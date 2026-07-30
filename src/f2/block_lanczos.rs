//! Montgomery's 64-way block-Lanczos recurrence over `GF(2)`.
//!
//! A long vector of `u64` values represents 64 vectors at once: bit `j` of
//! element `i` is coordinate `i` of vector `j`. The recurrence is the one in
//! Peter Montgomery, "A Block Lanczos Algorithm for Finding Dependencies over
//! GF(2)", EUROCRYPT '95. In particular, this is an iterative sparse
//! `Bᵀ(BV)` solver, not blocked Gaussian elimination.

use super::{DependencySet, SparseBinaryMatrix};

const BLOCK: usize = 64;
const ALL: u64 = u64::MAX;

type BlockMatrix = [u64; BLOCK];

impl SparseBinaryMatrix {
    /// Find up to `limit` dependencies with the 64-way Montgomery block-Lanczos
    /// recurrence. A few deterministic starting blocks are tried because the
    /// algorithm has a small, input-dependent probability of breakdown.
    pub(super) fn block_lanczos_dependencies(&self, limit: usize) -> DependencySet {
        if limit == 0 || self.columns() == 0 {
            return DependencySet::default();
        }
        for attempt in 0..4u64 {
            if let Some(mut dependencies) = Lanczos::new(
                self,
                0xd6e8_feb8_6659_fd93 ^ attempt.wrapping_mul(0x9e37_79b9_7f4a_7c15),
            )
            .solve()
            {
                dependencies.truncate(limit);
                if !dependencies.is_empty()
                    && dependencies
                        .iter()
                        .all(|dependency| self.verify_dependency(dependency))
                {
                    return DependencySet {
                        vectors: dependencies,
                    };
                }
            }
        }
        DependencySet::default()
    }

    /// `out = B * input`, where `input` and every output coordinate carry 64
    /// independent vectors in their bits.
    fn mul_b(&self, input: &[u64], out: &mut [u64]) {
        debug_assert_eq!(input.len(), self.columns());
        debug_assert_eq!(out.len(), self.rows());
        out.fill(0);
        for (column, &word) in input.iter().enumerate() {
            if word == 0 {
                continue;
            }
            let start = self.csc_offsets[column] as usize;
            let end = self.csc_offsets[column + 1] as usize;
            for &row in &self.csc_rows[start..end] {
                out[row as usize] ^= word;
            }
        }
    }

    /// `out = Bᵀ * input`.
    fn mul_bt(&self, input: &[u64], out: &mut [u64]) {
        debug_assert_eq!(input.len(), self.rows());
        debug_assert_eq!(out.len(), self.columns());
        for (column, value) in out.iter_mut().enumerate() {
            let start = self.csc_offsets[column] as usize;
            let end = self.csc_offsets[column + 1] as usize;
            *value = self.csc_rows[start..end]
                .iter()
                .fold(0, |word, &row| word ^ input[row as usize]);
        }
    }
}

struct Lanczos<'a> {
    matrix: &'a SparseBinaryMatrix,
    seed: u64,
    row_scratch: Vec<u64>,
}

impl<'a> Lanczos<'a> {
    fn new(matrix: &'a SparseBinaryMatrix, seed: u64) -> Self {
        Self {
            matrix,
            seed,
            row_scratch: vec![0; matrix.rows()],
        }
    }

    /// Apply the symmetric matrix `A = BᵀB` without materializing it.
    fn apply_a(&mut self, input: &[u64], out: &mut [u64]) {
        self.matrix.mul_b(input, &mut self.row_scratch);
        self.matrix.mul_bt(&self.row_scratch, out);
    }

    fn random_word(&mut self) -> u64 {
        // SplitMix64 is deterministic on every target and has no all-zero
        // absorbing state, which makes it suitable for starting blocks.
        self.seed = self.seed.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut value = self.seed;
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    }

    fn solve(&mut self) -> Option<Vec<Box<[u64]>>> {
        let n = self.matrix.columns();
        let mut x: Vec<u64> = (0..n).map(|_| self.random_word()).collect();
        let mut v = [x.clone(), vec![0; n], vec![0; n]];
        let mut initial_v = vec![0; n];
        self.apply_a(&v[0], &mut initial_v);
        v[0].copy_from_slice(&initial_v);

        let mut v_next = vec![0; n];
        let mut winv = [[0; BLOCK]; 3];
        let mut vt_a_v = [[0; BLOCK]; 2];
        let mut vt_a2_v = [[0; BLOCK]; 2];
        let mut vt_v0 = [[0; BLOCK]; 3];
        let mut s_previous = core::array::from_fn(|index| index);
        let mut previous_dim = BLOCK;
        let mut previous_mask = ALL;
        let mut dimensions_solved = 0usize;
        let mut last_dim = 0usize;

        // A healthy run advances almost 64 dimensions per iteration. This
        // generous bound catches a recurrence breakdown without allowing an
        // accidental non-progressing run to consume unbounded time.
        let max_iterations = self.matrix.rows().div_ceil(24).saturating_add(128);
        for iteration in 1..=max_iterations {
            self.apply_a(&v[0], &mut v_next);
            vt_a_v[0] = inner_product(&v[0], &v_next);
            vt_a2_v[0] = inner_product(&v_next, &v_next);
            if vt_a_v[0].iter().all(|&row| row == 0) {
                if last_dim == 0 {
                    return None;
                }
                let mut bx = vec![0; self.matrix.rows()];
                let mut bv = vec![0; self.matrix.rows()];
                self.matrix.mul_b(&x, &mut bx);
                self.matrix.mul_b(&v[0], &mut bv);
                return Some(combine_columns(n, self.matrix.rows(), &x, &v[0], &bx, &bv));
            }

            let (current_dim, s_current, inverse) =
                find_nonsingular_sub(&vt_a_v[0], &s_previous, previous_dim)?;
            if current_dim == 0 {
                return None;
            }
            last_dim = current_dim;
            winv[0] = inverse;
            let current_mask = s_current[..current_dim]
                .iter()
                .fold(0u64, |mask, &column| mask | (1u64 << column));

            // Montgomery's recurrence requires every block column to occur in
            // the current or immediately preceding nonsingular submatrix. The
            // final iteration is exempt; its rank commonly falls short.
            if dimensions_solved < self.matrix.rows().saturating_sub(BLOCK)
                && current_mask | previous_mask != ALL
            {
                return None;
            }
            dimensions_solved += current_dim;
            if current_mask != ALL {
                for word in &mut v_next {
                    *word &= current_mask;
                }
            }

            if iteration < 4 {
                vt_v0[0] = inner_product(&v[0], &initial_v);
            }

            let mut d =
                core::array::from_fn(|row| (vt_a2_v[0][row] & current_mask) ^ vt_a_v[0][row]);
            d = matrix_mul(&winv[0], &d);
            xor_identity(&mut d);
            mul_block_acc(&v[0], &d, &mut v_next);
            let mut vt_v0_next = matrix_mul(&transpose(&d), &vt_v0[0]);

            let mut e = matrix_mul(&winv[1], &vt_a_v[0]);
            for row in &mut e {
                *row &= current_mask;
            }
            mul_block_acc(&v[1], &e, &mut v_next);
            xor_matrix(&mut vt_v0_next, &matrix_mul(&transpose(&e), &vt_v0[1]));

            // If the preceding block was rank-deficient, the third recurrence
            // term is essential; omitting it is a common "almost Lanczos" bug.
            if previous_mask != ALL {
                let mut f = matrix_mul(&vt_a_v[1], &winv[1]);
                xor_identity(&mut f);
                f = matrix_mul(&winv[2], &f);
                let f2 = core::array::from_fn(|row| {
                    ((vt_a2_v[1][row] & previous_mask) ^ vt_a_v[1][row]) & current_mask
                });
                f = matrix_mul(&f, &f2);
                mul_block_acc(&v[2], &f, &mut v_next);
                xor_matrix(&mut vt_v0_next, &matrix_mul(&transpose(&f), &vt_v0[2]));
            }

            d = matrix_mul(&winv[0], &vt_v0[0]);
            mul_block_acc(&v[0], &d, &mut x);

            v.rotate_right(1);
            core::mem::swap(&mut v[0], &mut v_next);
            winv.rotate_right(1);
            vt_v0.rotate_right(1);
            vt_v0[0] = vt_v0_next;
            vt_a_v.swap(0, 1);
            vt_a2_v.swap(0, 1);
            s_previous = s_current;
            previous_dim = current_dim;
            previous_mask = current_mask;
        }
        None
    }
}

/// `transpose(x) * y` for two `N × 64` bit matrices.
fn inner_product(x: &[u64], y: &[u64]) -> BlockMatrix {
    debug_assert_eq!(x.len(), y.len());
    let mut buckets = [[0u64; 256]; 8];
    for (&left, &right) in x.iter().zip(y) {
        for (byte, bucket) in buckets.iter_mut().enumerate() {
            bucket[((left >> (byte * 8)) & 0xff) as usize] ^= right;
        }
    }
    let mut product = [0u64; BLOCK];
    for bit in 0..8 {
        for (byte, bucket) in buckets.iter().enumerate() {
            let mut value = 0u64;
            for (index, &word) in bucket.iter().enumerate() {
                if index & (1 << bit) != 0 {
                    value ^= word;
                }
            }
            product[byte * 8 + bit] = value;
        }
    }
    product
}

/// XOR `vectors * matrix` into `out`.
fn mul_block_acc(vectors: &[u64], matrix: &BlockMatrix, out: &mut [u64]) {
    debug_assert_eq!(vectors.len(), out.len());
    let mut tables = [[0u64; 256]; 8];
    for (byte, table) in tables.iter_mut().enumerate() {
        for index in 1..256usize {
            let bit = index.trailing_zeros() as usize;
            table[index] = table[index & (index - 1)] ^ matrix[byte * 8 + bit];
        }
    }
    for (&word, output) in vectors.iter().zip(out) {
        let mut value = 0u64;
        for (byte, table) in tables.iter().enumerate() {
            value ^= table[((word >> (byte * 8)) & 0xff) as usize];
        }
        *output ^= value;
    }
}

fn matrix_mul(left: &BlockMatrix, right: &BlockMatrix) -> BlockMatrix {
    let mut product = [0u64; BLOCK];
    for (row, &word) in left.iter().enumerate() {
        let mut bits = word;
        while bits != 0 {
            let bit = bits.trailing_zeros() as usize;
            product[row] ^= right[bit];
            bits &= bits - 1;
        }
    }
    product
}

fn transpose(matrix: &BlockMatrix) -> BlockMatrix {
    let mut transposed = [0u64; BLOCK];
    for (row, &word) in matrix.iter().enumerate() {
        let mut bits = word;
        while bits != 0 {
            let column = bits.trailing_zeros() as usize;
            transposed[column] |= 1u64 << row;
            bits &= bits - 1;
        }
    }
    transposed
}

fn xor_identity(matrix: &mut BlockMatrix) {
    for (index, row) in matrix.iter_mut().enumerate() {
        *row ^= 1u64 << index;
    }
}

fn xor_matrix(left: &mut BlockMatrix, right: &BlockMatrix) {
    for (left, right) in left.iter_mut().zip(right) {
        *left ^= *right;
    }
}

/// Select and invert the maximal nonsingular submatrix required by the
/// Montgomery recurrence. Previously selected columns are considered last,
/// which preserves the current/previous block coverage invariant.
fn find_nonsingular_sub(
    input: &BlockMatrix,
    previous: &[usize; BLOCK],
    previous_dim: usize,
) -> Option<(usize, [usize; BLOCK], BlockMatrix)> {
    let mut left = *input;
    let mut right = core::array::from_fn(|row| 1u64 << row);
    let previous_mask = previous[..previous_dim]
        .iter()
        .fold(0u64, |mask, &column| mask | (1u64 << column));
    let mut columns = [0usize; BLOCK];
    for (offset, &column) in previous[..previous_dim].iter().enumerate() {
        columns[BLOCK - 1 - offset] = column;
    }
    let mut next = 0;
    for column in 0..BLOCK {
        if previous_mask & (1u64 << column) == 0 {
            columns[next] = column;
            next += 1;
        }
    }

    let mut selected = [0usize; BLOCK];
    let mut dimension = 0;
    for index in 0..BLOCK {
        let column = columns[index];
        let mask = 1u64 << column;
        let pivot = (index..BLOCK).find(|&candidate| left[columns[candidate]] & mask != 0);
        if let Some(pivot) = pivot {
            let pivot_row = columns[pivot];
            left.swap(column, pivot_row);
            right.swap(column, pivot_row);
            for &other in &columns {
                if other != column && left[other] & mask != 0 {
                    left[other] ^= left[column];
                    right[other] ^= right[column];
                }
            }
            selected[dimension] = column;
            dimension += 1;
            continue;
        }

        let pivot = (index..BLOCK).find(|&candidate| right[columns[candidate]] & mask != 0)?;
        let pivot_row = columns[pivot];
        left.swap(column, pivot_row);
        right.swap(column, pivot_row);
        for &other in &columns {
            if other != column && right[other] & mask != 0 {
                left[other] ^= left[column];
                right[other] ^= right[column];
            }
        }
        left[column] = 0;
        right[column] = 0;
    }
    Some((dimension, selected, right))
}

/// Convert the two terminal Lanczos blocks into actual nullspace vectors by
/// applying column elimination to `[B*x | B*v]` and mirroring it in `[x | v]`.
fn combine_columns(
    columns: usize,
    rows: usize,
    x: &[u64],
    v: &[u64],
    bx: &[u64],
    bv: &[u64],
) -> Vec<Box<[u64]>> {
    let words = columns.div_ceil(64);
    let image_words = rows.div_ceil(64);
    let mut vectors: Vec<Vec<u64>> = (0..128).map(|_| vec![0; words]).collect();
    let mut images: Vec<Vec<u64>> = (0..128).map(|_| vec![0; image_words]).collect();
    transpose_vectors(x, &mut vectors[..64]);
    transpose_vectors(v, &mut vectors[64..]);
    transpose_vectors(bx, &mut images[..64]);
    transpose_vectors(bv, &mut images[64..]);

    let mut rank = 0usize;
    for bit in 0..rows {
        let word = bit / 64;
        let mask = 1u64 << (bit % 64);
        let Some(pivot) = (rank..128).find(|&candidate| images[candidate][word] & mask != 0) else {
            continue;
        };
        images.swap(rank, pivot);
        vectors.swap(rank, pivot);
        for candidate in rank + 1..128 {
            if images[candidate][word] & mask == 0 {
                continue;
            }
            let (before, after) = images.split_at_mut(candidate);
            for (target, &source) in after[0].iter_mut().zip(&before[rank]) {
                *target ^= source;
            }
            let (before, after) = vectors.split_at_mut(candidate);
            for (target, &source) in after[0].iter_mut().zip(&before[rank]) {
                *target ^= source;
            }
        }
        rank += 1;
        if rank == 128 {
            break;
        }
    }
    if rank > 64 {
        return Vec::new();
    }
    vectors
        .into_iter()
        .take(64)
        .skip(rank)
        .filter(|dependency| dependency.iter().any(|&word| word != 0))
        .map(Vec::into_boxed_slice)
        .collect()
}

fn transpose_vectors(input: &[u64], output: &mut [Vec<u64>]) {
    debug_assert_eq!(output.len(), BLOCK);
    for (coordinate, &word) in input.iter().enumerate() {
        let mut bits = word;
        while bits != 0 {
            let vector = bits.trailing_zeros() as usize;
            output[vector][coordinate / 64] |= 1u64 << (coordinate % 64);
            bits &= bits - 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn block_helpers_match_scalar_arithmetic() {
        let mut state = 0x1234_5678_9abc_def0u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let left = core::array::from_fn(|_| next());
        let right = core::array::from_fn(|_| next());
        let product = matrix_mul(&left, &right);
        for row in 0..BLOCK {
            let expected = (0..BLOCK)
                .filter(|&column| left[row] & (1u64 << column) != 0)
                .fold(0, |value, column| value ^ right[column]);
            assert_eq!(product[row], expected);
        }
        assert_eq!(transpose(&transpose(&left)), left);

        let x: Vec<u64> = (0..137).map(|_| next()).collect();
        let y: Vec<u64> = (0..137).map(|_| next()).collect();
        let inner = inner_product(&x, &y);
        for (row, &actual) in inner.iter().enumerate() {
            let expected = x
                .iter()
                .zip(&y)
                .filter(|(word, _)| **word & (1u64 << row) != 0)
                .fold(0, |value, (_, &word)| value ^ word);
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn true_block_lanczos_finds_sparse_dependencies() {
        let rows = 192;
        let columns = 272;
        let mut state = 0xa076_1d64_78bd_642fu64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let matrix_columns: Vec<Vec<u32>> = (0..columns)
            .map(|_| {
                let mut column: Vec<u32> =
                    (0..7).map(|_| (next() as usize % rows) as u32).collect();
                column.sort_unstable();
                column.dedup();
                column
            })
            .collect();
        let matrix = SparseBinaryMatrix::from_columns(rows, &matrix_columns).unwrap();
        let dependencies = matrix.block_lanczos_dependencies(16);
        assert!(!dependencies.is_empty());
        assert!(
            dependencies
                .iter()
                .all(|dependency| matrix.verify_dependency(dependency))
        );
    }
}
