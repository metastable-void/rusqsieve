//! Sparse binary matrices and verified dependencies.
use core::fmt;

/// A matrix stored in both row- and column-oriented sparse formats.
#[derive(Clone, Debug)]
pub struct SparseBinaryMatrix {
    rows: u32,
    columns: u32,
    csr_offsets: Box<[u32]>,
    csr_columns: Box<[u32]>,
    csc_offsets: Box<[u32]>,
    csc_rows: Box<[u32]>,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum MatrixError {
    DimensionOverflow,
    IndexOutOfRange,
    MalformedOffsets,
    ResourceLimit,
}
impl fmt::Display for MatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "binary matrix error: {self:?}")
    }
}
impl std::error::Error for MatrixError {}
impl SparseBinaryMatrix {
    pub fn from_columns(rows: usize, columns: &[Vec<u32>]) -> Result<Self, MatrixError> {
        let r = u32::try_from(rows).map_err(|_| MatrixError::DimensionOverflow)?;
        let c = u32::try_from(columns.len()).map_err(|_| MatrixError::DimensionOverflow)?;
        let mut csc_o = Vec::with_capacity(columns.len() + 1);
        let mut csc_r = Vec::new();
        let mut rowcols = vec![Vec::new(); rows];
        csc_o.push(0);
        for (col, rs) in columns.iter().enumerate() {
            let mut sorted = rs.clone();
            sorted.sort_unstable();
            sorted.dedup();
            for &row in &sorted {
                if row >= r {
                    return Err(MatrixError::IndexOutOfRange);
                }
                csc_r.push(row);
                rowcols[row as usize].push(col as u32)
            }
            csc_o.push(u32::try_from(csc_r.len()).map_err(|_| MatrixError::DimensionOverflow)?)
        }
        let mut csr_o = Vec::with_capacity(rows + 1);
        let mut csr_c = Vec::new();
        csr_o.push(0);
        for cs in rowcols {
            csr_c.extend(cs);
            csr_o.push(u32::try_from(csr_c.len()).map_err(|_| MatrixError::DimensionOverflow)?)
        }
        Ok(Self {
            rows: r,
            columns: c,
            csr_offsets: csr_o.into_boxed_slice(),
            csr_columns: csr_c.into_boxed_slice(),
            csc_offsets: csc_o.into_boxed_slice(),
            csc_rows: csc_r.into_boxed_slice(),
        })
    }
    pub fn rows(&self) -> usize {
        self.rows as usize
    }
    pub fn columns(&self) -> usize {
        self.columns as usize
    }
    pub fn nonzeros(&self) -> usize {
        self.csc_rows.len()
    }
    pub fn verify_dependency(&self, selected: &[u64]) -> bool {
        if selected.len() < self.columns().div_ceil(64) {
            return false;
        }
        for row in 0..self.rows() {
            let a = self.csr_offsets[row] as usize;
            let b = self.csr_offsets[row + 1] as usize;
            if self.csr_columns[a..b].iter().fold(false, |v, &c| {
                v ^ ((selected[c as usize / 64] >> (c % 64)) & 1 != 0)
            }) {
                return false;
            }
        }
        true
    }
    pub fn dense_dependencies(&self) -> DependencySet {
        let cols = self.columns();
        let words = cols.div_ceil(64);
        let parity_words = self.rows().div_ceil(64);
        let mut basis: Vec<Option<(Vec<u64>, Vec<u64>)>> = vec![None; self.rows()];
        let mut deps = Vec::new();
        // Two working vectors reused across columns. Allocating them inside the loop cost one
        // allocation pair per column — 41 816 of them on a 256-bit matrix — and they are only moved
        // out on the two paths that consume them, where a fresh pair is taken for the next column.
        let mut parity = vec![0u64; parity_words];
        let mut comb = vec![0u64; words];
        for col in 0..cols {
            parity.clear();
            parity.resize(parity_words, 0);
            comb.clear();
            comb.resize(words, 0);
            let a = self.csc_offsets[col] as usize;
            let b = self.csc_offsets[col + 1] as usize;
            for &r in &self.csc_rows[a..b] {
                parity[r as usize / 64] ^= 1 << (r % 64)
            }
            comb[col / 64] |= 1 << (col % 64);
            loop {
                let Some(pivot) = highest_bit(&parity) else {
                    if self.verify_dependency(&comb) {
                        deps.push(core::mem::take(&mut comb).into_boxed_slice())
                    }
                    break;
                };
                if let Some((p, c)) = &basis[pivot] {
                    xor(&mut parity, p);
                    xor(&mut comb, c)
                } else {
                    basis[pivot] = Some((core::mem::take(&mut parity), core::mem::take(&mut comb)));
                    break;
                }
            }
        }
        DependencySet { vectors: deps }
    }

    /// Compute a bounded nullspace basis by row-reducing the parity matrix.
    ///
    /// The column-oriented reference solver above carries both a parity vector
    /// and a full provenance vector through every elimination.  Once sparse
    /// filtering has made the residual matrix fairly dense, reducing rows uses
    /// half as much live bitset data.  The echelon rows themselves are equations
    /// in the original column variables, so dependencies can be recovered by
    /// back-substitution without maintaining provenance during elimination.
    fn row_echelon_dependencies(&self, limit: usize) -> DependencySet {
        let cols = self.columns();
        if cols == 0 || limit == 0 {
            return DependencySet::default();
        }
        let words = cols.div_ceil(64);
        let mut basis: Vec<Option<Box<[u64]>>> = vec![None; cols];

        for row in 0..self.rows() {
            let a = self.csr_offsets[row] as usize;
            let b = self.csr_offsets[row + 1] as usize;
            if a == b {
                continue;
            }
            let highest_column = self.csr_columns[a..b].iter().copied().max().unwrap() as usize;
            let mut equation = vec![0u64; highest_column / 64 + 1];
            for &column in &self.csr_columns[a..b] {
                equation[column as usize / 64] ^= 1 << (column % 64);
            }
            // Reduce against the basis one pivot at a time, XORing only up to the pivot word: the
            // basis row's higher words are all zero by construction, so a truncated XOR is exact.
            // A Gray-code (M4RI) table over 4-column blocks was implemented and measured here and
            // was 75% SLOWER (f2_dense 2.34 s -> 4.08 s at 256-bit): pivots are installed
            // incrementally, so every insertion invalidates its block's 16-row table and the
            // rebuild cost dominates the XORs it saves. See CHANGELOG 0.2.1.
            while let Some(pivot) = highest_bit(&equation) {
                if let Some(prior) = &basis[pivot] {
                    xor(&mut equation[..=pivot / 64], &prior[..=pivot / 64]);
                } else {
                    equation.truncate(pivot / 64 + 1);
                    basis[pivot] = Some(equation.into_boxed_slice());
                    break;
                }
            }
        }

        let mut dependencies = Vec::new();
        for free in (0..cols)
            .filter(|&column| basis[column].is_none())
            .take(limit)
        {
            let mut dependency = vec![0u64; words];
            dependency[free / 64] |= 1 << (free % 64);
            // A pivot row has no set bits above its pivot.  Ascending
            // substitution therefore has every right-hand-side value ready.
            for (pivot, equation) in basis.iter().enumerate() {
                let Some(equation) = equation else {
                    continue;
                };
                let last = pivot / 64;
                let odd = parity_dot(&equation[..=last], &dependency[..=last]);
                if odd != 0 {
                    dependency[last] ^= 1 << (pivot % 64);
                }
            }
            if self.verify_dependency(&dependency) {
                dependencies.push(dependency.into_boxed_slice());
            }
        }
        DependencySet {
            vectors: dependencies,
        }
    }

    /// Nullspace via SPEC §8 filtering — iterative elimination of every row of weight 1 through
    /// `MAX_STRUCTURED_WEIGHT` (6), with Markowitz-style pivot selection to limit fill-in, not merely
    /// the singleton rows an earlier version of this comment described — followed by dense
    /// elimination on the much smaller reduced matrix. Dependencies are returned in the ORIGINAL column space (eliminated
    /// columns are held at zero) and every one is re-verified against `self`.
    ///
    /// For quadratic-sieve matrices this removes the many low-weight rows before
    /// the O(n³) dense step, turning the linear-algebra phase from a bottleneck
    /// into a small fraction of the run at large input sizes.
    pub fn filtered_dependencies(&self) -> Result<DependencySet, MatrixError> {
        let nrows = self.rows();
        let ncols = self.columns();
        if ncols == 0 {
            return Ok(DependencySet::default());
        }
        let mut row_cols: Vec<Vec<usize>> = (0..nrows)
            .map(|r| {
                let a = self.csr_offsets[r] as usize;
                let b = self.csr_offsets[r + 1] as usize;
                self.csr_columns[a..b].iter().map(|&c| c as usize).collect()
            })
            .collect();
        let mut col_rows: Vec<Vec<usize>> = (0..ncols)
            .map(|c| {
                let a = self.csc_offsets[c] as usize;
                let b = self.csc_offsets[c + 1] as usize;
                self.csc_rows[a..b].iter().map(|&r| r as usize).collect()
            })
            .collect();
        let mut col_alive = vec![true; ncols];
        // Each eliminated pivot satisfies x[pivot] = XOR(x[other] for other in rhs).
        // Replaying these records backwards expands a dependency of the reduced matrix into the
        // original column space without carrying dense provenance through the sparse phase.
        let mut eliminations: Vec<(usize, Vec<usize>)> = Vec::new();
        const MAX_STRUCTURED_WEIGHT: usize = 6;
        let mut stack: Vec<usize> = (0..nrows)
            .filter(|&r| (1..=MAX_STRUCTURED_WEIGHT).contains(&row_cols[r].len()))
            .collect();

        while let Some(r) = stack.pop() {
            let weight = row_cols[r].len();
            if weight == 0 || weight > MAX_STRUCTURED_WEIGHT {
                continue;
            }
            let equation: Vec<usize> = row_cols[r].to_vec();
            // Markowitz-style choice: eliminate the column occurring in the fewest other rows,
            // minimizing fill. Ties are stable because the equation is sorted.
            let pivot = *equation
                .iter()
                .min_by_key(|&&c| (col_rows[c].len(), c))
                .unwrap();
            let rhs: Vec<usize> = equation.iter().copied().filter(|&c| c != pivot).collect();

            row_cols[r].clear();
            for &c in &equation {
                sorted_remove(&mut col_rows[c], r);
            }
            let affected: Vec<usize> = col_rows[pivot].to_vec();
            for rr in affected {
                sorted_remove(&mut row_cols[rr], pivot);
                sorted_remove(&mut col_rows[pivot], rr);
                for &c in &rhs {
                    if sorted_remove(&mut row_cols[rr], c) {
                        sorted_remove(&mut col_rows[c], rr);
                    } else {
                        sorted_insert(&mut row_cols[rr], c);
                        sorted_insert(&mut col_rows[c], rr);
                    }
                }
                if (1..=MAX_STRUCTURED_WEIGHT).contains(&row_cols[rr].len()) {
                    stack.push(rr);
                }
            }
            col_rows[pivot].clear();
            col_alive[pivot] = false;
            eliminations.push((pivot, rhs));
        }

        let alive_cols: Vec<usize> = (0..ncols).filter(|&c| col_alive[c]).collect();
        let mut reduced_rows = 0usize;
        let mut row_map = vec![u32::MAX; nrows];
        for r in 0..nrows {
            if !row_cols[r].is_empty() {
                row_map[r] = reduced_rows as u32;
                reduced_rows += 1;
            }
        }
        if alive_cols.len() == ncols || alive_cols.len() <= reduced_rows {
            let dense_bytes = ncols.saturating_mul(nrows.div_ceil(64)).saturating_mul(16);
            if dense_bytes > 256 * 1024 * 1024 {
                return Err(MatrixError::ResourceLimit);
            }
            return Ok(self.dense_dependencies());
        }
        let reduced_cols: Vec<Vec<u32>> = alive_cols
            .iter()
            .map(|&c| {
                col_rows[c]
                    .iter()
                    .filter_map(|&r| {
                        let m = row_map[r];
                        (m != u32::MAX).then_some(m)
                    })
                    .collect()
            })
            .collect();
        let Ok(reduced) = SparseBinaryMatrix::from_columns(reduced_rows, &reduced_cols) else {
            return Err(MatrixError::MalformedOffsets);
        };
        let words = ncols.div_ceil(64);
        let mut out = Vec::new();
        // Cap the nullspace basis at 64 dependencies. Each one has an independent ~1/2 chance of
        // yielding a nontrivial gcd, so 64 makes exhausting them without a factor negligible, while
        // back-substitution is O(cols²/64) per dependency and there is no reason to compute the
        // hundreds the residual matrix usually admits. (This bound was previously justified by what
        // "block solvers conventionally return"; this crate has no block solver.)
        for dep in reduced.row_echelon_dependencies(64).iter() {
            let mut full = vec![0u64; words];
            for (j, &original_col) in alive_cols.iter().enumerate() {
                if (dep[j / 64] >> (j % 64)) & 1 != 0 {
                    full[original_col / 64] |= 1 << (original_col % 64);
                }
            }
            for (pivot, rhs) in eliminations.iter().rev() {
                let value = rhs
                    .iter()
                    .fold(false, |v, &c| v ^ ((full[c / 64] >> (c % 64)) & 1 != 0));
                if value {
                    full[pivot / 64] |= 1 << (pivot % 64);
                }
            }
            if self.verify_dependency(&full) {
                out.push(full.into_boxed_slice());
            }
        }
        Ok(DependencySet { vectors: out })
    }
}

fn sorted_remove(values: &mut Vec<usize>, value: usize) -> bool {
    match values.binary_search(&value) {
        Ok(index) => {
            values.remove(index);
            true
        }
        Err(_) => false,
    }
}

fn sorted_insert(values: &mut Vec<usize>, value: usize) {
    if let Err(index) = values.binary_search(&value) {
        values.insert(index, value);
    }
}
fn highest_bit(v: &[u64]) -> Option<usize> {
    v.iter()
        .rposition(|&x| x != 0)
        .map(|i| i * 64 + 63 - v[i].leading_zeros() as usize)
}
#[inline]
fn parity_dot(a: &[u64], b: &[u64]) -> u32 {
    let mut parity = [0u32; 4];
    let mut index = 0;
    while index + 4 <= a.len().min(b.len()) {
        for lane in 0..4 {
            parity[lane] ^= (a[index + lane] & b[index + lane]).count_ones();
        }
        index += 4;
    }
    let mut result = parity.into_iter().fold(0, |value, lane| value ^ lane);
    while index < a.len().min(b.len()) {
        result ^= (a[index] & b[index]).count_ones();
        index += 1;
    }
    result & 1
}
#[cfg(not(all(feature = "wasm-simd128", target_arch = "wasm32")))]
fn xor(a: &mut [u64], b: &[u64]) {
    for (x, y) in a.iter_mut().zip(b) {
        *x ^= *y
    }
}

#[cfg(all(feature = "wasm-simd128", target_arch = "wasm32"))]
#[allow(unsafe_code)]
fn xor(a: &mut [u64], b: &[u64]) {
    // The feature explicitly opts the whole wasm artifact into the simd128
    // baseline, so calling the specialized function is valid.
    unsafe { xor_wasm_simd(a, b) }
}

#[cfg(all(feature = "wasm-simd128", target_arch = "wasm32"))]
#[allow(unsafe_code)]
#[target_feature(enable = "simd128")]
unsafe fn xor_wasm_simd(a: &mut [u64], b: &[u64]) {
    use core::arch::wasm32::{v128, v128_load, v128_store, v128_xor};
    let len = a.len().min(b.len());
    let mut i = 0;
    while i + 2 <= len {
        // SAFETY: the loop condition proves that both slices contain the two
        // u64 lanes loaded here. WebAssembly v128 loads/stores are unaligned.
        unsafe {
            let av = v128_load(a.as_ptr().add(i).cast::<v128>());
            let bv = v128_load(b.as_ptr().add(i).cast::<v128>());
            v128_store(a.as_mut_ptr().add(i).cast::<v128>(), v128_xor(av, bv));
        }
        i += 2;
    }
    if i < len {
        a[i] ^= b[i];
    }
}
#[derive(Clone, Debug, Default)]
pub struct DependencySet {
    vectors: Vec<Box<[u64]>>,
}
impl DependencySet {
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &[u64]> {
        self.vectors.iter().map(AsRef::as_ref)
    }
    pub fn len(&self) -> usize {
        self.vectors.len()
    }
    pub fn is_empty(&self) -> bool {
        self.vectors.is_empty()
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn dense_dep() {
        let m = SparseBinaryMatrix::from_columns(3, &[vec![0, 1], vec![1, 2], vec![0, 2]]).unwrap();
        let d = m.dense_dependencies();
        assert_eq!(d.len(), 1);
        assert!(m.verify_dependency(d.iter().next().unwrap()));
    }

    #[test]
    fn filtered_dependencies_are_valid_and_present() {
        // Deterministic pseudo-random sparse matrices with a nullspace (cols > rows)
        // and plenty of singleton rows. Every filtered dependency must verify, and
        // when a dependency exists it must be found.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut rng = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..50 {
            let rows = 30 + (rng() as usize % 40);
            let cols = rows + 8 + (rng() as usize % 20);
            let columns: Vec<Vec<u32>> = (0..cols)
                .map(|_| {
                    let weight = 1 + (rng() as usize % 5);
                    (0..weight)
                        .map(|_| (rng() as usize % rows) as u32)
                        .collect()
                })
                .collect();
            let m = SparseBinaryMatrix::from_columns(rows, &columns).unwrap();
            let filtered = m.filtered_dependencies().unwrap();
            for d in filtered.iter() {
                assert!(
                    m.verify_dependency(d),
                    "filtered produced an invalid dependency"
                );
            }
            let dense = m.dense_dependencies();
            let echelon = m.row_echelon_dependencies(64);
            assert_eq!(echelon.len(), dense.len());
            for d in echelon.iter() {
                assert!(m.verify_dependency(d));
            }
            // cols > rows guarantees a nontrivial nullspace, so both solvers find one.
            assert!(!dense.is_empty());
            assert!(
                !filtered.is_empty(),
                "filtered found no dependency though one exists"
            );
        }
    }
}
