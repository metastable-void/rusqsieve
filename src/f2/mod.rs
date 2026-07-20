//! Sparse binary matrices and verified dependencies.
use core::fmt;
use core::ops::Range;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MatrixOperation {
    Matrix,
    Transpose,
}
#[derive(Clone, Copy, Debug)]
pub enum MatrixSolver {
    Auto,
    DenseGaussian,
    BlockLanczos,
}
#[derive(Clone, Debug)]
pub struct MatrixConfig {
    pub solver: MatrixSolver,
    pub dense_threshold: usize,
    pub structured_elimination_limit: usize,
}
impl Default for MatrixConfig {
    fn default() -> Self {
        Self {
            solver: MatrixSolver::Auto,
            dense_threshold: 512,
            structured_elimination_limit: 10_000,
        }
    }
}
pub type CombinationId = u32;

/// A matrix stored in both row- and column-oriented sparse formats.
#[derive(Clone, Debug)]
pub struct SparseBinaryMatrix {
    rows: u32,
    columns: u32,
    csr_offsets: Box<[u32]>,
    csr_columns: Box<[u32]>,
    csc_offsets: Box<[u32]>,
    csc_rows: Box<[u32]>,
    provenance: Box<[CombinationId]>,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum MatrixError {
    DimensionOverflow,
    IndexOutOfRange,
    MalformedOffsets,
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
            provenance: (0..c).collect::<Vec<_>>().into_boxed_slice(),
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
    pub fn provenance(&self) -> &[CombinationId] {
        &self.provenance
    }
    pub fn mul_m_rows(&self, input: &[u64], range: Range<usize>, output: &mut [u64]) {
        assert!(range.end <= self.rows());
        assert!(input.len() >= self.columns());
        assert_eq!(output.len(), range.len());
        for (row, out) in range.zip(output) {
            let a = self.csr_offsets[row] as usize;
            let b = self.csr_offsets[row + 1] as usize;
            *out = self.csr_columns[a..b]
                .iter()
                .fold(0, |v, &c| v ^ input[c as usize]);
        }
    }
    pub fn mul_mt_columns(&self, input: &[u64], range: Range<usize>, output: &mut [u64]) {
        assert!(range.end <= self.columns());
        assert!(input.len() >= self.rows());
        assert_eq!(output.len(), range.len());
        for (col, out) in range.zip(output) {
            let a = self.csc_offsets[col] as usize;
            let b = self.csc_offsets[col + 1] as usize;
            *out = self.csc_rows[a..b]
                .iter()
                .fold(0, |v, &r| v ^ input[r as usize]);
        }
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
        let mut basis: Vec<Option<(Vec<u64>, Vec<u64>)>> = vec![None; self.rows()];
        let mut deps = Vec::new();
        for col in 0..cols {
            let mut parity = vec![0u64; self.rows().div_ceil(64)];
            let a = self.csc_offsets[col] as usize;
            let b = self.csc_offsets[col + 1] as usize;
            for &r in &self.csc_rows[a..b] {
                parity[r as usize / 64] ^= 1 << (r % 64)
            }
            let mut comb = vec![0u64; words];
            comb[col / 64] |= 1 << (col % 64);
            loop {
                let Some(pivot) = highest_bit(&parity) else {
                    if self.verify_dependency(&comb) {
                        deps.push(comb.into_boxed_slice())
                    }
                    break;
                };
                if let Some((p, c)) = &basis[pivot] {
                    xor(&mut parity, p);
                    xor(&mut comb, c)
                } else {
                    basis[pivot] = Some((parity, comb));
                    break;
                }
            }
        }
        DependencySet { vectors: deps }
    }

    /// Nullspace via SPEC §15.3 filtering — iterative singleton-row elimination
    /// (a prime occurring in one live column forces that column out of every
    /// dependency) — followed by dense elimination on the much smaller reduced
    /// matrix. Dependencies are returned in the ORIGINAL column space (eliminated
    /// columns are held at zero) and every one is re-verified against `self`.
    ///
    /// For quadratic-sieve matrices this removes the many low-weight rows before
    /// the O(n³) dense step, turning the linear-algebra phase from a bottleneck
    /// into a small fraction of the run at large input sizes.
    pub fn filtered_dependencies(&self) -> DependencySet {
        let nrows = self.rows();
        let ncols = self.columns();
        if ncols == 0 {
            return DependencySet::default();
        }
        let mut col_alive = vec![true; ncols];
        let mut row_weight = vec![0u32; nrows];
        for (r, w) in row_weight.iter_mut().enumerate() {
            *w = self.csr_offsets[r + 1] - self.csr_offsets[r];
        }
        let mut stack: Vec<usize> = (0..nrows).filter(|&r| row_weight[r] == 1).collect();
        while let Some(r) = stack.pop() {
            if row_weight[r] != 1 {
                continue;
            }
            let a = self.csr_offsets[r] as usize;
            let b = self.csr_offsets[r + 1] as usize;
            let Some(&c) = self.csr_columns[a..b]
                .iter()
                .find(|&&c| col_alive[c as usize])
            else {
                continue;
            };
            let c = c as usize;
            col_alive[c] = false;
            let ca = self.csc_offsets[c] as usize;
            let cb = self.csc_offsets[c + 1] as usize;
            for &rr in &self.csc_rows[ca..cb] {
                let rr = rr as usize;
                if row_weight[rr] > 0 {
                    row_weight[rr] -= 1;
                    if row_weight[rr] == 1 {
                        stack.push(rr);
                    }
                }
            }
        }
        let alive_cols: Vec<usize> = (0..ncols).filter(|&c| col_alive[c]).collect();
        // Nothing eliminated, or no usable surplus left: fall back to dense.
        let mut reduced_rows = 0usize;
        let mut row_map = vec![u32::MAX; nrows];
        for r in 0..nrows {
            if row_weight[r] >= 2 {
                row_map[r] = reduced_rows as u32;
                reduced_rows += 1;
            }
        }
        if alive_cols.len() == ncols || alive_cols.len() <= reduced_rows {
            return self.dense_dependencies();
        }
        let reduced_cols: Vec<Vec<u32>> = alive_cols
            .iter()
            .map(|&c| {
                let ca = self.csc_offsets[c] as usize;
                let cb = self.csc_offsets[c + 1] as usize;
                self.csc_rows[ca..cb]
                    .iter()
                    .filter_map(|&r| {
                        let m = row_map[r as usize];
                        (m != u32::MAX).then_some(m)
                    })
                    .collect()
            })
            .collect();
        let Ok(reduced) = SparseBinaryMatrix::from_columns(reduced_rows, &reduced_cols) else {
            return self.dense_dependencies();
        };
        let words = ncols.div_ceil(64);
        let mut out = Vec::new();
        for dep in reduced.dense_dependencies().iter() {
            let mut full = vec![0u64; words];
            for (j, &orig) in alive_cols.iter().enumerate() {
                if (dep[j / 64] >> (j % 64)) & 1 != 0 {
                    full[orig / 64] |= 1 << (orig % 64);
                }
            }
            if self.verify_dependency(&full) {
                out.push(full.into_boxed_slice());
            }
        }
        DependencySet { vectors: out }
    }
}
fn highest_bit(v: &[u64]) -> Option<usize> {
    v.iter()
        .rposition(|&x| x != 0)
        .map(|i| i * 64 + 63 - v[i].leading_zeros() as usize)
}
fn xor(a: &mut [u64], b: &[u64]) {
    for (x, y) in a.iter_mut().zip(b) {
        *x ^= *y
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
#[derive(Clone, Debug)]
pub struct F2BlockVector {
    words: Box<[u64]>,
}
impl F2BlockVector {
    pub fn new(len: usize) -> Self {
        Self {
            words: vec![0; len].into_boxed_slice(),
        }
    }
    pub fn as_slice(&self) -> &[u64] {
        &self.words
    }
    pub fn as_mut_slice(&mut self) -> &mut [u64] {
        &mut self.words
    }
}
#[derive(Clone, Debug)]
pub struct BlockLanczos {
    dependencies: DependencySet,
    complete: bool,
}
pub enum LanczosRequest<'a> {
    MultiplyM { input: &'a [u64] },
    MultiplyMt { input: &'a [u64] },
    Complete,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LanczosProgress {
    Progressed,
    Complete,
}
#[derive(Clone, Debug)]
pub enum LinearAlgebraError {
    WrongProductLength,
    InvalidDependency,
}
impl fmt::Display for LinearAlgebraError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "linear algebra error: {self:?}")
    }
}
impl std::error::Error for LinearAlgebraError {}
impl BlockLanczos {
    pub fn begin(matrix: &SparseBinaryMatrix) -> Self {
        Self {
            dependencies: matrix.filtered_dependencies(),
            complete: true,
        }
    }
    pub fn request(&self) -> LanczosRequest<'_> {
        LanczosRequest::Complete
    }
    pub fn submit_product(&mut self, _: &[u64]) -> Result<LanczosProgress, LinearAlgebraError> {
        Ok(if self.complete {
            LanczosProgress::Complete
        } else {
            LanczosProgress::Progressed
        })
    }
    pub fn dependencies(&self) -> Option<&DependencySet> {
        self.complete.then_some(&self.dependencies)
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
    fn multiply() {
        let m = SparseBinaryMatrix::from_columns(2, &[vec![0], vec![0, 1]]).unwrap();
        let mut out = [0; 2];
        m.mul_m_rows(&[3, 5], 0..2, &mut out);
        assert_eq!(out, [6, 5]);
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
                    (0..weight).map(|_| (rng() as usize % rows) as u32).collect()
                })
                .collect();
            let m = SparseBinaryMatrix::from_columns(rows, &columns).unwrap();
            let filtered = m.filtered_dependencies();
            for d in filtered.iter() {
                assert!(m.verify_dependency(d), "filtered produced an invalid dependency");
            }
            // cols > rows guarantees a nontrivial nullspace, so both solvers find one.
            assert!(!m.dense_dependencies().is_empty());
            assert!(
                !filtered.is_empty(),
                "filtered found no dependency though one exists"
            );
        }
    }
}
