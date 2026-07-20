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
            dependencies: matrix.dense_dependencies(),
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
}
