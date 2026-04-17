//! Core 2D field and periodic indexing for the thin viscous sheet solver.

/// Row-major nx-by-ny field of f64 values. Row stride is `nx` (number of
/// columns per row). Indexing convention: `data[j * nx + i]`.
#[derive(Clone)]
pub struct Field2D {
    data: Vec<f64>,
    nx: usize,
    ny: usize,
}

impl Field2D {
    pub fn new(nx: usize, ny: usize) -> Self {
        Self { data: vec![0.0; nx * ny], nx, ny }
    }

    pub fn filled(nx: usize, ny: usize, value: f64) -> Self {
        Self { data: vec![value; nx * ny], nx, ny }
    }

    #[inline(always)]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[j * self.nx + i]
    }

    #[inline(always)]
    pub fn set(&mut self, i: usize, j: usize, val: f64) {
        self.data[j * self.nx + i] = val;
    }

    #[inline(always)]
    pub fn data(&self) -> &[f64] {
        &self.data
    }

    #[inline(always)]
    pub fn data_mut(&mut self) -> &mut [f64] {
        &mut self.data
    }

    #[inline(always)]
    pub fn nx(&self) -> usize {
        self.nx
    }

    #[inline(always)]
    pub fn ny(&self) -> usize {
        self.ny
    }
}

/// Precomputed periodic (wrapping) index lookup for an N-sized dimension.
pub struct PeriodicIndex {
    n: usize,
    prev: Vec<usize>,
    next: Vec<usize>,
}

impl PeriodicIndex {
    pub fn new(n: usize) -> Self {
        let prev: Vec<usize> = (0..n).map(|i| (i + n - 1) % n).collect();
        let next: Vec<usize> = (0..n).map(|i| (i + 1) % n).collect();
        Self { n, prev, next }
    }

    #[inline(always)]
    pub fn prev(&self, i: usize) -> usize {
        self.prev[i]
    }

    #[inline(always)]
    pub fn next(&self, i: usize) -> usize {
        self.next[i]
    }

    #[inline(always)]
    pub fn n(&self) -> usize {
        self.n
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn periodic_prev_of_zero_wraps() {
        let idx = PeriodicIndex::new(16);
        assert_eq!(idx.prev(0), 15);
    }

    #[test]
    fn periodic_next_of_last_wraps() {
        let idx = PeriodicIndex::new(16);
        assert_eq!(idx.next(15), 0);
    }

    #[test]
    fn periodic_roundtrip() {
        let idx = PeriodicIndex::new(32);
        for i in 0..32 {
            assert_eq!(idx.prev(idx.next(i)), i);
            assert_eq!(idx.next(idx.prev(i)), i);
        }
    }

    #[test]
    fn field2d_basic() {
        let mut f = Field2D::new(4, 4);
        f.set(2, 3, 42.0);
        assert_eq!(f.get(2, 3), 42.0);
        assert_eq!(f.get(0, 0), 0.0);
        assert_eq!(f.nx(), 4);
        assert_eq!(f.ny(), 4);
        assert_eq!(f.data().len(), 16);
    }

    #[test]
    fn field2d_filled() {
        let f = Field2D::filled(8, 8, 1.5);
        for val in f.data() {
            assert_eq!(*val, 1.5);
        }
    }

    /// On a rectangular field the row stride is `nx`, so `(i, j)` lives at
    /// raw index `j * nx + i`. A bug that uses `ny` as the stride is
    /// invisible on square grids but corrupts data here.
    #[test]
    fn field2d_rectangular_stride() {
        let mut f = Field2D::new(8, 4);
        assert_eq!(f.nx(), 8);
        assert_eq!(f.ny(), 4);
        assert_eq!(f.data().len(), 32);

        f.set(7, 3, 42.0);
        f.set(5, 2, 7.5);
        assert_eq!(f.get(7, 3), 42.0);
        assert_eq!(f.get(5, 2), 7.5);

        // The stride must be nx=8, so (7, 3) is at flat index 3*8 + 7 = 31
        // and (5, 2) is at 2*8 + 5 = 21. The corresponding stride-swap
        // indices (7*4 + 3 = 31, 5*4 + 2 = 22) happen to coincide or
        // collide elsewhere — verifying the raw buffer catches both.
        assert_eq!(f.data()[3 * 8 + 7], 42.0);
        assert_eq!(f.data()[2 * 8 + 5], 7.5);
    }

    /// Coprime dimensions so that swapping the stride lands at a different
    /// (but still valid) buffer position. A stride-swap bug would put the
    /// value somewhere other than `j * nx + i` and this test would catch it.
    #[test]
    fn field2d_coprime_dimensions_stride() {
        let mut f = Field2D::new(7, 3);
        assert_eq!(f.data().len(), 21);

        // Writing to (5, 1) with correct stride lands at 1*7 + 5 = 12.
        // A buggy `i * ny + j` would compute 5*3 + 1 = 16, a distinct
        // index that is also in bounds.
        f.set(5, 1, 99.0);
        assert_eq!(f.get(5, 1), 99.0);
        assert_eq!(f.data()[12], 99.0);
        assert_eq!(f.data()[16], 0.0);
    }
}
