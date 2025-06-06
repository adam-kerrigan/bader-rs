use crate::grid::Grid;
use rustc_hash::FxHashSet;
use std::cell::UnsafeCell;
use std::ops::{Deref, DerefMut};
use std::sync::atomic::{AtomicBool, AtomicIsize, Ordering};

/// Describes the state of the voxel.
pub enum Voxel<'a> {
    /// Contians the position of the voxel's maxima.
    Maxima(usize),
    /// Contians a vector of the maxima the current voxel contributes to and
    /// their weights.
    Boundary(&'a [f64]),
    /// A voxel beneath the vacuum tolerance and not contributing to any maxima.
    Vacuum,
}

/// A lock guard for write access to [`VoxelMap.weight_map`].
pub struct Lock<'a> {
    data: &'a BlockingVoxelMap,
}

unsafe impl Sync for Lock<'_> {}

/// Deref only exposes the weight_map field of a [`VoxelMap`].
impl Deref for Lock<'_> {
    type Target = Vec<Box<[f64]>>;
    fn deref(&self) -> &Vec<Box<[f64]>> {
        unsafe { &*self.data.weight_map.get() }
    }
}

/// DerefMut only exposes the weight_map field of a [`VoxelMap`].
impl DerefMut for Lock<'_> {
    fn deref_mut(&mut self) -> &mut Vec<Box<[f64]>> {
        unsafe { &mut *self.data.weight_map.get() }
    }
}

/// Make sure to free the lock when the struct is dropped.
impl Drop for Lock<'_> {
    fn drop(&mut self) {
        self.data.lock.store(false, Ordering::SeqCst);
    }
}

/// A structure for building and processing the map between voxel and maxima.
/// Bader maxima are stored in the voxel_map whilst the contributing weights are
/// stored in the weight_map. The weight_map is only written to once by each
/// point and so once a value has been written it is safe to read by any thread.
/// To check it has been written to `weight_get` monitors the state of corresponding
/// voxel_map value. Writing to the map is acheived by acquiring the lock, noting
/// the length of the weight_map, pushing the weight vector for voxel p to the
/// weight_map, droping the write lock and then storing the index of the inserted
/// vector using `weight_store`.
///
/// # Examples
/// ```
/// use bader::voxel_map::BlockingVoxelMap;
///
/// for p in 0..1isize {
///     let voxel_map = BlockingVoxelMap::new(
///         [2, 5, 2],
///         [[2.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 2.0]],
///         [0.0, 0.0, 0.0],
///     );
///     let i = {
///         let mut weight = voxel_map.lock();
///         (*weight).push(Vec::with_capacity(0).into());
///         weight.len() - 1
///     };
///     voxel_map.weight_store(p, i)
/// }
/// ```
pub struct BlockingVoxelMap {
    weight_map: UnsafeCell<Vec<Box<[f64]>>>,
    voxel_map: Vec<AtomicIsize>,
    pub grid: Grid,
    lock: AtomicBool,
}

unsafe impl Sync for BlockingVoxelMap {}

impl BlockingVoxelMap {
    /// Initialises a [`BlockingVoxelMap`] and the [`Grid`] that will faciliate movemoment around the
    /// map.
    pub fn new(
        grid: [usize; 3],
        lattice: [[f64; 3]; 3],
        voxel_origin: [f64; 3],
    ) -> Self {
        let grid = Grid::new(grid, lattice, voxel_origin);
        let size = grid.size.total;
        // For mapping the the voxels
        let weight_map =
            UnsafeCell::new(Vec::<Box<[f64]>>::with_capacity(size));
        let mut voxel_map = Vec::with_capacity(size);
        voxel_map.resize_with(size, || AtomicIsize::new(-1));
        let lock = AtomicBool::new(false);
        // For post processing
        Self {
            weight_map,
            voxel_map,
            grid,
            lock,
        }
    }

    /// Retrieves the state of the voxel, p. This will lock until p has been stored
    /// in the VoxelMap and then return either a `Voxel::Maxima` or `Voxel::Weight`.
    /// Calling this on a voxel, p, that is below the vacuum_tolerance will deadlock
    /// as a voxel is considered stored once voxel_map\[p\] > -1.
    pub fn weight_get(&self, i: isize) -> &[f64] {
        let i = -2 - i;
        &(unsafe { &*self.weight_map.get() })[i as usize]
    }

    /// Atomic loading of voxel, p, from voxel_map blocks if maxima == -1
    pub fn maxima_get(&self, p: isize) -> isize {
        loop {
            match self.voxel_map[p as usize].load(Ordering::Relaxed) {
                -1 => (),
                x => break x,
            }
        }
    }

    /// Check if a maxima is stored
    pub fn maxima_check(&self, p: isize) -> Option<isize> {
        match self.voxel_map[p as usize].load(Ordering::Relaxed) {
            -1 => None,
            x => Some(x),
        }
    }

    /// Stores the maxima of voxel, p, in the voxel_map.
    pub fn maxima_store(&self, p: isize, maxima: isize) {
        self.voxel_map[p as usize].store(maxima, Ordering::Relaxed);
    }

    /// Stores the index of p's weight contributions in weight_map into the
    /// weight_index.
    pub fn weight_store(&self, p: isize, i: usize) {
        self.maxima_store(p, -2 - (i as isize));
    }

    /// Locks the structure for write access unlock occurs when the returned
    /// Lock is dropped.
    pub fn lock(&self) -> Lock {
        while self.lock.swap(true, Ordering::SeqCst) {}
        Lock { data: self }
    }

    /// Extract the voxel map data.
    pub fn into_inner(self) -> (Vec<isize>, Vec<Box<[f64]>>, Grid) {
        (
            self.voxel_map.into_iter().map(|x| x.into_inner()).collect(),
            self.weight_map.into_inner(),
            self.grid,
        )
    }
}

/// A VoxelMap for if the maxima stored are atomic indices.
pub struct VoxelMap {
    /// The vector mapping the voxel to a maxima.
    pub voxel_map: Vec<isize>,
    /// The vector containing the weights for boundary voxels.
    pub weight_map: Vec<Box<[f64]>>,
    /// The Grid used to navigate the VoxelMap.
    pub grid: Grid,
}

impl VoxelMap {
    /// Create a new [`VoxelMap`]
    pub fn new(
        voxel_map: Vec<isize>,
        weight_map: Vec<Box<[f64]>>,
        grid: Grid,
    ) -> Self {
        Self {
            voxel_map,
            weight_map,
            grid,
        }
    }

    /// Create a new [`VoxelMap`] from a [`BlockingVoxelMap`].
    pub fn from_blocking_voxel_map(voxel_map: BlockingVoxelMap) -> Self {
        let (voxel_map, weight_map, grid) = voxel_map.into_inner();
        Self::new(voxel_map, weight_map, grid)
    }

    /// Produce an Iter over the boundary voxels.
    pub fn weight_iter(&self) -> std::slice::Iter<'_, Box<[f64]>> {
        self.weight_map.iter()
    }

    /// Get the length of the weight_map.
    pub fn weight_len(&self) -> usize {
        self.weight_map.len()
    }

    /// Get a refernce to the grid used by the VoxelMap.
    pub fn grid_get(&self) -> &Grid {
        &self.grid
    }

    /// Returns the atom associated with the point.
    pub fn maxima_to_atom(&self, maxima: usize) -> usize {
        maxima
    }

    /// Retrieval of the state of the voxel, p.
    pub fn maxima_to_voxel(&self, maxima: isize) -> Voxel {
        match maxima.cmp(&-1) {
            std::cmp::Ordering::Equal => Voxel::Vacuum,
            std::cmp::Ordering::Greater => Voxel::Maxima(maxima as usize),
            std::cmp::Ordering::Less => {
                Voxel::Boundary(self.maxima_to_weight(maxima))
            }
        }
    }

    /// Return a reference to the weights from the given maxima, Note: maxima here must be < -1.
    pub fn maxima_to_weight(&self, maxima: isize) -> &[f64] {
        &self.weight_map[(-2 - maxima) as usize]
    }

    /// Return an Iter over the maxima stored in the VoxelMap.
    pub fn maxima_iter(&self) -> std::slice::Iter<'_, isize> {
        self.voxel_map.iter()
    }

    /// Get the length of the voxel_map.
    pub fn maxima_len(&self) -> usize {
        self.voxel_map.len()
    }

    /// Return a Chunk over the maxima stored in the VoxelMap.
    pub fn maxima_chunks(
        &self,
        chunk_size: usize,
    ) -> std::slice::Chunks<'_, isize> {
        self.voxel_map.chunks(chunk_size)
    }

    /// Retrieval of the state of the voxel, p.
    pub fn voxel_get(&self, p: isize) -> Voxel {
        self.maxima_to_voxel(self.maxima_get(p))
    }

    /// Return the stored maxima at point p.
    pub fn maxima_get(&self, p: isize) -> isize {
        self.voxel_map[p as usize]
    }

    /// Produce a mask for a specific volume number.
    pub fn volume_map(&self, volume_number: isize) -> Vec<Option<f64>> {
        self.maxima_iter()
            .map(|maxima| {
                if *maxima == volume_number {
                    Some(1.0)
                } else if *maxima < -1 {
                    let mut w = None;
                    for weight in self.maxima_to_weight(*maxima).iter() {
                        let m = *weight as isize;
                        if m == volume_number {
                            w = Some(weight - m as f64);
                            break;
                        }
                    }
                    w
                } else {
                    None
                }
            })
            .collect()
    }
    /// Produce a mask for a collection volume numbers.
    pub fn multi_volume_map(
        &self,
        volume_numbers: &FxHashSet<isize>,
    ) -> Vec<Option<f64>> {
        self.maxima_iter()
            .map(|maxima| {
                if volume_numbers.contains(maxima) {
                    Some(1.0)
                } else if *maxima < -1 {
                    let mut w = 0.0;
                    for weight in self.maxima_to_weight(*maxima).iter() {
                        let m = *weight as isize;
                        if volume_numbers.contains(&m) {
                            w += weight - m as f64;
                        }
                    }
                    Some(w)
                } else {
                    None
                }
            })
            .collect()
    }
}
