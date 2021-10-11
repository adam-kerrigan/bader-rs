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
    Weight(&'a Vec<f64>),
    /// A voxel beneath the vacuum tolerance and not contributing to any maxima.
    Vacuum,
}

/// A lock guard for write access to [`VoxelMap.weight_map`]
pub struct Lock<'a> {
    data: &'a BlockingVoxelMap,
}

unsafe impl<'a> Sync for Lock<'a> {}

/// Deref only exposes the weight_map field of a [`VoxelMap`].
impl<'a> Deref for Lock<'a> {
    type Target = Vec<Vec<f64>>;
    fn deref(&self) -> &Vec<Vec<f64>> {
        unsafe { &*self.data.weight_map.get() }
    }
}

/// DerefMut only exposes the weight_map field of a [`VoxelMap`].
impl<'a> DerefMut for Lock<'a> {
    fn deref_mut(&mut self) -> &mut Vec<Vec<f64>> {
        unsafe { &mut *self.data.weight_map.get() }
    }
}

/// Make sure to free the lock when the struct is dropped.
impl<'a> Drop for Lock<'a> {
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
/// use bader::voxel_map::BlockingVoxelMap as VoxelMap;
///
/// for p in 0..1isize {
///     let voxel_map = VoxelMap::new([2, 5, 2],
///                                   [[2.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 2.0]],
///                                   [0.0, 0.0, 0.0]);
///     let i = {
///         let mut weight = voxel_map.lock();
///         (*weight).push(Vec::with_capacity(0));
///         weight.len() - 1
///     };
///     voxel_map.weight_store(p, i)
/// }
/// ```
pub struct BlockingVoxelMap {
    weight_map: UnsafeCell<Vec<Vec<f64>>>,
    voxel_map: Vec<AtomicIsize>,
    pub grid: Grid,
    lock: AtomicBool,
}

unsafe impl Sync for BlockingVoxelMap {}

impl BlockingVoxelMap {
    /// Initialises a VoxelMap and the bader::grid::Grid that will faciliate movemoment around the
    /// map.
    pub fn new(grid: [usize; 3],
               lattice: [[f64; 3]; 3],
               voxel_origin: [f64; 3])
               -> Self {
        let grid = Grid::new(grid, lattice, voxel_origin);
        let size = grid.size.total;
        // For mapping the the voxels
        let weight_map = UnsafeCell::new(Vec::<Vec<f64>>::with_capacity(size));
        let mut voxel_map = Vec::with_capacity(size);
        voxel_map.resize_with(size, || AtomicIsize::new(-1));
        let lock = AtomicBool::new(false);
        // For post processing
        Self { weight_map,
               voxel_map,
               grid,
               lock }
    }

    /// Retrieves the state of the voxel, p. This will lock until p has been stored
    /// in the VoxelMap and then return either a `Voxel::Maxima` or `Voxel::Weight`.
    /// Calling this on a voxel, p, that is below the vacuum_tolerance will deadlock
    /// as a voxel is considered stored once voxel_map\[p\] > -1.
    pub fn weight_get(&self, i: isize) -> &Vec<f64> {
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
    pub fn into_inner(self) -> (Vec<isize>, Vec<Vec<f64>>, Grid) {
        (self.voxel_map.into_iter().map(|x| x.into_inner()).collect(),
         self.weight_map.into_inner(),
         self.grid)
    }
}

pub struct NonBlockingVoxelMap {
    pub voxel_map: Vec<isize>,
    pub weight_map: Vec<Vec<f64>>,
    pub grid: Grid,
}

impl NonBlockingVoxelMap {
    pub fn new(voxel_map: Vec<isize>,
               weight_map: Vec<Vec<f64>>,
               grid: Grid)
               -> Self {
        Self { voxel_map,
               weight_map,
               grid }
    }

    pub fn from_blocking_voxel_map(voxel_map: BlockingVoxelMap) -> Self {
        let (voxel_map, weight_map, grid) = voxel_map.into_inner();
        Self::new(voxel_map, weight_map, grid)
    }

    pub fn weight_get(&self, maxima: isize) -> &Vec<f64> {
        let i = -2 - maxima;
        &self.weight_map[i as usize]
    }
    /// Atomic loading of voxel, p, from voxel_map
    pub fn maxima_get(&self, p: isize) -> isize {
        let maxima = self.voxel_map[p as usize];
        match maxima.cmp(&-1) {
            std::cmp::Ordering::Equal => -1,
            std::cmp::Ordering::Greater => maxima,
            std::cmp::Ordering::Less => {
                let weight = self.weight_get(maxima);
                weight[0] as isize
            }
        }
    }

    /// A none locking retrieval of the state of voxel, p. This should only be
    /// used once the VoxelMap has been fully populated.
    pub fn voxel_get(&self, p: isize) -> Voxel {
        let maxima = self.voxel_map[p as usize];
        match maxima.cmp(&-1) {
            std::cmp::Ordering::Equal => Voxel::Vacuum,
            std::cmp::Ordering::Greater => Voxel::Maxima(maxima as usize),
            std::cmp::Ordering::Less => {
                let weight = self.weight_get(maxima);
                Voxel::Weight(weight)
            }
        }
    }

    pub fn volume_map(&self, volume_number: isize) -> Vec<Option<f64>> {
        self.voxel_map
            .iter()
            .map(|maxima| {
                if *maxima == volume_number {
                    Some(1.0)
                } else if *maxima < -1 {
                    let mut w = None;
                    for weight in self.weight_get(*maxima) {
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

    pub fn multi_volume_map(&self,
                            volume_numbers: &FxHashSet<isize>)
                            -> Vec<Option<f64>> {
        self.voxel_map
            .iter()
            .map(|maxima| {
                if volume_numbers.contains(maxima) {
                    Some(1.0)
                } else if *maxima < -1 {
                    let mut w = 0.0;
                    for weight in self.weight_get(*maxima) {
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
