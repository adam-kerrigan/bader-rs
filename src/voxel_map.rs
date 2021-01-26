use std::cell::UnsafeCell;
use std::collections::BTreeSet;
use std::ops::{Deref, DerefMut};
use std::sync::atomic::{AtomicBool, AtomicIsize, Ordering};

/// Describes the state of the voxel.
pub enum Voxel<'a> {
    /// Contians the position of the voxel's maxima.
    Maxima(usize),
    /// Contians a vector of the maxima the current voxel contributes to and
    /// their weights.
    Weight(&'a Vec<(usize, f64)>),
    /// A voxel beneath the vacuum tolerance and not contributing to any maxima.
    Vacuum,
}

/// A lock guard for write access to VoxelMap.weight_map
pub struct Lock<'a> {
    data: &'a VoxelMap,
}

unsafe impl<'a> Sync for Lock<'a> {}

/// Deref only exposes the weight_map field of a [`VoxelMap`].
impl<'a> Deref for Lock<'a> {
    type Target = Vec<Vec<(usize, f64)>>;
    fn deref(&self) -> &Vec<Vec<(usize, f64)>> {
        unsafe { &*self.data.weight_map.get() }
    }
}

/// DerefMut only exposes the weight_map field of a [`VoxelMap`].
impl<'a> DerefMut for Lock<'a> {
    fn deref_mut(&mut self) -> &mut Vec<Vec<(usize, f64)>> {
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
/// use bader::voxel_map::VoxelMap;
///
/// for p in 0..1isize {
///     let voxel_map = VoxelMap::new(10);
///     let i = {
///         let mut weight = voxel_map.lock();
///         (*weight).push(Vec::with_capacity(0));
///         weight.len() - 1
///     };
///     voxel_map.weight_store(p, i)
/// }
/// ```
pub struct VoxelMap {
    weight_map: UnsafeCell<Vec<Vec<(usize, f64)>>>,
    voxel_map: Vec<AtomicIsize>,
    lock: AtomicBool,
}

unsafe impl Sync for VoxelMap {}

impl VoxelMap {
    /// Initialises a VoxelMap of dimensions, size.
    pub fn new(size: usize) -> Self {
        // For mapping the the voxels
        let weight_map =
            UnsafeCell::new(Vec::<Vec<(usize, f64)>>::with_capacity(size));
        let mut voxel_map = Vec::with_capacity(size);
        voxel_map.resize_with(size, || AtomicIsize::new(-1));
        let lock = AtomicBool::new(false);
        // For post processing
        Self { weight_map,
               voxel_map,
               lock }
    }

    /// Retrieves the state of the voxel, p. This will lock until p has been stored
    /// in the VoxelMap and then return either a `Voxel::Maxima` or `Voxel::Weight`.
    /// Calling this on a voxel, p, that is below the vacuum_tolerance will deadlock
    /// as a voxel is considered stored once voxel_map\[p\] > -1.
    pub fn weight_get(&self, i: isize) -> &Vec<(usize, f64)> {
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

    /// Atomic loading of voxel, p, from voxel_map
    pub fn maxima_non_block_get(&self, p: isize) -> isize {
        let volume_number = self.voxel_map[p as usize].load(Ordering::Relaxed);
        match volume_number.cmp(&-1) {
            std::cmp::Ordering::Equal => -1,
            std::cmp::Ordering::Greater => volume_number,
            std::cmp::Ordering::Less => {
                let weight = self.weight_get(volume_number);
                weight[0].0 as isize
            }
        }
    }

    /// A none locking retrieval of the state of voxel, p. This should only be
    /// used once the VoxelMap has been fully populated.
    pub fn voxel_get(&self, p: isize) -> Voxel {
        let volume_number = self.voxel_map[p as usize].load(Ordering::Relaxed);
        match volume_number.cmp(&-1) {
            std::cmp::Ordering::Equal => Voxel::Vacuum,
            std::cmp::Ordering::Greater => {
                Voxel::Maxima(volume_number as usize)
            }
            std::cmp::Ordering::Less => {
                let weight = self.weight_get(volume_number);
                Voxel::Weight(weight)
            }
        }
    }

    /// Finds the unique values contained in voxel_map and returns them sorted by
    /// value.
    pub fn maxima_list(&self) -> Vec<usize> {
        let len = self.voxel_map.len() as isize;
        let maximas = (0..len).filter_map(|x| {
                                  if let Voxel::Maxima(m) = self.voxel_get(x) {
                                      Some(m)
                                  } else {
                                      None
                                  }
                              })
                              .collect::<BTreeSet<usize>>();
        maximas.into_iter().collect()
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn voxel_map_maxima_store() {
        let voxel_map = VoxelMap::new(10);
        for p in 0..10isize {
            voxel_map.maxima_store(p, p - 1);
        }
        assert_eq!(voxel_map.maxima_non_block_get(0), -1);
        assert_eq!(voxel_map.maxima_non_block_get(9), 8);
    }

    #[test]
    fn voxel_map_weight_store() {}
}
