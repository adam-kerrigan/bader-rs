use crate::atoms::Atoms;
use crate::density::Density;
use crate::progress::Bar;
use crate::utils;
use indicatif::ProgressBar;
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
    weight_index: Vec<AtomicIsize>,
    voxel_map: Vec<AtomicIsize>,
    lock: AtomicBool,
    pub assigned_atom: Vec<usize>,
    pub minimum_distance: Vec<f64>,
    pub bader_volume: Vec<f64>,
    pub surface_distance: Vec<f64>,
    pub bader_charge: Vec<Vec<f64>>,
    pub bader_maxima: Vec<isize>,
    maxima_index: Vec<usize>,
}

unsafe impl Sync for VoxelMap {}

impl VoxelMap {
    /// Initialises a VoxelMap of dimensions, size.
    pub fn new(size: usize) -> Self {
        // For mapping the the voxels
        let weight_map =
            UnsafeCell::new(Vec::<Vec<(usize, f64)>>::with_capacity(size));
        let mut weight_index = Vec::with_capacity(size);
        weight_index.resize_with(size, || AtomicIsize::new(-1));
        let mut voxel_map = Vec::with_capacity(size);
        voxel_map.resize_with(size, || AtomicIsize::new(-1));
        let lock = AtomicBool::new(false);
        // For post processing
        let assigned_atom = Vec::with_capacity(0);
        let maxima_index = Vec::with_capacity(0);
        let minimum_distance = Vec::with_capacity(0);
        let bader_maxima = Vec::with_capacity(0);
        let bader_charge = vec![Vec::with_capacity(0)];
        let bader_volume = Vec::with_capacity(0);
        let surface_distance = Vec::with_capacity(0);
        Self { weight_map,
               weight_index,
               voxel_map,
               lock,
               assigned_atom,
               maxima_index,
               minimum_distance,
               bader_maxima,
               surface_distance,
               bader_charge,
               bader_volume }
    }

    /// Retrieves the state of the voxel, p. This will lock until p has been stored
    /// in the VoxelMap and then return either a `Voxel::Maxima` or `Voxel::Weight`.
    /// Calling this on a voxel, p, that is below the vacuum_tolerance will deadlock
    /// as a voxel is considered stored once voxel_map\[p\] > -1.
    pub fn weight_get(&self, p: isize) -> Voxel {
        let volume_number = loop {
            match self.voxel_map[p as usize].load(Ordering::Relaxed) {
                -1 => (),
                x => break x as usize,
            }
        };
        let i = match self.weight_index[p as usize].load(Ordering::Relaxed) {
            -1 => return Voxel::Maxima(volume_number),
            i => i as usize,
        };
        let vec = &unsafe { &*self.weight_map.get() }[i];
        Voxel::Weight(vec)
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
        self.voxel_map[p as usize].load(Ordering::Relaxed)
    }

    /// A none locking retrieval of the state of voxel, p. This should only be
    /// used once the VoxelMap has been fully populated.
    pub fn voxel_get(&self, p: isize) -> Voxel {
        let volume_number =
            match self.voxel_map[p as usize].load(Ordering::Relaxed) {
                -1 => return Voxel::Vacuum,
                x => x as usize,
            };
        let i = match self.weight_index[p as usize].load(Ordering::Relaxed) {
            -1 => return Voxel::Maxima(volume_number),
            i => i as usize,
        };
        let vec = &unsafe { &*self.weight_map.get() }[i];
        Voxel::Weight(vec)
    }

    /// Finds the unique values contained in voxel_map and returns them sorted by
    /// value.
    pub fn maxima_list(&self) -> Vec<isize> {
        let mut maxima = BTreeSet::<isize>::new();
        self.voxel_map.iter().for_each(|x| {
                                 let bm = x.load(Ordering::Relaxed);
                                 maxima.insert(bm);
                             });
        maxima.into_iter().collect()
    }

    /// Stores the maxima of voxel, p, in the voxel_map.
    pub fn maxima_store(&self, p: isize, maxima: isize) {
        self.voxel_map[p as usize].store(maxima, Ordering::Relaxed);
    }

    /// Stores the index of p's weight contributions in weight_map into the
    /// weight_index and unlocks the structure.
    pub fn weight_store(&self, p: isize, i: usize) {
        self.weight_index[p as usize].store(i as isize, Ordering::Relaxed);
        self.lock.store(false, Ordering::SeqCst);
    }

    /// Locks the structure for write access.
    pub fn lock(&self) -> Lock {
        while self.lock.swap(true, Ordering::SeqCst) {}
        Lock { data: self }
    }

    /// Returns the index of p's maxima in the bader_maxima vector.
    fn index_get(&self, p: isize) -> Option<usize> {
        let maxima = self.maxima_non_block_get(p);
        if maxima.is_negative() {
            return None;
        }
        Some(self.maxima_index[maxima as usize])
    }

    /// Returns the atoms to which p is assigned.
    fn atom_get(&self, p: isize) -> Option<usize> {
        if let Some(i) = self.index_get(p) {
            return Some(self.assigned_atom[i]);
        }
        None
    }

    /// Checks to see if p is a boundary voxel between atoms.
    fn is_boundary(&self,
                   p: isize,
                   atom_number: usize,
                   density: &Density)
                   -> bool {
        for shift in density.voronoi.vectors.iter() {
            let pn = density.voronoi_shift(p, shift);
            match self.atom_get(pn) {
                Some(atom_num) => {
                    if atom_num != atom_number {
                        return true;
                    }
                }
                None => return true,
            }
        }
        false
    }

    /// Assigns each Bader maxima to an atom.
    pub fn assign_atoms(&mut self, atoms: &Atoms, density: &Density) {
        let mut index = vec![0usize; density.size.total];
        let bader_maxima = self.maxima_list();
        if bader_maxima[0] < 0 {
            for (i, maxima) in bader_maxima[1..].iter().enumerate() {
                index[*maxima as usize] = i;
            }
        } else {
            for (i, maxima) in bader_maxima.iter().enumerate() {
                index[*maxima as usize] = i;
            }
        }
        let (assigned_atom, minimum_distance) =
            atoms.assign_maxima(&bader_maxima, density);
        self.bader_maxima = bader_maxima;
        self.maxima_index = index;
        self.assigned_atom = assigned_atom;
        self.minimum_distance = minimum_distance;
    }

    /// Sums the densities for each bader volume.
    pub fn charge_sum(&mut self,
                      densities: &[Vec<f64>],
                      atoms: &Atoms,
                      density: &Density) {
        let mut bader_volume = vec![0.; self.bader_maxima.len()];
        let mut bader_charge =
            { vec![vec![0f64; self.bader_maxima.len()]; densities.len()] };
        let mut surface_distance = vec![0f64; atoms.positions.len()];
        let mut minimum_distance = {
            vec![
                atoms.lattice.distance_matrix[0].powi(2);
                atoms.positions.len()
            ]
        };
        let pbar = ProgressBar::new(density.size.total as u64);
        let pbar = Bar::new(pbar, 100, String::from("Summing Charge: "));
        'charge_sum: for p in (0..density.size.total).into_iter() {
            match self.voxel_get(p as isize) {
                Voxel::Weight(weights) => {
                    for (volume_number, weight) in weights.iter() {
                        if *volume_number == (-1isize) as usize {
                            println!("{:?}", weights);
                        }
                        let i = self.maxima_index[*volume_number];
                        bader_volume[i] += weight;
                        for j in 0..densities.len() {
                            bader_charge[j][i] += weight * densities[j][p];
                        }
                    }
                }
                Voxel::Maxima(volume_number) => {
                    let i = self.maxima_index[volume_number];
                    bader_volume[i] += 1.;
                    for j in 0..densities.len() {
                        bader_charge[j][i] += densities[j][p];
                    }
                }
                Voxel::Vacuum => {
                    let i = self.bader_maxima.len() - 1;
                    bader_volume[i] += 1.;
                    for j in 0..densities.len() {
                        bader_charge[j][i] += densities[j][p];
                    }
                    continue 'charge_sum;
                }
            }
            let atom_num = self.atom_get(p as isize).unwrap();
            if self.is_boundary(p as isize, atom_num, density) {
                let p_cartesian = density.to_cartesian(p as isize);
                let p_cartesian =
                    utils::dot(p_cartesian, density.voxel_lattice.to_cartesian);
                let mut p_lll_fractional = utils::dot(p_cartesian,
                                                      atoms.reduced_lattice
                                                           .to_fractional);
                for f in &mut p_lll_fractional {
                    *f = f.rem_euclid(1.);
                }
                let p_lll_cartesian = utils::dot(p_lll_fractional,
                                                 atoms.reduced_lattice
                                                      .to_cartesian);
                let atom = atoms.reduced_positions[atom_num];
                for atom_shift in
                    atoms.reduced_lattice.cartesian_shift_matrix.iter()
                {
                    let distance = {
                        (p_lll_cartesian[0] - (atom[0] + atom_shift[0])).powi(2)
                            + (p_lll_cartesian[1] - (atom[1] + atom_shift[1]))
                                .powi(2)
                            + (p_lll_cartesian[2] - (atom[2] + atom_shift[2]))
                                .powi(2)
                    };
                    if distance < minimum_distance[atom_num] {
                        surface_distance[atom_num] = distance.powf(0.5);
                        minimum_distance[atom_num] = distance;
                    }
                }
            }
            pbar.tick();
        }
        drop(pbar);
        for charge in &mut bader_charge {
            for charge in charge.iter_mut() {
                *charge *= density.voxel_lattice.volume;
            }
        }
        self.bader_charge = bader_charge;
        self.bader_volume =
            bader_volume.into_iter()
                        .map(|x| x * density.voxel_lattice.volume)
                        .collect();
        self.surface_distance = surface_distance;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rayon::prelude::*;
    use std::sync::Arc;

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
    fn voxel_map_weight_store() {
        let voxel_map = VoxelMap::new(10);
        let voxel_map = Arc::new(voxel_map);
        (0..10isize).into_par_iter().for_each(|p| {
                                        let i = {
                                            let mut weight = voxel_map.lock();
                                            (*weight).push(vec![(p as usize,
                                                                 p as f64
                                                                 - 1.)]);
                                            weight.len() - 1
                                        };
                                        voxel_map.weight_store(p, i);
                                        voxel_map.maxima_store(p, p - 1);
                                    });
        match voxel_map.weight_get(9) {
            Voxel::Weight(vec) => assert_eq!(vec[0], (9, 8.)),
            _ => panic!("not found weight in map"),
        }
        match voxel_map.voxel_get(9) {
            Voxel::Weight(vec) => assert_eq!(vec[0], (9, 8.)),
            _ => panic!("not found weight in map"),
        }
        assert!(matches!(voxel_map.voxel_get(0), Voxel::Vacuum));
    }
}
