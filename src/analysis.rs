use crate::atoms::Atoms;
use crate::grid::Grid;
use crate::methods::{laplacian, CriticalPoint, CriticalPointKind};
use crate::progress::{Bar, HiddenBar, ProgressBar};
use crate::utils::{cross, dot, norm, subtract, vdot};
use crate::voxel_map::{Voxel, VoxelMap};
use crossbeam_utils::thread;
use rustc_hash::{FxHashMap, FxHashSet};

/// Sums the densities of each Bader volume.
///
/// #Example:
/// ```
/// use bader::analysis::calculate_bader_density;
/// use bader::atoms::{Atoms, Lattice};
/// use bader::grid::Grid;
/// use bader::voxel_map::VoxelMap;
///
/// // Intialise Atoms and VoxelMap structs as well as a density to sum.
/// let lattice =
///     Lattice::new([[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]);
/// let atoms = Atoms::new(
///     lattice,
///     vec![[0.0, 0.0, 0.0], [1.5, 1.5, 1.5]],
///     String::from(""),
/// );
/// let grid = Grid::new(
///     [10, 10, 10],
///     [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]],
///     [0.0, 0.0, 0.0],
/// );
/// // each atom gets 500 voxels all of value 1
/// let mut voxel_map = (0..1000).map(|i| i / 500).collect::<Vec<isize>>();
/// // add some vacuum meaning atom 2 has 499 voxels
/// voxel_map[600] = -1;
/// // add a weighted voxel meaning atom 1 now has 499.7 voxels and atom 2 has 499.3
/// voxel_map[400] = -2;
/// let weight_map: Vec<Box<[f64]>> = vec![vec![0.7, 1.3].into(); 1];
/// let voxel_map = VoxelMap::new(voxel_map, weight_map, grid);
/// let density = vec![1.0; 1000];
///
/// let summed_density =
///     calculate_bader_density(&density, &voxel_map, &atoms, 1, false);
/// let volume = voxel_map.grid_get().voxel_lattice.volume;
/// assert_eq!(
///     summed_density,
///     vec![
///         499.7 * volume,
///         499.3 * volume,
///         volume
///     ]
///     .into()
/// );
/// ```
pub fn calculate_bader_density(
    density: &[f64],
    voxel_map: &VoxelMap,
    atoms: &Atoms,
    threads: usize,
    visible_bar: bool,
) -> Box<[f64]> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            density.len(),
            String::from("Summing Bader Density"),
        )),
    };
    let pbar = &progress_bar;
    let mut bader_density = vec![0.0; atoms.positions.len() + 1];
    let vm = &voxel_map;
    // Calculate the size of the vector to be passed to each thread.
    let chunk_size =
        (density.len() / threads) + (density.len() % threads).min(1);
    thread::scope(|s| {
        let spawned_threads = voxel_map
            .maxima_chunks(chunk_size)
            .enumerate()
            .map(|(index, chunk)| {
                s.spawn(move |_| {
                    let mut bd = vec![0.0; atoms.positions.len() + 1];
                    chunk.iter().enumerate().for_each(
                        |(voxel_index, maxima)| {
                            let p = index * chunk.len() + voxel_index;
                            match vm.maxima_to_voxel(*maxima) {
                                Voxel::Maxima(m) => {
                                    let (m, _) =
                                        voxel_map.grid.decode_maxima(m);
                                    bd[m] += density[p];
                                }
                                Voxel::Boundary(weights) => {
                                    for weight in weights.iter() {
                                        let m = *weight as usize;
                                        let w = weight - (m as f64);
                                        let (m, _) =
                                            voxel_map.grid.decode_maxima(m);
                                        bd[m] += w * density[p];
                                    }
                                }
                                Voxel::Vacuum => {
                                    bd[atoms.positions.len()] += density[p]
                                }
                            };
                            pbar.tick();
                        },
                    );
                    bd
                })
            })
            .collect::<Vec<_>>();
        // Join each thread and collect the results.
        // If one thread terminates before the other this is not operated on first.
        // Either use the sorted index to remove vacuum from the summation or
        // find a way to operate on finshed threads first (ideally both).
        for thread in spawned_threads {
            if let Ok(tmp_bd) = thread.join() {
                bader_density.iter_mut().zip(tmp_bd.into_iter()).for_each(
                    |(a, b)| {
                        *a += b;
                    },
                );
            } else {
                panic!("Unable to join thread in sum_bader_densities.")
            };
        }
    })
    .unwrap();
    // The final result needs to be converted to a charge rather than a density.
    bader_density.iter_mut().for_each(|a| {
        *a *= voxel_map.grid_get().voxel_lattice.volume;
    });
    bader_density.into()
}

/// Calculates the volume and radius of each Bader atom.
///
/// #Example:
/// ```
/// use bader::analysis::calculate_bader_volumes_and_radii;
/// use bader::atoms::{Atoms, Lattice};
/// use bader::grid::Grid;
/// use bader::voxel_map::VoxelMap;
///
/// // Intialise Atoms and VoxelMap structs as well as a density to sum.
/// let lattice =
///     Lattice::new([[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]);
/// let atoms = Atoms::new(
///     lattice,
///     vec![[0.0, 0.0, 0.0], [1.5, 1.5, 1.5]],
///     String::from(""),
/// );
/// let grid = Grid::new(
///     [10, 10, 10],
///     [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]],
///     [0.0, 0.0, 0.0],
/// );
/// // each atom gets 500 voxels all of value 1
/// let mut voxel_map = (0..1000).map(|i| i / 500).collect::<Vec<isize>>();
/// // add some vacuum meaning atom 2 has 499 voxels
/// voxel_map[600] = -1;
/// // add a weighted voxel meaning atom 1 now has 499.7 voxels and atom 2 has 499.3
/// // this is the only factor in determining the radius; 2 * (a + b + c) for atom 1
/// // and 3 * (a + b + c) for atom 2.
/// voxel_map[222] = -2;
/// let weight_map: Vec<Box<[f64]>> = vec![vec![0.7, 1.3].into(); 1];
/// let voxel_map = VoxelMap::new(voxel_map, weight_map, grid);
///
/// let (volumes, radii) = calculate_bader_volumes_and_radii(&voxel_map, &atoms, 1, false);
/// let volume = voxel_map.grid_get().voxel_lattice.volume;
/// let a_b_c = (0.3_f64.powi(2) + 0.3_f64.powi(2) + 0.3_f64.powi(2)).powf(0.5);
/// assert_eq!(
///     volumes,
///     vec![
///         499.7 * volume,
///         499.3 * volume,
///         volume
///     ]
///     .into()
/// );
/// assert!(radii[0] - (2.0 * a_b_c) <= f64::EPSILON);
/// assert!(radii[1] - (3.0 * a_b_c) <= f64::EPSILON);
/// ```
pub fn calculate_bader_volumes_and_radii(
    voxel_map: &VoxelMap,
    atoms: &Atoms,
    threads: usize,
    visible_bar: bool,
) -> (Box<[f64]>, Box<[f64]>) {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            voxel_map.maxima_len(),
            String::from("Calculating Volumes"),
        )),
    };
    let pbar = &progress_bar;
    let mut bader_radius = vec![f64::INFINITY; atoms.positions.len()];
    let mut bader_volume = vec![0.0; atoms.positions.len() + 1];
    let vm = &voxel_map;
    // Calculate the size of the vector to be passed to each thread.
    let chunk_size = (voxel_map.maxima_len() / threads)
        + (voxel_map.maxima_len() % threads).min(1);
    thread::scope(|s| {
        let spawned_threads = voxel_map
            .maxima_chunks(chunk_size)
            .enumerate()
            .map(|(index, chunk)| {
                s.spawn(move |_| {
                    let mut br = vec![f64::INFINITY; atoms.positions.len()];
                    let mut bv = vec![0.0; atoms.positions.len() + 1];
                    chunk.iter().enumerate().for_each(|(voxel_index, maxima)| {
                        let p = index * chunk.len() + voxel_index;
                        match vm.maxima_to_voxel(*maxima) {
                            Voxel::Boundary(weights) => {
                                for weight in weights.iter() {
                                    let m = *weight as usize;
                                    let w = weight - (m as f64);
                                    let (m, _) =
                                        voxel_map.grid.decode_maxima(m);
                                    bv[m] += w;
                                    let atom_number = vm.maxima_to_atom(m);
                                    let p_c = vm.grid.to_cartesian(p as isize);
                                    let p_lll_c = atoms.lattice.cartesian_to_reduced(p_c);
                                    let atom = atoms.reduced_positions[atom_number];
                                    br[atom_number] = atoms.lattice.minimum_distance(p_lll_c, atom, Some(br[atom_number]));
                                }
                            }
                            Voxel::Maxima(m) => {
                                    let (m, _) =
                                        voxel_map.grid.decode_maxima(m);
                                bv[m] += 1.0;
                            }
                            Voxel::Vacuum => {
                                bv[atoms.positions.len()] += 1.0;
                            }
                        };
                        pbar.tick();
                    });
                    (bv, br)
                })
            })
            .collect::<Vec<_>>();
        // Join each thread and collect the results.
        // If one thread terminates before the other this is not operated on first.
        // Either use the sorted index to remove vacuum from the summation or
        // find a way to operate on finshed threads first (ideally both).
        for thread in spawned_threads {
            if let Ok((tmp_bv, tmp_br)) = thread.join() {
                bader_volume
                    .iter_mut()
                    .zip(tmp_bv.into_iter())
                    .for_each(|(a, b)| {
                        *a += b;
                    });
                bader_radius
                    .iter_mut()
                    .zip(tmp_br.into_iter())
                    .for_each(|(a, b)| {
                        *a = a.min(b);
                    });
            } else {
                panic!("Unable to join thread in calculate_bader_volumes_and_radii.")
            };
        }
    })
    .unwrap();
    // The distance isn't square rooted in the calcation of distance to save time.
    // As we need to filter out the infinite distances (atoms with no assigned maxima)
    // we can square root here also.
    bader_volume.iter_mut().for_each(|a| {
        *a *= voxel_map.grid_get().voxel_lattice.volume;
    });
    bader_radius.iter_mut().for_each(|d| {
        match (*d).partial_cmp(&f64::INFINITY) {
            Some(std::cmp::Ordering::Less) => *d = d.powf(0.5),
            _ => *d = 0.0,
        }
    });
    (bader_volume.into(), bader_radius.into())
}

/// calcuate the error associated with each atom from the Laplacian
pub fn calculate_bader_error(
    density: &[f64],
    voxel_map: &VoxelMap,
    atoms: &Atoms,
    threads: usize,
    visible_bar: bool,
) -> Box<[f64]> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            voxel_map.maxima_len(),
            String::from("Calculating Errors"),
        )),
    };
    let pbar = &progress_bar;
    let mut bader_error = vec![0.0; atoms.positions.len() + 1];
    let vm = &voxel_map;
    // Calculate the size of the vector to be passed to each thread.
    let chunk_size =
        (density.len() / threads) + (density.len() % threads).min(1);
    thread::scope(|s| {
        let spawned_threads = voxel_map
            .maxima_chunks(chunk_size)
            .enumerate()
            .map(|(index, chunk)| {
                s.spawn(move |_| {
                    let mut bd = vec![0.0; atoms.positions.len() + 1];
                    chunk.iter().enumerate().for_each(
                        |(voxel_index, maxima)| {
                            let p = index * chunk.len() + voxel_index;
                            let lapl = laplacian(p, density, &vm.grid);
                            match vm.maxima_to_voxel(*maxima) {
                                Voxel::Maxima(m) => {
                                    let (m, _) =
                                        voxel_map.grid.decode_maxima(m);
                                    bd[m] += lapl;
                                }
                                Voxel::Boundary(weights) => {
                                    for weight in weights.iter() {
                                        let m = *weight as usize;
                                        let w = weight - (m as f64);
                                        let (m, _) =
                                            voxel_map.grid.decode_maxima(m);
                                        bd[m] += w * lapl;
                                    }
                                }
                                Voxel::Vacuum => {
                                    bd[atoms.positions.len()] += lapl
                                }
                            };
                            pbar.tick();
                        },
                    );
                    bd
                })
            })
            .collect::<Vec<_>>();
        // Join each thread and collect the results.
        // If one thread terminates before the other this is not operated on first.
        // Either use the sorted index to remove vacuum from the summation or
        // find a way to operate on finshed threads first (ideally both).
        for thread in spawned_threads {
            if let Ok(tmp_bd) = thread.join() {
                bader_error.iter_mut().zip(tmp_bd.into_iter()).for_each(
                    |(a, b)| {
                        *a += b;
                    },
                );
            } else {
                panic!("Unable to join thread in sum_bader_densities.")
            };
        }
    })
    .unwrap();
    bader_error.into()
}

pub fn nuclei_ordering(
    nuclei: Vec<CriticalPoint>,
    density: &[f64],
    atom_len: usize,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            nuclei.len(),
            String::from("Pruning Nucleus Critical Points"),
        )),
    };
    let pbar = &progress_bar;
    // Find all the nuclei with the same atom number and group them based on which is largest
    // charge density. All maxima will be in image (0, 0, 0).
    // We need to order the nuclei by atom so that we can get the position of them by knowing the
    // atom number.
    let mut ordered_nuclei =
        vec![
            CriticalPoint::new(0, CriticalPointKind::Blank, Box::new([]));
            atom_len
        ];
    nuclei.iter().for_each(|cp| {
        let p = cp.position;
        let rho = density[p as usize];
        let atom_num = cp.atoms[0];
        if let CriticalPointKind::Blank = ordered_nuclei[atom_num].kind {
            ordered_nuclei[atom_num] =
                CriticalPoint::new(cp.position, cp.kind, cp.atoms.clone());
        } else if rho > density[ordered_nuclei[atom_num].position as usize] {
            ordered_nuclei[atom_num] =
                CriticalPoint::new(cp.position, cp.kind, cp.atoms.clone());
        }
        pbar.tick();
    });
    ordered_nuclei
}

pub fn bond_pruning(
    bonds: &[CriticalPoint],
    density: &[f64],
    grid: &Grid,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            bonds.len(),
            String::from("Pruning Bond Critical Points"),
        )),
    };
    let pbar = &progress_bar;
    bonds
        .iter()
        .filter_map(|cp| {
            pbar.tick();
            let mut origin_flag = false;
            cp.atoms.iter().for_each(|a| {
                let (_, image) = grid.decode_maxima(*a);
                if image[0].abs() + image[1].abs() + image[2].abs() == 0 {
                    origin_flag = true;
                }
            });
            let rho = density[cp.position as usize];
            let atom_num = match origin_flag {
                true => FxHashSet::from_iter(vec![cp.atoms.to_vec()]),
                false => cp
                    .atoms
                    .iter()
                    .map(|a| {
                        let (_, image) = grid.decode_maxima(*a);
                        cp.atoms
                            .iter()
                            .map(|a| {
                                let (a, i) = grid.decode_maxima(*a);
                                grid.encode_maxima(
                                    a,
                                    i.into_iter()
                                        .zip(image.iter())
                                        .map(|(i, ii)| i - *ii)
                                        .collect::<Vec<i8>>()
                                        .try_into()
                                        .unwrap(),
                                )
                            })
                            .collect::<Vec<usize>>()
                    })
                    .collect(),
            };
            for cp_t in bonds.iter() {
                let pt = cp_t.position;
                let mut origin_flag_t = false;
                cp_t.atoms.iter().for_each(|a| {
                    let (_, image) = grid.decode_maxima(*a);
                    if image[0].abs() + image[1].abs() + image[2].abs() == 0 {
                        origin_flag_t = true;
                    }
                });
                let atom_num_t = match origin_flag_t {
                    true => FxHashSet::from_iter(vec![cp_t.atoms.to_vec()]),
                    false => FxHashSet::from_iter(cp_t.atoms.iter().map(|a| {
                        let (_, image) = grid.decode_maxima(*a);
                        cp_t.atoms
                            .iter()
                            .map(|a| {
                                let (a, i) = grid.decode_maxima(*a);
                                grid.encode_maxima(
                                    a,
                                    i.into_iter()
                                        .zip(image.iter())
                                        .map(|(i, ii)| i - *ii)
                                        .collect::<Vec<i8>>()
                                        .try_into()
                                        .unwrap(),
                                )
                            })
                            .collect::<Vec<usize>>()
                    })),
                };
                for an in atom_num.iter() {
                    let an = FxHashSet::from_iter(an);
                    for ant in atom_num_t.iter() {
                        let ant = FxHashSet::from_iter(ant);
                        if an == ant && rho < density[pt as usize] {
                            return None;
                        }
                    }
                }
            }
            Some(CriticalPoint::new(cp.position, cp.kind, cp.atoms.clone()))
        })
        .collect()
}

pub fn ring_pruning(
    rings: &[CriticalPoint],
    ordered_nuclei: &[CriticalPoint],
    density: &[f64],
    atoms: &Atoms,
    grid: &Grid,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            rings.len(),
            String::from("Pruning Ring Critical Points"),
        )),
    };
    let pbar = &progress_bar;
    let mut pbc_rings = vec![];
    let rings = rings
        .iter()
        .filter_map(|cp| {
            let mut folded_atom_nums =
                Vec::<usize>::with_capacity(cp.atoms.len());
            let mut atom_images = Vec::<[i8; 3]>::with_capacity(cp.atoms.len());
            let mut origin_flag = false;
            let positions = cp.atoms[..3]
                .iter()
                .map(|a| {
                    let (atom_num, image) = grid.decode_maxima(*a);
                    folded_atom_nums.push(atom_num);
                    atom_images.push(image);
                    if image[0].abs() + image[1].abs() + image[2].abs() == 0 {
                        origin_flag = true;
                    }
                    let image_shift = dot(
                        image
                            .into_iter()
                            .map(|i| i as f64)
                            .collect::<Vec<f64>>()
                            .try_into()
                            .unwrap(),
                        atoms.lattice.to_cartesian,
                    );
                    grid.to_cartesian(ordered_nuclei[atom_num].position)
                        .iter()
                        .zip(image_shift)
                        .map(|(p, s)| *p + s)
                        .collect::<Vec<f64>>()
                        .try_into()
                        .unwrap()
                })
                .collect::<Vec<[f64; 3]>>();
            // if has 3 atoms they must form a plane
            if cp.atoms.len() > 3 {
                // create a vector normal to the plane span by the first 3 atoms
                let vec_1 = subtract(positions[1], positions[0]);
                let vec_2 = subtract(positions[2], positions[0]);
                let plane = cross(vec_1, vec_2);
                let plane_normal = norm(plane);
                let plane: [f64; 3] = plane
                    .into_iter()
                    .map(|f| f / plane_normal)
                    .collect::<Vec<f64>>()
                    .try_into()
                    .unwrap();
                // check if every other atom falls on that plane
                for a in cp.atoms[3..].iter() {
                    let (atom_num, image) = grid.decode_maxima(*a);
                    if image[0].abs() + image[1].abs() + image[2].abs() == 0 {
                        origin_flag = true;
                    }
                    let image_shift = dot(
                        image
                            .into_iter()
                            .map(|i| i as f64)
                            .collect::<Vec<f64>>()
                            .try_into()
                            .unwrap(),
                        atoms.lattice.to_cartesian,
                    );
                    let position = grid
                        .to_cartesian(ordered_nuclei[atom_num].position)
                        .iter()
                        .zip(image_shift.into_iter())
                        .map(|(p, s)| *p + s)
                        .collect::<Vec<f64>>()
                        .try_into()
                        .unwrap();
                    let vec_2 = subtract(position, positions[0]);
                    let plane_t = cross(vec_1, vec_2);
                    let plane_normal = norm(plane_t);
                    let plane_t: [f64; 3] = plane_t
                        .into_iter()
                        .map(|f| f / plane_normal)
                        .collect::<Vec<f64>>()
                        .try_into()
                        .unwrap();
                    // TODO: make this a tolerance currently 5.73 degrees
                    if vdot(plane, plane_t).abs() < 0.995 {
                        pbar.tick();
                        return None;
                    //TODO: remove this, its is for testing
                    } else {
                        folded_atom_nums.push(atom_num);
                        atom_images.push(image);
                    }
                }
            }
            // filter down all the same rings into the one with the highest charge density
            let p = cp.position;
            let rho = density[p as usize];
            let atom_nums = FxHashSet::from_iter(cp.atoms.iter());
            // check if the ring is the same as others and get rid if it has lower density
            for cp_t in rings.iter() {
                let pt = cp_t.position;
                let atom_num_t = FxHashSet::from_iter(cp_t.atoms.iter());
                if atom_nums == atom_num_t && rho < density[pt as usize] {
                    pbar.tick();
                    return None;
                }
            }
            // the origin flag is still false if none of the atoms are within the periodic bounds.
            // need to add these to a seperate list so they can be wrapped in and checked
            if !origin_flag {
                // these are all the unique collections of atoms that are shifted by the images
                let pbc_atoms =
                    FxHashSet::from_iter(cp.atoms.iter().map(|a| {
                        let (_, image) = grid.decode_maxima(*a);
                        cp.atoms
                            .iter()
                            .map(|a| {
                                let (pbc_a, pbc_image) = grid.decode_maxima(*a);
                                grid.encode_maxima(
                                    pbc_a,
                                    pbc_image
                                        .into_iter()
                                        .zip(image.iter())
                                        .map(|(pbc_i, i)| pbc_i - i)
                                        .collect::<Vec<i8>>()
                                        .try_into()
                                        .unwrap(),
                                )
                            })
                            .collect::<Vec<usize>>()
                    }));
                pbc_rings.push((
                    CriticalPoint::new(cp.position, cp.kind, cp.atoms.clone()),
                    pbc_atoms,
                ));
                None
            } else {
                Some(CriticalPoint::new(cp.position, cp.kind, cp.atoms.clone()))
            }
        })
        .collect::<Vec<CriticalPoint>>();
    // now we prune the out of origin atoms
    let mut pbc_rings = pbc_rings
        .iter()
        .enumerate()
        .filter_map(|(i, (cp, pbc_atoms))| {
            // we want to check every pbc ring against all of its pbc reflections
            let p = cp.position;
            let rho = density[p as usize];
            // first we compare against all of the other non origin rings and remove if it has
            // lower density
            for (ii, (cp_t, pbc_t)) in pbc_rings.iter().enumerate() {
                if i != ii {
                    let pt = cp_t.position;
                    for atom_nums in pbc_atoms.iter() {
                        let atom_nums = FxHashSet::from_iter(atom_nums);
                        for atom_nums_t in pbc_t.iter() {
                            let atom_nums_t = FxHashSet::from_iter(atom_nums_t);
                            if atom_nums == atom_nums_t {
                                if rho < density[pt as usize] {
                                    pbar.tick();
                                    return None;
                                }
                            } else if atom_nums.is_subset(&atom_nums_t) {
                                pbar.tick();
                                return None;
                            }
                        }
                    }
                }
            }
            // then compare against the other rings. this is mor complicated, if they are the same
            // remove if the density is lower, else we need to add this one and later remove the
            // one that has lower density
            for cp_t in rings.iter() {
                let pt = cp_t.position;
                for atom_nums in pbc_atoms.iter() {
                    let atom_nums = FxHashSet::from_iter(atom_nums);
                    let atom_nums_t = FxHashSet::from_iter(cp_t.atoms.iter());
                    // if they have the exact same atoms check the density
                    if atom_nums == atom_nums_t {
                        if rho < density[pt as usize] {
                            pbar.tick();
                            return None;
                        } else {
                            return Some(CriticalPoint::new(
                                cp.position,
                                cp.kind,
                                cp_t.atoms.clone(),
                            ));
                        }
                    // if its a subset remove it
                    } else if atom_nums.is_subset(&atom_nums_t) {
                        pbar.tick();
                        return None;
                    }
                }
            }
            Some(CriticalPoint::new(cp.position, cp.kind, cp.atoms.clone()))
        })
        .collect::<Vec<CriticalPoint>>();
    pbc_rings.extend(rings);
    pbc_rings
        .iter()
        .filter_map(|cp| {
            pbar.tick();
            let atom_num = FxHashSet::from_iter(cp.atoms.iter());
            let rho = density[cp.position as usize];
            // filter down all the same rings into the one with the highest charge density
            for cp_t in pbc_rings.iter() {
                let pt = cp_t.position;
                let atom_num_t = FxHashSet::from_iter(cp_t.atoms.iter());
                if atom_num == atom_num_t {
                    if rho < density[pt as usize] {
                        return None;
                    }
                } else if atom_num.is_subset(&atom_num_t) {
                    return None;
                }
            }
            Some(CriticalPoint::new(cp.position, cp.kind, cp.atoms.clone()))
        })
        .collect()
}

pub fn cage_pruning(
    cages: &[CriticalPoint],
    ordered_nuclei: &[CriticalPoint],
    density: &[f64],
    atoms: &Atoms,
    grid: &Grid,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            cages.len(),
            String::from("Pruning Cage Critical Points"),
        )),
    };
    let pbar = &progress_bar;
    cages
        .iter()
        .filter_map(|cp| {
            pbar.tick();
            if cp.atoms.len() < 4 {
                // impossible to have a cage with 3 atoms
                return None;
            } else {
                // form a plane with the first 3 atoms
                let positions = cp.atoms[..3]
                    .iter()
                    .map(|a| {
                        let (atom_num, image) = grid.decode_maxima(*a);
                        let image_shift = dot(
                            image
                                .into_iter()
                                .map(|i| i as f64)
                                .collect::<Vec<f64>>()
                                .try_into()
                                .unwrap(),
                            atoms.lattice.to_cartesian,
                        );
                        // get the position of the nuclei and apply the image shift
                        grid.to_cartesian(ordered_nuclei[atom_num].position)
                            .iter()
                            .zip(image_shift)
                            .map(|(p, s)| *p + s)
                            .collect::<Vec<f64>>()
                            .try_into()
                            .unwrap()
                    })
                    .collect::<Vec<[f64; 3]>>();
                let vec_1 = subtract(positions[1], positions[0]);
                let vec_2 = subtract(positions[2], positions[0]);
                let plane = cross(vec_1, vec_2);
                let plane_normal = norm(plane);
                let plane: [f64; 3] = plane
                    .into_iter()
                    .map(|f| f / plane_normal)
                    .collect::<Vec<f64>>()
                    .try_into()
                    .unwrap();
                let mut flag = true;
                // check every other atom against this plane, at least one should not be in the
                // plane
                for a in cp.atoms[3..].iter() {
                    let (atom_num, image) = grid.decode_maxima(*a);
                    let image_shift = dot(
                        image
                            .into_iter()
                            .map(|i| i as f64)
                            .collect::<Vec<f64>>()
                            .try_into()
                            .unwrap(),
                        atoms.lattice.to_cartesian,
                    );
                    let position = grid
                        .to_cartesian(ordered_nuclei[atom_num].position)
                        .iter()
                        .zip(image_shift.into_iter())
                        .map(|(p, s)| *p + s)
                        .collect::<Vec<f64>>()
                        .try_into()
                        .unwrap();
                    let vec_2 = subtract(position, positions[0]);
                    let plane_t = cross(vec_1, vec_2);
                    let plane_normal = norm(plane_t);
                    let plane_t: [f64; 3] = plane_t
                        .into_iter()
                        .map(|f| f / plane_normal)
                        .collect::<Vec<f64>>()
                        .try_into()
                        .unwrap();
                    // TODO: make this a tolerance currently 5.73 degrees
                    if vdot(plane, plane_t).abs() < 0.995 {
                        flag = false;
                    }
                }
                // if the flag is still true then all were in the same plane and the cage is a
                // density fluctuation
                if flag {
                    return None;
                }
            }
            let p = cp.position;
            let rho = density[p as usize];
            let atom_num = FxHashSet::from_iter(cp.atoms.iter());
            // now we check against the other cages, note there is no check on folded atoms
            for cp_t in cages {
                let pt = cp_t.position;
                let atom_num_t = FxHashSet::from_iter(cp_t.atoms.iter());
                if atom_num_t == atom_num && rho > density[pt as usize] {
                    return None;
                }
            }
            Some(CriticalPoint::new(cp.position, cp.kind, cp.atoms.clone()))
        })
        .collect()
}

/// Calculates the Laplacian at each saddle point. This is currently basic analysis, atoms images
/// are associated by distance not gradient paths and ring points are just being ignored.
pub fn calculate_bond_strengths(
    saddles: &[isize],
    density: &[f64],
    atoms: &Atoms,
    voxel_map: &VoxelMap,
    visible_bar: bool,
) -> Vec<FxHashMap<(usize, usize), f64>> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            voxel_map.maxima_len(),
            String::from("Calculating Bond Strength"),
        )),
    };
    let pbar = &progress_bar;
    let lat = &atoms.lattice;
    let mut bonds = vec![
        FxHashMap::<(usize, usize), f64>::default();
        atoms.positions.len()
    ];
    saddles.iter().for_each(|p| {
        let lapl = laplacian(*p as usize, density, voxel_map.grid_get());
        let p_lll = lat.cartesian_to_reduced(voxel_map.grid.to_cartesian(*p));
        // We know that the saddle point is a weight not a maximum.
        let weights = voxel_map.maxima_to_weight(voxel_map.maxima_get(*p));
        // If the weight is two atoms it's a bond > 2 is a ring?
        if let std::cmp::Ordering::Equal = weights.len().cmp(&2) {
            let atom_nums = weights
                .iter()
                .map(|w| {
                    let n = *w as usize;
                    let atom = atoms.reduced_positions[n];
                    let mut min_distance = f64::INFINITY;
                    let mut atom_image = 0;
                    for (i, atom_shift) in
                        lat.reduced_cartesian_shift_matrix.iter().enumerate()
                    {
                        let distance = {
                            (p_lll[0] - (atom[0] + atom_shift[0])).powi(2)
                                + (p_lll[1] - (atom[1] + atom_shift[1])).powi(2)
                                + (p_lll[2] - (atom[2] + atom_shift[2])).powi(2)
                        };
                        if distance < min_distance {
                            min_distance = distance;
                            atom_image = i;
                        }
                    }
                    (atom_image, n)
                })
                .collect::<Vec<(usize, usize)>>();
            // compare images
            // crossing (out[0] * 9 + out[1] * 3 + out[2] + 13) as usize
            let (image_1, atom_1) = atom_nums[0];
            let (image_2, atom_2) = atom_nums[1];
            let x1 = (image_1 / 9) as isize - 1;
            let y1 = (image_1 / 3).rem_euclid(3) as isize - 1;
            let z1 = image_1.rem_euclid(3) as isize - 1;
            let x2 = (image_2 / 9) as isize - 1;
            let y2 = (image_2 / 3).rem_euclid(3) as isize - 1;
            let z2 = image_2.rem_euclid(3) as isize - 1;
            let image_1_adjust =
                ((x1 - x2) * 9 + (y1 - y2) * 3 + (z1 - z2) + 13) as usize;
            let image_2_adjust =
                ((x2 - x1) * 9 + (y2 - y1) * 3 + (z2 - z1) + 13) as usize;
            let bonds_1 =
                bonds[atom_1].entry((atom_2, image_2_adjust)).or_insert(0.0);
            if let std::cmp::Ordering::Greater =
                lapl.abs().partial_cmp(&bonds_1.abs()).unwrap()
            {
                *bonds_1 = lapl;
                let bonds_2 = bonds[atom_2]
                    .entry((atom_1, image_1_adjust))
                    .or_insert(0.0);
                *bonds_2 = lapl;
            };
        }
        pbar.tick();
    });
    bonds.iter_mut().for_each(|b| b.shrink_to_fit());
    bonds
}
