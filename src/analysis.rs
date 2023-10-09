use crate::atoms::Atoms;
use crate::errors::MaximaError;
use crate::grid::Grid;
use crate::methods::laplacian;
use crate::progress::{Bar, HiddenBar, ProgressBar};
use crate::voxel_map::{Voxel, VoxelMap};
use crossbeam_utils::thread;
use rustc_hash::FxHashMap;

/// Assign the Bader maxima to the nearest atom.
///
/// # Example
/// ```
/// use bader::analysis::assign_maxima;
/// use bader::atoms::{Atoms, Lattice};
/// use bader::grid::Grid;
///
/// // Intialise [`Atoms`] and [`Grid`] structs as well as a list of maxima
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
/// let maxima = vec![0, 555]; // maxima placed at the sites of the atoms
///
/// // Run with default maxima distance tolerance
/// let maximum_distance = 0.1;
/// let atom_list =
///     assign_maxima(&maxima, &atoms, &grid, &maximum_distance, 1, false);
/// assert!(atom_list.is_ok());
/// assert_eq!(atom_list.unwrap(), vec![0, 1]);
///
/// // If the maxima is too far away we get an error.
/// let maxima = vec![1, 555]; // maxima placed 0.3 Ang away from atom 0
/// let atom_list =
///     assign_maxima(&maxima, &atoms, &grid, &maximum_distance, 1, false);
/// assert!(atom_list.is_err());
/// ```
pub fn assign_maxima(
    maxima: &[isize],
    atoms: &Atoms,
    grid: &Grid,
    maximum_distance: &f64,
    threads: usize,
    visible_bar: bool,
) -> Result<Vec<usize>, MaximaError> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => {
            Box::new(Bar::new(maxima.len(), String::from("Assigning to Atoms")))
        }
    };
    let pbar = &progress_bar;
    let chunk_size = (maxima.len() / threads) + (maxima.len() % threads).min(1);
    thread::scope(|s| {
        let mut assigned_atom = Vec::with_capacity(maxima.len());
        let spawned_threads = maxima
            .chunks(chunk_size)
            .map(|chunk| {
                s.spawn(move |_| {
                    let chunk_size = chunk.len();
                    // create vectors for storing the assigned atom and distance for each maxima
                    let mut ass_atom = Vec::with_capacity(chunk_size);
                    for m in chunk.iter() {
                        // convert the point first to cartesian, then to the reduced basis
                        let m_cartesian = grid.to_cartesian(*m);
                        let m_reduced_cartesian =
                            atoms.lattice.cartesian_to_reduced(m_cartesian);
                        let mut atom_num = 0;
                        let mut min_distance = f64::INFINITY;
                        // go through each atom in the reduced basis and shift in each
                        // reduced direction, save the atom with the shortest distance
                        for (i, atom) in
                            atoms.reduced_positions.iter().enumerate()
                        {
                            for atom_shift in atoms
                                .lattice
                                .reduced_cartesian_shift_matrix
                                .iter()
                            {
                                let distance = {
                                    (m_reduced_cartesian[0]
                                        - (atom[0] + atom_shift[0]))
                                        .powi(2)
                                        + (m_reduced_cartesian[1]
                                            - (atom[1] + atom_shift[1]))
                                            .powi(2)
                                        + (m_reduced_cartesian[2]
                                            - (atom[2] + atom_shift[2]))
                                            .powi(2)
                                };
                                if distance < min_distance {
                                    min_distance = distance;
                                    atom_num = i;
                                }
                            }
                        }
                        if min_distance.powf(0.5) > *maximum_distance {
                            return Err(MaximaError {
                                maximum: m_cartesian,
                                atom: atom_num,
                                distance: min_distance.powf(0.5),
                            });
                        }
                        ass_atom.push(atom_num);
                        pbar.tick();
                    }
                    Ok(ass_atom)
                })
            })
            .collect::<Vec<_>>();
        // we go through in the order they are spawned so order is retained.
        for thread in spawned_threads {
            if let Ok(result) = thread.join() {
                let ass_atom = match result {
                    Ok(v) => v,
                    Err(e) => return Err(e),
                };
                assigned_atom.extend(ass_atom.into_iter());
            } else {
                panic!("Failed to join thread in assign maxima.")
            };
        }
        Ok(assigned_atom)
    })
    .unwrap()
}

/// Sums the densities of each Bader volume.
///
/// #Example:
/// ```
/// use bader::analysis::calculate_bader_density;
/// use bader::atoms::{Atoms, Lattice};
/// use bader::grid::Grid;
/// use bader::voxel_map::VoxelMap;
///
/// // Intialise [`Atoms`] and [`VoxelMap`] structs as well as a density to sum.
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
                                    bd[m] += density[p];
                                }
                                Voxel::Boundary(weights) => {
                                    for weight in weights.iter() {
                                        let m = *weight as usize;
                                        let w = weight - (m as f64);
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
/// // Intialise [`Atoms`] and [`VoxelMap`] structs as well as a density to sum.
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
                                    bv[m] += w;
                                    let atom_number = vm.maxima_to_atom(m);
                                    let mr = &mut br[atom_number];
                                    let p_c = vm.grid_get().to_cartesian(p as isize);
                                    let p_lll_c = atoms.lattice.cartesian_to_reduced(p_c);
                                    let atom = atoms.reduced_positions[atom_number];
                                    for atom_shift in
                                        atoms.lattice.reduced_cartesian_shift_matrix.iter()
                                    {
                                        let distance = {
                                            (p_lll_c[0] - (atom[0] + atom_shift[0])).powi(2)
                                                + (p_lll_c[1] - (atom[1] + atom_shift[1])).powi(2)
                                                + (p_lll_c[2] - (atom[2] + atom_shift[2])).powi(2)
                                        };
                                        if distance < *mr {
                                            *mr = distance;
                                        }
                                    }
                                }
                            }
                            Voxel::Maxima(m) => {
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
                                    bd[m] += lapl;
                                }
                                Voxel::Boundary(weights) => {
                                    for weight in weights.iter() {
                                        let m = *weight as usize;
                                        let w = weight - (m as f64);
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
