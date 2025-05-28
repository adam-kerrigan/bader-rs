use crate::atoms::Atoms;
use crate::errors::MaximaError;
use crate::grid::Grid;
use crate::progress::{Bar, HiddenBar, ProgressBar};
use crate::voxel_map::{BlockingVoxelMap, VoxelMap};
use crossbeam_utils::thread;
use rustc_hash::FxHashMap;
use std::sync::atomic::AtomicUsize;
use std::sync::Arc;

/// Result of a Weight step.
///
/// TODO: turn this into an actual result type?
pub enum WeightResult {
    /// Length of the Box dictates the type of Critical Point, 1 -> Maxima, 2 -> Saddle,
    /// 3+ -> Saddle or minima. Critical Points with >=2 will be on boundaries.
    Critical(Box<[f64]>),
    /// Entirely assigned to a single Bader atom.
    Interior(usize),
    /// Meeting point at the edge of 2 or more Bader atoms.
    Boundary(Box<[f64]>),
    /// Maximum
    Maximum,
}

#[derive(Clone)]
pub struct CriticalPoint {
    pub position: isize,
    pub kind: CriticalPointKind,
    pub atoms: Box<[usize]>,
}

impl CriticalPoint {
    pub fn new(
        position: isize,
        kind: CriticalPointKind,
        atoms: Box<[usize]>,
    ) -> Self {
        CriticalPoint {
            position,
            kind,
            atoms,
        }
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd, Debug, Clone, Copy)]
pub enum CriticalPointKind {
    Nuclei,
    Bond,
    Ring,
    Cage,
    Blank,
}

/// Steps in the density grid, from point p, following the gradient.
///
/// This should be called from [`weight()`].
///
/// Note: This function will deadlock if the points above it have no associated
/// maxima in [`VoxelMap.voxel_map`].
///
/// * `p`: The point from which to step.
/// * `density`: The reference [`Grid`].
/// * `voxel_map`: An [`Arc`] wrapped [`BlockingVoxelMap`] for tracking the maxima.
/// * `weight_tolerance`: Minimum percentage value to consider the weight significant.
///
/// ### Returns:
/// [`WeightResult`]: The type of point `p` is Critical, Interior or Boundary and
/// the relevant data for each type.
///
/// # Examples
/// ```
/// use bader::methods::{weight_step, WeightResult};
/// use bader::voxel_map::BlockingVoxelMap as VoxelMap;
///
/// // Intialise the reference density, setting index 34 to 0. for easy maths.
/// let density = (0..64)
///     .map(|rho| if rho != 34 { rho as f64 } else { 0. })
///     .collect::<Vec<f64>>();
/// let voxel_map = VoxelMap::new(
///     [4, 4, 4],
///     [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]],
///     [0.0, 0.0, 0.0],
/// );
/// // The highest gradient between point, p = 33, and it's neighbours, with
/// // periodic boundary conditions, is with point p = 61.
///
/// // to avoid deadlock let's store maxima for all the values above us and
/// // store as either 61 or 62 to make the current point a boundary.
/// for (i, p) in [37, 45, 49].iter().enumerate() {
///     voxel_map.maxima_store(*p, 62 - (i as isize) % 2);
/// }
/// let weight: Vec<f64> = match weight_step(33, &density, &voxel_map, 1E-8) {
///     WeightResult::Critical(weights) => weights
///         .iter()
///         .map(|f| {
///             let maxima = *f as usize;
///             let weight = f - maxima as f64;
///             let (decoded_maxima, _) = voxel_map.grid.decode_maxima(maxima);
///             weight + decoded_maxima as f64
///         })
///         .collect(),
///     _ => panic!("None Weight"),
/// };
/// assert_eq!(weight, vec![61.375, 62.625])
/// ```
pub fn weight_step(
    p: isize,
    density: &[f64],
    voxel_map: &BlockingVoxelMap,
    weight_tolerance: f64,
) -> WeightResult {
    let control = density[p as usize];
    let grid = &voxel_map.grid;
    let mut t_sum = 0.;
    let mut weights = FxHashMap::<usize, f64>::default();
    let mut weight_count = 0;
    // colllect the shift and distances and iterate over them.
    grid.voronoi_shifts(p)
        .into_iter()
        .for_each(|((pt, image), alpha)| {
            let charge_diff = density[pt as usize] - control;
            // density differences of zero should be ignored to avoid division by
            // zero errors.
            if charge_diff > 0. {
                // calculate the gradient and add any weights to the HashMap.
                let rho = charge_diff * alpha;
                let maxima = voxel_map.maxima_get(pt);
                match maxima.cmp(&-1) {
                    // feeds into already weighted voxel therefore not a saddle point
                    std::cmp::Ordering::Less => {
                        let point_weights = voxel_map.weight_get(maxima);
                        weight_count = point_weights.len().max(weight_count);
                        for maxima_weight in point_weights.iter() {
                            let mut maxima = *maxima_weight as usize;
                            let w = maxima_weight - maxima as f64;
                            if image[0].abs() + image[1].abs() + image[2].abs()
                                != 0
                            {
                                maxima = grid.encode_maxima(maxima, image);
                            }
                            let weight = weights.entry(maxima).or_insert(0.);
                            *weight += w * rho;
                        }
                    }
                    // interior point
                    std::cmp::Ordering::Greater => {
                        let mut maxima = maxima as usize;
                        if image[0].abs() + image[1].abs() + image[2].abs() != 0
                        {
                            maxima = grid.encode_maxima(maxima, image);
                        }
                        let weight = weights.entry(maxima).or_insert(0.);
                        *weight += rho;
                    }
                    // going into vacuum (this be impossible)
                    std::cmp::Ordering::Equal => (),
                }
                t_sum += rho;
            }
        });
    match weights.len().cmp(&1) {
        // more than one weight is a boundary or saddle (if the weight is weighty enough)
        std::cmp::Ordering::Greater => {
            let mut total = 0.;
            // remove weights below the tolerance
            let weights = weights
                .into_iter()
                .filter_map(|(maxima, weight)| {
                    let weight = weight / t_sum;
                    if weight > weight_tolerance {
                        total += weight;
                        Some((maxima, weight))
                    } else {
                        None
                    }
                })
                .collect::<Vec<(usize, f64)>>();
            // still more than one weight then readjust the weights so that they sum to 1
            if let std::cmp::Ordering::Greater = weights.len().cmp(&1) {
                let weights = weights
                    .iter()
                    .map(|(maxima, w)| *maxima as f64 + w / total)
                    .collect::<Box<[f64]>>();
                // check if new maxima has joined the weights -> Critical Point (saddle/ring/cage)
                if weights.len() > weight_count {
                    WeightResult::Critical(weights)
                } else {
                    WeightResult::Boundary(weights)
                }
            } else {
                WeightResult::Interior(weights[0].0)
            }
        }
        // only feeds one atom means interior voxel
        std::cmp::Ordering::Equal => {
            WeightResult::Interior(*weights.keys().next().unwrap())
        }
        // no flux out means maximum
        std::cmp::Ordering::Less => WeightResult::Maximum,
    }
}

/// Assigns a maxima to the points within index.
///
/// Note: This function will deadlock if the points above it have no associated
/// maxima in [`VoxelMap.voxel_map`]. As such make sure index is sorted.
pub fn weight(
    density: &[f64],
    voxel_map: &BlockingVoxelMap,
    index: &[usize],
    weight_tolerance: f64,
    visible_bar: bool,
    threads: usize,
) -> (Vec<CriticalPoint>, Vec<CriticalPoint>) {
    let counter = Arc::new(AtomicUsize::new(0));
    let mut critical_points = (vec![], vec![]);
    let pbar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => {
            Box::new(Bar::new(index.len(), String::from("Bader Partitioning")))
        }
    };
    thread::scope(|s| {
        // Assign the remaining voxels to Bader maxima
        let th = (0..threads)
            .map(|_| {
                s.spawn(|_| {
                    let mut c_ps = (vec![], vec![]);
                    loop {
                        let p = {
                            let i = counter.fetch_add(
                                1,
                                std::sync::atomic::Ordering::Relaxed,
                            );
                            if i >= index.len() {
                                break;
                            };
                            index[i] as isize
                        };
                        match weight_step(
                            p,
                            density,
                            voxel_map,
                            weight_tolerance,
                        ) {
                            WeightResult::Maximum => {}
                            WeightResult::Interior(maxima) => {
                                voxel_map.maxima_store(p, maxima as isize);
                            }
                            WeightResult::Boundary(weights) => {
                                let i = {
                                    let mut weight = voxel_map.lock();
                                    let i = weight.len();
                                    (*weight).push(weights);
                                    i
                                };
                                voxel_map.weight_store(p, i);
                            }
                            WeightResult::Critical(weights) => {
                                // length = 1 is a maxima and doesn't need storing.
                                let (i, atoms) = {
                                    let mut weight = voxel_map.lock();
                                    let i = weight.len();
                                    let atoms: Vec<usize> = weights
                                        .iter()
                                        .map(|w| *w as usize)
                                        .collect();
                                    (*weight).push(weights);
                                    (i, atoms)
                                };
                                voxel_map.weight_store(p, i);
                                if atoms.len() < 3 {
                                    c_ps.0.push(CriticalPoint::new(
                                        p,
                                        CriticalPointKind::Bond,
                                        atoms.into(),
                                    ));
                                } else {
                                    c_ps.1.push(CriticalPoint::new(
                                        p,
                                        CriticalPointKind::Ring,
                                        atoms.into(),
                                    ));
                                }
                            }
                        }
                        pbar.tick();
                    }
                    c_ps
                })
            })
            .collect::<Vec<_>>();
        for thread in th {
            if let Ok(c_ps) = thread.join() {
                critical_points.0.extend(c_ps.0);
                critical_points.1.extend(c_ps.1);
            }
        }
    })
    .unwrap();
    {
        let mut weights = voxel_map.lock();
        weights.shrink_to_fit();
    }
    critical_points.0.shrink_to_fit();
    critical_points.1.shrink_to_fit();
    critical_points
}

/// Find the maxima within the charge density
pub fn maxima_finder(
    index: &[usize],
    density: &[f64],
    voxel_map: &BlockingVoxelMap,
    maximum_distance: &f64,
    atoms: &Atoms,
    threads: usize,
    visible_bar: bool,
) -> Result<Vec<CriticalPoint>, MaximaError> {
    let grid = &voxel_map.grid;
    let mut bader_maxima = Vec::<CriticalPoint>::new();
    let pbar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(index.len(), String::from("Maxima Finding"))),
    };
    let index_len = index.len();
    let chunk_size = (index_len / threads) + (index_len % threads).min(1);
    thread::scope(|s| {
        // Identify all the maxima
        let th = index
            .chunks(chunk_size)
            .map(|chunk| {
                s.spawn(|_| {
                    chunk
                        .iter()
                        .filter_map(|p| {
                            // we have to tick first due to early return
                            pbar.tick();
                            let rho = density[*p];
                            for (pt, _) in voxel_map
                                .grid
                                .voronoi_shifts_nocheck(*p as isize)
                            {
                                if density[pt as usize] > rho {
                                    return None;
                                }
                            }
                            // if we made it this far we have a maxima
                            // change this index to a value it could
                            // never be and return it
                            Some(
                                assign_maximum(
                                    *p as isize,
                                    atoms,
                                    grid,
                                    maximum_distance,
                                )
                                .map(|atom| {
                                    CriticalPoint::new(
                                        *p as isize,
                                        CriticalPointKind::Nuclei,
                                        Box::new([atom]),
                                    )
                                }),
                            )
                        })
                        .collect::<Result<Vec<CriticalPoint>, MaximaError>>()
                })
            })
            .collect::<Vec<_>>();
        for thread in th {
            if let Ok(maxima_list) = thread.join() {
                match maxima_list {
                    Ok(bm) => bader_maxima.extend(bm),
                    Err(e) => return Err(e),
                }
            } else {
                panic!("Failed to join thread in manima finder.")
            };
        }
        Ok(())
    })
    .unwrap()?; // There is no panic option in the threads that isn't covered
    bader_maxima.shrink_to_fit();
    Ok(bader_maxima)
}

/// Find minima in the charge density
pub fn minima_finder(
    index: &[usize],
    density: &[f64],
    voxel_map: &VoxelMap,
    threads: usize,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let mut bader_minima = Vec::<CriticalPoint>::new();
    let pbar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(index.len(), String::from("Minima Finding"))),
    };
    let index_len = index.len();
    let chunk_size = (index_len / threads) + (index_len % threads).min(1);
    thread::scope(|s| {
        // Identify all the maxima
        let th = index
            .chunks(chunk_size)
            .map(|chunk| {
                s.spawn(|_| {
                    chunk
                        .iter()
                        .filter_map(|p| {
                            // we have to tick first due to early return
                            pbar.tick();
                            let rho = density[*p];
                            for (pt, _) in voxel_map
                                .grid
                                .voronoi_shifts_nocheck(*p as isize)
                            {
                                if density[pt as usize] < rho {
                                    return None;
                                }
                            }
                            // if we made it this far we have a maxima
                            // change this index to a value it could
                            // never be and return it
                            // TODO: This needs to check if the cage is actually a boundary and if
                            // not complain that the weight tolerance is too high
                            Some(CriticalPoint::new(
                                *p as isize,
                                CriticalPointKind::Cage,
                                voxel_map
                                    .maxima_to_weight(
                                        voxel_map.maxima_get(*p as isize),
                                    )
                                    .iter()
                                    .map(|f| *f as usize)
                                    .collect(),
                            ))
                        })
                        .collect::<Vec<CriticalPoint>>()
                })
            })
            .collect::<Vec<_>>();
        for thread in th {
            if let Ok(minima_list) = thread.join() {
                bader_minima.extend(minima_list);
            } else {
                panic!("Failed to join thread in manima finder.")
            };
        }
    })
    .unwrap(); // There is no panic option in the threads that isn't covered
    bader_minima.shrink_to_fit();
    bader_minima
}

/// Assign the Bader maxima to the nearest atom.
///
/// # Example
/// ```
/// use bader::atoms::{Atoms, Lattice};
/// use bader::grid::Grid;
/// use bader::methods::assign_maximum;
///
/// // Intialise Atoms and Grid structs as well as a list of maxima
/// let lattice =
///     Lattice::new([[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]);
/// // Place atoms at 0 and 555 in the grid
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
///
/// // Run with default maxima distance tolerance
/// let maximum_distance = 0.1;
/// let atom_list = assign_maximum(555, &atoms, &grid, &maximum_distance);
/// assert!(atom_list.is_ok());
/// assert_eq!(atom_list.unwrap(), 1);
///
/// // If the maxima is too far away we get an error.
/// let atom_list = assign_maximum(554, &atoms, &grid, &maximum_distance);
/// assert!(atom_list.is_err());
/// ```
pub fn assign_maximum(
    maximum: isize,
    atoms: &Atoms,
    grid: &Grid,
    maximum_distance: &f64,
) -> Result<usize, MaximaError> {
    // convert the point first to cartesian, then to the reduced basis
    let m_cartesian = grid.to_cartesian(maximum);
    let m_reduced_cartesian = atoms.lattice.cartesian_to_reduced(m_cartesian);
    let mut atom_num = 0;
    let mut min_distance = f64::INFINITY;
    // go through each atom in the reduced basis and shift in each
    // reduced direction, save the atom with the shortest distance
    for (i, atom) in atoms.reduced_positions.iter().enumerate() {
        for atom_shift in atoms.lattice.reduced_cartesian_shift_matrix.iter() {
            let distance = {
                (m_reduced_cartesian[0] - (atom[0] + atom_shift[0])).powi(2)
                    + (m_reduced_cartesian[1] - (atom[1] + atom_shift[1]))
                        .powi(2)
                    + (m_reduced_cartesian[2] - (atom[2] + atom_shift[2]))
                        .powi(2)
            };
            if distance < min_distance {
                min_distance = distance;
                atom_num = i;
            }
        }
    }
    if min_distance.powf(0.5) > *maximum_distance {
        Err(MaximaError {
            maximum: m_cartesian,
            atom: atom_num,
            distance: min_distance.powf(0.5),
        })
    } else {
        Ok(atom_num)
    }
}

/// Calculate the Laplacian of the density at a point in the grid
pub fn laplacian(p: usize, density: &[f64], grid: &Grid) -> f64 {
    let rho = density[p];
    grid.voronoi_shifts_nocheck(p as isize)
        .iter()
        .fold(0.0, |acc, (pt, alpha)| {
            acc + alpha * (density[*pt as usize] - rho)
        })
        / grid.voronoi.volume
}
