use crate::grid::Grid;
use crate::progress::{Bar, HiddenBar, ProgressBar};
use crate::voxel_map::BlockingVoxelMap;
use crossbeam_utils::thread;
use rustc_hash::FxHashMap;
use std::sync::atomic::AtomicUsize;
use std::sync::Arc;

pub enum WeightResult {
    Maxima,
    Interior(usize),
    Boundary(Box<[f64]>),
    Saddle(Box<[f64]>),
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
/// * `weight_map`: An [`Arc`] wrapped [`VoxelMap`] for tracking the maxima.
///
/// ### Returns:
/// [`WeightResult`]: The type of point `p` is (maxima, interior, boundary) and
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
/// let weight = match weight_step(33, &density, &voxel_map, 1E-8) {
///     WeightResult::Saddle(weights) => weights,
///     _ => Vec::with_capacity(0).into(),
/// };
/// assert_eq!(weight, vec![61.375, 62.625].into())
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
    let mut saddle_flag = true;
    // colllect the shift and distances and iterate over them.
    for (pt, alpha) in grid.voronoi_shifts(p) {
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
                    saddle_flag = false;
                    let point_weights = voxel_map.weight_get(maxima);
                    for maxima_weight in point_weights.iter() {
                        let maxima = *maxima_weight as usize;
                        let w = maxima_weight - maxima as f64;
                        let weight = weights.entry(maxima).or_insert(0.);
                        *weight += w * rho;
                    }
                }
                // interior point
                std::cmp::Ordering::Greater => {
                    let weight = weights.entry(maxima as usize).or_insert(0.);
                    *weight += rho;
                }
                // going into vacuum (this is be impossible)
                std::cmp::Ordering::Equal => (),
            }
            t_sum += rho;
        }
    }
    match weights.len().cmp(&1) {
        // more than one weight is a boundary or saddle (if the weight is weighty enough)
        std::cmp::Ordering::Greater => {
            let mut total = 0.;
            // remove weuights below the tolerance
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
                // check if saddle point
                if saddle_flag {
                    WeightResult::Saddle(weights)
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
        std::cmp::Ordering::Less => WeightResult::Maxima,
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
) -> Vec<isize> {
    let counter = Arc::new(AtomicUsize::new(0));
    let mut saddle_points = vec![];
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
                    let mut saddles = vec![];
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
                            // Maxima should already be stored
                            WeightResult::Maxima => {
                                if voxel_map.maxima_check(p).is_none() {
                                    panic!(
                                        "Found new maxima in voxel assigning."
                                    )
                                }
                            }
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
                            WeightResult::Saddle(weights) => {
                                let i = {
                                    let mut weight = voxel_map.lock();
                                    let i = weight.len();
                                    (*weight).push(weights);
                                    i
                                };
                                voxel_map.weight_store(p, i);
                                saddles.push(p);
                            }
                        }
                        pbar.tick();
                    }
                    saddles
                })
            })
            .collect::<Vec<_>>();
        for thread in th {
            if let Ok(saddles) = thread.join() {
                saddle_points.extend(saddles);
            }
        }
    })
    .unwrap();
    {
        let mut weights = voxel_map.lock();
        weights.shrink_to_fit();
    }
    saddle_points.shrink_to_fit();
    saddle_points
}

/// Find the maxima within the charge density
pub fn maxima_finder(
    index: &[usize],
    density: &[f64],
    voxel_map: &BlockingVoxelMap,
    threads: usize,
    visible_bar: bool,
) -> Vec<isize> {
    let mut bader_maxima = Vec::<isize>::new();
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
                            for (pt, _) in
                                voxel_map.grid.voronoi_shifts(*p as isize)
                            {
                                if density[pt as usize] > rho {
                                    return None;
                                }
                            }
                            // if we made it this far we have a maxima
                            // change this index to a value it could
                            // never be and return it
                            Some(*p as isize)
                        })
                        .collect::<Vec<isize>>()
                })
            })
            .collect::<Vec<_>>();
        for thread in th {
            if let Ok(mut maxima_list) = thread.join() {
                bader_maxima.append(&mut maxima_list);
            } else {
                panic!("Failed to join thread in maxima finder.")
            };
        }
    })
    .unwrap(); // There is no panic option in the threads that isn't covered
    bader_maxima.shrink_to_fit();
    bader_maxima
}

/// Calculate the Laplacian of the density at a point in the grid
pub fn laplacian(p: usize, density: &[f64], grid: &Grid) -> f64 {
    let rho = density[p];
    grid.voronoi_shifts(p as isize)
        .iter()
        .fold(0.0, |acc, (pt, alpha)| {
            acc + alpha * (density[*pt as usize] - rho)
        })
        / grid.voronoi.volume
}
