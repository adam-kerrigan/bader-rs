use crate::progress::Bar;
use crate::voxel_map::BlockingVoxelMap as VoxelMap;
use anyhow::Result;
use atomic_counter::{AtomicCounter, RelaxedCounter};
use crossbeam_utils::thread;
use rustc_hash::FxHashMap;

pub enum WeightResult {
    Maxima,
    Interier(usize),
    Boundary(Box<[f64]>),
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
/// use bader::voxel_map::BlockingVoxelMap as VoxelMap;
/// use bader::methods::{WeightResult, weight_step};
///
/// // Intialise the reference density, setting index 34 to 0. for easy maths.
/// let density = (0..64).map(|rho| if rho != 34 { rho  as f64 } else {0.})
///                      .collect::<Vec<f64>>();
/// let voxel_map = VoxelMap::new([4, 4, 4],
///                               [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]],
///                               [0.0, 0.0, 0.0]);
/// // The highest gradient between point, p = 33, and it's neighbours, with
/// // periodic boundary conditions, is with point p = 61.
///
/// // to avoid deadlock let's store maxima for all the values above us and
/// // store as either 61 or 62 to make the current point a boundary.
/// for (i, p) in [37, 45, 49].iter().enumerate() {
///     voxel_map.maxima_store(*p, 62 - (i as isize) % 2);
/// }
/// let weight = match weight_step(33, &density, &voxel_map, 1E-8) {
///     WeightResult::Boundary(weights) => weights,
///     _ => Vec::with_capacity(0).into(),
/// };
/// assert_eq!(weight, vec![62.625, 61.375].into())
/// ```
pub fn weight_step(p: isize,
                   density: &[f64],
                   voxel_map: &VoxelMap,
                   weight_tolerance: f64)
                   -> WeightResult {
    let control = density[p as usize];
    let grid = &voxel_map.grid;
    let mut t_sum = 0.;
    let mut weights = FxHashMap::<usize, f64>::default();
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
                std::cmp::Ordering::Less => {
                    let point_weights = voxel_map.weight_get(maxima);
                    for maxima_weight in point_weights.iter() {
                        let maxima = *maxima_weight as usize;
                        let w = maxima_weight - maxima as f64;
                        let weight = weights.entry(maxima).or_insert(0.);
                        *weight += w * rho;
                    }
                }
                std::cmp::Ordering::Greater => {
                    let weight = weights.entry(maxima as usize).or_insert(0.);
                    *weight += rho;
                }
                std::cmp::Ordering::Equal => (),
            }
            t_sum += rho;
        }
    }
    // Sort the weights, if they exist, by the most probable.
    match weights.len().cmp(&1) {
        std::cmp::Ordering::Greater => {
            let mut total = 0.;
            let mut weights = weights.into_iter()
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
            if weights.len() > 1 {
                weights.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
                // re-adjust the weights
                let weights =
                    weights.iter()
                           .map(|(maxima, w)| *maxima as f64 + w / total)
                           .collect::<Box<[f64]>>();
                WeightResult::Boundary(weights)
            } else {
                WeightResult::Interier(weights[0].0)
            }
        }
        std::cmp::Ordering::Equal => {
            WeightResult::Interier(*weights.keys().next().unwrap())
        }
        std::cmp::Ordering::Less => WeightResult::Maxima,
    }
}

/// Assigns a maxima to the points within index.
///
/// Note: This function will deadlock if the points above it have no associated
/// maxima in [`VoxelMap.voxel_map`]. As such make sure index is sorted.
pub fn weight(density: &[f64],
              voxel_map: &VoxelMap,
              index: &[usize],
              progress_bar: Bar,
              threads: usize,
              weight_tolerance: f64) {
    let counter = RelaxedCounter::new(0);
    thread::scope(|s| {
        // Assign the remaining voxels to Bader maxima
        for _ in 0..threads {
            s.spawn(|_| loop {
                 let p = {
                     let i = counter.inc();
                     if i >= index.len() {
                         break;
                     };
                     index[i] as isize
                 };
                 match weight_step(p, density, voxel_map, weight_tolerance) {
                     WeightResult::Maxima => {
                         panic!("Found new maxima in voxel assigning.")
                     }
                     WeightResult::Interier(maxima) => {
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
                 }
                 progress_bar.tick();
             });
        }
    }).unwrap();
    {
        let mut weights = voxel_map.lock();
        weights.shrink_to_fit();
    }
}

/// Find (and remove from a sorted index) the maxima within the charge density
pub fn maxima_finder(index: &mut Vec<usize>,
                     density: &[f64],
                     voxel_map: &VoxelMap,
                     threads: usize,
                     progress_bar: Bar)
                     -> Result<Vec<isize>> {
    let mut bader_maxima = Vec::<isize>::new();
    let pbar = &progress_bar;
    let index_len = index.len();
    let chunk_size = (index_len / threads) + (index_len % threads).min(1);
    thread::scope(|s| {
        // Identify all the maxima
        let th = index.chunks_mut(chunk_size)
                      .map(|chunk| {
                          s.spawn(move |_| {
                               chunk.iter_mut()
                                    .filter_map(|p| {
                                        // we have to tick first due to early return
                                        pbar.tick();
                                        let rho = density[*p];
                                        for (pt, _) in
                                            voxel_map.grid
                                                     .voronoi_shifts(*p
                                                                     as isize)
                                        {
                                            if density[pt as usize] > rho {
                                                return None;
                                            }
                                        }
                                        // if we made it this far we have a maxima
                                        // change this index to a value it could
                                        // never be and return it
                                        let maxima = Some(*p as isize);
                                        *p = index_len;
                                        maxima
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
    }).unwrap();
    // remove the index of the maxima
    index.retain(|&i| i != index_len);
    Ok(bader_maxima)
}
