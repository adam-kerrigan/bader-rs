use crate::atoms::Atoms;
use crate::critical::{CriticalPoint, CriticalPointKind};
use crate::errors::MaximaError;
use crate::grid::Grid;
use crate::progress::{Bar, HiddenBar, ProgressBar};
use crate::voxel_map::{
    BlockingVoxelMap, EncodedAtom, EncodedImage, EncodedWeight, Voxel, VoxelMap,
};
use rustc_hash::FxHashMap;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::thread;

/// Result of a Weight step.
///
/// TODO: turn this into an actual result type?
pub enum WeightResult {
    /// Length of the Box dictates the type of Critical Point, 1 -> Maxima, 2 -> Saddle,
    /// 3+ -> Saddle or minima. Critical Points with >=2 will be on boundaries.
    Critical(Box<[EncodedWeight]>),
    /// Entirely assigned to a single Bader atom.
    Interior(usize),
    /// Meeting point at the edge of 2 or more Bader atoms.
    Boundary(Box<[EncodedWeight]>),
    /// Maximum
    Maximum,
}

/// Performs a single step of gradient ascent from a voxel to determine its owner.
///
/// This function looks at the current voxel `p` and its neighbors. It calculates the flux
/// of charge density into neighboring voxels that have a higher density (following the gradient).
///
/// # Logic
/// * **Maximum**: If no neighbors have higher density, `p` is a local maximum.
/// * **Interior**: If `p` flows entirely into neighbors belonging to the *same* Atom/Maxima,
///   then `p` also belongs to that Atom.
/// * **Boundary**: If `p` flows into neighbors belonging to *different* Atoms, `p` is on a boundary.
///   Returns a weighted list of contributions.
///
/// # Arguments
/// * `p`: The index of the voxel to step from.
/// * `density`: The charge density array.
/// * `voxel_map`: The map storing the state (Maxima/Boundary) of processed voxels.
/// * `weight_tolerance`: Minimum weight fraction (0.0-1.0) to be considered significant.
pub fn weight_step(
    p: isize,
    density: &[f64],
    voxel_map: &BlockingVoxelMap,
    weight_tolerance: f64,
) -> WeightResult {
    let control = density[p as usize];
    let grid = &voxel_map.grid;
    let mut t_sum = 0.;
    let mut weights = FxHashMap::<EncodedAtom, f64>::default();
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
                match voxel_map.voxel_get(pt) {
                    // feeds into already weighted voxel therefore not a saddle point
                    Voxel::Boundary(weight_map) => {
                        weight_count = weight_map.len().max(weight_count);
                        weight_map.into_iter().for_each(|(maxima, weight)| {
                            let maxima = match image.is_zero() {
                                true => maxima,
                                false => maxima.image_add(image),
                            };
                            let w = weights.entry(maxima).or_insert(0.);
                            *w += weight as f64 * rho
                        });
                    },
                    // interior point
                    Voxel::Maxima(maxima) => {
                        let maxima = match image.is_zero() {
                            true => maxima,
                            false => maxima.image_add(image),
                        };
                        let w = weights.entry(maxima).or_insert(0.);
                        *w += rho
                    }
                    // going into vacuum (this be impossible)
                    Voxel::Vacuum => panic!("Vacuum voxel found with higher charge density than the control voxel.")
                };
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
                        Some((maxima, weight as f32))
                    } else {
                        None
                    }
                })
                .collect::<Vec<(EncodedAtom, f32)>>();
            // still more than one weight then readjust the weights so that they sum to 1
            if let std::cmp::Ordering::Greater = weights.len().cmp(&1) {
                let weights = weights
                    .into_iter()
                    .map(|(maxima, w)| EncodedWeight::new(maxima, w))
                    .collect::<Box<[EncodedWeight]>>();
                // check if new maxima has joined the weights -> Critical Point (saddle/ring/cage)
                if weights.len() > weight_count {
                    WeightResult::Critical(weights)
                } else {
                    WeightResult::Boundary(weights)
                }
            } else {
                WeightResult::Interior(weights[0].0.0 as usize)
            }
        }
        // only feeds one atom means interior voxel
        std::cmp::Ordering::Equal => {
            WeightResult::Interior(weights.keys().next().unwrap().0 as usize)
        }
        // no flux out means maximum
        std::cmp::Ordering::Less => WeightResult::Maximum,
    }
}

/// Assigns every voxel in the sorted `index` list to a Bader volume (atom) or boundary.
///
/// This function iterates through the density grid (following the order in `index`) and uses
/// gradient ascent to determine which atom(s) each voxel belongs to. It populates the
/// `voxel_map` in-place and identifies Critical Points (saddles) found during the process.
///
/// # Thread Safety & Deadlocks
/// This function spawns multiple threads to process the voxels.
/// **Crucial**: The `index` slice **must** be sorted in descending order of charge density.
/// The algorithm relies on the fact that when processing a voxel `p`, all voxels `q` with
/// $\rho(q) > \rho(p)$ (i.e., further up the gradient path) have already been processed and
/// assigned. Violating this order will result in deadlocks or incorrect assignments.
///
/// # Arguments
/// * `density`: The full flattened charge density array.
/// * `voxel_map`: The thread-safe map where results (Maxima/Weights) will be stored.
/// * `index`: A list of voxel indices sorted by density (Highest $\to$ Lowest).
/// * `weight_tolerance`: Minimum contribution required for a boundary weight (typically `1e-6` to `1e-3`).
/// * `visible_bar`: If `true`, displays a progress bar to `stdout`.
/// * `threads`: The number of worker threads to spawn.
///
/// # Returns
/// A tuple containing two vectors of [`CriticalPoint`]:
/// 1. **Bond Points**: Critical points connecting exactly 2 atoms.
/// 2. **Ring/Cage Points**: Critical points connecting 3 or more atoms.
///
/// # Example
/// ```
/// use bader::methods::weight;
/// use bader::voxel_map::{BlockingVoxelMap, VoxelMap, EncodedAtom};
///
/// // 1. Setup Data
/// let density = vec![2.0, 2.0, 12.0, 2.0, 2.0, 11.0, 1.0, 2.0, 6.0, 1.0, 1.0, 5.0]; // 1 atom at voxel 2
/// // Sorted indices: 3 (10.0) -> 2 (2.0) -> 1 (1.0) -> 0 (0.0)
/// let sorted_indices = vec![2, 5, 8, 11, 0, 1, 3, 4, 7, 6, 9, 10];
///
/// // 2. Setup Map
/// let map = BlockingVoxelMap::new(
///     [2, 2, 3],
///     [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
///     [0.0, 0.0, 0.0]
/// );
/// // Pre-assign the maxima (index 2) to Atom 0
/// map.maxima_store(2, EncodedAtom::new_zero_image(0).0 as isize);
///
/// // 3. Run Partitioning
/// // processing indices (skipping 2 as it's already done)
/// let (bonds, rings) = weight(
///     &density,
///     &map,
///     &sorted_indices[1..], // Skip the maxima itself
///     1e-6,
///     false, // No progress bar
///     1      // Single thread
/// );
///
/// // 4. Check Results
/// // Voxel 1 should have flowed up to Voxel 3 (Atom 0)
/// let map = VoxelMap::from_blocking_voxel_map(map);
/// use bader::voxel_map::Voxel;
/// if let Voxel::Maxima(atom) = map.voxel_get(1) {
///     assert_eq!(atom.atom_index(), 0);
/// } else {
///     panic!("Voxel 1 should be assigned to Atom 1");
/// }
/// ```
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
                s.spawn(|| {
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
                                voxel_map.weight_store(p, weights);
                            }
                            WeightResult::Critical(weights) => {
                                // length = 1 is a maxima and doesn't need storing.
                                let atoms: Vec<EncodedAtom> = weights
                                    .iter()
                                    .map(|ed| ed.decode().0)
                                    .collect();
                                voxel_map.weight_store(p, weights);
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
    });
    critical_points.0.shrink_to_fit();
    critical_points.1.shrink_to_fit();
    critical_points
}

/// Scans the grid to identify all local maxima in the charge density.
///
/// This is the first step of the Bader analysis. It iterates over all voxels (in parallel)
/// and checks if a voxel is strictly greater than all its neighbors.
///
/// Valid maxima are then assigned to atoms using [`assign_maximum`].
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
                s.spawn(|| {
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
                                        Box::new([EncodedAtom::new(
                                            atom as u32,
                                            EncodedImage::new([0; 3]),
                                        )]),
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
    })?;
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
                s.spawn(|| {
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
                            if let Voxel::Boundary(weights) =
                                voxel_map.voxel_get(*p as isize)
                            {
                                return Some(CriticalPoint::new(
                                    *p as isize,
                                    CriticalPointKind::Cage,
                                    weights.into_keys().collect(),
                                ));
                            }
                            None
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
    });
    bader_minima.shrink_to_fit();
    bader_minima
}

/// Associates a grid point (Maxima) with the nearest Atom.
///
/// This function calculates the Euclidean distance from the grid point `maximum` to all atoms
/// in the system, respecting Periodic Boundary Conditions (PBC).
///
/// # Arguments
/// * `maximum`: The 1D index of the grid point.
/// * `atoms`: The struct containing atomic positions.
/// * `grid`: The grid definition (used for coordinate conversion).
/// * `maximum_distance`: The cutoff distance. If the nearest atom is further than this, an error is returned.
///
/// # Returns
/// * `Ok(usize)`: The index of the assigned atom.
/// * `Err(MaximaError)`: If no atom is found within `maximum_distance`.
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

/// Calculates the discrete Laplacian of the density at a specific voxel.
///
/// The Laplacian ($\nabla^2 \rho$) measures the curvature of the density.
/// * **Negative**: The density is at a local peak (charge concentration).
/// * **Positive**: The density is at a local minimum (charge depletion).
///
/// This implementation uses a Voronoi-weighted finite difference method compatible with
/// the grid's neighbor connectivity.
///
/// # Arguments
/// * `p`: The voxel index.
/// * `density`: The charge density array.
/// * `grid`: The grid containing neighbor and Voronoi weight information.
pub fn laplacian(p: usize, density: &[f64], grid: &Grid) -> f64 {
    let rho = density[p];
    grid.voronoi_shifts_nocheck(p as isize)
        .iter()
        .fold(0.0, |acc, (pt, alpha)| {
            acc + alpha * (density[*pt as usize] - rho)
        })
        / grid.voronoi.volume
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::{Atoms, Lattice};
    use crate::grid::Grid;
    use crate::voxel_map::BlockingVoxelMap;

    // --- Helper to Setup Environment ---
    fn setup_env(dim: usize) -> (Grid, Atoms, BlockingVoxelMap) {
        let lattice_mat = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]];
        let origin = [0.0, 0.0, 0.0];

        let grid = Grid::new([dim, dim, dim], lattice_mat, origin);

        // Two atoms: one at [0,0,0] (corner), one at [1.5, 1.5, 1.5] (center)
        let atoms = Atoms::new(
            Lattice::new(lattice_mat),
            vec![[0.0, 0.0, 0.0], [1.5, 1.5, 1.5]],
            String::from("Test"),
        );

        let map = BlockingVoxelMap::new([dim, dim, dim], lattice_mat, origin);

        (grid, atoms, map)
    }

    // --- weight_step Tests ---

    #[test]
    fn test_weight_step_maximum() {
        let dim = 3;
        let (_, _, map) = setup_env(dim);

        // Setup density where center (index 13) is a peak
        let mut density = vec![0.0; 27];
        density[13] = 10.0; // Peak
        // Neighbors are 0.0 by default

        // Calling weight_step on the peak should return Maximum
        let result = weight_step(13, &density, &map, 1e-8);

        match result {
            WeightResult::Maximum => (),
            _ => panic!("Expected Maximum, got {:?}", result_name(&result)),
        }
    }

    #[test]
    fn test_weight_step_interior() {
        let dim = 3;
        let (_, _, map) = setup_env(dim);

        // Setup density gradient: 13 (Center) > 14 (Right Neighbor)
        let mut density = vec![0.0; 27];
        density[13] = 10.0;
        density[14] = 5.0;

        // Pre-populate the Maxima at 13 so 14 has somewhere to flow
        map.maxima_store(13, 100); // Atom 100

        // Step from 14. It should climb to 13, see Atom 100, and return Interior(100)
        let result = weight_step(14, &density, &map, 1e-8);

        match result {
            WeightResult::Interior(atom_idx) => assert_eq!(atom_idx, 100),
            _ => panic!("Expected Interior, got {:?}", result_name(&result)),
        }
    }

    #[test]
    fn test_weight_step_boundary() {
        let dim = 3;
        let (_, _, map) = setup_env(dim);

        // Voxel 13 is center.
        // Voxel 12 (Left) -> Atom 1
        // Voxel 14 (Right) -> Atom 2
        // We set 13 to be lower than both, so it flows "up" to both 12 and 14?
        // Actually weight_step goes UP gradient.
        // So let 13 be a Saddle between 12 and 14.

        let mut density = vec![0.0; 27];
        density[13] = 5.0; // Saddle
        density[12] = 10.0; // Peak A
        density[14] = 10.0; // Peak B

        map.maxima_store(12, EncodedAtom::new_zero_image(1).0 as isize);
        map.maxima_store(14, EncodedAtom::new_zero_image(2).0 as isize);

        // Step from 13. It should see both 12 and 14 as higher neighbors.
        let result = weight_step(13, &density, &map, 1e-8);

        match result {
            WeightResult::Boundary(weights)
            | WeightResult::Critical(weights) => {
                // Should have 2 weights
                assert_eq!(weights.len(), 2);
                // We expect ~50/50 split if geometry is symmetric
                let w1 = weights
                    .iter()
                    .find(|w| w.decode().0.atom_index() == 1)
                    .unwrap();
                let w2 = weights
                    .iter()
                    .find(|w| w.decode().0.atom_index() == 2)
                    .unwrap();

                let val1 = w1.decode().1;
                let val2 = w2.decode().1;

                assert!((val1 - 0.5).abs() < 0.1);
                assert!((val2 - 0.5).abs() < 0.1);
            }
            _ => panic!(
                "Expected Boundary/Critical, got {:?}",
                result_name(&result)
            ),
        }
    }

    // --- assign_maximum Tests ---

    #[test]
    fn test_assign_maximum_nearest() {
        let dim = 10;
        let (grid, atoms, _) = setup_env(dim); // 10x10x10 grid, Atoms at [0,0,0] and [1.5,1.5,1.5]

        // Voxel at [0,0,0] -> Index 0. Should match Atom 0.
        // Grid spacing is 3.0 / 10 = 0.3.
        let max_dist = 1.0;

        // Test Origin
        let atom_idx = assign_maximum(0, &atoms, &grid, &max_dist).unwrap();
        assert_eq!(atom_idx, 0);

        // Test Point closer to Atom 1 ([1.5, 1.5, 1.5] is at index ~555)
        // [5, 5, 5] -> 1.5, 1.5, 1.5
        let center_idx = 5 * 100 + 5 * 10 + 5;
        let atom_idx_2 =
            assign_maximum(center_idx, &atoms, &grid, &max_dist).unwrap();
        assert_eq!(atom_idx_2, 1);
    }

    #[test]
    fn test_assign_maximum_cutoff() {
        let dim = 10;
        let (grid, atoms, _) = setup_env(dim);

        // Atom 0 at [0,0,0]. Atom 1 at [1.5, 1.5, 1.5].
        // Try point at [2.9, 0, 0]. Closest to Atom 0 (distance 0.1 via PBC 3.0->0.0),
        // BUT let's test a point far from everything if possible,
        // or restrict max_dist very strictly.

        let tight_cutoff = 0.05; // Smaller than grid spacing (0.3)
        // Point [1, 0, 0] is at x=0.3. Dist to Atom 0 is 0.3.
        // 0.3 > 0.05 -> Should fail.
        let point_idx = 100; // x=1, y=0, z=0

        let result = assign_maximum(point_idx, &atoms, &grid, &tight_cutoff);
        assert!(result.is_err());
    }

    // --- laplacian Tests ---

    #[test]
    fn test_laplacian_uniform() {
        let dim = 3;
        let (grid, _, _) = setup_env(dim);

        // Uniform density -> Gradient is 0 -> Laplacian is 0
        let density = vec![1.0; 27];

        let lap = laplacian(13, &density, &grid);
        assert!(lap.abs() < 1e-9);
    }

    #[test]
    fn test_laplacian_peak() {
        let dim = 3;
        let (grid, _, _) = setup_env(dim);

        // Peak at center. Curvature should be negative (concave down).
        let mut density = vec![0.0; 27];
        density[13] = 10.0; // Center
        // Neighbors 0.0

        let lap = laplacian(13, &density, &grid);

        // The Laplacian formula involves sum of (neighbor - center).
        // (0 - 10) is negative. Sum is negative.
        assert!(lap < 0.0);
    }

    // Utility for debug printing enum variants
    fn result_name(w: &WeightResult) -> &'static str {
        match w {
            WeightResult::Maximum => "Maximum",
            WeightResult::Interior(_) => "Interior",
            WeightResult::Boundary(_) => "Boundary",
            WeightResult::Critical(_) => "Critical",
        }
    }
}
