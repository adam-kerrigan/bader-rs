use crate::atoms::Atoms;
use crate::critical::{CriticalPoint, CriticalPointKind};
use crate::errors::MaximaError;
use crate::grid::Grid;
use crate::hash::IntMap;
use crate::progress::progress_bar;
use crate::threading::parallel_map_reduce;
use crate::voxel_map::{
    BlockingVoxelMap, EncodedAtom, EncodedWeight, Voxel, VoxelMap,
};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

/// Result of a Weight step.
pub enum WeightResult {
    /// A boundary voxel which flows into a maximum that none of its higher
    /// neighbours do, making it a saddle (2 atoms), ring or cage candidate.
    Critical(Box<[EncodedWeight]>),
    /// Entirely assigned to a single Bader atom.
    Interior(EncodedAtom),
    /// Meeting point at the edge of 2 or more Bader atoms.
    Boundary(Box<[EncodedWeight]>),
    /// No higher neighbours, a local maximum.
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
/// * **Critical**: A boundary that flows into more atoms than any of its higher neighbours.
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
    let mut flux_total = 0.;
    let mut weights = IntMap::<EncodedAtom, f64>::default();
    let mut neighbour_atoms = 0;
    for ((pt, image), alpha) in voxel_map.grid.voronoi_shifts(p) {
        let charge_diff = density[pt as usize] - control;
        // only flow uphill, differences of zero are ignored to avoid division by zero
        if charge_diff <= 0. {
            continue;
        }
        let flux = charge_diff * alpha;
        flux_total += flux;
        match voxel_map.voxel_get(pt) {
            Voxel::Boundary(neighbour_weights) => {
                neighbour_atoms = neighbour_atoms.max(neighbour_weights.len());
                for (atom, weight) in
                    neighbour_weights.iter().map(|w| w.decode())
                {
                    *weights.entry(atom.image_add(image)).or_insert(0.) +=
                        weight as f64 * flux;
                }
            }
            Voxel::Maxima(atom) => {
                *weights.entry(atom.image_add(image)).or_insert(0.) += flux;
            }
            Voxel::Vacuum => panic!(
                "Vacuum voxel found with higher charge density than the control voxel."
            ),
        }
    }
    if weights.is_empty() {
        return WeightResult::Maximum;
    }
    // normalise the weights and remove any below the tolerance
    let weights = weights
        .into_iter()
        .map(|(atom, flux)| (atom, (flux / flux_total) as f32))
        .filter(|(_, weight)| *weight as f64 > weight_tolerance)
        .collect::<Vec<(EncodedAtom, f32)>>();
    match weights.as_slice() {
        [] => WeightResult::Maximum,
        [(atom, _)] => WeightResult::Interior(*atom),
        _ => {
            let encoded = weights
                .into_iter()
                .map(|(atom, w)| EncodedWeight::new(atom, w))
                .collect::<Box<[EncodedWeight]>>();
            // a new maximum has joined the weights -> Critical Point (saddle/ring/cage)
            if encoded.len() > neighbour_atoms {
                WeightResult::Critical(encoded)
            } else {
                WeightResult::Boundary(encoded)
            }
        }
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
/// * `visible_bar`: If `true`, displays a progress bar on `stderr`.
/// * `threads`: The number of worker threads to spawn.
///
/// # Returns
/// A tuple containing two vectors of [`CriticalPoint`]:
/// 1. **Bond Points**: Critical points connecting exactly 2 atoms.
/// 2. **Ring Points**: Critical points connecting 3 or more atoms.
///
/// # Example
/// ```
/// use bader::methods::weight;
/// use bader::voxel_map::{BlockingVoxelMap, VoxelMap, EncodedAtom, Voxel};
///
/// // 1. Setup Data: a 2x2x3 grid with a single maximum at voxel 2
/// let density = vec![2.0, 2.0, 12.0, 2.0, 2.0, 11.0, 1.0, 2.0, 6.0, 1.0, 1.0, 5.0];
/// let sorted_indices = vec![2, 5, 8, 11, 0, 1, 3, 4, 7, 6, 9, 10];
///
/// // 2. Setup Map
/// let map = BlockingVoxelMap::new(
///     [2, 2, 3],
///     [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
///     [0.0, 0.0, 0.0]
/// );
/// // Pre-assign the maxima (index 2) to Atom 0
/// map.maxima_store(2, EncodedAtom::new_zero_image(0));
///
/// // 3. Run Partitioning on the remaining voxels
/// let (bonds, rings) = weight(
///     &density,
///     &map,
///     &sorted_indices[1..],
///     1e-6,
///     false, // No progress bar
///     1      // Single thread
/// );
///
/// // 4. Check Results: everything flows to Atom 0, through the periodic
/// // boundary the atom is bonded to its own images so bonds are found
/// let map = VoxelMap::from_blocking_voxel_map(map);
/// assert!(matches!(map.voxel_get(1), Voxel::Maxima(atom) if atom.atom_index() == 0));
/// assert!(bonds.iter().chain(&rings).all(|cp| cp.atoms[0].atom_index() == 0));
/// ```
pub fn weight(
    density: &[f64],
    voxel_map: &BlockingVoxelMap,
    index: &[usize],
    weight_tolerance: f64,
    visible_bar: bool,
    threads: usize,
) -> (Vec<CriticalPoint>, Vec<CriticalPoint>) {
    let counter = AtomicUsize::new(0);
    let pbar = progress_bar(visible_bar, index.len(), "Bader Partitioning");
    let mut bonds = vec![];
    let mut rings = vec![];
    thread::scope(|s| {
        let handles = (0..threads)
            .map(|_| {
                s.spawn(|| {
                    let mut bonds = vec![];
                    let mut rings = vec![];
                    // voxels are claimed one at a time so that every thread is
                    // always working on the highest unassigned density
                    while let Some(p) = index
                        .get(counter.fetch_add(1, Ordering::Relaxed))
                        .map(|p| *p as isize)
                    {
                        match weight_step(
                            p,
                            density,
                            voxel_map,
                            weight_tolerance,
                        ) {
                            WeightResult::Maximum => {}
                            WeightResult::Interior(atom) => {
                                voxel_map.maxima_store(p, atom);
                            }
                            WeightResult::Boundary(weights) => {
                                voxel_map.weight_store(p, weights);
                            }
                            WeightResult::Critical(weights) => {
                                let atoms: Box<[EncodedAtom]> = weights
                                    .iter()
                                    .map(|w| w.decode().0)
                                    .collect();
                                voxel_map.weight_store(p, weights);
                                if atoms.len() == 2 {
                                    bonds.push(CriticalPoint::new(
                                        p,
                                        CriticalPointKind::Bond,
                                        atoms,
                                    ));
                                } else {
                                    rings.push(CriticalPoint::new(
                                        p,
                                        CriticalPointKind::Ring,
                                        atoms,
                                    ));
                                }
                            }
                        }
                        pbar.tick();
                    }
                    (bonds, rings)
                })
            })
            .collect::<Vec<_>>();
        for handle in handles {
            let (b, r) = handle.join().unwrap();
            bonds.extend(b);
            rings.extend(r);
        }
    });
    (bonds, rings)
}

/// Whether every Voronoi neighbour of `p` satisfies `pred(neighbour, rho)`
/// where `rho` is the density at `p`.
fn neighbours_all(
    p: usize,
    density: &[f64],
    grid: &Grid,
    pred: impl Fn(f64, f64) -> bool,
) -> bool {
    let rho = density[p];
    grid.voronoi_shifts_nocheck(p as isize)
        .all(|(pt, _)| pred(density[pt as usize], rho))
}

/// Scans the grid to identify all local maxima in the charge density.
///
/// This is the first step of the Bader analysis. It iterates over all voxels in `index`
/// (in parallel) and checks if a voxel is greater than or equal to all its neighbors.
///
/// Valid maxima are then assigned to atoms using [`assign_maximum`].
///
/// # Errors
/// Returns the first [`MaximaError`] from [`assign_maximum`] if a maximum is further
/// than `maximum_distance` from every atom.
pub fn maxima_finder(
    index: &[usize],
    density: &[f64],
    voxel_map: &BlockingVoxelMap,
    maximum_distance: f64,
    atoms: &Atoms,
    threads: usize,
    visible_bar: bool,
) -> Result<Vec<CriticalPoint>, MaximaError> {
    let grid = &voxel_map.grid;
    let pbar = progress_bar(visible_bar, index.len(), "Maxima Finding");
    parallel_map_reduce(
        index,
        || Ok(Vec::new()),
        |maxima: &mut Result<Vec<CriticalPoint>, MaximaError>, _, &p| {
            if let Ok(found) = maxima {
                if neighbours_all(p, density, grid, |neighbour, rho| {
                    neighbour <= rho
                }) {
                    match assign_maximum(
                        p as isize,
                        atoms,
                        grid,
                        maximum_distance,
                    ) {
                        Ok(atom) => found.push(CriticalPoint::new(
                            p as isize,
                            CriticalPointKind::Nuclei,
                            Box::new([EncodedAtom::new_zero_image(
                                atom as u32,
                            )]),
                        )),
                        Err(e) => *maxima = Err(e),
                    }
                }
            }
        },
        |global, local| match (global.as_mut(), local) {
            (Ok(g), Ok(l)) => g.extend(l),
            (Ok(_), Err(e)) => *global = Err(e),
            (Err(_), _) => {}
        },
        threads,
        pbar,
    )
}

/// Find the minima in the charge density.
///
/// A minimum that lies on a boundary between Bader volumes is a cage critical point
/// candidate and is returned with the atoms that share the voxel.
pub fn minima_finder(
    index: &[usize],
    density: &[f64],
    voxel_map: &VoxelMap,
    threads: usize,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let grid = &voxel_map.grid;
    let pbar = progress_bar(visible_bar, index.len(), "Minima Finding");
    parallel_map_reduce(
        index,
        Vec::new,
        |minima: &mut Vec<CriticalPoint>, _, &p| {
            // TODO: This needs to check if the cage is actually a boundary and if
            // not complain that the weight tolerance is too high
            if neighbours_all(p, density, grid, |neighbour, rho| {
                neighbour >= rho
            }) {
                if let Voxel::Boundary(weights) =
                    voxel_map.voxel_get(p as isize)
                {
                    minima.push(CriticalPoint::new(
                        p as isize,
                        CriticalPointKind::Cage,
                        weights.iter().map(|w| w.decode().0).collect(),
                    ));
                }
            }
        },
        |global, local| global.extend(local),
        threads,
        pbar,
    )
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
    maximum_distance: f64,
) -> Result<usize, MaximaError> {
    // convert the point first to cartesian, then to the reduced basis
    let m_cartesian = grid.to_cartesian(maximum);
    let m_reduced = atoms.lattice.cartesian_to_reduced(m_cartesian);
    let mut atom_num = 0;
    let mut min_distance = f64::INFINITY;
    for (i, atom) in atoms.reduced_positions.iter().enumerate() {
        let distance = atoms.lattice.minimum_distance(
            m_reduced,
            *atom,
            Some(min_distance),
        );
        if distance < min_distance {
            min_distance = distance;
            atom_num = i;
        }
    }
    let distance = min_distance.sqrt();
    if distance > maximum_distance {
        Err(MaximaError {
            maximum: m_cartesian,
            atom: atom_num,
            distance,
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
        .map(|(pt, alpha)| alpha * (density[pt as usize] - rho))
        .sum::<f64>()
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
        map.maxima_store(13, EncodedAtom::new_zero_image(100)); // Atom 100

        // Step from 14. It should climb to 13, see Atom 100, and return Interior(100)
        let result = weight_step(14, &density, &map, 1e-8);

        match result {
            WeightResult::Interior(atom) => assert_eq!(atom.atom_index(), 100),
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

        map.maxima_store(12, EncodedAtom::new_zero_image(1));
        map.maxima_store(14, EncodedAtom::new_zero_image(2));

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
        let atom_idx = assign_maximum(0, &atoms, &grid, max_dist).unwrap();
        assert_eq!(atom_idx, 0);

        // Test Point closer to Atom 1 ([1.5, 1.5, 1.5] is at index ~555)
        // [5, 5, 5] -> 1.5, 1.5, 1.5
        let center_idx = 5 * 100 + 5 * 10 + 5;
        let atom_idx_2 =
            assign_maximum(center_idx, &atoms, &grid, max_dist).unwrap();
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

        let result = assign_maximum(point_idx, &atoms, &grid, tight_cutoff);
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

    // --- maxima_finder / minima_finder Tests ---

    /// A 6x6x6 cubic cell with two Gaussian peaks on the two atoms.
    fn two_peak_density(grid: &Grid, atoms: &Atoms) -> Vec<f64> {
        (0..grid.size.total)
            .map(|p| {
                let r = grid.to_cartesian(p as isize);
                atoms
                    .positions
                    .iter()
                    .map(|a| {
                        let d = atoms.lattice.minimum_distance(
                            atoms.lattice.cartesian_to_reduced(r),
                            atoms.lattice.cartesian_to_reduced(*a),
                            None,
                        );
                        (-d / 0.5).exp()
                    })
                    .sum()
            })
            .collect()
    }

    #[test]
    fn test_maxima_finder_two_peaks() {
        let (grid, atoms, map) = setup_env(6);
        let density = two_peak_density(&grid, &atoms);
        let index = crate::utils::index_generator(&density, 0.0).unwrap();
        let maxima =
            maxima_finder(&index, &density, &map, 1.0, &atoms, 3, false)
                .unwrap();
        let mut found = maxima
            .iter()
            .map(|cp| (cp.atoms[0].atom_index(), cp.position))
            .collect::<Vec<_>>();
        found.sort_unstable();
        // atom 0 at the origin (voxel 0) and atom 1 at [1.5; 3] (voxel 3, 3, 3)
        assert_eq!(found, vec![(0, 0), (1, 3 * 36 + 3 * 6 + 3)]);
    }

    #[test]
    fn test_maxima_finder_too_far() {
        let (grid, atoms, map) = setup_env(6);
        let density = two_peak_density(&grid, &atoms);
        let index = crate::utils::index_generator(&density, 0.0).unwrap();
        // a cut-off smaller than the grid spacing still finds the on-atom
        // peaks, but moving the atoms makes every maximum too far away
        let moved = Atoms::new(
            Lattice::new(atoms.lattice.to_cartesian),
            vec![[0.25, 0.25, 0.25], [1.75, 1.75, 1.75]],
            String::new(),
        );
        assert!(
            maxima_finder(&index, &density, &map, 0.1, &moved, 2, false)
                .is_err()
        );
        assert!(
            maxima_finder(&index, &density, &map, 1.0, &moved, 2, false)
                .is_ok()
        );
    }

    #[test]
    fn test_weight_and_minima_two_peaks() {
        let (grid, atoms, map) = setup_env(6);
        let density = two_peak_density(&grid, &atoms);
        let index = crate::utils::index_generator(&density, 0.0).unwrap();
        let maxima =
            maxima_finder(&index, &density, &map, 1.0, &atoms, 2, false)
                .unwrap();
        for cp in &maxima {
            map.maxima_store(cp.position, cp.atoms[0]);
        }
        let (bonds, rings) = weight(&density, &map, &index, 1e-6, false, 3);
        // the two atoms share one bond (through the periodic images there is
        // only one unique pair) and there are no rings with two atoms
        assert!(!bonds.is_empty());
        assert!(bonds.iter().all(|cp| cp.atoms.len() == 2));
        assert!(rings.iter().all(|cp| cp.atoms.len() > 2));
        let map = VoxelMap::from_blocking_voxel_map(map);
        // every voxel is assigned
        assert!(map.voxel_map.iter().all(|&v| v != -1));
        let cages = minima_finder(&index, &density, &map, 3, false);
        assert!(!cages.is_empty());
        assert!(cages.iter().all(|cp| cp.kind == CriticalPointKind::Cage));
        // each cage is on a boundary so has more than one atom
        assert!(cages.iter().all(|cp| cp.atoms.len() > 1));
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
