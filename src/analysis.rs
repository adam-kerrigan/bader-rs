use crate::atoms::Atoms;
use crate::critical::CriticalPoint;
use crate::grid::Grid;
use crate::methods::laplacian;
use crate::progress::{Bar, HiddenBar, ProgressBar};
use crate::voxel_map::{Voxel, VoxelMap};
use std::thread;

type BaderResult = (Box<[Box<[f64]>]>, Box<[f64]>, Box<[f64]>);

/// Sums the densities of each Bader volume and calculates associated properties.
///
/// This function integrates the charge density over the volume assigned to each atom.
/// It returns four vectors containing the results for each atom (plus a final index for vacuum).
///
/// # Arguments
/// * `density` - A slice of density vectors (e.g., charge, spin).
/// * `voxel_map` - The map assigning each voxel to an atom or boundary.
/// * `atoms` - The atomic structure information.
/// * `threads` - Number of threads to use for parallel execution.
/// * `visible_bar` - Whether to show a progress bar.
///
/// # Returns
/// A tuple containing:
/// 1. **Density**: `Box<[Box<[f64]>]>` - Integrated charge for each atom.
/// 2. **Volume**: `Box<[f64]>` - Volume size for each atom.
/// 4. **Error**: `Box<[f64]>` - Integrated Laplacian (should be close to 0 for perfect partitioning).
///
/// # Example
/// ```no_run
/// use bader::analysis::calculate_bader_density;
/// use bader::voxel_map::VoxelMap;
/// use bader::atoms::Atoms;
///
/// // Assume we have a populated VoxelMap and Atoms struct
/// # let voxel_map: VoxelMap = unsafe { std::mem::zeroed() };
/// # let atoms: Atoms = unsafe { std::mem::zeroed() };
/// // And a calculated charge density flat array
/// let charge_density = vec![vec![0.1; 1000]]; // 1000 voxels
///
/// let (charges, volumes, errors) = calculate_bader_density(
///     &charge_density,
///     &voxel_map,
///     &atoms,
///     4,    // Run on 4 threads
///     true, // Show progress bar
/// );
///
/// println!("Atom 0 Charge: {}", charges[0][0]);
/// println!("Atom 0 Volume: {}", volumes[0]);
/// ```
pub fn calculate_bader_density(
    density: &[Vec<f64>],
    voxel_map: &VoxelMap,
    atoms: &Atoms,
    threads: usize,
    visible_bar: bool,
) -> BaderResult {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            density[0].len(),
            String::from("Summing Bader Density"),
        )),
    };
    let pbar = &progress_bar;
    let mut bader_density =
        vec![vec![0.0; density.len()]; atoms.positions.len() + 1];
    let mut bader_volume = vec![0.0; atoms.positions.len() + 1];
    let mut bader_error = vec![0.0; atoms.positions.len() + 1];
    let vm = &voxel_map;
    // Calculate the size of the vector to be passed to each thread.
    let chunk_size =
        (density[0].len() / threads) + (density[0].len() % threads).min(1);
    thread::scope(|s| {
        let spawned_threads = voxel_map
            .maxima_chunks(chunk_size)
            .enumerate()
            .map(|(index, chunk)| {
                s.spawn(move || {
                    let mut bd = vec![
                        vec![0.0; density.len()];
                        atoms.positions.len() + 1
                    ];
                    let mut bv = vec![0.0; atoms.positions.len() + 1];
                    let mut be = vec![0.0; atoms.positions.len() + 1];
                    chunk.iter().enumerate().for_each(
                        |(voxel_index, maxima)| {
                            let p = index * chunk.len() + voxel_index;
                            let lapl = laplacian(p, &density[0], &vm.grid);
                            match vm.maxima_to_voxel(*maxima) {
                                Voxel::Maxima(m) => {
                                    let m = m.atom_index() as usize;
                                    // Bader density
                                    for (i, rho) in density.iter().enumerate() {
                                        bd[m][i] += rho[p];
                                    }
                                    // Bader Volume
                                    bv[m] += 1.0;
                                    // Bader Error
                                    be[m] += lapl
                                }
                                Voxel::Boundary(weights) => {
                                    for (m, w) in weights.into_iter() {
                                        let m = m.atom_index() as usize;
                                        let w = w as f64;
                                        // Bader density
                                        for (i, rho) in
                                            density.iter().enumerate()
                                        {
                                            bd[m][i] += w * rho[p];
                                        }
                                        // Bader volume
                                        bv[m] += w;
                                        // Bader error
                                        be[m] += w * lapl;
                                    }
                                }
                                Voxel::Vacuum => {
                                    // Bader Density
                                    for (i, rho) in density.iter().enumerate() {
                                        bd[atoms.positions.len()][i] += rho[p];
                                    }
                                    // Bader Volume
                                    bv[atoms.positions.len()] += 1.0;
                                    // Bader Error
                                    be[atoms.positions.len()] += lapl
                                }
                            };
                            pbar.tick();
                        },
                    );
                    (bd, bv, be)
                })
            })
            .collect::<Vec<_>>();
        // Join each thread and collect the results.
        // If one thread terminates before the other this is not operated on first.
        // Either use the sorted index to remove vacuum from the summation or
        // find a way to operate on finshed threads first (ideally both).
        for thread in spawned_threads {
            if let Ok((tmp_bd, tmp_bv, tmp_be)) = thread.join() {
                bader_density.iter_mut().zip(tmp_bd.into_iter()).for_each(
                    |(a, b)| {
                        a.iter_mut().zip(b).for_each(|(c, d)| *c += d);
                    },
                );
                bader_volume.iter_mut().zip(tmp_bv.into_iter()).for_each(
                    |(a, b)| {
                        *a += b;
                    },
                );
                bader_error.iter_mut().zip(tmp_be.into_iter()).for_each(
                    |(a, b)| {
                        *a += b;
                    },
                );
            } else {
                panic!("Unable to join thread in sum_bader_densities.")
            };
        }
    });
    // The final result needs to be converted to a charge rather than a density.
    bader_density.iter_mut().for_each(|a| {
        a.iter_mut()
            .for_each(|b| *b *= voxel_map.grid_get().voxel_lattice.volume);
    });
    // The distance isn't square rooted in the calcation of distance to save time.
    // As we need to filter out the infinite distances (atoms with no assigned maxima)
    // we can square root here also.
    bader_volume.iter_mut().for_each(|a| {
        *a *= voxel_map.grid_get().voxel_lattice.volume;
    });
    (
        bader_density
            .into_iter()
            .map(|bd| bd.into_boxed_slice())
            .collect(),
        bader_volume.into(),
        bader_error.into(),
    )
}

/// Calculates the Bader radius for each atom in the system.
///
/// The Bader radius is defined as the minimum distance from an atom's nucleus to any of
/// its associated Bond Critical Points (BCPs).
///
/// # Logic
/// 1. Initializes the radius for all atoms to infinity.
/// 2. Iterates through the given list of bond critical points.
/// 3. For each critical point, calculates the shortest periodic distance to the parent atoms.
/// 4. Updates the atom's minimum radius if the current BCP is closer than previously checked ones.
/// 5. Finally, takes the square root of the minimum squared distances to get the true radii.
///    If an atom has no associated bonds (distance is still infinity), its radius is set to `0.0`.
///
/// # Arguments
/// * `bonds` - A slice containing the [`CriticalPoint`]s to be measured.
///   The outer index corresponds to the atom's ID in the `Atoms` struct.
/// * `atoms` - The atomic geometry, containing reduced positions and lattice vectors.
/// * `grid` - The voxel grid, used to translate 1D grid positions into Cartesian coordinates.
/// * `visible_bar` - A boolean flag to determine whether to output a progress bar to the console.
///
/// # Returns
/// A boxed slice of `f64` values (`Box<[f64]>`), where the value at index `i` is the calculated
/// Bader radius for Atom `i`.
pub fn calculate_bader_radius(
    bonds: &[CriticalPoint],
    atoms: &Atoms,
    grid: &Grid,
    visible_bar: bool,
) -> Box<[f64]> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            bonds.len(),
            String::from("Calculating Bader Radii"),
        )),
    };
    let pbar = &progress_bar;
    let mut radii = vec![f64::INFINITY; atoms.positions.len()];
    bonds.iter().for_each(|critical_point| {
        let p_c = grid.to_cartesian(critical_point.position);
        let p_lll_c = atoms.lattice.cartesian_to_reduced(p_c);
        critical_point.atoms.iter().for_each(|encoded_atom| {
            let atom_number = encoded_atom.atom_index() as usize;
            let atom = atoms.reduced_positions[atom_number];
            println!("{:?} {:?}", atom, p_lll_c);
            radii[atom_number] = atoms.lattice.minimum_distance(
                p_lll_c,
                atom,
                Some(radii[atom_number]),
            );
        });
        pbar.tick();
    });
    radii
        .iter_mut()
        .for_each(|d| match (*d).partial_cmp(&f64::INFINITY) {
            Some(std::cmp::Ordering::Less) => *d = d.powf(0.5),
            _ => *d = 0.0,
        });
    radii.into()
}

#[cfg(test)]
mod tests {
    use super::*; // Import functions from analysis.rs
    use crate::atoms::{Atoms, Lattice};
    use crate::critical::{CriticalPoint, CriticalPointKind};
    use crate::grid::Grid;
    use crate::voxel_map::{
        EncodedAtom, EncodedImage, EncodedWeight, VoxelMap,
    };

    // Helper to create a basic environment
    fn setup_env(grid_dim: usize) -> (Grid, Atoms) {
        let lattice_mat = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]];
        let origin = [0.0, 0.0, 0.0];

        // 1. Setup Grid
        let grid =
            Grid::new([grid_dim, grid_dim, grid_dim], lattice_mat, origin);

        // 2. Setup Atoms (2 atoms at opposite corners)
        let atoms = Atoms::new(
            Lattice::new(lattice_mat),
            vec![[0.0, 0.0, 0.0], [1.5, 1.5, 1.5], [2.0, 0.0, 0.0]], // Positions
            String::from("Test System"),
        );

        (grid, atoms)
    }

    #[test]
    fn test_calculate_bader_density_uniform() {
        let dim = 4;
        let size = dim * dim * dim;
        let (grid, atoms) = setup_env(dim);

        // CASE: Uniform density of 1.0 everywhere
        // All voxels belong to Atom 0
        let density_data = vec![vec![1.0; size]];

        let voxel_indices = vec![0isize; size]; // All point to Atom 0
        let weight_map = vec![]; // No boundaries

        // Construct VoxelMap manually
        let v_map = VoxelMap {
            voxel_map: voxel_indices,
            weight_map,
            grid,
        };

        let (b_rho, b_vol, _b_err) = calculate_bader_density(
            &density_data,
            &v_map,
            &atoms,
            1,     // 1 thread
            false, // no progress bar
        );

        // Atom 0 should have all the charge
        // Total Volume = 3.0 * 3.0 * 3.0 = 27.0
        let total_vol = 27.0;

        assert!(
            (b_rho[0][0] - total_vol).abs() < 1e-4,
            "Atom 0 should contain all density"
        );
        assert!(
            (b_vol[0] - total_vol).abs() < 1e-4,
            "Atom 0 should have full volume"
        );

        // Atom 1 should be empty
        assert_eq!(b_rho[1][0], 0.0);
        assert_eq!(b_vol[1], 0.0);
    }

    #[test]
    fn test_calculate_bader_density_boundary() {
        let dim = 2; // 8 voxels total
        let size = 8;
        let (grid, atoms) = setup_env(dim);

        // CASE: 7 voxels for Atom 0, 1 voxel on boundary (50/50 split)
        let density_data = vec![vec![1.0; size]];

        let mut voxel_indices = vec![0isize; size];

        // Create Weights: Voxel 7 is shared 50/50 between Atom 0 and Atom 1
        let atom0 = EncodedAtom::new(0, EncodedImage::new([0, 0, 0]));
        let atom1 = EncodedAtom::new(1, EncodedImage::new([0, 0, 0]));

        let w0 = EncodedWeight::new(atom0, 0.5);
        let w1 = EncodedWeight::new(atom1, 0.5);

        let weights = vec![w0, w1].into_boxed_slice();

        // Set Voxel 7 to be boundary (index -2) pointing to weight_map[0]
        voxel_indices[7] = -2;
        let weight_map = vec![weights];

        let v_map = VoxelMap {
            voxel_map: voxel_indices,
            weight_map,
            grid,
        };

        let (__rho, b_vol, b_err) =
            calculate_bader_density(&density_data, &v_map, &atoms, 1, false);

        // Voxel Volume
        let v_vol = v_map.grid.voxel_lattice.volume;

        // Atom 0: 7 full voxels + 0.5 voxel
        let expected_vol_0 = (7.0 + 0.5) * v_vol;
        // Atom 1: 0 full voxels + 0.5 voxel
        let expected_vol_1 = 0.5 * v_vol;

        assert!((b_vol[0] - expected_vol_0).abs() < f64::EPSILON);
        assert!((b_vol[1] - expected_vol_1).abs() < f64::EPSILON);
        assert!((b_err[0] - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_calculate_bader_radius() {
        let dim = 10; // 1000 voxels total
        let (grid, atoms) = setup_env(dim);

        // Create dummy Critical Points using the helper from earlier.
        // We will assign two BCPs to Atom 0, and ZERO BCPs to Atom 1.
        let cp_closest = CriticalPoint::new(
            100,
            CriticalPointKind::Bond,
            Box::new([
                EncodedAtom::new_zero_image(0),
                EncodedAtom::new_zero_image(2),
            ]),
            0.0,
            0.0,
        ); // Distance = 0.3, 2.7 (1.3 from pbc)
        let cp_further = CriticalPoint::new(
            800,
            CriticalPointKind::Bond,
            Box::new([
                EncodedAtom::new_zero_image(0),
                EncodedAtom::new_zero_image(2),
            ]),
            0.0,
            0.0,
        ); // Distance = 0.6, 0.4

        // Group the critical points by atom (outer index is atom number)
        let bonds = vec![cp_closest, cp_further];

        // Run the function silently (visible_bar = false)
        let radii = calculate_bader_radius(&bonds, &atoms, &grid, false);

        // Assertions
        assert_eq!(radii.len(), 3);

        // Atom 0 should have a radius corresponding to the closest BCP (distance = 0.3)
        assert!((radii[0] - 0.3).abs() < f64::EPSILON);

        // Atom 1 has no bonds, so its radius should default to 0.0
        assert!((radii[1] - 0.0).abs() < f64::EPSILON);
        // Atom 2 should be 0.4
        assert!((radii[2] - 0.4).abs() < f64::EPSILON);
    }
}
