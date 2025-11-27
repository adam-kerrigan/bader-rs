use crate::atoms::Atoms;
use crate::methods::laplacian;
use crate::progress::{Bar, HiddenBar, ProgressBar};
use crate::voxel_map::{Voxel, VoxelMap};
use std::thread;

type BaderResult = (Box<[Box<[f64]>]>, Box<[f64]>, Box<[f64]>, Box<[f64]>);

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
/// 3. **Radius**: `Box<[f64]>` - Distance from the nucleus to the furthest assigned voxel.
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
/// let (charges, volumes, radii, errors) = calculate_bader_density(
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
    let mut bader_radius = vec![f64::INFINITY; atoms.positions.len()];
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
                    let mut br = vec![f64::INFINITY; atoms.positions.len()];
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
                                        // Bader radius
                                        let atom_number = vm.maxima_to_atom(m);
                                        let p_c =
                                            vm.grid.to_cartesian(p as isize);
                                        let p_lll_c = atoms
                                            .lattice
                                            .cartesian_to_reduced(p_c);
                                        let atom = atoms.reduced_positions
                                            [atom_number];
                                        br[atom_number] =
                                            atoms.lattice.minimum_distance(
                                                p_lll_c,
                                                atom,
                                                Some(br[atom_number]),
                                            );
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
                    (bd, bv, br, be)
                })
            })
            .collect::<Vec<_>>();
        // Join each thread and collect the results.
        // If one thread terminates before the other this is not operated on first.
        // Either use the sorted index to remove vacuum from the summation or
        // find a way to operate on finshed threads first (ideally both).
        for thread in spawned_threads {
            if let Ok((tmp_bd, tmp_bv, tmp_br, tmp_be)) = thread.join() {
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
                bader_radius.iter_mut().zip(tmp_br.into_iter()).for_each(
                    |(a, b)| {
                        *a = a.min(b);
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
    bader_radius.iter_mut().for_each(|d| {
        match (*d).partial_cmp(&f64::INFINITY) {
            Some(std::cmp::Ordering::Less) => *d = d.powf(0.5),
            _ => *d = 0.0,
        }
    });
    (
        bader_density
            .into_iter()
            .map(|bd| bd.into_boxed_slice())
            .collect(),
        bader_volume.into(),
        bader_radius.into(),
        bader_error.into(),
    )
}

#[cfg(test)]
mod tests {
    use super::*; // Import functions from analysis.rs
    use crate::atoms::{Atoms, Lattice};
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
            vec![[0.0, 0.0, 0.0], [1.5, 1.5, 1.5]], // Positions
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

        let (b_rho, b_vol, _b_rad, _b_err) = calculate_bader_density(
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

        let (__rho, b_vol, b_rad, b_err) =
            calculate_bader_density(&density_data, &v_map, &atoms, 1, false);

        // Voxel Volume
        let v_vol = v_map.grid.voxel_lattice.volume;

        // Atom 0: 7 full voxels + 0.5 voxel
        let expected_vol_0 = (7.0 + 0.5) * v_vol;
        // Atom 1: 0 full voxels + 0.5 voxel
        let expected_vol_1 = 0.5 * v_vol;

        assert!((b_vol[0] - expected_vol_0).abs() < f64::EPSILON);
        assert!((b_vol[1] - expected_vol_1).abs() < f64::EPSILON);
        assert!(
            (b_rad[0] - (3.0 * 1.5f64.powi(2)).powf(0.5)).abs() < f64::EPSILON
        );
        assert!((b_err[0] - 0.0).abs() < f64::EPSILON);
    }
}
