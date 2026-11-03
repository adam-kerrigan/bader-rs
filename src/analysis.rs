use crate::atoms::Atoms;
use crate::methods::laplacian;
use crate::progress::progress_bar;
use crate::threading::parallel_map_reduce;
use crate::voxel_map::{Voxel, VoxelMap};

/// The integrated properties of each Bader volume.
///
/// Every field has one entry per atom followed by a final entry for the
/// vacuum, except `radius` which has no vacuum entry.
#[derive(Debug, Clone, PartialEq)]
pub struct Partition {
    /// The integrated charge of each density (e.g. charge, spin) for each atom.
    pub charge: Vec<Vec<f64>>,
    /// The volume of each Bader volume.
    pub volume: Vec<f64>,
    /// The distance from each atom to the closest boundary voxel of its
    /// volume, zero if the atom has no boundary.
    pub radius: Vec<f64>,
    /// The integrated Laplacian of each volume, a measure of the error in the
    /// partitioning that should be close to zero.
    pub error: Vec<f64>,
}

impl Partition {
    fn new(densities: usize, atoms: usize) -> Self {
        Self {
            charge: vec![vec![0.0; densities]; atoms + 1],
            volume: vec![0.0; atoms + 1],
            radius: vec![f64::INFINITY; atoms],
            error: vec![0.0; atoms + 1],
        }
    }

    /// The number of atoms, excluding the vacuum.
    pub fn atoms(&self) -> usize {
        self.radius.len()
    }

    /// Add the contribution `weight` of voxel `p` to the volume with index `atom`.
    fn add(
        &mut self,
        atom: usize,
        weight: f64,
        p: usize,
        density: &[Vec<f64>],
        laplacian: f64,
    ) {
        for (charge, rho) in self.charge[atom].iter_mut().zip(density) {
            *charge += weight * rho[p];
        }
        self.volume[atom] += weight;
        self.error[atom] += weight * laplacian;
    }

    fn merge(&mut self, other: Self) {
        for (a, b) in self.charge.iter_mut().zip(other.charge) {
            a.iter_mut().zip(b).for_each(|(a, b)| *a += b);
        }
        self.volume
            .iter_mut()
            .zip(other.volume)
            .for_each(|(a, b)| *a += b);
        self.radius
            .iter_mut()
            .zip(other.radius)
            .for_each(|(a, b)| *a = a.min(b));
        self.error
            .iter_mut()
            .zip(other.error)
            .for_each(|(a, b)| *a += b);
    }
}

/// Sums the densities of each Bader volume and calculates associated properties.
///
/// This function integrates the charge density over the volume assigned to each atom.
///
/// # Arguments
/// * `density` - A slice of density vectors (e.g., charge, spin).
/// * `voxel_map` - The map assigning each voxel to an atom or boundary.
/// * `atoms` - The atomic structure information.
/// * `threads` - Number of threads to use for parallel execution.
/// * `visible_bar` - Whether to show a progress bar.
///
/// # Returns
/// The [`Partition`] with the charge, volume, radius and error of each atom and the vacuum.
///
/// # Example
/// ```
/// use bader::analysis::calculate_bader_density;
/// use bader::atoms::{Atoms, Lattice};
/// use bader::grid::Grid;
/// use bader::voxel_map::VoxelMap;
///
/// let lattice = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
/// let atoms = Atoms::new(Lattice::new(lattice), vec![[0.0; 3]], String::new());
/// // a 2x2x2 grid with every voxel assigned to atom 0
/// let voxel_map = VoxelMap::new(
///     vec![0; 8],
///     vec![],
///     Grid::new([2, 2, 2], lattice, [0.0; 3]),
/// );
/// let charge_density = vec![vec![0.5; 8]];
///
/// let partition = calculate_bader_density(
///     &charge_density,
///     &voxel_map,
///     &atoms,
///     2,     // Run on 2 threads
///     false, // No progress bar
/// );
///
/// // the whole cell belongs to atom 0: volume 8 and charge 0.5 * 8
/// assert!((partition.volume[0] - 8.0).abs() < 1e-12);
/// assert!((partition.charge[0][0] - 4.0).abs() < 1e-12);
/// assert_eq!(partition.volume[1], 0.0); // vacuum
/// ```
pub fn calculate_bader_density(
    density: &[Vec<f64>],
    voxel_map: &VoxelMap,
    atoms: &Atoms,
    threads: usize,
    visible_bar: bool,
) -> Partition {
    let pbar =
        progress_bar(visible_bar, density[0].len(), "Summing Bader Density");
    let atom_len = atoms.positions.len();
    let grid = &voxel_map.grid;
    let mut partition = parallel_map_reduce(
        &voxel_map.voxel_map,
        || Partition::new(density.len(), atom_len),
        |partition, p, _| {
            let lapl = laplacian(p, &density[0], grid);
            match voxel_map.voxel_get(p as isize) {
                Voxel::Maxima(atom) => partition.add(
                    atom.atom_index() as usize,
                    1.0,
                    p,
                    density,
                    lapl,
                ),
                Voxel::Boundary(weights) => {
                    let p_reduced = atoms
                        .lattice
                        .cartesian_to_reduced(grid.to_cartesian(p as isize));
                    for (atom, weight) in weights.iter().map(|w| w.decode()) {
                        let atom = atom.atom_index() as usize;
                        partition.add(atom, weight as f64, p, density, lapl);
                        // the radius is the closest boundary voxel to the atom
                        partition.radius[atom] =
                            atoms.lattice.minimum_distance(
                                p_reduced,
                                atoms.reduced_positions[atom],
                                Some(partition.radius[atom]),
                            );
                    }
                }
                Voxel::Vacuum => partition.add(atom_len, 1.0, p, density, lapl),
            }
        },
        Partition::merge,
        threads,
        pbar,
    );
    // convert from density to charge and from voxel count to volume
    let voxel_volume = grid.voxel_lattice.volume;
    for charge in partition.charge.iter_mut() {
        charge.iter_mut().for_each(|c| *c *= voxel_volume);
    }
    partition.volume.iter_mut().for_each(|v| *v *= voxel_volume);
    // the distances are squared and atoms without a boundary are infinite
    for radius in partition.radius.iter_mut() {
        *radius = if radius.is_finite() {
            radius.sqrt()
        } else {
            0.0
        };
    }
    partition
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

        let partition = calculate_bader_density(
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
            (partition.charge[0][0] - total_vol).abs() < 1e-4,
            "Atom 0 should contain all density"
        );
        assert!(
            (partition.volume[0] - total_vol).abs() < 1e-4,
            "Atom 0 should have full volume"
        );

        // Atom 1 should be empty
        assert_eq!(partition.charge[1][0], 0.0);
        assert_eq!(partition.volume[1], 0.0);
        assert_eq!(partition.radius[1], 0.0);
        assert_eq!(partition.atoms(), 2);
    }

    #[test]
    fn test_calculate_bader_density_thread_independent() {
        // the result must not depend on the number of threads, including
        // thread counts that do not divide the number of voxels
        let dim = 5;
        let size = dim * dim * dim;
        let (grid, atoms) = setup_env(dim);
        let density_data = vec![
            (0..size)
                .map(|i| 1.0 + (i % 7) as f64)
                .collect::<Vec<f64>>(),
            (0..size)
                .map(|i| (i % 3) as f64 - 1.0)
                .collect::<Vec<f64>>(),
        ];
        let voxel_indices = (0..size as isize)
            .map(|i| if i % 11 == 0 { -1 } else { i % 2 })
            .collect();
        let v_map = VoxelMap::new(voxel_indices, vec![], grid);
        let reference =
            calculate_bader_density(&density_data, &v_map, &atoms, 1, false);
        for threads in [2, 3, 4, 7] {
            let partition = calculate_bader_density(
                &density_data,
                &v_map,
                &atoms,
                threads,
                false,
            );
            for (a, b) in partition.charge.iter().zip(&reference.charge) {
                for (a, b) in a.iter().zip(b) {
                    assert!((a - b).abs() < 1e-9, "{} threads", threads);
                }
            }
            for (a, b) in partition.error.iter().zip(&reference.error) {
                assert!((a - b).abs() < 1e-9, "{} threads", threads);
            }
            assert_eq!(partition.volume, reference.volume);
        }
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

        let partition =
            calculate_bader_density(&density_data, &v_map, &atoms, 1, false);

        // Voxel Volume
        let v_vol = v_map.grid.voxel_lattice.volume;

        // Atom 0: 7 full voxels + 0.5 voxel
        let expected_vol_0 = (7.0 + 0.5) * v_vol;
        // Atom 1: 0 full voxels + 0.5 voxel
        let expected_vol_1 = 0.5 * v_vol;

        assert!((partition.volume[0] - expected_vol_0).abs() < f64::EPSILON);
        assert!((partition.volume[1] - expected_vol_1).abs() < f64::EPSILON);
        assert!(
            (partition.radius[0] - (3.0 * 1.5f64.powi(2)).sqrt()).abs()
                < f64::EPSILON
        );
        assert!((partition.error[0] - 0.0).abs() < f64::EPSILON);
    }
}
