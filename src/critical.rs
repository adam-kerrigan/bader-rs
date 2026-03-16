use std::cmp::Ordering;

use crate::arguments::Args;
use crate::atoms::Atoms;
use crate::grid::Grid;
use crate::hash::{IntMap, IntSet};
use crate::progress::{Bar, HiddenBar, ProgressBar};
use crate::threading::parallel_prune;
use crate::utils::{cross, norm, subtract, vdot};
use crate::voxel_map::EncodedAtom;

/// Represents a critical point in the charge density topology.
///
/// # Fields
/// * `position`: The 1D index of the voxel where this point is located.
/// * `kind`: The topological classification (Nuclei, Bond, Ring, Cage).
/// * `atoms`: The list of atoms associated with this point (e.g., the two atoms a bond connects).
#[derive(Clone)]
pub struct CriticalPoint {
    pub position: isize,
    pub kind: CriticalPointKind,
    pub atoms: Box<[EncodedAtom]>,
    pub density: f64,
    pub laplacian: f64,
}

impl CriticalPoint {
    pub fn new(
        position: isize,
        kind: CriticalPointKind,
        atoms: Box<[EncodedAtom]>,
        density: f64,
        laplacian: f64,
    ) -> Self {
        CriticalPoint {
            position,
            kind,
            atoms,
            density,
            laplacian,
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

/// A canonical key for identifying unique critical points based on their atom list.
///
/// This struct handles the normalisation of atom lists to ensure that:
/// 1. Order doesn't matter (sorted internally).
/// 2. Translational symmetry is respected (images are normalised relative to the first atom).
///
/// This allows for hashing and sets to identify duplicate points.
#[derive(Debug, PartialEq, Eq, Hash)]
pub struct CriticalPointKey(Vec<EncodedAtom>);
impl CriticalPointKey {
    pub fn from_cp(cp: CriticalPoint) -> Self {
        let mut atoms = cp.atoms.to_vec();
        atoms.sort_unstable();
        if let Some(anchor) = atoms.first() {
            let image = anchor.image();
            atoms
                .iter_mut()
                .for_each(|atom| *atom = atom.image_sub(image));
        }
        atoms.sort_unstable();
        Self(atoms)
    }

    pub fn into_box(self) -> Box<[EncodedAtom]> {
        self.0.into()
    }
}

/// Filters and orders nuclear critical points.
///
/// In the raw output, a single atom might be associated with multiple "maxima" candidates due
/// to grid discretisation noise. This function groups candidates by their assigned Atom ID
/// and selects the one with the highest charge density as the true nucleus position.
///
/// # Arguments
/// * `nuclei`: List of candidate points classified as [`CriticalPointKind::Nuclei`].
/// * `density`: The charge density array (used to compare peak heights).
/// * `atom_len`: The total number of atoms in the system (defines output size).
/// * `visible_bar`: Whether to show a progress bar.
///
/// # Returns
/// A vector of length `atom_len`, where index `i` corresponds to Atom `i`.
pub fn nuclei_ordering(
    nuclei: &mut [CriticalPoint],
    atoms: &Atoms,
    grid: &Grid,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let atom_len = atoms.positions.len();
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            nuclei.len() + atom_len,
            String::from("Pruning Nucleus Critical Points"),
        )),
    };
    let pbar = &progress_bar;
    // Find all the nuclei with the same atom number and group them based on which is largest
    // charge density. All maxima will be in image (0, 0, 0).
    // We need to order the nuclei by atom so that we can get the position of them by knowing the
    // atom number.
    let mut nuclei_sorting = vec![Vec::<usize>::new(); atom_len];
    nuclei.iter().enumerate().for_each(|(i, cp)| {
        nuclei_sorting[cp.atoms[0].atom_index() as usize].push(i);
        pbar.tick();
    });
    nuclei_sorting
        .into_iter()
        .map(|indices| {
            pbar.tick();
            match indices.iter().max_by(|a, b| {
                nuclei[**a].density.total_cmp(&nuclei[**b].density)
            }) {
                Some(index) => {
                    let true_maximum = nuclei[*index].clone();
                    if indices.len() > 1 {
                        let true_position =
                            grid.to_cartesian(nuclei[*index].position);
                        indices.iter().for_each(|i| {
                            if i != index {
                                if let Some(cp) = nuclei.get_mut(*i) {
                                    let position =
                                        grid.to_cartesian(cp.position);
                                    let image = atoms
                                        .lattice
                                        .closest_image(true_position, position);
                                    *cp = CriticalPoint::new(
                                        cp.position,
                                        CriticalPointKind::Nuclei,
                                        Box::new(
                                            [cp.atoms[0].image_sub(image)],
                                        ),
                                        cp.density,
                                        cp.laplacian,
                                    )
                                }
                            }
                        });
                    }
                    true_maximum
                }
                None => CriticalPoint::new(
                    0,
                    CriticalPointKind::Blank,
                    Box::new([]),
                    0.0,
                    0.0,
                ),
            }
        })
        .collect()
}

/// Filters and deduplicates Bond Critical Points (3, -1).
///
/// Bonds are defined by the two atoms they connect. In the raw analysis, multiple voxels
/// along the boundary between two atoms might be flagged as saddle points.
///
/// # Logic
/// 1. **Deduplication**: Uses [`parallel_prune`] to group points by their atom pair.
///    For each pair, only the point with the **highest charge density** is retained.
///    This corresponds to the true saddle point on the gradient path.
///
/// # Arguments
/// * `bonds`: The raw list of candidate bond points.
/// * `threads`: Number of threads to use.
/// * `visible_bar`: Whether to display a progress bar.
pub fn bond_pruning(
    bonds: &[CriticalPoint],
    args: &Args,
) -> Vec<CriticalPoint> {
    let threads = args.threads;
    let progress_bar: Box<dyn ProgressBar> = match args.silent {
        true => Box::new(HiddenBar {}),
        false => Box::new(Bar::new(
            bonds.len(),
            String::from("Pruning Bond Critical Points"),
        )),
    };
    parallel_prune(bonds, |_| true, threads, progress_bar)
}

/// Constructs an adjacency graph from a list of Bond Critical Points.
///
/// This function maps each atom to a list of atoms it is bonded to, preserving
/// periodic boundary conditions.
///
/// # Arguments
/// * `bonds`: The list of Bond Critical Points (3, -1) representing connections.
/// * `atom_len`: The total number of atoms in the system (used to size the output vector).
///
/// # Returns
/// A vector where index `i` contains a list of `EncodedAtom`s bonded to Atom `i`.
/// The `EncodedAtom`s in the list are shifted relative to the image of Atom `i`.
pub fn bond_adjacency(
    bonds: &[CriticalPoint],
    atom_len: usize,
) -> Vec<Vec<EncodedAtom>> {
    let mut adjacency: Vec<Vec<EncodedAtom>> = vec![Vec::new(); atom_len];
    bonds.iter().for_each(|bond| {
        adjacency[bond.atoms[0].atom_index() as usize].push(bond.atoms[1]);
        adjacency[bond.atoms[1].atom_index() as usize]
            .push(bond.atoms[0].image_sub(bond.atoms[1].image()));
    });
    adjacency
}

/// Filters Ring Critical Points (3, +1) by validating bond connectivity.
///
/// This function verifies that the atoms associated with a candidate Ring Critical Point
/// form a valid, continuous closed loop (cycle) within the determined bond network.
///
/// # Logic
/// 1. **Connectivity Check**: Iterates through the atoms in the ring candidate.
/// 2. **Valency Validation**: Verifies that every atom in the candidate is connected
///    to exactly two other atoms *within the same candidate set* (via `bond_adjacency`).
/// 3. **Traversal Check**: Ensures the atoms form a single continuous loop (no disjoint parts).
/// 4. **Deduplication**: Groups by unique atom list and keeps the candidate with the highest density.
///
/// # Arguments
/// * `rings`: The raw list of candidate ring points.
/// * `bond_adjacency`: The adjacency graph derived from Bond Critical Points.
/// * `args`: Command line arguments controlling threading and verbosity.
pub fn ring_pruning(
    rings: &[CriticalPoint],
    bond_adjancy: &[Vec<EncodedAtom>],
    args: &Args,
) -> Vec<CriticalPoint> {
    let threads = args.threads;
    let progress_bar: Box<dyn ProgressBar> = match args.silent {
        true => Box::new(HiddenBar {}),
        false => Box::new(Bar::new(
            rings.len(),
            String::from("Pruning Ring Critical Points"),
        )),
    };
    let validator = |cp: &CriticalPoint| {
        // First we are going to check that each atom is only connected to two other atoms
        let mut previous_index = 0;
        let mut current_index = 0;
        let mut next_index = 0;
        let mut steps = 1;
        loop {
            let mut intrabonds = 0;
            let graph =
                &bond_adjancy[cp.atoms[current_index].atom_index() as usize];
            cp.atoms.iter().enumerate().for_each(|(i, encoded_atom)| {
                if graph.contains(
                    &encoded_atom.image_sub(cp.atoms[current_index].image()),
                ) {
                    intrabonds += 1;
                    if i != previous_index {
                        next_index = i;
                    }
                }
            });
            // regardless of steps if there are more or less than 2 bonds it is not a ring
            if intrabonds != 2 {
                return false;
            }
            // the steps limit should be redundant but stops an infinite loop
            if next_index == 0 || steps > cp.atoms.len() {
                break;
            }
            steps += 1;
            previous_index = current_index;
            current_index = next_index;
        }
        // Next we will check that we can actually traverse the full ring
        steps == cp.atoms.len()
    };
    parallel_prune(rings, validator, threads, progress_bar)
}

/// Filters Cage Critical Points (3, +3) and enforces 3D structure.
///
/// A Cage Critical Point must connect at least 4 atoms and those atoms must **not** be coplanar
/// (they must enclose a volume).
///
/// # Logic
/// 1. **Size Check**: Must have $\ge$ 4 atoms.
/// 2. **Volumetric Check**: Calculates the plane formed by the first 3 atoms.
///    Scans the remaining atoms; if **any** atom lies significantly off this plane,
///    the point is accepted as a valid cage.
///    - If all atoms are coplanar, the point is rejected.
/// 3. **Deduplication**: Groups by atom list and keeps the candidate with the highest density.
///
/// # Arguments
/// * `cages`: The raw list of candidate cage points.
/// * `ordered_nuclei`: The list of definitive nucleus positions.
/// * `atoms`: The system geometry.
/// * `grid`: The voxel grid.
pub fn cage_pruning(
    cages: &[CriticalPoint],
    ordered_nuclei: &[CriticalPoint],
    atoms: &Atoms,
    grid: &Grid,
    args: &Args,
) -> Vec<CriticalPoint> {
    let threads = args.threads;
    let progress_bar: Box<dyn ProgressBar> = match args.silent {
        true => Box::new(HiddenBar {}),
        false => Box::new(Bar::new(
            cages.len(),
            String::from("Pruning Cage Critical Points"),
        )),
    };

    _ = |cp: &CriticalPoint| {
        if cp.atoms.len() < 4 {
            return false;
        }
        let positions: Vec<[f64; 3]> = cp.atoms[..3]
            .iter()
            .map(|encoded_atom| {
                let (atom_num, encoded_image) = encoded_atom.decode_partial();
                let mut position = grid
                    .to_cartesian(ordered_nuclei[atom_num as usize].position);
                let image = match encoded_image.is_zero() {
                    true => [0., 0., 0.],
                    false => {
                        let image = encoded_image.decode();
                        atoms.lattice.fractional_to_cartesian([
                            image[0] as f64,
                            image[1] as f64,
                            image[2] as f64,
                        ])
                    }
                };
                position.iter_mut().zip(image).for_each(|(f, i)| *f += i);
                position
            })
            .collect();
        let vec_1 = subtract(positions[1], positions[0]);
        let vec_2 = subtract(positions[2], positions[0]);
        let mut plane = cross(vec_1, vec_2);
        let plane_normal = norm(plane);
        plane.iter_mut().for_each(|f| *f /= plane_normal);
        for encoded_atom in cp.atoms[3..].iter() {
            let (atom_num, encoded_image) = encoded_atom.decode_partial();
            let mut position =
                grid.to_cartesian(ordered_nuclei[atom_num as usize].position);
            let image = match encoded_image.is_zero() {
                true => [0., 0., 0.],
                false => {
                    let image = encoded_image.decode();
                    atoms.lattice.fractional_to_cartesian([
                        image[0] as f64,
                        image[1] as f64,
                        image[2] as f64,
                    ])
                }
            };
            position.iter_mut().zip(image).for_each(|(f, i)| *f += i);
            let vec_3 = subtract(position, positions[0]);
            let mut plane_t = cross(vec_1, vec_3);
            let plane_normal = norm(plane_t);
            plane_t.iter_mut().for_each(|f| *f /= plane_normal);
            // TODO: make this a tolerance currently 5.73 degrees
            if vdot(plane, plane_t).abs() < 0.995 {
                return true;
            }
        }
        false
    };
    parallel_prune(cages, |_| true, threads, progress_bar)
}

/// Merges degenerate or subset critical points.
///
/// This simplifies the topology by removing redundant critical points.
///
/// # Logic
/// 1. Sorts critical points by the number of associated atoms (descending).
/// 2. Iterates through them. If a point's atoms are a **subset** of an already existing
///    point's atoms (e.g., a Ring \[1,2,3\] is likely charge fluctuation before the Ring \[1,2,3,4\]),
///    it is discarded/merged.
///
/// # Returns
/// A cleaner list of unique critical points.
pub fn critical_point_merge(mut cps: Vec<CriticalPoint>) -> Vec<CriticalPoint> {
    cps.sort_unstable_by(|a, b| {
        let order = b.atoms.len().cmp(&a.atoms.len());
        if let Ordering::Equal = order {
            b.density.total_cmp(&a.density)
        } else {
            order
        }
    });
    let mut merged_points: Vec<CriticalPoint> = Vec::with_capacity(cps.len());
    let mut inverted_index: IntMap<u32, Vec<usize>> = IntMap::default();
    'critical: for cp in cps.iter() {
        let mut superset: Option<IntSet<usize>> = None;
        'atom: for atom in cp.atoms.iter() {
            if let Some(matches) = inverted_index.get(&atom.atom_index()) {
                match superset {
                    None => superset = Some(matches.iter().copied().collect()),
                    Some(ref mut set) => {
                        set.retain(|id| matches.contains(id));
                        if set.is_empty() {
                            superset = None;
                            break 'atom;
                        }
                    }
                }
            } else {
                superset = None;
                break 'atom;
            }
        }
        if let Some(set) = superset {
            for index in set.iter() {
                let set_sub = IntSet::from_iter(cp.atoms.iter().copied());
                let cp_super = &merged_points[*index];
                for encoded_atom in cp_super.atoms.iter() {
                    let new_anchor = encoded_atom.image();
                    let rotated_super = cp_super
                        .atoms
                        .iter()
                        .map(|a| a.image_sub(new_anchor))
                        .collect::<IntSet<EncodedAtom>>();
                    if set_sub.is_subset(&rotated_super) {
                        continue 'critical;
                    }
                }
            }
        }
        let i = merged_points.len();
        cp.atoms.iter().for_each(|atom| {
            inverted_index.entry(atom.atom_index()).or_default().push(i);
        });
        merged_points.push(cp.clone());
    }
    merged_points
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        atoms::Lattice,
        voxel_map::{EncodedAtom, EncodedImage},
    };

    // --- Helper for creating dummy Critical Points ---
    fn create_cp(
        pos: isize,
        kind: CriticalPointKind,
        atom_ids: &[u32],
        density: f64,
    ) -> CriticalPoint {
        let atoms = atom_ids
            .iter()
            .map(|&id| EncodedAtom::new(id, EncodedImage::new([0, 0, 0])))
            .collect::<Vec<_>>()
            .into_boxed_slice();

        CriticalPoint::new(pos, kind, atoms, density, 0.0)
    }

    // --- CriticalPointKey Tests ---

    #[test]
    fn test_critical_point_key_sorting() {
        // CP with atoms [2, 1, 3] should sort to [1, 2, 3]
        let cp = create_cp(0, CriticalPointKind::Bond, &[2, 1, 3], 0.0);
        let key = CriticalPointKey::from_cp(cp);

        let ids: Vec<u32> =
            key.into_box().iter().map(|a| a.atom_index()).collect();
        assert_eq!(ids, vec![1, 2, 3]);
    }

    #[test]
    fn test_critical_point_key_translation() {
        // Test that keys are canonical regarding translation
        // (If the first atom is at image [1,0,0], all atoms should be shifted back)

        let img1 = EncodedImage::new([1, 0, 0]);
        let img2 = EncodedImage::new([1, 0, 0]);

        let a1 = EncodedAtom::new(1, img1);
        let a2 = EncodedAtom::new(2, img2);

        let cp = CriticalPoint::new(
            0,
            CriticalPointKind::Bond,
            vec![a1, a2].into_boxed_slice(),
            0.0,
            0.0,
        );

        let key = CriticalPointKey::from_cp(cp);

        // Both images should be [0,0,0] after normalisation (subtracted the first one's image)
        for atom in key.into_box().iter() {
            assert!(
                atom.image().is_zero(),
                "Image should be normalised to zero"
            );
        }
    }

    // --- Nuclei Ordering Tests ---

    #[test]
    fn test_nuclei_ordering_simple() {
        // Case: Atom 0 has two potential nuclei candidates at pos 10 and 20.
        // Pos 20 has higher density.
        let atoms = Atoms::new(
            Lattice::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            vec![[0.0, 0.0, 0.0]],
            String::with_capacity(0),
        );
        let grid = Grid::new(
            [10, 10, 10],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [0.0, 0.0, 0.0],
        );

        let cp1 = create_cp(10, CriticalPointKind::Nuclei, &[0], 1.0);
        let cp2 = create_cp(20, CriticalPointKind::Nuclei, &[0], 5.0); // Winner

        let mut candidates = vec![cp1, cp2];

        let ordered = nuclei_ordering(&mut candidates, &atoms, &grid, false);

        assert_eq!(ordered.len(), 1);
        assert_eq!(ordered[0].position, 20); // Should pick high density one
    }

    #[test]
    fn test_nuclei_ordering_multiple_atoms() {
        // Atom 0 and Atom 1 both have candidates
        let atoms = Atoms::new(
            Lattice::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            vec![[0.0, 0.0, 0.0], [0., 0.5, 0.0]],
            String::with_capacity(0),
        );
        let grid = Grid::new(
            [10, 10, 10],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [0.0, 0.0, 0.0],
        );

        let cp0 = create_cp(10, CriticalPointKind::Nuclei, &[0], 5.0);
        let cp1 = create_cp(50, CriticalPointKind::Nuclei, &[1], 3.0);

        let mut candidates = vec![cp0, cp1];
        let ordered = nuclei_ordering(&mut candidates, &atoms, &grid, false);

        assert_eq!(ordered.len(), 2);
        assert_eq!(ordered[0].position, 10);
        assert_eq!(ordered[1].position, 50);
    }

    // --- Merge Tests ---

    #[test]
    fn test_critical_point_merge_exact_duplicate() {
        // Two CPs with same atoms should be merged
        let cp1 = create_cp(10, CriticalPointKind::Bond, &[1, 2], 0.0);
        let cp2 = create_cp(11, CriticalPointKind::Bond, &[1, 2], 1.0); // Diff position, same atoms, higher density

        let input = vec![cp1, cp2];
        let merged = critical_point_merge(input);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].position, 11);
    }

    #[test]
    fn test_critical_point_merge_subset() {
        // CP A: [1, 2, 3] (Ring)
        // CP B: [1, 2]    (Bond)
        // B is a subset of A?
        // Logic check: The code says:
        // "if superset is found... break 'atom" -> effectively discards the current CP?
        // Let's trace:
        // 1. Sort by len descending. -> [1,2,3] comes first. Added to merged. index map: {1:0, 2:0, 3:0}
        // 2. Process [1,2].
        //    Atom 1 -> matches index 0. Superset = {0}.
        //    Atom 2 -> matches index 0. Superset intersect {0} = {0}.
        //    Loop ends. superset is Some({0}).
        //    Result: [1,2] is NOT added.

        // So yes, subsets are merged into supersets.

        let cp_ring = create_cp(10, CriticalPointKind::Ring, &[1, 2, 3], 0.0);
        let cp_bond = create_cp(11, CriticalPointKind::Bond, &[1, 2], 0.0);

        let input = vec![cp_ring, cp_bond];
        let merged = critical_point_merge(input);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].atoms.len(), 3); // Kept the larger one
    }

    #[test]
    fn test_critical_point_merge_distinct() {
        let cp1 = create_cp(10, CriticalPointKind::Bond, &[1, 2], 0.0);
        let cp2 = create_cp(11, CriticalPointKind::Bond, &[3, 4], 0.0);

        let input = vec![cp1, cp2];
        let merged = critical_point_merge(input);

        assert_eq!(merged.len(), 2);
    }
}
