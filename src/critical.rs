use crate::atoms::Atoms;
use crate::grid::Grid;
use crate::hash::{IntMap, IntSet};
use crate::progress::progress_bar;
use crate::threading::parallel_prune;
use crate::voxel_map::EncodedAtom;

/// Represents a critical point in the charge density topology.
#[derive(Clone, Debug)]
pub struct CriticalPoint {
    /// The 1D index of the voxel where this point is located.
    pub position: isize,
    /// The topological classification (Nuclei, Bond, Ring, Cage).
    pub kind: CriticalPointKind,
    /// The atoms associated with this point (e.g., the two atoms a bond connects).
    pub atoms: Box<[EncodedAtom]>,
}

impl CriticalPoint {
    /// Create a new critical point.
    pub fn new(
        position: isize,
        kind: CriticalPointKind,
        atoms: Box<[EncodedAtom]>,
    ) -> Self {
        CriticalPoint {
            position,
            kind,
            atoms,
        }
    }
}

/// The pruned critical points of a system, grouped by kind.
#[derive(Debug, Default)]
pub struct CriticalPoints {
    /// One per atom, see [`nuclei_ordering`].
    pub nuclei: Vec<CriticalPoint>,
    /// The (3, -1) saddles, see [`bond_pruning`].
    pub bonds: Vec<CriticalPoint>,
    /// The (3, +1) saddles, see [`ring_pruning`].
    pub rings: Vec<CriticalPoint>,
    /// The (3, +3) minima, see [`cage_pruning`].
    pub cages: Vec<CriticalPoint>,
}

/// The topological classification of a [`CriticalPoint`].
#[derive(Eq, Ord, PartialEq, PartialOrd, Debug, Clone, Copy)]
pub enum CriticalPointKind {
    /// A (3, -3) maximum in the density.
    Nuclei,
    /// A (3, -1) saddle between two atoms.
    Bond,
    /// A (3, +1) saddle in a ring of atoms.
    Ring,
    /// A (3, +3) minimum enclosed by atoms.
    Cage,
    /// A placeholder for an atom with no critical point.
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
    /// Builds the canonical key for a list of atoms.
    pub fn new(atoms: &[EncodedAtom]) -> Self {
        let mut atoms = atoms.to_vec();
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

    /// The canonically ordered atoms.
    pub fn into_box(self) -> Box<[EncodedAtom]> {
        self.0.into()
    }
}

/// Selects the true maximum of each atom and relates the others to it.
///
/// In the raw output, a single atom might be associated with multiple "maxima" candidates due
/// to grid discretisation noise or genuine secondary maxima. This function groups candidates
/// by their assigned atom and selects the one with the highest charge density as the true
/// nucleus position. The remaining candidates of that atom are re-encoded, in place, with the
/// periodic image of the atom they are closest to so that an atom can later be found bonded
/// to its own image.
///
/// # Arguments
/// * `nuclei`: Candidate points classified as [`CriticalPointKind::Nuclei`], all in image
///   `[0, 0, 0]`. The images of the non-selected candidates are updated.
/// * `density`: The charge density array (used to compare peak heights).
/// * `atoms`: The atoms of the system.
/// * `grid`: The grid the density is stored on.
///
/// # Returns
/// A vector with one point per atom, where index `i` is the true maximum of atom `i` or a
/// [`CriticalPointKind::Blank`] point if the atom has no maximum.
pub fn nuclei_ordering(
    nuclei: &mut [CriticalPoint],
    density: &[f64],
    atoms: &Atoms,
    grid: &Grid,
) -> Vec<CriticalPoint> {
    let mut candidates = vec![Vec::<usize>::new(); atoms.positions.len()];
    for (i, cp) in nuclei.iter().enumerate() {
        candidates[cp.atoms[0].atom_index() as usize].push(i);
    }
    candidates
        .into_iter()
        .map(|indices| {
            let Some(&best) = indices.iter().max_by(|&&a, &&b| {
                density[nuclei[a].position as usize]
                    .total_cmp(&density[nuclei[b].position as usize])
            }) else {
                return CriticalPoint::new(
                    0,
                    CriticalPointKind::Blank,
                    Box::new([]),
                );
            };
            let true_maximum = nuclei[best].clone();
            let true_position = grid.to_cartesian(true_maximum.position);
            for &i in indices.iter().filter(|&&i| i != best) {
                let position = grid.to_cartesian(nuclei[i].position);
                let image =
                    atoms.lattice.closest_image(true_position, position);
                nuclei[i].atoms[0] = nuclei[i].atoms[0].image_sub(image);
            }
            true_maximum
        })
        .collect()
}

/// Validates and deduplicates critical points with a labelled progress bar.
fn prune(
    critical_points: &[CriticalPoint],
    density: &[f64],
    validator: impl Fn(&CriticalPoint) -> bool + Sync,
    label: &str,
    threads: usize,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let pbar = progress_bar(visible_bar, critical_points.len(), label);
    parallel_prune(critical_points, density, validator, threads, pbar)
}

/// Deduplicates Bond Critical Points (3, -1).
///
/// Bonds are defined by the two atoms they connect. In the raw analysis, multiple voxels
/// along the boundary between two atoms might be flagged as saddle points. Candidates are
/// grouped by their atom pair with [`parallel_prune`] and only the point with the **highest
/// charge density** is retained, this corresponds to the true saddle point on the gradient
/// path.
///
/// # Arguments
/// * `bonds`: The raw list of candidate bond points.
/// * `density`: The charge density grid.
/// * `threads`: Number of threads to use.
/// * `visible_bar`: Whether to display a progress bar.
pub fn bond_pruning(
    bonds: &[CriticalPoint],
    density: &[f64],
    threads: usize,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    prune(
        bonds,
        density,
        |_| true,
        "Pruning Bond Critical Points",
        threads,
        visible_bar,
    )
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
    for bond in bonds {
        let [a, b] = [bond.atoms[0], bond.atoms[1]];
        adjacency[a.atom_index() as usize].push(b.image_sub(a.image()));
        adjacency[b.atom_index() as usize].push(a.image_sub(b.image()));
    }
    adjacency
}

/// Filters Ring Critical Points (3, +1) by validating bond connectivity.
///
/// This function verifies that the atoms associated with a candidate Ring Critical Point
/// form a valid, continuous closed loop (cycle) within the determined bond network.
///
/// # Logic
/// 1. **Valency Validation**: Verifies that every atom in the candidate is connected
///    to exactly two other atoms *within the same candidate set* (via `bond_adjacency`).
/// 2. **Traversal Check**: Ensures the atoms form a single continuous loop (no disjoint parts).
/// 3. **Deduplication**: Groups by unique atom list and keeps the candidate with the highest density.
///
/// # Arguments
/// * `rings`: The raw list of candidate ring points.
/// * `bond_adjacency`: The adjacency graph derived from Bond Critical Points.
/// * `density`: The charge density grid.
/// * `threads`: Number of threads to use.
/// * `visible_bar`: Whether to display a progress bar.
pub fn ring_pruning(
    rings: &[CriticalPoint],
    bond_adjacency: &[Vec<EncodedAtom>],
    density: &[f64],
    threads: usize,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let is_ring = |cp: &CriticalPoint| {
        // walk the ring, at each atom the two bonds inside the candidate
        // must be to the atom we came from and the next atom
        let mut previous = 0;
        let mut current = 0;
        let mut steps = 1;
        loop {
            let graph =
                &bond_adjacency[cp.atoms[current].atom_index() as usize];
            let anchor = cp.atoms[current].image();
            let mut bonded =
                cp.atoms.iter().enumerate().filter(|(_, atom)| {
                    graph.contains(&atom.image_sub(anchor))
                });
            let next = match (bonded.next(), bonded.next(), bonded.next()) {
                (Some(a), Some(b), None) => {
                    if a.0 == previous {
                        b.0
                    } else {
                        a.0
                    }
                }
                // anything other than two bonds within the ring is not a ring
                _ => return false,
            };
            // back at the start, the steps limit only guards against an infinite loop
            if next == 0 || steps > cp.atoms.len() {
                break;
            }
            steps += 1;
            previous = current;
            current = next;
        }
        // a single loop visits every atom
        steps == cp.atoms.len()
    };
    prune(
        rings,
        density,
        is_ring,
        "Pruning Ring Critical Points",
        threads,
        visible_bar,
    )
}

/// Deduplicates Cage Critical Points (3, +3).
///
/// Candidates are grouped by their atom list and only the point with the highest
/// density is kept.
///
/// # Arguments
/// * `cages`: The raw list of candidate cage points.
/// * `density`: The charge density grid.
/// * `threads`: Number of threads to use.
/// * `visible_bar`: Whether to display a progress bar.
pub fn cage_pruning(
    cages: &[CriticalPoint],
    density: &[f64],
    threads: usize,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    prune(
        cages,
        density,
        |_| true,
        "Pruning Cage Critical Points",
        threads,
        visible_bar,
    )
}

/// Merges degenerate or subset critical points.
///
/// This simplifies the topology by removing redundant critical points.
///
/// # Logic
/// 1. Sorts critical points by the number of associated atoms (descending).
/// 2. Iterates through them. If a point's atoms are a **subset** of an already kept
///    point's atoms (e.g., a Ring \[1,2,3\] is likely charge fluctuation before the Ring
///    \[1,2,3,4\]), under any translation of the periodic images, it is discarded.
///
/// # Returns
/// A cleaner list of unique critical points.
pub fn critical_point_merge(mut cps: Vec<CriticalPoint>) -> Vec<CriticalPoint> {
    cps.sort_unstable_by_key(|a| std::cmp::Reverse(a.atoms.len()));
    let mut merged_points: Vec<CriticalPoint> = Vec::with_capacity(cps.len());
    // the kept points that contain each atom index
    let mut inverted_index: IntMap<u32, Vec<usize>> = IntMap::default();
    for cp in cps {
        // the kept points containing every atom index of this point
        let mut candidates = cp.atoms.iter().map(|atom| {
            inverted_index
                .get(&atom.atom_index())
                .map(|v| v.as_slice())
                .unwrap_or(&[])
        });
        let mut supersets: IntSet<usize> =
            candidates.next().unwrap_or(&[]).iter().copied().collect();
        for matches in candidates {
            supersets.retain(|i| matches.contains(i));
        }
        // the point is redundant if its atoms are a subset of a candidate's
        // atoms under some translation of the periodic images
        let atoms: IntSet<EncodedAtom> = cp.atoms.iter().copied().collect();
        let is_subset = supersets.iter().any(|&i| {
            let superset = &merged_points[i].atoms;
            superset.iter().any(|anchor| {
                let image = anchor.image();
                let shifted: IntSet<EncodedAtom> =
                    superset.iter().map(|a| a.image_sub(image)).collect();
                atoms.is_subset(&shifted)
            })
        });
        if is_subset {
            continue;
        }
        for atom in cp.atoms.iter() {
            inverted_index
                .entry(atom.atom_index())
                .or_default()
                .push(merged_points.len());
        }
        merged_points.push(cp);
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
    ) -> CriticalPoint {
        let atoms = atom_ids
            .iter()
            .map(|&id| EncodedAtom::new(id, EncodedImage::new([0, 0, 0])))
            .collect::<Vec<_>>()
            .into_boxed_slice();

        CriticalPoint::new(pos, kind, atoms)
    }

    // --- CriticalPointKey Tests ---

    #[test]
    fn test_critical_point_key_sorting() {
        // CP with atoms [2, 1, 3] should sort to [1, 2, 3]
        let cp = create_cp(0, CriticalPointKind::Bond, &[2, 1, 3]);
        let key = CriticalPointKey::new(&cp.atoms);

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
        );

        let key = CriticalPointKey::new(&cp.atoms);

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

        let mut density = vec![0.0; 30];
        density[10] = 1.0;
        density[20] = 5.0; // Winner

        let cp1 = create_cp(10, CriticalPointKind::Nuclei, &[0]);
        let cp2 = create_cp(20, CriticalPointKind::Nuclei, &[0]);

        let mut candidates = vec![cp1, cp2];

        let ordered = nuclei_ordering(&mut candidates, &density, &atoms, &grid);

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
        let mut density = vec![0.0; 100];
        density[10] = 5.0; // Atom 0
        density[50] = 3.0; // Atom 1

        let cp0 = create_cp(10, CriticalPointKind::Nuclei, &[0]);
        let cp1 = create_cp(50, CriticalPointKind::Nuclei, &[1]);

        let mut candidates = vec![cp0, cp1];
        let ordered = nuclei_ordering(&mut candidates, &density, &atoms, &grid);

        assert_eq!(ordered.len(), 2);
        assert_eq!(ordered[0].position, 10);
        assert_eq!(ordered[1].position, 50);
    }

    // --- Merge Tests ---

    #[test]
    fn test_critical_point_merge_exact_duplicate() {
        // Two CPs with same atoms should be merged
        let cp1 = create_cp(10, CriticalPointKind::Bond, &[1, 2]);
        let cp2 = create_cp(11, CriticalPointKind::Bond, &[1, 2]); // Diff position, same atoms

        let input = vec![cp1, cp2];
        let merged = critical_point_merge(input);

        assert_eq!(merged.len(), 1);
    }

    #[test]
    fn test_critical_point_merge_subset() {
        // CP A: [1, 2, 3] (Ring)
        // CP B: [1, 2]    (Bond)
        // B is a subset of A so only A is kept.
        let cp_ring = create_cp(10, CriticalPointKind::Ring, &[1, 2, 3]);
        let cp_bond = create_cp(11, CriticalPointKind::Bond, &[1, 2]);

        let input = vec![cp_ring, cp_bond];
        let merged = critical_point_merge(input);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].atoms.len(), 3); // Kept the larger one
    }

    #[test]
    fn test_critical_point_merge_translated_subset() {
        // [1(0,0,0), 2(1,0,0)] is a subset of [1(-1,0,0), 2(0,0,0), 3(0,0,0)]
        // once the images are shifted by the first atom
        let superset = CriticalPoint::new(
            10,
            CriticalPointKind::Ring,
            vec![
                EncodedAtom::new(1, EncodedImage::new([-1, 0, 0])),
                EncodedAtom::new_zero_image(2),
                EncodedAtom::new_zero_image(3),
            ]
            .into(),
        );
        let subset = CriticalPoint::new(
            11,
            CriticalPointKind::Bond,
            vec![
                EncodedAtom::new_zero_image(1),
                EncodedAtom::new(2, EncodedImage::new([1, 0, 0])),
            ]
            .into(),
        );
        let mut merged = critical_point_merge(vec![subset.clone(), superset]);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].atoms.len(), 3);
        // but a different image relationship is a distinct point
        let other = CriticalPoint::new(
            11,
            CriticalPointKind::Bond,
            vec![
                EncodedAtom::new_zero_image(1),
                EncodedAtom::new(2, EncodedImage::new([0, 1, 0])),
            ]
            .into(),
        );
        let merged = critical_point_merge(vec![other, merged.remove(0)]);
        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn test_critical_point_merge_distinct() {
        let cp1 = create_cp(10, CriticalPointKind::Bond, &[1, 2]);
        let cp2 = create_cp(11, CriticalPointKind::Bond, &[3, 4]);

        let input = vec![cp1, cp2];
        let merged = critical_point_merge(input);

        assert_eq!(merged.len(), 2);
    }

    // --- Adjacency and Pruning Tests ---

    #[test]
    fn test_bond_adjacency_images() {
        let bond = CriticalPoint::new(
            5,
            CriticalPointKind::Bond,
            vec![
                EncodedAtom::new_zero_image(0),
                EncodedAtom::new(1, EncodedImage::new([1, 0, 0])),
            ]
            .into(),
        );
        let adjacency = bond_adjacency(&[bond], 3);
        assert_eq!(
            adjacency[0],
            vec![EncodedAtom::new(1, EncodedImage::new([1, 0, 0]))]
        );
        assert_eq!(
            adjacency[1],
            vec![EncodedAtom::new(0, EncodedImage::new([-1, 0, 0]))]
        );
        assert!(adjacency[2].is_empty());
    }

    #[test]
    fn test_bond_pruning_keeps_highest_density() {
        let mut density = vec![0.0; 30];
        density[10] = 1.0;
        density[20] = 5.0;
        let bonds = vec![
            create_cp(10, CriticalPointKind::Bond, &[1, 2]),
            create_cp(20, CriticalPointKind::Bond, &[2, 1]),
            create_cp(15, CriticalPointKind::Bond, &[2, 3]),
        ];
        let mut pruned = bond_pruning(&bonds, &density, 2, false);
        pruned.sort_unstable_by_key(|cp| cp.position);
        assert_eq!(pruned.len(), 2);
        assert_eq!(pruned[0].position, 15);
        assert_eq!(pruned[1].position, 20);
        let cages = cage_pruning(&bonds, &density, 2, false);
        assert_eq!(cages.len(), 2);
    }

    #[test]
    fn test_ring_pruning_requires_cycle() {
        let density = vec![1.0; 30];
        // bonds 0-1, 1-2, 2-0 form a triangle, 3 hangs off atom 0
        let bonds = vec![
            create_cp(1, CriticalPointKind::Bond, &[0, 1]),
            create_cp(2, CriticalPointKind::Bond, &[1, 2]),
            create_cp(3, CriticalPointKind::Bond, &[2, 0]),
            create_cp(4, CriticalPointKind::Bond, &[0, 3]),
        ];
        let adjacency = bond_adjacency(&bonds, 4);
        let rings = vec![
            create_cp(10, CriticalPointKind::Ring, &[0, 1, 2]), // a cycle
            create_cp(11, CriticalPointKind::Ring, &[0, 1, 3]), // 1-3 not bonded
            create_cp(12, CriticalPointKind::Ring, &[0, 1, 2, 3]), // 0 has 3 bonds
        ];
        let pruned = ring_pruning(&rings, &adjacency, &density, 2, false);
        assert_eq!(pruned.len(), 1);
        assert_eq!(pruned[0].position, 10);
    }

    #[test]
    fn test_nuclei_ordering_secondary_maximum_image() {
        // atom 0 has a secondary maximum across the periodic boundary which
        // should be re-encoded with the image of the true maximum
        let atoms = Atoms::new(
            Lattice::new([
                [10.0, 0.0, 0.0],
                [0.0, 10.0, 0.0],
                [0.0, 0.0, 10.0],
            ]),
            vec![[0.0, 0.0, 0.0]],
            String::new(),
        );
        let grid = Grid::new(
            [10, 10, 10],
            [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            [0.0; 3],
        );
        let mut density = vec![0.0; 1000];
        density[0] = 5.0; // true maximum at the origin
        density[900] = 1.0; // x = 9, closest to the origin's image at x = 10
        let mut nuclei = vec![
            create_cp(900, CriticalPointKind::Nuclei, &[0]),
            create_cp(0, CriticalPointKind::Nuclei, &[0]),
        ];
        let ordered = nuclei_ordering(&mut nuclei, &density, &atoms, &grid);
        assert_eq!(ordered[0].position, 0);
        assert!(ordered[0].atoms[0].image().is_zero());
        // the secondary maximum belongs to the image of the atom at x = 10
        assert_eq!(nuclei[0].atoms[0].image().decode(), [1, 0, 0]);
        assert!(nuclei[1].atoms[0].image().is_zero());
    }
}
