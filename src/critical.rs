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
}

impl CriticalPoint {
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
    nuclei: Vec<CriticalPoint>,
    density: &[f64],
    atom_len: usize,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            nuclei.len(),
            String::from("Pruning Nucleus Critical Points"),
        )),
    };
    let pbar = &progress_bar;
    // Find all the nuclei with the same atom number and group them based on which is largest
    // charge density. All maxima will be in image (0, 0, 0).
    // We need to order the nuclei by atom so that we can get the position of them by knowing the
    // atom number.
    let mut ordered_nuclei =
        vec![
            CriticalPoint::new(0, CriticalPointKind::Blank, Box::new([]));
            atom_len
        ];
    nuclei.iter().for_each(|cp| {
        let p = cp.position;
        let rho = density[p as usize];
        let atom_num = cp.atoms[0].atom_index() as usize;
        if let CriticalPointKind::Blank = ordered_nuclei[atom_num].kind {
            ordered_nuclei[atom_num] =
                CriticalPoint::new(cp.position, cp.kind, cp.atoms.clone());
        } else if rho > density[ordered_nuclei[atom_num].position as usize] {
            ordered_nuclei[atom_num] =
                CriticalPoint::new(cp.position, cp.kind, cp.atoms.clone());
        }
        pbar.tick();
    });
    ordered_nuclei
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
/// * `density`: The charge density grid.
/// * `threads`: Number of threads to use.
/// * `visible_bar`: Whether to display a progress bar.
pub fn bond_pruning(
    bonds: &[CriticalPoint],
    density: &[f64],
    threads: usize,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            bonds.len(),
            String::from("Pruning Bond Critical Points"),
        )),
    };
    parallel_prune(bonds, density, |_| true, threads, progress_bar)
}

/// Filters Ring Critical Points (3, +1) and enforces planarity.
///
/// A Ring Critical Point must connect at least 3 atoms and those atoms must lie approximately
/// on a single plane. Points that fail this geometric check are discarded.
///
/// # Logic
/// 1. **Size Check**: Must have $\ge$ 3 atoms.
/// 2. **Planarity Check**: Calculates the normal vector of the plane formed by the first 3 atoms.
///    Then verifies that all subsequent atoms lie on this plane (tolerance ~5.7°).
///    - If atoms are **not** coplanar, the point is rejected.
/// 3. **Deduplication**: Groups by atom list and keeps the candidate with the highest density.
///
/// # Arguments
/// * `rings`: The raw list of candidate ring points.
/// * `ordered_nuclei`: The list of definitive nucleus positions (used to get atom coordinates).
/// * `density`: The charge density grid.
/// * `atoms`: The system geometry (lattice/positions).
/// * `grid`: The voxel grid (for coordinate conversion).
pub fn ring_pruning(
    rings: &[CriticalPoint],
    ordered_nuclei: &[CriticalPoint],
    density: &[f64],
    atoms: &Atoms,
    grid: &Grid,
    threads: usize,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            rings.len(),
            String::from("Pruning Ring Critical Points"),
        )),
    };
    parallel_prune(
        rings,
        density,
        |cp| {
            if cp.atoms.len() < 3 {
                return false;
            }
            let positions: Vec<[f64; 3]> = cp.atoms[..3]
                .iter()
                .map(|encoded_atom| {
                    let (atom_num, encoded_image) =
                        encoded_atom.decode_partial();
                    let mut position = grid.to_cartesian(
                        ordered_nuclei[atom_num as usize].position,
                    );
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
                let vec_3 = subtract(position, positions[0]);
                let mut plane_t = cross(vec_1, vec_3);
                let plane_normal = norm(plane_t);
                plane_t.iter_mut().for_each(|f| *f /= plane_normal);
                // TODO: make this a tolerance currently 5.73 degrees
                if vdot(plane, plane_t).abs() < 0.995 {
                    return false;
                }
            }
            true
        },
        threads,
        progress_bar,
    )
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
/// * `density`: The charge density grid.
/// * `atoms`: The system geometry.
/// * `grid`: The voxel grid.
pub fn cage_pruning(
    cages: &[CriticalPoint],
    ordered_nuclei: &[CriticalPoint],
    density: &[f64],
    atoms: &Atoms,
    grid: &Grid,
    threads: usize,
    visible_bar: bool,
) -> Vec<CriticalPoint> {
    let progress_bar: Box<dyn ProgressBar> = match visible_bar {
        false => Box::new(HiddenBar {}),
        true => Box::new(Bar::new(
            cages.len(),
            String::from("Pruning Cage Critical Points"),
        )),
    };
    parallel_prune(
        cages,
        density,
        |cp| {
            if cp.atoms.len() < 4 {
                return false;
            }
            let positions: Vec<[f64; 3]> = cp.atoms[..3]
                .iter()
                .map(|encoded_atom| {
                    let (atom_num, encoded_image) =
                        encoded_atom.decode_partial();
                    let mut position = grid.to_cartesian(
                        ordered_nuclei[atom_num as usize].position,
                    );
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
        },
        threads,
        progress_bar,
    )
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
    cps.sort_unstable_by(|a, b| b.atoms.len().cmp(&a.atoms.len()));
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
    use crate::voxel_map::{EncodedAtom, EncodedImage};

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

        let mut density = vec![0.0; 30];
        density[10] = 1.0;
        density[20] = 5.0; // Winner

        let cp1 = create_cp(10, CriticalPointKind::Nuclei, &[0]);
        let cp2 = create_cp(20, CriticalPointKind::Nuclei, &[0]);

        let candidates = vec![cp1, cp2];

        let ordered = nuclei_ordering(candidates, &density, 1, false);

        assert_eq!(ordered.len(), 1);
        assert_eq!(ordered[0].position, 20); // Should pick high density one
    }

    #[test]
    fn test_nuclei_ordering_multiple_atoms() {
        // Atom 0 and Atom 1 both have candidates
        let mut density = vec![0.0; 100];
        density[10] = 5.0; // Atom 0
        density[50] = 3.0; // Atom 1

        let cp0 = create_cp(10, CriticalPointKind::Nuclei, &[0]);
        let cp1 = create_cp(50, CriticalPointKind::Nuclei, &[1]);

        let candidates = vec![cp0, cp1];
        let ordered = nuclei_ordering(candidates, &density, 2, false);

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

        let cp_ring = create_cp(10, CriticalPointKind::Ring, &[1, 2, 3]);
        let cp_bond = create_cp(11, CriticalPointKind::Bond, &[1, 2]);

        let input = vec![cp_ring, cp_bond];
        let merged = critical_point_merge(input);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].atoms.len(), 3); // Kept the larger one
    }

    #[test]
    fn test_critical_point_merge_distinct() {
        let cp1 = create_cp(10, CriticalPointKind::Bond, &[1, 2]);
        let cp2 = create_cp(11, CriticalPointKind::Bond, &[3, 4]);

        let input = vec![cp1, cp2];
        let merged = critical_point_merge(input);

        assert_eq!(merged.len(), 2);
    }
}
