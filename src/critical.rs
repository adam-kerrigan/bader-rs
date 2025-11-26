use rustc_hash::{FxHashMap, FxHashSet};

use crate::atoms::Atoms;
use crate::grid::Grid;
use crate::progress::{Bar, HiddenBar, ProgressBar};
use crate::threading::parallel_prune;
use crate::utils::{cross, norm, subtract, vdot};
use crate::voxel_map::EncodedAtom;

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

pub fn critical_point_merge(mut cps: Vec<CriticalPoint>) -> Vec<CriticalPoint> {
    cps.sort_unstable_by(|a, b| b.atoms.len().cmp(&a.atoms.len()));
    let mut merged_points: Vec<CriticalPoint> = Vec::with_capacity(cps.len());
    let mut inverted_index: FxHashMap<u32, Vec<usize>> = FxHashMap::default();
    'critical: for cp in cps.iter() {
        let mut superset: Option<FxHashSet<usize>> = None;
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
                let set_sub = FxHashSet::from_iter(cp.atoms.iter().copied());
                let cp_super = &merged_points[*index];
                for encoded_atom in cp_super.atoms.iter() {
                    let new_anchor = encoded_atom.image();
                    let rotated_super = cp_super
                        .atoms
                        .iter()
                        .map(|a| a.image_sub(new_anchor))
                        .collect::<FxHashSet<EncodedAtom>>();
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
