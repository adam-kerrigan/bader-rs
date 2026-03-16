use crate::critical::CriticalPoint;
use crate::grid::Grid;
use crate::io::{FileFormat, FileType};
use crate::voxel_map::EncodedAtom;
use core::fmt;
use std::fs::File;
use std::io::Write;

type CriticalPoints = (
    Vec<CriticalPoint>,
    Vec<CriticalPoint>,
    Vec<CriticalPoint>,
    Vec<CriticalPoint>,
);

type CriticalPointInfo = ([f64; 3], Box<[EncodedAtom]>, f64, f64);

pub struct CriticalPointOutput {
    splash: String,
    positions: Vec<[f64; 3]>,
    atom_nucleus_map: Vec<Vec<CriticalPointInfo>>,
    atom_bond_map: Vec<Vec<CriticalPointInfo>>,
    atom_ring_map: Vec<Vec<CriticalPointInfo>>,
    atom_cage_map: Vec<Vec<CriticalPointInfo>>,
}

impl CriticalPointOutput {
    pub fn new(
        positions: Vec<[f64; 3]>,
        critical_points: CriticalPoints,
        grid: &Grid,
        file_type: FileType,
    ) -> Self {
        let atom_num = positions.len();
        let (nuclei, bonds, rings, cages) = critical_points;
        let pivot_critical_points =
            |cps: &[CriticalPoint]| -> Vec<Vec<CriticalPointInfo>> {
                let mut pivot = vec![vec![]; atom_num];
                cps.iter().for_each(|cp| {
                    'pivot: for (i, encoded_atom) in cp.atoms.iter().enumerate()
                    {
                        let atom_index = encoded_atom.atom_index();
                        let image = encoded_atom.image();
                        for cp_check in cp.atoms[..i].iter() {
                            if atom_index == cp_check.atom_index() {
                                continue 'pivot;
                            }
                        }
                        pivot[atom_index as usize].push((
                            file_type.coordinate_format(
                                grid.to_cartesian(cp.position),
                            ),
                            cp.atoms
                                .iter()
                                .map(|encoded_atom| {
                                    encoded_atom.image_sub(image)
                                })
                                .collect(),
                            cp.density,
                            cp.laplacian,
                        ));
                    }
                });
                pivot
            };
        let positions = positions
            .iter()
            .map(|p| file_type.coordinate_format([p[0], p[1], p[2]]))
            .collect();
        let atom_nucleus_map = pivot_critical_points(&nuclei);
        let atom_bond_map = pivot_critical_points(&bonds);
        let atom_ring_map = pivot_critical_points(&rings);
        let atom_cage_map = pivot_critical_points(&cages);
        let splash = format!(
            "## Topological Information\n\n * Total nuclei: {}\n * Total bonds: {}\n * Total rings: {}\n * Total cages: {}\n",
            nuclei.len(),
            bonds.len(),
            rings.len(),
            cages.len()
        );
        Self {
            splash,
            positions,
            atom_nucleus_map,
            atom_bond_map,
            atom_ring_map,
            atom_cage_map,
        }
    }
}

impl fmt::Display for CriticalPointOutput {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut output = self.splash.clone();
        self.positions.iter().enumerate().for_each(|(i, position)| {
            output.push_str(&format!(
                "\n### Atom: {}\n * Position: {:.6} {:.6} {:.6}\n * Coordination: {}\n * Critical Points:\n",
                i + 1, position[0], position[1], position[2], self.atom_bond_map[i].len(),
            ));
            let self_atom = EncodedAtom::new_zero_image(i as u32);
            let mut crit_type= Vec::<String>::with_capacity(self.atom_nucleus_map[i].len() + self.atom_bond_map[i].len() + self.atom_ring_map[i].len() + self.atom_cage_map[i].len());
            let mut max_type: usize = 4;
            let mut crit_pos= Vec::<Vec<String>>::with_capacity(crit_type.capacity());
            let mut max_pos: usize = 8;
            let mut crit_den= Vec::<String>::with_capacity(crit_type.capacity());
            let mut max_den: usize = 6;
            let mut crit_lap= Vec::<String>::with_capacity(crit_type.capacity());
            let mut max_lap: usize = 8;
            let mut crit_mem= Vec::<String>::with_capacity(crit_type.capacity());
            let mut max_mem: usize = 7;
            self.atom_nucleus_map[i].iter().for_each(|(position, _, density, laplacian)| {
                crit_type.push(String::from("Nucleus"));
                max_type = max_type.max(7);
                crit_pos.push(Vec::from_iter(position.iter().map(|p| {
                    let f = format!("{:.6}", p);
                    max_pos = max_pos.max(f.len());
                    f
                })));
                let d = format!("{:.6}", density);
                max_den = max_den.max(d.len());
                crit_den.push(d);
                let l = format!("{:.6}", laplacian);
                max_lap = max_lap.max(l.len());
                crit_lap.push(l);
                crit_mem.push(String::with_capacity(0));
            });
            self.atom_bond_map[i].iter().for_each(|(position, atoms, density, laplacian)| {
                if let Some(atom) = atoms.iter().find(|a| a != &&self_atom) {
                    crit_type.push(String::from("Bond"));
                    crit_pos.push(Vec::from_iter(position.iter().map(|p| {
                        let f = format!("{:.6}", p);
                        max_pos = max_pos.max(f.len());
                        f
                    })));
                    let d = format!("{:.6}", density);
                    max_den = max_den.max(d.len());
                    crit_den.push(d);
                    let l = format!("{:.6}", laplacian);
                    max_lap = max_lap.max(l.len());
                    crit_lap.push(l);
                    let m = match atom.image().is_zero() {
                        true => format!("{}", atom.atom_index()),
                        false => {
                            let other_image = atom.image().decode();
                            format!("{}({} {} {})", atom.atom_index(), other_image[0], other_image[1], other_image[2])
                        }
                    };
                    max_mem = max_mem.max(m.len());
                    crit_mem.push(m);

                }
            });
            self.atom_ring_map[i].iter().for_each(|(position, atoms, density, laplacian)| {
                crit_type.push(String::from("Ring"));
                crit_pos.push(Vec::from_iter(position.iter().map(|p| {
                    let f = format!("{:.6}", p);
                    max_pos = max_pos.max(f.len());
                    f
                })));
                let d = format!("{:.6}", density);
                max_den = max_den.max(d.len());
                crit_den.push(d);
                let l = format!("{:.6}", laplacian);
                max_lap = max_lap.max(l.len());
                crit_lap.push(l);
                let mut members = String::from("");
                atoms.iter().for_each(|atom| {
                    if atom != &self_atom {
                        let m = match atom.image().is_zero() {
                            true => format!("{}, ", atom.atom_index()),
                            false => {
                                let other_image = atom.image().decode();
                                format!("{}({} {} {}), ", atom.atom_index(), other_image[0], other_image[1], other_image[2])
                            }
                        };
                        members.push_str(&m);
                    }
                });
                members.pop();
                members.pop();
                max_mem = max_mem.max(members.len());
                crit_mem.push(members);
            });
            self.atom_cage_map[i].iter().for_each(|(position, atoms, density, laplacian)| {
                crit_type.push(String::from("Cage"));
                crit_pos.push(Vec::from_iter(position.iter().map(|p| {
                    let f = format!("{:.6}", p);
                    max_pos = max_pos.max(f.len());
                    f
                })));
                let d = format!("{:.6}", density);
                max_den = max_den.max(d.len());
                crit_den.push(d);
                let l = format!("{:.6}", laplacian);
                max_lap = max_lap.max(l.len());
                crit_lap.push(l);
                let mut members = String::from("");
                atoms.iter().for_each(|atom| {
                    if atom != &self_atom {
                        let m = match atom.image().is_zero() {
                            true => format!("{}, ", atom.atom_index()),
                            false => {
                                let other_image = atom.image().decode();
                                format!("{}({} {} {}), ", atom.atom_index(), other_image[0], other_image[1], other_image[2])
                            }
                        };
                        members.push_str(&m);
                    }
                });
                members.pop();
                members.pop();
                max_mem = max_mem.max(members.len());
                crit_mem.push(members);
            });
            let max_title_pos = max_pos * 3 + 2;
            output.push_str(&format!("| {:^max_type$} | {:^max_title_pos$} | {:^max_den$} | {:^max_lap$} | {:^max_mem$} |\n", "Type", "Position", "Density", "Laplacian", "Members"));
            output.push_str(&format!("|-{:-^max_type$}-|-{:-^max_title_pos$}-|-{:-^max_den$}-|-{:-^max_lap$}-|-{:-^max_mem$}-|\n", "-", "-", "-", "-", "-"));
            crit_type.into_iter().zip(crit_pos).zip(crit_den).zip(crit_lap).zip(crit_mem).for_each(|((((t, p), d), l), m)| {
                output.push_str(&format!("| {:^max_type$} | {:>max_pos$} {:>max_pos$} {:>max_pos$} | {:>max_den$} | {:>max_lap$} | {:<max_mem$} |\n", t, p[0], p[1], p[2], d, l, m));
            });
        });
        write!(f, "{}", output)
    }
}

/// Structure that contains and builds the table.
pub struct PartitionTable {
    /// How wide each column is.
    column_width: Vec<usize>,
    /// The number of charge and spin densities.
    density_num: usize,
    /// The rows of the table as a vector of strings.
    rows: Vec<Vec<String>>,
    footer: String,
}

impl PartitionTable {
    /// Creates a new structure and sets the minimum widths of each.
    pub fn new(
        partitioned_density: &[Box<[f64]>],
        partitioned_volume: &[f64],
        radius: &[f64],
        errors: &[f64],
        boundary_voxels: usize,
        total_voxels: usize,
    ) -> Self {
        let density_num = partitioned_density[0].len();
        let mut rows: Vec<Vec<String>> =
            vec![Vec::with_capacity(4 + density_num)];
        let mut column_width: Vec<usize> = vec![6; 4 + density_num];
        column_width[2 + density_num] = 8;
        // calculate the total density for each density supplied
        let total_density: Vec<f64> = partitioned_density.iter().fold(
            vec![0.0; partitioned_density[0].len()],
            |mut sum, d| {
                sum.iter_mut().zip(d).for_each(|(tpd, pd)| *tpd += pd);
                sum
            },
        );
        // the last value is is the vacuum and it has definitely been added
        let vacuum_density = partitioned_density.last().unwrap();
        let total_partitioned_density = total_density
            .iter()
            .zip(vacuum_density)
            .map(|(td, vd)| td - vd)
            .collect::<Vec<f64>>();
        // the volume is the same for all densities
        let total_volume: f64 = partitioned_volume.iter().sum();
        // the last value is is the vacuum and it has definitely been added
        let vacuum_volume = *partitioned_volume.last().unwrap();
        let total_partitioned_volume = total_volume - vacuum_volume;
        partitioned_density
            .iter()
            .zip(partitioned_volume)
            .zip(radius)
            .zip(errors)
            .for_each(|(((density, volume), radius), error)| {
                let mut row: Vec<String> =
                    Vec::with_capacity(rows[0].capacity());
                row.push(format!("{}", rows.len()));
                density.iter().for_each(|d| row.push(format!("{:.6}", d)));
                row.push(format!("{:.6}", volume));
                row.push(format!("{:.6}", radius));
                row.push(format!("{:.6}", error));
                for (i, col) in row.iter().enumerate() {
                    column_width[i] = column_width[i].max(col.len());
                }
                rows.push(row);
            });
        let mut footer = match density_num.cmp(&2) {
            std::cmp::Ordering::Less => format!(
                "\n * Vacuum Charge: {:>18.4}\n * Vacuum Volume: {:>18.4}\n * Partitioned Charge: {:>13.4}\n * Partitioned Volume: {:>13.4}",
                vacuum_density[0],
                vacuum_volume,
                total_partitioned_density[0],
                total_partitioned_volume,
            ),
            std::cmp::Ordering::Equal => format!(
                "\n * Vacuum Charge: {:>18.4}\n * Vacuum Spin: {:>20.4}\n * Vacuum Volume: {:>18.4}\n * Partitioned Charge: {:>13.4}\n * Partitioned Spin: {:>15.4}\n * Partitioned Volume: {:>13.4}",
                vacuum_density[0],
                vacuum_density[1],
                vacuum_volume,
                total_partitioned_density[0],
                total_partitioned_density[1],
                total_partitioned_volume,
            ),
            std::cmp::Ordering::Greater => format!(
                "\n * Vacuum Charge: {:>18.4}\n * Vacuum Spin X: {:>18.4}\n * Vacuum Spin Y: {:>18.4}\n * Vacuum Spin Z: {:>18.4}\n * Vacuum Volume: {:>18.4}\n * Partitioned Charge: {:>13.4}\n * Partitioned Spin X: {:>13.4}\n * Partitioned Spin Y: {:>13.4}\n * Partitioned Spin Z: {:>13.4}\n * Partitioned Volume: {:>13.4}",
                vacuum_density[0],
                vacuum_density[1],
                vacuum_density[2],
                vacuum_density[3],
                vacuum_volume,
                total_partitioned_density[0],
                total_partitioned_density[1],
                total_partitioned_density[2],
                total_partitioned_density[3],
                total_partitioned_volume,
            ),
        };
        footer.push_str(&format!(
            "\n * Boundary Voxels: {:>16.4}\n * Total Voxels: {:>19.4}",
            boundary_voxels, total_voxels
        ));
        Self {
            column_width,
            density_num,
            rows,
            footer,
        }
    }

    /// Creates and formats the header.
    fn format_header(&self) -> String {
        let mut header = String::new();
        let mut iter = self.column_width.iter();
        header.push_str(&format!(
            "| {:^width$} |",
            "Atom #",
            width = iter.next().unwrap()
        ));
        header.push_str(&format!(
            " {:^width$} |",
            "Charge",
            width = iter.next().unwrap()
        ));
        match self.density_num.cmp(&2) {
            std::cmp::Ordering::Equal => {
                header.push_str(&format!(
                    " {:^width$} |",
                    "Spin",
                    width = iter.next().unwrap()
                ));
            }
            std::cmp::Ordering::Greater => {
                header.push_str(&format!(
                    " {:^width$} |",
                    "Spin X",
                    width = iter.next().unwrap()
                ));
                header.push_str(&format!(
                    " {:^width$} |",
                    "Spin Y",
                    width = iter.next().unwrap()
                ));
                header.push_str(&format!(
                    " {:^width$} |",
                    "Spin Z",
                    width = iter.next().unwrap()
                ));
            }
            std::cmp::Ordering::Less => (),
        }
        header.push_str(&format!(
            " {:^width$} |",
            "Volume",
            width = iter.next().unwrap()
        ));
        header.push_str(&format!(
            " {:^width$} |",
            "Distance",
            width = iter.next().unwrap()
        ));
        header.push_str(&format!(
            " {:^width$} |\n",
            "Error",
            width = iter.next().unwrap()
        ));
        header
    }

    /// Creates and formats a separator.
    fn format_separator(&self) -> String {
        let mut separator = String::from("|");
        self.column_width.iter().for_each(|w| {
            separator.push_str(&format!("-{:-^width$}-|", "-", width = w));
        });
        separator
    }
}

impl fmt::Display for PartitionTable {
    /// Creates a String representation of the Table.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut table = String::from("## Partition Information\n");
        table.push_str(&self.format_header());
        self.rows.iter().for_each(|r| {
            if r.is_empty() {
                table.push_str(&self.format_separator());
            } else {
                let mut row = String::from("|");
                r.iter().zip(&self.column_width).for_each(|(s, w)| {
                    row.push_str(&format!(" {:>width$} |", s, width = w))
                });
                table.push_str(&row);
            }
            table.push('\n');
        });
        table.push_str(&self.footer);
        write!(f, "{}", table)
    }
}

/// Write a string to filename. Creates a new file regardless of what exists.
pub fn write(string: String, filename: String) -> std::io::Result<()> {
    let mut bader_file = File::create(filename)?;
    bader_file.write_all(string.as_bytes())?;
    Ok(())
}
