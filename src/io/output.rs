use crate::analysis::Analysis;
use crate::atoms::Atoms;
use crate::grid::Grid;
use crate::io::{FileFormat, WriteType};
use crate::progress::Bar;
use crate::utils;
use crate::voxel_map::VoxelMap;
use std::fs::File;
use std::io::Write;

/// Enum of available tables.
enum TableType {
    /// Table for the ACF file.
    AtomsCharge,
    /// Table for the BCF file.
    BaderCharge,
}

/// Structure that contains and builds the table.
struct Table {
    /// How wide each column is.
    column_width: Vec<usize>,
    /// The number of charge and spin densities.
    density_num: usize,
    /// The rows of the table as a vector of strings.
    rows: Vec<Vec<String>>,
    /// What type of table the structure is.
    table_type: TableType,
}

impl Table {
    /// Creates a new structure and sets the minimum widths of each.
    fn new(table_type: TableType, density_num: usize) -> Self {
        let rows = vec![Vec::with_capacity(0)];
        let mut column_width = Vec::with_capacity(6 + density_num);
        column_width.push(1);
        column_width.push(1);
        column_width.push(1);
        column_width.push(1);
        column_width.push(6);
        match density_num.cmp(&2) {
            std::cmp::Ordering::Greater => {
                column_width.push(6);
                column_width.push(6);
                column_width.push(6);
            }
            std::cmp::Ordering::Equal => column_width.push(6),
            std::cmp::Ordering::Less => (),
        };
        column_width.push(6);
        column_width.push(8);
        Self { column_width,
               density_num,
               rows,
               table_type }
    }

    /// Adds a row the table.
    #[allow(clippy::borrowed_box)]
    fn add_row(&mut self,
               index: usize,
               p: [f64; 3],
               density: &[f64],
               volume: f64,
               distance: f64,
               file_type: &Box<dyn FileFormat>) {
        let mut row: Vec<String> = Vec::with_capacity(6 + self.density_num);
        row.push(format!("{}", index));
        let coord = file_type.coordinate_format(p);
        row.push(coord.0);
        row.push(coord.1);
        row.push(coord.2);
        density.iter().for_each(|d| row.push(format!("{:.6}", d)));
        row.push(format!("{:.6}", volume));
        row.push(format!("{:.6}", distance));
        for (i, col) in row.iter().enumerate() {
            self.column_width[i] = self.column_width[i].max(col.len());
        }
        self.rows.push(row);
    }

    /// Adds a blank row to be a separator in the final string.
    fn add_separator(&mut self) {
        self.rows.push(Vec::with_capacity(0));
    }

    /// Creates and formats the footer.
    fn format_footer(&self, analysis: &Analysis) -> String {
        match self.table_type {
            TableType::AtomsCharge => {
                let mut separator = self.format_separator(0);
                let footer = match self.density_num.cmp(&2) {
                        std::cmp::Ordering::Less => format!(
                            "\n  Vacuum Charge: {:>18.4}\n  Vacuum Volume: {:>18.4}\n  Partitioned Charge: {:>13.4}",
                            analysis.vacuum_charge[0],
                            analysis.vacuum_volume,
                            analysis.total_charge[0],
                        ),
                        std::cmp::Ordering::Equal =>  format!(
                            "\n  Vacuum Charge: {:>18.4}\n  Vacuum Spin: {:>20.4}\n  Vacuum Volume: {:>18.4}\n  Partitioned Charge: {:>13.4}\n  Partitioned Spin: {:>15.4}",
                            analysis.vacuum_charge[0],
                            analysis.vacuum_charge[1],
                            analysis.vacuum_volume,
                            analysis.total_charge[0], analysis.total_charge[1]
                        ),
                        std::cmp::Ordering::Greater =>  format!(
                            "\n  Vacuum Charge: {:>18.4}\n  Vacuum Spin X: {:>18.4}\n  Vacuum Spin Y: {:>18.4}\n  Vacuum Spin Z: {:>18.4}\n  Vacuum Volume: {:>18.4}\n  Partitioned Charge: {:>13.4}\n  Partitioned Spin X: {:>13.4}\n  Partitioned Spin Y: {:>13.4}\n  Partitioned Spin Z: {:>13.4}",
                            analysis.vacuum_charge[0],
                            analysis.vacuum_charge[1],
                            analysis.vacuum_charge[2],
                            analysis.vacuum_charge[3],
                            analysis.vacuum_volume,
                            analysis.total_charge[0],
                            analysis.total_charge[1],
                            analysis.total_charge[2],
                            analysis.total_charge[3]
                        ),
                };
                separator.push_str(&footer);
                separator
            }
            TableType::BaderCharge => String::new(),
        }
    }

    /// Creates and formats the header.
    fn format_header(&self) -> String {
        let mut header = String::new();
        let mut iter = self.column_width.iter();
        header.push_str(&format!(" {:^width$} |",
                                 "#",
                                 width = iter.next().unwrap()));
        header.push_str(&format!(" {:^width$} |",
                                 "X",
                                 width = iter.next().unwrap()));
        header.push_str(&format!(" {:^width$} |",
                                 "Y",
                                 width = iter.next().unwrap()));
        header.push_str(&format!(" {:^width$} |",
                                 "Z",
                                 width = iter.next().unwrap()));
        header.push_str(&format!(" {:^width$} |",
                                 "Charge",
                                 width = iter.next().unwrap()));
        match self.table_type {
            TableType::AtomsCharge | TableType::BaderCharge => {
                match self.density_num.cmp(&2) {
                    std::cmp::Ordering::Equal => {
                        header.push_str(&format!(" {:^width$} |",
                                                 "Spin",
                                                 width = iter.next().unwrap()));
                    }
                    std::cmp::Ordering::Greater => {
                        header.push_str(&format!(" {:^width$} |",
                                                 "Spin X",
                                                 width = iter.next().unwrap()));
                        header.push_str(&format!(" {:^width$} |",
                                                 "Spin Y",
                                                 width = iter.next().unwrap()));
                        header.push_str(&format!(" {:^width$} |",
                                                 "Spin Z",
                                                 width = iter.next().unwrap()));
                    }
                    std::cmp::Ordering::Less => (),
                }
            }
        }
        header.push_str(&format!(" {:^width$} |",
                                 "Volume",
                                 width = iter.next().unwrap()));
        header.push_str(&format!(" {:^width$}\n",
                                 "Distance",
                                 width = iter.next().unwrap()));
        header
    }

    /// Creates and formats a separator.
    fn format_separator(&self, i: usize) -> String {
        let mut separator = String::new();
        self.column_width.iter().for_each(|w| {
            separator.push_str(&format!("-{:-^width$}-+", "-", width = w));
        });
        separator.pop();
        separator.pop();
        if let TableType::BaderCharge = self.table_type {
            let len = self.column_width[0];
            separator.replace_range(1..(len + 7),
                                    &format!("Atom: {:>width$}",
                                             i,
                                             width = len));
        }
        separator
    }

    /// Creates a String representation of the Table.
    fn to_string(&self, analysis: &Analysis) -> String {
        let mut i = 0;
        let mut table = String::new();
        table.push_str(&self.format_header());
        self.rows.iter().for_each(|r| {
                            if r.is_empty() {
                                i += 1;
                                table.push_str(&self.format_separator(i));
                            } else {
                                let mut row = String::new();
                                r.iter()
                                 .zip(&self.column_width)
                                 .for_each(|(s, w)| {
                                     row.push_str(&format!(" {:>width$} |",
                                                           s,
                                                           width = w))
                                 });
                                row.pop();
                                table.push_str(&row);
                            }
                            table.push('\n');
                        });
        table.push_str(&self.format_footer(analysis));
        table
    }
}

/// Writes the tables for outputting charge data.
///
/// * analysis: The [`Analysis`] structure to be tabulated.
/// * atoms: The [`Atoms`] referenced to the structure.
/// * grid: The [`Grid`] structure for moving around the analysed density.
/// * file_type: [`FileFormat`] for printing the correct coordinates.
///
/// ### Returns:
/// (String, String): The ACF and BCF as Strings.
#[allow(clippy::borrowed_box)]
pub fn charge_files(analysis: &Analysis,
                    atoms: &Atoms,
                    grid: &Grid,
                    file_type: &Box<dyn FileFormat>)
                    -> (String, String) {
    let mut bader_table =
        Table::new(TableType::BaderCharge, analysis.bader_charge.len());
    let mut atoms_table =
        Table::new(TableType::AtomsCharge, analysis.bader_charge.len());
    let mut index: Vec<usize> = (0..analysis.bader_maxima.len()).collect();
    index.sort_by(|a, b| {
             analysis.assigned_atom[*a].cmp(&analysis.assigned_atom[*b])
         });
    let mut atom_num = 0;
    atoms_table.add_row(atom_num + 1,
                        atoms.positions[atom_num],
                        &analysis.atoms_charge
                                 .iter()
                                 .map(|charge| charge[atom_num])
                                 .collect::<Vec<f64>>(),
                        analysis.atoms_volume[atom_num],
                        analysis.surface_distance[atom_num],
                        file_type);
    for i in index {
        atom_num = {
            let a = analysis.assigned_atom[i];
            if a != atom_num {
                bader_table.add_separator();
                atoms_table.add_row(a + 1,
                                    atoms.positions[a],
                                    &analysis.atoms_charge
                                             .iter()
                                             .map(|charge| charge[a])
                                             .collect::<Vec<f64>>(),
                                    analysis.atoms_volume[a],
                                    analysis.surface_distance[a],
                                    file_type);
            }
            a
        };
        let maxima_cartesian =
            { grid.to_cartesian(analysis.bader_maxima[i] as isize) };
        let maxima_cartesian =
            utils::dot(maxima_cartesian, grid.voxel_lattice.to_cartesian);
        if analysis.bader_charge[0][i] >= grid.maxima_tolerance {
            bader_table.add_row(i + 1,
                                maxima_cartesian,
                                &analysis.bader_charge
                                         .iter()
                                         .map(|charge| charge[i])
                                         .collect::<Vec<f64>>(),
                                analysis.bader_volume[i],
                                analysis.minimum_distance[i],
                                file_type);
        }
    }
    let bader_charge_file = bader_table.to_string(analysis);
    let atoms_charge_file = atoms_table.to_string(analysis);
    (atoms_charge_file, bader_charge_file)
}

/// Write the files
///
/// * `atoms_charge_file`: The contents, as a String, of the ACF.dat file.
/// * `bader_charge_file`: The contents, as a String, of the BCF.dat file.
pub fn write(atoms_charge_file: String,
             bader_charge_file: String)
             -> std::io::Result<()> {
    let mut bader_file = File::create("BCF.dat")?;
    bader_file.write_all(bader_charge_file.as_bytes())?;
    let mut atoms_file = File::create("ACF.dat")?;
    atoms_file.write_all(atoms_charge_file.as_bytes())?;
    Ok(())
}

/// Write the densities of either Bader atoms or volumes.
#[allow(clippy::borrowed_box)]
pub fn write_densities(atoms: &Atoms,
                       analysis: &Analysis,
                       densities: Vec<Vec<f64>>,
                       grid: &Grid,
                       output: WriteType,
                       voxel_map: &VoxelMap,
                       file_type: &Box<dyn FileFormat>)
                       -> std::io::Result<()> {
    let filename = match densities.len().cmp(&2) {
        std::cmp::Ordering::Less => vec![String::from("charge")],
        std::cmp::Ordering::Equal => {
            vec![String::from("charge"), String::from("spin")]
        }
        std::cmp::Ordering::Greater => vec![String::from("charge"),
                                            String::from("spin_x"),
                                            String::from("spin_y"),
                                            String::from("spin_z"),],
    };
    match output {
        WriteType::Atom(a) => {
            println!("Writing out charge densities for atoms:");
            let atom_iter = if a.is_empty() {
                (0..atoms.positions.len()).collect()
            } else {
                a
            };
            for atom in atom_iter {
                println!("Atom {}:", atom + 1);
                let pbar = Bar::visible(grid.size.total as u64,
                                        100,
                                        String::from("Building map:"));
                let map = analysis.output_atom_map(grid, voxel_map, atom, pbar);
                for (i, den) in densities.iter().enumerate() {
                    let fname = format!("atom_{}_{}", atom + 1, filename[i]);
                    let pbar =
                        Bar::visible(1, 100, format!("Writing {}:", fname));
                    let den =
                        den.iter()
                           .zip(&map)
                           .map(|(d, weight)| weight.as_ref().map(|w| d * w))
                           .collect::<Vec<Option<f64>>>();
                    file_type.write(atoms, den, fname, pbar)?;
                }
            }
        }
        WriteType::Volume(v) => {
            println!("Writing out charge densities for atoms:");
            let volume_iter = if v.is_empty() {
                (0..analysis.bader_maxima.len()).collect()
            } else {
                v
            };
            for volume in volume_iter {
                print!("Atom: {}.", volume + 1);
                let pbar = Bar::visible(grid.size.total as u64,
                                        100,
                                        String::from("Building map:"));
                let map =
                    analysis.output_volume_map(grid, voxel_map, volume, pbar);
                for (i, den) in densities.iter().enumerate() {
                    let fname =
                        format!("volume_{}_{}", volume + 1, filename[i]);
                    let pbar =
                        Bar::visible(1, 100, format!("Writing {}:", fname));
                    let den =
                        den.iter()
                           .zip(&map)
                           .map(|(d, weight)| weight.as_ref().map(|w| d * w))
                           .collect::<Vec<Option<f64>>>();
                    file_type.write(atoms, den, fname, pbar)?;
                }
                println!(" Done.");
            }
        }
        WriteType::None => (),
    }
    Ok(())
}
