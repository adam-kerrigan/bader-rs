use crate::analysis::Analysis;
use crate::atoms::Atoms;
use crate::grid::Grid;
use crate::io::FileFormat;
use crate::utils;
use prettytable::{cell, format, row, Row, Table};
use std::fs::File;
use std::io::Write;

/// Function type for the different add_row functions.
type AddRow =
    fn(&mut Table, String, (String, String, String), Vec<f64>, f64, f64);

/// Writes the tables for outputting charge data.
pub fn charge_files(analysis: &Analysis,
                    atoms: &Atoms,
                    grid: &Grid,
                    file_type: Box<dyn FileFormat>)
                    -> (String, String) {
    let mut bader_table = Table::new();
    let mut atoms_table = Table::new();
    bader_table.set_format(table_format());
    atoms_table.set_format(table_format());
    let add_row: AddRow = match analysis.bader_charge.len() {
        1 => {
            bader_table.set_titles(bader_no_spin());
            atoms_table.set_titles(atom_no_spin());
            add_row_no_spin
        }
        2 => {
            bader_table.set_titles(bader_spin());
            atoms_table.set_titles(atom_spin());
            add_row_spin
        }
        4 => {
            bader_table.set_titles(bader_ncl_spin());
            atoms_table.set_titles(atom_ncl_spin());
            add_row_ncl
        }
        _ => panic!(),
    };
    let mut index: Vec<usize> = (0..analysis.bader_maxima.len()).collect();
    index.sort_by(|a, b| {
             analysis.assigned_atom[*a].cmp(&analysis.assigned_atom[*b])
         });
    let mut count = 1;
    let mut atom_num = 0;
    for i in index {
        atom_num = {
            let a = analysis.assigned_atom[i];
            if a != atom_num {
                count = 1;
            } else {
                count += 1;
            }
            a
        };
        let index = format!("{}: {}", atom_num + 1, count);
        let maxima_cartesian =
            { grid.to_cartesian(analysis.bader_maxima[i] as isize) };
        let maxima_cartesian =
            utils::dot(maxima_cartesian, grid.voxel_lattice.to_cartesian);
        let (x, y, z) = file_type.coordinate_format(maxima_cartesian);
        if analysis.bader_charge[0][i] >= 1E-3 {
            count += 1;
            add_row(&mut bader_table,
                    index,
                    (x, y, z),
                    analysis.bader_charge
                            .iter()
                            .map(|charge| charge[i])
                            .collect(),
                    analysis.bader_volume[i],
                    analysis.minimum_distance[i]);
        }
    }
    for (i, s_dist) in analysis.surface_distance.iter().enumerate() {
        let (x, y, z) = file_type.coordinate_format(atoms.positions[i]);
        let index = format!("{}", i + 1);
        add_row(&mut atoms_table,
                index,
                (x, y, z),
                analysis.atoms_charge
                        .iter()
                        .map(|charge| charge[i])
                        .collect(),
                analysis.atoms_volume[i],
                *s_dist);
    }
    let footer = footer(analysis);
    let mut atoms_charge_file = atoms_table.to_string();
    atoms_charge_file.push_str(&footer);
    let bader_charge_file = bader_table.to_string();
    (atoms_charge_file, bader_charge_file)
}

/// Creates a format for the output tables.
pub fn table_format() -> format::TableFormat {
    let line_position =
        &[format::LinePosition::Title, format::LinePosition::Bottom];
    let line_separator = format::LineSeparator::new('-', '+', '+', '+');
    format::FormatBuilder::new().column_separator('|')
                                .separators(line_position, line_separator)
                                .padding(1, 1)
                                .build()
}
/// Part of the collection of headers for the tables.
pub fn bader_no_spin() -> Row {
    row![c =>"#", "X", "Y", "Z", "Charge", "Volume", "Distance"]
}
/// Part of the collection of headers for the tables.
pub fn atom_no_spin() -> Row {
    row![c =>"#", "X", "Y", "Z", "Charge", "Volume", "Min. Dist."]
}
/// Part of the collection of headers for the tables.
pub fn bader_spin() -> Row {
    row![c =>"#", "X", "Y", "Z", "Charge", "Spin", "Volume", "Distance"]
}
/// Part of the collection of headers for the tables.
pub fn atom_spin() -> Row {
    row![c =>"#", "X", "Y", "Z", "Charge", "Spin", "Volume", "Min. Dist."]
}
/// Part of the collection of headers for the tables.
pub fn bader_ncl_spin() -> Row {
    row![c =>"#", "X", "Y", "Z", "Charge", "X Spin", "Y Spin", "Z Spin", "Volume", "Distance"]
}
/// Part of the collection of headers for the tables.
pub fn atom_ncl_spin() -> Row {
    row![c =>"#", "X", "Y", "Z", "Charge", "X Spin", "Y Spin", "Z Spin", "Volume", "Min. Dist."]
}

/// Part of the collection of functions for adding a row to the tables.
///
/// * `i_str`: The number to display in the '#' column of the table as a String.
/// * `pos_str`: The 3d position of the element as a String.
/// * `charge`: The densities of the current element.
/// * `volume`: The volume of the current element.
/// * `distance`: The minimum_distance or surface_distance of the current element.
pub fn add_row_no_spin(t: &mut Table,
                       i_str: String,
                       pos_str: (String, String, String),
                       charge: Vec<f64>,
                       volume: f64,
                       distance: f64) {
    let (x_str, y_str, z_str) = pos_str;
    let c_str = format!("{:.6}", charge[0]);
    let v_str = format!("{:.6}", volume);
    let d_str = format!("{:.6}", distance);
    t.add_row(row![r => i_str, x_str, y_str, z_str, c_str, v_str, d_str]);
}

/// Part of the collection of functions for adding a row to the tables.
///
/// * `i_str`: The number to display in the '#' column of the table as a String.
/// * `pos_str`: The 3d position of the element as a String.
/// * `charge`: The densities of the current element.
/// * `volume`: The volume of the current element.
/// * `distance`: The minimum_distance or surface_distance of the current element.
pub fn add_row_spin(t: &mut Table,
                    i_str: String,
                    pos_str: (String, String, String),
                    charge: Vec<f64>,
                    volume: f64,
                    distance: f64) {
    let (x_str, y_str, z_str) = pos_str;
    let c_str = format!("{:.6}", charge[0]);
    let s_str = format!("{:.6}", charge[1]);
    let v_str = format!("{:.6}", volume);
    let d_str = format!("{:.6}", distance);
    t.add_row(
        row![r => i_str, x_str, y_str, z_str, c_str, s_str, v_str, d_str],
    );
}

/// Part of the collection of functions for adding a row to the tables.
///
/// * `i_str`: The number to display in the '#' column of the table as a String.
/// * `pos_str`: The 3d position of the element as a String.
/// * `charge`: The densities of the current element.
/// * `volume`: The volume of the current element.
/// * `distance`: The minimum_distance or surface_distance of the current element.
pub fn add_row_ncl(t: &mut Table,
                   i_str: String,
                   pos_str: (String, String, String),
                   charge: Vec<f64>,
                   volume: f64,
                   distance: f64) {
    let (x_str, y_str, z_str) = pos_str;
    let c_str = format!("{:.6}", charge[0]);
    let sx_str = format!("{:.6}", charge[1]);
    let sy_str = format!("{:.6}", charge[2]);
    let sz_str = format!("{:.6}", charge[3]);
    let v_str = format!("{:.6}", volume);
    let d_str = format!("{:.6}", distance);
    t.add_row(row![r => i_str, x_str, y_str, z_str, c_str, sx_str, sy_str, sz_str, v_str, d_str]);
}

/// Produces the footer for the atoms charge file.
///
/// * `voxel_map`: [`VoxelMap`] post [`VoxelMap::charge_sum()`].
/// * `charge_total`: The sum of all the partitioned densities.
pub fn footer(analysis: &Analysis) -> String {
    match analysis.bader_charge.len() {
            1 => format!(
                "  Vacuum Charge: {:>18.4}\n  Vacuum Volume: {:>18.4}\n  Partitioned Charge: {:>13.4}",
                analysis.vacuum_charge[0],
                analysis.vacuum_volume,
                analysis.total_charge[0],
            ),
            2 =>  format!(
                "  Vacuum Charge: {:>18.4}\n  Vacuum Spin: {:>20.4}\n  Vacuum Volume: {:>18.4}\n  Partitioned Charge: {:>13.4}\n  Partitioned Spin: {:>15.4}",
                analysis.vacuum_charge[0],
                analysis.vacuum_charge[1],
                analysis.vacuum_volume,
                analysis.total_charge[0], analysis.total_charge[1]
            ),
            _ =>  format!(
                "  Vacuum Charge: {:>18.4}\n  Vacuum Spin X: {:>18.4}\n  Vacuum Spin Y: {:>18.4}\n  Vacuum Spin Z: {:>18.4}\n  Vacuum Volume: {:>18.4}\n  Partitioned Charge: {:>13.4}\n  Partitioned Spin X: {:>13.4}\n  Partitioned Spin Y: {:>13.4}\n  Partitioned Spin Z: {:>13.4}",
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
        }
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
