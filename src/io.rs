use crate::arguments::{Args, Reference};
use crate::atoms::Atoms;
use crate::density::Density;
use crate::utils;
use crate::voxel_map::VoxelMap;
use prettytable::{cell, format, row, Row, Table};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, Write};

/// File I/O for the gaussian cube format.
pub mod cube;
/// Custom BufReader.
pub mod reader;
/// File I/O for the VASP file format.
pub mod vasp;

/// Indicates the available file types of the density file.
pub enum FileType {
    /// CHGCAR, CHG and PARCHG.
    Vasp,
    /// Guassian, CP2K etc.
    Cube,
}

/// Return type of the read function in FileFormat.
pub type ReadFunction =
    io::Result<([f64; 3], [usize; 3], Atoms, Vec<Vec<f64>>)>;
/// Function type for the different add_row functions.
type AddRow =
    fn(&mut Table, String, (String, String, String), Vec<f64>, f64, f64);
/// Return type of the init function in FileFormat.
type InitReturn = (Vec<Vec<f64>>, Vec<f64>, Atoms, [usize; 3], [f64; 3]);

/// FileFormat trait. Used for handling input from a file.
pub trait FileFormat {
    /// Returns the parts required to build [`Density`] and [`Atoms`] structures.
    ///
    /// * `args`: [`Args`] parsed from the command line.
    fn init(&self, args: &Args) -> InitReturn {
        let (voxel_origin, grid, atoms, mut densities) =
            match self.read(args.file.clone()) {
                Ok(x) => x,
                Err(e) => panic!("Error: Problem reading file.\n{}", e),
            };
        if let Some(x) = args.spin.clone() {
            match densities.len() {
                1 => {
                    let (_, g, _, d) = match self.read(x.clone()) {
                        Ok(r) => r,
                        Err(e) => panic!("{}", e),
                    };
                    if 1 != d.len() {
                        panic!(
                               "Number of densities in original file is not 1.
Ambiguous how to handle spin density when {} contains {} densities.",
                               x,
                               d.len()
                        );
                    }
                    assert_eq!(g, grid,
                               "Error: Spin density has different grid size.");
                    densities.push(d[0].clone());
                }
                x => panic!(
                            "Number of densities in original file is not 1.
Ambiguous how to handle new spin when {} already has {} spin densities.",
                            args.file,
                            x - 1
                ),
            }
        }
        let rho = match args.reference.clone() {
            Reference::None => Vec::with_capacity(0),
            Reference::One(f) => {
                let (_, g, _, densities) = match self.read(f) {
                    Ok(r) => r,
                    Err(e) => panic!("{}", e),
                };
                assert_eq!(g, grid,
                           "Error: Reference density has different grid size.");
                densities[0].clone()
            }
            Reference::Two(f1, f2) => {
                let (_, g, _, densities) = match self.read(f1) {
                    Ok(r) => r,
                    Err(e) => panic!("{}", e),
                };
                assert_eq!(g, grid,
                           "Error: Reference density has different grid size.");
                let (_, g2, _, densities2) = match self.read(f2) {
                    Ok(r) => r,
                    Err(e) => panic!("{}", e),
                };

                assert_eq!(g2, grid,
                           "Error: Reference density has different grid size.");
                densities[0].par_iter()
                            .zip(&densities2[0])
                            .map(|(a, b)| a + b)
                            .collect::<Vec<f64>>()
            }
        };
        (densities, rho, atoms, grid, voxel_origin)
    }

    /// Returns the contents of the atoms charge file and bader charge file as
    /// Strings.
    ///
    /// * `voxel_map`: [`VoxelMap`] post [`VoxelMap::assign_atoms()`] and
    /// [`VoxelMap::charge_sum()`].
    /// * `atoms`: The associated [`Atoms`] struct for the density.
    /// * `density`: The reference &[`Density`].
    fn results(&self,
               voxel_map: VoxelMap,
               atoms: Atoms,
               density: &Density)
               -> (String, String) {
        let mut bader_table = Table::new();
        let mut atoms_table = Table::new();
        bader_table.set_format(table_format());
        atoms_table.set_format(table_format());
        let add_row: AddRow = match voxel_map.bader_charge.len() {
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
        let mut count = 1;
        let mut charge_t = vec![0f64; voxel_map.bader_charge.len()];
        for (i, s_dist) in voxel_map.surface_distance.iter().enumerate() {
            let mut charge_a = vec![0f64; 4];
            let mut volume_a = 0f64;
            for (ii, atom_num) in voxel_map.assigned_atom.iter().enumerate() {
                if *atom_num == i {
                    let index = format!("{}: {}", atom_num + 1, count);
                    let maxima_cartesian =
                        { density.to_cartesian(voxel_map.bader_maxima[ii]) };
                    let maxima_cartesian = utils::dot(maxima_cartesian,
                                                      density.voxel_lattice
                                                             .to_cartesian);
                    let (x, y, z) = self.coordinate_format(maxima_cartesian);
                    let volume = voxel_map.bader_volume[ii];
                    volume_a += volume;
                    let mut charge = vec![0f64; voxel_map.bader_charge.len()];
                    for j in 0..voxel_map.bader_charge.len() {
                        charge[j] = voxel_map.bader_charge[j][ii];
                        charge_a[j] += charge[j];
                        charge_t[j] += charge[j];
                    }
                    if charge[0] >= 1E-3 {
                        count += 1;
                        add_row(&mut bader_table,
                                index,
                                (x, y, z),
                                charge,
                                volume,
                                voxel_map.minimum_distance[ii]);
                    }
                }
            }
            count = 1;
            let (x, y, z) = self.coordinate_format(atoms.positions[i]);
            let index = format!("{}", i + 1);
            add_row(&mut atoms_table,
                    index,
                    (x, y, z),
                    charge_a,
                    volume_a,
                    *s_dist);
        }
        let footer = footer(voxel_map, charge_t);
        let mut atoms_charge_file = atoms_table.to_string();
        atoms_charge_file.push_str(&footer);
        let bader_charge_file = bader_table.to_string();
        (atoms_charge_file, bader_charge_file)
    }

    /// Reads the file into a [`ReadFunction`] containing the information
    /// needed from the file to build a [`Density`].
    ///
    /// * `filename`: The name of the file to read.
    fn read(&self, filename: String) -> ReadFunction;

    /// Reads the non-density section of the file into an [`Atoms`] object.
    ///
    /// * `atom_text`: The full string of non-density information from the
    /// density file.
    fn to_atoms(&self, atom_text: String) -> Atoms;

    /// Writes a specific density, data, to tile in the correct format.
    ///
    /// * `atoms`: The associated &[`Atoms`] object for the density file.
    /// * `data`: The density to write to file.
    fn write(&self, atoms: &Atoms, data: Vec<Vec<f64>>);

    /// How the format the positions of maxima and atoms
    ///
    /// * `coords`: The 3d representation of the position.
    fn coordinate_format(&self, coords: [f64; 3]) -> (String, String, String);
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
pub fn footer(voxel_map: VoxelMap, charge_total: Vec<f64>) -> String {
    let mut vacuum_charge = Vec::<f64>::with_capacity(4);
    for i in 0..voxel_map.bader_charge.len() {
        if voxel_map.bader_maxima[0] < 0 {
            vacuum_charge.push(*voxel_map.bader_charge[i].last().unwrap());
        } else {
            vacuum_charge.push(0.);
        }
    }
    let vacuum_volume = if voxel_map.bader_maxima[0] < 0 {
        *voxel_map.bader_volume.last().unwrap()
    } else {
        0.
    };
    match voxel_map.bader_charge.len() {
            1 => format!(
                "  Vacuum Charge: {:>18.4}\n  Vacuum Volume: {:>18.4}\n  Partitioned Charge: {:>13.4}",
                vacuum_charge[0],
                vacuum_volume,
                charge_total[0],
            ),
            2 =>  format!(
                "  Vacuum Charge: {:>18.4}\n  Vacuum Spin: {:>20.4}\n  Vacuum Volume: {:>18.4}\n  Partitioned Charge: {:>13.4}\n  Partitioned Spin: {:>15.4}",
                vacuum_charge[0],
                vacuum_charge[1],
                vacuum_volume,
                charge_total[0], charge_total[1]
            ),
            _ =>  format!(
                "  Vacuum Charge: {:>18.4}\n  Vacuum Spin X: {:>18.4}\n  Vacuum Spin Y: {:>18.4}\n  Vacuum Spin Z: {:>18.4}\n  Vacuum Volume: {:>18.4}\n  Partitioned Charge: {:>13.4}\n  Partitioned Spin X: {:>13.4}\n  Partitioned Spin Y: {:>13.4}\n  Partitioned Spin Z: {:>13.4}",
                vacuum_charge[0],
                vacuum_charge[1],
                vacuum_charge[2],
                vacuum_charge[3],
                vacuum_volume,
                charge_total[0], charge_total[1], charge_total[2], charge_total[3]
            ),
        }
}

/// Write the files
///
/// * `atoms_charge_file`: The contents, as a String, of the ACF.dat file.
/// * `bader_charge_file`: The contents, as a String, of the BCF.dat file.
pub fn write(atoms_charge_file: String,
             bader_charge_file: String)
             -> io::Result<()> {
    let mut bader_file = File::create("BCF.dat")?;
    bader_file.write_all(bader_charge_file.as_bytes())?;
    let mut atoms_file = File::create("ACF.dat")?;
    atoms_file.write_all(atoms_charge_file.as_bytes())?;
    Ok(())
}
