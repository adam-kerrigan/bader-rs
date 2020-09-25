use crate::atoms::Atoms;
use crate::density::{Density, VoxelMap};
use crate::utils;
use prettytable::{cell, format, row, Table};
use std::fs::File;
use std::io::{self, Write};

pub mod cube;
pub mod reader;
pub mod vasp;

/// What is the coordinated system of the file
enum Coord {
    Fractional,
    Cartesian,
}

pub type ReadFunction =
    fn(String) -> io::Result<([f64; 3], [usize; 3], Atoms, Vec<Vec<f64>>)>;

/// Write the results to file
pub struct Results {
    atoms_charge_file: String,
    bader_charge_file: String,
}

impl Results {
    pub fn new(voxel_map: VoxelMap,
               bader_charge: Vec<Vec<f64>>,
               bader_volume: Vec<f64>,
               assigned_atom: Vec<usize>,
               assigned_distance: Vec<f64>,
               surface_distance: Vec<f64>,
               density: Density,
               atoms: Atoms,
               vasp_flag: bool)
               -> Self {
        let mut vacuum_charge = 0.;
        let mut vacuum_volume = 0.;
        let bader_maxima = if voxel_map.bader_maxima[0] == -1 {
            vacuum_charge = bader_charge[0][bader_charge[0].len() - 1];
            vacuum_volume = bader_volume[bader_volume.len() - 1];
            &voxel_map.bader_maxima[1..]
        } else {
            &voxel_map.bader_maxima[..]
        };
        let mut bader_table = Table::new();
        let mut atoms_table = Table::new();
        let format =
            format::FormatBuilder::new().column_separator('|')
                                        .separators(&[format::LinePosition::Title,
                                                      format::LinePosition::Bottom],
                                                    format::LineSeparator::new('-',
                                                                               '+',
                                                                               '+',
                                                                               '+'))
                                        .padding(1, 1)
                                        .build();
        bader_table.set_format(format);
        atoms_table.set_format(format);
        match bader_charge.len() {
            1 => {
                bader_table.set_titles(row![c =>"#", "X", "Y", "Z", "Charge", "Volume", "Distance"]);
                atoms_table.set_titles(row![c =>"#", "X", "Y", "Z", "Charge", "Volume", "Min. Dist."]);
            }
            2 => {
                bader_table
                    .set_titles(row![c =>"#", "X", "Y", "Z", "Charge", "Spin", "Volume", "Distance"]);
                atoms_table
                    .set_titles(row![c =>"#", "X", "Y", "Z", "Charge", "Spin", "Volume", "Min. Dist."]);
            }
            4 => {
                bader_table.set_titles(row![c =>"#", "X", "Y", "Z", "Charge", "X Spin", "Y Spin", "Z Spin", "Volume", "Distance"]);
                atoms_table.set_titles(row![c =>"#", "X", "Y", "Z", "Charge", "X Spin", "Y Spin", "Z Spin", "Volume", "Min. Dist."]);
            }
            _ => {}
        }

        fn add_row_no_spin(t: &mut Table,
                           x: String,
                           y: String,
                           z: String,
                           charge: Vec<f64>,
                           volume: f64,
                           distance: f64) {
            let i = format!("{}", t.len() + 1);
            let c = format!("{:.6}", charge[0]);
            let v = format!("{:.6}", volume);
            let d = format!("{:.6}", distance);
            t.add_row(row![r => i, x, y, z, c, v, d]);
        }

        fn add_row_spin(t: &mut Table,
                        x: String,
                        y: String,
                        z: String,
                        charge: Vec<f64>,
                        volume: f64,
                        distance: f64) {
            let i = format!("{}", t.len() + 1);
            let c = format!("{:.6}", charge[0]);
            let s = format!("{:.6}", charge[1]);
            let v = format!("{:.6}", volume);
            let d = format!("{:.6}", distance);
            t.add_row(row![r => i, x, y, z, c, s, v, d]);
        }

        fn add_row_ncl(t: &mut Table,
                       x: String,
                       y: String,
                       z: String,
                       charge: Vec<f64>,
                       volume: f64,
                       distance: f64) {
            let i = format!("{}", t.len() + 1);
            let c = format!("{:.6}", charge[0]);
            let sx = format!("{:.6}", charge[1]);
            let sy = format!("{:.6}", charge[2]);
            let sz = format!("{:.6}", charge[3]);
            let v = format!("{:.6}", volume);
            let d = format!("{:.6}", distance);
            t.add_row(row![r => i, x, y, z, c, sx, sy, sz, v, d]);
        }

        type AddRowType =
            fn(&mut Table, String, String, String, Vec<f64>, f64, f64);

        let add_row: AddRowType = match bader_charge.len() {
            1 => add_row_no_spin,
            2 => add_row_spin,
            _ => add_row_ncl,
        };

        fn vasp_format(coords: [f64; 3]) -> (String, String, String) {
            let x = format!("{:.6}", coords[2]);
            let y = format!("{:.6}", coords[1]);
            let z = format!("{:.6}", coords[0]);
            (x, y, z)
        }
        fn standard_format(coords: [f64; 3]) -> (String, String, String) {
            let x = format!("{:.6}", coords[0]);
            let y = format!("{:.6}", coords[1]);
            let z = format!("{:.6}", coords[2]);
            (x, y, z)
        }

        type PositionFormat = fn([f64; 3]) -> (String, String, String);

        let position_format: PositionFormat = if vasp_flag {
            vasp_format
        } else {
            standard_format
        };

        let mut charge_t = vec![0f64; bader_charge.len()];
        for i in 0..atoms.positions.len() {
            let mut charge_a = vec![0f64; 4];
            let mut volume_a = 0f64;
            for (ii, atom_num) in assigned_atom.iter().enumerate() {
                if *atom_num == i {
                    let bx = density.voxel_origin[0]
                             + (bader_maxima[ii]
                                / (density.size.y * density.size.z))
                               as f64;
                    let by = density.voxel_origin[1]
                             + (bader_maxima[ii] / density.size.z).rem_euclid(density.size.y)
                               as f64;
                    let bz = density.voxel_origin[2]
                             + bader_maxima[ii].rem_euclid(density.size.z)
                               as f64;
                    let maxima_cartesian = utils::dot([bx, by, bz],
                                                      density.voxel_lattice
                                                             .to_cartesian);
                    let (x, y, z) = position_format(maxima_cartesian);
                    let volume =
                        bader_volume[ii] * density.voxel_lattice.volume;
                    volume_a += volume;
                    let mut charge = vec![0f64; bader_charge.len()];
                    for j in 0..bader_charge.len() {
                        charge[j] =
                            bader_charge[j][ii] * density.voxel_lattice.volume;
                        charge_a[j] += charge[j];
                        charge_t[j] += charge[j];
                    }
                    if charge[0] >= 1E-3 {
                        add_row(&mut bader_table,
                                x,
                                y,
                                z,
                                charge,
                                volume,
                                assigned_distance[ii]);
                    }
                }
            }
            let (x, y, z) = position_format(atoms.positions[i]);
            add_row(&mut atoms_table,
                    x,
                    y,
                    z,
                    charge_a,
                    volume_a,
                    surface_distance[i]);
        }
        let footer = match bader_charge.len() {
            1 => format!(
                "  Vacuum Charge: {:>14.4}\n  Vacuum Volume: {:>14.4}\n  Total Charge: {:>15.4}",
                vacuum_charge * density.voxel_lattice.volume,
                vacuum_volume * density.voxel_lattice.volume,
                charge_t[0],
            ),
            2 =>  format!(
                "  Vacuum Charge: {:>14.4}\n  Vacuum Volume: {:>14.4}\n  Total Charge: {:>15.4}\n  Total Spin: {:>17.4}",
                vacuum_charge * density.voxel_lattice.volume,
                vacuum_volume * density.voxel_lattice.volume,
                charge_t[0], charge_t[1]
            ),
            _ =>  format!(
                "  Vacuum Charge: {:>14.4}\n  Vacuum Volume: {:>14.4}\n  Total Charge: {:>15.4}\n  Total Spin X: {:>15.4}\n  Total Spin Y: {:>15.4}\n  Total Spin Z: {:>15.4}",
                vacuum_charge * density.voxel_lattice.volume,
                vacuum_volume * density.voxel_lattice.volume,
                charge_t[0], charge_t[1], charge_t[2], charge_t[3]
            ),
        };

        let mut atoms_charge_file = atoms_table.to_string();
        atoms_charge_file.push_str(&footer);
        let bader_charge_file = bader_table.to_string();
        Self { atoms_charge_file,
               bader_charge_file }
    }

    /// Write the files
    pub fn write(self) -> io::Result<()> {
        let mut bader_file = File::create("BCF.dat")?;
        bader_file.write_all(self.bader_charge_file.as_bytes())?;
        let mut atoms_file = File::create("ACF.dat")?;
        atoms_file.write_all(self.atoms_charge_file.as_bytes())?;
        return Ok(());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arguments::Weight;
    use crate::atoms::Lattice;

    #[test]
    fn io_results_new_no_spin() {
        let mut map: Vec<isize> = vec![-1; 4usize.pow(3)];
        let mut densities: Vec<Vec<f64>> = vec![vec![0.; 4usize.pow(3)]];
        for i in [0usize, 1, 3, 4, 5, 7, 12, 13, 15, 16, 17, 19, 20, 21, 23,
                  28, 29, 31, 48, 49, 51, 52, 53, 55, 60, 61, 63].iter()
        {
            map[*i] = 0;
            densities[0][*i] = 1.;
        }
        let bader_maxima = vec![-1, 0];
        let voxel_map = VoxelMap::new(map, bader_maxima);
        let lattice = Lattice::new([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]);
        let atoms = Atoms::new(lattice, vec![[0., 0., 0.]], String::new());
        let lattice = Lattice::new([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]);
        let reference = Density::new(&densities[0],
                                     [4, 4, 4],
                                     lattice.to_cartesian,
                                     Some(1E-3),
                                     [0., 0., 0.1]);
        let (assigned_atom, assigned_distance) =
            atoms.assign_maxima(&voxel_map, &reference);
        let (bader_charge, bader_volume, surface_distance) = {
            voxel_map.charge_sum(&densities,
                                 &assigned_atom,
                                 &atoms,
                                 &reference,
                                 (0..64).collect(),
                                 Weight::None)
        };
        let min_len = reference.voxel_lattice.distance_matrix[4];
        let results = Results::new(voxel_map,
                                   bader_charge,
                                   bader_volume,
                                   assigned_atom,
                                   assigned_distance,
                                   surface_distance,
                                   reference,
                                   atoms,
                                   true);
        let acf = match results.atoms_charge_file.split("\n").nth(2) {
            Some(s) => s.split(" | ")
                        .map(|x| x.trim().parse::<f64>().unwrap())
                        .collect::<Vec<f64>>(),
            //Some(s) => s.split(" | ").collect::<Vec<_>>(),
            None => panic!("No acf table content"),
        };
        let bcf = match results.bader_charge_file.split("\n").nth(2) {
            Some(s) => s.split(" | ")
                        .map(|x| x.trim().parse::<f64>().unwrap())
                        .collect::<Vec<f64>>(),
            //Some(s) => s.split(" | ").collect::<Vec<_>>(),
            None => panic!("No bcf table content"),
        };
        assert_eq!(bcf, [1., 0.1, 0., 0., 27., 27., 0.1]);
        assert_eq!(acf, [1., 0., 0., 0., 27., 27., min_len - 0.1]);
    }

    #[test]
    fn io_results_new_spin() {
        let mut map: Vec<isize> = vec![-1; 4usize.pow(3)];
        let mut densities: Vec<Vec<f64>> = vec![vec![0.; 4usize.pow(3)]; 2];
        for i in [0usize, 1, 3, 4, 5, 7, 12, 13, 15, 16, 17, 19, 20, 21, 23,
                  28, 29, 31, 48, 49, 51, 52, 53, 55, 60, 61, 63].iter()
        {
            map[*i] = 0;
            densities[0][*i] = 1.;
            densities[1][*i] = 2.;
        }
        let bader_maxima = vec![-1, 0];
        let voxel_map = VoxelMap::new(map, bader_maxima);
        let lattice = Lattice::new([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]);
        let atoms = Atoms::new(lattice, vec![[0., 0., 0.]], String::new());
        let lattice = Lattice::new([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]);
        let reference = Density::new(&densities[0],
                                     [4, 4, 4],
                                     lattice.to_cartesian,
                                     Some(1E-3),
                                     [0., 0., 0.1]);
        let (assigned_atom, assigned_distance) =
            atoms.assign_maxima(&voxel_map, &reference);
        let (bader_charge, bader_volume, surface_distance) = {
            voxel_map.charge_sum(&densities,
                                 &assigned_atom,
                                 &atoms,
                                 &reference,
                                 (0..64).collect(),
                                 Weight::None)
        };
        let min_len = reference.voxel_lattice.distance_matrix[4];
        let results = Results::new(voxel_map,
                                   bader_charge,
                                   bader_volume,
                                   assigned_atom,
                                   assigned_distance,
                                   surface_distance,
                                   reference,
                                   atoms,
                                   false);
        let acf = match results.atoms_charge_file.split("\n").nth(2) {
            Some(s) => s.split(" | ")
                        .map(|x| x.trim().parse::<f64>().unwrap())
                        .collect::<Vec<f64>>(),
            None => panic!("No acf table content"),
        };
        let bcf = match results.bader_charge_file.split("\n").nth(2) {
            Some(s) => s.split(" | ")
                        .map(|x| x.trim().parse::<f64>().unwrap())
                        .collect::<Vec<f64>>(),
            None => panic!("No bcf table content"),
        };
        assert_eq!(bcf, [1., 0., 0., 0.1, 27., 27. * 2., 27., 0.1]);
        assert_eq!(acf, [1., 0., 0., 0., 27., 27. * 2., 27., min_len - 0.1]);
    }

    #[test]
    fn io_results_new_ncl_spin() {
        let mut map: Vec<isize> = vec![-1; 4usize.pow(3)];
        let mut densities: Vec<Vec<f64>> = vec![vec![0.; 4usize.pow(3)]; 4];
        for i in [0usize, 1, 3, 4, 5, 7, 12, 13, 15, 16, 17, 19, 20, 21, 23,
                  28, 29, 31, 48, 49, 51, 52, 53, 55, 60, 61, 63].iter()
        {
            map[*i] = 0;
            densities[0][*i] = 1.;
            densities[1][*i] = 2.;
            densities[2][*i] = 3.;
            densities[3][*i] = 4.;
        }
        let bader_maxima = vec![-1, 0];
        let voxel_map = VoxelMap::new(map, bader_maxima);
        let lattice = Lattice::new([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]);
        let atoms = Atoms::new(lattice, vec![[0., 0., 0.]], String::new());
        let lattice = Lattice::new([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]);
        let reference = Density::new(&densities[0],
                                     [4, 4, 4],
                                     lattice.to_cartesian,
                                     Some(1E-3),
                                     [0., 0., 0.]);
        let (assigned_atom, assigned_distance) =
            atoms.assign_maxima(&voxel_map, &reference);
        let (bader_charge, bader_volume, surface_distance) = {
            voxel_map.charge_sum(&densities,
                                 &assigned_atom,
                                 &atoms,
                                 &reference,
                                 (0..64).collect(),
                                 Weight::None)
        };
        let min_len = reference.voxel_lattice.distance_matrix[4];
        let results = Results::new(voxel_map,
                                   bader_charge,
                                   bader_volume,
                                   assigned_atom,
                                   assigned_distance,
                                   surface_distance,
                                   reference,
                                   atoms,
                                   true);
        let acf = match results.atoms_charge_file.split("\n").nth(2) {
            Some(s) => s.split(" | ")
                        .map(|x| x.trim().parse::<f64>().unwrap())
                        .collect::<Vec<f64>>(),
            None => panic!("No acf table content"),
        };
        let bcf = match results.bader_charge_file.split("\n").nth(2) {
            Some(s) => s.split(" | ")
                        .map(|x| x.trim().parse::<f64>().unwrap())
                        .collect::<Vec<f64>>(),
            None => panic!("No bcf table content"),
        };
        assert_eq!(bcf,
                   [1.,
                    0.,
                    0.,
                    0.,
                    27.,
                    27. * 2.,
                    27. * 3.,
                    27. * 4.,
                    27.,
                    0.]);
        assert_eq!(acf,
                   [1.,
                    0.,
                    0.,
                    0.,
                    27.,
                    27. * 2.,
                    27. * 3.,
                    27. * 4.,
                    27.,
                    min_len]);
    }
}
