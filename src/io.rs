use crate::atoms::Atoms;
use crate::density::{Density, VoxelMap};
use crate::utils;
use prettytable::{cell, format, row, Table};
use std::fs::File;
use std::io::{self, Write};

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
pub fn results(voxel_map: VoxelMap,
               bader_charge: Vec<Vec<f64>>,
               bader_volume: Vec<usize>,
               assigned_atom: Vec<usize>,
               assigned_distance: Vec<f64>,
               surface_distance: Vec<f64>,
               vacuum_charge: f64,
               vacuum_volume: usize,
               density: Density,
               atoms: Atoms)
               -> io::Result<()> {
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

    let mut charge_t = vec![0f64; bader_charge.len()];
    for i in 0..atoms.positions.len() {
        let mut charge_a = vec![0f64; 4];
        let mut volume_a = 0f64;
        for (ii, atom_num) in assigned_atom.iter().enumerate() {
            if *atom_num == i {
                let bx = density.voxel_origin[0]
                         + (voxel_map.bader_maxima[ii]
                            / (density.size.y * density.size.z))
                           as f64;
                let by = density.voxel_origin[1]
                         + (voxel_map.bader_maxima[ii] / density.size.z).rem_euclid(density.size.y)
                           as f64;
                let bz = density.voxel_origin[2]
                         + voxel_map.bader_maxima[ii].rem_euclid(density.size.z)
                           as f64;
                let maxima_cartesian = utils::dot([bx, by, bz],
                                                  density.voxel_lattice
                                                         .to_cartesian);
                let x = format!("{:.6}", maxima_cartesian[2]);
                let y = format!("{:.6}", maxima_cartesian[1]);
                let z = format!("{:.6}", maxima_cartesian[0]);
                let volume =
                    bader_volume[ii] as f64 * density.voxel_lattice.volume;
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
        let x = format!("{:.6}", atoms.positions[i][2]);
        let y = format!("{:.6}", atoms.positions[i][1]);
        let z = format!("{:.6}", atoms.positions[i][0]);
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
            "  Vacuum Charge: {:>14.4}\n  Vacuum Volume: {:>14.4}\n  No. Electrons: {:>14.4}",
            vacuum_charge * density.voxel_lattice.volume,
            vacuum_volume as f64 * density.voxel_lattice.volume,
            charge_t[0],
        ),
        2 =>  format!(
            "  Vacuum Charge: {:>14.4}\n  Vacuum Volume: {:>14.4}\n  No. Electrons: {:>14.4}\n  Total Spin: {:>17.4}",
            vacuum_charge * density.voxel_lattice.volume,
            vacuum_volume as f64 * density.voxel_lattice.volume,
            charge_t[0], charge_t[1]
        ),
        _ =>  format!(
            "  Vacuum Charge: {:>14.4}\n  Vacuum Volume: {:>14.4}\n  No. Electrons: {:>14.4}\n  Total Spin X: {:>15.4}\n  Total Spin Y: {:>15.4}\n  Total Spin Z: {:>15.4}",
            vacuum_charge * density.voxel_lattice.volume,
            vacuum_volume as f64 * density.voxel_lattice.volume,
            charge_t[0], charge_t[1], charge_t[2], charge_t[3]
        ),
    };

    let mut atoms_out = atoms_table.to_string();
    atoms_out.push_str(&footer);
    let mut bader_file = File::create("BCF.dat")?;
    bader_file.write_all(bader_table.to_string().as_bytes())?;
    let mut atoms_file = File::create("ACF.dat")?;
    atoms_file.write_all(atoms_out.as_bytes())?;
    return Ok(());
}
