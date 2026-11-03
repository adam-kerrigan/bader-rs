use crate::atoms::{Atoms, Lattice};
use crate::io::{
    FileFormat, FortranFormat, ReadFunction, invalid_data, parse_line,
    read_line,
};
use crate::progress::progress_bar;
use crate::utils;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};

/// Convert from Bohr to Angstrom.
const LENGTH_UNITS: f64 = 0.52917721067;
/// Convert from Bohr^3 to Angstrom^3.
const VOLUME_UNITS: f64 = LENGTH_UNITS * LENGTH_UNITS * LENGTH_UNITS;

/// Structure for reading/writing a cube file.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cube {}

impl FileFormat for Cube {
    /// Reads a cube file from filename.
    ///
    /// The header is two comment lines, the number of atoms and origin, three
    /// lines of grid dimensions and axis vectors and a line per atom. The
    /// density follows with the last axis varying fastest.
    fn read(&self, filename: String) -> ReadFunction {
        let mut reader = BufReader::new(File::open(&filename)?);
        let mut buffer = String::new();
        let mut header = String::new();
        let mut next_line = |header: &mut String| -> std::io::Result<String> {
            let line = read_line(&mut reader, &mut buffer)?
                .ok_or_else(|| invalid_data(&filename, "cube"))?;
            header.push_str(line);
            Ok(line.to_string())
        };
        // two comment lines
        next_line(&mut header)?;
        next_line(&mut header)?;
        // number of atoms and the origin of the data
        let origin = parse_line::<f64>(&next_line(&mut header)?)
            .filter(|v| v.len() >= 4)
            .ok_or_else(|| invalid_data(&filename, "cube"))?;
        if origin[0] < 0.0 || origin.get(4).is_some_and(|&n| n != 1.0) {
            panic!("Error(Unsupported): Multiple values per voxel.");
        }
        let natoms = origin[0] as usize;
        let origin =
            [origin[1], origin[2], origin[3]].map(|x| x * LENGTH_UNITS);
        let mut grid_pts = [0usize; 3];
        for gp in &mut grid_pts {
            *gp = parse_line::<f64>(&next_line(&mut header)?)
                .and_then(|v| v.first().copied())
                .ok_or_else(|| invalid_data(&filename, "cube"))?
                as usize;
        }
        for _ in 0..natoms {
            next_line(&mut header)?;
        }
        let atoms = self.to_atoms(header);
        // the voxel origin is half a voxel from the data origin, in voxels
        let voxel_origin = utils::dot(origin, atoms.lattice.to_fractional);
        let voxel_origin =
            [0, 1, 2].map(|i| voxel_origin[i] * grid_pts[i] as f64 + 0.5);
        // the rest of the file is the density, convert out of Bohr
        let total = grid_pts.iter().product();
        let mut density = Vec::with_capacity(total);
        while let Some(line) = read_line(&mut reader, &mut buffer)? {
            let values = parse_line::<f64>(line)
                .ok_or_else(|| invalid_data(&filename, "cube"))?;
            density.extend(values.into_iter().map(|x| x / VOLUME_UNITS));
        }
        if density.len() != total {
            return Err(invalid_data(&filename, "cube"));
        }
        Ok((voxel_origin, grid_pts, atoms, vec![density]))
    }

    /// Read atoms information from file header.
    fn to_atoms(&self, atoms_text: String) -> Atoms {
        let lines = atoms_text.lines().skip(3);
        let numbers = |line: &str| {
            parse_line::<f64>(line).expect("cube header is not numeric")
        };
        // the axis lines are: number of voxels then the voxel vector in Bohr
        let mut lattice = [[0.0; 3]; 3];
        let mut axes = lines.clone().take(3).map(numbers);
        for vector in lattice.iter_mut() {
            let axis = axes.next().expect("cube header is missing axes");
            for (v, a) in vector.iter_mut().zip(&axis[1..4]) {
                *v = a * axis[0] * LENGTH_UNITS;
            }
        }
        let lattice = Lattice::new(lattice);
        // atom lines are: atomic number, charge, x, y, z in Bohr
        let positions = lines
            .skip(3)
            .map(|line| {
                let pos = numbers(line);
                let cartesian =
                    [pos[2], pos[3], pos[4]].map(|x| x * LENGTH_UNITS);
                // wrap the positions into the cell
                let fractional = utils::dot(cartesian, lattice.to_fractional)
                    .map(|x| x - x.floor());
                lattice.fractional_to_cartesian(fractional)
            })
            .collect();
        Atoms::new(lattice, positions, atoms_text)
    }

    /// Write a cube file.
    fn write(
        &self,
        atoms: &Atoms,
        data: Vec<f64>,
        filename: String,
        visible_pbar: bool,
    ) -> std::io::Result<()> {
        let filename = format!("{}.cube", filename);
        let mut buffer = BufWriter::new(File::create(&filename)?);
        let pbar = progress_bar(
            visible_pbar,
            data.len().div_ceil(6),
            &format!("Writing file {}:", filename),
        );
        buffer.write_all(atoms.text.as_bytes())?;
        for line in data.chunks(6) {
            for f in line {
                write!(buffer, " {:.5}", FortranFormat::new(f * VOLUME_UNITS))?;
            }
            writeln!(buffer)?;
            pbar.tick();
        }
        buffer.flush()
    }

    /// Cube files store the density with the first axis slowest so the
    /// coordinates do not need reordering.
    fn coordinate_format(&self, coords: [f64; 3]) -> [f64; 3] {
        coords
    }
}
