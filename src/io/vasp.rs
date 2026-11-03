use crate::atoms::{Atoms, Lattice};
use crate::io::{
    FileFormat, FortranFormat, ReadFunction, invalid_data, parse_line,
    read_line,
};
use crate::progress::progress_bar;
use crate::utils;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};

/// The VASP file format for reading/writing CHG, PARCHG and CHGCARs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vasp {}

impl FileFormat for Vasp {
    /// Read a VASP density.
    ///
    /// The file is a POSCAR followed by a blank line, the grid dimensions and
    /// the density with the first axis varying fastest. Spin densities follow
    /// each preceded by the grid dimensions again, in CHGCARs each density is
    /// also followed by augmentation occupancies which are skipped.
    fn read(&self, filename: String) -> ReadFunction {
        // the voxel origin in VASP is (0, 0, 0)
        let voxel_origin = [0f64; 3];
        let mut reader = BufReader::new(File::open(&filename)?);
        let mut buffer = String::new();
        // the POSCAR runs until the first blank line
        let mut header = String::new();
        loop {
            let line = read_line(&mut reader, &mut buffer)?
                .ok_or_else(|| invalid_data(&filename, "VASP"))?;
            header.push_str(line);
            if line.trim().is_empty() {
                break;
            }
        }
        // this line is the grid dimensions and separates each density
        let grid_line = read_line(&mut reader, &mut buffer)?
            .ok_or_else(|| invalid_data(&filename, "VASP"))?
            .to_string();
        header.push_str(&grid_line);
        let grid = parse_line::<usize>(&grid_line)
            .filter(|g| g.len() == 3)
            .ok_or_else(|| invalid_data(&filename, "VASP"))?;
        let total: usize = grid.iter().product();
        let atoms = self.to_atoms(header);
        // there could be a maximum of 4 densities, 1 total and then 1 or 3 spin
        let mut densities: Vec<Vec<f64>> = vec![Vec::with_capacity(total)];
        let mut skipping = false;
        while let Some(line) = read_line(&mut reader, &mut buffer)? {
            if line == grid_line {
                densities.push(Vec::with_capacity(total));
                skipping = false;
            } else if line.starts_with("aug") {
                // augmentation occupancies run until the next density
                skipping = true;
            } else if !skipping {
                let values = parse_line::<f64>(line)
                    .ok_or_else(|| invalid_data(&filename, "VASP"))?;
                // VASP stores the density multiplied by the volume
                densities.last_mut().unwrap().extend(
                    values.into_iter().map(|x| x / atoms.lattice.volume),
                );
            }
        }
        if densities.iter().any(|d| d.len() != total) {
            return Err(invalid_data(&filename, "VASP"));
        }
        // flip the grid points as VASP outputs density[z, y, x]
        let grid_pts = [grid[2], grid[1], grid[0]];
        Ok((voxel_origin, grid_pts, atoms, densities))
    }

    /// Read atom information from the POSCAR section of the file.
    ///
    /// As VASP stores the density with the first axis fastest the axes of
    /// the lattice and the positions are reversed so that the stored density
    /// is `[x, y, z]` in the reversed basis.
    fn to_atoms(&self, atoms_text: String) -> Atoms {
        let mut lines = atoms_text.lines().skip(1);
        let numbers = |line: &str| {
            parse_line::<f64>(line).expect("POSCAR header is not numeric")
        };
        let scale = numbers(lines.next().unwrap());
        let mut lattice = [[0.0; 3]; 3];
        for vector in lattice.iter_mut() {
            *vector = numbers(lines.next().unwrap())[..3].try_into().unwrap();
        }
        // a single negative scale is the volume of the cell, a single
        // positive scale is applied to every vector and three scales are
        // applied to the x, y and z components
        let scale: [f64; 3] = match scale[..] {
            [s] if s < 0.0 => {
                let volume = utils::vdot(
                    lattice[0],
                    utils::cross(lattice[1], lattice[2]),
                )
                .abs();
                [(-s / volume).cbrt(); 3]
            }
            [s] => [s; 3],
            [x, y, z] => [x, y, z],
            _ => panic!("POSCAR scaling factor must be 1 or 3 values."),
        };
        for vector in lattice.iter_mut() {
            for (v, s) in vector.iter_mut().zip(scale) {
                *v *= s;
            }
        }
        // reverse the axes and the components of each
        let reversed = |v: [f64; 3]| [v[2], v[1], v[0]];
        let lattice = Lattice::new([
            reversed(lattice[2]),
            reversed(lattice[1]),
            reversed(lattice[0]),
        ]);
        // the species line is optional, the counts line is all integers
        let mut line = lines.next().unwrap();
        if parse_line::<usize>(line).is_none() {
            line = lines.next().unwrap();
        }
        let total_atoms: usize =
            parse_line::<usize>(line).unwrap().iter().sum();
        // selective dynamics is optional before the coordinate type
        let mut coordinates = lines.next().unwrap().trim_start().to_lowercase();
        if coordinates.starts_with('s') {
            coordinates = lines.next().unwrap().trim_start().to_lowercase();
        }
        let fractional = coordinates.starts_with('d');
        let positions = lines
            .take(total_atoms)
            .map(|line| {
                let p = numbers(line);
                let p = reversed([p[0], p[1], p[2]]);
                let p = if fractional {
                    p
                } else {
                    // cartesian coordinates are scaled like the lattice
                    let p = [p[0] * scale[2], p[1] * scale[1], p[2] * scale[0]];
                    utils::dot(p, lattice.to_fractional)
                };
                lattice.fractional_to_cartesian(p.map(|x| x - x.floor()))
            })
            .collect();
        Atoms::new(lattice, positions, atoms_text)
    }

    /// Write a CHGCAR.
    fn write(
        &self,
        atoms: &Atoms,
        data: Vec<f64>,
        filename: String,
        visible_pbar: bool,
    ) -> std::io::Result<()> {
        let filename = format!("{}_CHGCAR", filename);
        let mut buffer = BufWriter::new(File::create(&filename)?);
        let pbar = progress_bar(
            visible_pbar,
            data.len().div_ceil(5),
            &format!("Writing file {}:", filename),
        );
        buffer.write_all(atoms.text.as_bytes())?;
        for line in data.chunks(5) {
            for f in line {
                write!(
                    buffer,
                    " {:.11}",
                    FortranFormat::new(f * atoms.lattice.volume)
                )?;
            }
            writeln!(buffer)?;
            pbar.tick();
        }
        buffer.flush()
    }

    /// Reverses the axes to undo the reordering of [`Vasp::to_atoms`].
    fn coordinate_format(&self, coords: [f64; 3]) -> [f64; 3] {
        [coords[2], coords[1], coords[0]]
    }
}
