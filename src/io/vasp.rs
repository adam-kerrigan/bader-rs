use crate::atoms::{Atoms, Lattice};
use crate::io::reader::BufReader;
use crate::io::{FileFormat, FortranFormat, ReadFunction};
use crate::progress::{Bar, HiddenBar, ProgressBar};
use crate::utils;
use std::fs::File;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};

/// The coordinate system.
enum Coord {
    /// Fractional coordinates.
    Fractional,
    /// Cartesian coordinates.
    Cartesian,
}

/// The VASP file format for reading/writing CHG, PARCHG and CHGCARs.
pub struct Vasp {}

impl FileFormat for Vasp {
    /// Read a VASP density.
    fn read(&self, filename: String) -> ReadFunction {
        // the voxel origin in VASP is (0, 0, 0)
        let voxel_origin = [0f64; 3];
        // find the start and end points of the density as well as the total file size
        let (grid, aug, total) = {
            // open the file in a buffer reader
            let mut reader = BufReader::open(filename.clone())?;

            let mut buffer = String::new();
            let mut grid: Vec<[usize; 2]> = vec![];
            let mut aug: Vec<usize> = vec![];
            let mut pos = 0;
            // the first 7 lines are useless to us
            for _ in 0..8 {
                let size = match reader.read_line(&mut buffer) {
                    Some(line) => {
                        let (_, size) = line?;
                        size
                    }
                    None => 0,
                };
                pos += size;
            }
            // lets find the start
            while let Some(line) = reader.read_line(&mut buffer) {
                let (text, size) = line?;
                pos += size;
                // empty line before the grid spacing
                if text.trim().is_empty() {
                    break;
                }
            }
            // this line is the grid spacing
            let mut buffer = String::new();
            let grid_spacing = reader
                .read_line(&mut buffer)
                .map(|line| {
                    let (text, size) = line.unwrap();
                    let start = pos;
                    pos += size;
                    let end = pos;
                    grid.push([start, end]);
                    text
                })
                .unwrap();
            // lets fast forward a bit
            let grid_length = grid_spacing
                .split_whitespace()
                .fold(1, |acc, val| val.parse::<usize>().unwrap() * acc);
            let mut buffer = String::new();
            let per_row = reader
                .read_line(&mut buffer)
                .map(|line| {
                    let (text, size) = line.unwrap();
                    pos += size;
                    text
                })
                .unwrap()
                .split_whitespace()
                .count();
            let grid_length = if grid_length.rem_euclid(per_row) == 0 {
                grid_length / per_row
            } else {
                grid_length / per_row + 1
            };
            let mut buffer = String::new();
            for _ in 1..grid_length {
                pos += reader.read_line(&mut buffer).unwrap().unwrap().1;
            }
            // lets start trying to match
            let mut buffer = String::new();
            while let Some(line) = reader.read_line(&mut buffer) {
                let (text, size) = line?;
                if text == grid_spacing {
                    let start = pos;
                    pos += size;
                    let end = pos;
                    grid.push([start, end]);
                } else if text.starts_with("aug") && aug.len() < grid.len() {
                    aug.push(pos);
                    pos += size;
                } else {
                    pos += size;
                }
            }
            (grid, aug, pos)
        };
        // Now we know where everything is so let's work out what to do
        // Start by making vector of start and end points of the densities
        let mut start = Vec::with_capacity(4);
        let mut stop = Vec::with_capacity(4);
        for (i, start_stop) in grid.iter().enumerate() {
            start.push(start_stop[1]);
            let s = if !aug.is_empty() {
                aug[i * aug.len() / grid.len()]
            } else if grid.len() > (i + 1) {
                grid[i + 1][0]
            } else {
                total
            };
            stop.push(s);
        }
        let mut file = File::open(filename)?;
        // assign Vectos with the capacity of what it is to hold
        let mut poscar_b = Vec::with_capacity(grid[0][0]);
        let mut grid_pts_b = Vec::with_capacity(grid[0][1] - grid[0][0]);
        let mut density_b = Vec::with_capacity(stop[0] - start[0]);
        // there could be a maximum of 4 densities 1 total and then 1 or 3 spin
        let mut density: Vec<Vec<f64>> = Vec::with_capacity(4);
        // read the poscar information poscar_b
        let _ = <File as Read>::by_ref(&mut file)
            .take(grid[0][1] as u64)
            .read_to_end(&mut poscar_b)?;
        file.seek(SeekFrom::Current(grid[0][0] as i64 - grid[0][1] as i64))?;
        // read the grid line into grid_pts_b
        let _ = <File as Read>::by_ref(&mut file)
            .take((grid[0][1] - grid[0][0]) as u64)
            .read_to_end(&mut grid_pts_b)?;
        // read the total charge density into density_b
        let _ = <File as Read>::by_ref(&mut file)
            .take((stop[0] - start[0]) as u64)
            .read_to_end(&mut density_b)?;
        // convert the bytes we have read into a String and an Atoms struct
        let poscar = String::from_utf8(poscar_b).unwrap();
        let grid_vec: Vec<usize> = {
            String::from_utf8(grid_pts_b)
                .unwrap()
                .split_whitespace()
                .map(|x| x.parse::<usize>().unwrap())
                .collect()
        };
        let atoms = self.to_atoms(poscar);
        // convert out of VASP's strange units
        density.push(
            String::from_utf8(density_b)
                .unwrap()
                .split_whitespace()
                .map(|x| x.parse::<f64>().unwrap() / atoms.lattice.volume)
                .collect(),
        );
        for i in 1..start.len() {
            let mut spin_b = Vec::with_capacity(stop[i] - start[i]);
            let _ = file.seek(SeekFrom::Start(start[i] as u64));
            let _ = <File as Read>::by_ref(&mut file)
                .take((stop[i] - start[i]) as u64)
                .read_to_end(&mut spin_b)?;
            density.push(
                String::from_utf8(spin_b)
                    .unwrap()
                    .split_whitespace()
                    .map(|x| x.parse::<f64>().unwrap() / atoms.lattice.volume)
                    .collect(),
            );
        }
        // flip the grid points as VASP outputs density[z, y, x]
        let grid_pts: [usize; 3] = [grid_vec[2], grid_vec[1], grid_vec[0]];
        Ok((voxel_origin, grid_pts, atoms, density))
    }

    /// Read atom information.
    fn to_atoms(&self, atoms_text: String) -> Atoms {
        // create regex for matching the (C|K)artesian | Direct line
        // the last match is the one we want so we don't match carbon or the comment line
        let mut lines = atoms_text.lines();
        // skip the comment line  and then read the lattice information
        let _ = lines.next();
        let mut scale = {
            lines
                .next()
                .unwrap()
                .to_string()
                .split_whitespace()
                .map(|x| x.parse::<f64>().unwrap())
                .collect::<Vec<f64>>()
        };
        // density[z, y, x] so lets swap the c and a
        let mut c = {
            lines
                .next()
                .unwrap()
                .to_string()
                .split_whitespace()
                .map(|x| x.parse::<f64>().unwrap())
                .collect::<Vec<f64>>()
        };
        let mut b = {
            lines
                .next()
                .unwrap()
                .to_string()
                .split_whitespace()
                .map(|x| x.parse::<f64>().unwrap())
                .collect::<Vec<f64>>()
        };
        let mut a = {
            lines
                .next()
                .unwrap()
                .to_string()
                .split_whitespace()
                .map(|x| x.parse::<f64>().unwrap())
                .collect::<Vec<f64>>()
        };
        let volume = {
            (c[0] * (a[1] * b[2] - a[2] * b[1])
                + c[1] * (a[2] * b[0] - a[0] * b[2])
                + c[2] * (a[0] * b[1] - a[1] * b[0]))
                .abs()
        };
        // the scale can be negative and this means that it is the volume of the cell
        // it can also be 3 values which is a multiplier for each lattice
        if scale.len() == 1 {
            if scale[0] < 0f64 {
                scale[0] /= -volume;
            }
            scale.push(scale[0]);
            scale.push(scale[0]);
        }
        for i in 0..3 {
            c[i] *= scale[2 - i];
            b[i] *= scale[2 - i];
            a[i] *= scale[2 - i];
        }
        let lattice = Lattice::new([
            [a[2], a[1], a[0]],
            [b[2], b[1], b[0]],
            [c[2], c[1], c[0]],
        ]);
        // now lets find out what type of file we are dealing with
        let dubious = lines.next().unwrap().split_whitespace();
        let total_atoms = match dubious
            .clone()
            .fold(String::new(), |acc, val| format!("{}{}", acc, val))
            .parse::<usize>()
            .is_ok()
        {
            true => dubious,
            false => lines.next().unwrap().split_whitespace(),
        }
        .fold(0, |acc, val| acc + val.parse::<usize>().unwrap());
        let mut dubious = lines.next().unwrap().trim_start().to_lowercase();
        if dubious.starts_with('s') {
            dubious = lines.next().unwrap().trim_start().to_lowercase();
        }
        let coord = if dubious.starts_with('d') {
            Coord::Fractional
        } else {
            Coord::Cartesian
        };
        let pos: Vec<[f64; 3]> = (0..total_atoms)
            .map(|_| {
                lines
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .take(3)
                    .map(|f| f.parse::<f64>().unwrap())
                    .collect::<Vec<f64>>()
                    .try_into() // we only take 3 so safe
                    .unwrap()
            })
            .collect();
        // make the positions fractional and swap c and a
        let positions = match coord {
            Coord::Fractional => pos
                .into_iter()
                .map(|p| {
                    utils::dot(
                        [
                            p[2].rem_euclid(1f64),
                            p[1].rem_euclid(1f64),
                            p[0].rem_euclid(1f64),
                        ],
                        lattice.to_cartesian,
                    )
                })
                .collect(),
            Coord::Cartesian => pos
                .into_iter()
                .map(|p| {
                    let p =
                        utils::dot([p[2], p[1], p[0]], lattice.to_fractional);
                    utils::dot(
                        [
                            p[0].rem_euclid(1f64),
                            p[1].rem_euclid(1f64),
                            p[2].rem_euclid(1f64),
                        ],
                        lattice.to_cartesian,
                    )
                })
                .collect(),
        };
        Atoms::new(lattice, positions, atoms_text)
    }

    /// Write a CHGCAR from a vector of options where None will be written as zero.
    fn write(
        &self,
        atoms: &Atoms,
        data: Vec<Option<f64>>,
        filename: String,
        visible_pbar: bool,
    ) -> std::io::Result<()> {
        let filename = format!("{}_CHGCAR", filename);
        let mut buffer = BufWriter::new(File::create(filename.clone())?);
        let length = data.len() / 5 + (data.len() % 5 != 0) as usize;
        let pbar: Box<dyn ProgressBar> = match visible_pbar {
            true => Box::new(Bar::new(
                length,
                format!("Writing file {}:", filename),
            )),
            false => Box::new(HiddenBar {}),
        };
        buffer.write_all(atoms.text.as_bytes())?;
        data.chunks(5).for_each(|line| {
            if let Err(e) = line.iter().try_for_each(|f| {
                write!(
                    buffer,
                    " {:.11}",
                    FortranFormat {
                        float: *f,
                        mult: atoms.lattice.volume
                    }
                )
            }) {
                panic!("Error occured during write: {}", e)
            };
            if let Err(e) = writeln!(buffer) {
                panic!("Error occured during write: {}", e)
            };
            pbar.tick();
        });
        Ok(())
    }

    /// Deals with fortran indexing.
    fn coordinate_format(&self, coords: [f64; 3]) -> (String, String, String) {
        let z = format!("{:.6}", coords[0]);
        let y = format!("{:.6}", coords[1]);
        let x = format!("{:.6}", coords[2]);
        (x, y, z)
    }
}
