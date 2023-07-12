use crate::atoms::{Atoms, Lattice};
use crate::io::reader::BufReader;
use crate::io::{FileFormat, FortranFormat, ReadFunction};
use crate::progress::Bar;
use crate::utils;
use regex::{Regex, RegexSet};
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
        println!("Reading {} as VASP format:", filename);
        // find the start and end points of the density as well as the total file size
        let (grid, aug, total) = {
            // open the file in a buffer reader
            let mut reader = BufReader::open(filename.clone())?;

            let mut buffer = String::new();
            let mut grid: Vec<[usize; 2]> = vec![];
            let mut aug: Vec<usize> = vec![];
            let mut pos = 0;
            // search for the grid lines or augmentation that bound the densities
            let regex = RegexSet::new([r"^\s*\d+\s+\d+\s+\d+\s*$",
                                        r"^\s*aug"]).unwrap();
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
            // lets start trying to match
            while let Some(line) = reader.read_line(&mut buffer) {
                let (text, size) = line?;
                if regex.is_match(text) {
                    let matches: Vec<usize> =
                        regex.matches(text).into_iter().collect();
                    let matches = matches[0];
                    match matches {
                        0 => {
                            let start = pos;
                            pos += size;
                            let end = pos;
                            grid.push([start, end]);
                        }
                        1 => {
                            aug.push(pos);
                            pos += size;
                        }
                        _ => {}
                    }
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
        let _ = <File as Read>::by_ref(&mut file).take(grid[0][1] as u64)
                                                 .read_to_end(&mut poscar_b)?;
        file.seek(SeekFrom::Current(grid[0][0] as i64 - grid[0][1] as i64))?;
        // read the grid line into grid_pts_b
        let _ =
            <File as Read>::by_ref(&mut file).take((grid[0][1] - grid[0][0])
                                                   as u64)
                                             .read_to_end(&mut grid_pts_b)?;
        // read the total charge density into density_b
        let _ = <File as Read>::by_ref(&mut file).take((stop[0] - start[0])
                                                       as u64)
                                                 .read_to_end(&mut density_b)?;
        // convert the bytes we have read into a String and an Atoms struct
        let poscar = String::from_utf8(poscar_b).unwrap();
        let grid_vec: Vec<usize> = {
            String::from_utf8(grid_pts_b).unwrap()
                                         .split_whitespace()
                                         .map(|x| x.parse::<usize>().unwrap())
                                         .collect()
        };
        let atoms = self.to_atoms(poscar);
        // convert out of VASP's strange units
        density.push(String::from_utf8(density_b).unwrap()
                                                 .split_whitespace()
                                                 .map(|x| {
                                                     x.parse::<f64>().unwrap()
                                                     / atoms.lattice.volume
                                                 })
                                                 .collect());
        for i in 1..start.len() {
            let mut spin_b = Vec::with_capacity(stop[i] - start[i]);
            let _ = file.seek(SeekFrom::Start(start[i] as u64));
            let _ = <File as Read>::by_ref(&mut file).take((stop[i] - start[i])
                                                           as u64)
                                                     .read_to_end(&mut spin_b)?;
            density.push(String::from_utf8(spin_b).unwrap()
                                                  .split_whitespace()
                                                  .map(|x| {
                                                      x.parse::<f64>().unwrap()
                                                      / atoms.lattice.volume
                                                  })
                                                  .collect());
        }
        // flip the grid points as VASP outputs density[z, y, x]
        let grid_pts: [usize; 3] = [grid_vec[2], grid_vec[1], grid_vec[0]];
        println!("File read successfully.");
        Ok((voxel_origin, grid_pts, atoms, density))
    }

    /// Read atom information.
    fn to_atoms(&self, atoms_text: String) -> Atoms {
        // create regex for matching the (C|K)artesian | Direct line
        let coord_regex = Regex::new(r"(?m)^\s*(c|C|k|K|d|D)\w*").unwrap();
        // the last match is the one we want so we don't match carbon or the comment line
        let matches = coord_regex.find_iter(&atoms_text).last().unwrap();
        let coord_text =
            &atoms_text[matches.start()..matches.end()].trim_start();
        let coord = if coord_text.starts_with('d') | coord_text.starts_with('D')
        {
            Coord::Fractional
        } else {
            Coord::Cartesian
        };
        let mut pos: Vec<f64> = vec![];
        // push the floats to the pos vector so we don't get caught out by selective dynamics
        let _ = {
            &atoms_text[matches.end()..].to_string()
                                        .split_whitespace()
                                        .map(|x| {
                                            if let Ok(x) = x.parse::<f64>() {
                                                pos.push(x);
                                            }
                                        })
                                        .collect::<Vec<_>>()
        };
        for _ in 0..3 {
            pos.pop();
        }
        let mut lines = atoms_text.lines();
        // skip the comment line  and then read the lattice information
        let _ = lines.next();
        let mut scale = {
            lines.next()
                 .unwrap()
                 .to_string()
                 .split_whitespace()
                 .map(|x| x.parse::<f64>().unwrap())
                 .collect::<Vec<f64>>()
        };
        // density[z, y, x] so lets swap the c and a
        let mut c = {
            lines.next()
                 .unwrap()
                 .to_string()
                 .split_whitespace()
                 .map(|x| x.parse::<f64>().unwrap())
                 .collect::<Vec<f64>>()
        };
        let mut b = {
            lines.next()
                 .unwrap()
                 .to_string()
                 .split_whitespace()
                 .map(|x| x.parse::<f64>().unwrap())
                 .collect::<Vec<f64>>()
        };
        let mut a = {
            lines.next()
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
        let lattice = Lattice::new([[a[2], a[1], a[0]],
                                    [b[2], b[1], b[0]],
                                    [c[2], c[1], c[0]]]);
        let mut positions: Vec<[f64; 3]> = vec![];
        // make the positions fractional and swap c and a
        match coord {
            Coord::Fractional => {
                for i in (0..pos.len()).step_by(3) {
                    positions.push(utils::dot([pos[i + 2].rem_euclid(1f64),
                                               pos[i + 1].rem_euclid(1f64),
                                               pos[i].rem_euclid(1f64)],
                                              lattice.to_cartesian));
                }
            }
            Coord::Cartesian => {
                for i in (0..pos.len()).step_by(3) {
                    let p = utils::dot([pos[i + 2], pos[i + 1], pos[i]],
                                       lattice.to_fractional);
                    positions.push(utils::dot([p[0].rem_euclid(1f64),
                                               p[1].rem_euclid(1f64),
                                               p[2].rem_euclid(1f64)],
                                              lattice.to_cartesian));
                }
            }
        }
        Atoms::new(lattice, positions, atoms_text)
    }

    /// Write a CHGCAR from a vector of options where None will be written as zero.
    fn write(&self,
             atoms: &Atoms,
             data: Vec<Option<f64>>,
             filename: String,
             pbar: Bar)
             -> std::io::Result<()> {
        let filename = format!("{}_CHGCAR", filename);
        let mut buffer = BufWriter::new(File::create(filename)?);
        pbar.set_length(data.len() / 5 + (data.len() % 5 != 0) as usize);
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
                                                     })
                          {
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
