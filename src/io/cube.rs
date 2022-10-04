use crate::atoms::{Atoms, Lattice};
use crate::io::reader::BufReader;
use crate::io::{FileFormat, FortranFormat, ReadFunction};
use crate::progress::Bar;
use crate::utils;
use std::fs::File;
use std::io::{BufWriter, Read, Write};

/// Convert from chemists.
const LENGTH_UNITS: f64 = 0.52917721067;
/// Convert from chemists.
const VOLUME_UNITS: f64 = LENGTH_UNITS * LENGTH_UNITS * LENGTH_UNITS;

/// Structure for reading/writing a cube file.
pub struct Cube {}

impl FileFormat for Cube {
    /// reads a cube file from filename.
    fn read(&self, filename: String) -> ReadFunction {
        // the voxel origin in cube files is (0.5, 0.5, 0.5)
        let mut voxel_origin = [0.5f64; 3];

        println!("Reading {} as cube format:", filename);
        // find the start and end points of the density as well as the total file size
        let (start, grid_pts) = {
            // open the file in a buffer reader
            let mut reader = BufReader::open(filename.clone())?;

            let mut buffer = String::new();
            let mut pos = 0;
            // first two lines are comments
            for _ in 0..2 {
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
            let natoms = match reader.read_line(&mut buffer) {
                Some(line) => {
                    let (text, size) = line?;
                    pos += size;
                    let split = text
                        .split_whitespace()
                        .map(|x| x.parse::<f64>())
                        .collect::<Vec<Result<f64, std::num::ParseFloatError>>>(
                        );
                    if split.len() == 5 && split[4] != Ok(1.) {
                        panic!(
                            "Error(Unsuppoerted): Multiple values per voxel."
                        );
                    }
                    let natoms = match split[0] {
                        Ok(x) => x as isize,
                        Err(_) => panic!("Error: Cannot read {} as cube file.",
                                         filename),
                    };
                    for i in 0..3 {
                        voxel_origin[i] += match split[i + 1] {
                            Ok(x) => x,
                            Err(_) => {
                                panic!("Error: Cannot read {} as cube file.",
                                       filename)
                            }
                        };
                    }
                    natoms
                }
                None => panic!("Error: Cannot read {} as cube file.", filename),
            };
            if natoms < 0 {
                panic!("Error(Unsuppoerted): Multiple values per voxel.");
            }
            let mut grid_pts = [0usize; 3];
            for gp in &mut grid_pts {
                *gp = match reader.read_line(&mut buffer) {
                    Some(line) => {
                        let (text, size) = line?;
                        pos += size;
                        match text.split_whitespace().next() {
                            Some(x) => match x.parse::<usize>() {
                                Ok(x) => x,
                                Err(_) => panic!(
                                    "Error: Cannot read {} as cube file.",
                                    filename
                                ),
                            },
                            None => {
                                panic!("Error: Cannot read {} as cube file.",
                                       filename)
                            }
                        }
                    }
                    None => {
                        panic!("Error: Cannot read {} as cube file.", filename)
                    }
                }
            }
            for _ in 0..natoms.abs() {
                match reader.read_line(&mut buffer) {
                    Some(line) => {
                        let (_, size) = line?;
                        pos += size;
                    }
                    None => {
                        panic!("Error: Cannot read {} as cube file.", filename)
                    }
                }
            }
            (pos, grid_pts)
        };
        // Now we know where everything is so let's work out what to do
        // Start by making vector of start and end points of the densities
        let mut file = File::open(filename)?;
        let total = file.metadata()?.len();
        // assign Vectos with the capacity of what it is to hold
        let mut xyz_b = Vec::with_capacity(start);
        let mut density_b = Vec::with_capacity(total as usize - start);
        // read the xyz information into xyz_b
        let _ = <File as Read>::by_ref(&mut file).take(start as u64)
                                                 .read_to_end(&mut xyz_b)?;
        // read the total charge density into density_b
        let _ = <File as Read>::by_ref(&mut file).take(total - start as u64)
                                                 .read_to_end(&mut density_b)?;
        // convert the bytes we have read into a String and an Atoms struct
        let xyz = String::from_utf8(xyz_b).unwrap();
        let atoms = self.to_atoms(xyz);
        // convert out of Bohr
        let density = String::from_utf8(density_b).unwrap()
                                                  .split_whitespace()
                                                  .map(|x| {
                                                      x.parse::<f64>().unwrap()
                                                      / VOLUME_UNITS
                                                  })
                                                  .collect::<Vec<f64>>();
        println!("File read successfully.");
        Ok((voxel_origin, grid_pts, atoms, vec![density]))
    }

    /// Read atoms information from file header.
    fn to_atoms(&self, atoms_text: String) -> Atoms {
        let mut lines = atoms_text.lines();
        // skip the 2 comment lines + voxel info and then read the lattice information
        let _ = lines.next();
        let _ = lines.next();
        let _ = lines.next();
        let mut a = {
            lines.next()
                 .unwrap()
                 .to_string()
                 .split_whitespace()
                 .map(|x| x.parse::<f64>().unwrap())
                 .collect::<Vec<f64>>()
        };
        // density[z, y, x] so lets swap the c and a
        let mut b = {
            lines.next()
                 .unwrap()
                 .to_string()
                 .split_whitespace()
                 .map(|x| x.parse::<f64>().unwrap())
                 .collect::<Vec<f64>>()
        };
        let mut c = {
            lines.next()
                 .unwrap()
                 .to_string()
                 .split_whitespace()
                 .map(|x| x.parse::<f64>().unwrap())
                 .collect::<Vec<f64>>()
        };
        for i in 1..4 {
            c[i] *= c[0] * LENGTH_UNITS;
            b[i] *= b[0] * LENGTH_UNITS;
            a[i] *= a[0] * LENGTH_UNITS;
        }
        let lattice = Lattice::new([[a[1], a[2], a[3]],
                                    [b[1], b[2], b[3]],
                                    [c[1], c[2], c[3]]]);
        let mut positions: Vec<[f64; 3]> = vec![];
        // make the positions fractional and swap c and a
        for line in lines {
            let pos = line.split_whitespace()
                          .map(|x| x.parse::<f64>().unwrap() * LENGTH_UNITS)
                          .collect::<Vec<f64>>();
            let pos_frac =
                utils::dot([pos[2], pos[3], pos[4]], lattice.to_fractional)
                    .iter()
                    .map(|x| x.rem_euclid(1.))
                    .collect::<Vec<f64>>();
            let pos_cart = {
                utils::dot([pos_frac[0], pos_frac[1], pos_frac[2]],
                           lattice.to_cartesian)
            };
            positions.push(pos_cart);
        }
        Atoms::new(lattice, positions, atoms_text)
    }

    /// Write a cube file from a vector of options where None will be written as
    /// zero.
    fn write(&self,
             atoms: &Atoms,
             data: Vec<Option<f64>>,
             filename: String,
             pbar: Bar)
             -> std::io::Result<()> {
        let filename = format!("{}.cube", filename);
        let mut buffer = BufWriter::new(File::create(filename)?);
        pbar.set_length(data.len() / 6 + (data.len() % 6 != 0) as usize);
        buffer.write_all(atoms.text.as_bytes())?;
        data.chunks(6).for_each(|line| {
                          if let Err(e) = line.iter().try_for_each(|f| {
                                                         write!(
                    buffer,
                    " {:.5}",
                    FortranFormat {
                        float: *f,
                        mult: VOLUME_UNITS
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

    /// Coordinate format for dealing with fortran indexing (doesn't affect cube).
    fn coordinate_format(&self, coords: [f64; 3]) -> (String, String, String) {
        let x = format!("{:.6}", coords[0]);
        let y = format!("{:.6}", coords[1]);
        let z = format!("{:.6}", coords[2]);
        (x, y, z)
    }
}
