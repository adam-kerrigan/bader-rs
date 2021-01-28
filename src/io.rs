use crate::arguments::{Args, Reference};
use crate::atoms::Atoms;
use crate::progress::Bar;

/// File I/O for the gaussian cube format.
pub mod cube;
/// Write analysis files.
pub mod output;
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

/// What type of density to write.
pub enum WriteType {
    /// Write a Bader Atom.
    Atom(Vec<usize>),
    /// Write a Bader Volume.
    Volume(Vec<usize>),
    /// Don't write anything.
    None,
}

/// Turn a float into fortran "scientific" notation (leading digit is zero).
pub struct FortranFormat {
    /// The float to convert to a string. Wrapped in an option as we need to log
    /// so `0f64` should be stored as None.
    float: Option<f64>,
    /// A value to multiply the float by before printing, eg. a volume.
    mult: f64,
}

impl std::fmt::Display for FortranFormat {
    /// Format the structure into a fortran style exponential.
    fn fmt(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        let prec = if let Some(prec) = formatter.precision() {
            prec
        } else {
            6
        };
        match self.float {
            None => {
                write!(formatter, " 0.{:0<width$}E{:+03}", 0, 0, width = prec)
            }
            Some(f) => {
                let float = f * self.mult;
                let exponant = float.log10() as i32 + 1;
                let decimals = float.abs() * 10f64.powi(prec as i32 - exponant);
                let decimals = decimals.round() as usize;
                if float.is_sign_negative() {
                    write!(formatter,
                           "-0.{:0<width$}E{:+03}",
                           decimals,
                           exponant,
                           width = prec)
                } else {
                    write!(formatter,
                           " 0.{:0<width$}E{:+03}",
                           decimals,
                           exponant,
                           width = prec)
                }
            }
        }
    }
}

/// Return type of the read function in FileFormat.
pub type ReadFunction =
    std::io::Result<([f64; 3], [usize; 3], Atoms, Vec<Vec<f64>>)>;
/// Return type of the init function in FileFormat.
type InitReturn = (Vec<Vec<f64>>, Vec<f64>, Atoms, [usize; 3], [f64; 3]);

/// FileFormat trait. Used for handling input from a file.
pub trait FileFormat {
    /// Returns the parts required to build [`Grid`] and [`Atoms`] structures.
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
                densities[0].iter()
                            .zip(&densities2[0])
                            .map(|(a, b)| a + b)
                            .collect::<Vec<f64>>()
            }
        };
        (densities, rho, atoms, grid, voxel_origin)
    }

    /// Reads the file into a [`ReadFunction`] containing the information
    /// needed from the file to build a [`Grid`].
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
    /// * `data`: The density to write to file wrapped in options with None representing 0.
    /// * `filename`: Where to save the file, minus any suffix as this should
    /// be applied in the function.
    /// * `pbar`: A progress bar for monitoring the write.
    fn write(&self,
             atoms: &Atoms,
             data: Vec<Option<f64>>,
             filename: String,
             pbar: Bar)
             -> std::io::Result<()>;

    /// How the format the positions of maxima and atoms
    ///
    /// * `coords`: The 3d representation of the position.
    fn coordinate_format(&self, coords: [f64; 3]) -> (String, String, String);
}
