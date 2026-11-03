use crate::arguments::{Args, Reference};
use crate::atoms::Atoms;
use crate::io::cube::Cube;
use crate::io::vasp::Vasp;
use std::io::BufRead;

/// File I/O for the gaussian cube format.
pub mod cube;
/// Write analysis files.
pub mod output;
/// File I/O for the VASP file format.
pub mod vasp;

/// Return type of the read function in FileFormat: the voxel origin, grid
/// dimensions, atoms and densities read from a file.
pub type ReadFunction =
    std::io::Result<([f64; 3], [usize; 3], Atoms, Vec<Vec<f64>>)>;

/// Return type of the init function in FileFormat: the densities, reference
/// density (empty if there is none), atoms, grid dimensions and voxel origin.
pub type InitReturn = (Vec<Vec<f64>>, Vec<f64>, Atoms, [usize; 3], [f64; 3]);

/// FileFormat trait. Used for handling input from a file.
pub trait FileFormat {
    /// Reads the file into a [`ReadFunction`] containing the information
    /// needed from the file to build a [`Grid`](crate::grid::Grid).
    ///
    /// * `filename`: The name of the file to read.
    fn read(&self, filename: String) -> ReadFunction;

    /// Reads the non-density section of the file into an [`Atoms`] object.
    ///
    /// * `atom_text`: The full string of non-density information from the
    ///   density file.
    fn to_atoms(&self, atom_text: String) -> Atoms;

    /// Writes a density, `data`, to file in the correct format.
    ///
    /// * `atoms`: The associated &[`Atoms`] object for the density file.
    /// * `data`: The density to write to file.
    /// * `filename`: Where to save the file, minus any suffix as this should
    ///   be applied in the function.
    /// * `visible_pbar`: Whether to show a progress bar for the write.
    fn write(
        &self,
        atoms: &Atoms,
        data: Vec<f64>,
        filename: String,
        visible_pbar: bool,
    ) -> std::io::Result<()>;

    /// How to format the positions of maxima and atoms for the file format.
    ///
    /// * `coords`: The 3d representation of the position.
    fn coordinate_format(&self, coords: [f64; 3]) -> [f64; 3];

    /// Reads the density file and any spin or reference files from the
    /// arguments, checking they all share a grid.
    ///
    /// * `args`: [`Args`] parsed from the command line.
    ///
    /// # Panics
    /// Panics if a file cannot be read or the grids do not match.
    fn init(&self, args: &Args) -> InitReturn {
        let (voxel_origin, grid, atoms, mut densities) = self
            .read(args.file.clone())
            .unwrap_or_else(|e| panic!("Error: Problem reading file.\n{}", e));
        // reads an additional file, checking it has the same grid
        let read_extra = |filename: &String, kind: &str| {
            let (_, g, _, d) =
                self.read(filename.clone()).unwrap_or_else(|e| {
                    panic!("Error: Problem reading {}.\n{}", filename, e)
                });
            assert_eq!(
                g, grid,
                "Error: {} density has different grid size.",
                kind
            );
            d
        };
        if let Some(spin_file) = &args.spin {
            if densities.len() != 1 {
                panic!(
                    "Number of densities in original file is not 1.\nAmbiguous how to handle new spin when {} already has {} spin densities.",
                    args.file,
                    densities.len() - 1
                );
            }
            let mut spin = read_extra(spin_file, "Spin");
            if spin.len() != 1 {
                panic!(
                    "Number of densities in original file is not 1.\nAmbiguous how to handle spin density when {} contains {} densities.",
                    spin_file,
                    spin.len()
                );
            }
            densities.push(spin.remove(0));
        }
        let rho = match &args.reference {
            Reference::None => Vec::new(),
            Reference::One(f) => read_extra(f, "Reference").remove(0),
            Reference::Two(f1, f2) => {
                let mut rho = read_extra(f1, "Reference").remove(0);
                let rho2 = read_extra(f2, "Reference").remove(0);
                rho.iter_mut().zip(rho2).for_each(|(a, b)| *a += b);
                rho
            }
        };
        (densities, rho, atoms, grid, voxel_origin)
    }
}

/// Indicates the available file types of the density file.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FileType {
    /// CHGCAR, CHG and PARCHG.
    Vasp(Vasp),
    /// Gaussian, CP2K etc.
    Cube(Cube),
}

impl FileFormat for FileType {
    fn read(&self, filename: String) -> ReadFunction {
        match self {
            Self::Vasp(f) => f.read(filename),
            Self::Cube(f) => f.read(filename),
        }
    }

    fn to_atoms(&self, atom_text: String) -> Atoms {
        match self {
            Self::Vasp(f) => f.to_atoms(atom_text),
            Self::Cube(f) => f.to_atoms(atom_text),
        }
    }

    fn write(
        &self,
        atoms: &Atoms,
        data: Vec<f64>,
        filename: String,
        visible_pbar: bool,
    ) -> std::io::Result<()> {
        match self {
            Self::Vasp(f) => f.write(atoms, data, filename, visible_pbar),
            Self::Cube(f) => f.write(atoms, data, filename, visible_pbar),
        }
    }

    fn coordinate_format(&self, coords: [f64; 3]) -> [f64; 3] {
        match self {
            Self::Vasp(f) => f.coordinate_format(coords),
            Self::Cube(f) => f.coordinate_format(coords),
        }
    }
}

/// What type of density to write.
pub enum WriteType {
    /// Write the Bader volumes of these atoms, all atoms if empty.
    Atom(Vec<usize>),
    /// Don't write anything.
    None,
}

/// Reads a line from `reader` into `buffer`, including the newline.
///
/// Returns `None` at the end of the file and an error if the file ends
/// before `expected` lines have been read.
pub(crate) fn read_line<'a>(
    reader: &mut impl BufRead,
    buffer: &'a mut String,
) -> std::io::Result<Option<&'a str>> {
    buffer.clear();
    match reader.read_line(buffer)? {
        0 => Ok(None),
        _ => Ok(Some(buffer.as_str())),
    }
}

/// An error for a file that does not match the expected format.
pub(crate) fn invalid_data(filename: &str, format: &str) -> std::io::Error {
    std::io::Error::new(
        std::io::ErrorKind::InvalidData,
        format!("Cannot read {} as a {} file.", filename, format),
    )
}

/// Parses every whitespace separated number on a line.
pub(crate) fn parse_line<T: std::str::FromStr>(line: &str) -> Option<Vec<T>> {
    line.split_whitespace()
        .map(|x| x.parse::<T>().ok())
        .collect()
}

/// Turn a float into fortran "scientific" notation (leading digit is zero).
///
/// The precision of the formatter sets the number of digits in the mantissa.
///
/// # Example
/// ```
/// use bader::io::FortranFormat;
///
/// assert_eq!(format!("{:.5}", FortranFormat::new(12.5)), " 0.12500E+02");
/// assert_eq!(format!("{:.5}", FortranFormat::new(-0.05)), "-0.50000E-01");
/// assert_eq!(format!("{:.3}", FortranFormat::new(0.0)), " 0.000E+00");
/// ```
pub struct FortranFormat(f64);

impl FortranFormat {
    /// Format `float`.
    pub fn new(float: f64) -> Self {
        Self(float)
    }
}

impl std::fmt::Display for FortranFormat {
    fn fmt(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        let prec = formatter.precision().unwrap_or(6);
        let sign = if self.0.is_sign_negative() { '-' } else { ' ' };
        let float = self.0.abs();
        if float == 0.0 || !float.is_finite() {
            return write!(formatter, "{}0.{:0<prec$}E+00", sign, 0);
        }
        // normalise so the mantissa is 0.d1d2...dn
        let mut exponent = float.log10().floor() as i32 + 1;
        let mut mantissa =
            (float * 10f64.powi(prec as i32 - exponent)).round() as u64;
        // rounding can carry into the next power of ten
        if mantissa >= 10u64.pow(prec as u32) {
            mantissa /= 10;
            exponent += 1;
        }
        write!(formatter, "{}0.{:0>prec$}E{:+03}", sign, mantissa, exponent)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fortran_format() {
        let f = |x: f64, p: usize| format!("{:.p$}", FortranFormat::new(x));
        assert_eq!(f(152.46, 11), " 0.15246000000E+03");
        assert_eq!(f(1.0, 5), " 0.10000E+01");
        assert_eq!(f(0.5, 5), " 0.50000E+00");
        assert_eq!(f(0.05, 5), " 0.50000E-01");
        assert_eq!(f(-0.10283642961E-07, 11), "-0.10283642961E-07");
        assert_eq!(f(-3.0, 2), "-0.30E+01");
        assert_eq!(f(0.0, 4), " 0.0000E+00");
        assert_eq!(f(-0.0, 4), "-0.0000E+00");
        // rounding up to the next power of ten
        assert_eq!(f(9.99999999, 5), " 0.10000E+02");
        assert_eq!(f(0.099999999, 3), " 0.100E+00");
    }

    #[test]
    fn test_parse_line() {
        assert_eq!(
            parse_line::<f64>(" 1.0 2.5\t-3e2 \n"),
            Some(vec![1.0, 2.5, -300.0])
        );
        assert_eq!(parse_line::<usize>("1 2 x"), None);
        assert_eq!(parse_line::<usize>(""), Some(vec![]));
    }
}
