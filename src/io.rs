use crate::arguments::{Args, Reference};
use crate::atoms::Atoms;

/// File I/O for the gaussian cube format.
pub mod cube;
/// Custom BufReader.
pub mod reader;
/// File I/O for the VASP file format.
pub mod vasp;
// Write analysis files.
pub mod output;

/// Indicates the available file types of the density file.
pub enum FileType {
    /// CHGCAR, CHG and PARCHG.
    Vasp,
    /// Guassian, CP2K etc.
    Cube,
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
    /// * `data`: The density to write to file.
    fn write(&self, atoms: &Atoms, data: Vec<Vec<f64>>);

    /// How the format the positions of maxima and atoms
    ///
    /// * `coords`: The 3d representation of the position.
    fn coordinate_format(&self, coords: [f64; 3]) -> (String, String, String);
}
