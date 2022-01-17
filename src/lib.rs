//! An incredibly fast, multi-threaded, Bader charge partitioning binary and
//! library. Based on methods presented in
//! [Yu Min  and Trinkle Dallas R. 2011  J. Che.m Phys. 134 064111] and with
//! adaptions for multi-threading.
//!
//! ### Supported Platforms
//! - Linux
//! - Os X
//! - Windows
//!
//! ## Installing the binary
//! ### Cargo
//! ```sh
//! $ cargo install bader
//! ```
//! ### From Source
//! To check out the lastest features not in the binaries yet you can compile
//! from source. To do this run the following, which will create the
//! ./target/release/bca executable.
//! ```sh
//! $ git clone https://github.com/adam-kerrigan/bader-rs
//! $ cd bader-rs
//! $ cargo build --verbose --release
//! ```
//! From here you can either move or link the binary to folder in your path.
//! ```sh
//! $ mv ./target/release/bca ~/bin
//! ```
//!
//! ## Using the library
//! Add the following to your Cargo.toml:
//! `bader = "0.4.1"`
//!
//! ### Minimum Supported Rust Version (MSRV)
//! This crate is guaranteed to compile on stable Rust 1.54.0 and up. It *might*
//! compile with older versions but that may change in any new patch release.
//! ## Usage
//! The program takes a charge density file as input and performs Bader analysis
//! of the data. Currently it supports density in [VASP] or [cube] formats. It
//! is recommended to run VASP calculations with [LAECHG] = .TRUE. to print the
//! core density and self-consistent valence density. These can then be passed
//! as reference files to the program using the -r, --reference flag where they
//! will be summed.
//! ```sh
//! $ bca CHGCAR -r AECCAR0 -r AECCAR2
//! ```
//! VASP charge density files containing spin densities will output the the
//! partitioned spin also. To achieve this for cube files requires using the
//! --spin flag to pass a second file to treat as the spin density.
//! ```sh
//! $ bca charge-density.cube -s spin-density.cube
//! ```
//! For a detailed list of usage options run
//! ```sh
//! $ bca --help
//! ```
//! ## Output
//! The program outputs two files, ACF.dat & BCF.dat. The Atomic Charge File
//! (ACF.dat) contians the charge (and spin) information for each atom and the
//! Bader Charge File (BCF.dat) contains the information about each Bader volume.
//! The BCF file also includes the atom number in the number column formatted as
//! 'atom number: bader volume'.
//! ## License
//! MIT
//!
//! [//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)
//!
//! [VASP]: <https://www.vasp.at/>
//! [cube]: <https://gaussian.com/>
//! [LAECHG]: <https://www.vasp.at/wiki/index.php/LAECHG>
//! [Yu Min  and Trinkle Dallas R. 2011  J. Che.m Phys. 134 064111]: <https://doi.org/10.1063/1.3553716>
//! [cargo]: <https://doc.rust-lang.org/cargo/getting-started/installation.html>

/// Performs analysis of the VoxelMap to find the partitioned charge, assigned atom and other
/// relevant properties.
pub mod analysis;
/// Builds the [clap::App] and parses command-line arguments.
pub mod arguments;
/// Contains [Atoms](atoms::Atoms) for storing the relevant data on the atoms
/// in the calculation. Also contains [Lattice](atoms::Lattice) and
/// [ReducedLattice](atoms::ReducedLattice) for storing information about the
/// cell in which the density is stored.
pub mod atoms;
/// Contains [Grid](grid::Grid) for managing the movement around the grid on
/// which the density is stored.
pub mod grid;
/// Handles the File I/O for both the density file and result files.
/// Provides a [FileFormat](io::FileFormat) trait to be implemented by modules designed to
/// cover a specific file format of a density file.
pub mod io;
/// Contains the three methods for partioning the density, ([Ongrid](methods::ongrid),
/// [Neargrid](methods::neargrid), and [Weight](methods::weight)), and functions for
/// performing a step for in each.
pub mod methods;
/// Provides [Bar](progress::Bar): A quicker thread-safe version of the [indicatif::ProgressBar].
pub mod progress;
/// Misc functions mainly for vector and matrix manipulation.
pub mod utils;
/// Calculates the Voronoi vectors, and their alpha values for the weight method,
/// for lattices. Also useful for periodic minimum distances.
pub mod voronoi;
/// Provides the [VoxelMap](voxel_map::VoxelMap) for storing the maxima and weights of
/// partioned voxels.
pub mod voxel_map;
