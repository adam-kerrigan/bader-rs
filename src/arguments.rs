use crate::io::{self, ReadFunction};
use crate::methods::{self, StepMethod};
use clap::{crate_authors, App, Arg};

/// Indicates how many reference files are passed
#[derive(Clone)]
pub enum Reference {
    One(String),
    Two(String, String),
    None,
}

/// Holds the arguments passed to the program from the command-line
pub struct Args {
    pub file: String,
    pub read: ReadFunction,
    pub method: StepMethod,
    pub reference: Reference,
    pub vacuum_tolerance: Option<f64>,
}

impl Args {
    /// Initialises the structure from the command-line arguments.
    pub fn new() -> Self {
        let arguments = App::new("Multi-threaded Bader Charge Analysis")
            .author(crate_authors!())
            .version("0.1.0")
            .arg(Arg::with_name("file")
                .required(true)
                .index(1))
            .arg(Arg::with_name("method")
                .short('m')
                .long("method")
                .takes_value(true)
                .possible_value("ongrid")
                .possible_value("neargrid")
                .case_insensitive(false)
                .about("method by which to partition the charge density")
                .long_about(
"Use the \"near-grid\" or \"on-grid\" methods based on the algorithms presented
in W. Tang et al. A grid-based Bader analysis algorithm without lattice bias,
J. Phys.: Condens. Matter 21, 084204 (2009)"))
            .arg(Arg::with_name("file type")
                .short('t')
                .long("type")
                .takes_value(true)
                .possible_value("cube")
                .possible_value("chgcar")
                .case_insensitive(false)
                .about("the file type of the charge density")
                .long_about(
"The file type of the input file. If this is not supplied the type will attempt
to be infered from the filename"))
            .arg(Arg::with_name("reference")
                .short('r')
                .long("ref")
                .multiple(true)
                .max_values(2)
                .number_of_values(1)
                .about("file(s) containing reference charge")
                .long_about(
"A reference charge to do the partitioning upon. Two files can be passed
by using multiple flags (bader CHGCAR -r AECCAR0 -r AECCAR2). If two files are
passed they are summed together."))
            .arg(Arg::with_name("all electron")
                .short('a')
                .long("aec")
                .about("convience flag for reading both aeccars")
                .takes_value(false)
                .multiple(false)
                .conflicts_with("reference"))
            .arg(Arg::with_name("vacuum tolerance")
                .short('v')
                .long("vac")
                .takes_value(true)
                .about("cut-off at which charge is considered vacuum")
                .long_about(
"Values of density below the supplied value are considered vacuum and are not
included in the calculation. A value of \"auto\" can be passed to use 1E-3 C*m^-3"))
            .get_matches();

        let file = match arguments.value_of("file") {
            Some(f) => String::from(f),
            None => String::new(),
        };
        let file_type = match arguments.value_of("file type") {
            Some(f) => Some(String::from(f)),
            None => None,
        };
        let method = match arguments.value_of("method") {
            Some("neargrid") => methods::neargrid,
            Some("ongrid") => methods::ongrid,
            Some(_) => methods::neargrid,
            None => methods::neargrid,
        };
        let vacuum_tolerance = match arguments.value_of("vacuum tolerance") {
            Some(s) => {
                if s.eq("auto") {
                    Some(1E-3)
                } else {
                    match s.parse::<f64>() {
                        Ok(x) => Some(x),
                        Err(_) => {
                            panic!("Couldn't parse vacuum tolerance: {}", s)
                        }
                    }
                }
            }
            None => None,
        };
        let references: Vec<_> = if arguments.is_present("all electron") {
            vec!["AECCAR0", "AECCAR2"]
        } else {
            match arguments.values_of("reference") {
                Some(x) => x.collect(),
                None => Vec::with_capacity(0),
            }
        };
        let reference = match references.len() {
            2 => Reference::Two(String::from(references[0]),
                                String::from(references[1])),
            1 => Reference::One(String::from(references[0])),
            _ => Reference::None,
        };
        let read: ReadFunction = match file_type {
            Some(s) => {
                if s.eq("vasp") {
                    io::vasp::read
                } else {
                    io::vasp::read
                }
            }
            None => {
                if file.to_lowercase().contains("car") {
                    io::vasp::read
                } else {
                    println!("Error: File-type cannot be infered, attempting to read as VASP");
                    io::vasp::read
                }
            }
        };
        return Self { file,
                      read,
                      method,
                      reference,
                      vacuum_tolerance };
    }
}
