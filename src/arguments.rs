use crate::io::{self, ReadFunction};
use clap::{crate_authors, App, Arg, ArgMatches};

/// Indicates how many reference files are passed
#[derive(Clone)]
pub enum Reference {
    One(String),
    Two(String, String),
    None,
}

// Indicates which method to use
#[derive(Clone)]
pub enum Method {
    OnGrid,
    NearGrid,
    Weight,
}

#[derive(Clone)]
pub enum Weight {
    Atoms,
    Volumes,
    None,
}
/// Create a container for dealing with clap and being able to test arg parsing
pub enum ClapApp {
    App,
}

impl ClapApp {
    /// Create and return the clap::App
    pub fn get(&self) -> App {
        App::new("Multi-threaded Bader Charge Analysis")
            .author(crate_authors!())
            .version("0.1.0")
            .arg(Arg::new("file")
                .required(true)
                .index(1)
                .about("The file to analyse."))
            .arg(Arg::new("method")
                .short('m')
                .long("method")
                .takes_value(true)
                .possible_value("ongrid")
                .possible_value("neargrid")
                .possible_value("weight")
                .case_insensitive(false)
                .about("method by which to partition the charge density")
                .long_about(
"Use the \"near-grid\" or \"on-grid\" methods based on the algorithms presented
in W. Tang et al. A grid-based Bader analysis algorithm without lattice bias,
J. Phys.: Condens. Matter 21, 084204 (2009). Or the \"weight\" method (Default)
presented in Min Yu and Dallas R. Trinkle. Accurate and efficient algorithm for
Bader charge integration, J. Chem. Phys. 134, 064111 (2011)."))
            .arg(Arg::new("weight")
                .short('w')
                .long("weight")
                .takes_value(true)
                .possible_value("atoms")
                .possible_value("volumes")
                .possible_value("none")
                .default_value("volumes")
                .case_insensitive(false)
                .about("How to calculate the weight of the boundary voxels")
                .long_about(
"Determines the condition for a voxel to be considered at boundary. When multiple
Bader volumes are assigned to the same atom it can be faster to only consider
voxels that neighbour a different atom to be at a boundary. Selecting none
disregards the weighting of voxels and their full charge is assigned to their
Bader volume."))
            .arg(Arg::new("file type")
                .short('t')
                .long("type")
                .takes_value(true)
                .possible_value("cube")
                .possible_value("vasp")
                .case_insensitive(false)
                .about("The file type of the charge density.")
                .long_about(
"The file type of the input file. If this is not supplied the type will attempt
to be infered from the filename."))
            .arg(Arg::new("reference")
                .short('r')
                .long("ref")
                .multiple(true)
                .max_values(2)
                .number_of_values(1)
                .about("File(s) containing reference charge.")
                .long_about(
"A reference charge to do the partitioning upon. Two files can be passed
by using multiple flags (bader CHGCAR -r AECCAR0 -r AECCAR2). If two files are
passed they are summed together."))
            .arg(Arg::new("spin")
                .short('s')
                .long("spin")
                .number_of_values(1)
                .about("File containing spin density.")
                .long_about(
"A path to the spin density associated with the original file. This is primarily
for cube files as if spin density exists in a CHGCAR it will be read automatically.
If using with VASP outputs then the files for charge and spin density must only
contain a single density (ie. the original file has been split)."))
            .arg(Arg::new("all electron")
                .short('a')
                .long("aec")
                .about("Convience flag for reading both aeccars.")
                .takes_value(false)
                .multiple(false)
                .conflicts_with("reference"))
            .arg(Arg::new("vacuum tolerance")
                .short('v')
                .long("vac")
                .takes_value(true)
                .about("Cut-off at which charge is considered vacuum.")
                .long_about(
"Values of density below the supplied value are considered vacuum and are not
included in the calculation. A value of \"auto\" can be passed to use 1E-3 C*m^-3."))
            .arg(Arg::new("threads")
                .short('J')
                .long("threads")
                .takes_value(true)
                .default_value("0")
                .about("Number of threads to distribute the calculation over.")
                .long_about(
"The number of threads to be used by the program. A default value of 0 is used
to allow the program to best decide how to use the available hardware."))
    }
}

/// Holds the arguments passed to the program from the command-line
pub struct Args {
    pub file: String,
    pub read: ReadFunction,
    pub method: Method,
    pub weight: Weight,
    pub reference: Reference,
    pub spin: Option<String>,
    pub threads: usize,
    pub vacuum_tolerance: Option<f64>,
    pub zyx_format: bool,
}

impl Args {
    /// Initialises the structure from the command-line arguments.
    pub fn new(arguments: ArgMatches) -> Self {
        let file = match arguments.value_of("file") {
            Some(f) => String::from(f),
            None => String::new(),
        };
        let file_type = match arguments.value_of("file type") {
            Some(f) => Some(String::from(f)),
            None => None,
        };
        let mut zyx_format = false;
        let read: ReadFunction = match file_type {
            Some(ftype) => {
                if ftype.eq("cube") {
                    io::cube::read
                } else {
                    zyx_format = true;
                    io::vasp::read
                }
            }
            None => {
                if file.to_lowercase().contains("cube") {
                    io::cube::read
                } else if file.to_lowercase().contains("car") {
                    zyx_format = true;
                    io::vasp::read
                } else {
                    println!("Error: File-type cannot be infered, attempting to read as VASP");
                    zyx_format = true;
                    io::vasp::read
                }
            }
        };
        let weight = match arguments.value_of("weight") {
            Some("atoms") => Weight::Atoms,
            Some("volumes") => Weight::Volumes,
            Some("none") => Weight::None,
            _ => Weight::Volumes,
        };
        let method = match arguments.value_of("method") {
            Some("neargrid") => Method::NearGrid,
            Some("weight") => Method::Weight,
            Some("ongrid") => Method::OnGrid,
            _ => match arguments.value_of("weight") {
                Some("none") => Method::NearGrid,
                _ => Method::Weight,
            },
        };
        // safe to unwrap as threads has a default value of 0
        let threads =
            match arguments.value_of("threads").unwrap().parse::<usize>() {
                Ok(x) => x,
                Err(e) => panic!("Couldn't parse threads into integer:\n{}", e),
            };
        let vacuum_tolerance = match arguments.value_of("vacuum tolerance") {
            Some(s) => {
                if s.eq("auto") {
                    Some(1E-3)
                } else {
                    match s.parse::<f64>() {
                        Ok(x) => Some(x),
                        Err(e) => {
                            panic!("Couldn't parse vacuum tolerance into float:\n{}", e)
                        }
                    }
                }
            }
            None => None,
        };
        let references: Vec<_> = if arguments.is_present("all electron") {
            if read as usize != io::vasp::read as usize {
                panic!("Error: Cannot use AECCAR flag for non VASP file-types.")
            }
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
        let spin = match arguments.value_of("spin") {
            Some(x) => Some(String::from(x)),
            None => None,
        };
        Self { file,
               read,
               method,
               weight,
               reference,
               threads,
               spin,
               vacuum_tolerance,
               zyx_format }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clapapp_get() {
        let app = ClapApp::App.get();
        assert_eq!(app.get_name(), "Multi-threaded Bader Charge Analysis")
    }

    #[test]
    fn argument_file() {
        let app = ClapApp::App.get();
        let matches = app.get_matches_from(vec!["bader", "CHGCAR"]);
        let args = Args::new(matches);
        assert_eq!(args.file, String::from("CHGCAR"));
    }

    #[test]
    #[should_panic]
    fn argument_no_file() {
        let app = ClapApp::App.get();
        let _ = app.try_get_matches_from(vec!["bader"])
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_method_ongrid() {
        let app = ClapApp::App.get();
        let matches =
            app.get_matches_from(vec!["bader", "CHGCAR", "-m", "ongrid"]);
        let args = Args::new(matches);
        match args.method {
            Method::OnGrid => (),
            _ => panic!("OnGrid passed but didnt get OnGrid"),
        }
    }

    #[test]
    fn argument_method_neargrid() {
        let app = ClapApp::App.get();
        let matches = app.get_matches_from(vec!["bader", "CHGCAR",
                                                "--method", "neargrid"]);
        let args = Args::new(matches);
        match args.method {
            Method::NearGrid => (),
            _ => panic!("NearGrid passed but didnt get NearGrid"),
        }
    }

    #[test]
    fn argument_method_default_no_weight() {
        let app = ClapApp::App.get();
        let matches =
            app.get_matches_from(vec!["bader", "CHGCAR", "-w", "none"]);
        let args = Args::new(matches);
        match args.method {
            Method::NearGrid => (),
            _ => panic!("No argument passed, didnt get NearGrid"),
        }
    }

    #[test]
    fn argument_method_default() {
        let app = ClapApp::App.get();
        let matches = app.get_matches_from(vec!["bader", "CHGCAR"]);
        let args = Args::new(matches);
        match args.method {
            Method::Weight => (),
            _ => panic!("No argument passed, didnt get Weight"),
        }
    }

    #[test]
    #[should_panic]
    fn argument_method_not_method() {
        let app = ClapApp::App.get();
        let _ = app.try_get_matches_from(vec!["bader", "CHGCAR", "-m",
                                              "ngrid"])
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_file_type_default_vasp() {
        let app = ClapApp::App.get();
        let matches = app.get_matches_from(vec!["bader", "CHGCAR"]);
        let args = Args::new(matches);
        let read: ReadFunction = io::vasp::read;
        assert_eq!(args.read as usize, read as usize);
    }

    #[test]
    fn argument_file_type_default_unknown() {
        let app = ClapApp::App.get();
        let matches = app.get_matches_from(vec!["bader", "CHG"]);
        let args = Args::new(matches);
        let read: ReadFunction = io::vasp::read;
        assert_eq!(args.read as usize, read as usize);
    }

    #[test]
    fn argument_file_type_vasp() {
        let app = ClapApp::App.get();
        let matches =
            app.get_matches_from(vec!["bader", "CHGCAR", "-t", "vasp"]);
        let args = Args::new(matches);
        let read: ReadFunction = io::vasp::read;
        assert_eq!(args.read as usize, read as usize);
    }

    #[test]
    fn argument_file_type_default_cube() {
        let app = ClapApp::App.get();
        let matches = app.get_matches_from(vec!["bader", "charge.cube"]);
        let args = Args::new(matches);
        let read: ReadFunction = io::cube::read;
        assert_eq!(args.read as usize, read as usize);
    }

    #[test]
    fn argument_file_type_cube() {
        let app = ClapApp::App.get();
        let matches = app.get_matches_from(vec!["bader",
                                                "charge.cube",
                                                "--type",
                                                "cube"]);
        let args = Args::new(matches);
        let read: ReadFunction = io::cube::read;
        assert_eq!(args.read as usize, read as usize);
    }

    #[test]
    #[should_panic]
    fn argument_file_type_not_type() {
        let app = ClapApp::App.get();
        let _ = app.try_get_matches_from(vec!["bader", "CHGCAR", "-t", "basp"])
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_spin() {
        let app = ClapApp::App.get();
        let matches =
            app.get_matches_from(vec!["bader", "density.cube", "-s", "spin.cube"]);
        let args = Args::new(matches);
        assert_eq!(args.spin, Some(String::from("spin.cube")))
    }

    #[test]
    fn argument_reference_one() {
        let app = ClapApp::App.get();
        let matches =
            app.get_matches_from(vec!["bader", "CHGCAR", "-r", "CHGCAR_sum"]);
        let args = Args::new(matches);
        let flag = matches!(args.reference, Reference::One(_));
        assert!(flag)
    }

    #[test]
    fn argument_reference_two() {
        let app = ClapApp::App.get();
        let v = vec!["bader", "CHGCAR", "-r", "AECCAR0", "--ref", "AECCAR2"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        let flag = matches!(args.reference, Reference::Two(_, _));
        assert!(flag)
    }

    #[test]
    fn argument_reference_none() {
        let app = ClapApp::App.get();
        let v = vec!["bader", "CHGCAR"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        let flag = matches!(args.reference, Reference::None);
        assert!(flag)
    }

    #[test]
    fn argument_aeccar() {
        let app = ClapApp::App.get();
        let v = vec!["bader", "CHGCAR", "-a"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        let flag = match args.reference {
            Reference::Two(x, y) => (x == *"AECCAR0") && (y == *"AECCAR2"),
            _ => false,
        };
        assert!(flag)
    }

    #[test]
    #[should_panic]
    fn argument_aeccar_cube() {
        let app = ClapApp::App.get();
        let v = vec!["bader", "charge.cube", "-a"];
        let matches = app.get_matches_from(v);
        let _ = Args::new(matches);
    }

    #[test]
    fn argument_vacuum_tolerance_auto() {
        let app = ClapApp::App.get();
        let v = vec!["bader", "CHGCAR", "-v", "auto"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        assert_eq!(args.vacuum_tolerance, Some(1E-3))
    }

    #[test]
    fn argument_vacuum_tolerance_float() {
        let app = ClapApp::App.get();
        let v = vec!["bader", "CHGCAR", "--vac", "1E-4"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        assert_eq!(args.vacuum_tolerance, Some(1E-4))
    }

    #[test]
    #[should_panic]
    fn argument_vacuum_tolerance_not_float() {
        let app = ClapApp::App.get();
        let v = vec!["bader", "CHGCAR", "-v", "0.00.1"];
        let matches = app.get_matches_from(v);
        let _ = Args::new(matches);
    }

    #[test]
    fn argument_threads_default() {
        let app = ClapApp::App.get();
        let v = vec!["bader", "CHGCAR"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        assert_eq!(args.threads, 0)
    }

    #[test]
    fn argument_threads_int() {
        let app = ClapApp::App.get();
        let v = vec!["bader", "CHGCAR", "--threads", "1"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        assert_eq!(args.threads, 1)
    }

    #[test]
    #[should_panic]
    fn argument_threads_not_int() {
        let app = ClapApp::App.get();
        let v = vec!["bader", "CHGCAR", "-J", "0.1"];
        let matches = app.get_matches_from(v);
        let _ = Args::new(matches);
    }
}
