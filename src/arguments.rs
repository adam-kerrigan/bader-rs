use crate::io::{FileType, WriteType};
use clap::{crate_authors, App, Arg, ArgMatches};

/// Indicates how many reference files are passed
#[derive(Clone)]
pub enum Reference {
    /// One file as a reference, usually either a CHGCAR_sum or spin.cube.
    One(String),
    /// Two files as a reference, these files will be summed together.
    Two(String, String),
    /// No reference, just use the density file.
    None,
}

/// Create a container for dealing with clap and being able to test arg parsing.
pub enum ClapApp {}

impl<'a> ClapApp {
    /// Create and return the clap::App
    pub fn get() -> App<'a> {
        App::new("Multi-threaded Bader Charge Analysis")
            .author(crate_authors!())
            .version("0.3.1")
            .arg(Arg::new("file")
                .required(true)
                .index(1)
                .about("The file to analyse."))
            .arg(Arg::new("output")
                .short('o')
                .long("output")
                .takes_value(true)
                .possible_value("atoms")
                .possible_value("volumes")
                .case_insensitive(false)
                .about("Output the Bader atoms or volumes.")
                .long_about(
"Output the Bader atoms or the Bader volumes in the same file formtat as the
input density. This can be used in conjunction with the index flag to specify a
specific atoms or volumes. Without the index flag it will print all the atoms or
volumes."))
            .arg(Arg::new("index")
                .short('i')
                .long("index")
                .multiple(true)
                .number_of_values(1)
                .requires("output")
                .about("Index of Bader atoms or volumes to be written out.")
                .long_about(
"An index of a Bader atom or volume to be written out, starting at 1. This flag
requires the output flag to be set. Multiple atoms or volumes can be written by
repeating the flag ie. bca CHGCAR -o atoms -i 1 -i 2."))
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
by using multiple flags (bca CHGCAR -r AECCAR0 -r AECCAR2). If two files are
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
included in the calculation. A value of \"auto\" can be passed to use 1E-6 C*m^-3."))
            .arg(Arg::new("maxima tolerance")
                .short('m')
                .long("maxima")
                .takes_value(true)
                .about("Cut-off for charge at which a maxima is not printed.")
                .long_about(
"Values of charge for the Bader maxima below the supplied value are not written
to the Bader charge file (BCF.dat). A default value of 1E-6 is used."))
            .arg(Arg::new("weight tolerance")
                .short('w')
                .long("weight")
                .takes_value(true)
                .about("Cut-off at which contributions to the weighting will be ignored.")
                .long_about(
"Values of density below the supplied value are ignored from the weighting and
included in the calculation. A default vaule of 1E-6 is used. By raising the
tolerance the calculation speed can be increased but every ignored weight is
unaccounted for in the final partitions. Be sure to test this!"))
            .arg(Arg::new("threads")
                .short('J')
                .long("threads")
                .takes_value(true)
                .default_value("0")
                .about("Number of threads to distribute the calculation over.")
                .long_about(
"The number of threads to be used by the program. A default value of 0 is used
to allow the program to best decide how to use the available hardware. It does
this by using the minimum value out of the number cores available and 12."))
    }
}

/// Holds the arguments passed to the program from the command-line
pub struct Args {
    /// The filename.
    pub file: String,
    /// The file format.
    pub file_type: FileType,
    /// Tolerance to disregard weights at.
    pub weight_tolerance: f64,
    /// Tolerance to disregard maxima at.
    pub maxima_tolerance: f64,
    /// Output Writing
    pub output: WriteType,
    /// Is there a reference file.
    pub reference: Reference,
    /// Is there a spin density to include as well.
    pub spin: Option<String>,
    /// How many threads to use in the calculation.
    pub threads: usize,
    /// Is there a tolerance to consider a density vacuum.
    pub vacuum_tolerance: Option<f64>,
}

impl Args {
    /// Initialises the structure from the command-line arguments.
    pub fn new(arguments: ArgMatches) -> Self {
        // Collect file
        let file = match arguments.value_of("file") {
            Some(f) => String::from(f),
            None => String::new(),
        };

        // Collect write charge info
        let output = match arguments.value_of("output") {
            Some("atoms") => {
                let atoms = match arguments.values_of("index") {
                    Some(vec) => {
                        vec.map(|s| match s.parse::<usize>() {
                            Ok(u) => match u.checked_sub(1) {
                                Some(u) => u,
                                None => {
                                    panic!("Counting for index starts at 1.")
                                }
                            },
                            Err(_) => {
                                panic!("Unable to parse index, ({}) to usize.",
                                       s)
                            }
                        })
                        .collect::<Vec<usize>>()
                    }
                    None => Vec::with_capacity(0),
                };
                WriteType::Atom(atoms)
            }
            Some("volumes") => {
                let volumes = match arguments.values_of("index") {
                    Some(vec) => {
                        vec.map(|s| match s.parse::<usize>() {
                            Ok(u) => match u.checked_sub(1) {
                                Some(u) => u,
                                None => {
                                    panic!("Counting for index starts at 1.")
                                }
                            },
                            Err(_) => {
                                panic!("Unable to parse index, ({}) to usize.",
                                       s)
                            }
                        })
                        .collect::<Vec<usize>>()
                    }
                    None => Vec::with_capacity(0),
                };
                WriteType::Volume(volumes)
            }
            _ => WriteType::None,
        };

        // Collect file type
        let file_type = match arguments.value_of("file type") {
            Some(f) => Some(String::from(f)),
            None => None,
        };
        let file_type = match file_type {
            Some(ftype) => {
                if ftype.eq("cube") {
                    FileType::Cube
                } else {
                    FileType::Vasp
                }
            }
            None => {
                if file.to_lowercase().contains("cube") {
                    FileType::Cube
                } else if file.to_lowercase().contains("car") {
                    FileType::Vasp
                } else {
                    println!("Error: File-type cannot be infered, attempting to read as VASP");
                    FileType::Vasp
                }
            }
        };
        // Collect weight tolerance
        let weight_tolerance = match arguments.value_of("weight tolerance") {
            Some(x) => match x.parse::<f64>() {
                Ok(x) => x,
                Err(e) => {
                    panic!("Couldn't parse weight tolerance into float:\n{}", e)
                }
            },
            _ => 1E-6,
        };
        // Collect maxima tolerance
        let maxima_tolerance = match arguments.value_of("maxima tolerance") {
            Some(x) => match x.parse::<f64>() {
                Ok(x) => x,
                Err(e) => {
                    panic!("Couldn't parse maxima tolerance into float:\n{}", e)
                }
            },
            _ => 1E-6,
        };
        // Collect threads
        // safe to unwrap as threads has a default value of 0
        let threads = {
            match arguments.value_of("threads").unwrap().parse::<usize>() {
                Ok(0) => num_cpus::get().min(12),
                Ok(x) => x,
                Err(e) => panic!("Couldn't parse threads into integer:\n{}", e),
            }
        };
        // Collect vacuum tolerance
        let vacuum_tolerance = match arguments.value_of("vacuum tolerance") {
            Some(s) => {
                if s.eq("auto") {
                    Some(1E-6)
                } else {
                    match s.parse::<f64>() {
                        Ok(x) => Some(x),
                        Err(e) => panic!(
                            "Couldn't parse vacuum tolerance into float:\n{}",
                            e
                        ),
                    }
                }
            }
            None => None,
        };
        // Collect reference files
        let references: Vec<_> = if arguments.is_present("all electron") {
            if let FileType::Vasp = file_type {
                vec!["AECCAR0", "AECCAR2"]
            } else {
                panic!("Error: Cannot use AECCAR flag for non VASP file-types.")
            }
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
               file_type,
               weight_tolerance,
               maxima_tolerance,
               output,
               reference,
               threads,
               spin,
               vacuum_tolerance }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clapapp_get() {
        let app = ClapApp::get();
        assert_eq!(app.get_name(), "Multi-threaded Bader Charge Analysis")
    }

    #[test]
    fn argument_file() {
        let app = ClapApp::get();
        let matches = app.get_matches_from(vec!["bca", "CHGCAR"]);
        let args = Args::new(matches);
        assert_eq!(args.file, String::from("CHGCAR"));
    }

    #[test]
    #[should_panic]
    fn argument_no_file() {
        let app = ClapApp::get();
        let _ = app.try_get_matches_from(vec!["bca"])
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_file_type_default_vasp() {
        let app = ClapApp::get();
        let matches = app.get_matches_from(vec!["bca", "CHGCAR"]);
        let args = Args::new(matches);
        let flag = matches!(args.file_type, FileType::Vasp);
        assert!(flag);
    }

    #[test]
    fn argument_file_type_default_unknown() {
        let app = ClapApp::get();
        let matches = app.get_matches_from(vec!["bca", "CHG"]);
        let args = Args::new(matches);
        let flag = matches!(args.file_type, FileType::Vasp);
        assert!(flag);
    }

    #[test]
    fn argument_file_type_vasp() {
        let app = ClapApp::get();
        let matches = app.get_matches_from(vec!["bca", "CHGCAR", "-t", "vasp"]);
        let args = Args::new(matches);
        let flag = matches!(args.file_type, FileType::Vasp);
        assert!(flag);
    }

    #[test]
    fn argument_file_type_default_cube() {
        let app = ClapApp::get();
        let matches = app.get_matches_from(vec!["bca", "charge.cube"]);
        let args = Args::new(matches);
        let flag = matches!(args.file_type, FileType::Cube);
        assert!(flag);
    }

    #[test]
    fn argument_file_type_cube() {
        let app = ClapApp::get();
        let matches =
            app.get_matches_from(vec!["bca", "charge.cube", "--type", "cube",]);
        let args = Args::new(matches);
        let flag = matches!(args.file_type, FileType::Cube);
        assert!(flag);
    }

    #[test]
    #[should_panic]
    fn argument_file_type_not_type() {
        let app = ClapApp::get();
        let _ = app.try_get_matches_from(vec!["bca", "CHGCAR", "-t", "basp"])
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_output_atoms() {
        let app = ClapApp::get();
        let matches =
            app.get_matches_from(vec!["bca", "CHGCAR", "-o", "atoms"]);
        let args = Args::new(matches);
        match args.output {
            WriteType::Atom(v) => assert!(v.is_empty()),
            _ => panic!(),
        }
    }

    #[test]
    fn argument_output_volumes() {
        let app = ClapApp::get();
        let matches =
            app.get_matches_from(vec!["bca", "CHGCAR", "-o", "volumes"]);
        let args = Args::new(matches);
        match args.output {
            WriteType::Volume(v) => assert!(v.is_empty()),
            _ => panic!(),
        }
    }

    #[test]
    #[should_panic]
    fn argument_output_not_output() {
        let app = ClapApp::get();
        let _ = app.try_get_matches_from(vec!["bca", "CHGCAR", "-o", "Atoms"])
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_output_index() {
        let app = ClapApp::get();
        let matches = app.get_matches_from(vec!["bca", "CHGCAR", "-o",
                                                "volumes", "-i", "1"]);
        let args = Args::new(matches);
        match args.output {
            WriteType::Volume(v) => assert_eq!(v, vec![0]),
            _ => panic!(),
        }
    }

    #[test]
    fn argument_output_mult_index() {
        let app = ClapApp::get();
        let matches = app.get_matches_from(vec!["bca", "CHGCAR", "-o",
                                                "atoms", "--index", "1",
                                                "-i", "3"]);
        let args = Args::new(matches);
        match args.output {
            WriteType::Atom(v) => assert_eq!(v, vec![0, 2]),
            _ => panic!(),
        }
    }

    #[test]
    #[should_panic]
    fn argument_index_zero() {
        let app = ClapApp::get();
        let matches = app.get_matches_from(vec!["bca", "CHGCAR", "-o",
                                                "atoms", "-i", "0"]);
        let _ = Args::new(matches);
    }

    #[test]
    #[should_panic]
    fn argument_index_no_output() {
        let app = ClapApp::get();
        let _ = app.try_get_matches_from(vec!["bca", "CHGCAR", "-i", "1"])
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    #[should_panic]
    fn argument_index_not_parse() {
        let app = ClapApp::get();
        let _ = app.try_get_matches_from(vec!["bca", "CHGCAR", "-i", "1,6"])
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_spin() {
        let app = ClapApp::get();
        let matches = app.get_matches_from(vec!["bca",
                                                "density.cube",
                                                "-s",
                                                "spin.cube",]);
        let args = Args::new(matches);
        assert_eq!(args.spin, Some(String::from("spin.cube")))
    }

    #[test]
    fn argument_reference_one() {
        let app = ClapApp::get();
        let matches =
            app.get_matches_from(vec!["bca", "CHGCAR", "-r", "CHGCAR_sum"]);
        let args = Args::new(matches);
        let flag = matches!(args.reference, Reference::One(_));
        assert!(flag)
    }

    #[test]
    fn argument_reference_two() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "-r", "AECCAR0", "--ref", "AECCAR2"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        let flag = matches!(args.reference, Reference::Two(_, _));
        assert!(flag)
    }

    #[test]
    fn argument_reference_none() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        let flag = matches!(args.reference, Reference::None);
        assert!(flag)
    }

    #[test]
    fn argument_aeccar() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "-a"];
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
        let app = ClapApp::get();
        let v = vec!["bca", "charge.cube", "-a"];
        let matches = app.get_matches_from(v);
        let _ = Args::new(matches);
    }

    #[test]
    fn argument_vacuum_tolerance_auto() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "-v", "auto"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        assert_eq!(args.vacuum_tolerance, Some(1E-6))
    }

    #[test]
    fn argument_vacuum_tolerance_float() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "--vac", "1E-4"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        assert_eq!(args.vacuum_tolerance, Some(1E-4))
    }

    #[test]
    #[should_panic]
    fn argument_vacuum_tolerance_not_float() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "-v", "0.00.1"];
        let matches = app.get_matches_from(v);
        let _ = Args::new(matches);
    }

    #[test]
    fn argument_weight_tolerance_float() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "--weight", "1E-4"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        assert_eq!(args.weight_tolerance, 1E-4)
    }

    #[test]
    #[should_panic]
    fn argument_weight_tolerance_not_float() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "-w", "0.00.1"];
        let matches = app.get_matches_from(v);
        let _ = Args::new(matches);
    }

    #[test]
    fn argument_maxima_tolerance_float() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "--maxima", "1E-4"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        assert_eq!(args.maxima_tolerance, 1E-4)
    }

    #[test]
    #[should_panic]
    fn argument_maxima_tolerance_not_float() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "-m", "0.00.1"];
        let matches = app.get_matches_from(v);
        let _ = Args::new(matches);
    }

    #[test]
    fn argument_threads_default() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        let threads = num_cpus::get().min(12);
        assert_eq!(args.threads, threads)
    }

    #[test]
    fn argument_threads_int() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "--threads", "1"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        assert_eq!(args.threads, 1)
    }

    #[test]
    #[should_panic]
    fn argument_threads_not_int() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "-J", "0.1"];
        let matches = app.get_matches_from(v);
        let _ = Args::new(matches);
    }
}
