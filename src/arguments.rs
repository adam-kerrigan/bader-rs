use crate::io::{FileType, WriteType};
use clap::{crate_authors, value_parser, Arg, ArgAction, ArgMatches, Command};

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

impl ClapApp {
    /// Create and return the clap::App
    pub fn get() -> Command {
        Command::new("Multi-threaded Bader Charge Analysis")
            .author(crate_authors!())
            .version("0.4.2")
            .arg(Arg::new("file")
                .required(true)
                .index(1)
                .help("The file to analyse."))
            .arg(Arg::new("output")
                .short('o')
                .long("output")
                .action(ArgAction::SetTrue)
                .help("Output the Bader atoms.")
                .long_help(
"Output the Bader atoms in the same file format as the input density.
This can be used in conjunction with the index flag to specify a
specific atom. Without the index flag it will print all the atoms."))
            .arg(Arg::new("index")
                .short('i')
                .long("index")
                .action(ArgAction::Append)
                .requires("output")
                .value_parser(1..)
                .help("Index of Bader atoms to be written out.")
                .long_help(
"An index of a Bader atom to be written out, starting at 1. This flag
requires the output flag to be set. Multiple atoms can be written by
repeating the flag ie. bca CHGCAR -oi 1 -i 2."))
            .arg(Arg::new("file type")
                .short('t')
                .long("type")
                .num_args(1)
                .value_parser(["cube", "vasp"])
                .ignore_case(false)
                .help("The file type of the charge density.")
                .long_help(
"The file type of the input file. If this is not supplied the type will attempt
to be infered from the filename."))
            .arg(Arg::new("reference")
                .short('r')
                .long("ref")
                .action(ArgAction::Append)
                .num_args(1)
                .value_parser(value_parser!(String))
                .help("File(s) containing reference charge.")
                .long_help(
"A reference charge to do the partitioning upon. Two files can be passed
by using multiple flags (bca CHGCAR -r AECCAR0 -r AECCAR2). If two files are
passed they are summed together."))
            .arg(Arg::new("spin")
                .short('s')
                .long("spin")
                .num_args(1)
                .value_parser(value_parser!(String))
                .help("File containing spin density.")
                .long_help(
"A path to the spin density associated with the original file. This is primarily
for cube files as if spin density exists in a CHGCAR it will be read automatically.
If using with VASP outputs then the files for charge and spin density must only
contain a single density (ie. the original file has been split)."))
            .arg(Arg::new("all electron")
                .short('a')
                .long("aec")
                .help("Convience flag for reading both aeccars.")
                .action(ArgAction::SetTrue)
                .conflicts_with("reference"))
            .arg(Arg::new("vacuum tolerance")
                .short('v')
                .long("vac")
                .default_value("1E-6")
                .value_parser(value_parser!(f64))
                .num_args(1)
                .help("Cut-off at which charge is considered vacuum.")
                .long_help(
"Values of density below the supplied value are considered vacuum and are not
included in the calculation."))
            .arg(Arg::new("weight tolerance")
                .short('w')
                .long("weight")
                .default_value("1E-8")
                .value_parser(value_parser!(f64))
                .num_args(1)
                .help("Cut-off at which contributions to the weighting will be ignored.")
                .long_help(
"Values of density below the supplied value are ignored from the weighting and
included in the calculation. By raising the tolerance the calculation speed can
be increased but every ignored weight is unaccounted for in the final partitions.
Be sure to test this!"))
            .arg(Arg::new("maximum distance")
                .short('m')
                .long("max_dist")
                .default_value("0.1")
                .value_parser(value_parser!(f64))
                .num_args(1)
                .help("Cut-off after which an error will be thrown for distance of Bader maximum to atom.")
                .long_help(
"The distance allowed between Bader maximum and its associated atom, in angstrom, before
an error is thrown. This will cause a hard crash of the program, consider whether
increasing the cut-off or adding a \"ghost atom\" at the location of the Bader
maximum is more appropriate."))
            .arg(Arg::new("threads")
                .short('J')
                .long("threads")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Number of threads to distribute the calculation over.")
                .long_help(
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
    pub maximum_distance: f64,
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
        let file = arguments.get_one::<String>("file").unwrap().to_string();
        // Collect write charge info
        let output = if arguments.get_flag("output") {
            let atoms: Vec<isize> = match arguments.get_many::<i64>("index") {
                Some(vec) => vec.into_iter().map(|i| *i as isize - 1).collect(),
                None => Vec::with_capacity(0),
            };
            WriteType::Atom(atoms)
        } else {
            WriteType::None
        };

        // Collect file type
        let file_type = arguments.get_one::<String>("file type");
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
        // safe to unwrap as threads has a default value of 1E-8
        let weight_tolerance = arguments.get_one::<f64>("weight tolerance")
                                        .unwrap()
                                        .max(1E-13);
        // Collect maxima tolerance
        // safe to unwrap as threads has a default value of 0.1
        let maximum_distance =
            *arguments.get_one::<f64>("maximum distance").unwrap();
        // Collect threads
        // safe to unwrap as threads has a default value of 0
        let threads = match *arguments.get_one::<usize>("threads").unwrap() {
            0 => num_cpus::get().min(12),
            x => x,
        };
        // Collect vacuum tolerance
        let vacuum_tolerance =
            arguments.get_one::<f64>("vacuum tolerance").copied();
        // Collect reference files
        let references: Vec<_> = if arguments.get_flag("all electron") {
            if let FileType::Vasp = file_type {
                vec!["AECCAR0", "AECCAR2"]
            } else {
                panic!("Error: Cannot use AECCAR flag for non VASP file-types.")
            }
        } else {
            match arguments.get_many::<String>("reference") {
                Some(x) => x.map(|s| s.as_str()).collect(),
                None => Vec::with_capacity(0),
            }
        };
        let reference = match references.len() {
            2 => Reference::Two(String::from(references[0]),
                                String::from(references[1])),
            1 => Reference::One(String::from(references[0])),
            _ => Reference::None,
        };
        let spin = arguments.get_one::<String>("spin")
                            .map(|x| String::from(x.as_str()));
        Self { file,
               file_type,
               weight_tolerance,
               maximum_distance,
               output,
               reference,
               spin,
               threads,
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
            app.get_matches_from(vec!["bca", "charge.cube", "--type", "cube"]);
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
        let matches = app.get_matches_from(vec!["bca", "CHGCAR", "-o"]);
        let args = Args::new(matches);
        match args.output {
            WriteType::Atom(v) => assert!(v.is_empty()),
            _ => panic!(),
        }
    }

    #[test]
    fn argument_output_index() {
        let app = ClapApp::get();
        let matches =
            app.get_matches_from(vec!["bca", "CHGCAR", "-o", "-i", "1",]);
        let args = Args::new(matches);
        match args.output {
            WriteType::Atom(v) => assert_eq!(v, vec![0]),
            _ => panic!(),
        }
    }

    #[test]
    fn argument_output_mult_index() {
        let app = ClapApp::get();
        let matches = app.get_matches_from(vec!["bca", "CHGCAR", "-o",
                                                "--index", "1", "-i", "3",]);
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
        let _ = app.try_get_matches_from(vec!["bca", "CHGCAR", "-o", "-i",
                                              "0"])
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
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
        let v = vec!["bca", "CHGCAR", "--vac", "0.00.1"];
        let _ = app.try_get_matches_from(v)
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
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
        let _ = app.try_get_matches_from(v)
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_maximum_distance_float() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "--max_dist", "1E-4"];
        let matches = app.get_matches_from(v);
        let args = Args::new(matches);
        assert_eq!(args.maximum_distance, 1E-4)
    }

    #[test]
    #[should_panic]
    fn argument_maximum_distance_not_float() {
        let app = ClapApp::get();
        let v = vec!["bca", "CHGCAR", "-m", "0.00.1"];
        let _ = app.try_get_matches_from(v)
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
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
        let _ = app.try_get_matches_from(v)
                   .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }
}
