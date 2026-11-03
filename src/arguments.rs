use crate::{
    errors::ArgumentError,
    io::{FileType, WriteType, cube::Cube, vasp::Vasp},
};
use std::collections::HashMap;
use std::fmt::{Debug, Display, Write};
use std::str::FromStr;
use std::thread::available_parallelism;

/// Indicates how many reference files are passed
#[derive(Clone, Debug, PartialEq)]
pub enum Reference {
    /// One file as a reference, usually either a CHGCAR_sum or spin.cube.
    One(String),
    /// Two files as a reference, these files will be summed together.
    Two(String, String),
    /// No reference, just use the density file.
    None,
}

/// The definition of a command-line option.
#[derive(Debug)]
struct Arg {
    /// Name of the argument, used as the key for its values.
    name: &'static str,
    /// Text displayed as help for the argument when -h is passed.
    short_help: &'static str,
    /// Text displayed as help for the argument when --help is passed.
    long_help: &'static str,
    /// Single character for flag e.g. -h.
    short_flag: char,
    /// Longer flag e.g. --help.
    long_flag: &'static str,
    /// Whether the argument takes a value.
    takes_value: bool,
    /// Whether the argument can be passed multiple times.
    multiple_values: bool,
    /// The default value, if any, as it is parsed and displayed.
    default_value: Option<&'static str>,
    /// Whether to restrict the argument to specific values and if so what they are.
    allowed_values: &'static [&'static str],
}

/// The main configuration parser for the application.
///
/// This struct defines the available command-line arguments and handles the logic
/// for parsing the raw argument vector into a structured [`Args`] object.
///
/// # Logic
/// It supports:
/// * **Short Flags**: Single character (e.g., `-h`).
/// * **Long Flags**: Full name (e.g., `--help`).
/// * **Clusters**: Combining short flags (e.g., `-xo` equivalent to `-x -o`).
/// * **Positional Arguments**: Identifying the input filename amidst flags.
///
/// # Examples
/// ```
/// use bader::arguments::App;
///
/// let app = App::new();
/// let args = vec!["bca", "CHGCAR", "-t", "4"];
/// let parsed = app.parse_args(args).expect("Failed to parse args");
///
/// assert_eq!(parsed.file, "CHGCAR");
/// assert_eq!(parsed.threads, 4);
/// ```
pub struct App {
    /// What arguments it can take, [OPTIONS] in the help.
    options: Vec<Arg>,
    /// The help specific argument.
    help: Arg,
    /// The longest argument name's size.
    max_width: usize,
}

impl App {
    /// Defines the options of the program.
    pub fn new() -> Self {
        // Start options
        let options = vec![
            Arg {
                name: "aec",
                short_help: "Convenience flag for reading both AECCARs.",
                long_help: "
\tFlag for reading and summing both the AECCAR0 and AECCAR2 from a VASP calculation",
                short_flag: 'a',
                long_flag: "aec",
                takes_value: false,
                multiple_values: false,
                default_value: None,
                allowed_values: &[],
            },
            Arg {
                name: "file type",
                short_help: "The file format of the density file.",
                long_help: "
\tThe file format of the density file. If not set the format is guessed from the
\tfilename: a name containing \"cube\" is read as a cube file and anything else
\tis read as a VASP file.",
                short_flag: 'f',
                long_flag: "file_type",
                takes_value: true,
                multiple_values: false,
                default_value: None,
                allowed_values: &["cube", "vasp"],
            },
            Arg {
                name: "index",
                short_help: "Index of Bader atoms to be written out.",
                long_help: "
\tAn index of a Bader atom to be written out, starting at 1. This flag requires
\tthe output flag to be set. Multiple atoms can only be written by repeating the
\tflag ie. bca CHGCAR -oi 1 -i 2.",
                short_flag: 'i',
                long_flag: "index",
                takes_value: true,
                multiple_values: true,
                default_value: None,
                allowed_values: &[],
            },
            Arg {
                name: "maximum distance",
                short_help: "Cut-off after which an error will be thrown for distance of Bader maximum to atom.",
                long_help: "
\tThe distance allowed between Bader maximum and its associated atom, in angstrom, before
\tan error is thrown. This will cause a hard crash of the program, consider whether
\tincreasing the cut-off or adding a \"ghost atom\" at the location of the Bader
\tmaximum is more appropriate.",
                short_flag: 'm',
                long_flag: "max_dist",
                takes_value: true,
                multiple_values: false,
                default_value: Some("1"),
                allowed_values: &[],
            },
            Arg {
                name: "output",
                short_help: "Output the Bader atoms.",
                long_help: "
\tOutput the bader atoms in the same file format as the input density.
\tthis can be used in conjunction with the index flag to specify a
\tspecific atom. without the index flag it will print all the atoms.",
                short_flag: 'o',
                long_flag: "output",
                takes_value: false,
                multiple_values: false,
                default_value: None,
                allowed_values: &[],
            },
            Arg {
                name: "reference",
                short_help: "File(s) containing reference charge.",
                long_help: "
\tA reference charge to do the partitioning upon. Two files can be passed
\tby using multiple flags (bca CHGCAR -r AECCAR0 -r AECCAR2). If two files are
\tpassed they are summed together.",
                short_flag: 'r',
                long_flag: "ref",
                takes_value: true,
                multiple_values: true,
                default_value: None,
                allowed_values: &[],
            },
            Arg {
                name: "silent",
                short_help: "Whether to display any output.",
                long_help: "
\tRuns the program without displaying any output text or progress bars.",
                short_flag: 'x',
                long_flag: "silent",
                takes_value: false,
                multiple_values: false,
                default_value: None,
                allowed_values: &[],
            },
            Arg {
                name: "spin",
                short_help: "File containing spin density.",
                long_help: "
\tA path to the spin density associated with the original file. This is primarily
\tfor cube files as if spin density exists in a CHGCAR it will be read automatically.
\tIf using with VASP outputs then the files for charge and spin density must only
\tcontain a single density (ie. the original file has been split).",
                short_flag: 's',
                long_flag: "spin",
                takes_value: true,
                multiple_values: false,
                default_value: None,
                allowed_values: &[],
            },
            Arg {
                name: "threads",
                short_help: "Number of threads to distribute the calculation over.",
                long_help: "
\tThe number of threads to be used by the program. A default value of 0 is used
\tto allow the program to best decide how to use the available hardware. It does
\tthis by using the minimum value out of the number cores available and 12.",
                short_flag: 't',
                long_flag: "threads",
                takes_value: true,
                multiple_values: false,
                default_value: Some("0"),
                allowed_values: &[],
            },
            Arg {
                name: "vacuum tolerance",
                short_help: "Cut-off at which charge is considered vacuum.",
                long_help: "
\tValues of density below the supplied value are considered vacuum and are not
\tincluded in the calculation.",
                short_flag: 'v',
                long_flag: "vac",
                takes_value: true,
                multiple_values: false,
                default_value: Some("0.000001"),
                allowed_values: &[],
            },
            Arg {
                name: "weight tolerance",
                short_help: "Cut-off at which contributions to the weighting will be ignored.",
                long_help: "
\tWeights below the supplied value are ignored from the weighting and not
\tincluded in the calculation. By raising the tolerance the calculation speed can
\tbe increased but every ignored weight is unaccounted for in the final partitions.
\tBe sure to test this!",
                short_flag: 'w',
                long_flag: "weight",
                takes_value: true,
                multiple_values: false,
                default_value: Some("0.000001"),
                allowed_values: &[],
            },
        ];
        // End options
        let help = Arg {
            name: "help",
            short_help: "Print help (see more with '--help')",
            long_help: "Print help (see a summary with '-h')",
            short_flag: 'h',
            long_flag: "help",
            takes_value: false,
            multiple_values: false,
            default_value: None,
            allowed_values: &[],
        };
        let max_width =
            options.iter().map(|o| o.flags().len()).max().unwrap_or(0);
        Self {
            options,
            help,
            max_width,
        }
    }

    /// Get the text for all the options either in short or long format.
    fn get_options_text(&self, long: bool) -> String {
        self.options.iter().fold(String::new(), |mut output, o| {
            let help = if long {
                o.long_help.to_string()
            } else {
                format!("  {}", o.short_help)
            };
            let default = match (o.default_value, long) {
                (None, true) => String::from("\n"),
                (None, false) => String::new(),
                (Some(d), true) => format!("\n\n\t[default: {}]\n", d),
                (Some(d), false) => format!(" [default: {}]", d),
            };
            let allowed = match (o.allowed_values.is_empty(), long) {
                (true, _) => String::new(),
                (false, true) => {
                    format!("\n\t[possible values: {:?}]\n", o.allowed_values)
                }
                (false, false) => {
                    format!(" [possible values: {:?}]", o.allowed_values)
                }
            };
            let _ = writeln!(
                output,
                " {:<width$}{}{}{}",
                o.flags(),
                help,
                default,
                allowed,
                width = self.max_width
            );
            output
        })
    }

    /// Get the help text in long or short formats.
    fn get_help_text(&self, long: bool) -> String {
        format!(
            " {:<width$}{}",
            self.help.flags(),
            if long {
                format!("\n\t{}", self.help.long_help)
            } else {
                format!("  {}", self.help.short_help)
            },
            width = self.max_width
        )
    }

    /// Get the `Arg` from a supplied short flag.
    fn get_option_from_short_flag(
        &self,
        f: char,
    ) -> Result<&Arg, ArgumentError<'_>> {
        if f == self.help.short_flag {
            return Err(ArgumentError::ShortHelp(self));
        }
        self.options
            .iter()
            .find(|opt| opt.short_flag == f)
            .ok_or_else(|| ArgumentError::NotFlag(f.to_string()))
    }

    /// Get the `Arg` from a supplied long flag.
    fn get_option_from_long_flag(
        &self,
        f: &str,
    ) -> Result<&Arg, ArgumentError<'_>> {
        if f == self.help.long_flag {
            return Err(ArgumentError::LongHelp(self));
        }
        self.options
            .iter()
            .find(|opt| opt.long_flag == f)
            .ok_or_else(|| ArgumentError::NotFlag(f.to_string()))
    }

    /// Get the `Arg` with the given name.
    fn get_option(&self, name: &str) -> &Arg {
        self.options.iter().find(|opt| opt.name == name).unwrap()
    }

    /// Parse the single value of the option `name`, falling back to the
    /// default if it was not passed.
    fn parse_value<T: FromStr>(
        &self,
        values: &HashMap<&str, Vec<String>>,
        name: &'static str,
    ) -> Result<T, ArgumentError<'_>> {
        let option = self.get_option(name);
        let value = match values.get(name) {
            Some(v) => v.last().unwrap().as_str(),
            None => option.default_value.unwrap(),
        };
        value.parse().map_err(|_| {
            ArgumentError::Unparsable(
                name.to_string(),
                value.to_string(),
                std::any::type_name::<T>().to_string(),
            )
        })
    }

    /// Parse arguments from flags and values.
    ///
    /// The first argument is the name of the program and is ignored. Every
    /// other argument that does not start with `-` is the file to analyse.
    pub fn parse_args(
        &self,
        args: Vec<&str>,
    ) -> Result<Args, ArgumentError<'_>> {
        // the values passed for each option, keyed by the option's name
        let mut values: HashMap<&str, Vec<String>> = HashMap::new();
        let mut file = None;
        let mut args = args.into_iter().skip(1);
        // stores the option, taking the next argument as its value if needed
        let mut store = |option: &Arg,
                         flag: &str,
                         args: &mut dyn Iterator<Item = &str>|
         -> Result<(), ArgumentError<'_>> {
            let value = if option.takes_value {
                args.next()
                    .ok_or_else(|| ArgumentError::NoValue(flag.to_string()))?
            } else {
                ""
            };
            let entry = values.entry(option.name).or_default();
            if !option.multiple_values {
                entry.clear();
            }
            entry.push(value.to_string());
            Ok(())
        };
        while let Some(arg) = args.next() {
            if let Some(long_flag) = arg.strip_prefix("--") {
                let option = self.get_option_from_long_flag(long_flag)?;
                store(option, long_flag, &mut args)?;
            } else if let Some(short_flags) = arg.strip_prefix('-') {
                if short_flags.is_empty() {
                    return Err(ArgumentError::NotFlag(arg.to_string()));
                }
                // a cluster of short flags, each may take the next argument
                for flag in short_flags.chars() {
                    let option = self.get_option_from_short_flag(flag)?;
                    store(option, &flag.to_string(), &mut args)?;
                }
            } else {
                file = Some(arg.to_string());
            }
        }
        let file = file.ok_or(ArgumentError::NoFile(self))?;
        let file_type = match values.get("file type") {
            Some(ftype) => {
                match ftype.last().unwrap().to_lowercase().as_str() {
                    "cube" => FileType::Cube(Cube {}),
                    "vasp" => FileType::Vasp(Vasp {}),
                    other => {
                        return Err(ArgumentError::NotValidValue(
                            String::from("file type"),
                            other.to_string(),
                        ));
                    }
                }
            }
            None => parse_filetype(&file),
        };
        let weight_tolerance = self.parse_value(&values, "weight tolerance")?;
        let maximum_distance = self.parse_value(&values, "maximum distance")?;
        let vacuum_tolerance = self.parse_value(&values, "vacuum tolerance")?;
        let threads = match self.parse_value::<usize>(&values, "threads")? {
            0 => available_parallelism().map_or(1, |n| n.get().min(12)),
            t => t,
        };
        let output = match (values.contains_key("output"), values.get("index"))
        {
            (true, Some(indices)) => WriteType::Atom(
                indices
                    .iter()
                    .map(|i| match i.parse::<usize>() {
                        Ok(i) if i > 0 => Ok(i - 1),
                        _ => Err(ArgumentError::Unparsable(
                            String::from("index"),
                            i.to_string(),
                            String::from("positive integer"),
                        )),
                    })
                    .collect::<Result<Vec<usize>, _>>()?,
            ),
            (true, None) => WriteType::Atom(Vec::new()),
            (false, Some(_)) => {
                return Err(ArgumentError::MissingDependant(
                    String::from("index"),
                    String::from("output"),
                ));
            }
            (false, None) => WriteType::None,
        };
        let reference = match values.get("reference").map(|v| v.as_slice()) {
            Some([f]) => Reference::One(f.clone()),
            Some([f1, f2]) => Reference::Two(f1.clone(), f2.clone()),
            Some(v) => {
                return Err(ArgumentError::TooManyValues(
                    String::from("reference"),
                    2,
                    v.len(),
                ));
            }
            None if values.contains_key("aec") => match file_type {
                FileType::Vasp(_) => Reference::Two(
                    String::from("AECCAR0"),
                    String::from("AECCAR2"),
                ),
                FileType::Cube(_) => {
                    return Err(ArgumentError::WrongFileType(
                        String::from("aec"),
                        String::from("cube"),
                    ));
                }
            },
            None => Reference::None,
        };
        Ok(Args {
            file,
            file_type,
            weight_tolerance,
            maximum_distance,
            output,
            reference,
            silent: values.contains_key("silent"),
            spin: values.get("spin").map(|s| s.last().unwrap().clone()),
            threads,
            vacuum_tolerance,
        })
    }
}

impl Arg {
    /// The flags as displayed in the help, e.g. `-t,--threads [threads]`.
    fn flags(&self) -> String {
        if self.takes_value {
            format!("-{},--{} [{}]", self.short_flag, self.long_flag, self.name)
        } else {
            format!("-{},--{}", self.short_flag, self.long_flag)
        }
    }
}

impl Default for App {
    fn default() -> Self {
        Self::new()
    }
}

/// The short help, printed with `-h`.
impl Display for App {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Usage: bca [OPTIONS] <file>\n\nArguments:\n <file>  The file to analyse.\n\nOptions:\n{}{}",
            self.get_options_text(false),
            self.get_help_text(false)
        )
    }
}

/// The long help, printed with `--help`.
impl Debug for App {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Usage: bca [OPTIONS] <file>\n\nArguments:\n <file>  The file containing the charge density to analyse.\n\nOptions:\n{}{}",
            self.get_options_text(true),
            self.get_help_text(true)
        )
    }
}

/// The validated configuration for a Bader analysis run.
///
/// This struct is produced by [`App::parse_args`] and contains all the settings required
/// to execute the analysis, including file paths, thread counts, and physical tolerances.
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
    /// Should the program be ran silently.
    pub silent: bool,
    /// Is there a spin density to include as well.
    pub spin: Option<String>,
    /// How many threads to use in the calculation.
    pub threads: usize,
    /// Is there a tolerance to consider a density vacuum.
    pub vacuum_tolerance: f64,
}

/// Parse the file type from the filename.
pub fn parse_filetype(fname: &str) -> FileType {
    let f = fname.to_lowercase();
    if f.contains("cube") {
        FileType::Cube(Cube {})
    } else {
        if !f.contains("car") {
            eprintln!("Cannot detect file type, attempting to read as VASP.");
        }
        FileType::Vasp(Vasp {})
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn argument_file() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR"];
        let args = app.parse_args(v).unwrap();
        assert_eq!(args.file, String::from("CHGCAR"));
    }

    #[test]
    #[should_panic]
    fn argument_no_file() {
        let app = App::new();
        let _ = app
            .parse_args(vec!["bca"])
            .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_file_type_default_vasp() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.file_type, FileType::Vasp(_));
        assert!(flag);
    }

    #[test]
    fn argument_file_type_default_unknown() {
        let app = App::new();
        let v = vec!["bca", "CHG"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.file_type, FileType::Vasp(_));
        assert!(flag);
    }

    #[test]
    fn argument_file_type_vasp() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "-f", "vasp"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.file_type, FileType::Vasp(_));
        assert!(flag);
    }

    #[test]
    fn argument_file_type_default_cube() {
        let app = App::new();
        let v = vec!["bca", "charge.cube"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.file_type, FileType::Cube(_));
        assert!(flag);
    }

    #[test]
    fn argument_file_type_cube() {
        let app = App::new();
        let v = vec!["bca", "charge.cube", "--file_type", "cube"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.file_type, FileType::Cube(_));
        assert!(flag);
    }

    #[test]
    #[should_panic]
    fn argument_file_type_not_type() {
        let app = App::new();
        let _ = app
            .parse_args(vec!["bca", "CHGCAR", "-f", "basp"])
            .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_output_atoms() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "-o"];
        let args = app.parse_args(v).unwrap();
        match args.output {
            WriteType::Atom(v) => assert!(v.is_empty()),
            _ => panic!(),
        }
    }

    #[test]
    fn argument_output_index() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "-o", "-i", "1"];
        let args = app.parse_args(v).unwrap();
        match args.output {
            WriteType::Atom(v) => assert_eq!(v, vec![0]),
            _ => panic!(),
        }
    }

    #[test]
    fn argument_output_mult_index() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "-o", "--index", "1", "-i", "3"];
        let args = app.parse_args(v).unwrap();
        match args.output {
            WriteType::Atom(v) => assert_eq!(v, vec![0, 2]),
            _ => panic!(),
        }
    }

    #[test]
    #[should_panic]
    fn argument_index_no_output() {
        let app = App::new();
        let _ = app
            .parse_args(vec!["bca", "CHGCAR", "-i", "1"])
            .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    #[should_panic]
    fn argument_index_not_parse() {
        let app = App::new();
        let _ = app
            .parse_args(vec!["bca", "CHGCAR", "-oi", "1,6"])
            .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_spin() {
        let app = App::new();
        let v = vec!["bca", "density.cube", "-s", "spin.cube"];
        let args = app.parse_args(v).unwrap();
        assert_eq!(args.spin, Some(String::from("spin.cube")))
    }

    #[test]
    fn argument_reference_one() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "-r", "CHGCAR_sum"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.reference, Reference::One(_));
        assert!(flag)
    }

    #[test]
    fn argument_reference_two() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "-r", "AECCAR0", "--ref", "AECCAR2"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.reference, Reference::Two(_, _));
        assert!(flag)
    }

    #[test]
    fn argument_reference_none() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.reference, Reference::None);
        assert!(flag)
    }

    #[test]
    fn argument_aeccar() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "-a"];
        let args = app.parse_args(v).unwrap();
        let flag = match args.reference {
            Reference::Two(x, y) => (x == *"AECCAR0") && (y == *"AECCAR2"),
            _ => false,
        };
        assert!(flag)
    }

    #[test]
    #[should_panic]
    fn argument_aeccar_cube() {
        let app = App::new();
        let v = vec!["bca", "charge.cube", "-a"];
        let _ = app.parse_args(v).unwrap();
    }

    #[test]
    fn argument_vacuum_tolerance_float() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "--vac", "1E-4"];
        let args = app.parse_args(v).unwrap();
        assert_eq!(args.vacuum_tolerance, 1E-4)
    }

    #[test]
    fn argument_vacuum_tolerance_default() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR"];
        let args = app.parse_args(v).unwrap();
        assert_eq!(args.vacuum_tolerance, 1E-6)
    }

    #[test]
    #[should_panic]
    fn argument_vacuum_tolerance_not_float() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "--vac", "0.00.1"];
        let _ = app
            .parse_args(v)
            .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_weight_tolerance_float() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "--weight", "1E-4"];
        let args = app.parse_args(v).unwrap();
        assert_eq!(args.weight_tolerance, 1E-4)
    }

    #[test]
    #[should_panic]
    fn argument_weight_tolerance_not_float() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "-w", "0.00.1"];
        let _ = app
            .parse_args(v)
            .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_maximum_distance_float() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "--max_dist", "1E-4"];
        let args = app.parse_args(v).unwrap();
        assert_eq!(args.maximum_distance, 1E-4)
    }

    #[test]
    #[should_panic]
    fn argument_maximum_distance_not_float() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "-m", "0.00.1"];
        let _ = app
            .parse_args(v)
            .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn argument_threads_default() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR"];
        let args = app.parse_args(v).unwrap();
        let threads = available_parallelism().unwrap().get().min(12);
        assert_eq!(args.threads, threads)
    }

    #[test]
    fn argument_threads_int() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "--threads", "1"];
        let args = app.parse_args(v).unwrap();
        assert_eq!(args.threads, 1)
    }

    #[test]
    #[should_panic]
    fn argument_threads_not_int() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "-t", "0.1"];
        let _ = app
            .parse_args(v)
            .unwrap_or_else(|e| panic!("An error occurs: {}", e));
    }

    #[test]
    fn test_app_creation() {
        let app = App::new();
        assert!(!app.options.is_empty());
    }

    #[test]
    fn test_positional_file_parsing() {
        let app = App::new();
        // File at end
        let v1 = vec!["bca", "-t", "4", "CHGCAR"];
        let args1 = app.parse_args(v1).unwrap();
        assert_eq!(args1.file, "CHGCAR");
        assert_eq!(args1.threads, 4);

        // File in middle
        let v2 = vec!["bca", "-t", "2", "CHGCAR", "-x"];
        let args2 = app.parse_args(v2).unwrap();
        assert_eq!(args2.file, "CHGCAR");
        assert_eq!(args2.threads, 2);

        // File at start
        let v2 = vec!["bca", "CHGCAR", "-t", "2"];
        let args2 = app.parse_args(v2).unwrap();
        assert_eq!(args2.file, "CHGCAR");
        assert_eq!(args2.threads, 2);
    }

    #[test]
    fn test_help_flags() {
        let app = App::new();
        assert!(matches!(
            app.parse_args(vec!["bca", "-h"]),
            Err(ArgumentError::ShortHelp(_))
        ));
        assert!(matches!(
            app.parse_args(vec!["bca", "CHGCAR", "--help"]),
            Err(ArgumentError::LongHelp(_))
        ));
        assert!(matches!(
            app.parse_args(vec!["bca", "-t", "2"]),
            Err(ArgumentError::NoFile(_))
        ));
        // the help lists every option
        let short = format!("{}", app);
        let long = format!("{:?}", app);
        for option in &app.options {
            assert!(short.contains(&format!("--{}", option.long_flag)));
            assert!(long.contains(option.long_help));
        }
    }

    #[test]
    fn test_argument_errors() {
        let app = App::new();
        assert!(matches!(
            app.parse_args(vec!["bca", "CHGCAR", "-t"]),
            Err(ArgumentError::NoValue(_))
        ));
        assert!(matches!(
            app.parse_args(vec!["bca", "CHGCAR", "--nope"]),
            Err(ArgumentError::NotFlag(_))
        ));
        assert!(matches!(
            app.parse_args(vec!["bca", "CHGCAR", "-q"]),
            Err(ArgumentError::NotFlag(_))
        ));
        assert!(matches!(
            app.parse_args(vec![
                "bca", "CHGCAR", "-r", "a", "-r", "b", "-r", "c"
            ]),
            Err(ArgumentError::TooManyValues(_, 2, 3))
        ));
        assert!(matches!(
            app.parse_args(vec!["bca", "CHGCAR", "-oi", "0"]),
            Err(ArgumentError::Unparsable(_, _, _))
        ));
        // repeated single value options take the last value
        let args = app
            .parse_args(vec!["bca", "CHGCAR", "-t", "1", "-t", "3"])
            .unwrap();
        assert_eq!(args.threads, 3);
    }

    #[test]
    fn test_short_flag_cluster() {
        let app = App::new();
        // -x (silent), -a (aec) are booleans.
        // We can combine them: -xa
        let v = vec!["bca", "CHGCAR", "-xa"];

        let args = app.parse_args(v).unwrap();
        match args.reference {
            Reference::Two(a, b) => {
                assert_eq!(a, "AECCAR0");
                assert_eq!(b, "AECCAR2");
            }
            _ => panic!("Expected AECCAR reference from -a flag"),
        }
    }
}
