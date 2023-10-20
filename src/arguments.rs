use crate::{
    errors::ArgumentError,
    io::{FileType, WriteType},
};
use rustc_hash::FxHashMap;
use std::fmt::{Debug, Display, Write};

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

/// Defualt values for arguments.
enum DefaultValue {
    /// None
    None,
    /// Integer
    Int(usize),
    /// Float
    Float(f64),
}

/// Allowed values for arguments.
enum AllowedValue {
    /// Anything
    None,
    /// Specific list
    Strs(Vec<String>),
}

/// Everything an argument needs.
struct Arg {
    /// Name of the argument.
    name: String,
    /// Text displayed as help for the argument when -h is passed.
    short_help: String,
    /// Text displayed as help for the argument when --help is passed.
    long_help: String,
    /// Single character for flag e.g. -h.
    short_flag: char,
    /// Longer flag e.g. --help.
    long_flag: String,
    /// Whether the argument takes a value.
    takes_value: bool,
    /// Whether the argument takes multiple values
    multiple_values: bool,
    /// Whether the argument had a default value and if so what it is.
    default_value: DefaultValue,
    /// Whether to restrict the argument to specific values and if so what they are.
    allowed_values: AllowedValue,
}

/// Everything about the bca app.
pub struct App {
    /// What arguments it can take, [OPTIONS] in the help.
    options: Vec<Arg>,
    /// The help specific argument.
    help: Arg,
    /// The longest argument name's size.
    max_width: usize,
}

impl App {
    pub fn new() -> Self {
        // Start options
        let options = vec![
            Arg{
                name: String::from("aec"),
                short_help: String::from("Convience flag for reading both AECCARs."),
                long_help: String::from("
\tFlag for reading and summing both the AECCAR0 and AECCAR2 from a VASP calculation"),
                short_flag: 'a',
                long_flag: String::from("aec"),
                takes_value: false,
                multiple_values: false,
                default_value: DefaultValue::None,
                allowed_values: AllowedValue::None,
            },
            Arg {
                name: String::from("file type"),
                short_help: String::from("File containing spin density."),
                long_help: String::from("
\tA path to the spin density associated with the original file. This is primarily
\tfor cube files as if spin density exists in a CHGCAR it will be read automatically.
\tIf using with VASP outputs then the files for charge and spin density must only
\tcontain a single density (ie. the original file has been split)."),
                short_flag: 'f',
                long_flag: String::from("file_type"),
                takes_value: true,
                multiple_values: false,
                default_value: DefaultValue::None,
                allowed_values: AllowedValue::Strs(vec![String::from("cube"), String::from("vasp")]),
            },
            Arg {
                name: String::from("index"),
                short_help: String::from("Index of Bader atoms to be written out."),
                long_help: String::from("
\tAn index of a Bader atom to be written out, starting at 1. This flag requires
\tthe output flag to be set. Multiple atoms can only be written by repeating the
\tflag ie. bca CHGCAR -oi 1 -i 2."),
                short_flag: 'i',
                long_flag: String::from("index"),
                takes_value: true,
                multiple_values: true,
                default_value: DefaultValue::None,
                allowed_values: AllowedValue::None,
            },
            Arg {
                name: String::from("maximum distance"),
                short_help: String::from("Cut-off after which an error will be thrown for distance of Bader maximum to atom."),
                long_help: String::from("
\tThe distance allowed between Bader maximum and its associated atom, in angstrom, before
\tan error is thrown. This will cause a hard crash of the program, consider whether
\tincreasing the cut-off or adding a \"ghost atom\" at the location of the Bader
\tmaximum is more appropriate."),
                short_flag: 'm',
                long_flag: String::from("max_dist"),
                takes_value: true,
                multiple_values: false,
                default_value: DefaultValue::Float(0.1),
                allowed_values: AllowedValue::None,
            },
            Arg {
                name: String::from("output"),
                short_help: String::from("Output the Bader atoms."),
                long_help: String::from("
\tOutput the bader atoms in the same file format as the input density.
\tthis can be used in conjunction with the index flag to specify a
\tspecific atom. without the index flag it will print all the atoms."),
                short_flag: 'o',
                long_flag: String::from("output"),
                takes_value: false,
                multiple_values: false,
                default_value: DefaultValue::None,
                allowed_values: AllowedValue::None,
            },
            Arg {
                name: String::from("reference"),
                short_help: String::from("File(s) containing reference charge."),
                long_help: String::from("
\tA reference charge to do the partitioning upon. Two files can be passed
\tby using multiple flags (bca CHGCAR -r AECCAR0 -r AECCAR2). If two files are
\tpassed they are summed together."),
                short_flag: 'r',
                long_flag: String::from("ref"),
                takes_value: true,
                multiple_values: true,
                default_value: DefaultValue::None,
                allowed_values: AllowedValue::None,
            },
            Arg {
                name: String::from("silent"),
                short_help: String::from("Whether to display any output."),
                long_help: String::from("
\tRuns the program without displaying any output text or progress bars."),
                short_flag: 'x',
                long_flag: String::from("silent"),
                takes_value: false,
                multiple_values: false,
                default_value: DefaultValue::None,
                allowed_values: AllowedValue::None,
            },
            Arg {
                name: String::from("spin"),
                short_help: String::from("File containing spin density."),
                long_help: String::from("
\tA path to the spin density associated with the original file. This is primarily
\tfor cube files as if spin density exists in a CHGCAR it will be read automatically.
\tIf using with VASP outputs then the files for charge and spin density must only
\tcontain a single density (ie. the original file has been split)."),
                short_flag: 's',
                long_flag: String::from("spin"),
                takes_value: true,
                multiple_values: false,
                default_value: DefaultValue::None,
                allowed_values: AllowedValue::None,
            },
            Arg {
                name: String::from("threads"),
                short_help: String::from("Number of threads to distribute the calculation over."),
                long_help: String::from("
\tThe number of threads to be used by the program. A default value of 0 is used
\tto allow the program to best decide how to use the available hardware. It does
\tthis by using the minimum value out of the number cores available and 12."),
                short_flag: 't',
                long_flag: String::from("threads"),
                takes_value: true,
                multiple_values: false,
                default_value: DefaultValue::Int(0),
                allowed_values: AllowedValue::None,
            },
            Arg {
                name: String::from("vacuum tolerance"),
                short_help: String::from("Cut-off at which charge is considered vacuum."),
                long_help: String::from("
\tValues of density below the supplied value are considered vacuum and are not
\tincluded in the calculation."),
                short_flag: 'v',
                long_flag: String::from("vac"),
                takes_value: true,
                multiple_values: false,
                default_value: DefaultValue::Float(1E-6),
                allowed_values: AllowedValue::None,
            },
            Arg {
                name: String::from("weight tolerance"),
                short_help: String::from("Cut-off at which contributions to the weighting will be ignored."),
                long_help: String::from("
\tValues of density below the supplied value are ignored from the weighting and
\tincluded in the calculation. By raising the tolerance the calculation speed can
\tbe increased but every ignored weight is unaccounted for in the final partitions.
\tBe sure to test this!"),
                short_flag: 'w',
                long_flag: String::from("weight"),
                takes_value: true,
                multiple_values: false,
                default_value: DefaultValue::Float(1E-8),
                allowed_values: AllowedValue::None,
            }
        ];
        // End options
        let help = Arg {
            name: String::from("help"),
            short_help: String::from("Print help (see more with '--help')"),
            long_help: String::from("Print help (see a summary with '-h')"),
            short_flag: 'h',
            long_flag: String::from("help"),
            takes_value: false,
            multiple_values: false,
            default_value: DefaultValue::None,
            allowed_values: AllowedValue::None,
        };
        let max_width = options
            .iter()
            .map(|o| {
                1 + o.long_flag.len()
                    + if o.takes_value { 7 + o.name.len() } else { 4 }
            })
            .max()
            .unwrap_or(0);
        Self {
            options,
            help,
            max_width,
        }
    }

    /// Get the text for all the options either in short or long format.
    fn get_options_text(&self, long: bool) -> String {
        self.options.iter().fold(String::new(), |mut output, o| {
            let _ = writeln!(
                output,
                " {:<width$}{}{}{}",
                // the short and long flag
                if o.takes_value {
                    format!("-{},--{} [{}]", o.short_flag, o.long_flag, o.name)
                } else {
                    format!("-{},--{}", o.short_flag, o.long_flag)
                },
                // the help text
                if long {
                    o.long_help.to_string()
                } else {
                    format!("  {}", o.short_help)
                },
                // does it have a default value
                match o.default_value {
                    DefaultValue::None =>
                        if long {
                            String::from("\n")
                        } else {
                            String::with_capacity(0)
                        },
                    DefaultValue::Int(u) =>
                        if long {
                            format!("\n\n\t[default: {}]\n", u)
                        } else {
                            format!(" [default: {}]", u)
                        },
                    DefaultValue::Float(f) =>
                        if long {
                            format!("\n\n\t[default: {}]\n", f)
                        } else {
                            format!(" [default: {}]", f)
                        },
                },
                // does it have allowed values
                match &o.allowed_values {
                    AllowedValue::None => String::with_capacity(0),
                    AllowedValue::Strs(v) =>
                        if long {
                            format!("\n\t[possible values: {:?}]\n", v)
                        } else {
                            format!(" [possible values: {:?}]", v)
                        },
                },
                width = self.max_width
            );
            output
        })
    }

    /// Get the help text in long or short formats.
    fn get_help_text(&self, long: bool) -> String {
        format!(
            " {:<width$}{}",
            format!("-{},--{}", self.help.short_flag, self.help.long_flag),
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
    ) -> Result<&Arg, ArgumentError> {
        if f == self.help.short_flag {
            return Err(ArgumentError::ShortHelp(self));
        }
        for opt in self.options.iter() {
            if f == opt.short_flag {
                return Ok(opt);
            }
        }
        Err(ArgumentError::NotFlag(f.to_string()))
    }

    /// Get the `Arg` from a supplied long flag.
    fn get_option_from_long_flag(
        &self,
        f: String,
    ) -> Result<&Arg, ArgumentError> {
        if f == self.help.long_flag {
            return Err(ArgumentError::LongHelp(self));
        }
        for opt in self.options.iter() {
            if f == opt.long_flag {
                return Ok(opt);
            }
        }
        Err(ArgumentError::NotFlag(f))
    }

    /// Parse arguments from flags and values.
    pub fn parse_args(&self, args: Vec<&str>) -> Result<Args, ArgumentError> {
        let mut arguments = FxHashMap::<String, String>::default();
        let mut multi_arguments = FxHashMap::<String, Vec<String>>::default();
        args.iter().enumerate()
                   .try_for_each(|(i, flag)| {
                       match flag.strip_prefix('-') {
                           // starts with - so is arg
                           Some(stripped_flag) => {
                               // longer than 1 char after being stripped means either cluster
                               // of flags or long flag
                               match stripped_flag.len().cmp(&1) {
                                   std::cmp::Ordering::Greater => match stripped_flag.strip_prefix('-') {
                                       // matches - again is long flag
                                       Some(long_flag) => {
                                           match self.get_option_from_long_flag(long_flag.to_string()) {
                                               // if it's in the argument list we are good and
                                               // just need to check if it takes values
                                               Ok(o) => {if o.takes_value {
                                                   let value = match args.get(i + 1) {
                                                       Some(v) => v.to_string(),
                                                       None => return Err(ArgumentError::NoValue(long_flag.to_string())),
                                                   };
                                                   if o.multiple_values {
                                                       let key = multi_arguments.entry(o.name.to_string()).or_insert(vec![]);
                                                       key.push(value);
                                                   } else {
                                                       arguments.insert(o.name.to_string(), value);
                                                   }
                                               } else {
                                                   arguments.insert(o.name.to_string(), String::with_capacity(0));
                                               };
                                               Ok(())},
                                               // else it isn't an argument
                                               Err(e) => Err(e),
                                           }
                                       },
                                       // doesn't match again so is cluster of short flags
                                       None => {
                                           stripped_flag.chars().try_for_each(|f| match self.get_option_from_short_flag(f) {
                                               // if it's in the argument list we are good and
                                               // just need to check if it takes values
                                               Ok(o) => {if o.takes_value {
                                                   let value = match args.get(i + 1) {
                                                       Some(v) => v.to_string(),
                                                       None => return Err(ArgumentError::NoValue(f.to_string())),
                                                   };
                                                   if o.multiple_values {
                                                       let key = multi_arguments.entry(o.name.to_string()).or_insert(vec![]);
                                                       key.push(value);
                                                   } else {
                                                       arguments.insert(o.name.to_string(), value);
                                                   }
                                               } else {
                                                   arguments.insert(o.name.to_string(), String::with_capacity(0));
                                               };
                                               Ok(())},
                                               // else it isn't an argument
                                               Err(e) => Err(e),
                                           })
                                       },
                                   },
                               // is 1 in length so short flag
                               std::cmp::Ordering::Equal => {
                                   let f = stripped_flag.chars().nth(0).unwrap();
                                   match self.get_option_from_short_flag(f) {
                                       // if it's in the argument list we are good and just
                                       // need to check if it takes values
                                       Ok(o) => {if o.takes_value {
                                           let value = match args.get(i + 1) {
                                               Some(v) => v.to_string(),
                                               None => return Err(ArgumentError::NoValue(stripped_flag.to_string())),
                                           };
                                           if o.multiple_values {
                                               let key = multi_arguments.entry(o.name.to_string()).or_insert(vec![]);
                                               key.push(value);
                                           } else {
                                               arguments.insert(o.name.to_string(), value);
                                           }
                                       } else {
                                           arguments.insert(o.name.to_string(), String::with_capacity(0));
                                       };
                                       Ok(())},
                                       // else it isn't an argument
                                       Err(e) => Err(e),
                                   }
                               },
                               std::cmp::Ordering::Less => Err(ArgumentError::NotFlag("-".to_string()))
                               }
                           },
                           None => Ok(()),
                       }
                   })?;
        let mut file = String::with_capacity(0);
        args.windows(2).for_each(|window| {
            let flag = window[0].to_string();
            let arg = window[1].to_string();
            // doesn't start with - so could be out file
            if !arg.starts_with('-') {
                match flag.strip_prefix('-') {
                    // previous argument does start with - so we need to check if it is a passed value
                    Some(short_flag) => match short_flag.strip_prefix('-') {
                        Some(long_flag) => {
                            if !self
                                .get_option_from_long_flag(
                                    long_flag.to_string(),
                                )
                                .unwrap()
                                .takes_value
                            {
                                file = arg;
                            }
                        }
                        None => {
                            if short_flag.len() > 1 {
                                let mut takes_value_flag = false;
                                for f in short_flag.chars() {
                                    if self
                                        .get_option_from_short_flag(f)
                                        .unwrap()
                                        .takes_value
                                    {
                                        takes_value_flag = true;
                                        break;
                                    }
                                }
                                if !takes_value_flag {
                                    file = arg;
                                }
                            } else if !self
                                .get_option_from_short_flag(
                                    // this should have a flag bu if not pass something that won't
                                    // be matched by any flag
                                    short_flag.chars().nth(0).unwrap_or(','),
                                )
                                .unwrap()
                                .takes_value
                            {
                                file = arg;
                            }
                        }
                    },
                    None => file = arg,
                }
            }
        });
        if file.is_empty() {
            return Err(ArgumentError::NoFile(self));
        }
        let file_type = match arguments.get("file type") {
            Some(ftype) => {
                let ftype = ftype.to_lowercase();
                if ftype.eq("cube") {
                    FileType::Cube
                } else if ftype.eq("vasp") {
                    FileType::Vasp
                } else {
                    return Err(ArgumentError::NotValidValue(
                        String::from("file type"),
                        ftype,
                    ));
                }
            }
            None => {
                let f = file.to_lowercase();
                if f.contains("cube") {
                    FileType::Cube
                } else if f.contains("car") {
                    FileType::Vasp
                } else {
                    eprintln!(
                        "Cannot detect file type, attempting to read as VASP."
                    );
                    FileType::Vasp
                }
            }
        };
        let weight_tolerance = match arguments.get("weight tolerance") {
            Some(w) => match w.parse::<f64>() {
                Ok(f) => f,
                Err(_) => {
                    return Err(ArgumentError::Unparsable(
                        String::from("weight tolerance"),
                        w.to_string(),
                        String::from("f64"),
                    ))
                }
            },
            None => match self
                .get_option_from_short_flag('w')
                .unwrap()
                .default_value
            {
                DefaultValue::Float(f) => f,
                _ => panic!(""),
            },
        };
        let maximum_distance = match arguments.get("maximum distance") {
            Some(m) => match m.parse::<f64>() {
                Ok(f) => f,
                Err(_) => {
                    return Err(ArgumentError::Unparsable(
                        String::from("maximum distance"),
                        m.to_string(),
                        String::from("f64"),
                    ))
                }
            },
            None => match self
                .get_option_from_short_flag('m')
                .unwrap()
                .default_value
            {
                DefaultValue::Float(f) => f,
                _ => panic!(""),
            },
        };
        let vacuum_tolerance = match arguments.get("vacuum tolerance") {
            Some(m) => {
                if m.to_lowercase() == "none" {
                    None
                } else {
                    match m.parse::<f64>() {
                        Ok(f) => Some(f),
                        Err(_) => {
                            return Err(ArgumentError::Unparsable(
                                String::from("vacuum tolerance"),
                                m.to_string(),
                                String::from("f64"),
                            ))
                        }
                    }
                }
            }
            None => match self
                .get_option_from_short_flag('v')
                .unwrap()
                .default_value
            {
                DefaultValue::Float(f) => Some(f),
                _ => panic!(""),
            },
        };
        let mut threads = match arguments.get("threads") {
            Some(t) => match t.parse::<usize>() {
                Ok(f) => f,
                Err(_) => {
                    return Err(ArgumentError::Unparsable(
                        String::from("threads"),
                        t.to_string(),
                        String::from("usize"),
                    ))
                }
            },
            None => match self
                .get_option_from_short_flag('t')
                .unwrap()
                .default_value
            {
                DefaultValue::Int(u) => u,
                _ => panic!(""),
            },
        };
        if threads == 0 {
            threads = num_cpus::get().min(12);
        }
        let output = match arguments.get("output") {
            Some(_) => match multi_arguments.get("index") {
                Some(v) => {
                    let mut a = Vec::with_capacity(v.len());
                    for i in v.iter() {
                        match i.parse::<isize>() {
                            Ok(i) => a.push(i - 1),
                            Err(_) => {
                                return Err(ArgumentError::Unparsable(
                                    String::from("index"),
                                    i.to_string(),
                                    String::from("isize"),
                                ))
                            }
                        }
                    }
                    WriteType::Atom(a)
                }
                None => WriteType::Atom(Vec::with_capacity(0)),
            },
            None => match multi_arguments.get("index") {
                Some(_) => {
                    return Err(ArgumentError::MissingDependant(
                        String::from("index"),
                        String::from("output"),
                    ))
                }
                None => WriteType::None,
            },
        };
        let reference = match multi_arguments.get("reference") {
            Some(v) => match v.len().cmp(&2) {
                std::cmp::Ordering::Less => Reference::One(v[0].to_string()),
                std::cmp::Ordering::Equal => {
                    Reference::Two(v[0].to_string(), v[1].to_string())
                }
                std::cmp::Ordering::Greater => {
                    return Err(ArgumentError::TooManyValues(
                        String::from("reference"),
                        2,
                        v.len(),
                    ))
                }
            },
            None => match arguments.get("aec") {
                Some(_) => match file_type {
                    FileType::Vasp => Reference::Two(
                        String::from("AECCAR0"),
                        String::from("AECCAR2"),
                    ),
                    FileType::Cube => {
                        return Err(ArgumentError::WrongFileType(
                            String::from("aec"),
                            String::from("cube"),
                        ))
                    }
                },
                None => Reference::None,
            },
        };
        let spin = arguments.get("spin").cloned();
        let silent = arguments.get("silent").is_some();
        Ok(Args {
            file,
            file_type,
            weight_tolerance,
            maximum_distance,
            output,
            reference,
            silent,
            spin,
            threads,
            vacuum_tolerance,
        })
    }
}

// To make clippy happy.
impl Default for App {
    fn default() -> Self {
        Self::new()
    }
}

impl Display for App {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let usage = String::from("Usage: bca [OPTIONS] <file>");
        let arguments =
            String::from("Arguments:\n <file>  The file to analyse.");
        let options = self.get_options_text(false);
        let help = self.get_help_text(false);
        write!(
            f,
            "{}\n\n{}\n\nOptions:\n{}{}",
            usage, arguments, options, help
        )
    }
}

impl Debug for App {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let usage = String::from("Usage: bca [OPTIONS] <file>");
        let arguments =
            String::from("Arguments:\n <file>  The file containing the charge density to analyse.");
        let options = self.get_options_text(true);
        let help = self.get_help_text(true);
        write!(
            f,
            "{}\n\n{}\n\nOptions:\n{}{}",
            usage, arguments, options, help
        )
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
    /// Should the program be ran silently.
    pub silent: bool,
    /// Is there a spin density to include as well.
    pub spin: Option<String>,
    /// How many threads to use in the calculation.
    pub threads: usize,
    /// Is there a tolerance to consider a density vacuum.
    pub vacuum_tolerance: Option<f64>,
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
        let flag = matches!(args.file_type, FileType::Vasp);
        assert!(flag);
    }

    #[test]
    fn argument_file_type_default_unknown() {
        let app = App::new();
        let v = vec!["bca", "CHG"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.file_type, FileType::Vasp);
        assert!(flag);
    }

    #[test]
    fn argument_file_type_vasp() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR", "-f", "vasp"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.file_type, FileType::Vasp);
        assert!(flag);
    }

    #[test]
    fn argument_file_type_default_cube() {
        let app = App::new();
        let v = vec!["bca", "charge.cube"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.file_type, FileType::Cube);
        assert!(flag);
    }

    #[test]
    fn argument_file_type_cube() {
        let app = App::new();
        let v = vec!["bca", "charge.cube", "--file_type", "cube"];
        let args = app.parse_args(v).unwrap();
        let flag = matches!(args.file_type, FileType::Cube);
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
        assert_eq!(args.vacuum_tolerance, Some(1E-4))
    }

    #[test]
    fn argument_vacuum_tolerance_default() {
        let app = App::new();
        let v = vec!["bca", "CHGCAR"];
        let args = app.parse_args(v).unwrap();
        assert_eq!(args.vacuum_tolerance, Some(1E-6))
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
        let threads = num_cpus::get().min(12);
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
}
