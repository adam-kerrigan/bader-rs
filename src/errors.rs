use crate::arguments::App;
use std::fmt::{Debug, Display};

/// An error for not being able to assign Bader maxima to atoms
pub struct MaximaError {
    /// The position of the maxima
    pub maximum: [f64; 3],
    /// The index of the atom
    pub atom: usize,
    /// the distance between maximum and atom
    pub distance: f64,
}

/// Display should never be called explicitly as the Bader maximum needs to be formatted by the
/// file type
impl Display for MaximaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "")
    }
}

/// Debug should never be called explicitly as the Bader maximum needs to be formatted by the
/// file type
impl Debug for MaximaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "")
    }
}

/// An error for the vacuum tolerance being heigher than the highest density value
pub struct VacuumError {
    pub vacuum_tolerance: f64,
    pub density: f64,
}

impl Display for VacuumError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Vacuum tolerance ({}) is higher than maximum value of density ({}).",
            self.vacuum_tolerance, self.density
        )
    }
}

impl Debug for VacuumError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

/// Error for reading of Arguments.
pub enum ArgumentError<'a> {
    /// Passed an argument that isn't a flag.
    NotFlag(String),
    /// Passed a value that isn't parsable.
    /// Unparsable(flag, value, type)
    Unparsable(String, String, String),
    /// Didn't pass a value.
    /// NoValue(flag)
    NoValue(String),
    /// Passed too many values.
    /// TooManyValues(flag, max, supplied)
    TooManyValues(String, usize, usize),
    /// Passed an unvalid value for the flag, ie. filetype.
    /// NotValidValue(flag, value)
    NotValidValue(String, String),
    /// Missing a dependant flag.
    /// MissingDependant(flag, dependancy)
    MissingDependant(String, String),
    /// Passed a flag for a different filetype than the one given.
    /// WrongFileType(flag, filetype)
    WrongFileType(String, String),
    /// No file given.
    NoFile(&'a App),
    /// Asked for help with -h.
    ShortHelp(&'a App),
    /// Asked for help with --help.
    LongHelp(&'a App),
}

impl Display for ArgumentError<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NotFlag(flag) => {
                write!(f, "The flag: {} does not exist", flag)
            }
            Self::Unparsable(flag, value, typ) => write!(
                f,
                "The supplied value \"{}\" for the option \"{}\" is unparsable as a {}.",
                value, flag, typ
            ),
            Self::NoValue(flag) => write!(
                f,
                "The option \"{}\" requires a value to be supplied.",
                flag
            ),
            Self::TooManyValues(flag, max, supplied) => write!(
                f,
                "Option \"{}\" was supplied {} times, the maximum allowed is \"{}\".",
                flag, supplied, max
            ),
            Self::NotValidValue(flag, value) => write!(
                f,
                "The value \"{}\" is not valid input for the option \"{}\".",
                value, flag
            ),
            Self::MissingDependant(flag, dependant) => write!(
                f,
                "The option \"{}\" requires the option \"{}\" to also be set.",
                flag, dependant
            ),
            Self::WrongFileType(flag, file_type) => write!(
                f,
                "The option \"{}\" cannot be set for the file type \"{}\".",
                flag, file_type
            ),
            Self::NoFile(app) => write!(f, "No file supplied.\n\n{}", app),
            Self::ShortHelp(app) => write!(f, "{}", app),
            Self::LongHelp(app) => write!(f, "{:?}", app),
        }
    }
}

impl Debug for ArgumentError<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}
