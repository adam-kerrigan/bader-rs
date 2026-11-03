use crate::arguments::App;
use std::fmt::{Debug, Display};

/// An error for not being able to assign a Bader maximum to an atom.
///
/// The position is in the Cartesian frame of the [`Grid`](crate::grid::Grid),
/// the binary formats it for the file type before reporting it.
#[derive(Debug)]
pub struct MaximaError {
    /// The position of the maximum
    pub maximum: [f64; 3],
    /// The index of the nearest atom
    pub atom: usize,
    /// The distance between the maximum and the atom
    pub distance: f64,
}

impl Display for MaximaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Bader maximum at {:?} is too far away from nearest atom: {} with a distance of {} Ang.",
            self.maximum,
            self.atom + 1,
            self.distance
        )
    }
}

/// An error for the vacuum tolerance being higher than the highest density value
#[derive(Debug)]
pub struct VacuumError {
    /// The tolerance that was passed.
    pub vacuum_tolerance: f64,
    /// The highest value in the density.
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
    /// Passed an invalid value for the flag, ie. filetype.
    /// NotValidValue(flag, value)
    NotValidValue(String, String),
    /// Missing a flag this one depends on.
    /// MissingDependant(flag, dependency)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_messages() {
        let e = MaximaError {
            maximum: [1.0, 2.0, 3.0],
            atom: 3,
            distance: 1.5,
        };
        assert_eq!(
            e.to_string(),
            "Bader maximum at [1.0, 2.0, 3.0] is too far away from nearest atom: 4 with a distance of 1.5 Ang."
        );
        let e = VacuumError {
            vacuum_tolerance: 1.0,
            density: 0.5,
        };
        assert_eq!(
            e.to_string(),
            "Vacuum tolerance (1) is higher than maximum value of density (0.5)."
        );
        let e = ArgumentError::TooManyValues(String::from("ref"), 2, 3);
        assert_eq!(
            e.to_string(),
            "Option \"ref\" was supplied 3 times, the maximum allowed is \"2\"."
        );
        assert_eq!(format!("{:?}", e), e.to_string());
    }
}
