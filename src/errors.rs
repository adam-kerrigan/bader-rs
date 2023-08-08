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
        write!(f, "Vacuum tolerance ({}) is higher than maximum value of density ({}).", self.vacuum_tolerance, self.density)
    }
}

impl Debug for VacuumError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Vacuum tolerance ({}) is higher than maximum value of density ({}).", self.vacuum_tolerance, self.density)
    }
}
