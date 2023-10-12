use crate::utils;
use std::cmp::Ordering::Less;

/// struct for containing the information about the atoms.
pub struct Atoms {
    /// The lattice of the structure.
    pub lattice: Lattice,
    /// The positions of the atoms in cartesian coordinates.
    pub positions: Vec<[f64; 3]>,
    /// Text representation from the input file.
    pub text: String,
    /// The positions of the atoms in the LLL-reduced basis.
    pub reduced_positions: Vec<[f64; 3]>,
}

impl Atoms {
    /// Initialises the structure.
    pub fn new(
        lattice: Lattice,
        positions: Vec<[f64; 3]>,
        text: String,
    ) -> Self {
        let reduced_positions = positions
            .iter()
            .map(|p| lattice.cartesian_to_reduced(*p))
            .collect::<Vec<[f64; 3]>>();
        Self {
            lattice,
            positions,
            text,
            reduced_positions,
        }
    }
}

/// Lattice - structure for containing information on the cell
///
/// <pre class="rust">
/// shift matrix ordering:
///     0 -> (-1,-1,-1)   7 -> (-1, 1, 0)  14 -> (0, 0, 1)  21 -> (1, 0,-1)
///     1 -> (-1,-1, 0)   8 -> (-1, 1, 1)  15 -> (0, 1,-1)  22 -> (1, 0, 0)
///     2 -> (-1,-1, 1)   9 -> (0,-1,-1)   16 -> (0, 1, 0)  23 -> (1, 0, 1)
///     3 -> (-1, 0,-1)  10 -> (0,-1, 0)   17 -> (0, 1, 1)  24 -> (1, 1,-1)
///     4 -> (-1, 0, 0)  11 -> (0,-1, 1)   18 -> (1,-1,-1)  25 -> (1, 1, 0)
///     5 -> (-1, 0, 1)  12 -> (0, 0,-1)   19 -> (1,-1, 0)  26 -> (1, 1, 1)
///     6 -> (-1, 1,-1)  13 -> (0, 0, 0)   20 -> (1,-1, 1)
/// </pre>
pub struct Lattice {
    /// The cartesian vectors for every combination of lattice vector.
    pub cartesian_shift_matrix: [[f64; 3]; 27],
    /// Transformation matrix for converting to fractional coordinates.
    pub to_fractional: [[f64; 3]; 3],
    /// Transformation matrix for converting to cartesian coordinates.
    pub to_cartesian: [[f64; 3]; 3],
    /// The cartesian vectors for every combination of reduced lattice vector.
    pub reduced_cartesian_shift_matrix: [[f64; 3]; 27],
    /// The conversion of the reduced shift matrix to the individual steps in the
    /// [`crate::grid::Grid`]
    pub reduced_grid_shift_matrix: Vec<Vec<usize>>,
    /// Transformation matrix for converting to fractional coordinates.
    pub reduced_to_fractional: [[f64; 3]; 3],
    /// Transformation matrix for converting to cartesian coordinates.
    pub reduced_to_cartesian: [[f64; 3]; 3],
    /// Volume of the lattice.
    pub volume: f64,
}

impl Lattice {
    /// Initialises the structure. Builds all the fields of the lattice structure
    /// from a 2d vector in the form:
    ///
    /// <pre class="rust">
    /// [
    ///     [ax, ay, az],
    ///     [bx, by, bz],
    ///     [cx, cy, cz],
    /// ]
    /// </pre>
    pub fn new(lattice: [[f64; 3]; 3]) -> Self {
        let cartesian_shift_matrix =
            Lattice::create_cartesian_shift_matrix(&lattice);
        let to_fractional = match utils::invert_lattice(&lattice) {
            Some(inv) => inv,
            None => panic!("Supplied lattice does not span 3d space."),
        };
        let reduced_lattice = lll_lattice(lattice);
        let reduced_cartesian_shift_matrix =
            Lattice::create_cartesian_shift_matrix(&reduced_lattice);
        let reduced_to_fractional =
            match utils::invert_lattice(&reduced_lattice) {
                Some(inv) => inv,
                None => panic!("Supplied lattice does not span 3d space."),
            };
        let reduced_grid_shift_matrix = Lattice::create_grid_shift_matrix(
            &reduced_cartesian_shift_matrix,
            &reduced_to_fractional,
        );
        let volume =
            utils::vdot(lattice[0], utils::cross(lattice[1], lattice[2])).abs();
        let to_cartesian = lattice;
        let reduced_to_cartesian = reduced_lattice;
        Self {
            cartesian_shift_matrix,
            to_fractional,
            to_cartesian,
            reduced_cartesian_shift_matrix,
            reduced_grid_shift_matrix,
            reduced_to_fractional,
            reduced_to_cartesian,
            volume,
        }
    }

    /// Turn fractional coordinates into Cartesian coordinates in the reduced basis.
    pub fn fractional_to_reduced(&self, p: [f64; 3]) -> [f64; 3] {
        self.cartesian_to_reduced(utils::dot(p, self.to_cartesian))
    }

    /// Map Cartesian coordinates into the reduced basis.
    pub fn cartesian_to_reduced(&self, p: [f64; 3]) -> [f64; 3] {
        let pn = utils::dot(p, self.reduced_to_fractional)
            .iter()
            .map(|p| p.rem_euclid(1.0))
            .collect::<Vec<f64>>()
            .try_into()
            .unwrap();
        utils::dot(pn, self.reduced_to_cartesian)
    }

    /// Create the shift matrix from the lattice supplied.
    fn create_cartesian_shift_matrix(
        lattice: &[[f64; 3]; 3],
    ) -> [[f64; 3]; 27] {
        let x = lattice[0];
        let y = lattice[1];
        let z = lattice[2];
        [
            [
                -x[0] - y[0] - z[0],
                -x[1] - y[1] - z[1],
                -x[2] - y[2] - z[2],
            ],
            [-x[0] - y[0], -x[1] - y[1], -x[2] - y[2]],
            [
                -x[0] - y[0] + z[0],
                -x[1] - y[1] + z[1],
                -x[2] - y[2] + z[2],
            ],
            [-x[0] - z[0], -x[1] - z[1], -x[2] - z[2]],
            [-x[0], -x[1], -x[2]],
            [-x[0] + z[0], -x[1] + z[1], -x[2] + z[2]],
            [
                -x[0] + y[0] - z[0],
                -x[1] + y[1] - z[1],
                -x[2] + y[2] - z[2],
            ],
            [-x[0] + y[0], -x[1] + y[1], -x[2] + y[2]],
            [
                -x[0] + y[0] + z[0],
                -x[1] + y[1] + z[1],
                -x[2] + y[2] + z[2],
            ],
            [-y[0] - z[0], -y[1] - z[1], -y[2] - z[2]],
            [-y[0], -y[1], -y[2]],
            [-y[0] + z[0], -y[1] + z[1], -y[2] + z[2]],
            [-z[0], -z[1], -z[2]],
            [0.0, 0.0, 0.0],
            [z[0], z[1], z[2]],
            [y[0] - z[0], y[1] - z[1], y[2] - z[2]],
            [y[0], y[1], y[2]],
            [y[0] + z[0], y[1] + z[1], y[2] + z[2]],
            [x[0] - y[0] - z[0], x[1] - y[1] - z[1], x[2] - y[2] - z[2]],
            [x[0] - y[0], x[1] - y[1], x[2] - y[2]],
            [x[0] - y[0] + z[0], x[1] - y[1] + z[1], x[2] - y[2] + z[2]],
            [x[0] - z[0], x[1] - z[1], x[2] - z[2]],
            [x[0], x[1], x[2]],
            [x[0] + z[0], x[1] + z[1], x[2] + z[2]],
            [x[0] + y[0] - z[0], x[1] + y[1] - z[1], x[2] + y[2] - z[2]],
            [x[0] + y[0], x[1] + y[1], x[2] + y[2]],
            [x[0] + y[0] + z[0], x[1] + y[1] + z[1], x[2] + y[2] + z[2]],
        ]
    }

    /// Turn the shift matrix into a vector of all the required steps in the [`crate::grid::Grid`]
    /// required to move by the vector.
    fn create_grid_shift_matrix(
        shift_matrix: &[[f64; 3]; 27],
        to_fractional: &[[f64; 3]; 3],
    ) -> Vec<Vec<usize>> {
        shift_matrix
            .iter()
            .map(|c_shift| {
                let shift = utils::idot(*c_shift, *to_fractional);
                // how many times are we going to have to reduce the vector
                let max = shift.iter().map(|x| x.abs()).max().unwrap();
                (0..max)
                    .map(|i| {
                        let out = shift
                            .iter()
                            .map(|s| {
                                // if the value is 0 or below we have
                                // finshed reducing this axis
                                if let Less = (s.abs() - i).cmp(&1) {
                                    0
                                // if it is 1 or above then we need to
                                // add a 1 with the same sign as the value
                                } else {
                                    s.signum()
                                }
                            })
                            .collect::<Vec<isize>>();
                        (out[0] * 9 + out[1] * 3 + out[2] + 13) as usize
                    })
                    .collect()
            })
            .collect()
    }
}

/// Calculates the lll reduction of a lattice.
pub fn lll_lattice(lattice: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let delta = 0.75;
    let mut a = lattice;
    let (mut b, mut mu) = gram_schmidt(&a);
    let mut i = 1usize;
    while i <= 2 {
        for j in (0..i).rev() {
            match mu[i][j] {
                q if q.abs() <= 0.5 => (),
                q => {
                    for k in 0..3 {
                        a[i][k] -= q.round() * a[j][k];
                    }
                    let (b_temp, mu_temp) = gram_schmidt(&a);
                    b = b_temp;
                    mu = mu_temp;
                }
            }
        }
        if utils::vdot(b[i], b[i])
            >= (delta - mu[i][i - 1].powi(2)) * utils::vdot(b[i - 1], b[i - 1])
        {
            i += 1;
        } else {
            for j in 0..3 {
                b[0][0] = a[i][j];
                a[i][j] = a[i - 1][j];
                a[i - 1][j] = b[0][0];
            }
            let (b_temp, mu_temp) = gram_schmidt(&a);
            b = b_temp;
            mu = mu_temp;
            i = 1usize.max(i - 1);
        }
    }
    a
}

/// Calculates the Gram-Schmidt co-effecients for the lll-reduction.
fn gram_schmidt(v: &[[f64; 3]; 3]) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let mut u = [[0f64; 3]; 3];
    let mut mu = [[0f64; 3]; 3];
    u[0] = [v[0][0], v[0][1], v[0][2]];
    mu[1][0] = utils::vdot(v[1], u[0]) / utils::vdot(u[0], u[0]);
    for i in 0..3 {
        u[1][i] = v[1][i] - (mu[1][0] * u[0][i]);
    }
    mu[2][0] = utils::vdot(v[2], u[0]) / utils::vdot(u[0], u[0]);
    mu[2][1] = utils::vdot(v[2], u[1]) / utils::vdot(u[1], u[1]);
    for i in 0..3 {
        u[2][i] = v[2][i] - (mu[2][0] * u[0][i]) - (mu[2][1] * u[1][i]);
    }
    (u, mu)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atoms_new() {
        let positions = vec![[0.; 3]];
        let lattice = Lattice::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]);
        let text = String::new();
        let atoms = Atoms::new(lattice, positions, text);
        let positions = vec![[0.; 3]];
        let lattice = Lattice::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]);
        let text = String::new();
        assert_eq!(atoms.lattice.to_cartesian, lattice.to_cartesian);
        assert_eq!(atoms.positions, positions);
        assert_eq!(atoms.text, text);
    }

    #[test]
    #[should_panic]
    fn lattice_new_non_invert() {
        let _ = Lattice::new([[1., 0., 0.], [1., 0., 0.], [0., 0., 2.]]);
    }
}
