use crate::utils;
use std::cmp::Ordering::Less;

/// Stores the configuration of the atomic system.
///
/// This structure holds both the physical positions of atoms (in Cartesian coordinates)
/// and their positions in the reduced lattice basis (used for internal calculations).
///
/// # Examples
/// ```
/// use bader::atoms::{Atoms, Lattice};
///
/// let lattice = Lattice::new([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]);
/// let positions = vec![[15.0, 5.0, 5.0]];
///
/// let atoms = Atoms::new(lattice, positions, "".to_string());
///
/// // Original Cartesian positions
/// assert_eq!(atoms.positions[0], [15.0, 5.0, 5.0]);
///
/// // Reduced positions (wrapped into the LLL-reduced lattice)
/// assert_eq!(atoms.reduced_positions[0], [5.0, 5.0, 5.0]);
/// ```
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

/// Manages the crystal lattice, coordinate transformations, and Periodic Boundary Conditions (PBC).
///
/// Upon initialization, this struct automatically performs LLL reduction to find a
/// "better" (more orthogonal/shorter) basis set. This reduced basis is used for
/// internal algorithms to improve numerical stability and neighbor finding.
///
/// # Key Features
/// * **Shift Matrices**: Pre-calculated vectors for visiting 26 nearest neighbors (PBC aware).
/// * **LLL Reduction**: Automatic basis reduction.
/// * **Coordinate Transforms**: Cartesian <-> Fractional <-> Reduced.
///
/// # Shift Matrix Ordering
/// Neighbors are indexed 0-26. Center (0,0,0) is index 13.
/// <pre>
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
    /// Initializes a new Lattice from a 3x3 matrix.
    ///
    /// The input matrix `lattice` should be provided as `[a, b, c]` row vectors.
    ///
    /// # Panics
    /// Panics if the provided lattice vectors are linearly dependent (volume is zero),
    /// meaning they do not span 3D space.
    ///
    /// # Logic
    /// 1. Computes the inverse and volume of the input lattice.
    /// 2. Performs LLL reduction to create a `reduced_lattice`.
    /// 3. Generates shift matrices for finding neighbors in the reduced basis.
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

    pub fn fractional_to_cartesian(&self, p: [f64; 3]) -> [f64; 3] {
        utils::dot(p, self.to_cartesian)
    }

    /// Map Cartesian coordinates into the reduced basis.
    pub fn cartesian_to_reduced(&self, p: [f64; 3]) -> [f64; 3] {
        let pn = utils::dot(p, self.reduced_to_fractional)
            .iter()
            .map(|p| p - p.floor())
            .collect::<Vec<f64>>()
            .try_into()
            .unwrap();
        utils::dot(pn, self.reduced_to_cartesian)
    }

    /// Calculates the shortest distance between two points under Periodic Boundary Conditions.
    ///
    /// This method checks the distance between `a` and `b` as well as all 26 periodic
    /// images of `b` defined by the reduced shift matrix.
    ///
    /// # Arguments
    /// * `a` - Point A (Cartesian).
    /// * `b` - Point B (Cartesian).
    /// * `min_dist` - An initial upper bound.
    ///
    /// # Returns
    /// The scalar Euclidean distance.
    pub fn minimum_distance(
        &self,
        a: [f64; 3],
        b: [f64; 3],
        min_dist: Option<f64>,
    ) -> f64 {
        let mut min_dist = min_dist.unwrap_or(f64::INFINITY);
        for periodic_shift in self.reduced_cartesian_shift_matrix.iter() {
            let distance = {
                (a[0] - (b[0] + periodic_shift[0])).powi(2)
                    + (a[1] - (b[1] + periodic_shift[1])).powi(2)
                    + (a[2] - (b[2] + periodic_shift[2])).powi(2)
            };
            if distance < min_dist {
                min_dist = distance;
            }
        }
        min_dist
    }

    pub fn closest_image(&self, a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        let mut min_dist = f64::INFINITY;
        let mut position = [0.0; 3];
        for periodic_shift in self.reduced_cartesian_shift_matrix.iter() {
            let image_position = b
                .iter()
                .zip(periodic_shift)
                .map(|(f, image)| *f + *image)
                .collect::<Vec<f64>>();
            let distance = a
                .iter()
                .zip(&image_position)
                .fold(0.0, |acc, (f, p)| acc + (f - p).powi(2));
            if distance < min_dist {
                min_dist = distance;
                position = image_position.try_into().unwrap();
            }
        }
        position
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
    use crate::utils;

    // --- Lattice Tests ---

    #[test]
    fn test_lattice_cubic() {
        // Simple cubic lattice 10x10x10
        let matrix = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let lattice = Lattice::new(matrix);

        assert_eq!(lattice.volume, 1000.0);

        // In a cubic lattice, reduced basis is identical to input
        assert_eq!(lattice.to_cartesian, lattice.reduced_to_cartesian);
    }

    #[test]
    #[should_panic(expected = "Supplied lattice does not span 3d space")]
    fn test_lattice_singular_panic() {
        // Two vectors are identical -> Volume 0 -> Singular
        Lattice::new([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
    }

    #[test]
    fn test_coordinate_transformations() {
        let matrix = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
        let lattice = Lattice::new(matrix);

        // 1. Fractional to Cartesian
        let frac = [0.5, 0.5, 0.5];
        let cart = lattice.fractional_to_cartesian(frac);
        assert_eq!(cart, [1.0, 1.0, 1.0]);

        // 2. Cartesian to Reduced
        // For a cubic lattice, reduced basis == original basis.
        // Input: [1.0, 1.0, 1.0] (Center of box)
        // Reduced: Should be the same Cartesian vector [1.0, 1.0, 1.0]
        let p_cart = [1.0, 1.0, 1.0];
        let p_red = lattice.cartesian_to_reduced(p_cart);
        assert_eq!(p_red, [1.0, 1.0, 1.0]);

        // 3. Cartesian to Reduced (PBC Wrapping)
        // Input: [3.0, 1.0, 1.0] (Outside box in X)
        // Reduced: Should wrap back to [1.0, 1.0, 1.0] inside the 2x2x2 cell
        let p_outside = [3.0, 1.0, 1.0];
        let p_wrapped = lattice.cartesian_to_reduced(p_outside);
        assert!((p_wrapped[0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_minimum_distance_pbc() {
        // 10x10x10 box
        let matrix = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let lattice = Lattice::new(matrix);

        let p1 = [1.0, 1.0, 1.0];
        let p2 = [9.0, 1.0, 1.0];

        // Real distance is 8.0, but PBC distance is 2.0 (wrapping around x)
        let dist = lattice.minimum_distance(p1, p2, None).powf(0.5);
        assert!(
            (dist - 2.0).abs() < 1e-6,
            "PBC distance failed, got {}",
            dist
        );
    }

    // --- LLL Reduction Test (Wikipedia Example) ---

    #[test]
    fn test_lll_reduction_wikipedia() {
        // From Wikipedia: "Lenstra–Lenstra–Lovász lattice basis reduction algorithm"
        // Basis given by COLUMNS of:
        // [ 1, -1, 3 ]
        // [ 1,  0, 5 ]
        // [ 1,  2, 6 ]
        //
        // Transposing to Row-Major format for our struct:
        let input_lattice = [
            [1.0, 1.0, 1.0],  // Col 1
            [-1.0, 0.0, 2.0], // Col 2
            [3.0, 5.0, 6.0],  // Col 3
        ];

        // Expected Reduced Basis given by COLUMNS of:
        // [ 0, 1, -1 ]
        // [ 1, 0,  0 ]
        // [ 0, 1,  2 ]
        //
        // Transposing to Row-Major format:
        // Row 1: [0, 1, 0]
        // Row 2: [1, 0, 1]
        // Row 3: [-1, 0, 2]

        let lattice = Lattice::new(input_lattice);
        let reduced = lattice.reduced_to_cartesian;

        // Note: The LLL algorithm guarantees a specific reduced quality, but the
        // ordering of vectors might vary slightly depending on implementation details
        // (though usually standard LLL is deterministic).
        // We check if the rows match the expected set.

        let expected_rows =
            vec![[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [-1.0, 0.0, 2.0]];

        for row in reduced.iter() {
            let found = expected_rows.iter().any(|exp| {
                // Check for equality (allowing for sign flips on the whole vector)
                let diff_pos = utils::subtract(*row, *exp);
                let diff_neg =
                    utils::subtract(*row, [-exp[0], -exp[1], -exp[2]]);

                utils::vdot(diff_pos, diff_pos) < 1e-6
                    || utils::vdot(diff_neg, diff_neg) < 1e-6
            });

            assert!(found, "Unexpected reduced vector: {:?}", row);
        }

        // Check volumes match (determinant invariant)
        // Vol input = |1(0-10) - 1(-6-6) + 1(-5-0)| = |-10 + 12 - 5| = |-3| = 3.0
        // Vol output = 1.0 * 1.0 * sqrt(5) ... easier to just check lattice.volume
        assert!((lattice.volume - 3.0).abs() < 1e-6);
    }

    // --- Atoms Tests ---

    #[test]
    fn test_atoms_initialization() {
        let matrix = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let lattice = Lattice::new(matrix);

        // Atom at [15.0, 5.0, 5.0] (Outside the 10x10x10 box)
        let positions = vec![[15.0, 5.0, 5.0]];
        let text = "Test Atom".to_string();

        let atoms = Atoms::new(lattice, positions.clone(), text.clone());

        // 1. Check stored Cartesian positions (Unchanged)
        assert_eq!(atoms.positions, positions);

        // 2. Check Reduced Positions (Cartesian coordinates in reduced cell)
        // Should wrap 15.0 -> 5.0
        let red = atoms.reduced_positions[0];
        assert!((red[0] - 5.0).abs() < 1e-6);
        assert!((red[1] - 5.0).abs() < 1e-6);
        assert!((red[2] - 5.0).abs() < 1e-6);
    }
}
