use crate::errors::VacuumError;

/// compute the cross product between two vectors
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// subtract two \[f64; 3\] vectors
pub fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// compute the product of a row vector and a matrix, v.M
pub fn dot(v: [f64; 3], m: [[f64; 3]; 3]) -> [f64; 3] {
    [0, 1, 2].map(|i| v[0] * m[0][i] + v[1] * m[1][i] + v[2] * m[2][i])
}

/// compute v.M rounded to the nearest integers
pub fn idot(v: [f64; 3], m: [[f64; 3]; 3]) -> [isize; 3] {
    dot(v, m).map(|x| x.round() as isize)
}

/// compute the dot product between two vectors
pub fn vdot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b.iter()).map(|(a, b)| a * b).sum()
}

/// compute the squared norm of a vector
pub fn norm_squared(a: [f64; 3]) -> f64 {
    vdot(a, a)
}

/// compute the norm of a vector
pub fn norm(a: [f64; 3]) -> f64 {
    norm_squared(a).sqrt()
}

/// compute the sum of two vectors
pub fn vsum(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// calculates the inverse of a 3x3 lattice if it is invertible
pub fn invert_lattice(lattice: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let minor00 = lattice[1][1] * lattice[2][2] - lattice[1][2] * lattice[2][1];
    let minor01 = lattice[1][0] * lattice[2][2] - lattice[1][2] * lattice[2][0];
    let minor02 = lattice[1][0] * lattice[2][1] - lattice[1][1] * lattice[2][0];
    let determinant = lattice[0][0] * minor00 - lattice[0][1] * minor01
        + lattice[0][2] * minor02;
    // a determinant of zero is not invertible
    if determinant.abs() <= f64::EPSILON {
        None
    } else {
        Some([
            [
                minor00 / determinant,
                (lattice[0][2] * lattice[2][1] - lattice[2][2] * lattice[0][1])
                    / determinant,
                (lattice[0][1] * lattice[1][2] - lattice[1][1] * lattice[0][2])
                    / determinant,
            ],
            [
                -minor01 / determinant,
                (lattice[0][0] * lattice[2][2] - lattice[2][0] * lattice[0][2])
                    / determinant,
                (lattice[0][2] * lattice[1][0] - lattice[1][2] * lattice[0][0])
                    / determinant,
            ],
            [
                minor02 / determinant,
                (lattice[0][1] * lattice[2][0] - lattice[2][1] * lattice[0][0])
                    / determinant,
                (lattice[0][0] * lattice[1][1] - lattice[1][0] * lattice[0][1])
                    / determinant,
            ],
        ])
    }
}

/// Returns the indices of every voxel above the vacuum `tolerance`, sorted
/// from the highest density to the lowest.
///
/// # Errors
/// Returns a [`VacuumError`] if no voxel is above the tolerance.
pub fn index_generator(
    density: &[f64],
    tolerance: f64,
) -> Result<Vec<usize>, VacuumError> {
    let mut index: Vec<usize> = (0..density.len())
        .filter(|&i| density[i] > tolerance)
        .collect();
    if index.is_empty() {
        return Err(VacuumError {
            vacuum_tolerance: tolerance,
            density: density.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        });
    }
    index.sort_unstable_by(|a, b| density[*b].total_cmp(&density[*a]));
    Ok(index)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn utils_dot() {
        assert_eq!(
            dot([0.5, 0.5, 0.5], [[1., 0., 0.], [1., 2., 0.], [2., 2., 4.]]),
            [2., 2., 2.]
        )
    }

    #[test]
    fn utils_vdot() {
        assert_eq!(vdot([1., 2., 3.], [1., 2., 3.]), 14.)
    }

    #[test]
    fn utils_norm() {
        assert_eq!(norm([3., 4., 12.]), 13.)
    }

    #[test]
    fn utils_norm_squared_and_vsum() {
        assert_eq!(norm_squared([3., 4., 12.]), 169.);
        assert_eq!(vsum([1., 2., 3.], [4., 5., 6.]), [5., 7., 9.]);
        assert_eq!(subtract([4., 5., 6.], [1., 2., 3.]), [3., 3., 3.]);
    }

    #[test]
    fn utils_index_generator_high() {
        let data = (0..60).map(|x| x as f64).collect::<Vec<f64>>();
        assert!(index_generator(&data, 100.).is_err())
    }

    #[test]
    fn utils_index_generator_low() {
        let data = (0..60).map(|x| x as f64).collect::<Vec<f64>>();
        let index = index_generator(&data, 1.0).unwrap();
        assert_eq!(index[0], 59);
        assert_eq!(index.len(), 58);
    }

    #[test]
    fn utils_cross_product() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert_eq!(cross(a, b), [0.0, 0.0, 1.0]);
    }

    #[test]
    fn utils_matrix_dot() {
        let v = [1.0, 1.0, 0.0];
        let m = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
        // Vector * Matrix
        assert_eq!(dot(v, m), [2.0, 2.0, 0.0]);
    }

    #[test]
    fn utils_idot_rounding() {
        // Test that floating point coordinates round correctly to integer shifts
        let v_approx = [0.99, 1.01, -0.01];
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

        // Should round to [1, 1, 0]
        assert_eq!(idot(v_approx, identity), [1, 1, 0]);
    }

    #[test]
    fn utils_invert_lattice() {
        // 1. Identity
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let inv = invert_lattice(&identity).unwrap();
        assert_eq!(inv, identity);

        // 2. Scaling (Inverse of 2I is 0.5I)
        let scaled = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
        let inv_scaled = invert_lattice(&scaled).unwrap();
        assert!((inv_scaled[0][0] - 0.5).abs() < 1e-9);

        // 3. Singular (Determinant = 0)
        let singular = [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(invert_lattice(&singular).is_none());
    }
}
