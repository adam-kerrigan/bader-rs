use crate::errors::VacuumError;

/// compute the cross product between two vectors
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// subtract two [f64;3] vectors
pub fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// compute the dot product between a vector and a matrix
pub fn dot(v: [f64; 3], m: [[f64; 3]; 3]) -> [f64; 3] {
    (0..3)
        .map(|i| (v[0] * m[0][i] + v[1] * m[1][i] + v[2] * m[2][i]))
        .collect::<Vec<f64>>()
        .try_into()
        .unwrap() // safe to unwrap as is size 3
}

/// compute the integer dot product between a vector and a matrix
pub fn idot(v: [f64; 3], m: [[f64; 3]; 3]) -> [isize; 3] {
    (0..3)
        .map(|i| {
            (v[0] * m[0][i] + v[1] * m[1][i] + v[2] * m[2][i]).round() as isize
        })
        .collect::<Vec<isize>>()
        .try_into()
        .unwrap() // safe to unwrap as is size 3
}

/// compute the dot product between two vectors
pub fn vdot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b.iter()).map(|(a, b)| a * b).sum()
}

/// compute the norm of a vector
pub fn norm(a: [f64; 3]) -> f64 {
    a.iter().map(|a| a.powi(2)).sum::<f64>().powf(0.5)
}

/// compute the sum of two vectors
pub fn vsum(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// compute M.T * M
pub fn transpose_square(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [
            vdot([m[0][0], m[1][0], m[2][0]], [m[0][0], m[1][0], m[2][0]]),
            vdot([m[0][0], m[1][0], m[2][0]], [m[0][1], m[1][1], m[2][1]]),
            vdot([m[0][0], m[1][0], m[2][0]], [m[0][2], m[1][2], m[2][2]]),
        ],
        [
            vdot([m[0][1], m[1][1], m[2][1]], [m[0][0], m[1][0], m[2][0]]),
            vdot([m[0][1], m[1][1], m[2][1]], [m[0][1], m[1][1], m[2][1]]),
            vdot([m[0][1], m[1][1], m[2][1]], [m[0][2], m[1][2], m[2][2]]),
        ],
        [
            vdot([m[0][2], m[1][2], m[2][2]], [m[0][0], m[1][0], m[2][0]]),
            vdot([m[0][2], m[1][2], m[2][2]], [m[0][1], m[1][1], m[2][1]]),
            vdot([m[0][2], m[1][2], m[2][2]], [m[0][2], m[1][2], m[2][2]]),
        ],
    ]
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

/// returns the first index that is not vacuum from a sorted index list
pub fn vacuum_index(
    density: &[f64],
    index: &[usize],
    tolerance: Option<f64>,
) -> Result<usize, VacuumError> {
    match tolerance {
        Some(tol) => {
            for (i, p) in index.iter().rev().enumerate() {
                if density[*p] > tol {
                    return Ok(index.len() - i);
                }
            }
            Err(VacuumError {
                vacuum_tolerance: tol,
                density: density[index[0]],
            })
        }
        None => Ok(index.len()),
    }
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
    fn utils_transpose_square() {
        let matrix = [[3., 0., 0.], [2.5, 2., 0.], [0., 0., 5.]];
        let t_squared = [[15.25, 5., 0.], [5., 4., 0.], [0., 0., 25.]];
        assert_eq!(transpose_square(matrix), t_squared)
    }

    #[test]
    fn utils_vacuum_index_some_high() {
        let data = (0..60).map(|x| x as f64).collect::<Vec<f64>>();
        let index = (0..60).rev().collect::<Vec<usize>>();
        assert!(vacuum_index(&data, &index, Some(100.)).is_err())
    }

    #[test]
    fn utils_vacuum_index_some_low() {
        let data = (0..60).map(|x| x as f64).collect::<Vec<f64>>();
        let index = (0..60).rev().collect::<Vec<usize>>();
        let i = vacuum_index(&data, &index, Some(-1.)).unwrap();
        assert_eq!(i, 60)
    }

    #[test]
    fn utils_vacuum_index_some() {
        let data = (0..60).map(|x| x as f64).collect::<Vec<f64>>();
        let index = (0..60).rev().collect::<Vec<usize>>();
        let i = vacuum_index(&data, &index, Some(10.)).unwrap();
        assert_eq!(i, 49)
    }

    #[test]
    fn utils_vacuum_index_none() {
        let data = (0..60).map(|x| x as f64).collect::<Vec<f64>>();
        let index = (0..60).rev().collect::<Vec<usize>>();
        let i = vacuum_index(&data, &index, None).unwrap();
        assert_eq!(i, 60)
    }
}
