/// compute the dot product between a vector and a matrix
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1],
     a[2] * b[0] - a[0] * b[2],
     a[0] * b[1] - a[1] * b[0]]
}
/// compute the dot product between a vector and a matrix
pub fn dot(v: [f64; 3], m: [[f64; 3]; 3]) -> [f64; 3] {
    let mut out = [0f64; 3];
    for i in 0..3 {
        out[i] = v[0] * m[0][i] + v[1] * m[1][i] + v[2] * m[2][i]
    }
    return out;
}

/// compute the dot product between two vectors
pub fn vdot(a: [f64; 3], b: [f64; 3]) -> f64 {
    let mut out = 0f64;
    for i in 0..3 {
        out += a[i] * b[i]
    }
    return out;
}

/// compute the norm of a vector
pub fn norm(a: [f64; 3]) -> f64 {
    a.iter().map(|a| a.powi(2)).sum::<f64>().powf(0.5)
}

/// compute M.T * M
pub fn transpose_square(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [[vdot([m[0][0], m[1][0], m[2][0]], [m[0][0], m[1][0], m[2][0]]),
      vdot([m[0][0], m[1][0], m[2][0]], [m[0][1], m[1][1], m[2][1]]),
      vdot([m[0][0], m[1][0], m[2][0]], [m[0][2], m[1][2], m[2][2]])],
     [vdot([m[0][1], m[1][1], m[2][1]], [m[0][0], m[1][0], m[2][0]]),
      vdot([m[0][1], m[1][1], m[2][1]], [m[0][1], m[1][1], m[2][1]]),
      vdot([m[0][1], m[1][1], m[2][1]], [m[0][2], m[1][2], m[2][2]])],
     [vdot([m[0][2], m[1][2], m[2][2]], [m[0][0], m[1][0], m[2][0]]),
      vdot([m[0][2], m[1][2], m[2][2]], [m[0][1], m[1][1], m[2][1]]),
      vdot([m[0][2], m[1][2], m[2][2]], [m[0][2], m[1][2], m[2][2]])]]
}

pub fn invert_lattice(lattice: &[[f64; 3]; 3])
                      -> Result<[[f64; 3]; 3], String> {
    let minor00 = lattice[1][1] * lattice[2][2] - lattice[1][2] * lattice[2][1];
    let minor01 = lattice[1][0] * lattice[2][2] - lattice[1][2] * lattice[2][0];
    let minor02 = lattice[1][0] * lattice[2][1] - lattice[1][1] * lattice[2][0];
    let determinant = lattice[0][0] * minor00 - lattice[0][1] * minor01
                      + lattice[0][2] * minor02;
    if determinant.abs() < 1e-16 {
        Err(String::from("Lattice doesn't span 3D space"))
    } else {
        Ok([[minor00 / determinant,
             (lattice[0][2] * lattice[2][1] - lattice[2][2] * lattice[0][1])
             / determinant,
             (lattice[0][1] * lattice[1][2] - lattice[1][1] * lattice[0][2])
             / determinant],
            [-minor01 / determinant,
             (lattice[0][0] * lattice[2][2] - lattice[2][0] * lattice[0][2])
             / determinant,
             (lattice[0][2] * lattice[1][0] - lattice[1][2] * lattice[0][0])
             / determinant],
            [minor02 / determinant,
             (lattice[0][1] * lattice[2][0] - lattice[2][1] * lattice[0][0])
             / determinant,
             (lattice[0][0] * lattice[1][1] - lattice[1][0] * lattice[0][1])
             / determinant]])
    }
}

pub struct BTMap {
    map: Vec<(Vec<(isize, f64)>, usize)>,
}

impl BTMap {
    pub fn new(size: usize) -> Self {
        Self { map: vec![(vec![], 0); size] }
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    pub fn get(&mut self, key: isize) -> Option<Vec<(isize, f64)>> {
        let (value, count) = match self.map.get(key as usize) {
            Some(x) => x.clone(),
            None => (vec![], 0),
        };
        if count == 0 {
            return None;
        } else {
            self.insert(key, value.clone(), count - 1);
        }
        Some(value)
    }

    pub fn get_sly(&self, key: isize) -> Option<Vec<(isize, f64)>> {
        let (value, count) = match self.map.get(key as usize) {
            Some(x) => x.clone(),
            None => (vec![], 0),
        };
        if count == 0 {
            return None;
        }
        Some(value)
    }

    pub fn insert(&mut self,
                  key: isize,
                  value: Vec<(isize, f64)>,
                  count: usize) {
        if count == 0 {
            self.map[key as usize] = (vec![], 0);
        } else {
            self.map[key as usize] = (value, count);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn utils_dot() {
        assert_eq!(dot([1., 2., 3.],
                       [[1., 0., 0.], [0., 2., 0.], [0., 0., 3.]]),
                   [1., 4., 9.])
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
}
