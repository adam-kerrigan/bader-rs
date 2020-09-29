use crate::density::Density;
use crate::voxel_map::VoxelMap;
use crate::utils;

/// struct for containing the information about the atoms
///
/// > lattice: Lattice - the lattice of the structure
/// > positions: Vec<[f64; 3]> - the positions of the atoms in cartesian coordinates
/// > text: String - text representation from the input file
pub struct Atoms {
    pub lattice: Lattice,
    pub positions: Vec<[f64; 3]>,
    pub text: String,
    pub reduced_lattice: ReducedLattice,
    pub reduced_positions: Vec<[f64; 3]>,
}

impl Atoms {
    /// initialises the structure
    pub fn new(lattice: Lattice,
               positions: Vec<[f64; 3]>,
               text: String)
               -> Self {
        let reduced_lattice = ReducedLattice::from_lattice(&lattice);
        let reduced_positions =
            positions.iter()
                     .map(|v| {
                         let mut frac =
                             utils::dot(*v, reduced_lattice.to_fractional);
                         for f in &mut frac {
                             *f = f.rem_euclid(1.);
                         }
                         utils::dot(frac, reduced_lattice.to_cartesian)
                     })
                     .collect::<Vec<[f64; 3]>>();
        Self { lattice,
               positions,
               reduced_lattice,
               reduced_positions,
               text }
    }

    /// assigns bader volumes to their nearest atom
    pub fn assign_maxima(&self,
                         map: &VoxelMap,
                         density: &Density)
                         -> (Vec<usize>, Vec<f64>) {
        let max_distance =
            (self.lattice.a + self.lattice.b + self.lattice.c).powi(2);
        let mut assigned_atom =
            Vec::<usize>::with_capacity(map.bader_maxima.len());
        let mut assigned_distance =
            Vec::<f64>::with_capacity(map.bader_maxima.len());
        let maximas = if map.bader_maxima[0] == -1 {
            &map.bader_maxima[1..]
        } else {
            &map.bader_maxima[..]
        };
        for maxima in maximas.iter() {
            let bx = density.voxel_origin[0]
                     + (maxima / (density.size.y * density.size.z)) as f64;
            let by = density.voxel_origin[1]
                     + (maxima / density.size.z).rem_euclid(density.size.y)
                       as f64;
            let bz = density.voxel_origin[2]
                     + maxima.rem_euclid(density.size.z) as f64;
            let maxima_cartesian =
                utils::dot([bx, by, bz], density.voxel_lattice.to_cartesian);
            let mut maxima_lll_fractional = utils::dot(maxima_cartesian,
                                                       self.reduced_lattice
                                                           .to_fractional);
            for f in &mut maxima_lll_fractional {
                *f = f.rem_euclid(1.);
            }
            let maxima_lll_cartesian = utils::dot(maxima_lll_fractional,
                                                  self.reduced_lattice
                                                      .to_cartesian);
            let mut atom_num = 0;
            let mut min_distance = max_distance;
            for (i, atom) in self.reduced_positions.iter().enumerate() {
                for atom_shift in
                    self.reduced_lattice.cartesian_shift_matrix.iter()
                {
                    let distance = {
                        (maxima_lll_cartesian[0] - (atom[0] + atom_shift[0])).powi(2)
                        + (maxima_lll_cartesian[1] - (atom[1] + atom_shift[1])).powi(2)
                        + (maxima_lll_cartesian[2] - (atom[2] + atom_shift[2])).powi(2)
                    };
                    if distance < min_distance {
                        min_distance = distance;
                        atom_num = i;
                    }
                }
            }
            assigned_atom.push(atom_num);
            assigned_distance.push(min_distance.powf(0.5));
        }
        (assigned_atom, assigned_distance)
    }
}

/// Lattice - structure for containing information on the cell
///
/// > a: f64 - length of the a-vector
/// > b: f64 - length of the b-vector
/// > c: f64 - length of the c-vector
/// > distance_matrix: [f64; 26] - distance to move in each direction, follows the
/// >                              ordering listed below
/// > gradient_transform: [[f64; 3]; 3] - transformation matrix for converting central
/// >                                      difference gradients to the lattice basis
/// > to_fractional: [[f64; 3]; 3] - transformation matrix for converting to fractional
/// >                                coordinates
/// > to_cartesian: [[f64; 3]; 3] - transformation matrix for converting to cartesian
/// >                               coordinates
///
/// > distance_matrix ordering:
/// >     0 -> (-1,-1,-1)   7 -> (-1, 1, 0)  14 -> (0, 1,-1)  21 -> (1, 0, 0)
/// >     1 -> (-1,-1, 0)   8 -> (-1, 1, 1)  15 -> (0, 1, 0)  22 -> (1, 0, 1)
/// >     2 -> (-1,-1, 1)   9 -> (0,-1,-1)   16 -> (0, 1, 1)  23 -> (1, 1,-1)
/// >     3 -> (-1, 0,-1)  10 -> (0,-1, 0)   17 -> (1,-1,-1)  24 -> (1, 1, 0)
/// >     4 -> (-1, 0, 0)  11 -> (0,-1, 1)   18 -> (1,-1, 0)  25 -> (1, 1, 1)
/// >     5 -> (-1, 0, 1)  12 -> (0, 0,-1)   19 -> (1,-1, 1)
/// >     6 -> (-1, 1,-1)  13 -> (0, 0, 1)   20 -> (1, 0,-1)
pub struct Lattice {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub distance_matrix: [f64; 26],
    pub shift_matrix: [[f64; 3]; 27],
    pub gradient_transform: [[f64; 3]; 3],
    pub to_fractional: [[f64; 3]; 3],
    pub to_cartesian: [[f64; 3]; 3],
    pub volume: f64,
}

impl Lattice {
    /// Initialises the structure. Builds all the fields of the lattice structure
    /// from a 2d vector in the form:
    ///
    /// > [
    /// >     [ax, ay, az],
    /// >     [bx, by, bz],
    /// >     [cx, cy, cz],
    /// >  ]
    pub fn new(lattice: [[f64; 3]; 3]) -> Self {
        let a = utils::norm(lattice[0]);
        let b = utils::norm(lattice[1]);
        let c = utils::norm(lattice[2]);
        let x = lattice[0];
        let y = lattice[1];
        let z = lattice[2];
        let shift_matrix =
            [[-x[0] - y[0] - z[0],
              -x[1] - y[1] - z[1],
              -x[2] - y[2] - z[2]],
             [-x[0] - y[0], -x[1] - y[1], -x[2] - y[2]],
             [-x[0] - y[0] + z[0],
              -x[1] - y[1] + z[1],
              -x[2] - y[2] + z[2]],
             [-x[0] - z[0], -x[1] - z[1], -x[2] - z[2]],
             [-x[0], -x[1], -x[2]],
             [-x[0] + z[0], -x[1] + z[1], -x[2] + z[2]],
             [-x[0] + y[0] - z[0],
              -x[1] + y[1] - z[1],
              -x[2] + y[2] - z[2]],
             [-x[0] + y[0], -x[1] + y[1], -x[2] + y[2]],
             [-x[0] + y[0] + z[0],
              -x[1] + y[1] + z[1],
              -x[2] + y[2] + z[2]],
             [-y[0] - z[0], -y[1] - z[1], -y[2] - z[2]],
             [-y[0], -y[1], -y[2]],
             [-y[0] + z[0], -y[1] + z[1], -y[2] + z[2]],
             [-z[0], -z[1], -z[2]],
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
             [0., 0., 0.]];
        let mut distance_matrix = [0f64; 26];
        for i in 0..26 {
            distance_matrix[i] = utils::norm(shift_matrix[i]);
        }
        let to_fractional = match utils::invert_lattice(&lattice) {
            Ok(inv) => inv,
            Err(e) => panic!("{}", e),
        };
        let to_cartesian = lattice;
        let volume = {
            (lattice[0][0]
             * (lattice[1][1] * lattice[2][2] - lattice[1][2] * lattice[2][1])
             + lattice[0][1]
               * (lattice[1][2] * lattice[2][0]
                  - lattice[1][0] * lattice[2][2])
             + lattice[0][2]
               * (lattice[1][0] * lattice[2][1]
                  - lattice[1][1] * lattice[2][0]))
                                                   .abs()
        };
        let gradient_transform = utils::transpose_square(to_fractional);
        Self { a,
               b,
               c,
               distance_matrix,
               shift_matrix,
               gradient_transform,
               to_fractional,
               to_cartesian,
               volume }
    }
}

pub struct ReducedLattice {
    pub shift_matrix: Vec<Vec<usize>>,
    pub cartesian_shift_matrix: [[f64; 3]; 27],
    pub distance_matrix: [f64; 26],
    pub to_cartesian: [[f64; 3]; 3],
    pub to_fractional: [[f64; 3]; 3],
}

impl ReducedLattice {
    pub fn from_lattice(lattice: &Lattice) -> Self {
        let reduced_lattice =
            Lattice::new(ReducedLattice::lll_lattice(lattice.to_cartesian));
        let to_cartesian = reduced_lattice.to_cartesian;
        let to_fractional = reduced_lattice.to_fractional;
        let distance_matrix = reduced_lattice.distance_matrix;
        let cartesian_shift_matrix = reduced_lattice.shift_matrix;
        let mut shift_matrix = Vec::with_capacity(26);
        for c_shift in cartesian_shift_matrix.iter().take(26) {
            shift_matrix.push({
                            let mut shift_vec = Vec::<usize>::new();
                            let mut shift =
                    utils::dot(*c_shift, lattice.to_fractional).iter()
                                                     .map(|x| {
                                                         ((x * 1E14).round()
                                                          / 1E14)
                                                         as isize
                                                     })
                                                     .collect::<Vec<isize>>();
                            let max_shift =
                                shift.iter().map(|x| x.abs()).max().unwrap();
                            for _ in 0..max_shift {
                                let mut out = [0isize; 3];
                                for k in 0..3 {
                                    if shift[k].abs() > 1 {
                                        out[k] = shift[k] / shift[k].abs();
                                    } else {
                                        out[k] = shift[k];
                                    }
                                    shift[k] -= out[k];
                                }
                                shift_vec.push((out[0] * 9
                                                + out[1] * 3
                                                + out[2]
                                                + 13)
                                               as usize);
                            }
                            shift_vec
                        });
        }
        Self { shift_matrix,
               cartesian_shift_matrix,
               distance_matrix,
               to_cartesian,
               to_fractional }
    }

    fn lll_lattice(lattice: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let delta = 0.75;
        let mut a = lattice;
        let (mut b, mut mu) = ReducedLattice::gram_schmidt(&a);
        let mut i = 1usize;
        while i <= 2 {
            for j in (0..i).rev() {
                match mu[i][j] {
                    q if q.abs() <= 0.5 => (),
                    q => {
                        for k in 0..3 {
                            a[i][k] -= q.round() * a[j][k];
                        }
                        let (b_temp, mu_temp) =
                            ReducedLattice::gram_schmidt(&a);
                        b = b_temp;
                        mu = mu_temp;
                    }
                }
            }
            if utils::vdot(b[i], b[i])
               >= (delta - mu[i][i - 1].powi(2))
                  * utils::vdot(b[i - 1], b[i - 1])
            {
                i += 1;
            } else {
                for j in 0..3 {
                    b[0][0] = a[i][j];
                    a[i][j] = a[i - 1][j];
                    a[i - 1][j] = b[0][0];
                }
                let (b_temp, mu_temp) = ReducedLattice::gram_schmidt(&a);
                b = b_temp;
                mu = mu_temp;
                i = 1usize.max(i - 1);
            }
        }
        a
    }

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
    fn atoms_assign_maxima() {
        let positions = vec![[0.; 3]];
        let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
        let text = String::new();
        let atoms = Atoms::new(lattice, positions, text);
        let map = vec![63; 64];
        let bader_maxima = vec![63];
        let voxel_map = VoxelMap::new(map, bader_maxima);
        let data = (0..64).map(|x| x as f64).collect::<Vec<f64>>();
        let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
        let density = Density::new(&data,
                                   [4, 4, 4],
                                   lattice.to_cartesian,
                                   None,
                                   [0., 0., 0.]);
        let (atom, distance) = atoms.assign_maxima(&voxel_map, &density);
        assert_eq!(atom, vec![0]);
        assert_eq!(distance[0], density.voxel_lattice.distance_matrix[0])
    }

    #[test]
    fn atoms_assign_maxima_vacuum() {
        let positions = vec![[0.; 3]];
        let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
        let text = String::new();
        let atoms = Atoms::new(lattice, positions, text);
        let map = vec![63; 64];
        let bader_maxima = vec![-1, 63];
        let voxel_map = VoxelMap::new(map, bader_maxima);
        let data = (0..64).map(|x| x as f64).collect::<Vec<f64>>();
        let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
        let density = Density::new(&data,
                                   [4, 4, 4],
                                   lattice.to_cartesian,
                                   None,
                                   [0., 0., 0.]);
        let (atom, distance) = atoms.assign_maxima(&voxel_map, &density);
        assert_eq!(atom, vec![0]);
        assert_eq!(distance[0], density.voxel_lattice.distance_matrix[0])
    }

    #[test]
    fn lattice_new() {
        let lattice = Lattice::new([[1., 0., 0.], [0., 2., 0.], [0., 0., 2.]]);
        let distance_matrix = [3.,
                               5f64.powf(0.5),
                               3.,
                               5f64.powf(0.5),
                               1.,
                               5f64.powf(0.5),
                               3.,
                               5f64.powf(0.5),
                               3.,
                               8f64.powf(0.5),
                               2.,
                               8f64.powf(0.5),
                               2.,
                               2.,
                               8f64.powf(0.5),
                               2.,
                               8f64.powf(0.5),
                               3.,
                               5f64.powf(0.5),
                               3.,
                               5f64.powf(0.5),
                               1.,
                               5f64.powf(0.5),
                               3.,
                               5f64.powf(0.5),
                               3.];
        assert_eq!(lattice.distance_matrix, distance_matrix)
    }

    #[test]
    #[should_panic]
    fn lattice_new_non_invert() {
        let _ = Lattice::new([[1., 0., 0.], [1., 0., 0.], [0., 0., 2.]]);
    }

    #[test]
    fn reduced_lattice_lll_lattice() {
        let lat = [[0., 1., 0.], [1., 0., 1.], [-1., 0., 2.]];
        assert_eq!(lat,
                   ReducedLattice::lll_lattice([[1., 1., 1.],
                                                [-1., 0., 2.],
                                                [3., 5., 6.]]));
    }
}
