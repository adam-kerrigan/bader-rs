use crate::density::{Density, VoxelMap};
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
}

impl Atoms {
    /// initialises the structure
    pub fn new(lattice: Lattice,
               positions: Vec<[f64; 3]>,
               text: String)
               -> Self {
        return Self { lattice,
                      positions,
                      text };
    }

    /// assigns bader volumes to their nearest atom
    pub fn assign_maxima(&self,
                         map: &VoxelMap,
                         density: &Density)
                         -> (Vec<usize>, Vec<f64>) {
        let max_distance = self.lattice.distance_matrix[0].powi(2);
        let mut assigned_atom =
            Vec::<usize>::with_capacity(map.bader_maxima.len());
        let mut assigned_distance =
            Vec::<f64>::with_capacity(map.bader_maxima.len());
        let maximas = match density.vacuum_tolerance {
            Some(_) => &map.bader_maxima[1..],
            None => &map.bader_maxima[..],
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
            let mut atom_num = 0;
            let mut min_distance = max_distance;
            for (i, atom) in self.positions.iter().enumerate() {
                for atom_shift in self.lattice.shift_matrix.iter() {
                    let distance = {
                        (maxima_cartesian[0] - (atom[0] + atom_shift[0])).powi(2)
                        + (maxima_cartesian[1] - (atom[1] + atom_shift[1])).powi(2)
                        + (maxima_cartesian[2] - (atom[2] + atom_shift[2])).powi(2)
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
        return (assigned_atom, assigned_distance);
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
        let to_fractional = {
            let minor00 =
                lattice[1][1] * lattice[2][2] - lattice[1][2] * lattice[2][1];
            let minor01 =
                lattice[1][0] * lattice[2][2] - lattice[1][2] * lattice[2][0];
            let minor02 =
                lattice[1][0] * lattice[2][1] - lattice[1][1] * lattice[2][0];
            let determinant = lattice[0][0] * minor00 - lattice[0][1] * minor01
                              + lattice[0][2] * minor02;
            [[minor00 / determinant,
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
              / determinant]]
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
        return Self { a,
                      b,
                      c,
                      distance_matrix,
                      shift_matrix,
                      gradient_transform,
                      to_fractional,
                      to_cartesian,
                      volume };
    }
}
