use crate::atoms::{Lattice, ReducedLattice};
use crate::utils::{cross, invert_lattice, vdot};

pub struct Voronoi {
    pub vectors: Vec<Vec<usize>>,
    pub alphas: Vec<f64>,
    pub lll_lattice: ReducedLattice,
}

impl Voronoi {
    pub fn new(lattice: &Lattice) -> Self {
        let lll_lattice = ReducedLattice::from_lattice(lattice);
        let (vectors, alphas) = Voronoi::voronoi_vectors(&lll_lattice);
        Self { vectors,
               alphas,
               lll_lattice }
    }

    fn voronoi_vectors(lll: &ReducedLattice) -> (Vec<Vec<usize>>, Vec<f64>) {
        // allocate the storage for voronoi vectors and flux coefficients
        let mut vectors = Vec::<Vec<usize>>::with_capacity(14);
        let mut alphas = Vec::<f64>::with_capacity(14);
        // allocate the vertex storage and vector/matrix for calculating them
        let mut vertices = Vec::<[f64; 3]>::with_capacity(28);
        let mut vector_basis = [[0f64; 3]; 3];
        let mut vector_mag = [0f64; 3];
        // allocate the plane vectors for each voronoi vector
        let mut rx = [0f64; 3];
        // calculate the reduced lattice
        'vector: for vec_i in 0..26 {
            let c_shift = lll.cartesian_shift_matrix[vec_i];
            vector_basis[0][..3].clone_from_slice(&c_shift[..3]);
            vector_mag[0] = vdot(c_shift, c_shift) * 0.5;
            for neigh_a in 0..26 {
                let c_neigh_a = lll.cartesian_shift_matrix[neigh_a];
                vector_basis[1][..3].clone_from_slice(&c_neigh_a[..3]);
                vector_mag[1] = vdot(c_neigh_a, c_neigh_a) * 0.5;
                'neigh_b: for neigh_b in (neigh_a + 1)..26 {
                    let c_neigh_b = lll.cartesian_shift_matrix[neigh_b];
                    vector_basis[2][..3].clone_from_slice(&c_neigh_b[..3]);
                    vector_mag[2] = vdot(c_neigh_b, c_neigh_b) * 0.5;
                    match invert_lattice(&vector_basis) {
                        Ok(vector_inv) => {
                            let mut vertex = [0f64; 3];
                            for i in 0..3 {
                                vertex[i] = vdot(vector_mag, vector_inv[i])
                            }
                            for cell_check in lll.cartesian_shift_matrix.iter()
                            {
                                let vector_2 =
                                    0.5 * vdot(*cell_check, *cell_check);
                                if vdot(vertex, *cell_check) > vector_2 + 1E-8 {
                                    continue 'neigh_b;
                                }
                            }
                            let vertex_mag = vdot(vertex, vertex);
                            if (vertex_mag - 0.25 * vdot(c_shift, c_shift)).abs()
                               < 1E-8
                            {
                                vertices.clear();
                                continue 'vector;
                            } else if vertex_mag > 0. {
                                vertices.push(vertex);
                            }
                        }
                        Err(_) => continue 'neigh_b,
                    }
                }
            }
            if vertices.is_empty() {
                continue 'vector;
            }
            rx[..3].clone_from_slice(&vertices[0][..3]);
            let r_coeff = vdot(rx, c_shift) / vdot(c_shift, c_shift);
            for (i, c_shift) in c_shift.iter().enumerate() {
                rx[i] -= c_shift * r_coeff;
            }
            let r_coeff = vdot(rx, rx).powf(-0.5);
            for rx in &mut rx {
                *rx *= r_coeff;
            }
            let mut ry = cross(c_shift, rx);
            let r_coeff = vdot(ry, ry).powf(-0.5);
            for ry in &mut ry {
                *ry *= r_coeff;
            }
            vertices.sort_unstable_by({
                        |a, b| {
                            let c = vdot(*a, ry).atan2(vdot(*a, rx));
                            let d = vdot(*b, ry).atan2(vdot(*b, rx));
                            c.partial_cmp(&d).unwrap()
                        }
                    });
            let num_vertices = vertices.len();
            let alpha = vertices.iter()
                                .enumerate()
                                .map(|(i, v)| {
                                    vdot(*v,
                                         cross(vertices
                                                   [(i + 1) % num_vertices],
                                               c_shift))
                                })
                                .sum::<f64>()
                        / (2. * vdot(c_shift, c_shift));
            if alpha.abs() < 1E-8 {
                vertices.clear();
                continue 'vector;
            }

            vectors.push(lll.shift_matrix[vec_i].clone());
            alphas.push(alpha);
            vertices.clear();
        }
        vectors.shrink_to_fit();
        alphas.shrink_to_fit();
        (vectors, alphas)
    }
}
/*
#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_voronoi_vectors() {
        let lattice = Lattice::new([[8.39866536000000 * 2. / 240., 0., 0.], [0., 8.39866536000000 * 2. / 240., 0.], [0., 0., 8.39866536000000 * 3.286996 / 384.]]);
        let vecs = Vec::<Vec<usize>>::with_capacity(14);
        assert_eq!(vecs, voronoi_vectors(&lattice))
    }
}
*/
