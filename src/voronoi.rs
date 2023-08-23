use crate::atoms::{Lattice, ReducedLattice};
use crate::utils::{cross, invert_lattice, vdot};

/// Holds information on the Voronoi vectors.
pub struct Voronoi {
    /// Voronoi vectors as indices in a shift matrix.
    pub vectors: Vec<Vec<usize>>,
    /// The alphas associated with each Voronoi vector.
    /// alphas are used to multiply the charge difference by to calculate flux.
    pub alphas: Vec<f64>,
    /// The volume of the Voronoi cell.
    pub volume: f64,
    /// The LLL-reduced lattice for the voxel basis.
    pub lll_lattice: ReducedLattice,
}

impl Voronoi {
    /// Generates a Voronoi struct from a [`Lattice`].
    pub fn new(lattice: &Lattice) -> Self {
        let lll_lattice = ReducedLattice::from_lattice(lattice);
        let (vectors, alphas, volume) = Voronoi::voronoi_vectors(&lll_lattice);
        Self { vectors,
               alphas,
               volume,
               lll_lattice }
    }

    /// Calculates the Voronoi vectors and their alphas from a reduced basis.
    /// As the volume is identical for each cell alpha is effectively the area of the facet divided
    /// by the length of the vector.
    fn voronoi_vectors(lll: &ReducedLattice)
                       -> (Vec<Vec<usize>>, Vec<f64>, f64) {
        // allocate the storage for voronoi vectors and flux coefficients
        let mut vectors = Vec::<Vec<usize>>::with_capacity(14);
        let mut alphas = Vec::<f64>::with_capacity(14);
        let mut volume = 0f64;
        // allocate the vertex storage and vector/matrix for calculating them
        let mut vertices = Vec::<[f64; 3]>::with_capacity(28);
        let mut vector_mag = [0f64; 3];
        // allocate the plane vectors for each voronoi vector
        let mut rx = [0f64; 3];
        // find the vertices for each plane by solving 3-way intersections between the plane and
        // every other plane pair and then checking if it falls within the Voronoi volume
        'vector: for vec_i in 0..26 {
            let n = lll.cartesian_shift_matrix[vec_i];
            vector_mag[0] = vdot(n, n) * 0.5;
            for neigh_a in 0..26 {
                let r1 = lll.cartesian_shift_matrix[neigh_a];
                vector_mag[1] = vdot(r1, r1) * 0.5;
                'neigh_b: for neigh_b in (neigh_a + 1)..26 {
                    let r2 = lll.cartesian_shift_matrix[neigh_b];
                    vector_mag[2] = vdot(r2, r2) * 0.5;
                    // If not invertable then no crossing point
                    if let Some(vector_inv) = invert_lattice(&[n, r1, r2]) {
                        let mut vertex = [0f64; 3];
                        for i in 0..3 {
                            vertex[i] = vdot(vector_mag, vector_inv[i])
                        }
                        // Check that the vertex isn't on the orginal vector which would make all
                        // vertices apart from this one outside the Voronoi volume
                        let vertex_mag = vdot(vertex, vertex);
                        if (vertex_mag - 0.25 * vdot(n, n)).abs() < f64::EPSILON
                        {
                            vertices.clear();
                            continue 'vector;
                        } else if vertex_mag < f64::EPSILON {
                            continue 'neigh_b;
                        }
                        // is this vertex inside the Voronoi volume, project along every lll vector
                        // and compare to vector * vector / 2 if higher then outside Voronoi volume
                        for s in lll.cartesian_shift_matrix.iter() {
                            let ss2 = 0.5 * vdot(*s, *s);
                            if vdot(vertex, *s) > ss2 + f64::EPSILON {
                                continue 'neigh_b;
                            }
                        }
                        vertices.push(vertex);
                    }
                }
            }
            // if the current vector does not form the Voronoi bounding planes then skip
            if vertices.is_empty() {
                continue 'vector;
            }
            // order the vertices by projecting on to an orthogonal basis that is itself orthogonal
            // to n for calculating the area
            rx[..3].clone_from_slice(&vertices[0][..3]);
            let r_coeff = vdot(rx, n) / vdot(n, n);
            for (i, n) in n.iter().enumerate() {
                rx[i] -= n * r_coeff;
            }
            let r_coeff = vdot(rx, rx).powf(-0.5);
            for rx in &mut rx {
                *rx *= r_coeff;
            }
            let mut ry = cross(n, rx);
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
            // calculate the area of the facet and divide by the length
            // first calculate the volume of the tetrahedron of v[i], v[i+1] & n/2. This is the
            // scalar tripple product of the vectors divided by 6. The volume is for the Laplacian
            // calculation later
            let wedge_volume =
                vertices.iter()
                        .enumerate()
                        .map(|(i, v)| {
                            vdot(*v, cross(vertices[(i + 1) % num_vertices], n))
                        })
                        .sum::<f64>()
                / 12f64;
            if wedge_volume.abs() < f64::EPSILON {
                vertices.clear();
                continue 'vector;
            }
            // now we need to turn the volume into an area and divide by |n|
            // V = Ah/3 where h = n/2
            let alpha = 6f64 * wedge_volume / vdot(n, n);
            vectors.push(lll.shift_matrix[vec_i].clone());
            alphas.push(alpha);
            volume += wedge_volume;
            vertices.clear();
        }
        vectors.shrink_to_fit();
        alphas.shrink_to_fit();
        (vectors, alphas, volume)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_voronoi_vectors() {
        let lattice =
            Lattice::new([[1.0, 0., 0.], [0.707, 0.707, 0.], [0., 0., 1.]]);
        let voronoi = Voronoi::new(&lattice);
        let vecs = vec![vec![19],
                        vec![22],
                        vec![10],
                        vec![12],
                        vec![14],
                        vec![16],
                        vec![4],
                        vec![7],];
        assert_eq!(vecs, voronoi.vectors)
    }
}
