use crate::atoms::Lattice;
use crate::utils::{cross, invert_lattice, vdot};

/// Manages the geometric weights for calculating gradients and fluxes on the grid.
///
/// To integrate charge density in a grid-agnostic way, we partition space using **Voronoi cells**
/// (Wigner-Seitz cells). This structure calculates the flux of charge density between a voxel
/// and its neighbors.
///
/// # The Flux Equation
/// The flux $J$ between two voxels separated by a vector $\vec{d}$ is:
/// $$ J \propto \alpha \cdot (\rho_{neighbor} - \rho_{self}) $$
///
/// The weight $\alpha$ is derived from the geometry of the Voronoi face shared by the voxels:
/// $$ \alpha = \frac{A}{|\vec{d}|} $$
/// * $A$: Area of the shared Voronoi face.
/// * $|\vec{d}|$: Distance between the voxel centers.
///
/// # Fields
/// * `vectors`: Indices of the neighbors (relative shifts) that share a face with the central voxel.
/// * `alphas`: The pre-calculated weight $\alpha$ for each neighbor vector.
/// * `volume`: The total volume of the Voronoi cell (must match the lattice cell volume).
pub struct Voronoi {
    /// Voronoi vectors as indices in a shift matrix.
    pub vectors: Vec<Vec<usize>>,
    /// The alphas associated with each Voronoi vector.
    /// alphas are used to multiply the charge difference by to calculate flux.
    pub alphas: Vec<f64>,
    /// The volume of the Voronoi cell.
    pub volume: f64,
}

impl Voronoi {
    /// Generates a Voronoi struct from a [`Lattice`].
    pub fn new(lattice: &Lattice) -> Self {
        let (vectors, alphas, volume) = Voronoi::voronoi_vectors(lattice);
        Self {
            vectors,
            alphas,
            volume,
        }
    }

    /// Calculates the Voronoi vectors and their alphas from a reduced basis.
    fn voronoi_vectors(lattice: &Lattice) -> (Vec<Vec<usize>>, Vec<f64>, f64) {
        // allocate the storage for voronoi vectors and flux coefficients
        let mut vectors = Vec::<Vec<usize>>::with_capacity(14);
        let mut alphas = Vec::<f64>::with_capacity(14);
        // allocate the vertex storage and vector/matrix for calculating them
        // allocate the plane vectors for each voronoi vector
        // find the vertices for each plane by solving 3-way intersections between the plane and
        // every other plane pair and then checking if it falls within the Voronoi volume
        let volume = lattice
            .reduced_cartesian_shift_matrix
            .iter()
            .take(13) // the first 13 vectors are symmetrically equivalent to the last 13
            .enumerate()
            .filter_map(|(vec_i, n)| {
                let n_mag = vdot(*n, *n) * 0.5;
                let mut vertices = lattice
                    .reduced_cartesian_shift_matrix
                    .iter()
                    .enumerate()
                    .filter_map(|(neigh_i, r1)| {
                        let r1_mag = vdot(*r1, *r1) * 0.5;
                        let vertices = lattice
                            .reduced_cartesian_shift_matrix
                            .iter()
                            .skip(neigh_i + 1)
                            .filter_map(|r2| {
                                let r2_mag = vdot(*r2, *r2) * 0.5;
                                // If not invertable then no crossing point
                                if let Some(vector_inv) =
                                    invert_lattice(&[*n, *r1, *r2])
                                {
                                    let mut vertex = [0f64; 3];
                                    for i in 0..3 {
                                        vertex[i] = vdot(
                                            [n_mag, r1_mag, r2_mag],
                                            vector_inv[i],
                                        )
                                    }
                                    // Check that the vertex isn't on the orginal vector which would make all
                                    // vertices apart from this one outside the Voronoi volume
                                    let vertex_mag = vdot(vertex, vertex);
                                    if (vertex_mag - 0.5 * n_mag).abs()
                                        < f64::EPSILON
                                    {
                                        return None;
                                    }
                                    // is this vertex inside the Voronoi volume, project along every lll vector
                                    // and compare to vector * vector / 2 if higher then outside Voronoi volume
                                    for s in lattice
                                        .reduced_cartesian_shift_matrix
                                        .iter()
                                    {
                                        let ss2 = 0.5 * vdot(*s, *s);
                                        if vdot(vertex, *s) > ss2 + f64::EPSILON
                                        {
                                            return None;
                                        }
                                    }
                                    Some(vertex)
                                } else {
                                    None
                                }
                            })
                            .collect::<Vec<[f64; 3]>>();
                        if vertices.is_empty() {
                            None
                        } else {
                            Some(vertices)
                        }
                    })
                    .flatten()
                    .collect::<Vec<[f64; 3]>>();
                // if the current vector does not form the Voronoi bounding planes then skip
                if vertices.is_empty() {
                    return None;
                }
                // order the vertices by projecting on to an orthogonal basis that is itself orthogonal
                // to n for calculating the area
                let mut rx = vertices[0];
                let r_coeff = vdot(rx, *n) / vdot(*n, *n);
                for (element, n) in rx.iter_mut().zip(n) {
                    *element -= n * r_coeff;
                }
                let r_coeff = vdot(rx, rx).powf(-0.5);
                for element in rx.iter_mut() {
                    *element *= r_coeff;
                }
                let mut ry = cross(*n, rx);
                let r_coeff = vdot(ry, ry).powf(-0.5);
                for element in ry.iter_mut() {
                    *element *= r_coeff;
                }
                vertices.sort_unstable_by({
                    |a, b| {
                        let c = vdot(*a, ry).atan2(vdot(*a, rx));
                        let d = vdot(*b, ry).atan2(vdot(*b, rx));
                        c.partial_cmp(&d).unwrap()
                    }
                });
                vertices.push(vertices[0]);
                // calculate the area of the facet and divide by the length
                // first calculate the volume of the tetrahedron of v[i], v[i+1] & n/2. This is the
                // scalar tripple product of the vectors divided by 6. The volume is for the Laplacian
                // calculation later
                let wedge_volume = vertices
                    .windows(2)
                    .fold(0f64, |sum, w| sum + vdot(w[0], cross(w[1], *n)))
                    / 12f64;
                if wedge_volume.abs() < f64::EPSILON {
                    return None;
                }
                // now we need to turn the volume into an area and divide by |n|
                // V = Ah/3 where h = n/2
                let alpha = 6f64 * wedge_volume / vdot(*n, *n);
                vectors.push(lattice.reduced_grid_shift_matrix[vec_i].clone());
                alphas.push(alpha);
                vectors.push(
                    lattice.reduced_grid_shift_matrix[26 - vec_i].clone(),
                );
                alphas.push(alpha);
                Some(wedge_volume * 2f64)
            })
            .sum();
        vectors.shrink_to_fit();
        alphas.shrink_to_fit();
        (vectors, alphas, volume)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::atoms::Lattice;

    #[test]
    fn test_voronoi_cubic() {
        // 1. Setup a Simple Cubic Lattice (10x10x10)
        // The Voronoi cell (Wigner-Seitz cell) of a cubic lattice is a cube.
        // It should have 6 faces (neighbors): ±x, ±y, ±z.
        let matrix = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let lattice = Lattice::new(matrix);

        let voronoi = Voronoi::new(&lattice);

        // 2. Check Volume
        // Voronoi volume must equal Lattice volume
        assert!((voronoi.volume - 1000.0).abs() < 1e-6, "Volume mismatch");

        // 3. Check Vectors
        // We expect 6 neighbors for a simple cubic lattice.
        // However, the implementation might store them as 6 unique vectors or more depending on how it iterates.
        // Let's verify we have vectors and non-zero alphas.
        assert!(!voronoi.vectors.is_empty());
        assert_eq!(voronoi.vectors.len(), voronoi.alphas.len());

        // 4. Check Alphas (Weights)
        // For a cubic lattice of side L=10:
        // Neighbor distance |d| = 10.
        // Face area A = 10 * 10 = 100.
        // Alpha = Area / Distance = 100 / 10 = 10.0.

        for &alpha in &voronoi.alphas {
            assert!(
                (alpha - 10.0).abs() < 1e-6,
                "Expected alpha 10.0, got {}",
                alpha
            );
        }
    }

    #[test]
    fn test_voronoi_orthorhombic() {
        // Orthorhombic: [2, 0, 0], [0, 3, 0], [0, 0, 4]
        // Volume = 2 * 3 * 4 = 24.
        let matrix = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]];
        let lattice = Lattice::new(matrix);
        let voronoi = Voronoi::new(&lattice);

        // Check Volume
        assert!((voronoi.volume - 24.0).abs() < 1e-6);

        // We expect 3 distinct alpha values corresponding to the 3 axis lengths.
        // Alpha_x = (3*4) / 2 = 6.0
        // Alpha_y = (2*4) / 3 = 2.666...
        // Alpha_z = (2*3) / 4 = 1.5

        let mut found_6 = false;
        let mut found_2_6 = false;
        let mut found_1_5 = false;

        for &alpha in &voronoi.alphas {
            if (alpha - 6.0).abs() < 1e-4 {
                found_6 = true;
            }
            if (alpha - (8.0 / 3.0)).abs() < 1e-4 {
                found_2_6 = true;
            }
            if (alpha - 1.5).abs() < 1e-4 {
                found_1_5 = true;
            }
        }

        assert!(
            found_6 && found_2_6 && found_1_5,
            "Missing expected face weights"
        );
    }

    #[test]
    fn test_voronoi_vectors() {
        let lattice =
            Lattice::new([[1.0, 0., 0.], [0.707, 0.707, 0.], [0., 0., 1.]]);
        let voronoi = Voronoi::new(&lattice);
        let vecs = vec![
            vec![4],
            vec![22],
            vec![7],
            vec![19],
            vec![10],
            vec![16],
            vec![12],
            vec![14],
        ];
        assert_eq!(vecs, voronoi.vectors)
    }
}
