use crate::density::Density;
use crate::utils;
use crate::voxel_map::{Voxel, VoxelMap};
use std::collections::HashMap;
use std::sync::Arc;

/// A type for the maxima search function
pub type WeightMap = Arc<VoxelMap>;
pub type StepMethod =
    fn(usize, &Density, WeightMap) -> (isize, Vec<(usize, f64)>);

/// Steps in the density grid, from point p, following the gradient
pub fn ongrid_step(p: isize, density: &Density) -> isize {
    let mut pn = p;
    let mut pt = p;
    let control = density[pn];
    let mut max_val = density[pt];
    // colllect the shift and distances and iterate over them
    for (shift, d) in density.full_shift(p)
                             .iter()
                             .zip(&density.voxel_lattice.distance_matrix)
    {
        pt = p + shift;
        // calculate the gradient and if it is greater that the current max
        let rho = control + (density[pt] - control) / d;
        if rho > max_val {
            max_val = rho;
            pn = pt;
        }
    }
    pn
}

/// Finds the maxima associated with the current point, p.
pub fn ongrid(p: usize,
              density: &Density,
              weight_map: WeightMap)
              -> (isize, Vec<(usize, f64)>) {
    let mut pt = p as isize;
    let maxima = {
        pt = ongrid_step(pt, density);
        if pt == p as isize {
            pt
        } else {
            weight_map.maxima_get(pt)
        }
    };
    (maxima, Vec::with_capacity(0))
}

/// Step in the density grid following the gradient.
///
/// Use the current point, p, and step along density and use the residual
/// gradient, dr, to account for grid bias.
pub fn neargrid_step(p: isize, dr: &mut [f64; 3], density: &Density) -> isize {
    let mut den_t = [0f64; 2];
    let mut grad = [0f64; 3];
    let mut max_grad = -1f64;
    let mut pt = p;
    // retrive the [+x, -x, +y, -y, +z, -z] shifts
    let p_shift = density.reduced_shift(p);
    for i in 0..3 {
        den_t[0] = density[p + p_shift[i * 2]];
        den_t[1] = density[p + p_shift[i * 2 + 1]];
        // if not a maxima calculate the density
        if (density[p] < den_t[0]) | (density[p] < den_t[1]) {
            grad[i] = (den_t[0] - den_t[1]) / 2f64;
        }
    }
    // convert the gradient to the lattice basis and find the max value
    grad = utils::dot(grad, density.voxel_lattice.gradient_transform);
    for g in &grad {
        if g.abs() > max_grad {
            max_grad = g.abs();
        }
    }
    // if we arent at a maxima make a step in the gradient
    if max_grad > 0f64 {
        let mut dr_i = [0f64; 3];
        for i in 0..3 {
            grad[i] /= max_grad;
            dr[i] += grad[i] - grad[i].round();
            grad[i] = grad[i].round();
            dr_i[i] = dr[i].round();
            dr[i] -= dr[i].round();
        }
        pt += density.gradient_shift(pt, grad);
        pt += density.gradient_shift(pt, dr_i);
    }
    pt
}

/// Find the maxima assiociated with the current point using the neargrid method
pub fn neargrid(p: usize,
                density: &Density,
                _weight_map: WeightMap)
                -> (isize, Vec<(usize, f64)>) {
    let mut pt = p as isize;
    if let Some(x) = density.vacuum_tolerance {
        if density[pt] <= x {
            return (-1, Vec::with_capacity(0));
        }
    }
    let mut pn = pt;
    let mut dr = [0f64; 3];
    let mut path = vec![];
    let maxima = 'pathloop: loop {
        path.push(pn);
        pt = neargrid_step(pn, &mut dr, density);
        if path.contains(&pt) {
            pt = ongrid_step(pn, density);
            dr = [0f64; 3];
        }
        if pt == pn {
            break 'pathloop pn;
        } else {
            pn = pt;
        };
    };
    (maxima, Vec::with_capacity(0))
}

pub fn weight_step(p: isize,
                   density: &Density,
                   weight_map: WeightMap)
                   -> (isize, Vec<(usize, f64)>) {
    let control = density[p];
    let mut t_sum = 0.;
    let mut max_rho = 0.;
    let mut pn = p as isize;
    let mut weights = HashMap::<usize, f64>::new();
    // colllect the shift and distances and iterate over them
    for (shift, alpha) in
        density.voronoi.vectors.iter().zip(&density.voronoi.alphas)
    {
        let pt = density.voronoi_shift(p, shift);
        let charge_diff = density[pt] - control;
        if charge_diff > 0. {
            // calculate the gradient and if it is greater that the current max
            let rho = charge_diff * alpha;
            if rho > max_rho {
                pn = pt;
                max_rho = rho;
            }
            match weight_map.weight_get(pt) {
                Voxel::Weight(point_weights) => {
                    for (volume_number, w) in point_weights.iter() {
                        let weight =
                            weights.entry(*volume_number).or_insert(0.);
                        *weight += w * rho;
                    }
                }
                Voxel::Maxima(volume_number) => {
                    let weight = weights.entry(volume_number).or_insert(0.);
                    *weight += rho;
                }
                Voxel::Vacuum => (),
            }
            t_sum += rho;
        }
    }
    let weights = if weights.len() > 1 {
        weights.into_iter()
               .filter_map(|(volume_number, weight)| {
                   let weight = weight / t_sum;
                   if (weight * control).abs() > density.weight_tolerance {
                       Some((volume_number, weight))
                   } else {
                       None
                   }
               })
               .collect::<Vec<(usize, f64)>>()
    } else {
        Vec::with_capacity(0)
    };
    (pn, weights)
}

/// Finds the maxima associated with the current point, p.
pub fn weight(p: usize,
              density: &Density,
              weight_map: WeightMap)
              -> (isize, Vec<(usize, f64)>) {
    let pt = p as isize;
    let (pt, weights) = weight_step(pt, density, Arc::clone(&weight_map));
    if pt == p as isize {
        (pt, weights)
    } else {
        (weight_map.maxima_get(pt), weights)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::Lattice;

    #[test]
    fn methods_ongrid_step() {
        let data = (0..64).map(|x| x as f64).collect::<Vec<f64>>();
        let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
        let density = Density::new(&data,
                                   [4, 4, 4],
                                   lattice.to_cartesian,
                                   1E-8,
                                   Some(1E-3),
                                   [0., 0., 0.0]);
        assert_eq!(ongrid_step(61, &density), 62)
    }

    #[test]
    fn methods_neargrid_step() {
        let data = (0..64).map(|x| x as f64).collect::<Vec<f64>>();
        let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
        let density = Density::new(&data,
                                   [4, 4, 4],
                                   lattice.to_cartesian,
                                   1E-8,
                                   Some(1E-3),
                                   [0., 0., 0.0]);
        assert_eq!(neargrid_step(61, &mut [0., 0., -1.], &density), 61)
    }
}
