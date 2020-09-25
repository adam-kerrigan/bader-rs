use crate::density::Density;
use crate::utils;
use std::sync::atomic::{AtomicIsize, Ordering};
use std::sync::Arc;

/// A type for the maxima search function
pub type StepMethod = fn(usize, &Density, Arc<Vec<AtomicIsize>>) -> isize;

/// Steps in the density grid, from point p, following the gradient
pub fn ongrid_step(p: isize, density: &Density) -> isize {
    let mut pn = p.clone();
    let mut pt = p.clone();
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
    return pn;
}

/// Finds the maxima associated with the current point, p.
pub fn ongrid(p: usize,
              density: &Density,
              voxel_map: Arc<Vec<AtomicIsize>>)
              -> isize {
    let mut pt = p as isize;
    match density.vacuum_tolerance {
        Some(x) => {
            if density[pt] <= x {
                return -1;
            }
        }
        None => (),
    }
    let mut pn = pt;
    let maxima = 'pathloop: loop {
        pt = ongrid_step(pn, density);
        match voxel_map[pt as usize].load(Ordering::Relaxed) {
            -1 => (),
            x => {
                voxel_map[p].store(x, Ordering::Relaxed);
                return x;
            }
        }
        if pt == pn {
            break 'pathloop pn;
        } else {
            pn = pt;
        };
    };
    voxel_map[p].store(maxima, Ordering::Relaxed);
    return maxima;
}

/// Step in the density grid following the gradient.
///
/// Use the current point, p, and step along density and use the residual
/// gradient, dr, to account for grid bias.
pub fn neargrid_step(p: isize, dr: &mut [f64; 3], density: &Density) -> isize {
    let mut den_t = [0f64; 2];
    let mut grad = [0f64; 3];
    let mut max_grad = -1f64;
    let mut pt = p.clone();
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
    for i in 0..3 {
        if grad[i].abs() > max_grad {
            max_grad = grad[i].abs();
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
    return pt;
}

/// Find the maxima assiociated with the current point using the neargrid method
pub fn neargrid(p: usize,
                density: &Density,
                voxel_map: Arc<Vec<AtomicIsize>>)
                -> isize {
    let mut pt = p as isize;
    match density.vacuum_tolerance {
        Some(x) => {
            if density[pt] <= x {
                return -1;
            }
        }
        None => (),
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
    voxel_map[p].store(maxima, Ordering::Relaxed);
    return maxima;
}

pub fn weight_step(p: isize, density: &Density) -> isize {
    let mut pn = p.clone();
    let control = density[pn];
    let mut max_val = 0.;
    // colllect the shift and distances and iterate over them
    for (shift, alpha) in
        density.voronoi.vectors.iter().zip(&density.voronoi.alphas)
    {
        let pt = density.voronoi_shift(p, shift);
        // calculate the gradient and if it is greater that the current max
        let rho = (density[pt] - control) * alpha;
        if rho > max_val {
            max_val = rho;
            pn = pt;
        }
    }
    return pn;
}

/// Finds the maxima associated with the current point, p.
pub fn weight(p: usize,
              density: &Density,
              voxel_map: Arc<Vec<AtomicIsize>>)
              -> isize {
    let mut pt = p as isize;
    match density.vacuum_tolerance {
        Some(x) => {
            if density[pt] <= x {
                return -1;
            }
        }
        None => (),
    }
    let mut pn = pt;
    let maxima = 'pathloop: loop {
        pt = weight_step(pn, density);
        match voxel_map[pt as usize].load(Ordering::Relaxed) {
            -1 => (),
            x => {
                voxel_map[p].store(x, Ordering::Relaxed);
                return x;
            }
        }
        if pt == pn {
            break 'pathloop pn;
        } else {
            pn = pt;
        };
    };
    voxel_map[p].store(maxima, Ordering::Relaxed);
    return maxima;
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
                                   Some(1E-3),
                                   [0., 0., 0.0]);
        assert_eq!(ongrid_step(61, &density), 62)
    }

    #[test]
    fn methods_ongrid() {
        let data = (0..64).map(|x| x as f64).collect::<Vec<f64>>();
        let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
        let mut map = Vec::new();
        map.resize_with(64, || AtomicIsize::new(-1));
        let map = Arc::new(map);
        let density = Density::new(&data,
                                   [4, 4, 4],
                                   lattice.to_cartesian,
                                   Some(1E-3),
                                   [0., 0., 0.0]);
        assert_eq!(ongrid(61, &density, map), 63)
    }

    #[test]
    fn methods_neargrid_step() {
        let data = (0..64).map(|x| x as f64).collect::<Vec<f64>>();
        let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
        let density = Density::new(&data,
                                   [4, 4, 4],
                                   lattice.to_cartesian,
                                   Some(1E-3),
                                   [0., 0., 0.0]);
        assert_eq!(neargrid_step(61, &mut [0., 0., -1.], &density), 61)
    }

    #[test]
    fn methods_neargrid() {
        let data = (0..64).map(|x| x as f64).collect::<Vec<f64>>();
        let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
        let mut map = Vec::new();
        map.resize_with(64, || AtomicIsize::new(-1));
        let map = Arc::new(map);
        let density = Density::new(&data,
                                   [4, 4, 4],
                                   lattice.to_cartesian,
                                   Some(1E-3),
                                   [0., 0., 0.0]);
        assert_eq!(neargrid(61, &density, map), 63)
    }
}
