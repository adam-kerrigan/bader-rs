use crate::density::Density;
use crate::utils;

/// A type for the maxima search function
pub type StepMethod = fn(usize, &Density) -> isize;

/// Steps in the density grid, from point p, following the gradient
pub fn ongrid_step(p: isize, density: &Density) -> isize {
    let mut pn = p.clone();
    let mut pt = p.clone();
    let control = density[pn];
    let mut max_val = density[pt];
    let distance = density.voxel_lattice.distance_matrix;
    // colllect the shift and distances and iterate over them
    for (i, p_shift) in density.full_shift(p).iter().enumerate() {
        pt = p + p_shift;
        // calculate the gradient and if it is greater that the current max
        let rho = control + (density[pt] - control) / distance[i];
        if rho > max_val {
            max_val = rho;
            pn = pt;
        }
    }
    return pn;
}

/// Finds the maxima associated with the current point, p.
pub fn ongrid(p: usize, density: &Density) -> isize {
    let mut pt = p as isize;
    match density.vacuum_tolerance {
        Some(x) => {
            if density[pt] <= x {
                return -1isize;
            }
        }
        None => (),
    }
    let mut pn = pt;
    let maxima = 'pathloop: loop {
        pt = ongrid_step(pn, density);
        if pt == pn {
            break 'pathloop pn;
        } else {
            pn = pt;
        };
    };
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
pub fn neargrid(p: usize, density: &Density) -> isize {
    let mut pt = p as isize;
    match density.vacuum_tolerance {
        Some(x) => {
            if density[pt] <= x {
                return -1isize;
            }
        }
        None => (),
    }
    let mut pn = pt;
    let mut dr = [0f64; 3];
    let mut path = vec![];
    let vol_num = 'pathloop: loop {
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
    return vol_num;
}
