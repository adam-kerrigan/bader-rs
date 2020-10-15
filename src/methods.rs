use crate::density::Density;
use crate::utils;
use crate::voxel_map::{Voxel, VoxelMap};
use std::collections::HashMap;
use std::sync::Arc;

/// A type for the maxima search function
pub type WeightMap = Arc<VoxelMap>;
/// A type for monkey patching the method to available in the main.rs file.
pub type StepMethod =
    fn(usize, &Density, WeightMap) -> (isize, Vec<(usize, f64)>);

/// Indicates which method to use.
pub enum Method {
    /// The ongrid method.
    OnGrid,
    /// The neargrid method.
    NearGrid,
    /// The weighted method.
    Weight,
}

/// Steps in the density grid, from point p, following the gradient.
///
/// This should be called from [`ongrid()`] or [`neargrid()`].
///
/// * `p`: The point from which to step.
/// * `density`: The reference [`Density`].
///
/// ### Returns:
/// `isize`: The point with highest gradient.
///
/// # Examples
/// ```
/// use bader::density::Density;
/// use bader::atoms::Lattice;
/// use bader::methods::ongrid_step;
///
/// // Intialise the reference density.
/// let data = (0..64).map(|rho| rho as f64).collect::<Vec<f64>>();
/// let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
/// let density = Density::new(&data,
///                            [4, 4, 4],
///                            lattice.to_cartesian,
///                            1E-8,
///                            None,
///                            [0., 0., 0.]);
/// // The highest gradient between point, p = 33, and it's neighbours, with
/// // periodic boundary conditions, is with point p = 61.
/// assert_eq!(ongrid_step(33, &density), 61)
/// ```
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
///
/// * `p`: The point from which to step.
/// * `density`: The reference [`Density`].
/// * `weight_map`: An [`Arc`] wrapped [`VoxelMap`] for tracking the maxima.
///
/// ### Returns:
/// `(isize, Vec::with_capacity(0))`: The point current point's maxima and it's
/// weights. As this is not the weight method there are no weights.
///
/// # Examples
/// ```
/// use bader::density::Density;
/// use bader::atoms::Lattice;
/// use bader::voxel_map::VoxelMap;
/// use bader::methods::ongrid;
/// use std::sync::Arc;
///
/// // Intialise the reference density.
/// let data = (0..64).map(|rho| rho as f64).collect::<Vec<f64>>();
/// let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
/// let density = Density::new(&data,
///                            [4, 4, 4],
///                            lattice.to_cartesian,
///                            1E-8,
///                            None,
///                            [0., 0., 0.]);
/// let voxel_map = VoxelMap::new(64);
/// let voxel_map = Arc::new(voxel_map);
/// // The highest gradient between point, p = 33, and it's neighbours, with
/// // periodic boundary conditions, is with point p = 61 however as this hasn't
/// // been stored so this would block.
/// // assert_eq!(ongrid(33, &density, Arc::clone(&voxel_map)).0, -1);
/// // after stroing the maxima for the point, p = 61, ongrid will now return
/// // this value.
/// voxel_map.maxima_store(61, 63);
/// assert_eq!(ongrid(33, &density, Arc::clone(&voxel_map)).0, 63)
/// ```
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
///
/// This should be called from [`neargrid()`].
///
/// * `p`: The point from which to step.
/// * `dr`: A vector for storing the remainder of the unused step.
/// * `density`: The reference [`Density`].
///
/// ### Returns:
/// `isize`: The point with highest gradient.
///
/// # Examples
/// ```
/// use bader::density::Density;
/// use bader::atoms::Lattice;
/// use bader::methods::neargrid_step;
///
/// // Intialise the reference density.
/// let data = (0..64).map(|rho| rho as f64).collect::<Vec<f64>>();
/// let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
/// let density = Density::new(&data,
///                            [4, 4, 4],
///                            lattice.to_cartesian,
///                            1E-8,
///                            None,
///                            [0., 0., 0.]);
/// let mut dr = [0.; 3];
/// // The highest gradient between point, p = 33, and it's neighbours, with
/// // periodic boundary conditions, is with point p = 49.
/// assert_eq!(neargrid_step(33, &mut dr, &density), 49);
/// // with a dr showing the true gradient lies between p = 49, p = 61 and p = 62.
/// assert_eq!(dr, [0., -0.25, 0.0625]);
/// // performing this step another time, carrying forward dr, would result in a
/// // correction to the step from the rounding of the -0.5 from p = 49 to p = 61.
/// assert_eq!(neargrid_step(33, &mut dr, &density), 61);
/// // the new dr has the difference between the full step taken to correct
/// // and the value of dr when the correction was made.
/// assert_eq!(dr, [0., 0.5, 0.125]);
/// ```
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
///
/// * `p`: The point from which to step.
/// * `density`: The reference [`Density`].
/// * `_weight_map`: An unused [`Arc`] wrapped [`VoxelMap`] for tracking the maxima.
///
/// ### Returns:
/// `(isize, Vec::with_capacity(0))`: The point current point's maxima and it's
/// weights. As this is not the weight method there are no weights.
///
/// # Examples
/// ```
/// use bader::density::Density;
/// use bader::atoms::Lattice;
/// use bader::voxel_map::VoxelMap;
/// use bader::methods::neargrid;
/// use std::sync::Arc;
///
/// // Intialise the reference density.
/// let data = (0..64).map(|rho| rho as f64).collect::<Vec<f64>>();
/// let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
/// let density = Density::new(&data,
///                            [4, 4, 4],
///                            lattice.to_cartesian,
///                            1E-8,
///                            None,
///                            [0., 0., 0.]);
/// let voxel_map = VoxelMap::new(64);
/// let voxel_map = Arc::new(voxel_map);
/// // Unlike the other methods this follows the density gradient all the way to
/// // a maxima as the correction to the steps means that landing at a point
/// // doesn't mean the the original point's maxima will be the same.
/// assert_eq!(neargrid(22, &density, Arc::clone(&voxel_map)).0, 63)
/// ```
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

/// Steps in the density grid, from point p, following the gradient.
///
/// This should be called from [`weight()`].
///
/// Note: This function will deadlock if the points above it have no associated
/// maxima in [`VoxelMap.voxel_map`].
///
/// * `p`: The point from which to step.
/// * `density`: The reference [`Density`].
/// * `weight_map`: An [`Arc`] wrapped [`VoxelMap`] for tracking the maxima.
///
/// ### Returns:
/// `(isize, Vec::<f64>::with_capacity(14))`: The highest density gradient from
/// point, p, and it's associated weights.
///
/// # Examples
/// ```
/// use bader::density::Density;
/// use bader::atoms::Lattice;
/// use bader::voxel_map::VoxelMap;
/// use bader::methods::weight_step;
/// use std::sync::Arc;
///
/// // Intialise the reference density.
/// let data = (0..64).map(|rho| rho as f64).collect::<Vec<f64>>();
/// let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
/// let density = Density::new(&data,
///                            [4, 4, 4],
///                            lattice.to_cartesian,
///                            1E-8,
///                            None,
///                            [0., 0., 0.]);
/// let voxel_map = VoxelMap::new(64);
/// let voxel_map = Arc::new(voxel_map);
/// // The highest gradient between point, p = 33, and it's neighbours, with
/// // periodic boundary conditions, is with point p = 61.
///
/// // to avoid deadlock let's store maxima for all the values above us and
/// // store as either 61 or 62 to make the current point a boundary.
/// for i in [34, 37, 45, 49].iter() {
///     voxel_map.maxima_store(*i, 62 - i % 2);
/// }
/// let (pn, weights) = weight_step(33, &density, Arc::clone(&voxel_map));
/// let mut maxima = weights.into_iter().map(|(m, w)| m).collect::<Vec<usize>>();
/// maxima.sort_unstable();
/// assert_eq!(pn, 49);
/// assert_eq!(maxima, vec![61, 62]);
/// ```
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
///
/// Note: This function will deadlock if the points above it have no associated
/// maxima in [`VoxelMap.voxel_map`].
///
/// * `p`: The point from which to step.
/// * `density`: The reference [`Density`].
/// * `weight_map`: An [`Arc`] wrapped [`VoxelMap`] for tracking the maxima.
///
/// ### Returns:
/// `(isize, Vec::<f64>::with_capacity(14))`: The maxima associated with the
/// point, p, and it's associated weights.
///
/// # Examples
/// ```
/// use bader::density::Density;
/// use bader::atoms::Lattice;
/// use bader::voxel_map::VoxelMap;
/// use bader::methods::weight;
/// use std::sync::Arc;
///
/// // Intialise the reference density.
/// let data = (0..64).map(|rho| rho as f64).collect::<Vec<f64>>();
/// let lattice = Lattice::new([[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]]);
/// let density = Density::new(&data,
///                            [4, 4, 4],
///                            lattice.to_cartesian,
///                            1E-8,
///                            None,
///                            [0., 0., 0.]);
/// let voxel_map = VoxelMap::new(64);
/// let voxel_map = Arc::new(voxel_map);
/// // The highest gradient between point, p = 33, and it's neighbours, with
/// // periodic boundary conditions, is with point p = 61.
///
/// // to avoid deadlock let's store maxima for all the values above us and
/// // store as either 61 or 62 to make the current point a boundary.
/// for i in [34, 37, 45, 49].iter() {
///     voxel_map.maxima_store(*i, 62 - i % 2);
/// }
/// let (p_maxima, weights) = weight(33, &density, Arc::clone(&voxel_map));
/// let mut maxima = weights.into_iter().map(|(m, w)| m).collect::<Vec<usize>>();
/// maxima.sort_unstable();
/// assert_eq!(p_maxima, 61);
/// assert_eq!(maxima, vec![61, 62]);
/// ```
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
