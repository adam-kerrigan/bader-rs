use crate::arguments::Weight;
use crate::atoms::Atoms;
use crate::density::Density;
use crate::progress::Bar;
use crate::utils;
use indicatif::ProgressBar;
use std::ops::Index;

enum Boundary {
    Weight((Vec<(isize, f64)>, usize, bool)),
    None(bool),
}

/// VoxelMap - enum for mapping the voxels to bader volumes
///
/// > Just(Vec<usize>) - Just the voxel map (ongrid)
/// > Known { map: Vec<usize>, known: Vec<bool> } - Holds the map and a known
/// >                                               array
pub struct VoxelMap {
    pub map: Vec<isize>,
    pub bader_maxima: Vec<isize>,
    index: Vec<usize>,
}

impl Index<isize> for VoxelMap {
    type Output = isize;

    /// Index the reference charge inside the Density structure
    fn index(&self, i: isize) -> &Self::Output {
        &self.map[i as usize]
    }
}

impl Index<usize> for VoxelMap {
    type Output = isize;

    /// Index the reference charge inside the Density structure
    fn index(&self, i: usize) -> &Self::Output {
        &self.map[i]
    }
}

impl VoxelMap {
    /// Initialise the structure. Create a HashMap to index the volumes.
    pub fn new(map: Vec<isize>, bader_maxima: Vec<isize>) -> Self {
        let mut index = vec![0usize; map.len()];
        if bader_maxima[0] < 0 {
            for (i, maxima) in bader_maxima[1..].iter().enumerate() {
                index[*maxima as usize] = i;
            }
        } else {
            for (i, maxima) in bader_maxima.iter().enumerate() {
                index[*maxima as usize] = i;
            }
        }
        Self { map,
               bader_maxima,
               index }
    }

    pub fn index_get(&self, p: isize) -> usize {
        self.index[self[p] as usize]
    }

    /// Is the point known, use shifts to move around
    fn is_boundary_atoms(map: &Self,
                         p: isize,
                         atom_number: usize,
                         density: &Density,
                         assigned_atom: &Vec<usize>,
                         weight_map: &utils::BTMap)
                         -> Boundary {
        let mut t = Vec::<(isize, f64)>::with_capacity(14);
        let mut t_total = 0.;
        let rho = density[p];
        let mut is_boundary = false;
        let mut is_atom = false;
        let mut count = 0;
        for (shift, alpha) in
            density.voronoi.vectors.iter().zip(&density.voronoi.alphas)
        {
            let pn = density.voronoi_shift(p, shift);
            if density[pn] > rho {
                if !is_atom {
                    match weight_map.get_sly(pn) {
                        Some(volume_numbers) => {
                            for (vn, _) in volume_numbers.iter() {
                                if atom_number
                                   != assigned_atom[map.index[*vn as usize]]
                                {
                                    is_atom = true;
                                }
                            }
                        }
                        None => {
                            if atom_number != assigned_atom[map.index_get(pn)] {
                                is_atom = true;
                            }
                        }
                    }
                }
                t.push((pn, alpha * (density[pn] - rho)));
                t_total += t.last().unwrap().1;
            } else {
                count += 1;
                if map[pn] >= 0 {
                    if atom_number != assigned_atom[map.index_get(pn)] {
                        is_boundary = true;
                    }
                } else {
                    is_boundary = true;
                }
            }
        }
        if is_atom {
            Boundary::Weight((t.into_iter()
                               .map(|(i, x)| (i, x / t_total))
                               .collect(),
                              count,
                              true))
        } else {
            Boundary::None(is_boundary)
        }
    }

    fn is_boundary_volumes(map: &Self,
                           p: isize,
                           atom_number: usize,
                           density: &Density,
                           assigned_atom: &Vec<usize>,
                           weight_map: &utils::BTMap)
                           -> Boundary {
        let volume_number = map[p];
        let mut t = Vec::<(isize, f64)>::with_capacity(14);
        let mut t_total = 0.;
        let rho = density[p];
        let mut is_boundary = false;
        let mut is_volume = false;
        let mut count = 0;
        for (shift, alpha) in
            density.voronoi.vectors.iter().zip(&density.voronoi.alphas)
        {
            let pn = density.voronoi_shift(p, shift);
            if density[pn] > rho {
                if !is_boundary || !is_volume {
                    match weight_map.get_sly(pn) {
                        Some(volume_numbers) => {
                            for (vn, _) in volume_numbers.into_iter() {
                                if atom_number
                                   != assigned_atom[map.index[vn as usize]]
                                {
                                    is_boundary = true;
                                }
                            }
                            is_volume = true;
                        }
                        None => {
                            if map[pn] != volume_number {
                                if atom_number
                                   != assigned_atom[map.index_get(pn)]
                                {
                                    is_boundary = true;
                                }
                                is_volume = true
                            }
                        }
                    }
                }
                t.push((pn, alpha * (density[pn] - rho)));
                t_total += t.last().unwrap().1;
            } else {
                count += 1;
                if map[pn] >= 0 {
                    if atom_number != assigned_atom[map.index_get(pn)] {
                        is_boundary = true;
                    }
                } else {
                    is_boundary = true;
                }
            }
        }
        if is_volume {
            Boundary::Weight((t.into_iter()
                               .map(|(i, x)| (i, x / t_total))
                               .collect(),
                              count,
                              is_boundary))
        } else {
            Boundary::None(is_boundary)
        }
    }

    fn is_boundary(map: &VoxelMap,
                   p: isize,
                   atom_number: usize,
                   density: &Density,
                   assigned_atom: &Vec<usize>,
                   _weight_map: &utils::BTMap)
                   -> Boundary {
        for shift in density.voronoi.vectors.iter() {
            let pn = density.voronoi_shift(p, shift);
            if map[pn] >= 0 {
                if atom_number != assigned_atom[map.index_get(pn)] {
                    return Boundary::None(true);
                }
            } else {
                return Boundary::None(true);
            }
        }
        Boundary::None(false)
    }

    /// Sum the densities for each bader volume
    pub fn charge_sum(&self,
                      densities: &Vec<Vec<f64>>,
                      assigned_atom: &Vec<usize>,
                      atoms: &Atoms,
                      density: &Density,
                      index: Vec<usize>,
                      weight: Weight)
                      -> (Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
        type WeightMethod = fn(&VoxelMap,
                               isize,
                               usize,
                               &Density,
                               &Vec<usize>,
                               &utils::BTMap)
                               -> Boundary;
        let is_boundary: WeightMethod = match weight {
            Weight::Atoms => Self::is_boundary_atoms,
            Weight::Volumes => Self::is_boundary_volumes,
            Weight::None => Self::is_boundary,
        };
        let mut weight_map = utils::BTMap::new(density.size.total);
        let mut bader_volume = vec![0.; self.bader_maxima.len()];
        let mut bader_charge =
            { vec![vec![0f64; self.bader_maxima.len()]; densities.len()] };
        let mut surface_distance = vec![0f64; atoms.positions.len()];
        let mut minimum_distance = {
            vec![
                atoms.lattice.distance_matrix[0].powi(2);
                atoms.positions.len()
            ]
        };

        let pbar = ProgressBar::new(density.size.total as u64);
        let pbar = Bar::new(pbar, 100, String::from("Summing Charge: "));
        for p in index.into_iter() {
            let maxima = self[p];
            if maxima < 0 {
                for j in 0..densities.len() {
                    bader_charge[j][self.bader_maxima.len() - 1] +=
                        densities[j][p];
                }
                bader_volume[self.bader_maxima.len() - 1] += 1.;
                continue;
            }
            let i = self.index[maxima as usize];
            let atom_num = assigned_atom[i];
            let boundary = match is_boundary(&self,
                                             p as isize,
                                             atom_num,
                                             density,
                                             assigned_atom,
                                             &weight_map)
            {
                Boundary::Weight((probabilities, count, is_boundary)) => {
                    let weights = density.probability(probabilities,
                                                      self,
                                                      &mut weight_map);
                    for (maxima, w) in weights.iter() {
                        let i = self.index[*maxima as usize];
                        bader_volume[i] += w;
                        for j in 0..densities.len() {
                            bader_charge[j][i] += w * densities[j][p];
                        }
                    }
                    weight_map.insert(p as isize, weights, count);
                    is_boundary
                }
                Boundary::None(boundary) => {
                    bader_volume[i] += 1.;
                    for j in 0..densities.len() {
                        bader_charge[j][i] += densities[j][p];
                    }
                    boundary
                }
            };
            if boundary {
                let px = density.voxel_origin[0]
                         + (p as isize / (density.size.y * density.size.z))
                           as f64;
                let py =
                    density.voxel_origin[1]
                    + (p as isize / density.size.z).rem_euclid(density.size.y)
                      as f64;
                let pz = density.voxel_origin[2]
                         + (p as isize).rem_euclid(density.size.z) as f64;
                let p_cartesian = utils::dot([px, py, pz],
                                             density.voxel_lattice
                                                    .to_cartesian);
                let mut p_lll_fractional = utils::dot(p_cartesian,
                                                      atoms.reduced_lattice
                                                           .to_fractional);
                for f in &mut p_lll_fractional {
                    *f = f.rem_euclid(1.);
                }
                let p_lll_cartesian = utils::dot(p_lll_fractional,
                                                 atoms.reduced_lattice
                                                      .to_cartesian);
                let atom = atoms.reduced_positions[atom_num];
                for atom_shift in
                    atoms.reduced_lattice.cartesian_shift_matrix.iter()
                {
                    let distance = {
                        (p_lll_cartesian[0] - (atom[0] + atom_shift[0])).powi(2)
                        + (p_lll_cartesian[1] - (atom[1] + atom_shift[1])).powi(2)
                        + (p_lll_cartesian[2] - (atom[2] + atom_shift[2])).powi(2)
                    };
                    if distance < minimum_distance[atom_num] {
                        surface_distance[atom_num] = distance.powf(0.5);
                        minimum_distance[atom_num] = distance;
                    }
                }
            }
            pbar.tick();
        }
        drop(pbar);
        (bader_charge, bader_volume, surface_distance)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn voxel_map_new_vacuum() {
        let map: Vec<isize> = vec![1, 1, 1, 2, 2, 1, -1, 2, 2, -1, 2, 1];
        let bader_maxima = vec![-1, 1, 2];
        let voxel_map = VoxelMap::new(map, bader_maxima);
        assert_eq!(voxel_map.index[1], 0);
        assert_eq!(voxel_map.index[2], 1);
    }

    #[test]
    fn voxel_map_new_no_vacuum() {
        let map: Vec<isize> = vec![1, 1, 1, 2, 2, 1, 2, 2, 2, 1];
        let bader_maxima = vec![1, 2];
        let voxel_map = VoxelMap::new(map, bader_maxima);
        assert_eq!(voxel_map.index[1], 0);
        assert_eq!(voxel_map.index[2], 1);
    }
}
