use crate::atoms::Atoms;
use crate::grid::Grid;
use crate::progress::Bar;
use crate::utils;
use crate::voxel_map::{Voxel, VoxelMap};
use anyhow::{bail, Context, Result};
use crossbeam_utils::thread;

/// Calculates the distance between a maxima and its nearest atom.
/// Chunk represents a collection of bader maxima positions withing the density
/// array.
fn maxima_to_atom(chunk: &[isize],
                  atoms: &Atoms,
                  grid: &Grid,
                  maximum_distance: &f64,
                  progress_bar: &Bar)
                  -> Result<Vec<usize>> {
    let chunk_size = chunk.len();
    // create vectors for storing the assigned atom and distance for each maxima
    let mut ass_atom = Vec::with_capacity(chunk_size);
    for m in chunk.iter() {
        // convert the point first to cartesian, then to the reduced basis
        let m_cartesian = grid.to_cartesian(*m);
        let m_reduced_cartesian = atoms.reduced_lattice.to_reduced(m_cartesian);
        let mut atom_num = 0;
        let mut min_distance = f64::INFINITY;
        // go through each atom in the reduced basis and shift in each
        // reduced direction, save the atom with the shortest distance
        for (i, atom) in atoms.reduced_positions.iter().enumerate() {
            for atom_shift in
                atoms.reduced_lattice.cartesian_shift_matrix.iter()
            {
                let distance = {
                    (m_reduced_cartesian[0]
                                        - (atom[0] + atom_shift[0]))
                                                                    .powi(2)
                                       + (m_reduced_cartesian[1]
                                          - (atom[1] + atom_shift[1]))
                                                                      .powi(2)
                                       + (m_reduced_cartesian[2]
                                          - (atom[2] + atom_shift[2]))
                                                                      .powi(2)
                };
                if distance < min_distance {
                    min_distance = distance;
                    atom_num = i;
                }
            }
        }
        if min_distance.powf(0.5) > *maximum_distance {
            bail!(
                "Bader maximum ({}, {}, {}) is too far from nearest atom ({}): {} Ang",
                m_cartesian[0],
                m_cartesian[1],
                m_cartesian[2],
                atom_num + 1,
                min_distance.powf(0.5),
            )
        }
        // remember to square root the distance
        ass_atom.push(atom_num);
        progress_bar.tick()
    }
    Ok(ass_atom)
}

/// Assign the Bader maxima to the nearest atom.
/// Threading will split the slice of maxima into chunks and operate on each
/// chunk in parallel.
pub fn assign_maxima(maxima: &[isize],
                     atoms: &Atoms,
                     grid: &Grid,
                     maximum_distance: &f64,
                     threads: usize,
                     progress_bar: Bar)
                     -> Result<Vec<usize>> {
    let mut assigned_atom = vec![0; maxima.len()];
    let pbar = &progress_bar;
    // this is basically a thread handling function for running the
    // maxima_to_atom function
    match threads.cmp(&1) {
        std::cmp::Ordering::Greater => {
            let chunk_size =
                (maxima.len() / threads) + (maxima.len() % threads).min(1);
            thread::scope(|s| {
                let spawned_threads =
                    maxima.chunks(chunk_size)
                          .enumerate()
                          .map(|(index, chunk)| {
                              s.spawn(move |_| {
                                  match maxima_to_atom(chunk, atoms, grid, maximum_distance, pbar) {
                                      Ok(result) => (result, index),
                                      _ => panic!("Failed to match maxima to atom"),
                                  }
                               })
                          })
                          .collect::<Vec<_>>();
                for thread in spawned_threads {
                    if let Ok((ass_atom, chunk_index)) =
                        thread.join()
                    {
                        // is this required? is the collection of handles not
                        // already sorted like this, is it possible to join as
                        // they finish?
                        let i = chunk_index * chunk_size;
                        assigned_atom.splice(i..(i + ass_atom.len()), ass_atom);
                    } else {
                        panic!("Failed to join thread in assign maxima.")
                    };
                }
            }).unwrap();
        }
        _ => {
            let ass_atom =
                maxima_to_atom(maxima, atoms, grid, maximum_distance, pbar).context("Failed to assign maxima to atom.")?;
            assigned_atom = ass_atom;
        }
    }
    Ok(assigned_atom)
}

/// Sums the densities of each Bader volume.
pub fn calculate_bader_density(density: &[f64],
                               voxel_map: &impl VoxelMap,
                               atoms: &Atoms,
                               threads: usize,
                               progress_bar: Bar)
                               -> Box<[f64]> {
    let pbar = &progress_bar;
    let mut bader_density = vec![0.0; atoms.positions.len() + 1];
    let vm = &voxel_map;
    // Calculate the size of the vector to be passed to each thread.
    let chunk_size =
        (density.len() / threads) + (density.len() % threads).min(1);
    thread::scope(|s| {
        let spawned_threads =
            voxel_map.maxima_chunks(chunk_size)
                     .enumerate()
                     .map(|(index, chunk)| {
                         s.spawn(move |_| {
                              let mut bd = vec![0.0; atoms.positions.len() + 1];
                              chunk.iter()
                                   .enumerate()
                                   .for_each(|(voxel_index, maxima)| {
                                       let p =
                                           index * chunk.len() + voxel_index;
                                       match vm.maxima_to_voxel(*maxima) {
                                           Voxel::Maxima(m) => {
                                               bd[m] += density[p];
                                           }
                                           Voxel::Boundary(weights) => {
                                               for weight in weights.iter() {
                                                   let m = *weight as usize;
                                                   let w = weight - (m as f64);
                                                   bd[m] += w * density[p];
                                               }
                                           }
                                           Voxel::Vacuum => {
                                               bd[atoms.positions.len()] +=
                                                   density[p]
                                           }
                                       };
                                       pbar.tick();
                                   });
                              bd
                          })
                     })
                     .collect::<Vec<_>>();
        // Join each thread and collect the results.
        // If one thread terminates before the other this is not operated on first.
        // Either use the sorted index to remove vacuum from the summation or
        // find a way to operate on finshed threads first (ideally both).
        for thread in spawned_threads {
            if let Ok(tmp_bd) = thread.join() {
                bader_density.iter_mut()
                             .zip(tmp_bd.into_iter())
                             .for_each(|(a, b)| {
                                 *a += b;
                             });
            } else {
                panic!("Unable to join thread in sum_bader_densities.")
            };
        }
    }).unwrap();
    // The distance isn't square rooted in the calcation of distance to save time.
    // As we need to filter out the infinite distances (atoms with no assigned maxima)
    // we can square root here also.
    bader_density.iter_mut().for_each(|a| {
                                *a *= voxel_map.grid_get().voxel_lattice.volume;
                            });
    bader_density.into()
}

/// Calculates the volume and radius of each Bader atom.
pub fn calculate_bader_volume_radius(density: &[f64],
                                     voxel_map: &impl VoxelMap,
                                     atoms: &Atoms,
                                     threads: usize,
                                     progress_bar: Bar)
                                     -> (Box<[f64]>, Box<[f64]>) {
    let pbar = &progress_bar;
    let mut bader_radius = vec![f64::INFINITY; atoms.positions.len()];
    let mut bader_volume = vec![0.0; atoms.positions.len() + 1];
    let vm = &voxel_map;
    // Calculate the size of the vector to be passed to each thread.
    let chunk_size =
        (density.len() / threads) + (density.len() % threads).min(1);
    thread::scope(|s| {
        let spawned_threads =
            voxel_map.maxima_chunks(chunk_size)
                     .enumerate()
                     .map(|(index, chunk)| s.spawn(move |_| {
                         let mut br = vec![f64::INFINITY; atoms.positions.len()];
                         let mut bv = vec![0.0; atoms.positions.len() + 1];
                         chunk.iter()
                              .enumerate()
                              .for_each(|(voxel_index, maxima)| {
                                  let p = index * chunk.len() + voxel_index;
                                  match vm.maxima_to_voxel(*maxima) {
                                      Voxel::Boundary(weights) => {
                                          for weight in weights.iter() {
                                              let m = *weight as usize;
                                              let w = weight - (m as f64);
                                              bv[m] += w;
                                              let atom_number = vm.maxima_to_atom(m);
                                              let mr = &mut br[atom_number];
                                              let p_c = vm.grid_get().to_cartesian(p as isize);
                                              let mut p_lll_f = utils::dot(p_c, atoms.reduced_lattice.to_fractional);
                                              for f in &mut p_lll_f {
                                                  *f = f.rem_euclid(1.);
                                              }
                                              let p_lll_c = utils::dot(p_lll_f, atoms.reduced_lattice.to_cartesian);
                                              let atom = atoms.reduced_positions[atom_number];
                                              for atom_shift in
                                                  atoms.reduced_lattice.cartesian_shift_matrix.iter()
                                              {
                                                  let distance = {
                                                      (p_lll_c[0]
                                                           - (atom[0] + atom_shift[0]))
                                                                                       .powi(2)
                                                          + (p_lll_c[1]
                                                             - (atom[1] + atom_shift[1]))
                                                                                         .powi(2)
                                                          + (p_lll_c[2]
                                                             - (atom[2] + atom_shift[2]))
                                                                                         .powi(2)
                                                  };
                                                  if distance < *mr {
                                                      *mr = distance;
                                                  }
                                              }
                                          }
                                      },
                                      Voxel::Maxima(m) => {bv[m] += 1.0;},
                                      Voxel::Vacuum => {bv[atoms.positions.len()] += 1.0;},
                                  };
                                  pbar.tick();
                              });
                        (bv, br)
                     }))
                     .collect::<Vec<_>>();
        // Join each thread and collect the results.
        // If one thread terminates before the other this is not operated on first.
        // Either use the sorted index to remove vacuum from the summation or
        // find a way to operate on finshed threads first (ideally both).
        for thread in spawned_threads {
            if let Ok((tmp_bv, tmp_br)) = thread.join() {
                bader_volume.iter_mut()
                            .zip(tmp_bv.into_iter())
                            .for_each(|(a, b)| {
                                *a += b;
                            });
                bader_radius.iter_mut()
                                .zip(tmp_br.into_iter())
                                .for_each(|(a, b)| {
                                    *a = a.min(b);
                                });
            } else {
                panic!("Unable to join thread in sum_bader_densities.")
            };
        }
    }).unwrap();
    // The distance isn't square rooted in the calcation of distance to save time.
    // As we need to filter out the infinite distances (atoms with no assigned maxima)
    // we can square root here also.
    bader_volume.iter_mut().for_each(|a| {
                               *a *= voxel_map.grid_get().voxel_lattice.volume;
                           });
    bader_radius.iter_mut()
                .for_each(|d| match (*d).partial_cmp(&f64::INFINITY) {
                    Some(std::cmp::Ordering::Less) => *d = d.powf(0.5),
                    _ => *d = 0.0,
                });
    (bader_volume.into(), bader_radius.into())
}
