use anyhow::{bail, Context, Result};
use bader::analysis::{
    assign_maxima, sum_atoms_densities, sum_bader_densities,
};
use bader::arguments::{Args, ClapApp, Verbosity};
use bader::io::{self, FileFormat, FileType, WriteType};
use bader::methods::{maxima_finder, weight};
use bader::progress::Bar;
use bader::utils::vacuum_index;
use bader::voxel_map::{
    AtomVoxelMap, BaderVoxelMap, BlockingVoxelMap, VoxelMap,
};
use rustc_hash::FxHashSet;

fn main() -> Result<()> {
    // argument parsing
    let app = ClapApp::get();
    let args = Args::new(app.get_matches());
    // print splash
    println!("Multi-threaded Bader Charge Analysis ({})",
             env!("CARGO_PKG_VERSION"));
    // read the input files into a densities vector and a Grid struct
    let file_type: Box<dyn FileFormat> = match args.file_type {
        FileType::Vasp => Box::new(io::vasp::Vasp {}),
        FileType::Cube => Box::new(io::cube::Cube {}),
    };
    println!("Running on {} threads.", args.threads);
    let (densities, rho, atoms, grid, voxel_origin) = file_type.init(&args);
    let reference = if rho.is_empty() { &densities[0] } else { &rho };
    let voxel_map =
        BlockingVoxelMap::new(grid, atoms.lattice.to_cartesian, voxel_origin);
    let total_density = densities.iter()
                                 .map(|d| {
                                     d.iter().sum::<f64>()
                                     * voxel_map.grid.voxel_lattice.volume
                                 })
                                 .collect::<Vec<f64>>();
    // create the index list which will tell us in which order to evaluate the
    // voxels
    let mut index: Vec<usize> = (0..voxel_map.grid.size.total).collect();
    let pbar =
        Bar::visible(index.len() as u64, 100, String::from("Maxima Finding: "));
    let bader_maxima =
        maxima_finder(&mut index, reference, &voxel_map, args.threads, pbar)?;
    index.sort_unstable_by(|a, b| {
             reference[*b].partial_cmp(&reference[*a]).unwrap()
         });
    // remove from the indices any voxel that is below the vacuum limit
    index.truncate(vacuum_index(reference, &index, args.vacuum_tolerance)
         .context("Failed to apply vacuum tolerance")?);
    // Start a thread-safe progress bar and run the main calculation
    let pbar = Bar::visible(bader_maxima.len() as u64,
                            100,
                            String::from("Assigning to Atoms: "));
    let (atom_map, minimum_distance) = assign_maxima(&bader_maxima,
                                                     &atoms,
                                                     &voxel_map.grid,
                                                     args.threads,
                                                     pbar)?;
    let pbar = Bar::visible(index.len() as u64,
                            100,
                            String::from("Bader Partitioning: "));
    // input the maxima into the voxel map
    if let Verbosity::Atoms = args.verbosity {
        bader_maxima.iter().enumerate().for_each(|(i, maxima)| {
                                           voxel_map.maxima_store(*maxima,
                                                                  atom_map[i]
                                                                  as isize);
                                       });
    } else {
        bader_maxima.iter()
                    .enumerate()
                    .for_each(|(i, maxima)| {
                        voxel_map.maxima_store(*maxima, i as isize);
                    })
    }
    // calculate the weights
    weight(reference,
           &voxel_map,
           &index,
           pbar,
           args.threads,
           args.weight_tolerance);
    // convert into a NonBlockingVoxelMap as the map is filled
    let (maxima_n, voxel_map): (usize, Box<dyn VoxelMap>) = match args.verbosity
    {
        Verbosity::Atoms => {
            (atoms.positions.len(),
             Box::new(AtomVoxelMap::from_blocking_voxel_map(voxel_map)))
        }
        _ => (atom_map.len(),
              Box::new(BaderVoxelMap::from_blocking_voxel_map(voxel_map,
                                                              atom_map))),
    };
    let pbar = Bar::visible(index.len() as u64,
                            100,
                            String::from("Summing Densities: "));
    // sum the densities and then write the charge partition files
    let (bader_density, bader_volume, min_surf_dist) =
        sum_bader_densities(&densities,
                            &voxel_map,
                            &atoms,
                            args.threads,
                            maxima_n,
                            pbar)?;
    let (atoms_density, atoms_volume) = match args.verbosity {
        Verbosity::Atoms => (bader_density, bader_volume),
        _ => {
            let positions = bader_maxima.iter()
                                        .map(|p| {
                                            let p = voxel_map.grid_get()
                                                             .to_cartesian(*p);
                                            file_type.coordinate_format(p)
                                        })
                                        .collect();
            let bader_charge_file =
                io::output::partitions_file(positions,
                                            &bader_density,
                                            &bader_volume,
                                            &total_density,
                                            atoms.lattice.volume,
                                            &minimum_distance,
                                            voxel_map.atom_map())?;
            // check that the write was successfull
            io::output::write(bader_charge_file, String::from("BCF.dat"))?;
            sum_atoms_densities(&bader_density,
                                &bader_volume,
                                voxel_map.atom_map().unwrap(),
                                atoms.positions.len())?
        }
    };
    let positions = atoms.positions
                         .iter()
                         .map(|coords| file_type.coordinate_format(*coords))
                         .collect();
    let mut atoms_charge_file = io::output::partitions_file(positions,
                                                            &atoms_density,
                                                            &atoms_volume,
                                                            &total_density,
                                                            atoms.lattice
                                                                 .volume,
                                                            &min_surf_dist,
                                                            None)?;
    if let Verbosity::Full = args.verbosity {
        atoms_charge_file.push_str(&format!("\n  Bader Maxima: {}\n  Boundary Voxels: {}\n  Total Voxels: {}",
                                            bader_maxima.len(),
                                            voxel_map.boundary_iter().len(),
                                            reference.len())
                                   );
    }
    // check that the write was successfull
    io::output::write(atoms_charge_file, String::from("ACF.dat"))?;
    // Prepare to write any densities that have been requested.
    let filename = match densities.len().cmp(&2) {
        std::cmp::Ordering::Less => vec![String::from("charge")],
        std::cmp::Ordering::Equal => {
            vec![String::from("charge"), String::from("spin")]
        }
        std::cmp::Ordering::Greater => vec![String::from("charge"),
                                            String::from("spin_x"),
                                            String::from("spin_y"),
                                            String::from("spin_z"),],
    };
    // create a map that has an optional weight for each voxel and store that
    // with an id for each volume that is to be outputted. Save this as a lazy
    // iterator as to save memory? This is now a large part of the binary,
    // should it be moved?
    let write_map: Box<dyn Iterator<Item = (isize, Vec<Option<f64>>)>> =
        match (args.output, args.verbosity) {
            (WriteType::Volume(_), Verbosity::Atoms) => bail!(
"Unable to write Bader volumes at this level of verbosity, either increase
verbosity or export atoms."
                ),
            (WriteType::Volume(v), _) => {
                let volume_iter = if v.is_empty() {
                    (0..bader_maxima.len() as isize).collect()
                } else {
                    v
                };
                Box::new(volume_iter.into_iter()
                       .map(|volume_number| (volume_number, voxel_map.volume_map(volume_number))))
            }
            (WriteType::Atom(a), Verbosity::Atoms) => {
                let atom_iter = if a.is_empty() {
                    (0..atoms.positions.len() as isize).collect()
                } else {
                    a
                };
                Box::new(atom_iter.into_iter()
                     .map(|atom_number| {
                         let map = voxel_map.volume_map(atom_number);
                         (atom_number, map)
                     }))
            }
            (WriteType::Atom(a), _) => {
                let atom_iter: Vec<FxHashSet<isize>> = if a.is_empty() {
                    let mut a_i =
                        vec![FxHashSet::default(); atoms.positions.len()];
                    voxel_map.atom_map()
                             .unwrap()
                             .iter()
                             .enumerate()
                             .for_each(|(i, atom)| {
                                 a_i[*atom].insert(i as isize);
                             });
                    a_i
                } else {
                    a.iter()
                     .map(|atom_number| {
                         voxel_map.atom_map()
                                  .unwrap()
                                  .iter()
                                  .enumerate()
                                  .filter_map(|(i, atom)| {
                                      if (*atom as isize) == *atom_number {
                                          Some(i as isize)
                                      } else {
                                          None
                                      }
                                  })
                                  .collect::<FxHashSet<isize>>()
                     })
                     .collect()
                };
                Box::new(atom_iter.into_iter()
                     .filter_map(|volume_numbers| {
                         if !volume_numbers.is_empty() {
                             let map = voxel_map.multi_volume_map(&volume_numbers);
                             let vn = *volume_numbers.iter().next().unwrap() as usize;
                             let atom_number = voxel_map.maxima_to_atom(vn) as isize;
                             Some((atom_number, map))
                         } else {
                             None
                         }
                     }))
            }
            (WriteType::None, _) => Box::new(Vec::with_capacity(0).into_iter()),
        };
    write_map.for_each(|(id, weight_map)| {
        densities.iter()
                 .zip(&filename)
                 .for_each(|(rho, flnm)| {
                     let pbar = Bar::visible(index.len() as u64,
                                             100,
                                             format!("Writing file {}:", flnm));
                     if file_type.write(
                         &atoms,
                         weight_map.iter()
                                   .zip(rho)
                                   .map(|(weight, charge)| weight.map(|w| w * charge))
                                   .collect(),
                         format!("{}_{}", id + 1, flnm),
                         pbar).is_err()
                     {
                         panic!("Error in writing {}", flnm)
                     }
                 });
    });
    Ok(())
}
