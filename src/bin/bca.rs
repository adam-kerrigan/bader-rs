use anyhow::{bail, Context, Result};
use bader::analysis::{assign_maxima, sum_bader_densities};
use bader::arguments::{Args, ClapApp};
use bader::io::{self, FileFormat, FileType, WriteType};
use bader::methods::{maxima_finder, weight};
use bader::progress::Bar;
use bader::utils::vacuum_index;
use bader::voxel_map::{AtomVoxelMap, BlockingVoxelMap, VoxelMap};

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
    // find the maxima in the system and store them whilst removing them from
    // the index list
    let bader_maxima =
        maxima_finder(&mut index, reference, &voxel_map, args.threads, pbar)?;
    index.sort_unstable_by(|a, b| {
             reference[*b].partial_cmp(&reference[*a]).unwrap()
         });
    // remove from the indices any voxel that is below the vacuum limit
    index.truncate(vacuum_index(reference, &index, args.vacuum_tolerance)
         .context("Failed to apply vacuum tolerance")?);
    // Start a thread-safe progress bar and assign the maxima to atoms
    let pbar = Bar::visible(bader_maxima.len() as u64,
                            100,
                            String::from("Assigning to Atoms: "));
    let atom_map = assign_maxima(&bader_maxima,
                                 &atoms,
                                 &voxel_map.grid,
                                 &args.maximum_distance,
                                 args.threads,
                                 pbar)?;
    let pbar = Bar::visible(index.len() as u64,
                            100,
                            String::from("Bader Partitioning: "));
    // input the maxima as atoms into the voxel map
    bader_maxima.iter().enumerate().for_each(|(i, maxima)| {
                                       voxel_map.maxima_store(*maxima,
                                                              atom_map[i]
                                                              as isize);
                                   });
    // calculate the weights
    weight(reference,
           &voxel_map,
           &index,
           pbar,
           args.threads,
           args.weight_tolerance);
    // convert into a NonBlockingVoxelMap as the map is filled
    let (maxima_n, voxel_map): (usize, Box<dyn VoxelMap>) = {
        (atoms.positions.len(),
         Box::new(AtomVoxelMap::from_blocking_voxel_map(voxel_map)))
    };
    let pbar = Bar::visible(index.len() as u64,
                            100,
                            String::from("Summing Densities: "));
    // sum the densities and then write the charge partition files
    let (atoms_density, atoms_volume, min_surf_dist) =
        sum_bader_densities(&densities,
                            &voxel_map,
                            &atoms,
                            args.threads,
                            maxima_n,
                            pbar)?;
    // prepare the positions for writing out
    let positions = atoms.positions
                         .iter()
                         .map(|coords| file_type.coordinate_format(*coords))
                         .collect();
    // generate the output file
    let mut atoms_charge_file = io::output::partitions_file(positions,
                                                            &atoms_density,
                                                            &atoms_volume,
                                                            &total_density,
                                                            atoms.lattice
                                                                 .volume,
                                                            &min_surf_dist,
                                                            None)?;
    atoms_charge_file.push_str(&format!("\n  Bader Maxima: {}\n  Boundary Voxels: {}\n  Total Voxels: {}",
                                            bader_maxima.len(),
                                            voxel_map.boundary_iter().len(),
                                            reference.len()));
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
    // iterator as to save memory?
    let write_map: Box<dyn Iterator<Item = (isize, Vec<Option<f64>>)>> = {
        if let WriteType::Atom(a) = args.output {
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
        } else {
            bail!("BADBADNOTGOOD");
        }
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
