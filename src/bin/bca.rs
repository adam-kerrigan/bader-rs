use anyhow::{Context, Result};
use bader::analysis::{assign_maxima, sum_bader_densities, Analysis};
use bader::arguments::{Args, ClapApp};
use bader::io::{self, FileFormat, FileType};
use bader::methods::weight;
use bader::progress::Bar;
use bader::utils::vacuum_index;
use bader::voxel_map::{BlockingVoxelMap, NonBlockingVoxelMap};
use rustc_hash::FxHashMap;
use std::time::{Duration, Instant};

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

    let (densities, rho, atoms, grid, voxel_origin) = file_type.init(&args);
    let reference = if rho.is_empty() { &densities[0] } else { &rho };
    let voxel_map =
        BlockingVoxelMap::new(grid, atoms.lattice.to_cartesian, voxel_origin);
    let mut index: Vec<usize> = (0..voxel_map.grid.size.total).collect();
    index.sort_unstable_by(|a, b| {
             reference[*b].partial_cmp(&reference[*a]).unwrap()
         });
    let vacuum_index = vacuum_index(reference, &index, args.vacuum_tolerance).context("Failed to apply vacuum tolerance")?;
    // Start a thread-safe progress bar and run the main calculation
    let pbar = Bar::visible(vacuum_index as u64,
                            100,
                            String::from("Bader Partitioning: "));
    weight(reference,
           &voxel_map,
           &index,
           pbar,
           args.threads,
           vacuum_index,
           args.weight_tolerance);
    let voxel_map = NonBlockingVoxelMap::from_blocking_voxel_map(voxel_map);
    let mut analysis =
        Analysis::new(&voxel_map, densities.len(), atoms.positions.len());
    // find the nearest atom to each Bader maxima
    let pbar = Bar::visible(analysis.bader_maxima.len() as u64,
                            100,
                            String::from("Assigning to Atoms: "));
    let start = Instant::now();
    analysis.assign_atoms(&atoms, &voxel_map.grid, pbar);
    let duration = start.elapsed();
    println!("analysis assign: {:?}", duration);
    let pbar = Bar::visible(voxel_map.grid.size.total as u64,
                            100,
                            String::from("summing charge: "));
    // Sum the charge in each volume
    let start = Instant::now();
    if let Err(e) = analysis.charge_sum(&atoms, &densities, &voxel_map, pbar) {
        println!("{}", e);
    }
    let duration = start.elapsed();
    println!("analysis charge sum: {:?}", duration);
    // Sum the charge in each atom
    analysis.atoms_charge_sum();
    // build the results
    println!("Writing output files:");
    let (atoms_charge_file, bader_charge_file) =
        io::output::charge_files(&analysis,
                                 &atoms,
                                 &voxel_map.grid,
                                 &file_type,
                                 args.maxima_tolerance);
    // check that the write was successfull
    if let Err(e) = io::output::write(atoms_charge_file, bader_charge_file) {
        panic!("Error occured: {}", e);
    }
    println!("ACF.dat and BCF.dat written successfully.");

    // START TESTING
    let pbar = Bar::visible(voxel_map.maxima_map.len() as u64,
                            100,
                            String::from("Assigning to Atoms: "));
    let start = Instant::now();
    let bader_maxima = voxel_map.maxima_map
                                .keys()
                                .map(|maxima| *maxima as isize)
                                .collect::<Vec<isize>>();
    let (assigned_atom, minimum_distance) =
        assign_maxima(&bader_maxima, &atoms, &voxel_map.grid, 1, pbar)?;
    let atom_map: FxHashMap<usize, usize> =
        assigned_atom.iter()
                     .enumerate()
                     .map(|(i, atom)| (bader_maxima[i] as usize, *atom))
                     .collect();
    let duration = start.elapsed();
    println!("assign atoms 1 thread: {:?}", duration);
    let pbar = Bar::visible(voxel_map.maxima_map.len() as u64,
                            100,
                            String::from("Assigning to Atoms: "));
    let start = Instant::now();
    let bader_maxima = voxel_map.maxima_map
                                .keys()
                                .map(|maxima| *maxima as isize)
                                .collect::<Vec<isize>>();
    let (assigned_atom, minimum_distance) =
        assign_maxima(&bader_maxima, &atoms, &voxel_map.grid, 2, pbar)?;
    let atom_map: FxHashMap<usize, usize> =
        assigned_atom.iter()
                     .enumerate()
                     .map(|(i, atom)| (bader_maxima[i] as usize, *atom))
                     .collect();
    let duration = start.elapsed();
    println!("assign atoms 2 thread: {:?}", duration);
    let pbar = Bar::visible(voxel_map.maxima_map.len() as u64,
                            100,
                            String::from("Assigning to Atoms: "));
    let start = Instant::now();
    let bader_maxima = voxel_map.maxima_map
                                .keys()
                                .map(|maxima| *maxima as isize)
                                .collect::<Vec<isize>>();
    let (assigned_atom, minimum_distance) =
        assign_maxima(&bader_maxima, &atoms, &voxel_map.grid, 4, pbar)?;
    let atom_map: FxHashMap<usize, usize> =
        assigned_atom.iter()
                     .enumerate()
                     .map(|(i, atom)| (bader_maxima[i] as usize, *atom))
                     .collect();
    let duration = start.elapsed();
    println!("assign atoms 4 thread: {:?}", duration);

    let pbar = Bar::visible(voxel_map.grid.size.total as u64,
                            100,
                            String::from("summing charge: "));
    let start = Instant::now();
    let (bader_charge, bader_volume, min_surf_dist) =
        sum_bader_densities(&densities, &voxel_map, &atoms, &atom_map, 1,
                            pbar)?;
    let duration = start.elapsed();
    println!("charge sum 1 threads: {:?}", duration);
    let pbar = Bar::visible(voxel_map.grid.size.total as u64,
                            100,
                            String::from("summing charge: "));
    let start = Instant::now();
    let (bader_charge, bader_volume, min_surf_dist) =
        sum_bader_densities(&densities, &voxel_map, &atoms, &atom_map, 2,
                            pbar)?;
    let duration = start.elapsed();
    println!("charge sum 2 threads: {:?}", duration);
    let pbar = Bar::visible(voxel_map.grid.size.total as u64,
                            100,
                            String::from("summing charge: "));
    let start = Instant::now();
    let (bader_charge, bader_volume, min_surf_dist) =
        sum_bader_densities(&densities, &voxel_map, &atoms, &atom_map, 4,
                            pbar)?;
    let duration = start.elapsed();
    println!("charge sum 4 threads: {:?}", duration);

    // END TESTING

    if let Err(e) = io::output::write_densities(&atoms,
                                                &analysis,
                                                densities,
                                                args.output,
                                                &voxel_map,
                                                &file_type)
    {
        panic!("Error occured: {}", e);
    }
    Ok(())
}
