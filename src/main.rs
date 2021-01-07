use atomic_counter::{AtomicCounter, RelaxedCounter};
use bader::arguments::{Args, ClapApp};
use bader::density::Density;
use bader::io::{self, FileFormat, FileType};
use bader::methods::{self, Method};
use bader::progress::Bar;
use bader::utils::vacuum_tolerance;
use bader::voxel_map::VoxelMap;
use rayon::prelude::*;
use std::sync::Arc;

fn main() {
    // argument parsing
    let app = ClapApp::get();
    let args = Args::new(app.get_matches());
    // set up the threads
    rayon::ThreadPoolBuilder::new().num_threads(args.threads)
                                   .build_global()
                                   .unwrap();
    // print splash
    println!("Multi-threaded Bader Charge Analysis ({})",
             env!("CARGO_PKG_VERSION"));
    // read the input files into a densities vector and a Density struct
    let file_type: Box<dyn FileFormat> = match args.file_type {
        FileType::Vasp => Box::new(io::vasp::Vasp {}),
        FileType::Cube => Box::new(io::cube::Cube {}),
    };

    let (densities, rho, atoms, grid, voxel_origin) = file_type.init(&args);
    let reference = if rho.is_empty() {
        Density::new(&densities[0],
                     grid,
                     atoms.lattice.to_cartesian,
                     args.weight_tolerance,
                     args.vacuum_tolerance,
                     voxel_origin)
    } else {
        Density::new(&rho,
                     grid,
                     atoms.lattice.to_cartesian,
                     args.weight_tolerance,
                     args.vacuum_tolerance,
                     voxel_origin)
    };
    let mut index: Vec<usize> = (0..reference.size.total).collect();
    // Start a thread-safe progress bar and run the main calculation
    let method: methods::StepMethod = match args.method {
        Method::OnGrid => {
            println!("Sorting density.");
            index.par_sort_unstable_by(|a, b| {
                     reference.data[*b].partial_cmp(&reference.data[*a])
                                       .unwrap()
                 });
            methods::ongrid
        }
        Method::Weight => {
            println!("Sorting density.");
            index.par_sort_unstable_by(|a, b| {
                     reference.data[*b].partial_cmp(&reference.data[*a])
                                       .unwrap()
                 });
            methods::weight
        }
        Method::NearGrid => methods::neargrid,
    };
    let voxel_map = VoxelMap::new(reference.size.total);
    let voxel_map = Arc::new(voxel_map);
    let len = match args.method {
        Method::NearGrid => reference.size.total,
        _ => vacuum_tolerance(&reference, &index),
    };
    {
        let counter = RelaxedCounter::new(0);
        let pbar =
            Bar::new(len as u64, 100, String::from("Bader Partitioning: "));
        pbar.display();
        (0..len).into_par_iter().for_each(|_| {
                                    let p = {
                                        let i = counter.inc();
                                        index[i]
                                    };
                                    let (maxima, weights) =
                                        method(p,
                                               &reference,
                                               Arc::clone(&voxel_map));
                                    if !weights.is_empty() {
                                        let i = {
                                            let mut weight = voxel_map.lock();
                                            let i = weight.len();
                                            (*weight).push(weights);
                                            i
                                        };
                                        voxel_map.weight_store(p as isize, i);
                                    }
                                    voxel_map.maxima_store(p as isize, maxima);
                                    pbar.tick();
                                });
    }
    // find the nearest atom to each Bader maxima
    let mut voxel_map = match Arc::try_unwrap(voxel_map) {
        Ok(voxel_map) => voxel_map,
        _ => panic!(),
    };
    voxel_map.collect_maxima();
    let pbar = Bar::new(voxel_map.bader_maxima.len() as u64,
                        100,
                        String::from("Assigning to Atoms: "));
    pbar.display();
    voxel_map.assign_atoms(&atoms, &reference, pbar);
    let pbar = Bar::new(reference.size.total as u64,
                        100,
                        String::from("Summing Charge: "));
    pbar.display();
    voxel_map.charge_sum(&atoms, &densities, &reference, pbar);
    // build the results
    println!("Writing output files:");
    let (atoms_charge_file, bader_charge_file) =
        { file_type.results(voxel_map, atoms, &reference) };
    // check that the write was successfull
    match io::write(atoms_charge_file, bader_charge_file) {
        Ok(_) => {}
        Err(e) => println!("{}", e),
    }
    println!("ACF.dat and BCF.dat written successfully.");
}
