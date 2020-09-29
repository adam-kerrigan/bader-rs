use bader::arguments::{Args, ClapApp, Method, Reference, Weight};
use bader::density::Density;
use bader::voxel_map::VoxelMap;
use bader::io::{self, ReadFunction};
use bader::methods::{self, StepMethod};
use bader::progress::Bar;
use indicatif::ProgressBar;
use rayon::prelude::*;
use std::collections::BTreeSet;
use std::sync::atomic::AtomicIsize;
use std::sync::Arc;

fn main() {
    // argument parsing
    let app = ClapApp::App.get();
    let args = Args::new(app.get_matches());
    // set up the threads
    rayon::ThreadPoolBuilder::new().num_threads(args.threads)
                                   .build_global()
                                   .unwrap();
    // print splash
    println!("Multi-threaded Bader Charge Analysis ({})",
             env!("CARGO_PKG_VERSION"));
    // read the input files into a densities vector and a Density struct
    let read_function: ReadFunction = args.read;
    let (voxel_origin, grid, atoms, densities) = match read_function(args.file)
    {
        Ok(r) => r,
        Err(e) => panic!("{}", e),
    };
    let reference_d = match args.reference.clone() {
        Reference::None => vec![],
        Reference::One(f) => {
            let (_, g, _, densities) = match read_function(f) {
                Ok(r) => r,
                Err(e) => panic!("{}", e),
            };
            assert_eq!(g, grid);
            densities
        }
        Reference::Two(f1, f2) => {
            let (_, g, _, densities) = match read_function(f1) {
                Ok(r) => r,
                Err(e) => panic!("{}", e),
            };
            assert_eq!(g, grid);
            let (_, g2, _, densities2) = match read_function(f2) {
                Ok(r) => r,
                Err(e) => panic!("{}", e),
            };

            assert_eq!(g, g2);
            let d = densities[0].par_iter()
                                .zip(&densities2[0])
                                .map(|(a, b)| a + b)
                                .collect::<Vec<f64>>();
            vec![d]
        }
    };
    let reference = match args.reference {
        Reference::None => Density::new(&densities[0],
                                        grid,
                                        atoms.lattice.to_cartesian,
                                        args.vacuum_tolerance,
                                        voxel_origin),
        _ => Density::new(&reference_d[0],
                          grid,
                          atoms.lattice.to_cartesian,
                          args.vacuum_tolerance,
                          voxel_origin),
    };
    let mut index: Vec<usize> = (0..reference.size.total).collect();
    let method: StepMethod = match args.method {
        Method::OnGrid => methods::ongrid,
        Method::Weight => methods::weight,
        Method::NearGrid => methods::neargrid,
    };
    match args.weight {
        Weight::None => (),
        _ => match args.method {
            Method::NearGrid => (),
            _ => {
                println!("Sorting the density.");
                index.par_sort_unstable_by(|a, b| {
                         reference.data[*b].partial_cmp(&reference.data[*a])
                                           .unwrap()
                     });
            }
        },
    }
    // Start a thread-safe progress bar and run the main calculation
    let pbar = ProgressBar::new(reference.size.total as u64);
    let pbar = Bar::new(pbar, 100, String::from("Bader Calculation: "));
    let mut map = Vec::new();
    map.resize_with(reference.size.total, || AtomicIsize::new(-1));
    let map = Arc::new(map);
    let bader_maxima = index.par_iter()
                            .map(|p| {
                                pbar.tick();
                                method(*p, &reference, Arc::clone(&map))
                            })
                            .collect::<BTreeSet<isize>>();
    // each maxima is a unique value in the map with vacuum being -1
    drop(pbar);
    let map = match Arc::try_unwrap(map) {
        Ok(map) => map.into_par_iter()
                      .map(|bv| bv.into_inner())
                      .collect::<Vec<isize>>(),
        Err(_) => panic!(), // I don't know why this would panic
    };
    let voxel_map = VoxelMap::new(map, bader_maxima.into_iter().collect());
    // find the nearest atom to each Bader maxima
    println!("Assigning maxima to atoms.");
    let (assigned_atom, assigned_distance) =
        atoms.assign_maxima(&voxel_map, &reference);
    // find the nearest point to the atom and sum the atoms charge
    match args.weight {
        Weight::None => (),
        _ => {
            if let Method::NearGrid = args.method {
                println!("Sorting the density.");
                index.par_sort_unstable_by(|a, b| {
                         reference.data[*b].partial_cmp(&reference.data[*a])
                                           .unwrap()
                     })
            }
        }
    }
    let (bader_charge, bader_volume, surface_distance) = {
        voxel_map.charge_sum(&densities,
                             &assigned_atom,
                             &atoms,
                             &reference,
                             index,
                             args.weight)
    };
    // build the results
    println!("Writing output files:");
    let results = io::Results::new(voxel_map,
                                   bader_charge,
                                   bader_volume,
                                   assigned_atom,
                                   assigned_distance,
                                   surface_distance,
                                   reference,
                                   atoms,
                                   args.zyx_format);
    // check that the write was successfull
    match results.write() {
        Ok(_) => {}
        Err(e) => println!("{}", e),
    }
    println!("ACF.dat and BCF.dat written successfully.");
}
