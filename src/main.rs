use bader::arguments::{Args, ClapApp, Reference};
use bader::density::{Density, VoxelMap};
use bader::io::{self, ReadFunction};
use bader::methods::StepMethod;
use bader::progress::Bar;
use indicatif::ProgressBar;
use rayon::prelude::*;

fn main() {
    let app = ClapApp::App;
    let args = Args::new(app.get().get_matches());

    rayon::ThreadPoolBuilder::new().num_threads(args.threads)
                                   .build_global()
                                   .unwrap();

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

    let reference = match args.reference.clone() {
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

    let method: StepMethod = args.method;

    let pbar = ProgressBar::new(reference.size.total as u64);
    let pbar = Bar::new(pbar, 100, String::from("Bader Calculation:"));
    let map = (0..reference.size.total).into_par_iter()
                                       .map(|i| {
                                           pbar.tick();
                                           method(i, &reference)
                                       })
                                       .collect::<Vec<isize>>();
    let mut bader_maxima = map.clone();
    bader_maxima.par_sort();
    bader_maxima.dedup();
    let voxel_map = VoxelMap::new(map, bader_maxima);
    let (assigned_atom, assigned_distance) =
        atoms.assign_maxima(&voxel_map, &reference);
    let surface_distance =
        voxel_map.surface_distance(&assigned_atom, &atoms, &reference);
    let (bader_charge, bader_volume, vacuum_charge, vacuum_volume) =
        { voxel_map.charge_sum(&densities) };
    let results = io::Results::new(voxel_map,
                                   bader_charge,
                                   bader_volume,
                                   assigned_atom,
                                   assigned_distance,
                                   surface_distance,
                                   vacuum_charge,
                                   vacuum_volume,
                                   reference,
                                   atoms,
                                   args.zyx_format);
    match results.write() {
        Ok(_) => {}
        Err(e) => println!("{}", e),
    }
}
