use bader::analysis::{calculate_bader_density, calculate_bader_radius};
use bader::arguments::App;
use bader::critical::{
    bond_adjacency, bond_pruning, cage_pruning, critical_point_merge,
    nuclei_ordering, ring_pruning,
};
use bader::errors::ArgumentError;
use bader::io::{self, FileFormat, WriteType};
use bader::methods::{maxima_finder, minima_finder, weight};
use bader::utils::index_generator;
use bader::voxel_map::{BlockingVoxelMap, VoxelMap};

fn main() {
    // argument parsing
    let app = App::new();
    let env_args = std::env::args().collect::<Vec<String>>();
    let args =
        match app.parse_args(env_args.iter().map(|s| s.as_str()).collect()) {
            Ok(a) => a,
            Err(e) => match e {
                ArgumentError::ShortHelp(_)
                | ArgumentError::LongHelp(_)
                | ArgumentError::NoFile(_) => {
                    println!("{}", e);
                    return;
                }
                _ => panic!("{}", e),
            },
        };
    // print the splash
    let splash = format!(
        "{}: v{}\nRunning on {} threads.",
        env!("CARGO_PKG_DESCRIPTION"),
        env!("CARGO_PKG_VERSION"),
        args.threads
    );
    if !args.silent {
        println!("{}", splash);
    }
    // read the input files into a densities vector and a Grid struct
    let (densities, rho, atoms, grid, voxel_origin) =
        args.file_type.init(&args);
    let reference = if rho.is_empty() { &densities[0] } else { &rho };
    let voxel_map =
        BlockingVoxelMap::new(grid, atoms.lattice.to_cartesian, voxel_origin);
    // create the index list which will tell us in which order to evaluate the
    // voxels and remove from the indices any voxel that is below the vacuum limit
    let index = match index_generator(reference, args.vacuum_tolerance) {
        Ok(i) => i,
        Err(e) => panic!("{}", e),
    };
    // find the maxima in the system and store them whilst removing them from
    // the index list
    let mut nuclei = match maxima_finder(
        &index,
        reference,
        &voxel_map,
        &args.maximum_distance,
        &atoms,
        args.threads,
        !args.silent,
    ) {
        Ok(v) => v,
        Err(e) => panic!(
            "\nBader maximum at {:#?}\n is too far away from nearest atom: {} with a distance of {} Ang.",
            args.file_type.coordinate_format(e.maximum),
            e.atom + 1,
            e.distance,
        ),
    };
    // input the maxima as atoms into the voxel map
    let ordered_nuclei =
        nuclei_ordering(&mut nuclei, &atoms, &voxel_map.grid, !args.silent);
    nuclei.iter().for_each(|maximum| {
        voxel_map.maxima_store(maximum.position, maximum.atoms[0].0 as isize);
    });
    // calculate the weights leave the critical points for now
    let (bonds, rings) = weight(
        reference,
        &voxel_map,
        &index,
        args.weight_tolerance,
        !args.silent,
        args.threads,
    );
    // convert into a VoxelMap as the map is filled and no longer needs to block
    let voxel_map = VoxelMap::from_blocking_voxel_map(voxel_map);
    // Find the minima
    let cages = minima_finder(
        &index,
        reference,
        &voxel_map,
        args.threads,
        !args.silent,
    );
    // sum the densities and then write the charge partition files
    let (atoms_density, atoms_volume, atoms_error) = calculate_bader_density(
        &densities,
        &voxel_map,
        &atoms,
        args.threads,
        !args.silent,
    );
    let bonds = bond_pruning(&bonds, &args);
    let bond_adjacency = bond_adjacency(&bonds, atoms.positions.len());
    let rings =
        critical_point_merge(ring_pruning(&rings, &bond_adjacency, &args));
    let cages = critical_point_merge(cage_pruning(
        &cages,
        &ordered_nuclei,
        &atoms,
        voxel_map.grid_get(),
        &args,
    ));
    let atoms_radius =
        calculate_bader_radius(&bonds, &atoms, &voxel_map.grid, !args.silent);
    // generate the output file
    let partition_table = io::output::PartitionTable::new(
        &atoms_density,
        &atoms_volume,
        &atoms_radius,
        &atoms_error,
        voxel_map.weight_len(),
        reference.len(),
    );
    let critical_points = (ordered_nuclei, bonds, rings, cages);
    let critical_point_output = io::output::CriticalPointOutput::new(
        atoms.positions.clone(),
        critical_points,
        &voxel_map.grid,
        args.file_type,
    );
    // check that the write was successfull
    if io::output::write(
        format!(
            "# {}\n\n{}\n\n{}",
            splash, partition_table, critical_point_output
        ),
        String::from("bader_charge_analysis.md"),
    )
    .is_err()
    {
        panic!("Error in writing bader_charge_analysis.md")
    }
    // let bonds_file = io::output::bonds_file(&bonds);
    // if io::output::write(bonds_file, String::from("BF.dat")).is_err() {
    //     panic!("Error in writing BF.dat")
    // }
    // Prepare to write any densities that have been requested.
    let filename = match densities.len().cmp(&2) {
        std::cmp::Ordering::Less => vec![String::from("charge")],
        std::cmp::Ordering::Equal => {
            vec![String::from("charge"), String::from("spin")]
        }
        std::cmp::Ordering::Greater => vec![
            String::from("charge"),
            String::from("spin_x"),
            String::from("spin_y"),
            String::from("spin_z"),
        ],
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
            Box::new(atom_iter.into_iter().map(|atom_number| {
                let map = voxel_map.volume_map(atom_number);
                (atom_number, map)
            }))
        } else {
            Box::new(Vec::with_capacity(0).into_iter())
        }
    };
    write_map.for_each(|(id, weight_map)| {
        densities.iter().zip(&filename).for_each(|(rho, flnm)| {
            if args
                .file_type
                .write(
                    &atoms,
                    weight_map
                        .iter()
                        .zip(rho)
                        .map(|(weight, charge)| weight.map(|w| w * charge))
                        .collect(),
                    format!("{}_{}", id + 1, flnm),
                    !args.silent,
                )
                .is_err()
            {
                panic!("Error in writing {}", flnm)
            }
        });
    });
}
