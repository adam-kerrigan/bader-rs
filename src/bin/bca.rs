use bader::analysis::calculate_bader_density;
use bader::arguments::App;
use bader::critical::{
    CriticalPoints, bond_adjacency, bond_pruning, cage_pruning,
    critical_point_merge, nuclei_ordering, ring_pruning,
};
use bader::errors::ArgumentError;
use bader::io::output::{self, CriticalPointOutput, PartitionTable};
use bader::io::{FileFormat, WriteType};
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
            Err(
                e @ (ArgumentError::ShortHelp(_)
                | ArgumentError::LongHelp(_)
                | ArgumentError::NoFile(_)),
            ) => {
                println!("{}", e);
                return;
            }
            Err(e) => panic!("{}", e),
        };
    let visible = !args.silent;
    // print the splash
    let splash = format!(
        "{}: v{}\nRunning on {} threads.",
        env!("CARGO_PKG_DESCRIPTION"),
        env!("CARGO_PKG_VERSION"),
        args.threads
    );
    if visible {
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
    let index = index_generator(reference, args.vacuum_tolerance)
        .unwrap_or_else(|e| panic!("{}", e));
    // find the maxima in the system and assign them to atoms
    let mut nuclei = maxima_finder(
        &index,
        reference,
        &voxel_map,
        args.maximum_distance,
        &atoms,
        args.threads,
        visible,
    )
    .unwrap_or_else(|e| {
        panic!(
            "\nBader maximum at {:#?}\n is too far away from nearest atom: {} with a distance of {} Ang.",
            args.file_type.coordinate_format(e.maximum),
            e.atom + 1,
            e.distance,
        )
    });
    // pick the true maximum of each atom and store every maximum in the map
    let ordered_nuclei =
        nuclei_ordering(&mut nuclei, reference, &atoms, &voxel_map.grid);
    for maximum in &nuclei {
        voxel_map.maxima_store(maximum.position, maximum.atoms[0]);
    }
    // partition the density, collecting the saddle points on the way
    let (bonds, rings) = weight(
        reference,
        &voxel_map,
        &index,
        args.weight_tolerance,
        visible,
        args.threads,
    );
    // convert into a VoxelMap as the map is filled and no longer needs to block
    let voxel_map = VoxelMap::from_blocking_voxel_map(voxel_map);
    let cages =
        minima_finder(&index, reference, &voxel_map, args.threads, visible);
    // sum the densities
    let partition = calculate_bader_density(
        &densities,
        &voxel_map,
        &atoms,
        args.threads,
        visible,
    );
    // prune the critical points
    let bonds = bond_pruning(&bonds, reference, args.threads, visible);
    let adjacency = bond_adjacency(&bonds, atoms.positions.len());
    let critical_points = CriticalPoints {
        nuclei: ordered_nuclei,
        rings: critical_point_merge(ring_pruning(
            &rings,
            &adjacency,
            reference,
            args.threads,
            visible,
        )),
        cages: critical_point_merge(cage_pruning(
            &cages,
            reference,
            args.threads,
            visible,
        )),
        bonds,
    };
    // generate the output file
    let partition_table = PartitionTable::new(
        &partition,
        voxel_map.weight_map.len(),
        reference.len(),
    );
    let critical_point_output = CriticalPointOutput::new(
        &atoms,
        &critical_points,
        reference,
        &voxel_map.grid,
        args.file_type,
    );
    output::write(
        format!(
            "# {}\n\n{}\n\n{}",
            splash, partition_table, critical_point_output
        ),
        String::from("bader_charge_analysis.md"),
    )
    .unwrap_or_else(|e| {
        panic!("Error in writing bader_charge_analysis.md: {}", e)
    });
    // write any Bader volumes that have been requested
    let WriteType::Atom(atom_indices) = args.output else {
        return;
    };
    let atom_indices = if atom_indices.is_empty() {
        (0..atoms.positions.len()).collect()
    } else {
        atom_indices
    };
    let filenames = output::density_names(densities.len())
        .iter()
        .map(|name| name.to_lowercase().replace(' ', "_"))
        .collect::<Vec<String>>();
    for atom in atom_indices {
        let volume = voxel_map.volume_map(atom);
        for (density, name) in densities.iter().zip(&filenames) {
            let data =
                volume.iter().zip(density).map(|(w, rho)| w * rho).collect();
            args.file_type
                .write(&atoms, data, format!("{}_{}", atom + 1, name), visible)
                .unwrap_or_else(|e| panic!("Error in writing {}: {}", name, e));
        }
    }
}
