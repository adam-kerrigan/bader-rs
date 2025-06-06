use bader::analysis::{
    bond_pruning, cage_pruning, calculate_bader_density, calculate_bader_error,
    calculate_bader_volumes_and_radii, nuclei_ordering, ring_pruning,
};
use bader::arguments::App;
use bader::errors::ArgumentError;
use bader::io::{self, FileFormat, FileType, WriteType};
use bader::methods::{maxima_finder, minima_finder, weight};
use bader::utils::vacuum_index;
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
    if !args.silent {
        let version = env!("CARGO_PKG_VERSION");
        let description = env!("CARGO_PKG_DESCRIPTION");
        println!("{}: v{}", description, version);
        println!("Running on {} threads.", args.threads);
    }
    // read the input files into a densities vector and a Grid struct
    let file_type: Box<dyn FileFormat> = match args.file_type {
        FileType::Vasp => Box::new(io::vasp::Vasp {}),
        FileType::Cube => Box::new(io::cube::Cube {}),
    };
    let (densities, rho, atoms, grid, voxel_origin) = file_type.init(&args);
    let reference = if rho.is_empty() { &densities[0] } else { &rho };
    let voxel_map =
        BlockingVoxelMap::new(grid, atoms.lattice.to_cartesian, voxel_origin);
    // create the index list which will tell us in which order to evaluate the
    // voxels
    let mut index: Vec<usize> = (0..voxel_map.grid.size.total).collect();
    index.sort_unstable_by(|a, b| {
        reference[*b].partial_cmp(&reference[*a]).unwrap()
    });
    // remove from the indices any voxel that is below the vacuum limit
    let vacuum_i = match vacuum_index(reference, &index, args.vacuum_tolerance)
    {
        Ok(i) => i,
        Err(e) => panic!("{}", e),
    };
    index.truncate(vacuum_i);
    // find the maxima in the system and store them whilst removing them from
    // the index list
    let nuclei = match maxima_finder(
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
            file_type.coordinate_format(e.maximum),
            e.atom + 1,
            e.distance,
        ),
    };
    // input the maxima as atoms into the voxel map
    nuclei.iter().for_each(|maximum| {
        voxel_map.maxima_store(maximum.position, maximum.atoms[0] as isize);
    });
    let n_bader_maxima = nuclei.len();
    let nuclei =
        nuclei_ordering(nuclei, reference, atoms.positions.len(), !args.silent);
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
    let (atoms_volume, atoms_radius) = calculate_bader_volumes_and_radii(
        &voxel_map,
        &atoms,
        args.threads,
        !args.silent,
    );
    let mut atoms_density =
        vec![vec![0.0; densities.len()]; atoms_volume.len()];
    densities.iter().enumerate().for_each(|(i, density)| {
        atoms_density
            .iter_mut()
            .zip(
                calculate_bader_density(
                    density,
                    &voxel_map,
                    &atoms,
                    args.threads,
                    !args.silent,
                )
                .iter(),
            )
            .for_each(|(ad, bd)| ad[i] += bd);
    });
    let atoms_error = calculate_bader_error(
        reference,
        &voxel_map,
        &atoms,
        args.threads,
        !args.silent,
    );
    let bonds =
        bond_pruning(&bonds, reference, voxel_map.grid_get(), !args.silent);
    let rings = ring_pruning(
        &rings,
        &nuclei,
        reference,
        &atoms,
        voxel_map.grid_get(),
        !args.silent,
    );
    let cages = cage_pruning(
        &cages,
        &nuclei,
        reference,
        &atoms,
        voxel_map.grid_get(),
        !args.silent,
    );
    println!(
        "{} {} {} {}",
        nuclei.len(),
        bonds.len(),
        rings.len(),
        cages.len()
    );
    /*
    critical_points.0.iter().for_each(|cp| {
        let [x, y, z] = voxel_map.grid.to_3d(cp.position);
        let x = x as f64 / voxel_map.grid.size.x as f64;
        let y = y as f64 / voxel_map.grid.size.y as f64;
        let z = z as f64 / voxel_map.grid.size.z as f64;
        let (x, y, z) = file_type.coordinate_format([x, y, z]);
        println!("{} {} {}        {:?}", x, y, z, cp.atoms);
    });
    critical_points.1.iter().for_each(|cp| {
        let [x, y, z] = voxel_map.grid.to_3d(cp.position);
        let x = x as f64 / voxel_map.grid.size.x as f64;
        let y = y as f64 / voxel_map.grid.size.y as f64;
        let z = z as f64 / voxel_map.grid.size.z as f64;
        let (x, y, z) = file_type.coordinate_format([x, y, z]);
        println!("{} {} {}        {:?}", x, y, z, cp.atoms);
    });
    critical_points.2.iter().for_each(|cp| {
        let [x, y, z] = voxel_map.grid.to_3d(cp.position);
        let x = x as f64 / voxel_map.grid.size.x as f64;
        let y = y as f64 / voxel_map.grid.size.y as f64;
        let z = z as f64 / voxel_map.grid.size.z as f64;
        let (x, y, z) = file_type.coordinate_format([x, y, z]);
        println!("{} {} {}        {:?}", x, y, z, cp.atoms);
    });
    critical_points.3.iter().for_each(|cp| {
        let [x, y, z] = voxel_map.grid.to_3d(cp.position);
        let x = x as f64 / voxel_map.grid.size.x as f64;
        let y = y as f64 / voxel_map.grid.size.y as f64;
        let z = z as f64 / voxel_map.grid.size.z as f64;
        let (x, y, z) = file_type.coordinate_format([x, y, z]);
        println!("{} {} {}        {:?}", x, y, z, cp.atoms);
    });
    */
    // prepare the positions for writing out
    let positions = atoms
        .positions
        .iter()
        .map(|coords| file_type.coordinate_format(*coords))
        .collect();
    // generate the output file
    let mut atoms_charge_file = io::output::partitions_file(
        positions,
        &atoms_density,
        &atoms_volume,
        &atoms_radius,
        &atoms_error,
    );
    atoms_charge_file.push_str(&format!(
        "\n  Bader Maxima: {}\n  Boundary Voxels: {}\n  Total Voxels: {}",
        n_bader_maxima,
        voxel_map.weight_len(),
        reference.len()
    ));
    // check that the write was successfull
    if io::output::write(atoms_charge_file, String::from("ACF.dat")).is_err() {
        panic!("Error in writing ACF.dat")
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
            if file_type
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
