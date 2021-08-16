use bader::analysis::Analysis;
use bader::arguments::{Args, ClapApp};
use bader::io::{self, FileFormat, FileType};
use bader::progress::Bar;
use bader::utils::vacuum_index;
use bader::voxel_map::VoxelMap;

fn main() {
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
        VoxelMap::new(grid, atoms.lattice.to_cartesian, voxel_origin);
    let mut index: Vec<usize> = (0..voxel_map.grid.size.total).collect();
    index.sort_unstable_by(|a, b| {
             reference[*b].partial_cmp(&reference[*a]).unwrap()
         });
    let vacuum_index = vacuum_index(reference, &index, args.vacuum_tolerance);
    // Start a thread-safe progress bar and run the main calculation
    let pbar = Bar::visible(vacuum_index as u64,
                            100,
                            String::from("Bader Partitioning: "));
    voxel_map.calc(reference,
                   &index,
                   pbar,
                   args.threads,
                   vacuum_index,
                   args.weight_tolerance);
    let mut analysis =
        Analysis::new(&voxel_map, densities.len(), atoms.positions.len());
    // find the nearest atom to each Bader maxima
    let pbar = Bar::visible(analysis.bader_maxima.len() as u64,
                            100,
                            String::from("Assigning to Atoms: "));
    analysis.assign_atoms(&atoms, &voxel_map.grid, pbar);
    let pbar = Bar::visible(voxel_map.grid.size.total as u64,
                            100,
                            String::from("Summing Charge: "));
    // Sum the charge in each volume
    if let Err(e) = analysis.charge_sum(&atoms, &densities, &voxel_map, pbar) {
        println!("{}", e);
    }
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
    if let Err(e) = io::output::write_densities(&atoms,
                                                &analysis,
                                                densities,
                                                args.output,
                                                &voxel_map,
                                                &file_type)
    {
        panic!("Error occured: {}", e);
    }
}
