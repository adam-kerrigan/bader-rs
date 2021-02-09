use atomic_counter::{AtomicCounter, RelaxedCounter};
use bader::analysis::Analysis;
use bader::arguments::{Args, ClapApp};
use bader::grid::Grid;
use bader::io::{self, FileFormat, FileType};
use bader::methods::weight;
use bader::progress::Bar;
use bader::utils::vacuum_tolerance;
use bader::voxel_map::VoxelMap;
use crossbeam_utils::thread;

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
    let grid = Grid::new(grid,
                         atoms.lattice.to_cartesian,
                         args.weight_tolerance,
                         args.maxima_tolerance,
                         args.vacuum_tolerance,
                         voxel_origin);
    let voxel_map = VoxelMap::new(grid.size.total);
    {
        let mut index: Vec<usize> = (0..grid.size.total).collect();
        // Start a thread-safe progress bar and run the main calculation
        println!("Sorting density.");
        index.sort_unstable_by(|a, b| {
                 reference[*b].partial_cmp(&reference[*a]).unwrap()
             });
        let counter = RelaxedCounter::new(0);
        let vacuum_index =
            vacuum_tolerance(&reference, &index, grid.vacuum_tolerance);
        let pbar = Bar::visible(vacuum_index as u64,
                                100,
                                String::from("Bader Partitioning: "));
        thread::scope(|s| {
            for _ in 0..args.threads {
                s.spawn(|_| loop {
                     let p = {
                         let i = counter.inc();
                         if i >= vacuum_index {
                             break;
                         };
                         index[i]
                     };
                     weight(p, &grid, &reference, &voxel_map);
                     pbar.tick();
                 });
            }
        }).unwrap();
    }
    {
        let mut weights = voxel_map.lock();
        weights.shrink_to_fit();
    }
    let mut analysis =
        Analysis::new(&voxel_map, densities.len(), atoms.positions.len());
    // find the nearest atom to each Bader maxima
    let pbar = Bar::visible(analysis.bader_maxima.len() as u64,
                            100,
                            String::from("Assigning to Atoms: "));
    analysis.assign_atoms(&atoms, &grid, pbar);
    let pbar = Bar::visible(grid.size.total as u64,
                            100,
                            String::from("Summing Charge: "));
    // Sum the charge in each volume
    if let Err(e) =
        analysis.charge_sum(&atoms, &densities, &grid, &voxel_map, pbar)
    {
        println!("{}", e);
    }
    // Sum the charge in each atom
    analysis.atoms_charge_sum();
    // build the results
    println!("Writing output files:");
    let (atoms_charge_file, bader_charge_file) =
        io::output::charge_files(&analysis, &atoms, &grid, &file_type);
    // check that the write was successfull
    if let Err(e) = io::output::write(atoms_charge_file, bader_charge_file) {
        panic!("Error occured: {}", e);
    }
    println!("ACF.dat and BCF.dat written successfully.");
    if let Err(e) = io::output::write_densities(&atoms,
                                                &analysis,
                                                densities,
                                                &grid,
                                                args.output,
                                                &voxel_map,
                                                &file_type)
    {
        panic!("Error occured: {}", e);
    }
}
