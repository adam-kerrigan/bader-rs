//! Runs the whole analysis through the library on a small VASP file and checks
//! the totals, matching what the `bca` binary does.
use bader::analysis::calculate_bader_density;
use bader::critical::{
    bond_adjacency, bond_pruning, cage_pruning, critical_point_merge,
    nuclei_ordering, ring_pruning,
};
use bader::io::FileFormat;
use bader::io::vasp::Vasp;
use bader::methods::{maxima_finder, minima_finder, weight};
use bader::utils::index_generator;
use bader::voxel_map::{BlockingVoxelMap, VoxelMap};

#[test]
fn fcc_nickel_partition() {
    let (voxel_origin, grid, atoms, densities) = Vasp {}
        .read(String::from("tests/vasp/CHGCAR_spin"))
        .unwrap();
    let reference = &densities[0];
    let threads = 3;
    let voxel_map =
        BlockingVoxelMap::new(grid, atoms.lattice.to_cartesian, voxel_origin);
    let index = index_generator(reference, 1e-6).unwrap();
    assert_eq!(index.len(), reference.len(), "no vacuum in a metal");
    let mut nuclei = maxima_finder(
        &index, reference, &voxel_map, 1.0, &atoms, threads, false,
    )
    .unwrap();
    assert_eq!(nuclei.len(), 1);
    assert_eq!(nuclei[0].position, 0, "the maximum is on the atom");
    let ordered =
        nuclei_ordering(&mut nuclei, reference, &atoms, &voxel_map.grid);
    assert_eq!(ordered.len(), 1);
    for cp in &nuclei {
        voxel_map.maxima_store(cp.position, cp.atoms[0]);
    }
    let (bonds, rings) =
        weight(reference, &voxel_map, &index, 1e-6, false, threads);
    let voxel_map = VoxelMap::from_blocking_voxel_map(voxel_map);
    let cages = minima_finder(&index, reference, &voxel_map, threads, false);
    let partition =
        calculate_bader_density(&densities, &voxel_map, &atoms, threads, false);
    // one atom so it has every electron and the whole cell
    assert!(
        (partition.charge[0][0] - 10.0).abs() < 1e-4,
        "{:?}",
        partition.charge
    );
    assert!(partition.charge[0][1].abs() < 1e-4, "no net spin");
    assert!((partition.volume[0] - atoms.lattice.volume).abs() < 1e-6);
    assert_eq!(partition.charge[1], vec![0.0, 0.0], "no vacuum charge");
    assert_eq!(partition.volume[1], 0.0);
    assert!(partition.radius[0] > 1.0 && partition.radius[0] < 1.3);
    assert!(partition.error[0].abs() < 1e-2);
    // the boundary weights are f32 so the total is only approximate
    let cell =
        voxel_map.voxel_map.len() as f64 * voxel_map.grid.voxel_lattice.volume;
    assert!((partition.volume.iter().sum::<f64>() - cell).abs() < 1e-5);
    // the atom is bonded to its own images across the periodic boundary
    let bonds = bond_pruning(&bonds, reference, threads, false);
    assert!(!bonds.is_empty());
    assert!(bonds.iter().all(|cp| cp.atoms.len() == 2
        && cp.atoms.iter().all(|a| a.atom_index() == 0)));
    let adjacency = bond_adjacency(&bonds, 1);
    assert_eq!(adjacency[0].len(), 2 * bonds.len());
    let rings = critical_point_merge(ring_pruning(
        &rings, &adjacency, reference, threads, false,
    ));
    let cages =
        critical_point_merge(cage_pruning(&cages, reference, threads, false));
    assert!(rings.iter().all(|cp| cp.atoms.len() >= 3));
    assert!(cages.iter().all(|cp| cp.atoms.len() >= 2));
    // every atom in a cage is the nickel atom in some image
    assert!(
        cages
            .iter()
            .all(|cp| cp.atoms.iter().all(|a| a.atom_index() == 0))
    );
}
