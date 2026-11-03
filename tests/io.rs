//! Tests of reading and writing complete files through the `FileFormat` trait.
use bader::arguments::App;
use bader::io::cube::Cube;
use bader::io::vasp::Vasp;
use bader::io::{FileFormat, FileType};

/// A unique directory in the temporary directory for written files.
fn scratch(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bader-test-{}-{}",
        name,
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn vasp_write_round_trip() {
    let vasp = Vasp {};
    let (_, grid, atoms, densities) =
        vasp.read(String::from("tests/vasp/CHGCAR_spin")).unwrap();
    let dir = scratch("vasp");
    let filename = dir.join("round_trip").to_string_lossy().to_string();
    vasp.write(&atoms, densities[1].clone(), filename.clone(), false)
        .unwrap();
    let (_, grid_2, atoms_2, written) =
        vasp.read(format!("{}_CHGCAR", filename)).unwrap();
    assert_eq!(grid, grid_2);
    assert_eq!(atoms.positions, atoms_2.positions);
    assert_eq!(written.len(), 1);
    // the spin density is negative and small, it must survive the write
    for (a, b) in densities[1].iter().zip(&written[0]) {
        assert!((a - b).abs() <= 1e-10 * a.abs(), "{} != {}", a, b);
    }
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn cube_write_round_trip() {
    let cube = Cube {};
    let (origin, grid, atoms, densities) =
        cube.read(String::from("tests/cube/anatase.cube")).unwrap();
    let dir = scratch("cube");
    let filename = dir.join("round_trip").to_string_lossy().to_string();
    cube.write(&atoms, densities[0].clone(), filename.clone(), false)
        .unwrap();
    let (origin_2, grid_2, atoms_2, written) =
        cube.read(format!("{}.cube", filename)).unwrap();
    assert_eq!(origin, origin_2);
    assert_eq!(grid, grid_2);
    assert_eq!(atoms.positions, atoms_2.positions);
    for (a, b) in densities[0].iter().zip(&written[0]) {
        assert!((a - b).abs() <= 1e-4 * a.abs(), "{} != {}", a, b);
    }
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn init_reference_sums_files() {
    let app = App::new();
    let args = app
        .parse_args(vec![
            "bca",
            "tests/vasp/CHGCAR_no_spin",
            "-r",
            "tests/vasp/CHGCAR_no_spin",
            "-r",
            "tests/vasp/CHG_no_spin",
        ])
        .unwrap();
    let (densities, rho, _, grid, _) = args.file_type.init(&args);
    assert_eq!(grid, [32, 32, 32]);
    assert_eq!(densities.len(), 1);
    // the CHG is the CHGCAR to 5 significant figures
    for ((d, r), i) in densities[0].iter().zip(&rho).zip(0..) {
        assert!((2.0 * d - r).abs() < 1e-3 * d.abs(), "voxel {}", i);
    }
    let args = app
        .parse_args(vec![
            "bca",
            "tests/vasp/CHGCAR_no_spin",
            "-r",
            "tests/vasp/CHG_no_spin",
        ])
        .unwrap();
    let (_, rho_one, _, _, _) = args.file_type.init(&args);
    assert_eq!(rho_one.len(), rho.len());
    let args = app
        .parse_args(vec!["bca", "tests/vasp/CHGCAR_no_spin"])
        .unwrap();
    let (_, rho_none, _, _, _) = args.file_type.init(&args);
    assert!(rho_none.is_empty());
}

#[test]
fn init_spin_file_for_vasp() {
    let app = App::new();
    let args = app
        .parse_args(vec![
            "bca",
            "tests/vasp/CHG_no_spin",
            "-s",
            "tests/vasp/CHGCAR_no_spin",
        ])
        .unwrap();
    let (densities, _, _, _, _) = args.file_type.init(&args);
    assert_eq!(densities.len(), 2);
    assert_eq!(densities[1].len(), 32 * 32 * 32);
}

#[test]
#[should_panic(expected = "already has 1 spin densities")]
fn init_spin_file_rejects_spin_polarised() {
    let app = App::new();
    let args = app
        .parse_args(vec![
            "bca",
            "tests/vasp/CHG_spin",
            "-s",
            "tests/vasp/CHGCAR_no_spin",
        ])
        .unwrap();
    let _ = args.file_type.init(&args);
}

#[test]
fn file_type_dispatch() {
    assert_eq!(
        FileType::Vasp(Vasp {}).coordinate_format([1.0, 2.0, 3.0]),
        [3.0, 2.0, 1.0]
    );
    assert_eq!(
        FileType::Cube(Cube {}).coordinate_format([1.0, 2.0, 3.0]),
        [1.0, 2.0, 3.0]
    );
    assert!(
        FileType::Vasp(Vasp {})
            .read(String::from("missing"))
            .is_err()
    );
    assert!(
        FileType::Cube(Cube {})
            .read(String::from("tests/vasp/CHG_no_spin"))
            .is_err()
    );
}
