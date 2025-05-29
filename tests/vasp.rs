#[cfg(test)]
mod tests {
    use bader::io::FileFormat;
    use bader::io::vasp::Vasp;

    #[test]
    fn vasp_read_no_spin() {
        let filename = String::from("tests/vasp/CHGCAR_no_spin");
        let vasp = Vasp {};
        let (voxel_origin, grid, atoms, densities) = match vasp.read(filename) {
            Ok(r) => r,
            Err(e) => panic!("{}", e),
        };
        assert_eq!(voxel_origin, [0.; 3]);
        assert_eq!(grid, [32, 32, 32]);
        assert_eq!(atoms.positions, vec![[0., 0., 0.]]);
        assert_eq!(densities[0][0], 0.15246059033E+03 / atoms.lattice.volume);
        assert_eq!(
            densities[0][32767],
            0.13036296982E+03 / atoms.lattice.volume
        );
    }

    #[test]
    fn vasp_read_no_spin_chg() {
        let filename = String::from("tests/vasp/CHG_no_spin");
        let vasp = Vasp {};
        let (voxel_origin, grid, atoms, densities) = match vasp.read(filename) {
            Ok(r) => r,
            Err(e) => panic!("{}", e),
        };
        assert_eq!(voxel_origin, [0.; 3]);
        assert_eq!(grid, [32, 32, 32]);
        assert_eq!(atoms.positions, vec![[0., 0., 0.]]);
        assert_eq!(densities[0][0], 152.46 / atoms.lattice.volume);
        assert_eq!(densities[0][32767], 130.36 / atoms.lattice.volume);
    }

    #[test]
    fn vasp_read_spin() {
        let filename = String::from("tests/vasp/CHGCAR_spin");
        let vasp = Vasp {};
        let (voxel_origin, grid, atoms, densities) = match vasp.read(filename) {
            Ok(r) => r,
            Err(e) => panic!("{}", e),
        };
        assert_eq!(voxel_origin, [0.; 3]);
        assert_eq!(grid, [32, 32, 32]);
        assert_eq!(atoms.positions, vec![[0., 0., 0.]]);
        assert_eq!(densities[0][0], 0.15245934681E+03 / atoms.lattice.volume);
        assert_eq!(
            densities[0][32767],
            0.13036192086E+03 / atoms.lattice.volume
        );
        assert_eq!(densities[1][0], -0.10283642961E-07 / atoms.lattice.volume);
        assert_eq!(
            densities[1][32767],
            -0.87468511150E-08 / atoms.lattice.volume
        );
    }

    #[test]
    fn vasp_read_spin_chg() {
        let filename = String::from("tests/vasp/CHG_spin");
        let vasp = Vasp {};
        let (voxel_origin, grid, atoms, densities) = match vasp.read(filename) {
            Ok(r) => r,
            Err(e) => panic!("{}", e),
        };
        assert_eq!(voxel_origin, [0.; 3]);
        assert_eq!(grid, [32, 32, 32]);
        assert_eq!(atoms.positions, vec![[0., 0., 0.]]);
        assert_eq!(densities[0][0], 152.46 / atoms.lattice.volume);
        assert_eq!(densities[0][32767], 130.36 / atoms.lattice.volume);
        assert_eq!(densities[1][0], -0.10284E-07 / atoms.lattice.volume);
        assert_eq!(densities[1][32767], -0.87469E-08 / atoms.lattice.volume);
    }

    #[test]
    fn vasp_read_ncl() {
        let filename = String::from("tests/vasp/CHGCAR_ncl");
        let vasp = Vasp {};
        let (voxel_origin, grid, atoms, densities) = match vasp.read(filename) {
            Ok(r) => r,
            Err(e) => panic!("{}", e),
        };
        assert_eq!(voxel_origin, [0.; 3]);
        assert_eq!(grid, [32, 32, 32]);
        assert_eq!(atoms.positions, vec![[0., 0., 0.]]);
        assert_eq!(densities[0][0], 0.15229118148E+03 / atoms.lattice.volume);
        assert_eq!(
            densities[0][32767],
            0.13021559741E+03 / atoms.lattice.volume
        );
        assert_eq!(densities[1][0], -0.50501186231E-02 / atoms.lattice.volume);
        assert_eq!(
            densities[1][32767],
            -0.56304248048E-02 / atoms.lattice.volume
        );
        assert_eq!(densities[2][0], -0.89074011765E-03 / atoms.lattice.volume);
        assert_eq!(
            densities[2][32767],
            -0.95861710945E-03 / atoms.lattice.volume
        );
        assert_eq!(densities[3][0], 0.16139598297E+02 / atoms.lattice.volume);
        assert_eq!(
            densities[3][32767],
            0.13834498321E+02 / atoms.lattice.volume
        );
    }

    #[test]
    fn vasp_read_ncl_chg() {
        let filename = String::from("tests/vasp/CHG_ncl");
        let vasp = Vasp {};
        let (voxel_origin, grid, atoms, densities) = match vasp.read(filename) {
            Ok(r) => r,
            Err(e) => panic!("{}", e),
        };
        assert_eq!(voxel_origin, [0.; 3]);
        assert_eq!(grid, [32, 32, 32]);
        assert_eq!(atoms.positions, vec![[0., 0., 0.]]);
        assert_eq!(densities[0][0], 152.29 / atoms.lattice.volume);
        assert_eq!(densities[0][32767], 130.22 / atoms.lattice.volume);
        assert_eq!(densities[1][0], -0.50501E-02 / atoms.lattice.volume);
        assert_eq!(densities[1][32767], -0.56304E-02 / atoms.lattice.volume);
        assert_eq!(densities[2][0], -0.89074E-03 / atoms.lattice.volume);
        assert_eq!(densities[2][32767], -0.95862E-03 / atoms.lattice.volume);
        assert_eq!(densities[3][0], 16.140 / atoms.lattice.volume);
        assert_eq!(densities[3][32767], 13.834 / atoms.lattice.volume);
    }
}
