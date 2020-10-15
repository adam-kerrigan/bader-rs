#[cfg(test)]
mod tests {
    use bader::io::cube::Cube;
    use bader::io::FileFormat;

    const LENGTH_UNITS: f64 = 0.52917721067;
    const VOLUME_UNITS: f64 = LENGTH_UNITS * LENGTH_UNITS * LENGTH_UNITS;

    #[test]
    fn cube_read() {
        let filename = String::from("tests/cube/anatase.cube");
        let cube = Cube {};
        let (voxel_origin, grid, atoms, densities) = match cube.read(filename) {
            Ok(r) => r,
            Err(e) => panic!("{}", e),
        };
        assert_eq!(voxel_origin, [0.5; 3]);
        assert_eq!(grid, [96, 96, 180]);
        assert_eq!(atoms.positions.len(), 576);
        assert_eq!(densities[0][0], 0.13387E-02 / VOLUME_UNITS);
        assert_eq!(densities[0][1658879], 0.11782E+01 / VOLUME_UNITS);
    }
}
