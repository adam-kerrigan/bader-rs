use crate::atoms::Lattice;
use crate::utils::dot;
use crate::voronoi::Voronoi;

/// Structure for managing the movement within the reference density.
pub struct Grid {
    /// The [`Shift`] structure for movement.
    shift: Shift,
    /// The 3d size of the data.
    pub size: Size,
    /// The voronoi vectors and their alphas.
    pub voronoi: Voronoi,
    /// Information on the voxel basis.
    pub voxel_lattice: Lattice,
    /// The origin of each voxel.
    pub voxel_origin: [f64; 3],
}

impl Grid {
    /// Initialises a grid structure. Computes the voxel_lattice from the grid and lattice.
    pub fn new(
        grid: [usize; 3],
        lattice: [[f64; 3]; 3],
        voxel_origin: [f64; 3],
    ) -> Self {
        let size = Size::new(grid[0], grid[1], grid[2]);
        let shift = Shift::new(&size);
        let voxel_lattice = Lattice::new([
            [
                lattice[0][0] / grid[0] as f64,
                lattice[0][1] / grid[0] as f64,
                lattice[0][2] / grid[0] as f64,
            ],
            [
                lattice[1][0] / grid[1] as f64,
                lattice[1][1] / grid[1] as f64,
                lattice[1][2] / grid[1] as f64,
            ],
            [
                lattice[2][0] / grid[2] as f64,
                lattice[2][1] / grid[2] as f64,
                lattice[2][2] / grid[2] as f64,
            ],
        ]);
        let voronoi = Voronoi::new(&voxel_lattice);
        Self {
            shift,
            size,
            voronoi,
            voxel_lattice,
            voxel_origin,
        }
    }

    /// get the full shift to visit the surrounding 26 voxels
    pub fn full_shift(&self, p: isize) -> [isize; 26] {
        let (shift, _) = self.shift.get(p);
        shift
            .iter()
            .take(13)
            .chain(shift.iter().skip(14))
            .copied()
            .collect::<Vec<isize>>()
            .try_into()
            .unwrap() // safe to unwrap as we know that we will get 26 items
    }

    /// get the reduced shift to visit the surrounding 6 voxels
    pub fn reduced_shift(&self, p: isize) -> [isize; 6] {
        let (shift, _) = self.shift.get(p);
        // [ +x, -x, +y, -y, +z, -z ]
        [
            shift[22], shift[4], shift[16], shift[10], shift[14], shift[12],
        ]
    }

    /// Get the gradient shift from the current point
    pub fn gradient_shift(&self, p: isize, gradient: [f64; 3]) -> isize {
        let (shift, _) = self.shift.get(p);
        // [ [-x, 0, +x], [-y, 0, +y], [-z, 0, +z] ]
        let i = (gradient[0] * 9f64 + gradient[1] * 3f64 + gradient[2] + 13f64)
            as usize;
        shift[i]
    }

    /// Shifts a point, p, by a single voronoi vector.
    pub fn voronoi_shifts(&self, p: isize) -> Vec<((isize, [i8; 3]), f64)> {
        self.voronoi
            .vectors
            .iter()
            .map(|shifts| {
                shifts.iter().fold((p, [0; 3]), |(pn, pbc), p_shift| {
                    let (shift, pbc_t) = self.shift.get(pn);
                    let pt = pn + shift[*p_shift];
                    let [x, y, z] = pbc_t[*p_shift];
                    (pt, [pbc[0] + x, pbc[1] + y, pbc[2] + z])
                })
            })
            .zip(self.voronoi.alphas.iter().copied())
            .collect()
    }

    pub fn voronoi_shifts_nocheck(&self, p: isize) -> Vec<(isize, f64)> {
        self.voronoi
            .vectors
            .iter()
            .map(|shifts| {
                shifts
                    .iter()
                    .fold(p, |pn, p_shift| pn + self.shift.get(pn).0[*p_shift])
            })
            .zip(self.voronoi.alphas.iter().copied())
            .collect()
    }

    /// Converts a 1D index of the array into a 3D index
    pub fn to_3d(&self, p: isize) -> [isize; 3] {
        let x = p / (self.size.y * self.size.z);
        let y = (p / self.size.z).rem_euclid(self.size.y);
        let z = p.rem_euclid(self.size.z);
        [x, y, z]
    }

    /// Converts a point in the array to cartesian.
    pub fn to_cartesian(&self, p: isize) -> [f64; 3] {
        let [x, y, z] = self.to_3d(p);
        let p = [
            x as f64 + self.voxel_origin[0],
            y as f64 + self.voxel_origin[1],
            z as f64 + self.voxel_origin[2],
        ];
        dot(p, self.voxel_lattice.to_cartesian)
    }

    pub fn encode_maxima(&self, maxima: usize, image: [i8; 3]) -> usize {
        let (decoded_maxima, decoded_image) = self.decode_maxima(maxima);
        let [x, y, z] = if decoded_maxima != maxima {
            decoded_image
                .into_iter()
                .zip(image)
                .map(|(di, i)| (di + i) as u8)
                .collect::<Vec<u8>>()
                .try_into()
                .unwrap()
        } else {
            image
                .into_iter()
                .map(|i| i as u8)
                .collect::<Vec<u8>>()
                .try_into()
                .unwrap()
        };
        // starting with 2 zeros allows for 32,767 atoms
        // TODO: Add a check at the start for this?
        let maxima_shift = usize::from_le_bytes([0, 0, x, y, z, 0, 0, 0]);
        decoded_maxima + maxima_shift
    }

    pub fn decode_maxima(&self, maxima: usize) -> (usize, [i8; 3]) {
        // the true maxima number has been shifted but by an integer multiple of the number of
        // atoms
        let encoded_maxima = maxima.to_le_bytes();
        let maxima = usize::from_le_bytes([
            encoded_maxima[0],
            encoded_maxima[1],
            0,
            0,
            0,
            0,
            0,
            0,
        ]);
        let maxima_bytes: Vec<i8> = encoded_maxima
            .into_iter()
            // the shift is offset by 2 to allow for 32,767 atoms
            .skip(2)
            .take(3)
            .map(|u| u as i8)
            .collect();
        (maxima, maxima_bytes.try_into().unwrap())
    }
}

/// Structure for holding the periodic shifts
struct Shift {
    di: [[isize; 27]; 27],
    pbc: [[[i8; 3]; 27]; 27],
    index: Vec<u8>,
}

impl Shift {
    /// Gets the periodic boundary shift
    fn get(&self, i: isize) -> ([isize; 27], [[i8; 3]; 27]) {
        let ii = self.index[i as usize] as usize;
        (self.di[ii], self.pbc[ii])
    }

    /// Generate the indices for Shift.index
    fn index_gen(size: &Size) -> Vec<u8> {
        let mut index = vec![0u8; size.total];
        let mut i: usize = 0;
        let bound_x: usize = (size.x - 1) as usize;
        let bound_y: usize = (size.y - 1) as usize;
        let bound_z: usize = (size.z - 1) as usize;

        // first x
        // first y
        index[i] = 26;
        i += 1;
        for _ in 1..bound_z {
            index[i] = 24;
            i += 1;
        }
        index[i] = 25;
        i += 1;
        // bulk y
        for _ in 1..bound_y {
            index[i] = 20;
            i += 1;
            for _ in 1..bound_z {
                index[i] = 18;
                i += 1;
            }
            index[i] = 19;
            i += 1;
        }
        // last y
        index[i] = 23;
        i += 1;
        for _ in 1..bound_z {
            index[i] = 21;
            i += 1;
        }
        index[i] = 22;
        i += 1;

        // bulk x
        // first y
        for _ in 1..bound_x {
            index[i] = 8;
            i += 1;
            for _ in 1..bound_z {
                index[i] = 6;
                i += 1;
            }
            index[i] = 7;
            i += 1;
            // bulk y
            for _ in 1..bound_y {
                index[i] = 2;
                i += bound_z;
                index[i] = 1;
                i += 1;
            }
            // last y
            index[i] = 5;
            i += 1;
            for _ in 1..bound_z {
                index[i] = 3;
                i += 1;
            }
            index[i] = 4;
            i += 1;
        }

        // last x
        // first y
        index[i] = 17;
        i += 1;
        for _ in 1..bound_z {
            index[i] = 15;
            i += 1;
        }
        index[i] = 16;
        i += 1;
        // bulk y
        for _ in 1..bound_y {
            index[i] = 11;
            i += 1;
            for _ in 1..bound_z {
                index[i] = 9;
                i += 1;
            }
            index[i] = 10;
            i += 1;
        }
        // last y
        index[i] = 14;
        i += 1;
        for _ in 1..bound_z {
            index[i] = 12;
            i += 1;
        }
        index[i] = 13;
        index
    }

    /// Creates the shift matrix for moving in the array
    /// the outer index in di is unraveling the index with a periodic roll so the point
    /// is at (0, 0, 0) -> 0
    /// <pre class="rust">
    /// 0 -> (0, 0, 0)   7 -> (0,-1, 1)  14 -> (1, 1,-1)   21 -> (-1, 1, 0)
    /// 1 -> (0, 0, 1)   8 -> (0,-1,-1)  15 -> (1,-1, 0)   22 -> (-1, 1, 1)
    /// 2 -> (0, 0,-1)   9 -> (1, 0, 0)  16 -> (1,-1, 1)   23 -> (-1, 1,-1)
    /// 3 -> (0, 1, 0)  10 -> (1, 0, 1)  17 -> (1,-1,-1)   24 -> (-1,-1, 0)
    /// 4 -> (0, 1, 1)  11 -> (1, 0,-1)  18 -> (-1, 0, 0)  25 -> (-1,-1, 1)
    /// 5 -> (0, 1,-1)  12 -> (1, 1, 0)  19 -> (-1, 0, 1)  26 -> (-1,-1,-1)
    /// 6 -> (0,-1, 0)  13 -> (1, 1, 1)  20 -> (-1, 0,-1)
    /// </pre>
    ///
    /// Whereas the inner index does not have this roll
    /// <pre class="rust">
    /// 0 -> (-1,-1,-1)   7 -> (-1, 1, 0)  14 -> (0, 0, 1)  21 -> (1, 0,-1)
    /// 1 -> (-1,-1, 0)   8 -> (-1, 1, 1)  15 -> (0, 1,-1)  22 -> (1, 0, 0)
    /// 2 -> (-1,-1, 1)   9 -> (0,-1,-1)   16 -> (0, 1, 0)  23 -> (1, 0, 1)
    /// 3 -> (-1, 0,-1)  10 -> (0,-1, 0)   17 -> (0, 1, 1)  24 -> (1, 1,-1)
    /// 4 -> (-1, 0, 0)  11 -> (0,-1, 1)   18 -> (1,-1,-1)  25 -> (1, 1, 0)
    /// 5 -> (-1, 0, 1)  12 -> (0, 0,-1)   19 -> (1,-1, 0)  26 -> (1, 1, 1)
    /// 6 -> (-1, 1,-1)  13 -> (0, 0, 0)   20 -> (1,-1, 1)
    /// </pre>
    /// BEWARE: The code for this is god awful
    fn new(size: &Size) -> Self {
        let mut di = [[0isize; 27]; 27];
        let index = Shift::index_gen(size);

        let _x: isize = -(size.y * size.z);
        let _xx: isize = (size.x * size.y * size.z) + _x;
        let x: isize = size.y * size.z;
        let xx: isize = -(size.x * size.y * size.z) + x;
        let _y: isize = -size.z;
        let _yy: isize = (size.y * size.z) + _y;
        let y: isize = size.z;
        let yy: isize = -(size.y * size.z) + y;
        let _z: isize = -1;
        let _zz: isize = size.z + _z;
        let z: isize = 1;
        let zz: isize = -size.z + z;

        // GET READY FOR SOME FRESH HELL VECTOR CREATION
        // We need to make the 27 index shifts for the 27 different positions
        // we can be in the density grid

        // (0, 0, 0)
        di[0] = [
            _x + _y + _z,
            _x + _y,
            _x + _y + z,
            _x + _z,
            _x,
            _x + z,
            _x + y + _z,
            _x + y,
            _x + y + z,
            _y + _z,
            _y,
            _y + z,
            _z,
            0,
            z,
            y + _z,
            y,
            y + z,
            x + _y + _z,
            x + _y,
            x + _y + z,
            x + _z,
            x,
            x + z,
            x + y + _z,
            x + y,
            x + y + z,
        ];
        // (0, 0, 1)
        di[1] = [
            _x + _y + _z,
            _x + _y,
            _x + _y + zz,
            _x + _z,
            _x,
            _x + zz,
            _x + y + _z,
            _x + y,
            _x + y + zz,
            _y + _z,
            _y,
            _y + zz,
            _z,
            0,
            zz,
            y + _z,
            y,
            y + zz,
            x + _y + _z,
            x + _y,
            x + _y + zz,
            x + _z,
            x,
            x + zz,
            x + y + _z,
            x + y,
            x + y + zz,
        ];
        // (0, 0,-1)
        di[2] = [
            _x + _y + _zz,
            _x + _y,
            _x + _y + z,
            _x + _zz,
            _x,
            _x + z,
            _x + y + _zz,
            _x + y,
            _x + y + z,
            _y + _zz,
            _y,
            _y + z,
            _zz,
            0,
            z,
            y + _zz,
            y,
            y + z,
            x + _y + _zz,
            x + _y,
            x + _y + z,
            x + _zz,
            x,
            x + z,
            x + y + _zz,
            x + y,
            x + y + z,
        ];
        // (0, 1, 0)
        di[3] = [
            _x + _y + _z,
            _x + _y,
            _x + _y + z,
            _x + _z,
            _x,
            _x + z,
            _x + yy + _z,
            _x + yy,
            _x + yy + z,
            _y + _z,
            _y,
            _y + z,
            _z,
            0,
            z,
            yy + _z,
            yy,
            yy + z,
            x + _y + _z,
            x + _y,
            x + _y + z,
            x + _z,
            x,
            x + z,
            x + yy + _z,
            x + yy,
            x + yy + z,
        ];
        // (0, 1, 1)
        di[4] = [
            _x + _y + _z,
            _x + _y,
            _x + _y + zz,
            _x + _z,
            _x,
            _x + zz,
            _x + yy + _z,
            _x + yy,
            _x + yy + zz,
            _y + _z,
            _y,
            _y + zz,
            _z,
            0,
            zz,
            yy + _z,
            yy,
            yy + zz,
            x + _y + _z,
            x + _y,
            x + _y + zz,
            x + _z,
            x,
            x + zz,
            x + yy + _z,
            x + yy,
            x + yy + zz,
        ];
        // (0, 1,-1)
        di[5] = [
            _x + _y + _zz,
            _x + _y,
            _x + _y + z,
            _x + _zz,
            _x,
            _x + z,
            _x + yy + _zz,
            _x + yy,
            _x + yy + z,
            _y + _zz,
            _y,
            _y + z,
            _zz,
            0,
            z,
            yy + _zz,
            yy,
            yy + z,
            x + _y + _zz,
            x + _y,
            x + _y + z,
            x + _zz,
            x,
            x + z,
            x + yy + _zz,
            x + yy,
            x + yy + z,
        ];
        // (0,-1, 0)
        di[6] = [
            _x + _yy + _z,
            _x + _yy,
            _x + _yy + z,
            _x + _z,
            _x,
            _x + z,
            _x + y + _z,
            _x + y,
            _x + y + z,
            _yy + _z,
            _yy,
            _yy + z,
            _z,
            0,
            z,
            y + _z,
            y,
            y + z,
            x + _yy + _z,
            x + _yy,
            x + _yy + z,
            x + _z,
            x,
            x + z,
            x + y + _z,
            x + y,
            x + y + z,
        ];
        // (0,-1, 1)
        di[7] = [
            _x + _yy + _z,
            _x + _yy,
            _x + _yy + zz,
            _x + _z,
            _x,
            _x + zz,
            _x + y + _z,
            _x + y,
            _x + y + zz,
            _yy + _z,
            _yy,
            _yy + zz,
            _z,
            0,
            zz,
            y + _z,
            y,
            y + zz,
            x + _yy + _z,
            x + _yy,
            x + _yy + zz,
            x + _z,
            x,
            x + zz,
            x + y + _z,
            x + y,
            x + y + zz,
        ];
        // (0,-1,-1)
        di[8] = [
            _x + _yy + _zz,
            _x + _yy,
            _x + _yy + z,
            _x + _zz,
            _x,
            _x + z,
            _x + y + _zz,
            _x + y,
            _x + y + z,
            _yy + _zz,
            _yy,
            _yy + z,
            _zz,
            0,
            z,
            y + _zz,
            y,
            y + z,
            x + _yy + _zz,
            x + _yy,
            x + _yy + z,
            x + _zz,
            x,
            x + z,
            x + y + _zz,
            x + y,
            x + y + z,
        ];
        // (1, 0, 0)
        di[9] = [
            _x + _y + _z,
            _x + _y,
            _x + _y + z,
            _x + _z,
            _x,
            _x + z,
            _x + y + _z,
            _x + y,
            _x + y + z,
            _y + _z,
            _y,
            _y + z,
            _z,
            0,
            z,
            y + _z,
            y,
            y + z,
            xx + _y + _z,
            xx + _y,
            xx + _y + z,
            xx + _z,
            xx,
            xx + z,
            xx + y + _z,
            xx + y,
            xx + y + z,
        ];
        // (1, 0, 1)
        di[10] = [
            _x + _y + _z,
            _x + _y,
            _x + _y + zz,
            _x + _z,
            _x,
            _x + zz,
            _x + y + _z,
            _x + y,
            _x + y + zz,
            _y + _z,
            _y,
            _y + zz,
            _z,
            0,
            zz,
            y + _z,
            y,
            y + zz,
            xx + _y + _z,
            xx + _y,
            xx + _y + zz,
            xx + _z,
            xx,
            xx + zz,
            xx + y + _z,
            xx + y,
            xx + y + zz,
        ];
        // (1, 0,-1)
        di[11] = [
            _x + _y + _zz,
            _x + _y,
            _x + _y + z,
            _x + _zz,
            _x,
            _x + z,
            _x + y + _zz,
            _x + y,
            _x + y + z,
            _y + _zz,
            _y,
            _y + z,
            _zz,
            0,
            z,
            y + _zz,
            y,
            y + z,
            xx + _y + _zz,
            xx + _y,
            xx + _y + z,
            xx + _zz,
            xx,
            xx + z,
            xx + y + _zz,
            xx + y,
            xx + y + z,
        ];
        // (1, 1, 0)
        di[12] = [
            _x + _y + _z,
            _x + _y,
            _x + _y + z,
            _x + _z,
            _x,
            _x + z,
            _x + yy + _z,
            _x + yy,
            _x + yy + z,
            _y + _z,
            _y,
            _y + z,
            _z,
            0,
            z,
            yy + _z,
            yy,
            yy + z,
            xx + _y + _z,
            xx + _y,
            xx + _y + z,
            xx + _z,
            xx,
            xx + z,
            xx + yy + _z,
            xx + yy,
            xx + yy + z,
        ];
        // (1, 1, 1)
        di[13] = [
            _x + _y + _z,
            _x + _y,
            _x + _y + zz,
            _x + _z,
            _x,
            _x + zz,
            _x + yy + _z,
            _x + yy,
            _x + yy + zz,
            _y + _z,
            _y,
            _y + zz,
            _z,
            0,
            zz,
            yy + _z,
            yy,
            yy + zz,
            xx + _y + _z,
            xx + _y,
            xx + _y + zz,
            xx + _z,
            xx,
            xx + zz,
            xx + yy + _z,
            xx + yy,
            xx + yy + zz,
        ];
        // (1, 1,-1)
        di[14] = [
            _x + _y + _zz,
            _x + _y,
            _x + _y + z,
            _x + _zz,
            _x,
            _x + z,
            _x + yy + _zz,
            _x + yy,
            _x + yy + z,
            _y + _zz,
            _y,
            _y + z,
            _zz,
            0,
            z,
            yy + _zz,
            yy,
            yy + z,
            xx + _y + _zz,
            xx + _y,
            xx + _y + z,
            xx + _zz,
            xx,
            xx + z,
            xx + yy + _zz,
            xx + yy,
            xx + yy + z,
        ];
        // (1,-1, 0)
        di[15] = [
            _x + _yy + _z,
            _x + _yy,
            _x + _yy + z,
            _x + _z,
            _x,
            _x + z,
            _x + y + _z,
            _x + y,
            _x + y + z,
            _yy + _z,
            _yy,
            _yy + z,
            _z,
            0,
            z,
            y + _z,
            y,
            y + z,
            xx + _yy + _z,
            xx + _yy,
            xx + _yy + z,
            xx + _z,
            xx,
            xx + z,
            xx + y + _z,
            xx + y,
            xx + y + z,
        ];
        // (1,-1, 1)
        di[16] = [
            _x + _yy + _z,
            _x + _yy,
            _x + _yy + zz,
            _x + _z,
            _x,
            _x + zz,
            _x + y + _z,
            _x + y,
            _x + y + zz,
            _yy + _z,
            _yy,
            _yy + zz,
            _z,
            0,
            zz,
            y + _z,
            y,
            y + zz,
            xx + _yy + _z,
            xx + _yy,
            xx + _yy + zz,
            xx + _z,
            xx,
            xx + zz,
            xx + y + _z,
            xx + y,
            xx + y + zz,
        ];
        // (1,-1,-1)
        di[17] = [
            _x + _yy + _zz,
            _x + _yy,
            _x + _yy + z,
            _x + _zz,
            _x,
            _x + z,
            _x + y + _zz,
            _x + y,
            _x + y + z,
            _yy + _zz,
            _yy,
            _yy + z,
            _zz,
            0,
            z,
            y + _zz,
            y,
            y + z,
            xx + _yy + _zz,
            xx + _yy,
            xx + _yy + z,
            xx + _zz,
            xx,
            xx + z,
            xx + y + _zz,
            xx + y,
            xx + y + z,
        ];
        // (-1, 0, 0)
        di[18] = [
            _xx + _y + _z,
            _xx + _y,
            _xx + _y + z,
            _xx + _z,
            _xx,
            _xx + z,
            _xx + y + _z,
            _xx + y,
            _xx + y + z,
            _y + _z,
            _y,
            _y + z,
            _z,
            0,
            z,
            y + _z,
            y,
            y + z,
            x + _y + _z,
            x + _y,
            x + _y + z,
            x + _z,
            x,
            x + z,
            x + y + _z,
            x + y,
            x + y + z,
        ];
        // (-1, 0, 1)
        di[19] = [
            _xx + _y + _z,
            _xx + _y,
            _xx + _y + zz,
            _xx + _z,
            _xx,
            _xx + zz,
            _xx + y + _z,
            _xx + y,
            _xx + y + zz,
            _y + _z,
            _y,
            _y + zz,
            _z,
            0,
            zz,
            y + _z,
            y,
            y + zz,
            x + _y + _z,
            x + _y,
            x + _y + zz,
            x + _z,
            x,
            x + zz,
            x + y + _z,
            x + y,
            x + y + zz,
        ];
        // (-1, 0,-1)
        di[20] = [
            _xx + _y + _zz,
            _xx + _y,
            _xx + _y + z,
            _xx + _zz,
            _xx,
            _xx + z,
            _xx + y + _zz,
            _xx + y,
            _xx + y + z,
            _y + _zz,
            _y,
            _y + z,
            _zz,
            0,
            z,
            y + _zz,
            y,
            y + z,
            x + _y + _zz,
            x + _y,
            x + _y + z,
            x + _zz,
            x,
            x + z,
            x + y + _zz,
            x + y,
            x + y + z,
        ];
        // (-1, 1, 0)
        di[21] = [
            _xx + _y + _z,
            _xx + _y,
            _xx + _y + z,
            _xx + _z,
            _xx,
            _xx + z,
            _xx + yy + _z,
            _xx + yy,
            _xx + yy + z,
            _y + _z,
            _y,
            _y + z,
            _z,
            0,
            z,
            yy + _z,
            yy,
            yy + z,
            x + _y + _z,
            x + _y,
            x + _y + z,
            x + _z,
            x,
            x + z,
            x + yy + _z,
            x + yy,
            x + yy + z,
        ];
        // (-1, 1, 1)
        di[22] = [
            _xx + _y + _z,
            _xx + _y,
            _xx + _y + zz,
            _xx + _z,
            _xx,
            _xx + zz,
            _xx + yy + _z,
            _xx + yy,
            _xx + yy + zz,
            _y + _z,
            _y,
            _y + zz,
            _z,
            0,
            zz,
            yy + _z,
            yy,
            yy + zz,
            x + _y + _z,
            x + _y,
            x + _y + zz,
            x + _z,
            x,
            x + zz,
            x + yy + _z,
            x + yy,
            x + yy + zz,
        ];
        // (-1, 1,-1)
        di[23] = [
            _xx + _y + _zz,
            _xx + _y,
            _xx + _y + z,
            _xx + _zz,
            _xx,
            _xx + z,
            _xx + yy + _zz,
            _xx + yy,
            _xx + yy + z,
            _y + _zz,
            _y,
            _y + z,
            _zz,
            0,
            z,
            yy + _zz,
            yy,
            yy + z,
            x + _y + _zz,
            x + _y,
            x + _y + z,
            x + _zz,
            x,
            x + z,
            x + yy + _zz,
            x + yy,
            x + yy + z,
        ];
        // (-1,-1, 0)
        di[24] = [
            _xx + _yy + _z,
            _xx + _yy,
            _xx + _yy + z,
            _xx + _z,
            _xx,
            _xx + z,
            _xx + y + _z,
            _xx + y,
            _xx + y + z,
            _yy + _z,
            _yy,
            _yy + z,
            _z,
            0,
            z,
            y + _z,
            y,
            y + z,
            x + _yy + _z,
            x + _yy,
            x + _yy + z,
            x + _z,
            x,
            x + z,
            x + y + _z,
            x + y,
            x + y + z,
        ];
        // (-1,-1, 1)
        di[25] = [
            _xx + _yy + _z,
            _xx + _yy,
            _xx + _yy + zz,
            _xx + _z,
            _xx,
            _xx + zz,
            _xx + y + _z,
            _xx + y,
            _xx + y + zz,
            _yy + _z,
            _yy,
            _yy + zz,
            _z,
            0,
            zz,
            y + _z,
            y,
            y + zz,
            x + _yy + _z,
            x + _yy,
            x + _yy + zz,
            x + _z,
            x,
            x + zz,
            x + y + _z,
            x + y,
            x + y + zz,
        ];
        // (-1,-1,-1)
        di[26] = [
            _xx + _yy + _zz,
            _xx + _yy,
            _xx + _yy + z,
            _xx + _zz,
            _xx,
            _xx + z,
            _xx + y + _zz,
            _xx + y,
            _xx + y + z,
            _yy + _zz,
            _yy,
            _yy + z,
            _zz,
            0,
            z,
            y + _zz,
            y,
            y + z,
            x + _yy + _zz,
            x + _yy,
            x + _yy + z,
            x + _zz,
            x,
            x + z,
            x + y + _zz,
            x + y,
            x + y + z,
        ];
        let pbc = [
            // (0, 0, 0)
            [
                [0, 0, 0], //(-1, -1, -1)
                [0, 0, 0], //(-1, -1,  0)
                [0, 0, 0], //(-1, -1,  1)
                [0, 0, 0], //(-1,  0, -1)
                [0, 0, 0], //(-1,  0,  0)
                [0, 0, 0], //(-1,  0,  1)
                [0, 0, 0], //(-1,  1, -1)
                [0, 0, 0], //(-1,  1,  0)
                [0, 0, 0], //(-1,  1,  1)
                [0, 0, 0], //( 0, -1, -1)
                [0, 0, 0], //( 0, -1,  0)
                [0, 0, 0], //( 0, -1,  1)
                [0, 0, 0], //( 0,  0, -1)
                [0, 0, 0], //( 0,  0,  0)
                [0, 0, 0], //( 0,  0,  1)
                [0, 0, 0], //( 0,  1, -1)
                [0, 0, 0], //( 0,  1,  0)
                [0, 0, 0], //( 0,  1,  1)
                [0, 0, 0], //( 1, -1, -1)
                [0, 0, 0], //( 1, -1,  0)
                [0, 0, 0], //( 1, -1,  1)
                [0, 0, 0], //( 1,  0, -1)
                [0, 0, 0], //( 1,  0,  0)
                [0, 0, 0], //( 1,  0,  1)
                [0, 0, 0], //( 1,  1, -1)
                [0, 0, 0], //( 1,  1,  0)
                [0, 0, 0], //( 1,  1,  1)
            ],
            // (0, 0, 1)
            [
                [0, 0, 0], //(-1, -1, -1)
                [0, 0, 0], //(-1, -1,  0)
                [0, 0, 1], //(-1, -1,  1)
                [0, 0, 0], //(-1,  0, -1)
                [0, 0, 0], //(-1,  0,  0)
                [0, 0, 1], //(-1,  0,  1)
                [0, 0, 0], //(-1,  1, -1)
                [0, 0, 0], //(-1,  1,  0)
                [0, 0, 1], //(-1,  1,  1)
                [0, 0, 0], //( 0, -1, -1)
                [0, 0, 0], //( 0, -1,  0)
                [0, 0, 1], //( 0, -1,  1)
                [0, 0, 0], //( 0,  0, -1)
                [0, 0, 0], //( 0,  0,  0)
                [0, 0, 1], //( 0,  0,  1)
                [0, 0, 0], //( 0,  1, -1)
                [0, 0, 0], //( 0,  1,  0)
                [0, 0, 1], //( 0,  1,  1)
                [0, 0, 0], //( 1, -1, -1)
                [0, 0, 0], //( 1, -1,  0)
                [0, 0, 1], //( 1, -1,  1)
                [0, 0, 0], //( 1,  0, -1)
                [0, 0, 0], //( 1,  0,  0)
                [0, 0, 1], //( 1,  0,  1)
                [0, 0, 0], //( 1,  1, -1)
                [0, 0, 0], //( 1,  1,  0)
                [0, 0, 1], //( 1,  1,  1)
            ],
            // (0, 0,-1)
            [
                [0, 0, -1], //(-1, -1, -1)
                [0, 0, 0],  //(-1, -1,  0)
                [0, 0, 0],  //(-1, -1,  1)
                [0, 0, -1], //(-1,  0, -1)
                [0, 0, 0],  //(-1,  0,  0)
                [0, 0, 0],  //(-1,  0,  1)
                [0, 0, -1], //(-1,  1, -1)
                [0, 0, 0],  //(-1,  1,  0)
                [0, 0, 0],  //(-1,  1,  1)
                [0, 0, -1], //( 0, -1, -1)
                [0, 0, 0],  //( 0, -1,  0)
                [0, 0, 0],  //( 0, -1,  1)
                [0, 0, -1], //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 0],  //( 0,  0,  1)
                [0, 0, -1], //( 0,  1, -1)
                [0, 0, 0],  //( 0,  1,  0)
                [0, 0, 0],  //( 0,  1,  1)
                [0, 0, -1], //( 1, -1, -1)
                [0, 0, 0],  //( 1, -1,  0)
                [0, 0, 0],  //( 1, -1,  1)
                [0, 0, -1], //( 1,  0, -1)
                [0, 0, 0],  //( 1,  0,  0)
                [0, 0, 0],  //( 1,  0,  1)
                [0, 0, -1], //( 1,  1, -1)
                [0, 0, 0],  //( 1,  1,  0)
                [0, 0, 0],  //( 1,  1,  1)
            ],
            // (0, 1, 0)
            [
                [0, 0, 0], //(-1, -1, -1)
                [0, 0, 0], //(-1, -1,  0)
                [0, 0, 0], //(-1, -1,  1)
                [0, 0, 0], //(-1,  0, -1)
                [0, 0, 0], //(-1,  0,  0)
                [0, 0, 0], //(-1,  0,  1)
                [0, 1, 0], //(-1,  1, -1)
                [0, 1, 0], //(-1,  1,  0)
                [0, 1, 0], //(-1,  1,  1)
                [0, 0, 0], //( 0, -1, -1)
                [0, 0, 0], //( 0, -1,  0)
                [0, 0, 0], //( 0, -1,  1)
                [0, 0, 0], //( 0,  0, -1)
                [0, 0, 0], //( 0,  0,  0)
                [0, 0, 0], //( 0,  0,  1)
                [0, 1, 0], //( 0,  1, -1)
                [0, 1, 0], //( 0,  1,  0)
                [0, 1, 0], //( 0,  1,  1)
                [0, 0, 0], //( 1, -1, -1)
                [0, 0, 0], //( 1, -1,  0)
                [0, 0, 0], //( 1, -1,  1)
                [0, 0, 0], //( 1,  0, -1)
                [0, 0, 0], //( 1,  0,  0)
                [0, 0, 0], //( 1,  0,  1)
                [0, 1, 0], //( 1,  1, -1)
                [0, 1, 0], //( 1,  1,  0)
                [0, 1, 0], //( 1,  1,  1)
            ],
            // (0, 1, 1)
            [
                [0, 0, 0], //(-1, -1, -1)
                [0, 0, 0], //(-1, -1,  0)
                [0, 0, 1], //(-1, -1,  1)
                [0, 0, 0], //(-1,  0, -1)
                [0, 0, 0], //(-1,  0,  0)
                [0, 0, 1], //(-1,  0,  1)
                [0, 1, 0], //(-1,  1, -1)
                [0, 1, 0], //(-1,  1,  0)
                [0, 1, 1], //(-1,  1,  1)
                [0, 0, 0], //( 0, -1, -1)
                [0, 0, 0], //( 0, -1,  0)
                [0, 0, 1], //( 0, -1,  1)
                [0, 0, 0], //( 0,  0, -1)
                [0, 0, 0], //( 0,  0,  0)
                [0, 0, 1], //( 0,  0,  1)
                [0, 1, 0], //( 0,  1, -1)
                [0, 1, 0], //( 0,  1,  0)
                [0, 1, 1], //( 0,  1,  1)
                [0, 0, 0], //( 1, -1, -1)
                [0, 0, 0], //( 1, -1,  0)
                [0, 0, 1], //( 1, -1,  1)
                [0, 0, 0], //( 1,  0, -1)
                [0, 0, 0], //( 1,  0,  0)
                [0, 0, 1], //( 1,  0,  1)
                [0, 1, 0], //( 1,  1, -1)
                [0, 1, 0], //( 1,  1,  0)
                [0, 1, 1], //( 1,  1,  1)
            ],
            // (0, 1,-1)
            [
                [0, 0, -1], //(-1, -1, -1)
                [0, 0, 0],  //(-1, -1,  0)
                [0, 0, 0],  //(-1, -1,  1)
                [0, 0, -1], //(-1,  0, -1)
                [0, 0, 0],  //(-1,  0,  0)
                [0, 0, 0],  //(-1,  0,  1)
                [0, 1, -1], //(-1,  1, -1)
                [0, 1, 0],  //(-1,  1,  0)
                [0, 1, 0],  //(-1,  1,  1)
                [0, 0, -1], //( 0, -1, -1)
                [0, 0, 0],  //( 0, -1,  0)
                [0, 0, 0],  //( 0, -1,  1)
                [0, 0, -1], //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 0],  //( 0,  0,  1)
                [0, 1, -1], //( 0,  1, -1)
                [0, 1, 0],  //( 0,  1,  0)
                [0, 1, 0],  //( 0,  1,  1)
                [0, 0, -1], //( 1, -1, -1)
                [0, 0, 0],  //( 1, -1,  0)
                [0, 0, 0],  //( 1, -1,  1)
                [0, 0, -1], //( 1,  0, -1)
                [0, 0, 0],  //( 1,  0,  0)
                [0, 0, 0],  //( 1,  0,  1)
                [0, 1, -1], //( 1,  1, -1)
                [0, 1, 0],  //( 1,  1,  0)
                [0, 1, 0],  //( 1,  1,  1)
            ],
            // (0,-1, 0)
            [
                [0, -1, 0], //(-1, -1, -1)
                [0, -1, 0], //(-1, -1,  0)
                [0, -1, 0], //(-1, -1,  1)
                [0, 0, 0],  //(-1,  0, -1)
                [0, 0, 0],  //(-1,  0,  0)
                [0, 0, 0],  //(-1,  0,  1)
                [0, 0, 0],  //(-1,  1, -1)
                [0, 0, 0],  //(-1,  1,  0)
                [0, 0, 0],  //(-1,  1,  1)
                [0, -1, 0], //( 0, -1, -1)
                [0, -1, 0], //( 0, -1,  0)
                [0, -1, 0], //( 0, -1,  1)
                [0, 0, 0],  //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 0],  //( 0,  0,  1)
                [0, 0, 0],  //( 0,  1, -1)
                [0, 0, 0],  //( 0,  1,  0)
                [0, 0, 0],  //( 0,  1,  1)
                [0, -1, 0], //( 1, -1, -1)
                [0, -1, 0], //( 1, -1,  0)
                [0, -1, 0], //( 1, -1,  1)
                [0, 0, 0],  //( 1,  0, -1)
                [0, 0, 0],  //( 1,  0,  0)
                [0, 0, 0],  //( 1,  0,  1)
                [0, 0, 0],  //( 1,  1, -1)
                [0, 0, 0],  //( 1,  1,  0)
                [0, 0, 0],  //( 1,  1,  1)
            ],
            // (0,-1, 1)
            [
                [0, -1, 0], //(-1, -1, -1)
                [0, -1, 0], //(-1, -1,  0)
                [0, -1, 1], //(-1, -1,  1)
                [0, 0, 0],  //(-1,  0, -1)
                [0, 0, 0],  //(-1,  0,  0)
                [0, 0, 1],  //(-1,  0,  1)
                [0, 0, 0],  //(-1,  1, -1)
                [0, 0, 0],  //(-1,  1,  0)
                [0, 0, 1],  //(-1,  1,  1)
                [0, -1, 0], //( 0, -1, -1)
                [0, -1, 0], //( 0, -1,  0)
                [0, -1, 1], //( 0, -1,  1)
                [0, 0, 0],  //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 1],  //( 0,  0,  1)
                [0, 0, 0],  //( 0,  1, -1)
                [0, 0, 0],  //( 0,  1,  0)
                [0, 0, 1],  //( 0,  1,  1)
                [0, -1, 0], //( 1, -1, -1)
                [0, -1, 0], //( 1, -1,  0)
                [0, -1, 1], //( 1, -1,  1)
                [0, 0, 0],  //( 1,  0, -1)
                [0, 0, 0],  //( 1,  0,  0)
                [0, 0, 1],  //( 1,  0,  1)
                [0, 0, 0],  //( 1,  1, -1)
                [0, 0, 0],  //( 1,  1,  0)
                [0, 0, 1],  //( 1,  1,  1)
            ],
            // (0,-1,-1)
            [
                [0, -1, -1], //(-1, -1, -1)
                [0, -1, 0],  //(-1, -1,  0)
                [0, -1, 0],  //(-1, -1,  1)
                [0, 0, -1],  //(-1,  0, -1)
                [0, 0, 0],   //(-1,  0,  0)
                [0, 0, 0],   //(-1,  0,  1)
                [0, 0, -1],  //(-1,  1, -1)
                [0, 0, 0],   //(-1,  1,  0)
                [0, 0, 0],   //(-1,  1,  1)
                [0, -1, -1], //( 0, -1, -1)
                [0, -1, 0],  //( 0, -1,  0)
                [0, -1, 0],  //( 0, -1,  1)
                [0, 0, -1],  //( 0,  0, -1)
                [0, 0, 0],   //( 0,  0,  0)
                [0, 0, 0],   //( 0,  0,  1)
                [0, 0, -1],  //( 0,  1, -1)
                [0, 0, 0],   //( 0,  1,  0)
                [0, 0, 0],   //( 0,  1,  1)
                [0, -1, -1], //( 1, -1, -1)
                [0, -1, 0],  //( 1, -1,  0)
                [0, -1, 0],  //( 1, -1,  1)
                [0, 0, -1],  //( 1,  0, -1)
                [0, 0, 0],   //( 1,  0,  0)
                [0, 0, 0],   //( 1,  0,  1)
                [0, 0, -1],  //( 1,  1, -1)
                [0, 0, 0],   //( 1,  1,  0)
                [0, 0, 0],   //( 1,  1,  1)
            ],
            // (1, 0, 0)
            [
                [0, 0, 0], //(-1, -1, -1)
                [0, 0, 0], //(-1, -1,  0)
                [0, 0, 0], //(-1, -1,  1)
                [0, 0, 0], //(-1,  0, -1)
                [0, 0, 0], //(-1,  0,  0)
                [0, 0, 0], //(-1,  0,  1)
                [0, 0, 0], //(-1,  1, -1)
                [0, 0, 0], //(-1,  1,  0)
                [0, 0, 0], //(-1,  1,  1)
                [0, 0, 0], //( 0, -1, -1)
                [0, 0, 0], //( 0, -1,  0)
                [0, 0, 0], //( 0, -1,  1)
                [0, 0, 0], //( 0,  0, -1)
                [0, 0, 0], //( 0,  0,  0)
                [0, 0, 0], //( 0,  0,  1)
                [0, 0, 0], //( 0,  1, -1)
                [0, 0, 0], //( 0,  1,  0)
                [0, 0, 0], //( 0,  1,  1)
                [1, 0, 0], //( 1, -1, -1)
                [1, 0, 0], //( 1, -1,  0)
                [1, 0, 0], //( 1, -1,  1)
                [1, 0, 0], //( 1,  0, -1)
                [1, 0, 0], //( 1,  0,  0)
                [1, 0, 0], //( 1,  0,  1)
                [1, 0, 0], //( 1,  1, -1)
                [1, 0, 0], //( 1,  1,  0)
                [1, 0, 0], //( 1,  1,  1)
            ],
            // (1, 0, 1)
            [
                [0, 0, 0], //(-1, -1, -1)
                [0, 0, 0], //(-1, -1,  0)
                [0, 0, 1], //(-1, -1,  1)
                [0, 0, 0], //(-1,  0, -1)
                [0, 0, 0], //(-1,  0,  0)
                [0, 0, 1], //(-1,  0,  1)
                [0, 0, 0], //(-1,  1, -1)
                [0, 0, 0], //(-1,  1,  0)
                [0, 0, 1], //(-1,  1,  1)
                [0, 0, 0], //( 0, -1, -1)
                [0, 0, 0], //( 0, -1,  0)
                [0, 0, 1], //( 0, -1,  1)
                [0, 0, 0], //( 0,  0, -1)
                [0, 0, 0], //( 0,  0,  0)
                [0, 0, 1], //( 0,  0,  1)
                [0, 0, 0], //( 0,  1, -1)
                [0, 0, 0], //( 0,  1,  0)
                [0, 0, 1], //( 0,  1,  1)
                [1, 0, 0], //( 1, -1, -1)
                [1, 0, 0], //( 1, -1,  0)
                [1, 0, 1], //( 1, -1,  1)
                [1, 0, 0], //( 1,  0, -1)
                [1, 0, 0], //( 1,  0,  0)
                [1, 0, 1], //( 1,  0,  1)
                [1, 0, 0], //( 1,  1, -1)
                [1, 0, 0], //( 1,  1,  0)
                [1, 0, 1], //( 1,  1,  1)
            ],
            // (1, 0,-1)
            [
                [0, 0, -1], //(-1, -1, -1)
                [0, 0, 0],  //(-1, -1,  0)
                [0, 0, 0],  //(-1, -1,  1)
                [0, 0, -1], //(-1,  0, -1)
                [0, 0, 0],  //(-1,  0,  0)
                [0, 0, 0],  //(-1,  0,  1)
                [0, 0, -1], //(-1,  1, -1)
                [0, 0, 0],  //(-1,  1,  0)
                [0, 0, 0],  //(-1,  1,  1)
                [0, 0, -1], //( 0, -1, -1)
                [0, 0, 0],  //( 0, -1,  0)
                [0, 0, 0],  //( 0, -1,  1)
                [0, 0, -1], //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 0],  //( 0,  0,  1)
                [0, 0, -1], //( 0,  1, -1)
                [0, 0, 0],  //( 0,  1,  0)
                [0, 0, 0],  //( 0,  1,  1)
                [1, 0, -1], //( 1, -1, -1)
                [1, 0, 0],  //( 1, -1,  0)
                [1, 0, 0],  //( 1, -1,  1)
                [1, 0, -1], //( 1,  0, -1)
                [1, 0, 0],  //( 1,  0,  0)
                [1, 0, 0],  //( 1,  0,  1)
                [1, 0, -1], //( 1,  1, -1)
                [1, 0, 0],  //( 1,  1,  0)
                [1, 0, 0],  //( 1,  1,  1)
            ],
            // (1, 1, 0)
            [
                [0, 0, 0], //(-1, -1, -1)
                [0, 0, 0], //(-1, -1,  0)
                [0, 0, 0], //(-1, -1,  1)
                [0, 0, 0], //(-1,  0, -1)
                [0, 0, 0], //(-1,  0,  0)
                [0, 0, 0], //(-1,  0,  1)
                [0, 1, 0], //(-1,  1, -1)
                [0, 1, 0], //(-1,  1,  0)
                [0, 1, 0], //(-1,  1,  1)
                [0, 0, 0], //( 0, -1, -1)
                [0, 0, 0], //( 0, -1,  0)
                [0, 0, 0], //( 0, -1,  1)
                [0, 0, 0], //( 0,  0, -1)
                [0, 0, 0], //( 0,  0,  0)
                [0, 0, 0], //( 0,  0,  1)
                [0, 1, 0], //( 0,  1, -1)
                [0, 1, 0], //( 0,  1,  0)
                [0, 1, 0], //( 0,  1,  1)
                [1, 0, 0], //( 1, -1, -1)
                [1, 0, 0], //( 1, -1,  0)
                [1, 0, 0], //( 1, -1,  1)
                [1, 0, 0], //( 1,  0, -1)
                [1, 0, 0], //( 1,  0,  0)
                [1, 0, 0], //( 1,  0,  1)
                [1, 1, 0], //( 1,  1, -1)
                [1, 1, 0], //( 1,  1,  0)
                [1, 1, 0], //( 1,  1,  1)
            ],
            // (1, 1, 1)
            [
                [0, 0, 0], //(-1, -1, -1)
                [0, 0, 0], //(-1, -1,  0)
                [0, 0, 1], //(-1, -1,  1)
                [0, 0, 0], //(-1,  0, -1)
                [0, 0, 0], //(-1,  0,  0)
                [0, 0, 0], //(-1,  0,  1)
                [0, 1, 0], //(-1,  1, -1)
                [0, 1, 0], //(-1,  1,  0)
                [0, 1, 1], //(-1,  1,  1)
                [0, 0, 0], //( 0, -1, -1)
                [0, 0, 0], //( 0, -1,  0)
                [0, 0, 1], //( 0, -1,  1)
                [0, 0, 0], //( 0,  0, -1)
                [0, 0, 0], //( 0,  0,  0)
                [0, 0, 1], //( 0,  0,  1)
                [0, 1, 0], //( 0,  1, -1)
                [0, 1, 0], //( 0,  1,  0)
                [0, 1, 1], //( 0,  1,  1)
                [1, 0, 0], //( 1, -1, -1)
                [1, 0, 0], //( 1, -1,  0)
                [1, 0, 1], //( 1, -1,  1)
                [1, 0, 0], //( 1,  0, -1)
                [1, 0, 0], //( 1,  0,  0)
                [1, 0, 1], //( 1,  0,  1)
                [1, 1, 0], //( 1,  1, -1)
                [1, 1, 0], //( 1,  1,  0)
                [1, 1, 1], //( 1,  1,  1)
            ],
            // (1, 1,-1)
            [
                [0, 0, -1], //(-1, -1, -1)
                [0, 0, 0],  //(-1, -1,  0)
                [0, 0, 0],  //(-1, -1,  1)
                [0, 0, -1], //(-1,  0, -1)
                [0, 0, 0],  //(-1,  0,  0)
                [0, 0, 0],  //(-1,  0,  1)
                [0, 1, -1], //(-1,  1, -1)
                [0, 1, 0],  //(-1,  1,  0)
                [0, 1, 0],  //(-1,  1,  1)
                [0, 0, -1], //( 0, -1, -1)
                [0, 0, 0],  //( 0, -1,  0)
                [0, 0, 0],  //( 0, -1,  1)
                [0, 0, -1], //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 0],  //( 0,  0,  1)
                [0, 1, -1], //( 0,  1, -1)
                [0, 1, 0],  //( 0,  1,  0)
                [0, 1, 0],  //( 0,  1,  1)
                [1, 0, -1], //( 1, -1, -1)
                [1, 0, 0],  //( 1, -1,  0)
                [1, 0, 0],  //( 1, -1,  1)
                [1, 0, -1], //( 1,  0, -1)
                [1, 0, 0],  //( 1,  0,  0)
                [1, 0, 0],  //( 1,  0,  1)
                [1, 1, -1], //( 1,  1, -1)
                [1, 1, 0],  //( 1,  1,  0)
                [1, 1, 0],  //( 1,  1,  1)
            ],
            // (1,-1, 0)
            [
                [0, -1, 0], //(-1, -1, -1)
                [0, -1, 0], //(-1, -1,  0)
                [0, -1, 0], //(-1, -1,  1)
                [0, 0, 0],  //(-1,  0, -1)
                [0, 0, 0],  //(-1,  0,  0)
                [0, 0, 0],  //(-1,  0,  1)
                [0, 0, 0],  //(-1,  1, -1)
                [0, 0, 0],  //(-1,  1,  0)
                [0, 0, 0],  //(-1,  1,  1)
                [0, -1, 0], //( 0, -1, -1)
                [0, -1, 0], //( 0, -1,  0)
                [0, -1, 0], //( 0, -1,  1)
                [0, 0, 0],  //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 0],  //( 0,  0,  1)
                [0, 0, 0],  //( 0,  1, -1)
                [0, 0, 0],  //( 0,  1,  0)
                [0, 0, 0],  //( 0,  1,  1)
                [1, -1, 0], //( 1, -1, -1)
                [1, -1, 0], //( 1, -1,  0)
                [1, -1, 0], //( 1, -1,  1)
                [1, 0, 0],  //( 1,  0, -1)
                [1, 0, 0],  //( 1,  0,  0)
                [1, 0, 0],  //( 1,  0,  1)
                [1, 0, 0],  //( 1,  1, -1)
                [1, 0, 0],  //( 1,  1,  0)
                [1, 0, 0],  //( 1,  1,  1)
            ],
            // (1,-1, 1)
            [
                [0, -1, 0], //(-1, -1, -1)
                [0, -1, 0], //(-1, -1,  0)
                [0, -1, 1], //(-1, -1,  1)
                [0, 0, 0],  //(-1,  0, -1)
                [0, 0, 0],  //(-1,  0,  0)
                [0, 0, 1],  //(-1,  0,  1)
                [0, 0, 0],  //(-1,  1, -1)
                [0, 0, 0],  //(-1,  1,  0)
                [0, 0, 1],  //(-1,  1,  1)
                [0, -1, 0], //( 0, -1, -1)
                [0, -1, 0], //( 0, -1,  0)
                [0, -1, 1], //( 0, -1,  1)
                [0, 0, 0],  //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 1],  //( 0,  0,  1)
                [0, 0, 0],  //( 0,  1, -1)
                [0, 0, 0],  //( 0,  1,  0)
                [0, 0, 1],  //( 0,  1,  1)
                [1, -1, 0], //( 1, -1, -1)
                [1, -1, 0], //( 1, -1,  0)
                [1, -1, 1], //( 1, -1,  1)
                [1, 0, 0],  //( 1,  0, -1)
                [1, 0, 0],  //( 1,  0,  0)
                [1, 0, 1],  //( 1,  0,  1)
                [1, 0, 0],  //( 1,  1, -1)
                [1, 0, 0],  //( 1,  1,  0)
                [1, 0, 1],  //( 1,  1,  1)
            ],
            // (1,-1,-1)
            [
                [0, -1, -1], //(-1, -1, -1)
                [0, -1, 0],  //(-1, -1,  0)
                [0, -1, 0],  //(-1, -1,  1)
                [0, 0, -1],  //(-1,  0, -1)
                [0, 0, 0],   //(-1,  0,  0)
                [0, 0, 0],   //(-1,  0,  1)
                [0, 0, -1],  //(-1,  1, -1)
                [0, 0, 0],   //(-1,  1,  0)
                [0, 0, 0],   //(-1,  1,  1)
                [0, -1, -1], //( 0, -1, -1)
                [0, -1, 0],  //( 0, -1,  0)
                [0, -1, 0],  //( 0, -1,  1)
                [0, 0, -1],  //( 0,  0, -1)
                [0, 0, 0],   //( 0,  0,  0)
                [0, 0, 0],   //( 0,  0,  1)
                [0, 0, -1],  //( 0,  1, -1)
                [0, 0, 0],   //( 0,  1,  0)
                [0, 0, 0],   //( 0,  1,  1)
                [1, -1, -1], //( 1, -1, -1)
                [1, -1, 0],  //( 1, -1,  0)
                [1, -1, 0],  //( 1, -1,  1)
                [1, 0, -1],  //( 1,  0, -1)
                [1, 0, 0],   //( 1,  0,  0)
                [1, 0, 0],   //( 1,  0,  1)
                [1, 0, -1],  //( 1,  1, -1)
                [1, 0, 0],   //( 1,  1,  0)
                [1, 0, 0],   //( 1,  1,  1)
            ],
            // (-1, 0, 0)
            [
                [-1, 0, 0], //(-1, -1, -1)
                [-1, 0, 0], //(-1, -1,  0)
                [-1, 0, 0], //(-1, -1,  1)
                [-1, 0, 0], //(-1,  0, -1)
                [-1, 0, 0], //(-1,  0,  0)
                [-1, 0, 0], //(-1,  0,  1)
                [-1, 0, 0], //(-1,  1, -1)
                [-1, 0, 0], //(-1,  1,  0)
                [-1, 0, 0], //(-1,  1,  1)
                [0, 0, 0],  //( 0, -1, -1)
                [0, 0, 0],  //( 0, -1,  0)
                [0, 0, 0],  //( 0, -1,  1)
                [0, 0, 0],  //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 0],  //( 0,  0,  1)
                [0, 0, 0],  //( 0,  1, -1)
                [0, 0, 0],  //( 0,  1,  0)
                [0, 0, 0],  //( 0,  1,  1)
                [0, 0, 0],  //( 1, -1, -1)
                [0, 0, 0],  //( 1, -1,  0)
                [0, 0, 0],  //( 1, -1,  1)
                [0, 0, 0],  //( 1,  0, -1)
                [0, 0, 0],  //( 1,  0,  0)
                [0, 0, 0],  //( 1,  0,  1)
                [0, 0, 0],  //( 1,  1, -1)
                [0, 0, 0],  //( 1,  1,  0)
                [0, 0, 0],  //( 1,  1,  1)
            ],
            // (-1, 0, 1)
            [
                [-1, 0, 0], //(-1, -1, -1)
                [-1, 0, 0], //(-1, -1,  0)
                [-1, 0, 1], //(-1, -1,  1)
                [-1, 0, 0], //(-1,  0, -1)
                [-1, 0, 0], //(-1,  0,  0)
                [-1, 0, 1], //(-1,  0,  1)
                [-1, 0, 0], //(-1,  1, -1)
                [-1, 0, 0], //(-1,  1,  0)
                [-1, 0, 1], //(-1,  1,  1)
                [0, 0, 0],  //( 0, -1, -1)
                [0, 0, 0],  //( 0, -1,  0)
                [0, 0, 1],  //( 0, -1,  1)
                [0, 0, 0],  //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 1],  //( 0,  0,  1)
                [0, 0, 0],  //( 0,  1, -1)
                [0, 0, 0],  //( 0,  1,  0)
                [0, 0, 1],  //( 0,  1,  1)
                [0, 0, 0],  //( 1, -1, -1)
                [0, 0, 0],  //( 1, -1,  0)
                [0, 0, 1],  //( 1, -1,  1)
                [0, 0, 0],  //( 1,  0, -1)
                [0, 0, 0],  //( 1,  0,  0)
                [0, 0, 1],  //( 1,  0,  1)
                [0, 0, 0],  //( 1,  1, -1)
                [0, 0, 0],  //( 1,  1,  0)
                [0, 0, 1],  //( 1,  1,  1)
            ],
            // (-1, 0,-1)
            [
                [-1, 0, -1], //(-1, -1, -1)
                [-1, 0, 0],  //(-1, -1,  0)
                [-1, 0, 0],  //(-1, -1,  1)
                [-1, 0, -1], //(-1,  0, -1)
                [-1, 0, 0],  //(-1,  0,  0)
                [-1, 0, 0],  //(-1,  0,  1)
                [-1, 0, -1], //(-1,  1, -1)
                [-1, 0, 0],  //(-1,  1,  0)
                [-1, 0, 0],  //(-1,  1,  1)
                [0, 0, -1],  //( 0, -1, -1)
                [0, 0, 0],   //( 0, -1,  0)
                [0, 0, 0],   //( 0, -1,  1)
                [0, 0, -1],  //( 0,  0, -1)
                [0, 0, 0],   //( 0,  0,  0)
                [0, 0, 0],   //( 0,  0,  1)
                [0, 0, -1],  //( 0,  1, -1)
                [0, 0, 0],   //( 0,  1,  0)
                [0, 0, 0],   //( 0,  1,  1)
                [0, 0, -1],  //( 1, -1, -1)
                [0, 0, 0],   //( 1, -1,  0)
                [0, 0, 0],   //( 1, -1,  1)
                [0, 0, -1],  //( 1,  0, -1)
                [0, 0, 0],   //( 1,  0,  0)
                [0, 0, 0],   //( 1,  0,  1)
                [0, 0, -1],  //( 1,  1, -1)
                [0, 0, 0],   //( 1,  1,  0)
                [0, 0, 0],   //( 1,  1,  1)
            ],
            // (-1, 1, 0)
            [
                [-1, 0, 0], //(-1, -1, -1)
                [-1, 0, 0], //(-1, -1,  0)
                [-1, 0, 0], //(-1, -1,  1)
                [-1, 0, 0], //(-1,  0, -1)
                [-1, 0, 0], //(-1,  0,  0)
                [-1, 0, 0], //(-1,  0,  1)
                [-1, 1, 0], //(-1,  1, -1)
                [-1, 1, 0], //(-1,  1,  0)
                [-1, 1, 0], //(-1,  1,  1)
                [0, 0, 0],  //( 0, -1, -1)
                [0, 0, 0],  //( 0, -1,  0)
                [0, 0, 0],  //( 0, -1,  1)
                [0, 0, 0],  //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 0],  //( 0,  0,  1)
                [0, 1, 0],  //( 0,  1, -1)
                [0, 1, 0],  //( 0,  1,  0)
                [0, 1, 0],  //( 0,  1,  1)
                [0, 0, 0],  //( 1, -1, -1)
                [0, 0, 0],  //( 1, -1,  0)
                [0, 0, 0],  //( 1, -1,  1)
                [0, 0, 0],  //( 1,  0, -1)
                [0, 0, 0],  //( 1,  0,  0)
                [0, 0, 0],  //( 1,  0,  1)
                [0, 1, 0],  //( 1,  1, -1)
                [0, 1, 0],  //( 1,  1,  0)
                [0, 1, 0],  //( 1,  1,  1)
            ],
            // (-1, 1, 1)
            [
                [-1, 0, 0], //(-1, -1, -1)
                [-1, 0, 0], //(-1, -1,  0)
                [-1, 0, 1], //(-1, -1,  1)
                [-1, 0, 0], //(-1,  0, -1)
                [-1, 0, 0], //(-1,  0,  0)
                [-1, 0, 1], //(-1,  0,  1)
                [-1, 1, 0], //(-1,  1, -1)
                [-1, 1, 0], //(-1,  1,  0)
                [-1, 1, 1], //(-1,  1,  1)
                [0, 0, 0],  //( 0, -1, -1)
                [0, 0, 0],  //( 0, -1,  0)
                [0, 0, 1],  //( 0, -1,  1)
                [0, 0, 0],  //( 0,  0, -1)
                [0, 0, 0],  //( 0,  0,  0)
                [0, 0, 1],  //( 0,  0,  1)
                [0, 1, 0],  //( 0,  1, -1)
                [0, 1, 0],  //( 0,  1,  0)
                [0, 1, 1],  //( 0,  1,  1)
                [0, 0, 0],  //( 1, -1, -1)
                [0, 0, 0],  //( 1, -1,  0)
                [0, 0, 1],  //( 1, -1,  1)
                [0, 0, 0],  //( 1,  0, -1)
                [0, 0, 0],  //( 1,  0,  0)
                [0, 0, 1],  //( 1,  0,  1)
                [0, 1, 0],  //( 1,  1, -1)
                [0, 1, 0],  //( 1,  1,  0)
                [0, 1, 1],  //( 1,  1,  1)
            ],
            // (-1, 1,-1)
            [
                [-1, 0, -1], //(-1, -1, -1)
                [-1, 0, 0],  //(-1, -1,  0)
                [-1, 0, 0],  //(-1, -1,  1)
                [-1, 0, -1], //(-1,  0, -1)
                [-1, 0, 0],  //(-1,  0,  0)
                [-1, 0, 0],  //(-1,  0,  1)
                [-1, 1, -1], //(-1,  1, -1)
                [-1, 1, 0],  //(-1,  1,  0)
                [-1, 1, 0],  //(-1,  1,  1)
                [0, 0, -1],  //( 0, -1, -1)
                [0, 0, 0],   //( 0, -1,  0)
                [0, 0, 0],   //( 0, -1,  1)
                [0, 0, -1],  //( 0,  0, -1)
                [0, 0, 0],   //( 0,  0,  0)
                [0, 0, 0],   //( 0,  0,  1)
                [0, 1, -1],  //( 0,  1, -1)
                [0, 1, 0],   //( 0,  1,  0)
                [0, 1, 0],   //( 0,  1,  1)
                [0, 0, -1],  //( 1, -1, -1)
                [0, 0, 0],   //( 1, -1,  0)
                [0, 0, 0],   //( 1, -1,  1)
                [0, 0, -1],  //( 1,  0, -1)
                [0, 0, 0],   //( 1,  0,  0)
                [0, 0, 0],   //( 1,  0,  1)
                [0, 1, -1],  //( 1,  1, -1)
                [0, 1, 0],   //( 1,  1,  0)
                [0, 1, 0],   //( 1,  1,  1)
            ],
            // (-1,-1, 0)
            [
                [-1, -1, 0], //(-1, -1, -1)
                [-1, -1, 0], //(-1, -1,  0)
                [-1, -1, 0], //(-1, -1,  1)
                [-1, 0, 0],  //(-1,  0, -1)
                [-1, 0, 0],  //(-1,  0,  0)
                [-1, 0, 0],  //(-1,  0,  1)
                [-1, 0, 0],  //(-1,  1, -1)
                [-1, 0, 0],  //(-1,  1,  0)
                [-1, 0, 0],  //(-1,  1,  1)
                [0, -1, 0],  //( 0, -1, -1)
                [0, -1, 0],  //( 0, -1,  0)
                [0, -1, 0],  //( 0, -1,  1)
                [0, 0, 0],   //( 0,  0, -1)
                [0, 0, 0],   //( 0,  0,  0)
                [0, 0, 0],   //( 0,  0,  1)
                [0, 0, 0],   //( 0,  1, -1)
                [0, 0, 0],   //( 0,  1,  0)
                [0, 0, 0],   //( 0,  1,  1)
                [0, -1, 0],  //( 1, -1, -1)
                [0, -1, 0],  //( 1, -1,  0)
                [0, -1, 0],  //( 1, -1,  1)
                [0, 0, 0],   //( 1,  0, -1)
                [0, 0, 0],   //( 1,  0,  0)
                [0, 0, 0],   //( 1,  0,  1)
                [0, 0, 0],   //( 1,  1, -1)
                [0, 0, 0],   //( 1,  1,  0)
                [0, 0, 0],   //( 1,  1,  1)
            ],
            // (-1,-1, 1)
            [
                [-1, -1, 0], //(-1, -1, -1)
                [-1, -1, 0], //(-1, -1,  0)
                [-1, -1, 1], //(-1, -1,  1)
                [-1, 0, 0],  //(-1,  0, -1)
                [-1, 0, 0],  //(-1,  0,  0)
                [-1, 0, 1],  //(-1,  0,  1)
                [-1, 0, 0],  //(-1,  1, -1)
                [-1, 0, 0],  //(-1,  1,  0)
                [-1, 0, 1],  //(-1,  1,  1)
                [0, -1, 0],  //( 0, -1, -1)
                [0, -1, 0],  //( 0, -1,  0)
                [0, -1, 1],  //( 0, -1,  1)
                [0, 0, 0],   //( 0,  0, -1)
                [0, 0, 0],   //( 0,  0,  0)
                [0, 0, 1],   //( 0,  0,  1)
                [0, 0, 0],   //( 0,  1, -1)
                [0, 0, 0],   //( 0,  1,  0)
                [0, 0, 1],   //( 0,  1,  1)
                [0, -1, 0],  //( 1, -1, -1)
                [0, -1, 0],  //( 1, -1,  0)
                [0, -1, 1],  //( 1, -1,  1)
                [0, 0, 0],   //( 1,  0, -1)
                [0, 0, 0],   //( 1,  0,  0)
                [0, 0, 1],   //( 1,  0,  1)
                [0, 0, 0],   //( 1,  1, -1)
                [0, 0, 0],   //( 1,  1,  0)
                [0, 0, 1],   //( 1,  1,  1)
            ],
            // (-1,-1,-1)
            [
                [-1, -1, -1], //(-1, -1, -1)
                [-1, -1, 0],  //(-1, -1,  0)
                [-1, -1, 0],  //(-1, -1,  1)
                [-1, 0, -1],  //(-1,  0, -1)
                [-1, 0, 0],   //(-1,  0,  0)
                [-1, 0, 0],   //(-1,  0,  1)
                [-1, 0, -1],  //(-1,  1, -1)
                [-1, 0, 0],   //(-1,  1,  0)
                [-1, 0, 0],   //(-1,  1,  1)
                [0, -1, -1],  //( 0, -1, -1)
                [0, -1, 0],   //( 0, -1,  0)
                [0, -1, 0],   //( 0, -1,  1)
                [0, 0, -1],   //( 0,  0, -1)
                [0, 0, 0],    //( 0,  0,  0)
                [0, 0, 0],    //( 0,  0,  1)
                [0, 0, -1],   //( 0,  1, -1)
                [0, 0, 0],    //( 0,  1,  0)
                [0, 0, 0],    //( 0,  1,  1)
                [0, -1, -1],  //( 1, -1, -1)
                [0, -1, 0],   //( 1, -1,  0)
                [0, -1, 0],   //( 1, -1,  1)
                [0, 0, -1],   //( 1,  0, -1)
                [0, 0, 0],    //( 1,  0,  0)
                [0, 0, 0],    //( 1,  0,  1)
                [0, 0, -1],   //( 1,  1, -1)
                [0, 0, 0],    //( 1,  1,  0)
                [0, 0, 0],    //( 1,  1,  1)
            ],
        ];
        Self { di, index, pbc }
    }
}

/// Size of the density data in 3d
pub struct Size {
    /// Number of voxels in the x-direction.
    pub x: isize,
    /// Number of voxels in the y-direction.
    pub y: isize,
    /// Number of voxels in the z-direction.
    pub z: isize,
    /// Total number of voxels.
    pub total: usize,
}

impl Size {
    /// The length of the flattened array for the density data in 3d
    fn new(x: usize, y: usize, z: usize) -> Self {
        let x = x as isize;
        let y = y as isize;
        let z = z as isize;
        let total = match x.checked_mul(y) {
            Some(xy) => match xy.checked_mul(z) {
                Some(xyz) => xyz as usize,
                None => panic!("Grid size is too large!"),
            },
            None => panic!("Grid size is too large!"),
        };
        // VoxelMap stores the maxima of weighs as -2 - maxima the largest
        // value maxima can be is Size::total so must check
        if total as isize + 1 == isize::MAX {
            panic!("Grid size is too large!");
        }
        Self { x, y, z, total }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_new() {
        let lattice = Lattice::new([[3., 3., 0.], [-3., 3., 0.], [1., 1., 1.]]);
        let grid = Grid::new([4, 4, 4], lattice.to_cartesian, [0., 0., 0.0]);
        assert!(
            (grid.voxel_lattice.volume - (lattice.volume / 64.0f64)).abs()
                < f64::EPSILON
        )
    }

    #[test]
    #[should_panic]
    fn grid_new_bad_grid() {
        let lattice = Lattice::new([[3., 3., 0.], [-3., 3., 0.], [1., 1., 1.]]);
        let _ = Grid::new([1, 4, 4], lattice.to_cartesian, [0., 0., 0.0]);
    }

    #[test]
    fn grid_full_shift() {
        let lattice = Lattice::new([[3., 3., 0.], [-3., 3., 0.], [1., 1., 1.]]);
        let grid = Grid::new([3, 4, 5], lattice.to_cartesian, [0., 0., 0.0]);
        let shift = [
            -26, -25, -24, -21, -20, -19, -16, -15, -14, -6, -5, -4, -1, 1, 4,
            5, 6, 14, 15, 16, 19, 20, 21, 24, 25, 26,
        ];
        assert_eq!(shift, grid.full_shift(26))
    }

    #[test]
    fn grid_reduced_shift() {
        let lattice = Lattice::new([[3., 3., 0.], [-3., 3., 0.], [1., 1., 1.]]);
        let grid = Grid::new([3, 4, 5], lattice.to_cartesian, [0., 0., 0.0]);
        let shift = [20, -20, 5, -5, 1, -1];
        assert_eq!(shift, grid.reduced_shift(26))
    }

    #[test]
    fn grid_gradient_shift() {
        let lattice = Lattice::new([[3., 3., 0.], [-3., 3., 0.], [1., 1., 1.]]);
        let grid = Grid::new([3, 4, 5], lattice.to_cartesian, [0., 0., 0.0]);
        assert_eq!(1, grid.gradient_shift(26, [0., 0., 1.]))
    }

    #[test]
    fn shift_index_gen() {
        let index = Shift::index_gen(&Size::new(3, 4, 5));
        let test_index: [usize; 27] = [
            0, 1, 4, 5, 6, 9, 15, 16, 19, 20, 21, 24, 25, 26, 29, 35, 36, 39,
            40, 41, 44, 45, 46, 49, 55, 56, 59,
        ];
        let test_value: [u8; 27] = [
            26, 24, 25, 20, 18, 19, 23, 21, 22, 8, 6, 7, 2, 0, 1, 5, 3, 4, 17,
            15, 16, 11, 9, 10, 14, 12, 13,
        ];
        for (i, v) in test_index.iter().zip(&test_value) {
            assert_eq!(index[*i], *v)
        }
    }

    #[test]
    fn shift_new() {
        let shift = Shift::new(&Size::new(3, 4, 5));
        let di = [
            [
                -26, -25, -24, -21, -20, -19, -16, -15, -14, -6, -5, -4, -1, 0,
                1, 4, 5, 6, 14, 15, 16, 19, 20, 21, 24, 25, 26,
            ],
            [
                -26, -25, -29, -21, -20, -24, -16, -15, -19, -6, -5, -9, -1, 0,
                -4, 4, 5, 1, 14, 15, 11, 19, 20, 16, 24, 25, 21,
            ],
            [
                -21, -25, -24, -16, -20, -19, -11, -15, -14, -1, -5, -4, 4, 0,
                1, 9, 5, 6, 19, 15, 16, 24, 20, 21, 29, 25, 26,
            ],
            [
                -26, -25, -24, -21, -20, -19, -36, -35, -34, -6, -5, -4, -1, 0,
                1, -16, -15, -14, 14, 15, 16, 19, 20, 21, 4, 5, 6,
            ],
            [
                -26, -25, -29, -21, -20, -24, -36, -35, -39, -6, -5, -9, -1, 0,
                -4, -16, -15, -19, 14, 15, 11, 19, 20, 16, 4, 5, 1,
            ],
            [
                -21, -25, -24, -16, -20, -19, -31, -35, -34, -1, -5, -4, 4, 0,
                1, -11, -15, -14, 19, 15, 16, 24, 20, 21, 9, 5, 6,
            ],
            [
                -6, -5, -4, -21, -20, -19, -16, -15, -14, 14, 15, 16, -1, 0, 1,
                4, 5, 6, 34, 35, 36, 19, 20, 21, 24, 25, 26,
            ],
            [
                -6, -5, -9, -21, -20, -24, -16, -15, -19, 14, 15, 11, -1, 0,
                -4, 4, 5, 1, 34, 35, 31, 19, 20, 16, 24, 25, 21,
            ],
            [
                -1, -5, -4, -16, -20, -19, -11, -15, -14, 19, 15, 16, 4, 0, 1,
                9, 5, 6, 39, 35, 36, 24, 20, 21, 29, 25, 26,
            ],
            [
                -26, -25, -24, -21, -20, -19, -16, -15, -14, -6, -5, -4, -1, 0,
                1, 4, 5, 6, -46, -45, -44, -41, -40, -39, -36, -35, -34,
            ],
            [
                -26, -25, -29, -21, -20, -24, -16, -15, -19, -6, -5, -9, -1, 0,
                -4, 4, 5, 1, -46, -45, -49, -41, -40, -44, -36, -35, -39,
            ],
            [
                -21, -25, -24, -16, -20, -19, -11, -15, -14, -1, -5, -4, 4, 0,
                1, 9, 5, 6, -41, -45, -44, -36, -40, -39, -31, -35, -34,
            ],
            [
                -26, -25, -24, -21, -20, -19, -36, -35, -34, -6, -5, -4, -1, 0,
                1, -16, -15, -14, -46, -45, -44, -41, -40, -39, -56, -55, -54,
            ],
            [
                -26, -25, -29, -21, -20, -24, -36, -35, -39, -6, -5, -9, -1, 0,
                -4, -16, -15, -19, -46, -45, -49, -41, -40, -44, -56, -55, -59,
            ],
            [
                -21, -25, -24, -16, -20, -19, -31, -35, -34, -1, -5, -4, 4, 0,
                1, -11, -15, -14, -41, -45, -44, -36, -40, -39, -51, -55, -54,
            ],
            [
                -6, -5, -4, -21, -20, -19, -16, -15, -14, 14, 15, 16, -1, 0, 1,
                4, 5, 6, -26, -25, -24, -41, -40, -39, -36, -35, -34,
            ],
            [
                -6, -5, -9, -21, -20, -24, -16, -15, -19, 14, 15, 11, -1, 0,
                -4, 4, 5, 1, -26, -25, -29, -41, -40, -44, -36, -35, -39,
            ],
            [
                -1, -5, -4, -16, -20, -19, -11, -15, -14, 19, 15, 16, 4, 0, 1,
                9, 5, 6, -21, -25, -24, -36, -40, -39, -31, -35, -34,
            ],
            [
                34, 35, 36, 39, 40, 41, 44, 45, 46, -6, -5, -4, -1, 0, 1, 4, 5,
                6, 14, 15, 16, 19, 20, 21, 24, 25, 26,
            ],
            [
                34, 35, 31, 39, 40, 36, 44, 45, 41, -6, -5, -9, -1, 0, -4, 4,
                5, 1, 14, 15, 11, 19, 20, 16, 24, 25, 21,
            ],
            [
                39, 35, 36, 44, 40, 41, 49, 45, 46, -1, -5, -4, 4, 0, 1, 9, 5,
                6, 19, 15, 16, 24, 20, 21, 29, 25, 26,
            ],
            [
                34, 35, 36, 39, 40, 41, 24, 25, 26, -6, -5, -4, -1, 0, 1, -16,
                -15, -14, 14, 15, 16, 19, 20, 21, 4, 5, 6,
            ],
            [
                34, 35, 31, 39, 40, 36, 24, 25, 21, -6, -5, -9, -1, 0, -4, -16,
                -15, -19, 14, 15, 11, 19, 20, 16, 4, 5, 1,
            ],
            [
                39, 35, 36, 44, 40, 41, 29, 25, 26, -1, -5, -4, 4, 0, 1, -11,
                -15, -14, 19, 15, 16, 24, 20, 21, 9, 5, 6,
            ],
            [
                54, 55, 56, 39, 40, 41, 44, 45, 46, 14, 15, 16, -1, 0, 1, 4, 5,
                6, 34, 35, 36, 19, 20, 21, 24, 25, 26,
            ],
            [
                54, 55, 51, 39, 40, 36, 44, 45, 41, 14, 15, 11, -1, 0, -4, 4,
                5, 1, 34, 35, 31, 19, 20, 16, 24, 25, 21,
            ],
            [
                59, 55, 56, 44, 40, 41, 49, 45, 46, 19, 15, 16, 4, 0, 1, 9, 5,
                6, 39, 35, 36, 24, 20, 21, 29, 25, 26,
            ],
        ];
        assert_eq!(shift.di, di)
    }

    #[test]
    fn shift_get() {
        assert_eq!(Shift::new(&Size::new(3, 4, 5)).get(25).0[0], -21)
    }

    #[test]
    fn size_new() {
        assert_eq!(Size::new(3, 4, 5).total, 3 * 4 * 5)
    }

    #[test]
    #[should_panic]
    #[cfg(target_pointer_width = "32")]
    fn size_too_large_32() {
        let _ = Size::new(1300, 1300, 1300);
    }

    #[test]
    #[should_panic]
    fn size_too_large() {
        let _ = Size::new(2100000, 2100000, 2100000);
    }
}
