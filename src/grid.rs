use crate::atoms::Lattice;
use crate::utils::dot;
use crate::voronoi::Voronoi;
use crate::voxel_map::EncodedImage;

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
    ///
    /// # Panics
    /// Panics if any dimension of `grid` is less than 2 or the total size would overflow.
    pub fn new(
        grid: [usize; 3],
        lattice: [[f64; 3]; 3],
        voxel_origin: [f64; 3],
    ) -> Self {
        let size = Size::new(grid[0], grid[1], grid[2]);
        let shift = Shift::new(&size);
        let mut voxel_lattice = lattice;
        for (vector, n) in voxel_lattice.iter_mut().zip(grid) {
            vector.iter_mut().for_each(|x| *x /= n as f64);
        }
        let voxel_lattice = Lattice::new(voxel_lattice);
        let voronoi = Voronoi::new(&voxel_lattice);
        Self {
            shift,
            size,
            voronoi,
            voxel_lattice,
            voxel_origin,
        }
    }

    /// Iterates over the neighbours of `p` that share a Voronoi face with it.
    ///
    /// Each item is `((neighbour, image), alpha)` where `image` is the periodic
    /// image the neighbour was found in relative to `p` and `alpha` the Voronoi
    /// flux coefficient of the shared face.
    pub fn voronoi_shifts(
        &self,
        p: isize,
    ) -> impl Iterator<Item = ((isize, EncodedImage), f64)> + '_ {
        self.voronoi
            .vectors
            .iter()
            .map(move |steps| {
                steps.iter().fold(
                    (p, EncodedImage::new([0; 3])),
                    |(pn, image), step| {
                        let (shift, pbc) = self.shift.get(pn);
                        (pn + shift[*step], image.image_add(pbc[*step]))
                    },
                )
            })
            .zip(self.voronoi.alphas.iter().copied())
    }

    /// As [`Grid::voronoi_shifts`] but without tracking the periodic image.
    pub fn voronoi_shifts_nocheck(
        &self,
        p: isize,
    ) -> impl Iterator<Item = (isize, f64)> + '_ {
        self.voronoi
            .vectors
            .iter()
            .map(move |steps| {
                steps
                    .iter()
                    .fold(p, |pn, step| pn + self.shift.get(pn).0[*step])
            })
            .zip(self.voronoi.alphas.iter().copied())
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
}

/// Pre-computed index offsets for stepping to the 27 surrounding voxels.
///
/// A voxel can be in one of 27 "edge classes" depending on whether it sits
/// on the lower edge, in the bulk or on the upper edge of each axis. For each
/// class `di` holds the flat-index offset to every neighbour (wrapping
/// periodically where needed) and `pbc` the periodic image crossed by the step.
///
/// The outer (class) index is ordered with `[bulk, upper edge, lower edge]`
/// per axis so that `class = cx * 9 + cy * 3 + cz`:
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
/// The inner (neighbour) index is ordered `[-1, 0, 1]` per axis so that
/// `neighbour = (dx + 1) * 9 + (dy + 1) * 3 + (dz + 1)`, matching
/// [`Lattice`]'s shift matrices:
/// <pre class="rust">
/// 0 -> (-1,-1,-1)   7 -> (-1, 1, 0)  14 -> (0, 0, 1)  21 -> (1, 0,-1)
/// 1 -> (-1,-1, 0)   8 -> (-1, 1, 1)  15 -> (0, 1,-1)  22 -> (1, 0, 0)
/// 2 -> (-1,-1, 1)   9 -> (0,-1,-1)   16 -> (0, 1, 0)  23 -> (1, 0, 1)
/// 3 -> (-1, 0,-1)  10 -> (0,-1, 0)   17 -> (0, 1, 1)  24 -> (1, 1,-1)
/// 4 -> (-1, 0, 0)  11 -> (0,-1, 1)   18 -> (1,-1,-1)  25 -> (1, 1, 0)
/// 5 -> (-1, 0, 1)  12 -> (0, 0,-1)   19 -> (1,-1, 0)  26 -> (1, 1, 1)
/// 6 -> (-1, 1,-1)  13 -> (0, 0, 0)   20 -> (1,-1, 1)
/// </pre>
struct Shift {
    /// Flat-index offsets for each edge class.
    di: [[isize; 27]; 27],
    /// Periodic image crossed by each step for each edge class.
    pbc: [[[i8; 3]; 27]; 27],
    /// The edge class of every voxel.
    index: Vec<u8>,
}

impl Shift {
    /// Gets the offsets and periodic images for the neighbours of voxel `i`.
    fn get(&self, i: isize) -> (&[isize; 27], &[[i8; 3]; 27]) {
        let class = self.index[i as usize] as usize;
        (&self.di[class], &self.pbc[class])
    }

    /// The edge class of a coordinate along an axis of length `n`:
    /// 0 for the bulk, 1 for the upper edge and 2 for the lower edge.
    fn edge_class(i: isize, n: isize) -> u8 {
        if i == 0 {
            2
        } else if i == n - 1 {
            1
        } else {
            0
        }
    }

    /// Generate the edge class of every voxel.
    fn index_gen(size: &Size) -> Vec<u8> {
        let mut index = Vec::with_capacity(size.total);
        for x in 0..size.x {
            let cx = Self::edge_class(x, size.x) * 9;
            for y in 0..size.y {
                let cy = Self::edge_class(y, size.y) * 3;
                for z in 0..size.z {
                    index.push(cx + cy + Self::edge_class(z, size.z));
                }
            }
        }
        index
    }

    /// Creates the offset tables for every edge class.
    fn new(size: &Size) -> Self {
        let dims = [size.x, size.y, size.z];
        let strides = [size.y * size.z, size.z, 1];
        let mut di = [[0isize; 27]; 27];
        let mut pbc = [[[0i8; 3]; 27]; 27];
        for (class, (di, pbc)) in di.iter_mut().zip(pbc.iter_mut()).enumerate()
        {
            // the side of each axis this class sits on: 0 bulk, 1 upper, -1 lower
            let edges = [class / 9, (class / 3) % 3, class % 3]
                .map(|c| [0, 1, -1][c] as isize);
            for (neighbour, (di, pbc)) in
                di.iter_mut().zip(pbc.iter_mut()).enumerate()
            {
                let steps = [neighbour / 9, (neighbour / 3) % 3, neighbour % 3]
                    .map(|s| s as isize - 1);
                for axis in 0..3 {
                    // stepping off the edge we sit on wraps around the cell
                    let wrap = if steps[axis] == edges[axis] {
                        steps[axis]
                    } else {
                        0
                    };
                    *di += (steps[axis] - wrap * dims[axis]) * strides[axis];
                    pbc[axis] = wrap as i8;
                }
            }
        }
        Self {
            di,
            pbc,
            index: Shift::index_gen(size),
        }
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
    ///
    /// # Panics
    /// Panics if any dimension is less than 2 or the total size is too large
    /// to be stored in the voxel map.
    fn new(x: usize, y: usize, z: usize) -> Self {
        if x < 2 || y < 2 || z < 2 {
            panic!("Grid must be at least 2 voxels in every direction!");
        }
        // VoxelMap stores the index of weights as -2 - index so the total
        // must leave room for this encoding in an isize.
        let total = x
            .checked_mul(y)
            .and_then(|xy| xy.checked_mul(z))
            .filter(|&t| t < isize::MAX as usize - 1)
            .unwrap_or_else(|| panic!("Grid size is too large!"));
        Self {
            x: x as isize,
            y: y as isize,
            z: z as isize,
            total,
        }
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
    fn grid_to_3d_and_cartesian() {
        let grid = Grid::new(
            [3, 4, 5],
            [[3., 0., 0.], [0., 4., 0.], [0., 0., 5.]],
            [0.5, 0.5, 0.5],
        );
        // index = x * 20 + y * 5 + z
        assert_eq!(grid.to_3d(2 * 20 + 3 * 5 + 4), [2, 3, 4]);
        assert_eq!(grid.to_3d(0), [0, 0, 0]);
        // each voxel is 1 x 1 x 1 so cartesian = index + origin
        assert_eq!(grid.to_cartesian(2 * 20 + 3 * 5 + 4), [2.5, 3.5, 4.5]);
    }

    #[test]
    fn grid_voronoi_shifts_bulk() {
        // cubic lattice so the Voronoi neighbours are the 6 nearest voxels
        let grid = Grid::new(
            [4, 4, 4],
            [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]],
            [0.; 3],
        );
        // voxel (1, 1, 1) is in the bulk so no images are crossed
        let p = 16 + 4 + 1;
        let mut neighbours = grid
            .voronoi_shifts(p)
            .map(|((pt, image), alpha)| {
                assert!(image.is_zero());
                assert!((alpha - 1.0).abs() < 1e-10);
                pt
            })
            .collect::<Vec<isize>>();
        neighbours.sort_unstable();
        assert_eq!(
            neighbours,
            vec![p - 16, p - 4, p - 1, p + 1, p + 4, p + 16]
        );
    }

    #[test]
    fn grid_voronoi_shifts_periodic() {
        let grid = Grid::new(
            [4, 4, 4],
            [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]],
            [0.; 3],
        );
        // voxel (0, 0, 0) wraps in the negative direction on every axis
        let mut shifts = grid
            .voronoi_shifts(0)
            .map(|((pt, image), _)| (pt, image.decode()))
            .collect::<Vec<_>>();
        shifts.sort_unstable();
        assert_eq!(
            shifts,
            vec![
                (1, [0, 0, 0]),
                (3, [0, 0, -1]),
                (4, [0, 0, 0]),
                (12, [0, -1, 0]),
                (16, [0, 0, 0]),
                (48, [-1, 0, 0]),
            ]
        );
        // the nocheck variant visits the same voxels
        let mut nocheck = grid
            .voronoi_shifts_nocheck(0)
            .map(|(pt, _)| pt)
            .collect::<Vec<isize>>();
        nocheck.sort_unstable();
        assert_eq!(nocheck, vec![1, 3, 4, 12, 16, 48]);
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
        let shift = Shift::new(&Size::new(3, 4, 5));
        // voxel 25 = (1, 1, 0) is on the lower z edge
        let (di, pbc) = shift.get(25);
        assert_eq!(di[0], -21);
        assert_eq!(pbc[0], [0, 0, -1]);
        assert_eq!(pbc[13], [0, 0, 0]);
        assert_eq!(pbc[14], [0, 0, 0]);
        // voxel 59 = (2, 3, 4) is on the upper edge of every axis
        let (di, pbc) = shift.get(59);
        assert_eq!(di[26], -59);
        assert_eq!(pbc[26], [1, 1, 1]);
        assert_eq!(pbc[0], [0, 0, 0]);
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
