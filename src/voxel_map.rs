use crate::grid::Grid;
use rustc_hash::{FxHashMap, FxHashSet};
use std::mem::MaybeUninit;
use std::ops::{Add, Sub};
use std::sync::Arc;
use std::sync::atomic::{AtomicIsize, AtomicUsize, Ordering};

/// An [i4; 3] stored as a single u16
#[derive(Clone, Copy, Debug)]
pub struct EncodedImage(u16);

impl EncodedImage {
    const BIAS: i8 = 8;
    const BITS: usize = 4;
    const MASK: u16 = 2u16.pow(Self::BITS as u32) - 1; // 4 bits
    const ZERO: u16 = ((Self::BIAS as u16) & Self::MASK)
        | (((Self::BIAS as u16) & Self::MASK) << Self::BITS)
        | (((Self::BIAS as u16) & Self::MASK) << (2 * Self::BITS));

    /// Create a new encoded image from an [i8; 3]
    pub fn new(image: [i8; 3]) -> Self {
        let mut encoded: u16 = 0;
        image.iter().enumerate().for_each(|(i, img)| {
            let biased = (img + Self::BIAS) as u16;
            debug_assert!(biased <= Self::MASK, "Image out of encoding range");
            encoded |= (biased & Self::MASK) << (i * Self::BITS);
        });
        Self(encoded)
    }

    /// Decode to an [i8; 3]
    pub fn decode(self) -> [i8; 3] {
        let mut image = [0; 3];
        image.iter_mut().enumerate().for_each(|(i, img)| {
            let biased = (self.0 >> (i * Self::BITS)) & Self::MASK;
            *img = (biased as i8) - Self::BIAS;
        });
        return image;
    }

    pub fn image_add(self, b: [i8; 3]) -> Self {
        let a = self.decode();
        Self::new([a[0] + b[0], a[1] + b[1], a[2] + b[2]])
    }

    pub fn is_zero(&self) -> bool {
        self.0 == Self::ZERO
    }
}

impl Add for EncodedImage {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let a = self.decode();
        let b = rhs.decode();
        Self::new([a[0] + b[0], a[1] + b[1], a[2] + b[2]])
    }
}
impl Sub for EncodedImage {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let a = self.decode();
        let b = rhs.decode();
        Self::new([a[0] - b[0], a[1] - b[1], a[2] - b[2]])
    }
}

/// An encoded atom number and image.
///
/// [Atom Number: 20 bits] | [Image: 12 bits]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
#[repr(transparent)]
pub struct EncodedAtom(pub u32);

impl EncodedAtom {
    const SHIFT: usize = EncodedImage::BITS * 3;
    const MASK: u32 = 2u32.pow(Self::SHIFT as u32) - 1;
    const BITS: u32 = 20;
    const MAX_ATOM: u32 = 2u32.pow(Self::BITS) - 1;

    pub fn new(atom: u32, image: EncodedImage) -> Self {
        debug_assert!(
            atom < Self::MAX_ATOM,
            "Atom Number out of EncodedAtom Range"
        );
        Self((atom << Self::SHIFT) | (image.0) as u32)
    }

    pub fn atom_index(&self) -> u32 {
        self.0 >> Self::SHIFT
    }

    pub fn image(&self) -> EncodedImage {
        EncodedImage((self.0 & Self::MASK) as u16)
    }

    pub fn image_add(self, image: EncodedImage) -> Self {
        Self::new(self.atom_index(), self.image() + image)
    }

    pub fn image_sub(self, image: EncodedImage) -> Self {
        Self::new(self.atom_index(), self.image() - image)
    }

    pub fn decode_partial(self) -> (u32, EncodedImage) {
        (self.atom_index(), self.image())
    }

    pub fn decode_full(self) -> (u32, [i8; 3]) {
        (self.atom_index(), self.image().decode())
    }
}

#[derive(Clone, Copy, Debug)]
pub struct EncodedWeight(u64);

impl EncodedWeight {
    pub fn new(encoded_atom: EncodedAtom, weight: f32) -> Self {
        Self((encoded_atom.0 as u64) | ((weight.to_bits() as u64) << 32))
    }

    pub fn decode(self) -> (EncodedAtom, f32) {
        (
            EncodedAtom(self.0 as u32),
            f32::from_bits((self.0 >> 32) as u32),
        )
    }
}

/// Describes the state of the voxel.
pub enum Voxel {
    /// Contians the position of the voxel's maxima.
    Maxima(EncodedAtom),
    /// Contians a vector of the maxima the current voxel contributes to and
    /// their weights.
    Boundary(FxHashMap<EncodedAtom, f32>),
    /// A voxel beneath the vacuum tolerance and not contributing to any maxima.
    Vacuum,
}

/// A structure for building and processing the map between voxel and maxima.
/// Bader maxima are stored in the voxel_map whilst the contributing weights are
/// stored in the weight_map. The weight_map is only written to once by each
/// point and so once a value has been written it is safe to read by any thread.
/// To check it has been written to `weight_get` monitors the state of corresponding
/// voxel_map value. Writing to the map is acheived by acquiring the lock, noting
/// the length of the weight_map, pushing the weight vector for voxel p to the
/// weight_map, droping the write lock and then storing the index of the inserted
/// vector using `weight_store`.
///
/// # Examples
/// ```
/// use bader::voxel_map::BlockingVoxelMap;
///
/// for p in 0..1isize {
///     let voxel_map = BlockingVoxelMap::new(
///         [2, 5, 2],
///         [[2.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 2.0]],
///         [0.0, 0.0, 0.0],
///     );
///     let i = {
///         let mut weight = voxel_map.lock();
///         (*weight).push(Vec::with_capacity(0).into());
///         weight.len() - 1
///     };
///     voxel_map.weight_store(p, i)
/// }
/// ```
pub struct BlockingVoxelMap {
    weight_map: Arc<[MaybeUninit<Box<[EncodedWeight]>>]>,
    voxel_map: Arc<[AtomicIsize]>,
    pub grid: Grid,
    weight_counter: AtomicUsize,
}

impl BlockingVoxelMap {
    /// Initialises a [`BlockingVoxelMap`] and the [`Grid`] that will faciliate movemoment around the
    /// map.
    pub fn new(
        grid: [usize; 3],
        lattice: [[f64; 3]; 3],
        voxel_origin: [f64; 3],
    ) -> Self {
        let grid = Grid::new(grid, lattice, voxel_origin);
        let size = grid.size.total;
        // For mapping the the voxels
        let mut weight_map = Vec::with_capacity(size);
        weight_map.resize_with(size, MaybeUninit::uninit);
        let weight_map = Arc::from(weight_map.into_boxed_slice());
        let mut voxel_map = Vec::with_capacity(size);
        voxel_map.resize_with(size, || AtomicIsize::new(-1));
        let voxel_map = Arc::from(voxel_map.into_boxed_slice());
        let weight_counter = AtomicUsize::new(0);
        // For post processing
        Self {
            weight_map,
            voxel_map,
            grid,
            weight_counter,
        }
    }

    /// Retrieves the state of the voxel, p. This will lock until p has been stored
    /// in the VoxelMap and then return either a `Voxel::Maxima` or `Voxel::Weight`.
    /// Calling this on a voxel, p, that is below the vacuum_tolerance will deadlock
    /// as a voxel is considered stored once voxel_map\[p\] > -1.
    pub fn weight_get(&self, i: isize) -> FxHashMap<EncodedAtom, f32> {
        let i = -2 - i;
        (unsafe { self.weight_map.get_unchecked(i as usize).assume_init_ref() })
            .iter()
            .map(|u| u.decode())
            .collect()
    }

    /// Atomic loading of voxel, p, from voxel_map blocks if maxima == -1
    pub fn maxima_get(&self, p: isize) -> isize {
        loop {
            match self.voxel_map[p as usize].load(Ordering::Acquire) {
                -1 => std::thread::yield_now(),
                x => break x,
            }
        }
    }
    ///
    /// Atomic loading of voxel, p, from voxel_map blocks if maxima == -1
    pub fn voxel_get(&self, p: isize) -> Voxel {
        let i = self.maxima_get(p);
        match i.cmp(&-1) {
            std::cmp::Ordering::Less => Voxel::Boundary(self.weight_get(i)),
            std::cmp::Ordering::Equal => Voxel::Vacuum,
            std::cmp::Ordering::Greater => Voxel::Maxima(EncodedAtom(i as u32)),
        }
    }

    /// Check if a maxima is stored
    pub fn maxima_check(&self, p: isize) -> Option<isize> {
        match self.voxel_map[p as usize].load(Ordering::Relaxed) {
            -1 => None,
            x => Some(x),
        }
    }

    /// Stores the maxima of voxel, p, in the voxel_map.
    pub fn maxima_store(&self, p: isize, maxima: isize) {
        self.voxel_map[p as usize].store(maxima, Ordering::Release);
    }

    /// Stores the index of p's weight contributions in weight_map into the
    /// weight_index.
    pub fn weight_store(&self, p: isize, weights: Box<[EncodedWeight]>) {
        let i = self.weight_counter.fetch_add(1, Ordering::Relaxed);
        unsafe {
            let ptr: *mut Box<[EncodedWeight]> =
                self.weight_map.get_unchecked(i) as *const _ as *mut _;
            ptr.write(weights)
        }
        self.maxima_store(p, -2 - (i as isize));
    }

    /// Extract the voxel map data.
    pub fn into_inner(self) -> (Vec<isize>, Vec<Box<[EncodedWeight]>>, Grid) {
        (
            self.voxel_map
                .iter()
                .map(|x| x.load(Ordering::Relaxed))
                .collect(),
            self.weight_map
                .iter()
                .take(self.weight_counter.into_inner())
                .map(|mu| unsafe { mu.assume_init_read() })
                .collect(),
            self.grid,
        )
    }
}

/// A VoxelMap for if the maxima stored are atomic indices.
pub struct VoxelMap {
    /// The vector mapping the voxel to a maxima.
    pub voxel_map: Vec<isize>,
    /// The vector containing the weights for boundary voxels.
    pub weight_map: Vec<Box<[EncodedWeight]>>,
    /// The Grid used to navigate the VoxelMap.
    pub grid: Grid,
}

impl VoxelMap {
    /// Create a new [`VoxelMap`]
    pub fn new(
        voxel_map: Vec<isize>,
        weight_map: Vec<Box<[EncodedWeight]>>,
        grid: Grid,
    ) -> Self {
        Self {
            voxel_map,
            weight_map,
            grid,
        }
    }

    /// Create a new [`VoxelMap`] from a [`BlockingVoxelMap`].
    pub fn from_blocking_voxel_map(voxel_map: BlockingVoxelMap) -> Self {
        let (voxel_map, weight_map, grid) = voxel_map.into_inner();
        Self::new(voxel_map, weight_map, grid)
    }

    /// Produce an Iter over the boundary voxels.
    pub fn weight_iter(&self) -> std::slice::Iter<'_, Box<[EncodedWeight]>> {
        self.weight_map.iter()
    }

    /// Get the length of the weight_map.
    pub fn weight_len(&self) -> usize {
        self.weight_map.len()
    }

    /// Get a refernce to the grid used by the VoxelMap.
    pub fn grid_get(&self) -> &Grid {
        &self.grid
    }

    /// Returns the atom associated with the point.
    pub fn maxima_to_atom(&self, maxima: usize) -> usize {
        maxima
    }

    /// Retrieval of the state of the voxel, p.
    pub fn maxima_to_voxel(&self, maxima: isize) -> Voxel {
        match maxima.cmp(&-1) {
            std::cmp::Ordering::Equal => Voxel::Vacuum,
            std::cmp::Ordering::Greater => {
                Voxel::Maxima(EncodedAtom(maxima as u32))
            }
            std::cmp::Ordering::Less => {
                Voxel::Boundary(self.maxima_to_weight(maxima))
            }
        }
    }

    /// Return a reference to the weights from the given maxima, Note: maxima here must be < -1.
    pub fn maxima_to_weight(
        &self,
        maxima: isize,
    ) -> FxHashMap<EncodedAtom, f32> {
        self.weight_map[(-2 - maxima) as usize]
            .iter()
            .map(|ed| ed.decode())
            .collect()
    }

    /// Return an Iter over the maxima stored in the VoxelMap.
    pub fn maxima_iter(&self) -> std::slice::Iter<'_, isize> {
        self.voxel_map.iter()
    }

    /// Get the length of the voxel_map.
    pub fn maxima_len(&self) -> usize {
        self.voxel_map.len()
    }

    /// Return a Chunk over the maxima stored in the VoxelMap.
    pub fn maxima_chunks(
        &self,
        chunk_size: usize,
    ) -> std::slice::Chunks<'_, isize> {
        self.voxel_map.chunks(chunk_size)
    }

    /// Retrieval of the state of the voxel, p.
    pub fn voxel_get(&self, p: isize) -> Voxel {
        self.maxima_to_voxel(self.maxima_get(p))
    }

    /// Return the stored maxima at point p.
    pub fn maxima_get(&self, p: isize) -> isize {
        self.voxel_map[p as usize]
    }

    /// Produce a mask for a specific volume number.
    pub fn volume_map(&self, volume_number: isize) -> Vec<Option<f64>> {
        self.maxima_iter()
            .map(|maxima| {
                if *maxima == volume_number {
                    Some(1.0)
                } else if *maxima < -1 {
                    let mut w = None;
                    for (m, weight) in
                        self.maxima_to_weight(*maxima).into_iter()
                    {
                        if (m.atom_index() as isize) == volume_number {
                            w = Some(weight as f64);
                            break;
                        }
                    }
                    w
                } else {
                    None
                }
            })
            .collect()
    }
    /// Produce a mask for a collection volume numbers.
    pub fn multi_volume_map(
        &self,
        volume_numbers: &FxHashSet<isize>,
    ) -> Vec<Option<f64>> {
        self.maxima_iter()
            .map(|maxima| {
                if volume_numbers.contains(maxima) {
                    Some(1.0)
                } else if *maxima < -1 {
                    let mut w = 0.0;
                    for (m, weight) in
                        self.maxima_to_weight(*maxima).into_iter()
                    {
                        if volume_numbers.contains(&(m.atom_index() as isize)) {
                            w += weight as f64;
                        }
                    }
                    Some(w)
                } else {
                    None
                }
            })
            .collect()
    }
}
