use crate::grid::Grid;
use rustc_hash::{FxHashMap, FxHashSet};
use std::mem::MaybeUninit;
use std::ops::{Add, Sub};
use std::sync::Arc;
use std::sync::atomic::{AtomicIsize, AtomicUsize, Ordering};

/// A compressed 3D image offset stored as a single `u16`.
///
/// This structure packs three 4-bit integers into 16 bits using a bias of 8.
/// This allows for vector components in the range `[-8, 7]`.
///
/// # Layout
/// * **Bits 0-3**: X component
/// * **Bits 4-7**: Y component
/// * **Bits 8-11**: Z component
/// * **Bias**: +8 per component
///
/// # Examples
/// ```
/// use bader::voxel_map::EncodedImage;
///
/// let vec = [1, -1, 0];
/// let encoded = EncodedImage::new(vec);
/// assert_eq!(encoded.decode(), vec);
/// ```
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
        image
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

/// A packed identifier containing an Atom ID and an Image Offset.
///
/// Stored as a `u32` to efficiently handle around 1 million atoms.
///
/// # Layout
/// | Component | Bits | Description |
/// |-----------|------|-------------|
/// | Atom ID   | 20   | Unique identifier for the atom (max ~1 million) |
/// | Image     | 12   | Encoded image offset (via [`EncodedImage`]) |
///
/// # Examples
/// ```
/// use bader::voxel_map::{EncodedAtom, EncodedImage};
///
/// let atom_id = 42;
/// let img = EncodedImage::new([0, 0, 1]);
/// let encoded = EncodedAtom::new(atom_id, img);
///
/// assert_eq!(encoded.atom_index(), 42);
/// ```
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

    pub fn new_zero_image(atom: u32) -> Self {
        Self((atom << Self::SHIFT) | (EncodedImage::ZERO) as u32)
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

/// A compact representation of an atom and its associated weight.
///
/// This structure packs an [`EncodedAtom`] and a `f32` weight into a single `u64`
/// to minimise memory usage when storing millions of boundary weights.
///
/// # Layout
/// | Component | Bits | Description |
/// |-----------|------|-------------|
/// | Weight    | 32   | The weight as an IEEE 754 `f32` (stored in high bits) |
/// | Atom      | 32   | The [`EncodedAtom`] identifier (stored in low bits) |
///
/// # Examples
/// ```
/// use bader::voxel_map::{EncodedWeight, EncodedAtom, EncodedImage};
///
/// let atom = EncodedAtom::new(42, EncodedImage::new([0, 0, 0]));
/// let weight_val = 0.5f32;
///
/// // Pack
/// let encoded = EncodedWeight::new(atom, weight_val);
///
/// // Unpack
/// let (decoded_atom, decoded_weight) = encoded.decode();
/// assert_eq!(decoded_atom, atom);
/// assert_eq!(decoded_weight, 0.5);
/// ```
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

/// A thread-safe, write-optimised map for populating Bader volumes.
///
/// Designed for concurrent generation of voxel assignments. Threads can safely
/// store "Maxima" (integer IDs) or "Weights" (boundary contributions) without
/// global locking.
///
/// # Storage Logic
/// The `voxel_map` stores an `AtomicIsize` for every voxel:
/// * **`>= 0`**: The voxel belongs entirely to the Atom with this ID (Maxima).
/// * **`-1`**: The voxel is Vacuum or currently processing.
/// * **`<-1`**: The voxel is on a boundary. The value is `-2 - index`, where `index`
///   points to a slice of weights in `weight_map`.
///
/// # Examples
/// ```
/// use bader::voxel_map::{BlockingVoxelMap, EncodedWeight, EncodedAtom, EncodedImage};
///
/// // 1. Init
/// let map = BlockingVoxelMap::new(
///     [2, 2, 2],
///     [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
///     [0.0, 0.0, 0.0],
/// );
/// let atom = EncodedAtom::new(100, EncodedImage::new([0,0,0]));
///
/// // 2. Store a Maxima
/// map.maxima_store(0, atom.0 as isize); // Voxel 0 -> Atom 100 Image [0, 0, 0]
///
/// // 3. Store Weights (Boundary)
/// let w = EncodedWeight::new(atom, 0.5);
/// map.weight_store(1, vec![w].into_boxed_slice());
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

    /// Stores the maxima of voxel, p, in the voxel_map. Note: maximas should be stored in their
    /// encoded form.
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

/// A read-optimised, non-blocking map for analysing Bader partitions.
///
/// While [`BlockingVoxelMap`] is designed for concurrent *write* operations during the
/// partitioning phase, `VoxelMap` is designed for efficient *read* operations during
/// the analysis phase. It provides methods to calculate partial volumes and integrate
/// properties over atoms.
///
/// # Structure
/// * **`voxel_map`**: A flat array matching the grid size.
///   * `i >= 0`: The voxel belongs to the Maxima (Atom) with index `i`.
///   * `i == -1`: Vacuum.
///   * `i < -1`: Boundary voxel. Points to weights at index `(-2 - i)`.
/// * **`weight_map`**: A collection of weights for boundary voxels.
///
/// # Examples
/// ```
/// use bader::voxel_map::{BlockingVoxelMap, VoxelMap};
///
/// // 1. Construct and populate a BlockingVoxelMap (concurrently)
/// let blocking = BlockingVoxelMap::new(
///     [2, 2, 2],
///     [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
///     [0.0, 0.0, 0.0],
/// );
/// // ... (spawn threads to populate map) ...
///
/// // 2. Convert to VoxelMap for analysis (consumes the blocking map)
/// let map = VoxelMap::from_blocking_voxel_map(blocking);
///
/// // 3. ... (Perform charge summing and critical point analysis) ...
///
/// // 4. Analyse volumes
/// let atom_volume = map.volume_map(0); // Get contributions for Atom 0
/// ```
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

#[cfg(test)]
mod tests {
    use super::*;
    use rustc_hash::FxHashSet;

    // --- EncodedImage Tests ---

    #[test]
    fn test_encoded_image_packing() {
        // Test bounds and zero
        let zero = EncodedImage::new([0, 0, 0]);
        assert_eq!(zero.decode(), [0, 0, 0]);
        assert!(zero.is_zero());

        // Test max range (bias 8, 4 bits -> [-8, 7])
        let max = EncodedImage::new([7, 7, 7]);
        assert_eq!(max.decode(), [7, 7, 7]);

        let min = EncodedImage::new([-8, -8, -8]);
        assert_eq!(min.decode(), [-8, -8, -8]);
    }

    #[test]
    fn test_encoded_image_arithmetic() {
        let a = EncodedImage::new([1, 2, 3]);
        let b = EncodedImage::new([-1, 0, 1]);

        // Add
        let sum = a + b;
        assert_eq!(sum.decode(), [0, 2, 4]);

        // Sub
        let sub = a - b;
        assert_eq!(sub.decode(), [2, 2, 2]);

        // Helper image_add
        let helper = a.image_add([-1, -2, -3]);
        assert!(helper.is_zero());
    }

    #[test]
    #[should_panic(expected = "Image out of encoding range")]
    fn test_encoded_image_out_of_bounds() {
        EncodedImage::new([8, 0, 0]); // Max is 7
    }

    // --- EncodedAtom Tests ---

    #[test]
    fn test_encoded_atom_round_trip() {
        let atom_id = 12345;
        let offset = EncodedImage::new([1, -1, 0]);
        let encoded = EncodedAtom::new(atom_id, offset);

        assert_eq!(encoded.atom_index(), atom_id);
        assert_eq!(encoded.image().decode(), offset.decode());

        let (id, img) = encoded.decode_full();
        assert_eq!(id, atom_id);
        assert_eq!(img, [1, -1, 0]);
    }

    #[test]
    fn test_encoded_atom_zero_round_trip() {
        let atom_id = 12345;
        let encoded = EncodedAtom::new_zero_image(atom_id);

        let (id, img) = encoded.decode_full();
        assert_eq!(id, atom_id);
        assert_eq!(img, [0, 0, 0]);
    }

    #[test]
    fn test_encoded_atom_operations() {
        let start = EncodedAtom::new(1, EncodedImage::new([0, 0, 0]));
        let shift = EncodedImage::new([1, 1, 1]);

        let moved = start.image_add(shift);
        assert_eq!(moved.image().decode(), [1, 1, 1]);
        assert_eq!(moved.atom_index(), 1);

        let back = moved.image_sub(shift);
        assert_eq!(back.image().decode(), [0, 0, 0]);
    }

    // --- BlockingVoxelMap & VoxelMap Integration Tests ---

    #[test]
    fn test_voxel_map_full_flow() {
        // 1. SETUP
        // Create a 4x4x4 grid (64 voxels)
        let grid_dims = [4, 4, 4];
        let lattice = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
        let voxel_origin = [0.0, 0.0, 0.0];

        let b_map = BlockingVoxelMap::new(grid_dims, lattice, voxel_origin);

        // 2. POPULATION (Simulate threads)

        // We need to construct weights
        let atom1 = EncodedAtom::new(0, EncodedImage::new([0, 0, 0]));
        let atom2 = EncodedAtom::new(1, EncodedImage::new([0, 0, 0]));

        let w1 = EncodedWeight::new(atom1, 0.5);
        let w2 = EncodedWeight::new(atom2, 0.5);
        let weights = vec![w1, w2].into_boxed_slice();

        // Case A: Voxel 0 is a Maxima for Atom 1
        b_map.maxima_store(0, atom1.0 as isize);

        // Case B: Voxel 1 is a Boundary (50% Atom 1, 50% Atom 2)
        b_map.weight_store(1, weights);

        // Case C: Voxel 2 is Vacuum (Explicitly stored as -1 or just ignored)
        // In BlockingVoxelMap, default is -1 (Vacuum/Unset).
        // We won't touch it, or we can explicitly set it if we had a method for it.

        // 3. CONVERSION
        let v_map = VoxelMap::from_blocking_voxel_map(b_map);

        // 4. ASSERTIONS on VoxelMap

        // Check Voxel 0 (Maxima)
        match v_map.voxel_get(0) {
            Voxel::Maxima(a) => assert_eq!(a.atom_index(), 0),
            _ => panic!("Voxel 0 should be Maxima"),
        }

        // Check Voxel 1 (Boundary)
        match v_map.voxel_get(1) {
            Voxel::Boundary(weights) => {
                assert!(
                    (weights.get(&atom1).unwrap() - 0.5).abs() < f32::EPSILON
                );
                assert!(
                    (weights.get(&atom2).unwrap() - 0.5).abs() < f32::EPSILON
                );
            }
            _ => panic!("Voxel 1 should be Boundary"),
        }

        // Check Voxel 2 (Vacuum)
        // Since we never wrote to it, it should remain -1 (Vacuum)
        match v_map.voxel_get(2) {
            Voxel::Vacuum => (),
            _ => panic!("Voxel 2 should be Vacuum"),
        }
    }

    #[test]
    fn test_volume_map_generation() {
        // Setup a simple map manually to test the `volume_map` logic
        // Grid is irrelevant for this specific test, but required for struct
        let grid = Grid::new(
            [4, 4, 4],
            [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
            [0.0; 3],
        );

        // Construct Weights: Voxel 1 is shared between Atom 1 and Atom 2
        let atom1 = EncodedAtom::new(1, EncodedImage::new([0, 0, 0]));
        let atom2 = EncodedAtom::new(2, EncodedImage::new([0, 0, 0]));
        let w_box = vec![
            EncodedWeight::new(atom1, 0.25),
            EncodedWeight::new(atom2, 0.75),
        ]
        .into_boxed_slice();

        // Map:
        // 0 -> Maxima (Atom 1)
        // 1 -> Boundary (points to w_box)
        // 2 -> Vacuum
        let voxel_data = vec![1, -2, -1];
        let weight_data = vec![w_box];

        let map = VoxelMap::new(voxel_data, weight_data, grid);

        // Get volume map for Atom 1
        let volumes = map.volume_map(1);

        assert_eq!(volumes[0], Some(1.0)); // Fully Atom 1
        assert_eq!(volumes[1], Some(0.25)); // Partially Atom 1
        assert_eq!(volumes[2], None); // Vacuum

        // Test Multi-Volume Map (e.g. Atoms 1 and 2 combined)
        let mut set = FxHashSet::default();
        set.insert(1);
        set.insert(2);

        let multi_vol = map.multi_volume_map(&set);
        assert_eq!(multi_vol[0], Some(1.0)); // Atom 1 is in set
        assert_eq!(multi_vol[1], Some(1.0)); // 0.25 + 0.75 = 1.0
        assert_eq!(multi_vol[2], None);
    }
}
