use crate::grid::Grid;
use std::mem::MaybeUninit;
use std::ops::{Add, Sub};
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
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EncodedImage(u16);

impl EncodedImage {
    const BIAS: i8 = 8;
    const BITS: usize = 4;
    const MASK: u16 = (1 << Self::BITS) - 1;
    const ZERO: u16 = (Self::BIAS as u16)
        | ((Self::BIAS as u16) << Self::BITS)
        | ((Self::BIAS as u16) << (2 * Self::BITS));

    /// Create a new encoded image from an [i8; 3]
    ///
    /// # Panics
    /// In debug builds, panics if a component is outside `[-8, 7]`.
    pub fn new(image: [i8; 3]) -> Self {
        let mut encoded: u16 = 0;
        for (i, img) in image.iter().enumerate() {
            let biased = (img + Self::BIAS) as u16;
            debug_assert!(biased <= Self::MASK, "Image out of encoding range");
            encoded |= (biased & Self::MASK) << (i * Self::BITS);
        }
        Self(encoded)
    }

    /// Decode to an [i8; 3]
    pub fn decode(self) -> [i8; 3] {
        [0, 1, 2].map(|i| {
            ((self.0 >> (i * Self::BITS)) & Self::MASK) as i8 - Self::BIAS
        })
    }

    /// Add a raw image to this one.
    pub fn image_add(self, b: [i8; 3]) -> Self {
        let a = self.decode();
        Self::new([a[0] + b[0], a[1] + b[1], a[2] + b[2]])
    }

    /// Whether this is the image `[0, 0, 0]`.
    pub fn is_zero(&self) -> bool {
        self.0 == Self::ZERO
    }
}

impl Add for EncodedImage {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.image_add(rhs.decode())
    }
}

impl Sub for EncodedImage {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let b = rhs.decode();
        self.image_add([-b[0], -b[1], -b[2]])
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
/// assert_eq!(encoded.image(), img);
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
#[repr(transparent)]
pub struct EncodedAtom(u32);

impl EncodedAtom {
    const SHIFT: usize = EncodedImage::BITS * 3;
    const MASK: u32 = (1 << Self::SHIFT) - 1;
    const BITS: u32 = 20;
    const MAX_ATOM: u32 = (1 << Self::BITS) - 1;

    /// Create a new encoded atom from an atom index and its image.
    ///
    /// # Panics
    /// In debug builds, panics if `atom` cannot be stored in 20 bits.
    pub fn new(atom: u32, image: EncodedImage) -> Self {
        debug_assert!(
            atom < Self::MAX_ATOM,
            "Atom Number out of EncodedAtom Range"
        );
        Self((atom << Self::SHIFT) | (image.0) as u32)
    }

    /// Create a new encoded atom in the image `[0, 0, 0]`.
    pub fn new_zero_image(atom: u32) -> Self {
        Self::new(atom, EncodedImage(EncodedImage::ZERO))
    }

    /// The index of the atom.
    pub fn atom_index(&self) -> u32 {
        self.0 >> Self::SHIFT
    }

    /// The image of the atom.
    pub fn image(&self) -> EncodedImage {
        EncodedImage((self.0 & Self::MASK) as u16)
    }

    /// The same atom with `image` added to its image.
    pub fn image_add(self, image: EncodedImage) -> Self {
        Self::new(self.atom_index(), self.image() + image)
    }

    /// The same atom with `image` subtracted from its image.
    pub fn image_sub(self, image: EncodedImage) -> Self {
        Self::new(self.atom_index(), self.image() - image)
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
    /// Pack an atom and its weight.
    pub fn new(encoded_atom: EncodedAtom, weight: f32) -> Self {
        Self((encoded_atom.0 as u64) | ((weight.to_bits() as u64) << 32))
    }

    /// Unpack into the atom and its weight.
    pub fn decode(self) -> (EncodedAtom, f32) {
        (
            EncodedAtom(self.0 as u32),
            f32::from_bits((self.0 >> 32) as u32),
        )
    }
}

/// Describes the state of a voxel.
///
/// # Storage
/// Both voxel maps store an `isize` for every voxel:
/// * **`>= 0`**: The voxel belongs entirely to the [`EncodedAtom`] with this value.
/// * **`-1`**: The voxel is vacuum (or not yet assigned in a [`BlockingVoxelMap`]).
/// * **`< -1`**: The voxel is on a boundary. The value is `-2 - index`, where `index`
///   points to a slice of weights in the weight map.
pub enum Voxel<'a> {
    /// Entirely assigned to the maximum of this atom.
    Maxima(EncodedAtom),
    /// Shared between the atoms in the slice with the given weights.
    Boundary(&'a [EncodedWeight]),
    /// A voxel beneath the vacuum tolerance and not contributing to any maxima.
    Vacuum,
}

impl<'a> Voxel<'a> {
    /// Decode a stored value, `weights` looks up the slice for a boundary index.
    fn from_stored(
        stored: isize,
        weights: impl FnOnce(usize) -> &'a [EncodedWeight],
    ) -> Self {
        match stored {
            -1 => Voxel::Vacuum,
            i if i < -1 => Voxel::Boundary(weights((-2 - i) as usize)),
            i => Voxel::Maxima(EncodedAtom(i as u32)),
        }
    }
}

/// Append-only, lock-free storage for the weights of boundary voxels.
///
/// Slots are claimed with an atomic counter so each is written exactly once,
/// after which it may be read freely. The slot index is published to readers
/// through the voxel map's Release/Acquire ordering so a slot is never read
/// before it has been written.
struct WeightStore {
    slots: Box<[MaybeUninit<Box<[EncodedWeight]>>]>,
    len: AtomicUsize,
}

impl WeightStore {
    fn new(capacity: usize) -> Self {
        let mut slots = Vec::with_capacity(capacity);
        slots.resize_with(capacity, MaybeUninit::uninit);
        Self {
            slots: slots.into_boxed_slice(),
            len: AtomicUsize::new(0),
        }
    }

    /// Store the weights and return their slot index.
    fn push(&self, weights: Box<[EncodedWeight]>) -> usize {
        let i = self.len.fetch_add(1, Ordering::Relaxed);
        let slot = &self.slots[i] as *const _
            as *mut MaybeUninit<Box<[EncodedWeight]>>;
        // SAFETY: the counter hands out each index once so no other thread
        // writes this slot, and readers only access it after the index has
        // been published with Release ordering.
        unsafe { (*slot).write(weights) };
        i
    }

    /// The weights stored in slot `i`, which must already have been written.
    fn get(&self, i: usize) -> &[EncodedWeight] {
        // SAFETY: callers obtain `i` from the voxel map after an Acquire load
        // which synchronises with the Release store that followed the write.
        unsafe { self.slots[i].assume_init_ref() }
    }

    fn into_vec(mut self) -> Vec<Box<[EncodedWeight]>> {
        let len = *self.len.get_mut();
        let slots = std::mem::take(&mut self.slots);
        *self.len.get_mut() = 0;
        slots
            .into_vec()
            .into_iter()
            .take(len)
            // SAFETY: the first `len` slots have been written.
            .map(|slot| unsafe { slot.assume_init() })
            .collect()
    }
}

impl Drop for WeightStore {
    fn drop(&mut self) {
        for slot in self.slots.iter_mut().take(*self.len.get_mut()) {
            // SAFETY: the first `len` slots have been written and are
            // dropped exactly once here.
            unsafe { slot.assume_init_drop() };
        }
    }
}

/// A thread-safe, write-optimised map for populating Bader volumes.
///
/// Designed for concurrent generation of voxel assignments. Threads can safely
/// store "Maxima" (atoms) or "Weights" (boundary contributions) without
/// global locking. Reading a voxel blocks until it has been stored, see
/// [`Voxel`] for how the states are stored.
///
/// # Examples
/// ```
/// use bader::voxel_map::{BlockingVoxelMap, EncodedWeight, EncodedAtom, EncodedImage, Voxel};
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
/// map.maxima_store(0, atom); // Voxel 0 -> Atom 100 Image [0, 0, 0]
///
/// // 3. Store Weights (Boundary)
/// let w = EncodedWeight::new(atom, 0.5);
/// map.weight_store(1, vec![w].into_boxed_slice());
///
/// assert!(matches!(map.voxel_get(0), Voxel::Maxima(a) if a == atom));
/// assert!(matches!(map.voxel_get(1), Voxel::Boundary(w) if w.len() == 1));
/// ```
pub struct BlockingVoxelMap {
    weight_map: WeightStore,
    voxel_map: Box<[AtomicIsize]>,
    /// The Grid used to navigate the VoxelMap.
    pub grid: Grid,
}

impl BlockingVoxelMap {
    /// Initialises a [`BlockingVoxelMap`] and the [`Grid`] that will facilitate movement around the
    /// map.
    pub fn new(
        grid: [usize; 3],
        lattice: [[f64; 3]; 3],
        voxel_origin: [f64; 3],
    ) -> Self {
        let grid = Grid::new(grid, lattice, voxel_origin);
        let size = grid.size.total;
        let mut voxel_map = Vec::with_capacity(size);
        voxel_map.resize_with(size, || AtomicIsize::new(-1));
        Self {
            weight_map: WeightStore::new(size),
            voxel_map: voxel_map.into_boxed_slice(),
            grid,
        }
    }

    /// Retrieves the state of the voxel, `p`, blocking until it has been
    /// stored.
    ///
    /// Calling this on a voxel that is below the vacuum tolerance will block
    /// forever as such voxels are never stored.
    pub fn voxel_get(&self, p: isize) -> Voxel<'_> {
        let stored = loop {
            match self.voxel_map[p as usize].load(Ordering::Acquire) {
                -1 => std::thread::yield_now(),
                x => break x,
            }
        };
        Voxel::from_stored(stored, |i| self.weight_map.get(i))
    }

    /// Assigns the voxel, `p`, entirely to the maximum of `atom`.
    pub fn maxima_store(&self, p: isize, atom: EncodedAtom) {
        self.voxel_map[p as usize].store(atom.0 as isize, Ordering::Release);
    }

    /// Stores the weights of the boundary voxel, `p`.
    pub fn weight_store(&self, p: isize, weights: Box<[EncodedWeight]>) {
        let i = self.weight_map.push(weights);
        self.voxel_map[p as usize].store(-2 - (i as isize), Ordering::Release);
    }

    /// Extract the voxel map data.
    pub fn into_inner(self) -> (Vec<isize>, Vec<Box<[EncodedWeight]>>, Grid) {
        (
            self.voxel_map
                .iter()
                .map(|x| x.load(Ordering::Relaxed))
                .collect(),
            self.weight_map.into_vec(),
            self.grid,
        )
    }
}

/// A read-optimised, non-blocking map for analysing Bader partitions.
///
/// While [`BlockingVoxelMap`] is designed for concurrent *write* operations during the
/// partitioning phase, `VoxelMap` is designed for efficient *read* operations during
/// the analysis phase. See [`Voxel`] for how the states are stored.
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
/// // 3. Analyse volumes
/// let atom_volume = map.volume_map(0); // Get contributions for Atom 0
/// assert_eq!(atom_volume.len(), 8);
/// ```
pub struct VoxelMap {
    /// The stored state of each voxel, see [`Voxel`].
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

    /// Retrieval of the state of the voxel, `p`.
    pub fn voxel_get(&self, p: isize) -> Voxel<'_> {
        Voxel::from_stored(self.voxel_map[p as usize], |i| &self.weight_map[i])
    }

    /// Produce the fraction of each voxel that belongs to the atom with index
    /// `atom`, summed over every image of the atom.
    pub fn volume_map(&self, atom: usize) -> Vec<f64> {
        (0..self.voxel_map.len())
            .map(|p| match self.voxel_get(p as isize) {
                Voxel::Maxima(a) => {
                    (a.atom_index() as usize == atom) as u8 as f64
                }
                Voxel::Boundary(weights) => weights
                    .iter()
                    .map(|w| w.decode())
                    .filter(|(a, _)| a.atom_index() as usize == atom)
                    .map(|(_, w)| w as f64)
                    .sum(),
                Voxel::Vacuum => 0.0,
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert_eq!(encoded.image(), offset);
        assert_eq!(encoded.image().decode(), [1, -1, 0]);
    }

    #[test]
    fn test_encoded_atom_zero_round_trip() {
        let atom_id = 12345;
        let encoded = EncodedAtom::new_zero_image(atom_id);
        assert_eq!(encoded.atom_index(), atom_id);
        assert!(encoded.image().is_zero());
        assert_eq!(
            encoded,
            EncodedAtom::new(atom_id, EncodedImage::new([0; 3]))
        );
    }

    #[test]
    fn test_encoded_weight_round_trip() {
        let atom = EncodedAtom::new(7, EncodedImage::new([1, 0, -1]));
        let (a, w) = EncodedWeight::new(atom, 0.25).decode();
        assert_eq!(a, atom);
        assert_eq!(w, 0.25);
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
        b_map.maxima_store(0, atom1);

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
                assert_eq!(weights.len(), 2);
                assert_eq!(weights[0].decode(), (atom1, 0.5));
                assert_eq!(weights[1].decode(), (atom2, 0.5));
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
        let voxel_data = vec![atom1.0 as isize, -2, -1];
        let weight_data = vec![w_box];

        let map = VoxelMap::new(voxel_data, weight_data, grid);

        // Get volume map for Atom 1
        let volumes = map.volume_map(1);
        assert_eq!(volumes, vec![1.0, 0.25, 0.0]);
        // And for Atom 2
        let volumes = map.volume_map(2);
        assert_eq!(volumes, vec![0.0, 0.75, 0.0]);
    }

    #[test]
    fn test_volume_map_sums_images() {
        let grid = Grid::new(
            [4, 4, 4],
            [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
            [0.0; 3],
        );
        // an atom bonded to its own image shares the voxel with itself
        let w_box = vec![
            EncodedWeight::new(EncodedAtom::new_zero_image(1), 0.25),
            EncodedWeight::new(
                EncodedAtom::new(1, EncodedImage::new([1, 0, 0])),
                0.75,
            ),
        ]
        .into_boxed_slice();
        let map = VoxelMap::new(vec![-2], vec![w_box], grid);
        assert_eq!(map.volume_map(1), vec![1.0]);
    }

    #[test]
    fn test_blocking_voxel_map_drop_without_into_inner() {
        // weights stored in a BlockingVoxelMap that is dropped must not leak
        // or double free
        let b_map = BlockingVoxelMap::new(
            [2, 2, 2],
            [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
            [0.0; 3],
        );
        let atom = EncodedAtom::new_zero_image(0);
        b_map.weight_store(0, vec![EncodedWeight::new(atom, 1.0)].into());
        b_map.weight_store(1, vec![EncodedWeight::new(atom, 1.0)].into());
        drop(b_map);
    }
}
