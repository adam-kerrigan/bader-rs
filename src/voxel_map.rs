use crate::grid::Grid;
use rustc_hash::{FxHashMap, FxHashSet};
use std::cell::UnsafeCell;
use std::mem::MaybeUninit;
use std::ops::{Deref, DerefMut};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicIsize, AtomicUsize, Ordering};

pub struct EncodedData(u64);

impl EncodedData {
    pub fn new(encoded_maxima: u32, weight: f32) -> Self {
        Self(u64::from_le_bytes(
            [encoded_maxima.to_le_bytes(), weight.to_le_bytes()]
                .concat()
                .try_into()
                .unwrap(),
        ))
    }

    pub fn decode_maxima(encoded_maxima: u32) -> (u16, u16) {
        let bytes = encoded_maxima.to_le_bytes();
        (
            u16::from_le_bytes(bytes[0..2].try_into().unwrap()),
            u16::from_le_bytes(bytes[2..4].try_into().unwrap()),
        )
    }
    pub fn encode_maxima(maxima: u16, encoded_image: u16) -> u32 {
        u32::from_le_bytes(
            [maxima.to_le_bytes(), encoded_image.to_le_bytes()]
                .concat()
                .try_into()
                .unwrap(),
        )
    }

    pub fn add_image(encoded_maxima: u32, other_image: [i8; 3]) -> u32 {
        let (maxima, image) = Self::decode_maxima(encoded_maxima);
        let decoded_image = Self::decode_image(image);
        let encoded_image = Self::encode_image(
            decoded_image
                .into_iter()
                .zip(other_image)
                .map(|(i, ii)| i + ii)
                .collect::<Vec<i8>>()
                .try_into()
                .unwrap(),
        );
        Self::encode_maxima(maxima, encoded_image)
    }
    pub fn subtract_image(encoded_maxima: u32, other_image: [i8; 3]) -> u32 {
        let (maxima, image) = Self::decode_maxima(encoded_maxima);
        let decoded_image = Self::decode_image(image);
        let encoded_image = Self::encode_image(
            decoded_image
                .into_iter()
                .zip(other_image)
                .map(|(i, ii)| i - ii)
                .collect::<Vec<i8>>()
                .try_into()
                .unwrap(),
        );
        Self::encode_maxima(maxima, encoded_image)
    }

    pub fn add_encoded_image(encoded_maxima: u32, other_image: u16) -> u32 {
        let (maxima, image) = Self::decode_maxima(encoded_maxima);
        let decoded_image = Self::decode_image(image);
        let decoded_other_image = Self::decode_image(other_image);
        let encoded_image = Self::encode_image(
            decoded_image
                .into_iter()
                .zip(decoded_other_image)
                .map(|(i, ii)| i + ii)
                .collect::<Vec<i8>>()
                .try_into()
                .unwrap(),
        );
        Self::encode_maxima(maxima, encoded_image)
    }

    pub fn decode_image(encoded_image: u16) -> [i8; 3] {
        (0..15)
            .step_by(5)
            .map(|shift: u16| {
                (((encoded_image >> shift) & 0x1F) as i8) << 3 >> 3
            })
            .collect::<Vec<i8>>()
            .try_into()
            .unwrap()
    }

    pub fn encode_image(image: [i8; 3]) -> u16 {
        image
            .into_iter()
            .enumerate()
            .fold(0, |acc, (i, img)| acc | (((img as u16) & 0x1F) << i * 5))
    }

    pub fn decode_self(&self) -> (u32, f32) {
        let bytes: [u8; 8] = self.0.to_le_bytes();
        let maxima = u32::from_le_bytes(bytes[0..4].try_into().unwrap());
        let weight = f32::from_le_bytes(bytes[4..8].try_into().unwrap());
        (maxima, weight)
    }
}

/// Describes the state of the voxel.
pub enum Voxel {
    /// Contians the position of the voxel's maxima.
    Maxima(usize),
    /// Contians a vector of the maxima the current voxel contributes to and
    /// their weights.
    Boundary(FxHashMap<u32, f32>),
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
    weight_map: Arc<[MaybeUninit<Box<[EncodedData]>>]>,
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
        weight_map.resize_with(size, || MaybeUninit::uninit());
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
    pub fn weight_get(&self, i: isize) -> FxHashMap<u32, f32> {
        let i = -2 - i;
        (unsafe { self.weight_map.get_unchecked(i as usize).assume_init_ref() })
            .iter()
            .map(|u| u.decode_self())
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
            std::cmp::Ordering::Greater => Voxel::Maxima(i as usize),
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
    pub fn weight_store(&self, p: isize, weights: Box<[EncodedData]>) {
        let i = self.weight_counter.fetch_add(1, Ordering::Relaxed);
        unsafe {
            let ptr: *mut Box<[EncodedData]> =
                self.weight_map.get_unchecked(i) as *const _ as *mut _;
            ptr.write(weights)
        }
        self.maxima_store(p, -2 - (i as isize));
    }

    /// Extract the voxel map data.
    pub fn into_inner(self) -> (Vec<isize>, Vec<Box<[EncodedData]>>, Grid) {
        (
            self.voxel_map
                .into_iter()
                .map(|x| x.load(Ordering::Relaxed))
                .collect(),
            self.weight_map
                .into_iter()
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
    pub weight_map: Vec<Box<[EncodedData]>>,
    /// The Grid used to navigate the VoxelMap.
    pub grid: Grid,
}

impl VoxelMap {
    /// Create a new [`VoxelMap`]
    pub fn new(
        voxel_map: Vec<isize>,
        weight_map: Vec<Box<[EncodedData]>>,
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
    pub fn weight_iter(&self) -> std::slice::Iter<'_, Box<[EncodedData]>> {
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
            std::cmp::Ordering::Greater => Voxel::Maxima(maxima as usize),
            std::cmp::Ordering::Less => {
                Voxel::Boundary(self.maxima_to_weight(maxima))
            }
        }
    }

    /// Return a reference to the weights from the given maxima, Note: maxima here must be < -1.
    pub fn maxima_to_weight(&self, maxima: isize) -> FxHashMap<u32, f32> {
        self.weight_map[(-2 - maxima) as usize]
            .iter()
            .map(|ed| ed.decode_self())
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
                        if (m as isize) == volume_number {
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
                        if volume_numbers.contains(&(m as isize)) {
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
