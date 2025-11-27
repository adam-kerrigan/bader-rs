use std::hash::{BuildHasherDefault, Hasher};
use std::ops::BitXor;

const GOLDEN_RATIO: u64 = 0x9E3779B97F4A7C15;

/// A specialized, non-mixing hasher for single integers.
///
/// It ignores previous state and simply multiplies the input by the Golden Ratio.
/// This is the fastest possible hasher for single integer keys.
///
/// # Safety
/// DO NOT use this for composite keys (structs with >1 field or vectors).
/// It will overwrite the hash of previous fields, leading to collisions.
#[derive(Default)]
pub struct SingleIntHasher(u64);

impl Hasher for SingleIntHasher {
    #[inline(always)]
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, _: &[u8]) {
        panic!(
            "SingleIntHasher does not support generic byte streams. Use IntSliceHasher instead."
        );
    }

    #[inline(always)]
    fn write_u32(&mut self, i: u32) {
        self.0 = (i as u64).wrapping_mul(GOLDEN_RATIO);
    }

    #[inline(always)]
    fn write_u64(&mut self, i: u64) {
        self.0 = i.wrapping_mul(GOLDEN_RATIO);
    }

    #[inline(always)]
    fn write_usize(&mut self, i: usize) {
        self.0 = (i as u64).wrapping_mul(GOLDEN_RATIO);
    }
}

/// A specialized hasher for sequences of integers.
///
/// It uses a Rotate-XOR-Multiply (FxHash-like) algorithm to mix multiple integers
/// into a single hash state. Use this for `Vec<u32>`, slices, or composite structs.
#[derive(Default)]
pub struct IntSliceHasher(u64);

impl Hasher for IntSliceHasher {
    #[inline(always)]
    fn finish(&self) -> u64 {
        self.0
    }

    #[inline(always)]
    fn write(&mut self, bytes: &[u8]) {
        // Fallback mixing for generic bytes
        for &b in bytes {
            self.0 = self
                .0
                .rotate_left(5)
                .bitxor(b as u64)
                .wrapping_mul(GOLDEN_RATIO);
        }
    }

    #[inline(always)]
    fn write_u32(&mut self, i: u32) {
        self.mix(i as u64);
    }

    #[inline(always)]
    fn write_u64(&mut self, i: u64) {
        self.mix(i);
    }

    #[inline(always)]
    fn write_usize(&mut self, i: usize) {
        self.mix(i as u64);
    }
}

impl IntSliceHasher {
    #[inline(always)]
    fn mix(&mut self, i: u64) {
        // FxHash-style mixing with Golden Ratio
        self.0 = self.0.rotate_left(5).bitxor(i).wrapping_mul(GOLDEN_RATIO);
    }
}

/// Map for single integer keys (`u32`, `usize`, [`crate::voxel_map::EncodedAtom`]).
/// Fast "Overwrite" hashing.
pub type IntMap<K, V> =
    std::collections::HashMap<K, V, BuildHasherDefault<SingleIntHasher>>;

/// Set for single integer keys.
pub type IntSet<K> =
    std::collections::HashSet<K, BuildHasherDefault<SingleIntHasher>>;

/// Map for composite keys (`Vec<u32>`, [`crate::critical::CriticalPointKey`]).
/// "Mixing" hashing.
pub type SliceMap<K, V> =
    std::collections::HashMap<K, V, BuildHasherDefault<IntSliceHasher>>;

/// Set for composite keys.
pub type SliceSet<K> =
    std::collections::HashSet<K, BuildHasherDefault<IntSliceHasher>>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{HashMap, HashSet};
    use std::hash::{Hash, Hasher};

    // --- SingleIntHasher Tests ---

    #[test]
    fn test_single_int_golden_ratio() {
        let mut hasher = SingleIntHasher::default();
        hasher.write_u32(1);
        let result = hasher.finish();

        // Expect exact multiplication by Golden Ratio
        assert_eq!(result, 0x9E3779B97F4A7C15);
    }

    #[test]
    fn test_single_int_overwrite_behavior() {
        // This hasher is designed to OVERWRITE state, not mix it.
        // This test confirms that writing A then B is exactly the same as just writing B.
        let mut h1 = SingleIntHasher::default();
        h1.write_u32(100); // This should be discarded
        h1.write_u32(200);

        let mut h2 = SingleIntHasher::default();
        h2.write_u32(200); // Just write B

        assert_eq!(
            h1.finish(),
            h2.finish(),
            "SingleIntHasher should overwrite previous state"
        );
    }

    #[test]
    #[should_panic(
        expected = "SingleIntHasher does not support generic byte streams"
    )]
    fn test_single_int_panic_on_bytes() {
        let mut hasher = SingleIntHasher::default();
        hasher.write(&[1, 2, 3]);
    }

    #[test]
    #[should_panic(
        expected = "SingleIntHasher does not support generic byte streams"
    )]
    fn test_single_int_panic_on_string() {
        // Strings use .write(), which should panic
        let mut hasher = SingleIntHasher::default();
        "Hello".hash(&mut hasher);
    }

    #[test]
    fn test_single_int_map_integration() {
        // Verify it works as a BuildHasher for a standard HashMap
        let mut map: IntMap<u32, &str> = HashMap::default();
        map.insert(42, "Answer");
        assert_eq!(map.get(&42), Some(&"Answer"));
    }

    // --- IntSliceHasher Tests ---

    #[test]
    fn test_slice_mixing_order_dependence() {
        // For a slice/vector, [1, 2] must have a different hash than [2, 1].
        // If it used simple XOR without rotation, they would be identical (bad).

        let mut h1 = IntSliceHasher::default();
        h1.write_u32(1);
        h1.write_u32(2);

        let mut h2 = IntSliceHasher::default();
        h2.write_u32(2);
        h2.write_u32(1);

        assert_ne!(
            h1.finish(),
            h2.finish(),
            "IntSliceHasher must be order-dependent"
        );
    }

    #[test]
    fn test_slice_mixing_accumulation() {
        // Verify that writing A then B results in a mixed state, not just B.
        let mut h1 = IntSliceHasher::default();
        h1.write_u32(100);
        let hash_after_first = h1.finish();

        h1.write_u32(200);
        let hash_after_second = h1.finish();

        assert_ne!(hash_after_first, hash_after_second);

        // Also verify it doesn't just equal the hash of 200 alone
        let mut h2 = IntSliceHasher::default();
        h2.write_u32(200);
        assert_ne!(h1.finish(), h2.finish());
    }

    #[test]
    fn test_slice_bytes_support() {
        // Unlike SingleIntHasher, this one SHOULD support bytes (fallback mode)
        let mut h = IntSliceHasher::default();
        h.write(&[1, 2, 3]);
        // Should not panic
        let res = h.finish();
        assert_ne!(res, 0);
    }

    #[test]
    fn test_slice_set_integration() {
        // Verify it works with HashSet for a vector of integers
        let mut set: SliceSet<Vec<u32>> = HashSet::default();

        set.insert(vec![1, 2, 3]);
        set.insert(vec![4, 5, 6]);

        assert!(set.contains(&vec![1, 2, 3]));
        assert!(!set.contains(&vec![3, 2, 1])); // Order matters
    }
}
