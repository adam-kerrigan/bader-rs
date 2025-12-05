use crate::{
    critical::{CriticalPoint, CriticalPointKey},
    hash::{IntSet, SliceMap},
    progress::ProgressBar,
};
use std::thread;

/// A generic parallel Map-Reduce implementation using `std::thread::scope`.
///
/// This function splits the input `items` into chunks, processes each chunk in a separate thread
/// to build a local state (Map phase), and then combines these local states into a final result (Reduce phase).
///
/// # Type Parameters
/// * `State`: The type of the accumulator
/// * `Item`: The type of the input data items.
/// * `Init`: Function to initialise the thread-local state.
/// * `Map`: Function to process a single item and update the local state.
/// * `Reduce`: Function to merge a thread-local state into the global state.
///
/// # Arguments
/// * `items`: The slice of data to process.
/// * `init_state`: Generator for the initial state of each thread.
/// * `map_op`: The mapper function `fn(&mut State, &Item)`.
/// * `reduce_op`: The reducer function `fn(&mut State, State)`.
/// * `threads`: Number of parallel threads to spawn.
/// * `progress_bar`: [`ProgressBar`] to update as items are processed.
///
/// # Example
/// ```
/// use bader::threading::parallel_map_reduce;
/// use bader::progress::HiddenBar;
///
/// let numbers = vec![1, 2, 3, 4, 5];
/// let sum = parallel_map_reduce(
///     &numbers,
///     || 0,                    // Init
///     |acc, x| *acc += x,      // Map
///     |global, local| *global += local, // Reduce
///     2,                       // Threads
///     Box::new(HiddenBar {})
/// );
/// assert_eq!(sum, 15);
/// ```
pub fn parallel_map_reduce<State, Item, Init, Map, Reduce>(
    items: &[Item],
    init_state: Init,
    map_op: Map,
    reduce_op: Reduce,
    threads: usize,
    progress_bar: Box<dyn ProgressBar>,
) -> State
where
    State: Send,
    Init: Fn() -> State + Sync + Send,
    Map: Fn(&mut State, &Item) + Sync + Send + Copy,
    Reduce: Fn(&mut State, State) + Sync + Send + Copy,
    Item: Sync + Clone,
{
    let total_length = items.len();
    let chunk_size = (total_length / threads) + (total_length % threads).min(1);
    thread::scope(|s| {
        let mut handles = Vec::with_capacity(threads);
        items.chunks(chunk_size).for_each(|chunk| {
            let mut local_state = init_state();
            handles.push(s.spawn(|| {
                chunk.iter().for_each(|index| {
                    map_op(&mut local_state, index);
                    progress_bar.tick();
                });
                local_state
            }));
        });
        let mut handle_iter = handles.into_iter();
        let mut global_state = match handle_iter.next() {
            Some(h) => h.join().unwrap(),
            None => panic!(""),
        };
        for handle in handle_iter {
            reduce_op(&mut global_state, handle.join().unwrap());
        }
        global_state
    })
}

/// Filters and deduplicates Critical Points in parallel.
///
/// This function serves two purposes:
/// 1. **Validation**: Applies a user-defined filter (e.g., "is this a valid ring?") to all points.
/// 2. **Deduplication**: Groups valid points by their unique `CriticalPointKey` (Atom IDs).
///    If multiple points exist for the same key (e.g., multiple candidates for the same bond),
///    only the one with the highest charge density is kept.
///
/// # Logic
/// * **Map Phase**: Each thread processes a chunk of `critical_points`. Valid points are inserted
///   into a local `HashMap`. If a collision occurs (same key), the point with higher density overwrites.
/// * **Reduce Phase**: Local HashMaps are merged into a global HashMap, again preserving only the
///   highest density candidates.
///
/// # Arguments
/// * `critical_points`: List of candidate points.
/// * `density`: Charge density array for comparison.
/// * `validator`: Function returning `true` if a point should be kept.
/// * `threads`: Number of parallel threads to spawn.
/// * `progress_bar`: [`ProgressBar`] to update as items are processed.
///
/// # Returns
/// A vector of unique, validated [`CriticalPoint`]s.
pub fn parallel_prune<F>(
    critical_points: &[CriticalPoint],
    density: &[f64],
    validator: F,
    self_bonds: bool,
    threads: usize,
    progress_bar: Box<dyn ProgressBar>,
) -> Vec<CriticalPoint>
where
    F: Fn(&CriticalPoint) -> bool + Sync + Send + Copy,
{
    let final_map = parallel_map_reduce(
        // Items
        critical_points,
        // Init
        SliceMap::default,
        // Map
        |local_state, cp| {
            if validator(cp) {
                let key = CriticalPointKey::from_cp(cp.clone());
                local_state
                    .entry(key)
                    .and_modify(|existing: &mut CriticalPoint| {
                        let rho_new = density[cp.position as usize];
                        let rho_old = density[existing.position as usize];
                        if rho_new > rho_old {
                            *existing = cp.clone();
                        }
                    })
                    .or_insert(cp.clone());
            }
        },
        // Reduce
        |global_state, local_state| {
            local_state.into_iter().for_each(|(key, cp)| {
                global_state
                    .entry(key)
                    .and_modify(|existing: &mut CriticalPoint| {
                        let rho_new = density[cp.position as usize];
                        let rho_old = density[existing.position as usize];
                        if rho_new > rho_old {
                            *existing = cp.clone();
                        }
                    })
                    .or_insert(cp.clone());
            });
        },
        threads,
        progress_bar,
    );
    let mut critical_points: Vec<CriticalPoint> = final_map
        .into_iter()
        .map(|(k, v)| CriticalPoint::new(v.position, v.kind, k.into_box()))
        .collect();
    if !self_bonds {
        critical_points.retain(|cp| {
            let atom_set = IntSet::from_iter(
                cp.atoms
                    .iter()
                    .map(|encoded_atom| encoded_atom.atom_index()),
            );
            atom_set.len() == cp.atoms.len()
        });
    }
    critical_points
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::critical::{CriticalPoint, CriticalPointKind};
    use crate::progress::HiddenBar;
    use crate::voxel_map::EncodedAtom;

    // --- Helper to create a basic CriticalPoint ---
    fn create_cp(pos: isize, atoms: &[u32]) -> CriticalPoint {
        let atoms_enc = atoms
            .iter()
            .map(|&id| EncodedAtom::new_zero_image(id))
            .collect::<Vec<_>>()
            .into_boxed_slice();
        CriticalPoint::new(pos, CriticalPointKind::Bond, atoms_enc)
    }

    // --- Parallel Map Reduce Tests ---

    #[test]
    fn test_parallel_map_reduce_sum() {
        // Sum numbers from 1 to 100
        let data: Vec<i32> = (1..=100).collect();
        let threads = 4;
        let pbar = Box::new(HiddenBar {});

        let sum = parallel_map_reduce(
            &data,
            || 0,                             // Init local state
            |acc, item| *acc += item, // Map: add item to local accumulator
            |global, local| *global += local, // Reduce: add local sum to global
            threads,
            pbar,
        );

        assert_eq!(sum, 5050);
    }

    #[test]
    fn test_parallel_map_reduce_chunking() {
        // Ensure it works even if threads > items
        let data = vec![1, 2, 3];
        let threads = 8;
        let pbar = Box::new(HiddenBar {});

        let sum = parallel_map_reduce(
            &data,
            || 0,
            |acc, item| *acc += item,
            |global, local| *global += local,
            threads,
            pbar,
        );

        assert_eq!(sum, 6);
    }

    // --- Parallel Prune Tests ---

    #[test]
    fn test_parallel_prune_logic() {
        // Setup: Two Critical Points identifying the SAME bond (same atoms [1, 2])
        // CP1 at pos 10, Density = 1.0
        // CP2 at pos 20, Density = 5.0
        // Expected: Pruning should keep only CP2 (highest density).

        let cp1 = create_cp(10, &[1, 2]);
        let cp2 = create_cp(20, &[1, 2]);

        let cps = vec![cp1, cp2];
        let mut density = vec![0.0; 30];
        density[10] = 1.0;
        density[20] = 5.0;

        let threads = 2;
        let pbar = Box::new(HiddenBar {});

        // Validator always returns true (keep everything initially)
        let pruned =
            parallel_prune(&cps, &density, |_| true, false, threads, pbar);

        assert_eq!(pruned.len(), 1);
        assert_eq!(pruned[0].position, 20); // Should be the higher density one
    }

    #[test]
    fn test_parallel_prune_filtering() {
        // Setup: CP1 (Bond) and CP2 (Ring). Validator will reject Ring.
        let cp1 = create_cp(10, &[1, 2]); // Bond
        let mut cp2 = create_cp(20, &[3, 4]);
        cp2.kind = CriticalPointKind::Ring;

        let cps = vec![cp1.clone(), cp2];
        let density = vec![0.0; 30]; // Density irrelevant here as they are distinct keys

        let threads = 1;
        let pbar = Box::new(HiddenBar {});

        let pruned = parallel_prune(
            &cps,
            &density,
            |cp| matches!(cp.kind, CriticalPointKind::Bond), // Keep only Bonds
            false,
            threads,
            pbar,
        );

        assert_eq!(pruned.len(), 1);
        assert_eq!(pruned[0].kind, CriticalPointKind::Bond);
    }
}
