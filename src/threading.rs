use crate::{
    critical::{CriticalPoint, CriticalPointKey},
    progress::ProgressBar,
};
use rustc_hash::FxHashMap;
use std::thread;

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
        for handle in handle_iter.into_iter() {
            reduce_op(&mut global_state, handle.join().unwrap());
        }
        global_state
    })
}

pub fn parallel_prune<F>(
    critical_points: &[CriticalPoint],
    density: &[f64],
    validator: F,
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
        || FxHashMap::default(),
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
    final_map
        .into_iter()
        .map(|(k, v)| CriticalPoint::new(v.position, v.kind, k.into_box()))
        .collect()
}
