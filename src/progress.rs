use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::thread;
use std::time::Duration;

/// Trait for reporting progress.
///
/// This allows for different implementations of progress reporting, such as
/// a visible progress bar for interactive use or a silent no-op for scripts/logging.
pub trait ProgressBar: Send + Sync {
    /// Increment the progress counter by one.
    fn tick(&self);
}

/// A threaded, visible progress bar.
///
/// This struct spawns a background thread upon creation that periodically updates
/// the progress bar printed to `stderr`. It uses an `AtomicUsize` to share state
/// between the main worker threads (which call `tick`) and the display thread.
pub struct Bar {
    counter: Arc<AtomicUsize>,
}

impl Bar {
    /// Creates a new `Bar` and starts the display thread.
    ///
    /// # Arguments
    /// * `length` - The total number of items to process (100% progress).
    /// * `text` - A label to display alongside the bar.
    ///
    /// The display thread will exit automatically when the `Bar` struct is dropped
    /// (when `Arc::strong_count` drops to 1) or when the counter reaches `length`.
    pub fn new(length: usize, text: String) -> Self {
        let counter = Arc::new(AtomicUsize::new(0));
        let thread_counter = counter.clone();
        thread::spawn(move || {
            while Arc::strong_count(&thread_counter) > 1 {
                let count = thread_counter
                    .fetch_min(length, std::sync::atomic::Ordering::Relaxed);
                if let std::cmp::Ordering::Less = count.cmp(&length) {
                    let progress = (count * 40) / length;
                    let done = format!("{:=<width$}", "", width = progress);
                    let remain =
                        format!("{:-<width$}", "", width = 39 - progress);
                    eprint!("\r{}: [{}>{}]", text, done, remain);
                } else {
                    break;
                }
                thread::sleep(Duration::from_millis(20));
            }
            eprint!("\r{: <width$}\r", " ", width = text.len() + 44);
        });
        Self { counter }
    }
}

impl ProgressBar for Bar {
    fn tick(&self) {
        self.counter
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }
}

/// A silent progress bar (No-op).
///
/// Used when the user requests "silent" mode or when output is redirected.
/// All calls to `tick` are ignored.
pub struct HiddenBar {}
impl ProgressBar for HiddenBar {
    fn tick(&self) {}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn progress_new() {
        let bar = Bar::new(10, String::new());
        assert_eq!(bar.counter.load(std::sync::atomic::Ordering::Relaxed), 0)
    }

    #[test]
    fn progress_tick() {
        let bar = Bar::new(10, String::new());
        bar.tick();
        assert_eq!(bar.counter.load(std::sync::atomic::Ordering::Relaxed), 1)
    }

    #[test]
    fn test_bar_concurrent_ticks() {
        // Verify that the atomic counter handles concurrent updates correctly
        let bar = Arc::new(Bar::new(1000, "Concurrent".to_string()));
        let mut handles = vec![];

        for _ in 0..10 {
            let b = bar.clone();
            handles.push(thread::spawn(move || {
                for _ in 0..100 {
                    b.tick();
                }
            }));
        }

        for h in handles {
            h.join().unwrap();
        }

        assert_eq!(
            bar.counter.load(std::sync::atomic::Ordering::Relaxed),
            1000
        );
    }

    #[test]
    fn test_hidden_bar() {
        // Verify HiddenBar implements the trait and doesn't panic
        let bar = HiddenBar {};
        bar.tick();
        bar.tick();
    }
}
