use atomic_counter::{AtomicCounter, RelaxedCounter};
use indicatif::{ProgressBar, ProgressStyle};
use std::sync::Arc;
use std::thread;
use std::time::Duration;

/// Contains the indicatif progress bar and an atomic_counter RelaxedCounter
/// The counter is used to update the bar for use with a rayon iterator
pub struct Bar {
    counter: Arc<RelaxedCounter>,
    pub pbar: Arc<ProgressBar>,
}

impl Bar {
    /// Creates the Bar struct with a size, refresh_rate and prefix for the bar
    pub fn new(progress_bar: ProgressBar,
               refresh_rate: u64,
               prefix: String)
               -> Self {
        progress_bar.set_prefix(&prefix);
        progress_bar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{prefix}[{bar:40}] [{elapsed_precise}] {percent:>3}%",
                )
                .progress_chars("=>-"),
        );
        let pb = Arc::new(progress_bar);
        let pb2 = pb.clone();
        let counter = Arc::new(RelaxedCounter::new(0));
        let counter2 = counter.clone();
        thread::spawn(move || {
            while Arc::strong_count(&counter2) > 1 && !pb2.is_finished() {
                pb2.set_position(counter2.get() as u64);
                thread::sleep(Duration::from_millis(refresh_rate));
            }
        });
        Self { counter, pbar: pb }
    }

    /// tick the progress bar
    pub fn tick(&self) {
        self.counter.inc();
    }
}

impl Drop for Bar {
    /// make sure we clear bars when the object is dropped
    fn drop(&mut self) {
        if !self.pbar.is_finished() {
            let value = self.counter.get() as u64;
            self.pbar.set_position(value);
            self.pbar.finish_and_clear();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn progress_new() {
        let bar = ProgressBar::hidden();
        let bar = Bar::new(bar, 1, String::new());
        assert_eq!(bar.counter.get(), 0);
    }

    #[test]
    fn progress_tick() {
        let bar = ProgressBar::hidden();
        let bar = Bar::new(bar, 1, String::new());
        bar.tick();
        assert_eq!(bar.counter.get(), 1)
    }
}
