use std::sync::atomic::AtomicUsize;
use std::sync::Arc;
use std::thread;
use std::time::Duration;

pub struct Bar {
    counter: Arc<AtomicUsize>,
}

impl Bar {
    pub fn new(length: usize, text: String, hidden: bool) -> Self {
        let counter = Arc::new(AtomicUsize::new(0));
        let thread_counter = counter.clone();
        if let false = hidden {
            thread::spawn(move || {
                while Arc::strong_count(&thread_counter) > 1 {
                    let count = thread_counter.fetch_min(length, std::sync::atomic::Ordering::Relaxed);
                    if let std::cmp::Ordering::Less = count.cmp(&length) {
                        let progress = (count * 40) / length;
                        let done = format!("{:=<width$}", "", width = progress);
                        let remain =
                            format!("{:-<width$}", "", width = 39 - progress);
                        eprint!("\r{}: [{}>{}]", text, done, remain);
                    } else {
                        break;
                    }
                    thread::sleep(Duration::from_millis(100));
                }
                eprint!("\r{: <width$}\r", " ", width = text.len() + 44);
            });
        }
        Self { counter }
    }

    pub fn visible(length: usize, text: String) -> Self {
        Self::new(length, text, false)
    }

    pub fn tick(&self) {
        self.counter
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }

    pub fn set_length(&self, length: usize) {
        self.counter
            .store(length, std::sync::atomic::Ordering::Relaxed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn progress_new() {
        let bar = Bar::new(10, String::new(), false);
        assert_eq!(bar.counter.load(std::sync::atomic::Ordering::Relaxed), 0)
    }

    #[test]
    fn progress_tick() {
        let bar = Bar::new(10, String::new(), false);
        bar.tick();
        assert_eq!(bar.counter.load(std::sync::atomic::Ordering::Relaxed), 1)
    }
}
