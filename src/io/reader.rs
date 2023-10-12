use std::fs::File;
use std::io::{BufRead, BufReader as Reader, Result};

/// Read a file into a mutable buffer
pub struct BufReader {
    reader: Reader<File>,
}

impl BufReader {
    /// Opens the file from the path into a reader
    pub fn open(path: impl AsRef<std::path::Path>) -> Result<Self> {
        let file = File::open(path)?;
        let reader = Reader::new(file);

        Ok(Self { reader })
    }

    /// Reads a line from the buffer reader to mutable string
    pub fn read_line<'buf>(
        &mut self,
        buffer: &'buf mut String,
    ) -> Option<Result<(&'buf mut String, usize)>> {
        buffer.clear();

        self.reader
            .read_line(buffer)
            .map(|u| if u == 0 { None } else { Some((buffer, u)) })
            .transpose()
    }
}
