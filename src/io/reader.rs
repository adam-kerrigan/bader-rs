use std::fs::File;
use std::io::{self, prelude::*};

/// Read a file into a mutable buffer
pub struct BufReader {
    reader: io::BufReader<File>,
}

impl BufReader {
    /// Opens the file from the path into a reader
    pub fn open(path: impl AsRef<std::path::Path>) -> io::Result<Self> {
        let file = File::open(path)?;
        let reader = io::BufReader::new(file);

        Ok(Self { reader })
    }

    /// Reads a line from the buffer reader to mutable string
    pub fn read_line<'buf>(&mut self,
                           buffer: &'buf mut String)
                           -> Option<io::Result<(&'buf mut String, usize)>>
    {
        buffer.clear();

        self.reader
            .read_line(buffer)
            .map(|u| if u == 0 { None } else { Some((buffer, u)) })
            .transpose()
    }
}
