use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastqRecord {
    pub id: Vec<u8>,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

pub trait PairSource {
    fn for_each_pair_dyn(
        &self,
        max_pairs: Option<usize>,
        callback: &mut dyn FnMut(usize, FastqRecord, Option<FastqRecord>) -> Result<()>,
    ) -> Result<usize>;
}

pub struct FastqPairSource<'a> {
    pub read1: &'a Path,
    pub read2: Option<&'a Path>,
}

impl PairSource for FastqPairSource<'_> {
    fn for_each_pair_dyn(
        &self,
        max_pairs: Option<usize>,
        callback: &mut dyn FnMut(usize, FastqRecord, Option<FastqRecord>) -> Result<()>,
    ) -> Result<usize> {
        for_each_pair(self.read1, self.read2, max_pairs, |index, left, right| {
            callback(index, left, right)
        })
    }
}

pub struct FastqReader {
    reader: Box<dyn BufRead + Send>,
    line_number: usize,
    path: PathBuf,
}

impl FastqReader {
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path)
            .with_context(|| format!("failed to open FASTQ: {}", path.display()))?;
        let reader: Box<dyn Read + Send> = if path
            .extension()
            .and_then(|value| value.to_str())
            .is_some_and(|value| value.eq_ignore_ascii_case("gz"))
        {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        Ok(Self {
            reader: Box::new(BufReader::with_capacity(1 << 20, reader)),
            line_number: 0,
            path,
        })
    }

    pub fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        let Some(header) = self.read_line()? else {
            return Ok(None);
        };
        let sequence = self
            .read_line()?
            .ok_or_else(|| anyhow::anyhow!("truncated FASTQ sequence line"))?;
        let plus = self
            .read_line()?
            .ok_or_else(|| anyhow::anyhow!("truncated FASTQ plus line"))?;
        let quality = self
            .read_line()?
            .ok_or_else(|| anyhow::anyhow!("truncated FASTQ quality line"))?;

        if !header.starts_with(b"@") {
            bail!(
                "{}:{}: FASTQ header must start with @",
                self.path.display(),
                self.line_number.saturating_sub(3)
            );
        }
        if !plus.starts_with(b"+") {
            bail!(
                "{}:{}: FASTQ separator must start with +",
                self.path.display(),
                self.line_number.saturating_sub(1)
            );
        }
        if sequence.len() != quality.len() {
            bail!(
                "{}:{}: sequence/quality length mismatch: {} != {}",
                self.path.display(),
                self.line_number,
                sequence.len(),
                quality.len()
            );
        }

        Ok(Some(FastqRecord {
            id: header[1..].to_vec(),
            sequence,
            quality,
        }))
    }

    fn read_line(&mut self) -> Result<Option<Vec<u8>>> {
        let mut buffer = Vec::new();
        let bytes = self
            .reader
            .read_until(b'\n', &mut buffer)
            .with_context(|| format!("failed reading {}", self.path.display()))?;
        if bytes == 0 {
            return Ok(None);
        }
        self.line_number += 1;
        while buffer
            .last()
            .is_some_and(|value| *value == b'\n' || *value == b'\r')
        {
            buffer.pop();
        }
        Ok(Some(buffer))
    }
}

pub fn for_each_pair<F>(
    read1: &Path,
    read2: Option<&Path>,
    max_pairs: Option<usize>,
    mut callback: F,
) -> Result<usize>
where
    F: FnMut(usize, FastqRecord, Option<FastqRecord>) -> Result<()>,
{
    let mut reader1 = FastqReader::from_path(read1)?;
    let mut reader2 = read2.map(FastqReader::from_path).transpose()?;
    let mut index = 0_usize;

    loop {
        if max_pairs.is_some_and(|limit| index >= limit) {
            break;
        }
        let left = reader1.next_record()?;
        let right = match reader2.as_mut() {
            Some(reader) => reader.next_record()?,
            None => None,
        };

        match (left, right, reader2.is_some()) {
            (None, None, _) => break,
            (Some(left), Some(right), true) => callback(index, left, Some(right))?,
            (Some(left), None, false) => callback(index, left, None)?,
            _ => bail!("paired FASTQ files contain different record counts"),
        }
        index += 1;
    }

    if let Some(reader) = reader2.as_mut() {
        if !max_pairs.is_some_and(|limit| index >= limit) && reader.next_record()?.is_some() {
            bail!("paired FASTQ files contain different record counts");
        }
    }

    Ok(index)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_fastq() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "@r1\nACGT\n+\nIIII").unwrap();
        let mut reader = FastqReader::from_path(file.path()).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.id, b"r1");
        assert_eq!(record.sequence, b"ACGT");
        assert_eq!(record.quality, b"IIII");
        assert!(reader.next_record().unwrap().is_none());
    }
}
