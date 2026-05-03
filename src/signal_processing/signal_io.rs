use std::fs;
use std::io::{self, BufRead, Write};
use std::path::Path;

/// Read a signal from a text `.dat` file (one sample per line).
///
/// # Errors
///
/// Returns an I/O error if the file cannot be opened or read.
pub fn read_dat<P: AsRef<Path>>(path: P) -> io::Result<Vec<f64>> {
    let file = fs::File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut signal = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        match trimmed.parse::<f64>() {
            Ok(value) => signal.push(value),
            Err(_) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Could not parse '{trimmed}' as f64"),
                ));
            }
        }
    }
    Ok(signal)
}

/// Write a signal to a text `.dat` file (one sample per line).
///
/// # Errors
///
/// Returns an I/O error if the file cannot be created or written to.
pub fn write_dat<P: AsRef<Path>>(path: P, signal: &[f64]) -> io::Result<()> {
    let mut file = fs::File::create(path)?;
    for &sample in signal {
        writeln!(file, "{sample}")?;
    }
    Ok(())
}

/// Read a complex signal from two `.dat` files (real and imaginary parts).
///
/// # Errors
///
/// Returns an I/O error if either file cannot be opened or read.
pub fn read_complex_dat<P: AsRef<Path>>(
    real_path: P,
    imag_path: P,
) -> io::Result<(Vec<f64>, Vec<f64>)> {
    let real = read_dat(real_path)?;
    let imag = read_dat(imag_path)?;
    Ok((real, imag))
}

/// Write a complex signal to two `.dat` files (real and imaginary parts).
///
/// # Errors
///
/// Returns an I/O error if either file cannot be created or written to.
pub fn write_complex_dat<P: AsRef<Path>>(
    real_path: P,
    imag_path: P,
    real: &[f64],
    imag: &[f64],
) -> io::Result<()> {
    write_dat(real_path, real)?;
    write_dat(imag_path, imag)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn write_and_read_roundtrip() {
        let dir = std::env::temp_dir();
        let path = dir.join("test_signal_io_roundtrip.dat");

        let signal = vec![1.0, -2.5, 3.14159, 0.0, 42.0];
        write_dat(&path, &signal).expect("write failed");

        let loaded = read_dat(&path).expect("read failed");
        assert_eq!(loaded.len(), signal.len());
        for (a, b) in signal.iter().zip(loaded.iter()) {
            assert!((a - b).abs() < 1e-10);
        }

        let _ = fs::remove_file(&path);
    }

    #[test]
    fn read_dat_skips_blank_lines() {
        let dir = std::env::temp_dir();
        let path = dir.join("test_signal_io_blanks.dat");

        fs::write(&path, "1.0\n\n2.0\n  \n3.0\n").expect("write failed");

        let loaded = read_dat(&path).expect("read failed");
        assert_eq!(loaded, vec![1.0, 2.0, 3.0]);

        let _ = fs::remove_file(&path);
    }

    #[test]
    fn read_dat_nonexistent_file_errors() {
        let result = read_dat("nonexistent_file_12345.dat");
        assert!(result.is_err());
    }
}
