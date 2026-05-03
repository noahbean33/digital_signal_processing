use super::fft::{self, Complex};
use super::windowing;

/// Result of one STFT frame analysis.
#[derive(Debug, Clone)]
pub struct StftFrame {
    pub start_sample: usize,
    pub spectrum: Vec<Complex>,
    pub magnitude: Vec<f64>,
}

/// Compute the Short-Time Fourier Transform of a signal.
///
/// `window_size` – length of each analysis window.
/// `hop_size` – step between successive frames.
///
/// Returns a vector of `StftFrame`s.
#[must_use]
pub fn stft(signal: &[f64], window_size: usize, hop_size: usize) -> Vec<StftFrame> {
    if signal.is_empty() || window_size == 0 || hop_size == 0 {
        return Vec::new();
    }
    let window = windowing::hann(window_size);
    let mut frames = Vec::new();

    let mut start = 0;
    while start + window_size <= signal.len() {
        let segment = extract_segment(signal, start, window_size);
        let windowed = windowing::apply_window(&segment, &window);
        let spectrum = perform_fft(&windowed);
        let magnitude = fft::magnitude_spectrum(&spectrum);

        frames.push(StftFrame {
            start_sample: start,
            spectrum,
            magnitude,
        });
        start += hop_size;
    }
    frames
}

/// Extract a segment from the signal, zero-padding if necessary.
#[must_use]
pub fn extract_segment(signal: &[f64], start: usize, window_size: usize) -> Vec<f64> {
    let mut segment = vec![0.0; window_size];
    for i in 0..window_size {
        if start + i < signal.len() {
            segment[i] = signal[start + i];
        }
    }
    segment
}

/// Perform FFT on a real-valued segment (zero-pads to next power of 2).
#[must_use]
pub fn perform_fft(segment: &[f64]) -> Vec<Complex> {
    if segment.is_empty() {
        return Vec::new();
    }
    let mut data = fft::real_to_complex(segment);
    data = fft::zero_pad_to_power_of_2(&data);
    fft::fft(&mut data);
    data
}

/// Detect spectral peaks that exceed a threshold in the magnitude spectrum.
///
/// Returns a vector of bin indices where peaks were found.
#[must_use]
pub fn detect_peaks(magnitude: &[f64], threshold: f64) -> Vec<usize> {
    if magnitude.len() < 3 {
        return Vec::new();
    }
    let mut peak_indices = Vec::new();
    for i in 1..magnitude.len() - 1 {
        if magnitude[i] > threshold
            && magnitude[i] > magnitude[i - 1]
            && magnitude[i] > magnitude[i + 1]
        {
            peak_indices.push(i);
        }
    }
    peak_indices
}

/// Compute the spectrogram (magnitude² matrix) from STFT frames.
///
/// Returns a 2D vector: `[frame_index][frequency_bin]`.
#[must_use]
pub fn spectrogram(frames: &[StftFrame]) -> Vec<Vec<f64>> {
    frames
        .iter()
        .map(|f| f.magnitude.iter().map(|m| m * m).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stft_produces_frames() {
        let signal: Vec<f64> = (0..256)
            .map(|i| (2.0 * std::f64::consts::PI * 10.0 * i as f64 / 256.0).sin())
            .collect();
        let frames = stft(&signal, 64, 32);
        assert!(!frames.is_empty());
        assert!(frames[0].magnitude.len() > 0);
    }

    #[test]
    fn detect_peaks_finds_local_maxima() {
        let mag = vec![0.0, 1.0, 5.0, 3.0, 1.0, 8.0, 2.0];
        let peaks = detect_peaks(&mag, 2.0);
        assert!(peaks.contains(&2));
        assert!(peaks.contains(&5));
    }

    #[test]
    fn extract_segment_zero_pads() {
        let signal = vec![1.0, 2.0, 3.0];
        let seg = extract_segment(&signal, 1, 5);
        assert_eq!(seg, vec![2.0, 3.0, 0.0, 0.0, 0.0]);
    }
}
