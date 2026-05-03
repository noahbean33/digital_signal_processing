use std::f64::consts::PI;

use super::fft::Complex;

/// Calculate FIR low-pass filter coefficients using the windowed-sinc method
/// with a Hamming window.
///
/// `taps` – number of filter coefficients (filter order + 1).
/// `cutoff` – normalised cutoff frequency (0.0 … 0.5 of Nyquist).
#[must_use]
pub fn design_lowpass(taps: usize, cutoff: f64) -> Vec<f64> {
    if taps <= 1 {
        return Vec::new();
    }
    let m = taps - 1;
    let mut coeffs = Vec::with_capacity(taps);
    let mut sum = 0.0;
    for n in 0..taps {
        let pos = n as f64 - m as f64 / 2.0;
        let ideal = if pos.abs() < 1e-6 {
            1.0
        } else {
            (PI * cutoff * pos).sin() / (PI * cutoff * pos)
        };
        let window = 0.54 - 0.46 * (2.0 * PI * n as f64 / m as f64).cos();
        let c = ideal * window;
        coeffs.push(c);
        sum += c;
    }
    if sum.abs() > 1e-9 {
        for c in &mut coeffs {
            *c /= sum;
        }
    }
    coeffs
}

/// Calculate FIR high-pass filter coefficients via spectral inversion of a
/// windowed-sinc low-pass filter.
///
/// `taps` – number of filter coefficients (must be odd for a type-I high-pass).
/// `cutoff` – normalised cutoff frequency (0.0 … 0.5 of Nyquist).
#[must_use]
pub fn design_highpass(taps: usize, cutoff: f64) -> Vec<f64> {
    let mut lp = design_lowpass(taps, cutoff);
    if lp.is_empty() {
        return lp;
    }
    // Spectral inversion: negate all coefficients, add 1 to centre tap
    for c in lp.iter_mut() {
        *c = -*c;
    }
    let centre = (taps - 1) / 2;
    lp[centre] += 1.0;
    lp
}

/// Calculate FIR band-pass filter coefficients using the windowed-sinc method
/// with a Hamming window.
///
/// `taps` – number of filter coefficients (filter order + 1).
/// `low_cutoff` – normalised lower cutoff frequency (0.0 … 0.5).
/// `high_cutoff` – normalised upper cutoff frequency (low_cutoff … 0.5).
#[must_use]
pub fn design_bandpass(taps: usize, low_cutoff: f64, high_cutoff: f64) -> Vec<f64> {
    if taps <= 1 || low_cutoff >= high_cutoff {
        return Vec::new();
    }
    let m = taps - 1;
    let mut coeffs = Vec::with_capacity(taps);
    let mut sum = 0.0;
    for n in 0..taps {
        let pos = n as f64 - m as f64 / 2.0;
        let ideal = if pos.abs() < 1e-6 {
            2.0 * (high_cutoff - low_cutoff)
        } else {
            (PI * high_cutoff * pos).sin() / (PI * pos)
                - (PI * low_cutoff * pos).sin() / (PI * pos)
        };
        let window = 0.54 - 0.46 * (2.0 * PI * n as f64 / m as f64).cos();
        let c = ideal * window;
        coeffs.push(c);
        sum += c;
    }
    if sum.abs() > 1e-9 {
        for c in &mut coeffs {
            *c /= sum;
        }
    }
    coeffs
}

/// Calculate FIR band-stop (notch) filter coefficients using the windowed-sinc
/// method with a Hamming window.
///
/// `taps` – number of filter coefficients (filter order + 1).
/// `low_cutoff` – normalised lower cutoff frequency (0.0 … 0.5).
/// `high_cutoff` – normalised upper cutoff frequency (low_cutoff … 0.5).
#[must_use]
pub fn design_bandstop(taps: usize, low_cutoff: f64, high_cutoff: f64) -> Vec<f64> {
    if taps <= 1 || low_cutoff >= high_cutoff {
        return Vec::new();
    }
    let m = taps - 1;
    let mut coeffs = Vec::with_capacity(taps);
    let mut sum = 0.0;
    for n in 0..taps {
        let pos = n as f64 - m as f64 / 2.0;
        let ideal = if pos.abs() < 1e-6 {
            1.0 - 2.0 * (high_cutoff - low_cutoff)
        } else {
            (PI * low_cutoff * pos).sin() / (PI * pos)
                - (PI * high_cutoff * pos).sin() / (PI * pos)
        };
        let window = 0.54 - 0.46 * (2.0 * PI * n as f64 / m as f64).cos();
        let c = ideal * window;
        coeffs.push(c);
        sum += c;
    }
    if sum.abs() > 1e-9 {
        for c in &mut coeffs {
            *c /= sum;
        }
    }
    coeffs
}

/// Apply a FIR filter to an input signal via direct convolution (valid mode).
///
/// Output length = `signal.len() - coeffs.len() + 1`.
#[must_use]
pub fn apply(signal: &[f64], coeffs: &[f64]) -> Vec<f64> {
    if signal.len() < coeffs.len() || coeffs.is_empty() {
        return Vec::new();
    }
    let out_size = signal.len() - coeffs.len() + 1;
    let mut output = vec![0.0; out_size];
    for n in 0..out_size {
        let mut acc = 0.0;
        for k in 0..coeffs.len() {
            acc += signal[n + k] * coeffs[k];
        }
        output[n] = acc;
    }
    output
}

/// Apply a FIR filter in full-length mode (output length = signal length).
///
/// Past samples that fall before the signal start are treated as zero.
#[must_use]
pub fn apply_full(signal: &[f64], coeffs: &[f64]) -> Vec<f64> {
    if signal.is_empty() || coeffs.is_empty() {
        return Vec::new();
    }
    let mut output = vec![0.0; signal.len()];
    for i in 0..signal.len() {
        let mut acc = 0.0;
        for (j, &c) in coeffs.iter().enumerate() {
            if i >= j {
                acc += signal[i - j] * c;
            }
        }
        output[i] = acc;
    }
    output
}

/// Compute the complex frequency response H(e^{jω}) of a FIR filter at a
/// given frequency.
#[must_use]
pub fn frequency_response(coeffs: &[f64], frequency: f64, sampling_rate: f64) -> Complex {
    let mut response = Complex::zero();
    for (n, &c) in coeffs.iter().enumerate() {
        let angle = -2.0 * PI * frequency * n as f64 / sampling_rate;
        response += Complex::new(c, 0.0) * Complex::new(angle.cos(), angle.sin());
    }
    response
}

/// Quantise FIR coefficients to fixed-point representation.
#[must_use]
pub fn quantize_coefficients(coeffs: &[f64], q_factor: i32) -> Vec<i32> {
    coeffs
        .iter()
        .map(|&c| (c * f64::from(q_factor)).round() as i32)
        .collect()
}

// ─── Streaming / Stateful FIR Filter ─────────────────────────────────────────

/// A stateful FIR filter that processes one sample at a time using a circular
/// buffer, retaining state between calls.
pub struct StreamingFir {
    coefficients: Vec<f64>,
    buffer: Vec<f64>,
    index: usize,
}

impl StreamingFir {
    /// Create a new streaming FIR filter with the given coefficients.
    #[must_use]
    pub fn new(coefficients: &[f64]) -> Self {
        let len = coefficients.len();
        Self {
            coefficients: coefficients.to_vec(),
            buffer: vec![0.0; len],
            index: 0,
        }
    }

    /// Reset internal buffer to zero.
    pub fn reset(&mut self) {
        self.buffer.fill(0.0);
        self.index = 0;
    }

    /// Process a single input sample and return the filtered output.
    pub fn process_sample(&mut self, sample: f64) -> f64 {
        let len = self.coefficients.len();
        if len == 0 {
            return 0.0;
        }
        // Write sample into circular buffer
        self.buffer[self.index] = sample;

        // Convolve: walk backwards through the buffer using the circular index
        let mut acc = 0.0;
        let mut buf_idx = self.index;
        for k in 0..len {
            acc += self.coefficients[k] * self.buffer[buf_idx];
            // Wrap-around decrement
            if buf_idx == 0 {
                buf_idx = len - 1;
            } else {
                buf_idx -= 1;
            }
        }

        // Advance circular index
        self.index = (self.index + 1) % len;

        acc
    }

    /// Process an entire signal through the streaming filter.
    #[must_use]
    pub fn apply(&mut self, signal: &[f64]) -> Vec<f64> {
        signal.iter().map(|&s| self.process_sample(s)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lowpass_coefficients_sum_to_one() {
        let coeffs = design_lowpass(51, 0.25);
        let sum: f64 = coeffs.iter().sum();
        assert!((sum - 1.0).abs() < 1e-6);
    }

    #[test]
    fn apply_impulse_response() {
        let coeffs = vec![0.25, 0.5, 0.25];
        let signal = vec![0.0, 0.0, 1.0, 0.0, 0.0];
        let out = apply(&signal, &coeffs);
        assert_eq!(out.len(), 3);
        assert!((out[0] - 0.25).abs() < 1e-10);
        assert!((out[1] - 0.5).abs() < 1e-10);
        assert!((out[2] - 0.25).abs() < 1e-10);
    }

    #[test]
    fn quantize_roundtrip() {
        let coeffs = vec![0.25, 0.5, 0.25];
        let q = quantize_coefficients(&coeffs, 1000);
        assert_eq!(q, vec![250, 500, 250]);
    }

    #[test]
    fn streaming_fir_matches_batch() {
        let coeffs = vec![0.25, 0.5, 0.25];
        let signal = vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        let batch = apply_full(&signal, &coeffs);
        let mut fir = StreamingFir::new(&coeffs);
        let streaming = fir.apply(&signal);
        assert_eq!(batch.len(), streaming.len());
        for (a, b) in batch.iter().zip(streaming.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn streaming_fir_sample_by_sample() {
        let coeffs = vec![1.0, 0.0, 0.0];
        let mut fir = StreamingFir::new(&coeffs);
        assert!((fir.process_sample(5.0) - 5.0).abs() < 1e-10);
        assert!((fir.process_sample(3.0) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn streaming_fir_reset() {
        let coeffs = vec![0.5, 0.5];
        let mut fir = StreamingFir::new(&coeffs);
        fir.process_sample(10.0);
        fir.reset();
        // After reset, buffer should be zeroed
        let out = fir.process_sample(2.0);
        assert!((out - 1.0).abs() < 1e-10); // 0.5*2 + 0.5*0
    }

    #[test]
    fn highpass_attenuates_dc() {
        let coeffs = design_highpass(51, 0.25);
        assert_eq!(coeffs.len(), 51);
        // Sum of coefficients should be near zero (rejects DC)
        let sum: f64 = coeffs.iter().sum();
        assert!(sum.abs() < 0.05);
    }

    #[test]
    fn bandpass_correct_length() {
        let coeffs = design_bandpass(51, 0.1, 0.3);
        assert_eq!(coeffs.len(), 51);
    }

    #[test]
    fn bandstop_correct_length() {
        let coeffs = design_bandstop(51, 0.1, 0.3);
        assert_eq!(coeffs.len(), 51);
    }

    #[test]
    fn bandpass_invalid_cutoffs_empty() {
        let coeffs = design_bandpass(51, 0.3, 0.1);
        assert!(coeffs.is_empty());
    }
}
