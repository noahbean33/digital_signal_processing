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
}
