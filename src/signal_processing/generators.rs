use std::f64::consts::PI;

use super::fft::Complex;

/// Generate a sine wave signal.
///
/// `x[n] = amplitude * sin(2π * frequency * n / sample_rate + phase)`
#[must_use]
pub fn generate_sine(
    length: usize,
    amplitude: f64,
    frequency: f64,
    sample_rate: f64,
    phase: f64,
) -> Vec<f64> {
    (0..length)
        .map(|n| amplitude * (2.0 * PI * frequency * n as f64 / sample_rate + phase).sin())
        .collect()
}

/// Generate a cosine wave signal.
///
/// `x[n] = amplitude * cos(2π * frequency * n / sample_rate + phase)`
#[must_use]
pub fn generate_cosine(
    length: usize,
    amplitude: f64,
    frequency: f64,
    sample_rate: f64,
    phase: f64,
) -> Vec<f64> {
    (0..length)
        .map(|n| amplitude * (2.0 * PI * frequency * n as f64 / sample_rate + phase).cos())
        .collect()
}

/// Generate an impulse (delta) signal with a single non-zero sample.
#[must_use]
pub fn generate_impulse(length: usize, position: usize, amplitude: f64) -> Vec<f64> {
    let mut signal = vec![0.0; length];
    if position < length {
        signal[position] = amplitude;
    }
    signal
}

/// Generate a step (Heaviside) signal: zero before `position`, `amplitude` from
/// `position` onward.
#[must_use]
pub fn generate_step(length: usize, position: usize, amplitude: f64) -> Vec<f64> {
    let mut signal = vec![0.0; length];
    for s in signal.iter_mut().skip(position) {
        *s = amplitude;
    }
    signal
}

/// Generate uniform white noise in the range `[-amplitude, +amplitude]`.
///
/// Uses a simple linear congruential generator seeded by `seed`.
#[must_use]
pub fn generate_noise(length: usize, amplitude: f64, seed: u64) -> Vec<f64> {
    let mut state = seed.wrapping_add(1);
    (0..length)
        .map(|_| {
            // xorshift64
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let uniform = (state as f64) / (u64::MAX as f64) * 2.0 - 1.0;
            amplitude * uniform
        })
        .collect()
}

/// Generate a chirp (frequency-sweep) signal from `f0` to `f1` over `length`
/// samples at the given `sample_rate`.
#[must_use]
pub fn generate_chirp(
    length: usize,
    amplitude: f64,
    f0: f64,
    f1: f64,
    sample_rate: f64,
) -> Vec<f64> {
    if length == 0 {
        return Vec::new();
    }
    let duration = length as f64 / sample_rate;
    let k = (f1 - f0) / duration;
    (0..length)
        .map(|n| {
            let t = n as f64 / sample_rate;
            amplitude * (2.0 * PI * (f0 * t + 0.5 * k * t * t)).sin()
        })
        .collect()
}

/// Normalise a signal to the range `[-1, 1]`.
#[must_use]
pub fn normalize(signal: &[f64]) -> Vec<f64> {
    if signal.is_empty() {
        return Vec::new();
    }
    let max_abs = signal
        .iter()
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max);
    if max_abs < 1e-30 {
        return signal.to_vec();
    }
    signal.iter().map(|&x| x / max_abs).collect()
}

/// Scale every sample by a constant factor.
#[must_use]
pub fn scale(signal: &[f64], factor: f64) -> Vec<f64> {
    signal.iter().map(|&x| x * factor).collect()
}

/// Add a DC offset to every sample.
#[must_use]
pub fn add_offset(signal: &[f64], offset: f64) -> Vec<f64> {
    signal.iter().map(|&x| x + offset).collect()
}

/// Remove the DC component (subtract the mean).
#[must_use]
pub fn remove_dc(signal: &[f64]) -> Vec<f64> {
    if signal.is_empty() {
        return Vec::new();
    }
    let m = signal.iter().sum::<f64>() / signal.len() as f64;
    signal.iter().map(|&x| x - m).collect()
}

/// Convert rectangular (real, imag) representation to polar (magnitude, phase).
///
/// Returns `(magnitudes, phases)`.
#[must_use]
pub fn rect_to_polar(real: &[f64], imag: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let len = real.len().min(imag.len());
    let mut magnitudes = Vec::with_capacity(len);
    let mut phases = Vec::with_capacity(len);
    for i in 0..len {
        let c = Complex::new(real[i], imag[i]);
        magnitudes.push(c.norm());
        if real[i] == 0.0 && imag[i] == 0.0 {
            phases.push(0.0);
        } else {
            phases.push(c.arg());
        }
    }
    (magnitudes, phases)
}

/// Convert polar (magnitude, phase) representation to rectangular (real, imag).
///
/// Returns `(real, imag)`.
#[must_use]
pub fn polar_to_rect(magnitudes: &[f64], phases: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let len = magnitudes.len().min(phases.len());
    let mut real = Vec::with_capacity(len);
    let mut imag = Vec::with_capacity(len);
    for i in 0..len {
        real.push(magnitudes[i] * phases[i].cos());
        imag.push(magnitudes[i] * phases[i].sin());
    }
    (real, imag)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sine_starts_at_zero() {
        let s = generate_sine(100, 1.0, 10.0, 100.0, 0.0);
        assert!(s[0].abs() < 1e-10);
    }

    #[test]
    fn cosine_starts_at_amplitude() {
        let s = generate_cosine(100, 2.0, 10.0, 100.0, 0.0);
        assert!((s[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn impulse_has_one_nonzero() {
        let s = generate_impulse(10, 5, 3.0);
        assert!((s[5] - 3.0).abs() < 1e-10);
        let sum: f64 = s.iter().sum();
        assert!((sum - 3.0).abs() < 1e-10);
    }

    #[test]
    fn step_correct() {
        let s = generate_step(6, 3, 1.0);
        assert_eq!(s, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    }

    #[test]
    fn noise_has_correct_length() {
        let n = generate_noise(50, 1.0, 42);
        assert_eq!(n.len(), 50);
    }

    #[test]
    fn normalize_peaks_at_one() {
        let s = vec![2.0, -4.0, 1.0];
        let n = normalize(&s);
        assert!((n[1] + 1.0).abs() < 1e-10);
        assert!(n.iter().all(|&v| v.abs() <= 1.0 + 1e-10));
    }

    #[test]
    fn scale_doubles() {
        let s = vec![1.0, 2.0, 3.0];
        let scaled = scale(&s, 2.0);
        assert_eq!(scaled, vec![2.0, 4.0, 6.0]);
    }

    #[test]
    fn add_offset_shifts() {
        let s = vec![1.0, 2.0, 3.0];
        let shifted = add_offset(&s, 10.0);
        assert_eq!(shifted, vec![11.0, 12.0, 13.0]);
    }

    #[test]
    fn remove_dc_zeroes_mean() {
        let s = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let dc_removed = remove_dc(&s);
        let m: f64 = dc_removed.iter().sum::<f64>() / dc_removed.len() as f64;
        assert!(m.abs() < 1e-10);
    }

    #[test]
    fn rect_polar_roundtrip() {
        let re = vec![1.0, 0.0, -1.0];
        let im = vec![0.0, 1.0, 0.0];
        let (mag, phase) = rect_to_polar(&re, &im);
        let (re2, im2) = polar_to_rect(&mag, &phase);
        for i in 0..3 {
            assert!((re[i] - re2[i]).abs() < 1e-10);
            assert!((im[i] - im2[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn chirp_correct_length() {
        let c = generate_chirp(200, 1.0, 100.0, 1000.0, 8000.0);
        assert_eq!(c.len(), 200);
    }
}
