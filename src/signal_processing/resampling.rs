use std::f64::consts::PI;

/// Normalised sinc function: sinc(x) = sin(πx) / (πx).
#[must_use]
pub fn sinc(x: f64) -> f64 {
    if x.abs() < 1e-7 {
        return 1.0;
    }
    (PI * x).sin() / (PI * x)
}

/// Generate a FIR low-pass filter kernel using the windowed-sinc method
/// (Hamming window).
///
/// `num_taps` – filter length.
/// `cutoff` – normalised cutoff frequency (0.0 … 0.5).
#[must_use]
pub fn generate_fir_lowpass(num_taps: usize, cutoff: f64) -> Vec<f64> {
    if num_taps <= 1 {
        return Vec::new();
    }
    let m = num_taps - 1;
    let mut kernel = Vec::with_capacity(num_taps);
    for i in 0..num_taps {
        let n = i as f64 - m as f64 / 2.0;
        let h = 2.0 * cutoff * sinc(2.0 * cutoff * n);
        let w = 0.54 - 0.46 * (2.0 * PI * i as f64 / m as f64).cos();
        kernel.push(h * w);
    }
    kernel
}

/// Decimate a signal by keeping every `factor`-th sample.
#[must_use]
pub fn decimate(input: &[f64], factor: usize) -> Vec<f64> {
    if factor == 0 || input.is_empty() {
        return Vec::new();
    }
    input.iter().step_by(factor).copied().collect()
}

/// Interpolate (upsample) by inserting `factor - 1` zeros between each sample.
#[must_use]
pub fn interpolate(input: &[f64], factor: usize) -> Vec<f64> {
    if factor == 0 || input.is_empty() {
        return Vec::new();
    }
    let mut output = vec![0.0; input.len() * factor];
    for (i, &v) in input.iter().enumerate() {
        output[i * factor] = v;
    }
    output
}

/// Apply a FIR filter (causal convolution, output length = input length).
#[must_use]
pub fn apply_fir_filter(signal: &[f64], kernel: &[f64]) -> Vec<f64> {
    if signal.is_empty() || kernel.is_empty() {
        return Vec::new();
    }
    let signal_size = signal.len();
    let kernel_size = kernel.len();
    let mut filtered = vec![0.0; signal_size];
    for i in 0..signal_size {
        let mut sum = 0.0;
        for j in 0..kernel_size {
            if i >= j {
                sum += signal[i - j] * kernel[j];
            }
        }
        filtered[i] = sum;
    }
    filtered
}

/// Rational sample-rate conversion by a factor of L/M.
///
/// Upsamples by `l`, applies an anti-aliasing FIR filter, then decimates by `m`.
#[must_use]
pub fn sample_rate_conversion(
    input: &[f64],
    l: usize,
    m: usize,
    fir_kernel: &[f64],
) -> Vec<f64> {
    if l == 0 || m == 0 || input.is_empty() {
        return Vec::new();
    }
    let upsampled = interpolate(input, l);
    let filtered = apply_fir_filter(&upsampled, fir_kernel);
    decimate(&filtered, m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sinc_at_zero_is_one() {
        assert!((sinc(0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn decimate_by_2() {
        let input = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let result = decimate(&input, 2);
        assert_eq!(result, vec![1.0, 3.0, 5.0]);
    }

    #[test]
    fn interpolate_by_3() {
        let input = vec![1.0, 2.0];
        let result = interpolate(&input, 3);
        assert_eq!(result, vec![1.0, 0.0, 0.0, 2.0, 0.0, 0.0]);
    }

    #[test]
    fn src_preserves_approximate_length() {
        let signal: Vec<f64> = (0..100)
            .map(|n| (2.0 * PI * 50.0 * n as f64 / 1000.0).sin())
            .collect();
        let kernel = generate_fir_lowpass(29, 0.25);
        let converted = sample_rate_conversion(&signal, 3, 2, &kernel);
        // L/M = 3/2, so output ≈ 150 samples (minus filter transient)
        assert!(converted.len() > 100);
    }
}
