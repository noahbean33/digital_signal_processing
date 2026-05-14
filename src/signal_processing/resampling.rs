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

// ─── Spectral Interpolation (Gap Filling) ────────────────────────────────────

/// Fill a gap in a signal using spectral interpolation.
///
/// Averages the FFTs of the segments immediately before and after the gap,
/// adds a linear ramp to match boundary values, and inserts the result.
///
/// * `signal` – the input signal with a gap (values in the gap are ignored).
/// * `gap_start` – first index of the gap (inclusive).
/// * `gap_end` – one past the last index of the gap (exclusive).
///
/// Returns a new signal with the gap filled.  If the gap or context is
/// invalid, returns the original signal unchanged.
#[must_use]
pub fn spectral_interpolation(signal: &[f64], gap_start: usize, gap_end: usize) -> Vec<f64> {
    let n = signal.len();
    if gap_end <= gap_start || gap_start == 0 || gap_end >= n {
        return signal.to_vec();
    }

    let gap_len = gap_end - gap_start;

    // We need context segments of the same length as the gap on either side
    let pre_start = if gap_start >= gap_len {
        gap_start - gap_len
    } else {
        0
    };
    let post_end = (gap_end + gap_len).min(n);

    let pre_seg = &signal[pre_start..gap_start];
    let post_seg = &signal[gap_end..post_end];

    // Use the shorter segment length for FFT
    let seg_len = pre_seg.len().min(post_seg.len());
    if seg_len < 2 {
        return signal.to_vec();
    }

    // Simple DFT of pre and post segments (real-valued, no external FFT needed)
    let dft = |data: &[f64], len: usize| -> Vec<(f64, f64)> {
        let mut result = Vec::with_capacity(len);
        for k in 0..len {
            let mut re = 0.0;
            let mut im = 0.0;
            for (n_i, &x) in data.iter().take(len).enumerate() {
                let angle = 2.0 * PI * k as f64 * n_i as f64 / len as f64;
                re += x * angle.cos();
                im -= x * angle.sin();
            }
            result.push((re, im));
        }
        result
    };

    let idft = |freq: &[(f64, f64)], len: usize| -> Vec<f64> {
        let mut result = Vec::with_capacity(len);
        let n_f = len as f64;
        for n_i in 0..len {
            let mut val = 0.0;
            for (k, &(re, im)) in freq.iter().enumerate() {
                let angle = 2.0 * PI * k as f64 * n_i as f64 / n_f;
                val += re * angle.cos() - im * angle.sin();
            }
            result.push(val / n_f);
        }
        result
    };

    let fft_pre = dft(pre_seg, seg_len);
    let fft_post = dft(post_seg, seg_len);

    // Average the two spectra
    let avg_fft: Vec<(f64, f64)> = fft_pre
        .iter()
        .zip(fft_post.iter())
        .map(|(&(r1, i1), &(r2, i2))| ((r1 + r2) / 2.0, (i1 + i2) / 2.0))
        .collect();

    // Inverse FFT → mixed data
    let mixed = idft(&avg_fft, seg_len);

    // Detrend the mixed data (remove its own linear trend)
    let mixed_mean = mixed.iter().sum::<f64>() / mixed.len() as f64;
    let detrended: Vec<f64> = mixed.iter().map(|&v| v - mixed_mean).collect();

    // Linear ramp from boundary values
    let left_val = signal[gap_start - 1];
    let right_val = if gap_end < n { signal[gap_end] } else { left_val };

    // Build the interpolated segment (gap_len samples from seg_len data)
    let mut result = signal.to_vec();
    for i in 0..gap_len {
        let t = (i as f64 + 1.0) / (gap_len as f64 + 1.0);
        let ramp = left_val * (1.0 - t) + right_val * t;
        let mixed_idx = if seg_len >= gap_len {
            i
        } else {
            i % seg_len
        };
        let mixed_val = if mixed_idx < detrended.len() {
            detrended[mixed_idx]
        } else {
            0.0
        };
        result[gap_start + i] = mixed_val + ramp;
    }
    result
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

    // ─── Spectral interpolation tests ────────────────────────────────────────

    #[test]
    fn spectral_interpolation_fills_gap() {
        // Create a smooth low-frequency signal, punch a small hole, and fill it
        let n = 200;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * 1.0 * i as f64 / n as f64).sin())
            .collect();

        // Zero out a small gap
        let mut gapped = signal.clone();
        for i in 90..110 {
            gapped[i] = 0.0;
        }

        let filled = spectral_interpolation(&gapped, 90, 110);
        assert_eq!(filled.len(), n);

        // The filled region should have non-trivial values (not all zero)
        let filled_energy: f64 = (90..110)
            .map(|i| filled[i].powi(2))
            .sum::<f64>();
        assert!(filled_energy > 0.01);
    }

    #[test]
    fn spectral_interpolation_invalid_gap_returns_original() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        // gap_start == 0 → invalid
        let result = spectral_interpolation(&signal, 0, 3);
        assert_eq!(result, signal);
        // gap_end >= n → invalid
        let result2 = spectral_interpolation(&signal, 2, 5);
        assert_eq!(result2, signal);
    }

    #[test]
    fn spectral_interpolation_preserves_outside_gap() {
        let signal: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let filled = spectral_interpolation(&signal, 30, 50);
        // Values outside the gap should be unchanged
        for i in 0..30 {
            assert!((filled[i] - signal[i]).abs() < 1e-10);
        }
        for i in 50..100 {
            assert!((filled[i] - signal[i]).abs() < 1e-10);
        }
    }
}
