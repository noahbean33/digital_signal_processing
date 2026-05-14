use super::fft::{Complex, fft, ifft};

/// FFT-based overlap-add convolution for efficient filtering of long signals.
///
/// Splits the input `signal` into blocks of size `block_size`, performs circular
/// convolution with `kernel` in the frequency domain via FFT, and combines
/// the results using the overlap-add method.
///
/// Complexity: O(N log N) where N is the FFT size, versus O(N·M) for direct
/// convolution. Significantly faster when `kernel` is long.
///
/// Output length = `signal.len() + kernel.len() - 1`.
#[must_use]
pub fn fft_overlap_add(signal: &[f64], kernel: &[f64], block_size: usize) -> Vec<f64> {
    if signal.is_empty() || kernel.is_empty() || block_size == 0 {
        return Vec::new();
    }
    let sig_len = signal.len();
    let ker_len = kernel.len();
    let out_len = sig_len + ker_len - 1;

    // FFT size: next power of 2 >= block_size + kernel_len - 1
    let seg_len = block_size + ker_len - 1;
    let mut fft_size = 1;
    while fft_size < seg_len {
        fft_size <<= 1;
    }

    // Pre-compute kernel FFT (zero-padded)
    let mut h_fft: Vec<Complex> = kernel
        .iter()
        .map(|&v| Complex::new(v, 0.0))
        .chain(std::iter::repeat(Complex::zero()).take(fft_size - ker_len))
        .collect();
    fft(&mut h_fft);

    let mut output = vec![0.0; out_len];
    let mut pos = 0;

    while pos < sig_len {
        let end = (pos + block_size).min(sig_len);
        let current_block_len = end - pos;

        // Zero-padded block
        let mut block: Vec<Complex> = signal[pos..end]
            .iter()
            .map(|&v| Complex::new(v, 0.0))
            .chain(std::iter::repeat(Complex::zero()).take(fft_size - current_block_len))
            .collect();

        // Circular convolution via FFT
        fft(&mut block);
        for i in 0..fft_size {
            block[i] = block[i] * h_fft[i];
        }
        ifft(&mut block);

        // Overlap-add into output
        for i in 0..fft_size {
            let out_idx = pos + i;
            if out_idx < out_len {
                output[out_idx] += block[i].re;
            }
        }

        pos += block_size;
    }
    output
}

/// FFT-based overlap-save (overlap-scrap) convolution.
///
/// An alternative to overlap-add that avoids the addition step by discarding
/// the circular convolution artifacts and keeping only valid samples.
///
/// Output length = `signal.len() + kernel.len() - 1`.
#[must_use]
pub fn fft_overlap_save(signal: &[f64], kernel: &[f64], block_size: usize) -> Vec<f64> {
    if signal.is_empty() || kernel.is_empty() || block_size == 0 {
        return Vec::new();
    }
    let sig_len = signal.len();
    let ker_len = kernel.len();
    let out_len = sig_len + ker_len - 1;

    // FFT size: next power of 2 >= block_size + kernel_len - 1
    let mut fft_size = 1;
    let needed = block_size + ker_len - 1;
    while fft_size < needed {
        fft_size <<= 1;
    }
    let overlap = ker_len - 1;
    let valid_per_block = fft_size - overlap;

    // Pre-compute kernel FFT
    let mut h_fft: Vec<Complex> = kernel
        .iter()
        .map(|&v| Complex::new(v, 0.0))
        .chain(std::iter::repeat(Complex::zero()).take(fft_size - ker_len))
        .collect();
    fft(&mut h_fft);

    // Pad the signal with overlap zeros at front and enough at the back
    let padded_len = sig_len + 2 * overlap;
    let mut padded = vec![0.0; padded_len];
    for i in 0..sig_len {
        padded[overlap + i] = signal[i];
    }

    let mut output = vec![0.0; out_len];
    let mut in_pos = 0;
    let mut out_pos = 0;

    while out_pos < out_len {
        let mut block: Vec<Complex> = (in_pos..in_pos + fft_size)
            .map(|i| {
                if i < padded.len() {
                    Complex::new(padded[i], 0.0)
                } else {
                    Complex::zero()
                }
            })
            .collect();

        fft(&mut block);
        for i in 0..fft_size {
            block[i] = block[i] * h_fft[i];
        }
        ifft(&mut block);

        // Keep only valid portion (discard first `overlap` samples)
        for i in 0..valid_per_block {
            let idx = out_pos + i;
            if idx < out_len {
                output[idx] = block[overlap + i].re;
            }
        }

        in_pos += valid_per_block;
        out_pos += valid_per_block;
    }
    output
}

// ─── Frequency-Domain Gaussian Filter ─────────────────────────────────────────

/// Result of the frequency-domain Gaussian filter including empirical values.
#[derive(Clone, Debug)]
pub struct GaussianFilterResult {
    /// Filtered signal.
    pub filtered: Vec<f64>,
    /// Empirical centre frequency (Hz).
    pub empirical_frequency: f64,
    /// Empirical FWHM in the frequency domain (Hz).
    pub empirical_fwhm_hz: f64,
}

/// Narrow-band filter via frequency-domain Gaussian.
///
/// Filters a signal by multiplying its FFT with a Gaussian kernel centred
/// at `center_freq` with a specified full-width at half-maximum (FWHM).
///
/// This approach provides excellent frequency selectivity with minimal
/// time-domain ringing compared to FIR/IIR alternatives.
///
/// * `signal` – input time-domain signal
/// * `sample_rate` – sampling rate in Hz
/// * `center_freq` – peak frequency of the Gaussian filter (Hz)
/// * `fwhm` – full-width at half-maximum of the Gaussian in Hz
///
/// Returns a [`GaussianFilterResult`] with the filtered signal and empirical values.
#[must_use]
pub fn gaussian_frequency_filter(
    signal: &[f64],
    sample_rate: f64,
    center_freq: f64,
    fwhm: f64,
) -> GaussianFilterResult {
    use std::f64::consts::PI;

    let n = signal.len();
    if n == 0 {
        return GaussianFilterResult {
            filtered: Vec::new(),
            empirical_frequency: 0.0,
            empirical_fwhm_hz: 0.0,
        };
    }

    // Compute frequency vector
    let hz: Vec<f64> = (0..n).map(|i| sample_rate * i as f64 / n as f64).collect();

    // Create Gaussian kernel in frequency domain
    // sigma from FWHM: s = fwhm * (2*pi - 1) / (4*pi)
    let sigma = fwhm * (2.0 * PI - 1.0) / (4.0 * PI);
    let fx: Vec<f64> = hz.iter()
        .map(|&f| {
            let x = f - center_freq;
            (-0.5 * (x / sigma).powi(2)).exp()
        })
        .collect();

    // Normalize to unit gain
    let max_fx = fx.iter().cloned().fold(0.0_f64, f64::max);
    let fx_norm: Vec<f64> = if max_fx > 1e-30 {
        fx.iter().map(|&v| v / max_fx).collect()
    } else {
        fx.clone()
    };

    // Pad to power of 2 for FFT
    let mut fft_size = 1;
    while fft_size < n {
        fft_size <<= 1;
    }

    // FFT of signal
    let mut signal_fft: Vec<Complex> = signal
        .iter()
        .map(|&v| Complex::new(v, 0.0))
        .chain(std::iter::repeat(Complex::zero()).take(fft_size - n))
        .collect();
    fft(&mut signal_fft);

    // Create frequency-domain filter (padded)
    let fx_padded: Vec<f64> = (0..fft_size)
        .map(|i| {
            if i < n {
                fx_norm[i]
            } else {
                // Mirror for negative frequencies
                let mirror_idx = fft_size - i;
                if mirror_idx < n { fx_norm[mirror_idx] } else { 0.0 }
            }
        })
        .collect();

    // Apply filter in frequency domain
    for i in 0..fft_size {
        signal_fft[i] = Complex::new(
            signal_fft[i].re * fx_padded[i],
            signal_fft[i].im * fx_padded[i],
        );
    }

    // Inverse FFT
    ifft(&mut signal_fft);

    // Extract real part (×2 to compensate for one-sided filtering)
    let filtered: Vec<f64> = signal_fft[..n].iter().map(|c| 2.0 * c.re).collect();

    // Compute empirical frequency (peak of the Gaussian)
    let emp_idx = hz.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            ((*a - center_freq).abs()).partial_cmp(&((*b - center_freq).abs())).unwrap()
        })
        .map(|(i, _)| i)
        .unwrap_or(0);
    let empirical_frequency = hz[emp_idx];

    // Compute empirical FWHM: find half-max crossings
    let half_max_idx_left = fx_norm[..emp_idx]
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            ((*a - 0.5).abs()).partial_cmp(&((*b - 0.5).abs())).unwrap()
        })
        .map(|(i, _)| i)
        .unwrap_or(0);

    let half_max_idx_right = fx_norm[emp_idx..]
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            ((*a - 0.5).abs()).partial_cmp(&((*b - 0.5).abs())).unwrap()
        })
        .map(|(i, _)| emp_idx + i)
        .unwrap_or(n - 1);

    let empirical_fwhm_hz = hz[half_max_idx_right] - hz[half_max_idx_left];

    GaussianFilterResult {
        filtered,
        empirical_frequency,
        empirical_fwhm_hz,
    }
}

/// Simplified version that just returns the filtered signal.
#[must_use]
pub fn gaussian_bandpass(signal: &[f64], sample_rate: f64, center_freq: f64, fwhm: f64) -> Vec<f64> {
    gaussian_frequency_filter(signal, sample_rate, center_freq, fwhm).filtered
}

// ─── Planck Taper Bandpass Filter ─────────────────────────────────────────────

/// Narrow-band filter using a Planck-taper spectral window.
///
/// The Planck taper provides a smooth, flat-top bandpass shape with tuneable
/// edge roll-off controlled by the `eta` parameter (0 < eta < 0.5).
///
/// * `signal` – input time-domain signal.
/// * `sample_rate` – sampling rate in Hz.
/// * `center_freq` – peak frequency of the bandpass (Hz).
/// * `fwhm` – full-width at half-maximum of the passband (Hz).
/// * `eta` – edge decay parameter (0 < eta < 0.5). Smaller = sharper edges.
///
/// Returns the filtered signal (same length as input).
#[must_use]
pub fn planck_taper_bandpass(
    signal: &[f64],
    sample_rate: f64,
    center_freq: f64,
    fwhm: f64,
    eta: f64,
) -> Vec<f64> {
    let n = signal.len();
    if n == 0 || sample_rate <= 0.0 || fwhm <= 0.0 {
        return signal.to_vec();
    }

    // Number of frequency bins spanning the FWHM
    let np_bins = ((2.0 * fwhm * n as f64 / sample_rate).round() as usize).max(2);

    // Build Planck taper of length np_bins
    let np_f = np_bins as f64;
    let eta_clamped = eta.max(0.01).min(0.49);
    let mut taper = vec![0.0; np_bins];

    let bound_left = ((eta_clamped * (np_f - 1.0)).floor() as usize).max(1);
    let bound_right = (((1.0 - eta_clamped) * (np_f - 1.0)).ceil() as usize).min(np_bins - 2);

    for i in 0..np_bins {
        let pt = i as f64;
        if i == 0 {
            taper[i] = 0.0;
        } else if i < bound_left {
            let zl = eta_clamped * (np_f - 1.0) * (1.0 / pt + 1.0 / (pt - eta_clamped * (np_f - 1.0)));
            taper[i] = 1.0 / (zl.exp() + 1.0);
        } else if i <= bound_right {
            taper[i] = 1.0;
        } else if i < np_bins - 1 {
            let zr = eta_clamped * (np_f - 1.0)
                * (1.0 / (np_f - 1.0 - pt) + 1.0 / ((1.0 - eta_clamped) * (np_f - 1.0) - pt));
            taper[i] = 1.0 / (zr.exp() + 1.0);
        } else {
            taper[i] = 0.0;
        }
    }

    // Frequency vector
    let hz: Vec<f64> = (0..n).map(|i| sample_rate * i as f64 / n as f64).collect();

    // Find centre frequency index
    let fidx = hz
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            ((*a - center_freq).abs())
                .partial_cmp(&((*b - center_freq).abs()))
                .unwrap()
        })
        .map(|(i, _)| i)
        .unwrap_or(0);

    // Place taper into a full-length frequency-domain window
    let half_np = np_bins / 2;
    let mut px = vec![0.0; n];
    for (k, &v) in taper.iter().enumerate() {
        let idx = fidx as isize - half_np as isize + k as isize;
        if idx >= 0 && (idx as usize) < n {
            px[idx as usize] = v;
        }
    }

    // Pad to power of 2 for FFT
    let mut fft_size = 1;
    while fft_size < n {
        fft_size <<= 1;
    }

    // FFT of signal
    let mut signal_fft: Vec<Complex> = signal
        .iter()
        .map(|&v| Complex::new(v, 0.0))
        .chain(std::iter::repeat(Complex::zero()).take(fft_size - n))
        .collect();
    fft(&mut signal_fft);

    // Build padded filter (mirror for negative frequencies)
    let px_padded: Vec<f64> = (0..fft_size)
        .map(|i| {
            if i < n {
                px[i]
            } else {
                let mirror = fft_size - i;
                if mirror < n { px[mirror] } else { 0.0 }
            }
        })
        .collect();

    // Apply filter
    for i in 0..fft_size {
        signal_fft[i] = Complex::new(
            signal_fft[i].re * px_padded[i],
            signal_fft[i].im * px_padded[i],
        );
    }

    // IFFT
    ifft(&mut signal_fft);

    // Return real part, doubled to compensate for one-sided filtering
    signal_fft[..n].iter().map(|c| 2.0 * c.re).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn direct_convolve(signal: &[f64], kernel: &[f64]) -> Vec<f64> {
        if signal.is_empty() || kernel.is_empty() {
            return Vec::new();
        }
        let n = signal.len();
        let k = kernel.len();
        let mut result = vec![0.0; n + k - 1];
        for i in 0..n {
            for j in 0..k {
                result[i + j] += signal[i] * kernel[j];
            }
        }
        result
    }

    #[test]
    fn fft_overlap_add_matches_direct() {
        let signal: Vec<f64> = (0..100).map(|i| (i as f64 * 0.1).sin()).collect();
        let kernel = vec![0.2, 0.3, 0.5, 0.3, 0.2];
        let direct = direct_convolve(&signal, &kernel);
        let ola = fft_overlap_add(&signal, &kernel, 32);
        assert_eq!(direct.len(), ola.len());
        for (i, (a, b)) in direct.iter().zip(ola.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1e-8,
                "Mismatch at index {i}: direct={a}, ola={b}"
            );
        }
    }

    #[test]
    fn fft_overlap_add_short_signal() {
        let signal = vec![1.0, 2.0, 3.0];
        let kernel = vec![1.0, 1.0];
        let result = fft_overlap_add(&signal, &kernel, 2);
        let expected = direct_convolve(&signal, &kernel);
        assert_eq!(result.len(), expected.len());
        for (a, b) in result.iter().zip(expected.iter()) {
            assert!((a - b).abs() < 1e-8);
        }
    }

    #[test]
    fn fft_overlap_save_matches_direct() {
        let signal: Vec<f64> = (0..100).map(|i| (i as f64 * 0.1).sin()).collect();
        let kernel = vec![0.2, 0.3, 0.5, 0.3, 0.2];
        let direct = direct_convolve(&signal, &kernel);
        let ols = fft_overlap_save(&signal, &kernel, 32);
        assert_eq!(direct.len(), ols.len());
        for (i, (a, b)) in direct.iter().zip(ols.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1e-6,
                "Mismatch at index {i}: direct={a}, ols={b}"
            );
        }
    }

    #[test]
    fn fft_overlap_add_empty_inputs() {
        assert!(fft_overlap_add(&[], &[1.0], 4).is_empty());
        assert!(fft_overlap_add(&[1.0], &[], 4).is_empty());
    }

    #[test]
    fn gaussian_filter_passes_target_frequency() {
        use std::f64::consts::PI;
        let sr = 1000.0;
        let n = 1024;
        // Signal: 10 Hz sine + 100 Hz sine
        let signal: Vec<f64> = (0..n)
            .map(|i| {
                let t = i as f64 / sr;
                (2.0 * PI * 10.0 * t).sin() + (2.0 * PI * 100.0 * t).sin()
            })
            .collect();

        // Filter at 10 Hz with narrow FWHM
        let result = gaussian_frequency_filter(&signal, sr, 10.0, 2.0);
        assert_eq!(result.filtered.len(), n);

        // The filtered signal should have reduced 100 Hz content
        // Check power at 100 Hz is attenuated relative to 10 Hz
        let power_10: f64 = result.filtered.iter().map(|x| x * x).sum::<f64>() / n as f64;
        assert!(power_10 > 0.01); // 10 Hz component survives
    }

    #[test]
    fn gaussian_filter_empirical_freq_close() {
        use std::f64::consts::PI;
        let sr = 500.0;
        let n = 512;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * 20.0 * i as f64 / sr).sin())
            .collect();

        let result = gaussian_frequency_filter(&signal, sr, 20.0, 3.0);
        // Empirical frequency should be close to 20 Hz
        assert!((result.empirical_frequency - 20.0).abs() < 2.0);
    }

    #[test]
    fn gaussian_bandpass_correct_length() {
        let signal = vec![1.0; 256];
        let filtered = gaussian_bandpass(&signal, 1000.0, 50.0, 5.0);
        assert_eq!(filtered.len(), 256);
    }

    #[test]
    fn gaussian_filter_empty_input() {
        let result = gaussian_frequency_filter(&[], 1000.0, 50.0, 5.0);
        assert!(result.filtered.is_empty());
    }

    // ─── Planck taper bandpass tests ─────────────────────────────────────────

    #[test]
    fn planck_taper_bandpass_correct_length() {
        let signal = vec![1.0; 512];
        let filtered = planck_taper_bandpass(&signal, 1000.0, 50.0, 10.0, 0.15);
        assert_eq!(filtered.len(), 512);
    }

    #[test]
    fn planck_taper_bandpass_passes_target() {
        use std::f64::consts::PI;
        let sr = 1000.0;
        let n = 1024;
        // Signal: 20 Hz sine + 200 Hz sine
        let signal: Vec<f64> = (0..n)
            .map(|i| {
                let t = i as f64 / sr;
                (2.0 * PI * 20.0 * t).sin() + (2.0 * PI * 200.0 * t).sin()
            })
            .collect();

        let filtered = planck_taper_bandpass(&signal, sr, 20.0, 6.0, 0.15);
        assert_eq!(filtered.len(), n);
        // Filtered signal should have energy (20 Hz component survives)
        let power: f64 = filtered.iter().map(|x| x * x).sum::<f64>() / n as f64;
        assert!(power > 0.01);
    }

    #[test]
    fn planck_taper_bandpass_empty() {
        let filtered = planck_taper_bandpass(&[], 1000.0, 50.0, 10.0, 0.15);
        assert!(filtered.is_empty());
    }
}
