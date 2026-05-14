/// Single-level Discrete Wavelet Transform.
///
/// Returns `(approximation_coefficients, detail_coefficients)`.
#[must_use]
pub fn dwt_single_level(
    signal: &[f64],
    low_filter: &[f64],
    high_filter: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    if low_filter.is_empty()
        || high_filter.is_empty()
        || low_filter.len() != high_filter.len()
        || signal.is_empty()
    {
        return (Vec::new(), Vec::new());
    }
    let filter_len = low_filter.len();
    let mut approx = Vec::new();
    let mut detail = Vec::new();

    let mut i = 0;
    while i + filter_len <= signal.len() {
        let mut a = 0.0;
        let mut d = 0.0;
        for k in 0..filter_len {
            a += signal[i + k] * low_filter[k];
            d += signal[i + k] * high_filter[k];
        }
        approx.push(a);
        detail.push(d);
        i += 2;
    }
    (approx, detail)
}

/// Single-level Inverse Discrete Wavelet Transform.
#[must_use]
pub fn idwt_single_level(
    approx: &[f64],
    detail: &[f64],
    rec_low_filter: &[f64],
    rec_high_filter: &[f64],
) -> Vec<f64> {
    if approx.is_empty() || detail.is_empty() || rec_low_filter.is_empty() || rec_high_filter.is_empty() {
        return Vec::new();
    }
    let len = approx.len();
    let filter_len = rec_low_filter.len().max(rec_high_filter.len());
    let out_size = 2 * len + filter_len - 2;
    let mut reconstructed = vec![0.0; out_size];

    for i in 0..len {
        let index = 2 * i;
        for k in 0..rec_low_filter.len() {
            if index + k < reconstructed.len() {
                reconstructed[index + k] += approx[i] * rec_low_filter[k];
            }
        }
        for k in 0..rec_high_filter.len() {
            if index + k < reconstructed.len() {
                reconstructed[index + k] += detail[i] * rec_high_filter[k];
            }
        }
    }
    reconstructed
}

/// Multi-level DWT decomposition.
///
/// Returns a vector of `(approximation, detail)` pairs, one per level.
#[must_use]
pub fn dwt_multi_level(
    signal: &[f64],
    low_filter: &[f64],
    high_filter: &[f64],
    levels: usize,
) -> Vec<(Vec<f64>, Vec<f64>)> {
    if levels == 0 {
        return Vec::new();
    }
    let mut decomposition = Vec::with_capacity(levels);
    let mut current = signal.to_vec();
    for _ in 0..levels {
        let (a, d) = dwt_single_level(&current, low_filter, high_filter);
        current = a.clone();
        decomposition.push((a, d));
    }
    decomposition
}

/// Multi-level inverse DWT reconstruction.
#[must_use]
pub fn idwt_multi_level(
    decomposition: &[(Vec<f64>, Vec<f64>)],
    rec_low_filter: &[f64],
    rec_high_filter: &[f64],
) -> Vec<f64> {
    if decomposition.is_empty() {
        return Vec::new();
    }
    let last = decomposition.len() - 1;
    let mut current = idwt_single_level(
        &decomposition[last].0,
        &decomposition[last].1,
        rec_low_filter,
        rec_high_filter,
    );
    for i in (0..last).rev() {
        current = idwt_single_level(&current, &decomposition[i].1, rec_low_filter, rec_high_filter);
    }
    current
}

/// Soft thresholding of wavelet coefficients for denoising.
#[must_use]
pub fn soft_threshold(coeffs: &[f64], threshold: f64) -> Vec<f64> {
    coeffs
        .iter()
        .map(|&c| {
            let abs_val = c.abs();
            if abs_val < threshold {
                0.0
            } else {
                c.signum() * (abs_val - threshold)
            }
        })
        .collect()
}

/// Hard thresholding of wavelet coefficients.
#[must_use]
pub fn hard_threshold(coeffs: &[f64], threshold: f64) -> Vec<f64> {
    coeffs
        .iter()
        .map(|&c| if c.abs() < threshold { 0.0 } else { c })
        .collect()
}

/// Kahan-compensated convolution for improved numerical accuracy.
#[must_use]
pub fn kahan_convolve(signal: &[f64], filter: &[f64], start_index: usize) -> f64 {
    if filter.is_empty() {
        return 0.0;
    }
    let mut sum = 0.0;
    let mut compensation = 0.0;
    for (i, &f) in filter.iter().enumerate() {
        if start_index + i < signal.len() {
            let product = signal[start_index + i] * f;
            let y = product - compensation;
            let t = sum + y;
            compensation = (t - sum) - y;
            sum = t;
        }
    }
    sum
}

// ─── Wavelet Packet functions ───────────────────────────────────────────────

/// Single-level wavelet packet decomposition.
///
/// Returns `[approximation, detail]` coefficient vectors.
#[must_use]
pub fn wavelet_packet_decompose(
    signal: &[f64],
    low_pass_filter: &[f64],
    high_pass_filter: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    if low_pass_filter.len() != high_pass_filter.len() || low_pass_filter.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let filter_size = low_pass_filter.len();
    let length = signal.len();
    let num_coeffs = length / 2;
    let mut approx = vec![0.0; num_coeffs];
    let mut detail = vec![0.0; num_coeffs];

    for i in 0..num_coeffs {
        let mut a_val = 0.0;
        let mut d_val = 0.0;
        for j in 0..filter_size {
            let idx = (2 * i + j) % length;
            a_val += signal[idx] * low_pass_filter[j];
            d_val += signal[idx] * high_pass_filter[j];
        }
        approx[i] = a_val;
        detail[i] = d_val;
    }
    (approx, detail)
}

/// Compute energy of wavelet packet coefficients.
#[must_use]
pub fn packet_energy(coefficients: &[f64]) -> f64 {
    coefficients.iter().map(|c| c * c).sum()
}

/// Single-level wavelet packet reconstruction.
#[must_use]
pub fn wavelet_packet_reconstruct(
    approx: &[f64],
    detail: &[f64],
    synthesis_low_filter: &[f64],
    synthesis_high_filter: &[f64],
) -> Vec<f64> {
    if synthesis_low_filter.len() != synthesis_high_filter.len() || synthesis_low_filter.is_empty() {
        return Vec::new();
    }
    let num_coeffs = approx.len();
    let signal_length = num_coeffs * 2;
    let filter_size = synthesis_low_filter.len();
    let mut signal = vec![0.0; signal_length];

    for i in 0..num_coeffs {
        for j in 0..filter_size {
            let idx = (2 * i + j) % signal_length;
            signal[idx] += approx[i] * synthesis_low_filter[j]
                + detail[i] * synthesis_high_filter[j];
        }
    }
    signal
}

/// Recursive wavelet packet decomposition to obtain terminal (leaf) nodes.
#[must_use]
pub fn recursive_packet_decomposition(
    signal: &[f64],
    low_pass_filter: &[f64],
    high_pass_filter: &[f64],
    max_level: usize,
) -> Vec<Vec<f64>> {
    if max_level == 0 {
        return vec![signal.to_vec()];
    }
    let (approx, detail) = wavelet_packet_decompose(signal, low_pass_filter, high_pass_filter);
    let mut result = recursive_packet_decomposition(&approx, low_pass_filter, high_pass_filter, max_level - 1);
    let right = recursive_packet_decomposition(&detail, low_pass_filter, high_pass_filter, max_level - 1);
    result.extend(right);
    result
}

/// Compute RMSE between original and reconstructed signals.
#[must_use]
pub fn reconstruction_rmse(original: &[f64], reconstructed: &[f64]) -> f64 {
    if original.len() != reconstructed.len() || original.is_empty() {
        return 0.0;
    }
    let sq_error: f64 = original
        .iter()
        .zip(reconstructed.iter())
        .map(|(&o, &r)| (o - r) * (o - r))
        .sum();
    (sq_error / original.len() as f64).sqrt()
}

/// Haar wavelet filter coefficients.
pub const HAAR_LOW: [f64; 2] = [0.707_106_781_186_547_5, 0.707_106_781_186_547_5];
pub const HAAR_HIGH: [f64; 2] = [-0.707_106_781_186_547_5, 0.707_106_781_186_547_5];
pub const HAAR_REC_LOW: [f64; 2] = [0.707_106_781_186_547_5, 0.707_106_781_186_547_5];
pub const HAAR_REC_HIGH: [f64; 2] = [0.707_106_781_186_547_5, -0.707_106_781_186_547_5];

// ─── Wavelet Generation ──────────────────────────────────────────────────────

/// Generate a real-valued Morlet wavelet.
///
/// `morlet[t] = cos(2π·freq·t) · exp(−4·ln2·t² / fwhm²)`
///
/// * `freq`  – peak frequency in Hz.
/// * `fwhm`  – full-width at half-maximum of the Gaussian envelope **in seconds**.
/// * `sample_rate` – sampling rate in Hz.
/// * `half_len` – number of samples on each side of the centre
///   (total length = `2 * half_len + 1`).
///
/// Returns a vector of wavelet coefficients centred at index `half_len`.
#[must_use]
pub fn morlet_wavelet(freq: f64, fwhm: f64, sample_rate: f64, half_len: usize) -> Vec<f64> {
    use std::f64::consts::PI;
    let len = 2 * half_len + 1;
    let mut wav = Vec::with_capacity(len);
    for i in 0..len {
        let t = (i as f64 - half_len as f64) / sample_rate;
        let cosine = (2.0 * PI * freq * t).cos();
        let gaussian = (-4.0 * 2.0_f64.ln() * t * t / (fwhm * fwhm)).exp();
        wav.push(cosine * gaussian);
    }
    wav
}

/// Generate a complex Morlet wavelet as `(real_part, imaginary_part)` vectors.
///
/// `cmw[t] = exp(j·2π·freq·t) · exp(−4·ln2·t² / fwhm²)`
///
/// Useful for time-frequency analysis where both amplitude and phase are needed.
#[must_use]
pub fn complex_morlet_wavelet(
    freq: f64,
    fwhm: f64,
    sample_rate: f64,
    half_len: usize,
) -> (Vec<f64>, Vec<f64>) {
    use std::f64::consts::PI;
    let len = 2 * half_len + 1;
    let mut re = Vec::with_capacity(len);
    let mut im = Vec::with_capacity(len);
    for i in 0..len {
        let t = (i as f64 - half_len as f64) / sample_rate;
        let gaussian = (-4.0 * 2.0_f64.ln() * t * t / (fwhm * fwhm)).exp();
        let phase = 2.0 * PI * freq * t;
        re.push(phase.cos() * gaussian);
        im.push(phase.sin() * gaussian);
    }
    (re, im)
}

/// Generate a Haar wavelet of a given length.
///
/// The first half is +1, the second half is −1, normalised by `1/√N`.
///
/// * `length` – total number of samples (should be even for a symmetric wavelet).
#[must_use]
pub fn haar_wavelet(length: usize) -> Vec<f64> {
    if length == 0 {
        return Vec::new();
    }
    let half = length / 2;
    let norm = 1.0 / (length as f64).sqrt();
    let mut wav = Vec::with_capacity(length);
    for i in 0..length {
        if i < half {
            wav.push(norm);
        } else {
            wav.push(-norm);
        }
    }
    wav
}

/// Generate a Mexican hat (Ricker) wavelet.
///
/// `ψ(t) = (2 / (√3σ · π^¼)) · (1 − t²/σ²) · exp(−t² / (2σ²))`
///
/// * `sigma` – width parameter controlling the wavelet's time spread.
/// * `sample_rate` – sampling rate in Hz.
/// * `half_len` – number of samples on each side of the centre.
#[must_use]
pub fn mexican_hat_wavelet(sigma: f64, sample_rate: f64, half_len: usize) -> Vec<f64> {
    use std::f64::consts::PI;
    let len = 2 * half_len + 1;
    let norm = 2.0 / ((3.0 * sigma).sqrt() * PI.powf(0.25));
    let mut wav = Vec::with_capacity(len);
    for i in 0..len {
        let t = (i as f64 - half_len as f64) / sample_rate;
        let t_over_s = t / sigma;
        wav.push(norm * (1.0 - t_over_s * t_over_s) * (-t * t / (2.0 * sigma * sigma)).exp());
    }
    wav
}

/// Generate a Difference of Gaussians (DoG) wavelet.
///
/// Approximates the Laplacian of Gaussian by subtracting a broad Gaussian
/// from a narrow one.
///
/// * `sigma_pos` – standard deviation of the narrow (positive) Gaussian.
/// * `sigma_neg` – standard deviation of the broad (negative) Gaussian.
/// * `sample_rate` – sampling rate in Hz.
/// * `half_len` – number of samples on each side of the centre.
#[must_use]
pub fn dog_wavelet(
    sigma_pos: f64,
    sigma_neg: f64,
    sample_rate: f64,
    half_len: usize,
) -> Vec<f64> {
    use std::f64::consts::PI;
    let len = 2 * half_len + 1;
    let two_pi_sqrt = (2.0 * PI).sqrt();
    let mut wav = Vec::with_capacity(len);
    for i in 0..len {
        let t = (i as f64 - half_len as f64) / sample_rate;
        let g1 = (-t * t / (2.0 * sigma_pos * sigma_pos)).exp() / (sigma_pos * two_pi_sqrt);
        let g2 = (-t * t / (2.0 * sigma_neg * sigma_neg)).exp() / (sigma_neg * two_pi_sqrt);
        wav.push(g1 - g2);
    }
    wav
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dwt_single_level_basic() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let (a, d) = dwt_single_level(&signal, &HAAR_LOW, &HAAR_HIGH);
        assert_eq!(a.len(), 4);
        assert_eq!(d.len(), 4);
    }

    #[test]
    fn soft_threshold_zeros_small() {
        let coeffs = vec![0.1, 0.5, -0.3, 2.0, -1.5];
        let thresholded = soft_threshold(&coeffs, 0.5);
        assert!((thresholded[0]).abs() < 1e-10);
        assert!((thresholded[1]).abs() < 1e-10);
        assert!((thresholded[2]).abs() < 1e-10);
        assert!((thresholded[3] - 1.5).abs() < 1e-10);
        assert!((thresholded[4] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn packet_energy_correct() {
        let coeffs = vec![1.0, 2.0, 3.0];
        assert!((packet_energy(&coeffs) - 14.0).abs() < 1e-10);
    }

    // ─── Wavelet generation tests ────────────────────────────────────────────

    #[test]
    fn morlet_wavelet_length_and_peak() {
        let wav = morlet_wavelet(10.0, 0.5, 1000.0, 500);
        assert_eq!(wav.len(), 1001);
        // Peak should be at the centre (t=0 → cos(0)=1, gaussian=1)
        let center = 500;
        let peak_idx = wav
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;
        assert_eq!(peak_idx, center);
    }

    #[test]
    fn morlet_wavelet_decays_to_zero() {
        let wav = morlet_wavelet(10.0, 0.2, 1000.0, 500);
        // Far from centre, values should be near zero
        assert!(wav[0].abs() < 1e-6);
        assert!(wav[1000].abs() < 1e-6);
    }

    #[test]
    fn complex_morlet_wavelet_real_imag_orthogonal() {
        let (re, im) = complex_morlet_wavelet(10.0, 0.5, 1000.0, 500);
        assert_eq!(re.len(), 1001);
        assert_eq!(im.len(), 1001);
        // At t=0: cos(0)=1 (real peak), sin(0)=0 (imag zero)
        assert!((re[500] - 1.0).abs() < 1e-10);
        assert!(im[500].abs() < 1e-10);
    }

    #[test]
    fn haar_wavelet_structure() {
        let wav = haar_wavelet(100);
        assert_eq!(wav.len(), 100);
        let norm = 1.0 / 10.0; // 1/sqrt(100)
        // First half positive, second half negative
        assert!((wav[0] - norm).abs() < 1e-10);
        assert!((wav[49] - norm).abs() < 1e-10);
        assert!((wav[50] + norm).abs() < 1e-10);
        assert!((wav[99] + norm).abs() < 1e-10);
    }

    #[test]
    fn haar_wavelet_sums_to_zero() {
        let wav = haar_wavelet(64);
        let sum: f64 = wav.iter().sum();
        assert!(sum.abs() < 1e-10);
    }

    #[test]
    fn haar_wavelet_empty() {
        assert!(haar_wavelet(0).is_empty());
    }

    #[test]
    fn mexican_hat_wavelet_center_positive() {
        let wav = mexican_hat_wavelet(0.4, 1000.0, 500);
        assert_eq!(wav.len(), 1001);
        // Centre (t=0) should be the maximum positive value
        let peak_idx = wav
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;
        assert_eq!(peak_idx, 500);
        assert!(wav[500] > 0.0);
    }

    #[test]
    fn mexican_hat_wavelet_has_negative_lobes() {
        let wav = mexican_hat_wavelet(0.4, 1000.0, 500);
        // Should have negative values away from centre
        let has_negative = wav.iter().any(|&v| v < -1e-10);
        assert!(has_negative);
    }

    #[test]
    fn dog_wavelet_center_positive() {
        let wav = dog_wavelet(0.1, 0.5, 1000.0, 500);
        assert_eq!(wav.len(), 1001);
        // Narrow Gaussian peaks higher than broad at centre → positive peak
        assert!(wav[500] > 0.0);
    }

    #[test]
    fn dog_wavelet_has_negative_values() {
        let wav = dog_wavelet(0.1, 0.5, 1000.0, 500);
        // Broad Gaussian dominates far from centre → negative values exist
        let has_negative = wav.iter().any(|&v| v < -1e-10);
        assert!(has_negative);
    }

    #[test]
    fn dog_wavelet_approaches_zero_at_edges() {
        let wav = dog_wavelet(0.1, 0.5, 1000.0, 4000);
        assert!(wav[0].abs() < 1e-6);
        assert!(wav[8000].abs() < 1e-6);
    }
}
