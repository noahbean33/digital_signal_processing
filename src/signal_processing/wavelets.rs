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
}
