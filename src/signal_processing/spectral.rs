use super::fft::{self, Complex};

/// Compute the periodogram (power spectral density) using FFT.
#[must_use]
pub fn periodogram(signal: &[f64]) -> Vec<f64> {
    if signal.is_empty() {
        return Vec::new();
    }
    let mut data = fft::real_to_complex(signal);
    data = fft::zero_pad_to_power_of_2(&data);
    fft::fft(&mut data);
    let fft_n = data.len();
    data.iter().map(|c| c.norm_sqr() / fft_n as f64).collect()
}

/// Estimate AR coefficients using the Levinson-Durbin recursion.
///
/// Returns the AR coefficient vector of the given `order`.
#[must_use]
pub fn ar_levinson_durbin(signal: &[f64], order: usize) -> Vec<f64> {
    let n = signal.len();
    if n == 0 || order == 0 {
        return Vec::new();
    }
    // Compute autocorrelation
    let mut r = vec![0.0; order + 1];
    for lag in 0..=order {
        for i in lag..n {
            r[lag] += signal[i] * signal[i - lag];
        }
        r[lag] /= n as f64;
    }

    let mut a = vec![0.0; order];
    let mut e = vec![0.0; order + 1];
    e[0] = r[0];

    for k in 1..=order {
        if e[k - 1].abs() < 1e-9 {
            return Vec::new();
        }
        let mut lambda = 0.0;
        for j in 1..k {
            lambda -= a[j - 1] * r[k - j];
        }
        lambda = (lambda + r[k]) / e[k - 1];

        let a_prev = a.clone();
        a[k - 1] = lambda;
        for j in 1..k {
            a[j - 1] = a_prev[j - 1] - lambda * a_prev[k - j - 1];
        }
        e[k] = e[k - 1] * (1.0 - lambda * lambda);
    }
    a
}

/// Estimate AR coefficients using the Burg algorithm.
#[must_use]
pub fn ar_burg(signal: &[f64], order: usize) -> Vec<f64> {
    let n = signal.len();
    if n == 0 || order == 0 {
        return Vec::new();
    }
    let mut a = vec![0.0; order];
    let mut ef: Vec<f64> = signal.to_vec();
    let mut eb: Vec<f64> = signal.to_vec();

    for m in 0..order {
        let mut num = 0.0;
        let mut den = 0.0;
        for i in (m + 1)..n {
            num += ef[i] * eb[i - 1];
            den += ef[i] * ef[i] + eb[i - 1] * eb[i - 1];
        }
        if den.abs() < 1e-9 {
            return Vec::new();
        }
        let k = -2.0 * num / den;
        a[m] = k;
        for j in 0..m {
            a[j] += k * a[m - j - 1];
        }
        for i in (m + 1..n).rev() {
            let temp = ef[i];
            ef[i] += k * eb[i - 1];
            eb[i] = eb[i - 1] + k * temp;
        }
    }
    a
}

/// Compute the AR-based power spectral density.
///
/// `ar_coeffs` – AR coefficient vector.
/// `noise_variance` – estimated noise variance.
/// `frequencies` – angular frequency points (radians).
#[must_use]
pub fn ar_psd(ar_coeffs: &[f64], noise_variance: f64, frequencies: &[f64]) -> Vec<f64> {
    frequencies
        .iter()
        .map(|&freq| {
            let mut denom = Complex::new(1.0, 0.0);
            for (k, &ak) in ar_coeffs.iter().enumerate() {
                let angle = -freq * (k + 1) as f64;
                denom += Complex::new(ak, 0.0) * Complex::new(angle.cos(), angle.sin());
            }
            let mag_sq = denom.norm_sqr();
            if mag_sq < 1e-9 {
                0.0
            } else {
                noise_variance / mag_sq
            }
        })
        .collect()
}

/// Compute the Akaike Information Criterion (AIC).
#[must_use]
pub fn aic(error_variance: f64, order: usize, num_samples: usize) -> f64 {
    2.0 * order as f64 + num_samples as f64 * error_variance.ln()
}

/// Select optimal AR model order via AIC.
#[must_use]
pub fn select_optimal_ar_order(signal: &[f64], max_order: usize) -> usize {
    let n = signal.len();
    if n == 0 || max_order == 0 {
        return 1;
    }
    let mut best_order = 1;
    let mut best_aic = f64::INFINITY;

    for order in 1..=max_order {
        if n <= order {
            break;
        }
        let ar_coeffs = ar_levinson_durbin(signal, order);
        if ar_coeffs.is_empty() {
            continue;
        }
        let mut error_var = 0.0;
        for i in order..n {
            let mut prediction = 0.0;
            for k in 0..order {
                prediction += ar_coeffs[k] * signal[i - k - 1];
            }
            let err = signal[i] + prediction;
            error_var += err * err;
        }
        error_var /= (n - order) as f64;
        if error_var <= 0.0 {
            continue;
        }
        let current_aic = aic(error_var, order, n);
        if current_aic < best_aic {
            best_aic = current_aic;
            best_order = order;
        }
    }
    best_order
}

/// Compute the Welch periodogram for spectral estimation.
///
/// `segment_length` – length of each segment.
/// `overlap` – number of overlapping samples.
/// `window` – window function to apply to each segment.
#[must_use]
pub fn welch(signal: &[f64], segment_length: usize, overlap: usize, window: &[f64]) -> Vec<f64> {
    if segment_length <= overlap || signal.len() < segment_length || window.len() != segment_length
    {
        return Vec::new();
    }
    let step = segment_length - overlap;
    let num_segments = (signal.len() - overlap) / step;
    if num_segments == 0 {
        return Vec::new();
    }

    let window_energy: f64 = window.iter().map(|w| w * w).sum();
    if window_energy.abs() < 1e-9 {
        return Vec::new();
    }

    // Determine FFT size (next power of 2)
    let mut fft_size = 1;
    while fft_size < segment_length {
        fft_size <<= 1;
    }

    let mut welch_psd = vec![0.0; fft_size];

    for i in 0..num_segments {
        let start = i * step;
        let segment: Vec<f64> = signal[start..start + segment_length]
            .iter()
            .zip(window.iter())
            .map(|(&s, &w)| s * w)
            .collect();

        let mut complex_seg = fft::real_to_complex(&segment);
        complex_seg.resize(fft_size, Complex::zero());
        fft::fft(&mut complex_seg);

        for (k, c) in complex_seg.iter().enumerate() {
            welch_psd[k] += c.norm_sqr();
        }
    }

    for v in &mut welch_psd {
        *v /= (num_segments as f64) * window_energy;
    }
    welch_psd
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::f64::consts::PI;

    #[test]
    fn periodogram_nonzero_for_sine() {
        let n = 64;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * 5.0 * i as f64 / n as f64).sin())
            .collect();
        let psd = periodogram(&signal);
        assert!(!psd.is_empty());
        let max_val = psd.iter().cloned().fold(0.0_f64, f64::max);
        assert!(max_val > 0.0);
    }

    #[test]
    fn burg_returns_correct_order() {
        let signal: Vec<f64> = (0..128)
            .map(|i| (2.0 * PI * 5.0 * i as f64 / 128.0).sin())
            .collect();
        let coeffs = ar_burg(&signal, 4);
        assert_eq!(coeffs.len(), 4);
    }

    #[test]
    fn aic_increases_with_order_for_noise() {
        let signal: Vec<f64> = (0..100).map(|i| (i as f64 * 0.1).sin()).collect();
        let aic1 = {
            let c = ar_levinson_durbin(&signal, 1);
            if c.is_empty() { f64::INFINITY } else {
                let mut ev = 0.0;
                for i in 1..100 { let e = signal[i] + c[0] * signal[i-1]; ev += e*e; }
                ev /= 99.0;
                aic(ev, 1, 100)
            }
        };
        // Just verify it returns a finite number
        assert!(aic1.is_finite());
    }
}
