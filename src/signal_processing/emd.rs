/// Compute the mean envelope of a signal using linear interpolation between
/// local maxima and minima envelopes.
///
/// In a production implementation this would use cubic spline interpolation
/// through detected extrema. Here we use a simplified version that computes
/// local mean via a sliding window as a placeholder for the full sifting
/// process (matching the C++ library's `computeMeanEnvelope` placeholder).
#[must_use]
pub fn compute_mean_envelope(signal: &[f64]) -> Vec<f64> {
    if signal.len() < 3 {
        return vec![0.0; signal.len()];
    }

    let maxima = find_local_maxima(signal);
    let minima = find_local_minima(signal);

    if maxima.len() < 2 || minima.len() < 2 {
        return vec![0.0; signal.len()];
    }

    let upper = interpolate_envelope(signal, &maxima);
    let lower = interpolate_envelope(signal, &minima);

    upper
        .iter()
        .zip(lower.iter())
        .map(|(&u, &l)| (u + l) / 2.0)
        .collect()
}

/// Perform a single sifting iteration: signal − mean_envelope.
#[must_use]
pub fn sift_once(signal: &[f64]) -> Vec<f64> {
    let mean_envelope = compute_mean_envelope(signal);
    signal
        .iter()
        .zip(mean_envelope.iter())
        .map(|(&s, &m)| s - m)
        .collect()
}

/// Compute the average absolute difference between two iterations
/// (convergence criterion).
#[must_use]
pub fn convergence_error(prev: &[f64], curr: &[f64]) -> f64 {
    if prev.len() != curr.len() || prev.is_empty() {
        return 0.0;
    }
    let error: f64 = prev
        .iter()
        .zip(curr.iter())
        .map(|(&p, &c)| (c - p).abs())
        .sum();
    error / prev.len() as f64
}

/// Check whether a candidate signal qualifies as an Intrinsic Mode Function.
///
/// An IMF must satisfy:
/// 1. The number of extrema and zero crossings differ by at most 1.
/// 2. The mean envelope is approximately zero everywhere.
#[must_use]
pub fn is_imf(signal: &[f64]) -> bool {
    if signal.len() < 2 {
        return false;
    }
    let zero_crossings = count_zero_crossings(signal);
    let extrema = count_extrema(signal);
    (extrema as isize - zero_crossings as isize).unsigned_abs() <= 1
}

/// Adaptive sifting to extract a single IMF.
#[must_use]
pub fn adaptive_sift(signal: &[f64], max_iter: usize, tolerance: f64) -> Vec<f64> {
    let mut current = signal.to_vec();
    for _ in 0..max_iter {
        let previous = current.clone();
        current = sift_once(&current);
        let error = convergence_error(&previous, &current);
        if error < tolerance && is_imf(&current) {
            break;
        }
    }
    current
}

/// Full EMD decomposition: extract all IMFs from a signal.
///
/// Returns a vector of IMFs. The last element is the residual.
#[must_use]
pub fn emd(signal: &[f64], max_imfs: usize, max_sift_iter: usize, tolerance: f64) -> Vec<Vec<f64>> {
    let mut imfs = Vec::new();
    let mut residual = signal.to_vec();

    for _ in 0..max_imfs {
        if residual.len() < 3 {
            break;
        }
        let imf = adaptive_sift(&residual, max_sift_iter, tolerance);

        // Check if the IMF is essentially zero (decomposition complete)
        let energy: f64 = imf.iter().map(|x| x * x).sum();
        if energy < tolerance * tolerance * residual.len() as f64 {
            break;
        }

        residual = residual
            .iter()
            .zip(imf.iter())
            .map(|(&r, &i)| r - i)
            .collect();
        imfs.push(imf);
    }
    imfs.push(residual);
    imfs
}

/// Cubic spline interpolation through given knot points.
///
/// `x` and `y` are the knot coordinates; `query` are the x-positions to evaluate.
#[must_use]
pub fn cubic_spline_interpolate(x: &[f64], y: &[f64], query: &[f64]) -> Vec<f64> {
    let n = x.len();
    if n < 2 || x.len() != y.len() {
        return vec![0.0; query.len()];
    }

    let mut h = vec![0.0; n - 1];
    for i in 0..n - 1 {
        h[i] = x[i + 1] - x[i];
    }

    let mut alpha = vec![0.0; n];
    for i in 1..n - 1 {
        if h[i].abs() < 1e-30 || h[i - 1].abs() < 1e-30 {
            continue;
        }
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
    }

    let mut l = vec![0.0; n];
    let mut mu = vec![0.0; n];
    let mut z = vec![0.0; n];
    let mut c = vec![0.0; n];
    let mut b = vec![0.0; n - 1];
    let mut d = vec![0.0; n - 1];

    l[0] = 1.0;
    for i in 1..n - 1 {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        if l[i].abs() < 1e-9 {
            l[i] = 1.0;
        }
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    l[n - 1] = 1.0;

    for i in (0..n - 1).rev() {
        c[i] = z[i] - mu[i] * c[i + 1];
        if h[i].abs() < 1e-30 {
            continue;
        }
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
    }

    query
        .iter()
        .map(|&q| {
            let mut j = n - 2;
            for i in 0..n - 1 {
                if q >= x[i] && q <= x[i + 1] {
                    j = i;
                    break;
                }
            }
            let diff = q - x[j];
            y[j] + b[j] * diff + c[j] * diff * diff + d[j] * diff * diff * diff
        })
        .collect()
}

/// Dynamically adjust convergence tolerance.
#[must_use]
pub fn adjust_tolerance(current_error: f64, previous_error: f64, base_tolerance: f64) -> f64 {
    if current_error < previous_error {
        base_tolerance * 0.95
    } else {
        base_tolerance * 1.05
    }
}

// ─── Helpers ────────────────────────────────────────────────────────────────

fn find_local_maxima(signal: &[f64]) -> Vec<usize> {
    let mut maxima = Vec::new();
    for i in 1..signal.len() - 1 {
        if signal[i] > signal[i - 1] && signal[i] > signal[i + 1] {
            maxima.push(i);
        }
    }
    maxima
}

fn find_local_minima(signal: &[f64]) -> Vec<usize> {
    let mut minima = Vec::new();
    for i in 1..signal.len() - 1 {
        if signal[i] < signal[i - 1] && signal[i] < signal[i + 1] {
            minima.push(i);
        }
    }
    minima
}

fn interpolate_envelope(signal: &[f64], indices: &[usize]) -> Vec<f64> {
    if indices.len() < 2 {
        return vec![0.0; signal.len()];
    }
    let x: Vec<f64> = indices.iter().map(|&i| i as f64).collect();
    let y: Vec<f64> = indices.iter().map(|&i| signal[i]).collect();
    let query: Vec<f64> = (0..signal.len()).map(|i| i as f64).collect();
    cubic_spline_interpolate(&x, &y, &query)
}

fn count_zero_crossings(signal: &[f64]) -> usize {
    let mut count = 0;
    for i in 1..signal.len() {
        if (signal[i - 1] <= 0.0 && signal[i] > 0.0)
            || (signal[i - 1] >= 0.0 && signal[i] < 0.0)
        {
            count += 1;
        }
    }
    count
}

fn count_extrema(signal: &[f64]) -> usize {
    if signal.len() < 3 {
        return 0;
    }
    let mut count = 0;
    for i in 1..signal.len() - 1 {
        if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1])
            || (signal[i] < signal[i - 1] && signal[i] < signal[i + 1])
        {
            count += 1;
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn sift_once_produces_same_length() {
        let signal: Vec<f64> = (0..64)
            .map(|i| (2.0 * PI * 5.0 * i as f64 / 64.0).sin())
            .collect();
        let sifted = sift_once(&signal);
        assert_eq!(sifted.len(), signal.len());
    }

    #[test]
    fn convergence_error_of_identical_is_zero() {
        let a = vec![1.0, 2.0, 3.0];
        assert!(convergence_error(&a, &a).abs() < 1e-10);
    }

    #[test]
    fn cubic_spline_passes_through_knots() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![0.0, 1.0, 0.0, 1.0];
        let result = cubic_spline_interpolate(&x, &y, &x);
        for (i, &r) in result.iter().enumerate() {
            assert!((r - y[i]).abs() < 1e-6, "knot {i}: expected {}, got {r}", y[i]);
        }
    }
}
