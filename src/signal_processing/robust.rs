/// Huber weight function for robust estimation.
///
/// Returns 1.0 for small residuals, `threshold / |residual|` for large ones.
#[must_use]
pub fn huber_weight(residual: f64, threshold: f64) -> f64 {
    if residual.abs() <= threshold {
        1.0
    } else {
        threshold / residual.abs()
    }
}

/// Compute the median of a slice (non-destructive; clones internally).
#[must_use]
pub fn median(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

/// Compute the Median Absolute Deviation (MAD).
#[must_use]
pub fn mad(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let med = median(data);
    let deviations: Vec<f64> = data.iter().map(|&v| (v - med).abs()).collect();
    median(&deviations)
}

/// Compute a robust confidence interval using MAD.
///
/// Returns `(lower_bound, upper_bound)`.
#[must_use]
pub fn robust_confidence_interval(data: &[f64], confidence_factor: f64) -> (f64, f64) {
    if data.is_empty() {
        return (0.0, 0.0);
    }
    let med = median(data);
    let m = mad(data);
    (med - confidence_factor * m, med + confidence_factor * m)
}

/// Perform one iteration of Iteratively Reweighted Least Squares (IRLS)
/// for a simple linear model y = β·x.
///
/// Returns the updated β coefficient.
#[must_use]
pub fn irls_update(x: &[f64], y: &[f64], beta: f64, threshold: f64) -> f64 {
    if x.len() != y.len() || x.is_empty() {
        return beta;
    }
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for i in 0..x.len() {
        let residual = y[i] - beta * x[i];
        let w = huber_weight(residual, threshold);
        numerator += w * y[i] * x[i];
        denominator += w * x[i] * x[i];
    }
    if denominator > 1e-9 {
        numerator / denominator
    } else {
        beta
    }
}

/// Robust moving median filter.
///
/// `window_size` must be odd and ≥ 3.
#[must_use]
pub fn robust_moving_median(signal: &[f64], window_size: usize) -> Vec<f64> {
    if window_size % 2 == 0 || window_size < 3 || signal.len() < window_size {
        return Vec::new();
    }
    let half = window_size / 2;
    let mut filtered = Vec::with_capacity(signal.len() - window_size + 1);
    for i in half..signal.len() - half {
        let window = &signal[i - half..=i + half];
        filtered.push(median(window));
    }
    filtered
}

/// Robust weighted sum using Huber weights.
#[must_use]
pub fn robust_weighted_sum(residuals: &[f64], threshold: f64) -> f64 {
    residuals
        .iter()
        .map(|&r| huber_weight(r, threshold) * r)
        .sum()
}

/// Full IRLS iteration to convergence.
///
/// Returns the final β estimate.
#[must_use]
pub fn irls(
    x: &[f64],
    y: &[f64],
    initial_beta: f64,
    threshold: f64,
    tolerance: f64,
    max_iterations: usize,
) -> f64 {
    let mut beta = initial_beta;
    for _ in 0..max_iterations {
        let new_beta = irls_update(x, y, beta, threshold);
        if (new_beta - beta).abs() < tolerance {
            return new_beta;
        }
        beta = new_beta;
    }
    beta
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn huber_weight_small_residual() {
        assert!((huber_weight(0.5, 1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn huber_weight_large_residual() {
        let w = huber_weight(2.0, 1.0);
        assert!((w - 0.5).abs() < 1e-10);
    }

    #[test]
    fn median_odd_length() {
        assert!((median(&[3.0, 1.0, 4.0, 1.0, 5.0]) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn median_even_length() {
        assert!((median(&[1.0, 2.0, 3.0, 4.0]) - 2.5).abs() < 1e-10);
    }

    #[test]
    fn mad_constant_is_zero() {
        assert!((mad(&[5.0, 5.0, 5.0]) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn irls_converges() {
        let x: Vec<f64> = (1..=10).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| 2.0 * xi + 0.1).collect();
        let beta = irls(&x, &y, 1.0, 1.5, 1e-6, 100);
        assert!((beta - 2.0).abs() < 0.5);
    }
}
