/// Normalised Least-Mean-Squares (NLMS) weight update.
///
/// Returns the updated weight vector.
#[must_use]
pub fn nlms_update(
    weights: &[f64],
    input: &[f64],
    error: f64,
    mu: f64,
    epsilon: f64,
) -> Vec<f64> {
    if weights.len() != input.len() || weights.is_empty() {
        return weights.to_vec();
    }
    let input_energy: f64 = input.iter().map(|x| x * x).sum::<f64>() + epsilon;
    weights
        .iter()
        .zip(input.iter())
        .map(|(&w, &x)| w + mu * error * x / input_energy)
        .collect()
}

/// Standard LMS weight update.
pub fn lms_update(weights: &mut [f64], step_size: f64, error: f64, input: &[f64]) {
    for (i, w) in weights.iter_mut().enumerate() {
        if i < input.len() {
            *w += step_size * error * input[i];
        }
    }
}

/// LMS weight update with L2 (Tikhonov) regularisation.
pub fn lms_update_regularized(
    weights: &mut [f64],
    step_size: f64,
    error: f64,
    input: &[f64],
    lambda: f64,
) {
    for (i, w) in weights.iter_mut().enumerate() {
        if i < input.len() {
            *w += step_size * (error * input[i] - lambda * *w);
        }
    }
}

/// Dynamic adjustment of step-size based on instantaneous error.
#[must_use]
pub fn adaptive_step_size(current_error: f64, mu: f64, beta: f64) -> f64 {
    mu / (1.0 + beta * current_error.abs())
}

/// Compute the dot product of two equal-length slices.
#[must_use]
pub fn dot_product(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

/// Compute the error signal: e(n) = d(n) − y(n).
#[must_use]
pub fn compute_error(desired: f64, output: f64) -> f64 {
    desired - output
}

/// Compute mean squared error from an error signal.
#[must_use]
pub fn mean_squared_error(errors: &[f64]) -> f64 {
    if errors.is_empty() {
        return 0.0;
    }
    let sum: f64 = errors.iter().map(|e| e * e).sum();
    sum / errors.len() as f64
}

/// L2 norm difference between two weight vectors (convergence measure).
#[must_use]
pub fn weight_difference(prev: &[f64], curr: &[f64]) -> f64 {
    prev.iter()
        .zip(curr.iter())
        .map(|(&a, &b)| (b - a) * (b - a))
        .sum::<f64>()
        .sqrt()
}

/// Clamped weight update to restrict magnitude of change.
#[must_use]
pub fn clamped_weight_update(current_weight: f64, update_value: f64, max_update: f64) -> f64 {
    let applied = update_value.clamp(-max_update, max_update);
    current_weight + applied
}

/// Convergence metric: average absolute difference between consecutive errors.
#[must_use]
pub fn convergence_metric(error_history: &[f64]) -> f64 {
    if error_history.len() < 2 {
        return 0.0;
    }
    let sum: f64 = error_history
        .windows(2)
        .map(|w| (w[1] - w[0]).abs())
        .sum();
    sum / (error_history.len() - 1) as f64
}

/// RLS gain vector computation.
///
/// `p` – inverse correlation matrix (flattened row-major, n×n).
/// `x` – input vector (length n).
/// `lambda` – forgetting factor.
///
/// Returns the gain vector k of length n.
#[must_use]
pub fn rls_gain(p: &[f64], x: &[f64], lambda: f64, n: usize) -> Vec<f64> {
    let px = mat_vec_mul(p, x, n);
    let denom = lambda + dot_product(x, &px);
    if denom.abs() < 1e-9 {
        return vec![0.0; n];
    }
    px.iter().map(|&v| v / denom).collect()
}

/// Update the RLS inverse correlation matrix.
///
/// Returns the new P matrix (flattened row-major, n×n).
#[must_use]
pub fn rls_update_inverse_correlation(
    p: &[f64],
    x: &[f64],
    lambda: f64,
    n: usize,
) -> Vec<f64> {
    let px = mat_vec_mul(p, x, n);
    let denom = lambda + dot_product(x, &px);
    if denom.abs() < 1e-9 {
        return p.to_vec();
    }
    let mut new_p = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            new_p[i * n + j] = (p[i * n + j] - px[i] * px[j] / denom) / lambda;
        }
    }
    new_p
}

/// Run a full adaptive filter (LMS with regularisation and dynamic step size).
///
/// Returns `(final_weights, output_signal, error_signal)`.
pub fn adaptive_filter_process(
    input_signal: &[f64],
    desired_signal: &[f64],
    initial_weights: &[f64],
    initial_step_size: f64,
    beta: f64,
    lambda: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let filter_order = initial_weights.len();
    let signal_length = input_signal.len();
    let mut weights = initial_weights.to_vec();
    let mut output_signal = vec![0.0; signal_length];
    let mut error_signal = vec![0.0; signal_length];

    for n in (filter_order - 1)..signal_length {
        let input_segment: Vec<f64> = (0..filter_order)
            .map(|k| input_signal[n - k])
            .collect();

        let y = dot_product(&weights, &input_segment);
        output_signal[n] = y;

        let error = compute_error(desired_signal[n], y);
        error_signal[n] = error;

        let step_size = adaptive_step_size(error, initial_step_size, beta);
        lms_update_regularized(&mut weights, step_size, error, &input_segment, lambda);
    }

    (weights, output_signal, error_signal)
}

/// Matrix × vector multiplication for a row-major n×n matrix.
fn mat_vec_mul(mat: &[f64], vec: &[f64], n: usize) -> Vec<f64> {
    let mut result = vec![0.0; n];
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..n {
            s += mat[i * n + j] * vec[j];
        }
        result[i] = s;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nlms_zero_error_no_change() {
        let w = vec![1.0, 2.0, 3.0];
        let x = vec![0.5, 0.5, 0.5];
        let updated = nlms_update(&w, &x, 0.0, 0.05, 1e-6);
        for (a, b) in w.iter().zip(updated.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn mse_of_zero_error_is_zero() {
        let errors = vec![0.0, 0.0, 0.0];
        assert!(mean_squared_error(&errors).abs() < 1e-10);
    }

    #[test]
    fn adaptive_step_size_decreases_with_error() {
        let s1 = adaptive_step_size(0.1, 0.05, 0.1);
        let s2 = adaptive_step_size(10.0, 0.05, 0.1);
        assert!(s1 > s2);
    }
}
