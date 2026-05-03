/// Apply harmonic distortion (soft clipping) to a signal.
///
/// Samples exceeding `threshold` are compressed using `tanh(gain * x / threshold)`.
#[must_use]
pub fn harmonic_distortion(input: &[f64], threshold: f64, gain: f64) -> Vec<f64> {
    if threshold.abs() < 1e-9 {
        return input.to_vec();
    }
    input
        .iter()
        .map(|&sample| {
            if sample.abs() > threshold {
                threshold * (gain * (sample / threshold)).tanh()
            } else {
                sample
            }
        })
        .collect()
}

/// μ-law style dynamic range compression.
///
/// `compression_ratio` controls the amount of compression.
#[must_use]
pub fn dynamic_range_compression(input: &[f64], compression_ratio: f64) -> Vec<f64> {
    if compression_ratio <= 0.0 {
        return input.to_vec();
    }
    let log_cr = (1.0 + compression_ratio).ln();
    if log_cr.abs() < 1e-9 {
        return input.to_vec();
    }
    input
        .iter()
        .map(|&sample| {
            let abs_s = sample.abs();
            let compressed = (1.0 + compression_ratio * abs_s).ln() / log_cr;
            if sample >= 0.0 {
                compressed
            } else {
                -compressed
            }
        })
        .collect()
}

/// Non-linear FIR filter with sigmoid activation.
///
/// Each output is `sigmoid(steepness * Σ coeffs[k] * input[n-k])`.
#[must_use]
pub fn nonlinear_fir(input: &[f64], coeffs: &[f64], steepness: f64) -> Vec<f64> {
    let filter_size = coeffs.len();
    let mut output = vec![0.0; input.len()];
    for n in 0..input.len() {
        let mut sum = 0.0;
        for k in 0..filter_size {
            if n >= k {
                sum += coeffs[k] * input[n - k];
            }
        }
        output[n] = 1.0 / (1.0 + (-steepness * sum).exp());
    }
    output
}

/// Newton-Raphson method to invert a nonlinear function.
///
/// Finds `x` such that `f(x) ≈ target`.
#[must_use]
pub fn newton_raphson_invert(
    target: f64,
    f: impl Fn(f64) -> f64,
    df: impl Fn(f64) -> f64,
    initial_guess: f64,
    tolerance: f64,
    max_iterations: usize,
) -> f64 {
    let mut x = initial_guess;
    for _ in 0..max_iterations {
        let fx = f(x) - target;
        if fx.abs() < tolerance {
            return x;
        }
        let dfx = df(x);
        if dfx.abs() < 1e-12 {
            break;
        }
        x -= fx / dfx;
    }
    x
}

/// Damped Newton-Raphson solver for `f(x) = 0`.
#[must_use]
pub fn damped_newton_raphson(
    f: impl Fn(f64) -> f64,
    df: impl Fn(f64) -> f64,
    initial: f64,
    damping: f64,
    tolerance: f64,
    max_iterations: usize,
) -> f64 {
    let mut x = initial;
    for _ in 0..max_iterations {
        let val = f(x);
        if val.abs() < tolerance {
            break;
        }
        let derivative = df(x);
        if derivative.abs() < 1e-12 {
            break;
        }
        x -= damping * (val / derivative);
    }
    x
}

/// Adaptive nonlinear weight update using a sigmoid activation.
///
/// Returns the updated weight vector.
#[must_use]
pub fn adaptive_nonlinear_update(
    weights: &[f64],
    input: &[f64],
    desired: f64,
    learning_rate: f64,
    steepness: f64,
) -> Vec<f64> {
    if weights.len() != input.len() {
        return weights.to_vec();
    }
    let sum: f64 = weights.iter().zip(input.iter()).map(|(&w, &x)| w * x).sum();
    let output = 1.0 / (1.0 + (-steepness * sum).exp());
    let error = desired - output;
    let gradient = steepness * output * (1.0 - output);

    weights
        .iter()
        .zip(input.iter())
        .map(|(&w, &x)| w + learning_rate * error * gradient * x)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distortion_leaves_small_signals_unchanged() {
        let input = vec![0.1, 0.2, -0.1];
        let output = harmonic_distortion(&input, 0.5, 2.0);
        assert_eq!(output, input);
    }

    #[test]
    fn compression_preserves_sign() {
        let input = vec![-0.5, 0.5];
        let output = dynamic_range_compression(&input, 2.0);
        assert!(output[0] < 0.0);
        assert!(output[1] > 0.0);
    }

    #[test]
    fn newton_raphson_finds_log3() {
        let result = newton_raphson_invert(
            3.0,
            |x: f64| x.exp(),
            |x: f64| x.exp(),
            1.0,
            1e-6,
            100,
        );
        assert!((result - 3.0_f64.ln()).abs() < 1e-5);
    }

    #[test]
    fn nonlinear_fir_output_bounded() {
        let input = vec![1.0, -1.0, 0.5, -0.5];
        let coeffs = vec![0.5, 0.3];
        let output = nonlinear_fir(&input, &coeffs, 1.0);
        assert!(output.iter().all(|&v| v >= 0.0 && v <= 1.0));
    }
}
