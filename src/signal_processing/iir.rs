use std::f64::consts::PI;

use super::fft::Complex;

/// Process a single sample using Direct Form I structure.
///
/// `b` – feedforward (numerator) coefficients.
/// `a` – feedback (denominator) coefficients (a[0] is the normalisation term).
/// `x` – input buffer (all samples up to index `n`).
/// `y` – output buffer (all samples up to index `n - 1`).
/// `n` – current sample index.
#[must_use]
pub fn process_sample_df1(b: &[f64], a: &[f64], x: &[f64], y: &[f64], n: usize) -> f64 {
    if a.is_empty() || a[0].abs() < 1e-9 {
        return 0.0;
    }
    let mut acc = 0.0;
    for (i, &bi) in b.iter().enumerate() {
        if n >= i {
            acc += bi * x[n - i];
        }
    }
    for (i, &ai) in a.iter().enumerate().skip(1) {
        if n >= i {
            acc -= ai * y[n - i];
        }
    }
    acc / a[0]
}

/// Apply IIR filtering using a Direct Form II transposed structure.
///
/// `b` and `a` must have the same length. `a[0]` is assumed normalised to 1.
#[must_use]
pub fn apply_df2(input: &[f64], a: &[f64], b: &[f64]) -> Vec<f64> {
    if input.is_empty() || a.is_empty() || b.is_empty() || a.len() != b.len() {
        return Vec::new();
    }
    let order = a.len() - 1;
    let mut output = vec![0.0; input.len()];
    let mut w = vec![0.0; order];

    for (n, &xn) in input.iter().enumerate() {
        let mut w0 = xn;
        for i in 1..=order {
            w0 -= a[i] * w[i - 1];
        }
        let mut yn = b[0] * w0;
        for i in 1..=order {
            yn += b[i] * w[i - 1];
        }
        output[n] = yn;
        if order > 0 {
            for i in (1..order).rev() {
                w[i] = w[i - 1];
            }
            w[0] = w0;
        }
    }
    output
}

/// Check filter stability: all poles must lie inside the unit circle.
#[must_use]
pub fn is_stable(poles: &[Complex]) -> bool {
    poles.iter().all(|p| p.norm() < 1.0)
}

/// Compute the frequency response H(e^{jω}) of the IIR filter.
#[must_use]
pub fn frequency_response(
    b: &[f64],
    a: &[f64],
    frequency: f64,
    sampling_rate: f64,
) -> Complex {
    let mut numerator = Complex::zero();
    let mut denominator = Complex::zero();
    for (k, &bk) in b.iter().enumerate() {
        let angle = -2.0 * PI * frequency * k as f64 / sampling_rate;
        numerator += Complex::new(bk, 0.0) * Complex::new(angle.cos(), angle.sin());
    }
    for (k, &ak) in a.iter().enumerate() {
        let angle = -2.0 * PI * frequency * k as f64 / sampling_rate;
        denominator += Complex::new(ak, 0.0) * Complex::new(angle.cos(), angle.sin());
    }
    if denominator.norm() < 1e-9 {
        return Complex::zero();
    }
    numerator / denominator
}

/// Normalise IIR coefficients so that `a[0] == 1`.
#[must_use]
pub fn normalize_coefficients(b: &[f64], a: &[f64]) -> (Vec<f64>, Vec<f64>) {
    if a.is_empty() || a[0].abs() < 1e-9 {
        return (b.to_vec(), a.to_vec());
    }
    let scale = a[0];
    let b_norm: Vec<f64> = b.iter().map(|&v| v / scale).collect();
    let a_norm: Vec<f64> = a.iter().map(|&v| v / scale).collect();
    (b_norm, a_norm)
}

/// Quantise IIR coefficients to fixed-point representation.
#[must_use]
pub fn quantize_coefficients(coeffs: &[f64], q_factor: i32) -> Vec<i32> {
    coeffs
        .iter()
        .map(|&c| (c * f64::from(q_factor)).round() as i32)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn df2_impulse_response() {
        let b = vec![0.2929, 0.5858, 0.2929];
        let a = vec![1.0, 0.0, 0.1716];
        let mut input = vec![0.0; 10];
        input[0] = 1.0;
        let output = apply_df2(&input, &a, &b);
        assert_eq!(output.len(), 10);
        assert!((output[0] - 0.2929).abs() < 1e-4);
    }

    #[test]
    fn stable_poles() {
        let poles = vec![
            Complex::new(0.8, 0.1),
            Complex::new(0.7, -0.3),
        ];
        assert!(is_stable(&poles));
    }

    #[test]
    fn unstable_pole() {
        let poles = vec![Complex::new(1.0, 0.0)];
        assert!(!is_stable(&poles));
    }

    #[test]
    fn normalize_scales_correctly() {
        let b = vec![0.5, 1.0];
        let a = vec![2.0, 0.4];
        let (bn, an) = normalize_coefficients(&b, &a);
        assert!((an[0] - 1.0).abs() < 1e-10);
        assert!((bn[0] - 0.25).abs() < 1e-10);
    }
}
