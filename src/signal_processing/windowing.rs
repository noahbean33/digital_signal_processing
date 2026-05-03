use std::f64::consts::PI;

/// Apply a Hann (Hanning) window to a signal.
#[must_use]
pub fn hann(size: usize) -> Vec<f64> {
    if size <= 1 {
        return vec![1.0; size];
    }
    (0..size)
        .map(|n| 0.5 * (1.0 - (2.0 * PI * n as f64 / (size - 1) as f64).cos()))
        .collect()
}

/// Apply a Hamming window.
#[must_use]
pub fn hamming(size: usize) -> Vec<f64> {
    if size <= 1 {
        return vec![1.0; size];
    }
    (0..size)
        .map(|n| 0.54 - 0.46 * (2.0 * PI * n as f64 / (size - 1) as f64).cos())
        .collect()
}

/// Apply a Blackman window.
#[must_use]
pub fn blackman(size: usize) -> Vec<f64> {
    if size <= 1 {
        return vec![1.0; size];
    }
    (0..size)
        .map(|n| {
            let x = n as f64 / (size - 1) as f64;
            0.42 - 0.5 * (2.0 * PI * x).cos() + 0.08 * (4.0 * PI * x).cos()
        })
        .collect()
}

/// Apply a rectangular (boxcar) window — all ones.
#[must_use]
pub fn rectangular(size: usize) -> Vec<f64> {
    vec![1.0; size]
}

/// Apply a Bartlett (triangular) window.
#[must_use]
pub fn bartlett(size: usize) -> Vec<f64> {
    if size <= 1 {
        return vec![1.0; size];
    }
    let half = (size - 1) as f64 / 2.0;
    (0..size)
        .map(|n| 1.0 - ((n as f64 - half) / half).abs())
        .collect()
}

/// Apply a Blackman-Harris window.
#[must_use]
pub fn blackman_harris(size: usize) -> Vec<f64> {
    if size <= 1 {
        return vec![1.0; size];
    }
    let a0 = 0.355_768;
    let a1 = 0.487_396;
    let a2 = 0.144_232;
    let a3 = 0.012_604;
    (0..size)
        .map(|n| {
            let x = 2.0 * PI * n as f64 / (size - 1) as f64;
            a0 - a1 * x.cos() + a2 * (2.0 * x).cos() - a3 * (3.0 * x).cos()
        })
        .collect()
}

/// Apply a Flat-Top window (maximum amplitude accuracy).
#[must_use]
pub fn flat_top(size: usize) -> Vec<f64> {
    if size <= 1 {
        return vec![1.0; size];
    }
    let a0 = 0.215_578_95;
    let a1 = 0.416_631_58;
    let a2 = 0.277_263_16;
    let a3 = 0.083_578_947;
    let a4 = 0.006_947_368;
    (0..size)
        .map(|n| {
            let x = 2.0 * PI * n as f64 / (size - 1) as f64;
            a0 - a1 * x.cos() + a2 * (2.0 * x).cos() - a3 * (3.0 * x).cos()
                + a4 * (4.0 * x).cos()
        })
        .collect()
}

/// Apply a Kaiser window with parameter beta.
#[must_use]
pub fn kaiser(size: usize, beta: f64) -> Vec<f64> {
    if size <= 1 {
        return vec![1.0; size];
    }
    let denom = bessel_i0(beta);
    if denom.abs() < 1e-30 {
        return vec![1.0; size];
    }
    let half = (size - 1) as f64 / 2.0;
    (0..size)
        .map(|n| {
            let x = (n as f64 - half) / half;
            bessel_i0(beta * (1.0 - x * x).max(0.0).sqrt()) / denom
        })
        .collect()
}

/// Zeroth-order modified Bessel function of the first kind (series approx).
fn bessel_i0(x: f64) -> f64 {
    let mut sum = 1.0;
    let mut term = 1.0;
    let x_half_sq = (x / 2.0) * (x / 2.0);
    for k in 1..=25 {
        term *= x_half_sq / (k as f64 * k as f64);
        sum += term;
        if term.abs() < 1e-16 * sum.abs() {
            break;
        }
    }
    sum
}

/// Multiply a signal element-wise by a window.
#[must_use]
pub fn apply_window(signal: &[f64], window: &[f64]) -> Vec<f64> {
    signal
        .iter()
        .zip(window.iter())
        .map(|(&s, &w)| s * w)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hann_endpoints_are_zero() {
        let w = hann(64);
        assert!(w[0].abs() < 1e-10);
        assert!(w[63].abs() < 1e-10);
    }

    #[test]
    fn hamming_centre_is_one() {
        let w = hamming(65);
        assert!((w[32] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn rectangular_is_all_ones() {
        let w = rectangular(10);
        assert!(w.iter().all(|&v| (v - 1.0).abs() < 1e-10));
    }

    #[test]
    fn apply_window_scales_correctly() {
        let sig = vec![2.0; 4];
        let win = vec![0.5; 4];
        let out = apply_window(&sig, &win);
        assert!(out.iter().all(|&v| (v - 1.0).abs() < 1e-10));
    }
}
