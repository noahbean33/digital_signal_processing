/// A-weighting filter magnitude response per IEC 61672.
///
/// `f` – frequency in Hz.
///
/// Returns the A-weighting gain (linear, not dB) at the given frequency.
#[must_use]
pub fn a_weighting(f: f64) -> f64 {
    let f2 = f * f;
    let num = 12200.0_f64.powi(2) * f2 * f2;
    let den = (f2 + 20.6_f64.powi(2))
        * ((f2 + 107.7_f64.powi(2)) * (f2 + 737.9_f64.powi(2))).sqrt()
        * (f2 + 12200.0_f64.powi(2));
    if den.abs() < 1e-30 {
        return 0.0;
    }
    let unnorm = num / den;
    // Normalisation constant so that A(1000 Hz) ≈ 1
    const GAIN: f64 = 1.258_896_629_083_327_7;
    unnorm * GAIN
}

/// B-weighting filter magnitude response.
///
/// `f` – frequency in Hz.
#[must_use]
pub fn b_weighting(f: f64) -> f64 {
    let f2 = f * f;
    let num = 12200.0_f64.powi(2) * f2 * f.abs();
    let den = (f2 + 20.6_f64.powi(2))
        * (f2 + 158.5_f64.powi(2)).sqrt()
        * (f2 + 12200.0_f64.powi(2));
    if den.abs() < 1e-30 {
        return 0.0;
    }
    let unnorm = num / den;
    const GAIN: f64 = 1.019_718_247_837_232_6;
    unnorm * GAIN
}

/// C-weighting filter magnitude response per IEC 61672.
///
/// `f` – frequency in Hz.
#[must_use]
pub fn c_weighting(f: f64) -> f64 {
    let f2 = f * f;
    let num = 12200.0_f64.powi(2) * f2;
    let den = (f2 + 20.6_f64.powi(2)) * (f2 + 12200.0_f64.powi(2));
    if den.abs() < 1e-30 {
        return 0.0;
    }
    let unnorm = num / den;
    const GAIN: f64 = 1.007_145_835_141_091_1;
    unnorm * GAIN
}

/// Convert a linear weighting gain to dB.
#[must_use]
pub fn weighting_db(linear: f64) -> f64 {
    if linear.abs() < 1e-30 {
        return f64::NEG_INFINITY;
    }
    20.0 * linear.log10()
}

/// Compute the A-weighted magnitude spectrum from a set of frequency/magnitude pairs.
///
/// `frequencies` – frequencies in Hz.
/// `magnitudes` – linear magnitudes at each frequency.
///
/// Returns the A-weighted magnitudes.
#[must_use]
pub fn apply_a_weighting(frequencies: &[f64], magnitudes: &[f64]) -> Vec<f64> {
    frequencies
        .iter()
        .zip(magnitudes.iter())
        .map(|(&f, &m)| m * a_weighting(f))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_weighting_at_1khz_is_unity() {
        let w = a_weighting(1000.0);
        assert!((w - 1.0).abs() < 0.02);
    }

    #[test]
    fn a_weighting_low_freq_attenuated() {
        let w = a_weighting(20.0);
        assert!(w < 0.05); // heavily attenuated at 20 Hz
    }

    #[test]
    fn c_weighting_at_1khz_is_unity() {
        let w = c_weighting(1000.0);
        assert!((w - 1.0).abs() < 0.02);
    }

    #[test]
    fn b_weighting_at_1khz_near_unity() {
        let w = b_weighting(1000.0);
        assert!((w - 1.0).abs() < 0.05);
    }

    #[test]
    fn weighting_db_unity_is_zero() {
        assert!(weighting_db(1.0).abs() < 1e-10);
    }

    #[test]
    fn apply_a_weighting_length() {
        let freqs = vec![100.0, 1000.0, 10000.0];
        let mags = vec![1.0, 1.0, 1.0];
        let result = apply_a_weighting(&freqs, &mags);
        assert_eq!(result.len(), 3);
        // 1 kHz should be close to 1.0
        assert!((result[1] - 1.0).abs() < 0.02);
    }
}
