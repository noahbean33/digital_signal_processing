/// Compute the arithmetic mean of a signal.
#[must_use]
pub fn mean(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    signal.iter().sum::<f64>() / signal.len() as f64
}

/// Compute the variance of a signal.
///
/// Uses population variance (divide by N).
#[must_use]
pub fn variance(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    let m = mean(signal);
    signal.iter().map(|&x| (x - m) * (x - m)).sum::<f64>() / signal.len() as f64
}

/// Compute the sample variance of a signal (Bessel-corrected, divide by N-1).
#[must_use]
pub fn sample_variance(signal: &[f64]) -> f64 {
    if signal.len() < 2 {
        return 0.0;
    }
    let m = mean(signal);
    signal.iter().map(|&x| (x - m) * (x - m)).sum::<f64>() / (signal.len() - 1) as f64
}

/// Compute the standard deviation of a signal (population).
#[must_use]
pub fn std_dev(signal: &[f64]) -> f64 {
    variance(signal).sqrt()
}

/// Compute the sample standard deviation (Bessel-corrected).
#[must_use]
pub fn sample_std_dev(signal: &[f64]) -> f64 {
    sample_variance(signal).sqrt()
}

/// Find the minimum value in a signal.
#[must_use]
pub fn min(signal: &[f64]) -> f64 {
    signal
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min)
}

/// Find the maximum value in a signal.
#[must_use]
pub fn max(signal: &[f64]) -> f64 {
    signal
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max)
}

/// Compute the total energy of a signal: Σ x[n]².
#[must_use]
pub fn energy(signal: &[f64]) -> f64 {
    signal.iter().map(|&x| x * x).sum()
}

/// Compute the average power of a signal: (1/N) Σ x[n]².
#[must_use]
pub fn power(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    energy(signal) / signal.len() as f64
}

/// Compute the root mean square (RMS) of a signal.
#[must_use]
pub fn rms(signal: &[f64]) -> f64 {
    power(signal).sqrt()
}

/// Compute the peak-to-peak amplitude of a signal.
#[must_use]
pub fn peak_to_peak(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    max(signal) - min(signal)
}

/// Compute the crest factor (peak / RMS).
#[must_use]
pub fn crest_factor(signal: &[f64]) -> f64 {
    let r = rms(signal);
    if r < 1e-30 {
        return 0.0;
    }
    let peak = signal.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    peak / r
}

// ─── dB Conversion Utilities ──────────────────────────────────────────────────

/// Convert an amplitude (voltage) value to decibels: 20 * log10(|amp|).
///
/// Returns `f64::NEG_INFINITY` for zero amplitude.
#[must_use]
pub fn amp_to_db(amp: f64) -> f64 {
    let abs_amp = amp.abs();
    if abs_amp < 1e-30 {
        return f64::NEG_INFINITY;
    }
    20.0 * abs_amp.log10()
}

/// Convert a decibel value to amplitude (voltage): 10^(dB / 20).
#[must_use]
pub fn db_to_amp(db: f64) -> f64 {
    10.0_f64.powf(db / 20.0)
}

/// Convert a power value to decibels: 10 * log10(|power|).
///
/// Returns `f64::NEG_INFINITY` for zero power.
#[must_use]
pub fn power_to_db(power: f64) -> f64 {
    let abs_power = power.abs();
    if abs_power < 1e-30 {
        return f64::NEG_INFINITY;
    }
    10.0 * abs_power.log10()
}

/// Convert a decibel value to power: 10^(dB / 10).
#[must_use]
pub fn db_to_power(db: f64) -> f64 {
    10.0_f64.powf(db / 10.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_of_constant() {
        assert!((mean(&[5.0, 5.0, 5.0]) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn mean_of_empty() {
        assert!((mean(&[]) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn variance_of_constant_is_zero() {
        assert!(variance(&[3.0, 3.0, 3.0]).abs() < 1e-10);
    }

    #[test]
    fn std_dev_known_value() {
        let sig = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let sd = std_dev(&sig);
        assert!((sd - 2.0).abs() < 0.1);
    }

    #[test]
    fn sample_variance_bessel() {
        let sig = vec![2.0, 4.0, 6.0];
        let sv = sample_variance(&sig);
        // mean=4, deviations: 4,0,4 → sum=8, /(3-1)=4
        assert!((sv - 4.0).abs() < 1e-10);
    }

    #[test]
    fn min_max_correct() {
        let sig = vec![3.0, 1.0, 4.0, 1.0, 5.0];
        assert!((min(&sig) - 1.0).abs() < 1e-10);
        assert!((max(&sig) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn energy_correct() {
        let sig = vec![1.0, 2.0, 3.0];
        assert!((energy(&sig) - 14.0).abs() < 1e-10);
    }

    #[test]
    fn power_correct() {
        let sig = vec![1.0, 2.0, 3.0];
        assert!((power(&sig) - 14.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn rms_dc_signal() {
        assert!((rms(&[3.0, 3.0, 3.0]) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn peak_to_peak_range() {
        let sig = vec![-2.0, 0.0, 3.0];
        assert!((peak_to_peak(&sig) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn amp_db_roundtrip() {
        let amp = 2.5;
        let db = amp_to_db(amp);
        let recovered = db_to_amp(db);
        assert!((recovered - amp).abs() < 1e-10);
    }

    #[test]
    fn power_db_roundtrip() {
        let pwr = 100.0;
        let db = power_to_db(pwr);
        assert!((db - 20.0).abs() < 1e-10); // 10*log10(100) = 20
        let recovered = db_to_power(db);
        assert!((recovered - pwr).abs() < 1e-6);
    }

    #[test]
    fn amp_to_db_unity_is_zero() {
        assert!(amp_to_db(1.0).abs() < 1e-10);
    }

    #[test]
    fn amp_to_db_zero_is_neg_inf() {
        assert!(amp_to_db(0.0).is_infinite());
        assert!(amp_to_db(0.0) < 0.0);
    }

    #[test]
    fn db_to_amp_zero_is_one() {
        assert!((db_to_amp(0.0) - 1.0).abs() < 1e-10);
    }
}
