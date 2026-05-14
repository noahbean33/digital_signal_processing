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

// ─── Histogram Functions ──────────────────────────────────────────────────────

/// Compute a histogram of `signal` values into `num_bins` equally-spaced bins.
///
/// Bins span `[min_val, max_val]`. Values outside the range are clamped to the
/// nearest bin. Returns a vector of counts (as `f64` for easy normalisation).
#[must_use]
pub fn histogram(signal: &[f64], min_val: f64, max_val: f64, num_bins: usize) -> Vec<f64> {
    if signal.is_empty() || num_bins == 0 || max_val <= min_val {
        return Vec::new();
    }
    let mut bins = vec![0.0; num_bins];
    let range = max_val - min_val;
    for &x in signal {
        let frac = (x - min_val) / range;
        let idx = (frac * num_bins as f64).floor() as isize;
        let idx = idx.max(0).min(num_bins as isize - 1) as usize;
        bins[idx] += 1.0;
    }
    bins
}

/// Compute a cumulative histogram of `signal` values.
///
/// Each bin `i` contains the count of samples in bins `0..=i`.
#[must_use]
pub fn histogram_cumulative(signal: &[f64], min_val: f64, max_val: f64, num_bins: usize) -> Vec<f64> {
    let hist = histogram(signal, min_val, max_val, num_bins);
    if hist.is_empty() {
        return Vec::new();
    }
    let mut cum = Vec::with_capacity(hist.len());
    let mut sum = 0.0;
    for &h in &hist {
        sum += h;
        cum.push(sum);
    }
    cum
}

/// Histogram equalisation: remap signal values so the histogram is approximately
/// uniform over `[0, new_peak]`.
#[must_use]
pub fn histogram_equalize(signal: &[f64], new_peak: f64) -> Vec<f64> {
    if signal.is_empty() {
        return Vec::new();
    }
    let n = signal.len() as f64;
    let min_v = min(signal);
    let max_v = max(signal);
    if (max_v - min_v).abs() < 1e-30 {
        return vec![new_peak / 2.0; signal.len()];
    }
    let num_bins = 256;
    let cum = histogram_cumulative(signal, min_v, max_v, num_bins);
    let range = max_v - min_v;

    signal
        .iter()
        .map(|&x| {
            let frac = (x - min_v) / range;
            let idx = (frac * num_bins as f64).floor() as isize;
            let idx = idx.max(0).min(num_bins as isize - 1) as usize;
            cum[idx] / n * new_peak
        })
        .collect()
}

// ─── Peak-to-Average Ratio ────────────────────────────────────────────────────

/// Peak-to-average ratio (PAR).
///
/// `PAR = max(|x|) / mean(|x|)`
#[must_use]
pub fn peak_to_average_ratio(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    let peak = signal.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    let avg = signal.iter().map(|x| x.abs()).sum::<f64>() / signal.len() as f64;
    if avg < 1e-30 {
        return 0.0;
    }
    peak / avg
}

/// Peak-to-average power ratio (PAPR).
///
/// `PAPR = max(|x|²) / mean(|x|²)`
#[must_use]
pub fn peak_to_average_power_ratio(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    let peak_sq = signal.iter().map(|x| x * x).fold(0.0_f64, f64::max);
    let avg_sq = signal.iter().map(|x| x * x).sum::<f64>() / signal.len() as f64;
    if avg_sq < 1e-30 {
        return 0.0;
    }
    peak_sq / avg_sq
}

/// Peak-to-average power ratio in decibels.
///
/// `PAPR_dB = 10 * log10(PAPR)`
#[must_use]
pub fn peak_to_average_power_ratio_db(signal: &[f64]) -> f64 {
    let papr = peak_to_average_power_ratio(signal);
    if papr < 1e-30 {
        return f64::NEG_INFINITY;
    }
    10.0 * papr.log10()
}

/// dBm to voltage conversion.
///
/// `V = sqrt(Z * 10^((dBm - 30) / 10))`
///
/// * `zero_dbm_level` – impedance in ohms (typically 50 Ω)
#[must_use]
pub fn dbm_to_voltage(dbm: f64, impedance: f64) -> f64 {
    (impedance * 10.0_f64.powf((dbm - 30.0) / 10.0)).sqrt()
}

/// Voltage to dBm conversion.
///
/// `dBm = 10 * log10(V² / Z) + 30`
#[must_use]
pub fn voltage_to_dbm(voltage: f64, impedance: f64) -> f64 {
    if impedance < 1e-30 {
        return f64::NEG_INFINITY;
    }
    10.0 * (voltage * voltage / impedance).log10() + 30.0
}

// ─── Covariance Utilities ─────────────────────────────────────────────────────

/// Compute the sample covariance matrix of multi-channel data.
///
/// * `data` – flat row-major (channels × samples)
/// * `channels` – number of channels/features
/// * `samples` – number of observations
///
/// Returns flat row-major (channels × channels) covariance matrix.
#[must_use]
pub fn covariance_matrix(data: &[f64], channels: usize, samples: usize) -> Vec<f64> {
    assert_eq!(data.len(), channels * samples);
    let denom = if samples > 1 { (samples - 1) as f64 } else { 1.0 };

    // Compute means
    let mut means = vec![0.0; channels];
    for ch in 0..channels {
        let start = ch * samples;
        means[ch] = data[start..start + samples].iter().sum::<f64>() / samples as f64;
    }

    let mut cov = vec![0.0; channels * channels];
    for i in 0..channels {
        for j in i..channels {
            let mut sum = 0.0;
            for s in 0..samples {
                sum += (data[i * samples + s] - means[i]) * (data[j * samples + s] - means[j]);
            }
            let val = sum / denom;
            cov[i * channels + j] = val;
            cov[j * channels + i] = val;
        }
    }
    cov
}

/// Apply shrinkage regularization to a covariance matrix.
///
/// `C_reg = (1 - γ) * C + γ * mean(eigenvalues) * I`
///
/// This stabilizes ill-conditioned covariance matrices by shrinking them
/// toward a scaled identity matrix.
///
/// * `cov` – flat row-major (n × n) covariance matrix
/// * `n` – dimension
/// * `gamma` – shrinkage proportion (0.0 to 1.0, typical: 0.01)
#[must_use]
pub fn covariance_regularize(cov: &[f64], n: usize, gamma: f64) -> Vec<f64> {
    assert_eq!(cov.len(), n * n);
    let gamma = gamma.clamp(0.0, 1.0);

    // Estimate mean eigenvalue from the trace (sum of diagonal = sum of eigenvalues)
    let trace: f64 = (0..n).map(|i| cov[i * n + i]).sum();
    let mean_eval = trace / n as f64;

    let mut reg = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            reg[i * n + j] = (1.0 - gamma) * cov[i * n + j];
        }
        reg[i * n + i] += gamma * mean_eval;
    }
    reg
}

/// Compute the Frobenius distance between a single-trial covariance and
/// a reference covariance matrix.
///
/// `dist = sqrt(sum((A_flat · B_flat)))`
///
/// Used for artifact detection: trials with unusually large distance
/// from the average covariance are likely artifacts.
///
/// * `cov_trial` – flat (n × n) trial covariance
/// * `cov_ref` – flat (n × n) reference/average covariance
/// * `n` – dimension
#[must_use]
pub fn covariance_distance(cov_trial: &[f64], cov_ref: &[f64], n: usize) -> f64 {
    assert_eq!(cov_trial.len(), n * n);
    assert_eq!(cov_ref.len(), n * n);

    let dot: f64 = cov_trial.iter().zip(cov_ref.iter()).map(|(a, b)| a * b).sum();
    dot.abs().sqrt()
}

/// Compute z-scored covariance distances for a set of trials.
///
/// Each trial's covariance matrix is compared to the average covariance.
/// Returns z-scores that can be thresholded for artifact rejection.
///
/// * `trials` – vector of flat (n × n) covariance matrices, one per trial
/// * `n` – dimension
///
/// Returns z-scores for each trial.
#[must_use]
pub fn covariance_distances_z(trials: &[Vec<f64>], n: usize) -> Vec<f64> {
    if trials.is_empty() {
        return Vec::new();
    }

    // Compute average covariance
    let n_trials = trials.len();
    let mut avg_cov = vec![0.0; n * n];
    for trial in trials {
        assert_eq!(trial.len(), n * n);
        for (i, &v) in trial.iter().enumerate() {
            avg_cov[i] += v;
        }
    }
    for v in &mut avg_cov {
        *v /= n_trials as f64;
    }

    // Compute distances
    let distances: Vec<f64> = trials
        .iter()
        .map(|trial| covariance_distance(trial, &avg_cov, n))
        .collect();

    // Z-score
    let mean_dist = distances.iter().sum::<f64>() / n_trials as f64;
    let std_dist = (distances.iter().map(|&d| (d - mean_dist).powi(2)).sum::<f64>()
        / n_trials as f64)
        .sqrt();

    if std_dist < 1e-30 {
        return vec![0.0; n_trials];
    }

    distances.iter().map(|&d| (d - mean_dist) / std_dist).collect()
}

// ─── Permutation Testing ──────────────────────────────────────────────────────

/// Result of a permutation test.
#[derive(Clone, Debug)]
pub struct PermutationTestResult {
    /// The observed test statistic.
    pub observed: f64,
    /// The null distribution of test statistics.
    pub null_distribution: Vec<f64>,
    /// Z-score of the observed statistic relative to the null.
    pub z_score: f64,
    /// P-value (proportion of null >= observed).
    pub p_value: f64,
}

/// Permutation test for comparing two groups of values.
///
/// Tests whether the difference in means between `group_a` and `group_b`
/// is statistically significant by randomly permuting group labels.
///
/// * `group_a` – first group of observations
/// * `group_b` – second group of observations
/// * `n_permutations` – number of random permutations
/// * `seed` – random seed for reproducibility
#[must_use]
pub fn permutation_test_means(
    group_a: &[f64],
    group_b: &[f64],
    n_permutations: usize,
    seed: u64,
) -> PermutationTestResult {
    let na = group_a.len();
    let nb = group_b.len();
    let n_total = na + nb;

    // Observed statistic: difference in means
    let mean_a = if na > 0 { group_a.iter().sum::<f64>() / na as f64 } else { 0.0 };
    let mean_b = if nb > 0 { group_b.iter().sum::<f64>() / nb as f64 } else { 0.0 };
    let observed = (mean_a - mean_b).abs();

    // Pool all data
    let mut pooled: Vec<f64> = Vec::with_capacity(n_total);
    pooled.extend_from_slice(group_a);
    pooled.extend_from_slice(group_b);

    // Simple LCG random number generator
    let mut rng_state = seed;
    let mut next_rand = || -> u64 {
        rng_state = rng_state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1442695040888963407);
        rng_state
    };

    // Fisher-Yates shuffle helper
    let mut null_distribution = Vec::with_capacity(n_permutations);
    let mut shuffled = pooled.clone();

    for _ in 0..n_permutations {
        // Shuffle
        for i in (1..n_total).rev() {
            let j = (next_rand() as usize) % (i + 1);
            shuffled.swap(i, j);
        }

        // Compute statistic on shuffled data
        let perm_mean_a = shuffled[..na].iter().sum::<f64>() / na as f64;
        let perm_mean_b = shuffled[na..].iter().sum::<f64>() / nb as f64;
        null_distribution.push((perm_mean_a - perm_mean_b).abs());
    }

    // Compute z-score and p-value
    let null_mean = null_distribution.iter().sum::<f64>() / n_permutations as f64;
    let null_std = (null_distribution
        .iter()
        .map(|&x| (x - null_mean).powi(2))
        .sum::<f64>()
        / n_permutations as f64)
        .sqrt();

    let z_score = if null_std > 1e-30 {
        (observed - null_mean) / null_std
    } else {
        0.0
    };

    let p_value = null_distribution
        .iter()
        .filter(|&&x| x >= observed)
        .count() as f64
        / n_permutations as f64;

    PermutationTestResult {
        observed,
        null_distribution,
        z_score,
        p_value,
    }
}

/// Permutation test for a single statistic value against a null distribution.
///
/// Given an observed statistic and a null distribution (e.g., from shuffled data),
/// computes the z-score and p-value.
///
/// * `observed` – the observed test statistic
/// * `null_distribution` – vector of statistics from permuted data
#[must_use]
pub fn permutation_significance(observed: f64, null_distribution: &[f64]) -> PermutationTestResult {
    let n = null_distribution.len();
    if n == 0 {
        return PermutationTestResult {
            observed,
            null_distribution: Vec::new(),
            z_score: 0.0,
            p_value: 1.0,
        };
    }

    let null_mean = null_distribution.iter().sum::<f64>() / n as f64;
    let null_std = (null_distribution
        .iter()
        .map(|&x| (x - null_mean).powi(2))
        .sum::<f64>()
        / n as f64)
        .sqrt();

    let z_score = if null_std > 1e-30 {
        (observed - null_mean) / null_std
    } else {
        0.0
    };

    let p_value = null_distribution
        .iter()
        .filter(|&&x| x >= observed)
        .count() as f64
        / n as f64;

    PermutationTestResult {
        observed,
        null_distribution: null_distribution.to_vec(),
        z_score,
        p_value,
    }
}

// ─── Median & Rank ────────────────────────────────────────────────────────────

/// Compute the median of a signal.
///
/// For even-length signals, returns the average of the two middle values.
#[must_use]
pub fn median(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    let mut sorted = signal.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

/// Rank-transform a signal (average ranks for ties).
///
/// Assigns fractional ranks 1..N to each element, with ties receiving the
/// average of the ranks they would span.
#[must_use]
pub fn rank_transform(signal: &[f64]) -> Vec<f64> {
    let n = signal.len();
    if n == 0 {
        return Vec::new();
    }

    // Create (value, original_index) pairs and sort by value
    let mut indexed: Vec<(f64, usize)> = signal.iter().copied().enumerate().map(|(i, v)| (v, i)).collect();
    indexed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut ranks = vec![0.0; n];
    let mut i = 0;
    while i < n {
        // Find the end of this group of ties
        let mut j = i + 1;
        while j < n && (indexed[j].0 - indexed[i].0).abs() < 1e-15 {
            j += 1;
        }
        // Average rank for this group (1-indexed)
        let avg_rank = (i + 1 + j) as f64 / 2.0;
        for item in indexed.iter().take(j).skip(i) {
            ranks[item.1] = avg_rank;
        }
        i = j;
    }
    ranks
}

// ─── MAD & Modified Z-score ──────────────────────────────────────────────────

/// Median Absolute Deviation (MAD).
///
/// `MAD = median(|x - median(x)|)`
///
/// A robust measure of dispersion, resistant to outliers.
#[must_use]
pub fn mad(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    let med = median(signal);
    let deviations: Vec<f64> = signal.iter().map(|&x| (x - med).abs()).collect();
    median(&deviations)
}

/// Modified Z-score using MAD instead of standard deviation.
///
/// `Mz = 0.6745 * (x - median(x)) / MAD`
///
/// The constant 0.6745 is `Φ⁻¹(3/4)`, making the MAD consistent with
/// standard deviation for normally distributed data.
///
/// Useful for outlier detection in non-normal distributions.
#[must_use]
pub fn modified_z_score(signal: &[f64]) -> Vec<f64> {
    if signal.is_empty() {
        return Vec::new();
    }
    let med = median(signal);
    let mad_val = mad(signal);
    if mad_val < 1e-30 {
        return vec![0.0; signal.len()];
    }
    // 0.6745 ≈ inverse normal CDF at 0.75
    let k = 0.674_489_750_196_082;
    signal.iter().map(|&x| k * (x - med) / mad_val).collect()
}

// ─── Z-score & Normalizations ─────────────────────────────────────────────────

/// Z-score normalization (standardization).
///
/// `z = (x - mean) / std`
///
/// Transforms data to zero mean and unit variance.
///
/// * `ddof` – delta degrees of freedom for std (0 = population, 1 = sample)
#[must_use]
pub fn z_score(signal: &[f64], ddof: usize) -> Vec<f64> {
    if signal.len() < 2 {
        return vec![0.0; signal.len()];
    }
    let m = mean(signal);
    let n = signal.len();
    let denom = if ddof == 0 { n } else { n - ddof.min(n - 1) };
    let var = signal.iter().map(|&x| (x - m) * (x - m)).sum::<f64>() / denom as f64;
    let s = var.sqrt();
    if s < 1e-30 {
        return vec![0.0; signal.len()];
    }
    signal.iter().map(|&x| (x - m) / s).collect()
}

/// Min-max scaling to an arbitrary range `[new_min, new_max]`.
///
/// `x_scaled = (x - min) / (max - min) * (new_max - new_min) + new_min`
#[must_use]
pub fn min_max_scale(signal: &[f64], new_min: f64, new_max: f64) -> Vec<f64> {
    if signal.is_empty() {
        return Vec::new();
    }
    let min_v = min(signal);
    let max_v = max(signal);
    let range = max_v - min_v;
    if range < 1e-30 {
        return vec![(new_min + new_max) / 2.0; signal.len()];
    }
    let new_range = new_max - new_min;
    signal.iter().map(|&x| (x - min_v) / range * new_range + new_min).collect()
}

/// Fisher-Z transform: `F(r) = arctanh(r)`.
///
/// Transforms correlation coefficients to approximately normally distributed
/// values, enabling standard statistical tests on correlations.
///
/// Input values should be in (-1, 1).
#[must_use]
pub fn fisher_z(r: f64) -> f64 {
    r.atanh()
}

/// Inverse Fisher-Z transform: `r = tanh(z)`.
#[must_use]
pub fn fisher_z_inv(z: f64) -> f64 {
    z.tanh()
}

/// Apply Fisher-Z transform to a vector of correlation coefficients.
#[must_use]
pub fn fisher_z_vec(correlations: &[f64]) -> Vec<f64> {
    correlations.iter().map(|&r| r.atanh()).collect()
}

// ─── Entropy ──────────────────────────────────────────────────────────────────

/// Shannon entropy of a discrete probability distribution.
///
/// `H = -Σ p * log2(p)`
///
/// * `probabilities` – probability values (should sum to ~1)
///
/// Returns entropy in bits. Zero-probability entries are skipped.
#[must_use]
pub fn shannon_entropy(probabilities: &[f64]) -> f64 {
    let mut h = 0.0;
    for &p in probabilities {
        if p > 1e-30 {
            h -= p * p.log2();
        }
    }
    h
}

/// Compute Shannon entropy of a signal by binning into a histogram.
///
/// * `signal` – input data
/// * `num_bins` – number of histogram bins
///
/// Returns entropy in bits.
#[must_use]
pub fn signal_entropy(signal: &[f64], num_bins: usize) -> f64 {
    if signal.is_empty() || num_bins == 0 {
        return 0.0;
    }
    let min_v = min(signal);
    let max_v = max(signal);
    if (max_v - min_v) < 1e-30 {
        return 0.0;
    }
    let hist = histogram(signal, min_v, max_v, num_bins);
    let n = signal.len() as f64;
    let probs: Vec<f64> = hist.iter().map(|&c| c / n).collect();
    shannon_entropy(&probs)
}

// ─── Correlation Functions ────────────────────────────────────────────────────

/// Pearson correlation coefficient between two signals.
///
/// `r = Σ((x-μx)(y-μy)) / sqrt(Σ(x-μx)² * Σ(y-μy)²)`
///
/// Returns a value in [-1, 1]. Returns 0 if either signal has zero variance.
#[must_use]
pub fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    if n < 2 {
        return 0.0;
    }
    let mx = x[..n].iter().sum::<f64>() / n as f64;
    let my = y[..n].iter().sum::<f64>() / n as f64;

    let mut num = 0.0;
    let mut dx2 = 0.0;
    let mut dy2 = 0.0;
    for i in 0..n {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        num += dx * dy;
        dx2 += dx * dx;
        dy2 += dy * dy;
    }

    let denom = (dx2 * dy2).sqrt();
    if denom < 1e-30 {
        return 0.0;
    }
    num / denom
}

/// Spearman rank correlation coefficient.
///
/// Computes the Pearson correlation of the rank-transformed data.
/// Measures monotonic (not necessarily linear) relationships.
#[must_use]
pub fn spearman_correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    if n < 2 {
        return 0.0;
    }
    let rx = rank_transform(&x[..n]);
    let ry = rank_transform(&y[..n]);
    pearson_correlation(&rx, &ry)
}

/// Kendall tau-b correlation coefficient.
///
/// Measures ordinal association: the difference between the probability
/// of concordant and discordant pairs.
#[must_use]
pub fn kendall_correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    if n < 2 {
        return 0.0;
    }

    let mut concordant: i64 = 0;
    let mut discordant: i64 = 0;
    let mut ties_x: i64 = 0;
    let mut ties_y: i64 = 0;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = x[i] - x[j];
            let dy = y[i] - y[j];
            let product = dx * dy;

            if product.abs() < 1e-15 {
                // At least one tie
                if dx.abs() < 1e-15 && dy.abs() < 1e-15 {
                    // Tie in both
                    ties_x += 1;
                    ties_y += 1;
                } else if dx.abs() < 1e-15 {
                    ties_x += 1;
                } else {
                    ties_y += 1;
                }
            } else if product > 0.0 {
                concordant += 1;
            } else {
                discordant += 1;
            }
        }
    }

    let n0 = (n * (n - 1)) as f64 / 2.0;
    let denom = ((n0 - ties_x as f64) * (n0 - ties_y as f64)).sqrt();
    if denom < 1e-30 {
        return 0.0;
    }
    (concordant - discordant) as f64 / denom
}

/// Cosine similarity between two vectors.
///
/// `cos_sim = (x · y) / (||x|| * ||y||)`
///
/// Unlike correlation, cosine similarity is sensitive to the mean/offset.
/// Returns a value in [-1, 1].
#[must_use]
pub fn cosine_similarity(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    if n == 0 {
        return 0.0;
    }
    let mut dot = 0.0;
    let mut nx = 0.0;
    let mut ny = 0.0;
    for i in 0..n {
        dot += x[i] * y[i];
        nx += x[i] * x[i];
        ny += y[i] * y[i];
    }
    let denom = (nx * ny).sqrt();
    if denom < 1e-30 {
        return 0.0;
    }
    dot / denom
}

/// Compute a correlation matrix for multi-channel data.
///
/// * `data` – flat row-major (channels × samples)
/// * `channels` – number of channels/features
/// * `samples` – number of observations
///
/// Returns flat row-major (channels × channels) Pearson correlation matrix.
#[must_use]
pub fn correlation_matrix(data: &[f64], channels: usize, samples: usize) -> Vec<f64> {
    assert_eq!(data.len(), channels * samples);
    let mut corr = vec![0.0; channels * channels];
    for i in 0..channels {
        corr[i * channels + i] = 1.0;
        for j in (i + 1)..channels {
            let r = pearson_correlation(
                &data[i * samples..(i + 1) * samples],
                &data[j * samples..(j + 1) * samples],
            );
            corr[i * channels + j] = r;
            corr[j * channels + i] = r;
        }
    }
    corr
}

// ─── Signal Detection Theory ─────────────────────────────────────────────────

/// Compute d-prime (d') for signal detection theory.
///
/// `d' = Z(hit_rate) - Z(false_alarm_rate)`
///
/// where Z is the inverse normal CDF (probit function).
///
/// * `hit_rate` – proportion of hits (0 < hr < 1)
/// * `false_alarm_rate` – proportion of false alarms (0 < far < 1)
///
/// Rates of exactly 0 or 1 are clamped to 0.001 and 0.999.
#[must_use]
pub fn d_prime(hit_rate: f64, false_alarm_rate: f64) -> f64 {
    let hr = hit_rate.clamp(0.001, 0.999);
    let far = false_alarm_rate.clamp(0.001, 0.999);
    probit(hr) - probit(far)
}

/// Probit function (inverse normal CDF) using rational approximation.
///
/// Abramowitz & Stegun approximation for the inverse of the standard
/// normal cumulative distribution function.
fn probit(p: f64) -> f64 {
    // Ensure p is in valid range
    let p = p.clamp(1e-10, 1.0 - 1e-10);

    // Use symmetry: if p > 0.5, compute for 1-p and negate
    if p < 0.5 {
        -rational_approx((-2.0 * p.ln()).sqrt())
    } else {
        rational_approx((-2.0 * (1.0 - p).ln()).sqrt())
    }
}

/// Rational approximation helper for probit.
fn rational_approx(t: f64) -> f64 {
    // Coefficients from Peter Acklam's approximation
    let c0 = 2.515_517;
    let c1 = 0.802_853;
    let c2 = 0.010_328;
    let d1 = 1.432_788;
    let d2 = 0.189_269;
    let d3 = 0.001_308;

    t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t)
}

// ─── Bootstrap Confidence Intervals ──────────────────────────────────────────

/// Result of a bootstrap confidence interval computation.
#[derive(Clone, Debug)]
pub struct BootstrapCIResult {
    /// Sample statistic (e.g., mean of the original data).
    pub statistic: f64,
    /// Lower bound of the confidence interval.
    pub ci_lower: f64,
    /// Upper bound of the confidence interval.
    pub ci_upper: f64,
    /// Bootstrap distribution of the statistic.
    pub bootstrap_distribution: Vec<f64>,
}

/// Bootstrap confidence interval for the mean.
///
/// Resamples the data with replacement `n_bootstrap` times, computes the mean
/// of each resample, and finds the percentile-based confidence interval.
///
/// * `data` – input data
/// * `confidence` – confidence level (e.g., 95.0 for 95%)
/// * `n_bootstrap` – number of bootstrap resamples
/// * `seed` – random seed
#[must_use]
pub fn bootstrap_ci_mean(
    data: &[f64],
    confidence: f64,
    n_bootstrap: usize,
    seed: u64,
) -> BootstrapCIResult {
    let n = data.len();
    if n == 0 {
        return BootstrapCIResult {
            statistic: 0.0,
            ci_lower: 0.0,
            ci_upper: 0.0,
            bootstrap_distribution: Vec::new(),
        };
    }

    let statistic = mean(data);

    // LCG RNG
    let mut rng_state = seed;
    let mut next_rand = || -> u64 {
        rng_state = rng_state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        rng_state
    };

    let mut boot_means = Vec::with_capacity(n_bootstrap);
    for _ in 0..n_bootstrap {
        let mut sum = 0.0;
        for _ in 0..n {
            let idx = (next_rand() as usize) % n;
            sum += data[idx];
        }
        boot_means.push(sum / n as f64);
    }

    // Sort for percentile computation
    boot_means.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let alpha = (100.0 - confidence) / 200.0;
    let lower_idx = (alpha * n_bootstrap as f64).floor() as usize;
    let upper_idx = ((1.0 - alpha) * n_bootstrap as f64).ceil() as usize;

    let ci_lower = boot_means[lower_idx.min(n_bootstrap - 1)];
    let ci_upper = boot_means[upper_idx.min(n_bootstrap - 1)];

    BootstrapCIResult {
        statistic,
        ci_lower,
        ci_upper,
        bootstrap_distribution: boot_means,
    }
}

// ─── T-tests ──────────────────────────────────────────────────────────────────

/// Result of a t-test.
#[derive(Clone, Debug)]
pub struct TTestResult {
    /// The t-statistic.
    pub t_statistic: f64,
    /// Degrees of freedom.
    pub df: f64,
    /// Two-tailed p-value (approximate).
    pub p_value: f64,
}

/// One-sample t-test: test whether the mean of a sample differs from a value.
///
/// * `data` – sample data
/// * `population_mean` – hypothesized population mean (often 0)
#[must_use]
pub fn t_test_one_sample(data: &[f64], population_mean: f64) -> TTestResult {
    let n = data.len();
    if n < 2 {
        return TTestResult { t_statistic: 0.0, df: 0.0, p_value: 1.0 };
    }
    let m = mean(data);
    let s = sample_std_dev(data);
    let se = s / (n as f64).sqrt();
    let t = if se > 1e-30 { (m - population_mean) / se } else { 0.0 };
    let df = (n - 1) as f64;
    let p = t_to_p(t.abs(), df);

    TTestResult { t_statistic: t, df, p_value: p }
}

/// Independent two-sample t-test (Welch's t-test, unequal variances).
///
/// Tests whether the means of two independent samples differ.
#[must_use]
pub fn t_test_two_sample(group_a: &[f64], group_b: &[f64]) -> TTestResult {
    let na = group_a.len();
    let nb = group_b.len();
    if na < 2 || nb < 2 {
        return TTestResult { t_statistic: 0.0, df: 0.0, p_value: 1.0 };
    }

    let ma = mean(group_a);
    let mb = mean(group_b);
    let va = sample_variance(group_a);
    let vb = sample_variance(group_b);

    let se = (va / na as f64 + vb / nb as f64).sqrt();
    let t = if se > 1e-30 { (ma - mb) / se } else { 0.0 };

    // Welch-Satterthwaite degrees of freedom
    let num = (va / na as f64 + vb / nb as f64).powi(2);
    let den = (va / na as f64).powi(2) / (na - 1) as f64
        + (vb / nb as f64).powi(2) / (nb - 1) as f64;
    let df = if den > 1e-30 { num / den } else { 1.0 };
    let p = t_to_p(t.abs(), df);

    TTestResult { t_statistic: t, df, p_value: p }
}

/// Approximate two-tailed p-value from t-statistic using a normal approximation
/// for large df, or a simple beta incomplete function approximation.
fn t_to_p(t_abs: f64, df: f64) -> f64 {
    // For large df, t → z, use normal CDF approximation
    // Use the approximation: p ≈ 2 * (1 - Φ(|t| * sqrt(df/(df-2+t²)))) for rough estimate
    // More accurate: use regularized incomplete beta function
    let x = df / (df + t_abs * t_abs);
    // Approximate regularized incomplete beta I_x(df/2, 1/2) using series
    let p = regularized_incomplete_beta(x, df / 2.0, 0.5);
    p.clamp(0.0, 1.0)
}

/// Simple approximation of the regularized incomplete beta function I_x(a, b).
/// Used internally for p-value computation from t-statistics.
fn regularized_incomplete_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }

    // Use continued fraction (Lentz's algorithm) for better accuracy
    let mut f = 1.0;
    let mut c = 1.0;
    let mut d;

    let max_iter = 200;
    let eps = 1e-10;

    // First, compute ln(B(a,b)) using Stirling for the prefactor
    let log_prefix = a * x.ln() + b * (1.0 - x).ln()
        - (a.ln() + ln_beta(a, b));

    let prefix = log_prefix.exp();

    // Continued fraction for I_x(a, b) = prefix * 1/(1+ d1/(1+ d2/(1+ ...)))
    d = 1.0;
    for m in 0..max_iter {
        let m_f = m as f64;
        // Compute numerator coefficient
        let numerator = if m == 0 {
            1.0
        } else if m % 2 == 0 {
            let k = m_f / 2.0;
            k * (b - k) * x / ((a + 2.0 * k - 1.0) * (a + 2.0 * k))
        } else {
            let k = (m_f - 1.0) / 2.0;
            -(a + k) * (a + b + k) * x / ((a + 2.0 * k) * (a + 2.0 * k + 1.0))
        };

        d = 1.0 + numerator * d;
        if d.abs() < 1e-30 {
            d = 1e-30;
        }
        d = 1.0 / d;

        c = 1.0 + numerator / c;
        if c.abs() < 1e-30 {
            c = 1e-30;
        }

        let delta = c * d;
        f *= delta;

        if (delta - 1.0).abs() < eps {
            break;
        }
    }

    (prefix * f / a).clamp(0.0, 1.0)
}

/// Natural log of the Beta function using log-gamma (Stirling approximation).
fn ln_beta(a: f64, b: f64) -> f64 {
    ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
}

/// Stirling approximation for ln(Γ(x)).
fn ln_gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    // Lanczos approximation coefficients
    let g = 7.0;
    let coef = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];

    let x = x - 1.0;
    let mut sum = coef[0];
    for (i, &c) in coef.iter().enumerate().skip(1) {
        sum += c / (x + i as f64);
    }

    let t = x + g + 0.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + sum.ln()
}

// ─── Dispersion Measures ──────────────────────────────────────────────────────

/// Coefficient of Variation (CV).
///
/// `CV = std / mean`
///
/// A normalized measure of dispersion. Meaningful only for positive-valued data.
#[must_use]
pub fn coefficient_of_variation(signal: &[f64]) -> f64 {
    let m = mean(signal);
    if m.abs() < 1e-30 {
        return 0.0;
    }
    std_dev(signal) / m.abs()
}

/// Fano factor.
///
/// `F = variance / mean`
///
/// A measure of dispersion relative to a Poisson process (F=1 for Poisson).
/// Meaningful only for positive-valued (count) data.
#[must_use]
pub fn fano_factor(signal: &[f64]) -> f64 {
    let m = mean(signal);
    if m.abs() < 1e-30 {
        return 0.0;
    }
    variance(signal) / m.abs()
}

// ─── Mutual Information ──────────────────────────────────────────────────────

/// Result of a mutual information computation.
#[derive(Clone, Debug)]
pub struct MutualInfoResult {
    /// Mutual information in bits.
    pub mi: f64,
    /// Entropy of x in bits.
    pub entropy_x: f64,
    /// Entropy of y in bits.
    pub entropy_y: f64,
    /// Joint entropy of (x, y) in bits.
    pub joint_entropy: f64,
    /// Number of bins used for discretisation.
    pub num_bins: usize,
}

/// Compute the mutual information between two equal-length vectors.
///
/// Uses histogram binning to estimate the probability distributions.
/// If `num_bins` is 0, the number of bins is chosen automatically via the
/// Freedman-Diaconis rule.
///
/// MI(X;Y) = H(X) + H(Y) − H(X,Y)
#[must_use]
pub fn mutual_information(x: &[f64], y: &[f64], num_bins: usize) -> MutualInfoResult {
    let n = x.len();
    assert_eq!(n, y.len(), "x and y must have equal length");
    if n == 0 {
        return MutualInfoResult {
            mi: 0.0,
            entropy_x: 0.0,
            entropy_y: 0.0,
            joint_entropy: 0.0,
            num_bins: 0,
        };
    }

    // Determine bins via Freedman-Diaconis if num_bins == 0
    let bins = if num_bins == 0 {
        let fd = |data: &[f64]| -> usize {
            let mut sorted = data.to_vec();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let q1 = sorted[sorted.len() / 4];
            let q3 = sorted[3 * sorted.len() / 4];
            let iqr = q3 - q1;
            let range = sorted.last().unwrap() - sorted.first().unwrap();
            if iqr.abs() < 1e-30 || range.abs() < 1e-30 {
                return 10;
            }
            let bin_width = 2.0 * iqr * (n as f64).powf(-1.0 / 3.0);
            (range / bin_width).ceil() as usize
        };
        let b1 = fd(x);
        let b2 = fd(y);
        ((b1 + b2) / 2).max(2)
    } else {
        num_bins
    };

    // Bin assignment helper
    let bin_assign = |data: &[f64]| -> Vec<usize> {
        let min_val = data.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_val = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let range = max_val - min_val;
        if range.abs() < 1e-30 {
            return vec![0; data.len()];
        }
        data.iter()
            .map(|&v| {
                let b = ((v - min_val) / range * (bins as f64 - 1.0)).round() as usize;
                b.min(bins - 1)
            })
            .collect()
    };

    let bins_x = bin_assign(x);
    let bins_y = bin_assign(y);

    // Marginal distributions
    let mut px = vec![0.0; bins];
    let mut py = vec![0.0; bins];
    for i in 0..n {
        px[bins_x[i]] += 1.0;
        py[bins_y[i]] += 1.0;
    }
    let n_f = n as f64;
    for v in px.iter_mut() {
        *v /= n_f;
    }
    for v in py.iter_mut() {
        *v /= n_f;
    }

    // Joint distribution
    let mut pxy = vec![0.0; bins * bins];
    for i in 0..n {
        pxy[bins_x[i] * bins + bins_y[i]] += 1.0;
    }
    for v in pxy.iter_mut() {
        *v /= n_f;
    }

    let eps = f64::EPSILON;
    let entropy_x = -px
        .iter()
        .map(|&p| if p > eps { p * p.log2() } else { 0.0 })
        .sum::<f64>();
    let entropy_y = -py
        .iter()
        .map(|&p| if p > eps { p * p.log2() } else { 0.0 })
        .sum::<f64>();
    let joint_entropy = -pxy
        .iter()
        .map(|&p| if p > eps { p * p.log2() } else { 0.0 })
        .sum::<f64>();

    let mi = entropy_x + entropy_y - joint_entropy;

    MutualInfoResult {
        mi,
        entropy_x,
        entropy_y,
        joint_entropy,
        num_bins: bins,
    }
}

/// Mutual information with permutation testing.
///
/// Computes the MI and then runs `n_permutations` shuffles of `y` to build
/// a null distribution.  Returns the MI as a z-score relative to the null.
#[must_use]
pub fn mutual_information_z(
    x: &[f64],
    y: &[f64],
    num_bins: usize,
    n_permutations: usize,
    seed: u64,
) -> f64 {
    let observed = mutual_information(x, y, num_bins).mi;
    if n_permutations == 0 {
        return observed;
    }

    let n = y.len();
    let mut y_shuffled = y.to_vec();
    let mut perm_mi = Vec::with_capacity(n_permutations);

    // Simple LCG PRNG for reproducible shuffling
    let mut rng_state = seed;
    let lcg_next = |state: &mut u64| -> u64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *state
    };

    for _ in 0..n_permutations {
        // Fisher-Yates shuffle
        for i in (1..n).rev() {
            let j = (lcg_next(&mut rng_state) as usize) % (i + 1);
            y_shuffled.swap(i, j);
        }
        perm_mi.push(mutual_information(x, &y_shuffled, num_bins).mi);
    }

    let perm_mean = perm_mi.iter().sum::<f64>() / perm_mi.len() as f64;
    let perm_var = perm_mi
        .iter()
        .map(|&v| (v - perm_mean) * (v - perm_mean))
        .sum::<f64>()
        / perm_mi.len() as f64;
    let perm_std = perm_var.sqrt();

    if perm_std.abs() < 1e-30 {
        return 0.0;
    }
    (observed - perm_mean) / perm_std
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

    #[test]
    fn histogram_uniform_distribution() {
        // 100 values evenly spread across [0, 10)
        let signal: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();
        let hist = histogram(&signal, 0.0, 10.0, 10);
        assert_eq!(hist.len(), 10);
        // Each bin should have 10 samples
        for &h in &hist {
            assert!((h - 10.0).abs() < 1e-10);
        }
    }

    #[test]
    fn histogram_cumulative_last_equals_total() {
        let signal: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let cum = histogram_cumulative(&signal, 0.0, 50.0, 5);
        assert_eq!(cum.len(), 5);
        assert!((*cum.last().unwrap() - 50.0).abs() < 1e-10);
        // Should be non-decreasing
        for i in 1..cum.len() {
            assert!(cum[i] >= cum[i - 1]);
        }
    }

    #[test]
    fn histogram_equalize_range() {
        let signal: Vec<f64> = (0..100).map(|i| (i as f64).powi(2)).collect();
        let eq = histogram_equalize(&signal, 1.0);
        assert_eq!(eq.len(), 100);
        // All values should be in [0, 1]
        for &v in &eq {
            assert!(v >= -1e-10 && v <= 1.0 + 1e-10);
        }
    }

    #[test]
    fn histogram_empty_inputs() {
        assert!(histogram(&[], 0.0, 1.0, 10).is_empty());
        assert!(histogram_cumulative(&[], 0.0, 1.0, 10).is_empty());
        assert!(histogram_equalize(&[], 1.0).is_empty());
    }

    #[test]
    fn par_of_dc_is_one() {
        let signal = vec![3.0; 100];
        assert!((peak_to_average_ratio(&signal) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn papr_of_dc_is_one() {
        let signal = vec![2.0; 100];
        assert!((peak_to_average_power_ratio(&signal) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn papr_db_of_dc_is_zero() {
        let signal = vec![1.0; 100];
        assert!((peak_to_average_power_ratio_db(&signal) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn papr_of_sine_about_3db() {
        use std::f64::consts::PI;
        let signal: Vec<f64> = (0..1000)
            .map(|i| (2.0 * PI * i as f64 / 1000.0).sin())
            .collect();
        let papr_db = peak_to_average_power_ratio_db(&signal);
        // PAPR of a sine wave is 10*log10(2) ≈ 3.01 dB
        assert!((papr_db - 3.01).abs() < 0.1);
    }

    #[test]
    fn dbm_voltage_roundtrip() {
        let v = 1.0; // 1 volt
        let z = 50.0;
        let dbm = voltage_to_dbm(v, z);
        let v_back = dbm_to_voltage(dbm, z);
        assert!((v_back - v).abs() < 1e-6);
    }

    #[test]
    fn covariance_matrix_symmetric() {
        let data = vec![
            1.0, 2.0, 3.0, 4.0, 5.0,
            5.0, 4.0, 3.0, 2.0, 1.0,
        ];
        let cov = covariance_matrix(&data, 2, 5);
        assert_eq!(cov.len(), 4);
        // Should be symmetric
        assert!((cov[1] - cov[2]).abs() < 1e-10);
        // Diagonal should be positive (variance)
        assert!(cov[0] > 0.0);
        assert!(cov[3] > 0.0);
    }

    #[test]
    fn covariance_matrix_identical_channels() {
        let data = vec![
            1.0, 2.0, 3.0, 4.0,
            1.0, 2.0, 3.0, 4.0,
        ];
        let cov = covariance_matrix(&data, 2, 4);
        // Perfectly correlated: off-diagonal should equal diagonal
        assert!((cov[0] - cov[1]).abs() < 1e-10);
    }

    #[test]
    fn covariance_regularize_increases_diagonal() {
        let cov = vec![2.0, 0.5, 0.5, 1.0];
        let reg = covariance_regularize(&cov, 2, 0.1);
        // Trace is preserved approximately
        let trace_orig = cov[0] + cov[3];
        let trace_reg = reg[0] + reg[3];
        assert!((trace_orig - trace_reg).abs() < 1e-6);
        // Off-diagonal shrinks
        assert!(reg[1].abs() < cov[1].abs());
    }

    #[test]
    fn covariance_distance_self_is_norm() {
        let cov = vec![2.0, 0.5, 0.5, 1.0];
        let dist = covariance_distance(&cov, &cov, 2);
        // Should be sqrt(sum of squares)
        let expected = cov.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((dist - expected).abs() < 1e-10);
    }

    #[test]
    fn covariance_distances_z_mean_is_zero() {
        let trials = vec![
            vec![1.0, 0.0, 0.0, 1.0],
            vec![2.0, 0.5, 0.5, 2.0],
            vec![1.5, 0.2, 0.2, 1.5],
            vec![1.0, 0.1, 0.1, 1.0],
        ];
        let z = covariance_distances_z(&trials, 2);
        assert_eq!(z.len(), 4);
        let mean_z: f64 = z.iter().sum::<f64>() / z.len() as f64;
        assert!(mean_z.abs() < 1e-10);
    }

    #[test]
    fn permutation_test_significant_difference() {
        // Large difference should yield small p-value
        let group_a = vec![10.0, 11.0, 12.0, 13.0, 14.0];
        let group_b = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let result = permutation_test_means(&group_a, &group_b, 500, 42);
        assert!(result.p_value < 0.05);
        assert!(result.z_score > 2.0);
    }

    #[test]
    fn permutation_test_no_difference() {
        // Same distribution should yield large p-value
        let group_a = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let group_b = vec![1.5, 2.5, 3.5, 4.5, 5.5];
        let result = permutation_test_means(&group_a, &group_b, 500, 42);
        // p-value should not be extremely small
        assert!(result.p_value > 0.01);
    }

    #[test]
    fn permutation_significance_works() {
        let null = vec![1.0, 1.5, 2.0, 2.5, 3.0, 1.2, 1.8, 2.2, 2.8, 1.6];
        let result = permutation_significance(10.0, &null);
        // 10.0 is way above null → p should be 0
        assert!(result.p_value < 0.01);
        assert!(result.z_score > 3.0);
    }

    // ─── Median & Rank Tests ─────────────────────────────────────────────────

    #[test]
    fn median_odd() {
        assert!((median(&[3.0, 1.0, 2.0]) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn median_even() {
        assert!((median(&[4.0, 1.0, 3.0, 2.0]) - 2.5).abs() < 1e-10);
    }

    #[test]
    fn median_empty() {
        assert!((median(&[]) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn rank_transform_no_ties() {
        let r = rank_transform(&[10.0, 30.0, 20.0]);
        assert!((r[0] - 1.0).abs() < 1e-10);
        assert!((r[1] - 3.0).abs() < 1e-10);
        assert!((r[2] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn rank_transform_with_ties() {
        let r = rank_transform(&[1.0, 2.0, 2.0, 4.0]);
        assert!((r[0] - 1.0).abs() < 1e-10);
        assert!((r[1] - 2.5).abs() < 1e-10); // average rank
        assert!((r[2] - 2.5).abs() < 1e-10);
        assert!((r[3] - 4.0).abs() < 1e-10);
    }

    // ─── MAD & Modified Z-score Tests ────────────────────────────────────────

    #[test]
    fn mad_symmetric() {
        // MAD of {1,2,3,4,5}: median=3, deviations={2,1,0,1,2}, MAD=1
        assert!((mad(&[1.0, 2.0, 3.0, 4.0, 5.0]) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn modified_z_score_center_is_zero() {
        let mz = modified_z_score(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        // Median element (3.0) should have mz = 0
        assert!(mz[2].abs() < 1e-10);
    }

    // ─── Z-score & Normalization Tests ───────────────────────────────────────

    #[test]
    fn z_score_zero_mean_unit_var() {
        let data = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        let z = z_score(&data, 0);
        let m: f64 = z.iter().sum::<f64>() / z.len() as f64;
        let v: f64 = z.iter().map(|&x| (x - m) * (x - m)).sum::<f64>() / z.len() as f64;
        assert!(m.abs() < 1e-10);
        assert!((v - 1.0).abs() < 1e-10);
    }

    #[test]
    fn min_max_scale_bounds() {
        let data = vec![1.0, 5.0, 3.0];
        let scaled = min_max_scale(&data, 0.0, 1.0);
        assert!((scaled[0] - 0.0).abs() < 1e-10); // min maps to 0
        assert!((scaled[1] - 1.0).abs() < 1e-10); // max maps to 1
        assert!((scaled[2] - 0.5).abs() < 1e-10); // middle maps to 0.5
    }

    #[test]
    fn min_max_scale_custom_range() {
        let data = vec![0.0, 10.0];
        let scaled = min_max_scale(&data, -1.0, 1.0);
        assert!((scaled[0] - (-1.0)).abs() < 1e-10);
        assert!((scaled[1] - 1.0).abs() < 1e-10);
    }

    // ─── Fisher-Z Tests ─────────────────────────────────────────────────────

    #[test]
    fn fisher_z_roundtrip() {
        let r = 0.7;
        let z = fisher_z(r);
        let r_back = fisher_z_inv(z);
        assert!((r - r_back).abs() < 1e-10);
    }

    #[test]
    fn fisher_z_zero() {
        assert!(fisher_z(0.0).abs() < 1e-10);
    }

    #[test]
    fn fisher_z_vec_values() {
        let rs = vec![0.0, 0.5, -0.5];
        let zs = fisher_z_vec(&rs);
        assert!(zs[0].abs() < 1e-10);
        assert!((zs[1] - (-zs[2])).abs() < 1e-10); // antisymmetric
    }

    // ─── Entropy Tests ───────────────────────────────────────────────────────

    #[test]
    fn shannon_entropy_uniform() {
        // Uniform over 8 outcomes: H = log2(8) = 3
        let p = vec![0.125; 8];
        assert!((shannon_entropy(&p) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn shannon_entropy_certain() {
        // One certain outcome: H = 0
        let p = vec![1.0, 0.0, 0.0];
        assert!(shannon_entropy(&p).abs() < 1e-10);
    }

    #[test]
    fn signal_entropy_positive() {
        let data: Vec<f64> = (0..100).map(|i| (i as f64 * 0.1).sin()).collect();
        let h = signal_entropy(&data, 10);
        assert!(h > 0.0);
    }

    // ─── Correlation Tests ───────────────────────────────────────────────────

    #[test]
    fn pearson_perfect_positive() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        assert!((pearson_correlation(&x, &y) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn pearson_perfect_negative() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![10.0, 8.0, 6.0, 4.0, 2.0];
        assert!((pearson_correlation(&x, &y) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn spearman_monotonic() {
        // Spearman should be 1.0 for any monotonically increasing relationship
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![1.0, 8.0, 27.0, 64.0, 125.0]; // y = x³
        assert!((spearman_correlation(&x, &y) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn kendall_concordant() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((kendall_correlation(&x, &y) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn cosine_similarity_identical() {
        let x = vec![1.0, 2.0, 3.0];
        assert!((cosine_similarity(&x, &x) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn cosine_similarity_orthogonal() {
        let x = vec![1.0, 0.0];
        let y = vec![0.0, 1.0];
        assert!(cosine_similarity(&x, &y).abs() < 1e-10);
    }

    #[test]
    fn correlation_matrix_diagonal_ones() {
        let data = vec![
            1.0, 2.0, 3.0, 4.0,
            5.0, 6.0, 7.0, 8.0,
            9.0, 10.0, 11.0, 12.0,
        ];
        let corr = correlation_matrix(&data, 3, 4);
        assert!((corr[0] - 1.0).abs() < 1e-10); // [0,0]
        assert!((corr[4] - 1.0).abs() < 1e-10); // [1,1]
        assert!((corr[8] - 1.0).abs() < 1e-10); // [2,2]
        // Symmetric
        assert!((corr[1] - corr[3]).abs() < 1e-10); // [0,1] == [1,0]
    }

    // ─── d-prime Tests ───────────────────────────────────────────────────────

    #[test]
    fn d_prime_equal_rates() {
        // When hit_rate == false_alarm_rate, d' = 0
        assert!(d_prime(0.5, 0.5).abs() < 1e-6);
    }

    #[test]
    fn d_prime_good_detection() {
        // High hit rate, low FA → positive d'
        let dp = d_prime(0.9, 0.1);
        assert!(dp > 2.0);
    }

    #[test]
    fn d_prime_clamped_extremes() {
        // Rates of 0 and 1 should not panic
        let dp = d_prime(1.0, 0.0);
        assert!(dp > 0.0);
    }

    // ─── Bootstrap CI Tests ──────────────────────────────────────────────────

    #[test]
    fn bootstrap_ci_contains_mean() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let result = bootstrap_ci_mean(&data, 95.0, 1000, 12345);
        let m = mean(&data);
        assert!((result.statistic - m).abs() < 1e-10);
        assert!(result.ci_lower <= m);
        assert!(result.ci_upper >= m);
    }

    #[test]
    fn bootstrap_ci_narrow_with_constant() {
        let data = vec![5.0; 50];
        let result = bootstrap_ci_mean(&data, 95.0, 500, 42);
        assert!((result.ci_lower - 5.0).abs() < 1e-10);
        assert!((result.ci_upper - 5.0).abs() < 1e-10);
    }

    // ─── T-test Tests ────────────────────────────────────────────────────────

    #[test]
    fn t_test_one_sample_zero_mean() {
        // Data centered at 0 → not significant
        let data = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
        let result = t_test_one_sample(&data, 0.0);
        assert!(result.t_statistic.abs() < 1e-10);
        assert!((result.df - 4.0).abs() < 1e-10);
    }

    #[test]
    fn t_test_one_sample_shifted() {
        // Data far from 0 → should have large |t|
        let data = vec![100.0, 101.0, 102.0, 103.0, 104.0];
        let result = t_test_one_sample(&data, 0.0);
        assert!(result.t_statistic > 10.0);
        assert!(result.p_value < 0.01);
    }

    #[test]
    fn t_test_two_sample_same() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let b = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let result = t_test_two_sample(&a, &b);
        assert!(result.t_statistic.abs() < 1e-10);
    }

    #[test]
    fn t_test_two_sample_different() {
        let a = vec![10.0, 11.0, 12.0, 13.0, 14.0, 15.0];
        let b = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let result = t_test_two_sample(&a, &b);
        assert!(result.t_statistic > 3.0);
        assert!(result.p_value < 0.05);
    }

    // ─── Dispersion Measures Tests ───────────────────────────────────────────

    #[test]
    fn cv_positive_data() {
        let data = vec![10.0, 10.0, 10.0];
        assert!(coefficient_of_variation(&data).abs() < 1e-10);
    }

    #[test]
    fn cv_varying_data() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let cv = coefficient_of_variation(&data);
        // std/mean for {1..5}: std = sqrt(2), mean = 3 → ~0.4714
        assert!(cv > 0.0 && cv < 1.0);
    }

    #[test]
    fn fano_factor_poisson_like() {
        // For Poisson: variance ≈ mean, so F ≈ 1
        // Constant data: variance = 0, so F = 0
        let data = vec![5.0, 5.0, 5.0];
        assert!(fano_factor(&data).abs() < 1e-10);
    }

    #[test]
    fn fano_factor_positive() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!(fano_factor(&data) > 0.0);
    }

    // ─── Mutual information tests ────────────────────────────────────────────

    #[test]
    fn mutual_information_identical_signals() {
        // MI of a signal with itself should be equal to its entropy
        let x: Vec<f64> = (0..200).map(|i| (i % 10) as f64).collect();
        let result = mutual_information(&x, &x, 10);
        assert!(result.mi > 0.0);
        // MI(X;X) = H(X), so mi ≈ entropy_x
        assert!((result.mi - result.entropy_x).abs() < 0.5);
    }

    #[test]
    fn mutual_information_independent_signals() {
        // Two unrelated signals should have lower MI than correlated ones
        let x: Vec<f64> = (0..1000).map(|i| (i as f64 * 0.1).sin()).collect();
        let y: Vec<f64> = (0..1000).map(|i| ((i as f64 + 137.0) * 0.37).cos()).collect();
        let result_indep = mutual_information(&x, &y, 20);
        // Compare with perfectly correlated
        let y_corr: Vec<f64> = x.iter().map(|&v| v * 2.0 + 1.0).collect();
        let result_corr = mutual_information(&x, &y_corr, 20);
        assert!(result_indep.mi < result_corr.mi);
    }

    #[test]
    fn mutual_information_auto_bins() {
        let x: Vec<f64> = (0..300).map(|i| i as f64 / 30.0).collect();
        let y: Vec<f64> = x.iter().map(|&v| v * 2.0 + 1.0).collect();
        let result = mutual_information(&x, &y, 0); // auto bins
        assert!(result.mi > 0.0);
        assert!(result.num_bins >= 2);
    }

    #[test]
    fn mutual_information_z_correlated_high() {
        // Highly correlated → z-score should be large
        let x: Vec<f64> = (0..200).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|&v| v * 3.0 + 5.0).collect();
        let z = mutual_information_z(&x, &y, 15, 100, 42);
        assert!(z > 2.0); // significantly above null
    }
}
