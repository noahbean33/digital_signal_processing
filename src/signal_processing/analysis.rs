/// Moving average filter with a sliding window.
///
/// Output length = `input.len() - window_size + 1`.
#[must_use]
pub fn moving_average(input: &[f64], window_size: usize) -> Vec<f64> {
    if window_size == 0 || input.len() < window_size {
        return Vec::new();
    }
    let mut output = Vec::with_capacity(input.len() - window_size + 1);
    let mut sum: f64 = input[..window_size].iter().sum();
    output.push(sum / window_size as f64);
    for i in window_size..input.len() {
        sum += input[i] - input[i - window_size];
        output.push(sum / window_size as f64);
    }
    output
}

/// Finite difference (first derivative approximation).
///
/// `interval` is the sampling interval Δt.
#[must_use]
pub fn finite_difference(signal: &[f64], interval: f64) -> Vec<f64> {
    if signal.len() < 2 || interval.abs() < 1e-30 {
        return Vec::new();
    }
    signal
        .windows(2)
        .map(|w| (w[1] - w[0]) / interval)
        .collect()
}

/// Count zero crossings in a signal.
#[must_use]
pub fn count_zero_crossings(signal: &[f64]) -> usize {
    let mut count = 0;
    for i in 1..signal.len() {
        if (signal[i - 1] < 0.0 && signal[i] >= 0.0)
            || (signal[i - 1] > 0.0 && signal[i] <= 0.0)
        {
            count += 1;
        }
    }
    count
}

/// Adaptive peak detection: finds local maxima exceeding an adaptive
/// threshold of `mean + threshold_multiplier * stddev`.
#[must_use]
pub fn detect_adaptive_peaks(signal: &[f64], threshold_multiplier: f64) -> Vec<usize> {
    if signal.len() < 3 {
        return Vec::new();
    }
    let n = signal.len() as f64;
    let mean: f64 = signal.iter().sum::<f64>() / n;
    let variance: f64 = signal.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / n;
    let stddev = variance.sqrt();
    let adaptive_threshold = mean + threshold_multiplier * stddev;

    let mut peaks = Vec::new();
    for i in 1..signal.len() - 1 {
        if signal[i] > signal[i - 1]
            && signal[i] > signal[i + 1]
            && signal[i] > adaptive_threshold
        {
            peaks.push(i);
        }
    }
    peaks
}

/// Segment events based on an amplitude threshold.
///
/// Returns a list of `(start, end)` index pairs.
#[must_use]
pub fn segment_events(signal: &[f64], amplitude_threshold: f64) -> Vec<(usize, usize)> {
    let mut segments = Vec::new();
    let mut active = false;
    let mut start = 0;

    for (i, &v) in signal.iter().enumerate() {
        if !active && v.abs() > amplitude_threshold {
            active = true;
            start = i;
        } else if active && v.abs() <= amplitude_threshold {
            segments.push((start, i));
            active = false;
        }
    }
    if active {
        segments.push((start, signal.len() - 1));
    }
    segments
}

/// Sliding-window median filter for impulse noise suppression.
///
/// `window_size` must be odd.
#[must_use]
pub fn median_filter(signal: &[f64], window_size: usize) -> Vec<f64> {
    if window_size % 2 == 0 || window_size == 0 || signal.len() < window_size {
        return Vec::new();
    }
    let half = window_size / 2;
    let mut filtered = Vec::with_capacity(signal.len());
    for i in 0..signal.len() {
        let mut window = Vec::with_capacity(window_size);
        for j in 0..window_size {
            let idx = (i as isize + j as isize - half as isize)
                .max(0)
                .min(signal.len() as isize - 1) as usize;
            window.push(signal[idx]);
        }
        window.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        filtered.push(window[half]);
    }
    filtered
}

/// Linear interpolation at a fractional index.
#[must_use]
pub fn linear_interpolation(signal: &[f64], fractional_index: f64) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    let index = fractional_index as isize;
    if index < 0 {
        return signal[0];
    }
    if index + 1 >= signal.len() as isize {
        return *signal.last().unwrap_or(&0.0);
    }
    let alpha = fractional_index - index as f64;
    let i = index as usize;
    (1.0 - alpha) * signal[i] + alpha * signal[i + 1]
}

/// Estimate the synchronisation offset (delay) between two signals using
/// cross-correlation.
#[must_use]
pub fn synchronization_offset(stream1: &[f64], stream2: &[f64], max_delay: usize) -> isize {
    if stream1.is_empty() || stream2.is_empty() {
        return 0;
    }
    let length = stream1.len();
    let mut best_offset: isize = 0;
    let mut max_correlation = f64::NEG_INFINITY;

    let max_delay = max_delay as isize;
    for delay in -max_delay..=max_delay {
        let mut correlation = 0.0;
        for i in 0..length {
            let j = i as isize + delay;
            if j >= 0 && (j as usize) < length {
                correlation += stream1[i] * stream2[j as usize];
            }
        }
        if correlation > max_correlation {
            max_correlation = correlation;
            best_offset = delay;
        }
    }
    best_offset
}

/// Compute the signal energy: Σ x[n]².
#[must_use]
pub fn signal_energy(signal: &[f64]) -> f64 {
    signal.iter().map(|&x| x * x).sum()
}

/// Compute the RMS (root mean square) of a signal.
#[must_use]
pub fn rms(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    (signal_energy(signal) / signal.len() as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn moving_average_constant_signal() {
        let signal = vec![5.0; 10];
        let avg = moving_average(&signal, 3);
        assert!(avg.iter().all(|&v| (v - 5.0).abs() < 1e-10));
    }

    #[test]
    fn zero_crossings_sine() {
        let signal: Vec<f64> = (0..100)
            .map(|i| (2.0 * std::f64::consts::PI * 5.0 * i as f64 / 100.0).sin())
            .collect();
        let crossings = count_zero_crossings(&signal);
        assert!(crossings >= 8 && crossings <= 12);
    }

    #[test]
    fn median_filter_removes_spike() {
        let mut signal = vec![1.0; 11];
        signal[5] = 100.0;
        let filtered = median_filter(&signal, 3);
        assert!((filtered[5] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn linear_interp_midpoint() {
        let signal = vec![0.0, 10.0];
        assert!((linear_interpolation(&signal, 0.5) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn sync_offset_identity() {
        let signal = vec![0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0];
        let offset = synchronization_offset(&signal, &signal, 3);
        assert_eq!(offset, 0);
    }
}
