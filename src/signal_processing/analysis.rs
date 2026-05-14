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

/// Centered (symmetric) moving average filter.
///
/// For each sample, averages `window_size` points centered on that sample.
/// `window_size` must be odd. Samples near edges use clamped (replicated) boundary.
#[must_use]
pub fn centered_moving_average(signal: &[f64], window_size: usize) -> Vec<f64> {
    if window_size == 0 || window_size % 2 == 0 || signal.is_empty() {
        return Vec::new();
    }
    let half = window_size / 2;
    let n = signal.len();
    let mut output = Vec::with_capacity(n);
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..window_size {
            let idx_signed = i as isize + j as isize - half as isize;
            let idx = idx_signed.max(0).min(n as isize - 1) as usize;
            sum += signal[idx];
        }
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

// ─── DC Removal Filter ───────────────────────────────────────────────────────

/// Remove DC offset from a signal using a single-pole high-pass IIR filter.
///
/// Implements: `y[n] = x[n] - x[n-1] + alpha * y[n-1]`
/// where `alpha = (1 - sin(2π * cutoff_hz / fs)) / cos(2π * cutoff_hz / fs)`.
///
/// `cutoff_hz` – cutoff frequency in Hz (typically very low, e.g. 5–20 Hz).
/// `sample_rate` – sampling frequency in Hz.
#[must_use]
pub fn dc_remove(signal: &[f64], cutoff_hz: f64, sample_rate: f64) -> Vec<f64> {
    if signal.is_empty() || sample_rate <= 0.0 || cutoff_hz <= 0.0 {
        return signal.to_vec();
    }
    let omega = 2.0 * std::f64::consts::PI * cutoff_hz / sample_rate;
    let alpha = (1.0 - omega.sin()) / omega.cos();

    let mut output = vec![0.0; signal.len()];
    output[0] = signal[0];
    for n in 1..signal.len() {
        output[n] = signal[n] - signal[n - 1] + alpha * output[n - 1];
    }
    output
}

// ─── Envelope Detection ───────────────────────────────────────────────────────

/// Envelope follower with configurable attack and decay coefficients.
///
/// Tracks the amplitude envelope of a signal using a one-pole filter:
/// - When `|x[n]| > env[n-1]`: `env[n] = attack * |x[n]| + (1 - attack) * env[n-1]`
/// - Otherwise:                 `env[n] = decay  * |x[n]| + (1 - decay)  * env[n-1]`
///
/// `attack` and `decay` are in the range `(0, 1]`. Larger values track faster.
#[must_use]
pub fn envelope(signal: &[f64], attack: f64, decay: f64) -> Vec<f64> {
    if signal.is_empty() {
        return Vec::new();
    }
    let mut env = Vec::with_capacity(signal.len());
    let mut state = 0.0_f64;
    for &x in signal {
        let abs_x = x.abs();
        if abs_x > state {
            state = attack * abs_x + (1.0 - attack) * state;
        } else {
            state = decay * abs_x + (1.0 - decay) * state;
        }
        env.push(state);
    }
    env
}

/// RMS envelope follower.
///
/// Computes a smoothed RMS envelope using a one-pole filter on `x²`:
/// `rms²[n] = α * x[n]² + (1 - α) * rms²[n-1]`,  output = √(rms²[n]).
///
/// `alpha` controls the smoothing (higher = faster tracking).
#[must_use]
pub fn envelope_rms(signal: &[f64], alpha: f64) -> Vec<f64> {
    if signal.is_empty() {
        return Vec::new();
    }
    let mut env = Vec::with_capacity(signal.len());
    let mut state = 0.0_f64;
    for &x in signal {
        state = alpha * x * x + (1.0 - alpha) * state;
        env.push(state.sqrt());
    }
    env
}

/// Hilbert envelope (analytic signal magnitude).
///
/// Computes the instantaneous amplitude envelope by:
/// 1. Taking the FFT of the signal.
/// 2. Zeroing the negative-frequency components (creating the analytic signal).
/// 3. Taking the inverse FFT.
/// 4. Returning the magnitude of the complex analytic signal.
///
/// The signal is zero-padded to the next power of 2.
#[must_use]
pub fn envelope_hilbert(signal: &[f64]) -> Vec<f64> {
    use super::fft::{Complex, fft, ifft, real_to_complex, zero_pad_to_power_of_2};

    if signal.is_empty() {
        return Vec::new();
    }
    let orig_len = signal.len();
    let mut data = zero_pad_to_power_of_2(&real_to_complex(signal));
    let n = data.len();
    fft(&mut data);

    // Build analytic signal: keep DC and Nyquist, double positive freqs, zero negative
    // Bin 0 (DC) stays, bins 1..N/2-1 doubled, bin N/2 stays, bins N/2+1..N-1 zeroed
    if n > 1 {
        for i in 1..n / 2 {
            data[i] = data[i] * 2.0;
        }
        for i in (n / 2 + 1)..n {
            data[i] = Complex::zero();
        }
    }

    ifft(&mut data);

    data[..orig_len].iter().map(|c| c.norm()).collect()
}

// ─── Signal Integration ───────────────────────────────────────────────────────

/// Numerical integration (cumulative trapezoidal rule).
///
/// `interval` is the sampling interval Δt.
/// `y[n] = y[n-1] + Δt * (x[n] + x[n-1]) / 2`
#[must_use]
pub fn integrate(signal: &[f64], interval: f64) -> Vec<f64> {
    if signal.is_empty() {
        return Vec::new();
    }
    let mut output = Vec::with_capacity(signal.len());
    output.push(0.0);
    for i in 1..signal.len() {
        let prev = output[i - 1];
        output.push(prev + interval * (signal[i] + signal[i - 1]) / 2.0);
    }
    output
}

// ─── Threshold / Clamp / Limit ────────────────────────────────────────────────

/// Threshold mode for signal processing utilities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ThresholdMode {
    /// Zero values below the threshold.
    AboveZero,
    /// Zero values above the threshold.
    BelowZero,
    /// Replace values below the threshold with the threshold.
    AboveClamp,
    /// Replace values above the threshold with the threshold.
    BelowClamp,
}

/// Apply a threshold to a signal.
///
/// Depending on the mode, values that do not meet the threshold condition
/// are either zeroed out or clamped.
#[must_use]
pub fn threshold(signal: &[f64], thresh: f64, mode: ThresholdMode) -> Vec<f64> {
    signal
        .iter()
        .map(|&x| match mode {
            ThresholdMode::AboveZero => {
                if x >= thresh { x } else { 0.0 }
            }
            ThresholdMode::BelowZero => {
                if x <= thresh { x } else { 0.0 }
            }
            ThresholdMode::AboveClamp => {
                if x >= thresh { x } else { thresh }
            }
            ThresholdMode::BelowClamp => {
                if x <= thresh { x } else { thresh }
            }
        })
        .collect()
}

/// Soft threshold: shrink values towards zero.
///
/// `y = sign(x) * max(|x| - thresh, 0)`
#[must_use]
pub fn soft_threshold(signal: &[f64], thresh: f64) -> Vec<f64> {
    signal
        .iter()
        .map(|&x| {
            let abs = x.abs();
            if abs > thresh {
                x.signum() * (abs - thresh)
            } else {
                0.0
            }
        })
        .collect()
}

/// Clamp (limit) signal values to `[min_val, max_val]`.
#[must_use]
pub fn clamp(signal: &[f64], min_val: f64, max_val: f64) -> Vec<f64> {
    signal
        .iter()
        .map(|&x| x.max(min_val).min(max_val))
        .collect()
}

/// Count how many values exceed the threshold.
#[must_use]
pub fn count_over_threshold(signal: &[f64], thresh: f64) -> usize {
    signal.iter().filter(|&&x| x > thresh).count()
}

/// Count how many absolute values exceed the threshold.
#[must_use]
pub fn count_abs_over_threshold(signal: &[f64], thresh: f64) -> usize {
    signal.iter().filter(|&&x| x.abs() > thresh).count()
}

/// Element-wise maximum of two signals.
#[must_use]
pub fn select_max(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x.max(y)).collect()
}

/// Element-wise minimum of two signals.
#[must_use]
pub fn select_min(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x.min(y)).collect()
}

// ─── Teager-Kaiser Energy Operator (TKEO) ────────────────────────────────────

/// Teager-Kaiser Energy Operator for instantaneous energy estimation.
///
/// `tkeo[n] = x[n]² − x[n−1] · x[n+1]`
///
/// Output length = `signal.len()`.  The first and last samples are set to zero
/// because the operator requires one neighbour on each side.
#[must_use]
pub fn tkeo(signal: &[f64]) -> Vec<f64> {
    let n = signal.len();
    if n < 3 {
        return vec![0.0; n];
    }
    let mut out = vec![0.0; n];
    for i in 1..n - 1 {
        out[i] = signal[i] * signal[i] - signal[i - 1] * signal[i + 1];
    }
    out
}

// ─── Gaussian Smoothing Kernel ───────────────────────────────────────────────

/// Generate a normalised Gaussian smoothing kernel parameterised by FWHM.
///
/// * `fwhm` – full-width at half-maximum **in samples** (e.g. `fwhm_seconds * sample_rate`).
/// * `half_len` – number of samples on each side of the centre (kernel length = `2 * half_len + 1`).
///
/// The kernel sums to 1.0 so that convolution preserves signal amplitude.
#[must_use]
pub fn gaussian_kernel(fwhm: f64, half_len: usize) -> Vec<f64> {
    let len = 2 * half_len + 1;
    let mut kernel = Vec::with_capacity(len);
    for i in 0..len {
        let t = i as f64 - half_len as f64;
        kernel.push((-4.0 * 2.0_f64.ln() * t * t / (fwhm * fwhm)).exp());
    }
    let sum: f64 = kernel.iter().sum();
    if sum.abs() > 1e-30 {
        for v in &mut kernel {
            *v /= sum;
        }
    }
    kernel
}

/// Gaussian-smooth a signal in the time domain.
///
/// Applies a FWHM-parameterised Gaussian kernel via direct convolution.
/// Edge samples outside the kernel radius are left unchanged (copied from the original).
///
/// * `fwhm` – full-width at half-maximum **in samples**.
/// * `half_len` – number of samples on each side of the centre.
#[must_use]
pub fn gaussian_smooth(signal: &[f64], fwhm: f64, half_len: usize) -> Vec<f64> {
    let n = signal.len();
    if n == 0 {
        return Vec::new();
    }
    let kernel = gaussian_kernel(fwhm, half_len);
    let k = half_len;
    let mut out = signal.to_vec();
    for i in k..n.saturating_sub(k) {
        let mut acc = 0.0;
        for (j, &w) in kernel.iter().enumerate() {
            acc += signal[i + j - k] * w;
        }
        out[i] = acc;
    }
    out
}

// ─── Windowed Variance / RMS ─────────────────────────────────────────────────

/// Sliding-window variance.
///
/// For each sample, computes the variance of the surrounding `2 * half_win + 1` points.
/// Boundary samples use a smaller window (clamped to signal edges).
#[must_use]
pub fn windowed_variance(signal: &[f64], half_win: usize) -> Vec<f64> {
    let n = signal.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let lo = if i >= half_win { i - half_win } else { 0 };
        let hi = (i + half_win + 1).min(n);
        let seg = &signal[lo..hi];
        let count = seg.len() as f64;
        let mean = seg.iter().sum::<f64>() / count;
        let var = seg.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / count;
        out.push(var);
    }
    out
}

/// Sliding-window RMS (root mean square).
///
/// For each sample, computes the RMS of the surrounding `2 * half_win + 1` points.
/// Boundary samples use a smaller window.
#[must_use]
pub fn windowed_rms(signal: &[f64], half_win: usize) -> Vec<f64> {
    let n = signal.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let lo = if i >= half_win { i - half_win } else { 0 };
        let hi = (i + half_win + 1).min(n);
        let seg = &signal[lo..hi];
        let count = seg.len() as f64;
        let rms_val = (seg.iter().map(|&x| x * x).sum::<f64>() / count).sqrt();
        out.push(rms_val);
    }
    out
}

// ─── Windowed SNR ────────────────────────────────────────────────────────────

/// Sliding-window Signal-to-Noise Ratio: `mean / std` in each window.
///
/// Returns the raw (linear) SNR at each sample.  Convert to dB via `10 * log10(snr)`.
/// Windows at the boundaries are clamped.
#[must_use]
pub fn windowed_snr(signal: &[f64], half_win: usize) -> Vec<f64> {
    let n = signal.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let lo = if i >= half_win { i - half_win } else { 0 };
        let hi = (i + half_win + 1).min(n);
        let seg = &signal[lo..hi];
        let count = seg.len() as f64;
        let mean = seg.iter().sum::<f64>() / count;
        let var = seg.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / count;
        let std = var.sqrt();
        if std.abs() < 1e-30 {
            out.push(0.0);
        } else {
            out.push(mean / std);
        }
    }
    out
}

// ─── Dynamic Time Warping ────────────────────────────────────────────────────

/// Result of Dynamic Time Warping.
#[derive(Clone, Debug)]
pub struct DtwResult {
    /// Total accumulated distance along the optimal warping path.
    pub distance: f64,
    /// Warping path as `(index_in_x, index_in_y)` pairs, from start to end.
    pub path: Vec<(usize, usize)>,
}

/// Dynamic Time Warping between two sequences.
///
/// Computes the optimal alignment between `x` and `y` by minimising the
/// cumulative absolute-difference cost.  Returns the total distance and the
/// warping path.
#[must_use]
pub fn dtw(x: &[f64], y: &[f64]) -> DtwResult {
    let nx = x.len();
    let ny = y.len();
    if nx == 0 || ny == 0 {
        return DtwResult {
            distance: 0.0,
            path: Vec::new(),
        };
    }

    // Build cost matrix
    let mut dm = vec![vec![f64::INFINITY; ny]; nx];
    dm[0][0] = (x[0] - y[0]).abs();

    for i in 1..nx {
        dm[i][0] = (x[i] - y[0]).abs() + dm[i - 1][0];
    }
    for j in 1..ny {
        dm[0][j] = (x[0] - y[j]).abs() + dm[0][j - 1];
    }
    for i in 1..nx {
        for j in 1..ny {
            let cost = (x[i] - y[j]).abs();
            dm[i][j] = cost + dm[i - 1][j].min(dm[i][j - 1]).min(dm[i - 1][j - 1]);
        }
    }

    let distance = dm[nx - 1][ny - 1];

    // Trace-back
    let mut path = Vec::new();
    let mut i = nx - 1;
    let mut j = ny - 1;
    path.push((i, j));
    while i > 0 || j > 0 {
        if i == 0 {
            j -= 1;
        } else if j == 0 {
            i -= 1;
        } else {
            let diag = dm[i - 1][j - 1];
            let left = dm[i][j - 1];
            let up = dm[i - 1][j];
            if diag <= left && diag <= up {
                i -= 1;
                j -= 1;
            } else if left < up {
                j -= 1;
            } else {
                i -= 1;
            }
        }
        path.push((i, j));
    }
    path.reverse();
    DtwResult { distance, path }
}

// ─── FWHM Measurement ───────────────────────────────────────────────────────

/// Measure the Full-Width at Half-Maximum of the tallest peak in a signal.
///
/// Normalises the signal to its peak value, then locates the half-max crossings
/// on either side.  Returns the width in **samples** (fractional).
///
/// Returns `None` if the signal is too short or the peak is at an edge.
#[must_use]
pub fn measure_fwhm(signal: &[f64]) -> Option<f64> {
    if signal.len() < 3 {
        return None;
    }

    // Find peak
    let (peak_idx, &peak_val) = signal
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())?;

    if peak_val.abs() < 1e-30 {
        return None;
    }

    let half = peak_val / 2.0;

    // Pre-peak half-max crossing (scan left from peak)
    let mut pre = None;
    for i in (1..=peak_idx).rev() {
        if signal[i - 1] <= half && signal[i] >= half {
            // Linear interpolation
            let frac = (half - signal[i - 1]) / (signal[i] - signal[i - 1]);
            pre = Some((i - 1) as f64 + frac);
            break;
        }
    }

    // Post-peak half-max crossing (scan right from peak)
    let mut post = None;
    for i in peak_idx..signal.len() - 1 {
        if signal[i] >= half && signal[i + 1] <= half {
            let frac = (half - signal[i]) / (signal[i + 1] - signal[i]);
            post = Some(i as f64 + frac);
            break;
        }
    }

    match (pre, post) {
        (Some(p), Some(q)) => Some(q - p),
        _ => None,
    }
}

// ─── Area Under the Curve ────────────────────────────────────────────────────

/// A contiguous lobe (above-zero region) of a signal.
#[derive(Clone, Debug)]
pub struct SignalLobe {
    /// Start index (inclusive).
    pub start: usize,
    /// End index (exclusive).
    pub end: usize,
    /// Area under the curve (trapezoidal integration, or `sum * dt`).
    pub area: f64,
}

/// Compute the area under each positive lobe of a signal.
///
/// A lobe is a contiguous region where the signal is > 0.  The area is
/// computed as the sum of sample values times `dt` (the sampling interval).
///
/// * `dt` – sampling interval (1.0 / sample_rate).
#[must_use]
pub fn signal_lobe_areas(signal: &[f64], dt: f64) -> Vec<SignalLobe> {
    let n = signal.len();
    let mut lobes = Vec::new();
    let mut in_lobe = false;
    let mut start = 0;

    for i in 0..n {
        if !in_lobe && signal[i] > 0.0 {
            in_lobe = true;
            start = i;
        } else if in_lobe && signal[i] <= 0.0 {
            let area: f64 = signal[start..i].iter().sum::<f64>() * dt;
            lobes.push(SignalLobe {
                start,
                end: i,
                area,
            });
            in_lobe = false;
        }
    }
    if in_lobe {
        let area: f64 = signal[start..n].iter().sum::<f64>() * dt;
        lobes.push(SignalLobe {
            start,
            end: n,
            area,
        });
    }
    lobes
}

// ─── Polynomial Detrending ───────────────────────────────────────────────────

/// Fit a polynomial of a given `order` to the signal and return the coefficients.
///
/// Uses a least-squares Vandermonde approach.  Coefficients are returned in
/// descending power order: `[a_order, a_{order-1}, ..., a_1, a_0]`.
#[must_use]
pub fn polyfit(signal: &[f64], order: usize) -> Vec<f64> {
    let n = signal.len();
    if n == 0 || order >= n {
        return Vec::new();
    }
    let p = order + 1;
    // Build normal equations: V^T V c = V^T y
    let mut vtv = vec![0.0; p * p];
    let mut vty = vec![0.0; p];

    for i in 0..n {
        let t = i as f64;
        let mut powers = vec![1.0; p];
        for k in 1..p {
            powers[k] = powers[k - 1] * t;
        }
        for r in 0..p {
            for c in 0..p {
                vtv[r * p + c] += powers[r] * powers[c];
            }
            vty[r] += powers[r] * signal[i];
        }
    }

    // Solve via Gaussian elimination with partial pivoting
    let mut aug = vec![vec![0.0; p + 1]; p];
    for r in 0..p {
        for c in 0..p {
            aug[r][c] = vtv[r * p + c];
        }
        aug[r][p] = vty[r];
    }

    for col in 0..p {
        // Pivot
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in col + 1..p {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        aug.swap(col, max_row);

        let diag = aug[col][col];
        if diag.abs() < 1e-30 {
            return vec![0.0; p];
        }
        for c in col..=p {
            aug[col][c] /= diag;
        }
        for row in 0..p {
            if row == col {
                continue;
            }
            let factor = aug[row][col];
            for c in col..=p {
                aug[row][c] -= factor * aug[col][c];
            }
        }
    }

    // Coefficients in ascending power; reverse to descending
    let mut coeffs: Vec<f64> = (0..p).map(|r| aug[r][p]).collect();
    coeffs.reverse();
    coeffs
}

/// Evaluate a polynomial (coefficients in descending power order) at a point.
#[must_use]
pub fn polyval(coeffs: &[f64], x: f64) -> f64 {
    coeffs.iter().fold(0.0, |acc, &c| acc * x + c)
}

/// Remove a polynomial trend from a signal.
///
/// Returns the residual signal after subtracting the best-fit polynomial.
#[must_use]
pub fn polynomial_detrend(signal: &[f64], order: usize) -> Vec<f64> {
    let coeffs = polyfit(signal, order);
    if coeffs.is_empty() {
        return signal.to_vec();
    }
    (0..signal.len())
        .map(|i| signal[i] - polyval(&coeffs, i as f64))
        .collect()
}

/// Select the best polynomial detrending order using the Bayesian Information
/// Criterion (BIC).
///
/// Evaluates orders in `min_order..=max_order` and returns the order with the
/// smallest BIC value.
#[must_use]
pub fn best_detrend_order(signal: &[f64], min_order: usize, max_order: usize) -> usize {
    let n = signal.len();
    if n < 2 {
        return min_order;
    }
    let ln_n = (n as f64).ln();
    let n_f = n as f64;

    let mut best_order = min_order;
    let mut best_bic = f64::INFINITY;

    for order in min_order..=max_order {
        let coeffs = polyfit(signal, order);
        if coeffs.is_empty() {
            continue;
        }
        let sse: f64 = (0..n)
            .map(|i| {
                let residual = signal[i] - polyval(&coeffs, i as f64);
                residual * residual
            })
            .sum::<f64>()
            / n_f;

        if sse < 1e-30 {
            return order;
        }
        let bic = n_f * sse.ln() + (order as f64) * ln_n;
        if bic < best_bic {
            best_bic = bic;
            best_order = order;
        }
    }
    best_order
}

// ─── Template Projection (Artifact Removal) ──────────────────────────────────

/// Remove an artifact signal from a data signal via least-squares projection.
///
/// For each column (trial), fits `data = β₀ + β₁ · artifact` and returns the
/// residual.  This is the standard regression-based artifact removal used in
/// EEG/EMG processing (e.g. removing eye-movement artifacts).
///
/// * `data` – the contaminated signal.
/// * `artifact` – the reference artifact signal (same length as `data`).
///
/// Returns the cleaned (residual) signal.
#[must_use]
pub fn template_projection(data: &[f64], artifact: &[f64]) -> Vec<f64> {
    let n = data.len();
    if n == 0 || artifact.len() != n {
        return data.to_vec();
    }

    // X = [1, artifact]  → 2-column design matrix
    // Normal equations: (X^T X) b = X^T y
    let mut s1 = 0.0; // Σ artifact
    let mut s2 = 0.0; // Σ artifact²
    let mut sy = 0.0; // Σ data
    let mut say = 0.0; // Σ artifact * data
    let n_f = n as f64;

    for i in 0..n {
        let a = artifact[i];
        let d = data[i];
        s1 += a;
        s2 += a * a;
        sy += d;
        say += a * d;
    }

    // [n,  s1 ] [b0]   [sy ]
    // [s1, s2 ] [b1] = [say]
    let det = n_f * s2 - s1 * s1;
    if det.abs() < 1e-30 {
        return data.to_vec();
    }

    let b0 = (s2 * sy - s1 * say) / det;
    let b1 = (n_f * say - s1 * sy) / det;

    (0..n)
        .map(|i| data[i] - (b0 + b1 * artifact[i]))
        .collect()
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

    #[test]
    fn centered_moving_average_constant() {
        let signal = vec![5.0; 10];
        let avg = centered_moving_average(&signal, 3);
        assert_eq!(avg.len(), 10);
        assert!(avg.iter().all(|&v| (v - 5.0).abs() < 1e-10));
    }

    #[test]
    fn centered_moving_average_smooths_spike() {
        let mut signal = vec![0.0; 7];
        signal[3] = 7.0;
        let avg = centered_moving_average(&signal, 3);
        // Centre point: (0 + 7 + 0) / 3
        assert!((avg[3] - 7.0 / 3.0).abs() < 1e-10);
        // Neighbours: (0 + 0 + 7) / 3  and  (7 + 0 + 0) / 3
        assert!((avg[2] - 7.0 / 3.0).abs() < 1e-10);
        assert!((avg[4] - 7.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn centered_moving_average_even_window_returns_empty() {
        let signal = vec![1.0, 2.0, 3.0];
        assert!(centered_moving_average(&signal, 2).is_empty());
    }

    #[test]
    fn dc_remove_removes_offset() {
        // Signal with DC offset of 5.0
        let signal: Vec<f64> = (0..1000)
            .map(|i| 5.0 + (2.0 * std::f64::consts::PI * 50.0 * i as f64 / 1000.0).sin())
            .collect();
        let filtered = dc_remove(&signal, 10.0, 1000.0);
        // After settling, mean should be near zero
        let tail_mean: f64 = filtered[500..].iter().sum::<f64>() / 500.0;
        assert!(tail_mean.abs() < 0.5);
    }

    #[test]
    fn dc_remove_preserves_ac() {
        let signal: Vec<f64> = (0..1000)
            .map(|i| (2.0 * std::f64::consts::PI * 100.0 * i as f64 / 1000.0).sin())
            .collect();
        let filtered = dc_remove(&signal, 5.0, 1000.0);
        // AC energy should be mostly preserved
        let in_energy: f64 = signal[200..].iter().map(|x| x * x).sum();
        let out_energy: f64 = filtered[200..].iter().map(|x| x * x).sum();
        assert!((out_energy / in_energy - 1.0).abs() < 0.1);
    }

    #[test]
    fn dc_remove_empty_passthrough() {
        let empty: Vec<f64> = Vec::new();
        assert!(dc_remove(&empty, 10.0, 1000.0).is_empty());
    }

    #[test]
    fn envelope_tracks_amplitude() {
        let pi = std::f64::consts::PI;
        // AM signal: carrier modulated by a slow envelope
        let signal: Vec<f64> = (0..1000)
            .map(|i| {
                let t = i as f64 / 1000.0;
                let carrier = (2.0 * pi * 100.0 * t).sin();
                let mod_env = 0.5 + 0.5 * (2.0 * pi * 5.0 * t).sin();
                carrier * mod_env
            })
            .collect();
        let env = envelope(&signal, 0.1, 0.01);
        assert_eq!(env.len(), signal.len());
        // Envelope should be non-negative
        assert!(env.iter().all(|&v| v >= -1e-10));
    }

    #[test]
    fn envelope_rms_positive() {
        let pi = std::f64::consts::PI;
        let signal: Vec<f64> = (0..200)
            .map(|i| (2.0 * pi * 50.0 * i as f64 / 1000.0).sin())
            .collect();
        let env = envelope_rms(&signal, 0.05);
        assert_eq!(env.len(), signal.len());
        assert!(env.iter().all(|&v| v >= 0.0));
        // RMS should settle near 1/√2 ≈ 0.707 for a unit sine
        let tail_mean: f64 = env[100..].iter().sum::<f64>() / env[100..].len() as f64;
        assert!((tail_mean - 0.707).abs() < 0.2);
    }

    #[test]
    fn envelope_hilbert_of_sine() {
        let pi = std::f64::consts::PI;
        let n = 256;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * pi * 10.0 * i as f64 / n as f64).sin())
            .collect();
        let env = envelope_hilbert(&signal);
        assert_eq!(env.len(), n);
        // Hilbert envelope of a pure sine should be ~constant (~1.0)
        // Check the middle portion (edges have transient effects)
        let mid = &env[n / 4..3 * n / 4];
        let mean: f64 = mid.iter().sum::<f64>() / mid.len() as f64;
        assert!((mean - 1.0).abs() < 0.15);
    }

    #[test]
    fn integrate_constant() {
        // Integral of constant 2.0 with dt=0.1 over 10 samples
        let signal = vec![2.0; 10];
        let result = integrate(&signal, 0.1);
        assert_eq!(result.len(), 10);
        assert!((result[0] - 0.0).abs() < 1e-10);
        // After 9 intervals: 9 * 0.1 * 2.0 = 1.8
        assert!((result[9] - 1.8).abs() < 1e-10);
    }

    #[test]
    fn integrate_empty() {
        assert!(integrate(&[], 1.0).is_empty());
    }

    #[test]
    fn threshold_above_zero() {
        let signal = vec![-2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
        let out = threshold(&signal, 1.0, ThresholdMode::AboveZero);
        assert_eq!(out, vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn threshold_below_clamp() {
        let signal = vec![-2.0, 0.0, 5.0, 10.0];
        let out = threshold(&signal, 5.0, ThresholdMode::BelowClamp);
        assert_eq!(out, vec![-2.0, 0.0, 5.0, 5.0]);
    }

    #[test]
    fn soft_threshold_shrinks() {
        let signal = vec![-3.0, -1.0, 0.0, 1.0, 3.0];
        let out = soft_threshold(&signal, 2.0);
        assert!((out[0] - (-1.0)).abs() < 1e-10);
        assert!((out[1] - 0.0).abs() < 1e-10);
        assert!((out[2] - 0.0).abs() < 1e-10);
        assert!((out[3] - 0.0).abs() < 1e-10);
        assert!((out[4] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn clamp_limits_values() {
        let signal = vec![-5.0, -1.0, 0.0, 1.0, 5.0];
        let out = clamp(&signal, -2.0, 2.0);
        assert_eq!(out, vec![-2.0, -1.0, 0.0, 1.0, 2.0]);
    }

    #[test]
    fn count_over_threshold_works() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(count_over_threshold(&signal, 3.0), 2);
        assert_eq!(count_abs_over_threshold(&signal, 3.0), 2);
    }

    #[test]
    fn select_max_min() {
        let a = vec![1.0, 5.0, 3.0];
        let b = vec![4.0, 2.0, 6.0];
        assert_eq!(select_max(&a, &b), vec![4.0, 5.0, 6.0]);
        assert_eq!(select_min(&a, &b), vec![1.0, 2.0, 3.0]);
    }

    // ─── TKEO tests ──────────────────────────────────────────────────────────

    #[test]
    fn tkeo_pure_sine() {
        use std::f64::consts::PI;
        let n = 200;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * 10.0 * i as f64 / n as f64).sin())
            .collect();
        let energy = tkeo(&signal);
        assert_eq!(energy.len(), n);
        // First and last are zero
        assert!((energy[0]).abs() < 1e-10);
        assert!((energy[n - 1]).abs() < 1e-10);
        // Interior values should be positive for a sine
        assert!(energy[50] > 0.0);
    }

    #[test]
    fn tkeo_constant_is_zero() {
        let signal = vec![5.0; 20];
        let energy = tkeo(&signal);
        for &e in &energy[1..energy.len() - 1] {
            assert!(e.abs() < 1e-10);
        }
    }

    #[test]
    fn tkeo_short_signal() {
        assert_eq!(tkeo(&[1.0, 2.0]), vec![0.0, 0.0]);
        assert_eq!(tkeo(&[]).len(), 0);
    }

    // ─── Gaussian kernel / smooth tests ──────────────────────────────────────

    #[test]
    fn gaussian_kernel_sums_to_one() {
        let k = gaussian_kernel(10.0, 30);
        let sum: f64 = k.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
        assert_eq!(k.len(), 61);
    }

    #[test]
    fn gaussian_kernel_peak_at_center() {
        let k = gaussian_kernel(5.0, 15);
        let max_idx = k
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;
        assert_eq!(max_idx, 15);
    }

    #[test]
    fn gaussian_smooth_constant() {
        let signal = vec![3.0; 100];
        let smoothed = gaussian_smooth(&signal, 10.0, 20);
        // Constant signal should remain constant after smoothing
        for &v in &smoothed[20..80] {
            assert!((v - 3.0).abs() < 1e-10);
        }
    }

    #[test]
    fn gaussian_smooth_reduces_noise() {
        // Simple test: spike should be attenuated
        let mut signal = vec![0.0; 100];
        signal[50] = 100.0;
        let smoothed = gaussian_smooth(&signal, 10.0, 20);
        assert!(smoothed[50] < 100.0);
    }

    // ─── Windowed variance / RMS tests ───────────────────────────────────────

    #[test]
    fn windowed_variance_constant_is_zero() {
        let signal = vec![7.0; 50];
        let var = windowed_variance(&signal, 5);
        for &v in &var {
            assert!(v.abs() < 1e-10);
        }
    }

    #[test]
    fn windowed_variance_detects_variability() {
        let mut signal = vec![0.0; 100];
        // Add high-variance region
        for i in 40..60 {
            signal[i] = if i % 2 == 0 { 10.0 } else { -10.0 };
        }
        let var = windowed_variance(&signal, 5);
        // Variance in the noisy region should be much higher than quiet region
        let quiet_var = var[10];
        let noisy_var = var[50];
        assert!(noisy_var > quiet_var + 1.0);
    }

    #[test]
    fn windowed_rms_constant() {
        let signal = vec![3.0; 50];
        let r = windowed_rms(&signal, 5);
        for &v in &r {
            assert!((v - 3.0).abs() < 1e-10);
        }
    }

    #[test]
    fn windowed_rms_length() {
        let signal = vec![1.0; 20];
        assert_eq!(windowed_rms(&signal, 3).len(), 20);
    }

    // ─── Windowed SNR tests ──────────────────────────────────────────────────

    #[test]
    fn windowed_snr_constant_is_zero() {
        // Constant signal → std = 0 → SNR = 0 (guard)
        let signal = vec![5.0; 50];
        let snr = windowed_snr(&signal, 5);
        for &v in &snr {
            assert!(v.abs() < 1e-10);
        }
    }

    #[test]
    fn windowed_snr_positive_mean_positive_snr() {
        // Signal with positive mean and some noise
        let signal: Vec<f64> = (0..100)
            .map(|i| 10.0 + (i as f64 * 0.3).sin())
            .collect();
        let snr = windowed_snr(&signal, 10);
        // With a large positive mean, SNR should be positive
        let mid_snr = snr[50];
        assert!(mid_snr > 0.0);
    }

    // ─── DTW tests ───────────────────────────────────────────────────────────

    #[test]
    fn dtw_identical_signals() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let result = dtw(&x, &x);
        assert!((result.distance).abs() < 1e-10);
        assert_eq!(result.path.len(), 5);
        for (i, &(a, b)) in result.path.iter().enumerate() {
            assert_eq!(a, i);
            assert_eq!(b, i);
        }
    }

    #[test]
    fn dtw_different_lengths() {
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y = vec![1.0, 3.0, 4.0];
        let result = dtw(&x, &y);
        assert!(result.distance >= 0.0);
        assert!(!result.path.is_empty());
        // Path starts at (0,0) and ends at (3,2)
        assert_eq!(result.path[0], (0, 0));
        assert_eq!(*result.path.last().unwrap(), (3, 2));
    }

    #[test]
    fn dtw_empty() {
        let result = dtw(&[], &[1.0]);
        assert!(result.path.is_empty());
    }

    // ─── FWHM tests ─────────────────────────────────────────────────────────

    #[test]
    fn fwhm_gaussian_peak() {
        // Create a Gaussian with known FWHM
        let fwhm_expected = 20.0;
        let center = 50;
        let signal: Vec<f64> = (0..100)
            .map(|i| {
                let t = i as f64 - center as f64;
                (-4.0 * 2.0_f64.ln() * t * t / (fwhm_expected * fwhm_expected)).exp()
            })
            .collect();
        let fwhm = measure_fwhm(&signal);
        assert!(fwhm.is_some());
        assert!((fwhm.unwrap() - fwhm_expected).abs() < 1.0);
    }

    #[test]
    fn fwhm_flat_signal_is_none() {
        let signal = vec![0.0; 50];
        assert!(measure_fwhm(&signal).is_none());
    }

    // ─── AUC / signal lobe tests ─────────────────────────────────────────────

    #[test]
    fn signal_lobe_areas_basic() {
        let signal = vec![-1.0, 1.0, 2.0, 3.0, -1.0, 1.0, 1.0, -1.0];
        let lobes = signal_lobe_areas(&signal, 1.0);
        assert_eq!(lobes.len(), 2);
        // First lobe: indices 1..4, values [1,2,3], area = 6
        assert_eq!(lobes[0].start, 1);
        assert_eq!(lobes[0].end, 4);
        assert!((lobes[0].area - 6.0).abs() < 1e-10);
        // Second lobe: indices 5..7, values [1,1], area = 2
        assert_eq!(lobes[1].start, 5);
        assert_eq!(lobes[1].end, 7);
        assert!((lobes[1].area - 2.0).abs() < 1e-10);
    }

    #[test]
    fn signal_lobe_areas_all_negative() {
        let signal = vec![-1.0, -2.0, -3.0];
        let lobes = signal_lobe_areas(&signal, 0.001);
        assert!(lobes.is_empty());
    }

    // ─── Polynomial detrending tests ─────────────────────────────────────────

    #[test]
    fn polyfit_linear() {
        // y = 2x + 1
        let signal: Vec<f64> = (0..50).map(|i| 2.0 * i as f64 + 1.0).collect();
        let coeffs = polyfit(&signal, 1);
        // Descending: [a1, a0] = [2.0, 1.0]
        assert_eq!(coeffs.len(), 2);
        assert!((coeffs[0] - 2.0).abs() < 1e-6);
        assert!((coeffs[1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn polyval_quadratic() {
        // p(x) = x² + 2x + 3  →  coeffs = [1, 2, 3]
        let coeffs = vec![1.0, 2.0, 3.0];
        assert!((polyval(&coeffs, 0.0) - 3.0).abs() < 1e-10);
        assert!((polyval(&coeffs, 1.0) - 6.0).abs() < 1e-10);
        assert!((polyval(&coeffs, 2.0) - 11.0).abs() < 1e-10);
    }

    #[test]
    fn polynomial_detrend_removes_linear_trend() {
        let signal: Vec<f64> = (0..100)
            .map(|i| 0.5 * i as f64 + (i as f64 * 0.3).sin())
            .collect();
        let detrended = polynomial_detrend(&signal, 1);
        // Mean of detrended should be near zero
        let mean: f64 = detrended.iter().sum::<f64>() / detrended.len() as f64;
        assert!(mean.abs() < 1.0);
    }

    #[test]
    fn best_detrend_order_picks_linear_for_linear() {
        let signal: Vec<f64> = (0..200).map(|i| 3.0 * i as f64 - 10.0).collect();
        let order = best_detrend_order(&signal, 1, 10);
        // For a perfectly linear signal, any order >= 1 gives SSE ≈ 0,
        // but BIC penalises higher orders → should pick a low order
        assert!(order <= 3);
    }

    // ─── Template projection tests ───────────────────────────────────────────

    #[test]
    fn template_projection_removes_artifact() {
        let n = 200;
        // Clean signal
        let clean: Vec<f64> = (0..n)
            .map(|i| (i as f64 * 0.1).sin())
            .collect();
        // Artifact
        let artifact: Vec<f64> = (0..n)
            .map(|i| 5.0 * (i as f64 * 0.02).cos())
            .collect();
        // Contaminated = clean + 3 * artifact
        let contaminated: Vec<f64> = clean
            .iter()
            .zip(artifact.iter())
            .map(|(&c, &a)| c + 3.0 * a)
            .collect();

        let residual = template_projection(&contaminated, &artifact);
        assert_eq!(residual.len(), n);

        // Residual should be closer to clean than the contaminated signal
        let err_before: f64 = contaminated
            .iter()
            .zip(clean.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>();
        let err_after: f64 = residual
            .iter()
            .zip(clean.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>();
        assert!(err_after < err_before * 0.1);
    }

    #[test]
    fn template_projection_mismatched_lengths() {
        let data = vec![1.0, 2.0, 3.0];
        let artifact = vec![1.0, 2.0];
        let result = template_projection(&data, &artifact);
        assert_eq!(result, data); // Returns original when lengths mismatch
    }
}
