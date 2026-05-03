/// Discrete linear convolution: y[n] = Σ x[m] * h[n-m].
///
/// Output length = `signal.len() + kernel.len() - 1`.
#[must_use]
pub fn convolve(signal: &[f64], kernel: &[f64]) -> Vec<f64> {
    if signal.is_empty() || kernel.is_empty() {
        return Vec::new();
    }
    let n = signal.len();
    let k = kernel.len();
    let conv_size = n + k - 1;
    let mut result = vec![0.0; conv_size];
    for i in 0..conv_size {
        for j in 0..k {
            let si = i as isize - j as isize;
            if si >= 0 && (si as usize) < n {
                result[i] += signal[si as usize] * kernel[j];
            }
        }
    }
    result
}

/// Basic (valid) correlation without reversing the kernel.
///
/// Output length = `signal.len() - kernel.len() + 1`.
#[must_use]
pub fn correlate(signal: &[f64], kernel: &[f64]) -> Vec<f64> {
    if signal.is_empty() || kernel.is_empty() || signal.len() < kernel.len() {
        return Vec::new();
    }
    let n = signal.len();
    let k = kernel.len();
    let corr_size = n - k + 1;
    let mut result = vec![0.0; corr_size];
    for i in 0..corr_size {
        let mut sum = 0.0;
        for j in 0..k {
            sum += signal[i + j] * kernel[j];
        }
        result[i] = sum;
    }
    result
}

/// Full cross-correlation between two signals.
///
/// Output length = `x.len() + y.len() - 1`.
#[must_use]
pub fn cross_correlate(x: &[f64], y: &[f64]) -> Vec<f64> {
    if x.is_empty() || y.is_empty() {
        return Vec::new();
    }
    let n = x.len();
    let m = y.len();
    let corr_size = n + m - 1;
    let mut result = vec![0.0; corr_size];
    for lag in 0..corr_size {
        let mut sum = 0.0;
        for i in 0..n {
            let j = lag as isize - i as isize;
            if j >= 0 && (j as usize) < m {
                sum += x[i] * y[j as usize];
            }
        }
        result[lag] = sum;
    }
    result
}

/// Overlap-add method for efficient convolution of long signals.
#[must_use]
pub fn overlap_add(signal: &[f64], kernel: &[f64], block_size: usize) -> Vec<f64> {
    if signal.is_empty() || kernel.is_empty() || block_size == 0 {
        return Vec::new();
    }
    let signal_size = signal.len();
    let kernel_size = kernel.len();
    let output_size = signal_size + kernel_size - 1;
    let mut output = vec![0.0; output_size];

    let mut start = 0;
    while start < signal_size {
        let current_block_size = block_size.min(signal_size - start);
        let block = &signal[start..start + current_block_size];
        for i in 0..current_block_size {
            for j in 0..kernel_size {
                output[start + i + j] += block[i] * kernel[j];
            }
        }
        start += block_size;
    }
    output
}

/// Normalised cross-correlation (Pearson correlation coefficient).
///
/// Returns a single scalar for equal-length signals.
#[must_use]
pub fn normalized_cross_correlation(signal1: &[f64], signal2: &[f64]) -> f64 {
    if signal1.len() != signal2.len() || signal1.is_empty() {
        return 0.0;
    }
    let n = signal1.len() as f64;
    let mean1: f64 = signal1.iter().sum::<f64>() / n;
    let mean2: f64 = signal2.iter().sum::<f64>() / n;

    let mut num = 0.0;
    let mut den1 = 0.0;
    let mut den2 = 0.0;
    for i in 0..signal1.len() {
        let d1 = signal1[i] - mean1;
        let d2 = signal2[i] - mean2;
        num += d1 * d2;
        den1 += d1 * d1;
        den2 += d2 * d2;
    }
    if den1 < 1e-30 || den2 < 1e-30 {
        return 0.0;
    }
    num / (den1 * den2).sqrt()
}

/// Block convolution of a signal segment with an impulse response.
#[must_use]
pub fn block_convolution(block: &[f64], impulse_response: &[f64]) -> Vec<f64> {
    convolve(block, impulse_response)
}

/// Autocorrelation: cross-correlation of a signal with itself.
///
/// Output length = `2 * signal.len() - 1`.
#[must_use]
pub fn autocorrelate(signal: &[f64]) -> Vec<f64> {
    cross_correlate(signal, signal)
}

/// Running sum (cumulative integration).
///
/// `y[n] = Σ_{k=0}^{n} x[k]`
#[must_use]
pub fn running_sum(signal: &[f64]) -> Vec<f64> {
    if signal.is_empty() {
        return Vec::new();
    }
    let mut output = Vec::with_capacity(signal.len());
    let mut sum = 0.0;
    for &x in signal {
        sum += x;
        output.push(sum);
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn convolve_with_impulse() {
        let signal = vec![1.0, 2.0, 3.0];
        let kernel = vec![1.0];
        let result = convolve(&signal, &kernel);
        assert_eq!(result.len(), 3);
        assert!((result[0] - 1.0).abs() < 1e-10);
        assert!((result[1] - 2.0).abs() < 1e-10);
        assert!((result[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn overlap_add_matches_direct() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let kernel = vec![0.2, 0.5, 0.2];
        let direct = convolve(&signal, &kernel);
        let ola = overlap_add(&signal, &kernel, 3);
        assert_eq!(direct.len(), ola.len());
        for (a, b) in direct.iter().zip(ola.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn normalized_cross_correlation_identical() {
        let s = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let ncc = normalized_cross_correlation(&s, &s);
        assert!((ncc - 1.0).abs() < 1e-10);
    }

    #[test]
    fn autocorrelate_peak_at_centre() {
        let s = vec![1.0, 2.0, 3.0, 2.0, 1.0];
        let ac = autocorrelate(&s);
        assert_eq!(ac.len(), 9);
        // Peak is at the centre (lag 0)
        let centre = ac.len() / 2;
        for (i, &v) in ac.iter().enumerate() {
            assert!(v <= ac[centre] + 1e-10, "lag {i} exceeded centre");
        }
    }

    #[test]
    fn running_sum_correct() {
        let s = vec![1.0, 2.0, 3.0, 4.0];
        let rs = running_sum(&s);
        assert_eq!(rs, vec![1.0, 3.0, 6.0, 10.0]);
    }

    #[test]
    fn running_sum_empty() {
        assert!(running_sum(&[]).is_empty());
    }
}
