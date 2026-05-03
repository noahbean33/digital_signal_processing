use std::f64::consts::PI;

use super::fft::Complex;

/// Compute the sample covariance matrix from a data matrix.
///
/// `data` is stored in row-major order: `data[row * cols + col]`.
/// `rows` = number of samples (snapshots), `cols` = number of sensors.
///
/// Returns the covariance matrix (cols × cols, row-major).
#[must_use]
pub fn covariance_matrix(data: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    if rows <= 1 || cols == 0 || data.len() != rows * cols {
        return vec![0.0; cols * cols];
    }
    // Compute column means
    let mut mean = vec![0.0; cols];
    for r in 0..rows {
        for c in 0..cols {
            mean[c] += data[r * cols + c];
        }
    }
    for m in &mut mean {
        *m /= rows as f64;
    }

    // Compute covariance
    let mut cov = vec![0.0; cols * cols];
    for r in 0..rows {
        for i in 0..cols {
            let di = data[r * cols + i] - mean[i];
            for j in 0..cols {
                let dj = data[r * cols + j] - mean[j];
                cov[i * cols + j] += di * dj;
            }
        }
    }
    let scale = 1.0 / (rows - 1) as f64;
    for v in &mut cov {
        *v *= scale;
    }
    cov
}

/// Compute a regularised covariance matrix: C + α·I.
#[must_use]
pub fn regularized_covariance(data: &[f64], rows: usize, cols: usize, alpha: f64) -> Vec<f64> {
    let mut cov = covariance_matrix(data, rows, cols);
    for i in 0..cols {
        cov[i * cols + i] += alpha;
    }
    cov
}

/// Compute a steering vector for a Uniform Linear Array (ULA).
///
/// `frequency` – normalised spatial frequency.
/// `array_size` – number of sensors.
#[must_use]
pub fn steering_vector(frequency: f64, array_size: usize) -> Vec<Complex> {
    (0..array_size)
        .map(|i| {
            let angle = -2.0 * PI * frequency * i as f64;
            Complex::new(angle.cos(), angle.sin())
        })
        .collect()
}

/// Compute the MUSIC pseudo-spectrum value at a given steering vector.
///
/// `noise_subspace` is stored column-major: `noise_subspace[row * noise_cols + col]`.
/// `noise_rows` = array size, `noise_cols` = number of noise eigenvectors.
#[must_use]
pub fn music_pseudospectrum(
    noise_subspace: &[Complex],
    noise_rows: usize,
    noise_cols: usize,
    steering: &[Complex],
) -> f64 {
    if steering.len() != noise_rows || noise_subspace.len() != noise_rows * noise_cols {
        return 0.0;
    }
    // Compute P = E_n * E_n^H, then a^H * P * a
    // More efficiently: || E_n^H * a ||^2
    let mut sum_sq = 0.0;
    for col in 0..noise_cols {
        let mut dot = Complex::zero();
        for row in 0..noise_rows {
            // E_n^H[col, row] = conj(E_n[row, col])
            let en_conj = noise_subspace[row * noise_cols + col].conj();
            dot += en_conj * steering[row];
        }
        sum_sq += dot.norm_sqr();
    }
    if sum_sq < 1e-9 {
        return 0.0;
    }
    1.0 / sum_sq
}

/// Scan the MUSIC pseudo-spectrum over a range of normalised frequencies.
///
/// Returns `(frequencies, spectrum_values)`.
#[must_use]
pub fn music_spectrum_scan(
    noise_subspace: &[Complex],
    noise_rows: usize,
    noise_cols: usize,
    num_points: usize,
) -> (Vec<f64>, Vec<f64>) {
    let mut frequencies = Vec::with_capacity(num_points);
    let mut spectrum = Vec::with_capacity(num_points);
    for i in 0..num_points {
        let freq = i as f64 / num_points as f64 - 0.5;
        let sv = steering_vector(freq, noise_rows);
        let val = music_pseudospectrum(noise_subspace, noise_rows, noise_cols, &sv);
        frequencies.push(freq);
        spectrum.push(val);
    }
    (frequencies, spectrum)
}

/// Simple eigenvalue decomposition for a real symmetric matrix using the
/// Jacobi eigenvalue algorithm.
///
/// `matrix` is n×n row-major.
///
/// Returns `(eigenvalues, eigenvectors_column_major)`.
#[must_use]
pub fn symmetric_eigen(matrix: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let max_iterations = 100 * n * n;
    let mut a = matrix.to_vec();
    // Eigenvector matrix (starts as identity)
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }

    for _ in 0..max_iterations {
        // Find largest off-diagonal element
        let mut max_val = 0.0_f64;
        let mut p = 0;
        let mut q = 1;
        for i in 0..n {
            for j in (i + 1)..n {
                let val = a[i * n + j].abs();
                if val > max_val {
                    max_val = val;
                    p = i;
                    q = j;
                }
            }
        }
        if max_val < 1e-12 {
            break;
        }

        // Compute rotation
        let app = a[p * n + p];
        let aqq = a[q * n + q];
        let apq = a[p * n + q];
        let theta = if (app - aqq).abs() < 1e-30 {
            PI / 4.0
        } else {
            0.5 * (2.0 * apq / (app - aqq)).atan()
        };
        let cos_t = theta.cos();
        let sin_t = theta.sin();

        // Apply Jacobi rotation
        let mut new_a = a.clone();
        for i in 0..n {
            if i != p && i != q {
                let aip = a[i * n + p];
                let aiq = a[i * n + q];
                new_a[i * n + p] = cos_t * aip + sin_t * aiq;
                new_a[p * n + i] = new_a[i * n + p];
                new_a[i * n + q] = -sin_t * aip + cos_t * aiq;
                new_a[q * n + i] = new_a[i * n + q];
            }
        }
        new_a[p * n + p] = cos_t * cos_t * app + 2.0 * cos_t * sin_t * apq + sin_t * sin_t * aqq;
        new_a[q * n + q] = sin_t * sin_t * app - 2.0 * cos_t * sin_t * apq + cos_t * cos_t * aqq;
        new_a[p * n + q] = 0.0;
        new_a[q * n + p] = 0.0;
        a = new_a;

        // Update eigenvectors
        for i in 0..n {
            let vip = v[i * n + p];
            let viq = v[i * n + q];
            v[i * n + p] = cos_t * vip + sin_t * viq;
            v[i * n + q] = -sin_t * vip + cos_t * viq;
        }
    }

    let eigenvalues: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
    (eigenvalues, v)
}

/// Extract noise subspace from a covariance matrix.
///
/// `num_signals` – number of signals (sources); the remaining eigenvectors
/// form the noise subspace.
///
/// Returns the noise subspace as a `Complex` vector (column-major, rows = array_size,
/// cols = array_size - num_signals).
#[must_use]
pub fn extract_noise_subspace(
    cov_matrix: &[f64],
    array_size: usize,
    num_signals: usize,
) -> Vec<Complex> {
    let (eigenvalues, eigenvectors) = symmetric_eigen(cov_matrix, array_size);

    // Sort eigenvalues (ascending) and get indices
    let mut indices: Vec<usize> = (0..array_size).collect();
    indices.sort_by(|&a, &b| {
        eigenvalues[a]
            .partial_cmp(&eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let noise_dim = array_size.saturating_sub(num_signals);
    let mut noise_sub = vec![Complex::zero(); array_size * noise_dim];
    for (col, &idx) in indices.iter().take(noise_dim).enumerate() {
        for row in 0..array_size {
            noise_sub[row * noise_dim + col] = Complex::new(eigenvectors[row * array_size + idx], 0.0);
        }
    }
    noise_sub
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn covariance_of_constant_is_zero() {
        let data = vec![1.0, 2.0, 1.0, 2.0, 1.0, 2.0];
        let cov = covariance_matrix(&data, 3, 2);
        assert!((cov[0]).abs() < 1e-10); // var of column 0
        assert!((cov[3]).abs() < 1e-10); // var of column 1
    }

    #[test]
    fn steering_vector_first_element_is_one() {
        let sv = steering_vector(0.1, 4);
        assert!((sv[0].re - 1.0).abs() < 1e-10);
        assert!(sv[0].im.abs() < 1e-10);
    }

    #[test]
    fn eigen_identity_gives_ones() {
        let mat = vec![1.0, 0.0, 0.0, 1.0];
        let (vals, _) = symmetric_eigen(&mat, 2);
        assert!((vals[0] - 1.0).abs() < 1e-10);
        assert!((vals[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn eigen_diagonal_matrix() {
        let mat = vec![3.0, 0.0, 0.0, 5.0];
        let (mut vals, _) = symmetric_eigen(&mat, 2);
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((vals[0] - 3.0).abs() < 1e-10);
        assert!((vals[1] - 5.0).abs() < 1e-10);
    }
}
