use std::f64::consts::PI;

/// A complex number with double-precision real and imaginary parts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    #[must_use]
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    #[must_use]
    pub fn zero() -> Self {
        Self { re: 0.0, im: 0.0 }
    }

    #[must_use]
    pub fn from_polar(r: f64, theta: f64) -> Self {
        Self {
            re: r * theta.cos(),
            im: r * theta.sin(),
        }
    }

    #[must_use]
    pub fn norm(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }

    #[must_use]
    pub fn norm_sqr(&self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    #[must_use]
    pub fn conj(&self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    #[must_use]
    pub fn arg(&self) -> f64 {
        self.im.atan2(self.re)
    }

    #[must_use]
    pub fn exp(self) -> Self {
        let e = self.re.exp();
        Self {
            re: e * self.im.cos(),
            im: e * self.im.sin(),
        }
    }
}

impl std::ops::Add for Complex {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl std::ops::Sub for Complex {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl std::ops::Mul for Complex {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl std::ops::Div for Complex {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let denom = rhs.norm_sqr();
        if denom < 1e-30 {
            return Self::zero();
        }
        Self {
            re: (self.re * rhs.re + self.im * rhs.im) / denom,
            im: (self.im * rhs.re - self.re * rhs.im) / denom,
        }
    }
}

impl std::ops::Mul<f64> for Complex {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self {
            re: self.re * rhs,
            im: self.im * rhs,
        }
    }
}

impl std::ops::AddAssign for Complex {
    fn add_assign(&mut self, rhs: Self) {
        self.re += rhs.re;
        self.im += rhs.im;
    }
}

impl std::ops::MulAssign for Complex {
    fn mul_assign(&mut self, rhs: Self) {
        let re = self.re * rhs.re - self.im * rhs.im;
        let im = self.re * rhs.im + self.im * rhs.re;
        self.re = re;
        self.im = im;
    }
}

/// Compute the twiddle factor W_N^k = exp(-2πik/N).
#[must_use]
pub fn twiddle_factor(k: usize, n: usize) -> Complex {
    if n == 0 {
        return Complex::new(1.0, 0.0);
    }
    let angle = -2.0 * PI * k as f64 / n as f64;
    Complex::new(angle.cos(), angle.sin())
}

/// Perform bit-reversal permutation on `data` in-place.
pub fn bit_reversal_permutation(data: &mut [Complex]) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    let mut n_bits = 0u32;
    while (1usize << n_bits) < n {
        n_bits += 1;
    }
    for i in 0..n {
        let mut j = 0usize;
        for k in 0..n_bits {
            if i & (1 << k) != 0 {
                j |= 1 << (n_bits - 1 - k);
            }
        }
        if j > i {
            data.swap(i, j);
        }
    }
}

/// In-place iterative Cooley-Tukey radix-2 FFT.
///
/// `data` length **must** be a power of 2.
pub fn fft(data: &mut [Complex]) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    bit_reversal_permutation(data);

    let mut len = 2;
    while len <= n {
        let angle = -2.0 * PI / len as f64;
        let wlen = Complex::new(angle.cos(), angle.sin());
        let mut i = 0;
        while i < n {
            let mut w = Complex::new(1.0, 0.0);
            for j in 0..len / 2 {
                let u = data[i + j];
                let v = data[i + j + len / 2] * w;
                data[i + j] = u + v;
                data[i + j + len / 2] = u - v;
                w *= wlen;
            }
            i += len;
        }
        len <<= 1;
    }
}

/// In-place inverse FFT. Divides by N after the transform.
pub fn ifft(data: &mut [Complex]) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    // Conjugate input
    for c in data.iter_mut() {
        c.im = -c.im;
    }
    fft(data);
    // Conjugate and scale
    let scale = 1.0 / n as f64;
    for c in data.iter_mut() {
        c.re *= scale;
        c.im = -c.im * scale;
    }
}

/// Convert a real-valued signal to a complex vector (imaginary parts zero).
#[must_use]
pub fn real_to_complex(signal: &[f64]) -> Vec<Complex> {
    signal.iter().map(|&v| Complex::new(v, 0.0)).collect()
}

/// Zero-pad a complex vector to the next power of 2.
#[must_use]
pub fn zero_pad_to_power_of_2(data: &[Complex]) -> Vec<Complex> {
    let n = data.len();
    let mut size = 1;
    while size < n {
        size <<= 1;
    }
    let mut padded = Vec::with_capacity(size);
    padded.extend_from_slice(data);
    padded.resize(size, Complex::zero());
    padded
}

/// Compute the magnitude spectrum from a complex FFT result.
#[must_use]
pub fn magnitude_spectrum(fft_result: &[Complex]) -> Vec<f64> {
    fft_result.iter().map(Complex::norm).collect()
}

/// Compute the power spectrum (|X[k]|² / N).
#[must_use]
pub fn power_spectrum(fft_result: &[Complex]) -> Vec<f64> {
    let n = fft_result.len() as f64;
    fft_result.iter().map(|c| c.norm_sqr() / n).collect()
}

/// Compute the phase spectrum in radians.
#[must_use]
pub fn phase_spectrum(fft_result: &[Complex]) -> Vec<f64> {
    fft_result.iter().map(Complex::arg).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fft_of_dc_signal() {
        let mut data = vec![Complex::new(1.0, 0.0); 4];
        fft(&mut data);
        assert!((data[0].re - 4.0).abs() < 1e-10);
        for c in &data[1..] {
            assert!(c.norm() < 1e-10);
        }
    }

    #[test]
    fn fft_ifft_roundtrip() {
        let original = vec![
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
            Complex::new(4.0, 0.0),
        ];
        let mut data = original.clone();
        fft(&mut data);
        ifft(&mut data);
        for (a, b) in original.iter().zip(data.iter()) {
            assert!((a.re - b.re).abs() < 1e-10);
            assert!((a.im - b.im).abs() < 1e-10);
        }
    }

    #[test]
    fn bit_reversal_identity_for_len_1() {
        let mut data = vec![Complex::new(42.0, 0.0)];
        bit_reversal_permutation(&mut data);
        assert!((data[0].re - 42.0).abs() < 1e-10);
    }

    #[test]
    fn twiddle_w0_is_one() {
        let w = twiddle_factor(0, 8);
        assert!((w.re - 1.0).abs() < 1e-10);
        assert!(w.im.abs() < 1e-10);
    }
}
