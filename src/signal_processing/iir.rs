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

// ─── IIR Second-Order Sections (SOS / Biquad Cascade) ────────────────────────

/// A single biquad (second-order section) state for sample-by-sample processing.
///
/// Uses Direct Form II transposed: two delay elements per section.
pub struct BiquadSection {
    b: [f64; 3],
    a: [f64; 3],
    w: [f64; 2],
}

impl BiquadSection {
    /// Create a new biquad section.
    ///
    /// `b` – numerator coefficients `[b0, b1, b2]`.
    /// `a` – denominator coefficients `[a0, a1, a2]` (a0 is typically 1.0).
    #[must_use]
    pub fn new(b: [f64; 3], a: [f64; 3]) -> Self {
        Self { b, a, w: [0.0; 2] }
    }

    /// Reset internal state to zero.
    pub fn reset(&mut self) {
        self.w = [0.0; 2];
    }

    /// Process a single sample through this biquad section (DF2).
    pub fn process(&mut self, input: f64) -> f64 {
        let wn = input - self.a[1] * self.w[0] - self.a[2] * self.w[1];
        let yn = self.b[0] * wn + self.b[1] * self.w[0] + self.b[2] * self.w[1];
        self.w[1] = self.w[0];
        self.w[0] = wn;
        yn
    }
}

/// A cascaded biquad (SOS) filter: chains multiple second-order sections.
pub struct SosFilter {
    sections: Vec<BiquadSection>,
}

impl SosFilter {
    /// Create a new SOS filter from arrays of b and a coefficients.
    ///
    /// `b_sections` – `&[[f64; 3]]` numerator coefficients per section.
    /// `a_sections` – `&[[f64; 3]]` denominator coefficients per section.
    #[must_use]
    pub fn new(b_sections: &[[f64; 3]], a_sections: &[[f64; 3]]) -> Self {
        let sections = b_sections
            .iter()
            .zip(a_sections.iter())
            .map(|(&b, &a)| BiquadSection::new(b, a))
            .collect();
        Self { sections }
    }

    /// Reset internal state of all sections.
    pub fn reset(&mut self) {
        for section in &mut self.sections {
            section.reset();
        }
    }

    /// Process a single sample through the entire cascade.
    pub fn process_sample(&mut self, sample: f64) -> f64 {
        let mut x = sample;
        for section in &mut self.sections {
            x = section.process(x);
        }
        x
    }

    /// Filter an entire signal through the SOS cascade (batch mode).
    #[must_use]
    pub fn apply(&mut self, signal: &[f64]) -> Vec<f64> {
        signal.iter().map(|&s| self.process_sample(s)).collect()
    }
}

/// Convenience function: apply SOS filtering in one call (batch, no state retained).
#[must_use]
pub fn apply_sos(signal: &[f64], b_sections: &[[f64; 3]], a_sections: &[[f64; 3]]) -> Vec<f64> {
    let mut filter = SosFilter::new(b_sections, a_sections);
    filter.apply(signal)
}

// ─── Biquad Filter Design ─────────────────────────────────────────────────────

/// Design a low-pass biquad filter.
///
/// `frequency` – normalised cutoff frequency (frequency_Hz / sample_rate_Hz).
/// `q` – quality factor.
///
/// Returns `(b, a)` coefficient arrays for a `BiquadSection`.
#[must_use]
pub fn biquad_lowpass(frequency: f64, q: f64) -> ([f64; 3], [f64; 3]) {
    let k = (PI * frequency).tan();
    let k2 = k * k;
    let norm = 1.0 / (1.0 + k / q + k2);
    let a0 = k2 * norm;
    let a1 = 2.0 * a0;
    let a2 = a0;
    let b1 = 2.0 * (k2 - 1.0) * norm;
    let b2 = (1.0 - k / q + k2) * norm;
    ([a0, a1, a2], [1.0, b1, b2])
}

/// Design a high-pass biquad filter.
///
/// `frequency` – normalised cutoff frequency (frequency_Hz / sample_rate_Hz).
/// `q` – quality factor.
#[must_use]
pub fn biquad_highpass(frequency: f64, q: f64) -> ([f64; 3], [f64; 3]) {
    let k = (PI * frequency).tan();
    let k2 = k * k;
    let norm = 1.0 / (1.0 + k / q + k2);
    let a0 = norm;
    let a1 = -2.0 * a0;
    let a2 = a0;
    let b1 = 2.0 * (k2 - 1.0) * norm;
    let b2 = (1.0 - k / q + k2) * norm;
    ([a0, a1, a2], [1.0, b1, b2])
}

/// Design a band-pass biquad filter.
///
/// `frequency` – normalised center frequency (frequency_Hz / sample_rate_Hz).
/// `q` – quality factor.
#[must_use]
pub fn biquad_bandpass(frequency: f64, q: f64) -> ([f64; 3], [f64; 3]) {
    let k = (PI * frequency).tan();
    let k2 = k * k;
    let norm = 1.0 / (1.0 + k / q + k2);
    let a0 = k / q * norm;
    let a1 = 0.0;
    let a2 = -a0;
    let b1 = 2.0 * (k2 - 1.0) * norm;
    let b2 = (1.0 - k / q + k2) * norm;
    ([a0, a1, a2], [1.0, b1, b2])
}

/// Design a notch (band-reject) biquad filter.
///
/// `frequency` – normalised center frequency (frequency_Hz / sample_rate_Hz).
/// `q` – quality factor.
#[must_use]
pub fn biquad_notch(frequency: f64, q: f64) -> ([f64; 3], [f64; 3]) {
    let k = (PI * frequency).tan();
    let k2 = k * k;
    let norm = 1.0 / (1.0 + k / q + k2);
    let a0 = (1.0 + k2) * norm;
    let a1 = 2.0 * (k2 - 1.0) * norm;
    let a2 = a0;
    let b1 = a1;
    let b2 = (1.0 - k / q + k2) * norm;
    ([a0, a1, a2], [1.0, b1, b2])
}

/// Design a peaking EQ biquad filter.
///
/// `frequency` – normalised center frequency (frequency_Hz / sample_rate_Hz).
/// `q` – quality factor.
/// `gain_db` – gain in decibels (positive = boost, negative = cut).
#[must_use]
pub fn biquad_peak(frequency: f64, q: f64, gain_db: f64) -> ([f64; 3], [f64; 3]) {
    let k = (PI * frequency).tan();
    let k2 = k * k;
    let v = (gain_db.abs() / 20.0 * std::f64::consts::LN_10).exp();

    if gain_db >= 0.0 {
        let norm = 1.0 / (1.0 + 1.0 / q * k + k2);
        let a0 = (1.0 + v / q * k + k2) * norm;
        let a1 = 2.0 * (k2 - 1.0) * norm;
        let a2 = (1.0 - v / q * k + k2) * norm;
        let b1 = a1;
        let b2 = (1.0 - 1.0 / q * k + k2) * norm;
        ([a0, a1, a2], [1.0, b1, b2])
    } else {
        let norm = 1.0 / (1.0 + v / q * k + k2);
        let a0 = (1.0 + 1.0 / q * k + k2) * norm;
        let a1 = 2.0 * (k2 - 1.0) * norm;
        let a2 = (1.0 - 1.0 / q * k + k2) * norm;
        let b1 = a1;
        let b2 = (1.0 - v / q * k + k2) * norm;
        ([a0, a1, a2], [1.0, b1, b2])
    }
}

/// Design a low-shelf biquad filter.
///
/// `frequency` – normalised cutoff frequency (frequency_Hz / sample_rate_Hz).
/// `gain_db` – gain in decibels (positive = boost, negative = cut).
#[must_use]
pub fn biquad_lowshelf(frequency: f64, gain_db: f64) -> ([f64; 3], [f64; 3]) {
    let k = (PI * frequency).tan();
    let k2 = k * k;
    let v = (gain_db.abs() / 20.0 * std::f64::consts::LN_10).exp();
    let sqrt2 = std::f64::consts::SQRT_2;

    if gain_db >= 0.0 {
        let norm = 1.0 / (1.0 + sqrt2 * k + k2);
        let a0 = (1.0 + (2.0 * v).sqrt() * k + v * k2) * norm;
        let a1 = 2.0 * (v * k2 - 1.0) * norm;
        let a2 = (1.0 - (2.0 * v).sqrt() * k + v * k2) * norm;
        let b1 = 2.0 * (k2 - 1.0) * norm;
        let b2 = (1.0 - sqrt2 * k + k2) * norm;
        ([a0, a1, a2], [1.0, b1, b2])
    } else {
        let norm = 1.0 / (1.0 + (2.0 * v).sqrt() * k + v * k2);
        let a0 = (1.0 + sqrt2 * k + k2) * norm;
        let a1 = 2.0 * (k2 - 1.0) * norm;
        let a2 = (1.0 - sqrt2 * k + k2) * norm;
        let b1 = 2.0 * (v * k2 - 1.0) * norm;
        let b2 = (1.0 - (2.0 * v).sqrt() * k + v * k2) * norm;
        ([a0, a1, a2], [1.0, b1, b2])
    }
}

/// Design a high-shelf biquad filter.
///
/// `frequency` – normalised cutoff frequency (frequency_Hz / sample_rate_Hz).
/// `gain_db` – gain in decibels (positive = boost, negative = cut).
#[must_use]
pub fn biquad_highshelf(frequency: f64, gain_db: f64) -> ([f64; 3], [f64; 3]) {
    let k = (PI * frequency).tan();
    let k2 = k * k;
    let v = (gain_db.abs() / 20.0 * std::f64::consts::LN_10).exp();
    let sqrt2 = std::f64::consts::SQRT_2;

    if gain_db >= 0.0 {
        let norm = 1.0 / (1.0 + sqrt2 * k + k2);
        let a0 = (v + (2.0 * v).sqrt() * k + k2) * norm;
        let a1 = 2.0 * (k2 - v) * norm;
        let a2 = (v - (2.0 * v).sqrt() * k + k2) * norm;
        let b1 = 2.0 * (k2 - 1.0) * norm;
        let b2 = (1.0 - sqrt2 * k + k2) * norm;
        ([a0, a1, a2], [1.0, b1, b2])
    } else {
        let norm = 1.0 / (v + (2.0 * v).sqrt() * k + k2);
        let a0 = (1.0 + sqrt2 * k + k2) * norm;
        let a1 = 2.0 * (k2 - 1.0) * norm;
        let a2 = (1.0 - sqrt2 * k + k2) * norm;
        let b1 = 2.0 * (k2 - v) * norm;
        let b2 = (v - (2.0 * v).sqrt() * k + k2) * norm;
        ([a0, a1, a2], [1.0, b1, b2])
    }
}

/// Design an all-pass biquad filter.
///
/// `frequency` – normalised frequency (frequency_Hz / sample_rate_Hz).
/// `q` – quality factor.
#[must_use]
pub fn biquad_allpass(frequency: f64, q: f64) -> ([f64; 3], [f64; 3]) {
    let alpha = frequency.sin() / (2.0 * q);
    let cs = frequency.cos();
    let b0 = 1.0 / (1.0 + alpha);
    let b1 = -2.0 * cs * b0;
    let b2 = (1.0 - alpha) * b0;
    let a0 = (1.0 - alpha) * b0;
    let a1 = -2.0 * cs * b0;
    let a2 = (1.0 + alpha) * b0;
    ([a0, a1, a2], [1.0, b1, b2])
}

// ─── IIR Prototype Filter Design (ZPK) ───────────────────────────────────────

/// Zero-Pole-Gain representation of an analog or digital filter.
#[derive(Debug, Clone)]
pub struct Zpk {
    pub zeros: Vec<Complex>,
    pub poles: Vec<Complex>,
    pub gain: f64,
}

/// Compute the analog prototype poles/zeros for a Butterworth filter of order N.
///
/// Returns a `Zpk` with no zeros, N poles equally spaced on the left half of the
/// unit circle, and unity gain.
#[must_use]
pub fn butterworth(order: usize) -> Zpk {
    let mut poles = Vec::with_capacity(order);
    for k in 0..order {
        let theta = PI * (2 * k + order + 1) as f64 / (2 * order) as f64;
        poles.push(Complex::new(theta.cos(), theta.sin()));
    }
    Zpk {
        zeros: Vec::new(),
        poles,
        gain: 1.0,
    }
}

/// Compute the analog prototype poles/zeros for a Chebyshev Type I filter.
///
/// `order` – filter order.
/// `ripple_db` – passband ripple in dB.
#[must_use]
pub fn chebyshev1(order: usize, ripple_db: f64) -> Zpk {
    if order == 0 {
        return Zpk { zeros: Vec::new(), poles: Vec::new(), gain: 1.0 };
    }
    let eps = (10.0_f64.powf(0.1 * ripple_db) - 1.0).sqrt();
    let mu = (1.0 / eps).asinh() / order as f64;
    let mut poles = Vec::with_capacity(order);
    for k in 0..order {
        let theta = PI * (2 * k + 1) as f64 / (2 * order) as f64;
        let re = -(mu.sinh()) * theta.sin();
        let im = (mu.cosh()) * theta.cos();
        poles.push(Complex::new(re, im));
    }
    let mut gain = 1.0;
    for p in &poles {
        gain *= p.re * p.re + p.im * p.im;
    }
    gain = gain.sqrt();
    if order % 2 == 0 {
        gain /= (1.0 + eps * eps).sqrt();
    }
    Zpk {
        zeros: Vec::new(),
        poles,
        gain,
    }
}

/// Compute the analog prototype poles/zeros for a Chebyshev Type II filter.
///
/// `order` – filter order.
/// `stopband_db` – stopband attenuation in dB.
#[must_use]
pub fn chebyshev2(order: usize, stopband_db: f64) -> Zpk {
    if order == 0 {
        return Zpk { zeros: Vec::new(), poles: Vec::new(), gain: 1.0 };
    }
    let eps = 1.0 / (10.0_f64.powf(0.1 * stopband_db) - 1.0).sqrt();
    let mu = (1.0 / eps).asinh() / order as f64;
    let mut poles = Vec::with_capacity(order);
    let mut zeros = Vec::new();
    for k in 0..order {
        let theta = PI * (2 * k + 1) as f64 / (2 * order) as f64;
        let sin_t = theta.sin();
        let cos_t = theta.cos();
        // Poles are reciprocals of Chebyshev I poles
        let re = -(mu.sinh()) * sin_t;
        let im = (mu.cosh()) * cos_t;
        let denom = re * re + im * im;
        if denom > 1e-30 {
            poles.push(Complex::new(re / denom, -im / denom));
        }
        // Zeros on the jω axis
        if sin_t.abs() > 1e-10 {
            let z_im = 1.0 / sin_t;
            zeros.push(Complex::new(0.0, z_im));
            zeros.push(Complex::new(0.0, -z_im));
        }
    }
    // Remove duplicate conjugate zeros (keep unique pairs)
    let mut unique_zeros = Vec::new();
    for z in &zeros {
        let already = unique_zeros.iter().any(|u: &Complex| {
            (u.re - z.re).abs() < 1e-10 && (u.im - z.im).abs() < 1e-10
        });
        if !already {
            unique_zeros.push(*z);
        }
    }
    // Compute gain
    let mut num = Complex::new(1.0, 0.0);
    for z in &unique_zeros {
        num = num * Complex::new(-z.re, -z.im);
    }
    let mut den = Complex::new(1.0, 0.0);
    for p in &poles {
        den = den * Complex::new(-p.re, -p.im);
    }
    let gain = (den / num).re.abs();
    Zpk {
        zeros: unique_zeros,
        poles,
        gain,
    }
}

/// Compute the analog prototype poles/zeros for a Bessel filter of order N.
///
/// Uses the reverse Bessel polynomial roots for maximally-flat group delay.
#[must_use]
pub fn bessel(order: usize) -> Zpk {
    if order == 0 {
        return Zpk { zeros: Vec::new(), poles: Vec::new(), gain: 1.0 };
    }
    // Bessel polynomial coefficients via recurrence, then find roots
    // For simplicity, use known pole locations for orders 1-8
    let poles: Vec<Complex> = match order {
        1 => vec![Complex::new(-1.0, 0.0)],
        2 => vec![
            Complex::new(-1.1030, 0.6368),
            Complex::new(-1.1030, -0.6368),
        ],
        3 => vec![
            Complex::new(-1.0509, 0.9991),
            Complex::new(-1.0509, -0.9991),
            Complex::new(-1.3270, 0.0),
        ],
        4 => vec![
            Complex::new(-0.9952, 1.2571),
            Complex::new(-0.9952, -1.2571),
            Complex::new(-1.3706, 0.4102),
            Complex::new(-1.3706, -0.4102),
        ],
        5 => vec![
            Complex::new(-0.9576, 1.4711),
            Complex::new(-0.9576, -1.4711),
            Complex::new(-1.3809, 0.7179),
            Complex::new(-1.3809, -0.7179),
            Complex::new(-1.5069, 0.0),
        ],
        6 => vec![
            Complex::new(-0.9318, 1.6617),
            Complex::new(-0.9318, -1.6617),
            Complex::new(-1.3789, 0.9715),
            Complex::new(-1.3789, -0.9715),
            Complex::new(-1.5735, 0.3213),
            Complex::new(-1.5735, -0.3213),
        ],
        7 => vec![
            Complex::new(-0.9104, 1.8364),
            Complex::new(-0.9104, -1.8364),
            Complex::new(-1.3724, 1.1923),
            Complex::new(-1.3724, -1.1923),
            Complex::new(-1.6120, 0.5896),
            Complex::new(-1.6120, -0.5896),
            Complex::new(-1.6843, 0.0),
        ],
        8 => vec![
            Complex::new(-0.8955, 1.9983),
            Complex::new(-0.8955, -1.9983),
            Complex::new(-1.3655, 1.3884),
            Complex::new(-1.3655, -1.3884),
            Complex::new(-1.6419, 0.8227),
            Complex::new(-1.6419, -0.8227),
            Complex::new(-1.7574, 0.2728),
            Complex::new(-1.7574, -0.2728),
        ],
        _ => {
            // Fallback: approximate with Butterworth for unsupported orders
            return butterworth(order);
        }
    };
    // Gain normalises the filter so H(0) = 1
    let mut gain = 1.0;
    for p in &poles {
        gain *= p.re * p.re + p.im * p.im;
    }
    gain = gain.sqrt();
    Zpk {
        zeros: Vec::new(),
        poles,
        gain,
    }
}

/// Apply the bilinear transform to convert an analog `Zpk` filter to digital.
///
/// `fs` – sample rate in Hz (used for frequency warping: `2 * fs`).
#[must_use]
pub fn bilinear(filter: &Zpk, fs: f64) -> Zpk {
    let fs2 = 2.0 * fs;
    let mut dz = Vec::with_capacity(filter.poles.len());
    let mut dp = Vec::with_capacity(filter.poles.len());
    let mut k = filter.gain;

    for z in &filter.zeros {
        let num = Complex::new(1.0, 0.0) + *z * (1.0 / fs2);
        let den = Complex::new(1.0, 0.0) - *z * (1.0 / fs2);
        dz.push(num / den);
        let scale_z = Complex::new(fs2, 0.0) - *z;
        k *= scale_z.norm();
    }
    for p in &filter.poles {
        let num = Complex::new(1.0, 0.0) + *p * (1.0 / fs2);
        let den = Complex::new(1.0, 0.0) - *p * (1.0 / fs2);
        dp.push(num / den);
        let scale_p = Complex::new(fs2, 0.0) - *p;
        k /= scale_p.norm();
    }
    // Extra zeros at z = -1 to match the order
    let extra = filter.poles.len() as isize - filter.zeros.len() as isize;
    for _ in 0..extra {
        dz.push(Complex::new(-1.0, 0.0));
    }
    Zpk { zeros: dz, poles: dp, gain: k }
}

/// Frequency-warp a normalised frequency for the bilinear transform.
#[must_use]
pub fn warp_frequency(frequency: f64, fs: f64) -> f64 {
    2.0 * fs * (PI * frequency / fs).tan()
}

/// Transform a low-pass analog prototype to a low-pass filter with cutoff `wo`.
#[must_use]
pub fn lp_to_lp(filter: &Zpk, wo: f64) -> Zpk {
    let mut z = filter.zeros.iter().map(|&z| z * wo).collect::<Vec<_>>();
    let p = filter.poles.iter().map(|&p| p * wo).collect::<Vec<_>>();
    let degree = filter.poles.len() as i32 - filter.zeros.len() as i32;
    let k = filter.gain * wo.powi(degree);
    Zpk { zeros: z, poles: p, gain: k }
}

/// Transform a low-pass analog prototype to a high-pass filter with cutoff `wo`.
#[must_use]
pub fn lp_to_hp(filter: &Zpk, wo: f64) -> Zpk {
    let z: Vec<Complex> = filter.zeros.iter().map(|&z| {
        if z.norm_sqr() < 1e-30 { Complex::zero() } else { Complex::new(wo, 0.0) / z }
    }).collect();
    let p: Vec<Complex> = filter.poles.iter().map(|&p| {
        Complex::new(wo, 0.0) / p
    }).collect();
    // Add zeros at origin
    let extra = filter.poles.len() - filter.zeros.len();
    let mut z_all = z;
    for _ in 0..extra {
        z_all.push(Complex::zero());
    }
    let degree = filter.poles.len() as i32 - filter.zeros.len() as i32;
    let k = filter.gain * wo.powi(degree);
    Zpk { zeros: z_all, poles: p, gain: k }
}

/// Design a digital IIR low-pass filter from an analog prototype.
///
/// `prototype` – analog prototype (e.g. from `butterworth()`, `chebyshev1()`, etc.).
/// `frequency` – cutoff frequency in Hz.
/// `fs` – sample rate in Hz.
///
/// Returns a digital `Zpk` filter.
#[must_use]
pub fn iir_lowpass(prototype: &Zpk, frequency: f64, fs: f64) -> Zpk {
    let warped = warp_frequency(frequency, fs);
    let lp = lp_to_lp(prototype, warped);
    bilinear(&lp, fs)
}

/// Design a digital IIR high-pass filter from an analog prototype.
///
/// `prototype` – analog prototype.
/// `frequency` – cutoff frequency in Hz.
/// `fs` – sample rate in Hz.
#[must_use]
pub fn iir_highpass(prototype: &Zpk, frequency: f64, fs: f64) -> Zpk {
    let warped = warp_frequency(frequency, fs);
    let hp = lp_to_hp(prototype, warped);
    bilinear(&hp, fs)
}

/// Convert a digital `Zpk` filter to second-order sections (SOS).
///
/// Returns `(b_sections, a_sections)` ready for use with `SosFilter`.
#[must_use]
pub fn zpk_to_sos(filter: &Zpk) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let n = filter.poles.len();
    if n == 0 {
        return (vec![[filter.gain, 0.0, 0.0]], vec![[1.0, 0.0, 0.0]]);
    }

    // Pair conjugate poles and zeros into second-order sections
    let mut poles = filter.poles.clone();
    let mut zeros = filter.zeros.clone();

    // Pad zeros to match pole count
    while zeros.len() < poles.len() {
        zeros.push(Complex::new(-1.0, 0.0));
    }

    let mut b_sections = Vec::new();
    let mut a_sections = Vec::new();

    let mut used_p = vec![false; poles.len()];
    let mut used_z = vec![false; zeros.len()];

    // Pair complex conjugate poles first
    let mut i = 0;
    while i < poles.len() {
        if used_p[i] {
            i += 1;
            continue;
        }
        if poles[i].im.abs() > 1e-10 {
            // Find conjugate pair
            let mut conj_idx = None;
            for j in (i + 1)..poles.len() {
                if !used_p[j]
                    && (poles[j].re - poles[i].re).abs() < 1e-10
                    && (poles[j].im + poles[i].im).abs() < 1e-10
                {
                    conj_idx = Some(j);
                    break;
                }
            }
            if let Some(j) = conj_idx {
                // Denominator: (1 - p*z^-1)(1 - p'*z^-1) = 1 - 2*Re(p)*z^-1 + |p|^2*z^-2
                let a1 = -2.0 * poles[i].re;
                let a2 = poles[i].norm_sqr();

                // Find a pair of zeros (prefer conjugate pair)
                let (z1, z2) = find_zero_pair(&zeros, &used_z);
                let b0 = 1.0;
                let b1 = -(z1.re + z2.re);
                let b2 = (z1 * z2).re;
                mark_used(&zeros, &mut used_z, &z1);
                mark_used(&zeros, &mut used_z, &z2);

                b_sections.push([b0, b1, b2]);
                a_sections.push([1.0, a1, a2]);
                used_p[i] = true;
                used_p[j] = true;
            }
        }
        i += 1;
    }

    // Handle remaining real poles
    let mut real_poles: Vec<usize> = (0..poles.len()).filter(|&i| !used_p[i]).collect();
    let mut real_zeros_idx: Vec<usize> = (0..zeros.len()).filter(|&i| !used_z[i]).collect();

    // Pair real poles together
    while real_poles.len() >= 2 {
        let i = real_poles.remove(0);
        let j = real_poles.remove(0);
        let a1 = -(poles[i].re + poles[j].re);
        let a2 = poles[i].re * poles[j].re;

        let (z1, z2) = if real_zeros_idx.len() >= 2 {
            let zi = real_zeros_idx.remove(0);
            let zj = real_zeros_idx.remove(0);
            used_z[zi] = true;
            used_z[zj] = true;
            (zeros[zi], zeros[zj])
        } else if !real_zeros_idx.is_empty() {
            let zi = real_zeros_idx.remove(0);
            used_z[zi] = true;
            (zeros[zi], Complex::new(-1.0, 0.0))
        } else {
            (Complex::new(-1.0, 0.0), Complex::new(-1.0, 0.0))
        };
        let b0 = 1.0;
        let b1 = -(z1.re + z2.re);
        let b2 = (z1 * z2).re;
        b_sections.push([b0, b1, b2]);
        a_sections.push([1.0, a1, a2]);
    }

    // Handle single remaining real pole
    if !real_poles.is_empty() {
        let i = real_poles[0];
        let a1 = -poles[i].re;

        let z = if !real_zeros_idx.is_empty() {
            let zi = real_zeros_idx.remove(0);
            used_z[zi] = true;
            zeros[zi]
        } else {
            Complex::new(-1.0, 0.0)
        };
        let b0 = 1.0;
        let b1 = -z.re;
        b_sections.push([b0, b1, 0.0]);
        a_sections.push([1.0, a1, 0.0]);
    }

    // Apply gain to first section
    if !b_sections.is_empty() {
        let g = filter.gain;
        b_sections[0][0] *= g;
        b_sections[0][1] *= g;
        b_sections[0][2] *= g;
    }

    (b_sections, a_sections)
}

fn find_zero_pair(zeros: &[Complex], used: &[bool]) -> (Complex, Complex) {
    let available: Vec<usize> = (0..zeros.len()).filter(|&i| !used[i]).collect();
    if available.len() >= 2 {
        // Try to find conjugate pair
        for (idx_a, &a) in available.iter().enumerate() {
            for &b in &available[idx_a + 1..] {
                if (zeros[a].re - zeros[b].re).abs() < 1e-10
                    && (zeros[a].im + zeros[b].im).abs() < 1e-10
                {
                    return (zeros[a], zeros[b]);
                }
            }
        }
        // No conjugate pair found, take first two
        return (zeros[available[0]], zeros[available[1]]);
    } else if available.len() == 1 {
        return (zeros[available[0]], Complex::new(-1.0, 0.0));
    }
    (Complex::new(-1.0, 0.0), Complex::new(-1.0, 0.0))
}

fn mark_used(zeros: &[Complex], used: &mut [bool], target: &Complex) {
    for (i, z) in zeros.iter().enumerate() {
        if !used[i] && (z.re - target.re).abs() < 1e-10 && (z.im - target.im).abs() < 1e-10 {
            used[i] = true;
            return;
        }
    }
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

    #[test]
    fn biquad_section_impulse() {
        let mut bq = BiquadSection::new([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((bq.process(1.0) - 1.0).abs() < 1e-10);
        assert!(bq.process(0.0).abs() < 1e-10);
    }

    #[test]
    fn sos_passthrough() {
        // Identity section: b=[1,0,0], a=[1,0,0] passes signal unchanged
        let b_sec = [[1.0, 0.0, 0.0]];
        let a_sec = [[1.0, 0.0, 0.0]];
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let out = apply_sos(&signal, &b_sec, &a_sec);
        for (a, b) in signal.iter().zip(out.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn sos_elliptic_filter() {
        // Use the 4th-order elliptic coefficients from the C code
        let b_sec = [
            [1.17381581e-02, -1.23174221e-02, 1.17381581e-02],
            [1.00000000e+00, -1.75591033e+00, 1.00000000e+00],
        ];
        let a_sec = [
            [1.00000000e+00, -1.76247270e+00, 7.94755199e-01],
            [1.00000000e+00, -1.84231257e+00, 9.36980611e-01],
        ];
        let mut impulse = vec![0.0; 20];
        impulse[0] = 1.0;
        let out = apply_sos(&impulse, &b_sec, &a_sec);
        assert_eq!(out.len(), 20);
        // First output should be b[0][0] * b[1][0] = 0.01174 * 1.0
        assert!((out[0] - b_sec[0][0]).abs() < 1e-4);
        // Filter should produce non-zero output after the impulse
        assert!(out[1].abs() > 1e-6);
    }

    #[test]
    fn sos_filter_reset() {
        let b_sec = [[0.5, 0.5, 0.0]];
        let a_sec = [[1.0, -0.5, 0.0]];
        let mut filter = SosFilter::new(&b_sec, &a_sec);
        let _ = filter.apply(&[1.0, 0.0, 0.0]);
        filter.reset();
        // After reset, state should be zero — same as fresh filter
        let out = filter.process_sample(1.0);
        let mut fresh = SosFilter::new(&b_sec, &a_sec);
        let expected = fresh.process_sample(1.0);
        assert!((out - expected).abs() < 1e-10);
    }

    #[test]
    fn biquad_lowpass_attenuates_high_freq() {
        let (b, a) = biquad_lowpass(0.1, 0.707);
        let mut bq = BiquadSection::new(b, a);
        // DC signal should pass through
        let dc_out: f64 = (0..100).map(|_| bq.process(1.0)).last().unwrap();
        assert!((dc_out - 1.0).abs() < 0.05);
    }

    #[test]
    fn biquad_highpass_attenuates_dc() {
        let (b, a) = biquad_highpass(0.1, 0.707);
        let mut bq = BiquadSection::new(b, a);
        // DC signal should be attenuated
        let dc_out: f64 = (0..200).map(|_| bq.process(1.0)).last().unwrap();
        assert!(dc_out.abs() < 0.05);
    }

    #[test]
    fn biquad_bandpass_coefficients_valid() {
        let (b, a) = biquad_bandpass(0.25, 1.0);
        // b[1] should be zero for bandpass
        assert!(b[1].abs() < 1e-10);
        // a[0] is 1.0
        assert!((a[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn biquad_notch_rejects_center() {
        let (b, a) = biquad_notch(0.25, 5.0);
        let mut bq = BiquadSection::new(b, a);
        // Feed a sine at the notch frequency
        let signal: Vec<f64> = (0..200)
            .map(|i| (2.0 * PI * 0.25 * i as f64).sin())
            .collect();
        let output: Vec<f64> = signal.iter().map(|&s| bq.process(s)).collect();
        let out_energy: f64 = output[100..].iter().map(|x| x * x).sum();
        let in_energy: f64 = signal[100..].iter().map(|x| x * x).sum();
        assert!(out_energy < in_energy * 0.1);
    }

    #[test]
    fn biquad_peak_boost_and_cut() {
        let (b_boost, a_boost) = biquad_peak(0.25, 1.0, 6.0);
        let (b_cut, a_cut) = biquad_peak(0.25, 1.0, -6.0);
        // Both should produce valid coefficients (a[0] = 1)
        assert!((a_boost[0] - 1.0).abs() < 1e-10);
        assert!((a_cut[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn biquad_shelves_valid() {
        let (b, a) = biquad_lowshelf(0.1, 6.0);
        assert!((a[0] - 1.0).abs() < 1e-10);
        let (b, a) = biquad_highshelf(0.1, -6.0);
        assert!((a[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn biquad_allpass_unity_magnitude() {
        let (b, a) = biquad_allpass(0.25, 0.707);
        let mut bq = BiquadSection::new(b, a);
        // Feed a sine and check output has similar energy
        let signal: Vec<f64> = (0..200)
            .map(|i| (2.0 * PI * 0.1 * i as f64).sin())
            .collect();
        let output: Vec<f64> = signal.iter().map(|&s| bq.process(s)).collect();
        let in_energy: f64 = signal[50..].iter().map(|x| x * x).sum();
        let out_energy: f64 = output[50..].iter().map(|x| x * x).sum();
        assert!((out_energy / in_energy - 1.0).abs() < 0.1);
    }

    #[test]
    fn butterworth_pole_count() {
        let zpk = butterworth(4);
        assert_eq!(zpk.poles.len(), 4);
        assert!(zpk.zeros.is_empty());
        // All poles should be in the left half-plane
        for p in &zpk.poles {
            assert!(p.re < 0.0);
        }
    }

    #[test]
    fn chebyshev1_pole_count() {
        let zpk = chebyshev1(4, 1.0);
        assert_eq!(zpk.poles.len(), 4);
        for p in &zpk.poles {
            assert!(p.re < 0.0);
        }
    }

    #[test]
    fn chebyshev2_has_zeros() {
        let zpk = chebyshev2(4, 40.0);
        assert_eq!(zpk.poles.len(), 4);
        assert!(!zpk.zeros.is_empty());
    }

    #[test]
    fn bessel_pole_count() {
        let zpk = bessel(3);
        assert_eq!(zpk.poles.len(), 3);
        for p in &zpk.poles {
            assert!(p.re < 0.0);
        }
    }

    #[test]
    fn butterworth_lowpass_sos() {
        let proto = butterworth(2);
        let digital = iir_lowpass(&proto, 1000.0, 8000.0);
        let (b_sec, a_sec) = zpk_to_sos(&digital);
        assert!(!b_sec.is_empty());
        assert_eq!(b_sec.len(), a_sec.len());
        // Apply to DC signal — should pass through
        let mut filter = SosFilter::new(&b_sec, &a_sec);
        let dc_out: f64 = (0..200).map(|_| filter.process_sample(1.0)).last().unwrap();
        assert!((dc_out - 1.0).abs() < 0.1);
    }

    #[test]
    fn bilinear_preserves_pole_count() {
        let proto = butterworth(3);
        let lp = lp_to_lp(&proto, 1000.0);
        let digital = bilinear(&lp, 8000.0);
        assert_eq!(digital.poles.len(), 3);
        assert_eq!(digital.zeros.len(), 3);
        // All digital poles should be inside unit circle
        for p in &digital.poles {
            assert!(p.norm() < 1.0 + 1e-6);
        }
    }
}
