use std::f64::consts::PI;

use super::fft::Complex;

/// Map an integer symbol to its PSK constellation point.
#[must_use]
pub fn map_psk_symbol(symbol: usize, modulation_order: usize) -> Complex {
    let phase = 2.0 * PI * symbol as f64 / modulation_order as f64;
    Complex::from_polar(1.0, phase)
}

/// Map an integer symbol to its QAM constellation point.
///
/// `sqrt_mod_order` is the square root of the modulation order (e.g. 4 for 16-QAM).
#[must_use]
pub fn map_qam_symbol(symbol: usize, sqrt_mod_order: usize) -> Complex {
    let mod_order = sqrt_mod_order * sqrt_mod_order;
    if symbol >= mod_order || sqrt_mod_order == 0 {
        return Complex::zero();
    }
    let row = symbol / sqrt_mod_order;
    let col = symbol % sqrt_mod_order;
    let norm_factor = ((2.0 * (mod_order as f64 - 1.0)) / 3.0).sqrt();
    if norm_factor.abs() < 1e-9 {
        return Complex::zero();
    }
    let i_val = (2 * col as isize + 1 - sqrt_mod_order as isize) as f64 / norm_factor;
    let q_val = (2 * row as isize + 1 - sqrt_mod_order as isize) as f64 / norm_factor;
    Complex::new(i_val, q_val)
}

/// Generate an FSK signal for a given symbol.
#[must_use]
pub fn generate_fsk_signal(
    symbol: usize,
    symbol_duration: f64,
    sample_rate: f64,
    base_frequency: f64,
    frequency_deviation: f64,
) -> Vec<f64> {
    if sample_rate <= 0.0 {
        return Vec::new();
    }
    let num_samples = (symbol_duration * sample_rate) as usize;
    if num_samples == 0 {
        return Vec::new();
    }
    let frequency = base_frequency + symbol as f64 * frequency_deviation;
    (0..num_samples)
        .map(|n| {
            let time = n as f64 / sample_rate;
            (2.0 * PI * frequency * time).sin()
        })
        .collect()
}

/// Demodulate a PSK symbol from a received complex sample.
#[must_use]
pub fn demodulate_psk(received: Complex, modulation_order: usize) -> usize {
    if modulation_order == 0 {
        return 0;
    }
    let mut phase = received.arg();
    if phase < 0.0 {
        phase += 2.0 * PI;
    }
    let symbol_interval = 2.0 * PI / modulation_order as f64;
    let symbol = (phase / symbol_interval).round() as usize % modulation_order;
    symbol
}

/// Demodulate a QAM symbol from a received complex sample.
#[must_use]
pub fn demodulate_qam(received: Complex, sqrt_mod_order: usize) -> usize {
    if sqrt_mod_order == 0 {
        return 0;
    }
    let mod_order = sqrt_mod_order * sqrt_mod_order;
    let norm_factor = ((2.0 * (mod_order as f64 - 1.0)) / 3.0).sqrt();
    if norm_factor.abs() < 1e-9 {
        return 0;
    }
    let col = ((received.re * norm_factor + sqrt_mod_order as f64 - 1.0) / 2.0)
        .round()
        .clamp(0.0, (sqrt_mod_order - 1) as f64) as usize;
    let row = ((received.im * norm_factor + sqrt_mod_order as f64 - 1.0) / 2.0)
        .round()
        .clamp(0.0, (sqrt_mod_order - 1) as f64) as usize;
    row * sqrt_mod_order + col
}

/// Demodulate an FSK symbol based on zero-crossing rate.
#[must_use]
pub fn demodulate_fsk(
    segment: &[f64],
    sample_rate: f64,
    base_frequency: f64,
    frequency_deviation: f64,
    modulation_order: usize,
) -> usize {
    if segment.is_empty() || frequency_deviation.abs() < 1e-9 {
        return 0;
    }
    let num_samples = segment.len();
    let mut zero_crossings = 0usize;
    for n in 1..num_samples {
        if segment[n - 1] * segment[n] < 0.0 {
            zero_crossings += 1;
        }
    }
    let estimated_freq = (zero_crossings as f64 * sample_rate) / (2.0 * num_samples as f64);
    let symbol = ((estimated_freq - base_frequency) / frequency_deviation)
        .round()
        .clamp(0.0, (modulation_order - 1) as f64) as usize;
    symbol
}

/// Correct carrier phase offset in a received symbol.
#[must_use]
pub fn correct_carrier_phase(received: Complex, estimated_phase_offset: f64) -> Complex {
    let correction = Complex::from_polar(1.0, -estimated_phase_offset);
    received * correction
}

/// Find the nearest constellation symbol by Euclidean distance.
///
/// Returns the index of the nearest symbol in the constellation.
#[must_use]
pub fn find_nearest_symbol(received: Complex, constellation: &[Complex]) -> usize {
    let mut min_distance = f64::MAX;
    let mut best_index = 0;
    for (i, &sym) in constellation.iter().enumerate() {
        let diff = received - sym;
        let distance = diff.norm_sqr();
        if distance < min_distance {
            min_distance = distance;
            best_index = i;
        }
    }
    best_index
}

/// Compute symbol likelihoods based on a Gaussian noise model.
#[must_use]
pub fn symbol_likelihoods(
    received: Complex,
    constellation: &[Complex],
    sigma: f64,
) -> Vec<f64> {
    let sigma_sq = sigma * sigma;
    if sigma_sq < 1e-30 {
        return vec![0.0; constellation.len()];
    }
    constellation
        .iter()
        .map(|&sym| {
            let diff = received - sym;
            (-diff.norm_sqr() / (2.0 * sigma_sq)).exp()
        })
        .collect()
}

/// Equalise a received signal using zero-forcing with a channel estimate.
#[must_use]
pub fn equalize_signal(received: &[Complex], channel_estimate: &[Complex]) -> Vec<Complex> {
    if received.len() != channel_estimate.len() {
        return Vec::new();
    }
    received
        .iter()
        .zip(channel_estimate.iter())
        .map(|(&r, &h)| {
            if h.norm() > 1e-12 {
                r / h
            } else {
                r
            }
        })
        .collect()
}

/// Generate a full PSK constellation.
#[must_use]
pub fn psk_constellation(modulation_order: usize) -> Vec<Complex> {
    (0..modulation_order)
        .map(|i| map_psk_symbol(i, modulation_order))
        .collect()
}

/// Generate a full QAM constellation.
#[must_use]
pub fn qam_constellation(sqrt_mod_order: usize) -> Vec<Complex> {
    let mod_order = sqrt_mod_order * sqrt_mod_order;
    (0..mod_order)
        .map(|i| map_qam_symbol(i, sqrt_mod_order))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn psk_roundtrip() {
        for symbol in 0..8 {
            let mapped = map_psk_symbol(symbol, 8);
            let demod = demodulate_psk(mapped, 8);
            assert_eq!(demod, symbol, "PSK roundtrip failed for symbol {symbol}");
        }
    }

    #[test]
    fn qam_roundtrip() {
        for symbol in 0..16 {
            let mapped = map_qam_symbol(symbol, 4);
            let demod = demodulate_qam(mapped, 4);
            assert_eq!(demod, symbol, "QAM roundtrip failed for symbol {symbol}");
        }
    }

    #[test]
    fn fsk_generates_correct_length() {
        let sig = generate_fsk_signal(2, 0.01, 10000.0, 1000.0, 100.0);
        assert_eq!(sig.len(), 100);
    }

    #[test]
    fn phase_correction_restores_symbol() {
        let constellation = psk_constellation(8);
        let offset = 0.2;
        let received = constellation[3] * Complex::from_polar(1.0, offset);
        let corrected = correct_carrier_phase(received, offset);
        let idx = find_nearest_symbol(corrected, &constellation);
        assert_eq!(idx, 3);
    }
}
