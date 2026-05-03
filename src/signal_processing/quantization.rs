use std::f64::consts::PI;

/// Simulate sampling a continuous sine-wave signal at a given time.
#[must_use]
pub fn simulate_sample(time: f64, frequency: f64) -> f64 {
    (2.0 * PI * frequency * time).sin()
}

/// Quantise a normalised sample (assumed in [-1, 1]) to an N-bit integer.
///
/// Returns a value in `[0, 2^bits - 1]`.
#[must_use]
pub fn quantize_sample(sample: f64, bits: u32) -> u32 {
    if bits == 0 {
        return 0;
    }
    let levels = 1u32 << bits;
    let scaled = (sample + 1.0) / 2.0 * (levels - 1) as f64;
    scaled.round().clamp(0.0, (levels - 1) as f64) as u32
}

/// Simulate an ADC conversion: sample → normalise → quantise.
#[must_use]
pub fn simulate_adc(time: f64, frequency: f64, bits: u32) -> u32 {
    if bits == 0 {
        return 0;
    }
    let analog = (2.0 * PI * frequency * time).sin();
    let normalized = (analog + 1.0) / 2.0;
    let levels = 1u32 << bits;
    (normalized * (levels - 1) as f64).round().clamp(0.0, (levels - 1) as f64) as u32
}

/// Capture a series of ADC samples over a time interval.
#[must_use]
pub fn capture_samples(
    frequency: f64,
    bits: u32,
    start_time: f64,
    end_time: f64,
    interval: f64,
) -> Vec<u32> {
    if interval <= 0.0 || end_time <= start_time {
        return Vec::new();
    }
    let count = ((end_time - start_time) / interval) as usize;
    let mut samples = Vec::with_capacity(count);
    let mut t = start_time;
    while t < end_time {
        samples.push(simulate_adc(t, frequency, bits));
        t += interval;
    }
    samples
}

/// First-order IIR low-pass filter.
///
/// `alpha` controls the balance between the current sample and previous output.
#[must_use]
pub fn lowpass_filter_sample(current_sample: f64, previous_output: f64, alpha: f64) -> f64 {
    alpha * current_sample + (1.0 - alpha) * previous_output
}

/// Process a buffer by applying a gain factor.
#[must_use]
pub fn apply_gain(input: &[f64], gain: f64) -> Vec<f64> {
    input.iter().map(|&v| v * gain).collect()
}

/// Compute the quantisation error between the original analog value and its
/// digital approximation.
#[must_use]
pub fn quantization_error(analog_value: f64, quantized_value: u32, bits: u32) -> f64 {
    if bits == 0 {
        return analog_value;
    }
    let levels = 1u32 << bits;
    if levels <= 1 {
        return analog_value;
    }
    let normalized_quantized = quantized_value as f64 / (levels - 1) as f64;
    let reconverted = 2.0 * normalized_quantized - 1.0;
    analog_value - reconverted
}

/// Calibrate a raw ADC sample by subtracting an offset and applying a scale.
#[must_use]
pub fn calibrate_adc(raw_sample: f64, offset: f64, scale: f64) -> f64 {
    (raw_sample - offset) * scale
}

/// Convert a quantised ADC value back to the analog range [-1, 1].
#[must_use]
pub fn dac_reconstruct(quantized_value: u32, bits: u32) -> f64 {
    if bits == 0 {
        return 0.0;
    }
    let levels = 1u32 << bits;
    if levels <= 1 {
        return 0.0;
    }
    2.0 * quantized_value as f64 / (levels - 1) as f64 - 1.0
}

/// Compute the Signal-to-Quantisation-Noise Ratio (SQNR) in dB.
///
/// Theoretical SQNR ≈ 6.02 * bits + 1.76 dB for a full-scale sinusoid.
#[must_use]
pub fn theoretical_sqnr_db(bits: u32) -> f64 {
    6.02 * bits as f64 + 1.76
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantize_midpoint() {
        // Sample = 0.0 → midpoint of [0, 255] for 8 bits
        let q = quantize_sample(0.0, 8);
        assert_eq!(q, 128); // (0+1)/2 * 255 = 127.5 → rounds to 128
    }

    #[test]
    fn quantize_extremes() {
        assert_eq!(quantize_sample(-1.0, 8), 0);
        assert_eq!(quantize_sample(1.0, 8), 255);
    }

    #[test]
    fn quantization_error_small_for_8bit() {
        let sample = 0.5;
        let q = quantize_sample(sample, 8);
        let error = quantization_error(sample, q, 8);
        assert!(error.abs() < 0.01);
    }

    #[test]
    fn dac_reconstruct_roundtrip() {
        let original = 0.5;
        let q = quantize_sample(original, 16);
        let reconstructed = dac_reconstruct(q, 16);
        assert!((original - reconstructed).abs() < 0.001);
    }

    #[test]
    fn lowpass_filter_smooths() {
        let mut output = 0.0;
        for _ in 0..100 {
            output = lowpass_filter_sample(1.0, output, 0.1);
        }
        assert!((output - 1.0).abs() < 0.01);
    }

    #[test]
    fn theoretical_sqnr_8bit() {
        let sqnr = theoretical_sqnr_db(8);
        assert!((sqnr - 49.92).abs() < 0.1);
    }
}
