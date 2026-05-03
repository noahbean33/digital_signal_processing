#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <numeric>
#include <algorithm>
#include <random>

// Constant for PI
const double PI = std::acos(-1);

// Function to calculate FIR filter coefficients using the window method.
// It computes the ideal impulse response using the sinc function and applies a Hamming window.
std::vector<double> calculateFIRCoefficients(int taps, double cutoff) {
    if (taps <= 1) {
        return {};
    }
    std::vector<double> coeffs(taps);
    int M = taps - 1;
    double sum = 0.0;
    for (int n = 0; n < taps; ++n) {
        double pos = n - M / 2.0;
        double ideal = (std::abs(pos) < 1e-6)
                         ? 1.0
                         : std::sin(PI * cutoff * pos) / (PI * cutoff * pos);
        double window = 0.54 - 0.46 * std::cos(2 * PI * n / M);
        coeffs[n] = ideal * window;
        sum += coeffs[n];
    }
    // Normalizing the coefficients to ensure unity gain
    if (std::abs(sum) > 1e-9) {
        for (int n = 0; n < taps; ++n) {
            coeffs[n] /= sum;
        }
    }
    return coeffs;
}

// Function to apply FIR filter to an input signal using convolution
std::vector<double> applyFIRFilter(const std::vector<double>& signal, const std::vector<double>& coeffs) {
    if (signal.size() < coeffs.size() || coeffs.empty()) {
        return {};
    }
    int outputSize = signal.size() - coeffs.size() + 1;
    std::vector<double> output(outputSize, 0.0);
    for (std::size_t n = 0; n < output.size(); ++n) {
        for (std::size_t k = 0; k < coeffs.size(); ++k) {
            output[n] += signal[n + k] * coeffs[k];
        }
    }
    return output;
}

// Optimized FIR filtering using the standard library inner_product
std::vector<double> optimizedFIRFilter(const std::vector<double>& signal, const std::vector<double>& coeffs) {
    if (signal.size() < coeffs.size() || coeffs.empty()) {
        return {};
    }
    int outSize = signal.size() - coeffs.size() + 1;
    std::vector<double> output(outSize, 0.0);
    for (int i = 0; i < outSize; ++i) {
        output[i] = std::inner_product(coeffs.begin(), coeffs.end(), signal.begin() + i, 0.0);
    }
    return output;
}

// Function to compute the complex frequency response of the FIR filter at a given frequency.
std::complex<double> computeFIRFrequencyResponse(
    const std::vector<double>& coeffs,
    double frequency,
    double samplingRate) {
    std::complex<double> response(0.0, 0.0);
    for (std::size_t n = 0; n < coeffs.size(); ++n) {
        response += coeffs[n] * std::exp(std::complex<double>(0.0, -2.0 * PI * frequency * n / samplingRate));
    }
    return response;
}

// Function to quantize FIR filter coefficients using a fixed-point scaling factor (qFactor).
std::vector<int> quantizeFIRCoefficients(const std::vector<double>& coeffs, int qFactor) {
    std::vector<int> qCoeffs;
    qCoeffs.reserve(coeffs.size());
    for (double c : coeffs) {
        qCoeffs.push_back(static_cast<int>(std::round(c * qFactor)));
    }
    return qCoeffs;
}

// Function to compute Mean Squared Error (MSE) between two frequency response magnitudes.
double computeMSE(const std::vector<double>& theoretical, const std::vector<double>& computed) {
    if (theoretical.size() != computed.size() || theoretical.empty()) {
        return 0.0;
    }
    double mse = 0.0;
    size_t n = theoretical.size();
    for (size_t i = 0; i < n; ++i) {
        double diff = theoretical[i] - computed[i];
        mse += diff * diff;
    }
    return mse / n;
}

// Monte Carlo simulation to evaluate filter robustness under random perturbations.
void monteCarloSimulation(int trials, int taps, double cutoff, double samplingRate) {
    std::default_random_engine rng;
    std::normal_distribution<double> noise(0.0, 0.01); // Noise with standard deviation 0.01
    std::vector<double> mseResults;

    // Define frequency bins for analysis (from 0 to Nyquist frequency)
    int numFreqBins = 128;
    std::vector<double> freqs(numFreqBins);
    for (int i = 0; i < numFreqBins; ++i) {
        freqs[i] = i * samplingRate / (2.0 * numFreqBins);
    }

    // Compute ideal frequency response for the filter
    std::vector<double> coeffsIdeal = calculateFIRCoefficients(taps, cutoff);
    std::vector<double> idealResponse(numFreqBins, 0.0);
    for (int i = 0; i < numFreqBins; ++i) {
        std::complex<double> resp = computeFIRFrequencyResponse(coeffsIdeal, freqs[i], samplingRate);
        idealResponse[i] = std::abs(resp);
    }

    // Run simulation trials
    for (int t = 0; t < trials; ++t) {
        std::vector<double> coeffsPerturbed = coeffsIdeal;
        // Introduce random perturbations
        for (auto& c : coeffsPerturbed) {
            c += noise(rng);
        }

        // Compute frequency response for perturbed coefficients
        std::vector<double> perturbedResponse(numFreqBins, 0.0);
        for (int i = 0; i < numFreqBins; ++i) {
            std::complex<double> resp = computeFIRFrequencyResponse(coeffsPerturbed, freqs[i], samplingRate);
            perturbedResponse[i] = std::abs(resp);
        }
        double mse = computeMSE(idealResponse, perturbedResponse);
        mseResults.push_back(mse);
    }

    // Calculate average MSE across all trials
    if (!mseResults.empty()) {
        double averageMSE = std::accumulate(mseResults.begin(), mseResults.end(), 0.0) / mseResults.size();
        std::cout << "Monte Carlo Simulation Average MSE: " << averageMSE << std::endl;
    }
}

int main() {
    // Filter design parameters
    int taps = 51;
    double cutoff = 0.25; // Normalized cutoff frequency (0 to 0.5)
    double samplingRate = 1000.0; // Sampling rate in Hz
    int qFactor = 10000; // Quantization scaling factor

    // Generate a test signal (e.g., a sine wave)
    int signalLength = 500;
    std::vector<double> signal(signalLength);
    double freqSignal = 50.0; // Frequency of the sine wave in Hz
    for (int i = 0; i < signalLength; ++i) {
        signal[i] = std::sin(2 * PI * freqSignal * i / samplingRate);
    }

    // Calculate FIR filter coefficients using Hamming window method
    std::vector<double> coeffs = calculateFIRCoefficients(taps, cutoff);

    // Apply FIR filter using standard convolution
    std::vector<double> filteredSignal = applyFIRFilter(signal, coeffs);

    // Apply FIR filter using optimized convolution with inner_product
    std::vector<double> optimizedFilteredSignal = optimizedFIRFilter(signal, coeffs);

    // Compute and display frequency response at a test frequency (e.g., 100 Hz)
    double testFrequency = 100.0;
    std::complex<double> freqResp = computeFIRFrequencyResponse(coeffs, testFrequency, samplingRate);
    std::cout << "Frequency response at " << testFrequency << " Hz: " << freqResp << std::endl;

    // Quantize coefficients for fixed-point arithmetic
    std::vector<int> quantizedCoeffs = quantizeFIRCoefficients(coeffs, qFactor);
    std::cout << "Quantized coefficients:";
    for (int c : quantizedCoeffs) {
        std::cout << c << " ";
    }
    std::cout << std::endl;

    // Run Monte Carlo simulation to test filter robustness
    int trials = 100;
    monteCarloSimulation(trials, taps, cutoff, samplingRate);

    return 0;
}
