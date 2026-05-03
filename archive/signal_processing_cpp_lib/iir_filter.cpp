#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <numeric>

// Process a single sample using Direct Form I structure.
double processIIRSample(const std::vector<double>& b,
                        const std::vector<double>& a,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        int n) {
    if (a.empty() || std::abs(a[0]) < 1e-9) {
        return 0.0; // Avoid division by zero
    }
    double acc = 0.0;
    int order = b.size();
    // Feedforward portion: accumulate input contributions.
    for (int i = 0; i < order; ++i) {
        if (n - i >= 0)
            acc += b[i] * x[n - i];
    }
    // Feedback portion: subtract past output contributions.
    for (int i = 1; i < a.size(); ++i) {
        if (n - i >= 0)
            acc -= a[i] * y[n - i];
    }
    return acc / a[0];
}

// Check filter stability by ensuring all poles lie inside the unit circle.
bool isFilterStable(const std::vector<std::complex<double>>& poles) {
    for (const auto& pole : poles) {
        if (std::abs(pole) >= 1.0)
            return false;
    }
    return true;
}

// Apply IIR filtering using a Direct Form II transposed structure.
std::vector<double> applyIIRDirectFormII(const std::vector<double>& input,
                                       const std::vector<double>& a,
                                       const std::vector<double>& b) {
    if (input.empty() || a.empty() || b.empty() || a.size() != b.size()) {
        return {};
    }
    std::vector<double> output(input.size(), 0.0);
    // Filter order determined by denominator coefficients (a[0] is assumed normalized).
    int order = a.size() - 1;
    // State vector (delay elements) initialized to zero.
    std::vector<double> w(order, 0.0);
    for (size_t n = 0; n < input.size(); ++n) {
        // Compute intermediate state value.
        double w0 = input[n];
        for (int i = 1; i <= order; ++i) {
            w0 -= a[i] * w[i - 1];
        }
        // Compute output sample using feedforward coefficients.
        double y = b[0] * w0;
        for (int i = 1; i <= order; ++i) {
            y += b[i] * w[i - 1];
        }
        output[n] = y;
        // Update state vector (shift the delay line).
        if (order > 0) {
            for (int i = order - 1; i > 0; --i) {
                w[i] = w[i - 1];
            }
            w[0] = w0;
        }
    }
    return output;
}

// Compute the frequency response H(e^(j)) of the IIR filter at a given frequency.
std::complex<double> computeIIRFrequencyResponse(const std::vector<double>& b,
                                                 const std::vector<double>& a,
                                                 double frequency,
                                                 double samplingRate) {
    std::complex<double> numerator(0.0, 0.0);
    std::complex<double> denominator(0.0, 0.0);
    const double PI = std::acos(-1);
    for (std::size_t k = 0; k < b.size(); ++k) {
        numerator += b[k] * std::exp(std::complex<double>(0.0, -2.0 * PI * frequency * k / samplingRate));
    }
    for (std::size_t k = 0; k < a.size(); ++k) {
        denominator += a[k] * std::exp(std::complex<double>(0.0, -2.0 * PI * frequency * k / samplingRate));
    }
    if (std::abs(denominator) < 1e-9) {
        return std::complex<double>(0.0, 0.0);
    }
    return numerator / denominator;
}

// Quantize IIR coefficients by scaling and rounding.
std::vector<int> quantizeIIRCoefficients(const std::vector<double>& coeffs, int qFactor) {
    std::vector<int> qCoeffs;
    qCoeffs.reserve(coeffs.size());
    for (double c : coeffs) {
        qCoeffs.push_back(static_cast<int>(std::round(c * qFactor)));
    }
    return qCoeffs;
}

// Normalize the numerator and denominator coefficients so that a[0] is 1.
std::pair<std::vector<double>, std::vector<double>> normalizeIIRCoefficients(
    const std::vector<double>& b,
    const std::vector<double>& a) {
    if (a.empty() || std::abs(a[0]) < 1e-9) {
        return {b, a};
    }
    std::vector<double> b_norm(b.size());
    std::vector<double> a_norm(a.size());
    double scale = a[0];
    for (std::size_t i = 0; i < b.size(); ++i) {
        b_norm[i] = b[i] / scale;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        a_norm[i] = a[i] / scale;
    }
    return {b_norm, a_norm};
}

// Update the state buffer used in real-time recursive filtering.
void updateIIRStateBuffer(std::vector<double>& stateBuffer, double newValue) {
    for (std::size_t i = stateBuffer.size() - 1; i > 0; --i) {
        stateBuffer[i] = stateBuffer[i - 1];
    }
    stateBuffer[0] = newValue;
}

// Demonstration of the IIR filter algorithms.
int main() {
    // Example IIR filter coefficients (feedforward and feedback).
    std::vector<double> b = {0.2929, 0.5858, 0.2929}; // Numerator coefficients.
    std::vector<double> a = {1.0000, -0.0000, 0.1716}; // Denominator coefficients.

    // Normalize coefficients so that a[0] equals 1.
    auto normalizedCoeffs = normalizeIIRCoefficients(b, a);
    b = normalizedCoeffs.first;
    a = normalizedCoeffs.second;

    // Create a test input signal (e.g., an impulse response).
    std::vector<double> input(20, 0.0);
    input[0] = 1.0;

    // Apply IIR filter using Direct Form II Transposed implementation.
    std::vector<double> output = applyIIRDirectFormII(input, a, b);

    // Compute frequency response at 1000 Hz for a sampling rate of 8000 Hz.
    double frequency = 1000.0;
    double samplingRate = 8000.0;
    std::complex<double> freqResponse = computeIIRFrequencyResponse(b, a, frequency, samplingRate);

    // Display filtered output.
    std::cout << "Filtered Output:" << std::endl;
    for (double val : output) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    // Display frequency response at specified frequency.
    std::cout << "Frequency Response at " << frequency << " Hz: " << freqResponse << std::endl;

    // Quantize coefficients with a quantization factor (e.g., 1000).
    int quantFactor = 1000;
    std::vector<int> qB = quantizeIIRCoefficients(b, quantFactor);
    std::vector<int> qA = quantizeIIRCoefficients(a, quantFactor);
    std::cout << "Quantized Coefficients:" << std::endl;
    std::cout << "b coefficients: ";
    for (int coeff : qB) {
        std::cout << coeff << " ";
    }
    std::cout << std::endl;
    std::cout << "a coefficients: ";
    for (int coeff : qA) {
        std::cout << coeff << " ";
    }
    std::cout << std::endl;

    // Demonstrate updating the state buffer for real-time filtering.
    std::vector<double> stateBuffer(2, 0.0); // State buffer size equal to filter order.
    double newInput = 0.5;
    updateIIRStateBuffer(stateBuffer, newInput);
    std::cout << "Updated State Buffer: ";
    for (double state : stateBuffer) {
        std::cout << state << " ";
    }
    std::cout << std::endl;

    // Example stability check with sample poles.
    std::vector<std::complex<double>> poles = {{0.8, 0.1}, {0.7, -0.3}};
    if (isFilterStable(poles)) {
        std::cout << "Filter is stable." << std::endl;
    } else {
        std::cout << "Filter is unstable." << std::endl;
    }

    return 0;
}
