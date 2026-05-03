#include <iostream>
#include <cmath>
#include <vector>

using std::vector;
using std::cout;
using std::endl;
using std::sin;
using std::round;

// Define M_PI if not already defined.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Function to simulate sampling a continuous sine-wave signal at a given time instance.
double simulateSample(double time, double frequency) {
    // Computes the instantaneous amplitude of a sine wave.
    return std::sin(2 * M_PI * frequency * time);
}

// Quantize a normalized sample (assumed to be in [-1, 1]) into an integer based on the provided bit depth.
int quantizeSample(double sample, int bits) {
    if (bits <= 0) {
        return 0;
    }
    int levels = 1 << bits; // Total quantization levels: 2^bits.
    // Map sample from [-1, 1] to [0, levels-1] and round to nearest integer.
    double scaledSample = (sample + 1.0) / 2.0 * (levels - 1);
    return static_cast<int>(std::round(scaledSample));
}

// Simulate an Analog-to-Digital Conversion (ADC) by sampling, normalizing, and quantizing the signal.
int simulateADCConversion(double time, double frequency, int bits) {
    if (bits <= 0) {
        return 0;
    }
    double analogValue = std::sin(2 * M_PI * frequency * time);
    // Normalize the value from [-1, 1] to [0, 1].
    double normalized = (analogValue + 1.0) / 2.0;
    int levels = 1 << bits;
    return static_cast<int>(std::round(normalized * (levels - 1)));
}

// Capture a series of ADC samples over a specified time interval using a fixed sampling interval.
std::vector<int> captureSamples(double frequency, int bits, double startTime, double endTime, double interval) {
    if (interval <= 0) {
        return {};
    }
    std::vector<int> samples;
    int count = static_cast<int>((endTime - startTime) / interval);
    if (count > 0) {
        samples.reserve(count);
        for (double t = startTime; t < endTime; t += interval) {
            samples.push_back(simulateADCConversion(t, frequency, bits));
        }
    }
    return samples;
}

// Apply a simple first-order IIR low-pass filter to smooth the sampled signal.
double applyLowPassFilter(double currentSample, double previousOutput, double alpha) {
    // 'alpha' determines the weighting of the current sample vs. the previous output.
    return alpha * currentSample + (1.0 - alpha) * previousOutput;
}

// Process an array of samples by applying a gain factor using pointer arithmetic for optimization.
void processSampleBuffer(const double* input, double* output, int size, double gain) {
    const double* end = input + size;
    while (input < end) {
        *output = (*input) * gain;
        ++input;
        ++output;
    }
}

// Compute the quantization error between the original analog value and its digital approximation.
double computeQuantizationError(double analogValue, int quantizedValue, int bits) {
    if (bits <= 0) {
        return analogValue;
    }
    int levels = 1 << bits;
    if (levels <= 1) {
        return analogValue;
    }
    // Convert the quantized value back to a normalized value in [0, 1].
    double normalizedQuantized = static_cast<double>(quantizedValue) / (levels - 1);
    // Reconversion to the original analog range [-1, 1].
    double reconvertedValue = 2.0 * normalizedQuantized - 1.0;
    return analogValue - reconvertedValue;
}

// Calibrate a raw ADC sample by subtracting an offset and applying a scaling factor.
double calibrateADC(double rawSample, double offset, double scale) {
    return (rawSample - offset) * scale;
}

// Main function demonstrating the key digital signal processing algorithms.
int main() {
    // Parameters for signal simulation and conversion.
    double frequency = 5.0; // Signal frequency in Hz.
    double startTime = 0.0; // Start time for sample capture (seconds).
    double endTime = 1.0; // End time for sample capture (seconds).
    double interval = 0.01; // Sampling interval (seconds).
    int bits = 8; // Bit depth for ADC conversion.

    // Simulate a single sample capture.
    double currentTime = 0.05;
    double sampleValue = simulateSample(currentTime, frequency);
    std::cout << "Simulated continuous sample at time " << currentTime << "s: " << sampleValue << std::endl;

    // Quantize the sample and output the quantized value.
    int quantizedValue = quantizeSample(sampleValue, bits);
    std::cout << "Quantized sample (using " << bits << " bits): " << quantizedValue << std::endl;

    // Perform a full ADC conversion simulation.
    int adcValue = simulateADCConversion(currentTime, frequency, bits);
    std::cout << "Simulated ADC conversion value: " << adcValue << std::endl;

    // Capture multiple samples within the specified time window.
    std::vector<int> adcSamples = captureSamples(frequency, bits, startTime, endTime, interval);
    std::cout << "Captured " << adcSamples.size() << " ADC samples." << std::endl;

    // Demonstrate low-pass filtering over captured samples.
    double filteredValue = 0.0;
    const double LOW_PASS_ALPHA = 0.1;
double alpha = LOW_PASS_ALPHA; // Filter coefficient to control smoothing.
    for (size_t i = 0; i < adcSamples.size(); ++i) {
        // Convert quantized ADC value back to [-1, 1] range.
        double currentSample = 2.0 * adcSamples[i] / ((1 << bits) - 1) - 1.0;
        filteredValue = applyLowPassFilter(currentSample, filteredValue, alpha);
    }
    std::cout << "Filtered output after applying low-pass filter: " << filteredValue << std::endl;

    // Example of processing a small fixed buffer using pointer arithmetic.
    const int bufferSize = 5;
    double inputBuffer[bufferSize] = {0.1, 0.2, -0.1, 0.3, -0.2};
    double outputBuffer[bufferSize] = {0};
    double gain = 2.0;
    processSampleBuffer(inputBuffer, outputBuffer, bufferSize, gain);
    std::cout << "Processed sample buffer with gain " << gain << ":" << std::endl;
    for (int i = 0; i < bufferSize; ++i) {
        std::cout << outputBuffer[i] << " ";
    }
    std::cout << std::endl;

    // Analyze quantization error for the previously sampled and quantized value.
    double quantError = computeQuantizationError(sampleValue, quantizedValue, bits);
    std::cout << "Quantization error: " << quantError << std::endl;

    // Demonstrate ADC calibration using example offset and scaling parameters.
    double rawADC = static_cast<double>(adcValue);
    double offset = 0.05; // Example offset correction.
    double scale = 1.02; // Example scaling factor.
    double calibratedValue = calibrateADC(rawADC, offset, scale);
    std::cout << "Calibrated ADC value: " << calibratedValue << std::endl;

    return 0;
}
