#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <functional>
#include <future>
#include <algorithm>

// Compute the twiddle factor: W_N^k = exp(-2ik/N)
std::complex<double> computeTwiddleFactor(int k, int N) {
    if (N == 0) {
        return std::complex<double>(1.0, 0.0);
    }
    const double PI = std::acos(-1);
    return std::exp(std::complex<double>(0.0, -2.0 * PI * k / N));
}

// Butterfly operation: combining two complex numbers with a twiddle factor
std::pair<std::complex<double>, std::complex<double>> fftButterfly(
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> twiddle)
{
    std::complex<double> temp = b * twiddle;
    return {a + temp, a - temp};
}

// Iterative Cooley-Tukey FFT implementation (radix-2 FFT)
void iterativeFFT(std::vector<std::complex<double>>& data) {
    const size_t N = data.size();
    if (N <= 1) return;

    // Bit-reversal permutation
    bitReversalPermutation(data);

    // Iterative FFT
    for (size_t len = 2; len <= N; len <<= 1) {
        double angle = -2.0 * PI / len;
        std::complex<double> wlen(cos(angle), sin(angle));
        for (size_t i = 0; i < N; i += len) {
            std::complex<double> w(1);
            for (size_t j = 0; j < len / 2; j++) {
                std::complex<double> u = data[i + j];
                std::complex<double> v = data[i + j + len / 2] * w;
                data[i + j] = u + v;
                data[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

// Bit-reversal permutation: reorder data in-place for in-place FFT efficiency
void bitReversalPermutation(std::vector<std::complex<double>>& data) {
    const std::size_t N = data.size();
    std::size_t nBits = 0;
    while ((1u << nBits) < N)
        ++nBits;

    for (std::size_t i = 0; i < N; ++i) {
        std::size_t j = 0;
        for (std::size_t k = 0; k < nBits; ++k) {
            if (i & (1u << k))
                j |= 1u << (nBits - 1 - k);
        }
        if (j > i)
            std::swap(data[i], data[j]);
    }
}

// Apply a Hanning window to the time-domain signal to reduce spectral leakage
std::vector<double> applyHanningWindow(const std::vector<double>& signal) {
    const std::size_t N = signal.size();
    if (N <= 1) {
        return signal;
    }
    std::vector<double> windowed(N);
    const double PI = std::acos(-1);
    for (std::size_t n = 0; n < N; ++n) {
        double w = 0.5 * (1 - std::cos((2 * PI * n) / (N - 1)));
        windowed[n] = signal[n] * w;
    }
    return windowed;
}

// Measure execution time of an FFT function using high-resolution clock
std::chrono::duration<double> measureFFTExecutionTime(
    std::vector<std::complex<double>>& data,
    std::function<void(std::vector<std::complex<double>>&)> fftFunction)
{
    auto start = std::chrono::high_resolution_clock::now();
    fftFunction(data);
    auto end = std::chrono::high_resolution_clock::now();
    return end - start;
}

int main() {
    // Define a simple test signal: a combination of sine waves
    const std::size_t N = 16; // Signal length (must be a power of 2 for radix-2 FFT)
    const double PI = std::acos(-1);
    std::vector<double> signal(N);
    double frequency1 = 2.0; // Frequency for first sine wave (Hz)
    double frequency2 = 3.0; // Frequency for second sine wave (Hz)
    double samplingRate = 16.0; // Sampling rate (Hz)

    // Generate the composite signal with two sine components
    for (std::size_t n = 0; n < N; ++n) {
        signal[n] = std::sin(2 * PI * frequency1 * n / samplingRate)
                    + 0.5 * std::sin(2 * PI * frequency2 * n / samplingRate);
    }

    // Preprocess: Apply Hanning window for leakage mitigation
    auto windowedSignal = applyHanningWindow(signal);

    // Convert the real-valued signal to complex format (imaginary parts are zero)
    std::vector<std::complex<double>> complexSignal;
    for (const auto& value : windowedSignal) {
        complexSignal.push_back(std::complex<double>(value, 0.0));
    }

    // Compute FFT using the iterative FFT function
    std::vector<std::complex<double>> fftResult = complexSignal;
    iterativeFFT(fftResult);

    // Measure execution time of the FFT implementation
    auto duration = measureFFTExecutionTime(complexSignal, iterativeFFT);

    // Output the FFT results
    std::cout << "FFT Result (Iterative FFT):" << std::endl;
    for (const auto& value : fftResult) {
        std::cout << value << std::endl;
    }

    // Output measured execution time
    std::cout << "\nExecution Time (Iterative FFT): "
              << duration.count() << " seconds" << std::endl;

    // Demonstrate bit reversal permutation on the FFT result
    bitReversalPermutation(fftResult);
    std::cout << "\nFFT Result after Bit-Reversal Permutation:" << std::endl;
    for (const auto& value : fftResult) {
        std::cout << value << std::endl;
    }

    return 0;
}
