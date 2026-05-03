#include <iostream>
#include <std::vector>
#include <cmath>
#include <algorithm>
#include <numeric>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// Discrete convolution: y[n] = sum_{m=0}^{n} x[m]*h[n-m]
std::vector<double> convolve(const std::vector<double>& signal, const std::vector<double>& kernel) {
    if (signal.empty() || kernel.empty()) {
        return {};
    }
    int n = signal.size();
    int k = kernel.size();
    int convSize = n + k - 1;
    std::vector<double> result(convSize, 0.0);
    for (int i = 0; i < convSize; ++i) {
        for (int j = 0; j < k; ++j) {
            if (i - j >= 0 && i - j < n)
                result[i] += signal[i - j] * kernel[j];
        }
    }
    return result;
}

// Basic correlation without reversing the kernel
std::vector<double> correlate(const std::vector<double>& signal, const std::vector<double>& kernel) {
    if (signal.empty() || kernel.empty()) {
        return {};
    }
    int n = signal.size();
    int k = kernel.size();
    int corrSize = n - k + 1;
    std::vector<double> result(corrSize, 0.0);
    for (int i = 0; i < corrSize; ++i) {
        double sum = 0.0;
        for (int j = 0; j < k; ++j) {
            sum += signal[i + j] * kernel[j];
        }
        result[i] = sum;
    }
    return result;
}

// Full cross-correlation between two signals
std::vector<double> crossCorrelate(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.empty() || y.empty()) {
        return {};
    }
    int n = x.size();
    int m = y.size();
    int corrSize = n + m - 1;
    std::vector<double> result(corrSize, 0.0);
    for (int lag = 0; lag < corrSize; ++lag) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            int j = lag - i;
            if (j >= 0 && j < m)
                sum += x[i] * y[j];
        }
        result[lag] = sum;
    }
    return result;
}

// Overlap-add method for efficient convolution of long signals
std::vector<double> overlapAddConvolution(const std::vector<double>& signal, const std::vector<double>& kernel, int blockSize) {
    if (signal.empty() || kernel.empty()) {
        return {};
    }
    int signalSize = signal.size();
    int kernelSize = kernel.size();
    int outputSize = signalSize + kernelSize - 1;
    std::vector<double> output(outputSize, 0.0);
    for (int start = 0; start < signalSize; start += blockSize) {
        int currentBlockSize = std::min(blockSize, signalSize - start);
        std::vector<double> block(currentBlockSize, 0.0);
        for (int i = 0; i < currentBlockSize; ++i) {
            block[i] = signal[start + i];
        }
        for (int i = 0; i < currentBlockSize; ++i) {
            for (int j = 0; j < kernelSize; ++j) {
                output[start + i + j] += block[i] * kernel[j];
            }
        }
    }
    return output;
}

// Normalized cross-correlation: scales correlation by local energies
double normalizedCrossCorrelation(const std::vector<double>& signal1, const std::vector<double>& signal2) {
    if (signal1.size() != signal2.size() || signal1.empty()) {
        return 0.0;
    }
    double mean1 = std::accumulate(signal1.begin(), signal1.end(), 0.0) / signal1.size();
    double mean2 = std::accumulate(signal2.begin(), signal2.end(), 0.0) / signal2.size();

    double num = 0.0;
    double den1 = 0.0;
    double den2 = 0.0;

    for (size_t i = 0; i < signal1.size(); ++i) {
        num += (signal1[i] - mean1) * (signal2[i] - mean2);
        den1 += (signal1[i] - mean1) * (signal1[i] - mean1);
        den2 += (signal2[i] - mean2) * (signal2[i] - mean2);
    }

    if (den1 == 0 || den2 == 0) {
        return 0.0;
    }

    return num / std::sqrt(den1 * den2);
}

// Apply a Hanning window to reduce spectral leakage
std::vector<double> applyHanningWindow(const std::vector<double>& signal) {
    int N = signal.size();
    std::vector<double> windowedSignal(N, 0.0);
    for (int n = 0; n < N; ++n) {
        double multiplier = 0.5 * (1.0 - std::cos((2.0 * M_PI * n) / (N - 1)));
        windowedSignal[n] = signal[n] * multiplier;
    }
    return windowedSignal;
}

// Block convolution for processing a segment with an impulse response
std::vector<double> blockConvolution(const std::vector<double>& block, const std::vector<double>& impulseResponse) {
    int blockSize = block.size();
    int impulseSize = impulseResponse.size();
    int outputSize = blockSize + impulseSize - 1;
    std::vector<double> output(outputSize, 0.0);
    for (int m = 0; m < outputSize; ++m) {
        for (int n = 0; n < impulseSize; ++n) {
            if (m - n >= 0 && m - n < blockSize) {
                output[m] += block[m - n] * impulseResponse[n];
            }
        }
    }
    return output;
}

// Apply a Hamming window to further mitigate boundary artifacts
std::vector<double> applyHammingWindow(const std::vector<double>& signal) {
    int N = signal.size();
    std::vector<double> windowedSignal(N, 0.0);
    for (int n = 0; n < N; ++n) {
        double multiplier = 0.54 - 0.46 * std::cos((2.0 * M_PI * n) / (N - 1));
        windowedSignal[n] = signal[n] * multiplier;
    }
    return windowedSignal;
}

int main() {
    std::vector<double> signal = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> kernel = {0.2, 0.5, 0.2};
    std::vector<double> pattern = {2.0, 3.0, 4.0};

    std::vector<double> convResult = convolve(signal, kernel);
    std::vector<double> corrResult = correlate(signal, kernel);
    std::vector<double> crossCorrResult = crossCorrelate(signal, kernel);
    std::vector<double> overlapAddResult = overlapAddConvolution(signal, kernel, 3);
    std::vector<double> normCrossCorrResult = normalizedCrossCorrelation(signal, pattern);
    std::vector<double> hanningResult = applyHanningWindow(signal);
    std::vector<double> blockConvResult = blockConvolution(signal, kernel);
    std::vector<double> hamstd::mingResult = applyHamstd::mingWindow(signal);

    std::cout << "Convolution Result: ";
    for (double val : convResult) { std::cout << val << " "; }
    std::cout << std::endl;

    std::cout << "Correlation Result: ";
    for (double val : corrResult) { std::cout << val << " "; }
    std::cout << std::endl;

    std::cout << "Cross-Correlation Result: ";
    for (double val : crossCorrResult) { std::cout << val << " "; }
    std::cout << std::endl;

    std::cout << "Overlap-Add Convolution Result: ";
    for (double val : overlapAddResult) { std::cout << val << " "; }
    std::cout << std::endl;

    std::cout << "Normalized Cross-Correlation Result: ";
    for (double val : normCrossCorrResult) { std::cout << val << " "; }
    std::cout << std::endl;

    std::cout << "Hanning Window Result: ";
    for (double val : hanningResult) { std::cout << val << " "; }
    std::cout << std::endl;

    std::cout << "Block Convolution Result: ";
    for (double val : blockConvResult) { std::cout << val << " "; }
    std::cout << std::endl;

    std::cout << "Hamstd::ming Window Result: ";
    for (double val : hamstd::mingResult) { std::cout << val << " "; }
    std::cout << std::endl;

    return 0;
}
