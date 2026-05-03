#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

using std::vector;
using std::complex;
using std::cout;
using std::endl;
using std::sin;
using std::cos;
using std::abs;
using std::polar;
using std::size_t;

// Define M_PI if it is not already defined.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Generate the Hann window of a given size.
std::vector<double> generateHannWindow(std::size_t windowSize) {
  if (windowSize <= 1) {
      return {1.0};
  }
  std::vector<double> window(windowSize, 0.0);
  for (std::size_t n = 0; n < windowSize; ++n) {
    window[n] = 0.5 * (1 - std::cos((2 * M_PI * n) / (windowSize - 1)));
  }
  return window;
}

// Extract a segment from the signal starting at 'start' with length 'windowSize'.
std::vector<double> extractSegment(const std::vector<double>& signal, std::size_t start, std::size_t windowSize) {
  std::vector<double> segment(windowSize, 0.0);
  for (std::size_t i = 0; i < windowSize && (start + i) < signal.size(); ++i) {
    segment[i] = signal[start + i];
  }
  return segment;
}

// Apply a window to a signal segment via elementwise multiplication.
std::vector<double> applyWindow(const std::vector<double>& segment, const std::vector<double>& window) {
  if (segment.size() != window.size()) {
      throw std::invalid_argument("Segment and window must have the same size.");
  }
  std::vector<double> windowedSegment(segment.size());
  for (std::size_t i = 0; i < segment.size(); ++i) {
    windowedSegment[i] = segment[i] * window[i];
  }
  return windowedSegment;
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

// Iterative Cooley-Tukey FFT implementation (radix-2 FFT)
void iterativeFFT(std::vector<std::complex<double>>& data) {
    const size_t N = data.size();
    if (N <= 1) return;

    // Bit-reversal permutation
    bitReversalPermutation(data);

    // Iterative FFT
    for (size_t len = 2; len <= N; len <<= 1) {
        double angle = -2.0 * M_PI / len;
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

// Perform a naive Discrete Fourier Transform (DFT) on the input signal segment.
std::vector<std::complex<double>> performFFT(const std::vector<double>& input) {
  std::size_t N = input.size();
  if (N == 0) {
      return {};
  }
  std::vector<std::complex<double>> complexSignal;
  for (const auto& value : input) {
      complexSignal.push_back(std::complex<double>(value, 0.0));
  }
  iterativeFFT(complexSignal);
  return complexSignal;
}

// Compute the magnitude spectrum from the complex FFT result.
std::vector<double> computeMagnitudeSpectrum(const std::vector<std::complex<double>>& fftResult) {
  std::vector<double> magnitude(fftResult.size());
  for (std::size_t i = 0; i < fftResult.size(); ++i) {
    magnitude[i] = std::abs(fftResult[i]);
  }
  return magnitude;
}

// Detect peaks in the magnitude spectrum that exceed a given threshold.
std::vector<std::size_t> detectPeaks(const std::vector<double>& magnitude, double threshold) {
  std::vector<std::size_t> peakIndices;
  if (magnitude.size() < 3) {
      return peakIndices;
  }
  for (std::size_t i = 1; i < magnitude.size() - 1; ++i) {
    if (magnitude[i] > threshold && magnitude[i] > magnitude[i - 1] && magnitude[i] > magnitude[i + 1]) {
      peakIndices.push_back(i);
    }
  }
  return peakIndices;
}

// Main function: demonstrates the STFT processing pipeline.
int main() {
  const std::size_t signalLength = 1024;
  std::vector<double> signal(signalLength, 0.0);
  double sampleRate = 100.0;
  double frequency1 = 5.0;
  double frequency2 = 20.0;

  for (std::size_t n = 0; n < signalLength; ++n) {
    double t = static_cast<double>(n) / sampleRate;
    signal[n] = std::sin(2 * M_PI * frequency1 * t) + 0.5 * std::sin(2 * M_PI * frequency2 * t);
  }

  std::size_t windowSize = 128;
  std::size_t hopSize = windowSize / 2;
  std::vector<double> hannWindow = generateHannWindow(windowSize);

  for (std::size_t start = 0; start < signal.size(); start += hopSize) {
    std::vector<double> segment = extractSegment(signal, start, windowSize);
    std::vector<double> windowedSegment = applyWindow(segment, hannWindow);
    std::vector<std::complex<double>> fftResult = performFFT(windowedSegment);
    std::vector<double> magnitude = computeMagnitudeSpectrum(fftResult);
    const double PEAK_THRESHOLD = 10.0;
double threshold = PEAK_THRESHOLD;
    std::vector<std::size_t> peaks = detectPeaks(magnitude, threshold);

    std::cout << "Segment starting at sample " << start << ":\n";
    std::cout << "Detected peak frequency bins: ";
    for (std::size_t idx : peaks) {
      std::cout << idx << " ";
    }
    std::cout << "\n\n";
  }

  return 0;
}
