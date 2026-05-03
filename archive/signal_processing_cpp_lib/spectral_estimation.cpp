#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include <random>

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
        double angle = -2.0 * acos(-1.0) / len;
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

// Naive Discrete Fourier Transform (DFT) implementation
vector<complex<double>> performFFT(const vector<double>& signal) {
  size_t N = signal.size();
  if (N == 0) {
      return {};
  }
  vector<complex<double>> complexSignal;
  for (const auto& value : signal) {
      complexSignal.push_back(std::complex<double>(value, 0.0));
  }
  iterativeFFT(complexSignal);
  return complexSignal;
}

// Compute periodogram (Power Spectral Density) using FFT
vector<double> computePeriodogram(const vector<double>& signal) {
  size_t N = signal.size();
  if (N == 0) {
      return {};
  }
  vector<double> periodogram(N, 0.0);
  vector<complex<double>> X = performFFT(signal);
  for (size_t k = 0; k < N; ++k) {
    periodogram[k] = norm(X[k]) / static_cast<double>(N);
  }
  return periodogram;
}

// Estimate AR Coefficients using Levinson-Durbin recursion
vector<double> estimateARCoefficients(const vector<double>& signal, int order) {
  size_t N = signal.size();
  if (N == 0 || order <= 0) {
      return {};
  }
  vector<double> r(order + 1, 0.0);
  for (int lag = 0; lag <= order; ++lag) {
    for (size_t n = lag; n < N; ++n) {
      r[lag] += signal[n] * signal[n - lag];
    }
    r[lag] /= static_cast<double>(N);
  }
  vector<double> a(order, 0.0);
  vector<double> E(order + 1, 0.0);
  E[0] = r[0];
  for (int k = 1; k <= order; ++k) {
    double lambda = 0.0;
    if (std::abs(E[k - 1]) < 1e-9) {
        return {};
    }
    for (int j = 1; j < k; ++j) {
      lambda -= a[j - 1] * r[k - j];
    }
    lambda /= E[k - 1];

    a[k - 1] = lambda;
    for (int j = 1; j < k; ++j) {
      a[j - 1] = a[j - 1] - lambda * a[k - j - 1];
    }
    E[k] = E[k - 1] * (1.0 - lambda * lambda);
  }
  return a;
}

// Compute Akaike Information Criterion (AIC)
double computeAIC(double errorVariance, int order, int numSamples) {
  return 2.0 * order + numSamples * log(errorVariance);
}

// Compute AR-based Power Spectral Density (PSD)
vector<double> computeARPSD(const vector<double>& arCoeffs, double noiseVariance, const vector<double>& frequencies) {
  size_t order = arCoeffs.size();
  vector<double> psd(frequencies.size(), 0.0);
  for (size_t i = 0; i < frequencies.size(); ++i) {
    complex<double> denominator = 1.0;
    for (size_t k = 0; k < order; ++k) {
      denominator += arCoeffs[k] * exp(complex<double>(0, -frequencies[i] * (k + 1)));
    }
    double magnitudeSquared = norm(denominator);
    if (std::abs(magnitudeSquared) < 1e-9) {
        psd[i] = 0.0;
    } else {
        psd[i] = noiseVariance / magnitudeSquared;
    }
  }
  return psd;
}

// Compute Welch Periodogram for spectral estimation
vector<double> computeWelchPeriodogram(const vector<double>& signal, int segmentLength, int overlap, const vector<double>& window) {
  if (segmentLength <= overlap) {
      return {};
  }
  int step = segmentLength - overlap;
  if (step == 0) {
      return {};
  }
  int numSegments = (signal.size() - overlap) / step;
  if (numSegments <= 0) {
      return {};
  }
  int fftSize = segmentLength;
  vector<double> welchPSD(fftSize, 0.0);
  double windowEnergy = accumulate(window.begin(), window.end(), 0.0, [](double sum, double val) { return sum + val * val; });
  if (std::abs(windowEnergy) < 1e-9) {
      return {};
  }

  for (int i = 0; i < numSegments; ++i) {
    int start = i * step;
    vector<double> segment(signal.begin() + start, signal.begin() + start + segmentLength);
    for (int j = 0; j < segmentLength; ++j) {
      segment[j] *= window[j];
    }
    vector<complex<double>> spectrum = performFFT(segment);
    for (int k = 0; k < fftSize; ++k) {
      welchPSD[k] += norm(spectrum[k]);
    }
  }
  for (int k = 0; k < fftSize; ++k) {
    welchPSD[k] /= (numSegments * windowEnergy);
  }
  return welchPSD;
}

// Estimate AR Coefficients using Burg Algorithm
vector<double> estimateARBurg(const vector<double>& signal, int order) {
  int N = signal.size();
  if (N == 0 || order <= 0) {
      return {};
  }
  vector<double> a(order, 0.0);
  vector<double> ef(signal);
  vector<double> eb(signal);
  double error = 0.0;
  for (int i = 0; i < N; ++i) {
    error += signal[i] * signal[i];
  }
  error /= N;
  for (int m = 0; m < order; ++m) {
    double num = 0.0, den = 0.0;
    for (int n = m + 1; n < N; ++n) {
      num += ef[n] * eb[n - 1];
      den += ef[n] * ef[n] + eb[n - 1] * eb[n - 1];
    }
    if (std::abs(den) < 1e-9) {
        return {};
    }
    double k = -2.0 * num / den;
    a[m] = k;
    for (int j = 0; j < m; ++j) {
      a[j] = a[j] + k * a[m - j - 1];
    }
    for (int n = N - 1; n > m; --n) {
      double temp = ef[n];
      ef[n] = ef[n] + k * eb[n - 1];
      eb[n] = eb[n - 1] + k * temp;
    }
    error *= (1.0 - k * k);
  }
  return a;
}

// Select Optimal AR Model Order using AIC
int selectOptimalAROrder(const vector<double>& signal, int maxOrder) {
  int N = signal.size();
  if (N == 0 || maxOrder <= 0) {
      return 1;
  }
  int optimalOrder = 1;
  double bestAIC = numeric_limits<double>::infinity();
  for (int order = 1; order <= maxOrder; ++order) {
    if (N <= order) {
        break;
    }
    vector<double> arCoeffs = estimateARCoefficients(signal, order);
    if (arCoeffs.empty()) {
        continue;
    }
    double errorVariance = 0.0;
    for (int n = order; n < N; ++n) {
      double prediction = 0.0;
      for (int k = 0; k < order; ++k) {
        prediction += arCoeffs[k] * signal[n - k - 1];
      }
      double error = signal[n] + prediction;
      errorVariance += error * error;
    }
    errorVariance /= (N - order);
    if (errorVariance <= 0) {
        continue;
    }
    double currentAIC = 2.0 * order + N * log(errorVariance);
    if (currentAIC < bestAIC) {
      bestAIC = currentAIC;
      optimalOrder = order;
    }
  }
  return optimalOrder;
}

int main() {
      std::mt19937 rng(static_cast<unsigned int>(time(nullptr)));
    std::uniform_real_distribution<double> dist(-0.05, 0.05);

  const int N = 128;
  vector<double> signal(N, 0.0);
  const double PI = acos(-1.0);
  double frequency = 5.0;
  for (int n = 0; n < N; ++n) {
    signal[n] = sin(2 * PI * frequency * n / N) + dist(rng);
  }

  vector<double> periodogram = computePeriodogram(signal);
  std::cout << "Periodogram (first 10 values):" << std::endl;
  for (size_t i = 0; i < 10 && i < periodogram.size(); ++i) { std::cout << periodogram[i] << " "; }
  std::cout << "\n\n";

  int arOrder = 4;
  vector<double> arCoeffs = estimateARCoefficients(signal, arOrder);
  std::cout << "Estimated AR Coefficients (Levinson-Durbin, order " << arOrder << "):";
  for (size_t i = 0; i < arCoeffs.size(); ++i) { std::cout << arCoeffs[i] << " "; }
  std::cout << "\n\n";

  double errorVariance = 0.0;
  for (int n = arOrder; n < N; ++n) {
    double prediction = 0.0;
    for (int k = 0; k < arOrder; ++k) {
      prediction += arCoeffs[k] * signal[n - k - 1];
    }
    double error = signal[n] + prediction;
    errorVariance += error * error;
  }
  errorVariance /= (N - arOrder);
  std::cout << "Error Variance: " << errorVariance << "\n\n";

  double aicValue = computeAIC(errorVariance, arOrder, N);
  std::cout << "AIC Value: " << aicValue << "\n\n";

  int numFreqs = 65;
  vector<double> frequencies(numFreqs, 0.0);
  for (int i = 0; i < numFreqs; ++i) { frequencies[i] = i * PI / (numFreqs - 1); }
  vector<double> arPSD = computeARPSD(arCoeffs, errorVariance, frequencies);
  std::cout << "AR-based PSD (first 10 values):" << std::endl;
  for (size_t i = 0; i < 10 && i < arPSD.size(); ++i) { std::cout << arPSD[i] << " "; }
  std::cout << "\n\n";

  int segmentLength = 64;
  int overlap = 32;
  vector<double> window(segmentLength, 0.0);
  for (int n = 0; n < segmentLength; ++n) { window[n] = 0.54 - 0.46 * cos(2 * PI * n / (segmentLength - 1)); }
  vector<double> welchPSD = computeWelchPeriodogram(signal, segmentLength, overlap, window);
  std::cout << "Welch Periodogram (first 10 values):" << std::endl;
  for (size_t i = 0; i < 10 && i < welchPSD.size(); ++i) { std::cout << welchPSD[i] << " "; }
  std::cout << "\n\n";

  vector<double> burgCoeffs = estimateARBurg(signal, arOrder);
  std::cout << "Estimated AR Coefficients (Burg Algorithm, order " << arOrder << "):";
  for (size_t i = 0; i < burgCoeffs.size(); ++i) { std::cout << burgCoeffs[i] << " "; }
  std::cout << "\n\n";

  int maxOrder = 10;
  int optimalOrder = selectOptimalAROrder(signal, maxOrder);
  std::cout << "Optimal AR Model Order (via AIC): " << optimalOrder << std::endl;

  return 0;
}
