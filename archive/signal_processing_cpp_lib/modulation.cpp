#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>

// Define M_PI if not defined.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Function to map an integer symbol to its PSK representation.
std::complex<double> mapPSKSymbol(int symbol, int modulationOrder) {
  if (modulationOrder <= 0) {
      throw std::invalid_argument("Modulation order must be positive.");
  }
  if (symbol < 0 || symbol >= modulationOrder) {
    throw std::invalid_argument("Invalid symbol for PSK modulation");
  }
  double phase = 2.0 * M_PI * symbol / modulationOrder;
  return std::polar(1.0, phase);
}

// Function to map an integer symbol to its QAM constellation point.
std::complex<double> mapQAMSymbol(int symbol, int sqrtModOrder) {
  if (sqrtModOrder <= 0) {
      throw std::invalid_argument("Square root of modulation order must be positive.");
  }
  int modOrder = sqrtModOrder * sqrtModOrder;
  if (symbol < 0 || symbol >= modOrder) {
    throw std::invalid_argument("Invalid symbol for QAM modulation");
  }
  int row = symbol / sqrtModOrder;
  int col = symbol % sqrtModOrder;
  // Normalization factor to maintain unit average energy.
  double normalizationFactor = std::sqrt((2.0 * (modOrder - 1)) / 3.0);
  if (std::abs(normalizationFactor) < 1e-9) {
      return std::complex<double>(0.0, 0.0);
  }
  double I = (2 * col - sqrtModOrder + 1) / normalizationFactor;
  double Q = (2 * row - sqrtModOrder + 1) / normalizationFactor;
  return std::complex<double>(I, Q);
}

// Function to generate an FSK signal for a given symbol.
std::vector<double> generateFSKSignal(int symbol, double symbolDuration, double sampleRate,
    double baseFrequency, double frequencyDeviation) {
  if (sampleRate <= 0) {
      throw std::invalid_argument("Sample rate must be positive.");
  }
  int numSamples = static_cast<int>(symbolDuration * sampleRate);
  if (numSamples <= 0) {
      return {};
  }
  std::vector<double> signal(numSamples);
  // Calculate frequency for the given symbol.
  double frequency = baseFrequency + symbol * frequencyDeviation;
  for (int n = 0; n < numSamples; ++n) {
    double time = n / sampleRate;
    signal[n] = std::sin(2 * M_PI * frequency * time);
  }
  return signal;
}

// Function to demodulate a PSK symbol from a received complex sample.
int demodulatePSKSymbol(const std::complex<double>& received, int modulationOrder) {
  if (modulationOrder <= 0) {
      throw std::invalid_argument("Modulation order must be positive.");
  }
  // Compute the phase of the received signal.
  double phase = std::arg(received);
  if (phase < 0) {
    phase += 2 * M_PI;
  }
  double symbolInterval = 2 * M_PI / modulationOrder;
  int symbol = static_cast<int>(std::round(phase / symbolInterval)) % modulationOrder;
  return symbol;
}

// Function to demodulate a QAM symbol from a received complex sample.
int demodulateQAMSymbol(const std::complex<double>& received, int sqrtModOrder) {
  if (sqrtModOrder <= 0) {
      throw std::invalid_argument("Square root of modulation order must be positive.");
  }
  int modOrder = sqrtModOrder * sqrtModOrder;
  double normalizationFactor = std::sqrt((2.0 * (modOrder - 1)) / 3.0);
  if (std::abs(normalizationFactor) < 1e-9) {
      return 0;
  }
  double I = received.real();
  double Q = received.imag();
  int col = static_cast<int>(std::round((I * normalizationFactor + sqrtModOrder - 1) / 2));
  int row = static_cast<int>(std::round((Q * normalizationFactor + sqrtModOrder - 1) / 2));
  col = std::max(0, std::min(col, sqrtModOrder - 1));
  row = std::max(0, std::min(row, sqrtModOrder - 1));
  int symbol = row * sqrtModOrder + col;
  return symbol;
}

// Function to demodulate an FSK symbol based on zero crossing rate.
int demodulateFSKSymbol(const std::vector<double>& symbolSegment, double sampleRate,
    double baseFrequency, double frequencyDeviation, int modulationOrder) {
  if (symbolSegment.empty() || std::abs(frequencyDeviation) < 1e-9) {
      return 0;
  }
  int numSamples = symbolSegment.size();
  int zeroCrossings = 0;
  // Count zero crossings in the received segment.
  for (int n = 1; n < numSamples; ++n) {
    if (symbolSegment[n - 1] * symbolSegment[n] < 0) {
      ++zeroCrossings;
    }
  }
  // Estimate frequency from zero crossing count.
  double estimatedFrequency = (zeroCrossings * sampleRate) / (2.0 * numSamples);
  int symbol = static_cast<int>(std::round((estimatedFrequency - baseFrequency) / frequencyDeviation));
  symbol = std::max(0, std::min(symbol, modulationOrder - 1));
  return symbol;
}

// Function to correct carrier phase offsets in a received PSK symbol.
std::complex<double> correctCarrierPhase(const std::complex<double>& received, double estimatedPhaseOffset) {
  std::complex<double> correctionFactor = std::polar(1.0, -estimatedPhaseOffset);
  return received * correctionFactor;
}

// Function to identify the nearest symbol in a constellation based on Euclidean distance.
int findNearestSymbol(const std::complex<double>& received, const std::vector<std::complex<double>>& constellation) {
  double minDistance = std::numeric_limits<double>::max();
  int bestIndex = -1;
  for (std::size_t i = 0; i < constellation.size(); ++i) {
    double distance = std::norm(received - constellation[i]);
    if (distance < minDistance) {
      minDistance = distance;
      bestIndex = static_cast<int>(i);
    }
  }
  return bestIndex;
}

// Function to compute likelihoods for each symbol based on a Gaussian probability model.
std::vector<double> computeSymbolLikelihoods(const std::complex<double>& received, const std::vector<std::complex<double>>& constellation) {
  std::vector<double> likelihoods;
  likelihoods.reserve(constellation.size());
  double sigma = 1.0; // Assumed noise standard deviation.
  for (const auto &symbol : constellation) {
    double distance = std::norm(received - symbol);
    double likelihood = std::exp(-distance / (2 * sigma * sigma));
    likelihoods.push_back(likelihood);
  }
  return likelihoods;
}

// Function to equalize a received signal using a provided channel estimate.
std::vector<std::complex<double>> equalizeSignal(const std::vector<std::complex<double>>& received, const std::vector<std::complex<double>>& channelEstimate) {
  if (received.size() != channelEstimate.size()) {
      throw std::invalid_argument("Received signal and channel estimate must have the same size.");
  }
  std::vector<std::complex<double>> equalized(received.size());
  for (std::size_t i = 0; i < received.size(); ++i) {
    if (std::abs(channelEstimate[i]) > 1e-12) {
      equalized[i] = received[i] / channelEstimate[i];
    } else {
      equalized[i] = received[i];
    }
  }
  return equalized;
}

int main() {
  using namespace std;

  //=== PSK Modulation/Demodulation Demonstration ===//
  int modulationOrder = 8; // Example: 8-PSK
  cout << "PSK Constellation: " << endl;
  vector<complex<double>> pskConstellation;
  for (int i = 0; i < modulationOrder; ++i) {
    complex<double> symbol = mapPSKSymbol(i, modulationOrder);
    pskConstellation.push_back(symbol);
    cout << "Symbol " << i << ": " << symbol << endl;
  }

  complex<double> testPSK = pskConstellation[3];
  int demodulatedPSK = demodulatePSKSymbol(testPSK, modulationOrder);
  cout << "\nDemodulated PSK symbol: " << demodulatedPSK << endl;

  //=== QAM Modulation/Demodulation Demonstration ===//
  int sqrtModOrder = 4; // Example: 16-QAM (4x4 grid)
  cout << "\nQAM Constellation: " << endl;
  vector<complex<double>> qamConstellation;
  for (int i = 0; i < sqrtModOrder * sqrtModOrder; ++i) {
    complex<double> symbol = mapQAMSymbol(i, sqrtModOrder);
    qamConstellation.push_back(symbol);
    cout << "Symbol " << i << ": " << symbol << endl;
  }
  complex<double> testQAM = qamConstellation[10];
  int demodulatedQAM = demodulateQAMSymbol(testQAM, sqrtModOrder);
  cout << "\nDemodulated QAM symbol: " << demodulatedQAM << endl;

  //=== FSK Modulation/Demodulation Demonstration ===//
  int fskSymbol = 2;
  double symbolDuration = 0.01;
  double sampleRate = 10000;
  double baseFrequency = 1000;
  double frequencyDeviation = 100;
  vector<double> fskSignal = generateFSKSignal(fskSymbol, symbolDuration, sampleRate, baseFrequency, frequencyDeviation);
  int demodulatedFSK = demodulateFSKSymbol(fskSignal, sampleRate, baseFrequency, frequencyDeviation, 8);
  cout << "\nDemodulated FSK symbol: " << demodulatedFSK << endl;

  //=== Carrier Phase Correction Demonstration ===//
  double phaseOffset = 0.2;
  complex<double> receivedPSK = pskConstellation[5] * std::polar(1.0, phaseOffset);
  complex<double> correctedSymbol = correctCarrierPhase(receivedPSK, phaseOffset);
  cout << "\nOriginal PSK symbol: " << pskConstellation[5] << endl;
  cout << "Received symbol with phase offset: " << receivedPSK << endl;
  cout << "Corrected symbol: " << correctedSymbol << endl;

  //=== Euclidean Distance Based Detection Demonstration ===//
  int nearestSymbolIndex = findNearestSymbol(receivedPSK, pskConstellation);
  cout << "\nNearest PSK symbol index based on Euclidean distance: " << nearestSymbolIndex << endl;

  //=== Symbol Likelihood Computation Demonstration ===//
  vector<double> likelihoods = computeSymbolLikelihoods(receivedPSK, pskConstellation);
  cout << "\nSymbol Likelihoods: " << endl;
  for (size_t i = 0; i < likelihoods.size(); ++i) {
    cout << "Symbol " << i << ": " << likelihoods[i] << endl;
  }

  //=== Signal Equalization Demonstration ===//
  cout << "\nSignal Equalization Demonstration: " << endl;
  vector<complex<double>> receivedSignal = { {1.0, 0.5}, {0.8, 0.3}, {1.2, -0.2} };
  vector<complex<double>> channelEstimate = { {0.9, 0.1}, {0.95, 0.05}, {1.0, 0.0} };
  vector<complex<double>> equalizedSignal = equalizeSignal(receivedSignal, channelEstimate);
  for (size_t i = 0; i < equalizedSignal.size(); ++i) {
    cout << "Equalized sample " << i << ": " << equalizedSignal[i] << endl;
  }

  return 0;
}
