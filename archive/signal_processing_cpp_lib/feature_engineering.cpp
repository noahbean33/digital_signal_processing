#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <thread>
#include <stdexcept>

// Compute the total energy of a signal.
double computeSignalEnergy(const std::vector<double>& signal) {
  if (signal.empty()) {
      return 0.0;
  }
  double energy = 0.0;
  for (double sample : signal) {
    energy += sample * sample;
  }
  return energy;
}

// Calculate the Euclidean distance between two feature vectors.
double computeEuclideanDistance(const std::vector<double>& vec1, const std::vector<double>& vec2) {
  if (vec1.size() != vec2.size()) {
      throw std::invalid_argument("Vectors must have the same size.");
  }
  double sum = 0.0;
  for (std::size_t i = 0; i < vec1.size(); ++i) {
    double diff = vec1[i] - vec2[i];
    sum += diff * diff;
  }
  return std::sqrt(sum);
}

// Normalize a vector of features to the range [0, 1].
std::vector<double> normalizeFeatures(const std::vector<double>& features) {
  if (features.empty()) {
      return {};
  }
  double minVal = *std::min_element(features.begin(), features.end());
  double maxVal = *std::max_element(features.begin(), features.end());
  double range = maxVal - minVal;
  if (range < 1e-9) {
      return std::vector<double>(features.size(), 0.0);
  }
  std::vector<double> normalized;
  normalized.reserve(features.size());
  for (const double value : features) {
    normalized.push_back((value - minVal) / range);
  }
  return normalized;
}

// Fuse two feature vectors into a single vector.
std::vector<double> fuseFeatures(const std::vector<double>& timeFeatures, const std::vector<double>& frequencyFeatures) {
  std::vector<double> fused;
  fused.reserve(timeFeatures.size() + frequencyFeatures.size());
  fused.insert(fused.end(), timeFeatures.begin(), timeFeatures.end());
  fused.insert(fused.end(), frequencyFeatures.begin(), frequencyFeatures.end());
  return fused;
}

// Perform a single gradient descent step with L2 regularization.
double gradientDescentStep(double currentParam, double gradient, double learningRate, double regularization) {
  return currentParam - learningRate * (gradient + regularization * currentParam);
}

// Parallel normalization of features using multiple threads.
std::vector<double> parallelNormalizeFeatures(std::vector<double> features) {
  if (features.empty()) {
      return {};
  }
  std::size_t numThreads = std::thread::hardware_concurrency();
  std::size_t dataSize = features.size();
  if (numThreads == 0) numThreads = 2; // Fallback

  double minVal = *std::min_element(features.begin(), features.end());
  double maxVal = *std::max_element(features.begin(), features.end());
  double range = maxVal - minVal;
  if (range < 1e-9) {
      std::fill(features.begin(), features.end(), 0.0);
      return features;
  }

  std::vector<std::thread> threads;
  std::size_t blockSize = (dataSize + numThreads - 1) / numThreads;

  for (std::size_t i = 0; i < numThreads; ++i) {
    std::size_t start = i * blockSize;
    if (start >= dataSize) break;
    std::size_t end = std::min(start + blockSize, dataSize);
    threads.emplace_back([&features, minVal, range, start, end]() {
      for (std::size_t j = start; j < end; ++j) {
        features[j] = (features[j] - minVal) / range;
      }
    });
  }
  for (auto& t : threads) {
    t.join();
  }
  return features;
}

// Compute the Root Mean Square Error (RMSE).
double computeRMSE(const std::vector<double>& reference, const std::vector<double>& reconstructed) {
  if (reference.size() != reconstructed.size() || reference.empty()) {
      return 0.0;
  }
  double sumSqError = 0.0;
  std::size_t n = reference.size();
  for (std::size_t i = 0; i < n; ++i) {
    double error = reference[i] - reconstructed[i];
    sumSqError += error * error;
  }
  return std::sqrt(sumSqError / n);
}

// Compute the F1 score.
double computeF1Score(int truePos, int falsePos, int falseNeg) {
  double precision = (truePos + falsePos > 0) ? static_cast<double>(truePos) / (truePos + falsePos) : 0.0;
  double recall = (truePos + falseNeg > 0) ? static_cast<double>(truePos) / (truePos + falseNeg) : 0.0;
  if (precision + recall < 1e-9) {
      return 0.0;
  }
  return 2.0 * precision * recall / (precision + recall);
}

// Compute the variance of a dataset.
double computeVariance(const std::vector<double>& data) {
  if (data.empty()) {
      return 0.0;
  }
  double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
  double variance = 0.0;
  for (double value : data) {
    variance += (value - mean) * (value - mean);
  }
  return variance / data.size();
}

int main() {
  std::vector<double> signal = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> frequencySignal = {0.5, 1.5, 2.5, 3.5, 4.5};

  double energy = computeSignalEnergy(signal);
  std::cout << "Signal Energy: " << energy << std::endl;

  double distance = computeEuclideanDistance(signal, frequencySignal);
  std::cout << "Euclidean Distance: " << distance << std::endl;

  std::vector<double> normalizedSignal = normalizeFeatures(signal);
  std::cout << "Normalized Signal: ";
  for (double val : normalizedSignal) { std::cout << val << " "; }
  std::cout << std::endl;

  std::vector<double> fusedFeatures = fuseFeatures(signal, frequencySignal);
  std::cout << "Fused Features: ";
  for (double val : fusedFeatures) { std::cout << val << " "; }
  std::cout << std::endl;

  const double CURRENT_PARAM = 0.5;
  const double GRADIENT = 0.1;
  const double LEARNING_RATE = 0.05;
  const double REGULARIZATION = 0.01;
  double updatedParam = gradientDescentStep(CURRENT_PARAM, GRADIENT, LEARNING_RATE, REGULARIZATION);
  std::cout << "Updated Parameter after Gradient Descent: " << updatedParam << std::endl;

  std::vector<double> parallelNormalized = parallelNormalizeFeatures(signal);
  std::cout << "Parallel Normalized Signal: ";
  for (double val : parallelNormalized) { std::cout << val << " "; }
  std::cout << std::endl;

  double rmse = computeRMSE(signal, frequencySignal);
  std::cout << "RMSE: " << rmse << std::endl;

  const int TRUE_POS = 40, FALSE_POS = 10, FALSE_NEG = 5;
  double f1Score = computeF1Score(TRUE_POS, FALSE_POS, FALSE_NEG);
  std::cout << "F1 Score: " << f1Score << std::endl;

  double variance = computeVariance(signal);
  std::cout << "Variance: " << variance << std::endl;

  return 0;
}
