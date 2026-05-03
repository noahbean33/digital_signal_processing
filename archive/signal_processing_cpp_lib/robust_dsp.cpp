#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>
#include <thread>
#include <mutex>
#include <stdexcept>

// Huber weight function
double huberWeight(double residual, double threshold) {
  if (std::fabs(residual) <= threshold) {
    return 1.0;
  }
  return threshold / std::fabs(residual);
}

// Function to compute the median of a vector.
double computeMedian(std::vector<double> window) {
  if (window.empty()) {
      return 0.0;
  }
  std::sort(window.begin(), window.end());
  std::size_t n = window.size();
  if (n % 2 == 0) {
    return (window[n/2 - 1] + window[n/2]) / 2.0;
  }
  return window[n/2];
}

// Computes a robust confidence interval using the median absolute deviation (MAD).
std::pair<double, double> robustConfidenceInterval(const std::vector<double>& data, double confidenceFactor) {
  if (data.empty()) {
      return {0.0, 0.0};
  }
  std::vector<double> sortedData = data;
  std::sort(sortedData.begin(), sortedData.end());
  double median = sortedData[sortedData.size() / 2];
  std::vector<double> deviations;
  deviations.reserve(data.size());
  for (auto value : sortedData) {
    deviations.push_back(std::fabs(value - median));
  }
  std::sort(deviations.begin(), deviations.end());
  double mad = deviations[deviations.size() / 2];
  double lowerBound = median - confidenceFactor * mad;
  double upperBound = median + confidenceFactor * mad;
  return {lowerBound, upperBound};
}

// Utility function to compute the Median Absolute Deviation (MAD).
double computeMAD(const std::vector<double>& data) {
  if (data.empty()) {
      return 0.0;
  }
  std::vector<double> sortedData = data;
  std::sort(sortedData.begin(), sortedData.end());
  double median = sortedData[sortedData.size() / 2];
  std::vector<double> deviations;
  deviations.reserve(data.size());
  for (double value : sortedData) {
    deviations.push_back(std::fabs(value - median));
  }
  std::sort(deviations.begin(), deviations.end());
  return deviations[deviations.size() / 2];
}

// Implements one iteration update of an IRLS (Iterative Reweighted Least Squares) algorithm
double irlsUpdate(const std::vector<double>& x, const std::vector<double>& y, double beta, double threshold) {
  if (x.size() != y.size() || x.empty()) {
      throw std::invalid_argument("Input vectors must have the same size and not be empty.");
  }
  double numerator = 0.0;
  double denominator = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
      double residual = y[i] - beta * x[i];
      double weight = huberWeight(residual, threshold);
      numerator += weight * y[i] * x[i];
      denominator += weight * x[i] * x[i];
  }
  return (denominator > 1e-9) ? numerator / denominator : beta;
}

// Implements a robust moving median filter.
std::vector<double> robustMovingMedian(const std::vector<double>& signal, std::size_t windowSize) {
  if (windowSize % 2 == 0 || windowSize < 3 || signal.size() < windowSize) {
    return {};
  }
  std::vector<double> filtered;
  filtered.reserve(signal.size() - windowSize + 1);
  std::size_t halfWindow = windowSize / 2;
  for (std::size_t i = halfWindow; i < signal.size() - halfWindow; ++i) {
    std::vector<double> window(signal.begin() + i - halfWindow, signal.begin() + i + halfWindow + 1);
    filtered.push_back(computeMedian(window));
  }
  return filtered;
}

// Computes a robust sum of weighted residuals in a thread-safe manner.
double threadSafeRobustSum(const std::vector<double>& residuals) {
    double sum = 0.0;
    std::mutex mtx;
    auto worker = [&](std::size_t start, std::size_t end) {
      double localSum = 0.0;
      for (std::size_t i = start; i < end; ++i) {
        double weight = huberWeight(residuals[i], 1.0); // Using threshold of 1.0
        localSum += weight * residuals[i];
      }
      std::lock_guard<std::mutex> lock(mtx);
      sum += localSum;
    };
    
    std::size_t numThreads = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;
    std::size_t blockSize = (residuals.size() + numThreads - 1) / numThreads;
    
    for (std::size_t i = 0; i < numThreads; ++i) {
      std::size_t start = i * blockSize;
      if (start >= residuals.size()) break;
      std::size_t end = std::min(start + blockSize, residuals.size());
      threads.emplace_back(worker, start, end);
    }
    
    for (auto &t : threads) {
      t.join();
    }
    return sum;
}

int main() {
    std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<double> y = {2.1, 4.1, 6.0, 7.8, 10.2, 12.1, 14.2, 16.5, 18.1, 20.0};
    const double INITIAL_BETA = 2.0;
    const double IRLS_THRESHOLD = 1.5;
    double updatedBeta = irlsUpdate(x, y, INITIAL_BETA, IRLS_THRESHOLD);
    std::cout << "Updated beta from IRLS: " << updatedBeta << std::endl;

    std::vector<double> data = {5.0, 7.0, 6.5, 10.0, 4.0, 8.5, 6.0, 7.5, 9.0};
    const double CONFIDENCE_FACTOR = 1.4826; // For Gaussian-like data
    auto ci = robustConfidenceInterval(data, CONFIDENCE_FACTOR);
    std::cout << "Robust Confidence Interval: [" << ci.first << ", " << ci.second << "]" << std::endl;
      
    std::vector<double> window = {3, 1, 4, 1, 5};
    double medianValue = computeMedian(window);
    std::cout << "Median of window: " << medianValue << std::endl;
    
    double madValue = computeMAD(data);
    std::cout << "Median Absolute Deviation (MAD): " << madValue << std::endl;

    std::vector<double> signal = {1, 3, 2, 8, 7, 10, 9, 5, 6, 4, 3, 2};
    const std::size_t WINDOW_SIZE = 3;
    std::vector<double> filteredSignal = robustMovingMedian(signal, WINDOW_SIZE);
    std::cout << "Filtered Signal (Moving Median): ";
    for (double val : filteredSignal) { std::cout << val << " "; }
    std::cout << std::endl;

    std::vector<double> residuals = {0.5, 2.0, 3.0, 0.2, 1.5, 0.8};
    double robustSum = threadSafeRobustSum(residuals);
    std::cout << "Thread Safe Robust Sum: " << robustSum << std::endl;

    return 0;
}
