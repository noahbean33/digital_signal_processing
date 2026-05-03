#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <numeric>
#include <stdexcept>

// Predicts the next sample using past samples and autoregressive
// coefficients.
double predictNextSample(const std::vector<double>& pastSamples,
  const std::vector<double>& arCoeffs) {
  if (pastSamples.size() < arCoeffs.size()) {
      return 0.0; // Not enough data to make a prediction
  }
  double prediction = 0.0;
  int order = arCoeffs.size();
  for (int i = 0; i < order; ++i) {
    // Uses the most recent p samples weighted by their
    // coefficients.
    prediction += arCoeffs[i] * pastSamples[pastSamples.size() -
      i - 1];
  }
  return prediction;
}

double movingAverageFilter(const std::vector<double>& signal, int currentIndex, int windowSize) {
  double sum = 0.0;
  int count = 0;
  for (int i = currentIndex; i > currentIndex - windowSize && i >= 0; --i) {
    sum += signal[i];
    ++count;
  }
  return (count > 0) ? (sum / count) : 0.0;
}

double computeMeanSquaredError(const std::vector<double>& predictions, const std::vector<double>& actual) {
  if (predictions.size() != actual.size() || predictions.empty()) {
    return 0.0; // Return 0 if inputs are mismatched or empty
  }
  double mse = 0.0;
  int n = predictions.size();
  for (int i = 0; i < n; ++i) {
    double error = predictions[i] - actual[i];
    mse += error * error;
  }
  return mse / n;
}

std::vector<double> updateFilterCoefficients(const std::vector<double>& coeffs, const std::vector<double>& pastSamples, double error, double learningRate) {
  std::vector<double> updatedCoeffs(coeffs.size());
  for (size_t i = 0; i < coeffs.size(); ++i) {
    updatedCoeffs[i] = coeffs[i] - learningRate * error * pastSamples[pastSamples.size() - i - 1];
  }
  return updatedCoeffs;
}

std::vector<double> updateRLSCoefficients(const std::vector<double>& coeffs, const std::vector<double>& gain, double error) {
  std::vector<double> updatedCoeffs(coeffs.size());
  for (size_t i = 0; i < coeffs.size(); ++i) {
    updatedCoeffs[i] = coeffs[i] + gain[i] * error;
  }
  return updatedCoeffs;
}

bool hasConverged(const std::vector<double>& mseHistory, double tolerance) {
  double averageMSE = 0.0;
  if (mseHistory.empty()) return false;
  for (double mse : mseHistory) {
    averageMSE += mse;
  }
  averageMSE /= mseHistory.size();
  return averageMSE < tolerance;
}

// Computes the dot product of two vectors using OpenMP for
// parallelization.
double parallelDotProduct(const std::vector<double>& a, const
  std::vector<double>& b) {
  if (a.size() != b.size()) {
      return 0.0; // Vectors must have the same size
  }
  double result = 0.0;
  #pragma omp parallel for reduction(+:result)
  for (size_t i = 0; i < a.size(); ++i) {
    result += a[i] * b[i];
  }
  return result;
}

int main() {
  std::vector<double> signal = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  std::vector<double> arCoeffs = {0.5, -0.2};

  double prediction = predictNextSample(signal, arCoeffs);
  std::cout << "Predicted value: " << prediction << std::endl;

  double filteredValue = movingAverageFilter(signal, 5, 3);
  std::cout << "Moving Average Filter Output: " << filteredValue << std::endl;

  std::vector<double> predictions = {2.1, 3.9, 5.8};
  std::vector<double> actual = {2.0, 4.0, 6.0};
  double mse = computeMeanSquaredError(predictions, actual);
  std::cout << "Mean Squared Error: " << mse << std::endl;

  double error = actual.back() - prediction;
  const double LEARNING_RATE = 0.01;
double learningRate = LEARNING_RATE;
  std::vector<double> updatedCoeffs = updateFilterCoefficients(arCoeffs, signal, error, learningRate);
  std::cout << "Updated AR Coefficients (Gradient Descent): ";
  for (double coeff : updatedCoeffs) { std::cout << coeff << " "; }
  std::cout << std::endl;

  const std::vector<double> RLS_GAIN = {0.1, 0.1};
std::vector<double> gain = RLS_GAIN;
  std::vector<double> rlsCoeffs = updateRLSCoefficients(arCoeffs, gain, error);
  std::cout << "Updated AR Coefficients (RLS): ";
  for (double coeff : rlsCoeffs) { std::cout << coeff << " "; }
  std::cout << std::endl;

  std::vector<double> mseHistory = {0.05, 0.04, 0.03, 0.04};
  const double CONVERGENCE_TOLERANCE = 0.05;
double tolerance = CONVERGENCE_TOLERANCE;
  if (hasConverged(mseHistory, tolerance)) {
    std::cout << "The adaptive filter has converged." << std::endl;
  } else {
    std::cout << "The adaptive filter has not converged yet." << std::endl;
  }

  std::vector<double> vecA = {1.0, 2.0, 3.0};
  std::vector<double> vecB = {4.0, 5.0, 6.0};
  double dotProd = parallelDotProduct(vecA, vecB);
  std::cout << "Parallel Dot Product: " << dotProd << std::endl;

  return 0;
}
