#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

// Function to compute the excess kurtosis of a signal vector.
double computeKurtosis(const std::vector<double>& signal) {
  int n = signal.size();
  if (n == 0) {
      return 0.0;
  }
  double mean = 0.0;
  for (const double &val : signal) {
    mean += val;
  }
  mean /= n;

  double m2 = 0.0, m4 = 0.0;
  for (const double &val : signal) {
    double diff = val - mean;
    m2 += diff * diff;
    m4 += diff * diff * diff * diff;
  }
  m2 /= n;
  m4 /= n;

  if (m2 < 1e-9) {
      return 0.0;
  }
  return m4 / (m2 * m2) - 3.0;
}

// FastICA update: perform one iteration for a weight vector w
Eigen::VectorXd fastICAUpdate(const Eigen::MatrixXd &w, const Eigen::MatrixXd &X) {
  int n = X.cols();
  if (n == 0) {
      return w;
  }
  Eigen::RowVectorXd wx = w.transpose() * X;
  Eigen::RowVectorXd gwx = wx.array().tanh();
  Eigen::RowVectorXd g_wx_deriv = 1.0 - gwx.array().square();
  Eigen::VectorXd new_w = (X * gwx.transpose()) / n - g_wx_deriv.mean() * w;
  new_w.normalize();
  return new_w;
}

// Adaptive FastICA update that incorporates a learning rate.
Eigen::VectorXd adaptiveFastICAUpdate(const Eigen::VectorXd &w,
                                    const Eigen::MatrixXd &X,
                                    double learningRate) {
  int n = X.cols();
  if (n == 0) {
      return w;
  }
  Eigen::RowVectorXd wx = w.transpose() * X;
  Eigen::RowVectorXd gwx = wx.array().tanh();
  Eigen::RowVectorXd g_wx_deriv = 1.0 - gwx.array().square();
  Eigen::VectorXd updated_w = (X * gwx.transpose()) / n - g_wx_deriv.mean() * w;
  Eigen::VectorXd new_w = w + learningRate * (updated_w - w);
  new_w.normalize();
  return new_w;
}

// Unmix the observed signals X using the unmixing matrix W.
Eigen::MatrixXd unmixSignals(const Eigen::MatrixXd &W, const Eigen::MatrixXd &X) {
  return W * X;
}

// Orthogonalize the weight matrix using the Gram-Schmidt process.
Eigen::MatrixXd orthogonalizeWeightMatrix(const Eigen::MatrixXd &W) {
  Eigen::MatrixXd orthW = W;
  int n = orthW.cols();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      double proj = orthW.col(i).dot(orthW.col(j));
      orthW.col(i) -= proj * orthW.col(j);
    }
    orthW.col(i).normalize();
  }
  return orthW;
}

// Compute the relative change between two unmixing matrices W_old and W_new
double computeRelativeChange(const Eigen::MatrixXd &W_old, const Eigen::MatrixXd &W_new) {
  double normDiff = (W_new - W_old).norm();
  double normW_old = W_old.norm();
  if (normW_old < 1e-9) {
      return 0.0;
  }
  return normDiff / normW_old;
}

// Update the learning rate based on the current error relative to the previous error.
double updateLearningRate(double currentError, double previousError, double currentRate) {
  if (currentError > previousError) {
    return currentRate * 0.85;
  }
  return currentRate * 1.10;
}

int main() {
  const int numSources = 3;
  const int numSamples = 1000;

  Eigen::MatrixXd X = Eigen::MatrixXd::Random(numSources, numSamples);

  Eigen::MatrixXd W = Eigen::MatrixXd::Identity(numSources, numSources);

  double learningRate = 0.1;

  double relativeChange = 1.0;
  const double tolerance = 1e-6;
  const int maxIterations = 1000;
  int iteration = 0;
  double previousError = 1e10;

  while (relativeChange > tolerance && iteration < maxIterations) {
    Eigen::MatrixXd W_old = W;

    for (int i = 0; i < numSources; ++i) {
      Eigen::VectorXd currentW = W.row(i).transpose();
      Eigen::VectorXd newW = adaptiveFastICAUpdate(currentW, X, learningRate);
      W.row(i) = newW.transpose();
    }

    W = orthogonalizeWeightMatrix(W);

    relativeChange = computeRelativeChange(W_old, W);

    learningRate = updateLearningRate(relativeChange, previousError, learningRate);
    previousError = relativeChange;
    iteration++;
  }

  std::cout << "ICA converged in " << iteration << " iterations." << std::endl;
  std::cout << "Estimated Unmixing Matrix W:" << std::endl << W << std::endl;

  Eigen::MatrixXd S = unmixSignals(W, X);
  std::cout << "Recovered Source Signals (first 5 samples):" << std::endl;
  std::cout << S.leftCols(5) << "..." << std::endl;

  std::vector<double> signal(numSamples);
  for (int i = 0; i < numSamples; ++i) {
    signal[i] = S(0, i);
  }
  double kurtosisValue = computeKurtosis(signal);
  std::cout << "Kurtosis of the first source signal: " << kurtosisValue << std::endl;

  return 0;
}
