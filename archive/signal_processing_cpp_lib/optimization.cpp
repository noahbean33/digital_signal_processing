#include <iostream>
#include <vector>
#include <cmath>
#include <future>
#include <numeric>
#include <algorithm>
#include <Eigen/Dense>

double updateParameter(double current, double gradient, double learningRate) {
  return current - learningRate * gradient;
}

int selectModel(const std::vector<double>& costs) {
  if (costs.empty()) {
    throw std::invalid_argument("Costs vector is empty");
  }
  auto minIt = std::min_element(costs.begin(), costs.end());
  return std::distance(costs.begin(), minIt);
}

double computeLeastSquaresError(const std::vector<double>& predictions, const std::vector<double>& observations) {
  if (predictions.size() != observations.size() || predictions.empty()) {
    throw std::invalid_argument("Predictions and observations must have the same size and not be empty");
  }
  double sumError = 0.0;
  for (size_t i = 0; i < predictions.size(); ++i) {
    double error = predictions[i] - observations[i];
    sumError += error * error;
  }
  if (predictions.size() == 0) {
    throw std::runtime_error("Division by zero in computeLeastSquaresError");
  }
  return sumError / static_cast<double>(predictions.size());
}

double simulateNonlinearSystem(double state, double input) {
  return state + input - 0.2 * std::sin(state);
}

double computeMeanAbsoluteError(const std::vector<double>& predictions, const std::vector<double>& observations) {
  if (predictions.size() != observations.size() || predictions.empty()) {
    throw std::invalid_argument("Predictions and observations must have the same size and not be empty");
  }
  double totalError = 0.0;
  for (size_t i = 0; i < predictions.size(); ++i) {
    totalError += std::abs(predictions[i] - observations[i]);
  }
  if (predictions.size() == 0) {
    throw std::runtime_error("Division by zero in computeMeanAbsoluteError");
  }
  return totalError / static_cast<double>(predictions.size());
}

Eigen::VectorXd computeTrustRegionStep(const Eigen::MatrixXd& H, const Eigen::VectorXd& grad, double delta) {
  if (H.rows() != H.cols() || grad.size() != H.rows()) {
    throw std::invalid_argument("H must be a square matrix and grad must have the same size as H's rows");
  }
  double gradDotHGrad = grad.dot(H * grad);
  if (std::abs(gradDotHGrad) < 1e-9) { // Avoid division by zero
    return Eigen::VectorXd::Zero(grad.size());
  }
  Eigen::VectorXd p = -(grad.dot(grad) / gradDotHGrad) * grad;
  if (p.norm() > delta) {
    p = delta * p.normalized();
  }
  return p;
}

double computeL2RegularizedError(const std::vector<double>& predictions, const std::vector<double>& observations, const std::vector<double>& parameters, double lambda) {
  double mse = computeLeastSquaresError(predictions, observations);
  double regSum = 0.0;
  for (double param : parameters) {
    regSum += param * param;
  }
  return mse + lambda * regSum;
}

double parallelComputeMSE(const std::vector<double>& predictions, const std::vector<double>& observations) {
  if (predictions.size() != observations.size() || predictions.empty()) {
    throw std::invalid_argument("Predictions and observations must have the same size and not be empty");
  }
  size_t mid = predictions.size() / 2;
  auto computePartialError = [&predictions, &observations](size_t start, size_t end) -> double {
    double sumError = 0.0;
    for (size_t i = start; i < end; ++i) {
      double error = predictions[i] - observations[i];
      sumError += error * error;
    }
    return sumError;
  };

  std::future<double> futureError = std::async(std::launch::async, computePartialError, 0, mid);
  double secondHalfError = computePartialError(mid, predictions.size());
  double firstHalfError = futureError.get();
  double totalError = firstHalfError + secondHalfError;
  if (predictions.size() == 0) {
    throw std::runtime_error("Division by zero in parallelComputeMSE");
  }
  return totalError / static_cast<double>(predictions.size());
}

int main() {
  double param = 1.0;
  double grad = 0.1;
  double learningRate = 0.05;
  double updatedParam = updateParameter(param, grad, learningRate);
  std::cout << "Updated parameter: " << updatedParam << std::endl;

  std::vector<double> modelCosts = {10.5, 8.3, 12.7, 7.9 };
  int bestModelIndex = selectModel(modelCosts);
  std::cout << "Best model index: " << bestModelIndex << std::endl;

  std::vector<double> predictions = {1.2, 2.3, 3.7, 4.1 };
  std::vector<double> observations = {1.0, 2.5, 3.5, 4.0 };

  double mse = computeLeastSquaresError(predictions, observations);
  std::cout << "Mean Squared Error: " << mse << std::endl;

  double mae = computeMeanAbsoluteError(predictions, observations);
  std::cout << "Mean Absolute Error: " << mae << std::endl;

  double state = 0.5, input = 1.0;
  double nextState = simulateNonlinearSystem(state, input);
  std::cout << "Next state (simulated nonlinear system): " << nextState << std::endl;

  Eigen::MatrixXd H(2, 2);
  H << 4, 1, 1, 3;
  Eigen::VectorXd gradient(2);
  gradient << 2, -1;
  double delta = 0.5;
  Eigen::VectorXd trustStep = computeTrustRegionStep(H, gradient, delta);
  std::cout << "Trust region step: [" << trustStep[0] << ", " << trustStep[1] << "]" << std::endl;

  std::vector<double> parameters = { 0.5, -0.3, 0.8 };
  double lambda = 0.1;
  double regError = computeL2RegularizedError(predictions, observations, parameters, lambda);
  std::cout << "L2 Regularized Error: " << regError << std::endl;

  double parallelMSE = parallelComputeMSE(predictions, observations);
  std::cout << "Parallel Computed MSE: " << parallelMSE << std::endl;

  return 0;
}
