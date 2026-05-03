#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <Eigen/Dense>

// Normalized Least Mean Squares (LMS) Update Function
std::vector<double> normalizedLMSUpdate(const std::vector<double>& weights,
                                        const std::vector<double>& input,
                                        double error,
                                        double mu,
                                        double epsilon = 1e-6) {
    if (weights.size() != input.size() || weights.empty()) {
        throw std::invalid_argument("Weights and input vectors must have the same size and not be empty.");
    }
    std::vector<double> updatedWeights(weights.size());
    double inputEnergy = 0.0;
    for (auto x : input) {
        inputEnergy += x * x;
    }
    inputEnergy += epsilon; // Prevent division by zero in low energy conditions
    for (size_t i = 0; i < weights.size(); ++i) {
        updatedWeights[i] = weights[i] + mu * error * input[i] / inputEnergy;
    }
    return updatedWeights;
}

// Recursive Least Squares (RLS) Gain Computation using Eigen
Eigen::VectorXd computeRLSGain(const Eigen::MatrixXd& P,
                               const Eigen::VectorXd& x,
                               double lambda) {
    Eigen::VectorXd Px = P * x;
    double denominator = lambda + x.dot(Px);
    if (std::abs(denominator) < 1e-9) {
        return Eigen::VectorXd::Zero(x.size());
    }
    return Px / denominator;
}

// Update Inverse Correlation Matrix for RLS Algorithm
Eigen::MatrixXd updateInverseCorrelation(const Eigen::MatrixXd& P,
                                         const Eigen::VectorXd& x,
                                         double lambda) {
    Eigen::VectorXd Px = P * x;
    double denominator = lambda + x.dot(Px);
    if (std::abs(denominator) < 1e-9) {
        return P;
    }
    // Using the matrix inversion lemma for efficient update
    return (P - (Px * Px.transpose()) / denominator) / lambda;
}

// Clamped Weight Update Function to restrict weight change magnitude
double clampedWeightUpdate(double currentWeight,
                         double updateValue,
                         double maxUpdate) {
    double appliedUpdate = std::clamp(updateValue, -maxUpdate, maxUpdate);
    return currentWeight + appliedUpdate;
}

// Adaptive Step Size Adjustment Function based on instantaneous error
double computeAdaptiveStepSize(double errorMagnitude,
                             double previousStepSize,
                             double adjustmentFactor) {
    double newStep = (errorMagnitude > 1.0)
                     ? previousStepSize * (1.0 - adjustmentFactor)
                     : previousStepSize * (1.0 + adjustmentFactor);
    return newStep;
}

// Compute Convergence Metric from error history
double computeConvergenceMetric(const std::vector<double>& errorHistory) {
    if (errorHistory.size() < 2) return 0.0;
    double sumAbsDiff = 0.0;
    for (size_t i = 1; i < errorHistory.size(); ++i) {
        sumAbsDiff += std::abs(errorHistory[i] - errorHistory[i - 1]);
    }
    return sumAbsDiff / static_cast<double>(errorHistory.size() - 1);
}

// Clamp a given value within specified minimum and maximum bounds
double clampValue(double value, double minVal, double maxVal) {
    return std::clamp(value, minVal, maxVal);
}

int main() {
    // Example 1: Normalized LMS Update
    std::vector<double> weights = {0.0, 0.0, 0.0};
    std::vector<double> input = {1.2, -0.7, 0.5};
    double error = 0.3;
    double mu = 0.05;
    double epsilon = 1e-6;
    std::vector<double> updatedWeights =
        normalizedLMSUpdate(weights, input, error, mu, epsilon);
    std::cout << "Normalized LMS Updated Weights:" << std::endl;
    for (double w : updatedWeights) {
        std::cout << w << " ";
    }
    std::cout << std::endl << std::endl;

    // Example 2: RLS Algorithm Components
    int filterOrder = 3;
    Eigen::VectorXd rlsWeights = Eigen::VectorXd::Zero(filterOrder);
    Eigen::MatrixXd P = Eigen::MatrixXd::Identity(filterOrder, filterOrder) * 1000; // Initial large inverse correlation matrix
    Eigen::VectorXd x(filterOrder);
    x << 1.2, -0.7, 0.5;
    double lambda = 0.99;

    // Compute RLS gain vector k
    Eigen::VectorXd k = computeRLSGain(P, x, lambda);
    std::cout << "RLS Gain Vector:" << std::endl << k << std::endl << std::endl;

    // Update RLS filter coefficients: w(n) = w(n-1) + k * e(n)
    rlsWeights = rlsWeights + k * error;
    std::cout << "RLS Updated Filter Coefficients:" << std::endl << rlsWeights << std::endl << std::endl;

    // Update the inverse correlation matrix P using the new data
    P = updateInverseCorrelation(P, x, lambda);
    std::cout << "Updated Inverse Correlation Matrix P:" << std::endl << P << std::endl << std::endl;

    // Example 3: Adaptive Step Size adjustment based on error magnitude
    double previousStepSize = 0.05;
    double adjustmentFactor = 0.01;
    double newStepSize = computeAdaptiveStepSize(std::abs(error), previousStepSize, adjustmentFactor);
    std::cout << "Adaptive Step Size: " << newStepSize << std::endl << std::endl;

    // Example 4: Clamped Weight Update to limit sudden changes
    double currentWeight = 0.10;
    double proposedUpdate = 0.08;
    double maxAllowedUpdate = 0.05;
    double clampedWeight = clampedWeightUpdate(currentWeight, proposedUpdate, maxAllowedUpdate);
    std::cout << "Clamped Weight Updated Value: " << clampedWeight << std::endl << std::endl;

    // Example 5: Convergence Metric Calculation from error history
    std::vector<double> errorHistory = {0.5, 0.4, 0.35, 0.32, 0.31};
    double convergenceMetric =
        computeConvergenceMetric(errorHistory);
    std::cout << "Convergence Metric: " << convergenceMetric << std::endl << std::endl;

    // Example 6: Demonstration of value clamping to prevent undesirable update magnitudes
    double value = 1.2;
    double clampedVal = clampValue(value, 0.0, 1.0);
    std::cout << "Clamped Value (1.2 clamped to [0, 1]): " << clampedVal << std::endl;

    return 0;
}
