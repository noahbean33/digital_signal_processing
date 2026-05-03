#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Forward declarations
void adaptiveFilterProcess(const std::vector<double>& inputSignal,
                           const std::vector<double>& desiredSignal,
                           std::vector<double>& weights, double initialStepSize,
                           std::vector<double>& outputSignal,
                           std::vector<double>& errorSignal);
double computeMeanSquaredError(const std::vector<double>& errors);

// Function to update weights using the LMS algorithm.
void updateLMS(std::vector<double>& weights, double stepSize, double error,
             const std::vector<double>& input) {
    for (std::size_t i = 0; i < weights.size(); ++i) {
        weights[i] += stepSize * error * input[i];
    }
}

// Function to compute the error signal: e(n) = d(n) - y(n)
double computeError(double desired, double output) {
    return desired - output;
}

// Dynamic adjustment of step-size
double updateStepSize(double currentError, double mu, double beta) {
    return mu / (1 + beta * std::fabs(currentError));
}

// Update weights with an L2 regularization term.
void updateWeightsRegularized(std::vector<double>& weights, double stepSize,
                            double error, const std::vector<double>& inputSegment,
                            double lambda) {
    for (std::size_t i = 0; i < weights.size(); ++i) {
        weights[i] += stepSize * (error * inputSegment[i] - lambda * weights[i]);
    }
}

// Compute the dot product between two vectors.
double dotProduct(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    if (vec1.size() != vec2.size()) {
        return 0.0; // Vectors must be the same size
    }
    double product = 0.0;
    for (std::size_t i = 0; i < vec1.size(); ++i) {
        product += vec1[i] * vec2[i];
    }
    return product;
}

// Compute the Euclidean (L2) norm difference between two weight vectors.
double computeWeightDifference(const std::vector<double>& prevWeights,
                               const std::vector<double>& currWeights) {
    double normDiff = 0.0;
    for (std::size_t i = 0; i < prevWeights.size(); ++i) {
        double diff = currWeights[i] - prevWeights[i];
        normDiff += diff * diff;
    }
    return std::sqrt(normDiff);
}

int main() {
    std::srand(static_cast<unsigned int>(std::time(0)));

    const std::size_t signalLength = 100;
    const std::size_t filterOrder = 5;
    double initialStepSize = 0.05;

    std::vector<double> inputSignal(signalLength, 0.0);
    std::vector<double> desiredSignal(signalLength, 0.0);

    double frequency = 2 * M_PI / 20.0;
    for (std::size_t n = 0; n < signalLength; ++n) {
        double noise = ((std::rand() % 100) / 500.0 - 0.1);
        inputSignal[n] = std::sin(frequency * n) + noise;
        desiredSignal[n] = std::sin(frequency * n);
    }

    std::vector<double> weights(filterOrder, 0.1);
    std::vector<double> outputSignal;
    std::vector<double> errorSignal;

    std::vector<double> initialWeights = weights;

    adaptiveFilterProcess(inputSignal, desiredSignal, weights, initialStepSize, outputSignal, errorSignal);

    double mse = computeMeanSquaredError(errorSignal);

    std::cout << "Final Weights: ";
    for (double w : weights) {
        std::cout << w << " ";
    }
    std::cout << std::endl;

    std::cout << "Mean Squared Error: " << mse << std::endl;

    double weightDiff = computeWeightDifference(initialWeights, weights);
    std::cout << "Weight Difference (Convergence Measure): " << weightDiff << std::endl;

    return 0;
}

void adaptiveFilterProcess(const std::vector<double>& inputSignal,
                           const std::vector<double>& desiredSignal,
                           std::vector<double>& weights, double initialStepSize,
                           std::vector<double>& outputSignal,
                           std::vector<double>& errorSignal) {
    std::size_t filterOrder = weights.size();
    std::size_t signalLength = inputSignal.size();

    outputSignal.resize(signalLength, 0.0);
    errorSignal.resize(signalLength, 0.0);

    for (std::size_t n = filterOrder - 1; n < signalLength; ++n) {
        std::vector<double> inputSegment(filterOrder, 0.0);
        for (std::size_t k = 0; k < filterOrder; ++k) {
            inputSegment[k] = inputSignal[n - k];
        }

        double y = dotProduct(weights, inputSegment);
        outputSignal[n] = y;

        double error = computeError(desiredSignal[n], y);
        errorSignal[n] = error;

        const double BETA = 0.1;
double stepSize = updateStepSize(error, initialStepSize, BETA);

        const double LAMBDA = 0.01;
updateWeightsRegularized(weights, stepSize, error, inputSegment, LAMBDA);
    }
}

double computeMeanSquaredError(const std::vector<double>& errors) {
    if (errors.empty()) {
        return 0.0; // Avoid division by zero
    }
    double sum = 0.0;
    for (double err : errors) {
        sum += err * err;
    }
    return sum / errors.size();
}
