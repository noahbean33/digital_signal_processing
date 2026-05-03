#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <numeric>

// Function to perform one-dimensional valid convolution on a signal
std::vector<double> applyConvolution(const std::vector<double>& signal, const std::vector<double>& kernel) {
  if (signal.empty() || kernel.empty() || kernel.size() > signal.size()) {
      return {};
  }
  int signalSize = signal.size();
  int kernelSize = kernel.size();
  int outputSize = signalSize - kernelSize + 1;
  std::vector<double> output(outputSize, 0.0);
  for (int i = 0; i < outputSize; ++i) {
    for (int j = 0; j < kernelSize; ++j) {
      output[i] += signal[i + j] * kernel[j];
    }
  }
  return output;
}

// Function to compute the softmax activation for a given vector of logits.
std::vector<double> softmax(const std::vector<double>& logits) {
  if (logits.empty()) {
      return {};
  }
  std::vector<double> expValues;
  expValues.reserve(logits.size());
  double maxLogit = *std::max_element(logits.begin(), logits.end());
  double sumExp = 0.0;
  for (double value : logits) {
    double expVal = std::exp(value - maxLogit);
    expValues.push_back(expVal);
    sumExp += expVal;
  }
  
  if (sumExp < 1e-9) {
      return std::vector<double>(logits.size(), 1.0 / logits.size());
  }

  std::vector<double> probabilities;
  probabilities.reserve(logits.size());
  for (double expVal : expValues) {
    probabilities.push_back(expVal / sumExp);
  }
  return probabilities;
}

// Function to perform Xavier Initialization for a layer.
std::vector<double> xavierInitialization(int fanIn, int fanOut) {
  if (fanIn <= 0 || fanOut <= 0) {
      return {};
  }
  std::vector<double> weights(fanIn * fanOut);
  double limit = std::sqrt(6.0 / (fanIn + fanOut));
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(-limit, limit);
  for (auto &weight : weights) {
    weight = dis(gen);
  }
  return weights;
}

// Function to apply dropout to a set of activations.
std::vector<double> applyDropout(const std::vector<double>& activations, double dropoutProb) {
  if (dropoutProb < 0.0 || dropoutProb > 1.0) {
      throw std::invalid_argument("Dropout probability must be between 0.0 and 1.0.");
  }
  std::vector<double> output(activations.size());
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0.0, 1.0);
  for (std::size_t i = 0; i < activations.size(); ++i) {
    output[i] = (dis(gen) > dropoutProb) ? activations[i] : 0.0;
  }
  return output;
}

// Function to clip a gradient value to ensure it does not exceed a predefined threshold.
double clipGradient(double gradient, double threshold) {
  return std::max(-threshold, std::min(threshold, gradient));
}

// Function to apply Batch Normalization to the input activations.
std::vector<double> applyBatchNormalization(const std::vector<double>& inputs,
                                          double mean,
                                          double variance,
                                          double gamma,
                                          double beta,
                                          double epsilon) {
  if (inputs.empty()) {
      return {};
  }
  std::vector<double> outputs;
  outputs.reserve(inputs.size());
  double denom = std::sqrt(variance + epsilon);
  if (denom < 1e-9) {
      return inputs; // Avoid division by zero
  }
  for (double x : inputs) {
    outputs.push_back(gamma * ((x - mean) / denom) + beta);
  }
  return outputs;
}

int main() {
  // Example 1: Convolution
  std::vector<double> signal = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> kernel = {0.2, 0.5, 0.3};
  std::vector<double> convOutput = applyConvolution(signal, kernel);
  
  std::cout << "Convolution Output: ";
  for (double val : convOutput) {
    std::cout << val << " ";
  }
  std::cout << std::endl;

  // Example 2: Softmax Activation
  std::vector<double> logits = {1.2, 0.9, 3.0};
  std::vector<double> probabilities = softmax(logits);
  std::cout << "Softmax Probabilities: ";
  for (double prob : probabilities) {
    std::cout << prob << " ";
  }
  std::cout << std::endl;

  // Example 3: Xavier Initialization
  const int FAN_IN = 3;
  const int FAN_OUT = 2;
  std::vector<double> weights = xavierInitialization(FAN_IN, FAN_OUT);
  std::cout << "Xavier Initialized Weights: ";
  for (double w : weights) {
    std::cout << w << " ";
  }
  std::cout << std::endl;

  // Example 4: Applying Dropout
  std::vector<double> activations = {0.5, 0.8, 0.3, 0.9};
  const double DROPOUT_PROB = 0.25;
  std::vector<double> droppedActivations = applyDropout(activations, DROPOUT_PROB);
  std::cout << "Activations after Dropout: ";
  for (double act : droppedActivations) {
    std::cout << act << " ";
  }
  std::cout << std::endl;

  // Example 5: Gradient Clipping
  const double GRADIENT = 5.5;
  const double THRESHOLD = 3.0;
  double clippedGradient = clipGradient(GRADIENT, THRESHOLD);
  std::cout << "Clipped Gradient: " << clippedGradient << std::endl;

  // Example 6: Batch Normalization
  std::vector<double> inputs = {1.0, 2.0, 3.0, 4.0};
  const double MEAN = 2.5;
  const double VARIANCE = 1.25;
  const double GAMMA = 1.0;
  const double BETA = 0.0;
  const double EPSILON = 1e-5;
  
  std::vector<double> normalizedOutputs = applyBatchNormalization(inputs, MEAN, VARIANCE, GAMMA, BETA, EPSILON);
  std::cout << "Batch Normalized Outputs: ";
  for (double normVal : normalizedOutputs) {
    std::cout << normVal << " ";
  }
  std::cout << std::endl;

  return 0;
}
