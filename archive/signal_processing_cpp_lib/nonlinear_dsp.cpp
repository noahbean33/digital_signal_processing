#include <iostream>
#include <std::vector>
#include <cmath>
#include <std::functional>
#include <algorithm>
#include <std::execution>
#include <numeric>
#include <cstdlib>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


std::vector<double> applyHarmonicDistortion(const std::vector<double>& input, double threshold, double gain) {
  if (std::abs(threshold) < 1e-9) {
      return input; // or throw an exception
  }
  std::vector<double> output;
  output.reserve(input.size());
  for (double sample : input) {
    double distorted = sample;
    if (std::fabs(sample) > threshold) {
      distorted = threshold * std::tanh(gain * (sample / threshold));
    }
    output.push_back(distorted);
  }
  return output;
}

std::vector<double> compressDynamicRange(const std::vector<double>& input, double compressionRatio) {
  if (compressionRatio <= 0) {
      throw std::invalid_argument("Compression ratio must be positive.");
  }
  double logCompressionRatio = std::log1p(compressionRatio);
  if (std::abs(logCompressionRatio) < 1e-9) {
      return input; // or throw an exception
  }
  std::vector<double> output;
  output.reserve(input.size());
  for (double sample : input) {
    double absSample = std::fabs(sample);
    double compressed = std::log1p(compressionRatio * absSample) / logCompressionRatio;
    output.push_back(sample >= 0.0 ? compressed : -compressed);
  }
  return output;
}

std::vector<double> nonlinearFIRFilter(const std::vector<double>& input, const std::vector<double>& coeffs, double steepness) {
  std::vector<double> output(input.size(), 0.0);
  size_t filterSize = coeffs.size();
  for (size_t n = 0; n < input.size(); ++n) {
    double sum = 0.0;
    for (size_t k = 0; k < filterSize; ++k) {
      if (n >= k) {
        sum += coeffs[k] * input[n - k];
      }
    }
    output[n] = 1.0 / (1.0 + std::exp(-steepness * sum));
  }
  return output;
}

double newtonRaphsonInvert(double target, std::function<double(double)> f, std::function<double(double)> df, double initialGuess, double tolerance = 1e-6, int maxIterations = 100) {
  double x = initialGuess;
  for (int i = 0; i < maxIterations; ++i) {
    double fx = f(x) - target;
    if (std::fabs(fx) < tolerance) {
      return x;
    }
    double dfx = df(x);
    if (std::fabs(dfx) < 1e-12) {
      break;
    }
    x -= fx / dfx;
  }
  return x;
}

std::vector<double> adaptiveNonlinearUpdate(const std::vector<double>& weights, const std::vector<double>& input, double desired, double learningRate, double steepness) {
  if (weights.size() != input.size()) {
      throw std::invalid_argument("Weights and input std::vectors must have the same size.");
  }
  double sum = 0.0;
  for (size_t i = 0; i < weights.size(); ++i) {
    sum += weights[i] * input[i];
  }
  double output = 1.0 / (1.0 + std::exp(-steepness * sum));
  double error = desired - output;
  double gradient = steepness * output * (1.0 - output);

  std::vector<double> updatedWeights(weights.size());
  for (size_t i = 0; i < weights.size(); ++i) {
    updatedWeights[i] = weights[i] + learningRate * error * gradient * input[i];
  }
  return updatedWeights;
}

double monteCarloEstimate(std::function<double()> simulationFunc, int iterations) {
  if (iterations <= 0) {
      return 0.0;
  }
  double sum = 0.0;
  for (int i = 0; i < iterations; ++i) {
    sum += simulationFunc();
  }
  return sum / iterations;
}

double dampedNewtonRaphsonSolve(std::function<double(double)> f, std::function<double(double)> df, double initial, double damping, double tol = 1e-6, int maxIter = 100) {
  double x = initial;
  for (int i = 0; i < maxIter; ++i) {
    double val = f(x);
    if (std::fabs(val) < tol) break;
    double derivative = df(x);
    if (std::fabs(derivative) < 1e-12) break;
    x -= damping * (val / derivative);
  }
  return x;
}

std::vector<double> applyNonlinearTransformationParallel(const std::vector<double>& input, std::function<double(double)> nonlinearFunc) {
  std::vector<double> output(input.size());
  std::transform(std::execution::par, input.begin(), input.end(), output.begin(), nonlinearFunc);
  return output;
}

int main() {
  std::vector<double> signal = {-0.8, -0.5, 0.0, 0.5, 0.8, 1.0, -1.0};

  double threshold = 0.7;
  double gain = 2.0;
  std::vector<double> distortedSignal = applyHarmonicDistortion(signal, threshold, gain);
  std::cout << "Harmonic Distortion: " << std::endl;
  for (double val : distortedSignal) { std::cout << val << " "; }
  std::cout << std::endl;

  double compressionRatio = 2.0;
  std::vector<double> compressedSignal = compressDynamicRange(signal, compressionRatio);
  std::cout << "Dynamic Range Compression:" << std::endl;
  for (double val : compressedSignal) { std::cout << val << " "; }
  std::cout << std::endl;

  std::vector<double> coeffs = {0.2, 0.3, 0.5};
  double steepness = 1.0;
  std::vector<double> firOutput = nonlinearFIRFilter(signal, coeffs, steepness);
  std::cout << "Nonlinear FIR Filter Output:" << std::endl;
  for (double val : firOutput) { std::cout << val << " "; }
  std::cout << std::endl;

  auto f = [](double x) { return std::exp(x); };
  auto df = [](double x) { return std::exp(x); };
  double target = 3.0;
  double initialGuess = 1.0;
  const double TOLERANCE = 1e-6;
const int MAX_ITERATIONS = 100;
double inverted = newtonRaphsonInvert(target, f, df, initialGuess, TOLERANCE, MAX_ITERATIONS);
  std::cout << "Newton-Raphson Inversion (std::exp(x)=3) gives x = " << inverted << std::endl;

  std::vector<double> weights = {0.5, -0.3, 0.8};
  std::vector<double> input = {1.0, 0.5, -1.0};
  double desired = 0.7;
  double learningRate = 0.1;
  std::vector<double> updatedWeights = adaptiveNonlinearUpdate(weights, input, desired, learningRate, steepness);
  std::cout << "Adaptive Nonlinear Update Weights:" << std::endl;
  for (double w : updatedWeights) { std::cout << w << " "; }
  std::cout << std::endl;

  auto simulationFunc = []() -> double {
    double x = (static_cast<double>(rand()) / RAND_MAX) * M_PI;
    return sin(x);
  };
  double estimatedExpectation = monteCarloEstimate(simulationFunc, 100000);
  std::cout << "Monte Carlo Estimated Expectation of sin(x) over [0,pi]: " << estimatedExpectation << std::endl;

  auto f_cos = [](double x) { return cos(x); };
  auto df_cos = [](double x) { return -sin(x); };
  double initial = 1.5;
  double damping = 0.9;
  const double DAMPING_TOLERANCE = 1e-6;
const int DAMPING_MAX_ITERATIONS = 100;
double solvedX = dampedNewtonRaphsonSolve(f_cos, df_cos, initial, damping, DAMPING_TOLERANCE, DAMPING_MAX_ITERATIONS);
  std::cout << "Damped Newton-Raphson Solve for cos(x)=0 gives x = " << solvedX << std::endl;

  auto squareFunc = [](double x) { return x * x; };
  std::vector<double> parallelOutput = applyNonlinearTransformationParallel(signal, squareFunc);
  std::cout << "Parallel Nonlinear Transformation (Square):" << std::endl;
  for (double val : parallelOutput) { std::cout << val << " "; }
  std::cout << std::endl;

  return 0;
}
