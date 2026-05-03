#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Function to compute the Gaussian likelihood of an observation
double calculateLikelihood(double observed, double predicted, double variance) {
  if (variance <= 0) {
      throw std::invalid_argument("Variance must be positive.");
  }
  double diff = observed - predicted;
  double exponent = -0.5 * (diff * diff) / variance;
  double likelihood = 1.0 / std::sqrt(2.0 * M_PI * variance) * std::exp(exponent);
  return likelihood;
}

// Function implementing systematic resampling.
std::vector<int> systematicResample(const std::vector<double>& weights) {
  int N = weights.size();
  if (N == 0) {
      return {};
  }
  std::vector<int> indices;
  indices.reserve(N);
  double step = 1.0 / N;
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0.0, step);
  double r = dist(gen);

  double c = weights[0];
  int i = 0;
  for (int m = 0; m < N; ++m) {
    double U = r + m * step;
    while (U > c && i < N - 1) {
      i++;
      c += weights[i];
    }
    indices.push_back(i);
  }
  return indices;
}

// Function to predict the next state in a linear state-space model.
Eigen::VectorXd predictState(const Eigen::VectorXd &currentState,
                           const Eigen::MatrixXd &transitionMatrix,
                           const Eigen::VectorXd &controlInput,
                           const Eigen::VectorXd &processNoise) {
  return transitionMatrix * currentState + controlInput + processNoise;
}

// Function to compute the Effective Sample Size (ESS)
double computeEffectiveSampleSize(const std::vector<double> &normalizedWeights) {
  if (normalizedWeights.empty()) {
      return 0.0;
  }
  double sumSq = 0.0;
  for (double w : normalizedWeights) {
    sumSq += w * w;
  }
  if (sumSq < 1e-9) {
      return 0.0;
  }
  return 1.0 / sumSq;
}

// Function implementing residual resampling.
std::vector<int> residualResample(const std::vector<double>& normalizedWeights) {
  int N = normalizedWeights.size();
  if (N == 0) {
      return {};
  }
  std::vector<int> indices;
  indices.reserve(N);

  std::vector<int> deterministicCount(N, 0);
  std::vector<double> residualWeights;
  residualWeights.reserve(N);

  int sumCount = 0;
  for (int i = 0; i < N; ++i) {
    deterministicCount[i] = static_cast<int>(N * normalizedWeights[i]);
    sumCount += deterministicCount[i];
    residualWeights.push_back(N * normalizedWeights[i] - deterministicCount[i]);
  }

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < deterministicCount[i]; j++) {
      indices.push_back(i);
    }
  }

  double totalResidual = 0.0;
  for (double w : residualWeights) {
    totalResidual += w;
  }
  if (totalResidual > 1e-9) {
      for (int i = 0; i < N; ++i) {
        residualWeights[i] /= totalResidual;
      }
  }

  int remaining = N - sumCount;
  if (remaining > 0) {
      double step = 1.0 / remaining;
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<> dist(0.0, step);
      double r = dist(gen);
      double c = residualWeights[0];
      int index = 0;
      for (int m = 0; m < remaining; ++m) {
        double U = r + m * step;
        while (U > c && index < N - 1) {
          index++;
          c += residualWeights[index];
        }
        indices.push_back(index);
      }
  }
  
  return indices;
}

// Function to generate an adaptive proposal state for importance sampling.
Eigen::VectorXd generateAdaptiveProposal(const Eigen::VectorXd &currentState,
                                       const Eigen::MatrixXd &proposalCovariance) {
  Eigen::LLT<Eigen::MatrixXd> lltOfCov(proposalCovariance);
  if(lltOfCov.info() == Eigen::NumericalIssue)
  {
    throw std::runtime_error("Possibly non semi-positive definite matrix!");
  }
  Eigen::MatrixXd L = lltOfCov.matrixL();

  Eigen::VectorXd standardGaussian(currentState.size());
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> dist(0.0, 1.0);
  for (int i = 0; i < currentState.size(); ++i) {
    standardGaussian(i) = dist(gen);
  }

  Eigen::VectorXd proposalNoise = L * standardGaussian;
  return currentState + proposalNoise;
}

// Function to propagate a set of particles concurrently using OpenMP.
void propagateParticlesConcurrently(std::vector<Eigen::VectorXd> &particles,
                                  const Eigen::MatrixXd &transitionMatrix,
                                  const Eigen::VectorXd &controlInput) {
  #pragma omp parallel for
  for (size_t i = 0; i < particles.size(); ++i) {
    particles[i] = transitionMatrix * particles[i] + controlInput;
  }
}

int main() {
  double observed = 1.0;
  double predicted = 0.8;
  double variance = 0.1;
  double likelihood = calculateLikelihood(observed, predicted, variance);
  std::cout << "Likelihood: " << likelihood << std::endl;

  std::vector<double> weights = {0.2, 0.3, 0.1, 0.4};

  std::vector<int> systematicIndices = systematicResample(weights);
  std::cout << "Systematic Resampling Indices: ";
  for (int idx : systematicIndices) {
    std::cout << idx << " ";
  }
  std::cout << std::endl;

  double ESS = computeEffectiveSampleSize(weights);
  std::cout << "Effective Sample Size: " << ESS << std::endl;

  std::vector<int> residualIndices = residualResample(weights);
  std::cout << "Residual Resampling Indices: ";
  for (int idx : residualIndices) {
    std::cout << idx << " ";
  }
  std::cout << std::endl;

  Eigen::VectorXd currentState(2);
  currentState << 1.0, 2.0;
  Eigen::MatrixXd transitionMatrix(2, 2);
  transitionMatrix << 0.9, 0.1, 0.2, 0.8;
  Eigen::VectorXd controlInput(2);
  controlInput << 0.5, -0.3;
  Eigen::VectorXd processNoise(2);
  processNoise << 0.05, -0.05;

  Eigen::VectorXd nextState = predictState(currentState, transitionMatrix, controlInput, processNoise);
  std::cout << "Predicted Next State:\n" << nextState << std::endl;

  Eigen::MatrixXd proposalCovariance = Eigen::MatrixXd::Identity(2, 2) * 0.1;
  Eigen::VectorXd adaptiveProposal = generateAdaptiveProposal(currentState, proposalCovariance);
  std::cout << "Adaptive Proposal State:\n" << adaptiveProposal << std::endl;

  std::vector<Eigen::VectorXd> particles;
  for (int i = 0; i < 4; ++i) {
    Eigen::VectorXd particle = currentState;
    particle(0) += 0.1 * i;
    particle(1) += 0.1 * i;
    particles.push_back(particle);
  }

  propagateParticlesConcurrently(particles, transitionMatrix, controlInput);
  std::cout << "Particles after concurrent propagation:" << std::endl;
  for (const auto &p : particles) {
    std::cout << p.transpose() << std::endl;
  }

  return 0;
}
