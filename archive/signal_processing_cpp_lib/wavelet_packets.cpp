#include <iostream>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Eigen/Dense>

// Single-level wavelet packet decomposition:
std::vector<Eigen::VectorXd> waveletPacketDecompose(const Eigen::VectorXd &signal,
                                                  const Eigen::VectorXd &lowPassFilter,
                                                  const Eigen::VectorXd &highPassFilter) {
  if (lowPassFilter.size() != highPassFilter.size()) {
      throw std::invalid_argument("Filters must have the same size.");
  }
  int filterSize = lowPassFilter.size();
  int length = signal.size();
  int numCoeffs = length / 2;
  Eigen::VectorXd approx = Eigen::VectorXd::Zero(numCoeffs);
  Eigen::VectorXd detail = Eigen::VectorXd::Zero(numCoeffs);

  for (int i = 0; i < numCoeffs; ++i) {
    double aVal = 0.0, dVal = 0.0;
    for (int j = 0; j < filterSize; ++j) {
      int idx = (2 * i + j) % length;
      aVal += signal(idx) * lowPassFilter(j);
      dVal += signal(idx) * highPassFilter(j);
    }
    approx(i) = aVal;
    detail(i) = dVal;
  }
  return {approx, detail};
}

// Compute the energy of the wavelet packet coefficients.
double computePacketEnergy(const Eigen::VectorXd &coefficients) {
  return coefficients.squaredNorm();
}

// Adaptive basis selection using recursive decomposition:
Eigen::VectorXd adaptiveBasisSelection(const Eigen::VectorXd &signal,
                                     const Eigen::VectorXd &lowPassFilter,
                                     const Eigen::VectorXd &highPassFilter,
                                     double energyThreshold,
                                     int maxLevel,
                                     int currentLevel = 0) {
  if (currentLevel == maxLevel) {
    return signal;
  }
  
  std::vector<Eigen::VectorXd> packets = waveletPacketDecompose(signal, lowPassFilter, highPassFilter);
  double energyLow = computePacketEnergy(packets[0]);
  double energyHigh = computePacketEnergy(packets[1]);

  if (energyLow > energyThreshold) {
    packets[0] = adaptiveBasisSelection(packets[0], lowPassFilter, highPassFilter, energyThreshold, maxLevel, currentLevel + 1);
  }
  
  if (energyHigh > energyThreshold) {
    packets[1] = adaptiveBasisSelection(packets[1], lowPassFilter, highPassFilter, energyThreshold, maxLevel, currentLevel + 1);
  }

  Eigen::VectorXd reassembled(signal.size());
  int halfSize = reassembled.size() / 2;
  reassembled.head(halfSize) = packets[0];
  reassembled.tail(halfSize) = packets[1];
  return reassembled;
}

// Recursive packet decomposition to obtain terminal (leaf) nodes of the wavelet packet tree.
std::vector<Eigen::VectorXd> recursivePacketDecomposition(const Eigen::VectorXd &signal,
                                                        const Eigen::VectorXd &lowPassFilter,
                                                        const Eigen::VectorXd &highPassFilter,
                                                        int maxLevel) {
  if (maxLevel == 0) {
    return {signal};
  }
  
  std::vector<Eigen::VectorXd> packets = waveletPacketDecompose(signal, lowPassFilter, highPassFilter);
  std::vector<Eigen::VectorXd> result;
  std::vector<Eigen::VectorXd> leftNodes = recursivePacketDecomposition(packets[0], lowPassFilter, highPassFilter, maxLevel - 1);
  std::vector<Eigen::VectorXd> rightNodes = recursivePacketDecomposition(packets[1], lowPassFilter, highPassFilter, maxLevel - 1);
  result.insert(result.end(), leftNodes.begin(), leftNodes.end());
  result.insert(result.end(), rightNodes.begin(), rightNodes.end());
  return result;
}

// Single-level reconstruction of the original signal from approximation and detail coefficients.
Eigen::VectorXd waveletPacketReconstruct(const Eigen::VectorXd &approx,
                                       const Eigen::VectorXd &detail,
                                       const Eigen::VectorXd &synthesisLowFilter,
                                       const Eigen::VectorXd &synthesisHighFilter) {
  if (synthesisLowFilter.size() != synthesisHighFilter.size()) {
      throw std::invalid_argument("Synthesis filters must have the same size.");
  }
  int numCoeffs = approx.size();
  int signalLength = numCoeffs * 2;
  Eigen::VectorXd signal = Eigen::VectorXd::Zero(signalLength);
  int filterSize = synthesisLowFilter.size();
  
  for (int i = 0; i < numCoeffs; ++i) {
    for (int j = 0; j < filterSize; ++j) {
      int idx = (2 * i + j) % signalLength;
      signal(idx) += approx(i) * synthesisLowFilter(j) + detail(i) * synthesisHighFilter(j);
    }
  }
  return signal;
}

// Compute the Root Mean Square Error (RMSE) between the original and reconstructed signals.
double computeReconstructionRMSE(const Eigen::VectorXd &original,
                               const Eigen::VectorXd &reconstructed) {
    if (original.size() != reconstructed.size() || original.size() == 0) {
        return 0.0;
    }
  double squaredError = (original - reconstructed).squaredNorm();
  return std::sqrt(squaredError / original.size());
}

// Parallel convolution using OpenMP to accelerate filtering operations.
Eigen::VectorXd parallelConvolve(const Eigen::VectorXd &signal,
                             const Eigen::VectorXd &filter) {
  int signalSize = signal.size();
  int filterSize = filter.size();
  Eigen::VectorXd output = Eigen::VectorXd::Zero(signalSize);
  #pragma omp parallel for
  for (int i = 0; i < signalSize; ++i) {
    double sum = 0.0;
    for (int j = 0; j < filterSize; ++j) {
      int idx = (i + j) % signalSize;
      sum += signal(idx) * filter(j);
    }
    output(i) = sum;
  }
  return output;
}

// Decision function to determine if further decomposition of a node is required based on its energy.
bool shouldDecomposeFurther(double nodeEnergy, double energyThreshold) {
  return nodeEnergy > energyThreshold;
}

int main() {
  Eigen::VectorXd signal(8);
  signal << 1, 2, 3, 4, 5, 6, 7, 8;

  Eigen::VectorXd lowPassFilter(2);
  lowPassFilter << 0.5, 0.5;
  Eigen::VectorXd highPassFilter(2);
  highPassFilter << -0.5, 0.5;

  Eigen::VectorXd synthesisLowFilter = lowPassFilter;
  Eigen::VectorXd synthesisHighFilter = highPassFilter;

  std::vector<Eigen::VectorXd> packets = waveletPacketDecompose(signal, lowPassFilter, highPassFilter);
    
  double energyLow = computePacketEnergy(packets[0]);
  double energyHigh = computePacketEnergy(packets[1]);
  std::cout << "Energy of Approximation: " << energyLow << std::endl;
  std::cout << "Energy of Detail: " << energyHigh << std::endl;
  
  double energyThreshold = 10.0;
  int maxLevel = 3;
  Eigen::VectorXd adaptiveSignal = adaptiveBasisSelection(signal, lowPassFilter, highPassFilter, energyThreshold, maxLevel);
  std::cout << "Adaptive Basis Selected Signal: " << adaptiveSignal.transpose() << std::endl;

  std::vector<Eigen::VectorXd> leafNodes = recursivePacketDecomposition(signal, lowPassFilter, highPassFilter, maxLevel);
  std::cout << "Number of leaf nodes: " << leafNodes.size() << std::endl;

  Eigen::VectorXd reconstructedSignal = waveletPacketReconstruct(packets[0], packets[1], synthesisLowFilter, synthesisHighFilter);
  std::cout << "Reconstructed Signal: " << reconstructedSignal.transpose() << std::endl;
    
  double rmse = computeReconstructionRMSE(signal, reconstructedSignal);
  std::cout << "Reconstruction RMSE: " << rmse << std::endl;

  Eigen::VectorXd convResult = parallelConvolve(signal, lowPassFilter);
  std::cout << "Parallel Convolution Result: " << convResult.transpose() << std::endl;

  if (shouldDecomposeFurther(energyLow, energyThreshold)) {
    std::cout << "Further decomposition required for the approximation component." << std::endl;
  } else {
    std::cout << "No further decomposition needed for the approximation component." << std::endl;
  }

  return 0;
}
