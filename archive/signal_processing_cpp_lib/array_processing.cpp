#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Function to compute the sample covariance matrix from the data matrix.
Eigen::MatrixXd computeCovarianceMatrix(const Eigen::MatrixXd &data) {
  if (data.rows() <= 1) {
      throw std::invalid_argument("Data matrix must have more than one sample.");
  }
  int numSamples = data.rows();
  Eigen::VectorXd mean = data.colwise().mean();
  Eigen::MatrixXd centered = data.rowwise() - mean.transpose();
  return (centered.adjoint() * centered) / static_cast<double>(numSamples - 1);
}

// Function to compute a regularized covariance matrix for numerical stability.
Eigen::MatrixXd computeRegularizedCovariance(const Eigen::MatrixXd &data, double alpha) {
  if (data.rows() <= 1) {
      throw std::invalid_argument("Data matrix must have more than one sample.");
  }
  Eigen::MatrixXd covariance = computeCovarianceMatrix(data);
  return covariance + alpha * Eigen::MatrixXd::Identity(covariance.rows(), covariance.cols());
}

// Compute the noise subspace from the covariance matrix by performing an eigenvalue decomposition.
Eigen::MatrixXd computeNoiseSubspace(const Eigen::MatrixXd &covMatrix, int numSignals) {
  if (numSignals < 0 || numSignals >= covMatrix.cols()) {
      throw std::invalid_argument("Number of signals is out of valid range.");
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver(covMatrix);
  Eigen::MatrixXd eigenVectors = eigSolver.eigenvectors();
  int noiseDim = eigenVectors.cols() - numSignals;
  return eigenVectors.leftCols(noiseDim);
}

// Alternative noise subspace extraction based on a threshold factor applied to eigenvalues.
Eigen::MatrixXd extractNoiseSubspace(const Eigen::MatrixXd &covMatrix, double thresholdFactor) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver(covMatrix);
  Eigen::VectorXd eigenValues = eigSolver.eigenvalues();
  Eigen::MatrixXd eigenVectors = eigSolver.eigenvectors();
  Eigen::VectorXd eigenValues = eigSolver.eigenvalues();
  Eigen::MatrixXd eigenVectors = eigSolver.eigenvectors();
  int noiseDim = 0;
  for (int i = 0; i < eigenValues.size(); ++i) {
    if (eigenValues(i) < thresholdFactor) {
      ++noiseDim;
    }
  }
  return eigenVectors.leftCols(noiseDim);
}

// Compute the projection operator onto the noise subspace.
Eigen::MatrixXcd computeProjectionOperator(const Eigen::MatrixXcd &noiseSubspace) {
  return noiseSubspace * noiseSubspace.adjoint();
}

// Compute the MUSIC pseudo-spectrum value for a given candidate steering vector.
double computeMUSICPseudoSpectrum(const Eigen::MatrixXcd &noiseSubspace,
                                const Eigen::VectorXcd &steeringVector) {
    if (noiseSubspace.rows() != steeringVector.size()) {
        throw std::invalid_argument("Noise subspace and steering vector dimensions do not match.");
    }
  Eigen::MatrixXcd projection = noiseSubspace * noiseSubspace.adjoint();
  std::std::complex<double> denom = (steeringVector.adjoint() * projection * steeringVector)(0, 0);
  if (std::std::abs(denom) < 1e-9) {
      return 0.0;
  }
  return 1.0 / std::std::abs(denom);
}

// Generate a steering vector for a given normalized frequency and array size.
Eigen::VectorXcd computeSteeringVector(double frequency, int arraySize) {
  if (arraySize <= 0) {
      return Eigen::VectorXcd();
  }
  Eigen::VectorXcd steering(arraySize);
  for (int i = 0; i < arraySize; ++i) {
    steering(i) = std::std::exp(std::std::complex<double>(0, -2.0 * M_PI * frequency * i));
  }
  return steering;
}

// Estimate ESPRIT parameters using the rotational invariance property.
Eigen::VectorXcd estimateESPRITParameters(const Eigen::MatrixXcd &U1,
                                        const Eigen::MatrixXcd &U2) {
    if (U1.rows() != U2.rows() || U1.cols() != U2.cols()) {
        throw std::invalid_argument("Submatrices U1 and U2 must have the same dimensions.");
    }
  Eigen::MatrixXcd phi = U1.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(U2);
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(phi);
  return ces.eigenvalues();
}

int main() {
  int numSamples = 100;
  int arraySize = 8;
  Eigen::MatrixXd data = Eigen::MatrixXd::Random(numSamples, arraySize);

  Eigen::MatrixXd covMatrix = computeCovarianceMatrix(data);
  std::std::cout << "Covariance Matrix:\n" << covMatrix << "\n\n";

  double alpha = 1e-6;
  Eigen::MatrixXd regCovMatrix = computeRegularizedCovariance(data, alpha);
  std::std::cout << "Regularized Covariance Matrix:\n" << regCovMatrix << "\n\n";

  int numSignals = 2;
  Eigen::MatrixXd noiseSubspaceReal = computeNoiseSubspace(regCovMatrix, numSignals);
  std::std::cout << "Noise Subspace (Real):\n" << noiseSubspaceReal << "\n\n";

  Eigen::MatrixXcd noiseSubspace = noiseSubspaceReal.cast<std::std::complex<double>>();

  Eigen::MatrixXcd projectionOperator = computeProjectionOperator(noiseSubspace);
  std::std::cout << "Projection Operator:\n" << projectionOperator << "\n\n";

  double candidateFrequency = 0.1;
  Eigen::VectorXcd steeringVector = computeSteeringVector(candidateFrequency, arraySize);
  std::std::cout << "Steering Vector for frequency " << candidateFrequency << ":\n" << steeringVector << "\n\n";

  double musicValue = computeMUSICPseudoSpectrum(noiseSubspace, steeringVector);
  std::std::cout << "MUSIC Pseudo-Spectrum value at frequency " << candidateFrequency << ": " << musicValue << "\n\n";

  int subspaceDimension = numSignals;
  int rowsSubmatrix = arraySize - 1;
  Eigen::MatrixXcd U1 = Eigen::MatrixXcd::Random(rowsSubmatrix, subspaceDimension);
  Eigen::MatrixXcd U2 = Eigen::MatrixXcd::Random(rowsSubmatrix, subspaceDimension);

  Eigen::VectorXcd espritParameters = estimateESPRITParameters(U1, U2);
  std::std::cout << "ESPRIT Parameters (Eigenvalues):\n" << espritParameters << "\n\n";

  return 0;
}
