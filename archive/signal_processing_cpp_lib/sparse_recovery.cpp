#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

// Function: applyThreshold
Eigen::VectorXd applyThreshold(const Eigen::VectorXd &signal, double threshold) {
  Eigen::VectorXd sparse = signal;
  for (int i = 0; i < sparse.size(); ++i) {
    if (std::abs(sparse(i)) < threshold) {
      sparse(i) = 0.0;
    }
  }
  return sparse;
}

// Function: computeL1Norm
double computeL1Norm(const Eigen::VectorXd &vector) {
  return vector.lpNorm<1>();
}

// Function: softThreshold
Eigen::VectorXd softThreshold(const Eigen::VectorXd &x, double lambda) {
  Eigen::VectorXd result = x;
  for (int i = 0; i < x.size(); ++i) {
    if (x(i) > lambda)
      result(i) = x(i) - lambda;
    else if (x(i) < -lambda)
      result(i) = x(i) + lambda;
    else
      result(i) = 0.0;
  }
  return result;
}

// Function: computeResidualNorm
double computeResidualNorm(const Eigen::MatrixXd &A, const Eigen::VectorXd &x,
                         const Eigen::VectorXd &y) {
  return (A * x - y).norm();
}

// Function: iterativeSoftThresholding
Eigen::VectorXd iterativeSoftThresholding(const Eigen::MatrixXd &A,
                                        const Eigen::VectorXd &y,
                                        double lambda,
                                        double tau,
                                        int maxIter) {
  Eigen::VectorXd x = Eigen::VectorXd::Zero(A.cols());
  for (int iter = 0; iter < maxIter; ++iter) {
    Eigen::VectorXd grad = A.transpose() * (A * x - y);
    x = softThreshold(x - tau * grad, lambda);
  }
  return x;
}

// Function: orthogonalMatchingPursuit
Eigen::VectorXd orthogonalMatchingPursuit(const Eigen::MatrixXd &A,
                                        const Eigen::VectorXd &y,
                                        int sparsity) {
  int n = A.cols();
  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  std::vector<int> support;
  Eigen::VectorXd residual = y;
  for (int k = 0; k < sparsity; ++k) {
    Eigen::VectorXd correlations = A.transpose() * residual;
    int idx;
    correlations.cwiseAbs().maxCoeff(&idx);
    support.push_back(idx);
    Eigen::MatrixXd A_sub(y.size(), support.size());
    for (size_t i = 0; i < support.size(); ++i) {
      A_sub.col(i) = A.col(support[i]);
    }
    
    Eigen::VectorXd x_sub = A_sub.colPivHouseholderQr().solve(y);
    residual = y - A_sub * x_sub;
  }
  for (size_t i = 0; i < support.size(); ++i) {
    x(support[i]) = x_sub(i);
  }
  return x;
}

// Function: acceleratedIterativeSoftThresholding
Eigen::VectorXd acceleratedIterativeSoftThresholding(const Eigen::MatrixXd &A,
                                                     const Eigen::VectorXd &y,
                                                     double lambda,
                                                     double tau,
                                                     int maxIter) {
  Eigen::VectorXd x = Eigen::VectorXd::Zero(A.cols());
  Eigen::VectorXd z = x;
  double t = 1.0;
  for (int iter = 0; iter < maxIter; ++iter) {
    Eigen::VectorXd grad = A.transpose() * (A * z - y);
    Eigen::VectorXd x_new = softThreshold(z - tau * grad, lambda);
    double t_new = (1.0 + std::sqrt(1.0 + 4.0 * t * t)) / 2.0;
    z = x_new + ((t - 1.0) / t_new) * (x_new - x);
    x = x_new;
    t = t_new;
  }
  return x;
}

// Function: hybridSparseRecovery
Eigen::VectorXd hybridSparseRecovery(const Eigen::MatrixXd &A,
                                   const Eigen::VectorXd &y,
                                   int initSparsity,
                                   double lambda,
                                   double tau,
                                   int maxIter) {
  Eigen::VectorXd x_greedy = orthogonalMatchingPursuit(A, y, initSparsity);
  std::vector<int> support;
  for (int i = 0; i < x_greedy.size(); ++i) {
    if (std::abs(x_greedy(i)) > 1e-6)
      support.push_back(i);
  }
  
  Eigen::VectorXd x = x_greedy;
  for (int iter = 0; iter < maxIter; ++iter) {
    Eigen::VectorXd grad = A.transpose() * (A * x - y);
    for (int i = 0; i < grad.size(); ++i) {
      if (std::find(support.begin(), support.end(), i) == support.end())
        grad(i) = 0;
    }
    x = softThreshold(x - tau * grad, lambda);
  }
  return x;
}

// Function: parallelGradientComputation
Eigen::VectorXd parallelGradientComputation(const Eigen::MatrixXd &A,
                                          const Eigen::VectorXd &x,
                                          const Eigen::VectorXd &y) {
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(A.cols());
  Eigen::VectorXd residual = A * x - y;
  #pragma omp parallel for
  for (int i = 0; i < A.cols(); ++i) {
    double accum = 0.0;
    for (int j = 0; j < A.rows(); ++j) {
      accum += A(j, i) * residual(j);
    }
    grad(i) = accum;
  }
  return grad;
}

int main() {
  const int m = 50;
  const int n = 100;

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(m, n);

  Eigen::VectorXd x_true = Eigen::VectorXd::Zero(n);
  std::vector<int> true_support;
  std::mt19937 rng(static_cast<unsigned int>(time(nullptr)));
  std::uniform_int_distribution<int> dist(0, n - 1);
  std::uniform_real_distribution<double> val_dist(-1.0, 1.0);

  while (true_support.size() < 10) {
    int idx = dist(rng);
    if (std::find(true_support.begin(), true_support.end(), idx) == true_support.end()) {
      true_support.push_back(idx);
      x_true(idx) = val_dist(rng);
    }
  }

  Eigen::VectorXd y = A * x_true;

  double lambda = 0.1;
  double tau = 0.001;
  int maxIter = 100;
  int initSparsity = 10;

  Eigen::VectorXd x_ist = iterativeSoftThresholding(A, y, lambda, tau, maxIter);
  Eigen::VectorXd x_omp = orthogonalMatchingPursuit(A, y, initSparsity);
  Eigen::VectorXd x_acc = acceleratedIterativeSoftThresholding(A, y, lambda, tau, maxIter);
  Eigen::VectorXd x_hybrid = hybridSparseRecovery(A, y, initSparsity, lambda, tau, maxIter);

  std::cout << "Residual Norm (Iterative Soft Thresholding): " << computeResidualNorm(A, x_ist, y) << std::endl;
  std::cout << "Residual Norm (Orthogonal Matching Pursuit): " << computeResidualNorm(A, x_omp, y) << std::endl;
  std::cout << "Residual Norm (Accelerated IST): " << computeResidualNorm(A, x_acc, y) << std::endl;
  std::cout << "Residual Norm (Hybrid Recovery): " << computeResidualNorm(A, x_hybrid, y) << std::endl;

  std::cout << "L1 Norm of true sparse signal: " << computeL1Norm(x_true) << std::endl;

  return 0;
}
