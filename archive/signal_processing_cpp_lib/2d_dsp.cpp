#include <iostream>
#include <Eigen/Dense>
#include <cmath>

// Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 2D Convolution with kernel flipping
Eigen::MatrixXd convolve2D(const Eigen::MatrixXd &input, const Eigen::MatrixXd &kernel) {
  if (input.size() == 0 || kernel.size() == 0) {
      return Eigen::MatrixXd();
  }
  int rows = input.rows();
  int cols = input.cols();
  int kRows = kernel.rows();
  int kCols = kernel.cols();
  int padRows = kRows / 2;
  int padCols = kCols / 2;
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      double sum = 0.0;
      for (int m = 0; m < kRows; ++m) {
        for (int n = 0; n < kCols; ++n) {
          int rowIndex = i - m + padRows;
          int colIndex = j - n + padCols;
          if (rowIndex >= 0 && rowIndex < rows && colIndex >= 0 && colIndex < cols) {
            sum += input(rowIndex, colIndex) * kernel(kRows - 1 - m, kCols - 1 - n);
          }
        }
      }
      output(i, j) = sum;
    }
  }
  return output;
}

// Create a normalized Gaussian kernel
Eigen::MatrixXd createGaussianKernel(int kernelSize, double sigma) {
  if (kernelSize <= 0 || sigma <= 0) {
      return Eigen::MatrixXd();
  }
  Eigen::MatrixXd kernel(kernelSize, kernelSize);
  int halfSize = kernelSize / 2;
  double sum = 0.0;
  double s = 2.0 * sigma * sigma;
  
  for (int x = -halfSize; x <= halfSize; ++x) {
    for (int y = -halfSize; y <= halfSize; ++y) {
      double r = x * x + y * y;
      double value = std::exp(-r / s) / (M_PI * s);
      kernel(x + halfSize, y + halfSize) = value;
      sum += value;
    }
  }
  if (sum > 1e-9) {
    kernel /= sum; // Normalize the kernel
  }
  return kernel;
}

// 2D Correlation without kernel flipping
Eigen::MatrixXd correlate2D(const Eigen::MatrixXd &input, const Eigen::MatrixXd &kernel) {
  if (input.size() == 0 || kernel.size() == 0) {
      return Eigen::MatrixXd();
  }
  int rows = input.rows();
  int cols = input.cols();
  int kRows = kernel.rows();
  int kCols = kernel.cols();
  int padRows = kRows / 2;
  int padCols = kCols / 2;
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      double sum = 0.0;
      for (int m = 0; m < kRows; ++m) {
        for (int n = 0; n < kCols; ++n) {
          int rowIndex = i + m - padRows;
          int colIndex = j + n - padCols;
          if (rowIndex >= 0 && rowIndex < rows && colIndex >= 0 && colIndex < cols) {
            sum += input(rowIndex, colIndex) * kernel(m, n);
          }
        }
      }
      output(i, j) = sum;
    }
  }
  return output;
}

// Compute the Sobel gradient magnitude
Eigen::MatrixXd computeSobelGradient(const Eigen::MatrixXd &image) {
  if (image.size() == 0) {
      return Eigen::MatrixXd();
  }
  int rows = image.rows();
  int cols = image.cols();
  Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(rows, cols);
  
  Eigen::Matrix3d sobelX;
  sobelX << -1, 0, 1, -2, 0, 2, -1, 0, 1;
            
  Eigen::Matrix3d sobelY;
  sobelY << -1, -2, -1, 0,  0,  0, 1,  2,  1;
             
  int pad = 1;
  for (int i = pad; i < rows - pad; ++i) {
    for (int j = pad; j < cols - pad; ++j) {
      double gx = 0.0;
      double gy = 0.0;
      for (int m = -pad; m <= pad; ++m) {
        for (int n = -pad; n <= pad; ++n) {
          double pixel = image(i + m, j + n);
          gx += pixel * sobelX(m + pad, n + pad);
          gy += pixel * sobelY(m + pad, n + pad);
        }
      }
      gradient(i, j) = std::sqrt(gx * gx + gy * gy);
    }
  }
  return gradient;
}

// Downsampling image using nearest-neighbor interpolation
Eigen::MatrixXd downsampleImage(const Eigen::MatrixXd &image, int factor) {
  if (factor <= 0 || image.size() == 0) {
      return Eigen::MatrixXd();
  }
  int originalRows = image.rows();
  int originalCols = image.cols();
  int newRows = originalRows / factor;
  int newCols = originalCols / factor;
  Eigen::MatrixXd downsampled = Eigen::MatrixXd::Zero(newRows, newCols);
  for (int i = 0; i < newRows; ++i) {
    for (int j = 0; j < newCols; ++j) {
      downsampled(i, j) = image(i * factor, j * factor);
    }
  }
  return downsampled;
}

// Apply Laplacian filter for edge detection
Eigen::MatrixXd applyLaplacianFilter(const Eigen::MatrixXd &image) {
  if (image.size() == 0) {
      return Eigen::MatrixXd();
  }
  Eigen::Matrix3d laplacian;
  laplacian << 0,  1,  0, 1, -4,  1, 0,  1,  0;
  int rows = image.rows();
  int cols = image.cols();
  int pad = 1;
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(rows, cols);
  for (int i = pad; i < rows - pad; ++i) {
    for (int j = pad; j < cols - pad; ++j) {
      double sum = 0.0;
      for (int m = -pad; m <= pad; ++m) {
        for (int n = -pad; n <= pad; ++n) {
          sum += image(i + m, j + n) * laplacian(m + pad, n + pad);
        }
      }
      output(i, j) = sum;
    }
  }
  return output;
}

// Apply Gaussian blur for noise reduction
Eigen::MatrixXd applyGaussianBlur(const Eigen::MatrixXd &image, int kernelSize, double sigma) {
  if (image.size() == 0) {
      return Eigen::MatrixXd();
  }
  Eigen::MatrixXd kernel = createGaussianKernel(kernelSize, sigma);
  return correlate2D(image, kernel);
}

int main() {
  Eigen::MatrixXd image = Eigen::MatrixXd::Random(10, 10);
  image = (image + Eigen::MatrixXd::Constant(10, 10, 1.0)) / 2.0;
  std::cout << "Original Image:" << std::endl;
  std::cout << image << std::endl << std::endl;
  
  Eigen::MatrixXd avgKernel = Eigen::MatrixXd::Constant(3, 3, 1.0 / 9.0);
  Eigen::MatrixXd convResult = convolve2D(image, avgKernel);
  std::cout << "Result of 2D Convolution (Averaging Filter):" << std::endl;
  std::cout << convResult << std::endl << std::endl;
    
  Eigen::MatrixXd gaussianKernel = createGaussianKernel(5, 1.0);
  std::cout << "Gaussian Kernel (5x5, sigma=1.0):" << std::endl;
  std::cout << gaussianKernel << std::endl << std::endl;
  
  Eigen::MatrixXd corrResult = correlate2D(image, avgKernel);
  std::cout << "Result of 2D Correlation (with Averaging Kernel):" << std::endl;
  std::cout << corrResult << std::endl << std::endl;
    
  Eigen::MatrixXd sobelGradient = computeSobelGradient(image);
  std::cout << "Sobel Gradient Magnitude:" << std::endl;
  std::cout << sobelGradient << std::endl << std::endl;
  
  Eigen::MatrixXd downsampled = downsampleImage(image, 2);
  std::cout << "Downsampled Image (factor = 2):" << std::endl;
  std::cout << downsampled << std::endl << std::endl;
  
  Eigen::MatrixXd laplacianEdges = applyLaplacianFilter(image);
  std::cout << "Laplacian Filtered Image:" << std::endl;
  std::cout << laplacianEdges << std::endl << std::endl;
  
  Eigen::MatrixXd blurredImage = applyGaussianBlur(image, 5, 1.0);
  std::cout << "Gaussian Blurred Image:" << std::endl;
  std::cout << blurredImage << std::endl << std::endl;

  return 0;
}
