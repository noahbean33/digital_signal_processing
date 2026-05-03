#include <iostream>
#include <vector>
#include <cmath>

using std::vector;
using std::cout;
using std::endl;
using std::sin;
using std::cos;
using std::fabs;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double sinc(double x) {
  if (std::fabs(x) < 1e-7) {
    return 1.0;
  }
  return std::sin(M_PI * x) / (M_PI * x);
}

std::vector<double> generateFIRLowpass(int numTaps, double cutoff) {
  if (numTaps <= 1) {
      return {};
  }
  std::vector<double> kernel(numTaps, 0.0);
  int M = numTaps - 1;
  for (int i = 0; i < numTaps; ++i) {
    double n = i - M / 2.0;
    kernel[i] = 2 * cutoff * sinc(2 * cutoff * n);
    kernel[i] *= (0.54 - 0.46 * std::cos((2 * M_PI * i) / M));
  }
  return kernel;
}

std::vector<double> decimateSignal(const std::vector<double>& input, int factor) {
  if (factor <= 0) {
      return {};
  }
  std::vector<double> output;
  if (!input.empty()) {
      output.reserve(input.size() / factor);
      for (size_t i = 0; i < input.size(); i += factor) {
        output.push_back(input[i]);
      }
  }
  return output;
}

std::vector<double> interpolateSignal(const std::vector<double>& input, int factor) {
  if (factor <= 0) {
      return {};
  }
  std::vector<double> output(input.size() * factor, 0.0);
  for (size_t i = 0; i < input.size(); ++i) {
    output[i * factor] = input[i];
  }
  return output;
}

std::vector<double> applyFIRFilter(const std::vector<double>& signal, const std::vector<double>& kernel) {
  if (signal.empty() || kernel.empty()) {
      return {};
  }
  int signalSize = signal.size();
  int kernelSize = kernel.size();
  std::vector<double> filtered(signalSize, 0.0);
  for (int i = 0; i < signalSize; ++i) {
    double sum = 0.0;
    for (int j = 0; j < kernelSize; ++j) {
      if (i - j >= 0) {
        sum += signal[i - j] * kernel[j];
      }
    }
    filtered[i] = sum;
  }
  return filtered;
}

std::vector<double> sampleRateConversion(const std::vector<double>& input, int L, int M, const std::vector<double>& firKernel) {
  if (L <= 0 || M <= 0) {
      return {};
  }
  std::vector<double> upsampled = interpolateSignal(input, L);
  std::vector<double> filtered = applyFIRFilter(upsampled, firKernel);
  std::vector<double> output = decimateSignal(filtered, M);
  return output;
}

int main() {
  double f_s = 1000.0;
  int N = 100;
  double frequency = 50.0;

  std::vector<double> signal(N, 0.0);
  for (int n = 0; n < N; ++n) {
    signal[n] = std::sin(2 * M_PI * frequency * n / f_s);
  }

  int L = 3;
  int M = 2;

  int numTaps = 29;
  double normalizedCutoff = 0.25;
  std::vector<double> firKernel = generateFIRLowpass(numTaps, normalizedCutoff);

  std::vector<double> convertedSignal = sampleRateConversion(signal, L, M, firKernel);

  std::cout << "Original Signal Size: " << signal.size() << std::endl;
  std::cout << "Converted Signal Size: " << convertedSignal.size() << std::endl;
  std::cout << "New Sampling Frequency: " << (static_cast<double>(L) / M) * f_s << " Hz" << std::endl;

  std::cout << "\nFirst 10 samples of Converted Signal:" << std::endl;
  for (size_t i = 0; i < 10 && i < convertedSignal.size(); ++i) {
    std::cout << convertedSignal[i] << std::endl;
  }

  return 0;
}
