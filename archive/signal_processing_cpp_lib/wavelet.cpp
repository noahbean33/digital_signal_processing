#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>

using std::vector;
using std::pair;
using std::cout;
using std::endl;
using std::make_pair;
using std::max;
using std::fabs;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

std::pair<std::vector<double>, std::vector<double>> dwtSingleLevel(
    const std::vector<double>& signal,
    const std::vector<double>& lowFilter,
    const std::vector<double>& highFilter)
{
  if (lowFilter.empty() || highFilter.empty() || lowFilter.size() != highFilter.size()) {
      throw std::invalid_argument("Filters must be non-empty and have the same size.");
  }
  std::vector<double> approx;
  std::vector<double> detail;
  for (size_t i = 0; i + lowFilter.size() <= signal.size(); i += 2)
  {
    double a = 0.0;
    double d = 0.0;
    for (size_t k = 0; k < lowFilter.size(); ++k)
    {
      a += signal[i + k] * lowFilter[k];
      d += signal[i + k] * highFilter[k];
    }
    approx.push_back(a);
    detail.push_back(d);
  }
  return std::make_pair(approx, detail);
}

std::vector<double> idwtSingleLevel(
    const std::vector<double>& approx,
    const std::vector<double>& detail,
    const std::vector<double>& recLowFilter,
    const std::vector<double>& recHighFilter)
{
  if (approx.empty() || detail.empty() || recLowFilter.empty() || recHighFilter.empty()) {
      return {};
  }
  size_t len = approx.size();
  size_t filterLen = std::max(recLowFilter.size(), recHighFilter.size());
  size_t outSize = 2 * len + filterLen - 2;
  std::vector<double> reconstructed(outSize, 0.0);

  for (size_t i = 0; i < len; ++i)
  {
    size_t index = 2 * i;
    for (size_t k = 0; k < recLowFilter.size(); ++k)
    {
      if (index + k < reconstructed.size())
        reconstructed[index + k] += approx[i] * recLowFilter[k];
    }
    for (size_t k = 0; k < recHighFilter.size(); ++k)
    {
      if (index + k < reconstructed.size())
        reconstructed[index + k] += detail[i] * recHighFilter[k];
    }
  }
  return reconstructed;
}

std::vector<std::pair<std::vector<double>, std::vector<double>>>
  dwtMultiLevel(
    const std::vector<double>& signal,
    const std::vector<double>& lowFilter,
    const std::vector<double>& highFilter,
    int levels)
{
  if (levels <= 0) {
      return {};
  }
  std::vector<std::pair<std::vector<double>, std::vector<double>>> decomposition;
  std::vector<double> currentSignal = signal;
  for (int lvl = 0; lvl < levels; ++lvl)
  {
    auto result = dwtSingleLevel(currentSignal, lowFilter, highFilter);
    decomposition.push_back(result);
    currentSignal = result.first;
  }
  return decomposition;
}

std::vector<double> thresholdCoefficients(const std::vector<double>& coeffs, double threshold)
{
  std::vector<double> thresholded;
  thresholded.reserve(coeffs.size());
  for (double c : coeffs)
  {
    double absVal = std::fabs(c);
    if (absVal < threshold)
    {
      thresholded.push_back(0.0);
    }
    else
    {
      double sign = (c >= 0.0 ? 1.0 : -1.0);
      thresholded.push_back(sign * (absVal - threshold));
    }
  }
  return thresholded;
}

double kahanConvolve(const std::vector<double>& signal,
    const std::vector<double>& filter,
    size_t startIndex)
{
  if (filter.empty()) {
      return 0.0;
  }
  double sum = 0.0;
  double compensation = 0.0;
  for (size_t i = 0; i < filter.size(); ++i)
  {
    if (startIndex + i < signal.size())
    {
      double product = signal[startIndex + i] * filter[i];
      double y = product - compensation;
      double t = sum + y;
      compensation = (t - sum) - y;
      sum = t;
    }
  }
  return sum;
}

std::vector<double> multiLevelIDWT(
    const std::vector<std::pair<std::vector<double>, std::vector<double>>>& decomposition,
    const std::vector<double>& recLowFilter,
    const std::vector<double>& recHighFilter)
{
  if (decomposition.empty())
  {
    return std::vector<double>();
  }

  std::vector<double> current = idwtSingleLevel(
    decomposition.back().first,
    decomposition.back().second,
    recLowFilter,
    recHighFilter);

  for (int i = decomposition.size() - 2; i >= 0; --i)
  {
    // The size of the detail coefficients might not match the upsampled 'current' signal.
    // We need to handle this carefully. For simplicity, we assume they are compatible.
    current = idwtSingleLevel(current, decomposition[i].second, recLowFilter, recHighFilter);
  }
  return current;
}

int main()
{
  std::vector<double> signal = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
  std::vector<double> lowFilter = { 0.70710678, 0.70710678 };
  std::vector<double> highFilter = { -0.70710678, 0.70710678 };
  std::vector<double> recLowFilter = { 0.70710678, 0.70710678 };
  std::vector<double> recHighFilter = { 0.70710678, -0.70710678 };

  auto dwtResult = dwtSingleLevel(signal, lowFilter, highFilter);
  std::cout << "Single-Level DWT Approximation: ";
  for (double a : dwtResult.first) { std::cout << a << " "; }
  std::cout << std::endl;

  std::cout << "Single-Level DWT Detail: ";
  for (double d : dwtResult.second) { std::cout << d << " "; }
  std::cout << std::endl;

  int levels = 2;
  auto decomposition = dwtMultiLevel(signal, lowFilter, highFilter, levels);

  const double THRESHOLD = 0.5;
std::vector<double> thresholdedDetail = thresholdCoefficients(decomposition[0].second, THRESHOLD);
  std::cout << "Thresholded Detail Coefficients (Level 1): ";
  for (double t : thresholdedDetail) { std::cout << t << " "; }
  std::cout << std::endl;

  std::vector<double> reconstructed = multiLevelIDWT(decomposition, recLowFilter, recHighFilter);
  std::cout << "Reconstructed Signal: ";
  for (double r : reconstructed) { std::cout << r << " "; }
  std::cout << std::endl;

  double convResult = kahanConvolve(signal, lowFilter, 0);
  std::cout << "Kahan-compensated convolution (start at index 0): " << convResult << std::endl;

  return 0;
}
