#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <emmintrin.h>
#include <numeric>
#include <stdexcept>

// Template function to measure the execution time of any callable
template <typename Func>
double measureExecutionTime(Func f) {
  auto start = std::chrono::high_resolution_clock::now();
  f();
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;
  return duration.count();
}

// Optimized vector addition using pre-allocation
std::vector<double> optimizedVectorAddition(const std::vector<double>& A, const std::vector<double>& B) {
  if (A.size() != B.size()) {
      throw std::invalid_argument("Vectors must have the same size.");
  }
  std::vector<double> result;
  result.reserve(A.size());
  for (std::size_t i = 0; i < A.size(); ++i) {
    result.push_back(A[i] + B[i]);
  }
  return result;
}

// Loop unrolling for vector addition
std::vector<double> unrolledVectorAddition(const std::vector<double>& A, const std::vector<double>& B) {
  if (A.size() != B.size()) {
      throw std::invalid_argument("Vectors must have the same size.");
  }
  std::size_t n = A.size();
  std::vector<double> result(n);
  std::size_t i = 0;
  for (; i + 3 < n; i += 4) {
    result[i] = A[i] + B[i];
    result[i + 1] = A[i + 1] + B[i + 1];
    result[i + 2] = A[i + 2] + B[i + 2];
    result[i + 3] = A[i + 3] + B[i + 3];
  }
  for (; i < n; ++i) {
    result[i] = A[i] + B[i];
  }
  return result;
}

// Parallel vector addition using std::thread
std::vector<double> parallelVectorAddition(const std::vector<double>& A, const std::vector<double>& B) {
  if (A.size() != B.size()) {
      throw std::invalid_argument("Vectors must have the same size.");
  }
  std::size_t n = A.size();
  std::vector<double> result(n, 0.0);
  unsigned int numThreads = std::thread::hardware_concurrency();
  if (numThreads == 0) numThreads = 2;
  std::vector<std::thread> threads;
  std::size_t blockSize = (n + numThreads - 1) / numThreads;

  auto worker = [&](std::size_t start, std::size_t end) {
    for (std::size_t i = start; i < end; ++i) {
      result[i] = A[i] + B[i];
    }
  };

  for (unsigned int t = 0; t < numThreads; ++t) {
    std::size_t start = t * blockSize;
    if (start >= n) break;
    std::size_t end = std::min(start + blockSize, n);
    threads.emplace_back(worker, start, end);
  }

  for (auto& thread : threads) {
    thread.join();
  }
  return result;
}

// SIMD dot product computation using SSE2 intrinsics
double simdDotProduct(const double* a, const double* b, std::size_t n) {
  std::size_t i = 0;
  __m128d sumVec = _mm_setzero_pd();
  for (; i + 1 < n; i += 2) {
    __m128d vecA = _mm_loadu_pd(a + i);
    __m128d vecB = _mm_loadu_pd(b + i);
    sumVec = _mm_add_pd(sumVec, _mm_mul_pd(vecA, vecB));
  }
  double temp[2];
  _mm_storeu_pd(temp, sumVec);
  double dot = temp[0] + temp[1];
  for (; i < n; ++i) {
    dot += a[i] * b[i];
  }
  return dot;
}

// Static loop unrolling via template recursion
template <int N>
struct StaticUnrollSum {
  static double sum(const double* arr) {
    return arr[N - 1] + StaticUnrollSum<N - 1>::sum(arr);
  }
};

template <>
struct StaticUnrollSum<1> {
  static double sum(const double* arr) {
    return arr[0];
  }
};

// Inline multiplication function
inline double forceInlineMultiply(double a, double b) {
  return a * b;
}

// Dynamic frequency scaling function
double adjustFrequency(double load, double baseFrequency) {
  const double LOW_LOAD_THRESHOLD = 0.5;
  const double HIGH_LOAD_THRESHOLD = 0.8;
  const double LOW_FREQ_SCALE = 0.8;
  const double HIGH_FREQ_SCALE = 1.2;
  if (load < LOW_LOAD_THRESHOLD) {
    return baseFrequency * LOW_FREQ_SCALE;
  } else if (load < HIGH_LOAD_THRESHOLD) {
    return baseFrequency;
  } else {
    return baseFrequency * HIGH_FREQ_SCALE;
  }
}

// Function to decide if further optimization is needed
bool shouldOptimizeFurther(double execTime, double threshold) {
  return execTime > threshold;
}

int main() {
  std::vector<double> A(16), B(16);
  std::iota(A.begin(), A.end(), 0.0);
  std::iota(B.begin(), B.end(), 0.0);
  for(auto& val : B) { val *= 2; }

  auto unrolledResult = unrolledVectorAddition(A, B);
  std::cout << "Unrolled Vector Addition Result:" << std::endl;
  for (const auto& value : unrolledResult) { std::cout << value << " "; }
  std::cout << std::endl;

  auto parallelResult = parallelVectorAddition(A, B);
  std::cout << "Parallel Vector Addition Result:" << std::endl;
  for (const auto& value : parallelResult) { std::cout << value << " "; }
  std::cout << std::endl;

  double dot = simdDotProduct(A.data(), B.data(), A.size());
  std::cout << "SIMD Dot Product: " << dot << std::endl;

  double execTime = measureExecutionTime([&]() {
    auto result = unrolledVectorAddition(A, B);
    volatile double sum = 0.0;
    for (const auto& val : result) { sum += val; }
  });
  std::cout << "Execution Time (Unrolled Addition): " << execTime << " seconds" << std::endl;

  double arr[4] = {1.0, 2.0, 3.0, 4.0};
  double staticSum = StaticUnrollSum<4>::sum(arr);
  std::cout << "Static Unroll Sum: " << staticSum << std::endl;

  double multiplyResult = forceInlineMultiply(3.5, 2.0);
  std::cout << "Force Inline Multiply: " << multiplyResult << std::endl;

  const double BASE_FREQ = 2.0;
  const double LOAD = 0.85;
  double adjustedFreq = adjustFrequency(LOAD, BASE_FREQ);
  std::cout << "Adjusted Frequency: " << adjustedFreq << " (relative units)" << std::endl;

  const double THRESHOLD = 0.001;
  bool optimizeMore = shouldOptimizeFurther(execTime, THRESHOLD);
  std::cout << "Should Optimize Further: " << (optimizeMore ? "Yes" : "No") << std::endl;

  return 0;
}
