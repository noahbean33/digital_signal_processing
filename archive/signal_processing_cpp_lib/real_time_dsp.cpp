#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <mutex>
#include <atomic>
#include <thread>
#include <chrono>
#include <algorithm>

using std::vector;
using std::cout;
using std::endl;
using std::abs;
using std::sin;
using std::mutex;
using std::atomic;
using std::thread;
using std::chrono::milliseconds;
using std::this_thread::sleep_for;
using std::out_of_range;
using std::invalid_argument;
using std::exception;
using std::cerr;
using std::lock_guard;
using std::fill;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void addSampleToCircularBuffer(std::vector<double>& buffer, int& writeIndex, double sample) {
  if (buffer.empty()) {
      return;
  }
  buffer[writeIndex] = sample;
  writeIndex = (writeIndex + 1) % buffer.size();
}

bool isBufferOverloaded(int currentOccupancy, int bufferCapacity, double latencyThreshold, double sampleRate) {
  if (sampleRate <= 0) {
      return false;
  }
  double estimatedLatency = static_cast<double>(currentOccupancy) / sampleRate;
  return estimatedLatency > latencyThreshold;
}

double computeTimeOffset(int sampleIndex1, int sampleIndex2, double sampleRate) {
  if (sampleRate <= 0) {
      return 0.0;
  }
  return std::abs(sampleIndex1 - sampleIndex2) / sampleRate;
}

void threadSafeAddSample(std::vector<double>& buffer, std::atomic<int>& writeIndex, double sample) {
  if (buffer.empty()) {
      return;
  }
  int currentWriteIndex = writeIndex.fetch_add(1);
  buffer[currentWriteIndex % buffer.size()] = sample;
}

void resetCircularBuffer(std::vector<double>& buffer, int& writeIndex, double initialValue = 0.0) {
  std::fill(buffer.begin(), buffer.end(), initialValue);
  writeIndex = 0;
}

double safeGetSample(const std::vector<double>& buffer, int index) {
  if (index < 0 || index >= static_cast<int>(buffer.size())) {
    throw std::out_of_range("Index out of bounds in circular buffer access");
  }
  return buffer[index];
}

void processBufferSample(const std::vector<double>& buffer, int index) {
  try {
    double sample = safeGetSample(buffer, index);
    std::cout << "Processed sample at index " << index << ": " << sample << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "Error processing sample at index " << index << ": " << e.what() << std::endl;
  }
}

void resizeCircularBuffer(std::vector<double>& buffer, size_t newCapacity) {
  if (newCapacity < 1) {
    throw std::invalid_argument("Buffer capacity must be at least 1");
  }
  buffer.resize(newCapacity, 0.0);
}

int main() {
  size_t bufferSize = 10;
  std::vector<double> circularBuffer(bufferSize, 0.0);
  int writeIndex = 0;

  for (int i = 0; i < 15; ++i) {
    double sample = std::sin(2 * M_PI * i / bufferSize);
    addSampleToCircularBuffer(circularBuffer, writeIndex, sample);
  }

  int currentOccupancy = 7;
  double sampleRate = 1000.0;
  const double LATENCY_THRESHOLD_S = 0.005;
double latencyThreshold = LATENCY_THRESHOLD_S;
  if (isBufferOverloaded(currentOccupancy, circularBuffer.size(), latencyThreshold, sampleRate)) {
    std::cout << "Buffer overloaded! Adjust the processing pipeline." << std::endl;
  } else {
    std::cout << "Buffer latency is within acceptable limits." << std::endl;
  }

  int sampleIndex1 = 3;
  int sampleIndex2 = 7;
  double timeOffset = computeTimeOffset(sampleIndex1, sampleIndex2, sampleRate);
  std::cout << "Time offset between streams: " << timeOffset << " seconds" << std::endl;

  for (int i = 0; i < 15; ++i) {
    double sample = std::sin(2 * M_PI * i / bufferSize);
    addSampleToCircularBuffer(circularBuffer, writeIndex, sample);
  }

  std::atomic<int> atomicWriteIndex(writeIndex);
  std::thread t1([&]() {
    for (int i = 0; i < 5; ++i) {
      threadSafeAddSample(circularBuffer, atomicWriteIndex, i * 0.1);
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  });
  std::thread t2([&]() {
    for (int i = 0; i < 5; ++i) {
      threadSafeAddSample(circularBuffer, atomicWriteIndex, i * 0.2);
      std::this_thread::sleep_for(std::chrono::milliseconds(15));
    }
  });
  t1.join();
  t2.join();

  // update the non-atomic writeIndex after threads are done
  writeIndex = atomicWriteIndex.load();

  for (int i = 0; i < static_cast<int>(circularBuffer.size()); ++i) {
    processBufferSample(circularBuffer, i);
  }

  resetCircularBuffer(circularBuffer, writeIndex, 0.0);
  std::cout << "Circular buffer has been reset." << std::endl;

  try {
    resizeCircularBuffer(circularBuffer, 15);
    std::cout << "Buffer resized to capacity " << circularBuffer.size() << "." << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "Error resizing buffer: " << e.what() << std::endl;
  }

  return 0;
}
