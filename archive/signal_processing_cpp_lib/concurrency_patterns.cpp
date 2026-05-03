#include <iostream>
#include <vector>
#include <thread>
#include <future>
#include <mutex>
#include <atomic>
#include <functional>
#include <cmath>
#include <chrono>

// Simple DSP task: applies a constant gain to each signal sample.
void processDSPTask(const std::vector<double>& input,
                  std::vector<double>& output) {
  const double GAIN = 0.75;
  output.resize(input.size());
  for (std::size_t i = 0; i < input.size(); ++i) {
    output[i] = input[i] * GAIN;
  }
}

// DSP worker: processes a segment by modulating each element with its sine.
void dspWorker(const std::vector<double>& dataSegment,
             std::vector<double>& resultSegment) {
  resultSegment.resize(dataSegment.size());
  for (std::size_t i = 0; i < dataSegment.size(); ++i) {
    resultSegment[i] = dataSegment[i] * std::sin(dataSegment[i]);
  }
}

// Thread-safe update of a shared data buffer using mutex.
void updateSharedBuffer(std::vector<double>& sharedBuffer, double newValue, std::mutex& mtx) {
  std::lock_guard<std::mutex> lock(mtx);
  sharedBuffer.push_back(newValue);
}

// Performs a lock-free atomic increment on a shared counter.
void atomicIncrement(std::atomic<int>& counter) {
  counter.fetch_add(1, std::memory_order_relaxed);
}

// Schedules a DSP task by creating a new thread executing the provided callable.
void scheduleTask(std::function<void()> task, std::thread& worker) {
  worker = std::thread(task);
}

// Dynamically allocates tasks using an atomic counter to ensure unique task assignment.
bool tryAllocateTask(std::atomic<int>& currentTask, int totalTasks, int &allocatedTask) {
  int oldTask = currentTask.load(std::memory_order_relaxed);
  while (oldTask < totalTasks) {
    if (currentTask.compare_exchange_weak(oldTask, oldTask + 1, std::memory_order_release, std::memory_order_relaxed)) {
      allocatedTask = oldTask;
      return true;
    }
  }
  return false;
}

// Demonstrates non-blocking send using an atomic_flag to control access.
bool nonBlockingSend(std::atomic_flag &lock, std::function<void()> sendOperation) {
  if (!lock.test_and_set(std::memory_order_acquire)) {
    sendOperation();
    lock.clear(std::memory_order_release);
    return true;
  }
  return false;
}

// Performs parallel processing of a DSP operation by partitioning the data and aggregating results.
double parallelDSPCompute(const std::vector<double>& data, std::function<double(double)> process) {
  const std::size_t numThreads = std::thread::hardware_concurrency();
  if (numThreads == 0 || data.empty()) {
      return 0.0;
  }
  std::size_t blockSize = (data.size() + numThreads - 1) / numThreads;
  std::vector<std::future<double>> futures;
  for (std::size_t i = 0; i < numThreads; ++i) {
    auto start = data.begin() + i * blockSize;
    if (start >= data.end()) break; // Avoid creating threads for empty ranges
    auto end = (i == numThreads - 1) ? data.end() : start + blockSize;
    if (end > data.end()) end = data.end();

    futures.push_back(std::async(std::launch::async, [start, end, process]() {
      double partialSum = 0.0;
      for (auto it = start; it != end; ++it) {
        partialSum += process(*it);
      }
      return partialSum;
    }));
  }

  double total = 0.0;
  for (auto& f : futures) {
    total += f.get();
  }
  return total;
}

int main() {
  const int SIGNAL_SIZE = 1000;
  const double SIGNAL_FREQ = 0.01;
  std::vector<double> inputSignal(SIGNAL_SIZE);
  for (std::size_t i = 0; i < inputSignal.size(); ++i) {
    inputSignal[i] = std::sin(SIGNAL_FREQ * i);
  }

  std::vector<double> processedSignal;
  processDSPTask(inputSignal, processedSignal);

  std::vector<double> segmentResult;
  const std::size_t segmentSize = 200;
  std::vector<double> dataSegment(inputSignal.begin(), inputSignal.begin() + segmentSize);
  dspWorker(dataSegment, segmentResult);

  std::vector<double> sharedBuffer;
  std::mutex bufferMutex;
  updateSharedBuffer(sharedBuffer, 1.234, bufferMutex);
  updateSharedBuffer(sharedBuffer, 2.468, bufferMutex);

  std::atomic<int> atomicCounter(0);
  atomicIncrement(atomicCounter);
  std::cout << "Atomic counter value: " << atomicCounter.load() << std::endl;

  std::thread scheduledThread;
  scheduleTask([&]() {
    std::cout << "Scheduled DSP task executed in thread id: " << std::this_thread::get_id() << std::endl;
  }, scheduledThread);
  if (scheduledThread.joinable()) {
    scheduledThread.join();
  }

  std::atomic<int> currentTask(0);
  const int totalTasks = 10;
  const int NUM_WORKERS = 4;
  std::vector<std::thread> taskThreads;
  for (int i = 0; i < NUM_WORKERS; ++i) {
    taskThreads.emplace_back([&]() {
      int allocatedTask;
      while (tryAllocateTask(currentTask, totalTasks, allocatedTask)) {
        std::cout << "Thread " << std::this_thread::get_id() << " processing task " << allocatedTask << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
      }
    });
  }
  
  for (auto &t : taskThreads) {
    t.join();
  }

  std::atomic_flag sendLock = ATOMIC_FLAG_INIT;
  bool sendSuccess = nonBlockingSend(sendLock, []() {
    std::cout << "Non-blocking send operation executed." << std::endl;
  });
  if (!sendSuccess) {
    std::cout << "Send operation skipped due to active lock." << std::endl;
  }

  double result = parallelDSPCompute(processedSignal, [](double value) {
    return std::abs(value);
  });
  std::cout << "Parallel DSP compute (sum of absolute values): " << result << std::endl;

  return 0;
}
