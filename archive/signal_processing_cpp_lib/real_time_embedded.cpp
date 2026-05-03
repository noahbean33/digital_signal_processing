#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include <atomic>
#include <functional>
#include <future>
#include <string>
#include <stdexcept>

// Inline function for real-time signal processing
inline void processSignalFrame(const float* input, float* output, std::size_t frameSize) {
  const float GAIN = 0.5f;
  for (std::size_t i = 0; i < frameSize; ++i) {
    output[i] = input[i] * GAIN;
  }
}

// Function to read a value from a circular (ring) buffer
int readCircularBuffer(const std::vector<int>& buffer, std::size_t& head) {
  if (buffer.empty()) {
      throw std::out_of_range("Buffer is empty.");
  }
  int data = buffer[head];
  head = (head + 1) % buffer.size();
  return data;
}

// Scheduler function that executes a given task periodically
void schedulePeriodicTask(std::function<void()> task, std::chrono::milliseconds interval, int iterations) {
  auto next = std::chrono::high_resolution_clock::now() + interval;
  for (int i = 0; i < iterations; ++i) {
    task();
    std::this_thread::sleep_until(next);
    next += interval;
  }
}

// Jitter-compensated sleep routine
inline void jitterCompensatedSleep(std::chrono::steady_clock::time_point &scheduled, const std::chrono::milliseconds period) {
  auto now = std::chrono::steady_clock::now();
  if (now < scheduled) {
    std::this_thread::sleep_until(scheduled);
  }
  auto actual = std::chrono::steady_clock::now();
  scheduled += period + std::chrono::duration_cast<std::chrono::milliseconds>(actual - scheduled);
}

struct Event {
  int id;
  std::string message;
};

// Lock-free enqueue operation for a circular buffer
bool lockFreeEnqueue(std::atomic<std::size_t>& writeIndex,
                   std::atomic<std::size_t>& readIndex,
                   std::size_t capacity,
                   std::vector<Event>& buffer,
                   const Event& event) {
  std::size_t current = writeIndex.load(std::memory_order_relaxed);
  std::size_t next = (current + 1) % capacity;
  if (next == readIndex.load(std::memory_order_acquire)) {
    return false; // Buffer is full
  }
  buffer[current] = event;
  writeIndex.store(next, std::memory_order_release);
  return true;
}

// Asynchronous task dispatcher
template <typename Func>
auto dispatchTask(Func task) -> std::future<decltype(task())> {
  return std::async(std::launch::async, task);
}

// Adaptive fault handler
inline int adaptiveFaultHandler(double errorRate, int currentMode) {
  const double HIGH_ERROR_THRESHOLD = 0.1;
  const double LOW_ERROR_THRESHOLD = 0.05;
  if (errorRate > HIGH_ERROR_THRESHOLD) {
    return currentMode + 1;
  } else if (errorRate < LOW_ERROR_THRESHOLD && currentMode > 0) {
    return currentMode - 1;
  }
  return currentMode;
}

int main() {
  const std::size_t frameSize = 8;
  float inputFrame[frameSize] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
  float outputFrame[frameSize] = {};
  processSignalFrame(inputFrame, outputFrame, frameSize);
  std::cout << "Processed Signal Frame: ";
  for (std::size_t i = 0; i < frameSize; ++i) { std::cout << outputFrame[i] << " "; }
  std::cout << std::endl;

  std::vector<int> circBuffer = {10, 20, 30, 40, 50};
  std::size_t head = 0;
  std::cout << "Circular Buffer Read: ";
  for (std::size_t i = 0; i < circBuffer.size() * 2; ++i) {
    std::cout << readCircularBuffer(circBuffer, head) << " ";
  }
  std::cout << std::endl;

  std::cout << "Periodic Task Execution:" << std::endl;
  schedulePeriodicTask([&]() {
    auto now = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();
    std::cout << "Task executed at: " << ms << " ms" << std::endl;
  }, std::chrono::milliseconds(500), 5);

  std::chrono::steady_clock::time_point scheduled = std::chrono::steady_clock::now() + std::chrono::milliseconds(100);
  std::cout << "Jitter Compensated Sleep:" << std::endl;
  for (int i = 0; i < 5; ++i) {
    jitterCompensatedSleep(scheduled, std::chrono::milliseconds(100));
    auto wakeTime = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(wakeTime.time_since_epoch()).count();
    std::cout << "Woke up at: " << ms << " ms" << std::endl;
  }

  const std::size_t bufferCapacity = 10;
  std::vector<Event> eventBuffer(bufferCapacity);
  std::atomic<std::size_t> writeIdx(0), readIdx(0);
  std::cout << "Lock-Free Enqueue:" << std::endl;
  for (int i = 0; i < 12; ++i) {
    Event e = {i, "Event " + std::to_string(i)};
    if (lockFreeEnqueue(writeIdx, readIdx, bufferCapacity, eventBuffer, e)) {
      std::cout << "Enqueued: " << e.message << std::endl;
    } else {
      std::cout << "Buffer full, cannot enqueue: " << e.message << std::endl;
    }
  }

  auto asyncResult = dispatchTask([]() -> int {
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    return 42;
  });
  std::cout << "Asynchronous Task Result: " << asyncResult.get() << std::endl;

  int mode = 1;
  double errorRate = 0.12;
  mode = adaptiveFaultHandler(errorRate, mode);
  std::cout << "Adaptive Mode after error rate " << errorRate << ": " << mode << std::endl;

  return 0;
}
