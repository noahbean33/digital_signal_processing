#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <queue>
#include <functional>
#include <algorithm>
#include <numeric>
#include <chrono>

std::vector<double> processDspSegment(const std::vector<double>& segment) {
  std::vector<double> result(segment.size());
  for (size_t i = 0; i < segment.size(); ++i) {
    result[i] = 2.0 * segment[i];
  }
  return result;
}

void updateSharedParameter(std::atomic<double>& sharedParam, double newValue) {
  sharedParam.store(newValue, std::memory_order_relaxed);
}

std::pair<size_t, size_t> computeTaskChunk(size_t totalTasks, size_t totalThreads, size_t threadId) {
  if (totalThreads == 0) {
      return {0, 0};
  }
  size_t chunkSize = totalTasks / totalThreads;
  size_t remainder = totalTasks % totalThreads;
  size_t start = threadId * chunkSize + std::min(threadId, remainder);
  size_t end = start + chunkSize + (threadId < remainder ? 1 : 0);
  return {start, end};
}

void waitForProcessing(std::condition_variable& cv, std::mutex& mtx, bool& readyFlag) {
  std::unique_lock<std::mutex> lock(mtx);
  cv.wait(lock, [&readyFlag] { return readyFlag; });
}

void submitTask(std::queue<std::function<void()>>& taskQueue, std::mutex& queueMutex, std::condition_variable& cv, std::function<void()> task) {
  {
    std::lock_guard<std::mutex> lock(queueMutex);
    taskQueue.push(task);
  }
  cv.notify_one();
}

void workerThread(std::queue<std::function<void()>>& taskQueue, std::mutex& queueMutex, std::condition_variable& cv, std::atomic<bool>& stopFlag) {
  while (!stopFlag.load(std::memory_order_acquire)) {
    std::function<void()> task;
    {
      std::unique_lock<std::mutex> lock(queueMutex);
      cv.wait(lock, [&] { return stopFlag.load() || !taskQueue.empty(); });
      if (stopFlag.load() && taskQueue.empty()) {
        return;
      }
      task = taskQueue.front();
      taskQueue.pop();
    }
    task();
  }
}

std::vector<double> vectorizedTransform(const std::vector<double>& inputData) {
  std::vector<double> transformedData(inputData.size());
  std::transform(inputData.begin(), inputData.end(), transformedData.begin(), [](double sample) { return sample * 2.0; });
  return transformedData;
}

int main() {
  const size_t signalSize = 100;
  std::vector<double> signal(signalSize);
  std::iota(signal.begin(), signal.end(), 0.0);

  std::atomic<double> sharedParam(0.0);

  const size_t totalThreads = 4;
  std::vector<std::thread> threads;
  std::vector<std::vector<double>> processedSegments(totalThreads);
  std::vector<std::vector<double>> segments(totalThreads);
  for (size_t i = 0; i < totalThreads; ++i) {
    auto indices = computeTaskChunk(signal.size(), totalThreads, i);
    segments[i] = std::vector<double>(signal.begin() + indices.first, signal.begin() + indices.second);
  }

  for (size_t i = 0; i < totalThreads; ++i) {
    threads.emplace_back([&, i]() {
      processedSegments[i] = processDspSegment(segments[i]);
      double sum = std::accumulate(processedSegments[i].begin(), processedSegments[i].end(), 0.0);
      if (!processedSegments[i].empty()) {
          double average = sum / processedSegments[i].size();
          updateSharedParameter(sharedParam, average);
      }
    });
  }
  for (auto& t : threads) {
    t.join();
  }
  
  std::cout << "Updated Shared Parameter (from one segment average): " << sharedParam.load() << std::endl;

  std::queue<std::function<void()>> taskQueue;
  std::mutex queueMutex;
  std::condition_variable cv;
  std::atomic<bool> stopFlag(false);
  const size_t poolSize = 4;
  std::vector<std::thread> threadPool;
  for (size_t i = 0; i < poolSize; ++i) {
    threadPool.emplace_back(workerThread, std::ref(taskQueue), std::ref(queueMutex), std::ref(cv), std::ref(stopFlag));
  }

  for (size_t i = 0; i < totalThreads; ++i) {
    submitTask(taskQueue, queueMutex, cv, [=, &segments]() {
      if (!segments[i].empty()) {
          auto transformed = vectorizedTransform(segments[i]);
          std::cout << "Thread Pool Task " << i << " - First sample after vectorized transform: " << transformed.front() << std::endl;
      }
    });
  }

  std::this_thread::sleep_for(std::chrono::seconds(1));
  
  stopFlag.store(true);
  cv.notify_all();
  for (auto& worker : threadPool) {
    worker.join();
  }
  
  return 0;
}
