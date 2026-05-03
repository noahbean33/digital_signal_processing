#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <numeric>

// Moving Average Filter: computes the average over a sliding window.
std::vector<double> movingAverageFilter(const std::vector<double>& input, std::size_t windowSize) {
    if (windowSize == 0 || input.size() < windowSize) {
        return {};
    }
    std::vector<double> output;
    output.reserve(input.size() - windowSize + 1);
    double sum = 0.0;
    for (std::size_t i = 0; i < windowSize; ++i) {
        sum += input[i];
    }
    output.push_back(sum / windowSize);
    for (std::size_t i = windowSize; i < input.size(); ++i) {
        sum += input[i] - input[i - windowSize];
        output.push_back(sum / windowSize);
    }
    return output;
}

// Finite Difference: Approximates the derivative of the signal.
std::vector<double> computeFiniteDifference(const std::vector<double>& signal, double interval) {
    if (signal.size() < 2) {
        return {};
    }
    std::vector<double> derivative;
    derivative.reserve(signal.size() - 1);
    for (std::size_t i = 1; i < signal.size(); ++i) {
        derivative.push_back((signal[i] - signal[i - 1]) / interval);
    }
    return derivative;
}

// Zero Crossing Detection: Counts the number of sign changes in the signal.
int countZeroCrossings(const std::vector<double>& signal) {
    int zeroCrossings = 0;
    for (std::size_t i = 1; i < signal.size(); ++i) {
        if ((signal[i - 1] < 0 && signal[i] >= 0) || (signal[i - 1] > 0 && signal[i] <= 0)) {
            ++zeroCrossings;
        }
    }
    return zeroCrossings;
}

// Adaptive Peak Detection: Identifies local maxima that exceed an adaptive threshold.
std::vector<int> detectAdaptivePeaks(const std::vector<double>& signal, double thresholdMultiplier) {
    if (signal.size() < 3) {
        return {};
    }
    std::vector<int> peakIndices;
    // Compute global mean.
    double sum = 0.0;
    for (double val : signal) {
        sum += val;
    }
    double mean = sum / signal.size();
    // Compute standard deviation.
    double variance = 0.0;
    for (double val : signal) {
        variance += (val - mean) * (val - mean);
    }
    double stddev = std::sqrt(variance / signal.size());
    double adaptiveThreshold = mean + thresholdMultiplier * stddev;
    // Identify peaks based on local neighborhood and threshold.
    for (std::size_t i = 1; i < signal.size() - 1; ++i) {
        if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1] && signal[i] > adaptiveThreshold) {
            peakIndices.push_back(i);
        }
    }
    return peakIndices;
}

// Event Segmentation: Segments the signal based on an amplitude threshold.
std::vector<std::pair<int, int>> segmentEvents(const std::vector<double>& signal, double amplitudeThreshold) {
    std::vector<std::pair<int, int>> segments;
    bool segmentActive = false;
    int segmentStart = 0;
    for (std::size_t i = 0; i < signal.size(); ++i) {
        if (!segmentActive && std::abs(signal[i]) > amplitudeThreshold) {
            segmentActive = true;
            segmentStart = static_cast<int>(i);
        } else if (segmentActive && std::abs(signal[i]) <= amplitudeThreshold) {
            segments.push_back({segmentStart, static_cast<int>(i)});
            segmentActive = false;
        }
    }
    // Close any open segment at the end.
    if (segmentActive) {
        segments.push_back({segmentStart, static_cast<int>(signal.size() - 1)});
    }
    return segments;
}

// Robust Median Filter: Applies a sliding window median filter for impulse noise suppression.
std::vector<double> medianFilter(const std::vector<double>& signal, std::size_t windowSize) {
    if (windowSize % 2 == 0 || windowSize == 0 || signal.size() < windowSize) {
        return {};
    }
    std::vector<double> filtered;
    int halfWindow = windowSize / 2;
    filtered.reserve(signal.size());
    for (std::size_t i = 0; i < signal.size(); ++i) {
        std::vector<double> window;
        for (int j = -halfWindow; j <= halfWindow; ++j) {
            int index = std::min(std::max(static_cast<int>(i) + j, 0), static_cast<int>(signal.size() - 1));
            window.push_back(signal[index]);
        }
        std::sort(window.begin(), window.end());
        filtered.push_back(window[halfWindow]);
    }
    return filtered;
}

// Linear Interpolation: Estimates the signal at non-integer positions.
double linearInterpolation(const std::vector<double>& signal, double fractionalIndex) {
    if (signal.empty()) {
        return 0.0;
    }
    int index = static_cast<int>(fractionalIndex);
    double alpha = fractionalIndex - index;
    if (index < 0)
        return signal.front();
    if (index + 1 >= static_cast<int>(signal.size()))
        return signal.back();
    return (1.0 - alpha) * signal[index] + alpha * signal[index + 1];
}

// Synchronization Offset: Uses cross-correlation to estimate the delay between two signals.
int computeSynchronizationOffset(const std::vector<double>& stream1, const std::vector<double>& stream2, int maxDelay) {
    if (stream1.empty() || stream2.empty()) {
        return 0;
    }
    int bestOffset = 0;
    double maxCorrelation = -std::numeric_limits<double>::infinity();
    int length = static_cast<int>(stream1.size());
    for (int delay = -maxDelay; delay <= maxDelay; ++delay) {
        double correlation = 0.0;
        for (int i = 0; i < length; ++i) {
            int j = i + delay;
            if (j < 0 || j >= length)
                continue;
            correlation += stream1[i] * stream2[j];
        }
        if (correlation > maxCorrelation) {
            maxCorrelation = correlation;
            bestOffset = delay;
        }
    }
    return bestOffset;
}

int main() {
    // Example input signal for demonstration purposes.
    std::vector<double> signal = {0.1, 0.5, -0.2, 0.6, 1.2, 0.9, -0.4, -0.8, 0.3, 0.7, -0.1, 0.0};

    // 1. Apply Moving Average Filter with window size 3.
    std::vector<double> smoothed = movingAverageFilter(signal, 3);
    std::cout << "Smoothed Signal (Moving Average):" << std::endl;
    for (double val : smoothed) {
        std::cout << val << " ";
    }
    std::cout << std::endl << std::endl;

    // 2. Compute Finite Difference (Derivative) with interval = 1.0.
    std::vector<double> derivative = computeFiniteDifference(signal, 1.0);
    std::cout << "Finite Difference (Derivative):" << std::endl;
    for (double val : derivative) {
        std::cout << val << " ";
    }
    std::cout << std::endl << std::endl;

    // 3. Count the number of Zero Crossings in the signal.
    int zeroCross = countZeroCrossings(signal);
    std::cout << "Number of Zero Crossings: " << zeroCross << std::endl << std::endl;

    // 4. Detect Adaptive Peaks with a threshold multiplier of 1.0.
    std::vector<int> peaks = detectAdaptivePeaks(signal, 1.0);
    std::cout << "Adaptive Peak Indices:" << std::endl;
    for (int idx : peaks) {
        std::cout << idx << " ";
    }
    std::cout << std::endl << std::endl;

    // 5. Segment Events based on an amplitude threshold of 0.5.
    std::vector<std::pair<int, int>> segments = segmentEvents(signal, 0.5);
    std::cout << "Segmented Events (Start, End):" << std::endl;
    for (const auto& seg : segments) {
        std::cout << "(" << seg.first << ", " << seg.second << ") ";
    }
    std::cout << std::endl << std::endl;

    // 6. Apply Median Filter with window size 3.
    std::vector<double> medianFiltered = medianFilter(signal, 3);
    std::cout << "Median Filtered Signal:" << std::endl;
    for (double val : medianFiltered) {
        std::cout << val << " ";
    }
    std::cout << std::endl << std::endl;

    // 7. Estimate signal value at a fractional index (e.g., 5.5) using Linear Interpolation.
    double interpValue = linearInterpolation(signal, 5.5);
    std::cout << "Linear Interpolation at index 5.5: " << interpValue << std::endl << std::endl;

    // 8. Compute Synchronization Offset between two streams.
    std::vector<double> stream1 = signal;
    std::vector<double> stream2 = {0.0, 0.1, 0.5, -0.2, 0.6, 1.2, 0.9, -0.4, -0.8, 0.3, 0.7, -0.1};
    int syncOffset = computeSynchronizationOffset(stream1, stream2, 3);
    std::cout << "Synchronization Offset: " << syncOffset << std::endl;

    return 0;
}
