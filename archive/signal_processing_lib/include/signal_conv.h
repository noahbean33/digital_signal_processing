/**
 * @file signal_conv.h
 * @brief Signal convolution operations
 * 
 * This module provides functions for convolving signals, which is
 * fundamental to filtering and system analysis in signal processing.
 */

#ifndef SIGNAL_CONV_H
#define SIGNAL_CONV_H

#include "signal_io.h"

/**
 * @brief Perform convolution of two signals
 * 
 * Computes the convolution of input signal with impulse response:
 * y[n] = Σ x[k] * h[n-k]
 * 
 * Output length = input_length + impulse_length - 1
 * 
 * @param input Input signal
 * @param impulse Impulse response (filter kernel)
 * @param output Output signal (must be pre-allocated with correct length)
 * @return 0 on success, -1 on failure
 */
int convolve(const Signal *input, const Signal *impulse, Signal *output);

/**
 * @brief Perform correlation of two signals
 * 
 * Computes the cross-correlation of two signals:
 * r[n] = Σ x[k] * y[k+n]
 * 
 * @param sig1 First signal
 * @param sig2 Second signal
 * @param output Output correlation signal
 * @return 0 on success, -1 on failure
 */
int correlate(const Signal *sig1, const Signal *sig2, Signal *output);

/**
 * @brief Perform autocorrelation of a signal
 * 
 * Computes the autocorrelation (correlation with itself):
 * r[n] = Σ x[k] * x[k+n]
 * 
 * @param input Input signal
 * @param output Output autocorrelation signal
 * @return 0 on success, -1 on failure
 */
int autocorrelate(const Signal *input, Signal *output);

/**
 * @brief Running sum filter (integrator)
 * 
 * Computes cumulative sum:
 * y[n] = Σ (k=0 to n) x[k]
 * 
 * @param input Input signal
 * @param output Output signal
 * @return 0 on success, -1 on failure
 */
int running_sum(const Signal *input, Signal *output);

/**
 * @brief Moving average filter
 * 
 * Computes moving average with specified window size:
 * y[n] = (1/M) * Σ x[n-k] for k=0 to M-1
 * 
 * @param input Input signal
 * @param output Output signal
 * @param window_size Size of averaging window
 * @return 0 on success, -1 on failure
 */
int moving_average(const Signal *input, Signal *output, int window_size);

#endif /* SIGNAL_CONV_H */
