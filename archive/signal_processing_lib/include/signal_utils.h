/**
 * @file signal_utils.h
 * @brief Signal utility functions
 * 
 * This module provides utility functions for signal processing including
 * coordinate conversions, normalization, and signal generation.
 */

#ifndef SIGNAL_UTILS_H
#define SIGNAL_UTILS_H

#include "signal_io.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Convert rectangular to polar form
 * 
 * Converts complex signal from rectangular (real, imag) to polar (mag, phase):
 * Magnitude = sqrt(real² + imag²)
 * Phase = atan2(imag, real)
 * 
 * @param rectangular Input complex signal (rectangular form)
 * @param magnitude Output magnitude signal
 * @param phase Output phase signal (in radians)
 * @return 0 on success, -1 on failure
 */
int rect_to_polar(const ComplexSignal *rectangular, 
                  Signal *magnitude, Signal *phase);

/**
 * @brief Convert polar to rectangular form
 * 
 * Converts from polar (mag, phase) to rectangular (real, imag):
 * Real = magnitude * cos(phase)
 * Imag = magnitude * sin(phase)
 * 
 * @param magnitude Input magnitude signal
 * @param phase Input phase signal (in radians)
 * @param rectangular Output complex signal (rectangular form)
 * @return 0 on success, -1 on failure
 */
int polar_to_rect(const Signal *magnitude, const Signal *phase,
                  ComplexSignal *rectangular);

/**
 * @brief Normalize signal to range [-1, 1]
 * 
 * Scales signal so that maximum absolute value is 1.
 * 
 * @param signal Input/output signal (modified in place)
 * @return 0 on success, -1 on failure
 */
int signal_normalize(Signal *signal);

/**
 * @brief Scale signal by a constant factor
 * 
 * Multiplies all samples by the given factor.
 * 
 * @param signal Input/output signal (modified in place)
 * @param factor Scaling factor
 * @return 0 on success, -1 on failure
 */
int signal_scale(Signal *signal, double factor);

/**
 * @brief Add DC offset to signal
 * 
 * Adds a constant value to all samples.
 * 
 * @param signal Input/output signal (modified in place)
 * @param offset DC offset value
 * @return 0 on success, -1 on failure
 */
int signal_add_offset(Signal *signal, double offset);

/**
 * @brief Generate sine wave
 * 
 * Creates a sine wave signal:
 * x[n] = amplitude * sin(2π * frequency * n / sample_rate + phase)
 * 
 * @param output Output signal
 * @param length Number of samples
 * @param amplitude Amplitude of sine wave
 * @param frequency Frequency in Hz
 * @param sample_rate Sample rate in Hz
 * @param phase Initial phase in radians
 * @return 0 on success, -1 on failure
 */
int generate_sine(Signal *output, int length, double amplitude, 
                  double frequency, double sample_rate, double phase);

/**
 * @brief Generate cosine wave
 * 
 * Creates a cosine wave signal.
 * 
 * @param output Output signal
 * @param length Number of samples
 * @param amplitude Amplitude of cosine wave
 * @param frequency Frequency in Hz
 * @param sample_rate Sample rate in Hz
 * @param phase Initial phase in radians
 * @return 0 on success, -1 on failure
 */
int generate_cosine(Signal *output, int length, double amplitude,
                    double frequency, double sample_rate, double phase);

/**
 * @brief Generate impulse signal
 * 
 * Creates an impulse (delta function) at specified position.
 * 
 * @param output Output signal
 * @param length Number of samples
 * @param position Position of impulse
 * @param amplitude Amplitude of impulse
 * @return 0 on success, -1 on failure
 */
int generate_impulse(Signal *output, int length, int position, double amplitude);

/**
 * @brief Generate white noise
 * 
 * Creates random white noise with uniform distribution.
 * 
 * @param output Output signal
 * @param length Number of samples
 * @param amplitude Maximum amplitude (±amplitude)
 * @return 0 on success, -1 on failure
 */
int generate_noise(Signal *output, int length, double amplitude);

#endif /* SIGNAL_UTILS_H */
