/**
 * @file dft.h
 * @brief Discrete Fourier Transform operations
 * 
 * This module provides functions for performing DFT and IDFT on real
 * and complex signals. The DFT converts time-domain signals to frequency domain.
 */

#ifndef DFT_H
#define DFT_H

#include "signal_io.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Compute DFT of a real signal
 * 
 * Computes the Discrete Fourier Transform of a real-valued signal.
 * Output contains N/2 frequency bins (DC to Nyquist).
 * 
 * Formula: X[k] = Σ x[n] * exp(-j*2π*k*n/N)
 * 
 * @param input Input signal (time domain)
 * @param output Output complex signal (frequency domain)
 * @return 0 on success, -1 on failure
 */
int dft_real(const Signal *input, ComplexSignal *output);

/**
 * @brief Compute DFT of a complex signal
 * 
 * Computes the Discrete Fourier Transform of a complex-valued signal.
 * 
 * @param input Input complex signal (time domain)
 * @param output Output complex signal (frequency domain)
 * @return 0 on success, -1 on failure
 */
int dft_complex(const ComplexSignal *input, ComplexSignal *output);

/**
 * @brief Compute Inverse DFT (IDFT)
 * 
 * Converts frequency domain signal back to time domain.
 * 
 * Formula: x[n] = (1/N) * Σ X[k] * exp(j*2π*k*n/N)
 * 
 * @param input Input complex signal (frequency domain)
 * @param output Output signal (time domain, real part only)
 * @return 0 on success, -1 on failure
 */
int idft(const ComplexSignal *input, Signal *output);

/**
 * @brief Compute IDFT (complex output)
 * 
 * Converts frequency domain signal back to time domain complex signal.
 * 
 * @param input Input complex signal (frequency domain)
 * @param output Output complex signal (time domain)
 * @return 0 on success, -1 on failure
 */
int idft_complex(const ComplexSignal *input, ComplexSignal *output);

/**
 * @brief Compute magnitude spectrum
 * 
 * Calculates the magnitude of each frequency bin:
 * Magnitude[k] = sqrt(Real[k]² + Imag[k]²)
 * 
 * @param complex_sig Input complex signal
 * @param magnitude Output magnitude signal
 * @return 0 on success, -1 on failure
 */
int compute_magnitude(const ComplexSignal *complex_sig, Signal *magnitude);

/**
 * @brief Compute phase spectrum
 * 
 * Calculates the phase angle of each frequency bin:
 * Phase[k] = atan2(Imag[k], Real[k])
 * 
 * @param complex_sig Input complex signal
 * @param phase Output phase signal (in radians)
 * @return 0 on success, -1 on failure
 */
int compute_phase(const ComplexSignal *complex_sig, Signal *phase);

/**
 * @brief Compute power spectrum
 * 
 * Calculates the power of each frequency bin:
 * Power[k] = Real[k]² + Imag[k]²
 * 
 * @param complex_sig Input complex signal
 * @param power Output power signal
 * @return 0 on success, -1 on failure
 */
int compute_power(const ComplexSignal *complex_sig, Signal *power);

#endif /* DFT_H */
