/**
 * @file signal_stats.h
 * @brief Signal statistical analysis functions
 * 
 * This module provides functions for computing statistical properties
 * of signals including mean, variance, standard deviation, and more.
 */

#ifndef SIGNAL_STATS_H
#define SIGNAL_STATS_H

#include "signal_io.h"
#include <math.h>

/**
 * @brief Compute signal mean (average)
 * 
 * Calculates the arithmetic mean of signal samples:
 * Mean = (1/N) * Σ x[n]
 * 
 * @param sig Input signal
 * @return Mean value
 */
double signal_mean(const Signal *sig);

/**
 * @brief Compute signal variance
 * 
 * Calculates the variance (average squared deviation from mean):
 * Var = (1/N) * Σ (x[n] - mean)²
 * 
 * @param sig Input signal
 * @param mean Signal mean (pass 0 to compute automatically)
 * @return Variance value
 */
double signal_variance(const Signal *sig, double mean);

/**
 * @brief Compute signal standard deviation
 * 
 * Calculates the standard deviation (square root of variance):
 * StdDev = sqrt(Var)
 * 
 * @param sig Input signal
 * @param mean Signal mean (pass 0 to compute automatically)
 * @return Standard deviation value
 */
double signal_std_dev(const Signal *sig, double mean);

/**
 * @brief Compute signal RMS (Root Mean Square)
 * 
 * Calculates the RMS value:
 * RMS = sqrt((1/N) * Σ x[n]²)
 * 
 * @param sig Input signal
 * @return RMS value
 */
double signal_rms(const Signal *sig);

/**
 * @brief Find signal minimum value
 * 
 * @param sig Input signal
 * @return Minimum value in signal
 */
double signal_min(const Signal *sig);

/**
 * @brief Find signal maximum value
 * 
 * @param sig Input signal
 * @return Maximum value in signal
 */
double signal_max(const Signal *sig);

/**
 * @brief Compute signal energy
 * 
 * Calculates the total energy:
 * Energy = Σ x[n]²
 * 
 * @param sig Input signal
 * @return Energy value
 */
double signal_energy(const Signal *sig);

/**
 * @brief Compute signal power
 * 
 * Calculates the average power:
 * Power = (1/N) * Σ x[n]²
 * 
 * @param sig Input signal
 * @return Power value
 */
double signal_power(const Signal *sig);

#endif /* SIGNAL_STATS_H */
