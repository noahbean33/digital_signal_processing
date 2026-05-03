/**
 * @file signal_io.h
 * @brief Signal input/output operations
 * 
 * This module provides functions for reading and writing signal data
 * to/from text files (.dat format).
 */

#ifndef SIGNAL_IO_H
#define SIGNAL_IO_H

#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Signal data structure
 */
typedef struct {
    double *data;       /* Signal samples */
    int length;         /* Number of samples */
} Signal;

/**
 * @brief Complex signal data structure
 */
typedef struct {
    double *real;       /* Real part */
    double *imag;       /* Imaginary part */
    int length;         /* Number of samples */
} ComplexSignal;

/**
 * @brief Allocate memory for a signal
 * 
 * @param sig Pointer to Signal structure
 * @param length Number of samples
 * @return 0 on success, -1 on failure
 */
int signal_alloc(Signal *sig, int length);

/**
 * @brief Free signal memory
 * 
 * @param sig Pointer to Signal structure
 */
void signal_free(Signal *sig);

/**
 * @brief Allocate memory for a complex signal
 * 
 * @param sig Pointer to ComplexSignal structure
 * @param length Number of samples
 * @return 0 on success, -1 on failure
 */
int complex_signal_alloc(ComplexSignal *sig, int length);

/**
 * @brief Free complex signal memory
 * 
 * @param sig Pointer to ComplexSignal structure
 */
void complex_signal_free(ComplexSignal *sig);

/**
 * @brief Read signal from file
 * 
 * Reads signal data from a text file with one value per line.
 * 
 * @param filename Path to input file
 * @param sig Pointer to Signal structure
 * @param max_length Maximum number of samples to read
 * @return Number of samples read, or -1 on error
 */
int signal_read(const char *filename, Signal *sig, int max_length);

/**
 * @brief Write signal to file
 * 
 * Writes signal data to a text file with one value per line.
 * 
 * @param filename Path to output file
 * @param sig Pointer to Signal structure
 * @return 0 on success, -1 on failure
 */
int signal_write(const char *filename, const Signal *sig);

/**
 * @brief Read complex signal from files
 * 
 * Reads real and imaginary parts from separate files.
 * 
 * @param filename_real Path to real part file
 * @param filename_imag Path to imaginary part file
 * @param sig Pointer to ComplexSignal structure
 * @param max_length Maximum number of samples to read
 * @return Number of samples read, or -1 on error
 */
int complex_signal_read(const char *filename_real, const char *filename_imag, 
                        ComplexSignal *sig, int max_length);

/**
 * @brief Write complex signal to files
 * 
 * Writes real and imaginary parts to separate files.
 * 
 * @param filename_real Path to real part output file
 * @param filename_imag Path to imaginary part output file
 * @param sig Pointer to ComplexSignal structure
 * @return 0 on success, -1 on failure
 */
int complex_signal_write(const char *filename_real, const char *filename_imag, 
                         const ComplexSignal *sig);

/**
 * @brief Copy signal data
 * 
 * @param src Source signal
 * @param dst Destination signal (must be pre-allocated)
 * @return 0 on success, -1 on failure
 */
int signal_copy(const Signal *src, Signal *dst);

#endif /* SIGNAL_IO_H */
