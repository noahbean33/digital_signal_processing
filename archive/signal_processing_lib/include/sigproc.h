/**
 * @file sigproc.h
 * @brief Main header file for the Signal Processing Library
 * 
 * Include this single header to access all library functions.
 * 
 * Usage:
 *   #include "sigproc.h"
 */

#ifndef SIGPROC_H
#define SIGPROC_H

/* Include all module headers */
#include "signal_io.h"
#include "dft.h"
#include "signal_stats.h"
#include "signal_conv.h"
#include "signal_utils.h"

/* Library version */
#define SIGPROC_VERSION_MAJOR 1
#define SIGPROC_VERSION_MINOR 0
#define SIGPROC_VERSION_PATCH 0

#endif /* SIGPROC_H */
