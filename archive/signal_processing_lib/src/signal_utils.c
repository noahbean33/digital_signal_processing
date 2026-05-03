/**
 * @file signal_utils.c
 * @brief Implementation of signal utility functions
 */

#include "../include/signal_utils.h"
#include <stdlib.h>
#include <time.h>

/* Initialize random number generator */
static int rng_initialized = 0;

static void init_rng() {
    if (!rng_initialized) {
        srand((unsigned int)time(NULL));
        rng_initialized = 1;
    }
}

int rect_to_polar(const ComplexSignal *rectangular, 
                  Signal *magnitude, Signal *phase) {
    int k;
    double real_val, imag_val;
    
    /* Allocate outputs */
    if (signal_alloc(magnitude, rectangular->length) != 0) {
        return -1;
    }
    if (signal_alloc(phase, rectangular->length) != 0) {
        signal_free(magnitude);
        return -1;
    }
    
    /* Convert to polar */
    for (k = 0; k < rectangular->length; k++) {
        real_val = rectangular->real[k];
        imag_val = rectangular->imag[k];
        
        /* Calculate magnitude */
        magnitude->data[k] = sqrt(real_val * real_val + imag_val * imag_val);
        
        /* Calculate phase (handle zero case) */
        if (real_val == 0.0 && imag_val == 0.0) {
            phase->data[k] = 0.0;
        } else {
            phase->data[k] = atan2(imag_val, real_val);
        }
    }
    
    return 0;
}

int polar_to_rect(const Signal *magnitude, const Signal *phase,
                  ComplexSignal *rectangular) {
    int k;
    
    /* Check lengths match */
    if (magnitude->length != phase->length) {
        printf("Error: Magnitude and phase signals must have same length\n");
        return -1;
    }
    
    /* Allocate output */
    if (complex_signal_alloc(rectangular, magnitude->length) != 0) {
        return -1;
    }
    
    /* Convert to rectangular */
    for (k = 0; k < magnitude->length; k++) {
        rectangular->real[k] = magnitude->data[k] * cos(phase->data[k]);
        rectangular->imag[k] = magnitude->data[k] * sin(phase->data[k]);
    }
    
    return 0;
}

int signal_normalize(Signal *signal) {
    int i;
    double max_abs = 0.0;
    double abs_val;
    
    /* Find maximum absolute value */
    for (i = 0; i < signal->length; i++) {
        abs_val = fabs(signal->data[i]);
        if (abs_val > max_abs) {
            max_abs = abs_val;
        }
    }
    
    /* Avoid division by zero */
    if (max_abs == 0.0) {
        return 0;
    }
    
    /* Normalize */
    for (i = 0; i < signal->length; i++) {
        signal->data[i] /= max_abs;
    }
    
    return 0;
}

int signal_scale(Signal *signal, double factor) {
    int i;
    
    for (i = 0; i < signal->length; i++) {
        signal->data[i] *= factor;
    }
    
    return 0;
}

int signal_add_offset(Signal *signal, double offset) {
    int i;
    
    for (i = 0; i < signal->length; i++) {
        signal->data[i] += offset;
    }
    
    return 0;
}

int generate_sine(Signal *output, int length, double amplitude, 
                  double frequency, double sample_rate, double phase) {
    int i;
    
    /* Allocate output */
    if (signal_alloc(output, length) != 0) {
        return -1;
    }
    
    /* Generate sine wave */
    for (i = 0; i < length; i++) {
        output->data[i] = amplitude * sin(2.0 * M_PI * frequency * i / sample_rate + phase);
    }
    
    return 0;
}

int generate_cosine(Signal *output, int length, double amplitude,
                    double frequency, double sample_rate, double phase) {
    int i;
    
    /* Allocate output */
    if (signal_alloc(output, length) != 0) {
        return -1;
    }
    
    /* Generate cosine wave */
    for (i = 0; i < length; i++) {
        output->data[i] = amplitude * cos(2.0 * M_PI * frequency * i / sample_rate + phase);
    }
    
    return 0;
}

int generate_impulse(Signal *output, int length, int position, double amplitude) {
    int i;
    
    /* Allocate output */
    if (signal_alloc(output, length) != 0) {
        return -1;
    }
    
    /* Initialize to zero */
    for (i = 0; i < length; i++) {
        output->data[i] = 0.0;
    }
    
    /* Set impulse */
    if (position >= 0 && position < length) {
        output->data[position] = amplitude;
    }
    
    return 0;
}

int generate_noise(Signal *output, int length, double amplitude) {
    int i;
    
    init_rng();
    
    /* Allocate output */
    if (signal_alloc(output, length) != 0) {
        return -1;
    }
    
    /* Generate noise */
    for (i = 0; i < length; i++) {
        output->data[i] = amplitude * (2.0 * (double)rand() / RAND_MAX - 1.0);
    }
    
    return 0;
}
