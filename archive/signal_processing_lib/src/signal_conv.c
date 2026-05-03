/**
 * @file signal_conv.c
 * @brief Implementation of convolution operations
 */

#include "../include/signal_conv.h"

int convolve(const Signal *input, const Signal *impulse, Signal *output) {
    int i, j;
    int output_length = input->length + impulse->length - 1;
    
    /* Allocate output */
    if (signal_alloc(output, output_length) != 0) {
        return -1;
    }
    
    /* Initialize output to zero */
    for (i = 0; i < output_length; i++) {
        output->data[i] = 0.0;
    }
    
    /* Perform convolution */
    for (i = 0; i < input->length; i++) {
        for (j = 0; j < impulse->length; j++) {
            output->data[i + j] += input->data[i] * impulse->data[j];
        }
    }
    
    return 0;
}

int correlate(const Signal *sig1, const Signal *sig2, Signal *output) {
    int i, j, n;
    int max_len = (sig1->length > sig2->length) ? sig1->length : sig2->length;
    int output_length = sig1->length + sig2->length - 1;
    
    /* Allocate output */
    if (signal_alloc(output, output_length) != 0) {
        return -1;
    }
    
    /* Initialize output to zero */
    for (i = 0; i < output_length; i++) {
        output->data[i] = 0.0;
    }
    
    /* Perform correlation */
    for (n = 0; n < output_length; n++) {
        int shift = n - sig1->length + 1;
        for (i = 0; i < sig1->length; i++) {
            j = i + shift;
            if (j >= 0 && j < sig2->length) {
                output->data[n] += sig1->data[i] * sig2->data[j];
            }
        }
    }
    
    return 0;
}

int autocorrelate(const Signal *input, Signal *output) {
    /* Autocorrelation is just correlation with itself */
    return correlate(input, input, output);
}

int running_sum(const Signal *input, Signal *output) {
    int i;
    double sum = 0.0;
    
    /* Allocate output */
    if (signal_alloc(output, input->length) != 0) {
        return -1;
    }
    
    /* Compute running sum */
    for (i = 0; i < input->length; i++) {
        sum += input->data[i];
        output->data[i] = sum;
    }
    
    return 0;
}

int moving_average(const Signal *input, Signal *output, int window_size) {
    int i, j;
    double sum;
    int count;
    
    /* Allocate output */
    if (signal_alloc(output, input->length) != 0) {
        return -1;
    }
    
    /* Compute moving average */
    for (i = 0; i < input->length; i++) {
        sum = 0.0;
        count = 0;
        
        /* Average over window */
        for (j = 0; j < window_size; j++) {
            int idx = i - j;
            if (idx >= 0) {
                sum += input->data[idx];
                count++;
            }
        }
        
        output->data[i] = sum / count;
    }
    
    return 0;
}
