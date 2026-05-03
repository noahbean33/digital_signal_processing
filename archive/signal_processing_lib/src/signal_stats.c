/**
 * @file signal_stats.c
 * @brief Implementation of signal statistical functions
 */

#include "../include/signal_stats.h"

double signal_mean(const Signal *sig) {
    int i;
    double sum = 0.0;
    
    for (i = 0; i < sig->length; i++) {
        sum += sig->data[i];
    }
    
    return sum / sig->length;
}

double signal_variance(const Signal *sig, double mean) {
    int i;
    double sum = 0.0;
    double diff;
    
    /* Compute mean if not provided */
    if (mean == 0.0) {
        mean = signal_mean(sig);
    }
    
    /* Compute variance */
    for (i = 0; i < sig->length; i++) {
        diff = sig->data[i] - mean;
        sum += diff * diff;
    }
    
    return sum / sig->length;
}

double signal_std_dev(const Signal *sig, double mean) {
    return sqrt(signal_variance(sig, mean));
}

double signal_rms(const Signal *sig) {
    int i;
    double sum = 0.0;
    
    for (i = 0; i < sig->length; i++) {
        sum += sig->data[i] * sig->data[i];
    }
    
    return sqrt(sum / sig->length);
}

double signal_min(const Signal *sig) {
    int i;
    double min_val = sig->data[0];
    
    for (i = 1; i < sig->length; i++) {
        if (sig->data[i] < min_val) {
            min_val = sig->data[i];
        }
    }
    
    return min_val;
}

double signal_max(const Signal *sig) {
    int i;
    double max_val = sig->data[0];
    
    for (i = 1; i < sig->length; i++) {
        if (sig->data[i] > max_val) {
            max_val = sig->data[i];
        }
    }
    
    return max_val;
}

double signal_energy(const Signal *sig) {
    int i;
    double energy = 0.0;
    
    for (i = 0; i < sig->length; i++) {
        energy += sig->data[i] * sig->data[i];
    }
    
    return energy;
}

double signal_power(const Signal *sig) {
    return signal_energy(sig) / sig->length;
}
