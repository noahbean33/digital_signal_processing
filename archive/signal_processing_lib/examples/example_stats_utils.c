/**
 * @file example_stats_utils.c
 * @brief Example demonstrating signal statistics and utilities
 * 
 * This example shows how to:
 * - Compute statistical properties of signals
 * - Convert between rectangular and polar forms
 * - Generate various test signals
 * - Normalize and scale signals
 * 
 * Compile:
 *   gcc -o example_stats example_stats_utils.c ../src/*.c -I../include -lm
 * 
 * Run:
 *   ./example_stats
 */

#include <stdio.h>
#include "../include/sigproc.h"

int main() {
    Signal signal, magnitude, phase;
    ComplexSignal complex_sig, rect_sig;
    
    /* Initialize structures */
    signal.data = NULL;
    magnitude.data = NULL;
    phase.data = NULL;
    complex_sig.real = NULL;
    complex_sig.imag = NULL;
    rect_sig.real = NULL;
    rect_sig.imag = NULL;
    
    printf("=== Signal Statistics and Utilities Example ===\n\n");
    
    /* Generate test signal */
    int length = 100;
    printf("Generating test signal (sine + noise)...\n");
    
    Signal sine, noise;
    sine.data = NULL;
    noise.data = NULL;
    
    generate_sine(&sine, length, 2.0, 10.0, 1000.0, 0.0);
    generate_noise(&noise, length, 0.5);
    
    if (signal_alloc(&signal, length) == 0) {
        for (int i = 0; i < length; i++) {
            signal.data[i] = sine.data[i] + noise.data[i];
        }
    }
    
    signal_free(&sine);
    signal_free(&noise);
    
    /* Compute statistics */
    printf("\n=== Signal Statistics ===\n");
    double mean = signal_mean(&signal);
    double variance = signal_variance(&signal, mean);
    double std_dev = signal_std_dev(&signal, mean);
    double rms = signal_rms(&signal);
    double energy = signal_energy(&signal);
    double power = signal_power(&signal);
    double min_val = signal_min(&signal);
    double max_val = signal_max(&signal);
    
    printf("Mean:              %.6f\n", mean);
    printf("Variance:          %.6f\n", variance);
    printf("Std Deviation:     %.6f\n", std_dev);
    printf("RMS:               %.6f\n", rms);
    printf("Energy:            %.6f\n", energy);
    printf("Power:             %.6f\n", power);
    printf("Min value:         %.6f\n", min_val);
    printf("Max value:         %.6f\n", max_val);
    
    /* Normalize signal */
    printf("\n=== Signal Normalization ===\n");
    Signal normalized;
    normalized.data = NULL;
    if (signal_alloc(&normalized, signal.length) == 0) {
        for (int i = 0; i < signal.length; i++) {
            normalized.data[i] = signal.data[i];
        }
        signal_normalize(&normalized);
        printf("Signal normalized to [-1, 1] range\n");
        printf("New min: %.6f, max: %.6f\n", 
               signal_min(&normalized), signal_max(&normalized));
        signal_write("normalized_signal.dat", &normalized);
    }
    
    /* Test rectangular to polar conversion */
    printf("\n=== Rectangular to Polar Conversion ===\n");
    printf("Creating complex signal...\n");
    
    if (complex_signal_alloc(&complex_sig, 50) == 0) {
        /* Create complex signal: real=cos, imag=sin (unit circle) */
        for (int i = 0; i < 50; i++) {
            double angle = 2.0 * M_PI * i / 50.0;
            complex_sig.real[i] = 2.0 * cos(angle);
            complex_sig.imag[i] = 2.0 * sin(angle);
        }
        
        /* Convert to polar */
        if (rect_to_polar(&complex_sig, &magnitude, &phase) == 0) {
            printf("Converted to polar form\n");
            signal_write("magnitude.dat", &magnitude);
            signal_write("phase.dat", &phase);
            
            printf("Magnitude statistics:\n");
            printf("  Mean: %.6f\n", signal_mean(&magnitude));
            printf("  Min:  %.6f\n", signal_min(&magnitude));
            printf("  Max:  %.6f\n", signal_max(&magnitude));
            
            /* Convert back to rectangular */
            printf("\nConverting back to rectangular...\n");
            if (polar_to_rect(&magnitude, &phase, &rect_sig) == 0) {
                printf("Conversion successful\n");
                
                /* Check error */
                double error = 0.0;
                for (int i = 0; i < rect_sig.length; i++) {
                    double diff_r = complex_sig.real[i] - rect_sig.real[i];
                    double diff_i = complex_sig.imag[i] - rect_sig.imag[i];
                    error += diff_r * diff_r + diff_i * diff_i;
                }
                error = sqrt(error / rect_sig.length);
                printf("Conversion error (RMS): %e\n", error);
            }
        }
    }
    
    /* Generate various test signals */
    printf("\n=== Generating Test Signals ===\n");
    
    Signal sine_wave, cosine_wave, impulse;
    sine_wave.data = NULL;
    cosine_wave.data = NULL;
    impulse.data = NULL;
    
    generate_sine(&sine_wave, 100, 1.0, 5.0, 100.0, 0.0);
    signal_write("sine_wave.dat", &sine_wave);
    printf("Sine wave generated\n");
    
    generate_cosine(&cosine_wave, 100, 1.0, 5.0, 100.0, 0.0);
    signal_write("cosine_wave.dat", &cosine_wave);
    printf("Cosine wave generated\n");
    
    generate_impulse(&impulse, 100, 50, 1.0);
    signal_write("impulse.dat", &impulse);
    printf("Impulse generated at position 50\n");
    
    /* Clean up */
    signal_free(&signal);
    signal_free(&normalized);
    signal_free(&magnitude);
    signal_free(&phase);
    signal_free(&sine_wave);
    signal_free(&cosine_wave);
    signal_free(&impulse);
    complex_signal_free(&complex_sig);
    complex_signal_free(&rect_sig);
    
    printf("\nStatistics and utilities example complete!\n");
    printf("Files created:\n");
    printf("  - normalized_signal.dat\n");
    printf("  - magnitude.dat, phase.dat\n");
    printf("  - sine_wave.dat, cosine_wave.dat\n");
    printf("  - impulse.dat\n");
    
    return 0;
}
