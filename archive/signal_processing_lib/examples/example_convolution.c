/**
 * @file example_convolution.c
 * @brief Example demonstrating convolution operations
 * 
 * This example shows how to:
 * - Perform convolution with a low-pass filter
 * - Apply moving average filter
 * - Compute autocorrelation
 * 
 * Compile:
 *   gcc -o example_conv example_convolution.c ../src/*.c -I../include -lm
 * 
 * Run:
 *   ./example_conv
 */

#include <stdio.h>
#include "../include/sigproc.h"

int main() {
    Signal input, impulse, output, autocorr;
    
    /* Initialize structures */
    input.data = NULL;
    impulse.data = NULL;
    output.data = NULL;
    autocorr.data = NULL;
    
    printf("=== Convolution Example ===\n\n");
    
    /* Generate noisy test signal */
    Signal clean, noise;
    clean.data = NULL;
    noise.data = NULL;
    
    int length = 200;
    double sample_rate = 1000;  /* 1 kHz */
    
    printf("Generating noisy signal...\n");
    generate_sine(&clean, length, 1.0, 50.0, sample_rate, 0.0);
    generate_noise(&noise, length, 0.3);
    
    if (signal_alloc(&input, length) == 0) {
        for (int i = 0; i < length; i++) {
            input.data[i] = clean.data[i] + noise.data[i];
        }
        signal_write("noisy_signal.dat", &input);
        printf("Noisy signal saved\n\n");
    }
    
    signal_free(&clean);
    signal_free(&noise);
    
    /* Create simple low-pass filter (moving average kernel) */
    printf("Creating low-pass filter kernel...\n");
    int kernel_length = 11;
    if (signal_alloc(&impulse, kernel_length) == 0) {
        for (int i = 0; i < kernel_length; i++) {
            impulse.data[i] = 1.0 / kernel_length;  /* Normalized */
        }
        signal_write("filter_kernel.dat", &impulse);
        printf("Filter kernel saved\n\n");
    }
    
    /* Perform convolution */
    printf("Applying filter via convolution...\n");
    if (convolve(&input, &impulse, &output) == 0) {
        signal_write("filtered_signal.dat", &output);
        printf("Filtered signal saved to filtered_signal.dat\n\n");
        
        /* Calculate signal statistics */
        double mean = signal_mean(&output);
        double std_dev = signal_std_dev(&output, mean);
        double rms = signal_rms(&output);
        
        printf("Output signal statistics:\n");
        printf("  Mean:   %.6f\n", mean);
        printf("  StdDev: %.6f\n", std_dev);
        printf("  RMS:    %.6f\n\n", rms);
    }
    
    /* Apply moving average filter */
    Signal smoothed;
    smoothed.data = NULL;
    printf("Applying moving average filter (window=15)...\n");
    if (moving_average(&input, &smoothed, 15) == 0) {
        signal_write("smoothed_signal.dat", &smoothed);
        printf("Smoothed signal saved\n\n");
    }
    
    /* Compute autocorrelation */
    printf("Computing autocorrelation...\n");
    if (autocorrelate(&input, &autocorr) == 0) {
        signal_write("autocorrelation.dat", &autocorr);
        printf("Autocorrelation saved\n\n");
    }
    
    /* Running sum */
    Signal integrated;
    integrated.data = NULL;
    printf("Computing running sum (integration)...\n");
    if (running_sum(&input, &integrated) == 0) {
        signal_write("running_sum.dat", &integrated);
        printf("Running sum saved\n\n");
    }
    
    /* Clean up */
    signal_free(&input);
    signal_free(&impulse);
    signal_free(&output);
    signal_free(&autocorr);
    signal_free(&smoothed);
    signal_free(&integrated);
    
    printf("Convolution example complete!\n");
    printf("Files created:\n");
    printf("  - noisy_signal.dat\n");
    printf("  - filter_kernel.dat\n");
    printf("  - filtered_signal.dat\n");
    printf("  - smoothed_signal.dat\n");
    printf("  - autocorrelation.dat\n");
    printf("  - running_sum.dat\n");
    
    return 0;
}
