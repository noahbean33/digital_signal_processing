/**
 * @file example_dft.c
 * @brief Example demonstrating DFT and IDFT operations
 * 
 * This example shows how to:
 * - Generate a test signal (sum of two sine waves)
 * - Compute the DFT
 * - Calculate magnitude and phase spectra
 * - Perform IDFT to recover the original signal
 * 
 * Compile:
 *   gcc -o example_dft example_dft.c ../src/*.c -I../include -lm
 * 
 * Run:
 *   ./example_dft
 */

#include <stdio.h>
#include "../include/sigproc.h"

int main() {
    Signal input, output, magnitude, phase;
    ComplexSignal dft_output;
    
    /* Initialize structures */
    input.data = NULL;
    output.data = NULL;
    magnitude.data = NULL;
    phase.data = NULL;
    dft_output.real = NULL;
    dft_output.imag = NULL;
    
    printf("=== DFT Example ===\n\n");
    
    /* Generate test signal: sum of 1kHz and 5kHz sine waves */
    Signal sig1, sig2;
    sig1.data = NULL;
    sig2.data = NULL;
    
    int length = 320;
    double sample_rate = 48000;  /* 48 kHz */
    
    printf("Generating test signal (1kHz + 5kHz)...\n");
    generate_sine(&sig1, length, 1.0, 1000.0, sample_rate, 0.0);
    generate_sine(&sig2, length, 0.5, 5000.0, sample_rate, 0.0);
    
    /* Combine signals */
    if (signal_alloc(&input, length) == 0) {
        for (int i = 0; i < length; i++) {
            input.data[i] = sig1.data[i] + sig2.data[i];
        }
        signal_write("input_signal.dat", &input);
        printf("Input signal saved to input_signal.dat\n\n");
    }
    
    signal_free(&sig1);
    signal_free(&sig2);
    
    /* Compute DFT */
    printf("Computing DFT...\n");
    if (dft_real(&input, &dft_output) == 0) {
        /* Save DFT output */
        Signal temp_real, temp_imag;
        temp_real.data = dft_output.real;
        temp_real.length = dft_output.length;
        temp_imag.data = dft_output.imag;
        temp_imag.length = dft_output.length;
        
        signal_write("dft_real.dat", &temp_real);
        signal_write("dft_imag.dat", &temp_imag);
        printf("DFT output saved to dft_real.dat and dft_imag.dat\n\n");
        
        /* Compute magnitude spectrum */
        printf("Computing magnitude spectrum...\n");
        if (compute_magnitude(&dft_output, &magnitude) == 0) {
            signal_write("dft_magnitude.dat", &magnitude);
            printf("Magnitude spectrum saved to dft_magnitude.dat\n\n");
        }
        
        /* Compute phase spectrum */
        printf("Computing phase spectrum...\n");
        if (compute_phase(&dft_output, &phase) == 0) {
            signal_write("dft_phase.dat", &phase);
            printf("Phase spectrum saved to dft_phase.dat\n\n");
        }
        
        /* Perform IDFT to recover signal */
        printf("Performing IDFT to recover signal...\n");
        if (idft(&dft_output, &output) == 0) {
            signal_write("idft_output.dat", &output);
            printf("IDFT output saved to idft_output.dat\n\n");
            
            /* Compare with original */
            double error = 0.0;
            for (int i = 0; i < input.length && i < output.length; i++) {
                double diff = input.data[i] - output.data[i];
                error += diff * diff;
            }
            error = sqrt(error / input.length);
            printf("Reconstruction error (RMS): %e\n", error);
        }
    }
    
    /* Clean up */
    signal_free(&input);
    signal_free(&output);
    signal_free(&magnitude);
    signal_free(&phase);
    complex_signal_free(&dft_output);
    
    printf("\nDFT example complete!\n");
    printf("Files created:\n");
    printf("  - input_signal.dat\n");
    printf("  - dft_real.dat, dft_imag.dat\n");
    printf("  - dft_magnitude.dat, dft_phase.dat\n");
    printf("  - idft_output.dat\n");
    
    return 0;
}
