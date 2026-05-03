/**
 * @file dft.c
 * @brief Implementation of DFT operations
 */

#include "../include/dft.h"

int dft_real(const Signal *input, ComplexSignal *output) {
    int i, k;
    int N = input->length;
    int N_out = N / 2;
    double angle;
    
    /* Allocate output */
    if (complex_signal_alloc(output, N_out) != 0) {
        return -1;
    }
    
    /* Initialize output to zero */
    for (k = 0; k < N_out; k++) {
        output->real[k] = 0.0;
        output->imag[k] = 0.0;
    }
    
    /* Compute DFT */
    for (k = 0; k < N_out; k++) {
        for (i = 0; i < N; i++) {
            angle = 2.0 * M_PI * k * i / N;
            output->real[k] += input->data[i] * cos(angle);
            output->imag[k] -= input->data[i] * sin(angle);
        }
    }
    
    return 0;
}

int dft_complex(const ComplexSignal *input, ComplexSignal *output) {
    int i, k;
    int N = input->length;
    double angle, cos_val, sin_val;
    double temp_real, temp_imag;
    
    /* Allocate output */
    if (complex_signal_alloc(output, N) != 0) {
        return -1;
    }
    
    /* Initialize output to zero */
    for (k = 0; k < N; k++) {
        output->real[k] = 0.0;
        output->imag[k] = 0.0;
    }
    
    /* Compute DFT */
    for (k = 0; k < N; k++) {
        for (i = 0; i < N; i++) {
            angle = 2.0 * M_PI * k * i / N;
            cos_val = cos(angle);
            sin_val = sin(angle);
            
            /* Complex multiplication: (real + j*imag) * (cos - j*sin) */
            temp_real = input->real[i] * cos_val + input->imag[i] * sin_val;
            temp_imag = input->imag[i] * cos_val - input->real[i] * sin_val;
            
            output->real[k] += temp_real;
            output->imag[k] += temp_imag;
        }
    }
    
    return 0;
}

int idft(const ComplexSignal *input, Signal *output) {
    int i, k;
    int N_in = input->length;
    int N = N_in * 2;  /* Reconstruct full signal length */
    double angle;
    
    /* Allocate output */
    if (signal_alloc(output, N) != 0) {
        return -1;
    }
    
    /* Initialize output to zero */
    for (i = 0; i < N; i++) {
        output->data[i] = 0.0;
    }
    
    /* Compute IDFT */
    for (i = 0; i < N; i++) {
        for (k = 0; k < N_in; k++) {
            angle = 2.0 * M_PI * k * i / N;
            /* Real part only (assuming conjugate symmetry) */
            output->data[i] += input->real[k] * cos(angle) - input->imag[k] * sin(angle);
        }
        /* Normalize by N */
        output->data[i] /= N;
    }
    
    return 0;
}

int idft_complex(const ComplexSignal *input, ComplexSignal *output) {
    int i, k;
    int N = input->length;
    double angle, cos_val, sin_val;
    double temp_real, temp_imag;
    
    /* Allocate output */
    if (complex_signal_alloc(output, N) != 0) {
        return -1;
    }
    
    /* Initialize output to zero */
    for (i = 0; i < N; i++) {
        output->real[i] = 0.0;
        output->imag[i] = 0.0;
    }
    
    /* Compute IDFT */
    for (i = 0; i < N; i++) {
        for (k = 0; k < N; k++) {
            angle = 2.0 * M_PI * k * i / N;
            cos_val = cos(angle);
            sin_val = sin(angle);
            
            /* Complex multiplication: (real + j*imag) * (cos + j*sin) */
            temp_real = input->real[k] * cos_val - input->imag[k] * sin_val;
            temp_imag = input->imag[k] * cos_val + input->real[k] * sin_val;
            
            output->real[i] += temp_real;
            output->imag[i] += temp_imag;
        }
        /* Normalize by N */
        output->real[i] /= N;
        output->imag[i] /= N;
    }
    
    return 0;
}

int compute_magnitude(const ComplexSignal *complex_sig, Signal *magnitude) {
    int k;
    
    /* Allocate output */
    if (signal_alloc(magnitude, complex_sig->length) != 0) {
        return -1;
    }
    
    /* Calculate magnitude */
    for (k = 0; k < complex_sig->length; k++) {
        magnitude->data[k] = sqrt(complex_sig->real[k] * complex_sig->real[k] + 
                                  complex_sig->imag[k] * complex_sig->imag[k]);
    }
    
    return 0;
}

int compute_phase(const ComplexSignal *complex_sig, Signal *phase) {
    int k;
    double real_val, imag_val;
    
    /* Allocate output */
    if (signal_alloc(phase, complex_sig->length) != 0) {
        return -1;
    }
    
    /* Calculate phase */
    for (k = 0; k < complex_sig->length; k++) {
        real_val = complex_sig->real[k];
        imag_val = complex_sig->imag[k];
        
        /* Handle zero case */
        if (real_val == 0.0 && imag_val == 0.0) {
            phase->data[k] = 0.0;
        } else {
            phase->data[k] = atan2(imag_val, real_val);
        }
    }
    
    return 0;
}

int compute_power(const ComplexSignal *complex_sig, Signal *power) {
    int k;
    
    /* Allocate output */
    if (signal_alloc(power, complex_sig->length) != 0) {
        return -1;
    }
    
    /* Calculate power */
    for (k = 0; k < complex_sig->length; k++) {
        power->data[k] = complex_sig->real[k] * complex_sig->real[k] + 
                         complex_sig->imag[k] * complex_sig->imag[k];
    }
    
    return 0;
}
