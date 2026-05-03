/**
 * @file signal_io.c
 * @brief Implementation of signal I/O operations
 */

#include "../include/signal_io.h"
#include <string.h>

int signal_alloc(Signal *sig, int length) {
    sig->data = (double *)malloc(length * sizeof(double));
    if (sig->data == NULL) {
        printf("Error: Unable to allocate memory for signal\n");
        return -1;
    }
    sig->length = length;
    return 0;
}

void signal_free(Signal *sig) {
    if (sig->data != NULL) {
        free(sig->data);
        sig->data = NULL;
    }
    sig->length = 0;
}

int complex_signal_alloc(ComplexSignal *sig, int length) {
    sig->real = (double *)malloc(length * sizeof(double));
    sig->imag = (double *)malloc(length * sizeof(double));
    
    if (sig->real == NULL || sig->imag == NULL) {
        printf("Error: Unable to allocate memory for complex signal\n");
        if (sig->real) free(sig->real);
        if (sig->imag) free(sig->imag);
        return -1;
    }
    
    sig->length = length;
    return 0;
}

void complex_signal_free(ComplexSignal *sig) {
    if (sig->real != NULL) {
        free(sig->real);
        sig->real = NULL;
    }
    if (sig->imag != NULL) {
        free(sig->imag);
        sig->imag = NULL;
    }
    sig->length = 0;
}

int signal_read(const char *filename, Signal *sig, int max_length) {
    FILE *file;
    int count = 0;
    double value;
    
    file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error: Unable to open file '%s' for reading\n", filename);
        return -1;
    }
    
    /* Allocate memory */
    if (signal_alloc(sig, max_length) != 0) {
        fclose(file);
        return -1;
    }
    
    /* Read values */
    while (count < max_length && fscanf(file, "%lf", &value) == 1) {
        sig->data[count] = value;
        count++;
    }
    
    fclose(file);
    sig->length = count;
    
    printf("Read %d samples from %s\n", count, filename);
    return count;
}

int signal_write(const char *filename, const Signal *sig) {
    FILE *file;
    int i;
    
    file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error: Unable to open file '%s' for writing\n", filename);
        return -1;
    }
    
    /* Write values */
    for (i = 0; i < sig->length; i++) {
        fprintf(file, "%f\n", sig->data[i]);
    }
    
    fclose(file);
    printf("Wrote %d samples to %s\n", sig->length, filename);
    return 0;
}

int complex_signal_read(const char *filename_real, const char *filename_imag, 
                        ComplexSignal *sig, int max_length) {
    FILE *file_real, *file_imag;
    int count = 0;
    double real_val, imag_val;
    
    file_real = fopen(filename_real, "r");
    file_imag = fopen(filename_imag, "r");
    
    if (file_real == NULL || file_imag == NULL) {
        printf("Error: Unable to open input files\n");
        if (file_real) fclose(file_real);
        if (file_imag) fclose(file_imag);
        return -1;
    }
    
    /* Allocate memory */
    if (complex_signal_alloc(sig, max_length) != 0) {
        fclose(file_real);
        fclose(file_imag);
        return -1;
    }
    
    /* Read values */
    while (count < max_length && 
           fscanf(file_real, "%lf", &real_val) == 1 &&
           fscanf(file_imag, "%lf", &imag_val) == 1) {
        sig->real[count] = real_val;
        sig->imag[count] = imag_val;
        count++;
    }
    
    fclose(file_real);
    fclose(file_imag);
    sig->length = count;
    
    printf("Read %d complex samples\n", count);
    return count;
}

int complex_signal_write(const char *filename_real, const char *filename_imag, 
                         const ComplexSignal *sig) {
    FILE *file_real, *file_imag;
    int i;
    
    file_real = fopen(filename_real, "w");
    file_imag = fopen(filename_imag, "w");
    
    if (file_real == NULL || file_imag == NULL) {
        printf("Error: Unable to open output files\n");
        if (file_real) fclose(file_real);
        if (file_imag) fclose(file_imag);
        return -1;
    }
    
    /* Write values */
    for (i = 0; i < sig->length; i++) {
        fprintf(file_real, "%f\n", sig->real[i]);
        fprintf(file_imag, "%f\n", sig->imag[i]);
    }
    
    fclose(file_real);
    fclose(file_imag);
    
    printf("Wrote %d complex samples\n", sig->length);
    return 0;
}

int signal_copy(const Signal *src, Signal *dst) {
    if (dst->length != src->length) {
        printf("Error: Destination signal must have same length as source\n");
        return -1;
    }
    
    memcpy(dst->data, src->data, src->length * sizeof(double));
    return 0;
}
