# Signal Processing Library

A comprehensive C library for digital signal processing (DSP) operations. This library provides clean, modular implementations of fundamental DSP algorithms including DFT/IDFT, convolution, statistical analysis, and signal utilities.

## Features

### 📊 Signal I/O
- Read/write signals from/to text files (.dat format)
- Support for real and complex signals
- Dynamic memory management
- Signal copying and manipulation utilities

### 🔄 Fourier Transform (DFT/IDFT)
- **DFT** - Discrete Fourier Transform for real and complex signals
- **IDFT** - Inverse Discrete Fourier Transform
- **Magnitude spectrum** - Frequency domain amplitude
- **Phase spectrum** - Frequency domain phase
- **Power spectrum** - Frequency domain power

### 🔢 Signal Statistics
- **Mean** - Average value
- **Variance** - Spread of values
- **Standard deviation** - Square root of variance
- **RMS** - Root mean square value
- **Min/Max** - Extreme values
- **Energy** - Total signal energy
- **Power** - Average signal power

### ⚡ Convolution Operations
- **Convolution** - Linear convolution of two signals
- **Correlation** - Cross-correlation of signals
- **Autocorrelation** - Self-correlation
- **Running sum** - Cumulative integration
- **Moving average** - Smoothing filter

### 🛠️ Signal Utilities
- **Rectangular ↔ Polar** - Coordinate conversion
- **Signal normalization** - Scale to [-1, 1]
- **Signal scaling** - Multiply by constant
- **DC offset** - Add constant to signal
- **Signal generation** - Sine, cosine, impulse, noise

## Directory Structure

```
signal_processing_lib/
├── include/              # Header files
│   ├── signal_io.h      # Signal I/O
│   ├── dft.h            # DFT/IDFT operations
│   ├── signal_stats.h   # Statistical functions
│   ├── signal_conv.h    # Convolution operations
│   ├── signal_utils.h   # Utility functions
│   └── sigproc.h        # Master header
├── src/                 # Implementation files
│   ├── signal_io.c
│   ├── dft.c
│   ├── signal_stats.c
│   ├── signal_conv.c
│   └── signal_utils.c
├── examples/            # Example programs
│   ├── example_dft.c
│   ├── example_convolution.c
│   └── example_stats_utils.c
├── test_signals/        # Sample signal data
└── README.md            # This file
```

## Getting Started

### Prerequisites
- C compiler (GCC, Clang, or MSVC)
- Standard C library
- Math library (libm)

### Compilation

#### Compile a single example:
```bash
# Linux/macOS
gcc -o example_dft examples/example_dft.c src/*.c -Iinclude -lm

# Windows (MinGW)
gcc -o example_dft.exe examples/example_dft.c src/*.c -Iinclude
```

#### Compile all examples:
```bash
# Linux/macOS
make

# Windows
build.bat
```

### Running Examples

```bash
# Run DFT example
./example_dft

# Run convolution example
./example_conv

# Run statistics example
./example_stats
```

## Usage Examples

### Basic Signal I/O

```c
#include "sigproc.h"

int main() {
    Signal sig;
    sig.data = NULL;
    
    /* Generate a sine wave */
    generate_sine(&sig, 100, 1.0, 10.0, 1000.0, 0.0);
    
    /* Write to file */
    signal_write("output.dat", &sig);
    
    /* Clean up */
    signal_free(&sig);
    return 0;
}
```

### DFT Analysis

```c
#include "sigproc.h"

int main() {
    Signal input, magnitude;
    ComplexSignal dft_output;
    
    input.data = NULL;
    magnitude.data = NULL;
    dft_output.real = NULL;
    dft_output.imag = NULL;
    
    /* Read input signal */
    signal_read("signal.dat", &input, 1000);
    
    /* Compute DFT */
    dft_real(&input, &dft_output);
    
    /* Get magnitude spectrum */
    compute_magnitude(&dft_output, &magnitude);
    signal_write("spectrum.dat", &magnitude);
    
    /* Clean up */
    signal_free(&input);
    signal_free(&magnitude);
    complex_signal_free(&dft_output);
    
    return 0;
}
```

### Signal Filtering

```c
#include "sigproc.h"

int main() {
    Signal input, output, kernel;
    
    input.data = NULL;
    output.data = NULL;
    kernel.data = NULL;
    
    /* Read noisy signal */
    signal_read("noisy.dat", &input, 1000);
    
    /* Create low-pass filter */
    signal_alloc(&kernel, 11);
    for (int i = 0; i < 11; i++) {
        kernel.data[i] = 1.0 / 11.0;  /* Average filter */
    }
    
    /* Apply filter via convolution */
    convolve(&input, &kernel, &output);
    signal_write("filtered.dat", &output);
    
    /* Clean up */
    signal_free(&input);
    signal_free(&output);
    signal_free(&kernel);
    
    return 0;
}
```

### Statistical Analysis

```c
#include "sigproc.h"

int main() {
    Signal sig;
    sig.data = NULL;
    
    signal_read("data.dat", &sig, 1000);
    
    /* Compute statistics */
    double mean = signal_mean(&sig);
    double std_dev = signal_std_dev(&sig, mean);
    double rms = signal_rms(&sig);
    double energy = signal_energy(&sig);
    
    printf("Mean: %.6f\n", mean);
    printf("Std Dev: %.6f\n", std_dev);
    printf("RMS: %.6f\n", rms);
    printf("Energy: %.6f\n", energy);
    
    signal_free(&sig);
    return 0;
}
```

### Coordinate Conversion

```c
#include "sigproc.h"

int main() {
    ComplexSignal complex_sig;
    Signal magnitude, phase;
    
    complex_sig.real = NULL;
    complex_sig.imag = NULL;
    magnitude.data = NULL;
    phase.data = NULL;
    
    /* Read complex signal */
    complex_signal_read("real.dat", "imag.dat", &complex_sig, 1000);
    
    /* Convert to polar form */
    rect_to_polar(&complex_sig, &magnitude, &phase);
    
    signal_write("magnitude.dat", &magnitude);
    signal_write("phase.dat", &phase);
    
    /* Clean up */
    complex_signal_free(&complex_sig);
    signal_free(&magnitude);
    signal_free(&phase);
    
    return 0;
}
```

## API Reference

### Core Data Structures

#### Signal
```c
typedef struct {
    double *data;   /* Signal samples */
    int length;     /* Number of samples */
} Signal;
```

#### ComplexSignal
```c
typedef struct {
    double *real;   /* Real part */
    double *imag;   /* Imaginary part */
    int length;     /* Number of samples */
} ComplexSignal;
```

### Function Return Values
- **0** - Success
- **-1** - Failure (error message printed to stdout)

### Memory Management
- Always initialize pointers to NULL
- Call appropriate `_free()` functions to release memory
- Functions handle memory allocation internally

## File Format

### Signal Data Files (.dat)
- Text files with one value per line
- Real numbers (double precision)
- Example:
```
1.234
-0.567
2.890
...
```

## Performance Considerations

- **DFT complexity**: O(N²) where N is signal length
- **Convolution**: O(N×M) where N=signal length, M=kernel length
- **Memory**: All operations use dynamic allocation
- **Precision**: Double precision (64-bit floating point)

## Improvements from Original Code

This library provides significant improvements:

1. **✅ Eliminated redundancy** - Consolidated duplicate waveform data and I/O functions
2. **✅ Modular design** - Organized by functionality (DFT, stats, convolution, utils)
3. **✅ Comprehensive API** - Clean, documented interface
4. **✅ Dynamic memory** - No fixed-size arrays
5. **✅ Error handling** - Return codes and messages
6. **✅ Complete examples** - Working programs for each module
7. **✅ Better structure** - Header/source separation

## Consolidation Summary

### Original Structure (14 folders)
- DFT, IDFT, IDFT_ecg, ComplexDFT
- Convolution
- Signal_Mean, SignalVariance, Signal_Standard_deviation
- Rect_to_Polar
- Various waveform/signal data folders

### New Structure (5 modules)
- **signal_io** - Unified I/O for all signal types
- **dft** - All Fourier transform operations
- **signal_stats** - All statistical functions
- **signal_conv** - Convolution and related operations
- **signal_utils** - Utilities and signal generation

## Known Limitations

- DFT is O(N²); consider FFT for large signals
- No support for multi-dimensional signals
- File I/O is text-based (not binary)
- Limited to double precision

## Future Enhancements

Potential additions:
1. **FFT** - Fast Fourier Transform (O(N log N))
2. **Filter design** - FIR/IIR filter design functions
3. **Window functions** - Hamming, Hanning, Blackman, etc.
4. **Spectrograms** - Short-time Fourier Transform
5. **Resampling** - Interpolation and decimation
6. **Binary I/O** - Faster file operations

## References

- Oppenheim & Schafer, "Discrete-Time Signal Processing"
- Proakis & Manolakis, "Digital Signal Processing"
- Smith, "The Scientist and Engineer's Guide to Digital Signal Processing"

## License

This is educational code. Use freely for learning and non-commercial purposes.

---

**Note**: This library consolidates and modernizes the code from the `signal_processing/` folder, which has been preserved for reference.
