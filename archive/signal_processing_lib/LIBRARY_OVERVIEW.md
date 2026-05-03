# Signal Processing Library - Overview

## Purpose

This library is a cleaned-up, modularized version of the original `signal_processing/` codebase. It consolidates redundant code, improves structure, and provides a professional API for digital signal processing operations.

## Key Improvements

### 1. Eliminated Code Redundancy
**Before:** Every program had duplicate waveform data and separate implementations.
**After:** Centralized signal generation in `signal_utils.c`, unified I/O in `signal_io.c`.

### 2. Modular Organization
**Before:** 14+ separate folders, each with standalone programs
**After:** 5 functional modules with clear separation:
- `signal_io` - File I/O and memory management
- `dft` - All Fourier transform operations
- `signal_stats` - Statistical analysis functions
- `signal_conv` - Convolution and correlation
- `signal_utils` - Utilities and signal generation

### 3. Professional Structure
```
Before:                          After:
DFT/                            include/    (headers)
  DFT/                          src/        (implementations)
    main.c                      examples/   (sample programs)
    waveforms.c                 test_signals/ (test data)
    *.dat
```

### 4. Better Memory Management
**Before:** Fixed-size arrays (`double Output_REX[SIG_LENGTH/2]`)
**After:** Dynamic allocation with proper cleanup functions

### 5. Comprehensive Documentation
**Before:** Minimal comments
**After:** Complete documentation:
- Function headers with @brief, @param, @return
- Module-level documentation
- README with usage examples
- This overview document

## Module Details

### signal_io Module
**Purpose:** Core signal file handling and memory management
**Key Functions:**
- `signal_alloc()`, `signal_free()` - Memory management
- `signal_read()`, `signal_write()` - File I/O
- `complex_signal_alloc()`, `complex_signal_free()` - Complex signals
- `complex_signal_read()`, `complex_signal_write()` - Complex I/O

**Improvements:**
- Dynamic memory allocation
- Support for both real and complex signals
- Error checking and reporting
- No hardcoded signal lengths

### dft Module
**Purpose:** Discrete Fourier Transform operations
**Consolidates:** DFT, IDFT, IDFT_ecg, ComplexDFT (4 original folders)

**Key Functions:**
- `dft_real()` - DFT for real-valued signals
- `dft_complex()` - DFT for complex signals
- `idft()` - Inverse DFT (real output)
- `idft_complex()` - Inverse DFT (complex output)
- `compute_magnitude()` - Magnitude spectrum
- `compute_phase()` - Phase spectrum
- `compute_power()` - Power spectrum

**Improvements:**
- Unified DFT implementation for real and complex
- Separate functions for magnitude/phase calculation
- Proper normalization in IDFT
- Clean separation of concerns

### signal_stats Module
**Purpose:** Statistical analysis of signals
**Consolidates:** Signal_Mean, SignalVariance, Signal_Standard_deviation (3 original folders)

**Key Functions:**
- `signal_mean()` - Average value
- `signal_variance()` - Variance calculation
- `signal_std_dev()` - Standard deviation
- `signal_rms()` - RMS value
- `signal_min()`, `signal_max()` - Find extremes
- `signal_energy()`, `signal_power()` - Energy/power

**Improvements:**
- All statistics in one module
- Consistent interface
- Added RMS, energy, power functions
- Efficient single-pass algorithms

### signal_conv Module
**Purpose:** Convolution and related operations
**Consolidates:** Convolution, correlation operations (1 original folder + new functions)

**Key Functions:**
- `convolve()` - Linear convolution
- `correlate()` - Cross-correlation
- `autocorrelate()` - Autocorrelation
- `running_sum()` - Cumulative integration
- `moving_average()` - Smoothing filter

**Improvements:**
- Generic convolution implementation
- Added correlation functions
- Running sum and moving average
- Proper output length calculation

### signal_utils Module
**Purpose:** Utility functions and signal generation
**Consolidates:** Rect_to_Polar, waveform generation, test signals (3+ original folders)

**Key Functions:**
- `rect_to_polar()` - Rectangular to polar conversion
- `polar_to_rect()` - Polar to rectangular conversion
- `signal_normalize()` - Normalize to [-1, 1]
- `signal_scale()`, `signal_add_offset()` - Signal manipulation
- `generate_sine()`, `generate_cosine()` - Waveform generation
- `generate_impulse()`, `generate_noise()` - Test signals

**Improvements:**
- Proper atan2 handling in phase calculation
- Signal generation functions (no need for external data files)
- Normalization and scaling utilities
- Consistent interface

## Files Not Included

### Build Artifacts
- `*.exe`, `*.o` - Compiled binaries
- `bin/`, `obj/` - Build directories
- `*.cbp`, `*.layout` - IDE project files

### Redundant Data Files
- Multiple copies of waveform data in each folder
- Output `.dat` files from test runs

### Irrelevant Folders
- `_20hz_sig/`, `_ecg_320/`, `ecg/` - Just data files
- `waveforms/`, `impulse_response/` - Data moved to test_signals/

## Usage Philosophy

### Simple Operations
For single operations, use individual functions:
```c
Signal sig;
sig.data = NULL;

signal_read("input.dat", &sig, 1000);
double mean = signal_mean(&sig);
printf("Mean: %.6f\n", mean);
signal_free(&sig);
```

### Complex Pipelines
Chain operations for signal processing pipelines:
```c
Signal input, filtered, magnitude;
ComplexSignal dft_out;

signal_read("noisy.dat", &input, 1000);

/* Filter in time domain */
moving_average(&input, &filtered, 11);

/* Analyze in frequency domain */
dft_real(&filtered, &dft_out);
compute_magnitude(&dft_out, &magnitude);

signal_write("spectrum.dat", &magnitude);

/* Cleanup */
signal_free(&input);
signal_free(&filtered);
signal_free(&magnitude);
complex_signal_free(&dft_out);
```

## Design Decisions

### 1. Dynamic Memory Allocation
**Why:** Supports arbitrary signal lengths without recompilation
**Trade-off:** Slightly more complex memory management, but worth the flexibility

### 2. Separate Real and Complex Functions
**Why:** Clearer intent, better type safety, avoids confusion
**Trade-off:** More functions, but each is simpler and focused

### 3. Data Structure Approach
**Why:** Encapsulates signal properties (data + length)
**Trade-off:** Requires passing structs, but prevents length mismatches

### 4. File I/O Format
**Why:** Text files are human-readable and tool-independent
**Trade-off:** Slower than binary, but better for learning and debugging

### 5. No Global Variables
**Why:** Better encapsulation, reentrant code, thread-safe potential
**Trade-off:** More function parameters, but cleaner design

## Testing Strategy

### Example Programs
Three comprehensive examples demonstrate:
1. DFT/IDFT operations and spectral analysis
2. Convolution, filtering, and correlation
3. Statistics and utility functions

### Validation
- DFT/IDFT round-trip testing
- Known signal properties (sine wave statistics)
- Visual inspection of output files
- Compare with original implementations

## Performance Characteristics

### DFT: O(N²)
- Suitable for signals up to ~1000 samples
- Consider FFT for larger signals (future enhancement)

### Convolution: O(N×M)
- N = signal length, M = kernel length
- Efficient for short kernels (< 100 samples)

### Statistics: O(N)
- Linear time, single pass when possible
- Minimal overhead

### Memory: O(N)
- Linear memory usage
- All allocations are temporary and freed

## Future Enhancements

### High Priority
1. **FFT** - Fast Fourier Transform (O(N log N))
2. **Window functions** - Hamming, Hanning, Blackman
3. **Filter design** - FIR/IIR coefficient generation

### Medium Priority
4. **Binary I/O** - Faster file operations
5. **Resampling** - Interpolation and decimation
6. **Spectrogram** - Short-Time Fourier Transform

### Low Priority
7. **Multi-channel** - Support for stereo signals
8. **Fixed-point** - Integer arithmetic option
9. **SIMD** - Vectorized operations for speed

## Migration Guide

### From Original Code

**Old (DFT example):**
```c
#define SIG_LENGTH 320
extern double InputSignal_f32_1kHz_15kHz[320];
double Output_REX[SIG_LENGTH/2];
double Output_IMX[SIG_LENGTH/2];

calc_sig_dft((double *)&InputSignal_f32_1kHz_15kHz[0],
             (double *)&Output_REX[0],
             (double *)&Output_IMX[0],
             (int)SIG_LENGTH);
```

**New:**
```c
Signal input;
ComplexSignal dft_out;

input.data = NULL;
dft_out.real = NULL;
dft_out.imag = NULL;

generate_sine(&input, 320, 1.0, 1000.0, 48000.0, 0.0);
dft_real(&input, &dft_out);

signal_free(&input);
complex_signal_free(&dft_out);
```

### Benefits of Migration
- No hardcoded lengths
- Automatic memory management
- Better error handling
- Cleaner, more readable code

## Conclusion

This library represents a significant improvement over the original codebase while preserving all functionality. It serves as:
- **A learning resource** for DSP algorithms
- **A practical tool** for signal analysis
- **A foundation** for advanced DSP projects

The original `signal_processing/` folder remains available for reference and comparison.

## Credits

Based on original code from the digital signal processing notes repository. Refactored and enhanced for educational purposes.
