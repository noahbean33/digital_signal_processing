# Quick Start Guide

## Get Started in 5 Minutes

### Step 1: Build the Library

**Windows (with MinGW):**
```bash
build.bat
```

**Linux/macOS:**
```bash
make
```

### Step 2: Run Examples

**Windows:**
```bash
cd bin
example_dft.exe
example_convolution.exe
example_stats_utils.exe
```

**Linux/macOS:**
```bash
./bin/example_dft
./bin/example_convolution
./bin/example_stats_utils
```

### Step 3: View Results

Output signal files will be created in the `bin/` directory:
- `input_signal.dat` - Input time-domain signal
- `dft_magnitude.dat` - Frequency spectrum
- `filtered_signal.dat` - Filtered output
- And many more!

## Write Your First Program

Create `my_program.c`:

```c
#include "sigproc.h"

int main() {
    Signal sig, filtered;
    sig.data = NULL;
    filtered.data = NULL;
    
    /* Generate a noisy sine wave */
    Signal clean, noise;
    clean.data = NULL;
    noise.data = NULL;
    
    generate_sine(&clean, 200, 1.0, 50.0, 1000.0, 0.0);
    generate_noise(&noise, 200, 0.3);
    
    signal_alloc(&sig, 200);
    for (int i = 0; i < 200; i++) {
        sig.data[i] = clean.data[i] + noise.data[i];
    }
    
    /* Apply moving average filter */
    moving_average(&sig, &filtered, 11);
    
    /* Save results */
    signal_write("noisy.dat", &sig);
    signal_write("clean.dat", &filtered);
    
    /* Cleanup */
    signal_free(&sig);
    signal_free(&filtered);
    signal_free(&clean);
    signal_free(&noise);
    
    printf("Done! Check noisy.dat and clean.dat\n");
    return 0;
}
```

**Compile and run:**
```bash
# Windows
gcc -o my_program.exe my_program.c src/*.c -Iinclude

# Linux/macOS
gcc -o my_program my_program.c src/*.c -Iinclude -lm
```

## Common Tasks

### Task 1: Analyze Frequency Content
```c
Signal input, magnitude;
ComplexSignal dft_out;

input.data = NULL;
magnitude.data = NULL;
dft_out.real = NULL;
dft_out.imag = NULL;

signal_read("signal.dat", &input, 1000);
dft_real(&input, &dft_out);
compute_magnitude(&dft_out, &magnitude);
signal_write("spectrum.dat", &magnitude);

signal_free(&input);
signal_free(&magnitude);
complex_signal_free(&dft_out);
```

### Task 2: Filter a Signal
```c
Signal input, kernel, output;

input.data = NULL;
kernel.data = NULL;
output.data = NULL;

signal_read("noisy.dat", &input, 1000);

/* Create simple low-pass filter */
signal_alloc(&kernel, 11);
for (int i = 0; i < 11; i++) {
    kernel.data[i] = 1.0 / 11.0;
}

convolve(&input, &kernel, &output);
signal_write("filtered.dat", &output);

signal_free(&input);
signal_free(&kernel);
signal_free(&output);
```

### Task 3: Compute Statistics
```c
Signal sig;
sig.data = NULL;

signal_read("data.dat", &sig, 1000);

double mean = signal_mean(&sig);
double std_dev = signal_std_dev(&sig, 0);
double rms = signal_rms(&sig);

printf("Mean: %.4f, StdDev: %.4f, RMS: %.4f\n", 
       mean, std_dev, rms);

signal_free(&sig);
```

### Task 4: Generate Test Signals
```c
Signal sine, cosine, impulse, noise;

sine.data = NULL;
cosine.data = NULL;
impulse.data = NULL;
noise.data = NULL;

generate_sine(&sine, 100, 1.0, 10.0, 1000.0, 0.0);
generate_cosine(&cosine, 100, 1.0, 10.0, 1000.0, 0.0);
generate_impulse(&impulse, 100, 50, 1.0);
generate_noise(&noise, 100, 0.5);

signal_write("sine.dat", &sine);
signal_write("cosine.dat", &cosine);
signal_write("impulse.dat", &impulse);
signal_write("noise.dat", &noise);

signal_free(&sine);
signal_free(&cosine);
signal_free(&impulse);
signal_free(&noise);
```

## Available Functions

### DFT Operations
- `dft_real()` - DFT of real signal
- `dft_complex()` - DFT of complex signal
- `idft()` - Inverse DFT
- `idft_complex()` - Complex IDFT
- `compute_magnitude()` - Magnitude spectrum
- `compute_phase()` - Phase spectrum
- `compute_power()` - Power spectrum

### Convolution
- `convolve()` - Linear convolution
- `correlate()` - Cross-correlation
- `autocorrelate()` - Autocorrelation
- `running_sum()` - Cumulative sum
- `moving_average()` - Smoothing filter

### Statistics
- `signal_mean()` - Average
- `signal_variance()` - Variance
- `signal_std_dev()` - Standard deviation
- `signal_rms()` - RMS value
- `signal_min()`, `signal_max()` - Extremes
- `signal_energy()` - Total energy
- `signal_power()` - Average power

### Utilities
- `rect_to_polar()` - Convert to polar form
- `polar_to_rect()` - Convert to rectangular
- `signal_normalize()` - Scale to [-1, 1]
- `signal_scale()` - Multiply by factor
- `signal_add_offset()` - Add DC offset
- `generate_sine()`, `generate_cosine()` - Waveforms
- `generate_impulse()`, `generate_noise()` - Test signals

## Troubleshooting

### Build Errors
- **"gcc not found"**: Install MinGW (Windows) or GCC (Linux/macOS)
- **Math library errors**: Add `-lm` flag on Linux/macOS

### Runtime Errors
- **"Unable to open file"**: Check file path exists
- **"Unable to allocate memory"**: Signal too large or insufficient RAM
- **Segmentation fault**: Initialize pointers to NULL, call `_free()` functions

### Data Issues
- **Wrong results**: Check signal length matches file data
- **File format**: Ensure one value per line in .dat files

## Next Steps

1. **Read the full README.md** for detailed API documentation
2. **Study the examples/** to see complete programs
3. **Experiment** with different parameters and signal types
4. **Visualize** output .dat files with plotting tools (MATLAB, Python, gnuplot)

## Need Help?

- Check `README.md` for full documentation
- Look at `examples/` for working code
- Review the original code in `../signal_processing/` folder

Happy signal processing! 📊
