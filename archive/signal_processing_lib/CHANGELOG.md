# Changelog

## Version 1.0.0 - Initial Release

### Created
- **Core Modules**
  - `signal_io` - Signal I/O with dynamic memory management
  - `dft` - DFT, IDFT, and spectral analysis (magnitude, phase, power)
  - `signal_stats` - 8 statistical functions (mean, variance, std dev, RMS, min, max, energy, power)
  - `signal_conv` - Convolution, correlation, autocorrelation, running sum, moving average
  - `signal_utils` - Coordinate conversion, normalization, signal generation

### Features
- Support for real and complex signals
- Dynamic memory allocation (no fixed sizes)
- Comprehensive API documentation with doxygen-style comments
- Error handling with return codes and messages
- 3 complete example programs demonstrating library usage

### Documentation
- README.md - Comprehensive user guide with examples
- QUICKSTART.md - 5-minute getting started guide
- LIBRARY_OVERVIEW.md - Technical design documentation
- CHANGELOG.md - This file
- Inline code documentation

### Build System
- Makefile for Linux/macOS/MinGW
- build.bat for Windows
- Example compilation instructions

### Test Data
- test_signals/ directory for sample data files
- Examples generate their own test signals

### Improvements Over Original Code
1. ✅ Eliminated duplicate waveform data across folders
2. ✅ Consolidated 14+ folders into 5 modules
3. ✅ Replaced fixed-size arrays with dynamic allocation
4. ✅ Added comprehensive documentation
5. ✅ Removed IDE-specific files (*.cbp, *.layout, bin/, obj/)
6. ✅ Created professional library structure
7. ✅ Unified signal I/O functions
8. ✅ Added signal generation utilities
9. ✅ Provided working example programs
10. ✅ Created build automation scripts

### Files Consolidated

**DFT Operations** (4 folders → 1 module):
- DFT
- IDFT
- IDFT_ecg
- ComplexDFT

**Statistics** (3 folders → 1 module):
- Signal_Mean
- SignalVariance
- Signal_Standard_deviation

**Convolution** (1 folder → 1 module):
- Convolution

**Utilities** (1 folder + new → 1 module):
- Rect_to_Polar
- Signal generation functions (new)
- Normalization utilities (new)

**Signal I/O** (consolidated):
- Common waveform data
- Unified read/write functions
- Dynamic memory management

### Known Issues
None currently. This is a stable initial release.

### Future Enhancements
See LIBRARY_OVERVIEW.md for potential additions including:
- FFT (Fast Fourier Transform)
- Window functions (Hamming, Hanning, etc.)
- FIR/IIR filter design
- Spectrograms (STFT)
- Resampling functions

---

**Note**: The original `signal_processing/` folder has been preserved for reference.
