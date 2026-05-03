# Makefile for Signal Processing Library
# Works with GCC on Linux/macOS/MinGW

CC = gcc
CFLAGS = -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm

# Source files
SRC_DIR = src
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(SRCS:.c=.o)

# Example files
EXAMPLE_DIR = examples
EXAMPLES = $(wildcard $(EXAMPLE_DIR)/*.c)
EXAMPLE_BINS = $(EXAMPLES:.c=)

# Output directory for binaries
BIN_DIR = bin

.PHONY: all examples clean help

all: examples

# Create bin directory
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Build all examples
examples: $(BIN_DIR) $(EXAMPLE_BINS)

# Rule to build each example
$(EXAMPLE_DIR)/%: $(EXAMPLE_DIR)/%.c $(SRCS)
	$(CC) $(CFLAGS) -o $(BIN_DIR)/$(notdir $@) $< $(SRCS) $(LDFLAGS)
	@echo "Built: $(notdir $@)"

# Clean build artifacts
clean:
	rm -rf $(BIN_DIR)
	rm -f $(OBJS)
	rm -f *.dat
	rm -f $(EXAMPLE_DIR)/*.dat
	@echo "Cleaned build artifacts"

# Help message
help:
	@echo "Signal Processing Library - Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  make all       - Build all examples (default)"
	@echo "  make examples  - Build all example programs"
	@echo "  make clean     - Remove build artifacts"
	@echo "  make help      - Show this help message"
	@echo ""
	@echo "Example binaries will be in the bin/ directory"
	@echo ""
	@echo "To run examples:"
	@echo "  ./bin/example_dft"
	@echo "  ./bin/example_convolution"
	@echo "  ./bin/example_stats_utils"
