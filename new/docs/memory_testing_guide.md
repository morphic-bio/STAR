# Memory Testing Guide for STAR Tag Table Functionality

This guide explains how to compile STAR with AddressSanitizer (ASan) and run memory testing for the tag table serialization functionality.

## Overview

The memory testing system uses AddressSanitizer to detect memory errors, leaks, and other memory-related issues in STAR's tag table functionality. This is particularly important for validating the binary serialization implementation.

## Prerequisites

- GCC or Clang compiler with AddressSanitizer support
- Standard STAR build dependencies (htslib, etc.)
- Sufficient RAM (ASan increases memory usage ~3x)
- Test data (genome index, FASTQ files, whitelist)

## Quick Start

```bash
# 1. Build STAR with AddressSanitizer
cd /path/to/STAR/source
ASAN=1 make clean
ASAN=1 make STAR

# 2. Run memory test
cd ..
./mem_test_tags.sh

# 3. Rebuild production binary (optional)
cd source
make clean
make STAR
```

## Detailed Instructions

### Step 1: Enable AddressSanitizer Build

The STAR Makefile includes conditional ASan support. When `ASAN=1` is set, it adds the following flags:

```makefile
ifdef ASAN
CFLAGS   += -fsanitize=address -fno-omit-frame-pointer -O1
CXXFLAGS += -fsanitize=address -fno-omit-frame-pointer -O1
LDFLAGS  += -fsanitize=address
endif
```

### Step 2: Build Instrumented Binary

```bash
cd source/
ASAN=1 make clean    # Clean previous builds
ASAN=1 make STAR     # Build with AddressSanitizer
```

**Verify ASan linkage:**
```bash
ldd STAR | grep asan
# Should output: libasan.so.X => /path/to/libasan.so.X
```

### Step 3: Configure Test Environment

The `mem_test_tags.sh` script accepts environment variables for configuration:

```bash
# Required paths (defaults shown)
export STAR_BIN="./source/STAR"
export GENOME_DIR="/storage/scRNAseq_output/indices-98-32/star"
export FASTQ_DIR="/storage/downsampled/SC2300771"
export WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"

# Optional settings
export OUTDIR="mem_test_output"        # Output directory
export THREADS="4"                     # Number of threads
```

### Step 4: Run Memory Test

```bash
# Basic run
./mem_test_tags.sh

# With custom ASan options
ASAN_OPTIONS="detect_leaks=1:abort_on_error=1" ./mem_test_tags.sh

# Capture full log
ASAN_OPTIONS="detect_leaks=1" ./mem_test_tags.sh 2>&1 | tee mem_test_full.log
```

### Step 5: Analyze Results

The test produces several outputs:

**Success Indicators:**
- `✓ Binary tag table created successfully`
- `✓ BAM file created successfully`
- STAR completes with "finished successfully"

**Memory Analysis:**
- ASan reports appear after STAR completion
- Check for `ERROR: AddressSanitizer` messages
- Review leak summary at end of output

## ASan Options Reference

Common `ASAN_OPTIONS` environment variable settings:

```bash
# Basic leak detection
ASAN_OPTIONS="detect_leaks=1"

# Abort on first error (faster debugging)
ASAN_OPTIONS="detect_leaks=1:abort_on_error=1"

# Detailed stack traces
ASAN_OPTIONS="detect_leaks=1:symbolize=1:print_stacktrace=1"

# Suppress known third-party issues
ASAN_OPTIONS="detect_leaks=1:suppressions=asan_suppressions.txt"
```

## Test Script Details

The `mem_test_tags.sh` script performs these operations:

1. **Validation**: Checks for required files and directories
2. **FASTQ Discovery**: Automatically finds R1/R2 files (first 2 lanes for speed)
3. **STAR Execution**: Runs with minimal parameters to trigger tag table code:
   - `--soloWriteTagTable Default` (creates binary table)
   - `--soloAddTagsToUnsorted yes` (adds tags to BAM)
4. **Output Verification**: Confirms expected files are created

**Key STAR Parameters Used:**
```bash
--soloType CB_UMI_Simple
--soloCBlen 16 --soloUMIlen 12
--soloWriteTagTable Default        # Creates .bin file
--soloAddTagsToUnsorted yes        # Adds CB/UB tags to BAM
--outSAMtype BAM Unsorted
```

## Expected Output Files

After successful execution:

```
mem_test_output/
├── Aligned.out.bam              # BAM with CB/UB tags
├── Aligned.out.cb_ub.bin        # Binary tag table (compact format)
├── Log.final.out                # STAR completion log
├── Log.out                      # Detailed STAR log
└── Solo.out/                    # Solo matrices directory
```

## Troubleshooting

### Common Issues

**1. ASan Not Linked**
```bash
# Check linkage
ldd source/STAR | grep asan
# If no output, rebuild with ASAN=1
```

**2. Memory Exhaustion**
```bash
# Reduce thread count
export THREADS=2
# Or use smaller test dataset
```

**3. Missing Dependencies**
```bash
# Ensure paths exist
ls -la $GENOME_DIR $FASTQ_DIR $WHITELIST
```

**4. Compilation Errors**
```bash
# Check compiler version
gcc --version
# ASan requires GCC 4.8+ or Clang 3.1+
```

### Memory Error Types

**Critical Errors (require fixing):**
- Buffer overflows: `heap-buffer-overflow`
- Use after free: `heap-use-after-free`
- Double free: `attempting double-free`

**Non-Critical (typical for bioinformatics tools):**
- End-of-program leaks: `detected memory leaks`
- Uninitialized reads in third-party code

## Performance Impact

ASan introduces significant overhead:
- **Memory**: ~3x increase in RAM usage
- **Runtime**: ~2-5x slower execution
- **Disk**: Debug symbols increase binary size

**Recommendations:**
- Use smaller test datasets for routine testing
- Run on machines with adequate RAM (>32GB recommended)
- Consider using faster storage for temporary files

## Integration with CI/CD

For automated testing:

```bash
#!/bin/bash
# ci_memory_test.sh
set -euo pipefail

# Build with ASan
cd source
ASAN=1 make clean
ASAN=1 make STAR

# Run memory test with timeout
cd ..
timeout 1800 ./mem_test_tags.sh

# Check for critical errors (ignore end-of-program leaks)
if grep -q "ERROR: AddressSanitizer.*heap-\|use-after-free\|double-free" mem_test.log; then
    echo "CRITICAL memory errors detected!"
    exit 1
fi

echo "Memory test passed - no critical errors"
```

## Cleanup

After testing, rebuild production binary:

```bash
cd source/
make clean      # Remove ASan artifacts
make STAR       # Build optimized production binary

# Verify clean build
ldd STAR | grep asan || echo "Clean production build confirmed"
```

## Related Files

- `source/Makefile` - ASan build configuration
- `mem_test_tags.sh` - Main test script
- `asan_findings.md` - Example analysis report
- `tools/decode_tag_binary` - Binary format decoder for validation

## Additional Resources

- [AddressSanitizer Documentation](https://github.com/google/sanitizers/wiki/AddressSanitizer)
- [GCC Sanitizer Options](https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html)
- [Clang AddressSanitizer](https://clang.llvm.org/docs/AddressSanitizer.html)
