# Memory Testing Quick Reference

## Build with AddressSanitizer
```bash
cd source/
ASAN=1 make clean
ASAN=1 make STAR
```

## Run Memory Test
```bash
# Basic test
./mem_test_tags.sh

# With leak detection
ASAN_OPTIONS="detect_leaks=1" ./mem_test_tags.sh 2>&1 | tee mem_test.log
```

## Configuration (Environment Variables)
```bash
export STAR_BIN="./source/STAR"
export GENOME_DIR="/path/to/genome/index"
export FASTQ_DIR="/path/to/fastq/files"
export WHITELIST="/path/to/cb_whitelist.txt"
export THREADS="4"
```

## Expected Outputs
- `mem_test_output/Aligned.out.cb_ub.bin` - Binary tag table
- `mem_test_output/Aligned.out.bam` - BAM with CB/UB tags
- ASan reports in terminal/log file

## Rebuild Production Binary
```bash
cd source/
make clean
make STAR
```

## Critical vs Non-Critical Errors
**Critical (fix required):**
- `heap-buffer-overflow`
- `heap-use-after-free`
- `attempting double-free`

**Non-Critical (typical):**
- `detected memory leaks` (end-of-program)

For detailed instructions, see `docs/memory_testing_guide.md`
