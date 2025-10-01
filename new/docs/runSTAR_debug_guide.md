# runSTAR_debug.sh Usage Guide

## Overview

`runSTAR_debug.sh` is a debug-instrumented version of `runSTAR.sh` designed to track down memory errors and data corruption issues in STAR's tag table functionality. It enables comprehensive debug logging, memory debugging, and real-time analysis to help identify the root cause of issues like the trailer parsing bug discovered during testing.

## Key Differences from runSTAR.sh

### Debug Features Enabled
- **`STAR_DEBUG_TAG=1`** - Enables comprehensive debug logging
- **`MALLOC_CHECK_=2`** - Enables glibc malloc debugging  
- **`MALLOC_PERTURB_=42`** - Perturbs freed memory to catch use-after-free bugs
- **Real-time log analysis** - Automatically analyzes debug output for issues
- **Resource monitoring** - Tracks memory and CPU usage during run
- **Comprehensive error reporting** - Categorizes and reports all detected issues

### Configuration Changes
- **Threads reduced** from 24 to 8 for debug stability
- **Buffer sizes reduced** for faster debug iterations
- **Compression reduced** from 6 to 1 for faster I/O
- **Debug output directory** uses `_debug` suffix
- **Both tag modes enabled** - `--soloWriteTagTable` and `--soloAddTagsToUnsorted`

## Usage

### Basic Usage
```bash
./runSTAR_debug.sh
```

### Prerequisites
1. **Debug-instrumented STAR binary** in `source/STAR`
2. **Input data** available in `/storage/JAX_sequences/SC2300772/`
3. **Sufficient disk space** for debug logs and output
4. **Root access** (optional) for enhanced malloc debugging

### Monitoring Progress
```bash
# Monitor debug output in real-time
tail -f /tmp/star_debug_logs/star_debug_*.log

# Monitor resource usage
tail -f /tmp/star_debug_logs/star_debug_*.log.resources

# Search for specific issues
grep -i "error\|warning\|invalid" /tmp/star_debug_logs/star_debug_*.log
```

## Output Structure

### Debug Logs
- **Main log**: `/tmp/star_debug_logs/star_debug_YYYYMMDD_HHMMSS.log`
- **Resource log**: `/tmp/star_debug_logs/star_debug_YYYYMMDD_HHMMSS.log.resources`

### Output Files
- **Debug output directory**: `/storage/Alignments/SC2300772_debug/`
- **Binary tag table**: `Aligned.out.cb_ub.bin`
- **BAM with tags**: `Aligned.out.bam`
- **Solo matrices**: `Solo.out/` directory

## Debug Output Analysis

### Automatic Analysis
The script automatically analyzes debug output and reports:

#### ✅ Success Indicators
- Debug message counts
- Successful file creation
- Resource usage statistics
- Clean exit code (0)

#### ❌ Error Detection
- **Memory corruption**: `Invalid read index` messages
- **Allocation failures**: `Allocation failed` messages  
- **Data consistency**: Record count mismatches
- **File format**: Binary format validation errors

#### ⚠️ Warning Detection
- Record index ordering issues
- Buffer size warnings
- Resource usage alerts

### Manual Analysis Commands
```bash
# Count debug message types
grep -c "DEBUG_TAG" /tmp/star_debug_logs/star_debug_*.log
grep -c "ERROR:" /tmp/star_debug_logs/star_debug_*.log
grep -c "WARNING:" /tmp/star_debug_logs/star_debug_*.log

# Find specific issues
grep -A 5 -B 5 "Invalid read index" /tmp/star_debug_logs/star_debug_*.log
grep -A 3 "Allocation failed" /tmp/star_debug_logs/star_debug_*.log
grep "Record.*out of order" /tmp/star_debug_logs/star_debug_*.log

# Analyze resource usage
grep "Memory usage:" /tmp/star_debug_logs/star_debug_*.log.resources
```

## Common Issues and Solutions

### Issue: Invalid Read Index Error
```
[DEBUG_TAG] ERROR: Invalid iread=3573412790272 at record 0 (readInfo.size=20001)
```
**Cause**: Memory corruption in trailer parsing
**Investigation**: 
- Check tmp file integrity
- Verify readInfo array bounds
- Look for upstream indexing bugs

### Issue: Allocation Failures
```
ERROR: Allocation failed in BAMTagBuffer::entries (current size=1000000)
```
**Cause**: Insufficient memory
**Solutions**:
- Reduce `--runThreadN` further
- Increase system RAM
- Use smaller input dataset

### Issue: Record Index Out of Order
```
[DEBUG_TAG] WARNING: Record index out of order at entry 1500: current=1499, previous=1500
```
**Cause**: Threading race condition
**Investigation**:
- Check sorting logic
- Verify thread synchronization
- Review BAMTagBuffer locking

### Issue: File Size Mismatch
```
[DEBUG_TAG] File size validation: actual=1580 bytes, expected~1520 bytes
```
**Cause**: Binary format corruption
**Investigation**:
- Check bit packing logic
- Verify header calculations
- Review BAMTagBinaryWriter implementation

## Performance Expectations

### Debug Mode Overhead
- **Runtime**: 2-3x slower than production
- **Memory**: ~20% increase due to debug structures
- **Disk I/O**: Higher due to comprehensive logging

### Typical Debug Run Times
- **Small dataset** (mem_test): ~30 seconds
- **Full SC2300772**: ~2-4 hours (vs 1-2 hours production)

## Integration with Other Tools

### With AddressSanitizer
```bash
# Build with ASan first
cd source
ASAN=1 make clean
ASAN=1 make STAR

# Run debug script with ASan binary
./runSTAR_debug.sh
```

### With Valgrind
```bash
# Modify runSTAR_debug.sh to use valgrind
valgrind --tool=memcheck --leak-check=full --track-origins=yes \
  ./source/STAR [parameters] 2>&1 | tee valgrind_debug.log
```

### With GDB
```bash
# For interactive debugging
gdb --args ./source/STAR [parameters]
(gdb) set environment STAR_DEBUG_TAG=1
(gdb) run
```

## Customization

### Modify Thread Count
```bash
# Edit runSTAR_debug.sh
THREADS=4  # Further reduce for more stability
```

### Enable Additional Debug Features
```bash
# Add to runSTAR_debug.sh
export STAR_DEBUG_SOLO=1        # Additional Solo debugging
export STAR_DEBUG_MEMORY=1      # Memory allocation tracking
export MALLOC_TRACE=debug.mtrace # Malloc tracing
```

### Custom Log Analysis
```bash
# Add custom analysis function
analyze_custom_debug() {
    local log_file=$1
    echo "=== Custom Analysis ==="
    
    # Look for specific patterns
    grep "BAMTagBuffer" "$log_file" | wc -l
    grep "writeTagBinary" "$log_file" | head -5
    
    # Custom metrics
    local first_append=$(grep -m1 "First BAMTagBuffer::append" "$log_file")
    echo "First append: $first_append"
}
```

## Troubleshooting

### Script Won't Start
1. **Check permissions**: `chmod +x runSTAR_debug.sh`
2. **Verify paths**: Ensure input directories exist
3. **Check disk space**: Ensure sufficient space for logs
4. **Test STAR binary**: Run `./source/STAR --version`

### No Debug Output
1. **Verify instrumentation**: Run `./test_debug.sh`
2. **Check environment**: Ensure `STAR_DEBUG_TAG=1` is set
3. **Rebuild STAR**: `cd source && make clean && make STAR`

### High Memory Usage
1. **Reduce threads**: Set `THREADS=2` or `THREADS=1`
2. **Reduce buffers**: Lower `limitIObufferSize` values
3. **Use smaller dataset**: Test with subset of input files

### Log Files Too Large
1. **Filter debug output**: Modify logging levels in source code
2. **Rotate logs**: Implement log rotation in script
3. **Focus analysis**: Use `grep` to extract specific issues

## See Also

- [debug_instrumentation.md](debug_instrumentation.md) - Complete debug system guide
- [DEBUG_REFERENCE.md](../DEBUG_REFERENCE.md) - Quick debug reference
- [memory_testing_guide.md](memory_testing_guide.md) - Memory testing with ASan
- [runSTAR.sh](../runSTAR.sh) - Original production script
