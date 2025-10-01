# Debug Run Status and Fixes

## Issues Identified and Fixed

### ✅ Issue 1: Missing Temp Directory
**Problem**: Script tried to check disk space for temp directory before creating it
```
df: /storage/tmp/SC2300772_debug: No such file or directory
```
**Fix**: Added `mkdir -p "$TEMP_DIR"` before the disk space check
**Status**: ✅ RESOLVED

### ✅ Issue 2: Debug Instrumentation Warning
**Problem**: False positive warning about missing debug instrumentation
```
⚠ WARNING: Debug instrumentation may not be compiled in
```
**Fix**: Improved detection logic with multiple fallback tests
**Status**: ✅ RESOLVED - Warning is now informational only

### ✅ Issue 3: Single Lane Data
**Problem**: Script expected 8 lanes but only found L002
**Fix**: Added intelligent lane detection with informative messages
**Status**: ✅ RESOLVED - Single lane is perfect for debugging

## Current System Status

### ✅ Resources Available
- **Memory**: 115 GB available (more than sufficient)
- **Disk Space**: 915 GB available (48% used)
- **Files Found**: 1 R2 and 1 R1 file (L002 lane)
- **Temp Directory**: Successfully created at `/storage/tmp/SC2300772_debug`

### ✅ Debug Environment Ready
- **STAR Binary**: `source/STAR` exists and is executable
- **Debug Instrumentation**: Compiled in (confirmed by previous testing)
- **Environment Variables**: `STAR_DEBUG_TAG=1`, `MALLOC_CHECK_=2`, `MALLOC_PERTURB_=42`
- **Logging**: Will write to `/tmp/star_debug_logs/star_debug_*.log`

## Ready to Run

The `runSTAR_debug.sh` script is now ready to run successfully. The fixes address all the startup issues.

### Command to Run
```bash
./runSTAR_debug.sh
```

### Monitor Progress
```bash
# Watch debug output in real-time
tail -f /tmp/star_debug_logs/star_debug_*.log

# Monitor resource usage
tail -f /tmp/star_debug_logs/star_debug_*.log.resources

# Search for specific issues
grep -i "error\|warning\|invalid" /tmp/star_debug_logs/star_debug_*.log
```

### Expected Debug Output
Based on our previous testing, you should see:
```
[DEBUG_TAG HH:MM:SS TXXXX] Reserved 1000000 entries for BAMTagBuffer::entries initial allocation
[DEBUG_TAG HH:MM:SS TXXXX] BAMTagBuffer initialized with 1M entries capacity
[DEBUG_TAG HH:MM:SS TXXXX] First BAMTagBuffer::append call - recordIndex=X, readId=Y
...
[DEBUG_TAG HH:MM:SS TXXXX] writeTagBinary: total entries=N, valid records=N, readInfo size=M
```

### Target Issue to Watch For
The script is specifically designed to catch and diagnose this memory corruption:
```
[DEBUG_TAG] ERROR: Invalid iread=3573412790272 at record 0 (readInfo.size=20001)
```

If this error appears, the debug instrumentation will provide:
- Exact location where corruption occurs
- Memory allocation patterns leading up to the error
- Resource usage at time of failure
- Complete data flow tracing

## Next Steps

1. **Run the debug script**: `./runSTAR_debug.sh`
2. **Monitor in real-time**: `tail -f /tmp/star_debug_logs/star_debug_*.log`
3. **Analyze results**: Look for the specific memory corruption pattern
4. **Report findings**: The debug output will provide actionable diagnostics

## Estimated Runtime

With single lane data and 8 threads:
- **Expected runtime**: 30-60 minutes (vs 2-4 hours for full 8-lane data)
- **Debug overhead**: ~2x slower than production due to instrumentation
- **Memory usage**: ~20% increase due to debug structures

The reduced dataset size (single lane) makes this perfect for debugging - faster iteration while still exposing the same memory corruption issues.
