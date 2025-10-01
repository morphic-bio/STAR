# Debug Reference Card

## Quick Start
```bash
# Enable debug mode
export STAR_DEBUG_TAG=1

# Run with debug logging
./source/STAR [normal options]

# Disable debug mode
unset STAR_DEBUG_TAG
```

## Environment Variables
| Variable | Purpose | Example |
|----------|---------|---------|
| `STAR_DEBUG_TAG` | Enable debug logging & validation | `STAR_DEBUG_TAG=1` |

## Debug Output Format
```
[DEBUG_TAG HH:MM:SS TXXXX] message
```
- `HH:MM:SS`: Timestamp
- `TXXXX`: Thread ID (last 4 digits)
- `message`: Debug information

## Key Debug Points

### BAMTagBuffer Allocation
```
[DEBUG_TAG 04:52:42 T0688] Reserved 1000000 entries for BAMTagBuffer::entries initial allocation
[DEBUG_TAG 04:52:42 T0688] BAMTagBuffer initialized with 1M entries capacity
```

### Binary Tag Writer
```
[DEBUG_TAG 04:52:51 T0688] Binary tag header: statusBits=1, cbBits=20, umiBits=24, recordCount=21
[DEBUG_TAG 04:52:51 T0688] Record 0: status=1, cbIdx=452607, umiPacked=8333434, recordIndex=0
[DEBUG_TAG 04:52:51 T0688] File size validation: actual=158 bytes, expected~158 bytes
```

### BAM Tag Injection
```
[DEBUG_TAG 04:52:51 T0688] Allocated BAM buffers: bam0=100008 bytes, bam1=10000 bytes
[DEBUG_TAG 04:52:51 T0688] Record 0: iread=832, hasCB=true, hasUMI=true
[DEBUG_TAG 04:52:51 T0688] BAM processing completed: 21 records processed
```

## Error Types

### Critical Errors (Fatal)
- `ERROR: Allocation failed in [context]` → Reduce threads/increase memory
- `ERROR: Invalid read index X in trailer` → Corrupted data/indexing bug
- `ERROR: BAM record size X exceeds buffer` → Buffer overflow protection

### Warnings (Non-Fatal)
- `WARNING: Record index out of order` → Threading/sorting issue
- `WARNING: Expected to write X records but wrote Y` → Count mismatch
- `WARNING: Not all entries have valid readIds` → Data filtering

## Common Usage Patterns

### Development/Testing
```bash
STAR_DEBUG_TAG=1 ./source/STAR --runThreadN 2 [options] 2>&1 | tee debug.log
```

### Production Debug Run
```bash
./runSTAR_debug.sh  # Full production run with debug instrumentation
```

### Memory Testing
```bash
STAR_DEBUG_TAG=1 ./mem_test_tags.sh
```

### Production (No Debug)
```bash
./source/STAR --runThreadN 8 [options]  # No debug overhead
```

### CI/Automated Testing
```bash
if STAR_DEBUG_TAG=1 ./test.sh 2>&1 | grep -q "ERROR:"; then
  echo "Validation failed"
  exit 1
fi
```

## Performance Impact
- **Production**: 0% overhead (debug code compiled out)
- **Debug mode**: ~5-10% slower, minimal memory increase

## See Also
- [docs/debug_instrumentation.md](docs/debug_instrumentation.md) - Complete guide
- [docs/memory_testing_guide.md](docs/memory_testing_guide.md) - Memory testing
- [MEMORY_TESTING.md](MEMORY_TESTING.md) - Memory testing quick reference
