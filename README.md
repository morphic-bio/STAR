STAR 2.7.11b (Fork)
===================
Spliced Transcripts Alignment to a Reference  
© Alexander Dobin, 2009-2024  
https://www.ncbi.nlm.nih.gov/pubmed/23104886

**This is a fork** of the upstream STAR repository with extended STARsolo functionality for unsorted BAM workflows, tag table export, and enhanced gene annotation.

UPSTREAM AUTHOR/SUPPORT
========================
Alex Dobin, dobin@cshl.edu  
https://github.com/alexdobin/STAR/issues  
https://groups.google.com/d/forum/rna-star

FORK ENHANCEMENTS
==================

This fork adds three major enhancements to STAR 2.7.11b:

### 1. Two-Pass Unsorted BAM with CB/UB Tags (`--soloAddTagsToUnsorted`)
Enables corrected cell barcode (CB) and UMI (UB) tags in unsorted BAM output. When combined with `--outSAMtype BAM Unsorted`, STAR now captures alignments during pass 1, then replays them in pass 2 after Solo correction to inject accurate barcode tags.

**Quick Usage:**
```bash
STAR --outSAMtype BAM Unsorted \
     --soloAddTagsToUnsorted yes \
     --soloType CB_UMI_Simple \
     --soloFeatures Gene GeneFull \
     [other parameters...]
```

### 2. Binary Tag Table Export (`--soloWriteTagTable`)
Exports corrected CB/UB assignments to a compact binary sidecar file without rewriting the entire BAM. Useful for audit trails, secondary analysis, or lightweight barcode extraction.

**Quick Usage:**
```bash
STAR --soloWriteTagTable Default \
     --soloType CB_UMI_Simple \
     --soloFeatures Gene \
     [other parameters...]
```

### 3. Custom Gene Annotation Tags (ZG/ZX)
Two new BAM tags that provide comprehensive gene annotation beyond standard GX/GN tags:
- **ZG**: Comma-separated list of Ensembl gene IDs for all overlapping genes
- **ZX**: Genomic overlap classification (`exonic`, `intronic`, `none`, `spanning`)

**Quick Usage:**
```bash
STAR --outSAMattributes NH HI AS nM NM CR CY UR UY GX GN gx gn ZG ZX \
     --soloFeatures Gene GeneFull \
     --soloStrand Unstranded \
     [other parameters...]
```

**Key Benefits:**
- 28% more gene annotations than standard GX tags
- 99.8% read coverage vs 79.2% with GX tags
- Perfect 100% concordance with existing GX tags
- Production validated on large-scale datasets

DOCUMENTATION
==============

### Fork-Specific Documentation
Located in `new/docs/`:

- **[TECHNICAL_NOTES.md](new/docs/TECHNICAL_NOTES.md)** - Complete technical implementation details for all fork features
- **[ZG_ZX_Implementation_Summary.md](new/docs/ZG_ZX_Implementation_Summary.md)** - ZG/ZX tag technical documentation with code references
- **[two_pass_unsorted_usage.md](new/docs/two_pass_unsorted_usage.md)** - Usage guide for unsorted BAM workflows
- **[CHANGES_FORK.md](new/docs/CHANGES_FORK.md)** - Detailed changelog for fork modifications
- **[RELEASEnotes.md](new/docs/RELEASEnotes.md)** - Fork release notes and validation results
- **[memory_testing_guide.md](new/docs/memory_testing_guide.md)** - AddressSanitizer testing procedures
- **[debug_instrumentation.md](new/docs/debug_instrumentation.md)** - Debug logging and validation

### Upstream STAR Documentation
- **[STARmanual.pdf](doc/STARmanual.pdf)** - Original STAR manual
- **[CHANGES.md](CHANGES.md)** - Upstream STAR changes and release history

### Production Scripts
Located in `new/scripts/`:
- `runSTAR.sh` - Production script with comprehensive ZG/ZX configuration
- `runSTAR_debug.sh` - Debug-enabled run with monitoring
- `validate_zg_zx.py` - Validation script for ZG/ZX tag correctness

### Test Suites
Located in `new/tests/`:
- `emit_test.sh` - Binary tag stream and unsorted BAM validation
- `integration_test.sh` - End-to-end CB/UB patching tests
- `mem_test_tags.sh` - AddressSanitizer memory safety tests

HARDWARE/SOFTWARE REQUIREMENTS
================================
- x86-64 compatible processors
- 64 bit Linux or Mac OS X
- At least 16GB RAM for mammalian genomes (32GB recommended)

DIRECTORY CONTENTS
===================
```
/
├── source/          # All source files required for compilation
├── bin/             # Pre-compiled executables for Linux and Mac OS X
├── doc/             # Upstream STAR documentation
├── extras/          # Miscellaneous files and scripts
├── tools/           # Binary tag decoder and utilities
└── new/             # Fork-only material
    ├── docs/        # Fork documentation (see TECHNICAL_NOTES.md)
    ├── scripts/     # Production and validation scripts
    ├── tests/       # Test harnesses
    ├── plans/       # Development planning documents
    └── testing/     # Test data and outputs
```

COMPILING FROM SOURCE
======================

### Standard Compilation (Linux)
```bash
cd source
make STAR
```

For processors without AVX support:
```bash
make STAR CXXFLAGS_SIMD=sse
```

### Mac OS X Compilation
```bash
# Install brew and gcc
brew install gcc

# Build STAR (adjust gcc version path as needed)
cd source
make STARforMacStatic CXX=/usr/local/Cellar/gcc/8.2.0/bin/g++-8

# Install to system path
cp STAR /usr/local/bin
```

### Non-Standard GCC
If g++ is not on the path:
```bash
cd source
make STAR CXX=/path/to/g++
```

### Platform-Specific Optimization
```bash
# Platform-specific optimization
make CXXFLAGSextra=-march=native

# With link-time optimization
make LDFLAGSextra=-flto CXXFLAGSextra="-flto -march=native"
```

### AddressSanitizer Build (Memory Testing)
```bash
cd source
ASAN=1 make clean
ASAN=1 make STAR

# Run memory tests
cd ..
ASAN_OPTIONS="detect_leaks=1" ./new/tests/mem_test_tags.sh
```

**Note:** ASan builds require ~3x more RAM and run 2-5x slower than production builds.

USAGE EXAMPLES
===============

### Example 1: Unsorted BAM with CB/UB Tags
```bash
STAR --runThreadN 24 \
     --genomeDir /path/to/genome \
     --readFilesIn R2.fastq.gz R1.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM Unsorted \
     --outSAMattributes NH HI AS nM NM CR CY UR UY \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist /path/to/whitelist.txt \
     --soloFeatures Gene GeneFull \
     --soloAddTagsToUnsorted yes \
     --outFileNamePrefix output/
```

### Example 2: Tag Table Export Only
```bash
STAR --runThreadN 24 \
     --genomeDir /path/to/genome \
     --readFilesIn R2.fastq.gz R1.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM Unsorted \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist /path/to/whitelist.txt \
     --soloFeatures Gene \
     --soloWriteTagTable Default \
     --outFileNamePrefix output/
```

### Example 3: Full Feature Set (Unsorted BAM + Tag Table + ZG/ZX)
```bash
STAR --runThreadN 24 \
     --genomeDir /path/to/genome \
     --readFilesIn R2.fastq.gz R1.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM Unsorted \
     --outSAMattributes NH HI AS nM NM CR CY UR UY GX GN gx gn ZG ZX \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist /path/to/whitelist.txt \
     --soloFeatures Gene GeneFull \
     --soloStrand Unstranded \
     --soloAddTagsToUnsorted yes \
     --soloWriteTagTable Default \
     --quantMode GeneCounts \
     --outFileNamePrefix output/
```

### Example 4: Debug Mode with Full Instrumentation
```bash
# Use production debug script
./new/scripts/runSTAR_debug.sh

# Or enable debug logging manually
STAR_DEBUG_TAG=1 STAR [parameters...]
```

INSPECTING OUTPUT
==================

### Examine ZG/ZX Tags
```bash
# View first 10 reads with ZG/ZX tags
samtools view output/Aligned.out.bam | grep -E "ZG:Z:|ZX:Z:" | head -10

# Extract ZG/ZX tags for analysis
samtools view output/Aligned.out.bam | awk '{
    for(i=1;i<=NF;i++) 
        if($i~/^ZG:Z:/ || $i~/^ZX:Z:/) 
            print $1"\t"$i
}' > zg_zx_tags.txt

# Validate ZG/ZX tags
python3 new/scripts/validate_zg_zx.py output/Aligned.out.bam allowed_genes.txt
```

### Decode Binary Tag Table
```bash
# Compile decoder
cd tools
make

# Decode tag table
./decode_tag_binary ../output/Aligned.out.cb_ub.bin > decoded_tags.txt
```

LIMITATIONS
============
This release has been tested with default parameters for human and mouse genomes.
Mammalian genomes require at least 16GB of RAM, ideally 32GB.

**Fork-Specific Notes:**
- Two-pass unsorted BAM requires sufficient disk space for temporary files
- Tag table binary format is specific to this fork; use provided decoder tool
- ZG/ZX tags are BAM-only (not emitted in SAM format)
- GeneFull must be in `--soloFeatures` for ZG/ZX tags to populate

FREEBSD PORTS
==============
STAR can be installed on FreeBSD via the FreeBSD ports system:
```bash
pkg install star
```
**Note:** FreeBSD ports may not include fork-specific features. Compile from source for full functionality.

DOCKER
=======
Docker build scripts are available in `new/scripts/`:
```bash
./new/scripts/docker_build.sh
./new/scripts/runSTAR_docker.sh
```

CONTRIBUTING
=============
For issues or contributions related to:
- **Upstream STAR**: Use https://github.com/alexdobin/STAR/issues
- **Fork-specific features**: Contact the fork maintainer or create issues in the fork repository

LICENSE
========
See [LICENSE](LICENSE) file for details.

FUNDING
========
The development of STAR is supported by the National Human Genome Research Institute of
the National Institutes of Health under Award Number R01HG009318.
The content is solely the responsibility of the authors and does not necessarily represent 
the official views of the National Institutes of Health.
