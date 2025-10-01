#ifndef CODE_BAMoutput
#define CODE_BAMoutput

#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H
#include "Parameters.h"

// Forward declaration to avoid circular dependency
class BAMTagBuffer;

struct BAMRecordMeta {
    uint64_t recordIndex;
    uint64_t iReadAll;
    uint32_t mate;      // 0 or 1
    uint32_t alignIdx;  // multi-mapper ordinal
    const char* qname;  // pointer or nullptr
};

class BAMoutput {//
public:
    //sorted output
    BAMoutput (int iChunk, string tmpDir, Parameters &Pin);
    void coordOneAlign (char *bamIn, uint bamSize, uint iRead);
    void coordBins ();
    void coordFlush ();
    //unsorted output
    BAMoutput (BGZF *bgzfBAMin, Parameters &Pin);
    BAMoutput (BGZF *bgzfBAMin, Parameters &Pin, BAMTagBuffer* tagBufferIn);
    void unsortedOneAlign (char *bamIn, uint bamSize, uint bamSize2);
    void unsortedOneAlign (char *bamIn, uint bamSize, uint bamSize2, const BAMRecordMeta& meta);
    void unsortedFlush ();
    void coordUnmappedPrepareBySJout();

    uint32 nBins; //number of bins to split genome into
    uint* binTotalN; //total number of aligns in each bin
    uint* binTotalBytes;//total size of aligns in each bin
private:
    uint64 bamArraySize; //this size will be allocated
    char* bamArray; //large array to store the bam alignments, pre-sorted
    uint64 binSize, binSize1;//storage size of each bin
    uint64 binGlen;//bin genomic length
    char **binStart; //pointers to starts of the bins
    uint64 *binBytes, binBytes1;//number of bytes currently written to each bin
    ofstream **binStream;//output streams for each bin
    BGZF *bgzfBAM;
    Parameters &P;
    string bamDir;
    BAMTagBuffer* tagBuffer; // Optional tag buffer for metadata collection
};

#endif
