#ifndef CODE_BAMoutputSoloTmp
#define CODE_BAMoutputSoloTmp

#include "IncludeDefine.h"
#include "Parameters.h"
#include "BAMoutput.h" // for BAMRecordMeta

// Forward declaration to avoid circular dependency
class BAMTagBuffer;

class BAMoutputSoloTmp {
public:
    // Constructor for solo tmp mode
    BAMoutputSoloTmp(ofstream* tmpStreamIn, Parameters &Pin);
    BAMoutputSoloTmp(ofstream* tmpStreamIn, Parameters &Pin, BAMTagBuffer* tagBufferIn);
    ~BAMoutputSoloTmp();
    
    // Main method to write one alignment with iRead trailer
    void unsortedOneAlign(char *bamIn, uint bamSize, uint iReadAll);
    void unsortedOneAlign(char *bamIn, uint bamSize, uint iReadAll, const BAMRecordMeta& meta);
    
    // Flush buffer to file
    void flush();
    
    // Close and finalize
    void close();

private:
    ofstream* tmpStream;
    Parameters &P;
    
    uint64 bufferSize;      // size of the buffer
    char* buffer;           // buffer to store alignments before writing
    uint64 bufferUsed;      // bytes currently used in buffer
    BAMTagBuffer* tagBuffer; // Optional tag buffer for metadata collection
    
    static const uint64 defaultBufferSize = 1024*1024*64; // 64MB buffer
    // BAMoutput_oneAlignMaxBytes is defined in IncludeDefine.h
};

#endif
