#ifndef CODE_BAMTagBuffer
#define CODE_BAMTagBuffer

#include "IncludeDefine.h"
#include "BAMoutput.h"
#include "SoloCommon.h"
#include <string>
#include <vector>
#include <mutex>
#include <cstdint>

struct BAMTagEntry {
    uint32_t recordIndex; // BAM record index for ordering
    uint32_t readId;      // iReadAll truncated to 32 bits
    
    BAMTagEntry() : recordIndex(0), readId(0) {}
    BAMTagEntry(const BAMRecordMeta& meta) 
        : recordIndex(static_cast<uint32_t>(meta.recordIndex)), readId(static_cast<uint32_t>(meta.iReadAll)) {}
};

class BAMTagBuffer {
public:
    BAMTagBuffer();
    ~BAMTagBuffer();
    
    // Thread-safe append of alignment metadata
    void append(const BAMRecordMeta& meta);
    
    // Write binary tag stream to file (called at end of processing)
    void writeTagBinary(const std::string& outputPath,
                        const std::vector<readInfoStruct>& readInfo,
                        uint64_t cbBits,
                        uint64_t umiBits);
    
    // Clear buffer to free memory
    void clear();
    
    // Get number of entries
    size_t size() const { return entries.size(); }

private:
    std::vector<BAMTagEntry> entries;
    std::mutex entriesMutex; // Thread safety for append operations
};

#endif
