#ifndef CODE_BAMTagBuffer
#define CODE_BAMTagBuffer

#include "IncludeDefine.h"
#include "BAMoutput.h"
#include "SoloCommon.h"
#include <string>
#include <vector>
#include <mutex>
#include <cstdint>

#ifdef SOLO_USE_PACKED_READINFO
#include "PackedReadInfo.h"
#endif

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

    void append(const BAMRecordMeta& meta);

#ifndef SOLO_USE_PACKED_READINFO
    void writeTagBinary(const std::string& outputPath,
                        const std::vector<readInfoStruct>& readInfo,
                        uint64_t cbBits,
                        uint64_t umiBits);
#else
    void writeTagBinaryPacked(const std::string& outputPath,
                              const PackedReadInfo& packed,
                              uint64_t cbBits,
                              uint64_t umiBits);
#endif
    void clear();
    size_t size() const { return entries.size(); }
private:
    std::vector<BAMTagEntry> entries;
    std::mutex entriesMutex;
};
#endif
