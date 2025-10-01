#ifndef CODE_BAMTagBinaryWriter
#define CODE_BAMTagBinaryWriter

#include <iostream>
#include <fstream>
#include <cstdint>

// Binary tag stream format:
// - Header: 4x uint64_t (statusBits, cbBits, umiBits, recordCount)  
// - Records: bit-packed (status + CB index + packed UMI), byte-aligned per record
struct TagStreamHeader {
    uint64_t statusBits;  // Always 1 for current format
    uint64_t cbBits;      // Bits needed for CB indices (0 = missing)
    uint64_t umiBits;     // Bits needed for packed UMI values
    uint64_t recordCount; // Number of payload records that follow
    
    TagStreamHeader() : statusBits(1), cbBits(0), umiBits(0), recordCount(0) {}
    TagStreamHeader(uint64_t cb, uint64_t umi, uint64_t count) : statusBits(1), cbBits(cb), umiBits(umi), recordCount(count) {}
};

class BAMTagBinaryWriter {
public:
    BAMTagBinaryWriter();
    ~BAMTagBinaryWriter();
    
    // Write header to stream (4x uint64_t in little-endian)
    void writeHeader(std::ostream& out, const TagStreamHeader& header);
    
    // Write a single record to stream (bit-packed, byte-aligned)
    void writeRecord(std::ostream& out, uint8_t status, uint64_t cbIdx, uint64_t umiPacked, const TagStreamHeader& header);
    
    // Utility: compute ceil(log2(n)) for bit width calculation
    static uint64_t integerLog2ceil(uint64_t n);

private:
    // Helper to write little-endian uint64_t
    void writeLittleEndian64(std::ostream& out, uint64_t value);
    
    // Helper to pack bits into uint64_t buffer and write when byte-aligned
    void packAndWrite(std::ostream& out, uint64_t value, uint64_t bits, uint64_t& buffer, uint64_t& bufferBits);
    
    // Helper to flush remaining bits in buffer (pad to byte boundary)
    void flushBuffer(std::ostream& out, uint64_t buffer, uint64_t bufferBits);
};

#endif
