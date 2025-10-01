#include "BAMTagBinaryWriter.h"
#include <iostream>
#include <cstring>

BAMTagBinaryWriter::BAMTagBinaryWriter() {
}

BAMTagBinaryWriter::~BAMTagBinaryWriter() {
}

void BAMTagBinaryWriter::writeHeader(std::ostream& out, const TagStreamHeader& header) {
    // Write 4x uint64_t in little-endian format
    writeLittleEndian64(out, header.statusBits);
    writeLittleEndian64(out, header.cbBits);
    writeLittleEndian64(out, header.umiBits);
    writeLittleEndian64(out, header.recordCount);
}

void BAMTagBinaryWriter::writeRecord(std::ostream& out, uint8_t status, uint64_t cbIdx, uint64_t umiPacked, const TagStreamHeader& header) {
    // Pack bits: status (1 bit) + cbIdx (cbBits) + umiPacked (umiBits)
    uint64_t totalBits = header.statusBits + header.cbBits + header.umiBits;
    uint64_t buffer = 0;
    uint64_t bufferBits = 0;
    
    // Pack status bit (1 = OK, 0 = missing/invalid)
    packAndWrite(out, status & 1, header.statusBits, buffer, bufferBits);
    
    // Pack CB index (0 = missing when status=0)
    packAndWrite(out, cbIdx, header.cbBits, buffer, bufferBits);
    
    // Pack UMI (0 when status=0)
    packAndWrite(out, umiPacked, header.umiBits, buffer, bufferBits);
    
    // Flush any remaining bits to ensure byte alignment per record
    if (bufferBits > 0) {
        flushBuffer(out, buffer, bufferBits);
    }
}

uint64_t BAMTagBinaryWriter::integerLog2ceil(uint64_t n) {
    if (n <= 1) return 1;
    
    uint64_t bits = 0;
    uint64_t temp = n - 1;  // For exact powers of 2
    while (temp > 0) {
        temp >>= 1;
        bits++;
    }
    return bits;
}

void BAMTagBinaryWriter::writeLittleEndian64(std::ostream& out, uint64_t value) {
    // Write 8 bytes in little-endian order
    for (int i = 0; i < 8; i++) {
        uint8_t byte = static_cast<uint8_t>(value & 0xFF);
        out.write(reinterpret_cast<const char*>(&byte), 1);
        value >>= 8;
    }
}

void BAMTagBinaryWriter::packAndWrite(std::ostream& out, uint64_t value, uint64_t bits, uint64_t& buffer, uint64_t& bufferBits) {
    if (bits == 0) return;
    
    // Mask value to specified bit width
    uint64_t mask = (1ULL << bits) - 1;
    value &= mask;
    
    // Add to buffer
    buffer |= (value << bufferBits);
    bufferBits += bits;
    
    // Write full bytes
    while (bufferBits >= 8) {
        uint8_t byte = static_cast<uint8_t>(buffer & 0xFF);
        out.write(reinterpret_cast<const char*>(&byte), 1);
        buffer >>= 8;
        bufferBits -= 8;
    }
}

void BAMTagBinaryWriter::flushBuffer(std::ostream& out, uint64_t buffer, uint64_t bufferBits) {
    if (bufferBits == 0) return;
    
    // Pad to byte boundary and write
    while (bufferBits > 0) {
        uint8_t byte = static_cast<uint8_t>(buffer & 0xFF);
        out.write(reinterpret_cast<const char*>(&byte), 1);
        buffer >>= 8;
        bufferBits = (bufferBits >= 8) ? (bufferBits - 8) : 0;
    }
}
