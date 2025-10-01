#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

// Simple decoder for binary tag stream format
// Header: 4x uint64_t (statusBits, cbBits, umiBits, recordCount)
// Records: bit-packed (status + CB index + packed UMI)

struct TagStreamHeader {
    uint64_t statusBits;
    uint64_t cbBits;
    uint64_t umiBits;
    uint64_t recordCount;
};

class BinaryTagDecoder {
public:
    static uint64_t readLittleEndian64(std::istream& in) {
        uint64_t value = 0;
        for (int i = 0; i < 8; i++) {
            uint8_t byte;
            in.read(reinterpret_cast<char*>(&byte), 1);
            if (in.fail()) return 0;
            value |= (static_cast<uint64_t>(byte) << (8 * i));
        }
        return value;
    }
    
    static uint64_t readBits(std::istream& in, uint64_t bits, uint64_t& buffer, uint64_t& bufferBits) {
        if (bits == 0) return 0;
        
        uint64_t result = 0;
        uint64_t resultBits = 0;
        
        while (resultBits < bits) {
            // Fill buffer if needed
            while (bufferBits < 8 && !in.eof()) {
                uint8_t byte;
                in.read(reinterpret_cast<char*>(&byte), 1);
                if (in.fail()) break;
                buffer |= (static_cast<uint64_t>(byte) << bufferBits);
                bufferBits += 8;
            }
            
            if (bufferBits == 0) break; // EOF
            
            // Extract bits
            uint64_t needed = bits - resultBits;
            uint64_t available = std::min(needed, bufferBits);
            uint64_t mask = (1ULL << available) - 1;
            
            result |= ((buffer & mask) << resultBits);
            buffer >>= available;
            bufferBits -= available;
            resultBits += available;
        }
        
        return result;
    }
    
    static std::string unpackUmi(uint64_t packed, size_t umiLength) {
        if (packed == 0) return "-";
        
        std::string umi;
        umi.reserve(umiLength);
        
        const char bases[] = {'A', 'C', 'G', 'T'};
        // Read bits from high order end (MSB to LSB) to match STAR's packing order
        for (size_t i = 0; i < umiLength; ++i) {
            uint64_t base = (packed >> (2 * (umiLength - 1 - i))) & 3;
            umi += bases[base];
        }
        
        return umi;
    }
};

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <binary_file> <cb_whitelist> <umi_length>" << std::endl;
        return 1;
    }
    
    std::string binaryFile = argv[1];
    std::string whitelistFile = argv[2];
    size_t umiLength = std::stoul(argv[3]);
    
    // Read CB whitelist
    std::vector<std::string> cbWhitelist;
    std::ifstream wlFile(whitelistFile);
    std::string cb;
    while (std::getline(wlFile, cb)) {
        cbWhitelist.push_back(cb);
    }
    wlFile.close();
    
    // Open binary file
    std::ifstream binFile(binaryFile, std::ios::binary);
    if (!binFile.is_open()) {
        std::cerr << "ERROR: Cannot open binary file: " << binaryFile << std::endl;
        return 1;
    }
    
    // Read header (v2 format: 4x uint64_t)
    TagStreamHeader header;
    header.statusBits = BinaryTagDecoder::readLittleEndian64(binFile);
    header.cbBits = BinaryTagDecoder::readLittleEndian64(binFile);
    header.umiBits = BinaryTagDecoder::readLittleEndian64(binFile);
    header.recordCount = BinaryTagDecoder::readLittleEndian64(binFile);
    
    std::cerr << "Header: statusBits=" << header.statusBits 
              << ", cbBits=" << header.cbBits 
              << ", umiBits=" << header.umiBits 
              << ", recordCount=" << header.recordCount << std::endl;
    
    std::cout << "# bam_record_index\tiReadAll\tmate\talign_idx\tCB\tUB\tstatus" << std::endl;
    
    // Read records
    uint64_t buffer = 0;
    uint64_t bufferBits = 0;
    uint64_t recordIndex = 0;
    
    while (recordIndex < header.recordCount && !binFile.eof()) {
        // Read status bit
        uint64_t status = BinaryTagDecoder::readBits(binFile, header.statusBits, buffer, bufferBits);
        if (binFile.eof() && bufferBits == 0) break;
        
        // Read CB index
        uint64_t cbIdx = BinaryTagDecoder::readBits(binFile, header.cbBits, buffer, bufferBits);
        
        // Read UMI packed
        uint64_t umiPacked = BinaryTagDecoder::readBits(binFile, header.umiBits, buffer, bufferBits);
        
        // Convert to strings
        std::string cbStr = "-";
        if (status && cbIdx > 0 && (cbIdx - 1) < cbWhitelist.size()) {
            cbStr = cbWhitelist[cbIdx - 1]; // -1 because we stored +1
        }
        
        std::string ubStr = "-";
        if (status && umiPacked > 0) {
            ubStr = BinaryTagDecoder::unpackUmi(umiPacked, umiLength);
        }
        
        std::string statusStr = status ? "OK" : "NO_CB_UMI";
        
        // Output TSV format (simplified - missing some columns for compatibility)
        std::cout << recordIndex << "\t"
                  << recordIndex << "\t"  // Use recordIndex as readId placeholder
                  << "\t"                 // mate (empty)
                  << "0\t"               // alignIdx placeholder
                  << cbStr << "\t"
                  << ubStr << "\t"
                  << statusStr << std::endl;
        
        recordIndex++;
    }
    
    binFile.close();
    return 0;
}
