#include "BAMTagBuffer.h"
#include "BAMTagBinaryWriter.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"
#include <fstream>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <thread>

// Debug flag - controlled by STAR_DEBUG_TAG environment variable
static const bool g_debugTag = (std::getenv("STAR_DEBUG_TAG") != nullptr);

// Debug logging helper with timestamp and thread ID
inline void logDebugTag(const std::string &msg) {
    if (!g_debugTag) return;
    
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    std::ostringstream thread_id;
    thread_id << std::this_thread::get_id();
    
    std::cerr << "[DEBUG_TAG " << std::put_time(&tm, "%H:%M:%S") 
              << " T" << thread_id.str().substr(thread_id.str().length()-4) << "] " 
              << msg << std::endl;
}

// Safe allocation helper for vectors
template <typename T>
inline void safe_push_back(std::vector<T> &vec, const T &value, const char *context) {
    try {
        vec.push_back(value);
        if (g_debugTag && vec.size() % 100000 == 0) {
            logDebugTag("Vector " + std::string(context) + " size: " + std::to_string(vec.size()));
        }
    } catch (const std::bad_alloc &e) {
        std::cerr << "ERROR: Allocation failed in " << context
                  << " (current size=" << vec.size() << ")" << std::endl;
        std::cerr << "Consider reducing --runThreadN or checking system memory." << std::endl;
        std::cerr << "Exception: " << e.what() << std::endl;
        std::exit(1);
    } catch (const std::exception &e) {
        std::cerr << "ERROR: Exception in " << context 
                  << " (current size=" << vec.size() << "): " << e.what() << std::endl;
        std::exit(1);
    }
}

// Safe reserve helper for vectors
template <typename T>
inline void safe_reserve(std::vector<T> &vec, size_t capacity, const char *context) {
    try {
        vec.reserve(capacity);
        if (g_debugTag) {
            logDebugTag("Reserved " + std::to_string(capacity) + " entries for " + std::string(context));
        }
    } catch (const std::bad_alloc &e) {
        std::cerr << "ERROR: Reserve failed in " << context
                  << " (requested=" << capacity << ", current size=" << vec.size() << ")" << std::endl;
        std::cerr << "Consider reducing dataset size or checking system memory." << std::endl;
        std::cerr << "Exception: " << e.what() << std::endl;
        std::exit(1);
    }
}

BAMTagBuffer::BAMTagBuffer() {
    // Reserve some initial capacity to reduce reallocations
    safe_reserve(entries, 1000000, "BAMTagBuffer::entries initial allocation"); // Start with 1M entries capacity
    
    if (g_debugTag) {
        logDebugTag("BAMTagBuffer initialized with 1M entries capacity");
    }
}

BAMTagBuffer::~BAMTagBuffer() {
    clear();
}

void BAMTagBuffer::append(const BAMRecordMeta& meta) {
    std::lock_guard<std::mutex> lock(entriesMutex);
    
    if (g_debugTag && entries.size() == 0) {
        logDebugTag("First BAMTagBuffer::append call - recordIndex=" + std::to_string(meta.recordIndex) + 
                   ", readId=" + std::to_string(meta.iReadAll));
    }
    
    // Guard against overflow of recordIndex field
    if (meta.recordIndex > UINT32_MAX) {
        std::cerr << "ERROR: Record index " << meta.recordIndex << " exceeds UINT32_MAX limit. "
                  << "Dataset too large for current implementation." << std::endl;
        exit(1);
    }
    
    // Guard against overflow of readId field  
    if (meta.iReadAll > UINT32_MAX) {
        std::cerr << "ERROR: Read ID " << meta.iReadAll << " exceeds UINT32_MAX limit. "
                  << "Dataset too large for current implementation." << std::endl;
        exit(1);
    }
    
    // Create entry with recordIndex and readId
    BAMTagEntry entry(meta);
    safe_push_back(entries, entry, "BAMTagBuffer::entries");
}



void BAMTagBuffer::writeTagBinary(const std::string& outputPath,
                                   const std::vector<readInfoStruct>& readInfo,
                                   uint64_t cbBits,
                                   uint64_t umiBits) {
    std::lock_guard<std::mutex> lock(entriesMutex);
    
    std::ofstream outFile(outputPath, std::ios::binary);
    if (!outFile.is_open()) {
        std::cerr << "ERROR: Cannot open binary tag output file: " << outputPath << std::endl;
        return;
    }
    
    // First pass: count valid records
    size_t validRecords = 0;
    for (const auto& entry : entries) {
        if (entry.readId < readInfo.size()) {
            validRecords++;
        }
    }
    
    if (g_debugTag) {
        logDebugTag("writeTagBinary: total entries=" + std::to_string(entries.size()) + 
                   ", valid records=" + std::to_string(validRecords) + 
                   ", readInfo size=" + std::to_string(readInfo.size()));
        
        // Assert valid records should match entries when all readIds are valid
        if (validRecords != entries.size()) {
            logDebugTag("WARNING: Not all entries have valid readIds - " + 
                       std::to_string(entries.size() - validRecords) + " entries will be skipped");
        }
    }
    
    // Create binary writer and header with record count
    BAMTagBinaryWriter writer;
    TagStreamHeader header(cbBits, umiBits, validRecords); // statusBits=1 (default)
    
    if (g_debugTag) {
        logDebugTag("Binary tag header: statusBits=1, cbBits=" + std::to_string(cbBits) + 
                   ", umiBits=" + std::to_string(umiBits) + ", recordCount=" + std::to_string(validRecords));
    }
    
    // Write binary header (4x uint64_t)
    writer.writeHeader(outFile, header);
    
    // Sort entries by recordIndex to restore BAM record order
    std::sort(entries.begin(), entries.end(), 
              [](const BAMTagEntry& a, const BAMTagEntry& b) {
                  return a.recordIndex < b.recordIndex;
              });
    
    // Write binary records in BAM order
    size_t recordsWritten = 0;
    size_t skippedRecords = 0;
    uint32_t lastRecordIndex = 0;
    
    for (size_t i = 0; i < entries.size(); i++) {
        const auto& entry = entries[i];
        
        if (entry.readId >= readInfo.size()) {
            skippedRecords++;
            if (g_debugTag && skippedRecords <= 5) {
                logDebugTag("Skipping entry " + std::to_string(i) + " with invalid readId=" + 
                           std::to_string(entry.readId) + " (readInfo.size=" + std::to_string(readInfo.size()) + ")");
            }
            continue; // Skip invalid entries
        }
        
        // Phase 5: Check record index monotonicity
        if (g_debugTag && recordsWritten > 0 && entry.recordIndex < lastRecordIndex) {
            logDebugTag("WARNING: Record index out of order at entry " + std::to_string(i) + 
                       ": current=" + std::to_string(entry.recordIndex) + 
                       ", previous=" + std::to_string(lastRecordIndex));
        }
        lastRecordIndex = entry.recordIndex;
        
        const readInfoStruct& readData = readInfo[entry.readId];
        
        // Compute binary status: 1 if both CB and UMI are valid, 0 otherwise
        bool hasCB = (readData.cb != (uint64)-1);
        bool hasUMI = (readData.umi != (uint32)-1);
        uint8_t status = (hasCB && hasUMI) ? 1 : 0;
        
        // Compute CB index: shift by +1 if valid (0 reserved for missing), 0 if invalid/missing
        uint64_t cbIdx = 0;
        if (status && hasCB) {
            cbIdx = readData.cb + 1; // +1 because 0 is reserved for missing
        }
        
        // Use packed UMI directly from readInfo (0 if invalid/missing)
        uint64_t umiPacked = 0;
        if (status && hasUMI) {
            umiPacked = readData.umi;
        }
        
        // Phase 3: Validate record data
        if (g_debugTag) {
            if (status == 1 && cbIdx == 0) {
                logDebugTag("ERROR: Record " + std::to_string(recordsWritten) + 
                           " has status=OK but cbIdx=0 (reserved for missing)");
            }
            if (recordsWritten < 5) {
                logDebugTag("Record " + std::to_string(recordsWritten) + 
                           ": status=" + std::to_string(status) + 
                           ", cbIdx=" + std::to_string(cbIdx) + 
                           ", umiPacked=" + std::to_string(umiPacked) + 
                           ", recordIndex=" + std::to_string(entry.recordIndex));
            }
        }
        
        // Write binary record
        writer.writeRecord(outFile, status, cbIdx, umiPacked, header);
        recordsWritten++;
    }
    
    if (g_debugTag && skippedRecords > 0) {
        logDebugTag("Skipped " + std::to_string(skippedRecords) + " entries with invalid readIds");
    }
    
    // Sanity check
    if (recordsWritten != validRecords) {
        std::cerr << "WARNING: Expected to write " << validRecords 
                  << " records but actually wrote " << recordsWritten << std::endl;
    }
    
    // Phase 5: Verify final record index
    if (g_debugTag && recordsWritten > 0) {
        // Note: lastRecordIndex is the actual BAM record index, not necessarily recordsWritten-1
        logDebugTag("Final record index: " + std::to_string(lastRecordIndex) + 
                   ", total records written: " + std::to_string(recordsWritten));
    }
    
    outFile.close();
    
    // Get file size for reporting and Phase 3 validation
    std::ifstream sizeFile(outputPath, std::ios::binary | std::ios::ate);
    size_t fileSize = sizeFile.tellg();
    sizeFile.close();
    
    // Phase 3: Verify file size matches expected format
    if (g_debugTag) {
        size_t recordBits = 1 + cbBits + umiBits; // status + CB + UMI bits
        size_t recordBytes = (recordBits + 7) / 8; // Round up to bytes
        size_t expectedSize = 32 + (recordsWritten * recordBytes); // 32-byte header + records
        
        logDebugTag("File size validation: actual=" + std::to_string(fileSize) + 
                   " bytes, expected~" + std::to_string(expectedSize) + 
                   " bytes (32-byte header + " + std::to_string(recordsWritten) + 
                   " records * " + std::to_string(recordBytes) + " bytes each)");
                   
        if (fileSize < 32) {
            logDebugTag("ERROR: File size too small - missing header");
        }
    }
    
    std::cout << "Binary tag stream written to " << outputPath 
              << " with " << validRecords << " records"
              << " (" << fileSize << " bytes, 32-byte header + records, " 
              << "CB:" << cbBits << " bits, UMI:" << umiBits << " bits)" << std::endl;
}

void BAMTagBuffer::clear() {
    std::lock_guard<std::mutex> lock(entriesMutex);
    entries.clear();
    entries.shrink_to_fit();
}
