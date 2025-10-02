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

static const bool g_debugTag = (std::getenv("STAR_DEBUG_TAG") != nullptr);

inline void logDebugTag(const std::string &msg) {
    if (!g_debugTag) return;
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    std::ostringstream thread_id; thread_id << std::this_thread::get_id();
    std::cerr << "[DEBUG_TAG " << std::put_time(&tm, "%H:%M:%S") << " T"
              << thread_id.str().substr(thread_id.str().length()-4) << "] " << msg << std::endl;
}

template <typename T>
inline void safe_push_back(std::vector<T> &vec, const T &value, const char *context) {
    try {
        vec.push_back(value);
        if (g_debugTag && vec.size() % 100000 == 0) {
            logDebugTag(std::string(context) + " size=" + std::to_string(vec.size()));
        }
    } catch (...) {
        std::cerr << "ERROR: Allocation failed in " << context << std::endl; std::exit(1);
    }
}

template <typename T>
inline void safe_reserve(std::vector<T> &vec, size_t capacity, const char *context) {
    try { vec.reserve(capacity); } catch (...) {
        std::cerr << "ERROR: reserve failed in " << context << std::endl; std::exit(1);
    }
}

BAMTagBuffer::BAMTagBuffer() {
    safe_reserve(entries, 1000000, "BAMTagBuffer::entries");
    if (g_debugTag) logDebugTag("BAMTagBuffer init with capacity 1M");
}

BAMTagBuffer::~BAMTagBuffer() { clear(); }

void BAMTagBuffer::append(const BAMRecordMeta& meta) {
    std::lock_guard<std::mutex> lock(entriesMutex);
    if (meta.recordIndex > UINT32_MAX || meta.iReadAll > UINT32_MAX) {
        std::cerr << "ERROR: index overflow in BAMTagBuffer::append" << std::endl; std::exit(1);
    }
    safe_push_back(entries, BAMTagEntry(meta), "BAMTagBuffer::entries");
}

#ifndef SOLO_USE_PACKED_READINFO
void BAMTagBuffer::writeTagBinary(const std::string& outputPath,
                                   const std::vector<readInfoStruct>& readInfo,
                                   uint64_t cbBits, uint64_t umiBits) {
    std::lock_guard<std::mutex> lock(entriesMutex);
    std::ofstream out(outputPath, std::ios::binary);
    if (!out.is_open()) { std::cerr << "ERROR: cannot open " << outputPath << std::endl; return; }
    size_t valid=0; for (auto &e: entries) if (e.readId < readInfo.size()) valid++;
    BAMTagBinaryWriter writer; TagStreamHeader header(cbBits, umiBits, valid); writer.writeHeader(out, header);
    std::sort(entries.begin(), entries.end(), [](auto&a,auto&b){return a.recordIndex<b.recordIndex;});
    for (auto &e: entries) {
        if (e.readId >= readInfo.size()) continue;
        const readInfoStruct &ri = readInfo[e.readId];
        uint8_t status = (ri.cb != (uint64)-1 && ri.umi != (uint32)-1) ? 1 : 0;
        uint64_t cbIdx = status ? (ri.cb + 1) : 0; // 0 reserved
        uint64_t umi = status ? ri.umi : 0;
        writer.writeRecord(out, status, cbIdx, umi, header);
    }
}
#else
void BAMTagBuffer::writeTagBinaryPacked(const std::string& outputPath,
                                        const PackedReadInfo& packed,
                                        uint64_t cbBits, uint64_t umiBits) {
    std::lock_guard<std::mutex> lock(entriesMutex);
    std::ofstream out(outputPath, std::ios::binary);
    if (!out.is_open()) { std::cerr << "ERROR: cannot open " << outputPath << std::endl; return; }
    size_t valid=0; for (auto &e: entries) if (e.readId < packed.data.size()) valid++;
    BAMTagBinaryWriter writer; TagStreamHeader header(cbBits, umiBits, valid); writer.writeHeader(out, header);
    std::sort(entries.begin(), entries.end(), [](auto&a,auto&b){return a.recordIndex<b.recordIndex;});
    for (auto &e: entries) {
        if (e.readId >= packed.data.size()) continue;
        uint32_t cb = packed.getCB(e.readId);
        uint32_t umi = packed.getUMI(e.readId);
        uint8_t status = packed.getStatus(e.readId); // map multi-state → 1 bit for now
        uint8_t statusBit = (status == 1) ? 1 : 0;
        uint64_t cbIdx = statusBit ? (uint64_t)cb : 0; // assume cb already 0-based + sentinel shift handled elsewhere
        writer.writeRecord(out, statusBit, cbIdx, umi, header);
    }
}
#endif

void BAMTagBuffer::clear() {
    std::lock_guard<std::mutex> lock(entriesMutex);
    entries.clear();
    entries.shrink_to_fit();
}
