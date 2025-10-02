#include "PackedReadInfo.h"

uint8_t log2ceil_u32(uint32_t n) {
    if (n <= 1) return 1;
    uint32_t x = n - 1;
    uint8_t bits = 0;
    while (x) { x >>= 1; ++bits; }
    return bits;
}

void PackedReadInfo::init(uint32_t nReads, uint32_t wlSize, uint8_t umiLength) {
    L.whitelistSize = wlSize;
    L.umiLength = umiLength;
    L.umiBits = umiLength * 2;
    if (L.umiBits == 0 || L.umiBits > 32)
        throw std::runtime_error("PackedReadInfo: UMI length unsupported");
    L.cbBits = log2ceil_u32(wlSize + 1);
    if (L.cbBits > 32)
        throw std::runtime_error("PackedReadInfo: whitelist too large");
    L.statusBits = 3;
    L.reserved = 0;
    if (L.umiBits + L.cbBits + L.statusBits > 64)
        throw std::runtime_error("PackedReadInfo: bit overflow");
    data.assign(nReads, 0ULL);
}

inline uint64_t PackedReadInfo::maskUMI() const {
    return (L.umiBits==64)?~0ULL:((1ULL<<L.umiBits)-1ULL);
}

inline uint64_t PackedReadInfo::maskCB() const {
    return (L.cbBits==64)?~0ULL:((1ULL<<L.cbBits)-1ULL);
}

inline uint64_t PackedReadInfo::pack(uint32_t cbIdx, uint32_t umiPacked, uint8_t status) const {
    uint64_t w=0ULL;
    w |= (uint64_t)(umiPacked & maskUMI());
    w |= (uint64_t)(cbIdx & maskCB()) << L.umiBits;
    w |= (uint64_t)(status & ((1u<<L.statusBits)-1u)) << (L.umiBits + L.cbBits);
    return w;
}

inline void PackedReadInfo::unpack(uint64_t w, uint32_t &cbIdx, uint32_t &umiPacked, uint8_t &status) const {
    umiPacked = (uint32_t)(w & maskUMI());
    cbIdx = (uint32_t)((w >> L.umiBits) & maskCB());
    status = (uint8_t)((w >> (L.umiBits + L.cbBits)) & ((1u<<L.statusBits)-1u));
}

inline void PackedReadInfo::set(uint32_t readId, uint32_t cbIdx, uint32_t umiPacked, uint8_t status) {
    data[readId] = pack(cbIdx, umiPacked, status);
}

inline uint32_t PackedReadInfo::getCB(uint32_t readId) const {
    uint64_t w = data[readId];
    return (uint32_t)((w >> L.umiBits) & maskCB());
}

inline uint32_t PackedReadInfo::getUMI(uint32_t readId) const {
    uint64_t w = data[readId];
    return (uint32_t)(w & maskUMI());
}

inline uint8_t PackedReadInfo::getStatus(uint32_t readId) const {
    uint64_t w = data[readId];
    return (uint8_t)((w >> (L.umiBits + L.cbBits)) & ((1u<<L.statusBits)-1u));
}

inline void PackedReadInfo::setStatus(uint32_t readId, uint8_t status) {
    uint64_t w = data[readId];
    uint64_t mask = ((1ULL<<L.statusBits)-1ULL) << (L.umiBits + L.cbBits);
    w &= ~mask;
    w |= (uint64_t)(status & ((1u<<L.statusBits)-1u)) << (L.umiBits + L.cbBits);
    data[readId] = w;
}
