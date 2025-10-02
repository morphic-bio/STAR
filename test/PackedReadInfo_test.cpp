#include "../source/PackedReadInfo.h"
#include <cassert>
#include <stdexcept>

int main() {
    // Basic init
    PackedReadInfo pri;
    pri.init(10, /*wlSize*/ 1024, /*umiLength*/ 12);

    // Round-trip pack/unpack
    pri.set(0, /*cb*/ 123, /*umi*/ 0xABCDE, /*status*/ 1);
    assert(pri.getCB(0) == 123);
    assert(pri.getUMI(0) == 0xABCDE);
    assert(pri.getStatus(0) == 1);

    // Zero-UMI should be allowed
    PackedReadInfo pri0;
    pri0.init(5, 1024, 0);
    pri0.set(0, 55, 0, 1);
    assert(pri0.getCB(0) == 55);
    assert(pri0.getUMI(0) == 0);
    assert(pri0.getStatus(0) == 1);

    // Overflow guards: use a tiny whitelist and UMI bits
    PackedReadInfo small;
    small.init(1, 1, 1); // cbBits>=1, umiBits=2
    bool threw=false;
    try { small.set(0, /*cbIdx*/ 9999, /*umi*/ 0, /*status*/ 1); } catch (...) { threw=true; }
    assert(threw);

    threw=false;
    try { small.set(0, /*cbIdx*/ 1, /*umi*/ 0xFFFF, /*status*/ 1); } catch (...) { threw=true; }
    assert(threw);

    return 0;
}


