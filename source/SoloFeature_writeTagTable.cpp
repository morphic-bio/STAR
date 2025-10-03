#include "SoloFeature.h"
#include "BAMTagBuffer.h"
#include "BAMTagBinaryWriter.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include <algorithm>

void SoloFeature::writeTagTableIfRequested(bool filteredPass)
{
    if (!pSolo.writeTagTableEnabled || featureType != pSolo.samAttrFeature) {
        return;
    }
    
    if (!pSolo.bamTagBuffer) {
        P.inOut->logMain << "WARNING: Tag table export requested but BAMTagBuffer is not available" << endl;
        return;
    }
    
    time_t rawTime;
    time(&rawTime);
    
    uint64_t cbBits = std::max<uint64_t>(1, BAMTagBinaryWriter::integerLog2ceil(pSolo.cbWLstr.size() + 1));
    uint64_t umiBits = pSolo.umiL * 2; // 2 bits per base
    if (umiBits > 32) {
        P.inOut->logMain << "WARNING: UMI length " << pSolo.umiL << " requires " << umiBits
                         << " bits, truncating to 32" << endl;
        umiBits = 32;
    }
    
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Writing binary tag stream to " << pSolo.writeTagTablePath
                     << " (32-byte header + records, CB:" << cbBits << " bits, UMI:" << umiBits << " bits)" << endl;

    pSolo.bamTagBuffer->writeTagBinaryPacked(pSolo.writeTagTablePath, packedReadInfo, cbBits, umiBits);

    pSolo.bamTagBuffer->clear();
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Finished writing binary tag stream and cleared buffer" << endl;
}
