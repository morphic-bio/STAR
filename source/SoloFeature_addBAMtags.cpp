#include "SoloFeature.h"
#include "BAMfunctions.h"
#include "SequenceFuns.h"

void SoloFeature::addBAMtags(char *&bam0, uint32 &size0, char *bam1)
{//add extra tags to the BAM record
    
    uint64 iread = * ((uint64*) (bam0+size0));
    iread = iread >> 32; //iRead encoded in upper 32 bits

#ifndef SOLO_USE_PACKED_READINFO
    string cb="-";
    string umi="-";
    if (readInfo[iread].cb+1 != 0)             
        cb  = pSolo.cbWLstr[readInfo[iread].cb];
    if (readInfo[iread].umi+1 != 0)
        umi = convertNuclInt64toString(readInfo[iread].umi, pSolo.umiL);
#else
    uint32_t cbIdx = packedReadInfo.getCB((uint32_t)iread);
    uint32_t umiPacked = packedReadInfo.getUMI((uint32_t)iread);
    uint8_t status = packedReadInfo.getStatus((uint32_t)iread);
    string cb="-";
    string umi="-";
    if (status == 1 && cbIdx < pSolo.cbWLstr.size()) {
        cb = pSolo.cbWLstr[cbIdx];
    }
    if (status == 1) {
        umi.reserve(pSolo.umiL);
        uint32_t tmp = umiPacked;
        for (int i = pSolo.umiL - 1; i >= 0; --i) {
            uint8_t b = tmp & 0x3; tmp >>= 2;
            char c = 'N';
            switch (b) { case 0: c='A'; break; case 1: c='C'; break; case 2: c='G'; break; case 3: c='T'; }
            umi.push_back(c);
        }
        std::reverse(umi.begin(), umi.end());
    }
#endif

    memcpy(bam1, bam0, size0);

    size0 += bamAttrArrayWrite(cb,  "CB", bam1+size0);
    size0 += bamAttrArrayWrite(umi, "UB", bam1+size0);
    uint32 *bam1i = (uint32*) bam1;
    bam1i[0] = size0-sizeof(uint32);
    bam0=bam1;
};
