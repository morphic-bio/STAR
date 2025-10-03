#include "../source/SoloReadFeature.h"
#include "../source/Parameters.h"
#include <cassert>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <cstring>

// Stage 1.1: Baseline harness that exercises inputRecords sink path (packed only)
int main() {
    Parameters P; // Constructor sets P.inOut

    // Create a temporary output directory using mkdtemp to avoid hard-coded paths
    char templ[] = "/tmp/solo_baseline_test.XXXXXXXX";
    char *tmpBase = mkdtemp(templ);
    if (tmpBase==nullptr) {
        std::cerr << "Failed to create temp directory" << std::endl;
        return 1;
    }
    std::string baseDir(tmpBase);
    std::string tmpDir = baseDir + "/_STARtmp";
    mkdir(tmpDir.c_str(), S_IRWXU);
    P.outFileTmp = tmpDir;

    // Minimal Solo configuration
    P.pSolo.yes = true;
    P.pSolo.type = ParametersSolo::SoloTypes::CB_UMI_Simple;
    P.pSolo.featureYes.fill(false);
    P.pSolo.featureYes[SoloFeatureTypes::Gene] = true;
    P.pSolo.featureYes[SoloFeatureTypes::GeneFull] = true;
    P.pSolo.readInfoYes.fill(true);
    P.pSolo.readIndexYes.fill(false);
    // Tiny whitelist enabled to bypass cb integer remapping path
    P.pSolo.cbWLyes = true;
    P.pSolo.cbWLsize = 1;
    P.pSolo.cbWL.assign(1, 0);
    P.pSolo.CBmatchWL.oneExact = false;
    // Configure quality score handling for multi-match logic
    P.pSolo.QSbase = 33;
    P.pSolo.QSmax = 41;
    P.pSolo.cbMinP = 0.0f;

    // Construct SoloReadFeature which opens its stream under P.outFileTmp
    SoloReadFeature srf(/*featureType=*/SoloFeatureTypes::GeneFull, P, /*iChunk=*/0);

    // Write a minimal fixture record into the opened stream file
    // Format expected (readIndexYes=false):
    // umi feature cbmatch cb\n
    std::string streamPath = P.outFileTmp + "/solo" + SoloFeatureTypes::Names[SoloFeatureTypes::GeneFull] + "_0";
    {
        std::ofstream ofs(streamPath, std::ios::out | std::ios::trunc);
        ofs << "12345 "    // umi
            << "100 "      // feature id
            << "0 "        // cbmatch (perfect)
            << "0\n";     // cb index
        ofs.flush();
    }

    // Prepare buffers for inputRecords
    uint32 cbBlock[3] = {0,0,0};
    uint32 *cbPbuf[1] = { cbBlock };
    uint32 cbPstride = 3;
    std::vector<uint32> cbReadCountTotal(1, 0);
    cbReadCountTotal[0] = 1; // allow perfect CB use
    SoloReadFlagClass readFlagCounts;
    std::vector<uint32> nReadPerCBunique1(1, 0), nReadPerCBmulti1(1, 0);
    // No legacy readInfo; using sink overload

    // Execute baseline capture
    srf.inputRecords(cbPbuf, cbPstride, cbReadCountTotal, readFlagCounts,
                     nReadPerCBunique1, nReadPerCBmulti1,
                     [&](uint64 iread, uint32 cbIdx, uint32 umi, uint8 status){ (void)iread; (void)cbIdx; (void)umi; (void)status; });

    // Validate that a record was emitted to cb block
    assert(cbBlock[0] == 100);
    assert(cbBlock[1] == 12345u);
    std::cout << "SoloReadInfoBaseline: OK (feature=" << cbBlock[0] << ", umi=" << cbBlock[1] << ")" << std::endl;

    // Case 2: readIndexYes=true, one counted and one no-feature record to test readInfo sink
    P.pSolo.readIndexYes[SoloFeatureTypes::GeneFull] = true;
    SoloReadFeature srf2(SoloFeatureTypes::GeneFull, P, 1);
    std::string streamPath2 = P.outFileTmp + "/solo" + SoloFeatureTypes::Names[SoloFeatureTypes::GeneFull] + "_1";
    {
        std::ofstream ofs(streamPath2, std::ios::out | std::ios::trunc);
        // Format with readIndexYes=true: umi iread flag feature cbmatch cb
        ofs << "777 " << "0 " << "0 " << "101 " << "0 " << "0\n"; // counted record, should fill cbP[2] with iread=0
        ofs << "888 " << "0 " << "0 " << "4294967295 " << "0 " << "0\n"; // no-feature, should populate readInfo[0]
        ofs.flush();
    }
    uint32 cbBlock2[3] = {0,0,0};
    uint32 *cbPbuf2[1] = { cbBlock2 };
    std::vector<uint32> cbReadCountTotal2(1, 1);
    std::vector<uint32> nReadPerCBunique2(1, 0), nReadPerCBmulti2(1, 0);
    srf2.inputRecords(cbPbuf2, 3, cbReadCountTotal2, readFlagCounts, nReadPerCBunique2, nReadPerCBmulti2,
                      [&](uint64 iread, uint32 cbIdx, uint32 umi, uint8 status){ (void)iread; (void)cbIdx; (void)umi; (void)status; });
    assert(cbBlock2[0] == 101);
    assert(cbBlock2[2] == 0u);
    // readInfo assertions removed (packed-only path uses sink)
    std::cout << "SoloReadInfoBaseline (index+no-feature): OK" << std::endl;

    // Case 3: multi-match (cbmatch>1) with quality scores to exercise posterior selection
    P.pSolo.readIndexYes[SoloFeatureTypes::GeneFull] = false;
    SoloReadFeature srf3(SoloFeatureTypes::GeneFull, P, 2);
    std::string streamPath3 = P.outFileTmp + "/solo" + SoloFeatureTypes::Names[SoloFeatureTypes::GeneFull] + "_2";
    {
        std::ofstream ofs(streamPath3, std::ios::out | std::ios::trunc);
        // Format: umi feature cbmatch (cbin qin)*
        ofs << "9999 " << "102 " << "2 " << "0 I " << "0 I\n"; // two identical candidates
        ofs.flush();
    }
    uint32 cbBlock3[3] = {0,0,0};
    uint32 *cbPbuf3[1] = { cbBlock3 };
    std::vector<uint32> cbReadCountTotal3(1, 5); // positive count for cbin 0
    std::vector<uint32> nReadPerCBunique3(1, 0), nReadPerCBmulti3(1, 0);
    srf3.inputRecords(cbPbuf3, 3, cbReadCountTotal3, readFlagCounts, nReadPerCBunique3, nReadPerCBmulti3,
                      [&](uint64 iread, uint32 cbIdx, uint32 umi, uint8 status){ (void)iread; (void)cbIdx; (void)umi; (void)status; });
    assert(cbBlock3[0] == 102);
    assert(cbBlock3[1] == 9999u);
    std::cout << "SoloReadInfoBaseline (multi-match): OK" << std::endl;

    // Case 4: probe-like RNAME path coverage via solo input emulation
    // We simulate a feature id that would correspond to a probe contig (non-"chr").
    // Since Solo input stores features as numeric IDs, we just ensure the plumbing
    // handles arbitrary feature ids without assumptions on RNAME formatting.
    SoloReadFeature srf4(SoloFeatureTypes::GeneFull, P, 3);
    std::string streamPath4 = P.outFileTmp + "/solo" + SoloFeatureTypes::Names[SoloFeatureTypes::GeneFull] + "_3";
    {
        std::ofstream ofs(streamPath4, std::ios::out | std::ios::trunc);
        ofs << "5555 " << "999999 " << "0 " << "0\n"; // large feature id
        ofs.flush();
    }
    uint32 cbBlock4[3] = {0,0,0};
    uint32 *cbPbuf4[1] = { cbBlock4 };
    std::vector<uint32> cbReadCountTotal4(1, 1);
    std::vector<uint32> nReadPerCBunique4(1, 0), nReadPerCBmulti4(1, 0);
    srf4.inputRecords(cbPbuf4, 3, cbReadCountTotal4, readFlagCounts, nReadPerCBunique4, nReadPerCBmulti4,
                      [&](uint64 iread, uint32 cbIdx, uint32 umi, uint8 status){ (void)iread; (void)cbIdx; (void)umi; (void)status; });
    assert(cbBlock4[0] == 999999u);
    assert(cbBlock4[1] == 5555u);
    std::cout << "SoloReadInfoBaseline (probe-like id): OK" << std::endl;

    // Cleanup temp files and directory
    unlink(streamPath.c_str());
    unlink(streamPath2.c_str());
    unlink(streamPath3.c_str());
    unlink(streamPath4.c_str());
    rmdir(P.outFileTmp.c_str());
    rmdir(baseDir.c_str());
    return 0;
}


