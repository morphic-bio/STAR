#include "../source/SoloReadInfoLoader.h"
#include "../source/SoloReadInfoSink.h"
#include "../source/SoloFeature.h"
#include "../source/SoloReadFeature.h"
#include "../source/Parameters.h"
#include <cassert>
#include <sys/stat.h>
#include <unistd.h>
#include <map>
#include <tuple>
#include <iostream>

using RecTuple = std::tuple<uint32_t,uint32_t,uint8_t>; // cb, umi, status

int main(){
    Parameters P;
    char templ[] = "/tmp/solo_parity_test.XXXXXXXX";
    char *tmpBase = mkdtemp(templ);
    if (tmpBase==nullptr) return 1;
    std::string baseDir(tmpBase);
    std::string tmpDir = baseDir + "/_STARtmp";
    mkdir(tmpDir.c_str(), S_IRWXU);
    P.outFileTmp = tmpDir;

    P.pSolo.yes = true;
    P.pSolo.type = ParametersSolo::SoloTypes::CB_UMI_Simple;
    P.pSolo.featureYes.fill(false);
    P.pSolo.featureYes[SoloFeatureTypes::GeneFull] = true;
    P.pSolo.readInfoYes.fill(true);
    P.pSolo.readIndexYes.fill(false);
    P.pSolo.cbWLyes = true; P.pSolo.cbWLsize = 1; P.pSolo.cbWL.assign(1,0);
    P.pSolo.QSbase = 33; P.pSolo.QSmax = 41; P.pSolo.cbMinP = 0.0f;

    // SoloReadFeature with a small synthetic stream
    SoloReadFeature rf(SoloFeatureTypes::GeneFull, P, 0);
    std::string path = P.outFileTmp + "/solo" + SoloFeatureTypes::Names[SoloFeatureTypes::GeneFull] + "_0";
    {
        // Format: UMI feature cbmatch cb
        // Three records: two counted (feature!=UINT32_MAX) and one no-feature
        std::ofstream ofs(path, std::ios::out | std::ios::trunc);
        ofs << "111 100 0 0\n";
        ofs << "222 101 0 0\n";
        ofs << "333 4294967295 0 0\n"; // no-feature
        ofs.flush();
    }

    // SoloFeature to carry effects (not strictly required for loader equivalence)
    // Minimal transcriptome usage not needed since loader only reads temp stream here
    Transcriptome *transPtr = (Transcriptome*)nullptr;
    SoloFeature *soloFeatAll[1] = {nullptr};
    // Avoid dereferencing null transcriptome by using a dummy instance-less SoloFeature interface: not used in loader
    // Construct a local SoloFeature with a fake Transcriptome reference is not safe; we won't call into it.
    // We'll only use SoloReadFeature rf with loader below and ignore sinks' side effects.

    SoloReadInfoLoader loader;
    std::vector<uint32_t> cbCounts(1,1), nU(1,0), nM(1,0);
    SoloReadFlagClass flags;

    std::map<uint64_t, RecTuple> minimalByRead;
    std::map<uint64_t, RecTuple> countingByRead;

    // Minimal mode capture
    loader.load(rf, SoloReadInfoMode::Minimal,
                [&](const ReadInfoRecord &rec){ minimalByRead[rec.readId] = RecTuple{rec.cbIdx, rec.umi, rec.status}; },
                cbCounts, flags, nU, nM);

    // Counting mode capture; require good feature to match CountingSink semantics
    loader.load(rf, SoloReadInfoMode::Counting,
                [&](const ReadInfoRecord &rec){ if (rec.featureId!=(uint32_t)-1 && rec.status==1) countingByRead[rec.readId] = RecTuple{rec.cbIdx, rec.umi, rec.status}; },
                cbCounts, flags, nU, nM);

    // Parity: counting emissions must be a subset of minimal; tuples must match
    for (const auto &kv : countingByRead) {
        auto it = minimalByRead.find(kv.first);
        assert(it != minimalByRead.end());
        assert(kv.second == it->second);
    }

    unlink(path.c_str());
    rmdir(P.outFileTmp.c_str());
    rmdir(baseDir.c_str());
    std::cout << "SoloReadInfoParity_test: OK" << std::endl;
    return 0;
}


