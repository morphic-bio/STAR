#include "../source/SoloReadInfoLoader.h"
#include "../source/SoloReadInfoSink.h"
#include "../source/SoloFeature.h"
#include "../source/SoloReadFeature.h"
#include "../source/Parameters.h"
#include <cassert>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>

int main(){
    const char* skip = std::getenv("STAR_TEST_SKIP_TXOME");
    if (skip && std::strcmp(skip, "1")==0) {
        std::cout << "MinimalSink_test: skipped via STAR_TEST_SKIP_TXOME=1" << std::endl;
        return 0;
    }
    const char* txdir = std::getenv("STAR_TXOME_DIR");
    if (txdir==nullptr || access(txdir, R_OK)!=0) {
        std::cout << "MinimalSink_test: skipped (STAR_TXOME_DIR not set/readable)" << std::endl;
        return 0;
    }
    // Ensure required transcriptome metadata exists to avoid crashes
    {
        std::string gtab = std::string(txdir) + "/geneInfo.tab";
        if (access(gtab.c_str(), R_OK)!=0) {
            std::cout << "MinimalSink_test: skipped (" << gtab << " not found)" << std::endl;
            return 0;
        }
    }
    Parameters P;
    char templ[] = "/tmp/solo_sink_test.XXXXXXXX";
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

    SoloReadFeature rf(SoloFeatureTypes::GeneFull, P, 0);
    std::string path = P.outFileTmp + "/solo" + SoloFeatureTypes::Names[SoloFeatureTypes::GeneFull] + "_0";
    { std::ofstream ofs(path, std::ios::out | std::ios::trunc); ofs << "123 500 0 0\n"; ofs.flush(); }

    // Instantiate a SoloFeature to hold results (requires a Transcriptome). Use provided STAR_TXOME_DIR.
    P.quant.yes = true;
    // Point transcriptome loader to genomeDir directly (reads geneInfo.tab, etc.)
    P.pGe.sjdbGTFfile = "-";
    P.pGe.gDir = std::string(txdir);
    Transcriptome trans(P);
    SoloFeature *soloFeatAll[1] = {nullptr};
    SoloFeature feature(P, nullptr, trans, SoloFeatureTypes::GeneFull, nullptr, soloFeatAll);
    feature.resetPackedStorage(1);

    SoloReadInfoLoader loader;
    MinimalSink sink;
    std::vector<uint32_t> cbCounts(1,1), nU(1,0), nM(1,0);
    SoloReadFlagClass flags;

    loader.load(rf, SoloReadInfoMode::Minimal,
                [&](const ReadInfoRecord &rec){ sink.onRecord(feature, rec); },
                cbCounts, flags, nU, nM);
    sink.finalize(feature);

    assert(feature.getPackedUMI(0)==123u);
    assert(feature.getPackedCB(0)==0u);

    unlink(path.c_str());
    rmdir(P.outFileTmp.c_str());
    rmdir(baseDir.c_str());
    std::cout << "MinimalSink_test: OK" << std::endl;
    return 0;
}


