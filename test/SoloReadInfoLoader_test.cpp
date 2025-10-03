#include "../source/SoloReadInfoLoader.h"
#include "../source/SoloReadFeature.h"
#include "../source/Parameters.h"
#include <cassert>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>

int main() {
    Parameters P;
    char templ[] = "/tmp/solo_loader_test.XXXXXXXX";
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

    // Case A: simple single-match, no read index
    SoloReadFeature rf(SoloFeatureTypes::GeneFull, P, 0);
    std::string streamPath = P.outFileTmp + "/solo" + SoloFeatureTypes::Names[SoloFeatureTypes::GeneFull] + "_0";
    { std::ofstream ofs(streamPath, std::ios::out | std::ios::trunc); ofs << "111 200 0 0\n"; ofs.flush(); }
    SoloReadInfoLoader loader;
    std::vector<uint32_t> cbReadCountTotal(1,1), nU(1,0), nM(1,0);
    SoloReadFlagClass flags;
    ReadInfoRecord captured{}; bool got=false;
    loader.load(rf, SoloReadInfoMode::Minimal, [&](const ReadInfoRecord &rec){ captured=rec; got=true; }, cbReadCountTotal, flags, nU, nM);
    assert(got && captured.cbIdx==0u && captured.umi==111u && captured.featureId==200u && captured.readIndex==(uint32_t)-1 && captured.status==1);

    unlink(streamPath.c_str());

    // Case B: readIndexYes=true with counted and no-feature
    P.pSolo.readIndexYes[SoloFeatureTypes::GeneFull] = true;
    SoloReadFeature rf2(SoloFeatureTypes::GeneFull, P, 1);
    std::string streamPath2 = P.outFileTmp + "/solo" + SoloFeatureTypes::Names[SoloFeatureTypes::GeneFull] + "_1";
    { std::ofstream ofs(streamPath2, std::ios::out | std::ios::trunc); ofs << "777 0 0 101 0 0\n" << "888 0 0 4294967295 0 0\n"; ofs.flush(); }
    ReadInfoRecord rec2a{}, rec2b{}; int seen=0;
    std::vector<uint32_t> cbCounts2(1,1);
    loader.loadMinimal(rf2, [&](const ReadInfoRecord &rec){ if(seen==0) rec2a=rec; else rec2b=rec; seen++; }, cbCounts2);
    assert(seen==2);
    assert(rec2a.status==1 && rec2a.cbIdx==0 && rec2a.umi==777 && rec2a.featureId==101 && rec2a.readIndex==0);
    assert(rec2b.status==1 && rec2b.cbIdx==0 && rec2b.umi==888 && rec2b.featureId==4294967295u && rec2b.readIndex==0);
    unlink(streamPath2.c_str());

    // Case C: multi-match failure emits missing-CB status (0)
    P.pSolo.readIndexYes[SoloFeatureTypes::GeneFull] = false;
    SoloReadFeature rf3(SoloFeatureTypes::GeneFull, P, 2);
    std::string streamPath3 = P.outFileTmp + "/solo" + SoloFeatureTypes::Names[SoloFeatureTypes::GeneFull] + "_2";
    { std::ofstream ofs(streamPath3, std::ios::out | std::ios::trunc); ofs << "9999 102 2 0 I 0 I\n"; ofs.flush(); }
    ReadInfoRecord rec3{}; bool got3=false;
    // Force no accepted candidates by using zero cbReadCountTotal
    std::vector<uint32_t> cbZero2(1,0);
    loader.loadMinimal(rf3, [&](const ReadInfoRecord &rec){ rec3=rec; got3=true; }, cbZero2);
    assert(got3 && rec3.cbIdx==(uint32_t)-1 && rec3.umi==(uint32_t)-1);
    assert(rec3.status==0);
    unlink(streamPath3.c_str());

    // Case D: oneExact without exact whitelist count => missing CB (status=0)
    P.pSolo.CBmatchWL.oneExact = true;
    std::vector<uint32_t> cbZero(1,0), nUzero(1,0), nMzero(1,0);
    SoloReadFeature rf4(SoloFeatureTypes::GeneFull, P, 3);
    std::string streamPath4 = P.outFileTmp + "/solo" + SoloFeatureTypes::Names[SoloFeatureTypes::GeneFull] + "_3";
    { std::ofstream ofs(streamPath4, std::ios::out | std::ios::trunc); ofs << "222 103 1 0\n"; ofs.flush(); }
    ReadInfoRecord rec4{}; bool got4=false;
    loader.loadMinimal(rf4, [&](const ReadInfoRecord &rec){ rec4=rec; got4=true; }, cbZero);
    assert(got4 && rec4.cbIdx==(uint32_t)-1 && rec4.umi==(uint32_t)-1);
    assert(rec4.status==0);
    unlink(streamPath4.c_str());

    // Case E: invalid UMI sentinel -> status=2
    P.pSolo.readIndexYes[SoloFeatureTypes::GeneFull] = false;
    P.pSolo.CBmatchWL.oneExact = false;
    SoloReadFeature rf5(SoloFeatureTypes::GeneFull, P, 4);
    std::string streamPath5 = P.outFileTmp + "/solo" + SoloFeatureTypes::Names[SoloFeatureTypes::GeneFull] + "_4";
    { std::ofstream ofs(streamPath5, std::ios::out | std::ios::trunc); ofs << "4294967295 104 0 0\n"; ofs.flush(); }
    ReadInfoRecord rec5{}; bool got5=false;
    std::vector<uint32_t> cbCounts5(1,1);
    loader.loadMinimal(rf5, [&](const ReadInfoRecord &rec){ rec5=rec; got5=true; }, cbCounts5);
    assert(got5 && rec5.status==2 && rec5.cbIdx==0 && rec5.umi==(uint32_t)-1);
    unlink(streamPath5.c_str());
    rmdir(P.outFileTmp.c_str());
    rmdir(baseDir.c_str());
    std::cout << "SoloReadInfoLoader_test: OK" << std::endl;
    return 0;
}


