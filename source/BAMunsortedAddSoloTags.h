#ifndef CODE_BAMunsortedAddSoloTags
#define CODE_BAMunsortedAddSoloTags

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Genome.h"
#include "Solo.h"

void BAMunsortedAddSoloTags(const std::string &tmpPath,
                            const std::string &outBamPath,
                            Parameters &P,
                            Genome &genome,
                            Solo &solo);

#endif
