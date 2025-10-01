#include "BAMunsortedAddSoloTags.h"
#include "ErrorWarning.h"
#include "BAMfunctions.h"
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <thread>

// Debug flag - controlled by STAR_DEBUG_TAG environment variable
static const bool g_debugTag = (std::getenv("STAR_DEBUG_TAG") != nullptr);

// Debug logging helper with timestamp and thread ID
inline void logDebugTag(const std::string &msg) {
    if (!g_debugTag) return;
    
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    std::ostringstream thread_id;
    thread_id << std::this_thread::get_id();
    
    std::cerr << "[DEBUG_TAG " << std::put_time(&tm, "%H:%M:%S") 
              << " T" << thread_id.str().substr(thread_id.str().length()-4) << "] " 
              << msg << std::endl;
}

void BAMunsortedAddSoloTags(const std::string &tmpPath,
                            const std::string &outBamPath,
                            Parameters &P,
                            Genome &genome,
                            Solo &solo)
{
    // Guard: only proceed if CB/UB tags are requested
    if (!solo.pSolo.samAttrYes) {
        // Nothing to inject, just return
        return;
    }

    // Open the tmp file for reading (binary)
    std::ifstream tmpStream(tmpPath.c_str(), std::ios::binary);
    if (!tmpStream.is_open()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not open solo tmp file for reading: " << tmpPath << "\n";
        errOut << "SOLUTION: check that the file exists and has proper permissions\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    // Open the final unsorted BAM for writing
    BGZF *bgzfOut = bgzf_open(outBamPath.c_str(), ("w"+to_string((long long) P.outBAMcompression)).c_str());
    if (bgzfOut == NULL) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not create final unsorted BAM file: " << outBamPath << "\n";
        errOut << "SOLUTION: check that the directory exists and has proper permissions\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    // Write BAM header immediately
    outBAMwriteHeader(bgzfOut, P.samHeader, genome.chrNameAll, genome.chrLengthAll);

    // Allocate scratch buffers
    char *bam0 = nullptr;
    char *bam1 = nullptr;
    
    try {
        bam0 = new char[BAMoutput_oneAlignMaxBytes + sizeof(uint64)]; // space for record + trailer
        bam1 = new char[BAM_ATTR_MaxSize]; // workspace for addBAMtags
        
        if (g_debugTag) {
            logDebugTag("Allocated BAM buffers: bam0=" + std::to_string(BAMoutput_oneAlignMaxBytes + sizeof(uint64)) + 
                       " bytes, bam1=" + std::to_string(BAM_ATTR_MaxSize) + " bytes");
        }
    } catch (const std::bad_alloc &e) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: failed to allocate BAM buffers\n";
        errOut << "bam0 size: " << (BAMoutput_oneAlignMaxBytes + sizeof(uint64)) << " bytes\n";
        errOut << "bam1 size: " << BAM_ATTR_MaxSize << " bytes\n";
        errOut << "SOLUTION: reduce --runThreadN or increase system memory\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    // Loop until EOF: read records from tmp file and process them
    uint32 size_with_len;
    uint64 trailer;
    uint64 recordsProcessed = 0;
    
    if (g_debugTag) {
        logDebugTag("Starting BAM record processing loop");
    }
    
    while (tmpStream.read(reinterpret_cast<char*>(&size_with_len), sizeof(uint32))) {
        // Guardrail: Check record size before reading
        if (size_with_len > BAMoutput_oneAlignMaxBytes) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: BAM record size " << size_with_len << " exceeds maximum buffer size " << BAMoutput_oneAlignMaxBytes << "\n";
            errOut << "SOLUTION: increase BAMoutput_oneAlignMaxBytes or check for corrupted tmp file: " << tmpPath << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        
        // Read the payload (BAM record)
        if (!tmpStream.read(bam0, size_with_len)) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: could not read BAM record payload from solo tmp file: " << tmpPath << "\n";
            errOut << "SOLUTION: the tmp file may be corrupted\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }

        // Read the trailer (contains iReadAll)
        if (!tmpStream.read(reinterpret_cast<char*>(&trailer), sizeof(uint64))) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: could not read trailer from solo tmp file: " << tmpPath << "\n";
            errOut << "SOLUTION: the tmp file may be corrupted\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }

        // Phase 4: Validate the trailer contains a valid iread index
        // The trailer is a 64-bit packed structure: iReadAll in upper 32 bits, auxiliary data in lower 32 bits
        uint32_t iread = static_cast<uint32_t>(trailer >> 32); // Extract iReadAll from upper 32 bits
        uint32_t aux = static_cast<uint32_t>(trailer & 0xFFFFFFFFu); // Extract auxiliary data from lower 32 bits
        
        if (g_debugTag && recordsProcessed < 5) {
            logDebugTag("Record " + std::to_string(recordsProcessed) + 
                       " trailer: full=0x" + std::to_string(trailer) + 
                       ", iread=" + std::to_string(iread) + 
                       ", aux=0x" + std::to_string(aux));
        }
        
        if (iread >= solo.soloFeat[solo.pSolo.featureInd[solo.pSolo.samAttrFeature]]->readInfo.size()) {
            if (g_debugTag) {
                logDebugTag("ERROR: Invalid iread=" + std::to_string(iread) + 
                           " at record " + std::to_string(recordsProcessed) + 
                           " (readInfo.size=" + std::to_string(solo.soloFeat[solo.pSolo.featureInd[solo.pSolo.samAttrFeature]]->readInfo.size()) + 
                           ", full trailer=0x" + std::to_string(trailer) + ")");
            }
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: Invalid read index " << iread << " in trailer\n";
            errOut << "Record: " << recordsProcessed << ", readInfo.size: " << solo.soloFeat[solo.pSolo.featureInd[solo.pSolo.samAttrFeature]]->readInfo.size() << "\n";
            errOut << "Full trailer: 0x" << std::hex << trailer << std::dec << " (iread=" << iread << ", aux=0x" << std::hex << aux << std::dec << ")\n";
            errOut << "SOLUTION: check for corrupted tmp file or upstream indexing issues: " << tmpPath << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        
        // Phase 4: Debug logging for first few records
        if (g_debugTag && recordsProcessed < 5) {
            const auto& readData = solo.soloFeat[solo.pSolo.featureInd[solo.pSolo.samAttrFeature]]->readInfo[iread];
            bool hasCB = (readData.cb != (uint64)-1);
            bool hasUMI = (readData.umi != (uint32)-1);
            
            logDebugTag("Record " + std::to_string(recordsProcessed) + 
                       ": iread=" + std::to_string(iread) + 
                       ", hasCB=" + (hasCB ? "true" : "false") + 
                       ", hasUMI=" + (hasUMI ? "true" : "false"));
        }
        
        // Append the trailer to the BAM record (addBAMtags expects it there)
        memcpy(bam0 + size_with_len, &trailer, sizeof(uint64));
        
        uint32 size0 = size_with_len; // Initial size (will be updated by addBAMtags)
        
        // Call addBAMtags to inject CB/UB tags
        char *bam0_ptr = bam0; // addBAMtags takes char*& so we need a pointer to modify
        solo.soloFeat[solo.pSolo.featureInd[solo.pSolo.samAttrFeature]]->addBAMtags(bam0_ptr, size0, bam1);

        // Guardrail: Check final size after tag injection
        if (size0 > BAM_ATTR_MaxSize) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: BAM record size after tag injection " << size0 << " exceeds scratch buffer size " << BAM_ATTR_MaxSize << "\n";
            errOut << "SOLUTION: increase BAM_ATTR_MaxSize or check for excessive tag data\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }

        // Phase 4: Additional size validation after tag injection
        if (g_debugTag && size0 != size_with_len) {
            logDebugTag("Record " + std::to_string(recordsProcessed) + 
                       " size changed: " + std::to_string(size_with_len) + 
                       " -> " + std::to_string(size0) + " bytes (+tags)");
        }
        
        // Write the resulting record to the final BAM (without trailer)
        bgzf_write(bgzfOut, bam0_ptr, size0);
        recordsProcessed++;
    }

    // Check if we exited due to EOF or error
    if (!tmpStream.eof()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: error reading from solo tmp file: " << tmpPath << "\n";
        errOut << "SOLUTION: the tmp file may be corrupted\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    // Final debug logging
    if (g_debugTag) {
        logDebugTag("BAM processing completed: " + std::to_string(recordsProcessed) + " records processed");
    }
    
    // Cleanup: close files and free memory
    delete[] bam0;
    delete[] bam1;
    
    tmpStream.close();
    
    if (bgzf_close(bgzfOut) != 0) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not close final unsorted BAM file: " << outBamPath << "\n";
        errOut << "SOLUTION: check disk space and permissions\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    // Delete the tmp file on success
    if (remove(tmpPath.c_str()) != 0) {
        // Non-fatal warning - tmp file cleanup failed but processing succeeded
        P.inOut->logMain << "WARNING: could not delete temporary file: " << tmpPath << std::endl;
    }
}
