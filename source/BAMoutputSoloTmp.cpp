#include "BAMoutputSoloTmp.h"
#include "BAMTagBuffer.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"
#include <cstring>

BAMoutputSoloTmp::BAMoutputSoloTmp(ofstream* tmpStreamIn, Parameters &Pin) 
    : tmpStream(tmpStreamIn), P(Pin), tagBuffer(nullptr) {
    
    bufferSize = defaultBufferSize;
    buffer = new char[bufferSize];
    bufferUsed = 0;
    
    if (!tmpStream || !tmpStream->is_open()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: BAMoutputSoloTmp could not initialize with closed stream\n";
        errOut << "DEBUG: tmpStream=" << (tmpStream ? "not null" : "null") << "\n";
        if (tmpStream) {
            errOut << "DEBUG: tmpStream->is_open()=" << (tmpStream->is_open() ? "true" : "false") << "\n";
            errOut << "DEBUG: tmpStream->good()=" << (tmpStream->good() ? "true" : "false") << "\n";
        }
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
}

BAMoutputSoloTmp::BAMoutputSoloTmp(ofstream* tmpStreamIn, Parameters &Pin, BAMTagBuffer* tagBufferIn) 
    : tmpStream(tmpStreamIn), P(Pin), tagBuffer(tagBufferIn) {
    
    bufferSize = defaultBufferSize;
    buffer = new char[bufferSize];
    bufferUsed = 0;
    
    if (!tmpStream || !tmpStream->is_open()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: BAMoutputSoloTmp tmpStream is not open\n";
        if (tmpStream) {
            errOut << "DEBUG: tmpStream->is_open()=" << (tmpStream->is_open() ? "true" : "false") << "\n";
            errOut << "DEBUG: tmpStream->good()=" << (tmpStream->good() ? "true" : "false") << "\n";
        }
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
}

BAMoutputSoloTmp::~BAMoutputSoloTmp() {
    if (buffer) {
        delete[] buffer;
        buffer = nullptr;
    }
}

void BAMoutputSoloTmp::unsortedOneAlign(char *bamIn, uint bamSize, uint iReadAll) {
    if (bamSize == 0) return; // no output
    
    // Guardrail: Check bamSize is reasonable
    if (bamSize > BAMoutput_oneAlignMaxBytes) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: BAM record size " << bamSize << " exceeds maximum " << BAMoutput_oneAlignMaxBytes << " in tmp writer\n";
        errOut << "SOLUTION: increase BAMoutput_oneAlignMaxBytes or check upstream BAM generation\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    
    // Calculate total space needed: size field + BAM record + trailer
    uint32 totalSize = sizeof(uint32) + bamSize + sizeof(uint64);
    
    // Check if we need to flush buffer
    if (bufferUsed + totalSize > bufferSize) {
        flush();
    }
    
    // If single record is larger than buffer, write directly
    if (totalSize > bufferSize) {
        // Write directly to stream with mutex protection
        if (g_threadChunks.threadBool) pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
        
        // Write size field
        uint32 sizeWithLen = bamSize;
        tmpStream->write(reinterpret_cast<char*>(&sizeWithLen), sizeof(uint32));
        
        // Write BAM record
        tmpStream->write(bamIn, bamSize);
        
        // Write trailer: pack iReadAll into upper 32 bits
        uint64 trailer = ((uint64)iReadAll) << 32;
        tmpStream->write(reinterpret_cast<char*>(&trailer), sizeof(uint64));
        
        if (g_threadChunks.threadBool) pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);
        return;
    }
    
    // Buffer the record
    char* writePos = buffer + bufferUsed;
    
    // Write size field to buffer
    uint32 sizeWithLen = bamSize;
    memcpy(writePos, &sizeWithLen, sizeof(uint32));
    writePos += sizeof(uint32);
    
    // Write BAM record to buffer
    memcpy(writePos, bamIn, bamSize);
    writePos += bamSize;
    
    // Write trailer to buffer: pack iReadAll into upper 32 bits
    uint64 trailer = ((uint64)iReadAll) << 32;
    memcpy(writePos, &trailer, sizeof(uint64));
    
    bufferUsed += totalSize;
}

void BAMoutputSoloTmp::unsortedOneAlign(char *bamIn, uint bamSize, uint iReadAll, const BAMRecordMeta& meta) {
    // Store metadata in tag buffer if available
    if (tagBuffer) {
        tagBuffer->append(meta);
    }
    
    // Delegate to the original method for BAM output
    unsortedOneAlign(bamIn, bamSize, iReadAll);
}

void BAMoutputSoloTmp::flush() {
    if (bufferUsed == 0) return;
    
    // Write buffer to stream with mutex protection
    if (g_threadChunks.threadBool) pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
    tmpStream->write(buffer, bufferUsed);
    if (g_threadChunks.threadBool) pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);
    
    bufferUsed = 0; // reset buffer
}

void BAMoutputSoloTmp::close() {
    flush(); // flush any remaining data
    
    if (tmpStream && tmpStream->is_open()) {
        tmpStream->close();
    }
}
