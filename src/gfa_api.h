#ifndef GFA_API_H_
#define GFA_API_H_

#include "gfa.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GFA_Params {
    int minGQrep;
    int maxGQspacer;
    int minMRrep;
    int maxMRspacer;
    int minIRrep;
    int maxIRspacer;
    int shortIRcut;
    int shortIRspacer;
    int maxDRspacer;
    int minDRrep;
    int maxDRrep;
    int minATracts;
    int minATractSep;
    int maxATractSep;
    int maxAPRlen;
    int minAPRlen;
    int minZlen;
    int minSTR;
    int maxSTR;
    int minSTRreps;
    int minSTRbp;
    int minCruciformRep;
    int maxCruciformSpacer;
    int minTriplexYRpercent;
    int maxTriplexSpacer;
    int maxSlippedSpacer;
    int minKVscore;
    BOOLEAN findSTR;
    BOOLEAN findDR;
    BOOLEAN findMR;
    BOOLEAN findGQ;
    BOOLEAN findIR;
    BOOLEAN findZ;
    BOOLEAN findAPR;
    BOOLEAN doCruciform;
    BOOLEAN doTriplex;
    BOOLEAN doSlipped;
    BOOLEAN doKVzdna;
} GFA_Params;

typedef struct GFA_Results {
    int ireps;
    int mreps;
    int dreps;
    int greps;
    int zreps;
    int sreps;
    int areps;
} GFA_Results;

void init_default_params(GFA_Params *params);

// Explicit memory management
int allocate_buffers(int dna_len);
void free_buffers();

// Process a single sequence string
GFA_Results run_gfa_analysis(const char* sequence, GFA_Params *params);

// Core analysis logic assuming 'dna' global is already populated
GFA_Results run_gfa_core(GFA_Params *params, int total_bases);

#ifdef __cplusplus
}
#endif

#endif
