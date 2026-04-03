#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "gfa_api.h"

void init_default_params(GFA_Params *params) {
    params->minGQrep = 3;
    params->maxGQspacer = 7;
    params->minMRrep = 10;
    params->maxMRspacer = 100;
    params->minIRrep = 6;
    params->maxIRspacer = 100;
    params->shortIRcut = 9;
    params->shortIRspacer = 4;
    params->maxDRspacer = 10;
    params->minDRrep = 10;
    params->maxDRrep = 300;
    params->minATracts = 3;
    params->minATractSep = 10;
    params->maxATractSep = 11;
    params->maxAPRlen = 9;
    params->minAPRlen = 3;
    params->minZlen = 10;
    params->minSTR = 1;
    params->maxSTR = 9;
    params->minSTRreps = 3;
    params->minSTRbp = 10;
    params->minCruciformRep = 10;
    params->maxCruciformSpacer = 4;
    params->minTriplexYRpercent = 10;
    params->maxTriplexSpacer = 8;
    params->maxSlippedSpacer = 0;
    params->minKVscore = 33;
    params->findSTR = TRUE;
    params->findDR = TRUE;
    params->findMR = TRUE;
    params->findGQ = TRUE;
    params->findIR = TRUE;
    params->findZ = TRUE;
    params->findAPR = TRUE;
    params->doCruciform = TRUE;
    params->doTriplex = TRUE;
    params->doSlipped = TRUE;
    params->doKVzdna = TRUE;
}

// Internal function declarations
void rcdna(int ndna);
void cdna(int ndna);
int findIR(int minIRrep, int maxIRspacer, int shortIRcut, int shortIRspacer, int total_bases);
int findMR(int minMRrep, int maxMRspacer, int total_bases);
int findDR(int minDRrep, int maxDRrep, int maxDRspacer, int total_bases);
int findZDNA(int minZ, int total_bases);
int findSTR(int minSTR, int maxSTR, int minSTRbp, int minSTRreps, int total_bases);
int findAPR(int minAPRlen, int maxAPRlen, int minATracts, int total_bases);
int findGQ(int minGQrep, int maxGQspacer);
int getGislands(int minGQrep, int total_bases);
void is_subset(int nreps, char X, int max_loop, int limit);

int allocate_buffers(int dna_len) {
    free_buffers();

    dna = (char*)calloc(dna_len + 1, sizeof(char));
    dna2 = (char*)calloc(dna_len + 1, sizeof(char));
    dna3 = (char*)calloc(dna_len + 1, sizeof(char));
    
    irep = (REP*)calloc(MAX_REPS + 1, sizeof(REP));
    mrep = (REP*)calloc(MAX_REPS + 1, sizeof(REP));
    drep = (REP*)calloc(MAX_REPS + 1, sizeof(REP));
    grep = (REP*)calloc(MAX_REPS + 1, sizeof(REP));
    zrep = (REP*)calloc(MAX_REPS + 1, sizeof(REP));
    srep = (REP*)calloc(MAX_REPS + 1, sizeof(REP));
    arep = (REP*)calloc(MAX_REPS + 1, sizeof(REP));
    
    gisle = (G_Island*)calloc(5 * MAX_REPS + 1, sizeof(G_Island));
    rcgisle = (G_Island*)calloc(5 * MAX_REPS + 1, sizeof(G_Island));
    pAPRs = (potential_Bent_DNA*)calloc(5 * MAX_REPS + 1, sizeof(potential_Bent_DNA));

    if (!dna || !dna2 || !dna3 || !irep || !mrep || !drep || !grep || !zrep || !srep || !arep || !gisle || !rcgisle || !pAPRs) {
        free_buffers();
        return 0;
    }
    return 1;
}

void free_buffers() {
    if (dna) free(dna); dna = NULL;
    if (dna2) free(dna2); dna2 = NULL;
    if (dna3) free(dna3); dna3 = NULL;
    if (irep) free(irep); irep = NULL;
    if (mrep) free(mrep); mrep = NULL;
    if (drep) free(drep); drep = NULL;
    if (grep) free(grep); grep = NULL;
    if (zrep) free(zrep); zrep = NULL;
    if (srep) free(srep); srep = NULL;
    if (arep) free(arep); arep = NULL;
    if (gisle) free(gisle); gisle = NULL;
    if (rcgisle) free(rcgisle); rcgisle = NULL;
    if (pAPRs) free(pAPRs); pAPRs = NULL;
}

GFA_Results run_gfa_core(GFA_Params *params, int total_bases) {
    GFA_Results res = {0, 0, 0, 0, 0, 0, 0};
    if (total_bases <= 0) return res;

    if (params->findGQ || params->findIR) {
        cdna(total_bases);
    }
    if (params->findGQ || params->findAPR || params->findIR) {
        rcdna(total_bases);
    }

    if (params->findIR) {
        res.ireps = findIR(params->minIRrep, params->maxIRspacer, params->shortIRcut, params->shortIRspacer, total_bases);
        if (params->doCruciform) {
            is_subset(res.ireps, 'I', params->maxCruciformSpacer, params->minCruciformRep);
        }
    }

    if (params->findGQ) {
        getGislands(params->minGQrep, total_bases);
        res.greps = findGQ(params->minGQrep, params->maxGQspacer);
    }

    if (params->findMR) {
        res.mreps = findMR(params->minMRrep, params->maxMRspacer, total_bases);
        if (params->doTriplex) {
            is_subset(res.mreps, 'M', params->maxTriplexSpacer, params->minTriplexYRpercent);
        }
    }

    if (params->findDR) {
        res.dreps = findDR(params->minDRrep, params->maxDRrep, params->maxDRspacer, total_bases);
        if (params->doSlipped) {
            is_subset(res.dreps, 'D', params->maxSlippedSpacer, -999);
        }
    }

    if (params->findZ) {
        res.zreps = findZDNA(params->minZlen, total_bases);
        if (params->doKVzdna) {
            is_subset(res.zreps, 'Z', -999, params->minKVscore);
        }
    }

    if (params->findSTR) {
        res.sreps = findSTR(params->minSTR, params->maxSTR, params->minSTRbp, params->minSTRreps, total_bases);
    }

    if (params->findAPR) {
        res.areps = findAPR(params->minAPRlen, params->maxAPRlen, params->minATracts, total_bases);
    }

    return res;
}

GFA_Results run_gfa_analysis(const char* sequence, GFA_Params *params) {
    int i;
    int total_bases = 0;
    GFA_Results res = {0, 0, 0, 0, 0, 0, 0};

    int seq_len = strlen(sequence);
    if (!allocate_buffers(seq_len)) {
        return res;
    }

    i = 0;
    while (sequence[i] != '\0') {
        if (isalpha(sequence[i])) {
            dna[i] = (char)tolower(sequence[i]);
            i++;
        } else {
            // skip non-alpha
        }
    }
    total_bases = i;

    res = run_gfa_core(params, total_bases);
    
    return res;
}
