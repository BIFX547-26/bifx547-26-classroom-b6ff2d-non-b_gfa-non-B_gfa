#include <Rcpp.h>
#include "gfa_api.h"

using namespace Rcpp;

// Helper to convert REP array to DataFrame
DataFrame rep_to_df(REP* reps, int count) {
    IntegerVector start(count);
    IntegerVector end(count);
    IntegerVector len(count);
    IntegerVector loop(count);
    IntegerVector num(count);
    IntegerVector sub(count);
    IntegerVector strand(count);
    IntegerVector special(count);

    if (reps == NULL) return DataFrame::create();

    for (int i = 0; i < count; i++) {
        start[i] = reps[i].start;
        end[i] = reps[i].end;
        len[i] = reps[i].len;
        loop[i] = reps[i].loop;
        num[i] = reps[i].num;
        sub[i] = reps[i].sub;
        strand[i] = reps[i].strand;
        special[i] = reps[i].special;
    }

    return DataFrame::create(
        Named("start") = start,
        Named("end") = end,
        Named("len") = len,
        Named("loop") = loop,
        Named("num") = num,
        Named("sub") = sub,
        Named("strand") = strand,
        Named("special") = special
    );
}

// [[Rcpp::export]]
List find_motifs_cpp(std::string sequence, List params_list) {
    GFA_Params params;
    init_default_params(&params);

    // Override defaults from R list
    if (params_list.containsElementNamed("minGQrep")) params.minGQrep = params_list["minGQrep"];
    if (params_list.containsElementNamed("maxGQspacer")) params.maxGQspacer = params_list["maxGQspacer"];
    if (params_list.containsElementNamed("minMRrep")) params.minMRrep = params_list["minMRrep"];
    if (params_list.containsElementNamed("maxMRspacer")) params.maxMRspacer = params_list["maxMRspacer"];
    if (params_list.containsElementNamed("minIRrep")) params.minIRrep = params_list["minIRrep"];
    if (params_list.containsElementNamed("maxIRspacer")) params.maxIRspacer = params_list["maxIRspacer"];
    if (params_list.containsElementNamed("shortIRcut")) params.shortIRcut = params_list["shortIRcut"];
    if (params_list.containsElementNamed("shortIRspacer")) params.shortIRspacer = params_list["shortIRspacer"];
    if (params_list.containsElementNamed("maxDRspacer")) params.maxDRspacer = params_list["maxDRspacer"];
    if (params_list.containsElementNamed("minDRrep")) params.minDRrep = params_list["minDRrep"];
    if (params_list.containsElementNamed("maxDRrep")) params.maxDRrep = params_list["maxDRrep"];
    if (params_list.containsElementNamed("minATracts")) params.minATracts = params_list["minATracts"];
    if (params_list.containsElementNamed("minATractSep")) params.minATractSep = params_list["minATractSep"];
    if (params_list.containsElementNamed("maxATractSep")) params.maxATractSep = params_list["maxATractSep"];
    if (params_list.containsElementNamed("maxAPRlen")) params.maxAPRlen = params_list["maxAPRlen"];
    if (params_list.containsElementNamed("minAPRlen")) params.minAPRlen = params_list["minAPRlen"];
    if (params_list.containsElementNamed("minZlen")) params.minZlen = params_list["minZlen"];
    if (params_list.containsElementNamed("minSTR")) params.minSTR = params_list["minSTR"];
    if (params_list.containsElementNamed("maxSTR")) params.maxSTR = params_list["maxSTR"];
    if (params_list.containsElementNamed("minSTRreps")) params.minSTRreps = params_list["minSTRreps"];
    if (params_list.containsElementNamed("minSTRbp")) params.minSTRbp = params_list["minSTRbp"];
    if (params_list.containsElementNamed("minCruciformRep")) params.minCruciformRep = params_list["minCruciformRep"];
    if (params_list.containsElementNamed("maxCruciformSpacer")) params.maxCruciformSpacer = params_list["maxCruciformSpacer"];
    if (params_list.containsElementNamed("minTriplexYRpercent")) params.minTriplexYRpercent = params_list["minTriplexYRpercent"];
    if (params_list.containsElementNamed("maxTriplexSpacer")) params.maxTriplexSpacer = params_list["maxTriplexSpacer"];
    if (params_list.containsElementNamed("maxSlippedSpacer")) params.maxSlippedSpacer = params_list["maxSlippedSpacer"];
    if (params_list.containsElementNamed("minKVscore")) params.minKVscore = params_list["minKVscore"];

    GFA_Results res = run_gfa_analysis(sequence.c_str(), &params);

    List out = List::create(
        Named("IR") = rep_to_df(irep, res.ireps),
        Named("MR") = rep_to_df(mrep, res.mreps),
        Named("DR") = rep_to_df(drep, res.dreps),
        Named("GQ") = rep_to_df(grep, res.greps),
        Named("Z") = rep_to_df(zrep, res.zreps),
        Named("STR") = rep_to_df(srep, res.sreps),
        Named("APR") = rep_to_df(arep, res.areps)
    );

    free_buffers();
    return out;
}
