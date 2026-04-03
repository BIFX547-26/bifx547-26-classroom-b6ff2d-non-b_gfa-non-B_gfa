#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "gfa.h"

int getAtracts(int minAT, int maxAT, int total_bases) {
	int n = 0;
	int n_rc = 0;
	int nPATs = 0;
	int maxATlen = 0;
	int maxTlen = 0;
	int Alen = 0;
	int Tlen = 0;
	int ATlen = 0;
	int TAlen = 0;
	int maxATend = 0;
	int ATend = 0;
	int maxATlen_rc = 0;
	int maxTlen_rc = 0;
	int Alen_rc = 0;
	int Tlen_rc = 0;
	int ATlen_rc = 0;
	int TAlen_rc = 0;
	int maxATend_rc = 0;
	register int i;
	int nAs = 0;
	i = 0;
	int strt = 0;

	while (i < total_bases) {
		if (i < total_bases && (dna[i] == 'a' || dna[i] == 't')) {
			nAs++;
		} else {
			if ((nAs >= minAT) && (nAs <= maxAT)) {
				strt = i - nAs + 1;
				ATend = (strt + nAs);
				Alen = 0; Tlen = 0; ATlen = 0; maxATlen = 0; maxTlen = 0; TAlen = 0; maxATend = 0;
				Alen_rc = 0; Tlen_rc = 0; ATlen_rc = 0; maxATlen_rc = 0; maxTlen_rc = 0; TAlen_rc = 0; maxATend_rc = 0;
				n_rc = (total_bases - ATend);
				for (n = strt - 1; n < ATend - 1; n++) {
					if (n < 0 || n >= total_bases) continue;
					n_rc++;
					if (dna[n] == 'a') {
						Tlen = 0; TAlen = 0;
						if (n > 0 && dna[n - 1] == 't') { Alen = 0; ATlen = 0; }
						else { Alen++; ATlen++; }
					}
					if (n_rc >= 0 && n_rc < total_bases && dna2[n_rc] == 'a') {
						Tlen_rc = 0; TAlen_rc = 0;
						if (n_rc > 0 && dna2[n_rc - 1] == 't') { Alen_rc = 0; ATlen_rc = 0; }
						else { Alen_rc++; ATlen_rc++; }
					}
					if (dna[n] == 't') {
						if (TAlen < Alen) { TAlen++; ATlen++; }
						else { Tlen++; TAlen = 0; ATlen = 0; Alen = 0; }
					}
					if (n_rc >= 0 && n_rc < total_bases && dna2[n_rc] == 't') {
						if (TAlen_rc < Alen_rc) { TAlen_rc++; ATlen_rc++; }
						else { Tlen_rc++; TAlen_rc = 0; ATlen_rc = 0; Alen_rc = 0; }
					}
					if (maxATlen < ATlen) { maxATlen = ATlen; maxATend = n; }
					if (maxTlen < Tlen) maxTlen = Tlen;
					if (maxATlen_rc < ATlen_rc) { maxATlen_rc = ATlen_rc; maxATend_rc = n_rc; }
					if (maxTlen_rc < Tlen_rc) maxTlen_rc = Tlen_rc;
				}
				if (((maxATlen - maxTlen) >= minAT) || ((maxATlen_rc - maxTlen_rc) >= minAT)) {
					if (nPATs < 5 * MAX_REPS) {
						pAPRs[nPATs].end = strt + nAs;
						pAPRs[nPATs].strt = strt;
						if ((maxATlen - maxTlen) >= (maxATlen_rc - maxTlen_rc)) {
							pAPRs[nPATs].a_center = ((double) maxATend - (((double) maxATlen - 1) / 2)) + 1;
						} else {
							pAPRs[nPATs].a_center = total_bases - (((double) maxATend_rc - (((double) maxATlen_rc - 1) / 2)));
						}
						nPATs++;
					}
				}
			}
			nAs = 0;
		}
		i++;
	}
	fprintf(stderr, "n potential a tracts = %d\n", nPATs);
	return (nPATs);
}

int findAPR(int minAPR, int maxAPR, int minATracts, int total_bases) {
	register int i;
	int nProcessedATs;
	nProcessedATs = getAtracts(minAPR, maxAPR, total_bases);
	int tracts = 1;
	double distToNext = 0;
	int ndx = 0;
	for (i = 0; i < nProcessedATs - 1; i++) {
		distToNext = pAPRs[i + 1].a_center - pAPRs[i].a_center;
		if ((distToNext <= 11.1) && (distToNext >= 9.9)) {
			tracts++;
		} else {
			if (tracts >= minATracts) {
				if (ndx < MAX_REPS) {
					arep[ndx].start = pAPRs[(i - tracts) + 1].strt;
					arep[ndx].num = tracts;
					arep[ndx].strand = 0;
					arep[ndx].len = tracts;
					arep[ndx].end = pAPRs[i].end - 1;
					ndx++;
				}
			}
			tracts = 1;
		}
	}
	// check last one
	if (tracts >= minATracts && nProcessedATs > 0) {
        if (ndx < MAX_REPS) {
            arep[ndx].start = pAPRs[(nProcessedATs - tracts)].strt;
            arep[ndx].num = tracts;
            arep[ndx].strand = 0;
            arep[ndx].len = tracts;
            arep[ndx].end = pAPRs[nProcessedATs - 1].end - 1;
            ndx++;
        }
	}
	return (ndx);
}
