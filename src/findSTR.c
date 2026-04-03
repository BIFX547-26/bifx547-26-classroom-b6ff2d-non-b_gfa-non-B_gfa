#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include "gfa.h"

static int compar4(const void *a, const void *b) {
	int i = 0;
	i = ((REP *) a)->start - ((REP *) b)->start;
	if (i == 0) i = ((REP *) a)->len - ((REP *) b)->len;
	return (i);
}

static void removeSTR(int nreps, int toRemove) {
	int i = 0;
	for (i = toRemove; i < nreps; i++) {
		srep[i] = srep[i + 1];
	}
}

int nonBstr(int start, int len) {
	int code = 0;
	int j = 0;
	int i = 0;
	BOOLEAN isEven = (len % 2 == 0);
	BOOLEAN isSymetric = TRUE;
	BOOLEAN isPUPY = TRUE;
	BOOLEAN isComp = TRUE;
	
	j = start + len - 2;
	if (len >= 2) {
		for (i = 0; i <= (len / 2) - 1; i++) {
			if (dna[start + i - 1] != dna[j]) {
				isSymetric = FALSE;
			}
			if (isEven) {
				char c1 = dna[start + i - 1];
				char c2 = dna[j];
				if ((c1 == 'a' && c2 != 't') || (c1 == 't' && c2 != 'a') ||
                    (c1 == 'c' && c2 != 'g') || (c1 == 'g' && c2 != 'c')) {
					isComp = FALSE;
				}
			} else {
				isComp = FALSE;
			}
			j--;
		}
		for (i = start; i < (len + start - 1); i++) {
			char c1 = dna[i];
			char c2 = dna[i - 1];
			if (((c1 == 'a' || c1 == 'g') && (c2 == 'a' || c2 == 'g')) ||
                ((c1 == 't' || c1 == 'c') && (c2 == 't' || c2 == 'c'))) {
				isPUPY = FALSE;
			}
		}
	} else {
		isComp = FALSE;
		isSymetric = TRUE;
		isPUPY = FALSE;
	}

	if (isEven) code += 1;
	if (isPUPY) code += 2;
	if (isSymetric) code += 4;
	if (isComp) code += 8;
	return (code);
}

int filterSTRs(int nSTRs) {
	qsort(srep, nSTRs, sizeof(*srep), compar4);
	int i = 0;
	for (i = 1; i < nSTRs; i++) {
		if (srep[i].end <= srep[i - 1].end) {
			removeSTR(nSTRs, i);
			nSTRs--;
			i--;
		}
	}
	return (nSTRs);
}

int findSTR(int minSTR, int maxSTR, int minSTRlen, int minReps, int total_bases) {
	register int i, j;
	int ndx = 0;
	int rpsz = 0;
	int reps = 1;
	int remainder = 0;
	int rs = 0;
	int re = 0;
	for (i = 0; i < (total_bases - minSTRlen); i++) {
		while (i < total_bases && dna[i] == 'n' && i != (total_bases - 1)) {
			i++;
		}
		for (rpsz = minSTR; rpsz <= maxSTR; rpsz++) {
			reps = 1;
			j = i + rpsz;
			while (j + rpsz <= total_bases && strncmp(&dna[i], &dna[j], rpsz) == 0) {
				reps++;
				j = j + rpsz;
			}
			if (reps >= minReps) {
				remainder = 0;
				rs = i;
				re = j;
				while (re < total_bases && dna[rs] == dna[re]) {
					remainder++;
					rs++;
					re++;
				}
				if ((((reps * rpsz) + remainder) >= minSTRlen)) {
					if (ndx >= MAX_REPS) return ndx;
					if (ndx == 0 || srep[ndx - 1].end < re) {
						srep[ndx].start = i + 1;
						srep[ndx].end = re;
						srep[ndx].num = reps;
						srep[ndx].loop = nonBstr(i + 1, rpsz);
						srep[ndx].len = rpsz;
						srep[ndx].sub = remainder;
						srep[ndx].strand = 0;
						ndx++;
						i = re - minSTRlen + 1;
						break;
					}
				}
			}
		}
	}
	return (ndx);
}
