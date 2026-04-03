#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "gfa.h"

/*******************************************
 * Start looking for G_Islands *************
 *******************************************
 */
void getGislands(int minGQ, int total_bases) {
	register int i;
	int ngs;
	int ncs;

	nGisls = 0;
	nCisls = 0;
	i = 0;
	ngs = 0;
	ncs = 0;
	while (i <= total_bases) {
		if (i < total_bases && dna[i] == 'g') {
			ngs++;
		}
		else {
			if (ngs >= minGQ) {
				if (nGisls < 5 * MAX_REPS) {
					gisle[nGisls].strt = i - ngs + 1;
					gisle[nGisls].len = ngs;
					nGisls++;
				}
			}
			ngs = 0;
		}
		if (i < total_bases && dna[i] == 'c') {
			ncs++;
		}
		else {
			if (ncs >= minGQ) {
				if (nCisls < 5 * MAX_REPS) {
					rcgisle[nCisls].strt = i - ncs + 1;
					rcgisle[nCisls].len = ncs;
					nCisls++;
				}
			}
			ncs = 0;
		}
		i++;
	}
}

/*******************************************
 * Process G_Islands to look for ***********
 * G-quadruplex forming motifs *************
 *******************************************
 */
int findGQ(int minGQ, int maxGQspacer) {
	int ndx = 0;
	int nIls;
	G_Island *islands;
	register int i, i2;
	int npos, j, k, m;
	int nposMax, maxGQ;
	int conIls;
	int strand;
	for (strand = 0; strand < 2; strand++) {
		if (strand == 0) {
			nIls = nGisls;
			islands = gisle;
		} else {
			nIls = nCisls;
			islands = rcgisle;
		}
		for (i = 0; i < nIls; i++) {
			conIls = 1;
			npos = (int) (floor((islands[i].len + 1) / (minGQ + 1)));
			i2 = i + 1;
			while ((i2 < nIls) && ((islands[i2].strt - (islands[i2 - 1].strt
					+ islands[i2 - 1].len)) <= maxGQspacer)) {
				conIls++;
				npos += (int) (floor((islands[i2].len + 1) / (minGQ + 1)));
				i2++;
			}
			if (npos >= 4) {
				maxGQ = minGQ;
				for (j = i; j<i2; j++) {
					for (k = islands[j].len; k>maxGQ; k--) {
						nposMax = (int) (floor((islands[j].len + 1) / (k + 1)));
						for (m = j+1; m<i2; m++) {
							nposMax+= (int) (floor((islands[m].len + 1) / (k + 1)));
							if (nposMax>=4) {
								maxGQ = k;
								break;
							}
							if ((int) (floor((islands[m].len + 1) / (k + 1)))==0) {
								if((m+1 < i2) && (islands[m+1].strt>(islands[m-1].strt + islands[m-1].len + maxGQspacer))) {
									break;
								}
							}
						}
					}
				}
				
				if (ndx >= MAX_REPS) return ndx;

				grep[ndx].start = islands[i].strt;
				grep[ndx].num = npos;
				grep[ndx].sub = conIls;
				grep[ndx].len = maxGQ;
				grep[ndx].end = (islands[i2 - 1].strt + islands[i2 - 1].len) - 1;
				grep[ndx].strand = strand;
				ndx++;
			}
			i = i + conIls - 1;
		}
	}
	return (ndx);
}
