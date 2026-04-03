#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stddef.h>
#include <ctype.h>

#include "gfa.h"
#include "gfa_api.h"

/********************************************************************
 *  Genomic Feature Analyzer; non-B motifs in nucleic acid sequence *
 ********************************************************************/

int main(int argc, char *argv[]) {

	//initialize file name and other strings
	char dna_filename[181];

	char gffout_filenameI[181];
	char gffout_filenameM[181];
	char gffout_filenameD[181];
	char gffout_filenameG[181];
	char gffout_filenameZ[181];
	char gffout_filenameS[181];
	char gffout_filenameA[181];

	char tsvout_filenameI[181];
	char tsvout_filenameM[181];
	char tsvout_filenameD[181];
	char tsvout_filenameG[181];
	char tsvout_filenameZ[181];
	char tsvout_filenameS[181];
	char tsvout_filenameA[181];

	char fasta_title[MAX_FASTA_SIZE + 1];
	char chrom[80];
	char seq_title[80];
	char wget_string[130]; //char array for storing wget command to tell php that prog. is done
	char chmod_stringTSV[130]; //set all to read, set group and owner to write
	char chmod_stringGFF[130]; //set group and owner to write
	strcpy(
			wget_string,
			"wget -S - \"http://nonb.abcc.ncifcrf.gov/modules/nBMSTc/controllers/notify.php?file=");
	strcpy(chmod_stringGFF, "chmod 664 ");
	strcpy(chmod_stringTSV, "chmod 664 ");

	//registered loop counters
	register int i = 0;
	register int j = 0;

	//total bases in input file
	int total_bases = 0;
	int fasta_count = 0;
	int fasta = 0; //main loop counter

    GFA_Params params;
    init_default_params(&params);

	//counters for tsv printing, only want to print header once
	int Ic, Mc, Dc, Zc, Ac, Gc, Sc;
	Ic = Mc = Dc = Zc = Ac = Gc = Sc = 0;

	BOOLEAN DO_WGET = TRUE; //system call to wget for PHP?
	BOOLEAN DO_CHMOD = FALSE; //system call to change output file permissions, everyone read, owner and group write
	BOOLEAN FATAL = FALSE; //fatal error?
	BOOLEAN ISEQ = FALSE; //input file given?
	BOOLEAN CHROM = FALSE; //chromosome name given as com line argument?
	//BOOLEAN KEEP_TIME = TRUE; //for benchmarking etc.

	FILE *dna_file; //input file name
	//output files
	FILE *gffout_fileI;
	FILE *gffout_fileM;
	FILE *gffout_fileD;
	FILE *gffout_fileG;
	FILE *gffout_fileZ;
	FILE *gffout_fileS;
	FILE *gffout_fileA;

	FILE *tsvout_fileI;
	FILE *tsvout_fileM;
	FILE *tsvout_fileD;
	FILE *tsvout_fileG;
	FILE *tsvout_fileZ;
	FILE *tsvout_fileS;
	FILE *tsvout_fileA;

	/*****************************
	 * Procedure Declarations  ***
	 *****************************
	 */

	void print_usage(char pgmname[]); //prints help text to screen
	int read_mult_fasta(FILE *dna_file, int fasta, char fasta_title[]);
	int get_fasta_count(FILE *dna_file);
	void print_gff_file(FILE *gffout_file, int nreps, char chrom[], char X,
			int total_bases);
	void print_tsv_file(FILE *tsvout_file, int nreps, char chrom[], char X,
			int nFasta, int total_bases);

	// if program with no arguments
	if (argc == 1) {
		print_usage(argv[0]);
		exit(1);
	}

	/*******************************
	 * get command line parameters *
	 *******************************
	 */
	for (i = 1; i < argc; i++) {

		//input file name: REQUIRED
		if (strncmp(argv[i], "-seq", 4) == 0) {
			memset((char *) dna_filename, '\0', 180);
			if (argv[i + 1] != NULL) {
				strncpy(dna_filename, argv[i + 1], strlen(argv[i + 1]));
				fprintf(stderr, "-seq value = %s\n", dna_filename);
				ISEQ = TRUE;
			} else {
				fprintf(
						stderr,
						"FATAL ERROR: No argument for -seq switch (input file name)\n");
				FATAL = TRUE;
			}
		}
		//Chromosome name, truncated to 10 chars, REQUIRED if .gff output = TRUE
		if (strncmp(argv[i], "-chrom", 6) == 0) {
			memset((char *) chrom, '\0', 10);
			if (argv[i + 1] != NULL) {
				j = strlen(argv[i + 1]);
				if (j > 10) {
					j = 10;
				}
				strncpy(chrom, argv[i + 1], j);
				CHROM = TRUE;
				fprintf(stderr, "-chrom value = %s\n", chrom);
			} else {
				fprintf(stderr, "FATAL ERROR: No argument for -chrom switch\n");
				FATAL = TRUE;
			}
		}

		//command line overrides for main motif default values
		if (strncmp(argv[i], "-minGQrep", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.minGQrep);
				fprintf(stderr, "-minGQrep value = %d\n", params.minGQrep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minGQrep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxGQspacer", 12) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.maxGQspacer);
				fprintf(stderr, "-maxGQspacer value = %d\n", params.maxGQspacer);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxGQspacer switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minMRrep", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.minMRrep);
				fprintf(stderr, "-minMRrep value = %d\n", params.minMRrep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minMRrep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minIRrep", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.minIRrep);
				fprintf(stderr, "-minIRrep value = %d\n", params.minIRrep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minIRrep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxMRspacer", 12) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.maxMRspacer);
				fprintf(stderr, "-maxMRspacer value = %d\n", params.maxMRspacer);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxMRspacer switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxIRspacer", 12) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.maxIRspacer);
				fprintf(stderr, "-maxIRspacer value = %d\n", params.maxIRspacer);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxIRspacer switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-shortIRcut", 11) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.shortIRcut);
				fprintf(stderr, "-shortIRcut value = %d\n", params.shortIRcut);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -shortIRcut switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-shortIRspacer", 14) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.shortIRspacer);
				fprintf(stderr, "-shortIRspacer value = %d\n", params.shortIRspacer);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -shortIRspacer switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxDRspacer", 12) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.maxDRspacer);
				fprintf(stderr, "-maxDRspacer value = %d\n", params.maxDRspacer);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxDRspacer switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minDRrep", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.minDRrep);
				fprintf(stderr, "-minDRrep value = %d\n", params.minDRrep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minDRrep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxDRrep", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.maxDRrep);
				fprintf(stderr, "-maxDRrep value = %d\n", params.maxDRrep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxDRrep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minATracts", 11) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.minATracts);
				fprintf(stderr, "-minATracts value = %d\n", params.minATracts);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minATtracts switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minATractSep", 13) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.minATractSep);
				fprintf(stderr, "-minATractSep value = %d\n", params.minATractSep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minATtractSep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxATractSep", 13) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.maxATractSep);
				fprintf(stderr, "-maxATractSep value = %d\n", params.maxATractSep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxATractSep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxAPRlen", 10) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.maxAPRlen);
				fprintf(stderr, "-maxAPRlen value = %d\n", params.maxAPRlen);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxAPRlen switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minAPRlen", 10) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.minAPRlen);
				fprintf(stderr, "-minAPRlen value = %d\n", params.minAPRlen);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minAPRlen switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minZlen", 8) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.minZlen);
				fprintf(stderr, "-minZlen value = %d\n", params.minZlen);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minZlen switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minSTR", 7) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.minSTR);
				fprintf(stderr, "-minSTR value = %d\n", params.minSTR);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minSTR switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxSTR", 7) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.maxSTR);
				fprintf(stderr, "-maxSTR value = %d\n", params.maxSTR);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxSTR switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minSTRbp", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.minSTRbp);
				fprintf(stderr, "-minSTRbp value = %d\n", params.minSTRbp);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minSTRbp switch\n");
				FATAL = TRUE;
			}
		}

		//command line overrides for subset motif default values
		if (strncmp(argv[i], "-minCruciformRep", 16) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.minCruciformRep);
				fprintf(stderr, "-minCruciformRep value = %d\n",
						params.minCruciformRep);
			} else {
				fprintf(
						stderr,
						"FATAL ERROR: No argument for -maxCruciformRep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxCruciformSpacer", 19) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &params.maxCruciformSpacer);
				fprintf(stderr, "-maxCruciformSpacer value = %d\n",
						params.maxCruciformSpacer);
			} else {
				fprintf(
						stderr,
						"FATAL ERROR: No argument for -maxCruciformSpacer switch\n");
				FATAL = TRUE;
			}
		}

		//command line boolean overrides
		if (strncmp(argv[i], "-skipZ", 6) == 0) {
			params.findZ = FALSE;
			fprintf(stderr, "-skipZ = Z-DNA motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipDR", 7) == 0) {
			params.findDR = FALSE;
			fprintf(stderr,
					"-skipDR = Direct Repeat motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipIR", 7) == 0) {
			params.findIR = FALSE;
			fprintf(stderr,
					"-skipIR = Inverted Repeat motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipMR", 7) == 0) {
			params.findMR = FALSE;
			fprintf(stderr,
					"-skipMR = Mirror Repeat motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipGQ", 7) == 0) {
			params.findGQ = FALSE;
			fprintf(stderr,
					"-skipGQ = G-Quadruplex motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipSTR", 8) == 0) {
			params.findSTR = FALSE;
			fprintf(
					stderr,
					"-skipSTR = Short Tandem Repeat motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipAPR", 8) == 0) {
			params.findAPR = FALSE;
			fprintf(stderr,
					"-skipAPR = A-Phased Repeat motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipSlipped", 12) == 0) {
			params.doSlipped = FALSE;
			fprintf(
					stderr,
					"-skipSlipped = Slipped subset of Direct Repeats will not be found\n");
		}
		if (strncmp(argv[i], "-skipTriplex", 12) == 0) {
			params.doTriplex = FALSE;
			fprintf(
					stderr,
					"-skipTriplex = Triples subset of Mirror Repeats will not be found\n");
		}
		if (strncmp(argv[i], "-skipCruciform", 14) == 0) {
			params.doCruciform = FALSE;
			fprintf(
					stderr,
					"-skipCruciform = Cruciform subset of Inverted Repeats will not be found\n");
		}
		if (strncmp(argv[i], "-skipKVzdna", 11) == 0) {
			params.doKVzdna = FALSE;
			fprintf(
					stderr,
					"-skipKVzdna = Karen Vasquex subset of Z-DNA will not be found\n");
		}
		if (strncmp(argv[i], "-skipWGET", 9) == 0) {
			DO_WGET = FALSE;
			fprintf(stderr, "-skipWGET = wget PHP call will be skipped\n");
		}
		if (strncmp(argv[i], "-doCHMOD", 10) == 0) {
			DO_CHMOD = TRUE;
			fprintf(stderr,
					"-doCHMOD = output files will get chmod 664 permissions\n");
		}
		//gff
		if (strncmp(argv[i], "-out", 4) == 0) {
			if (argv[i + 1] == NULL) {
				fprintf(
						stderr,
						"FATAL ERROR: No argument for -out switch (output file name prefix)\n");
				FATAL = TRUE;
			} else { //out file name given, create output files
				if (params.findIR) {
					memset((char *) gffout_filenameI, '\0', 80);
					memset((char *) tsvout_filenameI, '\0', 80);
					strncpy(gffout_filenameI, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameI, argv[i + 1], strlen(argv[i + 1]));
				}
				if (params.findMR) {
					memset((char *) gffout_filenameM, '\0', 80);
					memset((char *) tsvout_filenameM, '\0', 80);
					strncpy(gffout_filenameM, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameM, argv[i + 1], strlen(argv[i + 1]));
				}
				if (params.findDR) {
					memset((char *) gffout_filenameD, '\0', 80);
					memset((char *) tsvout_filenameD, '\0', 80);
					strncpy(gffout_filenameD, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameD, argv[i + 1], strlen(argv[i + 1]));
				}
				if (params.findGQ) {
					memset((char *) gffout_filenameG, '\0', 80);
					memset((char *) tsvout_filenameG, '\0', 80);
					strncpy(gffout_filenameG, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameG, argv[i + 1], strlen(argv[i + 1]));
				}
				if (params.findZ) {
					memset((char *) gffout_filenameZ, '\0', 80);
					memset((char *) tsvout_filenameZ, '\0', 80);
					strncpy(gffout_filenameZ, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameZ, argv[i + 1], strlen(argv[i + 1]));
				}
				if (params.findSTR) {
					memset((char *) gffout_filenameS, '\0', 80);
					memset((char *) tsvout_filenameS, '\0', 80);
					strncpy(gffout_filenameS, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameS, argv[i + 1], strlen(argv[i + 1]));
				}
				if (params.findAPR) {
					memset((char *) gffout_filenameA, '\0', 80);
					memset((char *) tsvout_filenameA, '\0', 80);
					strncpy(gffout_filenameA, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameA, argv[i + 1], strlen(argv[i + 1]));
				}
				strncat(wget_string, argv[i + 1], strlen(argv[i + 1]));
				strncat(chmod_stringTSV, argv[i + 1], strlen(argv[i + 1]));
				strncat(chmod_stringGFF, argv[i + 1], strlen(argv[i + 1]));
				strncat(wget_string, "\"", 2);
				strncat(chmod_stringTSV, "*.tsv", 5);
				strncat(chmod_stringGFF, "*.gff", 5);

			}

		}

	}

	/*************************************************
	 * check for required command line parameters ****
	 *************************************************
	 */
	if (!ISEQ) {
		fprintf(stderr, " FATAL ERROR: no sequence file given (-seq)\n");
		FATAL = TRUE;
	}
	if (!CHROM) {
		fprintf(
				stderr,
				"-chrom not given. Sequence will be named using first word of fasta title line.\n");
	}
	if (FATAL) {
		fprintf(stderr, " FATAL Errors, exiting program\n");
		exit(20);
	}

	//open fasta input file
	if ((dna_file = fopen(dna_filename, "r")) == NULL) {
		fprintf(stderr, "\nERROR - Cannot open dna input file %s \n",
				dna_filename);
		exit(2);
	}

	if (params.findIR) {
		strncat(gffout_filenameI, "_IR.gff", 7);
		if ((gffout_fileI = fopen(gffout_filenameI, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameI);
			exit(3);
		}
	}
	if (params.findMR) {
		strncat(gffout_filenameM, "_MR.gff", 7);
		if ((gffout_fileM = fopen(gffout_filenameM, "w")) == NULL) {

			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameM);
			exit(3);
		}
	}
	if (params.findDR) {
		strncat(gffout_filenameD, "_DR.gff", 7);
		if ((gffout_fileD = fopen(gffout_filenameD, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameD);
			exit(3);
		}
	}
	if (params.findGQ) {
		strncat(gffout_filenameG, "_GQ.gff", 7);
		if ((gffout_fileG = fopen(gffout_filenameG, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameG);
			exit(3);
		}
	}
	if (params.findZ) {
		strncat(gffout_filenameZ, "_Z.gff", 6);
		if ((gffout_fileZ = fopen(gffout_filenameZ, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameZ);
			exit(3);
		}
	}
	if (params.findSTR) {
		strncat(gffout_filenameS, "_STR.gff", 8);
		if ((gffout_fileS = fopen(gffout_filenameS, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameS);
			exit(3);
		}
	}
	if (params.findAPR) {
		strncat(gffout_filenameA, "_APR.gff", 8);
		if ((gffout_fileA = fopen(gffout_filenameA, "w")) == NULL) {

			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameA);
			exit(3);
		}
	}

	if (params.findIR) {
		strncat(tsvout_filenameI, "_IR.tsv", 7);
		if ((tsvout_fileI = fopen(tsvout_filenameI, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameI);
			exit(3);
		}
	}
	if (params.findMR) {
		strncat(tsvout_filenameM, "_MR.tsv", 7);
		if ((tsvout_fileM = fopen(tsvout_filenameM, "w")) == NULL) {

			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameM);
			exit(3);
		}
	}
	if (params.findDR) {
		strncat(tsvout_filenameD, "_DR.tsv", 7);
		if ((tsvout_fileD = fopen(tsvout_filenameD, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameD);
			exit(3);
		}
	}
	if (params.findGQ) {
		strncat(tsvout_filenameG, "_GQ.tsv", 7);
		if ((tsvout_fileG = fopen(tsvout_filenameG, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameG);
			exit(3);
		}
	}
	if (params.findZ) {
		strncat(tsvout_filenameZ, "_Z.tsv", 6);
		if ((tsvout_fileZ = fopen(tsvout_filenameZ, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameZ);
			exit(3);
		}
	}
	if (params.findSTR) {
		strncat(tsvout_filenameS, "_STR.tsv", 8);
		if ((tsvout_fileS = fopen(tsvout_filenameS, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameS);
			exit(3);
		}
	}
	if (params.findAPR) {
		strncat(tsvout_filenameA, "_APR.tsv", 8);
		if ((tsvout_fileA = fopen(tsvout_filenameA, "w")) == NULL) {

			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameA);
			exit(3);
		}
	}

	fasta_count = get_fasta_count(dna_file);
	fprintf(stderr, "Fasta Sections = %d\n", fasta_count);

	/************************************
	 * Master Loop for each FASTA entry *
	 ************************************
	 */
	for (fasta = 1; fasta <= fasta_count; fasta++) {

		fclose(dna_file);
		dna_file = fopen(dna_filename, "r");

        // Simple way to get sequence length: scan until next '>' or EOF
        // For CLI, we can just use MAX_DNA if we are lazy, but let's be better.
        // Actually, for simplicity in the fix, I'll use a large buffer.
        // But the real fix is to make sure ALL array accesses are checked.
        if (!allocate_buffers(MAX_DNA)) {
            fprintf(stderr, "Failed to allocate memory for sequence!\n");
            exit(21);
        }

		total_bases = read_mult_fasta(dna_file, fasta, fasta_title);
		if (CHROM) {
			strcpy(seq_title, chrom);
		} else {
			int w = 1;
			while (fasta_title[w] != '\0' && !isspace(fasta_title[w])) {
				w++;
			}
			memset((char *) seq_title, '\0', 80);
			strncpy(seq_title, &fasta_title[0], w);
		}

		if (CHROM && fasta_count > 1) {
			char tmp_str[12];
			sprintf(tmp_str, "%d", fasta);
			strncat(seq_title, "_", 1);
			strncat(seq_title, tmp_str, 5);
		}

		fprintf(stderr, "Total Bases Read=%d \n", total_bases);
		if (total_bases <= 0) {
			fprintf(stderr, " Empty Sequence!. Skipping!\n");
            free_buffers();
			continue;
		}

        GFA_Results results = run_gfa_core(&params, total_bases);

		if (fasta_count == 1) {
            if (params.findIR && results.ireps > 0) {
                print_gff_file(gffout_fileI, results.ireps, seq_title, 'I', total_bases);
                fclose(gffout_fileI);
                print_tsv_file(tsvout_fileI, results.ireps, seq_title, 'I', ++Ic, total_bases);
                fclose(tsvout_fileI);
            }
            if (params.findGQ && results.greps > 0) {
                print_gff_file(gffout_fileG, results.greps, seq_title, 'G', total_bases);
                fclose(gffout_fileG);
                print_tsv_file(tsvout_fileG, results.greps, seq_title, 'G', ++Gc, total_bases);
                fclose(tsvout_fileG);
            }
            if (params.findMR && results.mreps > 0) {
                print_gff_file(gffout_fileM, results.mreps, seq_title, 'M', total_bases);
                fclose(gffout_fileM);
                print_tsv_file(tsvout_fileM, results.mreps, seq_title, 'M', ++Mc, total_bases);
                fclose(tsvout_fileM);
            }
            if (params.findDR && results.dreps > 0) {
                print_gff_file(gffout_fileD, results.dreps, seq_title, 'D', total_bases);
                fclose(gffout_fileD);
                print_tsv_file(tsvout_fileD, results.dreps, seq_title, 'D', ++Dc, total_bases);
                fclose(tsvout_fileD);
            }
            if (params.findZ && results.zreps > 0) {
                print_gff_file(gffout_fileZ, results.zreps, seq_title, 'Z', total_bases);
                fclose(gffout_fileZ);
                print_tsv_file(tsvout_fileZ, results.zreps, seq_title, 'Z', ++Zc, total_bases);
                fclose(tsvout_fileZ);
            }
            if (params.findSTR && results.sreps > 0) {
                print_gff_file(gffout_fileS, results.sreps, seq_title, 'S', total_bases);
                fclose(gffout_fileS);
                print_tsv_file(tsvout_fileS, results.sreps, seq_title, 'S', ++Sc, total_bases);
                fclose(tsvout_fileS);
            }
            if (params.findAPR && results.areps > 0) {
                print_gff_file(gffout_fileA, results.areps, seq_title, 'A', total_bases);
                fclose(gffout_fileA);
                print_tsv_file(tsvout_fileA, results.areps, seq_title, 'A', ++Ac, total_bases);
                fclose(tsvout_fileA);
            }
		} else {
			if (params.findIR && results.ireps > 0) print_gff_file(gffout_fileI, results.ireps, seq_title, 'I', total_bases);
			if (params.findMR && results.mreps > 0) print_gff_file(gffout_fileM, results.mreps, seq_title, 'M', total_bases);
			if (params.findDR && results.dreps > 0) print_gff_file(gffout_fileD, results.dreps, seq_title, 'D', total_bases);
			if (params.findZ && results.zreps > 0) print_gff_file(gffout_fileZ, results.zreps, seq_title, 'Z', total_bases);
			if (params.findGQ && results.greps > 0) print_gff_file(gffout_fileG, results.greps, seq_title, 'G', total_bases);
			if (params.findSTR && results.sreps > 0) print_gff_file(gffout_fileS, results.sreps, seq_title, 'S', total_bases);
			if (params.findAPR && results.areps > 0) print_gff_file(gffout_fileA, results.areps, seq_title, 'A', total_bases);

			if (params.findIR && results.ireps > 0) print_tsv_file(tsvout_fileI, results.ireps, seq_title, 'I', ++Ic, total_bases);
			if (params.findMR && results.mreps > 0) print_tsv_file(tsvout_fileM, results.mreps, seq_title, 'M', ++Mc, total_bases);
			if (params.findDR && results.dreps > 0) print_tsv_file(tsvout_fileD, results.dreps, seq_title, 'D', ++Dc, total_bases);
			if (params.findZ && results.zreps > 0) print_tsv_file(tsvout_fileZ, results.zreps, seq_title, 'Z', ++Zc, total_bases);
			if (params.findGQ && results.greps > 0) print_tsv_file(tsvout_fileG, results.greps, seq_title, 'G', ++Gc, total_bases);
			if (params.findSTR && results.sreps > 0) print_tsv_file(tsvout_fileS, results.sreps, seq_title, 'S', ++Sc, total_bases);
			if (params.findAPR && results.areps > 0) print_tsv_file(tsvout_fileA, results.areps, seq_title, 'A', ++Ac, total_bases);
		}
        free_buffers();
	}

    fclose(dna_file);
	if (fasta_count > 1) {
        if (params.findIR) { fclose(gffout_fileI); fclose(tsvout_fileI); }
        if (params.findMR) { fclose(gffout_fileM); fclose(tsvout_fileM); }
        if (params.findDR) { fclose(gffout_fileD); fclose(tsvout_fileD); }
        if (params.findZ) { fclose(gffout_fileZ); fclose(tsvout_fileZ); }
        if (params.findGQ) { fclose(gffout_fileG); fclose(tsvout_fileG); }
        if (params.findSTR) { fclose(gffout_fileS); fclose(tsvout_fileS); }
        if (params.findAPR) { fclose(gffout_fileA); fclose(tsvout_fileA); }
	}
	if (DO_WGET) system(wget_string);
	if (DO_CHMOD) {
		system(chmod_stringTSV);
		system(chmod_stringGFF);
	}

	exit(0);
}
