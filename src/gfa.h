#ifndef GFA_H_
#define GFA_H_

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_REPS 2500000
#define MAX_DNA 300000000
#define MAXCOL 101
#define MAX_FASTA_SIZE 80
#define MAX_LINE 256

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
typedef int BOOLEAN;

typedef struct REP {
	int start;
	int loop;
	int len;
	int num;
	int end;
	int sub;
	int strand; //0 for plus, 1 for minus
	int special;//1 for yes, 0 for no
} REP;

typedef struct A_Tract {
	int strt;
	short len;
} A_Tract;

typedef struct potential_Bent_DNA {
	double a_center;
	int strt;
	int end;
} potential_Bent_DNA;

typedef struct G_Island {
	int strt;
	int len;
} G_Island;

extern int nGisls;
extern int nCisls;

extern char *dna;
extern char *dna2; //reverse complement DNA
extern char *dna3; //complement DNA

extern G_Island *gisle;
extern G_Island *rcgisle;
extern const G_Island null_gisle;

extern potential_Bent_DNA *pAPRs;

extern REP *irep;
extern REP *mrep;
extern REP *drep;
extern REP *grep;
extern REP *zrep;
extern REP *srep;
extern REP *arep;
extern const REP null_rep;

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#ifdef __cplusplus
}
#endif

#endif /* GFA_H_ */
