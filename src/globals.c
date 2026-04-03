#include "gfa.h"
#include <stddef.h>

potential_Bent_DNA *pAPRs = NULL;

G_Island *gisle = NULL;
G_Island *rcgisle = NULL;
const G_Island null_gisle;

int nGisls;
int nCisls;

char *dna = NULL;
char *dna2 = NULL; //reverse complement strand
char *dna3 = NULL; //complement strand

REP *mrep = NULL; //mirror
REP *irep = NULL; //inverted
REP *drep = NULL; //direct
REP *grep = NULL; //g-quadraplex
REP *zrep = NULL; //z-dna
REP *srep = NULL; //str
REP *arep = NULL; //a-phased-repeat
const REP null_rep;
