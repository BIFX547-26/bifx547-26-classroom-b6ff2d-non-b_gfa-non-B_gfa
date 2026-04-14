/*
 * nulls.c
 * 
 * Initialize and clear global arrays
 * Part of the GFA (Genomic Feature Analyzer) package
 * 
 * For more information, see gfa.h and the main gfa.c file
 */

#include <stdlib.h>

void nulls(char line[], int n) {
    int i;
    for (i = 0; i < n + 1; i++)
        line[i] = '\00';
}
