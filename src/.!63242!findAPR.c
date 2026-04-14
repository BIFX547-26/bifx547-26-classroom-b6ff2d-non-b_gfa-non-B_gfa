#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "gfa.h"

/*******************************************
 * A tract definition:
 *
 *Find all AT stretches of 4 to 9,
 * regardless of A/T order.
 *Find Longest A or AnTn run (NAs).
 *Start at first A.  Restart count at TA steps
