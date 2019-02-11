#ifndef GTHMKVTREE_H
#define GTHMKVTREE_H

#include "gth/gthoutput.h"
#include "types.h"

int gthcallmkvtree(const char *filename, const char *progname,
                   const char *proteinsmap, bool creatednaindex, bool transdnax,
                   bool addsuffix, bool oistab, bool completeindex, GthOutput*,
                   GtError*);

#endif
