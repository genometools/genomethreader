#ifndef GTHPRE_H
#define GTHPRE_H

#include "gth/gthoutput.h"
#include "gth/gthalphatype.h"
#include "gth/input.h"

int gthmakesureindexexists(const char *filename, const char *progname,
                           const char *proteinsmap, bool isreferencefile,
                           bool dnasuffix, bool oistab, bool completeindex,
                           GtUword *overall_maxlength,
                           GthAlphatype alphatype, GthOutput*, GtError*);

/* The following function checks if all inputfiles given by <genomicfiles> and
  <referencefiles> are valid (and their index exists).
  It also sets the output width used for genomic positions ``widthforgenpos''
  (part of <out>). */
int gthpreprocessinputfiles(GthInput *input,
                            bool gthconsensus,
                            bool noautoindex,
                            bool skipindexcheck,
                            bool maskpolyAtails,
                            bool online,
                            bool inverse,
                            const char *progname,
                            unsigned int translationtable,
                            GthOutput *out, GtError*);

#endif
