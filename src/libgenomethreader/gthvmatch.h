#ifndef GTHVMATCH_H
#define GTHVMATCH_H

#include "gth/gthmatch.h"
#include "gth/input.h"
#include "gth/chaining.h"
#include "types.h"
#include "match.h"
#include "virtualdef.h"
#include "maxfiles.h"

#define NUMOFVMATCHARGS         17

typedef struct
{
  Argctype argc;               /* number of arguments for vmatch */
  char *argv[NUMOFVMATCHARGS]; /* argument vector for vmatch */
  Virtualtree virtualtree,     /* the virtual trees needed by vmatch */
              queryvirtualtree,
              sixframeofqueryvirtualtree,
              dnavirtualtree;
} GthVmatchInfo;

GthVmatchInfo* gth_vmatch_info_new(bool checksubstrspec,
                                   GthInput *input,
                                   const char *queryfilename,
                                   const char *indexfilename,
                                   bool directmatches,
                                   bool refseqisdna,
                                   const char *progname,
                                   char *proteinsmap,
                                   bool exact,
                                   bool edist,
                                   bool hamming,
                                   GtUword hammingdistance,
                                   GtUword minmatchlength,
                                   GtUword seedlength,
                                   GtUword exdrop,
                                   GtUword prminmatchlen,
                                   GtUword prseedlength,
                                   GtUword prhdist,
                                   GtUword translationtable,
                                   bool online,
                                   bool noautoindex,
                                   bool usepolyasuffix,
                                   bool dbmaskmatch);
void gth_vmatch_info_delete(GthVmatchInfo*);
void gth_vmatch_runner(GthVmatchInfo*, GthShowVerbose, GthShowVerboseVM,
                       GthMatchProcessorInfo *match_processor_info);

/* similar stuff with prototypes suitable as GthMatcher */
void* gth_vmatch_matcher_arguments_new(bool checksubstrspec,
                                       GthInput *input,
                                       const char *queryfilename,
                                       const char *indexfilename,
                                       bool directmatches,
                                       bool refseqisdna,
                                       const char *progname,
                                       char *proteinsmap,
                                       bool exact,
                                       bool edist,
                                       bool hamming,
                                       GtUword hammingdistance,
                                       GtUword minmatchlength,
                                       GtUword seedlength,
                                       GtUword exdrop,
                                       GtUword prminmatchlen,
                                       GtUword prseedlength,
                                       GtUword prhdist,
                                       GtUword translationtable,
                                       bool online,
                                       bool noautoindex,
                                       bool usepolyasuffix,
                                       bool dbmaskmatch);
void gth_vmatch_matcher_arguments_delete(void *matcher_arguments);
void gth_vmatch_matcher_runner(void *matcher_arguments, GthShowVerbose,
                               GthShowVerboseVM, void *match_processor_info);

#endif
