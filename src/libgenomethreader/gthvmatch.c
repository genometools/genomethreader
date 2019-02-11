#include "core/assert_api.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "gth/gthdef.h"
#include "gth/gthoutput.h"
#include "gth/chaining.h"
#include "libgenomethreader/gthvmatch.h"
#include "libgenomethreader/seq_con_multiseq.h"
#include "types.h"
#include "virtualdef.h"
#include "select.h"
#include "readvirt.pr"
#include "vmatch.pr"

#ifndef NOLICENSEMANAGER
#include "licensemanager.h"
#endif

#define VMATCHARGSIZE  512
#define QSPEEDUPARG    5

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
                                   bool dbmaskmatch)
{
  GthVmatchInfo *vmatch_info;
  GtUword i, offset = 0;
  Sint rval;
  char suffix[MAXSUFFIXLEN];

  gt_assert(!(noautoindex && usepolyasuffix));

  vmatch_info = gt_malloc(sizeof *vmatch_info);

  /* init argv */
  for (i = 0; i < NUMOFVMATCHARGS; i++) {
    vmatch_info->argv[i] = gt_malloc(sizeof (char) * VMATCHARGSIZE);
    vmatch_info->argv[i][0] = '\0';
  }

  /* save program name */
  gt_assert(offset < NUMOFVMATCHARGS);;
  rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "%s", progname);
  gt_assert(rval < VMATCHARGSIZE);;

  /* save direction */
  if (directmatches || !refseqisdna) {
  /* XXX:           ^^^^^^^^^^^^^^^^ remove this when vmatch can do -p */
    /* save flag for direct matches */
    gt_assert(offset < NUMOFVMATCHARGS);;
    rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "-d");
    gt_assert(rval < VMATCHARGSIZE);;
  }
  else {
    /* save flag for reverse complemented (palindromic) matches */
    gt_assert(offset < NUMOFVMATCHARGS);;
    rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "-p");
    gt_assert(rval < VMATCHARGSIZE);;
  }

  /* save verbose flag */
  gt_assert(offset < NUMOFVMATCHARGS);;
  rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "-v");
  gt_assert(rval < VMATCHARGSIZE);;

  /* save length flag... */
  gt_assert(offset < NUMOFVMATCHARGS);;
  rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "-l");
  gt_assert(rval < VMATCHARGSIZE);;

  /* ...and the actual length */
  gt_assert(offset < NUMOFVMATCHARGS);;
  if (refseqisdna) {
    rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, ""GT_WU"",
                    minmatchlength);
  }
  else {
    gt_assert(prminmatchlen != GT_UNDEF_ULONG);
    rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, ""GT_WU"",
                    prminmatchlen);
  }
  gt_assert(rval < VMATCHARGSIZE);;

  /* if -exact was not used save exdrop flag and seedlength */
  if (refseqisdna && !exact) { /* DNA matching */
    if (!hamming) {
      gt_assert(offset < NUMOFVMATCHARGS);;
      rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE,
                      "-seedlength");
      gt_assert(rval < VMATCHARGSIZE);
      gt_assert(offset < NUMOFVMATCHARGS);;
      rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, ""GT_WU"",
                      seedlength);
      gt_assert(rval < VMATCHARGSIZE);;

      if (edist) {
        gt_assert(offset < NUMOFVMATCHARGS);;
        rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "-e");
        gt_assert(rval < VMATCHARGSIZE);;
      }
      else {
        gt_assert(offset < NUMOFVMATCHARGS);;
        rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE,
                        "-exdrop");
        gt_assert(rval < VMATCHARGSIZE);;
      }

      gt_assert(offset < NUMOFVMATCHARGS);;
      rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, ""GT_WU"",
                      exdrop);
      gt_assert(rval < VMATCHARGSIZE);;
    }
    else
    {
      gt_assert(offset < NUMOFVMATCHARGS);;
      rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "-h");
      gt_assert(rval < VMATCHARGSIZE);;
      gt_assert(offset < NUMOFVMATCHARGS);;
      rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, ""GT_WU"",
                      hammingdistance);
      gt_assert(rval < VMATCHARGSIZE);;
    }
  }
  else if (!refseqisdna) { /* protein matching */
      gt_assert(prhdist != GT_UNDEF_ULONG && prseedlength != GT_UNDEF_ULONG);

      /* we have to state explicitly, that we want to match a protein index
         against the six frame translation of a DNA query sequence. */
      gt_assert(offset < NUMOFVMATCHARGS);;
      rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE,
                      "-dnavsprot");
      gt_assert(rval < VMATCHARGSIZE);;
      /* additionally we give the translation scheme number */
      gt_assert(offset < NUMOFVMATCHARGS);;
      rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, ""GT_WU"",
                      translationtable);
      gt_assert(rval < VMATCHARGSIZE);;

      /* the parameters used for the matching */
      gt_assert(offset < NUMOFVMATCHARGS);;
      rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "-h");
      gt_assert(rval < VMATCHARGSIZE);;
      gt_assert(offset < NUMOFVMATCHARGS);;
      rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE , ""GT_WU"",
                      prhdist);
      gt_assert(rval < VMATCHARGSIZE);;
      gt_assert(offset < NUMOFVMATCHARGS);;
      rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE,
                      "-seedlength");
      gt_assert(rval < VMATCHARGSIZE);
      gt_assert(offset < NUMOFVMATCHARGS);;
      rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, ""GT_WU"",
                      prseedlength);
      gt_assert(rval < VMATCHARGSIZE);;
  }

  /* set online flag if necessary */
  if (online) {
    gt_assert(offset < NUMOFVMATCHARGS);;
    rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "-online");
    gt_assert(rval < VMATCHARGSIZE);;
  }

  /* set dbmaskmatch if necessary */
  if (dbmaskmatch) {
    gt_assert(offset < NUMOFVMATCHARGS);;
    rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "-dbmaskmatch");
    gt_assert(rval < VMATCHARGSIZE);;
    gt_assert(offset < NUMOFVMATCHARGS);;
    rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "X");
    gt_assert(rval < VMATCHARGSIZE);;
  }

  /* save query file flag... */
  gt_assert(offset < NUMOFVMATCHARGS);;
  rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "-q");
  gt_assert(rval < VMATCHARGSIZE);;

  /* ...and the actual names */
  gt_assert(offset < NUMOFVMATCHARGS);;
  /* set suffix */
  if (refseqisdna) {
    if (usepolyasuffix)
      rval = snprintf(suffix, MAXSUFFIXLEN, ".%s", POLYASUFFIX);
    else if (noautoindex)
      rval = snprintf(suffix, MAXSUFFIXLEN, "%s", ""); /* XXX  */
    else
      rval = snprintf(suffix, MAXSUFFIXLEN, ".%s", DNASUFFIX);
    gt_assert(rval < MAXSUFFIXLEN);
  }
  else
    suffix[0] = '\0';
  rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE , "%s%s",
                  queryfilename , suffix);
  gt_assert(rval < VMATCHARGSIZE);;

  /* ...with substring specification, if one is defined */
  if (checksubstrspec && gth_input_use_substring_spec(input)) {
    gt_assert(offset < NUMOFVMATCHARGS);
    rval = snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE,
                    "("GT_WU","GT_WU")",
                    gth_input_genomic_substring_from(input),
                    gth_input_genomic_substring_to(input));
    gt_assert(rval < VMATCHARGSIZE);;
  }

  /* save index name */
  gt_assert(offset < NUMOFVMATCHARGS);;
  /* set suffix */
  if (refseqisdna) {
    if (usepolyasuffix)
      rval = snprintf(suffix, MAXSUFFIXLEN, ".%s", POLYASUFFIX);
    else if (noautoindex)
      rval = snprintf(suffix, MAXSUFFIXLEN, "%s", ""); /* change */
    else
      rval = snprintf(suffix, MAXSUFFIXLEN, ".%s", DNASUFFIX);
  }
  else
    rval = snprintf(suffix, MAXSUFFIXLEN, ".%s", proteinsmap);
  gt_assert(rval < MAXSUFFIXLEN);

  rval =  snprintf(vmatch_info->argv[offset++], VMATCHARGSIZE, "%s%s",
                   indexfilename, suffix);
  gt_assert(rval < VMATCHARGSIZE);;

  /* save argument counter (argc) */
  vmatch_info->argc = (Argctype) offset;

  /* init virtual trees */
  makeemptyvirtualtree(&vmatch_info->virtualtree);
  makeemptyvirtualtree(&vmatch_info->queryvirtualtree);
  makeemptyvirtualtree(&vmatch_info->sixframeofqueryvirtualtree);
  makeemptyvirtualtree(&vmatch_info->dnavirtualtree);

  return vmatch_info;
}

void gth_vmatch_info_delete(GthVmatchInfo *vmatch_info)
{
  GtUword i;

  /* free argv */
  for (i = 0; i < NUMOFVMATCHARGS; i++)
    gt_free(vmatch_info->argv[i]);

  /* free virtual trees */
  if (wrapvmatch(&vmatch_info->virtualtree,
                 &vmatch_info->queryvirtualtree,
                 &vmatch_info->sixframeofqueryvirtualtree,
                 &vmatch_info->dnavirtualtree)) {
    fprintf(stderr, "%s\n", messagespace());
    exit(EXIT_FAILURE);
  }

  gt_free(vmatch_info);
}

static Sint vmatch_match_processor(void *data, GT_UNUSED Multiseq *multiseq,
                                   GT_UNUSED Multiseq *querymultiseq,
                                   StoreMatch *storematch)
{
  GthMatchProcessorInfo *info = data;
  GthMatch match;
  /* XXX: remove this filter when -p in vmatch works */
  if (info->refseqisdna ||
      ( info->directmatches && !(storematch->Storeflag & FLAGPPRIGHTREVERSE)) ||
      (!info->directmatches &&  (storematch->Storeflag & FLAGPPRIGHTREVERSE))) {

    /* save match... */

    if (info->refseqisindex) {
      match.Storeseqnumgenomic     = storematch->Storeseqnum2;
      match.Storepositiongenomic   = storematch->Storeposition2;
      match.Storelengthgenomic     = storematch->Storelength2;
      match.Storeseqnumreference   = storematch->Storeseqnum1;
      match.Storepositionreference = storematch->Storeposition1;
      match.Storelengthreference   = storematch->Storelength1;
    }
    else {
      match.Storeseqnumgenomic     = storematch->Storeseqnum1;
      match.Storepositiongenomic   = storematch->Storeposition1;
      match.Storelengthgenomic     = storematch->Storelength1;
      match.Storeseqnumreference   = storematch->Storeseqnum2;
      match.Storepositionreference = storematch->Storeposition2;
      match.Storelengthreference   = storematch->Storelength2;
    }
    match.Storescore = EVALDISTANCE2SCORE(storematch->Storedistance,
                                          match.Storelengthreference,
                                          match.Storelengthgenomic);

    return gth_match_processor(data, info->gen_seq_con, info->ref_seq_con,
                               &match);
  }
  return 0;
}

void gth_vmatch_runner(GthVmatchInfo *vmatch_info,
                       GthShowVerbose showverbose,
                       GthShowVerboseVM showverboseVM,
                       GthMatchProcessorInfo *match_processor_info)
{
  gt_assert(vmatch_info && match_processor_info);

  if (match_processor_info->refseqisindex) {
    match_processor_info->gen_seq_con =
      gth_seq_con_multiseq_new_vstree(&vmatch_info->queryvirtualtree);
    match_processor_info->ref_seq_con =
      gth_seq_con_multiseq_new_vstree(&vmatch_info->virtualtree);
  }
  else {
    match_processor_info->gen_seq_con =
      gth_seq_con_multiseq_new_vstree(&vmatch_info->virtualtree);
    match_processor_info->ref_seq_con =
      gth_seq_con_multiseq_new_vstree(&vmatch_info->queryvirtualtree);
  }

  if (callvmatch(vmatch_info->argc,
                 (const char**) vmatch_info->argv,
                 match_processor_info,
                 NULL,
                 "vmatch_match_processor",
                 vmatch_match_processor,
                 showverboseVM,
                 showverbose == NULL ? NULL : stdout,
                 NULL,
                 &vmatch_info->virtualtree,
                 &vmatch_info->queryvirtualtree,
                 &vmatch_info->sixframeofqueryvirtualtree,
                 &vmatch_info->dnavirtualtree)) {
    fprintf(stderr,"%s\n", messagespace());
    exit(EXIT_FAILURE);
  }
}

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
                                       bool dbmaskmatch)
{
#ifndef NOLICENSEMANAGER
  lm_license_check_c();
#endif
  return gth_vmatch_info_new(checksubstrspec, input, queryfilename,
                             indexfilename, directmatches, refseqisdna,
                             progname, proteinsmap, exact, edist, hamming,
                             hammingdistance, minmatchlength, seedlength,
                             exdrop, prminmatchlen, prseedlength, prhdist,
                             translationtable, online, noautoindex,
                             usepolyasuffix, dbmaskmatch);
}

void gth_vmatch_matcher_arguments_delete(void *matcher_arguments)
{
#ifndef NOLICENSEMANAGER
  lm_license_check_d();
#endif
  gth_vmatch_info_delete(matcher_arguments);
}

void gth_vmatch_matcher_runner(void *matcher_arguments,
                               GthShowVerbose showverbose,
                               GthShowVerboseVM showverboseVM,
                               void *match_processor_info)
{
#ifndef NOLICENSEMANAGER
  lm_license_check_e();
#endif
  gth_vmatch_runner(matcher_arguments, showverbose, showverboseVM,
                    match_processor_info);
}
