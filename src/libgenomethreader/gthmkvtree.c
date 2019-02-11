#include <stdbool.h>
#include "core/assert_api.h"
#include "core/compat.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/ma.h"
#include "gth/gthdef.h"
#include "gth/gthoutput.h"
#include "libgenomethreader/gthmkvtree.h"
#include "types.h"
#include "virtualdef.h"
#include "readvirt.pr"

/*
  This file contains all functions necessary to interface with mkvtree.
*/

#define NUMOFMKVTREEARGS        17

#define MKVTREEARGSIZE          512

#define CHECKOFFSET\
        gt_assert(offset < NUMOFMKVTREEARGS);

#define CHECKRVAL\
        gt_assert(rval < MKVTREEARGSIZE);

#define POLYA_SMAPFILE                  "TransDNAX"
#define POLYA_PL_VALUE                  7

#define SAVE_SMAP_OPTION\
        CHECKOFFSET;\
        rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE,\
                         "-smap");\
        CHECKRVAL;\
        CHECKOFFSET;\
        rval =  snprintf(gthmkvtreeinfo->argv[offset++] , MKVTREEARGSIZE, "%s",\
                         smapfile);\
        CHECKRVAL

/*
  This structure stores the argument vector and the argument counter for
  the mkvtree call.
*/
typedef struct
{
  Argctype argc;                /* number of arguments for mkvtree */
  char *argv[NUMOFMKVTREEARGS]; /* argument vector for mkvtree */
} Gthmkvtreeinfo;

static int setsmapfile(char *smapfile, const char *mapping, GtError *err)
{
  int rval, had_err = 0;
  gt_error_check(err);
  if (!getenv(GTHDATAENVNAME)) {
    gt_error_set(err, "$%s not defined. Please set correctly", GTHDATAENVNAME);
    had_err = -1;
  }
  if (!had_err) {
    rval = snprintf(smapfile, PATH_MAX+1, "%s%c%s", getenv(GTHDATAENVNAME),
                    GT_PATH_SEPARATOR, mapping);
    gt_assert(rval < PATH_MAX + 1);
    if (!gt_file_exists(smapfile)) {
      gt_error_set(err, "cannot open smap file '%s'", smapfile);
      had_err = -1;
    }
  }
  if (had_err) {
    /* check for file relative to binary */
    GtStr *path;
    int new_err;
    path = gt_str_new();
    new_err = gt_file_find_exec_in_path(path, gt_error_get_progname(err), NULL);
    if (!new_err) {
      gt_str_append_char(path, GT_PATH_SEPARATOR);
      gt_str_append_cstr(path, GTHDATADIRNAME);
      gt_str_append_char(path, GT_PATH_SEPARATOR);
      gt_str_append_cstr(path, mapping);
      if (gt_file_exists(gt_str_get(path))) {
        rval = snprintf(smapfile, PATH_MAX+1, "%s", gt_str_get(path));
        gt_assert(rval < PATH_MAX + 1);
        gt_error_unset(err);
        had_err = 0;
      }
    }
    gt_str_delete(path);
  }
  return had_err;
}

static int initGthmkvtreeinfo(Gthmkvtreeinfo *gthmkvtreeinfo,
                              const char *filename, bool creatednaindex,
                              bool transdnax, bool addsuffix, Uint demand,
                              const char *progname, const char *proteinsmap,
                              GtError *err)
{
  char smapfile[PATH_MAX+1];
  Uint i, offset = 0;
  Sint rval;

  gt_error_check(err);
  gt_assert(!(!creatednaindex && transdnax));

  /* init argv */
  for (i = 0; i < NUMOFMKVTREEARGS; i++) {
    gthmkvtreeinfo->argv[i] = gt_malloc(sizeof (char) * MKVTREEARGSIZE);
    gthmkvtreeinfo->argv[i][0] = '\0';
  }

  /* save program name */
  CHECKOFFSET;
  rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "%s",
                   progname);
  CHECKRVAL;

  /* save alphabet type */
  if (creatednaindex) {
    /* set DNA alphabet */
    if (transdnax) {
      if (setsmapfile(smapfile, POLYA_SMAPFILE, err))
        return -1;
      SAVE_SMAP_OPTION;
    }
    else {
      CHECKOFFSET;
      rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "-dna");
      CHECKRVAL;
    }
  }
  else {
    /* set protein alphabet */
    if (strcmp(proteinsmap, "protein") == 0) {
      CHECKOFFSET;
      rval = snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE,
                      "-protein");
      CHECKRVAL;
    }
    else {
      if (setsmapfile(smapfile, proteinsmap, err))
        return -1;
      SAVE_SMAP_OPTION;
    }
  }

  /* save verbose flag */
  CHECKOFFSET;
  rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "-v");
  CHECKRVAL;

  /* set flags for different tables */
  if (demand & OISTAB) {
    CHECKOFFSET;
    rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "-ois");
    CHECKRVAL;
  }
  if (demand & TISTAB) {
    CHECKOFFSET;
    rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "-tis");
    CHECKRVAL;
  }

  if (demand & BCKTAB) {
    CHECKOFFSET;
    rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "-bck");
    CHECKRVAL;
    CHECKOFFSET;
    rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "-pl");
    CHECKRVAL;
  }

  if (demand & SUFTAB) {
    CHECKOFFSET;
    rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "-suf");
    CHECKRVAL;
  }
  if (demand & BWTTAB) {
    CHECKOFFSET;
    rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "-bwt");
    CHECKRVAL;
  }
  if (demand & LCPTAB) {
    CHECKOFFSET;
    rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "-lcp");
    CHECKRVAL;
  }

  if (demand & STI1TAB) {
    CHECKOFFSET;
    rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "-sti1");
    CHECKRVAL;
  }

  /* save filename for which index needs to be computed */
  CHECKOFFSET;
  rval = snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "-db");
  CHECKRVAL;
  CHECKOFFSET;
  rval = snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "%s",
                  filename);
  CHECKRVAL;

  /* save name of index */
  CHECKOFFSET;
  rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE ,
                   "-indexname");
  CHECKRVAL;

  if (addsuffix) {
    if (transdnax) {
      CHECKOFFSET;
      rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "%s.%s",
                       filename, POLYASUFFIX);
      CHECKRVAL;
    }
    else if (!creatednaindex) {
      CHECKOFFSET;
      rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "%s.%s",
                       filename, proteinsmap);
      CHECKRVAL;
    }
    else {
      CHECKOFFSET;
      rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "%s.%s",
                       filename, DNASUFFIX);
      CHECKRVAL;
    }
  }
  else {
    CHECKOFFSET;
    rval =  snprintf(gthmkvtreeinfo->argv[offset++], MKVTREEARGSIZE, "%s",
                     filename);
    CHECKRVAL;
  }

  /* save argument counter (argc) */
  gthmkvtreeinfo->argc = (Argctype) offset;

  return 0;
}

static void freeGthmkvtreeinfo(Gthmkvtreeinfo *gthmkvtreeinfo)
{
  unsigned int i;
  if (!gthmkvtreeinfo) return;
  /* free argv */
  for (i = 0; i < NUMOFMKVTREEARGS; i++)
    gt_free(gthmkvtreeinfo->argv[i]);
}

/* The following function is used to interface with mkvtree.
  <filename> is the name of the file for which an index needs to be
  constructred. */
int gthcallmkvtree(const char *filename, const char *progname,
                   const char *proteinsmap, bool creatednaindex, bool transdnax,
                   bool addsuffix, bool oistab, bool completeindex,
                   GthOutput *out, GtError *err)
{
  Gthmkvtreeinfo gthmkvtreeinfo;
  Virtualtree virtualtree;
  Uint demand = TISTAB;

  gt_error_check(err);

  if (oistab)
    demand |= OISTAB;
  if (completeindex)
    demand |= SUFTAB | BWTTAB | LCPTAB | BCKTAB | STI1TAB;

  /* init argc and argv for mkvtree */
  if (initGthmkvtreeinfo(&gthmkvtreeinfo, filename, creatednaindex, transdnax,
                         addsuffix, demand, progname, proteinsmap, err)) {
    return -1;
  }

  if (out && out->showverbose)
    out->showverbose("call mkvtree to compute index");

  /* call mkvtree */
  makeemptyvirtualtree(&virtualtree);
  if (callmkvtree(gthmkvtreeinfo.argc, (const char**) gthmkvtreeinfo.argv, true,
                  &virtualtree, true, out ? out->showverboseVM : NULL)) {
    fprintf(stderr,"%s\n", messagespace());
    exit(EXIT_FAILURE);
  }
  if (freevirtualtree(&virtualtree)) {
    fprintf(stderr,"%s\n", messagespace());
    exit(EXIT_FAILURE);
  }

  /* free space */
  freeGthmkvtreeinfo(&gthmkvtreeinfo);

  return 0;
}
