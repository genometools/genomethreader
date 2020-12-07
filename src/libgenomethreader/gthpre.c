#include <math.h>
#include <stdbool.h>
#include "core/assert_api.h"
#include "core/compat_api.h"
#include "core/fa_api.h"
#include "core/fileutils_api.h"
#include "core/str_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "gth/gthdef.h"
#include "gth/gthoutput.h"
#include "gth/input.h"
#include "libgenomethreader/findfile.h"
#include "libgenomethreader/gthmkvtree.h"
#include "libgenomethreader/gthpolyafunc.h"
#include "libgenomethreader/gthpre.h"
#include "libgenomethreader/gthvmatch.h"
#include "types.h"
#include "select.h"
#include "virtualdef.h"
#include "fhandledef.h"
#include "filehandle.pr"
#include "multiseq-adv.pr"
#include "readvirt.pr"
#include "scanpaths.pr"
#include "vmatch.pr"

/* copied from readvirt.c */
#define MAXSUFFIXSIZE 4

#define INFOSUFFIX      ".info"
#define INFOSUFFIXLEN   5

/* this definitions are saved in the polya .info file */
#define POLYA_INFOLINELEN         80
#define POLYA_VERSIONLINE         "polyaversion=0.9"
#define POLYA_MINMATCHLENLINE     "polyaminmatchlen="
#define POLYA_HAMMINGDISTANCELINE "polyahammingdistance="

#define POLYA_FASTAFILENAME       "polyatail"

#define POLYA_MINMATCHLEN         14
#define POLYA_HAMMINGDISTANCE     1

/* This function checks if <alpha> is a strict DNA alphabet.
   That is, that both the lowercase and the uppercase version of ACGTU are
   mapped to the corresponding integer. */
static bool isstrictdnaalphabet(Alphabet *alpha)
{
  if (alpha->symbolmap['a'] == 0 && alpha->symbolmap[ 'A'] == 0 &&
      alpha->symbolmap['c'] == 1 && alpha->symbolmap[ 'C'] == 1 &&
      alpha->symbolmap['g'] == 2 && alpha->symbolmap[ 'G'] == 2 &&
      alpha->symbolmap['t'] == 3 && alpha->symbolmap[ 'T'] == 3 &&
      alpha->symbolmap['u'] == 3 && alpha->symbolmap[ 'U'] == 3) {
    return true;
  }
  return false;
}

/*
  This function checks if the tables given by <demand> are mappable and have the
  correct length. That is, that their length equals <totallength>.
  This function is similar to the function mapvirtualtreeifyoucan from
  readvirt.c, but it directly frees a table after it has been mapped.
*/
static bool tablesaremappableandhavecorrectlength(char *indexname,
                                                  Uint demand, Uint numofcodes,
                                                  Uint totallength)
{
  char tmpfilename[PATH_MAX+MAXSUFFIXSIZE+1];
  void *table;
  size_t numofbytes;

  if (demand & BCKTAB) {
    sprintf(tmpfilename, "%s.bck", indexname);
    table = gt_fa_mmap_read(tmpfilename, &numofbytes, NULL);
    if (!table)
      return false;
    if (numofbytes != 2 * numofcodes *  sizeof (Uint)) {
      gt_fa_xmunmap(table);
      return false;
    }
    gt_fa_xmunmap(table);
  }
  if (demand & SUFTAB) {
    sprintf(tmpfilename, "%s.suf", indexname);
    table = gt_fa_mmap_read(tmpfilename, &numofbytes, NULL);
    if (!table)
      return false;
    if (numofbytes != (totallength + 1) *  sizeof (Uint)) {
      gt_fa_xmunmap(table);
      return false;
    }
    gt_fa_xmunmap(table);
  }
  if (demand & BWTTAB) {
    sprintf(tmpfilename, "%s.bwt", indexname);
    table = gt_fa_mmap_read(tmpfilename, &numofbytes, NULL);
    if (!table)
      return false;
    if (numofbytes != (totallength + 1) * sizeof (unsigned char)) {
      gt_fa_xmunmap(table);
      return false;
    }
    gt_fa_xmunmap(table);
  }
  if (demand & LCPTAB) {
    sprintf(tmpfilename, "%s.lcp", indexname);
    table = gt_fa_mmap_read(tmpfilename, &numofbytes, NULL);
    if (!table)
      return false;
    if (numofbytes != (totallength + 1) *  sizeof (unsigned char)) {
      gt_fa_xmunmap(table);
      return false;
    }
    gt_fa_xmunmap(table);
  }
  if (demand & STI1TAB) {
    sprintf(tmpfilename, "%s.sti1", indexname);
    table = gt_fa_mmap_read(tmpfilename, &numofbytes, NULL);
    if (!table)
      return false;
    if (numofbytes != (totallength + 1) *  sizeof (unsigned char)) {
      gt_fa_xmunmap(table);
      return false;
    }
    gt_fa_xmunmap(table);
  }

  /* demand contains only BCKTAB | SUFTAB | BWTTAB | LCPTAB | STI1TAB */
  gt_assert((demand & ~(BCKTAB | SUFTAB | BWTTAB | LCPTAB | STI1TAB)) == 0);

  return true;
}

static bool polyafileisdifferent(FILE *fp)
{
  GtStr *line;
  char polyainfoline[POLYA_INFOLINELEN];

  /* init */
  line = gt_str_new();

  /* check if version number is the same */
  if (gt_str_read_next_line(line, fp) ==  EOF) { /* read first line */
    gt_str_delete(line);
    return true;
  }
  if (strcmp(gt_str_get(line), POLYA_VERSIONLINE)) {
    gt_str_delete(line);
    return true;
  }

  /* check if minmatchlen is the same */
  gt_str_reset(line);
  if (gt_str_read_next_line(line, fp) ==  EOF) { /* second */
    gt_str_delete(line);
    return true;
  }
  sprintf(polyainfoline, "%s%u", POLYA_MINMATCHLENLINE, POLYA_MINMATCHLEN);
  if (strcmp(gt_str_get(line), polyainfoline))
  {
    gt_str_delete(line);
    return true;
  }

  /* check if hamming distance is the same */
  gt_str_reset(line);
  if (gt_str_read_next_line(line, fp) ==  EOF) { /* third */
    gt_str_delete(line);
    return true;
  }
  sprintf(polyainfoline, "%s%u", POLYA_HAMMINGDISTANCELINE,
          POLYA_HAMMINGDISTANCE);
  if (strcmp(gt_str_get(line), polyainfoline)) {
    gt_str_delete(line);
    return true;
  }

  /* free */
  gt_str_delete(line);

  return false;
}

static void writeinfofile(char *filename)
{
  char polyainfoline[POLYA_INFOLINELEN];
  FILE *fp;

  /* open file */
  fp = gt_fa_xfopen(filename, "w");
  gt_assert(fp);

  /* write version line */
  gt_xfwrite(POLYA_VERSIONLINE,  sizeof (char), strlen(POLYA_VERSIONLINE), fp);
  gt_xfwrite("\n", sizeof (char), 1, fp);

  /* write minmatchlen line */
  sprintf(polyainfoline, "%s%u", POLYA_MINMATCHLENLINE, POLYA_MINMATCHLEN);
  gt_xfwrite(polyainfoline,  sizeof (char),  strlen(polyainfoline), fp);
  gt_xfwrite("\n", sizeof (char), 1, fp);

  /* write hammingdistance line */
  sprintf(polyainfoline, "%s%u", POLYA_HAMMINGDISTANCELINE,
          POLYA_HAMMINGDISTANCE);
  gt_xfwrite(polyainfoline,  sizeof (char),  strlen(polyainfoline), fp);
  gt_xfwrite("\n", sizeof (char), 1, fp);

  /* close .info file */
  gt_fa_xfclose(fp);
}

static int maskpolyAtailsandcreateindex(const char *filename,
                                        const char *progname,
                                        GtStr *proteinsmap, bool online,
                                        bool inverse, Uint translationtable,
                                        GthOutput *out, GtError *err)
{
  bool createmaskedfile = false,
       createindex = false,
       tablesok;
  FILE *fp;
  char maskedfilename[PATH_MAX+MAXSUFFIXLEN+1],
       dnafilename[PATH_MAX+MAXSUFFIXLEN+1],
       maskedinfofilename[PATH_MAX+MAXSUFFIXLEN+INFOSUFFIXLEN+1];
  GthVmatchInfo *vmatch_info;
  Virtualtree virtualtree;
  Uint numofcodes, totallength;
  SelectBundle selectbundle;
  int had_err = 0;

  gt_error_check(err);

  /* set the functions in the selectbundle */
  selectbundle.selectmatchHeader     = NULL;
  selectbundle.selectmatchInit       = NULL;
  selectbundle.selectmatch           = gthpolyaselectmatch;
  selectbundle.selectmatchWrap       = NULL;
  selectbundle.selectmatchFinaltable = NULL;

  /* check if masked file exists */
  sprintf(maskedfilename, "%s.%s", filename, POLYASUFFIX);
  sprintf(maskedinfofilename, "%s.%s%s", filename, POLYASUFFIX, INFOSUFFIX);
  if (out->showverbose) {
    out->showverbose("check if the following file containing masked reference "
                     "sequences exists:");
    out->showverbose(maskedfilename);
  }
  fp = gt_fa_fopen(maskedfilename, "r", NULL);
  if (fp) {
    gt_fa_xfclose(fp);
    if (gt_file_is_newer(filename, maskedfilename)) {
      if (out->showverbose) {
        out->showverbose("file exists, but is too old => create masked file "
                         "again");
      }
      createmaskedfile = true;
    }
    else {
      /* masked file exists and is up-to-date */
      if (out->showverbose) {
        out->showverbose("file exists, check if the corresponding .info file "
                         "exists");
      }

      /* check if .info file exists */
      fp = gt_fa_fopen(maskedinfofilename, "r", NULL);
      if (!fp) {
        /* .info file not there */
        if (out->showverbose)
          out->showverbose(".info file not there => create masked file again");
        createmaskedfile = true;
      }
      else {
        createmaskedfile = polyafileisdifferent(fp);
        gt_fa_xfclose(fp);
        if (out->showverbose) {
          if (createmaskedfile) {
            out->showverbose(".info file is not as it supposed to be => "
                             "create masked file again");
          }
          else {
            out->showverbose(".info file exists and was made with the same "
                             "parameters");
          }
        }
      }
    }
  }
  else {
    /* masked file does not exist */
    if (out->showverbose)
      out->showverbose("masked file does not exist => create it");
    createmaskedfile = true;
  }

  /* create masked file if necessary */
  if (createmaskedfile) {
    int rval = 0;
    GtStr *path = gt_str_new();

    sprintf(dnafilename, "%s.%s", POLYA_FASTAFILENAME, DNASUFFIX);
    had_err = gth_find_file(dnafilename, GTHDATAENVNAME, GTHDATADIRNAME, path,
                            err);
    if (!had_err) {
       rval = snprintf(dnafilename, PATH_MAX+1, "%s", gt_str_get(path));
       gt_assert(rval < PATH_MAX + 1);
       /* clip off DNASUFFIX */
       gt_assert(strlen(dnafilename) > 5);
       dnafilename[strlen(dnafilename)-4] = '\0';
    }
    gt_str_delete(path);

    /* call vmatch and save masked file */
    if (!had_err) {
      /* open file pointer for masked reference file */
      fp = gt_fa_xfopen(maskedfilename, "w");
      gt_assert(fp);

      /* prepare argv for vmatch */
      vmatch_info =
        gth_vmatch_info_new(false,
                            NULL,
                            dnafilename,
                            filename,
                            true,
                            true,
                            progname,
                            gt_str_get(proteinsmap),
                            false,
                            false,
                            true,
                            POLYA_HAMMINGDISTANCE,
                            POLYA_MINMATCHLEN,
                            GT_UNDEF_ULONG,
                            GT_UNDEF_ULONG,
                            GT_UNDEF_ULONG,
                            GT_UNDEF_ULONG,
                            GT_UNDEF_ULONG,
                            translationtable,
                            true,
                            false,
                            false,
                            true);

      /* begin XML comment */
      if (out->xmlout)
        gt_file_xprintf(out->outfp, "<!--\n");

      if (callvmatch(vmatch_info->argc,
                     (const char**) vmatch_info->argv,
                     NULL,
                     NULL,
                     "maskpolyAtailsandcreateindex",
                     NULL,
                     out->showverboseVM,
                     fp,
                     &selectbundle,
                     &vmatch_info->virtualtree,
                     &vmatch_info->queryvirtualtree,
                     &vmatch_info->sixframeofqueryvirtualtree,
                     &vmatch_info->dnavirtualtree)) {
        fprintf(stderr,"%s\n", messagespace());
        exit(EXIT_FAILURE);
      }

      /* end XML comment */
      if (out->xmlout)
        gt_file_xprintf(out->outfp, "-->\n");

      /* masked reference file has been computed, close it */
      gt_fa_xfclose(fp);

      if (out->showverbose) {
        out->showverbose("the following masked reference sequence file has "
                         "been created:");
        out->showverbose(maskedfilename);
        out->showverbose("write the corresponding .info file");
      }

      /* free */
      gth_vmatch_info_delete(vmatch_info);

      /* write .info file */
      writeinfofile(maskedinfofilename);

      /* a new index for the masked sequence needs to be created */
      createindex = true;
    }
  }

  if (!had_err) {
    /* check if masked reference index has already been created */
    if (!createindex) {
      makeemptyvirtualtree(&virtualtree);
      if (mapvirtualtreeifyoucan(&virtualtree, maskedfilename, TISTAB | DESTAB))
      {
        createindex = true;
      }

      if (!createindex && !online && inverse) {
        /* mapping was successful, save numofcodes and totallength */
        numofcodes  = virtualtree.numofcodes;
        totallength = virtualtree.multiseq.totallength;
        tablesok = tablesaremappableandhavecorrectlength(maskedfilename,
                                                         SUFTAB | BWTTAB |
                                                         LCPTAB | BCKTAB |
                                                         STI1TAB, numofcodes,
                                                         totallength);
        if (!tablesok)
          createindex = true;
      }

      if (freevirtualtree(&virtualtree)) {
        fprintf(stderr,"%s\n", messagespace());
        exit(EXIT_FAILURE);
      }
    }
  }

  if (!had_err) {
    if (createindex)
    {
      if (out->showverbose) {
        out->showverbose("create index for the following file:");
        out->showverbose(maskedfilename);
      }

      /* create index with TransDNAX mapping */
      had_err = gthcallmkvtree(maskedfilename, progname,
                               gt_str_get(proteinsmap), true, true, false,
                               false, !online && inverse, out, err);

      if (!had_err && out->showverbose)
        out->showverbose("index created");
    }
    else {
      if (out->showverbose) {
        out->showverbose("the index for the following file exists already:");
        out->showverbose(maskedfilename);
      }
    }
  }

  return had_err;
}

static void determine_overall_maxlength(Multiseq *multiseq,
                                        GtUword *overall_maxlength)
{
  ExtremeAverageSequences extreme;

  gt_assert(multiseq);
  gt_assert(overall_maxlength);

  /* for determing the width which is used to format genomic sequence
     positions */
  if (calculateseqparm(multiseq, &extreme)) {
    fprintf(stderr,"%s\n", messagespace());
    exit(EXIT_FAILURE);
  }
  if (extreme.maxlength > *overall_maxlength)
    *overall_maxlength = extreme.maxlength;
}

int gthmakesureindexexists(const char *filename, const char *progname,
                           const char *proteinsmap, bool isreferencefile,
                           bool dnasuffix, bool oistab, bool completeindex,
                           GtUword *overall_maxlength,
                           GthAlphatype alphatype, GthOutput *out, GtError *err)
{
  GtStr *indexname,
        *tisfilename;
  Virtualtree virtualtree;
  bool createindex, tablesok;
  Uint numofcodes = 0, totallength = 0, demand = TISTAB | DESTAB;
  int had_err = 0;

  gt_error_check(err);
  /* valid alphabet type */
  gt_assert(alphatype == DNA_ALPHA ||alphatype == PROTEIN_ALPHA);
  /* proteinsmap is defined iff reference files are preprocessed */
  gt_assert((isreferencefile && proteinsmap) ||
            !(isreferencefile && proteinsmap));

  if (oistab)
    demand |= OISTAB;

  createindex = false;

  if (out && out->showverbose) {
    out->showverbose("check the following file for index:");
    out->showverbose(filename);
  }

  /* init virtual tree */
  makeemptyvirtualtree(&virtualtree);

  indexname = gt_str_new_cstr(filename);
  gt_str_append_char(indexname, '.');
  if (alphatype == DNA_ALPHA) {
    gt_str_append_cstr(indexname, dnasuffix ? DNASUFFIX : POLYASUFFIX);
  }
  else if (alphatype == PROTEIN_ALPHA)
    gt_str_append_cstr(indexname, proteinsmap);

  tisfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tisfilename, ".tis");

  if (!gt_file_exists(gt_str_get(tisfilename)) ||
      gt_file_is_newer(filename, gt_str_get(tisfilename)) ||
      mapvirtualtreeifyoucan(&virtualtree, gt_str_get(indexname), demand)) {
    createindex = true;
  }
  if (!createindex) {
    /* mapping was successful, save numofcodes and totallength */
    numofcodes  = virtualtree.numofcodes;
    totallength = virtualtree.multiseq.totallength;

    /* check if mapped sequence has DNA alphabet */
    if (alphatype == DNA_ALPHA && !isstrictdnaalphabet(&virtualtree.alpha)) {
      gt_error_set(err, "alphabet type of index %s has to be DNA",
                   gt_str_get(indexname));
      had_err = -1;
    }

    /* if this is a genomic index, store the extreme values */
    if (!had_err && overall_maxlength && !isreferencefile)
      determine_overall_maxlength(&virtualtree.multiseq, overall_maxlength);
  }
  if (freevirtualtree(&virtualtree)) {
    fprintf(stderr,"%s\n", messagespace());
    exit(EXIT_FAILURE);
  }

  if (!had_err && !createindex && completeindex) {
    /* in this case check if the additional tables also exist. */
    tablesok = tablesaremappableandhavecorrectlength(gt_str_get(indexname),
                                                     SUFTAB | BWTTAB | LCPTAB |
                                                     BCKTAB | STI1TAB,
                                                     numofcodes, totallength);
    if (!tablesok)
      createindex = true;
  }

  if (!had_err) {
    if (createindex) {
      /* index does not exist */
      if (out && out->showverbose) {
        out->showverbose("index does not exist, is corrupted, or not "
                         "up-to-date => create it");
      }

      if (alphatype == DNA_ALPHA) {
        if (out && out->showverbose)
          out->showverbose("create DNA index");

        /* create DNA index */
        if (gthcallmkvtree(filename, progname, proteinsmap, true, !dnasuffix,
                           true, oistab, completeindex, out, err)) {
          had_err = -1;
        }
      }
      else {
        if (out && out->showverbose)
          out->showverbose("create protein index");

        /* create protein index */
        if (gthcallmkvtree(filename, progname, proteinsmap, false, false, true,
                           oistab, completeindex, out, err)) {
          had_err = -1;
        }
      }

      if (!had_err && out && out->showverbose)
        out->showverbose("index created");

      /* if this is a genomic index, store the extreme values */
      if (!had_err && overall_maxlength && !isreferencefile) {
        if (mapvirtualtreeifyoucan(&virtualtree, gt_str_get(indexname),
                                   demand)) {
          fprintf(stderr,"%s\n", messagespace());
          exit(EXIT_FAILURE);
        }
        determine_overall_maxlength(&virtualtree.multiseq, overall_maxlength);
        if (freevirtualtree(&virtualtree)) {
          fprintf(stderr,"%s\n", messagespace());
          exit(EXIT_FAILURE);
        }
      }
    }
    else { /* index exists */
      if (out && out->showverbose)
        out->showverbose("index exists");
    }
  }

  gt_str_delete(tisfilename);
  gt_str_delete(indexname);

  return had_err;
}

static int preprocessgenomicfile(const char* genomicfilename,
                                 bool gthconsensus, bool noautoindex,
                                 bool maskpolyAtails, bool online,
                                 bool inverse, bool skipindexcheck,
                                 GthAlphatype overall_alphatype,
                                 const char *progname,
                                 GtUword *overall_maxlength,
                                 GthOutput *out, GtError *err)
{
  gt_error_check(err);
  /* we cannot skip this step even if <skipindexcheck> is true, because
     <overall_maxlength> needs to be determined */
  if (gthmakesureindexexists(genomicfilename, progname, NULL, false, true, true,
                             (!gthconsensus && !noautoindex &&
                              !maskpolyAtails && !online && !inverse &&
                              overall_alphatype != PROTEIN_ALPHA),
                             overall_maxlength , DNA_ALPHA, out, err)) {
    return -1;
  }
  if (!skipindexcheck) {
    if (maskpolyAtails) {
      if (gthmakesureindexexists(genomicfilename, progname, NULL, false, false,
                                 false, (!online && !inverse),
                                 overall_maxlength, DNA_ALPHA, out, err)) {
        return -1;
      }
    }
  }
  return 0;
}

static int preprocessreferencefile(const char *referencefilename,
                                   GtStr *proteinsmap, bool gthconsensus,
                                   bool noautoindex, bool maskpolyAtails,
                                   bool online, bool inverse,
                                   bool skipindexcheck, const char *progname,
                                   GthAlphatype alphatype,
                                   Uint translationtable, GthOutput *out,
                                   GtError *err)
{
  gt_error_check(err);
  if (!skipindexcheck) {
    if (gthmakesureindexexists(referencefilename, progname,
                               gt_str_get(proteinsmap), true, true, true,
                               ((!gthconsensus && !noautoindex && !online) &&
                                ((!maskpolyAtails && inverse &&
                                  alphatype == DNA_ALPHA) ||
                                 (alphatype == PROTEIN_ALPHA))), NULL,
                               alphatype, out, err)) {
      return -1;
    }
    if (alphatype == DNA_ALPHA && maskpolyAtails) {
      if (maskpolyAtailsandcreateindex(referencefilename, progname, proteinsmap,
                                       online, inverse, translationtable, out,
                                       err)) {
        return -1;
      }
    }
  }
  return 0;
}

static int preprocessinputfiles(GthInput *input,
                                bool gthconsensus, bool noautoindex,
                                bool maskpolyAtails, bool online, bool inverse,
                                bool skipindexcheck, const char *progname,
                                GtStr *proteinsmap,
                                GtUword *overall_maxlength,
                                Uint translationtable, GthOutput *out,
                                GtError *err)
{
  GtUword i;

  gt_error_check(err);

  /* make sure the indices of all genomic input files exist */
  if (out->showverbose) {
    out->showverbose("make sure the necessary indices of all genomic input "
                     "files exist");
  }
  for (i = 0; i < gth_input_num_of_gen_files(input); i++) {
    if (preprocessgenomicfile(gth_input_get_genomic_filename(input, i),
                              gthconsensus, noautoindex, maskpolyAtails, online,
                              inverse, skipindexcheck,
                              gth_input_overall_alphatype(input), progname,
                              overall_maxlength, out, err)) {
      return -1;
    }
  }

  /* make sure the indices of all reference input files exist */
  if (out->showverbose) {
    out->showverbose("make sure the necessary indices off all reference input "
                     "files exist");
  }
  for (i = 0; i < gth_input_num_of_ref_files(input); i++) {
    if (preprocessreferencefile(gth_input_get_reference_filename(input, i),
                                proteinsmap, gthconsensus, noautoindex,
                                maskpolyAtails, online, inverse, skipindexcheck,
                                progname,
                                gth_input_get_alphatype(input, i),
                                translationtable, out, err)) {
      return -1;
    }
  }

  return 0;
}

int gthpreprocessinputfiles(GthInput *input, bool gthconsensus,
                            bool noautoindex, bool skipindexcheck,
                            bool maskpolyAtails, bool online, bool inverse,
                            const char *progname,
                            unsigned int translationtable, GthOutput *out,
                            GtError *err)
{
  GtUword overall_maxlength = 0;

  gt_error_check(err);

  if (out->showverbose)
    out->showverbose("make sure all necessary indices exist");

  /* make sure all indices exist */
  if (preprocessinputfiles(input, gthconsensus, noautoindex, maskpolyAtails,
                           online, inverse, skipindexcheck, progname,
                           gth_input_proteinsmap(input),
                           &overall_maxlength, translationtable, out, err)) {
    return -1;
  }

  /* set the width which is used to format genomic sequence positions */
  if (out->gs2out) {
    if (overall_maxlength < 1000000)
      out->widthforgenpos = 6;
    else
      out->widthforgenpos = 10;
  }
  else if (gth_input_num_of_gen_files(input)) { /* we have at least one genomic
                                                   file */
    gt_assert(overall_maxlength); /* overall max length has been determined */
    out->widthforgenpos = floor(log10(overall_maxlength)) + 1;
  }

  return 0;
}
