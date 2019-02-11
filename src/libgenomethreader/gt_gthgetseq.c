#include "core/bool_matrix.h"
#include "core/cstr_api.h"
#include "core/fasta.h"
#include "core/file.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/unused_api.h"
#include "gth/default.h"
#include "gth/intermediate.h"
#include "gth/sa_filter.h"
#include "libgenomethreader/gthpre.h"
#include "libgenomethreader/gthversionfunc.h"
#include "libgenomethreader/gt_gthgetseq.h"
#include "libgenomethreader/seq_con_multiseq.h"
#include "alphadef.h"
#include "multidef.h"
#include "multiseq-adv.pr"

#define GETCDNA_OPT_CSTR               "getcdna"
#define GETCDNACOMPLEMENT_OPT_CSTR     "getcdnacomp"
#define GETPROTEIN_OPT_CSTR            "getprotein"
#define GETPROTEINCOMPLEMENT_OPT_CSTR  "getproteincomp"
#define GETGENOMIC_OPT_CSTR            "getgenomic"
#define GETGENOMICCOMPLEMENT_OPT_CSTR  "getgenomiccomp"

#define DEFAULT_WIDTH  60

typedef struct {
  bool get_cdna_seq,
       get_cdna_seq_complement,
       get_protein_seq,
       get_protein_seq_complement,
       get_genomic_seq,
       get_genomic_seq_complement,
       process_cdna_seqs,
       process_protein_seqs,
       process_reference_seqs_normal,
       process_reference_seqs_complement,
       process_genomic_seqs;
  GtFileMode genfilemode;
  GtStrArray *consensusfiles;
  GthSAFilter *sa_filter;
  char *progname;
} Gthgetseqinfo;

typedef struct {
  Gthgetseqinfo *gthgetseqinfo;
  GthInput *input;
  GtUword current_genomicfile,
                current_referencefile;
  GtBoolMatrix *genseq_boolmatrix,
               *refseq_boolmatrix;
} Extractseq_data;

static int getusedseq(void *data, GthSA *sa,
                      GT_UNUSED const char *outputfilename, GtError *err)
{
  Extractseq_data *getseq_data = (Extractseq_data*) data;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(sa);
  gt_assert(getseq_data->gthgetseqinfo->sa_filter);

  /* filter before we do any further processing */
  if (gth_sa_filter_filter_sa(getseq_data->gthgetseqinfo->sa_filter, sa)) {
    /* and free the spliced alignment afterwards */
    gth_sa_delete(sa);
    /* discard */
    return 0;
  }

  /* make sure the necessary indices exist */
  if (getseq_data->current_referencefile <
      gth_input_num_of_ref_files(getseq_data->input)) {
    /* exactly one missing index */
    gt_assert(gth_input_num_of_ref_files(getseq_data->input) -
              getseq_data->current_referencefile == 1);

    if ((gth_input_get_alphatype(getseq_data->input,
                                    getseq_data->current_referencefile)
         == DNA_ALPHA &&
         getseq_data->gthgetseqinfo->process_cdna_seqs) ||
        (gth_input_get_alphatype(getseq_data->input,
                                    getseq_data->current_referencefile)
         == PROTEIN_ALPHA &&
         getseq_data->gthgetseqinfo->process_protein_seqs)) {
      if (gthmakesureindexexists(
            gth_input_get_reference_filename(getseq_data->input,
                                            getseq_data->current_referencefile),
                                 getseq_data->gthgetseqinfo->progname,
                                 GTH_DEFAULT_PROTEINSMAP, true, true, true,
                                 false, NULL,
                                 gth_input_get_alphatype(getseq_data->input,
                                            getseq_data->current_referencefile),
                                 NULL, err)) {
        had_err = -1;
      }
    }
    getseq_data->current_referencefile++;
  }
  if (!had_err && getseq_data->current_genomicfile <
                  gth_input_num_of_gen_files(getseq_data->input)) {
    /* exactly one missing index */
    gt_assert(gth_input_num_of_gen_files(getseq_data->input) -
              getseq_data->current_genomicfile == 1);

    if (getseq_data->gthgetseqinfo->process_genomic_seqs) {
      if (gthmakesureindexexists(
            gth_input_get_genomic_filename(getseq_data->input,
                                           getseq_data->current_genomicfile),
                                 getseq_data->gthgetseqinfo->progname,
                                 GTH_DEFAULT_PROTEINSMAP, false, true, true,
                                 false, NULL, DNA_ALPHA, NULL, err)) {
        had_err = -1;
      }
    }
    getseq_data->current_genomicfile++;
  }

  if (!had_err &&
      (
       (getseq_data->gthgetseqinfo->process_cdna_seqs &&
        gth_input_get_alphatype(getseq_data->input,
                                getseq_data->current_referencefile-1)
        == DNA_ALPHA) ||
       (getseq_data->gthgetseqinfo->process_protein_seqs &&
        gth_input_get_alphatype(getseq_data->input,
                                getseq_data->current_referencefile-1)
        == PROTEIN_ALPHA))) {
    gt_bool_matrix_set(getseq_data->refseq_boolmatrix, gth_sa_ref_file_num(sa),
                       gth_sa_ref_seq_num(sa), true);
  }
  if (!had_err && getseq_data->gthgetseqinfo->process_genomic_seqs) {
    gt_bool_matrix_set(getseq_data->genseq_boolmatrix, gth_sa_gen_file_num(sa),
                       gth_sa_gen_seq_num(sa), true);
  }

  /* and free the spliced alignment afterwards */
  gth_sa_delete(sa);

  return had_err;
}

static void initGthgetseqinfo(Gthgetseqinfo *gthgetseqinfo,
                              const char *progname)
{
  gthgetseqinfo->genfilemode                = GT_FILE_MODE_UNCOMPRESSED;
  gthgetseqinfo->sa_filter                  = gth_sa_filter_new();
  gthgetseqinfo->consensusfiles             = gt_str_array_new();
  gthgetseqinfo->progname                   = gt_cstr_dup(progname);
}

static void freeGthgetseqinfo(Gthgetseqinfo *gthgetseqinfo)
{
  if (!gthgetseqinfo) return;
  gt_str_array_delete(gthgetseqinfo->consensusfiles);
  gth_sa_filter_delete(gthgetseqinfo->sa_filter);
  gt_free(gthgetseqinfo->progname);
}

static int gthgetseq_parse_options(int *parsed_args,
                                   Gthgetseqinfo *gthgetseqinfo,
                                   int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *optgetcdna, *optgetcdnacomplement, *optgetprotein,
         *optgetproteincomplement, *optgetgenomic, *optgetgenomiccomplement,
         *optgzip, *optbzip2;
  bool gzip, bzip2;
  GtOPrval oprval;

  gt_error_check(err);

  op = gt_option_parser_new("-getcdna | -getprotein | -getgenomic [option ...] "
                          "[file ...]", "Get FASTA sequences from "
                          "GenomeThreader files containing intermediate "
                          "results.\nThe sequences are shown on stdout.");

  /* specify all options with a corresponding help-text */
  optgetcdna = gt_option_new_bool(GETCDNA_OPT_CSTR, "get cDNA/EST sequences",
                               &gthgetseqinfo->get_cdna_seq, false);
  gt_option_parser_add_option(op, optgetcdna);

  optgetcdnacomplement = gt_option_new_bool(GETCDNACOMPLEMENT_OPT_CSTR,
                                         "get complement of cDNA/EST sequences",
                                         &gthgetseqinfo
                                         ->get_cdna_seq_complement, false);
  gt_option_parser_add_option(op, optgetcdnacomplement);

  optgetprotein = gt_option_new_bool(GETPROTEIN_OPT_CSTR,
                                     "get protein sequences",
                                     &gthgetseqinfo->get_protein_seq, false);
  gt_option_parser_add_option(op, optgetprotein);

  optgetproteincomplement = gt_option_new_bool(GETPROTEINCOMPLEMENT_OPT_CSTR,
                                            "get complement of protein "
                                            "sequences",
                                            &gthgetseqinfo
                                            ->get_protein_seq_complement,
                                            false);
  gt_option_parser_add_option(op, optgetproteincomplement);

  optgetgenomic = gt_option_new_bool(GETGENOMIC_OPT_CSTR,
                                     "get genomic sequences",
                                     &gthgetseqinfo->get_genomic_seq, false);
  gt_option_parser_add_option(op, optgetgenomic);

  optgetgenomiccomplement = gt_option_new_bool(GETGENOMICCOMPLEMENT_OPT_CSTR,
                                            "get complement of genomic "
                                            "sequences",
                                            &gthgetseqinfo
                                            ->get_genomic_seq_complement,
                                            false);
  gt_option_parser_add_option(op, optgetgenomiccomplement);

  /* add sa_filter options */
  gth_sa_filter_register_options(op, gthgetseqinfo->sa_filter, false);

  optgzip = gt_option_new_bool("gzip", "gzip compressed input file(s)", &gzip,
                               false);
  gt_option_parser_add_option(op, optgzip);

  optbzip2 = gt_option_new_bool("bzip2", "bzip2 compressed input file(s)",
                                &bzip2, false);
  gt_option_parser_add_option(op, optbzip2);

  /* option exclusions */
  gt_option_exclude(optgetcdna, optgetcdnacomplement);
  gt_option_exclude(optgetprotein, optgetproteincomplement);
  gt_option_exclude(optgetgenomic, optgetgenomiccomplement);
  gt_option_exclude(optgzip, optbzip2);

  gt_option_parser_set_mail_address(op, "<gordon@gremme.org>");
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gthversionfunc,
                                  err);

  if (oprval == GT_OPTION_PARSER_OK && gzip)
    gthgetseqinfo->genfilemode = GT_FILE_MODE_GZIP;
  if (oprval == GT_OPTION_PARSER_OK && bzip2)
    gthgetseqinfo->genfilemode = GT_FILE_MODE_BZIP2;

  if (oprval == GT_OPTION_PARSER_OK &&
      !(gt_option_is_set(optgetcdna)              ||
        gt_option_is_set(optgetcdnacomplement)    ||
        gt_option_is_set(optgetprotein)           ||
        gt_option_is_set(optgetproteincomplement) ||
        gt_option_is_set(optgetgenomic)           ||
        gt_option_is_set(optgetgenomiccomplement))) {
    gt_error_set(err, "either option -%s, -%s, -%s, -%s, -%s, or -%s is "
                  "mandatory", GETCDNA_OPT_CSTR, GETCDNACOMPLEMENT_OPT_CSTR,
                  GETPROTEIN_OPT_CSTR, GETPROTEINCOMPLEMENT_OPT_CSTR,
                  GETGENOMIC_OPT_CSTR, GETGENOMICCOMPLEMENT_OPT_CSTR);
    oprval = GT_OPTION_PARSER_ERROR;
  }

  if (oprval == GT_OPTION_PARSER_OK) {
    gthgetseqinfo->process_cdna_seqs = gthgetseqinfo->get_cdna_seq ||
                                        gthgetseqinfo->get_cdna_seq_complement;
    gthgetseqinfo->process_protein_seqs = gthgetseqinfo->get_protein_seq ||
                                          gthgetseqinfo
                                          ->get_protein_seq_complement;
    gthgetseqinfo->process_reference_seqs_normal = gthgetseqinfo
                                                   ->get_cdna_seq ||
                                                   gthgetseqinfo
                                                   ->get_protein_seq;
    gthgetseqinfo
    ->process_reference_seqs_complement = gthgetseqinfo
                                          ->get_cdna_seq_complement ||
                                          gthgetseqinfo
                                          ->get_protein_seq_complement;
    gthgetseqinfo->process_genomic_seqs = gthgetseqinfo->get_genomic_seq ||
                                          gthgetseqinfo
                                          ->get_genomic_seq_complement;
  }

  /* save consensus files */
  if (oprval == GT_OPTION_PARSER_OK) {
    while (*parsed_args < argc) {
      gt_str_array_add_cstr(gthgetseqinfo->consensusfiles, argv[*parsed_args]);
      (*parsed_args)++;
    }
  }

  if (oprval == GT_OPTION_PARSER_OK &&
      !gt_str_array_size(gthgetseqinfo->consensusfiles) &&
      (gt_option_is_set(optgzip) || gt_option_is_set(optbzip2))) {
    gt_error_set(err, "to use compression, at least on input file has to be "
                   "supplied");
    oprval = GT_OPTION_PARSER_ERROR;
  }

  gt_option_parser_delete(op);

  return oprval;
}

static int gthgetseq_process_files(Gthgetseqinfo *gthgetseqinfo, GtError *err)
{
  Extractseq_data getseq_data;
  GthInput *input;
  bool seq_status;
  GtUword i, j;
  GthSeqCon *seq_con;
  GtStr *desc;
  int had_err = 0;

  gt_error_check(err);

  /* initialization */
  desc = gt_str_new();
  input= gth_input_new(gthpreprocessinputfiles, gth_seq_con_multiseq_new);
  gth_input_set_forward_only(input);
  getseq_data.gthgetseqinfo         = gthgetseqinfo;
  getseq_data.input                 = input;
  getseq_data.current_genomicfile   = 0;
  getseq_data.current_referencefile = 0;
  getseq_data.genseq_boolmatrix     = gt_bool_matrix_new();
  getseq_data.refseq_boolmatrix     = gt_bool_matrix_new();

  /* get used sequences */
  if (!had_err) {
    had_err = gth_process_intermediate_files(input,
                                             gthgetseqinfo->consensusfiles,
                                             getusedseq, &getseq_data, NULL,
                                             err);
  }

  /* print reference sequences */
  if (!had_err) {
    if (gthgetseqinfo->process_reference_seqs_normal ||
        gthgetseqinfo->process_reference_seqs_complement) {
      /* iterate over reference files */
      for (i = 0;
           !had_err && i < gth_input_num_of_ref_files(getseq_data.input);
           i++) {
        if ((gthgetseqinfo->process_cdna_seqs &&
             gth_input_get_alphatype(getseq_data.input, i) == DNA_ALPHA) ||
            (gthgetseqinfo->process_protein_seqs &&
             gth_input_get_alphatype(getseq_data.input, i) == PROTEIN_ALPHA)) {
          gth_input_load_reference_file(input, i, false);
          seq_con = gth_input_current_ref_seq_con(input);
          /* iterate over reference sequences */
          for (j = 0;
               !had_err && j < gth_input_num_of_ref_seqs(input, i);
               j++) {
            seq_status =
              gt_bool_matrix_get(getseq_data.refseq_boolmatrix, i, j);
            if ((gthgetseqinfo->process_reference_seqs_normal && seq_status) ||
                (gthgetseqinfo->process_reference_seqs_complement &&
                 !seq_status)) {
              gt_str_reset(desc);
              gth_seq_con_get_description(seq_con, j, desc);
              gt_fasta_show_entry(gt_str_get(desc),
                                  (const char*)
                                  gth_seq_con_get_orig_seq(seq_con, j),
                                  gth_seq_con_get_length(seq_con, j),
                                  DEFAULT_WIDTH, NULL);
            }
          }
        }
      }
    }
  }

  if (!had_err) {
    if (gthgetseqinfo->process_genomic_seqs) {
      /* iterate over genomic files */
      for (i = 0;
           !had_err &&
           i < gth_input_num_of_gen_files(getseq_data.input);
           i++) {
        gth_input_load_genomic_file(input, i, false);
        seq_con = gth_input_current_gen_seq_con(input);
        /* iterate over genomic sequences */
        for (j = 0; !had_err && j < gth_input_num_of_gen_seqs(input, i); j++) {
          seq_status =
            gt_bool_matrix_get(getseq_data.genseq_boolmatrix, i, j);
          if ((gthgetseqinfo->get_genomic_seq && seq_status) ||
              (gthgetseqinfo->get_genomic_seq_complement && !seq_status)) {
            gt_str_reset(desc);
            gth_seq_con_get_description(seq_con, j, desc);
            gt_fasta_show_entry(gt_str_get(desc),
                                (const char*)
                                gth_seq_con_get_orig_seq(seq_con, j),
                                gth_seq_con_get_length(seq_con, j),
                                DEFAULT_WIDTH, NULL);
          }
        }
      }
    }
  }

  /* free */
  gth_input_delete_complete(input);
  gt_bool_matrix_delete(getseq_data.genseq_boolmatrix);
  gt_bool_matrix_delete(getseq_data.refseq_boolmatrix);
  gt_str_delete(desc);

  return had_err;
}

int gt_gthgetseq(int argc, const char **argv, GtError *err)
{
  Gthgetseqinfo gthgetseqinfo;
  int had_err, parsed_args;
  gt_error_check(err);

  /* init data structures */
  initGthgetseqinfo(&gthgetseqinfo, argv[0]);

  switch (gthgetseq_parse_options(&parsed_args, &gthgetseqinfo, argc, argv,
                                  err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR:
      freeGthgetseqinfo(&gthgetseqinfo);
      return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT:
      freeGthgetseqinfo(&gthgetseqinfo);
      return 0;
  }
  gt_assert(parsed_args == argc);

  /* process files */
  had_err = gthgetseq_process_files(&gthgetseqinfo, err);

  /* free */
  freeGthgetseqinfo(&gthgetseqinfo);

  return had_err;
}
