#include "core/ma_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "gth/align_dna.h"
#include "gth/input.h"
#include "gth/sa.h"
#include "gth/splice_site_model.h"
#include "libgenomethreader/gthpre.h"
#include "libgenomethreader/gt_align_dna.h"
#include "libgenomethreader/gthversionfunc.h"
#include "libgenomethreader/seq_con_multiseq.h"

typedef struct {
  bool gen_reverse;
  GtRange gen_range;
} AlignDNAArguments;

static void* gt_align_dna_arguments_new(void)
{
  return gt_calloc(1, sizeof (AlignDNAArguments));
}

static void gt_align_dna_arguments_delete(void *tool_arguments)
{
  AlignDNAArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_free(arguments);
}

static GtOptionParser* gt_align_dna_option_parser_new(void *tool_arguments)
{
  AlignDNAArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);
  op = gt_option_parser_new("[option ...] genomic_file cDNA_file",
                            "Spliced align genomic_file with cDNA_file.");
  option = gt_option_new_bool("genreverse", "align reverse strand of genomic "
                              "sequence", &arguments->gen_reverse, false);
  gt_option_parser_add_option(op, option);
  option = gt_option_new_range("genrange", "set aligned genomic range",
                               &arguments->gen_range, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_parser_set_min_max_args(op, 2, 2);
  gt_option_parser_set_version_func(op, gthversionfunc);
  return op;
}

static int gt_align_dna_runner(GT_UNUSED int argc, const char **argv,
                               int parsed_args, void *tool_arguments,
                               GtError *err)
{
  AlignDNAArguments *arguments = tool_arguments;
  const GtUword gen_file_num = 0, gen_seq_num = 0,
                      ref_file_num = 0, ref_seq_num = 0;
  GthSpliceSiteModel *splice_site_model = NULL;
  const char *genomic_file, *cdna_file;
  GthInput *input;
  GthSA *sa = NULL;
  int had_err;

  gt_error_check(err);
  gt_assert(arguments);

  genomic_file = argv[parsed_args];
  cdna_file = argv[parsed_args+1];

  input = gth_input_new(gthpreprocessinputfiles, gth_seq_con_multiseq_new);
  gth_input_add_genomic_file(input, genomic_file);
  gth_input_add_cdna_file(input, cdna_file);
  if (arguments->gen_range.start == GT_UNDEF_ULONG) {
    gth_input_load_genomic_file(input, gen_file_num, true);
    arguments->gen_range = gth_input_get_genomic_range(input, gen_file_num,
                                                       gen_seq_num);
  }

  had_err = gth_input_make_indices(input, argv[0], err);

  if (!had_err) {
    splice_site_model = gth_splice_site_model_new();
    /* XXX */
    had_err = gth_splice_site_model_load_bssm(splice_site_model, "human", err);
  }

  if (!had_err) {
    sa = gth_align_dna_simple(input, &arguments->gen_range, gen_file_num,
                              gen_seq_num, !arguments->gen_reverse,
                              ref_file_num, ref_seq_num, splice_site_model);
    if (!sa) {
      gt_error_set(err, "out of memory");
      had_err = -1;
    }
  }

  gth_input_delete_current(input);

  if (!had_err)
    gth_sa_show(sa, input, NULL);

  gth_sa_delete(sa);
  gth_splice_site_model_delete(splice_site_model);
  gth_input_delete_complete(input);

  return had_err;
}

GtTool* gt_align_dna(void)
{
  return gt_tool_new(gt_align_dna_arguments_new,
                     gt_align_dna_arguments_delete,
                     gt_align_dna_option_parser_new,
                     NULL,
                     gt_align_dna_runner);
}
