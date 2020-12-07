#include "gth_config.h"
#include "core/tooldriver_api.h"
#include "gth/gt_gth.h"
#include "libgenomethreader/gthpre.h"
#include "libgenomethreader/gthversionfunc.h"
#include "libgenomethreader/gthvmatch.h"
#include "libgenomethreader/jt.h"
#include "libgenomethreader/jt_align_dna.h"
#include "libgenomethreader/jt_align_protein.h"
#include "libgenomethreader/seq_con_multiseq.h"

static int gth(int argc, const char **argv, GtError *err)
{
  static const GthPlugins plugins = { gthpreprocessinputfiles,
                                      gth_seq_con_multiseq_new,
                                      gth_vmatch_matcher_arguments_new,
                                      gth_vmatch_matcher_arguments_delete,
                                      gth_vmatch_matcher_runner,
                                      GTH_VERSION,
                                      gthversionfunc,
                                      gth_jump_table_new,
                                      gth_jump_table_new_reverse,
                                      gth_jump_table_delete,
                                      gth_dna_complete_path_matrix_jt,
                                      gth_protein_complete_path_matrix_jt };
  gt_error_check(err);
  return gt_gth(argc, argv, &plugins, err);
}

int main(int argc, char *argv[])
{
  int rval;
  rval = gt_tooldriver(gth, argc, argv);
  return rval;
}
