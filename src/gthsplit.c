#include "gth_config.h"
#include "core/tooldriver_api.h"
#include "gth/gt_gthsplit.h"
#include "libgenomethreader/gthpre.h"
#include "libgenomethreader/gthversionfunc.h"
#include "libgenomethreader/gthvmatch.h"
#include "libgenomethreader/seq_con_multiseq.h"

static int gthsplit(int argc, const char **argv, GtError *err)
{
  static const GthPlugins plugins = { gthpreprocessinputfiles,
                                      gth_seq_con_multiseq_new,
                                      gth_vmatch_matcher_arguments_new,
                                      gth_vmatch_matcher_arguments_delete,
                                      gth_vmatch_matcher_runner,
                                      GTH_VERSION,
                                      gthversionfunc };
  gt_error_check(err);
  return gt_gthsplit(argc, argv, &plugins, err);
}

int main(int argc, char *argv[])
{
  return gt_tooldriver(gthsplit, argc, argv);
}
