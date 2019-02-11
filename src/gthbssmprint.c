#include "core/tooldriver.h"
#include "gth/gt_gthbssmprint.h"
#include "libgenomethreader/gthversionfunc.h"

static int gthbssmprint(int argc, const char **argv, GtError *err)
{
  gt_error_check(err);
  return gt_gthbssmprint_with_version_func(argc, argv, gthversionfunc, err);
}

int main(int argc, char *argv[])
{
  return gt_tooldriver(gthbssmprint, argc, argv);
}
