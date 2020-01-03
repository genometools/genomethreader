#include "core/tooldriver_api.h"
#include "gth/gt_gthbssmfileinfo.h"
#include "libgenomethreader/gthversionfunc.h"

static int gthbssmfileinfo(int argc, const char **argv, GtError *err)
{
  gt_error_check(err);
  return gt_gthbssmfileinfo_with_version_func(argc, argv, gthversionfunc, err);
}

int main(int argc, char *argv[])
{
  return gt_tooldriver(gthbssmfileinfo, argc, argv);
}
