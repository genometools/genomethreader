#include "core/tooldriver_api.h"
#include "gth/gt_gthbssmbuild.h"
#include "libgenomethreader/gthversionfunc.h"

static int gthbssmbuild(int argc, const char **argv, GtError *err)
{
  gt_error_check(err);
  return gt_gthbssmbuild_with_version_func(argc, argv, gthversionfunc, err);
}

int main(int argc, char *argv[])
{
  return gt_tooldriver(gthbssmbuild, argc, argv);
}
