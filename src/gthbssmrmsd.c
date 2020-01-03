#include "gth_config.h"
#include "core/init_api.h"
#include "core/tooldriver_api.h"
#include "gth/gt_gthbssmrmsd.h"
#include "libgenomethreader/gthversionfunc.h"

int main(int argc, char *argv[])
{
  return gt_toolobjdriver(gt_gthbssmrmsd, gthversionfunc, argc, argv);
}
