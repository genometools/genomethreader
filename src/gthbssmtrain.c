#include "gth_config.h"
#include "core/tooldriver.h"
#include "gth/gt_gthbssmtrain.h"
#include "libgenomethreader/gthversionfunc.h"

int main(int argc, char *argv[])
{
  return gt_toolobjdriver(gt_gthbssmtrain, gthversionfunc, argc, argv);
}
