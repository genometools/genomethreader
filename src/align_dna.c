#include "core/tooldriver_api.h"
#include "libgenomethreader/gt_align_dna.h"
#include "libgenomethreader/gthversionfunc.h"

int main(int argc, char *argv[])
{
  return gt_toolobjdriver(gt_align_dna, gthversionfunc, argc, argv);
}
