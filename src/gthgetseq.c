#include "core/tooldriver_api.h"
#include "libgenomethreader/gt_gthgetseq.h"

int main(int argc, char *argv[])
{
  return gt_tooldriver(gt_gthgetseq, argc, argv);
}
