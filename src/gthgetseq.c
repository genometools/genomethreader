#include "core/tooldriver.h"
#include "libgenomethreader/gt_gthgetseq.h"

#ifndef NOLICENSEMANAGER
#include "licensemanager.h"
#endif

int main(int argc, char *argv[])
{
#ifndef NOLICENSEMANAGER
  return gt_tooldriver_with_license(gt_gthgetseq, argc, argv,
                                    (GtLicense**) &lm_license,
                                    lm_license_new_gth, lm_license_delete_gt);
#else
  return gt_tooldriver(gt_gthgetseq, argc, argv);
#endif
}
