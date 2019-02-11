#include "core/tooldriver.h"
#include "libgenomethreader/gt_align_dna.h"
#include "libgenomethreader/gthversionfunc.h"

#ifndef NOLICENSEMANAGER
#include "licensemanager.h"
#endif

int main(int argc, char *argv[])
{
#ifndef NOLICENSEMANAGER
  return gt_toolobjdriver_with_license(gt_align_dna, gthversionfunc, argc, argv,
                                       (GtLicense**) &lm_license,
                                       lm_license_new_gth,
                                       lm_license_delete_gt);
#else
  return gt_toolobjdriver(gt_align_dna, gthversionfunc, argc, argv);
#endif
}
