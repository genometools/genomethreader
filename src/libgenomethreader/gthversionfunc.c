#include <stdio.h>
#include "core/version_api.h"
#include "gth_config.h"

#ifndef NOLICENSEMANAGER
#include "licensemanager.h"
#include "zlm.h"
#endif

void gthversionfunc(const char *progname)
{
  printf("%s (GenomeThreader) %s\n", progname, GTH_VERSION);
  printf("\n");
  printf("libgenometools-%s:\n", gt_version());
  printf("Copyright (c) 2003-2018 G. Gremme, S. Steinbiss, S. Kurtz, and "
         "CONTRIBUTORS\n");
  printf("Copyright (c) 2003-2018 Center for Bioinformatics, University of "
         "Hamburg\n\n");
  printf("libgenomethreader:\n");
  printf("Copyright (c) 2009-2018 Wikena GmbH\n\n");
  printf("libvmatch:\n");
  printf("Copyright (c) 2000-2017 LScSA-Software GmbH\n\n");
#ifndef NOLICENSEMANAGER
  printf("libzlm-%s:\n", zlm_version());
  printf("Copyright (c) 2013-2016 Wikena GmbH\n\n");
#endif
  printf("Email: gordon@gremme.org\n\n");
  printf("Used compiler: %s\n", GTH_CC);
  printf("Compile flags: %s\n", GTH_CFLAGS);
#ifndef NOLICENSEMANAGER
  if (lm_license) {
    putchar('\n');
    lm_license_show_info(lm_license);
  }
#endif
}
