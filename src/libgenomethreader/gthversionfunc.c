#include <stdio.h>
#include "core/version_api.h"
#include "gth_config.h"

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
  printf("Email: gordon@gremme.org\n\n");
  printf("Used compiler: %s\n", GTH_CC);
  printf("Compile flags: %s\n", GTH_CFLAGS);
}
