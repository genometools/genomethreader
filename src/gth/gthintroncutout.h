/*
  Copyright (c) 2004-2007 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2004-2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef GTHINTRONCUTOUT_H
#define GTHINTRONCUTOUT_H

#include <stdbool.h>

typedef struct
{
  bool introncutout;                 /* use intron cutouts */
  unsigned int autoicmaxmatrixsize,  /* maximum matrix size (in megabytes) for
                                        automatic intron cutout technique
                                        (enabled if > 0) */
               icinitialdelta,       /* initial intron cutout delta */
               iciterations,         /* number of intron cutout iterations */
               icdeltaincrease,      /* the delta increase during every
                                        iteration */
               icminremintronlength; /* intron cutout minimum remaining intron
                                        length */
} Introncutoutinfo;

#endif
