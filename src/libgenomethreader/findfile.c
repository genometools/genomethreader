/*
  Copyright (c) 2020 Sascha Steinbiss <sascha@steinbiss.name>

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

#include "findfile.h"
#include "core/fileutils_api.h"
#include "core/compat_api.h"

int gth_find_file(const char *filename, const char *envname,
                  const char *dirname, GtStr *out, GtError *err)
{
  int had_err = 0;
  const char **defaultpath;
  static const char *defaultpaths[] = {
    "/usr/share/genomethreader",
    "/usr/local/share/genomethreader",
    "/usr/lib/genomethreader",
    "/usr/local/lib/genomethreader",
    NULL
  };

  gt_str_reset(out);
  /* check in directory given by envname */
  if (getenv(envname)) {
    gt_str_append_cstr(out, getenv(envname));
    gt_str_append_char(out, GT_PATH_SEPARATOR);
    gt_str_append_cstr(out, filename);
    if (gt_file_exists(gt_str_get(out))) {
      return 0;
    }
  }

  /* check for file relative to binary */
  gt_str_reset(out);
  had_err = gt_file_find_exec_in_path(out, gt_error_get_progname(err), NULL);
  if (!had_err) {
    gt_str_append_char(out, GT_PATH_SEPARATOR);
    gt_str_append_cstr(out, dirname);
    gt_str_append_char(out, GT_PATH_SEPARATOR);
    gt_str_append_cstr(out, filename);
    if (gt_file_exists(gt_str_get(out))) {
      return 0;
    }
  }

  /* check for file in set of default paths */
  for (defaultpath = defaultpaths; *defaultpath; defaultpath++) {
    gt_str_reset(out);
    gt_str_append_cstr(out, *defaultpath);
    gt_str_append_char(out, GT_PATH_SEPARATOR);
    gt_str_append_cstr(out, dirname);
    gt_str_append_char(out, GT_PATH_SEPARATOR);
    gt_str_append_cstr(out, filename);
    if (gt_file_exists(gt_str_get(out))) {
      return 0;
    }
  }

  gt_error_set(err, "could not find file '%s' in any of the search paths",
               filename);
  return -1;
}
