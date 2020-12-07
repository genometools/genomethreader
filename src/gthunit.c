#include "core/hashmap_api.h"
#include "core/init_api.h"
#include "core/unit_testing_api.h"
#include "core/unused_api.h"
#include "libgenomethreader/jt.h"

static GtHashmap* gth_unit_tests(void)
{
  GtHashmap *unit_tests = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);

  gt_hashmap_add(unit_tests, "jump table class", gth_jump_table_unit_test);

  return unit_tests;
}

int main(GT_UNUSED int argc, GT_UNUSED char *argv[])
{
  int test_err = 0, had_err;
  GtHashmap *unit_tests;
  GtError *err;
  gt_lib_init();
  err = gt_error_new();
  unit_tests = gth_unit_tests();
  had_err = gt_hashmap_foreach_in_key_order(unit_tests, gt_unit_test_run,
                                            &test_err, err);
  gt_assert(!had_err); /* cannot happen, gt_unit_test_run() is sane */
  gt_hashmap_delete(unit_tests);
  gt_error_delete(err);
  if (gt_lib_clean())
    return GT_EXIT_PROGRAMMING_ERROR; /* programmer error */
  if (test_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
