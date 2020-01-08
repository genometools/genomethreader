#include "core/bittab_api.h"
#include "core/divmodmul_api.h"
#include "core/ensure_api.h"
#include "core/ma_api.h"
#include "core/minmax_api.h"
#include "core/range_api.h"
#include "libgenomethreader/jt.h"

#ifndef NOLICENSEMANAGER
#include "licensemanager.h"
#endif

typedef struct {
  GtRange gen_range,
          ref_range;
  bool is_line; /* line if true, box otherwise */
} JTElem;

struct GthJumpTable {
  GtArray *elems;
};

static int compare_matches_genseq(const GthJTMatch *match_a,
                                  const GthJTMatch *match_b)
{
  gt_assert(match_a && match_b);
  return gt_range_compare(&match_a->gen_range, &match_b->gen_range);
}

static bool element_conditions_are_met(const GtArray *elems, bool debug)
{
  JTElem *elem_a, *elem_b;
  GtUword i, j;
  gt_assert(elems && gt_array_size(elems));
  if (debug)
    printf("jt: check element conditions\n");
  for (i = 0; i < gt_array_size(elems) - 1; i++) {
    for (j = i + 1; j < gt_array_size(elems); j++) {
      elem_a = gt_array_get(elems, i);
      elem_b = gt_array_get(elems, j);
      if (gt_range_overlap(&elem_a->gen_range, &elem_b->gen_range) ||
          gt_range_overlap(&elem_a->ref_range, &elem_b->ref_range) ||
          gt_range_compare(&elem_a->ref_range, &elem_b->ref_range) > 0) {
        if (debug) {
          printf("jt: not ok!\n");
          printf("jt: i="GT_WU", j="GT_WU"\n", i, j);
          printf("jt: elem_a: gen: "GT_WU", "GT_WU", ref: "GT_WU", "GT_WU"\n",
                 elem_a->gen_range.start, elem_a->gen_range.end,
                 elem_a->ref_range.start, elem_a->ref_range.end);
          printf("jt: elem_b: gen: "GT_WU", "GT_WU", ref: "GT_WU", "GT_WU"\n",
                 elem_b->gen_range.start, elem_b->gen_range.end,
                 elem_b->ref_range.start, elem_b->ref_range.end);
        }
        return false;
      }
    }
  }
  if (debug)
    printf("jt: ok\n");
  return true;
}

static void proc_matches(GthJumpTable *jt, GthJTMatch *matches,
                         GtUword num_of_matches, bool debug)
{
  GtUword i, j, cur_bits, prev_bits;
  GtBittab *bt;
  JTElem elem;
  gt_assert(jt && matches && num_of_matches);
  /* sort jump table matches according to genomic position */
  qsort(matches, num_of_matches, sizeof (GthJTMatch),
        (GtCompare) compare_matches_genseq);
  /* proc matches */
  if (debug) {
    printf("jt: matches\n");
    for (i = 0; i < num_of_matches; i++) {
      printf("jt: i=%2"GT_WUS", gen_start="GT_WU", gen_end="GT_WU", "
             "ref_start="GT_WU", ref_end="GT_WU"\n", i,
             matches[i].gen_range.start, matches[i].gen_range.end,
             matches[i].ref_range.start, matches[i].ref_range.end);
    }
  }
  bt = gt_bittab_new(num_of_matches);
  for (i = 0; i < num_of_matches; i++) {
    if (!gt_bittab_bit_is_set(bt, i)) {
      if (debug)
        printf("jt: bit "GT_WU" not set\n", i);
      elem.gen_range = matches[i].gen_range;
      elem.ref_range = matches[i].ref_range;
      if (gt_range_length(&elem.gen_range) == gt_range_length(&elem.ref_range))
        elem.is_line = true;
      else
        elem.is_line = false;
      gt_bittab_set_bit(bt, i);

      /* XXX: a more elegant algorithm probably exists */
      prev_bits = 0;
      while ((cur_bits = gt_bittab_count_set_bits(bt)) != prev_bits) {
        if (debug)
          printf("jt: looking for overlaps\n");
        prev_bits = cur_bits;
        for (j = i + 1; j < num_of_matches; j++) {
          if (!gt_bittab_bit_is_set(bt, j)) {
            if (gt_range_overlap(&elem.gen_range, &matches[j].gen_range) ||
                gt_range_overlap(&elem.ref_range, &matches[j].ref_range) ||
                gt_range_compare(&elem.ref_range, &matches[j].ref_range) > 0) {
              elem.gen_range = gt_range_join(&elem.gen_range,
                                             &matches[j].gen_range);
              elem.ref_range = gt_range_join(&elem.ref_range,
                                             &matches[j].ref_range);
              elem.is_line = false;
              gt_bittab_set_bit(bt, j);
              if (debug)
                printf("jt: bit "GT_WU" set\n", j);
            }
          }
        }
      }
      gt_array_add(jt->elems, elem);
    }
  }
  gt_assert(gt_bittab_count_set_bits(bt) == num_of_matches);
  gt_bittab_delete(bt);
  if (debug) {
    printf("jt: elements\n");
    for (i = 0; i < gt_array_size(jt->elems); i++) {
      JTElem *elem = gt_array_get(jt->elems, i);
      printf("jt: gen_start="GT_WU", gen_end="GT_WU", ref_start="GT_WU", "
             "ref_end="GT_WU", is_line=%s\n", elem->gen_range.start,
             elem->gen_range.end, elem->ref_range.start, elem->ref_range.end,
             elem->is_line ? "true" : "false");
    }
  }
  gt_assert(element_conditions_are_met(jt->elems, debug));
}

GthJumpTable* gth_jump_table_new(GthJTMatch *matches,
                                 GtUword num_of_matches, bool debug)
{
  GthJumpTable *jt;
  gt_assert(matches && num_of_matches);
#ifndef NOLICENSEMANAGER
  lm_license_check_f();
#endif
  jt = gt_malloc(sizeof *jt);
  jt->elems = gt_array_new(sizeof (JTElem));
  proc_matches(jt, matches, num_of_matches, debug);
  return jt;
}

static void fill_reverse(GthJumpTable *rjt, const GthJumpTable *fjt,
                         GtUword gen_total_length,
                         GtUword gen_offset,
                         GtUword ref_total_length,
                         GtUword ref_offset)
{
  GtUword i;
  gt_assert(rjt && fjt);
  for (i = gt_array_size(fjt->elems); i > 0; i--) {
    JTElem *forward_elem, reverse_elem;
    forward_elem = gt_array_get(fjt->elems, i-1);

    reverse_elem.gen_range.start = gen_total_length - 1
                                   - (forward_elem->gen_range.end - gen_offset)
                                   + gen_offset;

    reverse_elem.gen_range.end = gen_total_length - 1
                                 - (forward_elem->gen_range.start - gen_offset)
                                 + gen_offset;

    reverse_elem.ref_range.start = ref_total_length - 1
                                   - (forward_elem->ref_range.end - ref_offset)
                                   + ref_offset;

    reverse_elem.ref_range.end = ref_total_length - 1
                                 - (forward_elem->ref_range.start - ref_offset)
                                 + ref_offset;

    reverse_elem.is_line = forward_elem->is_line;
    gt_array_add(rjt->elems, reverse_elem);
  }
}

GthJumpTable* gth_jump_table_new_reverse(const GthJumpTable *fjt,
                                         GtUword gen_total_length,
                                         GtUword gen_offset,
                                         GtUword ref_total_length,
                                         GtUword ref_offset)

{
  GthJumpTable *rjt;
  gt_assert(fjt && gt_array_size(fjt->elems));
#ifndef NOLICENSEMANAGER
  lm_license_check_g();
#endif
  rjt = gt_malloc(sizeof *rjt);
  rjt->elems = gt_array_new(sizeof (JTElem));
  fill_reverse(rjt, fjt, gen_total_length, gen_offset, ref_total_length,
               ref_offset);
  return rjt;
}

void gth_jump_table_delete(GthJumpTable *jt)
{
  if (!jt) return;
#ifndef NOLICENSEMANAGER
  lm_license_check_h();
#endif
  gt_array_delete(jt->elems);
  gt_free(jt);
}

static GthJumpTab* make_normal_tab(const GthJumpTable *jt,
                                   const GtRange *gen_range,
                                   GtUword ref_dp_length,
                                   GtUword ref_offset,
                                   GtUword overlap,
                                   bool debug)
{
  GtUword i, el, len, gen_idx, ref_idx, offset;
  JTElem *elem = NULL, *prev_elem;
  bool was_line = true;
  GthJumpTab *tab;
  gt_assert(jt);

  gen_idx = gen_range->start;
  ref_idx = ref_offset;
  offset = gen_range->start;
  len = gt_range_length(gen_range);
  tab = gt_malloc(sizeof *tab * len);

  if (debug) {
    printf("jt: genomic elements\n");
    printf("jt: gen_range->start="GT_WU", gen_range->end="GT_WU"\n",
           gen_range->start, gen_range->end);
    for (i = 0; i < gt_array_size(jt->elems); i++) {
      JTElem *elem = gt_array_get(jt->elems, i);
      printf("jt: gen_start="GT_WU", gen_end="GT_WU", is_line=%s\n",
             elem->gen_range.start - offset, elem->gen_range.end - offset,
             elem->is_line ? "true" : "false");
    }
    printf("jt: table\n");
    printf("jt: gen_range.start="GT_WU", gen_range.end="GT_WU"\n",
           gen_range->start, gen_range->end);
    printf("jt: gen_range.start-offset="GT_WU", gen_range.end-offset="GT_WU"\n",
           gen_range->start - offset, gen_range->end - offset);
    printf("jt: offset="GT_WU"\n", offset);
    printf("jt: ref_dp_length="GT_WU"\n", ref_dp_length);
    printf("jt: len="GT_WU"\n", len);
  }

  for (el = 0; el < gt_array_size(jt->elems); el++) {
    prev_elem = elem;
    elem = gt_array_get(jt->elems, el);

    /* process part before jump table element */
    for (i = 0; gen_idx < elem->gen_range.start; i++, gen_idx++) {
      tab[gen_idx - offset].from = elem->ref_range.start + overlap - 1
                                   - ref_offset;
      if (!was_line && i + 1 < overlap)
        tab[gen_idx - offset].to = prev_elem->ref_range.start - ref_offset;
      else
        tab[gen_idx - offset].to = ref_idx - ref_offset;
    }
    if (!overlap)
      tab[gen_idx - offset - 1].to = elem->ref_range.start - ref_offset;

    /* process jump table element itself */
    gt_assert(gen_idx == elem->gen_range.start);
    for (i = 0; gen_idx + overlap <= elem->gen_range.end; i++, gen_idx++) {
      gt_assert(gen_idx >= offset);
      gt_assert(gen_idx - offset < len);
      if (elem->is_line) {
        if (gen_idx + 1 < elem->gen_range.start + overlap) {
          tab[gen_idx - offset].from = elem->ref_range.start + overlap - 1
                                       - ref_offset;
          tab[gen_idx - offset].to = ref_idx - ref_offset;
        }
        else {
          tab[gen_idx - offset].from = elem->ref_range.start + i - ref_offset;
          tab[gen_idx - offset].to = elem->ref_range.start + i + 1 - ref_offset;
        }
        was_line = true;
      }
      else {
        if (overlap) {
          if (gen_idx + 1 == elem->gen_range.start + overlap)
            ref_idx = elem->ref_range.start;
        }
        else
          ref_idx = elem->ref_range.start;
        tab[gen_idx - offset].from = elem->ref_range.end - ref_offset;
        if (gen_idx == elem->gen_range.end) {
          tab[gen_idx - offset].to = elem->ref_range.end + 1 - overlap
                                     - ref_offset;
        }
        else
          tab[gen_idx - offset].to = ref_idx - ref_offset;
        was_line = false;
      }
    }
    if (elem->ref_range.end > overlap)
      ref_idx = elem->ref_range.end - overlap + 1;
    else
      ref_idx = ref_offset;
  }

  /* process part after last jump table element */
  for (i = 0; gen_idx <= gen_range->end; i++, gen_idx++) {
    tab[gen_idx - offset].from = ref_dp_length - 1;
    if (!was_line && i + 1 < overlap)
      tab[gen_idx - offset].to = elem->ref_range.start - ref_offset;
    else
      tab[gen_idx - offset].to = ref_idx - ref_offset;
  }

  if (debug) {
    bool dots = true;
    for (i = 0; i < len; i++) {
      if (!i || tab[i-1].from != tab[i].from || tab[i-1].to != tab[i].to) {
        printf("jt: tab["GT_WU"].from="GT_WU", tab["GT_WU"].to="GT_WU"\n",
               i, tab[i].from, i, tab[i].to);
        dots = true;
      }
      else if (dots) {
        printf("jt: ...\n");
        dots = false;
      }
    }
  }

  return tab;
}

static GthJumpTab* make_spliced_tab(GthJumpTab *ntab, GtArray *ranges)
{
  GtUword i, len, offset;
  GthJumpTab *tab, *tabptr, *ntabptr;
  gt_assert(ntab && ranges);
  len = gt_ranges_total_length(ranges);
  tab = gt_malloc(sizeof *tab * len);
  tabptr = tab;
  offset = ((GtRange*) gt_array_get_first(ranges))->start;
  for (i = 0; i < gt_array_size(ranges); i++) {
    for (ntabptr = ntab + ((GtRange*) gt_array_get(ranges, i))->start - offset;
         ntabptr <= ntab + ((GtRange*) gt_array_get(ranges, i))->end - offset;
         *tabptr++ = *ntabptr++);

  }
  gt_assert(tabptr == tab + len);
  return tab;
}

GthJumpTab* gth_jump_table_make_tab(const GthJumpTable *jt, GtArray *gen_ranges,
                                    GtUword ref_dp_length,
                                    GtUword ref_offset,
                                    GtUword overlap, bool debug)
{
  GthJumpTab *normal_tab, *spliced_tab;
  GtRange gen_range;
  gt_assert(jt && gen_ranges && gt_array_size(gen_ranges));

  gen_range.start = ((GtRange*) gt_array_get_first(gen_ranges))->start;
  gen_range.end = ((GtRange*) gt_array_get_last(gen_ranges))->end;
  normal_tab = make_normal_tab(jt, &gen_range, ref_dp_length, ref_offset,
                               overlap, debug);
  if (gt_array_size(gen_ranges) == 1)
    return normal_tab;
  else {
    spliced_tab = make_spliced_tab(normal_tab, gen_ranges);
    gt_free(normal_tab);
    return spliced_tab;
  }
}

GtRowInfo* gth_jump_tab_get_full_row_info(const GthJumpTab *jt,
                                          GtUword gen_dp_length)
{
  GtUword i;
  GtRowInfo *ri;
  gt_assert(jt);

  ri = gt_malloc((gen_dp_length + 1) * sizeof *ri);

  ri[0].offset = 0;
  ri[0].length = jt[0].from + 1 + 1;

  ri[1].offset = 0;
  ri[1].length = jt[0].from + 1 + 1;

  for (i = 1; i < gen_dp_length; i++) {
    ri[i+1].offset = jt[i-1].to;
    ri[i+1].length = jt[i].from - jt[i-1].to + 1;
    if (ri[i+1].offset == 0)
      ri[i+1].length++;
    else
      ri[i+1].offset++;
  }

  return ri;
}

static GtRowInfo* get_half_row_info(GtRowInfo *ri_full,
                                    GtUword gen_dp_length,
                                    GtUword *size)
{
  GtUword i, x, y, full_length, half_length;
  GtRowInfo *ri_half;
  gt_assert(ri_full && size);

  full_length = gen_dp_length + 1;
  half_length = GT_DIV2(full_length) + GT_MOD2(full_length);
  ri_half = gt_malloc(half_length * sizeof *ri_half);
  *size = 0;

  for (i = 0; i < half_length; i++) {
    x = 2 * i;
    y = x + 1;
    ri_half[i].offset = ri_full[x].offset;
    if (y < full_length) {
      if (ri_full[x].offset == ri_full[y].offset)
        ri_half[i].length = GT_MAX(ri_full[x].length, ri_full[y].length);
      else {
        gt_assert(ri_full[x].offset < ri_full[y].offset);
        if (ri_full[x].offset + ri_full[x].length > ri_full[y].offset) {
          /* length(x) + length(y) - overlap(x,y) */
          ri_half[i].length = ri_full[x].length + ri_full[y].length -
                              (ri_full[x].offset + ri_full[x].length -
                               ri_full[y].offset);
        }
        else
          ri_half[i].length = ri_full[x].length + ri_full[y].length;
      }
    }
    else
      ri_half[i].length = ri_full[x].length;
    *size += ri_half[i].length;
  }

  return ri_half;
}

GtRowInfo* gth_jump_tab_get_row_info(const GthJumpTab *jt,
                                     GtUword gen_dp_length,
                                     GtUword *size)
{
  GtRowInfo *ri_full, *ri_half;
  gt_assert(jt && size);
  ri_full = gth_jump_tab_get_full_row_info(jt, gen_dp_length);
  ri_half = get_half_row_info(ri_full, gen_dp_length, size);
  gt_free(ri_full);
  return ri_half;
}

static GthJumpTab *simple_balanced_jump_tab(GtUword overlap)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch match;
  GthJumpTab *jt;

  match.gen_range.start = 3;
  match.gen_range.end = 6;
  match.ref_range.start = 3;
  match.ref_range.end = 6;

  jump_table = gth_jump_table_new(&match, 1, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 9;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 10, 0, overlap, false);

  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return jt;
}

static int simple_balanced_no_overlap(GtError *err)
{
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  jt = simple_balanced_jump_tab(0);

  gt_ensure(jt[0].from == 2 && jt[0].to == 0);
  gt_ensure(jt[1].from == 2 && jt[1].to == 0);
  gt_ensure(jt[2].from == 2 && jt[2].to == 3);
  gt_ensure(jt[3].from == 3 && jt[3].to == 4);
  gt_ensure(jt[4].from == 4 && jt[4].to == 5);
  gt_ensure(jt[5].from == 5 && jt[5].to == 6);
  gt_ensure(jt[6].from == 6 && jt[6].to == 7);
  gt_ensure(jt[7].from == 9 && jt[7].to == 7);
  gt_ensure(jt[8].from == 9 && jt[8].to == 7);
  gt_ensure(jt[9].from == 9 && jt[9].to == 7);

  gt_free(jt);

  return had_err;
}

static int simple_unbalanced_no_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch match;
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  match.gen_range.start = 3;
  match.gen_range.end = 7;
  match.ref_range.start = 3;
  match.ref_range.end = 6;

  jump_table = gth_jump_table_new(&match, 1, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 9;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 10, 0, 0, false);

  gt_ensure(jt[0].from == 2 && jt[0].to == 0);
  gt_ensure(jt[1].from == 2 && jt[1].to == 0);
  gt_ensure(jt[2].from == 2 && jt[2].to == 3);
  gt_ensure(jt[3].from == 6 && jt[3].to == 3);
  gt_ensure(jt[4].from == 6 && jt[4].to == 3);
  gt_ensure(jt[5].from == 6 && jt[5].to == 3);
  gt_ensure(jt[6].from == 6 && jt[6].to == 3);
  gt_ensure(jt[7].from == 6 && jt[7].to == 7);
  gt_ensure(jt[8].from == 9 && jt[8].to == 7);
  gt_ensure(jt[9].from == 9 && jt[9].to == 7);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int simple_cluster_no_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch matches[2];
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  matches[0].gen_range.start = 3;
  matches[0].gen_range.end = 5;
  matches[0].ref_range.start = 3;
  matches[0].ref_range.end = 5;
  matches[1].gen_range.start = 5;
  matches[1].gen_range.end = 7;
  matches[1].ref_range.start = 4;
  matches[1].ref_range.end = 6;

  jump_table = gth_jump_table_new(matches, 2, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 9;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 10, 0, 0, false);

  gt_ensure(jt[0].from == 2 && jt[0].to == 0);
  gt_ensure(jt[1].from == 2 && jt[1].to == 0);
  gt_ensure(jt[2].from == 2 && jt[2].to == 3);
  gt_ensure(jt[3].from == 6 && jt[3].to == 3);
  gt_ensure(jt[4].from == 6 && jt[4].to == 3);
  gt_ensure(jt[5].from == 6 && jt[5].to == 3);
  gt_ensure(jt[6].from == 6 && jt[6].to == 3);
  gt_ensure(jt[7].from == 6 && jt[7].to == 7);
  gt_ensure(jt[8].from == 9 && jt[8].to == 7);
  gt_ensure(jt[9].from == 9 && jt[9].to == 7);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int simple_balanced_small_overlap(GtError *err)
{
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  jt = simple_balanced_jump_tab(1);

  gt_ensure(jt[0].from == 3 && jt[0].to == 0);
  gt_ensure(jt[1].from == 3 && jt[1].to == 0);
  gt_ensure(jt[2].from == 3 && jt[2].to == 0);
  gt_ensure(jt[3].from == 3 && jt[3].to == 4);
  gt_ensure(jt[4].from == 4 && jt[4].to == 5);
  gt_ensure(jt[5].from == 5 && jt[5].to == 6);
  gt_ensure(jt[6].from == 9 && jt[6].to == 6);
  gt_ensure(jt[7].from == 9 && jt[7].to == 6);
  gt_ensure(jt[8].from == 9 && jt[8].to == 6);
  gt_ensure(jt[9].from == 9 && jt[9].to == 6);

  gt_free(jt);

  return had_err;
}

static int simple_unbalanced_small_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch match;
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  match.gen_range.start = 3;
  match.gen_range.end = 7;
  match.ref_range.start = 3;
  match.ref_range.end = 6;

  jump_table = gth_jump_table_new(&match, 1, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 9;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 10, 0, 1, false);

  gt_ensure(jt[0].from == 3 && jt[0].to == 0);
  gt_ensure(jt[1].from == 3 && jt[1].to == 0);
  gt_ensure(jt[2].from == 3 && jt[2].to == 0);
  gt_ensure(jt[3].from == 6 && jt[3].to == 3);
  gt_ensure(jt[4].from == 6 && jt[4].to == 3);
  gt_ensure(jt[5].from == 6 && jt[5].to == 3);
  gt_ensure(jt[6].from == 6 && jt[6].to == 3);
  gt_ensure(jt[7].from == 9 && jt[7].to == 6);
  gt_ensure(jt[8].from == 9 && jt[8].to == 6);
  gt_ensure(jt[9].from == 9 && jt[9].to == 6);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int simple_cluster_small_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch matches[2];
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  matches[0].gen_range.start = 3;
  matches[0].gen_range.end = 5;
  matches[0].ref_range.start = 3;
  matches[0].ref_range.end = 5;
  matches[1].gen_range.start = 5;
  matches[1].gen_range.end = 7;
  matches[1].ref_range.start = 4;
  matches[1].ref_range.end = 6;

  jump_table = gth_jump_table_new(matches, 2, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 9;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 10, 0, 1, false);

  gt_ensure(jt[0].from == 3 && jt[0].to == 0);
  gt_ensure(jt[1].from == 3 && jt[1].to == 0);
  gt_ensure(jt[2].from == 3 && jt[2].to == 0);
  gt_ensure(jt[3].from == 6 && jt[3].to == 3);
  gt_ensure(jt[4].from == 6 && jt[4].to == 3);
  gt_ensure(jt[5].from == 6 && jt[5].to == 3);
  gt_ensure(jt[6].from == 6 && jt[6].to == 3);
  gt_ensure(jt[7].from == 9 && jt[7].to == 6);
  gt_ensure(jt[8].from == 9 && jt[8].to == 6);
  gt_ensure(jt[9].from == 9 && jt[9].to == 6);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int simple_balanced_large_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch match;
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  match.gen_range.start = 5;
  match.gen_range.end = 19;
  match.ref_range.start = 5;
  match.ref_range.end = 19;

  jump_table = gth_jump_table_new(&match, 1, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 24;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 25, 0, 5, false);

  gt_ensure(jt[0].from == 9 && jt[0].to == 0);
  gt_ensure(jt[1].from == 9 && jt[1].to == 0);
  gt_ensure(jt[2].from == 9 && jt[2].to == 0);
  gt_ensure(jt[3].from == 9 && jt[3].to == 0);
  gt_ensure(jt[4].from == 9 && jt[4].to == 0);
  gt_ensure(jt[5].from == 9 && jt[5].to == 0);
  gt_ensure(jt[6].from == 9 && jt[6].to == 0);
  gt_ensure(jt[7].from == 9 && jt[7].to == 0);
  gt_ensure(jt[8].from == 9 && jt[8].to == 0);
  gt_ensure(jt[9].from == 9 && jt[9].to == 10);
  gt_ensure(jt[10].from == 10 && jt[10].to == 11);
  gt_ensure(jt[11].from == 11 && jt[11].to == 12);
  gt_ensure(jt[12].from == 12 && jt[12].to == 13);
  gt_ensure(jt[13].from == 13 && jt[13].to == 14);
  gt_ensure(jt[14].from == 14 && jt[14].to == 15);
  gt_ensure(jt[15].from == 24 && jt[15].to == 15);
  gt_ensure(jt[16].from == 24 && jt[16].to == 15);
  gt_ensure(jt[17].from == 24 && jt[17].to == 15);
  gt_ensure(jt[18].from == 24 && jt[18].to == 15);
  gt_ensure(jt[19].from == 24 && jt[19].to == 15);
  gt_ensure(jt[20].from == 24 && jt[20].to == 15);
  gt_ensure(jt[21].from == 24 && jt[21].to == 15);
  gt_ensure(jt[22].from == 24 && jt[22].to == 15);
  gt_ensure(jt[23].from == 24 && jt[23].to == 15);
  gt_ensure(jt[24].from == 24 && jt[24].to == 15);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int simple_unbalanced_large_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch match;
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  match.gen_range.start = 5;
  match.gen_range.end = 20;
  match.ref_range.start = 5;
  match.ref_range.end = 19;

  jump_table = gth_jump_table_new(&match, 1, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 24;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 25, 0, 5, false);

  gt_ensure(jt[0].from == 9 && jt[0].to == 0);
  gt_ensure(jt[1].from == 9 && jt[1].to == 0);
  gt_ensure(jt[2].from == 9 && jt[2].to == 0);
  gt_ensure(jt[3].from == 9 && jt[3].to == 0);
  gt_ensure(jt[4].from == 9 && jt[4].to == 0);
  gt_ensure(jt[5].from == 19 && jt[5].to == 0);
  gt_ensure(jt[6].from == 19 && jt[6].to == 0);
  gt_ensure(jt[7].from == 19 && jt[7].to == 0);
  gt_ensure(jt[8].from == 19 && jt[8].to == 0);
  gt_ensure(jt[9].from == 19 && jt[9].to == 5);
  gt_ensure(jt[10].from == 19 && jt[10].to == 5);
  gt_ensure(jt[11].from == 19 && jt[11].to == 5);
  gt_ensure(jt[12].from == 19 && jt[12].to == 5);
  gt_ensure(jt[13].from == 19 && jt[13].to == 5);
  gt_ensure(jt[14].from == 19 && jt[14].to == 5);
  gt_ensure(jt[15].from == 19 && jt[15].to == 5);
  gt_ensure(jt[16].from == 24 && jt[16].to == 5);
  gt_ensure(jt[17].from == 24 && jt[17].to == 5);
  gt_ensure(jt[18].from == 24 && jt[18].to == 5);
  gt_ensure(jt[19].from == 24 && jt[19].to == 5);
  gt_ensure(jt[20].from == 24 && jt[20].to == 15);
  gt_ensure(jt[21].from == 24 && jt[21].to == 15);
  gt_ensure(jt[22].from == 24 && jt[22].to == 15);
  gt_ensure(jt[23].from == 24 && jt[23].to == 15);
  gt_ensure(jt[24].from == 24 && jt[24].to == 15);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int simple_cluster_large_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch matches[2];
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  matches[0].gen_range.start = 5;
  matches[0].gen_range.end = 14;
  matches[0].ref_range.start = 5;
  matches[0].ref_range.end = 14;
  matches[1].gen_range.start = 11;
  matches[1].gen_range.end = 20;
  matches[1].ref_range.start = 10;
  matches[1].ref_range.end = 19;

  jump_table = gth_jump_table_new(matches, 2, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 24;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 25, 0, 5, false);

  gt_ensure(jt[0].from == 9 && jt[0].to == 0);
  gt_ensure(jt[1].from == 9 && jt[1].to == 0);
  gt_ensure(jt[2].from == 9 && jt[2].to == 0);
  gt_ensure(jt[3].from == 9 && jt[3].to == 0);
  gt_ensure(jt[4].from == 9 && jt[4].to == 0);
  gt_ensure(jt[5].from == 19 && jt[5].to == 0);
  gt_ensure(jt[6].from == 19 && jt[6].to == 0);
  gt_ensure(jt[7].from == 19 && jt[7].to == 0);
  gt_ensure(jt[8].from == 19 && jt[8].to == 0);
  gt_ensure(jt[9].from == 19 && jt[9].to == 5);
  gt_ensure(jt[10].from == 19 && jt[10].to == 5);
  gt_ensure(jt[11].from == 19 && jt[11].to == 5);
  gt_ensure(jt[12].from == 19 && jt[12].to == 5);
  gt_ensure(jt[13].from == 19 && jt[13].to == 5);
  gt_ensure(jt[14].from == 19 && jt[14].to == 5);
  gt_ensure(jt[15].from == 19 && jt[15].to == 5);
  gt_ensure(jt[16].from == 24 && jt[16].to == 5);
  gt_ensure(jt[17].from == 24 && jt[17].to == 5);
  gt_ensure(jt[18].from == 24 && jt[18].to == 5);
  gt_ensure(jt[19].from == 24 && jt[19].to == 5);
  gt_ensure(jt[20].from == 24 && jt[20].to == 15);
  gt_ensure(jt[21].from == 24 && jt[21].to == 15);
  gt_ensure(jt[22].from == 24 && jt[22].to == 15);
  gt_ensure(jt[23].from == 24 && jt[23].to == 15);
  gt_ensure(jt[24].from == 24 && jt[24].to == 15);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int double_balanced_no_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch matches[2];
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  matches[0].gen_range.start = 5;
  matches[0].gen_range.end = 14;
  matches[0].ref_range.start = 5;
  matches[0].ref_range.end = 14;
  matches[1].gen_range.start = 20;
  matches[1].gen_range.end = 29;
  matches[1].ref_range.start = 20;
  matches[1].ref_range.end = 29;

  jump_table = gth_jump_table_new(matches, 2, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 34;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 35, 0, 0, false);

  gt_ensure(jt[0].from == 4 && jt[0].to == 0);
  gt_ensure(jt[1].from == 4 && jt[1].to == 0);
  gt_ensure(jt[2].from == 4 && jt[2].to == 0);
  gt_ensure(jt[3].from == 4 && jt[3].to == 0);
  gt_ensure(jt[4].from == 4 && jt[4].to == 5);
  gt_ensure(jt[5].from == 5 && jt[5].to == 6);
  gt_ensure(jt[6].from == 6 && jt[6].to == 7);
  gt_ensure(jt[7].from == 7 && jt[7].to == 8);
  gt_ensure(jt[8].from == 8 && jt[8].to == 9);
  gt_ensure(jt[9].from == 9 && jt[9].to == 10);
  gt_ensure(jt[10].from == 10 && jt[10].to == 11);
  gt_ensure(jt[11].from == 11 && jt[11].to == 12);
  gt_ensure(jt[12].from == 12 && jt[12].to == 13);
  gt_ensure(jt[13].from == 13 && jt[13].to == 14);
  gt_ensure(jt[14].from == 14 && jt[14].to == 15);
  gt_ensure(jt[15].from == 19 && jt[15].to == 15);
  gt_ensure(jt[16].from == 19 && jt[16].to == 15);
  gt_ensure(jt[17].from == 19 && jt[17].to == 15);
  gt_ensure(jt[18].from == 19 && jt[18].to == 15);
  gt_ensure(jt[19].from == 19 && jt[19].to == 20);
  gt_ensure(jt[20].from == 20 && jt[20].to == 21);
  gt_ensure(jt[21].from == 21 && jt[21].to == 22);
  gt_ensure(jt[22].from == 22 && jt[22].to == 23);
  gt_ensure(jt[23].from == 23 && jt[23].to == 24);
  gt_ensure(jt[24].from == 24 && jt[24].to == 25);
  gt_ensure(jt[25].from == 25 && jt[25].to == 26);
  gt_ensure(jt[26].from == 26 && jt[26].to == 27);
  gt_ensure(jt[27].from == 27 && jt[27].to == 28);
  gt_ensure(jt[28].from == 28 && jt[28].to == 29);
  gt_ensure(jt[29].from == 29 && jt[29].to == 30);
  gt_ensure(jt[30].from == 34 && jt[30].to == 30);
  gt_ensure(jt[31].from == 34 && jt[31].to == 30);
  gt_ensure(jt[32].from == 34 && jt[32].to == 30);
  gt_ensure(jt[33].from == 34 && jt[33].to == 30);
  gt_ensure(jt[34].from == 34 && jt[34].to == 30);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int double_cluster_no_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch matches[4];
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  matches[0].gen_range.start = 5;
  matches[0].gen_range.end = 11;
  matches[0].ref_range.start = 5;
  matches[0].ref_range.end = 11;
  matches[1].gen_range.start = 9;
  matches[1].gen_range.end = 15;
  matches[1].ref_range.start = 8;
  matches[1].ref_range.end = 14;
  matches[2].gen_range.start = 20;
  matches[2].gen_range.end = 27;
  matches[2].ref_range.start = 20;
  matches[2].ref_range.end = 27;
  matches[3].gen_range.start = 24;
  matches[3].gen_range.end = 30;
  matches[3].ref_range.start = 23;
  matches[3].ref_range.end = 29;

  jump_table = gth_jump_table_new(matches, 4, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 34;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 35, 0, 0, false);

  gt_ensure(jt[0].from == 4 && jt[0].to == 0);
  gt_ensure(jt[1].from == 4 && jt[1].to == 0);
  gt_ensure(jt[2].from == 4 && jt[2].to == 0);
  gt_ensure(jt[3].from == 4 && jt[3].to == 0);
  gt_ensure(jt[4].from == 4 && jt[4].to == 5);
  gt_ensure(jt[5].from == 14 && jt[5].to == 5);
  gt_ensure(jt[6].from == 14 && jt[6].to == 5);
  gt_ensure(jt[7].from == 14 && jt[7].to == 5);
  gt_ensure(jt[8].from == 14 && jt[8].to == 5);
  gt_ensure(jt[9].from == 14 && jt[9].to == 5);
  gt_ensure(jt[10].from == 14 && jt[10].to == 5);
  gt_ensure(jt[11].from == 14 && jt[11].to == 5);
  gt_ensure(jt[12].from == 14 && jt[12].to == 5);
  gt_ensure(jt[13].from == 14 && jt[13].to == 5);
  gt_ensure(jt[14].from == 14 && jt[14].to == 5);
  gt_ensure(jt[15].from == 14 && jt[15].to == 15);
  gt_ensure(jt[16].from == 19 && jt[16].to == 15);
  gt_ensure(jt[17].from == 19 && jt[17].to == 15);
  gt_ensure(jt[18].from == 19 && jt[18].to == 15);
  gt_ensure(jt[19].from == 19 && jt[19].to == 20);
  gt_ensure(jt[20].from == 29 && jt[20].to == 20);
  gt_ensure(jt[21].from == 29 && jt[21].to == 20);
  gt_ensure(jt[22].from == 29 && jt[22].to == 20);
  gt_ensure(jt[23].from == 29 && jt[23].to == 20);
  gt_ensure(jt[24].from == 29 && jt[24].to == 20);
  gt_ensure(jt[25].from == 29 && jt[25].to == 20);
  gt_ensure(jt[26].from == 29 && jt[26].to == 20);
  gt_ensure(jt[27].from == 29 && jt[27].to == 20);
  gt_ensure(jt[28].from == 29 && jt[28].to == 20);
  gt_ensure(jt[29].from == 29 && jt[29].to == 20);
  gt_ensure(jt[30].from == 29 && jt[30].to == 30);
  gt_ensure(jt[31].from == 34 && jt[31].to == 30);
  gt_ensure(jt[32].from == 34 && jt[32].to == 30);
  gt_ensure(jt[33].from == 34 && jt[33].to == 30);
  gt_ensure(jt[34].from == 34 && jt[34].to == 30);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int double_unbalanced_no_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch matches[2];
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  matches[0].gen_range.start = 5;
  matches[0].gen_range.end = 15;
  matches[0].ref_range.start = 5;
  matches[0].ref_range.end = 14;
  matches[1].gen_range.start = 20;
  matches[1].gen_range.end = 30;
  matches[1].ref_range.start = 20;
  matches[1].ref_range.end = 29;

  jump_table = gth_jump_table_new(matches, 2, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 34;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 35, 0, 0, false);

  gt_ensure(jt[0].from == 4 && jt[0].to == 0);
  gt_ensure(jt[1].from == 4 && jt[1].to == 0);
  gt_ensure(jt[2].from == 4 && jt[2].to == 0);
  gt_ensure(jt[3].from == 4 && jt[3].to == 0);
  gt_ensure(jt[4].from == 4 && jt[4].to == 5);
  gt_ensure(jt[5].from == 14 && jt[5].to == 5);
  gt_ensure(jt[6].from == 14 && jt[6].to == 5);
  gt_ensure(jt[7].from == 14 && jt[7].to == 5);
  gt_ensure(jt[8].from == 14 && jt[8].to == 5);
  gt_ensure(jt[9].from == 14 && jt[9].to == 5);
  gt_ensure(jt[10].from == 14 && jt[10].to == 5);
  gt_ensure(jt[11].from == 14 && jt[11].to == 5);
  gt_ensure(jt[12].from == 14 && jt[12].to == 5);
  gt_ensure(jt[13].from == 14 && jt[13].to == 5);
  gt_ensure(jt[14].from == 14 && jt[14].to == 5);
  gt_ensure(jt[15].from == 14 && jt[15].to == 15);
  gt_ensure(jt[16].from == 19 && jt[16].to == 15);
  gt_ensure(jt[17].from == 19 && jt[17].to == 15);
  gt_ensure(jt[18].from == 19 && jt[18].to == 15);
  gt_ensure(jt[19].from == 19 && jt[19].to == 20);
  gt_ensure(jt[20].from == 29 && jt[20].to == 20);
  gt_ensure(jt[21].from == 29 && jt[21].to == 20);
  gt_ensure(jt[22].from == 29 && jt[22].to == 20);
  gt_ensure(jt[23].from == 29 && jt[23].to == 20);
  gt_ensure(jt[24].from == 29 && jt[24].to == 20);
  gt_ensure(jt[25].from == 29 && jt[25].to == 20);
  gt_ensure(jt[26].from == 29 && jt[26].to == 20);
  gt_ensure(jt[27].from == 29 && jt[27].to == 20);
  gt_ensure(jt[28].from == 29 && jt[28].to == 20);
  gt_ensure(jt[29].from == 29 && jt[29].to == 20);
  gt_ensure(jt[30].from == 29 && jt[30].to == 30);
  gt_ensure(jt[31].from == 34 && jt[31].to == 30);
  gt_ensure(jt[32].from == 34 && jt[32].to == 30);
  gt_ensure(jt[33].from == 34 && jt[33].to == 30);
  gt_ensure(jt[34].from == 34 && jt[34].to == 30);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int double_balanced_with_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch matches[2];
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  matches[0].gen_range.start = 5;
  matches[0].gen_range.end = 14;
  matches[0].ref_range.start = 5;
  matches[0].ref_range.end = 14;
  matches[1].gen_range.start = 20;
  matches[1].gen_range.end = 29;
  matches[1].ref_range.start = 20;
  matches[1].ref_range.end = 29;

  jump_table = gth_jump_table_new(matches, 2, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 34;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 35, 0, 3, false);

  gt_ensure(jt[0].from == 7 && jt[0].to == 0);
  gt_ensure(jt[1].from == 7 && jt[1].to == 0);
  gt_ensure(jt[2].from == 7 && jt[2].to == 0);
  gt_ensure(jt[3].from == 7 && jt[3].to == 0);
  gt_ensure(jt[4].from == 7 && jt[4].to == 0);
  gt_ensure(jt[5].from == 7 && jt[5].to == 0);
  gt_ensure(jt[6].from == 7 && jt[6].to == 0);
  gt_ensure(jt[7].from == 7 && jt[7].to == 8);
  gt_ensure(jt[8].from == 8 && jt[8].to == 9);
  gt_ensure(jt[9].from == 9 && jt[9].to == 10);
  gt_ensure(jt[10].from == 10 && jt[10].to == 11);
  gt_ensure(jt[11].from == 11 && jt[11].to == 12);
  gt_ensure(jt[12].from == 22 && jt[12].to == 12);
  gt_ensure(jt[13].from == 22 && jt[13].to == 12);
  gt_ensure(jt[14].from == 22 && jt[14].to == 12);
  gt_ensure(jt[15].from == 22 && jt[15].to == 12);
  gt_ensure(jt[16].from == 22 && jt[16].to == 12);
  gt_ensure(jt[17].from == 22 && jt[17].to == 12);
  gt_ensure(jt[18].from == 22 && jt[18].to == 12);
  gt_ensure(jt[19].from == 22 && jt[19].to == 12);
  gt_ensure(jt[20].from == 22 && jt[20].to == 12);
  gt_ensure(jt[21].from == 22 && jt[21].to == 12);
  gt_ensure(jt[22].from == 22 && jt[22].to == 23);
  gt_ensure(jt[23].from == 23 && jt[23].to == 24);
  gt_ensure(jt[24].from == 24 && jt[24].to == 25);
  gt_ensure(jt[25].from == 25 && jt[25].to == 26);
  gt_ensure(jt[26].from == 26 && jt[26].to == 27);
  gt_ensure(jt[27].from == 34 && jt[27].to == 27);
  gt_ensure(jt[28].from == 34 && jt[28].to == 27);
  gt_ensure(jt[29].from == 34 && jt[29].to == 27);
  gt_ensure(jt[30].from == 34 && jt[30].to == 27);
  gt_ensure(jt[31].from == 34 && jt[31].to == 27);
  gt_ensure(jt[32].from == 34 && jt[32].to == 27);
  gt_ensure(jt[33].from == 34 && jt[33].to == 27);
  gt_ensure(jt[34].from == 34 && jt[34].to == 27);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int double_unbalanced_with_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch matches[2];
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  matches[0].gen_range.start = 5;
  matches[0].gen_range.end = 15;
  matches[0].ref_range.start = 5;
  matches[0].ref_range.end = 14;
  matches[1].gen_range.start = 20;
  matches[1].gen_range.end = 30;
  matches[1].ref_range.start = 20;
  matches[1].ref_range.end = 29;

  jump_table = gth_jump_table_new(matches, 2, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 34;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 35, 0, 3, false);

  gt_ensure(jt[0].from == 7 && jt[0].to == 0);
  gt_ensure(jt[1].from == 7 && jt[1].to == 0);
  gt_ensure(jt[2].from == 7 && jt[2].to == 0);
  gt_ensure(jt[3].from == 7 && jt[3].to == 0);
  gt_ensure(jt[4].from == 7 && jt[4].to == 0);
  gt_ensure(jt[5].from == 14 && jt[5].to == 0);
  gt_ensure(jt[6].from == 14 && jt[6].to == 0);
  gt_ensure(jt[7].from == 14 && jt[7].to == 5);
  gt_ensure(jt[8].from == 14 && jt[8].to == 5);
  gt_ensure(jt[9].from == 14 && jt[9].to == 5);
  gt_ensure(jt[10].from == 14 && jt[10].to == 5);
  gt_ensure(jt[11].from == 14 && jt[11].to == 5);
  gt_ensure(jt[12].from == 14 && jt[12].to == 5);
  gt_ensure(jt[13].from == 22 && jt[13].to == 5);
  gt_ensure(jt[14].from == 22 && jt[14].to == 5);
  gt_ensure(jt[15].from == 22 && jt[15].to == 12);
  gt_ensure(jt[16].from == 22 && jt[16].to == 12);
  gt_ensure(jt[17].from == 22 && jt[17].to == 12);
  gt_ensure(jt[18].from == 22 && jt[18].to == 12);
  gt_ensure(jt[19].from == 22 && jt[19].to == 12);
  gt_ensure(jt[20].from == 29 && jt[20].to == 12);
  gt_ensure(jt[21].from == 29 && jt[21].to == 12);
  gt_ensure(jt[22].from == 29 && jt[22].to == 20);
  gt_ensure(jt[23].from == 29 && jt[23].to == 20);
  gt_ensure(jt[24].from == 29 && jt[24].to == 20);
  gt_ensure(jt[25].from == 29 && jt[25].to == 20);
  gt_ensure(jt[26].from == 29 && jt[26].to == 20);
  gt_ensure(jt[27].from == 29 && jt[27].to == 20);
  gt_ensure(jt[28].from == 34 && jt[28].to == 20);
  gt_ensure(jt[29].from == 34 && jt[29].to == 20);
  gt_ensure(jt[30].from == 34 && jt[30].to == 27);
  gt_ensure(jt[31].from == 34 && jt[31].to == 27);
  gt_ensure(jt[32].from == 34 && jt[32].to == 27);
  gt_ensure(jt[33].from == 34 && jt[33].to == 27);
  gt_ensure(jt[34].from == 34 && jt[34].to == 27);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int double_cluster_with_overlap(GtError *err)
{
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtRange gen_range;
  GthJTMatch matches[4];
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  matches[0].gen_range.start = 5;
  matches[0].gen_range.end = 11;
  matches[0].ref_range.start = 5;
  matches[0].ref_range.end = 11;
  matches[1].gen_range.start = 9;
  matches[1].gen_range.end = 15;
  matches[1].ref_range.start = 8;
  matches[1].ref_range.end = 14;
  matches[2].gen_range.start = 20;
  matches[2].gen_range.end = 27;
  matches[2].ref_range.start = 20;
  matches[2].ref_range.end = 27;
  matches[3].gen_range.start = 24;
  matches[3].gen_range.end = 30;
  matches[3].ref_range.start = 23;
  matches[3].ref_range.end = 29;

  jump_table = gth_jump_table_new(matches, 4, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 34;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 35, 0, 3, false);

  gt_ensure(jt[0].from == 7 && jt[0].to == 0);
  gt_ensure(jt[1].from == 7 && jt[1].to == 0);
  gt_ensure(jt[2].from == 7 && jt[2].to == 0);
  gt_ensure(jt[3].from == 7 && jt[3].to == 0);
  gt_ensure(jt[4].from == 7 && jt[4].to == 0);
  gt_ensure(jt[5].from == 14 && jt[5].to == 0);
  gt_ensure(jt[6].from == 14 && jt[6].to == 0);
  gt_ensure(jt[7].from == 14 && jt[7].to == 5);
  gt_ensure(jt[8].from == 14 && jt[8].to == 5);
  gt_ensure(jt[9].from == 14 && jt[9].to == 5);
  gt_ensure(jt[10].from == 14 && jt[10].to == 5);
  gt_ensure(jt[11].from == 14 && jt[11].to == 5);
  gt_ensure(jt[12].from == 14 && jt[12].to == 5);
  gt_ensure(jt[13].from == 22 && jt[13].to == 5);
  gt_ensure(jt[14].from == 22 && jt[14].to == 5);
  gt_ensure(jt[15].from == 22 && jt[15].to == 12);
  gt_ensure(jt[16].from == 22 && jt[16].to == 12);
  gt_ensure(jt[17].from == 22 && jt[17].to == 12);
  gt_ensure(jt[18].from == 22 && jt[18].to == 12);
  gt_ensure(jt[19].from == 22 && jt[19].to == 12);
  gt_ensure(jt[20].from == 29 && jt[20].to == 12);
  gt_ensure(jt[21].from == 29 && jt[21].to == 12);
  gt_ensure(jt[22].from == 29 && jt[22].to == 20);
  gt_ensure(jt[23].from == 29 && jt[23].to == 20);
  gt_ensure(jt[24].from == 29 && jt[24].to == 20);
  gt_ensure(jt[25].from == 29 && jt[25].to == 20);
  gt_ensure(jt[26].from == 29 && jt[26].to == 20);
  gt_ensure(jt[27].from == 29 && jt[27].to == 20);
  gt_ensure(jt[28].from == 34 && jt[28].to == 20);
  gt_ensure(jt[29].from == 34 && jt[29].to == 20);
  gt_ensure(jt[30].from == 34 && jt[30].to == 27);
  gt_ensure(jt[31].from == 34 && jt[31].to == 27);
  gt_ensure(jt[32].from == 34 && jt[32].to == 27);
  gt_ensure(jt[33].from == 34 && jt[33].to == 27);
  gt_ensure(jt[34].from == 34 && jt[34].to == 27);

  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

static int row_info_simple_balanced_no_overlap(GtError *err)
{
  GtRowInfo *full_ri, *half_ri;
  GtUword size;
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  jt = simple_balanced_jump_tab(0);
  full_ri = gth_jump_tab_get_full_row_info(jt, 10);

  gt_ensure(full_ri[0].offset == 0 && full_ri[0].length == 4);
  gt_ensure(full_ri[1].offset == 0 && full_ri[1].length == 4);
  gt_ensure(full_ri[2].offset == 0 && full_ri[2].length == 4);
  gt_ensure(full_ri[3].offset == 0 && full_ri[3].length == 4);
  gt_ensure(full_ri[4].offset == 4 && full_ri[4].length == 1);
  gt_ensure(full_ri[5].offset == 5 && full_ri[5].length == 1);
  gt_ensure(full_ri[6].offset == 6 && full_ri[6].length == 1);
  gt_ensure(full_ri[7].offset == 7 && full_ri[7].length == 1);
  gt_ensure(full_ri[8].offset == 8 && full_ri[8].length == 3);
  gt_ensure(full_ri[9].offset == 8 && full_ri[9].length == 3);
  gt_ensure(full_ri[10].offset == 8 && full_ri[10].length == 3);

  half_ri = gth_jump_tab_get_row_info(jt, 10, &size);

  gt_ensure(half_ri[0].offset == 0 && half_ri[0].length == 4);
  gt_ensure(half_ri[1].offset == 0 && half_ri[1].length == 4);
  gt_ensure(half_ri[2].offset == 4 && half_ri[2].length == 2);
  gt_ensure(half_ri[3].offset == 6 && half_ri[3].length == 2);
  gt_ensure(half_ri[4].offset == 8 && half_ri[4].length == 3);
  gt_ensure(half_ri[5].offset == 8 && half_ri[5].length == 3);
  gt_ensure(size == 18);

  gt_free(half_ri);
  gt_free(full_ri);
  gt_free(jt);

  return had_err;
}

static int row_info_simple_unbalanced_with_overlap(GtError *err)
{
  GtRowInfo *full_ri, *half_ri;
  GthJumpTable *jump_table;
  GtArray *gen_ranges;
  GtUword size;
  GtRange gen_range;
  GthJTMatch match;
  GthJumpTab *jt;
  int had_err = 0;
  gt_error_check(err);

  match.gen_range.start = 5;
  match.gen_range.end = 20;
  match.ref_range.start = 5;
  match.ref_range.end = 18;

  jump_table = gth_jump_table_new(&match, 1, false);
  gen_ranges = gt_array_new(sizeof (GtRange));
  gen_range.start = 0;
  gen_range.end = 24;
  gt_array_add(gen_ranges, gen_range);
  jt = gth_jump_table_make_tab(jump_table, gen_ranges, 25, 0, 5, false);

  gt_ensure(jt[0].from == 9 && jt[0].to == 0);
  gt_ensure(jt[1].from == 9 && jt[1].to == 0);
  gt_ensure(jt[2].from == 9 && jt[2].to == 0);
  gt_ensure(jt[3].from == 9 && jt[3].to == 0);
  gt_ensure(jt[4].from == 9 && jt[4].to == 0);
  gt_ensure(jt[5].from == 18 && jt[5].to == 0);
  gt_ensure(jt[6].from == 18 && jt[6].to == 0);
  gt_ensure(jt[7].from == 18 && jt[7].to == 0);
  gt_ensure(jt[8].from == 18 && jt[8].to == 0);
  gt_ensure(jt[9].from == 18 && jt[9].to == 5);
  gt_ensure(jt[10].from == 18 && jt[10].to == 5);
  gt_ensure(jt[11].from == 18 && jt[11].to == 5);
  gt_ensure(jt[12].from == 18 && jt[12].to == 5);
  gt_ensure(jt[13].from == 18 && jt[13].to == 5);
  gt_ensure(jt[14].from == 18 && jt[14].to == 5);
  gt_ensure(jt[15].from == 18 && jt[15].to == 5);
  gt_ensure(jt[16].from == 24 && jt[16].to == 5);
  gt_ensure(jt[17].from == 24 && jt[17].to == 5);
  gt_ensure(jt[18].from == 24 && jt[18].to == 5);
  gt_ensure(jt[19].from == 24 && jt[19].to == 5);
  gt_ensure(jt[20].from == 24 && jt[20].to == 14);
  gt_ensure(jt[21].from == 24 && jt[21].to == 14);
  gt_ensure(jt[22].from == 24 && jt[22].to == 14);
  gt_ensure(jt[23].from == 24 && jt[23].to == 14);
  gt_ensure(jt[24].from == 24 && jt[24].to == 14);

  full_ri = gth_jump_tab_get_full_row_info(jt, 25);

  gt_ensure(full_ri[0].offset == 0 && full_ri[0].length == 11);
  gt_ensure(full_ri[1].offset == 0 && full_ri[1].length == 11);
  gt_ensure(full_ri[2].offset == 0 && full_ri[2].length == 11);
  gt_ensure(full_ri[3].offset == 0 && full_ri[3].length == 11);
  gt_ensure(full_ri[4].offset == 0 && full_ri[4].length == 11);
  gt_ensure(full_ri[5].offset == 0 && full_ri[5].length == 11);
  gt_ensure(full_ri[6].offset == 0 && full_ri[6].length == 20);
  gt_ensure(full_ri[7].offset == 0 && full_ri[7].length == 20);
  gt_ensure(full_ri[8].offset == 0 && full_ri[8].length == 20);
  gt_ensure(full_ri[9].offset == 0 && full_ri[9].length == 20);
  gt_ensure(full_ri[10].offset == 0 && full_ri[10].length == 20);
  gt_ensure(full_ri[11].offset == 6 && full_ri[11].length == 14);
  gt_ensure(full_ri[12].offset == 6 && full_ri[12].length == 14);
  gt_ensure(full_ri[13].offset == 6 && full_ri[13].length == 14);
  gt_ensure(full_ri[14].offset == 6 && full_ri[14].length == 14);
  gt_ensure(full_ri[15].offset == 6 && full_ri[15].length == 14);
  gt_ensure(full_ri[16].offset == 6 && full_ri[16].length == 14);
  gt_ensure(full_ri[17].offset == 6 && full_ri[17].length == 20);
  gt_ensure(full_ri[18].offset == 6 && full_ri[18].length == 20);
  gt_ensure(full_ri[19].offset == 6 && full_ri[19].length == 20);
  gt_ensure(full_ri[20].offset == 6 && full_ri[20].length == 20);
  gt_ensure(full_ri[21].offset == 6 && full_ri[21].length == 20);
  gt_ensure(full_ri[22].offset == 15 && full_ri[22].length == 11);
  gt_ensure(full_ri[23].offset == 15 && full_ri[23].length == 11);
  gt_ensure(full_ri[24].offset == 15 && full_ri[24].length == 11);
  gt_ensure(full_ri[25].offset == 15 && full_ri[25].length == 11);

  half_ri = gth_jump_tab_get_row_info(jt, 25, &size);

  gt_ensure(half_ri[0].offset == 0 && half_ri[0].length == 11);
  gt_ensure(half_ri[1].offset == 0 && half_ri[1].length == 11);
  gt_ensure(half_ri[2].offset == 0 && half_ri[2].length == 11);
  gt_ensure(half_ri[3].offset == 0 && half_ri[3].length == 20);
  gt_ensure(half_ri[4].offset == 0 && half_ri[4].length == 20);
  gt_ensure(half_ri[5].offset == 0 && half_ri[5].length == 20);
  gt_ensure(half_ri[6].offset == 6 && half_ri[6].length == 14);
  gt_ensure(half_ri[7].offset == 6 && half_ri[7].length == 14);
  gt_ensure(half_ri[8].offset == 6 && half_ri[8].length == 20);
  gt_ensure(half_ri[9].offset == 6 && half_ri[9].length == 20);
  gt_ensure(half_ri[10].offset == 6 && half_ri[10].length == 20);
  gt_ensure(half_ri[11].offset == 15 && half_ri[11].length == 11);
  gt_ensure(half_ri[12].offset == 15 && half_ri[12].length == 11);
  gt_ensure(size == 203);

  gt_free(half_ri);
  gt_free(full_ri);
  gt_free(jt);
  gt_array_delete(gen_ranges);
  gth_jump_table_delete(jump_table);

  return had_err;
}

int gth_jump_table_unit_test(GtError *err)
{
  int had_err;
  gt_error_check(err);
  had_err = simple_balanced_no_overlap(err);
  if (!had_err)
    had_err = simple_unbalanced_no_overlap(err);
  if (!had_err)
    had_err = simple_cluster_no_overlap(err);
  if (!had_err)
    had_err = simple_balanced_small_overlap(err);
  if (!had_err)
    had_err = simple_unbalanced_small_overlap(err);
  if (!had_err)
    had_err = simple_cluster_small_overlap(err);
  if (!had_err)
    had_err = simple_balanced_large_overlap(err);
  if (!had_err)
    had_err = simple_unbalanced_large_overlap(err);
  if (!had_err)
    had_err = simple_cluster_large_overlap(err);
  if (!had_err)
    had_err = double_balanced_no_overlap(err);
  if (!had_err)
    had_err = double_unbalanced_no_overlap(err);
  if (!had_err)
    had_err = double_cluster_no_overlap(err);
  if (!had_err)
    had_err = double_balanced_with_overlap(err);
  if (!had_err)
    had_err = double_unbalanced_with_overlap(err);
  if (!had_err)
    had_err = double_cluster_with_overlap(err);
  if (!had_err)
    had_err = row_info_simple_balanced_no_overlap(err);
  if (!had_err)
    had_err = row_info_simple_unbalanced_with_overlap(err);
  return had_err;
}

void gth_jt_show(GthJumpTab *jt, GtUword from, GtUword to)
{
  GtUword i;
  gt_assert(jt);
  putchar('\n');
  for (i = from; i <= to; i++)
    printf("jt["GT_WU"].from="GT_WU", jt["GT_WU"].to="GT_WU"\n",
           i, jt[i].from, i, jt[i].to);
}

void gth_ri_show(GtRowInfo *ri, GtUword from, GtUword to)
{
  GtUword i;
  gt_assert(ri);
  putchar('\n');
  for (i = from; i <= to; i++) {
    printf("ri["GT_WU"].offset="GT_WU", ri["GT_WU"].length="GT_WU"\n",
           i, ri[i].offset, i, ri[i].length);
  }
}
