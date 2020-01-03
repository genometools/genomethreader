#ifndef JT_H
#define JT_H

#include "core/array2dim_sparse_api.h"
#include "gth/jump_table.h"

typedef struct {
  GtUword from,
                to;
} GthJumpTab;

GthJumpTable* gth_jump_table_new(GthJTMatch *matches,
                                 GtUword num_of_matches, bool debug);
GthJumpTable* gth_jump_table_new_reverse(const GthJumpTable*,
                                         GtUword gen_total_length,
                                         GtUword gen_offset,
                                         GtUword ref_total_length,
                                         GtUword ref_offset);
void          gth_jump_table_delete(GthJumpTable*);
GthJumpTab*   gth_jump_table_make_tab(const GthJumpTable*, GtArray *gen_ranges,
                                      GtUword ref_dp_length,
                                      GtUword ref_offset,
                                      GtUword overlap, bool debug);
GtRowInfo*    gth_jump_tab_get_row_info(const GthJumpTab*,
                                        GtUword gen_dp_length,
                                        GtUword *size);
GtRowInfo*    gth_jump_tab_get_full_row_info(const GthJumpTab*,
                                             GtUword gen_dp_length);
int           gth_jump_table_unit_test(GtError*);

void          gth_jt_show(GthJumpTab* jt, GtUword from, GtUword to);
void          gth_ri_show(GtRowInfo* ri, GtUword from, GtUword to);

#endif
