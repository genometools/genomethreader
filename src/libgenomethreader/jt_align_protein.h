#ifndef JT_ALIGN_PROTEIN_H
#define JT_ALIGN_PROTEIN_H

void gth_protein_complete_path_matrix_jt(GthDPtables *dpm,
                                         GthAlignInputProtein *input,
                                         bool proteinexonpenal,
                                         const unsigned char
                                         *gen_seq_tran,
                                         GtUword gen_dp_length,
                                         GtUword ref_dp_length,
                                         GthDPParam *dp_param,
                                         GthDPOptionsCore
                                         *dp_options_core,
                                         GthDPScoresProtein
                                         *dp_scores_protein,
                                         GthJumpTable *jump_table,
                                         GtArray *gen_ranges,
                                         GtUword ref_offset);

#endif
