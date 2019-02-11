#ifndef JT_ALIGN_DNA_H
#define JT_ALIGN_DNA_H

void gth_dna_complete_path_matrix_jt(GthDPMatrix *dpm,
                                     const unsigned char *gen_seq_tran,
                                     const unsigned char *ref_seq_tran,
                                     GtUword genomic_offset,
                                     GtAlphabet *gen_alphabet,
                                     GthDPParam *dp_param,
                                     GthDPOptionsEST *dp_options_est,
                                     GthDPOptionsCore *dp_options_core,
                                     GthJumpTable *jump_table,
                                     GtArray *gen_ranges,
                                     GtUword ref_dp_length,
                                     GtUword ref_offset,
                                     GthPathMatrix **pm);

#endif
