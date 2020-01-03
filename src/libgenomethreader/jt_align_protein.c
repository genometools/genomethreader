#include "core/assert_api.h"
#include "core/codon_api.h"
#include "core/divmodmul_api.h"
#include "core/ma_api.h"
#include "core/safearith_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "gth/gthenum.h"
#include "gth/gtherror.h"
#include "gth/gthoutput.h"
#include "gth/gthstopcodon.h"
#include "gth/align_common.h"
#include "gth/align_protein_imp.h"
#include "gth/compute_scores.h"
#include "gth/dp_param.h"
#include "gth/gthtravalign.h"
#include "libgenomethreader/jt.h"
#include "libgenomethreader/jt_align_protein.h"

/* XXX: derived code from GenomeTools: src/gth/align_protein.c */
/* the following function evaluate the dynamic programming tables */
void gth_protein_complete_path_matrix_jt(GT_UNUSED GthDPtables *dpm,
                                         GT_UNUSED GthAlignInputProtein *input,
                                         GT_UNUSED bool proteinexonpenal,
                                         GT_UNUSED const unsigned char
                                                   *gen_seq_tran,
                                         GT_UNUSED GtUword gen_dp_length,
                                         GT_UNUSED GtUword ref_dp_length,
                                         GT_UNUSED GthDPParam *dp_param,
                                         GT_UNUSED GthDPOptionsCore
                                                   *dp_options_core,
                                         GT_UNUSED GthDPScoresProtein
                                                   *dp_scores_protein,
                                         GT_UNUSED GthJumpTable *jump_table,
                                         GT_UNUSED GtArray *gen_ranges,
                                         GT_UNUSED GtUword ref_offset)
{
#if 0
  GtUword n, m, modn, modnminus1, modnminus2, modnminus3;
  unsigned char origreferencechar;
  GthFlt value, maxvalue;
  bool jumped = false;
  GthPath retrace;
  GthJumpTab *jt;

  gt_assert(jump_table);

  jt = gth_jump_table_make_tab(jump_table, gen_ranges, ref_dp_length,
                               ref_offset, dp_options_core->jtoverlap,
                               dp_options_core->jtdebug);

  /* stepping along the genomic sequence */
  m = REFERENCEDPSTART;
  for (n = GENOMICDPSTART; n <= gen_dp_length; n++) {
    modn       = GT_MOD4(n),
    modnminus1 = GT_MOD4(n-1),
    modnminus2 = GT_MOD4(n-2),
    modnminus3 = GT_MOD4(n-3);

    PATH(E_STATE, n, 0)  = (GthPath) E_N1;
    PATH(IA_STATE, n, 0) = (GthPath) IA_N1;
    PATH(IB_STATE, n, 0) = (GthPath) IB_N1;
    PATH(IC_STATE, n, 0) = (GthPath) IC_N1;

    /* stepping along the protein sequence */
    for (; m <= ref_dp_length; m++) {
      origreferencechar = input->ref_seq_orig[m-1];

      /* evaluate E_nm */
      /* 0. */
      maxvalue = SCORE(E_STATE, modnminus3, m-1) +
                 /* XXX: why is here no extra condition? */
                 (dp_param->log_1minusPdonor[n-3] +
                  GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-3],
                              gen_seq_tran[n-2], gen_seq_tran[n-1],
                              origreferencechar));
      retrace  = (GthPath) E_N3M;

      /* 1. */
      value = SCORE(E_STATE, modnminus2, m-1);
      if (n < gen_dp_length || m < WSIZE_PROTEIN) {
        value += dp_param->log_1minusPdonor[n-2] +
                 GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-2],
                             gen_seq_tran[n-1], DASH, origreferencechar);
      }
      UPDATEMAX(E_N2M);

      /* 2. */
      value = SCORE(E_STATE, modnminus1, m-1);
      if (n < gen_dp_length || m < WSIZE_PROTEIN) {
        value += dp_param->log_1minusPdonor[n-1] +
                 GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-1], DASH, DASH,
                             origreferencechar);
      }
      UPDATEMAX(E_N1M);

      if (jumped) {
        gt_assert(m != REFERENCEDPSTART && n != GENOMICDPSTART);
      /* 3. */
      value = SCORE(E_STATE, modn, m-1);
      if (n < gen_dp_length || m < WSIZE_PROTEIN) {
        if (n == gen_dp_length) {
          /* in this case the value used in the 'else' branch below is not
             defined. */
          value += dp_param->log_1minusPdonor[n-1];
        }
        else {
          value += dp_param->log_1minusPdonor[n];
                                     /* XXX: ^^^  why n? */
        }
        value += GTHGETSCORE(dp_scores_protein, DASH, DASH, DASH,
                             origreferencechar);
      }
      UPDATEMAX(E_M);
      }
      else
        jumped = false;

      if (jt[n-2].from + 1 >= m || m == REFERENCEDPSTART ||
          n == GENOMICDPSTART) {
      /* 4. */
      value = SCORE(E_STATE, modnminus3, m);
      if (m < ref_dp_length || n < WSIZE_DNA) {
        value += dp_param->log_1minusPdonor[n-3] +
                 GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-3],
                             gen_seq_tran[n-2], gen_seq_tran[n-1], DASH);
      }
      UPDATEMAX(E_N3);

      /* 5. */
      value = SCORE(E_STATE, modnminus2, m);
      if (m < ref_dp_length || n < WSIZE_DNA) {
        value += dp_param->log_1minusPdonor[n-2] +
                 GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-2],
                             gen_seq_tran[n-1], DASH, DASH);
      }
      UPDATEMAX(E_N2);

      /* 6. */
      value = SCORE(E_STATE, modnminus1, m);
      if (m < ref_dp_length || n < WSIZE_DNA) {
        value += dp_param->log_1minusPdonor[n-1] +
                 GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-1], DASH, DASH,
                             DASH);
      }
      UPDATEMAX(E_N1);
      }

      /* 7. */
      value = SCORE(IA_STATE, modnminus3, m-1);
      if (n > GENOMICDPSTART) /* the value below is only defined in this case */
        value += dp_param->log_Pacceptor[n-4];
      value += GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-3],
                           gen_seq_tran[n-2], gen_seq_tran[n-1],
                           origreferencechar);
      if (n - 2 - dpm->intronstart_A[modnminus3][m-1] <
          dp_options_core->dpminintronlength) {
        value -= dp_options_core->shortintronpenalty;
      }
      UPDATEMAX(IA_N3M);

      /* 8.
         this recurrence is only used if an intron has already been introduced.
         (in this case "dpm->splitcodon_B[modnminus1][m-1]" is different from
         "UNSET". */
      if (dpm->splitcodon_B[modnminus1][m-1] != (unsigned char) UNSET) {
        value = SCORE(IB_STATE, modnminus2, m-1) +
                dp_param->log_Pacceptor[n-3] +
                GTHGETSCORE(dp_scores_protein,
                            dpm->splitcodon_B[modnminus2][m-1],
                            gen_seq_tran[n-2], gen_seq_tran[n-1],
                            origreferencechar);
        if (n - 1 - dpm->intronstart_B[modnminus2][m-1] <
            dp_options_core->dpminintronlength) {
          value -= dp_options_core->shortintronpenalty;
        }
        UPDATEMAX(IB_N2M);
      }

      /* 9.
         explanation for check see above.
         "dpm->splitcodon_C2[modnminus1][m-1]" needs not to be checked,
         because it is always set in conjunction with
         "dpm->splitcodon_C1[modnminus1][m-1]". */
      if (dpm->splitcodon_C1[modnminus1][m-1] != (unsigned char) UNSET) {
        value = SCORE(IC_STATE, modnminus1, m-1) +
                dp_param->log_Pacceptor[n-2] +
                GTHGETSCORE(dp_scores_protein,
                            dpm->splitcodon_C1[modnminus1][m-1],
                            dpm->splitcodon_C2[modnminus1][m-1],
                            gen_seq_tran[n-1], origreferencechar);
        if (n - dpm->intronstart_C[modnminus1][m-1] <
            dp_options_core->dpminintronlength) {
          value -= dp_options_core->shortintronpenalty;
        }
        UPDATEMAX(IC_N1M);
      }

      /* save maximum values */
      SCORE(E_STATE, modn, m) = maxvalue;
      PATH(E_STATE, n, m)     = retrace;

      if (proteinexonpenal) {
        switch (retrace) {
          case E_N3M:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus3][m-1];
            break;
          case E_N2M:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus2][m-1];
            break;
          case E_N1M:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus1][m-1];
            break;
          case E_M:
            dpm->exonstart[modn][m] = dpm->exonstart[modn][m-1];
            break;
          case E_N3:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus3][m];
            break;
          case E_N2:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus2][m];
            break;
          case E_N1:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus1][m];
            break;
          case IA_N3M:
          case IB_N2M:
          case IC_N1M:
            dpm->exonstart[modn][m] = n;
            break;
          case IC_N1:
            dpm->exonstart[modn][m] = n;
            break;
          default: gt_assert(0);
        }
      }

      if (jt[n-2].from + 1 >= m || m == REFERENCEDPSTART ||
          n == GENOMICDPSTART) {
      /* evaluate IA_nm */
      maxvalue = SCORE(IA_STATE, modnminus1, m);
      if (!dp_options_core->freeintrontrans)
        maxvalue += dp_param->log_1minusPacceptor[n-2];
      retrace  = (GthPath) IA_N1;

      value = SCORE(E_STATE, modnminus1, m) + dp_param->log_Pdonor[n-1];
      if (proteinexonpenal) {
        if (n - dpm->exonstart[modnminus1][m] <
            dp_options_core->dpminexonlength) {
          value -= dp_options_core->shortexonpenalty;
        }
      }
      UPDATEMAX(E_N1);

      /* save maximum values */
      SCORE(IA_STATE, modn, m) = maxvalue;
      PATH(IA_STATE, n, m) = retrace;

      switch (retrace) {
        case IA_N1:
          dpm->intronstart_A[modn][m] = dpm->intronstart_A[modnminus1][m];
          break;
        case E_N1:
          dpm->intronstart_A[modn][m] = n;
          break;
        default: gt_assert(0);
      }

      /* evaluate IB_nm */
      maxvalue = SCORE(IB_STATE, modnminus1, m);
      if (!dp_options_core->freeintrontrans)
        maxvalue += dp_param->log_1minusPacceptor[n-2];
      retrace  = (GthPath) IB_N1;

      value = SCORE(E_STATE, modnminus2, m) + dp_param->log_Pdonor[n-1];
      if (proteinexonpenal) {
        if (n - 1 - dpm->exonstart[modnminus2][m] <
            dp_options_core->dpminexonlength) {
          value -= dp_options_core->shortexonpenalty;
        }
      }
      UPDATEMAX(E_N2);

      /* save maximum values */
      SCORE(IB_STATE, modn, m) = maxvalue;
      PATH(IB_STATE, n, m) = retrace;

      switch (retrace) {
        case(IB_N1):
          dpm->intronstart_B[modn][m] = dpm->intronstart_B[modnminus1][m];
          dpm->splitcodon_B[modn][m]  = dpm->splitcodon_B[modnminus1][m];
          break;
        case(E_N2):
          dpm->intronstart_B[modn][m] = n;
          dpm->splitcodon_B[modn][m]  = gen_seq_tran[n-2];
          break;
        default: gt_assert(0);
      }

      /* evaluate IC_nm */
      maxvalue = SCORE(IC_STATE, modnminus1, m);
      if (!dp_options_core->freeintrontrans)
        maxvalue += dp_param->log_1minusPacceptor[n-2];
      retrace = (GthPath) IC_N1;

      value = SCORE(E_STATE, modnminus3, m) + dp_param->log_Pdonor[n-1];
      if (proteinexonpenal) {
        if (n - 2 - dpm->exonstart[modnminus3][m] <
            dp_options_core->dpminexonlength) {
          value -= dp_options_core->shortexonpenalty;
        }
      }
      UPDATEMAX(E_N3);

      /* save maximum values */
      SCORE(IC_STATE, modn, m) = maxvalue;
      PATH(IC_STATE, n, m) = retrace;

      switch (retrace) {
        case(IC_N1):
          dpm->intronstart_C[modn][m] = dpm->intronstart_C[modnminus1][m];
          dpm->splitcodon_C1[modn][m] = dpm->splitcodon_C1[modnminus1][m];
          dpm->splitcodon_C2[modn][m] = dpm->splitcodon_C2[modnminus1][m];
          break;
        case(E_N3):
          dpm->intronstart_C[modn][m] = n;
          dpm->splitcodon_C1[modn][m] = gen_seq_tran[n-3];
          dpm->splitcodon_C2[modn][m] = gen_seq_tran[n-2];
          break;
        default: gt_assert(0);
      }
      }

      if (jt[n-1].from == m-1 && m != REFERENCEDPSTART && n != GENOMICDPSTART) {
        m = jt[n-1].to + 1;
        jumped = true;
        break;
      }
    }
  }

  gt_free(jt);
#endif
}
