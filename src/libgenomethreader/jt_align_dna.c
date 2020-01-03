#include <math.h>
#include "core/array2dim_api.h"
#include "core/divmodmul_api.h"
#include "core/undef_api.h"
#include "gth/align_dna_imp.h"
#include "gth/compute_scores.h"
#include "libgenomethreader/jt.h"
#include "libgenomethreader/jt_align_dna.h"

#ifndef NOLICENSEMANAGER
#include "licensemanager.h"
#endif

/* XXX: duplicated code from GenomeTools: src/gth/align_dna.c */
/* the following function evaluates state E_1m for indices 1 and <m> and
   stores a backtrace reference */
static void E_1m(GthDPMatrix *dpm, unsigned char genomicchar,
                 const unsigned char *ref_seq_tran, GtUword m,
                 GtAlphabet *gen_alphabet, GthDbl log_probies,
                 GthDPOptionsEST *dp_options_est,
                 GthDPOptionsCore *dp_options_core)
{
  GthFlt value, maxvalue;
  GthPath retrace;
  GthDbl rval = 0.0, outputweight = 0.0;
  unsigned char referencechar = ref_seq_tran[m-1];
  unsigned int gen_alphabet_mapsize = gt_alphabet_size(gen_alphabet);

  /* 0. */
  rval = log_probies;
  ADDOUTPUTWEIGHT(rval,genomicchar,referencechar);
  if ((m < dp_options_est->wdecreasedoutput ||
       m > dpm->ref_dp_length - dp_options_est->wdecreasedoutput) &&
       genomicchar == referencechar) {
    ADDOUTPUTWEIGHT(outputweight, genomicchar, referencechar);
    rval -= (outputweight / 2.0);
  }
  maxvalue = (GthFlt) (dpm->score[DNA_E_STATE][0][m-1] + rval);
  retrace  = DNA_E_NM;

  /* 1. */
  outputweight = 0.0;
  rval = log_probies;
  ADDOUTPUTWEIGHT(rval,genomicchar,referencechar);
  if ((m < dp_options_est->wdecreasedoutput ||
       m > dpm->ref_dp_length - dp_options_est->wdecreasedoutput) &&
       genomicchar == referencechar) {
    ADDOUTPUTWEIGHT(outputweight, genomicchar, referencechar);
    rval -= (outputweight / 2.0);
  }
  value = (GthFlt) (dpm->score[DNA_I_STATE][0][m-1] + rval);
  /* intron from intronstart to n-1 => n-1 - intronstart + 1 */
  if (1 - dpm->intronstart[0][m - 1] < dp_options_core->dpminintronlength)
    value -= dp_options_core->shortintronpenalty;
  UPDATEMAX(DNA_I_NM);

  /* 2. */
  rval = log_probies;
  if (m < dpm->ref_dp_length) {
    ADDOUTPUTWEIGHT(rval, genomicchar, (unsigned char) DASH);
  }
  value = (GthFlt) (dpm->score[DNA_E_STATE][0][m] + rval);
  UPDATEMAX(DNA_E_N);

  /* 3. */
  rval = log_probies;
  if (m < dpm->ref_dp_length) {
    ADDOUTPUTWEIGHT(rval, genomicchar, (unsigned char) DASH);
  }
  value = (GthFlt) (dpm->score[DNA_I_STATE][0][m] + rval);
  /* intron from intronstart to n-1 => n-1 - intronstart + 1 */
  if (1 - dpm->intronstart[0][m] < dp_options_core->dpminintronlength)
    value -= dp_options_core->shortintronpenalty;
  UPDATEMAX(DNA_I_N);

  /* 4. */
  rval = log_probies;
  ADDOUTPUTWEIGHT(rval, (unsigned char) DASH, referencechar);
  value = (GthFlt) (dpm->score[DNA_E_STATE][1][m-1] + rval);
  UPDATEMAX(DNA_E_M);

  /* 5. */
  rval = dpm->score[DNA_I_STATE][1][m-1] + log_probies;
  ADDOUTPUTWEIGHT(rval, (unsigned char) DASH, referencechar);
  value = (GthFlt) (dpm->score[DNA_I_STATE][1][m-1] + rval);
  /* intron from intronstart to n => n - intronstart + 1 */
  if (1 - dpm->intronstart[1][m - 1] + 1 < dp_options_core->dpminintronlength)
    value -= dp_options_core->shortintronpenalty;
  UPDATEMAX(DNA_I_M);

  /* save maximum values */
  dpm->score[DNA_E_STATE][1][m] = maxvalue;
  dpm->path[0][m] |= (retrace << 4);
  if (dpm->path_jt) {
    dpm->path_jt[0][m] |= (retrace << 4);
    gt_assert(dpm->path[0][m] == dpm->path_jt[0][m]);
  }

  switch (retrace) {
    case DNA_I_NM:
    case DNA_I_N:
    case DNA_I_M:
      dpm->exonstart[1][m] = 1;
      break;
    case DNA_E_NM:
      dpm->exonstart[1][m] = dpm->exonstart[0][m - 1];
      break;
    case DNA_E_N:
      dpm->exonstart[1][m] = dpm->exonstart[0][m];
      break;
    case DNA_E_M:
      dpm->exonstart[1][m] = dpm->exonstart[1][m - 1];
      break;
    default: gt_assert(0);
  }
}

/* XXX: duplicated code from GenomeTools: src/gth/align_dna.c */
/* the following function evaluates state I_1m for indices 1 and <m> and
   stores a backtrace reference */
static void I_1m(GthDPMatrix *dpm, GtUword m, GthDbl log_1minusprobies)
{
  GthFlt value, maxvalue;
  GthPath retrace;

  /* 0. */
  maxvalue = (GthFlt) (dpm->score[DNA_E_STATE][0][m] + log_1minusprobies);
  retrace  = I_STATE_E_N;

  /* 1. */
  value = (GthFlt) (dpm->score[DNA_I_STATE][0][m] + log_1minusprobies);
  UPDATEMAX(I_STATE_I_N);

  /* save maximum values */
  dpm->score[DNA_I_STATE][1][m] = maxvalue;
  dpm->path[0][m] |= (retrace << 4);
  if (dpm->path_jt) {
    dpm->path_jt[0][m] |= (retrace << 4);
    gt_assert(dpm->path[0][m] == dpm->path_jt[0][m]);
  }

  switch (retrace) {
    case I_STATE_E_N:
      /* begin of a new intron */
      dpm->intronstart[1][m] = 1;
      break;
    case I_STATE_I_N:
      /* continue existing intron */
      dpm->intronstart[1][m] = dpm->intronstart[0][m];
      break;
    default: gt_assert(0);
  }
}

/* XXX: derived code from GenomeTools: src/gth/align_dna.c */
/* the following function evaluate the dynamic programming tables */
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
                                     GthPathMatrix **pm)
{
  GthFlt value, maxvalue;
  GthPath retrace;
  GtUword n, m, modn, modnminus1;
  GthDbl rval, outputweight, **outputweights,
         log_probies,          /* initial exon state probability */
         log_1minusprobies;    /* initial intron state probability */
  GthFlt log_probdelgen,       /* deletion in genomic sequence */
         log_1minusprobdelgen;
  bool jumped = false;
  unsigned char genomicchar, referencechar;
  unsigned int gen_alphabet_mapsize = gt_alphabet_size(gen_alphabet);
  GthJumpTab *jt;

  gt_assert(dpm->gen_dp_length > 1);
  gt_assert(jump_table);

#ifndef NOLICENSEMANAGER
  lm_license_check_i();
#endif

  jt = gth_jump_table_make_tab(jump_table, gen_ranges, ref_dp_length,
                               ref_offset, dp_options_core->jtoverlap,
                               dp_options_core->jtdebug);

  GtUword size;
  GtRowInfo *ri = gth_jump_tab_get_row_info(jt, dpm->gen_dp_length, &size);

  if (dp_options_core->jtdebug) {
    gth_jt_show(jt, 0, dpm->gen_dp_length-1);
    gth_ri_show(ri, 0, GT_DIV2(dpm->gen_dp_length+1) +
                       GT_MOD2(dpm->gen_dp_length+1) - 1);
  }

  gt_array2dim_sparse_calloc(dpm->path_jt,
                             GT_DIV2(dpm->gen_dp_length + 1) +
                             GT_MOD2(dpm->gen_dp_length + 1), size, ri);

  dpm->path_jt[0][0]  = DNA_E_NM;
  dpm->path_jt[0][0] |= I_STATE_E_N;
  for (m = 1; m < ri[0].length; m++) {
    dpm->path_jt[0][m]  = DNA_E_M;
    dpm->path_jt[0][m] |= I_STATE_I_N;
  }
  /* XXX */
  /* gt_free(ri); */

  log_probies = (GthDbl) log((double) dp_options_est->probies);
  log_1minusprobies = (GthDbl) log(1.0 - dp_options_est->probies);
  log_probdelgen = (GthFlt) log((double) dp_options_est->probdelgen);
  log_1minusprobdelgen = (GthFlt) log(1.0 - dp_options_est->probdelgen);

  /* precompute outputweights
     XXX: move this to somewhere else, maybe make it smaller */
  gt_array2dim_calloc(outputweights, UCHAR_MAX+1, UCHAR_MAX+1);
  for (n = 0; n <= UCHAR_MAX; n++) {
    for (m = 0; m <= UCHAR_MAX; m++) {
      ADDOUTPUTWEIGHT(outputweights[n][m], n, m);
    }
  }

  gt_assert(dpm->path[0][0] == dpm->path_jt[0][0]);

  if (!genomic_offset) {
    /* handle case for n equals 1 */
    dpm->path[0][0] |= UPPER_E_N;
    dpm->path[0][0] |= UPPER_I_STATE_I_N;
    if (dpm->path_jt) {
      dpm->path_jt[0][0] |= UPPER_E_N;
      dpm->path_jt[0][0] |= UPPER_I_STATE_I_N;
      gt_assert(dpm->path[0][0] == dpm->path_jt[0][0]);
    }

    /* stepping along the cDNA/EST sequence */
    for (m = 1; m <= dpm->ref_dp_length; m++) {
      E_1m(dpm, gen_seq_tran[0], ref_seq_tran, m, gen_alphabet, log_probies,
           dp_options_est, dp_options_core);
      I_1m(dpm, m, log_1minusprobies);
      /* XXX */
      if (jt[0].from == m-1)
        break;
    }
  }

  /* handle all other n's
     stepping along the genomic sequence */
  if (genomic_offset)
    n = genomic_offset + 1;
  else
    n = 2;

  m = 1;
  for (; n <= dpm->gen_dp_length; n++) {
    modn = GT_MOD2(n);
    modnminus1 = GT_MOD2(n-1);
    genomicchar = gen_seq_tran[n-1];

    if (modn) {
      dpm->path[GT_DIV2(n)][0] |= UPPER_E_N;
      dpm->path[GT_DIV2(n)][0] |= UPPER_I_STATE_I_N;
      if (dpm->path_jt && jt[n-2].to == 0) {
        dpm->path_jt[GT_DIV2(n)][0] |= UPPER_E_N;
        dpm->path_jt[GT_DIV2(n)][0] |= UPPER_I_STATE_I_N;
        gt_assert(dpm->path[GT_DIV2(n)][0] == dpm->path_jt[GT_DIV2(n)][0]);
      }
    }
    else {
      dpm->path[GT_DIV2(n)][0]  = DNA_E_N;
      dpm->path[GT_DIV2(n)][0] |= I_STATE_I_N;
      if (dpm->path_jt && jt[n-2].to == 0) {
        dpm->path_jt[GT_DIV2(n)][0]  = DNA_E_N;
        dpm->path_jt[GT_DIV2(n)][0] |= I_STATE_I_N;
        gt_assert(dpm->path[GT_DIV2(n)][0] == dpm->path_jt[GT_DIV2(n)][0]);
      }
    }

    /* stepping along the cDNA/EST sequence */
    for (; m <= dpm->ref_dp_length; m++) {
      referencechar = ref_seq_tran[m-1];

    /* XXX */
#if 1
      gt_assert(m >= ri[GT_DIV2(n)].offset);
      gt_assert(m < ri[GT_DIV2(n)].offset + ri[GT_DIV2(n)].length);
#endif

      /* evaluate E_nm */

      if (n <= 2 || (n > 2 && (jt[n-3].to == 0 || jt[n-3].to + 1 < m)
                           && (jt[n-2].from + 1 + 1 >= m))) {
        /* 0. */
        outputweight = 0.0;
        rval = (GthDbl) (log_1minusprobdelgen +
                         dp_param->log_1minusPdonor[n-1]);
        rval += outputweights[genomicchar][referencechar];
        if ((m < dp_options_est->wdecreasedoutput ||
             m > dpm->ref_dp_length - dp_options_est->wdecreasedoutput) &&
             genomicchar == referencechar) {
          outputweight += outputweights[genomicchar][referencechar];
          rval -= (outputweight / 2.0);
        }
        maxvalue = (GthFlt) (dpm->score[DNA_E_STATE][modnminus1][m-1] + rval);
        retrace  = DNA_E_NM;

        /* 1. */
        outputweight = 0.0;
        rval = (GthDbl) (dp_param->log_Pacceptor[n-2] + log_1minusprobdelgen);
        rval += outputweights[genomicchar][referencechar];
        if ((m < dp_options_est->wdecreasedoutput ||
             m > dpm->ref_dp_length - dp_options_est->wdecreasedoutput) &&
             genomicchar == referencechar) {
          outputweight += outputweights[genomicchar][referencechar];
          rval -= (outputweight / 2.0);
        }
        value = (GthFlt) (dpm->score[DNA_I_STATE][modnminus1][m-1] + rval);
        /* intron from intronstart to n-1 => n-1 - intronstart + 1 */
        if (n - dpm->intronstart[modnminus1][m - 1] <
            dp_options_core->dpminintronlength) {
          value -= dp_options_core->shortintronpenalty;
        }
        UPDATEMAX(DNA_I_NM);
      }
      else {
        maxvalue = GTH_MINUSINFINITY;
        retrace  = DNA_E_NM;
      }

      if (jt[n-2].from + 1 >= m) {
        /* 2. */
        rval = 0.0;
        if (m < dpm->ref_dp_length || n < dp_options_est->wzerotransition)
          rval += (log_1minusprobdelgen + dp_param->log_1minusPdonor[n-1]);
        if (m < dpm->ref_dp_length)
          rval += outputweights[genomicchar][DASH];
        value = (GthFlt) (dpm->score[DNA_E_STATE][modnminus1][m] + rval);
        UPDATEMAX(DNA_E_N);

        /* 3. */
        rval = (GthDbl) (dp_param->log_Pacceptor[n-2] + log_1minusprobdelgen);
        if (m < dpm->ref_dp_length)
          rval += outputweights[genomicchar][DASH];
        value = (GthFlt) (dpm->score[DNA_I_STATE][modnminus1][m] + rval);
        /* intron from intronstart to n-1 => n-1 - intronstart + 1 */
        if (n - dpm->intronstart[modnminus1][m] <
            dp_options_core->dpminintronlength) {
          value -= dp_options_core->shortintronpenalty;
        }
        UPDATEMAX(DNA_I_N);
      }

      if (!jumped) {
        /* 4. */
        rval = 0.0;
        if (n < dpm->gen_dp_length || m < dp_options_est->wzerotransition)
          rval = (GthDbl) log_probdelgen;
        if (n < dpm->gen_dp_length)
          rval += outputweights[DASH][referencechar];
        value = (GthFlt) (dpm->score[DNA_E_STATE][modn][m-1] + rval);
        UPDATEMAX(DNA_E_M);

        /* 5. */
        rval = 0.0;
        if (n < dpm->gen_dp_length)
         rval += (dp_param->log_Pacceptor[n-1] + log_probdelgen);
        if (n < dpm->gen_dp_length)
          rval += outputweights[DASH][referencechar];
        value = (GthFlt) (dpm->score[DNA_I_STATE][modn][m-1] + rval);
        /* intron from intronstart to n => n - intronstart + 1 */
        if (n - dpm->intronstart[modn][m - 1] + 1 <
            dp_options_core->dpminintronlength) {
          value -= dp_options_core->shortintronpenalty;
        }
        UPDATEMAX(DNA_I_M);
      }
      else
        jumped = false;

      /* save maximum values */
      dpm->score[DNA_E_STATE][modn][m] = maxvalue;
      if (modn) {
        dpm->path[GT_DIV2(n)][m] |= (retrace << 4);
        if (dpm->path_jt) {
          dpm->path_jt[GT_DIV2(n)][m] |= (retrace << 4);
          gt_assert(dpm->path[GT_DIV2(n)][m] == dpm->path_jt[GT_DIV2(n)][m]);
        }
      }
      else {
        dpm->path[GT_DIV2(n)][m]  = retrace;
        if (dpm->path_jt) {
          dpm->path_jt[GT_DIV2(n)][m]  = retrace;
          gt_assert(dpm->path[GT_DIV2(n)][m] == dpm->path_jt[GT_DIV2(n)][m]);
        }
      }

      switch (retrace) {
        case DNA_I_NM:
        case DNA_I_N:
        case DNA_I_M:
          dpm->exonstart[modn][m] = n;
          break;
        case DNA_E_NM:
          dpm->exonstart[modn][m] = dpm->exonstart[modnminus1][m - 1];
          break;
        case DNA_E_N:
          dpm->exonstart[modn][m] = dpm->exonstart[modnminus1][m];
          break;
        case DNA_E_M:
          dpm->exonstart[modn][m] = dpm->exonstart[modn][m - 1];
          break;
        default: gt_assert(0);
      }

      /* evaluate I_nm */

      if (jt[n-2].from + 1 >= m) {
        /* 0. */
        maxvalue = dpm->score[DNA_E_STATE][modnminus1][m] +
                   (log_1minusprobdelgen + dp_param->log_Pdonor[n-1]);
        if (n - dpm->exonstart[modnminus1][m]
            < dp_options_core->dpminexonlength) {
           maxvalue -= dp_options_core->shortexonpenalty;
        }
        retrace  = I_STATE_E_N;

        /* 1. */
        value = dpm->score[DNA_I_STATE][modnminus1][m];
        if (!dp_options_core->freeintrontrans && m < dpm->ref_dp_length)
          value += dp_param->log_1minusPacceptor[n-2];
        UPDATEMAX(I_STATE_I_N);

        /* save maximum values */
        dpm->score[DNA_I_STATE][modn][m] = maxvalue;
        if (modn) {
          dpm->path[GT_DIV2(n)][m] |= (retrace << 4);
          if (dpm->path_jt) {
            dpm->path_jt[GT_DIV2(n)][m] |= (retrace << 4);
            gt_assert(dpm->path[GT_DIV2(n)][m] == dpm->path_jt[GT_DIV2(n)][m]);
          }
        }
        else {
          dpm->path[GT_DIV2(n)][m] |= retrace;
          if (dpm->path_jt) {
            dpm->path_jt[GT_DIV2(n)][m] |= retrace;
            gt_assert(dpm->path[GT_DIV2(n)][m] == dpm->path_jt[GT_DIV2(n)][m]);
          }
        }

        switch (retrace) {
          case I_STATE_E_N:
            /* begin of a new intron */
            dpm->intronstart[modn][m] = n;
            break;
          case I_STATE_I_N:
            /* continue existing intron */
            dpm->intronstart[modn][m] = dpm->intronstart[modnminus1][m];
            break;
          default: gt_assert(0);
        }
      }

      if (jt[n-1].from == m-1) {
        m = jt[n-1].to + 1;
        jumped = true;
        break;
      }
    }
  }

  if (dp_options_core->btmatrixgenrange.start != GT_UNDEF_ULONG) {
    GtRowInfo *full_ri = gth_jump_tab_get_full_row_info(jt, dpm->gen_dp_length);
    *pm = gth_path_matrix_new(dpm->path_jt, dpm->gen_dp_length,
                              dpm->ref_dp_length,
                              &dp_options_core->btmatrixgenrange,
                              &dp_options_core->btmatrixrefrange, full_ri);
    gt_free(full_ri);
  }

  /* free space  */
  gt_array2dim_delete(outputweights);
  gt_free(jt);
  gt_free(ri);
}
