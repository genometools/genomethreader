#include <ctype.h>
#include "core/cstr_api.h"
#include "core/ma.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "gth/seq_con_rep.h"
#include "libgenomethreader/seq_con_multiseq.h"
#include "types.h"
#include "virtualdef.h"
#include "alphabet.pr"
#include "readmulti.pr"
#include "readvirt.pr"
#include "multiseq.pr"
#include "multiseq-adv.pr"

#ifndef NOLICENSEMANAGER
#include "licensemanager.h"
#endif

struct GthSeqConMultiseq {
  const GthSeqCon parent_instance;

  char *indexname;

  /* vstree stuff */
  Alphabet *alpha;
  BOOL hasspecialsymbols;
  Multiseq *multiseq;

  bool owns_vstree_stuff;
  GtAlphabet *alphabet;
};

#define gth_seq_con_multiseq_cast(SC)\
        gth_seq_con_cast(gth_seq_con_multiseq_class(), SC)

static void multiseq_find_boundaries(Multiseq *multiseq, GtUword idx,
                                     PairUint *range)
{
  gt_assert(multiseq && range);
  if (findboundaries(multiseq, idx, range)) {
    fprintf(stderr,"%s\n", messagespace());
    exit(EXIT_FAILURE);
  }
}

static void gth_seq_con_multiseq_demand_orig_seq(GthSeqCon *sc)
{
  GthSeqConMultiseq *scm = gth_seq_con_multiseq_cast(sc);
  GtStr *tmpfilename;
  Uint numofbytes;
  gt_assert(scm);
  gt_assert(scm->multiseq && !scm->multiseq->originalsequence);
  gt_assert(scm->indexname);
  tmpfilename = gt_str_new_cstr(scm->indexname);
  gt_str_append_cstr(tmpfilename, ".ois");
  scm->multiseq->originalsequence =
    (Uchar*) CREATEMEMORYMAP(gt_str_get(tmpfilename), False, &numofbytes);
  gt_assert(scm->multiseq->originalsequence);
  gt_assert(numofbytes == scm->multiseq->totallength * (Uint) sizeof (Uchar));
  gt_str_delete(tmpfilename);
}

static GtUchar* gth_seq_con_multiseq_get_orig_seq(GthSeqCon *sc,
                                                  GtUword seq_num)
{
  PairUint range;
  GthSeqConMultiseq *scm = gth_seq_con_multiseq_cast(sc);
  multiseq_find_boundaries(scm->multiseq, seq_num, &range);
  return scm->multiseq->originalsequence + range.uint0;
}

static GtUchar* gth_seq_con_multiseq_get_tran_seq(GthSeqCon *sc,
                                                  GtUword seq_num)
{
  PairUint range;
  GthSeqConMultiseq *scm = gth_seq_con_multiseq_cast(sc);
  multiseq_find_boundaries(scm->multiseq, seq_num, &range);
  return scm->multiseq->sequence + range.uint0;
}

static GtUchar* gth_seq_con_multiseq_get_orig_seq_rc(GthSeqCon *sc,
                                                     GtUword seq_num)
{
  PairUint range;
  GthSeqConMultiseq *scm = gth_seq_con_multiseq_cast(sc);
  multiseq_find_boundaries(scm->multiseq, seq_num, &range);
  return scm->multiseq->rcoriginalsequence + range.uint0;
}

static GtUchar* gth_seq_con_multiseq_get_tran_seq_rc(GthSeqCon *sc,
                                                     GtUword seq_num)
{
  PairUint range;
  GthSeqConMultiseq *scm = gth_seq_con_multiseq_cast(sc);
  multiseq_find_boundaries(scm->multiseq, seq_num, &range);
  return scm->multiseq->rcsequence + range.uint0;
}

static void gth_seq_con_multiseq_get_description(GthSeqCon *sc,
                                                 GtUword seq_num,
                                                 GtStr *desc)
{
  GthSeqConMultiseq *scm;
  GtUword desclen;
  char *descptr;
  gt_assert(desc);
  scm = gth_seq_con_multiseq_cast(sc);
  desclen = DESCRIPTIONLENGTH(scm->multiseq, seq_num);
  descptr = (char*) DESCRIPTIONPTR(scm->multiseq, seq_num);
  gt_assert(desclen);
  desclen--; /* -1 to ignore '\n' */
  /* if the description ends with "\r\n" we ignore the '\r', too */
  if (desclen > 0 && descptr[desclen-1] == '\r')
    desclen--;
  gt_str_append_cstr_nt(desc, descptr, desclen);
}

static void gth_seq_con_multiseq_echo_description(GthSeqCon *sc,
                                                  GtUword seq_num,
                                                  GtFile *outfp)
{
  GtUword i, len;
  char *desc;
  GthSeqConMultiseq *scm = gth_seq_con_multiseq_cast(sc);
  len  = DESCRIPTIONLENGTH(scm->multiseq, seq_num);
  desc = (char*) DESCRIPTIONPTR(scm->multiseq, seq_num);
  gt_assert(len);
  len--; /* -1 to ignore '\n' */
  /* if the description ends with "\r\n" we ignore the '\r', too */
  if (len > 0 && desc[len-1] == '\r')
    len--;
  for (i = 0; i < len; i++)
    gt_file_xfputc(desc[i],outfp);
}

static GtUword gth_seq_con_multiseq_num_of_seqs(GthSeqCon *sc)
{
  GthSeqConMultiseq *scm = gth_seq_con_multiseq_cast(sc);
  return scm->multiseq->numofsequences;
}

static GtUword gth_seq_con_multiseq_total_length(GthSeqCon *sc)
{
  GthSeqConMultiseq *scm = gth_seq_con_multiseq_cast(sc);
  return scm->multiseq->totallength;
}

static GtRange gth_seq_con_multiseq_get_range(GthSeqCon *sc,
                                              GtUword seq_num)
{
  GtRange absolute_range;
  PairUint range;
  GthSeqConMultiseq *scm = gth_seq_con_multiseq_cast(sc);
  gt_assert(seq_num < gth_seq_con_multiseq_num_of_seqs(sc));
  multiseq_find_boundaries(scm->multiseq, seq_num, &range);
  absolute_range.start = range.uint0;
  absolute_range.end = range.uint1;
  return absolute_range;
}

static GtAlphabet* fill_alphabet(Alphabet *alpha)
{
  GtAlphabet *alphabet = NULL;
  gt_assert(alpha);
  if (vm_isdnaalphabet(alpha))
    alphabet = gt_alphabet_new_dna();
  else if (vm_isproteinalphabet(alpha))
    alphabet = gt_alphabet_new_protein();
  else
    gt_assert(false);
  return alphabet;
}

static GtAlphabet* gth_seq_con_multiseq_get_alphabet(GthSeqCon *sc)
{
  GthSeqConMultiseq *scm = gth_seq_con_multiseq_cast(sc);
  if (!scm->alphabet)
    scm->alphabet = fill_alphabet(scm->alpha);
  return scm->alphabet;
}

static void gth_seq_con_multiseq_free(GthSeqCon *sc)
{
  GthSeqConMultiseq *scm;
  if (!sc) return;
  scm = gth_seq_con_multiseq_cast(sc);
  if (scm->owns_vstree_stuff) {
    freemultiseq(scm->multiseq);
    gt_free(scm->multiseq);
    gt_free(scm->alpha);
  }
  gt_alphabet_delete(scm->alphabet);
  gt_free(scm->indexname);
}

const GthSeqConClass* gth_seq_con_multiseq_class(void)
{
  static const GthSeqConClass *ssc = NULL;
  if (!ssc) {
    ssc = gth_seq_con_class_new(sizeof (GthSeqConMultiseq),
                                gth_seq_con_multiseq_demand_orig_seq,
                                gth_seq_con_multiseq_get_orig_seq,
                                gth_seq_con_multiseq_get_tran_seq,
                                gth_seq_con_multiseq_get_orig_seq_rc,
                                gth_seq_con_multiseq_get_tran_seq_rc,
                                gth_seq_con_multiseq_get_description,
                                gth_seq_con_multiseq_echo_description,
                                gth_seq_con_multiseq_num_of_seqs,
                                gth_seq_con_multiseq_total_length,
                                gth_seq_con_multiseq_get_range,
                                gth_seq_con_multiseq_get_alphabet,
                                gth_seq_con_multiseq_free);
  }
  return ssc;
}

GthSeqCon* gth_seq_con_multiseq_new(const char *indexname, bool assign_rc,
                                    bool orig_seq, bool tran_seq)
{
  GthSeqConMultiseq *scm;
  GthSeqCon *sc;
  Uint indexsize;

  DefinedUint longestptr;
  Uint prefixlengthptr, maxbranchdepthptr, sixframeindex;
  ArrayPairUint largelcpvalues;
  BOOL rcmindex;

  gt_assert(indexname);

#ifndef NOLICENSEMANAGER
  lm_license_check_b();
#endif

  sc = gth_seq_con_create(gth_seq_con_multiseq_class());
  scm = gth_seq_con_multiseq_cast(sc);
  scm->indexname = gt_cstr_dup(indexname);
  scm->alpha = gt_calloc(1, sizeof *scm->alpha);
  scm->multiseq = gt_calloc(1, sizeof *scm->multiseq);
  /* map alphabet */
  if (mapalphabetifyoucan(&scm->hasspecialsymbols, scm->alpha, indexname)) {
    /* XXX */
    fprintf(stderr,"%s\n", messagespace());
    exit(EXIT_FAILURE);
  }
  /* map multiseq */
  initmultiseq(scm->multiseq);
  INITARRAY(&largelcpvalues, PairUint);
  if (parseprojectfile(scm->multiseq, &longestptr, &prefixlengthptr,
                       &largelcpvalues, &maxbranchdepthptr, &rcmindex,
                       &sixframeindex, indexname)) {
    /* XXX */
    fprintf(stderr,"%s\n", messagespace());
    exit(EXIT_FAILURE);
  }
  FREEARRAY(&largelcpvalues, PairUint);
  if (mapmultiseqifyoucan(&indexsize, scm->multiseq, indexname,
                          tran_seq, /* TISTAB */
                          orig_seq, /* OISTAB */
                          true,     /* DESTAB */
                          true      /* SSPTAB */)) {
    /* XXX */
    fprintf(stderr,"%s\n", messagespace());
    exit(EXIT_FAILURE);
  }
  if (assign_rc) {
    if (tran_seq) {
      scm->multiseq->rcsequence = copymultiseqRC(scm->multiseq);
      if (!scm->multiseq->rcsequence) {
        /* XXX */
        fprintf(stderr,"%s\n", messagespace());
        exit(EXIT_FAILURE);
      }
    }
    if (orig_seq) {
      scm->multiseq->rcoriginalsequence = copymultiseqRCorig(scm->multiseq);
      if (!scm->multiseq->rcoriginalsequence) {
        /* XXX */
        fprintf(stderr,"%s\n", messagespace());
        exit(EXIT_FAILURE);
      }
    }
  }
  scm->owns_vstree_stuff = true;
  scm->alphabet = NULL;
  return sc;
}

GthSeqCon* gth_seq_con_multiseq_new_vstree(Virtualtree *vstree)
{
  GthSeqConMultiseq *scm;
  GthSeqCon *sc;
  gt_assert(vstree);
  sc = gth_seq_con_create(gth_seq_con_multiseq_class());
  scm = gth_seq_con_multiseq_cast(sc);
  scm->alpha = &vstree->alpha;
  scm->multiseq = &vstree->multiseq;
  scm->owns_vstree_stuff = false;
  scm->alphabet = NULL;
  return sc;
}
