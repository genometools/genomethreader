#ifndef SEQ_CON_MULTISEQ_H
#define SEQ_CON_MULTISEQ_H

#include "gth/seq_con.h"
#include "alphadef.h"
#include "multidef.h"
#include "virtualdef.h"

typedef struct GthSeqConMultiseq GthSeqConMultiseq;

const GthSeqConClass* gth_seq_con_multiseq_class(void);
GthSeqCon*            gth_seq_con_multiseq_new(const char *indexname,
                                               bool assign_rc, bool orig_seq,
                                               bool tran_seq);
GthSeqCon*            gth_seq_con_multiseq_new_vstree(Virtualtree *vstree);

#endif
