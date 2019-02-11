#include <stdlib.h>
#include "core/types_api.h"
#include "core/unused_api.h"
#include "select.h"

/*
  This module selects all matches which start with the first
  character of the database sequence or end with the last character
  of the database sequence. To this end, we have to find the boundaries
  of a match, using the following function.
*/

/*
  Suppose we have \(k\) sequences numbered \(0,\ldots,k-1\).
  Given a sequence number \texttt{snum} in the range \([0,k-1]\),
  the function \texttt{polya_find_boundaries} computes in the integer pair
  \texttt{range}, the first and last index of \(T_{snum}\), where
  \(T_{snum}\) is the sequence with number \(snum\):
  \texttt{range->uint0} is the first index, and
  \texttt{range->uint1} is the last index with respect to
  the concatenation of all sequences. If the given sequence
  number is not valid, then the function terminates with an
  exit code.
*/

static void polya_find_boundaries(Multiseq *multiseq,Uint snum,PairUint *range)
{
  if (snum >= multiseq->numofsequences)
  {
    fprintf(stderr,"polyafunc: sequence "GT_WU" does not exist\n",
            (GtUword) snum);
    exit(EXIT_FAILURE);
  }
  if (snum == 0)
  {
    range->uint0 = 0;
    if (multiseq->numofsequences == 1)
    {
      range->uint1 = multiseq->totallength - 1;
    } else
    {
      range->uint1 = multiseq->markpos.spaceUint[snum] - 1;
    }
  } else
  {
    range->uint0 = multiseq->markpos.spaceUint[snum-1] + 1;
    if (snum == multiseq->numofsequences - 1)
    {
      range->uint1 = multiseq->totallength - 1;
    } else
    {
      range->uint1 = multiseq->markpos.spaceUint[snum] - 1;
    }
  }
}

/*
  The selection function bundle.
*/

/*
  The following function selects the matches as described above.
*/

Sint gthpolyaselectmatch(GT_UNUSED Alphabet *alpha,
                         Multiseq *virtualmultiseq,
                         GT_UNUSED Multiseq *querymultiseq,
                         StoreMatch *storematch)
{
  PairUint range;
  Uint endpos;

  /* find the boundaries of current cDNA/EST */
  polya_find_boundaries(virtualmultiseq, storematch->Storeseqnum1, &range);

  if (storematch->Storeseqnum2 == 1 && /* is match of poly(T) head */
      range.uint0 == storematch->Storeposition1) {
    /* The match begins with the first position of the db-sequence
       (modulo VCTORLENGTH) */
    return 1;
  }

  endpos = storematch->Storeposition1 + storematch->Storelength1 - 1;
  if (storematch->Storeseqnum2 == 0 && /* is match of poly(A) tail */
      range.uint1 == endpos) {
    /* The match ends with the last position of the db-sequence
       (modulo VECTORLENGTH) */
    return 1;
  }
  return 0;
}
