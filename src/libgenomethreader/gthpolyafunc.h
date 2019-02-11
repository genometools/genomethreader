#ifndef GTHPOLYAFUNC_H
#define GTHPOLYAFUNC_H

#include "match.h"

Sint gthpolyaselectmatch(/*@unused@*/ Alphabet *alpha,
                         Multiseq *virtualmultiseq,
                         /*@unused@*/ Multiseq *querymultiseq,
                         StoreMatch *storematch);

#endif
