//
//  mfl_bitset.h
//  Morphy
//
//  Created by mbrazeau on 27/04/2016.
//
//

#ifndef mfl_bitset_h
#define mfl_bitset_h



#include <stdio.h>
#include "morphy.h"

typedef struct mfl_bitset_t {
    int bts_nfields;
    uint64_t* bts_bitfields;
} mfl_bitset_t;


/**/
bool            mfl_bts_AND(mfl_bitset_t* set1, mfl_bitset_t* set2, mfl_bitset_t* target);
bool            mfl_bts_OR(mfl_bitset_t* set1, mfl_bitset_t* set2, mfl_bitset_t* target);
bool            mfl_bts_XOR(mfl_bitset_t* set1, mfl_bitset_t* set2, mfl_bitset_t* target);
bool            mfl_bts_COMPLEMENT(mfl_bitset_t* bitset, mfl_bitset_t* target);
int             mfl_bts_calculate_n_bitfieds(int n_minbits);
mfl_bitset_t*   mfl_bts_create_bitset(int n_minbits);
bool            mfl_bts_destroy_bitset(mfl_bitset_t* oldbts);

#endif /* mfl_bitset_h */
