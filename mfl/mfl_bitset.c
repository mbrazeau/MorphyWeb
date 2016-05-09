//
//  mfl_bitset.c
//  Morphy
//
//  Created by mbrazeau on 27/04/2016.
//
//

#include "mfl_bitset.h"

void mfl_bts_setbit(mfl_bitset_t* bitset, mfl_bitfield_t set_to, int setposition)
{
    int i = 0;
    
    while ( (++i * MFL_BTS_IN_BITSET) < setposition) {
        printf("i: %i\n", i);
    };
    --i;
    
    bitset->bts_bitfields[i] |= ( set_to << ((setposition - (i * MFL_BTS_IN_BITSET)) - 1) );
}


bool mfl_bts_AND(mfl_bitset_t* set1, mfl_bitset_t* set2, mfl_bitset_t* target)
{
    bool ret = true;
    int i = 0;
    int n_bitfields = set1->bts_nfields;
    
    if (target) {
        for (i = 0; i < n_bitfields; ++i) {
            if (!(target->bts_bitfields[i] = (set1->bts_bitfields[i] & set2->bts_bitfields[i]))) {
                ret = false;
            }
        }
    } else {
        for (i = 0; i < n_bitfields; ++i) {
            if (!(set1->bts_bitfields[i] & set2->bts_bitfields[i])) {
                return false;
            }
        }
    }
    
    return ret;
}


bool mfl_bts_OR(mfl_bitset_t* set1, mfl_bitset_t* set2, mfl_bitset_t* target)
{
    bool ret = false;
    int i = 0;
    int n_bitfields = set1->bts_nfields;
    
    if (target) {
        for (i = 0; i < n_bitfields; ++i) {
            if ((target->bts_bitfields[i] = (set1->bts_bitfields[i] | set2->bts_bitfields[i]))) {
                ret = true;
            }
        }
    } else {
        for (i = 0; i < n_bitfields; ++i) {
            if ((set1->bts_bitfields[i] | set2->bts_bitfields[i])) {
                ret = true;
            }
        }
    }
    
    return ret;
}


bool mfl_bts_XOR(mfl_bitset_t* set1, mfl_bitset_t* set2, mfl_bitset_t* target)
{
    bool ret = false;
    int i = 0;
    int n_bitfields = set1->bts_nfields;
    
    if (target) {
        for (i = 0; i < n_bitfields; ++i) {
            if ((target->bts_bitfields[i] = (set1->bts_bitfields[i] ^ set2->bts_bitfields[i]))) {
                ret = true;
            }
        }
    } else {
        for (i = 0; i < n_bitfields; ++i) {
            if ((set1->bts_bitfields[i] ^ set2->bts_bitfields[i])) {
                ret = true;
            }
        }
    }
    
    return ret;
}


bool mfl_bts_COMPLEMENT(mfl_bitset_t* bitset, mfl_bitset_t* target)
{
    bool ret = false;
    int i = 0;
    int n_bitfields = bitset->bts_nfields;
    
    for (i = 0; i < n_bitfields; ++i) {
        if ((target->bts_bitfields[i] = ~(bitset->bts_bitfields[i]) )) {
            ret = true;
        }
    }

    return ret;
}


int mfl_bts_calculate_n_bitfieds(int n_minbits)
{
    int n_fields = 0;
    
    n_fields = n_minbits / 64;
    if (n_minbits % 64) {
        n_fields += 1;
    }
    
    return n_fields;
}


mfl_bitset_t* mfl_bts_create_bitset(int n_minbits)
{
    int n_fields = 0;
    mfl_bitset_t* newbitset = NULL;
    
    newbitset = (mfl_bitset_t*)malloc(sizeof(mfl_bitset_t));
    if (!newbitset) {
        dbg_eprintf("unable to allocate memory for new bitset");
        return NULL;
    }
    else {
        memset(newbitset, 0, sizeof(mfl_bitset_t));
    }
    
    // Calculate the number of fields required;
    n_fields = mfl_bts_calculate_n_bitfieds(n_minbits);
    
    // Allocate the bitsets
    newbitset->bts_bitfields = (uint64_t*)malloc(n_fields * sizeof(uint64_t));
    if (!newbitset->bts_bitfields) {
        dbg_eprintf("unable to allocate memory for new bitfields in bitset");
        free(newbitset);
        return NULL;
    }
    else {
        memset(newbitset->bts_bitfields, 0, n_fields * sizeof(uint64_t));
    }
    
    newbitset->bts_max_bitfields = n_fields;
    newbitset->bts_nfields = n_fields;
    newbitset->bts_max_bit = n_minbits;
    
    return newbitset;
}


bool mfl_bts_destroy_bitset(mfl_bitset_t* oldbts)
{
    bool ret = false;
    if (oldbts) {
        if (oldbts->bts_bitfields) {
            free(oldbts->bts_bitfields);
        }
        free(oldbts);
        ret = true;
    }
    
    return ret;
}

