//
//  mfl_bitset.c
//  Morphy
//
//  Created by mbrazeau on 27/04/2016.
//
//

#include "mfl_bitset.h"

/*!
 @discussion Sets a specified bit position in the bitset to the specified value 
 (either 1 or 0).
 @param bitset (mfl_bitset_t*) pointer to the bitset being manipulated
 @param set_to (mfl_bitfield_t) either 0 or 1, specifying the desired value at 
 the set position
 @param setposition (int) the bit position to be set.
 */
void mfl_bts_setbit(mfl_bitset_t* bitset, mfl_bitfield_t set_to, int setposition)
{
    int i = 0;
    
    while ((++i * MFL_BTS_IN_BITSET) < setposition);
    --i;
    
    bitset->bts_bitfields[i] |= (set_to << ((setposition - (i * MFL_BTS_IN_BITSET)) - 1));
}


/*!
 @discussion Simulates a bitwise AND on an mfl_bitset_t struct. If a target 
 bitset is supplied, that target will be set to the result of a bitwise AND 
 between the test sets. Otherwise, returns a true/false value
 @param set1 (mfl_bitset_t*) a bitset operand
 @param set2 (mfl_bitset_t*) the second bitset operand
 @param target (mfl_bitset_t*) an optional target bitset to the receive the 
 result of bitwise and between set1 and set2.
 @return true if bitwise AND exists between set1 and set2; otherwise false
 */
bool mfl_bts_AND(mfl_bitset_t* set1, mfl_bitset_t* set2, mfl_bitset_t* target)
{
    bool ret = false;
    int i = 0;
    int n_bitfields = set1->bts_nfields;
    
    if (target) {
        for (i = 0; i < n_bitfields; ++i) {
            if ((target->bts_bitfields[i] = (set1->bts_bitfields[i] & set2->bts_bitfields[i]))) {
                ret = true;
            }
        }
    } else {
        for (i = 0; i < n_bitfields; ++i) {
            if ((set1->bts_bitfields[i] & set2->bts_bitfields[i])) {
                return true;
            }
        }
    }
    
    return ret;
}


/*!
 @discussion Simulates a bitwise OR on an mfl_bitset_t struct. If a target
 bitset is supplied, that target will be set to the result of a bitwise OR
 between the test sets. Otherwise, returns a true/false value.
 @param set1 (mfl_bitset_t*) a bitset operand
 @param set2 (mfl_bitset_t*) the second bitset operand
 @param target (mfl_bitset_t*) an optional target bitset to the receive the
 result of bitwise OR between set1 and set2.
 @return true if bitwise OR exists between set1 and set2; otherwise false
 */
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


/*!
 @discussion Simulates a bitwise XOR on an mfl_bitset_t struct. If a target
 bitset is supplied, that target will be set to the result of a bitwise XOR
 between the test sets. Otherwise, returns a true/false value.
 @param set1 (mfl_bitset_t*) a bitset operand
 @param set2 (mfl_bitset_t*) the second bitset operand
 @param target (mfl_bitset_t*) an optional target bitset to the receive the
 result of bitwise XOR between set1 and set2.
 @return true if bitwise XOR exists between set1 and set2; otherwise false
 */
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


/*!
 @discussion Simulates the bitwise complement operation on an mfl_bitset_t
 @param bitset (mfl_bitset_t*) the source bitset that is operated on
 @param target (mfl_bitset_t*) the target bitset for storing the result
 @return true if the complement is non-null; otherwise false
 */
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


/*!
 @discussion Calculates the number of bitsets of size MFL_BTS_IN_BITSET required
 to accommodate bit fields wider than the machine's maximum integer width.
 @param n_minbits (int) the number of bits required by the bitfield (e.g. the 
 number of taxa in the dataset/trees for creating taxon bipartitions
 @return the number of bitfields required.
 */
int mfl_bts_calculate_n_bitfieds(int n_minbits)
{
    int n_fields = 0;
    
    n_fields = n_minbits / MFL_BTS_IN_BITSET;
    if (n_minbits % MFL_BTS_IN_BITSET) {
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

