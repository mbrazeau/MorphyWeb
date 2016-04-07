/*
 *
 * mfl_characters.c
 *
 * Functions for processing the character input. These functions assign
 * the character types to their relevant partitions as mfl_charstate type
 * vectors.
 *
 */

#include "morphy.h"

/*typedef struct mfl_nodedata_t {
    int dp_n_characters;                        // The number of characters within the datablock.
    int dp_n_taxa;                              // The number of taxa the datablock applies to.
    mfl_optimisation_t dp_optimisation_method;  // The optimisation method applied to all characters in this datablock.
    bool dp_inapplicables;                      // false: no inapplicables; true: has inapplicables.
    mfl_costs_t *dp_costmatrix;                 // Cost matrix associated with these characters.
    mfl_parsim_fn dp_downpass;                  // The downpass parsimony function
    mfl_parsim_fn dp_uppass;                    // The uppass parsimony function
    charstate *dp_prelim_set;                 // The characters to which the datablock applies.
    charstate *dp_final_set;
    charstate *dp_subtree_prelim_set;
    charstate *dp_subtree_final_set;
} mfl_nodedata_t;*/

mfl_nodedata_t* mfl_alloc_datapart(void)
{
    mfl_nodedata_t* newdatapart = NULL;
    
    newdatapart = (mfl_nodedata_t*)malloc(sizeof(mfl_nodedata_t));
    
    if (!newdatapart) {
        dbg_printf("ERROR in mfl_alloc_datapart(): unable to allocate memory for new data partition.\n");
        return NULL;
    }
    else {
        memset(newdatapart, 0, sizeof(mfl_nodedata_t));
    }
    
    return newdatapart;
}


void mfl_free_datapart(mfl_nodedata_t *olddata)
{
    if (olddata->dp_prelim_set) {
        free(olddata->dp_prelim_set);
        olddata->dp_prelim_set = NULL;
    }
    if (olddata->dp_final_set) {
        free(olddata->dp_final_set);
        olddata->dp_final_set = NULL;
    }
    if (olddata->dp_subtree_prelim_set) {
        free(olddata->dp_subtree_prelim_set);
        olddata->dp_subtree_prelim_set = NULL;
    }
    if (olddata->dp_subtree_final_set) {
        free(olddata->dp_subtree_final_set);
        olddata->dp_subtree_final_set = NULL;
    }
    
    // Depending on how cost matrices are handled, they might be freed here, too.
    
    free(olddata);
}


mfl_charstate* mfl_allocate_nodal_character_set(int num_characters)
{
    mfl_charstate *newcharset = NULL;
    
    newcharset = (mfl_charstate*)malloc(num_characters * sizeof(mfl_charstate));
    if (!newcharset) {
        dbg_printf("ERROR in mfl_allocate_nodal_character_set(): unable to allocate memory for character data.\n");
        return NULL;
    }
    else {
        memset(newcharset, 0, num_characters * sizeof(mfl_charstate));
    }
    
    return newcharset;
}
