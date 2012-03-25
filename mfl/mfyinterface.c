/*
 *  mfyinterface.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 3/1/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include "morphy.h"

bool mfl_heuristic           (mfl_handle_t mfl_handle)
{
    mfl_handle_s *mfl_struct = mfl_t2s(mfl_handle);

    return mfl_heuristic_search(mfl_struct);
}

void mfl_free_input_data(mfl_handle_s *mfl_struct)
{
    if (mfl_struct->input_data)
    {
        free(mfl_struct->input_data);
        mfl_struct->input_data = NULL;
    }
}

bool mfl_set_ntax(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_struct->n_taxa = (long int)(param_data);
    return true;
}

bool mfl_set_nchar(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_struct->n_chars = (long int)(param_data);
    return true;
}

bool mfl_set_searchtype(mfl_handle_s *mfl_struct, void *param_data)
{
    /* might want to put range checks in here... with error return... */
    mfl_struct->search_type = (mfl_search_t)(long int)(param_data);
    return true;
}

bool mfl_set_numiterations(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_struct->n_iterations = (long int)(param_data);
    return true;
}

bool mfl_set_treelimit(mfl_handle_s *mfl_struct, void *param_data)
{
    if ((long int)(param_data)) 
    {
        mfl_struct->n_treelimit = (long int)(param_data);
    }
    else 
    {
        mfl_struct->n_treelimit = MORPHY_DEFAULT_TREE_LIMIT;
    }

    return true;
}

bool mfl_set_branchswap_t(mfl_handle_s *mfl_struct, void *param_data)
{
    /* might want to put range checks in here... with error return... */
    mfl_struct->bswap_type = (mfl_branch_swap_t)(long int)(param_data);
    return true;
}

bool mfl_set_ratchet_status(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_struct->is_ratchet = (bool)(param_data);
    return true;
}

bool mfl_attach_inputdata(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_free_input_data(mfl_struct);
    mfl_struct->input_data = strdup((char*)param_data);
    return true;
}

bool mfl_set_addseq_t(mfl_handle_s *mfl_struct, void *param_data)
{
    /* might want to put range checks in here... with error return... */
    mfl_struct->addseq_type = (mfl_add_sequence_t)(long int)param_data;
    return true;
}

bool mfl_set_collapse(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_struct->collapse_nolen = (bool)(param_data);
    return true;
}

bool mfl_set_collapse_value(mfl_handle_s *mfl_struct, void *param_data)
{
    /* might want to put range checks in here... with error return... */
    mfl_struct->collapse_at = (mfl_set_collapse_at_t)(long int)param_data;
    return true;
}

bool mfl_set_gapormissing(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_struct->gap_as_missing = (bool)(param_data);
    return true;
}

mfl_handle_t mfl_create_handle()
{
    mfl_handle_s *mfl_struct;
    mfl_struct = (mfl_handle_s*)malloc(sizeof(mfl_handle_s));
    
    memset(mfl_struct, 0, sizeof(mfl_handle_s));

    /* setup reasonable defaults */
    mfl_set_branchswap_t(mfl_struct, (void*)MFL_BST_SPR);

    return mfl_s2t(mfl_struct);
}

void mfl_destroy_handle(mfl_handle_t mfl_handle)
{
    mfl_handle_s *mfl_struct = mfl_t2s(mfl_handle);
    mfl_free_input_data(mfl_struct);
    free(mfl_struct);
}

bool mfl_set_parameter(mfl_handle_t mfl_handle, mfl_param_t param_type, void *param_data)
{
    bool ret = false;
    mfl_handle_s *mfl_struct = mfl_t2s(mfl_handle);

    switch (param_type) {
        case MFL_PT_NUM_TAX:
            ret = mfl_set_ntax(mfl_struct, param_data);
            break;

        case MFL_PT_NUM_CHAR:
            ret = mfl_set_nchar(mfl_struct, param_data);
            break;

        case MFL_PT_SEARCH_TYPE:
            ret = mfl_set_searchtype(mfl_struct, param_data);
            break;

        case MFL_PT_NUM_ITERATIONS:
            ret = mfl_set_numiterations(mfl_struct, param_data);
            break;

        case MFL_PT_TREELIMIT:
            ret = mfl_set_treelimit(mfl_struct, param_data);
            break;

        case MFL_PT_BRANCH_SWAP_TYPE:
            ret = mfl_set_branchswap_t(mfl_struct, param_data);
            break;

        case MFL_PT_RATCHET_SEARCH:
            ret = mfl_set_ratchet_status(mfl_struct, param_data);
            break;

        case MFL_PT_INPUT_DATA:
            ret = mfl_attach_inputdata(mfl_struct, param_data);
            break;

        case MFL_PT_ADD_SEQUENCE_TYPE:
            ret = mfl_set_addseq_t(mfl_struct, param_data);
            break;

        case MFL_PT_COLLAPSE:
            ret = mfl_set_collapse(mfl_struct, param_data);
            break;

        case MFL_PT_COLLAP_AT:
            ret = mfl_set_collapse_value(mfl_struct, param_data);
            break;

        case MFL_PT_GAP:
            ret = mfl_set_gapormissing(mfl_struct, param_data);
            break;

        default:
            break;
    }
    
    return ret;
}

/*
** These may seem dumb now, but may be useful one day
*/
mfl_handle_t mfl_s2t(mfl_handle_s *mfl_handle)
{
    return (mfl_handle_t)mfl_handle;
}

mfl_handle_s *mfl_t2s(mfl_handle_t mfl_handle)
{
    return (mfl_handle_s*)mfl_handle;
}

