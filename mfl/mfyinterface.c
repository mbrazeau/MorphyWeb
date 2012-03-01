/*
 *  mfyinterface.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 3/1/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

mfl_handle_t* mfl_create_handle()
{
    mfl_handle_t *newhandl;
    
    return newhandl = (mfl_handle_t*)malloc(sizeof(mfl_handle_t));
}

void mfl_destroy_handle(mfl_handle_t *mfl_handle)
{
    free(mfl_handle);
}

bool mfl_set_parameter(mfl_handle_t *mfl_handle, mfl_param_t param_type, void *param_data)
{
    
    bool success;
    
    switch (param_type) {
        case MFL_PT_NUM_TAX:
            mfl_set_ntax(mfl_handle, param_data);
            success = true;
            break;
        case MFL_PT_NUM_CHAR:
            mfl_set_nchar(mfl_handle, param_data);
            success = true;
            break;
        case MFL_PT_SEARCH_TYPE:
            mfl_set_searchtype(mfl_handle, param_data);
            success = true;
            break;
        case MFL_PT_NUM_ITERATIONS:
            mfl_set_numiterations(mfl_handle, param_data);
            success = true;
            break;
        case MFL_PT_TREELIMIT:
            mfl_set_treelimit(mfl_handle, param_data);
            success = true;
            break;
        case MFL_PT_BRANCH_SWAP_TYPE:
            mfl_set_branchswap_t(mfl_handle, param_data);
            success = true;
            break;
        case MFL_PT_RATCHET_SEARCH:
            mfl_set_ratchet_status(mfl_handle, param_data);
            success = true;
            break;
        case MFL_PT_INPUT_DATA:
            mfl_attach_inputdata(mfl_handle, param_data);
            success = true;
            break;
        case MFL_PT_ADD_SEQUENCE_TYPE:
            mfl_set_addseq_t(mfl_handle, param_data);
            success = true;
            break;
        case MFL_PT_COLLAPSE:
            mfl_set_collapse(mfl_handle, param_data);
            success = true;
            break;
        case MFL_PT_COLLAP_AT:
            mfl_set_collapse_value(mfl_handle, param_data);
            success = true;
            break;
        default:
            success = false;
            break;
    }
    
    return success;
}
