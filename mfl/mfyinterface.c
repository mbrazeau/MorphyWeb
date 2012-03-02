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

mfl_searchrec * mfl_create_searchrec()
{
    mfl_searchrec *newsearchrec;
    
    return newsearchrec = (mfl_searchrec*)malloc(sizeof(mfl_searchrec));
}

void mfl_set_ntax(mfl_handle_t *mfl_handle, void *param_data)
{
    mfl_handle->n_taxa = (int)(*param_data);
}

void mfl_set_nchar(mfl_handle_t *mfl_handle, void *param_data)
{
    mfl_handle->n_chars = (int)(*param_data);
}

void mfl_set_searchtype(mfl_handle_t *mfl_handle, void *param_data)
{
    switch (*param_data) { //should *param_data be type-cast?
        case 0:
            // Exhaustive
            mfl_handle->search_type = 0;
            break;
        case 1:
            // Branch-and-bound
            mfl_handle->search_type = 1;
            break;
        case 2:
            // Heuristic
            mfl_handle->search_type = 2;
            break;
        default:
            //Return an error
            break;
    }
}

void mfl_set_numiterations(mfl_handle_t *mfl_handle, void *param_data)
{
    mfl_handle->n_iterations = (int)(*param_data);
}

void mfl_set_treelimit(mfl_handle_t *mfl_handle, void *param_data)
{
    mfl_handle->n_treelimit = (int)(*param_data);
}

void mfl_set_branchswap_t(mfl_handle_t *mfl_handle, void *param_data)
{
    switch (*param_data) { //should *param_data be type-cast?
        case 0:
            // TBR
            mfl_handle->bswap_type = 0;
            break;
        case 1:
            // SPR
            mfl_handle->bswap_type = 1;
            break;
        case 2:
            // NNI
            mfl_handle->bswap_type = 2;
            break;
        default:
            //Return an error
            break;
    }
}

void mfl_set_ratchet_status(mfl_handle_t *mfl_handle, void *param_data)
{
    mfl_handle->is_ratchet = (bool)(*param_data);
}

void mfl_attach_inputdata(mfl_handle_t *mfl_handle, void *param_data)
{
    mfl_handle->input_data = (char *)param_data;
}

void mfl_set_addseq_t(mfl_handle_t *mfl_handle, void *param_data)
{
    switch (*param_data) { //should *param_data be type-cast?
        case 0:
            // As-is
            mfl_handle->addseq_type = 0;
            break;
        case 1:
            // Random
            mfl_handle->addseq_type = 1;
            break;
        case 2:
            // Simple
            // Not implemented yet
            //mfl_handle->addseq_type = 2;
            break;
        case 3:
            // Closest
            // Not implemented yet
            // mfl_handle->addseq_type = 3;
            break;

        default:
            //Return an error
            break;
    }
}

void mfl_set_collapse(mfl_handle_t *mfl_handle, void *param_data)
{
    mfl_handle->collapse_nolen = (bool)(*param_data);
}

void mfl_set_collapse_value(mfl_handle_t *mfl_handle, void *param_data)
{
    switch (*param_data) { //should *param_data be type-cast?
        case 0:
            // Max length 0
            mfl_handle->collapse_at = 0;
            break;
        case 1:
            // Min length 0
            mfl_handle->collapse_at = 1;
            break;
        case 2:
            // Equal reconstruction sets
            mfl_handle->collapse_at = 3;
            break;
        default:
            //Return an error
            break;
    }
}

void mfl_set_gapormissing(mfl_handle_t *mfl_handle, void *param_data)
{
    mfl_handle->gap_as_missing = (bool)(*param_data);
}

bool mfl_set_parameter(mfl_handle_t *mfl_handle, mfl_param_t param_type, void *param_data)
{
    switch (param_type) {
        case MFL_PT_NUM_TAX:
            mfl_set_ntax(mfl_handle, param_data);
            return true;
        case MFL_PT_NUM_CHAR:
            mfl_set_nchar(mfl_handle, param_data);
            return true;
        case MFL_PT_SEARCH_TYPE:
            mfl_set_searchtype(mfl_handle, param_data);
            return true;
        case MFL_PT_NUM_ITERATIONS:
            mfl_set_numiterations(mfl_handle, param_data);
            return true;
        case MFL_PT_TREELIMIT:
            mfl_set_treelimit(mfl_handle, param_data);
            return true;
        case MFL_PT_BRANCH_SWAP_TYPE:
            mfl_set_branchswap_t(mfl_handle, param_data);
            return true;
        case MFL_PT_RATCHET_SEARCH:
            mfl_set_ratchet_status(mfl_handle, param_data);
            return true;
        case MFL_PT_INPUT_DATA:
            mfl_attach_inputdata(mfl_handle, param_data);
            return true;
        case MFL_PT_ADD_SEQUENCE_TYPE:
            mfl_set_addseq_t(mfl_handle, param_data);
            return true;
        case MFL_PT_COLLAPSE:
            mfl_set_collapse(mfl_handle, param_data);
            return true;
        case MFL_PT_COLLAP_AT:
            mfl_set_collapse_value(mfl_handle, param_data);
            return true;
        case MFL_PT_GAP:
            mfl_set_gapormissing(mfl_handle, param_data);
            return true;
        default:
            return false;
    }
    
    return success;
}
