/*
 *  mfyinterface.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 3/1/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include "morphy.h"
#include <sstream>

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
    mfl_struct->search_type = (mfl_search_t)(long int)(param_data);
    if (!((mfl_struct->search_type >= 0) && (mfl_struct->search_type < MFL_ST_MAX)))
    {
        throw mfl_exception(mfl_s2t(mfl_struct), "Invalid search type");
    }
    return true;
}

bool mfl_set_numiterations(mfl_handle_s *mfl_struct, void *param_data)
{
    /* what is the valid range for this param? */
    mfl_struct->n_iterations = (long int)(param_data);
    return true;
}

bool mfl_set_treelimit(mfl_handle_s *mfl_struct, void *param_data)
{
    /* what is the valid range for this param? */
    mfl_struct->n_treelimit = (long int)(param_data);
    return true;
}

bool mfl_set_branchswap_t(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_struct->bswap_type = (mfl_branch_swap_t)(long int)(param_data);
    if (!((mfl_struct->bswap_type >= 0) && (mfl_struct->bswap_type < MFL_BST_MAX)))
    {
        throw mfl_exception(mfl_s2t(mfl_struct), "Invalid branch swap type");
    }
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

void* mfl_get_addseq_t(mfl_handle_s *mfl_struct)
{
    void* ret = (void*)mfl_struct->addseq_type;
    return ret;
}

bool mfl_set_addseq_t(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_struct->addseq_type = (mfl_add_sequence_t)(long int)param_data;
    if (!((mfl_struct->addseq_type >= 0) && (mfl_struct->addseq_type < MFL_AST_MAX)))
    {
        throw mfl_exception(mfl_s2t(mfl_struct), "Invalid add seq type");
    }
    return true;
}

bool mfl_set_collapse(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_struct->collapse_nolen = (bool)(param_data);
    return true;
}

bool mfl_set_collapse_value(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_struct->collapse_at = (mfl_set_collapse_at_t)(long int)param_data;
    if (!((mfl_struct->collapse_at >= 0) && (mfl_struct->collapse_at < MFL_SC_MAX)))
    {
        throw mfl_exception(mfl_s2t(mfl_struct), "Invalid collapse at parameter");
    }
    return true;
}

bool mfl_set_gapormissing(mfl_handle_s *mfl_struct, void *param_data)
{
    mfl_struct->gap_as_missing = (mfl_gap_t)(long int)(param_data);
    return true;
}

mfl_handle_t mfl_create_handle()
{
    mfl_handle_s *mfl_struct;
    mfl_struct = (mfl_handle_s*)malloc(sizeof(mfl_handle_s));
    
    memset(mfl_struct, 0, sizeof(mfl_handle_s));

    /* setup reasonable defaults */
    mfl_set_branchswap_t(mfl_struct, (void*)MFL_BST_SPR);
    mfl_set_collapse_value(mfl_struct, (void*)MFL_SC_MAX_LEN);
    mfl_set_treelimit(mfl_struct, (void*)MORPHY_DEFAULT_TREE_LIMIT);
    mfl_set_gapormissing(mfl_struct, (void*)MFL_GAP_INAPPLICABLE);

    return mfl_s2t(mfl_struct);
}

char** mfl_get_saved_trees_newick(mfl_handle_t mfl_handle)
{
    mfl_handle_s *mfl_struct = mfl_t2s(mfl_handle);
    
    return mfl_struct->resultant_data->newicktrees;
}

void mfl_erase_trees_newick(mfl_handle_t mfl_handle)
{
    int i;
    mfl_handle_s *mfl_struct = mfl_t2s(mfl_handle);
    char **nwktrees = mfl_struct->resultant_data->newicktrees;
    
    for (i = 0; nwktrees[i]; ++i) {
        if (nwktrees[i]) {
            free(nwktrees[i]);
        }
    }
    free(nwktrees);
}

void mfl_destroy_resultant_data(mfl_handle_t mfl_handle)
{
    
    mfl_handle_s *mfl_struct = mfl_t2s(mfl_handle);
    if (mfl_struct->resultant_data->newicktrees) {
        mfl_erase_trees_newick(mfl_handle);
    }
    free(mfl_struct->resultant_data);
}

void mfl_destroy_handle(mfl_handle_t mfl_handle)
{
    mfl_handle_s *mfl_struct = mfl_t2s(mfl_handle);
    if (mfl_struct->resultant_data) {
        mfl_destroy_resultant_data(mfl_handle);
    }
    free(mfl_struct);
}

int mfl_get_island_count(mfl_handle_s *mfl_struct)
{
    return mfl_struct->resultant_data->num_islands;
}

int mfl_get_island_size(mfl_handle_s *mfl_struct, int island_number)
{
    return mfl_struct->resultant_data->island_sizes[island_number];
}

int mfl_get_island_length(mfl_handle_s *mfl_struct, int island_number)
{
    return mfl_struct->resultant_data->island_lengths[island_number];
}

long mfl_get_num_rearrangements(mfl_handle_s *mfl_struct)
{
    return mfl_struct->resultant_data->n_rearrangements;
}

int mfl_get_num_savedtrees(mfl_handle_s *mfl_struct)
{
    return mfl_struct->resultant_data->n_savetrees;
}

int mfl_get_shortest_treelen(mfl_handle_s *mfl_struct)
{
    return mfl_struct->resultant_data->bestlength;
}

int mfl_get_resultant_data(mfl_handle_t mfl_handle, mfl_resultant_data_t resultant_data, int param)
{
    int ret;
    mfl_handle_s *mfl_struct = mfl_t2s(mfl_handle);
    switch (resultant_data)
    {
        case MFL_RT_ISLAND_COUNT:
            ret = mfl_get_island_count(mfl_struct);
            break;
        case MFL_RT_ISLAND_SIZE:
            ret = mfl_get_island_size(mfl_struct, param);
            break;
        case MFL_RT_ISLAND_LENGTH:
            ret = mfl_get_island_length(mfl_struct, param);
            break;
        case MFL_RT_NUM_REARRANGMENTS:
            ret = mfl_get_num_rearrangements(mfl_struct);
            break;
        case MFL_RT_NUM_SAVED_TREES:
            ret = mfl_get_num_savedtrees(mfl_struct);
            break;
        case MFL_RT_SHORTEST_TREE_LEN:
            ret = mfl_get_shortest_treelen(mfl_struct);
            break;
        default:
            /*
            ** not all cases are handled yet, if you get this assert, 
            ** then you need to update this code to handle the case
            ** you need
            */
            assert(0);
            break;
    }
    return ret;
}

void* mfl_get_parameter(mfl_handle_t mfl_handle, mfl_param_t param_type)
{
    void *ret = NULL;
    mfl_handle_s *mfl_struct = mfl_t2s(mfl_handle);
    switch (param_type)
    {
        case MFL_PT_ADD_SEQUENCE_TYPE:
            ret = mfl_get_addseq_t(mfl_struct);
            break;

        default:
            /* add other cases as needed... */
            assert(0);
            break;
    }
    return ret;
}

bool mfl_set_parameter(mfl_handle_t mfl_handle, mfl_param_t param_type, void *param_data)
{
    bool ret = false;
    mfl_handle_s *mfl_struct = mfl_t2s(mfl_handle);

    switch (param_type)
    {
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
** See... I just added an exception throw to t2s... how useful :D
** I only had to do that once thanks to these little guys...
*/
mfl_handle_t mfl_s2t(mfl_handle_s *mfl_handle)
{
    return (mfl_handle_t)mfl_handle;
}

mfl_handle_s *mfl_t2s(mfl_handle_t mfl_handle)
{
    if (mfl_handle)
    {
        return (mfl_handle_s*)mfl_handle;
    }
    else
    {
        throw mfl_exception(mfl_handle, "Nexus file not open");
        return NULL;
    }
}

mfl_exception::mfl_exception(mfl_handle_t mfl_handle, string const& msg) :
    runtime_error(msg)
{
    m_mfl_handle = mfl_handle;
}

const char * mfl_exception::what() const throw()
{
    string ret = runtime_error::what();
    if (m_mfl_handle)
    {
#ifdef MFY_DEBUG
        mfl_handle_s *mfl_struct = mfl_t2s(m_mfl_handle);
        stringstream ss;
        ss<<endl<<"n_taxa:          "<<mfl_struct->n_taxa;
        ss<<endl<<"n_chars:         "<<mfl_struct->n_chars;
        ss<<endl<<"search_type:     "<<mfl_struct->search_type;
        ss<<endl<<"n_iterations:    "<<mfl_struct->n_iterations;
        ss<<endl<<"n_treelimit:     "<<mfl_struct->n_treelimit;
        ss<<endl<<"bswap_type:      "<<mfl_struct->bswap_type;
        ss<<endl<<"is_ratchet:      "<<mfl_struct->is_ratchet;
        ss<<endl<<"addseq_type:     "<<mfl_struct->addseq_type;
        ss<<endl<<"collapse_nolen:  "<<mfl_struct->collapse_nolen;
        ss<<endl<<"collapse_at:     "<<mfl_struct->collapse_at;
        ss<<endl<<"gap_as_missing:  "<<mfl_struct->gap_as_missing;
        ss<<endl;
        ret.append(ss.str());
#endif
    }
    else
    {
        ret.append(" - Morphy handle is NULL");
    }
    return ret.c_str();
}


