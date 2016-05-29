/*
 *  mfl_searchrec.c
 *
 *  THE MORPHY FUNCTION LIBRARY
 *  A library for phylogenetic analysis with emphasis on parsimony and
 *  morphology (but someday other methods)
 *
 *  Copyright (C) 2016  by Martin D. Brazeau, Thomas Guillerme,
 *  and Chris Desjardins
 *
 *  Created by Martin Brazeau and Chris Desjardins on 1/26/12.
 *  Some data structs, routines and ideas derived from:
 *      - The PHYLIP package by Joe Felsenstein
 *          <http://evolution.genetics.washington.edu/phylip.html>
 *      - MrBayes by John Huelsenbeck and Fredrik Ronquist
 *          <http://mrbayes.sourceforge.net/>
 *
 *  Any bugs, errors, inefficiences and general amateurish handling are our own
 *  and most likely the responsibility of MDB. We make no guarantees.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "morphy.h"

/* temporary place for prototypes */
void mfl_initialise_searchrec(mfl_searchrec_t* searchrec, const mfl_handle_s* handle);
void mfl_copy_row_from_partition_into_nodedata(mfl_charstate* target, mfl_datapartition_t* datapart, int row);
void mfl_copy_from_all_partitions_into_node_data(mfl_node_t* n, mfl_partition_set_t* partset);
void mfl_apply_characters_to_tips(mfl_tree_t* t, mfl_handle_s* handle, mfl_partition_set_t* parts);
/**/

void mfl_initialise_searchrec(mfl_searchrec_t* searchrec, const mfl_handle_s* handle)
{
    assert(handle);
    searchrec->sr_abort_rep = false;
    searchrec->sr_abort_swapping = false;
    searchrec->sr_best_length   = MORPHY_UINTMAX;
    searchrec->sr_num_taxa_included = handle->n_taxa;
    searchrec->sr_num_chars_included = handle->n_chars;
    // TODO: Handle included/excluded chars.
    
    // User/default values from the handle
    searchrec->sr_stepwise      = handle->addseq_type;
    searchrec->sr_searchtype    = handle->search_type;
    searchrec->sr_bswaptype     = handle->bswap_type;
    
    if (handle->n_iterations) {
        searchrec->sr_num_reps_stepwise = handle->n_iterations;
    }
    else {
        searchrec->sr_num_reps_stepwise = MORPHY_DEFAULT_ADDITION_SEQUENCE_REPS;
    }
    
    if (handle->maxtrees) {
        searchrec->sr_maxtrees = handle->maxtrees;
    }
    else {
        searchrec->sr_maxtrees = MORPHY_DEFAULT_TREE_LIMIT;
    }
    
    if (!handle->autoincrease) {
        searchrec->sr_increase_treebuffer = MORPHY_DEFAULT_TREEBUFFER_AUTOINCREASE_SWITCH;
    }
    else {
        if (!handle->autoinc_incr) {
            searchrec->sr_autoinc_increment = MORPHY_DEFAULT_TREEBUFFER_AUTOINCREASE_AMOUNT;
        }
    }
    
    /* Values that will have to have been calculated before the search but which
       require some initialised values */
    if (!searchrec->sr_num_partitions) {
        searchrec->sr_num_partitions = 1;
    }
}

mfl_searchrec_t* mfl_create_searchrec(mfl_handle_s* handle)
{
    
    mfl_searchrec_t* newrec = (mfl_searchrec_t*)mfl_malloc(sizeof(mfl_searchrec_t), 0);
    mfl_initialise_searchrec(newrec, handle);
    
    return newrec;
}

bool mfl_search_environment(mfl_handle_s* handle)
{
    bool ret = false;
    
    // Set up random number generator
    
    // Set up the searchrec
    
    // Enter the requested type of search
    
    // Append the resultant data to the handle
    
    return ret;
}
