/*
 *  mfl_starttree.c
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


// Data needed:
    // Type of addition sequence
    // Number of taxa
    // Character data
    // Is or is not rooted (depends on optimality used)
    // Outgroup specification


typedef struct {
    int as_num_to_hold;
    mfl_node_t* as_newbranch;
    mfl_nodearray_t as_optimal_sites;
} mfl_stepwise_t;


void mfl_copy_row_from_partition_into_nodedata(mfl_charstate* target, mfl_datapartition_t* datapart, int row)
{
    int i = 0;
    int numchars = datapart->part_n_chars_max;
    
    for (i = 0; i < numchars; ++i) {
        target[i] = datapart->part_matrix[i + row * numchars];
    }
}

void mfl_copy_from_all_partitions_into_tip_nodedata(mfl_node_t* n, mfl_partition_set_t* partset)
{
    
    
    int i;
    
    assert(n->nodet_tip);
    int rownumber = n->nodet_index;
    
    int numparts = partset->ptset_n_parts;
    assert(n->nodet_num_dat_partitions == numparts);
    
    mfl_datapartition_t* datapart = NULL;
    mfl_nodedata_t* ndata = NULL;
    
    for (i = 0; i < numparts; ++i) {
        // Stuff here.
    }
}

void mfl_apply_characters_to_tips(mfl_tree_t* t, mfl_handle_s* handle, mfl_partition_set_t* parts)
{
    int i = 0;
    int numtaxa = handle->n_taxa;
    mfl_nodearray_t nds = t->treet_treenodes;
    
    for (i = 0; i < numtaxa; ++i) {
        // For each partition, copy into the downpass set of the tip
        nds[i]->nodet_charstates;
    }
}

void mfl_randomise_array(int* array, int nelems)
{
    
    int i = 0;
    int j = 0;
    int t = 0;
    
    if (nelems > 1) {
        for (i = 0; i < nelems - 1; ++i) {
            j = i + random() / (RAND_MAX / (nelems - i) + 1);
            t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

int* mfl_create_default_taxon_array(int num_taxa)
{
    int i = 0;
    int* taxa = (int*)mfl_malloc(num_taxa * sizeof(int), 0);
    
    for (i = 0; i < num_taxa; ++i) {
        taxa[i] = i;
    }
    
    return taxa;
}



void mfl_calculate_advancement_index(mfl_node_t* t, const mfl_node_t* a)
{
    int i = 0;
    int j = 0;
    int ai = 0;
    
    assert(t->nodet_num_dat_partitions == a->nodet_num_dat_partitions);
    int num_partitions = t->nodet_num_dat_partitions;
    int num_chars = 0;
    
    for (i = 0; i < num_partitions; ++i) {
        num_chars = t->nodet_charstates[i]->nd_n_characters;
        for (j = 0; j < num_chars; ++j) {
            // mfl_distance function.
        }
    }
}

void mfl_order_array_by_advancement_index(int* taxa, mfl_datapartition_t* chardata /*some list of outgroup taxa*/)
{
    // If there's a pre-defined outgroup
    //  Select an outgroup taxon
    //      If there's more than one OG taxon, select one of them arbitrarily
    // Else:
    //  If there's a pre-speficied reference taxon, use that one
    // Else:
    //  Select an arbitrary taxon to be the root (input order).
    
    // Create an array of ints for advancement indices.
    // Set the selected outgroup/ancestor to be the reference taxon
    // For each taxon in the array, compute its advancement index
    // Now, sort the array by advancement index
    // *** Return the array of taxon numbers
}

int* mfl_addition_sequence_generator(mfl_handle_s* handle, mfl_searchrec_t* searchrec)
{
    int* addsequence = NULL;
    
    if (searchrec->sr_included_taxa) {
        addsequence = (int*)mfl_malloc(searchrec->sr_num_taxa_included, 0);
        memcpy(addsequence, searchrec->sr_included_taxa, searchrec->sr_num_taxa_included * sizeof(int));
    }
    else {
        addsequence = mfl_create_default_taxon_array(handle->n_taxa);
    }
    
    switch (handle->addseq_type) {
        case MFL_AST_ASIS:
            return addsequence;
        case MFL_AST_RANDOM:
            mfl_randomise_array(addsequence, searchrec->sr_num_taxa_included);
            return addsequence;
        case MFL_AST_SIMPLE:
            // mfl_order_array_by_advancement_index
            return addsequence;
        // Do something here for cases that aren't yet handled.
        default:
            return addsequence;
            break;
    }
}


void mfl_try_all_insertions(mfl_node_t* newbranch, mfl_node_t* entrynode)
{
    
}


bool mfl_generate_starting_trichotomy(mfl_tree_t* t, int* taxon_addition_sequence)
{
    
}

bool mfl_setup_outgroup(mfl_tree_t* t, int* outgroup_taxa, int num_outgroup_taxa)
{
    int i = 0;
    
    t->treet_outgroup_tips = (mfl_nodearray_t)mfl_malloc(num_outgroup_taxa * sizeof(mfl_node_t*), 0);
    
    if (!t->treet_outgroup_tips) {
        return false;
    }
    
    t->treet_num_og_tips = 0;
    for (i = 0; i < num_outgroup_taxa; ++i) {
        t->treet_outgroup_tips[i] = t->treet_treenodes[outgroup_taxa[i]];
        ++t->treet_num_og_tips;
    }
    
    if (t->treet_num_og_tips) {
        return true;
    }
    else {
        return false;
    }
}

mfl_tree_t* mfl_get_start_tree(mfl_partition_set_t* dataparts, mfl_handle_s* handle, mfl_searchrec_t* searchrec)
{
    int* taxon_addition_seq = mfl_addition_sequence_generator(handle, searchrec);
    mfl_tree_t* starttree = mfl_alloctree_with_nodes(handle->n_taxa);
    
    // Create a three-taxon tree (either according to add seq or honouring ingroup/outgroup partition)
    
    // Begin adding branches to this tree, one at a time, trying all positions
    
    // Clean up resources:
    //  Free int arrays
    //  Free nodesets
    
    return starttree;
}