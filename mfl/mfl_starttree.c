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
    int try_length;
    mfl_node_t* try_site;
} mfl_try_t;

typedef struct {
    int             sptadd_num_added;
    int             stpadd_num_toadd;
    mfl_nodearray_t stpadd_addedtips;
    mfl_nodearray_t stpadd_tipstoadd;
    int             stpadd_max_hold;
    
    mfl_nodearray_t* stpadd_holdthreads;    // The stepwise addition sites for the last tree in thread i
    int**           stpadd_thread_ids;      // The thread ids for the ancestral topology of the current topology
    
    int             stpadd_shortest_try;
    int             stpadd_longest_try;
    mfl_node_t*     stpadd_newbranch;
    mfl_node_t*     stpadd_lastnewbranch;
    int             stpadd_num_held_old;
    int             stpadd_num_held_new;
    mfl_try_t**     stpadd_newtries;
    mfl_try_t**     stpadd_oldtries;
    
} mfl_stepwise_addition_t;


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
    
    for (i = 0; i < numparts; ++i) {
        // Stuff here.
        mfl_copy_row_from_partition_into_nodedata(
                                                  n->nodet_charstates[i]->nd_prelim_set,
                                                  partset->ptset_partitions[i],
                                                  rownumber
                                                  );
    }
}


void mfl_apply_characters_to_tips(mfl_tree_t* t, mfl_handle_s* handle, mfl_partition_set_t* parts)
{
    int i = 0;
    int numtaxa = handle->n_taxa;
    mfl_nodearray_t nds = t->treet_treenodes;
    
    for (i = 0; i < numtaxa; ++i) {
        // For each partition, copy into the downpass set of the tip
        mfl_copy_from_all_partitions_into_tip_nodedata(nds[i], parts);
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

void mfl_push_tip_to_addseq_array(mfl_node_t* tip, mfl_stepwise_addition_t* sarec)
{
    sarec->stpadd_tipstoadd[sarec->stpadd_num_toadd] = tip;
    ++sarec->stpadd_num_toadd;
}

void mfl_set_random_addition_sequence(mfl_tree_t* t, mfl_stepwise_addition_t* sarec)
{
    int i = 0;
    int numtips = t->treet_num_taxa;
    int* addsequence = mfl_create_default_taxon_array(t->treet_num_taxa);
    mfl_randomise_array(addsequence, numtips);
    mfl_nodearray_t nds = t->treet_treenodes;
    
    for (i = numtips-1; i >= 0; --i) {
        mfl_push_tip_to_addseq_array(nds[addsequence[i]], sarec);
    }
    
    free(addsequence);
}

void mfl_set_addseq_as_is(mfl_tree_t* t, mfl_stepwise_addition_t* sarec)
{
    int *addsequence = mfl_create_default_taxon_array(t->treet_num_taxa);
    int i = 0;
    int numtips = t->treet_num_taxa;
    mfl_nodearray_t nds = t->treet_treenodes;
    
    for (i = numtips-1; i >= 0; --i) {
        mfl_push_tip_to_addseq_array(nds[addsequence[i]], sarec);
    }
    
    free(addsequence);
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

int mfl_compare_holdthreads(mfl_nodearray_t thread1, mfl_nodearray_t thread2, int n)
{
    int i = n;
    
    do {
        if (thread1[i] == thread2[i]) {
            return i;
        }
        --i;
    } while (i);
}

void mfl_destroy_stepwise_addition(mfl_stepwise_addition_t* sarec)
{
    int i = 0;
    
    free(sarec->stpadd_addedtips);
    free(sarec->stpadd_tipstoadd);
    free(sarec);
}

mfl_stepwise_addition_t* mfl_generate_stepwise_addition(mfl_tree_t* t, mfl_handle_s* handle, mfl_searchrec_t* searchrec)
{
    int i = 0;
    int hold = 0;
    
    mfl_stepwise_addition_t* sarec = (mfl_stepwise_addition_t*)mfl_malloc(sizeof(mfl_stepwise_addition_t), 0);
    sarec->stpadd_tipstoadd = (mfl_nodearray_t)mfl_malloc(t->treet_num_taxa * sizeof(mfl_node_t*), 0);
    sarec->stpadd_num_toadd = 0;
    sarec->stpadd_addedtips = (mfl_nodearray_t)mfl_malloc(t->treet_num_taxa * sizeof(mfl_node_t*), 0);
    sarec->sptadd_num_added = 0;
    
    if (!searchrec->sr_num_trees_held_stepwise) {
        searchrec->sr_num_trees_held_stepwise = 1;
        hold = 1;
    }
    else {
        hold = searchrec->sr_num_trees_held_stepwise;
    }
    
    sarec->stpadd_max_hold = hold;
    
    
    
    
    // Set the terminal node addition sequence in the SA record.
    
    switch (handle->addseq_type) {
        case MFL_AST_ASIS:
            mfl_set_addseq_as_is(t, sarec);
            return sarec;
        case MFL_AST_RANDOM:
            mfl_set_random_addition_sequence(t, sarec);
            return sarec;
        case MFL_AST_SIMPLE:
            // mfl_order_array_by_advancement_index
            return sarec;
        // Do something here for cases that aren't yet handled.
        default:
            return sarec;
            break;
    }
}


int mfl_compare_tries_by_length(const void* t1, const void* t2)
{
    mfl_try_t* try1 = *(mfl_try_t**)t1;
    mfl_try_t* try2 = *(mfl_try_t**)t2;
    
    return (int)(try1->try_length - try2->try_length);
}


/*!
 @discussion Rolls back the addition of terminals to the tree by any number of
 steps. The result is a tree denuded of all of the terminals added between the 
 current topology and the number of steps indicated.
 @param sarecord (mfl_stepwise_addition_t*) the stepwise addition record
 @param steps (int) the number of steps in the stepwise addition to roll back by
 */
void mfl_rollback_additions(mfl_stepwise_addition_t* sarecord, int steps)
{
    int i = 0;
    int lasti = sarecord->sptadd_num_added - 1;
    
    for (i = 0; i < steps; ++i) {
        // mfl_disconnect_branch(sarecord->stpadd_addedtips[lasti - i]);
    }
}


void mfl_restore_addition_to_last_step(mfl_stepwise_addition_t* sarecord, int head)
{
    int i = 0;
    
    // Loop over each tip, following the reinsertion prescribed in a record
}


void mfl_reset_stepwise_addition_record_for_new_branch(mfl_node_t* newbranch, mfl_stepwise_addition_t* sarec)
{
    int i = 0;
    
    sarec->stpadd_shortest_try = 0;
    sarec->stpadd_longest_try = 0;
    sarec->stpadd_lastnewbranch = sarec->stpadd_newbranch;
    sarec->stpadd_newbranch = newbranch;
    
    for (i = 0; i < sarec->stpadd_num_held_new; ++i) {
        sarec->stpadd_oldtries[i]->try_site = sarec->stpadd_newtries[i]->try_site;
        sarec->stpadd_oldtries[i]->try_length = sarec->stpadd_newtries[i]->try_length;
        sarec->stpadd_newtries[i]->try_site = NULL;
        sarec->stpadd_newtries[i]->try_length = 0;
    }
    
    sarec->stpadd_num_held_old = sarec->stpadd_num_held_new;
    sarec->stpadd_num_held_new = 0;
}


bool mfl_push_try_to_record(mfl_node_t* tgt, mfl_stepwise_addition_t* sarecord, int length, mfl_searchrec_t* searchrec)
{
    assert(length);
    bool ret = false;
    int i = 0;
    int num_equal = 0;
    int discardrec = 0;
    
    
    
//    if (sarecord->stpadd_num_held_new < sarecord->stpadd_max_hold) {
//        sarecord->stpadd_newtries[sarecord->stpadd_num_held_new]->try_site = tgt;
//        sarecord->stpadd_newtries[sarecord->stpadd_num_held_new]->try_length = length;
//        ++sarecord->stpadd_num_held_new;
//        ret = true;
//    }
//    else {
//        if (sarecord->stpadd_shortest_try < sarecord->stpadd_longest_try) {
//            sarecord->stpadd_newtries[sarecord->stpadd_num_held_new-1]->try_site = tgt;
//            sarecord->stpadd_newtries[sarecord->stpadd_num_held_new-1]->try_length = length;
//            ret = true;
//        }
//        else {
//            // Randomly choose a tree of equal length.
//            
//            // Count the number of equal tries.
//            for (i = 0; i < sarecord->stpadd_num_held_new; ++i) {
//                if (sarecord->stpadd_newtries[i]->try_length == length) {
//                    ++num_equal;
//                }
//            }
//            
//            // Choose one of these to discard
//            discardrec = gsl_rng_uniform_int(searchrec->sr_random_number, (unsigned long)num_equal);
//            
//            // If that includes the new try, just return
//            if (discardrec == i) {
//                return ret = false;
//            }
//            
//            // Otherwise, swap out the old try record info with the new site
//            num_equal = 0;
//            
//            // Find the try record.
//            for (i = 0; i < sarecord->stpadd_num_held_new; ++i) {
//                if (sarecord->stpadd_newtries[i]->try_length == length) {
//                    ++num_equal;
//                    if (num_equal == discardrec) {
//                        sarecord->stpadd_newtries[i]->try_site = tgt;
//                        break;
//                    }
//                }
//            }
//            
//            ret = true;
//        }
//    }
//
//    qsort(sarecord->stpadd_newtries, sarecord->stpadd_num_held_new, sizeof(mfl_try_t*), &mfl_compare_tries_by_length);
//    sarecord->stpadd_longest_try = sarecord->stpadd_newtries[sarecord->stpadd_num_held_new-1]->try_length;
    
    return ret;
}

void mfl_tryall_traversal(mfl_node_t* n, mfl_node_t* newbranch, mfl_stepwise_addition_t* sarecord, mfl_searchrec_t* searchrec)
{
    mfl_node_t* p = NULL;
    
    if (!n->nodet_tip) {
        p = n->nodet_next;
        
        do {
            mfl_tryall_traversal(p->nodet_edge, newbranch, sarecord, searchrec);
            p = p->nodet_next;
        } while (p != n);
    }
    
    // Try the insertion.
    if (newbranch->nodet_index == n->nodet_index + 1) { // Some temporary arbitrary push criterion
        // push the node ref to the list on the current hold thread number
    }
    
    // If it works, store it.
    
    // Move on.
}


void mfl_try_all_insertions(mfl_node_t* newbranch, mfl_tree_t* t, mfl_searchrec_t* searchrec)
{

    
    
}

mfl_node_t* mfl_get_next_terminal_in_addseq(mfl_stepwise_addition_t* sarec)
{
    --sarec->stpadd_num_toadd;
    
    sarec->stpadd_addedtips[sarec->sptadd_num_added] = sarec->stpadd_tipstoadd[sarec->stpadd_num_toadd];
    
    sarec->stpadd_tipstoadd[sarec->stpadd_num_toadd] = NULL;
    
    ++sarec->sptadd_num_added;
    
    return sarec->stpadd_addedtips[sarec->sptadd_num_added-1];
}

mfl_node_t* mfl_generate_starting_trichotomy(mfl_tree_t* t, mfl_stepwise_addition_t* sarec)
{
    mfl_node_t* a = NULL;
    mfl_node_t* l = NULL;
    mfl_node_t* r = NULL;
    mfl_node_t* ringnd = NULL;
    
    a = mfl_get_next_terminal_in_addseq(sarec);
    l = mfl_get_next_terminal_in_addseq(sarec);
    r = mfl_get_next_terminal_in_addseq(sarec);
    
    ringnd = mfl_make_new_n_ary_ring_node(2, t->treet_nodestack);
    
    mfl_join_node_edges(a, ringnd);
    mfl_join_node_edges(l, ringnd->nodet_next);
    mfl_join_node_edges(r, ringnd->nodet_next->nodet_next);
    
    return ringnd;
}

void mfl_append_tip_to_ringnode(mfl_node_t* tip, mfl_tree_t* t)
{
    mfl_node_t* ring = mfl_make_new_n_ary_ring_node(2, t->treet_nodestack);
    
    mfl_join_node_edges(tip, ring->nodet_next);
}

void mfl_setup_tree_for_stepwise_addition(mfl_tree_t* t, mfl_stepwise_addition_t* sarec)
{
    mfl_node_t* startpoint = mfl_generate_starting_trichotomy(t, sarec);
    
    // Attach all remaining tips to a ringnode
    int i = 0;
    int unattached = sarec->stpadd_num_toadd;
    
    // Give all remaining tips a ring base so they can be inserted to the starttree
    for (i = 0; i < unattached; ++i) {
        mfl_append_tip_to_ringnode(sarec->stpadd_tipstoadd[i], t);
    }
    
    t->treet_start = startpoint;
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


// Functions for allocating nodal data
void mfl_setup_nodedata(mfl_node_t* node, mfl_partition_set_t* dataparts, bool bottom)
{
    int i = 0;
    int numparts = dataparts->ptset_n_parts;
    
    node->nodet_num_dat_partitions = dataparts->ptset_n_parts;
    
    node->nodet_charstates = (mfl_nodedata_t**)mfl_malloc(node->nodet_num_dat_partitions * sizeof(mfl_nodedata_t*), 0);
    
    for (i = 0; i < numparts; ++i) {
        node->nodet_charstates[i] = (mfl_nodedata_t*)mfl_malloc(sizeof(mfl_nodedata_t), 0);
        node->nodet_charstates[i]->nd_prelim_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        node->nodet_charstates[i]->nd_subtree_prelim_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        
        if (bottom) {
            node->nodet_charstates[i]->nd_final_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
            node->nodet_charstates[i]->nd_subtree_final_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        }
    }
}

mfl_tree_t* mfl_generate_new_starting_tree(mfl_partition_set_t* dataparts, mfl_handle_s* handle, mfl_searchrec_t* searchrec)
{
    int i = 0;
    bool bottom = false;
    
    mfl_tree_t* t = mfl_alloctree_with_nodes(searchrec->sr_num_taxa_included);
    
    for (i = 0; i < searchrec->sr_num_taxa_included; ++i) {
        bottom = false;
        if (t->treet_treenodes[i]->nodet_tip) { // Just in case this ever gets used to allocate memory for ring nodes.
            bottom = true;
        }
        mfl_setup_nodedata(t->treet_treenodes[i], dataparts, bottom);
    }
    
    mfl_apply_characters_to_tips(t, handle, dataparts);

    return t;
}

mfl_treebuffer_t* mfl_get_start_trees(mfl_partition_set_t* dataparts, mfl_handle_s* handle, mfl_searchrec_t* searchrec)
{
    mfl_treebuffer_t* starttrees;
    
    // Building stuff
    mfl_tree_t* t = mfl_generate_new_starting_tree(dataparts, handle, searchrec);
    
    mfl_stepwise_addition_t * sarec = mfl_generate_stepwise_addition(t, handle, searchrec);
    
    mfl_setup_tree_for_stepwise_addition(t, sarec);
    
    char *showtree = mfl_convert_mfl_tree_t_to_newick(t, false);
    dbg_printf("the starting trichotomy: %s\n", showtree);
    free(showtree);
    
    

    
    return starttrees;
}