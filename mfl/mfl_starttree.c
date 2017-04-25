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
        mfl_copy_row_from_partition_into_nodedata(n->nodet_charstates[i]->nd_prelim_set,
                                                  partset->ptset_partitions[i],
                                                  rownumber);
        memcpy(n->nodet_charstates[i]->nd_initprelim, n->nodet_charstates[i]->nd_prelim_set, n->nodet_charstates[i]->nd_n_characters * sizeof(mfl_charstate));
        memcpy(n->nodet_charstates[i]->nd_initfinal, n->nodet_charstates[i]->nd_prelim_set, n->nodet_charstates[i]->nd_n_characters * sizeof(mfl_charstate));
        memcpy(n->nodet_charstates[i]->nd_final_set, n->nodet_charstates[i]->nd_prelim_set, n->nodet_charstates[i]->nd_n_characters * sizeof(mfl_charstate));
        
        
    }
}


void mfl_apply_characters_to_tips(mfl_tree_t* t, mfl_partition_set_t* parts)
{
    int i = 0;
    int numtaxa = t->treet_num_taxa;
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
    //int ai = 0;
    
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
    
    // TODO: This is temporary
    return -1;
}


void mfl_destroy_stepwise_addition(mfl_stepwise_addition_t* sarec)
{
    
    free(sarec->stpadd_addedtips);
    free(sarec->stpadd_tipstoadd);
    //free(sarec->stpadd_heldtrees);
    // TODO: Add more free calls
    // TODO: Destroy the oldtries record
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
        // TODO: set these as default macros
        searchrec->sr_num_trees_held_stepwise = 1;
        hold = 1;
    }
    else {
        hold = searchrec->sr_num_trees_held_stepwise;
    }
    
    sarec->stpadd_max_hold = hold;
    sarec->stpadd_num_held = 0;
    sarec->stpadd_random_number = searchrec->sr_random_number;
    
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
    
    return try1->try_length - try2->try_length;
}

void mfl_reset_try_buffers(mfl_stepwise_addition_t *sarec)
{
    int i = 0;
    mfl_treebuffer_t *tbfa = sarec->stpadd_newtries;
    mfl_treebuffer_t *tbfb = sarec->stpadd_oldtries;
    
    sarec->stpadd_oldtries = tbfa;
    sarec->stpadd_newtries = tbfb;
    
    for (i = 0; i < sarec->stpadd_num_held; ++i) {
        sarec->stpadd_newtries->tb_savedtrees[i]->treet_parsimonylength = 0;
    }
    
    sarec->stpadd_num_saved = sarec->stpadd_num_held;
}

void mfl_reset_stepwise_addition_record_for_new_branch(mfl_node_t* newbranch, mfl_stepwise_addition_t* sarec)
{
    sarec->stpadd_longest_try   = 0;
    sarec->stpadd_shortest_try  = 0;
    sarec->stpadd_current_try   = 0;
    
    mfl_reset_try_buffers(sarec);
    
    sarec->stpadd_num_held      = 0;
}


int mfl_longest_tree_length(mfl_treebuffer_t* trbuf)
{
    int i = 0;
    int longest = 0;
    
    for (i = 0; i < trbuf->tb_num_trees; ++i) {
        if (trbuf->tb_savedtrees[i]->treet_parsimonylength > longest) {
            longest = trbuf->tb_savedtrees[i]->treet_parsimonylength;
        }
    }
    
    return longest;
}


void mfl_attempt_replacement_stepadd(mfl_treebuffer_t* heldtrs, mfl_tree_t* t, mfl_stepwise_addition_t *sarec)
{
    int i = 0;
    int equals = 0;
    int longer = 0;
    unsigned long breakindex;
    mfl_tree_t* equaltrs[heldtrs->tb_num_trees];
    mfl_tree_t* longest[heldtrs->tb_num_trees];;
    
    if (sarec->stpadd_max_hold == 1) {
        breakindex = gsl_rng_uniform_int(sarec->stpadd_random_number, sarec->stpadd_max_hold + 1);
        
        if (breakindex < sarec->stpadd_max_hold) {
            mfl_update_stored_topology(t, heldtrs->tb_savedtrees[breakindex]);
        }
        
        return;
    }
    
    for (i = 0; i < heldtrs->tb_num_trees; ++i) {
        if (heldtrs->tb_savedtrees[i]->treet_parsimonylength == t->treet_parsimonylength) {
            equaltrs[equals] = heldtrs->tb_savedtrees[i];
            ++equals;
        }
        else if (heldtrs->tb_savedtrees[i]->treet_parsimonylength > t->treet_parsimonylength) {
            longest[longer] = heldtrs->tb_savedtrees[i];
            ++longer;
        }
    }
    
    if (longer) {
        
        int max_len = 0;
        
        for (i = 0; i < longer; ++i) {
            if (longest[i]->treet_parsimonylength > max_len) {
                max_len = longest[i]->treet_parsimonylength;
                breakindex = i;
            }
        }
        sarec->stpadd_longest_try = mfl_longest_tree_length(heldtrs);
        mfl_update_stored_topology(t, longest[breakindex]);
    }
    else {
        
        breakindex = gsl_rng_uniform_int(sarec->stpadd_random_number, equals);
        
        if (breakindex < equals + 1) {
            // Replace tree at breakindex
            mfl_update_stored_topology(t, equaltrs[breakindex]);
            
        }
    }
    
    sarec->stpadd_longest_try = mfl_longest_tree_length(heldtrs);
    
}

mfl_partition_set_t *testparts;

int mfl_test_expected_length(mfl_node_t* src, mfl_node_t* tgt, mfl_tree_t* t, mfl_searchrec_t* search)
{
    int origlength = t->treet_parsimonylength;
    t->treet_parsimonylength = 0;
    mfl_cliprec_t clip;
    char *tree = NULL;
    
    mfl_insert_branch_with_ring_base(src, tgt);
    src->nodet_edge->nodet_weight = 3;
    
    mfl_fullpass_tree_optimisation(t, testparts);
    
    if (origlength != t->treet_parsimonylength) {
        dbg_pfail("lengths unequal: ");
        dbg_printf("Indirect: %i; direct: %i\n", origlength, t->treet_parsimonylength);
        dbg_printf("Tip added: %i\n", src->nodet_tip);
        tree = mfl_convert_mfl_tree_t_to_newick(t, false);
        dbg_printf("On this topology: %s\n", tree);
        dbg_printf("At node: %i\n\n", tgt->nodet_index);
        free(tree);
    }
    
    mfl_clip_branch(src, &clip);
    
    mfl_fullpass_tree_optimisation(t, testparts);
//    t->treet_parsimonylength = origlength;
    
    return 1;
}

void mfl_tryall_traversal(mfl_node_t* n, mfl_node_t* newbranch, mfl_stepwise_addition_t* sarecord, mfl_searchrec_t* searchrec)
{
    int cost = 0;
    int length = 0;
    mfl_node_t* p = NULL;
    mfl_cliprec_t clip;
    
    if (!n->nodet_tip) {
        p = n->nodet_next;
        
        do {
            mfl_tryall_traversal(p->nodet_edge, newbranch, sarecord, searchrec);
            p = p->nodet_next;
        } while (p != n);
    }
    
    // TODO: You can probably optimise by passing the length in instead of -1

    length = searchrec->sr_swaping_on->treet_parsimonylength;
    
    mfl_local_add_cost(newbranch, n->nodet_edge, NULL, &cost);
    searchrec->sr_swaping_on->treet_parsimonylength += cost;
    
//    mfl_test_expected_length(newbranch, n, searchrec->sr_swaping_on, searchrec);
    
    if (sarecord->stpadd_num_held < sarecord->stpadd_max_hold) {
        mfl_insert_branch_with_ring_base(newbranch, n);
        newbranch->nodet_edge->nodet_weight = 3;
        // Push try to record.
        mfl_hold_new_tree(searchrec->sr_swaping_on, sarecord);
        // Update best/worst tries.
        mfl_clip_branch(newbranch, &clip);
        sarecord->stpadd_longest_try = mfl_longest_tree_length(sarecord->stpadd_newtries);
    }
    else {
        if (searchrec->sr_swaping_on->treet_parsimonylength < sarecord->stpadd_longest_try) {
            
            mfl_insert_branch_with_ring_base(newbranch, n);
            newbranch->nodet_edge->nodet_weight = 3;

            mfl_attempt_replacement_stepadd(sarecord->stpadd_newtries, searchrec->sr_swaping_on, sarecord);
            
            mfl_clip_branch(newbranch, &clip);
        }
    }
    
    searchrec->sr_swaping_on->treet_parsimonylength = length;
}


mfl_node_t* mfl_get_next_terminal_in_addseq(mfl_stepwise_addition_t* sarec)
{
    mfl_node_t* ret = NULL;
    
    if (sarec->stpadd_num_toadd) {
        
        --sarec->stpadd_num_toadd;
        
        if (sarec->stpadd_tipstoadd[sarec->stpadd_num_toadd]) {
            sarec->stpadd_addedtips[sarec->sptadd_num_added] = sarec->stpadd_tipstoadd[sarec->stpadd_num_toadd];
            sarec->stpadd_tipstoadd[sarec->stpadd_num_toadd] = NULL;
        }
#ifdef MFY_DEBUG
        else {
            dbg_eprintf("Attempt to retrieve tip from NULL pointer\n");
        }
#endif

        ++sarec->sptadd_num_added;
        ret = sarec->stpadd_addedtips[sarec->sptadd_num_added-1];
    }
#ifdef MFY_DEBUG
    else {
        //dbg_eprintf("Attempt to add negative tip number\n");
        dbg_printf("Warning: attempt to add negative tip number\n");
    }
#endif
    
    return ret;
}


void mfl_setup_nodedata(mfl_node_t* node, mfl_partition_set_t* dataparts, bool bottom)
{
    int i = 0;
    int numparts = dataparts->ptset_n_parts;
    
    node->nodet_num_dat_partitions = dataparts->ptset_n_parts;
    
    node->nodet_charstates = (mfl_nodedata_t**)mfl_malloc(node->nodet_num_dat_partitions * sizeof(mfl_nodedata_t*), 0);
    
    for (i = 0; i < numparts; ++i) {
        node->nodet_charstates[i] = (mfl_nodedata_t*)mfl_malloc(sizeof(mfl_nodedata_t), 0);
        
        node->nodet_charstates[i]->nd_n_characters = dataparts->ptset_partitions[i]->part_n_chars_included;
        node->nodet_charstates[i]->nd_completed = false;
        
        node->nodet_charstates[i]->nd_prelim_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        node->nodet_charstates[i]->nd_initprelim = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        node->nodet_charstates[i]->nd_prelim2_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        node->nodet_charstates[i]->nd_prelim3_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        node->nodet_charstates[i]->nd_subtree_prelim_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        node->nodet_charstates[i]->nd_subtree_activestates = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        node->nodet_charstates[i]->nd_region_activestates = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        
        node->nodet_charstates[i]->nd_parent_partition = dataparts->ptset_partitions[i];
        node->nodet_charstates[i]->nd_downpass_full = dataparts->ptset_partitions[i]->part_downpass_full;
        node->nodet_charstates[i]->nd_uppass_full = dataparts->ptset_partitions[i]->part_uppass_full;
        node->nodet_charstates[i]->nd_NAdownpass_full = dataparts->ptset_partitions[i]->part_NAdownpass_full;
        node->nodet_charstates[i]->nd_NAuppass_full = dataparts->ptset_partitions[i]->part_NAuppass_full;
        node->nodet_charstates[i]->nd_local = dataparts->ptset_partitions[i]->part_local;
                
        if (bottom) {
            node->nodet_charstates[i]->nd_final_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
            node->nodet_charstates[i]->nd_initfinal = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
            node->nodet_charstates[i]->nd_subtree_final_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        }
    }
}


void mfl_connect_uppass_sets(mfl_node_t* ringmmber, const mfl_node_t* ringbase, int num_parts)
{
    int i = 0;
    
    for (i = 0; i < num_parts; ++i) {
        ringmmber->nodet_charstates[i]->nd_final_set = ringbase->nodet_charstates[i]->nd_final_set;
        ringmmber->nodet_charstates[i]->nd_initfinal = ringbase->nodet_charstates[i]->nd_initfinal;
        ringmmber->nodet_charstates[i]->nd_subtree_final_set = ringbase->nodet_charstates[i]->nd_subtree_final_set;
    }
    
}


void mfl_setup_ringnode_data(mfl_node_t* ringentry, mfl_partition_set_t* dataparts)
{
    mfl_node_t* p = NULL;
    
    p = ringentry;
    mfl_setup_nodedata(p, dataparts, true);
    
    p = p->nodet_next;
    do {
        mfl_setup_nodedata(p, dataparts, false);
        mfl_connect_uppass_sets(p, ringentry, dataparts->ptset_n_parts);
        p = p->nodet_next;
    } while (p != ringentry);
}


mfl_node_t* mfl_make_n_ary_ring_with_nodedata(int num_branches, mfl_nodestack_t* ndstk, mfl_partition_set_t* dataparts)
{
    
    mfl_node_t* ring = mfl_make_new_n_ary_ring_node(num_branches, ndstk);

    mfl_setup_ringnode_data(ring, dataparts);
    
    return ring;
}


void mfl_setup_internal_nodedata_in_assembled_tree(mfl_node_t* n, mfl_partition_set_t* dataparts)
{
    mfl_node_t* p = NULL;
    
    if (n->nodet_tip) {
        return;
    }
    
    p = n->nodet_next;
    do {
        mfl_setup_internal_nodedata_in_assembled_tree(p->nodet_edge, dataparts);
        p = p->nodet_next;
    } while (p != n);
    
    mfl_setup_ringnode_data(n, dataparts);
    
}

// TODO: Incorporate this function (conditionally) into any rerooting function
void mfl_setup_starttree_root(mfl_tree_t* t, mfl_partition_set_t* dataparts)
{
    assert(t->treet_root);
    int i = 0;
    
    mfl_join_node_edges(&t->treet_dummynode, t->treet_root);
    t->treet_dummynode.nodet_num_dat_partitions = dataparts->ptset_n_parts;
    t->treet_dummynode.nodet_charstates = (mfl_nodedata_t**)mfl_malloc(dataparts->ptset_n_parts * sizeof(mfl_nodedata_t*), 0);
    
    for (i = 0; i < dataparts->ptset_n_parts; ++i) {
        
        t->treet_dummynode.nodet_charstates[i] = (mfl_nodedata_t*)mfl_malloc(sizeof(mfl_nodedata_t), 0);
        
        t->treet_dummynode.nodet_charstates[i]->nd_final_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        
        t->treet_dummynode.nodet_charstates[i]->nd_prelim2_set = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        
        t->treet_dummynode.nodet_charstates[i]->nd_initprelim = t->treet_dummynode.nodet_charstates[i]->nd_final_set;
        t->treet_dummynode.nodet_charstates[i]->nd_initfinal = t->treet_dummynode.nodet_charstates[i]->nd_final_set;
        
        t->treet_dummynode.nodet_charstates[i]->nd_prelim_set = t->treet_dummynode.nodet_charstates[i]->nd_final_set;
        
        t->treet_dummynode.nodet_charstates[i]->nd_subtree_activestates = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        
        t->treet_dummynode.nodet_charstates[i]->nd_region_activestates = (mfl_charstate*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
    }
}


void mfl_setup_input_tree_with_node_data(mfl_tree_t* t, mfl_partition_set_t* dataparts)
{
    int i = 0;
    bool bottom = false;
    int num_taxa = t->treet_num_taxa;
    mfl_node_t* entry;
    
    for (i = 0; i < num_taxa; ++i) {
        bottom = false;
        if (t->treet_treenodes[i]->nodet_tip) { // Just in case this ever gets used to allocate memory for ring nodes.
            bottom = true;
        }
        mfl_setup_nodedata(t->treet_treenodes[i], dataparts, bottom);
    }
    
    if (t->treet_root) {
        entry = t->treet_root;
    }
    else {
        entry = t->treet_start;
    }
    
    mfl_setup_internal_nodedata_in_assembled_tree(entry, dataparts);
    
    // TODO: Conditional rooting business
    
    
    // TODO: Centralise this process
    mfl_setup_starttree_root(t, dataparts);
    
    mfl_apply_characters_to_tips(t, dataparts);
}


mfl_node_t* mfl_generate_starting_trichotomy(mfl_tree_t* t, mfl_stepwise_addition_t* sarec, mfl_partition_set_t* dataparts)
{
    mfl_node_t* a = NULL;
    mfl_node_t* l = NULL;
    mfl_node_t* r = NULL;
    mfl_node_t* ringnd = NULL;
    
    a = mfl_get_next_terminal_in_addseq(sarec);
    l = mfl_get_next_terminal_in_addseq(sarec);
    r = mfl_get_next_terminal_in_addseq(sarec);
    
    ringnd = mfl_make_n_ary_ring_with_nodedata(2, t->treet_nodestack, dataparts);//mfl_make_new_n_ary_ring_node(2, t->treet_nodestack);
    
    mfl_join_node_edges(a, ringnd);
    mfl_join_node_edges(l, ringnd->nodet_next);
    mfl_join_node_edges(r, ringnd->nodet_next->nodet_next);
    
    t->treet_edges = (mfl_nodearray_t)mfl_malloc(t->treet_num_nodes * sizeof(mfl_node_t*), 0);
    //mfl_update_stored_topology(t, t);
    
    return ringnd;
}


void mfl_append_tip_to_ringnode(mfl_node_t* tip, mfl_tree_t* t, mfl_partition_set_t* dataparts)
{
    mfl_node_t* ring = mfl_make_n_ary_ring_with_nodedata(2, t->treet_nodestack, dataparts);//mfl_make_new_n_ary_ring_node(2, t->treet_nodestack);
    
    mfl_join_node_edges(tip, ring->nodet_next);
}



void mfl_setup_tree_for_stepwise_addition(mfl_tree_t* t, mfl_stepwise_addition_t* sarec, mfl_partition_set_t* dataparts)
{
    mfl_node_t* startpoint = mfl_generate_starting_trichotomy(t, sarec, dataparts);
    
    // Attach all remaining tips to a ringnode
    int i = 0;
    int unattached = sarec->stpadd_num_toadd;
    
    // Give all remaining tips a ring base so they can be inserted to the starttree
    for (i = 0; i < unattached; ++i) {
        mfl_append_tip_to_ringnode(sarec->stpadd_tipstoadd[i], t, dataparts);
    }
    
    // Create the root ring
    mfl_root_target_edge(t, startpoint);
    // TODO: Might combine the below
    mfl_setup_ringnode_data(t->treet_root, dataparts);
    mfl_setup_starttree_root(t, dataparts);
    
    // TODO: Decide if tree is rooted/unrooted
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
    
    mfl_apply_characters_to_tips(t, dataparts);

    return t;
}


void mfl_add_tips_stepwise(mfl_tree_t* t, mfl_partition_set_t* dataparts, mfl_searchrec_t* searchrec)
{
    // For every tip:
    //      For each topology in the buffer:
    //          Try all insertions
}

int mfl_hold_new_tree(mfl_tree_t* t, mfl_stepwise_addition_t *sarec)
{
    int ret = -1;
    
    if (sarec->stpadd_num_held <= sarec->stpadd_max_hold) {
        mfl_update_stored_topology(t, sarec->stpadd_newtries->tb_savedtrees[sarec->stpadd_num_held]);
        ++sarec->stpadd_num_held;
        ret = 0;
    }
#ifdef MFY_DEBUG
    else {
        dbg_eprintf("attempt to hold more trees than permitted by user\n");
        ret = 1;
    }
#endif

    return ret;
}


mfl_treebuffer_t* mfl_get_start_trees(mfl_partition_set_t* dataparts, mfl_handle_s* handle, mfl_searchrec_t* searchrec)
{
    int i = 0;
    int num_trees = MORPHY_DEFAULT_TREE_LIMIT;
    mfl_node_t *newbranch;
    
    if (handle->maxtrees) {
        num_trees = handle->maxtrees;
    }
    
    mfl_treebuffer_t* tbf1 = mfl_alloc_treebuffer(handle->n_to_hold);
    mfl_treebuffer_t* tbf2 = mfl_alloc_treebuffer(handle->n_to_hold);

    // Create starting tree
    mfl_tree_t* t = mfl_generate_new_starting_tree(dataparts, handle, searchrec);
    mfl_update_stored_topology(t, t);
    
//    mfl_nodearray_t current = t->treet_treenodes;
    mfl_stepwise_addition_t * sarec = mfl_generate_stepwise_addition(t, handle, searchrec);
    
    sarec->stpadd_oldtries = tbf1;
    sarec->stpadd_newtries = tbf2;
    
    for (i = 0; i < sarec->stpadd_max_hold; ++i) {
        mfl_append_tree_to_treebuffer(mfl_alloc_empty_tree(searchrec->sr_num_taxa_included), tbf1, handle);
        mfl_append_tree_to_treebuffer(mfl_alloc_empty_tree(searchrec->sr_num_taxa_included), tbf2, handle);
    }
    
    mfl_setup_tree_for_stepwise_addition(t, sarec, dataparts);
    
    mfl_fullpass_tree_optimisation(t, dataparts);
    mfl_hold_new_tree(t, sarec);
  
    testparts = dataparts;
    
#ifdef MFY_DEBUG
    char *showtree = mfl_convert_mfl_tree_t_to_newick(t, false);
    dbg_printf("the starting trichotomy: %s\n", showtree);
    free(showtree);
    int fail = 0;
    int oldlen;
#endif
    searchrec->sr_swaping_on = t;
    
    // While(new branches to add)
    while ((newbranch = mfl_get_next_terminal_in_addseq(sarec))) {
        
        dbg_printf("****\nAdding tip %i:\n***\n\n", newbranch->nodet_tip);
        // Reset bestlen thingy if this is the first loop
        mfl_reset_stepwise_addition_record_for_new_branch(newbranch, sarec);
        
        for (i = 0; i < sarec->stpadd_num_saved; ++i) {
            // FOR each tree held in old list {
            dbg_printf("**Hold %i in tip %i***\n", i, newbranch->nodet_tip);
            mfl_convert_from_stored_topol(sarec->stpadd_oldtries->tb_savedtrees[i], t);
            
            // Optimise the tree
            t->treet_parsimonylength = 0;
            mfl_fullpass_tree_optimisation(t, dataparts);
            oldlen = t->treet_parsimonylength;
            if (oldlen != sarec->stpadd_oldtries->tb_savedtrees[i]->treet_parsimonylength) {
                dbg_printf("DOH!\n");
                dbg_printf("Estimated: %i; direct: %i\n", sarec->stpadd_oldtries->tb_savedtrees[i]->treet_parsimonylength, oldlen);
                ++fail;
            }
            // Try all insertions for new branch
            // TODO: This should be performed on an unrooted tree
            mfl_tryall_traversal(t->treet_root->nodet_next->nodet_next->nodet_edge->nodet_next->nodet_next->nodet_edge, newbranch, sarec, searchrec);

        }
    }
    
#ifdef MFY_DEBUG
    if (fail) {
        dbg_pfail("lenth mismatches occurred.");
        dbg_printf("Function failed %i times\n", fail);
    }
    
    for (i = 0; i < sarec->stpadd_newtries->tb_num_trees; ++i) {
        mfl_convert_from_stored_topol(sarec->stpadd_newtries->tb_savedtrees[i], t);
        showtree = mfl_convert_mfl_tree_t_to_newick(t, false);
        dbg_printf("tree morphy_%i = %s\n", i, showtree);
        dbg_printf("[Length: %i]\n", sarec->stpadd_newtries->tb_savedtrees[i]->treet_parsimonylength);
        free(showtree);
    }
    tui_check_broken_tree(t, false);
#endif
    
    tbf1 = sarec->stpadd_newtries;
    mfl_destroy_treebuffer(sarec->stpadd_oldtries, false);
    mfl_destroy_stepwise_addition(sarec);
    
    return tbf1;
}
