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
    int             stpadd_max_hold;
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


int mfl_compare_tries_by_length(const void* t1, const void* t2)
{
    mfl_try_t* try1 = *(mfl_try_t**)t1;
    mfl_try_t* try2 = *(mfl_try_t**)t2;
    
    return (int)(try1->try_length - try2->try_length);
}


void mfl_roll_back_addition(mfl_stepwise_addition_t* sarecord, int position)
{
    
    // Excise last new branch, if necessary
    
    
    
    // Replace branch to its
    //mfl_insert_branch_with_ring_base(sarecord->stpadd_lastnewbranch, sarecord->stpadd_oldtries[position]);
    
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
    
    if (sarecord->stpadd_num_held_new < sarecord->stpadd_max_hold) {
        sarecord->stpadd_newtries[sarecord->stpadd_num_held_new]->try_site = tgt;
        sarecord->stpadd_newtries[sarecord->stpadd_num_held_new]->try_length = length;
        ++sarecord->stpadd_num_held_new;
        ret = true;
    }
    else {
        if (sarecord->stpadd_shortest_try < sarecord->stpadd_longest_try) {
            sarecord->stpadd_newtries[sarecord->stpadd_num_held_new-1]->try_site = tgt;
            sarecord->stpadd_newtries[sarecord->stpadd_num_held_new-1]->try_length = length;
            ret = true;
        }
        else {
            // Randomly choose a tree of equal length.
            
            // Count the number of equal tries.
            for (i = 0; i < sarecord->stpadd_num_held_new; ++i) {
                if (sarecord->stpadd_newtries[i]->try_length == length) {
                    ++num_equal;
                }
            }
            
            // Choose one of these to discard
            discardrec = gsl_rng_uniform_int(searchrec->sr_random_number, (unsigned long)num_equal);
            
            // If that includes the new try, just return
            if (discardrec == i) {
                return ret = false;
            }
            
            // Otherwise, swap out the old try record info with the new site
            num_equal = 0;
            
            // Find the try record.
            for (i = 0; i < sarecord->stpadd_num_held_new; ++i) {
                if (sarecord->stpadd_newtries[i]->try_length == length) {
                    ++num_equal;
                    if (num_equal == discardrec) {
                        sarecord->stpadd_newtries[i]->try_site = tgt;
                        break;
                    }
                }
            }
            
            ret = true;
        }
    }

    qsort(sarecord->stpadd_newtries, sarecord->stpadd_num_held_new, sizeof(mfl_try_t*), &mfl_compare_tries_by_length);
    sarecord->stpadd_longest_try = sarecord->stpadd_newtries[sarecord->stpadd_num_held_new-1]->try_length;
    
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
    
    // Try the insertion
    // Try pushing it to the record
    // Move on
}


void mfl_try_all_insertions(mfl_node_t* newbranch, mfl_tree_t* t, mfl_searchrec_t* searchrec)
{
    mfl_stepwise_addition_t stepaddrec = {0};
    stepaddrec.stpadd_max_hold = searchrec->sr_num_trees_held_stepwise;
    stepaddrec.stpadd_newtries = (mfl_try_t**)mfl_malloc(stepaddrec.stpadd_max_hold * sizeof(mfl_try_t*), 0);
    stepaddrec.stpadd_oldtries = (mfl_try_t**)mfl_malloc(stepaddrec.stpadd_max_hold * sizeof(mfl_try_t*), 0);
    // TODO: Probably need to set up a special try-setup function to allocate the remainig memory
    
    // Do
    
    // For each of the old insertion points;
    //  Restore old branch to that point
    //      Fetch next tip in addition sequence.
    //      Put new tip in the addition sequence record
    //
    
    // While tips are available
    
    free(stepaddrec.stpadd_newtries);
    free(stepaddrec.stpadd_oldtries);
}


mfl_node_t* mfl_generate_starting_trichotomy(mfl_tree_t* t, int* taxon_addition_sequence)
{
    mfl_node_t* a = NULL;
    mfl_node_t* l = NULL;
    mfl_node_t* r = NULL;
    mfl_node_t* ringn = NULL;
    
    // if (there's an outgroup) {
    // Take one outgroup and two ingroup taxa.
    // a = outgroup 1
    // l = ingroup 1
    // r = ingroup 2
    //}
    // else
    // just grab the first three from the addition sequence.
    //{
    a = t->treet_treenodes[taxon_addition_sequence[0]];
    l = t->treet_treenodes[taxon_addition_sequence[1]];
    r = t->treet_treenodes[taxon_addition_sequence[2]];
    //}
    
    ringn = mfl_make_new_n_ary_ring_node(2, t->treet_nodestack);
    
    mfl_join_node_edges(a, ringn);
    mfl_join_node_edges(l, ringn->nodet_next);
    mfl_join_node_edges(r, ringn->nodet_next->nodet_next);
    
    // Do this stuff after the return:
    
    // If there are any directed characters, then root the tree.
    
    // Calculate the length of this trichotomy. (Perhaps after return?)
    
    return ringn;
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


mfl_treebuffer_t* mfl_get_start_trees(mfl_partition_set_t* dataparts, mfl_handle_s* handle, mfl_searchrec_t* searchrec)
{
    mfl_treebuffer_t* holdbuffer1 = mfl_alloc_treebuffer(searchrec->sr_num_trees_held_stepwise);
    mfl_treebuffer_t* holdbuffer2 = mfl_alloc_treebuffer(searchrec->sr_num_trees_held_stepwise);
    
    int* taxon_addition_seq = mfl_addition_sequence_generator(handle, searchrec);
    
    
    // Create a three-taxon tree (either according to add seq or honouring ingroup/outgroup partition)
    
    // Begin adding branches to this tree, one at a time, trying all positions
    
    // Clean up resources:
    //  Free int arrays
    //  Free nodesets
    
    
    return holdbuffer2;
}