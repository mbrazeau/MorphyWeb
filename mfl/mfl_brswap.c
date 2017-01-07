/*
 *  mfl_brswap.c
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

#include <stdio.h>
#include "morphy.h"


void mfl_save_topology(mfl_tree_t* t, mfl_treebuffer_t* trbuf, mfl_searchrec_t* searchrec)
{
    mfl_tree_t* topolcopy = mfl_copy_tree_topology(t);
    mfl_append_tree_to_treebuffer(topolcopy, trbuf, searchrec->sr_handle_ptr);
}



bool mfl_clip_branch(mfl_node_t* n, mfl_cliprec_t* cliprec)
{

    if (n->nodet_edge->nodet_weight < 3) {
        return false;
    }
    
    cliprec->src1 = n->nodet_edge->nodet_next;
    cliprec->src2 = cliprec->src1->nodet_next;
    cliprec->tgt1 = cliprec->src1->nodet_edge;
    cliprec->tgt2 = cliprec->src2->nodet_edge;
    
    mfl_disconnect_node(cliprec->src1);
    mfl_disconnect_node(cliprec->src2);
    
    mfl_join_node_edges(cliprec->tgt1, cliprec->tgt2);
    
    return true;
}


bool mfl_restore_branching(mfl_cliprec_t* cliprec)
{
    if (!cliprec->src1) {
        return false;
    }
    if (!cliprec->src2) {
        return false;
    }
    if (!cliprec->tgt1) {
        return false;
    }
    if (!cliprec->tgt2) {
        return false;
    }
    
    mfl_join_node_edges(cliprec->src1, cliprec->tgt1);
    mfl_join_node_edges(cliprec->src2, cliprec->tgt2);
    
    return true;
}


inline void mfl_temp_rebranching(mfl_node_t* src, mfl_node_t* tgt, mfl_cliprec_t* regraft)
{
    regraft->tgt1 = tgt;
    regraft->tgt2 = tgt->nodet_edge;
    regraft->src1 = src;
    
    if (src->nodet_next->nodet_edge) {
        regraft->src2 = src->nodet_next->nodet_next;
    }
    else {
        regraft->src2 = src->nodet_next;
        assert(!regraft->src2->nodet_edge);
    }
    
    mfl_join_node_edges(regraft->src1, regraft->tgt1);
    mfl_join_node_edges(regraft->src2, regraft->tgt2);
}


inline void mfl_undo_temp_rebranching(mfl_cliprec_t* regraft)
{
    regraft->src1->nodet_edge = NULL;
    regraft->src2->nodet_edge = NULL;
    mfl_join_node_edges(regraft->tgt1, regraft->tgt2);
}


mfl_node_t* mfl_get_src_root(mfl_node_t *n)
{
    mfl_node_t* p = n->nodet_next;
    
    while (!p->nodet_edge) {
        p = p->nodet_next;
    }
    
    return p->nodet_edge;
}


void mfl_TBReconnect_traversal(mfl_node_t* src, mfl_node_t* tgt, mfl_searchrec_t* searchrec, bool isbackswap)
{
    mfl_node_t* p = NULL;
    
    // Preorder stuff
    
    if (tgt->nodet_tip) {
        return;
    }
    
    p = tgt->nodet_next;
    
    do {
        mfl_TBReconnect_traversal(src, p->nodet_edge, searchrec, isbackswap);
        p = p->nodet_next;
    } while (p != tgt);
    
    
}


void mfl_initialise_tbr_rerooting_record(int node_weight, mfl_subtree_edges_t* stedges)
{
    if (node_weight < 3) {
        return;
    }
    stedges->ste_max_edges = 2 * node_weight - 4;
    stedges->ste_num_edges = 0;
    stedges->ste_head = 0;
    stedges->ste_edges = (mfl_node_t**)mfl_malloc(stedges->ste_max_edges * sizeof(mfl_node_t*), 0);
}


void mfl_push_edge_ref_to_record(mfl_node_t* edge, mfl_subtree_edges_t* steges)
{
    steges->ste_edges[steges->ste_num_edges] = edge;
    ++steges->ste_num_edges;
    assert(steges->ste_num_edges <= steges->ste_max_edges);
}



int mfl_get_all_subtree_edges(mfl_node_t* n, mfl_subtree_edges_t* steges, mfl_node_t* start)
{
 
    mfl_node_t* p = NULL;
    
    if (n->nodet_tip) {
        return 0;
    }
    
    p = n->nodet_next;
    
    do {
        mfl_get_all_subtree_edges(p->nodet_edge, steges, start);
        if (n != start) {
            mfl_push_edge_ref_to_record(p->nodet_edge, steges);
        }
        p = p->nodet_next;
    } while (p != n);
    
    if (steges->ste_num_edges) {
        return steges->ste_num_edges;
    }
    else {
        return 0;
    }
}



bool mfl_reroot_subtree(mfl_node_t* subtr, mfl_subtree_edges_t* stedges, mfl_cliprec_t* initial)
{
    
    if (subtr->nodet_weight < 3) {
        return false;
    }

        
    if (!stedges->ste_head) {
        mfl_get_all_subtree_edges(subtr, stedges, subtr);
        initial->src1 = subtr->nodet_next;
        initial->src2 = subtr->nodet_next->nodet_next;
        initial->tgt1 = subtr->nodet_next->nodet_edge;
        initial->tgt2 = subtr->nodet_next->nodet_next->nodet_edge;
    }
    
    // Reroot the tree
    
    if (stedges->ste_head < stedges->ste_num_edges) {
        
        // Perform rerooting:
        // Remove the root from current position
        mfl_join_node_edges(subtr->nodet_next->nodet_edge, subtr->nodet_next->nodet_next->nodet_edge);
        
        // Insert the root between new edges
        mfl_join_node_edges(subtr->nodet_next->nodet_next, stedges->ste_edges[stedges->ste_head]->nodet_edge);
        mfl_join_node_edges(subtr->nodet_next, stedges->ste_edges[stedges->ste_head]);
        
        ++stedges->ste_head;
        return true;
    }
    
    // Restore the root to its original position.
    mfl_join_node_edges(initial->src2->nodet_edge, initial->src1->nodet_edge);
    mfl_restore_branching(initial);
    return false;
}



void mfl_regrafting_traversal(mfl_node_t* tgt, mfl_node_t* src, mfl_searchrec_t* searchrec, int startdistance, int trav)
{
    // Traverse the target tree, attempting the reinsertion
    
    mfl_node_t* p = NULL;
    mfl_cliprec_t regraft;

    ++trav;
    
    if (trav > startdistance) {
//        // Insert branch
//        mfl_temp_rebranching(src, tgt, &regraft);
//        
//        // Copy the tree, append it to the buffer
//        mfl_save_topology(searchrec->sr_swaping_on, searchrec->sr_treebuffer, searchrec);
//        
//        // Count the number of rearrangements
//        mfl_undo_temp_rebranching(&regraft);
        
        
        ++searchrec->sr_rearrangement_counter;
    }
    
    if (!tgt) {
        return;
    }
    
    if (tgt->nodet_tip) {
        return;
    }
    
    p = tgt->nodet_next;
    
    if (p->nodet_edge) {
        mfl_regrafting_traversal(p->nodet_edge, src, searchrec, startdistance, trav);
    }
    
    if (p->nodet_edge) {
        mfl_regrafting_traversal(p->nodet_next->nodet_edge, src, searchrec, startdistance, trav);
    }
    
    return;
}


void mfl_regraft_subtree(mfl_node_t* src, mfl_node_t* tgt, mfl_searchrec_t* searchrec, bool isbackswap)
{
    int startdistance = 0;
    int startskip = 1;
    
    if (isbackswap) {
        startdistance = 2;
        startskip = 0;
    }
    
    if (tgt) {
        mfl_regrafting_traversal(tgt, src, searchrec, startdistance, 0);
        mfl_regrafting_traversal(tgt->nodet_edge, src, searchrec, startdistance + startskip, 0);
    }
    else {
        if (!isbackswap) {
            mfl_regrafting_traversal(tgt, src, searchrec, 0, 1);
        }
    }
}


void mfl_break_branch(mfl_node_t* src, mfl_searchrec_t* searchrec)
{
    bool isbackswap = true;
    mfl_cliprec_t srcrootinit = {0};
    mfl_cliprec_t clippoint = {0};
    mfl_subtree_edges_t stedges = {0};
    mfl_initialise_tbr_rerooting_record(src->nodet_weight, &stedges);
    

    if (mfl_clip_branch(src, &clippoint)) {
        isbackswap = true;
    }
    
    do {
        mfl_regraft_subtree(src, clippoint.tgt1, searchrec, isbackswap);
        isbackswap = false;
    }
    while (mfl_reroot_subtree(src, &stedges, &srcrootinit));
    
    mfl_restore_branching(&clippoint);
    
    free(stedges.ste_edges);
}


void mfl_pruning_traversal(mfl_node_t* n, mfl_searchrec_t* searchrec)
{
    mfl_node_t* p = n->nodet_next;
    
    if (!n->nodet_tip) {
        
        do {
            mfl_pruning_traversal(p->nodet_edge, searchrec);
            n->nodet_weight += p->nodet_edge->nodet_weight;
            p = p->nodet_next;
        } while (p != n);
        assert(n->nodet_edge->nodet_weight != searchrec->sr_num_taxa_included);
    }
    
    n->nodet_edge->nodet_weight = searchrec->sr_num_taxa_included - n->nodet_weight;
    
    mfl_break_branch(n, searchrec);
    
}

bool mfl_heuristic_search(mfl_handle_s *mfl_handle)
{
    /* Eventually, this will parse the information in the handle and set up the 
     * heuristic search. For now, it just returns zero because we're not even 
     * close to that yet. */
    return 0;
}
