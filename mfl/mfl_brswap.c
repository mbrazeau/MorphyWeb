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


void mfl_regrafting_traversal(mfl_node_t* n, mfl_node_t* src, mfl_searchrec_t* searchrec)
{
    // Traverse the target tree, attempting the reinsertion
}

// TG: same function but different arguments (we might want to modify it)
void mfl_rebranching_traversal(mfl_treebuffer_t* tree_buffer, mfl_tree_t* rebranch_tree, mfl_tree_t* clip_tree, bool neighbour_rule)
{
    // Traverse the target tree (rebranch_tree), attempting the reinsertion (clip_tree) with or witout neighbour_rule
    
    // We might want to add a counter to know where in the buffer to store the new clipped tree?
}


void mfl_pruning_traversal(mfl_node_t* n, mfl_searchrec_t* searchrec)
{
    // Traverse the starting tree, making the requisit prunings
    
    // At each clipping:
    //      mfl_regrafting_traversal()
}

// TG: algorithmic node, will it be easier to do the whole algorithm as two nested traversal?
mfl_treebuffer_t* algorithm_idea(mfl_tree_t* starttree, mfl_searchrec_t* searchrec, mfl_node_t *start, int* num_taxa)
{
    //Variables
    int clade_counter = 0; //TG: this is the counter for applying the neigbhour rule: any subtree that is <= to clade_counter cannnot rebranch on the neigbhour edges.
    int total_num_trees = 0;
    mfl_tree_t* clipped_tree = NULL;
    mfl_tree_t* branching_tree = NULL;
    mfl_treebuffer_t* tree_buffer = NULL;
    
    // Allocate the right tree buffer size (which is (2*(n-3)*(2*(n-3)-1)) trees)
    total_num_trees = (2*(starttree->treet_num_taxa-3)*(2*(starttree->treet_num_taxa-3)-1));
    tree_buffer = mfl_alloc_treebuffer(total_num_trees);
    
    // Intialising the traversal on the first tip
    mfl_node_t* current_clip = starttree->treet_treenodes[0];
    // Checking if the first tip is actually a tip
    if(current_clip->nodet_tip != 0) {
        // Some error message?
    }

    // TRAVERSAL START
    
    
    // Create the two different trees (branching and clip)
    // Create the branching tree
    branching_tree = mfl_copy_tree_topology(starttree); //TG: needs some malloc here! Also: use void mfl_disconnect_node(mfl_node_t* n), not sure how to separate the two trees though?
    // Create the clipped tree
    clipped_tree = mfl_copy_tree_topology(starttree); //TG: needs some malloc here!
    
    if (clipped_tree->treet_num_taxa > clade_counter) {
        // rebranch everywhere
        mfl_rebranching_traversal(tree_buffer, branching_tree, clipped_tree, 0); //TG: don't apply neighbour rule
        // incrementing clade counter
        ++clade_counter;
    } else {
        // do the traversal with neigbouring rule
        mfl_rebranching_traversal(tree_buffer, branching_tree, clipped_tree, 1); //TG: don't apply neighbour rule
    }
    
    // Move to the next node
    current_clip = current_clip->nodet_next;
    
    do {
        // Go through the next nodes and print the next tip
//        algorithm_idea(); // TG: I think the difficulty is to make sure that the next node is never +2 taxa bigger than the previous one (to increment the clade_counter only one at the time).
        
        // Move to the next node
        current_clip = current_clip->nodet_next;
        
        // Stop once reached back the start tree
    } while (current_clip != starttree->treet_treenodes[0]);
    
    
    return tree_buffer;
}




void mfl_subtree_pruning_and_regrafting(mfl_tree_t* starttree, mfl_searchrec_t* searchrec)
{
    // Call mfl_pruning_traversal() on the start
}

void tui_spr_test_environment(void)
{
    // NOTE: This is a tui function and should be moved there.
}

bool mfl_heuristic_search(mfl_handle_s *mfl_handle)
{
    /* Eventually, this will parse the information in the handle and set up the 
     * heuristic search. For now, it just returns zero because we're not even 
     * close to that yet. */
    return 0;
}