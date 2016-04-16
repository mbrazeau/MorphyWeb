/* 
 * mfl_tree.c
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
 
 
 
/*
 * In this file:
 * Functions for making, destroying, and performing basic operations on trees. 
 * In this file you will find functions that allocate trees, allocate nodes for 
 * trees. Additionally, tools for basic operations like connecting nodes, 
 * separating nodes, inserting new branches, removing branches etc. are found 
 * here. Functions for freeing all memory used by mfl_tree_t's and mfl_node_t's 
 * are found in here. No functions for assembling trees are contained in this 
 * file. Such operations are handled in separate files depending on whether the 
 * tree is assembled according to a user-specified input (e.g. a Newick or 
 * PhyloXML input, as those become supported) or a criterion-based tree 
 * assembler, such as stepwise addition.
 * 
 */

#include "morphy.h"

int mfl_calculate_number_of_nodes_to_allocate(int num_taxa)
{
    /* This function allows us to size the node array in the tree and correctly 
     * allocate sufficient memory for a tree.
     */
    
    int num_nodes = 0;
    num_nodes = num_taxa + 3 * (num_taxa - 1);
    
    return num_nodes;
}


bool mfl_check_node_is_bottom(mfl_node_t *querynode)
{
    bool is_bottom = 0;
    
    if (querynode->nodet_isbottom) {
        is_bottom = 1;
    }
    
    return is_bottom;
}

void mfl_safe_reset_node_params(mfl_node_t* node)
{
    /* If a node is made available, reset safe values here */
    node->nodet_isbottom         = NULL;
    node->nodet_maxsteps         = NULL;
    node->nodet_minsteps         = NULL;
    node->nodet_downpass_visited = NULL;
    node->nodet_uppass_visited   = NULL;
}

void mfl_join_node_edges(mfl_node_t *node1, mfl_node_t *node2)
{
    node1->nodet_edge = node2;
    node2->nodet_edge = node1;
}


void mfl_disconnect_node_edges(mfl_node_t *node1, mfl_node_t *node2)
{
    node1->nodet_edge = NULL;
    node2->nodet_edge = NULL;
}


void mfl_make_node_available(mfl_node_t *node)
{
    node->nodet_next = NULL;
    node->nodet_edge = NULL;
    mfl_safe_reset_node_params(node);
}


bool mfl_node_is_available(mfl_node_t *node)
{
    bool is_available = false;
    
    if (node->nodet_next == NULL) {
        if (node->nodet_edge == NULL) {
            if (!node->nodet_tip) {
                is_available = true;
            }
        }
    }
    
    return is_available;
}


mfl_node_t * mfl_get_next_available_node(mfl_nodearray_t nodearray)
{
    int i = 0;
    bool success = false;
    mfl_node_t *available_node = NULL;
    
    for (i = 0; nodearray; ++i) {
        if (mfl_node_is_available(*nodearray)) {
            available_node = *nodearray;
            success = true;
            break;
        }
        else {
            ++nodearray;
        }
    }
    
#ifdef MFY_DEBUG
    if (!success) {
        dbg_printf("ERROR in mfl_get_next_available_node(): unable to find unused node\n");
        dbg_printf("All nodes either used or pointing to garbage\n\b");
    }
#endif
    
    return available_node;
}


mfl_node_t * mfl_remove_branch(mfl_node_t *free_node_bottom, mfl_node_t *free_node_top)
{
    mfl_node_t * source_branch_upper = NULL;
    mfl_node_t * source_branch_lower = NULL;
    mfl_node_t * ptr_to_removed_branch = NULL;
    
#ifdef MFY_DEBUG
    if (!free_node_bottom->nodet_isbottom) {
        dbg_printf("WARNING in mfl_remove_branch(): free_node_bottom has NULL nodet_isbottom\n");
    }
#endif
    /*Put the business code in */
    source_branch_lower = free_node_bottom->nodet_edge;
    source_branch_upper = free_node_top->nodet_edge;
    
    free_node_bottom->nodet_edge = NULL;
    free_node_top->nodet_edge = NULL; // Removes pointers to the source tree
    
    source_branch_upper->nodet_edge = source_branch_lower;
    source_branch_lower->nodet_edge = source_branch_upper;  // "Reconnects" the "broken" source tree
    
    ptr_to_removed_branch = free_node_bottom;
    
    return ptr_to_removed_branch;
}


void mfl_insert_branch(mfl_node_t *src_bottom_node, mfl_node_t *src_free_desendant_edge, mfl_node_t *tgt_branch_bottom)
{
    mfl_node_t * tgt_branch_top = NULL;
    tgt_branch_top = tgt_branch_bottom->nodet_edge;

#ifdef MFY_DEBUG
    if (!mfl_check_node_is_bottom(src_bottom_node)) {
        dbg_printf("WARNING in mfl_insert_branch(): src_bottom_node has NULL nodet_isbottom\n");
    }
    if (!mfl_check_node_is_bottom(tgt_branch_bottom)) {
        dbg_printf("WARNING in mfl_insert_branch(): tgt_branch_bottom has NULL nodet_isbottom\n");
    }
#endif
    
    // Point the target branches to the source node
    tgt_branch_bottom->nodet_edge   = src_free_desendant_edge;
    tgt_branch_top->nodet_edge      = src_bottom_node;
    
    // Point the source branches to the targe nodes
    src_bottom_node->nodet_edge         = tgt_branch_top;
    src_free_desendant_edge->nodet_edge = tgt_branch_top;
}


void mfl_make_ring(mfl_node_t *bottom_node, mfl_node_t *left_node, mfl_node_t *right_node)
{
    bottom_node->nodet_next = left_node;
    left_node->nodet_next   = right_node;
    right_node->nodet_next  = bottom_node;
}


mfl_node_t * mfl_alloc_node(void)
{
    mfl_node_t *newnode = NULL;
    
    newnode = (mfl_node_t*)malloc(sizeof(mfl_node_t));
    
    if (newnode == NULL)
    {
        dbg_printf("Error in mfl_alloc_node(void): failed to allocate new mfl_node_t.\n");
    }
    else
    {
        memset(newnode, 0, sizeof(mfl_node_t));
    }
    
    return newnode;
}


void mfl_free_node(mfl_node_t *node)
{
    /*
     * If any other memory allocation made in nodes, then calls to free that 
     * memory should be placed here.
     */
    
    free(node);
}


void mfl_free_treenodes(mfl_nodearray_t treenodes)
{
    int counter = 0;
    
    do {
        ++counter;
        if (*treenodes) {
            mfl_free_node(*treenodes);
        }
        ++treenodes;
    } while (*treenodes);
    
}


void mfl_allocate_nodes_in_array(mfl_nodearray_t nodearray, int num_nodes, int num_taxa)
{
    int i = 0;
    
    if (!num_nodes) {
        if (!num_taxa) {
            dbg_printf("Error in mfl_allocate_nodes_in_array: insufficient data for sizing node array\n");
        }
        else {
            num_nodes = mfl_calculate_number_of_nodes_to_allocate(num_taxa);
        }
    }
    
    for (i = 0; i < num_nodes; ++i) {
        *(nodearray + i) = mfl_alloc_node();
    }
    
    /* The last spot in the array is a NULL pointer. This just provides a bit of
     * extra safety in case of an attempt to read past the allocated nodes.
     */
    *(nodearray + num_nodes) = NULL;
    
}


void mfl_setup_nodearray(mfl_nodearray_t nodearray, int num_nodes, int num_taxa)
{
    int i = 0;
    
    for (i = 0; nodearray[i]; ++i) {
        nodearray[i]->nodet_index = i;
        if (i < num_taxa) {
            nodearray[i]->nodet_tip = i + 1;
        }
        else {
            nodearray[i]->nodet_tip = 0;
        }
    }
}


mfl_nodearray_t mfl_allocate_nodearray(int num_taxa, int num_nodes)
{

    mfl_nodearray_t new_nodearray = NULL;
    
    if (!num_nodes) {
        if (!num_taxa) {
            dbg_printf("Error in mfl_allocate_node_array: insufficient data for sizing node array\n");
        }
        else {
            num_nodes = mfl_calculate_number_of_nodes_to_allocate(num_taxa);
        }
    }
    
    
    new_nodearray = (mfl_node_t**)malloc( (num_nodes + 1) * sizeof(mfl_node_t*)); // The +1 allows for a NULL pointer at then end of the array to (hopefully) keep it read-safe.
    
    if (!new_nodearray) {
        dbg_printf("Error in mfl_allocate_node_array(int num_taxa, int num_nodes): failed to allocate new node array.\n");
    }
    else
    {
        memset(new_nodearray, 0, num_nodes * sizeof(mfl_node_t*));
    }
    
    mfl_allocate_nodes_in_array(new_nodearray, num_nodes, num_taxa);
    
    return new_nodearray;
}


void mfl_free_nodearray(mfl_nodearray_t nodearray)
{
    free(nodearray);
}


mfl_node_t *mfl_insert_node_in_ring(mfl_node_t *ring_start, mfl_node_t *new_node)
{
    mfl_node_t *last_in_ring = NULL;
    
    if (ring_start == new_node) {
        /* This prevents setting a node's nodet_next pointer to itself, thus 
         * ensuring that the nodet_next pointer will always point to either 
         * another node or NULL. It is not that this is necessarily bad, but 
         * this operation would seem to be unlikely. This rule will help avoid 
         * unexpected behaviour by other functions */
        dbg_printf("Warning in function calling mfl_insert_node_in_ring(): cannot point nodet_next of one node to itself\n");
        return NULL;
    }
    
    last_in_ring = ring_start;
    
    /* This conditional statement allows the use of this function to create a 
     * ring by inserting the nodes one at a time. It checks that both the next 
     * pointer in the starting node is available and that new_node is not, in 
     * fact, NULL. Otheriwise it would create a node to point its next pointer 
     * to itself and we would like to avoid that. */
    
    if (!ring_start->nodet_next && new_node) {
        ring_start->nodet_next = ring_start;
    }
    
    do {
        last_in_ring=last_in_ring->nodet_next;
    } while (last_in_ring->nodet_next != ring_start);
    
    last_in_ring->nodet_next = new_node;
    new_node->nodet_next = ring_start;
    
    return ring_start;
}


bool mfl_check_is_in_ring(mfl_node_t *start)
{
    bool is_ring = false;
    bool did_loop = false;
    mfl_node_t *node_p = NULL;
    
    node_p = start;
    
    do {
        if (node_p->nodet_next) {
            node_p = node_p->nodet_next;
            did_loop = true;
        }
        else {
            return false;
        }
    } while (node_p != start);
    
    if (node_p == start) {
        if (did_loop) {
            is_ring = true;
        }
        else {
            is_ring = false;
        }
    }
    
    return is_ring;
}


void mfl_initialise_ring_node(mfl_node_t *bottom_node)
{
    mfl_node_t *p = NULL;
    
    if (!bottom_node->nodet_isbottom) {
        dbg_printf("Error in function calling mfl_initialise_ring_node(): node passed is not a bottom node.\n");
        return;
    }
    
    do {
        p = bottom_node->nodet_next;
        p->nodet_isbottom = 0;
        p->nodet_index = bottom_node->nodet_index;
        
        /* Dependeing on how we manage character data at internal node rings, we
         * might choose to point all the charstate arrays to a single piece of 
         * data, rather than each node having its own character memory */
        
    } while (p != bottom_node);
    
}


mfl_node_t *mfl_make_new_n_ary_ring_node(mfl_node_t *bottom_node, int num_branches, mfl_nodearray_t nodes)
{
    int i = 0;
    mfl_node_t *new_node = NULL;
    mfl_node_t *node_p = NULL;
    
    if (!mfl_node_is_available(bottom_node)) {
        dbg_printf("Warning in function calling mfl_make_new_n_ary_ring_node(): selected node is unavailable. Return is NULL pointer.\n");
        return NULL;
    }
    
    node_p = bottom_node;
    
    for (i = 0; i < num_branches; ++i) {
        new_node = mfl_get_next_available_node(nodes);
        node_p->nodet_next = new_node;
        node_p = new_node;
    }
    
    node_p->nodet_next = bottom_node;
    

    if (!mfl_check_is_in_ring(bottom_node)) {
        bottom_node = NULL;
        dbg_printf("Warning in mfl_make_new_n_ary_ring_node(): did not succeed in making %i-branch ring. Returning NULL pointer\n", num_branches);
        return bottom_node;
    }
    
    return bottom_node;
}


void mfl_destroy_n_nary_ring(mfl_node_t *bottom_node)
{
    mfl_node_t *p = NULL;
    mfl_node_t *garbage_node = NULL;
    
    p = bottom_node->nodet_next;
    
    do {
        garbage_node = p;
        p = p->nodet_next;
        free(garbage_node);
        
        /* Other free() calls may go here when destroying whole trees after 
         * analysis*/
        
    } while (p != bottom_node);
    
    bottom_node->nodet_next = NULL;
    
}


void mfl_create_binary_fork(mfl_node_t *parent, mfl_node_t *child1, mfl_node_t *child2, mfl_nodearray_t nodes)
{
    if (!parent->nodet_next) {
        mfl_make_new_n_ary_ring_node(parent, 2, nodes);
    }
    else {
        dbg_printf("Warning in mfl_create_binary_fork(): parent node might be unavailable\n");
    }
    
    mfl_join_node_edges(parent->nodet_next, child1);
    mfl_join_node_edges(parent->nodet_next->nodet_next, child2);
}


void mfl_initialise_nodearray(mfl_nodearray_t nodearray, int num_taxa, int num_nodes)
{
    int i = 0;
    
    for (i = 0; i < num_nodes; ++i) {
        nodearray[i]->nodet_index = i;
        if (i < num_taxa) {
            nodearray[i]->nodet_tip = i + 1;
        }
    }
}


/*!
 Returns 0 if the node is binary, -1 if not defined (no branching) and otherwise 
 returns the number of branchings detected. 
 @param querynode (mfl_node_t*) the base of the internal node being queried for n-ariness
 @param test_n_branches (int) the expected number of descendant branches
 @return the number of descendant branches found if true, 0 if number of descendants does 
 not match the expected value.
 */
int mfl_node_is_n_ary(mfl_node_t *querynode, int test_n_branches)
{
    mfl_node_t *node_ptr = NULL;
    int num_branching = 0;
    
    
    if (!querynode->nodet_next) {
        dbg_printf("WARNING in mfl_node_is_n_ary(): query node has no valid nodet_next pointer\n");
        return -1;
    }
    
    node_ptr = querynode;
    do {
        node_ptr = node_ptr->nodet_next;
        ++num_branching;
    } while (node_ptr != querynode);
    
    if (num_branching == test_n_branches) {
        return num_branching;
    }
    else {
        return 0;
    }
    
}


mfl_node_t* mfl_find_rightmost_tip_in_tree(mfl_node_t* n)
{
    
    mfl_node_t *p;
    
    if (n->nodet_tip) {
        return n;
    }
    
    return mfl_find_rightmost_tip_in_tree(n->nodet_next->nodet_edge);
}


void mfl_unroot_tree(mfl_tree_t *tree)
{
    mfl_node_t *p = NULL;
    mfl_node_t *q = NULL;
    
    if (!tree->treet_root) {
        dbg_printf("WARNING in mfl_unroot_tree(): attempt to deroot tree with no valid pointer to root\n");
        return;
    }
    
    q = tree->treet_root->nodet_next;
    p = tree->treet_root;
    
    do {
        p = p->nodet_next;
    } while (p->nodet_next != tree->treet_root);
    
    if (p == q->nodet_next) {
        // The root is binary; perform binary unrooting
        mfl_join_node_edges(p->nodet_edge, q->nodet_edge);
        mfl_make_node_available(p);
        mfl_make_node_available(q);
    }
    else {
        // The root is polychotomous, perform polychotomous unrooting
        p->nodet_next = q;
    }
    
    mfl_make_node_available(tree->treet_root);
    
    // This bit sets the root as either the base of tip #1 or the base of some arbitary
    // tip. This can be changed, but the start pointer should always be at the base
    // of some tip. This keeps uppass calculations simple.
    if (tree->treet_treenodes[0]->nodet_edge) {
        tree->treet_start = tree->treet_treenodes[0]->nodet_edge;
    }
    else {
        q = mfl_find_rightmost_tip_in_tree(q);
        tree->treet_start = q->nodet_edge;
    }
    
    tree->treet_root = NULL;
    
}


void mfl_initialise_tree(mfl_tree_t *newtree, int num_taxa, int num_nodes)
{
    mfl_initialise_nodearray(newtree->treet_treenodes, num_taxa, num_nodes);
    
    newtree->treet_num_taxa = num_taxa;
}

mfl_tree_t * mfl_alloctree_with_nodes(int num_taxa)
{
    int num_nodes = 0;
    
    num_nodes = mfl_calculate_number_of_nodes_to_allocate(num_taxa);
    
    mfl_tree_t *newtree = NULL;
    
    newtree = (mfl_tree_t*)malloc(sizeof(mfl_tree_t));
    memset(newtree, 0, sizeof(mfl_tree_t));
    
    newtree->treet_treenodes = mfl_allocate_nodearray(num_taxa, num_nodes);
    
    mfl_initialise_tree(newtree, num_taxa, num_nodes);
    
    return newtree;
}


void mfl_free_tree(mfl_tree_t *tree_to_free)
{
    
    mfl_free_treenodes(tree_to_free->treet_treenodes);
    mfl_free_nodearray(tree_to_free->treet_treenodes);
    
    /*
     * Any other allocated memory in a tree should be freed here
     */
    
    // Free the tree
    free(tree_to_free);
}