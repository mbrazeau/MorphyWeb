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

/*!
 @discussion Gives the number of nodes required to build an mfl_tree_t with 
 3-node internal node ring cycles for a tree with num_taxa leaves.
 @param num_taxa (int) the number of terminal leaves
 @return the number of nodes required to build a completely binary, rooted 
 tree.
 */
int mfl_calculate_number_of_nodes_to_allocate(int num_taxa)
{
    int num_nodes = 0;
    num_nodes = num_taxa + 3 * (num_taxa - 1);
    
    return num_nodes;
}


/*!
 @discussion Verifies that a node points towards the root of the tree through 
 its nodet_edge pointer. Useful only in rooted trees or trees, or where knowing 
 where the root used to be is important.
 @param querynode (mfl_node_t*) the node to be checked.
 @return true if node points towards the root, otherwise false
 */
bool mfl_check_node_is_bottom(mfl_node_t *querynode)
{
    bool is_bottom = 0;
    
    if (querynode->nodet_isbottom) {
        is_bottom = 1;
    }
    
    return is_bottom;
}


/*!
 @discussion When a node is made available, this function is called to reset
 internal node parameters that are dependent on the topology. Any new topology-
 dependent parameters need to be added here.
 @param node (mfl_node_t*) the node to be reset.
 */
void mfl_safe_reset_node_params(mfl_node_t* node)
{
    /* If a node is made available, reset safe values here */
    node->nodet_isbottom         = NULL;
    node->nodet_maxsteps         = NULL;
    node->nodet_minsteps         = NULL;
    node->nodet_downpass_visited = NULL;
    node->nodet_uppass_visited   = NULL;
}


/*!
 @discussion Joins two nodes at their edge pointers
 @param node1 (mfl_node_t*) one side of the union
 @param node2 (mfl_node_t*) the other side of the union
 */
void mfl_join_node_edges(mfl_node_t *node1, mfl_node_t *node2)
{
    node1->nodet_edge = node2;
    node2->nodet_edge = node1;
}


/*!
 @discussion Disconnects the reciprocal edge pointers of two nodes.
 @note: this function provides no guarantees that the input nodes are 
 reciprocal. Thus, it is safer to call mfl_disconnect_node() to ensure that that
 @param node1 (mfl_node_t*) the first node of the union to break.
 @param node2 (mfl_node_t*) What *should* be the second node of the union to
 break.
 */
void mfl_disconnect_node_edges(mfl_node_t *node1, mfl_node_t *node2)
{
    node1->nodet_edge = NULL;
    node2->nodet_edge = NULL;
}


/*!
 @discussion Breaks a node from its reciprocal edge. This operation will result
 in breaking a tree into subtrees.
 @param n (mfl_node_t*) the node at the desired break point.
 */
void mfl_disconnect_node(mfl_node_t* n)
{
    mfl_node_t* p = n->nodet_edge;
    mfl_disconnect_node_edges(p, n);
}


/*!
 @discussion Makes a node available for re-use in another part of the tree by
 nullifying its connections to other nodes and
 @param node (mfl_node_t*) the node to be made available
 */
void mfl_make_node_available(mfl_node_t *node)
{
    node->nodet_next = NULL;
    node->nodet_edge = NULL;
    mfl_safe_reset_node_params(node);
}


/*!
 @discussion Checks whether an internal node is available by checking whether or
 not it has connections to other nodes.
 @param node (mfl_node_t*) the node to be queried
 @return true if available, false if unavailable (i.e. has connections).
 */
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


mfl_nodestack_t* mfl_create_empty_nodestack(int num_internal_nodes)
{
    mfl_nodestack_t* newndstk = NULL;
    
    newndstk = (mfl_nodestack_t*)malloc(sizeof(mfl_nodestack_t));
    if (!newndstk) {
        dbg_eprintf("unable to allocate memory for new nodestack");
    }
    else {
        memset(newndstk, 0, sizeof(mfl_nodestack_t));
    }
    
    newndstk->nstk_availbale_nds = (mfl_nodearray_t)malloc((num_internal_nodes + 1) * sizeof(mfl_node_t*));
    
    if (!newndstk->nstk_availbale_nds) {
        dbg_eprintf("unable to allocate memory for node array in new nodestack");
    }
    else {
        memset(newndstk->nstk_availbale_nds, 0, (num_internal_nodes + 1) * sizeof(mfl_node_t*));
    }
    
    newndstk->nstk_availbale_nds[num_internal_nodes] = NULL;
    
    newndstk->nstk_maxsize = num_internal_nodes;
    
    return newndstk;
}


void mfl_destroy_nodestack(mfl_nodestack_t* ndstk)
{
    free(ndstk->nstk_availbale_nds);
    free(ndstk);
}


void mfl_push_node_to_nodestack(mfl_node_t* n, mfl_nodestack_t* ndstk)
{
    mfl_nodestack_t* nds;
    
    if (!ndstk) {
         nds = n->nodet_ndstack;
    }
    else {
        nds = ndstk;
    }
    
    mfl_make_node_available(n);
    
    // !!!: Other cleanup operations might go here.
    
    if (nds->nstk_numnodes < nds->nstk_maxsize) {
        nds->nstk_availbale_nds[nds->nstk_numnodes-1] = n;
        ++nds->nstk_numnodes;
    }
    else {
        dbg_eprintf("insufficient space for new node on nodestack");
        return;
    }
}


mfl_node_t* mfl_get_node_from_nodestack(mfl_nodestack_t *nds)
{
    mfl_node_t* retnode = NULL;
    
    --nds->nstk_numnodes;
    retnode = nds->nstk_availbale_nds[nds->nstk_numnodes];
    nds->nstk_availbale_nds[nds->nstk_numnodes] = NULL;
    
    return retnode;
}


/*!
 @discussion Performs a linear search in a node array to find a node that is 
 not joined to any other nodes by either a next pointer or an edge pointer.
 @param nodearray (mfl_nodearray_t) the array of nodes to be searched.
 @return pointer to an available node.
 */
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


///*!
// @discussion Removes a branch from a target tree and returns a pointer to the 
// base of the excised node.
// @note I think there can be some big changes to this function and to generally
// how we handle removal and reinsertion of branches.
// @param free_node_bottom (mfl_node_t*) pointer to the base of the internal node 
// ring to be removed and which connects the node to the target tree.
// @param free_node_top (mfl_node_t*) pointer to the 'upper' node ring that 
// connects the source branch to the target tree.
// @return pointer to the ring forming the base of the excised branch.
// */
//mfl_node_t * mfl_remove_branch(mfl_node_t *free_node_bottom, mfl_node_t *free_node_top)
//{
//    mfl_node_t * tgt_branch_upper = NULL;
//    mfl_node_t * tgt_branch_bottom = NULL;
//    mfl_node_t * ptr_to_removed_branch = NULL;
//    
//#ifdef MFY_DEBUG
//    if (!free_node_bottom->nodet_isbottom) {
//        dbg_printf("WARNING in mfl_remove_branch(): free_node_bottom has NULL nodet_isbottom\n");
//    }
//#endif
//    /*Put the business code in */
//    tgt_branch_bottom = free_node_bottom->nodet_edge;
//    tgt_branch_upper = free_node_top->nodet_edge;
//    
//    free_node_bottom->nodet_edge = NULL;
//    free_node_top->nodet_edge = NULL; // Removes pointers to the source tree
//    
//    tgt_branch_upper->nodet_edge = tgt_branch_bottom;
//    tgt_branch_bottom->nodet_edge = tgt_branch_upper;  // "Reconnects" the
//                                                       // "broken" source tree
//    
//    ptr_to_removed_branch = free_node_bottom;
//    
//    return ptr_to_removed_branch;
//}

/*!
 @discussion Inserts a ring with only one edge connection into a target tree.
 @param src (mfl_node_t*) the base of the source node
 @param tgt (mfl_node_t*) a node straddling the target destination.
 */
void mfl_insert_branch(mfl_node_t *src, mfl_node_t *tgt)
{
    mfl_node_t* srcf = NULL;
    mfl_node_t* tgt_opp = tgt->nodet_edge;
    
    assert(!src->nodet_next->nodet_edge || !src->nodet_next->nodet_next->nodet_edge);
    
    if (src->nodet_next->nodet_edge) {
        srcf = src->nodet_next;
    } else {
        srcf = src->nodet_next;
    }
    
    mfl_join_node_edges(src, tgt);
    mfl_join_node_edges(srcf, tgt_opp);
}


/*!
 @discussion Sets three input nodes into a ring cycle for use as an internal
 node in a tree.
 @param bottom_node (mfl_node_t*) the node that seeds the ring and forms the 
 bottom "rootward" connection if directed
 @param left_node (mfl_node_t*) the next node after the root
 @param right_node (mfl_node_t*) the next node after the left node and joins the
 bottom node via its next pointer
 */
void mfl_make_ring(mfl_node_t *bottom_node, mfl_node_t *left_node, mfl_node_t *right_node)
{
    bottom_node->nodet_next = left_node;
    left_node->nodet_next   = right_node;
    right_node->nodet_next  = bottom_node;
}


/*!
 @discussion Wrapper function for node allocation and checking for any failures
 to allocate memory
 @return pointer to an mfl_node_t
 */
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


/*!
 @discussion Wrapper function on free() which is used to free node memory. 
 @note Any memory allocated within a node should be freed here as well.
 @param node (mfl_node_t*) pointer to the node that will be free.
 */
void mfl_free_node(mfl_node_t *node)
{
    /*
     * If any other memory allocation made in nodes, then calls to free that 
     * memory should be placed here.
     */
    mfl_bts_destroy_bitset(node->nodet_bipart);
    free(node);
}


/*!
 @discussion Loops through the node pointer array freeing each node.
 @param treenodes (mfl_nodearray_t) an array of pointers to nodes froma source
 tree.
 */
void mfl_free_treenodes(mfl_nodearray_t treenodes)
{

    do {
        if (*treenodes) {
            mfl_free_node(*treenodes);
        }
        ++treenodes;
    } while (*treenodes);
    
}


/*!
 @discussion After an array of pointers to nodes has been sized, this function
 is used to allocate memory for each node pointer, thereby allowing a tree to be
 created in memory. 
 @note the last node pointer in the array is a NULL pointer and signals the end
 of the array.
 @note a calling function must specify either num_nodes or num_taxa but does not
 need to specify both
 @param nodearray (mfl_nodearray_t) the array of empty node pointers.
 @param num_nodes (int) the number of nodes to allocate.
 @param num_taxa (int) the number of terminal taxa.
 */
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
        *(nodearray + i) = (mfl_node_t*)mfl_malloc(sizeof(mfl_node_t), 0, __FXN_NAME__);//mfl_alloc_node();
    }
    
    /* The last spot in the array is a NULL pointer. This just provides a bit of
     * extra safety in case of an attempt to read past the allocated nodes.
     */
    *(nodearray + num_nodes) = NULL;
    
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
    
    // The +1 allows for a NULL pointer at then end of the array to (hopefully)
    // keep it read-safe.
    new_nodearray = (mfl_node_t**)mfl_malloc((num_nodes + 1) * sizeof(mfl_node_t), 0, __FXN_NAME__);
    
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


void mfl_initialise_nodearray(mfl_tree_t* t, int num_taxa, int num_nodes)
{
    int i = 0;
    int j = 0;
    mfl_nodearray_t nodearray = t->treet_treenodes;
    
    for (i = 0; i < num_nodes; ++i) {
        nodearray[i]->nodet_index = i;
        nodearray[i]->nodet_ndstack = t->treet_nodestack;
        nodearray[i]->nodet_bipart = mfl_bts_create_bitset(num_taxa);
        if (i < num_taxa) {
            nodearray[i]->nodet_tip = i + 1;
            nodearray[i]->nodet_weight = 1;
            mfl_bts_setbit(nodearray[i]->nodet_bipart, 1, nodearray[i]->nodet_tip);
        }
        else {
            t->treet_nodestack->nstk_availbale_nds[j] = nodearray[i];
            ++j;
            t->treet_nodestack->nstk_numnodes = j;
        }
    }
}


/*!
 Returns 0 if the node is binary, -1 if not defined (no branching) and otherwise 
 returns the number of branchings detected. 
 @param querynode (mfl_node_t*) the base of the internal node being queried for 
 n-ariness
 @param test_n_branches (int) the expected number of descendant branches
 @return the number of descendant branches found if true, 0 if number of 
 descendants does not match the expected value.
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
    
    if (n->nodet_tip) {
        return n;
    }
    
    return mfl_find_rightmost_tip_in_tree(n->nodet_next->nodet_edge);
}


void mfl_unroot_tree(mfl_tree_t *tree)
{
    mfl_node_t *p = NULL;
    mfl_node_t *q = NULL;
    
    assert(tree->treet_root);
    
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
    mfl_initialise_nodearray(newtree, num_taxa, num_nodes);
    
    newtree->treet_num_taxa = num_taxa;
    newtree->treet_num_nodes = num_nodes;
    newtree->treet_dummynode.nodet_tip = -1;
}


mfl_tree_t * mfl_alloctree_with_nodes(int num_taxa)
{
    int num_nodes = 0;
    
    num_nodes = mfl_calculate_number_of_nodes_to_allocate(num_taxa);
    
    mfl_tree_t *newtree = NULL;
    
    newtree = (mfl_tree_t*)mfl_malloc(sizeof(mfl_tree_t), 0, __FXN_NAME__);
    
    newtree->treet_treenodes = mfl_allocate_nodearray(num_taxa, num_nodes);
    
    newtree->treet_nodestack = mfl_create_empty_nodestack(num_nodes - num_taxa);
    
    mfl_initialise_tree(newtree, num_taxa, num_nodes);
    
    return newtree;
}


void mfl_free_tree(mfl_tree_t *tree_to_free)
{
    
    mfl_free_treenodes(tree_to_free->treet_treenodes);
    mfl_free_nodearray(tree_to_free->treet_treenodes);
    mfl_destroy_nodestack(tree_to_free->treet_nodestack);
    /*
     * Any other allocated memory in a tree should be freed here
     */
    
    // Free the tree
    free(tree_to_free);
}


/*!
 Add a root to a node ring (creating a polytomy)
 @param input_tree a pointer to a mfl_tree_t to root
 @param target_node_ring_start a pointer to the node to root
 */
void mfl_root_target_node(mfl_tree_t *input_tree, mfl_node_t *target_node_ring_start)
{
    //The tree must not be rooted!
    if (input_tree->treet_root){
        if (!input_tree->treet_root->nodet_edge) {
            dbg_printf("ERROR in mfl_root_target_node(): the input tree is already rooted!"); // Addition to the message for the profane version: "the fuck you think you're doing with that root?"
            return;
        }
        else {
            dbg_eprintf("tree has root with unexpected nodet_edge pointer non-NULL value\n");
            return;
        }    } else {
        
        //Generate the new root node
        mfl_node_t *root_node;
        root_node = mfl_get_next_available_node(input_tree->treet_treenodes);
        
        //Connect the node to the ring
        mfl_insert_node_in_ring(target_node_ring_start, root_node);
        
        //Setting it's nodet_edge to NULL;
        root_node->nodet_edge = NULL;
        
        //Check if node is ring
        bool check_ring;
        check_ring = mfl_check_is_in_ring(target_node_ring_start);
        if(check_ring == false) {
            dbg_printf("ERROR in mfl_root_target_node(): Rooting target node broke the node ring!");
        }
        
        //Point the tree root to the root_node
        input_tree->treet_root = root_node;
    }
}

/*!
 Add a root to an edge
 @param input_tree a pointer to a mfl_tree_t to root
 @param target_node a pointer to the node linking to the next edge to root
 */
void mfl_root_target_edge(mfl_tree_t *input_tree, mfl_node_t *target_node)
{
    int i = 0;
    int num_branches = 2;
    mfl_node_t *root_node = NULL;
    mfl_node_t *n = NULL;
    mfl_nodearray_t nds = input_tree->treet_treenodes; // This is an optimisation trick, more than anything.
                                                       // Because this dereferencing operation is done a bunch
                                                       // of times in this function, we can just store the address
                                                       // of treet_treenodes in a local variable do avoid multiple
                                                       // dereferences
    
    //The tree must not be rooted!
    if (input_tree->treet_root){
        if (!input_tree->treet_root->nodet_edge) {
            dbg_printf("ERROR in mfl_root_target_node(): the input tree is already rooted!"); // Addition to the message for the profane version: "the fuck you think you're doing with that root?"
            return;
        }
        else {
            dbg_eprintf("tree has root with unexpected nodet_edge pointer non-NULL value\n");
            return;
        }
    } else {
        
        //Generate the new root node
        
        root_node = mfl_get_next_available_node(nds);
        // Temporary linking the root to itseld
        root_node->nodet_next = root_node;
        
        // MDB: I'd do it like this:
        do {
            n = mfl_get_next_available_node(nds);
            mfl_insert_node_in_ring(root_node, n);
            ++i;
        } while (i < num_branches);
        // This automates what you've done below. 
        
        i = 0; // Resent i to 0 in case we try to use it again.
        
        //Generate the two other nodes in the new root ring and pointing them to themselves
        /*mfl_node_t *root_ring_node_left;
        root_ring_node_left = mfl_get_next_available_node(input_tree->treet_treenodes);
        root_ring_node_left->nodet_next = root_ring_node_left;
        mfl_node_t *root_ring_node_right;
        root_ring_node_right = mfl_get_next_available_node(input_tree->treet_treenodes);
        root_ring_node_right->nodet_next = root_ring_node_right;*/

        
        //Generate the root node ring
        //mfl_make_ring(root_node, root_ring_node_left, root_ring_node_right); // bottom is pointing to NULL (root), left is pointing to target_node and right is pointing to target_node->edge
        
        
        //Check if node is ring //TG: this part might be over kill since mfl_make_ring is supose to work smoothly. Just to be sure though.
        /*
         MDB: Yep. This is a debugging check and will slow execution down quite a bit if rootings/
         rerootings are done often enough. If you're going to add these kinds of checks, I'd
         put them inside conditional debug statements, like so: */
        
#ifdef MFY_DEBUG
        bool check_ring;
        check_ring = mfl_check_is_in_ring(root_node);
        if(check_ring == false) {
            dbg_printf("ERROR in mfl_root_target_edge(): Rooting node ring is broken!");
        }
        /* Now this check won't be compiled into the release build and slow performance.
         * However, we can get rid of it altogether and let the work be done by our test
         * interface. */
#endif
        
        //Join the new edges
        mfl_join_node_edges(target_node->nodet_edge, root_node->nodet_next);
        mfl_join_node_edges(target_node, root_node->nodet_next->nodet_next);
        
        //Setting it's nodet_edge to NULL;
        root_node->nodet_edge = NULL;
        
        //Point the tree root to the root_node
        input_tree->treet_root = root_node;
        
        /*
         *  TG: It's probably safe to add a broken tree checker here no?
         *  MDB: TUI functions aren't exported to the library. They're for testing only.
         *  Given that rooting and unrooting might be performed many times during the search,
         *  I woudld not want to add the complete tree-check routine to the mix. If we did
         *  add it here, we should put it within conditional compilation statements so that 
         *  it's not included in the release builds.
         */
    }
}


/*!
 Increase the treebuffer size
 @param trbuf (*mfl_treebuffer_t) a pointer to a tree buffer
 @param addedlength (int) the number of spaces to add to the treebuffer
 */
void mfl_resize_treebuffer(mfl_treebuffer_t* trbuf, int addedlength)
{
    
    int newsize = trbuf->tb_max_buffersize + addedlength;
    
    mfl_tree_t** newtarray = (mfl_tree_t**)malloc(newsize * sizeof(mfl_tree_t*));
    if (!newtarray) {
        dbg_eprintf("unable to increase size of tree buffer");
        return;
    }
    else {
        memset(newtarray, 0, newsize * sizeof(mfl_tree_t*));
    }
    
    memcpy(newtarray, trbuf->tb_savedtrees, (trbuf->tb_num_trees * sizeof(mfl_tree_t*)) );
    
    free(trbuf->tb_savedtrees);
    trbuf->tb_savedtrees = newtarray;
    trbuf->tb_maxtrees = newsize;
}

/*!
 Adds a tree to a treebuffer
 @param newtree (*mfl_tree_t) a pointer to a mfl_tree_t
 @param trbuf (*mfl_treebuffer_t) a pointer to a treebuffer
 @param mfl_handle (*mfl_handle_s) a pointer to the mfl_handle options
 */
void mfl_append_tree_to_treebuffer(mfl_tree_t* newtree, mfl_treebuffer_t* trbuf, mfl_handle_s* mfl_handle)
{
    int addedlength = 0;
    
    if ((trbuf->tb_num_trees+1) > trbuf->tb_max_buffersize) {
        
        if (mfl_handle->autoincrease) {
            if (mfl_handle->autoinc_incr) {
                addedlength = mfl_handle->autoinc_incr;
            }
            else {
                addedlength = MORPHY_DEFAULT_TREEBUFFER_AUTOINCREASE_AMOUNT;
            }
            mfl_resize_treebuffer(trbuf, addedlength);
        }
        else {
            return;
        }
    }
    
    trbuf->tb_savedtrees[trbuf->tb_num_trees] = newtree;
    ++trbuf->tb_num_trees;
}

/*!
 Allocates memory for a treebuffer
 @param num_trees (int) the number of trees to allocate in the buffer
 @return a tree buffer (mfl_treebuffer_t)
 */
mfl_treebuffer_t* mfl_alloc_treebuffer(int num_trees)
{
    int trbufsize = 1;
    mfl_treebuffer_t* newtrbf = NULL;
    
    newtrbf = (mfl_treebuffer_t*)malloc(sizeof(mfl_treebuffer_t));
    if (!newtrbf) {
        dbg_eprintf("unable to allocate memory for new treebuffer");
        return NULL;
    }
    else {
        memset(newtrbf, 0, sizeof(mfl_treebuffer_t));
    }
    
    if (num_trees) {
        trbufsize = num_trees;
    }
    
    newtrbf->tb_savedtrees = (mfl_tree_t**)malloc(trbufsize * sizeof(mfl_tree_t*));
    if (!newtrbf->tb_savedtrees) {
        dbg_eprintf("unable to allocte memory for tree array in new treebuffer");
        free(newtrbf);
        return NULL;
    }
    else {
        memset(newtrbf->tb_savedtrees, 0, trbufsize * sizeof(mfl_tree_t*));
    }
    
    newtrbf->tb_max_buffersize = trbufsize;
    
    return newtrbf;
}

void mfl_reset_nodestack(mfl_nodestack_t* nstk)
{
    int i = 0;
    
    for (i = nstk->nstk_numnodes; i; --i) {
        if (!mfl_node_is_available(nstk->nstk_availbale_nds[i - 1])) {
            nstk->nstk_availbale_nds[i - 1] = NULL;
            --nstk->nstk_numnodes;
        }
    }
}

mfl_tree_t* mfl_copy_tree_topology(const mfl_tree_t* t)
{
    assert(t);
    int i = 0;
    int numnodes = 0;
    mfl_tree_t* trcopy = NULL;
    trcopy = mfl_alloctree_with_nodes(t->treet_num_taxa);
    
    mfl_nodearray_t srcnds = t->treet_treenodes;
    mfl_nodearray_t cpynds = trcopy->treet_treenodes;
    
    assert((numnodes = trcopy->treet_num_nodes) != 0);
    
    for (i = 0; i < numnodes ; ++i) {
        assert(srcnds[i]->nodet_index == cpynds[i]->nodet_index);
        if (srcnds[i]->nodet_edge) {
            cpynds[i]->nodet_edge = cpynds[srcnds[i]->nodet_edge->nodet_index];
        }
        if (srcnds[i]->nodet_next) {
            cpynds[i]->nodet_next = cpynds[srcnds[i]->nodet_next->nodet_index];
        }
    }
    
    // Copy safe variables
    trcopy->treet_num_taxa = t->treet_num_taxa;
    
    assert(!(t->treet_start && t->treet_root));
    
    if (t->treet_root) {
        trcopy->treet_root = cpynds[t->treet_root->nodet_index];
        trcopy->treet_start = NULL;
    }
    else {
        trcopy->treet_start = cpynds[t->treet_start->nodet_index];
        trcopy->treet_root = NULL;
    }
    
    // Reset the nodestack.
    mfl_reset_nodestack(trcopy->treet_nodestack);
    
    return trcopy;
}

/*!
 Frees memory from a treebuffer
 @param oldtreebuf (*mfl_treebuffer_t) the tree buffer to free
 @param cleartrees (bool) whether to also clear the trees linked from the buffer (true) or not (false)
 */
void mfl_destroy_treebuffer(mfl_treebuffer_t* oldtreebuf, bool cleartrees)
{
    if (cleartrees) {
        // Free all trees at pointers in savedtrees,
    }
        
    if (oldtreebuf->tb_savedtrees) {
        free(oldtreebuf->tb_savedtrees);
    }
    
    free(oldtreebuf);
}

//Allocate the memory to the edgetable
mfl_edgetable_t* mfl_initiate_edgetable_t(int num_tips)
{
    //Intialise the table
    mfl_edgetable_t* new_edgetable = NULL;
    
    //malloc the edgetable
    new_edgetable = (mfl_edgetable_t*)malloc(sizeof(mfl_edgetable_t));
    if (!new_edgetable) {
        dbg_eprintf("unable to allocate memory for new edgetable");
    }
    else {
        memset(new_edgetable, 0, sizeof(mfl_edgetable_t));
    }
    
    //malloc the size of the edge table
    new_edgetable->numentries = (2*num_tips - 2) * 2;
                                        
    //malloc the size of the edge table
    new_edgetable->edgetable = (int*)malloc((new_edgetable->numentries) * sizeof(int));
    
    return new_edgetable;
}


void mfl_edgetable_traversal(mfl_node_t* start, mfl_edgetable_t* edgetable, int* node_counter, int* tip_counter, int* counter)
{
    
    mfl_node_t* node = start;
    
    if(node->nodet_tip != 0) {
        // When node is a tip, always increment tip counter and print the tip.
        ++*tip_counter;
        edgetable->edgetable[*counter] = *tip_counter;
    } else {
        
        // Increment node counter only when next is not a tip
        if(node->nodet_next->nodet_tip != 1) {
            ++*node_counter;
        }
        
        edgetable->edgetable[*counter] = *node_counter;
    }
    
    ++*counter;
    node = node->nodet_next;
    
    do {
        mfl_edgetable_traversal(node->nodet_edge, edgetable, node_counter, tip_counter, counter);
        
        node = node->nodet_next;
        
    } while (node != start);
    
    
}

void mfl_get_edgetable(mfl_edgetable_t* edgetable, mfl_tree_t* tree)
{
    // variables before the traversal
    int counter = 0;
    int tip_counter = 0;
    int node_counter = tree->treet_num_taxa;     // node count starts at number of tips + 1
    
    mfl_edgetable_traversal(tree->treet_start, edgetable, &counter, &tip_counter, &node_counter);

}


// Function idea
bool mfl_compare_edge_tables(mfl_edgetable_t* t1, mfl_edgetable_t* t2)
{
    if (t1->numentries != t2->numentries) {
        return false;
    }
    
    if (!memcmp(t1, t2, t1->numentries * sizeof(int))) {
        return true; // memcmp returns 0 if there's a match
    }
    else {
        return false;
    }
}

void mfl_free_edgetable(mfl_edgetable_t* edgetable)
{
    //Free the edge table
}

void tui_test_edgetables(void)
{
    char* cliptesttree = NULL;
    
    //cliptesttree = "temp_examp6=[&R] ((1,2),(3,4));";
    cliptesttree = "temp_examp6=[&R] ((1,2),3);";
    //cliptesttree = "temp_examp6=[&R] ((1,2),(3,(4,5)));";
    //cliptesttree = "temp_examp6=[&R] (1,(2,(3,(4,(5,6)))));";
    //cliptesttree = "temp_examp6=[&R] (5,(4,(3,(2,1))));";
    //cliptesttree = "temp_examp6=[&R] ((1,(2,(6,7))),(3,(4,5)));";
    
    mfl_tree_t* testree = mfl_convert_newick_to_mfl_tree_t(cliptesttree, 0);
    
    // Allocate the table
    mfl_edgetable_t* test_edgetable = mfl_initiate_edgetable_t(testree->treet_num_taxa);
    
    //mfl_get_edgetable(test_edgetable, testree);
    
    
    
}
