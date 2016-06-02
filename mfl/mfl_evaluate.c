/*
 *  mfl_evaluate.c
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


//long double mfl_get_transformation_cost(mfl_costs_t weights, int from_state, int to_state)
//{
//    /* Just a 'prototype' for a crude stepmatrix calculation for rooted 
//     * characters. I suspect there's a cleverer way to do this. */
//    
//    int from_index = 0;
//    int to_index = 0;
//    long double return_weight = 0.0;
//    
//    from_index = from_state - 1;
//    to_index = to_state - 1;
//    
//    return_weight = *(weights + from_index + to_index);
//    
//    return return_weight;
//}

void mfl_fitch_downpass_binary_node(mfl_charstate* anc, mfl_charstate* left, mfl_charstate* right, mfl_datapartition_t* datapart, int* length)
{
    int i = 0;
    int* weights = datapart->part_int_weights;
    int num_chars = datapart->part_n_chars_included;
    mfl_charstate temp = 0;
    
    for (i = 0; i < num_chars; ++i) {
        if ((temp = left[i] & right[i])) {
            anc[i] = temp;
        }
        else {
            anc[i] = left[i] | right[i];
            if (length) {
                *length += weights[i];
            }
        }
    }
    
}

void mfl_fitch_downpass_INAPPLICABLE(mfl_node_t *node)
{
    int i = 0;
    mfl_node_t* lchild = NULL;
    mfl_node_t* rchild = NULL;
    lchild = node->nodet_next->nodet_edge;
    rchild = node->nodet_next->nodet_next->nodet_edge;
    mfl_charstate temp = NULL;
    
    //    mfl_charstate* parentchars = node->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
    //    mfl_charstate* leftchars   = lchild->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
    //    mfl_charstate* rightchars  = rchild->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
    //    int num_chars = node->nodet_dataparts[MFL_OPT_FITCH]->nd_n_characters;
    //
    //    for (i = 0; i < num_chars; ++i) {
    //        if ((temp = leftchars[i] & rightchars[i])) {
    //            parentchars[i] = temp;
    //        }
    //        else {
    //            parentchars[i] = leftchars[i] | rightchars[i];
    //            /* Increment the length of the tree, and maxsteps*/
    //        }
    //    }
}

void mfl_fitch_uppass_binary_node(mfl_node_t *node)
{
    int i = 0;
    mfl_node_t* lchild   = NULL;
    mfl_node_t* rchild   = NULL;
    mfl_node_t* ancestor = NULL;
    lchild = node->nodet_next->nodet_edge;
    rchild = node->nodet_next->nodet_next->nodet_edge;
    mfl_charstate temp   = NULL;
//    mfl_charstate* parent_prelim = node->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
//    mfl_charstate* parent_final  = node->nodet_dataparts[MFL_OPT_FITCH]->nd_final_set;
//    mfl_charstate* leftchars     = lchild->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
//    mfl_charstate* rightchars    = rchild->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
//    mfl_charstate* ancchars      = ancestor->nodet_dataparts[MFL_OPT_FITCH]->nd_final_set;
//    int num_chars = node->nodet_dataparts[MFL_OPT_FITCH]->nd_n_characters;
    
//    for (i = 0; i < num_chars; ++i) {
//        // Uppass business logic.
//    }
}

inline int mfl_wagner_stepcount(mfl_charstate leftchar, mfl_charstate rightchar, mfl_charstate* parentchar, int weight)
{
    /* Calculates the number of steps between non-overlapping state sets from two 
     * descendant branches at a binary node. There might be a better way to do this
     * and it might have unpredictable behaviour if the user supplies strange terminal
     * state sets like D1 = 0100010 and D2 = 0000101. For now, this should do for most
     * normal cases. */
    
    int length_increment = 0;
    mfl_charstate newset = 0;
    mfl_charstate big = 0;
    mfl_charstate small = 0;
    
    
    big = MAX(leftchar, rightchar);
    small = MIN(leftchar, rightchar);
    
    do {
        ++length_increment;
        newset = (big & (small << length_increment));
    } while (!newset);
    
    // Close the set between descendant sets
    do {
        newset = newset | (newset << 1);
    } while (!(big & newset));
    
    // Assign this new set to the parent set
    if (parentchar) {
        *parentchar = newset;
    }
    
    return length_increment * weight;
}

void mfl_wagner_downpass_binary_node(mfl_node_t *node)
{
    int i = 0;
    mfl_node_t* lchild = NULL;
    mfl_node_t* rchild = NULL;
    lchild = node->nodet_next->nodet_edge;
    rchild = node->nodet_next->nodet_next->nodet_edge;
    mfl_charstate temp = NULL;
    mfl_charstate* parentchars;// = node->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
    mfl_charstate* leftchars;//   = lchild->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
    mfl_charstate* rightchars;//  = rchild->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
    int num_chars;// = node->nodet_dataparts[MFL_OPT_FITCH]->nd_n_characters;
    
    for (i = 0; i < num_chars; ++i) {
        if ((temp = leftchars[i] & rightchars[i])) {
            
        }
        else {
            //parentchars[i] = leftchars[i] | rightchars[i];
            /*Lenght increase = */ mfl_wagner_stepcount(leftchars[i], rightchars[i], &parentchars[i],  NULL/* WEIGHT goes here*/);
        }
    }
}


void mfl_postorder_traversal(mfl_node_t *n, mfl_searchrec_t *search_rec)
{
    int i = 0;
    int num_dataparts;
    mfl_node_t *p = NULL;
    mfl_parsim_fn evaluator;
    mfl_node_t* left = n->nodet_next->nodet_edge;
    mfl_node_t* right = n->nodet_next->nodet_next->nodet_edge;
    
    if (n->nodet_tip) {
        return;
    }
    
    p = n;
    do {
        p = p->nodet_next;
        mfl_postorder_traversal(p->nodet_edge, search_rec);
    } while (p != n);
    
    num_dataparts = n->nodet_num_dat_partitions;
    
    // For each data partition at the node, set the correct type and evaluation
    for (i = 0; i < num_dataparts; ++i) {
        evaluator = n->nodet_charstates[i]->nd_downpass_full;
        evaluator(
                  n->nodet_charstates[i]->nd_prelim_set,
                  left->nodet_charstates[i]->nd_prelim_set,
                  right->nodet_charstates[i]->nd_prelim_set,
                  n->nodet_charstates[i]->nd_parent_partition,
                  &search_rec->sr_swaping_on->treet_parsimonylength
                  );
    }
    
    return;
}

void mfl_postorder_traversal_full(mfl_node_t *parent, mfl_searchrec_t *search_rec)
{
    int i = 0;
    int num_dataparts;
    mfl_node_t *p = NULL;
    mfl_parsim_fn evaluator;
    
    if (parent->nodet_tip) {
        return;
    }
    
    p = parent;
    do {
        p = p->nodet_next;
        mfl_postorder_traversal_full(p->nodet_edge, search_rec);
    } while (p != parent);
    
    num_dataparts = parent->nodet_num_dat_partitions;
    
    // For each data partition at the node, set the correct type and evaluation
    for (i = 0; i < num_dataparts; ++i) {
        evaluator = parent->nodet_charstates[i]->nd_downpass_full;
        //evaluator(parent);
    }
    
    return;
}


void mfl_preorder_traversal_partial(mfl_node_t *parent, mfl_searchrec_t *search_rec)
{
    mfl_node_t *p = NULL;
    
    if (parent->nodet_tip) {
        // Set uppass set for the tips
        return;
    }
    
    // Some function activity here
    // Possibly follow same procedure as in mfl_postorder_traversal().
    
    p = parent;
    do {
        p = p->nodet_next;
        mfl_preorder_traversal_partial(p->nodet_edge, search_rec);
    } while (p != parent);
    
    return;
}


int mfl_unordered_distance(mfl_charstate* t, const mfl_charstate* a, int* weights, int num_chars)
{
    int i = 0;
    int d = 0;
    
    for (i = 0; i < num_chars; ++i) {
        if (!(t[i] & a[i])) {
            d += weights[i];
        }
    }
}



int mfl_ordered_distance(mfl_charstate* t, const mfl_charstate* a, int* weights, int num_chars)
{
    int i = 0;
    int d = 0;
    
    for (i = 0; i < num_chars; ++i) {
        if (!(t[i] & a[i])) {
            d = mfl_wagner_stepcount(t[i], a[i], NULL, weights[i]);
        }
    }
}


/*!
 @discussion Defines if final set is different than preliminary set
 @param target_node (mfl_nodedata_t*) a target node
 @return a boolean, whether the final set is different (TRUE) or not (FALSE)
 */
bool mfl_is_final_different(mfl_nodedata_t *target_node)
{
    bool is_node_different = 1;     // Probably safer to initialise the bool as being different
    
    // Check if the prelim set is equal to the final set
    if(target_node->nd_prelim_set == target_node->nd_final_set){
        // Sets are not different
        is_node_different = 0;
        return is_node_different;
    }
    //Sets are different
    return is_node_different;
}

/*
//TG: next step could be to add that in the traversal and store the results in a array of bools:

 void mfl_my_traversal(bool *node_differences, int *count)
 {
    //...
    node_differences[*count] = mfl_is_final_different(target_node);
    //...
 }
 
//TG: and then use this in a character_to_change array
*/