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

void mfl_fitch_downpass_binary_node(mfl_nodedata_t* n_nd, mfl_nodedata_t* left_nd, mfl_nodedata_t* right_nd, mfl_nodedata_t* dummy, mfl_datapartition_t* datapart, int* length)
{
    int i = 0;
    int* weights = datapart->part_int_weights;
    int num_chars = datapart->part_n_chars_included;
    mfl_charstate* n_prelim = n_nd->nd_prelim_set;
    mfl_charstate* left = left_nd->nd_prelim_set;
    mfl_charstate* right = right_nd->nd_prelim_set;
    mfl_charstate temp = 0;
    
    for (i = 0; i < num_chars; ++i) {
        if ((temp = left[i] & right[i])) {
            n_prelim[i] = temp;
        }
        else {
            n_prelim[i] = left[i] | right[i];
            if (length) {
                *length += weights[i];
            }
        }
    }
}

void mfl_fitch_uppass_binary_node(mfl_nodedata_t* n_nd, mfl_nodedata_t* left_nd, mfl_nodedata_t* right_nd, mfl_nodedata_t* anc_nd, mfl_datapartition_t* datapart, int* length)
{
    int i = 0;
    int num_chars = datapart->part_n_chars_included;
    mfl_charstate* lft_char = NULL;
    mfl_charstate* rt_char = NULL;
    mfl_charstate* n_prelim = n_nd->nd_prelim_set;
    mfl_charstate* n_final = n_nd->nd_final_set;
    mfl_charstate* anc_char = anc_nd->nd_final_set;
    
    if (!left_nd) {
        assert(!right_nd);
        for (i = 0; i < num_chars; ++i) {
            
            if (n_prelim[i] & anc_char[i]) {
                n_final[i] = n_prelim[i] & anc_char[i];
            }
            else {
                n_final[i] = n_prelim[i];
            }
        
        }
        
        return;
    }
    
    lft_char = left_nd->nd_prelim_set;
    rt_char = right_nd->nd_prelim_set;
    
    for (i = 0; i < num_chars; ++i) {
        if ((anc_char[i] & n_prelim[i]) == anc_char[i]) {
            n_final[i] = anc_char[i] & n_prelim[i];
#ifdef MFY_DEBUG
            if (length) {
                if (!(lft_char[i] & rt_char[i])) {
                    //TODO: Change this to weights
                    *length = *length + 1;
                }
            }
#endif
        }
        else {
            if (lft_char[i] & rt_char[i]) {
                n_final[i] = ( n_prelim[i] | ( anc_char[i] & (lft_char[i] | rt_char[i])));
            } else {
                n_final[i] = n_prelim[i] | anc_char[i];
#ifdef MFY_DEBUG
                if (length) {
                    //TODO: Change this to weights
                    *length = *length + 1;
                }
#endif
            }
        }
        
        assert(n_final[i]);
    }
}

void mfl_fitch_downpass_inapplicables(mfl_nodedata_t*       n_nd,
                                      mfl_nodedata_t*       left_nd,
                                      mfl_nodedata_t*       right_nd,
                                      mfl_nodedata_t*       dummy,
                                      mfl_datapartition_t*  datapart,
                                      int*                  length)
{
    int i = 0;
    //int* weights = datapart->part_int_weights;
    int num_chars = datapart->part_n_chars_included;
    mfl_charstate* n_prelim = n_nd->nd_prelim_set;
    mfl_charstate* left = left_nd->nd_prelim_set;
    mfl_charstate* right = right_nd->nd_prelim_set;
    mfl_charstate temp = 0;
    //mfl_charstate* actives = datapart->part_activestates;
    
    for (i = 0; i < num_chars; ++i) {

        
        if ((temp =  left[i] & right[i])) {
            
            if ((left[i] & MORPHY_IS_APPLICABLE) && (right[i] & MORPHY_IS_APPLICABLE)) {
                
                if (temp & MORPHY_IS_APPLICABLE) {
                    n_prelim[i] = temp & MORPHY_IS_APPLICABLE;
                }
                else {
                    n_prelim[i] = (left[i] | right[i]) & MORPHY_IS_APPLICABLE;
                }
            }
            else {
                n_prelim[i] = MORPHY_INAPPLICABLE_BITPOS;
            }
            assert(n_prelim[i]);
        }
        else {
            
            n_prelim[i] = (left[i] | right[i]) & MORPHY_IS_APPLICABLE;

            assert(n_prelim[i]);
        }
        
        // Set all states active on this node so that the optimisation shortcut "knows" whether the state occurs upwards of this view.
    }
    
}


void mfl_fitch_uppass_inapplicables(mfl_nodedata_t*       n_nd,
                                    mfl_nodedata_t*       left_nd,
                                    mfl_nodedata_t*       right_nd,
                                    mfl_nodedata_t*       anc_nd,
                                    mfl_datapartition_t*  datapart,
                                    int*                  length)
{

    int* weights = datapart->part_int_weights;
    int i = 0;
    int num_chars = datapart->part_n_chars_included;
    mfl_charstate* lft_char = NULL;
    mfl_charstate* rt_char = NULL;
    mfl_charstate* n_prelim = n_nd->nd_prelim_set;
    mfl_charstate* n_final = n_nd->nd_final_set;
    mfl_charstate* anc_char = anc_nd->nd_final_set;
    mfl_charstate* actives = datapart->part_activestates;
    mfl_charstate temp = 0;
    
    if (!left_nd) {
        assert(!right_nd);
        for (i = 0; i < num_chars; ++i) {
            
            if (n_prelim[i] & anc_char[i]) {
                n_final[i] = n_prelim[i] & anc_char[i];
                assert(n_final[i]);
            }
            else {
                n_final[i] = n_prelim[i];
                assert(n_final[i]);
            }
            
            if (n_prelim[i] == 8) {
                dbg_printf("break\n");
            }
            
            if (length) {
                if ((n_final[i] & anc_char[i]) != anc_char[i]) {
                    if (n_final[i] == (n_final[i] & -n_final[i])) {
                        if (n_final[i] != MORPHY_MISSING_DATA_BITWISE) {
                            if (n_final[i] & actives[i]) {
                                *length += weights[i];
                                dbg_printf("n_prelim[i]: %llu\n", n_prelim[i]);
                                dbg_printf("n_final[i]:  %llu\n", n_final[i]);
                                dbg_printf("anc_char[i]: %llu\n\n", anc_char[i]);
                            }
                        }
                    }
                }
            }
            
            if (n_prelim[i] != MORPHY_MISSING_DATA_BITWISE) {
                actives[i] = actives[i] | (n_final[i] & MORPHY_IS_APPLICABLE);
            }
        }
        
        return;
    }
    
    lft_char = left_nd->nd_prelim_set;
    rt_char = right_nd->nd_prelim_set;
    
    assert(lft_char);
    assert(rt_char);
    
    for (i = 0; i < num_chars; ++i) {
        
        if (anc_char[i] == MORPHY_INAPPLICABLE_BITPOS) {
            if ((lft_char[i] | rt_char[i]) & MORPHY_INAPPLICABLE_BITPOS) {
                if ((lft_char[i] & MORPHY_IS_APPLICABLE) && (rt_char[i] & MORPHY_IS_APPLICABLE)) {
                    // If the intersection of applicable states exists, create it. Otherwise, union.
                    if (((temp = (lft_char[i] & MORPHY_IS_APPLICABLE) & (rt_char[i] & MORPHY_IS_APPLICABLE)))) {
                        n_final[i] = temp;
                    } else {
                        n_final[i] = (lft_char[i] | rt_char[i]) & MORPHY_IS_APPLICABLE;
                    }
                }
                else {
                    n_final[i] = MORPHY_INAPPLICABLE_BITPOS;
                }
            }
            else {
                n_final[i] = n_prelim[i];// & MORPHY_IS_APPLICABLE;
                assert(!(n_final[i] & MORPHY_INAPPLICABLE_BITPOS));
            }
        }
        else {
            if ((anc_char[i] & n_prelim[i]) == anc_char[i]) {
                n_final[i] = anc_char[i] & n_prelim[i];
                assert(n_final[i]);
            }
            else {
                if ((lft_char[i] & rt_char[i]) & MORPHY_IS_APPLICABLE) {
                    
                    if ((temp = (lft_char[i] | rt_char[i]) & anc_char[i])) {
                        n_final[i] = temp & MORPHY_IS_APPLICABLE;
                    }
                    else {
                        n_final[i] = (lft_char[i] | rt_char[i]) & MORPHY_IS_APPLICABLE;
                    }
                    assert(n_final[i]);
                } else {
                    if ((lft_char[i] | rt_char[i]) == MORPHY_INAPPLICABLE_BITPOS) {
                        n_final[i] = MORPHY_INAPPLICABLE_BITPOS;
                    }
                    else {
                            
                            if ((temp = (lft_char[i] | rt_char[i]) & anc_char[i])) {
                                n_final[i] = temp;
                            }
                            else {
                                n_final[i] = (lft_char[i] | rt_char[i]) & MORPHY_IS_APPLICABLE;
                            }
                        
                    }
                    assert(n_final[i]);
                }
            }
        }
        
        
        if (anc_char[i] != MORPHY_INAPPLICABLE_BITPOS) {
            //if (n_final[i] == (n_final[i] & anc_char[i])) {
                if (n_final[i] != anc_char[i]) {
                    if (n_final[i] != (n_final[i] & -n_final[i])) {
                        n_final[i] = n_final[i] ^ (n_final[i] & -n_final[i]);
                    }
                }
            //}
        }
        
        int x = n_final[i];
        x = x & -x;
        
        if (length) {
        
            if (!(n_final[i] & anc_char[i]) /*!= anc_char[i]*/) {
                if (n_final[i] & actives[i]) {
                    *length += weights[i];
                    dbg_printf("n_prelim[i]: %llu\n", n_prelim[i]);
                    dbg_printf("n_final[i]:  %llu\n", n_final[i]);
                    dbg_printf("lft_char[i]: %llu\n", lft_char[i]);
                    dbg_printf("rt_char[i]:  %llu\n", rt_char[i]);
                    dbg_printf("anc_char[i]: %llu\n\n", anc_char[i]);
                    
                }
            }
        }
        
        if (n_final[i] == x ) {
            if (n_final[i] != MORPHY_MISSING_DATA_BITWISE) {
                actives[i] = actives[i] | (n_final[i] & MORPHY_IS_APPLICABLE);
            }
        }
        assert(n_final[i]);
    }
}

void mfl_set_rootstates(mfl_node_t* dummyroot, mfl_node_t* rootnode, mfl_partition_set_t* dataparts)
{
    
    int i = 0;
    int j = 0;
    
    // TODO: This function can be optimised considerably by diminishing the amount of repeated indirection calls, but test it first.
    
    for (i = 0; i < dataparts->ptset_n_parts; ++i) {
        if (dataparts->ptset_partitions[i]->part_has_inapplicables) {
            // Do inapplicable root optimisation
            for (j = 0; j < dataparts->ptset_partitions[i]->part_n_chars_included; ++j) {
                if (rootnode->nodet_charstates[i]->nd_prelim_set[j] & MORPHY_IS_APPLICABLE) {
                    dummyroot->nodet_charstates[i]->nd_final_set[j] = rootnode->nodet_charstates[i]->nd_prelim_set[j] & MORPHY_IS_APPLICABLE;
                }
                else {
                    dummyroot->nodet_charstates[i]->nd_final_set[j] = MORPHY_INAPPLICABLE_BITPOS;
                }
            }
            
        }
        else {
            // Do regular root optimisation
            for (j = 0; j < dataparts->ptset_partitions[i]->part_n_chars_included; ++j) {
                dummyroot->nodet_charstates[i]->nd_final_set[j] = rootnode->nodet_charstates[i]->nd_prelim_set[j];
            }
        }
    }
    
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


void mfl_postorder_traversal(mfl_node_t *n, int* length)
{
    int i = 0;
    int num_dataparts;
    mfl_node_t *p = NULL;
    mfl_parsim_fn evaluator;
    mfl_node_t* left;
    mfl_node_t* right;
    
    if (n->nodet_tip) {
        return;
    }
    
    left = n->nodet_next->nodet_edge;
    right = n->nodet_next->nodet_next->nodet_edge;
    
    p = n->nodet_next;
    do {
        mfl_postorder_traversal(p->nodet_edge, length);
        p = p->nodet_next;
    } while (p != n);
    
    num_dataparts = n->nodet_num_dat_partitions;
    
    // For each data partition at the node, set the correct type and evaluation
    for (i = 0; i < num_dataparts; ++i) {
        evaluator = n->nodet_charstates[i]->nd_downpass_full;
        evaluator(
                  n->nodet_charstates[i],
                  left->nodet_charstates[i],
                  right->nodet_charstates[i],
                  NULL,
                  n->nodet_charstates[i]->nd_parent_partition,
                  length
                  );
    }
    
    return;
}


void mfl_preorder_traversal(mfl_node_t *n, int* length)
{
    int i = 0;
    int num_dataparts;
    mfl_node_t *p = NULL;
    mfl_parsim_fn evaluator;
    mfl_node_t* left;
    mfl_node_t* right;
    
    num_dataparts = n->nodet_num_dat_partitions;
    
    if (!n->nodet_tip) {
        left = n->nodet_next->nodet_edge;
        right = n->nodet_next->nodet_next->nodet_edge;
    }
    
    // For each data partition at the node, set the correct type and evaluation
    for (i = 0; i < num_dataparts; ++i) {
        evaluator = n->nodet_charstates[i]->nd_uppass_full;
        
        if (n->nodet_tip) {
            evaluator(
                      n->nodet_charstates[i],
                      NULL,
                      NULL,
                      n->nodet_edge->nodet_charstates[i],
                      n->nodet_charstates[i]->nd_parent_partition,
                      length
                      );
        }
        else {
            if (n->nodet_next->nodet_edge->nodet_tip == 15) {
                dbg_printf("break\n");
            }
            evaluator(
                      n->nodet_charstates[i],
                      left->nodet_charstates[i],
                      right->nodet_charstates[i],
                      n->nodet_edge->nodet_charstates[i],
                      n->nodet_charstates[i]->nd_parent_partition,
                      length
                      );
        }
        
        
    }
    
    if (n->nodet_tip) {
        return;
    }
    
    p = n->nodet_next;
    do {
        mfl_preorder_traversal(p->nodet_edge, length);
        p = p->nodet_next;
    } while (p != n);
    
    return;
}

void mfl_fullpass_tree_optimisation(mfl_tree_t* t, mfl_partition_set_t* dataparts)
{
    // TODO: finish this function
    // check is rooted; if not do the fucking root;
    
    mfl_postorder_traversal(t->treet_root, &t->treet_parsimonylength);
    dbg_printf("\nHere's the length after the downpass: %i\n", t->treet_parsimonylength);
    tui_check_broken_tree(t, false);
    
    mfl_set_rootstates(&t->treet_dummynode, t->treet_root, dataparts);
    tui_check_broken_tree(t, false);
    
    dbg_printf("\nResetting length to 0 before uppass\n");
    t->treet_parsimonylength = 0;
    
    tui_check_broken_tree(t, false);
    mfl_preorder_traversal(t->treet_root, &t->treet_parsimonylength);
    dbg_printf("\nHere's the length after the uppass: %i\n", t->treet_parsimonylength);
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