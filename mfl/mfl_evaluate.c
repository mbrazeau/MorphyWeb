//
//  mfl_evaluate.c
//  Morphy
//
//  Created by mbrazeau on 09/03/2016.
//  Copyright © 2016 Imperial College London. All rights reserved.
//

#include "morphy.h"


long double mfl_get_transformation_cost(mfl_costs_t weights, int from_state, int to_state)
{
    /* Just a 'prototype' for a crude stepmatrix calculation for rooted 
     * characters. I suspect there's a cleverer way to do this. */
    
    int from_index = 0;
    int to_index = 0;
    long double return_weight = 0.0;
    
    from_index = from_state - 1;
    to_index = to_state - 1;
    
    return_weight = *(weights + from_index + to_index);
    
    return return_weight;
}

void mfl_fitch_downpass_binary_node(mfl_node_t *node)
{
    int i = 0;
    mfl_node_t* lchild = NULL;
    mfl_node_t* rchild = NULL;
    lchild = node->nodet_next->nodet_edge;
    rchild = node->nodet_next->nodet_next->nodet_edge;
    mfl_charstate temp = NULL;
    mfl_charstate* parentchars = node->nodet_dataparts[MFL_IS_FITCH]->nd_prelim_set;
    mfl_charstate* leftchars   = lchild->nodet_dataparts[MFL_IS_FITCH]->nd_prelim_set;
    mfl_charstate* rightchars  = rchild->nodet_dataparts[MFL_IS_FITCH]->nd_prelim_set;
    int num_chars = node->nodet_dataparts[MFL_IS_FITCH]->nd_n_characters;
    
    for (i = 0; i < num_chars; ++i) {
        if ((temp = leftchars[i] & rightchars[i])) {
            parentchars[i] = temp;
        }
        else {
            parentchars[i] = leftchars[i] | rightchars[i];
            /* Increment the length of the tree, and maxsteps*/
        }
    }
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
    mfl_charstate* parent_prelim = node->nodet_dataparts[MFL_IS_FITCH]->nd_prelim_set;
    mfl_charstate* parent_final  = node->nodet_dataparts[MFL_IS_FITCH]->nd_final_set;
    mfl_charstate* leftchars     = lchild->nodet_dataparts[MFL_IS_FITCH]->nd_prelim_set;
    mfl_charstate* rightchars    = rchild->nodet_dataparts[MFL_IS_FITCH]->nd_prelim_set;
    mfl_charstate* ancchars      = ancestor->nodet_dataparts[MFL_IS_FITCH]->nd_final_set;
    int num_chars = node->nodet_dataparts[MFL_IS_FITCH]->nd_n_characters;
    
    for (i = 0; i < num_chars; ++i) {
        // Uppass business logic.
    }
}

mfl_parsim_fn mfl_fetch_downpass_parsimony_fxn(mfl_optimisation_t parsim_type)
{
    mfl_parsim_fn ret = NULL;
    
    switch (parsim_type)
    {
        case MFL_IS_FITCH:
            ret = mfl_fitch_downpass_binary_node;
            break;
        
        /*case MFL_IS_WAGNER:
            ret = mfl_wagner_downpass_binary_node;
            break;
        
        case MFL_IS_DOLLO:
            ret = mfl_dollo_downpass_binary_node;
            break;
        
        case MFL_IS_IRREVERSIBLE:
            ret = mfl_irreversible_downpass_binary_node;
            break;
            
        case MFL_IS_COST_MATRIX:
            ret = mfl_costmatrix_downpass_binary_node;
            break;*/
            
        default:
            break;
    }
    
    return  ret;
}

mfl_parsim_fn mfl_fetch_uppass_parsimony_fxn(mfl_optimisation_t parsim_type)
{
    mfl_parsim_fn ret = NULL;
    
    switch (parsim_type)
    {
        /*case MFL_IS_FITCH:
        ret = mfl_fitch_uppass_binary_node;
        break;*/
            
        /*case MFL_IS_WAGNER:
        ret = mfl_wagner_downpass_binary_node;
        break;*/
             
        /*case MFL_IS_DOLLO:
        ret = mfl_dollo_downpass_binary_node;
        break;*/
             
        /*case MFL_IS_IRREVERSIBLE:
        ret = mfl_irreversible_downpass_binary_node;
        break;*/
             
        /*case MFL_IS_COST_MATRIX:
        ret = mfl_costmatrix_downpass_binary_node;
        break;*/
            
        default:
            break;
    }
    
    return  ret;
}


void mfl_evaluate_downpass(mfl_node_t *node)
{
    int i = 0;
    int num_dataparts;
    mfl_parsim_fn evaluator;
    
    num_dataparts = node->nodet_num_partitions;
    
    // For each data partition at the node, set the correct type and evaluation
    for (i = 0; i < num_dataparts; ++i) {
        evaluator = mfl_fetch_downpass_parsimony_fxn(node->nodet_dataparts[i]->nd_optimisation_method);
        evaluator(node);
    }
    
    return;
}

void mfl_postorder_traversal(mfl_node_t *parent, mfl_searchrec_t *search_rec)
{
    
    mfl_node_t *p = NULL;
    
    if (parent->nodet_tip) {
        return;
    }
    
    p = parent;
    do {
        p = p->nodet_next;
        mfl_postorder_traversal(p->nodet_edge, search_rec);
    } while (p != parent);
    
    mfl_evaluate_downpass(parent);
    
    return;
}

void mfl_preorder_traversal(mfl_node_t *parent, mfl_searchrec_t *search_rec)
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
        mfl_preorder_traversal(p->nodet_edge, search_rec);
    } while (p != parent);
    
    return;
}