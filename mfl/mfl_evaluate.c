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

void mfl_postorder_traversal(mfl_node_t *parent, mfl_searchrec_t *search_rec)
{
    mfl_node_t *p = NULL;
    
    if (parent->nodet_tip) {
        return;
    }
    
    p = parent;
    do {
        p = parent->nodet_next;
        mfl_postorder_traversal(p->nodet_edge, search_rec);
    } while (p != parent);
    
    // Some function activity here
    // We could put a function pointer here for ancestral state reconstructions
    // The reconstruction algorithms would all have to take the same number and
    // and type of argument, which probably wouldn't be too difficult to do.
    // Function pointers are less efficient than a direct function call and
    // even less efficient than an inline function call. 
    
    return;
}

void mfl_preorder_traversal(mfl_node_t *parent, mfl_searchrec_t *search_rec)
{
    mfl_node_t *p = NULL;
    
    if (parent->nodet_tip) {
        return;
    }
    
    // Some function activity here
    // Possibly follow same procedure as in mfl_postorder_traversal().
    
    p = parent;
    do {
        p = parent->nodet_next;
        mfl_preorder_traversal(p->nodet_edge, search_rec);
    } while (p != parent);
    
    return;
}