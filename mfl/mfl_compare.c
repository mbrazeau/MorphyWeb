//
//  mfl_compare.c
//  Morphy
//
//  Created by mbrazeau on 05/05/2016.
//
//

#include "morphy.h"

/* This is mostly a temporary function, as it is likely that the setting of 
 * bipartitions will be handled simultaneously by other tree traversals. 
 */
void mfl_set_bipartitions(mfl_node_t* n)
{
    assert(n->nodet_bipart->bts_bitfields);
    
    mfl_node_t* p = NULL;
    
    if (n->nodet_tip) {
        return;
    }
    
    p = n->nodet_next;
    do {
        mfl_set_bipartitions(p->nodet_edge);
        mfl_bts_OR(p->nodet_edge->nodet_bipart, n->nodet_bipart, n->nodet_bipart);
        p = p->nodet_next;
    } while (p != n);
    
}
