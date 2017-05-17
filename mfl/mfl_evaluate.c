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

//int mfl_locreopt_cost(node *src, node *tgt1, node *tgt2, int nchar, int diff)
//{
//    /* Returns cost of inserting subtree src between tgt1 and tgt2 following
//     * the algorithms described by Ronquist (1998. Cladistics) and Goloboff
//     * (1993, 1996. Cladistics).*/
//    
//    int i;
//    int cost = 0;
//    charstate *srctemps = src->apomorphies;
//    charstate *tgt1apos = tgt1->apomorphies;
//    charstate *tgt2apos = tgt2->apomorphies;
//    
//    
//    for (i = 0; i < nchar; ++i, ++srctemps, ++tgt1apos, ++tgt2apos) {
//        if ( !(*srctemps & (*tgt1apos | *tgt2apos)) ) {
//            ++cost;
//            if (cost > diff) {
//                return cost;
//            }
//        }
//    }
//    
//    return cost;
//}

int mfl_test_fitch_local(const mfl_nodedata_t* src_nd,
                         const mfl_nodedata_t* tgt1_nd,
                         const mfl_nodedata_t* tgt2_nd,
                         mfl_datapartition_t* dataprt,
                         const int diff)
{
    int i = 0;
    int cost = 0;
    int *weights = dataprt->part_int_weights;
    int num_chars = dataprt->part_n_chars_included;
    mfl_charstate* src = src_nd->nd_prelim_set;
    mfl_charstate* tgt1 = tgt1_nd->nd_final_set;
    mfl_charstate* tgt2 = tgt2_nd->nd_final_set;
    
    // TODO: Optimise: increment pointers in loop head
    for (i = 0; i < num_chars; ++i) {
        if (!(src[i] & (tgt1[i] | tgt2[i]))) {
            cost += weights[i];
            // TODO: The outer condition permits writing one local check
            // TODO: function. It can be removed in a version of this function
            // TODO: that is used specifically during the search
//            if (!(diff < 0)) {
//                if (cost > diff) {
//                    return cost;
//                }
//            }
        }
    }
    
    return cost;
}


void mfl_add_to_change_list(int index, mfl_datapartition_t *dataprt)
{
    dataprt->part_char_changing[dataprt->nchanges] = index;
    ++dataprt->nchanges;
    assert(dataprt->nchanges <= dataprt->part_n_chars_included);
}


void mfl_reset_change_list(mfl_datapartition_t *dataprt)
{
    dataprt->nchanges = 0;
    memset(dataprt->part_char_changing, 0, dataprt->part_n_chars_max * sizeof(int));
}


int mfl_test_fitch_na_local(const mfl_nodedata_t* src_nd,
                            const mfl_nodedata_t* tgt1_nd,
                            const mfl_nodedata_t* tgt2_nd,
                            mfl_datapartition_t* dataprt,
                            const int diff)
{
    int i       = 0;
    int cost    = 0;
    int regions = 0;
    int *weights  = dataprt->part_int_weights;
    int num_chars = dataprt->part_n_chars_included;
    mfl_charstate* src   = src_nd->nd_prelim_set;
    mfl_charstate* tgt1f = tgt1_nd->nd_final_set;
    mfl_charstate* tgt2f = tgt2_nd->nd_final_set;
    mfl_charstate* tgt1p = tgt1_nd->nd_initprelim;
    mfl_charstate* tgt2p = tgt2_nd->nd_initprelim;
    mfl_charstate* tgt1a = tgt1_nd->nd_subtree_activestates;
    mfl_charstate* tgt2a = tgt2_nd->nd_subtree_activestates;
    
    
    mfl_reset_change_list(dataprt);
    
    // TODO: Optimise: increment pointers in loop head
    for (i = 0; i < num_chars; ++i) {
        
        if (!(src[i] & (tgt1f[i] | tgt2f[i]))) {
            if (src[i] & MORPHY_IS_APPLICABLE) {
                if ((tgt1f[i] | tgt2f[i]) & MORPHY_IS_APPLICABLE) {
                    cost += weights[i];
                    //  if (!(diff < 0)) {
                    //      if (cost > diff) {
                    //          return cost + regions;
                    //      }
                    //  }
                }
                else if ((tgt1p[i] | tgt2p[i]) & MORPHY_IS_APPLICABLE) {
                    // That's complicated!
                    // Push i to list; We'll do that one later!
                    mfl_add_to_change_list(i, dataprt);
                }
                else if (tgt1a[i] || tgt2a[i]) {
                    regions += weights[i];
                }
            }
            else {
                if (tgt1p[i] & tgt2p[i] & MORPHY_INAPPLICABLE_BITPOS) {
                    if (tgt1p[i] & tgt2p[i] & MORPHY_IS_APPLICABLE) {
                        regions += weights[i];
                    }
                }
                else if (tgt1p[i] == MORPHY_INAPPLICABLE_BITPOS || tgt2p[i] == MORPHY_INAPPLICABLE_BITPOS) {
                    if (tgt1a[i] && tgt2a[i]) {
                        if (tgt1f[i] == tgt2f[i]) {
                            regions += weights[i];
                        }
                    }
                }
            }
        }
    }
    
    return cost + regions;
}


void mfl_local_add_cost(mfl_node_t* src, mfl_node_t* tgt, const int diff, int *cost)
{
    int i = 0;
    int num_dataparts = 0;
    mfl_lparsim_fn evaluator;
    mfl_node_t *tgtopp = tgt->nodet_edge;
    
    num_dataparts = src->nodet_num_dat_partitions;
    assert(src->nodet_num_dat_partitions == tgt->nodet_num_dat_partitions);
    
    for (i = 0; i < num_dataparts; ++i) {
        
        evaluator = src->nodet_charstates[i]->nd_local;
        
        *cost += evaluator((const mfl_nodedata_t*)src->nodet_charstates[i],
                           (const mfl_nodedata_t*)tgt->nodet_charstates[i],
                           (const mfl_nodedata_t*)tgtopp->nodet_charstates[i],
                           src->nodet_charstates[i]->nd_parent_partition,
                           diff);
        
        if (src->nodet_charstates[i]->nd_parent_partition->nchanges) {
            dbg_printf("Doing it the hard way for ");
            dbg_printf("%i characters\n",
                       src->nodet_charstates[i]->nd_parent_partition->nchanges);
            
            
            // Insert branch for real.
            mfl_cliprec_t clip;
            mfl_insert_branch_with_ring_base(src, tgt);
            src->nodet_edge->nodet_weight = 3; // TODO: Fix this crap
            
            // Update the subtree state sets.

            // Set the NEW internal node to "marked"
            
            // Plumb down to the last node affected by the insertion
            mfl_node_t* p = tgt->nodet_edge;
            tgtopp->nodet_edge->nodet_isbottom = true;
//            assert(!tgtopp->nodet_isbottom);
            if (!p->nodet_isbottom && !p->nodet_tip) {
                do {
                    p = p->nodet_next;
                } while (!p->nodet_isbottom);
            }
            mfl_update_postorder(p, NULL);
            // Pass back up.
            
            // Then down again.
            
            // Then up again.
            
            // Restore everything to the way it was before
            src->nodet_charstates[i]->nd_parent_partition->nchanges = 0;
//            mfl_reset_change_list(src->nodet_charstates[i]->nd_parent_partition);
            
            // Remove the branch
            tgtopp->nodet_edge->nodet_isbottom = false;
            mfl_clip_branch(src, &clip);
        }
        else {
            dbg_printf("0\n");
        }
    }
}

void mfl_fitch_downpass_binary_node(mfl_nodedata_t* n_nd,
                                    mfl_nodedata_t* left_nd,
                                    mfl_nodedata_t* right_nd,
                                    mfl_nodedata_t* dummy,
                                    mfl_datapartition_t* datapart,
                                    int* length)
{
    int i = 0;
    int* weights = datapart->part_int_weights;
    int num_chars = datapart->part_n_chars_included;
    mfl_charstate* n_prelim = n_nd->nd_prelim_set;
    //mfl_charstate* st_prelim = n_nd->nd_subtree_prelim_set;
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
        
        //st_prelim[i] = n_prelim[i];
    }
}

void mfl_fitch_uppass_binary_node(mfl_nodedata_t* n_nd,
                                  mfl_nodedata_t* left_nd,
                                  mfl_nodedata_t* right_nd,
                                  mfl_nodedata_t* anc_nd,
                                  mfl_datapartition_t* datapart,
                                  int* length)
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
            assert(n_final[i]);
        }
        
        return;
    }
    
    lft_char = left_nd->nd_prelim_set;
    rt_char = right_nd->nd_prelim_set;
    
    for (i = 0; i < num_chars; ++i) {
        if ((anc_char[i] & n_prelim[i]) == anc_char[i]) {
            n_final[i] = anc_char[i] & n_prelim[i];
        }
        else {
            if (lft_char[i] & rt_char[i]) {
                n_final[i] = ( n_prelim[i] | ( anc_char[i] & (lft_char[i] | rt_char[i])));
            } else {
                n_final[i] = n_prelim[i] | anc_char[i];
            }
        }
        
        assert(n_final[i]);
    }
}

void mfl_first_fitch_na_downpass(mfl_nodedata_t*       n_nd,
                                 mfl_nodedata_t*       left_nd,
                                 mfl_nodedata_t*       right_nd,
                                 mfl_nodedata_t*       dummy,
                                 mfl_datapartition_t*  datapart,
                                 int*                  length)
{
    int i = 0;
    int num_chars = datapart->part_n_chars_included;
    mfl_charstate* left = left_nd->nd_initprelim;
    mfl_charstate* right = right_nd->nd_initprelim;
    mfl_charstate* n_initprelim = n_nd->nd_initprelim;
    mfl_charstate* lft_active = left_nd->nd_subtree_activestates;
    mfl_charstate* rt_active = right_nd->nd_subtree_activestates;
    mfl_charstate* subtreeactive = n_nd->nd_subtree_activestates;
    mfl_charstate temp;
    
    for (i = 0; i < num_chars; ++i) {

        temp = 0;
        
        if ((temp = (left[i] & right[i])) ) {
            
            n_initprelim[i] = temp;
            
            if (n_initprelim[i] == MORPHY_INAPPLICABLE_BITPOS) {
                
                if ((left[i] & MORPHY_IS_APPLICABLE) && (right[i] & MORPHY_IS_APPLICABLE)) {
                    n_initprelim[i] = left[i] | right[i];
                }
            }

        }
        else {
            
            n_initprelim[i] = left[i] | right[i];
            
            if ((left[i] & MORPHY_IS_APPLICABLE) && (right[i] & MORPHY_IS_APPLICABLE)) {
                n_initprelim[i] = n_initprelim[i] & MORPHY_IS_APPLICABLE;
            }
            
        }
        
        subtreeactive[i] = (lft_active[i] | rt_active[i]) & MORPHY_IS_APPLICABLE;
        
        assert(n_initprelim[i]);
    }
    
}


void mfl_first_fitch_na_uppass(mfl_nodedata_t*       n_nd,
                               mfl_nodedata_t*       left_nd,
                               mfl_nodedata_t*       right_nd,
                               mfl_nodedata_t*       anc_nd,
                               mfl_datapartition_t*  datapart,
                               int*                  length)
{

    int i = 0;
    int num_chars = datapart->part_n_chars_included;
    mfl_charstate* lft_char = NULL;
    mfl_charstate* rt_char = NULL;
    mfl_charstate* n_initprelim = n_nd->nd_initprelim;
    mfl_charstate* n_initfinal = n_nd->nd_initfinal;
    mfl_charstate* n_prelim = n_nd->nd_prelim_set;
    mfl_charstate* anc_char = anc_nd->nd_initfinal;
    mfl_charstate* subtreeactive = n_nd->nd_subtree_activestates;

    if (!left_nd) {
        assert(!right_nd);
        for (i = 0; i < num_chars; ++i) {
            
            if (n_initprelim[i] & anc_char[i]) {
                subtreeactive[i] = (n_initprelim[i] & anc_char[i]) & MORPHY_IS_APPLICABLE;
            }
            else {
                subtreeactive[i] |= n_initprelim[i] & MORPHY_IS_APPLICABLE;
            }
            
            n_initfinal[i] = n_initprelim[i];
            
            if (n_initfinal[i] & anc_char[i]) {
                if (anc_char[i] & MORPHY_IS_APPLICABLE) {
                    n_initfinal[i] &= MORPHY_IS_APPLICABLE;
                }
            }
            
            n_prelim[i] = n_initfinal[i];

            assert(n_initfinal[i]);
            
        }
        return;
    }
    
    
    lft_char = left_nd->nd_initprelim;
    rt_char = right_nd->nd_initprelim;
    
    for (i = 0; i < num_chars; ++i) {
        
        if (n_initprelim[i] & MORPHY_INAPPLICABLE_BITPOS) {
            
            if (n_initprelim[i] & MORPHY_IS_APPLICABLE) {
                
                if (anc_char[i] == MORPHY_INAPPLICABLE_BITPOS) {
                    n_initfinal[i] = MORPHY_INAPPLICABLE_BITPOS;
                }
                else {
                    n_initfinal[i] = n_initprelim[i] & MORPHY_IS_APPLICABLE;
                    
                    assert(!(anc_char[i] & MORPHY_INAPPLICABLE_BITPOS));
                }
            }
            else {
                if (anc_char[i] == MORPHY_INAPPLICABLE_BITPOS) {
                    n_initfinal[i] = MORPHY_INAPPLICABLE_BITPOS;
                }
                else {
                    if ((lft_char[i] | rt_char[i]) & MORPHY_IS_APPLICABLE) {
                        n_initfinal[i] = (lft_char[i] | rt_char[i]) & MORPHY_IS_APPLICABLE;
                    }
                    else {
                        n_initfinal[i] = MORPHY_INAPPLICABLE_BITPOS;
                    }
                }
            }
        }
        else {
            n_initfinal[i] = n_initprelim[i];
        }
        
//        nd_active[i] |= ractive[i];
        assert(n_initfinal[i]);
    }
    
    
}


void mfl_second_fitch_na_downpass(mfl_nodedata_t*       n_nd,
                                  mfl_nodedata_t*       left_nd,
                                  mfl_nodedata_t*       right_nd,
                                  mfl_nodedata_t*       anc_nd,
                                  mfl_datapartition_t*  datapart,
                                  int*                  length)
{
    int i = 0;
    int num_chars = datapart->part_n_chars_included;
    int* weights = datapart->part_int_weights;
    mfl_charstate* lft_char = left_nd->nd_prelim_set;
    mfl_charstate* rt_char = right_nd->nd_prelim_set;
    mfl_charstate* n_init = n_nd->nd_initfinal;
    mfl_charstate* n_prelim = n_nd->nd_prelim_set;
//    mfl_charstate* actives = datapart->part_activestates;
//    mfl_charstate* ractives = datapart->part_tempactives;
    mfl_charstate* lft_active = left_nd->nd_subtree_activestates;
    mfl_charstate* rt_active = right_nd->nd_subtree_activestates;
    mfl_charstate* subtreeactive = n_nd->nd_subtree_activestates;
//    mfl_charstate* lreg_active = left_nd->nd_region_activestates;
//    mfl_charstate* rreg_active = right_nd->nd_region_activestates;
//    mfl_charstate* nreg_active = n_nd->nd_region_activestates;
    mfl_charstate temp = 0;
    
    
    
    for (i = 0; i < num_chars; ++i) {
        
        temp = 0;
        
//        nreg_active[i] = (lreg_active[i] | rreg_active[i]);
        
        if (n_init[i] & MORPHY_IS_APPLICABLE) {
            
            if ((temp = (lft_char[i] & rt_char[i]))) {
                
//                if (temp & MORPHY_IS_APPLICABLE) {
                    n_prelim[i] = temp & MORPHY_IS_APPLICABLE;
//                }
//                else {
//                    n_prelim[i] = n_init[i]; // & IS_APPLIC?
//                    assert(!(n_prelim[i] & MORPHY_INAPPLICABLE_BITPOS));
//                }
                
            } else {
                n_prelim[i] = (lft_char[i] | rt_char[i]) & MORPHY_IS_APPLICABLE;
                
                // Add steps if necessary:
                if (lft_char[i] & MORPHY_IS_APPLICABLE && rt_char[i] & MORPHY_IS_APPLICABLE) {
                    if (length) {
                        *length += weights[i];
                    }
                }
                else if (lft_active[i] && rt_active[i]) {
                    if (length) {
                        *length += weights[i];
                    }
                }
            }
        }
        else {
            n_prelim[i] = n_init[i];
        }
        
        subtreeactive[i] = (lft_active[i] | rt_active[i]) & MORPHY_IS_APPLICABLE;
        assert(n_prelim[i]);

    }
}

void mfl_set_subtree_actives(mfl_nodedata_t*       n_nd,
                             mfl_nodedata_t*       left_nd,
                             mfl_nodedata_t*       right_nd,
                             mfl_nodedata_t*       anc_nd,
                             mfl_datapartition_t*  datapart)
{
    int num_chars = datapart->part_n_chars_included;
    mfl_charstate* lft_active = left_nd->nd_subtree_activestates;
    mfl_charstate* rt_active = right_nd->nd_subtree_activestates;
    mfl_charstate* subtreeactive = n_nd->nd_subtree_activestates;
    
    for (int i = 0; i < num_chars; ++i) {
        subtreeactive[i] = (lft_active[i] | rt_active[i]) & MORPHY_IS_APPLICABLE;
//        assert(subtreeactive[i] < MORPHY_MISSING_DATA_BITWISE - 1);
    }
}


void mfl_second_fitch_na_uppass(mfl_nodedata_t*       n_nd,
                                mfl_nodedata_t*       left_nd,
                                mfl_nodedata_t*       right_nd,
                                mfl_nodedata_t*       anc_nd,
                                mfl_datapartition_t*  datapart,
                                int*                  length)
{
    int i = 0;
    int num_chars = datapart->part_n_chars_included;
    int* weights = datapart->part_int_weights;
    mfl_charstate* n_prelim = n_nd->nd_prelim_set;
    mfl_charstate* n_final = n_nd->nd_final_set;
    mfl_charstate* anc_char = anc_nd->nd_final_set;
//    mfl_charstate* actives = datapart->part_activestates;
    mfl_charstate* subtreeactive = n_nd->nd_subtree_activestates;
//    mfl_charstate* regionactive = n_nd->nd_region_activestates;
    mfl_charstate temp = 0;
    
    if (!left_nd) {
        
        assert(!right_nd);
        
        for (i = 0; i < num_chars; ++i) {
            
            if (n_prelim[i] & anc_char[i]) {
                n_final[i] = n_prelim[i] & anc_char[i];
                
                if (n_final[i] & MORPHY_IS_APPLICABLE) {
                    n_final[i] &= MORPHY_IS_APPLICABLE;
                }
            }
            else {
                n_final[i] = n_prelim[i];
            }

            subtreeactive[i] = (n_final[i] & MORPHY_IS_APPLICABLE);
            
            assert(n_final[i]);
        }
        
        return;
    }
    
    mfl_charstate* lft_char = left_nd->nd_prelim_set; // Note left and right are pointing to final sets now, not prelim sets
    mfl_charstate* rt_char = right_nd->nd_prelim_set;
    mfl_charstate* lft_active = left_nd->nd_subtree_activestates;
    mfl_charstate* rt_active = right_nd->nd_subtree_activestates;

    
    for (i = 0; i < num_chars; ++i) {
        
        temp = 0;
        
        if (n_prelim[i] & MORPHY_IS_APPLICABLE) {
            
             if (anc_char[i] & MORPHY_IS_APPLICABLE) {
                
                if ((anc_char[i] & n_prelim[i]) == anc_char[i]) {
                    
                    n_final[i] = anc_char[i] & n_prelim[i];
                    
                }
                else {
                    if (lft_char[i] & rt_char[i]) { // TODO: Use temp for micro-optimisation
                        
                        n_final[i] = (n_prelim[i] | (anc_char[i] & (lft_char[i] & rt_char[i])));
    
                    }
                    else {
                        
                        
                        
                    }
                }
             }
             else {
                 n_final[i] = n_prelim[i];
             }
        }
        else {
            
            n_final[i] = n_prelim[i];
            assert(n_final[i] == MORPHY_INAPPLICABLE_BITPOS);
            
            if (lft_active[i] && rt_active[i]) {
                if (length) {
                    *length += weights[i];
                }
            }
        }
        
        assert(n_final[i]);
    }

}


void mfl_set_rootstates(mfl_node_t* dummyroot,
                        mfl_node_t* rootnode,
                        mfl_partition_set_t* dataparts)
{
    
    int i = 0;
    int j = 0;
    int n_char_in_parts = 0;
    
    for (i = 0; i < dataparts->ptset_n_parts; ++i) {
        
        n_char_in_parts = dataparts->ptset_partitions[i]->part_n_chars_included;
        
        if (!dataparts->ptset_partitions[i]->part_has_inapplicables) {
            for (j = 0; j < n_char_in_parts; ++j) {
                
                dummyroot->nodet_charstates[i]->nd_final_set[j] = rootnode->nodet_charstates[i]->nd_prelim_set[j];
                if (dummyroot->nodet_charstates[i]->nd_final_set[j] & MORPHY_IS_APPLICABLE) {
                    dummyroot->nodet_charstates[i]->nd_final_set[j] = dummyroot->nodet_charstates[i]->nd_final_set[j] & MORPHY_IS_APPLICABLE;
                }
                dummyroot->nodet_charstates[i]->nd_subtree_activestates[j] = dummyroot->nodet_charstates[i]->nd_prelim_set[j];
                
            }
        }
        else {
            for (j = 0; j < n_char_in_parts; ++j) {
                
                dummyroot->nodet_charstates[i]->nd_initfinal[j] = rootnode->nodet_charstates[i]->nd_initprelim[j];
                
                if (dummyroot->nodet_charstates[i]->nd_initfinal[j] & MORPHY_IS_APPLICABLE) {
                    dummyroot->nodet_charstates[i]->nd_initfinal[j] = dummyroot->nodet_charstates[i]->nd_initfinal[j] & MORPHY_IS_APPLICABLE;
                    dummyroot->nodet_charstates[i]->nd_final_set[j] = dummyroot->nodet_charstates[i]->nd_initfinal[j];
                }
                
                dummyroot->nodet_charstates[i]->nd_subtree_activestates[j] = dummyroot->nodet_charstates[i]->nd_prelim_set[j];
            }
        }
        
    }
    
}

//void mfl_set_inapplic_rootstates(mfl_node_t* dummyroot,
//                                 mfl_node_t* rootnode,
//                                 mfl_partition_set_t* dataparts)
//{
//    
//    int i = 0;
//    int j = 0;
//    int n_char_in_parts = 0;
//    
//    for (i = 0; i < dataparts->ptset_n_parts; ++i) {
//        
//        n_char_in_parts = dataparts->ptset_partitions[i]->part_n_chars_included;
//        
//        for (j = 0; j < n_char_in_parts; ++j) {
//            dummyroot->nodet_charstates[i]->nd_final_set[j] = rootnode->nodet_charstates[i]->nd_final_set[j];
//            if (dummyroot->nodet_charstates[i]->nd_final_set[j] & MORPHY_IS_APPLICABLE) {
//                dummyroot->nodet_charstates[i]->nd_final_set[j] = dummyroot->nodet_charstates[i]->nd_final_set[j] & MORPHY_IS_APPLICABLE;
//            }
//        }
//    }
//    
//}

inline int mfl_wagner_stepcount(mfl_charstate leftchar,
                                mfl_charstate rightchar,
                                mfl_charstate* parentchar,
                                int weight)
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

//void mfl_wagner_downpass_binary_node(mfl_node_t *node)
//{
//    int i = 0;
//    mfl_node_t* lchild = NULL;
//    mfl_node_t* rchild = NULL;
//    lchild = node->nodet_next->nodet_edge;
//    rchild = node->nodet_next->nodet_next->nodet_edge;
//    mfl_charstate temp = NULL;
//    mfl_charstate* parentchars;// = node->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
//    mfl_charstate* leftchars;//   = lchild->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
//    mfl_charstate* rightchars;//  = rchild->nodet_dataparts[MFL_OPT_FITCH]->nd_prelim_set;
//    int num_chars;// = node->nodet_dataparts[MFL_OPT_FITCH]->nd_n_characters;
//    
//    for (i = 0; i < num_chars; ++i) {
//        if ((temp = leftchars[i] & rightchars[i])) {
//            
//        }
//        else {
//            parentchars[i] = leftchars[i] | rightchars[i];
//            /*Lenght increase = */ mfl_wagner_stepcount(leftchars[i], rightchars[i], &parentchars[i],  NULL/* WEIGHT goes here*/);
//        }
//    }
//}


void mfl_postorder_traversal(mfl_node_t *n, int* length)
{
    
    if (!n->nodet_downpass_visited) {
    
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

        n->nodet_downpass_visited = true;
        
        if (length) {
            n->nodet_isbottom = true;
        }
    }

    return;
}

int mfl_update_first_fitch_na_uppass
(mfl_nodedata_t* n_nd, mfl_nodedata_t* left_nd, mfl_nodedata_t* right_nd, mfl_nodedata_t* anc_nd, mfl_datapartition_t* datapart, int* length)
{
    
    int i = 0;
    int j = 0;
    int num_chars = datapart->nchanges;
    int *indices = datapart->part_char_changing;
    int updates = 0;

    mfl_charstate* lft_char = NULL;
    mfl_charstate* rt_char = NULL;
    mfl_charstate* n_initprelim = n_nd->nd_subtree_initprelim;
    mfl_charstate* n_initfinal = n_nd->nd_subtree_initfinal;
    mfl_charstate* n_prelim = n_nd->nd_subtree_prelim_set;
    mfl_charstate* anc_char = anc_nd->nd_subtree_initfinal;
    mfl_charstate* subtreeactive = n_nd->nd_subtree_activestates;
    mfl_charstate* n_originitfinal = n_nd->nd_initfinal;
    
    if (!left_nd) {
        assert(!right_nd);
        for (j = 0; j < num_chars; ++j) {
            
            i = indices[j];
            if (n_initprelim[i] & anc_char[i]) {
                subtreeactive[i] = (n_initprelim[i] & anc_char[i]) & MORPHY_IS_APPLICABLE;
            }
            else {
                subtreeactive[i] |= n_initprelim[i] & MORPHY_IS_APPLICABLE;
            }
            
            n_initfinal[i] = n_initprelim[i];
            
            if (n_initfinal[i] & anc_char[i]) {
                if (anc_char[i] & MORPHY_IS_APPLICABLE) {
                    n_initfinal[i] &= MORPHY_IS_APPLICABLE;
                }
            }
            
            n_prelim[i] = n_initfinal[i];
            
            assert(n_initfinal[i]);
        }
        return 0;
    }
    
    // TODO: memcpy the original statesets into these fiels, otherwise they'll be empty
    lft_char = left_nd->nd_subtree_initprelim;
    rt_char = right_nd->nd_subtree_initprelim;;
    
    for (j = 0; j < num_chars; ++j) {
        
        i = indices[j];
        
        if (n_initprelim[i] & MORPHY_INAPPLICABLE_BITPOS) {
            
            if (n_initprelim[i] & MORPHY_IS_APPLICABLE) {
                
                if (anc_char[i] == MORPHY_INAPPLICABLE_BITPOS) {
                    n_initfinal[i] = MORPHY_INAPPLICABLE_BITPOS;
                }
                else {
                    n_initfinal[i] = n_initprelim[i] & MORPHY_IS_APPLICABLE;
                    
                    assert(!(anc_char[i] & MORPHY_INAPPLICABLE_BITPOS));
                }
            }
            else {
                if (anc_char[i] == MORPHY_INAPPLICABLE_BITPOS) {
                    n_initfinal[i] = MORPHY_INAPPLICABLE_BITPOS;
                }
                else {
                    if ((lft_char[i] | rt_char[i]) & MORPHY_IS_APPLICABLE) {
                        n_initfinal[i] = (lft_char[i] | rt_char[i]) & MORPHY_IS_APPLICABLE;
                    }
                    else {
                        n_initfinal[i] = MORPHY_INAPPLICABLE_BITPOS;
                    }
                }
            }
        }
        else {
            n_initfinal[i] = n_initprelim[i];
        }
        
        if (n_initfinal[j] != n_originitfinal[j]) {
            ++updates;
        }
        
        //        nd_active[i] |= ractive[i];
        assert(n_initfinal[i]);
    }
    
    return updates;
}


int mfl_update_first_fitch_na_downpass
(mfl_nodedata_t* n_nd, mfl_nodedata_t* left_nd, mfl_nodedata_t* right_nd,mfl_nodedata_t* dummy, mfl_datapartition_t*  datapart, int* length)
{
    int i = 0;
    int j = 0;
    int updates = 0;
    int num_chars = datapart->nchanges;
    int *indices = datapart->part_char_changing;
    
    mfl_charstate* left = left_nd->nd_subtree_initprelim;   // TODO: Change these three to subtree records?
    mfl_charstate* right = right_nd->nd_subtree_initprelim;
    mfl_charstate* n_initprelim = n_nd->nd_subtree_initprelim;
    mfl_charstate* n_originitprelim = n_nd->nd_initprelim;
    
    mfl_charstate* lft_active = left_nd->nd_subtree_activestates; // TODO: What about these?
    mfl_charstate* rt_active = right_nd->nd_subtree_activestates;
    mfl_charstate* subtreeactive = n_nd->nd_subtree_activestates;
    mfl_charstate temp;
    
    for (i = 0; i < num_chars; ++i) {
        
        temp = 0;
        j = indices[i];
        
        if ((temp = (left[j] & right[j])) ) {
            
            n_initprelim[j] = temp;
            
            if (n_initprelim[j] == MORPHY_INAPPLICABLE_BITPOS) {
                
                if ((left[j] & MORPHY_IS_APPLICABLE) && (right[j] & MORPHY_IS_APPLICABLE)) {
                    n_initprelim[j] = left[j] | right[j];
                }
            }
        }
        else {
            
            n_initprelim[j] = left[j] | right[j];
            
            if ((left[j] & MORPHY_IS_APPLICABLE) && (right[j] & MORPHY_IS_APPLICABLE)) {
                n_initprelim[j] = n_initprelim[j] & MORPHY_IS_APPLICABLE;
            }
            
        }
        
        subtreeactive[j] = (lft_active[j] | rt_active[j]) & MORPHY_IS_APPLICABLE;
        
        if (n_initprelim[j] != n_originitprelim[j]) {
            ++updates;
        }
        
        assert(n_initprelim[j]);
    }
    
    return updates; // Return the number of characters updated at this node
}


void mfl_update_preorder(mfl_node_t *n, int *length)
{
    int i = 0;
    int num_dataparts;
    mfl_node_t *p = NULL;
    mfl_parsim_fn evaluator;
    mfl_node_t* left;
    mfl_node_t* right;
    mfl_nodedata_t* leftchars;
    mfl_nodedata_t* rightchars;
    
    num_dataparts = n->nodet_num_dat_partitions;
    
    if (!n->nodet_tip) {
        left = n->nodet_next->nodet_edge;
        right = n->nodet_next->nodet_next->nodet_edge;
    }
    
    
    
    for (i = 0; i < num_dataparts; ++i) {
        if (n->nodet_charstates[i]->nd_parent_partition->part_has_inapplicables) {
            
        }
    }
    
    if (n->nodet_tip) {
        return;
    }
    
    p = n->nodet_next;
    do {
        mfl_first_preorder_traversal(p->nodet_edge, length);
        p = p->nodet_next;
    } while (p != n);
    
    
    for (i = 0; i < num_dataparts; ++i) {
        if (n->nodet_charstates[i]->nd_parent_partition->part_has_inapplicables) {

        }
    }
    
    return;
}

void mfl_update_postorder(mfl_node_t *n, int *length)
{
    
    int i = 0;
    int num_dataparts;
    mfl_node_t *p = NULL;
//    mfl_parsim_fn evaluator;

    dbg_printf("Plumbing...\n");
    
    if (n->nodet_tip) {
        return;
    }
    
    mfl_node_t* left;
    mfl_node_t* right;
    
    p = n;
    if (!p->nodet_isbottom) {
        do {
            p = p->nodet_next;
            assert(p != n);
        } while (!p->nodet_isbottom);
    }
    
    left    = p->nodet_next->nodet_edge;
    right   = p->nodet_next->nodet_next->nodet_edge;
    
    num_dataparts = p->nodet_num_dat_partitions;
    mfl_datapartition_t *part = NULL;
    
    for (i = 0; i < num_dataparts; ++i) {
        part = p->nodet_charstates[i]->nd_parent_partition;
        if (part->nchanges) {

            // Optimise downpass for all character affected by insertion
            // Return the node for the last node affected by the downpass
            // Elsewhere:
            //      Begin first uppass on that node;
            //      Go up 1 node past the last node affected by the insertion?
            //      Pass down until state sets no longer affected
            //      Pass up until 1 past last node affected
            // Where to add/subtract length?
            
            int updates = mfl_update_first_fitch_na_downpass
                          (p->nodet_charstates[i], left->nodet_charstates[i], right->nodet_charstates[i], NULL, p->nodet_charstates[i]->nd_parent_partition, length);
            if (!updates) {
                dbg_printf("Finished plumbing, returning...\n");
                // Do uppass
                return;
            }
        }
    }
    
    if (p->nodet_edge->nodet_tip > 0) {
        mfl_update_postorder(p->nodet_edge, NULL);
    }
    else {
        // Do uppass:
    }
    
    return;
}



void mfl_first_preorder_traversal(mfl_node_t *n, int* length)
{
    int i = 0;
    int num_dataparts;
    mfl_node_t *p = NULL;
    mfl_parsim_fn evaluator;
    mfl_node_t* left;
    mfl_node_t* right;
    mfl_nodedata_t* leftchars;
    mfl_nodedata_t* rightchars;
    
    num_dataparts = n->nodet_num_dat_partitions;
    
    if (!n->nodet_tip) {
        left = n->nodet_next->nodet_edge;
        right = n->nodet_next->nodet_next->nodet_edge;
    }
    
    // For each data partition at the node, set the correct type and evaluation
    
    for (i = 0; i < num_dataparts; ++i) {
        
        evaluator = n->nodet_charstates[i]->nd_uppass_full;
        
        if (n->nodet_tip) {
            leftchars = NULL;
            rightchars = NULL;
        }
        else {
            leftchars = left->nodet_charstates[i];
            rightchars = right->nodet_charstates[i];
        }
        evaluator(n->nodet_charstates[i],
                  leftchars,
                  rightchars,
                  n->nodet_edge->nodet_charstates[i],
                  n->nodet_charstates[i]->nd_parent_partition,
                  length);
    }
    
    if (n->nodet_tip) {
        return;
    }
    
    p = n->nodet_next;
    do {
        mfl_first_preorder_traversal(p->nodet_edge, length);
        p = p->nodet_next;
    } while (p != n);
    
    
    for (i = 0; i < num_dataparts; ++i) {
        if (n->nodet_charstates[i]->nd_parent_partition->part_has_inapplicables) {
            
            evaluator = n->nodet_charstates[i]->nd_NAdownpass_full;

            evaluator(n->nodet_charstates[i],
                      left->nodet_charstates[i],
                      right->nodet_charstates[i],
                      n->nodet_edge->nodet_charstates[i],
                      n->nodet_charstates[i]->nd_parent_partition,
                      length);
        }
    }
    
    return;
}

void mfl_second_preorder_traversal(mfl_node_t *n, int* length)
{
    int i = 0;
    int num_dataparts;
    mfl_node_t *p = NULL;
    mfl_parsim_fn evaluator;
    mfl_node_t* left = NULL;
    mfl_node_t* right = NULL;
    mfl_nodedata_t* leftchars;
    mfl_nodedata_t* rightchars;
    
    num_dataparts = n->nodet_num_dat_partitions;
    
    if (!n->nodet_tip) {
        left = n->nodet_next->nodet_edge;
        right = n->nodet_next->nodet_next->nodet_edge;
    }
    
    // For each data partition at the node, set the correct type and evaluation
    
    for (i = 0; i < num_dataparts; ++i) {

        evaluator = n->nodet_charstates[i]->nd_NAuppass_full;
        if (evaluator) {
            if (n->nodet_tip) {
                leftchars = NULL;
                rightchars = NULL;
            }
            else {
                leftchars = left->nodet_charstates[i];
                rightchars = right->nodet_charstates[i];
            }
            
            evaluator(n->nodet_charstates[i],
                      leftchars,
                      rightchars,
                      n->nodet_edge->nodet_charstates[i],
                      n->nodet_charstates[i]->nd_parent_partition,
                      length);
        }
        
    }

    if (n->nodet_tip) {
        return;
    }
    
    p = n->nodet_next;
    
    do {
        mfl_second_preorder_traversal(p->nodet_edge, length);
        p->nodet_uppass_visited = false;
        p = p->nodet_next;
    } while (p != n);
    
    for (i = 0; i < num_dataparts; ++i) {
        if (n->nodet_charstates[i]->nd_parent_partition->part_has_inapplicables) {
            mfl_set_subtree_actives(p->nodet_charstates[i],
                                    p->nodet_next->nodet_edge->nodet_charstates[i],
                                    p->nodet_next->nodet_next->nodet_edge->nodet_charstates[i],
                                    p->nodet_edge->nodet_charstates[i],
                                    n->nodet_charstates[i]->nd_parent_partition);
        }
    }
    
    n->nodet_downpass_visited = false;
    
    return;
}



void mfl_allviews_traversal(mfl_node_t* n)
{
    mfl_node_t *p = NULL;
    
    if (n->nodet_tip) {
//        dbg_printf("Tip visited: %i\n", n->nodet_tip);
        mfl_postorder_traversal(n->nodet_edge, NULL);
        return;
    }
    
    p = n->nodet_next;
    do {
        mfl_allviews_traversal(p->nodet_edge);
        p = p->nodet_next;
    } while (p != n);
    
}


bool mfl_simple_unroot(mfl_tree_t *t, mfl_cliprec_t* clip)
{
    bool ret = false;
    
    if (!t->treet_root) {
        return true;
    }
    
    mfl_node_t *start_opp = mfl_find_rightmost_tip_in_tree(t->treet_root);
    
    
    
    if (mfl_clip_branch(t->treet_root->nodet_edge, clip)) {
        ret = true;
    }
    
    t->treet_start = start_opp->nodet_edge;
    
    return ret;
}

bool mfl_simple_reroot(mfl_tree_t* t, mfl_cliprec_t* clip)
{
    bool ret = false;
    
    mfl_restore_branching(clip);
    
    return ret;
}

bool mfl_calculate_all_views(mfl_tree_t* t, mfl_partition_set_t* dataparts, int *length)
{
    int num_taxa = t->treet_num_taxa;
    mfl_cliprec_t orig_root;
    mfl_node_t *entry = NULL;
    
    if (t->treet_root) {
        entry = t->treet_root;
    }
    else {
        entry = t->treet_start;
    }
    // Perform the first traversal
    mfl_postorder_traversal(entry, length);
    dbg_printf("Length after first traversal: %i\n", *length);
    
    // Unroot the tree
    if (t->treet_root) {
        if (!t->treet_root->nodet_weight) {
            t->treet_root->nodet_weight = num_taxa;
        }

    }
    // TODO: This really needs to be generalised.
    // Perform a simple unrooting by
    if (mfl_simple_unroot(t, &orig_root)) {
        
        mfl_allviews_traversal(t->treet_start);

        mfl_simple_reroot(t, &orig_root);

#ifdef MFY_DEBUG
//        tui_check_broken_tree(t, false);
#endif
        return true;
    }
    
#ifdef MFY_DEBUG
//    tui_check_broken_tree(t, false);
#endif
    
    return false;
}

void mfl_clear_active_states(mfl_partition_set_t* dataparts)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < dataparts->ptset_n_parts; ++i) {
        for (j = 0; j < dataparts->ptset_partitions[i]->part_n_chars_max; ++j) {
            if (dataparts->ptset_partitions[i]->part_has_inapplicables) {
                memset(dataparts->ptset_partitions[i]->part_activestates, 0, dataparts->ptset_partitions[i]->part_n_chars_max* sizeof(mfl_charstate));
            }
        }
    }
}

void mfl_fullpass_tree_optimisation(mfl_tree_t* t, mfl_partition_set_t* dataparts)
{
    
    // TODO: Needs rooted/unrooted handling
    t->treet_parsimonylength = 0;
    
    for (int i = 0; i < t->treet_num_nodes; ++i) {
        t->treet_treenodes[i]->nodet_isbottom = false;
    }
    
    // Perform the downpass in all directions
    mfl_calculate_all_views(t, dataparts, &t->treet_parsimonylength);
    
    // Perform uppass and first inapplicable downpass
    mfl_set_rootstates(&t->treet_dummynode, t->treet_root, dataparts);
//    mfl_clear_active_states(dataparts);
    mfl_first_preorder_traversal(t->treet_root, &t->treet_parsimonylength);

//    mfl_calculate_all_views2(t, dataparts, NULL);
    // Perform the final uppass
    mfl_set_rootstates(&t->treet_dummynode, t->treet_root, dataparts);
//    mfl_clear_active_states(dataparts);
    mfl_second_preorder_traversal(t->treet_root, &t->treet_parsimonylength);
    
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
    
    return d;
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
    
    return d;
}


// MDB: No it doesn't. This function won't work as expected.
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
