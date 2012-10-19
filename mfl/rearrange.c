/*
 *  rearrange.c
 *  Morphy
 *
 *  Functions used in tree rearrangements for heuristic searches. 
 *
 */

#include "morphy.h"
#include <time.h>

/*void check_double_skip(nodearray nds, int numnodes)
{
    int i, skipcount = 0;
    node *p;
    bool fndskips = false;
    
    for (i = 0; i < numnodes; ++i) {
        if (nds[i]->skip) {
            ++skipcount;
        }
        if (nds[i]->next) {
            p = nds[i]->next;
            while (p != nds[i]) {
                if (p->skip) {
                    ++skipcount;
                }
                p = p->next;
            }
        }
    }
    
    if (skipcount > 1) {
        dbg_printf("doubleskip: %i\n", skipcount);
    }
}*/

mfl_searchrec * mfl_create_searchrec(void)
{
    mfl_searchrec *newsearchrec;
    
    newsearchrec = (mfl_searchrec*)malloc(sizeof(mfl_searchrec));
    
    newsearchrec->nextinbuffer = 0;
    newsearchrec->undertreelimit = true;
    newsearchrec->bestlength = 0;
    newsearchrec->bestinrep = 0;
    newsearchrec->trbufstart = 0;
    newsearchrec->foundbettertr = false;
    newsearchrec->success = false;
    newsearchrec->niter_total = 0;
    newsearchrec->niter_ontree = 0;
    newsearchrec->currentreplicate = 0;
    
    return newsearchrec;
}

void mfl_part_reset_searchrec(mfl_searchrec *searchrec)
{
    /*Partially resets the searchrec. Used when a new replicate is started. */
    
    searchrec->foundbettertr = false;
    searchrec->success = false;
    searchrec->niter_ontree = 0;
}

void mfl_destroy_searchrec(mfl_searchrec *searchrec)
{
    free(searchrec);
}

long long int mfl_rearr_num(bool reset)
{
    static long long int rn = 0;

    if (reset) {
        rn = 0;
    }
    
    return ++rn;
}


void mfl_bswap(node *p, node *q)
{
    /* Exchanges the position of two nodes */
    
    node *pbase, *qbase;
    
    pbase = p->edge;
    qbase = q->edge;
    
    p->edge = qbase;
    q->edge = pbase;
    pbase->edge = q;
    qbase->edge = p;
    
}

void mfl_remove_branch(node *n)
{
    node *p, *q, *nb;
    
    nb = n->edge;
    p = nb->next->edge;
    q = nb->next->next->edge;
    nb->next->edge = NULL;
    nb->next->next->edge = NULL;
    mfl_join_nodes(p, q);
}

void mfl_insert_branch(node *br, node *target, int ntax)
{
    // Inserts a branch with a ring base into another branch
    
    node *br1, *br2, *tdesc;
    
    tdesc = target->edge;
    
    if (br->tip) {
        br1 = br->edge->next;
        br2 = br1->next;
    }
    else {
        br1 = mfl_seek_ringnode(br, ntax);
        br2 = mfl_seek_ringnode(br1->next, ntax);
    }

    if (br1->edge || br2->edge) {
        dbg_printf("Error in branch insertion\n");
        return;
    }
    
    mfl_join_nodes(br1, target);
    mfl_join_nodes(br2, tdesc);
}

/* Nearest-neighbor interchange (NNI) */

void mfl_nni_traversal(node *n, tree *swapingon, tree **savedtrees, int ntax, 
                       int nchar, int numnodes, long int *current, 
                       charstate *tipdata, bool *undertreelimit, 
                       int *currentbesttree, bool *foundbettertree)
{
    int trlength = 0;
    node *p;
    
    if (n->start) {
        mfl_nni_traversal(n->edge, swapingon, savedtrees, ntax, nchar, numnodes, 
                          current, tipdata, undertreelimit, currentbesttree,
                          foundbettertree);
        return;
    }
    
    if (n->tip || !(*undertreelimit)) {
        return;
    }
    
    if (!n->edge->tip)
    {
        /* Nearly identical steps are conducted twice because all nearest-
         * neighbor interchanges produce two distinct tree topologies */
        
        if (*current + 1 <= TREELIMIT) {
            mfl_bswap(n->next->edge, n->edge->next->edge);
            mfl_root_tree(swapingon, 1, ntax);
            trlength = mfl_get_treelen(swapingon, ntax, nchar, currentbesttree);
            mfl_unroot(ntax, swapingon);
            if (trlength < *currentbesttree) {
                *foundbettertree = true;
                *currentbesttree = trlength;
                swapingon->length = trlength;
                mfl_reinit_treebuffer(savedtrees, swapingon, current, numnodes);
                trlength = 0;
                *current = *current + 1;
                return;
            }
            if (trlength == *currentbesttree) {
                *foundbettertree = false;
                savedtrees[*current] = mfl_copytree(swapingon, ntax, numnodes);
                savedtrees[*current]->index = *current;
                savedtrees[*current]->length = trlength;
                trlength = 0;
                *current = *current + 1;
            }
            mfl_bswap(n->next->edge, n->edge->next->edge);
            
            if (*current + 1 <= TREELIMIT){ 
                mfl_bswap(n->next->next->edge, n->edge->next->edge);
                mfl_root_tree(swapingon, 1, ntax);
                trlength = mfl_get_treelen(swapingon, ntax, nchar, currentbesttree);
                mfl_unroot(ntax, swapingon);
                if (trlength < *currentbesttree) {
                    *foundbettertree = true;
                    *currentbesttree = trlength;
                    swapingon->length = trlength;
                    mfl_reinit_treebuffer(savedtrees, swapingon, current, numnodes);
                    trlength = 0;
                    *current = *current + 1;
                    return;
                }
                if (trlength == *currentbesttree) {
                    *foundbettertree = false;
                    savedtrees[*current] = mfl_copytree(swapingon, ntax, numnodes);
                    savedtrees[*current]->index = *current;
                    savedtrees[*current]->length = trlength;
                    trlength = 0;
                    *current = *current + 1;
                }
                mfl_bswap(n->next->next->edge, n->edge->next->edge);
            }
        }
        else {
            dbg_printf("Hit tree limit\n");
            *undertreelimit = false;
            return;
        }
    }
    
    p = n->next;
    while (p != n) {
        mfl_nni_traversal(p->edge, swapingon, savedtrees, ntax, nchar, numnodes, 
                          current, tipdata, undertreelimit, currentbesttree,
                          foundbettertree);
        p = p->next;
    }
}


/* Quite a lot about this function will change. For now, it interfaces with the
 * test in main.c but that won't always be the case. */
void mfl_nni_search(int ntax, int nchar, int numnodes, charstate *tipdata, 
                    tree **savedtrees, int starttreelen)
{
    long int i = 0, j = 0; 
    int currentbest = starttreelen;
    long int trbufpos = 0;
    long int *nxtintrbuf = &trbufpos; 
    int *currentbest_p = &currentbest;
    long int numreps = 100;
    bool undertreelimit = true, *undertreelimit_p = &undertreelimit;
    bool foundbettertree = false, *foundbettertree_p = &foundbettertree;
    
    do {
        for (j = 0; j == 0 || j < *nxtintrbuf; ++j) {
            mfl_nni_traversal(savedtrees[j]->trnodes[0], savedtrees[j], savedtrees, 
                              ntax, nchar, numnodes, nxtintrbuf, tipdata, 
                              undertreelimit_p, currentbest_p, foundbettertree_p);
        }
        ++i;
    } while (i < numreps);
    
    dbg_printf("Next in trbuf: %li\n", *nxtintrbuf);
    
    dbg_printf("\nThe optimal tree(s) found by nearest neighbor interchange:\n");
    for (i = 0; i < *nxtintrbuf; ++i) {
        dbg_printf("Tree %li:\n", i + 1);
        mfl_root_tree(savedtrees[i], 0, ntax);
        /*printNewick(savedtrees[i]->root);*/
        dbg_printf(";\n");
        dbg_printf("Length: %i\n", savedtrees[i]->length);
    }
    dbg_printf("\n\n");
    mfl_clear_treebuffer(savedtrees, nxtintrbuf, numnodes);
    
    free(savedtrees);
}

long int mfl_spr_leftotry(int ntax)
{
    long int rleft = 0;
    rleft = 4 * (ntax - 3) * (ntax - 2);
    return rleft;
}

void mfl_regrafting_traversal(node *n, node *subtr, tree *swapingon, 
                              tree **savedtrees, int ntax, int nchar, int numnodes, mfl_searchrec *searchrec, int diff)
{
    /* Called from within any subtree pruning algorithm used in either SPR or
     * TBR branch swapping. Traverses a binary tree in preorder, inserting the 
     * subtree at each node and (hopefully) skipping a reinsertion at the original
     * site of pruning (flagged by the "visited" boolean values in those nodes) */
    
    if (searchrec->success) {
        return;
    }
    
    static long int counter = 0;
    int trlength = 0;
    int al = 0;
    node *up;
    static bool island = false;
    
    if (!(n->visited) && !(n->edge->visited)) {
        up = n->edge;
        
        al = mfl_locreopt_cost(subtr->next->edge, n, up, nchar, diff);
        trlength = searchrec->bestinrep - diff + al;
        
        searchrec->niter_total = searchrec->niter_total + 1;

        if (trlength < searchrec->bestinrep) {
            
            island = true;
            
            counter = 0;
            
            mfl_join_nodes(subtr->next->next, up);
            mfl_join_nodes(subtr, n);
            
            searchrec->foundbettertr = true;
            searchrec->success = true;
            swapingon->length = trlength;
            searchrec->bestinrep = trlength;
            
            if (searchrec->bestinrep < searchrec->bestlength) {
                searchrec->bestlength = searchrec->bestinrep;
                mfl_reinit_treebuffer(savedtrees, swapingon, &searchrec->nextinbuffer, numnodes);
            }
            else {
                mfl_reinit_tbinrange(savedtrees, swapingon, searchrec->trbufstart, &searchrec->nextinbuffer, numnodes);
            }
            free(swapingon->compressedtr);
            swapingon->compressedtr = mfl_compress_tree(swapingon, ntax, numnodes);
            searchrec->nextinbuffer = searchrec->nextinbuffer + 1;
            
            return;
        }
        
        searchrec->foundbettertr = false;
        searchrec->success = false;
        
        if (trlength == searchrec->bestinrep) 
        {
            ++counter;
            
            node *down, *top;
            if (subtr->tocalcroot) {
                down = subtr;
                top = subtr->next->next;
            }
            else {
                down = subtr->next->next;
                top = subtr;
            }
            
            mfl_join_nodes(top, up);
            mfl_join_nodes(down, n);

            
            if (!mfl_compare_alltrees(swapingon, savedtrees, ntax, numnodes, &searchrec->trbufstart, searchrec->nextinbuffer)) 
            {
                if (island) {
                    dbg_printf("length: %i\n", searchrec->bestinrep);
                    island = false;
                }
                if (searchrec->currentreplicate > 0) {
                    if (trlength == searchrec->bestlength) {
                        long int bstart = 0;
                        if (mfl_compare_alltrees(swapingon, savedtrees, ntax, numnodes, &bstart, searchrec->trbufstart)) {
                            mfl_reinit_tbinrange(savedtrees, swapingon, searchrec->trbufstart, &searchrec->nextinbuffer, numnodes);
                            searchrec->success = true;
                            return;
                        }
                    }
                }
                
                savedtrees[searchrec->nextinbuffer] = mfl_copytree(swapingon, ntax, numnodes);
                savedtrees[searchrec->nextinbuffer]->index = searchrec->nextinbuffer;
                savedtrees[searchrec->nextinbuffer]->length = trlength;
                savedtrees[searchrec->nextinbuffer]->swapped = false;
                searchrec->nextinbuffer = searchrec->nextinbuffer + 1;
                
            }
            trlength = 0;
        }
        
        mfl_join_nodes(n, up);
    }
    
    if (n->tip) {
        return;
    }

    mfl_regrafting_traversal(n->next->edge, subtr, swapingon, savedtrees, ntax,
                                 nchar, numnodes, searchrec, diff);
    if (searchrec->success) {
        return;
    }
    mfl_regrafting_traversal(n->next->next->edge, subtr, swapingon, savedtrees, ntax,
                             nchar, numnodes, searchrec, diff);

}

void mfl_set_updown(node *n, node **up, node **dn)
{
    if (n->next->tocalcroot) {
        *dn = n->next->edge;
        *up = n->next->next->edge;
    }
    else {
        *up = n->next->edge;
        *dn = n->next->next->edge;
    }
}


void mfl_subtree_pruning(node *n, tree *swapingon, tree **savedtrees, int ntax, 
                         int nchar, int numnodes, mfl_searchrec *searchrec)
{
    
    /* Indexes over the node array of a tree, pruning subtrees and calling
     * the regrafting function */
    
    int i, diff, *srcchanging, *tgtchanging;
    node *p, *q, *src, *tgt, *s_dn, *s_up, *t_dn, *t_up, *t_dn_N, *t_up_N;
    node *s_dn_N, *s_up_N;
    nodearray nds = swapingon->trnodes;
    
    
    for (i = ntax + 1; i < numnodes; ++i) {
        
        p = nds[i];
        q = nds[i];
        
        p->skip = true;
        
        do {
            if (!p->skip) {
                src = p->edge;
                tgt = p;
                
                if (src->tip) {
                    memcpy(src->apomorphies, src->tempapos, 
                           nchar * sizeof(charstate));
                }
                else if (src->tocalcroot) {
                    mfl_set_rootstates(src, nchar);
                }
                else {
                    
                    mfl_set_updown(src, &s_up, &s_dn);
                    srcchanging = mfl_get_tgt_changing(s_dn->edge, s_up, s_dn, 
                                                       nchar);
                    
                    s_up_N = s_up->edge;
                    s_dn_N = s_dn->edge;
                    mfl_join_nodes(s_dn, s_up);
                    mfl_partial_downpass(s_dn, swapingon, numnodes, ntax, nchar, 
                                         srcchanging);
                    mfl_join_nodes(s_up_N, s_up);
                    mfl_join_nodes(s_dn_N, s_dn);
                    mfl_reopt_subtr_root(src, nchar);
                    free(srcchanging);
                }
                
                if (!tgt->tocalcroot) {
                    
                    mfl_set_updown(tgt, &t_up, &t_dn);
                    
                    tgtchanging = mfl_get_tgt_changing(t_dn->edge, t_up, t_dn, 
                                                       nchar);
                    
                    // Pop out redundant node:
                    t_up_N = t_up->edge;
                    t_dn_N = t_dn->edge;
                    mfl_join_nodes(t_up, t_dn);
                    
                    mfl_partial_downpass(t_dn, swapingon, numnodes, ntax, nchar, 
                                         tgtchanging);
                    free(tgtchanging);
                }
                else {
                    mfl_set_updown(tgt, &t_up, &t_dn);
                    // Pop out the redundant node:
                    t_up_N = t_up->edge;
                    t_dn_N = t_dn->edge;
                    mfl_join_nodes(t_up, t_dn);
                }
                
                diff = mfl_subtr_reinsertion(src, t_up, t_dn, nchar);
                
                // Perform all reinsertions of SRC on TGT
                t_up->visited = true;
                t_dn->visited = true;
                mfl_regrafting_traversal(t_up, src->edge->next->next, swapingon, 
                                         savedtrees, ntax, nchar, numnodes, 
                                         searchrec, diff);
                
                mfl_regrafting_traversal(t_dn, src->edge->next->next, swapingon, 
                                         savedtrees, ntax, nchar, numnodes, 
                                         searchrec, diff);
                
                t_up->visited = false;
                t_dn->visited = false;
                p->skip = false;
                if (searchrec->success) {
                    return;
                }
                
                mfl_join_nodes(t_up, t_up_N);
                mfl_join_nodes(t_dn, t_dn_N);
                mfl_restore_origstates(swapingon, ntax, numnodes, nchar);
            }
            else {
                p->skip = false;
            }

            p = p->next;
            
        } while (p != q);
        
    }
}

node *mfl_find_atip(node *n)
{
    node *tip = NULL;
    
    if (n->tip) {
        tip = n;
        return tip;
    }
 
    tip = mfl_find_atip(n->next->edge);
    
    return tip;
}

void mfl_reroot_subtree(node *n, node *atip, node *subtr, node *base, node *up, 
                        node *dn, tree *swapingon, tree **savedtrees, int ntax, 
                        int nchar, int numnodes, mfl_searchrec *searchrec, 
                        int diff)
{
    /* Traverses the subtree and re-roots it at each branch in preorder and then
     * calls regrafting traversal (same as used in SPR).*/
    
    if (searchrec->success) {
        return;
    }

    mfl_join_nodes(base->next->next, n->edge);
    mfl_join_nodes(base->next, n);
    
    // Reoptimize the subtree base

    mfl_reopt_subtr_root(base, nchar);
    
    if (!up->origbase && !dn->origbase) {
        up->visited = 0;
        dn->visited = 0;
    }
    else {
        base->skip = false;
    }

    
    searchrec->niter_total = searchrec->niter_total + 1;
    
    mfl_regrafting_traversal(up, subtr, swapingon, savedtrees, ntax, nchar, 
                             numnodes, searchrec, diff);
    if (searchrec->success) {
        up->visited = 1;
        dn->visited = 1;
        base->skip = true;
        return;
    }
    mfl_regrafting_traversal(dn, subtr, swapingon, savedtrees, ntax, nchar, 
                             numnodes, searchrec, diff);
    
    up->visited = 1;
    dn->visited = 1;
    base->skip = true;
    
    if (searchrec->success) {
        return;
    }
    
    // Remove the base
    mfl_join_nodes(base->next->edge, base->next->next->edge);
    
    if (n->tip) {
        return;
    }
    
    mfl_reroot_subtree(n->next->edge, atip, subtr, base, up, dn, swapingon, 
                       savedtrees, ntax, nchar, numnodes, searchrec, diff);
    if (searchrec->success) {
        return;
    }
    mfl_reroot_subtree(n->next->next->edge, atip, subtr, base, up, dn, 
                       swapingon, savedtrees, ntax, nchar, numnodes, searchrec, 
                       diff);
    
}

bool mfl_subtr_isrerootable(node *n)
{
    bool rerootable = true;
    
    if(n->tip) {
        return false;
    }
    
    else {
        if (n->next->edge->tip) {
            if (n->next->next->edge->tip) {
                rerootable = false;
            }
        }
    }
    
    return rerootable;
}

void mfl_bisection_traversal(node *n, tree *swapingon, tree **savedtrees, 
                             int ntax, int nchar, int numnodes, 
                             mfl_searchrec *searchrec)
{
    /* Traverses a binary tree clipping out a subtree in postorder and passing 
     * a pointer to the subtree to mfl_reroot_subtree. */
    
    
    int i, diff, *srcchanging, *tgtchanging;
    node *p, *q, *src, *tgt, *s_dn, *s_up, *t_dn, *t_up, *t_dn_N, *t_up_N;
    node *s_dn_N, *s_up_N, *atip;
    nodearray nds = swapingon->trnodes;
    
    for (i = ntax + 1; i < numnodes; ++i) {
        
        p = nds[i];
        q = nds[i];
        p->skip = true;
        
        do {
            if (!p->skip) {
                src = p->edge;
                tgt = p;
                
                if (src->tip) {
                    memcpy(src->apomorphies, src->tempapos, 
                           nchar * sizeof(charstate));
                }
                else if (src->tocalcroot) {
                    mfl_set_updown(src, &s_up, &s_dn);
                    s_up_N = s_up->edge;
                    s_dn_N = s_dn->edge;
                    
                    srcchanging = mfl_get_subtr_changing(src, NULL, NULL, nchar);
                    mfl_reopt_subtr(src, swapingon, nchar, numnodes, srcchanging);
                }
                else {
                    
                    mfl_set_updown(src, &s_up, &s_dn);
                    srcchanging = mfl_get_tgt_changing(s_dn->edge, s_up, s_dn, 
                                                       nchar);
                    
                    s_up_N = s_up->edge;
                    s_dn_N = s_dn->edge;
                    mfl_join_nodes(s_dn, s_up);
                    
                    mfl_partial_downpass(s_dn, swapingon, numnodes, ntax, nchar, 
                                         srcchanging);
                    
                    mfl_join_nodes(s_up_N, s_up);
                    mfl_join_nodes(s_dn_N, s_dn);
                    
                    //mfl_reopt_subtr_root(src, nchar);
                    free(srcchanging);
                }
                
                if (!tgt->tocalcroot) {
                    
                    mfl_set_updown(tgt, &t_up, &t_dn);
                    
                    tgtchanging = mfl_get_tgt_changing(t_dn->edge, t_up, t_dn, 
                                                       nchar);
                    
                    // Pop out redundant node:
                    t_up_N = t_up->edge;
                    t_dn_N = t_dn->edge;
                    mfl_join_nodes(t_up, t_dn);
                    
                    mfl_partial_downpass(t_dn, swapingon, numnodes, ntax, nchar, 
                                         tgtchanging);
                    free(tgtchanging);
                }
                else {
                    mfl_set_updown(tgt, &t_up, &t_dn);
                    // Pop out the redundant node:
                    t_up_N = t_up->edge;
                    t_dn_N = t_dn->edge;
                    mfl_join_nodes(t_up, t_dn);
                }
                
                diff = mfl_subtr_reinsertion(src, t_up, t_dn, nchar);
                
                // Perform all reinsertions of SRC on TGT
                t_up->visited = true;
                t_dn->visited = true;
                
                if (mfl_subtr_isrerootable(src)) {
                    
                    atip = mfl_find_atip(src);
                    
                    mfl_join_nodes(s_up, s_dn);
                    
                    mfl_reroot_subtree(atip->edge, atip, src->edge->next->next, 
                                       src, t_up, t_dn, swapingon, savedtrees, 
                                       ntax, nchar, numnodes, searchrec, diff);
                    
                    if (searchrec->success) {
                        t_up->visited = false;
                        t_dn->visited = false;
                        p->skip = false;
                        return;
                    }
                    mfl_join_nodes(s_up_N, s_up);
                    mfl_join_nodes(s_dn_N, s_dn);
                    
                }
                else {
                    mfl_regrafting_traversal(t_up, src->edge->next->next, swapingon, 
                                             savedtrees, ntax, nchar, numnodes, 
                                             searchrec, diff);
                    
                    mfl_regrafting_traversal(t_dn, src->edge->next->next, swapingon, 
                                             savedtrees, ntax, nchar, numnodes, 
                                             searchrec, diff);
                }
                
                t_up->visited = false;
                t_dn->visited = false;
                p->skip = false;
                
                if (searchrec->success) {
                    return;
                }
                
                mfl_join_nodes(t_up, t_up_N);
                mfl_join_nodes(t_dn, t_dn_N);
                
                mfl_restore_origstates(swapingon, ntax, numnodes, nchar);
            }
            else {
                p->skip = false;
            }

            p = p->next;
            
        } while (p != q);
        
    }
    
    /*int diff = 0;
    node *p, *up, *dn, *up1, *dn1, *subtr, *atip, *base, *bc1, *bc2;
        
    if (n->start) {
        mfl_bisection_traversal(n->edge, swapingon, savedtrees, ntax,
                              nchar, numnodes, searchrec);
        return;
    }
    
    if (n->tip || searchrec->success) {
        return;
    }
    
    dn = n->edge;
    
    p = n->next;
    up = p->next->edge;
    
    while (p != n) {
        
        mfl_bisection_traversal(p->edge, swapingon, savedtrees, ntax,
                              nchar, numnodes, searchrec);
        
        subtr = p->next->next;
        base = subtr->next->edge;
        
        if (!base->skip && searchrec->tipstovisit) {
            
            if (searchrec->success) {
                return;
            }*/
            
            /* If the clipped subtree has only 1 or 2 terminals, then it 
             * cannot be rerooted. Therefore, only a normal SPR routine is
             * required. Otherwise, the procedure can proceed with rerooting
             * and reconnecting. */
            
            /*if (mfl_subtr_isrerootable(subtr)) {
                
                up->visited = 1;
                dn->visited = 1;
                
                up1 = up->edge;
                dn1 = dn->edge;
                
                base->skip = true;
                bc1 = base->next->edge;
                bc2 = base->next->next->edge;
                bc1->origbase = true;
                bc2->origbase = true;

                mfl_changing *changing = mfl_get_changing(base, n, up, dn, nchar);
                
                //Clip the tree
                mfl_join_nodes(up, dn);
                up->clip = true;
                
                //mfl_trav_allviews(swapingon->trnodes[0], swapingon, ntax, nchar, changing->tgtchanging);
                
                mfl_partial_downpass(dn, swapingon, numnodes, ntax, nchar, changing->tgtchanging);
                
                atip = mfl_find_atip(base);
                
                mfl_reopt_subtr(base, swapingon, nchar, numnodes, changing->srcchanging);
                
                diff = mfl_subtr_reinsertion(base, up, dn, nchar);
                
                // Unroot the source tree
                mfl_join_nodes(bc1, bc2);
                
                // Call rerooting function
                
                mfl_reroot_subtree(atip->edge, atip, subtr, base, up, dn, swapingon, savedtrees, ntax, nchar, numnodes, searchrec, diff);                    
                
                up->visited = 0;
                dn->visited = 0;
                bc1->origbase = false;
                bc2->origbase = false;
                up->clip = false;
                dn->clip = false;
                base->skip = false;
                
                free(changing->srcchanging);
                free(changing->tgtchanging);
                free(changing);
                if (searchrec->success) {
                    return;
                }
                
                // Reroot the source tree on its base
                mfl_join_nodes(bc1, base->next);
                mfl_join_nodes(bc2, base->next->next);
                
                // Reinsert the subtree
                mfl_join_nodes(up, up1);
                mfl_join_nodes(dn, dn1); 
                
                mfl_restore_origstates(swapingon, ntax, numnodes, nchar);
            }
            else {
                mfl_spr_cliptrees(n, up, dn, subtr, swapingon, savedtrees, ntax, nchar, numnodes, searchrec);
                if (searchrec->success) {
                    return;
                }
            }
        }
        else {
            base->skip = false;
        }

                
        up = n->next->edge;
        p = p->next;
    }   */
    
    
}

void mfl_store_results(mfl_handle_s *mfl_handle, tree **savedtrees, int ntax)
{
    long int i = 0;
    mfl_resultant_data_s *a_results = mfl_handle->resultant_data;
    a_results->newicktrees.reserve((long unsigned int)a_results->n_savetrees);
    string *tstring;
    
    for (i = 0; i < a_results->n_savetrees; ++i) {
        tstring = mfl_trstring(savedtrees[i], ntax);
        a_results->newicktrees.push_back(*tstring);
    }
    
}

void (*mfl_swap_controller(mfl_handle_s *mfl_handle)) (node*, tree*, tree**, int, int , int, mfl_searchrec*)
{
    switch (mfl_handle->bswap_type) {
        case MFL_BST_TBR:
            dbg_printf("Searching by TBR\n");
            return &mfl_bisection_traversal;
            break;
        case MFL_BST_SPR:
            dbg_printf("Searching by SPR\n");
            return &mfl_subtree_pruning;
        case MFL_BST_NNI:
            dbg_printf("Not implemented\n");
            //return; // Temporary. Would fail if called.
            break;
        default:
            dbg_printf("Not implemented\n");
            break;
    }
    return NULL;
}

bool mfl_heuristic_search(mfl_handle_s *mfl_handle)
{
    int ntax = mfl_handle->n_taxa, nchar = mfl_handle->n_chars;
    int numnodes = mfl_calc_numnodes(ntax);
    long int i = 0, j = 0;
    long int nreps = mfl_handle->n_iterations;
    double timein = 0;
    double timeout = 0;
    bool quit = false;
    
    void (*branch_swapper)(node*, tree*, tree**, int, int, int, mfl_searchrec*) = NULL;
    
    tree **savedtrees = (tree**) malloc(TREELIMIT * sizeof(tree*));
    
    tree *newreptree;
    
    charstate *tipdata = mfl_convert_tipdata(mfl_handle->input_data, mfl_handle->n_taxa, mfl_handle->n_chars, mfl_handle->gap_as_missing);
    
    mfl_searchrec *searchrec = mfl_create_searchrec();
        
    timein = (double)(clock() / (double)CLOCKS_PER_SEC);
    
    /* This outer loop makes it possible to do multiple replicates of random
     * addition sequence. The variable nreps is the number of times an initial
     * tree is generated using random addition sequence. */
    
    branch_swapper = mfl_swap_controller(mfl_handle);
    
    for (i = 0; i < nreps; ++i) {
        
        newreptree = mfl_addseq_randasis(ntax, nchar, numnodes, tipdata, mfl_handle->addseq_type, savedtrees);
        searchrec->currentreplicate = i;
        
        if (i == 0) {
            dbg_printf("Replicate: %li\n", i + 1);
            /* The particular addition sequence will have to be selected by the user */
            savedtrees[0] = newreptree;
            savedtrees[0]->compressedtr = mfl_compress_tree(savedtrees[0], ntax, numnodes);
            searchrec->bestinrep = mfl_all_views(savedtrees[0], ntax, nchar, &searchrec->bestinrep);
            searchrec->bestlength = searchrec->bestinrep;
            j = 0;
            dbg_printf("The length of the starting tree: %i steps\n\n", searchrec->bestinrep);
        }
        else {
            dbg_printf("Replicate: %li\n", i + 1);
            dbg_printf("next in buff at start of rep: %li\n", searchrec->nextinbuffer);
            savedtrees[searchrec->nextinbuffer] = newreptree;
            searchrec->bestinrep = mfl_all_views(savedtrees[searchrec->nextinbuffer], ntax, nchar, &searchrec->bestinrep);
            //dbg_printf("Best length in replicate: %i\n", searchrec->bestinrep);
            j = searchrec->nextinbuffer;
            searchrec->trbufstart = searchrec->nextinbuffer;
            searchrec->foundbettertr = false;
            //dbg_printf("j = %li\n", searchrec->nextinbuffer);
            quit = false;
            //break;
        }
        
        do {
            //if (i > 0) {
                //dbg_printf("swapping on tree: %li\n", j);
            //}
            mfl_part_reset_searchrec(searchrec);
            //searchrec->foundbettertr = false;
            //searchrec->success = false;
            mfl_apply_tipdata(savedtrees[j], tipdata, ntax, nchar);
            mfl_all_views(savedtrees[j], ntax, nchar, &searchrec->bestinrep);
            mfl_devisit_tree(savedtrees[j]->trnodes, numnodes);

            /* The branch swapper is the specific type of heuristic search 
             * routine: either TBR, SPR, or NNI. TBR by default */
            mfl_save_origstates(savedtrees[j], ntax, numnodes, nchar);
            branch_swapper(savedtrees[j]->trnodes[0], savedtrees[j], savedtrees, 
                                  ntax, nchar, numnodes, searchrec);
            
            if (searchrec->foundbettertr) {
                if (i > 0) {
                    //dbg_printf("trbuf start %li\n", searchrec->trbufstart);
                }
                j = searchrec->trbufstart;
            }
            else {
                savedtrees[j]->swapped = true;
                if (savedtrees[j]->trnodes) {
                    mfl_free_trnodes(savedtrees[j], numnodes);
                }
                ++j;
            }
            
            if (j >= searchrec->nextinbuffer || !searchrec->undertreelimit) {
                dbg_printf("number of rearrangements tried: %li\n", searchrec->niter_total);
                quit = true;
            }
            
        } while (!quit);
        
        dbg_printf("Next in buffer at end of rep: %li\n", searchrec->nextinbuffer);
        
        //dbg_printf("best in rep: %i\n", searchrec->bestinrep);
        //dbg_printf("best overall: %i\n", searchrec->bestlength);
        
        if (i != 0) {
            if (searchrec->bestinrep > searchrec->bestlength) {
                mfl_reinit_tbinrange(savedtrees, savedtrees[0], searchrec->trbufstart, &searchrec->nextinbuffer, numnodes);
            }
        }
    }
    
    timeout = (double)(clock() / (double)CLOCKS_PER_SEC);
    
    mfl_handle->resultant_data = (mfl_resultant_data_s*) malloc(sizeof(mfl_resultant_data_s));
    
    mfl_handle->resultant_data->bestlength = searchrec->bestlength;
    mfl_handle->resultant_data->n_rearrangements = searchrec->niter_total;
    mfl_handle->resultant_data->n_savetrees = searchrec->nextinbuffer;
    mfl_handle->resultant_data->searcht = (timeout - timein);
    
    //mfl_store_results(mfl_handle, savedtrees, ntax);
    
    /* TESTING ONLY. This is just for checking output as I build up the heuristic
     * search procedure. Eventually, all this stuff will be written to a struct
     * and handed over to the interface for outputting to screen. */
    
    dbg_printf("Total search time: %g\n", timeout - timein);
    
    dbg_printf("Number of saved trees: %li\n", searchrec->nextinbuffer);
    
    /*dbg_printf("\nThe optimal tree(s) found by heuristic search:\n");
    for (j = 0; j < searchrec->nextinbuffer; ++j) {
        //dbg_printf("TREE str_%li = [&U] ", j+1);
        //mfl_root_tree(savedtrees[j], 0, ntax);
        //printNewick(savedtrees[j]->root);
        cout << "TREE str_" << j+1 << " = [&U] " << mfl_handle->resultant_data->newicktrees[i] << endl;
        //dbg_printf(";\n");
    }
    dbg_printf("\n");*/
    
    /* END OF TESTING-ONLY SECTION */
    
    mfl_clear_treebuffer(savedtrees, &searchrec->nextinbuffer, numnodes);
    free(savedtrees);
    free(tipdata);
    mfl_destroy_searchrec(searchrec);
    
    return true;
}
