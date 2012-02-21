/*
 *  coptim.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 2/15/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "morphy.h"


void mfl_apply_tipdata(tree *currenttree, charstate *tipdata, int ntax, int nchar)
{   
    int i;
    
    for (i = 0; i < ntax; ++i) {
        currenttree->trnodes[i]->apomorphies = &tipdata[i * nchar];
        currenttree->trnodes[i]->tempapos = &tipdata[i * nchar];
    }
}

int mfl_locreopt_cost(node *src, node *tgt1, node *tgt2, int nchar)
{
    int i;
    int cost = 0;
    
    for (i = 0; i < nchar; ++i) {
        if (!(src->tempapos[i] & (tgt1->tempapos[i] | tgt2->tempapos[i]))) {
            ++cost;
        }
    }
    
    printf("cost: %i\n", cost);
    return cost;
}

void mfl_reopt_postorder(node *n, int nchar)
{
    node *p;
    
    if (n->tip || n->clip) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_reopt_postorder(p->outedge, nchar);
        p = p->next;
    }
    mfl_subtree_count(n->next->outedge, n->next->next->outedge, n, nchar, NULL);
}

void mfl_reopt_preorder(node *n, int nchar)
{
    node *p, *dl, *dr;
    
    if (n->tip || n->clip) {
        return;
    }
    
    dl = n->next->outedge;
    dr = n->next->next->outedge;
    
    if (!dl->tip) {
        mfl_combine_up(dl, n, nchar);
    }
    if (!dr->tip) {
        mfl_combine_up(dr, n, nchar);
    }
    
    p = n->next;
    while (p != n) {
        mfl_reopt_preorder(p->outedge, nchar);
        p = p->next;
    }
}

void mfl_subtree_count(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *trlength)
{
    int i;
    charstate lft_chars, rt_chars;
    
    for (i = 0; i < nchar; ++i) {
        if (leftdesc->tempapos[i] & rightdesc->tempapos[i]) 
        {
            ancestor->tempapos[i] = leftdesc->tempapos[i] & rightdesc->tempapos[i];
        }
        else
        {
            lft_chars = leftdesc->tempapos[i];
            rt_chars = rightdesc->tempapos[i];
            ancestor->tempapos[i] = lft_chars | rt_chars;
            if ((ancestor->tempapos[i] & IS_APPLIC) && 
                ( ((ancestor->tempapos[i] & IS_APPLIC) & lft_chars) && 
                 ((ancestor->tempapos[i] & IS_APPLIC) & rt_chars) )) {
                    ancestor->tempapos[i] = (ancestor->tempapos[i] & IS_APPLIC);
                    if (trlength) {
                        *trlength = *trlength + 1;
                    }
                }
        }
    }
}

int mfl_get_sttreelen(tree *testtree, charstate *tipdata, int ntax, int nchar, int *besttreelen)
{
    int treelen = 0;
    int *treelen_p = &treelen;
    
    mfl_apply_tipdata(testtree, tipdata, ntax, nchar);
    mfl_subtree_postorder(testtree->root, treelen_p, nchar);
    mfl_fitch_preorder(testtree->root, nchar);
    
    return *treelen_p;
}

void mfl_countsteps(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *trlength, int *besttreelen)
{
    int i;
    charstate lft_chars, rt_chars;
    /*if (ancestor->tip) {
        return;
    }*/
    
    for (i = 0; i < nchar; ++i) {
        if (leftdesc->tempapos[i] & rightdesc->tempapos[i]) 
        {
            ancestor->tempapos[i] = leftdesc->tempapos[i] & rightdesc->tempapos[i];
        }
        else
        {
            lft_chars = leftdesc->tempapos[i];
            rt_chars = rightdesc->tempapos[i];
            ancestor->tempapos[i] = lft_chars | rt_chars;
            if ((ancestor->tempapos[i] & IS_APPLIC) && ( ((ancestor->tempapos[i] & IS_APPLIC) & lft_chars) && ((ancestor->tempapos[i] & IS_APPLIC) & rt_chars) )) 
            {
                ancestor->tempapos[i] = (ancestor->tempapos[i] & IS_APPLIC);
                *trlength = *trlength + 1;
                /*if (*trlength > *besttreelen) 
                {
                    return;
                }*/
            }
        }
    }
}

void mfl_subtree_postorder(node *n, int *trlength, int nchar)
{
    node *p;
    
    if (n->tip) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_subtree_postorder(p->outedge, trlength, nchar);
        p = p->next;
    }
    if (!n->apomorphies) {
        n->apomorphies = (charstate*)malloc(nchar * sizeof(charstate));
    }
    if (!n->tempapos) {
        n->tempapos = (charstate*)malloc(nchar * sizeof(charstate));
    }
    /*if (!n->outedge) {
        n->apomorphies = n->tempapos;
    }*/
    mfl_subtree_count(n->next->outedge, n->next->next->outedge, n, nchar, trlength);
}

void mfl_fitch_postorder(node *n, int *trlength, int nchar, int *besttreelen)
{
    node *p;
    
    if (n->tip) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_fitch_postorder(p->outedge, trlength, nchar, besttreelen);
        p = p->next;
    }
    
    if (!n->apomorphies) {
        n->apomorphies = (charstate*)malloc(nchar * sizeof(charstate));
    }
    if (!n->tempapos) {
        n->tempapos = (charstate*)malloc(nchar * sizeof(charstate));
    }

    mfl_countsteps(n->next->outedge, n->next->next->outedge, n, nchar, trlength, besttreelen);
    n->nodelen = *trlength;
}

void mfl_fitch_allviews(node *n, int *trlength, int nchar, int *besttreelen)
{
    node *p;
    
    if (n->tip) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_fitch_allviews(p->outedge, trlength, nchar, besttreelen);
        p = p->next;
    }
    
    if (!n->apomorphies) {
        n->apomorphies = (charstate*)malloc(nchar * sizeof(charstate));
    }
    if (!n->tempapos) {
        n->tempapos = (charstate*)malloc(nchar * sizeof(charstate));
    }

    mfl_countsteps(n->next->outedge, n->next->next->outedge, n, nchar, trlength, besttreelen);
    n->nodelen = *trlength;
}

void mfl_combine_up(node *n, node *anc, int nchar)
{
    int i;
    charstate lft_chars, rt_chars;
    
    if (n->visited) {
        return;
    }
    
    for (i = 0; i < nchar; ++i) {
        
        if ((n->tempapos[i] & anc->apomorphies[i]) == anc->apomorphies[i]) 
        {
            n->apomorphies[i] = n->tempapos[i] & anc->apomorphies[i];            
        }
        else {
            lft_chars = n->next->outedge->tempapos[i];
            rt_chars = n->next->next->outedge->tempapos[i];
            if ( lft_chars & rt_chars ) { //III
                //V
                n->apomorphies[i] = (n->tempapos[i] |(anc->apomorphies[i] &(lft_chars | rt_chars))) & IS_APPLIC;
            }
            else {
                //IV
                n->apomorphies[i] = n->tempapos[i] | anc->apomorphies[i];
            }
        }
    }
    
    if (n->index == 14) {
        printf("state after: %i\n", n->apomorphies[20]);
    }
    n->visited = 1;
}

void mfl_fitch_preorder(node *n, int nchar)
{
    node *p, *dl, *dr;
    
    if (n->tip) {
        return;
    }
    
    dl = n->next->outedge;
    dr = n->next->next->outedge;
    
    if (!dl->tip) {
        mfl_combine_up(dl, n, nchar);
    }
    if (!dr->tip) {
        mfl_combine_up(dr, n, nchar);
    }
    
    p = n->next;
    while (p != n) {
        mfl_fitch_preorder(p->outedge, nchar);
        p = p->next;
    }
}

int mfl_get_subtreelen(node *n, int ntax, int nchar, int *besttreelen)
{
    int treelen = 0;
    int *treelen_p = &treelen;
    
    mfl_fitch_postorder(n, treelen_p, nchar, besttreelen);
    mfl_fitch_preorder(n, nchar);
    
    return *treelen_p;
}

int mfl_get_trlonepass(tree *testtree, charstate *tipdata, int ntax, int nchar, int *besttreelen)
{
    int treelen = 0;
    int *treelen_p = &treelen;
    
    mfl_fitch_postorder(testtree->root, treelen_p, nchar, besttreelen);
    
    return *treelen_p;
}

int mfl_get_treelen(tree *t, int ntax, int nchar, int *besttreelen)
{
    int i;
    int treelen = 0;
    int *treelen_p = &treelen;
    
    if (!t->root) {
        mfl_temproot(t, 1, ntax);
        mfl_fitch_postorder(t->root, treelen_p, nchar, besttreelen);
        mfl_fitch_preorder(t->root, nchar);
        mfl_undo_temproot(ntax, t);
    }
    else {
        mfl_fitch_postorder(t->root, treelen_p, nchar, besttreelen);
        mfl_fitch_preorder(t->root, nchar);
    }

    return *treelen_p;
}

void mfl_compute_allviews(node *n, tree *t, int ntax, int nchar, int *treelen, int *besttreelen)
{
    int j;
    node *p;
    
    joinNodes(n->outedge, t->trnodes[ntax]->next->next);
    joinNodes(n, t->trnodes[ntax]->next);
    t->root = t->trnodes[ntax];
    mfl_reopt_postorder(t->root, nchar);
    for (j = 0; j < nchar; ++j) {
        t->root->apomorphies[j] = (t->root->next->outedge->tempapos[j] | t->root->next->next->outedge->tempapos[j]);
    }
    t->root->visited = 0;
    mfl_reopt_preorder(t->root, nchar);
    t->root->visited = 0;
    joinNodes(t->trnodes[ntax]->next->next->outedge, n);
    t->root = NULL;
    
    
    /*if (n->start) {
        mfl_compute_allviews(n->outedge, t, ntax, nchar, treelen, besttreelen);
        return;
    }*/
    
    if (n->tip) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_compute_allviews(p->outedge, t, ntax, nchar, treelen, besttreelen);
        p = p->next;
    }
    
}

int mfl_all_views(tree *t, int ntax, int nchar, int *besttreelen)
{
    int i, j;
    int treelen = 0, fptreelen;
    int *treelen_p = &treelen;
    
    //mfl_compute_allviews(t->trnodes[0], t, ntax, nchar, treelen_p, besttreelen);
    
    mfl_devisit_tree(t->trnodes, 2 * ntax - 1);
    
    //int pos = 4;
    
    for (i = 0; i < ntax; ++i) {
        if (i != ntax) {
            *treelen_p = 0;
            mfl_temproot(t, i, ntax);
            mfl_fitch_allviews(t->root, treelen_p, nchar, besttreelen);
            for (j = 0; j < nchar; ++j) {
                if (t->root->next->outedge->tempapos[j] & t->root->next->next->outedge->tempapos[j]) {
                    t->root->apomorphies[j] = (t->root->next->outedge->tempapos[j] & t->root->next->next->outedge->tempapos[j]);
                }
                else {
                    t->root->apomorphies[j] = (t->root->next->outedge->tempapos[j] | t->root->next->next->outedge->tempapos[j]);
                }

            }
            t->root->visited = 0;
            mfl_fitch_preorder(t->root, nchar);
            t->root->visited = 0;
            mfl_undo_temproot(ntax, t);
            if (i == 0) {
                fptreelen = *treelen_p;
            }
            
        }
    }   
    return fptreelen;
}
