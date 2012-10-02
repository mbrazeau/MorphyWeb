/*
 *  coptim.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 2/15/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "morphy.h"

charstate * mfl_convert_tipdata(char *txtsrc, int ntax, int nchar, bool na_as_missing)
{
    
    /* One of Morphy's most important and unique features will be distinguishing
     * between a missing data entry "?" and a character-inapplicability entry "-".
     * Most existing programs just treat "?" and "-" as identical values (-1).
     * In Morphy, they can be treated differently, such that ancestral states
     * are not reconstructed for taxa which cannot logically have them */
    
    int i, j;
    
    charstate *tipdata = (charstate*) malloc(ntax * nchar * sizeof(charstate));
    
    if (!na_as_missing) {
        dbg_printf("Gap symbol ('-') treated as character inapplicability\n");
    }
    else {
        dbg_printf("Gap symbol ('-') treated as missing data\n");
    }

    
    for (i = 0, j = 0; txtsrc[i]; ++i, ++j) {
        if ((txtsrc[i] - '0') >= 0 && (txtsrc[i] - '0') <= 9) {
            tipdata[j] = 1 << (txtsrc[i] - '0' + 1);
            
        }
        else if (txtsrc[i] == '{' || txtsrc[i] == '(') {
            ++i;
            tipdata[j] = 0;
            while (txtsrc[i] != '}' && txtsrc[i] != ')') {
                if ((txtsrc[i] - '0') >= 0 && (txtsrc[i] - '0') <= 9) {
                    tipdata[j] = tipdata[j] | (1 << (txtsrc[i] - '0' + 1));
                    ++i;
                }
            }
        }
        else if (txtsrc[i] == '?') {
            tipdata[j] = IS_APPLIC;
        }
        else if (txtsrc[i] == '-') {
            if (na_as_missing) {
                tipdata[j] = IS_APPLIC;//-1;
            }
            else {
                tipdata[j] = 1;
            }

        }
        else if (txtsrc[i] == '\n') {
            --j;
        }
        else {
            --j;
        }
    }
    
    return tipdata;
}

void mfl_apply_tipdata(tree *currenttree, charstate *tipdata, int ntax, int nchar)
{   
    int i, j;
    
    for (i = 0; i < ntax; ++i) {
        //currenttree->trnodes[i]->apomorphies = &tipdata[i * nchar];
        currenttree->trnodes[i]->tempapos = &tipdata[i * nchar];
        if (!currenttree->trnodes[i]->apomorphies) {
            currenttree->trnodes[i]->apomorphies = (charstate*)malloc(nchar * sizeof(charstate));
        }
        for (j = 0; j < nchar; ++j) {
            currenttree->trnodes[i]->apomorphies[j] = currenttree->trnodes[i]->tempapos[j];
        }
    }
}

int mfl_locreopt_cost(node *src, node *tgt1, node *tgt2, int nchar, int diff)
{
    /* Returns cost of inserting subtree src between tgt1 and tgt2 following
     * the algorithms described by Ronquist (1998. Cladistics) and Goloboff 
     * (1993, 1996. Cladistics).*/
    
    int i;
    int cost = 0;
    charstate *srctemps = src->apomorphies;
    charstate *tgt1apos = tgt1->apomorphies;
    charstate *tgt2apos = tgt2->apomorphies;
    
    for (i = 0; i < nchar; ++i) {
        if ( !(srctemps[i] & (tgt1apos[i] | tgt2apos[i])) ) {
            ++cost;
            if (cost > diff) {
                return cost;
            }
        }
    }
    return cost;
}

int mfl_subtr_reinsertion(node *src, node *tgt1, node *tgt2, int nchar)
{
    /* Returns cost of reinserting the subtree src at the original place from
     * which it was clipped (tgt1 and tgt2). This score is used to compute the 
     * difference in length between the two subtrees so that neither the total 
     * tree length nor the length of the individual subtrees needs to be 
     * calculated (Ronquist 1998). */
    
    int i;
    int cost = 0;
    charstate *srctemps = src->apomorphies;
    charstate *tgt1apos = tgt1->apomorphies;
    charstate *tgt2apos = tgt2->apomorphies;
    
    for (i = 0; i < nchar; ++i) {
        if ( !(srctemps[i] & (tgt1apos[i] | tgt2apos[i])) ) {
            ++cost;
        }
    }
    return cost;
}

void mfl_subtree_count_ii(node *leftdesc, node *rightdesc, node *ancestor, int nchar)
{
    int i;
    charstate lft_chars, rt_chars;
    
    for (i = 0; i < nchar; ++i) {
        lft_chars = leftdesc->apomorphies[i];
        rt_chars = rightdesc->apomorphies[i];
        ancestor->apomorphies[i] = lft_chars | rt_chars;
    }
}

void mfl_reopt_subtr_root(node *n, int nchar)
{
    if (!n->tip) {
        mfl_subtree_count_ii(n->next->outedge, n->next->next->outedge, n, nchar);
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

            if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                ancestor->tempapos[i] = ancestor->tempapos[i] & IS_APPLIC;
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

int mfl_wagner_count(charstate lchar, charstate rchar)
{
    int length = 0;
    
    if (lchar > rchar) {
        while (!(lchar & rchar)) {
            lchar = lchar >> 1;
            ++length;
        }
    }
    else {
        while (!(lchar & rchar)) {
            rchar = rchar >> 1;
            ++length;
        }
    }
    return length;
}

void mfl_countsteps(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *trlength, int *besttreelen)
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

            if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                ancestor->tempapos[i] = ancestor->tempapos[i] & IS_APPLIC;
                *trlength = *trlength + 1;
            }
        }
    }
}

void mfl_combine_up(node *n, node *anc, int nchar)
{
    int i;
    charstate lft_chars, rt_chars;
    
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
                n->apomorphies[i] = (n->tempapos[i] | (anc->apomorphies[i] & (lft_chars | rt_chars))); //& IS_APPLIC;
            }
            else {
                //IV
                if ( (anc->apomorphies[i] & IS_APPLIC) && (n->tempapos[i] & IS_APPLIC)) {
                    n->apomorphies[i] = (n->tempapos[i] | anc->apomorphies[i]) & IS_APPLIC;
                }
                else {
                    n->apomorphies[i] = n->tempapos[i] | anc->apomorphies[i];
                }

            }
        }
    }
}

void mfl_subtree_postorder(node *n, int *trlength, int nchar)
{
    node *p;
    
    if (!n->apomorphies) {
        n->apomorphies = (charstate*)malloc(nchar * sizeof(charstate));
        if (n->next) {
            mfl_join_apomorphies(n);
        }
    }
    
    if (n->tip) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_subtree_postorder(p->outedge, trlength, nchar);
        p = p->next;
    }
    if (!n->tempapos) {
        n->tempapos = (charstate*)malloc(nchar * sizeof(charstate));
    }
    mfl_subtree_count(n->next->outedge, n->next->next->outedge, n, nchar, trlength);
}

void mfl_fitch_postorder(node *n, int *trlength, int nchar, int *besttreelen)
{
    node *p;
    
    if (!n->apomorphies) {
        n->apomorphies = (charstate*)malloc(nchar * sizeof(charstate));
        if (n->next) {
            mfl_join_apomorphies(n);
        }
    }
    if (!n->tempapos) {
        n->tempapos = (charstate*)malloc(nchar * sizeof(charstate));
    }
    
    if (n->tip) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_fitch_postorder(p->outedge, trlength, nchar, besttreelen);
        p = p->next;
    }

    mfl_countsteps(n->next->outedge, n->next->next->outedge, n, nchar, trlength, besttreelen);
    n->nodelen = *trlength;
}

void mfl_reopt_fitch(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *trlength)
{
    int i;
    charstate lft_chars, rt_chars, temp;
    bool allsame = false;
    
    charstate *ldtemps = leftdesc->tempapos;
    charstate *rdtemps = rightdesc->tempapos;
    charstate *antemps = ancestor->tempapos;
    
    for (i = 0; i < nchar; ++i) {
        if (ldtemps[i] & rdtemps[i]) 
        {
            temp = ldtemps[i] & rdtemps[i];
            if (temp != antemps[i]) {
                antemps[i] = temp;
                allsame = false;
            }
        }
        else
        {
            lft_chars = ldtemps[i];
            rt_chars = rdtemps[i];

            temp = lft_chars | rt_chars;

            if (lft_chars & IS_APPLIC && rt_chars & IS_APPLIC) {
                temp = temp & IS_APPLIC;
            }
            
            if (temp != antemps[i]) {
                antemps[i] = temp;
                allsame = false;
            }
        }
    }
    if (allsame) {
        ancestor->success = true;
    }
}

void mfl_reopt_postorder(node *n, int nchar)
{
    node *p;
    bool allsame = true;
        
    if (n->tip) {
        return;
    }
    
    if (n->clip) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_reopt_postorder(p->outedge, nchar);
        if (!p->outedge->success) {
            allsame = false;
        }
        //p->outedge->success = false;
        p = p->next;
    }
    
    if (!allsame) {
        mfl_reopt_fitch(n->next->outedge, n->next->next->outedge, n, nchar, NULL);
    }
    
}

void mfl_fitch_allviews(node *n, int *trlength, int nchar, int *besttreelen)
{
    node *p;
    int weight = 0;
    
    if (!n->apomorphies) {
        n->apomorphies = (charstate*)malloc(nchar * sizeof(charstate));
        if (n->next) {
            mfl_join_apomorphies(n);
        }
    }
    
    if (n->tip) {// || n->visited) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_fitch_allviews(p->outedge, trlength, nchar, besttreelen);
        weight = weight + p->outedge->vweight;
        p = p->next;
    }
    
    if (!n->tempapos) {
        n->tempapos = (charstate*)malloc(nchar * sizeof(charstate));
    }

    mfl_countsteps(n->next->outedge, n->next->next->outedge, n, nchar, trlength, besttreelen);
    n->visited = 1;
    n->vweight = weight;
}

void mfl_tip_apomorphies(node *tip, node *anc, int nchar)
{
    /* Reconstructs the tip set if it is polymorphic */
    
    int i;
    charstate *tiptemp = tip->tempapos;
    charstate *tipapos = tip->apomorphies;
    charstate *ancapos = anc->apomorphies;
    
    for (i = 0; i < nchar; ++i) {
        
        if (tiptemp[i] != 1) {
            if (tiptemp[i] & ancapos[i]) {
                tipapos[i] = tiptemp[i] & ancapos[i];
            }
        }
    }
}

void mfl_tip_reopt(tree *t, int ntax, int nchar)
{
    int i;
    for (i = 0; i < ntax; ++i) {
        mfl_tip_apomorphies(t->trnodes[i], t->trnodes[i]->outedge, nchar);
    }
}

void mfl_reopt_comb(node *n, node *anc, int nchar)
{
    int i;
    charstate lft_chars, rt_chars;
    charstate temp;
    charstate *ntemps = n->tempapos;
    charstate *napos = n->apomorphies;
    charstate *ancapos = anc->apomorphies;
    
    bool allsame = true;
    
    for (i = 0; i < nchar; ++i) {
        
        if ((ntemps[i] & ancapos[i]) == ancapos[i]) 
        {
            temp = ntemps[i] & ancapos[i]; 
            if (temp != napos[i]) {
                napos[i] = temp;
                allsame = false;
            }
        }
        else {
            lft_chars = n->next->outedge->tempapos[i];
            rt_chars = n->next->next->outedge->tempapos[i];
            
            if ( lft_chars & rt_chars ) { //III
                //V
                temp = (ntemps[i] |(ancapos[i] &(lft_chars | rt_chars)));// & IS_APPLIC;
                if (temp != napos[i]) {
                    napos[i] = temp;
                    allsame = false;
                }
            }
            else {
                //IV
                if ( (ancapos[i] & IS_APPLIC) && (ntemps[i] & IS_APPLIC)) {
                    temp = (ntemps[i] | ancapos[i]) & IS_APPLIC;
                }
                else {
                    temp = ntemps[i] | ancapos[i];
                }
                
                if (temp != napos[i]) {
                    napos[i] = temp;
                    allsame = false;
                }
            }
        }
    }
    if (allsame) {
        n->success = true;
        //n->finished = true;
    }
}


void mfl_set_rootstates(node *n, int nchar)
{
    int i;
    
    for (i = 0; i < nchar; ++i) {
        if (n->next->outedge->tempapos[i] & n->next->next->outedge->tempapos[i]) {
            n->apomorphies[i] = (n->next->outedge->tempapos[i] & n->next->next->outedge->tempapos[i]);
        }
        else {
            n->apomorphies[i] = (n->next->outedge->tempapos[i] | n->next->next->outedge->tempapos[i]);
            if (n->next->outedge->tempapos[i] & IS_APPLIC && n->next->next->outedge->tempapos[i] & IS_APPLIC) {
                n->apomorphies[i] = n->apomorphies[i] & IS_APPLIC;
            }
        }
    }
}

void mfl_reopt_preorder(node *n, int nchar)
{
    
    node *dl, *dr;
    
    if (n->tip) {
        return;
    }
    
    if (!n->outedge || n->isroot) {
        mfl_set_rootstates(n, nchar);
    }
    
    dl = n->next->outedge;
    dr = n->next->next->outedge;
    
    if (!dl->tip) {
        mfl_reopt_comb(dl, n, nchar);
    }
    else if (!dl->success) {
        mfl_tip_apomorphies(dl, n, nchar);
    }
    
    if (!dr->tip) {
        mfl_reopt_comb(dr, n, nchar);
    }
    else if (!dr->success) {
        mfl_tip_apomorphies(dr, n, nchar);
    }
    
    //n->success = false;
    
    mfl_reopt_preorder(n->next->outedge, nchar);
    mfl_reopt_preorder(n->next->next->outedge, nchar);
    
}

void mfl_fitch_preorder(node *n, int nchar)
{
    node *p, *dl, *dr;
    
    if (n->tip) {// || n->clip) {
        return;
    }
    
    if (!n->outedge || n->isroot) {
        mfl_set_rootstates(n, nchar);
    }
    
    dl = n->next->outedge;
    dr = n->next->next->outedge;
    
    if (!dl->tip) {
        mfl_reopt_comb(dl, n, nchar);
    }
    else {
        mfl_tip_apomorphies(dl, n, nchar);
    }
    
    if (!dr->tip) {
        mfl_reopt_comb(dr, n, nchar);
    }
    else {
        mfl_tip_apomorphies(dr, n, nchar);
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

int mfl_get_treelen(tree *t, int ntax, int nchar, int *besttreelen)
{
    int treelen = 0;
    int *treelen_p = &treelen;
    
    if (!t->root) {
        mfl_temproot(t, 0, ntax);
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

charstate *mfl_get_changing(node *n, int nchar)
{
    
    int i, j;
    
    charstate *prelims = n->tempapos;
    charstate *finals = n->apomorphies;
    charstate *changing = (charstate*)malloc((nchar + 1) * sizeof(charstate));
    
    j = 0;
    for (i = 0; i < nchar; ++i) {
        if (prelims[i] != finals[i]) {
            changing[j] = i;
            ++j;
        }
    }
    changing[j] = '\0';
    
    return changing;
}

void mfl_wipe_states(node *n, int nchar)
{
    memset(n->apomorphies, 0, nchar * sizeof(charstate));
    memset(n->tempapos, 0, nchar * sizeof(charstate));
}

void mfl_allviews_traversal(node *n, tree *t, int ntax, int nchar, int *treelen, int *besttreelen)
{

    
    if (n->start) {
        mfl_allviews_traversal(n->outedge, t, ntax, nchar, treelen, besttreelen);
        return;
    }
    
    
    if (n->tip || n->outedge->tip) {
    
        //dbg_printf("working on node %i\n", n->index);
        t->root = NULL;
        
        mfl_join_nodes(n->outedge, t->trnodes[ntax]->next->next);
        mfl_join_nodes(n, t->trnodes[ntax]->next);
        
        t->root = t->trnodes[ntax];
        
        mfl_reopt_postorder(t->root, nchar);
        
        t->root->visited = 0;
        
        mfl_join_nodes(t->trnodes[ntax]->next->next->outedge, n);
        
        if (n->tip) {
            return;
        }
    }

    mfl_allviews_traversal(n->next->outedge, t, ntax, nchar, treelen, besttreelen);
    mfl_allviews_traversal(n->next->next->outedge, t, ntax, nchar, treelen, besttreelen);
}

void mfl_subtr_allviews(node *base, tree *t, int nchar, charstate *changing)
{
    
    int numnodes = 2 * 47 -1;
    
    /* For subtree reoptimization only. */
    base->isroot = true;
    mfl_desuccess_tree(t, numnodes);
    //mfl_reopt_postorder_ii(base, nchar, changing);
    mfl_reopt_postorder(base, nchar);
    mfl_reopt_preorder(base, nchar);
    base->isroot = false;
    
}

void mfl_trav_allviews(node *n, tree *t, int ntax, int nchar, int *treelen, int *besttreelen)
{
    
    /* For subtree reoptimization only. */
    
    mfl_definish_tree(t, 2 * ntax - 1);
    mfl_temproot(t, 0, ntax);
    mfl_reopt_postorder(t->root, nchar);
    mfl_reopt_preorder(t->root, nchar);
    mfl_undo_temproot(ntax, t);
    //mfl_tip_reopt(t, ntax, nchar);
    
}

int mfl_all_views(tree *t, int ntax, int nchar, int *besttreelen)
{
    int i=0;
    int treelen = 0, fptreelen;
    int *treelen_p = &treelen;
    
    mfl_devisit_tree(t->trnodes, 2 * ntax - 1);
    mfl_definish_tree(t, 2 * ntax - 1);

    t->root = t->trnodes[ntax];
    for (i = 0; i < ntax; ++i) {
        *treelen_p = 0;
        mfl_temproot(t, i, ntax);
        mfl_fitch_allviews(t->root, treelen_p, nchar, besttreelen);
        t->root->visited = 0;
        mfl_undo_temproot(ntax, t);
        if (i == 0) {
            fptreelen = *treelen_p;
        }
    }
    
    mfl_temproot(t, 0, ntax);
    t->root->visited = 0;
    mfl_fitch_preorder(t->root, nchar);
    mfl_undo_temproot(ntax, t);
    
    return fptreelen;
}
