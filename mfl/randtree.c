/*
 *  random.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 11-07-25.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 *  This file contains code for generating random rooted (pre-specified) and unrooted trees
 *
 */

#include "morphy.h"

void mfl_init_taxarray(int *taxarray, int ntax)
{
    int i;
    
    for (i = 0; i < ntax; ++i) {
        taxarray[i] = i + 1;
    }
}

bool mfl_headtail (void)
{
    if (random() % 2 == 0)
    {
        return false;
    }
    else 
    {
        return true;
    }
}

/* Taxon shuffle using array shuffle by Ben Pfaff http://benpfaff.org/writings/clc/shuffle.html 
 * ****NOTE******: Will be modified to employ the GSL random number generator*/

void mfl_shuffle(int *taxarray, int ntax)
{
    int i, j, t;
    
    
    if (ntax > 1) {
        for (i = 0; i < ntax - 1; i++) {
            j = i + random() / (RAND_MAX / (ntax - i) + 1);
            t = taxarray[j];
            taxarray[j] = taxarray[i];
            taxarray[i] = t;
        }
    }
}

struct tree *randrooted(int ntax, int numnodes)
{
    /* Returns a random tree with an arbitrary root*/
    
    tree *randtree;
    randtree = randunrooted(ntax, numnodes);
    randtree->trnodes[0]->start = false;
    mfl_root_tree(randtree, 1, ntax);    // The second argument should 
                                         // eventually be replaced by a randomly 
                                         // drawn number between 0 and ntax-1
    return (randtree);
}

//struct tree *rand_w_outgroup (int root)
//{
    /* Returns a random ingroup topology on a given root
     * and will arbitrarily resolve the outgroup*/
//  return;
    
//}

struct tree *randunrooted(int ntax, int numnodes)
{
    /* Returns a random unrooted tree*/
    
    int i;
    int *taxarray;
    node *p, *q;
    tree *randtree;

    taxarray =(int*) malloc(ntax * sizeof(int));
    mfl_init_taxarray(taxarray, ntax);
    
    mfl_shuffle(taxarray, ntax);
    
    randtree = mfl_alloctree(ntax, numnodes);
    
    randtree->trnodes[0]->start = true;
    randtree->trnodes[0]->edge = randtree->trnodes[taxarray[0]];
    
    // Join all the internal nodes (except the root) together
    for (i = 1; i <= (ntax - 3); ++i) {
        p = randtree->trnodes[ntax + i]->next->next;
        q = randtree->trnodes[ntax + i + 1];
        mfl_join_nodes(p, q);
    }
    
    // Add all the tips to the appropriate internal nodes
    
    mfl_join_nodes(randtree->trnodes[ntax + 1], randtree->trnodes[taxarray[0] - 1]);
    
    for (i = 1; i < ntax - 1; ++i) {
        randtree->trnodes[taxarray[i] - 1]->edge = randtree->trnodes[ntax + i]->next;
        randtree->trnodes[ntax + i]->next->edge = randtree->trnodes[taxarray[i] - 1];
    } 
    
    randtree->trnodes[2 * ntax - 2]->next->next->edge = randtree->trnodes[taxarray[i] - 1];
    randtree->trnodes[taxarray[i] - 1]->edge = randtree->trnodes[2 * ntax - 2]->next->next;
    
    free(taxarray);
    
    /*printNewick(randtree->trnodes[0]);*/
    dbg_printf(";\n");
    
    return (randtree);
}

struct node *mfl_breaktie_intryall(node *bestpos, node *n)
{
    if (mfl_headtail()) {
        bestpos = n;
    }
    
    return bestpos;
}

struct node * mfl_tryall(node *n, node *newbranch, node *bestpos, int ntax, int nchar, 
                int numnodes, int *bestlen, tree *starttree, tree **savedtrees, charstate *tipdata)
{
    int trlen = 0;
    node *p;
    
    if (n->edge) {
        mfl_insert_branch(newbranch, n, ntax);
        trlen = mfl_get_sttreelen(starttree, tipdata, ntax, nchar, bestlen);
        starttree->length = trlen;
        //dbg_printf("tested length: %i\n", *bestlen);
        if (trlen < *bestlen) {
            //dbg_printf("Setting bestpos to current node\n");
            bestpos = n;
            *bestlen = trlen;
        } 
        if (trlen == *bestlen) 
        {
            bestpos = mfl_breaktie_intryall(n, bestpos);
        }
        mfl_remove_branch(newbranch);
    }
    
    if (n->tip) {
        return bestpos;
    }    
    
    p = n->next;
    while (p != n) {
        bestpos = mfl_tryall(p->edge, newbranch, bestpos, ntax, nchar, numnodes, 
                   bestlen, starttree, savedtrees, tipdata);
        p = p->next;
    }
    
    return bestpos;
}

tree *mfl_addseq_randasis(int ntax, int nchar, int numnodes,
                                 charstate *tipdata, bool addRandom,
                                 tree **savedtrees)
{
    int i, nbeslen = 0;
    int *bestlen = &nbeslen;
    int *taxarray;
    node *p, *bestpos;
    
    tree *newtree = mfl_alloc_noring(ntax, numnodes);
    
    taxarray = (int*)malloc(ntax * sizeof(int));
    memset(taxarray, 0, ntax * sizeof(int));  // This is to see if I can fix the problem
    mfl_init_taxarray(taxarray, ntax);
    
    
    if (addRandom) {
        dbg_printf("Joining taxa by random addition sequence\n");
        mfl_shuffle(taxarray, ntax);
    }
    else {
        dbg_printf("Joining taxa according to order in matrix\n");
    }
    
    p = newtree->trnodes[ntax + 1];
    mfl_newring(p, ntax);
    i = 0;
    do {
        mfl_join_nodes(newtree->trnodes[taxarray[i] - 1], p);
        p = p->next;
        ++i;
    } while (p != newtree->trnodes[ntax + 1]);

    mfl_temproot(newtree, taxarray[0] - 1, ntax);
    newtree->length = *bestlen;
    bestpos = newtree->trnodes[taxarray[0] - 1];
    
    for (i = 3; i < ntax; ++i) {
        mfl_newring(newtree->trnodes[ntax + i - 1], ntax);
        mfl_join_nodes(newtree->trnodes[taxarray[i] - 1], newtree->trnodes[ntax + i - 1]->next);
        mfl_insert_branch(newtree->trnodes[taxarray[i] - 1], newtree->trnodes[taxarray[0] - 1], ntax);
        *bestlen = mfl_get_sttreelen(newtree, tipdata, ntax, nchar, bestlen);
        //dbg_printf("Preliminary length: %i\n", *bestlen);
        mfl_remove_branch(newtree->trnodes[taxarray[i] - 1]);
        bestpos = mfl_tryall(newtree->root, newtree->trnodes[taxarray[i] - 1], bestpos, ntax, nchar, 
                   numnodes, bestlen, newtree, savedtrees, tipdata);
        //Join the nodes//
        mfl_insert_branch(newtree->trnodes[taxarray[i] - 1], bestpos, ntax);
        *bestlen = 0;
    }
    mfl_undo_temproot(ntax, newtree);
    //newtree->bipartitions = mfl_tree_biparts(newtree, ntax, numnodes);
    
    free(taxarray);
    return newtree;
}
