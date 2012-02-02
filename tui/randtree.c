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

bool mfl_headtail (void)
{
    if (rand()% 2 == 0)
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

void shuffle(int *taxarray, int ntax)
{
    int i, j, t;
    
    
    if (ntax > 1) {
        for (i = 0; i < ntax - 1; i++) {
            j = i + rand() / (RAND_MAX / (ntax - i) + 1);
            t = taxarray[j];
            taxarray[j] = taxarray[i];
            taxarray[i] = t;
        }
    }
}

void shuffle_w_outgroup(int *taxarray, int ntax, int outtaxa[])
{
    int i, j, t;
    
    
    if (ntax > 1) {
        for (i = 0; i < ntax - 1; i++) {
            j = i + rand() / (RAND_MAX / (ntax - i) + 1);
            t = taxarray[j];
            taxarray[j] = taxarray[i];
            taxarray[i] = t;
        }
    }
}

struct tree *randtrunk(tree *newtrunk, node *startn)
{
    /* Generates the 'trunk' of the tree which consists of the network of all
     * internal nodes in the tree */
    
    return (newtrunk);
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
    init_taxarray(taxarray, ntax);
    
    shuffle(taxarray, ntax);
    
    randtree = alloctree(ntax, numnodes);
    
    randtree->trnodes[0]->start = true;
    randtree->trnodes[0]->outedge = randtree->trnodes[taxarray[0]];
    
    // Join all the internal nodes (except the root) together
    for (i = 1; i <= (ntax - 3); ++i) {
        p = randtree->trnodes[ntax + i]->next->next;
        q = randtree->trnodes[ntax + i + 1];
        joinNodes(p, q);
    }
    
    // Add all the tips to the appropriate internal nodes
    
    joinNodes(randtree->trnodes[ntax + 1], randtree->trnodes[taxarray[0] - 1]);
    
    for (i = 1; i < ntax - 1; ++i) {
        randtree->trnodes[taxarray[i] - 1]->outedge = randtree->trnodes[ntax + i]->next;
        randtree->trnodes[ntax + i]->next->outedge = randtree->trnodes[taxarray[i] - 1];
    } 
    
    randtree->trnodes[2 * ntax - 2]->next->next->outedge = randtree->trnodes[taxarray[i] - 1];
    randtree->trnodes[taxarray[i] - 1]->outedge = randtree->trnodes[2 * ntax - 2]->next->next;
    
    free(taxarray);
    
    printNewick(randtree->trnodes[0]);
    printf(";\n");
    
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
    
    if (n->outedge) {
        mfl_insert_branch(newbranch, n, ntax);
        trlen = mfl_get_treelen(starttree, tipdata, ntax, nchar);
        starttree->length = trlen;
        //printf("tested length: %i\n", *bestlen);
        if (trlen < *bestlen) {
            //printf("Setting bestpos to current node\n");
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
        bestpos = mfl_tryall(p->outedge, newbranch, bestpos, ntax, nchar, numnodes, 
                   bestlen, starttree, savedtrees, tipdata);
        p = p->next;
    }
    
    return bestpos;
}

void mfl_addseq_randasis(int ntax, int nchar, int numnodes,
                                 charstate *tipdata, bool addRandom,
                                 tree **savedtrees)
{
    int i, nbeslen;
    int *bestlen = &nbeslen;
    int *taxarray;
    node *p, *bestpos;
    
    savedtrees[0] = alloc_noring(ntax, numnodes);
    
    taxarray = (int*)malloc(ntax * sizeof(int));
    memset(taxarray, 0, ntax * sizeof(int));  // This is to see if I can fix the problem
    init_taxarray(taxarray, ntax);
    
    
    if (addRandom) {
        printf("Joining taxa by random addition sequence\n");
        shuffle(taxarray, ntax);
    }
    else {
        printf("Joining taxa according to order in matrix\n");
    }
    
    p = savedtrees[0]->trnodes[ntax + 1];
    newring(p, ntax);
    i = 0;
    do {
        joinNodes(savedtrees[0]->trnodes[taxarray[i] - 1], p);
        p = p->next;
        ++i;
    } while (p != savedtrees[0]->trnodes[ntax + 1]);

    mfl_temproot(savedtrees[0], taxarray[0] - 1, ntax);
    savedtrees[0]->length = *bestlen;
    bestpos = savedtrees[0]->trnodes[taxarray[0] - 1];
    
    for (i = 3; i < ntax; ++i) {
        newring(savedtrees[0]->trnodes[ntax + i - 1], ntax);
        joinNodes(savedtrees[0]->trnodes[taxarray[i] - 1], savedtrees[0]->trnodes[ntax + i - 1]->next);
        mfl_insert_branch(savedtrees[0]->trnodes[taxarray[i] - 1], savedtrees[0]->trnodes[taxarray[0] - 1], ntax);
        *bestlen = mfl_get_treelen(savedtrees[0], tipdata, ntax, nchar);
        //printf("Preliminary length: %i\n", *bestlen);
        mfl_remove_branch(savedtrees[0]->trnodes[taxarray[i] - 1]);
        bestpos = mfl_tryall(savedtrees[0]->root, savedtrees[0]->trnodes[taxarray[i] - 1], bestpos, ntax, nchar, 
                   numnodes, bestlen, savedtrees[0], savedtrees, tipdata);
        //Join the nodes//
        mfl_insert_branch(savedtrees[0]->trnodes[taxarray[i] - 1], bestpos, ntax);
        *bestlen = 0;
    }
    mfl_undo_temproot(ntax, savedtrees[0]);
    
    free(taxarray);
}