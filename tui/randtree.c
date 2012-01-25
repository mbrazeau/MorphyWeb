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

/*bool isodd (float seedn)
{
    if (seedn % 2 == 0)
    {
        return false;
    }
    else 
    {
        return true;
    }
    
}*/

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
    mfl_root_tree(randtree, 1, ntax);    // The second argument should eventually be replaced by a randomly drawn number between 0 and ntax-1
    
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

    taxarray = malloc(ntax * sizeof(int));
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

struct tree *mfl_stree_asis(int ntax, int numnodes, tree *treebuffer, charstate *tipdata)
{
    int i;
    node *p;
    tree *asistree;
    
    asistree = alloctree(ntax, numnodes);
    
    /* create the base star */
    /* The 'magic number' 3 appears here because that is all the value can ever
     * be. The smallest non-trivial rooted bifurcating tree has three leaves */
    p = asistree->trnodes[ntax + 1];
    for (i = 0; i < 3; ++i) {
        joinNodes(asistree->trnodes[0], p);
        p = p->next;
    }
    
    for (i = 3; i < ntax; ++i) {
        printf("something will happen here\n");
    }
    
    return asistree;
}