/*
 *  compare.c
 *  Morphy
 *
 *  Created by Martin Brazeau and Chris Desjardins on 2/8/12.
 *  Copyright 2012. All rights reserved.
 *
 *  
 *
 */

#include "morphy.h"
#include <time.h>



bool mfl_compare_trees(int *t1, int *t2, int ntax)
{
    /* could probably change these errors/exits to some kind of exception thing. I'm just not sure how */
    //assert(t1 != NULL);
    
    //assert(t2 != NULL);
    
    int pos = 2 * ntax - 1;
    
    return memcmp(t1, t2, pos * sizeof(int));
    
}

void mfl_clear_cpindex(tree *t, int numnodes)
{
    int i;
    node *p;
    for (i = 0; i < numnodes; ++i) {
        t->trnodes[i]->cpindex = 0;
        if (t->trnodes[i]->next) {
            p = t->trnodes[i]->next;
            while (p != t->trnodes[i]) {
                p->cpindex = 0;
                p = p->next;
            }
        }
    }
}

int *mfl_compress_tree(tree *t, int ntax, int numnodes)
{
    int i, j;
    node *p, *n;
    bool wasrooted = true;
    
    int *storedtr = (int*)malloc(numnodes * sizeof(int));
    memset(storedtr, 0, numnodes * sizeof(int));
    
    if (!t->root) {
        mfl_root_tree(t, 0, ntax);
        wasrooted = false;
    }
    
    mfl_clear_cpindex(t, numnodes);
    
    t->root->cpindex = ntax;
    t->root->bottom = true;
    
    j = ntax+1;
    
    for (i = 0; i < ntax; ++i) {        
        
        p = t->trnodes[i]->edge;
        n = t->trnodes[i];
        n->cpindex = i;
        
        while (!p->cpindex) {
            p = p->next;
            if (p->bottom) {
                if (!p->cpindex) {
                    p->cpindex = j;
                    ++j;
                    storedtr[n->cpindex] = p->cpindex;                    
                    n = p;
                    p = p->edge;

                }
                else if (p->cpindex) {
                    storedtr[n->cpindex] = p->cpindex;
                    n = p;
                }
            }
        }
    }
    
    if (!wasrooted) {
        mfl_unroot(ntax, t);
    }
    
    //storedtr[numnodes] = '\0';
    
    /*for (i = 0; i <= numnodes; ++i) {
        dbg_printf("%i ", storedtr[i]);
    }
    dbg_printf("\n");*/
    
    return storedtr;
}

tree *mfl_decompress_tree(int *savedtr, int ntax, int numnodes)
{
    int i;
    tree *t;
    
    t = mfl_alloc_noring(ntax, numnodes);
    
    for (i = 0; i < numnodes; ++i) {
        if (i < ntax) {
            t->trnodes[i]->edge = mfl_allocnode();
            t->trnodes[i]->edge->edge = t->trnodes[i];
        }
    }
    
    /* Stuff will go here that reconstructs the tree from the integer array representation */
    
    return t;
}

bool mfl_compare_alltrees(tree *newtopol, tree **savedtrees, int ntax, int numnodes, long int *start, long int last)
{
    long int i = 0;
    int *newtr;
    bool foundtr = false;
    
    newtr = mfl_compress_tree(newtopol, ntax, numnodes);
    
    /* This search should be replaced by a binary search or some kind of hash 
     * function in order to speed it up. Otherwise it will be a major time hog 
     * if/when somebody loads a noisy dataset. */
    
    for (i = last; i > *start; --i) {
        if (!mfl_compare_trees(newtr, savedtrees[i-1]->compressedtr, ntax)) {
            foundtr = true;
            break;
        }
    }
    
    if (foundtr) {
        free(newtr);
    }
    else {
        newtopol->cmptrholder = newtr;
    }
    
    return foundtr;
}

