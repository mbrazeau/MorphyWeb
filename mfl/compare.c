/*
 *  compare.c
 *  Morphy
 *
 *  Created by Martin Brazeau and Chris Desjardins on 2/8/12.
 *  Copyright 2012. All rights reserved.
 *
 *  Tree comparison functions. These break down the tree into bipartition hash
 *  tables and search for tree equality by comparing the hash tables. This is
 *  used to check whether a new tree is already found in memory. 
 *
 *  The routine is both slow and crude, and probably has polynomial time 
 *  complexity. For now, it suffices for testing purposes.
 *
 */

#include "morphy.h"
#include <time.h>



bool mfl_compare_trees(int *t1, int *t2, int numnodes)
{
    /* could probably change these errors/exits to some kind of exception thing. I'm just not sure how */
    if (t1 == NULL) {
        dbg_printf("error in mfl_compare_trees: tree 1 is invalid\n");
        exit(1);
    }
    
    if (t2 == NULL) {
        dbg_printf("error in mfl_compare_trees: tree 2 is invalid\n");
        exit(2);
    }
    
    if (memcmp(t1, t2, numnodes * sizeof(int))) {
        return false;
    }
    else {
        return true;
    }

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
    
    if (!t->root) {
        mfl_root_tree(t, 0, ntax);
        wasrooted = false;
    }
    
    mfl_clear_cpindex(t, numnodes);
    
    t->root->cpindex = ntax;
    t->root->bottom = true;
    
    j = ntax+1;
    
    for (i = 0; i < ntax; ++i) {        
        
        p = t->trnodes[i]->outedge;
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
                    p = p->outedge;

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
    
    return storedtr;
}

tree *mfl_decompress_tree(int *savedtr, int ntax, int numnodes)
{
    int i;
    tree *t;
    
    t = mfl_alloc_noring(ntax, numnodes);
    
    for (i = 0; i < numnodes; ++i) {
        if (i < ntax) {
            t->trnodes[i]->outedge = mfl_allocnode();
            t->trnodes[i]->outedge->outedge = t->trnodes[i];
        }
    }
    
    /* Stuff will go here that reconstructs the tree from the integer array representation */
    
    return t;
}

bool mfl_compare_alltrees(tree *newtopol, tree **savedtrees, int ntax, int numnodes, long int *start, long int *last)
{
    int i = 0;
    int *newtr;
    bool foundtr = false;
    
    newtr = mfl_compress_tree(newtopol, ntax, numnodes);
    
    /* This search should be replaced by a binary search or some kind of hash 
     * function in order to speed it up. Otherwise it will be a major time hog 
     * if/when somebody loads a noisy dataset. */
    
    for (i = *start; i < *last; ++i) {
        if (savedtrees[i]->compressedtr) {
            if (mfl_compare_trees(newtr, savedtrees[i]->compressedtr, numnodes)) {
                foundtr = true;
                break;
            }
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

void test_tree_compress(void)
{
    
    char testtr1[] = "(((2,3),1),(4,5));";
    char testtr2[] = "(((3,2),1),(5,4));";
    int ntax = 5;
    int numnodes = 2 * ntax - 1;
    
    tree *t = readNWK(testtr1, true);
    mfl_unroot(ntax, t);
    
    int *compressed = mfl_compress_tree(t, ntax, numnodes);
    
    int i;
    
    dbg_printf("\n");

    for (i = 0; i < numnodes; ++i) {
        dbg_printf("%i ", i);
    }
    dbg_printf("\n");
    
    for (i = 0; i < numnodes; ++i) {
        dbg_printf("%i ", compressed[i]);
    }
    dbg_printf("\n");
    
    tree *t2 = readNWK(testtr2, true);
    mfl_unroot(ntax, t2);
    int *compressed2 = mfl_compress_tree(t2, ntax, numnodes);
    
    if (mfl_compare_trees(compressed, compressed2, numnodes)) {
        dbg_printf("\nthey're the same\n");
    }
    else {
        dbg_printf("\nthey're diff'rnt\n");
    }

    
}

