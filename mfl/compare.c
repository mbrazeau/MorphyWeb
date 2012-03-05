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

int mfl_compare_ints(const void * a, const void * b)
{
    
    return ( **(taxbipart**)a - **(taxbipart**)b );
}

int mfl_compare_ints2(const void * a, const void * b)
{
    
    int i1=*(int*)a;
    int i2=**(taxbipart**)b; 
    return i1 - i2;
    //return ( **(taxbipart**)a - **(taxbipart**)b );
}

int mfl_count_fields(int ntax)
{
    /*counts the number of int64_t needed to cover all taxa in the tree */
    
    int i, numfields;
    bool quit = false;
    
    if (ntax > 64) {
        // Find the number of 64-bit fields needed to cover all taxa.
        for (i = 2; !quit; ++i) {
            if (64 * i >= ntax) {
                quit = true;
                numfields = i;
            }
        }
    }
    else {
        numfields = 1;
    }
    
    return numfields;
}

bool mfl_comp_bipartition(taxbipart *bp1, taxbipart *bp2, int numfields)
{
    /* Compares one bipartition array against another*/
    
    int i = 0;
    bool same = true;
    do {
        if ((bp1[i] != bp2[i])) {
            same = false;
        }
    } while (i < numfields);
    return same;
}

bool mfl_compare_trees(taxbipart **t1, taxbipart **t2, int ntax, int numfields)
{
    /* Performs a binary search for each bipartition of t1 in t2.
     *
     **** This should be supplemented by a function that compares number
     **** of bipartitions first and returns false if they are different */
    
    int i, j;
    bool eqtrees = false;
    
    /*for (i = 0; i < ntax - 1; ++i) {
        eqtrees = false;
        //if (!(bsearch(&(*t1[i]), t2, ntax - 1, sizeof(taxbipart*), mfl_compare_ints2))) {
          //  eqtrees = false;
        /}
        for (j = 0; j < ntax - 1; ++j) {
            if (t1[i] == t2[j]) {
                eqtrees = true;
                continue;
            }
        }
        if (!eqtrees) {
            break;
        }
    }*/
    for (i = 0; i < ntax - 1; ++i) {
        eqtrees = false;
        for (j = 0; j < ntax - 1; ++j) {
            if (*t1[i] == *t2[j]) {
                eqtrees = true;
                break;
            }
        }
        if (!eqtrees) {
            break;
        }
    }
    return eqtrees;
}

void mfl_set_bipartition(node *n, node *d)
{
    if (d->tip) {
        if (d->index > 64) {    /* * * * Check this comparison * * * */
            //something
        }
        else {
            n->tipsabove[0] = n->tipsabove[0] | (1 << d->index);
        }
        return;
    }
    else {
        n->tipsabove[0] = n->tipsabove[0] | d->tipsabove[0];
    }

    /*while () {
        //something
    }*/
    
    return;
}

void mfl_set_tipsabove(node *n, int numfields, taxbipart **hashtab, int *bpcounter)
{
    /* Traverses an n-ary tree, allocating memory for an attendant hashcode for 
     * each node and calls mfl_set_bipartition to set the hashcode for the 
     * above that node. */
    node *p;
    
    if (n->tip) {
        return;
    }
    
    if (!n->tipsabove) {
        n->tipsabove = (taxbipart*)malloc(numfields * sizeof(taxbipart));
    }
    memset(n->tipsabove, 0, numfields * sizeof(taxbipart));

    p = n->next;
    while (p != n) {
        mfl_set_tipsabove(p->outedge, numfields, hashtab, bpcounter);
        p = p->next;
    }
    
    p = n->next;
    while (p != n) {
        mfl_set_bipartition(n, p->outedge);
        p = p->next;
    }
    
    memcpy(hashtab[*bpcounter], n->tipsabove, sizeof(taxbipart));
    *bpcounter = *bpcounter + 1;
}

void mfl_free_hashtab(taxbipart **hashtab, int numbiparts)
{
    int i;
    
    for (i = 0; i < numbiparts; ++i) {
        free(hashtab[i]);
    }
    
    free(hashtab);
}

taxbipart **mfl_tree_biparts(tree *t,int ntax, int numnodes)
{
    /* Creates a bipartition table describing the tree as a series of taxon 
     * bipartitions. These bipartitions are then sorted in ascending order
     * */
    int i;
    int numbiparts = ntax - 1;
    int bpcount = 0;
    int *bpcounter = &bpcount;
    
    int numfields = mfl_count_fields(ntax);
    
    taxbipart **hashtab;
    
    hashtab = (taxbipart**) malloc(numbiparts * sizeof(taxbipart*));
    for (i = 0; i < numbiparts; ++i) {
        hashtab[i] = (taxbipart*) malloc(numfields * sizeof(taxbipart));
    }
    
    if (!t->root) {
        mfl_temproot(t, 0, ntax);
        mfl_set_tipsabove(t->root, numfields, hashtab, bpcounter);
        mfl_undo_temproot(ntax, t);
    }
    else {
        mfl_set_tipsabove(t->root, numfields, hashtab, bpcounter);
    }

    //qsort(hashtab, ntax - 1, sizeof(taxbipart*), mfl_compare_ints);
    
    return hashtab;
}

bool mfl_compare_alltrees(tree *newtopol, tree **savedtrees, int ntax, int numnodes, long int *current)
{
    int i = 0;
    int numfields = mfl_count_fields(ntax);
    taxbipart **temphashtab;
    bool foundtr = false;
    double timein;
    double timeout;
    static double totaltime = 0;
    static double increm = 0;
    
    timein = (double)(clock() / (double)CLOCKS_PER_SEC);
    
    temphashtab = mfl_tree_biparts(newtopol, ntax, numnodes);
    
    for (i = 0; i < *current; ++i ) {
        if (mfl_compare_trees(temphashtab, savedtrees[i]->bipartitions, ntax, numfields)) {
            if (newtopol != savedtrees[i]) {
                foundtr = true;
                break;
            }
        }
    }
    
    if (!foundtr) {
        newtopol->hashtabholder = temphashtab;
    }
    else {
        mfl_free_hashtab(temphashtab, ntax - 1);
    }
    
    timeout = (double)(clock() / (double)CLOCKS_PER_SEC);
    
    totaltime = totaltime + (timeout - timein);
    increm = increm + (timeout - timein);
    
    if (increm > 0.025) {
        printf("time in tree comparison: %g\n", totaltime);
        increm = 0;
    }
    
    //printf("time in comparison: %g\n", totaltime);
    
    return foundtr;
}

void test_tree_comparison(void)
{
    
    
    char tree1[] = "(1,((((25,(16,(6,2))),(((24,14),17),15)),7),((21,(10,8)),((9,(27,(28,(5,(26,((23,(18,13)),3)))))),(20,((((4,22),19),12),11))))));"; /*"(1,((2,5),(3,4)));";*/
    char tree2[] = "(1,((((25,(16,(6,2))),(((24,14),17),15)),7),((21,(10,8)),((20,((((4,22),19),12),11)),(9,(27,(28,(5,(26,((23,(18,3)),13))))))))));"; /*"(1,((4,5),(3,2)));"; */
    /*char tree3[] = "(1,((3,4),(2,5)));";
    char tree4[] = "(1,(2,(3,(4,5))));";
    char tree5[] = "((2,((5,4),3)),1);";*/
    
    /* tree3 == tree1 */
    
    int ntax = 28;
    int numnodes = 2 * ntax - 1;
    int numfields = mfl_count_fields(ntax);
    
    tree *t1 = readNWK(tree1, 0);
    printNewick(t1->trnodes[0]);
    printf("\n");
    tree *t2 = readNWK(tree2, 0);
    printNewick(t2->trnodes[0]);
    printf("\n");
    /*tree *t3 = readNWK(tree3, 0);
    printNewick(t3->trnodes[0]);
    printf("\n");
    tree *t4 = readNWK(tree4, 0);
    printNewick(t4->trnodes[0]);
    printf("\n");
    tree *t5 = readNWK(tree5, 0);
    printNewick(t5->trnodes[0]);
    printf("\n");*/
    
    // Get bipartitions for t1
    t1->bipartitions = mfl_tree_biparts(t1, ntax, numnodes);
    print_hashtab(t1->bipartitions, ntax);
    
    // Get bipartitions for t2
    t2->bipartitions = mfl_tree_biparts(t2, ntax, numnodes);
    print_hashtab(t2->bipartitions, ntax);
    
    // Get bipartitions for t3
    /*t3->bipartitions = mfl_tree_biparts(t3, ntax, numnodes);
    
    // Get bipartitions for t4
    t4->bipartitions = mfl_tree_biparts(t4, ntax, numnodes);
    
    // Get bipartitions for t5
    t5->bipartitions = mfl_tree_biparts(t5, ntax, numnodes);*/
    
    // Compare all trees
    if (mfl_compare_trees(t1->bipartitions, t2->bipartitions, ntax, numfields))
    {
        printf("trees 1 and 2 are equal\n");
    }
    else {
        printf("trees 1 and 2 are different\n");
    }

    /*if (mfl_compare_trees(t1->bipartitions, t3->bipartitions, ntax, numfields))
    {
        printf("trees 1 and 3 are equal\n");
    }
    else {
        printf("trees 1 and 3 are different\n");
    }
    
    if (mfl_compare_trees(t2->bipartitions, t3->bipartitions, ntax, numfields))
    {
        printf("trees 3 and 2 are equal\n");
    }
    else {
        printf("trees 3 and 2 are different\n");
    }
    if (mfl_compare_trees(t4->bipartitions, t5->bipartitions, ntax, numfields))
    {
        printf("trees 4 and 5 are equal\n");
    }
    else {
        printf("trees 4 and 5 are different\n");
    }
    
    if (mfl_compare_trees(t1->bipartitions, t5->bipartitions, ntax, numfields))
    {
        printf("trees 1 and 5 are equal\n");
    }
    else {
        printf("trees 1 and 5 are different\n");
    }*/
    
    /*mfl_free_hashtab(t1->bipartitions, ntax - 1);
    mfl_free_hashtab(t2->bipartitions, ntax - 1);
    mfl_free_hashtab(t3->bipartitions, ntax - 1);
    mfl_free_hashtab(t4->bipartitions, ntax - 1);
    mfl_free_hashtab(t5->bipartitions, ntax - 1);*/
}