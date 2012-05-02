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
    /* Performs a binary search for each hashcode of t1 in t2.
     *
     **** This should be supplemented by a function that compares number
     **** of bipartitions first and returns false if they are different */
    
    int i;
    bool eqtrees = false;
    
    for (i = 0; i < ntax - 3; ++i) {
        eqtrees = false;
        if (bsearch(&(*t1[i]), t2, ntax - 1, sizeof(taxbipart*), mfl_compare_ints2)) {
            eqtrees = true;
            continue;
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
    /* Creates a hashtable describing the tree as a series of taxon 
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
    
    qsort(hashtab, ntax - 1, sizeof(taxbipart*), mfl_compare_ints);
    
    return hashtab;
}

/*bool mfl_compare_alltrees(tree *newtopol, tree **savedtrees, int ntax, int numnodes, long int *current)
{
    int i = 0;
    int numfields = mfl_count_fields(ntax);
    taxbipart **temphashtab;
    bool foundtr = false;
    double timein;
    double timeout;
    static double totaltime = 0;
    static double increm = 0;
    
    //dbg_printf("time in: %g\n", (double)(clock() / (double)CLOCKS_PER_SEC));
    
    timein = (double)(clock() / (double)CLOCKS_PER_SEC);
    
    temphashtab = mfl_tree_biparts(newtopol, ntax, numnodes);
    
    for (i = *current; i--; ) {
        if (mfl_compare_trees(temphashtab, savedtrees[i]->bipartitions, ntax, numfields)) {
            if (newtopol != savedtrees[i]) {
                foundtr = true;
                break;
            }
        }
        //mfl_free_hashtab(tph2, ntax - 1);
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
        dbg_printf("time in tree comparison: %g\n", totaltime);
        increm = 0;
    }
    
    //dbg_printf("time in comparison: %g\n", totaltime);
    
    return foundtr;
}*/

/*bool mfl_compare_alltrees(tree *newtopol, tree **savedtrees, int ntax, int numnodes, long int *start, long int *current)
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
    
    for (i = *start; i < *current; ++i ) {
        if (!mfl_compare_trees(temphashtab, savedtrees[i]->bipartitions, ntax, numfields)) {
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
        dbg_printf("time in tree comparison: %g\n", totaltime);
        increm = 0;
    }
    
    //dbg_printf("time in comparison: %g\n", totaltime);
    
    return foundtr;
}*/

bool mfl_compare_trees(int *t1, int *t2, int numnodes)
{
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
        if (mfl_compare_trees(newtr, savedtrees[i]->compressedtr, numnodes)) {
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
    /*printNewick(t1->trnodes[0]);*/
    dbg_printf("\n");
    tree *t2 = readNWK(tree2, 0);
    /*printNewick(t2->trnodes[0]);*/
    dbg_printf("\n");
    /*tree *t3 = readNWK(tree3, 0);
    printNewick(t3->trnodes[0]);
    dbg_printf("\n");
    tree *t4 = readNWK(tree4, 0);
    printNewick(t4->trnodes[0]);
    dbg_printf("\n");
    tree *t5 = readNWK(tree5, 0);
    printNewick(t5->trnodes[0]);
    dbg_printf("\n");*/
    
    // Get bipartitions for t1
    t1->bipartitions = mfl_tree_biparts(t1, ntax, numnodes);
    /*print_hashtab(t1->bipartitions, ntax);*/
    
    // Get bipartitions for t2
    t2->bipartitions = mfl_tree_biparts(t2, ntax, numnodes);
    /*print_hashtab(t2->bipartitions, ntax);*/
    
    // Get bipartitions for t3
    /*t3->bipartitions = mfl_tree_biparts(t3, ntax, numnodes);
    
    // Get bipartitions for t4
    t4->bipartitions = mfl_tree_biparts(t4, ntax, numnodes);
    
    // Get bipartitions for t5
    t5->bipartitions = mfl_tree_biparts(t5, ntax, numnodes);*/
    
    // Compare all trees
    if (mfl_compare_trees(t1->bipartitions, t2->bipartitions, ntax, numfields))
    {
        dbg_printf("trees 1 and 2 are equal\n");
    }
    else {
        dbg_printf("trees 1 and 2 are different\n");
    }

    /*if (mfl_compare_trees(t1->bipartitions, t3->bipartitions, ntax, numfields))
    {
        dbg_printf("trees 1 and 3 are equal\n");
    }
    else {
        dbg_printf("trees 1 and 3 are different\n");
    }
    
    if (mfl_compare_trees(t2->bipartitions, t3->bipartitions, ntax, numfields))
    {
        dbg_printf("trees 3 and 2 are equal\n");
    }
    else {
        dbg_printf("trees 3 and 2 are different\n");
    }
    if (mfl_compare_trees(t4->bipartitions, t5->bipartitions, ntax, numfields))
    {
        dbg_printf("trees 4 and 5 are equal\n");
    }
    else {
        dbg_printf("trees 4 and 5 are different\n");
    }
    
    if (mfl_compare_trees(t1->bipartitions, t5->bipartitions, ntax, numfields))
    {
        dbg_printf("trees 1 and 5 are equal\n");
    }
    else {
        dbg_printf("trees 1 and 5 are different\n");
    }*/
    
    /*mfl_free_hashtab(t1->bipartitions, ntax - 1);
    mfl_free_hashtab(t2->bipartitions, ntax - 1);
    mfl_free_hashtab(t3->bipartitions, ntax - 1);
    mfl_free_hashtab(t4->bipartitions, ntax - 1);
    mfl_free_hashtab(t5->bipartitions, ntax - 1);*/
}
