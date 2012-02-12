/*
 *  compare.c
 *  Morphy
 *
 *  Created by Martin Brazeau and Chris Desjardins on 2/8/12.
 *  Copyright 2012. All rights reserved.
 *
 */

#include "morphy.h"

int mfl_count_fields(int ntax)
{
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
    int i, j;
    bool eqtrees = false;
    
    for (i = 0; i < ntax - 1; ++i) {
        eqtrees = false;
        for (j = 0; j < ntax - 1; ++j) {
            /*if (mfl_comp_bipartition(t1[i], t2[j], numfields)) {
                eqtrees = true;
            }*/
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
    
    mfl_temproot(t, 0, ntax);
    mfl_set_tipsabove(t->root, numfields, hashtab, bpcounter);
    mfl_undo_temproot(ntax, t);
    
    return hashtab;
}

void print_hashtab(taxbipart **hashtab, int numbiparts)
{
    int i;
    taxbipart tempht;
    
    for (i = 0; i < numbiparts; ++i) {
        tempht = *hashtab[i];
        /*do {
         if (tempht & 1) {
         printf("*");
         }
         else {
         printf(".");
         }
         tempht >> 1;
         } while (tempht != 0);*/
        printf("%i\n", tempht);
    }
}

bool mfl_compare_alltrees(tree *newtopol, tree **savedtrees, int ntax, int numnodes, long int *current)
{
    int i = 0;
    int numfields = mfl_count_fields(ntax);
    taxbipart **temphashtab;
    taxbipart **tph2;
    bool foundtr = false;
    
    temphashtab = mfl_tree_biparts(newtopol, ntax, numnodes);
    
    for (i = 0; i < *current; ++i) {
        //tph2 = mfl_tree_biparts(savedtrees[i], ntax, numnodes);
        if (mfl_compare_trees(temphashtab, savedtrees[i]->bipartitions, ntax, numfields)) {
            if (newtopol != savedtrees[i]) {
                foundtr = true;
                break;
            }
        }
        //mfl_free_hashtab(tph2, ntax - 1);
    }
    
    mfl_free_hashtab(temphashtab, ntax - 1);
    
    return foundtr;
}

void test_tree_comparison(void)
{
    
    
    char tree1[] = "(1,((((25,(16,(6,2))),(((24,14),17),15)),7),((21,(10,8)),((9,(27,(28,(5,(26,((23,(18,3)),13)))))),(20,((((4,22),19),12),11))))));"/*"(1,((2,5),(3,4)));"*/;
    char tree2[] = "(1,((((25,(16,(6,2))),(((24,14),17),15)),7),((21,(10,8)),((20,((((4,22),19),12),11)),(9,(27,(28,(5,(26,((23,(18,3)),13))))))))));"/*"(1,((4,5),(3,2)));"*/;
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
    print_hashtab(t1->bipartitions, ntax - 1);
    
    // Get bipartitions for t2
    t2->bipartitions = mfl_tree_biparts(t2, ntax, numnodes);
    print_hashtab(t2->bipartitions, ntax - 1);
    
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