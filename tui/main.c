/*
 *  morphy.h
 *  Morphy
 *
 *  Created by Martin Brazeau on 11-04-26.
 *  Copyright 2011. All rights reserved.
 *  Morphy is provided as is with no warranty of any kind.
 *
 */


#include "morphy.h"

/*temporary defs for testing*/
#define MORPHY_NUM_ITERATIONS 15
//#define MAXSTATES 5
/**/

void call_index(node *n)
{
    node *p;
    
    if (n->next) {
        printf("index: %i\n", n->index);
        p = n->next;
        while (p != n) 
        {
            printf("index: %i\n", p->index);
            p = p->next;
        }
    }   
}

void dump_nodearray(nodearray nds, int ntax, int numnodes)
{
    int i;
    node *p; 
    
    for (i = 0; i < numnodes; ++i) {
        mfl_set_ring_to_n(nds[i]);
        printf("Index: %i, tip: %i, order: %i, address: %p, outedge %p", nds[i]->index, nds[i]->tip, nds[i]->order, nds[i], nds[i]->outedge);
        if (i >= ntax && nds[i]->next) {
            printf(", next: %p\n", nds[i]->next);
            p = nds[i]->next;
            while (p != nds[i]) {
                printf("In ring: index: %i, tip: %i, order: %i, address: %p, outedge: %p, next: %p\n", p->index, p->tip, p->order, p, p->outedge, p->next);
                p = p->next;
            }
        }
        else 
        {
            printf("\n");
        }

    }
}

void dump_connections(nodearray nds, int ntax, int numnodes)
{
    int i, nbi;
    node *p; 
    
    
    for (i = 0; i < numnodes; ++i) {
        //mfl_set_ring_to_n(nds[i]);
        if (nds[i]->outedge) {
            nbi = nds[i]->outedge->index;
        }
        else {
            nbi = 999;
        }

        printf("Index: %i, tip: %i, OE index: %i\n", nds[i]->index, nds[i]->tip, nbi);
        if (i >= ntax && nds[i]->next) {
            p = nds[i]->next;
            while (p != nds[i]) {
                if (p->outedge) {
                    nbi = p->outedge->index;
                }
                else {
                    nbi = 999;
                }
                printf("In ring: index: %i, tip: %i, OE index %i\n", p->index, p->tip, nbi);
                p = p->next;
            }
        }
    }
}


void dump_tree(tree *t, int ntax, int numnodes)
{
    printf("Tree %i:\n", t->index);
    printf("Root: %p\nLength: %i\n", t->root, t->length);
    dump_nodearray(t->trnodes, ntax, numnodes);
    printf("\n");
}

void dump_tree_connections(tree *t, int ntax, int numnodes)
{
    printf("Tree %i:\n", t->index);
    printf("Root: %p\nLength: %i\n", t->root, t->length);
    dump_connections(t->trnodes, ntax, numnodes);
    printf("\n");
}


void close_all_rings(nodearray nds, int ntax, int numnodes)
{
    int i;
  
    for (i = ntax; i < numnodes; ++i) {
        mfl_close_ring(nds[i]);
    }
}

int numberOfNodes(int ntax)
{
    int numnodes;
    
    return numnodes = 2 * ntax - 1;
}

void init_taxarray(int *taxarray, int ntax)
{
    int i;
    
    for (i = 0; i < ntax; ++i) {
        taxarray[i] = i + 1;
    }
}

void joinNodes(node *n, node *p)
{
    if (n->outedge) {
        n->outedge->outedge = NULL;
        //printf("terminating existing connection in n\n");
    }
    if (p->outedge) {
        p->outedge->outedge = NULL;
        //printf("terminating existing connection in p\n");
    }
    
    n->outedge = p;
    p->outedge = n;
}

struct node * allocnode(void)
{
    node *newNode;
    newNode = (node *)malloc(sizeof(node));
    if (newNode == NULL)
    {
        printf("Error: failed to allocate new node.\n");
    }
    else
    {
        memset(newNode, 0, sizeof(node));
    }
    return newNode;
}

struct tree *alloctree(int ntax, int numnodes)
{
    int i; //Loop counters
    tree *newtree;
    
    newtree = (tree*)malloc(sizeof(tree));
    memset(newtree, 0, sizeof(tree));
    
    if (newtree == NULL)
    {
        printf("Error: failed to allocate new tree.\n");
        return (struct tree*) 0;
    }
    
    newtree->trnodes = (node **)malloc( (numnodes) * sizeof(node*));
    
    for (i = 0; i < numnodes; ++i)
    {
        newtree->trnodes[i] = allocnode();
    }
    
    for (i = 0; i < numnodes; ++i)
    {
        /* First half of trnodes are initialized to this */
        if (i < ntax)
        {
            newtree->trnodes[i]->tip = i + 1;
            newtree->trnodes[i]->index = i;
            newtree->trnodes[i]->next = NULL;
            newtree->trnodes[i]->vweight = 1;
        }
        else /* second half are allocated as a newring */
        {
            newtree->trnodes[i]->index = i;
            newring(newtree->trnodes[i], ntax);
        }
        
        newtree->trnodes[i]->outedge = NULL;
    }
    
    newtree->root = NULL;
    
    return (newtree);
}

/*
 * freetree - deletes an entire tree, by doing the mirror image
 * of alloctree.
 */
void freetree(tree *newtree, int numnodes)
{
    int i;
    
    /* free all of the trnodes */
    for (i = 0; i < numnodes; ++i)
    {
        if (newtree->trnodes[i]->next)
        {
            mfl_close_ring(newtree->trnodes[i]);
            deletering(newtree->trnodes[i]);
        }
        free(newtree->trnodes[i]);
    }
    /* free the trnode list */
    free(newtree->trnodes);
    /* free the tree */
    free(newtree);
    
}

struct tree *alloc_noring(int ntax, int numnodes)
{
    int i;
    tree *newtree;
    
    newtree = (tree*)malloc(sizeof(tree));
    
    if (newtree == NULL)
    {
        printf("Error: failed to allocate new tree.\n");
        return (struct tree*) 0;
    }
    
    newtree->trnodes = (node **)malloc( (numnodes) * sizeof(node*));
    
    for (i = 0; i < numnodes; ++i)
    {
        newtree->trnodes[i] = allocnode();
        if (i < ntax) 
        {
            newtree->trnodes[i]->tip = i + 1;
            newtree->trnodes[i]->vweight = 1;
        }
        newtree->trnodes[i]->index = i;
    }
    
    newtree->root = NULL;
    return (newtree);
    
}

void printNewick(node *n)
{   
    /* Prints the tree in Newick format (i.e. using brackets plus commas to 
     * separate equally ranked objects). In an unrooted tree, function will 
     * be called on a terminal node that has its start variable set to 1. */
    
    node *p;
    
    if (n->start) {
        printf("(%i,", n->tip);
        printNewick(n->outedge);
        printf(")");
        return;
    }   
    
    if (n->tip && n->outedge->next->outedge) {
        if (n->outedge->next->outedge->tip && 
            !n->outedge->next->outedge->start) {
            printf("%i", n->tip);
            return;
        }
    }
    
    if (n->tip && !n->start) {
        printf("%i", n->tip);
        return;
    }
    
    printf("(");
    
    p = n->next;
    while (p != n) {
        printNewick(p->outedge);
        p = p->next;
        if (p != n) {
            printf(",");
        }
    }   
    printf(")");
}

void newring(node *r1, int ntax)
{
    /* Generates a new internal node composed of a ring of node structures
     * joined in a unidirectional manner by their next pointers. Used by any
     * function that might dynamically add branches at run-time. */
    
    node *r2, *r3;
    
    r2 = allocnode();
    r3 = allocnode();
    
    r1->next = r2;
    r2->next = r3;
    r3->next = r1;
    
    r1->tip = r2->tip = r3->tip = 0;
    r1->start = r2->start = r3->start = false;
    r1->skip = r2->skip = r3->skip = false;
    
    r1->apomorphies = r2->apomorphies = r3->apomorphies;
    r2->index = r1->index;
    r3->index = r1->index;
}

void newring_to_order(node *r1, int order, int ntax)
{
    int i = 1;
    node *p;
    
    if (order == 0 || order == 1) {
        printf("Error in newring_to_order: no order specified or order too small\n");
        return;
    }
    
    p = r1;
    do {
        p->next = allocnode();
        p = p->next;
        p->tip = r1->tip;
        p->start = r1->start;
        p->index = r1->index;
        p->apomorphies = r1->apomorphies;
        ++i;
    } while (i < order);
    p->next = r1;
    r1->order = i;
    p = r1;
    do {
        p->order = r1->order;
        p = p->next;
    } while (p != r1);
    p->order = r1->order;
}

/*
 * deletering - deletes a node (ring) by doing the mirror image of newring.
 */
void deletering(node *r1)
{
    if (r1->tipsabove) {
        free(r1->tipsabove);
    }
    if (r1->apomorphies) {
        free(r1->apomorphies);
    }
    
    node *p, *q;
    
    p = r1->next;
    while (p != r1) {
        q = p->next;
        if (p->tipsabove) {
            free(p->tipsabove);
        }
        free(p);
        p = q;
    }
}

void reIndex(node *n, int index_val)
{
    node *p;
    
    n->index = index_val;
    
    if (n->next) {
        p = n->next;
        while (p != n) {
            p->index = n->index;
            p = p->next;
        }
    }
}

struct tree * copytree(tree *origtr, int ntax, int numnodes)
{
    int i, tmpltorder, begin;
    tree *treecp; // Pointer to the tree copy
    node *p, *q, *r, *s;
    //bool isRooted = false;
    
    treecp = alloc_noring(ntax, numnodes);
    
    if (origtr->root)
    {
        mfl_set_order(origtr->trnodes[ntax]);
        tmpltorder = origtr->trnodes[ntax]->order;
        newring_to_order(treecp->trnodes[ntax], tmpltorder + 1, ntax);
        treecp->trnodes[ntax]->order = treecp->trnodes[ntax]->order - 1;
        treecp->root = treecp->trnodes[ntax];
        //isRooted = true;
        begin = ntax;
    }
    else 
    {
        treecp->trnodes[0]->start = true;
        begin = ntax + 1;
        newring(treecp->trnodes[ntax], ntax);
    }
    
    
    for (i = ntax + 1; i < numnodes; ++i) 
    {
        if (origtr->trnodes[i]->next) 
        {
            mfl_set_order(origtr->trnodes[i]);
            tmpltorder = origtr->trnodes[i]->order;
            newring_to_order(treecp->trnodes[i], tmpltorder, ntax);
        }
    }
    
    // Add the tips first
    for (i = begin; i < numnodes; ++i) 
    {
        p = origtr->trnodes[i];
        q = treecp->trnodes[i];
        q->initialized = p->initialized;
        q->skip = p->skip;
        if (p->next) 
        {
            if (p->order != q->order) 
            {
                printf("Error: template and copy have non-matching branch orders\n");
            }
            do 
            {
                if (p->outedge && p->outedge->tip) 
                {
                    joinNodes(q, treecp->trnodes[p->outedge->index]);
                }
                if (p->outedge && !p->outedge->tip && !q->outedge) 
                {
                    r = origtr->trnodes[p->outedge->index];
                    s = treecp->trnodes[p->outedge->index];
                    while (r->outedge->index != p->index) 
                    {
                        r = r->next;
                        s = s->next;
                    }
                    joinNodes(q, s);
                }
                p = p->next;
                q = q->next;
                q->initialized = p->initialized;
                q->skip = p->skip;
            } while (p != origtr->trnodes[i] && q != treecp->trnodes[i]);
        }
    }
    //dump_tree_connections(origtr, ntax, numnodes);
    //dump_tree_connections(treecp, ntax, numnodes);
    treecp->bipartitions = mfl_tree_biparts(treecp, ntax, numnodes);
    return treecp;
}

void mfl_point_bottom(node *n, node **nodes, int ntax, int *iteration)
{
    /* Re-sets the pointers to the internal nodes in the tree's
     * node array to point to the 'bottom' (rootward) node in the ring*/
    
    node *p;
    
    if (n->tip) {
        return;
    }
    
    if (n->outedge) {
        nodes[*iteration] = n;
        n->index = *iteration;
        mfl_set_ring_to_n(n);
        *iteration = *iteration + 1;
    }
    
    p = n->next;
    while (p != n) {
        mfl_point_bottom(p->outedge, nodes, ntax, iteration);
        p = p->next;        
    }
    
}

void mfl_root_tree(tree *trtoroot, int nRoot, int ntax)
{    
    /*Roots the tree between a terminal (leaf) and an internal node*/
 
    int counter = ntax + 1;
    int *count_ptr = &counter;
    
    node *nodeptr, *r2, *r3;
    
    if (!trtoroot->trnodes[ntax]->next) {
        newring(trtoroot->trnodes[ntax]->next, ntax);
    }
    
    nodeptr = trtoroot->trnodes[nRoot]->outedge;
    r2 = trtoroot->trnodes[ntax]->next;
    r3 = trtoroot->trnodes[ntax]->next->next;
    
    joinNodes(r2, trtoroot->trnodes[nRoot]);
    joinNodes(r3, nodeptr);
    
    trtoroot->root = trtoroot->trnodes[ntax];
    trtoroot->trnodes[ntax]->outedge = NULL;
        
    mfl_point_bottom(trtoroot->root, trtoroot->trnodes, ntax, count_ptr);
    trtoroot->trnodes[0]->start = false;
}

void collapseBiNode(node *n)
{
    /*collapses a binary node*/
    node *n2, *n3;
    node *an1, *an2;
    
    n2 = n->next;
    n3 = n2->next;
    
    an1 = n->outedge;
    an2 = an1->next;
    
    joinNodes(an1, n2->outedge);
    
    an1->next = n3;
    n3->next = an2;
    
    n2->index = an1->index;
    n2->index = an1->index;
    
    n->next = NULL;
    
}

void unroot(int ntax, tree *rootedtree)
{
    int lnumnodes = 2*ntax-1;
    node *proot, *leftdesc, *rightdesc, *subnode, *n, *m, *p, *q, *ftip;
    
    proot = rootedtree->root;
    
    mfl_set_order(proot);
    
    if (proot->order > 2) 
    {
        printf("derooting multifurcating node\n");
        
        subnode = mfl_seek_internal(ntax, lnumnodes, rootedtree->trnodes);
        
        /*grab the second branch in the ring*/
        
        p = proot->next->next;
        
        /*find the last node in the ring before the bottom node*/
        
        q = p;
        while (q->next != proot) {
            q = q->next;
        }
        
        q->next = NULL;
        
        if (subnode->next)
        {
            printf("Doing the trickier one\n");
            n = subnode->next;
            m = subnode;
            while (m->next != subnode) {
                m = m->next;
            }
            m->next = NULL;
            
            ftip = proot->next->outedge;
            joinNodes(ftip, subnode);
            subnode->next = p;
            q->next = subnode;
            proot->next = n;
            m->next = proot;
            
        }
        else 
        {
            ftip = proot->next->outedge;
            joinNodes(ftip, subnode);
            proot->next = NULL;
            subnode->next = p;
            q->next = subnode;
        }
    }
    else 
    {
        leftdesc = proot->next->outedge;
        rightdesc = proot->next->next->outedge;
        joinNodes(leftdesc, rightdesc);
    }
    
    rootedtree->root = NULL;
    rootedtree->trnodes[0]->start = true;

}

void pauseit(void)
{
    int c;
    
    do {
        printf("Enter c to continue: \n");
        c = getchar();
    } while (c != 'c');
    c = getchar();
}

void rand_tree (int ntax, int numnodes) 
{
    int i; //Loop counter
    
    /*generate a random tree*/
    
    tree **randtrees;
    randtrees = (tree **) malloc(MORPHY_NUM_ITERATIONS * sizeof(tree*));
    if (randtrees == NULL) {
        printf("Error in main(): failed malloc for randtrees\n");
        return;
    }
    
    for (i = 0; i < MORPHY_NUM_ITERATIONS; ++i) {
        randtrees[i] = randrooted(ntax, numnodes);
        printNewick(randtrees[i]->root);
        printf(";\n");
    }
    
    for (i = 0; i < MORPHY_NUM_ITERATIONS; ++i) {
        freetree(randtrees[i], numnodes);
    }
    
    free(randtrees);
    
}

void mfl_apply_tipdata(tree *currenttree, charstate *tipdata, int ntax, int nchar)
{   
    int i;
    
    for (i = 0; i < ntax; ++i) {
        currenttree->trnodes[i]->apomorphies = &tipdata[i * nchar];
    }
}

/*int mfl_pg_shortcut(node *src, node *t1, node *t2, node *subtr, tree *swapingon, int nchar)
{
    int i;
    int length;
    
    for (i = 0; i < nchar; ++i) {
        subtr->apomorphies[i] = t1->apomorphies[i] | t2->apomorphies[i];
        
        if (src->a) {
            //stuff
        }
        
    }
    
}*/

void mfl_countsteps(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *trlength, int *besttreelen)
{
    int i;
    charstate lft_chars, rt_chars;
    
    for (i = 0; i < nchar; ++i) {
        if (leftdesc->apomorphies[i] & rightdesc->apomorphies[i]) 
        {
            ancestor->apomorphies[i] = leftdesc->apomorphies[i] & rightdesc->apomorphies[i];
        }
        else
        {
            lft_chars = leftdesc->apomorphies[i];
            rt_chars = rightdesc->apomorphies[i];
            ancestor->apomorphies[i] = lft_chars | rt_chars;
            if ((ancestor->apomorphies[i] & IS_APPLIC) && ( ((ancestor->apomorphies[i] & IS_APPLIC) & lft_chars) && ((ancestor->apomorphies[i] & IS_APPLIC) & rt_chars) )) 
            {
                ancestor->apomorphies[i] = (ancestor->apomorphies[i] & IS_APPLIC);
                *trlength = *trlength + 1;
                if (*trlength > *besttreelen) 
                {
                    return;
                }
            }
        }
    }
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
    mfl_countsteps(n->next->outedge, n->next->next->outedge, n, nchar, trlength, besttreelen);
}

void mfl_combine_up(node *n, node *anc, int nchar)
{
    int i;
    charstate lft_chars, rt_chars;
    
    for (i = 0; i < nchar; ++i) {
        
        if ((n->apomorphies[i] & anc->apomorphies[i]) == anc->apomorphies[i]) 
        {
            n->apomorphies[i] = n->apomorphies[i] & anc->apomorphies[i];
        }
        else
        {
            lft_chars = n->next->outedge->apomorphies[i];
            rt_chars = n->next->next->outedge->apomorphies[i];
            
            if (lft_chars | rt_chars == n->apomorphies[i]) 
            {
                n->apomorphies[i] = n->apomorphies[i] | anc->apomorphies[i];
            }
            else 
            {
                n->apomorphies[i] = n->apomorphies[i] | (anc->apomorphies[i] & (lft_chars & rt_chars));
            }
        }
    }
}

void mfl_fitch_preorder(node *n, int *trlength, int nchar, int *besttreelen)
{
    node *p, *dl, *dr;
    
    if (n->tip) {
        return;
    }
    
    p = n->next;
    dl = p->outedge;
    dr = p->next->outedge;
    if (!dl->tip) {
        mfl_combine_up(dl, n, nchar);
    }
    if (!dr->tip) {
        mfl_combine_up(dr, n, nchar);
    }
    
    while (p != n) {
        mfl_fitch_preorder(p->outedge, trlength, nchar, besttreelen);
        p = p->next;
    }
    
}

int mfl_get_subtreelen(node *n, charstate *tipdata, int ntax, int nchar, int *besttreelen)
{
    int treelen = 0;
    int *treelen_p = &treelen;
    
    mfl_subtree_postorder(n, treelen_p, nchar, besttreelen);
    mfl_fitch_preorder(n, treelen_p, nchar, besttreelen);
    
    return *treelen_p;
}

int mfl_get_treelen(tree *testtree, charstate *tipdata, int ntax, int nchar, int *besttreelen)
{
    int treelen = 0;
    int *treelen_p = &treelen;
    
    mfl_fitch_postorder(testtree->root, treelen_p, nchar, besttreelen);
    mfl_fitch_preorder(testtree->root, treelen_p, nchar, besttreelen);
    
    return *treelen_p;
}

charstate * mfl_convert_tipdata(char *txtsrc, int ntax, int nchar)
{
    int i, j;
    
    charstate *tipdata = (charstate*) malloc(ntax * nchar * sizeof(charstate));
    
    for (i = 0, j = 0; txtsrc[i]; ++i) {
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
            tipdata[j] = -1;
        }
        else if (txtsrc[i] == '-') {
            tipdata[j] = 1;
        }
        else if (txtsrc[i] == '\n') {
            ++i;
        }
        else {
            ++i;
        }
        
        ++j;
    }
    
    return tipdata;
}

void mini_test_analysis(void)
{    
    /* The content of this function will model a simple analysis. Obviously,
     * of the content of this function will be moved to other functions or 
     * will be stored differently. For now, it's just one place for stuff that 
     * needs to be available to the search algorithms before they can find the 
     * shortest tree. The analysis should converge towards a tree like the one
     * in the Newick string below. */
    
    int i, j = 0;
    int ntax = 28, nchar = 95; // These values would be supplied by the datafile
    int treelimit = 10000; // Will always start its life as a default value and be 
                         // changed by the user (or automatically, as a 
                         // condition of user preference)
    int besttreelen = 0;
    int numnodes;
    
    numnodes = numberOfNodes(ntax);
    
    /* A completely balanced tree topology. 
     * This is simple enough to create a 'rigged' dataset for.*/
    //char aNewickTree[] = "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));";
    
    /* A search with this dataset rigged to favour the topology of aNewickTree */
    char usrTipdata[] = "00000000000000000000000---0000000-0000000000000000000-0-0000000000000000000--00000000000000000?\
11-10121001120200001001---0000100-1--10001000001101010320010000010100001000--30300000000000010?\
??????????????????????30--0010??0-????????1-?-??100?0-??1?--????0?00??10??1112??????100111-0-1?\
00??0?????????????????3111100?????????????011100??0???????11?1????????10?????2??????1?1-1??????\
000?0???????11??01111030--???1010-???0????1---0-10010-3?00--??010?000??0??111?1-?0?111011????1?\
11-101201011212?0001002---0001100-01010011000001101010320010000010101001000--10310000000000000?\
00?100(02)???01?1(23)0???10?0---00?00010002?001100000001000-31-010????-?30?000000--01-0???00?0000100?\
1011023111112131???1100---000000101--200010000101000103011?01-1-01001101010?-11-10??00?0?00010?\
?????????????????0????3000?10?010-0100?01100-0--100011300?10????0????0010110021-1???1001101000?\
10110231101111311001110---0000000-1---000100001010000-3011101-1-01311001010--11-001000000000100\
????????????????0?????3000?1?0??0-010?1?1100100?100?11??0?00????0??0?000111002??????10011010001\
????????????????0001??3101?100????0??????1001100100011??0?010100??????00??1002????1?10??????001\
????????????2????11?????--????????????????1-?-??1?????????--?????0?00?10?????2?????110011????10\
11-11?301?1111(23)?0010002---000000110121011100000111001011-010?0001?011000000--?021?0000000001001\
11-11?30012120000000001---000000110021011100000101000-11-000000010010001000--101110?00000000000\
0011012?1011?12??011002---00001010002?001100000?101010??1??0??001?10?00100??-???1??0???????0000\
1001??30??112???001001?---0000000-1--2011100000111000-2?-010000010211001000--0021?0000000000000\
?????????????????11???30--0010010-???????11-?-??10010-3010--0001?0000010?????2??????10011??0-10\
0???0?1??1?0?1??00????311111??010-0100?011001100100011??0?01??????0??00011???21-101?1001(01)010001\
00000?101110113000011030000100010-0100101100100010001?30001000000000010{01}1110021-101110?101-0001\
00100000101011101000110---0-000010011000010000101000103011100?0?00?0?101010--11-0?1?00000000000\
00000?1?1??0113?00011031111101010-01?0?0??011100100?0-3000010100?0?0??10?110(01)21-?011111-1010?11\
?00?????1?1?213?01?1????--0011????????????1-?-??10010-??0?--??0??0000??0??????????????????????0\
11-11?3???1?2???000000?---000000110111011100000111001011-010000010011000000--00211??00000001000\
00?00000?01021100000002---000010111--1001000000110101032001000001000?001000--10010??00?00000000\
??????????????????????30--????????????????1-?????0????????????0????0??????????????????????????0\
?????????????????1?1??3000?0?1????010????1001???100011??0??000??0??0?000??11021-1??1100111-0001\
?????????????????????????????????????????100?????????????0?000??????0?00?????21-???1100111---1?";
    
    /*real data*/ /*"0000000000?00000100000000000000000000000000000100000100010000000000000101000101000000000010201100-0?0100--0--0000000010000000001320-?0????\
2000000000?00000-00000000000000000000000000000000000100010000000000000000000101000010000000001100-0?0100--0--00000000000001000002002?01???\
??2?????11?00000?????????????????0?10?000???????211131?1????11?1????????????????????0???????????0-0?????0????1211????10????10????2????????\
0000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000101000-000000--0--00000000000000000000000??0???\
0000000000000000000000000000000000000000000000000000100010000000000000000000000000000000000000000-000100--0--00000000000000000010001001000\
1010030011200000100000001001000000000000000000101000201020210001110010000000000000000000000000000-1001----5--12110001121102000111002112111\
1010030011200000100000001001000000000000000000101000201020210001110010000000000000000000000000000-1001----5--12110001121102000011002112111\
1110040111200000110001002000000110001000010000100010100111210001210120100100001000100000000001010-010111--0--11110000201000000013012112111\
212111111012101122211211211211122111111201211111211111211111111122112110011111112021111111101011100112321110012111112201001121102112113000\
212111111012101122211211211211122111111201211111211111211111111122112110011111112021111111101011100112321110012111112201001121102112113000\
212111111122101122211211211211122101110111211211211111111111011122112010011111112020111101101001100112220010112111112201002121001102113000\
212112101012111122111211311311122111110222111111211131011111111122112021111101111130111111101111110112323111012111112310111110202102113000\
212112111012111122111211311311122111111222111211211111111111111122112021011101111030111111101001100112321021112111112300111100102212113000\
212112101011111122111211311311122111111222111211211131211111111122112021011101111130111121121111110112322132112111122300111111102122113000\
212112101012111122211211311311122111111222111211211111111111111122112021011101111030111111101111100112322031012101112310111100102212113000\
212112101011111122211211311311122111111222011211211131011111111122112021011101111130111121121111110112322243112111112310111120202202113000\
?0100?0?0??0?00??0??00???00?????0??0000?0100??0?11000000????????????????????????????00???000??????????????0??1011?000100000001003?????????";*/
    
    /*"00000000000000000000000---0000000-0000000000000000000-0-0000000000000000000--00000000000000000?\
     11-10121001120200001001---0000100-1--10001000001101010320010000010100001000--30300000000000010?\
     ??????????????????????30--0010??0-????????1-?-??100?0-??1?--????0?00??10??1112??????100111-0-1?\
     00??0?????????????????3111100?????????????011100??0???????11?1????????10?????2??????1?1-1??????\
     000?0???????11??01111030--???1010-???0????1---0-10010-3?00--??010?000??0??111?1-?0?111011????1?\
     11-101201011212?0001002---0001100-01010011000001101010320010000010101001000--10310000000000000?\
     00?100(02)???01?1(23)0???10?0---00?00010002?001100000001000-31-010????-?30?000000--01-0???00?0000100?\
     1011023111112131???1100---000000101--200010000101000103011?01-1-01001101010?-11-10??00?0?00010?\
     ?????????????????0????3000?10?010-0100?01100-0--100011300?10????0????0010110021-1???1001101000?\
     10110231101111311001110---0000000-1---000100001010000-3011101-1-01311001010--11-001000000000100\
     ????????????????0?????3000?1?0??0-010?1?1100100?100?11??0?00????0??0?000111002??????10011010001\
     ????????????????0001??3101?100????0??????1001100100011??0?010100??????00??1002????1?10??????001\
     ????????????2????11?????--????????????????1-?-??1?????????--?????0?00?10?????2?????110011????10\
     11-11?301?1111(23)?0010002---000000110121011100000111001011-010?0001?011000000--?021?0000000001001\
     11-11?30012120000000001---000000110021011100000101000-11-000000010010001000--101110?00000000000\
     0011012?1011?12??011002---00001010002?001100000?101010??1??0??001?10?00100??-???1??0???????0000\
     1001??30??112???001001?---0000000-1--2011100000111000-2?-010000010211001000--0021?0000000000000\
     ?????????????????11???30--0010010-???????11-?-??10010-3010--0001?0000010?????2??????10011??0-10\
     0???0?1??1?0?1??00????311111??010-0100?011001100100011??0?01??????0??00011???21-101?1001(01)010001\
     00000?101110113000011030000100010-0100101100100010001?30001000000000010{01}1110021-101110?101-0001\
     00100000101011101000110---0-000010011000010000101000103011100?0?00?0?101010--11-0?1?00000000000\
     00000?1?1??0113?00011031111101010-01?0?0??011100100?0-3000010100?0?0??10?110(01)21-?011111-1010?11\
     ?00?????1?1?213?01?1????--0011????????????1-?-??10010-??0?--??0??0000??0??????????????????????0\
     11-11?3???1?2???000000?---000000110111011100000111001011-010000010011000000--00211??00000001000\
     00?00000?01021100000002---000010111--1001000000110101032001000001000?001000--10010??00?00000000\
     ??????????????????????30--????????????????1-?????0????????????0????0??????????????????????????0\
     ?????????????????1?1??3000?0?1????010????1001???100011??0??000??0??0?000??11021-1??1100111-0001\
     ?????????????????????????????????????????100?????????????0?000??????0?00?????21-???1100111---1?";
     */    
    /*rigged data*/
    /*"11111111110000000000\
11111111110000000000\
00111111110000000000\
00001111110000000000\
00000011110000000000\
00000000110000000000\
00000000000000000011\
00000000000000001111\
00000000000000111111\
00000000000011111111\
00000000001111111111\
00000000001111111111";*/
        
        /* a noisier dataset */
        /*"11111111110000000011
00?1001???0000011011
00111111110000000000
00001111110000000000
00000011110000000000
00000000110000000000
00100000000000000000
00001100000000001110
00000001100000111110
001000001000111?1101
00000101001111111101
00000000001111111111";*/
        
        
    
    printf("User tip data:\n");
    for (i = 0; i < ntax; ++i) {
        for (j = 0; j < nchar; ++j) {
            printf("%c", usrTipdata[j + i * nchar]);
        }
        printf("\n");
    }
    printf("\n");
    
    charstate *morphyTipdata = mfl_convert_tipdata(usrTipdata, ntax, nchar);
    
    //int *currentchar = &j;
    
    /* Initialize the tree array that will store optimal trees */
    tree **savedtrees = (tree**) malloc(TREELIMIT * sizeof(tree*));
        
    // Start with a (random, in this case) starting tree. 
    // Different algorithms will be written for doing this (as there are better
    // ways to do it) but we'll use randunrooted for now.
    
    mfl_addseq_randasis(ntax, nchar, numnodes, morphyTipdata, 1, savedtrees);
    //Get a length for the starting tree.
    mfl_root_tree(savedtrees[0], 1, ntax);
    int *besttreelen_p = &besttreelen;

    *besttreelen_p = mfl_get_sttreelen(savedtrees[0], morphyTipdata, ntax, nchar, besttreelen_p);
    savedtrees[0]->length = *besttreelen_p;
    
    unroot(ntax, savedtrees[0]);
    
    // Now start making rearrangements to the starting tree and test them against
    // the random tree. If a shorter tree is found, discard the random tree and
    // keep only the shortest tree. Reset besttreelen_p to the length of the new
    // shortest tree. Continue swapping until a maximum number of iterations is
    // hit, or maxtrees is hit, or there are no more possible re-arrangments.

    printf("The starting tree:\n");
    printNewick(savedtrees[0]->trnodes[0]);
    printf("\n");
    printf("The length of the starting tree: %i steps\n\n", *besttreelen_p);
    
    // Choice of NNI or SPR heuristic search. Will add TBR.
    
    //mfl_nni_search(ntax, nchar, numnodes, morphyTipdata, savedtrees, besttreelen);
    
    mfl_spr_search(ntax, nchar, numnodes, morphyTipdata, savedtrees, besttreelen);
    
    //mfl_clear_treebuffer(savedtrees, , numnodes);   
}

int main(void)
{
    int ntax = 9;
    int numnodes;
    //bool isRooted = true;
    
    test_tree_comparison();
    
    numnodes = numberOfNodes(ntax);
    
    tree *anewtree;    
    tree *originaltree;
    tree *copiedtree;
    
    mini_test_analysis();
    
    //test_spr();
    
    /* This part is just for testing the collapseBiNode*/
    /*anewtree = randrooted(ntax, numnodes);
    printf("New tree: ");
    printNewick(anewtree->root);
    printf("\n");
    //dump_connections(anewtree->trnodes, ntax, numnodes);
    //dump_nodearray(anewtree->trnodes, ntax, numnodes);
    
    copiedtree = copytree(anewtree, ntax, numnodes);
    printf("A copied rooted tree: ");
    printNewick(copiedtree->root);
    printf("\n");
    
    mfl_bswap(anewtree->trnodes[1], anewtree->trnodes[2]);
    printf("Original tree w swapped br: ");
    printNewick(anewtree->root);
    printf("\n");
    
    copiedtree = copytree(anewtree, ntax, numnodes);
    printf("Copy with swapped branch: ");
    printNewick(copiedtree->root);
    printf("\n");
    
    freetree(copiedtree, numnodes);
    
    mfl_collapse(anewtree->trnodes[ntax + 2], anewtree->trnodes); // Magic number just for testing
    printf("With collapsed node: ");
    printNewick(anewtree->root);
    printf("\n");    
    //dump_connections(anewtree->trnodes, ntax, numnodes);
    
    copiedtree = copytree(anewtree, ntax, numnodes);
    //dump_nodearray(copiedtree->trnodes, ntax, numnodes);
    printf("Copying with collapsed node: ");
    printNewick(copiedtree->root);
    printf("\n");
    
    freetree(copiedtree, numnodes);
    
    mfl_arb_resolve(anewtree->trnodes[ntax + 1], anewtree->trnodes, ntax, numnodes); // Magic number just for testing
    //dump_connections(anewtree->trnodes, ntax, numnodes);
    printf("With resolved node: ");
    printNewick(anewtree->root);
    printf("\n");
    
    copiedtree = copytree(anewtree, ntax, numnodes);
    //dump_nodearray(copiedtree->trnodes, ntax, numnodes);
    printf("Copying with resolved node: ");
    printNewick(copiedtree->root);
    printf("\n");
    
    freetree(anewtree, numnodes);
    freetree(copiedtree, numnodes);
    
    originaltree = randunrooted(ntax, numnodes);
    printf("unrooted test: ");
    printNewick(originaltree->trnodes[0]);
    printf("\n");
    
    copiedtree = copytree(originaltree, ntax, numnodes);
    printf("Copying of unrooted tree: ");
    printNewick(copiedtree->trnodes[0]);
    printf("\n");
    
    freetree(originaltree, numnodes);
    freetree(copiedtree, numnodes);*/
    
    return 0;
}
