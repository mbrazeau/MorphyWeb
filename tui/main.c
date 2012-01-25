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
        }
        else /* second half are allocated as a newring */
        {
            newtree->trnodes[i]->index = i;
            newring(newtree->trnodes[i]);
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



void mfl_fitch_postorder(node *n, int *trlength)
{
    node *p;
    charstate lft_chars, rt_chars;
    charstate ancstate;
     
    if (n->tip) {
        return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_fitch_postorder(p->outedge, trlength);
        p = p->next;
    }
    
    if (n->next->outedge->apomorphies & n->next->next->outedge->apomorphies) 
    {
        ancstate = n->next->outedge->apomorphies & n->next->next->outedge->apomorphies;
    }
    else
    {
        lft_chars = n->next->outedge->apomorphies;
        rt_chars = n->next->next->outedge->apomorphies;
        ancstate = lft_chars | rt_chars;
        if ((ancstate & (-1 ^ 1)) && ( ((ancstate & (-1 ^ 1)) & lft_chars) && ((ancstate & (-1 ^ 1)) & rt_chars) )) {
            ancstate = (ancstate & (-1 ^ 1));
            *trlength = *trlength + 1;
        }
    }
    n->apomorphies = ancstate;
}


void newring(node *r1)
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
    r1->dummy = r2->dummy = r3->dummy = false;
    
    r1->apomorphies = r2->apomorphies = r3->apomorphies = 0;
    r2->index = r1->index;
    r3->index = r1->index;
}

void newring_to_order(node *r1, int order)
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
    /* Note: apomorphies is only allocated once, r2 and r3 both point to the r1 apomorphies. */
    //free(r1->apomorphies);
    
    node *p, *q;
    
    p = r1->next;
    while (p != r1) {
        q = p->next;
        free(p);
        p = q;
    }
}

void mfl_apply_tipdata(tree *currenttree, charstate *tipdata, int ntax, int nchar, int currentchar)
{   
    int i;
    
    for (i = 0; i < ntax; ++i) {
        currenttree->trnodes[i]->apomorphies = tipdata[currentchar + i * nchar];
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
        newring_to_order(treecp->trnodes[ntax], tmpltorder + 1);
        treecp->trnodes[ntax]->order = treecp->trnodes[ntax]->order - 1;
        treecp->root = treecp->trnodes[ntax];
        //isRooted = true;
        begin = ntax;
    }
    else 
    {
        treecp->trnodes[0]->start = true;
        begin = ntax + 1;
        newring(treecp->trnodes[ntax]);
    }
    
    
    for (i = ntax + 1; i < numnodes; ++i) {
        if (origtr->trnodes[i]->next) {
            mfl_set_order(origtr->trnodes[i]);
            tmpltorder = origtr->trnodes[i]->order;
            newring_to_order(treecp->trnodes[i], tmpltorder);
        }
    }
    
    // Add the tips first
    for (i = begin; i < numnodes; ++i) {
        p = origtr->trnodes[i];
        q = treecp->trnodes[i];
        if (p->next) {
            if (p->order != q->order) {
                printf("Error: template and copy have non-matching branch orders\n");
            }
            do {
                if (p->outedge && p->outedge->tip) {
                    joinNodes(q, treecp->trnodes[p->outedge->index]);
                }
                if (p->outedge && !p->outedge->tip && !q->outedge) {
                    r = origtr->trnodes[p->outedge->index];
                    s = treecp->trnodes[p->outedge->index];
                    while (r->outedge->index != p->index) {
                        r = r->next;
                        s = s->next;
                    }
                    joinNodes(q, s);
                }
                p = p->next;
                q = q->next;
            } while (p != origtr->trnodes[i] && q != treecp->trnodes[i]);
        }
    }
    //dump_tree_connections(origtr, ntax, numnodes);
    //dump_tree_connections(treecp, ntax, numnodes);
    
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
        newring(trtoroot->trnodes[ntax]->next);
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

int mfl_get_treelen(tree *testtree, charstate *tipdata, int ntax, int nchar)
{
    int i;
    int treelen = 0;
    int *treelen_p = &treelen;
    int currenchar = 0;
    int *currentchar_p = &currenchar;
    
    for (i = 0, *currentchar_p = 0; i < nchar; ++i) {
        mfl_apply_tipdata(testtree, tipdata, ntax, nchar, *currentchar_p);
        *currentchar_p = *currentchar_p + 1;        
        mfl_fitch_postorder(testtree->root, treelen_p);
    }
    return *treelen_p;
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
    int ntax = 12, nchar = 20; // These values would be supplied by the datafile
    int treelimit = 1000; // Will always start its life as a default value and be 
                         // changed by the user (or automatically, as a 
                         // condition of user preference)
    int treelength1 = 0, treelength2 = 0, besttreelen = 0;
    int numnodes;
    bool isRooted = true;
    
    int *treelength_p = &treelength1;
    int *treelength_q = &treelength2;
    
    numnodes = numberOfNodes(ntax);
    
    tree *anewtree;
    tree *arandomtree;
    
    /* A completely balanced tree topology. 
     * This is simple enough to create a 'rigged' dataset for.*/
    char aNewickTree[] = "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));";
    
    anewtree = readNWK(aNewickTree, isRooted);
    arandomtree = randrooted(ntax, numnodes);
    printf("\nThis is the target tree: \n");
    printNewick(anewtree->root);
    printf("\n\n");
    
    /* A search with this dataset rigged to favour the topology of aNewickTree */
    char usrTipdata[] = {  
        
        /* a 'clean' dataset */
        '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
        '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
        '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
        '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
        '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
        '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1',
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1',
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1',
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1',
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',

        /* a 'noisier' dataset */
        /*'1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1',
        '0', '0', '?', '1', '0', '0', '1', '?', '?', '?', '0', '0', '0', '0', '0', '1', '1', '0', '1', '1',
        '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
        '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
        '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
        '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
        '0', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
        '0', '0', '0', '0', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '0',
        '0', '0', '0', '0', '0', '0', '0', '1', '1', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '0',
        '0', '0', '1', '0', '0', '0', '0', '0', '1', '0', '0', '0', '1', '1', '1', '?', '1', '1', '0', '1',
        '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1',
        '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',*/
    };
    
    printf("User tip data:\n");
    for (i = 0; i < ntax; ++i) {
        for (j = 0; j < nchar; ++j) {
            printf("%c", usrTipdata[j + i * nchar]);
        }
        printf("\n");
    }
    printf("\n");
    
    charstate *morphyTipdata = (charstate*) malloc(ntax * nchar * sizeof(charstate));
    
    /* Convert the tip data to a form that mfl_fitch_postorder can use */
    for (i = 0; i < nchar; ++i) {
        for (j = 0; j < ntax; ++j) {
            if (usrTipdata[i + j * nchar] == '?') {
                morphyTipdata[i + j * nchar] = -1;
            }
            else if (usrTipdata[i + j * nchar] == '-') {
                morphyTipdata[i + j * nchar] = 1;
            }
            else {
                morphyTipdata[i + j * nchar] = 1 << (usrTipdata[i + j * nchar] - '0' + 1);
            }
        }
    }
    
    int *currentchar = &j;
    
    /* Initialize the tree array that will store optimal trees */
    tree **savedtrees = (tree**) malloc(treelimit * sizeof(tree*));
    
    // Start with a (random, in this case) starting tree. 
    // Different algorithms will be written for doing this (as there are better
    // ways to do it) but we'll use randunrooted for now.
    
    savedtrees[0] = randunrooted(ntax, numnodes);
    //Get a length for the starting tree.
    mfl_root_tree(savedtrees[0], 1, ntax);
    int *besttreelen_p = &besttreelen;

    *besttreelen_p = mfl_get_treelen(savedtrees[0], morphyTipdata, ntax, nchar);
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
    
    // The exact type of rearrangment alogrithm (NNI, SPR, or TBR) will be set by the user
    
    mfl_nni_search(ntax, nchar, numnodes, morphyTipdata, savedtrees, besttreelen);
    
    
    //mfl_clear_treebuffer(savedtrees, 1, numnodes);
    
    for (i = 0, *currentchar = 0; i < nchar; ++i) {
        mfl_apply_tipdata(anewtree, morphyTipdata, ntax, nchar, *currentchar);
        *currentchar = *currentchar + 1;
        mfl_fitch_postorder(anewtree->root, treelength_p);
    }
    
    for (i = 0, *currentchar = 0; i < nchar; ++i) {
        mfl_apply_tipdata(arandomtree, morphyTipdata, ntax, nchar, *currentchar);
        *currentchar = *currentchar + 1;
        mfl_fitch_postorder(arandomtree->root, treelength_q);
    }

    printf("The 'user' tree:\n");
    printNewick(anewtree->root);
    printf("\n");
    printf("The length of the user tree: %i steps\n\n", *treelength_p);
    printf("The random tree:\n");
    printNewick(arandomtree->root);
    printf("\n");
    printf("The length of the random tree: %i steps\n\n", *treelength_q);
}

void testNWKreading(void)
{
    bool isRooted = true;    
    char *nwktree;
    tree *anewtree;
    
    char newickTree1[] = "(2,(1,4,3,5,6));";  
    printf("The newick string: %s\n", newickTree1);
    char newickTree2[] = "(2,((1,4,3),(5,6)));";  
    printf("The newick string: %s\n", newickTree2);
    char newickTree3[] = "(((1,2,4),3),(5,6));";  
    printf("The newick string: %s\n", newickTree3);
    char newickTree4[] = "((((1,2),4),3),(5,6));";  
    printf("The newick string: %s\n", newickTree4);
    
    nwktree = newickTree1;
    anewtree = readNWK(nwktree, isRooted);
    freetree(anewtree, 2*6-1);
    nwktree = newickTree2;
    anewtree = readNWK(nwktree, isRooted);
    freetree(anewtree, 2*6-1);
    nwktree = newickTree3;
    anewtree = readNWK(nwktree, isRooted);
    freetree(anewtree, 2*6-1);
    nwktree = newickTree4;
    anewtree = readNWK(nwktree, isRooted);
    freetree(anewtree, 2*6-1);
}

int main(void)
{
    int ntax = 12;
    int numnodes;
    //bool isRooted = true;
    
    numnodes = numberOfNodes(ntax);
    
    tree *anewtree;    
    tree *originaltree;
    tree *copiedtree;
    
    mini_test_analysis();
    
    //test_nni(ntax, numnodes);
    //testNWKreading();
    
    /* This part is just for testing the collapseBiNode*/
    anewtree = randrooted(ntax, numnodes);
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
    
    mfl_collapse(anewtree->trnodes[ntax + 1], anewtree->trnodes); // Magic number just for testing
    printf("With collapsed node: ");
    printNewick(anewtree->root);
    printf("\n");    
    //dump_nodearray(anewtree->trnodes, ntax, numnodes);
    
    copiedtree = copytree(anewtree, ntax, numnodes);
    //dump_nodearray(copiedtree->trnodes, ntax, numnodes);
    printf("Copying with collapsed node: ");
    printNewick(copiedtree->root);
    printf("\n");
    
    freetree(copiedtree, numnodes);
    
    mfl_arb_resolve(anewtree->trnodes[ntax], anewtree->trnodes, ntax, numnodes); // Magic number just for testing
    printf("With resolved node: ");
    printNewick(anewtree->root);
    printf("\n");
    //dump_nodearray(anewtree->trnodes, ntax, numnodes);
    
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
    freetree(copiedtree, numnodes);
    
    return 0;
}
