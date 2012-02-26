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

void print_bipartition(taxbipart bipartition, int ntax)
{
    int i = 0;
    //printf("%i\n", bipartition);
    while (bipartition) {
        
        ++i;
        if (bipartition & 1) 
        {
            printf("*");
        }
        else 
        {
            printf(".");
        }
        bipartition = bipartition >> 1;
        if (bipartition == 0) {
            while (ntax - i) {
                ++i;
                printf(".");
            }
        }
    }
    printf("\n");
}

void print_hashtab(taxbipart **hashtab, int ntax)
{
    int i;
    int numbiparts = ntax - 1;
    taxbipart tempht;
    
    for (i = 0; i < numbiparts; ++i) {
        tempht = *hashtab[i];
        print_bipartition(tempht, ntax);
    }
    //printf("\n");
}

void print_f_states(node *n, int nchar)
{
    int i, j;
    charstate c;
    if (!n->tip) {
        printf("break\n");
    }
    
    for (i = 0; i < nchar; ++i) {
        printf("C%i: ", i + 1);
        //printf("as int % i, ", c);
        c = n->apomorphies[i];
        
        if (c == -1) {
            printf("?; ");
            continue;
        }
        else 
        {
            if (c & 1) {
                printf("-");
            }
            if ((c & IS_APPLIC)) {
                j = 0;
                c = c >> 1;
                while (c) 
                {
                    if (c & 1) {
                        printf("%i", j);
                    }
                    ++j;
                    c = c >> 1;
                }
            }
        }
        printf("; ");
    }
    printf("\n");
}

void print_final_allviews(tree *testtree, int ntax, int nchar, int numnodes)
{
    int i, j;
    
    for (i = ntax + 1; i < numnodes; ++i) {
        //print_f_states(testtree->trnodes[i], nchar);
        node *q = testtree->trnodes[i];
        printf("node %i\n", i);
        do {
            if (q->apomorphies) {
                for (j = 0; j < nchar; ++j) {
                    printf("%u ", q->apomorphies[j]);
                }
                printf("\n");
            }
            q = q->next;
        } while (q != testtree->trnodes[i]);
    }
}


void print_dp_states(node *n, int nchar)
{
    int i, j;
    charstate c;
    if (!n->tip) {
        printf("break\n");
    }
    
    for (i = nchar-1; i < nchar; ++i) {
        printf("C%i: ", i + 1);
        //printf("as int % i, ", c);
        c = n->tempapos[i];
        
        if (c == -1) {
            printf("?; ");
            continue;
        }
        else 
        {
            if (c & 1) {
                printf("-");
            }
            if ((c & IS_APPLIC)) {
                j = 0;
                c = c >> 1;
                while (c) 
                {
                    if (c & 1) {
                        printf("%i", j);
                    }
                    ++j;
                    c = c >> 1;
                }
            }
        }
        printf("; ");
    }
    printf("\n");
}


void print_nodedata(node *n, int nchar, int ntax)
{
    node *p;
    
    if (n->tip) {
        printf("Tip: %i\n", n->tip);
        print_f_states(n, nchar);
        printf("\n");
        return;
    }
    
    p = n->next;
    while (p != n) {
        print_nodedata(p->outedge, nchar, ntax);
        p = p->next;
    }
    
    printf("Node: %i\n", n->index);
    if (n->tipsabove) {
        print_bipartition(*n->tipsabove, ntax);
    }
    printf("nodelen %i\n", n->nodelen);
    print_dp_states(n, nchar);
    print_f_states(n, nchar);
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
    newtree->hashtabholder = NULL;
    
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
        if (i < ntax) 
        {
            newtree->trnodes[i]->tip = i + 1;
            newtree->trnodes[i]->vweight = 1;
        }
        newtree->trnodes[i]->index = i;
    }
    
    newtree->root = NULL;
    newtree->hashtabholder = NULL;
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
    r1->clip = r2->clip = r3->clip = false;
    
    r2->apomorphies = r3->apomorphies = r1->apomorphies;
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
    if (r1->tempapos) {
        free(r1->tempapos);
        r1->tempapos;
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
        if (p->tempapos) {
            free(p->tempapos);
        }
        /*if (p->apomorphies) {
            free(p->apomorphies);
        }*/
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
        q->visited = p->visited;
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
    if (origtr->hashtabholder) {
        treecp->bipartitions = origtr->hashtabholder;
        origtr->hashtabholder = NULL;
    }
    //treecp->bipartitions = mfl_tree_biparts(treecp, ntax, numnodes);
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

node *mv_onetwo_to_twelve(tree *t, int ntax, int numnodes)
{
    int i;
    node *p;
    node *onetwonode = NULL;
    
    /*mfl_remove_branch(t->trnodes[0]);
    mfl_insert_branch(t->trnodes[0], t->trnodes[11], ntax);
    mfl_remove_branch(t->trnodes[1]);
    mfl_insert_branch(t->trnodes[1], t->trnodes[0], ntax);*/
    
    //Find the (1,2) node
    for (i = ntax + 1; i < numnodes; ++i) {
        p = t->trnodes[i];
        do {
            if ((p->outedge->tip == 11) || (p->outedge->tip == 12)) {
                onetwonode = t->trnodes[i];
                printf("found: %i\n", p->outedge->tip);
                break;
            }
            p = p->next;
            if (onetwonode) {
                break;
            }
        } while (p != t->trnodes[i]);
    }
    
    if (onetwonode->outedge->tip == 11 || onetwonode->outedge->tip == 12) {
        printf("Must find 'bottom' node\n");
    }

    return onetwonode;    
}

void test_char_optimization(void)
{
    char trstring[] = "(7,((8,(9,10)),(((((1,2),3),(4,(11,12))),5),6)));";
    tree *testtree = readNWK(trstring, 1);
    
    printf("The test tree:\n");
    printNewick(testtree->root);
    printf("\n");
    
    int i, j, al=0;
    int ntax = 12;
    int nchar = 21;
    int trl = 0;
    int *trlength = &trl;
    int numnodes = numberOfNodes(ntax);
    charstate *tipdata;
    node *mvnode;
    
    //mfl_start_drawtree(testtree, ntax);
    
    char usrTipdata[] = "111111111100000000001\n\
111111111100000000002\n\
001111111100000000002\n\
000011111100000000000\n\
000000111100000000001\n\
000000001100000000002\n\
000000000000000000110\n\
000000000000000011110\n\
000000000000001111111\n\
000000000000111111111\n\
000000000011111111112\n\
000000000011111111110";
    
    //old c21 pattern: 122012001120
    
    tipdata = mfl_convert_tipdata(usrTipdata, ntax, nchar);
    mfl_apply_tipdata(testtree, tipdata, ntax, nchar);
    
    mvnode = mv_onetwo_to_twelve(testtree, ntax, numnodes);
    mfl_undo_temproot(ntax, testtree);
    trl = mfl_all_views(testtree, ntax, nchar, trlength);
    //mfl_root_tree(testtree, 0, ntax);
    mfl_temproot(testtree, 5, ntax);
    printf("The altered tree (if move made):\n");
    printNewick(testtree->root);
    printf("\n");
    //mfl_undo_temproot(ntax, testtree);
    //testtree->bipartitions = mfl_tree_biparts(testtree, ntax, numnodes);
    printf("tree length: %i\n", trl);
    //mfl_undo_temproot(ntax, testtree);
    
    print_final_allviews(testtree, ntax, nchar, numnodes);
    
    int pos = 3;
    int diff = 0;
    node *up, *dn;
    up = mvnode->outedge->next->outedge;
    dn = mvnode->outedge->next->next->outedge;
    
    printf("mvnode index: %i\n", mvnode->index);
    printf("mvnode base:  %i\n", mvnode->outedge->index);
    
    joinNodes(up, dn);
    printf("The pruned tree:\n");
    printNewick(testtree->root);
    printf("\n");
    up->clip = true;
    dn->clip = true;

    mfl_undo_temproot(ntax, testtree);
    //mfl_wipe_states(testtree->trnodes[ntax], nchar);
    mfl_trav_allviews(testtree->trnodes[0], testtree, ntax, nchar, NULL, NULL);
    
    diff = mfl_get_treelen(testtree, ntax, nchar, NULL);
    printf("length of three after clipping: %i\n", diff);
    
    diff = mfl_get_subtreelen(mvnode, ntax, nchar, NULL);
    printf("length of sourcetree: %i\n", diff);
    
    diff = mfl_subtr_reinsertion(mvnode, testtree->trnodes[pos], testtree->trnodes[pos]->outedge, nchar);
    printf("cost of rejoining to original place: %i\n", diff);
    
    diff = mfl_locreopt_cost(mvnode, testtree->trnodes[9], testtree->trnodes[9]->outedge, nchar, 10000);
    printf("cost of rejoining to a new place: %i\n", diff);
    
    printf("target index: %i\n", testtree->trnodes[pos]->outedge->index);
    printf("target outedge: %i\n", testtree->trnodes[pos]->index);
    
    print_final_allviews(testtree, ntax, nchar, numnodes);
    
    joinNodes(mvnode->outedge->next, testtree->trnodes[9]->outedge);
    joinNodes(mvnode->outedge->next->next, testtree->trnodes[9]);
    //mfl_undo_temproot(ntax, testtree);
    diff = mfl_all_views(testtree, ntax, nchar, NULL);
    printf("length of the tree after placement: %i\n", diff);
    
    print_final_allviews(testtree, ntax, nchar, numnodes);
    
    mfl_temproot(testtree, 6, ntax);
    printf("The altered tree (if move made):\n");
    printNewick(testtree->root);
    printf("\n");
    
    //freetree(testtree, ntax);
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

    char usrTipdata[] = "00000000000000000000000---0000000-0000000000000000000-0-0000000000000000000--00000000000000000?\n\
11-10121001120200001001---0000100-1--10001000001101010320010000010100001000--30300000000000010?\n\
??????????????????????30--0010??0-????????1-?-??100?0-??1?--????0?00??10??1112??????100111-0-1?\n\
00??0?????????????????3111100?????????????011100??0???????11?1????????10?????2??????1?1-1??????\n\
000?0???????11??01111030--???1010-???0????1---0-10010-3?00--??010?000??0??111?1-?0?111011????1?\n\
11-101201011212?0001002---0001100-01010011000001101010320010000010101001000--10310000000000000?\n\
00?100(02)???01?1(23)0???10?0---00?00010002?001100000001000-31-010????-?30?000000--01-0???00?0000100?\n\
1011023111112131???1100---000000101--200010000101000103011?01-1-01001101010?-11-10??00?0?00010?\n\
?????????????????0????3000?10?010-0100?01100-0--100011300?10????0????0010110021-1???1001101000?\n\
10110231101111311001110---0000000-1---000100001010000-3011101-1-01311001010--11-001000000000100\n\
????????????????0?????3000?1?0??0-010?1?1100100?100?11??0?00????0??0?000111002??????10011010001\n\
????????????????0001??3101?100????0??????1001100100011??0?010100??????00??1002????1?10??????001\n\
????????????2????11?????--????????????????1-?-??1?????????--?????0?00?10?????2?????110011????10\n\
11-11?301?1111(23)?0010002---000000110121011100000111001011-010?0001?011000000--?021?0000000001001\n\
11-11?30012120000000001---000000110021011100000101000-11-000000010010001000--101110?00000000000\n\
0011012?1011?12??011002---00001010002?001100000?101010??1??0??001?10?00100??-???1??0???????0000\n\
1001??30??112???001001?---0000000-1--2011100000111000-2?-010000010211001000--0021?0000000000000\n\
?????????????????11???30--0010010-???????11-?-??10010-3010--0001?0000010?????2??????10011??0-10\n\
0???0?1??1?0?1??00????311111??010-0100?011001100100011??0?01??????0??00011???21-101?1001(01)010001\n\
00000?101110113000011030000100010-0100101100100010001?30001000000000010{01}1110021-101110?101-0001\n\
00100000101011101000110---0-000010011000010000101000103011100?0?00?0?101010--11-0?1?00000000000\n\
00000?1?1??0113?00011031111101010-01?0?0??011100100?0-3000010100?0?0??10?110(01)21-?011111-1010?11\n\
?00?????1?1?213?01?1????--0011????????????1-?-??10010-??0?--??0??0000??0??????????????????????0\n\
11-11?3???1?2???000000?---000000110111011100000111001011-010000010011000000--00211??00000001000\n\
00?00000?01021100000002---000010111--1001000000110101032001000001000?001000--10010??00?00000000\n\
??????????????????????30--????????????????1-?????0????????????0????0??????????????????????????0\n\
?????????????????1?1??3000?0?1????010????1001???100011??0??000??0??0?000??11021-1??1100111-0001\n\
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
    for (i = 0; usrTipdata[i]; ++i) {
        printf("%c", usrTipdata[i]);
        /*for (j = 0; j < nchar; ++j) {
            printf("%c", usrTipdata[j + i * nchar]);
        }
        printf("\n");*/
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
    
    pauseit();
    
    //test_tree_comparison();
    test_char_optimization();
    
    numnodes = numberOfNodes(ntax);
    
    //tree *anewtree;    
    //tree *originaltree;
    //tree *copiedtree;
    
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
    
    pauseit();
    
    return 0;
}
