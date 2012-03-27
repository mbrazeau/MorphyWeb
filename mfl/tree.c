/*
 *  tree.c
 *  Morphy
 *
 *  Miscellaneous functions for trees
 *
 */

#include "morphy.h"

void mfl_close_all_rings(nodearray nds, int ntax, int numnodes)
{
    int i;
    
    for (i = ntax; i < numnodes; ++i) {
        mfl_close_ring(nds[i]);
    }
}

int mfl_calc_numnodes(int ntax)
{
    int numnodes;
    
    return numnodes = 2 * ntax - 1;
}


void mfl_join_nodes(node *n, node *p)
{
    if (n->outedge) {
        n->outedge->outedge = NULL;
        //dbg_printf("terminating existing connection in n\n");
    }
    if (p->outedge) {
        p->outedge->outedge = NULL;
        //dbg_printf("terminating existing connection in p\n");
    }
    
    n->outedge = p;
    p->outedge = n;
}

struct node * mfl_allocnode(void)
{
    node *newNode;
    newNode = (node *)malloc(sizeof(node));
    if (newNode == NULL)
    {
        dbg_printf("Error: failed to allocate new node.\n");
    }
    else
    {
        memset(newNode, 0, sizeof(node));
    }
    return newNode;
}

struct tree * mfl_alloctree(int ntax, int numnodes)
{
    int i; //Loop counters
    tree *newtree;
    
    newtree = (tree*)malloc(sizeof(tree));
    memset(newtree, 0, sizeof(tree));
    
    if (newtree == NULL)
    {
        dbg_printf("Error: failed to allocate new tree.\n");
        return (struct tree*) 0;
    }
    
    newtree->trnodes = (node **)malloc( (numnodes) * sizeof(node*));
    
    for (i = 0; i < numnodes; ++i)
    {
        newtree->trnodes[i] = mfl_allocnode();
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
            mfl_newring(newtree->trnodes[i], ntax);
        }
        
        newtree->trnodes[i]->outedge = NULL;
    }
    
    newtree->root = NULL;
    newtree->hashtabholder = NULL;
    
    return (newtree);
}

/*
 * mfl_freetree - deletes an entire tree, by doing the mirror image
 * of mfl_alloctree.
 */
void mfl_freetree(tree *newtree, int numnodes)
{
    int i;
    
    /* free all of the trnodes */
    for (i = 0; i < numnodes; ++i)
    {
        if (newtree->trnodes[i]->next)
        {
            mfl_close_ring(newtree->trnodes[i]);
            mfl_deletering(newtree->trnodes[i]);
        }
        free(newtree->trnodes[i]);
    }
    /* free the trnode list */
    free(newtree->trnodes);
    /* free the tree */
    free(newtree);
    
}

struct node * mfl_seek_internal(int ntax, int numnodes, node **nds)
{
    /* Searches for an unused internal node */
    /* NB: This function needs to return some kind of error msg
     * if no unused nodes are found. */
    
    int i;
    node *unused = NULL;
    node *p;
    bool isUsed = false;
    
    for (i = ntax + 1; i < numnodes; ++i) {
        if (!nds[i]->next /*&& !nds[i]->initialized*/ && !nds[i]->outedge) {
            unused = nds[i];
            i = 2 * ntax;
        }
    }
    
    if (!unused) 
    {
        for (i = ntax + 1; i < numnodes; ++i) 
        {
            isUsed = false;
            if (nds[i]->next) 
            {
                p = nds[i];                
                do
                {
                    if (p->outedge) 
                    {
                        isUsed = true;
                        p = nds[i];
                    }
                    else 
                    {
                        p = p->next;
                    }                
                } while (p != nds[i]);
                
                if (!isUsed) {
                    unused = nds[i];
                }
            }
        }
    }
    
    if (!unused) {
        dbg_printf("Error in tree memory allocation\n");
        return unused;
    }
    else {
        //unused->initialized = 1;
        return unused;
    }
}

struct node * mfl_seek_ringnode(node *n, int ntax)
{
    /*-Used by mfl_copytree-*/
    node *p;
    bool rootnode = false;
    
    if (n->index == ntax) {
        rootnode = true;
    }
    
    if (!rootnode && !n->outedge) {
        return n;
    } else {
        p = n->next;
        while (p != n) {
            if (!p->outedge) {
                return p;
            }
            p = p->next;
        }
    }
    
    dbg_printf("did not find an available node in ring\n");
    return n;
}

struct tree *mfl_alloc_noring(int ntax, int numnodes)
{
    int i;
    tree *newtree;
    
    newtree = (tree*)malloc(sizeof(tree));
    memset(newtree, 0, sizeof(tree));
    
    if (newtree == NULL)
    {
        dbg_printf("Error: failed to allocate new tree.\n");
        return (struct tree*) 0;
    }
    
    newtree->trnodes = (node **)malloc( (numnodes) * sizeof(node*));
    
    for (i = 0; i < numnodes; ++i)
    {
        newtree->trnodes[i] = mfl_allocnode();
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

void mfl_newring(node *r1, int ntax)
{
    /* Generates a new internal node composed of a ring of node structures
     * joined in a unidirectional manner by their next pointers. Used by any
     * function that might dynamically add branches at run-time. */
    
    node *r2, *r3;
    
    r2 = mfl_allocnode();
    r3 = mfl_allocnode();
    
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

void mfl_newring_to_order(node *r1, int order, int ntax)
{
    int i = 1;
    node *p;
    
    if (order == 0 || order == 1) {
        dbg_printf("Error in mfl_newring_to_order: no order specified or order too small\n");
        return;
    }
    
    p = r1;
    do {
        p->next = mfl_allocnode();
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
 * mfl_deletering - deletes a node (ring) by doing the mirror image of newring.
 */
void mfl_deletering(node *r1)
{
    if (r1->tipsabove) {
        free(r1->tipsabove);
    }
    if (r1->tempapos) {
        free(r1->tempapos);
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

struct tree * mfl_copytree(tree *origtr, int ntax, int numnodes)
{
    int i, tmpltorder, begin;
    tree *treecp; // Pointer to the tree copy
    node *p, *q, *r, *s;
    //bool isRooted = false;
    
    treecp = mfl_alloc_noring(ntax, numnodes);
    
    if (origtr->root)
    {
        mfl_set_order(origtr->trnodes[ntax]);
        tmpltorder = origtr->trnodes[ntax]->order;
        mfl_newring_to_order(treecp->trnodes[ntax], tmpltorder + 1, ntax);
        treecp->trnodes[ntax]->order = treecp->trnodes[ntax]->order - 1;
        treecp->root = treecp->trnodes[ntax];
        //isRooted = true;
        begin = ntax;
    }
    else 
    {
        treecp->trnodes[0]->start = true;
        begin = ntax + 1;
        mfl_newring(treecp->trnodes[ntax], ntax);
    }
    
    
    for (i = ntax + 1; i < numnodes; ++i) 
    {
        if (origtr->trnodes[i]->next) 
        {
            mfl_set_order(origtr->trnodes[i]);
            tmpltorder = origtr->trnodes[i]->order;
            mfl_newring_to_order(treecp->trnodes[i], tmpltorder, ntax);
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
                dbg_printf("Error: template and copy have non-matching branch orders\n");
            }
            do 
            {
                if (p->outedge && p->outedge->tip) 
                {
                    mfl_join_nodes(q, treecp->trnodes[p->outedge->index]);
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
                    mfl_join_nodes(q, s);
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

unsigned long long int mfl_subtree_id(bool reset)
{
    static unsigned long long subtreeid = 0;
    if (reset) {
        subtreeid = 0;
        return 0;
    }
    
    ++subtreeid;
    
    return subtreeid;
}


int mfl_tree_enumerator(void)
{
    static long long int treenum = 0;
    
    return ++treenum;
}


void mfl_set_vweight(node *n)
{
    node *p;
    p = n->next;
    while (p != n) {
        n->vweight = n->vweight + p->outedge->vweight;
    }
}

void mfl_close_ring(node *n)
{
    /* Makes sure there isn't a dangling next pointer*/
    node *p;
    
    p = n;
    
    do {
        if (p->next) {
            p = p->next;
        }
        if (!p->next) {
            dbg_printf("Report error: dangling next pointer at %p\n", p);
            p->next = n;
            p = p->next; 
        }
    }
    while (p != n);
}

void mfl_as_ring(node *n, int ntax)
{
    
    if (n->next) {
        mfl_close_ring(n);
        mfl_deletering(n);
    }
    
    mfl_newring(n, ntax);
}

void mfl_as_noring(node *n)
{
    if (n->next) {
        mfl_close_ring(n);
        mfl_deletering(n);
    }
}

void mfl_reindex_tree(nodearray nds, int ntax, int numnodes)
{
    int i;
    
    for (i = 0; i < numnodes; ++i)
    {
        /* First half of trnodes are initialized to this */
        if (i < ntax)
        {
            nds[i]->index = i;
        }
        else /* second half are allocated as a newring */
        {
            nds[i]->index = i;
            mfl_set_index(nds[i]);
        }        
    }
}

void mfl_join_apomorphies(node *n)
{
    node *p;
    p = n->next;
    while (p != n) {
        p->apomorphies = n->apomorphies;
        p = p->next;
    }
}

void mfl_set_ring_to_n(node *n)
{
    node *p;
    
    if (n->next) {
        p = n->next;
        while (p != n) 
        {
            p->index = n->index;
            //p->apomorphies = n->apomorphies;
            p->initialized = n->initialized;
            p->skip = n->skip;
            p = p->next;
        }
    }
    mfl_set_order(n);
    
}

void mfl_reset_ring(node *n)
{
    node *p;
    
    n->initialized = 0;
    n->clip = false;
    n->skip = 0;
    
    if (n->next) {
        p = n->next;
        while (p != n) 
        {
            p->index = n->index;
            p->initialized = n->initialized;
            p->clip = n->clip;
            p->skip = n->skip;
            p = p->next;
        }
    }
    mfl_set_order(n);    
}

void mfl_collapse(node *n1, nodearray nds)
{
    node *an1, *an2, *n2, *tmp;
    
    an1 = n1->outedge;
    an2 = an1->next;
    tmp = an2;
    while (an2->next != an1) {
        an2 = an2->next;
    }
    an2->next = n1->next;
    
    n2 = n1->next;
    while (n2->next != n1) {
        n2 = n2->next;
    }
    n2->next = an1->next;
    an1->next = NULL;
    n1->next = NULL;
    n1->outedge = NULL;
    an1->outedge = NULL;
    free(an1);
    nds[n1->index] = n1;
    mfl_set_order(tmp);
    mfl_set_vweight(tmp);
}

void mfl_arb_resolve(node *n, node **nds, int ntax, int numnodes)
{
    /* Arbitrarily resolves a non-binary node and leaves it binary*/
    
    int c, i, j, ord;
    node *p, *q, *in, *in2;
    
    // Make sure node there is a polytomy, otherwise exit resolve()
    mfl_set_order(n);
    ord = n->order;
    if (ord <= 3 && n->index != ntax) {
        dbg_printf("mfl_arb_resolve() called in error on node %i\n", n->index);
        return;
    }
    
    // Find available internal node(s) in nds
    in = mfl_seek_internal(ntax, numnodes, nds);
    in2 = in;
    if (in->next) {
        mfl_as_noring(in);
    }
    if (in->outedge) {
        if (in->outedge->outedge) {
            in->outedge->outedge = NULL;
        }
        in->outedge = NULL;
    }
    
    /* Randomly select branches to join to it. Cycles through the ring (skipping 
     * n itself) ntax % 10 times. This doesn't use a true random number 
     * generator, but should be both sufficiently arbitrary with respect to the
     * tree's topology, but sufficiently deterministic to be repeatable. */
    
    c = ntax % 10;
    
    for (i = 0; i <= (ord - 3); ++i)
    {
        p = n->next;
        for (j = 0; j <= c; ++j) 
        {
            q = p;
            p = p->next;
            if (p == n) 
            {
                q = p;
                p = p->next;
            }
        }
        in2->next = p;
        in2 = in2->next;
        q->next = p->next;
        p->next = NULL;
    }
    
    p->next = mfl_allocnode();
    if (p->next->outedge) {
        p->next->outedge = NULL;
    }
    p = p->next;
    p->next = in;
    
    if (ntax % 2) 
    {
        n = n->next;
    }
    else 
    {
        n = n->next->next;
    }
    
    mfl_insert_branch(in, n, ntax);
}

void mfl_collap_binode(node *n)
{
    /*collapses a binary node*/
    node *n2, *n3;
    node *an1, *an2;
    
    n2 = n->next;
    n3 = n2->next;
    
    an1 = n->outedge;
    an2 = an1->next;
    
    mfl_join_nodes(an1, n2->outedge);
    
    an1->next = n3;
    n3->next = an2;
    
    n2->index = an1->index;
    n2->index = an1->index;
    
    n->next = NULL;
    
}

int mfl_determ_order(node *n)
{
    /*determines the number of branchings in a node*/
    int i = 0;
    node *p;
    
    if (n->outedge) {
        i = 1;
    }
    
    if (n->tip || !n->next) {
        return i;
    }
    
    p = n->next;
    while (p != n) 
    {
        if (p->outedge) {
            ++i;
        }
        p = p->next;
    }
    
    return i;
}

void mfl_set_order(node *n)
{
    int ord;
    node *p;
    ord = mfl_determ_order(n);
    
    n->order = ord;
    if (n->next) 
    {
        p = n->next;
        while (p != n) {
            p->order = ord;
            p = p->next;
        }
    }
}

void mfl_clear_order(node *n)
{
    node *p;
    
    p = n->next;
    while (p != n) {
        p->order = 0;
        p = p->next;
    }
    
    n->order = 0;
}

void mfl_set_index(node *n)
{
    node *p;
    
    //dbg_printf("setting to index %i\n", n->index);
    
    if (n->next) {
        p = n->next;
        while (p != n) 
        {
            p->index = n->index;
            p = p->next;
        }
    }   
}

void mfl_devisit_tree(nodearray nds, int numnodes)
{
    int i;
    node *p;
    
    for (i = 0; i < numnodes; ++i) {
        nds[i]->visited = 0;
        if (nds[i]->next) {
            p = nds[i]->next;
            while (p != nds[i]) {
                p->visited = 0;
                p = p->next;
            }
        }
    }
}

void mfl_definish_tree(tree *t, int numnodes)
{
    int i;
    node *p;
    for (i = 0; i < numnodes; ++i) {
        t->trnodes[i]->finished = false;
        
        if (t->trnodes[i]->next) {
            p = t->trnodes[i]->next;
            while (p != t->trnodes[i]) {
                p->finished = false;
                p = p->next;
            }
        }
    }
}

void mfl_deinit_tree(tree *t, int numnodes)
{
    int i;
    node *p;
    
    for (i = 0; i < numnodes; ++i) 
    {
        t->trnodes[i]->initialized = 0;
        
        if (t->trnodes[i]->next) 
        {
            p = t->trnodes[i]->next;
            
            do {
                p->initialized = 0;
                if (p->next) 
                {
                    p = p->next;
                }
                else 
                {
                    p = t->trnodes[i]; /* Should prevent a crash in the event of a dangling next pointer*/
                }

            } while (p != t->trnodes[i]);
            
        }
    }
}

void mfl_put_branch_in_ring(node *n, node *rnode)
{
    /* Given a branch (two nodes joined by their outedge pointers), places the 
     * branch in a ring */
    
    node *rnode2 = NULL;
    
    if (rnode->next) {
        rnode2 = rnode->next;
    }
    
    rnode->next = n;
    
    if (rnode2) {
        n->next = rnode2;
    }
    else {
        n->next = rnode;
    }
    
    rnode->order = rnode->order + 1;
    mfl_set_order(rnode);
}


void mfl_temproot(tree *trtoroot, int root, int ntax)
{
    node *lftbr, *rtbr;
    
    lftbr = trtoroot->trnodes[root];
    rtbr = lftbr->outedge;
    
    if (!trtoroot->trnodes[ntax]->next) {
        mfl_newring(trtoroot->trnodes[ntax], ntax);
    }
    
    mfl_join_nodes(lftbr, trtoroot->trnodes[ntax]->next);
    mfl_join_nodes(rtbr, trtoroot->trnodes[ntax]->next->next);
    
    trtoroot->root = trtoroot->trnodes[ntax];
    trtoroot->trnodes[0]->start = false;
}

void mfl_undo_temproot(int ntax, tree *trtounroot)
{
    node *lftbr, *rtbr;
    
    lftbr = trtounroot->trnodes[ntax]->next->outedge;
    rtbr = trtounroot->trnodes[ntax]->next->next->outedge;
    
    mfl_join_nodes(lftbr, rtbr);
    
    trtounroot->root = NULL;
    trtounroot->trnodes[0]->start = true;
}

void mfl_point_bottom(node *n, node **nodes)
{
    /* Re-sets the pointers to the internal nodes in the tree's
     * node array to point to the 'bottom' (rootward) node in the ring*/
    
    node *p;
    
    if (n->tip) {
        n->bottom = false;
        return;
    }
    
    p = n->next;
    while (p != n) {
        p->bottom = false;
        mfl_point_bottom(p->outedge, nodes);
        p = p->next;        
    }
    n->bottom = true;
    nodes[n->index] = n;
}

void mfl_root_tree(tree *trtoroot, int nRoot, int ntax)
{    
    /*Roots the tree between a terminal (leaf) and an internal node*/
    
    node *nodeptr, *r2, *r3;
    
    if (!trtoroot->trnodes[ntax]->next) {
        mfl_newring(trtoroot->trnodes[ntax]->next, ntax);
    }
    
    nodeptr = trtoroot->trnodes[nRoot]->outedge;
    r2 = trtoroot->trnodes[ntax]->next;
    r3 = trtoroot->trnodes[ntax]->next->next;
    
    mfl_join_nodes(r2, trtoroot->trnodes[nRoot]);
    mfl_join_nodes(r3, nodeptr);
    
    trtoroot->root = trtoroot->trnodes[ntax];
    trtoroot->trnodes[ntax]->outedge = NULL;
    
    mfl_point_bottom(trtoroot->root, trtoroot->trnodes);
    trtoroot->trnodes[0]->start = false;
}

void mfl_unroot(int ntax, tree *rootedtree)
{
    int lnumnodes = 2*ntax-1;
    node *proot, *leftdesc, *rightdesc, *subnode, *n, *m, *p, *q, *ftip;
    
    proot = rootedtree->root;
    
    mfl_set_order(proot);
    
    if (proot->order > 2) 
    {
        dbg_printf("derooting multifurcating node\n");
        
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
            dbg_printf("Doing the trickier one\n");
            n = subnode->next;
            m = subnode;
            while (m->next != subnode) {
                m = m->next;
            }
            m->next = NULL;
            
            ftip = proot->next->outedge;
            mfl_join_nodes(ftip, subnode);
            subnode->next = p;
            q->next = subnode;
            proot->next = n;
            m->next = proot;
            
        }
        else 
        {
            ftip = proot->next->outedge;
            mfl_join_nodes(ftip, subnode);
            proot->next = NULL;
            subnode->next = p;
            q->next = subnode;
        }
    }
    else 
    {
        leftdesc = proot->next->outedge;
        rightdesc = proot->next->next->outedge;
        mfl_join_nodes(leftdesc, rightdesc);
    }
    
    rootedtree->root = NULL;
    rootedtree->trnodes[0]->start = true;
    
}

void mfl_save_newick(node *n, string *nwkstr)
{   
    /* Stores a treen in a string in Newick format. In an unrooted tree, 
     * function will be called on a terminal node that has its start variable 
     * set to 1. */
    
    node *p;
    char cbuffer[64];
    
    if (n->start) {
        sprintf(cbuffer, "(%i,", n->tip);
        *nwkstr += string(cbuffer);
        mfl_save_newick(n->outedge, nwkstr);
        *nwkstr += string(")");
        return;
    }   
    
    if (n->tip && n->outedge->next->outedge) {
        if (n->outedge->next->outedge->tip && 
            !n->outedge->next->outedge->start) {
            sprintf(cbuffer, "%i", n->tip);
            *nwkstr += string(cbuffer);
            return;
        }
    }
    
    if (n->tip && !n->start) {
        sprintf(cbuffer, "%i", n->tip);
        *nwkstr += string(cbuffer);
        return;
    }
    
    *nwkstr += string("(");
    
    p = n->next;
    while (p != n) {
        mfl_save_newick(p->outedge, nwkstr);
        p = p->next;
        if (p != n) {
            *nwkstr += string(",");
        }
    }   
    *nwkstr += string(")");
}

string mfl_trstring(tree *t, int ntax)
{
    string *trstring = new string;
    
    if (t->root) {
        mfl_save_newick(t->root, trstring);
    }
    else if (!t->root && t->trnodes[0]->start) {
        mfl_save_newick(t->trnodes[0], trstring);
    }
    else {
        dbg_printf("Error: tree has no starting node\n");
    }

    *trstring += string(";");
    
    return *trstring;
}

void mfl_resize_treebuffer(tree **treebuffer, int *treelimit, int sizeincrease)
{
    tree **newtreebuffer;
    newtreebuffer = (tree **)realloc(treebuffer, (*treelimit + sizeincrease) * sizeof(tree*));
    if (!newtreebuffer) {
        dbg_printf("Insufficient memory for new treelimit\n");
        newtreebuffer = treebuffer;
    }
    else {
        *treelimit = *treelimit + sizeincrease;
        dbg_printf("Treelimit changed to %i\n", *treelimit);
    }
}

void mfl_clear_treebuffer(tree **treebuffer, long int *numsavedtrees, int numnodes)
{
    int i;
    for (i = 0; i < *numsavedtrees; ++i) {
        if (treebuffer[i]) {
            mfl_freetree(treebuffer[i], numnodes);
        }
    }
    *numsavedtrees = 0;
}

void mfl_reinit_tbinrange(tree **treebuffer, tree *newbest, long int start, long int *endofrng, int numnodes)
{
    dbg_printf("Clearing tree buffer for replicate\n");
    dbg_printf("Best tree length: %i\n", newbest->length);
    int i;
    for (i = start; i < *endofrng; ++i) {
        if (treebuffer[i] && treebuffer[i] != newbest) {
            mfl_freetree(treebuffer[i], numnodes);
        }
        if (treebuffer[i] == newbest) {
            treebuffer[i] = NULL;
        }
    }
    treebuffer[start] = newbest;
    *endofrng = start;    
}

void mfl_reinit_treebuffer(tree **treebuffer, tree *newbest, long int *numsavedtrees, int numnodes)
{    
    dbg_printf("Found shorter tree. Reinitializing the treebuffer\n");
    dbg_printf("Best tree length: %i\n", newbest->length);
    int i;
    for (i = 0; i < *numsavedtrees; ++i) {
        if (treebuffer[i] && treebuffer[i] != newbest) {
            mfl_freetree(treebuffer[i], numnodes);
        }
        if (treebuffer[i] == newbest) {
            treebuffer[i] = NULL;
        }
    }
    
    newbest->index = 0;
    newbest->swapped = false;
    treebuffer[0] = newbest;
    *numsavedtrees = 0;
    mfl_subtree_id(true);
}
