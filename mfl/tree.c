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
    /* Loops over the internal nodes in the tree's node array, closing up any
     * rings which may have a dangling next pointer. */
    
    int i;
    
    for (i = ntax; i < numnodes; ++i) {
        mfl_close_ring(nds[i]);
    }
}

int mfl_calc_numnodes(int ntax)
{
    /* Returns the number of internal and terminal nodes possible for a fully 
     * binary tree with ntax number of terminal branches */
    
    int numnodes;
    
    return numnodes = 2 * ntax - 1;
}

void mfl_join_nodes(node *n, node *p)
{
    
    /* Joins two nodes via their edge pointers */
    
    if (n->edge) {
        n->edge->edge = NULL;
        //dbg_printf("terminating existing connection in n\n");
    }
    if (p->edge) {
        p->edge->edge = NULL;
        //dbg_printf("terminating existing connection in p\n");
    }
    
    n->edge = p;
    p->edge = n;
}

struct node * mfl_allocnode(void)
{
    
    /* Allocates a node */
    
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
    /* Allocates memory for a tree, but does not assign values to the edge
     * pointers. The returned structure is therefore not yet a true tree */
    
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
        
        newtree->trnodes[i]->edge = NULL;
    }
    
    newtree->root = NULL;
    newtree->hashtabholder = NULL;
    
    return (newtree);
}


void mfl_free_trnodes(tree *newtree, int numnodes)
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
        else {
            free(newtree->trnodes[i]->apomorphies);
            if (newtree->trnodes[i]->origtemps) {
                free(newtree->trnodes[i]->origtemps);
            }
            if (newtree->trnodes[i]->origfinals) {
                free(newtree->trnodes[i]->origfinals);
            }
        }
        
        free(newtree->trnodes[i]);
    }
    free(newtree->trnodes);
    newtree->trnodes = NULL;
}

void mfl_freetree(tree *newtree, int numnodes)
{
    /*
     * mfl_freetree - deletes an entire tree, by doing the mirror image
     * of mfl_alloctree.
     */
    
    if (newtree->trnodes) {
        mfl_free_trnodes(newtree, numnodes);
    }
    
    /* free the trnode list */
    if (newtree->compressedtr) {
        free(newtree->compressedtr);
    }
    if (newtree->cmptrholder) {
        free(newtree->cmptrholder);
    }
    if (newtree->newick_tree) {
        free(newtree->newick_tree);
    }
    
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
        if (!nds[i]->next /*&& !nds[i]->initialized*/ && !nds[i]->edge) {
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
                    if (p->edge) 
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
    /* Used by mfl_insert_branch() to find unoccupied nodes in the ringnode */
    
    /* Optimization: Lookit this crap. Why go through the process of checking if
     * it's the root, then checking again? */
    
    node *p;
    bool rootnode = false;
    
    if (n->index == ntax) {
        rootnode = true;
    }
    
    if (!rootnode && !n->edge) {
        return n;
    } else {
        p = n->next;
        while (p != n) {
            if (!p->edge) {
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
    
    /* Allocates memory for the node array in a tree struct, but does not create
     * the ring nodes for the internal nodes. Similar to mfl_alloctree() in that 
     * none of the nodes are joined. Used by mfl_copytree(), 
     * mfl_addseq_randasis(), and readNWK(). */
    
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
    
    /* Given a predicted number of branchings (order), returns an internal node
     * ring with that many branch points */
    
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

void mfl_deletering(node *r1)
{
    
    /*
     * mfl_deletering - deletes a node (ring) by doing the mirror image of newring.
     */
    
    if (r1->tempapos) {
        free(r1->tempapos);
    }
    if (r1->apomorphies) {
        free(r1->apomorphies);
    }
    if (r1->origtemps) {
        free(r1->origtemps);
    }
    if (r1->origfinals) {
        free(r1->origfinals);
    }
    
    node *p, *q;
    
    p = r1->next;
    while (p != r1) {
        q = p->next;
        if (p->tempapos) {
            free(p->tempapos);
        }
        if (p->origtemps) {
            free(p->origtemps);
        }
        if (p->origfinals) {
            free(p->origfinals);
        }
        free(p);
        p = q;
    }
}

struct tree * mfl_copytree(tree *origtr, int ntax, int numnodes)
{
    
    /* Given either a rooted or unrooted tree struct, origtree, copies the 
     * topology of that tree struct, but nothing else. However, the node indices
     * should be the same. This is not a guarantee, however, and should be 
     * tested before you try to make this assumption in any code */
    
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
        //q->visited = p->visited;
        q->skip = p->skip;
        if (p->next) 
        {
            if (p->order != q->order) 
            {
                dbg_printf("Error: template and copy have non-matching branch orders\n");
            }
            do 
            {
                if (p->edge && p->edge->tip) 
                {
                    mfl_join_nodes(q, treecp->trnodes[p->edge->index]);
                }
                if (p->edge && !p->edge->tip && !q->edge) 
                {
                    r = origtr->trnodes[p->edge->index];
                    s = treecp->trnodes[p->edge->index];
                    while (r->edge->index != p->index) 
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
    if (origtr->cmptrholder) {
        treecp->compressedtr = origtr->cmptrholder;
        origtr->cmptrholder = NULL;
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
        n->vweight = n->vweight + p->edge->vweight;
    }
}

void mfl_close_ring(node *n)
{
    /* Makes sure there isn't a dangling next pointer in an internal node ring */
    /*ADD AN ERROR CHECK*/
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
    
    /* Makes it so the all apomorphies pointers of each node in an internal ring 
     * node point to the same memory. */
    
    node *p;
    p = n->next;
    while (p != n) {
        p->apomorphies = n->apomorphies;
        p = p->next;
    }
}

void mfl_join_tempapos(node *n)
{
    
    /* Makes it so the all apomorphies pointers of each node in an internal ring 
     * node point to the same memory. */
    
    node *p;
    p = n->next;
    while (p != n) {
        p->tempapos = n->tempapos;
        p = p->next;
    }
}

void mfl_set_ring_to_n(node *n)
{
    /* Sets the initialized and skip values in the members of a ring to the 
     * same value as those in node n. */
    
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
    
    an1 = n1->edge;
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
    n1->edge = NULL;
    an1->edge = NULL;
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
    if (in->edge) {
        if (in->edge->edge) {
            in->edge->edge = NULL;
        }
        in->edge = NULL;
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
    if (p->next->edge) {
        p->next->edge = NULL;
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
    
    an1 = n->edge;
    an2 = an1->next;
    
    mfl_join_nodes(an1, n2->edge);
    
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
    
    if (n->edge) {
        i = 1;
    }
    
    if (n->tip || !n->next) {
        return i;
    }
    
    p = n->next;
    while (p != n) 
    {
        if (p->edge) {
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

void mfl_reset_nodes1(nodearray nds, int numnodes, int nchar)
{
    /* Resets all the flags used in a single round of optimization */
    
    int i;
    node *p;
    
    for (i = 0; i < numnodes; ++i) {
        
        p = nds[i];
        
        do {
            
            p->origbase = false;
            p->success = false;
            p->visited = false;
            p->finished = false;
            p->changed = false;
            p->clip = false;
            p->clippath = false;
            p->initialized = false;
            
            if (p->next) {
                p = p->next;
            }
        } while (p != nds[i]);
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

void mfl_undone_tree(nodearray nds, int numnodes)
{
    int i;
    node *p;
    
    for (i = 0; i < numnodes; ++i) {
        p = nds[i];
        do {
            p->done = false;
            if (p->next) {
                p = p->next;
            }
        } while (p != nds[i]);
    }
}

void mfl_definish_tree(tree *t, int numnodes)
{
    int i;
    node *p;
    for (i = 0; i < numnodes; ++i) {
        t->trnodes[i]->finished = false;
        t->trnodes[i]->changed = false;
        
        if (t->trnodes[i]->next) {
            p = t->trnodes[i]->next;
            while (p != t->trnodes[i]) {
                p->finished = false;
                p->changed = false;
                p = p->next;
            }
        }
    }
}

void mfl_desuccess_tree(tree *t, int numnodes)
{
    int i;
    node *p;
    for (i = 0; i < numnodes; ++i) {
        t->trnodes[i]->success = false;
        
        if (t->trnodes[i]->next) {
            p = t->trnodes[i]->next;
            while (p != t->trnodes[i]) {
                p->success = false;
                p = p->next;
            }
        }
    }
}

void mfl_erase_clippath(tree *t, int numnodes)
{
    int i;
    node *p;
    for (i = 0; i < numnodes; ++i) {
        t->trnodes[i]->clippath = false;
        
        if (t->trnodes[i]->next) {
            p = t->trnodes[i]->next;
            while (p != t->trnodes[i]) {
                p->clippath = false;
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
    /* Given a branch (two nodes joined by their edge pointers), places the 
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
    rtbr = lftbr->edge;
    
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
    
    lftbr = trtounroot->trnodes[ntax]->next->edge;
    rtbr = trtounroot->trnodes[ntax]->next->next->edge;
    
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
        mfl_point_bottom(p->edge, nodes);
        p = p->next;        
    }
    n->bottom = true;
    //nodes[n->index] = n;
}

void mfl_root_tree(tree *trtoroot, int nRoot, int ntax)
{    
    /*Roots the tree between a terminal (leaf) and an internal node*/
    
    node *nodeptr, *r2, *r3;
    
    if (!trtoroot->trnodes[ntax]->next) {
        mfl_newring(trtoroot->trnodes[ntax]->next, ntax);
    }
    
    nodeptr = trtoroot->trnodes[nRoot]->edge;
    r2 = trtoroot->trnodes[ntax]->next;
    r3 = trtoroot->trnodes[ntax]->next->next;
    
    mfl_join_nodes(r2, trtoroot->trnodes[nRoot]);
    mfl_join_nodes(r3, nodeptr);
    
    trtoroot->root = trtoroot->trnodes[ntax];
    trtoroot->trnodes[ntax]->edge = NULL;
    
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
            
            ftip = proot->next->edge;
            mfl_join_nodes(ftip, subnode);
            subnode->next = p;
            q->next = subnode;
            proot->next = n;
            m->next = proot;
            
        }
        else 
        {
            ftip = proot->next->edge;
            mfl_join_nodes(ftip, subnode);
            proot->next = NULL;
            subnode->next = p;
            q->next = subnode;
        }
    }
    else 
    {
        leftdesc = proot->next->edge;
        rightdesc = proot->next->next->edge;
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
        mfl_save_newick(n->edge, nwkstr);
        *nwkstr += string(")");
        return;
    }   
    
    if (n->tip && n->edge->next->edge) {
        if (n->edge->next->edge->tip && 
            !n->edge->next->edge->start) {
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
        mfl_save_newick(p->edge, nwkstr);
        p = p->next;
        if (p != n) {
            *nwkstr += string(",");
        }
    }   
    *nwkstr += string(")");
}

string mfl_trstring(tree *t, int ntax)
{
    string trstring;
    
    if (t->root) {
        trstring += string("= [&R] ");
        mfl_save_newick(t->root, &trstring);
    }
    else if (!t->root && t->trnodes[0]->start) {
        trstring += string("= [&U] ");
        mfl_save_newick(t->trnodes[0], &trstring);
    }
    else {
        dbg_printf("Error: tree has no starting node\n");
    }

    trstring += string(";");
    
    return trstring;
}

void mfl_save_newick_c(node *n, char *nwkstr)
{   
    /* Stores a treen in a c-style string in Newick format. In an unrooted tree, 
     * function will be called on a terminal node that has its start variable 
     * set to 1. */
    
    node *p;
    char cbuffer[64];
    
    if (n->start) {
        sprintf(cbuffer, "(%i,", n->tip);
        strcat(nwkstr, cbuffer);
        mfl_save_newick_c(n->edge, nwkstr);
        strcat(nwkstr, ")");
        return;
    }   
    
    if (n->tip && n->edge->next->edge) {
        if (n->edge->next->edge->tip && 
            !n->edge->next->edge->start) {
            sprintf(cbuffer, "%i", n->tip);
            strcat(nwkstr, cbuffer);
            return;
        }
    }
    
    if (n->tip && !n->start) {
        sprintf(cbuffer, "%i", n->tip);
        strcat(nwkstr, cbuffer);
        return;
    }
    
    strcat(nwkstr, "(");
    
    p = n->next;
    while (p != n) {
        mfl_save_newick_c(p->edge, nwkstr);
        p = p->next;
        if (p != n) {
            strcat(nwkstr, ",");
        }
    }   
    strcat(nwkstr, ")");
}

char *mfl_alloc_newick(int ntax, bool rooted)
{
    /* Returns the number of chars*/
    
    int i = 0;
    int n = ntax;
    int numnodes = 2 * ntax - 1;
    int nwksize = ntax;
    char *nwkstr;
    
    /* Get number of digits for ntax*/
    while (n) {
        nwksize = nwksize + ntax - (10*i-1);
        n = n/10;
        ++i;
    }
    --nwksize;
    
    /*there is one comma and two brackets for every node in the tree*/
    nwksize = nwksize + (numnodes * 3);
    
    nwkstr = (char*)malloc((nwksize + 6) * sizeof(char)); // 6 is the size of the rooted/unrooted headers appended below
    memset(nwkstr, 0, (nwksize + 1) * sizeof(char));
    
    if (rooted) {
        strcat(nwkstr, "[&R] ");
    } 
    else {
        strcat(nwkstr, "[&U] ");
    }

    return nwkstr;
    
}

char *mfl_newick_cstring(tree *t, int ntax)
{
    char *nwkstr;
    
    if (t->root) {
        nwkstr = mfl_alloc_newick(ntax, true);
        mfl_save_newick_c(t->root, nwkstr);
    }
    else if (!t->root && t->trnodes[0]->start) {
        nwkstr = mfl_alloc_newick(ntax, false);
        mfl_save_newick_c(t->trnodes[0], nwkstr);
    }
    else {
        dbg_printf("Error: tree has no starting node\n");
    }
    
    strcat(nwkstr, ";");
    
    return nwkstr;
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
    //dbg_printf("Clearing tree buffer for replicate\n");
    //dbg_printf("Best tree length: %i\n", newbest->length);
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
