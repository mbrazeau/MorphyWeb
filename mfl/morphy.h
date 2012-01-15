/*
 *  morphy.h
 *  Morphy
 *
 *  Created by Martin Brazeau on 11-04-26.
 *  Copyright 2011. All rights reserved.
 *  Morphy is provided as is with no warranty of any kind.
 *
 *  
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
/*#include <gsl/gsl_rng.h>*/

#define MAX_OG_SIZE 20
#define MAX_IG_SIZE 500

typedef char *statearray;

/* For node and tree structures, this program follows the format recommended by
 * Felsenstein (2004. Inferring Phylogenies. Sinauer, Mass.) and implemented in 
 * Felsenstein et al.'s Phylip package. This includes representing internal nodes 
 * as a ring of (minmally) three nodes joined by the next pointer. The outedge 
 * pointer joins the node to either a leaf or the nearest internal node. */

typedef struct node {
    struct node *outedge, *next;
    char *tipname;
    int tip;
    int index;
    int initialized;
    int order;
    bool start;
    bool dummy;
    int minsteps;
    int maxsteps;
    int numstates; //number of states of a character reconstructed at that node
    statearray apomorphies;
} node;

typedef node **nodearray;

typedef struct tree {
    nodearray trnodes;
    /* An array of pointers to nodes. Elements from 0 to ntax - 1 are the terminals. Element
     * ntax is generally reserved for the root, but this is just a convention. All other nodes
     * are internal nodes and not the root node. */
    node *root;
    int length;
    int index;
} tree;

/*Function prototypes*/

/*in main.c*/

void init_taxarray(int *taxarray, int ntax);
void joinNodes(node *n, node *p);
struct tree *alloctree(int ntax, int numnodes);
struct tree *alloc_noring(int ntax, int numnodes);
struct node * allocnode(void);
void printNewick(node *n);
void treelen(node *n, int *stepcount); // The traversal algorithm that calls fitchdown
void fitchdown(node *leftdesc, node *rightdesc, node *ancestor, int *stepcount); // The Fitch process for the downpass
struct tree * copytree(tree *origtree, int ntax, int numnodes); // Calls growcopy to copy a template tree
void growcopy(node *templ, node *target, tree *newtree, int *iter); // Called by copytree. Copies tree in preorder
void newring(node *r1);
void deletering(node *r1);
void detree(node *n);
void detree2(nodearray trnptr);
void point_bottom(node *n, node **nodes, int *counter);
void rootOnTerminal(tree *trtoroot, int root, int ntax);
void unroot(int ntax, tree *rootedtree);

/*in exhaustive.c*/ 
void allunrooted(void /*tree *treearray, int ntaxa*/);
void insert_allp(node *n, tree *origtree, int taxon, int calln, int *counter);
long long int factorial(long long int n);
long long int numtrees(int ntaxa);

/*in random.c*/
void allunrooted(void /*tree *treearray, int ntaxa*/);
void insert_allp(node *n, tree *origtree, int taxon, int calln, int *counter);
long long int factorial(long long int n);
long long int numtrees(int ntaxa);
struct tree *randrooted (int ntax, int numnodes);
struct tree *randunrooted (int ntax, int numnodes);

/*in taxpart*/
int strToInt (char string[]);
void wipe_Og(int outtaxa[], nodearray outgroup);
void wipe_Ig(int intaxa[], nodearray ingroup);
void defOutgroup(int ntax, int outtaxa[], nodearray outgroup, int intaxa[], nodearray ingroup, bool *OGdefined);

/*in tree.c*/
struct node * seekInternal(int ntax, node **nds);
void closeRing(node *n);
void asRing(node *n);
void asNoring(node *n);
void mfl_collapse(node *n, nodearray nds);
int mfl_determ_order(node *n);
void mfl_set_order(node *n);
void mfl_clear_order(node *n);
void mfl_set_index(node *n);
void putBranchInRing(node *n, node *rnode);
void insertBranch(node *br, node *target);
void mfl_arb_resolve(node *n, node **nds, int ntax, int numnodes);
void mfl_deinit_tree(tree *t);

/*End function prototypes*/
