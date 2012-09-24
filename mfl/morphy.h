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
#include <stdint.h>
#include <time.h>
#include <assert.h>
#include "mfl.h"
//#include <gsl/gsl_rng.h>


#include <iostream>
using namespace std;

#ifdef MFY_DEBUG
#define dbg_printf(...) printf(__VA_ARGS__)
#else
#define dbg_printf(...)
#endif

#define MORPHY_DEFAULT_TREE_LIMIT 200
#define IS_APPLIC ~1//(-1^1)
#define MORPHY_MAX_STATES 32


#define TREELIMIT 50000 //A temporary tree limit for testing.


/* For node and tree structures, this program follows the format recommended by
 * Felsenstein (2004. Inferring Phylogenies. Sinauer, Mass.) and implemented in 
 * Felsenstein et al.'s 2004. Phylip package. This includes representing 
 * internal nodes as a ring of (minmally) three nodes joined by the next 
 * pointer. The outedge pointer joins the node to either a leaf or the nearest 
 * internal node. */


typedef struct {
    long int n_rearrangements; // Number of tree topologies visited
    int n_savetrees;    // Number of trees found and/or saved
    int bestlength;     // Length of the best tree
    time_t searcht;     // Time in search
    string **newicktrees; // Array of trees in newick format
    string **texttrees;   // Array of trees drawn as text
#ifdef VERSION_1_5
    double consist_ind; // Consistency index
    double ret_ind;     // Retention index
    double resc_ci;     // Rescaled consistency index
    int num_islands;    // Number of parsimony islands hit
    int *times_hit      // List of times each island was hit
#endif    
} mfl_resultant_data_s;


typedef struct {
    int                     n_taxa;
    int                     n_chars;
    mfl_search_t            search_type;
    int                     n_iterations;
    int                     n_treelimit;
    mfl_branch_swap_t       bswap_type;
    bool                    is_ratchet;
    char                    *input_data;
    mfl_add_sequence_t      addseq_type;
    bool                    collapse_nolen;
    mfl_set_collapse_at_t   collapse_at;
    mfl_gap_t               gap_as_missing;
    mfl_resultant_data_s    *resultant_data;
} mfl_handle_s;

typedef int32_t charstate;
typedef int64_t taxbipart;

typedef struct node {
    struct node *outedge, *next;
    char *tipname;
    int tip;
    int index;
    int visited;
    int vweight;
    int initialized;
    int order;
    int nodelen;
    int cpindex;
    bool success;
    bool finished;
    bool start;
    bool skip;
    bool clip;
    bool bottom;
    bool isroot;
    int minsteps;
    int maxsteps;
    int charstates;
    taxbipart *tipsabove;
    charstate *tempapos;
    charstate *apomorphies;
} node;

typedef node **nodearray;

typedef struct tree {
    nodearray trnodes;
    /* An array of pointers to nodes. Elements from 0 to ntax - 1 are the 
     * terminals. Element ntax is generally reserved for the root, but this is 
     * just a convention. All other nodes are internal nodes and not the root 
     * node. */
    node *root;
    int length;
    int templen;
    bool swapped;
    int *compressedtr;
    int *cmptrholder;
    taxbipart **bipartitions;
    taxbipart **hashtabholder;
    int index;
} tree;

/*typedef struct chardata {
    int charnum;
    charstate *transeries;
    int optim_type;
    int maxvalue;
    int minvalue;
    bool included;
    bool informative;
    double cIndex;
    double rcIndex;
    double retIndex;
    double hIndex;
    int numstates;
    int *stepmatrix[MORPHY_MAX_STATES][MORPHY_MAX_STATES];
    void (*optimzation_algo)(node *n, int *trlength);
} chardata;*/

typedef struct {
    
    /* Keeps a record of the state of tree rearrangment.*/
    
    long int nextinbuffer;  // Next free position in the array of stored trees
    bool undertreelimit;
    int bestlength;         // Best tree length over all iterations
    int bestinrep;          // Best tree length in current iteration of search
    int trbufstart;         // Position in buffer of first tree to undergo rearrangement
    int tipscoll;
    bool foundbettertr;
    bool success;
    long int niter_total;
    long int niter_ontree;
} mfl_searchrec;

/*Function prototypes*/

/*in main.c*/
void call_index(node *n);
void dump_nodearray(nodearray nds, int ntax, int numnodes);
void dump_connections(nodearray nds, int ntax, int numnodes);
void dump_tree(tree *t, int ntax, int numnodes);
void print_bipartition(taxbipart bipartition, int ntax);
void print_hashtab(taxbipart **hashtab, int ntax);
void print_charstates(node *n, int nchar);
void print_final_allviews(tree *testtree, int ntax, int nchar, int numnodes);
void printNewick(node *n);
void treelen(node *n, int *stepcount); // The traversal algorithm that calls fitchdown

/*in compare.c*/
int mfl_compare_ints(const void * a, const void * b);
int mfl_compare_ints2(const void * a, const void * b);
int mfl_count_fields(int ntax);
bool mfl_comp_bipartition(taxbipart *bp1, taxbipart *bp2, int numfields);
void mfl_set_bipartition(node *n, node *d);
void mfl_set_tipsabove(node *n, int numfields, taxbipart **hashtab, int *bpcounter);
void mfl_free_hashtab(taxbipart **hashtab, int numbiparts);
taxbipart **mfl_tree_biparts(tree *t,int ntax, int numnodes);
bool mfl_compare_trees(taxbipart **t1, taxbipart **t2, int ntax, int numfields);
bool mfl_compare_alltrees(tree *newtopol, tree **savedtrees, int ntax, int numnodes, long int *start, long int *last);
void test_tree_comparison(void);
int *mfl_compress_tree(tree *t, int ntax, int numnodes);

/*in coptim.c*/
charstate * mfl_convert_tipdata(char *txtsrc, int ntax, int nchar, bool na_as_missing);
void mfl_apply_tipdata(tree *currenttree, charstate *tipdata, int ntax, int nchar);
void mfl_reopt_subtr(node *src, int nchar);
void mfl_reopt_subtr_root(node *n, int nchar);
void mfl_subtree_count(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *trlength);
void mfl_subtree_postorder(node *n, int *trlength, int nchar);
int mfl_reopt_shortcut(node *src, node *t1, node *t2, int nchar, int *trlength, int templen);
void mfl_countsteps(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *trlength, int *besttreelen);
int mfl_get_treelen(tree *testtree, int ntax, int nchar, int *bestreelen);
int mfl_get_subtreelen(node *n, int ntax, int nchar, int *besttreelen);
void mfl_fitch_postorder(node *n, int *trlength, int nchar, int *besttreelen);
void mfl_combine_up(node *n, node *anc, int nchar);
void mfl_fitch_preorder(node *n, int nchar);
int mfl_all_views(tree *t, int ntax, int nchar, int *besttreelen);
int mfl_get_sttreelen(tree *testtree, charstate *tipdata, int ntax, int nchar, int *besttreelen);
void mfl_reopt_postorder(node *n, int nchar);
void mfl_reopt_preorder(node *n, int nchar);
int mfl_locreopt_cost(node *src, node *tgt1, node *tgt2, int nchar, int diff);
int mfl_subtr_reinsertion(node *src, node *tgt1, node *tgt2, int nchar);
void mfl_wipe_states(node *n, int nchar);
void mfl_tip_apomorphies(node *tip, node *anc, int nchar);
void mfl_tip_reopt(tree *t, int ntax, int nchar);
void mfl_subtr_allviews(node *n, tree *t, int ntax, int nchar, int *treelen, int *besttreelen);
void mfl_allviews_traversal(node *n, tree *t, int ntax, int nchar, int *treelen, int *besttreelen);
void mfl_trav_allviews(node *n, tree *t, int ntax, int nchar, int *treelen, int *besttreelen);
void mfl_set_rootstates(node *n, int nchar);

/*in drawtree*/
void mfl_draw_tree(node *n, int *p_in_node, int *depth);
void mfl_start_drawtree(tree *t, int ntax);

/*in exhaustive.c*/ 
void allunrooted(void /*tree *treearray, int ntaxa*/);
void insert_allp(node *n, tree *origtree, int taxon, int calln, int *counter);
long long int factorial(long long int n);
long long int numtrees(int ntaxa);

/*in randtree.c*/
bool mfl_headtail (void);
void mfl_init_taxarray(int *taxarray, int ntax);
void mfl_shuffle(int *taxarray, int ntax);
void allunrooted(void /*tree *treearray, int ntaxa*/);
void insert_allp(node *n, tree *origtree, int taxon, int calln, int *counter);
long long int factorial(long long int n);
long long int numtrees(int ntaxa);
struct tree *randrooted (int ntax, int numnodes);
struct tree *randunrooted (int ntax, int numnodes);
struct node *mfl_breaktie_intryall(node *bestpos, node *n);
struct node * mfl_tryall(node *n, node *newbranch, node *bestpos, int ntax, int nchar, 
                         int numnodes, int *bestlen, tree *starttree, tree **savedtrees, charstate *tipdata);
tree * mfl_addseq_randasis(int ntax, int nchar, int numnodes, 
                                 charstate *tipdata, bool addRandom, tree** savedtrees);
/*in taxpart*/
int strToInt (char string[]);
void wipe_Og(int outtaxa[], nodearray outgroup);
void wipe_Ig(int intaxa[], nodearray ingroup);
void defOutgroup(int ntax, int outtaxa[], nodearray outgroup, int intaxa[], nodearray ingroup, bool *OGdefined);

/*in tree.c*/
void mfl_close_all_rings(nodearray nds, int ntax, int numnodes);
int mfl_calc_numnodes(int ntax);
void mfl_join_nodes(node *n, node *p);
struct tree *mfl_alloctree(int ntax, int numnodes);
void mfl_freetree(tree *newtree, int numnodes);
struct tree *mfl_alloc_noring(int ntax, int numnodes);
struct node * mfl_allocnode(void);
void mfl_countsteps(node *leftdesc, node *rightdesc, node *ancestor, int nchar, int *trlength, int *besttreelen);
struct tree * mfl_copytree(tree *origtree, int ntax, int numnodes);
void mfl_newring(node *r1, int ntax);
void mfl_newring_to_order(node *r1, int order, int ntax);
void mfl_deletering(node *r1);
struct node * mfl_seek_internal(int ntax,int numnodes, node **nds);
struct node * mfl_seek_ringnode(node *n, int ntax);
unsigned long long int mfl_subtree_id(bool reset);
void mfl_set_vweight(node *n);
void mfl_close_ring(node *n);
void mfl_as_ring(node *n, int ntax);
void mfl_as_noring(node *n);
void mfl_reindex_tree(nodearray nds, int ntax, int numnodes);
void mfl_set_ring_to_n(node *n);
void mfl_reset_ring(node *n);
void mfl_join_apomorphies(node *n);
void mfl_collapse(node *n, nodearray nds);
int mfl_determ_order(node *n);
void mfl_set_order(node *n);
void mfl_clear_order(node *n);
void mfl_set_index(node *n);
void mfl_devisit_tree(nodearray nds, int numnodes);
void mfl_put_branch_in_ring(node *n, node *rnode);
void mfl_insert_branch(node *br, node *target, int ntax);
void mfl_arb_resolve(node *n, node **nds, int ntax, int numnodes);
void mfl_collap_binode(node *n);
void mfl_definish_tree(tree *t, int numnodes);
void mfl_deinit_tree(tree *t, int numnodes);
void mfl_temproot(tree *trtoroot, int root, int ntax);
void mfl_undo_temproot(int ntax, tree *trtounroot);
int mfl_tree_enumerator(void);
void mfl_resize_treebuffer(tree **treebuffer, int *treelimit, int sizeincrease);
void mfl_clear_treebuffer(tree **treebuffer, long int *numsavedtrees, int numnodes);
void mfl_reinit_tbinrange(tree **treebuffer, tree *newbest, long int start, long int *endofrng, int numnodes);
void mfl_reinit_treebuffer(tree **treebuffer, tree *newbest, long int *numsavedtrees, int numnodes);
void mfl_point_bottom(node *n, node **nodes);
void mfl_root_tree(tree *trtoroot, int root, int ntax);
void mfl_unroot(int ntax, tree *rootedtree);
void mfl_collap_binode(node *n);
void mfl_save_newick(node *n, string *nwkstr);
string *mfl_trstring(tree *t, int ntax);

/*in readnewick.c*/
struct node * cpyfromNWK(char *nwktr, int nwklen, int ntax, int numnodes, int *pos, nodearray nds, bool isRooted);
void NWK_roothandl(char *nwktr, int nwklen, int ntax, int numnodes, tree *newtree, bool isRooted);
struct tree * readNWK (char *nwktr, bool isRooted);

/*in rearrange.c*/
mfl_searchrec * mfl_create_searchrec(void);
void mfl_part_reset_searchrec(mfl_searchrec *searchrec);
void mfl_destroy_searchrec(mfl_searchrec *searchrec);
long long int mfl_rearr_num(bool reset);
void mfl_bswap(node *p, node *q);
void mfl_remove_branch(node *n);
void mfl_insert_branch(node *br, node *target);
void mfl_nni_traversal(node *n, tree *swapingon, tree **savedtrees, int ntax, int nchar, int numnodes, long int *current, charstate *tipdata, bool *undertreelimit, int *currentbesttree, bool *foundbettertree);
void mfl_nni_search(int ntax, int nchar, int numnodes, charstate *tipdata, tree **savedtrees, int starttreelen);
void test_nni(int ntax, int numnodes);
long int mfl_spr_leftotry(int ntax);
void (*mfl_swap_controller(mfl_handle_s *mfl_handle)) (node*, tree*, tree**, int, int , int, mfl_searchrec*);
void mfl_regrafting_traversal(node *n, node *subtr, tree *swapingon, tree **savedtrees, int ntax, int nchar, int numnodes, mfl_searchrec *searchrec, int diff);
void mfl_pruning_traversal(node *n, tree *swapingon, tree **savedtrees, int ntax, int nchar, int numnodes, mfl_searchrec *searchrec);
node *mfl_find_atip(node *n);
void mfl_reroot_subtree(node *n, node *subtr, node *base, tree *swapingon, tree **savedtrees, int ntax, int nchar, int numnodes, mfl_searchrec *searchrec);
bool mfl_subtr_isrerootable(node *n);
bool mfl_heuristic_search(mfl_handle_s *mfl_handle/*int ntax, int nchar, int numnodes, char *txtsrcdata, tree **savedtrees, int starttreelen*/);
void mfl_spr_cliptrees(node *p, node *up, node *dn, node *subtr, tree *swapingon, tree **savedtrees, int ntax, int nchar, int numnodes, mfl_searchrec *searchrec);
void mfl_reroot_subtree(node *n, node *atip, node *subtr, node *base, node *up, node *dn, tree *swapingon, tree **savedtrees, int ntax, int nchar, int numnodes, mfl_searchrec *searchrec, int diff);
void mfl_bisection_traversal(node *n, tree *swapingon, tree **savedtrees, int ntax, int nchar, int numnodes, mfl_searchrec *searchrec);

/* in mfyinterface.c*/
bool mfl_set_ntax(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_set_nchar(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_set_searchtype(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_set_numiterations(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_set_treelimit(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_set_branchswap_t(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_set_ratchet_status(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_attach_inputdata(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_set_addseq_t(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_set_collapse(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_set_collapse_value(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_set_gapormissing(mfl_handle_s *mfl_struct, void *param_data);
bool mfl_set_parameter(mfl_handle_t *mfl_handle, mfl_param_t param_type, void *param_data);
mfl_resultant_data_s *mfl_morphy_controller(mfl_handle_s *mfl_handle);
mfl_handle_t mfl_s2t(mfl_handle_s *mfl_handle);
mfl_handle_s *mfl_t2s(mfl_handle_t mfl_handle);

/*End function prototypes*/
