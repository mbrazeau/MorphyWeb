/*
 *
 * THE MORPHY FUNCTION LIBRARY (MFL)
 *
 */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef MFY_DEBUG
#include <stdio.h>
#define dbg_printf(...) printf(__VA_ARGS__)
#else
#define dbg_printf(...)
#endif

/*
 *
 * Data type definitions
 *
 */


typedef uint64_t charstate; // Each character state is represented by a single unsigned 64-bit integer. Thus, one character may have 64 possible states.


typedef struct mfl_node_t {
	struct mfl_node_t *nodet_edge, *nodet_next; // Pointers to the neighboring node in the tree; the next node in the node ring.
	char *nodet_tipname;                        // Name of the tip from the dataset.
	int *nodet_tip_num;                         // 1-based identifier of terminal. Assigned 0 if node is internal.
	int nodet_index;                            // 0-based index of node in the node-array.
    int nodet_isbottom;                         // Indicates node points to (calculation) root.
    int nodet_downpass_visited;                 // Indicates successful visit from a downpass traversal.
	int nodet_uppass_visited;                   // Indicates successful visit from an uppass traversal.
    int nodet_minsteps;                         // Minimum number of transformations along branch represented by nodet_edge.
	int nodet_maxsteps;                         // Maximum number of transformations along branch represented by nodet_edge.
    double nodet_mean_branchlen;                // Mean branch length.
    int nodet_isroot;                           // Set non-zero if this node is the root.
    long long int nodet_tree_index;             // Identity of the tree to which this node belongs.
    charstate *nodet_original_firstpass;        // The downpass optimisations of the starting tree.
    charstate *nodet_original_secondpass;       // The uppass optimisations of the starting tree.
    charstate *nodet_subtree_firstpass_set;     // The downpass optimisations after the tree is broken.
	charstate *nodet_subtree_secondpass_set;    // The uppass optimisations after the tree is broken.
} mfl_node_t;


typedef mfl_node_t ** mfl_nodearray_t;  // A pointer of pointers to nodes used in trees.


typedef struct mfl_tree_t {
	mfl_nodearray_t treet_treenodes;        // The nodes of the tree. Never point these to another instance of mfl_tree.
	mfl_node_t *treet_root;                 // Pointer to the root of the tree.
    mfl_nodearray_t treet_outgroup_tips;    // Pointers to the outgroup tips.
	int treet_uw_parsimonylength;           // Unweighted number of steps under parsimony
    int treet_island_id;                    // An identification number for the parent start tree.
    int *treet_compressed_tree;             // The integer encoding of the tree for tree comparisons and saving memory
	int *treet_comptree_holder;             // Holder for comparisons.
	char *treet_newick_format;              // The tree stored as a Newick format string.
	int treet_index;                        // Identifier for the tree in the saved trees array.
} mfl_tree_t;


typedef enum mfl_optimisation_t {
    MFL_IS_FITCH,
    MFL_IS_WAGNER,
    MFL_IS_DOLLO,
    MFL_IS_IRREVERSIBLE,
    MFL_IS_COST_MATRIX,
} mfl_optimisation_t;


typedef long double *mfl_costs_t;       // A matrix of transition costs.


typedef struct mfl_datablock_t {
    int db_n_characters;                        // The number of characters within the datablock.
    int db_n_taxa;                              // The number of taxa the datablock applies to.
    mfl_optimisation_t db_optimisation_method;  // The optimisation method applied to all characters in this datablock.
    bool db_inapplicables;                      // 0 indicates no inapplicables; 1 means block has inapplicables.
    mfl_costs_t *db_costmatrix;                 // Cost matrix associated with these characters.
    charstate *db_characters;                   // The characters to which the datablock applies.
} mfl_datablock_t;


typedef mfl_datablock_t *mfl_chardata_t;    // Pointer to a set of different datablocks.

typedef struct mfl_searchrec_t {
    long int sr_num_reps;
    long int sr_best_length;
} mfl_searchrec_t;

/*
 *
 * Prototypes for all functions in the Morphy Function Library
 *
 */

/* In mfl_characters.c */

/* In mfl_starttree.c */

/* In mfl_tree.c */
mfl_node_t*     mfl_make_new_n_ary_ring_node(mfl_nodearray_t nodearray, mfl_node_t *bottom_node, int num_branches, int num_nodes);
bool            mfl_check_is_in_ring(mfl_node_t *start);
mfl_node_t *    mfl_insert_node_in_ring(mfl_node_t *ring_start, mfl_node_t *new_node);
mfl_node_t *    mfl_get_next_available_node(mfl_nodearray_t nodearray, int num_nodes);
bool            mfl_node_is_available(mfl_node_t *node);
void            mfl_disconnect_node_edges(mfl_node_t *node1, mfl_node_t *node2);
void            mfl_join_node_edges(mfl_node_t *node1, mfl_node_t *node2);
bool            mfl_check_node_is_bottom(mfl_node_t *querynode);
void            mfl_join_nodes(mfl_node_t *node1, mfl_node_t *node2);
mfl_node_t*     mfl_remove_branch(mfl_node_t *free_node_bottom, mfl_node_t *free_node_top);
void            mfl_insert_branch(mfl_node_t *src_bottom_node, mfl_node_t *src_free_desendant_edge, mfl_node_t *tgt_branch_bottom);
void            mfl_make_ring(mfl_node_t *bottom_node, mfl_node_t *left_node, mfl_node_t *right_node);
int             mfl_calculate_number_of_nodes_to_allocate(int num_taxa);
int             mfl_node_is_in_binary_ring(mfl_node_t *test_node);
void            mfl_set_internal_nodes_to_rings(mfl_nodearray_t nodearray, int num_taxa, int num_nodes);
mfl_node_t*     mfl_alloc_node(void);
void            mfl_free_node(mfl_node_t *node);
void            mfl_free_treenodes(mfl_nodearray_t treenodes, int num_nodes);
void            mfl_allocate_nodes_in_array(mfl_nodearray_t nodearray, int num_nodes, int num_taxa);
mfl_nodearray_t mfl_allocate_nodearray(int num_taxa, int num_nodes);
void            mfl_free_nodearray(mfl_nodearray_t nodearray);
struct          mfl_tree_t * mfl_alloctree_with_nodes(int num_taxa);
void            mfl_free_tree(mfl_tree_t *tree_to_free, int num_taxa, int num_nodes);

/* In mfl_initialise.c */

/* In mfl_starttree.c */

/* In mfl_newick.c */
int     mfl_check_invalid_newick(char *newick_input);
bool    mfl_newick_string_is_rooted(char *newick_string);
