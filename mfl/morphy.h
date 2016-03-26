/*
 *
 * THE MORPHY FUNCTION LIBRARY (MFL)
 *
 */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include "mfl.h"

#ifdef MFY_DEBUG
#include <stdio.h>
#define dbg_printf(...) printf(__VA_ARGS__)
#else
#define dbg_printf(...)
#endif


/* 
 *
 * Definitions of some default values
 *
 */

#define MORPHY_DEFAULT_TREE_LIMIT 200

/*
 *
 * Data type definitions
 *
 */


typedef uint64_t charstate; // Each character state is represented by a single unsigned 64-bit integer. Thus, one character may have 64 possible states.

typedef struct {
    long int n_rearrangements; // Number of tree topologies visited
    int n_savetrees;    // Number of trees found and/or saved
    int bestlength;     // Length of the best tree
    time_t searcht;     // Time in search
    char** newicktrees; // Vector of trees in newick format
    string **texttrees;   // Array of trees drawn as text
    int num_islands;    // Number of parsimony islands hit
    int *island_lengths;       // Lengths of each island
    int *island_sizes;
    int *times_hit_islands;      // List of times each island was hit
#ifdef VERSION_1_5
    double consist_ind; // Consistency index
    double ret_ind;     // Retention index
    double resc_ci;     // Rescaled consistency index
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


typedef struct mfl_node_t {
	struct mfl_node_t *nodet_edge, *nodet_next; // Pointers to the neighboring node in the tree; the next node in the node ring.
	char *nodet_tipname;                        // Name of the tip from the dataset.
	int nodet_tip;                              // 1-based identifier of terminal. Assigned 0 if node is internal.
	int nodet_index;                            // 0-based index of node in the node-array. In rings, this should be identical for all nodes.
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
void            mfl_setup_nodearray(mfl_nodearray_t nodearray, int num_nodes, int num_taxa);
void            mfl_initialise_ring_node(mfl_node_t *bottom_node);
mfl_node_t*     mfl_make_new_n_ary_ring_node(mfl_node_t *bottom_node, int num_branches);
void            mfl_create_binary_fork(mfl_node_t *parent, mfl_node_t *child1, mfl_node_t *child2);
void            mfl_destroy_n_nary_ring(mfl_node_t *bottom_node);
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

/* In mfl_brwap.c */
bool    mfl_heuristic_search(mfl_handle_s *mfl_handle);

/* in mfyinterface.c*/
void                    mfl_free_input_data(mfl_handle_s *mfl_struct);
bool                    mfl_set_ntax(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_set_nchar(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_set_searchtype(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_set_numiterations(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_set_treelimit(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_set_branchswap_t(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_set_ratchet_status(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_attach_inputdata(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_set_addseq_t(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_set_collapse(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_set_collapse_value(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_set_gapormissing(mfl_handle_s *mfl_struct, void *param_data);
bool                    mfl_set_parameter(mfl_handle_t *mfl_handle, mfl_param_t param_type, void *param_data);
mfl_resultant_data_s*   mfl_morphy_controller(mfl_handle_s *mfl_handle);
mfl_handle_t            mfl_s2t(mfl_handle_s *mfl_handle);
mfl_handle_s*           mfl_t2s(mfl_handle_t mfl_handle);
