/*
 *  morphy.h
 *
 *  THE MORPHY FUNCTION LIBRARY
 *  A library for phylogenetic analysis with emphasis on parsimony and
 *  morphology (but someday other methods)
 *
 *  Copyright (C) 2016  by Martin D. Brazeau, Thomas Guillerme,
 *  and Chris Desjardins
 *
 *  Some data structs, routines and ideas derived from:
 *      - The PHYLIP package by Joe Felsenstein
 *          <http://evolution.genetics.washington.edu/phylip.html>
 *      - MrBayes by John Huelsenbeck and Fredrik Ronquist
 *          <http://mrbayes.sourceforge.net/>
 *
 *  Any bugs, errors, inefficiences and general amateurish handling are our own
 *  and most likely the responsibility of MDB. We make no guarantees.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <ctype.h>
#include "mfl.h"

#ifdef MFY_DEBUG
#include <stdio.h>
#define dbg_printf(...) printf(__VA_ARGS__)
#define dbg_eprintf(...) printf("ERROR in %s(), %s\n\n", __FUNCTION__, __VA_ARGS__)
#define dbg_linerr(...) printf("ERROR at line: %i in %s. %s\n\n", __LINE__, __FILE__, __VA_ARGS__)
#define dbg_pfail(...) printf("== FAIL == %s(): %s\n", __FUNCTION__, __VA_ARGS__)
#define dbg_ppass(...) printf("== PASS == %s(): %s\n", __FUNCTION__, __VA_ARGS__)
#define dbg_pcall(...) printf("\t\tCalled by %s()\n\n", __VA_ARGS__)
#define __FXN_NAME__ (const char*)__FUNCTION__
#else
#define dbg_printf(...)
#define dbg_eprintf(...)
#define dbg_linerr(...)
#define dbg_pfail(...)
#define dbg_ppass(...)
#define dbg_pcall(...)
#endif

using namespace std;

/* 
 *
 * Definitions of some default values
 *
 */

#define MORPHY_SPECIAL_STATE_PAD 1 /* Bit width used to reserve a position for a special state*/
#define MORPHY_MAX_STATE_NUMBER (64 - MORPHY_SPECIAL_STATE_PAD) /* The reason this isn't 64 is because the
                                      the first bit position is reserved for 
                                      gap as a state or as logical impossibility
                                                                              */

//Defaults
#define MORPHY_DEFAULT_WEIGHT 1.0
#define MORPHY_DEFAULT_TREE_LIMIT 1000
#define MORPHY_DEFAULT_TREEBUFFER_AUTOINCREASE_DEFAULT 500
#define MORPHY_DEFAULT_STEPWISE_HOLD 1
#define MORPHY_DEFAULT_ADDITION_SEQUENCE_REPS 10
#define MORPHY_DEFAULT_CHARACTER_INCLUDE true
#define MORPHY_DEFAULT_MULTISTATE_HANDLE MFL_MULTSTATE_UNCERTAINTY

#define MORPHY_VALID_NONALPHA_STATES char* __MORPHY_NONALPHAS = {'+','-','@'};
#define MORPHY_NUM_VALID_NONALPHA  3


/*
 *
 * Data type definitions
 *
 */

typedef uint64_t mfl_charstate; // Each character state is represented by a single unsigned 64-bit integer. Thus, one character may have 64 possible states.

typedef struct {
    long int n_rearrangements;  // Number of tree topologies visited
    int n_savetrees;            // Number of trees found and/or saved
    int bestlength;             // Length of the best tree
    time_t searcht;             // Time in search
    char** newicktrees;         // Vector of trees in newick format
    string **texttrees;         // Array of trees drawn as text
    int num_islands;            // Number of parsimony islands hit
    int *island_lengths;        // Lengths of each island
    int *island_sizes;          // Number of trees in an island
    int *times_hit_islands;     // List of times each island was hit
#ifdef VERSION_1_5
    double consist_ind;         // Consistency index
    double ret_ind;             // Retention index
    double resc_ci;             // Rescaled consistency index
#endif
} mfl_resultant_data_s;


typedef struct {
    int                     n_taxa;
    int                     n_chars;
    mfl_inputformat_t       input_format;
    mfl_search_t            search_type;
    int                     n_iterations;
    int                     n_treelimit;
    int*                    fitch_list;
    int*                    wagner_list;
    int*                    dollo_list;
    int*                    irreversible_list;
    int*                    usertype_list;
    mfl_branch_swap_t       bswap_type;
    bool                    is_ratchet;
    char*                   format_symbols;
    char*                   input_data;
    mfl_add_sequence_t      addseq_type;
    bool                    collapse_nolen;
    mfl_set_collapse_at_t   collapse_at;
    mfl_gap_t               gap_method;
    mfl_resultant_data_s    *resultant_data;
} mfl_handle_s;


typedef enum {
    MFL_MST_UNCERTAIN,
    MFL_MST_POLYMORPH,
    MFL_MST_VARIABLE,
    
    MFL_MST_MAX
} mfl_multistate_t;

typedef long double *mfl_costs_t;       // A matrix of transition costs.


typedef void (*mfl_parsim_fn)(struct mfl_node_t* parent);       // Pointer for a function that performs parsimony calculations at a node.
typedef mfl_charstate (*mfl_char2bit_fn)(char *states);    // Pointer to conversion functions following conversion rules for a particular character type


typedef struct mfl_datapartition_t {
    int part_n_characters;
    mfl_optimisation_t part_optimisation_method;
    bool part_has_inapplicables;
    bool part_char_is_directed;
    int *part_char_indices;                         // The partitions can contain characters be non-sequentially and out of order, this allows them to be identified after a search.
    int part_weight;
    mfl_costs_t *part_costmatrix;
    mfl_charstate *part_matrix;
} mfl_datapartition_t;


typedef struct mfl_character_vector_t {
    int cv_num_colums;
    bool cv_has_gaps;
    mfl_multicell_t cv_multistate_method;
    mfl_char2bit_fn cv_conversion_rule;
    char** cv_character_cells;
    mfl_charstate* cv_chardata;
} mfl_character_vector_t;


typedef struct mfl_matrix_t {
    int mat_num_taxa;
    int mat_num_characters;
    int max_states;
    mfl_character_vector_t** mat_matrix;
} mfl_matrix_t;


typedef struct mfl_nodedata_t {
    int nd_n_characters;                        // The number of characters within the datablock.
    mfl_optimisation_t nd_optimisation_method;  // The optimisation method applied to all characters in this datablock.
    bool nd_has_inapplicables;                  // false: no inapplicables; true: has inapplicables.
    bool nd_char_is_directed;                   // Character depends on tree rooting or not.
    mfl_costs_t *nd_costmatrix;                 // Cost matrix associated with these characters.
    mfl_parsim_fn nd_downpass;                  // The downpass parsimony function.
    mfl_parsim_fn nd_uppass;                    // The uppass parsimony function.
    mfl_charstate *nd_prelim_set;               // The initial downpass set for the whole tree.
    mfl_charstate *nd_final_set;                // The final uppass set for the whole tree.
    mfl_charstate *nd_subtree_prelim_set;       // The initial downpass set of the subtree when the tree broken.
    mfl_charstate *nd_subtree_final_set;        // The final uppass set of the subtree when the tree is broken.
} mfl_nodedata_t;


typedef struct mfl_node_t {
	mfl_node_t *nodet_edge, *nodet_next; // Pointers to the neighboring node in the tree; the next node in the node ring.
	char *nodet_tipname;                        // Name of the tip from the dataset.
	int nodet_tip;                              // 1-based identifier of terminal. Assigned 0 if node is internal.
	int nodet_index;                            // 0-based index of node in the node-array. In rings, this should be identical for all nodes.
    bool nodet_isingroup;                       // Indicates if node is within the ingroup or not.
    int nodet_isbottom;                         // Indicates node points to (calculation) root.
    int nodet_downpass_visited;                 // Indicates successful visit from a downpass traversal.
	int nodet_uppass_visited;                   // Indicates successful visit from an uppass traversal.
    int nodet_minsteps;                         // Minimum number of transformations along branch represented by nodet_edge.
	int nodet_maxsteps;                         // Maximum number of transformations along branch represented by nodet_edge.
    double nodet_mean_branchlen;                // Mean branch length.
    int nodet_isroot;                           // Set non-zero if this node is the root.
    long long int nodet_tree_index;             // Identity of the tree to which this node belongs.
    int nodet_num_partitions;
    mfl_nodedata_t **nodet_dataparts;
} mfl_node_t;


typedef mfl_node_t ** mfl_nodearray_t;  // A pointer of pointers to nodes used in trees.


typedef struct mfl_tree_t {
	mfl_nodearray_t treet_treenodes;        // The nodes of the tree. Never point these to another instance of mfl_tree.
	mfl_node_t *treet_root;                 // Pointer to the root of the tree.
    mfl_node_t *treet_start;                // Starting node for operations on unrooted tree.
    mfl_nodearray_t treet_outgroup_tips;    // Pointers to the outgroup tips.
	int treet_uw_parsimonylength;           // Unweighted number of steps under parsimony.
    int treet_island_id;                    // An identification number for the parent start tree.
    int *treet_compressed_tree;             // The integer encoding of the tree for tree comparisons and saving memory
	int *treet_comptree_holder;             // Holder for comparisons.
	char *treet_newick_format;              // The tree stored as a Newick format string.
	int treet_index;                        // Identifier for the tree in the saved trees array.
} mfl_tree_t;


typedef struct mfl_taxon_partition_t {      // For outgroups or other clade constraints.
    int tp_num_taxa;
    mfl_nodearray_t tp_taxon_list;
} mfl_taxon_partition_t;


typedef struct mfl_treelist_t {
    int tl_maxtrees;
    int tl_maxlength;
    int tl_bestlength;
    int tl_num_islands;
    mfl_tree_t** tl_savedtrees;
} mfl_treelist_t;


typedef struct mfl_searchrec_t {
    long int sr_num_reps;
    long int sr_best_length;
    int      sr_num_partitions;
} mfl_searchrec_t;


typedef struct mfl_island_data_t {
    int isd_index;              // Each new, unique starting tree gets a new number.
    int isd_length;             // Length of tree(s) that define the island.
    int isd_size;               // Number of trees in the island.
    mfl_tree_t **isd_members;   // The list of tree members stored as mfl_tree_t's.
} mfl_island_data_t;




/*
 *
 * Prototypes for all functions in the Morphy Function Library
 *
 */

/* In mfl_characters.c */

//prob move these two to another file
void mfl_move_past_symbol(char** c, char symbol);
char *mfl_move_past_eq_sign(char *input);

bool            mfl_is_valid_morphy_ctype(char c);
void            mfl_populate_chartype_character_vector(mfl_matrix_t *matrix, char *input_data_matrix, int num_chars, int num_taxa);
void            mfl_copy_multistate_subtoken_to_substring(char** xstatetoken, char* substring);
void            mfl_copy_singleton_subtoken_to_substring(char singleton, char* substring);
void            mfl_move_in_nexus_multistate(char **col);
mfl_charstate*  mfl_allocate_nodal_character_set(int num_characters);
bool            mfl_check_nexus_matrix_dimensions(char *input_matrix, int input_num_taxa, int input_num_chars);
void            mfl_destroy_character_cells(char **char_cells, int num_states, int num_taxa);
void            mfl_destroy_mfl_matrix(mfl_matrix_t *oldmatrix, int num_states, int num_taxa, int num_chars);
mfl_matrix_t*   mfl_create_mfl_matrix(int num_taxa, int num_chars);
void            mfl_setup_new_empty_matrix(mfl_matrix_t *newmatrix, int num_states, int num_taxa, int num_chars);
void            mfl_free_nodedata(mfl_nodedata_t *olddata);
mfl_nodedata_t* mfl_alloc_datapart(void);
bool*           mfl_read_nexus_exset_subcmd(char *subcommand, int Nexus_NCHARS);
void            mfl_free_inclusion_list(bool *inclist);
bool*           mfl_alloc_character_inclusion_list(int num_chars);
void            mfl_set_inclusion_list(bool* includes, bool includeval, int listmax, char *subcommand);
void            mfl_move_current_to_digit(char** current);
void            mfl_set_include_range(int first, int last, bool includeval, bool* includes);
void            mfl_set_include_value(int vectornum, bool includeval, bool* includes);
int             mfl_get_numstates_from_matrix(char *inputmatrix);
int             mfl_read_nexus_type_int(char **current);
void            mfl_skip_spaces(char **current);
bool            mfl_is_nexus_stop_position(char a);
mfl_charstate   mfl_convert_nexus_multistate(char *xstates, char *datype_converter);

/* In mfl_starttree.c */

/* In mfl_tree.c */
void            mfl_setup_nodearray(mfl_nodearray_t nodearray, int num_nodes, int num_taxa);
void            mfl_initialise_ring_node(mfl_node_t *bottom_node);
void            mfl_create_binary_fork(mfl_node_t *parent, mfl_node_t *child1, mfl_node_t *child2, mfl_nodearray_t nodes);
mfl_node_t*     mfl_make_new_n_ary_ring_node(mfl_node_t *bottom_node, int num_branches, mfl_nodearray_t nodes);
void            mfl_destroy_n_nary_ring(mfl_node_t *bottom_node);
bool            mfl_check_is_in_ring(mfl_node_t *start);
mfl_node_t*     mfl_insert_node_in_ring(mfl_node_t *ring_start, mfl_node_t *new_node);
mfl_node_t*     mfl_get_next_available_node(mfl_nodearray_t nodearray);
bool            mfl_node_is_available(mfl_node_t *node);
void            mfl_disconnect_node_edges(mfl_node_t *node1, mfl_node_t *node2);
void            mfl_join_node_edges(mfl_node_t *node1, mfl_node_t *node2);
void            mfl_make_node_available(mfl_node_t *node);
void            mfl_safe_reset_node_params(mfl_node_t* node);
bool            mfl_check_node_is_bottom(mfl_node_t *querynode);
void            mfl_join_nodes(mfl_node_t *node1, mfl_node_t *node2);
mfl_node_t*     mfl_remove_branch(mfl_node_t *free_node_bottom, mfl_node_t *free_node_top);
void            mfl_insert_branch(mfl_node_t *src_bottom_node, mfl_node_t *src_free_desendant_edge, mfl_node_t *tgt_branch_bottom);
void            mfl_make_ring(mfl_node_t *bottom_node, mfl_node_t *left_node, mfl_node_t *right_node);
int             mfl_calculate_number_of_nodes_to_allocate(int num_taxa);
mfl_node_t*     mfl_alloc_node(void);
void            mfl_free_node(mfl_node_t *node);
void            mfl_free_treenodes(mfl_nodearray_t treenodes);
void            mfl_allocate_nodes_in_array(mfl_nodearray_t nodearray, int num_nodes, int num_taxa);
mfl_nodearray_t mfl_allocate_nodearray(int num_taxa, int num_nodes);
void            mfl_initialise_nodearray(mfl_nodearray_t nodearray, int num_taxa, int num_nodes);
void            mfl_free_nodearray(mfl_nodearray_t nodearray);
int             mfl_node_is_n_ary(mfl_node_t *querynode, int test_n_branches);
void            mfl_unroot_tree(mfl_tree_t *tree);
mfl_tree_t*     mfl_alloctree_with_nodes(int num_taxa);
void            mfl_free_tree(mfl_tree_t *tree_to_free);

/* In mfl_initialise.c */

/* In mfl_starttree.c */

/* In mfl_newick.c */
int         mfl_is_valid_newick(char *newick_input);
int         mfl_count_internal_nodes_in_newick(char *newick_string);
bool        mfl_newick_tree_is_rooted(char *newick_string);
int         mfl_read_newick_int(char **newick_position);
char*       mfl_find_next_opening_bracket_in_newick(char *newick_tree);
int         mfl_seek_largest_tip_number_newick(char *newick_string);
mfl_node_t* mfl_traverse_newick_recursively(char **newick_position, mfl_nodearray_t nodearray, int num_taxa);
mfl_tree_t* mfl_convert_newick_to_mfl_tree_t(char *newick_tree, int num_taxa);
int         mfl_number_of_digits_in_integer(int n);
int         mfl_power(int base, unsigned int exp);
int         mfl_number_of_digits_in_sequence(int n);
int         mfl_number_of_characters_in_newick(int num_taxa);
int         mfl_tree_t_number_of_taxa(mfl_tree_t *input_tree);
char*       mfl_convert_mfl_tree_t_to_newick(mfl_tree_t *root_start);

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
