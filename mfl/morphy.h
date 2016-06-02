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
#include <time.h>
#include <gsl/gsl_rng.h>
#include "mfl.h"

#ifdef MFY_DEBUG
#include <stdio.h>
#define dbg_printf(...) printf(__VA_ARGS__)
#define dbg_eprintf(...) printf("ERROR in %s(), %s\n\n", __FUNCTION__, __VA_ARGS__)
#define dbg_linerr(...) printf("ERROR at line: %i in %s. %s\n\n", __LINE__, __FILE__, __VA_ARGS__)
#define dbg_pfail(...) printf("== FAIL == %s(): %s\n", __FUNCTION__, __VA_ARGS__)
#define dbg_ppass(...) printf("== PASS == %s(): %s\n", __FUNCTION__, __VA_ARGS__)
#define dbg_pcall(...) printf("\t\tCalled by %s()\n\n", __VA_ARGS__)
#else
#define dbg_printf(...)
#define dbg_eprintf(...)
#define dbg_linerr(...)
#define dbg_pfail(...)
#define dbg_ppass(...)
#define dbg_pcall(...)
#endif

using namespace std;

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define __FXN_NAME__ (const char*)__FUNCTION__
#define mfl_malloc(size, memesetval) __MFL_MALLOC__(size, memesetval, __FXN_NAME__)

/* 
 *
 * Definitions of some default values
 *
 */
#define MORPHY_UINTMAX UINT64_MAX
#define MORPHY_INAPPLICABLE_BITPOS ((mfl_charstate)1)
#define MORPHY_IS_APPLICABLE (~MORPHY_INAPPLICABLE_BITPOS)
#define MORPHY_MISSING_DATA_BITWISE (~1)
#define MORPHY_VALID_NONALPHA_STATES char* __MORPHY_NONALPHAS = {'+','-','@'};
#define MORPHY_NUM_VALID_NONALPHA  3
#define MORPHY_SPECIAL_STATE_PAD 1 /* Bit width used to reserve a position for a special state*/
#define MORPHY_MAX_STATE_NUMBER (sizeof(mfl_charstate) - MORPHY_SPECIAL_STATE_PAD)
#define MORPHY_BTS_IN_BITSET (sizeof(mfl_bitfield_t) * CHAR_BIT)

//Defaults
#define MORPHY_DEFAULT_WEIGHT 1.0
#define MORPHY_DEFAULT_TREE_LIMIT 1000
#define MORPHY_DEFAULT_TREEBUFFER_AUTOINCREASE_SWITCH false
#define MORPHY_DEFAULT_TREEBUFFER_AUTOINCREASE_AMOUNT 500
#define MORPHY_DEFAULT_STEPWISE_HOLD 1
#define MORPHY_DEFAULT_ADDITION_SEQUENCE_REPS 50
#define MORPHY_DEFAULT_CHARACTER_INCLUDE true
#define MORPHY_DEFAULT_MULTISTATE_HANDLE MFL_MULTSTATE_UNCERTAINTY
#define MORPHY_DEFAULT_PARSIMONY_METHOD MFL_OPT_FITCH
#define MORPHY_DEFAULT_PARSIMONY_NAME "Fitch (unordered)"

/*
 *
 * Data type definitions
 *
 */

typedef uint64_t mfl_uint;
typedef mfl_uint mfl_charstate; // Each character state is represented by a single unsigned 64-bit integer. Thus, one character may have 64 possible states.

typedef mfl_uint mfl_bitfield_t;

typedef struct mfl_stepmatrix_t {
    int      sm_numstates;
    bool     sm_is_int;
    union {
        int*     sm_int_costs;
        double*  sm_flt_costs;
    };
} mfl_stepmatrix_t;

typedef struct mfl_joint_character_t {
    int jpt_num_allapplic;
    int jpt_applic_part_num;
    int jpt_num_winapplic;
    int jpt_inapplic_part_num;
} mfl_joint_character_t;

typedef struct mfl_bitset_t {
    int bts_nfields;
    int bts_max_bitfields;
    int bts_max_bit;
    mfl_bitfield_t* bts_bitfields;
} mfl_bitset_t;

typedef struct mfl_bitsetlist_t {
    int bsl_num_sets;
    int bsl_max_sets;
    mfl_bitset_t** bsl_bitsets;
} mfl_bitsetlist_t;

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
    bool*                   included_taxa;
    int                     n_chars;
    bool*                   included_chars;
    int                     n_outgroup_taxa;
    int*                    outgroup_taxa;
    mfl_inputformat_t       input_format;
    mfl_search_t            search_type;
    int                     n_iterations;
    int                     n_treelimit;
    bool                    autoincrease;
    int                     maxtrees;
    int                     autoinc_incr;
    int                     n_ctypes;
    char*                   ctypes_cmd[MFL_OPT_MAX];
    mfl_parsimony_t*        ctype_setters;
    int                     n_usertypes;
    mfl_stepmatrix_t*       usertypes;
    mfl_branch_swap_t       bswap_type;
    bool                    is_ratchet;
    int                     n_symbols;
    char*                   format_symbols;
    char*                   input_data;
    int                     n_input_newick_trees;
    char**                  input_newick_trees;
    int*                    int_weights;
    double*                 real_weights;
    int                     n_datatypes;
    mfl_datatype_t          datatypes[MFL_DATATYPE_MAX];
    mfl_add_sequence_t      addseq_type;
    int                     n_to_hold;
    int                     reference_taxon;
    bool                    collapse_nolen;
    mfl_set_collapse_at_t   collapse_at;
    mfl_gap_t               gap_method;
    mfl_resultant_data_s    *resultant_data;
    unsigned long int       rseed;
} mfl_handle_s;


typedef enum {
    MFL_MST_UNCERTAIN,
    MFL_MST_POLYMORPH,
    MFL_MST_VARIABLE,
    
    MFL_MST_MAX
} mfl_multistate_t;

//typedef long double *mfl_costs_t;       // A matrix of transition costs.


typedef void (*mfl_parsim_fn)(struct mfl_node_t* parent);       // Pointer for a function that performs parsimony calculations at a node.
typedef mfl_charstate (*mfl_char2bit_fn)(char *states, char* datype_converter, mfl_gap_t gaprule);    // Pointer to conversion functions following conversion rules for a particular character type


typedef struct mfl_datapartition_t {
    int part_index;
    int part_n_chars_max;
    int part_n_chars_included;
    mfl_parsimony_t part_optimisation_method;
    bool part_has_inapplicables;
    bool part_char_is_directed;
    int *part_char_indices;                         // The partitions can contain characters be non-sequentially and out of order, this allows them to be identified after a search.
    int* part_int_weights;
    double* part_real_weights;
    mfl_parsim_fn part_downpass_full;
    mfl_parsim_fn part_downpass_partial;
    mfl_parsim_fn part_uppass_full;
    mfl_parsim_fn part_uppass_partial;
    mfl_stepmatrix_t* part_stepmatrix;
    mfl_charstate *part_matrix;
} mfl_datapartition_t;


typedef struct {
    int ptset_n_parts;
    mfl_datapartition_t** ptset_partitions;
} mfl_partition_set_t;


typedef struct mfl_character_vector_t {
    long long int cv_col_number;
    int cv_num_gaps;
    mfl_parsimony_t cv_parsim_method;
    mfl_multicell_t cv_multistate_method;
    mfl_char2bit_fn cv_conversion_rule;
    int cv_num_states;
    int cv_partition_destination;
    char** cv_character_cells;
    mfl_charstate* cv_chardata;
} mfl_character_vector_t;


typedef struct mfl_matrix_t {
    int  mat_num_taxa;
    int  mat_num_characters;
    int  mat_max_states;
    bool* mat_included_characters;
    int* mat_int_weights;
    double* mat_real_weights;
    mfl_character_vector_t** mat_matrix;
} mfl_matrix_t;


typedef struct mfl_nodedata_t {
    int nd_n_characters;                        // The number of characters within the datablock.
    mfl_datapartition_t *nd_parent_partition;
    mfl_parsim_fn nd_downpass_full;                  // The downpass parsimony function.
    mfl_parsim_fn nd_downpass_partial;
    mfl_parsim_fn nd_uppass_full;                    // The uppass parsimony function.
    mfl_parsim_fn nd_uppass_partial;
    mfl_charstate *nd_prelim_set;               // The initial downpass set for the whole tree.
    mfl_charstate *nd_final_set;                // The final uppass set for the whole tree.
    mfl_charstate *nd_subtree_prelim_set;       // The initial downpass set of the subtree when the tree broken.
    mfl_charstate *nd_subtree_final_set;        // The final uppass set of the subtree when the tree is broken.
} mfl_nodedata_t;


typedef struct mfl_nodestack_t mfl_nodestack_t;


typedef struct mfl_node_t {
	mfl_node_t *nodet_edge, *nodet_next; // Pointers to the neighboring node in the tree; the next node in the node ring.
    mfl_nodestack_t* nodet_ndstack;           // The nodestack of the tree to which this node belongs.
	char *nodet_tipname;                        // Name of the tip from the dataset.
	int nodet_tip;                              // 1-based identifier of terminal. Assigned 0 if node is internal.
	int nodet_index;                            // 0-based index of node in the node-array. In rings, this should be identical for all nodes.
    mfl_bitset_t* nodet_bipart;
    int row;
    int col;
    int branchl_cdraw;
    bool nodet_isingroup;                       // Indicates if node is within the ingroup or not.
    int nodet_advancement_index;
    int nodet_weight;                           // The number of tips of this node; tips are 1.
    int nodet_isbottom;                         // Indicates node points to (calculation) root.
    int nodet_edge_ref;                         // Indicates the node number identifier for the edgetable (only when isbottom is true!).
    int nodet_downpass_visited;                 // Indicates successful visit from a downpass traversal.
	int nodet_uppass_visited;                   // Indicates successful visit from an uppass traversal.
    int nodet_minsteps;                         // Minimum number of transformations along branch represented by nodet_edge.
	int nodet_maxsteps;                         // Maximum number of transformations along branch represented by nodet_edge.
    double nodet_mean_branchlen;                // Mean branch length.
    long int nodet_tree_index;             // Identity of the tree to which this node belongs.
    int nodet_num_dat_partitions;
    mfl_nodedata_t **nodet_charstates;
} mfl_node_t;


typedef mfl_node_t ** mfl_nodearray_t;  // A pointer of pointers to nodes used in trees.


typedef struct {
    mfl_node_t* src1;
    mfl_node_t* src2;
    mfl_node_t* tgt1;
    mfl_node_t* tgt2;
} mfl_cliprec_t;

typedef struct {
    int ste_max_edges;
    int ste_num_edges;
    int ste_head;
    mfl_node_t** ste_edges;
} mfl_subtree_edges_t;


typedef struct mfl_nodestack_t {    // Manages memory for unused nodes and avoid allocation on the fly or searches in the array.
    int nstk_numnodes;
    int nstk_maxsize;
    mfl_nodearray_t nstk_availbale_nds;
} mfl_nodestack_t;


typedef struct mfl_tree_t {
	mfl_nodearray_t treet_treenodes;        // The nodes of the tree. Never point these to another instance of mfl_tree.
	mfl_node_t *treet_root;                 // Pointer to the root of the tree.
    mfl_node_t treet_dummynode;
    mfl_node_t *treet_start;                // Starting node for operations on unrooted tree.
    int treet_num_og_tips;
    mfl_nodearray_t treet_outgroup_tips;    // Pointers to the outgroup tips.
    mfl_nodestack_t* treet_nodestack;       // Unused nodes
    mfl_bitsetlist_t* treet_bipartitions;
    int treet_num_taxa;                     // Total number of terminals.
    int treet_num_nodes;
    int treet_uw_parsimonylength;           // Unweighted number of steps under parsimony.
    int treet_island_id;                    // An identification number for the parent start tree.
    int *treet_compressed_tree;             // The integer encoding of the tree for tree comparisons and saving memory
	int *treet_comptree_holder;             // Holder for comparisons.
	char *treet_newick_format;              // The tree stored as a Newick format string.
	int treet_index;                        // Identifier for the tree in the saved trees array.
} mfl_tree_t;


typedef struct mfl_treebuffer_t {
    int tb_num_trees;
    int tb_max_buffersize;
    int tb_maxtrees;
    int tb_max_steps;
    int tb_bestlength;
    int tb_num_islands;
    mfl_tree_t** tb_savedtrees;
} mfl_treebuffer_t;

/*!
 @discussion Contains the parameters of the tree search process.
 */
typedef struct mfl_searchrec_t {
    bool sr_abort_rep;
    bool sr_abort_swapping;
    int  sr_num_taxa_included;
    int* sr_included_taxa;
    int  sr_num_chars_included;
    int* sr_included_chars;
    int  sr_num_partitions;
    mfl_search_t sr_searchtype;
    mfl_branch_swap_t sr_bswaptype;
    mfl_add_sequence_t sr_stepwise;
    long int sr_num_reps_stepwise;
    mfl_uint sr_best_length;
    int  sr_num_trees_held_stepwise;
    bool sr_increase_treebuffer;
    int  sr_autoinc_increment;
    int  sr_maxtrees;
    long long int sr_rearrangement_counter;
    mfl_node_t* sr_swap_entry;
    mfl_handle_s* sr_handle_ptr;
    mfl_tree_t* sr_swaping_on;
    mfl_treebuffer_t* sr_treebuffer;
    gsl_rng* sr_random_number;
} mfl_searchrec_t;


typedef struct mfl_island_data_t {
    int isd_index;              // Each new, unique starting tree gets a new number.
    int isd_length;             // Length of tree(s) that define the island.
    int isd_size;               // Number of trees in the island.
    mfl_tree_t **isd_members;   // The list of tree members stored as mfl_tree_t's.
} mfl_island_data_t;


typedef struct {
    int numentries;             // The number of entries in the edge table
    int num_tips;               // The maximum number of tips
    int num_nodes;              // The maximum number of nodes
    int* edgetable;             // The edge table
    mfl_node_t** edge_tips;     // Pointer to the tips addresses
    mfl_node_t** edge_nodes;    // Pointer to the nodes addresses
} mfl_edgetable_t;



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
bool            mfl_is_valid_morphy_statesymbol(char c);
void            mfl_populate_chartype_character_vectors(mfl_matrix_t *matrix, char *input_data_matrix, int num_chars, int num_taxa);
void            mfl_copy_multistate_subtoken_to_substring(char** xstatetoken, char* substring);
void            mfl_copy_singleton_subtoken_to_substring(char singleton, char* substring);
void            mfl_move_in_nexus_multistate(char **col);
mfl_charstate*  mfl_allocate_nodal_character_set(int num_characters);
int             mfl_convert_nexus_symbol_to_shift_value(char in, char *datype_converter);
int             mfl_shift_value_DEFAULT_catdata(char in);
bool            mfl_char_is_DNA_base(char in);
int             mfl_shift_value_DNA_base(char in);
mfl_charstate   mfl_convert_nexus_multistate(char *xstates, char *datype_converter);
mfl_charstate   mfl_convert_using_default_rule(char *xstates);
mfl_charstate   mfl_convert_gap_character(mfl_gap_t gapmethod);
void            mfl_set_datatype_converter_from_symbol_list(char* datype_converter, char* symbols_list);
char*           mfl_search_char_in_chartypes_array(char* key, char* list, int *listend, int listmax);
mfl_parsimony_t*    mfl_alloc_chartype_list(int input_num_chars);
void            mfl_destroy_chartype_list(mfl_parsimony_t *ctype_list);
mfl_charstate   mfl_standard_conversion_rule(char *c, char* datype_converter, mfl_gap_t gaprule);
mfl_parsimony_t*    mfl_get_chartypes_list(const mfl_handle_s* mfl_handle);
mfl_char2bit_fn mfl_fetch_conversion_rule(mfl_parsimony_t parsimtype);
void            mfl_apply_conversion_rules(mfl_matrix_t *matrix);
void            mfl_convert_charcells_to_mfl_charstates(mfl_character_vector_t* cv, const mfl_handle_s* handle);
void            mfl_convert_all_characters_to_charstates(mfl_matrix_t* m, const mfl_handle_s* handle);
bool            mfl_check_nexus_matrix_dimensions(char *input_matrix, int input_num_taxa, int input_num_chars);
int             mfl_check_state_number_support(char *datatype_list);
void            mfl_destroy_character_cells(char **char_cells, int num_states, int num_taxa);
void            mfl_destroy_mfl_matrix(mfl_matrix_t *oldmatrix, int num_taxa, int num_chars);
mfl_matrix_t*   mfl_create_mfl_matrix(int num_taxa, int num_chars);
void            mfl_setup_new_empty_matrix(mfl_matrix_t *newmatrix, int num_states, int num_taxa, int num_chars);
void            mfl_free_nodedata(mfl_nodedata_t *olddata);
mfl_nodedata_t* mfl_alloc_datapart(void);
int*            mfl_alloc_set_list(int num_chars);
void            mfl_free_set_list(bool *inclist);
void            mfl_read_nexus_style_list_subcmd(char *subcommand, int setval, int* list, int nelems);
void            mfl_count_gaps_in_each_character(mfl_matrix_t* matrix);
void            mfl_set_cv_chartypes(mfl_matrix_t* matrix, const mfl_handle_s* mfl_handle, mfl_parsimony_t* chartypes);
void            mfl_set_include_value(int vectornum, int includeval, int* includes);
void            mfl_set_include_range(int first, int last, int includeval, int* includes);
void            mfl_set_inclusion_list(int* includes, int includeval, int listmax, char *subcommand);
void            mfl_move_current_to_digit(char** current);
int             mfl_get_numstates_from_matrix(char *inputmatrix);
int             mfl_read_nexus_type_int(char **current);
void            mfl_skip_spaces(char **current);
bool            mfl_is_nexus_stop_position(char a);
mfl_datapartition_t* mfl_alloc_empty_datapartition_t(void);
int             mfl_count_num_partitions_required(struct mfl_joint_character_t* jntchars, mfl_matrix_t* m, mfl_gap_t gaprule);
void            mfl_set_datapart_params_fitch(mfl_datapartition_t* d, bool has_inapplic);
void            mfl_set_datapart_params_wagner(mfl_datapartition_t* d, bool has_inapplic);
void            mfl_set_datapart_params_dollo(mfl_datapartition_t* d, bool has_inapplic, bool up);
void            mfl_set_datapart_params_irrev(mfl_datapartition_t* d, bool has_inapplic, bool up);
void            mfl_set_datapart_params_costmatrx(mfl_datapartition_t* d, mfl_handle_s* handle, bool has_inapplic);
void            mfl_set_datapart_params(mfl_datapartition_t* d, mfl_parsimony_t opt_t, bool has_inapplic, mfl_handle_s* handle);
void            mfl_populate_all_character_partitions(mfl_partition_set_t* ptset, mfl_matrix_t* m);
void            mfl_copy_column_into_partition(mfl_datapartition_t* prt, mfl_character_vector_t* cv, int intwt, int index, int num_rows);
int             mfl_compare_dataparts_by_ctype(const void* p1, const void* p2);
int             mfl_compare_dataparts_by_index(const void* p1, const void* p2);
mfl_partition_set_t* mfl_create_data_partitions_set(mfl_matrix_t* matrix, mfl_handle_s* handle);
void            mfl_destroy_partition_set(mfl_partition_set_t* ptset);
mfl_matrix_t*   mfl_create_internal_data_matrix(const mfl_handle_s* mfl_handle);
mfl_partition_set_t* mfl_create_data_partitions_set(mfl_matrix_t* matrix, mfl_handle_s* handle);

/* In mfl_starttree.c */

/* In mfl_tree.c */
void            mfl_initialise_ring_node(mfl_node_t *bottom_node);
void            mfl_create_binary_fork(mfl_node_t *parent, mfl_node_t *child1, mfl_node_t *child2, mfl_nodearray_t nodes);
mfl_node_t*     mfl_make_new_n_ary_ring_node(int num_branches, mfl_nodestack_t* ndstk);
void            mfl_dissolve_n_ary_ring(mfl_node_t* n, mfl_nodestack_t* ndstk);
void            mfl_destroy_n_nary_ring(mfl_node_t *bottom_node);
bool            mfl_check_is_in_ring(mfl_node_t *start);
mfl_node_t*     mfl_insert_node_in_ring(mfl_node_t *ring_start, mfl_node_t *new_node);
mfl_node_t*     mfl_get_node_from_nodestack(mfl_nodestack_t *nds);
void            mfl_destroy_nodestack(mfl_nodestack_t* ndstk);
mfl_nodestack_t* mfl_create_empty_nodestack(int num_internal_nodes);
void            mfl_push_node_to_nodestack(mfl_node_t* n, mfl_nodestack_t* ndstk);
mfl_node_t*     mfl_get_next_available_node(mfl_nodearray_t nodearray);
bool            mfl_node_is_available(mfl_node_t *node);
void            mfl_disconnect_node_edges(mfl_node_t *node1, mfl_node_t *node2);
void            mfl_disconnect_node(mfl_node_t* n);
void            mfl_join_node_edges(mfl_node_t *node1, mfl_node_t *node2);
void            mfl_make_node_available(mfl_node_t *node);
void            mfl_safe_reset_node_params(mfl_node_t* node);
bool            mfl_check_node_is_bottom(mfl_node_t *querynode);
mfl_node_t*     mfl_remove_branch(mfl_node_t *free_node_bottom, mfl_node_t *free_node_top);
void            mfl_insert_branch_with_ring_base(mfl_node_t *src, mfl_node_t *tgt);
void            mfl_make_ring(mfl_node_t *bottom_node, mfl_node_t *left_node, mfl_node_t *right_node);
int             mfl_calculate_number_of_nodes_to_allocate(int num_taxa);
mfl_node_t*     mfl_alloc_node(void);
void            mfl_free_node(mfl_node_t *node);
void            mfl_free_treenodes(mfl_nodearray_t treenodes);
void            mfl_allocate_nodes_in_array(mfl_nodearray_t nodearray, int num_nodes, int num_taxa);
mfl_nodearray_t mfl_allocate_nodearray(int num_taxa, int num_nodes);
void            mfl_initialise_nodearray(mfl_tree_t* t, int num_taxa, int num_nodes);
void            mfl_free_nodearray(mfl_nodearray_t nodearray);
int             mfl_node_is_n_ary(mfl_node_t *querynode, int test_n_branches);
mfl_node_t*     mfl_find_rightmost_tip_in_tree(mfl_node_t* n);
void            mfl_unroot_tree(mfl_tree_t *tree);
void            mfl_initialise_tree(mfl_tree_t *newtree, int num_taxa, int num_nodes);
mfl_tree_t*     mfl_alloctree_with_nodes(int num_taxa);
void            mfl_free_tree(mfl_tree_t *tree_to_free);
void            mfl_root_target_node(mfl_tree_t *input_tree, mfl_node_t *target_node_ring_start); //TG: WARNING: this creates a polytomy on the target node.
void            mfl_root_target_edge(mfl_tree_t *input_tree, mfl_node_t *target_node);
mfl_treebuffer_t* mfl_alloc_treebuffer(int num_trees);
void            mfl_append_tree_to_treebuffer(mfl_tree_t* newtree, mfl_treebuffer_t* trbuf, mfl_handle_s* mfl_handle);
void            mfl_resize_treebuffer(mfl_treebuffer_t* trbuf, int addedlength);
void            mfl_reset_nodestack(mfl_nodestack_t* nstk);
mfl_tree_t*     mfl_copy_tree_topology(const mfl_tree_t* t);
void            mfl_destroy_treebuffer(mfl_treebuffer_t* oldtreebuf, bool cleartrees);
void            mfl_assign_bottom_node(mfl_node_t* n);

//

/* In mfl_initialise.c */

/* In mfl_starttree.c */

mfl_treebuffer_t* mfl_get_start_trees(mfl_partition_set_t* dataparts, mfl_handle_s* handle, mfl_searchrec_t* searchrec);

/* In mfl_drawtree.c*/
char*   mfl_drawtree_create_virtual_grid(int num_taxa);
void    mfl_put_character_in_cell(char const ch, int row, int col, char* grid);
char    mfl_drawtree_get_character_in_cell(int row, int col, char* grid);
void    mfl_drawtree_write_into_tipfield(char* name, char* grid, int row, int col);
int     mfl_drawtree_calculate_branchlength(int num_taxa);
void    mfl_drawtree_set_node_coords(mfl_node_t *n, int row, int num_taxa);
void    mfl_drawtree_apply_subbranch(mfl_node_t* parent, mfl_node_t* desc, char *grid);
void    mfl_drawtree_set_coords_traversal(mfl_node_t* n, int* currentrow, char* grid, int num_taxa);
void    mfl_drawtree_add_nodebar(mfl_node_t* n, mfl_node_t* ldesc, mfl_node_t* rdesc, char* grid);
void    mfl_drawtree_draw_traversal(mfl_node_t* n, char *grid);
char*   mfl_drawtree(mfl_tree_t* t);

/* In mfl_newick.c */

// TODO: Update function prototypes as a bunch of these have been changed to have no return values
int         mfl_is_valid_newick(char *newick_input);
int         mfl_count_internal_nodes_in_newick(char *newick_string);
bool        mfl_newick_tree_is_rooted(char *newick_string);
int         mfl_read_newick_int(char **newick_position);
char*       mfl_find_next_opening_bracket_in_newick(char *newick_tree);
int         mfl_seek_largest_tip_number_newick(char *newick_string);
mfl_node_t* mfl_traverse_newick_recursively(char **newick_position, mfl_nodearray_t nodearray, int num_taxa);
mfl_tree_t* mfl_convert_newick_to_mfl_tree_t(char *newick_tree, int num_taxa);
int         mfl_number_of_digits_in_integer(int n);
int         mfl_traverse_tree_to_get_tip_char_length(mfl_node_t *start, int &tips_length);
int         mfl_traverse_mfl_tree_t_number_of_taxa(mfl_node_t *start, int &num_taxa);
int         mfl_number_of_characters_in_newick(int num_taxa, mfl_node_t *start);
char*       mfl_traverse_tree_to_print_newick_char_recursive(mfl_node_t *start, char *newick_tree_out, int &count);
char*       mfl_convert_mfl_tree_t_to_newick(mfl_tree_t *input_tree, bool root_polytomy);

void mfl_traverse_tree_to_get_tip_char_length(mfl_node_t *start, int *tips_length);
void mfl_traverse_mfl_tree_t_number_of_taxa(mfl_node_t *start, int* num_taxa);
void mfl_traverse_tree_to_print_newick_char_recursive(mfl_node_t *start, char *newick_tree_out, int* count);
char* mfl_alloc_empty_newick(int newick_string_length);
int mfl_get_num_active_taxa_in_tree(mfl_tree_t* input_tree);
int mfl_calculate_newick_length(mfl_tree_t* input_tree, int num_taxa);
void mfl_concatenate_newick_elements(char* destination, char* newick_substring, char* root_header);
char* mfl_get_newick_root_header(bool isrooted);
mfl_treebuffer_t *mfl_tree_from_newick_to_buffer(mfl_handle_s* mfl_handle);


/* In mfl_brwap.c */
inline void mfl_temp_rebranching(mfl_node_t* src, mfl_node_t* tgt, mfl_cliprec_t* regraft);
inline void mfl_undo_temp_rebranching(mfl_cliprec_t* regraft);
void mfl_regrafting_traversal(mfl_node_t* n, mfl_node_t* src, mfl_searchrec_t* searchrec, int startdistance, int trav);
void mfl_regraft_subtree(mfl_node_t* src, mfl_node_t* tgt, mfl_searchrec_t* searchrec, bool neighbor_rule);
bool    mfl_heuristic_search(mfl_handle_s *mfl_handle);
//inline mfl_node_t* mfl_clip_branch(mfl_node_t* n, mfl_cliprec_t* cliprec);
//inline void mfl_restore_branching(mfl_cliprec_t* cliprec);
void mfl_pruning_traversal(mfl_node_t* n, mfl_searchrec_t* searchrec);


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

/* in mfl_bitset.c*/
void            mfl_bts_setbit(mfl_bitset_t* bitset, mfl_uint set_to, int setposition);
bool            mfl_bts_AND(mfl_bitset_t* set1, mfl_bitset_t* set2, mfl_bitset_t* target);
bool            mfl_bts_OR(mfl_bitset_t* set1, mfl_bitset_t* set2, mfl_bitset_t* target);
bool            mfl_bts_XOR(mfl_bitset_t* set1, mfl_bitset_t* set2, mfl_bitset_t* target);
bool            mfl_bts_COMPLEMENT(mfl_bitset_t* bitset, mfl_bitset_t* target);
int             mfl_bts_calculate_n_bitfieds(int n_minbits);
mfl_bitset_t*   mfl_bts_create_bitset(int n_minbits);
bool            mfl_bts_destroy_bitset(mfl_bitset_t* oldbts);

/* in mfl_compare.c*/
void            mfl_set_bipartitions(mfl_node_t* n);
mfl_edgetable_t* mfl_initiate_edgetable_t(int num_tips, bool is_rooted);
void            mfl_get_edgetable_tips(mfl_edgetable_t* edgetable, mfl_tree_t* tree);
void            mfl_destroy_edgetable(mfl_edgetable_t* edgetable);
mfl_node_t*     mfl_find_bottom_node_in_ring(mfl_node_t* node);
void            mfl_set_edge_ref_in_ring(mfl_node_t* node, int reference);
int             mfl_get_edge_ref_from_ring(mfl_node_t* node);
void            mfl_get_edgetable(mfl_edgetable_t* edgetable, mfl_tree_t* tree);

bool            mfl_compare_edge_tables(mfl_edgetable_t* t1, mfl_edgetable_t* t2);

void            tui_print_edgetable(mfl_edgetable_t* edgetable); //TODO: move this function to tui
void            tui_test_edgetables(void); //TODO: move this function to tui

/* in mfl_searchrec.c*/
void                mfl_initialise_searchrec(mfl_searchrec_t* searchrec, const mfl_handle_s* handle);
mfl_searchrec_t*    mfl_create_searchrec(mfl_handle_s* handle);
/* in mfl_morphy.c*/
void*           __MFL_MALLOC__(size_t size, int memsetval, const char* fn_name);
