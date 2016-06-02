//
//  tuimfy.h
//  Morphy
//
//  Created by mbrazeau on 13/04/2016.
//
//

#ifndef tuimfy_h
#define tuimfy_h

#include "ncl/ncl.h"

#define TUI_SILENT 0
#define TUI_VERBOSE 1

typedef struct tui_testrec {
  
    unsigned long long int tr_num_err;   // Number of errors counted
    unsigned long long int tr_err_overflow; // If the number of errors overflows ULLONG_MAX
    unsigned long long int tr_num_warn;         // As above but for warnings
    unsigned long long int tr_warn_overflow;
    int tr_verbosity;
    
} tui_testrec;


typedef enum tui_testcmd_t {
    TUI_CMD_MATRIX_INPUT_BASIC,
    TUI_CMD_NEWICK_INPUT_BASIC,
    TUI_CMD_COMMAND_BASIC,
    
    TUI_CMD_MAX
} tui_testcmd_t;


/*
 * Function prototypes for testing
 */
void tui_print_out_converted_matrix(mfl_matrix_t *matrix, int num_taxa, int num_chars, bool printbits);
void tui_test_tree_printing();
void tui_spr_test_environment(void);
void tui_test_addition_sequence(void);
int tui_build_destroy_tree_from_binary_newick(mfl_handle_t *mfl_handle, char *input_newick);
int tui_getting_numstates_test(void);
void tui_test_character_stuff(void);
void tui_test_bipartition_setting(void);
void tui_nexus_reader(char* argv1);
void tui_test_newick_stuff(void);
void tui_test_tree_copying(void);
void tui_test_nary_ring_creation(void);

/*  tui_utilities.c
 */
int             tui_check_all_node_ring_circularity(const mfl_tree_t *t, int num_nodes, int *verbose);
int             tui_check_for_anastomosis(mfl_tree_t* t, int num_nodes, int *verbose);
int             tui_check_tree_for_connection_errors(mfl_tree_t* t, int num_nodes, int *verbose);
bool            tui_check_reciprocal_edge(mfl_node_t *n, mfl_tree_t* t, int num_nodes, int *verbose);
int             tui_count_num_taxa(mfl_tree_t *t);
void            mfl_tree_check_traversal(mfl_node_t *node);
void            tui_print_node_data(mfl_node_t* p, const char* calling_fxn);
int             tui_tip_check(mfl_node_t* n, const char* calling_fxn, const int* verbose);
void*           tui_check_binary_traversal(mfl_node_t *p, const int* verbose, const char* calling_fxn);
int             tui_check_node_is_root(mfl_node_t* p, const int* verbose, const char* calling_fxn);
int             tui_check_all_binary(mfl_tree_t *querytree, const int *verbose);
int             tui_check_broken_tree(mfl_tree_t *t, int *verbose);
void            tui_print_newick_recursive(mfl_node_t *start);
void            tui_print_newick(mfl_node_t *n);
void            tui_print_edgetable(mfl_edgetable_t* edgetable);

/*  tui_io.c
 */
void            tui_ignore_nexus_comment(char **current);
void            tui_print_node_bipartition(mfl_node_t* n);
void            tui_partition_print_traversal(mfl_node_t* n);
void            tui_print_charstate_bits(const mfl_charstate cell, const int max_states);
char*           tui_readfile_to_str(FILE *input);
void            tui_destroy_file_string(char* oldinput);
void            tui_parse_test_infile(char *infile);
int             tui_check_simple_table_formatted(const char* input_table);
mfl_handle_s*   tui_parse_test_file(const char* arg1, const char* arg2);

/*  tmtarixproc.c
 */
void            tui_simple_table_parser(const char* input_table, mfl_handle_s* test_handle);
int             tui_check_simple_table_dimensions(const char* table, int rows, int cols);
void            tui_get_simple_table_dimensions(const char*table, int* rows, int* cols);
char*           tui_get_simple_table_matrix(const char* input_table);
int tui_test_matrix_processing(mfl_handle_s *mfl_handle);

/* to clean up*/
void            mfl_test_newick_stuff();
void            tui_test_character_stuff();
void            tui_test_checktree_(void);

#endif /* tuimfy_h */

