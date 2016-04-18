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

/*
 * Function prototypes for testing
 */

void tui_test_tree_printing();

/*  tui_utilities.c
 */
int tui_check_all_node_ring_circularity(const mfl_tree_t *t, int num_nodes, int *verbose);
int tui_check_for_anastomosis(mfl_tree_t* t, int num_nodes, int *verbose);
int tui_check_tree_for_connection_errors(mfl_tree_t* t, int num_nodes, int *verbose);
bool tui_check_reciprocal_edge(mfl_node_t *n, mfl_tree_t* t, int num_nodes, int *verbose);
int tui_count_num_taxa(mfl_tree_t *t);
void mfl_tree_check_traversal(mfl_node_t *node);
void tui_print_node_data(mfl_node_t* p, const char* calling_fxn);
int tui_tip_check(mfl_node_t* n, const char* calling_fxn, const int* verbose);
void* tui_check_binary_traversal(mfl_node_t *p, const int* verbose, const char* calling_fxn);
int tui_check_node_is_root(mfl_node_t* p, const int* verbose, const char* calling_fxn);
int tui_check_all_binary(mfl_tree_t *querytree, const int *verbose);
int tui_check_broken_tree(mfl_tree_t *t, int *verbose);
void tui_print_newick_recursive(mfl_node_t *start);
void tui_print_newick(mfl_node_t *n);

/*  tui_io.c
 */
int tui_check_simple_table_formatted(const char* input_table);
void tui_simple_table_parser(const char* input_table, mfl_handle_t test_handle);
int tui_check_simple_table_dimensions(const char* table, int rows, int cols);
void tui_get_simple_table_dimensions(const char*table, int* rows, int* cols);
int tui_parse_test_file(const char* arg1, const char* arg2);

/*  tmtarixproc.c
 */
int tui_test_matrix_processing(void);

/* to clean up*/
void mfl_test_newick_stuff();
void tui_test_character_stuff();
void tui_test_checktree_(void);

#endif /* tuimfy_h */

