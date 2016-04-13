//
//  tuimfy.h
//  Morphy
//
//  Created by mbrazeau on 13/04/2016.
//
//

#ifndef tuimfy_h
#define tuimfy_h


/*
 * Function prototypes for testing
 */

/*  tui_utilities.c
 */
void tui_print_node_data(mfl_node_t* p, const char* calling_fxn);
int tui_tip_check(mfl_node_t* n, const char* calling_fxn, const int* verbose);
void* tui_check_binary_traversal(mfl_node_t *p, const int* verbose, const char* calling_fxn);
int tui_check_node_is_root(mfl_node_t* p, const int* verbose, const char* calling_fxn);
int tui_check_all_binary(mfl_tree_t *querytree, const int *verbose);

/*  tui_io.c
 */
int tui_check_simple_table_formatted(const char* input_table);
void tui_simple_table_parser(const char* input_table, mfl_handle_t test_handle);
int tui_check_simple_table_dimensions(const char* table, int rows, int cols);
void tui_get_simple_table_dimensions(const char*table, int* rows, int* cols);

/*  tmtarixproc.c
 */
int tui_test_matrix_processing(void);

#endif /* tuimfy_h */

