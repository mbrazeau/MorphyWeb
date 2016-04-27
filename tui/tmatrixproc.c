/*
    tmatrixproc.c
    Morphy

 
    Tests for matrix processing
 */

#include "morphy.h"
#include "tuimfy.h"

int tui_check_simple_table_dimensions(const char* table, int rows, int cols)
{
    int expected = rows * cols;
    int actual = 0;
    char* c = (char*)table;
    
    mfl_move_past_symbol(&c, ';');
    
    do {
        if (mfl_is_valid_morphy_ctype(*c)) {
            ++actual;
        }
        else if (*c == '(' || *c == '{') {
            mfl_move_in_nexus_multistate(&c);
            //++c;
            ++actual;
        }
        ++c;
    } while (*c != ';');
    
    if (expected == actual) {
        return 0;
    }
    else {
        dbg_eprintf("input matrix has unexpected dimensions: ");
        dbg_printf("%i\n\n", actual);
        exit(actual);
    }
}


void tui_get_simple_table_dimensions(const char*table, int* rows, int* cols)
{
    char *c = (char*)table;
    
    tui_check_simple_table_formatted(table);
    
    mfl_move_current_to_digit(&c);
    *rows = mfl_read_nexus_type_int(&c);
    dbg_printf("\nRows: %i; ", *rows);
    
    mfl_move_current_to_digit(&c);
    *cols = mfl_read_nexus_type_int(&c);
    dbg_printf("\nColumns: %i.\n\n", *cols);
}


char* tui_get_simple_table_matrix(const char* input_table)
{
    int i = 0;
    char *ch = (char*)input_table;
    char *matrix = NULL;
    
    mfl_move_past_symbol(&ch, ';');
    ++ch;
    
    do {
        ++ch;
        if (!isspace(*ch)) {
            ++i;
        }
    } while (*ch);
    
    // Reset ch
    ch = (char*)input_table;
    mfl_move_past_symbol(&ch, ';');
    ++ch;
    
    matrix = (char*)malloc((i+1)*sizeof(char));
    
    i = 0;
    do {
        if (!isspace(*ch)) {
            matrix[i] = *ch;
            ++i;
        }
        ++ch;
    } while (*ch);
    
    matrix[i] = *ch;
    
    return matrix;
    
}

void tui_simple_table_parser(const char* input_table, mfl_handle_s* test_handle)
{
    /*  Reads a simple table. Table still requires some formatting, consistent with
     *  but without all the extra Nexus overhead. The basic format will be:
     *
     *  r c;
     *  eeeeee...ee(e_1_c)
     *  eeeeee...ee(e_2_c)
     *  .
     *  .
     *  .
     *  eeeeee...ee(e_r_c); [<- terminal semicolon]
     *
     *  Square brackets will be ignored.
     *
     */
    
    int num_rows = 0;
    int num_cols = 0;
    
    if (tui_check_simple_table_formatted(input_table)) {
        dbg_eprintf("unable to verify formatting");
        return;
    }
    
    tui_get_simple_table_dimensions(input_table, &num_rows, &num_cols);
    
    if (tui_check_simple_table_dimensions(input_table, num_rows, num_cols)) {
        dbg_eprintf("unable to verify dimensions");
        return;
    }
    
    test_handle->n_chars = num_cols;
    test_handle->n_taxa = num_rows;
    
    // Append the table into the handle
    test_handle->input_data = tui_get_simple_table_matrix(input_table);
    
    dbg_printf("\nChecking the stored matrix:\n");
    dbg_printf("%s\n\n", test_handle->input_data);
}


int tui_test_matrix_processing(mfl_handle_s *mfl_handle)
{
    char wagners[] = "1-3 6";
    char costmat[] = "8-10";
    char symbols[] = "012734";
    mfl_handle->ctypes_cmd[MFL_OPT_WAGNER] = wagners;
    mfl_handle->ctypes_cmd[MFL_OPT_COST_MATRIX] = costmat;
    mfl_handle->n_ctypes = 2;
    mfl_handle->format_symbols = symbols;
    mfl_handle->n_symbols = 6;
    
    mfl_matrix_t* internmatrix = NULL;
    
    internmatrix = mfl_create_internal_data_matrix((const mfl_handle_s*)mfl_handle);
    
    tui_print_out_converted_matrix(internmatrix, mfl_handle->n_taxa, mfl_handle->n_chars, true);
    
    mfl_destroy_mfl_matrix(internmatrix, mfl_handle->n_taxa, mfl_handle->n_chars);
    dbg_printf("\nEnd test for matrix conversion\n\n");
    
    
    // Print the contents of the data matrix.
    
    
}