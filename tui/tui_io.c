//
//  tui_io.c
//  Morphy
//
//  Created by mbrazeau on 13/04/2016.
//
//

//#include <stdio.h>
#include "morphy.h"
#include "tuimfy.h"

void tui_ignore_nexus_comment(char **current)
{
    if (**current != '[') {
        dbg_eprintf("not on valid Nexus comment flag");
    }
    do {
        ++(*current);
    } while (**current != ']');
}

int tui_check_simple_table_formatted(const char* input_table)
{
    char *c = (char*)input_table;
    int semicolons = 0;
    int dimensions = 0;
    bool innum = false;
    
    do {
        if (*c == ';') {
            ++semicolons;
        }
        ++c;
    } while (*c);
 
    if (semicolons != 2) {
        dbg_eprintf("expected semicolon.");
        exit(semicolons);
    }
    
    // Move back to the start of the file
    c = (char*)input_table;
    
    do {
        if (isdigit(*c)) {
            innum = true;
        }
        if (innum) {
            if (isspace(*c)) {
                ++dimensions;
                innum = false;
            }
        }
        ++c;
    } while (*c != ';');
    
    if (innum) {
        if (dimensions || c != input_table) {
            ++dimensions;
        }
    }
    
    if (dimensions != 2) {
        dbg_eprintf("expected table dimensions.");
        exit(dimensions);
    }
    
    return 0;
}


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

char* tui_get_simple_table()
{
    
}

void tui_simple_table_parser(const char* input_table, mfl_handle_t test_handle)
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
    
    tui_check_simple_table_formatted(input_table);
    tui_get_simple_table_dimensions(input_table, &num_rows, &num_cols);
    tui_check_simple_table_dimensions(input_table, num_rows, num_cols);
    
}

char* tui_readfile_to_str(FILE *input)
{
    int i = 0, filesize = 0;
    char c = 0;
    char* inputfilestr = NULL;
    
    while (c != EOF) {
        ++filesize;
        c = fgetc(input);
    }

    dbg_printf("Input file is %i characters long.\n\n", filesize);
    rewind(input);
    
    inputfilestr = (char*)malloc(filesize * sizeof(char));
    
    
    do {
        inputfilestr[i] = fgetc(input);
        ++i;
    } while (i < (filesize-1));
    inputfilestr[filesize-1] = '\0';inputfilestr[filesize] = '\0';
    
    dbg_printf("The input file: \n");
    dbg_printf("%s\n\n", inputfilestr);
    
    return inputfilestr;
}

void tui_destroy_file_string(char* oldinput)
{
    free(oldinput);
}


int tui_parse_test_file(const char* arg1, const char* arg2)
{
    
    FILE* inputfile;
    
    char* filstr = NULL;
    
    if (!(inputfile = fopen(arg1, "r"))) {
        dbg_eprintf("file does not exist.\n");
        return(1);
    }
    else {
        dbg_printf("\nProcessing file % . . .\n", arg1);
    }
    
    
    filstr = tui_readfile_to_str(inputfile);
    
    mfl_handle_t testhandle;
    testhandle = mfl_create_handle();
    
    // Do stuff to the file here.
    /*  Options for arg2
     *
     *  matrix;
     *  newick
     */
    
    //Processing as a matrix
    tui_simple_table_parser((const char*)filstr, testhandle);
    
    tui_destroy_file_string(filstr);
    
    mfl_destroy_handle(testhandle);
    fclose(inputfile);
    
}

/*
 * I think we can incorporate the NCL for I/O, I just need to figure out 
 * how to do that. Will work on that below.
 */

void tui_Nexus_Newick_reader(char *filename)
{
    NxsTaxaBlock *taxa = new NxsTaxaBlock();
    NxsTreesBlock *trees = new NxsTreesBlock(taxa);
}

void tui_Nexus_matrix_reader(char *filename)
{
    
    NxsTaxaBlock *taxa = new NxsTaxaBlock();
    NxsAssumptionsBlock *assumptions = new NxsAssumptionsBlock(taxa);
    NxsCharactersBlock *characters = new NxsCharactersBlock(taxa, assumptions);
    
}
