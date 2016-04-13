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
        ++c;
    } while (*c != ';');
    
    if (expected == actual) {
        return 0;
    }
    else {
        dbg_eprintf("input matrix has unexpected dimensions");
        exit(actual);
    }
}

void tui_get_simple_table_dimensions(const char*table, int* rows, int* cols)
{
    tui_check_simple_table_formatted(table);
    
    
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
     *  eeeeee...ee(e1_c)
     *  eeeeee...ee(e2_c)
     *  .
     *  .
     *  .
     *  eeeeee...ee(er_c); [<- terminal semicolon]
     *
     *  Square brackets will be ignored.
     *
     */
    
    char* inputtable =  "4 10;"
                        "0012300001"
                        "0012300001"
                        "0012300001"
                        "0012300001;";
    
    tui_check_simple_table_formatted(inputtable);
    tui_check_simple_table_dimensions(inputtable, 21, 10);
    
    
}