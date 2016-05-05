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


void tui_print_node_bipartition(mfl_node_t* n)
{
    int i = 0;
    int j = 0;
    mfl_bitfield_t bitfield = 1;
    
    for (i = 0; i < n->nodet_bipart->bts_nfields; ++i) {
        for (j = 0; j < MFL_BTS_IN_BITSET; ++j) {
            if (n->nodet_bipart->bts_bitfields[i] & (bitfield << j)) {
                dbg_printf("*");
            }
            else {
                dbg_printf(".");
            }
        }
    }
    
    dbg_printf("\n");
}

void tui_partition_print_traversal(mfl_node_t* n)
{
    mfl_node_t* p = NULL;
    
    if (n->nodet_tip) {
        return;
    }
    
    p = n->nodet_next;
    do {
        tui_partition_print_traversal(p->nodet_edge);
        p = p->nodet_next;
    } while (p != n);
    
    tui_print_node_bipartition(n);
}

void tui_print_charstate_bits(const mfl_charstate cell, const int max_states)
{
    int i = 1;
    
    i = i << max_states;
    
    do {
        if (i & cell) {
            dbg_printf("1");
        }
        else {
            dbg_printf("0");
        }
        i = i >> 1;
    } while (i);
}

void tui_print_out_converted_matrix(mfl_matrix_t *matrix, int num_taxa, int num_chars, bool printbits)
{
    int i = 0;
    int j = 0;
    
    dbg_printf("Printing this obfuscated matrix:\n");
    
    dbg_printf("Char numbers:\n");
    dbg_printf("\t");
    for (i = 0; i < num_chars; ++i) {
        if (matrix->mat_matrix[i]) {
            dbg_printf("%i ", matrix->mat_matrix[i]->cv_num_gaps);
        } else {
            dbg_printf("# ");
        }
    }
    dbg_printf("\n");
    dbg_printf("\t");
    for (i = 0; i < num_chars; ++i) {
        if (matrix->mat_matrix[i]) {
            dbg_printf("%i ", matrix->mat_matrix[i]->cv_parsim_method);
        } else {
            dbg_printf("# ");
        }
    }
    
    dbg_printf("\nMatrix:\n");
    for (i = 0; i < num_taxa; ++i) {
        dbg_printf("\t");
        for (j = 0; j < num_chars; ++j) {
            dbg_printf("%s", matrix->mat_matrix[j]->cv_character_cells[i]);
            if (printbits) {
                dbg_printf("      ");
            }
        }
        if (printbits) {
            dbg_printf("\n\t");
            for (j = 0; j < num_chars; ++j) {
                tui_print_charstate_bits(matrix->mat_matrix[j]->cv_chardata[i], 5);
                dbg_printf(" ");
            }
            dbg_printf("\n\n");
        }
    }
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
    
    inputfilestr = (char*)malloc((filesize+1) * sizeof(char));
    
    
    do {
        inputfilestr[i] = fgetc(input);
        ++i;
    } while (i < (filesize-1));
    inputfilestr[filesize-1] = '\0';
    //inputfilestr[filesize] = '\0';
    
    //dbg_printf("The input file: \n");
    //dbg_printf("%s\n\n", inputfilestr);
    
    return inputfilestr;
}


void tui_destroy_file_string(char* oldinput)
{
    free(oldinput);
}


tui_testcmd_t tui_read_cmdline_argv(const char* arg2)
{
    //tui_testcmd_t ret = TUI_CMD_MAX;
    
    if (strcmp(arg2, "matrix")) {
        return TUI_CMD_MATRIX_INPUT_BASIC;
    }
    else if (strcmp(arg2, "newick")) {
        return TUI_CMD_NEWICK_INPUT_BASIC;
    }
    else if (strcmp(arg2, "command")) {
        return TUI_CMD_COMMAND_BASIC;
    }
    else {
        return TUI_CMD_MAX;
    }
}


mfl_handle_s* tui_parse_test_file(const char* arg1, const char* arg2)
{
    
    FILE* inputfile;
    
    char* filstr = NULL;
    
    if (!(inputfile = fopen(arg1, "r"))) {
        dbg_eprintf("file does not exist.\n");
        exit(-1);
    }
    else {
        dbg_printf("\nProcessing file %s . . .\n", arg1);
    }
    
    
    filstr = tui_readfile_to_str(inputfile);
    
    mfl_handle_s* testhandle;
    testhandle = mfl_t2s(mfl_create_handle());
    
    // Do stuff to the file here.
    /*  Options for arg2
     *
     *  matrix
     *  newick
     *  commands
     *
     */
    
    //Processing as a matrix
    tui_simple_table_parser((const char*)filstr, testhandle);
    
    dbg_printf("Test matrix processing:\n\n");
    tui_test_matrix_processing(testhandle);
    
    
    tui_destroy_file_string(filstr);
    mfl_destroy_handle(testhandle);
    fclose(inputfile);
    
    return testhandle;
    
}


/*
 * I think we can incorporate the NCL for I/O, I just need to figure out 
 * how to do that. Will work on that below.
 */
void tui_parse_test_infile(char *infile)
{
    int i = 0;
    unsigned int nchar = 0;
    
    MultiFormatReader infile_reader(-1);
    //infile_reader.ReadFilepath(*infile, MultiFormatReader::NEXML_FORMAT);
    
    NxsTaxaBlock* taxa = infile_reader.GetTaxaBlock(0);
    
    NxsCharactersBlock* nex_characters = infile_reader.GetCharactersBlock(taxa, 0);
 
    nchar = nex_characters->GetNChar();
    
}