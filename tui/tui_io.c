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


mfl_handle_s* tui_parse_test_file(const char* arg1, const char* arg2)
{
    
    FILE* inputfile;
    
    char* filstr = NULL;
    
    if (!(inputfile = fopen(arg1, "r"))) {
        dbg_eprintf("file does not exist.\n");
        exit(-1);
    }
    else {
        dbg_printf("\nProcessing file % . . .\n", arg1);
    }
    
    
    filstr = tui_readfile_to_str(inputfile);
    
    mfl_handle_s* testhandle;
    testhandle = mfl_t2s(mfl_create_handle());
    
    // Do stuff to the file here.
    /*  Options for arg2
     *
     *  matrix
     *  newick
     */
    
    //Processing as a matrix
    tui_simple_table_parser((const char*)filstr, testhandle);
    
    tui_destroy_file_string(filstr);
    
    //mfl_destroy_handle(testhandle);
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