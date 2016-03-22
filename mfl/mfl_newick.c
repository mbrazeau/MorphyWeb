//
//  mfl_newick.c
//  Morphy
//
//  Functions for reading and writing trees in Newick format
//

#include "morphy.h"

int mfl_check_invalid_newick(char *newick_input)
{
    int i = 0;
    int is_valid = 1;
    
    int n_open_brackets = 0;
    int n_closed_brackets = 0;
    int n_commas = 0;
    bool terminal_semicolon = false;
    
    for (i = 0; *(newick_input + i) != '\0'; ++i) {
        if (*(newick_input + i) == '(') {
            ++n_open_brackets;
        }
        else if (*(newick_input + i) == ')') {
            ++n_closed_brackets;
        }
        else if (*(newick_input + i) == ',') {
            ++n_commas;
        }
        if (*(newick_input + i) == ';') {
            terminal_semicolon = true;
        }
    }

#ifdef MFY_DEBUG
    if (!terminal_semicolon) {
        dbg_printf("Error in mfl_check_valid_newick(): no terminal semicolon\n");
        //return is_valid = 0;
    }
    
    if (n_open_brackets != n_closed_brackets) {
        dbg_printf("Error in mfl_check_valid_newick(): bracket number mismatch\n");
        //return is_valid = 0;
    }
    
    if (n_commas != n_closed_brackets) {
        dbg_printf("Error in mfl_check_valid_newick(): node-comma mismatch\n");
        //return is_valid = 0;
    }
#endif
    
    return is_valid;
}


bool mfl_newick_string_is_rooted(char *newick_string)
{
    bool is_rooted = false;
    
    while (*newick_string != '&') {
        ++newick_string;
    }
    ++newick_string;
    if (*newick_string == 'R' || *newick_string == 'r') {
        is_rooted = true;
    }
#ifdef MFY_DEBUG
    else if (*newick_string != 'R' && *newick_string != 'u') {
        is_rooted = false;
        dbg_printf("Error in mfl_newick_string_is_rooted(): Newick data has no rooting specification\n");
    }
#endif
    return is_rooted;
}


char * mfl_traverse_newick_recursively(char * newick_position)
{
    
    if (*newick_position == ';') {
        return newick_position;
    }
    
    dbg_printf("%c", *newick_position);
    ++newick_position;
    mfl_traverse_newick_recursively(newick_position);
    
    return newick_position;
}

mfl_tree_t *mfl_convert_newick_to_mfl_tree_t(char *newick_tree, int num_taxa)
{
    int num_taxa_local = 0;
    mfl_tree_t *tree_from_newick = NULL;
    
    if (!num_taxa) {
        //num_taxa_local = mfl_count_taxa_in_newick_string();
    }
    else {
        num_taxa_local = num_taxa;
    }
    
    tree_from_newick = mfl_alloctree_with_nodes(num_taxa_local);
    
    // PRINCIPAL BUSINESS LOGIC;
    // Get pointer to start of the newick string [first opening bracket?]
    // Begin recursion on
    
    
    return tree_from_newick;
}

void mfl_test_newick_stuff()
{
    char temp_example_newick_for_writing1[] = "temp_examp=[&R] (2,((3,4),(5,1)));";
    char temp_example_newick_for_writing2[] = "temp_examp=[&R] (2,(6,((3,4),(5,1))));";
    char temp_example_newick_for_writing3[] = "temp_examp=[&R] (2,(6,((3,4),5),1));";
    
    mfl_traverse_newick_recursively(temp_example_newick_for_writing1);
}
