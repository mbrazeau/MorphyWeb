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
    bool rooting_specifier = false;
    
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

    if (!terminal_semicolon) {
        dbg_printf("Error in mfl_check_valid_newick(): no terminal semicolon\n");
        return is_valid = 0;
    }
    
    if (!rooting_specifier) {
        dbg_printf("Error in mfl_check_valid_newick(): no rooting specifier\n");
        return is_valid = 0;
    }
    
    if (n_open_brackets != n_closed_brackets) {
        dbg_printf("Error in mfl_check_valid_newick(): bracket number mismatch: ");
        if (n_open_brackets > n_closed_brackets) {
            dbg_printf("open brackets outnumber closed brackets\n");
        }
        else {
            dbg_printf("closed brackets outnumber opening brackets\n");
        }
        return is_valid = 0;
    }
    
    if (n_commas != n_closed_brackets) {
        dbg_printf("Error in mfl_check_valid_newick(): node-comma mismatch\n");
        return is_valid = 0;
    }
    
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
        dbg_printf("newick string is rooted\n");
    }
#ifdef MFY_DEBUG
    else if (*newick_string != 'R' && *newick_string != 'U') {
        is_rooted = false;
        dbg_printf("Error in mfl_newick_string_is_rooted(): Newick data has no rooting specification\n");
    }
#endif
    return is_rooted;
}


int mfl_read_newick_int(char **newick_position)
{
    int tip_number = 0;
    
    // Read char until ',' or ')' character, converting an int.
    do {
        tip_number = 10 * tip_number;
        tip_number = tip_number + (**(newick_position) - '0'); // This should convert the value at the current newick position to an int.
        ++(*newick_position);
    } while (**newick_position != ',' && **newick_position != ')');
    
    dbg_printf("Tip to retrieve: %i\n", tip_number);
    
    return tip_number;
}


mfl_node_t * mfl_traverse_newick_recursively(char **newick_position, mfl_nodearray_t nodearray)
{
    
    mfl_node_t *new_parent = NULL;  // A new parent node that will be returned from this function
    mfl_node_t *new_child = NULL;   // Child nodes returned from calls on this function
    mfl_node_t *node_ptr = NULL;    //
    int tip_number = 0;
    
    if (**newick_position == ';') {
        dbg_printf("ERROR in mfl_traverse_newick_recursively(): no defined operation for terminal semicolon\n");
        return NULL;
    }
    
    do {
        if (**newick_position == '(') {
            new_parent = mfl_get_next_available_node(NULL, NULL); // FIX THESE NULLS BEFORE CONTINUING
            new_child = mfl_traverse_newick_recursively(newick_position, nodearray);
            
            // Put new child into ring
            node_ptr = mfl_alloc_node();
            mfl_join_node_edges(node_ptr, new_child);
            mfl_insert_node_in_ring(new_parent, node_ptr);  // MAKE SURE THIS FUNCTION IS COMPATIBLE WITH THIS USAGE
            node_ptr = NULL;
            
            return  new_parent;
        }
        else if (**newick_position >= '1' || **newick_position <= '9') {
            tip_number = mfl_read_newick_int(newick_position);
            return nodearray[tip_number-1];
        }
        
        //++newick_position;
        
    } while (**newick_position != ')');
    
    
    
    return new_parent;
}

char *mfl_find_next_opening_bracket_in_newick(char *newick_tree)
{
    if (*newick_tree == '(') {
        return newick_tree;
    }
    
    do {
        ++newick_tree;
    } while (*newick_tree != '(');
    
    return newick_tree;
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
    if (*newick_tree != '(') {
        newick_tree = mfl_find_next_opening_bracket_in_newick(newick_tree);
    }
    
    // Begin recursion on the string, processing the nodearray as it proceeds
    
    
    
    return tree_from_newick;
}

void mfl_test_newick_stuff()
{
    char temp_example_newick_for_writing1[] = "temp_examp1=[&R] (2,((3,4),(5,1)));";
    char temp_example_newick_for_writing2[] = "temp_examp2=[&R] (2,(6,((3,4),(5,1))));";
    char temp_example_newick_for_writing3[] = "temp_examp3=[&R] (2,(6,((3,4),5),1));";
    char temp_example_newick_for_writing4[] = "temp_examp3=[&R] (2,((3,4) , (5,6) , 7));";
    
    char *sample_newick = NULL;
    
    mfl_newick_string_is_rooted(temp_example_newick_for_writing1);
    mfl_newick_string_is_rooted(temp_example_newick_for_writing2);
    mfl_newick_string_is_rooted(temp_example_newick_for_writing3);
    
    sample_newick = temp_example_newick_for_writing1;
    
    mfl_traverse_newick_recursively(&sample_newick, NULL);
    
    mfl_tree_t *anewtree = mfl_alloctree_with_nodes(5);
    
    dbg_printf("\n\nTesting function to find start of newick string:\n");
    
    sample_newick = temp_example_newick_for_writing1;
    sample_newick = mfl_find_next_opening_bracket_in_newick(sample_newick);

    ++sample_newick;
    dbg_printf("\nsample_newick before mfl_read_newick_int(): %c\n\n", *sample_newick);
    mfl_read_newick_int(&sample_newick);
    dbg_printf("sample_newick after mfl_read_newick_int(): %c\n\n", *sample_newick);
    
    mfl_traverse_newick_recursively(&sample_newick, NULL);
    dbg_printf("\n\n");
}
