//
//  mfl_newick.c
//  Morphy
//
//  Functions for reading and writing trees in Newick format
//

#include "morphy.h"

int mfl_is_valid_newick(char *newick_input)
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

int mfl_count_internal_nodes_in_newick(char *newick_string)
{
    
    int num_open = 0;
    int num_closed = 0;
    int num_comma = 0;
    
    if (*newick_string != '(') {
        dbg_printf("ERROR in function calling mfl_count_internal_nodes_in_newick(): pointer is not at start of valid Newick string\n");
        return 0;
    }
    
    while (*newick_string != ';') {
        if (*newick_string == '(') {
            ++num_open;
        }
        if (*newick_string == ')') {
            ++num_closed;
        }
        if (*newick_string == ',') {
            ++num_comma;
        }
        ++newick_string;
    }
    
    return num_comma;
    
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


int mfl_seek_largest_tip_number_newick(char *newick_string)
{
    int retrieved_tip = 0;
    int largest_tip_number = 0;
    
    dbg_printf("The string to analyse: %s\n", newick_string);
    if (*newick_string != '(') {
        newick_string = mfl_find_next_opening_bracket_in_newick(newick_string);
    }
    
    while (*newick_string) {
        if (isdigit(*newick_string)) {
            retrieved_tip = mfl_read_newick_int(&newick_string);
            if (retrieved_tip > largest_tip_number) {
                largest_tip_number = retrieved_tip;
            }
        }
        ++newick_string;
    }
    
    return largest_tip_number;
}


mfl_node_t * mfl_traverse_newick_recursively(char **newick_position, mfl_nodearray_t nodearray, int num_taxa, int num_nodes)
{
    
    mfl_node_t *new_parent = NULL;  // A new parent node that will be returned from this function
    mfl_node_t *new_child = NULL;   // Child nodes returned from calls on this function
    mfl_node_t *node_ptr = NULL;    // For new node bases
    int tip_number = 0;
    
    if (**newick_position == ';') {
        dbg_printf("ERROR in mfl_traverse_newick_recursively(): no defined operation for terminal semicolon\n");
        return NULL;
    }
    if (**newick_position != '(') {
        dbg_printf("ERROR in mfl_traverse_newick_recursively(): called on invalid symbol\n");
        return NULL;
    }
    
    ++(*newick_position);
    
    new_parent = mfl_get_next_available_node(&nodearray[num_taxa], num_nodes);
    
    do {
        if (**newick_position == '(') {
            node_ptr = mfl_alloc_node();
            mfl_insert_node_in_ring(new_parent, node_ptr);
            new_child = mfl_traverse_newick_recursively(newick_position, nodearray, num_taxa, num_nodes);
            mfl_join_node_edges(node_ptr, new_child);
            
        }
        if (isdigit(**newick_position)) {
            node_ptr = mfl_alloc_node();
            tip_number = mfl_read_newick_int(newick_position);
            mfl_insert_node_in_ring(new_parent, node_ptr);
            new_child = nodearray[tip_number-1];
            mfl_join_node_edges(node_ptr, new_child);
        }
        if (**newick_position == ',') {
            ++(*newick_position);
        }
        if (isspace(**newick_position)) {
            ++(*newick_position);
        }
        
    } while (**newick_position != ')');
    
    ++(*newick_position);
    
    return new_parent;
}


mfl_tree_t *mfl_convert_newick_to_mfl_tree_t(char *newick_tree, int num_taxa)
{
    int num_taxa_local = 0; // If the number of taxa is not supplied by the user, it is easy to calculate it from the Newick string
    int num_nodes = 0;
    mfl_tree_t *tree_from_newick = NULL;
    char **newick_position = NULL;
    
    if (*newick_tree != '(') {
        newick_tree = mfl_find_next_opening_bracket_in_newick(newick_tree);
    }
    
    if (!num_taxa) {
        /* Ideally, Newick trees will never be supplied without an NTAX value. 
         * However, if that fails, the safest operation is to find the largest
         * tip number in the string and allocate at least that many nodes. Then 
         * the tip array will always be large enough to accommodate at least
         * many terminals. */
        num_taxa_local = mfl_seek_largest_tip_number_newick(newick_tree);
        num_taxa_local = num_taxa_local + 1;
    }
    else {
        num_taxa_local = num_taxa;
    }
    
    num_nodes = mfl_calculate_number_of_nodes_to_allocate(num_taxa_local);
    
    tree_from_newick = mfl_alloctree_with_nodes(num_taxa_local);
    
    
    // Begin recursion on the string, processing the nodearray as it proceeds
    newick_position = &newick_tree;
    tree_from_newick->treet_root = mfl_traverse_newick_recursively(newick_position, tree_from_newick->treet_treenodes, num_taxa_local, num_nodes);
    
    return tree_from_newick;
}

void mfl_test_newick_stuff()
{
    /* This function will be eliminated from the library. */
    
    char temp_example_newick_for_writing1[] = "temp_examp1=[&R] (2,((3,4),(5,1)));";
    char temp_example_newick_for_writing2[] = "temp_examp2=[&R] (2,(6,((3,4),(5,1))));";
    char temp_example_newick_for_writing3[] = "temp_examp3=[&R] (2,(6,((3,4),5),1));";
    char temp_example_newick_for_writing4[] = "temp_examp3=[&R] (2,((3,4),(5,20),1));"; // Polytomy and multi-digit tip number not in sequence
    
    char *sample_newick = NULL;
    
    int num_taxa = 0;
    int num_nodes = 0;
    
    mfl_tree_t *tree_from_newick = NULL;
    
    int largest = 0;
    
    largest = mfl_seek_largest_tip_number_newick(temp_example_newick_for_writing1);
    dbg_printf("Largest in example 1: %i\n", largest);
    largest = mfl_seek_largest_tip_number_newick(temp_example_newick_for_writing4);
    dbg_printf("Largest in example 4: %i\n", largest);
    
    
    sample_newick = mfl_find_next_opening_bracket_in_newick(temp_example_newick_for_writing1);
    dbg_printf("The string after finding the opening bracket:\n");
    dbg_printf("%s\n", sample_newick);
    
    tree_from_newick = mfl_convert_newick_to_mfl_tree_t(temp_example_newick_for_writing4, 0);
    
    mfl_free_tree(tree_from_newick, 5, (2*5-1));
    
}
