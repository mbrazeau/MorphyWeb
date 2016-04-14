/*
 *  mfl_newick.c
 *
 *  THE MORPHY FUNCTION LIBRARY
 *  A library for phylogenetic analysis with emphasis on parsimony and
 *  morphology (but someday other methods)
 *
 *  Copyright (C) 2016  by Martin D. Brazeau, Thomas Guillerme,
 *  and Chris Desjardins
 *
 *  Created by Martin Brazeau and Chris Desjardins on 1/26/12.
 *  Some data structs, routines and ideas derived from:
 *      - The PHYLIP package by Joe Felsenstein
 *          <http://evolution.genetics.washington.edu/phylip.html>
 *      - MrBayes by John Huelsenbeck and Fredrik Ronquist
 *          <http://mrbayes.sourceforge.net/>
 *
 *  Any bugs, errors, inefficiences and general amateurish handling are our own
 *  and most likely the responsibility of MDB. We make no guarantees.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


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
        if (*(newick_input + i) == '&') {
            rooting_specifier = true;
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

bool mfl_newick_tree_is_rooted(char *newick_string)
{
    bool is_rooted = false;
    
    while (*newick_string != '&') {
        ++newick_string;
    }
    
    ++newick_string;
    
    if (*newick_string != 'R' && *newick_string != 'U' && *newick_string != 'r' && *newick_string != 'u') {
        is_rooted = false;
        dbg_printf("Error in mfl_newick_string_is_rooted(): Newick data has no rooting specification\n");
    }
    else if (*newick_string == 'R' || *newick_string == 'r') {
        is_rooted = true;
        dbg_printf("Newick string is rooted\n");
    }
    else {
        dbg_printf("Newick string is unrooted\n");
    }
    
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


mfl_node_t * mfl_traverse_newick_recursively(char **newick_position, mfl_nodearray_t nodearray, int num_taxa)
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
    
    new_parent = mfl_get_next_available_node(nodearray);
    new_parent->nodet_next = new_parent;
    
    if (new_parent->nodet_index == 7) {
        dbg_printf("Just wait\n");
    }
    
    do {
        if (**newick_position == '(') {
            node_ptr = mfl_get_next_available_node(nodearray);
            mfl_insert_node_in_ring(new_parent, node_ptr);
            new_child = mfl_traverse_newick_recursively(newick_position, nodearray, num_taxa);
            mfl_join_node_edges(node_ptr, new_child);
            
        }
        if (isdigit(**newick_position)) {
            node_ptr = mfl_get_next_available_node(nodearray);
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


void mfl_test_newick_setup(mfl_node_t *node)
{
    /* RENAME AND MOVE TO TUI ALONG WITH OTHER TESTS */
    mfl_node_t *p = NULL;
    int i = 0;
    
    if (node->nodet_tip) {
        dbg_printf("tip: %i\n", node->nodet_tip);
        return;
    }
    
    p = node->nodet_next;
    
    do {
        mfl_test_newick_setup(p->nodet_edge);
        p = p->nodet_next;
    } while (p != node);
    
    dbg_printf("Node %i: ", node->nodet_index);
    
    p = node->nodet_next;
    
    do {
        ++i;
        dbg_printf("child %i: %i ", i, p->nodet_edge->nodet_index);
        p = p->nodet_next;
    } while (p != node);
    dbg_printf("\n");
    
    return;
}

mfl_tree_t *mfl_convert_newick_to_mfl_tree_t(char *newick_tree, int num_taxa)
{
    int num_taxa_local = 0; // If the number of taxa is not supplied by the user, it is easy to calculate it from the Newick string.
    mfl_tree_t* tree_from_newick = NULL;
    char *newicktr_copy = NULL;
    char **newick_position = NULL; // A pointer to the Newick string that can be incremented during the recursion on the string.
    
    /* Do safety checks on the Newick string. */
    if (!mfl_is_valid_newick(newick_tree)) {
        dbg_printf("ERROR in function calling mfl_convert_newick_to_mfl_tree_t(): string passed to function is not in valid Newick format\n");
        return NULL;
    }
    
    newicktr_copy = newick_tree;
    
    /* Ensure that the recursion begins on the opening bracket */
    if (*newicktr_copy != '(') {
        newicktr_copy = mfl_find_next_opening_bracket_in_newick(newicktr_copy);
    }
    
    /* Ideally, Newick trees will never be supplied without an NTAX value.
     * However, if that fails, the safest operation I can come up with is
     * to find the largest tip number in the string and allocate at least
     * that many nodes. Then the tip array will always be large enough to
     * accommodate at leastmany terminals. MDB. */
    if (!num_taxa) {
        num_taxa_local = mfl_seek_largest_tip_number_newick(newicktr_copy);
    }
    else {
        num_taxa_local = num_taxa;
    }
    
    newick_position = &newicktr_copy;
    
    tree_from_newick = mfl_alloctree_with_nodes(num_taxa_local);
    
    tree_from_newick->treet_root = mfl_traverse_newick_recursively(newick_position, tree_from_newick->treet_treenodes, num_taxa_local);
    
    dbg_printf("The newick string processed: %s\n", newick_tree);
    
    mfl_test_newick_setup(tree_from_newick->treet_root);
    
    /* Process the rooting options*/
    if (!mfl_newick_tree_is_rooted(newick_tree)) {
        /* FINISH: Unroot the tree. */
    }
    
    return tree_from_newick;
}

// Finding the number of digits in an integer
int mfl_number_of_digits_in_integer(int n)
{
    //Variables
    int count=0;

    //Counting the digits
    while(n!=0)
    {
        n/=10;
        ++count;
    }
  return(count);
}

// Getting the power of an integer (equivalent to pow from math but does not output a double!)
int mfl_power(const int base, const unsigned int exp) {

    //Variables
    int i, result = 1;

    //Calculating the power as an integer
    for (i = 0; i < exp; i++)
        result *= base;
    return result;
}

//Getting the total number of digit in a sequence from 1 to n
int mfl_number_of_digits_in_sequence(const int n)
{
    
    //Variables
    int digits_in_n, digits_in_sequence;

    // Get the number of digits in n
    digits_in_n = mfl_number_of_digits_in_integer(n);

    // Get the total number of digits in the sequence from 1 to n conditional on n
    /* TG: For n being an integer and x being the number of digit in n, the
    * number of digits in the sequence from 1 to n is equal to:
    * "x*n - ( 10^(x-1) + 10^(x-2), ..., + 10^(x-x) - (x-1) )"
    * or, more conveniently:
    * "x*n - ( ((10^x) - 1)/(10 - 1) - x )"" */

    digits_in_sequence = digits_in_n * n - ( (mfl_power(10,digits_in_n) - 1)/(10 - 1) - digits_in_n); 

    return(digits_in_sequence);
}

//Getting the maximum number of characters in newick string
int mfl_number_of_characters_in_newick(const int num_taxa)
{
    //Variables
    int num_brackets = 2*(num_taxa - 1), num_commas = num_taxa - 1, newick_size;

    //Get the newick string size
    newick_size = mfl_number_of_digits_in_sequence(num_taxa) + num_brackets + num_commas + 6; // The final +6 includes the closing semi-colon (1), the root header (4) and a space between the root header and the start of the newick string (1)

    return newick_size;
}

//Getting number of taxa from a mfl_tree_t
int mfl_tree_t_max_number_of_taxa(const mfl_tree_t *input_tree)
{
    /* MDB says: This function doesn't really return the number 
     * num_taxa_active, but a total num_taxa at time of tree allocation.
     * The number num_taxa_active should be determined by a traversal
     * on the tree, since it should refer exclusively to the number of
     * taxa active in the assembled tree, rather than in the allocated tree
     * 
     */
    
    int num_taxa_active = 0, count = 0;
    
    while(input_tree->treet_treenodes[count] != '\0') {
        
        //Count the number of nodet_tip
        if(input_tree->treet_treenodes[count]->nodet_tip != 0) {
            ++num_taxa_active;
        }
        
        ++count;
    }
    
    return num_taxa_active;
}


//char *source = "This is an example.";
//int i;
//
//for (i = 0; i < strlen(source); i++){
//    
//    printf("%c", source[i]);
//    
//}


// Traversing the tree to spit out the newich characters
char* mfl_traverse_tree_to_print_newick_char_recursive(mfl_node_t *start, char *newick_tree_out, int &count)
{
    int i = 0;
    int tip_length;
    
    // Set the node to start
    mfl_node_t *node = start;
    
    if (node->nodet_tip != 0) {
        // convert the tip into a character
        tip_length = mfl_number_of_digits_in_integer(node->nodet_tip);
        char tip_num[tip_length];
        sprintf(tip_num, "%d", node->nodet_tip);
        // print a tip (depending on the number of integer in there)
        i = mfl_number_of_digits_in_integer(node->nodet_tip);
        for (i = 0; i < tip_length; ++i) {
            newick_tree_out[count] = tip_num[i];
            count++;
        }
        
        return newick_tree_out;
        
    } else {
        // open a clade
        newick_tree_out[count] = '(';
        count++;
    }
    
    // Move to the next node
    node = node->nodet_next;
    
    do {
        if (node != start->nodet_next) {
            newick_tree_out[count] = ',';
            count++;
        }
        // Go through the next nodes and print the next tip
        mfl_traverse_tree_to_print_newick_char_recursive(node->nodet_edge, newick_tree_out, count);
        
        // close a clade if the next edge is not a tip
        if(node->nodet_edge->nodet_tip == 0) {
            newick_tree_out[count] = ')';
            count++;
        }
        
        // Move to the next node
        node = node->nodet_next;
    
    // Stop once reached back the start tree
    } while (node != start);
    
    return newick_tree_out;
}


// Write newick string from tree_t structure
char* mfl_convert_mfl_tree_t_to_newick(mfl_tree_t *input_tree, int num_taxa_active)
{
    //Variables
    char *newick_tree_out = NULL;
    char *root_command;
    char *tree_end = ");\0";
    int newick_string_length = 0;
    int count = 0; // This is the counter for the recursive loop
    int i;
    
    //Infering max number of taxa
    if(num_taxa_active == 0) {
        //Infer number of taxa from mfl_tree_t
        num_taxa_active = mfl_tree_t_max_number_of_taxa(input_tree);
    }

    //Get the newick string length
    newick_string_length = mfl_number_of_characters_in_newick(num_taxa_active) + 2;  // the + 2 is for the terminal '\0'
    
    //Allocating memory to the newick
    newick_tree_out = (char*)malloc(newick_string_length * sizeof(char));
    
    if (!newick_tree_out) {
        dbg_printf("ERROR in mfl_convert_mfl_tree_t_to_newick(): unable to allocate memory for Newick string\n");
        return NULL;
    }
    else {
        memset(newick_tree_out, 0, newick_string_length * sizeof(char));
    }
    
    //Adding starting with the root command
    if (input_tree->treet_root->nodet_isroot != 0) {
        root_command = "[&R] ";
    } else {
        root_command = "[&U] ";
    }
    for (i = 0; i < 5; ++i) {
        newick_tree_out[i] = root_command[i];
        count++;
    }
    
    //Creating the newick string out
    newick_tree_out = mfl_traverse_tree_to_print_newick_char_recursive(input_tree->treet_root, newick_tree_out, count);
    
    //Closing the tree
    strcat(newick_tree_out, tree_end);
    
    return newick_tree_out;

}


void mfl_test_newick_stuff()
{
    /* This function will be eliminated from the library. */
    
    char temp_example_newick_for_writing1[] = "temp_examp1=[&R] (2,((3,4),(5,1)));";
    char temp_example_newick_for_writing2[] = "temp_examp2=[&R] (2,(6,((3,4),(5,1))));";
    char temp_example_newick_for_writing3[] = "temp_examp3=[&R] (2,(6,((3,4),5),1));";
    char temp_example_newick_for_writing4[] = "temp_examp4=[&U] (2,((3,4),(5,20),1));"; // Polytomy and multi-digit tip number not in sequence
    char temp_example_newick_for_writing5[] = "temp_examp5=[&R] (((((1,4),5),3),2),6);";
    char temp_example_newick_for_writing6[] = "temp_examp6=[&R] (((((1,4),5),3),2),6,(7,8));";
    char temp_example_newick_for_writing7[] = "temp_examp7=[&R] ((1000,856),(2,3),(56,4));";
    
    char *sample_newick = NULL;
    
    int num_taxa = 0;
    int num_nodes = 0;
    
    mfl_tree_t *tree_from_newick = NULL;
    
    int largest = 0;
    
    largest = mfl_seek_largest_tip_number_newick(temp_example_newick_for_writing1);
    dbg_printf("Largest in example 1: %i\n", largest);
    largest = mfl_seek_largest_tip_number_newick(temp_example_newick_for_writing4);
    dbg_printf("Largest in example 4: %i\n", largest);
    
    sample_newick = mfl_find_next_opening_bracket_in_newick(temp_example_newick_for_writing4);
    dbg_printf("The string after finding the opening bracket:\n");
    dbg_printf("%s\n", sample_newick);

    //tree_from_newick =  mfl_convert_newick_to_mfl_tree_t(temp_example_newick_for_writing6, 0);
    
    dbg_printf("Imported tree is: %s\n", temp_example_newick_for_writing7);
    tree_from_newick =  mfl_convert_newick_to_mfl_tree_t(temp_example_newick_for_writing7, 0);

    dbg_printf("\n\n\nTesting convert to newick\n");
    
    dbg_printf("%s\n", temp_example_newick_for_writing7);
    char *newick_string;
    newick_string = mfl_convert_mfl_tree_t_to_newick(tree_from_newick, 0);
    
    dbg_printf("The output newick string is:\n%s", newick_string);
    
    dbg_printf("\n\n\n");
    
    
    
    mfl_free_tree(tree_from_newick); // Some bug here
    
}