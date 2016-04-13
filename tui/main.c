#include "ncl/ncl.h"
#include "morphy.h"
#include "tuimfy.h"
void mfl_test_newick_stuff();
void tui_test_character_stuff();


/* 
 *
 *   MORPHY UNIT TESTS
 *
 */
 
// == mfl_characters

 
// == mfl_evaluate

 
// == mfl_newick

 
// == mfl_starttree

 
// == mfl_tree
//---------------------------

 
// == mfl_brswap
//-----------------------------

// == funtion interactions

int tui_build_destroy_tree_from_binary_newick(mfl_handle_t *mfl_handle, char *input_newick)
{
    int num_taxa = 0;
    
    if (!mfl_handle && !input_newick) {
        
        input_newick = "((( (1,4) ( (5,) (2,6), 3))));";
        
    }
}
 
/*
 *
 *   END MORPHY UNIT TESTS
 *
 */


int tui_getting_numstates_test()
{
    int num_taxa = 5;
    int num_chars = 16;
    int num_states = 0;
    
    /*Too many states test. */
    
    //                                        1111111
    //                               1234567890123456
    char too_many_states_matrix[] = "1234567890-AABCD"
                                    "EFGHIJKLMNOPQRST"
                                    "UVWXYZabcdefghij"
                                    "klmnopqrstuvwxyz"
                                    "+0001203--1234@0;";
    
    
    num_states = mfl_get_numstates_from_matrix(too_many_states_matrix);
    
    if (num_states != 0) {
        dbg_printf("ERROR in call to mfl_get_numstates_from_matrix() by %s, line: %i", __FILE__, __LINE__);
        dbg_printf("\t Failed to recognize state proposal overload\n");
        return 1;
    }
    else {
        dbg_printf("== PASSED == %s, line: %i too many states check\n", __FILE__, __LINE__);
    }
    
    char just_enough_states_matrix[] =  "1234567890-AABCD"
                                        "EFGHIJKLMNOPQRST"
                                        "UVWXYZabcdefghij"
                                        "klmnopqrstuvwxyz"
                                        "0001203--1234@0A;";
    
    
    num_states = mfl_get_numstates_from_matrix(just_enough_states_matrix);
    
    if (num_states == 0) {
        dbg_printf("ERROR in call to mfl_get_numstates_from_matrix() by %s, line: %i", __FILE__, __LINE__);
        return 1;
    }
    else if (num_states != MORPHY_MAX_STATE_NUMBER){
        dbg_printf("ERROR in call to mfl_get_numstates_from_matrix():by %s, line: %i", __FILE__, __LINE__);
        dbg_printf("\t Input states does not match max states or function failed to count correct number of states\n");
    }
    else {
        dbg_printf("== PASSED == %s, line: %i just enough states check\n", __FILE__, __LINE__);
    }

    
}

int main (int argc, char *argv[])
{
    dbg_printf("\n\t****************************************\n\n");
    dbg_printf("\t  Welcome to the Morphy Test Interface\n\n");
    dbg_printf("\t****************************************\n\n\n");
    
    
    /* Begin the new bit */
    
    if (argc == 3) {
        dbg_printf("Processing file %s as %s . . .\n\n", argv[1], argv[2]);
        tui_parse_test_file(argv[1], argv[2]);
    }
    else if (argc > 3) {
        dbg_eprintf("too many arguments.\n");
        exit(1);
    }
    else if (argc == 1) {
        dbg_printf("argv[0] == %s...\n\m", argv[0]);
        //exit(0);
    }
    
    
    /* Everything below here should be written into or replaced by a real test. */
    dbg_printf("Testing include set values:\n");
    tui_test_character_stuff();
    tui_getting_numstates_test();
    dbg_printf("\n\n");
    
    dbg_printf("Generate a new tree\n\n");
    
    mfl_test_newick_stuff();
    
    mfl_tree_t *newtree;
    int num_taxa = 5;
    int num_nodes = 0;
    
    newtree = mfl_alloctree_with_nodes(num_taxa);
    
    /* Setting up a binary fork to start testing some of the new functions. Currently hard-coded*/
    mfl_create_binary_fork(newtree->treet_treenodes[num_taxa], newtree->treet_treenodes[0], newtree->treet_treenodes[1], newtree->treet_treenodes);
    
    /* Adding a hard-coded branch to the tree */
    mfl_make_new_n_ary_ring_node(newtree->treet_treenodes[num_taxa+1], 2, newtree->treet_treenodes);
    mfl_join_node_edges(newtree->treet_treenodes[2], newtree->treet_treenodes[num_taxa + 1]->nodet_next);
    mfl_insert_branch(newtree->treet_treenodes[num_taxa + 1], newtree->treet_treenodes[num_taxa + 1]->nodet_next->nodet_next, newtree->treet_treenodes[1]);
    
    dbg_printf("Free the tree\n\n");
    
    mfl_free_tree(newtree);
    
    dbg_printf("Exit the program\n\n");
    
	return 0;
}
