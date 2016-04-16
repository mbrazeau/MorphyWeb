
#include "morphy.h"
#include "tuimfy.h"



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

void tui_print_out_converted_matrix(mfl_matrix_t *matrix, int num_taxa, int num_chars)
{
    int i = 0;
    int j = 0;
    
    dbg_printf("Printing this obfuscated matrix:\n");
    
    for (i = 0; i < num_taxa; ++i) {
        for (j = 0; j < num_chars; ++j) {
            dbg_printf("%s ", matrix->mat_matrix[j]->cv_character_cells[i]);
        }
        dbg_printf("\n\n");
    }
}

void tui_test_character_stuff()
{
    int i = 0;
    int num_taxa = 0;
    int num_chars = 0;
    char subcmd1[] = "ExSet * Exclude= 1-5 8 17;";
    char subcmd2[] = "Exclude = 1-5, 8 17;";
    char subcmd3[] = "18-51 100";
    char subcmd4[] = "10-15 2-6";
    
    //                        1111111
    //               1234567890123456
    char *matrix =  "120{12}000-000??001"
    "3001110001101120"
    "00(12?3)0001110010030"
    "{123}0(123)0001110010030;";
    
    num_chars = 16;
    num_taxa = 4;
    int num_states = 0;
    
    dbg_printf("Doing matrixy stuff...\n\n");
    dbg_printf("This is the matrix to convert (without linebreaks):\n");
    dbg_printf("%s\n\n", matrix);
    
    num_states = mfl_get_numstates_from_matrix(matrix);
    
    mfl_check_nexus_matrix_dimensions(matrix, num_taxa, num_chars);
    
    mfl_matrix_t *testmatrix = NULL;
    
    testmatrix = mfl_create_mfl_matrix(num_taxa, num_chars);
    
    mfl_setup_new_empty_matrix(testmatrix, num_states, num_taxa, num_chars);
    mfl_populate_chartype_character_vector(testmatrix, matrix, num_chars, num_taxa);
    
    tui_print_out_converted_matrix(testmatrix, num_taxa, num_chars);
    
    /* Need to calculate num_states!*/
    
    mfl_destroy_mfl_matrix(testmatrix, num_states, num_taxa, num_chars);
    
    
    /*char *subcmd = NULL;
     subcmd = subcmd4;
     dbg_printf("Processing the following subcommand: %s\n\n", subcmd);
     
     bool* includelist = mfl_read_nexus_exset_subcmd(subcmd, num_chars);
     
     dbg_printf("\nIncluding characters:\n");
     for (i = 0; i < num_chars; ++i) {
     if (includelist[i]) {
     dbg_printf("%i ", i + 1);
     }
     }
     dbg_printf("\n");
     
     dbg_printf("\nExcluding characters:\n");
     for (i = 0; i < num_chars; ++i) {
     if (!includelist[i]) {
     dbg_printf("%i ", i + 1);
     }
     }
     dbg_printf("\n");
     
     mfl_free_inclusion_list(includelist);*/
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
        dbg_printf("argv[0] == %s...\n\n", argv[0]);
        //exit(0);
    }
    
    
    /* Everything below here should be written into or replaced by a real test. */
    dbg_printf("Testing include set values:\n");
    tui_test_character_stuff();
    tui_getting_numstates_test();
    dbg_printf("\n\n");
    
    dbg_printf("Test Newick stuff\n\n");
    mfl_test_newick_stuff();
    
    dbg_printf("Test treecheck stuff\n\n");
    tui_test_checktree_();
    
    dbg_printf("\n\nGoodbye!\n\n");
    
	return 0;
}
