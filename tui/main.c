
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
    "00(12?3)000?110010030"
    "{123}0(24)0001110010030;";
    
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
    mfl_populate_chartype_character_vectors(testmatrix, matrix, num_chars, num_taxa);
    
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


void tui_test_tree_printing()
{
    int num_taxa = 78;
    char *grid = mfl_drawtree_create_virtual_grid(num_taxa);
    
    char newick[] = "[&R] = tree PAUP_1 = [&U] (1,(2,(((((((((((((((((((((3,39),12),(11,(53,64))),30),(42,62)),48),(25,32)),74),21),((((((6,61),76),17),67),(8,45)),((((((((((14,22),38),(16,18)),((37,58),75)),(59,73)),15),26),68),(51,56)),36))),((((13,(40,((46,55),54))),49),((((((29,34),(33,63)),72),57),65),35)),23)),70),44),27),(31,43)),(((9,((19,41),(20,28))),24),47)),71),((4,10),69)),((50,78),52)),7),(((5,66),77),60))));";
    
    mfl_tree_t* printme = mfl_convert_newick_to_mfl_tree_t(newick, num_taxa);
    
    printme->treet_treenodes[0]->nodet_tipname = (char*)"Borosaurus jibberstoni";
    printme->treet_treenodes[1]->nodet_tipname = (char*)"Wayne";
    printme->treet_treenodes[2]->nodet_tipname = (char*)"Ur mom";
    printme->treet_treenodes[3]->nodet_tipname = (char*)"Taxon 1";
    printme->treet_treenodes[4]->nodet_tipname = (char*)"Taxon Two";
    printme->treet_treenodes[5]->nodet_tipname = (char*)"Foo";
    printme->treet_treenodes[6]->nodet_tipname = (char*)"Bar";
    printme->treet_treenodes[7]->nodet_tipname = (char*)"Fubarus";
    
    int firstrow = 0;
    mfl_drawtree_set_coords_traversal(printme->treet_root, &firstrow, grid, num_taxa);
    mfl_drawtree_draw_traversal(printme->treet_root, grid);
    
    grid = mfl_drawtree(printme);
    
    //dbg_printf("________\n");
    dbg_printf("          1.........2.........3.........4.........5.........6.........7.........8\n");
    dbg_printf("012345678901234567890123456789012345678901234567890123456789012345678901234567890\n");
    dbg_printf("%s", grid);
    
    return;
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
    
    dbg_printf("Testing the tree printing:\n");
    tui_test_tree_printing();
    dbg_printf("\nEnd tree print test\n");
    
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
