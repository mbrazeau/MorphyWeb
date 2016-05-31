
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
    
    tui_print_out_converted_matrix(testmatrix, num_taxa, num_chars, false);
    
    /* Need to calculate num_states!*/
    
    mfl_destroy_mfl_matrix(testmatrix, num_taxa, num_chars);
    
    
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
    char *grid = NULL;
    
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
    
    grid = mfl_drawtree(printme);
    
    //dbg_printf("________\n");
    dbg_printf("          1.........2.........3.........4.........5.........6.........7.........8\n");
    dbg_printf("012345678901234567890123456789012345678901234567890123456789012345678901234567890\n");
    dbg_printf("%s", grid);
    
    mfl_free_tree(printme);
    free(grid);
}

void tui_test_bipartition_setting()
{
    int num_taxa = 78;
    char nwktree[] = "tree1=[&R] (1,(2,(((((((((((((((((((((3,39),12),(11,(53,64))),30),(42,62)),48),(25,32)),74),21),((((((6,61),76),17),67),(8,45)),((((((((((14,22),38),(16,18)),((37,58),75)),(59,73)),15),26),68),(51,56)),36))),((((13,(40,((46,55),54))),49),((((((29,34),(33,63)),72),57),65),35)),23)),70),44),27),(31,43)),(((9,((19,41),(20,28))),24),47)),71),((4,10),69)),((50,78),52)),7),(((5,66),77),60))));";
    mfl_tree_t* t = mfl_convert_newick_to_mfl_tree_t(nwktree, num_taxa);
    
    mfl_set_bipartitions(t->treet_root);
    tui_partition_print_traversal(t->treet_root);
    mfl_free_tree(t);
}


void tui_nexus_reader(char* argv1)
{
    NxsReader *myreader = new NxsReader;
    myreader->ReadFilepath(argv1);
    
    NxsTaxaBlock *taxa = new NxsTaxaBlock();
    NxsAssumptionsBlock *pAssumptions = new NxsAssumptionsBlock(taxa);
    const NxsTransformationManager transmanager = pAssumptions->GetNxsTransformationManagerRef();
//    string stepmatrixname = transmanager.GetUserTypeNames();
//    const NxsIntStepMatrix inmatrix = transmanager.GetIntType(<#const std::string &name#>)
}


void tui_test_newick_stuff()
{
    /* This function will be eliminated from the library. */
    
    char temp_example_newick_for_writing1[] = "temp_examp1=[&U] (2,((3,4),(5,1)));";
    char temp_example_newick_for_writing2[] = "temp_examp2=[&R] (2,(6,((3,4),(5,1))));";
    char temp_example_newick_for_writing3[] = "temp_examp3=[&R] (2,(6,((3,4),5),1));";
    char temp_example_newick_for_writing4[] = "temp_examp4=[&U] (2,((3,4),(5,20),1));"; // Polytomy and multi-digit tip number not in sequence
    char temp_example_newick_for_writing5[] = "temp_examp5=[&R] (((((1,4),5),3),2),6);";
    char temp_example_newick_for_writing6[] = "temp_examp6=[&R] (((((1,4),5),3),2),6,(7,8));";
    char temp_example_newick_for_writing7[] = "temp_examp7=[&U] ((1000,856),(2,3),(56,4));";
    char temp_example_newick_for_writing8[] = "temp_examp8=[&R] (1,(2,(3,(4,5))));";
    
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
    
    dbg_printf("\n\n\n");
    tree_from_newick =  mfl_convert_newick_to_mfl_tree_t(temp_example_newick_for_writing2, 0);
    dbg_printf("Testing convert to newick\n");
    char *newick_string;
    newick_string = mfl_convert_mfl_tree_t_to_newick(tree_from_newick, true);
    dbg_printf("The output newick string is: %s\n", newick_string);
    dbg_printf("Unrooting tree:\n");
    mfl_unroot_tree(tree_from_newick);
    newick_string = mfl_convert_mfl_tree_t_to_newick(tree_from_newick, true);
    dbg_printf("The output newick string is (with polytomy): %s\n", newick_string);
    tree_from_newick =  mfl_convert_newick_to_mfl_tree_t(temp_example_newick_for_writing2, 0);
    mfl_unroot_tree(tree_from_newick);
    newick_string = mfl_convert_mfl_tree_t_to_newick(tree_from_newick, false);
    dbg_printf("The output newick string is (without polytomy): %s\n", newick_string);
    
    mfl_free_tree(tree_from_newick); // Some bug here
    
    
    dbg_printf("Putting the following trees in the buffer:\n%s\n%s\n%s\n", temp_example_newick_for_writing1, temp_example_newick_for_writing5, temp_example_newick_for_writing6);
    
    //Create a char** to trees
    char* list_of_newicks[] = {temp_example_newick_for_writing1, temp_example_newick_for_writing5, temp_example_newick_for_writing6};

    //Create a handle
    mfl_handle_s* testhandle;
    testhandle = mfl_t2s(mfl_create_handle());
    
    //Add the trees to the handle
    testhandle->n_input_newick_trees = 3;
    testhandle->input_newick_trees = list_of_newicks;
    
    //Create the buffer
    mfl_treebuffer_t *my_buffer = NULL;
    my_buffer = mfl_tree_from_newick_to_buffer(testhandle);
    
    //Print what's in the buffer
    dbg_printf("Trees in the buffer are:\n");
    newick_string = mfl_convert_mfl_tree_t_to_newick(my_buffer->tb_savedtrees[0], true);
    dbg_printf("%s\n", newick_string);
    newick_string = mfl_convert_mfl_tree_t_to_newick(my_buffer->tb_savedtrees[1], true);
    dbg_printf("%s\n", newick_string);
    newick_string = mfl_convert_mfl_tree_t_to_newick(my_buffer->tb_savedtrees[2], true);
    dbg_printf("%s\n", newick_string);
//
    //Cleaning
//    free(newick_string);
//    mfl_destroy_treebuffer(my_buffer, false);
    
    //TODO: clear the handle as well!
}

void tui_test_tree_copying(void)
{
    char cpytarget[] = "temp_examp6=[&R] (((((1,4),5),3),2),6,(7,8));";
    char* copiedtr = NULL;
    mfl_tree_t* tmplt = mfl_convert_newick_to_mfl_tree_t(cpytarget, 8);
    dbg_printf("The tree to be copied:\n %s\n", cpytarget);
    
    mfl_tree_t* cpytr = mfl_copy_tree_topology(tmplt);
    copiedtr = mfl_convert_mfl_tree_t_to_newick(cpytr, false);
    dbg_printf("The copied tree:\n %s\n", copiedtr);
    
    mfl_free_tree(cpytr);
    mfl_free_tree(tmplt);
    free(copiedtr);
    
}

void tui_test_nary_ring_creation(void)
{
    int num_taxa = 10;
    int num_desc = 2;
    mfl_tree_t* t = mfl_alloctree_with_nodes(num_taxa);
    
    mfl_node_t* ring = mfl_make_new_n_ary_ring_node(num_desc, t->treet_nodestack);
    
    
    dbg_printf("You made a ring with %i descendants\n", mfl_node_is_n_ary(ring, num_desc));
    
    dbg_printf("now put the node back:\n");
    mfl_dissolve_n_ary_ring(ring, t->treet_nodestack);
    
    return;
}

void tui_test_addition_sequence(void)
{
    
    int num_taxa = 5;
    int num_chars = 0;
    int num_og_tax = 0;
    
    char* matrix = NULL;
    
    mfl_handle_s* handle = mfl_t2s(mfl_create_handle());
    
    // Setup the handle:
    handle->n_taxa = num_taxa;
    handle->n_chars = num_chars;
    handle->n_outgroup_taxa = num_og_tax;
    handle->addseq_type = MFL_AST_ASIS;
    
    mfl_searchrec_t* searchrec = mfl_create_searchrec(handle);
    
    mfl_get_start_trees(NULL, handle, searchrec);
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
    
    dbg_printf("Testing the n-ary ring creation:\n");
    tui_test_nary_ring_creation();
    dbg_printf("\nEnd n-ary ring test\n");
    
    dbg_printf("Testing addition sequence:\n");
    tui_test_addition_sequence();
    dbg_printf("\nEnd addition sequence test\n");
    
    dbg_printf("Testing the tree printing:\n");
    tui_test_tree_printing();
    dbg_printf("\nEnd tree print test\n");
    
    
    dbg_printf("Testing bipartition setting:\n");
    tui_test_bipartition_setting();
    dbg_printf("\nEnd bipartition test\n");
    
    dbg_printf("Testing tree copying:\n");
    tui_test_tree_copying();
    dbg_printf("\nEnd tree copy test\n\n");
    
    dbg_printf("Testing branchbreaking:\n");
    tui_spr_test_environment();
    dbg_printf("\nEnd bbreak test\n");

    dbg_printf("Testing edge tables:\n");
    tui_test_edgetables();
    dbg_printf("\nEnd edge tables test\n");

    
    //-----
    
    /* Everything below here should be written into or replaced by a real test. */
    /*dbg_printf("Testing include set values:\n");
    tui_test_character_stuff();
    tui_getting_numstates_test();
    dbg_printf("\n\n"); */
    
//    dbg_printf("Test Newick stuff\n\n");
//    tui_test_newick_stuff();
    
    /*
    dbg_printf("Test treecheck stuff\n\n");
    tui_test_checktree_();
    
    dbg_printf("\n\nGoodbye!\n\n");
    tui_test_newick_stuff();
    dbg_printf("Test treecheck stuff\n\n");
    tui_test_checktree_();
    
    dbg_printf("\n\nGoodbye!\n\n");*/
    
    
	return 0;
}
