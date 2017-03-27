
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

/*
 *
 *   END MORPHY UNIT TESTS
 *
 */


int tui_getting_numstates_test(void)
{
    //int num_taxa = 5;
    //int num_chars = 16;
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
        return 1;
    }
    else {
        dbg_printf("== PASSED == %s, line: %i just enough states check\n", __FILE__, __LINE__);
        return 0;
    }

    
}


void tui_test_character_stuff(void)
{
    int num_taxa = 0;
    int num_chars = 0;
    //char subcmd1[] = "ExSet * Exclude= 1-5 8 17;";
    //char subcmd2[] = "Exclude = 1-5, 8 17;";
    //char subcmd3[] = "18-51 100";
    //char subcmd4[] = "10-15 2-6";
    
    //                        1111111
    //               1234567890123456
    char *matrix =  (char*)"120{12}000-000??001"
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
    //char temp_example_newick_for_writing3[] = "temp_examp3=[&R] (2,(6,((3,4),5),1));";
    char temp_example_newick_for_writing4[] = "temp_examp4=[&U] (2,((3,4),(5,20),1));"; // Polytomy and multi-digit tip number not in sequence
    char temp_example_newick_for_writing5[] = "temp_examp5=[&R] (((((1,4),5),3),2),6);";
    char temp_example_newick_for_writing6[] = "temp_examp6=[&R] (((((1,4),5),3),2),6,(7,8));";
    //char temp_example_newick_for_writing7[] = "temp_examp7=[&U] ((1000,856),(2,3),(56,4));";
    //char temp_example_newick_for_writing8[] = "temp_examp8=[&R] (1,(2,(3,(4,5))));";
    
    char *sample_newick = NULL;
    
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

void tui_spr_test_environment(void)
{
    // NOTE: This is a tui function and should be moved there.
    
    int temp_max_trees = 50000;
    
    mfl_handle_s* testhandle = mfl_t2s(mfl_create_handle());
    
    char* cliptesttree = NULL;
    
    cliptesttree = (char*)"tree1=[&R] (1,(2,(((((((((((((((((((((3,39),12),(11,(53,64))),30),(42,62)),48),(25,32)),74),21),((((((6,61),76),17),67),(8,45)),((((((((((14,22),38),(16,18)),((37,58),75)),(59,73)),15),26),68),(51,56)),36))),((((13,(40,((46,55),54))),49),((((((29,34),(33,63)),72),57),65),35)),23)),70),44),27),(31,43)),(((9,((19,41),(20,28))),24),47)),71),((4,10),69)),((50,78),52)),7),(((5,66),77),60))));";
    //cliptesttree = (char*)"tree1=[&R] (1,((2,(79,(80,((81,((82,88),(86,87))),(83,(84,85)))))),(((((((((((((((((((((3,39),12),(11,(53,64))),30),(42,62)),48),(25,32)),74),21),((((((6,61),76),17),67),(8,45)),((((((((((14,22),38),(16,18)),((37,58),75)),(59,73)),15),26),68),(51,56)),36))),((((13,(40,((46,55),54))),49),((((((29,34),(33,63)),72),57),65),35)),23)),70),44),27),(31,43)),(((9,((19,41),(20,28))),24),47)),71),((4,10),69)),((50,78),52)),7),(((5,(66,((89,90),((91,((92,98),(96,97))),(93,(94,(95,(99,100)))))))),77),60))));";
    //cliptesttree = (char*)"594=[&R] (((((1,4),5),3),2),(6,(7,((8,9),(10,(11,12))))));";
    //    cliptesttree = (char*)"tree_857 = [&R] (1,((4,(2,((6,(7,((8,9),(10,(11,12))))),3))),5));";
    //cliptesttree = (char*)"tree_838 = [&R] ((((12,(((((3,2),6),7),(8,9)),10)),11),1),(4,5));";
    //cliptesttree = (char*)"tree_795 = [&R] (1,(((((7,11),(10,((2,((4,8),6)),3))),5),9),12));";
    //cliptesttree = (char*)"temp_examp6=[&R] (1,(4,(5,(3,(2,(6,(7,(8,(9,(10,(11,12)))))))))));";
    cliptesttree = (char*)"sixonefour = [&R] (1,(4,(5,(2,(6,(7,(8,(9,((10,3),(11,12))))))))));";
    //cliptesttree = (char*)"fivethirty = [&R] (1,((2,((6,((7,4),((8,9),(10,(11,12))))),3)),5));";
    //cliptesttree = (char*)"temp_examp6=[&R] ((1,2),(3,4));";
    //cliptesttree = (char*)"temp_examp6=[&R] ((1,2),3);";
    //cliptesttree = (char*)"temp_examp6=[&R] ((1,2),(3,(4,5)));";
    //cliptesttree = (char*)"temp_examp6=[&R] (1,(2,(3,(4,(5,6)))));";
    //cliptesttree = (char*)"temp_examp6=[&R] (5,(4,(3,(2,1))));";
    //cliptesttree = (char*)"temp_examp6=[&R] ((1,(2,(6,7))),(3,(4,5)));";
    //cliptesttree = (char*)"temp_examp6=[&R] ((1,((2,8),(6,7))),(3,(4,5)));";
    //cliptesttree = "UNTITLED = [&R] ((((((((1,2),3),((((4,5),(((((6,7),8),(9,10)),(11,12)),((13,(14,15)),(16,17)))),18),((19,20),21))),((((((22,23),(24,25)),((((((26,27),28),29),30),(31,(32,33))),((34,35),(((36,37),((38,39),40)),41)))),(42,((43,((44,45),46)),(47,(48,(49,50)))))),51),((((52,53),54),((55,56),57)),(58,((59,60),61))))),(((62,((63,64),(65,66))),(67,(68,(69,(70,71))))),(((72,73),((74,(((75,76),77),(78,79))),80)),81))),((((82,83),(((84,85),86),(87,(88,89)))),(90,91)),((92,93),((94,(95,((96,97),98))),(99,100))))),(((((((101,102),(103,(((104,105),106),107))),(((108,109),(110,111)),(((112,(113,114)),(115,116)),(117,(118,119))))),(((120,((121,((122,(((123,((124,125),126)),127),(128,129))),(130,(131,132)))),(((133,((134,135),(136,((137,138),139)))),(((140,141),142),(((143,144),145),((146,147),148)))),(149,(((150,151),152),153))))),((154,155),156)),((157,(158,(159,160))),((((161,162),163),(164,165)),((166,167),168))))),(((((169,(170,171)),172),((173,174),175)),(((176,177),(178,(179,((180,181),(182,183))))),184)),(((185,186),((187,(188,189)),190)),((((((191,192),193),(((194,195),196),(197,198))),(199,((200,201),(202,203)))),204),(((205,206),207),(208,(209,210))))))),(((((211,212),(((213,214),(215,216)),(((217,218),((219,220),221)),222))),((223,224),(((225,226),(227,228)),(((229,(230,((231,232),233))),(((234,235),(236,237)),238)),(((239,240),((241,242),243)),(244,245)))))),(((((246,247),248),249),250),251)),((((((252,253),(254,255)),((256,((257,258),259)),(260,261))),(((262,263),(264,((265,266),(267,268)))),((((269,(270,271)),((272,273),(((274,275),276),277))),((278,(279,(280,281))),(((282,283),(284,285)),286))),(287,288)))),(((289,((290,((291,292),293)),(294,295))),(296,(297,((298,299),(300,301))))),(((302,303),304),((305,((306,307),((308,(309,310)),(311,312)))),(((313,(314,315)),316),317))))),((318,319),(((320,321),(322,323)),((324,325),(326,327))))))),(328,((((329,330),331),332),(333,334))))),(((((335,336),((337,((338,339),(340,341))),((342,343),344))),(((((345,346),347),348),349),(350,(351,352)))),(((353,354),((((355,356),357),((358,(359,360)),(361,362))),(((363,((364,365),366)),(((367,368),369),(370,371))),((((372,(373,374)),(375,376)),(((377,(378,379)),(380,(381,(382,383)))),384)),(((385,386),387),((388,389),390)))))),((((((391,392),(393,((394,395),((396,397),398)))),(399,400)),(((401,(((402,403),404),(405,(406,407)))),((((408,409),(410,(411,(412,413)))),414),415)),(((416,((417,(418,419)),420)),(((421,422),423),(424,425))),(((426,427),428),((429,430),431))))),(((432,433),((434,435),((436,437),438))),((439,(440,441)),(442,(443,(((444,445),(446,447)),448)))))),((449,450),(451,452))))),(((((453,((454,455),456)),(457,458)),((((((459,460),461),((462,(((463,(464,465)),(466,467)),468)),(469,(470,(471,472))))),((473,(((474,475),476),477)),((478,(((479,480),481),482)),((((483,484),485),(486,487)),(((488,489),(490,491)),((492,493),((494,495),(496,497)))))))),(((498,(499,((500,501),502))),(((503,504),(505,506)),(((507,508),509),(((510,511),(512,513)),(514,515))))),((516,517),((518,519),(520,((521,522),523)))))),((524,(525,526)),((((527,((528,(529,530)),531)),532),(((533,(534,(535,(536,537)))),(((((538,539),540),541),((542,543),(((544,(545,546)),547),548))),549)),((((550,551),552),(553,(554,555))),((556,557),558)))),(((((559,560),561),(562,(563,564))),((565,566),567)),((568,569),(570,(571,(572,573))))))))),((((((574,((575,(576,577)),578)),(((579,580),581),582)),((583,584),(585,(586,587)))),(((((588,((589,(590,(591,592))),(593,(594,(595,596))))),(597,(598,599))),(((600,601),((((602,603),604),605),(606,607))),((608,(609,610)),611))),(((612,(613,614)),615),(616,(((617,(618,619)),(620,(621,(622,623)))),(((624,625),626),627))))),((628,(629,(630,631))),632))),((633,((634,635),(636,637))),(((((638,((639,640),641)),(((642,643),644),(((645,(646,647)),(648,649)),650))),(651,(((652,653),654),(655,656)))),(657,((658,659),660))),(661,662)))),((((((((663,664),665),(666,667)),(((668,(669,670)),(671,672)),(673,674))),(675,676)),((((((677,678),679),((680,((681,(682,683)),((684,685),686))),((687,(688,689)),(690,691)))),((692,693),694)),(((695,696),((697,(698,699)),(700,(701,702)))),(703,((704,(705,706)),707)))),708)),(((((((709,710),711),(((712,(713,714)),715),((716,(717,718)),(719,(720,721))))),(((722,(723,724)),((725,726),727)),728)),(((729,(730,(((731,(732,733)),734),735))),((736,((737,(738,739)),((740,741),742))),((743,(744,745)),((746,747),748)))),((((749,(750,751)),((752,753),(754,(755,(((756,757),(758,759)),760))))),761),(762,763)))),((((764,765),(766,767)),((((768,769),(770,771)),(772,773)),((774,775),776))),(((777,778),((779,((780,781),782)),(783,784))),((785,((786,787),((788,789),790))),((((791,(792,793)),794),((795,796),797)),(798,799)))))),((((((800,(801,802)),803),804),(805,((806,807),(808,809)))),((810,811),((812,((813,814),(815,(816,817)))),818))),(((819,820),((821,(822,(823,824))),(((825,826),827),(828,829)))),((((830,831),(832,(833,834))),((835,836),((837,(838,839)),(840,841)))),(842,843)))))),((844,845),((((846,((847,(848,849)),(((850,(851,852)),(853,854)),((855,((856,857),858)),859)))),(((860,((861,862),(((863,(864,(865,866))),867),(868,869)))),(((870,871),(872,(873,(874,(875,876))))),(877,878))),(879,880))),((((881,882),(883,(884,885))),((886,887),(888,889))),890)),((891,(892,893)),(((894,895),896),(((897,(898,899)),900),901)))))))),((((((902,903),904),905),(906,((907,(908,909)),((910,911),(912,(913,914)))))),(((((915,916),917),918),((919,(((920,(921,922)),(923,924)),925)),(926,(927,928)))),(929,(930,931)))),((((932,(933,(934,935))),936),((937,938),939)),((((940,941),(((942,943),((944,945),946)),(947,(948,((949,950),951))))),((((952,953),((954,(955,((956,957),(958,959)))),960)),(961,962)),(((963,964),(965,966)),((((967,968),(969,970)),(971,972)),973)))),((((974,975),(976,977)),(((978,979),((980,981),(982,983))),((984,(985,986)),((987,((988,989),(990,991))),(992,993))))),((994,995),(996,((997,(998,999)),1000))))))))));";
    
    
    mfl_tree_t* testree = mfl_convert_newick_to_mfl_tree_t(cliptesttree, 0);
    
    char* treedraw = mfl_drawtree(testree);
    dbg_printf("%s", treedraw);
    
    testhandle->n_taxa = testree->treet_num_taxa;
    
    mfl_searchrec_t* searchrec = mfl_create_searchrec(testhandle);
    
    mfl_unroot_tree(testree);
    
    if (testree->treet_root) {
        searchrec->sr_swap_entry = testree->treet_root;
    } else {
        searchrec->sr_swap_entry = testree->treet_start;
    }
    
    searchrec->sr_treebuffer = mfl_alloc_treebuffer(temp_max_trees);
    searchrec->sr_swaping_on = testree;
    searchrec->sr_handle_ptr = testhandle;
    //searchrec->sr_bswaptype = MFL_BST_SPR;
    
    mfl_append_tree_to_treebuffer(testree, searchrec->sr_treebuffer, testhandle);
    
    time_t before = 0;
    time_t after = 0;
    
    //mfl_assign_nodeweights(searchrec->sr_swap_entry);
    
    before = clock();
    mfl_pruning_traversal(searchrec->sr_swap_entry, searchrec);
    after = clock();
    
    dbg_printf("\nNumber of rearrangements attempted for %i taxa: %lli\n", testree->treet_num_taxa, searchrec->sr_rearrangement_counter);
    //    dbg_printf("\nAnd here are the trees:\n");
    dbg_printf("Time: %f\n\n", (float)(after-before)/CLOCKS_PER_SEC);
    
    char* newicktr = NULL;
    int i = 0;
    for (i = 0; i < searchrec->sr_treebuffer->tb_num_trees; ++i) {
        newicktr = mfl_convert_mfl_tree_t_to_newick(searchrec->sr_treebuffer->tb_savedtrees[i], 0);
        dbg_printf("TREE tree_%i = ", (i+1));
        dbg_printf("%s\n", newicktr);
        free(newicktr);
    }
    
    mfl_free_tree(testree);
    
}

// TODO: Move to utilities
bool tui_check_values(int actual, int expected)
{
    
    
    if (actual == expected) {
        return false;
    } else {
        dbg_printf("===================\n");
        dbg_printf(" Expected:   %i\n", expected);
        dbg_printf(" Calculated: %i\n", actual);
        dbg_printf("===================\n");
        return true;
    }
}

int tui_test_treelength_calculation(mfl_handle_s* handle, int* expectedlengths)
{
    int i = 0;
    int numtax = 0;
    int numintrees = 0;
    int testfails = 0;
    char *intree = NULL;
    mfl_tree_t* t = NULL;
    mfl_partition_set_t* dataparts = NULL;

    
    numtax = handle->n_taxa;
    numintrees = handle->n_input_newick_trees;
    
    // Set up the input data
    dataparts = mfl_generate_search_data(handle);
    
    for (i = 0; i < numintrees; ++i) {
        // Create the tree
        intree = handle->input_newick_trees[i];
        t = mfl_convert_newick_to_mfl_tree_t(intree, numtax);
        t->treet_parsimonylength = 0;
        
        // Prepare the data and tree
        mfl_setup_input_tree_with_node_data(t, dataparts);
        
        // Output the length
        mfl_fullpass_tree_optimisation(t, dataparts);
        
        // Compare the length to the expected value (if given)
        if (tui_check_values(t->treet_parsimonylength, expectedlengths[i])) {
            ++testfails;
        }
        
        // Destroy the tree
        //tui_check_broken_tree(t, false);
        mfl_free_tree(t);
        t = NULL;
    }
    
    // Destroy the input data
    
    // Summarise the test
    dbg_printf("\n");
    dbg_printf("SUMMARY:\n");
    dbg_printf("Test of all topologies for matrix %s concluded\n", handle->input_data);
    if (testfails) {
        dbg_printf("Test FAILED %i times.\n", testfails);
    }
    else {
        dbg_printf("All tests PASSED\n");
    }
    dbg_printf("\n");
    
    return testfails;
}

void tui_test_basic_character_optimisation(void)
{
    int num_taxa = 12; //12;//78;
    int num_chars = 1;//236;
    int num_og_tax = 0;
                   //.........111
                   //123456789012

char matrix[] = "023-???1--32;";

//    "10300000000000000000"
//    "33-00000000000000000"
//    "0--40511110000000000"
//    "0--40511111111000000"
//    "0-300011111111111000;";
    //"5315-35135-151--;";//"43132525--4--4--;";//"03422--50555----;";//"0402153-1-50-40-";//////"0100?011-1?0-???10110-???1?111-000;";//"-330355233-131-5";//"4-24-330-11305-2;";////"05--344-1000205-;";//
     char* testnewick;
    testnewick = (char*)"[&R] (1,(2,(3,(4,(5,(6,(7,(8,(9,(10,(11,12)))))))))));";
    
    //dbg_printf("Here is the data matrix:\n%s\n", matrix);
    //dbg_printf("And here is the tree:\n%s\n", testnewick);
    
    time_t before = 0;
    time_t after = 0;
    
    mfl_handle_s* handle = mfl_t2s(mfl_create_handle());
    
    // Setup the handle:
    handle->n_taxa = num_taxa;
    handle->n_chars = num_chars;
    handle->n_outgroup_taxa = num_og_tax;
    handle->n_to_hold = 3;
    handle->input_data = matrix;
    handle->addseq_type = MFL_AST_ASIS;
    //handle->gap_method = MFL_GAP_MISSING_DATA;
    //    handle->addseq_type = MFL_AST_RANDOM;
    
    mfl_searchrec_t* searchrec = mfl_create_searchrec(handle);
    
    before = clock();
    mfl_partition_set_t* dataparts = mfl_generate_search_data(handle);
    after = clock();
    dbg_printf("Time to process matrix: %f\n", (float)(after-before)/CLOCKS_PER_SEC);
    
    mfl_tree_t* testtree = mfl_convert_newick_to_mfl_tree_t(testnewick, num_taxa);
    
    mfl_setup_input_tree_with_node_data(testtree, dataparts);
    
    tui_check_broken_tree(testtree, false);
    
    before = clock();
    //mfl_postorder_traversal(testtree->treet_root, &testtree->treet_parsimonylength);
    mfl_fullpass_tree_optimisation(testtree, dataparts);
    after = clock();
    
    dbg_printf("And here is it's length: %i\n", testtree->treet_parsimonylength);
    
    dbg_printf("Time: %f\n\n", (float)(after-before)/CLOCKS_PER_SEC);
    tui_check_broken_tree(testtree, false);
    
    //
    mfl_free_tree(testtree);
    
}


void tui_test_local_reoptimisation(void)
{
    /*
     * A seriously crufty test of LOCAL REOPTIMISATION, AT PRESENT
     */
    
    int num_taxa = 7;
    int num_chars = 20;
    int num_og_tax = 0;
    int num_trees = 4;
    int pass = 0;
    int fail = 0;
    
    char *testnwk[] = { (char*)"[&R] (1,(2,(3,(4,(5,(6,7))))));",
                        (char*)"[&R] (1,(2,(7,((6,5),(4,3)))));",
                        (char*)"[&R] (1,(2,((6,(7,4)),(5,3))));",
                        (char*)"[&R] (1,(2,(7,(6,(5,(4,3))))));"};
    
    //                        1         2
    //               1...5....0....5....0
    char matrix[] =
//    "10000000000000000000"
//    "-3--0000000000000000"
//    "----0511110000000000"
//    "----0511111111000000"
//    "--2-0011111111111000"
//    "-2300011111111111111"
//    "-1100011111111111111;";
    
//    "10000000000000000000"
//        "-3-00000000000000000"
//        "---40511110000000000"
//        "---40511111111000000"
//        "--200011111111111000"
//        "-2300011111111111111"
//        "-1100011111111111111;";
    "1030000000000000000-"
    "33--0000-0-0000-----"
    "--------------------"
    "----0511-1--------00"
    "--3-----1111111110--"
    "3330----1-1-11-11111"
    "11100011---111111111;";
    
//    "103000"
//    "33-000"
//    "---405"
//    "---405"
//    "--3000"
//    "333000"
//    "111000;";
//    
    mfl_handle_s* handle = mfl_t2s(mfl_create_handle());
    
    // Setup the handle:
    handle->n_taxa = num_taxa;
    handle->n_chars = num_chars;
    handle->n_outgroup_taxa = num_og_tax;
    handle->n_to_hold = 3;
    handle->input_data = matrix;
//    handle->gap_method = MFL_GAP_MISSING_DATA;
//    handle->addseq_type = MFL_AST_ASIS;
//    handle->addseq_type = MFL_AST_RANDOM;
    
    //mfl_searchrec_t* searchrec = mfl_create_searchrec(handle);
    mfl_partition_set_t* dataparts = mfl_generate_search_data(handle);
    
    for (int i = 0; i < num_trees; ++i) {
        
        dbg_printf("Testing local reopt in tree: %i\n", i);
        dbg_printf("===============================\n\n");
        mfl_tree_t* testtree = mfl_convert_newick_to_mfl_tree_t(testnwk[i], num_taxa);
        
        mfl_setup_input_tree_with_node_data(testtree, dataparts);
        
        mfl_fullpass_tree_optimisation(testtree, dataparts);
        
        
        int oldlen = testtree->treet_parsimonylength;
        
        // Prune tip 3;
        for (int j = 2; j < testtree->treet_num_taxa; ++j) {
            dbg_printf("Removing branch %i in tree %i:\n", j + 1, i);
            dbg_printf("==============================\n");
            dbg_printf("==============================\n");
            dbg_printf("%s\n\n", testnwk[i]);
            int diff = 0;
            mfl_cliprec_t clip;
            mfl_node_t* tgt1 = testtree->treet_treenodes[j]->nodet_edge->nodet_next->nodet_edge;
//            mfl_node_t* tgt2 = testtree->treet_treenodes[j]->nodet_edge->nodet_next->nodet_next->nodet_edge;
            
            testtree->treet_treenodes[j]->nodet_edge->nodet_weight = 3;
            testtree->treet_treenodes[j]->nodet_edge->nodet_downpass_visited = true;
            
            mfl_clip_branch(testtree->treet_treenodes[j], &clip);
            
            testtree->treet_parsimonylength = 0;
            
            mfl_fullpass_tree_optimisation(testtree, dataparts);
            
            mfl_local_add_cost(testtree->treet_treenodes[j], tgt1, -1, &diff);
            
            // Compare length and sum of lengths;
            printf("\nDiff: %i\n\n", diff);
            if (oldlen == (testtree->treet_parsimonylength + diff)) {
                dbg_ppass("lengths are equal.");
                ++pass;
            }
            else {
                dbg_pfail("lengths do not sum to full length as expected.\n");
                dbg_printf("Estimated: %i, directly calculated: %i\n", oldlen, (testtree->treet_parsimonylength + diff));
                ++fail;
            }
            dbg_printf("\n");
            mfl_restore_branching(&clip);
        }
        
        mfl_free_tree(testtree);
    }
    
    if (fail) {
        dbg_pfail("Test of local reoptimisation failed ");
        dbg_printf("%i times out of %i.\n", fail, num_trees * (num_taxa-2));
    }
    else {
        dbg_ppass("all tests of local reoptimisation succeeded\n");
    }
    
    return;
}


void tui_test_compare_replace(void)
{
    int i = 0;
    int num_taxa = 7;
    int num_chars = 20;
    int num_og_tax = 1;
    int num_trees = 20;
    
    //                        1         2
    //               1...5....0....5....0
    char matrix[] =     "--111111111111111110"
                        "--111111111111111110"
                        "--111111111111111100"
                        "--111111111111111101"
                        "-1111111111111111101"
                        "11111111111111111101"
                        "11000000000000000001;";

//    "--111111111111111111"
//    "--111111111111111111"
//    "--111111111111111111"
//    "--111111111111111111"
//    "--111111111111111111"
//    "--111111111111111111"
//    "11000000000000000000;";
    
//    "10300000000000000000"
//                     "33-00000000000000000"
//                     "---40511110000000000"
//                     "---40511111111000000"
//                     "--300011111111111000"
//                     "33300011111111111111"
//                     "11100011111111111111;";
    

    

    
    char *testnwk[] = { (char*)"[&R] (1,(2,(3,(4,(5,(6,7))))));",
                        (char*)"[&R] (1,(2,(3,(4,(5,(6,7))))));",
                        (char*)"[&R] (1,(2,(3,(4,(5,(6,7))))));",
                        (char*)"[&R] (1,(2,(3,(4,(5,(6,7))))));",
                        (char*)"[&R] (1,(2,(3,(4,(5,(6,7))))));",};
//    { (char*)"[&R] (1,(5,(3,(4,(2,(6,7))))));",
//                        (char*)"[&R] (1,((2,3),(4,(5,(6,7)))));",
//                        (char*)"[&R] (1,(7,((5,((6,3),2)),4)));",
//                        (char*)"[&R] (1,((7,5),(2,(3,(6,4)))));",
//                        (char*)"[&R] (7,(2,(3,(4,(5,(6,1))))));",
//                        (char*)"[&R] (6,(2,(3,(4,(5,(1,7))))));",
//                        (char*)"[&R] (5,(2,((3,7),(1,(6,4)))));",
//                        (char*)"[&R] (4,(2,((3,7),(1,(6,5)))));",
//                        (char*)"[&R] (4,(2,((3,7),(1,(6,5)))));",
//                        (char*)"[&R] (5,(2,((3,7),(1,(6,4)))));",
//                        (char*)"[&R] (1,((7,5),(2,(3,(6,4)))));" };
    
    char *bestnwk = (char*) "[&R] (1,(2,(3,(4,(5,(6,7))))));";
    
    mfl_tree_t* besttr = mfl_convert_newick_to_mfl_tree_t(bestnwk, num_taxa);
    besttr->treet_edges = (mfl_nodearray_t)mfl_malloc(besttr->treet_num_nodes * sizeof(mfl_node_t*), 0);
    
    
//    mfl_convert_from_stored_topol(besttr, besttr);
//    char *reprint = mfl_convert_mfl_tree_t_to_newick(besttr, false);
//    dbg_printf("After converting original tree to itself:\n");
//    dbg_printf("%s\n\n", reprint);
    
    
    mfl_handle_s* handle = mfl_t2s(mfl_create_handle());
    
    // Setup the handle:
    handle->n_taxa = num_taxa;
    handle->n_chars = num_chars;
    handle->n_outgroup_taxa = num_og_tax;
    handle->n_to_hold = num_trees;
    handle->input_data = matrix;
    //    handle->gap_method = MFL_GAP_MISSING_DATA;
    handle->addseq_type = MFL_AST_ASIS;
    //    handle->addseq_type = MFL_AST_RANDOM;
    
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, handle->rseed);
    
    mfl_searchrec_t* searchrec = mfl_create_searchrec(handle);
    
    
    // Set up the searchrec
    searchrec->sr_random_number = r;
    
    mfl_partition_set_t* dataparts = mfl_generate_search_data(handle);
    mfl_treebuffer_t* trbuf = mfl_alloc_treebuffer(num_trees);
    mfl_stepwise_addition_t *sarec = mfl_generate_stepwise_addition(besttr, handle, searchrec);
    
    for (i = 0; i < num_trees; ++i) {
        mfl_append_tree_to_treebuffer(mfl_convert_newick_to_mfl_tree_t(bestnwk, num_taxa), trbuf, handle);
    }
    
    mfl_setup_input_tree_with_node_data(besttr, dataparts);
    mfl_fullpass_tree_optimisation(besttr, dataparts);
    mfl_update_stored_topology(besttr, besttr);
    dbg_printf("The trees after storing and pushing to buffer:\n");
    char *nwk = NULL;
    
    for (i = 0; i < trbuf->tb_num_trees; ++i) {
        
        trbuf->tb_savedtrees[i]->treet_edges = (mfl_nodearray_t)mfl_malloc(trbuf->tb_savedtrees[i]->treet_num_nodes * sizeof(mfl_node_t*), 0);
        
        mfl_shuffle_tree(trbuf->tb_savedtrees[i], searchrec, 10);
        
        nwk = mfl_convert_mfl_tree_t_to_newick(trbuf->tb_savedtrees[i], false);
        
        mfl_update_stored_topology(trbuf->tb_savedtrees[i], trbuf->tb_savedtrees[i]);
        mfl_setup_input_tree_with_node_data(trbuf->tb_savedtrees[i], dataparts);
        mfl_fullpass_tree_optimisation(trbuf->tb_savedtrees[i], dataparts);
        
        dbg_printf("%s   length: %i\n", nwk, trbuf->tb_savedtrees[i]->treet_parsimonylength);
        free(nwk);
        ++sarec->stpadd_num_held;
    }
    
    mfl_attempt_replacement_stepadd(trbuf, besttr, sarec);

    
    dbg_printf("\nThe trees after searching and replacing:\n");
    for (i = 0; i < trbuf->tb_num_trees; ++i) {
        mfl_convert_from_stored_topol(trbuf->tb_savedtrees[i], trbuf->tb_savedtrees[i]);
        nwk = mfl_convert_mfl_tree_t_to_newick(trbuf->tb_savedtrees[i], false);
        dbg_printf("%s   length: %i\n", nwk, trbuf->tb_savedtrees[i]->treet_parsimonylength);
        free(nwk);
    }
    
    //trbuf->tb_savedtrees[10]->treet_parsimonylength = 0;
    //mfl_fullpass_tree_optimisation(trbuf->tb_savedtrees[10], dataparts);
}

void tui_test_stepwise_addition(void)
{
    int num_taxa = 7;
    int num_chars = 20;
    int num_og_tax = 1;
    
    //                        1         2
    //               1...5....0....5....0
    char matrix[] =
//    "031000"
//                    "3-3000"
//                    "---1-0"
//                    "---110"
//                    "-3---0"
//                    "333--1"
//                    "111-11;";

//    "10000000000000000000"
//    "-3-00000000000000000"
//    "---40511110000000000"
//    "---40511111111000000"
//    "--200011111111111000"
//    "-2300011111111111111"
//    "-1100011111111111111;";

//        "10300000000000000000"
//        "33300000000000000000"
//        "0-340511110000000000"
//        "0-340511111111000000"
//        "0-300011111111111000"
//        "33300011111111111111"
//        "11100011111111111111;";

//    "10300012101110100001"
//    "33210000101000022101"
//    "00101010001111001111"
//    "11020511211102102030"
//    "00311211111111111032"
//    "33303500121211211111"
//    "11100011112111100111;";
    
//    "10300000000000000000"
//    "33-00000000000000000"
//    "02-40511110000000000"
//    "02-40511111111000000"
//    "02300011111111111000"
//    "33300011111111111111"
//    "11100011111111111111;";
//
//                        1         2
//               1...5....0....5....0
//                "10300000000000000000"
//                "33-00000000000000000"
//                "0--40511110000000000"
//                "0--40511111111000000"
//                "0-300011111111111000"
//                "33300011111111111111"
//                "11100011111111111111;";
    "1030000000000000000-"
    "33--0000-0-0000-----"
    "--------------------"
    "----0511-1--------00"
    "--3-----1111111110--"
    "3330----1-1-11-11111"
    "11100011---111111111;";
    
//    "031000"
//                    "3-3000"
//                    "--0100"
//                    "--0110"
//                    "-30110"
//                    "333111"
//                    "111111";
//

    mfl_handle_s* handle = mfl_t2s(mfl_create_handle());
    
    // Setup the handle:
    handle->n_taxa = num_taxa;
    handle->n_chars = num_chars;
    handle->n_outgroup_taxa = num_og_tax;
    handle->n_to_hold = 10;
    handle->input_data = matrix;
//    handle->gap_method = MFL_GAP_MISSING_DATA;
    handle->addseq_type = MFL_AST_ASIS;
//    handle->addseq_type = MFL_AST_RANDOM;
    
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, handle->rseed);
    
    mfl_searchrec_t* searchrec = mfl_create_searchrec(handle);
    mfl_partition_set_t* dataparts = mfl_generate_search_data(handle);
    searchrec->sr_random_number = r;
    // TODO: Needs to pass return to variable
    // TODO: Return type needs designing
    /**/ mfl_get_start_trees(dataparts, handle, searchrec);
    
}


void tui_test_edgetables(void)
{
    char* cliptesttree1 = NULL;
    char* cliptesttree2 = NULL;
    bool compare = false;
    
    //cliptesttree = (char*)"temp_examp6=[&U] ((1,2),(3,4));";
    //cliptesttree = (char*)"temp_examp6=[&U] ((1,2),(3,(4,5)));";
    //cliptesttree = (char*)"temp_examp6=[&U] (1,(2,(3,(4,(5,6)))));";
    //cliptesttree = (char*)"temp_examp6=[&U] (5,(4,(3,(2,1))));";
    //cliptesttree = (char*)"temp_examp6=[&U] ((1,(2,(6,7))),(3,(4,5)));";
    cliptesttree1 = (char*)"equal_test=[&U] ((1,(2,3)), (4,(5,6)));";
    //cliptesttree = (char*)"equal_test=[&U] ((4,(5,6)), (1,(2,3)));";
    //cliptesttree1 = (char*)"equal_test=[&U] (2, ((4,7), ((1,(3,5)), (8,(6,9)))));";
    cliptesttree2 = (char*)"equal_test=[&U] (2, ((4,7), ((8,(6,9)), (1,(3,5)))));";
    //cliptesttree = (char*)"equal_test=[&U] (2, (((8,(6,9)), (1,(3,5))), (4,7)));";
    //cliptesttree = (char*)"tree1=[&U] (1,(2,(((((((((((((((((((((3,39),12),(11,(53,64))),30),(42,62)),48),(25,32)),74),21),((((((6,61),76),17),67),(8,45)),((((((((((14,22),38),(16,18)),((37,58),75)),(59,73)),15),26),68),(51,56)),36))),((((13,(40,((46,55),54))),49),((((((29,34),(33,63)),72),57),65),35)),23)),70),44),27),(31,43)),(((9,((19,41),(20,28))),24),47)),71),((4,10),69)),((50,78),52)),7),(((5,66),77),60))));";
    
    mfl_edgetable_t* test_edgetable1 = NULL;
    mfl_edgetable_t* test_edgetable2 = NULL;
    mfl_tree_t* testree1 = mfl_convert_newick_to_mfl_tree_t(cliptesttree1, 0);
    mfl_tree_t* testree2 = mfl_convert_newick_to_mfl_tree_t(cliptesttree2, 0);
    
    if(!testree1->treet_root) {
        mfl_assign_bottom_node(testree1->treet_start);
        test_edgetable1 = mfl_initiate_edgetable_t(testree1->treet_num_taxa, 0);
    } else {
        mfl_assign_bottom_node(testree1->treet_root);
        test_edgetable1 = mfl_initiate_edgetable_t(testree1->treet_num_taxa, 1);
    }
    
    if(!testree2->treet_root) {
        mfl_assign_bottom_node(testree2->treet_start);
        test_edgetable2 = mfl_initiate_edgetable_t(testree2->treet_num_taxa, 0);
    } else {
        mfl_assign_bottom_node(testree2->treet_root);
        test_edgetable2 = mfl_initiate_edgetable_t(testree2->treet_num_taxa, 1);
    }
    
    
    mfl_get_edgetable(test_edgetable1, testree1);
    mfl_get_edgetable(test_edgetable2, testree2);
    
    tui_print_edgetable(test_edgetable1);
    tui_print_edgetable(test_edgetable2);
    
    compare = mfl_compare_edge_tables(test_edgetable1, test_edgetable2);
    if(compare == true){
        dbg_printf("Trees are the same!\n");
    } else {
        dbg_printf("Trees are different!\n");
    }
    
    
    //
    //free(test_edgetable);
    mfl_destroy_edgetable(test_edgetable2);
    mfl_free_tree(testree2);
    
    // Destroy the table
    mfl_destroy_edgetable(test_edgetable1);
    mfl_free_tree(testree1);
}
//
//void tui_test_bipartition_tables(void)
//{
//    char* cliptesttree1 = NULL;
//    char* cliptesttree2 = NULL;
//    char* testtree1 = NULL;
//    char* testtree2 = NULL;
//    char* testtree3 = NULL;
//    char* testtree4 = NULL;
//    char* testtree5 = NULL;
//    
//    cliptesttree1 = (char*)"temp_examp6=[&R] ((1,2),(3,4));";
//    cliptesttree2 = (char*)"temp_examp6=[&R] ((1,3),(2,4));";
//    testtree1 = (char*)"testtree1 = [&R] ((((1,2),3),((4,5),6)),(((7,8),((9,(10,11)),12)),(13,(14,15))));";
//    testtree2 = (char*)"testtree2 = [&R] ((((4,5),6),((1,2),3)),(((10,11),((9,(7,8)),12)),(13,(14,15))));";
//    testtree3 = (char*)"testtree3 = [&R] ((((13,(14,15),((4,5),6)),(((10,11),((9,(7,8)),12)),(1,2),3))));";
//    testtree4 = (char*)"testtree4 = [&R] ((((1,2),3),(13,(14,15))),(((7,8),(((4,5),6),12)),(9,(10,11))));";
//    testtree5 = (char*)"testtree5 = [&R] ((((13,2),7),((14,6),11)),(((3,9),((8,(10,5)),12)),(1,(4,15))));";
//    
//    //Get the trees into mfl_tree_t
//    mfl_tree_t* simpletree1 = mfl_convert_newick_to_mfl_tree_t(cliptesttree1, 0);
//    mfl_tree_t* simpletree2 = mfl_convert_newick_to_mfl_tree_t(cliptesttree2, 0);
//    mfl_tree_t* testree1 = mfl_convert_newick_to_mfl_tree_t(testtree1, 0);
//    mfl_tree_t* testree2 = mfl_convert_newick_to_mfl_tree_t(testtree2, 0);
//    mfl_tree_t* testree3 = mfl_convert_newick_to_mfl_tree_t(testtree3, 0);
//    mfl_tree_t* testree4 = mfl_convert_newick_to_mfl_tree_t(testtree4, 0);
//    mfl_tree_t* testree5 = mfl_convert_newick_to_mfl_tree_t(testtree5, 0);
//    
//    //Setting the biparititions
//    mfl_set_bipartitions(simpletree1->treet_root);
//    mfl_set_bipartitions(simpletree2->treet_root);
//    mfl_set_bipartitions(testree1->treet_root);
//    mfl_set_bipartitions(testree2->treet_root);
//    mfl_set_bipartitions(testree3->treet_root);
//    mfl_set_bipartitions(testree4->treet_root);
//    mfl_set_bipartitions(testree5->treet_root);
//    
////    Initialising the tables
//    mfl_bipartition_table* simple_bipar_table = NULL;
//    simple_bipar_table = mfl_initialise_bipartition_table(simpletree1->treet_num_taxa);
//    
//    //Simple example (2 4 taxa trees)
//    //Generating the bipartition table for the first tree
//    mfl_get_bipartition_traversal(simpletree1->treet_root, simple_bipar_table);
//    mfl_get_bipartition_traversal(simpletree2->treet_root, simple_bipar_table);
//    tui_print_bipartition_tables(simple_bipar_table);
//
//    //Simple example 2 (5 15 taxa trees)
//    //Initialising the tables
////    mfl_bipartition_table* simple_bipar_table2 = NULL;
////    simple_bipar_table2 = mfl_initialise_bipartition_table(testree1->treet_num_taxa);
//    
//    //Initialising the tables
//    mfl_bipartition_table* comple_bipar_table = NULL;
//    comple_bipar_table = mfl_initialise_bipartition_table(testree1->treet_num_taxa);
//    
//    //Simple example (2 4 taxa trees)
//    //Generating the bipartition table for the first tree
//    mfl_get_bipartition_traversal(testree1->treet_root, comple_bipar_table);
//    mfl_get_bipartition_traversal(testree2->treet_root, comple_bipar_table);
//    mfl_get_bipartition_traversal(testree4->treet_root, comple_bipar_table);
////    mfl_get_bipartition_traversal(testree5->treet_root, comple_bipar_table);
//    tui_print_bipartition_tables(comple_bipar_table);
//    
//    //Destroying the table
////    mfl_destroy_bipartition_table(simple_bipar_table);
////    mfl_destroy_bipartition_table(comple_bipar_table);
//
//    
//    //Destroying the trees
//    mfl_free_tree(simpletree1);
//    mfl_free_tree(simpletree2);
//    mfl_free_tree(testree1);
//    mfl_free_tree(testree2);
//    mfl_free_tree(testree3);
//    mfl_free_tree(testree4);
//    mfl_free_tree(testree5);
//}


void tui_test_consensus_trees(void)
{
    char* testtree1 = NULL;
    char* testtree2 = NULL;
    char* testtree4 = NULL;
    
    testtree1 = (char*)"testtree1 = [&R] ((((1,2),3),((4,5),6)),(((7,8),((9,(10,11)),12)),(13,(14,15))));";
    testtree2 = (char*)"testtree2 = [&R] ((((4,5),6),((1,2),3)),(((10,11),((9,(7,8)),12)),(13,(14,15))));";
    testtree4 = (char*)"testtree4 = [&R] ((((1,2),3),(13,(14,15))),(((7,8),(((4,5),6),12)),(9,(10,11))));";
    
    //Get the trees into mfl_tree_t
    mfl_tree_t* testree1 = mfl_convert_newick_to_mfl_tree_t(testtree1, 0);
    mfl_tree_t* testree2 = mfl_convert_newick_to_mfl_tree_t(testtree2, 0);
    mfl_tree_t* testree4 = mfl_convert_newick_to_mfl_tree_t(testtree4, 0);
    
    //Setting the biparititions
    mfl_set_bipartitions(testree1->treet_root);
    mfl_set_bipartitions(testree2->treet_root);
    mfl_set_bipartitions(testree4->treet_root);
    
    //Initialising the bipartition table
    mfl_bipartition_table* bipar_table = NULL;
    
//    //Filling the bipartition table
//    bipar_table = mfl_initialise_bipartition_table(testree1->treet_num_taxa);
//    mfl_get_bipartition_traversal(testree1->treet_root, bipar_table);
//    mfl_get_bipartition_traversal(testree2->treet_root, bipar_table);
//    mfl_get_bipartition_traversal(testree4->treet_root, bipar_table);
////    tui_print_bipartition_tables(comple_bipar_table);
    
    mfl_tree_t* testing;
    
    testing = mfl_rake_tree(8);
    dbg_printf("Starting tree (rake):\n");
    tui_print_newick(testing->treet_root);
    
    int* clade1 = (int*)mfl_malloc(sizeof(int), 4);
    clade1[0] = 1;
    clade1[1] = 7;
    clade1[2] = 3;
    clade1[3] = 2;
    
//    int* clade1[4] = {1,7,3,2};
    
    mfl_set_tips_in_clade(testing, testing->treet_root, clade1);
    
    dbg_printf("\nOne clade solved:\n");
    tui_print_newick(testing->treet_root);
    
    int* clade2 = (int*)mfl_malloc(sizeof(int), 2);
    clade2[0] = 1;
    clade2[1] = 2;

    
    mfl_set_tips_in_clade(testing, testing->treet_treenodes[19], clade2);
    
    dbg_printf("\nAnother clade solved:\n");
    tui_print_newick(testing->treet_root);
    
    //Destroying the table
    //mfl_destroy_bipartition_table(bipar_table);
    
    //Destroying the trees
    mfl_free_tree(testree1);
    mfl_free_tree(testree2);
    mfl_free_tree(testree4);
}


void tui_test_counts(void)
{
    int i = 0;
    char* intree = (char*)"UNTITLED = [&R] ((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));";
    intree =       (char*)"UNTITLED = [&R] ((((4,((1,2),3)),5),6),(7,(8,(9,((11,12),10)))));";
    intree =       (char*)"UNTITLED = [&R] ((((4,((1,2),3)),5),6),((8,(9,((11,12),10))),7));";
    intree =       (char*)"UNTITLED = [&R] (((8,(9,((11,12),10))),7),(((4,((1,2),3)),5),6));";
    intree =       (char*)"UNTITLED = [&R] (((8,(9,((11,12),10))),7),(((4,(3,(1,2))),5),6));";
    intree =       (char*)"UNTITLED = [&R] (1,(2,(3,(4,(5,(6,(7,(8,(9,(10,(11,12)))))))))));";
//    intree =       (char*)"UNTITLED = [&R] (3,((1,2),(4,(5,(6,(7,(8,(9,(10,(11,12))))))))));";
    
    char* matrices[] = {(char*)"23--1??--032;", // 0
                        (char*)"1---1111---1;", // 1
                        (char*)"1100----1100;", // 2
                        (char*)"11-------100;", // 3
                        (char*)"----1111---1;", // 4
                        (char*)"01----010101;", // 5
                        (char*)"01---1010101;", // 6
                        (char*)"1??--??--100;", // 7
                        (char*)"21--3??--032;", // 8
                        (char*)"11--1??--111;", // 9
                        (char*)"11--1000001-;", // 10
                        (char*)"01------0101;", // 11
                        (char*)"110--?---100;", // 12
                        (char*)"11--1??--111;", // 13
                        (char*)"210--100--21;", // 14
                        (char*)"????----1???;", // 15
                        (char*)"23--1----032;", // 16
                        (char*)"1----1----1-;", // 17
                        (char*)"-1-1-1--1-1-;", // 18
                        (char*)"23--1??--032;", // 19
                        (char*)"--------0101;", // 20
                        (char*)"10101-----01;", // 21
                        (char*)"011--?--0011;", // 22
                        (char*)"110--??--100;", // 23
                        (char*)"11--1000001-;", // 24
                        (char*)"21--1----012;", // 25
                        (char*)"11----111111;", // 26
                        (char*)"10101-----01;", // 27
                        (char*)"210210------;", // 28
                        (char*)"----1111----;", // 29
                        (char*)"230--??1--32;", // 30
                        (char*)"023--??1--32;", // 31
                        (char*)"023-???1--32;", // 32
                        (char*)"23--1?1--023;", // 33
                        (char*)"----1010----;", // 34
                        (char*)"------11---1;", // 35
                        (char*)"10----11---1;", // 36
                        };
    
    int num_matrices = 37;
    int expected[] = {2, 2, 2, 1, 1, 4, 4, 1, 2, 2, 1, 3, 2, 2, 3, 0, 2, 2, 4, 2, 1, 3, 2, 2, 1, 3, 1, 3, 2, 0, 2, 2, 1, 2, 1, 1, 2};
    int testfails = 0;
    
    mfl_handle_s* handle = mfl_t2s(mfl_create_handle());
    
    handle->n_taxa = 12;
    handle->n_chars = 1;
    handle->n_input_newick_trees = 1;
    handle->input_newick_trees = (char**)mfl_malloc(1 * sizeof(char*), 0);
    handle->input_newick_trees[0] = intree;
    
    for (i = 0; i < num_matrices; ++i) {
        if (i == 33) {
            dbg_printf("break\n");
        }
        handle->input_data = matrices[i];
        if (tui_test_treelength_calculation(handle, &expected[i])) {
            ++testfails;
        }
    }
    
    dbg_printf("\n");
    dbg_printf("SUMMARY:\n");
    dbg_printf("Test of all topologies for all matrices concluded\n");
    if (testfails) {
        dbg_printf("Test FAILED %i times.\n", testfails);
    }
    else {
        dbg_printf("All tests PASSED\n");
    }
    dbg_printf("\n");
}


int main (int argc, char *argv[])
{
    dbg_printf("\n\t****************************************\n\n");
    dbg_printf("\t  Welcome to the Morphy Test Interface\n\n");
    dbg_printf("\t****************************************\n\n\n");
    
    dbg_printf("Missing data ");
    if (MORPHY_MISSING_DATA_BITWISE & 1) {
        dbg_printf("includes lowest order bit (1...1111).\n");
    }
    else {
        dbg_printf("DOES NOT include lowest order bit (1...1110).\n");
    }
    dbg_printf("\n");
    
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
    
//    dbg_printf("Testing the n-ary ring creation:\n");
//    tui_test_nary_ring_creation();
//    dbg_printf("\nEnd n-ary ring test\n");
//
    dbg_printf("Testing character optimisation (yay!):\n");
    tui_test_basic_character_optimisation();
    dbg_printf("\nEnd character optimisation test\n");
    
    dbg_printf("Testing local reoptimisation functions:\n");
    tui_test_local_reoptimisation();
    dbg_printf("\nEnd local optimisation test\n");
    
//    dbg_printf("Testing addition sequence:\n");
//    tui_test_addition_sequence();
//    dbg_printf("\nEnd addition sequence test\n");
    
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

    dbg_printf("Testing matrices:\n");
    tui_test_counts();
    dbg_printf("\nEnd counts test\n");
    
    dbg_printf("Testing compare and replace:\n");
    tui_test_compare_replace();
    dbg_printf("\nEnd compare and replace test\n");
    
    dbg_printf("Testing stepwise addition:\n");
    tui_test_stepwise_addition();
    dbg_printf("\nEnd stepwise addition test\n");
    
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
