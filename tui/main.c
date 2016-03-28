#include "morphy.h"

void mfl_test_newick_stuff();

int main (void)
{
    dbg_printf("\n\t****************************************\n\n");
    dbg_printf("\t  Welcome to the Morphy Test Interface\n\n");
    dbg_printf("\t****************************************\n\n\n");
    
    dbg_printf("Generate a new tree\n\n");
    
    mfl_test_newick_stuff();
    
    mfl_tree_t *newtree;
    int num_taxa = 5;
    int num_nodes = 0;
    
    newtree = mfl_alloctree_with_nodes(num_taxa);
    
    /* Setting up a binary fork to start testing some of the new functions. Currently hard-coded*/
    mfl_create_binary_fork(newtree->treet_treenodes[num_taxa], newtree->treet_treenodes[0], newtree->treet_treenodes[1]);
    
    /* Adding a hard-coded branch to the tree */
    mfl_make_new_n_ary_ring_node(newtree->treet_treenodes[num_taxa+1], 2);
    mfl_join_node_edges(newtree->treet_treenodes[2], newtree->treet_treenodes[num_taxa + 1]->nodet_next);
    mfl_insert_branch(newtree->treet_treenodes[num_taxa + 1], newtree->treet_treenodes[num_taxa + 1]->nodet_next->nodet_next, newtree->treet_treenodes[1]);
    
    dbg_printf("Free the tree\n\n");
    
    mfl_free_tree(newtree, num_taxa, num_nodes);
    
    dbg_printf("Exit the program\n\n");
    
	return 0;
}
