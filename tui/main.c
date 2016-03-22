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
    
    dbg_printf("Free the tree\n\n");
    
    mfl_free_tree(newtree, num_taxa, num_nodes);
    
    dbg_printf("Exit the program\n\n");
    
	return 0;
}
