//
//  mfl_drawtree.c
//  Morphy
//
//  Created by mbrazeau on 17/04/2016.
//
//
#include "morphy.h"
#include "tuimfy.h"
#define MAX_BUFFER_WIDTH 80
#define DEFAULT_TIP_COLUMN 63

char* mfl_drawtree_create_virtual_grid(int num_taxa)
{
    int i, j;
    //int max_buffer_width = 80;
    int depth = 0;
    char *virtual_grid = NULL;
    
    depth = 2 * num_taxa;
    
    virtual_grid = (char*)malloc(((MAX_BUFFER_WIDTH+1) * depth + 1) * sizeof(char)); // +1 to buffer for a newline character at buffer end of each row.
                                                                                     // +1 to sum for the terminal null character. 
    if (!virtual_grid) {
        dbg_eprintf("unable to allocate memory for tree drawing. Returning NULL");
    }
    else {
        memset(virtual_grid, ' ', ((MAX_BUFFER_WIDTH+1) * depth + 1) * sizeof(char));
    }
    
    for (i = 0; i < depth; ++i) {
        virtual_grid[(MAX_BUFFER_WIDTH * (i + 1)) + i] = '\n';
    }
    virtual_grid[(MAX_BUFFER_WIDTH+1) * depth] = '\0';
    
    return virtual_grid;
    
}

void tui_move_print_head_down(char **printhead, const int max_buffer_width)
{
    
    *printhead = *printhead + 2* (max_buffer_width + 1);// * sizeof(char*);
}


void mfl_put_character_in_cell(char const ch, int row, int col, char* grid)
{
    grid[ row * (MAX_BUFFER_WIDTH + 1) + col ] = ch;
}

char mfl_get_character_in_cell(int row, int col, char* grid)
{
    return grid[ row * (MAX_BUFFER_WIDTH + 1) + col ];
}

void mfl_write_into_tipfield(char* name, char* grid, int row, int col)
{
    int i = 0;
    
    while ((mfl_get_character_in_cell(row, col+i+1, grid) != '\n') && name[i] != '\n' && name[i]) {
        mfl_put_character_in_cell(name[i], row, col+i+1, grid);
        ++i;
    }
}




void mfl_printree_set_node_coordinates(mfl_node_t *n, int row)
{
    int branchlen = 8; // TODO: Get rid of this magic number; calculate estimated tree depth
    mfl_node_t* p = NULL;
    mfl_node_t* leftdesc = NULL;
    mfl_node_t* rightdesc = NULL;
    
    if (n->nodet_tip) {
        n->row = row;
        n->col = DEFAULT_TIP_COLUMN;
    }
    else {
        
        p = n->nodet_next;
        leftdesc = p->nodet_edge;
        
        do {
            rightdesc = p->nodet_edge;
            p = p->nodet_next;
        } while (p != n);
        
        // Set the row level:
        // This next part should get replaced by something smarter that knows about how many descendant tips there are
        // and what the largest and smallest row numbers are.
        
        if (rightdesc->row > leftdesc->row) {
            n->row = leftdesc->row + ((rightdesc->row - leftdesc->row) / 2);
        }
        else {
            n->row = rightdesc->row + ((leftdesc->row - rightdesc->row) / 2);
        }
        
        // Now set the new column level. Replace as above.
        if (leftdesc->col < rightdesc->col) {
            n->col = leftdesc->col - branchlen;
        }
        else {
            n->col = rightdesc->col - branchlen;
        }
        
    }
}

void mfl_drawtree_apply_subbranch(mfl_node_t* parent, mfl_node_t* desc, char *grid)
{
    int i = 0;
    int brlen = desc->branchl_cdraw;
    
    for (i = 1; i < desc->col - parent->col; ++i) {
        mfl_put_character_in_cell('-', desc->row, parent->col + i, grid);
    }
}

void mfl_printtree_set_coords_traversal(mfl_node_t* n, int* currentrow, char* grid)
{
    mfl_node_t* p;
    
    if (n->nodet_tip) {
        mfl_printree_set_node_coordinates(n, *currentrow);
        *currentrow = *currentrow + 2;
        return;
    }
    
    p = n->nodet_next;
    do {
        mfl_printtree_set_coords_traversal(p->nodet_edge, currentrow, grid);

        p = p->nodet_next;
    } while (p != n);
    
    
    mfl_printree_set_node_coordinates(n, *currentrow);
}

void mfl_drawtree_add_nodebar(mfl_node_t* n, mfl_node_t* ldesc, mfl_node_t* rdesc, char* grid)
{
    int i = 0;
    
    // Right descendant
    do {
        ++i;
        mfl_put_character_in_cell('|', n->row + i, n->col, grid);
    } while (n->row + i < rdesc->row);
    
    mfl_put_character_in_cell('\\', n->row + i, n->col, grid);
    
    // Left descendant
    i = 0;
    do {
        ++i;
        mfl_put_character_in_cell('|', n->row - i, n->col, grid);
    } while (n->row - i > ldesc->row);
    
    mfl_put_character_in_cell('/', n->row - i, n->col, grid);
    
}

void mfl_drawtree_draw_traversal(mfl_node_t* n, char *grid)
{
    mfl_node_t* p = NULL;
    mfl_node_t* ldesc = NULL;
    mfl_node_t* rdesc = NULL;
    
    if (n->nodet_tip) {
        if (n->nodet_tipname) {
            mfl_write_into_tipfield(n->nodet_tipname, grid, n->row, n->col);
        }
        return;
    }
    
    p = n->nodet_next;
    ldesc = p->nodet_edge;
    
    do {
        mfl_drawtree_draw_traversal(p->nodet_edge, grid);
        mfl_drawtree_apply_subbranch(n, p->nodet_edge, grid);
        rdesc = p->nodet_edge;
        p = p->nodet_next;
    } while (p != n);
    
    dbg_printf("Vising internal node: %i\n", n->nodet_index);
    
    mfl_put_character_in_cell('+', n->row, n->col, grid);
    mfl_drawtree_add_nodebar(n, ldesc, rdesc, grid);
    return;
}

void tui_test_tree_printing()
{
    
    char *grid = mfl_drawtree_create_virtual_grid(7);
    
    char newick[] = "[&R] = ((2,1,3),4,(7,(5,6)));";
    
    mfl_tree_t* printme = mfl_convert_newick_to_mfl_tree_t(newick, 7);
    
    printme->treet_treenodes[0]->nodet_tipname = (char*)"Borosaurus jibberstoni";
    printme->treet_treenodes[1]->nodet_tipname = (char*)"Wayne";
    printme->treet_treenodes[2]->nodet_tipname = (char*)"Ur mom";
    printme->treet_treenodes[3]->nodet_tipname = (char*)"Taxon 1";
    printme->treet_treenodes[4]->nodet_tipname = (char*)"Taxon Two";
    printme->treet_treenodes[5]->nodet_tipname = (char*)"Foo";
    printme->treet_treenodes[6]->nodet_tipname = (char*)"Bar";
    printme->treet_treenodes[7]->nodet_tipname = (char*)"Fubarus";
    
    int firstrow = 0;
    mfl_printtree_set_coords_traversal(printme->treet_root, &firstrow, grid);
    mfl_drawtree_draw_traversal(printme->treet_root, grid);
    //dbg_printf("________\n");
    dbg_printf("          1.........2.........3.........4.........5.........6.........7.........8\n");
    dbg_printf("012345678901234567890123456789012345678901234567890123456789012345678901234567890\n");
    dbg_printf("%s", grid);
    
    return;
}