//
//  mfl_drawtree.c
//  Morphy
//
//  Created by mbrazeau on 17/04/2016.
//
//
#include "morphy.h"
#include "tuimfy.h"
#include "mfl_drawtree.h"

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


void mfl_put_character_in_cell(char const ch, int row, int col, char* grid)
{
    grid[ row * (MAX_BUFFER_WIDTH + 1) + col ] = ch;
}


char mfl_drawtree_get_character_in_cell(int row, int col, char* grid)
{
    return grid[ row * (MAX_BUFFER_WIDTH + 1) + col ];
}


void mfl_drawtree_write_into_tipfield(char* name, char* grid, int row, int col)
{
    int i = 0;
    
    while ((mfl_drawtree_get_character_in_cell(row, col+i+1, grid) != '\n') && name[i] != '\n' && name[i]) {
        mfl_put_character_in_cell(name[i], row, col+i+1, grid);
        ++i;
    }
}


void mfl_drawtree_set_coords_traversal(mfl_node_t *n, int row, int num_taxa)
{
    int branchlen = (2*63) / num_taxa; // TODO: Get rid of this magic number; calculate estimated tree depth
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


void mfl_drawtree_set_coords_traversal(mfl_node_t* n, int* currentrow, char* grid, int num_taxa)
{
    mfl_node_t* p;
    
    if (n->nodet_tip) {
        mfl_drawtree_set_coords_traversal(n, *currentrow, num_taxa);
        *currentrow = *currentrow + 2;
        return;
    }
    
    p = n->nodet_next;
    do {
        mfl_drawtree_set_coords_traversal(p->nodet_edge, currentrow, grid, num_taxa);

        p = p->nodet_next;
    } while (p != n);
    
    
    mfl_drawtree_set_coords_traversal(n, *currentrow, num_taxa);
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
            mfl_drawtree_write_into_tipfield(n->nodet_tipname, grid, n->row, n->col);
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
    
    mfl_put_character_in_cell('+', n->row, n->col, grid);
    mfl_drawtree_add_nodebar(n, ldesc, rdesc, grid);
    return;
}


char *mfl_drawtree(mfl_tree_t* t)
{
    char *treedrawing = mfl_drawtree_create_virtual_grid(t->treet_num_taxa);
    int firstrow = 0;
    bool wasrooted = false;
    
    if (!t->treet_root) {
        // TODO: Root the tree.
        wasrooted = true;
    }
    
    mfl_drawtree_set_coords_traversal(t->treet_root, &firstrow, treedrawing, t->treet_num_taxa);
    mfl_drawtree_draw_traversal(t->treet_root, treedrawing);
    
    if (wasrooted) {
        // TODO: Unroot the tree.
    }
    
    return treedrawing;
}