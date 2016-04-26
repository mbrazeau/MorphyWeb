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


/*!
 @discussion Allocates memory for an array of characters with newline characters
 embedded in each row at the 81st character. This forms a coordinate space where
 the tree can be drawn.
 @param num_taxa (int) the number of taxa in the tree to be drawn
 @return pointer to character array with spaces and newlines embedded.
 */
char* mfl_drawtree_create_virtual_grid(int num_taxa)
{
    int i;
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


/*!
 @discussion Puts the specified character into the virtual grid at the coordinates
 specified so that you don't need to think about the math every time.
 @param ch (const char) the character to place in the matrix
 @param row (int) the destination row for ch
 @param col (int) the destination column for ch
 @param grid (char*) the destination virtual grid for ch
 */
void mfl_put_character_in_cell(char const ch, int row, int col, char* grid)
{
    grid[ row * (MAX_BUFFER_WIDTH + 1) + col ] = ch;
}


/*!
 @discussion Returns the character found in the virtual grid at a specified
 set of coordinates
 @param row (int) the source row
 @param col (int) the source column
 @param grid (char*) the source virtual grid being queried
 @return the character found at the input coordinates
 */
char mfl_drawtree_get_character_in_cell(int row, int col, char* grid)
{
    return grid[ row * (MAX_BUFFER_WIDTH + 1) + col ];
}


/*!
 @discussion Writes a given character string into the space in the grid reserved
 for tip information (i.e. tip name or tip number or both). Without calling this 
 function, the tree will have no tip labels.
 @param name (char*) the string corresponding to the tip's name or number
 @param grid (char*) the virtual grid where the tree is being drawn
 @param row (int) the tip's row
 @param col (the) the tip's colum
 */
void mfl_drawtree_write_into_tipfield(char* name, char* grid, int row, int col)
{
    int i = 0;
    
    while ((mfl_drawtree_get_character_in_cell(row, col+i+1, grid) != '\n') && name[i] != '\n' && name[i]) {
        mfl_put_character_in_cell(name[i], row, col+i+1, grid);
        ++i;
    }
}

int mfl_drawtree_calculate_branchlength(int num_taxa)
{
    // TODO: There will be a more sophisticated way of doing this.
    
    int num_units = 0;
    float branchlen;
    
    branchlen = 3 * 63 / (num_taxa - 1);
    
    if (branchlen < 1) {
        num_units = 1;
    }
    else {
        num_units = (int)branchlen;
        if (num_units < 1) {
            num_units = 1;
        }
    }
    
    return num_units;
}

/*!
 @discussion Sets the spatial coordinates for ASCII text drawing for a tip or 
 internal node based on either the input row or the row and column levels of 
 the descendant branches.
 @param n (mfl_node_t*) target node
 @param row (int) the row number for drawing a tip
 @param num_taxa (int) the number of terminal taxa
 */
void mfl_drawtree_set_node_coords(mfl_node_t *n, int row, int num_taxa)
{
    int branchlen = mfl_drawtree_calculate_branchlength(num_taxa);
    
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


/*!
 @discussion Applies the branch subtending a descendant of the input parent node
 into the gridspace for the drawing.
 @param parent (mfl_node_t*) the parent node
 @param desc (mfl_node_t*) the descendant of the parent
 @param grid (char*) the character 'gridspace' for the drawing.
 */
void mfl_drawtree_apply_subbranch(mfl_node_t* parent, mfl_node_t* desc, char *grid)
{
    int i = 0;
    
    for (i = 1; i < desc->col - parent->col; ++i) {
        mfl_put_character_in_cell('-', desc->row, parent->col + i, grid);
    }
}


/*!
 @discussion Traverses the tree in postorder setting the coordinates in the 
 gridspace for each node
 @param n (mfl_node_t*) the node to be visited on the traversal
 @param currentrow (int*) pointer to the value for the current row for tips. 
 Because the traversal is left-handed, the initial value for the int at this
 pointer will be set to zero and incremented as tips are added to the tree 
 drawing
 @param grid (char*) the virtual gridspace into which the ASCII text drawing 
 is applied.
 @param num_taxa (int) the number of active terminal taxa in the input tree.
 */
void mfl_drawtree_set_coords_traversal(mfl_node_t* n, int* currentrow, char* grid, int num_taxa)
{
    mfl_node_t* p;
    
    if (n->nodet_tip) {
        mfl_drawtree_set_node_coords(n, *currentrow, num_taxa);
        *currentrow = *currentrow + 2;
        return;
    }
    
    p = n->nodet_next;
    do {
        mfl_drawtree_set_coords_traversal(p->nodet_edge, currentrow, grid, num_taxa);

        p = p->nodet_next;
    } while (p != n);
    
    
    mfl_drawtree_set_node_coords(n, *currentrow, num_taxa);
}


/*!
 @discussion Called on internal nodes only, this
 @note this function could be updated to find ldesc and rdesc itself, and thus
 be a bit safer.
 @param n (mfl_node_t*) the parent node
 @param ldesc (mfl_node_t*) the left-most descendant of n
 @param rdesc (mfl_node_t*) the right-most descendant of n
 @param grid (char*) the character gridspace for the tree drawing.
 */
void mfl_drawtree_add_nodebar(mfl_node_t* n, mfl_node_t* ldesc, mfl_node_t* rdesc, char* grid)
{
    // TODO: Update this function to do its own job of finding ldesc and rdesc.
    int i = 0;
    
    if (n->nodet_tip) {
        dbg_eprintf("function cannot be called on a terminal node.");
        return;
    }
    
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


/*!
 @discussion Traverses the tree inputting characters into the gridspace based
 on the coordinates created in mfl_drawtree_set_coords_traversal.
 @param n (mfl_node_t*) the node being visited
 @param grid (char*) the character gridspace for the drawing.
 */
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


/*!
 @discussion Creates a drawing of an input tree as an 80-character wide ASCII 
 text representation.
 @param t (mfl_tree_t*) any Morphy-type input tree
 @return a two-dimensional C-type string drawing of the tree.
 */
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