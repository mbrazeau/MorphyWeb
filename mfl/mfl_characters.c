/*
 *  mfl_characters.c
 *
 *
 *  THE MORPHY FUNCTION LIBRARY
 *  A library for phylogenetic analysis with emphasis on parsimony and
 *  morphology (but someday other methods)
 *
 *  Copyright (C) 2016  by Martin D. Brazeau, Thomas Guillerme,
 *  and Chris Desjardins
 *
 *  Created by Martin Brazeau and Chris Desjardins on 1/26/12.
 *  Some data structs, routines and ideas derived from:
 *      - The PHYLIP package by Joe Felsenstein
 *          <http://evolution.genetics.washington.edu/phylip.html>
 *      - MrBayes by John Huelsenbeck and Fredrik Ronquist
 *          <http://mrbayes.sourceforge.net/>
 *
 *  Any bugs, errors, inefficiences and general amateurish handling are our own
 *  and most likely the responsibility of MDB. We make no guarantees.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


/* In this file:
 * Functions for processing the character input. These functions assign
 * the character types to their relevant partitions as mfl_charstate type
 * vectors in an mfl_matrix_t.*/

#include "morphy.h"

/*!
 ## bool mfl_is_valid_morphy_ctype(char c)
 Checks that a symbol is a valid type used by Morphy for conversion
 to bit code as an mfl_charstate type. This is a critical function
 in Morphy upon which many other functions might depend. If any
 new symbols are to be added/recognized, they should be added here.
 However, with a maximum of 64-1 states prossible, 65 different
 symbols (+'?') should be sufficient. Valid symbols include:
 
 - Alphanumeric characters
 - '?'
 - '-'
 - '+'
 - at symbol
 
 @param c (int) the input character
 @returns (bool) 0 for false, 1 for true
 */
bool mfl_is_valid_morphy_ctype(char c)
{
    if (isalnum(c)) {
        return true;
    }
    else if (c == '?') {
        return true;
    }
    else if (c == '-') {
        return true;
    }
    else if (c == '+') {
        return true;
    }
    else if (c == '@') {
        return true;
    }
    else {
        return false;
    }
}

/*!
 ## mfl_nodedata_t* mfl_alloc_datapart(void)
 Allocates a data partition.

 @return pointer to a datapartition (mfl_nodedata_t)
 */
mfl_nodedata_t* mfl_alloc_datapart(void)
{
    mfl_nodedata_t* newdatapart = NULL;
    
    newdatapart = (mfl_nodedata_t*)malloc(sizeof(mfl_nodedata_t));
    
    if (!newdatapart) {
        dbg_printf("ERROR in mfl_alloc_datapart(): unable to allocate memory for new data partition.\n");
        return NULL;
    }
    else {
        memset(newdatapart, 0, sizeof(mfl_nodedata_t));
    }
    
    return newdatapart;
}

/*!
 ## void mfl_free_nodedata(mfl_nodedata_t *olddata)
 Frees all of the data used in a SINGLE data partition from within a
 node. Any other associated memory allocated with a datapartition should
 be freed in here before returning.
 @param olddata (mfl_nodedata_t*) corresponding to the partition being freed
 */
void mfl_free_nodedata(mfl_nodedata_t *olddata)
{
    
    if (olddata->nd_prelim_set) {
        free(olddata->nd_prelim_set);
        olddata->nd_prelim_set = NULL;
    }
    if (olddata->nd_final_set) {
        free(olddata->nd_final_set);
        olddata->nd_final_set = NULL;
    }
    if (olddata->nd_subtree_prelim_set) {
        free(olddata->nd_subtree_prelim_set);
        olddata->nd_subtree_prelim_set = NULL;
    }
    if (olddata->nd_subtree_final_set) {
        free(olddata->nd_subtree_final_set);
        olddata->nd_subtree_final_set = NULL;
    }
    
    // Depending on how cost matrices are handled, they might be freed here, too.
    
    free(olddata);
}

/*!
 # mfl_charstate* mfl_allocate_nodal_character_set(int num_characters)
 Allocates a nodal character set used as either the downpass or uppass set
 within a particular partition for a particular node.
 @param num_characters (int) num_characters the number of characters
 @return mfl_charstate* the vector of characters.
 */
mfl_charstate* mfl_allocate_nodal_character_set(int num_characters)
{
    mfl_charstate *newcharset = NULL;
    
    newcharset = (mfl_charstate*)malloc(num_characters * sizeof(mfl_charstate));
    if (!newcharset) {
        dbg_printf("ERROR in mfl_allocate_nodal_character_set(): unable to allocate memory for character data.\n");
        return NULL;
    }
    else {
        memset(newcharset, 0, num_characters * sizeof(mfl_charstate));
    }
    
    return newcharset;
}

/*!
 #int mfl_convert_nexus_symbol_to_shift_value(char in, char *datype_converter)
 Looks up the symbol held by 'in' in the datype_converter array. The index
 within the array where the symbols match is found, plus the number of special 
 states reserved in Morphy determines the number of shifts to be used to convert 
 the symbol. Each character partition could have its own datype_converter and its 
 own way of shuffling the symbols. Also, if the Nexus standard is followed and 
 sybols list is supplied, then the order in which input symbols are declared 
 determines their ordering in asymmetric characters.
 @param in (char) the input symbol
 @param datype_converter (char*) the array of datatypes in their
 */
int mfl_convert_nexus_symbol_to_shift_value(char in, char *datype_converter)
{
    int i = 0;
    
    do {
        if (in == datype_converter[i]) {
            break;
        }
        ++i;
    }
    while (datype_converter[i]);
    
    return i + MORPHY_SPECIAL_STATE_PAD;
}


/*!
 # int mfl_shift_value_DEFAULT_catdata(char in)
 Takes an alphanumeric value on the scale 0-9ABC...xyz and returns an int
 representing the number of shifts needed for setting its corresponding
 bit in an mfl_charstate. This function is used as a DEFAULT if there is no 
 pre-specified list of state symbols. In effect, Morphy will convert them to 
 a character state corresponding to their order as a default.
 @param in (char) input character
 @returns the shift value of the alphanumeric character
 */
int mfl_shift_value_DEFAULT_catdata(char in)
{
    if (isdigit(in)) {
        return in - '0' + MORPHY_SPECIAL_STATE_PAD;
    }
    else if (isalpha(in))
    {
        if (isupper(in)) {
            return in - 'A' + 11 + MORPHY_SPECIAL_STATE_PAD;
        }
        else {
            return in - 'a' + 11 + 26 + MORPHY_SPECIAL_STATE_PAD;
        }
    }
    else if (in == '@') {
        return in - 'a' + 11 + 26 + MORPHY_SPECIAL_STATE_PAD + 1;
    }
    else if (in == '+') {
        return in - 'a' + 11 + 26 + MORPHY_SPECIAL_STATE_PAD + 2;
    }
    else {
        dbg_printf("ERROR in mfl_convert_alphanum_to_shift_value(): input is not alphanumeric\n");
        return 0;
    }
}


/**
 # bool mfl_char_is_DNA_base(char in)
 Checks whether an input character is a valid single-base DNA symbol
 @param in (char) in the input character
 @return bool true/false
 */
bool mfl_char_is_DNA_base(char in)
{
    if (in == 'a' || in == 'A') {
        return true;
    }
    else if (in == 'c' || in == 'C') {
        return true;
    }
    else if (in == 'g' || in == 'G') {
        return true;
    }
    else if (in == 't' || in == 'T') {
        return true;
    }
    else {
        return false;
    }
}

/**
 */
int mfl_shift_value_DNA_base(char in)
{
    if (in == 'a' || in == 'A') {
        return 2;
    }
    else if (in == 'c' || in == 'C') {
        return 3;
    }
    else if (in == 'g' || in == 'G') {
        return 4;
    }
    else if (in == 't' || in == 'T') {
        return 5;
    }
    else {
        dbg_eprintf("input not valid DNA base");
        return 0;
    }
}


/**
 */
mfl_charstate mfl_set_DNA_bits(char in)
{
    // Check is single-state symbol
    //  Get conversion number
    // If IUPAC symbol, then do the bitwise operations needed
}


int mfl_shift_value_amino_acid(char state_symbol)
{
    
}


/*!
 #mfl_charstate mfl_convert_nexus_multistate(char *xstates, char *datype_converter)
 Converts a string of character states corresponding to a matrix cell containing
 multistate values. For each state symbole, the shift value of the data type converter 
 is calculated and that bit 'shifted in' by that amount and then added to the output
 with a bitwise OR.
 @param xstates (char*) the string corresponding to the multistate cell.
 @param datype_converter (char*) the ordered array of symbols used to calculated the number of shifts
 @returns a mfl_charstate value
 */
mfl_charstate mfl_convert_nexus_multistate(char *xstates, char *datype_converter)
{
    u_int64_t inbit = 0;
    int shift_val;
    mfl_charstate xstate_bits = 0;
    
    do {
        inbit = 1;
        shift_val = mfl_convert_nexus_symbol_to_shift_value(*xstates, datype_converter);
        inbit = inbit << shift_val;
        xstate_bits = xstate_bits | inbit;
        ++xstates;
    } while (xstates && *xstates != ')' && *xstates != '}');
    
    return xstate_bits;
}


/**
 # mfl_charstate mfl_convert_gap_character(mfl_optimisation_t opt_method)
 Sets the bits in an mfl_charstate variable to the appropriate value and returns
 the result.
 @param gapmethod (mfl_optimisation_t) the optimality method applied to the gap character.
 @return mfl_charstate bits set to the corresponding treatment.
 */
mfl_charstate mfl_convert_gap_character(mfl_gap_t gapmethod)
{
    if (gapmethod == MFL_GAP_INAPPLICABLE || gapmethod == MFL_GAP_NEWSTATE) {
        return MORPHY_INAPPLICABLE_BITPOS;
    }
    else {
        return MORPHY_MISSING_DATA_BITWISE;
    }
}


/*!
 If the typeset list is supplied as an input string rather than an array, this
 sets the new char array to be a space-less vector of symbols use as character
 states in a Nexus-derived "Format Symbols" command. A lookup into this will return
 the required shift value.
 @param datype_converter (char*) the set of symbols to convert presented as a string
 @param symbols_list (char*) the subtoken of symbols taken from input
 */
void mfl_set_datatype_converter_from_symbol_list(char* datype_converter, char* symbols_list)
{
    
    int i = 0;
    char *dtype_ptr = NULL;
    
    dtype_ptr = symbols_list;
    
    while (dtype_ptr) {
        if (!isspace(*dtype_ptr)) {
            datype_converter[i] = *dtype_ptr;
            ++i;
        }
        ++dtype_ptr;
    }
}


char* mfl_search_char_in_chartypes_array(char* key, char* list, int *listend, int listmax)
{
    int i = 0;
    char *ptr = NULL;
    
    for (i = 0; i < *listend; ++i) {
        if (*key == *(list + i)) {
            ptr = (list + i);
            break;
        }
    }
    
    if (!ptr) {
        if (*listend < listmax) {
            list[*listend] = *key;
            ++(*listend);
            list[*listend] = '\0';
        }
        else {
            dbg_printf("ERROR in mfl_search_in_chartypes_array(): %s: line: %i\n", __FILE__, __LINE__);
            dbg_printf("\t Insufficient space in array for %i member \'%c\'. Exceeding max states: %i.\n", *listend, *key, MORPHY_MAX_STATE_NUMBER);
            dbg_printf("Returning NULL\n\n");
            ++(*listend); // Exceed list max for error reporting.
            return NULL;
        }
    }
    
    return ptr;
}


int mfl_get_numstates_from_matrix(char *inputmatrix)
{
    /* When no state symbols are specified and default reading is in effect,
     * the number of unique symbols used can, with some assumptions, be
     * counted directly from the matrix.
     * */
    
    int count = 0;
    char *current = NULL;
    int listmax = MORPHY_MAX_STATE_NUMBER+MORPHY_SPECIAL_STATE_PAD; // +1 for terminal null.
    char statesymbols[listmax];
    int dbg_loopcount = 0;
    
    statesymbols[0] = NULL;
    current = inputmatrix;
    
    do {
        if (*current != '?' && mfl_is_valid_morphy_ctype(*current)) {
            mfl_search_char_in_chartypes_array(current, statesymbols, &count, listmax);
        }
        ++current;
        ++dbg_loopcount;
    } while (*current);
    
    
    //CHECK FOR ERROR HERE
    if (count > listmax) {
        dbg_printf("ERROR in %s, %s, line: %i\n", __FUNCTION__, __FILE__, __LINE__);
        dbg_printf("\t State symbols outnumber MORPHY_MAX_STATE_NUMBER. Returning NULL\n\n");
        return NULL;
    }
    
    
    return count-1;
}


mfl_parsimony_t* mfl_alloc_chartype_list(int input_num_chars)
{
    mfl_parsimony_t* chars_by_ctype = (mfl_parsimony_t*)malloc(input_num_chars * sizeof(mfl_parsimony_t));
    
    if (!chars_by_ctype) {
        dbg_eprintf("unable to allocate memory for new chars_by_ctype array\n");
    }
    else {
        memset(chars_by_ctype, MORPHY_DEFAULT_PARSIMONY_METHOD, input_num_chars * sizeof(mfl_parsimony_t));
    }
    
    return chars_by_ctype;
}


void mfl_destroy_chartype_list(mfl_parsimony_t *ctype_list)
{
    free(ctype_list);
}

mfl_parsimony_t* mfl_get_chartypes_list(mfl_handle_s* mfl_handle)
{
    int i;
    mfl_parsimony_t* ctype_setters = mfl_alloc_chartype_list(mfl_handle->n_chars);
    
    if (mfl_handle->ctypes) {
        memccpy(ctype_setters, mfl_handle->ctypes, mfl_handle->n_chars, sizeof(mfl_parsimony_t));
    }
    else if (mfl_handle->ctypes_cmd) {
        if (!mfl_handle->n_ctypes) {
            dbg_printf("All characters are of type %s\n", MORPHY_DEFAULT_PARSIMONY_NAME);
        }
        else {
            // Parse each ctype command
            for (i = 0; i < MFL_OPT_MAX; ++i) {
                // THE IMPORTANT BITS HERE
                if (mfl_handle->ctypes_cmd[i]) {
                    mfl_read_nexus_style_list_subcmd(mfl_handle->ctypes_cmd[i], i, (int*)ctype_setters, mfl_handle->n_chars);
                }
            }
        }
    }
    else {
        dbg_printf("All characters are of type %s\n", MORPHY_DEFAULT_PARSIMONY_NAME);
    }
    
    return ctype_setters;
}


mfl_charstate mfl_standard_conversion_rules(char *c, mfl_gap_t gaprule)
{
    
}


/*mfl_charstate mfl_CUSTOM_CONVERSION_RULES(char *c, mfl_gap_t gaprule)
{
    __TEMPLATE__FOR__CUSTOM__CONVERSION__RULE
}*/


/*!
 An interface for the conversion rule used to turn input characters into
 set bits in an mfl_charstate variable. If ever there is a need to write 
 different conversion rules, they can be set in here.
 @param parsimtyme (mfl_parsimony_t) the type of parsimony function used
 @return Pointer to the wrapper on the conversion rules to be applied to
 the data.
 */
mfl_char2bit_fn mfl_fetch_conversion_rule(mfl_parsimony_t parsimtype)
{
    mfl_char2bit_fn ret = NULL;
    
    switch (parsimtype) {
        
        case MFL_OPT_FITCH:
            ret = mfl_standard_conversion_rules;
            break;
        case MFL_OPT_WAGNER:
            ret = mfl_standard_conversion_rules;
            break;
        case MFL_OPT_DOLLO_UP:
            ret = mfl_standard_conversion_rules;
            break;
        case MFL_OPT_DOLLO_DN:
            ret = mfl_standard_conversion_rules;
            break;
        case MFL_OPT_IRREVERSIBLE_UP:
            ret = mfl_standard_conversion_rules;
            break;
        case MFL_OPT_IRREVERSIBLE_DN:
            ret = mfl_standard_conversion_rules;
            break;
        case MFL_OPT_COST_MATRIX:
            ret = mfl_standard_conversion_rules;
            break;
        /*case ___CUSTOM_TYPE___:
             ret = __POINTER_TO_CUSTOM_CONVERSION_ROUTINE;
             break;*/
        default:
            break;
    }
    
    return ret;
}

void mfl_apply_conversion_rule(mfl_character_vector_t *character)
{
    character->cv_conversion_rule = mfl_fetch_conversion_rule(character->cv_parsim_method);
}

/*
 Converting the columns
 */

// Has a gap? What's the gap conversion rule? Convert for that.

// Did a datatype command come in?
//  If Yes:
//      Is it a mixed set (e.g. Std. and DNA)?
//      Did a symbols list come in?
//          If yes:
//              Set converter according to symbols list
//          If no:
//              Set converter according to defaults
//
//  If no:
//      Set a default or try to detect type?
//          For now: use default.
//
//

void mfl_setup_convert_charcells_to_mfl_charstates(mfl_matrix_t* m, mfl_handle_t handle)
{
    
}

/*
 *
 * Functions for parsing the input matrix string into substrings. 
 *
 */

void mfl_copy_singleton_subtoken_to_substring(char singleton, char* substring)
{
    *substring = singleton;
    substring[1] = '\0';
}


void mfl_copy_multistate_subtoken_to_substring(char** xstatetoken, char* substring)
{
    int i = 0;
    if (**xstatetoken == '{' || **xstatetoken == '(') {
        ++(*xstatetoken);
    }
    
    do {
        if (isalnum(**xstatetoken)) {   // Morphy ignores '?' and '-' in multistate taxa
            substring[i] = **xstatetoken;
            ++(*xstatetoken);
            ++i;
        }
        else if (**xstatetoken == '?' || **xstatetoken == '-') {
            ++(*xstatetoken);
        }
    } while (**xstatetoken != ')' && **xstatetoken != '}');
    substring[i] = '\0';
}


void mfl_populate_chartype_character_vector(mfl_matrix_t *matrix, char *input_data_matrix, int num_chars, int num_taxa)
{
    int column = 0;
    int row = 0;
    char *current = NULL;
    char *substring = NULL;
    
    current = input_data_matrix;
    
    substring = (char*)malloc((matrix->max_states + 1) * sizeof(char));
    if (!substring) {
        dbg_printf("ERROR in mfl_populate_chartype_character_vector(): unable to allocate memory for character substring\n");
        return;
    }
    else {
        memset(substring, 0, (matrix->max_states + 1) * sizeof(char));
    }
    
    do {
        // Collect the input matrix substring
        column = 0;
        do {
            if (mfl_is_valid_morphy_ctype(*current)) {
                mfl_copy_singleton_subtoken_to_substring(*current, substring);
                strcpy(matrix->mat_matrix[column]->cv_character_cells[row], substring);
                
                if (*current == '-') {                                  // This does obfuscate a process involving special states.
                    matrix->mat_matrix[column]->cv_has_gaps = true;     // Might need a better way to handle this.
                }
                
                //memset(substring, 0, (matrix->max_states + 1) * sizeof(char));
                ++column;
            }
            else if (*current == '(' || *current == '{') {
                mfl_copy_multistate_subtoken_to_substring(&current, substring);
                strcpy(matrix->mat_matrix[column]->cv_character_cells[row], substring);
                //memset(substring, 0, (matrix->max_states + 1) * sizeof(char));
                ++column;
            }
            ++current;
        } while (column < num_chars && *current != '\n');
        //++current;
        ++row;
    } while (row < num_taxa);
    
    free(substring);
}

// Find and list the columns that have gaps in them (can be done simultaneously), setting the vector's cv_has_gaps flag
//      This now gives the data needed to set a list of characters with gaps (the first partitioning)

// Go through each cv_character_cell and convert its entry to the corresponding matrix cell.
//      This is where the conversion rule needs to be defined.

// Get the number of character types: unord, ord, dollo, irreversible, usertypes/sankoff

// Set up the correct number of partitions

// Copy the data into the partitions

// Apply the data to the tips of the tree.

/**
 Allocates all the memory required to store character data (pre- and post-conversion
 in the appropriate vectors. Does not populate the vectors with data--only creates
 the empty matrix, ready for populating with data.
 
 @param mfl_matrix_t* newmatrix, the Morphy-type matrix to be set up
 @param num_states the number of character states
 @param int the number of taxa
 @param int the number of characters
 @returns void
 */
void mfl_setup_new_empty_matrix(mfl_matrix_t *newmatrix, int num_states, int num_taxa, int num_chars)
{
    int i = 0;
    int j = 0;
    int charcells_size = num_states+1; // +1 for endline character;
    
    if (!num_states) {
        dbg_printf("ERROR in mfl_setup_new_empty_matrix(): cannot setup matrix without number of states\n");
        return;
    }
    
    for (i = 0; i < num_chars; ++i) {
        
        // Allocate memory for the vectors
        newmatrix->mat_matrix[i] = (mfl_character_vector_t*)malloc(sizeof(mfl_character_vector_t));
        if (!newmatrix->mat_matrix[i]) {
            dbg_printf("ERROR in mfl_setup_new_empty_matrix(): unable to allocate memory for mfl_character_vector_t's\n");
        }
        else {
            memset(newmatrix->mat_matrix[i], 0, sizeof(mfl_character_vector_t));
        }
        
        
        // Allocate and set up the char-type transformation series vectors
        newmatrix->mat_matrix[i]->cv_character_cells = (char**)malloc(num_taxa * sizeof(char*));
        
        if (!newmatrix->mat_matrix[i]->cv_character_cells) {
            dbg_printf("ERROR in mfl_setup_new_empty_matrix(): unable to allocate memory for cv_character_cells\n");
        }
        else {
            memset(newmatrix->mat_matrix[i]->cv_character_cells, 0, num_taxa * sizeof(char));
        }
        
        for (j = 0; j < num_taxa; ++j) {
            newmatrix->mat_matrix[i]->cv_character_cells[j] = (char*)malloc(charcells_size * sizeof(char));
            if (!newmatrix->mat_matrix[i]->cv_character_cells[j]) {
                dbg_printf("ERROR in mfl_setup_new_empty_matrix(): unable to allocate memory for a character cell\n");
            }
            else {
                memset(newmatrix->mat_matrix[i]->cv_character_cells[j], 0, charcells_size * sizeof(char));
            }
        }
        
        
        // Allocate the vector of bitwise state representations.
        newmatrix->mat_matrix[i]->cv_chardata = (mfl_charstate*)malloc(num_taxa * sizeof(sizeof(mfl_charstate)));
        if (!newmatrix->mat_matrix[i]->cv_chardata) {
            dbg_printf("ERROR in mfl_setup_new_empty_matrix(): unable to allocate memory for cv_chardata\n");
        }
        else {
            memset(newmatrix->mat_matrix[i]->cv_chardata, 0, num_taxa * sizeof(mfl_charstate));
        }
    }
    
    newmatrix->max_states = num_states;
}


/**
 Allocates memory for a new matrix and allocates memory for the list of
 character vectors pointed to by the mat_matrix variable. The returned matrix
 has no storage for the content of the vectors, just a list for pointers to
 that eventual data. Further setups are handled by mfl_setup_new_empty_matrix()
 @param int num_taxa: the number of input taxa
 @param int num_chars: the number of input characters
 @return mfl_matrix_t* the resulting empty matrix.
 */
mfl_matrix_t* mfl_create_mfl_matrix(int num_taxa, int num_chars)
{
    
    mfl_matrix_t *newmatrix = NULL;
    
    newmatrix = (mfl_matrix_t*)malloc(sizeof(mfl_matrix_t));
    if (!newmatrix) {
        dbg_printf("ERROR in mfl_create_mfl_matrix(): unable to allocate memory for new matrix\n");
        return NULL;
    }
    else {
        memset(newmatrix, 0, sizeof(mfl_matrix_t));
    }
    
    newmatrix->mat_num_taxa = num_taxa;
    newmatrix->mat_num_characters = num_chars;
    
    
    newmatrix->mat_matrix = (mfl_character_vector_t**)malloc(num_chars * sizeof(mfl_character_vector_t*));
    if (!newmatrix->mat_matrix) {
        dbg_printf("ERROR in mfl_create_mfl_matrix(): unable to allocate memory for new matrix\n");
        free(newmatrix);
        return NULL;
    }
    else {
        memset(newmatrix->mat_matrix, 0, num_chars * sizeof(mfl_character_vector_t*));
    }
    
    
    return newmatrix;
}


void mfl_destroy_character_cells(char **char_cells, int num_states, int num_taxa)
{
    int i = 0;
    
    for (i = 0; i < num_taxa; ++i) {
        free(char_cells[i]);
    }
}


void mfl_destroy_mfl_matrix(mfl_matrix_t *oldmatrix, int num_states, int num_taxa, int num_chars)
{
    int i = 0;
    
    for (i = 0; i < num_chars; ++i) {
        mfl_destroy_character_cells(oldmatrix->mat_matrix[i]->cv_character_cells, num_states, num_taxa);
        free(oldmatrix->mat_matrix[i]->cv_character_cells);
        free(oldmatrix->mat_matrix[i]->cv_chardata);
        free(oldmatrix->mat_matrix[i]);
    }
    
    free(oldmatrix->mat_matrix);
    free(oldmatrix);
}


/* Secondary input processing input */
/* Data and higher-level instructions might be passed to Morphy requiring further
 * processing. */

char *mfl_move_past_eq_sign(char *input)
{
    /* Move through a line and return the position of the first character after
     * an '=' sign; return original input position if no such symbol found. */
    char *eq_ptr = NULL;
    eq_ptr = input;
    
    do {
        if (*eq_ptr == '=') {
            ++eq_ptr;
            return eq_ptr;
        }
        ++eq_ptr;
    } while (*eq_ptr);
    
    return input;
}


void mfl_move_past_symbol(char** c, char symbol)
{
    
    do {
        if (**c == symbol) {
            break;
        }
        ++(*c);
    } while (**c);
}


bool mfl_check_nexus_matrix_dimensions(char *input_matrix, int input_num_taxa, int input_num_chars)
{
    /* An input matrix should not have inline taxon names. This function 
     * scans each row of the input matrix to determine whether or not the 
     * number of places in the row corresponds to input number of 
     * of characters. If the number exceeds the expected number of data
     * columns (num_input_chars), then it is inferred that taxon names or
     * other extraneous info are included in the matrix. */
    
    char* current = NULL;
    int matrix_size = 0;
    int expected_size = 0;
    bool dimensionsOK = true;
    
    expected_size = input_num_chars * input_num_taxa;
    
    current = input_matrix;
    
    do {
        if (mfl_is_valid_morphy_ctype(*current)) {
            ++matrix_size;
        }
        else if (*current == '(' || *current == '{') {
            mfl_move_in_nexus_multistate(&current);
            ++matrix_size;
        }
        
        if (*current == '\n') {
            mfl_move_past_symbol(&current, '\n');
        }
        ++current;
    } while (*current != ';');
    
    if (matrix_size > expected_size) {
        dimensionsOK = false;
        dbg_printf("ERROR in mfl_check_nexus_matrix_dimensions(): matrix dimensions exceed expected value based on NTAX and NCHAR input\n");
    }
    else if (matrix_size < expected_size) {
        dimensionsOK = false;
        dbg_printf("ERROR in mfl_check_nexus_matrix_dimensions(): matrix dimensions less than expected value based on NTAX and NCHAR input\n");
    }
    
    return dimensionsOK;
}


int mfl_check_state_number_support(char *datatype_list)
{
    
    /* Checks if the number of input states is supported by the default maximum
     * number of states. Returns the */
    
    int num_symbols = 0;
    char *dtype_ptr = NULL;
    
    // In case the input contains the tokens for Format and DataType as used in
    // Nexus files, this will move past the equals sign that begins the
    // assignment.
    dtype_ptr = mfl_move_past_eq_sign(datatype_list);
    
    while (dtype_ptr) {
        if (!isspace(*dtype_ptr)) {
            ++num_symbols;
        }
        ++dtype_ptr;
    }
    
    assert(num_symbols != 0);
    
    if (num_symbols > MORPHY_MAX_STATE_NUMBER) {
        dbg_printf("ERROR in mfl_check_number_support(): number of input states exceeds maximum allowed (%i).\n", MORPHY_MAX_STATE_NUMBER);
        dbg_printf("State list will be truncated\n");
        return MORPHY_MAX_STATE_NUMBER;
    }
    else {
        return num_symbols;
    }
    
}


void mfl_move_in_nexus_multistate(char **col)
{
    if (**col != '(' && **col != '{') {
        dbg_printf("WARNING in mfl_move_in_nexus_multistate(): *col is not start of multistate entry\n");
        return;
    }
    
    if (**col == '(') {
        do {
            ++(*col);
        } while (**col != ')');
    }
    else {
        do {
            ++(*col);
        } while (**col != '}');
    }
}


bool mfl_is_nexus_stop_position(char a)
{
    bool is_stop = false;
    
    if (a) {
        if (isspace(a) || a == '-' || a == ';' || a == ',') {
            is_stop = true;
        }
    }
    else {
        is_stop = true;
    }
    
    return is_stop;
}


void mfl_skip_spaces(char **current)
{
    if (**current == ' ') {
        do {
            ++(*current);
        } while (isspace(**current));
    }
}


int mfl_read_nexus_type_int(char **current)
{
    /* Read digit chars as an int until Nexus-type stop value is found. Return
     * the corresponding integer. */
    
    int nexus_int = 0;
    
    if (!isdigit(**current)) {
        mfl_move_current_to_digit(current);
    }
    
    // Read char until whitespace, a dash '-', semicolon, or end of string
    do {
        nexus_int = 10 * nexus_int;
        nexus_int = nexus_int + (**(current) - '0');
        ++(*current);
    } while (!mfl_is_nexus_stop_position(**current));
    
    mfl_skip_spaces(current);
    
    return nexus_int;
}


void mfl_set_include_value(int vectornum, int includeval, int* includes)
{
    /* In any set of bool include array, this sets the entry's value to the
     * desired true/false value */
    
    
    if (includeval) {
        includes[vectornum-1] = true;
    }
    else {
        includes[vectornum-1] = false;
    }
}


void mfl_set_include_range(int first, int last, int includeval, int* includes)
{
    /* Set all boolean values in includes to the true/false value indicated in
     * by includeval */
    
    int i = 0;
    
    for (i = first; i <= last; ++i) {
        mfl_set_include_value(i, includeval, includes);
    }
}


void mfl_move_current_to_digit(char** current)
{
    /* Increment the pointer to a digit character */
    
    if (!isdigit(**current)) {
        do {
            ++(*current);
        } while (!isdigit(**current));
    }
}


void mfl_set_inclusion_list(int* includes, int includeval, int listmax, char *subcommand)
{
    /* A generic function for setting a list of boolean values that is used
     * for including or excluding items from a list supplied by the user.
     * The list needs to be formatted in the Nexus style and must have 
     * numerical values.
     *
     * Reads the given subcommand, setting the positing in includes to the
     * true/false value specified by includeval. Converts the numeric tokens to
     * integers (1-based) which are used to index the include array (0-based). 
     * If a range of values is specified by the '-' character, then all 
     * positions in that range are set to includeval. Attempts to index out of 
     * list maximum are disallowed; ranges out of list maximum are truncated. */
    
    int first = 0;
    int last = 0; // first and last for identifier ranges (e.g. taxa 1-4 or
                  //characters 11-17)
    char *current = NULL;
    current = subcommand;
    
    if (!isdigit(*current)) {
        mfl_move_current_to_digit(&current);
    }
    
    do {
        first = 0;
        last = 0;
        if (isdigit(*current)) {
            // Read token as either a number or value range.
            first = mfl_read_nexus_type_int(&current);
        }
        else {
            ++current;
        }
        
        if (*current == '-') {
            mfl_move_current_to_digit(&current);
            last = mfl_read_nexus_type_int(&current);
        }
        
        if (first) {
            if (last) {
                if (last > listmax) {
                    last = listmax;
                    
                    dbg_printf("WARNING in mfl_set_inclusion_list(): Proposed range is\
                               outside of list maximum. Proposed range will be truncated at maximum: %i\n", listmax);
                }
                mfl_set_include_range(first, last, includeval, includes);
                
            }
            else {
                if (first <= listmax) {
                    mfl_set_include_value(first, includeval, includes);
                }
                else {
                    dbg_printf("WARNING in mfl_set_inclusion_list(): value will attempt to index outside of range. Ignorning proposed index\n");
                }
            }
        }
    } while (*current && *current != ';');
}


int* mfl_alloc_set_list(int num_chars)
{
    int *new_inclusion_list = NULL;
    
    new_inclusion_list = (int*)malloc(num_chars *sizeof(int));
    if (!new_inclusion_list) {
        dbg_printf("ERROR in mfl_alloc_character_inclusion_list(): unable to allocate memory\n");
    }
    else {
        memset(new_inclusion_list, MORPHY_DEFAULT_CHARACTER_INCLUDE, num_chars * sizeof(int));
    }
    
    return new_inclusion_list;
}


void mfl_free_set_list(bool *inclist)
{
    free(inclist);
}

/*!
 Takes a Nexus-formatted integer list and sets corresponding positions in the 
 input list to the value indicated by setval
 */
void mfl_read_nexus_style_list_subcmd(char *subcommand, int setval, int* list, int nelems)
{
    //int* list = NULL;
    
    //list = mfl_alloc_set_list(nelems);
    subcommand = mfl_move_past_eq_sign(subcommand);
    mfl_set_inclusion_list(list, setval, nelems, subcommand);
}

/* It may be possible to get rid of this function in favour of a more
 * general one. */
/*int * mfl_read_nexus_exset_subcmd(char *subcommand, int Nexus_NCHARS)
{
    bool* includelist = NULL;
    
    includelist = mfl_alloc_character_inclusion_list(Nexus_NCHARS);
    subcommand = mfl_move_past_eq_sign(subcommand);
    mfl_set_inclusion_list(includelist, false, Nexus_NCHARS, subcommand);
    
    return includelist;
}*/


int mfl_calculate_data_partitions_required(mfl_handle_t mfl_handle)
{
    
}