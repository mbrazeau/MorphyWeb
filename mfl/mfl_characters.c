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
 @discussion Checks that an input symbol is an allowed character state symbol, 
 not including '?' or the '-' (gap) character.
 @param c (char) an input symbol in the data matrix
 @return true or false
 */
bool mfl_is_valid_morphy_statesymbol(char c)
{
    if (isalnum(c)) {
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
 Looks up the symbol held by 'in' in the datype_converter array. The index
 within the array where the symbols match is found, plus the number of special 
 states reserved in Morphy determines the number of shifts to be used to convert 
 the symbol. Each character partition could have its own datype_converter and 
 its own way of shuffling the symbols. Also, if the Nexus standard is followed 
 and sybols list is supplied, then the order in which input symbols are declared
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


/*!
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


/*!
 @discussion Calculates the number of shifts required to convert a DNA base 
 symbol in an mfl_charstate variable.
 @param in (char) the input alpha character for a DNA base
 @return The number of left-shifts required to set the bit
 */
int mfl_shift_value_DNA_base(char in)
{
    if (in == 'a' || in == 'A') {
        return 1 + MORPHY_SPECIAL_STATE_PAD;
    }
    else if (in == 'c' || in == 'C') {
        return 2 + MORPHY_SPECIAL_STATE_PAD;
    }
    else if (in == 'g' || in == 'G') {
        return 3 + MORPHY_SPECIAL_STATE_PAD;
    }
    else if (in == 't' || in == 'T') {
        return 4 + MORPHY_SPECIAL_STATE_PAD;
    }
    else {
        dbg_eprintf("input not valid DNA base");
        return 0;
    }
}


/*!
 Converts a string of character states corresponding to a matrix cell containing
 multistate values. For each state symbole, the shift value of the data type 
 converter is calculated and that bit 'shifted in' by that amount and then added 
 to the output with a bitwise OR.
 @param xstates (char*) the string corresponding to the multistate cell.
 @param datype_converter (char*) the ordered array of symbols used to calculated 
 the number of shifts
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
    } while (*xstates);
    
    return xstate_bits;
}


/*!
 @discussion When a symbols list is not supplied and Morphy's default symbol 
 translation is in effect, this function sets the corresponding bit shift value
 required for setting the bit(s) corresponding to the symbol's value.
 @param xstates (char*) the string corresponding to symbolic (char-type) 
 representation of the state in the matrix.
 @return An mfl_charstate variable with the corresponding bit(s) set.
 */
mfl_charstate mfl_convert_using_default_rule(char *xstates)
{
    u_int64_t inbit = 0;
    int shift_val;
    mfl_charstate xstate_bits = 0;
    
    do {
        inbit = 1;
        shift_val = mfl_shift_value_DEFAULT_catdata(*xstates);
        inbit = inbit << shift_val;
        xstate_bits = xstate_bits | inbit;
        ++xstates;
    } while (*xstates);
    
    return xstate_bits;
}


/*!
 Sets the bits in an mfl_charstate variable to the appropriate value and returns
 the result.
 @param gapmethod (mfl_optimisation_t) the optimality method applied to the gap 
 character.
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
 states in a Nexus-derived "Format Symbols" command. A lookup into this will 
 return the required shift value.
 @param datype_converter (char*) the set of symbols to convert presented as a 
 string
 @param symbols_list (char*) the subtoken of symbols taken from input
 */
void mfl_set_datatype_converter_from_symbol_list(char* datype_converter, char* symbols_list)
{
    
    int i = 0;
    char *dtype_ptr = NULL;
    
    dtype_ptr = symbols_list;
    
    while (*dtype_ptr) {
        // There is a check here in case of accidental inclusion of unrecognised symbols
        if (!isspace(*dtype_ptr) && *dtype_ptr != '-' && *dtype_ptr != '?') {
            datype_converter[i] = *dtype_ptr;
            ++i;
        }
        ++dtype_ptr;
    }
}


/*!
 @discussion Does a linear search for a character in a list.
 @param key (char*) the character to be searched in the array.
 @param list (char*) the array to search
 @param listend (int *) the index of the last entry in the list
 @param listmax (int) the maximum size the list can have.
 @return Pointer to the key's match in the array, NULL if no match
 */
char* mfl_search_char_in_chartypes_array(char* key, char* list, int *listend, int listmax)
{
    // !!!: This function (and those it depends on) need testing for off-by-1 errors
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
            dbg_printf("\t Insufficient space in array for %i member \'%c\'. Exceeding max states: %lu.\n", *listend, *key, MORPHY_MAX_STATE_NUMBER);
            dbg_printf("Returning NULL\n\n");
            ++(*listend); // Exceed list max for error reporting.
            return NULL;
        }
    }
    
    return ptr;
}


/*!
 When no state symbols are specified by the user or input file and default 
 reading is in effect, the number of unique symbols (including gap) used can, 
 with some assumptions, be counted directly from the matrix. This function 
 attempts to do this.
 @param inputmatrix (char*) the states-only string corresponding to the data 
 matrix (can include braces or brackets for multistate taxa)
 @return the number of states in the data matrix
 */
int mfl_get_numstates_from_matrix(char *inputmatrix)
{
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


/*!
 @discussion Allocates an array used to assign the parsimony type for each 
 column (transformation series)in the data matrix. Sets the default parsimony
 type for all characters.
 @param input_num_chars (int) the number of characters in the dataset
 @return a pointer to an array of integers (enum type mfl_parsimony_t) all set
 to the default parsimony method.
 */
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


/*!
 @discussion Wrapper on the free function for freeing the chartype list.
 @param ctype_list (mfl_parsimony_t*) the character type-by transformation 
 series array
 */
void mfl_destroy_chartype_list(mfl_parsimony_t *ctype_list)
{
    free(ctype_list);
}


/*!
 @discussion Sets up a new chartypes list by parsing the list of character 
 (parsimony) types in the mfl_handle.
 @param mfl_handle (mfl_handle_s*) pointer to the mfl_handle that was set up 
 by the interface functions
 @return an array of mfl_parsimony_t values set to the parsimony type indicated
 in the character types command in the handle.
 */
mfl_parsimony_t* mfl_get_chartypes_list(const mfl_handle_s* mfl_handle)
{
    int i;
    mfl_parsimony_t* ctype_setters = mfl_alloc_chartype_list(mfl_handle->n_chars);
    
    if (mfl_handle->ctype_setters) {
        memccpy(ctype_setters, mfl_handle->ctype_setters, mfl_handle->n_chars, sizeof(mfl_parsimony_t));
    }
    else if (mfl_handle->ctypes_cmd) {
        if (!mfl_handle->n_ctypes) {
            dbg_printf("All characters are of type %s\n", MORPHY_DEFAULT_PARSIMONY_NAME);
        }
        else {
            // Parse each ctype command
            for (i = 0; i < MFL_OPT_MAX; ++i) {
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


/*!
 @discussion Converts a symbol to a set bit in an mfl_charstate variable
 @param c (char*) the input character from the data matrix
 @param datype_converter (char*) the optional symbols list used to determine the
 ordering of symbol values
 @param gaprule (mfl_gap_t) the parameter for how the gap character will be 
 handled
 @return an mfl_charstate variable with the corresponding bit(s) set
 */
mfl_charstate mfl_standard_conversion_rule(char *c, char* datype_converter, mfl_gap_t gaprule)
{
    if (mfl_is_valid_morphy_statesymbol(*c)) {
        if (datype_converter) {
            return mfl_convert_nexus_multistate(c, datype_converter);
        }
        else {
            return mfl_convert_using_default_rule(c);
        }
    }
    else if (*c == '?') {
        if (gaprule == MFL_GAP_INAPPLICABLE || gaprule == MFL_GAP_MISSING_DATA) {
            return MORPHY_MISSING_DATA_BITWISE;
        }
        else {
            return MORPHY_MISSING_DATA_BITWISE | 1;
        }
    }
    else if (*c == '-') {
        if (gaprule == MFL_GAP_INAPPLICABLE) {
            // Gaps are translated into 1's. Missing data is inverse.
            return mfl_convert_gap_character(gaprule);
        }
        else if (gaprule == MFL_GAP_MISSING_DATA) {
            // Gaps are translated into missing data.
            return mfl_convert_gap_character(gaprule);
        }
        else if (gaprule == MFL_GAP_NEWSTATE) {
            // Gaps are translated into 1's. Missing data is not inverse.
            return mfl_convert_gap_character(gaprule);
        }
    }
    
    dbg_eprintf("unable to convert input character");
    return 0;
}


/*mfl_charstate mfl_CUSTOM_CONVERSION_RULES(char *c, mfl_gap_t gaprule)
{
    __TEMPLATE__FOR__CUSTOM__CONVERSION__RULE
}*/


/*!
 An interface for the conversion rule used to turn input characters into set 
 bits in an mfl_charstate variable. If ever there is a need to write different
 conversion rules, they can be set in here.
 @param parsimtype (mfl_parsimony_t) the type of parsimony function used
 @return Pointer to the wrapper on the conversion rules to be applied to the 
 data.
 */
mfl_char2bit_fn mfl_fetch_conversion_rule(mfl_parsimony_t parsimtype)
{
    mfl_char2bit_fn ret = NULL;
    
    switch (parsimtype) {
        
        case MFL_OPT_FITCH:
            ret = mfl_standard_conversion_rule;
            break;
        case MFL_OPT_WAGNER:
            ret = mfl_standard_conversion_rule;
            break;
        case MFL_OPT_DOLLO_UP:
            ret = mfl_standard_conversion_rule;
            break;
        case MFL_OPT_DOLLO_DN:
            ret = mfl_standard_conversion_rule;
            break;
        case MFL_OPT_IRREVERSIBLE_UP:
            ret = mfl_standard_conversion_rule;
            break;
        case MFL_OPT_IRREVERSIBLE_DN:
            ret = mfl_standard_conversion_rule;
            break;
        case MFL_OPT_COST_MATRIX:
            ret = mfl_standard_conversion_rule;
            break;
        /*case ___CUSTOM_TYPE___:
             ret = __POINTER_TO_CUSTOM_CONVERSION_ROUTINE;
             break;*/
        default:
            break;
    }
    
    return ret;
}


/*!
 @discussion Sets the conversion rule function pointer in each column vector in 
 the matrix according to the specified parsimony method. Mostly, this step is 
 not necessary, but provides an interface in case a programmer wished to modify
 Morphy's rules for converting characters given a custom parsimony/optimatlity 
 type.
 @param matrix (mfl_matrix_t*) The internal matrix derived from input data.
 */
void mfl_apply_conversion_rules(mfl_matrix_t *matrix)
{
    int i = 0;
    int num_chars = matrix->mat_num_characters;
    
    for (i = 0; i < num_chars; ++i) {
        matrix->mat_matrix[i]->cv_conversion_rule = mfl_fetch_conversion_rule(matrix->mat_matrix[i]->cv_parsim_method);
    }
}


/*!
 @discussion Converts the symbolic char-type representations of the characters
 to mfl_charstates in a single character column vector.
 @param cv (mfl_character_vector_t*) a single column in the data matrix.
 @param handle (const mfl_handle_s*) the interface handle that contains all of 
 the analysis parameters.
 */
void mfl_convert_charcells_to_mfl_charstates(mfl_character_vector_t* cv, const mfl_handle_s* handle)
{
    int i = 0;
    int num_taxa = 0;
    char* dataconverter = NULL;
    mfl_gap_t localgapmethod = handle->gap_method;
    
    if (cv->cv_num_gaps < 3 && handle->gap_method == MFL_GAP_INAPPLICABLE) {
        localgapmethod = MFL_GAP_MISSING_DATA;
    }
    
    if (handle->format_symbols) {
        dataconverter = (char*)malloc((handle->n_symbols + 1) * sizeof(char));
        if (!dataconverter) {
            dbg_eprintf("unable to allocate memory for new data converter");
        }
        else {
            memset(dataconverter, 0, (handle->n_symbols + 1) * sizeof(char));
        }
        mfl_set_datatype_converter_from_symbol_list(dataconverter, handle->format_symbols);
    }
    
    num_taxa = handle->n_taxa;
    
    for (i = 0; i < num_taxa; ++i) {
        cv->cv_chardata[i] = cv->cv_conversion_rule(cv->cv_character_cells[i], dataconverter, localgapmethod);
    }
    
    free(dataconverter);
}


/*!
 @discussion Converts all of the character column vectors from their char-type
 symbolic representations (i.e. the input) to their corresponding set bits in
 their respective mfl_charstate arrays in a matrix.
 @param m (mfl_matrix_t*) the internal matrix based on input data.
 @param handle (const mfl_handle_s*) the interface handle that contains all of 
 the analysis parameters.
 */
void mfl_convert_all_characters_to_charstates(mfl_matrix_t* m, const mfl_handle_s* handle)
{
    int i = 0;
    int num_chars = m->mat_num_characters;
    
    for (i = 0; i <num_chars; ++i) {
        mfl_convert_charcells_to_mfl_charstates(m->mat_matrix[i], handle);
    }
}



/*
 *
 * Functions for parsing the input matrix string into substrings. 
 *
 */

/*!
 @discussion Copies the state symbol of a single-state cell in the input matrix
 to the substring supplied and appends a terminal null. The result is a single-
 character string. Used for creating the internal data matrix.
 @param singleton (char) the single character state symbole
 @param substring (char*) the substring into which the character will be 
 copied at the first position. Substring must be pre-sized to fit the input.
 */
void mfl_copy_singleton_subtoken_to_substring(char singleton, char* substring)
{
    *substring = singleton;
    substring[1] = '\0';
}


/*!
 @discussion Copies the state symbols in a multistate cell in the input matrix
 into the supplied substring. Used for creating the internal data matrix.
 @param xstatetoken (char**) pointer to the substring within the input data
 string corresponding to the start of the multistate cell
 @param substring (char*) the substring into which the character symbols will be
 copied at the first position. Substring must be pre-sized to fit the input.
 */
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


/*!
 @discussion Populates the char-type array of strings for the state symbols 
 within the internal data matrix. These later get translated into mfl_charstate
 variables.
 @param matrix (mfl_matrix_t*) the internal matrix being set up when this
 routine is called.
 @param input_data_matrix (char*) the data matrix supplied as a string of
 character state symbols only (no taxon names).
 @param num_chars (int) the number of characters (columns) in the matrix
 @param num_taxa (int) the number of taxa (rows) in the matrix
 */
void mfl_populate_chartype_character_vectors(mfl_matrix_t *matrix, char *input_data_matrix, int num_chars, int num_taxa)
{
    int column = 0;
    int row = 0;
    char *current = NULL;
    char *substring = NULL;
    
    current = input_data_matrix;
    
    substring = (char*)malloc((matrix->mat_max_states + 1) * sizeof(char));
    if (!substring) {
        dbg_printf("ERROR in mfl_populate_chartype_character_vector(): unable to allocate memory for character substring\n");
        return;
    }
    else {
        memset(substring, 0, (matrix->mat_max_states + 1) * sizeof(char));
    }
    
    do {
        // Collect the input matrix substring
        column = 0;
        do {
            if (mfl_is_valid_morphy_ctype(*current)) {
                mfl_copy_singleton_subtoken_to_substring(*current, substring);
                strcpy(matrix->mat_matrix[column]->cv_character_cells[row], substring);
                ++column;
            }
            else if (*current == '(' || *current == '{') {
                mfl_copy_multistate_subtoken_to_substring(&current, substring);
                strcpy(matrix->mat_matrix[column]->cv_character_cells[row], substring);
                ++column;
            }
            ++current;
        } while (column < num_chars && *current != '\n');
        ++row;
    } while (row < num_taxa);
    
    free(substring);
}


/*!
 Allocates all the memory required to store character data (pre- and post-
 conversion in the appropriate vectors. Does not populate the vectors with 
 data--only creates the empty matrix, ready for populating with data.
 
 @param newmatrix (mfl_matrix_t*) the Morphy-type matrix to be set up
 @param num_states (int) the number of character states
 @param num_taxa (int) the number of taxa
 @param num_chars (int) the number of characters
 */
void mfl_setup_new_empty_matrix(mfl_matrix_t *newmatrix, int num_states, int num_taxa, int num_chars)
{
    int i = 0;
    int j = 0;
    int charcells_size = num_states+1; // +1 for terminal null character;
    
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
        newmatrix->mat_matrix[i]->cv_num_states = num_states;
        
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
    
    newmatrix->mat_max_states = num_states;
}


/*!
 Allocates memory for a new matrix and allocates memory for the list of
 character vectors pointed to by the mat_matrix variable. The returned matrix
 has no storage for the content of the vectors, just a list for pointers to
 that eventual data. Further setups are handled by mfl_setup_new_empty_matrix()
 @param num_taxa (int) the number of input taxa
 @param num_chars (int) the number of input characters
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
    
    
    newmatrix->mat_matrix = (mfl_character_vector_t**)mfl_malloc(num_chars * sizeof(mfl_character_vector_t*), 0);
    /*if (!newmatrix->mat_matrix) {
        dbg_printf("ERROR in mfl_create_mfl_matrix(): unable to allocate memory for new matrix\n");
        free(newmatrix);
        return NULL;
    }
    else {
        memset(newmatrix->mat_matrix, 0, num_chars * sizeof(mfl_character_vector_t*));
    }*/
    
    
    return newmatrix;
}


void mfl_destroy_character_cells(char **char_cells, int num_states, int num_taxa)
{
    int i = 0;
    
    for (i = 0; i < num_taxa; ++i) {
        free(char_cells[i]);
    }
}


void mfl_destroy_mfl_matrix(mfl_matrix_t *oldmatrix, int num_taxa, int num_chars)
{
    int i = 0;
    
    for (i = 0; i < num_chars; ++i) {
        mfl_destroy_character_cells(oldmatrix->mat_matrix[i]->cv_character_cells, oldmatrix->mat_matrix[i]->cv_num_states, num_taxa);
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
        dbg_printf("ERROR in mfl_check_number_support(): number of input states exceeds maximum allowed (%lu).\n", MORPHY_MAX_STATE_NUMBER);
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
        includes[vectornum-1] = includeval;
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


/*!
 A generic function for setting a list of integer values that is used for 
 applying enumerated rules to a list supplied by the user. The list needs to be 
 formatted in the Nexus style and must have numerical values.
 
 Reads the given subcommand, setting the positing in includes to the true/false 
 value specified by includeval. Converts the numeric tokens to integers 
 (1-based) which are used to index the include array (0-based). If a range of 
 values is specified by the '-' character, then all positions in that range are 
 set to includeval. Attempts to index out of list maximum are disallowed; ranges 
 out of list maximum are truncated.
 @param includes (int*) the array with each position corresponding to the 
 members of some other list (e.g. the input characters or taxa).
 @param includeval (int) the value to be set according to the command.
 @param listmax (int) the size of the list.
 @param subcommand (char*) the input string that is to be parsed.
 */
void mfl_set_inclusion_list(int* includes, int includeval, int listmax, char *subcommand)
{
    
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
                    dbg_printf("WARNING in mfl_set_inclusion_list(): value will attempt to index outside of range. Ignorning proposed index.\n");
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
 Takes a Nexus-formatted integer list and parses it to set the corresponding 
 positions in the input list to the value indicated by setval parameter. Thus, 
 setval should correspond to the instructions of the command and the subcommand 
 is the character string giving the numeric list of characters affected by the
 command.
 @param subcommand (char*) the string supplying the list of character affected
 by the command.
 @param setval (int) the corresponding subcommand corresponding to an enumerated
 data type.
 @param list (int*) the array of set values that will be setto the corresponding
 input setval.
 @param nelems (int) the length of the list parameter.
 */
void mfl_read_nexus_style_list_subcmd(char *subcommand, int setval, int* list, int nelems)
{
    subcommand = mfl_move_past_eq_sign(subcommand);
    mfl_set_inclusion_list(list, setval, nelems, subcommand);
}



void mfl_count_gaps_in_each_character(mfl_matrix_t* matrix)
{
    int i = 0;
    int j = 0;
    int num_gaps = 0;
    int numchar = matrix->mat_num_characters;
    int numtax = matrix->mat_num_taxa;
    
    for (i = 0; i < numchar; ++i) {
        num_gaps = 0;
        for (j = 0; j < numtax; ++j) {
            if (strchr(matrix->mat_matrix[i]->cv_character_cells[j], '-')) {
                ++num_gaps;
            }
        }
        matrix->mat_matrix[i]->cv_num_gaps = num_gaps;
    }
}



void mfl_set_cv_chartypes(mfl_matrix_t* matrix, const mfl_handle_s* mfl_handle, mfl_parsimony_t* chartypes)
{
    int i = 0;
    
    for (i = 0; i < mfl_handle->n_chars; ++i) {
        matrix->mat_matrix[i]->cv_parsim_method = chartypes[i];
    }
}


/*!
 Creates the array of character column vectors corresponding to the
 transformation series of the input dataset, sets all character type values and
 the pointers to the conversion rule for that character type. All conversion
 rules are applied and the appropriate treatment of gaps as either missing data 
 or a (special) state
 @param mfl_handle (mfl_handle_s*) the analysis parameters set by the user 
 interface
 @returns pointer to a new matrix.
 */
mfl_matrix_t* mfl_create_internal_data_matrix(const mfl_handle_s* mfl_handle)
{
    int i = 0;
    int num_states = 0;
    char* input_chardata = NULL;
    mfl_parsimony_t* chartypes = NULL;
    
    mfl_matrix_t* new_inmatrix = mfl_create_mfl_matrix(mfl_handle->n_taxa, mfl_handle->n_chars);
    
    // Get the matrix string from the handle.
    input_chardata = mfl_handle->input_data;
    
    // Perform checks on input data.
    // TODO: Wrap more checks into this function
    mfl_check_nexus_matrix_dimensions(input_chardata, mfl_handle->n_taxa, mfl_handle->n_chars);
    
    if (mfl_handle->n_symbols && mfl_handle->format_symbols) {
        num_states = mfl_handle->n_symbols;
        // TODO: do input checks.
    }
    else {
        dbg_printf("No symbols list supplied. Attempting to estimate state number from input matrix.\n");
        // If no value supplied, then attempt to get number of states
        // from the input matrix directly.
        num_states = mfl_get_numstates_from_matrix(input_chardata);
    }
    
    mfl_setup_new_empty_matrix(new_inmatrix, num_states, mfl_handle->n_taxa, mfl_handle->n_chars);
    
    mfl_populate_chartype_character_vectors(new_inmatrix, input_chardata, mfl_handle->n_chars, mfl_handle->n_taxa);
    
    if (mfl_handle->gap_method == MFL_GAP_INAPPLICABLE) {
        mfl_count_gaps_in_each_character(new_inmatrix);
    }
    
    chartypes = mfl_get_chartypes_list(mfl_handle);
    mfl_set_cv_chartypes(new_inmatrix, mfl_handle, chartypes);
    
    mfl_apply_conversion_rules(new_inmatrix);
    
    mfl_convert_all_characters_to_charstates(new_inmatrix, mfl_handle);
    
    // TODO: This needs to be updated to be modified to function on user input. This is temporary default.
    new_inmatrix->mat_included_characters = (bool*)mfl_malloc(new_inmatrix->mat_num_characters, true);
    
    
    if (!mfl_handle->int_weights) {
        if (!mfl_handle->real_weights) {
            new_inmatrix->mat_int_weights = (int*)mfl_malloc(new_inmatrix->mat_num_characters * sizeof(int), 0);
            for (i = 0; i < new_inmatrix->mat_num_characters; ++i) {
                new_inmatrix->mat_int_weights[i] = (int)MORPHY_DEFAULT_WEIGHT;
            }
        }
        else {
            new_inmatrix->mat_real_weights = mfl_handle->real_weights;
        }
    }
    else {
        new_inmatrix->mat_int_weights = mfl_handle->int_weights;
    }
    
    /*int i;
    dbg_printf("Printing chartypes array:\n");
    for (i = 0; i < mfl_handle->n_chars; ++i) {
        dbg_printf("%i ", chartypes[i]);
    }

    dbg_printf("\n");*/
    
    free(chartypes);
    
    return new_inmatrix;
}



mfl_datapartition_t* mfl_alloc_empty_datapartition_t(void)
{
    mfl_datapartition_t* newdatapart = (mfl_datapartition_t*)malloc(sizeof(mfl_datapartition_t));
    if (!newdatapart) {
        dbg_eprintf("unable to allocate memory for data partition");
        return NULL;
    }
    else {
        memset(newdatapart, 0, sizeof(mfl_datapartition_t));
    }
    
    return newdatapart;
}


int mfl_count_num_partitions_required(mfl_joint_character_t* jntchars, mfl_matrix_t* m, mfl_gap_t gaprule)
{
    int i = 0;
    int numdataparts = 0;
    
    for (i = 0; i < m->mat_num_characters; ++i) {
        if (m->mat_matrix[i]->cv_num_gaps > 2 && gaprule == MFL_GAP_INAPPLICABLE) {
            if (jntchars[m->mat_matrix[i]->cv_parsim_method].jpt_num_winapplic == 0) {
                jntchars[m->mat_matrix[i]->cv_parsim_method].jpt_inapplic_part_num = numdataparts;
                ++numdataparts;
            }
            
            m->mat_matrix[i]->cv_partition_destination = jntchars[m->mat_matrix[i]->cv_parsim_method].jpt_inapplic_part_num;
            
            jntchars[m->mat_matrix[i]->cv_parsim_method].jpt_num_winapplic += 1;
        }
        else {
            if (jntchars[m->mat_matrix[i]->cv_parsim_method].jpt_num_allapplic == 0) {
                jntchars[m->mat_matrix[i]->cv_parsim_method].jpt_applic_part_num = numdataparts;
                ++numdataparts;
            }
            
            m->mat_matrix[i]->cv_partition_destination = jntchars[m->mat_matrix[i]->cv_parsim_method].jpt_applic_part_num;
            
            jntchars[m->mat_matrix[i]->cv_parsim_method].jpt_num_allapplic += 1;
        }
    }
    
    dbg_printf("Using the joint ctypes: %i\n", numdataparts);
    
    return numdataparts;
}


void mfl_set_datapart_params_fitch(mfl_datapartition_t* d, bool has_inapplic)
{
    d->part_optimisation_method = MFL_OPT_FITCH;
    d->part_has_inapplicables = has_inapplic;
    if (has_inapplic) {
        d->part_char_is_directed = true;
    }
    else {
        d->part_char_is_directed = false;
    }
}


void mfl_set_datapart_params_wagner(mfl_datapartition_t* d, bool has_inapplic)
{
    d->part_optimisation_method = MFL_OPT_WAGNER;
    d->part_has_inapplicables = has_inapplic;
    if (has_inapplic) {
        d->part_char_is_directed = true;
    }
    else {
        d->part_char_is_directed = false;
    }
}


void mfl_set_datapart_params_dollo(mfl_datapartition_t* d, bool has_inapplic, bool up)
{
    if (up) {
        d->part_optimisation_method = MFL_OPT_DOLLO_UP;
    }
    else {
        d->part_optimisation_method = MFL_OPT_DOLLO_DN;
    }
    d->part_has_inapplicables = has_inapplic;
    d->part_char_is_directed = true;
}


void mfl_set_datapart_params_irrev(mfl_datapartition_t* d, bool has_inapplic, bool up)
{
    if (up) {
        d->part_optimisation_method = MFL_OPT_IRREVERSIBLE_UP;
    }
    else {
        d->part_optimisation_method = MFL_OPT_IRREVERSIBLE_DN;
    }
    d->part_has_inapplicables = has_inapplic;
    d->part_char_is_directed = true;
}


void mfl_set_datapart_params_costmatrx(mfl_datapartition_t* d, const mfl_handle_s* handle, bool has_inapplic)
{
    d->part_optimisation_method = MFL_OPT_COST_MATRIX;
    d->part_has_inapplicables = has_inapplic;
}


void mfl_set_datapart_params(mfl_datapartition_t* d, mfl_parsimony_t opt_t, bool has_inapplic, const mfl_handle_s* handle)
{
    switch (opt_t) {
        case MFL_OPT_FITCH:
            mfl_set_datapart_params_fitch(d, has_inapplic);
            break;
        case MFL_OPT_WAGNER:
            mfl_set_datapart_params_wagner(d, has_inapplic);
            break;
        case MFL_OPT_DOLLO_UP:
            mfl_set_datapart_params_dollo(d, has_inapplic, true);
            break;
        case MFL_OPT_DOLLO_DN:
            mfl_set_datapart_params_dollo(d, has_inapplic, false);
            break;
        case MFL_OPT_IRREVERSIBLE_UP:
            mfl_set_datapart_params_irrev(d, has_inapplic, true);
            break;
        case MFL_OPT_IRREVERSIBLE_DN:
            mfl_set_datapart_params_irrev(d, has_inapplic, false);
            break;
        case MFL_OPT_COST_MATRIX:
            mfl_set_datapart_params_costmatrx(d, handle, has_inapplic);
            break;
        case MFL_OPT_MAX:
            dbg_eprintf("unsupported optimality type");
            assert(0);
        default:
            break;
    }
}


int mfl_fetch_parsimony_fxn(mfl_datapartition_t* part, mfl_parsimony_t parsim_type, bool gapasinapplic)
{
    int ret = 0;
    
    switch (parsim_type)
    {
        case MFL_OPT_FITCH:
            if (gapasinapplic) {
                // TODO: Fix this before doing any work on inapplicable characters
                part->part_downpass_full    = &mfl_first_fitch_na_downpass;
                part->part_uppass_full      = &mfl_first_fitch_na_uppass;
                part->part_NAdownpass_full  = &mfl_second_fitch_na_downpass;
                part->part_NAuppass_full    = &mfl_second_fitch_na_uppass;
                part->part_NAlocal          = &mfl_test_fitch_na_local;
                ret = 0;
            }
            else {
                part->part_downpass_full    = &mfl_fitch_downpass_binary_node;
                part->part_uppass_full      = &mfl_fitch_uppass_binary_node;
                part->part_local            = &mfl_test_fitch_local;
                ret = 0;
            }
            break;
            
            /*case MFL_OPT_WAGNER:
             ret = mfl_wagner_downpass_binary_node;
             break;
             
             case MFL_OPT_DOLLO:
             ret = mfl_dollo_downpass_binary_node;
             break;
             
             case MFL_OPT_IRREVERSIBLE:
             ret = mfl_irreversible_downpass_binary_node;
             break;
             
             case MFL_OPT_COST_MATRIX:
             ret = mfl_costmatrix_downpass_binary_node;
             break;*/
            
        default:
            ret = -1;
            break;
    }
    
    return  ret;
}

void mfl_copy_column_into_partition(mfl_datapartition_t* prt, mfl_character_vector_t* cv, int intwt, int index, int num_rows)
{
    int i = 0;
    int n = 0;
    
    n = prt->part_n_chars_included;
    
    for (i = 0; i < num_rows; ++i) {
        
        prt->part_matrix[n + (i * prt->part_n_chars_max)] = cv->cv_chardata[i];
    }
    
    prt->part_char_indices[prt->part_n_chars_included] = index;
    // TODO: some more elegant way of handling weights, especially where real values might be used.
    prt->part_int_weights[prt->part_n_chars_included] = intwt;
    ++prt->part_n_chars_included;
}


void mfl_populate_all_character_partitions(mfl_partition_set_t* ptset, mfl_gap_t gapmethod, mfl_matrix_t* m)
{
    int i = 0;
    //int j = 0;
    int numparts = ptset->ptset_n_parts;
    bool inapplicable = false;
    mfl_parsimony_t ptype;
    
    for (i = 0; i < m->mat_num_characters; ++i) {
        if (m->mat_included_characters[i]) {
            mfl_copy_column_into_partition(ptset->ptset_partitions[m->mat_matrix[i]->cv_partition_destination],
                                           m->mat_matrix[i], m->mat_int_weights[i], i,
                                           m->mat_num_taxa);
        }
    }
    
    for (i = 0; i < numparts; ++i) {
        
        inapplicable = false;
        ptype = ptset->ptset_partitions[i]->part_optimisation_method;
        
        if (ptset->ptset_partitions[i]->part_has_inapplicables) {
            if (gapmethod == MFL_GAP_INAPPLICABLE) {
                inapplicable = true;
            }
            ptset->ptset_partitions[i]->part_activestates = (mfl_charstate*)mfl_malloc(ptset->ptset_partitions[i]->part_n_chars_included * sizeof(mfl_charstate), 0);
        }
    
        mfl_fetch_parsimony_fxn(ptset->ptset_partitions[i], ptype, inapplicable);
        
        /*
        ptset->ptset_partitions[i]->part_downpass_full = mfl_fetch_downpass_parsimony_fxn(ptype, inapplicable, true);
        ptset->ptset_partitions[i]->part_downpass_partial = mfl_fetch_downpass_parsimony_fxn(ptype, inapplicable, false);
        ptset->ptset_partitions[i]->part_uppass_full = mfl_fetch_uppass_parsimony_fxn(ptype, inapplicable, true);
        ptset->ptset_partitions[i]->part_uppass_partial = mfl_fetch_uppass_parsimony_fxn(ptype, inapplicable, false);*/
    }
    
}



int mfl_compare_dataparts_by_index(const void* p1, const void* p2)
{
    
    mfl_datapartition_t* p_1 = *(mfl_datapartition_t**)p1;
    mfl_datapartition_t* p_2 = *(mfl_datapartition_t**)p2;
    
    return (int)(p_1->part_index - p_2->part_index);
}


int mfl_compare_dataparts_by_ctype(const void* p1, const void* p2)
{
    
    mfl_datapartition_t* p_1 = *(mfl_datapartition_t**)p1;
    mfl_datapartition_t* p_2 = *(mfl_datapartition_t**)p2;
    
    if (p_1->part_optimisation_method == p_2->part_optimisation_method) {
        if (p_1->part_has_inapplicables || p_2->part_has_inapplicables) {
            if (p_1->part_has_inapplicables) {
                return 1;
            }
            else {
                return -1;
            }
        }
        else {
            return 0;
        }
    }
    else {
        return (int)(p_1->part_optimisation_method - p_2->part_optimisation_method);
    }
}


mfl_partition_set_t* mfl_create_data_partitions_set(mfl_matrix_t* matrix, const mfl_handle_s* handle)
{
    int i = 0;
    int part_i = 0;
    int numparts = 0;

    mfl_joint_character_t jntchars[MFL_OPT_MAX] = {{0}};
    
    numparts = mfl_count_num_partitions_required(jntchars, matrix, handle->gap_method);
    
    mfl_partition_set_t* dataparts = (mfl_partition_set_t*)mfl_malloc(numparts * sizeof(mfl_partition_set_t), 0);
    
    dataparts->ptset_n_parts = numparts;
    
    dataparts->ptset_partitions = (mfl_datapartition_t**)mfl_malloc(numparts * sizeof(mfl_datapartition_t*), 0);
    
    // Allocate the dataparts
    for (i = 0; i < numparts; ++i) {
        dataparts->ptset_partitions[i] = mfl_alloc_empty_datapartition_t();
    }
    
    // OPTIMISATION: Can probably try to make cost matrix characters come last out of all
    // Set the params for the partitions without inapplicables
    for (i = 0; i < MFL_OPT_MAX; ++i) {
        if (jntchars[i].jpt_num_allapplic) {
            mfl_set_datapart_params(dataparts->ptset_partitions[part_i], (mfl_parsimony_t)i, false, handle);
            dataparts->ptset_partitions[part_i]->part_n_chars_max = jntchars[i].jpt_num_allapplic;
            dataparts->ptset_partitions[part_i]->part_char_indices = (int*)mfl_malloc(dataparts->ptset_partitions[part_i]->part_n_chars_max * sizeof(int), 0);
            dataparts->ptset_partitions[part_i]->part_index = jntchars[i].jpt_applic_part_num;
            ++part_i;
        }
    }
    
    // Set the params for the partitions with inapplicables
    for (i = 0; i < MFL_OPT_MAX; ++i) {
        if (jntchars[i].jpt_num_winapplic) {
            mfl_set_datapart_params(dataparts->ptset_partitions[part_i], (mfl_parsimony_t)i, true, handle);
            dataparts->ptset_partitions[part_i]->part_n_chars_max = jntchars[i].jpt_num_winapplic;
            dataparts->ptset_partitions[part_i]->part_char_indices = (int*)mfl_malloc(dataparts->ptset_partitions[part_i]->part_n_chars_max * sizeof(int), 0);
            dataparts->ptset_partitions[part_i]->part_index = jntchars[i].jpt_inapplic_part_num;
            ++part_i;
        }
    }
    
    
    for (i = 0; i < numparts; ++i) {
        dataparts->ptset_partitions[i]->part_int_weights = (int*)mfl_malloc(dataparts->ptset_partitions[i]->part_n_chars_max * sizeof(int), 0);
        dataparts->ptset_partitions[i]->part_matrix = (mfl_charstate*)mfl_malloc((dataparts->ptset_partitions[i]->part_n_chars_max * matrix->mat_num_taxa)* sizeof(mfl_charstate), 0);
    }
    
    
    qsort(dataparts->ptset_partitions, dataparts->ptset_n_parts, sizeof(mfl_datapartition_t*), &mfl_compare_dataparts_by_index);
    
    mfl_populate_all_character_partitions(dataparts, handle->gap_method, matrix);
    
    qsort(dataparts->ptset_partitions, dataparts->ptset_n_parts, sizeof(mfl_datapartition_t*), &mfl_compare_dataparts_by_ctype);
    
    mfl_destroy_mfl_matrix(matrix, handle->n_taxa, handle->n_chars);
    
    return dataparts;

}

mfl_partition_set_t* mfl_generate_search_data(const mfl_handle_s* handle)
{
    mfl_matrix_t* intmatrix = mfl_create_internal_data_matrix(handle);
    return mfl_create_data_partitions_set(intmatrix, handle);
}

void mfl_destroy_partition_set(mfl_partition_set_t* ptset)
{
    int i = 0;
    int part_i = 0;
    int numparts = 0;

    numparts = ptset->ptset_n_parts;
    
    // Dellocate the dataparts
    for (i = 0; i < numparts; ++i) {
        free(ptset->ptset_partitions[i]->part_char_indices);
        free(ptset->ptset_partitions[i]->part_int_weights);
        free(ptset->ptset_partitions[i]->part_matrix);
        free(ptset->ptset_partitions[i]);
        ptset->ptset_partitions[i] = NULL;
    }
    
    free(ptset->ptset_partitions);
    ptset->ptset_partitions = NULL;
    free(ptset);
}
