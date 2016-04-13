/*
    tmatrixproc.c
    Morphy

 
    Tests for matrix processing
 */

#include "morphy.h"
#include "tuimfy.h"


int tui_test_matrix_processing(const char *filename)
{
    FILE* inputfile;
    
    inputfile = fopen(filename, "r");
    
    mfl_handle_t testhandle;
    testhandle = mfl_create_handle();
    
    
    
    mfl_destroy_handle(testhandle);
}