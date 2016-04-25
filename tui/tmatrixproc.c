/*
    tmatrixproc.c
    Morphy

 
    Tests for matrix processing
 */

#include "morphy.h"
#include "tuimfy.h"



int tui_test_matrix_processing(mfl_handle_s *mfl_handle)
{
    char wagners[] = "2-6 10 21-25";
    char costmat[] = "11-15";
    mfl_handle->ctypes_cmd[MFL_OPT_WAGNER] = wagners;
    mfl_handle->ctypes_cmd[MFL_OPT_COST_MATRIX] = costmat;
    
    // Get the partition
    
}