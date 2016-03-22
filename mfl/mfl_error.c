//
//  mfl_error.c
//  Morphy
//
//  Created by mbrazeau on 07/03/2016.
//  Copyright © 2016 Imperial College London. All rights reserved.
//

#include "morphy.h"

int mfl_error_report_to_stderr(char *error_message, int exit_code); // Move to morphy.h

int mfl_error_report_to_stderr(char *error_message, int exit_code)
{
    
    int fputs_return = 0;
    
    fputs_return = fputs(error_message, stderr);
    
    if (fputs_return >= 0) {
        return 0;
    }
    else {
        return 1; // Make this something else?
    }
    
}

