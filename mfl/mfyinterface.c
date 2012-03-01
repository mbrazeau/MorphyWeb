/*
 *  mfyinterface.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 3/1/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */


#include "morphy.h"

mfl_handle_t* mfl_create_handle()
{
    mfl_handle_t *newhandl;
    
    return newhandl = (mfl_handle_t*)malloc(sizeof(mfl_handle_t));;
}

void mfl_destroy_handle(mfl_handle_t *mfl_handle)
{
    free(mfl_handle);
}

