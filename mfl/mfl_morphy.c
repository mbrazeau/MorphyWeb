//
//  mfl_morphy.c
//  Morphy
//
//  Created by mbrazeau on 16/05/2016.
//
//

#include "morphy.h"

void* mfl_malloc(size_t size, int memsetval, const char* fn_name)
{
    void *ret = malloc(size);
    if (!ret) {
        dbg_eprintf("unable to allocate memory for call from", fn_name);
        return NULL;
    }
    else {
        memset(ret, memsetval, size);
    }
    return ret;
}