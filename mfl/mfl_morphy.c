//
//  mfl_morphy.c
//  Morphy
//
//  Created by mbrazeau on 16/05/2016.
//
//

#include "morphy.h"

void* __MFL_MALLOC__(size_t size, int memsetval, const char* fn_name)
{
    void *ret = malloc(size);
    if (!ret) {
        dbg_printf("Error in __MFL_MALLOC__(): unable to allocate memory for call from %s()\n\n", fn_name);
        return NULL;
    }
    else {
        memset(ret, memsetval, size);
    }
    return ret;
}