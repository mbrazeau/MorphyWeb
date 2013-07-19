/*
 *  cdata.c
 *  Morphy
 *
 *  Created by Martin Brazeau on 2/25/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "morphy.h"

chardata *mfl_new_chardata(void)
{
    return (chardata*)malloc(sizeof(chardata));
}

void mfl_clear_chardata(chardata *cd)
{
    if (cd->cd_nogaps) {
        free(cd->cd_nogaps);
    }
    if (cd->cd_wgaps) {
        free(cd->cd_wgaps);
    }
}

void mfl_delete_chardata(chardata *cd)
{
    mfl_clear_chardata(cd);
    free(cd);
}
        