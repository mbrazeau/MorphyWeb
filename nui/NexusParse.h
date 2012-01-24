#pragma once

#include "ncl/ncl.h"
#include "NexusReader.h"

class CNexusParse
{
public:
    CNexusParse(string *infname, string *outfname);
    ~CNexusParse();
    bool ReadNexusFile();
    void Report();

private:
    NxsCharactersBlock  *m_cChars;
    NxsTaxaBlock        *m_cTaxa;
    NxsTreesBlock       *m_cTrees;
    NxsDataBlock        *m_cData;
    CNexusReader        *m_cNexus;
};
