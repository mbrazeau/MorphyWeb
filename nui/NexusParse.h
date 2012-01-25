#pragma once

#include "ncl/ncl.h"
#include "NexusReader.h"

class CNexusParse
{
public:
    CNexusParse();
    ~CNexusParse();
    bool ReadNexusFile(string *infname, string *outfname);
    void Report();
    CNexusReader* GetNexusReader() { return m_cNexusReader; }
private:
    NxsCharactersBlock  *m_cChars;
    NxsTaxaBlock        *m_cTaxa;
    NxsTreesBlock       *m_cTrees;
    NxsDataBlock        *m_cData;

    CNexusReader        *m_cNexusReader;
};
