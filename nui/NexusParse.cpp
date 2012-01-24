#include "NexusParse.h"

CNexusParse::CNexusParse(string *infname, string *outfname)
{
    m_cNexus = new CNexusReader(infname, outfname);
    m_cTaxa  = new NxsTaxaBlock();
    m_cChars = new NxsCharactersBlock(0, 0);
    m_cTrees = new NxsTreesBlock(m_cTaxa);
    m_cData  = new NxsDataBlock(0, 0);

    m_cNexus->Add(m_cTaxa);
    m_cNexus->Add(m_cChars);
    m_cNexus->Add(m_cTrees);
    m_cNexus->Add(m_cData);
}

CNexusParse::~CNexusParse()
{
    delete m_cTaxa;
    delete m_cChars;
    delete m_cTrees;
    delete m_cData;
    delete m_cNexus;
}

bool CNexusParse::ReadNexusFile()
{
    bool bRet = false;
    istream &iStream = m_cNexus->GetInStream();
    /* This needs to be improved... */
    if (iStream != cin)
    {
        CNexusToken token(iStream, m_cNexus->GetOutStream());
        m_cNexus->Execute(token);
        bRet = true;
    }
    return bRet;
}

void CNexusParse::Report()
{
    ostream &oStream = m_cNexus->GetOutStream();
    m_cTaxa->Report (oStream);
    m_cChars->Report(oStream);
    m_cTrees->Report(oStream);
    m_cData->Report (oStream);
}

