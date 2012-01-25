#include "NexusParse.h"

CNexusParse::CNexusParse()
{
}

CNexusParse::~CNexusParse()
{
    if (m_cTaxa)
    {
        delete m_cTaxa;
    }
    if (m_cChars)
    {
        delete m_cChars;
    }
    if (m_cTrees)
    {
        delete m_cTrees;
    }
    if (m_cData)
    {
        delete m_cData;
    }
    if (m_cNexus)
    {
        delete m_cNexus;
    }
}

bool CNexusParse::ReadNexusFile(string *infname, string *outfname)
{
    bool bRet = false;

    /*
     * if a new fails, the the destructor will take care of freeing
     * any memory that happened to work...
     */
    m_cNexus = new CNexusReader(infname, outfname);
    m_cTaxa  = new NxsTaxaBlock();
    m_cChars = new NxsCharactersBlock(0, 0);
    m_cData  = new NxsDataBlock(0, 0);
    if (m_cNexus && m_cTaxa && m_cChars && m_cData)
    {
        m_cTrees = new NxsTreesBlock(m_cTaxa);
        if (m_cTrees)
        {
            m_cNexus->Add(m_cTaxa);
            m_cNexus->Add(m_cChars);
            m_cNexus->Add(m_cTrees);
            m_cNexus->Add(m_cData);
           
            istream &iStream = m_cNexus->GetInStream();
            /* This needs to be improved... */
            if (iStream != cin)
            {
                CNexusToken token(iStream, m_cNexus->GetOutStream());
                m_cNexus->Execute(token);
                bRet = true;
            }
        }
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

