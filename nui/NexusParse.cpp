#include "NexusParse.h"

CNexusParse::CNexusParse()
{
    m_cChars = NULL;
    m_cTaxa = NULL;
    m_cTrees = NULL;
    m_cData = NULL;
    m_cNexusReader = NULL;
}

CNexusParse::~CNexusParse()
{
    if (m_cTaxa)
    {
        delete m_cTaxa;
        m_cTaxa = NULL;
    }
    if (m_cChars)
    {
        delete m_cChars;
        m_cChars = NULL;
    }
    if (m_cTrees)
    {
        delete m_cTrees;
        m_cTrees = NULL;
    }
    if (m_cData)
    {
        delete m_cData;
        m_cData = NULL;
    }
    if (m_cNexusReader)
    {
        delete m_cNexusReader;
        m_cNexusReader = NULL;
    }
}

bool CNexusParse::ReadNexusFile(string *infname, string *outfname)
{
    bool bRet = false;

    /*
     * if a new fails, the the destructor will take care of freeing
     * any memory that happened to work...
     */
    m_cNexusReader = new CNexusReader(infname, outfname);
    m_cTaxa  = new NxsTaxaBlock();
    m_cChars = new NxsCharactersBlock(0, 0);
    m_cData  = new NxsDataBlock(m_cTaxa, 0);
    
    if (m_cNexusReader && m_cTaxa && m_cChars && m_cData)
    {
        m_cTrees = new NxsTreesBlock(m_cTaxa);
        if (m_cTrees)
        {
            m_cNexusReader->Add(m_cTaxa);
            m_cNexusReader->Add(m_cChars);
            m_cNexusReader->Add(m_cTrees);
            m_cNexusReader->Add(m_cData);
           
            istream &iStream = m_cNexusReader->GetInStream();
            if (iStream)
            {
                CNexusToken token(iStream, m_cNexusReader->GetOutStream());
                m_cNexusReader->Execute(token);
                bRet = true;
            }
        }
    }
    if (HasDataBlock())
    {
        m_cData->TransferTo(*m_cChars);
    }
    return bRet;
}

bool CNexusParse::HasDataBlock()
{
    BlockReaderList b1 = m_cNexusReader->GetUsedBlocksInOrder();
    for (std::list<NxsBlock *>::iterator i = b1.begin(); i != b1.end(); ++i) 
    {
        NxsBlock * b = (*i);
        if (b->GetID().compare("DATA") == 0)
        {
            return true;
        }
    }
    return false;
}

void CNexusParse::Report()
{
    if (m_cNexusReader && m_cTaxa && m_cChars && m_cData && m_cTrees)
    {
        ostream &oStream = m_cNexusReader->GetOutStream();
        m_cTaxa->Report (oStream);
        m_cChars->Report(oStream);
        m_cTrees->Report(oStream);
        m_cData->Report (oStream);
    }
}

