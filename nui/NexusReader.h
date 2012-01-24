#pragma once

#include "ncl/ncl.h"

class CNexusReader : public NxsReader
{
public:
    CNexusReader(char *infname, char *outfname);
    ~CNexusReader();
    void ExecuteStarting();
    void ExecuteStopping();
    bool EnteringBlock(NxsString blockName);
    void SkippingBlock(NxsString blockName);
    void SkippingDisabledBlock(NxsString blockName);
    void OutputComment(const NxsString &msg);
    void NexusError(NxsString msg, file_pos pos, unsigned line, unsigned col);
    ostream &GetOutStream();
    istream &GetInStream();
    void statusMessage(const std::string & m) const;

private:
    ofstream m_fOut;
    ifstream m_fIn;
};

class CNexusToken : public NxsToken
{
public:

    CNexusToken(istream &is, ostream &os) : NxsToken(is), m_fOut(os){}

    void OutputComment(const NxsString &msg)
    {
        m_fOut << msg << endl;
    }

private:
    ostream &m_fOut;
};
