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
    void statusMessage(const std::string & m) const;

    ifstream m_fIn;
   
private:
    ofstream m_fOut;
};

class MyToken : public NxsToken
{
public:

    MyToken(istream &is, ostream &os) : NxsToken(is), out(os){}

    void OutputComment(const NxsString &msg)
    {
        out << msg << endl;
    }

private:
    ostream &out;
};
