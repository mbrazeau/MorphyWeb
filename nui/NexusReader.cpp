#include "NexusReader.h"

CNexusReader::CNexusReader(char *infname, char *outfname) : NxsReader()
{
    if (infname)
    {
        m_fIn.open(infname, ios::binary);
    }
    if (outfname)
    {
        m_fOut.open(outfname);
    }
    if (!m_fIn)
    {
        GetOutStream()<<"No input file supplied, using stdin"<<endl;
    }
}

CNexusReader::~CNexusReader()
{
    m_fIn.close();
    m_fOut.close();
}

void CNexusReader::ExecuteStarting()
{
}

void CNexusReader::ExecuteStopping()
{
}

bool CNexusReader::EnteringBlock(NxsString blockName)
{
    GetOutStream() << "Reading \"" << blockName << "\" block..." << endl;

    // Returning true means it is ok to delete any data associated with 
    // blocks of this type read in previously
    return true;    
}

void CNexusReader::SkippingBlock(NxsString blockName)
{
    GetOutStream() << "Skipping unknown block (" << blockName << ")..." << endl;
}

void CNexusReader::SkippingDisabledBlock(NxsString blockName) 
{
}

void CNexusReader::OutputComment(const NxsString &msg)
{
    GetOutStream() << msg;
}

void CNexusReader::NexusError(NxsString msg, file_pos pos, unsigned line, unsigned col)
{
    GetOutStream() << endl;
    GetOutStream() << "Error found at line " << line;
    GetOutStream() << ", column " << col;
    GetOutStream() << " (file position " << pos << "):" << endl;
    GetOutStream() << msg << endl;

    exit(0);
}

ostream &CNexusReader::GetOutStream()
{
    if (m_fOut)
    {
        return m_fOut;
    }
    return cout;
}

istream &CNexusReader::GetInStream() 
{
    if (m_fIn)
    {
        return m_fIn;
    }
    return cin;
}

void CNexusReader::statusMessage(const std::string & m) const
{
    /* 
     * cant call non-const member, but this function must be const
     * because it overrides a const function in nxsreader...
     * of course the version in nxsreader just prints to stderr
     * which we dont want...
     */
    /*
	if (alwaysReportStatusMessages || currentWarningLevel == UNCOMMON_SYNTAX_WARNING) 
    {
	    GetOutStream() << m << endl;
	}
    */
}

