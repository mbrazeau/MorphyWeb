#include "NexusReader.h"

CNexusReader::CNexusReader(char *infname, char *outfname) : NxsReader()
{
    m_fIn.open(infname, ios::binary);
    m_fOut.open(outfname);
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
    m_fOut << "Reading \"" << blockName << "\" block..." << endl;

    // Returning true means it is ok to delete any data associated with 
    // blocks of this type read in previously
    //
    return true;    
}

void CNexusReader::SkippingBlock(NxsString blockName)
{
    m_fOut << "Skipping unknown block (" << blockName << ")..." << endl;
}

void CNexusReader::SkippingDisabledBlock(NxsString blockName) 
{
}

void CNexusReader::OutputComment(const NxsString &msg)
{
    m_fOut << msg;
}

void CNexusReader::NexusError(NxsString msg, file_pos pos, unsigned line, unsigned col)
{
    m_fOut << endl;
    m_fOut << "Error found at line " << line;
    m_fOut << ", column " << col;
    m_fOut << " (file position " << pos << "):" << endl;
    m_fOut << msg << endl;

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

void CNexusReader::statusMessage(const std::string & m) const
{
    /* 
     * cant call non-const member, but this function must be const
     * because it overrides a const function in nxsreader...
     * of course the version in nxsreader just prints to stderr
     * which we dont want...
     */
    /*
    ostream &oStream = GetOutStream();
    
	if (alwaysReportStatusMessages || currentWarningLevel == UNCOMMON_SYNTAX_WARNING) 
    {
	    oStream << m << std::endl;
	}
    */
}

