#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "NexusMenuBase.h"

CNexusMenuBase::CNexusMenuBase(const char * strCommand, const char * strHelpText)
{
    InitMenu(strCommand, strHelpText);
}

CNexusMenuBase::CNexusMenuBase(const char * strCommand, const char * strHelpText, map<const char*, int, ltstr> assignments)
{
    InitMenu(strCommand, strHelpText);
    m_mapAssignments = assignments;
}

CNexusMenuBase::CNexusMenuBase(const char * strCommand, const char * strHelpText, vector<int> assignments)
{
    vector<int>::iterator it;
    InitMenu(strCommand, strHelpText);
    for (it = assignments.begin(); it < assignments.end(); it++)
    {
        m_intAssignments.push_back(*it);
    }
}

CNexusMenuBase::~CNexusMenuBase()
{
}

void CNexusMenuBase::InitMenu(const char * strCommand, const char * strHelpText)
{
    if (strCommand)
    {
        m_strCommand      = strCommand;
        m_strCommandLower = strCommand;
        transform(m_strCommandLower.begin(), m_strCommandLower.end(), m_strCommandLower.begin(), ::tolower);
    }
    m_strHelpText = SplitToMaxLen(strHelpText, MAX_HELP_WIDTH);
}

void CNexusMenuBase::GetMenuOutput(string *output)
{
    ostringstream strStream;
    if (m_strCommand.length() > 0)
    {
        strStream<<setw(COMMAND_WIDTH)<<m_strCommand<<") ";
        strStream<<FormatHelpText();
    }
    else
    {
        strStream<<"\n=== ";
        strStream<<FormatHelpText();
    }
    *output = strStream.str();
}

void CNexusMenuBase::GetValidParams(string *params)
{
    map<const char*, int, ltstr>::iterator itm;
    ostringstream strStream;
    for (itm = m_mapAssignments.begin(); itm != m_mapAssignments.end(); ++itm)
    {
        strStream<<(*itm).first<<" ";
    }
    *params = strStream.str();
}

bool CNexusMenuBase::IsSelection(string strInput)
{
    int index;
    if (strInput.length() > 0)
    {
        transform(strInput.begin(), strInput.end(), strInput.begin(), ::tolower);
        index = m_strCommandLower.find(strInput);
        return (index == 0);
    }
    return false;
}

void CNexusMenuBase::PrintCurrentCommand(string fullAssignment)
{
    cout<<"Running: "<<m_strCommand;
    if (fullAssignment.length() > 0)
    {
        cout<<" = "<<fullAssignment;
    }
    cout<<endl;
}

ENexusMenuCommandStatus CNexusMenuBase::RunCommand(CNexusUserInterface *pNexusUserInterface, string value)
{
    ENexusMenuCommandStatus eRet;
    bool bStatus;
    int nMappedVal;
    string fullAssignment = value;
    eRet = ValidateIntInput(value, &nMappedVal);

    if (eRet == ENXS_MCS_OK)
    {
        eRet = ValidateMapInput(value, &fullAssignment, &nMappedVal);
    }

    if (eRet == ENXS_MCS_OK)
    {
        PrintCurrentCommand(fullAssignment);
        bStatus = MenuFunction(pNexusUserInterface, &value, nMappedVal);
        if (bStatus == false)
        {
            eRet = ENXS_MCS_COMMAND_FAIL;
        }
    }
    return eRet;
}

string CNexusMenuBase::FormatHelpText()
{
    string ret;
    vector<string>::iterator it;
    for (it = m_strHelpText.begin(); it < m_strHelpText.end(); it++)
    {
        ret += *it;
    }
    return ret;
}

vector<string> CNexusMenuBase::SplitToMaxLen(string text, size_t max)
{
    vector<string> ret;
    vector<string> tokens;
    istringstream iss(text);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));
    vector<string>::iterator it;
    ostringstream line;
    for (it = tokens.begin(); it < tokens.end(); it++)
    {
        line<<*it<<" ";
        if (line.str().size() >= max)
        {
            ret.push_back(line.str());
            line.str("");
            line<<"\n"<<setw(COMMAND_WIDTH + 2)<<" ";
        }
    }
    ret.push_back(line.str());

    return ret;
}

ENexusMenuCommandStatus CNexusMenuBase::ValidateMapInput(string value, string *fullAssignment, int *nMappedVal)
{
    ENexusMenuCommandStatus eRet = ENXS_MCS_OK;
    if (m_mapAssignments.size() > 0)
    {
        int ret = FindValueInAssignmentMap(value, fullAssignment);
        if (ret >= 0)
        {
            *nMappedVal = ret;
        }
        else
        {
            eRet = ENXS_MCS_INVALID_PARAM;
        }
    }
    return eRet;
}

ENexusMenuCommandStatus CNexusMenuBase::ValidateIntInput(string value, int *nMappedVal)
{
    ENexusMenuCommandStatus eRet = ENXS_MCS_OK;
    if (m_intAssignments.size() > 0)
    {
        int v;
        istringstream(value)>>v;
        eRet = ENXS_MCS_INVALID_PARAM;
        if ((v >= m_intAssignments[0]) && (v <= m_intAssignments[1]))
        {
            *nMappedVal = v;
            eRet = ENXS_MCS_OK;
        }
    }
    return eRet;
}

int CNexusMenuBase::FindValueInAssignmentMap(string value, string *fullAssignment)
{
    int nFound = 0;
    map<const char*, int>::const_iterator iFound;
    map<const char*, int>::const_iterator itr;
    string possibleValue;
    int index;
    transform(value.begin(), value.end(), value.begin(), ::tolower);
    for (itr = m_mapAssignments.begin(); itr != m_mapAssignments.end(); itr++)
    {
        possibleValue = (*itr).first;
        transform(possibleValue.begin(), possibleValue.end(), possibleValue.begin(), ::tolower);
        index = possibleValue.find(value);
        if (index == 0)
        {
            nFound++;
            iFound = itr;
        }
    }
    if (nFound == 1)
    {
        *fullAssignment = (*iFound).first;
        return (*iFound).second;
    }
    return -1;
}
