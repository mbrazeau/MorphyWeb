#pragma once
#include <iterator>
#include <string.h>
#include <map>
#include <vector>
class CNexusUserInterface;

using namespace std;

enum ENexusMenuCommandStatus
{
    ENXS_MCS_OK,
    ENXS_MCS_INVALID_PARAM,
    ENXS_MCS_COMMAND_FAIL,
};

struct ltstr
{
    bool operator()(const char* s1, const char* s2) const
    {
        return strcmp(s1, s2) < 0;
    }
};

#define COMMAND_WIDTH 12
#define MAX_HELP_WIDTH 65

class CNexusMenuBase
{
public:
    CNexusMenuBase(const char * strCommand, const char * strHelpText);
    CNexusMenuBase(const char * strCommand, const char * strHelpText, map<const char*, int, ltstr> assignments);
    CNexusMenuBase(const char * strCommand, const char * strHelpText, vector<int> assignments);
    void InitMenu(const char * strCommand, const char * strHelpText);
    virtual ~CNexusMenuBase();
    void GetMenuOutput(string *output);
    void GetValidParams(string *params);
    bool IsSelection(string strInput);
    ENexusMenuCommandStatus RunCommand(CNexusUserInterface *pNexusUserInterface, string value);
    virtual bool MenuFunction(CNexusUserInterface *pNexusUserInterface, string *value, int nMappedVal) = 0;
protected:
    
    string FormatHelpText();
    vector<string> SplitToMaxLen(string text, size_t max);
    ENexusMenuCommandStatus ValidateMapInput(string value, string *fullAssignment, int *nMappedVal);
    ENexusMenuCommandStatus ValidateIntInput(string value, int *nMappedVal);
    int FindValueInAssignmentMap(string value, string *fullAssignment);
    void PrintCurrentCommand(string fullAssignment);

    string m_strCommand;
    string m_strCommandLower;
    vector<string> m_strHelpText;
    vector<int> m_intAssignments;
    map<const char*, int, ltstr> m_mapAssignments;
};

