#pragma once
#include <iterator>
class CNexusUserInterface;

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
    CNexusMenuBase(const char * strCommand, const char * strHelpText)
    {
        InitMenu(strCommand, strHelpText);
    }

    CNexusMenuBase(const char * strCommand, const char * strHelpText, map<const char*, int, ltstr> assignments)
    {
        InitMenu(strCommand, strHelpText);
        m_mapAssignments = assignments;
    }

    CNexusMenuBase(const char * strCommand, const char * strHelpText, vector<int> assignments)
    {
        vector<int>::iterator it;
        InitMenu(strCommand, strHelpText);
        for (it = assignments.begin(); it < assignments.end(); it++)
        {
            m_intAssignments.push_back(*it);
        }
    }

    void InitMenu(const char * strCommand, const char * strHelpText)
    {
        if (strCommand)
        {
            m_strCommand      = strCommand;
            m_strCommandLower = strCommand;
            transform(m_strCommandLower.begin(), m_strCommandLower.end(), m_strCommandLower.begin(), ::tolower);
        }
        m_strHelpText = SplitToMaxLen(strHelpText, MAX_HELP_WIDTH);
    }
    virtual ~CNexusMenuBase()
    {
    }

    void GetMenuOutput(string *output)
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
    
    void GetValidParams(string *params)
    {
        map<const char*, int, ltstr>::iterator itm;
        ostringstream strStream;
        for (itm = m_mapAssignments.begin(); itm != m_mapAssignments.end(); ++itm)
        {
            strStream<<(*itm).first<<" ";
        }
        *params = strStream.str();
    }

    bool IsSelection(string strInput)
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

    ENexusMenuCommandStatus RunCommand(CNexusUserInterface *pNexusUserInterface, string value)
    {
        ENexusMenuCommandStatus eRet;
        bool bStatus;
        int nMappedVal;
        eRet = ValidateIntInput(value, &nMappedVal);

        if (eRet == ENXS_MCS_OK)
        {
            eRet = ValidateMapInput(value, &nMappedVal);
        }

        if (eRet == ENXS_MCS_OK)
        {
            bStatus = MenuFunction(pNexusUserInterface, &value, nMappedVal);
            if (bStatus == false)
            {
                eRet = ENXS_MCS_COMMAND_FAIL;
            }
        }
        return eRet;
    }

    virtual bool MenuFunction(CNexusUserInterface *pNexusUserInterface, string *value, int nMappedVal) = 0;
protected:
    
    string FormatHelpText()
    {
        string ret;
        vector<string>::iterator it;
        for (it = m_strHelpText.begin(); it < m_strHelpText.end(); it++)
        {
            ret += *it;
        }
        return ret;
    }

    vector<string> SplitToMaxLen(string text, size_t max)
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
    
    ENexusMenuCommandStatus ValidateMapInput(string value, int *nMappedVal)
    {
        ENexusMenuCommandStatus eRet = ENXS_MCS_OK;
        if (m_mapAssignments.size() > 0)
        {
            map<const char*, int>::const_iterator itr;
            itr = m_mapAssignments.find(value.c_str());
            if (itr == m_mapAssignments.end())
            {
                eRet = ENXS_MCS_INVALID_PARAM;
            }
            else
            {
                *nMappedVal = (*itr).second;
            }
        }
        return eRet;
    }

    ENexusMenuCommandStatus ValidateIntInput(string value, int *nMappedVal)
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
    string m_strCommand;
    string m_strCommandLower;
    vector<string> m_strHelpText;
    vector<int> m_intAssignments;
    map<const char*, int, ltstr> m_mapAssignments;
};

class CNexusMenuData
{
public:
    CNexusMenuData(string strMenuTitle)
    {
        bMenuShown = false;
        m_strMenuTitle = strMenuTitle;
    }
    ~CNexusMenuData()
    {
        vector<CNexusMenuBase*>::iterator it;
        CNexusMenuBase* pMenuItem;
       
        for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
        {
            pMenuItem = *it;
            if (pMenuItem)
            {
                delete(pMenuItem);
            }
        }
        m_vMenu.clear();
    }

    void AddMenuItem(CNexusMenuBase* pMenuItem)
    {
        m_vMenu.push_back(pMenuItem);
    }
    
    bool RunSelection(string strInput, CNexusUserInterface *pNexusUserInterface)
    {
        vector<CNexusMenuBase*> pMenuItems;
        string command;
        string value;
        SplitInput(strInput, &command, &value);
        pMenuItems = GetMenuSelection(command);
        if (pMenuItems.size() == 1)
        {
            bool bRet = true;
            cout<<endl;
            try
            {
                switch (pMenuItems[0]->RunCommand(pNexusUserInterface, value))
                {
                case ENXS_MCS_COMMAND_FAIL:
                    bRet = false;
                    break;

                case ENXS_MCS_INVALID_PARAM:
                    PrintError(strInput, pMenuItems);
                    break;

                default:
                    break;
                }
            }
            catch (const char *e)
            {
                cout<<"NUI Error: "<<e<<endl;
            }
            catch (mfl_exception e)
            {
                cout<<"MFL Error: "<<e.what()<<endl;
            }
            cout<<endl;
            return bRet;
        }
        PrintError(strInput, pMenuItems);

        return true;
    }

    bool Help(bool bForceShowMenu)
    {
        if ((bMenuShown == false) || (bForceShowMenu == true))
        {
            vector<CNexusMenuBase*>::iterator it;
            CNexusMenuBase* pMenu;
            string output;
           
            for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
            {
                pMenu = *it;
                if (pMenu)
                {
                    pMenu->GetMenuOutput(&output);
                    cout<<output<<endl;
                }
            }
            cout<<endl;
            bMenuShown = true;
        }
        return true;
    }

    string GetPrompt()
    {
        return m_strMenuTitle;
    }

protected:

    std::string &ltrim(std::string &s) 
    {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }

    std::string &rtrim(std::string &s) 
    {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    std::string &trim(std::string &s) 
    {
        return ltrim(rtrim(s));
    }

    void SplitInput(string strInput, string *command, string *value)
    {
        *command = strInput;
        size_t index = strInput.find('=');
        if (index != string::npos)
        {
            *command = strInput.substr(0, index);
            *value = strInput.substr(index+1);
            trim(*command);
            trim(*value);
        }
    }

    void PrintError(string strInput, vector<CNexusMenuBase*> pMenuItems)
    {
        vector<CNexusMenuBase*>::iterator it;
        cout<<" Unknown command: '"<<strInput<<"'"<<endl<<endl;
        CNexusMenuBase* pMenuItem;
        string params;
        string output;
        for (it = pMenuItems.begin(); it < pMenuItems.end(); it++)
        {
            pMenuItem = *it;
            pMenuItem->GetMenuOutput(&output);
            pMenuItem->GetValidParams(&params);
            cout<<output<<endl;
            if (params.length() > 0)
            {
                cout<<setw(COMMAND_WIDTH - 2)<<" "<<"Possible values:"<<endl;
                cout<<setw(COMMAND_WIDTH)<<" "<<params<<endl;
            }
        }
        cout<<endl;
    }

    vector<CNexusMenuBase*> GetMenuSelection(string strInput)
    {
        vector<CNexusMenuBase*>::iterator it;
        CNexusMenuBase* pMenuItem;
        vector<CNexusMenuBase*> ret;
        for (it = m_vMenu.begin(); it < m_vMenu.end(); it++)
        {
            pMenuItem = *it;
            if ((pMenuItem) && (pMenuItem->IsSelection(strInput)))
            {
                ret.push_back(pMenuItem);
            }
        }
        return ret;
    }

    bool bMenuShown;
    string m_strMenuTitle;
    vector <CNexusMenuBase*> m_vMenu;
};
