#pragma once

class myofstream : public ofstream
{
public:
    void open (const char * filename, ios_base::openmode mode = ios_base::out) 
    {
        ofstream::open(filename, mode);
        if (is_open())
        {
            m_strFileName = filename;
        }
    }
    void close ()
    {
        ofstream::close();
        m_strFileName.clear();
    }
    string &GetFileName()
    {
        if (!is_open())
        {
            m_strFileName = "File not open";
        }
        return m_strFileName;
    }
private:
    string m_strFileName;
};

class myifstream : public ifstream
{
public:
    void open (const char * filename, ios_base::openmode mode = ios_base::in)
    {
        ifstream::open(filename, mode);
        if (is_open())
        {
            m_strFileName = filename;
        }
    }
    void close ()
    {
        ifstream::close();
        m_strFileName.clear();
    }
    string &GetFileName()
    {
        if (!is_open())
        {
            m_strFileName = "File not open";
        }
        return m_strFileName;
    }

private:
    string m_strFileName;
};

