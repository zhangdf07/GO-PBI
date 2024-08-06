#include "FlInputFile.h"
#include "tinyfuns.h"
#include <boost/spirit/include/classic.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>

using namespace nkchem;
using namespace boost::spirit::classic;
using namespace std;

FlInputFile::FlInputFile(void) 
{
    /* Nothing to do */ 
}

FlInputFile::FlInputFile(const string& fn) 
{
    // Read data.
    ifstream fd(fn.c_str(), ios::in);
    if(!fd.is_open()) 
    { 
        ErrorExit("Cannot open file [ %s ] to load data.", fn.c_str());
    }
    const string filecontent((istreambuf_iterator<char>(fd)), istreambuf_iterator<char>());
    fd.close();    
    string tstr;
    istringstream iss(filecontent);

    // Parse pristine.
    int tnx, tny;
    double trCC, tz;
    const rule<> RULE_pristine = (*blank_p)>>(int_p)[assign(tnx)]>>(*blank_p)
        >>(int_p)[assign(tny)]>>(*blank_p)
        >>(real_p)[assign(trCC)]>>(*blank_p)
        >>(real_p)[assign(tz)]>>(*blank_p)
        >>(+graph_p)[assign(resname)]>>(*blank_p);
    getline(iss, tstr, '\n'); tstr = triminput(tstr);
    if(!parse(tstr.c_str(), RULE_pristine).full)
    {
        ErrorExit("The pristine graphene information syntax in [ %s ] is wrong. It should be: nx ny rCC z resname.",
            fn.c_str());
    }
    nx = tnx;
    ny = tny;
    rCC = trCC;
    z = tz;
    if(nx <= 0) { ErrorExit("nx (%d) should be a positive integer.", nx); }
    if(ny <= 0) { ErrorExit("ny (%d) should be a positive integer.", ny); }
    if(rCC <= 0) { ErrorExit("rCC (%15.8f) should be a positive integer.", rCC); }

    // Parse functionalization.
    getline(iss, tstr, '\n'); tstr = triminput(tstr); ParseFunctionlization("-Oe",   2, tstr, num_edge_Oe,       edge_Oe_indices);
    getline(iss, tstr, '\n'); tstr = triminput(tstr); ParseFunctionlization("-Oe",   2, tstr, num_internal_Oe,   internal_Oe_indices);
    getline(iss, tstr, '\n'); tstr = triminput(tstr); ParseFunctionlization("-OH",   1, tstr, num_edge_OH,       edge_OH_indices);
    getline(iss, tstr, '\n'); tstr = triminput(tstr); ParseFunctionlization("-OH",   1, tstr, num_internal_OH,   internal_OH_indices);
    getline(iss, tstr, '\n'); tstr = triminput(tstr); ParseFunctionlization("-COOH", 1, tstr, num_edge_COOH,     edge_COOH_indices);
    getline(iss, tstr, '\n'); tstr = triminput(tstr); ParseFunctionlization("-COOH", 1, tstr, num_internal_COOH, internal_COOH_indices);    
}

FlInputFile::~FlInputFile(void) 
{
    /* Nothing to do */ 
}

void FlInputFile::Export(const string& fn) const
{
    FILE* fd = fopen(fn.c_str(), "w");
    if(fd == NULL)
    {
        ErrorExit("Cannot write to file [ %s ].", fn.c_str());
    }
    fprintf(fd, "%4d %4d %15.8f %15.8f %s\n", nx, ny, rCC, z, resname.c_str());
    fprintf(fd, "%4d ", num_edge_Oe); for(int i = 0; i < num_edge_Oe*2; ++i) { fprintf(fd, "%d ", edge_Oe_indices[i]); } fprintf(fd, " # Edge Oe\n");
    fprintf(fd, "%4d ", num_internal_Oe); for(int i = 0; i < num_internal_Oe*2; ++i) { fprintf(fd, "%d ", internal_Oe_indices[i]); } fprintf(fd, " # Internal Oe\n");
    fprintf(fd, "%4d ", num_edge_OH); for(int i = 0; i < num_edge_OH; ++i) { fprintf(fd, "%d ", edge_OH_indices[i]); } fprintf(fd, " # Edge OH\n");
    fprintf(fd, "%4d ", num_internal_OH); for(int i = 0; i < num_internal_OH; ++i) { fprintf(fd, "%d ", internal_OH_indices[i]); } fprintf(fd, " # Internal OH\n");
    fprintf(fd, "%4d ", num_edge_COOH); for(int i = 0; i < num_edge_COOH; ++i) { fprintf(fd, "%d ", edge_COOH_indices[i]); } fprintf(fd, " # Edge COOH\n");
    fprintf(fd, "%4d ", num_internal_COOH); for(int i = 0; i < num_internal_COOH; ++i) { fprintf(fd, "%d ", internal_COOH_indices[i]); } fprintf(fd, " # Internal COOH\n");    
    fclose(fd);
}

void FlInputFile::ParseFunctionlization(const string& tag, int dup, const string& str, int& num, vector<int>& indices) const
{
    vector<string> buffer;
    string tstr;
    istringstream reader(str);
    while(reader >> tstr) 
    {
        if(tstr.length()) 
        {
            buffer.push_back(tstr);
        }
    }

    if(buffer.size() == 0) { ErrorExit("The number of functionlization groups (%s) is needed.", tag.c_str()); }
    int tn;
    const rule<> RULE_index = (int_p)[assign(tn)];
    if(!parse(buffer[0].c_str(), RULE_index).full || tn < 0)
    {
        ErrorExit("A non-negative integer is needed for the number of functionalization groups (%s) but not %s.",
            tag.c_str(), buffer[0].c_str());
    }
    num = tn;
    indices.clear();
    if(num > 0)
    {
        if(buffer.size() == 1) { return; }
        if(buffer.size() != num*dup+1) { ErrorExit("Excatly %d functionalization indices (%s) are needed.", num*dup, tag.c_str()); }
        for(int i = 1; i <= num*dup; ++i)
        {
            if(!parse(buffer[i].c_str(), RULE_index).full || tn == 0)
            {
                ErrorExit("A nonzero integer is needed for the functionalization index (%s) but not %s.", tag.c_str(),  buffer[i].c_str());
            }
            indices.push_back(tn);
        }
    }
}

string FlInputFile::triminput(const string& str) const
{
    // Delete the comment.
    string tstr = str.substr(0, str.find('#'));
    // Delete the spaces.
    tstr = tstr.substr(0, tstr.find_last_not_of(" \n\r\t")+1);
    if(!tstr.empty())
    {
    tstr = tstr.substr(tstr.find_first_not_of(" \n\r\t"));
    }
    return tstr;
}