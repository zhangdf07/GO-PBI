#include "GrInputFile.h"
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

GrInputFile::GrInputFile(const string& fn) 
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
    getline(iss, tstr, '\n'); tstr = triminput(tstr); ParseFunctionlization("-OH",   1, tstr, num_OH,   num_OH_up,   style_OH,   OH_indices);
    getline(iss, tstr, '\n'); tstr = triminput(tstr); ParseFunctionlization("-O-",   2, tstr, num_Oe,   num_Oe_up,   style_Oe,   Oe_indices);
    getline(iss, tstr, '\n'); tstr = triminput(tstr); ParseFunctionlization("-H",    1, tstr, num_H,    num_H_up,    style_H,    H_indices);
    getline(iss, tstr, '\n'); tstr = triminput(tstr); ParseFunctionlization("-COOH", 1, tstr, num_COOH, num_COOH_up, style_COOH, COOH_indices);
}

GrInputFile::~GrInputFile(void) 
{
    /* Nothing to do */ 
}

void GrInputFile::Export(const string& fn) const
{
    FILE* fd = fopen(fn.c_str(), "w");
    if(fd == NULL)
    {
        ErrorExit("Cannot write to file [ %s ].", fn.c_str());
    }
    fprintf(fd, "%4d %4d %15.8f %15.8f %s\n", nx, ny, rCC, z, resname.c_str());
    fprintf(fd, "%4d ", num_OH); for(int i = 0; i < num_OH; ++i) { fprintf(fd, "%d ", OH_indices[i]); } fprintf(fd, " # OH\n");
    fprintf(fd, "%4d ", num_Oe); for(int i = 0; i < num_Oe*2; ++i) { fprintf(fd, "%d ", Oe_indices[i]); } fprintf(fd, " # Oe\n");
    fprintf(fd, "%4d ", num_H);  for(int i = 0; i < num_H; ++i)  { fprintf(fd, "%d ", H_indices[i]); }  fprintf(fd, " # H\n");
    fprintf(fd, "%4d ", num_COOH);  for(int i = 0; i < num_COOH; ++i)  { fprintf(fd, "%d ", COOH_indices[i]); }  fprintf(fd, " # COOH\n");

    fclose(fd);
}

void GrInputFile::ParseFunctionlization(const string& tag, int dup, const string& str, int& num, int& num_up, Style& sty, vector<int>& indices) const
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
        if(buffer.size() < 2) { ErrorExit("At least one character/integer is needed for the functionalization (%s).", tag.c_str()); }
        if(buffer[1] == RandomChar)
        {
            sty = Random;
            if(buffer.size() < 3) { ErrorExit("An integer is needed for the number of upper functionalizatin groups (%s).", tag.c_str()); }
            if(!parse(buffer[2].c_str(), RULE_index).full || tn > num  || tn < 0)
            {
                ErrorExit("An integer between 0 and %d is needed for the number of upper functionalization groups (%s) but not %s.",
                    num, tag.c_str(), buffer[2].c_str());            
            }
            num_up = tn;           
        }
        else
        {
            if(buffer[1] == EvenlyChar)
            {
                sty = Evenly;
            }
            else
            {
                sty = Assigned;
                if(buffer.size() != num*dup+1) { ErrorExit("Excatly %d functionalization indices (%s) are needed, or give \"R\" or \"E\".", num*dup, tag.c_str()); }
                for(int i = 1; i <= num*dup; ++i)
                {
                    if(!parse(buffer[i].c_str(), RULE_index).full || tn == 0)
                    {
                        ErrorExit("A nonzero integer is needed for the functionalization index (%s) but not %s.", buffer[i].c_str(), tag.c_str());
                    }
                    indices.push_back(tn);
                }
            }
        }
    }
}

string GrInputFile::triminput(const string& str) const
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