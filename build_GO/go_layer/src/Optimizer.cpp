#include "Optimizer.h"
#include "GrFlake.h"
#include "FlInputFile.h"
#include "FlBuilder.h"
#include "tinyfuns.h"
#include <boost/spirit/include/classic.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/stat.h>
#endif

using namespace nkchem;
using namespace boost::spirit::classic;
using namespace std;

Optimizer::Optimizer(void) 
{
    /* Nothing to do */ 
}

Optimizer::~Optimizer(void) 
{
    /* Nothing to do */ 
}

void Optimizer::Do(const string& fn)
{
    printf("==========================================\n");
    printf("|        Graphene Flake Optimizer        |\n");
    printf("==========================================\n");
    printf("\n");
    printf("Input file: %s\n", fn.c_str());
    ParseInput(fn);
    printf("\n");

    printf("Graphene flake:  %d x %d, rCC = %.2f\n", nx, ny, rCC);
    printf("Flake name:      %s\n", resname.c_str());
    GrFlake flake0(nx, ny, rCC, z, GrFlake::Parallel);

    printf("Derocation:\n");
    printf(" Number of Oe:   %d\n", num_Oe);
    printf(" Number of OH:   %d\n", num_OH);
    printf(" Number of COOH: %d\n", num_COOH);

    printf("Result folder:           %s\n", result_folder.c_str());
    if(!makedir(result_folder)) { ErrorExit("Cannot create the folder \"%s\".", result_folder.c_str()); }
    printf("Maximum steps of search: %d\n", max_step);
    printf("xTB command:             %s\n", cmd.c_str());

    char buffer[4098];
    const char* Separator = "-------------------------------------------------------------\n";
    printf(Separator);
    printf("%6s   %s\n", "Step", "Energy");
    printf(Separator);
    for(int i = 0; i < max_step; ++i)
    {
        // Initialization.
        GrFlake flake(flake0);
        FlInputFile inp;
        inp.nx = nx; inp.ny = ny; inp.rCC = rCC; inp.z = z; inp.resname = resname;
        inp.num_edge_Oe = 0; inp.num_internal_Oe = num_Oe;
        inp.num_edge_OH = 0; inp.num_internal_OH = num_OH;
        inp.num_edge_COOH = num_COOH; inp.num_internal_COOH = 0;

        FlBuilder::Build(inp, flake);
        // Export.
        sprintf(buffer, "%s/%s-%d", result_folder.c_str(), inp.resname.c_str(), i);
        flake.Export(buffer);
        sprintf(buffer, "%s/%s-%d.inp", result_folder.c_str(), inp.resname.c_str(), i);
        inp.Export(buffer);

        // Calculate and Rank.
        sprintf(buffer, "%s %s/%s-%d.xyz -o > %s/%s-%d.xtbout 2>/dev/null ; cp xtbopt.xyz %s/%s-%d-opt.xyz", cmd.c_str(), result_folder.c_str(), inp.resname.c_str(), i, 
            result_folder.c_str(), inp.resname.c_str(), i, result_folder.c_str(), inp.resname.c_str(), i);
        system(buffer);
        FILE* fp = popen("awk 'NR==2{print $2}' xtbopt.xyz", "r");
        fgets(buffer, 4098, fp);
        fclose(fp);       
        printf("%6d   %5s", i, buffer);
    }
    printf(Separator);
}

void Optimizer::ParseInput(const string& fn) 
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
    
    int tnx, tny, tnz;
    double trCC, tz;
    string te;

    // Parse pristine.
    const rule<> RULE_pristine = (*blank_p)>>(int_p)[assign(tnx)]>>(*blank_p)
        >>(int_p)[assign(tny)]>>(*blank_p)
        >>(real_p)[assign(trCC)]>>(*blank_p)
        >>(real_p)[assign(tz)]>>(*blank_p)
        >>(+graph_p)[assign(resname)]>>(*blank_p);
    getline(iss, tstr, '\n'); tstr = triminput(tstr);
    if(!parse(tstr.c_str(), RULE_pristine).full)
    {
        ErrorExit("The pristine graphene information syntax in [ %s ] (%s) is wrong. It should be: nx ny rCC z resname.",
            fn.c_str(), tstr.c_str());
    }
    nx = tnx;
    ny = tny;
    rCC = trCC;
    z = tz;
    if(nx <= 0) { ErrorExit("nx (%d) should be a positive integer.", nx); }
    if(ny <= 0) { ErrorExit("ny (%d) should be a positive integer.", ny); }
    if(rCC <= 0) { ErrorExit("rCC (%15.8f) should be a positive integer.", rCC); }

    // Parse decorations.
    const rule<> RULE_decorations = (*blank_p)>>(int_p)[assign(tnx)]>>(*blank_p)
        >>(int_p)[assign(tny)]>>(*blank_p)
        >>(int_p)[assign(tnz)]>>(*blank_p);
    getline(iss, tstr, '\n'); tstr = triminput(tstr);
    if(!parse(tstr.c_str(), RULE_decorations).full || tnx < 0 || tny < 0 || tnz < 0)
    {
        ErrorExit("The decoration information syntax in [ %s ] (%s) is wrong. It should be three non-negative integers.",
            fn.c_str(), tstr.c_str());
    }
    num_Oe = tnx;
    num_OH = tny;
    num_COOH = tnz;

    // Parse working stuff.
    getline(iss, tstr, '\n'); tstr = triminput(tstr);
    if(!parse(tstr.c_str(), (*blank_p)>>(+graph_p)[assign(result_folder)]>>(*blank_p)).full)
    {
        ErrorExit("The result folder syntax in [ %s ] (%s) is wrong. It should be a string.",
            fn.c_str(), tstr.c_str());
    }

    getline(iss, tstr, '\n'); tstr = triminput(tstr);
    if(!parse(tstr.c_str(), (*blank_p)>>(int_p)[assign(max_step)]>>(*blank_p)).full || max_step <= 0)
    {
        ErrorExit("The max step syntax in [ %s ] (%s) is wrong. It should be a non-negative integer.",
            fn.c_str(), tstr.c_str());
    }

    // Parse cmd.
    getline(iss, tstr, '\n'); tstr = triminput(tstr);
    if(!parse(tstr.c_str(), (*blank_p)>>(+graph_p)[assign(cmd)]>>(*blank_p)).full)
    {
        ErrorExit("The command syntax in [ %s ] (%s) is wrong. It should be a string.",
            fn.c_str(), tstr.c_str());
    }
}

string Optimizer::triminput(const string& str) const
{
    // Delete the comment.
    string tstr = str.substr(0, str.find('#'));
    // Delete the spaces.
    tstr = tstr.substr(0, tstr.find_last_not_of(" \n\r\t")+1);
    if(!tstr.empty()) { tstr = tstr.substr(tstr.find_first_not_of(" \n\r\t")); }
    return tstr;
}

bool Optimizer::makedir(const string& str) const
{
#if defined(_WIN32)
    return (CreateDirectory(str.c_str(), NULL) == 1) ? true : false;
#else
    return (mkdir(str.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0) ? true : false;
#endif
}
