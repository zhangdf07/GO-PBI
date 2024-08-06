#include "Builder.h"
#include "Graphene.h"
#include "GrInputFile.h"
#include "tinyfuns.h"
#include <algorithm>
#include <vector>
#include <string>
#include <cstdio>

using namespace nkchem;
using namespace std;

Builder::Builder(void) 
{
    /* Nothing to do */ 
}

Builder::~Builder(void) 
{
    /* Nothing to do */ 
}

void Builder::PrintHelp(void)
{
    printf("grbuilder: A program for generating coordinates and topologies of functionalized graphene.\n");
    printf("Author: Jun Zhang; zhangjunqcc@gmail.com\n");
    printf("\n");
    printf("Usage: grbuilder type inputfilename\n");
    printf("\n");
    printf("An example of the input file is:\n");
    printf("\n"
           "7 4 1.4 0 GRA      # pristine: nx, ny, R(C-C), Z\n"
           "12 R 6             # OH\n"
           "3 1 2 6 9 -11 -12  # Oe\n"
           "3 E                # H\n"
           "0                  # COOH\n");    
}

void Builder::AddGroups(Graphene& gr, int num, int num_up, GrInputFile::Style style, vector<int>& indices, const string& group)
{
    const double RoundAngle = 360.;
    const int num_grapheneC = gr.GetNumGrapheneC();
    vector<int> Cindices;
    for(int i = 0; i < num_grapheneC; ++i) { Cindices.push_back(i); }

    if(num > 0)
    {
        if(style == GrInputFile::Random) // Random distribution.
        {
            indices.clear();
            random_shuffle(Cindices.begin(), Cindices.end());
            for(int i = 0, iC = 0; i < num; ++i)
            {
                while(true)
                {
                    if(group != "Oe")
                    {
                        const int idx = Cindices[iC];
                        if(gr.CanBeFunctionalized(idx))
                        {
                            if(i < num_up)
                            {
                                if(group == "OH")   { gr.AddOH(idx, Graphene::Up, RoundAngle*rand01()); }
                                if(group == "H")    { gr.AddH(idx, Graphene::Up); }
                                if(group == "COOH") { gr.AddCOOH(idx, Graphene::Up, RoundAngle*rand01()); }
                                indices.push_back(idx+1);
                            }
                            else
                            {
                                if(group == "OH")   { gr.AddOH(idx, Graphene::Down, RoundAngle*rand01()); }
                                if(group == "H")    { gr.AddH(idx, Graphene::Down); }
                                if(group == "COOH") { gr.AddCOOH(idx, Graphene::Down, RoundAngle*rand01()); }
                                indices.push_back(-(idx+1));
                            }                        
                            break;
                        }
                        else
                        {
                            if(iC == num_grapheneC-1)
                            {
                                ErrorExit("No available graphene carbon atoms for functionalization.");
                            }
                            ++iC;
                        }
                    }
                    else // group == "Oe"
                    {
                        const int idx = Cindices[iC];
                        if(gr.CanBeFunctionalized(idx))
                        {
                            vector<int> Cindices1; gr.GetConnectedC(idx, Cindices1);
                            random_shuffle(Cindices1.begin(), Cindices1.end());
                            bool succeed = false;
                            for(int j = 0; j < Cindices1.size(); ++j)
                            {
                                const int idx1 = Cindices1[j];
                                if(gr.CanBeFunctionalized(idx1))
                                {
                                    if(i < num_up)
                                    {
                                        gr.AddOe(idx, idx1, Graphene::Up);
                                        indices.push_back(idx+1);
                                        indices.push_back(idx1+1);
                                    }
                                    else
                                    {
                                        gr.AddOe(idx, idx1, Graphene::Down);
                                        indices.push_back(-(idx+1));
                                        indices.push_back(-(idx1+1));
                                    }
                                    succeed = true;
                                    break;
                                }
                            }
                            ++iC;                            
                            if(succeed) { break; }
                        }
                        else
                        {
                            if(iC == num_grapheneC-1)
                            {
                                ErrorExit("No available graphene carbon atoms for functionalization.");
                            }
                            ++iC;
                        }    
                    }
                }
            }
        }
        else
        {
            if(style == GrInputFile::Evenly) // Evenly distribution.   
            {
                indices.clear();
                const int v = num_grapheneC/num;
                const int i0 = (num_grapheneC-v*(num-1))/2;
                for(int i = 0; i < num; ++i)
                {                
                    const int idx = i*v+i0;
                    if(group == "OH")   { gr.AddOH(idx, Graphene::Up, RoundAngle*rand01()); }
                    if(group == "Oe")   { gr.AddOe(idx, idx+1, Graphene::Up); indices.push_back(idx+2); }
                    if(group == "H")    { gr.AddH(idx, Graphene::Up); }
                    if(group == "COOH") { gr.AddCOOH(idx, Graphene::Up, RoundAngle*rand01()); }
                    indices.push_back(idx+1);
                }
            }
            else // Assigned ones.
            {
                for(int i = 0; i < num; ++i)
                {    
                    int first_index;
                    int idx, idx1;
                    if(group != "Oe")
                    {
                        first_index = indices[i];
                        idx = abs(indices[i])-1;
                    }
                    else
                    {
                        first_index = indices[i*2];
                        idx = abs(indices[i*2])-1;   
                        idx1 = abs(indices[i*2+1])-1;
                        if(gr.GetBondOrder(idx, idx1) == 0)
                        {
                            ErrorExit("Cannot Functionalize atom %d and %d with epoxy oxygen: they are not connected.",
                                idx+1, idx1+1);
                        }
                    }
                    if(first_index > 0)
                    {
                        if(group == "OH")   { gr.AddOH(idx, Graphene::Up, RoundAngle*rand01()); }
                        if(group == "Oe")   { gr.AddOe(idx, idx1, Graphene::Up); }
                        if(group == "H")    { gr.AddH(idx, Graphene::Up); }
                        if(group == "COOH") { gr.AddCOOH(idx, Graphene::Up, RoundAngle*rand01()); }
                    }
                    else
                    {
                        if(group == "OH")   { gr.AddOH(idx, Graphene::Down, RoundAngle*rand01()); }
                        if(group == "Oe")   { gr.AddOe(idx, idx1, Graphene::Down); }
                        if(group == "H")    { gr.AddH(idx, Graphene::Down); }
                        if(group == "COOH") { gr.AddCOOH(idx, Graphene::Down, RoundAngle*rand01()); }
                    }
                }
            }
        }
    }
}

