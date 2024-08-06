#include "FlBuilder.h"
#include "GrFlake.h"
#include "FlInputFile.h"
#include "tinyfuns.h"
#include <algorithm>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>

using namespace nkchem;
using namespace std;

FlBuilder::FlBuilder(void) 
{
    /* Nothing to do */ 
}

FlBuilder::~FlBuilder(void) 
{
    /* Nothing to do */ 
}

void FlBuilder::Build(FlInputFile& inp, GrFlake& flake)
{
    const bool Edge = true;
    const bool Internal = false;
    FlBuilder::AddGroups("Oe", Edge, flake, inp.num_edge_Oe, inp.edge_Oe_indices);
    FlBuilder::AddGroups("Oe", Internal, flake, inp.num_internal_Oe, inp.internal_Oe_indices);
    FlBuilder::AddGroups("OH", Edge, flake, inp.num_edge_OH, inp.edge_OH_indices);
    FlBuilder::AddGroups("OH", Internal, flake, inp.num_internal_OH, inp.internal_OH_indices);
    FlBuilder::AddGroups("COOH", Edge, flake, inp.num_edge_COOH, inp.edge_COOH_indices);
    FlBuilder::AddGroups("COOH", Internal, flake, inp.num_internal_COOH, inp.internal_COOH_indices);
    FlBuilder::CloseWithH(flake);
}

void FlBuilder::CloseWithH(GrFlake& flake)
{
    const vector<int> edge_atom_indices = flake.GetEdgeAtomIndices();
    const int num_edge_atom_indices = edge_atom_indices.size();
    for(int i = 0; i < num_edge_atom_indices; ++i)
    { 
        flake.AddH(edge_atom_indices[i], GrFlake::Up);
    }  
}

void FlBuilder::AddGroups(const string& group_name, bool is_edge, GrFlake& flake, int num_groups, vector<int>& indices)
{
    if(indices.empty()) // Random
    {
        if(group_name != "Oe")
        {
            vector<int> all_indices = (is_edge) ? flake.GetAvailableEdgeAtoms() : flake.GetAvailableInternalAtoms();        
            random_shuffle(all_indices.begin(), all_indices.end());
            indices.clear();    
            if(all_indices.size() < num_groups) { ErrorExit("Insufficient number of available carbons."); }
            for(int i = 0; i < num_groups; ++i)
            {
                indices.push_back(all_indices[i]+1);
            }
        }
        else
        {
            vector<GrFlake::Bond> all_bonds = (is_edge) ? flake.GetAvailableEdgeBonds() : flake.GetAvailableInternalBonds();
            random_shuffle(all_bonds.begin(), all_bonds.end());
            if(all_bonds.size() < num_groups) { ErrorExit("Insufficient number of available carbons."); }
            vector<GrFlake::Bond> tbonds;
            for(int i = 0, p = 0; i < num_groups; ++i)
            {
                while(true)
                {
                    if(p >= all_bonds.size()) { ErrorExit("Insufficient number of available carbons."); }
                    const int idx0 = all_bonds[p].idx0;
                    const int idx1 = all_bonds[p].idx1;
                    bool bond_exists = false;
                    for(int j = 0; j < tbonds.size(); ++j)
                    {
                        if(idx0 == tbonds[j].idx0 || idx0 == tbonds[j].idx1 || idx1 == tbonds[j].idx0 || idx1 == tbonds[j].idx1)
                        {
                            bond_exists = true;
                            break;
                        }
                    }             
                    if(bond_exists)       
                    {
                        ++p;
                        continue;
                    }
                    else
                    {
                        tbonds.push_back({idx0, idx1});
                        ++p;
                        break;
                    }
                }                
            }            
            indices.clear();
            for(int i = 0; i < tbonds.size(); ++i)
            {
                indices.push_back(tbonds[i].idx0+1);
                indices.push_back(tbonds[i].idx1+1);
            }
        }
    }


    AddGroupsCore(group_name, is_edge, flake, num_groups, indices);
}

void FlBuilder::AddGroupsCore(const string& group_name, bool is_edge, GrFlake& flake, int num_groups, vector<int>& indices)
{
    const double PI2 = 6.28318530718;
    for(int i = 0; i < num_groups; ++i)
    {
        if(group_name != "Oe")
        {            
            const double phi = rand01()*PI2;
            const int idx = abs(indices[i])-1;
            if(is_edge && !flake.IsEdge(idx)) { ErrorExit("Requiring edge groups but it is not an edge one (%d)\n", indices[i]); }
            if(!is_edge && flake.IsEdge(idx)) { ErrorExit("Requiring internal groups but it is not an internal one (%d)\n", indices[i]); }
            const GrFlake::StereoType st = (indices[i] > 0) ? GrFlake::Up : GrFlake::Down ;
            if(group_name == "OH") { const bool succeed = flake.AddOH(idx, st, phi); if(!succeed) { ErrorExit("Cannot add OH at %d.", indices[i]); } }
            if(group_name == "COOH") { const bool succeed = flake.AddCOOH(idx, st, phi); if(!succeed) { ErrorExit("Cannot add OH at %d.", indices[i]); } }
        }
        else
        {
            const int idx0 = abs(indices[i*2])-1;
            const int idx1 = abs(indices[i*2+1])-1;
            if(is_edge && !(flake.IsEdge(idx0) && flake.IsEdge(idx1))) { ErrorExit("Requiring edge groups but it is not an edge one (%d, %d)\n", indices[i*2], indices[i*2+1]); }
            if(!is_edge && (flake.IsEdge(idx0) && flake.IsEdge(idx1))) { ErrorExit("Requiring internal groups but it is not an internal one (%d, %d)\n", indices[i*2], indices[i*2+1]); }
            const GrFlake::StereoType st = (indices[i*2] > 0 && indices[i*2+1] > 0) ? GrFlake::Up : GrFlake::Down ;
            const bool succeed = flake.AddOe(idx0, idx1, st);
            if(!succeed) { ErrorExit("Cannot add epoxy O at %d and %d.", indices[i*2], indices[i*2+1]); }
        }
    }
}