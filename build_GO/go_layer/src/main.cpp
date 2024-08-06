#include "Optimizer.h"
#include "Graphene.h"
#include "GrFlake.h"
#include "GrInputFile.h"
#include "FlInputFile.h"
#include "FlBuilder.h"
#include "Builder.h"
#include "tinyfuns.h"
#include <string>
#include <cstdlib>
#include <cstdio>

using namespace nkchem;
using namespace std;

#include "Molecule.h"

int main(int argc, char** argv)
{
    nkchem::srand();    
    const string fn(argv[1]);
    //Optimizer opter;
    //opter.Do(fn);

    const string option("PBC");





    if(option == "PBC")
    {
        GrInputFile inp(fn);
        // Build pristine graphene.
        Graphene gr(inp.nx, inp.ny, inp.rCC, inp.z);
        // Hydroxyl OH.
        Builder::AddGroups(gr, inp.num_OH,   inp.num_OH_up,   inp.style_OH,   inp.OH_indices,   string("OH"));
        Builder::AddGroups(gr, inp.num_Oe,   inp.num_Oe_up,   inp.style_Oe,   inp.Oe_indices,   string("Oe"));
        Builder::AddGroups(gr, inp.num_H,    inp.num_H_up,    inp.style_H,    inp.H_indices,    string("H"));
        Builder::AddGroups(gr, inp.num_COOH, inp.num_COOH_up, inp.style_COOH, inp.COOH_indices, string("COOH"));
        // Export them.
        gr.Export(inp.resname.c_str());
        inp.Export((inp.resname+".inp").c_str());
    }
#if 0
    if(option == "Flake")
    {
        FlInputFile inp(fn);
        // Build pristine flake.
        printf("Graphene flake: %d x %d, rCC = %.2f\n", inp.nx, inp.ny, inp.rCC);
        GrFlake flake(inp.nx, inp.ny, inp.rCC, inp.z, GrFlake::Parallel);
        vector<int> internal_atom_indices = flake.GetInternalAtomIndices();
        printf("Derocation:\n");
        printf(" Edge Oe:       %d\n", inp.num_edge_Oe);
        printf(" Internal Oe:   %d\n", inp.num_internal_Oe);
        printf(" Edge OH:       %d\n", inp.num_edge_OH);
        printf(" Internal OH:   %d\n", inp.num_internal_OH);
        printf(" Edge COOH:     %d\n", inp.num_edge_COOH);
        printf(" Internal COOH: %d\n", inp.num_internal_COOH);       



        // The seed groups.        
        FlBuilder::Build(inp, flake);
        // Export them.
        flake.Export(inp.resname.c_str());
        inp.Export((inp.resname+".inp").c_str());
        // Evaluate them.
        system("");
    }
#endif    

    return 0;
}