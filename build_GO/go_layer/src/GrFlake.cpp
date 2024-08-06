#include "GrFlake.h"
#include "Molecule.h"
#include "Vec3D.h"
#include "tinyfuns.h"
#include <algorithm>
#include <string>
#include <vector>
#include <utility>
#include <cstdio>
#include <cmath>
#include <cassert>

using namespace nkchem;
using namespace std;

GrFlake::GrFlake(int tnx, int tny, double trCC, double tz, SheetType tsheet_type) : 
    nx(tnx), ny(tny), rCC(trCC), z(tz), sheet_type(tsheet_type)
{
    if(tsheet_type == Parallel) { BuildParallelSheet(); }
}

GrFlake::~GrFlake(void)
{
    // Nothing to do.
}

void GrFlake::Export(const string& prefix) const
{
    ExportGRO(prefix);
    ExportXYZ(prefix);
    ExportITP(prefix);    
}

vector<int> GrFlake::GetAvailableEdgeAtoms(void) const
{
    vector<int> indices = edge_atom_indices;
    for(vector<int>::iterator iter = indices.begin(); iter != indices.end(); /* Nothing to do. */ )
    {
        if(GetBondingAtoms(*iter).size() >= 3) { iter = indices.erase(iter); }
        else { ++iter; }
    }
    return indices;
}

vector<int> GrFlake::GetAvailableInternalAtoms(void) const
{
    vector<int> indices = internal_atom_indices;
    for(vector<int>::iterator iter = indices.begin(); iter != indices.end(); /* Nothing to do. */ )
    {
        if(GetBondingAtoms(*iter).size() >= 4) { iter = indices.erase(iter); }
        else { ++iter; }
    }
    return indices;
}

vector<GrFlake::Bond> GrFlake::GetAvailableEdgeBonds(void) const
{
    vector<Bond> bonds;
    vector<int> indices = GetAvailableEdgeAtoms();
    const int num_available_atoms = indices.size();
    for(int i = 0; i < num_available_atoms; ++i)
    {
        const int idx0 = indices[i];
        vector<int> bonding_atoms = GetBondingAtoms(idx0);
        const int num_bonding_atoms = bonding_atoms.size();
        for(int j = 0; j < num_bonding_atoms; ++j)
        {
            const int idx1 = bonding_atoms[j];
            if(IsEdge(idx1) && GetBondingAtoms(idx1).size() <= 2)
            {
                bonds.push_back({idx0, idx1});
            }
        }
    }
    sort(bonds.begin(), bonds.end());
    bonds.erase(unique(bonds.begin(), bonds.end() ), bonds.end());
    return bonds;    
}

vector<GrFlake::Bond> GrFlake::GetAvailableInternalBonds(void) const
{
    vector<Bond> bonds;
    vector<int> indices = GetAvailableInternalAtoms();
    const int num_available_atoms = indices.size();
    for(int i = 0; i < num_available_atoms; ++i)
    {
        const int idx0 = indices[i];
        vector<int> bonding_atoms = GetBondingAtoms(idx0);
        const int num_bonding_atoms = bonding_atoms.size();
        for(int j = 0; j < num_bonding_atoms; ++j)
        {
            const int idx1 = bonding_atoms[j];
            if((IsEdge(idx1) && GetBondingAtoms(idx1).size() <= 2) || (!IsEdge(idx1) && GetBondingAtoms(idx1).size() <= 3))
            {
                bonds.push_back({idx0, idx1});
            }
        }
    }
    sort(bonds.begin(), bonds.end());
    bonds.erase(unique(bonds.begin(), bonds.end() ), bonds.end());
    return bonds;
}

const vector<int>& GrFlake::GetEdgeAtomIndices(void) const
{
    return edge_atom_indices;
}


const vector<int>& GrFlake::GetInternalAtomIndices(void) const
{
    return internal_atom_indices;
}

void GrFlake::BuildParallelSheet(void)
{
    // Add geometry information.
    atoms.clear();
    edge_atom_indices.clear();
    const double rCCCos30 = 0.8660254037844387*rCC;

    // Add atoms.
    // Low edge.
    int atom_idx = 0;
    for(int i = 0; i < nx+1; ++i)
    {
        atoms.push_back(Graphene_C); atoms.back().coord = Vec3D(i*2*rCCCos30, 0., z);
        ++atom_idx;
        if(i != nx)
        {
            atoms.push_back(Graphene_C); atoms.back().coord = Vec3D((i*2+1)*rCCCos30, -0.5*rCC, z);
            ++atom_idx;
        }
    }    
    // Internal edges.
    for(int i = 0; i < ny-1; ++i)
    {
        for(int j = 0; j < nx+1; ++j)
        {
            atoms.push_back(Graphene_C); atoms.back().coord = Vec3D((j*2+i)*rCCCos30, (1+1.5*i)*rCC, z);
            ++atom_idx;
            atoms.push_back(Graphene_C); atoms.back().coord = Vec3D((j*2+1+i)*rCCCos30, (1.5+1.5*i)*rCC, z);
            ++atom_idx;
        }
    }
    // High edge.
    for(int i = 0; i < nx+1; ++i)
    {
        atoms.push_back(Graphene_C); atoms.back().coord = Vec3D((i*2+ny-1)*rCCCos30, (1.5*ny-0.5)*rCC, z);
        ++atom_idx;
        if(i != nx)
        {
            atoms.push_back(Graphene_C); atoms.back().coord = Vec3D((i*2+ny)*rCCCos30, 1.5*ny*rCC, z);
            ++atom_idx;
        }
    }

    // Bonds.
    bond_matrix.assign(atom_idx, vector<BondOrder>(atom_idx, None));
    for(int i = 0; i < atom_idx; ++i)
    {
        for(int j = i+1; j < atom_idx; ++j)
        {
            if(atoms[i].coord.norm(atoms[j].coord) < rCC*1.01)
            {
                bond_matrix[i][j] = Aromatic;
                bond_matrix[j][i] = Aromatic;
            }
        }
    }
    
    // Number of atoms.
    num_graphene_C = atom_idx;    

    // Refine edge atoms.
    // for(vector<int>::iterator iter = edge_atom_indices.begin(); iter != edge_atom_indices.end(); )
    // {
    //     if(GetBondingAtoms(*iter).size() > 2) { iter = edge_atom_indices.erase(iter); }
    //     else { ++iter; }
    // }
    edge_atom_indices.clear();
    for(int i = 0; i < num_graphene_C; ++i)
    {
        if(GetBondingAtoms(i).size() == 2) { edge_atom_indices.push_back(i); }
    }
    internal_atom_indices.clear();
    for(int i = 0; i < num_graphene_C; ++i)
    {
        if(GetBondingAtoms(i).size() >= 3) { internal_atom_indices.push_back(i); }
    }
}

vector<int> GrFlake::GetBondingAtoms(int idx) const
{
    if(idx >= num_graphene_C)
    {
        ErrorExit("Invalid idx = %d since only %d graphene C atoms are available.", idx, num_graphene_C);
    }

    vector<int> indices;
    const int num_atoms = bond_matrix[0].size();
    for(int i = 0; i < num_atoms; ++i)
    {
        if(bond_matrix[i][idx] != None)
        {
            indices.push_back(i);
        }
    }    
    return indices;
}

bool GrFlake::IsEdge(int idx) const
{
    return find(edge_atom_indices.begin(), edge_atom_indices.end(), idx) != edge_atom_indices.end();
}

bool GrFlake::AddH(int idx, StereoType st)
{
    // Get positions.
    vector<int> indices = GetBondingAtoms(idx);
    const int num_bonds = indices.size();
    if(num_bonds >= 4) { return false; } // No groups can be added.
    if(num_bonds >= 3 && IsEdge(idx)) { return false; } // No groups can be added.
    Vec3D a1, u1;
    if(num_bonds == 3) { u1 = Vec3D(0., 0., (st == Up) ? +1. : -1.); }
    if(num_bonds == 2) { u1 = (atoms[idx].coord-(atoms[indices[0]].coord+atoms[indices[1]].coord)*0.5).normI(); }
    a1 = atoms[idx].coord+u1*rCH;

    // Set up groups, calculate coordinates.
    Molecule mol;
    vector<Vec3D> c = { Vec3D(0., 0., 0.) };
    vector<string> s = { string("H") };
    mol.SetCoordsAndSymbols(c, s);
    vector<Vec3D> coords = mol.CalcBondingCoords(0, 0, a1, u1);
    
    Atom atomH = Alkane_H;
    atomH.coord = coords[0];
    atoms.push_back(atomH);

    // Origin Carbon atoms.    
    Atom atomC = Alkane_C;
    atomC.coord = atoms[idx].coord;
    atoms[idx] = atomC;

    // Add bondings.
    const int num_atoms = bond_matrix.size();
    for(int i = 0; i < num_atoms; ++i)
    {
        bond_matrix[i].push_back((i == idx) ? Single : None);
    }
    bond_matrix.push_back(vector<BondOrder>(num_atoms+1, None));
    for(int i = 0; i < num_atoms+1; ++i)
    {
        bond_matrix[num_atoms][i] = bond_matrix[i][num_atoms];
    }

    // Done.
    return true;
}

bool GrFlake::AddOH(int idx, StereoType st, double phi)
{
    // Get positions.
    vector<int> indices = GetBondingAtoms(idx);
    const int num_bonds = indices.size();
    if(num_bonds >= 4) { return false; } // No groups can be added.
    if(num_bonds >= 3 && IsEdge(idx)) { return false; } // No groups can be added.
    Vec3D a1, u1;
    if(num_bonds == 3) { u1 = Vec3D(0., 0., (st == Up) ? +1. : -1.); }
    if(num_bonds == 2) { u1 = (atoms[idx].coord-(atoms[indices[0]].coord+atoms[indices[1]].coord)*0.5).normI(); }
    a1 = atoms[idx].coord+u1*rCO;

    // Set up groups, calculate coordinates.    
    Molecule mol;
    vector<Vec3D> c = { Vec3D(0., 0., 0.), Vec3D(0., rOH*cos(aHOC), rOH*sin(aHOC)), Vec3D(0., 0., 1.) };
    vector<string> s = { string("O"), string("H"), string("X") };    
    mol.SetCoordsAndSymbols(c, s);
    mol.RotZ(phi);
    vector<Vec3D> coords = mol.CalcBondingCoords(0, 2, a1, u1);
    
    
    Atom atomO = Hydroxyl_O;
    atomO.coord = coords[0];
    atoms.push_back(atomO);
    Atom atomH = Hydroxyl_H;
    atomH.coord = coords[1];
    atoms.push_back(atomH);

    // Origin Carbon atoms.    
    Atom atomC = Hydroxyl_C;
    atomC.coord = atoms[idx].coord;
    atoms[idx] = atomC;

    // Add bondings.
    const int num_atoms = bond_matrix.size();
    for(int i = 0; i < num_atoms; ++i)
    {
        bond_matrix[i].push_back((i == idx) ? Single : None);        
        bond_matrix[i].push_back(None);        
    }
    bond_matrix.push_back(vector<BondOrder>(num_atoms+2, None));
    bond_matrix.push_back(vector<BondOrder>(num_atoms+2, None));
    bond_matrix[num_atoms][num_atoms+1] = Single;
    bond_matrix[num_atoms+1][num_atoms] = Single;
    for(int i = 0; i < num_atoms+2; ++i)
    {
        bond_matrix[num_atoms][i] = bond_matrix[i][num_atoms];
        bond_matrix[num_atoms+1][i] = bond_matrix[i][num_atoms+1];
    }

    // Done.
    return true;
}

bool GrFlake::AddCOOH(int idx, StereoType st, double phi)
{
    // Get positions.
    vector<int> indices = GetBondingAtoms(idx);
    const int num_bonds = indices.size();
    if(num_bonds >= 4) { return false; } // No groups can be added.
    if(num_bonds >= 3 && IsEdge(idx)) { return false; } // No groups can be added.
    Vec3D a1, u1;
    if(num_bonds == 3) { u1 = Vec3D(0., 0., (st == Up) ? +1. : -1.); }
    if(num_bonds == 2) { u1 = (atoms[idx].coord-(atoms[indices[0]].coord+atoms[indices[1]].coord)*0.5).normI(); }
    a1 = atoms[idx].coord+u1*rCO;

    // Set up groups, calculate coordinates.    
    Molecule mol;
    vector<Vec3D> c = { Vec3D(0., 0., 0.),
        Vec3D(0., -rCO*cos(aOCC), rCO*sin(aOCC)),
        Vec3D(0., rCO*cos(aOCC), rCO*sin(aOCC)),        
        Vec3D(0., rOH*cos(aOCC), rOH+rCO*sin(aOCC)),
        Vec3D(0., 0., 1.) };
    vector<string> s = { string("C"), string("O"), string("O"), string("H"), string("X") };
    mol.SetCoordsAndSymbols(c, s);
    mol.RotZ(phi);
    vector<Vec3D> coords = mol.CalcBondingCoords(0, 4, a1, u1);
        
    Atom atomC = COOH_C;
    atomC.coord = coords[0];
    atoms.push_back(atomC);
    Atom atomO1 = Carboxyl_O;
    atomO1.coord = coords[1];
    atoms.push_back(atomO1);
    Atom atomO2 = Carboxyl_Oh;
    atomO2.coord = coords[2];
    atoms.push_back(atomO2);    
    Atom atomH = Hydroxyl_H;
    atomH.coord = coords[3];
    atoms.push_back(atomH);

    // Origin Carbon atoms.    
    atomC = Carboxyl_C;
    atomC.coord = atoms[idx].coord;
    atoms[idx] = atomC;

    // Add bondings.
    const int num_atoms = bond_matrix.size();
    for(int i = 0; i < num_atoms; ++i)
    {
        bond_matrix[i].push_back((i == idx) ? Single : None);        
        bond_matrix[i].push_back(None);        
        bond_matrix[i].push_back(None);        
        bond_matrix[i].push_back(None);        
    }
    bond_matrix.push_back(vector<BondOrder>(num_atoms+4, None));
    bond_matrix.push_back(vector<BondOrder>(num_atoms+4, None));
    bond_matrix.push_back(vector<BondOrder>(num_atoms+4, None));
    bond_matrix.push_back(vector<BondOrder>(num_atoms+4, None));

    bond_matrix[num_atoms][num_atoms+1] = Single;
    bond_matrix[num_atoms+1][num_atoms] = Single;
    bond_matrix[num_atoms][num_atoms+2] = Single;
    bond_matrix[num_atoms+2][num_atoms] = Single;
    bond_matrix[num_atoms+2][num_atoms+3] = Single;
    bond_matrix[num_atoms+3][num_atoms+2] = Single;    
    for(int i = 0; i < num_atoms+4; ++i)
    {
        bond_matrix[num_atoms][i] = bond_matrix[i][num_atoms];
        bond_matrix[num_atoms+1][i] = bond_matrix[i][num_atoms+1];
        bond_matrix[num_atoms+2][i] = bond_matrix[i][num_atoms+2];
        bond_matrix[num_atoms+3][i] = bond_matrix[i][num_atoms+3];
    }

    // Done.
    return true;
}

bool GrFlake::AddOe(int idx0, int idx1, StereoType st)
{
    // Get positions.
    if(!IsEdge(idx0) && GetBondingAtoms(idx0).size() >= 4) { return false; }
    if(!IsEdge(idx1) && GetBondingAtoms(idx1).size() >= 4) { return false; }
    if(IsEdge(idx0) && GetBondingAtoms(idx0).size() >= 3) { return false; }
    if(IsEdge(idx1) && GetBondingAtoms(idx1).size() >= 3) { return false; }
    if(bond_matrix[idx0][idx1] == None) { return false; }

    Atom atomO = Epoxy_O;
    atomO.coord = (atoms[idx0].coord+atoms[idx1].coord)/2;
    atomO.coord.z = (st == Up) ? +hOeG : -hOeG ;
    atoms.push_back(atomO);

    // Origin Carbon atoms.    
    Atom atomC = Epoxy_C;
    atomC.coord = atoms[idx0].coord;
    atoms[idx0] = atomC;
    atomC.coord = atoms[idx1].coord;
    atoms[idx1] = atomC;    

    // Add bondings.
    const int num_atoms = bond_matrix.size();
    for(int i = 0; i < num_atoms; ++i)
    {
        bond_matrix[i].push_back((i == idx0 || i == idx1) ? Single : None);
    }
    bond_matrix.push_back(vector<BondOrder>(num_atoms+1, None));
    for(int i = 0; i < num_atoms+1; ++i)
    {
        bond_matrix[num_atoms][i] = bond_matrix[i][num_atoms];
    }

    // Done.
    return true;
}

void GrFlake::ExportXYZ(const string& prefix) const
{
    const string fn = prefix+".xyz";
    FILE* fd = fopen(fn.c_str(), "w");
    if(fd == NULL) 
    {
        ErrorExit("Cannot write file %s.\n", fn.c_str());
    }
    const int num_atoms = atoms.size();
    fprintf(fd, "%d\n", num_atoms);
    fprintf(fd, "fragment: %d x %d, rCC = %.4f, %s\n", nx, ny, rCC, (sheet_type == Parallel) ? "Parallel" : "Zigzag");
    for(int i = 0; i < num_atoms; ++i)
    {
        fprintf(fd, "%5s  %18.10f %18.10f %18.10f\n", atoms[i].symbol.c_str(), atoms[i].coord.x, atoms[i].coord.y, atoms[i].coord.z);
    }    
    fclose(fd);
}

void GrFlake::ExportGRO(const string& prefix) const
{
    const string fn = prefix+".gro";
    FILE* fd = fopen(fn.c_str(), "w");
    if(fd == NULL) 
    {
        ErrorExit("Cannot write file %s.\n", fn.c_str());
    }
    const int num_atoms = atoms.size();

    fprintf(fd, "Generated by GRapheneBuilder.\n");
    fprintf(fd, "%d\n", num_atoms);
    const double AA2nm = 0.1;
    for(int i = 0; i < num_atoms; ++i)
    {
        fprintf(fd, "%5d%3s%7s%5d%8.3f%8.3f%8.3f\n", 1, prefix.c_str(), atoms[i].symbol.c_str(), i+1, atoms[i].coord.x*AA2nm, atoms[i].coord.y*AA2nm, atoms[i].coord.z*AA2nm);
    }
    fprintf(fd, "%15.8f %15.8f %15.8f\n", 10., 10., 10.);
    fclose(fd);
}

void GrFlake::ExportITP(const string& prefix) const
{
    const string fn = prefix+".itp";
    FILE* fd = fopen(fn.c_str(), "w");
    if(fd == NULL) 
    {
        ErrorExit("Cannot write file %s.\n", fn.c_str());
    }
    
    const int num_atoms = atoms.size();
    fprintf(fd, "; Generated by GRapheneBuilder.\n");
    fprintf(fd, "\n");

    fprintf(fd, "[ moleculetype ]\n");
    fprintf(fd, "; %-20s %10s\n", "Name", "nrexcl");
    fprintf(fd, "  %-20s %10d\n", prefix.c_str(), 3);
    fprintf(fd, "\n");

    fprintf(fd, "[ atoms ]\n");
    fprintf(fd, "; %8s %15s %8s %8s %8s %8s %8s %8s\n",
        "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass");
    for(int i = 0; i < num_atoms; ++i)
    {
        fprintf(fd, "  %8d %15s %8d %8s %8s %8d %8.3f %8.3f\n", 
            i+1, atoms[i].type.c_str(), 1, prefix.c_str(), atoms[i].symbol.c_str(), i+1, atoms[i].charge, atoms[i].mass);
    }
    fprintf(fd, "\n");

    fprintf(fd, "[ bonds ]\n");
    fprintf(fd, "; %5s %5s %5s\n", "ai", "aj", "funct");
    for(int i = 0; i < num_atoms; ++i)
    {
        for(int j = i+1; j < num_atoms; ++j)
        {
            if(bond_matrix[i][j] != None)
            {
                fprintf(fd, "  %5d %5d %5d    ; %4s %4s\n", i+1, j+1, 1, atoms[i].symbol.c_str(), atoms[j].symbol.c_str());
            }
        }        
    }    
    fprintf(fd, "\n");

    fprintf(fd, "[ pairs ]\n");
    fprintf(fd, "; %5s %5s %5s\n", "ai", "aj", "funct");
    for(int k = 0; k < num_atoms; ++k)
    {
        for(int l = k+1; l < num_atoms; ++l)
        {
            for(int i = 0; i < num_atoms; ++i)
            {
                if(i != k && i != l)
                {
                    for(int j = i+1; j < num_atoms; ++j)
                    {
                        if(j != k && j != l)
                        {
                            if(bond_matrix[k][l] != None && bond_matrix[k][i] != None && bond_matrix[l][j] != None)
                            {
                                fprintf(fd, "  %5d %5d %5d    ; %4s %4s\n", i+1, j+1, 1, atoms[i].symbol.c_str(), atoms[j].symbol.c_str());
                            }
                            if(bond_matrix[k][l] != None && bond_matrix[k][j] != None && bond_matrix[l][i] != None)
                            {
                                fprintf(fd, "  %5d %5d %5d    ; %4s %4s\n", i+1, j+1, 1, atoms[i].symbol.c_str(), atoms[j].symbol.c_str());
                            }
                        }
                    }
                }
            }
        }
    }    
    fprintf(fd, "\n");

    fprintf(fd, "[ angles ]\n");
    fprintf(fd, "; %5s %5s %5s %5s\n", "ai", "aj", "ak", "funct");
    for(int i = 0; i < num_atoms; ++i)
    {
        for(int j = 0; j < num_atoms; ++j)
        {
            for(int k = j+1; k < num_atoms; ++k)
            {
                if(bond_matrix[i][j] != None && bond_matrix[i][k] != None)
                {
                    fprintf(fd, "  %5d %5d %5d %5d    ; %4s %4s %4s\n", j+1, i+1, k+1, 1, atoms[j].symbol.c_str(), atoms[i].symbol.c_str(), atoms[k].symbol.c_str());
                }
            }
        }
    }
    fprintf(fd, "\n");

    fprintf(fd, "[ dihedrals ]\n");
    fprintf(fd, "; %5s %5s %5s %5s %5s\n", "ai", "aj", "ak", "al", "funct");
    for(int i = 0; i < num_atoms; ++i)
    {
        for(int j = i+1; j < num_atoms; ++j)
        {
            for(int k = 0; k < num_atoms; ++k)
            {
                if(k != i && k != j)
                {
                    for(int l = k+1; l < num_atoms; ++l)
                    {
                        if(l != i && l != j)
                        {
                            if(bond_matrix[i][j] != None && bond_matrix[i][k] != None && bond_matrix[j][l] != None)
                            {
                                fprintf(fd, "  %5d %5d %5d %5d %5d    ; %4s %4s %4s %4s\n", k+1, i+1, j+1, l+1, 3,
                                    atoms[k].symbol.c_str(), atoms[i].symbol.c_str(), atoms[j].symbol.c_str(), atoms[l].symbol.c_str());
                            }
                            if(bond_matrix[i][j] != None && bond_matrix[i][l] != None && bond_matrix[j][k] != None)
                            {
                                fprintf(fd, "  %5d %5d %5d %5d %5d    ; %4s %4s %4s %4s\n", l+1, i+1, j+1, k+1, 3,
                                    atoms[l].symbol.c_str(), atoms[i].symbol.c_str(), atoms[j].symbol.c_str(), atoms[k].symbol.c_str());
                            }
                        }
                    }
                }
            }
        }
    }    
    fprintf(fd, "\n");

    fclose(fd);
}
