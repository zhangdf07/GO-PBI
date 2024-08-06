#include "Graphene.h"
#include "Vec3D.h"
#include "tinyfuns.h"
#include <algorithm>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <iostream>

using namespace nkchem;
using namespace std;

Graphene::Graphene(int tnx, int tny, double trCC, double tz)
{
    BuildSheet(tnx, tny, trCC, tz);
}

Graphene::~Graphene(void)
{
    // Nothing to do.
}

int Graphene::GetNumGrapheneC(void) const
{
    return num_graphene_C;
}

void Graphene::Export(const string& prefix) const
{
    ExportGRO(prefix);
    ExportXYZ(prefix);
    ExportITP(prefix);
}

bool Graphene::CanBeFunctionalized(int idx) const
{
    if(idx >= num_graphene_C || idx < 0 || atoms[idx].type != Graphene_C.type) 
    {
        return false;
    }
    return true;
}

void Graphene::CheckAddValidity(int idx) const
{
    if(idx >= num_graphene_C) 
    {
        ErrorExit("idx(%d) >= num_graphene_C(%d).", idx+1, num_graphene_C);
    }
    if(idx < 0)
    {
        ErrorExit("idx(%d) < 0.", idx+1);
    }
    if(atoms[idx].type != Graphene_C.type)
    {
        ErrorExit("idx(%d) is not Graphene_C(%s) but already functionalized (%s).", 
            idx+1, Graphene_C.type.c_str(), atoms[idx].type.c_str());
    }    
}

int Graphene::GetBondOrder(int idx1, int idx2) const
{
    CheckAddValidity(idx1);
    CheckAddValidity(idx2);
    return (int)(bond_matrix[idx1][idx2]);
}

void Graphene::GetConnectedC(int idx, vector<int>& indices) const
{
    CheckAddValidity(idx);
    indices.clear();
    for(int i = 0; i < num_graphene_C; ++i)
    {
        if(bond_matrix[idx][i] != None) { indices.push_back(i); }
    }
}

void Graphene::AddH(int idx, StereoType st)
{
    const double dCH = 1.08; 
    CheckAddValidity(idx);

    const int prev_num_atoms = atoms.size(); // Get number of atoms before addition of atoms!     
    const int num_group_atoms = 1;
    // Add H, go! 
    Atom graphene_atm = Alkane_C; graphene_atm.coord = atoms[idx].coord;
    Atom connected_H = Alkane_H;
    if(st == Up) { connected_H.coord = Vec3D(atoms[idx].coord.x, atoms[idx].coord.y, atoms[idx].coord.z+dCH); }
    else { connected_H.coord = Vec3D(atoms[idx].coord.x, atoms[idx].coord.y, atoms[idx].coord.z-dCH); }
    AttachAtomToGraphene(idx, graphene_atm, connected_H);
    // Update topology.
    for(int i = 0; i < prev_num_atoms; ++i) { for(int j = 0; j < num_group_atoms; ++j) { bond_matrix[i].push_back(None); } }
    for(int i = 0; i < num_group_atoms; ++i) { bond_matrix.push_back(vector<BondOrder>(prev_num_atoms+num_group_atoms, None)); }
    bond_matrix[idx][prev_num_atoms+0] = Single; bond_matrix[prev_num_atoms+0][idx] = Single;
}

void Graphene::AddOH(int idx, StereoType st, double rot_ang)
{
    const double dCO = 1.40;
    const double dOH = 1.;
    const double aHOC = 120.*Ang2Rad;
    CheckAddValidity(idx);

    const int prev_num_atoms = atoms.size(); // Get number of atoms before addition of atoms!     
    const int num_group_atoms = 2;
    // Add O, go! 
    Atom graphene_atm = Hydroxyl_C; graphene_atm.coord = atoms[idx].coord;
    Atom connected_O = Hydroxyl_O;
    if(st == Up) { connected_O.coord = Vec3D(atoms[idx].coord.x, atoms[idx].coord.y, atoms[idx].coord.z+dCO); }
    else { connected_O.coord = Vec3D(atoms[idx].coord.x, atoms[idx].coord.y, atoms[idx].coord.z-dCO); }
    AttachAtomToGraphene(idx, graphene_atm, connected_O);
    // Add H, go!
    Atom connected_H = Hydroxyl_H;
    const double rOHsinaHOC = dOH*sin(aHOC);
    const double rOHcosaHOC = dOH*cos(aHOC);
    rot_ang *= Ang2Rad;
    if(st == Up) { connected_H.coord = Vec3D(atoms[idx].coord.x+rOHsinaHOC*cos(rot_ang), atoms[idx].coord.y+rOHsinaHOC*sin(rot_ang), atoms[idx].coord.z+dCO-rOHcosaHOC); }
    else { connected_H.coord = Vec3D(atoms[idx].coord.x+rOHsinaHOC*cos(rot_ang), atoms[idx].coord.y+rOHsinaHOC*sin(rot_ang), atoms[idx].coord.z-dCO+rOHcosaHOC); }    
    atoms.push_back(connected_H);
    // Update topology.
    for(int i = 0; i < prev_num_atoms; ++i) { for(int j = 0; j < num_group_atoms; ++j) { bond_matrix[i].push_back(None); } }
    for(int i = 0; i < num_group_atoms; ++i) { bond_matrix.push_back(vector<BondOrder>(prev_num_atoms+num_group_atoms, None)); }
    bond_matrix[idx][prev_num_atoms+0] = Single; bond_matrix[prev_num_atoms+0][idx] = Single;
    bond_matrix[prev_num_atoms+0][prev_num_atoms+1] = Single; bond_matrix[prev_num_atoms+1][prev_num_atoms+0] = Single;
}


void Graphene::AddCOOH(int idx, StereoType st, double rot_ang)
{
    const double dCC = 1.5;
    const double dCO = 1.5;
    const double dOH = 1.;
    const double aOCC = 120.*Ang2Rad;
    CheckAddValidity(idx);

    const int prev_num_atoms = atoms.size(); // Get number of atoms before addition of atoms!     
    const int num_group_atoms = 4;
    // Add C, go! 
    Atom graphene_atm = COOH_C; graphene_atm.coord = atoms[idx].coord; 
    Atom connected_C = Carboxyl_C;
    if(st == Up) { connected_C.coord = Vec3D(atoms[idx].coord.x, atoms[idx].coord.y, atoms[idx].coord.z+dCC); }
    else { connected_C.coord = Vec3D(atoms[idx].coord.x, atoms[idx].coord.y, atoms[idx].coord.z-dCC); }
    AttachAtomToGraphene(idx, graphene_atm, connected_C);    
    // Add O, go!
    Atom connected_O = Carboxyl_O;
    const double rCOsinaOCC = dCO*sin(aOCC);
    const double rCOcosaOCC = dCO*cos(aOCC);
    rot_ang *= Ang2Rad;
    if(st == Up) { connected_O.coord = Vec3D(atoms[idx].coord.x+rCOsinaOCC*cos(rot_ang), atoms[idx].coord.y+rCOsinaOCC*sin(rot_ang), atoms[idx].coord.z+dCC-rCOcosaOCC); }
    else { connected_O.coord = Vec3D(atoms[idx].coord.x+rCOsinaOCC*cos(rot_ang), atoms[idx].coord.y+rCOsinaOCC*sin(rot_ang), atoms[idx].coord.z-dCC+rCOcosaOCC); }    
    atoms.push_back(connected_O);
    // Add Oh, go!
    const double Pi = 3.1415926535;
    Atom connected_Oh = Carboxyl_Oh;
    if(st == Up) { connected_Oh.coord = Vec3D(atoms[idx].coord.x+rCOsinaOCC*cos(Pi+rot_ang), atoms[idx].coord.y+rCOsinaOCC*sin(Pi+rot_ang), atoms[idx].coord.z+dCC-rCOcosaOCC); }
    else { connected_Oh.coord = Vec3D(atoms[idx].coord.x+rCOsinaOCC*cos(Pi+rot_ang), atoms[idx].coord.y+rCOsinaOCC*sin(Pi+rot_ang), atoms[idx].coord.z-dCC+rCOcosaOCC); }    
    atoms.push_back(connected_Oh);
    // Add H, go!
    Atom connected_H = Carboxyl_H;
    if(st == Up) { connected_H.coord = Vec3D(connected_Oh.coord.x, connected_Oh.coord.y, connected_Oh.coord.z+dOH); }
    else { connected_H.coord = Vec3D(connected_Oh.coord.x, connected_Oh.coord.y, connected_Oh.coord.z-dOH); }
    atoms.push_back(connected_H);
    //printf("Add H: %f\n", connected_H)
    // Update topology.
    for(int i = 0; i < prev_num_atoms; ++i) { for(int j = 0; j < num_group_atoms; ++j) { bond_matrix[i].push_back(None); } }
    for(int i = 0; i < num_group_atoms; ++i) { bond_matrix.push_back(vector<BondOrder>(prev_num_atoms+num_group_atoms, None)); }
    bond_matrix[idx][prev_num_atoms+0] = Single; bond_matrix[prev_num_atoms+0][idx] = Single;
    bond_matrix[prev_num_atoms+0][prev_num_atoms+1] = Double; bond_matrix[prev_num_atoms+1][prev_num_atoms+0] = Double;
    bond_matrix[prev_num_atoms+0][prev_num_atoms+2] = Single; bond_matrix[prev_num_atoms+2][prev_num_atoms+0] = Single;
    bond_matrix[prev_num_atoms+2][prev_num_atoms+3] = Single; bond_matrix[prev_num_atoms+3][prev_num_atoms+2] = Single;
}

void Graphene::AddOe(int idx1, int idx2, StereoType st)
{
    const double dCO = 1.45; 
    CheckAddValidity(idx1);
    CheckAddValidity(idx2);

    const int prev_num_atoms = atoms.size(); // Get number of atoms before addition of atoms!     
    const int num_group_atoms = 1;
    // Add O, go! 
    Atom graphene_atm1 = Epoxy_C; graphene_atm1.coord = atoms[idx1].coord;
    Atom graphene_atm2 = Epoxy_C; graphene_atm2.coord = atoms[idx2].coord;
    Atom connected_O = Epoxy_O;
    Vec3D xyO = (atoms[idx1].coord+atoms[idx2].coord)/2;
    const Vec3D dxy = atoms[idx1].coord-atoms[idx2].coord;   // PBC
    if(abs(dxy.x) > cell_vec.x/2) { xyO.x += cell_vec.x/2; } // PBC
    if(abs(dxy.y) > cell_vec.y/2) { xyO.y += cell_vec.y/2; } // PBC
    const double hO = sqrt(dCO*dCO-rCC*rCC/4);
    if(st == Up) { connected_O.coord = Vec3D(xyO.x, xyO.y, xyO.z+hO); }
    else  { connected_O.coord = Vec3D(xyO.x, xyO.y, xyO.z-hO); }
    AttachAtomToGraphene(idx1, graphene_atm1, idx2, graphene_atm2, connected_O);
    // Update topology.
    for(int i = 0; i < prev_num_atoms; ++i) { for(int j = 0; j < num_group_atoms; ++j) { bond_matrix[i].push_back(None); } }
    for(int i = 0; i < num_group_atoms; ++i) { bond_matrix.push_back(vector<BondOrder>(prev_num_atoms+num_group_atoms, None)); }
    bond_matrix[idx1][prev_num_atoms+0] = Single; bond_matrix[prev_num_atoms+0][idx1] = Single;   
    bond_matrix[idx2][prev_num_atoms+0] = Single; bond_matrix[prev_num_atoms+0][idx2] = Single;   
}

void Graphene::AttachAtomToGraphene(int idx, const Atom& graphene_atm, const Atom& connected_atm)
{
    atoms.push_back(connected_atm);
    atoms[idx]= graphene_atm;
}

void Graphene::AttachAtomToGraphene(int idx1, const Atom& graphene_atm1,
    int idx2, const Atom& graphene_atm2, const Atom& connected_atm)
{
    atoms.push_back(connected_atm);
    atoms[idx1]= graphene_atm1;
    atoms[idx2]= graphene_atm2;
}

void Graphene::BuildSheet(int tnx, int tny, double trCC, double tz)
{
    // Usually this is the first function to call, so we erase all previous possible data.    
    nx = tnx;
    ny = tny;
    rCC = trCC;
    
    // nx = 3
    // ny = 1
    //       C2     C2     C2
    //     /    \  /   \  /
    //   C1      C1     C1
    //    |      |      |
    //   C3      C3     C3
    //     \    /  \   /  \
    //       C4     C4     C4
    //
    //  i = 0, 1, ..., nx-1
    //  j = 0, 1, ..., ny-1
    //  C1: 4nx*j+2i
    //  C2: 4nx*j+2i+1
    //  C3: 4nx*j+2nx+2i
    //  C4: 4nx*j+2nx+2i+1

    // Add geometry information.
    atoms.clear();
    const double rCCCos30 = 0.8660254037844387*rCC;
    
    // Cell.
    cell_vec.x = nx*2*rCCCos30;
    cell_vec.y = ny*3*rCC;
    cell_vec.z = (cell_vec.x+cell_vec.y); // Intersheet distance > 3*wide.
    const double z = tz;
    num_graphene_C = 4*nx*ny;

    // Add atoms.
    for(int j = 0; j < ny; ++j)
    {
        for(int i = 0; i < nx; ++i)
        {
            atoms.push_back(Graphene_C); atoms.back().coord = Vec3D(i*2*rCCCos30, j*3*rCC, z);
            atoms.push_back(Graphene_C); atoms.back().coord = Vec3D((i*2+1)*rCCCos30, (j*3-0.5)*rCC, z);
        }
        for(int i = 0; i < nx; ++i)
        {
            atoms.push_back(Graphene_C); atoms.back().coord = Vec3D(i*2*rCCCos30 , (j*3+1)*rCC, z);
            atoms.push_back(Graphene_C); atoms.back().coord = Vec3D((i*2+1)*rCCCos30, (j*3+1.5)*rCC, z);
        }        
    }

    // Add topology information.
    bond_matrix.clear();
    const int num_atoms = atoms.size();
    bond_matrix.assign(num_atoms, vector<BondOrder>(num_atoms, None));
    int num_bonds = 0;
    for(int j = 0; j < ny; ++j)
    {
        // Bonds in the cell.
        for(int i = 0; i < nx*2; ++i)
        {
            const int idx = j*nx*4+i*2;
            bond_matrix[idx][idx+1] = Aromatic; bond_matrix[idx+1][idx] = Aromatic;
            ++num_bonds;
        }
        for(int i = 0; i < nx; ++i)
        {
            const int idx = j*nx*4+i*2;
            bond_matrix[idx][idx+nx*2] = Aromatic; bond_matrix[idx+nx*2][idx] = Aromatic;
            ++num_bonds;            
        }
        for(int i = 0; i < nx-1; ++i)
        {
            const int idx = j*nx*4+i*2+1;
            bond_matrix[idx][idx+1] = Aromatic; bond_matrix[idx+1][idx] = Aromatic;
            ++num_bonds;            
        }
        for(int i = 0; i < nx-1; ++i)
        {
            const int idx = j*nx*4+nx*2+i*2+1;
            bond_matrix[idx][idx+1] = Aromatic; bond_matrix[idx+1][idx] = Aromatic;
            ++num_bonds;            
        }
        if(j < ny-1)
        {
            for(int i = 0; i < nx; ++i)
            {
                const int idx = j*nx*4+i*2+1;
                bond_matrix[idx+nx*2][idx+nx*4] = Aromatic; bond_matrix[idx+nx*4][idx+nx*2] = Aromatic;
                ++num_bonds;                
            }            
        }
        // Bonds across the cell.
        #if 1
        {
            const int idx = j*nx*4;
            bond_matrix[idx][idx+nx*2-1] = Aromatic; bond_matrix[idx+nx*2-1][idx] = Aromatic;
            bond_matrix[idx+nx*2][idx+nx*4-1] = Aromatic; bond_matrix[idx+nx*4-1][idx+nx*2] = Aromatic;
            num_bonds += 2;
        }
        if(j == ny-1)
        {
            for(int i = 0; i < nx; ++i)
            {
                bond_matrix[i*2+1][(ny-1)*nx*4+nx*2+i*2+1] = Aromatic; bond_matrix[(ny-1)*nx*4+nx*2+i*2+1][i*2+1] = Aromatic;
                ++num_bonds;                
            }        
        }    
        #endif    
    }
}

void Graphene::ExportGRO(const string& prefix) const
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
    fprintf(fd, "%15.8f %15.8f %15.8f\n", cell_vec.x*AA2nm, cell_vec.y*AA2nm, cell_vec.z*AA2nm);
    fclose(fd);
}

void Graphene::ExportXYZ(const string& prefix) const
{
    std::cout << "Starting exporting XYZ file \n";
    const string fn = prefix+".xyz";
    FILE* fd = fopen(fn.c_str(), "w");
    if(fd == NULL) 
    {
        ErrorExit("Cannot write file %s.\n", fn.c_str());
    }
    const int num_atoms = atoms.size();

    fprintf(fd, "%d\n", num_atoms);

    fprintf(fd, "Cell: %15.8f %15.8f %15.8f\n", cell_vec.x, cell_vec.y, cell_vec.z);

    for(int i = 0; i < num_atoms; ++i)
    {
        fprintf(fd, "%5s  %18.10f %18.10f %18.10f\n", atoms[i].symbol.c_str(), atoms[i].coord.x, atoms[i].coord.y, atoms[i].coord.z);
    }
    
    fclose(fd);
}

void Graphene::ExportITP(const string& prefix) const
{
    std::cout << "Starting exporting ITP file \n";
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
