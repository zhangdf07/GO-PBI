#ifndef    __GRFLAKE_H__
#define    __GRFLAKE_H__

#include "Vec3D.h"
#include <string>
#include <vector>
#include <utility>

namespace nkchem {

using namespace std;

class GrFlake {
public:
    struct Bond {
        int idx0;
        int idx1;
        
        bool operator < (const Bond& b0) const
        {
            const int p00 = min(idx0, idx1);
            const int p01 = max(idx0, idx1);
            const int p10 = min(b0.idx0, b0.idx1);
            const int p11 = max(b0.idx0, b0.idx1);
            if(p00 < p10) { return true; }
            if(p00 == p10) { if(p01 < p11) { return true;} }
            return false;
        }
        
        bool operator == (const Bond& b0) const
        {
            return (b0.idx0 == idx0) && ((b0.idx1 == idx1)) ||
                (b0.idx0 == idx1) && ((b0.idx1 == idx0));            
        }
    };


public:
    enum SheetType { Parallel, Zigzag };
    enum StereoType { Up, Down };


public:
    GrFlake(int tnx, int tny, double trCC, double tz, SheetType tsheet_type);
    ~GrFlake(void);

    vector<int> GetAvailableEdgeAtoms(void) const;
    vector<int> GetAvailableInternalAtoms(void) const;
    vector<Bond> GetAvailableEdgeBonds(void) const;
    vector<Bond> GetAvailableInternalBonds(void) const;

    
    void Export(const string& prefix) const;
    const vector<int>& GetEdgeAtomIndices(void) const;
    const vector<int>& GetInternalAtomIndices(void) const;
    vector<int> GetBondingAtoms(int idx) const;
    bool IsEdge(int idx) const;

    bool AddH(int idx, StereoType st);
    bool AddOH(int idx, StereoType st, double phi);
    bool AddCOOH(int idx, StereoType st, double phi);
    bool AddOe(int idx0, int idx1, StereoType st);

private:
    // Unit: Angstrom everywhere!!!
    struct Atom {
        string symbol;
        string type;
        double mass;
        double charge;
        Vec3D  coord;
    };
    enum BondOrder { None = 0, Single, Double, Triple, Aromatic };
    const double Ang2Rad = 3.1415926535/180.;

    // Topology.
    const Atom Graphene_C =   { string("C"), string("opls_145"),  12.011,     0., Vec3D(0, 0, 0) };
    const Atom Alkane_C =     { string("C"), string("opls_137"),  12.011,  -0.06, Vec3D(0, 0, 0) };
    const Atom COOH_C =       { string("C"), string("opls_137"),  12.011,     0., Vec3D(0, 0, 0) };
    const Atom Carboxyl_C =   { string("C"), string("opls_267"),  12.011,  +0.52, Vec3D(0, 0, 0) };
    const Atom Hydroxyl_C =   { string("C"), string("opls_166"),  12.011,  +0.15, Vec3D(0, 0, 0) };
    const Atom Epoxy_C =      { string("C"), string("opls_184"),  12.011,  +0.20, Vec3D(0, 0, 0) };
    
    const Atom Alkane_H =   { string("H"), string("opls_140"),   1.008,  +0.06, Vec3D(0, 0, 0) };
    const Atom Hydroxyl_H = { string("H"), string("opls_163"),   1.008, +0.435, Vec3D(0, 0, 0) };
    const Atom Carboxyl_H = { string("H"), string("opls_270"),   1.008,  +0.45, Vec3D(0, 0, 0) };

    const Atom Hydroxyl_O =  { string("O"), string("opls_162"), 15.9994, -0.585, Vec3D(0, 0, 0) };
    const Atom Carboxyl_O =  { string("O"), string("opls_269"), 15.9994, -0.44, Vec3D(0, 0, 0) };
    const Atom Carboxyl_Oh = { string("O"), string("opls_268"), 15.9994, -0.53, Vec3D(0, 0, 0) };
    const Atom Epoxy_O =     { string("O"), string("opls_180"), 15.9994, -0.40, Vec3D(0, 0, 0) };

    // Geometry.
    const double rCH = 1.08;
    const double rCO = 1.41;
    const double rOH = 0.97;
    const double hOeG = 1.0;
    const double aHOC = 17.7*Ang2Rad;
    const double aOCC = 32.8*Ang2Rad;

private:
    int                        nx;
    int                        ny;
    double                     rCC;
    double                     z;
    int                        num_graphene_C;
    SheetType                  sheet_type;
    vector<Atom>               atoms;
    vector<vector<BondOrder> > bond_matrix;
    vector<int>                edge_atom_indices;
    vector<int>                internal_atom_indices;
    
    void BuildParallelSheet(void);
    void BuildZigZagSheet(void);
    void ExportGRO(const string& prefix) const;
    void ExportXYZ(const string& prefix) const;
    void ExportITP(const string& prefix) const;

};    

}

#endif //  __GRFLAKE_H__