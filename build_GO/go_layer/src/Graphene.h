#ifndef    __GRAPHENE_H__
#define    __GRAPHENE_H__

#include "Vec3D.h"
#include <string>
#include <vector>

namespace nkchem {

using namespace std;

class Graphene {
public:
    enum StereoType { Up, Down };

public:
    Graphene(int tnx, int tny, double trCC, double tz);
    ~Graphene(void);

    int GetNumGrapheneC(void) const;
    void Export(const string& prefix) const;
    void CheckAddValidity(int idx) const;
    bool CanBeFunctionalized(int idx) const;
    int GetBondOrder(int idx1, int idx2) const;
    void GetConnectedC(int idx, vector<int>& indices) const;

    // Add an H
    void AddH(int idx, StereoType st);
    // Add an OH, OH rotates rot_ang. 
    void AddOH(int idx, StereoType st, double rot_ang);
    // Add a COOH, COOH rotates rot_ang. 
    void AddCOOH(int idx, StereoType st, double rot_ang);
    // Add an O as an epoxy.
    void AddOe(int idx1, int idx2, StereoType st);

private:
    // Unit: Angstrom everywhere!!!
    struct Atom {
        string symbol;
        string type;
        double mass;
        double charge;
        Vec3D  coord;
    };
    
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
    
    enum BondOrder { None = 0, Single, Double, Triple, Aromatic };  

    const double Ang2Rad = 3.1415926535/180.;

private:
    int nx;
    int ny;
    double rCC;
    int num_graphene_C;

    vector<Atom>               atoms;
    vector<vector<BondOrder> > bond_matrix;
    Vec3D                      cell_vec;

    void BuildSheet(int tnx, int tny, double trCC, double tz);
    void AttachAtomToGraphene(int idx, const Atom& graphene_atm, const Atom& connected_atm);
    void AttachAtomToGraphene(int idx1, const Atom& graphene_atm1,
        int idx2, const Atom& graphene_atm2, const Atom& connected_atm);

    void ExportGRO(const string& prefix) const;
    void ExportXYZ(const string& prefix) const;
    void ExportITP(const string& prefix) const;
};

}

#endif //  __GRAPHENE_H__