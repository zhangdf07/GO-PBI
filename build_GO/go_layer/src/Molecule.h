#ifndef    __MOLECULE_H__
#define    __MOLECULE_H__

#include "Vec3D.h"
#include <string>
#include <vector>

namespace nkchem {

using namespace std;

class Molecule {
public:
    Molecule(void);
    ~Molecule(void);

    void SetCoordsAndSymbols(const vector<Vec3D>& tcoords, const vector<string>& tsymbols);
    
    const vector<Vec3D>& GetCoords(void) const;
    const vector<string>& GetSymbols(void) const;
    vector<Vec3D> CalcBondingCoords(int bonding_idx, int directing_idx, const Vec3D& a1, const Vec3D& u1) const;
    void RotZ(double phi);
    
private:
    vector<Vec3D> coords;
    vector<string> symbols;

};

}

#endif //  __MOLECULE_H__
