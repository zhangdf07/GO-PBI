#include "Molecule.h"
#include "Vec3D.h"
#include "lbfgs.h"
#include <string>
#include <vector>
#include <cmath>

using namespace nkchem;
using namespace std;

Molecule::Molecule(void)
{
    /* Nothing to do. */
}

Molecule::~Molecule(void)
{
    /* Nothing to do. */
}

void Molecule::SetCoordsAndSymbols(const vector<Vec3D>& tcoords, const vector<string>& tsymbols)
{
    coords = tcoords;
    symbols = tsymbols;
}

const vector<Vec3D>& Molecule::GetCoords(void) const
{
    return coords;
}

const vector<string>& Molecule::GetSymbols(void) const
{
    return symbols;
}

vector<Vec3D> Molecule::CalcBondingCoords(int bonding_idx, int directing_idx, const Vec3D& u0, const Vec3D& un) const
{
    struct OptArgs {
        Vec3D a0;
        Vec3D a1;
        Vec3D u0;
        Vec3D u1;
    };

    auto progress = [](void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls)
    {
#if 0
        printf("%.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n", );
        printf("step: %5d energy: %15.8f gradient: %15.8f change: %15.8f\n", 
                k, fx, gnorm, step); 
#endif    
        return 0;
    };
 
    auto evaluate = [](void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, 
        const lbfgsfloatval_t step)
    {
        // Initialization.
        OptArgs* oa = (OptArgs*)instance;
        const int NumAtoms = 2;
        double iCoords[NumAtoms*3] = { oa->a0.x, oa->a0.y, oa->a0.z, oa->a1.x, oa->a1.y, oa->a1.z };
        double eCoords[NumAtoms*3] = { oa->u0.x, oa->u0.y, oa->u0.z, oa->u1.x, oa->u1.y, oa->u1.z };
        double coords[NumAtoms*3];
        double derivatives[NumAtoms*18];

        // Calculate ua0 and ua1.
        const double c1 = cos(x[0]);
        const double s1 = sin(x[0]);
        const double c2 = cos(x[1]);
        const double s2 = sin(x[1]);
        const double c3 = cos(x[2]);
        const double s3 = sin(x[2]);
        const double RM[18] = {
             c1*c2*c3-s1*s3,  s1*c2*c3+c1*s3,  s2*c3, 
            -s1*c3-c1*c2*s3,  c1*c3-s1*c2*s3, -s2*s3, 
                     -c1*s2,          -s1*s2,     c2,
                  -c1*s2*c3,       -s1*s2*c3,  c2*c3,
                   c1*s2*s3,        s1*s2*s3, -c2*s3,
                     -c1*c2,          -s1*c2,    -s2
        };
        const double RX = x[3];
        const double RY = x[4];
        const double RZ = x[5];
        for(int i = 0, i3 = 0, i18 = 0; i < NumAtoms; ++i, i3 += 3, i18 += 18)
        {
            const double x0 = iCoords[i3];
            const double y0 = iCoords[i3+1];
            const double z0 = iCoords[i3+2];
            const double t1 = RM[3]*x0+RM[4]*y0+RM[5]*z0;
            const double t2 = RM[0]*x0+RM[1]*y0+RM[2]*z0;
            // Coordinates
            coords[i3]   = t2+RX;
            coords[i3+1] = t1+RY;
            coords[i3+2] = RM[6]*x0+RM[7]*y0+RM[8]*z0+RZ;
            // Derivatives
            derivatives[i18]  = -RM[1]*x0+RM[0]*y0;
            derivatives[i18+1]  = RM[9]*x0+RM[10]*y0+RM[11]*z0;
            derivatives[i18+2]  = t1;
            derivatives[i18+3]  = 1.;
            derivatives[i18+4]  = 0.;
            derivatives[i18+5]  = 0.;
            derivatives[i18+6]  = -RM[4]*x0+RM[3]*y0;
            derivatives[i18+7]  = RM[12]*x0+RM[13]*y0+RM[14]*z0;
            derivatives[i18+8]  = -t2;
            derivatives[i18+9]  = 0.;
            derivatives[i18+10]  = 1.;
            derivatives[i18+11]  = 0.;        
            derivatives[i18+12]  = -RM[7]*x0+RM[6]*y0;
            derivatives[i18+13]  = RM[15]*x0+RM[16]*y0+RM[17]*z0;
            derivatives[i18+14]  = 0.;
            derivatives[i18+15]  = 0.;
            derivatives[i18+16]  = 0.;
            derivatives[i18+17]  = 1.;
        }

        ///////////////////////////////////////////
        double delta = 0.;
        for(int i = 0; i < n; ++i) { g[i] = 0.; }
        for(int i = 0, i3 = 0, i18 = 0; i < NumAtoms; ++i, i3 += 3, i18 += 18)
        {
            const double dx = coords[i3]-eCoords[i3];
            const double dy = coords[i3+1]-eCoords[i3+1];
            const double dz = coords[i3+2]-eCoords[i3+2];
            delta += dx*dx+dy*dy+dz*dz;
            for(int j = 0; j < 6; ++j) { g[j] += dx*derivatives[i18+j]; }
            for(int j = 0; j < 6; ++j) { g[j] += dy*derivatives[i18+6+j]; }
            for(int j = 0; j < 6; ++j) { g[j] += dz*derivatives[i18+12+j]; }
        }
        for(int i = 0; i < n; ++i) { g[i] *= 2; }
        return delta;
    };

    const Vec3D a0 = coords[bonding_idx];
    const Vec3D an = (directing_idx == bonding_idx) ? Vec3D(0., 0., 1.) : (coords[directing_idx]-coords[bonding_idx]);

    double delta;
    const int dim = 6;
    double eCoords[dim];
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    // Search solutions.
    OptArgs args = {a0, a0+an.normI(), u0, u0+un.normI()};
    const int ret = lbfgs(dim, eCoords, &delta, evaluate, progress, (void*)(&args), &param);

    // Rotate.
    const double c1 = cos(eCoords[0]);
    const double s1 = sin(eCoords[0]);
    const double c2 = cos(eCoords[1]);
    const double s2 = sin(eCoords[1]);
    const double c3 = cos(eCoords[2]);
    const double s3 = sin(eCoords[2]);
    const double RM[18] = {
         c1*c2*c3-s1*s3,  s1*c2*c3+c1*s3,  s2*c3, 
        -s1*c3-c1*c2*s3,  c1*c3-s1*c2*s3, -s2*s3, 
                 -c1*s2,          -s1*s2,     c2,
              -c1*s2*c3,       -s1*s2*c3,  c2*c3,
               c1*s2*s3,        s1*s2*s3, -c2*s3,
                 -c1*c2,          -s1*c2,    -s2
    };
    const double RX = eCoords[3];
    const double RY = eCoords[4];
    const double RZ = eCoords[5];
    vector<Vec3D> bonding_coords = coords;
    const int num_atoms = bonding_coords.size();
    for(int i = 0; i < num_atoms; ++i)
    {
        const double x0 = bonding_coords[i].x;
        const double y0 = bonding_coords[i].y;
        const double z0 = bonding_coords[i].z;
        const double t1 = RM[3]*x0+RM[4]*y0+RM[5]*z0;
        const double t2 = RM[0]*x0+RM[1]*y0+RM[2]*z0;
        // Coordinates
        bonding_coords[i].x = t2+RX;
        bonding_coords[i].y = t1+RY;
        bonding_coords[i].z = RM[6]*x0+RM[7]*y0+RM[8]*z0+RZ;
    }    
    return bonding_coords;
}

void Molecule::RotZ(double phi)
{
    const double sp = sin(phi);
    const double cp = cos(phi);
    const int num_atoms = coords.size();
    for(int i = 0; i < num_atoms; ++i)
    {        
        const double x = coords[i].x;
        const double y = coords[i].y;
        const double z = coords[i].z;
        coords[i] = Vec3D(cp*x+sp*y, -sp*x+cp*y, z);
    }
}

