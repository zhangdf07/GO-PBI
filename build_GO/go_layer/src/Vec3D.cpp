#include "Vec3D.h"
#include <cmath>

using namespace nkchem;
using namespace std;

// Threshold below which two double precision numbers are equal.
const double Vec3D::Threshold = 1.00E-10;

Vec3D::Vec3D(void) : x(0.), y(0.), z(0.) 
{ 
    /* Nothing to do */ 
}

Vec3D::Vec3D(double tx, double ty, double tz) : x(tx), y(ty), z(tz)
{ 
    /* Nothing to do */ 
}

Vec3D::Vec3D(const Vec3D& v0) : x(v0.x), y(v0.y), z(v0.z) 
{ 
    /* Nothing to do */ 
}

Vec3D::~Vec3D(void)
{ 
    /* Nothing to do */ 
}

Vec3D Vec3D::operator + (const Vec3D& v0) const    
{
    return Vec3D(x+v0.x, y+v0.y, z+v0.z); 
}

Vec3D& Vec3D::operator += (const Vec3D& v0)
{ 
    x += v0.x; y += v0.y; z += v0.z;
    return *this;
}

Vec3D Vec3D::operator - (const Vec3D& v0) const
{ 
    return Vec3D(x-v0.x, y-v0.y, z-v0.z); 
}

Vec3D& Vec3D::operator -= (const Vec3D& v0)
{ 
    x -= v0.x; y -= v0.y; z -= v0.z;
    return *this;
}

Vec3D Vec3D::operator * (double t) const
{ 
    return Vec3D(x*t, y*t, z*t); 
}

Vec3D& Vec3D::operator *= (double t)
{ 
    x *= t; y *= t; z *= t;
    return *this;
}

Vec3D Vec3D::operator / (double t) const
{ 
    const double rec_t = 1./t;
    return Vec3D(x*rec_t, y*rec_t, z*rec_t); 
}

Vec3D& Vec3D::operator /= (double t)
{ 
    const double rec_t = 1./t;
    x *= rec_t; y *= rec_t; z *= rec_t;
    return *this;
}

double Vec3D::operator & (const Vec3D& v0) const
{ 
    return x*v0.x+y*v0.y+z*v0.z; 
}

Vec3D Vec3D::operator * (const Vec3D& v0) const
{ 
    return Vec3D(y*v0.z-z*v0.y, z*v0.x-x*v0.z, x*v0.y-y*v0.x); 
}

Vec3D& Vec3D::operator *= (const Vec3D& v0)
{
    const double tx = y*v0.z-z*v0.y; 
    const double ty = z*v0.x-x*v0.z; 
    const double tz = x*v0.y-y*v0.x;
    x = tx; y = ty; z = tz;
    return *this;
}

Vec3D& Vec3D::operator = (const Vec3D& v0)
{
    x = v0.x; y = v0.y; z = v0.z;
    return *this;    
}

bool Vec3D::operator == (const Vec3D& v0) const
{
    return (fabs(x-v0.x) < Threshold && fabs(y-v0.y) < Threshold &&
		fabs(z-v0.z) < Threshold);
}

bool Vec3D::operator != (const Vec3D& v0) const
{
    return !(*this == v0);
}

bool Vec3D::operator < (const Vec3D& v0) const
{
    if(x < v0.x) return true;
    if(fabs(x-v0.x) < Threshold && y < v0.y) return true;
    if(fabs(x-v0.x) < Threshold && fabs(y-v0.y) < Threshold && z < v0.z) 
	return true;
    return false;
}

bool Vec3D::operator <= (const Vec3D& v0) const
{ 
    return (*this < v0) || (*this == v0); 
}

bool Vec3D::operator > (const Vec3D& v0) const
{ 
    return !(*this <= v0); 
}

bool Vec3D::operator >= (const Vec3D& v0) const
{
    return !(*this < v0); 
}

double Vec3D::norm2(void) const
{
    return x*x+y*y+z*z; 
}

double Vec3D::norm2(const Vec3D& v0) const
{
    return (*this-v0).norm2(); 
}

double Vec3D::norm(void) const
{ 
    return sqrt(this->norm2()); 
}

double Vec3D::norm(const Vec3D& v0) const    
{
    return sqrt(this->norm2(v0)); 
}

Vec3D Vec3D::normI(void) const
{
    const double t = this->norm();
    return (*this)/t;
}
