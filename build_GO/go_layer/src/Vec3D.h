#ifndef    __VEC3D_H__
#define    __VEC3D_H__

namespace nkchem {

using namespace std;

/**
 * \brief Vector operations.
 */
class Vec3D {
public:
    /**
     * \brief Construct a zero vector \f$(0,0,0)\f$.
     */
    Vec3D(void); 
    /**
     * \brief Construct a vector from its components: `x = tx; y = ty; z = tz`.
     */
    Vec3D(double tx, double ty, double tz);
    /**
     * \brief A copy function.
     */
    Vec3D(const Vec3D& v0);
    ~Vec3D(void);
    
    /**
     * \name Overload operators.
     *
     * \{
     */    
    /**
     * \brief Vector addition.
     *
     * \return 
     * \f[
     * \mathbf{v}+\mathbf{v}_0
     * \f]
     */
    Vec3D operator + (const Vec3D& v0) const;
    /**
     * \brief Vector addition.
     *
     * \return 
     * \f[
     * \mathbf{v} = \mathbf{v}+\mathbf{v}_0
     * \f]
     */
    Vec3D& operator += (const Vec3D& v0);
    /**
     * \brief Vector substraction.
     *
     * \return 
     * \f[
     * \mathbf{v}-\mathbf{v}_0
     * \f]
     */
    Vec3D operator - (const Vec3D& v0) const;
    /**
     * \brief Vector substraction.
     *
     * \return 
     * \f[
     * \mathbf{v} = \mathbf{v}-\mathbf{v}_0
     * \f]
     */
    Vec3D& operator -= (const Vec3D& v0);
    /**
     * \brief Vector scalar multiplication.
     *
     * \return 
     * \f[
     * t\mathbf{v}
     * \f]
     */
    Vec3D operator * (double t) const;
    /**
     * \brief Vector scalar multiplication.
     *
     * \return 
     * \f[
     * \mathbf{v} = t\mathbf{v}
     * \f]
     */
    Vec3D& operator *= (double t);
    /**
     * \brief Vector scalar division.
     *
     * \return 
     * \f[
     * \mathbf{v} = \frac{1}{t}\mathbf{v}
     * \f]
     *
     * \warning For the same calculation, multiplication is preferred to division since 
     * it is more efficient.
     */
    Vec3D operator / (double t) const;
    /**
     * \brief Vector scalar division.
     *
     * \return 
     * \f[
     * \mathbf{v} = \frac{1}{t}\mathbf{v}
     * \f]
     *
     * \warning For the same calculation, multiplication is preferred to division since 
     * it is more efficient.
     */
    Vec3D& operator /= (double t);
    /**
     * \brief Vector inner product.
     *
     * \return 
     * \f[
     * \mathbf{v} \cdot \mathbf{v0}  \equiv x x_0+y y_0+z z_0
     * \f]
     *
     */
    double operator & (const Vec3D& v0) const;
    /**
     * \brief Vector outer product.
     *
     * \return 
     * \f[
     * \mathbf{v} \times \mathbf{v0}  \equiv 
     * \left\{y z_0-z y_0, z x_0-x z_0, x y_0-y x_0\right\}
     * \f]
     *
     */
    Vec3D operator * (const Vec3D& v0) const;
    /**
     * \brief Vector outer Product.
     *
     * \return 
     * \f[
     * \mathbf{v}  = \mathbf{v} \times \mathbf{v0}
     * \f]
     *
     */
    Vec3D& operator *= (const Vec3D& v0);
    /**
     *
     * \brief Assignment.
     */
    Vec3D& operator = (const Vec3D& v0);
    /**
     * \}
     */

    /**
     * \name Order.
     *
     * Two vectors, we define \f$ \mathbf{v} < \mathbf{v}_0 \f$ if any of the following
     * conditions is satisfied:
     * * if \f$x < x_0\f$;
     * * if \f$x = x_0\f$ and \f$y < y_0\f$;
     * * if \f$x = x_0\f$ and \f$y = y_0\f$ and \f$z < z_0\f$.
     *
     * Two components, e.g. \f$x\f$ and \f$x_0\f$ are considered to be euqal if \f$|x-x_0|<\f$ nkchem::Vec3D::Threshold.
     *
     * \{
     */
    bool operator == (const Vec3D& v0) const;
    bool operator != (const Vec3D& v0) const;
    bool operator <  (const Vec3D& v0) const;
    bool operator <= (const Vec3D& v0) const;
    bool operator >  (const Vec3D& v0) const;
    bool operator >= (const Vec3D& v0) const;
    /**
     * \}
     */

    /**
     * \name Norms
     */
    /**
     * \brief Compute the square of the norm.
     * 
     * \return
     * \f[
     *  x^2+y^2+z^2
     * \f]
     * 
     */
    double norm2(void) const;
    /**
     * \brief Compute the square of the norm with respect to another vector.
     * 
     * \return
     * \f[
     *  (x-x_0)^2+(y-y_0)^2+(z-z_0)^2
     * \f]
     * 
     */
    double norm2(const Vec3D& v0) const;
    /**
     * \brief Compute the norm.
     * 
     * \return
     * \f[
     *  \sqrt{x^2+y^2+z^2}
     * \f]
     * 
     */
    double norm(void) const;    
    /**
     * \brief Compute the square of the norm with respect to another vector.
     * 
     * \return
     * \f[
     *  \sqrt{(x-x_0)^2+(y-y_0)^2+(z-z_0)^2}
     * \f]
     * 
     */
    double norm(const Vec3D& v0) const;
    /**
     * \brief Compute the unit vector.
     * 
     * \return
     * \f[
     *  \frac{\mathbf{v}}{|\mathbf{v}|}
     * \f]
     * 
     */
    Vec3D normI(void) const;
    /**
     * \}
     */
    
    static const double Threshold; //!< Threshold below which two double precision numbers are equal.

    double x; //!< X component
    double y; //!< Y component
    double z; //!< Z component
};

}

#endif //  __VEC3D_H__
