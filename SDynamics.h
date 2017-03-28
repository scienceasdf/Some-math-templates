#ifndef SDYNAMICS_H
#define SDYNAMICS_H

#include<cmath>

class mat33;

class vec3{

    friend mat33 rotationMatrix(const vec3& vec, double alpha);
    friend mat33 crossProductMat(const vec3& vec);

private:
    double x,y,z;


public:
    vec3() :x(.0), y(.0), z(.0) {}
    vec3(double xpos, double ypos, double zpos)
        :x(xpos), y(ypos), z(zpos) {}

    ~vec3() {}

    double length() const;
    double lengthSquared() const;
    void normalize();
    vec3 normalized() const;

    inline void setX(double xpos);
    void setY(double ypos);
    void setZ(double zpos);

    double getX() const;
    double getY() const;
    double getZ() const;
};

inline double vec3::length() const {return pow(x*x+y*y+z*z,.5);}
inline double vec3::lengthSquared() const {return (x*x+y*y+z*z);}
inline double vec3::getX() const {return x;}
inline double vec3::getY() const {return y;}
inline double vec3::getZ() const {return z;}
inline void vec3::setX(double xpos){x=xpos;}
inline void vec3::setY(double ypos){y=ypos;}
inline void vec3::setZ(double zpos){z=zpos;}

class quaternion{

    //friend const quaternion operator*(const quaternion& q1, const quaternion& q2) const;
    friend const quaternion operator*(const quaternion& q1, const quaternion& q2);
private:
    double x,y,z,s;
public:
    quaternion(double scalar, double xpos, double ypos, double zpos)
        :s(scalar), x(xpos), y(ypos), z(zpos) {}
    quaternion(double scalar, const vec3& vector);
    quaternion(): s(.0),x(.0),y(.0),z(.0) {}
    ~quaternion() {}

    double length() const;
    double lengthSquared() const;

    void normalize();
    quaternion normalized() const;
    quaternion conjugated() const;

    double getScalar() const;
    double getX() const;
    double getY() const;
    double getZ() const;

    void setScalar(double scalar);
    void setX(double xpos);
    void setY(double ypos);
    void setZ(double zpos);

    vec3 toEulerAngles() const;
    mat33 toRotationMatrix() const;


};

inline double quaternion::lengthSquared() const {return (x*x+y*y+z*z+s*s);}
inline double quaternion::length() const {return pow(x*x+y*y+z*z+s*s,.5);}
inline double quaternion::getX() const {return x;}
inline double quaternion::getY() const {return y;}
inline double quaternion::getZ() const {return z;}
inline double quaternion::getScalar() const {return s;}
inline void quaternion::setX(double xpos) {x=xpos;}
inline void quaternion::setY(double ypos) {y=ypos;}
inline void quaternion::setZ(double zpos) {z=zpos;}
inline void quaternion::setScalar(double scalar) {s=scalar;}

class mat33{

    friend mat33 quaternion::toRotationMatrix() const;
    friend mat33 rotationMatrix(const vec3& vec, double alpha);
    friend mat33 crossProductMat(const vec3& vec);
    friend const mat33 operator*(const mat33& m1, const mat33& m2);
private:
    double a[3][3];
    bool orthogonal;
public:
    mat33() {}
    ~mat33() {}

    mat33 inverse();
    mat33 transpose();

};

const mat33 operator*(const mat33& m1, const mat33& m2);

class mat66{
private:
    double a[6][6];
public:
    mat66() {}
};

class rigidBody
{
    mat66 massTensor;
};

mat33 crossProductMat(const vec3& vec);
mat33 rotationMatrix(const vec3& vec, double alpha);

const quaternion operator*(const quaternion& q1, const quaternion& q2);

#endif // SDYNAMICS_H
