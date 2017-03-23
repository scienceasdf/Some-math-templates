#ifndef SDYNAMICS_H
#define SDYNAMICS_H

class vec3{
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

    void setX(double xpos);
    void setY(double ypos);
    void setZ(double zpos);

    double x() const;
    double y() const;
    double z() const;
};

inline double vec3::length() const {return pow(x*x+y*y+z*z,.5);}
inline double vec3::lengthSquared() const {return (x*x+y*y+z*z);}
inline double vec3::x() const {return x;}
inline double vec3::y() const {return y;}
inline double vec3::z() const {return z;}
inline void vec3::setX(double xpos){x=xpos;}
inline void vec3::setY(double ypos){y=ypos;}
inline void vec3::setZ(double zpos){z=zpos;}

class quaternion{
    friend quaternion operator*(const quaternion& Q) const;
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

    double scalar() const;
    double x() const;
    double y() const;
    double z() const;

    void setScalar(double scalar);
    void setX(double xpos);
    void setY(double ypos);
    void setZ(double zpos);

    vec3 toEulerAngles() const;
    mat33 toRotationMatrix() const;


};

inline double quaternion::lengthSquared() const {return (x*x+y*y+z*z+s*s);}
inline double quaternion::length() const {return pow(x*x+y*y+z*z+s*s,.5);}
inline double quaternion::x() const {return x;}
inline double quaternion::y() const {return y;}
inline double quaternion::z() const {return z;}
inline double quaternion::scalar() const {return s;}
inline void quaternion::setX(double xpos) {x=xpos;}
inline void quaternion::setY(double ypos) {y=ypos;}
inline void quaternion::setZ(double zpos) {z=zpos;}
inline void quaternion::setScalar(double scalar) {s=scalar;}

class mat33{
    friend mat33 quaternion::toRotationMatrix() const;

private:
    double a[3][3];
public:
    mat33() {}

};

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

#endif // SDYNAMICS_H
