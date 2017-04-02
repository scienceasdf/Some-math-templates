#ifndef SQUATERNION_H
#define SQUATERNION_H

#include<cmath>



#include<QQuaternion>

/*#include<E:\eigen\Eigen/Dense>
#include<E:\eigen\Eigen/LU>
#include<E:\eigen\Eigen/Eigenvalues>

typedef Eigen::Matrix<double,3,1> vec3;
typedef Eigen::Matrix<double,3,3> mat33;*/

class mat33;
class vec3;

void  getAStart(double arcs[3][3],int n,double ans[3][3]);
double getA(double arcs[3][3],int n);
bool GetMatrixInverse(double src[3][3],int n,double des[3][3]);

class vec3r{
public:
    double x,y,z;
    double& operator()(int i);
    const double& operator()(int i) const;
};

class vec3{

    friend mat33 rotationMatrix(const vec3& vec, double alpha);
    friend mat33 crossProductMat(const vec3& vec);
    friend vec3 operator*(mat33 mat, vec3 vec);
    friend bool GetMatrixInverse(double src[3][3],int n,double des[3][3]);

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
    vec3r transpose() {vec3r v1; v1.x=x;v1.y=y;v1.z=z; return v1;}

    inline void setX(double xpos);
    void setY(double ypos);
    void setZ(double zpos);

    double getX() const;
    double getY() const;
    double getZ() const;

    double& operator[] (int i);
    double& operator()(int i);
    const double& operator()(int i) const;
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

const double radPerDeg=3.141592653589793/180.0;
const double degPerRad=180.0/3.141592653583297;

mat33 crossProductMat(const vec3& vec);
mat33 getMatrix(vec3 vec, double alpha);
double length(vec3& vec);

QQuaternion fromMat33(mat33& mat);

class mat33{

    //friend mat33 quaternion::toRotationMatrix();
    friend mat33 rotationMatrix(const vec3& vec, double alpha);
    friend mat33 crossProductMat(const vec3& vec);
    friend const mat33 operator*(const mat33& m1, const mat33& m2);
    friend mat33 operator*(double factor, mat33 mat);
    friend vec3 operator*(mat33 mat, vec3 vec);
    friend mat33 crossProductMat(const vec3& vec);
private:
    double a[3][3];
    bool orthogonal;
public:
    mat33() {}
    ~mat33() {}

    mat33 inverse();
    mat33 transpose();

    double& operator() (int i, int t);

};

const mat33 operator*(const mat33& m1, const mat33& m2);

class mat66{
private:
    double a[6][6];
public:
    mat66() {}
};


mat33 crossProductMat(const vec3& vec);
mat33 rotationMatrix(const vec3& vec, double alpha);

const quaternion operator*(const quaternion& q1, const quaternion& q2);
mat33 operator*(double factor, mat33 mat);
vec3 operator*(mat33 mat, vec3 vec);
vec3 operator *(double factor, vec3 vec);
vec3 operator +(vec3 veca,vec3 vecb);
double operator *(vec3r v1, vec3 v2);
mat33 operator *(vec3 v1,vec3r v2);
mat33 operator +(mat33 m1, mat33 m2);
mat33 operator -(mat33 m1, mat33 m2);

void print(vec3& vec);
void print(mat33& mat);


#endif // SQUATERNION_H
