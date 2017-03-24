#include<SDynamics.h>
#include<cmath>



const quaternion operator*(const quaternion& q1, const quaternion& q2) const
{
    
}

mat33 quaternion::toRotationMatrix() const
{
    mat33 res;

    res.a[0][0]=x*x-y*y-z*z+s*s;
    res.a[0][1]=2.0*(x*y+z*s);
    res.a[0][2]=2.0*(x*z-y*s);
    res.a[1][0]=2.0*(x*y-z*s);
    res.a[1][1]=-x*x+y*y-z*z+s*s;
    res.a[1][2]=2.0*(y*z+x*s);
    res.a[2][0]=2.0*(x*z+y*s);
    res.a[2][1]=2.0*(y*z-x*s);
    res.a[2][2]=-x*x-y*y+z*z+s*s;
    res.orthogonal=true;

    return res;
}

quaternion quaternion::conjugated() const
{
    quaternion res=quaternion(s,-x,-y,-z);
    return res;
}

mat33 mat33::transpose()
{
    mat33 res;
    res.a[0][0]=a[0][0];
    res.a[0][1]=a[1][0];
    res.a[0][2]=a[2][0];
    res.a[1][0]=a[0][1];
    res.a[1][1]=a[1][1];
    res.a[1][2]=a[2][1];
    res.a[2][0]=a[0][2];
    res.a[2][1]=a[1][2];
    res.a[2][2]=a[2][2];
    
    return res;
}

mat33 rotationMatrix(const vec3& vec, double alpha)
{
    double e1=vec.x,e2=vec.y,e3=vec.z;
    double s=pow(e1*e1+e2*e2+e3*e3,.5);
    e1/=s;
    e2/=s;
    e3/=s;
    
    double ca=cos(alpha);
    double sa=sin(alpha);
    
    mat33 res;
    
    res.a[0][0]=ca+e1*e1*(1.0-ca);
    res.a[0][1]=e1*e2*(1.0-ca)+e3*sa;
    res.a[0][2]=e1*e3*(1.0-ca)-e2*sa;
    res.a[1][0]=e1*e2*(1.0-ca)-e3*sa;
    res.a[1][1]=ca+e2*e2*(1.0-ca);
    res.a[1][2]=e2*e3*(1.0-ca)+e1*sa;
    res.a[2][0]=e1*e3*(1.0-ca)+e2*sa;
    res.a[2][1]=e2*e3*(1.0-ca)-e1*sa;
    res.a[2][2]=ca+e3*e3*(1.0-ca);
    
    return res;
}

mat33 crossProductMat(const vec3& vec)
{
    mat33 res;
    res.a[0][0]=res.a[1][1]=res.a[2][2]=.0;
    res.a[0][1]=-vec.z;
    res.a[0][2]=vec.y;
    res.a[1][0]=vec.z;
    res.a[1][2]=-vec.x;
    res.a[2][0]=-vec.y;
    res.a[2][2]=vec.x;
    return res;
}
