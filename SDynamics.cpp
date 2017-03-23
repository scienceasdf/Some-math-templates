#include<SDynamics.h>
#include<cmath>



quaternion operator*(const quaternion& Q)
{
    quaternion res;
    res.x=

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

    return res;
}

quaternion quaternion::conjugated() const
{
    quaternion res=quaternion(s,-x,-y,-z);
    return res;
}


