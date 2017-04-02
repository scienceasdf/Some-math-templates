#include"squaternion.h"
#include<QDebug>


//得到给定矩阵src的逆矩阵保存到des中。
bool GetMatrixInverse(double src[3][3],int n,double des[3][3])
{
    double flag=getA(src,n);
    double t[3][3];
    if(flag==0){
        return false;
    }
    else{
        getAStart(src,n,t);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                des[i][j]=t[i][j]/flag;
            }

        }
    }
    return true;
}

//按第一行展开计算|A|
double getA(double arcs[3][3],int n)
{
    if(n==1)
    {
        return arcs[0][0];
    }
    double ans = 0;
    double temp[3][3]={0.0};
    int i,j,k;
    for(i=0;i<n;i++){
        for(j=0;j<n-1;j++){
            for(k=0;k<n-1;k++){
                temp[j][k] = arcs[j+1][(k>=i)?k+1:k];

            }
        }
        double t = getA(temp,n-1);
        if(i%2==0){
            ans += arcs[0][i]*t;
        }
        else{
            ans -=  arcs[0][i]*t;
        }
    }
    return ans;
}

//计算每一行每一列的每个元素所对应的余子式，组成A*
void  getAStart(double arcs[3][3],int n,double ans[3][3])
{
    if(n==1){
        ans[0][0] = 1;
        return;
    }
    int i,j,k,t;
    double temp[3][3];
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            for(k=0;k<n-1;k++){
                for(t=0;t<n-1;t++){
                    temp[k][t] = arcs[k>=i?k+1:k][t>=j?t+1:t];
                }
            }

            ans[j][i]  =  getA(temp,n-1);
            if((i+j)%2 == 1){
                ans[j][i] = - ans[j][i];
            }
        }
    }
}

const quaternion operator*(const quaternion& q1, const quaternion& q2)
{
    double s,i,j,k,s1,s2,x1,x2,y1,y2,z1,z2;
    s1=q1.s; s2=q2.s;
    x1=q1.x; x2=q2.x;
    y1=q1.y; y2=q2.y;
    z1=q1.z; z2=q2.z;
    s=(s1*s2-x1*x2-y1*y2-z1*z2);
    i=(s1*x2+s2*x1+y1*z2-y2*z1);
    j=(s1*y2+s2*y1+z1*x2-x1*z2);
    k=(s1*z2+s2*z1+x1*y2-x2*y1);
    quaternion res(s,i,j,k);
    return res;
}


quaternion quaternion::conjugated() const
{
    quaternion res=quaternion(s,-x,-y,-z);
    return res;
}

void quaternion::normalize()
{
    double l=length();
    s/=l;
    x/=l;
    y/=l;
    z/=l;
}

void vec3::normalize()
{
    double l=length();
    x/=l;
    y/=l;
    z/=l;
}

mat33 crossProductMat(const vec3& vec)
{
    mat33 mat;

    mat(0,0)=.0;
    mat(1,1)=.0;
    mat(2,2)=.0;
    mat(0,1)=-vec(2);
    mat(0,2)=vec(1);
    mat(1,0)=vec(2);
    mat(1,2)=-vec(0);
    mat(2,0)=-vec(1);
    mat(2,1)=vec(0);
    return mat;
}

mat33 getMatrix(vec3 vec, double alpha)
{
    mat33 res;

    mat33 I;
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            I(i,j)=(i==j)?1.0:.0;
        }
    }

    vec.normalize();
    mat33 E=crossProductMat(vec);
    //print(E);
    //print(I);
    //std::cout<<vec.transpose()<<"\n";
    res=cos(alpha)*I+(1-cos(alpha))*vec*vec.transpose()-sin(alpha)*E;
    return res;
}

double length(vec3& vec)
{
    return pow(vec(0)*vec(0)+vec(1)*vec(1)+vec(2)*vec(2),.5);
}

QQuaternion fromMat33(mat33& mat)
{
    double q1=.5*pow(1.0+mat(0,0)-mat(1,1)-mat(2,2),.5);
    double q2=.25*(mat(0,1)+mat(1,0))/q1;
    double q3=.25*(mat(0,2)+mat(2,0))/q1;
    double q4=.25*(mat(1,2)-mat(2,1))/q1;

    QQuaternion res(q4,q1,q2,q3);
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
    //A matrix map a vector in the origin frame to the vector
    //in the body frame( or the rotated frame)

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

/*mat33 crossProductMat(const vec3& vec)
{
    mat33 res;
    res.a[0][0]=.0;
    res.a[1][1]=.0;
    res.a[2][2]=.0;
    res.a[0][1]=-vec.z;
    res.a[0][2]=vec.y;
    res.a[1][0]=vec.z;
    res.a[1][2]=-vec.x;
    res.a[2][0]=-vec.y;
    res.a[2][2]=vec.x;
    return res;
}*/

const mat33 operator*(const mat33& m1, const mat33& m2)
{
    mat33 res;
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            res.a[i][j]=.0;
            for(int k=0;k<3;++k)
                res.a[i][j]+=m1.a[i][k]*m2.a[k][j];
        }
    }

    return res;
}

mat33 operator*(double factor, mat33 mat)
{
    mat33 res;
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            res.a[i][j]=mat.a[i][j]*factor;
        }
    }

    return res;
}

vec3 operator*(mat33 mat, vec3 vec)
{
    vec3 res;

    res.x=mat.a[0][0]*vec.x+mat.a[0][1]*vec.y+mat.a[0][2]*vec.z;
    res.y=mat.a[1][0]*vec.x+mat.a[1][1]*vec.y+mat.a[1][2]*vec.z;
    res.z=mat.a[2][0]*vec.x+mat.a[2][1]*vec.y+mat.a[2][2]*vec.z;

    return res;
}

mat33 mat33::inverse()
{
    mat33 res;
    GetMatrixInverse(a,3,res.a);
    return res;
}

double& mat33::operator ()(int i,int j)
{
    return a[i][j];
}


double& vec3::operator [](int i)
{
    return (i==0)?x:((i==1)?y:z);
}

double& vec3::operator()(int i)
{
    return (i==0)?x:((i==1)?y:z);
}

const double& vec3::operator()(int i) const
{
    return (i==0)?x:((i==1)?y:z);
}

vec3 operator *(double factor, vec3 vec)
{
    vec3 res(factor*vec.getX(),factor*vec.getY(),factor*vec.getZ());
    return res;
}

vec3 operator +(vec3 veca,vec3 vecb)
{
    vec3 res(veca.getX()+vecb.getX(),veca.getY()+vecb.getY(),veca.getZ()+vecb.getZ());
    return res;
}

double& vec3r::operator()(int i)
{
    return (i==0)?x:((i==1)?y:z);
}

const double& vec3r::operator()(int i) const
{
    return (i==0)?x:((i==1)?y:z);
}

mat33 operator *(vec3 v1, vec3r v2)
{
    mat33 res;
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            res(i,j)=v1(i)*v2(j);
        }
    }
    return res;
}

double operator *(vec3r v1,vec3 v2)
{
    double res=.0;
    res=v1(0)*v2(0)+v1(1)*v2(1)+v1(2)*v2(2);
    return res;
}

mat33 operator +(mat33 m1, mat33 m2)
{
    mat33 res;
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            res(i,j)=m1(i,j)+m2(i,j);
        }
    }
    return res;
}

mat33 operator -(mat33 m1, mat33 m2)
{
    mat33 res;
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            res(i,j)=m1(i,j)-m2(i,j);
        }
    }
    return res;
}

void print(vec3& vec)
{
    qDebug() <<vec(0)<<"i+"<<vec(1)<<"j+"<<vec(2)<<"k\t";
}

void print(mat33& mat)
{
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            qDebug()<<mat(i,j)<<"\t";
        }
        qDebug()<<"\n";
    }

}
