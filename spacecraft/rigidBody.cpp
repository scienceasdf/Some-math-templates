#include "rigidbody.h"


rigidBody::rigidBody()
{

}

/*vec3 rigidBody::moment(vec3 &omega, double t)
{
    //vec3 res;   //all zeroes
    //return res;
    vec3 res(-.5*omega(0),-.5*omega(1),-.5*omega(2));
    return res;

}*/

void rigidBody::func(double t, vec3& vec, mat33& cosMat, vec3& resVec, mat33& resMat)
{


    vec3 M=moment(vec,t);
    double dwx=M.getX()/m_Ix-vec.getY()*vec.getZ()/m_Ix*(m_Iz-m_Iy);
    double dwy=M.getY()/m_Iy-vec.getX()*vec.getZ()/m_Iy*(m_Ix-m_Iz);
    double dwz=M.getZ()/m_Iz-vec.getX()*vec.getY()/m_Iz*(m_Iy-m_Ix);
    resVec=vec3(dwx,dwy,dwz);
    //resMat=zero-crossProductMat(vec)*cosMat;
    resMat=cosMat*crossProductMat3(vec);



}

void rigidBody::do_step(double dt)
{
    //using RK4 algorithm

    static vec3 vk1,vk2,vk3,vk4;
    static mat33 mk1,mk2,mk3,mk4;
    static vec3 resV;
    static mat33 resM;

    double t=m_time;

    func(t,m_omega,m_cosMat,vk1,mk1);
    mk1=dt*mk1;
    vk1=dt*vk1;
    resV=m_omega+.5*vk1;
    resM=m_cosMat+.5*mk1;

    func(t+.5*dt,resV,resM,vk2,mk2);
    mk2=dt*mk2;
    vk2=dt*vk2;
    resV=m_omega+.5*vk2;
    resM=m_cosMat+.5*mk2;

    func(t+.5*dt,resV,resM,vk3,mk3);
    mk3=dt*mk3;
    vk3=dt*vk3;
    resV=m_omega+vk3;
    resM=m_cosMat+mk3;

    func(t+dt,resV,resM,vk4,mk4);
    mk4=dt*mk4;
    vk4=dt*vk4;

    m_time+=dt;
    m_omega=m_omega+.16666666666666666666666666666666666666666666666666666666*(vk1+2.0*vk2+2.0*vk3+vk4);
    m_cosMat=m_cosMat+.16666666666666666666666666666666666666666666666666666666666*(mk1+2.0*mk2+2.0*mk3+mk4);

}

double rigidBody::getRotKineticEnergy()
{
    return m_Ix*m_omega(0)*m_omega(0)+m_Iy*m_omega(1)*m_omega(1)+m_Iz*m_omega(2)*m_omega(2);
}

vec3 rigidBody::getAngularMomentum()
{
    return m_cosMat*(m_inertia*m_omega);
}

vec3 rigidBody::getOmega()
{
    return m_cosMat*m_omega;
}

mat33 rigidBody::getInertiaTensor()
{
    return m_cosMat*m_inertia*m_cosMat.transpose();
}

mat33 rigidBody::getCosineMat()
{
    return m_cosMat;
}
