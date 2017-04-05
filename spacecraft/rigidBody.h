#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include<functional>

#include"squaternion.h"



class rigidBody
{
public:
    rigidBody();
    rigidBody(double Ix, double Iy, double Iz, mat33& cosineMat, vec3& omega,std::function<vec3(vec3&,double)> mome):
        m_Ix(Ix),m_Iy(Iy),m_Iz(Iz), m_inertia(mat33::fromDiag(Ix,Iy,Iz)),
        m_cosMat(cosineMat), m_omega(omega),m_time(.0), moment(mome) {}

    void do_step(double dt);    //calculate the state after time dt

    double getRotKineticEnergy();
    vec3 getAngularMomentum();     //expressed in the inertial frame
    vec3 getOmega();            //expressed in the inertial frame
    mat33 getInertiaTensor();   //expressed in the inertial frame
    mat33 getCosineMat();



private:
    //double m_mass;
    double m_Ix,m_Iy,m_Iz;
    mat33 m_inertia;    //expressed in the body frame
    mat33 m_cosMat;
    //vec3 m_pos;
    vec3 m_omega;       //expressed in the body frame
    //vec3 m_speed;
    double m_time;
    std::function<vec3(vec3&,double)> moment; //suppose M=M(omega,t), expressed in the body frame
    //vec3 force(parameters);

    void func(double t, vec3& omega, mat33& cosMat, vec3& resVec, mat33& resMat);

};

#endif // RIGIDBODY_H
