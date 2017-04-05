#include"transformcontroller.h"
#include"squaternion.h"

class damp
{
public:
    damp() {}
    vec3 operator()(vec3 &omega, double t)
    {
        vec3 res(-.5*omega(0),-.5*omega(1),-.5*omega(2));
        return res;
    }
};

class freeObj
{
public:
    freeObj() {};
    vec3 operator()(vec3 &omega, double t)
    {
        vec3 res;
        return res;
    }
};

TransformController::TransformController(QObject *parent)
    : QObject(parent)
    , m_target(nullptr)
    , m_quat()
    , m_time(0.0f)
{
    mat33 tensor,cosineMat;
    vec3 angularM,omega;

    double Ix=20.0,Iy=40.0,Iz=20.0;
    tensor=mat33::fromDiag(Ix,Iy,Iz);

    angularM[0]=.0;
    angularM[1]=130.0;
    angularM[2]=0.;

    vec3 vec;
    vec[0]=1.0;
    vec[1]=.0;
    vec[2]=.0;

    mat33 temp=getMatrix(vec,5.0*radPerDeg);

    tensor=temp.transpose()*tensor*temp;
    omega=tensor.inverse()*angularM;
    omega=temp*omega;
    vec3 am=tensor*(temp.transpose()*omega);
    //print(am);
    print(tensor);
    cosineMat=temp.transpose();

    m_body=rigidBody(Ix,Iy,Iz,cosineMat,omega,freeObj());
    am=m_body.getAngularMomentum();
    print(am);
    tensor=m_body.getInertiaTensor();
    print(tensor);

}

void TransformController::setTarget(Qt3DCore::QTransform *target)
{
    if (m_target != target) {
        m_target = target;
        emit targetChanged();
    }
}

Qt3DCore::QTransform *TransformController::target() const
{
    return m_target;
}

void TransformController::setTime(float time)
{
    if (!qFuzzyCompare(time, m_time)) {
        dt=(time-m_time)/12000000.0;//1000.0;
        //dt=.01;
        qDebug()<<dt<<"\t";
        m_time=time;
        updateQuat();
        emit timeChanged();
    }
}

float TransformController::time() const
{
    return m_time;
}

void TransformController::updateQuat()
{
    m_body.do_step(dt);
    m_quat=fromMat33(m_body.getCosineMat());
    m_target->setRotation(m_quat);


    static vec3 am,omega,z_;
    static vec3 z(.0,1.0,.0);
    z_=m_body.getCosineMat()*z;
    //print(z_);
    am=m_body.getAngularMomentum();
    print(omega);
    omega=m_body.getOmega();
    print(am);
    double T=m_body.getRotKineticEnergy(),theta=angle(z,z_)*degPerRad;
    qDebug()<<angle(am,omega)*degPerRad<<theta<<T;


}

