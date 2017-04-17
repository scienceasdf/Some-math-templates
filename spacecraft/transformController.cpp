#include"transformcontroller.h"




class damp
{
public:
    damp() {}
    vec3 operator()(vec3 &omega, mat33& cosMat, double t)
    {
        vec3 res(-.5*omega(0),-.5*omega(1),-.5*omega(2));
        return res;
    }
};

class freeObj
{
public:
    freeObj() {};
    vec3 operator()(vec3 &omega, mat33& cosMat, double t)
    {
        vec3 res;
        return res;
    }
};

class disturb
{
public:
    vec3 operator()(vec3 &omega, mat33& cosMat, double t)
    {
        vec3 res(10.0,.0,.0);
        return res;
    }
};

class ANCsystem
{
public:
    vec3 operator()(vec3 &omega, mat33& cosMat, double t)
    {
        double T=(omega(2)>.04)?-150.0:.0;//:(omega(2)<-1e-3)?6.0:.0;
        double S=.0;//=(omega(0)>.04)?-150.0:.0;
        vec3 res(-.0*omega(0)+S,-.0*omega(1),-.0*omega(2)+T);
        qDebug()<<res(2);
        return res;
    }
};

class PDcontroller
{
private:
    mat33 m_target;
    QQuaternion m_q;
public:
    PDcontroller():m_target(rotationMatrix(vec3(1.0,.0,.0),1.57))
    {
        m_target=m_target.transpose();
        m_q=fromMat33(m_target);
    }

    vec3 operator()(vec3 &omega, mat33& cosMat, double t)
    {
        static QQuaternion quat;
        static QQuaternion qe;
        quat=fromMat33(cosMat);
        qe=quat.conjugate()*m_q;
        double s=qe.scalar();


        double Tcx=22.0*(qe.x()*s)-24.0*omega(0);
        double Tcy=22.0*(qe.y()*s)-24.0*omega(1);
        double Tcz=22.0*(qe.z()*s)-24.0*omega(2);

        vec3 res(Tcx,Tcy,Tcz);
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

    double Ix=100.0,Iy=40.0,Iz=100.0;
    tensor=mat33::fromDiag(Ix,Iy,Iz);

    angularM[0]=.0;
    angularM[1]=430.0;
    angularM[2]=0.;

    vec3 vec;
    vec[0]=1.0;
    vec[1]=.0;
    vec[2]=.0;

    mat33 temp=getMatrix(vec,23.0*radPerDeg);

    tensor=temp.transpose()*tensor*temp;
    omega=tensor.inverse()*angularM;
    omega=temp*omega;
    vec3 am=tensor*(temp.transpose()*omega);
    //print(am);
    //print(tensor);
    cosineMat=temp.transpose();

    //cosineMat=mat33::Identity();
    //cosineMat=rotationMatrix(vec3(1.0,.0,.0),.0).transpose();
    //omega=vec3(.1,2.0,0.01);
    //omega=vec3(.01,10.0,.10);


    //m_body=rigidBody(Ix,Iy,Iz,cosineMat,omega,ANCsystem());
    am=m_body.getAngularMomentum();
    //print(am);
    tensor=m_body.getInertiaTensor();
    //print(tensor);

}

void TransformController::setType(int type)
{
    mat33 tensor,cosineMat;
    vec3 angularM,omega;
    double Ix,Iy,Iz;

    vec3 vec;
    vec[0]=1.0;
    vec[1]=.0;
    vec[2]=.0;

    mat33 temp;

    vec3 am;

    switch (type) {
    case 0:
        //short bold rotation dynamics simulation
        Ix=100;
        Iy=60;
        Iz=100;
        tensor=mat33::fromDiag(Ix,Iy,Iz);

        temp=getMatrix(vec,20.0*radPerDeg);

        tensor=temp.transpose()*tensor*temp;

        angularM[0]=.0;
        angularM[1]=230.0;
        angularM[2]=0.;
        omega=tensor.inverse()*angularM;

        omega=temp*omega;

        cosineMat=temp.transpose();
        m_body=rigidBody(Ix,Iy,Iz,cosineMat,omega,freeObj());

        break;
    case 1:
        //slender body rotation dynamics simulation
        Ix=60;
        Iy=100;
        Iz=60;
        tensor=mat33::fromDiag(Ix,Iy,Iz);

        temp=getMatrix(vec,20.0*radPerDeg);

        tensor=temp.transpose()*tensor*temp;

        angularM[0]=.0;
        angularM[1]=230.0;
        angularM[2]=0.;
        omega=tensor.inverse()*angularM;

        omega=temp*omega;
        am=tensor*(temp.transpose()*omega);
        print(am);
        qDebug()<<"aaaaooo";

        cosineMat=temp.transpose();
        m_body=rigidBody(Ix,Iy,Iz,cosineMat,omega,freeObj());

        break;
    case 2:
        //active nutation control simulation
        Ix=200.0;
        Iy=80.0;
        Iz=200.0;
        tensor=mat33::fromDiag(Ix,Iy,Iz);

        temp=getMatrix(vec,20.0*radPerDeg);
        tensor=temp.transpose()*tensor*temp;

        angularM[0]=.0;
        angularM[1]=430.0;
        angularM[2]=0.;
        omega=tensor.inverse()*angularM;

        omega=temp*omega;
        am=tensor*(temp.transpose()*omega);

        print(tensor);
        cosineMat=temp.transpose();

        m_body=rigidBody(Ix,Iy,Iz,cosineMat,omega,ANCsystem());
        break;

    case 3:
        Ix=100.0;
        Iy=40.0;
        Iz=100.0;

        cosineMat=getMatrix(vec3(1.0,.0,.0),10.0*radPerDeg).transpose();
        omega=vec3(.1,2.0,0.01);

        m_body=rigidBody(Ix,Iy,Iz,cosineMat,omega,PDcontroller());
        break;
    case 4:
        Ix=110.0;
        Iy=100.0;
        Iz=90.0;
        tensor=mat33::fromDiag(Ix,Iy,Iz);

        temp=getMatrix(vec,7.0*radPerDeg);
        tensor=temp.transpose()*tensor*temp;

        angularM[0]=.0;
        angularM[1]=430.0;
        angularM[2]=0.;
        omega=tensor.inverse()*angularM;

        omega=temp*omega;
        am=tensor*(temp.transpose()*omega);

        print(tensor);
        cosineMat=temp.transpose();

        m_body=rigidBody(Ix,Iy,Iz,cosineMat,omega,freeObj());
        break;
    default:
        break;
    }

    am=m_body.getAngularMomentum();
    qDebug()<<"eeeuuuu";
    print(am);
    tensor=m_body.getInertiaTensor();
    //print(tensor);
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


