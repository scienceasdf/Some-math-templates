#include"transformcontroller.h"
#include"squaternion.h"

TransformController::TransformController(QObject *parent)
    : QObject(parent)
    , m_target(nullptr)
    , m_quat()
    , m_time(0.0f)
{
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            tensor(i,j)=.0;
        }
    }
    tensor(0,0)=80.0;
    tensor(1,1)=12.0;
    tensor(2,2)=80.0;

    angularM[0]=.0;
    angularM[1]=30.0;
    angularM[2]=0.;

    vec3 vec;
    vec[0]=1.0;
    vec[1]=.0;
    vec[2]=.0;

    mat33 temp=getMatrix(vec,20*radPerDeg);

    tensor=temp.transpose()*tensor*temp;
    omega=tensor.inverse()*angularM;
    print(omega);
    cosineMat=temp;

    double T=omega.transpose()*(tensor*omega);
    qDebug()<<T<<"\n";

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
        dt=(time-m_time)/5000000.0;//1000.0;
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
    //m_quat=QQuaternion::fromAxisAndAngle(1,2,3,m_time);
    static mat33 tensor_;
    static vec3 ome;
    static vec3 omega_;
    static mat33 A;
    const static vec3 v1(.0,.0,1.0);
    static vec3 v2;

    ome=omega;

    A=getMatrix(ome,length(ome)*dt);
    tensor_=A.transpose()*tensor*A;
    omega_=tensor_.inverse()*angularM;
    ome=.5*ome+.5*omega_;

    A=getMatrix(ome,length(ome)*dt);
    tensor_=A.transpose()*tensor*A;
    omega_=tensor_.inverse()*angularM;
    ome=.5*ome+.5*omega_;

    A=getMatrix(ome,length(ome)*dt);
    tensor=A.transpose()*tensor*A;
    omega=tensor.inverse()*angularM;
    cosineMat=A*cosineMat;
    v2=cosineMat.transpose()*v1;


    print(omega);
    print(v2);
    double T=omega.transpose()*(tensor*omega);
    //double theta=std::atan2(omega.getY(),pow(omega.getX()*omega.getX()+omega.getZ()*omega.getZ(),.5));
    double theta=std::atan2(v2.getY(),pow(v2.getX()*v2.getX()+v2.getZ()*v2.getZ(),.5));
    qDebug()<<theta<<"\t"<<T<<"\n";

    m_quat=fromMat33(cosineMat);
    m_target->setRotation(m_quat);
}

