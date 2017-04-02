#ifndef TRANSFORMCONTROLLER_H
#define TRANSFORMCONTROLLER_H

#include<QObject>
#include<Qt3DCore/qtransform.h>
#include<QQuaternion>

#include"squaternion.h"

class TransformController: public QObject
{
    Q_OBJECT
    Q_PROPERTY(Qt3DCore::QTransform* target READ target WRITE setTarget NOTIFY targetChanged)
    Q_PROPERTY(float time READ time WRITE setTime NOTIFY timeChanged)

public:
    TransformController(QObject *parent = 0);

    void setTarget(Qt3DCore::QTransform *target);
    Qt3DCore::QTransform *target() const;

    void setTime(float time);
    float time() const;


signals:
    void targetChanged();
    void timeChanged();

protected:
    void updateQuat();

private:
    Qt3DCore::QTransform *m_target;
    QQuaternion m_quat;
    float m_time;
    float dt;

    mat33 tensor;
    vec3 omega;
    vec3 angularM;
    mat33 cosineMat;

};

#endif // TRANSFORMCONTROLLER_H
