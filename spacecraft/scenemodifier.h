


#ifndef SCENEMODIFIER_H
#define SCENEMODIFIER_H

#include <QtCore/QObject>

#include <Qt3DCore/qentity.h>
#include <Qt3DCore/qtransform.h>

#include <Qt3DExtras/QTorusMesh>
#include <Qt3DExtras/QConeMesh>
#include <Qt3DExtras/QCylinderMesh>
#include <Qt3DExtras/QCuboidMesh>
#include <Qt3DExtras/QPlaneMesh>
#include <Qt3DExtras/QSphereMesh>
#include <Qt3DExtras/QPhongMaterial>
#include<Qt3DExtras/QGoochMaterial>

class SceneModifier : public QObject
{
    Q_OBJECT

public:
    explicit SceneModifier(Qt3DCore::QEntity *rootEntity);
    ~SceneModifier();

    void up_date();

public slots:
    void spin1();
    void spin2();
    void pd_simulate();
    void unstabAxisSimulate();
    void up();

private:
    Qt3DCore::QEntity *m_rootEntity;
    Qt3DExtras::QTorusMesh *m_torus;
    Qt3DCore::QEntity *m_coneEntity;
    Qt3DCore::QEntity *m_satEntity;
    Qt3DCore::QEntity *m_axisEntity;
    Qt3DCore::QEntity *m_cuboidEntity;
    Qt3DCore::QEntity *m_planeEntity;
    Qt3DCore::QEntity *m_sphereEntity;

    void setSome(int type);
};

#endif // SCENEMODIFIER_H


