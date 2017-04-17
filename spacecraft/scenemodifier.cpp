

#include "scenemodifier.h"

#include <QtCore/QDebug>
#include<QPropertyAnimation>
#include<windows.h>

#include"transformcontroller.h"


SceneModifier::SceneModifier(Qt3DCore::QEntity *rootEntity)
    : m_rootEntity(rootEntity)
{

    Qt3DExtras::QCylinderMesh *c1 = new Qt3DExtras::QCylinderMesh();
    c1->setRadius(.1);
    c1->setLength(3);
    c1->setRings(100);
    c1->setSlices(20);
    Qt3DCore::QTransform *c1Transform = new Qt3DCore::QTransform();
    c1Transform->setScale(1.5f);
    c1Transform->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), .0f));
    c1Transform->setTranslation(QVector3D(-5.0f, 4.0f, -1.5));


    //! [2]
    Qt3DExtras::QPhongMaterial *torusMaterial = new Qt3DExtras::QPhongMaterial();
    torusMaterial->setDiffuse(QColor(QRgb(0xbeb32b)));
    //! [2]

    // Torus
    //! [3]
    m_axisEntity = new Qt3DCore::QEntity(m_rootEntity);
    m_axisEntity->addComponent(c1);
    m_axisEntity->addComponent(torusMaterial);
    m_axisEntity->addComponent(c1Transform);
    //! [3]

    // Cone shape data
    Qt3DExtras::QConeMesh *cone = new Qt3DExtras::QConeMesh();
    cone->setTopRadius(0.5);
    cone->setBottomRadius(1);
    cone->setLength(3);
    cone->setRings(50);
    cone->setSlices(20);

    // ConeMesh Transform
    Qt3DCore::QTransform *coneTransform = new Qt3DCore::QTransform();
    coneTransform->setScale(1.5f);
    coneTransform->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), 45.0f));
    coneTransform->setTranslation(QVector3D(0.0f, 4.0f, -1.5));

    Qt3DExtras::QPhongMaterial *coneMaterial = new Qt3DExtras::QPhongMaterial();
    coneMaterial->setDiffuse(QColor(QRgb(0x928327)));

    // Cone
    m_coneEntity = new Qt3DCore::QEntity(m_rootEntity);
    m_coneEntity->addComponent(cone);
    m_coneEntity->addComponent(coneMaterial);
    m_coneEntity->addComponent(coneTransform);

    // Cylinder shape data
    Qt3DExtras::QCylinderMesh *cylinder = new Qt3DExtras::QCylinderMesh();
    cylinder->setRadius(1);
    cylinder->setLength(3);
    cylinder->setRings(100);
    cylinder->setSlices(20);

    // CylinderMesh Transform
    Qt3DCore::QTransform *cylinderTransform = new Qt3DCore::QTransform();
    cylinderTransform->setScale(1.5f);
    cylinderTransform->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), 20.0f));
    cylinderTransform->setTranslation(QVector3D(-5.0f, 4.0f, -1.5));

    Qt3DExtras::QPhongMaterial *cylinderMaterial = new Qt3DExtras::QPhongMaterial();
    cylinderMaterial->setDiffuse(QColor(QRgb(0x928327)));

    // Cylinder
    m_satEntity = new Qt3DCore::QEntity(m_rootEntity);
    m_satEntity->addComponent(cylinder);
    m_satEntity->addComponent(cylinderMaterial);
    m_satEntity->addComponent(cylinderTransform);

    // Cuboid shape data
    Qt3DExtras::QCuboidMesh *cuboid = new Qt3DExtras::QCuboidMesh();
    cuboid->setXExtent(1.1);

    // CuboidMesh Transform
    Qt3DCore::QTransform *cuboidTransform = new Qt3DCore::QTransform();
    cuboidTransform->setScale(4.0f);
    cuboidTransform->setTranslation(QVector3D(5.0f, -4.0f, 0.0f));

    Qt3DExtras::QPhongMaterial *cuboidMaterial = new Qt3DExtras::QPhongMaterial();
    cuboidMaterial->setDiffuse(QColor(QRgb(0x665423)));

    //Cuboid
    m_cuboidEntity = new Qt3DCore::QEntity(m_rootEntity);
    m_cuboidEntity->addComponent(cuboid);
    m_cuboidEntity->addComponent(cuboidMaterial);
    m_cuboidEntity->addComponent(cuboidTransform);

    // Plane shape data
    Qt3DExtras::QPlaneMesh *planeMesh = new Qt3DExtras::QPlaneMesh();
    planeMesh->setWidth(2);
    planeMesh->setHeight(2);

    // Plane mesh transform
    Qt3DCore::QTransform *planeTransform = new Qt3DCore::QTransform();
    planeTransform->setScale(1.3f);
    planeTransform->setRotation(QQuaternion::fromAxisAndAngle(QVector3D(1.0f, 0.0f, 0.0f), 45.0f));
    planeTransform->setTranslation(QVector3D(0.0f, -4.0f, 0.0f));

    Qt3DExtras::QPhongMaterial *planeMaterial = new Qt3DExtras::QPhongMaterial();
    planeMaterial->setDiffuse(QColor(QRgb(0xa69929)));

    // Plane
    m_planeEntity = new Qt3DCore::QEntity(m_rootEntity);
    m_planeEntity->addComponent(planeMesh);
    m_planeEntity->addComponent(planeMaterial);
    m_planeEntity->addComponent(planeTransform);

    // Sphere shape data
    Qt3DExtras::QSphereMesh *sphereMesh = new Qt3DExtras::QSphereMesh();
    sphereMesh->setRings(20);
    sphereMesh->setSlices(20);
    sphereMesh->setRadius(2);

    // Sphere mesh transform
    Qt3DCore::QTransform *sphereTransform = new Qt3DCore::QTransform();

    sphereTransform->setScale(1.3f);
    sphereTransform->setTranslation(QVector3D(-5.0f, -4.0f, 0.0f));

    Qt3DExtras::QPhongMaterial *sphereMaterial = new Qt3DExtras::QPhongMaterial();
    sphereMaterial->setDiffuse(QColor(QRgb(0xa69929)));

    // Sphere
    m_sphereEntity = new Qt3DCore::QEntity(m_rootEntity);
    m_sphereEntity->addComponent(sphereMesh);
    m_sphereEntity->addComponent(sphereMaterial);
    m_sphereEntity->addComponent(sphereTransform);

}

SceneModifier::~SceneModifier()
{
}

void SceneModifier::setSome(int type)
{
    delete m_satEntity;
    m_satEntity = new Qt3DCore::QEntity(m_rootEntity);

    Qt3DExtras::QPhongMaterial *satMaterial=new Qt3DExtras::QPhongMaterial();
    satMaterial->setDiffuse(QColor(QRgb(0x928327)));



    // satMesh Transform
    Qt3DCore::QTransform *satTransform = new Qt3DCore::QTransform();
    satTransform->setScale(1.5f);

    satTransform->setTranslation(QVector3D(-5.0f, 4.0f, -1.5));
    TransformController *ctrl=new TransformController(satTransform);
    ctrl->setType(type);
    ctrl->setTarget(satTransform);
    ctrl->setTime(.0f);

    QPropertyAnimation *animation=new QPropertyAnimation(satTransform);
    animation->setTargetObject(ctrl);
    animation->setPropertyName("time");
    animation->setDuration(100000);
    animation->setStartValue(QVariant::fromValue(.0));
    animation->setEndValue(QVariant::fromValue(800000000));
    animation->setLoopCount(1);
    animation->start();

    // sat
    m_satEntity->addComponent(satTransform);
    m_satEntity->addComponent(satMaterial);

}

void SceneModifier::up_date()
{
    setSome(2);

    Qt3DExtras::QCuboidMesh *sat = new Qt3DExtras::QCuboidMesh();
    sat->setYExtent(.6);

    m_satEntity->addComponent(sat);

}

void SceneModifier::spin1()
{
    setSome(1);

    Qt3DExtras::QCuboidMesh *sat = new Qt3DExtras::QCuboidMesh();
    sat->setYExtent(.6);

    m_satEntity->addComponent(sat);

}

void SceneModifier::spin2()
{
    setSome(0);

    Qt3DExtras::QCuboidMesh *sat = new Qt3DExtras::QCuboidMesh();
    sat->setYExtent(1.2);

    m_satEntity->addComponent(sat);

}

void SceneModifier::pd_simulate()
{
    setSome(3);

    Qt3DExtras::QCuboidMesh *sat = new Qt3DExtras::QCuboidMesh();
    sat->setYExtent(.6);

    m_satEntity->addComponent(sat);
}

void SceneModifier::unstabAxisSimulate()
{
    setSome(4);

    Qt3DExtras::QCuboidMesh *sat = new Qt3DExtras::QCuboidMesh();
    sat->setYExtent(.6);

    m_satEntity->addComponent(sat);

}

void SceneModifier::up()
{
    up_date();
}



