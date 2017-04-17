#include "scenemodifier.h"

#include <QGuiApplication>

#include <Qt3DRender/qcamera.h>
#include <Qt3DCore/qentity.h>
#include <Qt3DRender/qcameralens.h>

#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QCommandLinkButton>
#include <QtWidgets/QPushButton>
#include <QtGui/QScreen>

#include <Qt3DInput/QInputAspect>

#include <Qt3DExtras/qtorusmesh.h>
#include <Qt3DRender/qmesh.h>
#include <Qt3DRender/qtechnique.h>
#include <Qt3DRender/qmaterial.h>
#include <Qt3DRender/qeffect.h>
#include <Qt3DRender/qtexture.h>
#include <Qt3DRender/qrenderpass.h>
#include <Qt3DRender/qsceneloader.h>
#include <Qt3DRender/qpointlight.h>

#include <Qt3DCore/qtransform.h>
#include <Qt3DCore/qaspectengine.h>

#include <Qt3DRender/qrenderaspect.h>
#include <Qt3DExtras/qforwardrenderer.h>

#include <Qt3DExtras/qt3dwindow.h>
#include <Qt3DExtras/qfirstpersoncameracontroller.h>

#include<squaternion.h>

int main(int argc, char **argv)
{
    /*mat33 m1,m2;
    vec3 v1(1,1,0);
    m1=getMatrix(v1,.2);
    m2=m1*m1.inverse();
    print(m2);*/

    QApplication app(argc, argv);
    Qt3DExtras::Qt3DWindow *view = new Qt3DExtras::Qt3DWindow();
    view->defaultFrameGraph()->setClearColor(QColor(QRgb(0x4d4d4f)));
    QWidget *container = QWidget::createWindowContainer(view);
    QSize screenSize = view->screen()->size();
    container->setMinimumSize(QSize(200, 100));
    container->setMaximumSize(screenSize);

    QWidget *widget = new QWidget;
    QHBoxLayout *hLayout = new QHBoxLayout(widget);
    QVBoxLayout *vLayout = new QVBoxLayout();
    vLayout->setAlignment(Qt::AlignTop);
    hLayout->addWidget(container, 1);
    hLayout->addLayout(vLayout);

    widget->setWindowTitle(QStringLiteral("Basic shapes"));

    Qt3DInput::QInputAspect *input = new Qt3DInput::QInputAspect;
    view->registerAspect(input);

    // Root entity
    Qt3DCore::QEntity *rootEntity = new Qt3DCore::QEntity();

    // Camera
    Qt3DRender::QCamera *cameraEntity = view->camera();

    cameraEntity->lens()->setPerspectiveProjection(45.0f, 16.0f/9.0f, .1f, 1000.0f);
    //cameraEntity->lens()->setProjectionType(Qt3DRender::QCameraLens::OrthographicProjection);
    cameraEntity->setPosition(QVector3D(-5.0f, 4.0f, 20.0f));
    cameraEntity->setUpVector(QVector3D(0, 1, 0));
    cameraEntity->setViewCenter(QVector3D(0, 0, 0));

    Qt3DCore::QEntity *lightEntity = new Qt3DCore::QEntity(rootEntity);
    Qt3DRender::QPointLight *light = new Qt3DRender::QPointLight(lightEntity);
    light->setColor("white");
    light->setIntensity(1);
    lightEntity->addComponent(light);
    Qt3DCore::QTransform *lightTransform = new Qt3DCore::QTransform(lightEntity);
    lightTransform->setTranslation(cameraEntity->position());
    lightEntity->addComponent(lightTransform);

    // For camera controls
    Qt3DExtras::QFirstPersonCameraController *camController = new Qt3DExtras::QFirstPersonCameraController(rootEntity);
    camController->setCamera(cameraEntity);

    // Scenemodifier
    SceneModifier *modifier = new SceneModifier(rootEntity);

    // Set root object of the scene
    view->setRootEntity(rootEntity);

    QPushButton *spinButton1=new QPushButton(widget);
    spinButton1->setText("短粗体卫星绕最大惯量轴自旋");

    QPushButton *spinButton2=new QPushButton(widget);
    spinButton2->setText("细长体卫星绕最大惯量轴自旋");

    QPushButton *pdSimulateButton=new QPushButton(widget);
    pdSimulateButton->setText("基于四元数的PD三轴姿态控制");

    QPushButton *unstabAxisButton=new QPushButton(widget);
    unstabAxisButton->setText("不对称卫星绕中间惯量轴的不稳定自旋");

    QPushButton *ANC_button=new QPushButton(widget);
    ANC_button->setText("主动章动控制仿真");

    vLayout->addWidget(spinButton1);
    vLayout->addWidget(spinButton2);
    vLayout->addWidget(pdSimulateButton);
    vLayout->addWidget(unstabAxisButton);
    vLayout->addWidget(ANC_button);

    QObject::connect(spinButton1,&QPushButton::clicked,
                     modifier,&SceneModifier::spin1);
    QObject::connect(spinButton2,&QPushButton::clicked,
                     modifier,&SceneModifier::spin2);
    QObject::connect(pdSimulateButton,&QPushButton::clicked,
                     modifier,&SceneModifier::pd_simulate);
    QObject::connect(unstabAxisButton,&QPushButton::clicked,
                     modifier,&SceneModifier::unstabAxisSimulate);
    QObject::connect(ANC_button,&QPushButton::clicked,
                     modifier,&SceneModifier::up);

    // Show window
    widget->show();
    widget->resize(1200, 800);



    return app.exec();
}

