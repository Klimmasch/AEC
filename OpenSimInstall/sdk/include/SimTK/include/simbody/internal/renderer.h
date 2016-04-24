#ifndef SimTK_SIMBODY_RENDERER_H
#define SimTK_SIMBODY_RENDERER_H

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "../src/VisualizerGeometry.h"

#include <cstdlib>
#include <cstdio>
#include <pthread.h>
#include <string>
#include <ctime>
#include <iostream>
#include <limits>

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include <GL/glut.h>

#include "SimTKcommon.h"
#include "gllibheader.h"
#include "renderedmesh.h"
#include "renderedline.h"
#include "mesh.h"
#include "scene.h"
//#include "lodepng.h"

using namespace std;

class SimTK_SIMBODY_EXPORT Renderer
{
public:
    Renderer(MultibodySystem& system);
    Renderer();
    void drawFrameNow(const State& state);

    void backgroundColorChange(fVec3& color);    
    void setKeepAlive(bool enable);    
    
    void forceActiveRedisplay();

    void setClearColorToBackgroundColor();

    void addMesh(vector<float> vertices, vector<float> faces);

    void addVec(vector<float>& data, float x, float y, float z);
    void addVec(vector<unsigned short>& data, int x, int y, int z);

    Mesh* makeBox();
    Mesh* makeSphere(unsigned short resolution);
    Mesh* makeCylinder(unsigned short resolution);
    Mesh* makeCircle(unsigned short resolution);

    void computeSceneBounds(Scene* scene, float& radius, fVec3& center);
    void zoomCameraToShowWholeScene();
    void setCameraTransform(fVec3 R, fVec3 p);

    void drawSkyVertex(fVec3 position, float texture);
    void drawGroundAndSky(float farClipDistance);
    void renderScene();
    void redrawDisplay();

    void saveImageTask(const string &filename, int width, int height, Array_<unsigned char> data);
    void writeImage(const string& filename);

    void addDecorationGenerator(DecorationGenerator* generator);
    void addDecoration(MobilizedBodyIndex mobodIx, const Transform& X_BD, const DecorativeGeometry& geom);
    MultibodySystem& getSystem();

    void setCameraClippingPlanes(Real nearPlane, Real farPlane);
    void setCameraFieldOfView(Real fov);

    //void setSystemUpDirection(const CoordinateDirection& upDir);

    //CoordinateDirection getSystemUpDirection();

    void setGroundHeight(Real height);

    Real getGroundHeight();



    //Variables
    // This is the transform giving the pose of the camera's local frame in the
    // model's ground frame. The camera local frame has Y as the up direction,
    // -Z as the "look at" direction, and X to the right. We can't know a good
    // default transform for the camera until we know what the SimTK::System
    // we're viewing considers to be its "up" and "look at" directions.
    fTransform X_GC;
    Scene* scene;      // This is the front scene.
    bool passiveRedisplayRequested;    
    
    int viewWidth, viewHeight;
    GLfloat fieldOfView;
    GLfloat nearClip;
    GLfloat farClip;
    GLfloat groundHeight;
    //CoordinateDirection groundNormal;
    bool showGround, showShadows;    
    fVec3 backgroundColor;    
    int nextMeshIndex;
        
    vector<vector<Mesh*> > meshes;
    

    //Initialize geometry
    MultibodySystem system;
    Array_<DecorativeGeometry>              addedGeometry;
//    Array_<RubberBandLine>                  lines;
    Array_<DecorationGenerator*>            generators;
//    Array_<Visualizer::InputListener*>      listeners;
//    Array_<Visualizer::FrameController*>    controllers;
    //CoordinateDirection                     upDirection;



};
#endif // RENDERER_H
