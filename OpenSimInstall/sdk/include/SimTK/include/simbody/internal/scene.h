#ifndef SCENE_H
#define SCENE_H

#include "renderedmesh.h"
#include "renderedline.h"

class Scene {
public:
    Scene();

    float simTime; // simulated time associated with this frame

    vector<RenderedMesh> drawnMeshes;
    vector<RenderedMesh> solidMeshes;
    vector<RenderedMesh> transparentMeshes;
    vector<RenderedLine> lines;
    bool sceneHasBeenDrawn;
};

#endif // SCENE_H
