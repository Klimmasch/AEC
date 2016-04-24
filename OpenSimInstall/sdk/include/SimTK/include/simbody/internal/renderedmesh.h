#ifndef RENDEREDMESH_H
#define RENDEREDMESH_H

#include "gllibheader.h"
#include "mesh.h"

class RenderedMesh
{
public:
    RenderedMesh(const fTransform& transform, const fVec3& scale, const fVec4& color, short representation, unsigned short meshIndex, unsigned short resolution, vector<vector<Mesh*> > meshes);
    void draw(bool setColor);
    fTransform& getTransform();
    void computeBoundingSphere(float& radius, fVec3& center);
private:
    fTransform transform;
    fVec3 scale;
    GLfloat color[4];
    short representation;
    unsigned short meshIndex, resolution;
    vector<vector<Mesh*> > meshes;
};

#endif // RENDEREDMESH_H
