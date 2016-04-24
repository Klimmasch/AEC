#ifndef MESH_H
#define MESH_H

#include "gllibheader.h"

class Mesh
{
public:
    Mesh(vector<float>& vertices, vector<float>& normals, vector<GLushort>& faces);
    void draw(short representation, int meshIndex);
    void getBoundingSphere(float& radius, fVec3& center);
private:
    int numVertices;
    GLuint vertBuffer, normBuffer, texBuffer;
    vector<GLushort> edges, faces;
    fVec3 center;
    float radius;
};

#endif // MESH_H

