#ifndef RENDEREDLINE_H
#define RENDEREDLINE_H

#include "gllibheader.h"

class RenderedLine
{
public:
    RenderedLine(const fVec3& color, float thickness);
    void draw(bool setColor = true);
    vector<GLfloat>& getLines();
    fVec3& getColor();
    float getThickness();
    void computeBoundingSphere(float& radius, fVec3& center);
private:
    fVec3 color;
    float thickness;
    vector<GLfloat> lines;
};

#endif // RENDEREDLINE_H
