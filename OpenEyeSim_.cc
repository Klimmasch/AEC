/** Demonstration of mexplus library.
 *
 * In this example, we create MEX APIs for the hypothetical Database class in
 * Matlab.
 *
 */
#include <mexplus.h>
#include "math.h"
#include <OpenSim/OpenSim.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <cstring>
#include <cstdio>
#include <vector>

#include <GL/glut.h>

using namespace OpenSim;
using namespace SimTK;
using namespace std;
using namespace mexplus;

static GLfloat fieldOfView = 0.8726664626;
static GLfloat nearClip = 1;
static GLfloat farClip = 10;
static GLfloat groundHeight = 0;
static CoordinateDirection groundNormal = YAxis; // the +Y direction
static const int DefaultWindowWidth  = 320;
static const int DefaultWindowHeight = 240;

GLfloat light_diffuse[]  = {0.8f, 0.8f, 0.8f, 1.0f};
GLfloat light_position[] = {1.0f, 1000.0f, -100.0f, 0.0f};
GLfloat light_ambient[]  = {0.2f, 0.2f, 0.2f, 1.0f};

GLfloat n[6][3] = { /* Normals for the 6 faces of a cube. */
    {-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0},
    {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0} };
GLint faces[6][4] = { /* Vertex indices for the 6 faces of a cube. */
    {0, 1, 2, 3}, {3, 2, 6, 7}, {7, 6, 5, 4},
    {4, 5, 1, 0}, {5, 6, 2, 1}, {7, 4, 0, 3} };
GLfloat v[8][3]; /* Will be filled in with X,Y,Z vertexes. */
static unsigned int widthTex, heightTex;
static int width = ((DefaultWindowWidth + 3) / 4) * 4; // must be a multiple of 4 pixels
static int height = DefaultWindowHeight;
static int i = 0;
static bool initialized = false;
static GLuint skyTexture;
static GLuint groundTexture;
static GLuint textureID;

// Hypothetical database class to be MEXed. This example is a proxy to C++ map.
class OpenEyeSim
{
public:
    // Database constructor. This is a stub.
    OpenEyeSim()
    {
        mexPrintf("Simulator is running");
        init = false;
        angle = 0.;
        strangle = 0.;
        texture = "1.bmp";
        distance = 0.f;
        planeScale = 1.4f;
    }

    // Database destructor.
    virtual ~OpenEyeSim() {}

    // Set texture,angle and distance
    void
    set_params(int texture_number, double angle_input, float distance_input, double strabismusAngle, float planeScaling)
    {
        this->angle = angle_input;
        this->distance = distance_input;
        this->texture_number = texture_number;
        this->strangle = strabismusAngle;
        this->planeScale = planeScaling;
    }

    void
    add_texture(int number, const string& texture_input)
    {
        loadBMPRaw(texture_input.c_str(), widthTex, heightTex, number);
    }

    void
    loadBMPRaw(const char * imagepath, unsigned int& outWidth, unsigned int& outHeight, int number)
    {
        bool flipY = false;
        // Data read from the header of the BMP file
        unsigned char header[54];
        unsigned int dataPos;
        unsigned int imageSize;
        // Actual RGB data
        unsigned char * data;

        // Open the file
        FILE * file = fopen(imagepath, "rb");
        if (!file)
        {
            printf("Image could not be opened\n");
            //return NULL;
        }

        // Read the header, i.e. the 54 first bytes

        // If less than 54 byes are read, problem
        if (fread(header, 1, 54, file) != 54)
        {
            printf("Not a correct BMP file\n");
            // return NULL;
        }


        // A BMP files always begins with "BM"
        if ((header[0] != 'B') || (header[1] != 'M'))
        {
            printf("Not a correct BMP file\n");
            // return NULL;
        }


        // Make sure this is a 24bpp file
        if (*(int*)&(header[0x1E]) != 0)
        {
            printf("Not a correct BMP file\n");
            //return NULL;
        }
        if (*(int*)&(header[0x1C]) != 24)
        {
            printf("Not a correct BMP file\n");
            //return NULL;
        }

        // Read the information about the image
        dataPos     = *(int*)&(header[0x0A]);
        imageSize   = *(int*)&(header[0x22]);
        outWidth    = *(int*)&(header[0x12]);
        outHeight   = *(int*)&(header[0x16]);

        // Some BMP files are misformatted, guess missing information
        if (imageSize == 0)
            imageSize = outWidth * outHeight * 3; // 3 : one byte for each Red, Green and Blue component
        if (dataPos == 0)
            dataPos = 54; // The BMP header is done that way

        // Create a buffer
        //data = new unsigned char [imageSize];

        // Read the actual data from the file into the buffer
        fread(&dataTex[0][number], 1, imageSize, file);

        // Everything is in memory now, the file wan be closed
        fclose (file);

        // if (flipY)
        // {
        //     // swap y-axis
        //     unsigned char * tmpBuffer = new unsigned char[outWidth*3];

        //     int size = outWidth*3;

        //     for (int i=0;i<outHeight/2;i++){

        //     // copy row i to tmp
        //     //memcpy_s(tmpBuffer,size,data+outWidth*3*i,size);
        //     memcpy(tmpBuffer,data+outWidth*3*i,size);

        //     // copy row h-i-1 to i
        //     //memcpy_s(data+outWidth*3*i, size, data+outWidth*3*(outHeight-i-1), size);
        //     memcpy(data+outWidth*3*i, data+outWidth*3*(outHeight-i-1), size);

        //     // copy tmp to row h-i-1
        //     //memcpy_s(data+outWidth*3*(outHeight-i-1), size,tmpBuffer, size);
        //     memcpy(data+outWidth*3*(outHeight-i-1),tmpBuffer, size);

        //     }
        //     delete [] tmpBuffer;
        // }

        //return data;
    }

    void
    applyTexture(void)
    {
        vector<float> texture1(8);
        texture1[0] = 0;
        texture1[1] = 0;
        texture1[2] = 1;
        texture1[3] = 0;
        texture1[4] = 0;
        texture1[5] = 1;
        texture1[6] = 1;
        texture1[7] = 1;

        glBindBuffer(GL_ARRAY_BUFFER, textureBuffer);
        glBufferData(GL_ARRAY_BUFFER, texture1.size() * sizeof(float), &texture1[0], GL_STATIC_DRAW);
        //glBindBuffer(GL_ARRAY_BUFFER, textureBuffer);
        glTexCoordPointer(2, GL_FLOAT, 0, 0);
        glEnable(GL_TEXTURE_2D);                          // Enable Texture Mapping ( NEW )
        glGenTextures(1, &textureID);
        glBindTexture(GL_TEXTURE_2D, textureID);
        //mexPrintf("WidthTex:%d:", widthTex);
        //mexPrintf("HeightTex:%d:", heightTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, widthTex, heightTex, 0, GL_BGR, GL_UNSIGNED_BYTE, &dataTex[0][texture_number]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    }

    static void
    drawSkyVertex(fVec3 position, float texture)
    {
        glTexCoord1f(texture);
        glVertex3fv(&position[0]);
    }

    static void
    drawGroundAndSky(float farClipDistance)
    {
        glDisable(GL_TEXTURE_2D);
        if (!initialized)
        {
            initialized = true;

            // Create a texture to use afor (int i=0;i<99999;i++) tic;result = Simulator.set_params('1.bmp',15,1);toc; ends the sky.
            glGenTextures(1, &skyTexture);
            glBindTexture(GL_TEXTURE_1D, skyTexture);
            glTexParameterf(GL_TEXTURE_1D, GL_GENERATE_MIPMAP, GL_TRUE);
            const int width = 256;
            float skyImage[3*width];
            for (int i = 0; i < width; i++)
            {
                float fract = pow(i / (float) width, 1.8f);
                skyImage[3 * i] = fract;
                skyImage[3 * i + 1] = fract;
                skyImage[3 * i + 2] = 1;
            }
            glTexImage1D(GL_TEXTURE_1D, 0, 3, width, 0, GL_RGB, GL_FLOAT, skyImage);

            // Create a texture to use as the ground.
            srand(0);
            glGenTextures(1, &groundTexture);
            glBindTexture(GL_TEXTURE_2D, groundTexture);
            glTexParameterf(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
            float groundImage[width*width];
            for (int i = 0; i < width; i++)
            {
                float x = i / (float) width;
                for (int j = 0; j < width; j++)
                {
                    float y = j / (float) width;
                    double line = min(min(min(x, y), 1.0f - x), 1.0f - y);
                    float noise = (rand() % 255) / 255.0f - 0.5f;
                    groundImage[i * width + j] = pow(line, 0.1) * (0.35f + noise);
                }
            }
            glTexImage2D(GL_TEXTURE_2D, 0, 1, width, width, 0, GL_RED, GL_FLOAT, groundImage);
        }

        // Draw the box to represent the sky.
        float viewDistance = farClipDistance * 0.5f;
        fVec3 center = fVec3(0, 0, 0);
        const float sign = (float)groundNormal.getDirection(); // 1 or -1
        const CoordinateAxis axis = groundNormal.getAxis();
        float top = center[axis] + sign * viewDistance;
        center[axis] = 0;
        fVec3 offset1 = fVec3(0, 0, 3);
        fVec3 offset2 = (offset1 % fUnitVec3(axis));
        fVec3 corner1 = center + (-offset1 - offset2);
        fVec3 corner2 = center + (offset1 - offset2);
        fVec3 corner3 = center + (offset1 + offset2);
        fVec3 corner4 = center + (-offset1 + offset2);
        glActiveTexture(GL_TEXTURE0);
        glEnable(GL_TEXTURE_1D);
        glBindTexture(GL_TEXTURE_1D, skyTexture);
        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glDisable(GL_DEPTH_TEST);
        glDepthMask(GL_FALSE);
        glBegin(GL_QUAD_STRIP);
        fVec3 offset3 = sign * fUnitVec3(axis);
        drawSkyVertex(corner1 + top * offset3, 0);
        drawSkyVertex(corner1 + groundHeight * offset3, 1);
        drawSkyVertex(corner2 + top * offset3, 0);
        drawSkyVertex(corner2 + groundHeight * offset3, 1);
        drawSkyVertex(corner3 + top * offset3, 0);
        drawSkyVertex(corner3 + groundHeight * offset3, 1);
        drawSkyVertex(corner4 + top * offset3, 0);
        drawSkyVertex(corner4 + groundHeight * offset3, 1);
        drawSkyVertex(corner1 + top * offset3, 0);
        drawSkyVertex(corner1 + groundHeight * offset3, 1);
        glEnd();
        glDisable(GL_TEXTURE_1D);
        glColor3f(0, 0, 1);
        glBegin(GL_QUADS);
        drawSkyVertex(corner1 + 0.99f * top * offset3, 0);
        drawSkyVertex(corner2 + 0.99f * top * offset3, 0);
        drawSkyVertex(corner3 + 0.99f * top * offset3, 0);
        drawSkyVertex(corner4 + 0.99f * top * offset3, 0);
        glEnd();
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

        // Draw the ground plane.
        center = (0, 0, 0);
        center[1] = 0;
        corner1 = center + fVec3(-farClipDistance, 0, -farClipDistance);
        corner2 = center + fVec3(farClipDistance, 0, -farClipDistance);
        corner3 = center + fVec3(farClipDistance, 0, farClipDistance);
        corner4 = center + fVec3(-farClipDistance, 0, farClipDistance);
        // We need to calculate the GL transform T_GP that gives the ground
        // plane's coordinate frame P in the ground frame G.
        Mat<4, 4, GLfloat> T_GP(1);
        fVec3::updAs(&T_GP(YAxis)[0]) = fVec3(fUnitVec3(groundNormal)); // signed
        fVec3::updAs(&T_GP(XAxis)[0]) = fVec3(fUnitVec3(axis.getPreviousAxis()));
        fVec3::updAs(&T_GP(ZAxis)[0]) = fVec3(fUnitVec3(axis.getNextAxis()));
        T_GP[axis][3] = sign * groundHeight;
        glActiveTexture(GL_TEXTURE0);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, groundTexture);
        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_BLEND);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glDisable(GL_CULL_FACE);
        glPushMatrix();
        glMultMatrixf(&T_GP[0][0]);
        glDepthRange(0.01, 1.0);
        float color2[] = {1.0f, 0.8f, 0.7f};
        glBindTexture(GL_TEXTURE_2D, groundTexture);
        glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, color2);
        glBegin(GL_QUADS);
        glColor3f(0.3f, 0.2f, 0.0f);
        glTexCoord2d(2.0f * corner1[0], 2.0f * corner1[2]);
        glVertex3f(corner1[0], corner1[1], corner1[2]);
        glTexCoord2d(2.0f * corner2[0], 2.0f * corner2[2]);
        glVertex3f(corner2[0], corner2[1], corner2[2]);
        glTexCoord2d(2.0f * corner3[0], 2.0f * corner3[2]);
        glVertex3f(corner3[0], corner3[1], corner3[2]);
        glTexCoord2d(2.0f * corner4[0], 2.0f * corner4[2]);
        glVertex3f(corner4[0], corner4[1], corner4[2]);
        glEnd();
        glDepthRange(0.0, 1.0);
        glEnable(GL_CULL_FACE);
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);
    }

    void
    drawBox(void)
    {
        glPushMatrix();
        //glRotatef(angle, 0.0, 1.0, 0.0);
        glTranslatef(0.0, 1.578, 2.0 - this->distance);
        glScalef(this->planeScale, this->planeScale, this->planeScale);
        applyTexture();

        for (i = 1; i < 6; i++)
        {
            glBegin(GL_QUADS);
            glEnable(GL_TEXTURE_2D);
            glNormal3fv(&n[i][0]);
            glTexCoord2f(1.0f, 0.0f);
            glVertex3fv(&v[faces[i][0]][0]);
            glTexCoord2f(1.0f, 1.0f);
            glVertex3fv(&v[faces[i][1]][0]);
            glTexCoord2f(0.0f, 1.0f);
            glVertex3fv(&v[faces[i][2]][0]);
            glTexCoord2f(0.0f, 0.0f);
            glVertex3fv(&v[faces[i][3]][0]);
            glEnd();
        }
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);
    }


    void
    display(int eye)
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        if (eye == 1)
        {
            // left eye
            gluLookAt(-0.028, 1.578, 2.0, /* eye is at (0,0,5) */
            -0.028 + sin(M_PI * this->angle / 180.), 1.578, 2.0 - cos(M_PI * this->angle / 180.), /* center is at (0,0,0) */
            0.0, 1.0, 0.); /* up is in positive Y direction */
        }
        else
        {
            // right eye with strabismus
            gluLookAt(0.028, 1.578, 2.0, /* eye is at (0,0,5) */
            0.028 + sin(-M_PI * (this->angle + this->strangle) / 180.), 1.578, 2.0 - cos(-M_PI * (this->angle + this->strangle) / 180.), /* center is at (0,0,0) */
            0.0, 1.0, 0.); /* up is in positive Y direction */
        }
        drawGroundAndSky(16.01);
        drawBox();
        glutSwapBuffers();
    }


    void
    request(int eye)
    {
        int width = ((DefaultWindowWidth + 3) / 4) * 4; // must be a multiple of 4 pixels
        int height = DefaultWindowHeight;

        GLuint frameBuffer, colorBuffer, depthBuffer;

        // Create offscreen buffers for rendering the image.
        glGenFramebuffersEXT(1, &frameBuffer);
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameBuffer);
        glGenRenderbuffersEXT(1, &colorBuffer);
        glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, colorBuffer);
        glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGB8, width, height);
        glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_RENDERBUFFER_EXT, colorBuffer);
        glGenRenderbuffersEXT(1, &depthBuffer);
        glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, depthBuffer);
        glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24, width, height);
        glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, depthBuffer);

        display(eye);

        glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, &bufferData[0]);
        glDeleteRenderbuffersEXT(1, &colorBuffer);
        glDeleteRenderbuffersEXT(1, &depthBuffer);
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
        glDeleteFramebuffersEXT(1, &frameBuffer);

        const int rowLength = 3 * width;
        for (int row = 0; row < height / 2; ++row)
        {
            const int base1 = row * rowLength;
            const int base2 = (height - 1 - row) * rowLength;
            for (int i = 0; i < rowLength; i++)
            {
                unsigned char temp = bufferData[base1 + i];
                bufferData[base1 + i] = bufferData[base2 + i];
                bufferData[base2 + i] = temp;
            }
        }
    }

    void
    initRenderer()
    {
        int argc = 1;
        char *argv[1] = {(char*)"Something"};
        if (init == false)
        {
            glutInit(&argc, argv);
            glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
            glutInitWindowSize(DefaultWindowWidth, DefaultWindowHeight);
            init = true;
        }
        glutCreateWindow("OpenEyeSim");
        glutHideWindow();
        /* Setup cube vertex data. */
        v[0][0] = v[1][0] = v[2][0] = v[3][0] = -0.5;
        v[4][0] = v[5][0] = v[6][0] = v[7][0] = 0.5;
        v[0][1] = v[1][1] = v[4][1] = v[5][1] = -0.5;
        v[2][1] = v[3][1] = v[6][1] = v[7][1] = 0.5;
        v[0][2] = v[3][2] = v[4][2] = v[7][2] = 0;
        v[1][2] = v[2][2] = v[5][2] = v[6][2] = 0;

        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
        glLightfv(GL_LIGHT0, GL_POSITION, light_position);
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, light_ambient);
        glClearColor(1, 1, 1, 1);
        glEnable(GL_LIGHT0);

        /* Use depth buffering for hidden surface elimination. */
        glEnable(GL_DEPTH_TEST);

        /* Setup the view of the cube. */
        glMatrixMode(GL_PROJECTION);
        gluPerspective(50, 1.33333333, (GLfloat)0.1, 10);
        glMatrixMode(GL_MODELVIEW);
        gluLookAt(0.0, 1.578, 2.0,  /* eye is at (0,0,5) */
        0.0, 1.578, 0.0,            /* center is at (0,0,0) */
        0.0, 1.0, 0.);              /* up is in positive Y direction */

        request(1);
        display(1);
    }

    void
    reinitRenderer()
    {
        if (init == false)
        {
            glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
            glutInitWindowSize(DefaultWindowWidth, DefaultWindowHeight);
            init = true;
        }
        glutCreateWindow("OpenEyeSim");
        glutHideWindow();
        /* Setup cube vertex data. */
        v[0][0] = v[1][0] = v[2][0] = v[3][0] = -0.5;
        v[4][0] = v[5][0] = v[6][0] = v[7][0] = 0.5;
        v[0][1] = v[1][1] = v[4][1] = v[5][1] = -0.5;
        v[2][1] = v[3][1] = v[6][1] = v[7][1] = 0.5;
        v[0][2] = v[3][2] = v[4][2] = v[7][2] = 0;
        v[1][2] = v[2][2] = v[5][2] = v[6][2] = 0;

        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
        glLightfv(GL_LIGHT0, GL_POSITION, light_position);
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, light_ambient);
        glClearColor(1, 1, 1, 1);
        glEnable(GL_LIGHT0);

        /* Use depth buffering for hidden surface elimination. */
        glEnable(GL_DEPTH_TEST);

        /* Setup the view of the cube. */
        glMatrixMode(GL_PROJECTION);
        gluPerspective(50, 1.33333333, (GLfloat)0.1, 10);
        glMatrixMode(GL_MODELVIEW);
        gluLookAt(0.0, 1.578, 2.0,  /* eye is at (0,0,5) */
        0.0, 1.578, 0.0,            /* center is at (0,0,0) */
        0.0, 1.0, 0.);              /* up is in positive Y direction */

        request(1);
        display(1);
    }

public:
    unsigned char bufferData[230400];
    // unsigned char dataTex[3 * 480 * 480][10]; // Celine's image resolution
    unsigned char dataTex[1024 * 1024 * 3][100]; // vanHateren image resolution

private:
    // Database implementation.
    GLuint textureBuffer;
    string texture;
    int texture_number;
    double angle;
    float distance;
    bool init;
    double strangle;
    float planeScale;
};

// Instance manager for OpenEyeSim.
template class mexplus::Session<OpenEyeSim>;

namespace {

// Defines MEX API for new.
MEX_DEFINE(new) (int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    InputArguments input(nrhs, prhs, 1);
    OutputArguments output(nlhs, plhs, 1);
    output.set(0, Session<OpenEyeSim>::create(new OpenEyeSim()));
}

// Defines MEX API for delete.
MEX_DEFINE(delete) (int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[])
{
    InputArguments input(nrhs, prhs, 1);
    OutputArguments output(nlhs, plhs, 0);
    Session<OpenEyeSim>::destroy(input.get(0));
}

// Defines MEX API for set_params (non const method).
MEX_DEFINE(set_params) (int nlhs, mxArray* plhs[],
                        int nrhs, const mxArray* prhs[])
{
    InputArguments input(nrhs, prhs, 6);
    OutputArguments output(nlhs, plhs, 0);
    OpenEyeSim* osim = Session<OpenEyeSim>::get(input.get(0));
    osim->set_params(input.get<int>(1), input.get<double>(2), input.get<float>(3), input.get<double>(4), input.get<float>(5));
}

// Defines MEX API for set_params (non const method).
MEX_DEFINE(add_texture) (int nlhs, mxArray* plhs[],
                         int nrhs, const mxArray* prhs[])
{
    InputArguments input(nrhs, prhs, 3);
    OutputArguments output(nlhs, plhs, 0);
    OpenEyeSim* osim = Session<OpenEyeSim>::get(input.get(0));
    osim->add_texture(input.get<int>(1), input.get<string>(2));
}

// Defines MEX API for generate left (non const method).
MEX_DEFINE(generate_left) (int nlhs, mxArray* plhs[],
                           int nrhs, const mxArray* prhs[]) {
    InputArguments input(nrhs, prhs, 1);
    OutputArguments output(nlhs, plhs, 1);
    OpenEyeSim* osim = Session<OpenEyeSim>::get(input.get(0));
    osim->request(1);
    std::vector<unsigned char> v(std::begin(osim->bufferData), std::end(osim->bufferData));
    plhs[0] =  MxArray::from(v);
    v.clear();
    glDeleteTextures(1, &textureID);
}

MEX_DEFINE(generate_right) (int nlhs, mxArray* plhs[],
                            int nrhs, const mxArray* prhs[])
{
    InputArguments input(nrhs, prhs, 1);
    OutputArguments output(nlhs, plhs, 1);
    OpenEyeSim* osim = Session<OpenEyeSim>::get(input.get(0));
    osim->request(2);
    std::vector<unsigned char> v2(std::begin(osim->bufferData), std::end(osim->bufferData));
    plhs[0] =  MxArray::from(v2);
    v2.clear();
    glDeleteTextures(1, &textureID);
}

// Defines MEX API for request (const method).
MEX_DEFINE(request) (int nlhs, mxArray* plhs[],
                     int nrhs, const mxArray* prhs[])
{
    InputArguments input(nrhs, prhs, 1);
    OutputArguments output(nlhs, plhs, 1);
    OpenEyeSim* osim = Session<OpenEyeSim>::get(input.get(0));
    osim->request(1);
}

// Defines MEX API for request (const method).
MEX_DEFINE(initRenderer) (int nlhs, mxArray* plhs[],
                          int nrhs, const mxArray* prhs[])
{
    InputArguments input(nrhs, prhs, 1);
    OutputArguments output(nlhs, plhs, 1);
    OpenEyeSim* osim = Session<OpenEyeSim>::get(input.get(0));
    osim->initRenderer();
}

// Defines MEX API for request (const method).
MEX_DEFINE(reinitRenderer) (int nlhs, mxArray* plhs[],
                            int nrhs, const mxArray* prhs[])
{
    InputArguments input(nrhs, prhs, 1);
    OutputArguments output(nlhs, plhs, 1);
    OpenEyeSim* osim = Session<OpenEyeSim>::get(input.get(0));
    osim->reinitRenderer();
}

} // namespace

MEX_DISPATCH // Don't forget to add this if MEX_DEFINE() is used.
