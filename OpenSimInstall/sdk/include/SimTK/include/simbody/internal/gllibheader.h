#ifndef GLLIBHEADER_H
#define GLLIBHEADER_H

#include "SimTKcommon.h"

#include <cstdlib>
#include <string>
#include <algorithm>
#include <set>
#include <vector>
#include <utility>
#include <limits>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <pthread.h>
#include <sys/stat.h>

#ifdef _WIN32
    #include <direct.h>
#endif

// Get gl and glut using the appropriate platform-dependent incantations.
#if defined(__APPLE__)
    // OSX comes with a good glut implementation.
    #include <GLUT/glut.h>
#elif defined(_WIN32)
    #include "glut32/glut.h"    // we have our own private headers
    #include "glut32/glext.h"

    // A Windows-only extension for disabling vsync, allowing unreasonably
    // high frame rates.
    PFNWGLSWAPINTERVALFARPROC wglSwapIntervalEXT;

    // These will hold the dynamically-determined function addresses.

    // These functions are needed for basic Visualizer functionality.
    PFNGLGENBUFFERSPROC glGenBuffers;
    PFNGLBINDBUFFERPROC glBindBuffer;
    PFNGLBUFFERDATAPROC glBufferData;
    PFNGLACTIVETEXTUREPROC glActiveTexture;

    // These are needed only for saving images and movies.
    // Use old EXT names for these so we only require OpenGL 2.0.
    PFNGLGENFRAMEBUFFERSEXTPROC glGenFramebuffersEXT;
    PFNGLGENRENDERBUFFERSEXTPROC glGenRenderbuffersEXT;
    PFNGLBINDFRAMEBUFFEREXTPROC glBindFramebufferEXT;
    PFNGLBINDRENDERBUFFEREXTPROC glBindRenderbufferEXT;
    PFNGLRENDERBUFFERSTORAGEEXTPROC glRenderbufferStorageEXT;
    PFNGLFRAMEBUFFERRENDERBUFFEREXTPROC glFramebufferRenderbufferEXT;
    PFNGLDELETERENDERBUFFERSEXTPROC glDeleteRenderbuffersEXT;
    PFNGLDELETEFRAMEBUFFERSEXTPROC glDeleteFramebuffersEXT;

    // see initGlextFuncPointerIfNeeded() at end of this file
#else
    // Linux: assume we have a good OpenGL 2.0 and working glut or freeglut.
    #define GL_GLEXT_PROTOTYPES
    #include <GL/gl.h>
    #include <GL/glu.h>
    #include <GL/glext.h>
    #include <GL/glut.h>
#endif

// Next, get the functions necessary for reading from and writing to pipes.
#ifdef _WIN32
    #include <io.h>
    #define READ _read
#else
    #include <unistd.h>
    #define READ read
#endif

using namespace SimTK;
using namespace std;

#endif // GLLIBHEADER_H
