#ifndef VIEWER_H
#define VIEWER_H 

#include "mh/base/defs.h"
#include "mh/base/imports.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

//#include "mh/ext/imgui/imgui.h"
//#include "mh/ext/imgui/imgui_impl_glfw.h"

#include "mh/gpu/shader.h"

namespace mh
{

inline void glfwErrorCallback(int error, const char * desc)
{
    fprintf(stderr, "Error %d: %s\n", error, desc);    
}

class Viewer
{
public:
    typedef std::map<std::string, std::shared_ptr<Shader> > ShaderMap;
    
    // constructors & destructors
    Viewer(int width=MH_DEFAULT_VIEWER_WIDTH, int height=MH_DEFAULT_VIEWER_HEIGHT, const std::string & title=MH_DEFAULT_VIEWER_TITLE)
      : m_width(width), m_height(height), m_title(title), m_shouldClose(false)
    { baseInit(); }

    virtual     ~Viewer() {}

    // baseInit, called by the constructors
    int         baseInit       (void);

    // initialization function for derived viewers
    virtual int init           (void) = 0;

    // main loop, this is the function that should get called in the main
    // loop of the app
    virtual void mainLoop      (void) = 0;
    
    // makes the GL context of this viewer current
    void         makeCurrent   (void);

    // swaps the back and front buffers
    void         swap          (void)
    { glfwSwapBuffers(m_window); }

    // shouldClose, used to propagate window-close events
    void         setShouldClose(bool shouldClose)
    { m_shouldClose = shouldClose; }
    bool         shouldClose   (void)             const
    { return m_shouldClose || glfwWindowShouldClose(m_window); }

    // add shader to collection
    void         addShader     (std::string name, std::shared_ptr<Shader> shader)
    { m_shaders.insert(std::pair<std::string, std::shared_ptr<Shader> >(name, shader)); }

    int          getWidth      (void)             const
    { return m_width; }
    int          getHeight     (void)             const
    { return m_height; }

protected:
    void         pollSize();

    GLFWwindow * m_window;

    int          m_width;
    int          m_height;

    std::string  m_title;

    ShaderMap    m_shaders;

    bool         m_shouldClose;

    void *       m_imguiState;

private:

}; // class Viewer



} // namespace mh

#endif /* VIEWER_H */
