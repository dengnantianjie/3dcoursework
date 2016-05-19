//
// Created by bontius on 25/01/16.
//

#ifndef MH_LIB_MESHVIEWER_H
#define MH_LIB_MESHVIEWER_H

#include "mh/base/defs.h"
#include "mh/base/imports.h"

#include "mh/3d/camera.h"
#include "mh/gui/sceneviewer.h"

extern bool ICP_triger;

namespace mh
{

class MeshViewer : public SceneViewer
{
    public:
        MeshViewer(int width=MH_DEFAULT_VIEWER_WIDTH, int height=MH_DEFAULT_VIEWER_HEIGHT)
            : SceneViewer(width, height), m_button1Down(false), m_button2Down(false), m_lastCursorPos(0.,0.), m_showNormals(true)
        { init(); }

        int  init        (void);

        void mainLoop    (void);

        void loadMesh(std::shared_ptr<Mesh>&mesh_,const std::string& path );

        void handleKeyPress( GLFWwindow* window, int key, int scancode, int action, int mods );
        void handleMouseButton(GLFWwindow* window,  int button, int action, int mods);
        void handleMouseScroll( GLFWwindow* window, double xoffset, double yoffset );
        void handleMouseMove( GLFWwindow* window, double xpos, double ypos );

    protected:
        bool            m_button1Down, m_button2Down;
        Eigen::Vector2d m_lastCursorPos;
    private:
        std::shared_ptr<Camera> m_camera;

        bool                    m_showNormals;

}; // class MeshViewer

inline void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    MeshViewer* meshViewer = reinterpret_cast<MeshViewer*>(glfwGetWindowUserPointer(window)); // this is setup in init
    std::cout << "[meshviewer.h::key_callback] key \'" << key << "\' pressed" << std::endl;
    meshViewer->handleKeyPress(window, key, scancode, action, mods);
}

inline void mouse_callback(GLFWwindow* window, int button, int action, int mods)
{
    MeshViewer* meshViewer = reinterpret_cast<MeshViewer*>(glfwGetWindowUserPointer(window)); // this is setup in init
    std::cout << "[meshviewer.h::mouse_callback] button \'" << button << "\' pressed" << std::endl;
    meshViewer->handleMouseButton(window, button, action, mods);
}

inline void mouse_cursor_callback(GLFWwindow* window, double xpos, double ypos )
{
    MeshViewer* meshViewer = reinterpret_cast<MeshViewer*>(glfwGetWindowUserPointer(window)); // this is setup in init
    meshViewer->handleMouseMove( window, xpos, ypos );
}

inline void mouse_scroll_callback(GLFWwindow* window, double xoffset, double yoffset )
{
    MeshViewer* meshViewer = reinterpret_cast<MeshViewer*>(glfwGetWindowUserPointer(window)); // this is setup in init
    meshViewer->handleMouseScroll( window, xoffset, yoffset);
}

inline void char_callback(GLFWwindow* window, unsigned int key)
{
    std::cout << "[viewer.cpp::char_callback] key " << key << " pressed" << std::endl;
}


} // namespace mh

#endif //MH_LIB_MESHVIEWER_H
