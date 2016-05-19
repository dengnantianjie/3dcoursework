//
// Created by bontius on 25/01/16.
//

#include "gui/meshviewer.h"
#include "mh/3d/camera.h"
//#include "mh/gui/dialogs.h"
#include "mh/io/meshio.h"
#include "mh/gpu/shader.h"
#include <vector>
#include "Algorithms.h"
#include <ctime>

bool ICP_triger = false;
int laplace_flag = 1;
namespace
{
    
using namespace mh;


std::shared_ptr<Shader> compileDefaultShader()
{
    auto shader = std::make_shared<Shader>();
    shader->vertexShader(
        GLSL(330,
            layout(location = 0) in vec3 vertexPosition;
            layout(location = 1) in vec3 vertexNormal;
            layout(location = 2) in vec3 vertexColor;

            uniform mat4 cameraToClip;
            uniform mat4 worldToCamera;
            uniform mat4 modelToWorld;

            out vec3 position;
            out vec3 normal;
            //out vec3 color;
            out vec3 color ;

            void main()
            {
                position = (worldToCamera * modelToWorld * vec4(vertexPosition, 1.0)).xyz;
                normal   = (worldToCamera * modelToWorld * vec4(vertexNormal,   0.0)).xyz;
                color    = vertexColor;
                //color = vec4(1.0,0.0,0.0,0.3);

                gl_Position = cameraToClip * vec4(position, 1.0);
            }
        ));

    shader->fragmentShader(
        GLSL(330,
            in vec3 position;
            in vec3 normal;
            in vec3 color;
            

            out vec4 fragColor;
             
            
            void main()
            {
                
                vec3 pixelNormal = normalize(normal);
                vec3 viewDirection = normalize(-position);
                float clampedCosine = max(dot(pixelNormal, viewDirection), 0.0);
                //fragColor = vec4(clampedCosine * color, 1.0);
                //fragColor = vec4(1.0,0.0,0.0,0.3);
                //fragColor = vec4(clampedCosine/3,0.3);
                vec4 diffuseColor = vec4(clampedCosine * color, 1.0);
                vec4 ambientColor = vec4(color,1.0);
                fragColor = 0.6 * diffuseColor + 0.4 * ambientColor;

                

                
            }
        ));

    shader->link();

    return shader;
}
    


std::shared_ptr<Shader> compileNormalShader()
{
    auto shader = std::make_shared<Shader>();

    shader->vertexShader(
        GLSL(330,
            layout(location = 0) in vec3 vertexPosition;
            layout(location = 1) in vec3 vertexNormal;
            layout(location = 2) in vec3 vertexColor;

            uniform mat4 cameraToClip;
            uniform mat4 worldToCamera;
            uniform mat4 modelToWorld;

            out vec3 position;
            out vec3 normal;
            out vec3 color;

            out VS_OUT {
            vec3 normal;
            vec3 normal_color;
        } vs_out;

            void main()
            {
                position = (worldToCamera * modelToWorld * vec4(vertexPosition, 1.0)).xyz;
                normal   = normalize((worldToCamera * modelToWorld * vec4(vertexNormal, 0.0)).xyz);
                color    = vertexColor;
                //color = vec3(1.0,0.0,0.0);

                gl_Position = cameraToClip * vec4(position, 1.0);

                //mat4 normalMatrix   = transpose(inverse(worldToCamera * modelToWorld));
                /*vec3 tiNormal = (normalMatrix * vec4(vertexNormal, 0.0)).xyz;
                     tiNormal = normalize(tiNormal);*/
                /*vec3 trNormal = (worldToCamera * modelToWorld * vec4(vertexNormal, 0.0)).xyz;
                     trNormal = normalize(trNormal);*/

                vs_out.normal = (cameraToClip * vec4(normal, 0.0)).xyz;
                vs_out.normal_color = abs(normalize(vertexNormal));
                //vs_out.normal = normalize((cameraToClip * vec4(normalMatrix * vertexNormal, 1.0)));
            }
        ));

    shader->geometryShader(
        GLSL(330,
            layout (points) in;
            layout (line_strip, max_vertices = 2) out;

            in VS_OUT {
            vec3 normal;
            vec3 normal_color;
        } gs_in[];

            out vec3 vertex_color;

            const float MAGNITUDE = 0.01;

            void GenerateLine(int index)
            {
                gl_Position = gl_in[index].gl_Position;
                vertex_color = abs(gs_in[index].normal_color);
                EmitVertex();
                gl_Position = gl_in[index].gl_Position + vec4(gs_in[index].normal, 0.0f) * MAGNITUDE;
                vertex_color = abs(gs_in[index].normal_color);
                EmitVertex();
                EndPrimitive();
            }

            void main()
            {
                GenerateLine(0);
                //GenerateLine(1);
                //GenerateLine(2);
            }
        ));

    shader->fragmentShader(
        GLSL(330,
            out vec4 fragColor;

            in vec3 vertex_color;

            void main()
            {
                fragColor = vec4(vertex_color, 1.0);
            }
        ));

    shader->link();

    return shader;
}

}

namespace mh
{

int MeshViewer::init(void)
{
    //// OpenGL setup
    glEnable    (GL_DEPTH_TEST);
    glDepthFunc (GL_LESS);

    glDisable   (GL_CULL_FACE);
    glCullFace  (GL_BACK);
    glFrontFace (GL_CCW);

    glEnable    (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    

    glPointSize (1.5);
    
    //// Default camera setup
    float near = 0.01f;
    float far  = 100.0f;
    float fov  = MH_CAM_DEFAULT_FOV;
    float speed = MH_CAM_DEFAULT_SPEED;
    float sensitivity = 1.0f / 1000;

    m_camera = std::make_shared<Camera>(near, far, fov, m_width / m_height, speed, sensitivity);
    std::cout<<"nanan"<<m_camera->m_forward(2,0)<<std::endl;

    //// Default shader setup
    
    addShader("default", compileDefaultShader());
    addShader("normals", compileNormalShader());
    

    // setup glfw interaction
    glfwSetKeyCallback(m_window,key_callback);
    glfwSetMouseButtonCallback(m_window, mouse_callback );
    glfwSetCursorPosCallback(m_window, mouse_cursor_callback);
    glfwSetScrollCallback(m_window, mouse_scroll_callback);
    //glfwSetCharCallback(m_window,char_callback);
    glfwSetWindowUserPointer( m_window, this ); // communicates this Viewer object to glfw, so it can be used in key_callback

    return 0;
}

void MeshViewer::mainLoop()
{
    
    glfwPollEvents();

    //ImGui_ImplGlfw_NewFrame();

    pollSize();
    m_camera->setAspect((float) m_width / (float) m_height);

    //ImGuiIO& io = ImGui::GetIO();
    //if (io.KeysDown[io.KeyMap[ImGuiKey_Escape]])
    //{
    //    setShouldClose(true);
    //}

    //if (m_scene)
    //{
    //    updateCameraWithImgui(*m_camera, io, m_scene->getCenter());
    //    updateCameraWithGlfw( *m_camera, m_scene->getCenter() );
    //}

//    if (ImGui::Button("Load mesh"))
//    {
//        std::string filename = openFileDialog("obj,ply");
//        if (filename != "")
//        {
//            std::vector<std::shared_ptr<Mesh> > meshes = loadMeshesFromFile(filename);
//
//            m_scene->addMeshes(meshes);
//            m_scene->centerToCenterOfMass();
//
//            setCameraLookatScene(*m_camera, *m_scene);
//        }
//    }

    //ImGui::Checkbox("Show normals", &m_showNormals);

    //// Draw calls
    glViewport(0, 0, m_width, m_height);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    

    glLineWidth(2.0f);
    

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // wireframe
    if (m_scene)
    {
        
        m_scene->draw(m_shaders["default"], m_camera);
        if (m_showNormals)
        {
            
            m_scene->draw(m_shaders["normals"], m_camera);
        }
    }

    glDisable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    //ImGui::Render();
}


void MeshViewer::loadMesh(std::shared_ptr<Mesh>&mesh_,const std::string &path)
{
    
    if ( !path.empty() )
    {

        std::vector<std::shared_ptr<Mesh> > meshes = loadMeshesFromFile(path);
        mesh_ = meshes[0];
        m_scene->addMeshes(meshes);
        m_scene->centerToCenterOfMass();
        setCameraLookatScene(*m_camera, *m_scene);
        
    }
    //std::cout << "myvector stores " << int(meshes_vector.size()) << std::endl;
} //...loadMesh()

void MeshViewer::handleKeyPress(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if ( !m_scene || !m_camera )
        return;

    Eigen::Vector3f position_delta(0.0f, 0.0f, 0.0f);
    Eigen::Vector3f right;
    
    right = m_camera->getForward().cross(m_camera->getUp()).normalized();
    
    if ( action == GLFW_PRESS )
    {
        if ( key == GLFW_KEY_ESCAPE )
        {
            this->setShouldClose(true);
        }
        else if ( key == GLFW_KEY_W )
        {
            position_delta += m_camera->getForward() * m_camera->getSpeed();
        }
        else if ( key == GLFW_KEY_S )
        {
            position_delta -= m_camera->getForward() * m_camera->getSpeed();
        }
        else if ( key == GLFW_KEY_D )
        {
            position_delta += right * m_camera->getSpeed();
        }
        else if ( key == GLFW_KEY_A )
        {
            position_delta -= right * m_camera->getSpeed();
        }
        else if ( key == GLFW_KEY_R )
        {
            //position_delta += m_camera->getUp() * m_camera->getSpeed();
            //rotateMesh(30, 1, 0, 0);
            
        }
        else if ( key == GLFW_KEY_F )
        {
            position_delta -= m_camera->getUp() * m_camera->getSpeed();
        }

//////////////////////////////////////////////////////Assignment 2 ///////////////
        else if ( key == GLFW_KEY_1 )
        {
            color_to_show(om1, uniformcurvature_n,mesh11);
        }
         else if ( key == GLFW_KEY_2 )
        {
            color_to_show(om1, uniformcurvature_n_toshow,mesh11);
        }
        else if ( key == GLFW_KEY_3 )
        {
            color_to_show(om1, gausscurvature_n,mesh11);
        }
        else if ( key == GLFW_KEY_4 )
        {
            color_to_show(om1, gausscurvature_n_toshow,mesh11);
        }

        else if ( key == GLFW_KEY_5 )
        {
            color_to_show(om1, lbcurvature_n,mesh11);
        }
        else if ( key == GLFW_KEY_6 )
        {
            color_to_show(om1, lbcurvature_n_toshow,mesh11);
        }

         else if ( key == GLFW_KEY_7 )
        {
            color_to_show(om1, principle_max_uni,mesh11);
        }

         else if ( key == GLFW_KEY_8 )
        {
            color_to_show(om1, principle_min_uni,mesh11);
        }

         else if ( key == GLFW_KEY_9 )
        {
            color_to_show(om1, principle_max_lb,mesh11);
        }

         else if ( key == GLFW_KEY_0 )
        {
            color_to_show(om1, principle_min_lb,mesh11);
        }

         else if ( key == GLFW_KEY_Z )
        {
            color_to_show(om1, principle_max_uni_toshow,mesh11);
        }

         else if ( key == GLFW_KEY_X )
        {
            color_to_show(om1, principle_min_uni_toshow,mesh11);
        }

         else if ( key == GLFW_KEY_C )
        {
            color_to_show(om1, principle_max_lb_toshow,mesh11);
        }

         else if ( key == GLFW_KEY_V )
        {
            color_to_show(om1, principle_min_lb_toshow,mesh11);
        }
/////////////// part 2 mesh smoothing //////////
        else if ( key == GLFW_KEY_ENTER)
        {
            reset_meshs(mesh11,om1);
        }
        
        else if ( key == GLFW_KEY_J )
        {
             om_explicit_laplace_20 = explicit_laplace_smooth(om1, mesh11,20, 0.5);
             color_to_show(om_explicit_laplace_20, uniformcurvature_n_toshow,mesh11);
        }

        else if ( key == GLFW_KEY_K)
        {
            om_explicit_laplace_20 = explicit_laplace_smooth(om1, mesh11,10, 0.5);
            color_to_show(om_explicit_laplace_20, uniformcurvature_n_toshow,mesh11);
            
        }

        else if ( key == GLFW_KEY_UP)
        {
            const clock_t begin_time = clock();
            om1 = explicit_laplace_smooth(om1, mesh11, 1, 0.5);
            std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC;
            color_to_show(om1, uniformcurvature_n_toshow,mesh11);

        }

        else if ( key == GLFW_KEY_N)
        {
            om1 = implicit_laplace_smooth(om1,mesh11, 20,0.0001);
            color_to_show(om1, uniformcurvature_n_toshow,mesh11);
        }

        else if ( key == GLFW_KEY_M)
        {
            om1 = implicit_laplace_smooth(om1,mesh11, 10, 0.0001);
            color_to_show(om1, uniformcurvature_n_toshow,mesh11);
        }

        else if ( key == GLFW_KEY_DOWN)
        {
            const clock_t begin_time = clock();


            om1 = implicit_laplace_smooth(om1,mesh11, 1, 0.0002);
            std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC;
            color_to_show(om1, uniformcurvature_n_toshow,mesh11);
        }





////////////////////////////////////////////////////////////////////////////////////        
        else if ( key == GLFW_KEY_I )
        {
            for (int i = 0; i < 25; i++)
            {
                meshAlign(10);
            }
        }
        else if ( key == GLFW_KEY_L )
        {
            meshAlignAll(10);
        }
    }

    m_camera->setPosition(m_camera->getPosition() + position_delta);
} //...handleKeyPress()


void MeshViewer::handleMouseScroll( GLFWwindow* window, double xoffset, double yoffset)
{
    if ( !m_scene || !m_camera )
        return;

    std::cout << "[MeshViewer::handleMouseScroll] " << xoffset << "," << yoffset << std::endl;
    m_camera->setFOV(m_camera->getFOV() + (yoffset * MH_DEFAULT_ZOOM_SPEED));
} //...handleMouseScroll()

void MeshViewer::handleMouseButton(GLFWwindow* window, int button, int action, int mods)
{
    if ( !m_scene || !m_camera )
        return;

    if ( action == GLFW_PRESS )
    {
        if ( button == GLFW_MOUSE_BUTTON_1 || button == GLFW_MOUSE_BUTTON_2 )
        {
            if ( button == GLFW_MOUSE_BUTTON_1 )
                m_button1Down = true;
            else
                m_button2Down = true;

            double xpos(0.), ypos(0.);
            glfwGetCursorPos(window, &xpos, &ypos);
            m_lastCursorPos(0) = xpos;
            m_lastCursorPos(1) = ypos;
        }
    }
    else if ( action == GLFW_RELEASE )
    {
        if ( button == GLFW_MOUSE_BUTTON_1 )
        {
            m_button1Down = false;
        }
        else if ( button == GLFW_MOUSE_BUTTON_2 )
        {
            m_button2Down = false;
        }
    }
}

void MeshViewer::handleMouseMove( GLFWwindow* window, double xpos, double ypos )
{
    if ( !m_scene || !m_camera )
        return;

    bool updatePos(false);

    if ( m_button1Down )
    {
        Eigen::Vector2f mouseDelta = (Eigen::Vector2d(xpos,ypos) - m_lastCursorPos).cast<float>();
        float yawAngle   = -mouseDelta(0) * m_camera->getSensitivity();
        float pitchAngle = -mouseDelta(1) * m_camera->getSensitivity();

        Eigen::Quaternionf rotation = Eigen::AngleAxisf(pitchAngle, m_camera->getForward().cross(m_camera->getUp()).normalized()) *
                                      Eigen::AngleAxisf(yawAngle, m_camera->getUp());
        Eigen::Vector3f center = m_scene->getCenter();
        Eigen::Vector3f newCameraPos = Eigen::Translation3f(center) * rotation * Eigen::Translation3f(-center) * m_camera->getPosition();

        m_camera->setForward(rotation * m_camera->getForward());
        m_camera->setUp(m_camera->getForward().cross(m_camera->getUp()).normalized().cross(m_camera->getForward()).normalized());
        m_camera->setPosition(newCameraPos);

        updatePos = true;
    }
    else if ( m_button2Down )
    {
        Eigen::Vector2f mouseDelta = (Eigen::Vector2d(xpos,ypos) - m_lastCursorPos).cast<float>() * m_camera->getSensitivity();
        Eigen::Vector3f translation = m_camera->getUp() * mouseDelta(1) + m_camera->getForward().cross(m_camera->getUp()).normalized() * -mouseDelta(0);
        Eigen::Vector3f newCameraPos = Eigen::Translation3f(translation) * m_camera->getPosition();
        m_camera->setPosition(newCameraPos);
        std::cout << "translation: "<< translation.transpose() << std::endl;

        updatePos = true;
    }

    if ( updatePos )
    {
        m_lastCursorPos(0) = xpos;
        m_lastCursorPos(1) = ypos;
    }
} //...handleMouseMove

} //...namespace

void pause(int dur)
{
    int temp = time(NULL) + dur;
    
    while(temp > time(NULL));
}


