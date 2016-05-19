//
//  Algorithms.h
//  
//
//  Created by 邓楠天婕 on 16-2-9.
//
//

#ifndef ____Algorithms__
#define ____Algorithms__
#include <iostream>
#include <stdio.h>
#include <random>
#include <math.h>
#include "gui/meshviewer.h"


#include "mh/base/defs.h"
#include "mh/base/imports.h"

#include "mh/3d/camera.h"
#include "mh/gui/sceneviewer.h"
#include <ctime>

#include <fstream>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"

#include "mh/base/defs.h"
#include "mh/base/imports.h"


#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>


#include <GL/glew.h>



#endif /* defined(____Algorithms__) */

using namespace mh;
struct MyTraits : public OpenMesh::DefaultTraits
{
    typedef OpenMesh::Vec4f Color;
    VertexAttributes (
                      OpenMesh::Attributes::Normal |
                      OpenMesh::Attributes::Color);
    FaceAttributes (
                    OpenMesh::Attributes::Normal |
                    OpenMesh::Attributes::Color);
};

typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> MyMesh;


void meshAlign(int ratio); 
void readPt(Mesh &mesh, ANNpoint &p, int downsample);
double ICPAlgo(std::shared_ptr<Mesh>& mesh11,std::shared_ptr<Mesh>& mesh22, int downsample);
void rotateMesh(std::shared_ptr<Mesh>& mesh11,double degree, double x, double y, double z);
void translateMesh(std::shared_ptr<Mesh>& mesh11);
void pause(int dur);
void addNoise(double sd);

void meshAlignAll(int ratio);
void setmeshcolor(std::shared_ptr<Mesh> mesh11,Eigen::Vector4f colors);
void setNormal(std::shared_ptr<Mesh> mesh11);
Mesh transform_structure(MyMesh om1);
void reset_meshs(std::shared_ptr<Mesh> mesh11,MyMesh om1);
extern std::shared_ptr<Mesh> mesh11,mesh22,mesh33,mesh44,mesh55,mesh66;
extern MyMesh om1,om_explicit_laplace_1,om_explicit_laplace_20,om_explicit_laplace_10;
extern std::shared_ptr<Mesh> mesh_reset;
extern MyMesh om_original;
extern MyMesh om_implicit_laplace_1,om_implicit_laplace_20,om_implicit_laplace_10;


extern OpenMesh::VPropHandleT<OpenMesh::Vec3f> uniformcurvature;
extern OpenMesh::VPropHandleT<MyMesh::Scalar> uniformcurvature_n;
extern OpenMesh::VPropHandleT<MyMesh::Scalar> uniformcurvature_n_toshow;

extern OpenMesh::EPropHandleT<MyMesh::Scalar>  weight;
extern OpenMesh::VPropHandleT<OpenMesh::Vec3f> gausscurvature;
extern OpenMesh::VPropHandleT<MyMesh::Scalar> gausscurvature_n;
extern OpenMesh::VPropHandleT<MyMesh::Scalar> gausscurvature_n_toshow;
extern OpenMesh::VPropHandleT<OpenMesh::Vec3f>  lbcurvature;
extern OpenMesh::VPropHandleT<MyMesh::Scalar>  lbcurvature_n;
extern OpenMesh::VPropHandleT<MyMesh::Scalar>  lbcurvature_n_toshow;


extern OpenMesh::VPropHandleT<MyMesh::Scalar> principle_max_uni;
extern OpenMesh::VPropHandleT<MyMesh::Scalar> principle_min_uni;
extern OpenMesh::VPropHandleT<MyMesh::Scalar> principle_max_lb;
extern OpenMesh::VPropHandleT<MyMesh::Scalar> principle_min_lb;
extern OpenMesh::VPropHandleT<MyMesh::Scalar> principle_max_uni_toshow;
extern OpenMesh::VPropHandleT<MyMesh::Scalar> principle_min_uni_toshow;
extern OpenMesh::VPropHandleT<MyMesh::Scalar> principle_max_lb_toshow;
extern OpenMesh::VPropHandleT<MyMesh::Scalar> principle_min_lb_toshow;
extern OpenMesh::VPropHandleT<float> centerarea;
extern OpenMesh::VPropHandleT<float> sumdegree;

extern OpenMesh::EPropHandleT<MyMesh::Scalar>  eweight_;
extern OpenMesh::VPropHandleT<MyMesh::Scalar>  eweightSum_;


void setmeshcolor_new(std::shared_ptr<Mesh> mesh11,Eigen::Vector4f colors);
void color_to_show(MyMesh &mesh, OpenMesh::VPropHandleT<MyMesh::Scalar> prop, std::shared_ptr<Mesh> mesh11);
Eigen::Vector4f calculate_color(MyMesh::Scalar curve, MyMesh::Scalar min, MyMesh::Scalar max);

MyMesh uniform_curvature(MyMesh mesh);
MyMesh gauss_curvature(MyMesh mesh);
MyMesh lb_mean_curvature(MyMesh mesh_);
MyMesh get_principles(MyMesh mesh, OpenMesh::VPropHandleT<MyMesh::Scalar> h, OpenMesh::VPropHandleT<MyMesh::Scalar> k,int flag);
MyMesh explicit_laplace_smooth(MyMesh mesh, std::shared_ptr<Mesh> mesh11, int iterations, double lamda);
MyMesh implicit_laplace_smooth(MyMesh mesh,  std::shared_ptr<Mesh> mesh11, int step, double lamda);
extern int laplace_flag;

//float calc_area(MyMesh mesh,MyMesh::VertexHandle center, MyMesh::VertexHandle pre, MyMesh::VertexHandle next);
//float calc_degree(MyMesh mesh,MyMesh::VertexHandle center, MyMesh::VertexHandle pre, MyMesh::VertexHandle next);
//
//void color_code(MyMesh &mesh, OpenMesh::VPropHandleT<MyMesh::Scalar> prop);
//OpenMesh::Vec3uc getcolor(MyMesh::Scalar curve, MyMesh::Scalar min, MyMesh::Scalar max);
//MyMesh calc_mean_curvature(MyMesh mesh_);
//MyMesh calc_weights(MyMesh mesh_);
//MyMesh calc_principle(MyMesh mesh, OpenMesh::VPropHandleT<MyMesh::Scalar> h, OpenMesh::VPropHandleT<MyMesh::Scalar> k);
float calc_theta(OpenMesh::Vec3f normal, OpenMesh::Vec3f target_vector);
MyMesh calc_weights(MyMesh mesh_);
float calc_degree(MyMesh mesh, MyMesh::VertexHandle center, MyMesh::VertexHandle pre, MyMesh::VertexHandle next);
float calc_area(MyMesh mesh, MyMesh::VertexHandle center, MyMesh::VertexHandle pre, MyMesh::VertexHandle next);
MyMesh calc_degree_area(MyMesh mesh);
//MyMesh laplace_smooth(MyMesh mesh, int iterations, float lamda);
//MyMesh implicit_laplace_smooth(MyMesh mesh, float lamda, int step);
//Eigen::VectorXd Conjugate(Eigen::VectorXd B, Eigen::VectorXd X, Eigen::MatrixXd A,int n);
//
//
//
//MyMesh::Color value_to_color(MyMesh::Scalar value, MyMesh::Scalar min, MyMesh::Scalar max);
//extern OpenMesh::VPropHandleT<OpenMesh::Vec3f>  vunicurvature1_, vcurvature1_, vgausscurvature1_;
//extern OpenMesh::VPropHandleT<MyMesh::Scalar>  vunicurvature_, vcurvature_, vgausscurvature_;
//extern OpenMesh::EPropHandleT<MyMesh::Scalar>  eweight_;
//extern OpenMesh::VPropHandleT<MyMesh::Scalar>  eweightSum_;

