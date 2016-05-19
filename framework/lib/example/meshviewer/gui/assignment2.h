#include <iostream>
#include <stdio.h>
#include <random>
#include <math.h>
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


//extern OpenMesh::VPropHandleT<OpenMesh::Vec3f>  uniformcurvature;
//extern OpenMesh::VPropHandleT<OpenMesh::TriMesh_ArrayKernelT<>::Scalar>  uniformcurvature_n;
//extern OpenMesh::VPropHandleT<OpenMesh::Vec3f> meancurvature;
//extern OpenMesh::VPropHandleT<MyMesh::Scalar> meancurvature_n;
//extern OpenMesh::EPropHandleT<MyMesh::Scalar>  weight;
//extern OpenMesh::VPropHandleT<OpenMesh::Vec3f> gausscurvature;
//extern OpenMesh::VPropHandleT<MyMesh::Scalar> gausscurvature_n;
//extern OpenMesh::VPropHandleT<MyMesh::Scalar> principle_max;
//extern OpenMesh::VPropHandleT<MyMesh::Scalar> principle_min;
//extern OpenMesh::VPropHandleT<float> centerarea;
//extern OpenMesh::VPropHandleT<float> sumdegree;


OpenMesh::TriMesh_ArrayKernelT<> uniform_curvature(OpenMesh::TriMesh_ArrayKernelT<> mesh);
//MyMesh mean_curvature(MyMesh mesh);
//MyMesh gauss_curvature(MyMesh mesh);
//MyMesh calc_degree_area(MyMesh mesh);
//
//
//float calc_area(MyMesh mesh,MyMesh::VertexHandle center, MyMesh::VertexHandle pre, MyMesh::VertexHandle next);
//float calc_degree(MyMesh mesh,MyMesh::VertexHandle center, MyMesh::VertexHandle pre, MyMesh::VertexHandle next);
//
//void color_code(MyMesh &mesh, OpenMesh::VPropHandleT<MyMesh::Scalar> prop);
//OpenMesh::Vec3uc getcolor(MyMesh::Scalar curve, MyMesh::Scalar min, MyMesh::Scalar max);
//MyMesh calc_mean_curvature(MyMesh mesh_);
//MyMesh calc_weights(MyMesh mesh_);
//MyMesh calc_principle(MyMesh mesh, OpenMesh::VPropHandleT<MyMesh::Scalar> h, OpenMesh::VPropHandleT<MyMesh::Scalar> k);
float calc_theta(OpenMesh::Vec3f normal, OpenMesh::Vec3f target_vector);
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
