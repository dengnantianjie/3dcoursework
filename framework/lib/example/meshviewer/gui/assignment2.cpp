#include "meshviewer.h"
#include "Algorithms.h"

#include "mh/base/defs.h"
#include "mh/base/imports.h"

#include "mh/3d/camera.h"
#include "mh/gui/sceneviewer.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "mh/gui/sceneviewer.h"
#include "mh/io/meshio.h"
#include <stdlib.h>

using namespace mh;
using namespace std;

std::shared_ptr<Mesh> mesh11,mesh22,mesh33,mesh44,mesh55,mesh66,mesh_reset;
std::vector<Eigen::Vector3f> vertData1,vertData2;
MyMesh om1,om_explicit_laplace_1,om_explicit_laplace_20,om_explicit_laplace_10;
MyMesh om_implicit_laplace_1,om_implicit_laplace_20,om_implicit_laplace_10;



#define Min(a, b) (((a) < (b)) ? (a) : (b))
#define Max(a, b) (((a) > (b)) ? (a) : (b))
#define EPSILON 1e-5
#define LARGE_VALUE 1e5
const float PI = 3.1415926;
using namespace OpenMesh;
typedef Eigen::Triplet<double> T;

void color_to_show(MyMesh &mesh, OpenMesh::VPropHandleT<MyMesh::Scalar> prop, std::shared_ptr<Mesh> mesh11){
    MyMesh::VertexIter  v_it;
    MyMesh::Scalar      curve, min(FLT_MAX), max(FLT_MAX);
    Eigen::Vector4f colors;

    std::vector<MyMesh::Scalar> scalar_v;
    scalar_v.reserve(mesh.n_vertices());
    for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
        scalar_v.push_back(mesh.property(prop, v_it));
    }
    unsigned int n = scalar_v.size() - 1;
    unsigned int discard = n / 5;
    std::sort(scalar_v.begin(), scalar_v.end());
    min = scalar_v[1 + discard];
    max = scalar_v[n-1 - discard];

    Mesh mesh_temp = Mesh(*mesh11);
    std::vector<Eigen::Vector3f> vertexData;
    vertexData = mesh_temp.getVertData();
    unsigned i_ = 0;
    for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {
        curve = mesh.property(prop, v_it);
        colors = calculate_color(curve, min, max);
        mesh_temp.getVertices().at(i_)->setColor(colors);
        ++i_;
    }
    swap(*mesh11,mesh_temp);

}

Eigen::Vector4f calculate_color(MyMesh::Scalar curve, MyMesh::Scalar min, MyMesh::Scalar max) {
    MyMesh::Scalar seg0,seg1, seg2, seg3, seg4,v0,v1,v2,v3,v4;
    Eigen::Vector4f col = Eigen::Vector4f(1.0f,1.0f,1.0f,0.3f);
    float R = 0;
    float G = 0;
    float B = 0;

    seg0 = min + 0.0 / 4.0 * (max - min);
    seg1 = min + 1.0 / 4.0 * (max - min);
    seg2 = min + 2.0 / 4.0 * (max - min);
    seg3 = min + 3.0 / 4.0 * (max - min);
    seg4 = min + 4.0 / 4.0 * (max - min);
    
    if (curve == 0){
        R = 0;
        G = 0;
        B = 0;
    }
    
    else{
        if (curve < seg0){ R = 0; G = 0; B = 255; }
        if (curve > seg4){ R = 255; G = 0; B = 0; }

        R = Max(0, (curve - seg2) / (seg3 - seg2));
        if (R > 1){ R = 1; }
        R = (float)(255.0 * R);

        G = (curve - seg0) / (seg1 - seg0);
        if (G > 3){ G = 1 - (curve - seg3) / (seg4 - seg3); }
        if (G > 1){ G = 1; }
        G = (float)(255 * G);


        B = Max(0, 1 - (curve - seg1) / (seg2 - seg1));
        if (B > 1){ B = 1; }
        B = (float)(255.0 * B);

    }

    col = Eigen::Vector4f(R/255, G/255, B/255 ,0.2f);
    
    return col;
}



MyMesh uniform_curvature(MyMesh mesh){
    MyMesh::VertexIter vIt;
    std::ofstream myfile("/Users/DNTJ/Desktop/3D_coursework/bunny/output/uniform.csv");
    for (vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt){
        Vec3f v(0, 0, 0);
        int num = 0;
        for (MyMesh::VertexVertexIter vv_it = mesh.vv_iter(vIt.handle()); vv_it; ++vv_it)
        {
            v = v + mesh.point(vv_it);
            num = num + 1;
        }
        v = v/ num;
        v = v - mesh.point(vIt);
        mesh.property(uniformcurvature, vIt) = v;
        Vec3f normal = mesh.normal(vIt);
        float cos_theta = calc_theta(normal, v);
        float normal_length = normal.length()*cos_theta;

        mesh.property(uniformcurvature_n, vIt) = (v.norm()/(- 2*normal_length));
        mesh.property(uniformcurvature_n_toshow, vIt) = abs(v.norm()/(- 2*normal_length));

        float output = (v.norm() / (-2 * normal_length));
        if (myfile.is_open())
        {
        	myfile << output <<","<<"\n";
            //myfile << normal <<","<<"\n";
        
        }
    }
    
    myfile.close();
    return mesh;
}

MyMesh gauss_curvature(MyMesh mesh){
    mesh = calc_degree_area(mesh);
    std::ofstream myfile("/Users/DNTJ/Desktop/3D_coursework/bunny/output/gauss.csv");
    MyMesh::VertexIter vIt;
    MyMesh::VertexVertexIter vvIt;
    for (vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt){
        float sum_degree = mesh.property(sumdegree, vIt.handle());
        float centric_sum_area = mesh.property(centerarea, vIt.handle());

        mesh.property(gausscurvature_n, vIt.handle()) = ((2 * PI - sum_degree)/centric_sum_area);
        mesh.property(gausscurvature_n_toshow, vIt.handle()) = abs((2 * PI - sum_degree)/centric_sum_area);

        float output = abs((2 * PI - sum_degree) / centric_sum_area);
        if (myfile.is_open())
        {
            myfile << output << "," << "\n";
        }
        
    }
    myfile.close();
    
    return mesh;
    
    
}

MyMesh lb_mean_curvature(MyMesh mesh_)
{
    mesh_ = calc_weights(mesh_);
    mesh_ = calc_degree_area(mesh_);
    std::ofstream myfile("/Users/DNTJ/Desktop/3D_coursework/bunny/output/mean.csv");
    for (MyMesh::VertexIter vIt = mesh_.vertices_begin(); vIt != mesh_.vertices_end(); ++vIt)
    {
        Vec3f L(0, 0, 0);
        float sumW = 0;
        float area;

        OpenMesh::HalfedgeHandle currentE, StartE;
        StartE = currentE = mesh_.halfedge_handle(vIt);
        Vec3f v = mesh_.point(vIt);

        do {
            Vec3f vi = mesh_.point(mesh_.to_vertex_handle(currentE));
            OpenMesh::EdgeHandle E = mesh_.edge_handle(currentE);
            float Weight = mesh_.property(eweight_, E);
            sumW += Weight;
            L = L + (Weight)*(vi - v);
            OpenMesh::HalfedgeHandle twinHEH = mesh_.opposite_halfedge_handle(currentE);
            currentE = mesh_.next_halfedge_handle(twinHEH);
        } while (currentE != StartE);

        area = mesh_.property(centerarea, vIt);
        mesh_.property(lbcurvature, vIt) = L;
        mesh_.property(lbcurvature_n, vIt) = L.norm()/(2*area);
        mesh_.property(lbcurvature_n_toshow, vIt) = abs(L.norm()/(2*area));
        mesh_.property(eweightSum_, vIt) = sumW;
        float output = L.norm() / (2 * area);
        if (myfile.is_open())
        {
            myfile << output << "," << "\n";
        }
    }
    myfile.close();
    return mesh_;
}

float calc_theta(Vec3f normal, Vec3f target_vector){
    float cos_theta = normal[0] * target_vector[0] + normal[1] * target_vector[1] + normal[2] * target_vector[2];
    cos_theta = cos_theta / (normal.length()*target_vector.length());
    return cos_theta;
}

MyMesh calc_weights(MyMesh mesh_)
{

    for (MyMesh::VertexIter vIt = mesh_.vertices_begin(); vIt != mesh_.vertices_end(); ++vIt)
    {
        OpenMesh::HalfedgeHandle currentE, StartE;

        StartE = currentE = mesh_.halfedge_handle(vIt);

        do {
            Vec3f v1 = mesh_.point(mesh_.from_vertex_handle(currentE));
            Vec3f v2 = mesh_.point(mesh_.to_vertex_handle(currentE));
            Vec3f cotangent_v_a = mesh_.point(mesh_.to_vertex_handle(mesh_.next_halfedge_handle(currentE)));
            Vec3f cotangent_v_b = mesh_.point(mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(currentE))));

            Vec3f a_v1 = (v1 - cotangent_v_a).normalize();
            Vec3f a_v2 = (v2 - cotangent_v_a).normalize();
            Vec3f b_v1 = (v1 - cotangent_v_b).normalize();
            Vec3f b_v2 = (v2 - cotangent_v_b).normalize();
            float cosA = a_v1 | a_v2;
            float cosB = b_v1 | b_v2;
            float sinA = sqrt(1 - cosA*cosA);
            float sinB = sqrt(1 - cosB*cosB);
            float cotA = sinA < EPSILON ? LARGE_VALUE : cosA / sinA;
            float cotB = sinB < EPSILON ? LARGE_VALUE : cosB / sinB;

            OpenMesh::EdgeHandle E = mesh_.edge_handle(currentE);
            mesh_.property(eweight_, E) = abs(cotA) + abs(cotB);
            OpenMesh::HalfedgeHandle twinHEH = mesh_.opposite_halfedge_handle(currentE);
            currentE = mesh_.next_halfedge_handle(twinHEH);
        } while (currentE != StartE);
    }
    return mesh_;

}

float calc_degree(MyMesh mesh, MyMesh::VertexHandle center, MyMesh::VertexHandle pre, MyMesh::VertexHandle next){
    MyMesh::Point o = mesh.point(center);
    MyMesh::Point a = mesh.point(pre);
    MyMesh::Point b = mesh.point(next);
    OpenMesh::Vec3f ao = o - a;
    OpenMesh::Vec3f bo = o - b;
    float dotProduct = dot(ao, bo) / (ao.length() * bo.length());
    float degree = acos(dotProduct);
    return degree;
}

float calc_area(MyMesh mesh, MyMesh::VertexHandle center, MyMesh::VertexHandle pre, MyMesh::VertexHandle next){
    MyMesh::Point o = mesh.point(center);
    MyMesh::Point a = mesh.point(pre);
    MyMesh::Point b = mesh.point(next);
    OpenMesh::Vec3f ao = o - a;
    OpenMesh::Vec3f bo = o - b;
    OpenMesh::Vec3f ab = a - b;
    float s = (ao.length() + bo.length() + ab.length()) / 2;
    float area = sqrt(s*(s - ao.length())*(s - bo.length())*(s - ab.length()));
    area = area / 3;
    return area;
}


MyMesh calc_degree_area(MyMesh mesh){
    std::ofstream myfile("/Users/DNTJ/Desktop/3D_coursework/bunny/output/area.csv");
    for (MyMesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt){
        float sum_degree = 0.0;
        float centric_sum_area = 0.0;
        std::vector<MyMesh::VertexHandle> neighbour_vertex;
        MyMesh::VertexHandle center = vIt.handle();

        MyMesh::VertexHandle start;
        MyMesh::VertexHandle prev;
        MyMesh::VertexHandle next;
        MyMesh::VertexHandle last;

        for (MyMesh::VertexVertexIter vvIt = mesh.vv_iter(vIt.handle()); vvIt; ++vvIt) {
            neighbour_vertex.push_back(vvIt.handle());
        }

        start = neighbour_vertex[0];
        neighbour_vertex.push_back(start);
        for (size_t vertex_id = 0; vertex_id < neighbour_vertex.size() - 1; ++vertex_id)
        {
            prev = neighbour_vertex[vertex_id];
            next = neighbour_vertex[vertex_id + 1];
            MyMesh::Point o = mesh.point(center);
            MyMesh::Point a = mesh.point(prev);
            MyMesh::Point b = mesh.point(next);
            OpenMesh::Vec3f ao = o - a;
            OpenMesh::Vec3f bo = o - b;
            OpenMesh::Vec3f ab = a - b;
            float dotProduct = dot(ao, bo) / (ao.length() * bo.length());
            float degree = acos(dotProduct);
            sum_degree += degree;
            float s = (ao.length() + bo.length() + ab.length()) / 2;
            float area = sqrt(s*(s - ao.length())*(s - bo.length())*(s - ab.length()));
            area = area / 3;
            centric_sum_area += area;
        }
        //std::cout << centric_sum_area << std::endl;
        if (myfile.is_open())
        {
            myfile << centric_sum_area << "," << "\n";
        }
        mesh.property(centerarea, vIt.handle()) = centric_sum_area;
        mesh.property(sumdegree, vIt.handle()) = sum_degree;
    }
    myfile.close();
    return mesh;
}


MyMesh get_principles(MyMesh mesh, OpenMesh::VPropHandleT<MyMesh::Scalar> h, OpenMesh::VPropHandleT<MyMesh::Scalar> k,int flag){
    MyMesh::VertexIter  v_it;
    MyMesh::Scalar mean_curve;
    MyMesh::Scalar gauss_curve;
    string path1;
    string path2;
    if (flag == 1){
         path1 = "/Users/DNTJ/Desktop/3D_coursework/bunny/output/k1_uni.csv";
         path2 = "/Users/DNTJ/Desktop/3D_coursework/bunny/output/k2_uni.csv";
    }

    if (flag == 2){
         path1 = "/Users/DNTJ/Desktop/3D_coursework/bunny/output/k1_lb.csv";
         path2 = "/Users/DNTJ/Desktop/3D_coursework/bunny/output/k2_lb.csv";
        
    }
    std::ofstream myfile1 (path1);
    std::ofstream myfile2 (path2);
    
    
    for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {
        mean_curve = mesh.property(h, v_it);
        gauss_curve = mesh.property(k, v_it);
        MyMesh::Scalar k1 = 0;
        MyMesh::Scalar k2 = 0;
        //::cout << mean_curve << std::endl;
        //std::cout << gauss_curve << std::endl;
        
        if (mean_curve*mean_curve - gauss_curve >= 0){
            k1 = mean_curve + sqrt(mean_curve*mean_curve - gauss_curve);
            k2 = mean_curve - sqrt(mean_curve*mean_curve - gauss_curve);
        }
        else{
            k1 = 0;
            k2 = 0;
        }
        if (myfile1.is_open())
        {
            myfile1 << k1 << "," << "\n";
        }
        if (myfile2.is_open())
        {
            myfile2 << k2 << "," << "\n";
        }
        if (flag == 1){
            mesh.property(principle_max_uni, v_it) = k1;
            mesh.property(principle_min_uni, v_it) = k2;
            mesh.property(principle_max_uni_toshow, v_it) = abs(k1);
            mesh.property(principle_min_uni_toshow, v_it) = abs(k2);
        }
        if (flag == 2){
            mesh.property(principle_max_lb, v_it) = k1;
            mesh.property(principle_min_lb, v_it) = k2;
            mesh.property(principle_max_lb_toshow, v_it) = abs(k1);
            mesh.property(principle_min_lb_toshow, v_it) = abs(k2);

        }
        

    }
    myfile1.close();
    myfile2.close();
    return mesh;

}

MyMesh explicit_laplace_smooth(MyMesh mesh, std::shared_ptr<Mesh> mesh11, int iterations, double lamda){
    MyMesh origin = mesh;
    for (int i = 0; i < iterations; i++){
        std::vector<MyMesh::Point> after_laplace;
        
        for (MyMesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt){
            MyMesh::Point p_new;
            OpenMesh::Vec3f laplace_operator;
            MyMesh::Point laplace;
            laplace_operator = mesh.property(uniformcurvature, vIt);
            p_new = mesh.point(vIt.handle()) + laplace_operator * lamda;
            after_laplace.push_back(p_new);
        }
        int k = 0;
        for (MyMesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt, ++k) {
            //std::cout << k << std::endl;
            MyMesh::Point p_new_laplace = after_laplace[k];
            mesh.set_point(vIt.handle(), p_new_laplace);
        }

        mesh.update_normals();
        mesh = uniform_curvature(mesh);
        mesh.update_normals();
    }
    float sum_v = 0;
    for (MyMesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt){
        MyMesh::Point p_1 = mesh.point(vIt.handle());
        MyMesh::Point p_2 = origin.point(vIt.handle());
        float sum = ((p_1[0] - p_2[0])*(p_1[0] - p_2[0]) + (p_1[1] - p_2[1])*(p_1[1] - p_2[1]) + (p_1[2] - p_2[2])*(p_1[2] - p_2[2]));
        sum = sqrt(sum);
        sum_v = sum_v + sum;
    }
    std::ofstream myfile("/Users/DNTJ/Desktop/3D_coursework/bunny/output/explicitlaplace.csv");
    sum_v = sum_v / mesh.n_vertices();
    myfile << sum_v <<","<<"\n";

    Mesh mesh_temp20 = transform_structure(mesh);
    swap(*mesh11,mesh_temp20);
    
    return mesh;
}

MyMesh implicit_laplace_smooth(MyMesh mesh,  std::shared_ptr<Mesh> mesh11, int step, double lamda){

    MyMesh origin = mesh;
    mesh = calc_weights(mesh);
    mesh = calc_degree_area(mesh);
    mesh = lb_mean_curvature(mesh);
    int n = mesh.n_vertices();
    Eigen::MatrixXd A_dense = Eigen::MatrixXd::Zero(n, n);
    Eigen::SparseMatrix<double> A(n, n);
    Eigen::VectorXd B0 = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd B1 = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd B2 = Eigen::VectorXd::Zero(n);

    Eigen::VectorXd X0 = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd X1 = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd X2 = Eigen::VectorXd::Zero(n);

    std::vector<T> coefficients;

    for (MyMesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt)
    {
        float sumW = 0;
        Vec3f v = mesh.point(vIt);

        OpenMesh::HalfedgeHandle currentE, StartE;
        StartE = currentE = mesh.halfedge_handle(vIt);
        size_t i_index = vIt.handle().idx();
        do {
            MyMesh::VertexHandle vi = mesh.to_vertex_handle(currentE);
            OpenMesh::EdgeHandle E = mesh.edge_handle(currentE);
            float Weight = mesh.property(eweight_, E);
            size_t j_index = vi.idx();
            
            A_dense(i_index, j_index) =  -lamda* Weight;
            OpenMesh::HalfedgeHandle twinHEH = mesh.opposite_halfedge_handle(currentE);
            currentE = mesh.next_halfedge_handle(twinHEH);
        } while (currentE != StartE);

        float area = mesh.property(centerarea, vIt);
        sumW = mesh.property(eweightSum_, vIt);
        A_dense(i_index, i_index) = (2 * area) + lamda * sumW;
        
        
        B0(i_index) = v[0] * (2 * area);
        B1(i_index) = v[1] * (2 * area);
        B2(i_index) = v[2] * (2 * area);
    }

    for (int j = 0; j<n; ++j)
    {
        for (int i = 0; i<n; ++i)
        {
            if (A_dense(j, i)*A_dense(j, i)>0) {
                //std::cout << A_dense(j, i) << std::endl;
                coefficients.push_back(T(j, i, A_dense(j, i)));
            }
        }
    }

    
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver0;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver1;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver2;

    solver0.compute(A);
    solver1.compute(A);
    solver2.compute(A);

    for (int count = 0; count < step; count++) 
    { 
        
        X0 = solver0.solve(B0); 
        
    }
    for (int count = 0; count < step; count++)
    {
        
        X1 = solver1.solve(B1);
        
    }
    for (int count = 0; count < step; count++) {
        
        X2 = solver2.solve(B2);
        
    }
    
    for (MyMesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) {
        MyMesh::Point p_new;
        size_t idx = vIt.handle().idx();
        p_new[0] = X0(idx);
        p_new[1] = X1(idx);
        p_new[2] = X2(idx);
        mesh.set_point(vIt.handle(), p_new);
        
    }
    mesh.update_normals();
    OpenMesh::IO::write_mesh(mesh, "/Users/DNTJ/Desktop/3D_coursework/bunny/data/dragon_update.obj");
    float sum_v = 0;
    for (MyMesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt){
        MyMesh::Point p_1 = mesh.point(vIt.handle());
        MyMesh::Point p_2 = origin.point(vIt.handle());
        float sum = ((p_1[0] - p_2[0])*(p_1[0] - p_2[0]) + (p_1[1] - p_2[1])*(p_1[1] - p_2[1]) + (p_1[2] - p_2[2])*(p_1[2] - p_2[2]));
        sum = sqrt(sum);
        sum_v = sum_v + sum;
    }
    std::ofstream myfile("/Users/DNTJ/Desktop/3D_coursework/bunny/output/implicitlaplace.csv");
    sum_v = sum_v / mesh.n_vertices();
    myfile << sum_v <<","<<"\n";

    Mesh mesh_temp20 = transform_structure(mesh);
    swap(*mesh11,mesh_temp20);

    return mesh;
    
}

