//
// Created by bontius on 25/01/16.
//

#include "gui/meshviewer.h"
#include "gui/Algorithms.h"

#include "mh/base/defs.h"
#include "mh/base/imports.h"

#include "mh/3d/camera.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "mh/gui/sceneviewer.h"
#include "mh/io/meshio.h"
#include <stdlib.h>

using namespace mh;
using namespace std;

double thresh = 0.0001;
int k = 1;
int dim = 3;
double eps = 0;

//std::shared_ptr<Mesh> mesh11,mesh22,mesh33,mesh44,mesh55,mesh66;
//std::vector<Eigen::Vector3f> vertData1,vertData2;
//MyMesh om1;


//#define Min(a, b) (((a) < (b)) ? (a) : (b))
//#define Max(a, b) (((a) > (b)) ? (a) : (b))
//#define EPSILON 1e-5
//#define LARGE_VALUE 1e5
//const float PI = 3.1415926;
//using namespace OpenMesh;
typedef Eigen::Triplet<double> T;


OpenMesh::VPropHandleT<OpenMesh::Vec3f> uniformcurvature;
OpenMesh::VPropHandleT<MyMesh::Scalar> uniformcurvature_n;
OpenMesh::VPropHandleT<MyMesh::Scalar> uniformcurvature_n_toshow;

OpenMesh::EPropHandleT<MyMesh::Scalar>  weight;
OpenMesh::VPropHandleT<OpenMesh::Vec3f> gausscurvature;
OpenMesh::VPropHandleT<MyMesh::Scalar> gausscurvature_n;
OpenMesh::VPropHandleT<MyMesh::Scalar> gausscurvature_n_toshow;

OpenMesh::VPropHandleT<OpenMesh::Vec3f>  lbcurvature;
OpenMesh::VPropHandleT<MyMesh::Scalar>  lbcurvature_n;
OpenMesh::VPropHandleT<MyMesh::Scalar>  lbcurvature_n_toshow;

OpenMesh::VPropHandleT<MyMesh::Scalar> principle_max_uni;
OpenMesh::VPropHandleT<MyMesh::Scalar> principle_min_uni;
OpenMesh::VPropHandleT<MyMesh::Scalar> principle_max_lb;
OpenMesh::VPropHandleT<MyMesh::Scalar> principle_min_lb;

OpenMesh::VPropHandleT<MyMesh::Scalar> principle_max_uni_toshow;
OpenMesh::VPropHandleT<MyMesh::Scalar> principle_min_uni_toshow;
OpenMesh::VPropHandleT<MyMesh::Scalar> principle_max_lb_toshow;
OpenMesh::VPropHandleT<MyMesh::Scalar> principle_min_lb_toshow;

OpenMesh::VPropHandleT<float> centerarea;
OpenMesh::VPropHandleT<float> sumdegree;
OpenMesh::EPropHandleT<MyMesh::Scalar>  eweight_;
OpenMesh::VPropHandleT<MyMesh::Scalar>  eweightSum_;



inline void print_usage( int argc, char** argv )
{
    std::cout << "Usage: " << argv[0] << " first.ply second.ply ..." << std::endl;
}

int main(int argc, char* argv[])
{
    MeshViewer viewer;
    string meshpath;
    meshpath = "/Users/DNTJ/Desktop/3D_coursework/bunny/data/casting.obj";
    
    viewer.loadMesh(mesh11,meshpath);
    Eigen::Vector4f color1 = Eigen::Vector4f(1.0f,1.0f,1.0f,0.3f);
    setmeshcolor(mesh11,color1);
    
    OpenMesh::IO::Options opt;
    opt += OpenMesh::IO::Options::VertexColor;
    opt += OpenMesh::IO::Options::VertexNormal;
    opt += OpenMesh::IO::Options::VertexTexCoord;
    opt += OpenMesh::IO::Options::FaceColor;
    opt += OpenMesh::IO::Options::FaceNormal;
    opt += OpenMesh::IO::Options::FaceTexCoord;
    om1.request_face_normals();
    om1.request_face_colors();
    om1.request_vertex_normals();
    om1.request_vertex_colors();
   
    om1.add_property(uniformcurvature);
    om1.add_property(uniformcurvature_n);
    om1.add_property(uniformcurvature_n_toshow);
    om1.add_property(weight);
    om1.add_property(gausscurvature);
    om1.add_property(gausscurvature_n);
    om1.add_property(gausscurvature_n_toshow);
    om1.add_property(principle_max_uni);
    om1.add_property(principle_min_uni);
    om1.add_property(principle_max_lb);
    om1.add_property(principle_min_lb);
    om1.add_property(principle_max_uni_toshow);
    om1.add_property(principle_min_uni_toshow);
    om1.add_property(principle_max_lb_toshow);
    om1.add_property(principle_min_lb_toshow);
    om1.add_property(centerarea);
    om1.add_property(sumdegree);


    om1.add_property(lbcurvature);
    om1.add_property(lbcurvature_n);
    om1.add_property(lbcurvature_n_toshow);
    om1.add_property(eweight_);
    om1.add_property(eweightSum_);
    
    if (!OpenMesh::IO::read_mesh(om1, meshpath, opt))
    {
        std::cerr << "read error\n";
        exit(1);
    }
    if(! opt.check(OpenMesh::IO::Options::VertexNormal))
    {
        om1.update_vertex_normals();
    }

///////////////   part 1 Discrete Curvature ///////////////////////////
    om1 = uniform_curvature(om1);
    om1 = gauss_curvature(om1);
    om1 = lb_mean_curvature(om1);

    om1 = get_principles(om1, uniformcurvature_n, gausscurvature_n,1);
    om1 = get_principles(om1, lbcurvature_n, gausscurvature_n,2);

///////////////   part 1 Discrete Curvature ///////////////////////////
    om_explicit_laplace_20 = explicit_laplace_smooth(om1, mesh11,20, 0.5) ;  //////////// iteration 20 times
    om_explicit_laplace_10 = explicit_laplace_smooth(om1, mesh11,10, 0.5) ; //////////// iteration 10 times


///////////////   part 2 Laplacian Mesh Smoothing ///////////////////////////



///////////////   part 2 Laplacian Mesh Smoothing ///////////////////////////


    std::ofstream myfile("/Users/DNTJ/Desktop/3D_coursework/bunny/output/test_openmesh.csv");
    for (MyMesh::VertexIter v_it=om1.vertices_begin(); v_it!=om1.vertices_end(); ++v_it) {
        if (myfile.is_open())
        {
            myfile << om1.point(v_it) <<"\n";
        }
    }


    
    //Eigen::Vector4f color1 = Eigen::Vector4f(1.0f,0.0f,0.0f,0.3f);

    Mesh mesh_temp1 = transform_structure(om1);
    swap(*mesh11,mesh_temp1);


    while (!viewer.shouldClose())
    {
        viewer.makeCurrent();
        viewer.mainLoop();
        viewer.swap();
    }
    
    glfwTerminate();
    
    return 0;
}

//set color
void setmeshcolor(std::shared_ptr<Mesh> mesh11,Eigen::Vector4f colors){
    Mesh mesh_temp = Mesh(*mesh11);
    std::vector<Eigen::Vector3f> vertexData;
    vertexData = mesh_temp.getVertData();
    std::ofstream myfile("/Users/DNTJ/Desktop/3D_coursework/bunny/output/test_mesh.csv");
    unsigned i_ = 0;
    for (Eigen::Vector3f point: vertexData){
        if (myfile.is_open())
        {
            myfile << point.transpose() << "\n";
            //myfile << normal <<","<<"\n";
        }
    }
    swap(*mesh11,mesh_temp);
}


void readPt(int downsample, Mesh &mesh, ANNpointArray &p){
    int i = 0;
    int count = 0;
    std::vector<Eigen::Vector3f> vertexData;
    
    vertexData = mesh.getVertData();
    
    for (Eigen::Vector3f point: vertexData)
    {
        count++;
        if (count == downsample)
        {
            p[i][0] = point[0];
            p[i][1] = point[1];
            p[i][2] = point[2];
            i++;
            count = 0;
        }
    }
}

double ICPAlgo(std::shared_ptr<Mesh>& mesh11,std::shared_ptr<Mesh>& mesh22, int downsample){
    
    ANNpoint sourcePT;
    ANNpoint targetPT;
    ANNpoint newsourcePT;
    ANNpointArray ann_source;
    ANNpointArray ann_target;
    ANNpointArray ann_newsource;
    ANNidxArray nnindex;
    ANNdistArray dists;
    ANNkd_tree* kdTree;
    sourcePT = annAllocPt(dim);
    targetPT = annAllocPt(dim);
    nnindex = new ANNidx[k];
    dists = new ANNdist[k];
    
    double mean_newsource[3] = {};
    double mean_source[3] = {};
    double sum_dis = 0.0;
    
    Mesh mesh_source = Mesh(*mesh11);
    Mesh mesh_target = Mesh(*mesh22);
    
    int num_meshsource = mesh_source.nVerts();
    int num_meshtarget = mesh_target.nVerts();
    int num_sampled_source = num_meshsource / downsample;
    int num_sampled_target = num_meshtarget / downsample;
    
    
    ann_source = annAllocPts(num_sampled_source, dim);
    ann_target = annAllocPts(num_sampled_target, dim);
    ann_newsource = annAllocPts(num_sampled_source, dim);
    
    readPt(downsample, mesh_source, ann_source);
    readPt(downsample, mesh_target, ann_target);
    
    kdTree = new ANNkd_tree(
                            ann_target,
                            num_sampled_target,
                            dim);    //set up kd-tree for target mesh
    for (int n_s = 0; n_s < num_sampled_source; n_s++){
        sourcePT = ann_source[n_s];
        kdTree->annkSearch(
                           sourcePT,
                           k,
                           nnindex,
                           dists,
                           eps);
        sum_dis += *dists/num_sampled_source;
        newsourcePT = ann_target[*nnindex];
        ann_newsource[n_s] = newsourcePT;
        for (int d = 0; d < dim; d++){
            mean_source[d] += sourcePT[d];
            mean_newsource[d] += newsourcePT[d];
            //sum of x,y,z
        }
    }
    
    for (int d = 0; d < dim; d++){
        double x = mean_source[d];
        double y = mean_newsource[d];
        mean_source[d] = x / num_sampled_source;
        mean_newsource[d] = y / num_sampled_source;
    }//mean of source and newsource
    
    
    gsl_matrix * h = gsl_matrix_alloc(dim, dim);
    double inter_diff_array[3][3] = {};
    
    
    for (int n = 0; n < num_sampled_source; n++)
    {
        double diff_source[3] = {};
        double diff_newsource[3] = {};
        sourcePT = ann_source[n];
        newsourcePT = ann_newsource[n];
        for (int d = 0; d<dim; d++)
        {
            diff_newsource[d] = newsourcePT[d] - mean_newsource[d];
            diff_source[d] = sourcePT[d] - mean_source[d];
        }
        
        for (int i = 0; i < dim; i++){
            for (int j = 0; j < dim; j++){
                inter_diff_array[i][j] += diff_newsource[i] * diff_source[j];
            }
        }
        
    }//sum of diff
    
    for (int ii = 0; ii < 3; ii++){
        for (int jj = 0; jj < 3; jj++){
            gsl_matrix_set(h, ii, jj, inter_diff_array[ii][jj]);
            
        }
    }
    delete[] nnindex;
    delete[] dists;
    delete kdTree;
    annClose();
    
    //gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0, diff_source, diff_newsource,0.0, h);
    //std::cout << "start to calculate Rotation and Translation Matrix" << std::endl;
    gsl_matrix * Rotation = gsl_matrix_alloc(dim, dim);
    gsl_matrix * translation = gsl_matrix_alloc(dim, 1);;
    
    gsl_matrix * V = gsl_matrix_alloc(dim, dim);
    gsl_vector * S = gsl_vector_alloc(dim);
    gsl_vector * work = gsl_vector_alloc(dim);
    gsl_matrix * u = gsl_matrix_alloc(dim, dim);
    gsl_matrix *center_source = gsl_matrix_alloc(dim, 1);
    gsl_matrix *center_newsource = gsl_matrix_alloc(dim, 1);
    gsl_linalg_SV_decomp(h, V, S, work);//
    
    
    for (int row = 0; row<dim; row++)
    {
        gsl_matrix_set(center_newsource, row, 0, mean_newsource[row]);
        gsl_matrix_set(center_source, row, 0, mean_source[row]);
    }
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, h, V, 0, Rotation);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Rotation, center_source, 0, translation);
    gsl_matrix_sub(center_newsource, translation);
    translation = center_newsource;
    
    std::vector<Eigen::Vector3f> vertexData;
    vertexData = mesh_source.getVertData();
    unsigned i_ = 0;
    
    for (Eigen::Vector3f point: vertexData)
    {
        
        //std::cout<<point[0]<<","<<point[1]<<","<<point[2]<<std::endl;
        gsl_matrix *sPt = gsl_matrix_alloc(dim, 1);
        gsl_matrix *nsPt = gsl_matrix_alloc(dim, 1);
        
        for (int d = 0; d<dim; d++)
            gsl_matrix_set(sPt, d, 0, point[d]);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Rotation, sPt, 1.0, nsPt);
        Eigen::Vector3f newpoint;
        for (int d = 0; d < dim; d++)
            newpoint[d] = gsl_matrix_get(nsPt, d, 0) + gsl_matrix_get(translation, d, 0);
        //std::cout<<"new    "<<newpoint[0]<<","<<newpoint[1]<<","<<newpoint[2]<<std::endl;
        //std::count<<mesh_source.getVertices().at(i_)<<std::endl;
        mesh_source.getVertices().at(i_)->setPosition(newpoint);
        //std::count<<"new  "<<mesh_source.getVertices().at(i_)<<std::endl;
        gsl_matrix_free(sPt);
        gsl_matrix_free(nsPt);
        ++i_;
    }// get new point position
    gsl_matrix_free(h);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_matrix_free(u);
    gsl_matrix_free(Rotation);
    gsl_matrix_free(translation);
    swap(*mesh11,mesh_source);
    swap(*mesh22,mesh_target);
    return sum_dis;
    
    
    
}

void rotateMesh(std::shared_ptr<Mesh>& mesh11,double degree, double x, double y, double z)
{
    gsl_matrix *rotation = gsl_matrix_alloc(dim, dim);
    Mesh mesh = Mesh(*mesh11);
    double translation[3] = {};
    double Pt_mean[3] = {};
    double Pt_sum[3] = {};
    
    double r[3][3] = {};
    double c = cos(degree*3.1415926 / 180);
    double s = sin(degree*3.1415926 / 180);
    
    r[0][0] = c + (1 - c)*x*x;
    r[0][1] = (1 - c)*x*y - s*z;
    r[0][2] = (1 - c)*x*z + s*y;
    r[1][0] = (1 - c)*y*x + s*z;
    r[1][1] = c + (1 - c)*y*y;
    r[1][2] = (1 - c)*y*z - s*x;
    r[2][0] = (1 - c)*z*x - s*y;
    r[2][1] = (1 - c)*z*y + s*x;
    r[2][2] = c + (1 - c)*z*z;
    
    for (int row = 0; row<dim; row++)
    {
        for (int col = 0; col<dim; col++)
        {
            gsl_matrix_set(rotation, row, col, r[row][col]);
        }
    }
    std::vector<Eigen::Vector3f> vertexData;
    vertexData = mesh.getVertData();
    
    for (Eigen::Vector3f point: vertexData)
    {
        double Pt[3] = {};
        for (int d = 0; d<dim; d++)
        {
            Pt[d] += point[d];
            Pt_sum[d] += Pt[d];
        }
    }
    
    for (int d = 0; d<dim; d++)
    {
        Pt_mean[d] = Pt_sum[d] / mesh.nVerts();
        
        translation[d] = Pt_mean[d];
    }
    int i_ = 0;
    for (Eigen::Vector3f point: vertexData)
    {
        gsl_matrix *oldPt = gsl_matrix_alloc(dim, 1);
        gsl_matrix *newPt = gsl_matrix_alloc(dim, 1);
        double oldPtarray[3] = {};
        for (int d = 0; d<dim; d++)
        {
            oldPtarray[d] = point[d];
            oldPtarray[d] -= translation[d];
            gsl_matrix_set(oldPt, d, 0, oldPtarray[d]);
        }
        
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rotation, oldPt, 1.0, newPt);
        Eigen::Vector3f newpoint;
        for (int d = 0; d<dim; d++)
            newpoint[d] = gsl_matrix_get(newPt, d, 0);
        
        mesh.getVertices().at(i_)->setPosition(newpoint);
        gsl_matrix_free(oldPt);
        gsl_matrix_free(newPt);
        i_ ++;
    }
    swap(*mesh11,mesh);
}



/*void translateMesh(std::shared_ptr<Mesh>& mesh11)
 {
 Mesh mesh = Mesh(*mesh11);
 double translation[3] = {};
 double Pt_mean[3] = {};
 double Pt_sum[3] = {};
 
 std::vector<Eigen::Vector3f> vertexData;
 vertexData = mesh.getVertData();
 
 for (Eigen::Vector3f point: vertexData)
 {
 double Pt[3] = {};
 for (int d = 0; d<dim; d++)
 {
 Pt[d] += point[d];
 Pt_sum[d] += Pt[d];
 }
 }
 
 for (int d = 0; d<dim; d++)
 {
 Pt_mean[d] = Pt_sum[d] / mesh.nVerts();
 
 translation[d] = Pt_mean[d];
 }
 int i_ = 0;
 for (Eigen::Vector3f point: vertexData)
 {
 gsl_matrix *oldPt = gsl_matrix_alloc(dim, 1);
 gsl_matrix *newPt = gsl_matrix_alloc(dim, 1);
 double oldPtarray[3] = {};
 for (int d = 0; d<dim; d++)
 {
 oldPtarray[d] = point[d];
 oldPtarray[d] -= translation[d];
 gsl_matrix_set(oldPt, d, 0, oldPtarray[d]);
 }
 
 Eigen::Vector3f newpoint;
 for (int d = 0; d<dim; d++)
 newpoint[d] = gsl_matrix_get(newPt, d, 0);
 mesh.getVertices().at(i_)->setPosition(newpoint);
 mesh.setPosition(newpoint);
 gsl_matrix_free(oldPt);
 gsl_matrix_free(newPt);
 ++ i_;
 }
 swap(*mesh11,mesh);
 }
 */





void meshAlign(int ratio)
{
    
    double err = ICPAlgo(mesh22, mesh11,ratio);
    std::cout << "Mean_Distance is: " << err << std::endl;
    
}




void meshAlignAll(int ratio)
{
    
    for (int i = 0; i < 25; i++){
        
        double err1 = ICPAlgo(mesh22, mesh11,ratio);
        std::cout << "Mean_Distance is: " << err1 << std::endl;
    }
    std::cout << "finish" <<std::endl;
    for (int i = 0; i < 25; i++){
        double err2 = ICPAlgo(mesh33, mesh22, ratio);
        std::cout << "Mean_Distance is: " << err2 << std::endl;
    }
    std::cout << "finish" <<std::endl;
    for (int i = 0; i < 25; i++){
        double err3 = ICPAlgo(mesh44, mesh33, ratio);
        std::cout << "Mean_Distance is: " << err3 << std::endl;
    }
    std::cout << "finish" <<std::endl;
    for (int i = 0; i < 25; i++){
        double err4 = ICPAlgo(mesh55, mesh44, ratio);
        std::cout << "Mean_Distance is: " << err4 << std::endl;
    }
    std::cout << "finish" <<std::endl;
    
    for (int i = 0; i < 25; i++){
        double err5 = ICPAlgo(mesh66, mesh55, ratio);
        std::cout << "Mean_Distance is: " << err5 << std::endl;
    }
    std::cout << "finish" <<std::endl;
    
}




void addNoise(double sd)
{
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0, sd);
    Mesh mesh = Mesh(*mesh22);
    std::vector<Eigen::Vector3f> vertexData;
    
    vertexData = mesh.getVertData();
    int i_ = 0;
    for (Eigen::Vector3f point: vertexData)
    {
        Eigen::Vector3f newpoint;
        for (int d = 0; d<3; d++)
        {
            double noise = distribution(generator);
            //std::cout<<noise<<std::endl;
            //if ((noise >= -0.5) && (noise <= 0.5))
            
            newpoint[d] = point[d]*(1 + noise/30);
            //else
            //    newpoint[d] = point[d];
        }
        mesh.getVertices().at(i_)->setPosition(newpoint);
        ++ i_;
    }
    swap(*mesh22,mesh);
    
}



//  Q7
void setNormal(std::shared_ptr<Mesh> mesh11){
    int kk = 10;
    Mesh mesh_source = Mesh(*mesh11);
    std::vector<Eigen::Vector3f> vertexData;
    vertexData = mesh_source.getVertData();
    
    ANNpoint sourcePT;
    ANNpoint newsourcePT;
    ANNpointArray ann_source;
    ANNpointArray ann_newsource;
    ANNidxArray nnindex;
    ANNdistArray dists;
    ANNkd_tree* kdTree;
    sourcePT = annAllocPt(dim);
    nnindex = new ANNidx[kk];
    dists = new ANNdist[kk];
    
    
    int num_sampled_source = mesh_source.nVerts();
    
    ann_source = annAllocPts(num_sampled_source, dim);
    ann_newsource = annAllocPts(num_sampled_source, dim);
    
    readPt(1, mesh_source, ann_source);
    
    kdTree = new ANNkd_tree(
                            ann_source,
                            num_sampled_source,
                            dim);
    
    for (int n_s = 0; n_s < num_sampled_source; n_s++){
        sourcePT = ann_source[n_s];
        kdTree->annkSearch(
                           sourcePT,
                           kk,
                           nnindex,
                           dists,
                           eps);
        
        double data[9] = {};
        data[0] = 0;
        data[1] = 0;
        data[2] = 0;
        data[3] = 0;
        data[4] = 0;
        data[5] = 0;
        data[6] = 0;
        data[7] = 0;
        data[8] = 0;
        
        for (int i = 0; i < kk; i++){
            newsourcePT = ann_source[nnindex[i]];
            double x = newsourcePT[0] - sourcePT[0];
            double y = newsourcePT[1] - sourcePT[1];
            double z = newsourcePT[2] - sourcePT[2];
            
            data[0]+= x*x;
            data[1]+= x*y;
            data[2]+= x*z;
            data[3]+= y*x;
            data[4]+= y*y;
            data[5]+= y*z;
            data[6]+= z*x;
            data[7]+= z*y;
            data[8]+= z*z;
        }
        
        
        // eigen value
        gsl_matrix_view m
        = gsl_matrix_view_array (data, 3, 3);
        
        gsl_vector *eval = gsl_vector_alloc (3);
        gsl_matrix *evec = gsl_matrix_alloc (3, 3);
        
        gsl_eigen_symmv_workspace * w =
        gsl_eigen_symmv_alloc (3);
        
        gsl_eigen_symmv (&m.matrix, eval, evec, w);
        
        gsl_eigen_symmv_free (w);
        
        gsl_eigen_symmv_sort (eval, evec,
                              GSL_EIGEN_SORT_ABS_ASC);
        
        {
            int i;
            Eigen::Vector3f newpoint;
            double min_eval_i = gsl_vector_get (eval, 0);
            gsl_vector_view min_evec_i = gsl_matrix_column (evec, 0);
            for (i = 0; i < 3; i++)
            {
                double eval_i
                = gsl_vector_get (eval, i);
                gsl_vector_view evec_i
                = gsl_matrix_column (evec, i);
                if (eval_i < min_eval_i){
                    min_evec_i
                    = evec_i;
                    min_eval_i = eval_i;
                }
                
                //printf ("eigenvalue = %g\n", eval_i);
                //printf ("eigenvector = \n");
                //gsl_vector_fprintf (stdout,
                //                    &evec_i.vector, "%g");
            }
            
            //std::cout<<gsl_vector_get(&min_evec_i.vector,1)<<std::endl;
            //gsl_vector_fprintf (stdout,
            //                                       &min_evec_i.vector, "%g");
            for (i = 0; i < 3; i ++){
                newpoint[i] = gsl_vector_get(&min_evec_i.vector,i);
            }
            //newpoint = - newpoint;
            //if(newpoint[1] * 1 > 0){
            
            //   newpoint = - newpoint;
            //}
            
            gsl_vector_free (eval);
            gsl_matrix_free (evec);
            mesh_source.getVertices().at(n_s)->setNormal(newpoint.cast<float>());
            
        }
        
        
        
        
        
    }
    
    delete[] nnindex;
    delete[] dists;
    delete kdTree;
    annClose();
    
    swap(*mesh11,mesh_source);
    
}

/////////////////////////////////////////////////////////Assignment 2 part ///////////////////////////



Mesh transform_structure(MyMesh om1)
    {
        std::vector<Eigen::Vector3f> vertData;
        std::vector<Eigen::Vector3f> normalData;
        std::vector<Eigen::Vector3i> faceData;

        for (MyMesh::VertexIter v_it=om1.vertices_begin(); v_it!=om1.vertices_end(); ++v_it)
        {

            vertData.push_back(Eigen::Vector3f(om1.point(v_it)[0],om1.point(v_it)[1],om1.point(v_it)[2]));
            normalData.push_back(Eigen::Vector3f(om1.normal(v_it)[0],om1.normal(v_it)[1],om1.normal(v_it)[2]));
            
        }
        for(MyMesh::FaceIter f_it = om1.faces_begin(); f_it != om1.faces_end(); ++f_it) 
        {    
            std::vector<int> faceint;
            for (MyMesh::FaceVertexIter fv_it = om1.fv_iter(f_it); fv_it; ++fv_it)
            {
                 faceint.push_back(fv_it.handle().idx());
            }
            faceData.push_back(Eigen::Vector3i(faceint[0],faceint[1],faceint[2]));
        }
       

        Mesh mesh(vertData, normalData, faceData);

        return mesh;
    } //...convertFromAssimp()

void setmeshcolor_new(std::shared_ptr<Mesh> mesh11,Eigen::Vector4f colors){
    Mesh mesh_temp = Mesh(*mesh11);
    std::vector<Eigen::Vector3f> vertexData;
    vertexData = mesh_temp.getVertData();
    
    unsigned i_ = 0;
    for (Eigen::Vector3f point: vertexData){
        mesh_temp.getVertices().at(i_)->setColor(colors);
        ++i_;
    }
    swap(*mesh11,mesh_temp);
}

void reset_meshs(std::shared_ptr<Mesh> mesh11,MyMesh om1){



    Mesh mesh_temp1 = transform_structure(om1);
    swap(*mesh11,mesh_temp1);
}






    











