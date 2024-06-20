
#define __FUNCSIG__ __func__ //替代Visual Studio的获取函数签名的宏

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>

#include "m_util.h"

#include "MyRenderFunction.h"
#include "GeometryFactory.h"

#include "open3d/Open3D.h"

using namespace open3d::visualization;
using namespace open3d::geometry;
using namespace myvisualization;
using std::cout;
using std::endl;
using std::vector;

#define up(l,r,i) for(int i=l;i<=r;i++)

inline double getDistance(Eigen::Vector3d o1, Eigen::Vector3d o2){
    auto del = o1 - o2;
    return std::sqrt(del.dot(del));
}

void AddFunction(MyVisualizer* vis, const std::function<double(const Eigen::Vector2d&)> func,
                 const AxisAlignedBoundingBox& boundingbox,  const double& resolution, 
                 const Eigen::Vector3d& color = {1,0,0},
                 const Eigen::Vector3d& trans = {0,0,0},
                 const Eigen::Vector3d& rotateAxis = {0,0,1}, const Eigen::Vector3d& center = {0,0,0},
                 const double& angle = 0
                 ){
        auto surface = mygeometry::CreateSurface(func, boundingbox, resolution);
        auto nrotate = rotateAxis.normalized();
        Eigen::AngleAxisd rotation(angle, nrotate);
        surface->Rotate(rotation.matrix(), center);
        surface->Translate(trans);
        surface->ComputeVertexNormals();
        surface->ComputeTriangleNormals();
        surface->PaintUniformColor(color);
        vis->AddGeometry(surface);

        auto cloud = std::make_shared<PointCloud>(surface->vertices_);
        vis->AddGeometry(cloud);
}

int main(int argc, char *argv[]) {

    std::mt19937 random_generator(42);


    MyVisualizer scene;

    if(!scene.CreateVisualizerWindow("Scene", 1600, 900)){
        std:: cout<<"[Init Error] Window Create Failed" << endl;
        return 1;
    }

    cout << "[Init] Window Created" << endl;

    glewExperimental = true;
    if (glewInit() != GLEW_OK) {
        utility::LogWarning("Failed to initialize GLEW.");
        return 1;
    }
    std :: cout << "[Init] InitOpenGL Succeed" << endl;

    auto coordinate_frame = TriangleMesh::CreateCoordinateFrame(5);
    scene.AddGeometry(coordinate_frame, true, false);

    cout << "[Init] AddGeometry Succeed " << endl;

    auto lineset = mygeometry::CreateCoordinateFrame({1,1,1},AxisAlignedBoundingBox({-100,-100,0},{100,100,0}));
    scene.AddGeometry(lineset,false,false);

    cout << "[Init] Add Frame Corrdinate Succeed " << endl;

    std::function<double(const Eigen::Vector2d&)> f0 = [](const Eigen::Vector2d& coor){
        auto x = coor(0), y = coor(1);
        double ret = 8/(1+exp(5-0.5*sqrt(x*x+y*y))-0.01*((x-y)*(x-y)));
        if(y <= 5*log(x+1) && y >= exp(x/5)-1) 
            return ret>0?ret:-11.0;
        return -11.0;
    };

    std::function<double(const Eigen::Vector2d&)> f0_5 = [](const Eigen::Vector2d& coor){
        auto x = coor(0)*1.25, y = coor(1)*1.25;
        double ret = 8/(1+exp(4-0.5*sqrt(x*x+y*y))-0.01*((x-y)*(x-y)));
        if(y <= 4.5*log(x+1) && y >= exp(x/4.5)-1) 
            return ret>0?ret:-11.0;
        return -11.0;
    };

    std::function<double(const Eigen::Vector2d&)> f1 = [](const Eigen::Vector2d& coor){
        auto x = coor(0)*1.5, y = coor(1)*1.5;
        double ret = 8/(1+exp(3-0.5*sqrt(x*x+y*y))-0.01*((x-y)*(x-y)));
        if(y <= 4.5*log(x+1) && y >= exp(x/4.5)-1) 
            return ret>0?ret:-11.0;
        return -11.0;
    };


    std::function<double(const Eigen::Vector2d&)> f2 = [](const Eigen::Vector2d& coor){
        auto x = coor(0)*3, y = coor(1)*3;
        double ret = 8/(1+exp(5-0.5*sqrt(x*x+y*y))+(-(30-x)*(30-x)/50)*sqrt((x-y)*(x-y)));
        if(y <= 4.5*log(x+1) && y >= exp(x/4.5)-1) 
            return ret>0?ret:-11.0;
        return -11.0;
    };
    //abs(sin(x*x+2*x*y)) - sin(x-2*y)
    //2/(1+pow(m_E, 5-2*sqrt(x*x+y*y)));
    //if(y < exp(x)-1 || y > log(x+1)) return -11.0;
    //    return 2/(1+exp(5-2*sqrt(x*x+y*y)));
    //best 彼岸花
    //return 8/(1+exp(5-0.5*sqrt(x*x+y*y))+(-(30-x)*(30-x)/50)*sqrt((x-y)*(x-y)));

    //
    //if(y <= 4.5*log(x+1) && y >= exp(x/4.5)-1) 
    //        return 8/(1+exp(5-0.5*sqrt(x*x+y*y))+(-(30-x)*(30-x)/50)*sqrt((x-y)*(x-y)));

    for(double i = 0; i < m_PI*2; i += m_PI/4)
        AddFunction(&scene, f1, {{0,0,-10},{20,20,10}}, 0.1,{1,0,0},{0,0,-0.5},{0,0,1},{0,0,0},i);

    for(double i = m_PI/3; i < m_PI*2+m_PI/3; i += m_PI/4)
        AddFunction(&scene, f0_5, {{0,0,-10},{20,20,10}}, 0.1,{1,0,0},{0,0,0},{0,0,1},{0,0,0},i);

    for(double i = 7*m_PI/15; i < m_PI*2+7*m_PI/15; i += m_PI/4)
        AddFunction(&scene, f0, {{0,0,-10},{20,20,10}}, 0.1,{1,0,0},{0,0,0},{0,0,1},{0,0,0},i);

    for(double i = 0; i < m_PI*2; i += m_PI/4)
        AddFunction(&scene, f2, {{0,0,-10},{20,20,10}}, 0.05,{1,0,0},{0,0,0},{0,0,1},{0,0,0},i);

    
    
    //AddFunction(&scene, f, {{-20,-20,-10},{20,20,10}}, 0.1,{1,0,0},{0,0,0},{0,0,1},{0,0,0},0);

    scene.Run();
    scene.DestroyVisualizerWindow();

    return 0;
}
