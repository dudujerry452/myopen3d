
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
                 const Eigen::Vector3d& trans = {0,0,0},  const Eigen::Vector3d& color = {1,0,0}){
        auto surface = mygeometry::CreateSurface(func, boundingbox, resolution);
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

    //scene.RegisterAnimationCallback(update);
    //cout << "[Init] Animation Function Register Succeed " << endl;

    std::function<double(const Eigen::Vector2d&)> f = [](const Eigen::Vector2d& coor){
        auto x = coor(0), y = coor(1);
        if(x*x + y*y >= 7) return -1.0;
        return pow(m_E, sin(x) + cos(y));
    };
    std::function<double(const Eigen::Vector2d&)> f2 = [](const Eigen::Vector2d& coor){
        auto x = coor(0), y = coor(1);
        if(x*x + y*y >= 7) return -1.0;
        return pow(m_E, x);
    };

    AddFunction(&scene, f, {{-10,-10,0},{10,10,10}}, 0.1);

    AddFunction(&scene, f2, {{-10,-10,0},{10,10,10}}, 0.1, {10,0,0});

    scene.Run();
    scene.DestroyVisualizerWindow();

    return 0;
}
