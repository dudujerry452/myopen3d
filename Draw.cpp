
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

#include "open3d/Open3D.h"

using namespace open3d::visualization;
using namespace open3d::geometry;
using namespace myvisualization;
using std::cout;
using std::endl;
using std::vector;

#define up(l,r,i) for(int i=l;i<=r;i++)

std::shared_ptr<LineSet> getFrameCorrdinate(Eigen::Vector3d dn_bound, Eigen::Vector3d up_bound, Eigen::Vector3d step){
    std::vector<Eigen::Vector3d> points;
    std::vector<Eigen::Vector2i> lines;

    for(int dim = 0; dim < 3; dim ++){
        for(int dim2 = 0; dim2 < 3; dim2 ++){
            if(dim == dim2) continue;
            for(double x = dn_bound(dim); x <= up_bound(dim); x += step(dim)){
                auto db = dn_bound(dim2),ub = up_bound(dim2);
                Eigen::Vector3d p({0,0,0});
                p(dim) = x; p(dim2) = db;
                points.push_back(p);
                p(dim2) = ub;
                points.push_back(p);
                lines.emplace_back(points.size()-2, points.size()-1);
            }
        }
    }
    auto lineset = std::make_shared<LineSet>(points,lines);
    return lineset;
}


inline double getDistance(Eigen::Vector3d o1, Eigen::Vector3d o2){
    auto del = o1 - o2;
    return std::sqrt(del.dot(del));
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

    auto lineset = getFrameCorrdinate({0,0,0},{100,100,100},{1,1,1});
    scene.AddGeometry(lineset,false,false);

    cout << "[Init] Add Frame Corrdinate Succeed " << endl;

    //scene.RegisterAnimationCallback(update);
    //cout << "[Init] Animation Function Register Succeed " << endl;

    auto sphere = TriangleMesh::CreateSphere(10,20);
    sphere->Translate({10,10,10});
    sphere->ComputeVertexNormals();
    sphere->ComputeTriangleNormals();
    sphere->PaintUniformColor({0,1,0});
    scene.AddGeometry(sphere);

    scene.Run();
    scene.DestroyVisualizerWindow();

    return 0;
}
