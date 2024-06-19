
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

std::shared_ptr<LineSet> getFrameCorrdinate(double boundx, double boundy, double stepx, double stepy){
    std::vector<Eigen::Vector3d> points;
    std::vector<Eigen::Vector2i> lines;
    int timesx = boundx/stepx;
    int timesy = boundy/stepy;
    up(0,timesx,i){
        points.push_back({stepx*i,0,0});
        points.push_back({stepx*i,boundy,0});
        lines.push_back({i*2,i*2+1});
    }
    up(0,timesy,i){
        points.push_back({0,stepy*i,0});
        points.push_back({boundx,stepy*i,0});
        lines.push_back({(timesx+i+1)*2,(timesx+i+1)*2+1});
    }
    auto lineset = std::make_shared<LineSet>(points, lines);
    return lineset;
}


inline double getDistance(Eigen::Vector3d o1, Eigen::Vector3d o2){
    auto del = o1 - o2;
    return std::sqrt(del.dot(del));
}

int main(int argc, char *argv[]) {

    std::mt19937 random_generator(42);


    MyVisualizer scene;

    auto coordinate_frame = TriangleMesh::CreateCoordinateFrame();

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

    /*if (glewInit() != GLEW_OK) {
        std:: cout << "GLEW Intilization Failed" << std::endl;
        return false;
    }
    cout << "Window Created " << endl;*/

    scene.AddGeometry(coordinate_frame, true, false);

    cout << "[Init] AddGeometry Succeed " << endl;

    auto lineset = getFrameCorrdinate(1000,1000,100,100);
    scene.AddGeometry(lineset,false,false);

    cout << "[Init] Add Frame Corrdinate Succeed " << endl;

    //scene.RegisterAnimationCallback(update);
    cout << "[Init] Animation Function Register Succeed " << endl;
    scene.Run();
    scene.DestroyVisualizerWindow();

    return 0;
}
