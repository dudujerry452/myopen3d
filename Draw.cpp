
#define __FUNCSIG__ __func__ //替代Visual Studio的获取函数签名的宏

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>

#include "MyRenderFunction.h"

#include "open3d/Open3D.h"

using namespace open3d::visualization;
using namespace open3d::geometry;
using namespace myvisualization;
using std::cout;
using std::endl;
using std::vector;

#define up(l,r,i) for(int i=l;i<=r;i++)

#define chnow() std::chrono::steady_clock::now()
#define mildiff(x) std::chrono::duration_cast<std::chrono::milliseconds>(x).count()

const int64_t animation_delay = 10;
const int resolution = 50;
const double minError = 0.05;
const double init_radius = 0.01;

bool getData(const char * path, std::string * data){
    std::ifstream file(path);
    if(!file.is_open()) return false;
    file.seekg(0, std::ios::end);
    data->resize(file.tellg());
    file.seekg(0, std::ios::beg);
    file.read(&(*data)[0],data->size());
    file.close();
    return true;
}

float getDistanceFromEye(Eigen::Vector3d eye_position, Eigen::Vector3d target_position, double radius){
    auto delta = target_position - eye_position;
    double dis = std::sqrt(delta.dot(delta));
    if(dis < radius) dis = dis + radius;
    else dis = dis - radius;
    return dis;
}

Eigen::Vector3d eye_pos;

bool sphereCmp(const Eigen::Vector4d& s1, const Eigen::Vector4d& s2){
    return getDistanceFromEye(eye_pos, s1.head<3>(), s1(3)) > getDistanceFromEye(eye_pos, s2.head<3>(), s2(3));
}

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

bool getInteractivePart(const std::vector<Eigen::Vector3d>& geo, 
                        Eigen::Vector3d o1, Eigen::Vector3d o2, double r1, double r2, 
                         std::vector<Eigen::Vector3d>& result){
    bool flag = false;
    std::vector<Eigen::Vector3d> points;
    std::vector<Eigen::Vector2i> lines;
    for(const auto vertex : geo) {
        auto err = std::abs(getDistance(o1,vertex)-r1) + std::abs(getDistance(o2,vertex)-r2);
        if(err <= minError) {points.push_back(vertex);flag = true;}
    }

    result = points;
    /*if(points.empty()) return flag;
    auto first = points[0];
    size_t len = points.size();
    up(0,len-1,i){
        

    }*/
    return flag;
}

int main(int argc, char *argv[]) {

    std::mt19937 random_generator(42);

    std::string data;
    std::vector< Eigen::Vector4d> spheres_coor; 

    if(!getData("./data.dat", &data)){
        cout <<"File invalid !" << endl;
        return 1;
    }
    std::stringstream ss;
    ss << data;
    int n;
    double step;
    double k[5];//经纬高程声速单位km
    ss>>n>>step;//输入时间增加步长，以秒为单位
    double biggest = 0;
    up(0,3,i) ss>>k[i];
    k[4] = 1/k[3];
    Eigen::Vector4d bias;
    ss>>bias(0)>>bias(1)>>bias(2); // 输入坐标的偏移量
    bias(3) = 0;
    up(0,n-1,i){
        Eigen::Vector4d v4d;
        up(0,3,j) {
            ss>>v4d(j);
            v4d(j) += bias(j);
            v4d[j] *= k[j];
        }
        spheres_coor.push_back(v4d);
        biggest = (v4d(3)>biggest)?v4d(3):biggest;
    }
    for (auto &v : spheres_coor) v(3) -= biggest; 
    //cube->ComputeVertexNormals();
    //cube->PaintUniformColor({1.0, 0.0, 0.0});

    std::vector< Eigen::Vector4d> points_coor;
    int n2;
    ss>>n2;
    points_coor.resize(n2+1);
    up(0, n2-1, i){
        up(0,3,j) {
            ss>>points_coor[i](j);
            points_coor[i](j) += bias[j];
            points_coor[i](j) *= k[j];
        }
        points_coor[i](3) -= biggest;
    }


    vector<std::shared_ptr<TriangleMesh> > spheres;
    up(0,n-1,i){
        spheres.emplace_back(TriangleMesh::CreateSphere(init_radius,resolution));
        for (auto& v : spheres[i]->vertices_){
            v += spheres_coor[i].head<3>();
        }
    }

    MyVisualizer scene;

    std::vector<std::shared_ptr<PointCloud> > inter_points;

    std::vector<std::shared_ptr<PointCloud> > draw_points_group;
    std::vector<std::shared_ptr<PointCloud> > temp_draw_points_group;
    // Construct points to draw
    std::uniform_real_distribution<double> dis(0.0, 2.0*M_PI);
    for(auto& v : points_coor){
        std::vector<Eigen::Vector3d> v_group;
        v_group.emplace_back(v.head<3>());
        up(0,99,i){
            double th = dis(random_generator);
            double ph = dis(random_generator);
            double co = cos(ph);
            v_group.emplace_back(v(0)+co*sin(th),v(1)+co*cos(th),v(2)+sin(ph));
        }
        draw_points_group.emplace_back(std::make_shared<PointCloud>(v_group));
        temp_draw_points_group.emplace_back(std::make_shared<PointCloud>(v_group));
    }

    std::function<bool(MyVisualizer *)> update = [&points_coor, &temp_draw_points_group, &draw_points_group, &biggest, &inter_points, &spheres,step,k,&spheres_coor,n](MyVisualizer* vis){    
        if(!vis) return false;

        double current_time = biggest*k[4] - vis->GetAnimationFrame()*step;
        if(vis->PushKeyDown()){
            printf("[Press \"Down\" Key] Current Time : %lf s (x%lf)\n", current_time, step);
        }
        if(!vis->IsNextFrame()) return false;

        auto view = vis->GetViewControl();
        eye_pos = view.GetEye().cast<double>();
        std::sort(spheres_coor.begin(), spheres_coor.end(), sphereCmp);

        up(0,n-1,i){

            double radius = spheres_coor[i](3) + vis->GetAnimationFrame()*step*k[3];
            if(radius <= 0) radius = init_radius;
            
            auto sph = TriangleMesh::CreateSphere(radius,resolution);
            spheres[i]->vertices_ = sph->vertices_;
            spheres[i]->triangles_ = sph->triangles_;
            for(auto& v : spheres[i]->vertices_){
                v += spheres_coor[i].head<3>();
            }
            spheres[i]->ComputeVertexNormals();
            spheres[i]->ComputeTriangleNormals();
            spheres[i]->PaintUniformColor({0,1,0});
        }

        for(const auto & remov : inter_points){ // 移除上次渲染的交线
            vis->RemoveGeometry(remov);
        }
        inter_points.clear();

        up(0,n-1,i){
            double radius1 = spheres_coor[i](3) + vis->GetAnimationFrame()*step*k[3];

            up(i+1,n-1,j){
                double radius2 = spheres_coor[j](3) + vis->GetAnimationFrame()*step*k[3];

                if(getDistance(spheres_coor[i].head<3>(),spheres_coor[j].head<3>()) > radius1 + radius2) continue;
                std::vector<Eigen::Vector3d> result;
                if(getInteractivePart(spheres[i]->vertices_, 
                            spheres_coor[i].head<3>(), spheres_coor[j].head<3>(),
                            radius1, radius2, 
                             result)){
                    inter_points.emplace_back(std::make_shared<PointCloud>(result));
                    inter_points.back()->PaintUniformColor({0,0,0});
                    vis->AddGeometryA(inter_points.back(), false, '0', false);
                    
                }
            }
        }

        int index = 0;
        for(auto& v : draw_points_group){
            v->points_ = temp_draw_points_group[index]->points_;
            v->Scale(10*sin((vis->GetAnimationFrame())/5.0),v->points_[0]);
            if(abs(current_time - points_coor[index](3)) < step){
                vis->AnimationPause();
                printf("[Preset Solution Time Reached] Target Time : %lf s\n", points_coor[index](3));
                points_coor[index](3) = 9999.0;
            }
            index++;
        }
        
        vis->UpdateGeometry();
        return true;
    };


    auto coordinate_frame = TriangleMesh::CreateCoordinateFrame();

    if(!scene.CreateVisualizerWindow("Scene", 1600, 900)){
        std:: cout<<"[Init Error] Window Create Failed" << endl;
        return 1;
    }

    //auto renderer_ptr = std::make_shared<glsl::TriangleMeshRenderer>();

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

    up(0,n-1,i){
        scene.AddGeometry(spheres[i]);
    }
    cout << "[Init] Add Custom Balls Succeed " << endl;

    up(0,n2-1,i){
        scene.AddGeometryA(draw_points_group[i], false, '0', true);
    }

    scene.RegisterAnimationCallback(update);
    cout << "[Init] Animation Function Register Succeed " << endl;
    scene.Run();
    scene.DestroyVisualizerWindow();

    return 0;
}
