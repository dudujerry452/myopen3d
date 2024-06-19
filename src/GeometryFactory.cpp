#define __FUNCSIG__ __func__ //替代Visual Studio的获取函数签名的宏

#include "GeometryFactory.h"

using namespace open3d::visualization;
using namespace open3d::geometry;
using namespace myvisualization;

namespace mygeometry{

std::shared_ptr<LineSet> CreateCoordinateFrame(const Eigen::Vector3d step, 
                 const AxisAlignedBoundingBox& boundingbox){

    std::vector<Eigen::Vector3d> points;
    std::vector<Eigen::Vector2i> lines;

    for(int dim = 0; dim < 3; dim ++){
        for(int dim2 = 0; dim2 < 3; dim2 ++){
            if(dim == dim2) continue;
            for(double x = boundingbox.GetMinBound()(dim); x <= boundingbox.GetMaxBound()(dim); x += step(dim)){
                auto db = boundingbox.GetMinBound()(dim2),ub = boundingbox.GetMaxBound()(dim2);
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

std::shared_ptr<TriangleMesh> CreateSurface(const std::function<double(const Eigen::Vector2d&)> func,
                 const AxisAlignedBoundingBox& boundingbox,  const double& resolution){
    
    auto mesh = std::make_shared<TriangleMesh>();

    auto& points_ = mesh->vertices_;
    auto& triangles_ = mesh->triangles_;

    std::vector<int> prefix;
    std::vector<bool> invis;

    auto minbd = boundingbox.GetMinBound();
    auto maxbd = boundingbox.GetMaxBound();

    int base1 = 0,base2 = 0;
    for(double x = minbd(0); x <= maxbd(0); x += resolution) base1++;
    for(double y = minbd(1); y <= maxbd(1); y += resolution) base2++;
    prefix.resize(base1*base2+1);
    prefix[0] = 0;
    invis.resize(base1*base2);
    invis.assign(base1*base2, false);

    if(base1 <= 1 || base2 <= 1)  return mesh;

    int counter = 1;
    for(double x = minbd(0); x <= maxbd(0); x += resolution){
        for(double y = minbd(1); y <= maxbd(1); y += resolution){
            double ret = func({x,y});
            if(ret > maxbd(2) || ret < minbd(2)){ 
                prefix[counter] = 1;
                invis[counter-1] = true;
            }
            else
                points_.push_back({x,y,ret});
            prefix[counter] += prefix[counter-1];
            counter++;
        }
    }

    //id = base2*x + y

    for(int i = 0; i < base1-1; i ++){
        for(int j = 0; j < base2-1; j++){

            int a = base2*i+j, b = base2*i+j+1, c = base2*(i+1)+j, d = base2*(i+1)+j+1;
            int a1 = a-prefix[a+1], b1 = b-prefix[b+1], c1 = c-prefix[c+1], d1 = d-prefix[d+1];

            if(invis[a] || invis[b] || invis[c]){
                if(!invis[b] && !invis[c] && !invis[d]){
                    triangles_.push_back({d1, c1, b1});
                    triangles_.push_back({d1, b1, c1});
                }
            }
            else{
                triangles_.push_back({a1, b1, c1});
                triangles_.push_back({a1, c1, b1});
                if(!invis[d]) {
                    triangles_.push_back({d1, c1, b1});
                    triangles_.push_back({d1, b1, c1});
                }
            }
        }
    }

    return mesh;



}


};