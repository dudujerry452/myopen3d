#define __FUNCSIG__ __func__ //替代Visual Studio的获取函数签名的宏

#ifndef _GEOMETRY_FACTORY_H_
#define _GEOMETRY_FACTORY_H_

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

namespace mygeometry{

    std::shared_ptr<LineSet> CreateCoordinateFrame(const Eigen::Vector3d step, 
                 const AxisAlignedBoundingBox& boundingbox);

    std::shared_ptr<TriangleMesh> CreateSurface(const std::function<double(const Eigen::Vector2d&)> func,
                 const AxisAlignedBoundingBox& boundingbox,  const double& resolution);


};

#endif