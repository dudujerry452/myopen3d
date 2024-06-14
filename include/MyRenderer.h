#ifndef _MYRENDERER_H_
#define _MYRENDERER_H_

#include "open3d/Open3D.h"
#include "MyShaders.h"
#include "MyGeometryRenderer.h"

#define __FUNCSIG__ __func__ //替代Visual Studio的获取函数签名的宏

using namespace open3d::visualization::glsl;
using namespace open3d::visualization;
using namespace open3d;

using namespace myglsl;



class MyTriangleMeshRenderer : public MyGeometryRenderer {
public:
    ~MyTriangleMeshRenderer() override {}

public:
    bool Render(const RenderOption &option, const MyViewControl &view) override;
    bool AddGeometry(
            std::shared_ptr<const open3d::geometry::Geometry> geometry_ptr) override;
    bool UpdateGeometry() override;

protected:
    
    MyPhongShaderForTriangleMesh phong_mesh_shader_; // Target to Modify (Add "My")
};

class MyLineSetRenderer : public MyGeometryRenderer {
public:
    ~MyLineSetRenderer() override {}

public:
    bool Render(const RenderOption &option, const MyViewControl &view) override;
    bool AddGeometry(
            std::shared_ptr<const geometry::Geometry> geometry_ptr) override;
    bool UpdateGeometry() override;

protected:
    MySimpleShaderForLineSet simple_lineset_shader_;
};

class MyPointCloudRenderer : public MyGeometryRenderer {
public:
    ~MyPointCloudRenderer() override {}

public:
    bool Render(const RenderOption &option, const MyViewControl &view) override;
    bool AddGeometry(
            std::shared_ptr<const geometry::Geometry> geometry_ptr) override;
    bool UpdateGeometry() override;

protected:
    MySimpleShaderForPointCloud simple_point_shader_;
    //PhongShaderForPointCloud phong_point_shader_;
    //NormalShaderForPointCloud normal_point_shader_;
    //SimpleBlackShaderForPointCloudNormal simpleblack_normal_shader_;
};

#endif 