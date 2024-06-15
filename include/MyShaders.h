#ifndef _MYSHADERS_H_
#define _MYSHADERS_H_

#include "open3d/Open3D.h"
#include "MySpecificShaders.h"
#include "MyViewControl.h"

#define __FUNCSIG__ __func__ //替代Visual Studio的获取函数签名的宏

#define chnow() std::chrono::steady_clock::now()
#define mildiff(x) std::chrono::duration_cast<std::chrono::milliseconds>(x).count()

using namespace open3d::visualization::glsl;
using namespace open3d::visualization;
using namespace open3d;

using namespace myvisualization;

// Work flow : Compile -> PrepareBinding -> BindGeometry -> PrepareRendering -> RenderGeometry

namespace myglsl{

    class ShaderWrapper {
public:
    virtual ~ShaderWrapper() {}
    ShaderWrapper(const ShaderWrapper &) = delete;
    ShaderWrapper &operator=(const ShaderWrapper &) = delete;

protected:
    ShaderWrapper(const std::string &name) : shader_name_(name) {}

public:
    /// Function to render geometry under condition of mode and view
    /// The geometry is updated in a passive manner (bind only when needed).
    /// Thus this function compiles shaders if not yet, binds geometry if not
    /// yet, then do the rendering.
    bool Render(const geometry::Geometry &geometry,
                const RenderOption &option,
                const MyViewControl &view);

    /// Function to invalidate the geometry (set the dirty flag and release
    /// geometry resource)
    void InvalidateGeometry();

    const std::string &GetShaderName() const { return shader_name_; }

    void PrintShaderWarning(const std::string &message) const;

protected:
    /// Function to compile shader
    /// In a derived class, this must be declared as final, and called from
    /// the constructor.
    virtual bool Compile() = 0;

    /// Function to release resource
    /// In a derived class, this must be declared as final, and called from
    /// the destructor.
    virtual void Release() = 0;

    virtual bool BindGeometry(const geometry::Geometry &geometry,
                              const RenderOption &option,
                              const MyViewControl &view) = 0;
    virtual bool RenderGeometry(const geometry::Geometry &geometry,
                                const RenderOption &option,
                                const MyViewControl &view) = 0;
    virtual void UnbindGeometry() = 0;

protected:
    bool ValidateShader(GLuint shader_index);
    bool ValidateProgram(GLuint program_index);
    bool CompileShaders(const char *const vertex_shader_code,
                        const char *const geometry_shader_code,
                        const char *const fragment_shader_code);
    void ReleaseProgram();

protected:
    GLuint vertex_shader_ = 0;
    GLuint geometry_shader_ = 0;
    GLuint fragment_shader_ = 0;
    GLuint program_ = 0;
    GLenum draw_arrays_mode_ = GL_POINTS;
    GLsizei draw_arrays_size_ = 0;
    bool compiled_ = false;
    bool bound_ = false;

    void SetShaderName(const std::string &shader_name) {
        shader_name_ = shader_name;
    }

private:
    std::string shader_name_ = "ShaderWrapper";
};

class MySimpleShader : public ShaderWrapper {
public:
    ~MySimpleShader() override { Release(); }

protected:
    MySimpleShader(const std::string &name) : ShaderWrapper(name) { Compile(); }

protected:
    bool Compile() final;
    void Release() final;
    bool BindGeometry(const geometry::Geometry &geometry,
                      const RenderOption &option,
                      const MyViewControl &view) final;
    bool RenderGeometry(const geometry::Geometry &geometry,
                        const RenderOption &option,
                        const MyViewControl &view) final;
    void UnbindGeometry() final;

protected:
    virtual bool PrepareRendering(const geometry::Geometry &geometry,
                                  const RenderOption &option,
                                  const MyViewControl &view) = 0;
    virtual bool PrepareBinding(const geometry::Geometry &geometry,
                                const RenderOption &option,
                                const MyViewControl &view,
                                std::vector<Eigen::Vector3f> &points,
                                std::vector<Eigen::Vector3f> &colors) = 0;

protected:
    GLuint vertex_position_;
    GLuint vertex_position_buffer_;
    GLuint vertex_color_;
    GLuint vertex_color_buffer_;
    GLuint MVP_;
};

}

class MySimpleShaderForPointCloud : public myglsl::MySimpleShader {
public:
    MySimpleShaderForPointCloud() : myglsl::MySimpleShader("SimpleShaderForPointCloud") {}

protected:
    bool PrepareRendering(const geometry::Geometry &geometry,
                          const RenderOption &option,
                          const MyViewControl &view) final;
    bool PrepareBinding(const geometry::Geometry &geometry,
                        const RenderOption &option,
                        const MyViewControl &view,
                        std::vector<Eigen::Vector3f> &points,
                        std::vector<Eigen::Vector3f> &colors) final;
};

class MySimpleShaderForLineSet : public myglsl::MySimpleShader {
public:
    MySimpleShaderForLineSet() : myglsl::MySimpleShader("SimpleShaderForLineSet") {}

protected:
    bool PrepareRendering(const geometry::Geometry &geometry,
                          const RenderOption &option,
                          const MyViewControl &view) final;
    bool PrepareBinding(const geometry::Geometry &geometry,
                        const RenderOption &option,
                        const MyViewControl &view,
                        std::vector<Eigen::Vector3f> &points,
                        std::vector<Eigen::Vector3f> &colors) final;
};

class MyPhongShader : public myglsl::ShaderWrapper {
public:
    ~MyPhongShader() override { Release(); }

protected:
    MyPhongShader(const std::string &name) : myglsl::ShaderWrapper(name) {std::cout<<"MyPhongShader" << std::endl;Compile(); }

protected:
    bool Compile() final; // 1 Create and Link program, Compile shaders and Get glAttribution from shaders
    void Release() final;
    bool BindGeometry(const open3d::geometry::Geometry &geometry, // 3 Generate glBuffers and Input the data into GPU
                      const RenderOption &option,
                      const MyViewControl &view) final;
    bool RenderGeometry(const open3d::geometry::Geometry &geometry, // 5 Use the data to call glDrawArray to Render picture to screen
                        const RenderOption &option,
                        const MyViewControl &view) final;
    void UnbindGeometry() final;

protected:
    virtual bool PrepareRendering(const open3d::geometry::Geometry &geometry, // 4 Do the glSettings
                                  const RenderOption &option,
                                  const MyViewControl &view) = 0;
    virtual bool PrepareBinding(const open3d::geometry::Geometry &geometry, // 2 Copy vertex data and Calculate matrix
                                const RenderOption &option,
                                const MyViewControl &view,
                                std::vector<Eigen::Vector3f> &points,
                                std::vector<Eigen::Vector3f> &normals,
                                std::vector<Eigen::Vector3f> &colors) = 0;

protected:
    void SetLighting(const MyViewControl &view, const RenderOption &option);

protected:
    GLuint vertex_position_;
    GLuint vertex_position_buffer_;
    GLuint vertex_color_;
    GLuint vertex_color_buffer_;
    GLuint vertex_normal_;
    GLuint vertex_normal_buffer_;
    GLuint MVP_;
    GLuint V_;
    GLuint M_;
    GLuint light_position_world_;
    GLuint light_color_;
    GLuint light_diffuse_power_;
    GLuint light_specular_power_;
    GLuint light_specular_shininess_;
    GLuint light_ambient_;

    // At most support 4 lights
    gl_util::GLMatrix4f light_position_world_data_;
    gl_util::GLMatrix4f light_color_data_;
    gl_util::GLVector4f light_diffuse_power_data_;
    gl_util::GLVector4f light_specular_power_data_;
    gl_util::GLVector4f light_specular_shininess_data_;
    gl_util::GLVector4f light_ambient_data_;
};

class MyPhongShaderForTriangleMesh : public MyPhongShader {
public:
    MyPhongShaderForTriangleMesh() : MyPhongShader("MyPhongShaderForTriangleMesh") {}

protected:
    bool PrepareRendering(const geometry::Geometry &geometry,
                          const RenderOption &option,
                          const MyViewControl &view) final;
    bool PrepareBinding(const geometry::Geometry &geometry,
                        const RenderOption &option,
                        const MyViewControl &view,
                        std::vector<Eigen::Vector3f> &points,
                        std::vector<Eigen::Vector3f> &normals,
                        std::vector<Eigen::Vector3f> &colors) final;
};


#endif