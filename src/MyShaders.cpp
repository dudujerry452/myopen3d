#include "MyShaders.h"

bool MySimpleShaderForLineSet::PrepareRendering(
        const geometry::Geometry &geometry,
        const RenderOption &option,
        const MyViewControl &view) {
    if (geometry.GetGeometryType() !=
        geometry::Geometry::GeometryType::LineSet) {
        PrintShaderWarning("Rendering type is not geometry::LineSet.");
        return false;
    }
    glLineWidth(GLfloat(option.line_width_));
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GLenum(option.GetGLDepthFunc()));
    return true;
}

bool MySimpleShaderForLineSet::PrepareBinding(
        const geometry::Geometry &geometry,
        const RenderOption &option,
        const MyViewControl &view,
        std::vector<Eigen::Vector3f> &points,
        std::vector<Eigen::Vector3f> &colors) {
    if (geometry.GetGeometryType() !=
        geometry::Geometry::GeometryType::LineSet) {
        PrintShaderWarning("Rendering type is not geometry::LineSet.");
        return false;
    }
    const geometry::LineSet &lineset = (const geometry::LineSet &)geometry;
    if (!lineset.HasLines()) {
        PrintShaderWarning("Binding failed with empty geometry::LineSet.");
        return false;
    }
    points.resize(lineset.lines_.size() * 2);
    colors.resize(lineset.lines_.size() * 2);
    for (size_t i = 0; i < lineset.lines_.size(); i++) {
        const auto point_pair = lineset.GetLineCoordinate(i);
        points[i * 2] = point_pair.first.cast<float>();
        points[i * 2 + 1] = point_pair.second.cast<float>();
        Eigen::Vector3d color;
        if (lineset.HasColors()) {
            color = lineset.colors_[i];
        } else {
            color = Eigen::Vector3d::Zero();
        }
        colors[i * 2] = colors[i * 2 + 1] = color.cast<float>();
    }
    draw_arrays_mode_ = GL_LINES;
    draw_arrays_size_ = GLsizei(points.size());
    return true;
}

bool MyPhongShader::Compile() {
    std:: cout << "Init Compile()" << std::endl;
    if (!CompileShaders(MyShaders::PhongVertexShader, NULL, MyShaders::PhongFragmentShader)) {
        std:: cout << "Compile Failed." << std::endl;
        PrintShaderWarning("Compiling shaders failed.");
        return false;
    }
    std:: cout << "Compile finished." << std::endl;
    std:: cout << "_program ID = " << program_ << std::endl;
    std:: cout << "vertex_shader ID = " << vertex_shader_ << std::endl;

    if(glGetShaderiv == NULL) { 
        std:: cout << "glew Initilized Failed" << std::endl;
    }

    GLint compiled = GL_FALSE;
    glGetShaderiv(vertex_shader_, GL_COMPILE_STATUS, &compiled);

    if (compiled == GL_FALSE) {
        // 着色器编译失败，获取编译日志
        GLint logLength = 0;
        glGetShaderiv(vertex_shader_, GL_INFO_LOG_LENGTH, &logLength);
        GLchar* log = (GLchar*)malloc(logLength);
        glGetShaderInfoLog(vertex_shader_, logLength, NULL, log);
        printf("Shader compilation log: %s\n", log);
        free(log);
    }
    else{
        std:: cout << "vertex shader compile succeed " << vertex_shader_ << std::endl;
    }



    vertex_position_ = glGetAttribLocation(program_, "vertex_position");
    vertex_normal_ = glGetAttribLocation(program_, "vertex_normal");
    vertex_color_ = glGetAttribLocation(program_, "vertex_color");
    MVP_ = glGetUniformLocation(program_, "MVP");
    V_ = glGetUniformLocation(program_, "V");
    M_ = glGetUniformLocation(program_, "M");
    light_position_world_ =
            glGetUniformLocation(program_, "light_position_world_4");
    light_color_ = glGetUniformLocation(program_, "light_color_4");
    light_diffuse_power_ =
            glGetUniformLocation(program_, "light_diffuse_power_4");
    light_specular_power_ =
            glGetUniformLocation(program_, "light_specular_power_4");
    light_specular_shininess_ =
            glGetUniformLocation(program_, "light_specular_shininess_4");
    light_ambient_ = glGetUniformLocation(program_, "light_ambient");
    std:: cout << "After Compile" << std::endl;
    return true;
}

void MyPhongShader::Release() {
    UnbindGeometry();
    ReleaseProgram();
}

bool MyPhongShader::BindGeometry(const geometry::Geometry &geometry,
                               const RenderOption &option,
                               const MyViewControl &view) {
    // If there is already geometry, we first unbind it.
    // We use GL_STATIC_DRAW. When geometry changes, we clear buffers and
    // rebind the geometry. Note that this approach is slow. If the geometry is
    // changing per frame, consider implementing a new ShaderWrapper using
    // GL_STREAM_DRAW, and replace UnbindGeometry() with Buffer Object
    // Streaming mechanisms.
    UnbindGeometry();

    // Prepare data to be passed to GPU
    std::vector<Eigen::Vector3f> points;
    std::vector<Eigen::Vector3f> normals;
    std::vector<Eigen::Vector3f> colors;
    if (!PrepareBinding(geometry, option, view, points, normals, colors)) {
        PrintShaderWarning("Binding failed when preparing data.");
        return false;
    }

    // Create buffers and bind the geometry
    glGenBuffers(1, &vertex_position_buffer_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_position_buffer_);
    glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(Eigen::Vector3f),
                 points.data(), GL_STATIC_DRAW);
    glGenBuffers(1, &vertex_normal_buffer_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_normal_buffer_);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(Eigen::Vector3f),
                 normals.data(), GL_STATIC_DRAW);
    glGenBuffers(1, &vertex_color_buffer_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_color_buffer_);
    glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(Eigen::Vector3f),
                 colors.data(), GL_STATIC_DRAW);
    bound_ = true;
    return true;
}

bool MyPhongShader::RenderGeometry(const geometry::Geometry &geometry,
                                 const RenderOption &option,
                                 const MyViewControl &view) {
    if (!PrepareRendering(geometry, option, view)) {
        PrintShaderWarning("Rendering failed during preparation.");
        return false;
    }
    glUseProgram(program_);
    glUniformMatrix4fv(MVP_, 1, GL_FALSE, view.GetMVPMatrix().data());
    glUniformMatrix4fv(V_, 1, GL_FALSE, view.GetViewMatrix().data());
    glUniformMatrix4fv(M_, 1, GL_FALSE, view.GetModelMatrix().data());
    glUniformMatrix4fv(light_position_world_, 1, GL_FALSE,
                       light_position_world_data_.data());
    glUniformMatrix4fv(light_color_, 1, GL_FALSE, light_color_data_.data());
    glUniform4fv(light_diffuse_power_, 1, light_diffuse_power_data_.data());
    glUniform4fv(light_specular_power_, 1, light_specular_power_data_.data());
    glUniform4fv(light_specular_shininess_, 1,
                 light_specular_shininess_data_.data());
    glUniform4fv(light_ambient_, 1, light_ambient_data_.data());
    glEnableVertexAttribArray(vertex_position_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_position_buffer_);
    glVertexAttribPointer(vertex_position_, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(vertex_normal_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_normal_buffer_);
    glVertexAttribPointer(vertex_normal_, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(vertex_color_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_color_buffer_);
    glVertexAttribPointer(vertex_color_, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glDrawArrays(draw_arrays_mode_, 0, draw_arrays_size_);
    glDisableVertexAttribArray(vertex_position_);
    glDisableVertexAttribArray(vertex_normal_);
    glDisableVertexAttribArray(vertex_color_);
    return true;
}

void MyPhongShader::UnbindGeometry() {
    if (bound_) {
        glDeleteBuffers(1, &vertex_position_buffer_);
        glDeleteBuffers(1, &vertex_normal_buffer_);
        glDeleteBuffers(1, &vertex_color_buffer_);
        bound_ = false;
    }
}

void MyPhongShader::SetLighting(const MyViewControl &view,
                              const RenderOption &option) {
    const auto &box = view.GetBoundingBox();
    light_position_world_data_.block<3,1>(0, 0) = view.GetEye();
    //light_position_world_data_.setOnes();
    //light_position_world_data_ *= 100;
    light_color_data_.setOnes();
    for (int i = 0; i < 4; i++) {
        light_position_world_data_.block<3, 1>(0, i) =
                box.GetCenter().cast<GLfloat>() +
                (float)box.GetMaxExtent() *
                        ((float)option.light_position_relative_[i](0) *
                                 view.GetRight() +
                         (float)option.light_position_relative_[i](1) *
                                 view.GetUp() +
                         (float)option.light_position_relative_[i](2) *
                                 view.GetFront());
        light_color_data_.block<3, 1>(0, i) =
                option.light_color_[i].cast<GLfloat>();
    }
    if (option.light_on_) {
        light_diffuse_power_data_ =
                Eigen::Vector4d(option.light_diffuse_power_).cast<GLfloat>();
        light_specular_power_data_ =
                Eigen::Vector4d(option.light_specular_power_).cast<GLfloat>();
        light_specular_shininess_data_ =
                Eigen::Vector4d(option.light_specular_shininess_)
                        .cast<GLfloat>();
        light_ambient_data_.block<3, 1>(0, 0) =
                option.light_ambient_color_.cast<GLfloat>() * 100;
        light_ambient_data_(3) = 1.0f;
    } else {
        light_diffuse_power_data_ = gl_util::GLVector4f::Zero();
        light_specular_power_data_ = gl_util::GLVector4f::Zero();
        light_specular_shininess_data_ = gl_util::GLVector4f::Ones();
        light_ambient_data_ = gl_util::GLVector4f(1.0f, 1.0f, 1.0f, 1.0f);
    }
}

bool MyPhongShaderForTriangleMesh::PrepareRendering(
        const geometry::Geometry &geometry,
        const RenderOption &option,
        const MyViewControl &view) {
    if (geometry.GetGeometryType() !=
                geometry::Geometry::GeometryType::TriangleMesh &&
        geometry.GetGeometryType() !=
                geometry::Geometry::GeometryType::HalfEdgeTriangleMesh) {
        PrintShaderWarning("Rendering type is not geometry::TriangleMesh.");
        return false;
    }
    if (option.mesh_show_back_face_) {
        glDisable(GL_CULL_FACE);
    } else {
        glEnable(GL_CULL_FACE);
    }
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GLenum(option.GetGLDepthFunc()));
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    if (option.mesh_show_wireframe_) {
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0, 1.0);
    } else {
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
    SetLighting(view, option);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        //glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);

    return true;
}

bool MyPhongShaderForTriangleMesh::PrepareBinding(
        const geometry::Geometry &geometry,
        const RenderOption &option,
        const MyViewControl &view,
        std::vector<Eigen::Vector3f> &points,
        std::vector<Eigen::Vector3f> &normals,
        std::vector<Eigen::Vector3f> &colors) {
    if (geometry.GetGeometryType() !=
                geometry::Geometry::GeometryType::TriangleMesh &&
        geometry.GetGeometryType() !=
                geometry::Geometry::GeometryType::HalfEdgeTriangleMesh) {
        PrintShaderWarning("Rendering type is not geometry::TriangleMesh.");
        return false;
    }
    const geometry::TriangleMesh &mesh =
            (const geometry::TriangleMesh &)geometry;
    if (!mesh.HasTriangles()) {
        PrintShaderWarning("Binding failed with empty triangle mesh.");
        return false;
    }
    if (!mesh.HasTriangleNormals() || !mesh.HasVertexNormals()) {
        PrintShaderWarning("Binding failed because mesh has no normals.");
        PrintShaderWarning("Call ComputeVertexNormals() before binding.");
        return false;
    }
    const ColorMap &global_color_map = *GetGlobalColorMap();
    points.resize(mesh.triangles_.size() * 3);
    normals.resize(mesh.triangles_.size() * 3);
    colors.resize(mesh.triangles_.size() * 3);

    for (size_t i = 0; i < mesh.triangles_.size(); i++) {
        const auto &triangle = mesh.triangles_[i];
        for (size_t j = 0; j < 3; j++) {
            size_t idx = i * 3 + j;
            size_t vi = triangle(j);
            const auto &vertex = mesh.vertices_[vi];
            points[idx] = vertex.cast<float>();

            Eigen::Vector3d color;
            switch (option.mesh_color_option_) {
                case RenderOption::MeshColorOption::XCoordinate:
                    color = global_color_map.GetColor(
                            view.GetBoundingBox().GetXPercentage(vertex(0)));
                    break;
                case RenderOption::MeshColorOption::YCoordinate:
                    color = global_color_map.GetColor(
                            view.GetBoundingBox().GetYPercentage(vertex(1)));
                    break;
                case RenderOption::MeshColorOption::ZCoordinate:
                    color = global_color_map.GetColor(
                            view.GetBoundingBox().GetZPercentage(vertex(2)));
                    break;
                case RenderOption::MeshColorOption::Color:
                    if (mesh.HasVertexColors()) {
                        color = mesh.vertex_colors_[vi];
                        break;
                    }
                    // fallthrough
                case RenderOption::MeshColorOption::Default:
                default:
                    color = option.default_mesh_color_;
                    break;
            }
            colors[idx] = color.cast<float>();

            if (option.mesh_shade_option_ ==
                RenderOption::MeshShadeOption::FlatShade) {
                normals[idx] = mesh.triangle_normals_[i].cast<float>();
            } else {
                normals[idx] = mesh.vertex_normals_[vi].cast<float>();
            }
        }
    }
    draw_arrays_mode_ = GL_TRIANGLES;
    draw_arrays_size_ = GLsizei(points.size());
    return true;
}

bool MySimpleShaderForPointCloud::PrepareRendering(
        const geometry::Geometry &geometry,
        const RenderOption &option,
        const MyViewControl &view) {
    if (geometry.GetGeometryType() !=
        geometry::Geometry::GeometryType::PointCloud) {
        PrintShaderWarning("Rendering type is not geometry::PointCloud.");
        return false;
    }
    glPointSize(GLfloat(option.point_size_));
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GLenum(option.GetGLDepthFunc()));
    return true;
}

bool MySimpleShaderForPointCloud::PrepareBinding(
        const geometry::Geometry &geometry,
        const RenderOption &option,
        const MyViewControl &view,
        std::vector<Eigen::Vector3f> &points,
        std::vector<Eigen::Vector3f> &colors) {
    if (geometry.GetGeometryType() !=
        geometry::Geometry::GeometryType::PointCloud) {
        PrintShaderWarning("Rendering type is not geometry::PointCloud.");
        return false;
    }
    const geometry::PointCloud &pointcloud =
            (const geometry::PointCloud &)geometry;
    if (!pointcloud.HasPoints()) {
        PrintShaderWarning("Binding failed with empty pointcloud.");
        return false;
    }
    const ColorMap &global_color_map = *GetGlobalColorMap();
    points.resize(pointcloud.points_.size());
    colors.resize(pointcloud.points_.size());
    for (size_t i = 0; i < pointcloud.points_.size(); i++) {
        const auto &point = pointcloud.points_[i];
        points[i] = point.cast<float>();
        Eigen::Vector3d color;
        switch (option.point_color_option_) {
            case RenderOption::PointColorOption::XCoordinate:
                color = global_color_map.GetColor(
                        view.GetBoundingBox().GetXPercentage(point(0)));
                break;
            case RenderOption::PointColorOption::YCoordinate:
                color = global_color_map.GetColor(
                        view.GetBoundingBox().GetYPercentage(point(1)));
                break;
            case RenderOption::PointColorOption::ZCoordinate:
                color = global_color_map.GetColor(
                        view.GetBoundingBox().GetZPercentage(point(2)));
                break;
            case RenderOption::PointColorOption::Color:
            case RenderOption::PointColorOption::Default:
            default:
                if (pointcloud.HasColors()) {
                    color = pointcloud.colors_[i];
                } else {
                    color = global_color_map.GetColor(
                            view.GetBoundingBox().GetZPercentage(point(2)));
                }
                break;
        }
        colors[i] = color.cast<float>();
    }
    draw_arrays_mode_ = GL_POINTS;
    draw_arrays_size_ = GLsizei(points.size());
    return true;
}

namespace myglsl{

const static std::vector<Eigen::Vector3i> cuboid_vertex_offsets{
        Eigen::Vector3i(0, 0, 0), Eigen::Vector3i(1, 0, 0),
        Eigen::Vector3i(0, 1, 0), Eigen::Vector3i(1, 1, 0),
        Eigen::Vector3i(0, 0, 1), Eigen::Vector3i(1, 0, 1),
        Eigen::Vector3i(0, 1, 1), Eigen::Vector3i(1, 1, 1),
};

// Vertex indices of 12 triangles in a cuboid, for right-handed manifold mesh
const static std::vector<Eigen::Vector3i> cuboid_triangles_vertex_indices{
        Eigen::Vector3i(0, 2, 1), Eigen::Vector3i(0, 1, 4),
        Eigen::Vector3i(0, 4, 2), Eigen::Vector3i(5, 1, 7),
        Eigen::Vector3i(5, 7, 4), Eigen::Vector3i(5, 4, 1),
        Eigen::Vector3i(3, 7, 1), Eigen::Vector3i(3, 1, 2),
        Eigen::Vector3i(3, 2, 7), Eigen::Vector3i(6, 4, 7),
        Eigen::Vector3i(6, 7, 2), Eigen::Vector3i(6, 2, 4),
};

// Vertex indices of 12 lines in a cuboid
const static std::vector<Eigen::Vector2i> cuboid_lines_vertex_indices{
        Eigen::Vector2i(0, 1), Eigen::Vector2i(0, 2), Eigen::Vector2i(0, 4),
        Eigen::Vector2i(3, 1), Eigen::Vector2i(3, 2), Eigen::Vector2i(3, 7),
        Eigen::Vector2i(5, 1), Eigen::Vector2i(5, 4), Eigen::Vector2i(5, 7),
        Eigen::Vector2i(6, 2), Eigen::Vector2i(6, 4), Eigen::Vector2i(6, 7),
};

bool MySimpleShader::Compile() {
    if (!CompileShaders(MyShaders::SimpleVertexShader, NULL, MyShaders::SimpleFragmentShader)) {
        PrintShaderWarning("Compiling shaders failed.");
        return false;
    }
    vertex_position_ = glGetAttribLocation(program_, "vertex_position");
    vertex_color_ = glGetAttribLocation(program_, "vertex_color");
    MVP_ = glGetUniformLocation(program_, "MVP");
    return true;
}

void MySimpleShader::Release() {
    UnbindGeometry();
    ReleaseProgram();
}

bool MySimpleShader::BindGeometry(const geometry::Geometry &geometry,
                                const RenderOption &option,
                                const MyViewControl &view) {
    // If there is already geometry, we first unbind it.
    // We use GL_STATIC_DRAW. When geometry changes, we clear buffers and
    // rebind the geometry. Note that this approach is slow. If the geometry is
    // changing per frame, consider implementing a new ShaderWrapper using
    // GL_STREAM_DRAW, and replace InvalidateGeometry() with Buffer Object
    // Streaming mechanisms.
    UnbindGeometry();

    // Prepare data to be passed to GPU
    std::vector<Eigen::Vector3f> points;
    std::vector<Eigen::Vector3f> colors;
    if (!PrepareBinding(geometry, option, view, points, colors)) {
        PrintShaderWarning("Binding failed when preparing data.");
        return false;
    }

    // Create buffers and bind the geometry
    glGenBuffers(1, &vertex_position_buffer_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_position_buffer_);
    glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(Eigen::Vector3f),
                 points.data(), GL_STATIC_DRAW);
    glGenBuffers(1, &vertex_color_buffer_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_color_buffer_);
    glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(Eigen::Vector3f),
                 colors.data(), GL_STATIC_DRAW);
    bound_ = true;
    return true;
}

bool MySimpleShader::RenderGeometry(const geometry::Geometry &geometry,
                                  const RenderOption &option,
                                  const MyViewControl &view) {
    if (!PrepareRendering(geometry, option, view)) {
        PrintShaderWarning("Rendering failed during preparation.");
        return false;
    }
    glUseProgram(program_);
    glUniformMatrix4fv(MVP_, 1, GL_FALSE, view.GetMVPMatrix().data());
    glEnableVertexAttribArray(vertex_position_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_position_buffer_);
    glVertexAttribPointer(vertex_position_, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(vertex_color_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_color_buffer_);
    glVertexAttribPointer(vertex_color_, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glDrawArrays(draw_arrays_mode_, 0, draw_arrays_size_);
    glDisableVertexAttribArray(vertex_position_);
    glDisableVertexAttribArray(vertex_color_);
    return true;
}

void MySimpleShader::UnbindGeometry() {
    if (bound_) {
        glDeleteBuffers(1, &vertex_position_buffer_);
        glDeleteBuffers(1, &vertex_color_buffer_);
        bound_ = false;
    }
}

};