#include "MyRenderer.h"
#include "MyViewControl.h"

using namespace myvisualization;


bool MyTriangleMeshRenderer::AddGeometry(
        std::shared_ptr<const open3d::geometry::Geometry> geometry_ptr) {

    if (geometry_ptr->GetGeometryType() !=
                open3d::geometry::Geometry::GeometryType::TriangleMesh &&
        geometry_ptr->GetGeometryType() !=
                open3d::geometry::Geometry::GeometryType::HalfEdgeTriangleMesh) {
        return false;
    }
    geometry_ptr_ = geometry_ptr;
    return UpdateGeometry();
}

bool MyTriangleMeshRenderer::Render(const RenderOption &option,
                                  const MyViewControl &view) {
    if (!is_visible_ || geometry_ptr_->IsEmpty()) return true;
    const auto &mesh = (const open3d::geometry::TriangleMesh &)(*geometry_ptr_);
    bool success = true;
    if (mesh.HasTriangleNormals() && mesh.HasVertexNormals()) {
            // Make sure there is no color option or lacking normal data
        success &= phong_mesh_shader_.Render(mesh, option, view);
    } 
    return success;
}

bool MyTriangleMeshRenderer::UpdateGeometry() {
    phong_mesh_shader_.InvalidateGeometry();
    return true;
}

bool MyLineSetRenderer::Render(const RenderOption &option,
                             const MyViewControl &view) {
    if (!is_visible_ || geometry_ptr_->IsEmpty()) return true;
    return simple_lineset_shader_.Render(*geometry_ptr_, option, view);
}

bool MyLineSetRenderer::AddGeometry(
        std::shared_ptr<const geometry::Geometry> geometry_ptr) {
    if (geometry_ptr->GetGeometryType() !=
        geometry::Geometry::GeometryType::LineSet) {
        return false;
    }
    geometry_ptr_ = geometry_ptr;
    return UpdateGeometry();
}

bool MyLineSetRenderer::UpdateGeometry() {
    simple_lineset_shader_.InvalidateGeometry();
    return true;
}

bool MyPointCloudRenderer::Render(const RenderOption &option,
                                const MyViewControl &view) {
    if (!is_visible_ || geometry_ptr_->IsEmpty()) return true;
    const auto &pointcloud = (const geometry::PointCloud &)(*geometry_ptr_);
    bool success = true;
    success &= simple_point_shader_.Render(pointcloud, option, view);
    
    return success;
}

bool MyPointCloudRenderer::AddGeometry(
        std::shared_ptr<const geometry::Geometry> geometry_ptr) {
    if (geometry_ptr->GetGeometryType() !=
        geometry::Geometry::GeometryType::PointCloud) {
        return false;
    }
    geometry_ptr_ = geometry_ptr;
    return UpdateGeometry();
}

bool MyPointCloudRenderer::UpdateGeometry() {
    simple_point_shader_.InvalidateGeometry();
    return true;
}

namespace myglsl {

bool ShaderWrapper::Render(const geometry::Geometry &geometry,
                           const RenderOption &option,
                           const MyViewControl &view) {
    if (!compiled_) {
        Compile();
    }
    if (!bound_) {
        BindGeometry(geometry, option, view);
    }
    if (!compiled_ || !bound_) {
        PrintShaderWarning("Something is wrong in compiling or binding.");
        return false;
    }
    return RenderGeometry(geometry, option, view);
}

void ShaderWrapper::InvalidateGeometry() {
    if (bound_) {
        UnbindGeometry();
    }
}

void ShaderWrapper::PrintShaderWarning(const std::string &message) const {
    utility::LogWarning("[{}] {}", GetShaderName(), message);
}

bool ShaderWrapper::CompileShaders(const char *const vertex_shader_code,
                                   const char *const geometry_shader_code,
                                   const char *const fragment_shader_code) {
    if (compiled_) {
        return true;
    }

    if (vertex_shader_code != NULL) {
        vertex_shader_ = glCreateShader(GL_VERTEX_SHADER);
        const GLchar *vertex_shader_code_buffer = vertex_shader_code;
        glShaderSource(vertex_shader_, 1, &vertex_shader_code_buffer, NULL);
        glCompileShader(vertex_shader_);
        if (!ValidateShader(vertex_shader_)) {
            return false;
        }
    }

    if (geometry_shader_code != NULL) {
        geometry_shader_ = glCreateShader(GL_GEOMETRY_SHADER);
        const GLchar *geometry_shader_code_buffer = geometry_shader_code;
        glShaderSource(geometry_shader_, 1, &geometry_shader_code_buffer, NULL);
        glCompileShader(geometry_shader_);
        if (!ValidateShader(geometry_shader_)) {
            return false;
        }
    }

    if (fragment_shader_code != NULL) {
        fragment_shader_ = glCreateShader(GL_FRAGMENT_SHADER);
        const GLchar *fragment_shader_code_buffer = fragment_shader_code;
        glShaderSource(fragment_shader_, 1, &fragment_shader_code_buffer, NULL);
        glCompileShader(fragment_shader_);
        if (!ValidateShader(fragment_shader_)) {
            return false;
        }
    }

    program_ = glCreateProgram();
    if (vertex_shader_code != NULL) {
        glAttachShader(program_, vertex_shader_);
    }
    if (geometry_shader_code != NULL) {
        glAttachShader(program_, geometry_shader_);
    }
    if (fragment_shader_code != NULL) {
        glAttachShader(program_, fragment_shader_);
    }
    glLinkProgram(program_);
    if (!ValidateProgram(program_)) {
        return false;
    }

    // Mark shader objects as deletable.
    // They will be released as soon as program is deleted.
    if (vertex_shader_code != NULL) {
        glDeleteShader(vertex_shader_);
    }
    if (geometry_shader_code != NULL) {
        glDeleteShader(geometry_shader_);
    }
    if (fragment_shader_code != NULL) {
        glDeleteShader(fragment_shader_);
    }

    compiled_ = true;
    return true;
}

void ShaderWrapper::ReleaseProgram() {
    if (compiled_) {
        glDeleteProgram(program_);
        compiled_ = false;
    }
}

bool ShaderWrapper::ValidateShader(GLuint shader_index) {
    GLint result = GL_FALSE;
    int info_log_length;
    glGetShaderiv(shader_index, GL_COMPILE_STATUS, &result);
    if (result == GL_FALSE) {
        glGetShaderiv(shader_index, GL_INFO_LOG_LENGTH, &info_log_length);
        if (info_log_length > 0) {
            std::vector<char> error_message(info_log_length + 1);
            glGetShaderInfoLog(shader_index, info_log_length, NULL,
                               &error_message[0]);
            utility::LogWarning("Shader error: {}", &error_message[0]);
        }
        return false;
    }
    return true;
}

bool ShaderWrapper::ValidateProgram(GLuint program_index) {
    GLint result = GL_FALSE;
    int info_log_length;
    glGetProgramiv(program_index, GL_LINK_STATUS, &result);
    if (result == GL_FALSE) {
        glGetProgramiv(program_index, GL_INFO_LOG_LENGTH, &info_log_length);
        if (info_log_length > 0) {
            std::vector<char> error_message(info_log_length + 1);
            glGetShaderInfoLog(program_index, info_log_length, NULL,
                               &error_message[0]);
            utility::LogWarning("Shader error: {}", &error_message[0]);
        }
        return false;
    }
    return true;
}




}