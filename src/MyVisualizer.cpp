#include "MyVisualizer.h"
#include "m_util.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <Eigen/Core>

#include <math.h>
#include <limits.h>
#include <algorithm>

#define NOMINMAX

using namespace open3d::visualization;
using namespace open3d;

namespace Temp{
class GLFWContext {
private:
    GLFWContext() {
        utility::LogDebug("GLFW init.");

#if defined(__APPLE__)
        // On macOS, GLFW changes the directory to the resource directory,
        // which will cause an unexpected change of directory if using a
        // framework build version of Python.
        glfwInitHint(GLFW_COCOA_CHDIR_RESOURCES, GLFW_FALSE);
#endif
        init_status_ = glfwInit();
    }

    GLFWContext(const GLFWContext &) = delete;
    GLFWContext &operator=(const GLFWContext &) = delete;

public:
    ~GLFWContext() {
        if (init_status_ == GLFW_TRUE) {
            glfwTerminate();
            utility::LogDebug("GLFW destruct.");
        }
    }

    /// \brief Get the glfwInit status.
    inline int InitStatus() const { return init_status_; }

    /// \brief Get a shared instance of the GLFW context.
    static std::shared_ptr<GLFWContext> GetInstance() {
        static std::weak_ptr<GLFWContext> singleton;

        auto res = singleton.lock();
        if (res == nullptr) {
            res = std::shared_ptr<GLFWContext>(new GLFWContext());
            singleton = res;
        }
        return res;
    }

    static void GLFWErrorCallback(int error, const char *description) {
        utility::LogWarning("GLFW Error: {}", description);
    }

private:
    /// \brief Status of the glfwInit call.
    int init_status_ = GLFW_FALSE;
};
}

namespace myvisualization{

bool MyVisualizer::AddGeometry(std::shared_ptr<open3d::geometry::Geometry> geometry_ptr
,bool reset_bounding_box, bool update_geometry) {
    return AddGeometryA(geometry_ptr, reset_bounding_box, highest_priority+1, update_geometry);
}

bool MyVisualizer::AddGeometryA(
            std::shared_ptr<geometry::Geometry> geometry_ptr,
            bool reset_bounding_box, char target_priority, bool update_geometry){

    std::string ori_name = geometry_ptr->GetName();
    if(ori_name.empty()) ori_name.resize(1);
    ori_name[0] = target_priority;
    geometry_ptr->SetName(ori_name);

    highest_priority = (m_max)(highest_priority, target_priority);

    if (!is_initialized_) {
            return false;
        }
    if (!geometry_ptr.get()) {
        open3d::utility::LogWarning(
                "[AddGeometry] Invalid pointer. Possibly a null pointer or "
                "None was passed in.");
        return false;
    }

    glfwMakeContextCurrent(window_);

    std::shared_ptr<myglsl::MyGeometryRenderer> renderer_ptr;
    if(!current_renderer_name.empty()){
        renderer_ptr = PickRenderer(current_renderer_name);
        if(!renderer_ptr->AddGeometry(geometry_ptr)){
            std:: cout << "Add Custom Geometry Failed" << std::endl;
            return false;
        }
    }
    else if (geometry_ptr->GetGeometryType() ==
        open3d::geometry::Geometry::GeometryType::Unspecified) {
        return false;
    } else if (geometry_ptr->GetGeometryType() ==
               open3d::geometry::Geometry::GeometryType::PointCloud) {
        renderer_ptr = std::make_shared<MyPointCloudRenderer>();
        if (!renderer_ptr->AddGeometry(geometry_ptr)) {
            return false;
        }
    } /*else if (geometry_ptr->GetGeometryType() ==
               open3d::geometry::Geometry::GeometryType::VoxelGrid) {
        renderer_ptr = std::make_shared<glsl::VoxelGridRenderer>();
        if (!renderer_ptr->AddGeometry(geometry_ptr)) {
            return false;
        }
    } else if (geometry_ptr->GetGeometryType() ==
               open3d::geometry::Geometry::GeometryType::Octree) {
        renderer_ptr = std::make_shared<glsl::OctreeRenderer>();
        if (!renderer_ptr->AddGeometry(geometry_ptr)) {
            return false;
        }*/
    else if (geometry_ptr->GetGeometryType() ==
               open3d::geometry::Geometry::GeometryType::LineSet) {
        renderer_ptr = std::make_shared<MyLineSetRenderer>();
        if (!renderer_ptr->AddGeometry(geometry_ptr)) {
            return false;
        }
    } else if (geometry_ptr->GetGeometryType() ==
                       open3d::geometry::Geometry::GeometryType::TriangleMesh ||
               geometry_ptr->GetGeometryType() ==
                       open3d::geometry::Geometry::GeometryType::HalfEdgeTriangleMesh) {
                        
        renderer_ptr = std::make_shared<MyTriangleMeshRenderer>();
        if (!renderer_ptr->AddGeometry(geometry_ptr)) {
            return false;
        }
    } /*else if (geometry_ptr->GetGeometryType() ==
               open3d::geometry::Geometry::GeometryType::Image) {
        renderer_ptr = std::make_shared<glsl::ImageRenderer>();
        if (!renderer_ptr->AddGeometry(geometry_ptr)) {
            return false;
        }
    } else if (geometry_ptr->GetGeometryType() ==
               open3d::geometry::Geometry::GeometryType::RGBDImage) {
        renderer_ptr = std::make_shared<glsl::RGBDImageRenderer>();
        if (!renderer_ptr->AddGeometry(geometry_ptr)) {
            return false;
        }
    } else if (geometry_ptr->GetGeometryType() ==
               open3d::geometry::Geometry::GeometryType::TetraMesh) {
        renderer_ptr = std::make_shared<glsl::TetraMeshRenderer>();
        if (!renderer_ptr->AddGeometry(geometry_ptr)) {
            return false;
        }
    } else if (geometry_ptr->GetGeometryType() ==
               open3d::geometry::Geometry::GeometryType::OrientedBoundingBox) {
        renderer_ptr = std::make_shared<glsl::OrientedBoundingBoxRenderer>();
        if (!renderer_ptr->AddGeometry(geometry_ptr)) {
            return false;
        }
    } else if (geometry_ptr->GetGeometryType() ==
               open3d::geometry::Geometry::GeometryType::AxisAlignedBoundingBox) {
        renderer_ptr = std::make_shared<glsl::AxisAlignedBoundingBoxRenderer>();
        if (!renderer_ptr->AddGeometry(geometry_ptr)) {
            return false;
        }
    } */else {
        return false;
    }

    if (reset_bounding_box) {
        renderer_ptr->SetResetBoundingBox(true);
        view_control_ptr_->FitInGeometry(*geometry_ptr);
    }
    else{
        renderer_ptr->SetResetBoundingBox(false);
    }
    if(update_geometry) {
        renderer_ptr->SetUpdate(true);
    }
    else{
        renderer_ptr->SetUpdate(false);
    }

    geometry_renderer_ptrs_.insert(renderer_ptr);
    geometry_ptrs_.insert(geometry_ptr);
    open3d::utility::LogDebug(
            "Add geometry and update bounding box to {}",
            view_control_ptr_->GetBoundingBox().GetPrintInfo().c_str());
    return UpdateGeometry(geometry_ptr);

}

bool MyVisualizer::SetCurrentRenderer(std::string name){
    current_renderer_name = name;
    return true;
}
bool MyVisualizer::ResetCurrentRenderer(){
    current_renderer_name.clear();
    return true;
}

std::shared_ptr<myglsl::MyGeometryRenderer> MyVisualizer::PickRenderer(std::string name){
    if(name == "MyTriangleRenderer"){
        return std::make_shared<MyTriangleMeshRenderer>();
    }
    return std::shared_ptr<myglsl::MyGeometryRenderer>();
}

bool MyVisualizer::InitOpenGL() {
#if defined(BUILD_GUI)
    // With the current link strategy the OpenGL functions are bound to
    // Filament's BlueGL internal stubs which are initially null. BlueGL loads
    // the 'real' OpenGL functions dynamically. In new visualizer, Filament
    // automatically initializes BlueGL for us, but here we have to do manually
    // otherwise the OpenGL functions will point to null functions and crash.
    if (bluegl::bind()) {
        utility::LogWarning("Visualizer::InitOpenGL: bluegl::bind() error.");
    }
#endif

    glewExperimental = true;
    if (glewInit() != GLEW_OK) {
        utility::LogWarning("Failed to initialize GLEW.");
        return false;
    }

    render_fbo_ = 0;

    glGenVertexArrays(1, &vao_id_);
    glBindVertexArray(vao_id_);

    // depth test
    glEnable(GL_DEPTH_TEST);
    glClearDepth(1.0f);

    // pixel alignment
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    // polygon rendering
    glEnable(GL_CULL_FACE);

    // glReadPixels always read front buffer
    glReadBuffer(GL_FRONT);

    return true;
}

void MyVisualizer::Render(bool render_screen) {

    glfwMakeContextCurrent(window_);

    view_control_ptr_->SetViewMatrices();

    if (render_screen) {

        if (render_fbo_ != 0) {
            utility::LogWarning("Render framebuffer is not released.");
        }

        glGenFramebuffers(1, &render_fbo_);
        glBindFramebuffer(GL_FRAMEBUFFER, render_fbo_);

        int tex_w = view_control_ptr_->GetWindowWidth();
        int tex_h = view_control_ptr_->GetWindowHeight();

        glGenTextures(1, &render_rgb_tex_);
        glBindTexture(GL_TEXTURE_2D, render_rgb_tex_);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex_w, tex_h, 0, GL_RGB,
                     GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                               GL_TEXTURE_2D, render_rgb_tex_, 0);

        glGenRenderbuffers(1, &render_depth_stencil_rbo_);
        glBindRenderbuffer(GL_RENDERBUFFER, render_depth_stencil_rbo_);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, tex_w,
                              tex_h);
        glBindRenderbuffer(GL_RENDERBUFFER, 0);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT,
                                  GL_RENDERBUFFER, render_depth_stencil_rbo_);
    }

    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND); // Translucent
    auto &background_color = render_option_ptr_->background_color_;
    glClearColor((GLclampf)background_color(0), (GLclampf)background_color(1),
                 (GLclampf)background_color(2), 1.0f);
    glClearDepth(1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    for (const auto &renderer_ptr : geometry_renderer_ptrs_) {
        renderer_ptr->Render(*render_option_ptr_, *view_control_ptr_);
    }



    for (const auto &renderer_ptr : utility_renderer_ptrs_) {
        RenderOption *opt = render_option_ptr_.get();
        auto optIt = utility_renderer_opts_.find(renderer_ptr);
        if (optIt != utility_renderer_opts_.end()) {
            opt = &optIt->second;
        }
        renderer_ptr->Render(*opt, *view_control_ptr_);
    }

    glfwSwapBuffers(window_);
}

void MyVisualizer::ResetViewPoint(bool reset_bounding_box /* = false*/) {
    if (reset_bounding_box) {
        view_control_ptr_->ResetBoundingBox();
        for (const auto &geometry_ptr : geometry_ptrs_) {
            view_control_ptr_->FitInGeometry(*(geometry_ptr));
        }
        if (coordinate_frame_mesh_ptr_ && coordinate_frame_mesh_renderer_ptr_) {
            const auto &boundingbox = view_control_ptr_->GetBoundingBox();
            *coordinate_frame_mesh_ptr_ =
                    *geometry::TriangleMesh::CreateCoordinateFrame(
                            boundingbox.GetMaxExtent() * 0.2,
                            boundingbox.min_bound_);
            coordinate_frame_mesh_renderer_ptr_->UpdateGeometry();
        }
    }
    view_control_ptr_->Reset();
    is_redraw_required_ = true;
}





std::shared_ptr<geometry::Image> MyVisualizer::CaptureScreenFloatBuffer(
        bool do_render /* = true*/) {
    geometry::Image screen_image;
    screen_image.Prepare(view_control_ptr_->GetWindowWidth(),
                         view_control_ptr_->GetWindowHeight(), 3, 4);
    if (do_render) {
        Render(true);
        is_redraw_required_ = false;
    }
    glFinish();
    glReadPixels(0, 0, view_control_ptr_->GetWindowWidth(),
                 view_control_ptr_->GetWindowHeight(), GL_RGB, GL_FLOAT,
                 screen_image.data_.data());

    if (render_fbo_ != 0) {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glDeleteFramebuffers(1, &render_fbo_);
        glDeleteRenderbuffers(1, &render_depth_stencil_rbo_);
        glDeleteTextures(1, &render_rgb_tex_);
        render_fbo_ = 0;
    }

    // glReadPixels get the screen in a vertically flipped manner
    // Thus we should flip it back.
    auto image_ptr = std::make_shared<geometry::Image>();
    image_ptr->Prepare(view_control_ptr_->GetWindowWidth(),
                       view_control_ptr_->GetWindowHeight(), 3, 4);
    int bytes_per_line = screen_image.BytesPerLine();
    for (int i = 0; i < screen_image.height_; i++) {
        memcpy(image_ptr->data_.data() + bytes_per_line * i,
               screen_image.data_.data() +
                       bytes_per_line * (screen_image.height_ - i - 1),
               bytes_per_line);
    }
    return image_ptr;
}

void MyVisualizer::CaptureScreenImage(const std::string &filename /* = ""*/,
                                    bool do_render /* = true*/) {
    std::string png_filename = filename;
    std::string camera_filename;
    if (png_filename.empty()) {
        std::string timestamp = utility::GetCurrentTimeStamp();
        png_filename = "ScreenCapture_" + timestamp + ".png";
        camera_filename = "ScreenCamera_" + timestamp + ".json";
    }
    geometry::Image screen_image;
    screen_image.Prepare(view_control_ptr_->GetWindowWidth(),
                         view_control_ptr_->GetWindowHeight(), 3, 1);
    if (do_render) {
        Render(true);
        is_redraw_required_ = false;
    }
    glFinish();
    glReadPixels(0, 0, view_control_ptr_->GetWindowWidth(),
                 view_control_ptr_->GetWindowHeight(), GL_RGB, GL_UNSIGNED_BYTE,
                 screen_image.data_.data());

    if (render_fbo_ != 0) {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glDeleteFramebuffers(1, &render_fbo_);
        glDeleteRenderbuffers(1, &render_depth_stencil_rbo_);
        glDeleteTextures(1, &render_rgb_tex_);
        render_fbo_ = 0;
    }

    // glReadPixels get the screen in a vertically flipped manner
    // Thus we should flip it back.
    geometry::Image png_image;
    png_image.Prepare(view_control_ptr_->GetWindowWidth(),
                      view_control_ptr_->GetWindowHeight(), 3, 1);
    int bytes_per_line = screen_image.BytesPerLine();
    for (int i = 0; i < screen_image.height_; i++) {
        memcpy(png_image.data_.data() + bytes_per_line * i,
               screen_image.data_.data() +
                       bytes_per_line * (screen_image.height_ - i - 1),
               bytes_per_line);
    }

    utility::LogDebug("[Visualizer] Screen capture to {}",
                      png_filename.c_str());
    io::WriteImage(png_filename, png_image);
    if (!camera_filename.empty()) {
        utility::LogDebug("[Visualizer] Screen camera capture to {}",
                          camera_filename.c_str());
    }
}

std::shared_ptr<geometry::Image> MyVisualizer::CaptureDepthFloatBuffer(
        bool do_render /* = true*/) {
    geometry::Image depth_image;
    depth_image.Prepare(view_control_ptr_->GetWindowWidth(),
                        view_control_ptr_->GetWindowHeight(), 1, 4);
    if (do_render) {
        Render();
        is_redraw_required_ = false;
    }
    glFinish();

#if __APPLE__
    // On OSX with Retina display and glfw3, there is a bug with glReadPixels().
    // When using glReadPixels() to read a block of depth data. The data is
    // horizontally stretched (vertically it is fine). This issue is related
    // to GLFW_SAMPLES hint. When it is set to 0 (anti-aliasing disabled),
    // glReadPixels() works fine. See this post for details:
    // http://stackoverflow.com/questions/30608121/glreadpixel-one-pass-vs-looping-through-points
    // The reason of this bug is unknown. The current workaround is to read
    // depth buffer column by column. This is 15~30 times slower than one block
    // reading glReadPixels().
    std::vector<float> float_buffer(depth_image.height_);
    float *p = (float *)depth_image.data_.data();
    for (int j = 0; j < depth_image.width_; j++) {
        glReadPixels(j, 0, 1, depth_image.height_, GL_DEPTH_COMPONENT, GL_FLOAT,
                     float_buffer.data());
        for (int i = 0; i < depth_image.height_; i++) {
            p[i * depth_image.width_ + j] = float_buffer[i];
        }
    }
#else   //__APPLE__
    // By default, glReadPixels read a block of depth buffer.
    glReadPixels(0, 0, depth_image.width_, depth_image.height_,
                 GL_DEPTH_COMPONENT, GL_FLOAT, depth_image.data_.data());
#endif  //__APPLE__

    // glReadPixels get the screen in a vertically flipped manner
    // We should flip it back, and convert it to the correct depth value
    auto image_ptr = std::make_shared<geometry::Image>();
    double z_near = view_control_ptr_->GetZNear();
    double z_far = view_control_ptr_->GetZFar();

    image_ptr->Prepare(view_control_ptr_->GetWindowWidth(),
                       view_control_ptr_->GetWindowHeight(), 1, 4);
    for (int i = 0; i < depth_image.height_; i++) {
        float *p_depth = (float *)(depth_image.data_.data() +
                                   depth_image.BytesPerLine() *
                                           (depth_image.height_ - i - 1));
        float *p_image = (float *)(image_ptr->data_.data() +
                                   image_ptr->BytesPerLine() * i);
        for (int j = 0; j < depth_image.width_; j++) {
            if (p_depth[j] == 1.0) {
                continue;
            }
            double z_depth =
                    2.0 * z_near * z_far /
                    (z_far + z_near -
                     (2.0 * (double)p_depth[j] - 1.0) * (z_far - z_near));
            p_image[j] = (float)z_depth;
        }
    }
    return image_ptr;
}

void MyVisualizer::CaptureDepthImage(const std::string &filename /* = ""*/,
                                   bool do_render /* = true*/,
                                   double depth_scale /* = 1000.0*/) {
    std::string png_filename = filename;
    std::string camera_filename;
    if (png_filename.empty()) {
        std::string timestamp = utility::GetCurrentTimeStamp();
        png_filename = "DepthCapture_" + timestamp + ".png";
        camera_filename = "DepthCamera_" + timestamp + ".json";
    }
    geometry::Image depth_image;
    depth_image.Prepare(view_control_ptr_->GetWindowWidth(),
                        view_control_ptr_->GetWindowHeight(), 1, 4);

    if (do_render) {
        Render();
        is_redraw_required_ = false;
    }
    glFinish();

#if __APPLE__
    // On OSX with Retina display and glfw3, there is a bug with glReadPixels().
    // When using glReadPixels() to read a block of depth data. The data is
    // horizontally stretched (vertically it is fine). This issue is related
    // to GLFW_SAMPLES hint. When it is set to 0 (anti-aliasing disabled),
    // glReadPixels() works fine. See this post for details:
    // http://stackoverflow.com/questions/30608121/glreadpixel-one-pass-vs-looping-through-points
    // The reason of this bug is unknown. The current workaround is to read
    // depth buffer column by column. This is 15~30 times slower than one block
    // reading glReadPixels().
    std::vector<float> float_buffer(depth_image.height_);
    float *p = (float *)depth_image.data_.data();
    for (int j = 0; j < depth_image.width_; j++) {
        glReadPixels(j, 0, 1, depth_image.width_, GL_DEPTH_COMPONENT, GL_FLOAT,
                     float_buffer.data());
        for (int i = 0; i < depth_image.height_; i++) {
            p[i * depth_image.width_ + j] = float_buffer[i];
        }
    }
#else   //__APPLE__
    // By default, glReadPixels read a block of depth buffer.
    glReadPixels(0, 0, depth_image.width_, depth_image.height_,
                 GL_DEPTH_COMPONENT, GL_FLOAT, depth_image.data_.data());
#endif  //__APPLE__

    // glReadPixels get the screen in a vertically flipped manner
    // We should flip it back, and convert it to the correct depth value
    geometry::Image png_image;
    double z_near = view_control_ptr_->GetZNear();
    double z_far = view_control_ptr_->GetZFar();

    png_image.Prepare(view_control_ptr_->GetWindowWidth(),
                      view_control_ptr_->GetWindowHeight(), 1, 2);
    for (int i = 0; i < depth_image.height_; i++) {
        float *p_depth = (float *)(depth_image.data_.data() +
                                   depth_image.BytesPerLine() *
                                           (depth_image.height_ - i - 1));
        uint16_t *p_png = (uint16_t *)(png_image.data_.data() +
                                       png_image.BytesPerLine() * i);
        for (int j = 0; j < depth_image.width_; j++) {
            if (p_depth[j] == 1.0) {
                continue;
            }
            double z_depth =
                    2.0 * z_near * z_far /
                    (z_far + z_near -
                     (2.0 * (double)p_depth[j] - 1.0) * (z_far - z_near));
            p_png[j] = (uint16_t)(m_min)(std::round(depth_scale * z_depth),
                                          (double)INT16_MAX);
        }
    }

    utility::LogDebug("[Visualizer] Depth capture to {}", png_filename.c_str());
    io::WriteImage(png_filename, png_image);
    if (!camera_filename.empty()) {
        utility::LogDebug("[Visualizer] Depth camera capture to {}",
                          camera_filename.c_str());
        camera::PinholeCameraParameters parameter;
    }
}

void MyVisualizer::CaptureDepthPointCloud(
        const std::string &filename /* = ""*/,
        bool do_render /* = true*/,
        bool convert_to_world_coordinate /* = false*/) {
    std::string ply_filename = filename;
    std::string camera_filename;
    if (ply_filename.empty()) {
        std::string timestamp = utility::GetCurrentTimeStamp();
        ply_filename = "DepthCapture_" + timestamp + ".ply";
        camera_filename = "DepthCamera_" + timestamp + ".json";
    }
    geometry::Image depth_image;
    depth_image.Prepare(view_control_ptr_->GetWindowWidth(),
                        view_control_ptr_->GetWindowHeight(), 1, 4);

    if (do_render) {
        Render();
        is_redraw_required_ = false;
    }
    glFinish();

#if __APPLE__
    // On OSX with Retina display and glfw3, there is a bug with glReadPixels().
    // When using glReadPixels() to read a block of depth data. The data is
    // horizontally stretched (vertically it is fine). This issue is related
    // to GLFW_SAMPLES hint. When it is set to 0 (anti-aliasing disabled),
    // glReadPixels() works fine. See this post for details:
    // http://stackoverflow.com/questions/30608121/glreadpixel-one-pass-vs-looping-through-points
    // The reason of this bug is unknown. The current workaround is to read
    // depth buffer column by column. This is 15~30 times slower than one block
    // reading glReadPixels().
    std::vector<float> float_buffer(depth_image.height_);
    float *p = (float *)depth_image.data_.data();
    for (int j = 0; j < depth_image.width_; j++) {
        glReadPixels(j, 0, 1, depth_image.width_, GL_DEPTH_COMPONENT, GL_FLOAT,
                     float_buffer.data());
        for (int i = 0; i < depth_image.height_; i++) {
            p[i * depth_image.width_ + j] = float_buffer[i];
        }
    }
#else   //__APPLE__
    // By default, glReadPixels read a block of depth buffer.
    glReadPixels(0, 0, depth_image.width_, depth_image.height_,
                 GL_DEPTH_COMPONENT, GL_FLOAT, depth_image.data_.data());
#endif  //__APPLE__

    gl_util::GLMatrix4f mvp_matrix;
    if (convert_to_world_coordinate) {
        mvp_matrix = view_control_ptr_->GetMVPMatrix();
    } else {
        mvp_matrix = view_control_ptr_->GetProjectionMatrix();
    }

    // glReadPixels get the screen in a vertically flipped manner
    // We should flip it back, and convert it to the correct depth value
    geometry::PointCloud depth_pointcloud;
    for (int i = 0; i < depth_image.height_; i++) {
        float *p_depth = (float *)(depth_image.data_.data() +
                                   depth_image.BytesPerLine() * i);
        for (int j = 0; j < depth_image.width_; j++) {
            if (p_depth[j] == 1.0) {
                continue;
            }
            depth_pointcloud.points_.push_back(gl_util::Unproject(
                    Eigen::Vector3d(j + 0.5, i + 0.5, p_depth[j]), mvp_matrix,
                    view_control_ptr_->GetWindowWidth(),
                    view_control_ptr_->GetWindowHeight()));
        }
    }

    utility::LogDebug("[Visualizer] Depth point cloud capture to {}",
                      ply_filename.c_str());
    io::WritePointCloud(ply_filename, depth_pointcloud);
    if (!camera_filename.empty()) {
        utility::LogDebug("[Visualizer] Depth camera capture to {}",
                          camera_filename.c_str());
    }
}

void MyVisualizer::CaptureRenderOption(const std::string &filename /* = ""*/) {
    std::string json_filename = filename;
    if (json_filename.empty()) {
        std::string timestamp = utility::GetCurrentTimeStamp();
        json_filename = "RenderOption_" + timestamp + ".json";
    }
    utility::LogDebug("[Visualizer] Render option capture to {}",
                      json_filename.c_str());
    io::WriteIJsonConvertible(json_filename, *render_option_ptr_);
}

MyVisualizer::MyVisualizer() :highest_priority('0'){}

MyVisualizer::~MyVisualizer() {
    DestroyVisualizerWindow();

#if defined(__APPLE__) && defined(BUILD_GUI)
    bluegl::unbind();
#endif
}

bool MyVisualizer::CreateVisualizerWindow(
        const std::string &window_name /* = "Open3D"*/,
        const int width /* = 640*/,
        const int height /* = 480*/,
        const int left /* = 50*/,
        const int top /* = 50*/,
        const bool visible /* = true*/) {
    window_name_ = window_name;
    if (window_) {  // window already created
        utility::LogDebug("[Visualizer] Reusing window.");
        UpdateWindowTitle();
        glfwSetWindowPos(window_, left, top);
        glfwSetWindowSize(window_, width, height);
#ifdef __APPLE__
        glfwSetWindowSize(window_,
                          std::round(width * pixel_to_screen_coordinate_),
                          std::round(height * pixel_to_screen_coordinate_));
        glfwSetWindowPos(window_,
                         std::round(left * pixel_to_screen_coordinate_),
                         std::round(top * pixel_to_screen_coordinate_));
#endif  //__APPLE__
        return true;
    }

    utility::LogDebug("[Visualizer] Creating window.");
    glfwSetErrorCallback(Temp::GLFWContext::GLFWErrorCallback);
    glfw_context_ = Temp::GLFWContext::GetInstance();
    if (glfw_context_->InitStatus() != GLFW_TRUE) {
        utility::LogWarning("Failed to initialize GLFW");
        return false;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
#ifndef HEADLESS_RENDERING
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_VISIBLE, visible ? 1 : 0);

    window_ = glfwCreateWindow(width, height, window_name_.c_str(), NULL, NULL);
    if (!window_) {
        utility::LogWarning("Failed to create window");
        return false;
    }
    glfwSetWindowPos(window_, left, top);
    glfwSetWindowUserPointer(window_, this);

#ifdef __APPLE__
    // Some hacks to get pixel_to_screen_coordinate_
    glfwSetWindowSize(window_, 100, 100);
    glfwSetWindowPos(window_, 100, 100);
    int pixel_width_in_osx, pixel_height_in_osx;
    glfwGetFramebufferSize(window_, &pixel_width_in_osx, &pixel_height_in_osx);
    if (pixel_width_in_osx > 0) {
        pixel_to_screen_coordinate_ = 100.0 / (double)pixel_width_in_osx;
    } else {
        pixel_to_screen_coordinate_ = 1.0;
    }
    glfwSetWindowSize(window_, std::round(width * pixel_to_screen_coordinate_),
                      std::round(height * pixel_to_screen_coordinate_));
    glfwSetWindowPos(window_, std::round(left * pixel_to_screen_coordinate_),
                     std::round(top * pixel_to_screen_coordinate_));
#endif  //__APPLE__

    auto window_refresh_callback = [](GLFWwindow *window) {
        static_cast<MyVisualizer *>(glfwGetWindowUserPointer(window))
                ->WindowRefreshCallback(window);
    };
    glfwSetWindowRefreshCallback(window_, window_refresh_callback);

    auto window_resize_callback = [](GLFWwindow *window, int w, int h) {
        static_cast<MyVisualizer *>(glfwGetWindowUserPointer(window))
                ->WindowResizeCallback(window, w, h);
    };
    glfwSetFramebufferSizeCallback(window_, window_resize_callback);

    auto mouse_move_callback = [](GLFWwindow *window, double x, double y) {
        static_cast<MyVisualizer *>(glfwGetWindowUserPointer(window))
                ->MouseMoveCallback(window, x, y);
    };
    glfwSetCursorPosCallback(window_, mouse_move_callback);

    auto mouse_scroll_callback = [](GLFWwindow *window, double x, double y) {
        static_cast<MyVisualizer *>(glfwGetWindowUserPointer(window))
                ->MouseScrollCallback(window, x, y);
    };
    glfwSetScrollCallback(window_, mouse_scroll_callback);

    auto mouse_button_callback = [](GLFWwindow *window, int button, int action,
                                    int mods) {
        static_cast<MyVisualizer *>(glfwGetWindowUserPointer(window))
                ->MouseButtonCallback(window, button, action, mods);
    };
    glfwSetMouseButtonCallback(window_, mouse_button_callback);

    auto key_press_callback = [](GLFWwindow *window, int key, int scancode,
                                 int action, int mods) {
        static_cast<MyVisualizer *>(glfwGetWindowUserPointer(window))
                ->KeyPressCallback(window, key, scancode, action, mods);
    };
    glfwSetKeyCallback(window_, key_press_callback);

    auto window_close_callback = [](GLFWwindow *window) {
        static_cast<MyVisualizer *>(glfwGetWindowUserPointer(window))
                ->WindowCloseCallback(window);
    };
    glfwSetWindowCloseCallback(window_, window_close_callback);

    glfwMakeContextCurrent(window_);
    glfwSwapInterval(1);

    if (!InitOpenGL()) {
        return false;
    }

    if (!InitViewControl()) {
        return false;
    }

    if (!InitRenderOption()) {
        return false;
    }

    int window_width, window_height;
    glfwGetFramebufferSize(window_, &window_width, &window_height);
    WindowResizeCallback(window_, window_width, window_height);

    UpdateWindowTitle();

    is_initialized_ = true;
    animate_timer = chnow();
    return true;
}

void MyVisualizer::DestroyVisualizerWindow() {
    if (!is_initialized_) {
        return;
    }

    utility::LogDebug("[Visualizer] Destroying window.");
    is_initialized_ = false;
    glDeleteVertexArrays(1, &vao_id_);
    vao_id_ = 0;
    glfwDestroyWindow(window_);
    window_ = nullptr;
    glfw_context_.reset();
}

void MyVisualizer::RegisterAnimationCallback(
        std::function<bool(MyVisualizer *)> callback_func) {
    animation_callback_func_ = callback_func;
}

bool MyVisualizer::InitViewControl() {
    view_control_ptr_ = std::unique_ptr<MyViewControl>(new MyViewControl);
    ResetViewPoint();
    return true;
}

bool MyVisualizer::InitRenderOption() {
    render_option_ptr_ = std::unique_ptr<RenderOption>(new RenderOption);
    return true;
}

void MyVisualizer::UpdateWindowTitle() {
    if (window_ != NULL) {
        glfwSetWindowTitle(window_, window_name_.c_str());
    }
}

void MyVisualizer::BuildUtilities() {
    glfwMakeContextCurrent(window_);

    std::cout << "Build";

    // 0. Build coordinate frame
    const auto boundingbox = GetViewControl().GetBoundingBox();
    double extent = (m_max)(0.01, boundingbox.GetMaxExtent() * 0.2);
    coordinate_frame_mesh_ptr_ = geometry::TriangleMesh::CreateCoordinateFrame(
            extent, boundingbox.min_bound_);
            std::cout << "Build";
    coordinate_frame_mesh_renderer_ptr_ =
            std::make_shared<glsl::CoordinateFrameRenderer>();
            std::cout << "Build";
    if (!coordinate_frame_mesh_renderer_ptr_->AddGeometry(
                coordinate_frame_mesh_ptr_)) {
        return;
    }
    utility_ptrs_.push_back(coordinate_frame_mesh_ptr_);
    //utility_renderer_ptrs_.push_back(coordinate_frame_mesh_renderer_ptr_);
}

void MyVisualizer::Run() {
    //BuildUtilities();

    //std::cout<<"1\n";

    UpdateWindowTitle();

    //std::cout<<"2\n";
    
    while (bool(animation_callback_func_) ? PollEvents() : WaitEvents()) {
        if (bool(animation_callback_func_in_loop_)) {
            if (animation_callback_func_in_loop_(this)) {
                //std::cout<<"3\n";
                UpdateGeometry();
            }
            // Set render flag as dirty anyways, because when we use callback
            // functions, we assume something has been changed in the callback
            // and the redraw event should be triggered.
            UpdateRender();
        }
    }

    
}

void MyVisualizer::Close() {
    glfwSetWindowShouldClose(window_, GL_TRUE);
    utility::LogDebug("[Visualizer] Window closing.");
}

bool MyVisualizer::WaitEvents() {
    if (!is_initialized_) {
        return false;
    }
    glfwMakeContextCurrent(window_);
    if (is_redraw_required_) {
        WindowRefreshCallback(window_);
    }
    animation_callback_func_in_loop_ = animation_callback_func_;
    glfwWaitEvents();
    return !glfwWindowShouldClose(window_);
}

bool MyVisualizer::PollEvents() {
    if (!is_initialized_) {
        return false;
    }

    glfwMakeContextCurrent(window_);

    if (is_redraw_required_) {
        WindowRefreshCallback(window_);
    }

    animation_callback_func_in_loop_ = animation_callback_func_;
    
    glfwPollEvents();

    

    return !glfwWindowShouldClose(window_);
}


bool MyVisualizer::RemoveGeometry(
        std::shared_ptr<const geometry::Geometry> geometry_ptr,
        bool reset_bounding_box) {
    if (!is_initialized_) {
        return false;
    }
    glfwMakeContextCurrent(window_);
    std::shared_ptr<myglsl::MyGeometryRenderer> geometry_renderer_delete = NULL;
    for (auto &geometry_renderer_ptr : geometry_renderer_ptrs_) {
        if (geometry_renderer_ptr->GetGeometry() == geometry_ptr)
            geometry_renderer_delete = geometry_renderer_ptr;
    }
    if (geometry_renderer_delete == NULL) return false;
    geometry_renderer_ptrs_.erase(geometry_renderer_delete);
    geometry_ptrs_.erase(geometry_ptr);
    utility::LogDebug(
            "Remove geometry and update bounding box to {}",
            view_control_ptr_->GetBoundingBox().GetPrintInfo().c_str());
    return UpdateGeometry(geometry_ptr);
}

bool MyVisualizer::ClearGeometries() {
    if (!is_initialized_) {
        return false;
    }
    glfwMakeContextCurrent(window_);
    geometry_renderer_ptrs_.clear();
    geometry_ptrs_.clear();
    return UpdateGeometry();
}

bool MyVisualizer::UpdateGeometry(
        std::shared_ptr<const geometry::Geometry> geometry_ptr) {
    glfwMakeContextCurrent(window_);
    bool success = true;

    if (geometry_ptr == nullptr){
        for (const auto &renderer_ptr : geometry_renderer_ptrs_) {
            if(renderer_ptr->NeedUpdate()){
                success = (success && renderer_ptr->UpdateGeometry());  
            }      
        }
    }
    else
        for (const auto &renderer_ptr : geometry_renderer_ptrs_) {
            
                if(renderer_ptr->HasGeometry(geometry_ptr) && renderer_ptr->NeedUpdate()) 
                    success = (success && renderer_ptr->UpdateGeometry());
            
        }
    UpdateRender();
    return success;
}

void MyVisualizer::UpdateRender() { is_redraw_required_ = true; }

bool MyVisualizer::HasGeometry() const { return !geometry_ptrs_.empty(); }

void MyVisualizer::SetFullScreen(bool fullscreen) {
    if (!fullscreen) {
        glfwSetWindowMonitor(window_, NULL, saved_window_pos_(0),
                             saved_window_pos_(1), saved_window_size_(0),
                             saved_window_size_(1), GLFW_DONT_CARE);
    } else {
        glfwGetWindowSize(window_, &saved_window_size_(0),
                          &saved_window_size_(1));
        glfwGetWindowPos(window_, &saved_window_pos_(0), &saved_window_pos_(1));
        GLFWmonitor *monitor = glfwGetPrimaryMonitor();
        if (const GLFWvidmode *mode = glfwGetVideoMode(monitor)) {
            glfwSetWindowMonitor(window_, monitor, 0, 0, mode->width,
                                 mode->height, mode->refreshRate);
        } else {
            utility::LogError(
                    "Internal error: glfwGetVideoMode returns nullptr.");
        }
    }
}

void MyVisualizer::ToggleFullScreen() {
    if (IsFullScreen()) {
        SetFullScreen(false);
    } else {
        SetFullScreen(true);
    }
}

bool MyVisualizer::IsFullScreen() {
    return glfwGetWindowMonitor(window_) != nullptr;
}

void MyVisualizer::PrintVisualizerHelp() {
    // clang-format off
    utility::LogInfo("  -- Mouse view control --");
    utility::LogInfo("    Left button + drag         : Rotate.");
    utility::LogInfo("    Ctrl + left button + drag  : Translate.");
    utility::LogInfo("    Wheel button + drag        : Translate.");
    utility::LogInfo("    Shift + left button + drag : Roll.");
    utility::LogInfo("    Wheel                      : Zoom in/out.");
    utility::LogInfo("");
    utility::LogInfo("  -- Keyboard view control --");
    utility::LogInfo("    [/]          : Increase/decrease field of view.");
    utility::LogInfo("    R            : Reset view point.");
    utility::LogInfo("    Ctrl/Cmd + C : Copy current view status into the clipboard.");
    utility::LogInfo("    Ctrl/Cmd + V : Paste view status from clipboard.");
    utility::LogInfo("");
    utility::LogInfo("  -- General control --");
    utility::LogInfo("    Q, Esc       : Exit window.");
    utility::LogInfo("    H            : Print help message.");
    utility::LogInfo("    P, PrtScn    : Take a screen capture.");
    utility::LogInfo("    D            : Take a depth capture.");
    utility::LogInfo("    O            : Take a capture of current rendering settings.");
    utility::LogInfo("    Alt + Enter  : Toggle between full screen and windowed mode.");
    utility::LogInfo("");
    utility::LogInfo("  -- Render mode control --");
    utility::LogInfo("    L            : Turn on/off lighting.");
    utility::LogInfo("    +/-          : Increase/decrease point size.");
    utility::LogInfo("    Ctrl + +/-   : Increase/decrease width of geometry::LineSet.");
    utility::LogInfo("    N            : Turn on/off point cloud normal rendering.");
    utility::LogInfo("    S            : Toggle between mesh flat shading and smooth shading.");
    utility::LogInfo("    W            : Turn on/off mesh wireframe.");
    utility::LogInfo("    B            : Turn on/off back face rendering.");
    utility::LogInfo("    I            : Turn on/off image zoom in interpolation.");
    utility::LogInfo("    T            : Toggle among image render:");
    utility::LogInfo("                   no stretch / keep ratio / freely stretch.");
    utility::LogInfo("");
    utility::LogInfo("  -- Color control --");
    utility::LogInfo("    0..4,9       : Set point cloud color option.");
    utility::LogInfo("                   0 - Default behavior, render point color.");
    utility::LogInfo("                   1 - Render point color.");
    utility::LogInfo("                   2 - x coordinate as color.");
    utility::LogInfo("                   3 - y coordinate as color.");
    utility::LogInfo("                   4 - z coordinate as color.");
    utility::LogInfo("                   9 - normal as color.");
    utility::LogInfo("    Ctrl + 0..4,9: Set mesh color option.");
    utility::LogInfo("                   0 - Default behavior, render uniform gray color.");
    utility::LogInfo("                   1 - Render point color.");
    utility::LogInfo("                   2 - x coordinate as color.");
    utility::LogInfo("                   3 - y coordinate as color.");
    utility::LogInfo("                   4 - z coordinate as color.");
    utility::LogInfo("                   9 - normal as color.");
    utility::LogInfo("    Shift + 0..4 : Color map options.");
    utility::LogInfo("                   0 - Gray scale color.");
    utility::LogInfo("                   1 - JET color map.");
    utility::LogInfo("                   2 - SUMMER color map.");
    utility::LogInfo("                   3 - WINTER color map.");
    utility::LogInfo("                   4 - HOT color map.");
    utility::LogInfo("");
    // clang-format on
}

void MyVisualizer::WindowRefreshCallback(GLFWwindow *window) {
    if (is_redraw_required_) {
        Render();
        is_redraw_required_ = false;
    }
}

void MyVisualizer::WindowResizeCallback(GLFWwindow *window, int w, int h) {
    view_control_ptr_->ChangeWindowSize(w, h);
    is_redraw_required_ = true;
}

void MyVisualizer::MouseMoveCallback(GLFWwindow *window, double x, double y) {
#ifdef __APPLE__
    x /= pixel_to_screen_coordinate_;
    y /= pixel_to_screen_coordinate_;
#endif
    if (mouse_control_.is_mouse_left_button_down) {
        if (mouse_control_.is_control_key_down) {
            view_control_ptr_->Translate(x - mouse_control_.mouse_position_x,
                                         y - mouse_control_.mouse_position_y,
                                         mouse_control_.mouse_position_x,
                                         mouse_control_.mouse_position_y);
        } else if (mouse_control_.is_shift_key_down) {
            view_control_ptr_->Roll(x - mouse_control_.mouse_position_x);
        } else if (mouse_control_.is_alt_key_down) {
            view_control_ptr_->CameraLocalRotate(
                    x - mouse_control_.mouse_position_x,
                    y - mouse_control_.mouse_position_y,
                    mouse_control_.mouse_position_x,
                    mouse_control_.mouse_position_y);
        } else {
            view_control_ptr_->Rotate(x - mouse_control_.mouse_position_x,
                                      y - mouse_control_.mouse_position_y,
                                      mouse_control_.mouse_position_x,
                                      mouse_control_.mouse_position_y);
        }
        is_redraw_required_ = true;
    }
    if (mouse_control_.is_mouse_middle_button_down) {
        view_control_ptr_->Translate(x - mouse_control_.mouse_position_x,
                                     y - mouse_control_.mouse_position_y,
                                     mouse_control_.mouse_position_x,
                                     mouse_control_.mouse_position_y);
        is_redraw_required_ = true;
    }
    mouse_control_.mouse_position_x = x;
    mouse_control_.mouse_position_y = y;
}

void MyVisualizer::MouseScrollCallback(GLFWwindow *window, double x, double y) {
    view_control_ptr_->Scale(y);
    is_redraw_required_ = true;
}

void MyVisualizer::MouseButtonCallback(GLFWwindow *window,
                                     int button,
                                     int action,
                                     int mods) {
    double x, y;
    glfwGetCursorPos(window, &x, &y);
#ifdef __APPLE__
    x /= pixel_to_screen_coordinate_;
    y /= pixel_to_screen_coordinate_;
#endif
    mouse_control_.mouse_position_x = x;
    mouse_control_.mouse_position_y = y;
    if (action == GLFW_PRESS) {
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            view_control_ptr_->ResetCameraLocalRotate();

            mouse_control_.is_mouse_left_button_down = true;
            mouse_control_.is_control_key_down = (mods & GLFW_MOD_CONTROL) != 0;
            mouse_control_.is_shift_key_down = (mods & GLFW_MOD_SHIFT) != 0;
            mouse_control_.is_alt_key_down = (mods & GLFW_MOD_ALT) != 0;
            mouse_control_.is_super_key_down = (mods & GLFW_MOD_SUPER) != 0;
        } else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
            mouse_control_.is_mouse_middle_button_down = true;
        }
    } else {
        mouse_control_.is_mouse_left_button_down = false;
        mouse_control_.is_mouse_middle_button_down = false;
        mouse_control_.is_control_key_down = false;
        mouse_control_.is_shift_key_down = false;
        mouse_control_.is_alt_key_down = false;
        mouse_control_.is_super_key_down = false;
    }
}

void MyVisualizer::KeyPressCallback(
        GLFWwindow *window, int key, int scancode, int action, int mods) {
    bool is_release = false;

    if (action == GLFW_RELEASE) {
        animate_back = false;
        if(key == GLFW_KEY_DOWN){
            key_down_press = true;
        }
        return;
    }

    switch (key) {
        case GLFW_KEY_SPACE:
            AnimationPause();
            break;
        case GLFW_KEY_LEFT:
            animate_back = true;
            SetAnimationFrame(GetAnimationFrame()-1);
            if (bool(animation_callback_func_in_loop_)) {
                if (animation_callback_func_in_loop_(this)) {
                    UpdateGeometry();
                }
                UpdateRender();
            }
            break;
        case GLFW_KEY_RIGHT:
            animate_back = true;
            SetAnimationFrame(GetAnimationFrame()+1);
            if (bool(animation_callback_func_in_loop_)) {
                if (animation_callback_func_in_loop_(this)) {
                    UpdateGeometry();
                }
                UpdateRender();
            }
            break;
        case GLFW_KEY_UP:
            view_control_ptr_->ResetBoundingBox();
            for (const auto &renderer_ptr : geometry_renderer_ptrs_) {
                if(renderer_ptr->NeedResetBoundingBox())
                    view_control_ptr_->FitInGeometry(*renderer_ptr->GetGeometry()); // 每个循环重新计算摄像机位置
            }
            break;

        case GLFW_KEY_ENTER:
            if (mods & GLFW_MOD_ALT) {
                if (IsFullScreen()) {
                    SetFullScreen(false);
                } else {
                    SetFullScreen(true);
                };
            }
            break;
        case GLFW_KEY_LEFT_BRACKET:
            view_control_ptr_->ChangeFieldOfView(-1.0);
            utility::LogDebug("[Visualizer] Field of view set to {:.2f}.",
                              view_control_ptr_->GetFieldOfView());
            break;
        case GLFW_KEY_RIGHT_BRACKET:
            view_control_ptr_->ChangeFieldOfView(1.0);
            utility::LogDebug("[Visualizer] Field of view set to {:.2f}.",
                              view_control_ptr_->GetFieldOfView());
            break;
        case GLFW_KEY_R:
            ResetViewPoint();
            utility::LogDebug("[Visualizer] Reset view point.");
            break;
        case GLFW_KEY_ESCAPE:
        case GLFW_KEY_Q:
            Close();
            break;
        case GLFW_KEY_H:
            PrintVisualizerHelp();
            break;
        case GLFW_KEY_P:
        case GLFW_KEY_PRINT_SCREEN:
            CaptureScreenImage();
            break;
        case GLFW_KEY_D:
            CaptureDepthImage();
            break;
        case GLFW_KEY_O:
            CaptureRenderOption();
            break;
        case GLFW_KEY_L:
            render_option_ptr_->ToggleLightOn();
            utility::LogDebug("[Visualizer] Lighting {}.",
                              render_option_ptr_->light_on_ ? "ON" : "OFF");
            break;
        case GLFW_KEY_EQUAL:
            if (mods & GLFW_MOD_SHIFT) {
                render_option_ptr_->ChangeLineWidth(1.0);
                utility::LogDebug("[Visualizer] Line width set to {:.2f}.",
                                  render_option_ptr_->line_width_);
            } else {
                render_option_ptr_->ChangePointSize(1.0);
                if (render_option_ptr_->point_show_normal_) {
                    UpdateGeometry();
                }
                utility::LogDebug("[Visualizer] Point size set to {:.2f}.",
                                  render_option_ptr_->point_size_);
            }
            break;
        case GLFW_KEY_MINUS:
            if (mods & GLFW_MOD_SHIFT) {
                render_option_ptr_->ChangeLineWidth(-1.0);
                utility::LogDebug("[Visualizer] Line width set to {:.2f}.",
                                  render_option_ptr_->line_width_);
            } else {
                render_option_ptr_->ChangePointSize(-1.0);
                if (render_option_ptr_->point_show_normal_) {
                    UpdateGeometry();
                }
                utility::LogDebug("[Visualizer] Point size set to {:.2f}.",
                                  render_option_ptr_->point_size_);
            }
            break;
        case GLFW_KEY_N:
            render_option_ptr_->TogglePointShowNormal();
            if (render_option_ptr_->point_show_normal_) {
                UpdateGeometry();
            }
            utility::LogDebug(
                    "[Visualizer] Point normal rendering {}.",
                    render_option_ptr_->point_show_normal_ ? "ON" : "OFF");
            break;
        case GLFW_KEY_S:
            render_option_ptr_->ToggleShadingOption();
            UpdateGeometry();
            utility::LogDebug(
                    "[Visualizer] Mesh shading mode is {}.",
                    render_option_ptr_->mesh_shade_option_ ==
                                    RenderOption::MeshShadeOption::FlatShade
                            ? "FLAT"
                            : "SMOOTH");
            break;
        case GLFW_KEY_W:
            render_option_ptr_->ToggleMeshShowWireframe();
            utility::LogDebug(
                    "[Visualizer] Mesh wireframe rendering {}.",
                    render_option_ptr_->mesh_show_wireframe_ ? "ON" : "OFF");
            break;
        case GLFW_KEY_B:
            render_option_ptr_->ToggleMeshShowBackFace();
            utility::LogDebug(
                    "[Visualizer] Mesh back face rendering {}.",
                    render_option_ptr_->mesh_show_back_face_ ? "ON" : "OFF");
            break;
        case GLFW_KEY_I:
            render_option_ptr_->ToggleInterpolationOption();
            UpdateGeometry();
            utility::LogDebug(
                    "[Visualizer] geometry::Image interpolation mode is {}.",
                    render_option_ptr_->interpolation_option_ ==
                                    RenderOption::TextureInterpolationOption::
                                            Nearest
                            ? "NEARST"
                            : "LINEAR");
            break;
        case GLFW_KEY_T:
            render_option_ptr_->ToggleImageStretchOption();
            utility::LogDebug(
                    "[Visualizer] geometry::Image stretch mode is #{}.",
                    int(render_option_ptr_->image_stretch_option_));
            break;
        case GLFW_KEY_0:
            if (mods & GLFW_MOD_CONTROL) {
                render_option_ptr_->mesh_color_option_ =
                        RenderOption::MeshColorOption::Default;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Mesh color set to DEFAULT.");
            } else if (mods & GLFW_MOD_SHIFT) {
                SetGlobalColorMap(ColorMap::ColorMapOption::Gray);
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Color map set to GRAY.");
            } else {
                render_option_ptr_->point_color_option_ =
                        RenderOption::PointColorOption::Default;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Point color set to DEFAULT.");
            }
            break;
        case GLFW_KEY_1:
            if (mods & GLFW_MOD_CONTROL) {
                render_option_ptr_->mesh_color_option_ =
                        RenderOption::MeshColorOption::Color;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Mesh color set to COLOR.");
            } else if (mods & GLFW_MOD_SHIFT) {
                SetGlobalColorMap(ColorMap::ColorMapOption::Jet);
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Color map set to JET.");
            } else {
                render_option_ptr_->point_color_option_ =
                        RenderOption::PointColorOption::Color;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Point color set to COLOR.");
            }
            break;
        case GLFW_KEY_2:
            if (mods & GLFW_MOD_CONTROL) {
                render_option_ptr_->mesh_color_option_ =
                        RenderOption::MeshColorOption::XCoordinate;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Mesh color set to X.");
            } else if (mods & GLFW_MOD_SHIFT) {
                SetGlobalColorMap(ColorMap::ColorMapOption::Summer);
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Color map set to SUMMER.");
            } else {
                render_option_ptr_->point_color_option_ =
                        RenderOption::PointColorOption::XCoordinate;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Point color set to X.");
            }
            break;
        case GLFW_KEY_3:
            if (mods & GLFW_MOD_CONTROL) {
                render_option_ptr_->mesh_color_option_ =
                        RenderOption::MeshColorOption::YCoordinate;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Mesh color set to Y.");
            } else if (mods & GLFW_MOD_SHIFT) {
                SetGlobalColorMap(ColorMap::ColorMapOption::Winter);
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Color map set to WINTER.");
            } else {
                render_option_ptr_->point_color_option_ =
                        RenderOption::PointColorOption::YCoordinate;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Point color set to Y.");
            }
            break;
        case GLFW_KEY_4:
            if (mods & GLFW_MOD_CONTROL) {
                render_option_ptr_->mesh_color_option_ =
                        RenderOption::MeshColorOption::ZCoordinate;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Mesh color set to Z.");
            } else if (mods & GLFW_MOD_SHIFT) {
                SetGlobalColorMap(ColorMap::ColorMapOption::Hot);
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Color map set to HOT.");
            } else {
                render_option_ptr_->point_color_option_ =
                        RenderOption::PointColorOption::ZCoordinate;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Point color set to Z.");
            }
            break;
        case GLFW_KEY_9:
            if (mods & GLFW_MOD_CONTROL) {
                render_option_ptr_->mesh_color_option_ =
                        RenderOption::MeshColorOption::Normal;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Mesh color set to NORMAL.");
            } else if (mods & GLFW_MOD_SHIFT) {
            } else {
                render_option_ptr_->point_color_option_ =
                        RenderOption::PointColorOption::Normal;
                UpdateGeometry();
                utility::LogDebug("[Visualizer] Point color set to NORMAL.");
            }
            break;
        default:
            break;
    }

    is_redraw_required_ = true;
}

void MyVisualizer::WindowCloseCallback(GLFWwindow *window) {
    // happens when user click the close icon to close the window
}

}