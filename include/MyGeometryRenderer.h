// ----------------------------------------------------------------------------
// -                        Open3D: www.open3d.org                            -
// ----------------------------------------------------------------------------
// Copyright (c) 2018-2023 www.open3d.org
// SPDX-License-Identifier: MIT
// ----------------------------------------------------------------------------

#pragma once

#include "open3d/geometry/Geometry.h"
#include "open3d/visualization/shader/ImageMaskShader.h"
#include "open3d/visualization/shader/ImageShader.h"
#include "open3d/visualization/shader/NormalShader.h"
#include "open3d/visualization/shader/PhongShader.h"
#include "open3d/visualization/shader/PickingShader.h"
#include "open3d/visualization/shader/RGBDImageShader.h"
#include "open3d/visualization/shader/Simple2DShader.h"
#include "open3d/visualization/shader/SimpleBlackShader.h"
#include "open3d/visualization/shader/SimpleShader.h"
#include "open3d/visualization/shader/TexturePhongShader.h"
#include "open3d/visualization/shader/TextureSimpleShader.h"
#include "MyViewControl.h"

using namespace myvisualization;

using namespace open3d::visualization;

namespace myglsl {

class MyGeometryRenderer {
public:
    virtual ~MyGeometryRenderer() {}

public:
    virtual bool Render(const RenderOption &option,
                        const MyViewControl &view) = 0;

    /// Function to add open3d::geometry to the renderer
    /// 1. After calling the function, the renderer owns the open3d::geometry object.
    /// 2. This function returns FALSE if the open3d::geometry type is not matched to
    /// the renderer.
    /// 3. If an added open3d::geometry is changed, programmer must call
    /// UpdateGeometry() to notify the renderer.
    virtual bool AddGeometry(
            std::shared_ptr<const open3d::geometry::Geometry> geometry_ptr) = 0;

    /// Function to update open3d::geometry
    /// Programmer must call this function to notify a change of the open3d::geometry
    virtual bool UpdateGeometry() = 0;

    bool HasGeometry() const { return bool(geometry_ptr_); }
    std::shared_ptr<const open3d::geometry::Geometry> GetGeometry() const {
        return geometry_ptr_;
    }

    bool HasGeometry(
            std::shared_ptr<const open3d::geometry::Geometry> geometry_ptr) const {
        return geometry_ptr_ == geometry_ptr;
    }

    inline bool NeedUpdate(){
        if(auto_update) return true;
        if(require_update){
            require_update = false;
            return true;
        }
        return false;
    }
    inline void SetUpdate(bool option){
        auto_update = option;
    }

    inline void SetResetBoundingBox(bool option){
        reset_bouding_box = option;
    }
    inline bool NeedResetBoundingBox(){
        return reset_bouding_box;
    }
    

    bool IsVisible() const { return is_visible_; }
    void SetVisible(bool visible) { is_visible_ = visible; };

protected:
    std::shared_ptr<const open3d::geometry::Geometry> geometry_ptr_;
    bool is_visible_ = true;

    bool auto_update = true;
    bool require_update = false;
    bool reset_bouding_box = false;
};


}  // namespace myglsl
