cmake_minimum_required(VERSION 3.18)

project(MyOpen3D LANGUAGES C CXX)

if(__APPLE__)

    set(Open3D_DIR "~/Codes/open3d/open3d-darwin/lib/cmake/Open3D")


    find_package(Open3D REQUIRED)

    set(GLEW_DIR "/opt/homebrew/Cellar/glew/2.2.0_1/lib/cmake/glew/")
    set(GLFW_DIR "/opt/homebrew/Cellar/glfw/3.4/lib/cmake/glfw3/")

    find_package(GLEW REQUIRED)
    find_package(GLFW3 REQUIRED)

    add_executable(Draw)

    file(GLOB_RECURSE SOURCE_H "${PROJECT_SOURCE_DIR}/include/*.h")
    file(GLOB_RECURSE SOURCE_S "${PROJECT_SOURCE_DIR}/src/*.cpp")

    message("Found ${SOURCE_H}")
    message("Found ${SOURCE_S}")

    target_sources(Draw PRIVATE 
        ${SOURCE_H}
        ${SOURCE_S}
        Draw.cpp
        )
    target_include_directories(Draw PRIVATE "${PROJECT_SOURCE_DIR}/include/")


    if(GLEW_FOUND)
        include_directories(${GLEW_INCLUDE_DIRS})
        target_link_libraries(Draw PRIVATE GLEW::GLEW)
    else()
        message(FATAL_ERROR "GLEW not found")
    endif()

    if(GLFW3_FOUND)
        target_link_libraries(Draw PRIVATE glfw)
    else()
        message(FATAL_ERROR "GLFW3 not found")
    endif()

    target_link_libraries(Draw PRIVATE Open3D::Open3D)

    message(STATUS "lib to${Open3D_LIBRARIES}")

endif()

if(WIN32)

    option(STATIC_WINDOWS_RUNTIME "Use static (MT/MTd) Windows runtime" ON)
    if(STATIC_WINDOWS_RUNTIME)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    else()
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
    endif()

    set(BUILD_SHARED_LIBS ON)
    set(Open3D_ROOT "${PROJECT_SOURCE_DIR}/install")
    set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/install")



    find_package(Open3D REQUIRED)

    add_executable(Draw)

    file(GLOB_RECURSE SOURCE_H "${PROJECT_SOURCE_DIR}/include/*.h")
    file(GLOB_RECURSE SOURCE_S "${PROJECT_SOURCE_DIR}/src/*.cpp")

    message("Found ${SOURCE_H}")
    message("Found ${SOURCE_S}")

    target_sources(Draw PRIVATE 
        ${SOURCE_H}
        ${SOURCE_S}
        Draw.cpp
        )
    target_include_directories(Draw PRIVATE "${CMAKE_INSTALL_PREFIX}/include/")
    target_include_directories(Draw PRIVATE "${PROJECT_SOURCE_DIR}/include/")

    target_link_libraries(Draw PRIVATE Open3D::Open3D)

    target_link_libraries(Draw PRIVATE "${CMAKE_INSTALL_PREFIX}/lib/glew32s.lib")
    target_link_libraries(Draw PRIVATE "${CMAKE_INSTALL_PREFIX}/lib/OpenGL32.lib")
    target_link_libraries(Draw PRIVATE "${CMAKE_INSTALL_PREFIX}/lib/glfw3.lib")
    #target_link_libraries(Draw PRIVATE "${CMAKE_INSTALL_PREFIX}/lib/glew32.lib")
    #target_link_libraries(Draw PRIVATE "${CMAKE_INSTALL_PREFIX}/lib/Open3D_3rdparty_glew.lib")




    get_target_property(open3d_type Open3D::Open3D TYPE)
    if(open3d_type STREQUAL "SHARED_LIBRARY")
        message(STATUS "Copying Open3D.dll to ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>")
        add_custom_command(TARGET Draw POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy
                                ${CMAKE_INSTALL_PREFIX}/bin/Open3D.dll
                                ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>)
    endif()
endif()
