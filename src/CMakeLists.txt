# Gibbs Sewater toolbox for fortran
# build only if it hasn't alredy been built by cmake
if( NOT TARGET gsw)
  file(GLOB GSW_TB_SRC "GSW-Fortran/toolbox/*.[f,F]90")
  file(GLOB GSW_MD_SRC "GSW-Fortran/modules/*.[f,F]90")
  add_library(gsw STATIC ${GSW_TB_SRC} ${GSW_MD_SRC})
  set_target_properties(gsw PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/gsw")
  set(gsw_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/gsw
    CACHE INTERNAL "gsw-fortran include directories for modules" FORCE)
endif()

set(KDTREE_SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/geoKdTree/kdtree.f90
  PARENT_SCOPE)
add_subdirectory(module)
add_subdirectory(driver)
