project(3dvar_module Fortran)

set( LIB_NAME 3dvar)


# required packages
#------------------------------------------------------------
include_directories(${NETCDF_INCLUDES})
include_directories(${MPI_Fortran_INCLUDE_PATH})


# build target
#------------------------------------------------------------
file( GLOB SRCS "*.[f,F]90" )
# add_library( ${LIB_NAME}_shared SHARED ${SRCS} )
# set_target_properties( ${LIB_NAME}_shared
#   PROPERTIES
#   OUTPUT_NAME ${LIB_NAME})

file( GLOB KDTREE_SRC "../geoKdTree/*.[f,F]90" )
add_library( ${LIB_NAME}_static STATIC ${SRCS} ${KDTREE_SRC})
set_target_properties( ${LIB_NAME}_static
  PROPERTIES
  OUTPUT_NAME ${LIB_NAME})
		       

