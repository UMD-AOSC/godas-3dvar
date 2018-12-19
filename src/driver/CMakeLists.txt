project(3dvar_driver Fortran)

set( EXE godas_3dvar )


# required packages
#------------------------------------------------------------
#include_directories(${NETCDF_INCLUDE_DIR})
#include_directories(${CMAKE_BINARY_DIR}/3dvar_module)


#TODO, allow selection of static vs shared
set(3DVAR_LIB_STATIC ${CMAKE_BINARY_DIR}/lib/lib3dvar.a)
set(3DVAR_LIB_SHARED ${CMAKE_BINARY_DIR}/lib/lib3dvar.so)
set(3DVAR_LIB ${3DVAR_LIB_STATIC})


# build target
# ------------------------------------------------------------
file( GLOB SRCS "*.[f,F]90")
add_executable(${EXE} ${SRCS})
target_link_libraries(${EXE} ${NETCDF_LIBRARIES})
target_link_libraries(${EXE} ${3DVAR_LIB})
target_link_libraries(${EXE} ${MPI_Fortran_LIBRARIES})
target_link_libraries(${EXE} ${LAPACK_LIBRARIES})
