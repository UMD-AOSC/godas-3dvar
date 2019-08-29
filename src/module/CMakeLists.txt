#--------------------------------------------------------------------------------
# Hybrid-GODAS 3dvar add_library
#--------------------------------------------------------------------------------

# source files
set(SRCS
  g3dv_bgcov.f90
  g3dv_datatable.f90
  g3dv_grid.f90
  g3dv_mpi.F90
  g3dv_obs.F90
  g3dv_obs_dat.f90
  g3dv_obs_nc.f90
  g3dv_solver.F90
  g3dv.F90
  qsort.f90
  timing.f90
  ../geoKdTree/kdtree.f90
)

# statically linked library
add_library( 3dvar  ${SRCS} ${KDTREE_SRC})

# exernal dependencies
#-------------------------------------------------------------------------------
# MPI
include_directories(${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(3dvar PRIVATE ${MPI_Fortran_LIBRARIES})

#NetCDF
include_directories(${NETCDF_INCLUDE_DIRS})
target_link_libraries(3dvar PRIVATE ${NETCDF_LIBRARIES})

# LAPACK
target_link_libraries(3dvar  PRIVATE ${LAPACK_LIBRARIES})
