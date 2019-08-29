#--------------------------------------------------------------------------------
# Hybrid-GODAS 3dvar main driver
#--------------------------------------------------------------------------------

# source files
set(SRCS
  godas_driver.f90
)

set( EXE 3dvardriver )
add_executable(${EXE} ${SRCS})
target_link_libraries(${EXE} 3dvar)
target_link_libraries(${EXE} ${LAPACK_LIBRARIES})
target_link_libraries(${EXE} gsw)
