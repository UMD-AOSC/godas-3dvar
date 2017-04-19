module g3dv_grid
  !! author: Travis Sluka
  !!
  !! g3dv_grid description
  !!
  !! g3dv_grid detailed description
  !!
  use kdtree
  use timing
  use g3dv_datatable, only : datatable_get
  use g3dv_mpi
  use netcdf


  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: grid_init
  public :: grid_scatter
  public :: grid_write


  ! public module variables
  !------------------------------------------------------------
  integer, public, protected :: grid_nx
  integer, public, protected :: grid_ny
  integer, public, protected :: grid_nz
  integer, public, protected :: grid_ns 
  !! total number of slabs = [(grid_nz * num_3d_vars) + num_2d_vars]

  integer, public, protected :: grid_var_t
  integer, public, protected :: grid_var_s
  integer, public, protected :: grid_var_u
  integer, public, protected :: grid_var_v
  integer, public, protected :: grid_var_ssh

  ! these are only valid on the root proc
  real, public, protected, allocatable :: grid_lat(:,:)
  real, public, protected, allocatable :: grid_lon(:,:)
  real, public, protected, allocatable :: grid_dpth(:)
  real, public, protected, allocatable :: grid_lvls(:,:)
  real, public, protected, allocatable :: grid_ssh(:,:)
  real, public, protected, allocatable :: grid_dens(:,:,:)

  ! the scattered grid parameters
  real, public, protected, allocatable :: grid_local_lat(:) 
  real, public, protected, allocatable :: grid_local_lon(:)
  real, public, protected, allocatable :: grid_local_lvls(:)
  real, public, protected, allocatable :: grid_local_ssh(:) 
  real, public, protected, allocatable :: grid_local_dens(:,:)


  ! private module variables
  !------------------------------------------------------------
  ! these are only valid on the root proc
  type(kd_root) :: ll_kdtree
  integer, allocatable :: ll_kdtree_x(:)
  integer, allocatable :: ll_kdtree_y(:)

  logical :: isroot



contains


  !================================================================================
  !================================================================================



  subroutine grid_init(isroot0, nml)
    logical, intent(in) :: isroot0
    character(len=*), intent(in) :: nml

    real, allocatable :: tree_lons(:), tree_lats(:)

    integer :: x, y
    integer :: unit
    integer :: timer

    namelist /g3dv_grid/ grid_nx, grid_ny, grid_nz
    

    isroot = isroot0


    timer = timer_init(' grid / state', TIMER_SYNC)
    call timer_start(timer)

    if(isroot) then
       print *,new_line('a'),&
            new_line('a'), '------------------------------------------------------------',&
            new_line('a'), 'grid_init() : ',&
            new_line('a'), '------------------------------------------------------------'
    end if

    ! read in parameters from namelist
    open(newunit=unit, file=nml)
    read(unit, nml=g3dv_grid)
    close(unit)
    if(isroot) then
       print g3dv_grid
       print *,""
    end if


    ! TODO, 
    grid_ns = grid_nz * 4 + 1
    grid_var_t = 1
    grid_var_s = grid_nz + grid_var_t
    grid_var_u = grid_nz + grid_var_u
    grid_var_v = grid_nz + grid_var_v
    grid_var_ssh = grid_var_v + 1


    ! some things to be done only on the root node
    !------------------------------------------------------------
    if(isroot) then
       ! read in the grid specifications
       print *, "Reading grid definition parameters..."
       allocate(grid_lat(grid_nx, grid_ny))
       call datatable_get('grid_y', grid_lat)
       allocate(grid_lon(grid_nx, grid_ny))
       call datatable_get('grid_x', grid_lon)
       allocate(grid_dpth(grid_nz))
       call datatable_get('grid_z', grid_dpth)
       allocate(grid_lvls(grid_nx, grid_ny))
       call datatable_get('grid_lvls', grid_lvls)

       print *, ""
       print '(A,I0,A,I0,A,I0)', " Grid size is ",grid_nx," x ",grid_ny," x ",grid_nz
       print *, " lat  range: ",minval(grid_lat),maxval(grid_lat)
       print *, " lon  range: ",minval(grid_lon),maxval(grid_lon)
       print *, "depth range: ",minval(grid_dpth),maxval(grid_dpth)
       
       ! make sure longitude is 0 to 360
!       where (grid_lon < 0) grid_lon = grid_lon + 360

       ! generate kd-tree for fast lookup of grid points given a lat/lon
       print *,""
       print *, "Generating kd-tree for fast lookup of grid x/y given lat/lon..."
       allocate(tree_lons(grid_nx*grid_ny))
       allocate(tree_lats(grid_nx*grid_ny))
       allocate(ll_kdtree_x(grid_nx*grid_ny))
       allocate(ll_kdtree_y(grid_nx*grid_ny))
       do x = 1, grid_nx
          do y = 1, grid_ny
          tree_lons((y-1)*grid_nx + x) = grid_lon(x,y)
          tree_lats((y-1)*grid_nx + x) = grid_lat(x,y)
          ll_kdtree_x((y-1)*grid_nx + x) = x
          ll_kdtree_y((y-1)*grid_nx + x) = y
          end do
       end do
       call kd_init(ll_kdtree, tree_lons, tree_lats)
       deallocate(tree_lons)
       deallocate(tree_lats)


       ! read in state fields
       print *, ""
       print *, "Reading background fields..."
       allocate(grid_ssh(grid_nx, grid_ny))
       call datatable_get('bg_ssh', grid_ssh)
       allocate(grid_dens(grid_nx, grid_ny, grid_nz))
       call datatable_get('bg_dens',grid_dens)
       
       print *, ""
       print *, " ssh      range: ",minval(grid_ssh), maxval(grid_ssh)
       print *, " density  range: ",minval(grid_dens),maxval(grid_dens)
    end if



    call timer_stop(timer)
  end subroutine grid_init




  !================================================================================
  !================================================================================


  
  subroutine grid_scatter()
    integer :: i

    if(isroot) then
       print *,new_line('a'),&
            new_line('a'), '------------------------------------------------------------',&
            new_line('a'), 'grid_scatter() : ',&
            new_line('a'), '------------------------------------------------------------'
    end if

    if (isroot) print *, "scattering LAT ..."
    allocate(grid_local_lat(g3dv_mpi_ijcount))
    call g3dv_mpi_grd2ij_real(grid_lat, grid_local_lat)
    
    if (isroot) print *, "scattering LON ..."
    allocate(grid_local_lon(g3dv_mpi_ijcount))
    call g3dv_mpi_grd2ij_real(grid_lon, grid_local_lon)

    if (isroot) print *, "scattering LVLS ..."
    allocate(grid_local_lvls(g3dv_mpi_ijcount))
    call g3dv_mpi_grd2ij_real(grid_lvls, grid_local_lvls)

    if (isroot) print *, "scattering SSH ..."
    allocate(grid_local_ssh(g3dv_mpi_ijcount))
    call g3dv_mpi_grd2ij_real(grid_ssh, grid_local_ssh)
    if (isroot) deallocate(grid_ssh)

    if (isroot) print *, "scattering DENS ..."
    allocate(grid_local_dens(grid_nz, g3dv_mpi_ijcount))
    do i = 1, grid_nz
       call g3dv_mpi_grd2ij_real(grid_dens(:,:,i), grid_local_dens(i,:))
    end do
    if (isroot) deallocate(grid_dens)

  end subroutine grid_scatter





  !================================================================================
  !================================================================================

  subroutine grid_write(grd, filename)
    real, intent(in) :: grd(grid_nx, grid_ny, grid_ns)
    character(len=*), intent(in) :: filename

    integer :: ncid, vid
    integer :: d_x, d_y, d_z
    integer :: v_x, v_y, v_z
    integer :: v_t, v_s, v_u, v_v, v_ssh 


    print *, "Saving analysis grid to",trim(filename)

    ! setup the output file
    call check(nf90_create(filename, NF90_WRITE, ncid))

    call check(nf90_def_dim(ncid, "grid_x", grid_nx, d_x))
    call check(nf90_def_var(ncid, "grid_x", nf90_real, (/d_x/), v_x))

    call check(nf90_def_dim(ncid, "grid_y", grid_ny, d_y))
    call check(nf90_def_var(ncid, "grid_y", nf90_real, (/d_y/), v_y))

    call check(nf90_def_dim(ncid, "grid_z", grid_nz, d_z))
    call check(nf90_def_var(ncid, "grid_z", nf90_real, (/d_z/), v_z))

    call check(nf90_def_var(ncid, "t", nf90_real, (/d_x, d_y, d_z/), v_t))
    call check(nf90_def_var(ncid, "s", nf90_real, (/d_x, d_y, d_z/), v_s))
    call check(nf90_def_var(ncid, "u", nf90_real, (/d_x, d_y, d_z/), v_u))
    call check(nf90_def_var(ncid, "v", nf90_real, (/d_x, d_y, d_z/), v_v))
    call check(nf90_def_var(ncid, "ssh", nf90_real, (/d_x, d_y/), v_ssh))

    call check(nf90_enddef(ncid))

    ! write out data
    call check(nf90_put_var(ncid, v_x, grid_lon(:,1)))
    call check(nf90_put_var(ncid, v_y, grid_lat(1,:)))
    call check(nf90_put_var(ncid, v_t, grd(:,:,grid_var_t:grid_var_t+grid_nz-1)))
    call check(nf90_put_var(ncid, v_s, grd(:,:,grid_var_s:grid_var_s+grid_nz-1)))
    call check(nf90_put_var(ncid, v_u, grd(:,:,grid_var_u:grid_var_u+grid_nz-1)))
    call check(nf90_put_var(ncid, v_v, grd(:,:,grid_var_v:grid_var_v+grid_nz-1)))
    call check(nf90_put_var(ncid, v_ssh, grd(:,:,grid_var_ssh)))


    call check(nf90_close(ncid))
    
    
  end subroutine grid_write


  subroutine check(status)
    integer, intent(in) :: status
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop 1
    end if
  end subroutine check

end module  g3dv_grid
