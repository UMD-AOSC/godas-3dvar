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

  use mpi

  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: grid_init

  public :: grid_genweights
  public :: grid_interp2d
  public :: grid_interp3d

  public :: grid_interpweights
  type grid_interpweights
     integer :: x(4)
     integer :: y(4)
     integer :: z
     real    :: w_h(4)
     real    :: w_v
  end type grid_interpweights



  ! public module variables
  !------------------------------------------------------------
  integer, public, protected :: grid_mpi_interpweights

  integer, public, protected :: grid_nx
  integer, public, protected :: grid_ny
  integer, public, protected :: grid_nz
  integer, public, protected :: grid_ns 
  !! total number of slabs = [(grid_nz * num_3d_vars) + num_2d_vars]

  integer, public, protected :: grid_var_t
  integer, public, protected :: grid_var_s
!  integer, public, protected :: grid_var_u
!  integer, public, protected :: grid_var_v
!  integer, public, protected :: grid_var_ssh

  real, public, protected, allocatable :: grid_dpth(:)

  ! the scattered grid parameters
  real, public, protected, allocatable :: grid_local_lat(:) 
  real, public, protected, allocatable :: grid_local_lon(:)
  real, public, protected, allocatable :: grid_local_mask(:)
  real, public, protected, allocatable :: grid_local_ssh(:) 
  real, public, protected, allocatable :: grid_local_coastdist(:) 
  real, public, protected, allocatable :: grid_local_dens(:,:)
  real, public, protected, allocatable :: grid_local_D(:)

  real, public, allocatable :: grid_local_diag3D_2(:,:)



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


  pure function grid_interp2d(grd, w) result(v)
    real, intent(in) :: grd(:,:)
    type(grid_interpweights), intent(in) :: w
    real :: v
    integer :: i
    v = 0.0
    do i = 1, 1
       if (w%w_h(i) == 0) cycle
       v = v + grd(w%x(i), w%y(i)) * w%w_h(i)
    end do
  end function grid_interp2d


  
  !================================================================================
  !================================================================================


  pure function grid_interp3d(grd, w) result(v)
    real, intent(in) :: grd(:,:,:)
    type(grid_interpweights), intent(in) :: w
    real :: v

    integer :: i
    v = 0.0
    do i = 1, 1
       if (w%w_h(i) == 0) cycle
       v = v + grd(w%x(i), w%y(i), w%z) * w%w_h(i)
    end do

    !TODO, this doesn't do z properly
  end function grid_interp3d


  !================================================================================
  !================================================================================



  subroutine grid_init(isroot0, nml)
    logical, intent(in) :: isroot0
    character(len=*), intent(in) :: nml

    real, allocatable :: tree_lons(:), tree_lats(:)
    integer :: x, y, i
    integer :: unit
    integer :: timer

    real, allocatable :: tmp2d(:,:)
    real, allocatable :: tmp3d(:,:,:)
    real, allocatable :: tmpij(:)
    
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

    ! initialize mpi types
    call g3dv_mpi_setgrid(grid_nx, grid_ny, grid_nz)
    call grid_mpi_init()
    

    ! TODO, 
    grid_ns = grid_nz * 2! + 1
    grid_var_t = 1
    grid_var_s = grid_nz + grid_var_t
!    grid_var_u = grid_nz + grid_var_s
!    grid_var_v = grid_nz + grid_var_u
!    grid_var_ssh = grid_nz +  grid_var_v



    !------------------------------------------------------------
    ! read in the grid specifications
    !------------------------------------------------------------
    allocate(grid_dpth(grid_nz))
    allocate(grid_local_lat(g3dv_mpi_ijcount))
    allocate(grid_local_lon(g3dv_mpi_ijcount))
    allocate(grid_local_mask(g3dv_mpi_ijcount))
    allocate(grid_local_D(g3dv_mpi_ijcount))
    allocate(grid_local_ssh(g3dv_mpi_ijcount))
    allocate(grid_local_coastdist(g3dv_mpi_ijcount))    
    allocate(grid_local_dens(grid_nz, g3dv_mpi_ijcount))
    allocate(grid_local_diag3d_2(grid_ns, g3dv_mpi_ijcount))
    allocate(tmpij(g3dv_mpi_ijcount))
    if(isroot) then
       allocate(tmp2d(grid_nx, grid_ny))
       allocate(tmp3d(grid_nx, grid_ny, grid_nz))       
       allocate(tree_lons(grid_nx*grid_ny))
       allocate(tree_lats(grid_nx*grid_ny))

       print *, ""
       print *, "Reading grid definition parameters..."

       print '(A,I0,A,I0,A,I0)', "  Grid size is ",grid_nx," x ",grid_ny," x ",grid_nz
    else
       !TODO, get rid of this
       allocate(tmp2d(1,1))
       allocate(tmp3d(1,1,1))
    end if
    
    
    ! depth
    if(isroot) then
       call datatable_get('grid_z', grid_dpth)
       call datatable_get('grid_D', tmp2d)
    end if
    call g3dv_mpi_grd2ij_real(tmp2d, grid_local_D)
    call mpi_bcast(grid_dpth, size(grid_dpth), mpi_real, g3dv_mpi_root, g3dv_mpi_comm, i)

    
    ! latitude
    if(isroot) then
       call datatable_get('grid_y', tmp2d)       
       do x = 1, grid_nx
          do y = 1, grid_ny
             tree_lats((y-1)*grid_nx + x) = tmp2d(x,y)
          end do
       end do
    end if
    call g3dv_mpi_grd2ij_real(tmp2d, grid_local_lat)

    
    ! longitude
    if(isroot) then
       call datatable_get('grid_x', tmp2d)
       do x = 1, grid_nx
          do y = 1, grid_ny
             tree_lons((y-1)*grid_nx + x) = tmp2d(x,y)
          end do
       end do
    end if
    call g3dv_mpi_grd2ij_real(tmp2d, grid_local_lon)

    
    ! mask
    if(isroot) then
       call datatable_get('grid_mask', tmp2d)
       print '(A,I0,A,I0,A,F5.1,A)', "  ocean points: ", nint(sum(tmp2d))," of ", &
            grid_nx*grid_ny, " (",sum(tmp2d)/grid_nx/grid_ny*100,"%)"
    end if
    call g3dv_mpi_grd2ij_real(tmp2d, grid_local_mask)

    !coast distance
    if(isroot) then
       call datatable_get('grid_coast', tmp2d)
    end if
    call g3dv_mpi_grd2ij_real(tmp2d, grid_local_coastdist)
       
    ! generate kd-tree for fast lookup of grid points given a lat/lon
    !------------------------------
    if(isroot) then
       print *,""
       print *, "Generating kd-tree for fast lookup of grid x/y given lat/lon..."
       allocate(ll_kdtree_x(grid_nx*grid_ny))
       allocate(ll_kdtree_y(grid_nx*grid_ny))
       do x = 1, grid_nx
          do y = 1, grid_ny
          ll_kdtree_x((y-1)*grid_nx + x) = x
          ll_kdtree_y((y-1)*grid_nx + x) = y
          end do
       end do
       call kd_init(ll_kdtree, tree_lons, tree_lats)
       deallocate(tree_lons)
       deallocate(tree_lats)
    end if

    
    ! read in state fields
    ! TODO, scatter density and MLD instantly (they are only needed for
    ! the vertical scale calculations)
    !------------------------------
    if(isroot) then
       print *, ""
       print *, "Reading background fields..."
    end if

    
    ! ssh
    if(isroot) then
       call datatable_get('bg_ssh', tmp2d)
    end if
    call g3dv_mpi_grd2ij_real(tmp2d, grid_local_ssh)

    
    ! density
    ! TODO, save on memory by reading in a 2d slice at a time and
    ! scattering immediately?
    if(isroot) then
       call datatable_get('bg_dens',tmp3d)
    end if
    do i = 1, grid_nz
       call g3dv_mpi_grd2ij_real(tmp3d(:,:,i), tmpij)
       grid_local_dens(i,:) = tmpij
    end do


    ! done, cleanup
    !------------------------------------------------------------
    deallocate(tmp2d)
    deallocate(tmp3d)
    deallocate(tmpij)
    call timer_stop(timer)
  end subroutine grid_init




  !================================================================================
  !================================================================================



  pure function grid_genweights(lat, lon, depth) result(w)
    real, intent(in) :: lat, lon, depth
    type(grid_interpweights) :: w


    integer :: z
    real    :: zdist, zdist2
    integer :: rpoints(1), rnum
    real    :: rdist(1)

    ! this currently only generates the closest grid points, not full linear interpolation
    ! TODO, do proper linear interpolation of all nearby ocean points
    
    ! get x/y points
    call kd_search_nnearest(ll_kdtree, lon, lat, 1, rpoints, rdist, rnum, .false.)
    w%x   = 0
    w%y   = 0
    w%w_h = 0
    w%x(1) = ll_kdtree_x(rpoints(1))
    w%y(1) = ll_kdtree_y(rpoints(1))
    w%w_h(1) = 1.0

    ! ! if the closest point is on land, abort
    ! ! TODO, find the closest OCEAN point
    ! if(grid_mask(w%x(1), w%y(1)) < 1) then
    !    w%x(1) = -1
    !    return
    ! end if

    ! get z interp
    ! find closest z level    
    zdist = 1e10
    do z = 1, grid_nz
       zdist2 = abs(depth - grid_dpth(z))
       if(zdist2 > zdist) exit
       zdist = zdist2
    end do
    z = z - 1
    w%z = z
    w%w_v = 0

  end function grid_genweights



  !================================================================================
  !================================================================================



  subroutine grid_mpi_init()
    integer, parameter :: n = 5
    integer :: blocklen(n)
    integer :: type(n)
    integer(kind=mpi_address_kind) :: disp(n), base
    type(grid_interpweights) :: grd
    integer :: ierr, i

    call mpi_get_address(grd%x,     disp(1), ierr)
    call mpi_get_address(grd%y,     disp(2), ierr)
    call mpi_get_address(grd%z,     disp(3), ierr)
    call mpi_get_address(grd%w_h,   disp(4), ierr)
    call mpi_get_address(grd%w_v,   disp(5), ierr)
    base = disp(1)
    do i =1,n
       disp(i) = disp(i) - base
    end do

    blocklen = 4
    blocklen(3) = 1
    blocklen(5) = 1

    type(1:3) = mpi_integer
    type(4:5)   = mpi_real
    call mpi_type_create_struct(n, blocklen, disp, type, grid_mpi_interpweights, ierr)
    call mpi_type_commit(grid_mpi_interpweights, ierr)
  end subroutine grid_mpi_init



end module  g3dv_grid
