module g3dv_obs
  !! author: Travis Sluka
  !!
  !! g3dv_obs description
  !!
  !! g3dv_obs detailed description
  !!

  use mpi

  use timing
  use kdtree
  use qsort

  use g3dv_grid
  use g3dv_mpi

#ifdef __INTEL_COMPILER
  use ieee_arithmetic
#endif

  implicit none
  private


  
  !------------------------------------------------------------
  ! public module methods
  !------------------------------------------------------------
  public :: obs_init
  public :: obs_search
  public :: obs_obsio_reg


  
  !------------------------------------------------------------
  ! public module variables
  !------------------------------------------------------------
  type, public :: observation
     ! observation type / location
     integer :: id = -1
     real :: lat
     real :: lon
     real :: dpth
     real :: hr

     ! observation increment and error
     real :: inc
     real :: err

     ! grid parameters interpolated to obs location     
     type(grid_interpweights) :: grd_w     
     real    :: grd_ssh
     real    :: grd_var
     real    :: grd_vtloc
     real    :: grd_coast
  end type observation

  
  integer, public, protected :: obs_num = 0
  !! number of valid observations in **obs**

  type(observation), public, target, allocatable :: obs(:) ! TODO, add protected back
  !! array of all observations that are to be assimilated
  !! **obs[1:obs_num]** have passed QC and are to be used
  !! (array is bigger than obs_num if some obs have failed QC)

  integer, public, protected, allocatable :: obs_local_idx(:)
  !! list of indices pointing to obs that this proc is reponsible
  !! for handling in the solver code


  
  !--------------------
  ! observation blocks
  integer, public, protected, allocatable :: obs_block_size(:)
  integer, public, protected, allocatable :: obs_block_start(:)
  integer, public, protected, allocatable :: obs_block_end(:)
  integer, public, protected, allocatable :: obs_block_proc(:)

  !--------------------
  ! observation ids  
  integer, public, protected :: obs_id_t
  integer, public, protected :: obs_id_s
  integer, public, protected :: obs_id_u
  integer, public, protected :: obs_id_v
  integer, public, protected :: obs_id_ssh


  
  ! obsio interface
  type, abstract, public :: obsio
     integer :: i
   contains
     procedure(I_obsio_getstr), deferred :: get_name
     procedure(I_obsio_getstr), deferred :: get_desc
     procedure(I_obsio_write),  deferred :: write
     procedure(I_obsio_read),   deferred :: read
  end type obsio
  abstract interface
     function I_obsio_getstr(self)
       import obsio
       class(obsio) :: self
       character(:), allocatable :: I_obsio_getstr
     end function I_obsio_getstr
     subroutine I_obsio_init(self)
       import obsio
       class(obsio) :: self
     end subroutine I_obsio_init
     subroutine I_obsio_write(self, file, obs)
       import obsio
       import observation
       class(obsio) :: self
       character(len=*),  intent(in) :: file
       type(observation), intent(in) :: obs(:)
     end subroutine I_obsio_write
     subroutine I_obsio_read(self, file, obs)
       import obsio
       import observation
       class(obsio) :: self
       character(len=*), intent(in) :: file
       type(observation), allocatable, intent(out) :: obs(:)
     end subroutine I_obsio_read
  end interface

  
  !------------------------------------------------------------
  ! private module variables
  !------------------------------------------------------------
  logical :: isroot
  type(kd_root)    :: obs_kdtree

  integer :: obs_mpi_observation

  integer :: obs_block_nx, obs_block_ny
  integer :: obs_block_dx=1, obs_block_dy=1
  
  ! obsio class registration and selection
  type obsioptr
     class(obsio), pointer ::p
  end type obsioptr 
  integer, parameter    :: obsio_reg_max = 100
  integer               :: obsio_reg_num = 0
  type(obsioptr)        :: obsio_reg(obsio_reg_max)
  class(obsio), pointer :: obsio_class


  
contains



  !================================================================================
  !================================================================================



  subroutine obs_init(root, nml, local_vtloc, local_var)
    !! TODO, describe what this does... this does a lot
    
    logical, intent(in) :: root
    character(len=*), intent(in) :: nml
    real, intent(in) :: local_vtloc(:,:)
    real, intent(in) :: local_var(:,:)

    character(len=:), allocatable :: ioclass
    
    integer :: timer
    integer :: unit, i, j, i1

    ! parameters read in from namelist
    logical :: test_obs = .false.
    integer :: test_obs_max = 20
    character(len=:), allocatable :: obsfile 
    real :: depth_max = -1
    logical :: use_t = .true.
    logical :: use_s = .true.

    real :: inc_max_t, inc_max_s
    real, allocatable :: kd_lons(:), kd_lats(:)

    real :: stats_s(2), stats_t(2)
    real :: stats_err_t(2), stats_err_s(2)

    real, allocatable :: tmp2d(:,:)
    real, allocatable :: tmp3d(:,:,:)
    real, allocatable :: tmpij(:)
    
    namelist /g3dv_obs/ test_obs, test_obs_max, ioclass, obsfile, depth_max,&
         obs_id_t, obs_id_s, use_t, use_s, inc_max_t, inc_max_s, obs_block_dx, obs_block_dy


    isroot = root
    timer = timer_init(" observations", TIMER_SYNC)
    call timer_start(timer)

    if(isroot) then
       print *, new_line('a'),&
            new_line('a'), "------------------------------------------------------------",&
            new_line('a'), "obs_init() : ",&
            new_line('a'), "------------------------------------------------------------"
    end if


    ! read in our section of the namelist
    allocate(character(1024) :: ioclass)
    allocate(character(1024) :: obsfile)
    open(newunit=unit, file=nml)
    read(unit, nml=g3dv_obs)    
    close(unit)
    ioclass = trim(ioclass)
    obsfile = trim(obsfile)
    if(isroot) then 
       print g3dv_obs
       print *, ""
    end if
    obs_block_nx = ceiling(1.0*grid_nx / obs_block_dx)
    obs_block_ny = ceiling(1.0*grid_ny / obs_block_dy)


    ! initialize MPI structure
    call obs_mpi_init()

    
    ! print a list of all the obsio classes that have been registered
    if(isroot) then
       print *, ""
       print *, "List of obsio classes registered:"
       do i=1,obsio_reg_num
          print "(A,A,3A)", " * ", trim(obsio_reg(i)%p%get_name()), &
               " (",obsio_reg(i)%p%get_desc(), ")"
       end do
       print *, ""
    end if

    
    !select the obsio class
    if(.not. test_obs) then
       nullify(obsio_class)
       do i = 1, obsio_reg_num
          if(trim(obsio_reg(i)%p%get_name()) == trim(ioClass)) then
             obsio_class => obsio_reg(i)%p
             exit
          end if
       end do
       if(.not. associated(obsio_class)) then
          if(isroot) &
               print*, 'ERROR: obsio class "', trim(ioclass),&
               '"not found. check that hte name is in the list of registered classes'
          stop 1
       end if
    end if


    ! ------------------------------
    allocate(tmpij(g3dv_mpi_ijcount))
    if (isroot) then
       allocate(tmp2d(grid_nx, grid_ny))
       allocate(tmp3d(grid_nx, grid_ny, grid_nz))
    else
       ! fortran complains in debug mode about tmp2d and tmp3d
       ! being used uninitialized (which they aren't), this is just
       ! to keep Fortran from complaining
       allocate(tmp2d(1,1))
       allocate(tmp3d(1,1,1))
    end if


    ! load/generate obs
    ! NOTE: this is performed only on the root proc
    !------------------------------
    if (isroot) then
       ! get the observerations by reading in
       ! or generating from test obs
       if (test_obs) then
          allocate(obs(test_obs_max))
          call genTestObs()
       else
          call obsio_class%read(obsfile, obs)
          obs_num = size(obs)
       end if       
       print *, "observations read in:",obs_num
       
       ! Do QC checks
       ! ------------------------------------------------------------
       print *, ""
       print *, "Performing QC..."
       
       ! remove obs below the max depth
       i1 = obs_num
       do i = obs_num, 1, -1
          if(depth_max >0 .and. obs(i)%dpth > depth_max) then
             obs(i) = obs(obs_num)
             obs_num = obs_num - 1
          end if
       end do
       print *, (i1-obs_num), " removed for BELOW MAX DEPTH"
       i1 = obs_num
       
       ! remove obs types we dont want to assimilate
       do i = obs_num, 1, -1
          if( (.not. use_t .and. obs(i)%id == obs_id_t ) .or. &
              (.not. use_s .and. obs(i)%id == obs_id_s ) ) then
             obs(i) = obs(obs_num)
             obs_num = obs_num - 1
          end if
       end do
       print *, (i1-obs_num), " removed for TYPE NOT ASSIMILATED"
       i1 = obs_num
       
       ! remove obs that fail the gross QC checks
       do i = obs_num, 1, -1
          if( (obs(i)%id == obs_id_t .and. abs(obs(i)%inc) > inc_max_t) .or. &
               (obs(i)%id == obs_id_s .and. abs(obs(i)%inc) > inc_max_s)) then
             obs(i) = obs(obs_num)
             obs_num = obs_num -1
          end if
       end do
       print *, (i1-obs_num), " removed for FAILING GROSS QC"
       i1 = obs_num

       ! remove obs that contain NaNs
       do i = obs_num, 1, -1
#ifdef __INTEL_COMPILER
          if(ieee_is_nan(obs(i)%inc) .or. ieee_is_nan(obs(i)%err)) then
             obs(i) = obs(obs_num)
             obs_num = obs_num -1
          end if
#else
          if(isnan(obs(i)%inc) .or. isnan(obs(i)%err)) then
             obs(i) = obs(obs_num)
             obs_num = obs_num -1
          end if
#endif
       end do
       print *, (i1-obs_num), " removed for NAN FOUND"
       i1 = obs_num
    end if

    ! generate interpolation weights / values    
    call g3dv_mpi_ij2grd_real(grid_local_mask, tmp2d)
    if(isroot) then
       i1 = obs_num       
       do i = obs_num, 1, -1
          ! generate interpolation weights          
          obs(i)%grd_w = grid_genweights(obs(i)%lat, obs(i)%lon, obs(i)%dpth)

          ! remove obs where interpolation weights place it on land
          ! TODO, fix the gen_weights function so this never happens
          if(tmp2d(obs(i)%grd_w%x(1), obs(i)%grd_w%y(1)) < 1) then
             obs(i) = obs(obs_num)
             obs_num = obs_num - 1
             cycle
          end if
       end do
       print *, (i1-obs_num), " removed for ON LAND"
    end if
    

    
    !interpolate certain other fields to the observation locations
    !------------------------------------------------------------

    ! ocean depth at obs loc
    call g3dv_mpi_ij2grd_real(grid_local_D, tmp2d)
    if (isroot) then
       !TODO, this check should really be done in the obsop
       i1 = obs_num
       do i = obs_num, 1, -1          
          if(grid_dpth(obs(i)%grd_w%z) >= tmp2d(obs(i)%grd_w%x(1),obs(i)%grd_w%y(1))) then
             obs(i) = obs(obs_num)
             obs_num = obs_num - 1
          end if
       end do
       print *, (i1-obs_num), " removed for BELOW MODEL DEPTH"
       i1 = obs_num
    end if 


    ! SSH 
    call g3dv_mpi_ij2grd_real(grid_local_ssh, tmp2d)
    if (isroot) then
       do i = 1, obs_num
          obs(i)%grd_ssh = grid_interp2d(tmp2d, obs(i)%grd_w)
       end do
    end if


    ! vertical localization distance
    do i = 1, grid_nz
       tmpij = local_vtloc(i,:)
       call g3dv_mpi_ij2grd_real(tmpij, tmp2d)
       if (isroot) tmp3d(:,:,i) = tmp2d
    end do
    if(isroot) then
       do i = 1, obs_num
          obs(i)%grd_vtloc = grid_interp3d(tmp3d, obs(i)%grd_w)
       end do
    end if

    ! coastline tensor
    call g3dv_mpi_ij2grd_real(grid_local_coastdist, tmp2d)
    if (isroot) then
       do i = 1, obs_num
          obs(i)%grd_coast = grid_interp2d(tmp2d, obs(i)%grd_w)
       end do
    end if 
   

    ! background variance
    ! TODO, do this correctly from the 3D field
    if(isroot) then
       do i = 1, obs_num
          if(obs(i)%id == obs_id_t) obs(i)%grd_var = 1.0
          if(obs(i)%id == obs_id_s) obs(i)%grd_var = 0.1
       end do
    end if


    ! other stats to print
    !------------------------------
    if(isroot) then
       stats_t(1) = 1e10
       stats_s(1) = 1e10
       stats_t(2) = -1e10
       stats_s(2) = -1e10
       stats_err_t(1) = 1e10
       stats_err_s(1) = 1e10
       stats_err_t(2) = -1e10
       stats_err_s(2) = -1e10
       do i = 1, obs_num
          if(obs(i)%id == obs_id_t) then
             if(obs(i)%inc > stats_t(2)) stats_t(2) = obs(i)%inc
             if(obs(i)%inc < stats_t(1)) stats_t(1) = obs(i)%inc
             if(obs(i)%err > stats_err_t(2)) stats_err_t(2) = obs(i)%err
             if(obs(i)%err < stats_err_t(1)) stats_err_t(1) = obs(i)%err
          end if
          if(obs(i)%id == obs_id_s) then
             if(obs(i)%inc > stats_s(2)) stats_s(2) = obs(i)%inc
             if(obs(i)%inc < stats_s(1)) stats_s(1) = obs(i)%inc
             if(obs(i)%err > stats_err_s(2)) stats_err_s(2) = obs(i)%err
             if(obs(i)%err < stats_err_s(1)) stats_err_s(1) = obs(i)%err
          end if          
       end do

       print '(A,I0,A)', "  ",obs_num," obs kepts after QC"       
       print *, "T increment min/max: ",stats_t
       print *, "T error     min/max: ",stats_err_t
       print *, "S increment min/max: ",stats_s            
       print *, "S error     min/max: ",stats_err_s
    end if    

        
    
    ! organize into blocks, obs_block_* variables will be set after this, 
    ! broadcast the blocks
    !------------------------------
    if(isroot) call blockdiv()
    
    if(isroot) i = size(obs_block_size)
    call mpi_bcast(i, 1, mpi_integer, g3dv_mpi_root, g3dv_mpi_comm, j)
    if(.not. isroot) then
       allocate(obs_block_size(i))
       allocate(obs_block_start(i))
       allocate(obs_block_end(i))
       allocate(obs_block_proc(i))
    end if
    call mpi_bcast(obs_block_size,  i, mpi_integer, g3dv_mpi_root, g3dv_mpi_comm, j)
    call mpi_bcast(obs_block_start, i, mpi_integer, g3dv_mpi_root, g3dv_mpi_comm, j)
    call mpi_bcast(obs_block_end,   i, mpi_integer, g3dv_mpi_root, g3dv_mpi_comm, j)
    call mpi_bcast(obs_block_proc,  i, mpi_integer, g3dv_mpi_root, g3dv_mpi_comm, j)

    
    ! Give every proc a copy of the obs
    !------------------------------
    call mpi_bcast(obs_num, 1, mpi_integer, g3dv_mpi_root, g3dv_mpi_comm, i)
    if(.not. isroot) allocate(obs(obs_num))    
    call mpi_bcast(obs, obs_num, obs_mpi_observation, g3dv_mpi_root, g3dv_mpi_comm, i)   

    
    ! determine proc distribution for observations
    ! for now this is just a round-robin
    !------------------------------------------------------------
    i = obs_num / g3dv_mpi_size
    j = mod(obs_num, g3dv_mpi_size)
    if( j > 0 .and. g3dv_mpi_rank < j) i = i + 1
    allocate(obs_local_idx(i))
    do i =1, size(obs_local_idx)
       obs_local_idx(i) = g3dv_mpi_rank + (i-1)*g3dv_mpi_size + 1
    end do
    i = size(obs_local_idx)
    if(isroot) then
       print *, ""
       print *, "Observation mpi distribution:"
       if (j == 0) then
          print '(A,I0,A)', "   each proc responsible for ", i, " obs"
       else
          print '(A,I0,A,I0,A)', "   each proc responsible for ", i-1, " to ",i," obs"
       end if
    end if


    ! create KD tree
    !------------------------------
    allocate(kd_lons(obs_num))
    allocate(kd_lats(obs_num))
    do i=1, obs_num
       kd_lons(i) = obs(i)%lon
       kd_lats(i) = obs(i)%lat
    end do
    call kd_init(obs_kdtree, kd_lons, kd_lats)
    deallocate(kd_lons)
    deallocate(kd_lats)

    
    ! All done, cleanup
    !------------------------------------------------------------
    deallocate(tmp2d)
    deallocate(tmp3d)
    deallocate(tmpij)
    call timer_stop(timer)

  end subroutine



  !================================================================================
  !================================================================================



  subroutine obs_search(lat, lon, dist, r_points, r_dist, r_num)
    real, intent(in) :: lat, lon ,dist
    integer, intent(inout) :: r_points(:)
    real,    intent(inout) :: r_dist(:)
    integer, intent(out)   :: r_num

    call kd_search_radius(obs_kdtree, lon, lat, dist, r_points, r_dist, r_num, .false.)

  end subroutine obs_search



  !================================================================================
  !================================================================================



  subroutine obs_obsio_reg(cls)
    class(obsio), pointer :: cls
    integer :: i

    if(obsio_reg_num == obsio_reg_max) then
       print *, "ERROR: too many obsio classes registered"
       stop 1
    end if

    !ensure there is not a class of the same name already registered
    do i = 1, obsio_reg_num
       if(obsio_reg(i)%p%get_name() == cls%get_name()) then
          print *, "ERROR: can't register obsio class '", trim(cls%get_name()), &
               "', a class by that name has already been registered"
          stop 1
       end if
    end do

    obsio_reg_num = obsio_reg_num + 1
    obsio_reg(obsio_reg_num)%p => cls
  end subroutine obs_obsio_reg

  


  !================================================================================
  !================================================================================

  
  ! subroutine obs_setbgcov(idx, var)
  !   integer, intent(in) :: idx
  !   real, intent(in) :: var

  !   obs(idx)%grd_var = var
  ! end subroutine obs_setbgcov


  
  subroutine genTestObs()
    character(len=1024) :: filename = "test_obs.3dvar"
    logical :: ex
    integer :: unit, iostat, i
    character(len=1024) :: line

    type testobs_line
       integer :: id
       real    :: lat
       real    :: lon
       real    :: dpth
       real    :: hr
       real    :: inc
       real    :: err
    end type testobs_line
    type(testobs_line) ob_tmp



    if(isroot) print *, "Generating TEST OBSERVATION locations..."

    ! read in the config file
    if (isroot) print *, "reading configuration from ",trim(filename),"..."
    inquire(file=filename, exist=ex)
    if (.not. ex) then
       print *, "ERROR: unable to open file ",trim(filename)
       stop 1
    end if
    open(newunit=unit, file=filename, action='read')
    do while(.true.)
       ! read in a new line
       read(unit, '(A)', iostat=iostat) line
       if (iostat <0) exit
       if (iostat > 0) then
          print *, "ERROR: problem reading file. Error code: ", iostat
          stop 1
       end if

       ! convert tabs to spaces
       do i = 1, len(line)
          if (line(i:i) == char(9)) line(i:i) = ' '
       end do

       ! ignore comments and empty lines
       line = adjustl(line)
       if (line(1:1) == '#') cycle
       if (len(trim(adjustl(line))) == 0) cycle

       ! parse the line
       read(line, *, iostat=iostat) ob_tmp
       if (iostat > 0) then
          print *, "ERROR: problem reading data file, line: ", trim(line)
          stop 1
       end if

       obs_num = obs_num + 1
       obs(obs_num)%id   = ob_tmp%id
       obs(obs_num)%lat  = ob_tmp%lat
       obs(obs_num)%lon  = ob_tmp%lon
       obs(obs_num)%dpth = ob_tmp%dpth
       obs(obs_num)%hr   = ob_tmp%hr
       obs(obs_num)%inc  = ob_tmp%inc
       obs(obs_num)%err  = ob_tmp%err
       if (isroot)   print ('(A, I0,6F8.1)'), "  ",ob_tmp
    end do
    close(unit)
    if(isroot) print *, ""

  end subroutine genTestObs



  !================================================================================
  !================================================================================

  subroutine  blockdiv()
    integer :: block_idx(obs_num)
    integer :: obs_idx(obs_num)
    integer :: i, j, k, blocks
    type(observation) :: obs_tmp(obs_num)

    ! variable used for block load balancing
    real :: load(g3dv_mpi_size)
    integer, allocatable :: blockdiv_idx(:)
    integer, allocatable :: blockdiv_load(:)
    
    ! for each observation, determine which block it belongs in
    do i = 1, obs_num
       block_idx(i) = ((obs(i)%grd_w%y(1)-1)/obs_block_dy) * obs_block_nx &
                    + ((obs(i)%grd_w%x(1)-1)/obs_block_dx)+1
    end do

    ! sort the obs list by block number
    do i = 1, obs_num
       obs_idx(i) = i
    end do
    call qsort_2i(block_idx, obs_idx)
    do i = 1, obs_num
       obs_tmp(i) = obs(obs_idx(i))
    end do
    obs(1:obs_num) = obs_tmp

    ! determine total number of blocks
    blocks = 0
    j = 0
    do i = 1, obs_num
       if(block_idx(i) /= j) then
          j = block_idx(i)
          blocks = blocks + 1
       end if
    end do

    ! determine block start/stop/size
    allocate(obs_block_size(blocks))
    allocate(obs_block_start(blocks))
    allocate(obs_block_end(blocks))
    obs_block_size = 0
    k = 0
    j = 0
    do i = 1, obs_num
       if(block_idx(i) /= j) then
          k = k+1
          j = block_idx(i)
          obs_block_start(k) = i
       end if
       obs_block_end(k) = i
       obs_block_size(k) = obs_block_end(k)-obs_block_start(k)+1
    end do 

    
    ! determine block distribution across procs
    !------------------------------
    allocate(obs_block_proc(blocks))
    allocate(blockdiv_idx(blocks))
    allocate(blockdiv_load(blocks))
    ! determine complexity of each block ( O(n^2) )    
    do i = 1, blocks
       blockdiv_idx(i)  = i
       blockdiv_load(i) = obs_block_size(i)**2
    end do
    ! temporarily sort the blocks in order of complexity    
    call qsort_2i(blockdiv_load, blockdiv_idx)
    ! starting with the largest blocks, assign each block
    ! to the proc with currently the lightest load
    load = 0
    do i = blocks, 1, -1
       j = minloc(load,1)
       load(j) = load(j) + blockdiv_load(i)
       obs_block_proc(blockdiv_idx(i)) = j-1
    end do
    !obs_block_proc now lists the proc responsible for each block

    
    ! print out block statistics
    print *, ""
    print *, "Observation blocks (for preconditioner):"
    print *, "  blocks  = ",size(obs_block_start)
    print *, "  max obs = ", maxval(obs_block_size)
    print *, "  min obs = ", minval(obs_block_size)
    print *, "  ave obs = ", sum(obs_block_size)*1.0 / size(obs_block_start)
  end subroutine blockdiv

  !================================================================================
  !================================================================================


  
  subroutine obs_mpi_init()
    !! define a derived type to MPI
    integer, parameter :: n = 12
    integer :: blocklen(n)
    integer :: type(n)
    integer(kind=mpi_address_kind) :: disp(n), base
    type(observation) :: ob
    integer :: ierr, i

    call mpi_get_address(ob%id,        disp( 1),  ierr)
    call mpi_get_address(ob%lat,       disp( 2),  ierr)
    call mpi_get_address(ob%lon,       disp( 3),  ierr)
    call mpi_get_address(ob%dpth,      disp( 4),  ierr)
    call mpi_get_address(ob%hr,        disp( 5),  ierr)
    call mpi_get_address(ob%inc,       disp( 6),  ierr)
    call mpi_get_address(ob%err,       disp( 7),  ierr)
    call mpi_get_address(ob%grd_w,     disp( 8), ierr)
    call mpi_get_address(ob%grd_ssh,   disp( 9), ierr)    
    call mpi_get_address(ob%grd_var,   disp(10), ierr)
    call mpi_get_address(ob%grd_vtloc, disp(11), ierr)
    call mpi_get_address(ob%grd_coast, disp(12), ierr)
    base = disp(1)
    do i=1,n
       disp(i) = disp(i) - base
    end do

    blocklen = 1

    type(1)     = mpi_integer
    type(2:7)   = mpi_real
    type(8)     = grid_mpi_interpweights
    type(9:12)  = mpi_real
    
    call mpi_type_create_struct(n, blocklen, disp, type, obs_mpi_observation, ierr)
    call mpi_type_commit(obs_mpi_observation, ierr)
  end subroutine obs_mpi_init

  
end module g3dv_obs
