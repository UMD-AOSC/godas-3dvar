module g3dv_obs
  !! author: Travis Sluka
  !!
  !! g3dv_obs description
  !!
  !! g3dv_obs detailed description
  !!

  use timing
  use kdtree
  use g3dv_grid
  use g3dv_obs_types

  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: obs_init  
  public :: obs_search


  ! public module variables
  !------------------------------------------------------------
  type(observation), public, protected, target, allocatable :: obs(:)
  !! array of all observations that are to be assimilated
  !! **obs[1:obs_num]** have passed QC and are to be used
  !! (array is bigger than obs_num if some obs have failed QC)

  integer, public, protected :: obs_num = 0
  !! number of valid observations in **obs**

  integer, public, protected, allocatable :: obs_block_size(:)
  integer, public, protected, allocatable :: obs_block_start(:)
  integer, public, protected, allocatable :: obs_block_end(:)
  integer, public, protected, allocatable :: obs_block_proc(:)


  ! private module variables
  !------------------------------------------------------------
  logical :: isroot
  type(kd_root)    :: obs_kdtree



contains



  !================================================================================
  !================================================================================




  subroutine obs_init(root, nml)
    logical, intent(in) :: root
    character(len=*), intent(in) :: nml
    
    integer :: timer
    integer :: unit, i
    
    logical :: test_obs
    integer :: test_obs_max

    real, allocatable :: kd_lons(:), kd_lats(:)

    namelist /g3dv_obs/ test_obs, test_obs_max


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
    !------------------------------
    open(newunit=unit, file=nml)
    read(unit, nml=g3dv_obs)
    close(unit)
    if(isroot) then 
       print g3dv_obs
       print *, ""
    end if


    ! load in the observations
    !------------------------------
    if (test_obs) then
       allocate(obs(test_obs_max))
       call genTestObs(nml)
    else
       print *, "ERROR: obs not being read in yet."
       stop 1
    end if

    if (isroot) print *, "observations read in:",obs_num

    ! process the observations that have been read in
    ! TODO, do QC checks
    ! TODO remove bad obs
    ! TODO determine intepolation weights
    if (isroot) then
       print *, ""
       print *, "Performing QC (TODO, actually DO QC !)..."
    end if


    ! separate into blocks
    !------------------------------
    ! TODO, do a clustering algorithm?
    ! TODO, this is all temporary here
    allocate(obs_block_size(1))
    allocate(obs_block_start(1))
    allocate(obs_block_end(1))
    obs_block_size(1) = obs_num
    obs_block_start(1) = 1
    obs_block_end(1) = obs_num
    if (isroot) then
       print *, ""
       print *, "Observation clustering:"
       print *, "  blocks         = ", size(obs_block_size)
       print *, "  min block size = ", minval(obs_block_size)
       print *, "  max block size = ", maxval(obs_block_size)
    end if

    ! determine load balancing of blocks
    ! ------------------------------
    ! TODO, this is temporary
    allocate(obs_block_proc(1))
    obs_block_proc(1) = 0

    ! create KD tree
    !------------------------------
    allocate(kd_lons(obs_num))
    allocate(kd_lats(obs_num))
    do i=1, obs_num
       kd_lons(i) = obs(i)%lon
       kd_lats(i) = obs(i)%lat
    end do
    if(isroot) then
       print *, ""
       print *, "Creating observation KD-tree..."
    end if
    call kd_init(obs_kdtree, kd_lons, kd_lats)
    deallocate(kd_lons)
    deallocate(kd_lats)

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



  subroutine genTestObs(nml)
    character(len=*), intent(in) :: nml

    character(len=1024) :: filename = "test_obs.3dvar"
    logical :: ex
    integer :: unit, iostat, i
    character(len=1024) :: line

    type testobs_line
       character(len=16) :: id
       real  :: lat
       real  :: lon
       real  :: dpth
       real  :: hr
       real  :: inc
       real  :: err
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
       obs(obs_num)%id   = 1 !TODO, do ID lookup
       obs(obs_num)%lat  = ob_tmp%lat
       obs(obs_num)%lon  = ob_tmp%lon
       obs(obs_num)%dpth = ob_tmp%dpth
       obs(obs_num)%hr   = ob_tmp%hr
       obs(obs_num)%inc  = ob_tmp%inc
       obs(obs_num)%err  = ob_tmp%err
       if (isroot)   print ('(A, A6,6F8.1)'), "  ",ob_tmp
    end do
    close(unit)
    if(isroot) print *, ""

  end subroutine genTestObs

end module  g3dv_obs
