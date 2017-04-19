module godas_3dvar
  use g3dv_mpi
  use g3dv_datatable
  use g3dv_grid
  use g3dv_obs
  use g3dv_solver
  use g3dv_bgcov
  use timing 


  implicit none
  private

  ! public module methods
  !------------------------------------------------------------
  public :: g3dv_init
  public :: g3dv_run
  public :: g3dv_final


  ! private module variables
  !------------------------------------------------------------
  logical :: isroot
  logical :: initialized = .false.
  integer :: timer_total



contains



  !================================================================================
  !================================================================================



  subroutine g3dv_init()
    !! Initialize the 3DVar module, this needs to be called before anything else

    character(len=1024) :: nml_filename="namelist.3dvar"
    logical :: ex
    integer :: timer_ini

    ! initialize some of the timers we are using
    timer_total = timer_init('TOTAL')
    timer_ini   = timer_init('(init)')
    call timer_start(timer_total)
    call timer_start(timer_ini)

    !TODO, allow for passing of mpi_comm from outside
    ! initialize mpi library
    call g3dv_mpi_init()
    call timing_init(g3dv_mpi_comm, g3dv_mpi_root)
    isroot = g3dv_mpi_isroot

    ! print header
    if(isroot) then
       print *,new_line('a'),"============================================================",&
            new_line('a'),"HYBRID-GODAS 3DVAR system",&
            new_line('a')," preconditioned conjugate gradient solver in obs space",&
            new_line('a')," Travis Sluka (travis.sluka@noaa.gov, tsluka@umd.edu)",&
            new_line('a'),"============================================================",&
            new_line('a')
    end if

    ! make sure the namelist is there
    inquire(file=nml_filename, exist=ex)
    if (.not. ex) then
       print*, "ERROR: unable to open namelist: ", trim(nml_filename)
       stop 1
    end if
    if (isroot) then
       print *, 'Using namelist file: ', trim(nml_filename)
    end if

    ! initialize grid / state variables
    call datatable_init(isroot, 'data_table.3dvar')
    call grid_init(isroot, nml_filename)
    call g3dv_mpi_setgrid(grid_nx, grid_ny, grid_nz)
    call grid_scatter()

    ! initizlie observations
    call obs_init(isroot, nml_filename)

    ! initialize other modules 
    call solver_init(isroot, nml_filename)
    call bgcov_init(isroot, nml_filename)


    ! all done
    initialized = .true.
    call timer_stop(timer_ini)
  end subroutine g3dv_init



  !================================================================================
  !================================================================================



  subroutine g3dv_run()
    real, allocatable :: ai(:,:,:)

    integer :: timer_output

    
    !make sure the module has been initialized first
    if( .not. initialized) then
       print *, "ERROR: calling g3dv_run() before g3dv_init()"
       stop 1
    end if

    ! run the preconditioned conjugate gradient solver in obs space
    if(isroot) then
       print *, new_line('a'),&
            new_line('a'), "============================================================",&
            new_line('a'), "Running preconditioned conjugate gradient solver",&
            new_line('a'), "============================================================"
    end if

    call solver_run(ai)

    ! TODO, make output to file or retreival from  module API configurable
    timer_output = timer_init('(Output)', TIMER_SYNC)
    call timer_start(timer_output)
    if(isroot) then
       call grid_write(ai, "ai.nc")
    end if
    call timer_stop(timer_output)

  end subroutine g3dv_run



  !================================================================================
  !================================================================================



  subroutine g3dv_final()
    ! make sure the module has been initialized first
    if( .not. initialized) then
       print *, "ERROR: calling g3dv_final() before g3dv_init()"
       stop 1
    end if

    ! end of run timing
    if (isroot) then
       print *,new_line('a'), &
            new_line('a'), "============================================================",&
            new_line('a'), "End of run statistics",&
            new_line('a'), "============================================================"
    end if
    call timer_stop(timer_total)
    call timer_print()

    ! cleanup MPI
    call g3dv_mpi_final()

  end subroutine g3dv_final


end module godas_3dvar
