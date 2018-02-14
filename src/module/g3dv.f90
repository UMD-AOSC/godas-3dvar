module godas_3dvar
  use g3dv_mpi
  use g3dv_datatable
  use g3dv_grid
  use g3dv_obs
  use g3dv_obs_dat
  use g3dv_obs_nc 
  use g3dv_solver
  use g3dv_bgcov
  use timing 
  use netcdf

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
    class(obsio),  pointer :: obsio_ptr

    
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
     call datatable_init(isroot, 'datatable.3dvar')
     call grid_init(isroot, nml_filename)
    
    
     ! initialize background error covariance model
     call bgcov_init(isroot, nml_filename)

    
     ! initialize observations
     allocate(obsio_dat :: obsio_ptr)
     call obs_obsio_reg(obsio_ptr)
     allocate(obsio_nc :: obsio_ptr)
     call obs_obsio_reg(obsio_ptr)    
     call obs_init(isroot, nml_filename, bgcov_local_vtloc, bgcov_local_var_t, bgcov_local_var_s)

    
     ! initialize other modules 
     call solver_init(isroot, nml_filename)


    ! all done
    initialized = .true.
    call timer_stop(timer_ini)
  end subroutine g3dv_init



  !================================================================================
  !================================================================================



  subroutine g3dv_run()
    real :: local_ai(grid_ns, g3dv_mpi_ijcount)

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

    call solver_run(local_ai)
    
    ! TODO, make output to file or retreival from  module API configurable
    timer_output = timer_init('(Output)', TIMER_SYNC)
    call timer_start(timer_output)
    call write_output(local_ai)
    call timer_stop(timer_output)

  end subroutine g3dv_run



  !================================================================================
  !================================================================================



  subroutine write_output(grd)
    real, intent(in) :: grd(grid_ns, g3dv_mpi_ijcount)

    integer :: i
    integer :: ncid1, ncid2
    integer :: d_x, d_y, d_z
    integer :: v_x, v_y, v_z
    integer :: v_ai_t, v_ai_s
    integer :: v_maxobcor_t, v_maxobcor_s

    real, allocatable :: tmp2d(:,:)
    real, allocatable :: tmp3d(:,:,:)
    real :: tmpij(g3dv_mpi_ijcount)
    

    ! setup the output file
    if(isroot) then
       print*, 'Saving analysis increments'

       call check(nf90_create("ana_inc.nc", NF90_WRITE, ncid1))       
       call check(nf90_def_dim(ncid1, "grid_x", grid_nx, d_x))
       call check(nf90_def_var(ncid1, "grid_x", nf90_real, (/d_x/), v_x))
       call check(nf90_put_att(ncid1, v_x, "units", "degrees_east"))
       call check(nf90_def_dim(ncid1, "grid_y", grid_ny, d_y))
       call check(nf90_def_var(ncid1, "grid_y", nf90_real, (/d_y/), v_y))
       call check(nf90_put_att(ncid1, v_y, "units", "degrees_north"))
       call check(nf90_def_dim(ncid1, "grid_z", grid_nz, d_z))
       call check(nf90_def_var(ncid1, "grid_z", nf90_real, (/d_z/), v_z))
       call check(nf90_put_att(ncid1, v_z, "units", "meters"))
       call check(nf90_def_var(ncid1, "Temp", nf90_real, (/d_x, d_y, d_z/), v_ai_t))
       call check(nf90_def_var(ncid1, "Salt", nf90_real, (/d_x, d_y, d_z/), v_ai_s))
       call check(nf90_enddef(ncid1))

       ! other optional diagnostics
       call check(nf90_create("ana_diag.nc", NF90_WRITE, ncid2))       
       call check(nf90_def_dim(ncid2, "grid_x", grid_nx, d_x))
       call check(nf90_def_var(ncid2, "grid_x", nf90_real, (/d_x/), v_x))
       call check(nf90_put_att(ncid2, v_x, "units", "degrees_east"))
       call check(nf90_def_dim(ncid2, "grid_y", grid_ny, d_y))
       call check(nf90_def_var(ncid2, "grid_y", nf90_real, (/d_y/), v_y))
       call check(nf90_put_att(ncid2, v_y, "units", "degrees_north"))
       call check(nf90_def_dim(ncid2, "grid_z", grid_nz, d_z))
       call check(nf90_def_var(ncid2, "grid_z", nf90_real, (/d_z/), v_z))
       call check(nf90_put_att(ncid2, v_z, "units", "meters"))
       call check(nf90_def_var(ncid2, "maxobcorr_t", nf90_real, (/d_x, d_y, d_z/), v_maxobcor_t))
       call check(nf90_def_var(ncid2, "maxobcorr_s", nf90_real, (/d_x, d_y, d_z/), v_maxobcor_s))
       call check(nf90_enddef(ncid2))
   
    end if
    if(isroot) then
       allocate(tmp2d(grid_nx, grid_ny))
       allocate(tmp3d(grid_nx, grid_ny, grid_nz))
    else
       allocate(tmp2d(1,1))
       allocate(tmp3d(1,1,1))
    end if
    

    ! gather data from the procs as needed and save to file
    !------------------------------------------------------------
    ! depth
    if(isroot) call check(nf90_put_var(ncid1, v_z, grid_dpth))
    if(isroot) call check(nf90_put_var(ncid2, v_z, grid_dpth))

    ! lat
    call g3dv_mpi_ij2grd_real(grid_local_lat, tmp2d)
    if(isroot) call check(nf90_put_var(ncid1, v_y, maxval(tmp2d, 1, tmp2d < 100)))
    if(isroot) call check(nf90_put_var(ncid2, v_y, maxval(tmp2d, 1, tmp2d < 100)))

    ! lon
    call g3dv_mpi_ij2grd_real(grid_local_lon, tmp2d)
    if(isroot) call check(nf90_put_var(ncid1, v_x, tmp2d(:,1)))
    if(isroot) call check(nf90_put_var(ncid2, v_x, tmp2d(:,1)))


    ! temp AI
    do i = 1, grid_nz
       tmpij = grd(i+grid_var_t-1,:)
       call g3dv_mpi_ij2grd_real(tmpij, tmp3d(:,:,i))
    end do
    if(isroot) call check(nf90_put_var(ncid1, v_ai_t, tmp3d))

    ! salt AI
    do i = 1, grid_nz
       tmpij = grd(i+grid_var_s-1,:)              
       call g3dv_mpi_ij2grd_real(tmpij, tmp3d(:,:,i))
    end do
    if(isroot) call check(nf90_put_var(ncid1, v_ai_s, tmp3d))

    ! other optional diagnostics
    do i = 1, grid_nz
       tmpij = solver_local_maxobcor(grid_var_t+i-1,:)
       call g3dv_mpi_ij2grd_real(tmpij, tmp3d(:,:,i))
    end do
    if(isroot) call check(nf90_put_var(ncid2, v_maxobcor_t, tmp3d))
    do i = 1, grid_nz
       tmpij = solver_local_maxobcor(grid_var_s+i-1,:)
       call g3dv_mpi_ij2grd_real(tmpij, tmp3d(:,:,i))
    end do
    if(isroot) call check(nf90_put_var(ncid2, v_maxobcor_s, tmp3d))


    ! all done, cleanup
    if(isroot) then
       call check(nf90_close(ncid1))
       call check(nf90_close(ncid2))

    end if
    deallocate(tmp2d)
    deallocate(tmp3d)
  end subroutine write_output


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



  !================================================================================
  !================================================================================
  subroutine check(status)
    integer, intent(in) :: status
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop 1
    end if
  end subroutine check

end module godas_3dvar
