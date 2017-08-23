module g3dv_bgcov
  !! author: Travis Sluka
  !!
  !! g3dv_bgcov description
  !!
  !! g3dv_bgcov detailed description
  !!

  use g3dv_mpi
  use g3dv_obs
  use g3dv_grid
  use g3dv_datatable, only : datatable_get
  use timing

  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: bgcov_init
  public :: bgcov_HBH
  public :: bgcov_BH
  public :: bgcov_hzdist

  
  ! public module variables
  !------------------------------------------------------------
  real, public, allocatable :: bgcov_local_vtloc(:,:)
  real, public, allocatable :: bgcov_local_var_t(:,:)
  real, public, allocatable :: bgcov_local_var_s(:,:)



  ! private module variables
  !------------------------------------------------------------  
  real, parameter :: pi = 4*atan(1.0)
  real, parameter :: re = 6371d3
  real, parameter :: omega = 7.29e-5
  
  ! variables read in from namelist
  real :: hz_loc(2) = -1
  real :: vt_loc_max = -1
  real :: vt_loc_min = -1
  real :: vt_loc_diff_scale = -1
  real :: time_loc = -1
  real :: tnsr_surf  = -1.0   ! surface (SSH) gradient tensor
  real :: tnsr_coast_dist = -1.0       ! coast gradient tensor
  real :: tnsr_coast_min  = 0.0
  real :: bgvar_t = -1
  real :: bgvar_s = -1
  real :: hz_loc_scale = 2.0

  
  interface loc
     module procedure loc_gc
  end interface


contains


  !================================================================================
  !================================================================================



  subroutine bgcov_init(isroot, nml)
    logical, intent(in) :: isroot
    character(len=*), intent(in) :: nml

    integer :: unit, i
    real :: r

    integer :: timer
    real, allocatable :: tmp3d(:,:,:)
    real, allocatable :: tmpij(:)

    namelist /g3dv_hzloc/ hz_loc, hz_loc_scale

    namelist /g3dv_bgcov/ vt_loc_min, vt_loc_max,&
         vt_loc_diff_scale, time_loc, tnsr_surf, tnsr_coast_dist, tnsr_coast_min, bgvar_t, bgvar_s


    if(isroot) then
       print *, new_line('a'), &
            new_line('a'), "------------------------------------------------------------",&
            new_line('a'), "g3dv_bgcov_init() : background error covariance model",&
            new_line('a'), "------------------------------------------------------------"
    end if

    !read in our section of the namelist
    open(newunit=unit, file=nml)
    read(unit, nml=g3dv_bgcov)
    rewind(unit)
    read(unit, nml=g3dv_hzloc)
    close(unit)
    if(isroot) print g3dv_hzloc
    if(isroot) print g3dv_bgcov

    
    ! initialize some variables
    allocate(tmpij(g3dv_mpi_ijcount))

    
    ! determine horizontal correlation lengths
    if(isroot) then
       print *, ""
       print *, " Horizontal correlation length scales summary:"
       print *, " Latitude         scale (km)"
       r = 0
       do while(r <= 90)
          print *, r, bgcov_hzdist(r)
          r = r + 10.0
       end do
    end if



    ! load vertical localization distance field
    ! ------------------------------------------------------------
    timer = timer_init(" bgcov_vtloc", TIMER_SYNC)
    call timer_start(timer)


    if(isroot) then
       print *,""
       print *, "Loading vertical localization values..."
    end if
    allocate(bgcov_local_vtloc(grid_nz, g3dv_mpi_ijcount))
    bgcov_local_vtloc = 0.0
    if(isroot) then
       allocate(tmp3d(grid_nx, grid_ny, grid_nz))       
       call datatable_get('vtloc', tmp3d)
    end if
    do i=1,grid_nz
       call g3dv_mpi_grd2ij_real(tmp3d(:,:,i), tmpij)
       bgcov_local_vtloc(i,:) = tmpij
    end do
    
    
    ! Load the background error variance
    !------------------------------------------------------------
    if (isroot) then
       print *, ""
       print *, "Loading background error variance..."
    end if
    allocate(bgcov_local_var_t(grid_nz, g3dv_mpi_ijcount))
    allocate(bgcov_local_var_s(grid_nz, g3dv_mpi_ijcount))

    ! temperature
    if(bgvar_t >= 0) then
       bgcov_local_var_t = bgvar_t
       if(isroot) print *, "Using fixed TEMP bg err var: ", bgvar_t
    else
       if(isroot)  call datatable_get('bgvar_t', tmp3d)
       do i = 1, grid_nz
          call g3dv_mpi_grd2ij_real(tmp3d(:,:,i), tmpij)
          bgcov_local_var_t(i, :) = tmpij
       end do
    end if

    ! salinity
    if(bgvar_s >= 0) then
       bgcov_local_var_s = bgvar_s
       if(isroot) print *, "Using fixed SALT bg err var: ", bgvar_s
    else
       if(isroot) call datatable_get('bgvar_s', tmp3d)
       do i = 1, grid_nz
          call g3dv_mpi_grd2ij_real(tmp3d(:,:,i), tmpij)
          bgcov_local_var_s(i,:) = tmpij
       end do
    end if

    if(allocated(tmp3d)) deallocate(tmp3d)
    deallocate(tmpij)
   


    ! TODO, density is no longer needed, have the grid module delete that data
    
    !------------------------------------------------------------
    call timer_stop(timer)
    
  end subroutine bgcov_init


  
  !================================================================================
  !================================================================================



  pure function bgcov_HBH(ob1, ob2, dist) result(cov)
    type(observation), intent(in) :: ob1
    type(observation), intent(in) :: ob2
    real, optional, intent(in) :: dist
    real :: cov

    real :: dist0, r, vtloc_d
    real :: hz_cor, vt_cor, time_cor, surf_tensor, coast_tensor

    cov = 0.0
    
    ! univariate for now
    if (ob1%id /= ob2%id) return


    ! gradient tensor localizations
    surf_tensor   = merge( &
                loc(abs(ob1%grd_ssh - ob2%grd_ssh), tnsr_surf),&
                1.0, tnsr_surf > 0.0)
    if(surf_tensor <= 0) return
   

    !calculate distance if it hasn't been given
    ! TODO, move this outside the function to prevent
    ! accidentally calling without a distance, which is slower.
    if(present(dist)) then
       dist0 = dist
    else
       r = sin(ob1%lat*pi/180.0)*sin(ob2%lat*pi/180.0) + &
           cos(ob1%lat*pi/180.0)*cos(ob2%lat*pi/180.0)*&
           cos((ob1%lon-ob2%lon)*pi/180.0)
       if(r > 1)  r =  1
       if(r < -1) r = -1
       dist0 = re*acos(r)
    end if

    ! horizontal localization
    hz_cor = loc(dist0, (bgcov_hzdist(ob1%lat)+ bgcov_hzdist(ob2%lat))/2.0)
    if(hz_cor <= 0) return


    ! vertical localization distance is the average of vt_loc of the two points 
    vtloc_d = sqrt(ob1%grd_vtloc*ob2%grd_vtloc)
!    vtloc_d = (ob1%grd_vtloc+ob2%grd_vtloc)/2.0
    vt_cor = loc(abs(ob1%dpth-ob2%dpth),  vtloc_d)


    !TODO: remove this part? it was probably only needed before vtloc smoothing was done.
    ! vtloc changing too much over a short hz/vt distance can cause the covariance matrix to
    ! not be positive definite, and cholesky decomp will fail. This modulates the vt_cor
    ! by a function of the difference of the two vt_locs.
    if(vt_loc_diff_scale > 0) then
       vt_cor = vt_cor * loc(abs(ob1%grd_vtloc-ob2%grd_vtloc), vtloc_d*vt_loc_diff_scale)
    end if
    if(vt_cor <= 0) return


    ! temporal localization
    time_cor = merge( &
                loc(abs(ob1%hr - ob2%hr), time_loc),&
                1.0, time_loc > 0.0)


    ! coast distance tensor
    coast_tensor = merge( &
        max(tnsr_coast_min, &
           1.0 - abs( min(ob1%grd_coast, tnsr_coast_dist) - &
                      min(ob2%grd_coast, tnsr_coast_dist))&
           / tnsr_coast_dist), &
        1.0, tnsr_coast_dist > 0.0)

    ! add it all up to get the final covariance
    cov = hz_cor * vt_cor * time_cor * surf_tensor * coast_tensor * ob1%grd_var * ob2%grd_var

  end function bgcov_HBH



  !================================================================================
  !================================================================================


  pure subroutine bgcov_BH(ob, dist, ij, cor, var)
    ! return correlation and variance separately
    ! final covariance will just be cor*var

    type(observation), intent(in) :: ob    
    real,    intent(in)  :: dist
    integer, intent(in)  :: ij
    real,    intent(out) :: cor(grid_ns), var(grid_ns)

    real :: hz_cor, time_cor, surf_tensor, coast_tensor, vt_cor(grid_nz)
    real :: vtloc_d

    integer :: i, idx_start
    real :: r


    cor = 0.0
    var = 0.0


    ! gradient tensor localizations
    surf_tensor = merge( &
               loc(abs(ob%grd_ssh - grid_local_ssh(ij)), tnsr_surf),&
               1.0, tnsr_surf > 0.0)
    if(surf_tensor <= 0.0) return


    ! horizontal localization    
    hz_cor = loc(dist, (bgcov_hzdist(grid_local_lat(ij))+ bgcov_hzdist(ob%lat))/2.0)
    if (hz_cor <= 0.0) return


    ! vertical localization distance is the average of vt_loc of the two points 
    vt_cor = 0.0    
    do i = 1, grid_nz
       if(bgcov_local_vtloc(i,ij) <= 0.0) exit
       vtloc_d=sqrt(ob%grd_vtloc*bgcov_local_vtloc(i,ij))
!       vtloc_d=(ob%grd_vtloc+bgcov_local_vtloc(i,ij))/2.0
       r = loc(abs(grid_dpth(i)-ob%dpth),  vtloc_d)

       !TODO: remove this part? it was probably only needed before vtloc smoothing was done.
       ! vtloc changing too much over a short hz/vt distance can cause the covariance matrix to
       ! not be positive definite, and cholesky decomp will fail. This modulates the vt_cor
       ! by a function of the difference of the two vt_locs.
       if(vt_loc_diff_scale > 0) then
          r = r * loc(abs(ob%grd_vtloc-bgcov_local_vtloc(i,ij)), vtloc_d*vt_loc_diff_scale)
       end if      
       vt_cor(i) = r
    end do
    

    ! temporal localization
    time_cor = merge( &
               loc(abs(ob%hr), time_loc),&
               1.0, time_loc > 0.0)


    ! coast distance gradient tensor
    coast_tensor = merge( &
        max(tnsr_coast_min, &
           1.0 - abs( min(ob%grd_coast, tnsr_coast_dist) - &
                      min(grid_local_coastdist(ij), tnsr_coast_dist))&
           / tnsr_coast_dist), &
        1.0, tnsr_coast_dist > 0.0)


    !variance 
    var(grid_var_t:grid_var_t+grid_nz-1) = bgcov_local_var_t(:,ij)*ob%grd_var
    var(grid_var_s:grid_var_s+grid_nz-1) = bgcov_local_var_s(:,ij)*ob%grd_var


    ! univariate correlation for now
    if(ob%id == obs_id_t) then       
       idx_start = grid_var_t
    else if(ob%id == obs_id_s) then
       idx_start = grid_var_s
    else
       return
    end if
    cor(idx_start:idx_start+grid_nz-1) = 1.0

    ! add up all the correlations
    cor(grid_var_t:grid_var_t+grid_nz-1) = cor(grid_var_t:grid_var_t+grid_nz-1) * vt_cor
    cor(grid_var_s:grid_var_s+grid_nz-1) = cor(grid_var_s:grid_var_s+grid_nz-1) * vt_cor
    cor = cor * hz_cor * time_cor * surf_tensor * coast_tensor

  end subroutine bgcov_BH



  !================================================================================
  !================================================================================



  pure function bgcov_hzdist(lat) result(cor)
    real, intent(in) :: lat
    real :: cor
    ! linear interpolation from EQ to Pole
    !    cor = hz_loc(2) + (hz_loc(1)-hz_loc(2))*(1.0-(abs(lat)/90.0))

    ! more accurate rossby radius calculation
    if ( abs(lat) < 0.1) then
       cor = hz_loc(1)
    else
       cor = max(hz_loc(2), min(hz_loc(1), &
            hz_loc_scale*2.6/(2*omega*abs(sin(lat*pi/180.0))) ))
    end if
  end function bgcov_hzdist



  !================================================================================
  !================================================================================



  pure function loc_gaus(z, L)
    !! gaussian localization function
    
    real, intent(in) :: z
    real, intent(in) :: L
    real :: loc_gaus
    loc_gaus = exp( -0.5 * z*z / (L*L))
  end function loc_gaus



  !================================================================================
  !================================================================================



  pure function loc_gc(z, L)
    real, intent(in) :: z
    real, intent(in) :: L
    real :: loc_gc

    real :: c
    real :: abs_z, z_c

    c = L / sqrt(0.3)
    abs_z = abs(z)
    z_c = abs_z / c

    if (abs_z >= 2*c) then
       loc_gc = 0.0
    elseif (abs_z < 2*c .and. abs_z > c) then
       loc_gc = &
            (1.0/12.0)*z_c**5 - 0.5*z_c**4 + &
            (5.0/8.0)*z_c**3 + (5.0/3.0)*z_c**2 &
            - 5.0*z_c + 4 - (2.0/3.0)*c/abs_z
    else
       loc_gc = &
            -0.25*z_c**5 + 0.5*z_c**4 + &
            (5.0/8.0)*z_c**3 - (5.0/3.0)*z_c**2 + 1
    end if
  end function loc_gc


  !================================================================================
  !================================================================================

end module  g3dv_bgcov
