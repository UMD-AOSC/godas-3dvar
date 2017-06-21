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
  real :: vt_loc = 0.15
  real :: vt_loc_max = -1
  real :: vt_loc_min = -1
  real :: vt_loc_pow = 1
  real :: vt_loc_diff_scale = 2.0
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

    integer :: unit, i, z
    integer ::  btm_lvl
    real :: r

    integer :: timer
    real, allocatable :: tmp3d(:,:,:)
    real, allocatable :: tmpij(:)

    namelist /g3dv_bgcov/ hz_loc, hz_loc_scale, vt_loc, vt_loc_min, vt_loc_max, vt_loc_pow,&
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
    close(unit)
    if(isroot) print g3dv_bgcov

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


    ! Compute vertical localization distance field
    ! ------------------------------------------------------------
    timer = timer_init(" bgcov_vtloc", TIMER_SYNC)
    call timer_start(timer)

    if(isroot) then
       print *,""
       print *, "Calculating vertical localization values..."
    end if
    allocate(bgcov_local_vtloc(grid_nz, g3dv_mpi_ijcount))
    bgcov_local_vtloc = 0.0
    do i = 1, g3dv_mpi_ijcount
       if(grid_local_mask(i) <= 0.0) cycle

       ! find the bottom level
       do z = 1, grid_nz
          if(grid_dpth(z) > grid_local_D(i)) then
             btm_lvl = z - 1
             exit
          end if
       end do
       ! calculate vertical localization distances
       bgcov_local_vtloc(1:btm_lvl, i) = col_vtloc(&
            grid_local_dens(1:btm_lvl, i), grid_dpth(1:btm_lvl))
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
       allocate(tmpij(g3dv_mpi_ijcount))
       if(isroot) then
          allocate(tmp3d(grid_nx, grid_ny, grid_nz))
          call datatable_get('bgvar_t', tmp3d)
       end if
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
       if(.not. allocated(tmpij)) allocate(tmpij(g3dv_mpi_ijcount))
       if(isroot) then
          if(.not. allocated(tmp3d)) allocate(tmp3d(grid_nx, grid_ny, grid_nz))
          call datatable_get('bgvar_s', tmp3d)
       end if
       do i = 1, grid_nz
          call g3dv_mpi_grd2ij_real(tmp3d(:,:,i), tmpij)
          bgcov_local_var_s(i,:) = tmpij
       end do
    end if

    if(allocated(tmp3d)) then
       deallocate(tmp3d)
       deallocate(tmpij)
    end if
    


    ! TODO, density is no longer needed, have the grid module delete that data
    
    !------------------------------------------------------------
    call timer_stop(timer)
    
  end subroutine bgcov_init


  
  !================================================================================
  !================================================================================


  
  function col_vtloc(dens, dpth) result(vtloc)
    real, intent(in) :: dens(:)
    real, intent(in) :: dpth(:)
    real :: vtloc(size(dens))

    real :: r
    real :: loc_u(size(dens))
    real :: loc_d(size(dens))
    
    integer :: z, z1, z2

    vtloc = 0.0

    !TODO, currently uses linear interpolation, should switch to cubic spline

    
    ! initial pass through the column
    ! to determine localization distances in the up/down directions
    !------------------------------
    loc_u = -1
    loc_d = -1
    do z = 1, size(dens)
       ! calculate UPWARD localization distance
       do z2 = z-1, 1 ,-1          
          if(dens(z)-dens(z2) >= vt_loc) then
             loc_u(z) = dpth(z) - dpth(z2+1) - &
                  (dens(z) - vt_loc-dens(z2+1))&
                  * (dpth(z2) - dpth(z2+1)) &
                  / (dens(z2) - dens(z2+1))
             r = (dpth(z)-dpth(z-1))
             if(loc_u(z) <  r) loc_u(z) = r
             exit
          end if
       end do
       
       ! calculate DOWNWARD localization distance
       do z2 = z+1, size(dens)
          if(dens(z2)-dens(z) >= vt_loc) then
             loc_d(z) = -dpth(z) + dpth(z2-1) + &
                  (dens(z) + vt_loc-dens(z2-1))&
                  * (dpth(z2) - dpth(z2-1)) &
                  / (dens(z2) - dens(z2-1))
             r = (dpth(z+1)-dpth(z))
             if(loc_d(z) <  r) loc_d(z) = r
             exit
          end if
       end do       
    end do

    
  
    ! set the lengths at the top/bottom boundary
    loc_u(1) = loc_d(1)
    loc_d(size(dens)) = loc_u(size(dens))

    
    ! for shallow/stable layers that have no localization lengths set at all:
    if(loc_d(1) < 0 .or. loc_u(size(dens)) < 0) then
       loc_d(1) = dpth(size(dens))*sqrt(0.3)/2.0
       loc_u(1) = loc_d(1)
       loc_u(size(dens)) = dpth(size(dens))*sqrt(0.3)/2.0
       loc_d(size(dens)) = loc_u(size(dens))
!       loc_d = dpth(size(dens))*sqrt(0.3)/2.0
!       loc_u = dpth(size(dens))*sqrt(0.3)/2.0
    end if

    
    ! clip to some min/max value
    if(vt_loc_max > 0) then
       loc_u = min(loc_u, vt_loc_max)
       loc_d = min(loc_d, vt_loc_max)
    end if
    if(vt_loc_min > 0) then
       where(loc_u > 0) loc_u = max(loc_u, vt_loc_min)
       where(loc_d > 0) loc_d = max(loc_d, vt_loc_min)
    end if
    

    
    ! fill in gaps
    !------------------------------

    ! upward lengths
    z1 = -1
    do z2 = 2, size(dens)
       if(loc_u(z2) > 0) then
          z1 = z2
          exit
       end if
    end do
    if(z1 > 2) then
       do z=2,z1-1
          loc_u(z) = loc_u(1) + &
               (dpth(z) - dpth(1)) &
               * (loc_u(z1) - loc_u(1))&
               / (dpth(z1) - dpth(1))
       end do
    end if

    
    ! downward lengths
    z1 = size(dens)+1
    do z2 = size(dens)-1, 1, -1
       if(loc_d(z2) > 0) then
          z1 = z2
          exit
       end if
    end do
    if(z1 < size(dens) - 1) then
       do z = z1+1, size(dens)-1
          loc_d(z) = loc_d(size(dens)) + &
               (dpth(z) - dpth(size(dens))) &
               * (loc_d(z1) - loc_d(size(dens))) &
               / (dpth(z1) - dpth(size(dens)))
       end do
    end if

    ! any small gaps remaining are due to density decreasing slightly with depth
    ! just fill these in with the previous level's value
    do z = 2, size(dens)       
       if(loc_u(z) <= 0) loc_u(z) = loc_u(z-1)
       if(loc_d(z) <= 0) loc_d(z) = loc_d(z-1)
    end do
    
    ! calculate the final correlation lengths
    vtloc = (loc_u * loc_d * min(loc_u,loc_d)**vt_loc_pow)**(1/(2.0+vt_loc_pow))
    
  end function col_vtloc

  
  
  !================================================================================
  !================================================================================



  pure function bgcov_HBH(ob1, ob2, dist) result(cov)
    type(observation), intent(in) :: ob1
    type(observation), intent(in) :: ob2
    real, optional, intent(in) :: dist
    real :: cov

    real :: dist0, r
    real :: hz_cor, vt_cor, time_cor, surf_tensor, coast_tensor

    cov = 0.0
    
    ! univariate for now
    if (ob1%id /= ob2%id) return

    
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


    ! vertical localization distance is the average of vt_loc of the two points 
    ! BUT, this is also modulated by a function of the difference of the two vt_locs.
    ! It makes the algorithm stable... trust me
    vt_cor = loc(abs(ob1%dpth-ob2%dpth),  (ob1%grd_vtloc+ob2%grd_vtloc)/2.0)
    vt_cor = vt_cor *&
         loc(abs(ob1%grd_vtloc-ob2%grd_vtloc), (ob1%grd_vtloc+ob2%grd_vtloc)/(2.0*vt_loc_diff_scale))
    if(vt_cor <= 0) return
    
    ! horizontal localization
    hz_cor = loc(dist0, (bgcov_hzdist(ob1%lat)+ bgcov_hzdist(ob2%lat))/2.0)

    ! temporal localization
    time_cor = merge( &
                loc(abs(ob1%hr - ob2%hr), time_loc),&
                1.0, time_loc > 0.0)

    ! gradient tensor localizations
    surf_tensor   = merge( &
                loc(abs(ob1%grd_ssh - ob2%grd_ssh), tnsr_surf),&
                1.0, tnsr_surf > 0.0)

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

    integer :: i, idx_start
    real :: r


    ! univariate correlation for now
    if(ob%id == obs_id_t) then       
       idx_start = grid_var_t
    else if(ob%id == obs_id_s) then
       idx_start = grid_var_s
    else
       cor =0
       var = 0
       return
    end if
    cor = 0.0
    cor(idx_start:idx_start+grid_nz-1) = 1.0


    ! vertical localization distance is the average of vt_loc of the two points 
    ! BUT, this is also modulated by a function of the difference of the two vt_locs.
    ! It makes the algorithm stable... trust me
    vt_cor = 0.0    
    do i = 1, grid_nz
       if(bgcov_local_vtloc(i,ij) <= 0.0) exit
       r = loc(abs(grid_dpth(i)-ob%dpth),  (ob%grd_vtloc+bgcov_local_vtloc(i,ij))/2.0)
       r = r * loc(abs(ob%grd_vtloc-bgcov_local_vtloc(i,ij)), &
            (ob%grd_vtloc+bgcov_local_vtloc(i,ij))/(2.0*vt_loc_diff_scale))
       vt_cor(i) = r
    end do
    
     
    !variance 
    var(grid_var_t:grid_var_t+grid_nz-1) = bgcov_local_var_t(:,ij)*ob%grd_var
    var(grid_var_s:grid_var_s+grid_nz-1) = bgcov_local_var_s(:,ij)*ob%grd_var


    ! horizontal localization    
    hz_cor = loc(dist, (bgcov_hzdist(grid_local_lat(ij))+ bgcov_hzdist(ob%lat))/2.0)


    ! temporal localization
    time_cor = merge( &
               loc(abs(ob%hr), time_loc),&
               1.0, time_loc > 0.0)


    ! gradient tensor localizations
    surf_tensor = merge( &
               loc(abs(ob%grd_ssh - grid_local_ssh(ij)), tnsr_surf),&
               1.0, tnsr_surf > 0.0)


    coast_tensor = merge( &
        max(tnsr_coast_min, &
           1.0 - abs( min(ob%grd_coast, tnsr_coast_dist) - &
                      min(grid_local_coastdist(ij), tnsr_coast_dist))&
           / tnsr_coast_dist), &
        1.0, tnsr_coast_dist > 0.0)


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
