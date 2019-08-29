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
  real, parameter :: d2r = pi / 180.0
  
  ! variables read in from namelist
  real :: vt_loc_diff_scale = -1
  real :: time_loc = -1
  real :: tnsr_surf  = -1.0   ! surface (SSH) gradient tensor
  real :: tnsr_coast_dist = -1.0       ! coast gradient tensor
  real :: tnsr_coast_min  = 0.0

  real :: vtcor_min = 1.0
  real :: vtcor_max = 250.0
  real :: rho_delta = 0.125
  
  real :: hzcor_min = 50e3
  real :: hzcor_max = 150e3
  real :: hzcor_cspd = 2.7
  real :: hzcor_stretch_lat = 10
  real :: hzcor_stretch_mul = 3.5

  real :: t_varmin_len        = -1
  real :: t_varmin_do         = -1
  real :: t_varmin_surf_const = -1
  real :: t_varmax            = -1
  real :: t_dz                = -1

  real :: s_varmin_len        = -1
  real :: s_varmin_do         = -1
  real :: s_varmin_surf_const = -1
  real :: s_varmax            = -1
  real :: s_dz                = -1

  integer :: gauss_iter = 2
  
  interface loc
     module procedure loc_gc
  end interface


contains


  !================================================================================
  !================================================================================



  subroutine bgcov_init(isroot, nml)
    logical, intent(in) :: isroot
    character(len=*), intent(in) :: nml

    integer :: unit, i, j,d
    real :: r

    integer :: timer

    real :: work_xyk(grid_nx, grid_ny, g3dv_mpi_kcount)
    
    logical :: global_mask(grid_nx, grid_ny)
    real :: global_depth(grid_nx, grid_ny)
    real :: global_lat(grid_nx, grid_ny)
    real :: global_lon(grid_nx, grid_ny)
    

    namelist /g3dv_hzcor/ hzcor_min, hzcor_max, hzcor_cspd, hzcor_stretch_lat, hzcor_stretch_mul

    namelist /g3dv_vtcor/ rho_delta, vtcor_min, vtcor_max
    
    namelist /g3dv_bgcov/ gauss_iter, vt_loc_diff_scale, time_loc,&
         tnsr_surf, tnsr_coast_dist, tnsr_coast_min, &
         t_varmin_len, t_varmin_do, t_varmin_surf_const, t_varmax, t_dz, &
         s_varmin_len, s_varmin_do, s_varmin_surf_const, s_varmax, s_dz

    if(isroot) then
       print *, new_line('a'), &
            new_line('a'), "------------------------------------------------------------",&
            new_line('a'), "g3dv_bgcov_init() : background error covariance model",&
            new_line('a'), "------------------------------------------------------------"
    end if

    !read in our section of the namelist
    open(newunit=unit, file=nml)
    read(unit, nml=g3dv_bgcov)
    if(isroot) print g3dv_bgcov

    rewind(unit)
    read(unit, nml=g3dv_hzcor)
    if(isroot) print g3dv_hzcor
    
    rewind(unit)
    read(unit, nml=g3dv_vtcor)
    if(isroot) print g3dv_vtcor

    close(unit)

    


    ! determine horizontal correlation lengths
    if(isroot) then
       print *, ""
       print *, " Horizontal correlation length scales summary:"
       print *, " Latitude         x-scale(km)     y-scale(km)"
       r = 0
       do while(r <= 90)
          print *, r, bgcov_hzdist(r)
          r = r + 2.5
       end do
    end if



    ! calculate the vertical correlation lengths
    ! ------------------------------------------------------------
    if(isroot) then
       print *,""
       print *, "calculating  vertical localization values..."
    end if
    allocate(bgcov_local_vtloc(grid_nz, g3dv_mpi_ijcount))
    bgcov_local_vtloc = 0.0
    do i=1,g3dv_mpi_ijcount
       if (grid_local_mask(i) <= 0) CYCLE
       do j=1,grid_nz
          if(grid_dpth(j) > grid_local_D(i)) EXIT
       end do
       bgcov_local_vtloc(:,i) = col_vtloc_mld(grid_local_rho(:,i), grid_dpth(:), j-1)
    end do

    
    ! TODO load the t/s surf min variance
    
    ! calculate  the background error variance
    !------------------------------------------------------------    
    if (isroot) then
       print *, ""
       print *, "Loading background error variance..."
    end if
    allocate(bgcov_local_var_t(grid_nz, g3dv_mpi_ijcount))
    allocate(bgcov_local_var_s(grid_nz, g3dv_mpi_ijcount))
    bgcov_local_var_t = 0.0
    bgcov_local_var_s = 0.0

    do i=1,g3dv_mpi_ijcount
       if (grid_local_mask(i) <= 0) CYCLE
       do j=1,grid_nz
          if(grid_dpth(j) > grid_local_D(i)) EXIT
       end do
       call col_bgvar(grid_local_temp(:,i), grid_dpth(:), j-1, t_varmin_surf_const, s_varmin_surf_const,&
            bgcov_local_var_t(:,i), bgcov_local_var_s(:,i))
    end do


    ! smooth the background error variance and vertical correlation lengths
    !--------------------------------------------------------------------------------
    timer = timer_init(" bgcov_smooth", TIMER_SYNC)
    call timer_start(timer)

    call g3dv_mpi_ij2grd_real(grid_local_D, global_depth, -1)
    call g3dv_mpi_ij2grd_real(grid_local_lat, global_lat, -1)
    call g3dv_mpi_ij2grd_real(grid_local_lon, global_lon, -1)
    if(isroot) print *, "SMoothing fields..."
    do i =1, gauss_iter
       ! vt smooth
       call smooth_vt(bgcov_local_vtloc, bgcov_local_vtloc)
              
       !hz smooth
       call g3dv_mpi_zij2xyk_real(bgcov_local_vtloc, work_xyk)       
       do j=1,g3dv_mpi_kcount
          d=g3dv_mpi_rank+1 + (j-1)*g3dv_mpi_size          
          global_mask = (global_depth < grid_dpth(d))
          call smooth_hz(work_xyk(:,:,j), global_mask, global_lat, global_lon)
       end do       
       call g3dv_mpi_xyk2zij_real(work_xyk, bgcov_local_vtloc)
    end do


    do i =1, gauss_iter
       ! vt smooth
       call smooth_vt(bgcov_local_var_t, bgcov_local_vtloc)              

       ! hz  smooth
       call g3dv_mpi_zij2xyk_real(bgcov_local_var_t, work_xyk)       
       do j=1,g3dv_mpi_kcount
          d=g3dv_mpi_rank+1 + (j-1)*g3dv_mpi_size          
          global_mask = (global_depth < grid_dpth(d))
          call smooth_hz(work_xyk(:,:,j), global_mask, global_lat, global_lon)
       end do
       call g3dv_mpi_xyk2zij_real(work_xyk, bgcov_local_var_t)

       !vt smooth
       call smooth_vt(bgcov_local_var_s, bgcov_local_vtloc)       

       ! hz  smooth
       call g3dv_mpi_zij2xyk_real(bgcov_local_var_s, work_xyk)       
       do j=1,g3dv_mpi_kcount
          d=g3dv_mpi_rank+1 + (j-1)*g3dv_mpi_size          
          global_mask = (global_depth < grid_dpth(d))
          call smooth_hz(work_xyk(:,:,j), global_mask, global_lat, global_lon)
       end do
       call g3dv_mpi_xyk2zij_real(work_xyk, bgcov_local_var_s)

    end do
    


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

    real :: dist0, r, vtloc_d, r2(2), distx, disty, lat0
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
       r = sin(ob1%lat*d2r)*sin(ob2%lat*d2r) + &
           cos(ob1%lat*d2r)*cos(ob2%lat*d2r)*&
           cos((ob1%lon-ob2%lon)*d2r)
       if(r > 1)  r =  1
       if(r < -1) r = -1
       dist0 = re*acos(r)
    end if

    ! horizontal localization
    lat0=(ob1%lat + ob2%lat)/2.0
    r2=bgcov_hzdist(lat0)
    if(r2(1) > r2(2)) then
       distx=abs(ob1%lon - ob2%lon)
       if (distx > 180) distx = 360.0 - distx
       distx = distx * d2r * re * cos(abs(lat0)*d2r)
       disty = abs(ob1%lat - ob2%lat)
       disty = disty * d2r * re
       hz_cor = loc(distx, r2(1)) * loc(disty, r2(2))
    else
       hz_cor = loc(dist0, r2(2))
    end if

    if(hz_cor <= 0) return



    ! vertical localization distance is the average of vt_loc of the two points 
!    vtloc_d = sqrt(ob1%grd_vtloc*ob2%grd_vtloc)
    vtloc_d = (ob1%grd_vtloc+ob2%grd_vtloc)/2.0
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
    real :: vtloc_d, distx,disty

    integer :: i, idx_start
    real :: r, r2(2), lat0


    cor = 0.0
    var = 0.0


    ! gradient tensor localizations
    surf_tensor = merge( &
               loc(abs(ob%grd_ssh - grid_local_ssh(ij)), tnsr_surf),&
               1.0, tnsr_surf > 0.0)
    if(surf_tensor <= 0.0) return


    ! horizontal localization    
    lat0=(grid_local_lat(ij) + ob%lat)/2.0
    r2=bgcov_hzdist(lat0)
    if(r2(1) > r2(2)) then
       distx=abs(grid_local_lon(ij) - ob%lon)
       if (distx > 180) distx = 360.0 - distx
       distx = distx * d2r * re * cos(abs(lat0)*d2r)
       disty = abs(grid_local_lat(ij) - ob%lat)
       disty = disty * d2r * re
       hz_cor = loc(distx, r2(1)) * loc(disty, r2(2))
    else
       hz_cor = loc(dist, r2(2))
    end if
    if (hz_cor <= 0.0) return


    ! vertical localization distance is the average of vt_loc of the two points 
    vt_cor = 0.0    
    do i = 1, grid_nz
       if(bgcov_local_vtloc(i,ij) <= 0.0) exit
!       vtloc_d=sqrt(ob%grd_vtloc*bgcov_local_vtloc(i,ij))
       vtloc_d=(ob%grd_vtloc+bgcov_local_vtloc(i,ij))/2.0
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
    real :: cor(2)
    
    ! correlation length based on rossby radius calculation
    if ( abs(lat) < 0.1) then
       cor = hzcor_max
    else
       cor = max(hzcor_min, min(hzcor_max, hzcor_cspd/(2*omega*abs(sin(lat*d2r))) ))
    end if

    ! stretching along equator
    cor(1) = cor(1)*(1 + (hzcor_stretch_mul-1.0)*loc_gc(lat, hzcor_stretch_lat * 0.5* sqrt(0.3)))

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


  pure function col_vtloc_mld(dens, dpth, bottom) result(vtloc)
    real,    intent(in) :: dens(:)
    real,    intent(in) :: dpth(:)
    integer, intent(in) :: bottom

    real :: vtloc(size(dens))

    real    :: mld_max 
    real    :: mld
    integer :: mld_z
    integer :: z2, z
    real    :: r

    vtloc=0.0
    mld_max = vtcor_max*2

    ! calculate vt loc distances from layer thickensses
    vtloc(1) = (dpth(2)-dpth(1))
    do z=2,size(dpth)-1
       vtloc(z) = (dpth(z+1)-dpth(z-1))/2.0
    end do
    vtloc(size(dpth)) = vtloc(size(dpth)-1)
    vtloc(bottom+1:size(dpth)) = 0.0


    ! calculate the mixed layer depth
    mld = -1
    mld_z = bottom
    do z2 = 2, bottom
       if(dens(z2)-dens(1) >= rho_delta) then                    
          r = (dens(1)+rho_delta-dens(z2-1)) / (dens(z2)-dens(z2-1))
          mld = dpth(z2-1) + ( dpth(z2)-dpth(z2-1))*r
          mld_z = z2-1
          exit
       end if       
       if(mld_max > 0 .and. dpth(z2) > mld_max) then
          mld_z = z2 -1
          mld = mld_max
          exit
       end if
    end do
    if(mld < 0) mld = dpth(bottom)
    if(mld_max > 0)  mld = min(mld_max, mld)

    ! linear interpolation of MLD from top to bottom of MLD
    ! value at top is *half* the MLD so that decorrelates at based of ML
    ! TODO check this logic ^
    do z = 1, mld_z
       r = 1.0 - (dpth(z) / mld)
       vtloc(z) = (mld/2.0)*r + vtloc(mld_z)*(1.0-r)
    end do    
  end function col_vtloc_mld


  !================================================================================
  !================================================================================

  pure subroutine col_bgvar(temp, dpth,  bottom, t_varmin, s_varmin, var_t, var_s)
    real,    intent(in) :: temp(:)
    real,    intent(in) :: dpth(:)
    integer, intent(in) :: bottom
    real,    intent(in) :: t_varmin, s_varmin
    real,    intent(out) :: var_t(:), var_s(:)

    integer :: z

    real :: t_r1(grid_nz), s_r1(grid_nz)
    real :: r2(grid_nz)
    
    ! calculate profile minimums
    do z=1,grid_nz
       t_r1(z) = t_varmin_do + (t_varmin - t_varmin_do)*exp((dpth(1) - dpth(z))/ t_varmin_len)
       s_r1(z) = s_varmin_do + (s_varmin - s_varmin_do)*exp((dpth(1) - dpth(z))/ s_varmin_len)
    end do

    !  calculate dT/dZ
    ! r2 = error calculated from dt/dz
    r2 = 0.0
    r2(1) = (temp(2)-temp(1)) / (dpth(2)-dpth(1))
    do z = 2, bottom-1
       r2(z) = (temp(z+1)-temp(z-1)) / (dpth(z+1)-dpth(z-1))
    end do
    r2(bottom) = (temp(bottom)-temp(bottom-1)) / (dpth(bottom)-dpth(bottom-1))

    !combine calculated dt/dz with the calculated minimums
    t_r1(1:bottom) = max(min(abs(r2(1:bottom))*t_dz,t_varmax), t_r1(1:bottom))
    s_r1(1:bottom) = max(min(abs(r2(1:bottom))*s_dz,s_varmax), s_r1(1:bottom))

    var_t(1:bottom)=t_r1(1:bottom)
    var_s(1:bottom)=s_r1(1:bottom)
  end subroutine col_bgvar


  !================================================================================
  !================================================================================
  

  subroutine smooth_vt(data, vtcor)
    real,    intent(inout) :: data(grid_nz, g3dv_mpi_ijcount)
    real,    intent(in)    :: vtcor(grid_nz, g3dv_mpi_ijcount)

    real :: data2(grid_nz, g3dv_mpi_ijcount)

    integer :: ij, z, i, z2
    real :: gb_s, gb_v, r, d


    data2 = data
    data = 0.0

    do ij=1,g3dv_mpi_ijcount
       do z= 1, grid_nz
          if(grid_local_D(ij) < grid_dpth(z)) exit

          gb_s=1.0
          gb_v=data2(z,ij)

          ! do the following loop twice, once for up (i==1) and down (i==-1)
          do i=1,-1,-2 
             z2 = z + i
             do while(z2 <= grid_nz .and. z2 >= 1)
                d = abs(grid_dpth(z)-grid_dpth(z2))
                r=loc(d, (vtcor(z,ij)+vtcor(z2,ij))/2.0)
!                if (r==0) exit !ERR, should never happen?
                if(grid_local_D(ij) < grid_dpth(z2)) exit
                gb_v=gb_v+data2(z2,ij)*r
                gb_s=gb_s+r                   
                z2 = z2 + i
             end do
          end do
          
          data(z,ij) = gb_v/gb_s
       end do
    end do
  end subroutine smooth_vt

  !================================================================================
  !================================================================================
  

  subroutine smooth_hz(data, mask, lat, lon)
    real,    intent(inout) :: data(:,:)
    logical, intent(in)    :: mask(:,:)
    real,    intent(in)    :: lat(:,:)
    real,    intent(in)    :: lon(:,:)

    real :: data2(grid_nx,grid_ny)

    integer :: y, x, i, x2, y2
    real :: gb_s, gb_v, r, r2(2),r2_2(2),  d

    
    data2 = data
    data = 0.0

    ! 1 iteration of recursive filter in the x direction
    do y=1, grid_ny
       do x=1, grid_nx
          if(mask(x,y)) cycle
          gb_s = 1.0
          gb_v = data2(x,y)

          ! do the following loop twice, once for right (i==1) and left (i==-1)
          do i=1,-1,-2
             x2 = x + i
             do while(x2 /= x)
                ! wrap in the x direction
                if(i==1) then
                   if (x2 > grid_nx) x2 = 1
                else
                   if (x2 < 1) x2 = grid_nx
                end if

                r2=bgcov_hzdist(lat(x,y))
                d=abs(lon(x,y) - lon(x2,y))
                if(d > 180) d = 360.0 -d
                d = d*d2r*re*cos(abs(lat(x,y)*d2r))
                r = loc(d, r2(1))
                if (r == 0) exit
                if(.not. mask(x2,y)) then
                   gb_v = gb_v + data2(x2,y)*r
                   gb_s = gb_s + r
                end if
                x2 = x2 + i
             end do
             
             ! if we looped back to the original point, dont bother doing the loop
             ! in the other direction
             if(x2 == x) exit
             
          end do
          data(x,y) = gb_v/gb_s
       end do
    end do
 

    ! 1 iteration of recursive filter in y direction
    data2 = data
    data = 0.0
    do y=1, grid_ny
       do x=1, grid_nx
          if(mask(x,y)) cycle
          gb_s=1.0
          gb_v=data2(x,y)

          ! do the following loop twice, once for up (i==1) and down (i==-1)
          do i=1,-1,-2 
             y2 = y + i
             do while(y2 <= size(data,2) .and. y2 >= 1)
                ! TODO, averaging of lat is wrong?
!                lat0=(lat(x,y)+lat(x,y2))/2.0
                !                r2=bgcov_hzdist(lat0)
                r2=bgcov_hzdist(lat(x,y))
                r2_2=bgcov_hzdist(lat(x,y2))
                r2=(r2+r2_2)/2.0
                d=abs(lat(x,y) - lat(x,y2))
                d=d*d2r*re
                r=loc(d,r2(2))                
                if (r == 0) exit
                if(.not. mask(x,y2)) then
                   gb_v = gb_v + data2(x,y2)*r
                   gb_s = gb_s + r
                end if
                y2 = y2 + i
             end do
          end do
          data(x,y) = gb_v/gb_s
        end do
     end do

  end subroutine smooth_hz
  
end module  g3dv_bgcov
