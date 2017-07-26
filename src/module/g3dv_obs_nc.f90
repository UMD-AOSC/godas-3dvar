module g3dv_obs_nc
  use g3dv_obs
  use netcdf

  implicit none
  private


  type, public, extends(obsio) :: obsio_nc
   contains
     procedure :: get_name  => obsio_get_name
     procedure :: get_desc  => obsio_get_desc
     procedure :: write     => obsio_nc_write
     procedure :: writeqc   => obsio_nc_writeqc
     procedure :: read      => obsio_nc_read
  end type obsio_nc


contains


  function obsio_get_name(self)
    class(obsio_nc) :: self
    character(:), allocatable :: obsio_get_name
    obsio_get_name = "OBSIO_NC"
    self%i = self%i
  end function obsio_get_name



  function obsio_get_desc(self)
    class(obsio_nc) :: self
    character(:), allocatable :: obsio_get_desc
    obsio_get_desc = "netCDF observation I/O format"
    self%i = self%i
  end function obsio_get_desc



  subroutine obsio_nc_write(self, file, obs, obs_qc)
    class(obsio_nc) :: self
    character(len=*),  intent(in)  :: file
    type(observation), intent(in)  :: obs(:)
    integer, intent(in) :: obs_qc(:)

    self%i = self%i
    if(size(obs) == size(obs)) continue

    print *, "writing to file ",trim(file)
    print *, "ERROR: obsio write method not yet implemented"

    stop 1
  end subroutine obsio_nc_write


  subroutine obsio_nc_writeqc(self, file, obs_qc)
    class(obsio_nc) :: self
    character(len=*),  intent(in) :: file
    integer,           intent(in) :: obs_qc(:)

    integer :: ncid, dimid, varid

    call check( nf90_create(file, nf90_clobber, ncid))
    call check( nf90_def_dim(ncid, "obs", nf90_unlimited, dimid))
    call check( nf90_def_var(ncid, "qc",  nf90_int, dimid, varid))
    call check( nf90_put_att(ncid, varid, "descrtiption", "quality control flag returned  by 3dvar"))
    call check( nf90_enddef(ncid))
    call check( nf90_put_var(ncid,varid, obs_qc))
    call check( nf90_close(ncid))   
  end subroutine obsio_nc_writeqc



  subroutine obsio_nc_read(self, file, obs)
    class(obsio_nc) :: self
    character(len=*),  intent(in) :: file
    type(observation), intent(out), allocatable :: obs(:)

    integer :: i
    integer :: ncid,dim_id
    integer :: vid_obid, vid_lon, vid_lat, vid_depth, vid_inc, vid_err, vid_hr
    integer :: nobs
    integer, allocatable :: tmp_i(:)
    real,    allocatable :: tmp_r(:)    

    self%i = self%i

    call check(nf90_open(file, nf90_nowrite, ncid))
    call check(nf90_inq_dimid(ncid, "obs",   dim_id))
    call check(nf90_inquire_dimension(ncid, dim_id, len=nobs))
    call check(nf90_inq_varid(ncid, "obid",  vid_obid))
    call check(nf90_inq_varid(ncid, "lon",   vid_lon))
    call check(nf90_inq_varid(ncid, "lat",   vid_lat))
    call check(nf90_inq_varid(ncid, "depth", vid_depth))
    call check(nf90_inq_varid(ncid, "val",   vid_inc))
    call check(nf90_inq_varid(ncid, "err",   vid_err))
    call check(nf90_inq_varid(ncid, "hr",   vid_hr))

    print *, "Observations in file: ", nobs
    allocate(obs(nobs))
    allocate(tmp_i(nobs))
    allocate(tmp_r(nobs))

    call check(nf90_get_var(ncid, vid_obid, tmp_i))
    do i=1,nobs;   obs(i)%id = tmp_i(i);   end do

    call check(nf90_get_var(ncid, vid_lon, tmp_r))
    do i=1,nobs;   obs(i)%lon = tmp_r(i);   end do

    call check(nf90_get_var(ncid, vid_lat, tmp_r))
    do i=1,nobs;   obs(i)%lat = tmp_r(i);   end do

    call check(nf90_get_var(ncid, vid_depth, tmp_r))
    do i=1,nobs;   obs(i)%dpth = tmp_r(i);   end do

    call check(nf90_get_var(ncid, vid_inc, tmp_r))
    do i=1,nobs;   obs(i)%inc = tmp_r(i);   end do

    call check(nf90_get_var(ncid, vid_err, tmp_r))
    do i=1,nobs;   obs(i)%err = tmp_r(i);   end do

    call check(nf90_get_var(ncid, vid_hr, tmp_r))
    do i=1,nobs;   obs(i)%hr = tmp_r(i);   end do


    deallocate(tmp_i)
    deallocate(tmp_r)
    call check(nf90_close(ncid))
  end subroutine obsio_nc_read


  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check 
end module g3dv_obs_nc
