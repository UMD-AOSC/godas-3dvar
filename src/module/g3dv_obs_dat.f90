module g3dv_obs_dat
  use g3dv_obs

  implicit none
  private


  type, public, extends(obsio) :: obsio_dat
   contains
     procedure :: get_name  => obsio_get_name
     procedure :: get_desc  => obsio_get_desc
     procedure :: write     => obsio_dat_write
     procedure :: read      => obsio_dat_read
  end type obsio_dat


contains


  function obsio_get_name(self)
    class(obsio_dat) :: self
    character(:), allocatable :: obsio_get_name
    obsio_get_name = "OBSIO_DAT"
    self%i = self%i
  end function obsio_get_name



  function obsio_get_desc(self)
    class(obsio_dat) :: self
    character(:), allocatable :: obsio_get_desc
    obsio_get_desc = "raw observation I/O format"
    self%i = self%i
  end function obsio_get_desc



  subroutine obsio_dat_write(self, file, obs, obs_qc)
    class(obsio_dat) :: self
    character(len=*),  intent(in)  :: file
    type(observation), intent(in)  :: obs(:)
    integer,           intent(in)  :: obs_qc(:)

    self%i = self%i
    if(size(obs) == size(obs)) continue

    print *, "writing to file ",trim(file)
    print *, "ERROR: obsio write method not yet implemented"

    stop 1
  end subroutine obsio_dat_write


  subroutine obsio_dat_read(self, file, obs)
    class(obsio_dat) :: self
    character(len=*),  intent(in) :: file
    type(observation), intent(out), allocatable :: obs(:)

    integer :: filesize, unit, i
    real(kind=4) :: record(7)
    
    self%i = self%i

    ! determien the number of observatiosn that will be read in
    inquire(file=file, size=filesize)
    allocate( obs(filesize/4/9))

    ! open the file and read in observations
    open(newunit=unit, file=file, form='unformatted', access='sequential', action='read')
    do i = 1, size(obs)
       read(unit) record
       obs(i)%id   = int(record(1))
       obs(i)%lon  = record(2)
       obs(i)%lat  = record(3)
       obs(i)%dpth = record(4)
       obs(i)%inc  = record(5)
       obs(i)%err  = record(6)
    end do
    
  end subroutine obsio_dat_read
end module g3dv_obs_dat
