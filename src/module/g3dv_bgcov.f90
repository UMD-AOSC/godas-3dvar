module g3dv_bgcov
  !! author: Travis Sluka
  !!
  !! g3dv_bgcov description
  !!
  !! g3dv_bgcov detailed description
  !!

  use g3dv_obs_types,    only : observation
  use g3dv_grid


  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: bgcov_init
  public :: bgcov_HBH
  public :: bgcov_BH


  ! public module variables
  !------------------------------------------------------------


  ! private module variables
  !------------------------------------------------------------
  real :: hz_cor(2) = -1
  real :: time_sigma = -1
  real :: ssh_grad_tensor = 0.0



contains


  !================================================================================
  !================================================================================



  subroutine bgcov_init(isroot, nml)
    logical, intent(in) :: isroot
    character(len=*), intent(in) :: nml

    integer :: unit

    namelist /g3dv_bgcov/ hz_cor, time_sigma,  ssh_grad_tensor


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

    ! read in background error variance
    ! TODO

    ! computer vertical localization distance field
    ! TODO
  end subroutine bgcov_init




  !================================================================================
  !================================================================================



  pure function bgcov_HBH(ob1, ob2, dist) result(cov)
    type(observation), intent(in) :: ob1
    type(observation), intent(in) :: ob2
    real, optional, intent(in) :: dist
    real :: cov

    cov = 1.0
  end function bgcov_HBH



  !================================================================================
  !================================================================================


  pure function bgcov_BH(ob, dist, ij) result(cov)
    type(observation), intent(in) :: ob
    
    real, intent(in) :: dist

    integer, intent(in) :: ij
    real :: cov(grid_ns)

    cov = 1.0
  end function bgcov_BH
  !================================================================================
  !================================================================================

end module  g3dv_bgcov
