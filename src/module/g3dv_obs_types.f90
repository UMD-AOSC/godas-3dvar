module g3dv_obs_types
  !! author: Travis Sluka
  !!
  !! g3dv_obs_types description
  !!
  !! g3dv_obs_types detailed description
  !!


  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: observation


  ! public module variables
  !------------------------------------------------------------


  ! private module variables
  !------------------------------------------------------------

  type observation
     integer :: id = -1
     real :: lat
     real :: lon
     real :: dpth
     real :: hr
    
     real :: inc
     real :: err

     integer :: qc
  end type observation


contains


  !================================================================================
  !================================================================================



!  subroutine
!  end subroutine



  !================================================================================
  !================================================================================


end module  g3dv_obs_types
