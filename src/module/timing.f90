!! author: Travis Sluka
!! category: support
!! multiprocessor MPI timing statistics

module timing
  !! Provides timing statistics for segements of code across MPI processes.
  !!
  !! Various timers can be created using [[timing_init]], these timers are then
  !! started with [[timer_start]] and stopped with [[timer_stop]]. Timers can be
  !! started and stopped multiple times, and the cumulative times are stored.
  !! At the end of a program's run, [[timer_print]] can be called to display statistics
  !! for min/max/mean/stddev of run times across the mpi processes.

  use mpi

  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: timing_init
  public :: timer_init
  public :: timer_start
  public :: timer_stop
  public :: timer_print

  ! public module variables
  !------------------------------------------------------------
  integer, parameter, public :: TIMER_SYNC = 1


  ! private module variables
  !------------------------------------------------------------
  type timer_obj
     integer(kind=8)     :: total_ticks
     integer(kind=8)     :: tick
     character(len=1024) :: name
     integer :: flags
     integer :: grain
  end type timer_obj

  integer, parameter :: max_timers = 1024
  integer            :: active_timers = 0
  type(timer_obj)    :: timer_objs(max_timers)
  integer :: mp_comm, mp_root, mp_rank, mp_size



contains




  !================================================================================
  !================================================================================



  subroutine timing_init(comm, root)
    integer, intent(in) :: comm
    integer, intent(in) :: root
    integer :: ierr
    mp_root = root
    mp_comm = comm
    call mpi_comm_rank(mp_comm, mp_rank, ierr)
    if(ierr /= 0) then
       print *,"ERROR in timing_init"
       stop 1
    end if
    call mpi_comm_size(mp_comm, mp_size, ierr)
    if(ierr /= 0) then
       print *,"ERROR in timing_init"
       stop 1
    end if
  end subroutine timing_init




  !================================================================================
  !================================================================================



  function gettimer(timer) result(id)
    character(len=*), intent(in) :: timer
    integer :: id

    ! see if the timer already exists
    id = 1
    do while (id <= active_timers)
       if (trim(timer_objs(id)%name) == trim(timer)) exit
       id = id + 1
    end do

    if (id > active_timers) id = -1
  end function gettimer



  !================================================================================
  !================================================================================



  function timer_init(name, flags) result(id)
    character(len=*), intent(in) :: name
    integer, optional, intent(in) :: flags
    integer :: id

    id = gettimer(name)

    if ( id < 0 ) then
       active_timers = active_timers + 1
       id = active_timers

       if (active_timers > max_timers) then
          print *, "ERROR, too many timers have been created, increase max_timers"
          stop 1
       end if
       timer_objs(id)%name = trim(name)
       timer_objs(id)%total_ticks = 0
       timer_objs(id)%tick = 0
       timer_objs(id)%flags = 0
       if (present(flags)) timer_objs(id)%flags = flags
    end if
  end function timer_init



  !================================================================================
  !================================================================================



  subroutine timer_start(id)
    integer, intent(in) :: id
    integer :: ierr

    if (id <= 0) return

    if (timer_objs(id)%tick > 0) then
       print *, "ERROR: trying to start a timer that is already started :" , trim(timer_objs(id)%name)
       stop 1
    end if

    if ( iand(timer_objs(id)%flags, TIMER_SYNC) > 0)then
       call mpi_barrier(mp_comm, ierr)
    end if

    call system_clock(timer_objs(id)%tick)
  end subroutine timer_start



  !================================================================================
  !================================================================================



  subroutine timer_stop(id)
    integer, intent(in) :: id

    integer(kind=8):: tick

    if (id <= 0) return

    if (timer_objs(id)%tick <= 0) then
       print *, "ERROR: trying to stop a timer that is not started: ", trim(timer_objs(id)%name)
       stop 1
    end if

    call system_clock(tick)
    timer_objs(id)%total_ticks = timer_objs(id)%total_ticks + &
         (tick - timer_objs(id)%tick)
    timer_objs(id)%tick = 0
  end subroutine timer_stop



  !================================================================================
  !================================================================================



  ! subroutine timer_gather(id, times)
  !   integer, intent(in) :: id
  !   real,intent(inout) :: times(:)
  !   integer(kind=8) :: timer_rate
  !   real :: t0
  !   integer :: ierr

  !   !! todo, error checking on valid timer id
  !   call system_clock(count_rate=timer_rate)
  !   t0 = real(timer_objs(id)%total_ticks) / real(timer_rate)
  !   call mpi_allgather(t0, 1, mpi_real, times, 1, mpi_real, mp_comm, ierr)
  ! end subroutine timer_gather



  !================================================================================
  !================================================================================



  subroutine timer_print
    integer :: i
    real :: t, tmin, tmax, tave, tdif2, tdif2g

    integer(kind=8) :: timer_rate
    integer :: ierr

    if (mp_root == mp_rank) then
       print *,""
       print *,"Timing:"
       print *,"============================================================"
       print '(A,I4,A)',"calculating timer statistics across ", mp_size," processes..."
       print *, ""
       print '(A18,4A10)',"","ave","min","max", "std"
    end if

    call system_clock(count_rate=timer_rate)
    i = 1
    do while (i <= active_timers)
       t = real(timer_objs(i)%total_ticks)/real(timer_rate)

       call mpi_reduce   (t, tmin, 1, mpi_real, mpi_min, mp_root, mp_comm, ierr)
       call mpi_reduce   (t, tmax, 1, mpi_real, mpi_max, mp_root, mp_comm, ierr)
       call mpi_allreduce(t, tave, 1, mpi_real, mpi_sum,          mp_comm, ierr)
       tave = tave / mp_size
       tdif2 = (t-tave)*(t-tave)
       call mpi_reduce (tdif2, tdif2g, 1, mpi_real, mpi_sum, mp_root, mp_comm, ierr)

       if (mp_root == mp_rank) then
          tdif2g = sqrt(tdif2g/mp_size)
          print '(A,A18,4F10.1)', " ",timer_objs(i)%name, tave, tmin, tmax, tdif2g
       end if

       i = i +1
    end do
  end subroutine timer_print




  !================================================================================


end module timing
