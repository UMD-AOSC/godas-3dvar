module g3dv_mpi
  !! author: Travis Sluka
  !!
  !! g3dv_mpi description
  !!
  !! g3dv_mpi detailed description
  !!

  use mpi
  use iso_fortran_env
#ifdef __INTEL_COMPILER
  use ifport
#endif

  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: g3dv_mpi_init
  public :: g3dv_mpi_setgrid
  public :: g3dv_mpi_final

  public :: g3dv_mpi_barrier

  public :: g3dv_mpi_bcstflag
  public :: g3dv_mpi_grd2ij_real
  public :: g3dv_mpi_ij2grd_real



  ! public module variables
  !------------------------------------------------------------
  logical, public, protected :: g3dv_mpi_isroot
  integer, public, protected :: g3dv_mpi_root
  integer, public, protected :: g3dv_mpi_size
  integer, public, protected :: g3dv_mpi_rank
  integer, public, protected :: g3dv_mpi_comm

  integer, public, protected :: g3dv_mpi_ijcount



  ! private module variables
  !------------------------------------------------------------
  integer, allocatable :: scatterv_count(:)
  integer, allocatable :: scatterv_displ(:)

  integer, allocatable :: ij_list(:)

  integer :: grid_nx
  integer :: grid_ny
  integer :: grid_nz

  logical, parameter :: interleave = .true.




contains


  !================================================================================
  !================================================================================



  subroutine g3dv_mpi_init(m_comm)
    !! initialize MPI if it has not already been done

    integer, intent(in), optional :: m_comm
    !! if MPI has been initialized outside this module already, use the given
    !! MPI_COMM, otherwise MPI_COMM_WORLD will be used

    integer :: ierr
    logical :: initialized


    ! check to see if mpi has already been initialized
    call MPI_Initialized(initialized, ierr)

    ! it it has not been initialized, do so
    if (.not. initialized) then
       if(present(m_comm)) then
          print *, "ERROR: an MPI_COMM has been passed to g3dv_mpi_init", &
               " but MPI has not been initialized"
          stop 1
       end if
       call mpi_init(ierr)
       g3dv_mpi_root = 0
       g3dv_mpi_comm = mpi_comm_world
    else
       if(present(m_comm)) then
          g3dv_mpi_comm = m_comm
       else
          g3dv_mpi_comm = mpi_comm_world
       end if
       ! TODO, allow specification of root proc if MPI has already been initialized
       g3dv_mpi_root = 0
    end if

    call mpi_comm_size(g3dv_mpi_comm, g3dv_mpi_size, ierr)
    call mpi_comm_rank(g3dv_mpi_comm, g3dv_mpi_rank, ierr)

    g3dv_mpi_isroot = g3dv_mpi_root == g3dv_mpi_rank

    if(g3dv_mpi_isroot) then
       if (initialized) then
          print *, "MPI already initialized"
       else
          print *, "initializing MPI..."
       end if
       print *,"mpi_size =",g3dv_mpi_size
       print *,"mpi_root =",g3dv_mpi_root
       if (g3dv_mpi_comm == mpi_comm_world) then
          print *, "mpi_comm = MPI_COMM_WORLD"
       else
          print *, "mpi_comm =", g3dv_mpi_comm
       end if
    end if

  end subroutine




  !================================================================================
  !================================================================================



  subroutine g3dv_mpi_setgrid(nx, ny, nz)
    integer, intent(in) :: nx, ny, nz
    integer :: i,j
    integer :: count, prev


    if(g3dv_mpi_isroot) then
       print *, new_line('a'),&
            new_line('a'),"------------------------------------------------------------",&
            new_line('a'), "g3dv_mpi_setgrid", &
            new_line('a'), "------------------------------------------------------------"
       print '(A,I0,A,I0,A,I0)', "  MPI setting grid to ",nx," x ",ny," x ",nz
    end if

    grid_nx = nx
    grid_ny = ny
    grid_nz = nz


    ! calculate the number of gridpoints to use for each proc

    allocate(scatterv_count(g3dv_mpi_size))
    allocate(scatterv_displ(g3dv_mpi_size))
    prev = 0
    do i = 0, g3dv_mpi_size-1
       count = nint(nx*ny*(1.0/g3dv_mpi_size))
       if (i == g3dv_mpi_size-1) count = nx*ny - prev
       if (i == g3dv_mpi_rank) then
          g3dv_mpi_ijcount = count
          allocate(ij_list(count))
          do j = 1, count
             ij_list(j) = prev + j
          end do
       end if
       scatterv_count(i+1) = count
       scatterv_displ(i+1) = prev
       prev = prev + count
    end do

    if (g3dv_mpi_isroot) then
       print *, ""
       i = minval(scatterv_count)
       j = maxval(scatterv_count)
       if(i==j) then
          print *, " all procs assigned ",i,"gridpoints"
       else
          print *, " all procs assigned between ",i,"and",j,"gridpoints"
       end if
    end if

  end subroutine g3dv_mpi_setgrid



  !================================================================================
  !================================================================================



  subroutine g3dv_mpi_grd2ij_real(grd, ij)
    !! takes a single grid on the root process and distributes portions of it to 
    !! worker processes

    real, intent(in) :: grd(grid_nx*grid_ny)
    real, intent(inout) :: ij(g3dv_mpi_ijcount)

    integer :: ierr, i, j 
    real :: wrk(grid_nx*grid_ny)

    if (.not. interleave) then
       call mpi_scatterv(grd, scatterv_count, scatterv_displ, mpi_real, &
            ij, g3dv_mpi_ijcount, mpi_real, g3dv_mpi_root, g3dv_mpi_comm, ierr)
    else
       if (g3dv_mpi_isroot) then
          do i = 1, g3dv_mpi_size
             do j = 1, scatterv_count(i)
                wrk(scatterv_displ(i)+j) = grd((j-1)*g3dv_mpi_size+i)
             end do
          end do
       end if
       call mpi_scatterv(wrk, scatterv_count, scatterv_displ, mpi_real, &
            ij, g3dv_mpi_ijcount, mpi_real, g3dv_mpi_root, g3dv_mpi_comm, ierr)
    end if
  end subroutine g3dv_mpi_grd2ij_real



  !================================================================================
  !================================================================================



  subroutine g3dv_mpi_ij2grd_real(ij, grd, rank)
    real, intent(in)    :: ij(g3dv_mpi_ijcount)
    real, intent(inout) :: grd(grid_nx*grid_ny)
    integer, intent(in), optional :: rank
    
    integer :: rank0
    integer :: ierr, i , j
    real :: wrk(grid_nx*grid_ny)

    if (present(rank)) then
       rank0 = rank
    else
       rank0 = g3dv_mpi_root
    end if

    if (.not. interleave) then
       call mpi_gatherv(ij, g3dv_mpi_ijcount, mpi_real, &
            grd, scatterv_count, scatterv_displ, mpi_real, rank0, g3dv_mpi_comm, ierr)
    else
       call mpi_gatherv(ij, g3dv_mpi_ijcount, mpi_real, &
            wrk, scatterv_count, scatterv_displ, mpi_real, rank0, g3dv_mpi_comm, ierr)
       if(rank0 == g3dv_mpi_rank) then
          do i = 1, g3dv_mpi_size
             do j = 1, scatterv_count(i)
                grd((j-1)*g3dv_mpi_size+i) = wrk(scatterv_displ(i)+j)
             end do
          end do
       end if
    end if
  end subroutine g3dv_mpi_ij2grd_real



  !================================================================================
  !================================================================================



  subroutine g3dv_mpi_final
    integer :: ierr
    real    :: mem, mem_avg, mem_min, mem_max, mem_sum

    ! calculate memory hiwater mark across all procs
    call getMaxMem(mem)
    call mpi_reduce(mem, mem_avg, 1, mpi_real, mpi_sum, g3dv_mpi_root, g3dv_mpi_comm, ierr)
    call mpi_reduce(mem, mem_min, 1, mpi_real, mpi_min, g3dv_mpi_root, g3dv_mpi_comm, ierr)
    call mpi_reduce(mem, mem_max, 1, mpi_real, mpi_max, g3dv_mpi_root, g3dv_mpi_comm, ierr)
    call mpi_reduce(mem, mem_sum, 1, mpi_real, mpi_sum, g3dv_mpi_root, g3dv_mpi_comm, ierr)
    mem_avg = mem_avg / g3dv_mpi_size
    if(g3dv_mpi_isroot) then
       print *, new_line('A'),&
            new_line('A'), "------------------------------------------------------------",&
            new_line('A'), "Memory hiwater mark: "
       print *, "per core:"
       print '(A,F5.2,A)', "   avg: ",mem_avg," GB"
       print '(A,F5.2,A)', "   min: ",mem_min," GB"
       print '(A,F5.2,A)', "   max: ",mem_max," GB"
       print *,"total:"
       print '(A,f5.2,A)', "  ", mem_sum," GB"
    end if

    !TODO, don't finalize if MPI was initialized BEFORE calling g3dv_mpi_init
    call mpi_finalize(ierr)

  end subroutine g3dv_mpi_final



  !================================================================================
  !================================================================================



  subroutine getMaxMem(mem)
    !! return memory hiwater mark, in gigabytes.
    !! this relies on the VmHWM field set in the /proc/(pid)/status file
    !! if this file or field can't be found, -1 is returned

    real, intent(out)   :: mem
    character(len=30)   :: pid_char
    character(len=1024) :: proc_file, line
    integer             :: unit, iostat, i
    logical             :: ex

    ! get process id
    write(pid_char,'(I10)') getpid()
    pid_char = adjustl(pid_char)

    ! make sure the file exists that contains memory hi usage info in linux
    proc_file="/proc/"//trim(pid_char)//"/status"
    inquire(file=proc_file, exist=ex)
    if (.not. ex) then
       print *, "can't open ",trim(proc_file)
       mem = -1
       return
    end if

    ! read in the VmHWM line of the file
    open(newunit=unit, file=proc_file, action='read')
    iostat = 0
    do while (iostat==0)
       read(unit,'(A)',iostat=iostat) line
       i = scan(line, ':')
       if(i<=0) cycle
       if(line(:i) == "VmHWM:") then
          ! the VmHWM line has been found, parse out the number of kilobytes
          line = line(i+1:)
          do i=1,len(line)
             if(line(i:i) == char(9)) line(i:i) = ' '
          end do
          line = adjustl(line)
          read(line,*) mem
          ! return the number of gigabytes
          mem = mem/1024/1024
          close(unit)
          return
       end if
    end do
    mem = -1
    close(unit)
  end subroutine getMaxMem



  !================================================================================
  !================================================================================



  subroutine g3dv_mpi_barrier(syncio)
    integer :: ierr
    logical, optional :: syncio

    integer :: i
    call mpi_barrier(g3dv_mpi_comm, ierr)
    
    if(present(syncio) .and. syncio) then
       flush(output_unit)
       i = system('usleep 1')
       call mpi_barrier(g3dv_mpi_comm, ierr)
    end if

  end subroutine g3dv_mpi_barrier

  !================================================================================
  !================================================================================


  function g3dv_mpi_bcstflag(flag) result(res)
    integer, intent(in) :: flag
    integer :: res
    integer :: err

    call mpi_allreduce(flag, res, 1, mpi_integer, mpi_sum, mpi_comm_world, err)
    if (res > 0) res = 1
  end function g3dv_mpi_bcstflag

  !================================================================================
  !================================================================================

end module  g3dv_mpi
