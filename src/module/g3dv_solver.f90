module g3dv_solver
  !! author: Travis Sluka
  !!
  !! g3dv_solver description
  !!
  !! g3dv_solver detailed description
  !!

  use timing
  use g3dv_obs
  use g3dv_mpi
  use g3dv_grid
  use g3dv_bgcov

  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: solver_init
  public :: solver_run


  ! public module variables
  !------------------------------------------------------------


  ! private module variables
  !------------------------------------------------------------

  logical :: isroot

  ! timers
  integer :: timer
  integer :: timer_cholesky
  integer :: timer_pcg
  integer :: timer_precond
  integer :: timer_bgcov_HBH
  integer :: timer_pcg_mpi
  integer :: timer_postx
  integer :: timer_bgcov_BH
  integer :: timer_postx_mpi

  ! variables read in from namelist
  integer :: maxitr = 10
  real    :: conv_ratio = 1e2


contains


  !================================================================================
  !================================================================================



  subroutine solver_init(root, nml)
    logical, intent(in) :: root
    character(len=*), intent(in) :: nml

    integer :: unit
    namelist /g3dv_solver/ maxitr, conv_ratio

    isroot = root

    if(isroot) then
       print *,new_line('a'),&
            new_line('a'), "------------------------------------------------------------",&
            new_line('a'), "g3dv_solveros_init() : preconditioned congjugate gradient solver",&
            new_line('a'), "------------------------------------------------------------"
    end if

    ! read in our section of the namelist
    open(newunit=unit, file=nml)
    read(unit, nml=g3dv_solver)
    close(unit)
    if (isroot) print g3dv_solver

    !setup timer
    timer           = timer_init('(Solver)', TIMER_SYNC)
    timer_cholesky  = timer_init('  cholesky decomp')
    timer_pcg       = timer_init('  PCG')
    timer_precond   = timer_init('    preconditioner')
    timer_bgcov_HBH = timer_init('    bgcov (ob/ob)')
    timer_pcg_mpi   = timer_init('    pcg mpi', TIMER_SYNC)
    timer_postx     = timer_init('  post multiply', TIMER_SYNC)
    timer_bgcov_BH  = timer_init('    bgcov (ob/grd)')
    timer_postx_mpi = timer_init('    postx mpi', TIMER_SYNC)

  end subroutine solver_init



  !================================================================================
  !================================================================================


  subroutine solver_run(ai)
    real, intent(out), allocatable :: ai(:,:,:)

    integer :: i, j, k, m, n, itr
    integer :: ob1, ob2
    integer :: err, c_error

    ! main variables needed by the conjugate gradient algorithm
    real :: cg_z(obs_num)
    real :: cg_r(obs_num)
    real :: cg_p(obs_num)
    real :: cg_q(obs_num)
    real :: cg_s(obs_num)
    real :: cg_beta, cg_alpha, cg_rs, cg_rs_prev
    
    real :: resid0, resid

    ! local segment of observation blocks
    type local_obs_block_T
       real, allocatable :: l(:)
       logical :: l_valid
       integer :: idx
    end type local_obs_block_T
    integer :: local_obs_block_num
    type(local_obs_block_T), allocatable :: local_obs_block(:)

    ! local segment of analysis grid
    real, allocatable :: local_ai(:,:)

    ! variables needed by the post multiply step
    integer :: obs_points(obs_num)
    real    :: obs_dist(obs_num)
    real    :: r, d
    integer :: num
    real    :: cov(grid_ns)

    call timer_start(timer)

    ! determine the number of observation blocks this proc is responsible for,
    ! and initialize the local obs block information
    local_obs_block_num = count(obs_block_proc == g3dv_mpi_rank)
    call g3dv_mpi_barrier(.true.)
    print ('(A,I4,A,I6,A)'), "  PROC",g3dv_mpi_rank," responsible for",&
         local_obs_block_num," observation blocks"
    call g3dv_mpi_barrier(.true.)
    allocate(local_obs_block(local_obs_block_num))
    n = 1
    do i = 1, size(obs_block_proc)
       if (obs_block_proc(i) .ne. g3dv_mpi_rank) cycle
       local_obs_block(n)%idx = i
       n = n+1
    end do
    ! local_obs_block, now contains an array, the number of obs blocks this proc is to handle
    ! with the %idx field pointing to the block element in the master obs_block arrays
    

    ! ------------------------------------------------------------
    ! cholesky decompositions
    ! ------------------------------------------------------------
    if (isroot) print *, new_line('a'), &
         " Performing Cholesky decomposition of observation blocks...", new_line('a')
    c_error = 0
    call timer_start(timer_cholesky)
    do i = 1, size(local_obs_block)
       ! create a compact array of size n*(n+1)/2
       ! where n = number of observations in this block
       ! in order to store the lower triangular matrix
       j = obs_block_size(local_obs_block(i)%idx)
       j = j*(j+1) / 2
       allocate(local_obs_block(i)%l(j))
       
       ! generate the values of the matrix [HBH^T+R] for this block
       n = obs_block_size(local_obs_block(i)%idx)
       do j = 1, n
          ob1  = obs_block_start(local_obs_block(i)%idx) + j -1
          do k = 1, j
             ob2 = obs_block_start(local_obs_block(i)%idx) + k -1
             m = j + (k-1)*(2*n-k)/2
             local_obs_block(i)%l(m) = bgcov_HBH(obs(ob1), obs(ob2))

             if (j == k) then
                local_obs_block(i)%l(m) = 1.0
                local_obs_block(i)%l(m) = local_obs_block(i)%l(m) + obs(ob1)%err**2
             end if
          end do
       end do
       local_obs_block(i)%l_valid = .true.

       ! calculate the Cholesky decomposition
       call spptrf('L', obs_block_size(local_obs_block(i)%idx), local_obs_block(i)%l, err)
       if (err /= 0) then
          local_obs_block(i)%l_valid = .false.
          c_error = c_error + 1
          print *, "ERROR in cholesky decomposition, (likely not positive definite matrix):",err
       end if
    end do
    c_error = g3dv_mpi_bcstflag(c_error)
    if(isroot .and. c_error > 0) then
       print *, ""
       print *, "ERROR: there were",c_error,"observation block(s) that failed ",&
            "the Cholesky decomposition. Trying to continue by disabling the ",&
            "preconditioner for those block(s), but PCG will likely not converge ",&
            "and results will be bad! Check observations for errors (innsufficient ",&
            "thinning, error values too low, etc.)"
       print *, ""
    end if
    call timer_stop(timer_cholesky)


    ! ------------------------------------------------------------
    ! preconditioned conjugate gradient sover
    !------------------------------------------------------------
    if(isroot) print *, "Running preconditined congjugate gradient iterations..."
    call timer_start(timer_pcg)

    ! initialize
    !------------------------------
    cg_z = 0
    do i = 1, size(cg_r)
       cg_r(i) = obs(i)%inc
    end do


    ! call the preconditioner
    ! TODO, implement this
    cg_s = cg_r

    ! compute r^T * r
    cg_rs = dot_product(cg_r, cg_s)
    resid0 = sqrt(dot_product(cg_r, cg_r))

    if(isroot) print'(A,I4,A,ES12.4)',' iteration',0,'  residual:',resid0

    ! begin the iterations
    !------------------------------
    do itr = 1, maxitr
       ! p_k, \beta
       if (itr == 1) then
          cg_p = cg_s
       else
          cg_beta = cg_rs / cg_rs_prev
          cg_p = cg_s + cg_beta*cg_p
       end if

       ! apply bg cov between each point
       ! ------------------------------
       call timer_start(timer_bgcov_HBH)
       cg_q = 0
       do i =1, size(local_obs_block)
          do j = obs_block_start(local_obs_block(i)%idx), obs_block_end(local_obs_block(i)%idx)
             cg_q(j) = cg_p(j)*obs(j)%err**2

             ! TODO, get proper localization distance
             r = 3000e5
             call obs_search(obs(j)%lat, obs(j)%lon, r,&
                  obs_points, obs_dist, num)
             do k = 1, num
                cg_q(j) = cg_q(j) + &
                     cg_p(obs_points(k))*bgcov_HBH(obs(j), obs(obs_points(k)), obs_dist(j))
             end do
          end do
       end do
       call timer_stop(timer_bgcov_HBH)

       ! mpi broadcast of cg_q
       ! TODO, only share segments this process edited

       !
       cg_alpha = cg_rs / dot_product(cg_p, cg_q)
       cg_z = cg_z + cg_alpha * cg_p
       cg_r = cg_r - cg_alpha * cg_q
       !TODO, call the preconditioner
       cg_s = cg_r
       cg_rs_prev = cg_rs
       cg_rs = dot_product(cg_r, cg_s)
       resid = sqrt(dot_product(cg_r, cg_r))

       if(isroot) print '(A,I4,A,ES12.4)',' iteration',itr,'  residual:',resid
       if ( resid == 0 .or. (conv_ratio > 0 .and. resid0/resid > conv_ratio)) exit

    end do

    call timer_stop(timer_pcg)


    ! ------------------------------------------------------------
    ! post multiply to get the state correction
    ! ------------------------------------------------------------
    if (isroot) then
       print *, ""
       print *, "Applying post multiplication to get AI in state space..."
    end if
    call timer_start(timer_postx)

    ! initialize the local grid analysis increments to 0
    allocate(local_ai(g3dv_mpi_ijcount, grid_ns))
    local_ai = 0

    ! for each gridpoint this proc is to handle...
    do i = 1, g3dv_mpi_ijcount
       ! find all nearby observations 
       ! TODO use actual localization values
       r = 500e3
       call obs_search(grid_local_lat(i), grid_local_lon(i), r, &
            obs_points, obs_dist, num)

       ! for each observation found...
       do j = 1, num
          k = obs_points(j)

          ! calculate the covariance between the observation and each variable / 
          ! level of the grid column
          call timer_start(timer_bgcov_BH)
          cov = bgcov_BH(obs(k), obs_dist(j), i)
          call timer_stop(timer_bgcov_BH)

          ! TODO, more efficient to swap array order?
          ! apply observation effect on grid
          local_ai(i,:) = local_ai(i,:) + cg_z(k) * cov
       end do
    end do


    ! MPI gather of analysis increment
    call timer_start(timer_postx_mpi)
    local_ai = g3dv_mpi_rank
    if(isroot) allocate(ai(grid_nx, grid_ny, grid_ns))
    do i = 1, grid_ns
       call g3dv_mpi_ij2grd_real(local_ai(:,i), ai(:,:,i))
    end do
    call timer_stop(timer_postx_mpi)


    if(isroot) then
       print *, ""
       print *, "**** Sover finished ****"
       print *, ""
    end if

    call timer_stop(timer)
  end subroutine solver_run

end module  g3dv_solver
