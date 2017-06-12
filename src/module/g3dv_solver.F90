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
  use mpi


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
  integer :: timer_obsearch_HBH  
  integer :: timer_bgcov_HBH
  integer :: timer_pcg_mpi
  integer :: timer_postx
  integer :: timer_bgcov_BH
  integer :: timer_postx_mpi
  integer :: timer_obsearch_BH
  
  ! variables read in from namelist
  integer :: maxitr = 10
  real    :: conv_ratio = 1e2

  !#define SOLVER_PREC_SINGLE
#ifdef SOLVER_PREC_SINGLE
  integer, parameter :: r_p = kind(1.e0)
  integer, parameter :: r_p_mpi = mpi_real
#else
  integer, parameter :: r_p = kind(1.d0)
  integer, parameter :: r_p_mpi = mpi_double
#endif

  ! other output
  real, public, protected, allocatable :: solver_local_maxobcor(:,:)

  
  type local_obs_block_T
     real(r_p), allocatable :: l(:)
     logical :: l_valid
     integer :: idx
  end type local_obs_block_T




    
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
            new_line('a'), "solver_init() : preconditioned congjugate gradient solver",&
            new_line('a'), "------------------------------------------------------------"
    end if
    
    

    ! read in our section of the namelist
    open(newunit=unit, file=nml)
    read(unit, nml=g3dv_solver)
    close(unit)
    if (isroot) print g3dv_solver

    if(isroot) then
       print *, ""
#ifdef SOLVER_PREC_SINGLE
       print *, "  Running in SINGLE precision mode."
#else
       print *, "  Running in DOUBLE precision mode."
#endif
    end if

    
    !setup timer
    timer              = timer_init('(Solver)', TIMER_SYNC)
    timer_cholesky     = timer_init('  cholesky decomp',TIMER_SYNC)
    timer_pcg          = timer_init('  PCG',TIMER_SYNC)
    timer_precond      = timer_init('    preconditioner')
    timer_obsearch_HBH = timer_init('    obsearch_HBH')    
    timer_bgcov_HBH    = timer_init('    bgcov (ob/ob)')
    timer_pcg_mpi      = timer_init('    pcg mpi', TIMER_SYNC)
    timer_postx        = timer_init('  post multiply', TIMER_SYNC)
    timer_obsearch_BH  = timer_init('    obsearch_BH')
    timer_bgcov_BH     = timer_init('    bgcov (ob/grd)')
    timer_postx_mpi    = timer_init('    postx mpi', TIMER_SYNC)

  end subroutine solver_init



  !================================================================================
  !================================================================================


  subroutine solver_run(local_ai)
    !    real, intent(out), allocatable :: ai(:,:,:)
    real, intent(inout) :: local_ai(grid_ns, g3dv_mpi_ijcount)
    
   
    integer :: i, j, k, m, n, itr
    integer :: ob1, ob2
    integer :: err, c_error
    real :: prev_lon, prev_lat
    
    ! main variables needed by the conjugate gradient algorithm
    real(r_p) :: cg_z(obs_num)
    real(r_p) :: cg_r(obs_num)
    real(r_p) :: cg_p(obs_num)
    real(r_p) :: cg_q(obs_num)
    real(r_p) :: cg_s(obs_num)
    real(r_p) :: cg_beta, cg_alpha, cg_rs, cg_rs_prev
    real(r_p) :: cg_q_j
    real(r_p) :: resid0, resid

    ! local segment of observation blocks
    integer :: local_obs_block_num
    type(local_obs_block_T), allocatable :: local_obs_block(:)

    ! variables needed by the post multiply step
    integer :: obs_points(obs_num)
    real    :: obs_dist(obs_num)
    real    :: r
    integer :: num
    real    :: cov(grid_ns), cor(grid_ns), var(grid_ns), cor_max(grid_ns)

    call timer_start(timer)
    
    allocate(solver_local_maxobcor(grid_ns, g3dv_mpi_ijcount))
    solver_local_maxobcor = 0.0

    ! determine the number of observation blocks this proc is responsible for,
    ! and initialize the local obs block information
    local_obs_block_num = count(obs_block_proc == g3dv_mpi_rank)
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
             
             if(ob1 == ob2)  then
                local_obs_block(i)%l(m) = 1.0 + obs(ob1)%err**2                
             else
                local_obs_block(i)%l(m) = bgcov_HBH(obs(ob1), obs(ob2))
             end if
          end do
       end do
       local_obs_block(i)%l_valid = .true.

       ! calculate the Cholesky decomposition
#ifdef SOLVER_PREC_SINGLE       
       call spptrf('L', obs_block_size(local_obs_block(i)%idx), local_obs_block(i)%l, err)
#else
       call dpptrf('L', obs_block_size(local_obs_block(i)%idx), local_obs_block(i)%l, err)
#endif
       
       if (err /= 0) then
          local_obs_block(i)%l_valid = .false.
          c_error = c_error + 1
          print *, "ERROR in cholesky decomposition, (likely not positive definite matrix):",err
          n = obs_block_size(local_obs_block(i)%idx)
          do j = 1, n
             ob1  = obs_block_start(local_obs_block(i)%idx) + j -1
             print *, j, obs(ob1)%id, obs(ob1)%lat, obs(ob1)%lon, obs(ob1)%dpth, obs(ob1)%grd_w%y(1), obs(ob1)%grd_vtloc
          end do
          stop 1          
       end if
    end do
    call timer_stop(timer_cholesky)
    
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




    ! ------------------------------------------------------------
    ! preconditioned conjugate gradient sover
    !------------------------------------------------------------
    if(isroot) print *, "Running preconditined congjugate gradient iterations..."
    call timer_start(timer_pcg)

    ! initialize
    !------------------------------
    cg_z = 0.0
    do i = 1, size(cg_r)
       cg_r(i) = obs(i)%inc
    end do

    ! call the preconditioner
    call precondition(local_obs_block, cg_r, cg_s)

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
       cg_q = 0.0
       prev_lat = 1e30
       prev_lon = 1e30
       do i = 1, size(obs_local_idx)
          j = obs_local_idx(i)

          if(prev_lat /= obs(j)%lat .or. prev_lon /= obs(j)%lon) then
             !TODO, ensure profile obs are kept together, in order to
             ! accelerate this segment
             ! find nearby points
             call timer_start(timer_obsearch_HBH)
             r = bgcov_hzdist(obs(j)%lat)*1.1 * 2.0/sqrt(0.3) ! expanded by 10% to account for point
              ! with small radius in range of point with larger radius
             call obs_search(obs(j)%lat, obs(j)%lon, r,&
                  obs_points, obs_dist, num)
             call timer_stop(timer_obsearch_HBH)
             prev_lon = obs(j)%lon
             prev_lat = obs(j)%lat
          end if
          
          ! calculate the covariance with the surrounding points
          !------------------------------
          call timer_start(timer_bgcov_HBH)

          ! obs error ( R )
          cg_q_j = cg_p(j)*obs(j)%err**2

          ! covariance ( HBH )
          do k = 1, num                 
             cg_q_j = cg_q_j + cg_p(obs_points(k)) * &
                  bgcov_HBH(obs(j), obs(obs_points(k)), obs_dist(k))
          end do
          call timer_stop(timer_bgcov_HBH)
          cg_q(j) = cg_q_j
       end do



       ! mpi broadcast of cg_q
       ! TODO, move this to g3dv_mpi module
       ! TODO, do this more efficiently
       call timer_start(timer_pcg_mpi)
       call mpi_allreduce(mpi_in_place, cg_q, size(cg_q), r_p_mpi, mpi_sum, g3dv_mpi_comm, i)
       call timer_stop(timer_pcg_mpi)

       ! calculate other things
       cg_alpha = cg_rs / dot_product(cg_p, cg_q)
       cg_z = cg_z + cg_alpha * cg_p
       cg_r = cg_r - cg_alpha * cg_q
       call precondition(local_obs_block, cg_r, cg_s)
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
    local_ai = 0

    
    ! for each gridpoint this proc is to handle...
    do i = 1, g3dv_mpi_ijcount

       if (grid_local_mask(i) <= 0) cycle

       ! find all nearby observations 
       ! TODO use actual localization values
       call timer_start(timer_obsearch_BH)       
       r = bgcov_hzdist(grid_local_lat(i)) * 2.0/sqrt(0.3)
       call obs_search(grid_local_lat(i), grid_local_lon(i), r, &
            obs_points, obs_dist, num)
       call timer_stop(timer_obsearch_BH)


       
       ! for each observation found...
       cov = 0.0
       cor_max = 0.0
       call timer_start(timer_bgcov_BH)       
       do j = 1, num
          k = obs_points(j)

          ! calculate the covariance between the observation and each variable / 
          ! level of the grid column
          call bgcov_BH(obs(k), obs_dist(j), i, cor, var)
          cov = cov + cg_z(k)*cor*var
          do k = 1, grid_ns
             cor_max(k) = max(cor_max(k), cor(k))
          end do
       end do
       call timer_stop(timer_bgcov_BH)
       ! apply observation effect on grid       
       local_ai(:,i) = cov
       solver_local_maxobcor(:,i) = cor_max

    end do

    call timer_stop(timer_postx)

    
    !------------------------------
    if(isroot) then
       print *, ""
       print *, "**** Sover finished ****"
       print *, ""
    end if

    call timer_stop(timer)
  end subroutine solver_run


  !================================================================================
  !================================================================================


  
  subroutine precondition(blocks, r,s)
    real(r_p), intent(in)  :: r(:)
    real(r_p), intent(out) :: s(:)
    type(local_obs_block_T), intent(in) :: blocks(:)

    integer :: i, err, vs, ve

    call timer_start(timer_precond)
    s = 0
    do i = 1, size(blocks)
       vs = obs_block_start(blocks(i)%idx)
       ve = obs_block_end(  blocks(i)%idx)
       s(vs:ve) = r(vs:ve)          
       if(blocks(i)%l_valid) then
#ifdef SOLVER_PREC_SINGLE          
          call spptrs('L', obs_block_size(blocks(i)%idx), 1, blocks(i)%L, &
               s(vs:ve), obs_block_size(blocks(i)%idx), err)
#else
          call dpptrs('L', obs_block_size(blocks(i)%idx), 1, blocks(i)%L, &
               s(vs:ve), obs_block_size(blocks(i)%idx), err)
#endif        
          if(err/=0)then
             print *, "ERROR applying preconditioner"
             stop 1
          end if
       end if
    end do

    ! mpi broadcast of cg_q
    ! TODO, move this to g3dv_mpi module
    ! TODO, do this more efficiently
    call timer_start(timer_pcg_mpi)           
    call mpi_allreduce(mpi_in_place, s, size(s), r_p_mpi, mpi_sum, g3dv_mpi_comm, i)
    call timer_stop(timer_pcg_mpi)
    
    call timer_stop(timer_precond)

  end subroutine precondition



  !================================================================================
  !================================================================================



end module  g3dv_solver
