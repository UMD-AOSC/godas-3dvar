!! author: Travis Sluka
!! category: support
!! Quick sort implementation

module qsort
  !! Basic quicksort implementation.
  !!
  !! Currently only handles a double  list of integers
  !!
  !! @note Algorithm derived from Numerical Recipes, 2007

  implicit none
  private


  ! public module methods
  !------------------------------------------------------------
  public :: qsort_2i



contains



  !================================================================================
  !================================================================================



  subroutine qsort_2i(arr1, arr2)
    !! sort two arrays of integers, using the first array as the sorting key

    integer, intent(inout) :: arr1(:)
    !! Arrays will be sorted according to the values of this array. After calling
    !! the subroutine arr1 will be sorted

    integer, intent(inout) :: arr2(:)
    !! This array will be sorted according to the values of arr1. After calling
    !! the subroutine the respective members of arr1 and arr2 will remain matched together

    integer :: i, j, a, b
    integer, parameter :: M = 7
    integer, parameter :: NSTACK = 64
    integer :: jstack, ir, k, l, istack(NSTACK)

    !------------------------------
    jstack = 0
    l = 1
    ir = size(arr1)
    do while(.true.)
       if(ir-l < M) then
          ! insertion sort if the subarray is small enough
          do j = l+1, ir
             a = arr1(j)
             b = arr2(j)
             do i=j-1,l,-1
                if(arr1(i) < a)  exit
                arr1(i+1) = arr1(i)
                arr2(i+1) = arr2(i)
             end do
             arr1(i+1) = a
             arr2(i+1) = b
          end do
          ! if done sorting the list...
          if(jstack == 0) exit
          ! otherwise pop the stack and begin a new roudn of partitioning
          ir = istack(jstack)
          l  = istack(jstack-1)
          jstack = jstack - 2
       else
          ! quick sort
          k = ishft(l+ir,-1)
          ! choose median of left, center, and right elements as partitioning element
          ! while ensuring those elements are in order
          ! arr1[l] <= arr1[l+1] <= arr1[ir]
          call swap(arr1, arr2, k, l+1)
          if(arr1(l)   > arr1(ir))   call swap(arr1, arr2, l,   ir)
          if(arr1(l+1) > arr1(ir))   call swap(arr1, arr2, l+1, ir)
          if(arr1(l)   > arr1(l+1))  call swap(arr1, arr2, l,   l+1)
          !initialize pointers for for partitioning
          i = l+1
          j = ir
          a = arr1(l+1)  !partitioning element
          b = arr2(l+1)
          do while(.true.)
             j = j-1
             i = i+1
             ! scan up to find element > a
             do while(arr1(i) < a)
                i = i + 1
             end do
             ! scan down to find element < a
             do while(arr1(j) > a)
                j = j - 1
             end do
             if (j < i) exit
             call swap(arr1, arr2, i, j)
          end do
          ! insert patitioning element into array
          arr1(l+1) = arr1(j)
          arr2(l+1) = arr2(j)
          arr1(j) = a
          arr2(j) = b
          ! push pointers to larger subarray on stack; parocess smaller subarray now
          jstack = jstack + 2
          if(jstack >= NSTACK) then
             print *, "ERROR: stack error"
             stop 1
          end if
          if(ir-i+1 >= j-1) then
             istack(jstack) = ir
             istack(jstack-1)=i
             ir = j-1
          else
             istack(jstack)=j-1
             istack(jstack-1)=l
             l=i
          end if
       end if
    end do
  end subroutine qsort_2i



  !================================================================================
  !================================================================================



  pure subroutine swap(a, b, i, j)
    !! Convenience function to swap two elements in two arrays.
    integer, intent(inout) :: a(:), b(:)
    integer, intent(in)    :: i, j
    integer :: t_a, t_b

    t_a = a(i)
    a(i) = a(j)
    a(j) = t_a
    t_b = b(i)
    b(i) = b(j)
    b(j) = t_b
  end subroutine swap



  !================================================================================

end module qsort
