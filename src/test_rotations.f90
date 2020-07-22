program test_rotations
  use rotations

  implicit none
  real(pReal), dimension(4,52) :: setOfQuaternions = &
    reshape([1.0, 0.0, 0.0, 0.0, &
             !------------------
             0.0, 1.0, 0.0, 0.0, &
             0.0, 0.0, 1.0, 0.0, &
             0.0, 0.0, 0.0, 1.0, &
             0.0,-1.0, 0.0, 0.0, &
             0.0, 0.0,-1.0, 0.0, &
             0.0, 0.0, 0.0,-1.0, &
             !------------------
             1.0, 1.0, 0.0, 0.0, &
             1.0, 0.0, 1.0, 0.0, &
             1.0, 0.0, 0.0, 1.0, &
             0.0, 1.0, 1.0, 0.0, &
             0.0, 1.0, 0.0, 1.0, &
             0.0, 0.0, 1.0, 1.0, &
             !------------------
             1.0,-1.0, 0.0, 0.0, &
             1.0, 0.0,-1.0, 0.0, &
             1.0, 0.0, 0.0,-1.0, &
             0.0, 1.0,-1.0, 0.0, &
             0.0, 1.0, 0.0,-1.0, &
             0.0, 0.0, 1.0,-1.0, &
             !------------------
             0.0, 1.0,-1.0, 0.0, &
             0.0, 1.0, 0.0,-1.0, &
             0.0, 0.0, 1.0,-1.0, &
             !------------------
             0.0,-1.0,-1.0, 0.0, &
             0.0,-1.0, 0.0,-1.0, &
             0.0, 0.0,-1.0,-1.0, &
             !------------------
             1.0, 1.0, 1.0, 0.0, &
             1.0, 1.0, 0.0, 1.0, &
             1.0, 0.0, 1.0, 1.0, &
             1.0,-1.0, 1.0, 0.0, &
             1.0,-1.0, 0.0, 1.0, &
             1.0, 0.0,-1.0, 1.0, &
             1.0, 1.0,-1.0, 0.0, &
             1.0, 1.0, 0.0,-1.0, &
             1.0, 0.0, 1.0,-1.0, &
             1.0,-1.0,-1.0, 0.0, &
             1.0,-1.0, 0.0,-1.0, &
             1.0, 0.0,-1.0,-1.0, &
             !------------------
             0.0, 1.0, 1.0, 1.0, &
             0.0, 1.0,-1.0, 1.0, &
             0.0, 1.0, 1.0,-1.0, &
             0.0,-1.0, 1.0, 1.0, &
             0.0,-1.0,-1.0, 1.0, &
             0.0,-1.0, 1.0,-1.0, &
             0.0,-1.0,-1.0,-1.0, &
             !------------------
             1.0, 1.0, 1.0, 1.0, &
             1.0,-1.0, 1.0, 1.0, &
             1.0, 1.0,-1.0, 1.0, &
             1.0, 1.0, 1.0,-1.0, &
             1.0,-1.0,-1.0, 1.0, &
             1.0,-1.0, 1.0,-1.0, &
             1.0, 1.0,-1.0,-1.0, &
             1.0,-1.0,-1.0,-1.0],[4,52])
  
  real(pReal), dimension(:,:), allocatable :: testData
  integer :: i,N
  
  N = size(setOfQuaternions,2)
  allocate(testData(4,N))
  testData(:,:size(setOfQuaternions,2)) = setOfQuaternions

  do i = 1, size(testData,2)
    testData(:,i) = testData(:,i)/norm2(testData(:,i))
  enddo

  do i = 1, size(testData,2)
    print*, i
    call quaternion(testData(:,i))
  enddo

contains 

!--------------------------------------------------------------------------------------------------
! Quaternion forward/backward
!--------------------------------------------------------------------------------------------------
subroutine quaternion(qu)
  real(pReal), dimension(4) :: qu
  
  if(.not. quaternion_equal(qu,om2qu(qu2om(qu)))) stop 'om2qu/qu2om'

end subroutine quaternion

function quaternion_equal(qu1,qu2) result(ok)
  
  real(pReal), intent(in), dimension(4) :: qu1,qu2
  logical :: ok

  ok = all(dEq(qu1,qu2,1.0e-8_pReal))
  if(dEq0(qu1(1),1.0e-12_pReal)) &
    ok = ok .or. all(dEq(-1.0_pReal*qu1,qu2,1.0e-8_pReal))
  if(.not. ok) print*, qu1,new_line(''),qu2

end function quaternion_equal



end program test_rotations
