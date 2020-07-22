program test_rotations
  use rotations

  implicit none
  integer, parameter :: N = 1100
  real(pReal), parameter :: scatter = 1.e-2_pReal
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
  real(pReal), dimension(size(setOfQuaternions,1),size(setOfQuaternions,2)) :: noise
  integer :: i 
  
  allocate(testData(4,N))
  call random_number(testData)
  testData(:,:size(setOfQuaternions,2)) = setOfQuaternions
  call random_number(noise)
  testData(:,size(setOfQuaternions,2)+1:size(SetOfQuaternions,2)*2) = setOfQuaternions &
                                                                    + (noise*2.0_pReal -1.0_pReal) * scatter

  do i = 1, size(testData,2)
    testData(:,i) = testData(:,i)/norm2(testData(:,i))
    if(testData(1,i) < 0.0_pReal) testData(:,i) = -1.0_pReal * testData(:,i)
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
  if(.not. quaternion_equal(qu,eu2qu(qu2eu(qu)))) stop 'eu2qu/qu2eu'
  if(.not. quaternion_equal(qu,ax2qu(qu2ax(qu)))) stop 'ax2qu/qu2ax'
  if(.not. quaternion_equal(qu,ro2qu(qu2ro(qu)))) stop 'ro2qu/qu2ro'
  if(.not. quaternion_equal(qu,ho2qu(qu2ho(qu)))) stop 'ho2qu/qu2ho'
  if(.not. quaternion_equal(qu,cu2qu(qu2cu(qu)))) stop 'cu2qu/qu2cu'

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
