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
    !print*, i
    call quaternion(testData(:,i))
    call matrix(qu2om(testData(:,i)))
    call Eulers(qu2eu(testData(:,i)))
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


!--------------------------------------------------------------------------------------------------
! Rotation matrix forward/backward
!--------------------------------------------------------------------------------------------------
subroutine matrix(om)
  real(pReal), dimension(3,3) :: om

  if(.not. matrix_equal(om,qu2om(om2qu(om)))) stop 'qu2om/om2qu'
  if(.not. matrix_equal(om,eu2om(om2eu(om)))) stop 'eu2om/om2eu'
  if(.not. matrix_equal(om,ax2om(om2ax(om)))) stop 'ax2om/om2ax'
  if(.not. matrix_equal(om,ro2om(om2ro(om)))) stop 'ro2om/om2ro'
  if(.not. matrix_equal(om,ho2om(om2ho(om)))) stop 'ho2om/om2ho'
  if(.not. matrix_equal(om,cu2om(om2cu(om)))) stop 'cu2om/om2cu'

end subroutine matrix

function matrix_equal(om1,om2) result(ok)

  real(pReal), intent(in), dimension(3,3) :: om1,om2
  logical :: ok

  ok = all(dEq(om1,om2,5.0e-8_pReal))
  if(.not. ok) print*, om1,new_line(''),om2

end function matrix_equal


!--------------------------------------------------------------------------------------------------
! Euler angles forward/backward
!--------------------------------------------------------------------------------------------------
subroutine Eulers(eu)
  real(pReal), dimension(3) :: eu

  if(.not. Eulers_equal(eu,qu2eu(eu2qu(eu)))) stop 'qu2eu/eu2qu'
  if(.not. Eulers_equal(eu,om2eu(eu2om(eu)))) stop 'om2eu/eu2om'
  if(.not. Eulers_equal(eu,ax2eu(eu2ax(eu)))) stop 'ax2eu/eu2ax'
  if(.not. Eulers_equal(eu,ro2eu(eu2ro(eu)))) stop 'ro2eu/eu2ro'
  if(.not. Eulers_equal(eu,ho2eu(eu2ho(eu)))) stop 'ho2eu/eu2ho'
  if(.not. Eulers_equal(eu,cu2eu(eu2cu(eu)))) stop 'cu2eu/eu2cu'

end subroutine Eulers

function Eulers_equal(eu1,eu2) result(ok)

  real(pReal), intent(in), dimension(3) :: eu1,eu2
  logical :: ok

  ! low tolerance needed for ho/cu
  ok = all(dEq(eu1,eu2,1.0e-5_pReal))
  if(.not. ok) print*, eu1,new_line(''),eu2

end function Eulers_equal

end program test_rotations
