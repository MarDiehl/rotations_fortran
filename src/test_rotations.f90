! Copyright 2011-20 Max-Planck-Institut für Eisenforschung GmbH
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------------------------------
!> @brief check correctness of rotation conversions
!--------------------------------------------------------------------------------------------------

program test_rotations
  use rotations

  implicit none
  integer, parameter :: N = 50000
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

  real(pReal), dimension(4,N) :: testData
  real(pReal), dimension(size(setOfQuaternions,1),size(setOfQuaternions,2)) :: noise
  integer :: i

  ! test data: special, special+noise, random
  call random_number(testData)
  call random_number(noise)
  do i=0,10
    testData(:,i*size(setOfQuaternions,2)+1:size(SetOfQuaternions,2)*(i+1)) &
      = setOfQuaternions + (noise*2.0_pReal -1.0_pReal) * real(i,pReal)*10._pReal**(-real(i,pReal))
  enddo

  do i = 1, size(testData,2)
    testData(:,i) = testData(:,i)/norm2(testData(:,i))
    if(testData(1,i) < 0.0_pReal) testData(:,i) = -1.0_pReal * testData(:,i)
  enddo

  do i = 1, size(testData,2)
    call quaternion(      testData(:,i))
    call matrix    (qu2om(testData(:,i)))
    call Eulers    (qu2eu(testData(:,i)))
    call axisAngle (qu2ax(testData(:,i)))
    call Rodrigues (qu2ro(testData(:,i)))
  enddo
  print*, 'All fine'

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

  ok = all(dEq(qu1,qu2,1.0e-4_pReal))
  if(dEq0(qu1(1),1.0e-12_pReal)) &
    ok = ok .or. all(dEq(-1.0_pReal*qu1,qu2,1.0e-4_pReal))
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

  ok = all(dEq(om1,om2,5.0e-3_pReal))
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
  real(pReal) sum_phi1,sum_phi2

  ok = all(dEq(eu1,eu2,1.0e-3_pReal))
  if(dEq0(eu1(2),1.0e-3_pReal) .or. dEq(eu1(2),PI,1.0e-3_pReal)) then
    sum_phi1 = mod(eu1(1)+eu1(3)+2.0_pReal*PI,2.0_pReal*PI)
    sum_phi2 = mod(eu2(1)+eu2(3)+2.0_pReal*PI,2.0_pReal*PI)
    ok = dEq(sum_phi1,sum_phi2,5.0e-3_pReal) .or. &
         dEq(sum_phi1+2._pReal*PI,sum_phi2,5.0e-3_pReal) .or. &
         dEq(sum_phi1-2._pReal*PI,sum_phi2,5.0e-3_pReal) .or. &
         dEq(sum_phi1,sum_phi2+2._pReal*PI,5.0e-3_pReal) .or. &
         dEq(sum_phi1,sum_phi2-2._pReal*PI,5.0e-3_pReal)
  endif
  if(.not. ok) print*, eu1,new_line(''),eu2

end function Eulers_equal


!--------------------------------------------------------------------------------------------------
! Axis angle forward/backward
!--------------------------------------------------------------------------------------------------
subroutine axisAngle(ax)
  real(pReal), dimension(4) :: ax

  if(.not. axisAngle_equal(ax,qu2ax(ax2qu(ax)))) stop 'qu2ax/ax2qu'
  if(.not. axisAngle_equal(ax,om2ax(ax2om(ax)))) stop 'om2ax/ax2om'
  if(.not. axisAngle_equal(ax,eu2ax(ax2eu(ax)))) stop 'eu2ax/ax2eu'
  if(.not. axisAngle_equal(ax,ro2ax(ax2ro(ax)))) stop 'ro2ax/ax2ro'
  !ToDo: ho
  !ToDo: cu

end subroutine axisAngle

function axisAngle_equal(ax1,ax2) result(ok)

  real(pReal), intent(in), dimension(4) :: ax1,ax2
  logical :: ok

  ok = all(dEq(ax1,ax2,1.0e-3_pReal))
  if(dEq(ax1(4),PI,5.0e-3_pReal)) ok = ok .or. all(dEq(ax1*real([-1,-1,-1,1],pReal),ax2,1.0e-3_pReal))
  if(dEq(ax2(4),PI,5.0e-3_pReal)) ok = ok .or. all(dEq(ax1*real([-1,-1,-1,1],pReal),ax1,1.0e-3_pReal))
  ok = ok .or. dEq0(ax1(4),1.0e-3_pReal) .or. dEq0(ax2(4),1.0e-3_pReal)
  if(.not. ok) print*, ax1,new_line(''),ax2

end function axisAngle_equal


!--------------------------------------------------------------------------------------------------
! Rodrigues vector forward/backward
!--------------------------------------------------------------------------------------------------
subroutine Rodrigues(ro)
  real(pReal), dimension(4) :: ro
  
  if(.not. Rodrigues_equal(ro,qu2ro(ro2qu(ro)))) stop 'qu2ro/ro2qu'
  if(.not. Rodrigues_equal(ro,om2ro(ro2om(ro)))) stop 'om2ro/ro2om'
  if(.not. Rodrigues_equal(ro,eu2ro(ro2eu(ro)))) stop 'eu2ro/ro2eu'
  if(.not. Rodrigues_equal(ro,ax2ro(ro2ax(ro)))) stop 'ax2ro/ro2ax'

end subroutine Rodrigues

function Rodrigues_equal(ro1,ro2) result(ok)

  real(pReal), intent(in), dimension(4) :: ro1,ro2
  logical :: ok
  real(pReal) :: cutoff = tan(PI*.5_pReal*(1.0-1.e-5_pReal))
  
  ok = all(dEq(math_clip(ro1,right=cutoff),math_clip(ro2,right=cutoff),1.0e-6_pReal)) \
     .or. dEq0(ro1(4),1e-7_pReal)
  if(.not. ok) print*, ro1,new_line(''),ro2

end function Rodrigues_equal


end program test_rotations
