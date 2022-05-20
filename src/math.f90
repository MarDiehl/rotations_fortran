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
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Mathematical library
!--------------------------------------------------------------------------------------------------
module math
  use prec

  implicit none
  real(pReal), parameter :: &
    PI = acos(-1.0_pReal), &                                                                        !< ratio of a circle's circumference to its diameter
    TAU = 2.0_pReal*PI                                                                              !< ratio of a circle's circumference to its radius

contains


!--------------------------------------------------------------------------------------------------
!> @brief trace of a 3x3 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_trace33(m)

  real(pReal), dimension(3,3), intent(in) :: m

  math_trace33 = m(1,1) + m(2,2) + m(3,3)

end function math_trace33


!--------------------------------------------------------------------------------------------------
!> @brief determinant of a 3x3 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_det33(m)

  real(pReal), dimension(3,3), intent(in) :: m

  math_det33 = m(1,1)* (m(2,2)*m(3,3)-m(2,3)*m(3,2)) &
             - m(1,2)* (m(2,1)*m(3,3)-m(2,3)*m(3,1)) &
             + m(1,3)* (m(2,1)*m(3,2)-m(2,2)*m(3,1))

end function math_det33


!--------------------------------------------------------------------------------------------------
!> @brief limits a scalar value to a certain range (either one or two sided)
! Will return NaN if left > right
!--------------------------------------------------------------------------------------------------
real(pReal) pure elemental function math_clip(a, left, right)

  real(pReal), intent(in) :: a
  real(pReal), intent(in), optional :: left, right

  math_clip = a
  if (present(left))  math_clip = max(left,math_clip)
  if (present(right)) math_clip = min(right,math_clip)
  if (present(left) .and. present(right)) &
    math_clip = merge (IEEE_value(1.0_pReal,IEEE_quiet_NaN),math_clip, left>right)

end function math_clip

end module math
