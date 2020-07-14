!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Mathematical library, including random number generation and tensor representations
!--------------------------------------------------------------------------------------------------
module math
  use prec

  implicit none
  real(pReal),    parameter :: PI = acos(-1.0_pReal)                                                !< ratio of a circle's circumference to its diameter

contains
!--------------------------------------------------------------------------------------------------
!> @brief trace of a 3x3 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_trace33(m)

  real(pReal), dimension(3,3), intent(in) :: m

  math_trace33 = m(1,1) + m(2,2) + m(3,3)

end function math_trace33


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
