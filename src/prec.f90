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
!> @brief  precision settings and floating point comparison
!--------------------------------------------------------------------------------------------------
module prec
  use, intrinsic :: IEEE_arithmetic

  implicit none
  public

  ! https://software.intel.com/en-us/blogs/2017/03/27/doctor-fortran-in-it-takes-all-kinds
  integer,              parameter :: pReal = IEEE_selected_real_kind(15,307)                        !< number with 15 significant digits, up to 1e+-307 (typically 64 bit)

  real(pReal), private, parameter :: PREAL_EPSILON = epsilon(0.0_pReal)                             !< minimum positive number such that 1.0 + EPSILON /= 1.0.
  real(pReal), private, parameter :: PREAL_MIN     = tiny(0.0_pReal)                                !< smallest normalized floating point number

contains


!--------------------------------------------------------------------------------------------------
!> @brief equality comparison for float with double precision
! replaces "==" but for certain (relative) tolerance. Counterpart to dNeq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! AlmostEqualRelative
!--------------------------------------------------------------------------------------------------
logical elemental pure function dEq(a,b,tol)

  real(pReal), intent(in)           :: a,b
  real(pReal), intent(in), optional :: tol
  real(pReal)                       :: eps

  if (present(tol)) then
    eps = tol
  else
    eps = PREAL_EPSILON * maxval(abs([a,b]))
  endif

  dEq = merge(.True.,.False.,abs(a-b) <= eps)

end function dEq


!--------------------------------------------------------------------------------------------------
!> @brief inequality comparison for float with double precision
! replaces "!=" but for certain (relative) tolerance. Counterpart to dEq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! AlmostEqualRelative NOT
!--------------------------------------------------------------------------------------------------
logical elemental pure function dNeq(a,b,tol)

  real(pReal), intent(in)           :: a,b
  real(pReal), intent(in), optional :: tol

  dNeq = .not. dEq(a,b,tol)

end function dNeq


!--------------------------------------------------------------------------------------------------
!> @brief equality to 0 comparison for float with double precision
! replaces "==0" but everything not representable as a normal number is treated as 0. Counterpart to dNeq0
! https://de.mathworks.com/help/matlab/ref/realmin.html
! https://docs.oracle.com/cd/E19957-01/806-3568/ncg_math.html
!--------------------------------------------------------------------------------------------------
logical elemental pure function dEq0(a,tol)

  real(pReal), intent(in)           :: a
  real(pReal), intent(in), optional :: tol
  real(pReal)                       :: eps

  if (present(tol)) then
    eps = tol
  else
    eps = PREAL_MIN * 10.0_pReal
  endif

  dEq0 = merge(.True.,.False.,abs(a) <= eps)

end function dEq0


!--------------------------------------------------------------------------------------------------
!> @brief inequality to 0 comparison for float with double precision
! replaces "!=0" but everything not representable as a normal number is treated as 0. Counterpart to dEq0
! https://de.mathworks.com/help/matlab/ref/realmin.html
! https://docs.oracle.com/cd/E19957-01/806-3568/ncg_math.html
!--------------------------------------------------------------------------------------------------
logical elemental pure function dNeq0(a,tol)

  real(pReal), intent(in)           :: a
  real(pReal), intent(in), optional :: tol

  dNeq0 = .not. dEq0(a,tol)

end function dNeq0


!--------------------------------------------------------------------------------------------------
!> @brief equality comparison for complex with double precision
! replaces "==" but for certain (relative) tolerance. Counterpart to cNeq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! probably a component wise comparison would be more accurate than the comparsion of the absolute
! value
!--------------------------------------------------------------------------------------------------
logical elemental pure function cEq(a,b,tol)

  complex(pReal), intent(in)           :: a,b
  real(pReal),    intent(in), optional :: tol
  real(pReal)                          :: eps

  if (present(tol)) then
    eps = tol
  else
    eps = PREAL_EPSILON * maxval(abs([a,b]))
  endif

  cEq = merge(.True.,.False.,abs(a-b) <= eps)

end function cEq


!--------------------------------------------------------------------------------------------------
!> @brief inequality comparison for complex with double precision
! replaces "!=" but for certain (relative) tolerance. Counterpart to cEq
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
! probably a component wise comparison would be more accurate than the comparsion of the absolute
! value
!--------------------------------------------------------------------------------------------------
logical elemental pure function cNeq(a,b,tol)

  complex(pReal), intent(in)           :: a,b
  real(pReal),    intent(in), optional :: tol

  cNeq = .not. cEq(a,b,tol)

end function cNeq

end module prec
