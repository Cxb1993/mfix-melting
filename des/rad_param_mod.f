!!*****************************************************************
!!** module rad_param
!!*****************************************************************
      module rad_param
! Parameters will be used to control radiation calculations.
! Occationally they may need to be adjusted.
! They could be converted to keywords if necessary

      real(8), parameter :: BETAMIN = 1.d-5
! minimum beta (extinction coefficient) for radiation calculation
! If the maximum extinction coefficient is below this number,
! the problem is considered as optically thin. Then absorption is
! negligible, and radiative heat source comes from emission only.

      real(8), parameter :: BETALB = 1.d-6
! the lower bound value of extinction coefficient. The extinction
! coefficient in the flow field will be lower bounded by this
! number. This is because 1/Beta is used in P1 calculation and
! a small number of beta causes numerical instability

      logical, parameter :: isSootGray=.true.  
! logical constant for gray soot in spectral calculation
      end module