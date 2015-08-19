MODULE energy

  USE fldvar
  USE param
  USE param1
  USE rte_do

!   Gas-phase heat of reaction
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  HOR_g

!   Solids-phase heat of reaction
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  HOR_s

!   Gas-solids heat transfer coefficient
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  GAMA_gs

!   Gas-phase radiation coefficient
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  GAMA_Rg

!   Solids-phase radiation coefficient
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  GAMA_Rs

!   Gas-phase radiation temperature
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  T_Rg

!   Solids-phase radiation temperature
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  T_Rs

CONTAINS

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
  !  Source terms for the energy equations.  The source term is linearized as
  !  S = S_c - S_p * T,  where S_c and S_p must be positive.
  !
  !  By default the source terms have been coded for radiation sources.

  !
  !  S_p for gas phase at i, j, k
  DOUBLE PRECISION FUNCTION S_Rpg(IJK)
    IMPLICIT NONE
    INTEGER IJK
    S_Rpg = 16.0d0 * EP_G(IJK) * KAPPA_G * SIGMA_SB * T_g(IJK)**3
  END FUNCTION S_Rpg

  !  S_c for gas phase at i, j, k
  DOUBLE PRECISION FUNCTION S_Rcg(IJK)
    IMPLICIT NONE
    INTEGER IJK
    S_Rcg = S_CONT_G(IJK) &
      + 12.0d0 * EP_G(IJK) * KAPPA_G * SIGMA_SB * T_g(IJK)**4
  END FUNCTION S_Rcg

  !  S_p for solids phase at i, j, k
  DOUBLE PRECISION FUNCTION S_Rps(IJK, M)
    IMPLICIT NONE
    INTEGER IJK, M
    S_Rps = 16.0d0 * EP_S(IJK,M) * KAPPA_S(M) * SIGMA_SB * T_s(IJK,M)**3
  END FUNCTION S_Rps

  !  S_c for solids phase at i, j, k
  DOUBLE PRECISION FUNCTION S_Rcs(IJK, M)
    IMPLICIT NONE
    INTEGER IJK, M
    S_Rcs = S_CONT_S(IJK,M) & 
      + 12.0d0 * EP_S(IJK,M) * KAPPA_S(M) * SIGMA_SB * T_s(IJK,M)**4
  END FUNCTION S_Rcs

END MODULE energy
