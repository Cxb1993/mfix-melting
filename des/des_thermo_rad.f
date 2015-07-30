!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_RADIATION                                          !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 25-Jun-10  !
!                                                                      !
!  Commen:                                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_RADIATION(I, iM, iIJK, FOCUS)

      USE constant
	  USE geometry
      USE des_thermo
      USE discretelement
      USE fldvar
      USE param1
      USE physprop, only: SMAX
      USE toleranc
      USE compar
      USE functions
      USE indices
      USE run

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
! Global index of particle
      INTEGER, INTENT(IN) :: I
! Solid phase of particle I
      INTEGER, INTENT(IN) :: iM
! Fluid cell index containing particle I
      INTEGER, INTENT(IN) :: iIJK
! Logical used for debugging
      LOGICAL, INTENT(IN) :: FOCUS

! Local variables
!---------------------------------------------------------------------//
! Surface area of particle
      DOUBLE PRECISION :: A_S
! Radiative heat transfer
      DOUBLE PRECISION :: Qrd
! Environment temperature
      DOUBLE PRECISION :: Tenv
! Particle Emmisivity
      DOUBLE PRECISION :: lEm
	  INTEGER :: IJK, J
	  
	  DO J = 2, PART_CELLS(I,1) + 1 
		IJK = PART_CELLS(I,J)
		IF(.NOT.FLUID_AT(IJK)) CYCLE
        IF(PART_VOL_INTERSEC(IJK,I) == 0.0d0) CYCLE
        Q_SOURCE(I) = Q_SOURCE(I) + &
			PART_VOL_INTERSEC(IJK,I)/TOT_VOL_INTERSEC(IJK)*&
			S_RC_DES(IJK)*VOL(IJK)
      END DO

      RETURN

      END SUBROUTINE DES_RADIATION
	  
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: DES_RAD                                                 !
!                                                                      !
!  Purpose: This routine is called from the continuum phase and        !
!  calculates the source term from the P-1 Radiation model             !
!                                                                      !
!  Author: J.Musser                                   Date: 15-Jan-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE  DES_RAD(S_C,M)

      USE compar
      Use constant
      Use des_thermo
      Use discretelement
      Use fldvar
      USE geometry
      USE indices
      Use interpolation
      Use param1
      Use physprop

      use run, only: ODT
      use functions

      IMPLICIT NONE

! Passed Variables
!---------------------------------------------------------------------//
! Source term on RHS
      DOUBLE PRECISION, INTENT(INOUT) :: S_C(DIMENSION_3)
! Solids phase
	  INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! IJK value of cell containing particle NP
      INTEGER :: IJK
!---------------------------------------------------------------------//

! Loop over fluid cells.
      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE IJK_LP
! Redistribute the energy over the fluid time step. Note that by the
! time this routine is called, S_C and S_P have already been multiplied
! by the fluid cell volume. Thus, the mapping should result in units
! of energy per time.
         S_C(IJK) = S_C(IJK) + DES_ENERGY_SOURCE_S(IJK,M)*ODT
      ENDDO IJK_LP ! End loop over fluid cells

      RETURN
      END SUBROUTINE  DES_RAD
