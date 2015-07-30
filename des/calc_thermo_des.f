!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_THERMO_DES                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_THERMO_DES

      use physprop, only: SMAX
      use physprop, only: K_s0

      USE compar
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE interpolation
      USE param1
      USE run

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! Index of neighbor particle of particle I such that I < J
      INTEGER IJK
! Loop index for particles.
      INTEGER NP, lNP
! Phase index for particle NP
      INTEGER M
! Identifies that the indicated particle is of interest for debugging
      LOGICAL FOCUS
! Variables needed for calculating new interpolation quantities for
! species and energy equations
      INTEGER INTERP_IJK(2**3)
      DOUBLE PRECISION INTERP_WEIGHTS(2**3)
	  DOUBLE PRECISION S_RC_CONT(DIMENSION_3)
	  DOUBLE PRECISION ALPHA_TOT
	  INTEGER J, pIJK

! Functions
!---------------------------------------------------------------------//

! This is a quick work-around to keep the thermo routines from causes
! issues while the "check_data" routines are rewritten. Moving forward
! this routine should be split apart to avoid the particle loops for
! cold-flow, non-reacting cases.
      IF(.NOT.ENERGY_EQ .AND. .NOT.ANY_SPECIES_EQ) RETURN

      IF(ANY(CALC_RADT_DES(:))) THEN 
         S_RC_DES(:) = ZERO
		 S_RC_CONT(:) = ZERO
         CALL CALC_RAD_HEAT_SRC_DES(S_RC_DES,S_RC_CONT)
      END IF	  
	  
! Loop over fluid cells.
!---------------------------------------------------------------------//
      IJK_LP: DO IJK = IJKSTART3, IJKEND3	  
         IF(.NOT.FLUID_AT(IJK)) CYCLE IJK_LP
		 
		 IF(ANY(CALC_RADT_DES(:))) THEN
			ALPHA_TOT = 0.0d0
			DO M = 1, SMAX
				ALPHA_TOT = ALPHA_TOT + ALPHA_S(M)*EP_S(IJK,M)
			END DO
			DO M = 1, SMAX
				IF(EP_S(IJK,M) > 0.0d0) THEN 
					DES_ENERGY_SOURCE_S(IJK,M) = DES_ENERGY_SOURCE_S(IJK,M) + &
						S_RC_CONT(IJK)*VOL(IJK)*ALPHA_S(M)*EP_S(IJK,M) / &
						ALPHA_TOT*DTSOLID
				END IF
			END DO
		 END IF
		 
         IF(PINC(IJK) == 0) CYCLE IJK_LP

! Interpolation: Removed J.Musser 11/8/2012
!---------------------------------------------------------------------//
!     IF(DES_INTERP_ON .AND. (ANY_SPECIES_EQ .OR. DES_CONV_EQ)) THEN
!         INTERP_IJK(:) = -1
!         INTERP_WEIGHTS(:) = ZERO
!         CALL INTERPOLATE_CC(NP, INTERP_IJK, INTERP_WEIGHTS, FOCUS)
!      ENDIF

! Preform user-defined calculations from fluid grid.
         IF(CALL_USR) CALL USR4_DES(IJK)

! Loop over all particles in cell IJK.
!---------------------------------------------------------------------//
         lNP_LP: DO lNP = 1, PINC(IJK)
            NP = PIC(IJK)%p(lNP)

! Skip indices that do not represent particles
            IF(.NOT.PEA(NP,1)) CYCLE lNP_LP
            IF(any(PEA(NP,2:4))) CYCLE lNP_LP

! Reset the debug flag
            FOCUS = .FALSE.

! Calculate time dependent physical properties
            CALL DES_PHYSICAL_PROP(NP, FOCUS)

! Identify the solid phases of each particle
            M = PIJK(NP,5) + SMAX

! calculate heat transfer via convection
            IF(CALC_CONV_DES) CALL DES_CONVECTION(NP, M, IJK, &
               INTERP_IJK, INTERP_WEIGHTS, FOCUS)

! calculate heat transfer via radiation
            IF(CALC_RADT_DES(M)) CALL DES_RADIATION(NP, M, IJK, FOCUS)
			
!---Begin Dan Moser Changes------------------------------!
!if particle will go over melting temp, set temp to melting 
!temp and decrement source accordingly
            IF(ANY_SPECIES_EQ .AND. ENERGY_EQ) THEN
               IF(DES_T_S_NEW(NP) + DTSOLID*(Q_Source(NP) / &
               (PMASS(NP) * DES_C_ps(NP))) > 463.0d0) THEN
                  Q_Source(NP)=Q_Source(NP)-(463.0d0-DES_T_S_NEW(NP))/&
                  DTSOLID*PMASS(NP)*DES_C_ps(NP)
                  !DES_T_S_OLD(NP) = 463.0d0
                  DES_T_S_NEW(NP) = 463.0d0
               END IF
            END IF
!---End Dan Moser Changes--------------------------------!

! Calculate reaction rates and interphase mass transfer
            IF(ANY_SPECIES_EQ) THEN
				DO J = 2, PART_CELLS(NP,1) + 1
					pIJK = PART_CELLS(NP,J)
					CALL DES_RRATES0(NP, M, pIJK, &
						INTERP_IJK, INTERP_WEIGHTS, FOCUS)
				END DO
			END IF
			
         ENDDO lNP_LP ! End loop over all particles
      ENDDO IJK_LP ! End loop over fluid cells

      END SUBROUTINE CALC_THERMO_DES
	  
	  
!------------------P-1 Radiation Model Functions/Modules--------!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
      subroutine calc_rad_heat_src_des(RAD_SRC,RAD_SRC_CONT)
!----------------------------------------------------------------
! Modules
!----------------------------------------------------------------
      use rte_p1, only : solveP1Eqn=>solve_p1_eqn
      Use constant
      Use des_thermo
      Use discretelement
      Use bc
      Use fldvar
      USE geometry
      USE indices
      Use interpolation
      Use param1
      Use run
      USE compar 
      
!----------------------------------------------------------------
! Passed variables
!----------------------------------------------------------------
      real(8), intent(out) :: RAD_SRC(DIMENSION_3)
	  real(8), intent(out) :: RAD_SRC_CONT(DIMENSION_3)

!---------------------------------------------------------------- 
! Local variables 
!----------------------------------------------------------------
! Indices
      INTEGER IJK, I, J, K, L
! Inverse of total extinction coefficient
      real(8) :: beta(DIMENSION_3) ! beta is extinction coef
      integer :: erflag
      real(8) :: G(DIMENSION_3) ! incident radiation field
      real(8) :: dQ_g(DIMENSION_3) ! del.q_g
      real(8) :: WallTemp(DIMENSION_BC)
      real(8) :: WallEmis(DIMENSION_BC)
      real(8) :: EmisW(DIMENSION_BC) ! wall emission
      real(8) :: EMIS_DES(DIMENSION_3)
	  real(8) :: EMIS_CONT(DIMENSION_3)
	  real(8) :: K_CONT(DIMENSION_3)
      real(8) :: K_DES(DIMENSION_3)
	  real(8) :: K_TOT(DIMENSION_3)
      real(8) :: SCAT_DES(DIMENSION_3)
      real(8) :: Emis(DIMENSION_3)
      real(8) :: omega(DIMENSION_3), w(DIMENSION_3)
      real(8) :: S_RC(DIMENSION_3) ! TOTAL heat source
      real start,finish
      real(8) :: coord
      integer :: step
      real(8) :: G_WALL
!----------------------------------------------------------------
!----------------------------------------------------------------
!     print*, "Radiation Calc"
!     CALL CPU_TIME(start)
      call calc_des_rad_prop(K_DES,K_CONT,EMIS_DES,EMIS_CONT,SCAT_DES)
!     CALL CPU_TIME(finish)
!     print *,'calc_des_rad_prop',finish-start
!     solve p1 for each g values
      CALL ExtrapolateGhostCells(K_DES)
	  CALL ExtrapolateGhostCells(K_CONT)
      CALL ExtrapolateGhostCells(SCAT_DES)
      K_TOT(:)=max(K_DES(:)+K_CONT(:), 1.d-9) ! bound k
      BETA(:)=SCAT_DES(:)+K_TOT(:)
      OMEGA(:) = SCAT_DES(:)/BETA(:)
      W(:) = 3.0d0*(1.d0-OMEGA(:)) ! isotropic scattering 
      EMIS(:) = 4.0d0*pi*(EMIS_DES(:)+EMIS_CONT(:))    
!     calculate wall emission. Using BC_TW_S for the first solid
!     BC for all solids for now. Can improve on this later
      !DO I=1,DIMENSION_BC
      !   IF (.NOT.BC_DEFINED(I)) THEN
      !      WallTemp(I) = 0.0d0
      !      WallEmis(I) = 0.0d0
      !   ELSE IF (BC_TW_S(I,2) == UNDEFINED) THEN
      !      WallTemp(I) = 0.0d0
      !      WallEmis(I) = 0.0d0
      !   ELSE
      !      WallTemp(I) = BC_TW_S(I,2)
! Assuming walls have emissivity of 1 for now. Can improve this later
      !      WallEmis(I) = 1.0d0
      !   END IF
      !END DO
	  WallEmis(:) = 1.0d0
	  WallTemp(:) = 0.0d0
	  WallTemp(4) = 5000.0d0
      EMISW(:) = 4.0d0 * SB_CONST * WallTemp(:)**4
!     CALL CPU_TIME(start)
      call SolveP1Eqn(dq_g, G, beta, w,omega, Emis, WallEmis, Emisw)
!     CALL CPU_TIME(finish)
!     print *,'p1 solve',finish-start
      RAD_SRC(:) = K_DES(:)*G(:) - 4.0d0*pi*EMIS_DES(:)
	  RAD_SRC_CONT(:) = K_CONT(:)*G(:) - 4.0d0*pi*EMIS_CONT(:)

      end subroutine calc_rad_heat_src_des

      subroutine calc_des_rad_prop(K_DES,K_CONT,EMIS_DES,EMIS_CONT,SCAT_DES)

      USE constant
      Use des_thermo
      Use discretelement
      Use fldvar
      USE geometry
      USE indices
      Use interpolation
      Use param1
      Use run
      USE compar
	  use functions
	  use physprop, only: SMAX
      
      IMPLICIT NONE
!----------------------------------------------------------------
! Passed variables
!----------------------------------------------------------------
      real(8),intent(out),DIMENSION(DIMENSION_3) :: K_DES,EMIS_DES,&
      SCAT_DES,K_CONT,EMIS_CONT

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle index
      INTEGER L
! accounted for particles
      INTEGER PC
! solids phase index
      INTEGER M 
! ijk indices      
      INTEGER IJK, I, J, K
! 1 over volume of fluid cell      
      REAL(8) :: OVOL
! percent of particle in a cell and projected area of particle
      REAL(8) :: P_IN_CELL, A_P
! interim variables for testing full volume averaging calc
!      REAL(8) :: K,EMIS,SCAT
!----------------------------------------------
	

      K_DES(:) = ZERO
      EMIS_DES(:) = ZERO
      SCAT_DES(:) = ZERO
	  K_CONT(:) = ZERO
	  EMIS_CONT(:) = ZERO

      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         OVOL = 1.0d0/VOL(IJK)
         DO L = 1, MAX_PIP
            IF(.NOT.PEA(L,1)) CYCLE
            P_IN_CELL = PART_VOL_INTERSEC(IJK,L)/(4.0d0/3.0d0*Pi*&
            DES_RADIUS(L)**3)
            M = PIJK(L,5) + SMAX
            A_P = Pi*DES_RADIUS(L)**2
            K_DES(IJK) = K_DES(IJK) + DES_EM(M)*A_P*OVOL*P_IN_CELL
            EMIS_DES(IJK) = EMIS_DES(IJK) + DES_EM(M)*A_P*&
                SB_CONST*DES_T_S_NEW(L)**4/Pi*OVOL*P_IN_CELL
            SCAT_DES(IJK) = SCAT_DES(IJK)+(1.0d0-DES_EM(M))*A_P*&
                OVOL*P_IN_CELL
         END DO
		 DO M = 1, SMAX
			K_CONT(IJK) = K_CONT(IJK) + ALPHA_S(M)*EP_S(IJK,M)
			EMIS_CONT(IJK) = EMIS_CONT(IJK) + ALPHA_S(M)*SB_CONST*&
				T_S(IJK,M)**4/Pi*EP_S(IJK,M)
		 END DO
      END DO
!      K = 0
!      EMIS = 0
!      SCAT = 0
!      DO L = 1, MAX_PIP
!         IF(.NOT.PEA(L,1)) CYCLE
!         M = PIJK(L,5)
!         A_P = Pi*DES_RADIUS(L)**2
!         K = K + DES_EM(M)*A_P
!         EMIS = EMIS + DES_EM(M)*A_P*SB_CONST*DES_T_S_NEW(L)**4/Pi
!         SCAT = SCAT + (1.0d0-DES_EM(M))*A_P
!      END DO
!      OVOL = 1/(XLENGTH*YLENGTH*ZLENGTH)
!      K_DES(:) = K*OVOL
!      EMIS_DES(:) = EMIS*OVOL
!      SCAT_DES(:) = SCAT*OVOL

      end subroutine calc_des_rad_prop

      subroutine ExtrapolateGhostCells(Phi)
      Use fldvar
      USE geometry
      USE indices
      Use interpolation
      Use run
      USE compar
      use param
      use fldvar
	  use functions
	  
	  IMPLICIT NONE
      real(8), dimension(DIMENSION_3), intent(inout) :: Phi
      
      integer :: ijk            ! loop variable over cells
      integer :: M		! loop variable over solid species
!     ideally i should have dependency on species and temperature but will put them back later
! try uniform here first
      integer :: IMJK, IPJK, IJMK, IJPK, IJKM, IJKP ! neighbor indices
      integer :: IJK_1, IJK_2	! two neighbor indices for linear extrapolation
      integer :: ind_1, ind_2, ind_0 ! directional indices of two real cells and ghost cell
! double precision :: coor_x(0:IMAX2), coor_y(0:JMAX2), coor_z(0:KMAX2) ! coordinate
! double precision :: YMIN, ZMIN
double precision :: dl_1, dl_2, dl_0 ! cell sizes at boundary for linear extrapolation
                                  ! dl_0 is ghost cell, dl_1 next to ghost, dl_2 next to dl_1

! calculate gas phase absorption coefficients

      do IJK = 1, DIMENSION_3   ! loop to search ghost cells
         if (.NOT.FLUID_AT(IJK)) then ! ghost cells
            if (DO_I) then
               IMJK = IM_OF(IJK)
               IPJK = IP_OF(IJK)
            else 
               IMJK = IJK
               IPJK = IJK
            end if

            if (DO_J) then
               IJMK = JM_OF(IJK)
               IJPK = JP_OF(IJK)
            else 
               IJMK = IJK
               IJPK = IJK
            end if

            if (DO_K) then
               IJKM = KM_OF(IJK)
               IJKP = KP_OF(IJK)
            else 
               IJKM = IJK
               IJKP = IJK
            end if

    ! search all neighbors, it is assumed that a ghost cell has only one real 
    ! cell neighbor
            if (DO_I.AND.FLUID_AT(IMJK)) then
       ! ghost cell at east boundary
               IJK_1 = IMJK
               IJK_2 = IM_OF(IJK_1) ! search toward west for second cell
               ind_0 = I_OF(IJK)
               ind_1 = I_OF(IJK_1)
               ind_2 = I_OF(IJK_2)
               dl_0 = DX(ind_0)
               dl_1 = DX(ind_1)
               dl_2 = DX(ind_2)
            else if (DO_I.AND.FLUID_AT(IPJK)) then
       ! ghost cell at west boundary
               IJK_1 = IPJK
               IJK_2 = IP_OF(IJK_1) ! search toward east for second cell
               ind_0 = I_OF(IJK)
               ind_1 = I_OF(IJK_1)
               ind_2 = I_OF(IJK_2)
               dl_0 = DX(ind_0)
               dl_1 = DX(ind_1)
               dl_2 = DX(ind_2)
            else if (DO_J.AND.FLUID_AT(IJMK)) then
       ! ghost cell at south boundary
               IJK_1 = IJMK
               IJK_2 = JM_OF(IJK_1) ! search toward north for second cell
               ind_0 = J_OF(IJK)
               ind_1 = J_OF(IJK_1)
               ind_2 = J_OF(IJK_2)
               dl_0 = DY(ind_0)
               dl_1 = DY(ind_1)
               dl_2 = DY(ind_2)
            else if (DO_J.AND.FLUID_AT(IJPK)) then
       ! ghost cell at north boundary
               IJK_1 = IJPK
               IJK_2 = JP_OF(IJK_1) ! search toward south for second cell
               ind_0 = J_OF(IJK)
               ind_1 = J_OF(IJK_1)
               ind_2 = J_OF(IJK_2)
               dl_0 = DY(ind_0)
               dl_1 = DY(ind_1)
               dl_2 = DY(ind_2)
            else if (DO_K.AND.FLUID_AT(IJKM)) then
       ! ghost cell at bottom boundary
               IJK_1 = IJKM
               IJK_2 = KM_OF(IJK_1) ! search toward top for second cell
               ind_0 = K_OF(IJK)
               ind_1 = K_OF(IJK_1)
               ind_2 = K_OF(IJK_2)
               dl_0 = DZ(ind_0)
               dl_1 = DZ(ind_1)
               dl_2 = DZ(ind_2)
            else if (DO_K.AND.FLUID_AT(IJKP)) then
       ! ghost cell at top boundary
               IJK_1 = IJKP
               IJK_2 = KP_OF(IJK_1) ! search toward bottom for second cell
               ind_0 = K_OF(IJK)
               ind_1 = K_OF(IJK_1)
               ind_2 = K_OF(IJK_2)
               dl_0 = DZ(ind_0)
               dl_1 = DZ(ind_1)
               dl_2 = DZ(ind_2)
            else		! corners
               IJK_1 = IJK
               IJK_2 = IJK
            end if
        
            if (FLUID_AT(IJK_1)) then
               if (FLUID_AT(IJK_2)) then
! regular ghost cells found two boundary cells
                  Phi(IJK) = linearExtrapolation(Phi(IJK_1), Phi(IJK_2), dl_0, dl_1, dl_2)
               else
        ! search out of boundary
                  print*, 'search out of boundary'
                  stop
               end if
            else
               if (FLUID_AT(IJK_2)) then
        ! double layer of ghost cells?
                  print*, 'Double layer of ghost cells?'
                  stop
               else
        !corner cells
                  Phi(IJK)= ZERO ! value does not matter
               end if
            end if


         end if                 ! if ghost cells
      end do			! loop to search ghost cells

      contains

      double precision function linearExtrapolation(var_1, var_2, &
      dl_0, dl_1, dl_2)

      double precision :: var_1, var_2, dl_0, dl_1, dl_2

      linearExtrapolation = var_1 - (var_2-var_1)*&
      (dl_0+dl_1)/(dl_1+dl_2)

      end function linearExtrapolation

      end subroutine ExtrapolateGhostCells


