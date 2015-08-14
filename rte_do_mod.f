!!****************************************************************
!!**  Module RTE_DO
!!****************************************************************
      MODULE RTE_DO
! Module to solve the Radiative Transport Equation (RTE) via Discrete Ordinates (DO)
      USE PARAM
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OMEGA_DO, &
         THETA_DO,PHI_DO
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S_I, N_I
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FLUX_E_DO, &
        FLUX_N_DO, FLUX_T_DO
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U_DO, V_DO, W_DO
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: I_DO
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S_CONT_S
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: S_DES, S_CONT_G
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: KAPPA, SIGMA, EMIS
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: KAPPA_S, SIGMA_S
      DOUBLE PRECISION :: KAPPA_G, SIGMA_G
      DOUBLE PRECISION :: D_THETA, D_PHI
      DOUBLE PRECISION :: EMIS_W = 1.0d0
      DOUBLE PRECISION :: SIGMA_SB
      INTEGER :: N_THETA = 5
      INTEGER :: N_PHI = 10
      DOUBLE PRECISION :: DO_TOL = 1.0D-6
      INTEGER :: MAX_DO_ITS = 100
      INTEGER :: N_DO

      CONTAINS

      SUBROUTINE INIT_RTE_DO

      USE CONSTANT, ONLY: PI
      USE COMPAR
      USE PARAM1
      USE GEOMETRY
      USE RUN

      IMPLICIT NONE

      INTEGER :: I, J, K, IJK

      N_DO = N_THETA * N_PHI

      !Set value of stefan-boltzman constant
      IF(UNITS == 'SI') THEN
         SIGMA_SB = 5.6704d0*(10.0d0**(-8)) ! W/((m^2).K^4)
      ELSE
         SIGMA_SB = 1.355282d0*(10.0d0**(-12)) ! cal/((cm^2).sec.K^4)
      ENDIF

      ALLOCATE( OMEGA_DO(N_DO))
      ALLOCATE( THETA_DO(N_DO))
      ALLOCATE( PHI_DO(N_DO))
      ALLOCATE( S_I(N_DO,3))
      ALLOCATE( N_I(N_DO,3))
      ALLOCATE( FLUX_E_DO(DIMENSION_3,N_DO))
      ALLOCATE( FLUX_N_DO(DIMENSION_3,N_DO))
      ALLOCATE( FLUX_T_DO(DIMENSION_3,N_DO))
      ALLOCATE( U_DO(DIMENSION_3,N_DO))
      ALLOCATE( V_DO(DIMENSION_3,N_DO))
      ALLOCATE( W_DO(DIMENSION_3,N_DO))
      ALLOCATE( I_DO(DIMENSION_3,N_DO))
      ALLOCATE( KAPPA(DIMENSION_3))
      ALLOCATE( SIGMA(DIMENSION_3))
      ALLOCATE( EMIS(DIMENSION_3))
      ALLOCATE( KAPPA_S(DIMENSION_M))
      ALLOCATE( SIGMA_S(DIMENSION_M))
      ALLOCATE( S_CONT_S(DIMENSION_3,DIMENSION_M))
      ALLOCATE( S_CONT_G(DIMENSION_3))
      ALLOCATE( S_DES(DIMENSION_3))

      KAPPA_S(:) = ZERO
      SIGMA_S(:) = ZERO
      KAPPA_S(1) = 30.0d0
      KAPPA_G = 10.0d0
      SIGMA_G = 1.0d0

      I_DO(:,:) = ONE

      D_THETA = PI/DBLE(N_THETA)
      D_PHI = 2.0d0*PI/DBLE(N_PHI)

      K = 1
      DO I = 0, N_THETA - 1
         DO J = 0, N_PHI - 1
           THETA_DO(K) = D_THETA/2.0d0 + DBLE(I)*D_THETA
           PHI_DO(K) = D_PHI/2.0d0 + DBLE(J)*D_PHI

           S_I(K,1) = SIN(PHI_DO(K))*SIN(0.5d0*D_PHI)*(D_THETA - &
             COS(2.0d0*THETA_DO(K))*SIN(D_THETA))
           S_I(K,2) = COS(PHI_DO(K))*SIN(0.5d0*D_PHI)*(D_THETA - &
             COS(2.0d0*THETA_DO(K))*SIN(D_THETA))
           S_I(K,3) = 0.5d0*D_PHI*SIN(2.0d0*THETA_DO(K))*SIN(D_THETA)

           N_I(K,1) = SIN(THETA_DO(K))*SIN(PHI_DO(K))
           N_I(K,2) = COS(PHI_DO(K))*SIN(THETA_DO(K))
           N_I(K,3) = COS(THETA_DO(K))

           
           FLUX_E_DO(:,K) = AYZ(:)*S_I(K,1)
           FLUX_N_DO(:,K) = AXZ(:)*S_I(K,2)
           FLUX_T_DO(:,K) = AXY(:)*S_I(K,3)
           

           U_DO(:,K) = N_I(K,1)
           V_DO(:,K) = N_I(K,2)
           W_DO(:,K) = N_I(K,3)

           OMEGA_DO(K) = 2.0d0*SIN(THETA_DO(K))*SIN(D_THETA/2.0d0)*D_PHI
           K = K + 1
        END DO
      END DO

      END SUBROUTINE

      SUBROUTINE SOLVE_RTE_DO(SOURCE, DEBUG)
      USE PARAM1
      USE AMBM
      USE GEOMETRY
      USE INDICES
      USE COMPAR
      USE RUN
      USE LEQSOL
      USE BC
      USE TOLERANC
      USE OUTPUT
      USE RESIDUAL
      USE FUNCTIONS
      USE FUN_AVG
      USE CONSTANT
      USE FLDVAR
      USE PHYSPROP, ONLY: SMAX

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(DIMENSION_3,N_DO),INTENT(IN) :: SOURCE
      LOGICAL, INTENT(IN) :: DEBUG
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: S_C
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: S_P
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: K_DES
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: K_CONT
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: Null_Vec, Unity_Vec
      DOUBLE PRECISION, DIMENSION(DIMENSION_3,N_DO) :: TMP_I_DO
      INTEGER :: IJK, ORD, IER, M
      DOUBLE PRECISION :: RES1, MAX_RES1, NUM_RES, DEN_RES, SUM_I
      INTEGER :: IJK_RES
      LOGICAL, DIMENSION(N_DO) :: CONVERGED
      INTEGER :: ITS
      INTEGER :: LEQM, LEQI

      Null_Vec = ZERO
      Unity_Vec = ONE

      S_DES(:) = ZERO
      S_CONT_G(:) = ZERO
      S_CONT_S(:,:) = ZERO

      CONVERGED(:) = .FALSE.
      IF(.NOT.DEBUG) THEN
         CALL CALC_RAD_PROP(K_DES, K_CONT, EMIS, SIGMA)
         KAPPA(:) = K_DES(:) + K_CONT(:)
      END IF

      DO ITS = 1, MAX_DO_ITS
         IF(ALL(CONVERGED)) EXIT

         IF(.NOT. DEBUG) CALL CALC_I_W(0)

         DO ORD = 1, N_DO
            DO IJK = IJKSTART3, IJKEND3
               S_P(IJK) = (KAPPA(IJK)+SIGMA(IJK))*OMEGA_DO(ORD)*VOL(IJK)
            !Isentropic scattering
               S_C(IJK) = OMEGA_DO(ORD)*VOL(IJK)*(KAPPA(IJK)*EMIS(IJK) + &
                  SIGMA(IJK)/4.0d0/PI*SUM(I_DO(IJK,:)*OMEGA_DO(:))+ &
                  SOURCE(IJK,ORD))
            END DO
            TMP_I_DO(:,ORD) = I_DO(:,ORD)
            CALL LOCK_AMBM
            CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER)
            CALL CONV_DIF_PHI (TMP_I_DO(1,ORD),Null_Vec,DISCRETIZE(6), &
               U_DO(1,ORD), V_DO(1,ORD), W_DO(1,ORD), Flux_E_DO(1,ORD),& 
               Flux_N_DO(1,ORD), Flux_T_DO(1,ORD), 0, A_M, B_M, IER)
         
            CALL BC_DO(TMP_I_DO(1,ORD), ORD, 0, A_M, B_M, IER)
            CALL SOURCE_PHI (S_P, S_C, Unity_Vec, TMP_I_DO(1,ORD), 0, &
               A_M, B_M, IER)

            CALL CALC_RESID_S (TMP_I_DO(1,ORD), A_M, B_M, 0, NUM_RES, & 
               DEN_RES, RES1, MAX_RES1, IJK_RES, ZERO)
            IF(RES1 <= DO_TOL) THEN
               CONVERGED(ORD) = .TRUE.
               CALL UNLOCK_AMBM
               CYCLE
            ELSE
               CONVERGED(ORD) = .FALSE.
            END IF

            CALL UNDER_RELAX_S (TMP_I_DO(1,ORD),A_M,B_M,0,1.0D0,IER) 
            CALL ADJUST_LEQ (RES1, LEQ_IT(9), LEQ_METHOD(9), &
               LEQI, LEQM, IER) 

            CALL SOLVE_LIN_EQ ('I_DO', 9,TMP_I_DO(1,ORD), A_M, B_M, 0, &
               LEQI, LEQM, LEQ_SWEEP(9), LEQ_TOL(9), LEQ_PC(9), IER) 

            CALL UNLOCK_AMBM
         END DO
         I_DO(:,:) = TMP_I_DO(:,:)
      END DO

      IF(.NOT. ALL(CONVERGED)) THEN
         PRINT *, 'RTE DO Iterations did not converge!'
         CALL MFIX_EXIT(myPE)
      END IF

      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         SUM_I = SUM(I_DO(IJK,:)*OMEGA_DO(:))
         S_DES(IJK) = SUM_I*K_DES(IJK)
         S_CONT_G(IJK) = SUM_I*EP_G(IJK)*KAPPA_G
         DO M = 1, SMAX
            S_CONT_S(IJK,M) = SUM_I*EP_S(IJK,M)*KAPPA_S(M)
         END DO
      END DO

      END SUBROUTINE

      SUBROUTINE MANUFACTURE_SOLUTION
      USE PARAM1
      USE AMBM
      USE GEOMETRY, ONLY: DX, DY, DZ, XLENGTH, YLENGTH, ZLENGTH, &
         IMIN3, IMAX3, JMIN3, JMAX3, KMIN3, KMAX3
      USE INDICES
      USE COMPAR
      USE RUN
      USE LEQSOL
      USE BC
      USE TOLERANC
      USE OUTPUT
      USE RESIDUAL
      USE FUNCTIONS
      USE FUN_AVG
      USE CONSTANT, ONLY: PI
      USE PHYSPROP, ONLY: SMAX

      DOUBLE PRECISION, DIMENSION(DIMENSION_3,N_DO) :: I_SOL, SOURCE
      DOUBLE PRECISION :: X, Y, Z
      DOUBLE PRECISION :: TEMP = 463.0d0
      INTEGER :: I, J, K, IJK, ORD
      DOUBLE PRECISION, DIMENSION(3) :: GRAD_I
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: INT_I
      DOUBLE PRECISION, DIMENSION(N_DO) :: NUM_RES, DEN_RES 
      DOUBLE PRECISION :: NUM_RES_TOT, DEN_RES_TOT, SUM_I
      
      !Set up parameters
      I_DO(:,:) = ONE
      KAPPA(:) = ONE
      EMIS(:) = SIGMA_SB*TEMP**4/PI
      SIGMA(:) = ZERO
      

      DO I = IMIN3, IMAX3
         DO J = JMIN3, JMAX3
            DO K = KMIN3, KMAX3
               !Set up X, Y, and Z
               IJK = FUNIJK(I,J,K)
               IF(I .EQ. IMIN3) THEN
                  X = ZERO
               ELSE IF (I .EQ. IMAX3) THEN
                  X = XLENGTH
               ELSE
                  X = HALF*DX(1)+DBLE(I-2)*DX(1)
               END IF

               IF(J .EQ. JMIN3) THEN
                  Y = ZERO
               ELSE IF (J .EQ. JMAX3) THEN
                  Y = YLENGTH
               ELSE
                  Y = HALF*DY(1)+DBLE(J-2)*DY(1)
               END IF

               IF(K .EQ. KMIN3) THEN
                  Z = ZERO
               ELSE IF (K .EQ. KMAX3) THEN
                  Z = ZLENGTH
               ELSE
                  Z = HALF*DZ(1)+DBLE(K-2)*DZ(1)
               END IF
               
               DO ORD = 1, N_DO
                  !set up the solution, this needs to be done first in
                  !case the solution is angle dependent
                  I_SOL(IJK,ORD)=(X**3+3*Y**2+Z**2) * &
                     (COS(THETA_DO(ORD))**2 + 1.0d0)
               END DO

               !define integral of the solution over all solid angles
               INT_I(IJK) = (X**3+3*Y**2+Z**2)*16.0d0/3.0d0*PI

               DO ORD = 1, N_DO
                  !set up the source term
                  !define gradient
                  GRAD_I(1) = (3.0d0*X**2) * &
                     (COS(THETA_DO(ORD))**2 + 1.0d0)
                  GRAD_I(2) = (6.0d0*Y) * &
                     (COS(THETA_DO(ORD))**2 + 1.0d0)
                  GRAD_I(3) = (2.0d0*Z) * &
                     (COS(THETA_DO(ORD))**2 + 1.0d0)

                  SOURCE(IJK,ORD) = (KAPPA(IJK)+SIGMA(IJK)) * &
                  I_SOL(IJK,ORD)-KAPPA(IJK)*SIGMA_SB*TEMP**4/PI - &
                  SIGMA(IJK)/4.0d0/PI*INT_I(IJK) + &
                  DOT_PRODUCT(GRAD_I,N_I(ORD,:))
               END DO

               !Set up boundary conditions
               IF(.NOT.FLUID_AT(IJK)) I_DO(IJK,:) = I_SOL(IJK,:)

            END DO
         END DO
      END DO
      !Solve
      CALL SOLVE_RTE_DO(SOURCE, .TRUE.)

      !Calculate and print residuals
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         SUM_I = SUM(I_DO(IJK,:)*OMEGA_DO(:))
         NUM_RES_TOT = NUM_RES_TOT + (SUM_I-INT_I(IJK))**2
         DEN_RES_TOT = DEN_RES_TOT + (INT_I(IJK))**2
         DO ORD = 1, N_DO
            NUM_RES(ORD)=NUM_RES(ORD)+(I_DO(IJK,ORD)-I_SOL(IJK,ORD))**2
            DEN_RES(ORD)=DEN_RES(ORD)+(I_SOL(IJK,ORD))**2
         END DO
      END DO

      DO ORD = 1, N_DO
         PRINT *,'Residual for ordinate', &
           ORD,THETA_DO(ORD),PHI_DO(ORD),SQRT(NUM_RES(ORD)/DEN_RES(ORD))
      END DO

      PRINT *,'Total residual',SQRT(NUM_RES_TOT/DEN_RES_TOT)
      CALL MFIX_EXIT(myPE)
      
      END SUBROUTINE

      SUBROUTINE CHECK_MATRIX
      USE CONSTANT
      USE PARAM1
      USE AMBM
      USE GEOMETRY
      USE INDICES
      USE COMPAR
      USE RUN
      USE LEQSOL
      USE BC
      USE TOLERANC
      USE OUTPUT
      USE RESIDUAL
      USE FUNCTIONS
      USE FUN_AVG

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(DIMENSION_3,N_DO) :: TMP_I_DO
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: K_DES
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: K_CONT
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: S_C
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: S_P
      INTEGER :: ORD, IJK, IER, I_WRITE, I, J, K
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: Null_Vec, Unity_Vec

      Null_Vec = ZERO
      Unity_Vec = ONE

      !Check that solid angles sum to 4*PI
      IF(ABS(SUM(OMEGA_DO(:))-4.0d0*PI) > 0.0001d0) PRINT *, &
         'Sum of solid angles does not sum to 4*PI'

      !For debugging
      TMP_I_DO(:,:) = I_DO(:,:)
      !End debugging

      CALL CALC_I_W(0)

      !Sanity checks for debugging
      DO IJK = IJKSTART3, IJKEND3
         IF(FLUID_AT(IJK)) THEN
            IF(.NOT.ALL(ABS(TMP_I_DO(IJK,:)-I_DO(IJK,:))<=0.0001d0))&
               PRINT *,'CALC_I_W modified fluid cell'
         END IF
      END DO
      !end debugging

      CALL CALC_RAD_PROP(K_DES, K_CONT, EMIS, SIGMA)
      KAPPA(:) = K_DES(:) + K_CONT(:)

      IF(ANY((KAPPA(:) <= 0.0d0).AND.(FLAG(:).EQ.1))) &
         PRINT *,'Zero or negative absorptivity'
      IF(ANY((EMIS(:) <= 0.0d0).AND.(FLAG(:).EQ.1))) &
         PRINT *,'Zero or negative emissivity'

      OPEN (unit = 2, file = "rad_prop")
      DO I_WRITE = 1, DIMENSION_3
         WRITE (2,*) I_WRITE, FLUID_AT(I_WRITE), KAPPA(I_WRITE), &
         K_DES(I_WRITE), K_CONT(I_WRITE), EMIS(I_WRITE), &
         SIGMA(I_WRITE)
      END DO
      CLOSE(2)

      DO ORD = 1, N_DO
         DO IJK = IJKSTART3, IJKEND3
            S_P(IJK) = (KAPPA(IJK)+SIGMA(IJK))*OMEGA_DO(ORD)*VOL(IJK)
            !Isentropic scattering
            S_C(IJK) = OMEGA_DO(ORD)*VOL(IJK)*(EMIS(IJK) + &
               SIGMA(IJK)/4.0d0/PI*SUM(I_DO(IJK,:)*OMEGA_DO(:)))
         END DO
         TMP_I_DO(:,ORD) = I_DO(:,ORD)
         CALL LOCK_AMBM
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER)
         CALL CONV_DIF_PHI (TMP_I_DO(1,ORD),Null_Vec,DISCRETIZE(6), &
            U_DO(1,ORD), V_DO(1,ORD), W_DO(1,ORD), Flux_E_DO(1,ORD), & 
            Flux_N_DO(1,ORD), Flux_T_DO(1,ORD), 0, A_M, B_M, IER)
         
         CALL BC_DO(TMP_I_DO(1,ORD), ORD, 0, A_M, B_M, IER)
         CALL SOURCE_PHI (S_P, S_C, Unity_Vec, TMP_I_DO(1,ORD), 0, &
            A_M, B_M, IER)
         CALL UNLOCK_AMBM

         !Sanity checks for debugging:
         IF(ANY(A_M(:,0,0) >= ZERO)) THEN 
            PRINT *,'Diagonal coefficient positive or zero',ORD
            OPEN (unit = 2, file = "do_matrix")
            OPEN (unit = 3, file = "do_flux_vel")
            DO I_WRITE = 1, DIMENSION_3
               WRITE (2,*) A_M(I_WRITE,-3,0), A_M(I_WRITE,-2,0), & 
               A_M(I_WRITE,-1,0), A_M(I_WRITE,0,0), A_M(I_WRITE,1,0), &
               A_M(I_WRITE,2,0), A_M(I_WRITE,3,0), B_M(I_WRITE,0)
               WRITE (3,*) FLUX_E_DO(I_WRITE,ORD), & 
               FLUX_N_DO(I_WRITE,ORD), FLUX_T_DO(I_WRITE,ORD), &
               U_DO(I_WRITE,ORD), V_DO(I_WRITE,ORD), W_DO(I_WRITE,ORD)
            END DO
            CLOSE(2)
            CLOSE(3)
         END IF

      END DO

      CALL MFIX_EXIT(myPE)

      END SUBROUTINE

      SUBROUTINE CALC_RAD_PROP(K_DES,K_CONT,EMIS,SIGMA)

      USE CONSTANT
      USE DISCRETELEMENT
      USE FLDVAR
      USE GEOMETRY
      USE INDICES
      USE INTERPOLATION
      USE PARAM1
      USE RUN
      USE COMPAR
      USE FUNCTIONS
      USE DES_THERMO
      USE PHYSPROP, ONLY: SMAX
      
      IMPLICIT NONE
!----------------------------------------------------------------
! Passed variables
!----------------------------------------------------------------
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(OUT) :: K_DES, &
        K_CONT, EMIS, SIGMA
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
      DOUBLE PRECISION :: OVOL, OPVOL
! percent of particle in a cell and projected area of particle
      DOUBLE PRECISION :: P_IN_CELL, A_P
!----------------------------------------------
	

      K_DES(:) = ZERO
      K_CONT(:) = ZERO
      EMIS(:) = ZERO
      SIGMA(:) = ZERO

      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         OVOL = ONE/VOL(IJK)
         IF(DISCRETE_ELEMENT) THEN
            DO L = 1, MAX_PIP
               IF(.NOT.PEA(L,1)) CYCLE
               IF(PART_VOL_INTERSEC(IJK,L) == ZERO) CYCLE
               P_IN_CELL = PART_VOL_INTERSEC(IJK,L)/PVOL(L)
               M = PIJK(L,5) + SMAX
               A_P = PI*DES_RADIUS(L)**2
               K_DES(IJK) = K_DES(IJK) + DES_EM(M)*A_P*OVOL*P_IN_CELL
               EMIS(IJK) = EMIS(IJK) + DES_EM(M)*A_P* &
                  SIGMA_SB*DES_T_S_NEW(L)**4/PI*OVOL*P_IN_CELL
               SIGMA(IJK) = SIGMA(IJK)+(ONE-DES_EM(M))*A_P*  &
                  OVOL*P_IN_CELL
            END DO
         END IF

         K_CONT(IJK) = K_CONT(IJK) + KAPPA_G*EP_G(IJK)
         EMIS(IJK) = EMIS(IJK) + KAPPA_G*SIGMA_SB*  &
            T_G(IJK)**4/PI*EP_G(IJK) 
         SIGMA(IJK) = SIGMA(IJK) + SIGMA_G*EP_G(IJK)

         DO M = 1, SMAX
            K_CONT(IJK) = K_CONT(IJK) + KAPPA_S(M)*EP_S(IJK,M)
            EMIS(IJK) = EMIS(IJK) + KAPPA_S(M)*SIGMA_SB*  &
               T_S(IJK,M)**4/PI*EP_S(IJK,M)
            SIGMA(IJK) = SIGMA(IJK) + SIGMA_S(M)*EP_S(IJK,M)
         END DO
      END DO

      END SUBROUTINE CALC_RAD_PROP

      SUBROUTINE CALC_I_W(M)
      USE PARAM1
      USE MATRIX
      USE GEOMETRY
      USE INDICES
      USE BC
      USE COMPAR
      USE FUN_AVG
      USE FUNCTIONS
      USE FLDVAR
      USE CONSTANT

      IMPLICIT NONE

! Phase of wall BC to use
      INTEGER, INTENT(IN) :: M
! Indices
      INTEGER :: I1, J1, K1, K2, IJK, IJK_F, ORD
! Wall temp
      DOUBLE PRECISION :: Tw, SUM

      IF (DO_K) THEN
! bottom xy plane
         K1 = 1
         DO J1 = jmin3, jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
               IJK_F = KP_OF(IJK)
               IF(M == 0) THEN
                  Tw = (T_G(IJK) + T_G(IJK_F))/2.0d0
               ELSE
                  Tw = (T_S(IJK,M) + T_S(IJK_F,M))/2.0d0
               END IF
               SUM = 0.0d0
               DO ORD = 1, N_DO
                  IF(-N_I(ORD,3)>ZERO) SUM=SUM-I_DO(IJK_F,ORD)*S_I(ORD,3)
               ENDDO
               I_DO(IJK,:) = (ONE-EMIS_W)/PI*SUM+EMIS_W*SIGMA_SB*Tw**4/PI
               
            ENDDO
         ENDDO

! top xy plane
         K1 = KMAX2
         DO J1 = jmin3, jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
               IJK_F = KM_OF(IJK)
               IF(M == 0) THEN
                  Tw = (T_G(IJK) + T_G(IJK_F))/2.0d0
               ELSE
                  Tw = (T_S(IJK,M) + T_S(IJK_F,M))/2.0d0
               END IF
               SUM = 0.0d0
               DO ORD = 1, N_DO
                  IF(N_I(ORD,3)>ZERO) SUM=SUM+I_DO(IJK_F,ORD)*S_I(ORD,3)
               ENDDO
               I_DO(IJK,:) = (ONE-EMIS_W)/PI*SUM+EMIS_W*SIGMA_SB*Tw**4/PI
            ENDDO
         ENDDO
      ENDIF

! south xz plane
      J1 = 1
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IJK_F = JP_OF(IJK)
            IF(M == 0) THEN
               Tw = (T_G(IJK) + T_G(IJK_F))/2.0d0
            ELSE
               Tw = (T_S(IJK,M) + T_S(IJK_F,M))/2.0d0
            END IF
            SUM = 0.0d0
            DO ORD = 1, N_DO
               IF(-N_I(ORD,2)>ZERO) SUM=SUM-I_DO(IJK_F,ORD)*S_I(ORD,2)
            ENDDO
            I_DO(IJK,:) = (ONE-EMIS_W)/PI*SUM+EMIS_W*SIGMA_SB*Tw**4/PI
         ENDDO
      ENDDO

! north xz plane
      J1 = JMAX2
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IJK_F = JM_OF(IJK)
            IF(M == 0) THEN
               Tw = (T_G(IJK) + T_G(IJK_F))/2.0d0
            ELSE
               Tw = (T_S(IJK,M) + T_S(IJK_F,M))/2.0d0
            END IF
            SUM = 0.0d0
            DO ORD = 1, N_DO
               IF(N_I(ORD,2)>ZERO) SUM=SUM+I_DO(IJK_F,ORD)*S_I(ORD,2)
            ENDDO
            I_DO(IJK,:) = (ONE-EMIS_W)/PI*SUM+EMIS_W*SIGMA_SB*Tw**4/PI
         ENDDO
      ENDDO

! west yz plane
      I1 = imin2
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IJK_F = IP_OF(IJK)
            IF(M == 0) THEN
               Tw = (T_G(IJK) + T_G(IJK_F))/2.0d0
            ELSE
               Tw = (T_S(IJK,M) + T_S(IJK_F,M))/2.0d0
            END IF
            SUM = 0.0d0
            DO ORD = 1, N_DO
               IF(-N_I(ORD,1)>ZERO) SUM=SUM-I_DO(IJK_F,ORD)*S_I(ORD,1)
            ENDDO
            I_DO(IJK,:) = (ONE-EMIS_W)/PI*SUM+EMIS_W*SIGMA_SB*Tw**4/PI
         ENDDO
      ENDDO

! east yz plane
      I1 = IMAX2
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IJK_F = IM_OF(IJK)
            IF(M == 0) THEN
               Tw = (T_G(IJK) + T_G(IJK_F))/2.0d0
            ELSE
               Tw = (T_S(IJK,M) + T_S(IJK_F,M))/2.0d0
            END IF
            SUM = 0.0d0
            DO ORD = 1, N_DO
               IF(N_I(ORD,1)>ZERO) SUM=SUM+I_DO(IJK_F,ORD)*S_I(ORD,1)
            ENDDO
            I_DO(IJK,:) = (ONE-EMIS_W)/PI*SUM+EMIS_W*SIGMA_SB*Tw**4/PI
         ENDDO
      ENDDO    

      END SUBROUTINE


      SUBROUTINE BC_DO(VAR, ORD, M, A_M, B_M, IER)
! Modules
!--------------------------------------------------------------------//
      USE PARAM1
      USE MATRIX
      USE GEOMETRY
      USE INDICES
      USE BC
      USE COMPAR
      USE FUN_AVG
      USE FUNCTIONS
      IMPLICIT NONE

! Dummy arguments
!--------------------------------------------------------------------//
! The field variable being solved for:
      DOUBLE PRECISION, INTENT(IN) :: VAR(DIMENSION_3)
! Ordinate being solved for
      INTEGER, INTENT(IN) :: ORD
! Phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local variables
!--------------------------------------------------------------------//
! Indices
      INTEGER :: I1, J1, K1, IJK
!--------------------------------------------------------------------//

      IF (DO_K) THEN
! bottom xy plane
         K1 = 1
         DO J1 = jmin3, jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
! first set the flow boundary cell value equal to zero
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
! then set cell value according to inflow/outflow
               IF(-N_I(ORD,3) < ZERO) THEN
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = -VAR(IJK)
               ELSE
                  A_M(IJK,T,M) = ONE
               END IF
            ENDDO
         ENDDO

! top xy plane
         K1 = KMAX2
         DO J1 = jmin3, jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
! first set the flow boundary cell value equal to zero
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
!then set cell value according to inflow/outflow
               IF(N_I(ORD,3) < ZERO) THEN
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = -VAR(IJK)
               ELSE
                  A_M(IJK,B,M) = ONE
               END IF
            ENDDO
         ENDDO
      ENDIF

! south xz plane
      J1 = 1
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
! first set the flow boundary cell value equal to zero
            A_M(IJK,E,M) = ZERO
            A_M(IJK,W,M) = ZERO
            A_M(IJK,N,M) = ZERO
            A_M(IJK,S,M) = ZERO
            A_M(IJK,T,M) = ZERO
            A_M(IJK,B,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO
! then set cell value according to inflow/outflow
            IF(-N_I(ORD,2) < ZERO) THEN
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = -VAR(IJK)
            ELSE
               A_M(IJK,N,M) = ONE
            END IF
         ENDDO
      ENDDO

! north xz plane
      J1 = JMAX2
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
! first set the flow boundary cell value equal to zero
            A_M(IJK,E,M) = ZERO
            A_M(IJK,W,M) = ZERO
            A_M(IJK,N,M) = ZERO
            A_M(IJK,S,M) = ZERO
            A_M(IJK,T,M) = ZERO
            A_M(IJK,B,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO
! then set cell value according to inflow/outflow
            IF(N_I(ORD,2) < ZERO) THEN
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = -VAR(IJK)
            ELSE
               A_M(IJK,S,M) = ONE
            END IF
         ENDDO
      ENDDO

! west yz plane
      I1 = imin2
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
! first set the flow boundary cell value equal to zero
            A_M(IJK,E,M) = ZERO
            A_M(IJK,W,M) = ZERO
            A_M(IJK,N,M) = ZERO
            A_M(IJK,S,M) = ZERO
            A_M(IJK,T,M) = ZERO
            A_M(IJK,B,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO
! then set cell value according to inflow/outflow
            IF(-N_I(ORD,1) < ZERO) THEN
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = -VAR(IJK)
            ELSE
               A_M(IJK,E,M) = ONE
            END IF
         ENDDO
      ENDDO

! east yz plane
      I1 = IMAX2
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
! first set the flow boundary cell value equal to zero
            A_M(IJK,E,M) = ZERO
            A_M(IJK,W,M) = ZERO
            A_M(IJK,N,M) = ZERO
            A_M(IJK,S,M) = ZERO
            A_M(IJK,T,M) = ZERO
            A_M(IJK,B,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO
! then set cell value according to inflow/outflow
            IF(N_I(ORD,1) < ZERO) THEN
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = -VAR(IJK)
            ELSE
               A_M(IJK,W,M) = ONE
            END IF
         ENDDO
      ENDDO

      END SUBROUTINE

      END MODULE
