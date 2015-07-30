!!****************************************************************
!!** module rte_p1
!!****************************************************************
      module rte_p1
! module to solve the P1 Radiative Transfer Equation (RTE)
      use param
      use param1
      use ambm
      use geometry
      use indices
      use compar
      use run
      use leqsol
      use bc
      use toleranc
      use output
      use residual
      use rad_param
	  use functions
	  use fun_avg
      implicit none

      contains
!------------------------------------------------------------------
! subroutine to solve generic P1-like PDE with boundary conditions
      subroutine solve_p1_eqn(Sq,G, beta, w, omega, S, ew, Sw)
! input and output
      real(8), intent(out), dimension(DIMENSION_3) :: Sq, G
! Sq heat source = -del.q
      real(8), intent(in), dimension(DIMENSION_3) :: beta, w, S, omega
! beta: extinction coefficient, omega: scattering albedo
! w: general P1 equation proportional constant, S: source term in general P1 equation
      real(8), intent(in), dimension(DIMENSION_BC) :: ew, Sw
! ew: wall emissivity, Sw: wall emission source 
! local variables
      real(8), dimension(DIMENSION_3) :: Null_Vec,Unity_Vec, invBeta
      real(8), dimension(DIMENSION_3) :: S_C, S_P
      real(8), dimension(DIMENSION_BC) :: BC_HW, BC_C, BC_W, BC_F
      integer(8) :: M
      character(len=8) :: Vname
!                      Indices
      INTEGER          IJK, I, J, K, L, TOTAL_BC
!                      linear equation solver method and iterations
      INTEGER          LEQM, LEQI
!                      temporary variables in residual computation 
      real(8) 	:: res1, mres1, num_res, den_res
      INTEGER          ires1 , ier
! initialize variables
      G=ZERO
      Sq = ZERO
      Null_Vec = ZERO
      Unity_Vec = ONE
      Vname = "G"
	
      if (maxval(beta)<BETAMIN) then
	! optically thin bypass P1 calculation 
         Sq = (omega-1.d0)*beta*S ! -del_dot_q
!	write(*,*) "P1 cacluation bypassed"
         return
      endif 

      call lock_ambm

      CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER)

! it is assumed that extrapolation has been done before entering this subroutine
! variables need extrapolation to ghost cells include beta, w, omega 
      invBeta = one/beta

      CALL CONV_DIF_G (G, invBeta, DISCRETIZE(6), Null_Vec, &
      Null_Vec, Null_Vec, Null_Vec, Null_Vec, Null_Vec, 0, A_M&
      , B_M, IER)               ! convection diffusion term discretization

!C-------------------------------------------------------------
!C  calculate boundary wall transfer
!c-------------------------------------------------------------
      BC_HW = ZERO              ! default values


      do L = 1,DIMENSION_BC
         if (.not.BC_DEFINED(L)) then
            BC_C(L) = UNDEFINED
            BC_HW(L) = UNDEFINED
         else
            BC_C(L) = ew(L)/(2d0-ew(L))*1.5d0*Sw(L)
            BC_HW(L) =  ew(L)/(2d0-ew(L))*1.5d0
         end if
      end do
      BC_W = ZERO
      BC_F = ZERO

      CALL BC_G2 (G, BC_F, BC_W, BC_HW, BC_C,beta, 0, A_M, B_M, IER)
	
      S_C = 3.0d0*VOL*S 
      S_P = w*beta*VOL
	
      CALL SOURCE_PHI (S_P, S_C, Unity_Vec, G, 0, A_M, B_M, IER)

      M= 0
! linear solver
      CALL CALC_RESID_S (G, A_M, B_M, M, num_res, den_res, res1, &
      mres1, ires1, ZERO, IER) 
      CALL UNDER_RELAX_S (G, A_M, B_M, M, 1.0D0, IER) 
      CALL ADJUST_LEQ (res1, LEQ_IT(9), LEQ_METHOD(9), &
      LEQI, LEQM, IER) 

      CALL SOLVE_LIN_EQ (Vname, 9, G, A_M, B_M, 0, LEQI, LEQM, &
      LEQ_SWEEP(9), LEQ_TOL(9), LEQ_PC(9), IER) 

      call unlock_ambm
	
!radiative heat source is absorption - emssion
      Sq = (1-omega)*beta*(G-S) ! -del_dot_q
      end subroutine solve_P1_eqn


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BC_G2(VAR, BC_phif, BC_Phiw, BC_hw_Phi, BC_C_Phi,BETA, M,  C
!                              A_m, B_m, IER)                          C
!  Purpose: Set up the phi boundary conditions                         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-APR-04  C
!  Reviewer:                                          Date:            C
!  Purpose: include the variable (VAR) in the interface                C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE BC_G2(VAR, BC_PHIF,BC_PHIW,BC_HW_PHI,BC_C_PHI,BETA,&
      M,A_M,B_M,IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Added BETA for extinction coefficient variations at walls J.CAI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE toleranc 
      USE run
      USE physprop
      USE fldvar
      !USE visc_s
      USE geometry
      USE output
      USE indices
      USE bc
      USE compar    
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Boundary condition
      INTEGER          L
!
!                      Indices
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK, &
                      IM, JM, KM
!
!                      Solids phase
      INTEGER          M
!
!                      The field variable being solved for
      DOUBLE PRECISION VAR(DIMENSION_3)
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Boundary conditions
      DOUBLE PRECISION BC_phif(DIMENSION_BC), BC_Phiw(DIMENSION_BC), &
                      BC_hw_Phi(DIMENSION_BC), BC_C_Phi(DIMENSION_BC)

      DOUBLE PRECISION BETA(DIMENSION_3)

	DOUBLE PRECISION BetaW ! extinction coefficient at the wall face
!-----------------------------------------------

!
!  Set up the default walls as non-conducting.
!
!	PRINT*, 'BC_G.F LOADED'
      IF (DO_K) THEN 
         K1 = 1 
!$omp    parallel do private(IJK, J1, I1)
         DO J1 = jmin3, jmax3 
            DO I1 = imin3, imax3 
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IJK = FUNIJK(I1,J1,K1) 
               IF (DEFAULT_WALL_AT(IJK)) THEN 
!
                  A_M(KP_OF(IJK),B,M) = ZERO 
!
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
         K1 = KMAX2 
!$omp    parallel do private(IJK, J1, I1)
         DO J1 = jmin3, jmax3 
            DO I1 = imin3, imax3 
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IJK = FUNIJK(I1,J1,K1) 
               IF (DEFAULT_WALL_AT(IJK)) THEN 
!
                  A_M(KM_OF(IJK),T,M) = ZERO 
!
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
      ENDIF 
!
      J1 = 1 
!$omp    parallel do private(IJK, K1, I1)
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3 
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) THEN 
!
               A_M(JP_OF(IJK),S,M) = ZERO 
!
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      
      J1 = JMAX2 
!$omp    parallel do private(IJK, K1, I1)
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3 
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) THEN 
!
               A_M(JM_OF(IJK),N,M) = ZERO 
!
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 

      I1 = imin2 
!$omp    parallel do private(IJK, K1, J1)
      DO K1 = kmin3, kmax3 
         DO J1 = jmin3, jmax3 
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) THEN 
!
               A_M(IP_OF(IJK),W,M) = ZERO 
!
               A_M(IJK,E,M) = ONE 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 

      I1 = IMAX2 
!$omp    parallel do private(IJK, K1, J1)
      DO K1 = kmin3, kmax3 
         DO J1 = jmin3, jmax3 
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) THEN 
!
               A_M(IM_OF(IJK),E,M) = ZERO 
!
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ONE 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      
      !first set the bc for walls then overwrite where ever inflow/outflows are
      !defined so that the order in which the bcs are defined in the data file
      !does not matter.  Here set wall bcs . . .
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
!            IF (BC_TYPE(L)=='NO_SLIP_WALL' .OR. BC_TYPE(L)=='FREE_SLIP_WALL'&
!                .OR. BC_TYPE(L)=='PAR_SLIP_WALL') THEN 
            IF (.TRUE.) THEN
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!$omp    parallel do private(IJK, K, J, I, IM, JM, KM)
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 		     
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IJK = FUNIJK(I,J,K) 
                        IM = IM1(I) 
                        JM = JM1(J) 
                        KM = KM1(K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = VAR(IJK) 
                        IF (FLUID_AT(EAST_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,E,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
					BetaW = HALF*(BETA(IJK)+BETA(EAST_OF(IJK)))
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)*BetaW+ODX_E(I)) 
                              A_M(IJK,E,M) = -(HALF*BC_HW_PHI(L)*BetaW-ODX_E(I)) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 ))*BetaW 
                           ENDIF 
                        ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,W,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
					BetaW = HALF*(BETA(IJK)+BETA(WEST_OF(IJK)))
                              A_M(IJK,W,M) = -(HALF*BC_HW_PHI(L)*BetaW-ODX_E(IM)) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)*BetaW+ODX_E(IM)) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 )) *BetaW
                           ENDIF 
                        ELSE IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,N,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
					BetaW = HALF*(BETA(IJK)+BETA(NORTH_OF(IJK)))
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)*BetaW+ODY_N(J)) 
                              A_M(IJK,N,M) = -(HALF*BC_HW_PHI(L)*BetaW-ODY_N(J)) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 )) *BetaW
                           ENDIF 
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,S,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
					BetaW = HALF*(BETA(IJK)+BETA(SOUTH_OF(IJK)))
                              A_M(IJK,S,M) = -(HALF*BC_HW_PHI(L)*BetaW-ODY_N(JM)) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)*BetaW+ODY_N(JM)) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 )) *BetaW
                           ENDIF 
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,T,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
					BetaW = HALF*(BETA(IJK)+BETA(TOP_OF(IJK)))
                              A_M(IJK,0,M)=-(HALF*BC_HW_PHI(L)*BetaW +OX(I)*ODZ_T(K)) 
                              A_M(IJK,T,M)=-(HALF*BC_HW_PHI(L)*BetaW -OX(I)*ODZ_T(K)) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 )) *BetaW
                           ENDIF 
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,B,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
					BetaW = HALF*(BETA(IJK)+BETA(BOTTOM_OF(IJK)))
                              A_M(IJK,B,M) = -(HALF*BC_HW_PHI(L)*BetaW -OX(I)*ODZ_T(KM&
                                 )) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)*BetaW +OX(I)*ODZ_T(KM&
                                 )) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 )) *BetaW
                           ENDIF 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      


      RETURN  
      END SUBROUTINE BC_G2 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
!// 350 1206 change do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:
!  CONV_DIF_Phi(Phi, Dif, Disc, Uf, Vf, Wf, Flux_E, Flux_N, Flux_T, M, A_m, B_m, IER)    C
!  Purpose: Determine convection diffusion terms for a sclar phi       C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C

!  The diffusion at the flow boundaries is prevented by setting the 
!  diffusion coefficients at boundary cells to zero and then using a 
!  harmonic average to calculate the boundary diffusivity.  The value
!  diffusivities at the boundaries are checked in check_data_30.  Ensure
!  that harmonic avergaing is used in this routine. 
!  See source_phi                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_G(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1
      USE run 
      USE geometry
      USE compar
      USE sendrecv
!      Use xsi_array
      USE mpi_utility
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Discretizationindex
      INTEGER          Disc
!
!                      Velocity components
      DOUBLE PRECISION Uf(DIMENSION_3), Vf(DIMENSION_3), Wf(DIMENSION_3) 
!
!                      Mass flux components
      DOUBLE PRECISION Flux_E(DIMENSION_3), Flux_N(DIMENSION_3), Flux_T(DIMENSION_3) 
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Error index
      INTEGER          IER

!

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	IF DEFERRED CORRECTION IS USED WITH THE SCALAR TRANSPORT EQN.
!
!print*, 'conv_dif_phi overloaded'
	IF(DEF_COR)THEN
		print*, 'No defered correction for incident radiation'
		stop
	ELSE
!
!	NO DEFERRED CORRECTION IS USED WITH THE SCALAR TRANSPORT EQN.
!
	  IF (DISC == 0) THEN                        
            CALL CONV_DIF_G0(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER)
	  ELSE
           	print*, 'DISC must be 0'
		stop
          ENDIF
	ENDIF 
	
        CALL DIF_PHI_IS (DIF, A_M, B_M, M, IER)

        RETURN  
      END SUBROUTINE CONV_DIF_G 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:
!   CONV_DIF_Phi0(Phi, Dif, Disc, Uf, Vf, Wf, Flux_E,Flux_N,Flux_T, M, A_m, B_m, IER)
!  Purpose: Determine convection diffusion terms for Phi balance       C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_phi                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_G0(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE toleranc 
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Discretizationindex
      INTEGER          Disc
!
!                      Velocity components
      DOUBLE PRECISION Uf(DIMENSION_3), Vf(DIMENSION_3), Wf(DIMENSION_3) 
!
!                      Mass flux components
      DOUBLE PRECISION Flux_E(DIMENSION_3), Flux_N(DIMENSION_3), Flux_T(DIMENSION_3) 
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IJK, IPJK, IJPK, IJKE, IJKN,&
                       IJKP, IJKT

      INTEGER          IMJK, IM, IJKW
      INTEGER          IJMK, JM, IJKS
      INTEGER          IJKM, KM, IJKB
!
!                      Face velocity
      DOUBLE PRECISION V_f
!
!                      Difusion parameter
      DOUBLE PRECISION D_f
!
!-----------------------------------------------
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do                                              &
!$omp&     private(I,  J, K,  IJK,  IPJK, IJPK, IJKE, IJKN,         &
!$omp&             IJKP, IJKT,  V_f, D_f,                    &
!$omp&             IMJK, IM, IJKW,                                  &
!$omp&             IJMK, JM, IJKS,                                  &
!$omp&             IJKM, KM,  IJKB)                     
      DO IJK = ijkstart3, ijkend3
!
       I = I_OF(IJK)
       J = J_OF(IJK)
       K = K_OF(IJK)
!
         IF (FLUID_AT(IJK)) THEN 
!
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
!print*, 'IPJK & IJKE', IPJK, IJKE 
!
!
!           East face (i+1/2, j, k)
            V_F = UF(IJK) 
!C----------------------------
!C OVERWRITE INDICES TO GHOST CELLS, JCAI
!            D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*AYZ(IJK) 
            D_F = AVG_X_H(DIF(IJK),DIF(IPJK),I)*ODX_E(I)*AYZ(IJK) 
!C----------------------------
            IF (V_F >= ZERO) THEN 
               A_M(IJK,E,M) = D_F 
               A_M(IPJK,W,M) = D_F + FLUX_E(IJK) 
            ELSE 
               A_M(IJK,E,M) = D_F - FLUX_E(IJK) 
               A_M(IPJK,W,M) = D_F 
            ENDIF 
!
!
!           North face (i, j+1/2, k)
            V_F = VF(IJK) 
!C----------------------------
!C OVERWRITE INDICES TO GHOST CELLS, JCAI
!            D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*AXZ(IJK) 
            D_F = AVG_Y_H(DIF(IJK),DIF(IJPK),J)*ODY_N(J)*AXZ(IJK) 
!C----------------------------
            IF (V_F >= ZERO) THEN 
               A_M(IJK,N,M) = D_F 
               A_M(IJPK,S,M) = D_F + FLUX_N(IJK) 
            ELSE 
               A_M(IJK,N,M) = D_F - FLUX_N(IJK) 
               A_M(IJPK,S,M) = D_F 
            ENDIF 
!
!           Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               V_F = WF(IJK) 
!C----------------------------
!C OVERWRITE INDICES TO GHOST CELLS, JCAI
!               D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*AXY(IJK) 
               D_F = AVG_Z_H(DIF(IJK),DIF(IJKP),K)*OX(I)*ODZ_T(K)*AXY(IJK) 
!C----------------------------
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,T,M) = D_F 
                  A_M(IJKP,B,M) = D_F + FLUX_T(IJK) 
               ELSE 
                  A_M(IJK,T,M) = D_F - FLUX_T(IJK) 
                  A_M(IJKP,B,M) = D_F 
               ENDIF 
            ENDIF 
!
!
!           West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               V_F = UF(IMJK) 
!C----------------------------
!C OVERWRITE INDICES TO GHOST CELLS, JCAI
!               D_F = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*ODX_E(IM)*AYZ(IMJK) 
               D_F = AVG_X_H(DIF(IMJK),DIF(IJK),IM)*ODX_E(IM)*AYZ(IMJK) 
!C----------------------------
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,W,M) = D_F + FLUX_E(IMJK) 
               ELSE 
                  A_M(IJK,W,M) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
               V_F = VF(IJMK) 
!C----------------------------
!C OVERWRITE INDICES TO GHOST CELLS, JCAI
!               D_F = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*ODY_N(JM)*AXZ(IJMK) 
               D_F = AVG_Y_H(DIF(IJMK),DIF(IJK),JM)*ODY_N(JM)*AXZ(IJMK) 
!C----------------------------
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,S,M) = D_F + FLUX_N(IJMK) 
               ELSE 
                  A_M(IJK,S,M) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  KM = KM1(K) 
                  IJKB = BOTTOM_OF(IJK) 
                  V_F = WF(IJKM) 
!C----------------------------
!C OVERWRITE INDICES TO GHOST CELLS, JCAI
!                  D_F = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*OX(I)*ODZ_T(KM)*AXY(&
!                     IJKM) 
                  D_F = AVG_Z_H(DIF(IJKM),DIF(IJK),KM)*OX(I)*ODZ_T(KM)*AXY(&
                     IJKM) 
!C----------------------------
                  IF (V_F >= ZERO) THEN 
                     A_M(IJK,B,M) = D_F + FLUX_T(IJKM) 
                  ELSE 
                     A_M(IJK,B,M) = D_F 
                  ENDIF 
               ENDIF 
            ENDIF 
!
         ENDIF
      END DO 
!
      RETURN  
      END SUBROUTINE CONV_DIF_G0 

end module