!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RRATES0(IER)                                           C
!  Purpose: Calculate reaction rates for various reactions in cell ijk C
!           using information from the data file                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 3-10-98    C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose:Replaced routines with new proceedures for automated        C
!          reaction rate calculations.                                 C
!  Author: J. Musser                                  Date: 10-Oct-12  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX, IJK, T_g, T_s1, D_p, X_g, X_s, EP_g,    C
!            P_g, HOR_g, HOR_s                                         C
!                                                                      C
!                                                                      C
!  Variables modified: M, N, R_gp, R_sp, RoX_gc, RoX_sc, SUM_R_g,      C
!                      SUM_R_s                                         C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DES_RRATES0(NP, pM, IJK, INTERP_IJK, INTERP_WEIGHTS, &
         FOCUS)

      USE compar
      USE constant
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE energy
      USE fldvar
      USE funits
      USE geometry
      USE indices
      USE parallel
      USE param
      USE param1
      USE physprop
      USE run
      USE rxns
      USE sendrecv
      USE usr
      Use parse
      use functions
      use toleranc, only: ZERO_X_gs, COMPARE

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: NP   ! particle index
      INTEGER, INTENT(IN) :: pM   ! Global Solids Phase index
      INTEGER, INTENT(IN) :: IJK  ! fluid cell index
! Variables needed for calculating new interpolation quantities for
! species and energy equations
      INTEGER, INTENT(IN) :: INTERP_IJK(2**3)
      DOUBLE PRECISION, INTENT(IN) :: INTERP_WEIGHTS(2**3)
! Identifies that the indicated particle is of interest for debugging
      LOGICAL, INTENT(IN) :: FOCUS

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: H    ! Reaction loop counter
      INTEGER :: M    ! Global Phase index loop counter
      INTEGER :: N    ! Global species index
      INTEGER :: lN   ! Local reaction speices index/loop counter
      INTEGER :: LM   !

      INTEGER :: mXfr ! Global phase index for mass transfer

! User-defined reaction rates returned from USR_RATES
      DOUBLE PRECISION :: DES_RATES(NO_OF_DES_RXNS)

      DOUBLE PRECISION :: lRate
      DOUBLE PRECISION :: lTp
      DOUBLE PRECISION :: lHoRs

      DOUBLE PRECISION :: RxH
	  DOUBLE PRECISION :: RxHs

! Local gas phase values.
      DOUBLE PRECISION :: lRgp(NMAX(0)) ! Rate of species production
      DOUBLE PRECISION :: lRgc(NMAX(0)) ! Rate of species consumption
      DOUBLE PRECISION :: lHoRg, llHoRg ! Heat of reaction
      DOUBLE PRECISION :: SUMlRg

! Local continuum solid phase values
      DOUBLE PRECISION :: lRsp(NMAX(1)) ! Rate of species production
      DOUBLE PRECISION :: lRsc(NMAX(1)) ! Rate of species consumption
	  DOUBLE PRECISION :: lHoRs_TFM, llHoRs_TFM ! Heat of reaction
      DOUBLE PRECISION :: SUMlRs

! Interphase mass transfer
      DOUBLE PRECISION :: lRPhase(DIMENSION_LM+DIMENSION_M-1)

! Reaction limiters. If a species mass fraction is less than this
! value, then the reaction is suppressed.
      DOUBLE PRECISION :: speciesLimiter
	  DOUBLE PRECISION :: p_vol

! External functions
!---------------------------------------------------------------------//
! Enthalpy calculations (cal/gram)
      DOUBLE PRECISION, EXTERNAL :: CALC_H

! Alias particle temperature.
      lTp = DES_T_s_NEW(NP)
! Particle volume
	  p_vol = 4.0d0/3.0d0*pi*DES_RADIUS(NP)**3
! Initialize storage arrays
      DES_RATES(:) = ZERO
      lRgp(:) = ZERO
      lRgc(:) = ZERO
	  lRsp(:) = ZERO
	  lRsc(:) = ZERO
      lHoRg = ZERO

! Set the species limiter:
      speciesLimiter = ZERO_X_gs

! Calculate user defined reaction rates.
      CALL USR_RATES_DES(NP, pM, IJK, DES_RATES)

! Loop over reactions.
      RXN_LP: DO H = 1, NO_OF_DES_RXNS

! Skip empty reactions
         IF(DES_Reaction(H)%nSpecies == 0) CYCLE RXN_LP
         IF(COMPARE(DES_RATES(H),ZERO)) CYCLE RXN_LP

! Initialize local loop arrays
         llHoRg = ZERO
		 lHoRs = ZERO
		 llHoRs_TFM = ZERO
         RxH = ZERO
		 RxHs = ZERO

! Calculate the rate of formation/consumption for each species.
!---------------------------------------------------------------------//
         DO lN = 1, DES_Reaction(H)%nSpecies
! Global phase index.
            M = DES_Reaction(H)%Species(lN)%pMap
! Global species index.
            N = DES_Reaction(H)%Species(lN)%sMap
! Index for interphase mass transfer. For a gas/solid reaction, the
! index is stored with the gas phase.
            mXfr = DES_Reaction(H)%Species(lN)%mXfr
            lRate = PART_VOL_INTERSEC(IJK,NP)/p_vol * DES_RATES(H) * DES_Reaction(H)%Species(lN)%MWxStoich
! Gas Phase:
            IF(M == 0) THEN
! Consumption of gas phase species.
               IF(lRate < ZERO) THEN
                  IF(X_g(IJK,N) > speciesLimiter) THEN
                     lRgc(N) = lRgc(N) - lRate
! Enthalpy transfer associated with mass transfer. (gas/solid)
                     IF(M /= mXfr) RxH = RxH +                         &
                        lRate*CALC_H(T_g(IJK),0,N)
                  ELSE
! There is an insignificant amount of reactant. Skip this reaction.
                     DES_RATES(H) = ZERO
                     CYCLE RXN_LP
                  ENDIF
               ELSE
! Formation of gas phase species.
                  lRgp(N) = lRgp(N) + lRate
! Enthalpy transfer associated with mass transfer. (gas/solid)
                  IF(M /= mXfr) RxH = RxH + lRate*CALC_H(lTp,0,N)
               ENDIF
! Continuous Solids Phase:
            ELSE IF(M == 1) THEN
               IF(lRate < ZERO) THEN
                  IF(X_s(IJK,M,N) > speciesLimiter) THEN
                     lRsc(N) = lRsc(N) - lRate
					 !IF(M /= mXfr) RxHs = RxHs +                         &
                     !   lRate*CALC_H(T_S(IJK,M),M,N)
					 !This is a hack and may not work for reactions involving
					 !multiple species or reactions where the continuous phase
					 !is the source phase
					 RxHs = RxHs + lRate*CALC_H(T_S(IJK,M),M,N)
                  ELSE
                     DES_RATES(H) = ZERO
                     CYCLE RXN_LP
                  ENDIF
               ELSE
                  lRsp(N) = lRsp(N) + lRate
				  !IF(M /= mXfr) RxHs = RxHs + lRate*CALC_H(lTp,M,N)
				  RxHs = RxHs + lRate*CALC_H(lTp,M,N)
               ENDIF
! Discrete Solids Phase:
            ELSE
! Consumption of solids phase species.
               IF(lRate < ZERO) THEN
                  DES_R_sc(NP,N) = DES_R_sc(NP,N) - lRate
               ELSE
! Formation of solids phase species.
                  DES_R_sp(NP,N) = DES_R_sp(NP,N) + lRate
               ENDIF
            ENDIF
         ENDDO ! Loop of species


! Calculate and store the heat of reaction.
!---------------------------------------------------------------------//
         IF(ENERGY_EQ) THEN
! Automated heat of reaction calculations
            IF(DES_Reaction(H)%Calc_DH) THEN
! Loop over reaction species.
               DO lN = 1, DES_Reaction(H)%nSpecies
! Global phase index.
                  M = DES_Reaction(H)%Species(lN)%pMap
! Global species index.
                  N = DES_Reaction(H)%Species(lN)%sMap
! Rate of formation/consumption for speices N
                  lRate = DES_RATES(H) * &
                     DES_Reaction(H)%Species(lN)%MWxStoich * &
					 PART_VOL_INTERSEC(IJK,NP)/p_vol
! Gas phase enthalpy chnage from energy equation derivation.
                  IF(M == 0) THEN
                     llHORg = llHORg + CALC_H(T_g(IJK),0,N) * lRate
				  ELSE IF(M == 1) THEN
					 llHoRs_TFM = llHoRs_TFM + CALC_H(T_s(IJK,M),M,N) * lRate
! Solid phase enthalpy change from energy equation derivation.
                  ELSE
                     lHORs = lHORs + CALC_H(lTp,M,N) * lRate
                  ENDIF
               ENDDO

! Apply enthalpy transfer associated with mass transfer to get the
! complete heat of reaction for Reaction H.
               llHORg = llHORg - RxH
               lHORs = lHORs + RxH
			   llHoRs_TFM = llHoRs_TFM - RxHs
			   lHORs = lHORs + RxHs

! Convert the heat of reaction to the appropriate units (if SI), and
! store in the global array.
               IF(UNITS == 'SI') THEN
                  lHORg = lHORg + 4.183925d3*llHORg
				  lHoRs_TFM = lHoRs_TFM + 4.183925d3*llHoRs_TFM
                  Q_Source(NP) = Q_Source(NP) - 4.183925d3*lHORs
               ELSE
                  lHORg = lHORg + llHORg
				  lHoRs_TFM = lHoRs_TFM + llHoRs_TFM
                  Q_Source(NP) = Q_Source(NP) - lHORs
               ENDIF
            ELSE
! User-defined heat of reaction.
               HOR_g(IJK) = HOR_g(IJK) + &
                  DES_Reaction(H)%HoR(0) * DES_RATES(H)
			   HOR_s(IJK,1) = HOR_s(IJK,1) + &
			      DES_Reaction(H)%HoR(1) * DES_RATES(H)
               Q_Source(NP) = Q_Source(NP) - &
                  DES_Reaction(H)%HoR(pM) * DES_RATES(H)
            ENDIF
         ENDIF

! Update rate of interphase mass transfer.
!---------------------------------------------------------------------//
		 DO LM=1, (DIMENSION_LM+DIMENSION_M-1)
             lRPhase(LM) = lRPhase(LM) +                       &
                DES_RATES(H) * DES_Reaction(H)%rPHASE(LM)
         ENDDO

         !LM = 1 + (pM - 1)*pM/2
         !lRPhase(LM) = lRPhase(LM) + &
         !   DES_RATES(H) * DES_Reaction(H)%rPHASE(LM)

      ENDDO RXN_LP ! Loop over reactions.

! Calculate the toal rate of formation and consumption for each species.
!---------------------------------------------------------------------//
      IF(SPECIES_EQ(0)) THEN
         SUMlRg = SUM(lRgp(:NMAX(0)) - lRgc(:NMAX(0)))
      ELSE
         DO H=1, NO_OF_DES_RXNS
            IF(DES_Reaction(H)%nPhases <= 0) CYCLE
            LM = 1 + ((pM-1)*pM)/2
            SUMlRg = SUMlRg + &
               DES_RATES(H) * DES_Reaction(H)%rPHASE(LM)
         ENDDO
      ENDIF

      IF(SPECIES_EQ(1)) THEN
         SUMlRs = SUM(lRsp(:NMAX(1)) - lRsc(:NMAX(1)))
      ENDIF

! Integrate over solids time step and store in global array. This
! needs updated when interpolation is reintroduced into thermo code.
!---------------------------------------------------------------------//
      DES_R_gp(IJK,:) = DES_R_gp(IJK,:) + lRgp(:) * DTSOLID
      DES_R_gc(IJK,:) = DES_R_gc(IJK,:) + lRgc(:) * DTSOLID
      DES_R_PHASE(IJK,:) = DES_R_PHASE(IJK,:) + lRPhase(:) * DTSOLID
      DES_HOR_G(IJK) = DES_HOR_G(IJK) + lHoRg * DTSOLID
	  DES_HOR_S(IJK) = DES_HOR_S(IJK) + lHoRs_TFM * DTSOLID
      DES_SUM_R_g(IJK) = DES_SUM_R_g(IJK) + SUMlRg * DTSOLID

      DES_R_sp_TFM(IJK,:) = DES_R_sp_TFM(IJK,:) + lRsp(:) * DTSOLID
      DES_R_sc_TFM(IJK,:) = DES_R_sc_TFM(IJK,:) + lRsc(:) * DTSOLID
      DES_SUM_R_s_TFM(IJK) = DES_SUM_R_s_TFM(IJK) + SUMlRs * DTSOLID

      !print *, SUMlRs
      !print *, DES_SUM_R_s_TFM(:)


      RETURN
      END SUBROUTINE DES_RRATES0
