!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS1_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms have been calculated but before they are     !
!  applied. The user may insert code in this routine or call user      !
!  defined subroutines.                                                !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR1_DES

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE


! Local variables
!---------------------------------------------------------------------//
! index to track accounted for particles
      INTEGER PC
! dummy index value
      INTEGER L


! Keep the particles stationary.
      !FC(:,:) = ZERO
      !TOW(:,:) = ZERO
      !GRAV(:) = ZERO
	  
	  pc = 1
      do l = 1,max_pip
         if(pc.gt.pip) exit
         if(.not.pea(l,1)) cycle 
         pc = pc+1
         if(pea(l,4)) cycle 
		 Q_Source(l) = Q_Source(l) + 0.0d0
	  end do

      RETURN

      END SUBROUTINE USR1_DES
