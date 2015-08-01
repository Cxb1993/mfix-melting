!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: COMP_MEAN_FIELDS                                        !
!  Author: J.Musser                                   Date: 11-NOV-14  !
!                                                                      !
!  Purpose: Driver routine for calculating field variables (DES_ROP_s, !
!   DES_U_S, DES_V_S, DES_W_S) from particle data.                     !
!                                                                      !
!  o The diffusion filter is only applied to the the solids bulk       !
!    density because DEM simulations do not utilize the other field    !
!    variables within a time loop.                                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE COMP_MEAN_FIELDS

      use particle_filter, only: DES_INTERP_MEAN_FIELDS
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use particle_filter, only: DES_INTERP_GARG

      use discretelement, only: DES_MMAX
      use discretelement, only: DES_ROP_S

! Flag: Diffuse DES field variables.
      use particle_filter, only: DES_DIFFUSE_MEAN_FIELDS

      IMPLICIT NONE

! Loop counter.
      INTEGER :: M

!......................................................................!

! Calculate field variables from particle data:
      IF(DES_INTERP_MEAN_FIELDS) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE(DES_INTERP_NONE) ; CALL COMP_MEAN_FIELDS_ZERO_ORDER
         CASE(DES_INTERP_GARG) ; CALL COMP_MEAN_FIELDS0
         CASE DEFAULT; CALL COMP_MEAN_FIELDS1
         END SELECT
      ELSE
         CALL COMP_MEAN_FIELDS_ZERO_ORDER
      ENDIF

! Apply the diffusion filter.
      IF(DES_DIFFUSE_MEAN_FIELDS) THEN
         DO M=1, DES_MMAX
            CALL DIFFUSE_MEAN_FIELD(DES_ROP_S(:,M),'DES_ROP_S')
         ENDDO
      ENDIF

! Calculate the gas phase volume fraction from DES_ROP_s.
      CALL CALC_EPG_DES

      RETURN
      END SUBROUTINE COMP_MEAN_FIELDS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      use desgrid
      use desmpi
      USE mfix_pic
      USE functions
      use indices, only: I_OF, J_OF, K_OF

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters: partciles, filter cells, phases
      INTEGER NP, M
! Fluid cell index
      INTEGER IJK
! Total Mth solids phase volume in IJK
      DOUBLE PRECISION :: SOLVOLINC(DIMENSION_3,DES_MMAX)
! One divided by the total solids volume.
      DOUBLE PRECISION :: OoSOLVOL
! PVOL times statistical weight
      DOUBLE PRECISION :: VOL_WT, VOL_INT
      INTEGER :: I, J, K, IJK_INDEX
!Min and max extents for cell and particle center
      DOUBLE PRECISION, DIMENSION(3):: BMIN,BMAX,X_C
      DOUBLE PRECISION :: VOL_TOL
      DOUBLE PRECISION :: VOL_SUM
      DOUBLE PRECISION :: FRAC_INT

      VOL_TOL = 0.0001d0
      SOLVOLINC(:,:) = ZERO
	  
      PART_VOL_INTERSEC(:,:) = ZERO
      TOT_VOL_INTERSEC(:) = ZERO
      PART_CELLS(:,:) = ZERO

      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO

      DO NP=1,MAX_PIP
         IF(.NOT.PEA(NP,1)) CYCLE
         IF(PEA(NP,4)) CYCLE
         M = PIJK(NP,5)
         X_C(:) = DES_POS_NEW(:,NP)
      
         DO IJK = IJKSTART3, IJKEND3
            IF(.NOT.FLUID_AT(IJK)) CYCLE
		 
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
		 
            BMIN = (/XE(I-1),YN(J-1),ZT(K-1)/)
            BMAX = (/XE(I),YN(J),ZT(K)/)
		 
            CALL VOL_INTERSECTION(BMIN-X_C,BMAX-X_C,DES_RADIUS(NP), & 
              VOL(IJK)*VOL_TOL,VOL_INT)
            IF(VOL_INT > 0.0d0) THEN
               PART_CELLS(NP,1) = PART_CELLS(NP,1) + 1
               PART_CELLS(NP,PART_CELLS(NP,1)+1) = IJK
               PART_VOL_INTERSEC(IJK,NP)=PART_VOL_INTERSEC(IJK,NP)+VOL_INT
            END IF 
         ENDDO
         VOL_SUM = 0.0d0
!Normalize volume fractions so they sum to total volume for each particle
         DO IJK_INDEX = 2,PART_CELLS(NP,1)+1
            IJK = PART_CELLS(NP,IJK_INDEX)
            VOL_SUM = VOL_SUM + PART_VOL_INTERSEC(IJK,NP)
         END DO
         DO IJK_INDEX = 2,PART_CELLS(NP,1)+1
            IJK = PART_CELLS(NP,IJK_INDEX)
            FRAC_INT = PART_VOL_INTERSEC(IJK,NP)/VOL_SUM
            VOL_INT = FRAC_INT*PVOL(NP)
            PART_VOL_INTERSEC(IJK,NP) = VOL_INT
            TOT_VOL_INTERSEC(IJK) = TOT_VOL_INTERSEC(IJK)+VOL_INT
            SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M)+VOL_INT

            DES_U_S(IJK,M) = DES_U_S(IJK,M) +                         &
              DES_VEL_NEW(1,NP)*VOL_INT
            DES_V_S(IJK,M) = DES_V_S(IJK,M) +                         &
              DES_VEL_NEW(2,NP)*VOL_INT
            DES_W_S(IJK,M) = DES_W_S(IJK,M) +                         &
              DES_VEL_NEW(3,NP)*VOL_INT
         END DO
      END DO
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         DO M = 1, DES_MMAX
            IF(SOLVOLINC(IJK,M).GT.ZERO) THEN
               OoSOLVOL = ONE/SOLVOLINC(IJK,M)
               DES_U_s(IJK,M) = DES_U_s(IJK,M)*OoSOLVOL
               DES_V_s(IJK,M) = DES_V_s(IJK,M)*OoSOLVOL
               DES_W_s(IJK,M) = DES_W_s(IJK,M)*OoSOLVOL
            ENDIF

! calculating the bulk density of solids phase m based on the total
! number of particles having their center in the cell
            DES_ROP_S(IJK,M) = DES_RO_S(M)*SOLVOLINC(IJK,M)/VOL(IJK)

         ENDDO
		 
      ENDDO

! Halo exchange of solids volume fraction data.
      CALL SEND_RECV(DES_ROP_S,2)

      RETURN
      END SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER
	  
      RECURSIVE SUBROUTINE VOL_INTERSECTION(BMIN,BMAX,R,TOL,VOL)
      !D.Moser 07/21/2014
      !Subroutine to calculate the volume of intersection between
      !a sphere centered at 0,0,0 with radius R and an AAB with max
      !and min coordinates given by BMAX and BMIN. For spheres not
      !centered at 0,0,0, subtract the coordinates of the sphere center
      !from BMAX and BMIN before passing to this function

      !Minimum and maximum coordinates of box
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: BMIN,BMAX
      !Sphere radius and minimum volume to subdivide
      DOUBLE PRECISION, INTENT(IN) :: R,TOL
      !Output volume of intersection
      DOUBLE PRECISION, INTENT(OUT) :: VOL

      !Local variables
      !Volume of box
      DOUBLE PRECISION BVOL
      !Arrays to store all x,y, and z coordinates for easy traversing
      DOUBLE PRECISION, DIMENSION(2) :: BX, BY, BZ
      !Coordinates of box midpoint and new min/max for recusive calling
      DOUBLE PRECISION, DIMENSION(3) :: MIDPT,BMIN_NEW,BMAX_NEW
      !Variable for checking if sphere and box intersect
      DOUBLE PRECISION DMIN
      !Temp volume variable for recursion
      DOUBLE PRECISION TMP_VOL
      !Indicies for loops
      INTEGER I,J,K
      !Pi
      DOUBLE PRECISION PI
      !Logical for breaking a loop
      LOGICAL BREAK

!set up variables
      bvol = (bmax(1)-bmin(1))*(bmax(2)-bmin(2))*(bmax(3)-bmin(3))
      midpt(:) = bmin(:) + (bmax(:)-bmin(:))/2.0d0
      pi = 3.14159265359d0
!volume is already smaller than prescribed by tolerance
      IF(bvol < tol) THEN
         IF(DOT_PRODUCT(midpt,midpt) <= r*r) THEN
            vol = bvol
            RETURN
         ELSE
            vol = 0.0d0
            RETURN
         END IF
      END IF
!sphere is entirely outside of box
      dmin = 0.0d0
      DO i=1,3
        IF (0.0d0 < bmin(i)) THEN
           dmin = dmin + (bmin(i))*(bmin(i))
        ELSE IF (0.0d0 > bmax(i)) THEN
           dmin = dmin + (bmax(i))*(bmax(i))
        END IF
      END DO
      IF (dmin >= r*r) THEN
         vol = 0.0d0
         RETURN
      END IF
!sphere is entirely inside box
      IF (ALL(r <= bmax(:)).AND.ALL(-r >= bmin(:))) THEN
         vol =  4.0d0/3.0d0*pi*r*r*r
         RETURN
      END IF
!box is entirely inside sphere
      bx = (/bmin(1),bmax(1)/)
      by = (/bmin(2),bmax(2)/)
      bz = (/bmin(3),bmax(3)/)
      IF (ALL(bmax(:) <= r) .AND. ALL(bmin(:) >= -r)) THEN
         break = .FALSE.
         DO i=1,2
            DO j=1,2
               DO k=1,2
                  IF (bx(i)*bx(i)+by(j)*by(j)+bz(k)*bz(k) > r*r) THEN
                     break = .TRUE.
                     EXIT
                  END IF
               END DO
               IF(break) EXIT
            END DO
            IF(break) EXIT
         END DO
         IF(.NOT.break) THEN
            vol = bvol
            RETURN
         END IF
      END IF
!box contains partial intersection, bisect and call recursively
      vol = 0.0d0
      DO i=1,2
         DO j=1,2
            DO k=1,2
               bmin_new = (/MIN(midpt(1),bx(i)),MIN(midpt(2),by(j)),&
               MIN(midpt(3),bz(k))/)
               bmax_new = (/MAX(midpt(1),bx(i)),MAX(midpt(2),by(j)),&
               MAX(midpt(3),bz(k))/)
               CALL VOL_INTERSECTION(BMIN_NEW,BMAX_NEW,R,TOL,TMP_VOL)
               vol = vol + tmp_vol
            END DO
         END DO
      END DO
      RETURN
      END SUBROUTINE VOL_INTERSECTION
