!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PAKOUT           PAcK OUTput converts real to character
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PACK OUTPUT CONVERTS A REAL*8 ARRAY TO CHARACTER*1 PACKED ARRAY
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 14 Feb 1997 - RRD
!
! USAGE:  CALL PAKOUT(RVAR,CVAR,NX,NY,NXY,PREC,NEXP,VAR1,KSUM)
!   INPUT ARGUMENT LIST:
!     RVAR  - real      input real array of dimensions NX,NY
!   OUTPUT ARGUMENT LIST:
!     CVAR  - char      packed char*1 output array of length NX*NY
!     PREC  - real      precision of packed data array
!     NEXP  - int packing scaling exponent
!     VAR1  - real      value of real array at position (1,1)
!     KSUM  - int rotating int checksum of packed data array
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!$$$

      SUBROUTINE PAKOUT(RVAR,CVAR,NX,NY,NXY,PREC,NEXP,VAR1,KSUM)

      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 RVAR(NX,NY)
      CHARACTER CVAR(NXY)

      VAR1=RVAR(1,1)

      ROLD= VAR1
      RMAX= 0.0
!     find the maximum difference between adjacent elements
      DO J=1,NY
         DO I=1,NX
!           compute max difference between elements along row
            RMAX=DMAX1(ABS(RVAR(I,J)-ROLD), RMAX)
            ROLD=RVAR(I,J)
         END DO
!        row element 1 difference always from previous row
         ROLD=RVAR(1,J)
      END DO

      SEXP=0.0
!     compute the required scaling exponent
      IF(RMAX.NE.0.0)SEXP=DLOG(RMAX)/DLOG(DBLE(2.))
      NEXP=INT(SEXP)
!     positive or whole number scaling round up for lower precision
      IF(SEXP.GE.0.0.OR.DMOD(SEXP,DBLE(1.0)).EQ.0.0)NEXP=NEXP+1
!     precision range is -127 to 127 or 254
      PREC=(2.0**NEXP)/254.0
      SCALE=2.0**(7-NEXP)

!     initialize checksum
      KSUM=0
!     set column1 value
      RCOL=VAR1

      K=0
!     pack the array from low to high
      DO J=1,NY
         ROLD=RCOL
         DO I=1,NX
            K=K+1

!           packed integer at element
            ICVAL=(RVAR(I,J)-ROLD)*SCALE+127.5
!           previous element as it would appear unpacked
            ROLD=(ICVAL-127.0)/SCALE+ROLD
!           save the first column element for next row
            IF(I.EQ.1)RCOL=ROLD
!           convert to character
            CVAR(K)=CHAR(ICVAL)

!           maintain rotating checksum
            KSUM=KSUM+ICVAL
!           if sum carries over the eighth bit add one
            IF(KSUM.GE.256)KSUM=KSUM-255

         END DO
      END DO

      RETURN
      END
