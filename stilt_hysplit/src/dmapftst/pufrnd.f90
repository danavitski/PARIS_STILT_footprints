!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFRND           PUFf RouNDs the puff position for sortting
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF ROUND ROUNDS THE PUFF POSITION ACCORDING TO DISPERSION SO
!   IF TWO POSITIONS ARE WITHIN THE FRACTION*SIGMA THEN
!   THEY ARE ASSUMED TO BE IN THE SAME POSITION FOR SORTING
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 - RRD
!
! USAGE:  CALL PUFRND(
!      ZMDL,FRACH,FRACV,FRACT,FMASS,HGD,XPOS1,YPOS1,
!      ZPOS1,SIGH1,SIGV1,HDWP1,PTYP1,MASS1,TMASS1,XPOS2,YPOS2,
!      ZPOS2,SIGH2,SIGV2,HDWP2,PTYP2,MASS2,TMASS2,SORTF,EQUAL)
!   INPUT ARGUMENT LIST:
!     ZMDL        - real model top scaling height
!     FRACH       - real horizontal position rounding fraction
!     FRACV       - real vertical position rounding fraction
!     FRACT       - real travel-time (sigma) rounding fraction
!     FMASS       - real particle mass sort cutoff value
!     HGD         - real horizontal grid distance (m)
!     XPOS,YPOS,ZPOS    - real puff center positions (grid units)
!     SIGH,SIGV         - real horiz (meters) and vert sigma (sigma)
!     HDWP,PTYP         - int  distribution, pollutant type
!     MASS        - real pollutant mass
!     TMASS       - real mass summation over elements
!   OUTPUT ARGUMENT LIST:
!     SORTF       - log  ascending sort flag for indicies 1,2
!     EQUAL       - log  all test variables equal
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

      SUBROUTINE PUFRND(                                                &
     &   ZMDL,FRACH,FRACV,FRACT,FMASS,HGD,XPOS1,YPOS1,                  &
     &   ZPOS1,SIGH1,SIGV1,HDWP1,PTYP1,MASS1,TMASS1,XPOS2,YPOS2,        &
     &   ZPOS2,SIGH2,SIGV2,HDWP2,PTYP2,MASS2,TMASS2,SORTF,EQUAL)

      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
      INCLUDE 'DEFSIZE.INC'

      INTEGER HDWP1,HDWP2,PTYP1,PTYP2
      REAL*8   MASS1(MAXDIM),MASS2(MAXDIM)

!     return flags and equality flag
      LOGICAL SORTF, EQUAL

!     determine any mass sorting options
      IF(FMASS.GT.0.0)THEN
         XMASS1=MASS1(1)
         XMASS2=MASS2(1)

         MM=MAXDIM
         DO WHILE(MM.GT.1)
            XMASS1=XMASS1+MASS1(MM)
            XMASS2=XMASS2+MASS2(MM)
            MM=MM-1
         END DO

!        set mass values to 0 or 1 depending on fmass
         TMASS1=0.0
         TMASS2=0.0
         IF(XMASS1.GT.FMASS)TMASS1=1.0
         IF(XMASS2.GT.FMASS)TMASS2=1.0
      ELSE
         TMASS1=0.0
         TMASS2=0.0
      END IF

!     horizontal sigma precision
      SPREC=FRACT*DMIN1(SIGH1,SIGH2)
!     horizontal position precision
      HPREC=FRACH*DMIN1(SIGH1,SIGH2)/HGD

!     vertical precision depends upon distribution
      VPREC=0
      IF(HDWP1.EQ.HDWP2)THEN
         IF(HDWP1.EQ.1.OR.HDWP1.EQ.2)THEN
!           vertical top-hat
            VPREC=FRACV*DMIN1(SIGV1,SIGV2)
         ELSE
!           vertical particle
            IF(SIGV1.EQ.0.0.OR.SIGV2.EQ.0.0)THEN
!              if either puff has just split then merge
!              for particles assume to be 10. of sigma-h
               VPREC=0.10*DMIN1(SIGH1,SIGH2)/ZMDL
!              with a limit of 10. of the model depth
               VPREC=FRACV*DMIN1(DBLE(0.1),VPREC)
            ELSE
!              no merge until just after horizontal split
               VPREC=0.0
            END IF
         END IF
      END IF

!     ascending order sort rules
      SORTF=.FALSE.
      EQUAL=.FALSE.

!     horizontal distribution type
      IF(HDWP2.LT.HDWP1)THEN
         SORTF=.TRUE.
      ELSEIF(HDWP2.EQ.HDWP1)THEN

!     pollutant type
      IF(PTYP2.LT.PTYP1)THEN
         SORTF=.TRUE.
      ELSEIF(PTYP2.EQ.PTYP1)THEN

!     horizontal sigma same as age
      IF(SIGH2-SIGH1+SPREC.LT.0.0)THEN
         SORTF=.TRUE.
      ELSEIF(ABS(SIGH2-SIGH1).LT.SPREC)THEN

!     particle mass
      IF(TMASS2.LT.TMASS1)THEN
         SORTF=.TRUE.
      ELSEIF(TMASS2.EQ.TMASS1)THEN

!     vertical position
      IF(ZPOS2-ZPOS1+VPREC.LT.0.0)THEN
         SORTF=.TRUE.
      ELSEIF(ABS(ZPOS2-ZPOS1).LT.VPREC)THEN

!     east-west position
      IF(XPOS2-XPOS1+HPREC.LT.0.0)THEN
         SORTF=.TRUE.
      ELSEIF(ABS(XPOS2-XPOS1).LT.HPREC)THEN

!     north-south position
      IF(YPOS2-YPOS1+HPREC.LT.0.0)THEN
         SORTF=.TRUE.
      ELSEIF(ABS(YPOS2-YPOS1).LT.HPREC)THEN
!        total equality
         EQUAL=.TRUE.
         RETURN
      END IF

      END IF
      END IF
      END IF
      END IF
      END IF
      END IF

      RETURN
      END
