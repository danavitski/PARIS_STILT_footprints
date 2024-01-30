!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONINI           CONcentration array INItialization routine
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONCENTRATION ARRAY INITIALIZATION ROUTINE IS A SINGLE CALL
!   TO OPEN CONCENTRATION OUTPUT FILES AND ZERO OUT THE
!   ELEMENTS OF THE CONCENTRATION MATRIX, AND WRITE THE INITIAL
!   INDEX RECORD DATA TO THE CONCENTRATION OUTPUT FILES.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 19 Jun 1997 - RRD
!
! USAGE:  CALL CONINI(NLOC,NUMGRD,NUMTYP,CSUM)
!   INPUT ARGUMENT LIST:
!     NLOC      - int   number of sources
!     NUMGRD    - int   number of concentration grids
!     NUMTYP    - int   number of pollutants
!     CSUM      - real  concentration matrix
!   OUTPUT ARGUMENT LIST:
!     NONE
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     UNITS 21,22,etc according to names set in input CONTROL file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: conini.f90,v 1.4 2007-02-16 17:53:46 tnehrkor Exp $
!
!$$$

      SUBROUTINE CONINI(NLOC,NUMGRD,NUMTYP,CSUM)

      use module_defgrid
      use module_defconc
      use module_defspot

      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     meteorology grid and file
!      INCLUDE 'DEFGRID.INC'
!     pollutant and concentration grid
!      INCLUDE 'DEFCONC.INC'
!     multiple source structure
!      INCLUDE 'DEFSPOT.INC'

!     master concentration array (x,y,z,grids,species)
      REAL*8   CSUM(MAXXP,MAXYP,MAXZP,MAXTYP,MAXGRD)
!     dummy character string array
      CHARACTER LABEL*40

!      COMMON /GBLGRD/ GRID, DREC, FILE
!      COMMON /GBLCON/ CONC, DIRT
!      COMMON /GBLSPT/ SPOT

!     default labels to first defined meteo grid
      KM=1

!     go through each grid
      DO KG=1,NUMGRD

!        set loop indicies
         NXP=CONC(KG)%NUMB_LON
         NYP=CONC(KG)%NUMB_LAT
         NZP=CONC(KG)%LEVELS

!        open the files for output
         KUNIT=20+KG
         CONC(KG)%UNIT=KUNIT
         LABEL=CONC(KG)%DIR
         KLEN=INDEX(LABEL,' ')-1

! JCL:  open conc output file as FORMATTED to generate text file
         OPEN(KUNIT,FILE=LABEL(1:KLEN)//CONC(KG)%FILE,                  &
     &      FORM='FORMATTED',ACCESS='SEQUENTIAL')

!         OPEN(KUNIT,FILE=LABEL(1:KLEN)//CONC(KG)%FILE,
!     :      FORM='UNFORMATTED',ACCESS='SEQUENTIAL')


!        rec#1  meteo file information and number of source locations
         WRITE(KUNIT,*) GRID(KM)%MODEL_ID, FILE(KM,1)%FIRST%YR,         &
     &      FILE(KM,1)%FIRST%MO,FILE(KM,1)%FIRST%DA,                    &
     &      FILE(KM,1)%FIRST%HR,FILE(KM,1)%FIRST%IC,                    &
     &      NLOC

!        rec#2->nloc  source information
!        emission start currently same for all sources
         DO N=1,NLOC
            WRITE(KUNIT,*)                                              &
     &         DIRT(1)%START%YR,  DIRT(1)%START%MO,                     &
     &         DIRT(1)%START%DA,  DIRT(1)%START%HR,                     &
     &         SPOT(N)%OLAT,SPOT(N)%OLON,SPOT(N)%OLVL
         END DO

!        horizontal grid index information
         WRITE(KUNIT,*)                                                 &
     &      CONC(KG)%NUMB_LAT,CONC(KG)%NUMB_LON,CONC(KG)%DELT_LAT,      &
     &      CONC(KG)%DELT_LON,CONC(KG)%X1Y1_LAT,CONC(KG)%X1Y1_LON

!        vertical grid information
         WRITE(KUNIT,*) NZP,(CONC(KG)%HEIGHT(KK),KK=1,NZP)

!        pollutant identification
         WRITE(KUNIT,*) NUMTYP, (DIRT(KK)%IDENT,KK=1,NUMTYP)

!        zero out the array
         DO KT=1,NUMTYP
         DO KL=1,NZP
         DO JJ=1,NYP
         DO II=1,NXP
            CSUM(II,JJ,KL,KT,KG)=0.0
         END DO
         END DO
         END DO
         END DO

      END DO
      RETURN
      END
