!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PAKNDX           PAcK iNDeX writes index record
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PACK INDEX - AFTER ALL THE RECORDS FOR A PARTICULAR TIME
!   PERIOD HAVE BEEN WRITTEN TO A FILE, THIS ROUTINE WRITES THE
!   INDEX RECORD FOR THAT TIME GROUP.  THE INDEX RECORD IS ALWAYS
!   THE FIRST RECORD OF THE TIME GROUP.  IT INCLUDES GRID DEFINITION
!   VARIABLES, AND CHECKSUM INFORMATION.
!
! PROGRAM HISTORY LOG:
!   Last Revised: 14 Feb 1997 - RRD
!
! USAGE:  CALL PAKNDX(LUNIT)
!   INPUT ARGUMENT LIST:
!     LUNIT - int file unit number
!     COMMON PAKCOM
!   OUTPUT ARGUMENT LIST:
!     NONE
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     UNIT as defined in COMMON PAKCOM
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!$$$

      SUBROUTINE PAKNDX(LUNIT)

      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'DEFPACK.INC'

!     standard record label, extender header
      CHARACTER LABEL*50, HEADER*(MLEN)

!     pass structure between routines
      COMMON / PAKCOM / GV, NG

!     pass diagnostic variables to other routines
      COMMON / DIAG01 / VAR1, PREC, NEXP

!==>determine which grid

      KG=0
      DO KK=1,NG
         IF(LUNIT.EQ.GV(KK)%KUNIT)KG=KK
      END DO
      IF(KG.EQ.0)THEN
         WRITE(*,*)'ERROR pakndx: Requesting uninitialized unit'
         STOP
      END IF

!==>conventional 50 byte label

      WRITE(LABEL,'(5I2,2I2,A4,I4,2E14.7)')                             &
     &   GV(KG)%IY0,GV(KG)%IM0,GV(KG)%ID0,GV(KG)%IH0,GV(KG)%IC0,        &
     &   0,GV(KG)%IG,'INDX',0,0.0,0.0

!==>first part of header: 1 -> 108

      WRITE(HEADER(1:108),'(A4,I3,I2,12F7.2,3I3,I2,I4)')                &
     &   GV(KG)%MODEL,GV(KG)%ICX,GV(KG)%MN0,GV(KG)%GRIDS,               &
     &   GV(KG)%NXG,GV(KG)%NYG,GV(KG)%NLVL,GV(KG)%KSYS,GV(KG)%LENH

!==>loop through remainder of the extended header

      KOL=109
      NLVL=GV(KG)%NLVL

      DO L=1,NLVL
         ZL=GV(KG)%HEIGHT(L)

!        precision depends upon the height coordinate
         IF(ZL.GE.10000.0)THEN
            WRITE(HEADER(KOL:KOL+7),'(F6.0,I2)')ZL,GV(KG)%NVAR(L)
         ELSEIF(ZL.GE.1000.0)THEN
            WRITE(HEADER(KOL:KOL+7),'(F6.1,I2)')ZL,GV(KG)%NVAR(L)
         ELSEIF(ZL.GE.100.0.AND.ZL.LT.1000.0)THEN
            WRITE(HEADER(KOL:KOL+7),'(F6.2,I2)')ZL,GV(KG)%NVAR(L)
         ELSEIF(ZL.GE.10.0.AND.ZL.LT.100.0)THEN
            WRITE(HEADER(KOL:KOL+7),'(F6.3,I2)')ZL,GV(KG)%NVAR(L)
         ELSEIF(ZL.GE.1.0.AND.ZL.LT.10.0)THEN
            WRITE(HEADER(KOL:KOL+7),'(F6.4,I2)')ZL,GV(KG)%NVAR(L)
         ELSE
            WRITE(HEADER(KOL:KOL+7),'(F6.5,I2)')ZL,GV(KG)%NVAR(L)
         END IF

!        add variable id's and checksums
         KOL=KOL+8
         NVAR=GV(KG)%NVAR(L)
         DO K=1,NVAR
            WRITE(HEADER(KOL:KOL+7),'(A4,I3)')                          &
     &         GV(KG)%VARB(K,L), GV(KG)%CHKS(K,L)
            KOL=KOL+8
         END DO
      END DO

!==>write extended header to disk

      NHL1=1
!     number of index records
      NREC=GV(KG)%NHREC
!     point to first index record
      JREC=GV(KG)%MREC

      DO K=1,NREC
!        byte count for each index
         NHL2=NHL1+GV(KG)%LREC-1
         IF(K.EQ.NREC)NHL2=NHL1+GV(KG)%NHBYT-1

         WRITE(GV(KG)%KUNIT,REC=JREC)LABEL,HEADER(NHL1:NHL2)
         JREC=JREC+1
         NHL1=NHL2+1
      END DO

!==>clear flags

!     checksum table
      DO J=1,MLVL
      DO I=1,MVAR
         GV(KG)%CHKS(I,J)=0
      END DO
      END DO

!     new time flag
      GV(KG)%NEWT=.TRUE.

      RETURN
      END
