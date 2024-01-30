!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PAKREC           PAcK a RECord writes one meteo record
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PACK A RECORD WRITES A SPECIFIC RECORD OF DATA FIELD TO UNIT KUN
!   PREVIOUSLY OPENED FOR DIRECT UNFORMATTED ACCESS
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 16 Oct 1997 - RRD
!
! USAGE:  CALL PAKREC(LUNIT,RVAR,CVAR,NX,NY,NXY,KVAR,
!              IY,IM,ID,IH,MN,IC,LL,KINI)
!   INPUT ARGUMENT LIST:
!     LUNIT - int output unit number
!     RVAR  - real      variable of input data to be packed
!     CVAR  - char      packed data array
!     NX,NY - int dimensions of RVAR
!     NXY   - int dimensions of CVAR
!     KVAR  - char      4 character descriptor of variable being written
!     IY,IM,ID    - int date identification
!     IH,MN       - int time identification (MN-minutes)
!     IC    - int forecast hour, use ICX extended hour for >99
!     LL    - int level indicator used to compute record (1 to NLVL)
!     KINI  - int indicator for initialization (0-no 1-yes)
!   OUTPUT ARGUMENT LIST:
!     NONE
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     UNIT as defined by LUNIT
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!$$$

      SUBROUTINE PAKREC(LUNIT,RVAR,CVAR,NX,NY,NXY,KVAR,                 &
     &             IY,IM,ID,IH,MN,IC,LL,KINI)

      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'DEFPACK.INC'

!     real and packed data array
      REAL*8 RVAR(NX,NY)
      CHARACTER CVAR(NXY)*1

!     descriptor of variable in output record
      CHARACTER KVAR*4

!     standard record label
      CHARACTER LABEL*50

!     pass structure to other routines
      COMMON / PAKCOM / GV, NG

!     pass diagnostic variables to other routines
      COMMON / DIAG01 / VAR1, PREC, NEXP

!==>check if requested unit number properly opened

      KG=0
      DO KK=1,NG
         IF(LUNIT.EQ.GV(KK)%KUNIT)KG=KK
      END DO
      IF(KG.EQ.0)THEN
         WRITE(*,*)'ERROR pakrec: Requesting uninitialized unit'
         WRITE(*,*)'Require initial call to PAKSET for unit: ',LUNIT
         STOP
      END IF

!==>test grid dimensions for consistency with initialization

      NXYG=GV(KG)%NXG*GV(KG)%NYG
      IF(GV(KG)%NXG.NE.NX.OR.GV(KG)%NYG.NE.NY.OR.NXYG.NE.NXY)THEN
         WRITE(*,*)'ERROR pakrec: file dimensions do not match'
         WRITE(*,*)'Initialization: ',GV(KG)%NXG,GV(KG)%NYG, NXYG
         WRITE(*,*)'Argument list : ',NX,NY,NXY
         STOP
      END IF

!==>standard forecast hour to write cannot exceed two digits

      ICW=MIN(IC,99)

!==>set all base variables with first entry at each new time

      IF(GV(KG)%NEWT)THEN

!        increment internal counter to next index record
         GV(KG)%MREC=GV(KG)%MREC+GV(KG)%NRPT

!        extended forecast hour (3 digit)
         GV(KG)%ICX=IC

!        save initial times for headers
         GV(KG)%IY0=IY
         GV(KG)%IM0=IM
         GV(KG)%ID0=ID
         GV(KG)%IH0=IH
         GV(KG)%MN0=MN
         GV(KG)%IC0=ICW

!        set switch, reset by pakndx
         GV(KG)%NEWT=.FALSE.

!        initialize all records in this time group to NULL
         IF(KINI.EQ.1)CALL PAKINI(KG,CVAR,NXY)

      ELSE

!==>check current if current record consistent with first

         IF(IY.NE.GV(KG)%IY0.OR.IM.NE.GV(KG)%IM0.OR.                    &
     &      ID.NE.GV(KG)%ID0.OR.IH.NE.GV(KG)%IH0)THEN

            WRITE(*,*)'ERROR pakrec - at index: ',GV(KG)%MREC
            WRITE(*,*)'Argument list times    : ',IY,IM,ID,IH
            WRITE(*,*)'Do not match initial   : ',GV(KG)%IY0,           &
     &         GV(KG)%IM0, GV(KG)%ID0, GV(KG)%IH0
            STOP

         END IF
      END IF

!==>when no data is supplied just do a normal return
!   normally used in conjunction with the initialization flag

      IF(KVAR.EQ.'NULL')RETURN

!==>check vertical index

      IF(LL.LT.1.OR.LL.GT.GV(KG)%NLVL)THEN
         WRITE(*,*)'ERROR pakrec  : Level indicator out of range'
         WRITE(*,*)'Argument level: ',LL
         WRITE(*,*)'Valid Range   : 1 to ',GV(KG)%NLVL
         STOP
      ELSE
!        level indicator should =0 at the surface
         IL=LL-1
      END IF

!==>compute the record offset based upon variable match

      NV=0
      NVAR=GV(KG)%NVAR(LL)
      DO K=1,NVAR
         IF(KVAR.EQ.GV(KG)%VARB(K,LL))NV=K
      END DO

      IF(NV.EQ.0)THEN
         WRITE(*,*)'ERROR pakrec: Variable not in CFG file'
         WRITE(*,*)'Argument list variable: ',KVAR
         STOP
      END IF

!==>pack data and write

!     convert real to packed character
      CALL PAKOUT(RVAR,CVAR,NX,NY,NXY,PREC,NEXP,VAR1,KSUM)

!     save checksum in table
      GV(KG)%CHKS(NV,LL)=KSUM

!     write index portion of record
      WRITE(LABEL,'(5I2,2I2,A4,I4,2E14.7)')                             &
     &   IY,IM,ID,IH,ICW,IL,GV(KG)%IG,KVAR,NEXP,PREC,VAR1

!     compute record based upon variable offset
      JREC=GV(KG)%MREC+GV(KG)%NREC(LL)+NV-1
      IF(JREC.LE.1)THEN
         WRITE(*,*)'ERROR pakrec: output record <=1'
         WRITE(*,*)'Index record: ',GV(KG)%MREC
         WRITE(*,*)'Level offset: ',GV(KG)%NREC(LL)
         WRITE(*,*)'Varbl offset: ',NV
         STOP
      ELSE
         WRITE(GV(KG)%KUNIT,REC=JREC)LABEL,CVAR
      END IF

      RETURN
      END
