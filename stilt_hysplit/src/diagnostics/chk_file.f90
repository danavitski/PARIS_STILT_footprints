      PROGRAM CHK_FILE
!======================================================================
!     full analysis of meteo file structure using standard
!     hysplit4 source code library routines
!----------------------------------------------------------------------
!     Last Revision: 31 Mar 1998 (RRD)
!                    10 May 1999 (RRD) - stop for input
!----------------------------------------------------------------------
! $Id: chk_file.f90,v 1.4 2008-06-03 15:21:54 skoerner Exp $
!=======================================================================

      USE module_defsize
      USE module_defgrid
      IMPLICIT REAL*8 (A-H,O-Z)


      LOGICAL BACK, PRMPT
      CHARACTER HEADER*10000, LABEL*50


!---------------------------------------------------------------------------------------------------
      PRMPT=.FALSE.

!=>check for command line arguments

!     win95: nargs()	unix: iargc()
      NARG = IARGC()
      DO WHILE(NARG > 0)
         IF (NARG > 1) THEN
            WRITE(*,*)'Usage: chk_file [-i{nteractive}]'
            STOP
         ELSE
            CALL GETARG(NARG,LABEL)
            IF (LABEL(1:2) == '-i') PRMPT=.TRUE.
         END IF
         NARG = NARG-1
      END DO

!=>read and initialize data file

      WRITE(*,*)'Enter meteorological directory ...'
      READ(5,'(A)') FILE(1,1)%DIR
      WRITE(*,*)'Enter meteorological file name ...'
      READ(5,'(A)') FILE(1,1)%METEO

      KVEL = 0
      IBMO = 0
      IBDA = 0
      IBHR = 0
      CALL METINI(HEADER,1,0d0,IBYR,IBMO,IBDA,IBHR,.FALSE.,KVEL)

!=>dump file and grid information

      WRITE(*,*)' '
      WRITE(*,'(A20,5I5)')'File start time:',    &
         FILE(1,1)%FIRST%YR,                     &
         FILE(1,1)%FIRST%MO, FILE(1,1)%FIRST%DA, &
         FILE(1,1)%FIRST%HR, FILE(1,1)%FIRST%MN
      WRITE(*,'(A20,5I5)')'File ending time:',   &
         FILE(1,1)%LAST%YR,                      &
         FILE(1,1)%LAST%MO, FILE(1,1)%LAST%DA,   &
         FILE(1,1)%LAST%HR, FILE(1,1)%LAST%MN

      WRITE(*,*)' '
      WRITE(*,'(A20,1X,I0)')'Last Record #:',FILE(1,1)%ENDREC
      WRITE(*,'(A20,1X,I0)')'Record length bytes:',FILE(1,1)%REC_LEN
      WRITE(*,'(A20,1X,A4)')'Meteo data model:',GRID(1)%MODEL_ID
      WRITE(*,'(A20,1X,3I5)')'Grid size x,y,z:', GRID(1)%NX,GRID(1)%NY,GRID(1)%NZ
      WRITE(*,'(A20,1X,I5)')'Vertical coordinate:',DREC(1)%Z_FLAG
      WRITE(*,'(A20,1X,I5)')'First forecast hr:',FILE(1,1)%FIRST%IC
      WRITE(*,'(A20,1X,I5)')'Last forecast hr:',FILE(1,1)%LAST%IC
      WRITE(*,'(A20,1X,I5)')'Records per time:',DREC(1)%REC_PER
      WRITE(*,'(A20,1X,I5)')'Minutes between:',DREC(1)%DELTA

      WRITE(*,*)' '
      WRITE(*,'(A20,2F10.4)')'Pole lat/lon:', &
         GRID(1)%POLE_LAT, GRID(1)%POLE_LON
      WRITE(*,'(A20,2F10.4)')'Reference lat/lon:', &
         GRID(1)%REF_LAT, GRID(1)%REF_LON
      WRITE(*,'(A20,F10.4)')'Grid size (km):',GRID(1)%SIZE
      WRITE(*,'(A20,F10.4)')'Orientation   :',GRID(1)%ORIENT
      WRITE(*,'(A20,F10.4)')'Tang Lat /Cone:',GRID(1)%TANG_LAT
      WRITE(*,'(A20,2F10.4)')'Synch pnt x,y:', &
         GRID(1)%SYNC_XP, GRID(1)%SYNC_YP
      WRITE(*,'(A20,2F10.4)')'Synch pnt lat/lon:', &
         GRID(1)%SYNC_LAT, GRID(1)%SYNC_LON

      WRITE(*,*)' '
      IF (GRID(1)%LATLON) THEN
        CALL GBL2LL(1,1d0,1d0,CLAT,CLON)
      ELSE
        CALL CXY2LL(GRID(1)%GBASE,1d0,1d0,CLAT,CLON)
      END IF
      WRITE(*,'(A20,2F10.4)')'Lower left corner:',CLAT,CLON
      XP = GRID(1)%NX
      YP = GRID(1)%NY
      IF (GRID(1)%LATLON) THEN
        CALL GBL2LL(1,XP,YP,CLAT,CLON)
      ELSE
        CALL CXY2LL(GRID(1)%GBASE,XP,YP,CLAT,CLON)
      END IF
      WRITE(*,'(A20,2F10.4)')'Upper right corner:',CLAT,CLON
      IF (GRID(1)%LATLON) THEN
        CALL GBL2XY(1,CLAT,CLON,XP,YP)
      ELSE
        CALL CLL2XY(GRID(1)%GBASE,CLAT,CLON,XP,YP)
      END IF
      WRITE(*,'(A20,2F10.4)')'Tangent point x,y:',XP,YP

      WRITE(*,*)' '
      IF (PRMPT) THEN
         WRITE(*,*)'Enter to continue ...'
         READ(*,*)
      END IF

!=>dump variable structure

      NLVL = GRID(1)%NZ
      DO L = NLVL,1,-1
         NVAR = DREC(1)%NUM_VARB(L)
         IF (NVAR <= 5) THEN
            WRITE(*, '(1X,F10.4,I3,1X, 5(2X,A4,1X,I5))')            &
               DREC(1)%HEIGHT(L), NVAR, (DREC(1)%VARB_ID(K,L),DREC(1)%CHK_SUM(K,L),K = 1,NVAR)
         ELSE
            WRITE(*, '(1X,F10.4,I3,1X, 5(2X,A4,1X,I5),100(/,15X,5(2X,A4,1X,I5)) )') &
               DREC(1)%HEIGHT(L), NVAR, (DREC(1)%VARB_ID(K,L),DREC(1)%CHK_SUM(K,L),K = 1,NVAR)
         END IF
      END DO

      WRITE(*,*)' '
      IF (PRMPT) THEN
         WRITE(*,*)'Enter to continue ...'
         READ(*,*)
      END IF

!=>dump header field of each data record

      IF (FILE(1,1)%ENDREC == 0) FILE(1,1)%ENDREC = 9999
      IREC = 1
      DO WHILE (IREC <= FILE(1,1)%ENDREC)
         JREC = IREC+DREC(1)%REC_PER-1
         DO KREC = IREC,JREC
            READ(FILE(1,1)%KUNIT,REC = KREC) LABEL
            WRITE(*,*) KREC,'  ',LABEL
         END DO
         IREC = JREC+1

         WRITE(*,*)' '
         IF (PRMPT) THEN
            WRITE(*,*)'Enter to continue ...'
            READ(*,*)
         END IF
      END DO

      END
