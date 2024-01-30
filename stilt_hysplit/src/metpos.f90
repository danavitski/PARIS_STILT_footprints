!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METPOS           METeorological POSitioning finds record 
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL POSITIONING ROUTINE USES THE CURRENT CALCULATION 
!   POSITION AND TIME TO CHECK IF THE POSITION FALLS WITHIN THE METEO
!   DATA ALREADY LOADED INTO MEMORY. IF NOT IS IS DETERMINED IF DATA
!   IN A NEW LOCATION ARE REQUIRED (SUBGRID MOVES) OR AT A NEW TIME,
!   BOTH.  THE ROUTINE RETURNS A POSITIVE RECORD NUMBER IF DATA ARE 
!   REQUIRED TO BE READ FROM THE INPUT FILE.                        
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 15 May 1998 (RRD)                                
!                  08 Mar 1999 (RRD) - optimized subgrid selection
!                                    - removed grid number dependence
!                  17 Jun 1999 (RRD) - corrected newx flag at first entry
!                  08 Nov 1999 (RRD) - do not set newx at times between files
!
! USAGE:  CALL METPOS(BACK,XP,YP,JET,KG,NXS,NYS,LX1,LY1,LXC,LYC,
!            LXR,LYR,KREC1,KREC2,OFFG,NEWX,MTIME,FTIME,KGRID) 
!   INPUT ARGUMENT LIST:
!     BACK      - log   define backward integration                          
!     XP        - real  x grid position                                       
!     YP        - real  y grid position                                       
!     JET       - int   elapsed time (minutes from 1 jan)                     
!     KG        - int   grid number being loaded                               
!     KGRID     - int   current open meteorological grid index              
!     LXC,LYC   - int   subgrid center position   
!     LXR,LYR   - int   maximum subgrid size
!   OUTPUT ARGUMENT LIST:
!     LX1,LY1   - int   lower left corner of meteo sub-grid               
!     NXS,NYS   - int   grid limits of dynamic meteo sub-grid             
!     KREC1     - int   record number of data at last time (1)              
!     KREC2     - int   record number of data at next time (2)          
!     OFFG      - log   off entire grid flag                                 
!     NEWX      - log   off sub-grid (read more data)                        
!     MTIME     - int   accumulated min of requested input data             
!     FTIME     - int   accumulated min of current meteo data in array      
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     UNIT 30 diagnostic ouput MESSAGE file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!
! $Id: metpos.f90,v 1.11 2007-05-03 13:09:12 skoerner Exp $
!
!$$$

SUBROUTINE METPOS (BACK,XP,YP,JET,KG,NGRD,NXS,NYS,LX1,LY1,LXC,LYC, &
                   LXR,LYR,KREC1,KREC2,OFFG,NEWX,MTIME,FTIME,KGRID)

      USE module_defgrid
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER, INTENT(IN) :: NGRD

!     file request and position in minutes for time 1 and time 2
      INTEGER   MTIME(2), FTIME(2)

!     flags for new position, time, direction, offgrid
      LOGICAL NEWX, NEWT, BACK, OFFG


!---------------------------------------------------------------------------------------------------
!     start default assumption that advection point is on the grid
      OFFG = .FALSE.

!     convert position to index units
      II = INT(XP)
      JJ = INT(YP)

!     zero return indicates no data read required
!     first time or sub-grid changes - both times must be read
!     when time changes only new time (krec2) requires input
      KREC1 = 0
      KREC2 = 0

!     flags for data required at new time or new sub-grid        
      NEWT = .FALSE.
      NEWX = .FALSE.
     
!     truncate elapsed time to previous data time
!     only upon first entry (ie model times not yet defined)
      IF (FTIME(2) == 0 .OR. KGRID == 0) THEN
         MTIME(2) = INT(JET/DREC(KG)%DELTA)*DREC(KG)%DELTA
         IF (BACK .AND. MOD(JET,DREC(KG)%DELTA) /= 0)                      &
            MTIME(2) = MTIME(2)+DREC(KG)%DELTA
!        set temporal flag
         NEWT = .TRUE.
      END IF

!     at first entry or when grid changes usually need to read both     
!     times and a new sub-grid
      IF (KG > KGRID) THEN
!        negative values indicates data read for time period 1 or 2
!        the new period (#2) is always read
         KREC2=-1          

!        if meteo grid is the same do not re-read old time period
!        as the files are probably continuous in time, a change in grid
!        requires old time data as well as new sub-grid. The ftime test
!        checks to see if the old file was actually read (setup error)
         IF (FTIME(1) == 0 .OR. KGRID == 0) THEN
            KREC1=-1
            NEWX = .TRUE.
         ELSEIF((GRID(KGRID)%MODEL_ID) /= (GRID(KG)%MODEL_ID) .OR.       &
                (GRID(KGRID)%SIZE) /= (GRID(KG)%SIZE)) THEN
            KREC1=-1
            NEWX = .TRUE.
!           changes in temporal interval require recomputation 
!           of the data request time, ie data may not be available
!           at the time indicated using the old parameters 
            IF (DREC(KGRID)%DELTA /= DREC(KG)%DELTA) THEN
               MTIME(2) = INT(JET/DREC(KG)%DELTA)*DREC(KG)%DELTA
               IF (BACK .AND. MOD(JET,DREC(KG)%DELTA) /= 0)                &
                  MTIME(2) = MTIME(2)+DREC(KG)%DELTA
            END IF
         END IF

!        set the current grid index number only the initial time
         IF (KGRID == 0) KGRID = KG

         WRITE (30,*) 'NOTICE metpos: initial data load for grid -',KG
         WRITE (30,*) 'Request time:',MTIME(1),MTIME(2)
         WRITE (30,*) 'Loaded times:',FTIME(1),FTIME(2)

      
!     current grid matches calculation point
      ELSEIF(KG == KGRID) THEN
         IF (.NOT.BACK .AND. JET > FTIME(2) .OR.                           &
             .NOT.BACK .AND. JET < FTIME(1) .OR.                           &
                  BACK .AND. JET < FTIME(2) .OR.                           &
                  BACK .AND. JET > FTIME(1)) THEN

!           calculations out of sequence
            WRITE (30,*) 'ERROR metpos: positioning out of sequence'
            WRITE (30,*) 'Current time:',JET
            WRITE (30,*) 'File   times:',FTIME(1),FTIME(2)
            WRITE (30,*) 'Request time:',MTIME(1),MTIME(2)
            STOP
         END IF

!     calculation point still on old grid, calling routine will switch
      ELSE
         OFFG = .TRUE.
         RETURN
      END IF

!     determine if data are required at a new time
      IF (JET == FTIME(2)) NEWT = .TRUE.
      IF (MTIME(2) /= FTIME(2)) THEN
         NEWT = .TRUE.
         KREC2=-1
      END IF

!     determine if new data are required because position
!     has moved off the sub-grid, a change in files, on the same grid
!     can lead to trouble if subgrid data are required at the old time
!     and there is no temporal overlap of the files  
      max_offset = 2
      if (GRID(KG)%MODEL_ID(2:4) .eq. 'WRF') max_offset = max_offset+1
      IF (II <= LX1+1 .OR. II >= (LX1+NXS-max_offset) .OR.                           &
         JJ <= LY1+1 .OR. JJ >= (LY1+NYS-max_offset)) THEN
!        only set flag if calculation not between grids
!        avoid reading subgrid at time 1 when file no longer available
         IF (.NOT.(KG > 1 .AND. MTIME(2) == FILE(KG,2)%FIRST%MACC))       &
            NEWX = .TRUE.

! JCL:(02/18/2005)
!        subgrid shifts not required for global grid
         IF (GRID(KG)%GLOBAL .AND. FTIME(2) /= 0) NEWX = .FALSE.

      END IF

!     new data required
      IF (NEWX .OR. NEWT) THEN

!        flip indicies at first time (ftime=0) or when the request   
!        time equals the file time
         IF (FTIME(2) == 0 .OR. (NEWT .AND. MTIME(2) == FTIME(2))) THEN
!           reverse index positions each time to avoid
!           copying data from new time to old time array
            FILE(KG,1)%PERIOD = MOD(FILE(KG,1)%PERIOD,2)+1
            FILE(KG,2)%PERIOD = MOD(FILE(KG,2)%PERIOD,2)+1

                          
!           index position 2 indicates new time
            MTIME(1) = MTIME(2)
            FTIME(1) = FTIME(2)
            IF (.NOT. BACK) THEN
               MTIME(2) = MTIME(1)+DREC(KG)%DELTA
            ELSE
               MTIME(2) = MTIME(1)-DREC(KG)%DELTA
            ENDIF

!           set new data read indicator
            KREC2=-1
         END IF

         IF (NEWX) THEN
! JCL:(02/18/2005)
            IF (.NOT.GRID(KG)%GLOBAL) THEN
!           check if position off grid entirely
            max_offset = 1
            if (GRID(KG)%MODEL_ID(2:4) .eq. 'WRF') max_offset = max_offset+1
            IF (II <= 2 .OR. II >= GRID(KG)%NX-max_offset .OR.                       &
                JJ <= 2 .OR. JJ >= GRID(KG)%NY-max_offset) THEN
               OFFG = .TRUE.
               NEWX = .FALSE.
!              only continue in this routine if data required
!              for a new time
               IF (.NOT.NEWT) RETURN
!              staying in routine set flag back (6/17/99)
               NEWX = .TRUE.  
            END IF 
            END IF

!           check if subgrid center position set, if not use default
            IF (LXC <= 0 .OR. LYC <= 0) THEN
!              compute new sub-grid corner based upon current position
               LX1 = MIN(MAX(1,II-NXM/2),GRID(KG)%NX)
               LY1 = MIN(MAX(1,JJ-NYM/2),GRID(KG)%NY)

            ELSE
!              subgrid optimization values set in main call advrng
!              set new lower left corner for sub-grid based on center
!              position of all computational points from last hour
               LX1 = MIN(MAX(1,LXC-LXR/2),GRID(KG)%NX)
               LY1 = MIN(MAX(1,LYC-LYR/2),GRID(KG)%NY)
            END IF

!           check upper right corner against maximum limits
!           because lower left at 1,1 -> upper right is size
            NXS = MIN(LX1+LXR-1,GRID(KG)%NX)-LX1+1
            NYS = MIN(LY1+LYR-1,GRID(KG)%NY)-LY1+1
            IF (NXS > NXM .OR. NYS > NYM) THEN              ! bounds exeeded
               PRINT 100, NXS, NYS, NXM, NYM, MAXVAL(GRID(1:NGRD)%NX), MAXVAL(GRID(1:NGRD)%NY)
               STOP
            ENDIF

!           set new data read indicator
!           subgrid shift requires new data at old and new time
            KREC1=-1
            KREC2=-1
         END IF

!        compute record numbers to index record for last time
         IF (KREC1 < 0) THEN
            KREC1 = DREC(KG)%REC_PER*(MTIME(1)-FILE(KG,1)%FIRST%MACC)     &
            /DREC(KG)%DELTA+1

            
!           are the desired records within the file
            IF (BACK .AND. KREC1 > (FILE(KG,1)%ENDREC)) THEN
               WRITE (*,*) 'ERROR metpos: Positioning beyond EOF'
               IF (FTIME(1) == 0 .AND. FTIME(2) == 0) THEN
                  WRITE (*,*) 'Model starting time between files'
               ELSE
                  WRITE (*,*) 'Between files - compile larger subgrid'
               END IF
               STOP
            END IF
           
            IF (.NOT.BACK .AND. KREC1 < 1) THEN
               WRITE (*,*) 'ERROR metpos: Positioning before file start'
               IF (FTIME(1) == 0 .AND. FTIME(2) == 0) THEN
                  WRITE (*,*) 'Model starting time between files'
               ELSE
                  WRITE (*,*) 'Between files - compile larger subgrid'
               END IF
               STOP
            END IF
         END IF

!        record numbers for next time require file end check
         IF (KREC2 < 0) THEN
            KREC2 = DREC(KG)%REC_PER*(MTIME(2)-FILE(KG,2)%FIRST%MACC)     &
            /DREC(KG)%DELTA+1

!           are the desired records within the file
            IF ( (BACK .AND. KREC2 < 1) .OR.                              &
                (.NOT.BACK .AND. KREC2 > (FILE(KG,2)%ENDREC))) THEN       
                OFFG = .TRUE.
                RETURN
            END IF

         END IF
      END IF

!---------------------------------------------------------------------------------------------------
100   FORMAT('ERROR metpos: subgrid dimension(s) NXS and/or NYS exceed(s) compiled limit:'/ &
             '   NXS = ', i0, ', NYS = ', i0/                                               &
             '   NXM = ', i0, ', NYM = ', i0/                                               &
             'Recompile after increasing NXM and/or NYM in module_defsize.'/                &
             'Best is to set it to max. met-file extents, NX=', i0, ', NY=', i0, '.'/       &
             'Stop.')


END SUBROUTINE METPOS
