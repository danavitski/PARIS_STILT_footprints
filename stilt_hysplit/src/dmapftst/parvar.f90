!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PARVAR           PARticle VARiance returns random component
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PARTICLE VARIANCE RETURNS RANDOM COMPONENT GIVEN VELOCITY SIGMA.
!   RESULTS WITH A GAUSSIAN DISTRIBUTION AND A STANDARD
!   DEVIATION OF THE INPUT SIGMA AND A MEAN OF ZERO.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 19 Dec 1998 (RRD) - replace iterative solution with table
!                 20 Jan 1999 (MDC) - fixed array bound test failure
!
! USAGE:  CALL PARVAR(SIGMA,VALOUT)
!   INPUT ARGUMENT LIST:
!     SIGMA     - real  input desired standard deviation
!   OUTPUT ARGUMENT LIST:
!     VALOUT    - real  output value
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY (see notes in code for platform conversion)
!
!$$$

! JCL:add 'RSEED', the random seed for random seed generator, as additional argument
      SUBROUTINE PARVAR(SIGMA,VALOUT,RSEED)
!      SUBROUTINE PARVAR(SIGMA,VALOUT)

! JCL:(5/8/00)have NSIG value that is 2 orders of magnitude larger than original
!        value of 200
!     range of sigma values over which values are returned (0-3)
!      PARAMETER (NSIG=30)
      PARAMETER (NSIG=3000)

! JCL:include the FSUBLIB.FI include file in order
!    to use GETTIM to get current time for random seed
!      INCLUDE 'FSUBLIB.FI'

! JCL:the random seed that was determined in HYMODELC
      INTEGER RSEED

! jcl:arguments to GETTIM--HSECS will be used for seed
!      INTEGER*2 HRS,MINS,SECS,HSECS

!_sgi REAL*8 XVAL,RAND
!ray  REAL*8   XVAL,RANF
! jcl:not need RAND in revised version
      REAL*8 XVAL
!      REAL*8   XVAL,RAND


      LOGICAL NEGV
! JCL:reduce 1-D array by 1
!     Gaussian cumulative distribution function
      REAL*8 CDF(0:NSIG)
!      REAL*8 CDF(0:NSIG+1)
!     square root of 2 pi
      DATA SQR2PI/2.5066283/
!     seed initialization if required
      DATA KVAL/0/

      SAVE SQR2PI,CDF

!     first entry into subroutine compute distribution
      IF(SQR2PI.GT.0.0)THEN
!        initialize positive side of distribution
         CDF(0)=0.5
         XSUM=CDF(0)*SQR2PI

!        integrate Gaussian from 0 to NSIG/10 sigma at 0.1 increments
         DO KNDX=1,NSIG
! JCL:(5/8/00)make integration step 100 times smaller to reduce
!           numerical error
            XVAL=KNDX/1000.0
!            XVAL=KNDX/10.0
! JCL:(5/8/00)also make integration step 100 times smaller
            XSUM=XSUM+0.001*EXP(-0.5*XVAL*XVAL)
!            XSUM=XSUM+0.1*EXP(-0.5*XVAL*XVAL)
            CDF(KNDX)=XSUM/SQR2PI
! JCL:
!            WRITE(45,*)KNDX,CDF(KNDX)
         END DO
!        set value to zero to skip this code section
         SQR2PI=0.0
      END IF


!     compute random value
!_PC      CALL RANDOM(RVAL)
!ray  RVAL=RANF()
!_Sun RVAL=RAND(KVAL)
!      RVAL=RAND()
! jcl:since Watcom doesn't support the 'RAND' function,
!       use the URAND function instead.  Use the current time
!       in hundredths of seconds for seed in URAND
! jcl:call GETTIM to get current time in hundredths of secs
!      CALL GETTIM(HRS, MINS, SECS, HSECS)
!      RVAL=URAND(HSECS)
! JCL:call RAN3, an algorithm from Press et al., to generate random number
!     use RSEED as random seed for random number generator
      RVAL=RAN3(RSEED)


!     check for negative side of distribution
      IF(RVAL.LT.0.5)THEN
         NEGV=.TRUE.
!         XVAL=RVAL+0.5
! JCL:   new XVAL to get rid of jumps
         XVAL=1-RVAL
      ELSE
         NEGV=.FALSE.
         XVAL=RVAL
      END IF

!     we assume that the random number generator returns
!     a value of the cumulative distribution (0->1)
!     then the table lookup returns the sigma that equals the CDF
      KNDX=0
      DO WHILE (CDF(KNDX).LT.XVAL.AND.KNDX.LE.NSIG)
         KNDX=KNDX+1
      END DO

!     the sigma value that falls within the sigma distribution
!     is computed from the index (additionally randomized at 0.01 level)
! JCL:(5/8/00)make integration step 100 times smaller, so divide by 1000
!        instead of 10
! JCL:turn off 0.01 level noise
      VALOUT=((KNDX-1)/1000.0)*SIGMA
!      VALOUT=((KNDX-1)/1000.0+(RVAL-0.5)/100.0)*SIGMA
!      VALOUT=((KNDX-1)/10.0+(RVAL-0.5)/100.0)*SIGMA
      IF(NEGV)VALOUT=-VALOUT

! JCL:
!      WRITE(45,*)SIGMA,VALOUT

      RETURN
      END
