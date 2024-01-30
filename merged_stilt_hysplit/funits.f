!-----------------------------------------------------------
! Module to define FORTRAN I/O Unit Numbers 
!-----------------------------------------------------------
! Last Revised: 01 Apr 2004 (RRD) - initial version 
!               15 Jun 2005 (RRD) - split control & message
!               25 Jul 2008 (RRD) - gem routines
!-----------------------------------------------------------

MODULE funits
   
  INTEGER, PARAMETER :: KF01 = 10 ! Meteorological base 
  INTEGER, PARAMETER :: KF02 = 165 ! Meteorological maximum

  INTEGER, PARAMETER :: KF11 = 168 ! Conc or Traj base    
  INTEGER, PARAMETER :: KF12 = 169 ! Concentration maximum

  INTEGER, PARAMETER :: KF21 = 170 ! MESSAGE 
  INTEGER, PARAMETER :: KF22 = 171 ! STARTUP  
  INTEGER, PARAMETER :: KF23 = 172 ! Particle Dump Input
  INTEGER, PARAMETER :: KF24 = 173 ! Particle Dump Output
  INTEGER, PARAMETER :: KF25 = 174 ! CONTROL 
  INTEGER, PARAMETER :: KF26 = 175 ! namelist: SETUP.CFG (in), CONC.CFG (out); ZICONTROL (in), WINDERR (in), ZIERR (in)
  INTEGER, PARAMETER :: KFPARDAT = 176 ! PARTICLE.DAT
  INTEGER, PARAMETER :: KFJCLMSG = 180 ! JCLMESSAGE (also, before JCLMESSAGE is opened: ZSG_LEVS.IN)

  INTEGER, PARAMETER :: KF27 = 187 ! GEM concentration output
  INTEGER, PARAMETER :: KF28 = 188 ! GEM concentration dump
  INTEGER, PARAMETER :: KF29 = 189 ! GEM integrated profile 

  INTEGER, PARAMETER :: KF41 = 177 ! Landuse & ASCDATA.CFG
  INTEGER, PARAMETER :: KF42 = 178 ! Roughness length
  INTEGER, PARAMETER :: KF43 = 179 ! Terrain height

  INTEGER, PARAMETER :: KF31 = 180 ! Gridded emission array
  INTEGER, PARAMETER :: KF32 = 181 ! Temporal emission values
  INTEGER, PARAMETER :: KF33 = 182 ! IER | CHEM | PRCHEM

  INTEGER, PARAMETER :: KF50 = 190 ! lagrangian: LAGSET.CFG
  INTEGER, PARAMETER :: KF51 = 191 ! lagrangian: min out file
  INTEGER, PARAMETER :: KF52 = 199 ! lagrangian: max out file

END MODULE funits 
