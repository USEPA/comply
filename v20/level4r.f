       PROGRAM COMPLYR

c ========================================================================
c COMPLY-r v2.0
c
c PURPOSE: Radiological dose assessment model calculating air 
c          concentrations and inhalation doses from stack emissions 
c          of radon
c
c MAJOR SECTIONS:
c   - Declaration of variables and array allocations
c   - File I/O operations and input data reading
c   - Gaussian plume model calculations
c   - Inhalation dose computations
c   - Output report generation
c ========================================================================

       IMPLICIT NONE

c ========================================================================
c SECTION 1: DECLARE DYNAMIC ARRAYS FOR STACK AND ENVIRONMENTAL DATA
c ========================================================================
c Stack parameters allocated per stack

       REAL*8, ALLOCATABLE :: releaseht(:)           !stack height (m)
       REAL*8, ALLOCATABLE :: diameter(:)            !stack diameter (m)
       REAL*8, ALLOCATABLE :: flowrate(:)            !volumetric flow rate (m3/s)
       REAL*8, ALLOCATABLE :: distance(:)            !distance between stack and receptor (m)
       REAL*8, ALLOCATABLE :: stacktemp(:)           !stack temperature (F)
       REAL*8, ALLOCATABLE :: latitude(:)            !latitude (decimal degrees)
       REAL*8, ALLOCATABLE :: longitude(:)           !longitude (decimal degrees)
       REAL*8, ALLOCATABLE :: inhalation(:)          !inhalation dose (mrem/yr)

c Radionuclide and dose arrays (per stack, per radionuclide)
       REAL*8, ALLOCATABLE :: rademis(:)           !radionuclide emissions (Ci/s) for each stack
       REAL*8, ALLOCATABLE :: maxconc(:)              !radionuclide air concentration (Bq/m3)
       REAL*8, ALLOCATABLE :: dose(:,:)              !radionuclide dose (mrem/yr)
       REAL*8, ALLOCATABLE :: ihdose(:)

c Wind direction sector arrays (meteorological data)
       REAL*8, ALLOCATABLE :: windbinfreq(:,:)       !Sector (wind dir) frequency from that direction
       REAL*8, ALLOCATABLE :: windbinspeed(:,:)      !Sector (wind dir) average speed from that direction (m/s)
       REAL*8, ALLOCATABLE :: windbinsrdist(:,:)     !Sector (wind dir) source-receptor distance (m)
       REAL*8, ALLOCATABLE ::  pctcalm(:)            !percentage of calm days at wind data location
       REAL*8, ALLOCATABLE ::  airtemp(:)            !air temperature (F)
       INTEGER, ALLOCATABLE :: nwindbins(:)          !number of wind direction bins at each stack

c Character string arrays for metadata
       CHARACTER*8,  ALLOCATABLE :: radname(:,:)       !radionuclide names at each stack
       CHARACTER*3,  ALLOCATABLE :: windbinname(:,:)   !wind direction sector 
       CHARACTER*3,  ALLOCATABLE :: sectorname(:,:)    !sector name for receptors
       CHARACTER*3,  ALLOCATABLE :: srcsector(:)       !
       CHARACTER*3,  ALLOCATABLE :: verticalmomentum(:)!vertical momentum (YES or NO)
       CHARACTER*3,  ALLOCATABLE :: meatsector(:)      !Sector (wind dir) location of meat
       CHARACTER*8,  ALLOCATABLE :: decaychain(:)      !Decay chain list
       CHARACTER*8,  ALLOCATABLE :: decaychainname(:)  !Decay chain radionuclide names

       CHARACTER*200,  ALLOCATABLE :: windsource(:)    !Source of wind data
       CHARACTER*200,  ALLOCATABLE :: winddates(:)     !Date coverage of wind data
       CHARACTER*200,  ALLOCATABLE :: windloc(:)       !Location of wind data
       CHARACTER*200,  ALLOCATABLE :: winddist(:)      !Distance from wind data to source


c ========================================================================
c SECTION 2: DECLARE SCALAR VARIABLES FOR CALCULATIONS
c ========================================================================
c Conversion factors and working variables for Gaussian plume model

       REAL*8 totihdose

       REAL*8 FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF,FACTORG
       REAL*8 FACTORH,FACTORI,FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN
       REAL*8 area,xcrit,SIGMAZ,CAPSIG,pp,alldose,idose,dist,altconc
       REAL*8 tempdose,maxdose,maxdist,totmilkdose,totvegdose,itotvegdose
       REAL*8 totmeatdose,itotcowdose
       REAL*8 capp,capn0,outprod,inprod,z1,term1,term2,sum
       REAL*8 transittime                            !transit time of air from emission release to receptor (s)
       REAL*8 fdux                                   !not sure what this is-units are 1/(m2/s) or s/m2
       REAL*8 wetdepadj                              !wet deposition adjustment factor
       REAL*8 xtime
       REAL*8 xremove
       REAL*8 RHOAIR,RHOSTK,STKVEL,R,FM,DELHMO,DELHBU,HEFF,FB,zz,INHFAC
       REAL*8 TINIT,T,EQILOUT,EQILIN,INHDOS,FINDOOR

c Mathematical and precomputed constants for efficiency
       REAL*8 PI_CONST                               !mathematical constant pi
       REAL*8 INV_PI                                 !1/pi for efficiency
       REAL*8 reciprocal_windspeed                   !precomputed 1/windbinspeed
       REAL*8 sigmaz_sqr                             !sigmaz squared
       REAL*8 one_minus_pctcalm                      !precomputed 1-pctcalm(n)

c Loop counters and program control variables
       INTEGER nstacks                               !number of unique stacks
       INTEGER receptor, N, M, I, J, K, L, II, JJ, P, Q !loop counters
       INTEGER nchain                                !number of nuclides in the decay series
       INTEGER maxsect
       INTEGER RECEP(16)


       CHARACTER*200 title,company,facility,address,poc,phone,inline
       CHARACTER*8 NUC
       CHARACTER*3 DIAGFLAG

       integer date_time(8)                          !for date stamp function
       character*10 b(3)                             !for date stamp function

c Wind direction and downwind receptor relationships
!       DATA FROM   / 'N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW' /
!       DATA TOWARD / 'S','SSW','SW','WSW','W','WNW','NW','NNW','N','NNE','NE','ENE','E','ESE','SE','SSE' /
!       DATA WINDS  /   1,    2,   3,    4,  5,    6,   7,    8,  9,   10,  11,   12, 13,   14,  15,  16 /
       DATA RECEP  /    9,   10,  11,   12, 13,   14,  15,   16,  1,    2,   3,    4,  5,    6,   7,   8 /

c ========================================================================
c SECTION 3: INITIALIZE MATHEMATICAL CONSTANTS
c ========================================================================
! Constants
       PI_CONST = ACOS(-1.0d0)                       !most accurate pi
       INV_PI = 1.0d0 / PI_CONST
       FINDOOR = 0.75d0
       TINIT = 11.8d0
       INHFAC = 0.0092d0 * 100.0d0 * 1000.0d0

c ========================================================================
c SECTION 4: OPEN INPUT AND OUTPUT FILES
c ========================================================================
! Open input and output files 

       write(*,*)
       write(*,*)'COMPLY-r v2.0'
       write(*,*)
  
       open(11,file='input.level4r.txt',status='old',form='formatted')
       open(12,file='output.level4r.txt',status='unknown',form='formatted')
       open(13,file='diag.level4r.txt',status='unknown',form='formatted')

c ========================================================================
c SECTION 5: READ HEADER INFORMATION AND ALLOCATE STACK ARRAYS
c ========================================================================
c Read input namelist file

        read(11,'(38x,a)') title
        read(11,'(38x,a)') company
        read(11,'(38x,a)') facility
        read(11,'(38x,a)') address
        read(11,'(38x,a)') poc
        read(11,'(38x,a)') phone
        read(11,'(38x,a3)') diagflag
        read(11,'(38x,a)') inline
        read(inline,*) nstacks

c ============================================================================
c SECTION 16: MEMORY ALLOCATION - ALLOCATE STACK-LEVEL PARAMETER ARRAYS
c ============================================================================
c Allocate arrays to store stack-specific parameters (heights, distances,
c temperatures, etc.) based on the number of stacks

        allocate(latitude(nstacks))
        allocate(longitude(nstacks))
        allocate(releaseht(nstacks))
        allocate(flowrate(nstacks))
        allocate(diameter(nstacks))
        allocate(distance(nstacks))
        allocate(stacktemp(nstacks))
        allocate(verticalmomentum(nstacks))
        allocate(srcsector(nstacks))
        allocate(inhalation(nstacks))
        allocate(pctcalm(nstacks))
        allocate(airtemp(nstacks))
        allocate(nwindbins(nstacks))
        allocate(windsource(nstacks))
        allocate(winddates(nstacks))
        allocate(windloc(nstacks))
        allocate(winddist(nstacks))

        allocate(radname(nstacks, 1))
        allocate(rademis(nstacks))
        allocate(maxconc(nstacks))
        allocate(dose(nstacks, 1))

        allocate(windbinname(nstacks,16))
        allocate(windbinfreq(nstacks,16))
        allocate(windbinspeed(nstacks,16))
        allocate(windbinsrdist(nstacks,16))
        allocate(sectorname(nstacks,16))


       do n = 1, nstacks


       read(11,'(a)') !blank spacer to improve readability of the input file

       read(11,'(38x,a)') windsource(n)
       read(11,'(38x,a)') winddates(n)
       read(11,'(38x,a)') windloc(n)
       read(11,'(38x,a)') winddist(n)
       read(11,'(38x,a)') inline
       read(inline,*) airtemp(n)
       read(11,'(38x,a)') inline
       read(inline,*) nwindbins(n)

       sum = 0.0
       do m = 1, nwindbins(n)
       read(11,'(38x,a)') inline
       read(inline,*) windbinname(n,m),windbinfreq(n,m),windbinspeed(n,m)
       sum = sum + windbinfreq(n,m)
       enddo !end loop over reading in wind data and source-receptor distances

       airtemp(n) = airtemp(n) * (9.D0/5.D0) + 32.D0 !convert C to F

       !check that calms make sense with the wind rose information
       if(sum.lt.1.0.and.nwindbins(n).gt.1) then
        pctcalm(n) = 1 - sum
       elseif(sum.gt.1.0) then
        write(*,*)'ERROR! Sum of wind frequencies exceed 1.0'
        stop
       else
        pctcalm(n) = 0.0
       endif
        write(*,*)'Percent calms = ',pctcalm(n)


       read(11,'(a)') !blank spacer to improve readability of the input file

       if(nwindbins(n).eq.1) then !simple conservative estimate of winds no wind rose
       read(11,'(38x,a)') inline
       read(inline,*) sectorname(n,1),windbinsrdist(n,1)
       else !full 16 direction wind rose
       do m = 1, nwindbins(n)
       read(11,'(38x,a)') inline
       read(inline,*) sectorname(n,RECEP(m)),windbinsrdist(n,RECEP(m))
       enddo !end loop over reading in wind data and source-receptor distances
       endif


        read(11,'(a)') !blank spacer to improve readability of the input file 
        read(11,'(38x,a)') inline
        read(inline,*) latitude(n)
        read(11,'(38x,a)') inline
        read(inline,*) longitude(n)
        read(11,'(38x,a)') inline
        read(inline,*) releaseht(n)
        read(11,'(38x,a)') inline
        read(inline,*) flowrate(n)
        read(11,'(38x,a)') inline
        read(inline,*) diameter(n)
        read(11,'(38x,a)') inline
        read(inline,*) stacktemp(n)
        read(11,'(38x,a)') verticalmomentum(n)

c ========================================================================
c SECTION 8: INPUT PARSING - READ RADIONUCLIDE EMISSION DATA
c ========================================================================
        
        read(11,'(38x,a)') inline
        read(inline,*) rademis(n)
        rademis(n) = rademis(n) * (1./365.)*(1./24.)*(1./60.)*(1./60.) !year to seconds

        enddo !end loop over stacks on input file

c ========================================================================
c SECTION 9: MAIN CALCULATION LOOP - AIR CONCENTRATIONS & INHALATION DOSE
c ========================================================================
c Calculate air concentrations using Gaussian plume model
c Emission rate assumption is Ci/s


        alldose = 0.0d0

        do n = 1, nstacks

        maxdose = 0.0d0
        inhalation(n) = 0.0d0

        do j = 1, nwindbins(n)

        dist = windbinsrdist(n,j)


       !Air concentration
       CALL calc_air_concentration( dist, pctcalm(n), releaseht(n), diameter(n),
     & flowrate(n), stacktemp(n), verticalmomentum(n),
     & windbinfreq(n,j), windbinspeed(n,j), rademis(n), altconc, airtemp(n), transittime)


       CALL calc_workinglevel(dist, INHFAC, TINIT,
     & windbinfreq(n,j), windbinspeed(n,j), rademis(n), altconc, totihdose)


c ========================================================================
c SECTION 12: REPORT MAXIMUM DOSE AT RECEPTOR LOCATION
c ========================================================================
c Save dose information at key input receptor locations for report

       write(*,*)'Air concentration:',altconc,sectorname(n,j),dist,totihdose

         if(totihdose.gt.maxdose) then 
             maxdose = totihdose
             maxsect = j
             maxdist = dist
             inhalation(n) = totihdose
             distance(n) = maxdist
             srcsector(n) = sectorname(n,j)
             maxconc(n) = altconc
         endif


        enddo !end loop over number of wind direction bins (j)

        alldose = alldose + inhalation(n)

        write(*,*)'Stack ',n,' max dose (',maxdose,') and conc 
     &  (',1.0E12*maxconc(n),') 
     & receptor at ',windbinsrdist(n,maxsect),sectorname(n,maxsect)
        write(*,*)
        write(*,2999)'Inhalation dose is ',inhalation(n),' from stack',n
        write(*,*)'----------' 

        enddo !end loop over stacks (n)

        write(*,*)
        write(*,3500)'Total dose (mrem/yr):',alldose
        write(*,*)'Total dose limit is 10.0 mrem/yr'
        write(*,*)

 2999  FORMAT(a,ES9.2,a,i3)
 3500  FORMAT(a,ES9.2)

c ========================================================================
c SECTION 13: WRITE STANDARDIZED OUTPUT REPORT
c ========================================================================


 4000 FORMAT('COMPLY-r version 2.0 level 4 Dose Summary Report')
 2001 FORMAT('------------------------------------------------------------------------------')
 2002 FORMAT(' ')

 1000 FORMAT('          Project:',a)
 1001 FORMAT('          Company:',a)
 1002 FORMAT('         Facility:',a)
 1003 FORMAT('          Address:',a)
 1004 FORMAT(' Point of Contact:',a)
 1005 FORMAT('     Phone Number:',a)

       write(12,2001)
       write(12,2001)
       write(12,4000)
       write(12,2001)
       write(12,2001)
       write(12,2002)
       write(12,1000) trim(title)
       write(12,1001) trim(company)
       write(12,1002) trim(facility)
       write(12,1003) trim(address)
       write(12,1004) trim(poc)
       write(12,1005) trim(phone)
       write(12,2002)
       write(12,2001)


 1102 FORMAT('DOSE ESTIMATES')
 1100 FORMAT('Total Effective Dose Equivalent all stacks (mrem/year):',ES9.2)
       write(12,1102)
       write(12,2002)
       write(12,1100) alldose

       write(12,2002)

 1101 FORMAT('Maximum dose receptor was located at ',f8.0,' m downwind in sector ',a)

       write(12,2002)
       write(12,2002)


       do n = 1, nstacks
       write(12,1006) n
       write(12,1101) distance(n),srcsector(n)
       write(12,2002)
       write(12,1998)

 1998 FORMAT('Emissions(Ci/y),Conc(Ci/m3),Inhalation Dose')
 1999 FORMAT(ES11.3,',',ES11.3,',',ES11.3)
!       !Dose units are mrem/yr and concentration units pCi/m3 (convert Ci/m3 to pCi/m3 here)
       write(12,1999) rademis(n)*60.*60.*24.*365.,maxconc(n),inhalation(n)

       write(12,2002)
       enddo




 1030 FORMAT('STACK INFORMATION')
 1006 FORMAT('Stack number                   :',i3)
 1007 FORMAT('Stack latitude (degrees)       :',f10.5)
 1008 FORMAT('Stack longitude (degrees)      :',f10.5)
 1009 FORMAT('Release height (m)             :',f10.1)
 1010 FORMAT('Stack diameter (m)             :',f10.1)
 1011 FORMAT('Flow rate (m3/s)               :',f10.1)
 1012 FORMAT('Stack temperature (C)          :',f10.1)
 1014 FORMAT('Vertical momentum (YES/NO)     :',a)
 1015 FORMAT('Building height (m)            :',f10.1)
 1016 FORMAT('Building width (m)             :',f10.1)
 1019 FORMAT('Number of radionuclides        :',I2)
 1033 FORMAT('Veg distance(m),downwind sector:',f10.1,1x,a3)
 1034 FORMAT('Milk dist(m),downwind sector   :',f10.1,1x,a3)
 1035 FORMAT('Meat dist(m),downwind sector   :',f10.1,1x,a3)

       write(12,2001)
       write(12,1030)
       write(12,2002)
       do n = 1,nstacks
       write(12,1006) n
       write(12,1007) latitude(n)
       write(12,1008) longitude(n)
       write(12,1009) releaseht(n)
       write(12,1010) diameter(n)
       write(12,1011) flowrate(n)
       write(12,1012) stacktemp(n)
       write(12,1014) trim(verticalmomentum(n))
       write(12,2002)
       enddo




 1025 FORMAT('METEOROLOGY')
 1013 FORMAT('Air temperature (C)          :',f10.1)
 1022 FORMAT('Calm winds (fraction)        :',f10.3)
 1023 FORMAT('Sector winds from, frequency(fraction), speed(m/s)')
 1024 FORMAT(a3,14x,',',3x,f10.3,4x,',',f10.3)
 1026 FORMAT('Sector winds to  , downwind receptor(m)')
 1027 FORMAT(a3,14x,',',f10.3)

       write(12,2002)
       write(12,2001)
       write(12,1025)
       write(12,2002)
       
     
        do n = 1, nstacks
        write(12,2002) 
        write(12,1006) n
        write(12,1013) airtemp(n)
         write(12,1022) pctcalm(n)
         write(12,1023)
         do m = 1, nwindbins(n)
        write(12,1024)windbinname(n,m),windbinfreq(n,m),windbinspeed(n,m)
        enddo

        write(12,2002)
        write(12,1026)

        do m = 1, nwindbins(n)
        write(12,1027)sectorname(n,m),windbinsrdist(n,m)
        enddo
        enddo



c ========================================================================
c EXECUTION TIMESTAMP AND FORMAT OUTPUT
c ========================================================================
c Add time stamp of model simulation to output file

        call date_and_time(b(1), b(2), b(3), date_time)


 2100 FORMAT('Date/Time (LST) stamp of code execution (MM/DD/YYYY HH:MM):',i2.2,'/',i2.2,'/',i4,' ',i2.2,':',i2.2)

       write(12,2001)
       write(12,2002)
       write(12,2100) date_time(2),date_time(3),date_time(1),
     &           date_time(5),date_time(6)

       write(*,2002)
       write(*,2100) date_time(2),date_time(3),date_time(1),
     &           date_time(5),date_time(6)



c ===========================================================================
c Diagnostic output file
c ===========================================================================

       if(diagflag.eq.'YES'.or.diagflag.eq.'TRUE'.or.
     &    diagflag.eq.'yes'.or.diagflag.eq.'true') then

        write(13,3002) 
        write(13,3001) !write header to file

        alldose = 0.0d0

        do n = 1, nstacks

        maxdose = 0.0d0
        inhalation(n) = 0.0d0

        do j = 1, nwindbins(n)

        do receptor = 25, 50000, 25

        dist = dble(receptor)

       !Air concentration
       CALL calc_air_concentration( dist, pctcalm(n), releaseht(n), diameter(n),
     & flowrate(n), stacktemp(n), verticalmomentum(n),
     & windbinfreq(n,j), windbinspeed(n,j), rademis(n), altconc, airtemp(n), transittime)


       CALL calc_workinglevel(dist, INHFAC, TINIT,
     & windbinfreq(n,j), windbinspeed(n,j), rademis(n), altconc, totihdose)


        write(13,3000) n,sectorname(n,j),dist,altconc,totihdose

        enddo
        enddo
        enddo
 
        endif

 3000  FORMAT(i2,',',a,',',f7.1,',',ES10.3,1(',',ES10.3))
 3001  FORMAT('STACK,SECTOR,DISTANCE,CONC,DOSE')
 3002  FORMAT('# Distance as m; concentration as Ci/m3; dose as mrem/year')

c ===========================================================================

! End of program

       END PROGRAM

c ============================================================================
c START OF SUBROUTINES
c ============================================================================

       SUBROUTINE calc_air_concentration(dist, xpctcalm,
     & xreleaseht, xdiameter, xflowrate, xstacktemp, xverticalmomentum,
     & xwindbinfreq, xwindbinspeed, xrademis, xconc, xairtemp, transittime)

       IMPLICIT NONE

       REAL*8  :: dist, xpctcalm, xcrit, xairtemp
       REAL*8  :: xreleaseht, xdiameter, xflowrate, xstacktemp
       REAL*8  :: xbldheight, xbldwidth, xwindbinfreq, xwindbinspeed
       REAL*8  :: xrademis
       CHARACTER*3 :: xverticalmomentum

       REAL*8 :: xconc
       REAL*8 :: HEFF, SIGMAZ, CAPSIG, transittime, pp

       ! Local variables
       REAL*8 :: RHOAIR, RHOSTK, STKVEL, R, FM, DELHMO, DELHBU
       REAL*8 :: FB, AREA, zz

       REAL*8, PARAMETER :: PI = 3.14159265358979D0


       ! ------------------------------------------------------------
       ! Effective stack height
       ! ------------------------------------------------------------

        if(xflowrate.le.0.0d0) then
        xflowrate = 0.3d0
        endif

       HEFF = xreleaseht

       IF (xreleaseht .GT. 0.0) THEN

       ! Temperature and density
       RHOAIR = 0.0735D0 * (80.0D0 + 460.0D0) / (xairtemp + 460.0D0)
       RHOAIR = RHOAIR * 454.0D0 / (0.3048D0 ** 3)
       RHOSTK = RHOAIR * (xairtemp + 460.0D0) / (xstacktemp + 460.0D0)

       ! Momentum rise
       AREA   = PI * (xdiameter ** 2.0D0) / 4.0D0
       STKVEL = xflowrate / AREA
       R      = STKVEL / xwindbinspeed
       FM     = (RHOSTK / RHOAIR) * STKVEL * xflowrate / PI

       IF (xverticalmomentum .EQ. 'YES') THEN
        DELHMO = 0.93D0 * (75.0D0 * FM / (((1.0D0 + 3.0D0 / R) * xwindbinspeed) ** 2.0D0)) ** (3.0D0 / 7.0D0) *
     &            (xreleaseht + 7.79D0 * SQRT(FM) / ((1.0D0 + 3.0D0 / R) * xwindbinspeed)) ** (1.0D0 / 7.0D0)
       ELSE
        DELHMO = 0.0D0
       END IF

       ! Buoyancy rise
       FB = 9.81D0 * (RHOAIR - RHOSTK) * xflowrate / (PI * RHOAIR)
       DELHBU = 23.67D0 * ((FB / (xwindbinspeed ** 3.0D0)) ** 0.6D0) *
     &        (xreleaseht + 195.1D0 * FB / (xwindbinspeed ** 3.0D0)) ** 0.4D0

       ! Use larger of buoyancy or momentum rise
       IF (DELHMO .GT. DELHBU) THEN
         HEFF = xreleaseht + DELHMO
       ELSE
         HEFF = xreleaseht + DELHBU
       END IF

       END IF

       ! ------------------------------------------------------------
       ! Gaussian plume concentration term
       ! ------------------------------------------------------------
       SIGMAZ = (0.06D0 * dist) / SQRT(1.0D0 + 0.0015D0 * dist)
       CAPSIG = SQRT(SIGMAZ * SIGMAZ)

       zz = 0.5D0 * ((HEFF / CAPSIG) ** 2.0D0)

       IF (zz .GT. 0.001D0 ) THEN
        pp = 2.032D0 * EXP(-zz) / (dist * CAPSIG)
       ELSE
        pp = 2.032D0 / (dist * CAPSIG)
       END IF

       ! Transit time from source to receptor (days)
       transittime = (dist / xwindbinspeed) * (1.0D0 / 3600.0D0) * (1.0D0 / 24.0D0)

       ! ------------------------------------------------------------
       ! Near-field / far-field concentration
       ! ------------------------------------------------------------

       xconc = ((xwindbinfreq / (1.0D0 - xpctcalm)) * xrademis * pp) / xwindbinspeed

       IF (xconc .LT. 1.0D-32) xconc = 1.0D-32

       END SUBROUTINE calc_air_concentration

!=================================================================================================

       SUBROUTINE calc_workinglevel(dist, INHFAC, TINIT,
     & xwindbinfreq, xwindbinspeed, xrademis, xconc, totihdose)

       IMPLICIT NONE

       REAL*8  :: dist
       REAL*8  :: xwindbinfreq, xwindbinspeed
       REAL*8  :: xrademis
       REAL*8 :: xconc, totihdose

       ! Local variables
       REAL*8 :: INHFAC
       REAL*8 :: TINIT,T,EQILOUT,EQILIN,INHDOS,FINDOOR

       FINDOOR = 0.75d0

c Calculate inhalation dose accounting for transit time and equilibrium

c ---- Transit Time Calculation ----
c Calculate time for air parcel to travel from emission to receptor
             T = TINIT + dist / ( xwindbinspeed * 60.0d0)

c ---- Outdoor Equilibrium Fraction ----
c Account for decay of short-lived progeny during transit
           EQILOUT = 1.0d0 - 0.0479d0 * EXP(-T / 4.39d0) -
     &              2.1963d0 * EXP(-T / 38.6d0) +
     &              1.2442d0 * EXP(-T / 28.4d0)

             IF( EQILOUT .GT. 1.0d0 ) EQILOUT = 1.0d0

c ---- Indoor Equilibrium Fraction ----
c Account for reduced progeny contribution indoors
           EQILIN =  0.35d0 * ( 1.0d0 + EQILOUT )
             IF( EQILIN .GT. 0.7d0 ) EQILIN = 0.7d0

c ---- Total Inhalation Dose ----
c Convert concentration to dose with unit conversions
c Ci/m**3   pCi/Ci    m**3/liter   WL-liter/pCi  WLmonths/year  mrem/WLM
           INHDOS = xconc * 1.0d12 * (1.0d0/1000.0d0) *
     &             (1.0d0/100.0d0) * (8760.0d0/170.0d0) * INHFAC


c ---- Indoor/Outdoor Weighted Dose ----
c Combine indoor and outdoor doses weighted by time spent in each
       totihdose = INHDOS * ( FINDOOR * EQILIN +
     &                (1.0d0 - FINDOOR) * EQILOUT )

       END SUBROUTINE calc_workinglevel
