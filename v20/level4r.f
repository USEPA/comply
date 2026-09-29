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
       REAL*8, ALLOCATABLE :: conc(:,:)              !radionuclide air concentration (Bq/m3)
       REAL*8, ALLOCATABLE :: dose(:,:)              !radionuclide dose (mrem/yr)
       REAL*8, ALLOCATABLE :: ihdose(:)
       REAL*8, ALLOCATABLE :: totihdose(:,:)

c Wind direction sector arrays (meteorological data)
       REAL*8, ALLOCATABLE :: windbinfreq(:)       !Sector (wind dir) frequency from that direction
       REAL*8, ALLOCATABLE :: windbinspeed(:)      !Sector (wind dir) average speed from that direction (m/s)
       REAL*8, ALLOCATABLE :: windbinsrdist(:)     !Sector (wind dir) source-receptor distance (m)

c Character string arrays for metadata
       CHARACTER*8,  ALLOCATABLE :: radname(:,:)     !radionuclide names at each stack
       CHARACTER*3,  ALLOCATABLE :: windbinname(:) !wind direction sector 
       CHARACTER*3,  ALLOCATABLE :: sectorname(:)  !sector name for receptors
       CHARACTER*3,  ALLOCATABLE :: srcsector(:)     !
       CHARACTER*3,  ALLOCATABLE :: verticalmomentum(:)    !vertical momentum (YES or NO)
       CHARACTER*3,  ALLOCATABLE :: meatsector(:)    !Sector (wind dir) location of meat
       CHARACTER*8,  ALLOCATABLE :: decaychain(:)    !Decay chain list
       CHARACTER*8,  ALLOCATABLE :: decaychainname(:)!Decay chain radionuclide names

c ========================================================================
c SECTION 2: DECLARE SCALAR VARIABLES FOR CALCULATIONS
c ========================================================================
c Conversion factors and working variables for Gaussian plume model

       INTEGER nwindbins          !number of wind direction bins at each stack
       REAL*8 pctcalm             !percentage of calm days at wind data location
       REAL*8 rainfall            !rainfall rate (m/yr)
       REAL*8 airtemp             !air temperature (F)


       REAL*8 FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF,FACTORG
       REAL*8 FACTORH,FACTORI,FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN
       REAL*8 area,xcrit,SIGMAZ,CAPSIG,pp,alldose,idose,dist,altconc,maxconc
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

       CHARACTER*200 windsource    !Source of wind data
       CHARACTER*200 winddates     !Date coverage of wind data
       CHARACTER*200 windloc       !Location of wind data
       CHARACTER*200 winddist      !Distance from wind data to source
       CHARACTER*200 title,company,facility,address,poc,phone,inline
       CHARACTER*8 NUC
       CHARACTER*2 NN

       integer date_time(8)                          !for date stamp function
       character*10 b(3)                             !for date stamp function

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
       write(*,*)'COMPLYr v2.0'
       write(*,*)
  
       open(11,file='input.level4r.txt',status='old',form='formatted')
       open(12,file='output.level4r.txt',status='unknown',form='formatted')
       open(13,file='diag.level4r.txt',status='unknown',form='formatted')

c ========================================================================
c SECTION 5: READ HEADER INFORMATION AND ALLOCATE STACK ARRAYS
c ========================================================================
c Read input namelist file

        read(11,'(54x,a)') title
        read(11,'(54x,a)') company
        read(11,'(54x,a)') facility
        read(11,'(54x,a)') address
        read(11,'(54x,a)') poc
        read(11,'(54x,a)') phone
        read(11,'(54x,i3)') nstacks

       read(11,'(a)') !blank spacer to improve readability of the input file

       read(11,'(54x,a)') windsource
       read(11,'(54x,a)') winddates
       read(11,'(54x,a)') windloc
       read(11,'(54x,a)') winddist
       read(11,'(54x,f10.5)') pctcalm
       read(11,'(54x,i3)') nwindbins

       allocate(windbinname(nwindbins))
       allocate(windbinfreq(nwindbins))
       allocate(windbinspeed(nwindbins))
       allocate(windbinsrdist(nwindbins))
       allocate(sectorname(nwindbins))

       do m = 1, nwindbins
       read(11,'(54x,a)') inline
       read(inline,*) windbinname(m),windbinfreq(m),windbinspeed(m),windbinsrdist(m),sectorname(m)
       enddo !end loop over reading in wind data and source-receptor distances

       read(11,'(54x,f10.5)') airtemp
       read(11,'(54x,f10.5)') rainfall

       airtemp = airtemp * (9.D0/5.D0) + 32.D0 !convert C to F
       rainfall = rainfall / 100.D0 !convert cm to m

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

        allocate(radname(nstacks, 1))
        allocate(rademis(nstacks))
        allocate(conc(nstacks, 1))
        allocate(dose(nstacks, 1))
        allocate(totihdose(nstacks, 1))

c ========================================================================
c SECTION 6: FIRST PASS - DETERMINE MAXIMUM ARRAY DIMENSIONS
c ========================================================================
        
        do n = 1, nstacks 
        read(11,'(a)') !blank spacer to improve readability of the input file 
        read(11,'(54x,f10.5)') latitude(n)
        read(11,'(54x,f10.5)') longitude(n)
        read(11,'(54x,f10.5)') releaseht(n)
        read(11,'(54x,f10.5)') flowrate(n)
        read(11,'(54x,f10.5)') diameter(n)
        read(11,'(54x,f10.5)') stacktemp(n)
        read(11,'(54x,a)') verticalmomentum(n)
c        read(11,'(54x,f10.5)') distance(n)
c        read(11,'(54x,a)') srcsector(n)

c ========================================================================
c SECTION 8: INPUT PARSING - READ RADIONUCLIDE EMISSION DATA
c ========================================================================
        
        read(11,'(54x,a)') inline
        read(inline,*) rademis(n)
        rademis(n) = rademis(n) * (1./365.)*(1./24.)*(1./60.)*(1./60.) !year to seconds

        enddo !end loop over stacks on input file

c ========================================================================
c SECTION 9: MAIN CALCULATION LOOP - AIR CONCENTRATIONS & INHALATION DOSE
c ========================================================================
c Calculate air concentrations using Gaussian plume model
c Emission rate assumption is Ci/s

        write(13,3001) !write header to file

        alldose = 0.0d0

        do n = 1, nstacks

        maxdose = 0.0d0
        inhalation(n) = 0.0d0

        do j = 1, nwindbins

c ---- OPTIMIZATION: Calculate stack-sector independent parameters once ----
c These depend only on n and j, not on receptor distance

        if(flowrate(n).le.0.0d0) then
        flowrate(n) = 0.3d0
        endif

        HEFF = releaseht(n)

        IF ( releaseht(n) .GT. 0.0d0 ) THEN

c ---- Temperature and Density Calculations ----
c Account for temperature-dependent air density corrections
         RHOAIR = 0.0735d0 * ( 80.0d0 + 460.0d0 ) / ( airtemp + 460.0d0 )
         RHOAIR = RHOAIR * 454.0d0 / ( 0.3048d0 ** 3 )
         RHOSTK = RHOAIR * ( airtemp + 460.0d0 ) / ( stacktemp(n) + 460.0d0 )

c ---- Momentum Rise Calculations ----
c Calculate buoyant and momentum rise effects on effective stack height
         AREA   = PI_CONST * ( diameter(n) ** 2.0d0 ) * 0.25d0
         STKVEL = flowrate(n) / AREA
         R      = STKVEL / windbinspeed(j)
         FM     = ( RHOSTK / RHOAIR ) * STKVEL * flowrate(n) * INV_PI
         TERM1  = 75.0d0 * FM / ( ( ( 1.0d0 + 3.0d0 / R ) * windbinspeed(j) ) ** 2 )
         TERM2  = releaseht(n) + 7.79d0 * SQRT( FM ) / ( ( 1.0d0 + 3.0d0 / R ) * 
     &           windbinspeed(j) )

c ---- Momentum Rise (if vertical stack) ----
c If the site has a vertical stack, then the momentum rise value applies
         IF ( verticalmomentum(n).eq.'YES' ) THEN
          DELHMO = 0.93d0 * ( TERM1 ** ( 3.0d0 / 7.0d0 ) ) *
     &        ( TERM2 ** ( 1.0d0 / 7.0d0 ) )
         ELSE
          DELHMO = 0.0d0
         ENDIF

c ---- Buoyant Rise Calculation ----
c Calculate rise due to temperature difference between stack and ambient air
         FB = 9.81d0 * ( RHOAIR - RHOSTK ) * flowrate(n) * INV_PI / RHOAIR
         DELHBU = 23.67d0 * ( ( FB / ( windbinspeed(j) ** 3 ) ) ** 0.6d0 ) *
     &           ( releaseht(n) + 195.1d0 * FB / ( windbinspeed(j) ** 3 ) ) ** 0.4d0

c ---- Effective Stack Height ----
c Use maximum of buoyant or momentum rise
        if(DELHMO.gt.DELHBU) then
         HEFF = releaseht(n) + DELHMO
        else
         HEFF = releaseht(n) + DELHBU
        endif

        ENDIF !end condition where H > 0

c ---- OPTIMIZATION: Precompute stack-independent factors ----
        reciprocal_windspeed = 1.0d0 / windbinspeed(j)
        one_minus_pctcalm = 1.0d0 - pctcalm

c ========================================================================
c SECTION 10: DISTANCE LOOP - GAUSSIAN PLUME CONCENTRATION & DOSE
c ========================================================================

        do receptor = 10, 51000, 10

        dist = dble(receptor)

c ---- Gaussian Plume Model Dispersion Parameter ----
c Calculate SIGMAZ and plume concentration factor (pp)
        SIGMAZ = (0.06d0 * dist) / SQRT( 1.0d0 + 0.0015d0 * dist)
        sigmaz_sqr = SIGMAZ * SIGMAZ
        CAPSIG = SQRT( sigmaz_sqr )  !area contribution set to 0 as in original
        zz = 0.5d0 * ( ( HEFF / CAPSIG ) ** 2 )

        if(zz.gt.0.001d0) then
          pp = 2.032d0 * EXP( -zz ) / ( dist * CAPSIG )
        else
          pp = 2.032d0 / (dist * CAPSIG)
        endif

c ---- Air Concentration Calculation ----
c Compute radionuclide air concentration with wind direction weighting
        conc(n,m) = ((windbinfreq(j) / one_minus_pctcalm) * 
     &             rademis(n) * pp)/ windbinspeed(j)

        if(conc(n,m).lt. 1.0d-32 ) then
          conc(n,m) = 1.0d-32
        endif

c ========================================================================
c SECTION 11: INHALATION DOSE CALCULATION (RADON + DTR MODEL)
c ========================================================================
c Calculate inhalation dose accounting for transit time and equilibrium

c ---- Transit Time Calculation ----
c Calculate time for air parcel to travel from emission to receptor
             T = TINIT + dist / ( windbinspeed(j) * 60.0d0)

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
           INHDOS = conc(n,m) * 1.0d12 * (1.0d0/1000.0d0) * 
     &             (1.0d0/100.0d0) * (8760.0d0/170.0d0) * INHFAC

c ---- Indoor/Outdoor Weighted Dose ----
c Combine indoor and outdoor doses weighted by time spent in each
       totihdose(n,m) = INHDOS * ( FINDOOR * EQILIN + 
     &                (1.0d0 - FINDOOR) * EQILOUT )

c ========================================================================
c SECTION 12: REPORT MAXIMUM DOSE AT RECEPTOR LOCATION
c ========================================================================
c Save dose information at key input receptor locations for report
        if(dist.eq.windbinsrdist(j)) then 
        write(*,*)'Air concentration:',conc(n,m),sectorname(j),dist,totihdose(n,m)
             tempdose = totihdose(n,m)
         if(tempdose.gt.maxdose) then 
             maxdose = tempdose
             maxsect = j
             maxdist = dist
             inhalation(n) = totihdose(n,m)
             distance(n) = maxdist
             srcsector(n) = sectorname(j)
             maxconc = conc(n,m)
         endif
        endif

        write(13,3000) n,sectorname(j),dist,1E+12*conc(n,m),totihdose(n,m)

        enddo !end loop over distances between source and receptor (receptor)
        enddo !end loop over number of wind direction bins (j)

        alldose = alldose + inhalation(n)

        write(*,*)'Stack ',n,' max dose (',maxdose,') and conc 
     &  (',1.0E12*maxconc,') 
     & receptor at ',windbinsrdist(maxsect),sectorname(maxsect)
        write(*,*)
        write(*,2999)'Inhalation dose is ',inhalation(n),' from stack',n
        write(*,*)'----------' 

        enddo !end loop over stacks (n)

        write(*,*)
        write(*,3500)'Total dose (mrem/yr):',alldose
        write(*,*)'Total dose limit is 10.0 mrem/yr'
        write(*,*)

 2999  FORMAT(a,f20.15,a,i3)
 3000  FORMAT(i2,',',a,',',f7.1,',',f30.20,1(',',f20.10))
 3001  FORMAT('STACK,SECTOR,DISTANCE,CONC,DOSE')
 3500  FORMAT(a,f10.4)

c ========================================================================
c SECTION 13: WRITE STANDARDIZED OUTPUT REPORT
c ========================================================================
c Write standardized output to output file

        write(12,4000)
        write(12,2001)
        write(12,1000) title
        write(12,1001) company
        write(12,1002) facility 
        write(12,1003) address
        write(12,1004) poc
        write(12,1005) phone

        do n = 1, nstacks
        write(12,2001)
        write(12,1006) n
        write(12,1007) latitude(n)
        write(12,1008) longitude(n)
        write(12,1009) releaseht(n)
        write(12,1010) diameter(n)
        write(12,1011) flowrate(n)
        write(12,1012) stacktemp(n)
        write(12,1013) airtemp
        write(12,1014) verticalmomentum(n)
        write(12,1015) distance(n)
        write(12,1019) srcsector(n)
        write(12,1016) rainfall
        write(12,1017) rademis(n)
        write(12,1018) inhalation(n)
        enddo !end loop over stacks

        write(12,2001)
        write(12,2000) alldose
        write(12,2001)
        write(12,2001)
        write(12,4000)

c ========================================================================
c SECTION 14: ADD TIMESTAMP AND FORMAT OUTPUT
c ========================================================================
c Add time stamp of model simulation to output file

        call date_and_time(b(1), b(2), b(3), date_time)

        print *,'Program execution date and time information:'
        print *,'year=',date_time(1)
        print *,'month_of_year=',date_time(2)
        print *,'day_of_month=',date_time(3)
        print *,'time difference in minutes=',date_time(4)
        print *,'hour of day=',date_time(5)
        print *,'minutes of hour=',date_time(6)
        print *,'seconds of minute=',date_time(7)
        print *,'milliseconds of second=',date_time(8)
        print *, 'DATE=',b(1)
        print *, 'TIME=',b(2)
        print *, 'ZONE=',b(3)

        write(12,*)
        write(12,2100) date_time(2),date_time(3),date_time(1),
     &           date_time(5),date_time(6)


c ========================================================================
c SECTION 15: FORMAT SPECIFICATIONS FOR OUTPUT REPORT
c ========================================================================
c Formatting for output file

 1000 FORMAT('Project:',a)
 1001 FORMAT('Company:',a)
 1002 FORMAT('Facility:',a)
 1003 FORMAT('Address:',a)
 1004 FORMAT('Point of Contact:',a)
 1005 FORMAT('Phone Number:',a)
 1006 FORMAT('Stack number:',i3)
 1007 FORMAT('Latitude (degrees):',f10.5)
 1008 FORMAT('Longitude (degrees):',f10.5)
 1009 FORMAT('Release height (m):',f10.1)
 1010 FORMAT('Stack diameter (m):',f10.1)
 1011 FORMAT('Flow rate (m3/s):',f10.1)
 1012 FORMAT('Stack temperature (F):',f10.1)
 1013 FORMAT('Air temperature (F):',f10.1)
 1014 FORMAT('Vertial momentum (YES/NO):',a)
 1015 FORMAT('Distance from source to receptor (m):',f10.1)
 1019 FORMAT('Receptor sector:',a3)
 1016 FORMAT('Rainfall rate (m/yr):',f10.1)
 1017 FORMAT('Emission rate (Ci/s):',f25.20)
 1018 FORMAT('Effective dose equivalent (mrem/yr):',f25.20)
 2000 FORMAT('Total effective dose equivalent all stacks (mrem/yr):',
     &       f10.5)
 2001 FORMAT('--------------------------------------------------------------')
 2100 FORMAT('Date/Time (LST) stamp of code execution (MM/DD/YYYY HH:MM):',
     &       i2,'/',i2,'/',i4,' ',i2,':',i2)

 4000 FORMAT('COMPLY-r version 2.0beta distributed through GitHub')

! End of program

       END 
