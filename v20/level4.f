      PROGRAM COMPLY

      IMPLICIT NONE

c ============================================================================
c COMPLY v2.0 - LEVEL 4 RADIONUCLIDE DOSE ASSESSMENT MODEL
c ============================================================================
c This program models atmospheric dispersion and radiological dose 
c calculations for radionuclide emissions from stack sources. It computes 
c air concentrations, environmental transfer, and resulting human doses 
c via multiple pathways (inhalation, immersion, ingestion, surface deposition).
c
c Key features:
c   - Gaussian plume dispersion model with buoyancy/momentum rise corrections
c   - Radioactive decay chain calculations for parent/progeny nuclides
c   - Multiple dose pathways: inhalation, immersion, surface, ingestion
c   - Ingestion via vegetables, milk, and meat
c   - Sector-based wind direction analysis
c ============================================================================

c ============================================================================
c SECTION 1: DYNAMIC ARRAY DECLARATIONS - METEOROLOGICAL & STACK PARAMETERS
c ============================================================================
c These arrays store input parameters for each stack source and associated
c meteorological/geographic data

      REAL*8, ALLOCATABLE :: releaseht(:)           !stack height (m)
      REAL*8, ALLOCATABLE :: diameter(:)            !stack diameter (m)
      REAL*8, ALLOCATABLE :: flowrate(:)            !volumetric flow rate (m3/s)
      REAL*8, ALLOCATABLE :: stacktemp(:)           !stack temperature (C)
      REAL*8, ALLOCATABLE :: bldheight(:)           !building height (m)
      REAL*8, ALLOCATABLE :: bldwidth(:)            !building width (m)
      REAL*8, ALLOCATABLE :: distancemei(:)         !distance between source and max exposed individual (m)
      REAL*8, ALLOCATABLE :: distanceveg(:)         !distance between source and garden receptor (m)
      REAL*8, ALLOCATABLE :: distancemilk(:)        !distance between source and cows/milk receptor (m)
      REAL*8, ALLOCATABLE :: distancemeat(:)        !distance between source and cows/meat receptor (m)
      REAL*8, ALLOCATABLE :: latitude(:)            !latitude (decimal degrees)
      REAL*8, ALLOCATABLE :: longitude(:)           !longitude (decimal degrees)

c ============================================================================
c SECTION 2: DYNAMIC ARRAY DECLARATIONS - EMISSION & CONCENTRATION DATA
c ============================================================================
c These arrays store radionuclide-specific emission rates and computed
c air/dose concentrations for each stack/radionuclide combination

      REAL*8, ALLOCATABLE :: rademis(:,:)           !radionuclide emissions (Ci/yr) for each stack
      REAL*8, ALLOCATABLE :: maxconc(:,:)              !radionuclide air concentration (Bq/m3)

c ============================================================================
c SECTION 3: DYNAMIC ARRAY DECLARATIONS - RADIOACTIVE DECAY CHAIN PARAMETERS
c ============================================================================
c These arrays track source terms, initial amounts, and activities for
c decay chains (parent nuclide + radioactive progeny)

      REAL*8, ALLOCATABLE :: sourceterm(:)          !source term for decay calculation
      REAL*8, ALLOCATABLE :: initialamt(:)          !initial amount of nuclide
      REAL*8, ALLOCATABLE :: activity(:)            !computed activity (Ci/m3)
      REAL*8, ALLOCATABLE :: surfaceactivity(:)     !activity on ground surface
      REAL*8, ALLOCATABLE :: ambientactivity(:)     !activity in air
      REAL*8, ALLOCATABLE :: normsourceactivity(:)  !normalized surface activity

c ============================================================================
c SECTION 4: DYNAMIC ARRAY DECLARATIONS - INGESTION PATHWAY TRANSFER FACTORS
c ============================================================================
c These arrays store intermediate calculated activities used in computing
c radionuclide concentrations in vegetables, milk, and meat

      REAL*8, ALLOCATABLE :: aveglf(:)              !vegetable leaf activity (short-term weathering)
      REAL*8, ALLOCATABLE :: avegrt(:)              !vegetable root activity (long-term deposition)
      REAL*8, ALLOCATABLE :: amilklf(:)             !milk forage leaf activity (short-term)
      REAL*8, ALLOCATABLE :: amilkrt(:)             !milk forage root activity (long-term)
      REAL*8, ALLOCATABLE :: ameatlf(:)             !meat forage leaf activity (short-term)
      REAL*8, ALLOCATABLE :: ameatrt(:)             !meat forage root activity (long-term)
      REAL*8, ALLOCATABLE :: cfveg(:)               !concentration factor in vegetables (Ci/kg)
      REAL*8, ALLOCATABLE :: cfmilk(:)              !concentration factor in milk (Ci/L)
      REAL*8, ALLOCATABLE :: cfmeat(:)              !concentration factor in meat (Ci/kg)
      REAL*8, ALLOCATABLE :: dveg(:)                !dose from vegetable ingestion (mrem/yr)
      REAL*8, ALLOCATABLE :: dmilk(:)               !dose from milk ingestion (mrem/yr)
      REAL*8, ALLOCATABLE :: dmeat(:)               !dose from meat ingestion (mrem/yr)
      REAL*8, ALLOCATABLE :: ingdose(:)             !total ingestion dose

c ============================================================================
c SECTION 5: DYNAMIC ARRAY DECLARATIONS - DEPOSITION & DECAY FACTORS
c ============================================================================
c These arrays store parameters controlling deposition rates and radioactive
c decay of nuclides and their environmental transfer

      REAL*8, ALLOCATABLE :: kvalue(:)              !decay constant (1/day)
      REAL*8, ALLOCATABLE :: ktot(:)                !total removal constant (decay + deposition)
      REAL*8, ALLOCATABLE :: xbveg(:,:,:)           !root uptake factor for vegetables
      REAL*8, ALLOCATABLE :: xbfor(:,:,:)           !root uptake factor for forage
      REAL*8, ALLOCATABLE :: xfmilk(:,:,:)          !transfer factor to milk
      REAL*8, ALLOCATABLE :: xfmeat(:,:,:)          !transfer factor to meat
      REAL*8, ALLOCATABLE :: ybveg(:)           !root uptake factor for vegetables
      REAL*8, ALLOCATABLE :: ybfor(:)           !root uptake factor for forage
      REAL*8, ALLOCATABLE :: yfmilk(:)          !transfer factor to milk
      REAL*8, ALLOCATABLE :: yfmeat(:)          !transfer factor to meat

c ============================================================================
c SECTION 6: DYNAMIC ARRAY DECLARATIONS - RADIONUCLIDE DOSE FACTORS
c ============================================================================
c These arrays store dose conversion factors (DCF) from NCRP/ICRP references
c that convert nuclide activity to effective dose equivalent

      REAL*8, ALLOCATABLE :: radchainlen(:,:)       !decay chain length 
      REAL*8, ALLOCATABLE :: halflife(:,:,:)        !radionuclide half life (years)
      REAL*8, ALLOCATABLE :: xhalflife(:)           !radionuclide half life (years)
      REAL*8, ALLOCATABLE :: inhalationfac(:,:,:)   !inhalation dose factor
      REAL*8, ALLOCATABLE :: immersionfac(:,:,:)    !immersion dose factor
      REAL*8, ALLOCATABLE :: surfacefac(:,:,:)      !surface dose factor
      REAL*8, ALLOCATABLE :: ingestionfac(:,:,:)    !ingestion dose factor
      REAL*8, ALLOCATABLE :: drydepvel(:,:,:)       !dry deposition velocity (m/day)
      REAL*8, ALLOCATABLE :: dosefactor(:)          !temp array for dose factors
      REAL*8, ALLOCATABLE :: xdrydepvel(:)          !dry deposition velocity (m/day)

c ============================================================================
c SECTION 7: DYNAMIC ARRAY DECLARATIONS - CUMULATIVE DOSE ARRAYS
c ============================================================================
c These arrays accumulate total doses across decay chains and receptor
c locations for immersion, inhalation, surface deposition, and ingestion

      REAL*8, ALLOCATABLE :: imdose(:)              !immersion dose component
      REAL*8, ALLOCATABLE :: ihdose(:)              !inhalation dose component
      REAL*8, ALLOCATABLE :: surfdose(:)            !surface dose component
      REAL*8, ALLOCATABLE :: totimdose(:,:)         !total immersion dose per stack/radionuclide
      REAL*8, ALLOCATABLE :: totihdose(:,:)         !total inhalation dose per stack/radionuclide
      REAL*8, ALLOCATABLE :: totsurfdose(:,:)       !total surface dose per stack/radionuclide
      REAL*8, ALLOCATABLE :: totingdose(:,:)        !total ingestion dose per stack/radionuclide
      REAL*8, ALLOCATABLE :: sectordose(:)          !total dose by sector bin

c ============================================================================
c SECTION 8: DYNAMIC ARRAY DECLARATIONS - PATHWAY-SPECIFIC DOSES
c ============================================================================
c These arrays store final calculated doses for specific ingestion pathways
c at specific receptor locations

      REAL*8, ALLOCATABLE :: milkdose(:,:)          !radionuclide dose by milk (mrem/yr)
      REAL*8, ALLOCATABLE :: vegdose(:,:)           !radionuclide dose by vegetables (mrem/yr)
      REAL*8, ALLOCATABLE :: meatdose(:,:)          !radionuclide dose by meat (mrem/yr)

c ============================================================================
c SECTION 9: DYNAMIC ARRAY DECLARATIONS - METEOROLOGICAL DATA BY SECTOR
c ============================================================================
c These arrays store wind direction sector information (frequency, speed,
c source-receptor distance) for computing sector-dependent air concentrations

      REAL*8, ALLOCATABLE :: windbinfreq(:)       !Sector (wind dir) frequency from that direction
      REAL*8, ALLOCATABLE :: windbinspeed(:)      !Sector (wind dir) average speed from that direction (m/s)
      REAL*8, ALLOCATABLE :: windbinsrdist(:)     !Sector (wind dir) source-receptor distance (m)

c ============================================================================
c SECTION 10: DYNAMIC ARRAY DECLARATIONS - INDEXING & METADATA
c ============================================================================

      INTEGER, ALLOCATABLE :: nrads(:)                  !number of emitted radionucludes at each stack
      INTEGER RECEP(16)

      CHARACTER*8,  ALLOCATABLE :: radname(:,:)         !radionuclide names at each stack
      CHARACTER*3,  ALLOCATABLE :: windbinname(:)       !wind direction sector 
      CHARACTER*3,  ALLOCATABLE :: sectorname(:)        !sector name for receptors
      CHARACTER*3,  ALLOCATABLE :: vegsector(:)         !Sector (wind dir) location of vegetables
      CHARACTER*3,  ALLOCATABLE :: milksector(:)        !Sector (wind dir) location of cow/milk
      CHARACTER*3,  ALLOCATABLE :: meatsector(:)        !Sector (wind dir) location of meat
      CHARACTER*8,  ALLOCATABLE :: decaychain(:)        !Decay chain list
      CHARACTER*8,  ALLOCATABLE :: decaychainname(:,:,:)!Decay chain radionuclide names
      CHARACTER*8,  ALLOCATABLE :: xdecaychainname(:)   !Decay chain radionuclide names
      CHARACTER*3,  ALLOCATABLE :: verticalmomentum(:)  !vertical momentum (YES or NO)
      CHARACTER*3,  ALLOCATABLE :: maxsector(:)         !maximum downwind sector of max. conc.

c ============================================================================
c SECTION 11: STATIC VARIABLE DECLARATIONS - DISPERSION CALCULATION TERMS
c ============================================================================
c These variables are used in the Gaussian plume model calculations for
c atmospheric dispersion and effective stack height corrections

      CHARACTER*200 title,company,facility,address,poc,phone !facility metadata
      CHARACTER*200 windsource     !Source of wind data
      CHARACTER*200 winddates      !Date coverage of wind data
      CHARACTER*200 windloc        !Location of wind data
      CHARACTER*200 winddist       !Distance from wind data to source

      INTEGER nwindbins            !number of wind direction bins at each stack

      REAL*8 rainfall              !rainfall rate (m/yr)
      REAL*8 airtemp               !air temperature (C)
      REAL*8 pctcalm               !frequency of calm winds 

      REAL*8 FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF
      REAL*8 FACTORG,FACTORH,FACTORI,FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN
      REAL*8 area,xcrit,SIGMAZ,CAPSIG,pp,alldose,idose,dist,altconc
      REAL*8 tempdose,maxdose,maxdist,totmilkdose,totvegdose,itotvegdose,totmeatdose,itotcowdose
      REAL*8 capp,capn0,outprod,inprod,z1,term1,term2,sum
      REAL*8 transittime                              !transit time of air from emission release to receptor (s)
      REAL*8 wetdepadj                                !wet deposition adjustment factor
      REAL*8 xtime, xremove
      REAL*8 ingestion, immersion, inhalation, surface
      REAL*8 RHOAIR,RHOSTK,STKVEL,R,FM,DELHMO,DELHBU,HEFF,FB,zz,INHFAC,TINIT,T
      REAL*8 xmax, pmax 

c STATIC VARIABLE DECLARATIONS - CONTROL & LOOP INDICES

      INTEGER nstacks                               !number of unique stacks
      INTEGER receptor, N, M, I, J, K, L, II, JJ, P, Q                        !loop counters
      INTEGER nchain                                !number of nuclides in the decay series
      INTEGER maxsect
      INTEGER nrads_max, progeny_max

      CHARACTER*200 inline
      CHARACTER*8 NUC
      CHARACTER*3 DIAGFLAG

      integer date_time(8)                      !for date stamp function
      character*10 b(3)                         !for date stamp function

c ============================================================================
c SECTION 12: PARAMETER DEFINITIONS
c ============================================================================

      DOUBLE PRECISION, PARAMETER :: BREATH = 8000.D0 !m3/yr ! breathing rate
      DOUBLE PRECISION, PARAMETER :: DEPTIM = 100.D0  !years ! period of long-term buildup in soil
      DOUBLE PRECISION, PARAMETER :: REMCON = 0.015D0 !1/year! removal constant
      DOUBLE PRECISION, PARAMETER :: TWEATH = 12.D0   !days  ! weathering half-life
      DOUBLE PRECISION, PARAMETER :: DELVEG = 11.D0   !days  ! delay time harvest to consumption
      DOUBLE PRECISION, PARAMETER :: DELMILK = 2.D0   !days  ! delay time harvest to consumption
      DOUBLE PRECISION, PARAMETER :: DELMEAT = 17.D0  !days  ! delay time harvest to consumption
      DOUBLE PRECISION, PARAMETER :: CONVEG = 70.D0   !kg/y  ! human vegetable consumption
      DOUBLE PRECISION, PARAMETER :: CONMILK = 160.D0 !kg/y  ! human milk consumption
      DOUBLE PRECISION, PARAMETER :: CONMEAT = 75.D0  !kg/y  ! human meat consumption
      DOUBLE PRECISION, PARAMETER :: TEVEG = 60.0D0   !days  ! period of above ground exposure vegetables
      DOUBLE PRECISION, PARAMETER :: TEMILK = 30.0D0  !days  ! period of above ground exposure milk
      DOUBLE PRECISION, PARAMETER :: TEMEAT = 30.0D0  !days  ! period of above ground exposure meat
      DOUBLE PRECISION, PARAMETER :: FRVEG = 0.1D0    !frac  ! fraction of activity intercepted & retained by vegetables
      DOUBLE PRECISION, PARAMETER :: FRMILK = 1.8D0   !frac  ! fraction of activity intercepted & retained on forage or feed
      DOUBLE PRECISION, PARAMETER :: FRMEAT = 1.8D0   !frac  ! fraction of activity intercepted & retained on for
      DOUBLE PRECISION, PARAMETER :: YVEG = 1.0D0     !kg/m2 ! edible crop
      DOUBLE PRECISION, PARAMETER :: YMILK = 1.0D0    !kg/m2 ! edible crop
      DOUBLE PRECISION, PARAMETER :: YMEAT = 1.0D0    !kg/m2 ! edible crop
      DOUBLE PRECISION, PARAMETER :: PVEG = 220.D0    !kg/m2 ! edible crop
      DOUBLE PRECISION, PARAMETER :: PMILK = 220.D0   !kg/m2 ! edible crop
      DOUBLE PRECISION, PARAMETER :: PMEAT = 220.D0   !kg/m2 ! edible crop
      DOUBLE PRECISION, PARAMETER :: QMILK = 16.D0    !kg/day! feed or forage consumption rate for dairy
      DOUBLE PRECISION, PARAMETER :: QMEAT = 12.D0    !kg/day! feed or forate consumption rate for meat
      DOUBLE PRECISION, PARAMETER :: QWMILK = 60.D0   !kg/day! water consumption by dairy cow
      DOUBLE PRECISION, PARAMETER :: QWMEAT = 50.D0   !kg/day! water consumption by beef cattle
      DOUBLE PRECISION, PARAMETER :: FHVEG = 0.1D0    !frac  ! fraction of hydrogen in vegetables
      DOUBLE PRECISION, PARAMETER :: FHMILK = 0.1D0   !frac  ! fraction of hydrogen in milk
      DOUBLE PRECISION, PARAMETER :: FHMEAT = 0.11D0  !frac  ! fraction of hydrogen in meat
      DOUBLE PRECISION, PARAMETER :: FHFEED = 0.07D0  !frac  ! fraction of hydrogen in feed
      DOUBLE PRECISION, PARAMETER :: FCVEG = 0.090    !frac  ! fraction of carbon in vegetables
      DOUBLE PRECISION, PARAMETER :: FCMILK = 0.24D0  !frac  ! fraction of carbon in milk
      DOUBLE PRECISION, PARAMETER :: FCMEAT = 0.07D0  !frac  ! fraction of carbon in meat
      DOUBLE PRECISION, PARAMETER :: CH2O = 0.008D0   !kg/m3 ! concentration of water in air
      DOUBLE PRECISION, PARAMETER :: CCARB = 0.00016D0!kg/m3 ! concentration of carbon in air

c Physical constants

       REAL*8, PARAMETER :: PI = 3.14159265358979d0
c
c Wind and stack parameter limits

       REAL*8, PARAMETER :: WIND_MIN = 0.01d0
       REAL*8, PARAMETER :: FLOWRATE_DEFAULT = 0.3d0

c Wind direction and downwind receptor relationships
!       DATA FROM   / 'N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW' /
!       DATA TOWARD / 'S','SSW','SW','WSW','W','WNW','NW','NNW','N','NNE','NE','ENE','E','ESE','SE','SSE' /
!       DATA WINDS  /   1,    2,   3,    4,  5,    6,   7,    8,  9,   10,  11,   12, 13,   14,  15,  16 /
       DATA RECEP  /    9,   10,  11,   12, 13,   14,  15,   16,  1,    2,   3,    4,  5,    6,   7,   8 /

c ============================================================================
c SECTION 14: FILE OPERATIONS - OPEN INPUT AND OUTPUT FILES
c ============================================================================
c Open the necessary input files for model parameters and output files for
c results reporting

      write(*,*)
      write(*,*)'COMPLY v2.0'
      write(*,*)
 
      open(11,file='input.level4.txt',status='old',form='formatted')
      open(12,file='output.level4.txt',status='unknown',form='formatted')
      open(13,file='diag.level4.txt',status='unknown',form='formatted')

      open(9,file='comply171_original_input_data.csv',status='old',form='formatted')
      open(10,file='comply171_original_decaychain_data.csv',status='old',form='formatted')

c ============================================================================
c SECTION 15: INPUT PARSING - READ FACILITY IDENTIFICATION AND STACK COUNT
c ============================================================================
c Read project metadata and determine the number of stacks to process

       read(11,'(38x,a)') title
       read(11,'(38x,a)') company
       read(11,'(38x,a)') facility
       read(11,'(38x,a)') address
       read(11,'(38x,a)') poc
       read(11,'(38x,a)') phone
       read(11,'(38x,a3)') diagflag
       read(11,'(38x,i3)') nstacks

       read(11,'(a)') !blank spacer to improve readability of the input file

       read(11,'(38x,a)') windsource
       read(11,'(38x,a)') winddates
       read(11,'(38x,a)') windloc
       read(11,'(38x,a)') winddist
       read(11,'(38x,f10.5)') airtemp
       read(11,'(38x,f10.5)') rainfall
       read(11,'(38x,i3)') nwindbins

       allocate(windbinname(nwindbins))
       allocate(windbinfreq(nwindbins))
       allocate(windbinspeed(nwindbins))
       allocate(windbinsrdist(nwindbins))
       allocate(sectorname(nwindbins))
       allocate(sectordose(nwindbins))

       sum = 0.0
       do m = 1, nwindbins
       read(11,'(38x,a)') inline
       read(inline,*) windbinname(m),windbinfreq(m),windbinspeed(m)
       sum = sum + windbinfreq(m)
       enddo !end loop over reading in wind data and source-receptor distances

       airtemp = airtemp * (9.D0/5.D0) + 32.D0 !convert C to F
       rainfall = rainfall / 100.D0 !convert cm to m

       !QA checks on the wind inputs
       if(nwindbins.eq.1) then
        windbinfreq(1) = 0.25
       endif

       if(sum.lt.1.0.and.nwindbins.gt.1) then
        pctcalm = 1 - sum
       elseif(sum.gt.1.0.and.nwindbins.gt.1) then
        write(*,*)'ERROR! Sum of wind frequencies exceed 1.0'
        stop
       else
        pctcalm = 0.0
       endif
        write(*,*)'Percent calms = ',pctcalm  

!        do m = 1, nwindbins
!        write(*,*)windbinname(m),windbinfreq(m),windbinspeed(m),windbinsrdist(m),sectorname(m)
!        enddo

       read(11,'(a)') !blank spacer to improve readability of the input file

       if(nwindbins.eq.1) then !simple conservative estimate of winds no wind rose
       read(11,'(38x,a)') inline
       read(inline,*) sectorname(1),windbinsrdist(1)
       else !full 16 direction wind rose
       do m = 1, nwindbins
       read(11,'(38x,a)') inline
       read(inline,*) sectorname(RECEP(m)),windbinsrdist(RECEP(m))
       enddo !end loop over reading in wind data and source-receptor distances
       endif

c ============================================================================
c SECTION 16: MEMORY ALLOCATION - ALLOCATE STACK-LEVEL PARAMETER ARRAYS
c ============================================================================
c Allocate arrays to store stack-specific parameters (heights, distances,
c temperatures, etc.) based on the number of stacks

       allocate(latitude(nstacks))
       allocate(longitude(nstacks))
       allocate(releaseht(nstacks))
       allocate(bldheight(nstacks))
       allocate(bldwidth(nstacks))
       allocate(flowrate(nstacks))
       allocate(diameter(nstacks))
       allocate(stacktemp(nstacks))
       allocate(verticalmomentum(nstacks))
       allocate(distancemei(nstacks))
       allocate(distanceveg(nstacks))
       allocate(distancemilk(nstacks))
       allocate(distancemeat(nstacks))
       allocate(vegsector(nstacks))
       allocate(milksector(nstacks))
       allocate(meatsector(nstacks))
       allocate(nrads(nstacks))
       allocate(maxsector(nstacks))

       nrads_max = 3000 !maximum number of radionuclides possible per stack

       allocate(radname(nstacks,nrads_max))
       allocate(rademis(nstacks,nrads_max))
       allocate(radchainlen(nstacks,nrads_max))
       allocate(maxconc(nstacks,nrads_max))
       allocate(totimdose(nstacks,nrads_max))
       allocate(totihdose(nstacks,nrads_max))
       allocate(totsurfdose(nstacks,nrads_max))
       allocate(totingdose(nstacks,nrads_max))
       allocate(milkdose(nstacks,nrads_max))
       allocate(vegdose(nstacks,nrads_max))
       allocate(meatdose(nstacks,nrads_max))

       progeny_max = 19 !maximum number of progeny for any radionuclide

       allocate(decaychainname(nstacks,nrads_max,progeny_max))
       allocate(halflife(nstacks,nrads_max,progeny_max))
       allocate(immersionfac(nstacks,nrads_max,progeny_max))
       allocate(inhalationfac(nstacks,nrads_max,progeny_max))
       allocate(ingestionfac(nstacks,nrads_max,progeny_max))
       allocate(surfacefac(nstacks,nrads_max,progeny_max))
       allocate(drydepvel(nstacks,nrads_max,progeny_max))
       allocate(xbveg(nstacks,nrads_max,progeny_max))
       allocate(xbfor(nstacks,nrads_max,progeny_max))
       allocate(xfmilk(nstacks,nrads_max,progeny_max))
       allocate(xfmeat(nstacks,nrads_max,progeny_max))

       allocate(decaychain(19))
       allocate(sourceterm(19))
       allocate(initialamt(19))
       allocate(activity(19))
       allocate(surfaceactivity(19))
       allocate(ambientactivity(19))
       allocate(normsourceactivity(19))
       allocate(kvalue(19))
       allocate(ktot(19))
       allocate(imdose(19))
       allocate(ihdose(19))
       allocate(surfdose(19))
       allocate(aveglf(19))
       allocate(avegrt(19))
       allocate(amilklf(19))
       allocate(amilkrt(19))
       allocate(ameatlf(19))
       allocate(ameatrt(19))
       allocate(cfveg(19))
       allocate(cfmilk(19))
       allocate(cfmeat(19))
       allocate(dveg(19))
       allocate(dmilk(19))
       allocate(dmeat(19))
       allocate(ingdose(19))
       allocate(xhalflife(19))
       allocate(dosefactor(19))
       allocate(xdrydepvel(19))
       allocate(xdecaychainname(19))
       allocate(ybveg(19))
       allocate(ybfor(19))
       allocate(yfmilk(19))
       allocate(yfmeat(19))

c ============================================================================
c SECTION 17: INPUT PARSING LOOP - READ STACK-SPECIFIC PARAMETERS
c ============================================================================
c Main loop to read input parameters for each stack source. Includes
c geometrical parameters, temperatures, distances to receptors, and
c number of emitted radionuclides

       do n = 1, nstacks 

       read(11,'(a)') !blank spacer to improve readability of the input file
       read(11,'(38x,f10.5)') latitude(n)
       read(11,'(38x,f10.5)') longitude(n)
       read(11,'(38x,f10.5)') releaseht(n)
       read(11,'(38x,f10.5)') bldheight(n)
       read(11,'(38x,f10.5)') bldwidth(n)
       read(11,'(38x,f10.5)') flowrate(n)
       read(11,'(38x,f10.5)') diameter(n)
       read(11,'(38x,f10.5)') stacktemp(n)
       read(11,'(38x,a)') verticalmomentum(n)
       read(11,'(38x,f10.5)') distanceveg(n)
       read(11,'(38x,f10.5)') distancemilk(n)
       read(11,'(38x,f10.5)') distancemeat(n)
       read(11,'(38x,a)') vegsector(n)
       read(11,'(38x,a)') milksector(n)
       read(11,'(38x,a)') meatsector(n)

       read(11,'(38x,i5)') nrads(n)

       stacktemp(n) = stacktemp(n) * (9.D0/5.D0) + 32.D0 !convert C to F

       if(nwindbins.eq.1) then !if simple test make sure all downwind receptor sectors are the same
        vegsector(1) = sectorname(1)
        milksector(1) = sectorname(1)
        meatsector(1) = sectorname(1)
       endif

c ============================================================================
c SECTION 19: INPUT PARSING - READ RADIONUCLIDE EMISSION DATA
c ============================================================================
c For each stack, read the list of radionuclides emitted and their
c emission rates (in Ci/s)

       do m = 1, nrads(n)
       read(11,'(38x,a)') inline
       read(inline,*) NUC,rademis(n,m)
       rademis(n,m) = rademis(n,m) * (1./365.)*(1./24.)*(1./60.)*(1./60.) !year to seconds

       !convert nuclide name to all uppercase
       do i = 1, len(NUC)
        j = iachar(NUC(i:i))
        if(j.ge.iachar("a").and.j.le.iachar("z")) then
                radname(n,m)(i:i) = achar(iachar(NUC(i:i))-32)
        else
                radname(n,m)(i:i) = NUC(i:i)
        endif
       enddo

        radchainlen(n,m) = 1 !This can change to a higher number if a decay chain is found for this nuclide
        !Fill arrays with zeros in case there is no nuclide specific data
        do l = 1, 19
        decaychainname(n,m,1) = "NONE"
        halflife(n,m,1) = 0.0
        inhalationfac(n,m,1) = 0.0
        immersionfac(n,m,1) = 0.0
        surfacefac(n,m,1) = 0.0
        ingestionfac(n,m,1) = 0.0
        xbveg(n,m,1) = 0.0
        xbfor(n,m,1) = 0.0
        xfmilk(n,m,1) = 0.0
        xfmeat(n,m,1) = 0.0
        drydepvel(n,m,1) = 0.0
        enddo
        decaychainname(n,m,1) = radname(n,m)

       rewind(9)
       rewind(10)
       read(9,'(a)') inline

       !Find decay chain and nuclide specific data
        do ii = 1, 879
         read(10,'(a)') inline
         read(inline,*) NUC,k,(decaychain(l),l=1,k)
         if(NUC.eq.radname(n,m)) then !found decay chain, get nuclide specific inputs

          do jj = 1, 826
          read(9,'(a)') inline
          read(inline,*) NUC,FACTOR,FACTORB,FACTORC,
     &      FACTORD,FACTORE,FACTORF,FACTORG,FACTORH,FACTORI,FACTORJ,
     &      FACTORK,FACTORL,FACTORM,FACTORN

           do l = 1, k
           if(NUC.eq.decaychain(l)) then
            decaychainname(n,m,l) = NUC
            halflife(n,m,l) = FACTOR*365. !! convert years to days
            inhalationfac(n,m,l) = FACTORF
            immersionfac(n,m,l) = FACTORG
            surfacefac(n,m,l) = FACTORH
            ingestionfac(n,m,l) = FACTORI
            xbveg(n,m,l) = FACTORJ
            xbfor(n,m,l) = FACTORK
            xfmilk(n,m,l) = FACTORL
            xfmeat(n,m,l) = FACTORM
            drydepvel(n,m,l) = FACTORN
           endif
           enddo !loop over decay chain
          enddo !loop over all nuclides on input reference table
          radchainlen(n,m) = k
         endif
         enddo !loop over all nuclides on decay chain input table

       enddo !end loop over reading in radionuclide emissions

       enddo !end loop over reading in stacks


!       do n = 1, nstacks
!       do m = 1, nrads(n)
!       write(*,*)'===================================='
!       write(*,*)'STACK:',n,' NUCLIDE #:',m
!       write(*,*)'NUCLIDE:',radname(n,m)
!       write(*,*)'EMISSIONS:',rademis(n,m)
!       write(*,*)'CHAIN LENGTH:',radchainlen(n,m)
!       do l = 1, radchainlen(n,m)
!       write(*,*)'  ',decaychainname(n,m,l),halflife(n,m,l),inhalationfac(n,m,l),drydepvel(n,m,l)
!       write(*,*)'          ',surfacefac(n,m,l),ingestionfac(n,m,l),xfmeat(n,m,l)
!       enddo
!       enddo
!       enddo

c ============================================================================
c SECTION 20: MAIN CALCULATION LOOP - AIR CONCENTRATIONS & DOSES
c ============================================================================
c This is the core calculation loop that computes atmospheric dispersion,
c radioactive decay, and dose calculations for all combinations of:
c   - Stacks (n)
c   - Radionuclides (m)
c   - Wind direction sectors (j)
c   - Receptor distances (receptor loop: 10 to 3010 m in 10 m increments)
c ============================================================================


       do n = 1, nstacks

       alldose = 0.
       idose = 0.

       do j = 1, nwindbins
       sectordose(j) = 0.0
       enddo

c ============================================================================
c SECTION 21: SURFACE RELEASE SPECIAL CASE
c ============================================================================
c If the release height is zero (ground-level/area source), clear the
c building dimensions since COMPLY assumes no release shape geometry for
c surface releases

       if(releaseht(n).eq.0.0) then !surface release
       bldwidth(n) = 0.0 !odd that comply assumes no release shape geometry for a surface (area) release
       bldheight(n) = 0.0
       endif

c ============================================================================
c SECTION 22: INITIALIZE DISPERSION PARAMETERS
c ============================================================================
c Calculate the critical distance (xcrit) which represents the crossover
c point between near-field (simple) and far-field (Gaussian plume) models

       xcrit = 2.5 * sqrt(bldwidth(n) * bldheight(n))
       write(*,*)'Critical distance (m) is ',xcrit

       do m = 1, nrads(n)

c ============================================================================
c SECTION 23: INITIALIZE DOSE ACCUMULATORS FOR CURRENT RADIONUCLIDE
c ============================================================================

       maxdose = 0.0
       inhalation = 0.
       immersion = 0.
       ingestion = 0.
       surface = 0.
       totmilkdose = 0.
       totvegdose = 0.
       totmeatdose = 0.

c ============================================================================
c SECTION 24: WIND SECTOR AND DISTANCE ITERATION LOOPS
c ============================================================================
c Loop over wind direction sectors and receptor distances. For each combination,
c calculate air concentration and dose contributions across all decay pathways


       do j = 1, nwindbins

       do receptor = 1, 4, 1

       if(receptor.eq.1) then
       dist = windbinsrdist(j)
       elseif(receptor.eq.2) then
       dist = distanceveg(n)
       elseif(receptor.eq.3) then
       dist = distancemilk(n)
       elseif(receptor.eq.4) then
       dist = distancemeat(n)
       else
       endif


       !Air concentration
       CALL calc_air_concentration( dist, pctcalm, xcrit, releaseht(n), diameter(n), 
     & flowrate(n), stacktemp(n), bldheight(n), bldwidth(n), verticalmomentum(n),
     & windbinfreq(j), windbinspeed(j), rademis(n,m), altconc, airtemp, transittime)


       !Wet deposition
       CALL wetdep_calc(dist,pctcalm,rainfall,windbinfreq(j),windbinspeed(j),rademis(n,m),altconc,wetdepadj)


c ============================================================================
c DOSE CALCULATIONS
c ============================================================================

         nchain = int(radchainlen(n,m))

         do l = 1, nchain
         imdose(l) = 0.0
         ihdose(l) = 0.0
         surfdose(l) = 0.0
         xdrydepvel(l) = drydepvel(n,m,l)
         xhalflife(l) = halflife(n,m,l)
         xdecaychainname(l) = decaychainname(n,m,l)
         enddo



         !Air immersion
         do l = 1, nchain
         dosefactor(l) = immersionfac(n,m,l)
         enddo


         CALL immersion_calc(progeny_max,xhalflife,dosefactor,altconc,transittime,nchain,imdose,ambientactivity)


         immersion = 0.0
         do l = 1, nchain
         immersion = immersion + imdose(l)
         enddo


         !Inhalation
         do l = 1, nchain
         dosefactor(l) = inhalationfac(n,m,l)
         enddo

         CALL inhalation_calc(BREATH,progeny_max,xhalflife,dosefactor,altconc,transittime,nchain,ihdose,ambientactivity)

         inhalation = 0.0
         do l = 1, nchain
         inhalation = inhalation + ihdose(l)
         enddo


         !Surface
         do l = 1, nchain 
         dosefactor(l) = surfacefac(n,m,l)
         enddo

         CALL surface_calc(wetdepadj,progeny_max,xdrydepvel,xhalflife,xdecaychainname,
     &        dosefactor,REMCON,DEPTIM,nchain,surfdose,ambientactivity)

         surface = 0.0
         do l = 1, nchain
         surface = surface + surfdose(l)
         enddo


         !Ingestion
         do l = 1, nchain
         dosefactor(l) = ingestionfac(n,m,l)
         ybveg(l) = xbveg(n,m,l)
         ybfor(l) = xbfor(n,m,l)
         yfmilk(l) = xfmilk(n,m,l)
         yfmeat(l) = xfmeat(n,m,l)
         enddo

         CALL ingestion_calc(wetdepadj,transittime,progeny_max,xdrydepvel,xhalflife,
     &    xdecaychainname,dosefactor,REMCON,DEPTIM,nchain,altconc,dveg,dmeat,dmilk,
     &    ybveg,ybfor,yfmilk,yfmeat,
     &    TWEATH,TEVEG,TEMILK,TEMEAT,FHVEG,FHFEED,CH2O,QMILK,QMEAT,FHMILK,FHMEAT,FCVEG,FCMILK,
     &    FCMEAT,CCARB,FRVEG,FRMILK,FRMEAT,YVEG,YMILK,YMEAT,PVEG,PMILK,PMEAT,QWMILK,QWMEAT,
     &    DELVEG,DELMILK,DELMEAT,CONVEG,CONMILK,CONMEAT)

          totvegdose = 0.0
          totmilkdose = 0.0
          totmeatdose = 0.0
          do l = 1, nchain
          totvegdose = totvegdose + dveg(l)
          totmilkdose = totmilkdose + dmilk(l)
          totmeatdose = totmeatdose + dmeat(l)
          enddo


c ============================================================================
c REPORT KEY RECEPTOR LOCATIONS
c ============================================================================
c For specific receptor locations defined in input (vegetable garden,
c milk cow, meat animal locations), report the computed doses. Also
c accumulate total dose across all pathways and locations.
c ============================================================================


       if(receptor.eq.2.and.vegsector(n).eq.sectorname(j)) then !vegetable receptor
        alldose = alldose + totvegdose
        vegdose(n,m) = totvegdose
        write(*,*)'Vegetable dose:',j,radname(n,m),vegdose(n,m),dist,sectorname(j)
        endif
       
       if(receptor.eq.3.and.milksector(n).eq.sectorname(j)) then !milk receptor
        alldose = alldose + totmilkdose
        milkdose(n,m) = totmilkdose
        write(*,*)'Milk dose:',j,radname(n,m),milkdose(n,m),dist,sectorname(j)
       endif

       if(receptor.eq.4.and.meatsector(n).eq.sectorname(j)) then !meat receptor
        alldose = alldose + totmeatdose
        meatdose(n,m) = totmeatdose
        write(*,*)'Meat dose:',j,radname(n,m),meatdose(n,m),dist,sectorname(j)
       endif
 
c ============================================================================
c CALCULATE MAXIMUM DOSE RECEPTOR
c ============================================================================
c Find the receptor location (distance and sector) receiving the maximum
c direct exposure dose (immersion + inhalation + surface deposition).
c This is the critical dose location for regulatory comparison.
c ============================================================================


       if(receptor.eq.1) then 
            tempdose = immersion + inhalation + surface
        if(tempdose.gt.maxdose) then 
            maxdose = tempdose
            maxsector(n) = sectorname(j)
            distancemei(n) = dist
            totsurfdose(n,m) = surface
            totimdose(n,m) = immersion
            totihdose(n,m) = inhalation
            maxconc(n,m) = altconc
        endif
       endif

       enddo !end loop over distances between source and receptor (receptor)

       !make a sector total for output report
       sectordose(j) = sectordose(j) + 
     &  immersion+inhalation+surface+totvegdose+totmilkdose+totmeatdose

       enddo !end loop over number of wind direction bins (j)

c ============================================================================
c CONSOLE OUTPUT - SUMMARY OF DOSES FOR CURRENT RADIONUCLIDE
c ============================================================================
c Print dose summaries by pathway (inhalation, immersion, surface, ingestion)
c and total dose at the maximum dose receptor location.
c ============================================================================

 2999  FORMAT(a30,a10,f30.15)
 3100  FORMAT(a,f15.4)

       if(radname(n,m)(1:2).eq.'I-') then
       idose = idose + totihdose(n,m)+totimdose(n,m)+totsurfdose(n,m)+
     &  meatdose(n,m)+milkdose(n,m)+vegdose(n,m)
       endif

       alldose = alldose + totihdose(n,m)+totimdose(n,m)+totsurfdose(n,m)
       write(*,*)'Max dose (',maxdose,') receptor at ',distancemei(n),maxsector(n)
       write(*,*)
       write(*,2999)'Maximum air concentration:',radname(n,m),maxconc(n,m)
       write(*,2999)'Inhalation dose:',radname(n,m),totihdose(n,m)
       write(*,2999)'Immersion dose:',radname(n,m),totimdose(n,m)
       write(*,2999)'Surface dose:',radname(n,m),totsurfdose(n,m)
       write(*,2999)'Ingestion meat dose:',radname(n,m),meatdose(n,m)
       write(*,2999)'Ingestion milk dose:',radname(n,m),milkdose(n,m)
       write(*,2999)'Ingestion vegetables dose:',radname(n,m),vegdose(n,m)
       write(*,2999)'Total ingestion dose:',radname(n,m),meatdose(n,m)+milkdose(n,m)+vegdose(n,m)
       write(*,2999)'Total dose:',radname(n,m),(totihdose(n,m)+totimdose(n,m)+totsurfdose(n,m)+
     &  meatdose(n,m)+milkdose(n,m)+vegdose(n,m))
       write(*,*)'--------------------------------------------------------------------------' 

       enddo !end loop over radionuclides (m)
       enddo !end loop over stacks (n)

       write(*,*)
       write(*,3100)'Total dose (mrem/yr):',alldose
       write(*,*)'Total dose limit is 10.0 mrem/yr'
       write(*,*)
       write(*,3100)'Iodine total dose (mrem/yr):',idose
       write(*,*)'Iodine limit is 3.0 mrem/yr'
       write(*,*)

c ============================================================================
c SUMMARY OUTPUT
c ============================================================================

 4000 FORMAT('COMPLY version 2.0 level 4 Dose Summary Report')
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
 1100 FORMAT('Effective dose equivalent all stacks (mrem/year):',f10.4)
       write(12,1102)
       write(12,2002)
       write(12,1100) alldose

 1111 FORMAT('Effective dose equivalent all stacks from iodine (mrem/year):',f10.4)
       write(12,1111) idose

       write(12,2002)

 1101 FORMAT('Maximum dose receptor was located at ',f8.0,' m downwind in sector ',a)

       write(12,2002)
       write(12,2002)

       do n = 1, nstacks 
       write(12,1006) n
       write(12,1101) distancemei(n),trim(maxsector(n))
       write(12,2002)
       write(12,1997)
       write(12,1998)
       do m = 1, nrads(n)

 1997 FORMAT('          ..........................Dose reported in mrem/year..............................    (Ci/yr)   (pCi/m3)')
 1998 FORMAT('Nuclide , Total Dose, Inhalation,  Immersion,    Surface,       Meat,       Milk,  Vegetable, Emissions ,  Concen')
 1999 FORMAT(a,9(',',ES11.3))

       write(12,1999) radname(n,m),
     &    (totimdose(n,m)+totihdose(n,m)+totsurfdose(n,m)+meatdose(n,m)+milkdose(n,m)+vegdose(n,m)),
     &    totihdose(n,m),totimdose(n,m),totsurfdose(n,m),meatdose(n,m),milkdose(n,m),vegdose(n,m),
     &    rademis(n,m)*60.*60.*24.*365.,
     &    1.0E12*maxconc(n,m)

       enddo
       write(12,2002)
       enddo
       
       
 1030 FORMAT('STACK INFORMATION')
 1006 FORMAT('Stack number                   :',i3)
 1007 FORMAT('Stack latitude (degrees)       :',f10.5)
 1008 FORMAT('Stack longitude (degrees)      :',f10.5)
 1009 FORMAT('Release height (m)             :',f10.1)
 1010 FORMAT('Stack diameter (m)             :',f10.1)
 1011 FORMAT('Flow rate (m3/s)               :',f10.1)
 1012 FORMAT('Stack temperature (F)          :',f10.1)
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
       write(12,1015) bldheight(n)
       write(12,1016) bldwidth(n)
       write(12,1019) nrads(n)
       write(12,1033) distanceveg(n),vegsector(n)
       write(12,1034) distancemilk(n),milksector(n)
       write(12,1035) distancemeat(n),meatsector(n)
       write(12,2002)
       enddo


 1130 FORMAT('RADIONUCLIDE AND DECAY CHAIN DATA')
 1131 FORMAT('Nuclide and chain length:',a10,1x,f3.0)
 1132 FORMAT(a10,6(1x,ES9.1))
 1133 FORMAT(10x,1x,f10.2,1x,f10.2,1x,f10.2)
 1134 FORMAT('  Nuclide   Halflife  DryDepVel  Inhalat.  Immers.   Surface  Ingestion')
 1135 FORMAT('Halflife expressed in days and dry deposition in m/day')
 1136 FORMAT('Immersion dose coefficients expressed as mrem/year per Ci/m3')
 1137 FORMAT('Surface dose coefficients expressed as mrem/year per Ci/m2')
 1138 FORMAT('Inhalation and ingestion dose coefficients expressed as mrem/Ci')

      write(12,2001)
      write(12,1130)
      write(12,2002)
      write(12,1135)
      write(12,1136)
      write(12,1137)
      write(12,1138)
      write(12,2002)

       do n = 1, nstacks
       write(12,2002)
       do m = 1, nrads(n)
        write(12,1131) radname(n,m),radchainlen(n,m)
        write(12,1134)
       do l = 1, radchainlen(n,m)
        write(12,1132) decaychainname(n,m,l),halflife(n,m,l),drydepvel(n,m,l),
     &    inhalationfac(n,m,l),immersionfac(n,m,l),surfacefac(n,m,l),ingestionfac(n,m,l)
       enddo
       write(12,2002)
       enddo
       enddo
 

 1025 FORMAT('METEOROLOGY')
 1018 FORMAT('Rainfall rate (m/yr)         :',f10.1)
 1013 FORMAT('Air temperature (F)          :',f10.1)
 1022 FORMAT('Calm winds (fraction)        :',f10.3)
 1023 FORMAT('Sector winds from,frequency(fraction),speed(m/s)')
 1024 FORMAT(a3,14x,',',4x,f10.3,5x,',',f10.3)
 1026 FORMAT('Sector winds to  ,downwind receptor(m),sector total dose(mrem/yr)')
 1027 FORMAT(a3,14x,',',5x,f10.3,5x,',',ES9.2)

       write(12,2002)
       write(12,2001)
       write(12,1025)
       write(12,2002)
       write(12,1018) rainfall
       write(12,1013) airtemp
       write(12,1022) pctcalm
       write(12,2002)
       write(12,1023)

        do m = 1, nwindbins
        write(12,1024)windbinname(m),windbinfreq(m),windbinspeed(m)
        enddo

        write(12,2002)
        write(12,1026)

        do m = 1, nwindbins
        write(12,1027)sectorname(m),windbinsrdist(m),sectordose(m)
        enddo


c ============================================================================
c EXECUTION TIMESTAMP
c ============================================================================
c Capture and write execution date and time for traceability and
c regulatory documentation purposes
c ============================================================================

       call date_and_time(b(1), b(2), b(3), date_time)

 2100 FORMAT('Date/Time (LST) stamp of code execution (MM/DD/YYYY HH:MM) --> ',i2.2,'/',i2.2,'/',i4,' ',i2.2,':',i2.2)

       write(12,2001)
       write(12,2002)
       write(12,2100) date_time(2),date_time(3),date_time(1),
     &           date_time(5),date_time(6)

       write(*,2002)
       write(*,2100) date_time(2),date_time(3),date_time(1),
     &           date_time(5),date_time(6)


c ============================================================================
c DIAGNOSTIC OUTPUT
c ============================================================================
c Write comma-separated values to diagnostic file for all distance/sector
c combinations, including concentrations and all dose components
c
c Concentration in pCi/m3 and dose in mrem/year
c ============================================================================

       if(diagflag.eq.'YES'.or.diagflag.eq.'TRUE'.or.
     &    diagflag.eq.'yes'.or.diagflag.eq.'true') then

       write(13,3002) !write units to header
       write(13,3001) !write header to file

       do n = 1, nstacks
       do m = 1, nrads(n)

       inhalation = 0.
       immersion = 0.
       ingestion = 0.
       surface = 0.
       totmilkdose = 0.
       totvegdose = 0.
       totmeatdose = 0.

       xcrit = 0.0
       if(releaseht(n).eq.0.0) then !surface release
       bldwidth(n) = 0.0 !odd that comply assumes no release shape geometry for a surface (area) release
       bldheight(n) = 0.0
       endif

       do j = 1, nwindbins
       do receptor = 25, 50000, 25
       dist = float(receptor)

       !Air concentration
       CALL calc_air_concentration( dist, pctcalm, xcrit, releaseht(n), diameter(n),
     & flowrate(n), stacktemp(n), bldheight(n), bldwidth(n), verticalmomentum(n),
     & windbinfreq(j), windbinspeed(j), rademis(n,m), altconc, airtemp, transittime)


       !Wet deposition
       CALL wetdep_calc(dist,pctcalm,rainfall,windbinfreq(j),windbinspeed(j),rademis(n,m),altconc,wetdepadj)

         nchain = int(radchainlen(n,m))

         do l = 1, nchain
         imdose(l) = 0.0
         ihdose(l) = 0.0
         surfdose(l) = 0.0
         xdrydepvel(l) = drydepvel(n,m,l)
         xhalflife(l) = halflife(n,m,l)
         xdecaychainname(l) = decaychainname(n,m,l)
         enddo

         !Air immersion
         do l = 1, nchain
         dosefactor(l) = immersionfac(n,m,l)
         enddo

         CALL immersion_calc(progeny_max,xhalflife,dosefactor,altconc,transittime,nchain,imdose,ambientactivity)

         immersion = 0.0
         do l = 1, nchain
         immersion = immersion + imdose(l)
         enddo


         !Inhalation
         do l = 1, nchain
         dosefactor(l) = inhalationfac(n,m,l)
         enddo

         CALL inhalation_calc(BREATH,progeny_max,xhalflife,dosefactor,altconc,transittime,nchain,ihdose,ambientactivity)

         inhalation = 0.0
         do l = 1, nchain
         inhalation = inhalation + ihdose(l)
         enddo

         !Surface
         do l = 1, nchain
         dosefactor(l) = surfacefac(n,m,l)
         enddo


         CALL surface_calc(wetdepadj,progeny_max,xdrydepvel,xhalflife,xdecaychainname,
     &        dosefactor,REMCON,DEPTIM,nchain,surfdose,ambientactivity)

         surface = 0.0
         do l = 1, nchain
         surface = surface + surfdose(l)
         enddo

         !Ingestion
         do l = 1, nchain
         dosefactor(l) = ingestionfac(n,m,l)
         ybveg(l) = xbveg(n,m,l)
         ybfor(l) = xbfor(n,m,l)
         yfmilk(l) = xfmilk(n,m,l)
         yfmeat(l) = xfmeat(n,m,l)
         enddo

         CALL ingestion_calc(wetdepadj,transittime,progeny_max,xdrydepvel,xhalflife,
     &    xdecaychainname,dosefactor,REMCON,DEPTIM,nchain,altconc,dveg,dmeat,dmilk,
     &    ybveg,ybfor,yfmilk,yfmeat,
     &    TWEATH,TEVEG,TEMILK,TEMEAT,FHVEG,FHFEED,CH2O,QMILK,QMEAT,FHMILK,FHMEAT,FCVEG,FCMILK,
     &    FCMEAT,CCARB,FRVEG,FRMILK,FRMEAT,YVEG,YMILK,YMEAT,PVEG,PMILK,PMEAT,QWMILK,QWMEAT,
     &    DELVEG,DELMILK,DELMEAT,CONVEG,CONMILK,CONMEAT)

          totvegdose = 0.0
          totmilkdose = 0.0
          totmeatdose = 0.0
          do l = 1, nchain
          totvegdose = totvegdose + dveg(l)
          totmilkdose = totmilkdose + dmilk(l)
          totmeatdose = totmeatdose + dmeat(l)
          enddo

 3000  FORMAT(i2,',',a,',',a,',',f9.1,',',ES10.3,7(',',ES10.3))
 3001  FORMAT('STACK,NUCLIDE,SECTOR,DISTANCE,CONC,TOTALDOSE,IMMERSION,INHALATION,SURFACE,MEAT,MILK,VEG')
 3002  FORMAT('# Distance as m; Concentration as Ci/m3; all dose units are mrem/year')

       write(13,3000) n,radname(n,m),sectorname(j),dist,
     &     altconc,
     &     (immersion+inhalation+surface+totmeatdose+totmilkdose+totvegdose),
     &     immersion,inhalation,surface,totmeatdose,totmilkdose,totvegdose

        enddo !receptors
        enddo !wind bins
        enddo !nuclides
        enddo !stacks

        endif !end condition for writing out diagnostic output

       END PROGRAM

c ============================================================================
c START OF SUBROUTINES
c ============================================================================

       SUBROUTINE calc_air_concentration(dist, pctcalm, xcrit, 
     & xreleaseht, xdiameter, xflowrate, xstacktemp, xbldheight, xbldwidth, xverticalmomentum, 
     & xwindbinfreq, xwindbinspeed, xrademis, xconc, airtemp, transittime)

       IMPLICIT NONE

       REAL*8  :: dist, pctcalm, xcrit, airtemp
       REAL*8  :: xreleaseht, xdiameter, xflowrate, xstacktemp
       REAL*8  :: xbldheight, xbldwidth, xwindbinfreq, xwindbinspeed
       REAL*8  :: xrademis
       CHARACTER*3 :: xverticalmomentum

       REAL*8 :: xconc
       REAL*8 :: HEFF, SIGMAZ, CAPSIG, transittime, pp

       ! Local variables
       REAL*8 :: RHOAIR, RHOSTK, STKVEL, R, FM, DELHMO, DELHBU
       REAL*8 :: FB, AREA, zz, TERM1, TERM2
       REAL*8 :: xmax, pmax
        
       REAL*8, PARAMETER :: PI = 3.14159265358979D0


       IF (xflowrate .le. 0.0) xflowrate = 0.3 !default flow rate

       ! ------------------------------------------------------------
       ! Plume rise
       ! ------------------------------------------------------------

       HEFF = xreleaseht

       IF (xreleaseht .GT. 2.5D0 * xbldheight .or. xreleaseht.eq. 0.0D0) THEN

       ! Temperature and density
       RHOAIR = 0.0735D0 * (80.0D0 + 460.0D0) / (airtemp + 460.0D0)
       RHOAIR = RHOAIR * 454.0D0 / (0.3048D0 ** 3)
       RHOSTK = RHOAIR * (airtemp + 460.0D0) / (xstacktemp + 460.0D0)

       ! Momentum rise
       AREA   = PI * (xdiameter ** 2.0D0) / 4.0D0
       STKVEL = xflowrate / AREA
       R      = STKVEL / xwindbinspeed
       FM     = (RHOSTK / RHOAIR) * STKVEL * xflowrate / PI

       IF (xverticalmomentum .EQ. 'YES') THEN
         TERM1  = (75.0D0 * FM / (((1.0D0 + 3.0D0 / R) * xwindbinspeed) ** 2.0D0))
         TERM2  = (xreleaseht + 7.79D0 * SQRT(FM) / ((1.0D0 + 3.0D0 / R) * xwindbinspeed))

        DELHMO = 0.93D0 * (TERM1) ** (3.0D0 / 7.0D0) * 
     &                    (TERM2) ** (1.0D0 / 7.0D0)
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


       IF (xreleaseht .GT. 2.5D0 * xbldheight .or. xreleaseht.eq. 0.0D0) THEN
       !pp = P = PP(dist,HEFF,0.0) for tall stack or stack = 0 (surface release)

        xmax = 0.0
        if (HEFF.gt.0) then
         call MAXP(HEFF,pmax,xmax)
        endif

        if (dist.lt.xmax) then !dist = xmax
         SIGMAZ = (0.06D0 * xmax) / SQRT(1.0D0 + 0.0015D0 * xmax)
         CAPSIG = SQRT(SIGMAZ**2 + (0.0/ PI))
         zz = 0.5D0 * ((HEFF / CAPSIG) ** 2.0D0)
         IF (zz .GT. 0.001D0 ) THEN
         pp = 2.032D0 * EXP(-zz) / (xmax * CAPSIG)
         write(*,*)'a1',dist,xmax,HEFF,CAPSIG,pp
         ELSE
         pp = 2.032D0 / (xmax * CAPSIG)
         write(*,*)'a2',dist,xmax,HEFF
         END IF
        else
         SIGMAZ = (0.06D0 * dist) / SQRT(1.0D0 + 0.0015D0 * dist)
         CAPSIG = SQRT(SIGMAZ**2 + (0.0/ PI))
         zz = 0.5D0 * ((HEFF / CAPSIG) ** 2.0D0)
         IF (zz .GT. 0.001D0 ) THEN
         pp = 2.032D0 * EXP(-zz) / (dist * CAPSIG)
              write(*,*)'b1',dist,xmax,HEFF,CAPSIG,pp
         ELSE
         pp = 2.032D0 / (dist * CAPSIG)
         write(*,*)'b2',dist,xmax,HEFF
         END IF
        endif

       ELSE !pp = B = PP(dist,0.0,AG) for other situations
        SIGMAZ = (0.06D0 * dist) / SQRT(1.0D0 + 0.0015D0 * dist)
        CAPSIG = SQRT(SIGMAZ**2 + ((xbldheight*xbldwidth) / PI))
        zz = 0.5D0 * ((0.0 / CAPSIG) ** 2.0D0)
        IF (zz .GT. 0.001D0) then
        pp = 2.032D0 * EXP(-zz) / (dist * CAPSIG)
        ELSE
        pp = 2.032D0 / (dist * CAPSIG)
        ENDIF

       ENDIF

       ! ------------------------------------------------------------
       ! Near-field / far-field concentration
       ! ------------------------------------------------------------

       IF (dist .LE. 0.0) THEN
        xconc = 0.0
       ELSE
       IF (dist .LE. xcrit) THEN
       IF (xbldwidth .LE. xbldheight) THEN
         xconc = ((xwindbinfreq / (1.0D0 - pctcalm)) * xrademis) / 
     &                 (PI * xwindbinspeed * 1.0D0 * xbldwidth)
       ELSE
         xconc = ((xwindbinfreq / (1.0D0 - pctcalm)) * xrademis) /   
     &                 (PI * xwindbinspeed * 1.0D0 * xbldheight)
       END IF !condition where building width le building height
       ELSE
       xconc = ((xwindbinfreq / (1.0D0 - pctcalm)) * xrademis * pp) / xwindbinspeed
       END IF !condition where dist lt xcrit
       END IF !condition where dist le 0.0 

       IF (xconc .LT. 1.0D-32) xconc = 1.0D-32


       ! Transit time from source to receptor (days)
       transittime = (dist / xwindbinspeed) * (1.0D0 / 3600.0D0) * (1.0D0 / 24.0D0)


       END SUBROUTINE calc_air_concentration

!=================================================================================================

       SUBROUTINE wetdep_calc(dist, pctcalm, rainfall,
     & xwindbinfreq, xwindbinspeed, xrademis, xconc, wetdepadj)

c ============================================================================
c WET DEPOSITION ADJUSTMENT FACTOR
c ============================================================================
c Account for removal of radioactive material via precipitation (wet
c deposition) combined with dry deposition velocity. This adjusts the
c ground-level deposition for rainfall conditions.
c
c Rainfall rate is in m/yr, emissions in Ci/s, concentration in Ci/m3
c ============================================================================

       IMPLICIT NONE

       REAL*8  :: dist, pctcalm, wetdepadj
       REAL*8  :: xwindbinfreq, xwindbinspeed
       REAL*8  :: xrademis, xconc, rainfall, fdux

       fdux = (xwindbinfreq/(1.0-pctcalm)) / (xwindbinspeed * dist)

     
     
!       wetdepadj = (8.0 * fdux * rainfall * 1.0E-05 * xrademis * 3600. * 24.) / (3.1416*xconc/xrademis*3600.*24.)
       wetdepadj = (8.0 * fdux * rainfall * 1.0E-05 * xrademis * 3600. * 24.) / (3.1416*xconc)

           !FDUX     = 0.25 / ( U * X )
           ! VWET = 8.0 * FDUX * RAINRATE * Cr * Q * SecPerDay / ( Pi * C )


       END SUBROUTINE wetdep_calc

!=================================================================================================

       SUBROUTINE immersion_calc(progeny_max,xhalflife,dosefactor,xconc,transittime,nchain,imdose,ambientactivity)

       IMPLICIT NONE

       INTEGER progeny_max

       REAL*8 xhalflife(progeny_max)              !radionuclide half life (years)
       REAL*8 imdose(progeny_max)                 !immersion dose
       REAL*8 sourceterm(progeny_max)             !source term for decay calculation
       REAL*8 initialamt(progeny_max)             !initial amount of nuclide
       REAL*8 activity(progeny_max)               !computed activity (Ci/m3)
       REAL*8 ambientactivity(progeny_max)        !activity in air
       REAL*8 dosefactor(progeny_max)             !dose factors 

       REAL*8 :: transittime,xconc
       REAL*8 :: xtime, xremove
       INTEGER :: nchain,l


         do l = 1, progeny_max !nchain !loop over radionuclide and progeny
         sourceterm(l) = 0.0
         initialamt(l) = 0.0
         activity(l) = 0.0
         ambientactivity(l) = 0.0
         imdose(l) = 0.0
         enddo

         initialamt(1) = xconc
         xremove = 0.0
         xtime = transittime

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

         do l = 1, nchain
         ambientactivity(l) = activity(l)
         imdose(l) = activity(l) * dosefactor(l)
         enddo

         END SUBROUTINE immersion_calc

!=======================================================================================

       SUBROUTINE inhalation_calc(breath,progeny_max,xhalflife,dosefactor,xconc,transittime,nchain,ihdose,ambientactivity)

       IMPLICIT NONE

       INTEGER progeny_max

       REAL*8 xhalflife(progeny_max)              !radionuclide half life (years)
       REAL*8 ihdose(progeny_max)                 !inhalation dose
       REAL*8 sourceterm(progeny_max)             !source term for decay calculation
       REAL*8 initialamt(progeny_max)             !initial amount of nuclide
       REAL*8 activity(progeny_max)               !computed activity (Ci/m3)
       REAL*8 ambientactivity(progeny_max)        !activity in air
       REAL*8 dosefactor(progeny_max)             !dose factors
 
       DOUBLE PRECISION breath

       REAL*8 :: transittime,xconc
       REAL*8 :: xtime, xremove
       INTEGER :: nchain,l


         do l = 1, progeny_max !nchain !loop over radionuclide and progeny
         sourceterm(l) = 0.0
         initialamt(l) = 0.0
         activity(l) = 0.0
         ambientactivity(l) = 0.0
         ihdose(l) = 0.0
         enddo

         initialamt(1) = xconc
         xremove = 0.0
         xtime = transittime

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

         do l = 1, nchain
         ambientactivity(l) = activity(l)
         ihdose(l) = activity(l) * dosefactor(l) * breath
         enddo

         END SUBROUTINE inhalation_calc

!=============================================================================================
       SUBROUTINE surface_calc(wetdepadj,progeny_max,xdrydepvel,xhalflife,xdecaychainname,
     &    dosefactor,REMCON,DEPTIM,nchain,surfdose,ambientactivity)

c ============================================================================
c SURFACE DEPOSITION DOSE CALCULATION
c ============================================================================
c Calculate doses from external radiation exposure to radionuclides
c deposited on the ground surface. Accounts for:
c   - Dry deposition (particle settling)
c   - Wet deposition (rain removal)
c   - Noble gases (no deposition): HE, NE, AR, KR, XE, RN
c ============================================================================


       IMPLICIT NONE

       INTEGER progeny_max

       REAL*8 xhalflife(progeny_max)              !radionuclide half life (years)
       REAL*8 surfdose(progeny_max)               !surface dose
       REAL*8 sourceterm(progeny_max)             !source term for decay calculation
       REAL*8 initialamt(progeny_max)             !initial amount of nuclide
       REAL*8 activity(progeny_max)               !computed activity (Ci/m3)
       REAL*8 ambientactivity(progeny_max)        !activity in air
       REAL*8 dosefactor(progeny_max)             !dose factors
       REAL*8 xdrydepvel(progeny_max)             !dry deposition velocities

       REAL*8 :: REMCON,DEPTIM
       REAL*8 :: xtime, xremove, wetdepadj
       INTEGER :: nchain,l

       CHARACTER*8 xdecaychainname(progeny_max)


         do l = 1, progeny_max !nchain !loop over radionuclide and progeny
         sourceterm(l) = ambientactivity(l)
         initialamt(l) = 0.0
         activity(l) = 0.0
         surfdose(l) = 0.0
         enddo

         xremove = REMCON / 365. !from original comply
         xtime = DEPTIM * 365.  !100 days is the dry dep time from original comply

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

         do l = 1, nchain
          if(xdecaychainname(l)(1:2).eq.'HE'.or.
     &       xdecaychainname(l)(1:2).eq.'NE'.or.
     &       xdecaychainname(l)(1:2).eq.'AR'.or.
     &       xdecaychainname(l)(1:2).eq.'KR'.or.
     &       xdecaychainname(l)(1:2).eq.'XE'.or.
     &       xdecaychainname(l)(1:2).eq.'RN') then
          surfdose(l) = activity(l) * dosefactor(l) * xdrydepvel(l)
          else
          surfdose(l) = activity(l) * dosefactor(l) * (xdrydepvel(l)+wetdepadj)
          endif
         enddo

         END SUBROUTINE surface_calc

!=======================================================================================

         SUBROUTINE ingestion_calc(wetdepadj,transittime,progeny_max,xdrydepvel,xhalflife,
     &    xdecaychainname,dosefactor,REMCON,DEPTIM,nchain,xconc,dveg,dmeat,dmilk,
     &    ybveg,ybfor,yfmilk,yfmeat,
     &    TWEATH,TEVEG,TEMILK,TEMEAT,FHVEG,FHFEED,CH2O,QMILK,QMEAT,FHMILK,FHMEAT,FCVEG,FCMILK,
     &    FCMEAT,CCARB,FRVEG,FRMILK,FRMEAT,YVEG,YMILK,YMEAT,PVEG,PMILK,PMEAT,QWMILK,QWMEAT,
     &    DELVEG,DELMILK,DELMEAT,CONVEG,CONMILK,CONMEAT)

       IMPLICIT NONE

       INTEGER progeny_max

       REAL*8 xhalflife(progeny_max)              !radionuclide half life (years)
       REAL*8 ingdose(progeny_max)                !ingestion dose
       REAL*8 sourceterm(progeny_max)             !source term for decay calculation
       REAL*8 initialamt(progeny_max)             !initial amount of nuclide
       REAL*8 activity(progeny_max)               !computed activity (Ci/m3)
       REAL*8 ambientactivity(progeny_max)        !activity in air
       REAL*8 dosefactor(progeny_max)             !dose factors
       REAL*8 xdrydepvel(progeny_max)             !dry deposition velocities
       REAL*8 normsourceactivity(progeny_max)
       REAL*8 aveglf(progeny_max)
       REAL*8 avegrt(progeny_max)
       REAL*8 amilklf(progeny_max)
       REAL*8 amilkrt(progeny_max)
       REAL*8 ameatlf(progeny_max)
       REAL*8 ameatrt(progeny_max)
       REAL*8 cfveg(progeny_max)
       REAL*8 cfmilk(progeny_max)
       REAL*8 cfmeat(progeny_max)
       REAL*8 ybveg(progeny_max)
       REAL*8 ybfor(progeny_max)
       REAL*8 yfmilk(progeny_max)
       REAL*8 yfmeat(progeny_max)
       REAL*8 dveg(progeny_max)
       REAL*8 dmilk(progeny_max)
       REAL*8 dmeat(progeny_max)

       REAL*8 :: REMCON,DEPTIM,TWEATH,TEVEG,TEMILK,TEMEAT,PMILK,PMEAT,QWMILK,QWMEAT,DELVEG
       REAL*8 :: FHVEG,FHFEED,CH2O,QMILK,QMEAT,FHMILK,FCVEG,FCMILK,FHMEAT,CONVEG,DELMEAT
       REAL*8 :: FCMEAT,CCARB,FRVEG,FRMILK,FRMEAT,YVEG,YMILK,YMEAT,PVEG,CONMILK,CONMEAT,DELMILK
       REAL*8 :: xtime, xremove, wetdepadj, transittime,xconc
       INTEGER :: nchain,l

       CHARACTER*8 xdecaychainname(progeny_max)

c ============================================================================
c INGESTION PATHWAY NORMALIZED ACTIVITY CALCULATION
c ============================================================================
c Calculate normalized surface activity (for unit initial concentration)
c to use as basis for ingestion pathway calculations. This is then
c multiplied by sector-dependent transfer factors to get pathway doses.
c ============================================================================

         do l = 1, nchain !nchain !loop over radionuclide and progeny
         sourceterm(l) = 0.0
         initialamt(l) = 0.0
         activity(l) = 0.0
         enddo

         initialamt(1) = 1.0
         xremove = 0.0
         xtime = transittime


      CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

         do l = 1, nchain
          if(xdecaychainname(l)(1:2).eq.'HE'.or.
     &       xdecaychainname(l)(1:2).eq.'NE'.or.
     &       xdecaychainname(l)(1:2).eq.'AR'.or.
     &       xdecaychainname(l)(1:2).eq.'KR'.or.
     &       xdecaychainname(l)(1:2).eq.'XE'.or.
     &       xdecaychainname(l)(1:2).eq.'RN') then
          normsourceactivity(l) = activity(l) * xdrydepvel(l)
          else
          normsourceactivity(l) = activity(l) * (xdrydepvel(l)+wetdepadj)
          endif
         enddo

c ============================================================================
c VEGETABLE INGESTION PATHWAY
c ============================================================================
c Calculate activity in vegetables via:
c   - Leaf deposition and weathering (short-term, TWEATH=0.693/12 day-1)
c   - Root uptake from soil (long-term, residence time ~100 years)
c Equations follow NCRP Commentary 3 recommendations.
c ============================================================================

         do l = 1, nchain !loop over radionuclide and progeny
         sourceterm(l) = normsourceactivity(l)
         initialamt(l) = 0.0D0
         activity(l) = 0.0D0
         enddo

         xremove = 0.693 / TWEATH
         xtime = TEVEG

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

         do l = 1, nchain
         aveglf(l) = activity(l)
         enddo

         do l = 1, nchain
         sourceterm(l) = normsourceactivity(l)
         initialamt(l) = 0.0D0
         activity(l) = 0.0D0
         enddo

         xremove = REMCON / 365.D0  ! REMCON variable convert 1/yr to 1/day
         xtime = DEPTIM * 365.D0 ! DEPTIM variable convert yr to day

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

         do l = 1, nchain
         avegrt(l) = activity(l)
         enddo

c ============================================================================
c MILK INGESTION PATHWAY (VIA FORAGE)
c ============================================================================
c Calculate activity in milk via:
c   - Forage leaf deposition and weathering (short-term)
c   - Forage root uptake from soil (long-term)
c   - Transfer to milk via cow consumption and biological transfer factor
c ============================================================================

         do l = 1, nchain !loop over radionuclide and progeny
         sourceterm(l) = normsourceactivity(l)
         initialamt(l) = 0.0
         activity(l) = 0.0
         enddo

         xremove = 0.693 / TWEATH
         xtime = TEMILK

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

         do l = 1, nchain
         amilklf(l) = activity(l)
         enddo

         do l = 1, nchain
         sourceterm(l) = normsourceactivity(l)
         initialamt(l) = 0.0
         activity(l) = 0.0
         enddo

         xremove = REMCON / 365. ! REMCON variable
         xtime = DEPTIM * 365. ! DEPTIM variable

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

         do l = 1, nchain
         amilkrt(l) = activity(l)
         enddo


c ============================================================================
c MEAT INGESTION PATHWAY (VIA FORAGE)
c ============================================================================
c Calculate activity in meat via:
c   - Forage leaf deposition and weathering (short-term)
c   - Forage root uptake from soil (long-term)
c   - Transfer to meat via animal consumption and biological transfer factor
c ============================================================================

         do l = 1, nchain !loop over radionuclide and progeny
         sourceterm(l) = normsourceactivity(l)
         initialamt(l) = 0.0
         activity(l) = 0.0
         enddo

         xremove = 0.693 / TWEATH
         xtime = TEMEAT

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

         do l = 1, nchain
         ameatlf(l) = activity(l)
         enddo

         do l = 1, nchain
         sourceterm(l) = normsourceactivity(l)
         initialamt(l) = 0.0
         activity(l) = 0.0
         enddo

         xremove = REMCON / 365.  ! REMCON variable
         xtime = DEPTIM * 365. ! DEPTIM variable

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

         do l = 1, nchain
         ameatrt(l) = activity(l)
         enddo

c ============================================================================
c FOOD TRANSFER FACTOR CALCULATIONS
c ============================================================================
c Calculate nuclide concentration in food media (vegetables, milk, meat)
c using element-specific transfer models:
c   - Tritium (H-3): water-based model (high transfer to biota)
c   - Carbon-14: carbon cycle model (high transfer to biota)
c   - Other nuclides: NCRP Commentary 3 Equation B2 model
c ============================================================================

         do l = 1, nchain

          if(xdecaychainname(l)(1:2).eq.'H-') then !H-3 (tritium) -  D. A. Baker, et al, EPA 600/9-76-016, July 1976.
          cfveg(l) = (xconc*9.0*FHVEG/CH2O)
          cfmilk(l) = (xconc*9.0*FHFEED/CH2O)*QMILK*FHMILK/(FHFEED*QMILK + QWMILK/9.0)
          cfmeat(l) = (xconc*9.0*FHFEED/CH2O)*QMEAT*FHMEAT/(FHFEED*QMEAT + QWMEAT/9.0)

          elseif(xdecaychainname(l)(1:2).eq.'C-') then !C14 - D. A. Baker, et al, EPA 600/9-76-016, July 1976.
          cfveg(l) = xconc*(FCVEG/CCARB)
          cfmilk(l) = xconc*(FCMILK/CCARB)
          cfmeat(l) = xconc*(FCMEAT/CCARB)

          else !Equation B2 of NCRP Commentary 3
          cfveg(l) =  xconc * ((FRVEG/YVEG)*aveglf(l)+(ybveg(l)/PVEG)*avegrt(l))
          cfmilk(l) = xconc * ((FRMILK/YMILK)*amilklf(l)+(ybfor(l)/PMILK)*amilkrt(l))*yfmilk(l)*QMILK
          cfmeat(l) = xconc * ((FRMEAT/YMEAT)*ameatlf(l)+(ybfor(l)/PMEAT)*ameatrt(l))*yfmeat(l)*QMEAT
           endif

        enddo


c ============================================================================
c INGESTION DOSE CALCULATIONS - VEGETABLES, MILK, MEAT
c ============================================================================
c For each food pathway:
c   - Decay concentration in food from harvest to consumption
c   - Apply ingestion dose conversion factor
c   - Multiply by consumption rate (ingestion factor)
c ============================================================================

           do l = 1, nchain !loop over radionuclide and progeny
           sourceterm(l) = 0.0
           initialamt(l) = cfveg(l)
           activity(l) = 0.0
           enddo

           xremove = 0.0
           xtime = DELVEG

           CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

            do l = 1, nchain
            dveg(l) = activity(l) * dosefactor(l) * CONVEG
            enddo



           do l = 1, nchain !loop over radionuclide and progeny
           sourceterm(l) = 0.0
           initialamt(l) = cfmilk(l)
           activity(l) = 0.0
           enddo

           xremove = 0.0
           xtime = DELMILK

           CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

            do l = 1, nchain
            dmilk(l) = activity(l) * dosefactor(l) * CONMILK
            enddo


           do l = 1, nchain !loop over radionuclide and progeny
           sourceterm(l) = 0.0
           initialamt(l) = cfmeat(l)
           activity(l) = 0.0
           enddo

           xremove = 0.0
           xtime = DELMEAT

           CALL skrable(activity,sourceterm,xtime,nchain,xremove,xhalflife,initialamt)

            do l = 1, nchain
            dmeat(l) = activity(l) * dosefactor(l) * CONMEAT
            enddo


       END SUBROUTINE ingestion_calc

!=========================================================================================

       SUBROUTINE MAXP( H, PMAX, XMAX )

!**********************************************************************
!*                                                                    *
!*  SUBROUTINE MAXP FINDS THE DISTANCE XMAX AT WHICH THE CONCENTRA-   *
!*  TION REACHES ITS MAXIMUM VALUE                                    *
!*                                                                    *
!**********************************************************************

      REAL*8     H
      REAL*8     PMAX
      REAL*8     XMAX
C
      REAL*8     DELX
      REAL*8    ESTDIS
      REAL*8     HH( 129 )
      INTEGER  I
      REAL*8     P1
      REAL*8     P2
      REAL*8     P3
      REAL*8     XLOW
      REAL*8     XUP
      REAL*8     XX( 129 )


!                            H                X @ PMAX FOR H
!                     ----------             ----------------
       DATA HH(   1 ) /      1.0/, XX(   1 ) /      8.76750374/
       DATA HH(   2 ) /      2.0/, XX(   2 ) /     18.52621080/
       DATA HH(   3 ) /      3.0/, XX(   3 ) /     30.23082350/
       DATA HH(   4 ) /      4.0/, XX(   4 ) /     42.14927290/
       DATA HH(   5 ) /      5.0/, XX(   5 ) /     53.96775820/
       DATA HH(   6 ) /      6.0/, XX(   6 ) /     69.07896420/
       DATA HH(   7 ) /      7.0/, XX(   7 ) /     80.74008180/
       DATA HH(   8 ) /      8.0/, XX(   8 ) /     96.74093630/
       DATA HH(   9 ) /      9.0/, XX(   9 ) /    113.46733100/
       DATA HH(  10 ) /     10.0/, XX(  10 ) /    124.42367600/
       DATA HH(  11 ) /     11.0/, XX(  11 ) /    141.56265300/
       DATA HH(  12 ) /     12.0/, XX(  12 ) /    151.58261100/
       DATA HH(  13 ) /     13.0/, XX(  13 ) /    160.57955900/
       DATA HH(  14 ) /     14.0/, XX(  14 ) /    177.52858000/
       DATA HH(  15 ) /     15.0/, XX(  15 ) /    194.91195700/
       DATA HH(  16 ) /     16.0/, XX(  16 ) /    212.71072400/
       DATA HH(  17 ) /     17.0/, XX(  17 ) /    230.90789800/
       DATA HH(  18 ) /     18.0/, XX(  18 ) /    237.51951600/
       DATA HH(  19 ) /     19.0/, XX(  19 ) /    242.86346400/
       DATA HH(  20 ) /     20.0/, XX(  20 ) /    260.33010900/
       DATA HH(  21 ) /     21.0/, XX(  21 ) /    278.10879500/
       DATA HH(  22 ) /     22.0/, XX(  22 ) /    296.18981900/
       DATA HH(  23 ) /     23.0/, XX(  23 ) /    314.56417800/
       DATA HH(  24 ) /     24.0/, XX(  24 ) /    333.22366300/
       DATA HH(  25 ) /     25.0/, XX(  25 ) /    352.16043100/
       DATA HH(  26 ) /     26.0/, XX(  26 ) /    371.36731000/
       DATA HH(  27 ) /     27.0/, XX(  27 ) /    390.83758500/
       DATA HH(  28 ) /     28.0/, XX(  28 ) /    410.56481900/
       DATA HH(  29 ) /     29.0/, XX(  29 ) /    409.99679600/
       DATA HH(  30 ) /     30.0/, XX(  30 ) /    407.83245800/
       DATA HH(  31 ) /     31.0/, XX(  31 ) /    426.34710700/
       DATA HH(  32 ) /     32.0/, XX(  32 ) /    445.07446300/
       DATA HH(  33 ) /     33.0/, XX(  33 ) /    464.01013200/
       DATA HH(  34 ) /     34.0/, XX(  34 ) /    483.15008500/
       DATA HH(  35 ) /     35.0/, XX(  35 ) /    502.49038700/
       DATA HH(  36 ) /     36.0/, XX(  36 ) /    522.02722200/
       DATA HH(  37 ) /     37.0/, XX(  37 ) /    541.75720200/
       DATA HH(  38 ) /     38.0/, XX(  38 ) /    561.67688000/
       DATA HH(  39 ) /     39.0/, XX(  39 ) /    581.78308100/
       DATA HH(  40 ) /     40.0/, XX(  40 ) /    602.07257100/
       DATA HH(  41 ) /     41.0/, XX(  41 ) /    622.54241900/
       DATA HH(  42 ) /     42.0/, XX(  42 ) /    643.18981900/
       DATA HH(  43 ) /     43.0/, XX(  43 ) /    664.01196300/
       DATA HH(  44 ) /     44.0/, XX(  44 ) /    685.00622600/
       DATA HH(  45 ) /     45.0/, XX(  45 ) /    706.17010500/
       DATA HH(  46 ) /     46.0/, XX(  46 ) /    727.50109900/
       DATA HH(  47 ) /     47.0/, XX(  47 ) /    748.99700900/
       DATA HH(  48 ) /     48.0/, XX(  48 ) /    770.65533400/
       DATA HH(  49 ) /     49.0/, XX(  49 ) /    792.47406000/
       DATA HH(  50 ) /     50.0/, XX(  50 ) /    814.45105000/
       DATA HH(  51 ) /     51.0/, XX(  51 ) /    836.58410600/
       DATA HH(  52 ) /     52.0/, XX(  52 ) /    858.87133800/
       DATA HH(  53 ) /     53.0/, XX(  53 ) /    881.31079100/
       DATA HH(  54 ) /     54.0/, XX(  54 ) /    903.90075700/
       DATA HH(  55 ) /     55.0/, XX(  55 ) /    926.63916000/
       DATA HH(  56 ) /     56.0/, XX(  56 ) /    949.52453600/
       DATA HH(  57 ) /     57.0/, XX(  57 ) /    972.55493200/
       DATA HH(  58 ) /     58.0/, XX(  58 ) /    995.72876000/
       DATA HH(  59 ) /     59.0/, XX(  59 ) /   1019.04468000/
       DATA HH(  60 ) /     60.0/, XX(  60 ) /   1042.50073000/
       DATA HH(  61 ) /     61.0/, XX(  61 ) /   1066.09570000/
       DATA HH(  62 ) /     62.0/, XX(  62 ) /   1089.82788000/
       DATA HH(  63 ) /     63.0/, XX(  63 ) /   1113.69604000/
       DATA HH(  64 ) /     64.0/, XX(  64 ) /   1137.69873000/
       DATA HH(  65 ) /     65.0/, XX(  65 ) /   1161.83447000/
       DATA HH(  66 ) /     66.0/, XX(  66 ) /   1186.10205000/
       DATA HH(  67 ) /     67.0/, XX(  67 ) /   1210.50024000/
       DATA HH(  68 ) /     68.0/, XX(  68 ) /   1235.02759000/
       DATA HH(  69 ) /     69.0/, XX(  69 ) /   1259.68286000/
       DATA HH(  70 ) /     70.0/, XX(  70 ) /   1284.46521000/
       DATA HH(  71 ) /     71.0/, XX(  71 ) /   1309.37305000/
       DATA HH(  72 ) /     72.0/, XX(  72 ) /   1334.40540000/
       DATA HH(  73 ) /     73.0/, XX(  73 ) /   1359.56116000/
       DATA HH(  74 ) /     74.0/, XX(  74 ) /   1384.83911000/
       DATA HH(  75 ) /     75.0/, XX(  75 ) /   1410.23840000/
       DATA HH(  76 ) /     76.0/, XX(  76 ) /   1435.75781000/
       DATA HH(  77 ) /     77.0/, XX(  77 ) /   1461.39636000/
       DATA HH(  78 ) /     78.0/, XX(  78 ) /   1487.15308000/
       DATA HH(  79 ) /     79.0/, XX(  79 ) /   1513.02686000/
       DATA HH(  80 ) /     80.0/, XX(  80 ) /   1539.01709000/
       DATA HH(  81 ) /     81.0/, XX(  81 ) /   1565.12231000/
       DATA HH(  82 ) /     82.0/, XX(  82 ) /   1591.34204000/
       DATA HH(  83 ) /     83.0/, XX(  83 ) /   1617.67529000/
       DATA HH(  84 ) /     84.0/, XX(  84 ) /   1644.12085000/
       DATA HH(  85 ) /     85.0/, XX(  85 ) /   1670.67822000/
       DATA HH(  86 ) /     86.0/, XX(  86 ) /   1697.34644000/
       DATA HH(  87 ) /     87.0/, XX(  87 ) /   1724.12463000/
       DATA HH(  88 ) /     88.0/, XX(  88 ) /   1751.01196000/
       DATA HH(  89 ) /     89.0/, XX(  89 ) /   1778.00781000/
       DATA HH(  90 ) /     90.0/, XX(  90 ) /   1805.11108000/
       DATA HH(  91 ) /     91.0/, XX(  91 ) /   1832.32129000/
       DATA HH(  92 ) /     92.0/, XX(  92 ) /   1859.63745000/
       DATA HH(  93 ) /     93.0/, XX(  93 ) /   1887.05908000/
       DATA HH(  94 ) /     94.0/, XX(  94 ) /   1914.58521000/
       DATA HH(  95 ) /     95.0/, XX(  95 ) /   1942.21521000/
       DATA HH(  96 ) /     96.0/, XX(  96 ) /   1969.94849000/
       DATA HH(  97 ) /     97.0/, XX(  97 ) /   1997.78394000/
       DATA HH(  98 ) /     98.0/, XX(  98 ) /   2025.72119000/
       DATA HH(  99 ) /     99.0/, XX(  99 ) /   2053.75977000/
       DATA HH( 100 ) /    100.0/, XX( 100 ) /   2081.89868000/
       DATA HH( 101 ) /    110.0/, XX( 101 ) /   2368.67383000/
       DATA HH( 102 ) /    120.0/, XX( 102 ) /   2664.83911000/
       DATA HH( 103 ) /    130.0/, XX( 103 ) /   3282.53271000/
       DATA HH( 104 ) /    140.0/, XX( 104 ) /   3629.00122000/
       DATA HH( 105 ) /    150.0/, XX( 105 ) /   3984.34863000/
       DATA HH( 106 ) /    160.0/, XX( 106 ) /   4348.18701000/
       DATA HH( 107 ) /    170.0/, XX( 107 ) /   4720.16992000/
       DATA HH( 108 ) /    180.0/, XX( 108 ) /   5366.60059000/
       DATA HH( 109 ) /    190.0/, XX( 109 ) /   6065.17432000/
       DATA HH( 110 ) /    200.0/, XX( 110 ) /   6501.37939000/
       DATA HH( 111 ) /    220.0/, XX( 111 ) /   7396.92480000/
       DATA HH( 112 ) /    240.0/, XX( 112 ) /   9198.33008000/
       DATA HH( 113 ) /    260.0/, XX( 113 ) /  10251.25200000/
       DATA HH( 114 ) /    280.0/, XX( 114 ) /  11920.05860000/
       DATA HH( 115 ) /    300.0/, XX( 115 ) /  13087.25590000/
       DATA HH( 116 ) /    350.0/, XX( 116 ) /  16997.97270000/
       DATA HH( 117 ) /    400.0/, XX( 117 ) /  21517.41020000/
       DATA HH( 118 ) /    450.0/, XX( 118 ) /  25237.74220000/
       DATA HH( 119 ) /    500.0/, XX( 119 ) /  30820.03710000/
       DATA HH( 120 ) /    550.0/, XX( 120 ) /  37208.56640000/
       DATA HH( 121 ) /    600.0/, XX( 121 ) /  41860.91410000/
       DATA HH( 122 ) /    650.0/, XX( 122 ) /  49608.53520000/
       DATA HH( 123 ) /    700.0/, XX( 123 ) /  54844.67190000/
       DATA HH( 124 ) /    750.0/, XX( 124 ) /  64161.60160000/
       DATA HH( 125 ) /    800.0/, XX( 125 ) /  70020.64060000/
       DATA HH( 127 ) /    900.0/, XX( 127 ) /  87683.97660000/
       DATA HH( 128 ) /    950.0/, XX( 128 ) /  94343.85160000/
       DATA HH( 129 ) /   1000.0/, XX( 129 ) / 108178.82800000/

        DO I = 1, 128

           IF ( HH( I ) .LE. H  .AND.  H .LE. HH( I + 1 ) ) THEN

              XLOW = XX( I )
              XUP  = XX( I + 1 )

              GOTO 2

           ENDIF

        ENDDO

   2     CONTINUE


!     Having extablished lower and upper bounds, now find point
!     of maximum P which must lie between them.
 
  20  CONTINUE

      DELX = ( XLOW - XUP  ) / 2.0

!     Find P at lower, mid-point, and upper values of x.

      P1 = PP( XUP ,       H, 0.0 )
      P2 = PP( XUP + DELX, H, 0.0 )
      P3 = PP( XLOW,       H, 0.0 )

       IF ( ( P1 .LT. 1.0E-40 ) .OR.
     &     ( P2 .LT. 1.0E-40 ) .OR.
     &     ( P3 .LT. 1.0E-40 ) ) THEN

         PMAX = PP( XUP + DELX, H, 0.0 )

         XMAX = XUP + DELX

         GO TO 30

       ENDIF

!     If ratios of P at middle x to P's at upper and lower x's
!     are within 1 %, stop iteration.   Maximum P will then be
!     known to within 1 %.  XMAX is taken as halfway between
!     XUP and XLOW.

       IF ( P2/P1 .LE. 1.01 .AND. P2/P3 .LE. 1.01 ) THEN

         PMAX = PP( XUP + DELX, H, 0.0 )
         XMAX = XUP + DELX

         GO TO 30

       ENDIF

!     If test was failed, increase XUP  and decrease XLOW.

       XUP  = XUP  + DELX

       XLOW = XLOW - DELX
       GO TO 20

  30   CONTINUE


       RETURN
       END SUBROUTINE MAXP

!========================================================================

      FUNCTION PP( X, H, AG )

      REAL*8  SIGMAZ, CAPSIG, ZZ
      REAL*8  X, H, AG
C
C     Equation A-5 of NCRP
C

      SIGMAZ = 0.06 * X / SQRT( 1.0 + 0.0015 * X )

C
C     Equation A-11 of NCRP. Note that if AG=0, CAPSIG=SIGMAZ.
C

      CAPSIG = SQRT( SIGMAZ ** 2 + AG / 3.14 )

C
C     Equation A-7 with SIGMAZ replaced by CAPSIG. When AG=0.
C     this is equivalent to Equation A-7 of NCRP.
C
      ZZ = 0.5 * ( ( H / CAPSIG ) ** 2 )

      IF ( ZZ .GT. 0.001 ) THEN

         PP = 2.032 * EXP( -ZZ ) / ( X * CAPSIG )

      ELSE

         PP = 2.032 / ( X * CAPSIG )

      ENDIF

C     write(*,'(a)')' x,h,capsig,zz,pp'
C     write( *,'( 1x,5(1pe9.2,1x) )' )x,h,capsig,zz,pp

      RETURN
      END

