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
      REAL*8, ALLOCATABLE :: distance(:)            !distance between stack and receptor (m)
      REAL*8, ALLOCATABLE :: distanceveg(:)         !distance between garden and receptor (m)
      REAL*8, ALLOCATABLE :: distancecow(:)         !distance between cows/milk and receptor (m)
      REAL*8, ALLOCATABLE :: distancemeat(:)        !distance between meat and receptor (m)
      REAL*8, ALLOCATABLE :: latitude(:)            !latitude (decimal degrees)
      REAL*8, ALLOCATABLE :: longitude(:)           !longitude (decimal degrees)

c ============================================================================
c SECTION 2: DYNAMIC ARRAY DECLARATIONS - EMISSION & CONCENTRATION DATA
c ============================================================================
c These arrays store radionuclide-specific emission rates and computed
c air/dose concentrations for each stack/radionuclide combination

      REAL*8, ALLOCATABLE :: rademis(:,:)           !radionuclide emissions (Ci/yr) for each stack
      REAL*8, ALLOCATABLE :: conc(:,:)              !radionuclide air concentration (Bq/m3)
      REAL*8, ALLOCATABLE :: dose(:,:)              !radionuclide dose (mrem/yr)

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
      REAL*8, ALLOCATABLE :: xbveg(:)               !root uptake factor for vegetables
      REAL*8, ALLOCATABLE :: xbfor(:)               !root uptake factor for forage
      REAL*8, ALLOCATABLE :: xfmilk(:)              !transfer factor to milk
      REAL*8, ALLOCATABLE :: xfmeat(:)              !transfer factor to meat

c ============================================================================
c SECTION 6: DYNAMIC ARRAY DECLARATIONS - RADIONUCLIDE DOSE FACTORS
c ============================================================================
c These arrays store dose conversion factors (DCF) from NCRP/ICRP references
c that convert nuclide activity to effective dose equivalent

      REAL*8, ALLOCATABLE :: halflife(:)            !radionuclide half life (years)
      REAL*8, ALLOCATABLE :: wbfac(:)               !whole body dose factor
      REAL*8, ALLOCATABLE :: xwbfac(:)              !whole body dose factor
      REAL*8, ALLOCATABLE :: vegfac(:)              !vegetable dose factor
      REAL*8, ALLOCATABLE :: cowfac(:)              !cow dose factor
      REAL*8, ALLOCATABLE :: inhalationfac(:)       !inhalation dose factor
      REAL*8, ALLOCATABLE :: immersionfac(:)        !immersion dose factor
      REAL*8, ALLOCATABLE :: surfacefac(:)          !surface dose factor
      REAL*8, ALLOCATABLE :: ingestionfac(:)        !ingestion dose factor
      REAL*8, ALLOCATABLE :: drydepvel(:)           !dry deposition velocity (m/day)

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

      INTEGER, ALLOCATABLE :: nrads(:)                !number of emitted radionucludes at each stack

      CHARACTER*8,  ALLOCATABLE :: radname(:,:)       !radionuclide names at each stack
      CHARACTER*3,  ALLOCATABLE :: windbinname(:)     !wind direction sector 
      CHARACTER*3,  ALLOCATABLE :: sectorname(:)      !sector name for receptors
      CHARACTER*3,  ALLOCATABLE :: vegsector(:)       !Sector (wind dir) location of vegetables
      CHARACTER*3,  ALLOCATABLE :: milksector(:)      !Sector (wind dir) location of cow/milk
      CHARACTER*3,  ALLOCATABLE :: meatsector(:)      !Sector (wind dir) location of meat
      CHARACTER*8,  ALLOCATABLE :: decaychain(:)      !Decay chain list
      CHARACTER*8,  ALLOCATABLE :: decaychainname(:)  !Decay chain radionuclide names
      CHARACTER*3,  ALLOCATABLE :: verticalmomentum(:)!vertical momentum (YES or NO)

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

      REAL*8 pctcalm               !percentage of calm days at wind data location
      REAL*8 rainfall              !rainfall rate (m/yr)
      REAL*8 airtemp               !air temperature (C)


      REAL*8 FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF,FACTORG,FACTORH,FACTORI
      REAL*8 FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN
      REAL*8 area,xcrit,SIGMAZ,CAPSIG,pp,alldose,idose,dist,altconc,maxconc
      REAL*8 tempdose,maxdose,maxdist,totmilkdose,totvegdose,itotvegdose,totmeatdose,itotcowdose
      REAL*8 capp,capn0,outprod,inprod,z1,term1,term2,sum
      REAL*8 transittime                              !transit time of air from emission release to receptor (s)
      REAL*8 wetdepadj                                !wet deposition adjustment factor
      REAL*8 xtime, xremove
      REAL*8 ingestion, immersion, inhalation, surface
      REAL*8 RHOAIR,RHOSTK,STKVEL,R,FM,DELHMO,DELHBU,HEFF,FB,zz,INHFAC,TINIT,T

c STATIC VARIABLE DECLARATIONS - CONTROL & LOOP INDICES

      INTEGER nstacks                               !number of unique stacks
      INTEGER receptor, N, M, I, J, K, L, II, JJ, P, Q                        !loop counters
      INTEGER nchain                                !number of nuclides in the decay series
      INTEGER maxsect
      INTEGER nrads_max

      CHARACTER*200 inline
      CHARACTER*8 NUC
      CHARACTER*2 NN

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


c ============================================================================
c SECTION 13: MEMORY ALLOCATION - DECAY CHAIN (N=19) RELATED
c ============================================================================

       allocate(halflife(19))
       allocate(immersionfac(19))
       allocate(inhalationfac(19))
       allocate(ingestionfac(19))
       allocate(surfacefac(19))
       allocate(xbveg(19))
       allocate(xbfor(19))
       allocate(xfmilk(19))
       allocate(xfmeat(19))
       allocate(wbfac(19))
       allocate(xwbfac(19))
       allocate(vegfac(19))
       allocate(cowfac(19))
       allocate(decaychain(19))
       allocate(drydepvel(19))
       allocate(decaychainname(19))
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


c Write to diagnostic output file
       write(12,4000)
       write(12,2001)
       write(12,1000) title
       write(12,1001) company
       write(12,1002) facility
       write(12,1003) address
       write(12,1004) poc
       write(12,1005) phone
       write(12,1013) airtemp
       write(12,1018) rainfall

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
       allocate(distance(nstacks))
       allocate(distanceveg(nstacks))
       allocate(distancecow(nstacks))
       allocate(distancemeat(nstacks))
       allocate(vegsector(nstacks))
       allocate(milksector(nstacks))
       allocate(meatsector(nstacks))
       allocate(nrads(nstacks))

       nrads_max = 3000 !maximum number of radionuclides possible per stack

       allocate(radname(nstacks,nrads_max))
       allocate(rademis(nstacks,nrads_max))
       allocate(conc(nstacks,nrads_max))
       allocate(dose(nstacks,nrads_max))
       allocate(totimdose(nstacks,nrads_max))
       allocate(totihdose(nstacks,nrads_max))
       allocate(totsurfdose(nstacks,nrads_max))
       allocate(totingdose(nstacks,nrads_max))
       allocate(milkdose(nstacks,nrads_max))
       allocate(vegdose(nstacks,nrads_max))
       allocate(meatdose(nstacks,nrads_max))

c ============================================================================
c SECTION 17: INPUT PARSING LOOP - READ STACK-SPECIFIC PARAMETERS
c ============================================================================
c Main loop to read input parameters for each stack source. Includes
c geometrical parameters, temperatures, distances to receptors, and
c number of emitted radionuclides

       do n = 1, nstacks 

       read(11,'(a)') !blank spacer to improve readability of the input file
       read(11,'(54x,f10.5)') latitude(n)
       read(11,'(54x,f10.5)') longitude(n)
       read(11,'(54x,f10.5)') releaseht(n)
       read(11,'(54x,f10.5)') bldheight(n)
       read(11,'(54x,f10.5)') bldwidth(n)
       read(11,'(54x,f10.5)') flowrate(n)
       read(11,'(54x,f10.5)') diameter(n)
       read(11,'(54x,f10.5)') stacktemp(n)
       read(11,'(54x,a)') verticalmomentum(n)
       read(11,'(54x,f10.5)') distance(n)
       read(11,'(54x,f10.5)') distanceveg(n)
       read(11,'(54x,f10.5)') distancecow(n)
       read(11,'(54x,f10.5)') distancemeat(n)
       read(11,'(54x,a)') vegsector(n)
       read(11,'(54x,a)') milksector(n)
       read(11,'(54x,a)') meatsector(n)

       read(11,'(54x,i5)') nrads(n)

       stacktemp(n) = stacktemp(n) * (9.D0/5.D0) + 32.D0 !convert C to F


c write to diagnostic output file

       write(12,2001)
       write(12,1006) n
       write(12,1007) latitude(n)
       write(12,1008) longitude(n)
       write(12,1009) releaseht(n)
       write(12,1010) diameter(n)
       write(12,1011) flowrate(n)
       write(12,1012) stacktemp(n)
       write(12,1014) verticalmomentum(n)
       write(12,1015) bldheight(n)
       write(12,1016) bldwidth(n)
       write(12,1017) distance(n)
       write(12,1018) rainfall
       write(12,1019) nrads(n)

c ============================================================================
c SECTION 19: INPUT PARSING - READ RADIONUCLIDE EMISSION DATA
c ============================================================================
c For each stack, read the list of radionuclides emitted and their
c emission rates (in Ci/s)

       do m = 1, nrads(n)
       read(11,'(54x,a)') inline
       read(inline,*) radname(n,m),rademis(n,m)
       rademis(n,m) = rademis(n,m) * (1./365.)*(1./24.)*(1./60.)*(1./60.) !year to seconds
       enddo !end loop over reading in radionuclide emissions


       enddo !end loop over reading in stacks

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

       write(13,3001) !write header to file

       do n = 1, nstacks

       alldose = 0.

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

       maxconc = 0.0
       maxdose = 0.0
       inhalation = 0.
       immersion = 0.
       ingestion = 0.
       surface = 0.
       totmilkdose = 0.
       totvegdose = 0.
       totmeatdose = 0.

         do l = 1, 19
         halflife(l) = 1.
         inhalationfac(l) = 0.
         immersionfac(l) = 0.
         surfacefac(l) = 0.
         ingestionfac(l) = 0.
         xbveg(l) = 0.
         xbfor(l) = 0.
         xfmilk(l) = 0.
         xfmeat(l) = 0.
         drydepvel(l) = 0.
         enddo


c ============================================================================
c SECTION 24: WIND SECTOR AND DISTANCE ITERATION LOOPS
c ============================================================================
c Loop over wind direction sectors and receptor distances. For each combination,
c calculate air concentration and dose contributions across all decay pathways

       do j = 1, nwindbins

       do receptor = 20, 6000, 20

       dist = float(receptor) 

c ============================================================================
c SECTION 25: STACK PLUME RISE CALCULATIONS
c ============================================================================
c Calculate the effective release height (HEFF) by accounting for:
c   - Buoyancy rise: hot exhaust plume rises due to thermal buoyancy
c   - Momentum rise: vertical velocity of exhaust provides additional rise
c   - Wind shear: ambient wind speed reduces plume rise
c
c Only applied if release height > 2.5 * building height (elevated release)
c ============================================================================

       if(flowrate(n).le.0.) then
       flowrate(n) = 0.3
       endif

       HEFF = releaseht(n)

       IF ( releaseht(n) .GT. 2.5*bldheight(n) ) THEN

c --- Temperature and density calculations ---
        RHOAIR = 0.0735 * ( 80.0 + 460.0 ) / ( airtemp + 460.0 )
        RHOAIR = RHOAIR * 454.0 / ( 0.3048 ** 3 )
        RHOSTK = RHOAIR * ( airtemp + 460.0 ) / ( stacktemp(n) + 460.0 )

c --- Momentum rise calculations ---
        AREA   = PI * ( diameter(n) ** 2.0 ) / 4.0
        STKVEL = flowrate(n) / AREA
        R      = STKVEL / windbinspeed(j)
        FM     = ( RHOSTK / RHOAIR ) * STKVEL * flowrate(n) / PI
        TERM1  = 75.0 * FM / ( ( ( 1.0 + 3.0 / R ) * windbinspeed(j) ) ** 2.0 )
        TERM2  = releaseht(n) + 7.79 * SQRT( FM ) / ( ( 1.0 + 3.0 / R ) * windbinspeed(j) )

c --- Apply momentum rise if vertical momentum enabled ---
        IF ( verticalmomentum(n).eq.'YES' ) THEN
         DELHMO = 0.93 * ( TERM1 ** ( 3.0 / 7.0 ) ) *
     &        ( TERM2 ** ( 1.0 / 7.0 ) )
        ELSE
         DELHMO = 0.0
        ENDIF

c --- Buoyancy rise calculations ---
        FB = 9.81 * ( RHOAIR - RHOSTK ) * flowrate(n) / ( PI * RHOAIR )
        DELHBU = 23.67 * ( ( FB / ( windbinspeed(j) ** 3.0 ) ) ** 0.6 ) *
     &           ( releaseht(n) + 195.1 * FB / ( windbinspeed(j) ** 3.0 ) ) ** 0.4

c --- Use maximum of buoyancy or momentum rise ---
       if(DELHMO.gt.DELHBU) then
       HEFF = releaseht(n) + DELHMO
       else
       HEFF = releaseht(n) + DELHBU 
       endif

       ENDIF !end condition where H > 0

c ============================================================================
c SECTION 26: GAUSSIAN PLUME MODEL - CONCENTRATION AT RECEPTOR
c ============================================================================
c Calculate air concentration at receptor distance using Gaussian plume
c model with lateral (sigma_z) and source area (CAPSIG) corrections.
c Includes near-field correction for distances less than xcrit.
c ============================================================================

       SIGMAZ = (0.06 * dist) / SQRT( 1.0 + 0.0015 * dist)
       CAPSIG = SQRT( SIGMAZ ** 2 + (bldwidth(n) * bldheight(n)) / PI )
       zz = 0.5 * ( ( HEFF / CAPSIG ) ** 2 )
        if(zz.gt.0.001.and.releaseht(n).gt.2.5*bldheight(n)) then
        pp = 2.032 * EXP( -zz ) / ( dist * CAPSIG )
        else
        pp = 2.032 / ( dist * CAPSIG )
        endif

       transittime = (dist / windbinspeed(j))*(1./3600.)*(1./24.) ! seconds to days
 
c --- Near-field vs far-field concentration calculation ---
       if(dist.lt.xcrit) then

        if(bldwidth(n).le.bldheight(n)) then        
       conc(n,m) = ( (windbinfreq(j)/(1.0-pctcalm)) *rademis(n,m)) /
     &          (PI*windbinspeed(j)*1.0*bldwidth(n))
        else
       conc(n,m) = ( (windbinfreq(j)/(1.0-pctcalm)) *rademis(n,m)) /
     &          (PI*windbinspeed(j)*1.0*bldheight(n))                 
        endif

       else
       conc(n,m) = ((windbinfreq(j)/(1.0-pctcalm)) *rademis(n,m)*pp) / windbinspeed(j)
       endif

       if(conc(n,m).lt. 1.0E-32 ) then
         conc(n,m) = 1.0E-32
       endif

       altconc = (windbinfreq(j)*rademis(n,m)*pp) / windbinspeed(j) !just for QA and not being used for anything else here

c ============================================================================
c SECTION 27: WET DEPOSITION ADJUSTMENT FACTOR
c ============================================================================
c Account for removal of radioactive material via precipitation (wet
c deposition) combined with dry deposition velocity. This adjusts the
c ground-level deposition for rainfall conditions.
c
c Rainfall rate is in m/yr, emissions in Ci/s, concentration in Ci/m3
c ============================================================================

       wetdepadj = (16.0/(2.0*3.416))*(1/dist)*((windbinfreq(j)/(1.0-pctcalm)) / (windbinspeed(j)))
     &    * rainfall * 1.0E-05 * (rademis(n,m) * 3600. * 24.) * (1./(conc(n,m)))
!       wetdepadj = (8.0 * fdux * rainfall(n) * 1.0E-05 * (rademis(n,m) * 3600. * 24.))/(3.1416*conc(n,m)/rademis(n,m) )

c ============================================================================
c SECTION 28: DECAY CHAIN LOOKUP - READ NUCLIDE DATA FROM REFERENCE TABLE
c ============================================================================
c For each emitted radionuclide, search the reference CSV file to find:
c   - Decay chain (parent + progeny nuclides)
c   - Half-lives for all nuclides in chain
c   - Dose conversion factors for all pathways
c   - Environmental transfer factors (root uptake, animal transfer, etc.)
c ============================================================================

       rewind(9)
       rewind(10)
       read(9,'(a)') inline


         do ii = 1, 879
         read(10,'(a)') inline
         read(inline,*) NUC,k,(decaychain(l),l=1,k)
         if(NUC.eq.radname(n,m)) then
          read(inline,*) NUC,k,(decaychain(l),l=1,k)

          do jj = 1, 826
          read(9,'(a)') inline
          read(inline,*) NUC,FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF,FACTORG,
     &      FACTORH,FACTORI,FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN

           do l = 1, k
           if(NUC.eq.decaychain(l)) then
            decaychainname(l) = NUC
            halflife(l) = FACTOR*365. !! convert years to days
            wbfac(l) = FACTORB
            xwbfac(l) = FACTORC
            vegfac(l) = FACTORD
            cowfac(l) = FACTORE
            inhalationfac(l) = FACTORF
            immersionfac(l) = FACTORG
            surfacefac(l) = FACTORH
            ingestionfac(l) = FACTORI
            xbveg(l) = FACTORJ
            xbfor(l) = FACTORK
            xfmilk(l) = FACTORL
            xfmeat(l) = FACTORM
            drydepvel(l) = FACTORN
           endif
           enddo !loop over decay chain
          enddo !loop over all nuclides on input reference table 
           nchain = k
         endif 
         enddo !loop over all nuclides on decay chain input table

!         do l = 1, nchain
!         write(*,*)k,l,decaychainname(l),halflife(l),inhalationfac(l),immersionfac(l),
!     &    surfacefac(l),ingestionfac(l),xbveg(l), xbfor(l),xfmilk(l),xfmeat(l),drydepvel(l) 
!         enddo

c ============================================================================
c SECTION 29: INHALATION & IMMERSION DOSE CALCULATION
c ============================================================================
c Calculate doses from inhalation of contaminated air and immersion in
c contaminated air. These are pathway doses based on air concentration
c and are independent of ground deposition.
c ============================================================================

         do l = 1, 16 !nchain !loop over radionuclide and progeny
         sourceterm(l) = 0.0
         initialamt(l) = 0.0
         activity(l) = 0.0
         ambientactivity(l) = 0.0 
         imdose(l) = 0.0
         ihdose(l) = 0.0
         enddo

         initialamt(1) = conc(n,m)
         xremove = 0.0
         xtime = transittime

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

 !Calculate immersion and inhalation dose
         totimdose(n,m) = 0.0
         totihdose(n,m) = 0.0

         do l = 1, nchain
         ambientactivity(l) = activity(l)
         imdose(l) = ambientactivity(l) * immersionfac(l)
         ihdose(l) = ambientactivity(l) * inhalationfac(l) * BREATH
         totimdose(n,m) = totimdose(n,m) + imdose(l)
         totihdose(n,m) = totihdose(n,m) + ihdose(l)
         enddo


c ============================================================================
c SECTION 30: SURFACE DEPOSITION DOSE CALCULATION
c ============================================================================
c Calculate doses from external radiation exposure to radionuclides
c deposited on the ground surface. Accounts for:
c   - Dry deposition (particle settling)
c   - Wet deposition (rain removal)
c   - Noble gases (no deposition): HE, NE, AR, KR, XE, RN
c ============================================================================

         do l = 1, nchain !loop over radionuclide and progeny
         sourceterm(l) = ambientactivity(l)
         initialamt(l) = 0.0
         activity(l) = 0.0
         enddo

         xremove = REMCON / 365. !from original comply
         xtime = DEPTIM * 365.  !100 days is the dry dep time from original comply

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

 !Calculate surface dose

         totsurfdose(n,m) = 0.0

         do l = 1, nchain
          if(decaychainname(l)(1:2).eq.'HE'.or.
     &       decaychainname(l)(1:2).eq.'NE'.or.
     &       decaychainname(l)(1:2).eq.'AR'.or.
     &       decaychainname(l)(1:2).eq.'KR'.or.
     &       decaychainname(l)(1:2).eq.'XE'.or.
     &       decaychainname(l)(1:2).eq.'RN') then
          surfaceactivity(l) = activity(l) * (drydepvel(l))
          else 
          surfaceactivity(l) = activity(l) * (drydepvel(l)+wetdepadj)
          endif
         surfdose(l) = surfaceactivity(l) * surfacefac(l)
         totsurfdose(n,m) = totsurfdose(n,m) + surfdose(l)
         enddo

c ============================================================================
c SECTION 31: INGESTION PATHWAY NORMALIZED ACTIVITY CALCULATION
c ============================================================================
c Calculate normalized surface activity (for unit initial concentration)
c to use as basis for ingestion pathway calculations. This is then
c multiplied by sector-dependent transfer factors to get pathway doses.
c ============================================================================

         do l = 1, nchain !loop over radionuclide and progeny
         sourceterm(l) = 0.0
         initialamt(l) = 0.0
         activity(l) = 0.0
         enddo

         initialamt(1) = 1.0
         xremove = 0.0
         xtime = transittime

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

         do l = 1, nchain

          if(decaychainname(l)(1:2).eq.'HE'.or.
     &       decaychainname(l)(1:2).eq.'NE'.or.
     &       decaychainname(l)(1:2).eq.'AR'.or.
     &       decaychainname(l)(1:2).eq.'KR'.or.
     &       decaychainname(l)(1:2).eq.'XE'.or.
     &       decaychainname(l)(1:2).eq.'RN') then
          normsourceactivity(l) = activity(l) * (drydepvel(l))
          else
          normsourceactivity(l) = activity(l) * (drydepvel(l)+wetdepadj)
          endif

         enddo

c ============================================================================
c SECTION 32: VEGETABLE INGESTION PATHWAY
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

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

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

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

         do l = 1, nchain
         avegrt(l) = activity(l)
         enddo

c ============================================================================
c SECTION 33: MILK INGESTION PATHWAY (VIA FORAGE)
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

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

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

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

         do l = 1, nchain
         amilkrt(l) = activity(l)
         enddo

c ============================================================================
c SECTION 34: MEAT INGESTION PATHWAY (VIA FORAGE)
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

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

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

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

         do l = 1, nchain
         ameatrt(l) = activity(l)
         enddo

c ============================================================================
c SECTION 35: FOOD TRANSFER FACTOR CALCULATIONS
c ============================================================================
c Calculate nuclide concentration in food media (vegetables, milk, meat)
c using element-specific transfer models:
c   - Tritium (H-3): water-based model (high transfer to biota)
c   - Carbon-14: carbon cycle model (high transfer to biota)
c   - Other nuclides: NCRP Commentary 3 Equation B2 model
c ============================================================================

         read(radname(n,m),'(a2)') NN

         do l = 1, nchain

          if(NN.eq.'H-') then !H-3 (tritium) -  D. A. Baker, et al, EPA 600/9-76-016, July 1976.
          cfveg(l) = (conc(n,m)*9.0*FHVEG/CH2O)
          cfmilk(l) = (conc(n,m)*9.0*FHFEED/CH2O)*QMILK*FHMILK/(FHFEED*QMILK + QWMILK/9.0)
          cfmeat(l) = (conc(n,m)*9.0*FHFEED/CH2O)*QMEAT*FHMEAT/(FHFEED*QMEAT + QWMEAT/9.0)

          elseif(NN.eq.'C-') then !C14 - D. A. Baker, et al, EPA 600/9-76-016, July 1976.
          cfveg(l) = conc(n,m)*(FCVEG/CCARB)
          cfmilk(l) = conc(n,m)*(FCMILK/CCARB)
          cfmeat(l) = conc(n,m)*(FCMEAT/CCARB)

          else !Equation B2 of NCRP Commentary 3
          cfveg(l) =  conc(n,m) * ((FRVEG/YVEG)*aveglf(l)+(xbveg(l)/PVEG)*avegrt(l))
          cfmilk(l) = conc(n,m) * ((FRMILK/YMILK)*amilklf(l)+(xbfor(l)/PMILK)*amilkrt(l))*xfmilk(l)*QMILK
          cfmeat(l) = conc(n,m) * ((FRMEAT/YMEAT)*ameatlf(l)+(xbfor(l)/PMEAT)*ameatrt(l))*xfmeat(l)*QMEAT
          endif

         enddo

c ============================================================================
c SECTION 36: INGESTION DOSE CALCULATIONS - VEGETABLES, MILK, MEAT
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

           CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

            do l = 1, nchain
            dveg(l) = activity(l) * ingestionfac(l) * CONVEG 
            enddo


           do l = 1, nchain !loop over radionuclide and progeny
           sourceterm(l) = 0.0
           initialamt(l) = cfmilk(l)
           activity(l) = 0.0
           enddo

           xremove = 0.0
           xtime = DELMILK 

           CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

            do l = 1, nchain
            dmilk(l) = activity(l) * ingestionfac(l) * CONMILK 
            enddo

 
           do l = 1, nchain !loop over radionuclide and progeny
           sourceterm(l) = 0.0
           initialamt(l) = cfmeat(l)
           activity(l) = 0.0
           enddo

           xremove = 0.0
           xtime = DELMEAT 

           CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

            do l = 1, nchain
            dmeat(l) = activity(l) * ingestionfac(l) * CONMEAT 
            enddo

c ============================================================================
c SECTION 37: ACCUMULATE TOTAL INGESTION PATHWAY DOSES
c ============================================================================

          vegdose(n,m) = 0.0
          milkdose(n,m) = 0.0
          meatdose(n,m) = 0.0
          do l = 1, nchain
          vegdose(n,m) = vegdose(n,m) + dveg(l)
          milkdose(n,m) = milkdose(n,m) + dmilk(l)
          meatdose(n,m) = meatdose(n,m) + dmeat(l)
          enddo

c ============================================================================
c SECTION 38: REPORT KEY RECEPTOR LOCATIONS
c ============================================================================
c For specific receptor locations defined in input (vegetable garden,
c milk cow, meat animal locations), report the computed doses. Also
c accumulate total dose across all pathways and locations.
c ============================================================================

       if(dist.eq.distanceveg(n).and.vegsector(n).eq.sectorname(j)) then !vegetable receptor
        alldose = alldose + vegdose(n,m)
        totvegdose = vegdose(n,m)
        write(*,*)'Vegetable dose:',radname(n,m),vegdose(n,m),dist,sectorname(j)
        endif
       
       if(dist.eq.distancecow(n).and.milksector(n).eq.sectorname(j)) then !milk receptor
        alldose = alldose + milkdose(n,m)
        totmilkdose = milkdose(n,m)
        write(*,*)'Milk dose:',radname(n,m),milkdose(n,m),dist,sectorname(j)
       endif

       if(dist.eq.distancemeat(n).and.meatsector(n).eq.sectorname(j)) then !meat receptor
        alldose = alldose + meatdose(n,m)
        totmeatdose = meatdose(n,m)
        write(*,*)'Meat dose:',radname(n,m),meatdose(n,m),dist,sectorname(j)
       endif
 
c ============================================================================
c SECTION 39: CALCULATE MAXIMUM DOSE RECEPTOR
c ============================================================================
c Find the receptor location (distance and sector) receiving the maximum
c direct exposure dose (immersion + inhalation + surface deposition).
c This is the critical dose location for regulatory comparison.
c ============================================================================

       if(dist.eq.windbinsrdist(j)) then 
            tempdose = totimdose(n,m) + totihdose(n,m) + totsurfdose(n,m)
        if(tempdose.gt.maxdose) then 
            maxdose = tempdose
            maxsect = j
            maxdist = dist
            surface = totsurfdose(n,m)
            immersion = totimdose(n,m)
            inhalation = totihdose(n,m)
            maxconc = conc(n,m)
        endif
       endif

c ============================================================================
c SECTION 40: WRITE DETAILED DIAGNOSTIC OUTPUT
c ============================================================================
c Write comma-separated values to diagnostic file for all distance/sector
c combinations, including concentrations and all dose components
c
c Concentration in pCi/m3 and dose in mrem/year
c ============================================================================

       write(13,3000) n,radname(n,m),sectorname(j),dist,
     &     1.0E+12*conc(n,m),totimdose(n,m),totihdose(n,m),totsurfdose(n,m),
     &     meatdose(n,m),milkdose(n,m),vegdose(n,m)

       enddo !end loop over distances between source and receptor (receptor)
       enddo !end loop over number of wind direction bins (j)

c ============================================================================
c SECTION 41: CONSOLE OUTPUT - SUMMARY OF DOSES FOR CURRENT RADIONUCLIDE
c ============================================================================
c Print dose summaries by pathway (inhalation, immersion, surface, ingestion)
c and total dose at the maximum dose receptor location.
c ============================================================================

       alldose = alldose + surface + immersion + inhalation
       write(*,*)'Max dose (',maxdose,') receptor at ',windbinsrdist(maxsect),sectorname(maxsect)
       write(*,*)
       write(*,2999)'Maximum air concentration:',radname(n,m),maxconc
       write(*,2999)'Inhalation dose:',radname(n,m),inhalation
       write(*,2999)'Immersion dose:',radname(n,m),immersion
       write(*,2999)'Surface dose:',radname(n,m),surface
       write(*,2999)'Ingestion meat dose:',radname(n,m),totmeatdose
       write(*,2999)'Ingestion milk dose:',radname(n,m),totmilkdose
       write(*,2999)'Ingestion vegetables dose:',radname(n,m),totvegdose
       write(*,2999)'Total ingestion dose:',radname(n,m),totmeatdose+totmilkdose+totvegdose
       write(*,2999)'Total dose:',radname(n,m),(inhalation+immersion+surface+totmeatdose+totmilkdose+totvegdose)
       write(*,*)'--------------------------------------------------------------------------' 

c write to diagnostic output file

       if(m.eq.1) then
       write(12,2002)
       write(12,1800)'Max dose receptor ',windbinsrdist(maxsect),' m downwind in sector ',sectorname(maxsect)
 1800  FORMAT(a,f8.0,a,a)
       write(12,2002)

       write(12,1998)
       endif

       !Dose units are mrem/yr and concentration units pCi/m3 (convert Ci/m3 to pCi/m3 here)
       write(12,1999) radname(n,m),rademis(n,m),1.0E12*maxconc,
     &     (immersion+inhalation+surface+totmeatdose+totmilkdose+totvegdose),
     &     immersion,inhalation,surface,totmeatdose,totmilkdose,totvegdose

       enddo !end loop over radionuclides (m)

c ============================================================================
c SECTION 42: CONSOLE OUTPUT - TOTAL FACILITY DOSE AND REGULATORY LIMITS
c ============================================================================
c Report total cumulative effective dose equivalent across all stacks and
c radionuclides. Compare to EPA dose limits:
c   - General public: 10 mrem/yr
c   - Iodine: 3 mrem/yr
c ============================================================================

       write(*,*)
       write(*,3100)'Total dose (mrem/yr):',alldose
       write(*,*)'Total dose limit is 10.0 mrem/yr'
       write(*,*)'Iodine limit is 3.0 mrem/yr'
       write(*,*)

       write(12,2001)
       write(12,2050) alldose
       write(12,2001)

       enddo !end loop over stacks (n)

c ============================================================================
c SECTION 43: PROGRAM EXECUTION TIMESTAMP
c ============================================================================
c Capture and write execution date and time for traceability and
c regulatory documentation purposes
c ============================================================================

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

       write(12,2002)
       write(12,2001)
       write(12,2002)
       write(12,2100) date_time(2),date_time(3),date_time(1),
     &           date_time(5),date_time(6)


c ============================================================================
c SECTION 44: FORMAT STATEMENTS FOR REPORT OUTPUT
c ============================================================================

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
 1014 FORMAT('Vertical momentum (YES/NO):',a)
 1015 FORMAT('Building height (m):',f10.1)
 1016 FORMAT('Building width (m):',f10.1)
 1017 FORMAT('Distance from source to receptor (m):',f10.1)
 1018 FORMAT('Rainfall rate (m/yr):',f10.1)
 1019 FORMAT('Number of radionuclides:',I2)

 1998 FORMAT('Radionuclide,Emissions,Concentration,Total Dose,Immersion,Inhalation,Surface,Meat,Milk,Vegetable')
 1999 FORMAT(a,',',E10.3,',',E10.3,7(',',E10.3))

 2001 FORMAT('-----------------------------------------------------------------------')
 2002 FORMAT(' ')

 2050 FORMAT('Total effective dose equivalent all stacks (mrem/yr):',E10.3)
 2100 FORMAT('Date/Time (LST) stamp of code execution (MM/DD/YYYY HH:MM):',i2,'/',i2,'/',i4,' ',i2,':',i2)

 2999  FORMAT(a30,a10,f30.15)
 3000  FORMAT(i2,',',a,',',a,',',f9.1,',',f30.20,6(',',f20.10))
 3001  FORMAT('STACK,RADNAME,SECTOR,DISTANCE,CONC,IMMERSION,INHALATION,SURFACE,MEAT,MILK,VEG')
 3100  FORMAT(a,f15.4)

 4000 FORMAT('COMPLY version 2.0 level 4 distributed through GitHub')

c ============================================================================
c END OF PROGRAM
c ============================================================================

      END 

