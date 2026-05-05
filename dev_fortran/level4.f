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

      REAL*4, ALLOCATABLE :: releaseht(:)           !stack height (m)
      REAL*4, ALLOCATABLE :: diameter(:)            !stack diameter (m)
      REAL*4, ALLOCATABLE :: exitvelocity(:)        !exit velocity (m/s)
      REAL*4, ALLOCATABLE :: airtemp(:)             !air temperature (F)
      REAL*4, ALLOCATABLE :: stacktemp(:)           !stack temperature (F)
      REAL*4, ALLOCATABLE :: bldheight(:)           !building height (m)
      REAL*4, ALLOCATABLE :: bldwidth(:)            !building width (m)
      REAL*4, ALLOCATABLE :: distance(:)            !distance between stack and receptor (m)
      REAL*4, ALLOCATABLE :: distanceveg(:)         !distance between garden and receptor (m)
      REAL*4, ALLOCATABLE :: distancecow(:)         !distance between cows/milk and receptor (m)
      REAL*4, ALLOCATABLE :: distancemeat(:)        !distance between meat and receptor (m)
      REAL*4, ALLOCATABLE :: latitude(:)            !latitude (decimal degrees)
      REAL*4, ALLOCATABLE :: longitude(:)           !longitude (decimal degrees)
      REAL*4, ALLOCATABLE :: pctcalm(:)             !percentage of calm days at wind data location
      REAL*4, ALLOCATABLE :: rainfall(:)            !rainfall rate (m/yr)

c ============================================================================
c SECTION 2: DYNAMIC ARRAY DECLARATIONS - EMISSION & CONCENTRATION DATA
c ============================================================================
c These arrays store radionuclide-specific emission rates and computed
c air/dose concentrations for each stack/radionuclide combination

      REAL*8, ALLOCATABLE :: rademis(:,:)           !radionuclide emissions (Ci/s) for each stack
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

      REAL*4, ALLOCATABLE :: halflife(:)            !radionuclide half life (years)
      REAL*4, ALLOCATABLE :: wbfac(:)               !whole body dose factor from NCRP Commentary 3 table 2
      REAL*4, ALLOCATABLE :: xwbfac(:)              !whole body dose factor from NCRP Commentary 3 table 3
      REAL*4, ALLOCATABLE :: vegfac(:)              !vegetable dose factor from NCRP Commentary 3 table 7
      REAL*4, ALLOCATABLE :: cowfac(:)              !cow dose factor from NCRP Commentary 3
      REAL*8, ALLOCATABLE :: inhalationfac(:)       !inhalation dose factor from ICRP table 2.1  
      REAL*8, ALLOCATABLE :: immersionfac(:)        !immersion dose factor from NCRP Commentary Table C5
      REAL*8, ALLOCATABLE :: surfacefac(:)          !surface dose factor from NCRP Commentary Table C6
      REAL*8, ALLOCATABLE :: ingestionfac(:)        !ingestion dose factor from ICRP Table 2.2
      REAL*8, ALLOCATABLE :: drydepvel(:)           !dry deposition velocity (m/day) from NCRP Commentary Table L4

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

      REAL*8, ALLOCATABLE :: windbinfreq(:,:)       !Sector (wind dir) frequency from that direction
      REAL*8, ALLOCATABLE :: windbinspeed(:,:)      !Sector (wind dir) average speed from that direction (m/s)
      REAL*8, ALLOCATABLE :: windbinsrdist(:,:)     !Sector (wind dir) source-receptor distance (m)

c ============================================================================
c SECTION 10: DYNAMIC ARRAY DECLARATIONS - INDEXING & METADATA
c ============================================================================

      INTEGER, ALLOCATABLE :: nrads(:)              !number of emitted radionucludes at each stack
      INTEGER, ALLOCATABLE :: nwindbins(:)          !number of wind direction bins at each stack

      CHARACTER*8,  ALLOCATABLE :: radname(:,:)     !radionuclide names at each stack
      CHARACTER*3,  ALLOCATABLE :: windbinname(:,:) !wind direction sector 
      CHARACTER*3,  ALLOCATABLE :: sectorname(:,:)  !sector name for receptors
      CHARACTER*3,  ALLOCATABLE :: vegsector(:)     !Sector (wind dir) location of vegetables
      CHARACTER*3,  ALLOCATABLE :: milksector(:)    !Sector (wind dir) location of cow/milk
      CHARACTER*3,  ALLOCATABLE :: meatsector(:)    !Sector (wind dir) location of meat
      CHARACTER*8,  ALLOCATABLE :: decaychain(:)    !Decay chain list
      CHARACTER*8,  ALLOCATABLE :: decaychainname(:)!Decay chain radionuclide names
      CHARACTER*3,  ALLOCATABLE :: verticalmomentum(:)    !vertical momentum (YES or NO)

      CHARACTER*200,ALLOCATABLE :: windsource(:)    !Source of wind data
      CHARACTER*200,ALLOCATABLE :: winddates(:)     !Date coverage of wind data
      CHARACTER*200,ALLOCATABLE :: windloc(:)       !Location of wind data
      CHARACTER*200,ALLOCATABLE :: winddist(:)      !Distance from wind data to source

c ============================================================================
c SECTION 11: STATIC VARIABLE DECLARATIONS - DISPERSION CALCULATION TERMS
c ============================================================================
c These variables are used in the Gaussian plume model calculations for
c atmospheric dispersion and effective stack height corrections

      REAL FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF,FACTORG,FACTORH,FACTORI
      REAL FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN
      REAL area,xcrit,SIGMAZ,CAPSIG,pp,alldose,idose,dist,altconc
      REAL tempdose,maxdose,maxdist,totmilkdose,totvegdose,itotvegdose,totmeatdose,itotcowdose
      REAL capp,capn0,outprod,inprod,z1,term1,term2,sum
      REAL transittime                              !transit time of air from emission release to receptor (s)
      REAL fdux                                     !flux term units are 1/(m2/s) or s/m2
      REAL wetdepadj,wetdepadjbad                   !wet deposition adjustment factor
      REAL xtime
      REAL xremove
      REAL ingestion, immersion, inhalation, surface
      REAL RHOAIR,RHOSTK,STKVEL,R,FM,DELHMO,DELHBU,HEFF,FB,zz,INHFAC,TINIT,T

c ============================================================================
c SECTION 12: STATIC VARIABLE DECLARATIONS - CONTROL & LOOP INDICES
c ============================================================================

      INTEGER nstacks                               !number of unique stacks
      INTEGER receptor, N, M, I, J, K, L, II, JJ, P, Q                        !loop counters
      INTEGER nchain                                !number of nuclides in the decay series
      INTEGER maxsect

      CHARACTER*200 title,company,facility,address,poc,phone,inline
      CHARACTER*8 NUC
      CHARACTER*2 NN

      integer date_time(8)                      !for date stamp function
      character*10 b(3)                         !for date stamp function


c ============================================================================
c SECTION 13: FILE OPERATIONS - OPEN INPUT AND OUTPUT FILES
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
c SECTION 14: INPUT PARSING - READ FACILITY IDENTIFICATION AND STACK COUNT
c ============================================================================
c Read project metadata and determine the number of stacks to process

       read(11,'(20x,a)') title
       read(11,'(20x,a)') company
       read(11,'(20x,a)') facility
       read(11,'(20x,a)') address
       read(11,'(20x,a)') poc
       read(11,'(20x,a)') phone

       read(11,'(20x,i3)') nstacks

c ============================================================================
c SECTION 15: MEMORY ALLOCATION - ALLOCATE STACK-LEVEL PARAMETER ARRAYS
c ============================================================================
c Allocate arrays to store stack-specific parameters (heights, distances,
c temperatures, etc.) based on the number of stacks

       allocate(latitude(nstacks))
       allocate(longitude(nstacks))
       allocate(releaseht(nstacks))
       allocate(bldheight(nstacks))
       allocate(bldwidth(nstacks))
       allocate(exitvelocity(nstacks))
       allocate(diameter(nstacks))
       allocate(stacktemp(nstacks))
       allocate(airtemp(nstacks))
       allocate(verticalmomentum(nstacks))
       allocate(distance(nstacks))
       allocate(distanceveg(nstacks))
       allocate(distancecow(nstacks))
       allocate(distancemeat(nstacks))
       allocate(vegsector(nstacks))
       allocate(milksector(nstacks))
       allocate(meatsector(nstacks))
       allocate(nrads(nstacks))
       allocate(windsource(nstacks))
       allocate(winddates(nstacks))
       allocate(windloc(nstacks))
       allocate(winddist(nstacks))
       allocate(pctcalm(nstacks))
       allocate(nwindbins(nstacks))
       allocate(rainfall(nstacks))

c ============================================================================
c SECTION 16: INPUT PARSING LOOP - READ STACK-SPECIFIC PARAMETERS
c ============================================================================
c Main loop to read input parameters for each stack source. Includes
c geometrical parameters, temperatures, distances to receptors, and
c number of emitted radionuclides

       do n = 1, nstacks 

       read(11,'(a)') !blank spacer to improve readability of the input file
       read(11,'(20x,f10.5)') latitude(n)
       read(11,'(20x,f10.5)') longitude(n)
       read(11,'(20x,f10.5)') releaseht(n)
       read(11,'(20x,f10.5)') bldheight(n)
       read(11,'(20x,f10.5)') bldwidth(n)
       read(11,'(20x,f10.5)') exitvelocity(n)
       read(11,'(20x,f10.5)') diameter(n)
       read(11,'(20x,f10.5)') stacktemp(n)
       read(11,'(20x,f10.5)') airtemp(n)
       read(11,'(20x,a)') verticalmomentum(n)
       read(11,'(20x,f10.5)') distance(n)
       read(11,'(20x,f10.5)') distanceveg(n)
       read(11,'(20x,f10.5)') distancecow(n)
       read(11,'(20x,f10.5)') distancemeat(n)

       read(11,'(20x,a)') vegsector(n)
       read(11,'(20x,a)') milksector(n)
       read(11,'(20x,a)') meatsector(n)

       read(11,'(20x,f10.5)') rainfall(n)

       read(11,'(20x,i3)') nrads(n)

c ============================================================================
c SECTION 17: MEMORY ALLOCATION - RADIONUCLIDE & DOSE FACTOR ARRAYS
c ============================================================================
c Allocate arrays for radionuclide emissions, concentrations, doses, and
c dose factors. Note: These are allocated within the stack loop but assume
c a maximum of 19 nuclides in decay chain

       allocate(radname(nstacks,nrads(n)))
       allocate(rademis(nstacks,nrads(n)))
       allocate(conc(nstacks,nrads(n)))
       allocate(dose(nstacks,nrads(n)))
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

       allocate(totimdose(nstacks,nrads(n)))
       allocate(totihdose(nstacks,nrads(n)))
       allocate(totsurfdose(nstacks,nrads(n)))
       allocate(totingdose(nstacks,nrads(n)))

       allocate(milkdose(nstacks,nrads(n)))
       allocate(vegdose(nstacks,nrads(n)))
       allocate(meatdose(nstacks,nrads(n)))

c ============================================================================
c SECTION 18: INPUT PARSING - READ RADIONUCLIDE EMISSION DATA
c ============================================================================
c For each stack, read the list of radionuclides emitted and their
c emission rates (in Ci/s)

       do m = 1, nrads(n)
       read(11,'(20x,a)') inline
       read(inline,*) radname(n,m),rademis(n,m)
       enddo !end loop over reading in radionuclide emissions

c ============================================================================
c SECTION 19: INPUT PARSING - READ METEOROLOGICAL DATA SOURCE & WIND SECTORS
c ============================================================================
c Read metadata about wind data source and sector-based wind statistics
c (frequency, speed, source-receptor distance) for each wind direction bin

       read(11,'(20x,a)') windsource(n)
       read(11,'(20x,a)') winddates(n)
       read(11,'(20x,a)') windloc(n)
       read(11,'(20x,a)') winddist(n)

       read(11,'(20x,f10.5)') pctcalm(n)
       read(11,'(20x,i3)') nwindbins(n)

       allocate(windbinname(nstacks,nwindbins(n)))
       allocate(windbinfreq(nstacks,nwindbins(n)))
       allocate(windbinspeed(nstacks,nwindbins(n)))
       allocate(windbinsrdist(nstacks,nwindbins(n)))
       allocate(sectorname(nstacks,nwindbins(n)))

       do m = 1, nwindbins(n)
       read(11,'(20x,a)') inline
       read(inline,*) windbinname(n,m),windbinfreq(n,m),windbinspeed(n,m),windbinsrdist(n,m),sectorname(n,m)
       enddo !end loop over reading in wind data and source-receptor distances

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

       area = bldwidth(n) * bldheight(n)
       xcrit = 2.5 * sqrt(area)

       do m = 1, nrads(n)

c ============================================================================
c SECTION 23: INITIALIZE DOSE ACCUMULATORS FOR CURRENT RADIONUCLIDE
c ============================================================================

       maxdose = 0.0
       inhalation = 0.
       immersion = 0.
       ingestion = 0.
       totmilkdose = 0.
       totvegdose = 0.
       totmeatdose = 0.

c ============================================================================
c SECTION 24: WIND SECTOR AND DISTANCE ITERATION LOOPS
c ============================================================================
c Loop over wind direction sectors and receptor distances. For each combination,
c calculate air concentration and dose contributions across all decay pathways

       do j = 1, nwindbins(n)

       do receptor = 10, 3010, 10

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

       if(exitvelocity(n).le.0.) then
       exitvelocity(n) = 0.3
       endif

       HEFF = releaseht(n)

       IF ( releaseht(n) .GT. 2.5*bldheight(n) ) THEN

c --- Temperature and density calculations ---
        RHOAIR = 0.0735 * ( 80.0 + 460.0 ) / ( airtemp(n) + 460.0 )
        RHOAIR = RHOAIR * 454.0 / ( 0.3048 ** 3 )
        RHOSTK = RHOAIR * ( airtemp(n) + 460.0 ) / ( stacktemp(n) + 460.0 )

c --- Momentum rise calculations ---
        AREA   = 3.14 * ( diameter(n) ** 2.0 ) / 4.0
        STKVEL = exitvelocity(n) / AREA
        R      = STKVEL / windbinspeed(n,j)
        FM     = ( RHOSTK / RHOAIR ) * STKVEL * exitvelocity(n) / 3.14
        TERM1  = 75.0 * FM / ( ( ( 1.0 + 3.0 / R ) * windbinspeed(n,j) ) ** 2.0 )
        TERM2  = releaseht(n) + 7.79 * SQRT( FM ) / ( ( 1.0 + 3.0 / R ) * windbinspeed(n,j) )

c --- Apply momentum rise if vertical momentum enabled ---
        IF ( verticalmomentum(n).eq.'YES' ) THEN
         DELHMO = 0.93 * ( TERM1 ** ( 3.0 / 7.0 ) ) *
     &        ( TERM2 ** ( 1.0 / 7.0 ) )
        ELSE
         DELHMO = 0.0
        ENDIF

c --- Buoyancy rise calculations ---
        FB = 9.81 * ( RHOAIR - RHOSTK ) * exitvelocity(n) / ( 3.14 * RHOAIR )
        DELHBU = 23.67 * ( ( FB / ( windbinspeed(n,j) ** 3.0 ) ) ** 0.6 ) *
     &           ( releaseht(n) + 195.1 * FB / ( windbinspeed(n,j) ** 3.0 ) ) ** 0.4

c --- Use maximum of buoyancy or momentum rise ---
       HEFF = releaseht(n) + AMAX1( DELHMO, DELHBU )

       ENDIF !end condition where H > 0

c ============================================================================
c SECTION 26: GAUSSIAN PLUME MODEL - CONCENTRATION AT RECEPTOR
c ============================================================================
c Calculate air concentration at receptor distance using Gaussian plume
c model with lateral (sigma_z) and source area (CAPSIG) corrections.
c Includes near-field correction for distances less than xcrit.
c ============================================================================

       SIGMAZ = (0.06 * dist) / SQRT( 1.0 + 0.0015 * dist)
       CAPSIG = SQRT( SIGMAZ ** 2 + area / 3.14 )
       zz = 0.5 * ( ( HEFF / CAPSIG ) ** 2 )
        if(zz.gt.0.001.and.releaseht(n).gt.2.5*bldheight(n)) then
        pp = 2.032 * EXP( -zz ) / ( dist * CAPSIG )
        else
        pp = 2.032 / ( dist * CAPSIG )
        endif

       transittime = (dist / windbinspeed(n,j))*(1./3600.)*(1./24.) ! seconds to days
       fdux = (windbinfreq(n,j)/(1.0-pctcalm(n))) / (windbinspeed(n,j)*dist) !s/m2 units

c --- Near-field vs far-field concentration calculation ---
       if(dist.lt.xcrit) then
       conc(n,m) = ( (windbinfreq(n,j)/(1.0-pctcalm(n))) *rademis(n,m)) /
     &          (3.14*windbinspeed(n,j)*1.0*amin1(bldwidth(n),bldheight(n)))
       else
       conc(n,m) = ((windbinfreq(n,j)/(1.0-pctcalm(n))) *rademis(n,m)*pp) / windbinspeed(n,j)
       endif

       if(conc(n,m).lt. 1.0E-32 ) then
         conc(n,m) = 1.0E-32
       endif

       altconc = (windbinfreq(n,j)*rademis(n,m)*pp) / windbinspeed(n,j) !just for QA and not being used for anything else here

c ============================================================================
c SECTION 27: WET DEPOSITION ADJUSTMENT FACTOR
c ============================================================================
c Account for removal of radioactive material via precipitation (wet
c deposition) combined with dry deposition velocity. This adjusts the
c ground-level deposition for rainfall conditions.
c
c Rainfall rate is in m/yr, emissions in Ci/s, concentration in Ci/m3
c ============================================================================

       wetdepadj = (8.0 * fdux * rainfall(n) * 1.0E-05 * (rademis(n,m) * 3600. * 24.))/(3.1416*conc(n,m) )
       wetdepadjbad = (8.0 * fdux * rainfall(n) * 1.0E-05 * (rademis(n,m) * 3600. * 24.))/(3.1416*conc(n,m)/rademis(n,m) )

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

c ============================================================================
c SECTION 29: INHALATION & IMMERSION DOSE CALCULATION
c ============================================================================
c Calculate doses from inhalation of contaminated air and immersion in
c contaminated air. These are pathway doses based on air concentration
c and are independent of ground deposition.
c ============================================================================

         do l = 1, nchain !loop over radionuclide and progeny
         sourceterm(l) = 0.0
         initialamt(l) = 0.0
         activity(l) = 0.0
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
         ihdose(l) = ambientactivity(l) * inhalationfac(l) * 8000. !BREATH variable = 8000 m3/yr breathing rate
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

         xremove = 1.5000000E-02 / 365. !from original comply
         xtime = 100.0 * 365.  !100 days is the dry dep time from original comply

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
         initialamt(l) = 0.0
         activity(l) = 0.0
         enddo

         xremove = 0.693 / 12. ! TWEATH variable 
         xtime = 60. ! TEVEG variable

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

         do l = 1, nchain
         aveglf(l) = activity(l) 
         enddo


         do l = 1, nchain
         sourceterm(l) = normsourceactivity(l)
         initialamt(l) = 0.0
         activity(l) = 0.0
         enddo

         xremove = 1.5000000E-02 / 365.  ! REMCON variable convert 1/yr to 1/day
         xtime = 100. * 365 ! DEPTIM variable convert yr to day

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

         xremove = 0.693 / 12. ! TWEATH variable
         xtime = 30. ! TEMILK variable

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

         do l = 1, nchain
         amilklf(l) = activity(l)
         enddo

         do l = 1, nchain
         sourceterm(l) = normsourceactivity(l)
         initialamt(l) = 0.0
         activity(l) = 0.0
         enddo

         xremove = 1.5000000E-02 / 365. ! REMCON variable
         xtime = 100. * 365. ! DEPTIM variable

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

         xremove = 0.693 / 12. ! TWEATH variable
         xtime = 30. ! TEMEAT variable

      CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

         do l = 1, nchain
         ameatlf(l) = activity(l)
         enddo

         do l = 1, nchain
         sourceterm(l) = normsourceactivity(l)
         initialamt(l) = 0.0
         activity(l) = 0.0
         enddo
 
         xremove = 1.5000000E-02 / 365.  ! REMCON variable
         xtime = 100. * 365. ! DEPTIM variable

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
          cfveg(l) = conc(n,m)*9.0*(0.1/0.8E-02) ! (FHVEG / CH2O )
          cfmilk(l) = (conc(n,m)*9.0*7.0E-02/8.0E-03)*16.0*0.11/(7.0E-02*16. + 220.0/9.0)
          cfmeat(l) = (conc(n,m)*9.0*7.0E-02/8.0E-03)*12.0*0.10/(7.0E-02*12. + 50.0/9.0) 

          elseif(NN.eq.'C-') then !C14 - D. A. Baker, et al, EPA 600/9-76-016, July 1976.
          cfveg(l) = conc(n,m)*(9.000000E-02/0.16E-03)
          cfmilk(l) = conc(n,m)*(7.0000000E-02/0.16E-03)
          cfmeat(l) = conc(n,m)*(0.2400000/0.16E-03)

          else !Equation B2 of NCRP Commentary 3
          cfveg(l) =  conc(n,m) * ((0.1/1.0)*aveglf(l)+(xbveg(l)/220.)*avegrt(l))
          cfmilk(l) = conc(n,m) * ((1.8/1.0)*amilklf(l)+(xbfor(l)/220.)*amilkrt(l))*xfmilk(l)*16.0
          cfmeat(l) = conc(n,m) * ((1.8/1.0)*ameatlf(l)+(xbfor(l)/220.)*ameatrt(l))*xfmeat(l)*12.0
          endif

          enddo

c ============================================================================
c SECTION 36: INGESTION DOSE CALCULATIONS - VEGETABLES, MILK, MEAT
c ============================================================================
c For each food pathway:
c   - Decay concentration in food from harvest to consumption
c   - Apply ingestion dose conversion factor
c   - Multiply by consumption rate (ingestion factor)
c
c Time delays account for:
c   - Vegetables: ~11 days (harvest to consumption)
c   - Milk: ~2 days (production to consumption)
c   - Meat: ~17 days (slaughter to consumption)
c ============================================================================

           do l = 1, nchain !loop over radionuclide and progeny
           sourceterm(l) = 0.0
           initialamt(l) = cfveg(l)
           activity(l) = 0.0
           enddo

           xremove = 0.0
           xtime = 11. !DELVEG variable

           CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

            do l = 1, nchain
            dveg(l) = activity(l) * ingestionfac(l) * 70. !CONVEG = 70
            enddo


           do l = 1, nchain !loop over radionuclide and progeny
           sourceterm(l) = 0.0
           initialamt(l) = cfmilk(l)
           activity(l) = 0.0
           enddo

           xremove = 0.0
           xtime = 2. !DELMILK variable

           CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

            do l = 1, nchain
            dmilk(l) = activity(l) * ingestionfac(l) * 160. !CONMILK = 160
            enddo

 
           do l = 1, nchain !loop over radionuclide and progeny
           sourceterm(l) = 0.0
           initialamt(l) = cfmeat(l)
           activity(l) = 0.0
           enddo

           xremove = 0.0
           xtime = 17. !DELMEAT variable

           CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

            do l = 1, nchain
            dmeat(l) = activity(l) * ingestionfac(l) * 75. !CONMEAT = 75
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

       if(dist.eq.distanceveg(n).and.vegsector(n).eq.sectorname(n,j)) then !vegetable receptor
        alldose = alldose + vegdose(n,m)
        totvegdose = vegdose(n,m)
        write(*,*)'Vegetable dose:',radname(n,m),vegdose(n,m),dist,sectorname(n,j)
        endif
       
       if(dist.eq.distancecow(n).and.milksector(n).eq.sectorname(n,j)) then !milk receptor
        alldose = alldose + milkdose(n,m)
        totmilkdose = milkdose(n,m)
        write(*,*)'Milk dose:',radname(n,m),milkdose(n,m),dist,sectorname(n,j)
       endif

       if(dist.eq.distancemeat(n).and.meatsector(n).eq.sectorname(n,j)) then !meat receptor
        alldose = alldose + meatdose(n,m)
        totmeatdose = meatdose(n,m)
        write(*,*)'Meat dose:',radname(n,m),meatdose(n,m),dist,sectorname(n,j)
       endif
 
c ============================================================================
c SECTION 39: CALCULATE MAXIMUM DOSE RECEPTOR
c ============================================================================
c Find the receptor location (distance and sector) receiving the maximum
c direct exposure dose (immersion + inhalation + surface deposition).
c This is the critical dose location for regulatory comparison.
c ============================================================================

       if(dist.eq.windbinsrdist(n,j)) then 
       write(*,*)'Air concentration:',radname(n,m),conc(n,m),sectorname(n,j),dist
            tempdose = totimdose(n,m) + totihdose(n,m) + totsurfdose(n,m)
        if(tempdose.gt.maxdose) then 
            maxdose = tempdose
            maxsect = j
            maxdist = dist
            surface = totsurfdose(n,m)
            immersion = totimdose(n,m)
            inhalation = totihdose(n,m)
        endif
       endif

c ============================================================================
c SECTION 40: WRITE DETAILED DIAGNOSTIC OUTPUT
c ============================================================================
c Write comma-separated values to diagnostic file for all distance/sector
c combinations, including concentrations and all dose components
c ============================================================================

       write(13,3000) n,radname(n,m),sectorname(n,j),dist,
     &     conc(n,m),totimdose(n,m),totihdose(n,m),totsurfdose(n,m),
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
       write(*,*)'Max dose (',maxdose,') receptor at ',windbinsrdist(n,maxsect),sectorname(n,maxsect)
       write(*,*)
       write(*,2999)'Inhalation dose:',radname(n,m),inhalation
       write(*,2999)'Immersion dose:',radname(n,m),immersion
       write(*,2999)'Surface dose:',radname(n,m),surface
       write(*,2999)'Ingestion meat dose:',radname(n,m),totmeatdose
       write(*,2999)'Ingestion milk dose:',radname(n,m),totmilkdose
       write(*,2999)'Ingestion vegetables dose:',radname(n,m),totvegdose
       write(*,2999)'Total ingestion dose:',radname(n,m),totmeatdose+totmilkdose+totvegdose
       write(*,2999)'Total dose:',radname(n,m),(inhalation+immersion+surface+totmeatdose+totmilkdose+totvegdose)
       write(*,*)'--------------------------------------------------------------------------' 
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

       enddo !end loop over stacks (n)

c ============================================================================
c SECTION 43: FORMAT STATEMENTS FOR OUTPUT
c ============================================================================

 2999  FORMAT(a30,a10,f20.15)
 3000  FORMAT(i2,',',a,',',a,',',f6.1,',',f30.20,6(',',f20.10))
 3001  FORMAT('STACK,RADNAME,SECTOR,DISTANCE,CONC,IMMERSION,INHALATION,SURFACE,MEAT,MILK,VEG')
 3100  FORMAT(a,f10.4)

c ============================================================================
c SECTION 44: WRITE STANDARDIZED OUTPUT REPORT
c ============================================================================
c Generate human-readable output report with facility identification,
c stack parameters, dose results by radionuclide and pathway, and
c execution timestamp
c ============================================================================

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
       write(12,1011) exitvelocity(n)
       write(12,1012) stacktemp(n)
       write(12,1013) airtemp(n)
       write(12,1014) verticalmomentum(n)
       write(12,1015) bldheight(n)
       write(12,1016) bldwidth(n)
       write(12,1017) distance(n)
       write(12,1018) rainfall(n)
       write(12,1019) nrads(n)

       write(12,2002)
       write(12,1800)'Max dose receptor ',windbinsrdist(n,maxsect),' m downwind in sector ',sectorname(n,maxsect)
 1800  FORMAT(a,f5.0,a,a)
       write(12,2002)

       write(12,1998)
       do m=1, nrads(n)
       write(12,1999) radname(n,m),rademis(n,m),conc(n,m),
     &     (totimdose(n,m)+totihdose(n,m)+totsurfdose(n,m)+meatdose(n,m)+milkdose(n,m)+vegdose(n,m)),
     &     totimdose(n,m),totihdose(n,m),totsurfdose(n,m),meatdose(n,m),milkdose(n,m),vegdose(n,m)

       enddo !end loop over radionuclides
       enddo !end loop over stacks

       write(12,2002)
       write(12,2001)
       write(12,2001)
       write(12,2002)
       write(12,2050) alldose
       write(12,2002)

c ============================================================================
c SECTION 45: PROGRAM EXECUTION TIMESTAMP
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
c SECTION 46: FORMAT STATEMENTS FOR REPORT OUTPUT
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
 1011 FORMAT('Exit velocity (m3/s):',f10.1)
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
 4000 FORMAT('COMPLY version 2.0 level 4 distributed through GitHub')

c ============================================================================
c END OF PROGRAM
c ============================================================================

      END 

