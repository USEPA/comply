        PROGRAM COMPLY
c
c ============================================================================
c  COMPLY v2.0 (Optimized) - Radionuclide Compliance Dose Assessment Model
c ============================================================================
c
c PURPOSE:
c   This program calculates radiation doses from atmospheric radionuclide
c   emissions for compliance with 10 CFR 61 and NESHAP regulations.
c   It models air concentrations, deposition, and ingestion pathways.
c
c METHODOLOGY:
c   - Uses Gaussian plume model for atmospheric transport
c   - Incorporates stack parameters (height, diameter, temperature)
c   - Accounts for buoyant and momentum rise of exhaust plumes
c   - Models decay chain transformations during transport
c   - Calculates inhalation, immersion, surface, and ingestion pathways
c
c INPUT FILES:
c   - input.level4.txt: Main input parameters (stacks, receptors, emissions)
c   - comply171_original_input_data.csv: Radionuclide dose factors
c   - comply171_original_decaychain_data.csv: Decay chain definitions
c
c OUTPUT FILES:
c   - output.level4.txt: Formatted compliance dose report
c   - diag.level4.txt: Diagnostic results at all distances/sectors
c
c AUTHOR: U.S. EPA
c ============================================================================
c
        IMPLICIT NONE
c
c ============================================================================
c SECTION 1: DYNAMIC ARRAY DECLARATIONS FOR STACK PARAMETERS
c ============================================================================
c
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
c
c ============================================================================
c SECTION 2: DYNAMIC ARRAY DECLARATIONS FOR CONCENTRATION AND DOSE DATA
c ============================================================================
c
        REAL*8, ALLOCATABLE :: rademis(:,:)           !radionuclide emissions (Ci/s) for each stack
        REAL*8, ALLOCATABLE :: conc(:,:)              !radionuclide air concentration (Bq/m3)
        REAL*8, ALLOCATABLE :: dose(:,:)              !radionuclide dose (mrem/yr)
c
c ============================================================================
c SECTION 3: DYNAMIC ARRAY DECLARATIONS FOR DECAY CHAIN AND SOURCE TERMS
c ============================================================================
c These arrays track activity levels at different stages of the calculation
c
        REAL*8, ALLOCATABLE :: sourceterm(:)          !source term activity
        REAL*8, ALLOCATABLE :: initialamt(:)          !initial amount of radionuclide
        REAL*8, ALLOCATABLE :: activity(:)            !activity after decay calculations
        REAL*8, ALLOCATABLE :: surfaceactivity(:)     !activity on ground surface
        REAL*8, ALLOCATABLE :: ambientactivity(:)     !activity in air
        REAL*8, ALLOCATABLE :: normsourceactivity(:)  !normalized source activity for ingestion
c
c ============================================================================
c SECTION 4: DYNAMIC ARRAY DECLARATIONS FOR INGESTION PATHWAY FACTORS
c ============================================================================
c These arrays store intermediate calculations for vegetables, milk, and meat
c pathways. The 'lf' suffix means "leaf" and 'rt' means "root"
c
        REAL*8, ALLOCATABLE :: aveglf(:)              !vegetation activity (leaf component)
        REAL*8, ALLOCATABLE :: avegrt(:)              !vegetation activity (root component)
        REAL*8, ALLOCATABLE :: amilklf(:)             !milk activity (leaf-based feed component)
        REAL*8, ALLOCATABLE :: amilkrt(:)             !milk activity (root-based feed component)
        REAL*8, ALLOCATABLE :: ameatlf(:)             !meat activity (leaf-based feed component)
        REAL*8, ALLOCATABLE :: ameatrt(:)             !meat activity (root-based feed component)
        REAL*8, ALLOCATABLE :: cfveg(:)               !concentration factor for vegetables
        REAL*8, ALLOCATABLE :: cfmilk(:)              !concentration factor for milk
        REAL*8, ALLOCATABLE :: cfmeat(:)              !concentration factor for meat
        REAL*8, ALLOCATABLE :: dveg(:)                !dose from vegetable ingestion
        REAL*8, ALLOCATABLE :: dmilk(:)               !dose from milk ingestion
        REAL*8, ALLOCATABLE :: dmeat(:)               !dose from meat ingestion
        REAL*8, ALLOCATABLE :: ingdose(:)             !total ingestion dose
c
c ============================================================================
c SECTION 5: DYNAMIC ARRAY DECLARATIONS FOR DOSE FACTOR COEFFICIENTS
c ============================================================================
c These coefficients are read from reference databases and used to calculate
c dose from activity (dose = activity * dose_factor)
c
        REAL*8, ALLOCATABLE :: kvalue(:)              !K coefficient for decay calculations
        REAL*8, ALLOCATABLE :: ktot(:)                !total K coefficient
        REAL*8, ALLOCATABLE :: xbveg(:)               !bioaccumulation factor vegetables
        REAL*8, ALLOCATABLE :: xbfor(:)               !bioaccumulation factor forage
        REAL*8, ALLOCATABLE :: xfmilk(:)              !transfer factor milk
        REAL*8, ALLOCATABLE :: xfmeat(:)              !transfer factor meat
c
        REAL*4, ALLOCATABLE :: halflife(:)            !radionuclide half life (years)
        REAL*4, ALLOCATABLE :: wbfac(:)               !whole body dose factor from NCRP Commentary 3 table 2
        REAL*4, ALLOCATABLE :: xwbfac(:)              !whole body dose factor from NCRP Commentary 3 table 3
        REAL*4, ALLOCATABLE :: vegfac(:)              !vegetable dose factor from NCRP Commentary 3 table 7
        REAL*4, ALLOCATABLE :: cowfac(:)              !cow dose factor from NCRP Commentary 3
        REAL*8, ALLOCATABLE :: inhalationfac(:)       !inhalation dose factor from ICRP table 2.1  
        REAL*8, ALLOCATABLE :: immersionfac(:)        !immersion dose factor from NCRP Commentary Table C5
        REAL*8, ALLOCATABLE :: surfacefac(:)          !surface dose factor from NCRP Commentary Table C6
        REAL*8, ALLOCATABLE :: ingestionfac(:)        !ingestion dose factor from ICRP Table 2.2
        REAL*8, ALLOCATABLE :: drydepvel(:)           !dry deposition velocity (meters/day) from NCRP Commentary Table L4
c
c ============================================================================
c SECTION 6: DYNAMIC ARRAY DECLARATIONS FOR DOSE ACCUMULATION BY PATHWAY
c ============================================================================
c These arrays accumulate doses from different exposure pathways
c
        REAL*8, ALLOCATABLE :: imdose(:)              !immersion dose array
        REAL*8, ALLOCATABLE :: ihdose(:)              !inhalation dose array
        REAL*8, ALLOCATABLE :: surfdose(:)            !surface dose array
        REAL*8, ALLOCATABLE :: totimdose(:,:)         !total immersion dose by stack/radionuclide
        REAL*8, ALLOCATABLE :: totihdose(:,:)         !total inhalation dose by stack/radionuclide
        REAL*8, ALLOCATABLE :: totsurfdose(:,:)       !total surface dose by stack/radionuclide
        REAL*8, ALLOCATABLE :: totingdose(:,:)        !total ingestion dose by stack/radionuclide
c
c ============================================================================
c SECTION 7: DYNAMIC ARRAY DECLARATIONS FOR INGESTION DOSE RESULTS
c ============================================================================
c
        REAL*8, ALLOCATABLE :: milkdose(:,:)          !radionuclide dose by milk pathway (mrem/yr)
        REAL*8, ALLOCATABLE :: vegdose(:,:)           !radionuclide dose by vegetable pathway (mrem/yr)
        REAL*8, ALLOCATABLE :: meatdose(:,:)          !radionuclide dose by meat pathway (mrem/yr)
c
c ============================================================================
c SECTION 8: DYNAMIC ARRAY DECLARATIONS FOR WIND DATA BY SECTOR
c ============================================================================
c Wind data is binned by direction sectors (16-point compass typically)
c
        REAL*8, ALLOCATABLE :: windbinfreq(:,:)       !Sector frequency (proportion of time wind blows from that direction)
        REAL*8, ALLOCATABLE :: windbinspeed(:,:)      !Sector average wind speed (m/s)
        REAL*8, ALLOCATABLE :: windbinsrdist(:,:)     !Sector source-receptor distance (m)
c
c ============================================================================
c SECTION 9: DYNAMIC ARRAY DECLARATIONS FOR METADATA
c ============================================================================
c
        INTEGER, ALLOCATABLE :: nrads(:)              !number of radionuclides emitted at each stack
        INTEGER, ALLOCATABLE :: nwindbins(:)          !number of wind direction bins at each stack
c
        CHARACTER*8,  ALLOCATABLE :: radname(:,:)     !radionuclide names at each stack
        CHARACTER*3,  ALLOCATABLE :: windbinname(:,:) !wind direction sector abbreviation (N, NNE, etc)
        CHARACTER*3,  ALLOCATABLE :: sectorname(:,:)  !sector name for receptors
        CHARACTER*3,  ALLOCATABLE :: vegsector(:)     !wind sector where vegetables are located
        CHARACTER*3,  ALLOCATABLE :: milksector(:)    !wind sector where cow/milk is located
        CHARACTER*3,  ALLOCATABLE :: meatsector(:)    !wind sector where meat animals are located
        CHARACTER*8,  ALLOCATABLE :: decaychain(:)    !list of radionuclides in decay chain
        CHARACTER*8,  ALLOCATABLE :: decaychainname(:)!names of decay chain radionuclides
        CHARACTER*3,  ALLOCATABLE :: verticalmomentum(:) !flag for vertical momentum (YES or NO)
c
        CHARACTER*200,ALLOCATABLE :: windsource(:)    !source of wind data (e.g., weather station name)
        CHARACTER*200,ALLOCATABLE :: winddates(:)     !date coverage of wind data
        CHARACTER*200,ALLOCATABLE :: windloc(:)       !location of wind data collection
        CHARACTER*200,ALLOCATABLE :: winddist(:)      !distance from wind data to source
c
c ============================================================================
c SECTION 10: CONSTANTS AND PARAMETERS
c ============================================================================
c Maximum array sizes for decay chains and nuclide database
c
        INTEGER, PARAMETER :: MAX_CHAIN = 19                !max length of decay chain
        INTEGER, PARAMETER :: MAX_NUCLIDES = 826            !max number of nuclides in database
        INTEGER, PARAMETER :: MAX_DECAYCHAIN_INPUT = 879    !max decay chains in input file
c
c Physical and biological constants for dose calculations
c
        REAL, PARAMETER :: BREATH_RATE = 8000.0             !breathing rate (m3/yr)
        REAL, PARAMETER :: CONVEG = 70.0                    !vegetable consumption (kg/yr)
        REAL, PARAMETER :: CONMILK = 160.0                  !milk consumption (liters/yr)
        REAL, PARAMETER :: CONMEAT = 75.0                   !meat consumption (kg/yr)
        REAL, PARAMETER :: REMCON = 1.5000000E-02           !removal constant (1/yr) - weathering removal
        REAL, PARAMETER :: DEPTIM = 100.0 * 365.0           !deposition time (days)
        REAL, PARAMETER :: TWEATH = 0.693 / 12.0            !weathering half-life (months)
        REAL, PARAMETER :: TEVEG = 60.0                     !vegetable transit time from deposition to consumption (days)
        REAL, PARAMETER :: TEMILK = 30.0                    !milk transit time from animal feed to consumption (days)
        REAL, PARAMETER :: TEMEAT = 30.0                    !meat transit time from animal feed to consumption (days)
        REAL, PARAMETER :: DELVEG = 11.0                    !vegetable delivery time (days)
        REAL, PARAMETER :: DELMILK = 2.0                    !milk delivery time (days)
        REAL, PARAMETER :: DELMEAT = 17.0                   !meat delivery time (days)
        REAL, PARAMETER :: PI = 3.14159265                  !mathematical constant
c
c ============================================================================
c SECTION 11: LOCAL VARIABLES FOR CALCULATIONS
c ============================================================================
c
        REAL FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF,FACTORG,FACTORH,FACTORI
        REAL FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN
        REAL area                                      !cross-sectional area for plume dispersion
        REAL xcrit                                     !critical distance for transition between near/far field
        REAL SIGMAZ                                   !vertical standard deviation of plume (m)
        REAL CAPSIG                                   !combined sigma including building effects
        REAL pp                                       !Gaussian plume concentration multiplier
        REAL alldose                                  !cumulative dose from all pathways and radionuclides
        REAL idose                                    !intermediate dose variable
        REAL dist                                     !distance from source to receptor (m)
        REAL altconc                                  !alternative concentration calculation (QA check)
        REAL tempdose                                 !temporary dose variable for comparisons
        REAL maxdose                                  !maximum dose found in sector/radionuclide loop
        REAL maxdist                                  !distance at which maximum dose occurs
        REAL totmilkdose                              !total milk ingestion dose
        REAL totvegdose                               !total vegetable ingestion dose
        REAL itotvegdose                              !intermediate total vegetable dose
        REAL totmeatdose                              !total meat ingestion dose
        REAL itotcowdose                              !intermediate total cow/milk dose
        REAL capp                                     !Pasquill parameter
        REAL capn0                                    !initial condition parameter
        REAL outprod                                  !output product variable
        REAL inprod                                   !input product variable
        REAL z1                                       !height parameter
        REAL term1                                    !term 1 of momentum rise equation
        REAL term2                                    !term 2 of momentum rise equation
        REAL sum                                      !summation variable
        REAL transittime                              !transit time of air from emission release to receptor (days)
        REAL fdux                                     !flux term for deposition (units: s/m2)
        REAL wetdepadj                                !wet deposition adjustment factor
        REAL wetdepadjbad                             !alternative wet deposition adjustment (legacy)
        REAL xtime                                    !time variable for decay calculations
        REAL xremove                                  !removal rate for environmental transfers
        REAL ingestion                                !ingestion pathway dose component
        REAL immersion                                !immersion pathway dose component
        REAL inhalation                               !inhalation pathway dose component
        REAL surface                                  !surface contamination dose component
        REAL RHOAIR                                   !air density (kg/m3)
        REAL RHOSTK                                   !stack gas density (kg/m3)
        REAL STKVEL                                   !stack velocity (m/s)
        REAL R                                        !ratio parameter for plume rise
        REAL FM                                       !momentum flux parameter
        REAL DELHMO                                   !plume rise due to momentum (m)
        REAL DELHBU                                   !plume rise due to buoyancy (m)
        REAL HEFF                                     !effective release height accounting for plume rise (m)
        REAL FB                                       !buoyancy flux parameter
        REAL zz                                       !vertical Gaussian factor
        REAL INHFAC                                   !inhalation factor
        REAL TINIT                                    !initial temperature variable
        REAL T                                        !temperature variable
c
c ============================================================================
c SECTION 12: LOOP COUNTERS AND METADATA VARIABLES
c ============================================================================
c
        INTEGER nstacks                               !number of unique stacks
        INTEGER receptor                              !loop counter for receptors
        INTEGER N, M, I, J, K, L, II, JJ, P, Q        !loop counters (N=stacks, M=rads, J=wind sectors, L=chain)
        INTEGER nchain                                !number of nuclides in the current decay series
        INTEGER maxsect                               !sector number of maximum dose location
c
        CHARACTER*200 title                           !project title
        CHARACTER*200 company                         !company name
        CHARACTER*200 facility                        !facility name
        CHARACTER*200 address                         !facility address
        CHARACTER*200 poc                             !point of contact name
        CHARACTER*200 phone                           !phone number
        CHARACTER*200 inline                          !buffer for reading input lines
        CHARACTER*8 NUC                               !current nuclide name
        CHARACTER*2 NN                                !nuclide prefix for pathway selection (e.g., H-, C-, I-)
c
        integer date_time(8)                          !system date/time array from intrinsic function
        character*10 b(3)                             !character date/time strings
c
c ============================================================================
c MAIN PROGRAM BEGIN
c ============================================================================
c
        write(*,*)
        write(*,*)'COMPLY v2.0 (Optimized)'
        write(*,*)
c
c Open input and output files
c File 11 = Main input parameters
c File 12 = Output report
c File 13 = Diagnostic output
c File 9,10 = Reference databases
c
        open(11,file='input.level4.txt',status='old',form='formatted')
        open(12,file='output.level4.txt',status='unknown',form='formatted')
        open(13,file='diag.level4.txt',status='unknown',form='formatted')
c
        open(9,file='comply171_original_input_data.csv',status='old',form='formatted')
        open(10,file='comply171_original_decaychain_data.csv',status='old',form='formatted')
c
c ============================================================================
c READ INPUT NAMELIST HEADER INFORMATION
c ============================================================================
c
         read(11,'(20x,a)') title
         read(11,'(20x,a)') company
         read(11,'(20x,a)') facility
         read(11,'(20x,a)') address
         read(11,'(20x,a)') poc
         read(11,'(20x,a)') phone
c
c Read number of stacks - this is a key parameter that determines array sizes
c
         read(11,'(20x,i3)') nstacks
c
c ============================================================================
c ALLOCATE ARRAYS THAT DEPEND ONLY ON NSTACKS
c ============================================================================
c These arrays have one dimension per stack and must be allocated after
c reading nstacks but before the main loop
c
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
         allocate(nwindbins(nstacks))
         allocate(windsource(nstacks))
         allocate(winddates(nstacks))
         allocate(windloc(nstacks))
         allocate(winddist(nstacks))
         allocate(pctcalm(nstacks))
         allocate(rainfall(nstacks))
c
c ============================================================================
c ALLOCATE FIXED-SIZE ARRAYS BEFORE STACK LOOP (BUG FIX)
c ============================================================================
c These arrays have constant size (MAX_CHAIN) and are used in calculations
c within the stack loop. Allocating once outside the loop improves performance
c and prevents memory leaks from repeated allocation/deallocation
c
         allocate(halflife(MAX_CHAIN))
         allocate(immersionfac(MAX_CHAIN))
         allocate(inhalationfac(MAX_CHAIN))
         allocate(ingestionfac(MAX_CHAIN))
         allocate(surfacefac(MAX_CHAIN))
         allocate(xbveg(MAX_CHAIN))
         allocate(xbfor(MAX_CHAIN))
         allocate(xfmilk(MAX_CHAIN))
         allocate(xfmeat(MAX_CHAIN))
         allocate(wbfac(MAX_CHAIN))
         allocate(xwbfac(MAX_CHAIN))
         allocate(vegfac(MAX_CHAIN))
         allocate(cowfac(MAX_CHAIN))
         allocate(decaychain(MAX_CHAIN))
         allocate(drydepvel(MAX_CHAIN))
         allocate(decaychainname(MAX_CHAIN))
         allocate(sourceterm(MAX_CHAIN))
         allocate(initialamt(MAX_CHAIN))
         allocate(activity(MAX_CHAIN))
         allocate(surfaceactivity(MAX_CHAIN))
         allocate(ambientactivity(MAX_CHAIN))
         allocate(normsourceactivity(MAX_CHAIN))
         allocate(kvalue(MAX_CHAIN))
         allocate(ktot(MAX_CHAIN))
         allocate(imdose(MAX_CHAIN))
         allocate(ihdose(MAX_CHAIN))
         allocate(surfdose(MAX_CHAIN))
         allocate(aveglf(MAX_CHAIN))
         allocate(avegrt(MAX_CHAIN))
         allocate(amilklf(MAX_CHAIN))
         allocate(amilkrt(MAX_CHAIN))
         allocate(ameatlf(MAX_CHAIN))
         allocate(ameatrt(MAX_CHAIN))
         allocate(cfveg(MAX_CHAIN))
         allocate(cfmilk(MAX_CHAIN))
         allocate(cfmeat(MAX_CHAIN))
         allocate(dveg(MAX_CHAIN))
         allocate(dmilk(MAX_CHAIN))
         allocate(dmeat(MAX_CHAIN))
         allocate(ingdose(MAX_CHAIN))
c
c ============================================================================
c MAIN LOOP OVER STACKS
c ============================================================================
c Read detailed parameters for each stack and all its radionuclide emissions
c
         do n = 1, nstacks 
c
         read(11,'(a)') !blank spacer to improve readability of the input file
c
c Read stack location and physical parameters
c
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
c
c Read receptor location parameters for this stack
c (distances and sectors for special receptors like gardens, dairy farms)
c
         read(11,'(20x,f10.5)') distance(n)
         read(11,'(20x,f10.5)') distanceveg(n)
         read(11,'(20x,f10.5)') distancecow(n)
         read(11,'(20x,f10.5)') distancemeat(n)
c
         read(11,'(20x,a)') vegsector(n)
         read(11,'(20x,a)') milksector(n)
         read(11,'(20x,a)') meatsector(n)
c
         read(11,'(20x,f10.5)') rainfall(n)
         read(11,'(20x,i3)') nrads(n)
c
c ============================================================================
c ALLOCATE STACK AND RADIONUCLIDE-DEPENDENT ARRAYS
c ============================================================================
c These 2D arrays depend on both number of stacks and radionuclides per stack.
c They are allocated once on first stack iteration after reading nrads data.
c This is done outside the innermost loop to avoid repeated allocations.
c
         if(n .eq. 1) then
           allocate(radname(nstacks,maxval(nrads)))
           allocate(rademis(nstacks,maxval(nrads)))
           allocate(conc(nstacks,maxval(nrads)))
           allocate(dose(nstacks,maxval(nrads)))
           allocate(totimdose(nstacks,maxval(nrads)))
           allocate(totihdose(nstacks,maxval(nrads)))
           allocate(totsurfdose(nstacks,maxval(nrads)))
           allocate(totingdose(nstacks,maxval(nrads)))
           allocate(milkdose(nstacks,maxval(nrads)))
           allocate(vegdose(nstacks,maxval(nrads)))
           allocate(meatdose(nstacks,maxval(nrads)))
         endif
c
c ============================================================================
c ALLOCATE WIND SECTOR-DEPENDENT ARRAYS
c ============================================================================
c Wind data arrays depend on number of wind bins per stack (typically 16)
c These can vary by stack so allocation must occur per stack
c
         allocate(windbinname(nstacks,nwindbins(n)))
         allocate(windbinfreq(nstacks,nwindbins(n)))
         allocate(windbinspeed(nstacks,nwindbins(n)))
         allocate(windbinsrdist(nstacks,nwindbins(n)))
         allocate(sectorname(nstacks,nwindbins(n)))
c
c ============================================================================
c READ RADIONUCLIDE EMISSIONS FOR THIS STACK
c ============================================================================
c
         do m = 1, nrads(n)
         read(11,'(20x,a)') inline
         read(inline,*) radname(n,m),rademis(n,m)
         enddo !end loop over reading in radionuclide emissions
c
c ============================================================================
c READ WIND DATA AND SOURCE-RECEPTOR DISTANCES
c ============================================================================
c
         read(11,'(20x,a)') windsource(n)
         read(11,'(20x,a)') winddates(n)
         read(11,'(20x,a)') windloc(n)
         read(11,'(20x,a)') winddist(n)
c
         read(11,'(20x,f10.5)') pctcalm(n)
         read(11,'(20x,i3)') nwindbins(n)
c
         do m = 1, nwindbins(n)
         read(11,'(20x,a)') inline
         read(inline,*) windbinname(n,m),windbinfreq(n,m),
     &                  windbinspeed(n,m),windbinsrdist(n,m),sectorname(n,m)
         enddo !end loop over reading in wind data and source-receptor distances
c
         enddo !end loop over reading in stacks
c
c ============================================================================
c BEGIN DOSE CALCULATION SECTION
c ============================================================================
c Emission rate assumption is Ci/s for all calculations below
c
         write(13,3001) !write header to file
c
c ============================================================================
c LOOP OVER ALL STACKS FOR DOSE CALCULATION
c ============================================================================
c
         do n = 1, nstacks
c
c Initialize dose accumulators for this stack
c
         alldose = 0.
         maxsect = 0
         maxdist = 0.0
c
c ============================================================================
c HANDLE SURFACE RELEASES (HEIGHT = 0)
c ============================================================================
c For ground-level releases, building dimensions are set to zero
c (no additional plume rise from building effects)
c
         if(releaseht(n).eq.0.0) then !surface release
         bldwidth(n) = 0.0
         bldheight(n) = 0.0
         endif
c
c Calculate building cross-sectional area and critical distance
c The critical distance marks the transition between near-field and far-field
c dispersion regimes. Beyond this distance, wind speed variation becomes important.
c
         area = bldwidth(n) * bldheight(n)
         xcrit = 2.5 * sqrt(area)
c
c ============================================================================
c LOOP OVER RADIONUCLIDES AT THIS STACK
c ============================================================================
c
         do m = 1, nrads(n)
c
c Initialize dose components for this radionuclide
c
         maxdose = 0.0
         inhalation = 0.
         immersion = 0.
         ingestion = 0.
         totmilkdose = 0.
         totvegdose = 0.
         totmeatdose = 0.
         surface = 0.
c
c ============================================================================
c LOOP OVER WIND DIRECTION SECTORS
c ============================================================================
c Wind data is typically divided into 16 sectors (22.5 degrees each)
c Calculations are performed for each sector independently
c
         do j = 1, nwindbins(n)
c
c ============================================================================
c LOOP OVER DISTANCES (RECEPTORS)
c ============================================================================
c Receptors are evaluated at distances from 10m to 3010m in 10m increments
c This provides a detailed dose profile with distance
c
         do receptor = 10, 3010, 10
c
         dist = float(receptor)
c
c ============================================================================
c ATMOSPHERIC DISPERSION AND PLUME RISE CALCULATIONS
c ============================================================================
c This section calculates the air concentration using the Gaussian plume model,
c accounting for plume rise due to buoyancy and momentum
c
         if(exitvelocity(n).le.0.) then
         exitvelocity(n) = 0.3  !Set minimum exit velocity if not specified
         endif
c
         HEFF = releaseht(n)
c
c Check if stack release height exceeds 2.5 times building height
c If so, calculate plume rise; otherwise use release height as-is
c
         IF ( releaseht(n) .GT. 2.5*bldheight(n) ) THEN
c
c ============================================================================
c CALCULATE ATMOSPHERIC PROPERTIES AND PLUME RISE
c ============================================================================
c
c Calculate air and stack gas densities based on temperature
c RHOAIR: Reference density at 80F converted to kg/m^3
c RHOSTK: Stack gas density accounting for temperature difference
c
          RHOAIR = 0.0735 * ( 80.0 + 460.0 ) / ( airtemp(n) + 460.0 )
          RHOAIR = RHOAIR * 454.0 / ( 0.3048 ** 3 )
          RHOSTK = RHOAIR * ( airtemp(n) + 460.0 ) / ( stacktemp(n) + 460.0 )
c
c ============================================================================
c CALCULATE MOMENTUM RISE PARAMETERS
c ============================================================================
c
c Stack velocity is calculated from exit velocity and area
c R is the ratio of stack velocity to wind speed
c FM is the momentum flux parameter
c TERM1 and TERM2 are components of the momentum rise equation
c
          AREA   = PI * ( diameter(n) ** 2.0 ) / 4.0
          STKVEL = exitvelocity(n) / AREA
          R      = STKVEL / windbinspeed(n,j)
          FM     = ( RHOSTK / RHOAIR ) * STKVEL * exitvelocity(n) / PI
          TERM1  = 75.0 * FM / ( ( ( 1.0 + 3.0 / R ) * windbinspeed(n,j) ) ** 2.0 )
          TERM2  = releaseht(n) + 7.79 * SQRT( FM ) / ( ( 1.0 + 3.0 / R ) * windbinspeed(n,j) )
c
c ============================================================================
c APPLY VERTICAL MOMENTUM FLAG
c ============================================================================
c If the stack has vertical momentum (YES), apply momentum rise calculation
c Otherwise, set momentum rise to zero (stack blows horizontally into wind)
c
          IF ( verticalmomentum(n).eq.'YES' ) THEN
           DELHMO = 0.93 * ( TERM1 ** ( 3.0 / 7.0 ) ) *
      &        ( TERM2 ** ( 1.0 / 7.0 ) )
          ELSE
           DELHMO = 0.0
          ENDIF
c
c ============================================================================
c CALCULATE BUOYANT RISE
c ============================================================================
c Buoyant rise occurs because hot stack gas is less dense than surrounding air
c FB is the buoyancy flux parameter
c DELHBU is the buoyant rise contribution
c
          FB = 9.81 * ( RHOAIR - RHOSTK ) * exitvelocity(n) / ( PI * RHOAIR )
          DELHBU = 23.67 * ( ( FB / ( windbinspeed(n,j) ** 3.0 ) ) ** 0.6 ) *
      &           ( releaseht(n) + 195.1 * FB / ( windbinspeed(n,j) ** 3.0 ) ) ** 0.4
c
c ============================================================================
c EFFECTIVE RELEASE HEIGHT
c ============================================================================
c Use the maximum of buoyant or momentum rise to get effective height.
c The effective height is used in Gaussian plume concentration calculations.
c
         HEFF = releaseht(n) + AMAX1( DELHMO, DELHBU )
c
         ENDIF !end condition where H > 0
c
c ============================================================================
c GAUSSIAN PLUME DISPERSION PARAMETERS
c ============================================================================
c
c SIGMAZ is the vertical standard deviation of the concentration distribution
c It increases with distance from the source (plume spreading)
c CAPSIG combines SIGMAZ with building effects via CAPSIG = sqrt(SIGMAZ^2 + area/pi)
c
         SIGMAZ = (0.06 * dist) / SQRT( 1.0 + 0.0015 * dist)
         CAPSIG = SQRT( SIGMAZ ** 2 + area / PI )
         zz = 0.5 * ( ( HEFF / CAPSIG ) ** 2 )
c
c ============================================================================
c GAUSSIAN PLUME CONCENTRATION FACTOR
c ============================================================================
c pp is the vertical Gaussian factor. For large zz, this becomes
c proportional to exp(-zz). For small zz, approximate as constant 2.032/(dist*CAPSIG)
c
          if(zz.gt.0.001.and.releaseht(n).gt.2.5*bldheight(n)) then
          pp = 2.032 * EXP( -zz ) / ( dist * CAPSIG )
          else
          pp = 2.032 / ( dist * CAPSIG )
          endif
c
c ============================================================================
c ATMOSPHERIC TRANSIT TIME AND DEPOSITION PARAMETERS
c ============================================================================
c
c transittime: Time for air parcel to travel from source to receptor (days)
c fdux: Flux parameter for wet/dry deposition calculations (s/m2)
c
         transittime = (dist / windbinspeed(n,j))*(1./3600.)*(1./24.) ! days
         fdux = (windbinfreq(n,j)/(1.0-pctcalm(n))) / (windbinspeed(n,j)*dist) !s/m2 units
c
c ============================================================================
c AIR CONCENTRATION CALCULATION
c ============================================================================
c Two regimes are used depending on whether we're in near-field or far-field:
c 1. Near-field (dist < xcrit): Use building dimension as mixing height
c 2. Far-field (dist >= xcrit): Use Gaussian plume equation with pp factor
c
c Correction for calm days is applied via (1.0 - pctcalm)
c
         if(dist.lt.xcrit) then
         conc(n,m) = ( (windbinfreq(n,j)/(1.0-pctcalm(n))) *rademis(n,m)) /
      &          (PI*windbinspeed(n,j)*1.0*amin1(bldwidth(n),bldheight(n)))
         else
         conc(n,m) = ((windbinfreq(n,j)/(1.0-pctcalm(n))) *rademis(n,m)*pp) / windbinspeed(n,j)
         endif
c
c Set minimum concentration threshold to avoid numerical issues
c
         if(conc(n,m).lt. 1.0E-32 ) then
           conc(n,m) = 1.0E-32
         endif
c
c Calculate alternative concentration for quality assurance checking
c
         altconc = (windbinfreq(n,j)*rademis(n,m)*pp) / windbinspeed(n,j) !QA only
c
c ============================================================================
c WET DEPOSITION ADJUSTMENT FACTOR
c ============================================================================
c This factor accounts for removal of radionuclides from air to ground via
c precipitation. The adjustment depends on rainfall rate and other parameters.
c
         wetdepadj = (8.0 * fdux * rainfall(n) * 1.0E-05 * 
     &                (rademis(n,m) * 3600. * 24.))/(PI*conc(n,m))
c
c ============================================================================
c READ DECAY CHAIN AND DOSE FACTOR DATA
c ============================================================================
c This section looks up the decay chain for the current radionuclide
c and reads all associated dose factors from the reference databases
c
         rewind(9)
         rewind(10)
         read(9,'(a)') inline
c
           do ii = 1, MAX_DECAYCHAIN_INPUT
           read(10,'(a)') inline
           read(inline,*) NUC,k,(decaychain(l),l=1,k)
           if(NUC.eq.radname(n,m)) then
            read(inline,*) NUC,k,(decaychain(l),l=1,k)
            do jj = 1, MAX_NUCLIDES
            read(9,'(a)') inline
            read(inline,*) NUC,FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF,FACTORG,
      &      FACTORH,FACTORI,FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN
             do l = 1, k
             if(NUC.eq.decaychain(l)) then
              decaychainname(l) = NUC
              halflife(l) = FACTOR*365.               !Convert half-life to days
              wbfac(l) = FACTORB                       !Whole body dose factor
              xwbfac(l) = FACTORC                      !External whole body dose factor
              vegfac(l) = FACTORD                      !Vegetable dose factor
              cowfac(l) = FACTORE                      !Cow dose factor
              inhalationfac(l) = FACTORF              !Inhalation dose factor
              immersionfac(l) = FACTORG               !Immersion dose factor
              surfacefac(l) = FACTORH                 !Surface dose factor
              ingestionfac(l) = FACTORI               !Ingestion dose factor
              xbveg(l) = FACTORJ                      !Bioaccumulation veg
              xbfor(l) = FACTORK                      !Bioaccumulation forage
              xfmilk(l) = FACTORL                     !Transfer factor milk
              xfmeat(l) = FACTORM                     !Transfer factor meat
              drydepvel(l) = FACTORN                  !Dry deposition velocity (m/day)
             endif
             enddo !loop over decay chain
            enddo !loop over all nuclides on input reference table 
            nchain = k                               !Store number of nuclides in chain
           endif 
           enddo !loop over all nuclides on decay chain input table
c
c ============================================================================
c CALCULATE IMMERSION AND INHALATION DOSES
c ============================================================================
c These doses result from direct exposure to radionuclides in air
c
           call init_array(activity, nchain)
           call init_array(sourceterm, nchain)
           call init_array(initialamt, nchain)
           
           initialamt(1) = conc(n,m)
           xremove = 0.0
           xtime = transittime
c
c Call decay chain subroutine to calculate activity at receptor location
c
        CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
           totimdose(n,m) = 0.0
           totihdose(n,m) = 0.0
c
c Sum doses from all decay chain members (parent + progeny)
c
           do l = 1, nchain
           ambientactivity(l) = activity(l)
           imdose(l) = ambientactivity(l) * immersionfac(l)
           ihdose(l) = ambientactivity(l) * inhalationfac(l) * BREATH_RATE
           totimdose(n,m) = totimdose(n,m) + imdose(l)
           totihdose(n,m) = totihdose(n,m) + ihdose(l)
           enddo
c
c ============================================================================
c CALCULATE SURFACE CONTAMINATION DOSE
c ============================================================================
c This dose results from radionuclides deposited on ground and building surfaces
c where people spend time.
c
           call init_array(sourceterm, nchain)
           call init_array(initialamt, nchain)
           call init_array(activity, nchain)
c
           do l = 1, nchain
           sourceterm(l) = ambientactivity(l)
           enddo
c
           xremove = REMCON / 365.
           xtime = DEPTIM
c
        CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
           totsurfdose(n,m) = 0.0
c
           do l = 1, nchain
c
c Noble gases (He, Ne, Ar, Kr, Xe, Rn) do not deposit, so only use dry deposition
c All other nuclides deposit via both dry and wet processes
c
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
c
c ============================================================================
c CALCULATE INGESTION PATHWAY DOSES
c ============================================================================
c
c Step 1: Calculate normalized surface activity for ingestion pathway
c This is the baseline activity from which crop and animal uptake is calculated
c
           call init_array(sourceterm, nchain)
           call init_array(initialamt, nchain)
           call init_array(activity, nchain)
c
           initialamt(1) = 1.0
           xremove = 0.0
           xtime = transittime
c
        CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
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
c
c ============================================================================
c VEGETABLES PATHWAY CALCULATION
c ============================================================================
c
c Step 2a: Calculate leaf uptake of deposited activity
c Leaves directly intercept deposited material and retain it during growing season
c
           call init_array(sourceterm, nchain)
           call init_array(initialamt, nchain)
           call init_array(activity, nchain)
c
           do l = 1, nchain
           sourceterm(l) = normsourceactivity(l)
           enddo
c
           xremove = TWEATH
           xtime = TEVEG
c
        CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
           do l = 1, nchain
           aveglf(l) = activity(l) 
           enddo
c
c Step 2b: Calculate root uptake from soil residue (long-term deposition)
c Roots slowly take up activity from contaminated soil
c
           call init_array(sourceterm, nchain)
           call init_array(initialamt, nchain)
           call init_array(activity, nchain)
c
           do l = 1, nchain
           sourceterm(l) = normsourceactivity(l)
           enddo
c
           xremove = REMCON / 365.
           xtime = DEPTIM
c
        CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
           do l = 1, nchain
           avegrt(l) = activity(l)
           enddo
c
c ============================================================================
c MILK PATHWAY CALCULATION
c ============================================================================
c
c Step 3a: Calculate activity in forage from animal consumption of contaminated plants
c
           call init_array(sourceterm, nchain)
           call init_array(initialamt, nchain)
           call init_array(activity, nchain)
c
           do l = 1, nchain
           sourceterm(l) = normsourceactivity(l)
           enddo
c
           xremove = TWEATH
           xtime = TEMILK
c
        CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
           do l = 1, nchain
           amilklf(l) = activity(l)
           enddo
c
c Step 3b: Calculate long-term residue uptake by animals
c
           call init_array(sourceterm, nchain)
           call init_array(initialamt, nchain)
           call init_array(activity, nchain)
c
           do l = 1, nchain
           sourceterm(l) = normsourceactivity(l)
           enddo
c
           xremove = REMCON / 365.
           xtime = DEPTIM
c
        CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
           do l = 1, nchain
           amilkrt(l) = activity(l)
           enddo
c
c ============================================================================
c MEAT PATHWAY CALCULATION
c ============================================================================
c
c Step 4a: Calculate activity in animal diet from contaminated forage
c
           call init_array(sourceterm, nchain)
           call init_array(initialamt, nchain)
           call init_array(activity, nchain)
c
           do l = 1, nchain
           sourceterm(l) = normsourceactivity(l)
           enddo
c
           xremove = TWEATH
           xtime = TEMEAT
c
        CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
           do l = 1, nchain
           ameatlf(l) = activity(l)
           enddo
c
c Step 4b: Calculate residual activity in meat animals from long-term pasture contamination
c
           call init_array(sourceterm, nchain)
           call init_array(initialamt, nchain)
           call init_array(activity, nchain)
c
           do l = 1, nchain
           sourceterm(l) = normsourceactivity(l)
           enddo
   
           xremove = REMCON / 365.
           xtime = DEPTIM
c
        CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
           do l = 1, nchain
           ameatrt(l) = activity(l)
           enddo
c
c ============================================================================
c ESTIMATE CONCENTRATION IN FOOD MEDIA
c ============================================================================
c This section uses known relationships between air/soil contamination
c and contamination levels in food products
c
           read(radname(n,m),'(a2)') NN
c
           do l = 1, nchain
c
c Special cases for tritium (H-3) and carbon-14 (C-14)
c These isotopes have unique biogeochemical behavior
c
            if(NN.eq.'H-') then 
c H-3 (tritium) - D. A. Baker, et al, EPA 600/9-76-016, July 1976.
c Tritium is highly mobile and concentrates in water-containing food products
c
            cfveg(l) = conc(n,m)*9.0*(0.1/0.8E-02)
            cfmilk(l) = (conc(n,m)*9.0*7.0E-02/8.0E-03)*16.0*0.11/(7.0E-02*16. + 220.0/9.0)
            cfmeat(l) = (conc(n,m)*9.0*7.0E-02/8.0E-03)*12.0*0.10/(7.0E-02*12. + 50.0/9.0) 
c
            elseif(NN.eq.'C-') then 
c C-14 - D. A. Baker, et al, EPA 600/9-76-016, July 1976.
c Carbon-14 is incorporated into organic matter throughout the food chain
c
            cfveg(l) = conc(n,m)*(9.000000E-02/0.16E-03)
            cfmilk(l) = conc(n,m)*(7.0000000E-02/0.16E-03)
            cfmeat(l) = conc(n,m)*(0.2400000/0.16E-03)
c
            else 
c General case: Equation B2 of NCRP Commentary 3
c Uses weighted combination of leaf and root uptake pathways
c
            cfveg(l) =  conc(n,m) * ((0.1/1.0)*aveglf(l)+(xbveg(l)/220.)*avegrt(l))
            cfmilk(l) = conc(n,m) * ((1.8/1.0)*amilklf(l)+(xbfor(l)/220.)*amilkrt(l))*xfmilk(l)*16.0
            cfmeat(l) = conc(n,m) * ((1.8/1.0)*ameatlf(l)+(xbfor(l)/220.)*ameatrt(l))*xfmeat(l)*12.0
            endif
c
            enddo
c
c ============================================================================
c VEGETABLES INGESTION DOSE CALCULATION
c ============================================================================
c
             call init_array(sourceterm, nchain)
             call init_array(initialamt, nchain)
             call init_array(activity, nchain)
c
             do l = 1, nchain
             sourceterm(l) = 0.0
             initialamt(l) = cfveg(l)
             enddo
c
             xremove = 0.0
             xtime = DELVEG                          !Days from harvest to consumption
c
             CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
              do l = 1, nchain
              dveg(l) = activity(l) * ingestionfac(l) * CONVEG
              enddo
c
c ============================================================================
c MILK INGESTION DOSE CALCULATION
c ============================================================================
c
             call init_array(sourceterm, nchain)
             call init_array(initialamt, nchain)
             call init_array(activity, nchain)
c
             do l = 1, nchain
             sourceterm(l) = 0.0
             initialamt(l) = cfmilk(l)
             enddo
c
             xremove = 0.0
             xtime = DELMILK                         !Days from production to consumption
c
             CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
              do l = 1, nchain
              dmilk(l) = activity(l) * ingestionfac(l) * CONMILK
              enddo
c
c ============================================================================
c MEAT INGESTION DOSE CALCULATION
c ============================================================================
c
             call init_array(sourceterm, nchain)
             call init_array(initialamt, nchain)
             call init_array(activity, nchain)
c
             do l = 1, nchain
             sourceterm(l) = 0.0
             initialamt(l) = cfmeat(l)
             enddo
c
             xremove = 0.0
             xtime = DELMEAT                         !Days from slaughter to consumption
c
             CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)
c
              do l = 1, nchain
              dmeat(l) = activity(l) * ingestionfac(l) * CONMEAT
              enddo
c
c ============================================================================
c TOTAL INGESTION DOSE BY PATHWAY
c ============================================================================
c Sum doses from all decay chain members for each food pathway
c
            vegdose(n,m) = 0.0
            milkdose(n,m) = 0.0
            meatdose(n,m) = 0.0
            do l = 1, nchain
            vegdose(n,m) = vegdose(n,m) + dveg(l)
            milkdose(n,m) = milkdose(n,m) + dmilk(l)
            meatdose(n,m) = meatdose(n,m) + dmeat(l)
            enddo
c
c ============================================================================
c SAVE DOSE INFORMATION AT KEY RECEPTOR LOCATIONS FOR FINAL REPORT
c ============================================================================
c
c Check if this distance and sector match key receptor locations.
c If so, save the dose for output reporting.
c
         if(dist.eq.distanceveg(n).and.vegsector(n).eq.sectorname(n,j)) then
          alldose = alldose + vegdose(n,m)
          totvegdose = vegdose(n,m)
          write(*,*)'Vegetable dose:',radname(n,m),vegdose(n,m),dist,sectorname(n,j)
          endif
         
         if(dist.eq.distancecow(n).and.milksector(n).eq.sectorname(n,j)) then
          alldose = alldose + milkdose(n,m)
          totmilkdose = milkdose(n,m)
          write(*,*)'Milk dose:',radname(n,m),milkdose(n,m),dist,sectorname(n,j)
         endif

         if(dist.eq.distancemeat(n).and.meatsector(n).eq.sectorname(n,j)) then
          alldose = alldose + meatdose(n,m)
          totmeatdose = meatdose(n,m)
          write(*,*)'Meat dose:',radname(n,m),meatdose(n,m),dist,sectorname(n,j)
         endif
   
c Check if this distance matches the wind bin source-receptor distance
c If so, this is the "critical receptor" for which we report maximum dose
c
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
c
c Write diagnostic output for all distances and sectors
c
         write(13,3000) n,radname(n,m),sectorname(n,j),dist,
      &     conc(n,m),totimdose(n,m),totihdose(n,m),totsurfdose(n,m),
      &     meatdose(n,m),milkdose(n,m),vegdose(n,m)
c
         enddo !end loop over distances between source and receptor
         enddo !end loop over wind direction bins
c
c ============================================================================
c OUTPUT DOSE SUMMARY FOR THIS RADIONUCLIDE
c ============================================================================
c
         alldose = alldose + surface + immersion + inhalation
         if(maxsect .gt. 0) then
           write(*,*)'Max dose (',maxdose,') receptor at ',
     &               windbinsrdist(n,maxsect),sectorname(n,maxsect)
         endif
         write(*,*)
         write(*,2999)'Inhalation dose:',radname(n,m),inhalation
         write(*,2999)'Immersion dose:',radname(n,m),immersion
         write(*,2999)'Surface dose:',radname(n,m),surface
         write(*,2999)'Ingestion meat dose:',radname(n,m),totmeatdose
         write(*,2999)'Ingestion milk dose:',radname(n,m),totmilkdose
         write(*,2999)'Ingestion vegetables dose:',radname(n,m),totvegdose
         write(*,2999)'Total ingestion dose:',radname(n,m),
     &                totmeatdose+totmilkdose+totvegdose
         write(*,2999)'Total dose:',radname(n,m),
     &     (inhalation+immersion+surface+totmeatdose+totmilkdose+totvegdose)
         write(*,*)'--------------------------------------------------------------------------' 
         enddo !end loop over radionuclides
c
         write(*,*)
         write(*,3100)'Total dose (mrem/yr):',alldose
         write(*,*)'Total dose limit is 10.0 mrem/yr'
         write(*,*)'Iodine limit is 3.0 mrem/yr'
         write(*,*)
c
         enddo !end loop over stacks
c
c ============================================================================
c FORMAT STATEMENTS FOR DIAGNOSTIC OUTPUT
c ============================================================================
c
   2999  FORMAT(a30,a10,f20.15)
   3000  FORMAT(i2,',',a,',',a,',',f6.1,',',f30.20,6(',',f20.10))
   3001  FORMAT('STACK,RADNAME,SECTOR,DISTANCE,CONC,IMMERSION,INHALATION,SURFACE,MEAT,MILK,VEG')
   3100  FORMAT(a,f10.4)
c
c ============================================================================
c WRITE STANDARDIZED COMPLIANCE OUTPUT REPORT
c ============================================================================
c
         write(12,4000)
         write(12,2001)
c
         write(12,1000) title
         write(12,1001) company
         write(12,1002) facility 
         write(12,1003) address
         write(12,1004) poc
         write(12,1005) phone
c
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
c
         write(12,2002)
         if(maxsect .gt. 0) then
           write(12,1800)'Max dose receptor ',windbinsrdist(n,maxsect),
     &                   ' m downwind in sector ',sectorname(n,maxsect)
         endif
   1800  FORMAT(a,f5.0,a,a)
         write(12,2002)
c
         write(12,1998)
         do m=1, nrads(n)
         write(12,1999) radname(n,m),rademis(n,m),conc(n,m),
      &     (totimdose(n,m)+totihdose(n,m)+totsurfdose(n,m)+
      &      meatdose(n,m)+milkdose(n,m)+vegdose(n,m)),
      &     totimdose(n,m),totihdose(n,m),totsurfdose(n,m),
      &     meatdose(n,m),milkdose(n,m),vegdose(n,m)
c
         enddo !end loop over radionuclides
         enddo !end loop over stacks
c
         write(12,2002)
         write(12,2001)
         write(12,2001)
         write(12,2002)
         write(12,2050) alldose
         write(12,2002)
c
c ============================================================================
c ADD TIMESTAMP TO OUTPUT FILE
c ============================================================================
c
         call date_and_time(b(1), b(2), b(3), date_time)
c
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
c
         write(12,2002)
         write(12,2001)
         write(12,2002)
         write(12,2100) date_time(2),date_time(3),date_time(1),
      &           date_time(5),date_time(6)
c
c ============================================================================
c FORMATTING FOR OUTPUT FILE
c ============================================================================
c
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
c
   1998 FORMAT('Radionuclide,Emissions,Concentration,Total Dose,Immersion,Inhalation,Surface,Meat,Milk,Vegetable')
   1999 FORMAT(a,',',E10.3,',',E10.3,7(',',E10.3))
c
   2001 FORMAT('-----------------------------------------------------------------------')
   2002 FORMAT(' ')
c
   2050 FORMAT('Total effective dose equivalent all stacks (mrem/yr):',E10.3)
   2100 FORMAT('Date/Time (LST) stamp of code execution (MM/DD/YYYY HH:MM):',i2,'/',i2,'/',i4,' ',i2,':',i2)
   4000 FORMAT('COMPLY version 2.0 level 4 distributed through GitHub')
c
c ============================================================================
c END OF MAIN PROGRAM
c ============================================================================
c
        END 
c
c ============================================================================
c HELPER SUBROUTINE: init_array
c ============================================================================
c
c PURPOSE:
c   Efficiently initialize an array to zero
c
c PARAMETERS:
c   arr(n) - Array to be initialized
c   n      - Number of elements in array
c
        SUBROUTINE init_array(arr, n)
        IMPLICIT NONE
        INTEGER n, i
        REAL*8 arr(n)
        do i = 1, n
          arr(i) = 0.0
        enddo
        END
