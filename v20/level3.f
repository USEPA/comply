       PROGRAM COMPLY
c
c==============================================================================
c PROGRAM: COMPLY v2.0 - EPA Radiation Dose Model (Level 3)
c
c PURPOSE:
c   Calculates effective radiation dose to human receptors from airborne 
c   radionuclide emissions. The model applies Gaussian plume dispersion 
c   theory to estimate ground-level concentrations at specified distances 
c   from emission stacks, then applies dose conversion factors to estimate
c   individual and collective dose impacts across multiple exposure pathways:
c   - Direct inhalation (whole body dose)
c   - Vegetable ingestion pathway
c   - Meat and milk ingestion pathway
c
c INPUT FILES:
c   - input.level3.txt: Facility and emission stack parameters
c   - comply171_original_input_data.csv: Radionuclide-specific dose factors
c
c OUTPUT FILES:
c   - output.level3.txt: Formatted dose calculation results
c   - diag.level3.txt: Detailed diagnostic output for all distances/radionuclides
c
c==============================================================================

      IMPLICIT NONE

c ============================================================================
c SECTION 1: DECLARE DYNAMIC ARRAYS
c ============================================================================
c
c Stack-specific input parameters (allocated by number of stacks)

       REAL*8, ALLOCATABLE :: releaseht(:)           !stack height (m)      
       REAL*8, ALLOCATABLE :: diameter(:)            !stack diameter (m)
       REAL*8, ALLOCATABLE :: flowrate(:)            !volumetric flow rate (m3/s)
       REAL*8, ALLOCATABLE :: airtemp(:)             !air temperature (C)
       REAL*8, ALLOCATABLE :: stacktemp(:)           !stack temperature (C)
       REAL*8, ALLOCATABLE :: bldheight(:)           !building height (m)
       REAL*8, ALLOCATABLE :: bldwidth(:)            !building width (m)

c Receptor location parameters (distance from stack to exposure point)

       REAL*8, ALLOCATABLE :: distance(:)        !distance between stack and receptor (m)
       REAL*8, ALLOCATABLE :: distanceveg(:)     !distance between garden and receptor (m)
       REAL*8, ALLOCATABLE :: distancecow(:)     !distance between cows and receptor (m)
       REAL*8, ALLOCATABLE :: windspeed(:)       !wind speed (m/s)
       REAL*8, ALLOCATABLE :: latitude(:)        !latitude (decimal degrees)
       REAL*8, ALLOCATABLE :: longitude(:)       !longitude (decimal degrees)

c Radionuclide emission and dose arrays (2D: nstacks x max_radionuclides)

       REAL*8, ALLOCATABLE :: rademis(:,:)       !radionuclide emissions (Ci/yr) for each stack
       REAL*8, ALLOCATABLE :: dose(:,:)          !radionuclide dose (mrem/yr)
       REAL*8, ALLOCATABLE :: wbfac(:,:)         !whole body dose factor from Oak Ridge table
       REAL*8, ALLOCATABLE :: vegfac(:,:)        !vegetable dose factor from Oak Ridge table
       REAL*8, ALLOCATABLE :: cowfac(:,:)        !cow dose factor from Oak Ridge table
       REAL*8, ALLOCATABLE :: vegdose(:,:)       !radionuclide dose by vegetables (mrem/yr)
       REAL*8, ALLOCATABLE :: cowdose(:,:)       !radionuclide dose by meat and milk (mrem/yr)

c Metadata arrays

       INTEGER, ALLOCATABLE :: nrads(:)          !number of emitted radionuclides at each stack

       CHARACTER*8, ALLOCATABLE :: radname(:,:)  !radionuclide names at each stack
       CHARACTER*3, ALLOCATABLE :: verticalmomentum(:)    !vertical momentum (YES or NO)

c
c ============================================================================
c SECTION 2: DECLARE WORKING VARIABLES AND CALCULATED TERMS
c ============================================================================
c
c Gaussian plume and dispersion parameters
c   xcrit = critical distance for near-field effects
c   SIGMAZ = vertical dispersion parameter
c   CAPSIG = effective cross-section parameter
c   pp = plume deposition/dispersion factor
c   zz = normalized plume height parameter

       REAL*8 area, xcrit, SIGMAZ, CAPSIG, pp, alldose, idose, dist, altconc, zz

c Concentration and dose working variables

       REAL*8 tempconc, tempdose, tempvegdose, tempcowdose
       REAL*8 altdose, altdoseb, altdosec
       REAL*8 FACTOR, FACTORB, FACTORC, FACTORD, FACTORE
       REAL*8 altdose_total, totvegdose, itotvegdose, totcowdose, itotcowdose
c
c Stack plume rise and atmospheric parameters:
c   RHOAIR = air density (kg/m3)
c   RHOSTK = stack gas density (kg/m3)
c   STKVEL = stack gas velocity (m/s)
c   R = velocity ratio parameter
c   FM = momentum flux parameter
c   DELHMO = plume rise due to momentum (m)
c   DELHBU = plume rise due to buoyancy (m)
c   HEFF = effective release height (m)
c   FB = buoyancy flux parameter

       REAL*8 RHOAIR, RHOSTK, STKVEL, R, FM, DELHMO, DELHBU, HEFF, FB
       REAL*8 TERM1, TERM2
       REAL*8 windspeed_cubed, bldheight_threshold, fac_scaling

c Primary loop indices and control variables

       INTEGER nstacks
       INTEGER N, M, I, j, jj, receptor, ios
       INTEGER date_time(8)

c Facility metadata and temporary input variables

       CHARACTER*200 title, company, facility, address, poc, phone, inline
       CHARACTER*8 NUC, radname_trimmed
       CHARACTER*2 NN
       CHARACTER*10 b(3)

c
c ============================================================================
c SECTION 3: DEFINE MATHEMATICAL CONSTANTS AND TUNING PARAMETERS
c ============================================================================
c

       INTEGER, PARAMETER :: nrads_max = 3000

c Physical constants

       REAL*8, PARAMETER :: PI = 3.14159265358979d0
c
c Wind and stack parameter limits

       REAL*8, PARAMETER :: WIND_MIN = 0.01d0
       REAL*8, PARAMETER :: FLOWRATE_DEFAULT = 0.3d0
c
c Temperature constants for density calculations

       REAL*8, PARAMETER :: CONST_80 = 80.0d0
       REAL*8, PARAMETER :: CONST_460 = 460.0d0
c
c Air density scaling factors (metric to English units conversion)

       REAL*8, PARAMETER :: RHOAIR_SCALE1 = 0.0735d0
       REAL*8, PARAMETER :: RHOAIR_SCALE2 = 454.0d0
       REAL*8, PARAMETER :: UNIT_CONV = 0.3048d0
c
c Plume rise coefficients per Briggs equations

       REAL*8, PARAMETER :: CRITICAL_HEIGHT_FACTOR = 2.5d0
       REAL*8, PARAMETER :: MOMENTUM_CONST1 = 0.93d0
       REAL*8, PARAMETER :: MOMENTUM_CONST2 = 7.79d0
       REAL*8, PARAMETER :: BUOYANCY_CONST1 = 23.67d0
       REAL*8, PARAMETER :: BUOYANCY_CONST2 = 0.6d0
       REAL*8, PARAMETER :: BUOYANCY_CONST3 = 195.1d0
       REAL*8, PARAMETER :: BUOYANCY_CONST4 = 0.4d0
c
c Wind frequency and dispersion model parameters

       REAL*8, PARAMETER :: WIND_FREQ_FACTOR = 0.25d0
       REAL*8, PARAMETER :: ZZ_THRESHOLD = 0.001d0
c
c Distance sampling grid for dose calculations

       REAL*8, PARAMETER :: DIST_START = 10.0d0
       REAL*8, PARAMETER :: DIST_END = 9000.0d0
       REAL*8, PARAMETER :: DIST_STEP = 10.0d0
c
c Gaussian plume dispersion model coefficients

       REAL*8, PARAMETER :: SIGMA_COEFF = 0.06d0
       REAL*8, PARAMETER :: SIGMA_DENOM_COEFF = 0.0015d0
       REAL*8, PARAMETER :: PP_CONST = 2.032d0

c ============================================================================
c SECTION 4: INITIALIZE PROGRAM AND OPEN I/O FILES
c ============================================================================

       WRITE(*,*)
       WRITE(*,*)'COMPLY v2.0'
       WRITE(*,*)

c Open input file with facility and stack parameters

       OPEN(11, FILE='input.level3.txt', STATUS='old', 
     &      FORM='formatted', IOSTAT=ios)
       IF (ios /= 0) THEN
           WRITE(*,*) 'ERROR: Cannot open input.level3.txt'
           STOP
       ENDIF

c Open output files for results and diagnostics

       OPEN(12, FILE='output.level3.txt', STATUS='unknown', 
     &      FORM='formatted')
       OPEN(13, FILE='diag.level3.txt', STATUS='unknown', 
     &      FORM='formatted')
c
c Open radionuclide dose factor lookup table from Oak Ridge database

       OPEN(9, FILE='comply171_original_input_data.csv', 
     &      STATUS='old', FORM='formatted', IOSTAT=ios)
       IF (ios /= 0) THEN
           WRITE(*,*) 'ERROR: Cannot open comply171_original_input_data.csv'
           STOP
       ENDIF

c ============================================================================
c SECTION 5: READ FACILITY METADATA AND INITIAL STACK CONFIGURATION
c ============================================================================

        READ(11, '(38x,a)', IOSTAT=ios) title
        IF (ios /= 0) STOP 'Error reading title'
        READ(11, '(38x,a)') company
        READ(11, '(38x,a)') facility
        READ(11, '(38x,a)') address
        READ(11, '(38x,a)') poc
        READ(11, '(38x,a)') phone
        READ(11, '(38x,i3)') nstacks

        IF (nstacks <= 0) THEN
            WRITE(*,*) 'ERROR: nstacks must be > 0'
            STOP
        ENDIF

c ============================================================================
c SECTION 6: ALLOCATE ARRAYS AND READ INITIAL STACK DATA
c ============================================================================
c
c Pre-allocate arrays based on number of stacks; sized to determine maximum
c number of radionuclides before allocating 2D arrays

        ALLOCATE(latitude(nstacks))
        ALLOCATE(longitude(nstacks))
        ALLOCATE(releaseht(nstacks))
        ALLOCATE(flowrate(nstacks))
        ALLOCATE(diameter(nstacks))
        ALLOCATE(stacktemp(nstacks))
        ALLOCATE(airtemp(nstacks))
        ALLOCATE(verticalmomentum(nstacks))
        ALLOCATE(bldheight(nstacks))
        ALLOCATE(bldwidth(nstacks))
        ALLOCATE(distance(nstacks))
        ALLOCATE(distanceveg(nstacks))
        ALLOCATE(distancecow(nstacks))
        ALLOCATE(windspeed(nstacks))
        ALLOCATE(nrads(nstacks))
        ALLOCATE(radname(nstacks, nrads_max))
        ALLOCATE(rademis(nstacks, nrads_max))
        ALLOCATE(dose(nstacks, nrads_max))
        ALLOCATE(wbfac(nstacks, nrads_max))
        ALLOCATE(vegfac(nstacks, nrads_max))
        ALLOCATE(cowfac(nstacks, nrads_max))
        ALLOCATE(vegdose(nstacks, nrads_max))
        ALLOCATE(cowdose(nstacks, nrads_max))

        wbfac = 0.0d0
        vegfac = 0.0d0
        cowfac = 0.0d0
        dose = 0.0d0
        vegdose = 0.0d0
        cowdose = 0.0d0
        rademis = 0.0d0

c Read stack data from input file and determine maximum number of radionuclides
c This allows proper dimensioning of 2D radionuclide arrays

        DO n = 1, nstacks
            READ(11, '(a)', IOSTAT=ios)  !blank spacer
            IF (ios /= 0) STOP 'Error reading blank line'
            READ(11, '(38x,f10.5)') latitude(n)
            READ(11, '(38x,f10.5)') longitude(n)
            READ(11, '(38x,f10.5)') releaseht(n)
            READ(11, '(38x,f10.5)') bldheight(n)
            READ(11, '(38x,f10.5)') bldwidth(n)
            READ(11, '(38x,f10.5)') flowrate(n)
            READ(11, '(38x,f10.5)') diameter(n)
            READ(11, '(38x,f10.5)') stacktemp(n)
            READ(11, '(38x,f10.5)') airtemp(n)
            READ(11, '(38x,a)') verticalmomentum(n)
            READ(11, '(38x,f10.5)') distance(n)
            READ(11, '(38x,f10.5)') distanceveg(n)
            READ(11, '(38x,f10.5)') distancecow(n)
            READ(11, '(38x,f10.5)') windspeed(n)
            READ(11, '(38x,i3)') nrads(n)

c ============================================================================
c SECTION 8: REWIND AND READ RADIONUCLIDE EMISSION DATA
c ============================================================================
c
            airtemp(n) = airtemp(n) * (9./5.) + 32. !convert C to F
            stacktemp(n) = stacktemp(n) * (9./5.) + 32. !convert C to F

            DO M = 1, nrads(n)
                READ(11, '(38x,a)') inline
                READ(inline, *) NUC, rademis(n,m)
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


                ENDDO
         
        ENDDO !end loop over stacks

c ============================================================================
c SECTION 9: READ DOSE CONVERSION FACTORS FROM LOOKUP TABLE
c ============================================================================
c
c Load radionuclide-specific dose conversion factors from Oak Ridge database
c These factors convert air concentration (Ci/m3) to dose equivalent (mrem/yr)
c for three major pathways: inhalation, vegetable ingestion, meat/milk ingestion
c
c File format: NUC | FACTOR | FACTORC | FACTORD | FACTORE
c   NUC = radionuclide identifier (e.g., "C-14", "I-131")
c   FACTORC = whole body dose factor (for inhalation pathway)
c   FACTORD = vegetable dose factor (for ingestion pathway)
c   FACTORE = meat/milk dose factor (for ingestion pathway)

        WRITE(*,*)
        READ(9, '(a)', IOSTAT=ios) inline  !skip header

        DO jj = 1, 826
            READ(9, '(a)', IOSTAT=ios) inline
            IF (ios /= 0) EXIT  !stop if end of file reached
            READ(inline, *, IOSTAT=ios) NUC, FACTOR, FACTORB, 
     &               FACTORC, FACTORD, FACTORE 
            IF (ios /= 0) CYCLE  !skip malformed lines

            radname_trimmed = ADJUSTL(NUC)
            
c Match dose factors to radionuclides - exit early when found

            DO N = 1, nstacks
                DO M = 1, nrads(n)
                    IF (TRIM(radname_trimmed) .EQ. TRIM(radname(n,m))) 
     &              THEN
                        wbfac(n,m) = FACTORC
                        vegfac(n,m) = FACTORD
                        cowfac(n,m) = FACTORE
                        WRITE(*,'(A,A,E15.6)') 
     &                  'WB dose factor: ', NUC, FACTORC
                        WRITE(*,'(A,A,E15.6)') 
     &                  'Veg dose factor: ', NUC, FACTORD
                        WRITE(*,'(A,A,E15.6)') 
     &                  'Cow dose factor: ', NUC, FACTORE
                        WRITE(*,*)
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

        WRITE(*,*)

c ============================================================================
c SECTION 10: INITIALIZE DOSE ACCUMULATORS
c ============================================================================

        alldose = 0.0d0
        idose = 0.0d0
        totvegdose = 0.0d0
        itotvegdose = 0.0d0
        totcowdose = 0.0d0
        itotcowdose = 0.0d0

        WRITE(13, 3001)  !write header to diagnostic file

c ============================================================================
c SECTION 11: MAIN CALCULATION LOOP - DISPERSION AND DOSE ESTIMATION
c ============================================================================
c
c For each stack, calculate plume rise due to momentum and buoyancy.
c Then iterate through distances from the stack for each radionuclide,
c computing Gaussian plume concentrations and corresponding dose estimates.
c
c Key assumptions:
c   - Emission rate input units: Ci/yr 
c   - WIND_FREQ_FACTOR represents default wind direction frequency (0.25 = 25%)
c

        DO N = 1, nstacks

c Validate wind speed; skip stack if invalid

            IF (windspeed(n) < WIND_MIN) THEN
                WRITE(*,'(A,I3,A)') 
     &          'WARNING: Stack ', N, ' has invalid windspeed'
                CYCLE  !skip this stack
            ENDIF

c Calculate building-related parameters for near-field effects

            xcrit = CRITICAL_HEIGHT_FACTOR * SQRT(bldwidth(n) * bldheight(n))
            bldheight_threshold = CRITICAL_HEIGHT_FACTOR * bldheight(n)

c Set default exit velocity if not specified

            IF (flowrate(n) <= 0.0d0) THEN
                flowrate(n) = FLOWRATE_DEFAULT
            ENDIF

c ===================================================================
c SUBSECTION 11a: PLUME RISE CALCULATIONS (BRIGGS MODEL)
c ===================================================================
c
c Calculate effective release height accounting for plume rise due to:
c   1. Momentum flux (forced plume rise)
c   2. Buoyancy flux (thermal plume rise)
c The larger of the two effects is used.
c

            IF (releaseht(n) > bldheight_threshold) THEN
                
c Calculate air density as a function of ambient temperature

                RHOAIR = RHOAIR_SCALE1 * (CONST_80 + CONST_460) / 
     &                  (airtemp(n) + CONST_460)
                RHOAIR = RHOAIR * RHOAIR_SCALE2 / 
     &                  (UNIT_CONV ** 3)
c
c Calculate stack gas density (accounting for temperature difference)

                RHOSTK = RHOAIR * (airtemp(n) + CONST_460) / 
     &                  (stacktemp(n) + CONST_460)

c Calculate stack velocity parameters

                area = PI * (diameter(n) ** 2) / 4.0d0
                STKVEL = flowrate(n) / area
                R = STKVEL / windspeed(n)
                FM = (RHOSTK / RHOAIR) * STKVEL * flowrate(n) / PI
                
c Scaling factor for wind speed effects on plume rise

                fac_scaling = 1.0d0 + 3.0d0 / R

c Calculate intermediate terms for momentum rise equation

                TERM1 = 75.0d0 * FM / ((fac_scaling * windspeed(n)) 
     &                  ** 2)
                TERM2 = releaseht(n) + MOMENTUM_CONST2 * 
     &                 SQRT(FM) / (fac_scaling * windspeed(n))

c Compute plume rise due to momentum (if flag set to YES)

                IF (TRIM(verticalmomentum(n)) .EQ. 'YES') THEN
                    DELHMO = MOMENTUM_CONST1 * (TERM1 ** (3.0d0/7.0d0)) 
     &                      * (TERM2 ** (1.0d0/7.0d0))
                ELSE
                    DELHMO = 0.0d0
                ENDIF

c Compute plume rise due to buoyancy

                windspeed_cubed = windspeed(n) ** 3
                FB = 9.81d0 * (RHOAIR - RHOSTK) * flowrate(n) / 
     &               (PI * RHOAIR)
                DELHBU = BUOYANCY_CONST1 * 
     &                  ((FB / windspeed_cubed) ** BUOYANCY_CONST2) *
     &                  (releaseht(n) + BUOYANCY_CONST3 * FB / 
     &                   windspeed_cubed) ** BUOYANCY_CONST4

c Determine effective release height (stack height + larger plume rise)

                HEFF = releaseht(n) + MAX(DELHMO, DELHBU)
            ELSE
                HEFF = releaseht(n)
            ENDIF

c ===================================================================
c SUBSECTION 11b: GAUSSIAN PLUME CONCENTRATION AND DOSE CALCULATIONS
c ===================================================================
c
c For each radionuclide at this stack, compute concentration and dose
c at grid of distances from the stack

            DO M = 1, nrads(n)


            do receptor = 1, 3, 1

             if(receptor.eq.1) then
             dist = distance(n)
             elseif(receptor.eq.2) then
             dist = distanceveg(n)
             elseif(receptor.eq.3) then
             dist = distancecow(n)
             else
             endif


c === Compute Gaussian plume dispersion parameters ===

                    SIGMAZ = (SIGMA_COEFF * dist) / 
     &                       SQRT(1.0d0 + SIGMA_DENOM_COEFF * dist)
                    CAPSIG = SQRT(SIGMAZ ** 2 + (bldwidth(n) * bldheight(n)) / PI)
                    zz = 0.5d0 * ((HEFF / CAPSIG) ** 2)

c Compute ground-level concentration factor accounting for plume 
c geometry and whether plume height is above building

                    IF (zz > ZZ_THRESHOLD .AND. 
     &                  releaseht(n) > bldheight_threshold) THEN
                        pp = PP_CONST * EXP(-zz) / (dist * CAPSIG)
                    ELSE
                        pp = PP_CONST / (dist * CAPSIG)
                    ENDIF

c === Compute radionuclide air concentrations ===
c
c Two estimates: tempconc (near-field aware) and altconc (far-field)

                    IF (dist < xcrit) THEN
                        tempconc = (WIND_FREQ_FACTOR * rademis(n,m)) / 
     &                             (PI * windspeed(n) * 
     &                              MIN(bldwidth(n), bldheight(n)))
                    ELSE
                        tempconc = (WIND_FREQ_FACTOR * rademis(n,m) * 
     &                             pp) / windspeed(n)
                    ENDIF

                    altconc = (WIND_FREQ_FACTOR * rademis(n,m) * pp) / 
     &                        windspeed(n)

c === Compute dose estimates for three pathways ===
c
c Apply dose conversion factors (Sv or mrem per unit concentration-time)
c to the air concentration to estimate annual effective dose equivalent

                    tempdose = wbfac(n,m) * tempconc
                    tempvegdose = vegfac(n,m) * tempconc
                    tempcowdose = cowfac(n,m) * tempconc
                    altdose = wbfac(n,m) * altconc
                    altdoseb = vegfac(n,m) * altconc
                    altdosec = cowfac(n,m) * altconc

c === Write diagnostic output for all distances ===
c
c Concentration output as Ci/m3 and dose as mrem/year

                    WRITE(13, 3000) N, TRIM(radname(n,m)), dist,
     &                  tempconc, tempdose, tempvegdose, tempcowdose
   

c ===================================================================
c SUBSECTION 11c: ACCUMULATE DOSES AT SPECIFIED RECEPTOR DISTANCES
c ===================================================================
c
c Check if current distance matches one of the three specified receptor
c distances (meat/milk, vegetable, or general population). When distance
c matches (within 0.1m tolerance), accumulate dose totals and track
c iodine-specific doses separately.

                    IF (receptor.eq.3) THEN
                        totcowdose = totcowdose + tempcowdose
                        cowdose(n,m) = tempcowdose
                        IF (radname(n,m)(1:2) .EQ. 'I-') THEN
                            itotcowdose = itotcowdose + tempcowdose
                        ENDIF
                    ENDIF

                    IF (receptor.eq.2) THEN
                        totvegdose = totvegdose + tempvegdose
                        vegdose(n,m) = tempvegdose
                        IF (radname(n,m)(1:2) .EQ. 'I-') THEN
                            itotvegdose = itotvegdose + tempvegdose
                        ENDIF
                    ENDIF

                    IF (receptor.eq.1) THEN
                        alldose = alldose + tempdose
                        dose(n,m) = tempdose
                        IF (radname(n,m)(1:2) .EQ. 'I-') THEN
                            idose = idose + tempdose
                        ENDIF
                    ENDIF

                ENDDO  !end loop over distances

            ENDDO  !end loop over radionuclides

        ENDDO  !end loop over stacks

c ============================================================================
c SECTION 12: CALCULATE AND REPORT TOTAL DOSE
c ============================================================================
c
c Combine doses from all stacks and all pathways. Compare against EPA
c regulatory limits:
c   - Total effective dose: 10 mrem/yr (all pathways)
c   - Iodine-131 dose: 3 mrem/yr (thyroid critical organ)

        WRITE(*,*)
        altdose_total = alldose + totvegdose + totcowdose
        IF (altdose_total > 0.0d0) THEN
            WRITE(*,'(A,F12.2)') 'Total dose (mrem/yr): ', altdose_total
            WRITE(*,*)'Total dose limit is 10.0 mrem/yr'
            WRITE(*,'(A,F12.2)') 'Total iodine dose (mrem/yr): ',
     &          idose + itotvegdose + itotcowdose
            WRITE(*,*)'Iodine limit is 3.0 mrem/yr'
            WRITE(*,*)
            WRITE(*,'(A,F8.4)') 
     &      'Total dose vegetable fraction: ',
     &      totvegdose / altdose_total
            WRITE(*,'(A,F8.4)') 
     &      'Total dose meat/milk fraction: ',
     &      totcowdose / altdose_total
            WRITE(*,'(A,F8.4)') 
     &      'Total dose other pathways fraction: ',
     &      alldose / altdose_total
        ELSE
            WRITE(*,*) 'WARNING: Total dose is zero'
        ENDIF
        WRITE(*,*)

c ============================================================================
c SECTION 13: WRITE STANDARDIZED OUTPUT REPORT
c ============================================================================

        WRITE(12, 4000)
        WRITE(12, 2001)
        WRITE(12, 1000) TRIM(title)
        WRITE(12, 1001) TRIM(company)
        WRITE(12, 1002) TRIM(facility)
        WRITE(12, 1003) TRIM(address)
        WRITE(12, 1004) TRIM(poc)
        WRITE(12, 1005) TRIM(phone)

        WRITE(12, 2001)
        WRITE(12, 2000) altdose_total
        WRITE(12, 2003) idose + itotvegdose + itotcowdose
        WRITE(12, 2002)


        DO N = 1, nstacks
            WRITE(12, 2001)
            WRITE(12, 1006) N
            WRITE(12, 1007) latitude(n)
            WRITE(12, 1008) longitude(n)
            WRITE(12, 1009) releaseht(n)
            WRITE(12, 1010) diameter(n)
            WRITE(12, 1011) flowrate(n)
            WRITE(12, 1012) stacktemp(n)
            WRITE(12, 1013) airtemp(n)
            WRITE(12, 1014) TRIM(verticalmomentum(n))
            WRITE(12, 1015) distance(n)
            WRITE(12, 1016) distanceveg(n)
            WRITE(12, 1017) distancecow(n)
            WRITE(12, 1018) windspeed(n)
            
            DO M = 1, nrads(n)
                WRITE(12, 2002)
                WRITE(12, 1019) TRIM(radname(n,m))
                WRITE(12, 1020) rademis(n,m)*60.*60.*24.*365.
                WRITE(12, 1021) 
     &              (dose(n,m) + vegdose(n,m) + cowdose(n,m))
                WRITE(12, 1022) dose(n,m)
                WRITE(12, 1023) vegdose(n,m)
                WRITE(12, 1024) cowdose(n,m)
            ENDDO
        ENDDO


!        WRITE(12, 2001)
        WRITE(12, 2001)
        WRITE(12, 4000)

c ============================================================================
c SECTION 14: ADD EXECUTION TIMESTAMP
c ============================================================================

        CALL DATE_AND_TIME(b(1), b(2), b(3), date_time)

        PRINT *, 'Program execution date and time information:'
        PRINT *, 'Year=', date_time(1)
        PRINT *, 'Month=', date_time(2)
        PRINT *, 'Day=', date_time(3)
        PRINT *, 'Hour=', date_time(5)
        PRINT *, 'Minute=', date_time(6)
        PRINT *, 'Second=', date_time(7)

        WRITE(12, *)
        WRITE(12, 2100) date_time(2), date_time(3), date_time(1),
     &                   date_time(5), date_time(6)

c ============================================================================
c SECTION 15: FORMAT STATEMENTS FOR OUTPUT FILE
c ============================================================================

 1000 FORMAT('Project: ', A)
 1001 FORMAT('Company: ', A)
 1002 FORMAT('Facility: ', A)
 1003 FORMAT('Address: ', A)
 1004 FORMAT('Point of Contact: ', A)
 1005 FORMAT('Phone Number: ', A)
 1006 FORMAT('Stack number: ', I3)
 1007 FORMAT('Latitude (degrees): ', F10.5)
 1008 FORMAT('Longitude (degrees): ', F10.5)
 1009 FORMAT('Release height (m): ', F10.1)
 1010 FORMAT('Stack diameter (m): ', F10.1)
 1011 FORMAT('Flow rate (m3/s): ', F10.1)
 1012 FORMAT('Stack temperature (F): ', F10.1)
 1013 FORMAT('Air temperature (F): ', F10.1)
 1014 FORMAT('Vertical momentum (YES/NO): ', A)
 1015 FORMAT('Distance from source to receptor (m): ', F10.1)
 1016 FORMAT('Distance from source to vegetable receptor (m): ', F10.1)
 1017 FORMAT('Distance from source to meat/milk receptor (m): ', F10.1)
 1018 FORMAT('Assumed wind speed (m/s) when no wind rose: ', F10.1)
 1019 FORMAT('Radionuclide: ', A)
 1020 FORMAT('Emission rate (Ci/y): ', ES9.2)
 1021 FORMAT('Effective dose equivalent all pathways (mrem/yr): ',
     &       ES9.2)
 1022 FORMAT('Effective dose equivalent non-ingestion (mrem/yr): ',
     &       ES9.2)
 1023 FORMAT('Effective dose equivalent vegetables (mrem/yr): ',
     &       ES9.2)
 1024 FORMAT('Effective dose equivalent meat/milk (mrem/yr): ',
     &       ES9.2)

 2000 FORMAT('Effective dose equivalent all stacks (mrem/yr): ',
     &       ES9.2)
 2003 FORMAT('Effective dose equivalent all stacks from iodine (mrem/yr): ',ES9.2)
 2001 FORMAT('-------------------------------------------'//
     &       '-------------------------------------------')
 2002 FORMAT(' ')
 2005 FORMAT('Vegetable fraction of total: ', F7.4)
 2006 FORMAT('Meat/milk fraction of total: ', F7.4)
 2007 FORMAT('Other pathway fraction of total: ', F7.4)

 2100 FORMAT('Execution date/time (MM/DD/YYYY HH:MM): ',
     &       I2.2, '/', I2.2, '/', I4, ' ', I2.2, ':', I2.2)

 3000 FORMAT(I3, ',', A, ',', F8.1, 4(',', ES16.8))
 3001 FORMAT('STACK,NUCLIDE,DISTANCE(m),CONC(Ci/m3),TOTALDOSE(mrem/yr),'
     &    'VEGDOSE(mrem/yr),COWDOSE(mrem/yr)')

 4000 FORMAT('COMPLY version 2.0 level 3')

c ============================================================================
c SECTION 16: PROGRAM TERMINATION
c ============================================================================

        CLOSE(9)
        CLOSE(11)
        CLOSE(12)
        CLOSE(13)
        
        STOP

       END
