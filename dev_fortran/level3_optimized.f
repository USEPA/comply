       PROGRAM COMPLY

c Declare dynamic arrays

       REAL*4, ALLOCATABLE :: releaseht(:)           !stack height (m)      
       REAL*4, ALLOCATABLE :: diameter(:)            !stack diameter (m)
       REAL*4, ALLOCATABLE :: exitvelocity(:)        !exit velocity (m/s)
       REAL*4, ALLOCATABLE :: airtemp(:)             !air temperature (F)
       REAL*4, ALLOCATABLE :: stacktemp(:)           !stack temperature (F)
       REAL*4, ALLOCATABLE :: bldheight(:)           !building height (m)
       REAL*4, ALLOCATABLE :: bldwidth(:)            !building width (m)

       REAL*4, ALLOCATABLE :: distance(:)        !distance between stack and receptor (m)
       REAL*4, ALLOCATABLE :: distanceveg(:)     !distance between garden and receptor (m)
       REAL*4, ALLOCATABLE :: distancecow(:)     !distance between cows and receptor (m)
       REAL*4, ALLOCATABLE :: windspeed(:)       !wind speed (m/s)
       REAL*4, ALLOCATABLE :: latitude(:)        !latitude (decimal degrees)
       REAL*4, ALLOCATABLE :: longitude(:)       !longitude (decimal degrees)

       REAL*8, ALLOCATABLE :: rademis(:,:)       !radionuclide emissions (Ci/s) for each stack
       REAL*8, ALLOCATABLE :: dose(:,:)          !radionuclide dose (mrem/yr)
       REAL*8, ALLOCATABLE :: wbfac(:,:)         !whole body dose factor from Oak Ridge table
       REAL*8, ALLOCATABLE :: vegfac(:,:)        !vegetable dose factor from Oak Ridge table
       REAL*8, ALLOCATABLE :: cowfac(:,:)        !cow dose factor from Oak Ridge table
       REAL*8, ALLOCATABLE :: vegdose(:,:)       !radionuclide dose by vegetables (mrem/yr)
       REAL*8, ALLOCATABLE :: cowdose(:,:)       !radionuclide dose by meat and milk (mrem/yr)

       INTEGER, ALLOCATABLE :: nrads(:)          !number of emitted radionuclides at each stack

       CHARACTER*8, ALLOCATABLE :: radname(:,:)  !radionuclide names at each stack
       CHARACTER*3, ALLOCATABLE :: verticalmomentum(:)    !vertical momentum (YES or NO)

c Declare other variables and calculated terms

       REAL*8 tempconc, tempdose, tempvegdose, tempcowdose
       REAL*8 altdose, altdoseb, altdosec
       REAL*8 FACTOR, FACTORC, FACTORD, FACTORE
       REAL*8 area, xcrit, SIGMAZ, CAPSIG, pp, alldose, idose, dist, altconc
       REAL*8 altdose_total, totvegdose, itotvegdose, totcowdose, itotcowdose
       REAL*8 RHOAIR, RHOSTK, STKVEL, R, FM, DELHMO, DELHBU, HEFF, FB, zz
       REAL*8 TERM1, TERM2, THRESHOLD_HEIGHT
       REAL*8 windspeed_cubed, bldheight_threshold, fac_scaling

       INTEGER nstacks, nrads_max
       INTEGER N, M, J, receptor, jj, ios
       INTEGER date_time(8)

       CHARACTER*200 title, company, facility, address, poc, phone, inline
       CHARACTER*8 NUC, radname_trimmed
       CHARACTER*2 NN
       CHARACTER*10 b(3)

       REAL*8, PARAMETER :: PI = 3.14159265358979d0
       REAL*8, PARAMETER :: WIND_MIN = 0.01d0
       REAL*8, PARAMETER :: EXIT_VEL_DEFAULT = 0.3d0
       REAL*8, PARAMETER :: CONST_80 = 80.0d0
       REAL*8, PARAMETER :: CONST_460 = 460.0d0
       REAL*8, PARAMETER :: RHOAIR_SCALE1 = 0.0735d0
       REAL*8, PARAMETER :: RHOAIR_SCALE2 = 454.0d0
       REAL*8, PARAMETER :: UNIT_CONV = 0.3048d0
       REAL*8, PARAMETER :: CRITICAL_HEIGHT_FACTOR = 2.5d0
       REAL*8, PARAMETER :: WIND_FREQ_FACTOR = 0.25d0
       REAL*8, PARAMETER :: ZZ_THRESHOLD = 0.001d0
       REAL*8, PARAMETER :: MOMENTUM_CONST1 = 0.93d0
       REAL*8, PARAMETER :: MOMENTUM_CONST2 = 7.79d0
       REAL*8, PARAMETER :: BUOYANCY_CONST1 = 23.67d0
       REAL*8, PARAMETER :: BUOYANCY_CONST2 = 0.6d0
       REAL*8, PARAMETER :: BUOYANCY_CONST3 = 195.1d0
       REAL*8, PARAMETER :: BUOYANCY_CONST4 = 0.4d0
       REAL*8, PARAMETER :: DIST_START = 10.0d0
       REAL*8, PARAMETER :: DIST_END = 2000.0d0
       REAL*8, PARAMETER :: DIST_STEP = 10.0d0
       REAL*8, PARAMETER :: SIGMA_COEFF = 0.06d0
       REAL*8, PARAMETER :: SIGMA_DENOM_COEFF = 0.0015d0
       REAL*8, PARAMETER :: PP_CONST = 2.032d0

c Open input and output files 

       WRITE(*,*)
       WRITE(*,*)'COMPLY v2.0 for internal testing only'
       WRITE(*,*)

       OPEN(11, FILE='input.level3.txt', STATUS='old', 
     &      FORM='formatted', IOSTAT=ios)
       IF (ios /= 0) THEN
           WRITE(*,*) 'ERROR: Cannot open input.level3.txt'
           STOP
       ENDIF

       OPEN(12, FILE='output.level3.txt', STATUS='unknown', 
     &      FORM='formatted')
       OPEN(13, FILE='diag.level3.txt', STATUS='unknown', 
     &      FORM='formatted')
       OPEN(9, FILE='comply171_original_input_data.csv', 
     &      STATUS='old', FORM='formatted', IOSTAT=ios)
       IF (ios /= 0) THEN
           WRITE(*,*) 'ERROR: Cannot open comply171_original_input_data.csv'
           STOP
       ENDIF

c Read input namelist file

        READ(11, '(20x,a)', IOSTAT=ios) title
        IF (ios /= 0) STOP 'Error reading title'
        READ(11, '(20x,a)') company
        READ(11, '(20x,a)') facility
        READ(11, '(20x,a)') address
        READ(11, '(20x,a)') poc
        READ(11, '(20x,a)') phone
        READ(11, '(20x,i3)') nstacks

        IF (nstacks <= 0) THEN
            WRITE(*,*) 'ERROR: nstacks must be > 0'
            STOP
        ENDIF

c Allocate arrays for stack properties

        ALLOCATE(latitude(nstacks))
        ALLOCATE(longitude(nstacks))
        ALLOCATE(releaseht(nstacks))
        ALLOCATE(exitvelocity(nstacks))
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

c Read stack data and determine max radionuclides

        nrads_max = 0
        DO N = 1, nstacks
            READ(11, '(a)', IOSTAT=ios)  !blank spacer
            IF (ios /= 0) STOP 'Error reading blank line'
            READ(11, '(20x,f10.5)') latitude(n)
            READ(11, '(20x,f10.5)') longitude(n)
            READ(11, '(20x,f10.5)') releaseht(n)
            READ(11, '(20x,f10.5)') bldheight(n)
            READ(11, '(20x,f10.5)') bldwidth(n)
            READ(11, '(20x,f10.5)') exitvelocity(n)
            READ(11, '(20x,f10.5)') diameter(n)
            READ(11, '(20x,f10.5)') stacktemp(n)
            READ(11, '(20x,f10.5)') airtemp(n)
            READ(11, '(20x,a)') verticalmomentum(n)
            READ(11, '(20x,f10.5)') distance(n)
            READ(11, '(20x,f10.5)') distanceveg(n)
            READ(11, '(20x,f10.5)') distancecow(n)
            READ(11, '(20x,f10.5)') windspeed(n)
            READ(11, '(20x,i3)') nrads(n)

            IF (nrads(n) > nrads_max) nrads_max = nrads(n)
        ENDDO

c Allocate 2D arrays with correct dimensions

        ALLOCATE(radname(nstacks, nrads_max))
        ALLOCATE(rademis(nstacks, nrads_max))
        ALLOCATE(dose(nstacks, nrads_max))
        ALLOCATE(wbfac(nstacks, nrads_max))
        ALLOCATE(vegfac(nstacks, nrads_max))
        ALLOCATE(cowfac(nstacks, nrads_max))
        ALLOCATE(vegdose(nstacks, nrads_max))
        ALLOCATE(cowdose(nstacks, nrads_max))

c Initialize dose factor arrays
        wbfac = 0.0d0
        vegfac = 0.0d0
        cowfac = 0.0d0

c Rewind to read radionuclide data
        REWIND(11)
        READ(11, '(20x,a)') title
        READ(11, '(20x,a)') company
        READ(11, '(20x,a)') facility
        READ(11, '(20x,a)') address
        READ(11, '(20x,a)') poc
        READ(11, '(20x,a)') phone
        READ(11, '(20x,i3)') nstacks

        DO N = 1, nstacks
            READ(11, '(a)')
            READ(11, '(20x,f10.5)') latitude(n)
            READ(11, '(20x,f10.5)') longitude(n)
            READ(11, '(20x,f10.5)') releaseht(n)
            READ(11, '(20x,f10.5)') bldheight(n)
            READ(11, '(20x,f10.5)') bldwidth(n)
            READ(11, '(20x,f10.5)') exitvelocity(n)
            READ(11, '(20x,f10.5)') diameter(n)
            READ(11, '(20x,f10.5)') stacktemp(n)
            READ(11, '(20x,f10.5)') airtemp(n)
            READ(11, '(20x,a)') verticalmomentum(n)
            READ(11, '(20x,f10.5)') distance(n)
            READ(11, '(20x,f10.5)') distanceveg(n)
            READ(11, '(20x,f10.5)') distancecow(n)
            READ(11, '(20x,f10.5)') windspeed(n)
            READ(11, '(20x,i3)') nrads(n)

            DO M = 1, nrads(n)
                READ(11, '(20x,a)') inline
                READ(inline, *) radname(n,m), rademis(n,m)
            ENDDO
        ENDDO

c Read in dose factor data

        WRITE(*,*)
        REWIND(9)
        READ(9, '(a)', IOSTAT=ios) inline  !skip header

        DO jj = 1, 826
            READ(9, '(a)', IOSTAT=ios) inline
            IF (ios /= 0) EXIT  !stop if end of file reached
            READ(inline, *, IOSTAT=ios) NUC, FACTOR, FACTORC, FACTORD, 
     &                                   FACTORE
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

c Initialize dose accumulators

        alldose = 0.0d0
        idose = 0.0d0
        totvegdose = 0.0d0
        itotvegdose = 0.0d0
        totcowdose = 0.0d0
        itotcowdose = 0.0d0

        WRITE(13, 3001)  !write header to diagnostic file

c Main calculation loop - Calculate air quality impacts and dose
c Emission rate assumption is Ci/s
c WIND_FREQ_FACTOR is default wind direction frequency fraction

        DO N = 1, nstacks

            IF (windspeed(n) < WIND_MIN) THEN
                WRITE(*,'(A,I3,A)') 
     &          'WARNING: Stack ', N, ' has invalid windspeed'
                CYCLE  !skip this stack
            ENDIF

            area = bldwidth(n) * bldheight(n)
            xcrit = CRITICAL_HEIGHT_FACTOR * SQRT(area)
            bldheight_threshold = CRITICAL_HEIGHT_FACTOR * bldheight(n)

            IF (exitvelocity(n) <= 0.0d0) THEN
                exitvelocity(n) = EXIT_VEL_DEFAULT
            ENDIF

c Precalculate temperature and momentum terms outside inner loop

            IF (releaseht(n) > bldheight_threshold) THEN
                
                RHOAIR = RHOAIR_SCALE1 * (CONST_80 + CONST_460) / 
     &                  (airtemp(n) + CONST_460)
                RHOAIR = RHOAIR * RHOAIR_SCALE2 / 
     &                  (UNIT_CONV ** 3)
                RHOSTK = RHOAIR * (airtemp(n) + CONST_460) / 
     &                  (stacktemp(n) + CONST_460)

                area = PI * (diameter(n) ** 2) / 4.0d0
                STKVEL = exitvelocity(n) / area
                R = STKVEL / windspeed(n)
                FM = (RHOSTK / RHOAIR) * STKVEL * exitvelocity(n) / PI
                
                fac_scaling = 1.0d0 + 3.0d0 / R

                TERM1 = 75.0d0 * FM / ((fac_scaling * windspeed(n)) 
     &                  ** 2)
                TERM2 = releaseht(n) + MOMENTUM_CONST2 * 
     &                 SQRT(FM) / (fac_scaling * windspeed(n))

                IF (TRIM(verticalmomentum(n)) .EQ. 'YES') THEN
                    DELHMO = MOMENTUM_CONST1 * (TERM1 ** (3.0d0/7.0d0)) 
     &                      * (TERM2 ** (1.0d0/7.0d0))
                ELSE
                    DELHMO = 0.0d0
                ENDIF

                windspeed_cubed = windspeed(n) ** 3
                FB = 9.81d0 * (RHOAIR - RHOSTK) * exitvelocity(n) / 
     &               (PI * RHOAIR)
                DELHBU = BUOYANCY_CONST1 * 
     &                  ((FB / windspeed_cubed) ** BUOYANCY_CONST2) *
     &                  (releaseht(n) + BUOYANCY_CONST3 * FB / 
     &                   windspeed_cubed) ** BUOYANCY_CONST4

                HEFF = releaseht(n) + MAX(DELHMO, DELHBU)
            ELSE
                HEFF = releaseht(n)
            ENDIF

            DO M = 1, nrads(n)

                DO receptor = INT(DIST_START), INT(DIST_END), 
     &                        INT(DIST_STEP)

                    dist = DBLE(receptor)

c Calculate dispersion parameters

                    SIGMAZ = (SIGMA_COEFF * dist) / 
     &                       SQRT(1.0d0 + SIGMA_DENOM_COEFF * dist)
                    CAPSIG = SQRT(SIGMAZ ** 2 + area / PI)
                    zz = 0.5d0 * ((HEFF / CAPSIG) ** 2)

                    IF (zz > ZZ_THRESHOLD .AND. 
     &                  releaseht(n) > bldheight_threshold) THEN
                        pp = PP_CONST * EXP(-zz) / (dist * CAPSIG)
                    ELSE
                        pp = PP_CONST / (dist * CAPSIG)
                    ENDIF

c Calculate concentration

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

c Calculate doses

                    tempdose = wbfac(n,m) * tempconc
                    tempvegdose = vegfac(n,m) * tempconc
                    tempcowdose = cowfac(n,m) * tempconc
                    altdose = wbfac(n,m) * altconc
                    altdoseb = vegfac(n,m) * altconc
                    altdosec = cowfac(n,m) * altconc

c Write diagnostic output

                    WRITE(13, 3000) N, TRIM(radname(n,m)), dist,
     &                  tempconc, tempdose, tempvegdose, tempcowdose,
     &                  altconc, altdose, altdoseb, altdosec

c Check distance matches and accumulate totals

                    IF (ABS(dist - distancecow(n)) < 0.1d0) THEN
                        totcowdose = totcowdose + tempcowdose
                        cowdose(n,m) = tempcowdose
                        IF (radname(n,m)(1:2) .EQ. 'I-') THEN
                            itotcowdose = itotcowdose + tempcowdose
                        ENDIF
                    ENDIF

                    IF (ABS(dist - distanceveg(n)) < 0.1d0) THEN
                        totvegdose = totvegdose + tempvegdose
                        vegdose(n,m) = tempvegdose
                        IF (radname(n,m)(1:2) .EQ. 'I-') THEN
                            itotvegdose = itotvegdose + tempvegdose
                        ENDIF
                    ENDIF

                    IF (ABS(dist - distance(n)) < 0.1d0) THEN
                        alldose = alldose + tempdose
                        dose(n,m) = tempdose
                        IF (radname(n,m)(1:2) .EQ. 'I-') THEN
                            idose = idose + tempdose
                        ENDIF
                    ENDIF

                ENDDO  !end loop over distances

            ENDDO  !end loop over radionuclides

        ENDDO  !end loop over stacks

c Calculate and report total dose

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

c Write standardized output

        WRITE(12, 4000)
        WRITE(12, 2001)
        WRITE(12, 1000) TRIM(title)
        WRITE(12, 1001) TRIM(company)
        WRITE(12, 1002) TRIM(facility)
        WRITE(12, 1003) TRIM(address)
        WRITE(12, 1004) TRIM(poc)
        WRITE(12, 1005) TRIM(phone)

        DO N = 1, nstacks
            WRITE(12, 2001)
            WRITE(12, 1006) N
            WRITE(12, 1007) latitude(n)
            WRITE(12, 1008) longitude(n)
            WRITE(12, 1009) releaseht(n)
            WRITE(12, 1010) diameter(n)
            WRITE(12, 1011) exitvelocity(n)
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
                WRITE(12, 1020) rademis(n,m)
                WRITE(12, 1021) 
     &              (dose(n,m) + vegdose(n,m) + cowdose(n,m))
                WRITE(12, 1022) dose(n,m)
                WRITE(12, 1023) vegdose(n,m)
                WRITE(12, 1024) cowdose(n,m)
            ENDDO
        ENDDO

        WRITE(12, 2002)
        WRITE(12, 2001)
        WRITE(12, 2001)
        WRITE(12, 2000) altdose_total
        WRITE(12, 2002)
        
        IF (altdose_total > 0.0d0) THEN
            WRITE(12, 2005) totvegdose / altdose_total
            WRITE(12, 2006) totcowdose / altdose_total
            WRITE(12, 2007) alldose / altdose_total
        ENDIF

        WRITE(12, 2001)
        WRITE(12, 2001)
        WRITE(12, 4000)

c Add time stamp of model simulation to output file

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

c Formatting for output file

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
 1011 FORMAT('Exit velocity (m/s): ', F10.1)
 1012 FORMAT('Stack temperature (F): ', F10.1)
 1013 FORMAT('Air temperature (F): ', F10.1)
 1014 FORMAT('Vertical momentum (YES/NO): ', A)
 1015 FORMAT('Distance from source to receptor (m): ', F10.1)
 1016 FORMAT('Distance from source to vegetable receptor (m): ', F10.1)
 1017 FORMAT('Distance from source to meat/milk receptor (m): ', F10.1)
 1018 FORMAT('Assumed wind speed (m/s) when no wind rose: ', F10.1)
 1019 FORMAT('Radionuclide: ', A)
 1020 FORMAT('Emission rate (Ci/s): ', E25.16)
 1021 FORMAT('Effective dose equivalent all pathways (mrem/yr): ',
     &       E25.16)
 1022 FORMAT('Effective dose equivalent non-ingestion (mrem/yr): ',
     &       E25.16)
 1023 FORMAT('Effective dose equivalent vegetables (mrem/yr): ',
     &       E25.16)
 1024 FORMAT('Effective dose equivalent meat/milk (mrem/yr): ',
     &       E25.16)

 2000 FORMAT('Total effective dose equivalent all stacks (mrem/yr): ',
     &       F10.2)
 2001 FORMAT('-------------------------------------------'//
     &       '-------------------------------------------')
 2002 FORMAT(' ')
 2005 FORMAT('Vegetable fraction of total: ', F7.4)
 2006 FORMAT('Meat/milk fraction of total: ', F7.4)
 2007 FORMAT('Other pathway fraction of total: ', F7.4)

 2100 FORMAT('Execution date/time (MM/DD/YYYY HH:MM): ',
     &       I2.2, '/', I2.2, '/', I4, ' ', I2.2, ':', I2.2)

 3000 FORMAT(I2, ',', A, ',', F8.1, 8(',', E16.8))
 3001 FORMAT('STACK,RADNAME,DISTANCE,CONC,WBDOSE,VEGDOSE,COWDOSE,',
     &       'ALTCONC,ALTWB,ALTVEG,ALTCOW')

 4000 FORMAT('COMPLY version 2.0 level 3 - EPA Radiation Dose Model')

c End of program

        CLOSE(9)
        CLOSE(11)
        CLOSE(12)
        CLOSE(13)
        
        STOP

       END
