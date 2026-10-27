      PROGRAM COMPLY
c
c ============================================================================
c COMPLY v2.0 - Level 1 Radionuclide Possession Fraction Calculator
c ============================================================================
c
c PURPOSE:
c   This program calculates the fraction of maximum permitted radionuclide
c   possession for EPA compliance analysis. It processes facility emission
c   data for multiple stacks and compares against regulatory possession limits.
c
c INPUT FILES:
c   - input.level1.txt: Facility information and radionuclide emissions
c   - comply171_original_possession_data.csv: Reference possession limits
c
c OUTPUT FILES:
c   - output.level1.txt: Formatted compliance report
c
c APPROACH:
c   1. Read facility metadata (title, company, location, contact)
c   2. Read stack-specific data (latitude, longitude, radionuclide emissions)
c   3. Load regulatory possession limits from reference database
c   4. Calculate possession fractions (emissions / limit) by physical state
c   5. Generate formatted compliance report
c
c KEY VARIABLES:
c   - latitude(nstacks): Geographic latitude of each emission stack
c   - longitude(nstacks): Geographic longitude of each emission stack
c   - rademis(nstacks,100): Radionuclide emission rates (Ci/year)
c   - fract(nstacks,100): Calculated possession fractions (unitless)
c   - posmax(nstacks,4): Regulatory limits by physical state (Ci):
c       Index 1: Gaseous form
c       Index 2: Liquid form
c       Index 3: Solid form
c       Index 4: concentration
c   - radtype(nstacks,100): Physical state of radionuclide (G/L/S/C)
c
c ============================================================================

c Declare dynamic arrays
c Dynamic arrays allow program to handle variable numbers of stacks and
c radionuclides without recompilation. Memory allocated once based on nstacks.

      REAL*4, ALLOCATABLE :: latitude(:)
      ! Array of latitude coordinates for each stack
      
      REAL*4, ALLOCATABLE :: longitude(:)
      ! Array of longitude coordinates for each stack
      
      REAL*4, ALLOCATABLE :: rademis(:,:)
      ! 2D array: (stack, radionuclide) -> emission rate in Ci/year
      ! Dimension 2 limited to 100 radionuclides per stack (practical limit)
      
      REAL*4, ALLOCATABLE :: fract(:,:)
      ! 2D array: (stack, radionuclide) -> possession fraction
      ! Calculated as: fract = emission / regulatory_limit
      ! Values > 1.0 indicate potential compliance issues
      
      REAL*4, ALLOCATABLE :: posmax(:,:)
      ! 2D array: (nstacks, 4) -> regulatory possession limits by state
      ! Column 1: Gaseous maximum possession (Ci)
      ! Column 2: Liquid maximum possession (Ci)
      ! Column 3: Solid maximum possession (Ci)
      ! Column 4: Particulate/other maximum possession (Ci)

      INTEGER, ALLOCATABLE :: nrads(:)
      ! 1D array: number of radionuclides at each stack
      ! nrads(n) = count of distinct radionuclides for stack n

      CHARACTER*8, ALLOCATABLE :: radname(:,:)
      ! 2D array: (stack, radionuclide) -> radionuclide name (e.g., 'I-131')
      ! Limited to 8 characters per name
      
      CHARACTER*1, ALLOCATABLE :: radtype(:,:)
      ! 2D array: (stack, radionuclide) -> physical state code
      ! Valid values: 'G' (gas), 'L' (liquid), 'S' (solid), 'C' (particulate)

c Declare other variables and calculated terms

      REAL alldose, idose
      ! alldose = sum of all possession fractions across all stacks/radnuclides
      ! idose = sum of possession fractions for iodine isotopes only
      ! Note: Despite the name 'dose', these are actually possession fractions
      
      INTEGER nstacks, n, m, i, j
      ! nstacks = total number of emission stacks at facility
      ! n = loop counter for stacks (outer loop)
      ! m = loop counter for radionuclides (inner loop)
      ! i = loop counter for reference database records
      ! j = loop counter for physical state columns
      
      CHARACTER*200 title, company, facility, address, poc, phone, inline
      ! title = project name/title from input file
      ! company = company/organization name
      ! facility = facility name (e.g., 'Reactor A')
      ! address = facility physical address
      ! poc = point of contact person name
      ! phone = contact phone number
      ! inline = temporary buffer for reading CSV lines from files
      
      CHARACTER*8 NUC
      ! NUC = radionuclide name read from reference database
      
      CHARACTER*1 radtype_upper
      ! radtype_upper = normalized physical state code (uppercase)
      ! Used for case-insensitive comparison of G/L/S/C codes

      integer date_time(8)
      ! Array populated by date_and_time() intrinsic:
      ! (1)=year, (2)=month, (3)=day, (4)=timezone, (5)=hour,
      ! (6)=minute, (7)=second, (8)=millisecond
      
      character*10 b(3)
      ! Buffer array for date_and_time() (rarely used in this context)

c ============================================================================
c SECTION 1: FILE INITIALIZATION AND USER OUTPUT
c ============================================================================

c Open input and output files 

      write(*,*)
      write(*,*)'COMPLY v2.0'
      write(*,*)
      ! Display banner to console to confirm program started
 
      open(11,file='input.level1.txt',status='old',form='formatted')
      ! Unit 11: Input file with facility metadata and emission data
      ! status='old' requires file must exist; causes error if not found
      
      open(12,file='output.level1.txt',status='unknown',form='formatted')
      ! Unit 12: Output file for compliance report
      ! status='unknown' creates file if it doesn't exist
      
      open(13,file='comply171_original_possession_data.csv',
     &     status='old',form='formatted')
      ! Unit 13: Reference database with regulatory possession limits
      ! Contains 823 radionuclide records with state-specific limits

c ============================================================================
c SECTION 2: READ INPUT FACILITY METADATA
c ============================================================================

c Read input namelist file
c File format: "Label: Value" with label right-justified in 20 columns

       read(11,'(49x,a)') title
       ! Skip first 20 characters (label), read remaining as project title
       
       read(11,'(49x,a)') company
       ! Read company/organization name
       
       read(11,'(49x,a)') facility
       ! Read facility identifier/name
       
       read(11,'(49x,a)') address
       ! Read facility street address
       
       read(11,'(49x,a)') poc
       ! Read point of contact person name
       
       read(11,'(49x,a)') phone
       ! Read contact phone number (kept as string to preserve formatting)
       
       read(11,'(49x,i3)') nstacks
       ! Read number of emission stacks as 3-digit integer
       ! This value determines size for all stack-related arrays

c ============================================================================
c SECTION 3: DYNAMIC MEMORY ALLOCATION
c ============================================================================

c Allocate all arrays ONCE, before loop
c IMPORTANT: Allocation happens after nstacks is known, allowing
c program to handle different facility sizes without recompilation.
c Each array is sized based on nstacks and max expected radionuclides.

       allocate(latitude(nstacks))
       ! Allocate latitude array for all stacks
       
       allocate(longitude(nstacks))
       ! Allocate longitude array for all stacks
       
       allocate(nrads(nstacks))
       ! Allocate counter array for radionuclides per stack
       
       allocate(radname(nstacks,100))
       ! Allocate 2D array for radionuclide names
       ! Assume maximum 100 radionuclides per stack (conservative estimate)
       ! If facility has > 100 radionuclides, this limit must be increased
       
       allocate(rademis(nstacks,100))
       ! Allocate 2D array for emission rates (Ci/year)
       
       allocate(fract(nstacks,100))
       ! Allocate 2D array for possession fractions (unitless)
       
       allocate(radtype(nstacks,100))
       ! Allocate 2D array for physical state codes (G/L/S/C)
       
       allocate(posmax(nstacks,4))
       ! Allocate 2D array for regulatory limits by physical state
       ! Dimensions: (nstacks, 4 physical states)

c ============================================================================
c SECTION 4: READ STACK AND RADIONUCLIDE DATA
c ============================================================================

       do n = 1, nstacks 
       ! Loop over each stack at the facility
       
       read(11,'(a)') 
       ! Skip separator line (e.g., '===...')
       
       read(11,'(49x,f10.5)') latitude(n)
       ! Read latitude in decimal degrees (format: skip 20 chars, read float)
       
       read(11,'(49x,f10.5)') longitude(n)
       ! Read longitude in decimal degrees
       
       read(11,'(49x,i3)') nrads(n)
       ! Read number of radionuclides for this stack

       do m = 1, nrads(n)
       ! Loop over each radionuclide at this stack
       
       read(11,'(49x,a)') inline
       ! Read line containing radionuclide data
       
       read(inline,*) NUC, radtype(n,m), rademis(n,m)
       ! Parse comma-separated values from line:
       !   radname: isotope identifier (e.g., 'I-131')
       !   radtype: physical state (G=gas, L=liquid, S=solid, C=aerosol)
       !   rademis: emission rate in Ci/year

       !convert nuclide name to all uppercase
       do i = 1, len(NUC)
        j = iachar(NUC(i:i))
        if(j.ge.iachar("a").and.j.le.iachar("z")) then
                radname(n,m)(i:i) = achar(iachar(NUC(i:i))-32)
        else
                radname(n,m)(i:i) = NUC(i:i)
        endif
       enddo


       enddo
       ! End loop over radionuclides
       
       enddo
       ! End loop over stacks

c ============================================================================
c SECTION 5: READ AND PROCESS REFERENCE POSSESSION LIMIT DATABASE
c ============================================================================

c Read in ancillary files - store in lookup table FIRST
c This reference file contains maximum possession quantities
c for each radionuclide, organized by physical state.

       write(*,*)
       ! Add blank line to console output
       
       read(13,'(a)') inline  
       ! Skip header line from CSV file
       
c Create hash map: radionuclide name -> index (or use direct search if only 823 lines)

       do i = 1, 823
       ! Loop through all 823 radionuclide records in reference database
       
       read(13,'(a)') inline
       ! Read next line from reference file
       
       read(inline,*) NUC, (posmax(1,j),j=1,4)
       ! Parse reference record:
       !   NUC = radionuclide name (must match names in input file)
       !   posmax(1,1) = gaseous maximum possession (Ci)
       !   posmax(1,2) = liquid maximum possession (Ci)
       !   posmax(1,3) = solid maximum possession (Ci)
       !   posmax(1,4) = concentration maximum possession (Ci)
       ! NOTE: posmax stored in row 1 for all stacks (limit doesn't vary by stack)
       
c Now search only once per radionuclide, not per stack/radionuclide combo
c This nested loop structure (i inside n/m loops) is inefficient but clear.
c Could be optimized by pre-loading all 823 records into arrays before searching.

       do n = 1, nstacks
       ! For each stack...
       
        do m = 1, nrads(n)
        ! ...and each radionuclide at that stack...
        
        if(radname(n,m) .eq. NUC) then
         ! If radionuclide name matches record from reference file:
        
         write(*,*)  
         write(*,*)'Match:', NUC, radtype(n,m)
         ! Print to console: matched radionuclide name and physical state
         
c Normalize case once
c Input file may contain lowercase letters; standardize to uppercase for
c consistent case-insensitive matching against SELECT CASE statement.

         radtype_upper = radtype(n,m)
         ! Copy physical state code
         
         if(radtype_upper .eq. 'g') radtype_upper = 'G'
         ! Convert lowercase 'g' to uppercase 'G'
         
         if(radtype_upper .eq. 'l') radtype_upper = 'L'
         ! Convert lowercase 'l' to uppercase 'L'
         
         if(radtype_upper .eq. 's') radtype_upper = 'S'
         ! Convert lowercase 's' to uppercase 'S'
         
         if(radtype_upper .eq. 'c') radtype_upper = 'C'
         ! Convert lowercase 'c' to uppercase 'C' (for particulate/aerosol)
         
c Use SELECT CASE for efficiency
c SELECT CASE is more efficient than nested IF statements for this lookup.

         select case(radtype_upper)
         ! Branch based on normalized physical state code
         
         case('G')
          ! Gaseous form: use gaseous possession limit (column 1)
          fract(n,m) = rademis(n,m) / posmax(1,1)
          write(*,*)'WARNING: These results are not valid if:'
          write(*,*)'1) A receptor is within 10 m of any release point; or'
          write(*,*)'2) A meat, milk, or vegetable receptor is within 100 m of any release point.'

         case('L')
          ! Liquid form: use liquid possession limit (column 2)
          fract(n,m) = rademis(n,m) / posmax(1,2)
          write(*,*)'WARNING: These results are not valid if:'
          write(*,*)'1) A receptor is within 10 m of any release point; or'
          write(*,*)'2) A meat, milk, or vegetable receptor is within 100 m of any release point.'

         case('S')
          ! Solid form: use solid possession limit (column 3)
          fract(n,m) = rademis(n,m) / posmax(1,3)
          write(*,*)'WARNING: These results are not valid if:'
          write(*,*)'1) A receptor is within 10 m of any release point; or'
          write(*,*)'2) A meat, milk, or vegetable receptor is within 100 m of any release point.'

         case('C')
          ! Concentration (column 4)
          fract(n,m) = rademis(n,m) / posmax(1,4)
          write(*,*)'WARNING: If anyone lives within 3 stack diameters of any release point these results are not valid.'

         case default
          ! Unrecognized physical state code
          write(*,*)'Unrecognized physical state:', radtype(n,m)
          fract(n,m) = 9999.
          ! Set to sentinel value (9999) to flag data quality issue
         
         end select
         ! End case statement
         
        endif
        ! End if-match block
        
        enddo
       ! End inner loop over radionuclides
       
       enddo
       ! End middle loop over stacks
       
       enddo
       ! End outer loop over reference database records

c ============================================================================
c SECTION 6: CALCULATE TOTAL POSSESSION AND IODINE SUBTOTALS
c ============================================================================

c Calculate possession - avoid string I/O in loop
c Sum total possession fraction and separately sum iodine fractions.
c Note: These sums across all radionuclides and all stacks represent
c cumulative possession relative to limits (unitless metric).

       alldose = 0.0
       ! Initialize total possession fraction accumulator
       
       idose = 0.0
       ! Initialize iodine possession fraction accumulator

       do n = 1, nstacks
       ! Loop over all stacks
       
       do m = 1, nrads(n)
       ! Loop over all radionuclides at this stack
       
       alldose = alldose + fract(n,m)
       ! Add this radionuclide's possession fraction to total
       
c Check for iodine directly with INDEX
c Iodine isotopes (I-131, I-125, I-133, etc.) tracked separately.
c INDEX function returns position of substring in string, or 0 if not found.

       if(index(radname(n,m), 'I-') > 0) then
        ! If radionuclide name contains 'I-' (iodine isotope marker):
        
        idose = idose + fract(n,m)
        ! Add to iodine subtotal
       
       endif
       ! End iodine check
       
       enddo
       ! End inner loop over radionuclides
       
       enddo
       ! End outer loop over stacks

       write(*,*)
       write(*,*)
       write(*,*)'Total possession (fraction):', alldose
       ! Display total to console
       
       write(*,*)'Total iodine possession (fraction):', idose
       ! Display iodine subtotal to console

c ============================================================================
c SECTION 7: GENERATE FORMATTED OUTPUT REPORT
c ============================================================================

c Write standardized output (unchanged)
c Output report structured with header, per-stack details, and summary.

       write(12,2001)
       ! Write separator line

       write(12,4000)
       ! Write program identification banner
       
       write(12,2001)
       ! Write separator line
       
       write(12,1000) trim(title)
       ! Write project title
       
       write(12,1001) trim(company)
       ! Write company name
       
       write(12,1002) trim(facility) 
       ! Write facility identifier
       
       write(12,1003) trim(address)
       ! Write facility address
       
       write(12,1004) trim(poc)
       ! Write point of contact
       
       write(12,1005) trim(phone)
       ! Write contact phone number

       write(12,2001)
       ! Write separator line

       write(12,2001)
       ! Write second separator line

       write(12,2000) alldose
       ! Write total possession fraction

       write(12,1999) idose
       ! Write iodine possession fraction subtotal

       write(12,2001)
       ! Write separator line

       write(12,2001)
       ! Write separator line

       do n = 1, nstacks
       ! Loop over all stacks to write per-stack details
       
       write(12,2002)
       ! Write blank line for spacing

       write(12,1008) n
       ! Write stack number
       
       write(12,1006) latitude(n)
       ! Write stack latitude
       
       write(12,1007) longitude(n)
       ! Write stack longitude
       
       write(12,1013) nrads(n)
       ! Write count of radionuclides at this stack
       
       write(12,2002)
       ! Write blank line for spacing
       
       write(12,2003)
       ! Write column headers for radionuclide table

       do m = 1, nrads(n)
       ! Loop over radionuclides at this stack to write detailed data
       
       write(12,1014) radname(n,m), radtype(n,m), 
     &                rademis(n,m), fract(n,m)
       ! Write radionuclide: name, state, emission rate, possession fraction
       
       enddo
       ! End inner loop over radionuclides
       
       enddo
       ! End outer loop over stacks

       write(12,2002)
       ! Write blank line
       
       write(12,2001)
       ! Write separator line
       
       write(12,2001)
       ! Write separator line
       
c ============================================================================
c SECTION 8: ADD EXECUTION TIMESTAMP
c ============================================================================

c Add time stamp
c Record when this calculation was performed for audit trail.

       write(*,*)
       ! Blank line on console
       
       call date_and_time(b(1), b(2), b(3), date_time)
       ! Call system intrinsic to populate date_time array with current date/time
       
       print *, 'Program execution date and time information:'
       ! Console label
       
       print *, 'year=', date_time(1)
       ! Display year (e.g., 2026)
       
       print *, 'month_of_year=', date_time(2)
       ! Display month (1-12)
       
       print *, 'day_of_month=', date_time(3)
       ! Display day (1-31)
       
       print *, 'hour of day=', date_time(5)
       ! Display hour (0-23)
       
       print *, 'minutes of hour=', date_time(6)
       ! Display minute (0-59)
       
       print *, 'seconds of minute=', date_time(7)
       ! Display second (0-59)

       write(12,*)
       ! Blank line in output file
       
       write(12,2100) date_time(2), date_time(3), date_time(1),
     &           date_time(5), date_time(6)
       ! Write formatted timestamp to output file: MM/DD/YYYY HH:MM

c ============================================================================
c SECTION 9: OUTPUT FORMAT STATEMENTS
c ============================================================================

c Formatting for output file
c These are Fortran format statements used by previous WRITE statements.
c Format label numbers (1000, 1001, etc.) referenced in WRITE statements.
c 'a' = string, 'f10.5' = float with 5 decimal places, 'i3' = 3-digit integer

 1000 FORMAT('Project:',a)
      ! Format for project title output
 
 1001 FORMAT('Company:',a)
      ! Format for company name output
 
 1002 FORMAT('Facility:',a)
      ! Format for facility name output
 
 1003 FORMAT('Address:',a)
      ! Format for facility address output
 
 1004 FORMAT('Point of Contact:',a)
      ! Format for point of contact output
 
 1005 FORMAT('Phone Number:',a)
      ! Format for phone number output
 
 1006 FORMAT('Latitude (degrees):',f10.5)
      ! Format for latitude: fixed width, 5 decimal places
 
 1007 FORMAT('Longitude (degrees):',f10.5)
      ! Format for longitude: fixed width, 5 decimal places
 
 1008 FORMAT('Stack number:',i3)
      ! Format for stack number: 3-digit integer
 
 1013 FORMAT('Number of nuclides:',I2)
      ! Format for radionuclide count: 2-digit integer
 
 1014 FORMAT(a,',',a,',',ES9.2,',',F5.1)
      ! Format for radionuclide data row (CSV-style):
      ! radname (string), radtype (string), rademis (float with 20 decimals),
      ! fract (float with 20 decimals)
      ! High precision (20 decimals) preserves calculation accuracy
 
 1999 FORMAT('Total possession iodine fraction:',f5.1)
      ! Format for iodine subtotal: float with 2 decimal places
 
 2000 FORMAT('Total possession fraction:',f5.1)
      ! Format for total possession: float with 2 decimal places
 
 2001 FORMAT('-----------------------------------------------------------------------')
      ! Format for separator line (fixed 71-character underline)
 
 2002 FORMAT(' ')
      ! Format for blank line
 
 2003 FORMAT('Nuclide, Form, Amount(Ci/yr), Possession Fraction')
      ! Format for column headers in detailed radionuclide table
 
 2100 FORMAT('Date/Time (LST) stamp of code execution (MM/DD/YYYY HH:MM):',
     &       i2.2,'/',i2.2,'/',i4,' ',i2.2,':',i2.2)
      ! Format for timestamp: MM/DD/YYYY HH:MM (Local Solar Time)
      ! Uses 2-digit integers with leading zeros for month/day/hour/minute
      ! Uses 4-digit integer for year
 
 4000 FORMAT('COMPLY version 2.0 level 1')
      ! Format for program identification banner

      END 
c ============================================================================
c END OF PROGRAM COMPLY
c ============================================================================
