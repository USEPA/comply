      PROGRAM COMPLY

c Declare dynamic arrays

      REAL*4, ALLOCATABLE :: latitude(:)
      REAL*4, ALLOCATABLE :: longitude(:)
      REAL*4, ALLOCATABLE :: rademis(:,:)
      REAL*4, ALLOCATABLE :: fract(:,:)
      REAL*4, ALLOCATABLE :: posmax(:,:)        ! Now 2D: (nstacks, 4)

      INTEGER, ALLOCATABLE :: nrads(:)

      CHARACTER*8, ALLOCATABLE :: radname(:,:)
      CHARACTER*1, ALLOCATABLE :: radtype(:,:)

c Declare other variables and calculated terms

      REAL alldose, idose
      INTEGER nstacks, n, m, i, j
      CHARACTER*200 title, company, facility, address, poc, phone, inline
      CHARACTER*8 NUC
      CHARACTER*1 radtype_upper

      integer date_time(8)
      character*10 b(3)

c Open input and output files 

      write(*,*)
      write(*,*)'COMPLY v2.0'
      write(*,*)
 
      open(11,file='input.level1.txt',status='old',form='formatted')
      open(12,file='output.level1.txt',status='unknown',form='formatted')
      open(13,file='comply171_original_possession_data.csv',status='old',form='formatted')

c Read input namelist file

       read(11,'(20x,a)') title
       read(11,'(20x,a)') company
       read(11,'(20x,a)') facility
       read(11,'(20x,a)') address
       read(11,'(20x,a)') poc
       read(11,'(20x,a)') phone
       read(11,'(20x,i3)') nstacks

c Allocate all arrays ONCE, before loop
       allocate(latitude(nstacks))
       allocate(longitude(nstacks))
       allocate(nrads(nstacks))
       allocate(radname(nstacks,100))         ! Assume max 100 radionuclides
       allocate(rademis(nstacks,100))
       allocate(fract(nstacks,100))
       allocate(radtype(nstacks,100))
       allocate(posmax(nstacks,4))

       do n = 1, nstacks 
       read(11,'(a)') 
       read(11,'(20x,f10.5)') latitude(n)
       read(11,'(20x,f10.5)') longitude(n)
       read(11,'(20x,i3)') nrads(n)

       do m = 1, nrads(n)
       read(11,'(20x,a)') inline
       read(inline,*) radname(n,m), radtype(n,m), rademis(n,m)
       enddo
       enddo

c Read in ancillary files - store in lookup table FIRST

       write(*,*)
       read(13,'(a)') inline  ! Skip header
       
c Create hash map: radionuclide name -> index (or use direct search if only 823 lines)
       do i = 1, 823
       read(13,'(a)') inline
       read(inline,*) NUC, (posmax(1,j),j=1,4)
       
c Now search only once per radionuclide, not per stack/radionuclide combo
       do n = 1, nstacks
        do m = 1, nrads(n)
        if(radname(n,m) .eq. NUC) then
         write(*,*)'Match:', NUC, radtype(n,m)
         
c Normalize case once
         radtype_upper = radtype(n,m)
         if(radtype_upper .eq. 'g') radtype_upper = 'G'
         if(radtype_upper .eq. 'l') radtype_upper = 'L'
         if(radtype_upper .eq. 's') radtype_upper = 'S'
         if(radtype_upper .eq. 'c') radtype_upper = 'C'
         
c Use SELECT CASE for efficiency
         select case(radtype_upper)
         case('G')
          fract(n,m) = rademis(n,m) / posmax(1,1)
         case('L')
          fract(n,m) = rademis(n,m) / posmax(1,2)
         case('S')
          fract(n,m) = rademis(n,m) / posmax(1,3)
         case('C')
          fract(n,m) = rademis(n,m) / posmax(1,4)
         case default
          write(*,*)'Unrecognized physical state:', radtype(n,m)
          fract(n,m) = 9999.
         end select
        endif
        enddo
       enddo
       enddo

c Calculate possession - avoid string I/O in loop
       alldose = 0.0
       idose = 0.0

       do n = 1, nstacks
       do m = 1, nrads(n)
       alldose = alldose + fract(n,m)
c Check for iodine directly with INDEX
       if(index(radname(n,m), 'I-') > 0) then
        idose = idose + fract(n,m)
       endif
       enddo
       enddo

       write(*,*)'Total possession (fraction):', alldose
       write(*,*)'Total iodine possession (fraction):', idose

c Write standardized output (unchanged)

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
       write(12,1008) n
       write(12,1006) latitude(n)
       write(12,1007) longitude(n)
       write(12,1013) nrads(n)
       write(12,2002)
       write(12,2003)

       do m = 1, nrads(n)
       write(12,1014) radname(n,m), radtype(n,m), rademis(n,m), fract(n,m)
       enddo
       enddo

       write(12,2002)
       write(12,2001)
       write(12,2001)
       write(12,2000) alldose
       write(12,1999) idose
       write(12,2001)
       write(12,2001)
       write(12,4000)

c Add time stamp

       write(*,*)
       call date_and_time(b(1), b(2), b(3), date_time)
       print *, 'Program execution date and time information:'
       print *, 'year=', date_time(1)
       print *, 'month_of_year=', date_time(2)
       print *, 'day_of_month=', date_time(3)
       print *, 'hour of day=', date_time(5)
       print *, 'minutes of hour=', date_time(6)
       print *, 'seconds of minute=', date_time(7)

       write(12,*)
       write(12,2100) date_time(2), date_time(3), date_time(1),
     &           date_time(5), date_time(6)

c Formatting for output file

 1000 FORMAT('Project:',a)
 1001 FORMAT('Company:',a)
 1002 FORMAT('Facility:',a)
 1003 FORMAT('Address:',a)
 1004 FORMAT('Point of Contact:',a)
 1005 FORMAT('Phone Number:',a)
 1006 FORMAT('Latitude (degrees):',f10.5)
 1007 FORMAT('Longitude (degrees):',f10.5)
 1008 FORMAT('Stack number:',i3)
 1013 FORMAT('Number of radionuclides:',I2)
 1014 FORMAT(a,',',a,',',f25.20,',',f25.20)
 1999 FORMAT('Total possession iodine fraction:',f10.2)
 2000 FORMAT('Total possession fraction:',f10.2)
 2001 FORMAT('-----------------------------------------------------------------------')
 2002 FORMAT(' ')
 2003 FORMAT('Radionuclide, Phase, Amount, Possession Fraction')
 2100 FORMAT('Date/Time (LST) stamp of code execution (MM/DD/YYYY HH:MM):',i2,'/',i2,'/',i4,' ',i2,':',i2)
 4000 FORMAT('COMPLY version 2.0 level 1 distributed through GitHub')

      END 
