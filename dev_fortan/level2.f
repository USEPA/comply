      PROGRAM COMPLY

c Declare dynamic arrays

      REAL*4, ALLOCATABLE :: releaseht(:)       !stack height (m)
      REAL*4, ALLOCATABLE :: bldheight(:)       !building height (m)
      REAL*4, ALLOCATABLE :: bldwidth(:)        !building width (m)
      REAL*4, ALLOCATABLE :: distance(:)        !distance between stack and receptor (m)
      REAL*4, ALLOCATABLE :: windspeed(:)       !wind speed (m/s)
      REAL*4, ALLOCATABLE :: latitude(:)        !latitude (decimal degrees)
      REAL*4, ALLOCATABLE :: longitude(:)       !longitude (decimal degrees)

      REAL*8, ALLOCATABLE :: rademis(:,:)       !radionuclide emissions (Ci/s) for each stack
      REAL*8, ALLOCATABLE :: conc(:,:)          !radionuclide air concentration (Bq/m3)
      REAL*8, ALLOCATABLE :: dose(:,:)          !radionuclide dose (mrem/yr)
      REAL*8, ALLOCATABLE :: wbfac(:,:)         !whole body dose factor from Oak Ridge table

      INTEGER, ALLOCATABLE :: nrads(:)          !number of emitted radionucludes at each stack

      CHARACTER*8, ALLOCATABLE :: radname(:,:)  !radionuclide names at each stack

c Declare other variables and calculated terms

      REAL FACTOR,area,xcrit,SIGMAZ,CAPSIG,pp,alldose,idose,dist,altconc
      REAL altdose

      INTEGER nstacks                           !number of unique stacks
      INTEGER N, M, I, J                        !loop counters

      CHARACTER*200 title,company,facility,address,poc,phone,inline
      CHARACTER*8 NUC
      CHARACTER*2 NN

      integer date_time(8)                      !for date stamp function
      character*10 b(3)                         !for date stamp function

c Open input and output files 

      write(*,*)
      write(*,*)'COMPLY v2.0 for internal testing only'
      write(*,*)
 
      open(11,file='input.level2.txt',status='old',form='formatted')
      open(12,file='output.level2.txt',status='unknown',form='formatted')
      open(13,file='diag.level2.txt',status='unknown',form='formatted')
      open(10,file='input.level2.table2.csv',status='unknown',form='formatted')

c Read input namelist file

       read(11,'(20x,a)') title
       read(11,'(20x,a)') company
       read(11,'(20x,a)') facility
       read(11,'(20x,a)') address
       read(11,'(20x,a)') poc
       read(11,'(20x,a)') phone

       read(11,'(20x,i3)') nstacks

       allocate(latitude(nstacks))
       allocate(longitude(nstacks))
       allocate(releaseht(nstacks))
       allocate(bldheight(nstacks))
       allocate(bldwidth(nstacks))
       allocate(distance(nstacks))
       allocate(windspeed(nstacks))
       allocate(nrads(nstacks))

       do n = 1, nstacks 

       read(11,'(20x,f10.5)') latitude(n)
       read(11,'(20x,f10.5)') longitude(n)
       read(11,'(20x,f10.5)') releaseht(n)
       read(11,'(20x,f10.5)') bldheight(n)
       read(11,'(20x,f10.5)') bldwidth(n)
       read(11,'(20x,f10.5)') distance(n)
       read(11,'(20x,f10.5)') windspeed(n)

       read(11,'(20x,i3)') nrads(n)

       allocate(radname(nstacks,nrads(n)))
       allocate(rademis(nstacks,nrads(n)))
       allocate(conc(nstacks,nrads(n)))
       allocate(dose(nstacks,nrads(n)))
       allocate(wbfac(nstacks,nrads(n)))

       do m = 1, nrads(n)
       read(11,'(20x,a)') inline
       read(inline,*) radname(n,m),rademis(n,m)
       enddo !end loop over reading in radionuclide emissions
   
       enddo !end loop over reading in stacks

c Read in ancillary files

       write(*,*)
       write(10,700)

       OPEN( 14, FILE='TABLE2.DAT', STATUS='OLD', IOSTAT=NERR,
     &         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=3 )

        do i = 1, 878
         READ( 14 , REC = i, err=999 ) NUC, FACTOR

c Make a text version of the originally distributed dose table
         write(10,701)  NUC, FACTOR
  700 FORMAT('Nuclide,Dose (mrem/yr)/(Ci/m3)')
  701 FORMAT(a,1(',',E10.3))

         do n = 1, nstacks
           do m = 1, nrads(n)
           if(radname(n,m).eq.NUC) then
           wbfac(n,m) = FACTOR
           write(*,*) 'Whole body dose factor from table:',NUC,FACTOR,wbfac(n,m)
           endif
           enddo !end loop over radionuclides
          enddo !end loop over stacks
        enddo !end loop over lines on the input file 
 999    continue

        write(*,*)

c Calculate air quality impacts and dose for each stack
c Emission rate assumption is Ci/s
c 0.25 term used in following equations is a default wind direction frequency fraction

       do n = 1, nstacks

       area = bldwidth(n) * bldheight(n)
       xcrit = 2.5 * sqrt(area)

       do m = 1, nrads(n)
       if(distance(n).lt.xcrit) then
       conc =   (0.25*rademis(n,m)) /
     &          (3.14*windspeed(n)*1.0*amin1(bldwidth(n),bldheight(n))) 
       else
       SIGMAZ = (0.06 * distance(n)) / SQRT( 1.0 + 0.0015 * distance(n))
       CAPSIG = SQRT( SIGMAZ ** 2 + area / 3.14 )
       pp = 2.032 / ( distance(n) * CAPSIG )
       conc(n,m) = (0.25*rademis(n,m)*pp) / windspeed(n)
       endif

       dose(n,m) = wbfac(n,m) * conc(n,m)

       write(*,*)'Concentration and dose:',n,radname(n,m),rademis(n,m),
     & distance(n),conc(n,m),dose(n,m)

       enddo !end loop over radionuclides
       enddo !end loop over stacks

       write(*,*)

c Calculate total dose

       alldose = 0.
       idose = 0.

       do n = 1, nstacks
       do m = 1, nrads(n)
       alldose = alldose + dose(n,m)
         read(radname(n,m),'(a2)') NN
         if(NN.eq.'I-') then
         idose = idose + dose(n,m)
         endif
       enddo !end loop over radionuclides
       enddo !end loop over stacks

       write(*,*)'Total dose (mrem/yr):',alldose
       write(*,*)'Total dose limit is 10.0 mrem/yr'
       write(*,*)'Total iodine dose (mrem/yr):',idose
       write(*,*)'Iodine limit is 3.0 mrem/yr'

c Write stardardized output

       write(12,1000) title
       write(12,1001) company
       write(12,1002) facility 
       write(12,1003) address
       write(12,1004) poc
       write(12,1005) phone

       do n = 1, nstacks
       write(12,1006) latitude(n)
       write(12,1007) longitude(n)
       write(12,1008) releaseht(n)
       write(12,1009) bldheight(n)
       write(12,1010) bldwidth(n)
       write(12,1011) distance(n)
       write(12,1012) windspeed(n)

       write(12,1013) nrads(n)

       do m=1, nrads(n)
       write(12,1014) radname(n,m),rademis(n,m),conc(n,m),dose(n,m)
       enddo !end loop over radionuclides

       enddo !end loop over stacks

       write(12,2000) alldose
       write(12,2001) idose

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


c Formatting for output file

 1000 FORMAT('Project:',a)
 1001 FORMAT('Company:',a)
 1002 FORMAT('Facility:',a)
 1003 FORMAT('Address:',a)
 1004 FORMAT('Point of Contact:',a)
 1005 FORMAT('Phone Number:',a)
 1006 FORMAT('Latitude (degrees):',f10.5)
 1007 FORMAT('Longitude (degrees):',f10.5)
 1008 FORMAT('Release height (m):',f10.1)
 1009 FORMAT('Building height (m):',f10.1)
 1010 FORMAT('Building width (m):',f10.1)
 1011 FORMAT('Distance from source to receptor (m):',f10.1)
 1012 FORMAT('Assumed wind speed (m/s):',f10.1)
 1013 FORMAT('Number of radionuclides:',I2)
 1014 FORMAT(a,',',f25.20,',',f25.20,',',f10.6)
 2000 FORMAT('Total effective dose equivalent (mrem/yr):',f10.2)
 2001 FORMAT('Total effective dose equivalent iodine only (mrem/yr):',f10.2)
 2100 FORMAT('Date/Time (LST) stamp of simulation (MM/DD/YYYY HH:MM):',i2,'/',i2,'/',i4,' ',i2,':',i2)

c Extra output for QA purposes only

       write(13,3001) !write header to file
 
       do n = 1, nstacks

       area = bldwidth(n) * bldheight(n)
       xcrit = 2.5 * sqrt(area)

       do m = 1, nrads(n)

       do dist = 5, 2000, 5

       SIGMAZ = (0.06 * dist) / SQRT( 1.0 + 0.0015 * dist)
       CAPSIG = SQRT( SIGMAZ ** 2 + area / 3.14 )
       pp = 2.032 / ( dist * CAPSIG )

       if(dist.lt.xcrit) then
       conc =   (0.25*rademis(n,m)) /
     &          (3.14*windspeed(n)*1.0*amin1(bldwidth(n),bldheight(n)))
       else
       conc(n,m) = (0.25*rademis(n,m)*pp) / windspeed(n)
       endif
        
       altconc = (0.25*rademis(n,m)*pp) / windspeed(n) 

       dose(n,m) = wbfac(n,m) * conc(n,m)
       altdose = wbfac(n,m) * altconc

       write(13,3000) n,radname(n,m),dist,
     &     conc(n,m),dose(n,m),altconc,altdose

       enddo !end loop over radionuclides
       enddo !end loop over distances between source and receptor
       enddo !end loop over stacks

 3000  FORMAT(i2,',',a,',',f6.1,',',f25.20,',',f25.20,',',f25.20,',',f25.20)
 3001  FORMAT('STACK,RADNAME,DISTANCE,CONC,DOSE,ALTCONC,ALTDOSE')

c End of program

      END 

