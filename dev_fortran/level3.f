      PROGRAM COMPLY

c Declare dynamic arrays

      REAL*4, ALLOCATABLE :: releaseht(:)           !stack height (m)      
      REAL*4, ALLOCATABLE :: diameter(:)            !stack diameter (m)
      REAL*4, ALLOCATABLE :: exitvelocity(:)        !exit velocity ()
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
      REAL*8, ALLOCATABLE :: conc(:,:)          !radionuclide air concentration (Bq/m3)
      REAL*8, ALLOCATABLE :: dose(:,:)          !radionuclide dose (mrem/yr)
      REAL*8, ALLOCATABLE :: wbfac(:,:)         !whole body dose factor from Oak Ridge table
      REAL*8, ALLOCATABLE :: vegfac(:,:)        !vegetable dose factor from Oak Ridge table
      REAL*8, ALLOCATABLE :: cowfac(:,:)        !cow dose factor from Oak Ridge table
      REAL*8, ALLOCATABLE :: vegdose(:,:)       !radionuclide dose by vegetables (mrem/yr)
      REAL*8, ALLOCATABLE :: cowdose(:,:)       !radionuclide dose by meat and milk (mrem/yr)

      INTEGER, ALLOCATABLE :: nrads(:)          !number of emitted radionucludes at each stack

      CHARACTER*8, ALLOCATABLE :: radname(:,:)  !radionuclide names at each stack
      CHARACTER*3,  ALLOCATABLE :: verticalmomentum(:)    !vertical momentum (YES or NO)

c Declare other variables and calculated terms

      REAL*8 tempconc,tempdose,tempvegconc,tempvegdose,tempcowconc,tempcowdose
      REAL FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF,FACTORG,FACTORH,FACTORI
      REAL FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN
      REAL area,xcrit,SIGMAZ,CAPSIG,pp,alldose,idose,dist,altconc

      REAL altdose,totvegdose,itotvegdose,totcowdose,itotcowdose
      REAL RHOAIR,RHOSTK,STKVEL,R,FM,DELHMO,DELHBU,HEFF,FB,zz,INHFAC,TINIT,T,EQILOUT,EQILIN,INHDOS,FINDOOR

      INTEGER nstacks                           !number of unique stacks
      INTEGER N, M, I, J                        !loop counters
      INTEGER receptor

      CHARACTER*200 title,company,facility,address,poc,phone,inline
      CHARACTER*8 NUC
      CHARACTER*2 NN

      integer date_time(8)                      !for date stamp function
      character*10 b(3)                         !for date stamp function

c Open input and output files 

      write(*,*)
      write(*,*)'COMPLY v2.0 for internal testing only'
      write(*,*)
 
      open(11,file='input.level3.txt',status='old',form='formatted')
      open(12,file='output.level3.txt',status='unknown',form='formatted')
      open(13,file='diag.level3.txt',status='unknown',form='formatted')

      open(9,file='comply171_original_input_data.csv',status='old',form='formatted')

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
       allocate(exitvelocity(nstacks))
       allocate(diameter(nstacks))
       allocate(stacktemp(nstacks))
       allocate(airtemp(nstacks))
       allocate(verticalmomentum(nstacks))
       allocate(bldheight(nstacks))
       allocate(bldwidth(nstacks))
       allocate(distance(nstacks))
       allocate(distanceveg(nstacks))
       allocate(distancecow(nstacks))
       allocate(windspeed(nstacks))
       allocate(nrads(nstacks))

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
       read(11,'(20x,f10.5)') windspeed(n)

       read(11,'(20x,i3)') nrads(n)

       if(n.eq.1) then
       allocate(radname(nstacks,nrads(n)))
       allocate(rademis(nstacks,nrads(n)))
       allocate(conc(nstacks,nrads(n)))
       allocate(dose(nstacks,nrads(n)))
       allocate(wbfac(nstacks,nrads(n)))
       allocate(vegfac(nstacks,nrads(n)))
       allocate(cowfac(nstacks,nrads(n)))
       allocate(vegdose(nstacks,nrads(n)))
       allocate(cowdose(nstacks,nrads(n)))
       endif

       do m = 1, nrads(n)
       read(11,'(20x,a)') inline
       read(inline,*) radname(n,m),rademis(n,m)
       enddo !end loop over reading in radionuclide emissions
   
       enddo !end loop over reading in stacks

c Read in ancillary files

       write(*,*)
       rewind(9)
       read(9,'(a)') inline

          do jj = 1, 826
          read(9,'(a)') inline
          read(inline,*) NUC,FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF,FACTORG,
     &      FACTORH,FACTORI,FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN

          do n = 1, nstacks
          do m = 1, nrads(n)
           if(NUC.eq.radname(n,m)) then
            wbfac(n,m) = FACTORC
            vegfac(n,m) = FACTORD
            cowfac(n,m) = FACTORE
            write(*,*) 'Whole body minus ingestion pathway dose factor:',NUC,FACTORC
            write(*,*) 'Vegetable dose factor from table:',NUC,FACTORD
            write(*,*) 'Cow dose factor from table:',NUC,FACTORE
            write(*,*)
           endif
          enddo
          enddo

          enddo

        write(*,*)

c Calculate air quality impacts and dose for each stack
c Emission rate assumption is Ci/s
c 0.25 term used in following equations is a default wind direction frequency fraction

       write(13,3001) !write header to file

       alldose = 0.
       idose = 0.
       totvegdose = 0.
       itotvegdose = 0.
       totcowdose = 0.
       itotcowdose = 0.


       do n = 1, nstacks

       area = bldwidth(n) * bldheight(n)
       xcrit = 2.5 * sqrt(area)

       do m = 1, nrads(n)

       do receptor = 10, 2000, 10

       dist = float(receptor)

       if(exitvelocity(n).le.0.) then
       exitvelocity(n) = 0.3
       endif

       HEFF = releaseht(n)

       IF ( releaseht(n) .GT. 2.5*bldheight(n) ) THEN

! Temperature
        RHOAIR = 0.0735 * ( 80.0 + 460.0 ) / ( airtemp(n) + 460.0 )
        RHOAIR = RHOAIR * 454.0 / ( 0.3048 ** 3 )
        RHOSTK = RHOAIR * ( airtemp(n) + 460.0 ) / ( stacktemp(n) + 460.0 )

! Momentum
        AREA   = 3.14 * ( diameter(n) ** 2.0 ) / 4.0
        STKVEL = exitvelocity(n) / AREA
        R      = STKVEL / windspeed(n)
        FM     = ( RHOSTK / RHOAIR ) * STKVEL * exitvelocity(n) / 3.14
        TERM1  = 75.0 * FM / ( ( ( 1.0 + 3.0 / R ) * windspeed(n) ) ** 2.0 )
        TERM2  = releaseht(n) + 7.79 * SQRT( FM ) / ( ( 1.0 + 3.0 / R ) * windspeed(n) ) !is this stack height or HEFF??

! If the site has a vertical stack, then the momentum rise value is zero, otherwise compute the value.
        IF ( verticalmomentum(n).eq.'YES' ) THEN
         DELHMO = 0.93 * ( TERM1 ** ( 3.0 / 7.0 ) ) *
     &        ( TERM2 ** ( 1.0 / 7.0 ) )
        ELSE
         DELHMO = 0.0
        ENDIF

! Buoyant
        FB = 9.81 * ( RHOAIR - RHOSTK ) * exitvelocity(n) / ( 3.14 * RHOAIR )
        DELHBU = 23.67 * ( ( FB / ( windspeed(n) ** 3.0 ) ) ** 0.6 ) *
     &           ( releaseht(n) + 195.1 * FB / ( windspeed(n) ** 3.0 ) ) ** 0.4 !stack height or HEFF??

! Use maximum of buoyant or momentum rise
       HEFF = releaseht(n) + AMAX1( DELHMO, DELHBU )
!!       write(*,*)'WW',releaseht(n),HEFF,DELHMO, DELHBU

       ENDIF !end condition where H > 0


! Calculate PP
       SIGMAZ = (0.06 * dist) / SQRT( 1.0 + 0.0015 * dist)
       CAPSIG = SQRT( SIGMAZ ** 2 + area / 3.14 )
       zz = 0.5 * ( ( HEFF / CAPSIG ) ** 2 )
        if(zz.gt.0.001.and.releaseht(n).gt.2.5*bldheight(n)) then
        pp = 2.032 * EXP( -zz ) / ( dist * CAPSIG )
        else
        pp = 2.032 / ( dist * CAPSIG )
        endif


       if(dist.lt.xcrit) then
       tempconc =   (0.25*rademis(n,m)) /
     &          (3.14*windspeed(n)*1.0*amin1(bldwidth(n),bldheight(n)))
       else
       tempconc = (0.25*rademis(n,m)*pp) / windspeed(n)
       endif

       altconc = (0.25*rademis(n,m)*pp) / windspeed(n)

       tempdose = wbfac(n,m) * tempconc
       tempvegdose = vegfac(n,m) * tempconc
       tempcowdose = cowfac(n,m) * tempconc
       altdose = wbfac(n,m) * altconc
       altdoseb = vegfac(n,m) * altconc
       altdosec = cowfac(n,m) * altconc

       write(13,3000) n,radname(n,m),dist,
     &     tempconc,tempdose,
     &     tempvegdose,tempcowdose,
     &     altconc,altdose,
     &     altdoseb,altdosec

       if(dist.eq.distancecow(n)) then
       totcowdose = totcowdose + tempcowdose
       cowdose(n,m) = tempcowdose
         read(radname(n,m),'(a2)') NN
         if(NN.eq.'I-') then
         itotcowdose = itotcowdose + tempcowdose
         endif
       endif

       if(dist.eq.distanceveg(n)) then
       totvegdose = totvegdose + tempvegdose
       vegdose(n,m) = tempvegdose
         read(radname(n,m),'(a2)') NN
         if(NN.eq.'I-') then
         itotvegdose = itotvegdose + tempvegdose
         endif
       endif

       if(dist.eq.distance(n)) then
       alldose = alldose + tempdose
       dose(n,m) = tempdose
         read(radname(n,m),'(a2)') NN
         if(NN.eq.'I-') then
         idose = idose + tempdose
         endif
       endif

       enddo !end loop over radionuclides
       enddo !end loop over distances between source and receptor
       enddo !end loop over stacks

c Calculate total dose

       write(*,*)
       write(*,*)'Total dose (mrem/yr):',alldose+totvegdose+totcowdose
       write(*,*)'Total dose limit is 10.0 mrem/yr'
       write(*,*)'Total iodine dose (mrem/yr):',idose+itotvegdose+itotcowdose
       write(*,*)'Iodine limit is 3.0 mrem/yr'
       write(*,*)
       write(*,*)'Total dose vegetable fraction:',totvegdose/(alldose+totvegdose+totcowdose)
       write(*,*)'Total dose meat/milk fraction:',totcowdose/(alldose+totvegdose+totcowdose)
       write(*,*)'Total dose other pathways fraction:',alldose/(alldose+totvegdose+totcowdose)
       write(*,*)

c Write stardardized output

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
       write(12,1015) distance(n)
       write(12,1016) distanceveg(n)
       write(12,1017) distancecow(n)
       write(12,1018) windspeed(n)
       do m = 1, nrads(n)
       write(12,2002)
       write(12,1019) radname(n,m)
       write(12,1020) rademis(n,m)
       write(12,1021) (dose(n,m)+vegdose(n,m)+cowdose(n,m))
       write(12,1022) dose(n,m)
       write(12,1023) vegdose(n,m)
       write(12,1024) cowdose(n,m)
       enddo
       enddo !end loop over stacks

       write(12,2002)
       write(12,2001)
       write(12,2001)
       write(12,2000) alldose+totvegdose+totcowdose
       write(12,2002)
       write(12,2005) totvegdose/(alldose+totvegdose+totcowdose)
       write(12,2006) totcowdose/(alldose+totvegdose+totcowdose)
       write(12,2007) alldose/(alldose+totvegdose+totcowdose)

       write(12,2001)
       write(12,2001)
       write(12,4000)


c Add time stamp of model simulation to output file

       call date_and_time(b(1), b(2), b(3), date_time)

       print *,'Program execution date and time information:'
       print *,'year=',date_time(1)
       print *,'month_of_year=',date_time(2)
       print *,'day_of_month=',date_time(3)
!       print *,'time difference in minutes=',date_time(4)
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
 1006 FORMAT('Stack number:',i3)
 1007 FORMAT('Latitude (degrees):',f10.5)
 1008 FORMAT('Longitude (degrees):',f10.5)
 1009 FORMAT('Release height (m):',f10.1)
 1010 FORMAT('Stack diameter (m):',f10.1)
 1011 FORMAT('Exit velocity (m3/s):',f10.1)
 1012 FORMAT('Stack temperature (F):',f10.1)
 1013 FORMAT('Air temperature (F):',f10.1)
 1014 FORMAT('Vertical momentum (YES/NO):',a)
 1015 FORMAT('Distance from source to receptor (m):',f10.1)
 1016 FORMAT('Distance from source to vegetable receptor (m):',f10.1)
 1017 FORMAT('Distance from source to meat/milk receptor (m):',f10.1)
 1018 FORMAT('Assumed wind speed (m/s) when no wind rose provided:',f10.1)
 1019 FORMAT('Radionuclide:',a10)
 1020 FORMAT('Emission rate (Ci/s):',f25.20)
 1021 FORMAT('Effective dose equivalent all pathways (mrem/yr):',f25.20)
 1022 FORMAT('Effective dose equivalent non-ingestion (mrem/yr):',f25.20)
 1023 FORMAT('Effective dose equivalent vegetables (mrem/yr):',f25.20)
 1024 FORMAT('Effective dose equivalent meat/milk (mrem/yr):',f25.20)

 2000 FORMAT('Total effective dose equivalent all stacks (mrem/yr):',f10.2)
 2001 FORMAT('-----------------------------------------------------------------------')
 2002 FORMAT(' ')
 2005 FORMAT('Vegetable fraction of total effective dose equivalent:',f5.2)
 2006 FORMAT('Meat/milk fraction of total effective dose equivalent:',f5.2)
 2007 FORMAT('All other pathway fraction of total effective dose equivalent:',f5.2)


 2100 FORMAT('Date/Time (LST) stamp of code execution (MM/DD/YYYY HH:MM):',i2,'/',i2,'/',i4,' ',i2,':',i2)

 3000 FORMAT(i2,',',a,',',f6.1,8(',',f25.20))
 3001 FORMAT('STACK,RADNAME,DISTANCE,CONC,WBDOSE,VEGDOSE,COWDOSE,ALTCONC,ALTWB,ALTVEG,ALTCOW')

 4000 FORMAT('COMPLY version 2.0 level 3 distributed through GitHub') 

c End of program

      END 

