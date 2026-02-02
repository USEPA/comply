      PROGRAM COMPLY

      IMPLICIT NONE
c Declare dynamic arrays

      REAL*4, ALLOCATABLE :: releaseht(:)           !stack height (m)
      REAL*4, ALLOCATABLE :: diameter(:)            !stack diameter (m)
      REAL*4, ALLOCATABLE :: exitvelocity(:)        !exit velocity ()
      REAL*4, ALLOCATABLE :: distance(:)            !distance between stack and receptor (m)
      REAL*4, ALLOCATABLE :: airtemp(:)             !air temperature (F)
      REAL*4, ALLOCATABLE :: stacktemp(:)           !stack temperature (F)
      REAL*4, ALLOCATABLE :: windspeed(:)           !wind speed (m/s)
      REAL*4, ALLOCATABLE :: latitude(:)            !latitude (decimal degrees)
      REAL*4, ALLOCATABLE :: longitude(:)           !longitude (decimal degrees)
      REAL*4, ALLOCATABLE :: pctcalm(:)             !percentage of calm days at wind data location

      REAL*8, ALLOCATABLE :: rademis(:,:)           !radionuclide emissions (Ci/s) for each stack
      REAL*8, ALLOCATABLE :: conc(:,:)              !radionuclide air concentration (Bq/m3)
      REAL*8, ALLOCATABLE :: dose(:,:)              !radionuclide dose (mrem/yr)

      REAL*8, ALLOCATABLE :: ihdose(:)
      REAL*8, ALLOCATABLE :: totihdose(:,:)

      REAL*8, ALLOCATABLE :: windbinfreq(:,:)       !Sector (wind dir) frequency from that direction
      REAL*8, ALLOCATABLE :: windbinspeed(:,:)      !Sector (wind dir) average speed from that direction (m/s)
      REAL*8, ALLOCATABLE :: windbinsrdist(:,:)     !Sector (wind dir) source-receptor distance (m)

      INTEGER, ALLOCATABLE :: nrads(:)              !number of emitted radionucludes at each stack
      INTEGER, ALLOCATABLE :: nwindbins(:)          !number of wind direction bins at each stack

      CHARACTER*8,  ALLOCATABLE :: radname(:,:)     !radionuclide names at each stack
      CHARACTER*3,  ALLOCATABLE :: windbinname(:,:) !wind direction sector 
      CHARACTER*3,  ALLOCATABLE :: sectorname(:,:)  !sector name for receptors
      CHARACTER*3,  ALLOCATABLE :: srcsector(:)     !
      CHARACTER*3,  ALLOCATABLE :: verticalmomentum(:)    !vertical momentum (YES or NO)
      CHARACTER*3,  ALLOCATABLE :: meatsector(:)    !Sector (wind dir) location of meat
      CHARACTER*8,  ALLOCATABLE :: decaychain(:)    !Decay chain list
      CHARACTER*8,  ALLOCATABLE :: decaychainname(:)!Decay chain radionuclide names

      CHARACTER*200,ALLOCATABLE :: windsource(:)    !Source of wind data
      CHARACTER*200,ALLOCATABLE :: winddates(:)     !Date coverage of wind data
      CHARACTER*200,ALLOCATABLE :: windloc(:)       !Location of wind data
      CHARACTER*200,ALLOCATABLE :: winddist(:)      !Distance from wind data to source

! Declare other variables and calculated terms

      REAL FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF,FACTORG,FACTORH,FACTORI
      REAL FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN
      REAL area,xcrit,SIGMAZ,CAPSIG,pp,alldose,idose,dist,altconc
      REAL tempdose,maxdose,maxdist,totmilkdose,totvegdose,itotvegdose,totmeatdose,itotcowdose
      REAL capp,capn0,outprod,inprod,z1,term1,term2,sum
      REAL transittime                              !transit time of air from emission release to receptor (s)
      REAL fdux                                     !not sure what this is-units are 1/(m2/s) or s/m2
      REAL wetdepadj,wetdepadjbad                                !wet deposition adjustment factor
      REAL xtime
      REAL xremove
      REAL ingestion, immersion, inhalation, surface
      REAL RHOAIR,RHOSTK,STKVEL,R,FM,DELHMO,DELHBU,HEFF,FB,zz,INHFAC,TINIT,T,EQILOUT,EQILIN,INHDOS,FINDOOR

      INTEGER nstacks                               !number of unique stacks
      INTEGER receptor, N, M, I, J, K, L, II, JJ, P, Q                        !loop counters
      INTEGER nchain                                !number of nuclides in the decay series
      INTEGER maxsect

      CHARACTER*200 title,company,facility,address,poc,phone,inline
      CHARACTER*8 NUC
      CHARACTER*2 NN

      integer date_time(8)                      !for date stamp function
      character*10 b(3)                         !for date stamp function


! Open input and output files 

      write(*,*)
      write(*,*)'COMPLYr v2.0 for internal testing only'
      write(*,*)
 
      open(11,file='input.level4r.txt',status='old',form='formatted')
      open(12,file='output.level4r.txt',status='unknown',form='formatted')
      open(13,file='diag.level4r.txt',status='unknown',form='formatted')

!      open(9,file='comply171_original_input_data.csv',status='old',form='formatted')
!      open(10,file='output.level4.tree.csv',status='old',form='formatted')

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
       allocate(distance(nstacks))
       allocate(stacktemp(nstacks))
       allocate(airtemp(nstacks))
       allocate(verticalmomentum(nstacks))
       allocate(srcsector(nstacks))
       allocate(windspeed(nstacks))
       allocate(nrads(nstacks))
       allocate(windsource(nstacks))
       allocate(winddates(nstacks))
       allocate(windloc(nstacks))
       allocate(winddist(nstacks))
       allocate(pctcalm(nstacks))
       allocate(nwindbins(nstacks))

       do n = 1, nstacks 

       read(11,'(20x,f10.5)') latitude(n)
       read(11,'(20x,f10.5)') longitude(n)
       read(11,'(20x,f10.5)') releaseht(n)
       read(11,'(20x,f10.5)') exitvelocity(n)
       read(11,'(20x,f10.5)') diameter(n)
       read(11,'(20x,f10.5)') stacktemp(n)
       read(11,'(20x,f10.5)') airtemp(n)
       read(11,'(20x,a)') verticalmomentum(n)
       read(11,'(20x,f10.5)') distance(n)

       read(11,'(20x,a)') srcsector(n)
       read(11,'(20x,f10.5)') windspeed(n)

       read(11,'(20x,i3)') nrads(n)

       if(n.eq.1) then
       allocate(radname(nstacks,nrads(n)))
       allocate(rademis(nstacks,nrads(n)))
       allocate(conc(nstacks,nrads(n)))
       allocate(dose(nstacks,nrads(n)))
       allocate(totihdose(nstacks,nrads(n)))
       endif

       do m = 1, nrads(n)
       read(11,'(20x,a)') inline
       read(inline,*) rademis(n,m)
       enddo !end loop over reading in radionuclide emissions

       read(11,'(20x,a)') windsource(n)
       read(11,'(20x,a)') winddates(n)
       read(11,'(20x,a)') windloc(n)
       read(11,'(20x,a)') winddist(n)

       read(11,'(20x,f10.5)') pctcalm(n)
       read(11,'(20x,i3)') nwindbins(n)

       if(n.eq.1) then
       allocate(windbinname(nstacks,nwindbins(n)))
       allocate(windbinfreq(nstacks,nwindbins(n)))
       allocate(windbinspeed(nstacks,nwindbins(n)))
       allocate(windbinsrdist(nstacks,nwindbins(n)))
       allocate(sectorname(nstacks,nwindbins(n)))
       endif

       do m = 1, nwindbins(n)
       read(11,'(20x,a)') inline
       read(inline,*) windbinname(n,m),windbinfreq(n,m),windbinspeed(n,m),windbinsrdist(n,m),sectorname(n,m)
       enddo !end loop over reading in wind data and source-receptor distances

       enddo !end loop over reading in stacks



!Calculate air concentrations
!Emission rate assumption is Ci/s

       write(13,3001) !write header to file

       alldose = 0.

       do n = 1, nstacks

       do m = 1, nrads(n)

       maxdose = 0.0
       inhalation = 0.
       immersion = 0.
       ingestion = 0.
       totmilkdose = 0.
       totvegdose = 0.
       totmeatdose = 0.

       do j = 1, nwindbins(n)

       do receptor = 100, 40000, 100

       dist = float(receptor)

       if(exitvelocity(n).le.0.) then
       exitvelocity(n) = 0.3 
       endif

       HEFF = releaseht(n)

       IF ( releaseht(n) .GT. 0.0 ) THEN

! Temperature
        RHOAIR = 0.0735 * ( 80.0 + 460.0 ) / ( airtemp(n) + 460.0 )
        RHOAIR = RHOAIR * 454.0 / ( 0.3048 ** 3 )
        RHOSTK = RHOAIR * ( airtemp(n) + 460.0 ) / ( stacktemp(n) + 460.0 )

! Momentum
        AREA   = 3.14 * ( diameter(n) ** 2.0 ) / 4.0
        STKVEL = exitvelocity(n) / AREA
        R      = STKVEL / windbinspeed(n,j)
        FM     = ( RHOSTK / RHOAIR ) * STKVEL * exitvelocity(n) / 3.14
        TERM1  = 75.0 * FM / ( ( ( 1.0 + 3.0 / R ) * windbinspeed(n,j) ) ** 2.0 )
        TERM2  = releaseht(n) + 7.79 * SQRT( FM ) / ( ( 1.0 + 3.0 / R ) * windbinspeed(n,j) ) !is this stack height or HEFF??

!!        write(*,*)'YY! ',verticalmomentum(n),diameter(n),TERM1,TERM2,R,FM,AREA,STKVEL

! If the site has a vertical stack, then the momentum rise value is zero, otherwise compute the value.
        IF ( verticalmomentum(n).eq.'YES' ) THEN
         DELHMO = 0.93 * ( TERM1 ** ( 3.0 / 7.0 ) ) *
     &        ( TERM2 ** ( 1.0 / 7.0 ) )
        ELSE
         DELHMO = 0.0
        ENDIF

! Buoyant
        FB = 9.81 * ( RHOAIR - RHOSTK ) * exitvelocity(n) / ( 3.14 * RHOAIR )
        DELHBU = 23.67 * ( ( FB / ( windbinspeed(n,j) ** 3.0 ) ) ** 0.6 ) *
     &           ( releaseht(n) + 195.1 * FB / ( windbinspeed(n,j) ** 3.0 ) ) ** 0.4 !stack height or HEFF??

! Use maximum of buoyant or momentum rise
       HEFF = releaseht(n) + AMAX1( DELHMO, DELHBU )
!!       write(*,*)'WW',releaseht(n),HEFF,DELHMO, DELHBU

       ENDIF !end condition where H > 0

! Subroutine MAXP finds XMAX (point of maximum ground-level concentration) and PMAX (the maximum value of P).
!     XMAX = 0.0
!      IF ( HEFF .GT. 0.0 ) THEN
!        CALL MAXP( HEFF, PMAX, XMAX )
!      ENDIF

!   If XX < XMAX, then the maximum value of P is used. (See Figure 3 of NCRP)
!      IF ( XX .LT. XMAX ) THEN
!        P = PMAX
!      ELSE
!        P = PP( XX, HEFF, 0.0 )
!      ENDIF

! Section A.3.1, Equation A-6
!      C = F * Q * P / U

!Section A.3.1, Equation A-6 (except Q, done in LEVEL4WR)
!      C = FF * P / UU


! Calculate pp
       SIGMAZ = (0.06 * dist) / SQRT( 1.0 + 0.0015 * dist)
       !CAPSIG = SQRT( SIGMAZ ** 2 + area / 3.14159265 ) 
       CAPSIG = SQRT( SIGMAZ ** 2 + 0.0 / 3.14159265 ) !area seems to be set to 0. in complyr
       zz = 0.5 * ( ( HEFF / CAPSIG ) ** 2 )
        if(zz.gt.0.001) then
        pp = 2.032 * EXP( -zz ) / ( dist * CAPSIG )
        else
        pp = 2.032 / ( dist * CAPSIG )
        endif


       conc(n,m) = ((windbinfreq(n,j)/(1.0-pctcalm(n))) *rademis(n,m)*pp) / windbinspeed(n,j)

       if(conc(n,m).lt. 1.0E-32 ) then
         conc(n,m) = 1.0E-32
       endif

       altconc = (windbinfreq(n,j)*rademis(n,m)*pp) / windbinspeed(n,j) !just for QA and not being used for anything else here

       !emissions should be in Ci/s
       !concentration should be in Ci/m3


!Calculate inhalation dose

!   Set value of dose factor for radon + dtr inhalation
!    NOTE-- this is a temporary value (not sure why this is "temporary")

!                   sv/wlm    rem/sv   mrem/rem
            INHFAC = .0092 *   100.   * 1000.

!   TINIT is time to attain 0.2 equil fract (outlet of mine)
            TINIT = 11.8

!   Calculate equil fract during transit

!           min =  min     m    m/s  s/min 
            T = TINIT + dist / ( windbinspeed(n,j) * 60. )

!   Calculate outdoor equil fract
          EQILOUT = 1. - .0479 * EXP(-T / 4.39 ) - 2.1963 * EXP(-T / 38.6 )
     &         +  1.2442 * EXP(-T / 28.4 )

            IF( EQILOUT .GT. 1.0 ) EQILOUT = 1.0

!   Calculate indoor equil fract

          EQILIN =  0.35 * ( 1. + EQILOUT )

            IF( EQILIN .GT. .7 ) EQILIN = .7

!   Calculate the inhalation dose
!      Ci/m**3   pCi/Ci    m**3/liter   WL-liter/pCi  WLmonths/year  mrem/WLM     = mrem/year
          INHDOS =    conc(n,m)  *  1.E12  *  (1./1000.)  * (1./100.)  *
     &       (8760./170.)   *  INHFAC

!   FINDOOR is the percentage of time the receptor is indoors.
      FINDOOR = .75

!   EQILIN is the indoor equilibrium fraction
!   (1.0 - FINDOOR) is the percentage of time the receptor is outdoors.

!   EQILOUT is the outdoor equilibrium fraction

      totihdose(n,m) = INHDOS * ( FINDOOR * EQILIN + (1.0 - FINDOOR) * EQILOUT )

! Save dose information at key input receptor locations for report

       if(dist.eq.windbinsrdist(n,j)) then 
       write(*,*)'Air concentration:',conc(n,m),sectorname(n,j),dist
            tempdose = totihdose(n,m)
        if(tempdose.gt.maxdose) then 
            maxdose = tempdose
            maxsect = j
            maxdist = dist
            inhalation = totihdose(n,m)
        endif
       endif

       write(13,3000) n,dist,conc(n,m),totihdose(n,m)

       enddo !end loop over distances between source and receptor (receptor)
       enddo !end loop over number of wind direction bins (j)

       alldose = alldose + inhalation

       write(*,*)'Stack ',n,' max dose (',maxdose,') receptor at ',windbinsrdist(n,maxsect),sectorname(n,maxsect)
       write(*,*)
       write(*,2999)'Inhalation dose is',inhalation,' from stack',n
       write(*,*)'--------------------------------------------------------------------------' 
       enddo !end loop over radionuclides (m)

       enddo !end loop over stacks (n)

       write(*,*)
       write(*,4000)'Total dose (mrem/yr):',alldose
       write(*,*)'Total dose limit is 10.0 mrem/yr'
       write(*,*)

!       enddo !end loop over stacks (n)

 2999  FORMAT(a,f20.15,a,i3)
 3000  FORMAT(i2,',',f6.1,',',f30.20,6(',',f20.10))
 3001  FORMAT('STACK,RADNAME,DISTANCE,CONC,IMMERSION,INHALATION,SURFACE,MEAT,MILK,VEG')
 4000  FORMAT(a,f10.4)

c Write stardardized output

       write(12,1000) title
       write(12,1001) company
       write(12,1002) facility 
       write(12,1003) address
       write(12,1004) poc
       write(12,1005) phone

       write(12,2001)
       do n = 1, nstacks
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
       write(12,1016) windspeed(n)
       write(12,1017) rademis(n,1)
       write(12,1018) totihdose(n,1)
       write(12,2001)
       enddo !end loop over stacks

       write(12,2001)
       write(12,2000) alldose
       write(12,2001)
       write(12,2001)

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
 1006 FORMAT('Stack number:',i3)
 1007 FORMAT('Latitude (degrees):',f10.5)
 1008 FORMAT('Longitude (degrees):',f10.5)
 1009 FORMAT('Release height (m):',f10.1)
 1010 FORMAT('Stack diameter (m):',f10.1)
 1011 FORMAT('Exit velocity (m3/s):',f10.1)
 1012 FORMAT('Stack temperature (F):',f10.1)
 1013 FORMAT('Air temperature (F):',f10.1)
 1014 FORMAT('Vertial momentum (YES/NO):',a)
 1015 FORMAT('Distance from source to receptor (m):',f10.1)
 1016 FORMAT('Assumed wind speed (m/s) when no wind rose provided:',f10.1)
 1017 FORMAT('Emission rate (Ci/s):',f25.20)
 1018 FORMAT('Effective dose equivalent (mrem/yr):',f25.20)
 2000 FORMAT('Total effective dose equivalent all stacks (mrem/yr):',f10.2)
 2001 FORMAT('-----------------------------------------------------------------------')
 2100 FORMAT('Date/Time (LST) stamp of simulation (MM/DD/YYYY HH:MM):',i2,'/',i2,'/',i4,' ',i2,':',i2)

! End of program

      END 

