       PROGRAM COMPLY

       IMPLICIT NONE
c Declare dynamic arrays

       REAL*4, ALLOCATABLE :: releaseht(:)           !stack height (m)
       REAL*4, ALLOCATABLE :: diameter(:)            !stack diameter (m)
       REAL*4, ALLOCATABLE :: exitvelocity(:)        !exit velocity ()
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

       REAL*8, ALLOCATABLE :: rademis(:,:)           !radionuclide emissions (Ci/s) for each stack
       REAL*8, ALLOCATABLE :: conc(:,:)              !radionuclide air concentration (Bq/m3)
       REAL*8, ALLOCATABLE :: dose(:,:)              !radionuclide dose (mrem/yr)

       REAL*8, ALLOCATABLE :: sourceterm(:)
       REAL*8, ALLOCATABLE :: initialamt(:)
       REAL*8, ALLOCATABLE :: activity(:)
       REAL*8, ALLOCATABLE :: surfaceactivity(:)
       REAL*8, ALLOCATABLE :: ambientactivity(:)
       REAL*8, ALLOCATABLE :: normsourceactivity(:)

       REAL*8, ALLOCATABLE :: aveglf(:)
       REAL*8, ALLOCATABLE :: avegrt(:)
       REAL*8, ALLOCATABLE :: amilklf(:)
       REAL*8, ALLOCATABLE :: amilkrt(:)
       REAL*8, ALLOCATABLE :: ameatlf(:)
       REAL*8, ALLOCATABLE :: ameatrt(:)
       REAL*8, ALLOCATABLE :: cfveg(:)
       REAL*8, ALLOCATABLE :: cfmilk(:)
       REAL*8, ALLOCATABLE :: cfmeat(:)
       REAL*8, ALLOCATABLE :: dveg(:)
       REAL*8, ALLOCATABLE :: dmilk(:)
       REAL*8, ALLOCATABLE :: dmeat(:)
       REAL*8, ALLOCATABLE :: ingdose(:)

       REAL*8, ALLOCATABLE :: kvalue(:)
       REAL*8, ALLOCATABLE :: ktot(:)
       REAL*8, ALLOCATABLE :: xbveg(:)
       REAL*8, ALLOCATABLE :: xbfor(:)
       REAL*8, ALLOCATABLE :: xfmilk(:)
       REAL*8, ALLOCATABLE :: xfmeat(:)

       REAL*4, ALLOCATABLE :: halflife(:)          !radionuclide half life (years)
       REAL*4, ALLOCATABLE :: wbfac(:)             !whole body dose factor from NCRP Commentary 3 table 2
       REAL*4, ALLOCATABLE :: xwbfac(:)            !whole body dose factor from NCRP Commentary 3 table 3
       REAL*4, ALLOCATABLE :: vegfac(:)            !vegetable dose factor from NCRP Commentary 3 table 7
       REAL*4, ALLOCATABLE :: cowfac(:)            !cow dose factor from NCRP Commentary 3
       REAL*8, ALLOCATABLE :: inhalationfac(:)     !inhalation dose factor from ICRP table 2.1  
       REAL*8, ALLOCATABLE :: immersionfac(:)      !immersion dose factor from NCRP Commentary Table C5
       REAL*8, ALLOCATABLE :: surfacefac(:)        !surface dose factor from NCRP Commentary Table C6
       REAL*8, ALLOCATABLE :: ingestionfac(:)      !ingestion dose factor from ICRP Table 2.2
       REAL*8, ALLOCATABLE :: drydepvel(:)         !dry deposition velocity (meters/day) from NCRP Commentary Table L4

       REAL*8, ALLOCATABLE :: imdose(:)
       REAL*8, ALLOCATABLE :: ihdose(:)
       REAL*8, ALLOCATABLE :: surfdose(:)
       REAL*8, ALLOCATABLE :: totimdose(:,:)
       REAL*8, ALLOCATABLE :: totihdose(:,:)
       REAL*8, ALLOCATABLE :: totsurfdose(:,:)
       REAL*8, ALLOCATABLE :: totingdose(:,:)

       REAL*8, ALLOCATABLE :: milkdose(:,:)          !radionuclide dose by milk (mrem/yr)
       REAL*8, ALLOCATABLE :: vegdose(:,:)           !radionuclide dose by vegetables (mrem/yr)
       REAL*8, ALLOCATABLE :: meatdose(:,:)          !radionuclide dose by meat (mrem/yr)
       REAL*8, ALLOCATABLE :: windbinfreq(:,:)       !Sector (wind dir) frequency from that direction
       REAL*8, ALLOCATABLE :: windbinspeed(:,:)      !Sector (wind dir) average speed from that direction (m/s)
       REAL*8, ALLOCATABLE :: windbinsrdist(:,:)     !Sector (wind dir) source-receptor distance (m)

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

c Define constants
       INTEGER, PARAMETER :: MAX_CHAIN = 19
       INTEGER, PARAMETER :: MAX_NUCLIDES = 826
       INTEGER, PARAMETER :: MAX_DECAYCHAIN_INPUT = 879
       REAL, PARAMETER :: BREATH_RATE = 8000.0      !m3/yr breathing rate
       REAL, PARAMETER :: CONVEG = 70.0               !vegetable consumption
       REAL, PARAMETER :: CONMILK = 160.0             !milk consumption
       REAL, PARAMETER :: CONMEAT = 75.0              !meat consumption
       REAL, PARAMETER :: REMCON = 1.5000000E-02      !removal constant (1/yr)
       REAL, PARAMETER :: DEPTIM = 100.0 * 365.0      !deposition time (days)
       REAL, PARAMETER :: TWEATH = 0.693 / 12.0        !weathering half-life
       REAL, PARAMETER :: TEVEG = 60.0                 !vegetable transit time (days)
       REAL, PARAMETER :: TEMILK = 30.0                !milk transit time (days)
       REAL, PARAMETER :: TEMEAT = 30.0                !meat transit time (days)
       REAL, PARAMETER :: DELVEG = 11.0                !vegetable delivery time
       REAL, PARAMETER :: DELMILK = 2.0                !milk delivery time
       REAL, PARAMETER :: DELMEAT = 17.0               !meat delivery time
       REAL, PARAMETER :: PI = 3.14159265

c Declare other variables and calculated terms

       REAL FACTOR,FACTORB,FACTORC,FACTORD,FACTORE,FACTORF,FACTORG,FACTORH,FACTORI
       REAL FACTORJ,FACTORK,FACTORL,FACTORM,FACTORN
       REAL area,xcrit,SIGMAZ,CAPSIG,pp,alldose,idose,dist,altconc
       REAL tempdose,maxdose,maxdist,totmilkdose,totvegdose,itotvegdose,totmeatdose,itotcowdose
       REAL capp,capn0,outprod,inprod,z1,term1,term2,sum
       REAL transittime                              !transit time of air from emission release to receptor (s)
       REAL fdux                                     !not sure what this is-units are 1/(m2/s) or s/m2
       REAL wetdepadj,wetdepadjbad                   !wet deposition adjustment factor
       REAL xtime
       REAL xremove
       REAL ingestion, immersion, inhalation, surface
       REAL RHOAIR,RHOSTK,STKVEL,R,FM,DELHMO,DELHBU,HEFF,FB,zz,INHFAC,TINIT,T

       INTEGER nstacks                               !number of unique stacks
       INTEGER receptor, N, M, I, J, K, L, II, JJ, P, Q                        !loop counters
       INTEGER nchain                                !number of nuclides in the decay series
       INTEGER maxsect

       CHARACTER*200 title,company,facility,address,poc,phone,inline
       CHARACTER*8 NUC
       CHARACTER*2 NN

       integer date_time(8)                      !for date stamp function
       character*10 b(3)                         !for date stamp function

c Open input and output files 

       write(*,*)
       write(*,*)'COMPLY v2.0 (Optimized)'
       write(*,*)
  
       open(11,file='input.level4.txt',status='old',form='formatted')
       open(12,file='output.level4.txt',status='unknown',form='formatted')
       open(13,file='diag.level4.txt',status='unknown',form='formatted')

       open(9,file='comply171_original_input_data.csv',status='old',form='formatted')
       open(10,file='comply171_original_decaychain_data.csv',status='old',form='formatted')

c Read input namelist file

        read(11,'(20x,a)') title
        read(11,'(20x,a)') company
        read(11,'(20x,a)') facility
        read(11,'(20x,a)') address
        read(11,'(20x,a)') poc
        read(11,'(20x,a)') phone

        read(11,'(20x,i3)') nstacks

c Allocate arrays that depend only on nstacks
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

c Allocate fixed-size arrays once BEFORE the stack loop (BUG FIX)
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

c Allocate stack and radionuclide-dependent arrays (moved here - BUG FIX)
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

        allocate(windbinname(nstacks,nwindbins(n)))
        allocate(windbinfreq(nstacks,nwindbins(n)))
        allocate(windbinspeed(nstacks,nwindbins(n)))
        allocate(windbinsrdist(nstacks,nwindbins(n)))
        allocate(sectorname(nstacks,nwindbins(n)))

        do m = 1, nrads(n)
        read(11,'(20x,a)') inline
        read(inline,*) radname(n,m),rademis(n,m)
        enddo !end loop over reading in radionuclide emissions

        read(11,'(20x,a)') windsource(n)
        read(11,'(20x,a)') winddates(n)
        read(11,'(20x,a)') windloc(n)
        read(11,'(20x,a)') winddist(n)

        read(11,'(20x,f10.5)') pctcalm(n)
        read(11,'(20x,i3)') nwindbins(n)

        do m = 1, nwindbins(n)
        read(11,'(20x,a)') inline
        read(inline,*) windbinname(n,m),windbinfreq(n,m),windbinspeed(n,m),windbinsrdist(n,m),sectorname(n,m)
        enddo !end loop over reading in wind data and source-receptor distances

        enddo !end loop over reading in stacks

c Calculate air concentrations
c Emission rate assumption is Ci/s

        write(13,3001) !write header to file

        do n = 1, nstacks

        alldose = 0.
        maxsect = 0
        maxdist = 0.0

        if(releaseht(n).eq.0.0) then !surface release
        bldwidth(n) = 0.0
        bldheight(n) = 0.0
        endif

        area = bldwidth(n) * bldheight(n)
        xcrit = 2.5 * sqrt(area)

        do m = 1, nrads(n)

        maxdose = 0.0
        inhalation = 0.
        immersion = 0.
        ingestion = 0.
        totmilkdose = 0.
        totvegdose = 0.
        totmeatdose = 0.
        surface = 0.

        do j = 1, nwindbins(n)

        do receptor = 10, 3010, 10

        dist = float(receptor)

        if(exitvelocity(n).le.0.) then
        exitvelocity(n) = 0.3
        endif

        HEFF = releaseht(n)

        IF ( releaseht(n) .GT. 2.5*bldheight(n) ) THEN

c Temperature
         RHOAIR = 0.0735 * ( 80.0 + 460.0 ) / ( airtemp(n) + 460.0 )
         RHOAIR = RHOAIR * 454.0 / ( 0.3048 ** 3 )
         RHOSTK = RHOAIR * ( airtemp(n) + 460.0 ) / ( stacktemp(n) + 460.0 )

c Momentum
         AREA   = PI * ( diameter(n) ** 2.0 ) / 4.0
         STKVEL = exitvelocity(n) / AREA
         R      = STKVEL / windbinspeed(n,j)
         FM     = ( RHOSTK / RHOAIR ) * STKVEL * exitvelocity(n) / PI
         TERM1  = 75.0 * FM / ( ( ( 1.0 + 3.0 / R ) * windbinspeed(n,j) ) ** 2.0 )
         TERM2  = releaseht(n) + 7.79 * SQRT( FM ) / ( ( 1.0 + 3.0 / R ) * windbinspeed(n,j) )

c If the site has a vertical stack, then momentum rise value is zero, otherwise compute the value.
         IF ( verticalmomentum(n).eq.'YES' ) THEN
          DELHMO = 0.93 * ( TERM1 ** ( 3.0 / 7.0 ) ) *
     &        ( TERM2 ** ( 1.0 / 7.0 ) )
         ELSE
          DELHMO = 0.0
         ENDIF

c Buoyant rise
         FB = 9.81 * ( RHOAIR - RHOSTK ) * exitvelocity(n) / ( PI * RHOAIR )
         DELHBU = 23.67 * ( ( FB / ( windbinspeed(n,j) ** 3.0 ) ) ** 0.6 ) *
     &           ( releaseht(n) + 195.1 * FB / ( windbinspeed(n,j) ** 3.0 ) ) ** 0.4

c Use maximum of buoyant or momentum rise
        HEFF = releaseht(n) + AMAX1( DELHMO, DELHBU )

        ENDIF !end condition where H > 0


        SIGMAZ = (0.06 * dist) / SQRT( 1.0 + 0.0015 * dist)
        CAPSIG = SQRT( SIGMAZ ** 2 + area / PI )
        zz = 0.5 * ( ( HEFF / CAPSIG ) ** 2 )
         if(zz.gt.0.001.and.releaseht(n).gt.2.5*bldheight(n)) then
         pp = 2.032 * EXP( -zz ) / ( dist * CAPSIG )
         else
         pp = 2.032 / ( dist * CAPSIG )
         endif


        transittime = (dist / windbinspeed(n,j))*(1./3600.)*(1./24.) ! days
        fdux = (windbinfreq(n,j)/(1.0-pctcalm(n))) / (windbinspeed(n,j)*dist) !s/m2 units

        if(dist.lt.xcrit) then
        conc(n,m) = ( (windbinfreq(n,j)/(1.0-pctcalm(n))) *rademis(n,m)) /
     &          (PI*windbinspeed(n,j)*1.0*amin1(bldwidth(n),bldheight(n)))
        else
        conc(n,m) = ((windbinfreq(n,j)/(1.0-pctcalm(n))) *rademis(n,m)*pp) / windbinspeed(n,j)
        endif

        if(conc(n,m).lt. 1.0E-32 ) then
          conc(n,m) = 1.0E-32
        endif

        altconc = (windbinfreq(n,j)*rademis(n,m)*pp) / windbinspeed(n,j) !QA only

c Wet deposition adjustment
        wetdepadj = (8.0 * fdux * rainfall(n) * 1.0E-05 * (rademis(n,m) * 3600. * 24.))/(PI*conc(n,m))

c Find decay chain and relevant half-life values (moved from original nested loop)
        rewind(9)
        rewind(10)
        read(9,'(a)') inline

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
             halflife(l) = FACTOR*365.
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

c Calculate immersion and inhalation dose
          call init_array(activity, nchain)
          call init_array(sourceterm, nchain)
          call init_array(initialamt, nchain)
          
          initialamt(1) = conc(n,m)
          xremove = 0.0
          xtime = transittime

       CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

          totimdose(n,m) = 0.0
          totihdose(n,m) = 0.0

          do l = 1, nchain
          ambientactivity(l) = activity(l)
          imdose(l) = ambientactivity(l) * immersionfac(l)
          ihdose(l) = ambientactivity(l) * inhalationfac(l) * BREATH_RATE
          totimdose(n,m) = totimdose(n,m) + imdose(l)
          totihdose(n,m) = totihdose(n,m) + ihdose(l)
          enddo

c Calculate surface dose

          call init_array(sourceterm, nchain)
          call init_array(initialamt, nchain)
          call init_array(activity, nchain)

          do l = 1, nchain
          sourceterm(l) = ambientactivity(l)
          enddo

          xremove = REMCON / 365.
          xtime = DEPTIM

       CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

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

c Calculate ingestion pathways - normalized surface activity

          call init_array(sourceterm, nchain)
          call init_array(initialamt, nchain)
          call init_array(activity, nchain)

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

c Vegetables pathway
          call init_array(sourceterm, nchain)
          call init_array(initialamt, nchain)
          call init_array(activity, nchain)

          do l = 1, nchain
          sourceterm(l) = normsourceactivity(l)
          enddo

          xremove = TWEATH
          xtime = TEVEG

       CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

          do l = 1, nchain
          aveglf(l) = activity(l) 
          enddo

          call init_array(sourceterm, nchain)
          call init_array(initialamt, nchain)
          call init_array(activity, nchain)

          do l = 1, nchain
          sourceterm(l) = normsourceactivity(l)
          enddo

          xremove = REMCON / 365.
          xtime = DEPTIM

       CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

          do l = 1, nchain
          avegrt(l) = activity(l)
          enddo

c Milk pathway
          call init_array(sourceterm, nchain)
          call init_array(initialamt, nchain)
          call init_array(activity, nchain)

          do l = 1, nchain
          sourceterm(l) = normsourceactivity(l)
          enddo

          xremove = TWEATH
          xtime = TEMILK

       CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

          do l = 1, nchain
          amilklf(l) = activity(l)
          enddo

          call init_array(sourceterm, nchain)
          call init_array(initialamt, nchain)
          call init_array(activity, nchain)

          do l = 1, nchain
          sourceterm(l) = normsourceactivity(l)
          enddo

          xremove = REMCON / 365.
          xtime = DEPTIM

       CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

          do l = 1, nchain
          amilkrt(l) = activity(l)
          enddo

c Meat pathway
          call init_array(sourceterm, nchain)
          call init_array(initialamt, nchain)
          call init_array(activity, nchain)

          do l = 1, nchain
          sourceterm(l) = normsourceactivity(l)
          enddo

          xremove = TWEATH
          xtime = TEMEAT

       CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

          do l = 1, nchain
          ameatlf(l) = activity(l)
          enddo

          call init_array(sourceterm, nchain)
          call init_array(initialamt, nchain)
          call init_array(activity, nchain)

          do l = 1, nchain
          sourceterm(l) = normsourceactivity(l)
          enddo
  
          xremove = REMCON / 365.
          xtime = DEPTIM

       CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

          do l = 1, nchain
          ameatrt(l) = activity(l)
          enddo

c Estimate concentration in media

          read(radname(n,m),'(a2)') NN

          do l = 1, nchain

           if(NN.eq.'H-') then !H-3 (tritium) - D. A. Baker, et al, EPA 600/9-76-016, July 1976.
           cfveg(l) = conc(n,m)*9.0*(0.1/0.8E-02)
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

c Vegetables ingestion dose
            call init_array(sourceterm, nchain)
            call init_array(initialamt, nchain)
            call init_array(activity, nchain)

            do l = 1, nchain
            sourceterm(l) = 0.0
            initialamt(l) = cfveg(l)
            enddo

            xremove = 0.0
            xtime = DELVEG

            CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

             do l = 1, nchain
             dveg(l) = activity(l) * ingestionfac(l) * CONVEG
             enddo

c Milk ingestion dose
            call init_array(sourceterm, nchain)
            call init_array(initialamt, nchain)
            call init_array(activity, nchain)

            do l = 1, nchain
            sourceterm(l) = 0.0
            initialamt(l) = cfmilk(l)
            enddo

            xremove = 0.0
            xtime = DELMILK

            CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

             do l = 1, nchain
             dmilk(l) = activity(l) * ingestionfac(l) * CONMILK
             enddo

c Meat ingestion dose
            call init_array(sourceterm, nchain)
            call init_array(initialamt, nchain)
            call init_array(activity, nchain)

            do l = 1, nchain
            sourceterm(l) = 0.0
            initialamt(l) = cfmeat(l)
            enddo

            xremove = 0.0
            xtime = DELMEAT

            CALL skrable(activity,sourceterm,xtime,nchain,xremove,halflife,initialamt)

             do l = 1, nchain
             dmeat(l) = activity(l) * ingestionfac(l) * CONMEAT
             enddo

c Total ingestion doses
           vegdose(n,m) = 0.0
           milkdose(n,m) = 0.0
           meatdose(n,m) = 0.0
           do l = 1, nchain
           vegdose(n,m) = vegdose(n,m) + dveg(l)
           milkdose(n,m) = milkdose(n,m) + dmilk(l)
           meatdose(n,m) = meatdose(n,m) + dmeat(l)
           enddo

c Save dose information at key input receptor locations for report

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

        write(13,3000) n,radname(n,m),sectorname(n,j),dist,
     &     conc(n,m),totimdose(n,m),totihdose(n,m),totsurfdose(n,m),
     &     meatdose(n,m),milkdose(n,m),vegdose(n,m)

        enddo !end loop over distances between source and receptor
        enddo !end loop over wind direction bins

        alldose = alldose + surface + immersion + inhalation
        if(maxsect .gt. 0) then
          write(*,*)'Max dose (',maxdose,') receptor at ',windbinsrdist(n,maxsect),sectorname(n,maxsect)
        endif
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
        enddo !end loop over radionuclides

        write(*,*)
        write(*,3100)'Total dose (mrem/yr):',alldose
        write(*,*)'Total dose limit is 10.0 mrem/yr'
        write(*,*)'Iodine limit is 3.0 mrem/yr'
        write(*,*)

        enddo !end loop over stacks

  2999  FORMAT(a30,a10,f20.15)
  3000  FORMAT(i2,',',a,',',a,',',f6.1,',',f30.20,6(',',f20.10))
  3001  FORMAT('STACK,RADNAME,SECTOR,DISTANCE,CONC,IMMERSION,INHALATION,SURFACE,MEAT,MILK,VEG')
  3100  FORMAT(a,f10.4)

c Write standardized output

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
        if(maxsect .gt. 0) then
          write(12,1800)'Max dose receptor ',windbinsrdist(n,maxsect),' m downwind in sector ',sectorname(n,maxsect)
        endif
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

        write(12,2002)
        write(12,2001)
        write(12,2002)
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

c End of program

       END 

c Helper subroutine to initialize arrays efficiently
       SUBROUTINE init_array(arr, n)
       IMPLICIT NONE
       INTEGER n, i
       REAL*8 arr(n)
       do i = 1, n
         arr(i) = 0.0
       enddo
       END
