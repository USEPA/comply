      PROGRAM COMPLY

      IMPLICIT NONE
c Declare dynamic arrays

      REAL*4, ALLOCATABLE :: releaseht(:)           !stack height (m)
      REAL*4, ALLOCATABLE :: bldheight(:)           !building height (m)
      REAL*4, ALLOCATABLE :: bldwidth(:)            !building width (m)
      REAL*4, ALLOCATABLE :: distance(:)            !distance between stack and receptor (m)
      REAL*4, ALLOCATABLE :: distanceveg(:)         !distance between garden and receptor (m)
      REAL*4, ALLOCATABLE :: distancecow(:)         !distance between cows/milk and receptor (m)
      REAL*4, ALLOCATABLE :: distancemeat(:)        !distance between meat and receptor (m)
      REAL*4, ALLOCATABLE :: windspeed(:)           !wind speed (m/s)
      REAL*4, ALLOCATABLE :: latitude(:)            !latitude (decimal degrees)
      REAL*4, ALLOCATABLE :: longitude(:)           !longitude (decimal degrees)
      REAL*4, ALLOCATABLE :: pctcalm(:)             !percentage of calm days at wind data location

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
      write(*,*)'COMPLY v2.0 for internal testing only'
      write(*,*)
 
      open(11,file='input.level4.txt',status='old',form='formatted')
      open(12,file='output.level4.txt',status='unknown',form='formatted')
      open(13,file='diag.level4.txt',status='unknown',form='formatted')

      open(9,file='comply171_original_input_data.csv',status='old',form='formatted')
      open(10,file='output.level4.tree.csv',status='old',form='formatted')

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
       allocate(distanceveg(nstacks))
       allocate(distancecow(nstacks))
       allocate(distancemeat(nstacks))
       allocate(vegsector(nstacks))
       allocate(milksector(nstacks))
       allocate(meatsector(nstacks))
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
       read(11,'(20x,f10.5)') bldheight(n)
       read(11,'(20x,f10.5)') bldwidth(n)
       read(11,'(20x,f10.5)') distance(n)
       read(11,'(20x,f10.5)') distanceveg(n)
       read(11,'(20x,f10.5)') distancecow(n)
       read(11,'(20x,f10.5)') distancemeat(n)

       read(11,'(20x,a)') vegsector(n)
       read(11,'(20x,a)') milksector(n)
       read(11,'(20x,a)') meatsector(n)

       read(11,'(20x,f10.5)') windspeed(n)

       read(11,'(20x,i3)') nrads(n)

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



!Calculate air concentrations
!Emission rate assumption is Ci/s

       write(13,3001) !write header to file

       do n = 1, nstacks

       alldose = 0.

       if(releaseht(n).eq.0.0) then !surface release
       bldwidth(n) = 0.0 !odd that comply assumes no release shape geometry for a surface (area) release
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

!move reading in decay chains to here??



       do j = 1, nwindbins(n)

       do receptor = 10, 3010, 5

       dist = float(receptor)

       SIGMAZ = (0.06 * dist) / SQRT( 1.0 + 0.0015 * dist)
       CAPSIG = SQRT( SIGMAZ ** 2 + area / 3.14159265 )
       pp = 2.032 / ( dist * CAPSIG )
       transittime = (dist / windbinspeed(n,j))*(1./3600.)*(1./24.) ! seconds to days
       fdux = (windbinfreq(n,j)/(1.0-pctcalm(n))) / (windbinspeed(n,j)*dist) !s/m2 units

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

       !rainfall rate supposedly m/yr
       !emissions should be in Ci/s
       !concentration should be in Ci/m3
       !dry dep velocity in m/day
       wetdepadj = (8.0 * fdux * 1.0 * 1.0E-05 * (rademis(n,m) * 3600. * 24.))/(3.1416*conc(n,m) )
       wetdepadjbad = (8.0 * fdux * 1.0 * 1.0E-05 * (rademis(n,m) * 3600. * 24.))/(3.1416*conc(n,m)/rademis(n,m) )
       !!! error in comply171       wetdepadj = (8.0 * fdux * 1.0 * 1.0E-05 * (rademis(n,m) * 3600. * 24.))/(3.1416*conc(n,m)/rademis(n,m) )

!       if(dist.eq.240) then
!               write(*,*)'AIR CONC:',radname(n,m),conc(n,m),windbinfreq(n,j),windbinspeed(n,j),wetdepadj
!       endif

!Find decay chain and relevant half-life values

       rewind(9)
       rewind(10)
       read(9,'(a)') inline

         do ii = 1, 879
         read(10,'(a)') inline
         read(inline,*) NUC,k,(decaychain(l),l=1,k)
         if(NUC.eq.radname(n,m)) then
          read(inline,*) NUC,k,(decaychain(l),l=1,k)
!          write(*,*)'Decay chain:', NUC,k,(decaychain(l),l=1,k)
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
!            inhalationfac(l) = 0.0
!            immersionfac(l) = 0.0
!            surfacefac(l) = 0.0
!            ingestionfac(l) = 1.0
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

!         if(dist.eq.10) then
!          do l = 1, k
!          write(*,*)'Decay chain:', k,decaychainname(l),halflife(l),wbfac(l),inhalationfac(l),immersionfac(l),surfacefac(l),
!     &      ingestionfac(l),drydepvel(l)
!            write(*,*) 'Decay chain:',decaychainname(l),xbveg(l),xbfor(l),xfmilk(l),xfmeat(l)
!          enddo
!          endif

!Calculate parent and progeny ambient activity

!         write(*,*)'Calculating activity for nuclides in the decay chain (N=',nchain,')'

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

!Calculate parent and progeny activity at the ground level

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
!         if(dist.eq.240) then
!         write(*,*)'!',l,decaychainname(l),activity(l),drydepvel(l),wetdepadj,surfaceactivity(l),surfacefac(l),surfdose(l)
!         endif
         totsurfdose(n,m) = totsurfdose(n,m) + surfdose(l)
         enddo

!Calculate parent and progeny activity for ingestion pathways

!!! Normalized surface activity 

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

!!! Vegetables

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

!!! milk

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

!!! meat

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

!Estimate concentration in media

         read(radname(n,m),'(a2)') NN

         do l = 1, nchain

          if(NN.eq.'H-') then !H-3 (tritium) -  D. A. Baker, et al, EPA 600/9-76-016, July 1976.
          cfveg(l) = conc(n,m)*9.0*(0.1/0.8E-02) ! (FHVEG / CH2O )
          !CFVEG(  KDTR ) = ( VEGAIR * 9.0 * FHVEG  / CH2O )
          cfmilk(l) = (conc(n,m)*9.0*7.0E-02/8.0E-03)*16.0*0.11/(7.0E-02*16. + 220.0/9.0)
          !CFMILK( KDTR ) = ( MILKAIR * 9.0 * FHFEED / CH2O ) * QMILK * FHMILK / ( FHFEED * QMILK) + (QWMILK / 9.0 )
          cfmeat(l) = (conc(n,m)*9.0*7.0E-02/8.0E-03)*12.0*0.10/(7.0E-02*12. + 50.0/9.0) 
          !CFMEAT( KDTR ) = ( MEATAIR * 9.0 * FHFEED / CH2O ) * QMEAT * FHMEAT / ( FHFEED * QMEAT + QWMEAT / 9.0 )

          elseif(NN.eq.'C-') then !C14 - D. A. Baker, et al, EPA 600/9-76-016, July 1976.
          cfveg(l) = conc(n,m)*(9.000000E-02/0.16E-03)
          !CFVEG(  KDTR ) = VEGAIR  * FCVEG  / CCARB
          cfmilk(l) = conc(n,m)*(7.0000000E-02/0.16E-03)
          !CFMILK( KDTR ) = MILKAIR * FCMILK / CCARB
          cfmeat(l) = conc(n,m)*(0.2400000/0.16E-03)
          !CFMEAT( KDTR ) = MEATAIR * FCMEAT / CCARB

          else !Equation B2 of NCRP Commentary 3
          cfveg(l) =  conc(n,m) * ((0.1/1.0)*aveglf(l)+(xbveg(l)/220.)*avegrt(l))
          !CFVEG( KDTR )  = VEGAIR * (( FRVEG / YVEG ) * AVEGLF( KDTR ) +( BVEG  / PVEG ) * AVEGRT( KDTR ) )
          cfmilk(l) = conc(n,m) * ((1.8/1.0)*amilklf(l)+(xbfor(l)/220.)*amilkrt(l))*xfmilk(l)*16.0
          !CFMILK(KDTR)= MILKAIR*((FRMILK/YMILK)*AMILKLF(KDTR)+(BFOR/PMILK)*AMILKRT(KDTR))*FMILK*QMILK
          cfmeat(l) = conc(n,m) * ((1.8/1.0)*ameatlf(l)+(xbfor(l)/220.)*ameatrt(l))*xfmeat(l)*12.0
          !CFMEAT(KDTR)=MEATAIR*((FRMEAT/YMEAT)*AMEATLF(KDTR)+(BFOR/PMEAT)*AMEATRT(KDTR))*FMEAT*QMEAT
          endif

          enddo


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


          vegdose(n,m) = 0.0
          milkdose(n,m) = 0.0
          meatdose(n,m) = 0.0
          do l = 1, nchain
          vegdose(n,m) = vegdose(n,m) + dveg(l)
          milkdose(n,m) = milkdose(n,m) + dmilk(l)
          meatdose(n,m) = meatdose(n,m) + dmeat(l)
          enddo

! Save dose information at key input receptor locations for report

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
!         read(radname(n,m),'(a2)') NN
!         if(NN.eq.'I-') then
!         idose = idose + totimdose(n,m) + totihdose(n,m) + totsurfdose(n,m)
!         endif
        endif
       endif

       write(13,3000) n,radname(n,m),dist,
     &     conc(n,m),totimdose(n,m),totihdose(n,m),totsurfdose(n,m),
     &     meatdose(n,m),milkdose(n,m),vegdose(n,m)

       enddo !end loop over distances between source and receptor (receptor)
       enddo !end loop over number of wind direction bins (j)

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

       write(*,*)
       write(*,4000)'Total dose (mrem/yr):',alldose
       write(*,*)'Total dose limit is 10.0 mrem/yr'
       write(*,*)'Iodine limit is 3.0 mrem/yr'
       write(*,*)

       enddo !end loop over stacks (n)

 2999  FORMAT(a30,a10,f20.15)
 3000  FORMAT(i2,',',a,',',f6.1,',',f30.20,6(',',f20.10))
 3001  FORMAT('STACK,RADNAME,DISTANCE,CONC,IMMERSION,INHALATION,SURFACE,MEAT,MILK,VEG')
 4000  FORMAT(a,f10.4)

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
       write(12,1014) radname(n,m),rademis(n,m),conc(n,m),(dose(n,m)+vegdose(n,m)+meatdose(n,m)),dose(n,m),vegdose(n,m)
       enddo !end loop over radionuclides

       enddo !end loop over stacks

       write(12,2000) alldose
       write(12,2002) totvegdose/alldose
       write(12,2004) totmilkdose/alldose


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
 1014 FORMAT(a,',',f25.20,',',f25.20,4(',',f10.6))
 2000 FORMAT('Total effective dose equivalent (mrem/yr):',f10.2)
 2001 FORMAT('Total effective dose equivalent iodine only (mrem/yr):',f10.2)

 2002 FORMAT('Vegetable fraction of total effective dose equivalent:',f5.2)
 2003 FORMAT('Vegetable fraction of total effective dose equivalent iodine only:',f5.2)
 2004 FORMAT('Meat/milk fraction of total effective dose equivalent:',f5.2)
 2005 FORMAT('Meat/milk fraction of total effective dose equivalent iodine only:',f5.2)


 2100 FORMAT('Date/Time (LST) stamp of simulation (MM/DD/YYYY HH:MM):',i2,'/',i2,'/',i4,' ',i2,':',i2)

! End of program

      END 

