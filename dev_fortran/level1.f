      PROGRAM COMPLY

c Declare dynamic arrays

      REAL*4, ALLOCATABLE :: latitude(:)        !latitude (decimal degrees)
      REAL*4, ALLOCATABLE :: longitude(:)       !longitude (decimal degrees)
      REAL*4, ALLOCATABLE :: rademis(:,:)       !radionuclide levels for each stack
      REAL*4, ALLOCATABLE :: fract(:,:)         !radionuclide possession fraction
      REAL*4, ALLOCATABLE :: posmax(:)          !maximum possession limit for solid, liquid, and gas forms

      INTEGER, ALLOCATABLE :: nrads(:)          !number of emitted radionucludes at each stack

      CHARACTER*8, ALLOCATABLE :: radname(:,:)  !radionuclide names at each stack
      CHARACTER*1, ALLOCATABLE :: radtype(:,:)  !radionuclide type: gas, liquid, solid


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

       allocate(latitude(nstacks))
       allocate(longitude(nstacks))
       allocate(nrads(nstacks))

       do n = 1, nstacks 

       read(11,'(20x,f10.5)') latitude(n)
       read(11,'(20x,f10.5)') longitude(n)

       read(11,'(20x,i3)') nrads(n)

       allocate(radname(nstacks,nrads(n)))
       allocate(rademis(nstacks,nrads(n)))
       allocate(fract(nstacks,nrads(n)))
       allocate(radtype(nstacks,nrads(n)))
       allocate(posmax(4))

       do m = 1, nrads(n)
       read(11,'(20x,a)') inline
       read(inline,*) radname(n,m),radtype(n,m),rademis(n,m)
       enddo !end loop over reading in radionuclide emissions
   
       enddo !end loop over reading in stacks

c Read in ancillary files

       write(*,*)
       read(13,'(a)') inline
       
        do i = 1, 823
        read(13,'(a)' ) inline
        read(inline,*) NUC, (POSMAX(j),j=1,4)
         ! 1=gas, 2=liquid, and 3=solid
          do n = 1, nstacks
           do m = 1, nrads(n)
           if(radname(n,m).eq.NUC) then
            write(*,*)'Match',NUC,radtype(n,m)
            if(radtype(n,m).eq.'G'.or.radtype(n,m).eq.'g') then !gas
            fract(n,m) = rademis(n,m)/POSMAX(1)
            write(*,*) 'Possession,Limit,Fraction:',NUC,rademis(n,m),POSMAX(1),fract(n,m)
            else if(radtype(n,m).eq.'L'.or.radtype(n,m).eq.'l') then !liquid
            fract(n,m) = rademis(n,m)/POSMAX(2)
            write(*,*)'Possession,Limit,Fraction:',NUC,rademis(n,m),POSMAX(2),fract(n,m)
            else if(radtype(n,m).eq.'S'.or.radtype(n,m).eq.'s') then !solid
            fract(n,m) = rademis(n,m)/POSMAX(3)
            write(*,*)'Possession,Limit,Fraction:',NUC,rademis(n,m),POSMAX(3),fract(n,m)
            else if(radtype(n,m).eq.'C'.or.radtype(n,m).eq.'c') then !concentration comparison
            fract(n,m) = rademis(n,m)/POSMAX(4)
            write(*,*)'Possession,Limit,Fraction:',NUC,rademis(n,m),POSMAX(4),fract(n,m)
            else
            write(*,*)'Select C for concentration table or G, L, or S for possession table.'
            write(*,*)'Unrecognized physical state (should be G, L, or S):',radtype(n,m)
            fract(n,m) = 9999.
            endif
           endif
           enddo !end loop over radionuclides
          enddo !end loop over stacks

          enddo !end loop over lines on the input file 
 998    continue

        write(*,*)

c Calculate possession

       alldose = 0.
       idose = 0.

       do n = 1, nstacks
       do m = 1, nrads(n)
       alldose = alldose + fract(n,m)
         read(radname(n,m),'(a2)') NN
         if(NN.eq.'I-') then
         idose = idose + fract(n,m)
         else
         endif
       enddo !end loop over radionuclides
       enddo !end loop over stacks

       write(*,*)'Total possession (fraction):',alldose
       write(*,*)'Total iodine possession (fraction):',idose

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
       write(12,1008) n
       write(12,1006) latitude(n)
       write(12,1007) longitude(n)

       write(12,1013) nrads(n)

       do m=1, nrads(n)
       write(12,1014) radname(n,m),radtype(n,m),rademis(n,m),fract(n,m)
       enddo !end loop over radionuclides

       enddo !end loop over stacks

       write(12,2002)
       write(12,2001)
       write(12,2001)

       write(12,2000) alldose
       write(12,1999) idose
       write(12,2001)
       write(12,2001)
       write(12,4000)


c Add time stamp of model simulation to output file

       write(*,*)
       call date_and_time(b(1), b(2), b(3), date_time)
       print *,'Program execution date and time information:'
       print *,'year=',date_time(1)
       print *,'month_of_year=',date_time(2)
       print *,'day_of_month=',date_time(3)
       print *,'hour of day=',date_time(5)
       print *,'minutes of hour=',date_time(6)
       print *,'seconds of minute=',date_time(7)

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
 1008 FORMAT('Stack number:',i3)
 1013 FORMAT('Number of radionuclides:',I2)
 1014 FORMAT(a,',',a,',',f25.20,',',f25.20)
 1999 FORMAT('Total possession iodine fraction:',f10.2)
 2000 FORMAT('Total possession fraction:',f10.2)
 2001 FORMAT('-----------------------------------------------------------------------')
 2002 FORMAT(' ')
 2100 FORMAT('Date/Time (LST) stamp of simulation (MM/DD/YYYY HH:MM):',i2,'/',i2,'/',i4,' ',i2,':',i2)
 4000 FORMAT('COMPLY version 2.0 level 1 distributed through GitHub')

c End of program

      END 

