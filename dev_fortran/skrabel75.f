
      SUBROUTINE skrable( activity, sourceterm, xtime, nchain,
     &                   xremove, halflife, initialamt )

      IMPLICIT NONE

!      REAL*4, ALLOCATABLE :: sourceterm(:)
!      REAL*4, ALLOCATABLE :: initialamt(:)
!      REAL*4, ALLOCATABLE :: activity(:)
!      REAL*4, ALLOCATABLE :: kvalue(:)
!      REAL*4, ALLOCATABLE :: ktot(:)

!      REAL, ALLOCATABLE :: halflife(:)          !radionuclide half life (years)
      REAL*8 activity(19)
      REAL halflife(19)
      REAL*8 sourceterm(19)
      REAL*8 initialamt(19)
      REAL*8 kvalue(19)
      REAL*8 ktot(19)

      INTEGER N, M, I, J, K, L, II, JJ, P, Q                        !loop counters
      INTEGER nchain                                !number of nuclides in the decay series
      REAL xtime
      REAL xremove
      REAL*8 capp
      REAL*8 capn0
      REAL*8 outprod
      REAL*8 xsum
      REAL*8 inprod
      REAL*8 z1
      REAL*8 term1
      REAL*8 term2
!      allocate(kvalue(19))
!      allocate(ktot(19))
!      allocate(sourceterm(19))
!      allocate(initialamt(19))
!      allocate(activity(19))
!      allocate(halflife(19))

!Calculate parent and progeny activity at the ground level


         do l = 1, nchain !loop over radionuclide and progeny

         kvalue(l) = 0.693417D0/halflife(l)         
         ktot(l) = kvalue(l) + xremove
         activity(l) = 0.0

         do jj = 1, l

          capp = sourceterm(jj) / kvalue(jj)
          capn0 = initialamt(jj) / kvalue(jj)

          if(capp.eq.0.and.capn0.eq.0) then
          !do nothing
          else
          outprod = 1.D0
           if (jj.lt.l) then
            do ii = jj, l-1
            outprod = outprod * kvalue(ii)
            enddo
           endif
          endif

          if(outprod.gt.1.0D36) outprod = 1.D36


          xsum = 0.0D0

           do p = jj, l

           inprod = 1.0D0

            do q = jj, l
            if(q.ne.p) then
             inprod = inprod * (ktot(q)-ktot(p))
            endif
            enddo !end loop q

           z1 = ktot(p) * xtime
           if(z1.gt.40.D0) z1 = 40.D0
           term1 = capn0 * exp(-z1)/inprod
           term2 = capp * ( 1.0D0 - exp(-z1))/(inprod*ktot(p))
           xsum = xsum + term1 + term2

           enddo !end loop p

          if(xsum.lt.1.0D-36) xsum = 0.0D0
          activity(l) = activity(l) + outprod * xsum

          enddo !end loop jj

         activity(l) = activity(l) * kvalue(l)
         if(activity(l).lt.1.0D-36) activity(l) = 0.0D0

         enddo !end loop l



      RETURN

      END


