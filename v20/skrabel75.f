      SUBROUTINE skrable( activity, sourceterm, xtime, nchain,
     &                    xremove, halflife, initialamt )

      IMPLICIT NONE
 
      ! Use consistent double precision for all reals here
      REAL*8 activity(19)
      REAL*8 halflife(19)
      REAL*8 sourceterm(19)
      REAL*8 initialamt(19)
      REAL*8 kvalue(19)
      REAL*8 ktot(19)

      INTEGER nchain
      INTEGER l, jj, ii, p, q
      REAL*8 xtime
      REAL*8 xremove
      REAL*8 capp, capn0, outprod, xsum, inprod, z1, term1, term2

      REAL*8, PARAMETER :: LOG2 = 0.693147180559945309417232121458176568D0
      REAL*8, PARAMETER :: EPS_DEN = 1.0D-30          ! small regularization for denominator
      REAL*8, PARAMETER :: EXP_CAP = 40.0D0           ! cap for exponent to avoid underflow
      REAL*8, PARAMETER :: TINY_ACT = 1.0D-36         ! activity floor

      ! This solver only works if none of the half lives are exactly the same
      do l = 1, nchain
      do q = 1, nchain
       if(l.ne.q) then !same index same nuclide
       if(halflife(l).eq.halflife(q)) then
       halflife(l) = halflife(l)*1.0001 !make a minor change to the halflife 
       else
       endif
       else
       endif
      enddo
      enddo

      ! Precompute decay and total removal constants
      do l = 1, nchain
       kvalue(l)   = LOG2 / halflife(l)
       ktot(l)     = kvalue(l) + xremove
       activity(l) = 0.0D0
      enddo

      ! Generalized Bateman with constant sources
      do l = 1, nchain

      do jj = 1, l
        
        capp  = sourceterm(jj) / kvalue(jj)
        capn0 = initialamt(jj) / kvalue(jj)

        ! If both source and initial are zero, no contribution from nuclide jj -> skip safely
        if (capp == 0.0D0 .and. capn0 == 0.0D0) cycle

        ! Product of decay constants from jj to l-1 (empty product = 1)
        outprod = 1.0D0
        if (jj < l) then
          do ii = jj, l-1
            outprod = outprod * kvalue(ii)
          enddo
        endif
        if (outprod > 1.0D36) outprod = 1.0D36

        xsum = 0.0D0

        do p = jj, l

          ! Denominator product Π_{q=jj..l, q≠p} (ktot(q) - ktot(p))
          inprod = 1.0D0
          do q = jj, l
            if (q /= p) then
              inprod = inprod * (ktot(q) - ktot(p))
            endif
          enddo

          ! Guard against exact or near-equal ktot leading to zero/ill-conditioned denominator
          if (abs(inprod) < EPS_DEN) then
            inprod = sign(EPS_DEN, inprod)
          endif

          z1 = ktot(p) * xtime
          if (z1 > EXP_CAP) z1 = EXP_CAP

          term1 = capn0 * exp(-z1) / inprod
          term2 = capp  * (1.0D0 - exp(-z1)) / (inprod * ktot(p))

          xsum = xsum + term1 + term2

        enddo ! p

        if (abs(xsum) < TINY_ACT) xsum = 0.0D0
        activity(l) = activity(l) + outprod * xsum

      enddo ! jj

      activity(l) = activity(l) * kvalue(l)
      if (abs(activity(l)) < TINY_ACT) activity(l) = 0.0D0

      enddo ! l

      RETURN
      END SUBROUTINE skrable

