! $Id: ab_ratio.F 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#ifdef LMD_MIXING

      subroutine ab_ratio_slice (ratio, Istr,Iend, j)
!
!--------------------------------------------------------------------
!   This subroutine calculates the ratio of the thermodynamic
!  expansion coefficients for potential temperature and salinity,
!  alpha/beta, at horizontal and vertical W-points from a polynomial
!  expression (Jackett and McDougall,  1992).
!   The polynomial expression was found from fitting to248 values in
!  the oceanographic ranges of salinity,potential temperature, and
!  pressure. It assumes no pressure variation along geopotential
!  surfaces, that is, depth and pressure are interchangeable. The
!  RMS error of this polynomial is 0.000894.
!
!  On Output:
!    ab_ratio  ratio of expansion coefficients alpha (temperature)
!              over beta (salinity), (Celsius/PSU).
!
!  Check Value:
!     ab_ratio=0.34763 (Celsius/PSU)  at s=40.0, t=10.0, p=4000.0
!--------------------------------------------------------------------
!
      implicit none
# include "param.h"
      integer Istr,Iend, j, i,k
      real ratio(PRIVATE_1D_SCRATCH_ARRAY,0:N)
# ifndef NONLIN_EOS
     &   , cff
# endif
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"
# ifdef NONLIN_EOS
      real A0, A1, A2, A3, A4, B0, B1, CO, D0, D1, D2,
     &                         E0, F0, G0, H0, Smean,  Tt, Ts, Tp
      parameter(A0=+0.665157E-01, A1=+0.170907E-01, A2=-0.203814E-03,
     &          A3=+0.298357E-05, A4=-0.255019E-07, B0=+0.378110E-02,
     &          B1=-0.846960E-04, CO=-0.678662E-05, D0=+0.380374E-04,
     &          D1=-0.933746E-06, D2=+0.791325E-08, E0=-0.164759E-06,
     &          F0=-0.251520E-11, G0=+0.512857E-12, H0=-0.302285E-13,
     &          Smean=35.0)
#  define tind nstp
!
!  Compute the ratio of thermal expansion
!  and saline contraction coefficients.
!
!  Nonlinear Equation of state.  The units are as follows:
!------------------------------
!     Ts        salinity (PSU) anomaly from Smean.
!     Tt        potential temperature (degC).
!     Tp        pressure/depth, (depth in meters and positive).
!
      do k=1,N-1
        do i=Istr,Iend
#  ifdef TEMPERATURE
          Tt=0.5*(t(i,j,k,tind,itemp)+t(i,j,k+1,tind,itemp))
#  else
          Tt=25.0                                    !  OF STATE
#  endif
#  ifdef SALINITY
          Ts=0.5*(t(i,j,k,tind,isalt)+t(i,j,k+1,tind,isalt))
     &                                                -Smean
#  else
          Ts=0.
#  endif
          Tp=-z_w(i,j,k)
!-
          ratio(i,k)=A0+Tt*(A1+Tt*(A2+Tt*(A3+Tt*A4)))
     &                 +Ts*(B0+Tt*B1+Ts*CO)
     &                 +Tp*(D0+Tt*(D1+Tt*D2)+Ts*E0
     &                        +Tp*(Ts*F0+Tt*Tt*G0+Tp*H0))
        enddo
      enddo    /* R8000: 16 clock cycles; 56%(75%) of peak */
# else                        /* Linear Equation of state. */
#  ifdef SALINITY             /* ------------------------- */
      if (Scoef.ne.0.) then
        cff=Tcoef/Scoef
      else
        cff=Tcoef
      endif
#  else
      cff=Tcoef
#  endif
      do k=1,N-1
        do i=Istr,Iend
          ratio(i,k)=cff
        enddo
      enddo
# endif /* NONLIN_EOS */
#else
      subroutine ab_ratio_empt
#endif /* LMD_MIXING */
      return
      end

