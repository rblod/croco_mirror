! $Id: wetdry.F 1408 2013-12-20 12:41:14Z marchesiello $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!
! This routine was adapted from John Warner's code
! by P. Marchesiello, IRD 2013
!======================================================================
!
#include "cppdefs.h"
#ifdef WET_DRY
!***********************************************************************
      SUBROUTINE wetdry_tile (Istr,Iend,Jstr,Jend)
!***********************************************************************
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "coupling.h"
!
      real   cff, eps, cff1,cff2,cff3
      parameter (eps=1.e-10)
# ifdef AGRIF
#  include "zoom.h"
      external rmask_wetinterp,umask_wetinterp,vmask_wetinterp
      real :: tinterp
# endif
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
!
# include "compute_auxiliary_bounds.h"
!
!-----------------------------------------------------------------------
! If wet/drying, compute new masks for cells with depth < Dcrit.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
      IF (iif.le.nfast) THEN
# endif
!
!  Wet/dry mask at RHO-points.
!
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
            wetdry(i,j)=1.
# ifdef MASKING
            wetdry(i,j)=wetdry(i,j)*rmask(i,j)
# endif
            IF ((zeta(i,j,knew)+h(i,j)).le.Dcrit(i,j)+eps) THEN
              wetdry(i,j)=0.
            END IF
# if defined AGRIF && defined AGRIF_2WAY
            IF ((wetdry(i,j)==1) .AND. (rmask_childs(i,j)==0.)) THEN
              wetdry(i,j) = 0.
            END IF
# endif
          END DO
        END DO

# ifdef AGRIF0
        if (.not.agrif_root()) then
C$OMP BARRIER
C$OMP MASTER
          tinterp = 1.
#  ifdef MASKING
          Agrif_UseSpecialValue = .false.
#  endif
          call agrif_bc_variable(rmask_wetid,calledweight=tinterp,
     &                           procname = rmask_wetinterp)
C$OMP END MASTER
C$OMP BARRIER
        endif
# endif

        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rmask_wet(i,j)=wetdry(i,j)
          END DO
        END DO
!
!  Wet/dry mask at PSI-points.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            pmask_wet(i,j)=wetdry(i-1,j  )*wetdry(i  ,j  )*
     &                     wetdry(i-1,j-1)*wetdry(i  ,j-1)
          END DO
        END DO
!
!  Wet/dry mask at U-points.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            umask_wet(i,j)=wetdry(i-1,j)+wetdry(i,j)
            IF (umask_wet(i,j).eq.1.) THEN
              umask_wet(i,j)=wetdry(i-1,j)-wetdry(i,j)
            END IF
          END DO
        END DO

# ifdef AGRIF0
        if (.not.agrif_root()) then
C$OMP BARRIER
C$OMP MASTER
          tinterp = 1.
#  ifdef MASKING
        Agrif_UseSpecialValue = .false.
#  endif
          call agrif_bc_variable(umask_wetid,calledweight=tinterp,
     &                           procname = umask_wetinterp)
C$OMP END MASTER
C$OMP BARRIER
        endif
# endif

!
!  Wet/dry mask at V-points.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vmask_wet(i,j)=wetdry(i,j-1)+wetdry(i,j)
            IF (vmask_wet(i,j).eq.1.) THEN
              vmask_wet(i,j)=wetdry(i,j-1)-wetdry(i,j)
            END IF
          END DO
        END DO

# ifdef AGRIF0
        if (.not.agrif_root()) then
C$OMP BARRIER
C$OMP MASTER
          tinterp = 1.
#  ifdef MASKING
        Agrif_UseSpecialValue = .false.
#  endif
          call agrif_bc_variable(vmask_wetid,calledweight=tinterp,
     &                           procname = vmask_wetinterp)
C$OMP END MASTER
C$OMP BARRIER
        endif
# endif

# ifdef SOLVE3D
      END IF
# endif
!
!--------------------------------------------------------------------
!  Averaged wet/dry mask for 3D slow mode if M2 filter
!--------------------------------------------------------------------
!
# if defined SOLVE3D && !defined M2FILTER_NONE
!
!  Wet/dry mask at RHO-points, averaged over all fast time-steps.
!
      IF (iif.le.nfast) THEN
        IF (FIRST_2D_STEP) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              rmask_wet_avg(i,j)=wetdry(i,j)
            END DO
          END DO
        ELSE
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              rmask_wet_avg(i,j)=rmask_wet_avg(i,j)+wetdry(i,j)
            END DO
          END DO
        END IF
      END IF
!
!  If done fast time-stepping, scale mask by 2 nfast.
!
      IF (iif.eq.nfast+1) THEN
        cff=1./REAL(nfast)
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
            wetdry(i,j)=AINT(rmask_wet_avg(i,j)*cff)
          END DO
        END DO
!
!  Wet/dry mask at RHO-points, averaged over all fast time-steps.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rmask_wet(i,j)=wetdry(i,j)
          END DO
        END DO
!
!  Wet/dry mask at PSI-points, averaged over all fast time-steps.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            pmask_wet(i,j)=wetdry(i-1,j  )*wetdry(i  ,j  )*
     &                     wetdry(i-1,j-1)*wetdry(i  ,j-1)
          END DO
        END DO
!
!  Wet/dry mask at U-points, averaged over all fast time-steps.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff1=wetdry(i-1,j)+wetdry(i,j)
            IF (cff1.eq.1.) THEN
              cff1=wetdry(i-1,j)-wetdry(i,j)
            END IF
            cff2=ABS(ABS(cff1)-1.)
            cff3=0.5+SIGN(0.5,DU_avg1(i,j,nnew))*cff1
            umask_wet(i,j)=0.5*cff1*cff2+cff3*(1.-cff2)
!                                                    catch lone ponds
            IF (DU_avg1(i,j,nnew).eq.0. .and.
     &          wetdry(i-1,j)+wetdry(i,j).le.1.) THEN
              umask_wet(i,j)=0.
            END IF
          END DO
        END DO

#  ifdef AGRIF0
        if (.not.agrif_root()) then
C$OMP BARRIER
C$OMP MASTER
          tinterp = 1.
#   ifdef MASKING
        Agrif_UseSpecialValue = .false.
#   endif
          call agrif_bc_variable(umask_wetid,calledweight=tinterp,
     &                           procname = umask_wetinterp)
C$OMP END MASTER
C$OMP BARRIER
        endif
#  endif

!
!  Wet/dry mask at V-points, averaged over all fast time-steps.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff1=wetdry(i,j-1)+wetdry(i,j)
            IF (cff1.eq.1.) THEN
              cff1=wetdry(i,j-1)-wetdry(i,j)
            END IF
            cff2=ABS(ABS(cff1)-1.)
            cff3=0.5+SIGN(0.5,DV_avg1(i,j,nnew))*cff1
            vmask_wet(i,j)=0.5*cff1*cff2+cff3*(1.-cff2)
!                                                    catch lone ponds
            IF (DV_avg1(i,j,nnew).eq.0. .and.
     &          wetdry(i,j-1)+wetdry(i,j).le.1.) THEN
              vmask_wet(i,j)=0.
            END IF
          END DO
        END DO

#  ifdef AGRIF0
        if (.not.agrif_root()) then
C$OMP BARRIER
C$OMP MASTER
          tinterp = 1.
#   ifdef MASKING
        Agrif_UseSpecialValue = .false.
#   endif
          call agrif_bc_variable(vmask_wetid,calledweight=tinterp,
     &                           procname = vmask_wetinterp)
C$OMP END MASTER
C$OMP BARRIER
        endif
#  endif

      END IF
# endif
!
!  Exchange boundary data.
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,  rmask_wet)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,  rmask_wet_avg)
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,  umask_wet)
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,  vmask_wet)
      call exchange_p2d_tile (Istr,Iend,Jstr,Jend,  pmask_wet)
# endif

# ifdef WETDRY_IO
!
!-----------------------------------------------------------------------
!  Set masks for I/O purposes.
!-----------------------------------------------------------------------
!
#  ifdef SOLVE3D
      IF (iif.gt.nfast) THEN
#  else
      IF (iif.eq.nfast) THEN
#  endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rmask_io(i,j)=rmask_wet(i,j) SWITCH rmask(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=Istr,IendR
            pmask_io(i,j)=pmask_wet(i,j) SWITCH pmask(i,j)
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            umask_io(i,j)=umask_wet(i,j) SWITCH umask(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vmask_io(i,j)=vmask_wet(i,j) SWITCH vmask(i,j)
          END DO
        END DO

#  ifdef PSOURCE
!
!  Insure that masks at mass point source locations are set to water
!  to avoid writting output with FillValue at those locations.
!
        DO is=1,Nsrc
          i=Isrc(is)
          j=Jsrc(is)
          IF (((IstrR.le.i).and.(i.le.IendR)).and.
     &        ((JstrR.le.j).and.(j.le.JendR))) THEN
            IF (INT(Dsrc(is)).eq.0) THEN
              umask_io(i,j)=1.
            ELSE
              vmask_io(i,j)=1.
            END IF
          END IF
        END DO
#  endif

      END IF
# endif /* WETDRY_IO */

      RETURN
      END SUBROUTINE wetdry_tile

# ifdef AGRIF
      subroutine rmask_wetinterp(tabres,i1,i2,j1,j2,before)
      implicit none
#  include "param.h"
#  include "ocean2d.h"
#  include "scalars.h"
#  include "zoom.h"
#  include "grid.h"

       integer i1,i2,j1,j2
       real tabres(i1:i2,j1:j2)
       logical before
       real, parameter :: eps=1.1e-10
       integer :: i,j

       if (before) then
          tabres(i1:i2,j1:j2) = wetdry(i1:i2,j1:j2)
        else
          wetdry(i1:i2,j1:j2) = tabres
          do j=j1,j2
          do i=i1,i2
            if (wetdry(i,j) == 1) then
            zeta(i,j,knew) = max(zeta(i,j,knew),Dcrit(i,j)+eps-h(i,j))
            endif
          enddo
          enddo
       endif
      return
      end

      subroutine umask_wetinterp(tabres,i1,i2,j1,j2,before)
      implicit none
#  include "param.h"
#  include "ocean2d.h"
#  include "scalars.h"
#  include "zoom.h"
#  include "grid.h"

       integer i1,i2,j1,j2
       real tabres(i1:i2,j1:j2)
       logical before

       if (before) then
          tabres(i1:i2,j1:j2) = umask_wet(i1:i2,j1:j2)
        else
          umask_wet(i1:i2,j1:j2) = tabres
       endif
      return
      end

      subroutine vmask_wetinterp(tabres,i1,i2,j1,j2,before)
      implicit none
#  include "param.h"
#  include "ocean2d.h"
#  include "scalars.h"
#  include "zoom.h"
#  include "grid.h"

       integer i1,i2,j1,j2
       real tabres(i1:i2,j1:j2)
       logical before

       if (before) then
          tabres(i1:i2,j1:j2) = vmask_wet(i1:i2,j1:j2)
        else
          vmask_wet(i1:i2,j1:j2) = tabres
       endif
      return
      end
# endif

#endif

