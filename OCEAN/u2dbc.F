! $Id: u2dbc.F 1409 2014-01-06 16:34:49Z marchesiello $
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
#ifndef CHILD
!
# include "cppdefs.h"
      subroutine u2dbc_tile(Istr,Iend,Jstr,Jend,grad)
# ifdef AGRIF
      use AGRIF_Util
      integer Istr,Iend,Jstr,Jend
      real grad(PRIVATE_2D_SCRATCH_ARRAY)

      if (AGRIF_Root()) then
        call u2dbc_parent_tile(Istr,Iend,Jstr,Jend,grad)
      else
        call u2dbc_child_tile(Istr,Iend,Jstr,Jend,grad)
      endif
      return
      end
!
! PARENT
!
      subroutine u2dbc_parent_tile (Istr,Iend,Jstr,Jend,grad)
# endif
!
! Set lateral boundary conditions for the barotropic (i.e.
! vertically integrated) XI-component velocity ubar(:,:,knew)
! for the parent grid.
!
#else
!
! CHILD
!
      subroutine u2dbc_child_tile (Istr,Iend,Jstr,Jend,grad)
!
! Set lateral boundary conditions for the barotropic (i.e.
! vertically integrated) XI-component velocity ubar(:,:,knew)
! for the child grid.
!
#endif /* CHILD */
!
! Common Code
!
#include "set_obc_definitions.h"
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j
      real    grad(PRIVATE_2D_SCRATCH_ARRAY)
#ifdef AGRIF
# include "private_scratch.h"
      integer trd
C$    integer omp_get_thread_num
#endif
      real    eps,cff, cx,cy,u_str,
     &        dft,dfx,dfy, tau,tau_in,tau_out,hx,zx
#ifdef OBC_REDUCED_PHYSICS
     &        ,bry_val,bry_pgr,bry_str,bry_cor
#endif
      parameter (eps=1.D-20)
#if defined WET_DRY || defined OBC_REDUCED_PHYSICS
      real cff1,cff2
#endif
#include "param.h"
#include "boundary.h"
#include "climat.h"
#include "grid.h"
#include "ocean2d.h"
#include "scalars.h"
#include "forces.h"
#include "zoom.h"
#ifdef NBQ
# include "nbq.h"
#endif
!
#include "compute_auxiliary_bounds.h"
!
!  If child grid, interpolate ubar from parent
!  AGRIF does not work yet with #undef SOLVE3D
!
#ifdef CHILD
# ifdef SOLVE3D
#  ifdef NBQ
      if (M2bc_nbq_flag) then
#  endif
      trd=0
C$    trd=omp_get_thread_num()
      call u2dbc_interp_tile(Istr,Iend,Jstr,Jend
     &                    ,A1dETA(1,5+3*NWEIGHT)
     &                    ,A1dETA(1,6+4*NWEIGHT)
     &                    ,A1dETA(1,7+5*NWEIGHT)
     &                    ,A1dETA(1,8+5*NWEIGHT)
     &                    ,A1dXI(1,5+3*NWEIGHT)
     &                    ,A1dXI(1,6+4*NWEIGHT)
     &                    ,A1dXI(1,7+5*NWEIGHT)
     &                    ,A1dXI(1,8+5*NWEIGHT))
#  ifdef NBQ
       endif
#  endif
# else
      Some coding is needed to have AGRIF without SOLVE3D.....
# endif
#endif
!
#if defined M2_FRC_BRY || defined M2CLIMATOLOGY
      tau_in=dtfast*tauM_in
      tau_out=dtfast*tauM_out
#endif
!
#if defined OBC_COM_M2CHARACT || defined OBC_COM_M2CHARACT_EAST  || defined OBC_COM_M2CHARACT_WEST
      grad = 1.
# if defined OBC_COM_SOUTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        grad(Istr,Jstr) = 0.5
      endif
# endif
# if defined OBC_COM_SOUTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        grad(Iend+1,Jstr) = 0.5
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        grad(Istr,Jend) = 0.5
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        grad(Iend+1,Jend) = 0.5
      endif
# endif
#endif

!
!====================================================================
!                            WESTERN BC
!====================================================================
#ifndef EW_COM_PERIODIC
      if (WESTERN_EDGE) then
# ifdef OBC_COM_WEST
#  ifdef OBC_COM_M2ORLANSKI
!                                           Western edge radiation BC
!                                           ======= ==== ========= ==
        do j=Jstr,Jend+1
          grad(Istr  ,j)=(ubar(Istr  ,j,kstp)-ubar(Istr  ,j-1,kstp))
#   ifdef MASKING
     &                                                *pmask(Istr,j)
#   endif
          grad(Istr+1,j)=(ubar(Istr+1,j,kstp)-ubar(Istr+1,j-1,kstp))
#   ifdef MASKING
     &                                              *pmask(Istr+1,j)
#   endif
        enddo
        do j=Jstr,Jend
          dft=ubar(Istr+1,j,kstp)-ubar(Istr+1,j,knew)
          dfx=ubar(Istr+1,j,knew)-ubar(Istr+2,j,knew)

          if (dfx*dft .lt. 0.) then
            dft=0.                       ! <-- cancel cx, if inflow
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
            tau=tau_in
          else
            tau=tau_out
#   endif
          endif

          if (dft*(grad(Istr+1,j)+grad(Istr+1,j+1)) .gt. 0.) then
            dfy=grad(Istr+1,j  )
          else
            dfy=grad(Istr+1,j+1)
          endif

#   ifdef OBC_COM_RAD_NORMAL
          dfy=0.
#   endif
          cff=max(dfx*dfx+dfy*dfy, eps)
          cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
          cy=0.
#   else
          cy=min(cff,max(dft*dfy,-cff))
#   endif

          ubar(Istr,j,knew)=( cff*ubar(Istr  ,j,kstp)
     &                        +cx*ubar(Istr+1,j,knew)
     &                     -max(cy,0.)*grad(Istr,j  )
     &                     -min(cy,0.)*grad(Istr,j+1)
     &                                     )/(cff+cx)
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
          ubar(Istr,j,knew)=(1.-tau)*ubar(Istr,j,knew)
#    ifdef M2_FRC_BRY
     &                            +tau*ubarbry_west(j)
#    else
     &                            +tau*ubclm(Istr,j)
#    endif
#   endif
#   ifdef MASKING
          ubar(Istr,j,knew)=ubar(Istr,j,knew)*umask(Istr,j)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT
!                                            Western edge Characteristic BC
!                                            ======= ==== ======= ==
        do j=Jstr,Jend
          cff=0.5*(h(Istr-1,j)+zeta(Istr-1,j,kstp)+
     &             h(Istr  ,j)+zeta(Istr  ,j,kstp))
          hx=sqrt(g/cff)
          cx=dtfast*cff*hx*0.5*(pm(Istr-1,j)+pm(Istr,j))
          zx=(0.5+cx)*zeta(istr,j,kstp)+(0.5-cx)*zeta(istr-1,j,kstp)
          if (cx .gt. 0.292893218813452) then
            zx=zx + ( zeta(istr,j,knew) +cx*zeta(istr-1,j,kstp)
     &                             -(1.+cx)*zeta(istr  ,j,kstp)
     &                           )*(1.-0.292893218813452/cx)**2
          endif
#   ifdef OBC_REDUCED_PHYSICS
          cff1=1.0/(0.5*(h(Istr-1,j)+zeta(Istr-1,j,kstp)+
     &                   h(Istr  ,j)+zeta(Istr  ,j,kstp)))
#    ifdef STRESS_AT_RHO_POINTS
          u_str=0.5*(sustr(Istr,j)+sustr(Istr-1,j))
#    else
          u_str=sustr(Istr,j)
#    endif

          bry_str=cff1*(u_str-bustr(Istr,j))
#    ifdef UV_COR
          bry_cor=0.125*(vbar(Istr-1,j  ,kstp)+
     &                   vbar(Istr-1,j+1,kstp)+
     &                   vbar(Istr  ,j  ,kstp)+
     &                   vbar(Istr  ,j+1,kstp))*
     &                  (f(Istr-1,j)+f(Istr,j))
#    else
          bry_cor=0.
#    endif
          bry_pgr=-g*(zeta(Istr,j,kstp)-zeta(Istr-1,j,kstp))*
     &                          0.5*(pm(Istr-1,j)+pm(Istr,j))

          Cx=1.0/SQRT(g*0.5*(h(Istr-1,j)+zeta(Istr-1,j,kstp)+
     &                       h(Istr  ,j)+zeta(Istr  ,j,kstp)))
          cff2=om_u(Istr,j)*Cx
          bry_val=ubar(Istr+1,j,kstp)+
     &                cff2*(bry_pgr+
     &                      bry_cor+
     &                      bry_str)
#   endif
          ubar(Istr,j,knew)= 0.5*( (1.-cx)*ubar(Istr  ,j,kstp)
     &                                 +cx*ubar(Istr+1,j,kstp)
#   ifdef M2_FRC_BRY
#    ifdef OBC_REDUCED_PHYSICS
     &                                     +bry_val
#    else
     &                                     +ubarbry_west(j)
#    endif
#   elif defined M2CLIMATOLOGY
     &                                     +ubclm(Istr,j)
#   endif
     &                                     -hx*( zx
#   ifdef Z_FRC_BRY
     &                                        -zetabry_west(j)
#   elif defined ZCLIMATOLOGY
     &                                        -ssh(Istr-1,j)
#   endif
     &                           ))
#   ifdef MASKING
          ubar(Istr,j,knew)=ubar(Istr,j,knew)*umask(Istr,j)
#   endif
        enddo
#  elif defined OBC_COM_M2SPECIFIED || defined OBC_COM_M2SPECIFIED_WEST
!                                            Western edge Specified BC
!                                            ======= ==== ========= ==
        do j=Jstr,Jend
#   ifdef M2_FRC_BRY
          ubar(Istr,j,knew)=ubarbry_west(j)
#   else
          ubar(Istr,j,knew)=ubclm(Istr,j)
#   endif
#   ifdef MASKING
     &                      *umask(Istr,j)
#   endif
        enddo
#  else
!                                            Western edge gradient BC
!                                            ======= ==== ======== ==
        do j=Jstr,Jend
          ubar(Istr,j,knew)=ubar(Istr+1,j,knew)
#   ifdef MASKING
     &                      *umask(Istr,j)
#   endif
        enddo

#  endif /* OBC_COM_ORLASNKI */
!
# else
                                                ! Western edge closed
        do j=jstr,jend                          ! ======= ==== ======
#  ifdef MRL_WCI
          ubar(istr,j,knew)=-ust2d(Istr,j)      ! no Lagrangian flux
#   ifdef MASKING
     &                      *umask(Istr,j)
#   endif
#  else
          ubar(istr,j,knew)=0.                  ! (no-flux, default)
#  endif
        enddo

# endif           /* OBC_COM_WEST */

      endif     !<-- WESTERN_EDGE
!
!====================================================================
!                            EASTERN BC
!====================================================================
      if (EASTERN_EDGE) then
# ifdef OBC_COM_EAST
#  ifdef OBC_COM_M2ORLANSKI
!                                           Eastern edge radiation BC
!                                           ======= ==== ========= ==
        do j=Jstr,Jend+1
          grad(Iend  ,j)=(ubar(Iend  ,j,kstp)-ubar(Iend  ,j-1,kstp))
#   ifdef MASKING
     &                                                *pmask(Iend,j)
#   endif
          grad(Iend+1,j)=(ubar(Iend+1,j,kstp)-ubar(Iend+1,j-1,kstp))
#   ifdef MASKING
     &                                              *pmask(Iend+1,j)
#   endif
        enddo
        do j=Jstr,Jend
          dft=ubar(Iend,j,kstp)-ubar(Iend  ,j,knew)
          dfx=ubar(Iend,j,knew)-ubar(Iend-1,j,knew)

          if (dfx*dft .lt. 0.) then
            dft=0.                       ! <-- cancel cx, if inflow
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
            tau=tau_in
          else
            tau=tau_out
#   endif
          endif

          if (dft*(grad(Iend,j)+grad(Iend,j+1)) .gt. 0.) then
            dfy=grad(Iend,j)
          else
            dfy=grad(Iend,j+1)
          endif

#   ifdef OBC_COM_RAD_NORMAL
          dfy=0.
#   endif
          cff=max(dfx*dfx+dfy*dfy, eps)
          cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
          cy=0.
#   else
          cy=min(cff,max(dft*dfy,-cff))
#   endif

          ubar(Iend+1,j,knew)=( cff*ubar(Iend+1,j,kstp)
     &                          +cx*ubar(Iend  ,j,knew)
     &                     -max(cy,0.)*grad(Iend+1,j  )
     &                     -min(cy,0.)*grad(Iend+1,j+1)
     &                                       )/(cff+cx)
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
          ubar(Iend+1,j,knew)=(1.-tau)*ubar(Iend+1,j,knew)
#    ifdef M2_FRC_BRY
     &                                +tau*ubarbry_east(j)
#    else
     &                                +tau*ubclm(Iend+1,j)
#    endif
#   endif
#   ifdef MASKING
          ubar(Iend+1,j,knew)=ubar(Iend+1,j,knew)*umask(Iend+1,j)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT || defined OBC_COM_M2CHARACT_EAST
!                                            Eastern edge Characteristic BC
!                                            ======= ==== ======= ==
        do j=Jstr,Jend
          cff=0.5*(h(Iend  ,j)+zeta(Iend  ,j,kstp)+
     &             h(Iend+1,j)+zeta(Iend+1,j,kstp))
          hx=sqrt(g/cff)
          cx=dtfast*cff*hx*0.5*(pm(Iend,j)+pm(Iend+1,j))
          zx=(0.5+cx)*zeta(iend,j,kstp)+(0.5-cx)*zeta(iend+1,j,kstp)
          if (cx .gt. 0.292893218813452) then
            zx=zx + ( zeta(iend,j,knew) +cx*zeta(iend+1,j,kstp)
     &                             -(1.+cx)*zeta(iend  ,j,kstp)
     &                           )*(1.-0.292893218813452/cx)**2
          endif
#   ifdef OBC_REDUCED_PHYSICS
          cff1=1.0/(0.5*(h(Iend  ,j)+zeta(Iend  ,j,kstp)+
     &                   h(Iend+1,j)+zeta(Iend+1,j,kstp)))
#    ifdef STRESS_AT_RHO_POINTS
          u_str=0.5*(sustr(Iend,j)+sustr(Iend+1,j))
#    else
          u_str=sustr(Iend+1,j)
#    endif
          bry_str=cff1*(u_str-bustr(Iend+1,j))
#    ifdef UV_COR
          bry_cor=0.125*(vbar(Iend  ,j  ,kstp)+
     &                   vbar(Iend  ,j+1,kstp)+
     &                   vbar(Iend+1,j  ,kstp)+
     &                   vbar(Iend+1,j+1,kstp))*
     &                   (f(Iend,j)+f(Iend+1,j))
#    else
          bry_cor=0.
#    endif
          bry_pgr=-g*(zeta(Iend+1,j,kstp)-zeta(Iend,j,kstp))*
     &                          0.5*(pm(Iend,j)+pm(Iend+1,j))
          Cx=1.0/SQRT(g*0.5*(h(Iend  ,j)+zeta(Iend  ,j,kstp)+
     &                       h(Iend+1,j)+zeta(Iend+1,j,kstp)))
          cff2=om_u(Iend+1,j)*Cx
          bry_val=ubar(Iend,j,kstp)+
     &                cff2*(bry_pgr+
     &                      bry_cor+
     &                      bry_str)
#   endif
#   ifdef MASKING
          ubar(Iend+1,j,knew)=ubar(Iend+1,j,knew)*umask(Iend+1,j)
#   endif
          ubar(Iend+1,j,knew)= 0.5*( (1.-cx)*ubar(Iend+1,j,kstp)
     &                                   +cx*ubar(Iend  ,j,kstp)
#   ifdef M2_FRC_BRY
#    ifdef OBC_REDUCED_PHYSICS
     &                                       +bry_val
#    else
     &                                       +ubarbry_east(j)
#    endif
#   elif defined M2CLIMATOLOGY
     &                                       +ubclm(Iend+1,j)
#   endif
     &                                       +hx*( zx
#   ifdef Z_FRC_BRY
     &                                         -zetabry_east(j)
#   elif defined ZCLIMATOLOGY
     &                                         -ssh(Iend,j)
#   endif
     &                             ))
        enddo
#  elif defined OBC_COM_M2SPECIFIED
!                                            Eastern edge Specified BC
!                                            ======= ==== ========= ==
        do j=Jstr,Jend
#   ifdef M2_FRC_BRY
          ubar(Iend+1,j,knew)=ubarbry_east(j)
#   else
          ubar(Iend+1,j,knew)=ubclm(Iend+1,j)
#   endif
#   ifdef MASKING
     &                        *umask(Iend+1,j)
#   endif
        enddo
#  else
!                                            Eastern edge gradient BC
!                                            ======= ==== ======== ==
        do j=Jstr,Jend
          ubar(Iend+1,j,knew)=ubar(Iend,j,knew)
#   ifdef MASKING
     &                        *umask(Iend+1,j)
#   endif
        enddo

#  endif /* OBC_COM_ORLANKSI */
!
# else
                                             ! Eastern edge closed BC
        do j=Jstr,Jend                       ! ======= ==== ====== ==
#  ifdef MRL_WCI
          ubar(Iend+1,j,knew)=-ust2d(Iend+1,j)     ! no Lagrangian flux
#   ifdef MASKING
     &                        *umask(Iend+1,j)
#   endif
#  else
          ubar(Iend+1,j,knew)=0.
#  endif
        enddo

# endif /* OBC_COM_EAST */

      endif         !<-- EASTERN_EDGE

#endif              /* !EW_COM_PERIODIC */


#ifndef NS_COM_PERIODIC
!
!====================================================================
!                            SOUTHERN BC
!====================================================================
      if (SOUTHERN_EDGE) then
# ifdef OBC_COM_SOUTH
#  if defined OBC_COM_M2SPECIFIED || defined OBC_COM_M2SPECIFIED_SOUTH
!                                           Southern edge Specified BC
!                                           ======== ==== ========= ==
        do i=IstrU,Iend
#   ifdef M2_FRC_BRY
          ubar(i,Jstr-1,knew)=ubarbry_south(i)
#   else
          ubar(i,Jstr-1,knew)=ubclm(i,Jstr-1)
#   endif
#   ifdef MASKING
     &                         *umask(i,Jstr-1)
#   endif
        enddo
#  elif defined OBC_COM_M2ORLANSKI
!                                          Southern edge radiation BC
!                                          ======== ==== ========= ==
        do i=IstrU-1,Iend
          grad(i,Jstr-1)=ubar(i+1,Jstr-1,kstp)-ubar(i,Jstr-1,kstp)
          grad(i,Jstr  )=ubar(i+1,Jstr  ,kstp)-ubar(i,Jstr  ,kstp)
        enddo
        do i=IstrU,Iend
          dft=ubar(i,Jstr,kstp)-ubar(i,Jstr  ,knew)
          dfx=ubar(i,Jstr,knew)-ubar(i,Jstr+1,knew)

          if (dfx*dft .lt. 0.) then
            dft=0.                        ! <-- cancel cx, if inflow
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
            tau=tau_in
          else
            tau=tau_out
#   endif
          endif

          if (dft*(grad(i-1,Jstr)+grad(i,Jstr)) .gt. 0.) then
            dfy=grad(i-1,Jstr)
          else
            dfy=grad(i  ,Jstr)
          endif

#   ifdef OBC_COM_RAD_NORMAL
          dfy=0.
#   endif
          cff=max(dfx*dfx+dfy*dfy, eps)
          cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
          cy=0.
#   else
          cy=min(cff,max(dft*dfy,-cff))
#   endif

          ubar(i,Jstr-1,knew)=( cff*ubar(i,Jstr-1,kstp)
     &                          +cx*ubar(i,Jstr  ,knew)
     &                     -max(cy,0.)*grad(i-1,Jstr-1)
     &                     -min(cy,0.)*grad(i  ,Jstr-1)
     &                                       )/(cff+cx)
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
          ubar(i,Jstr-1,knew)=(1.-tau)*ubar(i,Jstr-1,knew)
#    ifdef M2_FRC_BRY
     &                               +tau*ubarbry_south(i)
#    else
     &                               +tau*ubclm(i,Jstr-1)
#    endif
#   endif
#   ifdef MASKING
          ubar(i,Jstr-1,knew)=ubar(i,Jstr-1,knew)*umask(i,Jstr-1)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT
!                                           Southern edge Chapman BC
!                                           ======== ==== ======= ==
        do i=IstrU,Iend
          cff=sqrt(0.5*g*(h(i-1,Jstr-1)+zeta(i-1,Jstr-1,kstp)+
     &                    h(i  ,Jstr-1)+zeta(i  ,Jstr-1,kstp)))
          cx=dtfast*0.5*cff*(pn(i-1,Jstr-1)+pn(i,Jstr-1))
          ubar(i,Jstr-1,knew)=(ubar(i,Jstr-1,kstp)
     &                     +cx*ubar(i,Jstr  ,knew) )/(1.+cx)
#   ifdef MASKING
     &                         *umask(i,Jstr-1)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT0
!                                      Southern edge Characteristic BC
!                                      ======== ==== ============== ==
        do i=IstrU,Iend
          if (vbar(i,Jstr+1,knew).lt.0.) then
            ubar(i,Jstr-1,knew)=ubar(i,Jstr,knew)
          else
#   ifdef M2_FRC_BRY
            ubar(i,Jstr-1,knew)=ubarbry_south(i)
#   else
            ubar(i,Jstr-1,knew)=ubclm(i,Jstr-1)
#   endif
          endif
#   ifdef MASKING
          ubar(i,Jstr-1,knew)=ubar(i,Jstr-1,knew)*umask(i,Jstr-1)
#   endif
        enddo
#  else
!                                           Southern edge gradient BC
!                                           ======== ==== ======== ==
        do i=IstrU,Iend
          ubar(i,Jstr-1,knew)=ubar(i,Jstr,knew)
#   ifdef MASKING
     &                         *umask(i,Jstr-1)
#   endif
        enddo
#  endif
# else
#  ifdef EW_COM_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
                                        ! Wall: free-slip (gamma2=+1)
        do i=I_RANGE                    ! ====    no-slip (gamma2=-1)
          ubar(i,Jstr-1,knew)=gamma2*ubar(i,Jstr,knew)
#  ifdef MASKING
     &                                *umask(i,Jstr-1)
#  endif
        enddo
#  undef I_RANGE
# endif              /* OBC_COM_SOUTH */
      endif          !<-- SOUTHERN_EDGE
!
!====================================================================
!                            NORTHERN BC
!====================================================================
      if (NORTHERN_EDGE) then
# ifdef OBC_COM_NORTH
#  if defined OBC_COM_M2SPECIFIED || defined OBC_COM_M2SPECIFIED_NORTH
!                                           Northern edge Specified BC
!                                           ======== ==== ========= ==
        do i=IstrU,Iend
#   ifdef M2_FRC_BRY
          ubar(i,Jend+1,knew)=ubarbry_north(i)
#   else
          ubar(i,Jend+1,knew)=ubclm(i,Jend+1)
#   endif
#   ifdef MASKING
     &                        *umask(i,Jend+1)
#   endif
        enddo
#  elif defined OBC_COM_M2ORLANSKI
!                                          Northern edge radiation BC
!                                          ======== ==== ========= ==
        do i=IstrU-1,Iend
          grad(i,Jend  )=ubar(i+1,Jend  ,kstp)-ubar(i,Jend,kstp  )
          grad(i,Jend+1)=ubar(i+1,Jend+1,kstp)-ubar(i,Jend+1,kstp)
        enddo
        do i=IstrU,Iend
          dft=ubar(i,Jend,kstp)-ubar(i,Jend  ,knew)
          dfx=ubar(i,Jend,knew)-ubar(i,Jend-1,knew)

          if (dfx*dft .lt. 0.) then
            dft=0.                        ! <-- cancel cx, if inflow
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
            tau=tau_in
          else
            tau=tau_out
#   endif
          endif

          if (dft*(grad(i-1,Jend)+grad(i,Jend)) .gt. 0.) then
            dfy=grad(i-1,Jend)
          else
            dfy=grad(i  ,Jend)
          endif

#   ifdef OBC_COM_RAD_NORMAL
          dfy=0.
#   endif
          cff=max(dfx*dfx+dfy*dfy, eps)
          cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
          cy=0.
#   else
          cy=min(cff,max(dft*dfy,-cff))
#   endif

          ubar(i,Jend+1,knew)=( cff*ubar(i,Jend+1,kstp)
     &                          +cx*ubar(i,Jend  ,knew)
     &                     -max(cy,0.)*grad(i-1,Jend+1)
     &                     -min(cy,0.)*grad(i  ,Jend+1)
     &                                       )/(cff+cx)
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
          ubar(i,Jend+1,knew)=(1.-tau)*ubar(i,Jend+1,knew)
#    ifdef M2_FRC_BRY
     &                               +tau*ubarbry_north(i)
#    else
     &                               +tau*ubclm(i,Jend+1)
#    endif
#   endif
#   ifdef MASKING
          ubar(i,Jend+1,knew)=ubar(i,Jend+1,knew)*umask(i,Jend+1)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT
!                                           Northern edge Chapman BC
!                                           ======== ==== ======= ==
        do i=IstrU,Iend
          cff=sqrt(0.5*g*(h(i-1,Jend+1)+zeta(i-1,Jend+1,kstp)+
     &                    h(i  ,Jend+1)+zeta(i  ,Jend+1,kstp)))
          cx=dtfast*0.5*cff*(pn(i-1,Jend+1)+pn(i,Jend+1))
          ubar(i,Jend+1,knew)=(ubar(i,Jend+1,kstp)
     &                     +cx*ubar(i,Jend  ,knew))/(1.+cx)
#   ifdef MASKING
     &                        *umask(i,Jend+1)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT0
!                                    Northern edge Characteristic BC
!                                    ======== ==== ============== ==
        do i=IstrU,Iend
          if (vbar(i,Jend-1,knew).gt.0.) then
            ubar(i,Jend+1,knew)=ubar(i,Jend,knew)
          else
#   ifdef M2_FRC_BRY
            ubar(i,Jend+1,knew)=ubarbry_north(i)
#   else
            ubar(i,Jend+1,knew)=ubclm(i,Jend+1)
#   endif
          endif
#   ifdef MASKING
          ubar(i,Jend+1,knew)=ubar(i,Jend+1,knew)*umask(i,Jend+1)
#   endif
        enddo
#  else
!                                           Northern edge gradient BC
!                                           ======== ==== ======== ==
        do i=IstrU,Iend
          ubar(i,Jend+1,knew)=ubar(i,Jend,knew)
#   ifdef MASKING
     &                        *umask(i,Jend+1)
#   endif
        enddo

#  endif
# else
#  ifdef EW_COM_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
                                       !  Wall: free-slip (gamma2=+1)
        do i=I_RANGE                   !  =====   no-slip (gamma2=-1)
          ubar(i,Jend+1,knew)=gamma2*ubar(i,Jend,knew)
#  ifdef MASKING
     &                                *umask(i,Jend+1)
#  endif
        enddo
#  undef I_RANGE
# endif               /* OBC_COM_NORTH */
      endif          !<-- NORTHERN_EDGE
#endif                /* !NS_COM_PERIODIC */

!=======================================================================
                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========

#if defined OBC_COM_SOUTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        ubar(Istr,Jstr-1,knew)=0.5*( ubar(Istr+1,Jstr-1,knew)
     &                                  +ubar(Istr,Jstr,knew))
# ifdef MASKING
     &                        *umask(Istr,Jstr-1)
# endif
      endif
#endif
#if defined OBC_COM_SOUTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        ubar(Iend+1,Jstr-1,knew)=0.5*( ubar(Iend,Jstr-1,knew)
     &                                +ubar(Iend+1,Jstr,knew))
# ifdef MASKING
     &                        *umask(Iend+1,Jstr-1)
# endif
      endif
#endif
#if defined OBC_COM_NORTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        ubar(Istr,Jend+1,knew)=0.5*( ubar(Istr+1,Jend+1,knew)
     &                                  +ubar(Istr,Jend,knew))
# ifdef MASKING
     &                        *umask(Istr,Jend+1)
# endif
      endif
#endif
#if defined OBC_COM_NORTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        ubar(Iend+1,Jend+1,knew)=0.5*( ubar(Iend,Jend+1,knew)
     &                                +ubar(Iend+1,Jend,knew))
# ifdef MASKING
     &                        *umask(Iend+1,Jend+1)
# endif
      endif
#endif

!------------------------------------------------------------
#ifdef NBQ
       if (M2bc_nbq_flag) then
#endif
!------------------------------------------------------------

#ifdef WET_DRY
!=======================================================================
!  Impose wetting and drying conditions.
!=======================================================================
!
# ifndef EW_COM_PERIODIC
      if (WESTERN_EDGE) then
        DO j=Jstr,Jend
          cff1=ABS(ABS(umask_wet(Istr,j))-1.)
          cff2=0.5+SIGN(0.5,ubar(Istr,j,knew))*umask_wet(Istr,j)
          umask_wet(Istr,j)=0.5*umask_wet(Istr,j)*cff1
     &                                           +cff2*(1.-cff1)
          ubar(Istr,j,knew)=ubar(Istr,j,knew)*umask_wet(Istr,j)
        END DO
      END IF
      if (EASTERN_EDGE) then
        DO j=Jstr,Jend
          cff1=ABS(ABS(umask_wet(Iend+1,j))-1.)
          cff2=0.5+SIGN(0.5,ubar(Iend+1,j,knew))*umask_wet(Iend+1,j)
          umask_wet(Iend+1,j)=0.5*umask_wet(Iend+1,j)*cff1
     &                                               +cff2*(1.-cff1)
          ubar(Iend+1,j,knew)=ubar(Iend+1,j,knew)*umask_wet(Iend+1,j)
        END DO
      END IF
# endif
# ifndef NS_COM_PERIODIC
      if (SOUTHERN_EDGE) then
        DO i=IstrU,Iend
          cff1=ABS(ABS(umask_wet(i,Jstr-1))-1.)
          cff2=0.5+SIGN(0.5,ubar(i,Jstr-1,knew))*umask_wet(i,Jstr-1)
          umask_wet(i,Jstr-1)=0.5*umask_wet(i,Jstr-1)*cff1
     &                                               +cff2*(1.-cff1)
          ubar(i,Jstr-1,knew)=ubar(i,Jstr-1,knew)*umask_wet(i,Jstr-1)
        END DO
      END IF
      if (NORTHERN_EDGE) then
        DO i=Istr,Iend
          cff1=ABS(ABS(umask_wet(i,Jend+1))-1.)
          cff2=0.5+SIGN(0.5,ubar(i,Jend+1,knew))*umask_wet(i,Jend+1)
          umask_wet(i,Jend+1)=0.5*umask_wet(i,Jend+1)*cff1
     &                                               +cff2*(1.-cff1)
          ubar(i,Jend+1,knew)=ubar(i,Jend+1,knew)*umask_wet(i,Jend+1)
        END DO
      END IF
# endif
# if !defined EW_COM_PERIODIC && !defined NS_COM_PERIODIC
      if (SOUTHERN_EDGE .and. WESTERN_EDGE) then
        cff1=ABS(ABS(umask_wet(Istr,Jstr-1))-1.)
        cff2=0.5+SIGN(0.5,ubar(Istr,Jstr-1,knew))*umask_wet(Istr,Jstr-1)
        umask_wet(Istr,Jstr-1)=0.5*umask_wet(Istr,Jstr-1)*cff1
     &                                                   +cff2*(1.-cff1)
        ubar(Istr,Jstr-1,knew)=ubar(Istr,Jstr-1,knew)
     &                                        *umask_wet(Istr,Jstr-1)
      END IF
      if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
        cff1=ABS(ABS(umask_wet(Iend+1,Jstr-1))-1.)
        cff2=0.5+SIGN(0.5,ubar(Iend+1,Jstr-1,knew))
     &                                         *umask_wet(Iend+1,Jstr-1)
        umask_wet(Iend+1,Jstr-1)=0.5*umask_wet(Iend+1,Jstr-1)*cff1
     &                                                   +cff2*(1.-cff1)
        ubar(Iend+1,Jstr-1,knew)=ubar(Iend+1,Jstr-1,knew)
     &                                        *umask_wet(Iend+1,Jstr-1)
      END IF
      if (NORTHERN_EDGE .and. WESTERN_EDGE) then
        cff1=ABS(ABS(umask_wet(Istr,Jend+1))-1.)
        cff2=0.5+SIGN(0.5,ubar(Istr,Jend+1,knew))
     &                                       *umask_wet(Istr,Jend+1)
        umask_wet(Istr,Jend+1)=0.5*umask_wet(Istr,Jend+1)*cff1
     &                                                   +cff2*(1.-cff1)
        ubar(Istr,Jend+1,knew)=ubar(Istr,Jend+1,knew)
     &                                        *umask_wet(Istr,Jend+1)
      END IF
      if (NORTHERN_EDGE .and. EASTERN_EDGE) then
        cff1=ABS(ABS(umask_wet(Iend+1,Jend+1))-1.)
        cff2=0.5+SIGN(0.5,ubar(Iend+1,Jend+1,knew))
     &                                       *umask_wet(Iend+1,Jend+1)
        umask_wet(Iend+1,Jend+1)=0.5*umask_wet(Iend+1,Jend+1)*cff1
     &                                              +cff2*(1.-cff1)
        ubar(Iend+1,Jend+1,knew)=ubar(Iend+1,Jend+1,knew)
     &                                        *umask_wet(Iend+1,Jend+1)
      END IF
# endif
#endif  /* WET_DRY */

#ifdef NBQ
!=======================================================================
!  Computes normal boundary values of DU_nbq for NBQ OBCs
!=======================================================================
!
# ifndef EW_COM_PERIODIC
      if (WESTERN_EDGE) then
#  if defined OBC_COM_WEST || defined MRL_WCI
        do j=Jstr,Jend
          DU_nbq(Istr,j)=ubar(Istr,j,knew)
     &                     *0.5*(h(Istr  ,j)+zeta(Istr  ,j,knew)
     &                         + h(Istr-1,j)+zeta(Istr-1,j,knew))
#   ifdef NBQ_MASS
     &                           *0.5*(rhobar_nbq(Istr  ,j,knew)
     &                                +rhobar_nbq(Istr-1,j,knew))
#   endif
        enddo
#  else
        do j=Jstr,Jend
          DU_nbq(Istr,j)=0.  ! closed boundary
        enddo
#  endif
      endif
!
      if (EASTERN_EDGE) then
#  if defined OBC_COM_EAST || defined MRL_WCI
        do j=Jstr,Jend
          DU_nbq(Iend+1,j)=ubar(Iend+1,j,knew)
     &                     *0.5*(h(Iend+1,j)+zeta(Iend+1,j,knew)
     &                          +h(Iend  ,j)+zeta(Iend  ,j,knew))
#   ifdef NBQ_MASS
     &                           *0.5*(rhobar_nbq(Iend+1,j,knew)
     &                                +rhobar_nbq(Iend  ,j,knew))
#   endif
        enddo
#  else
        do j=Jstr,Jend
          DU_nbq(Iend+1,j)=0.  ! closed boundary
        enddo
#  endif
      endif
# endif /* EW_COM_PERIODIC */
!
# ifndef NS_COM_PERIODIC
      if (SOUTHERN_EDGE) then
#  if defined OBC_COM_SOUTH || defined MRL_WCI
        do i=IstrU,Iend
          DU_nbq(i,Jstr-1)=ubar(i,Jstr-1,knew)
     &                 *0.5*(h(i  ,Jstr-1)+zeta(i  ,Jstr-1,knew)
     &                      +h(i-1,Jstr-1)+zeta(i-1,Jstr-1,knew))
#   ifdef NBQ_MASS
     &                         *0.5*(rhobar_nbq(i  ,Jstr-1,knew)
     &                              +rhobar_nbq(i-1,Jstr-1,knew))
#   endif
        enddo
#  else
        do i=IstrU,Iend
          DU_nbq(i,Jstr-1)=gamma2*DU_nbq(i,Jstr)  ! slip/noslip boundary
        enddo
#  endif
      endif
!
      if (NORTHERN_EDGE) then
#  if defined OBC_COM_NORTH  || defined MRL_WCI
        do i=IstrU,Iend
          DU_nbq(i,Jend+1)=ubar(i,Jend+1,knew)
     &                 *0.5*(h(i  ,Jend+1)+zeta(i  ,Jend+1,knew)
     &                      +h(i-1,Jend+1)+zeta(i-1,Jend+1,knew))
#   ifdef NBQ_MASS
     &                         *0.5*(rhobar_nbq(i  ,Jend+1,knew)
     &                              +rhobar_nbq(i-1,Jend+1,knew))
#   endif
        enddo
#  else
        do i=IstrU,Iend
          DU_nbq(i,Jend+1)=gamma2*DU_nbq(i,Jend)
        enddo
#  endif
      endif
# endif /* NS_COM_PERIODIC */
!
      endif ! M2bc_nbq_flag
!
#endif /* NBQ */

      return
      end
#ifndef CHILD
# define CHILD
# ifdef AGRIF
#  include "u2dbc.F"
# endif
# undef CHILD
#endif  /* !CHILD */
