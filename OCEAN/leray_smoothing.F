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
#ifdef LERAY_TURB
!      
!Computation of the filter weights to apply to velocity components in
!advection operators (Petersen, 2008)
      subroutine calc_smoothing_weights()
      implicit none    
# include "param.h"
# include "leray_smoothing.h"

# ifdef LERAY_FILTER_9PTS     
      !j-4
      filter_weights(-4,-4) = evs*evs
      filter_weights(-3,-4) = evs*dvs
      filter_weights(-2,-4) = evs*cvs
      filter_weights(-1,-4) = evs*bvs
      filter_weights( 0,-4) = evs
      filter_weights( 1,-4) = evs*bvs
      filter_weights( 2,-4) = evs*cvs
      filter_weights( 3,-4) = evs*dvs
      filter_weights( 4,-4) = evs*evs
      !j-3
      filter_weights(-4,-3) = dvs*evs
      filter_weights( 4,-3) = dvs*evs
      !j-2
      filter_weights(-4,-2) = cvs*evs
      filter_weights( 4,-2) = cvs*evs
      !j-1
      filter_weights(-4,-1) = bvs*evs
      filter_weights( 4,-1) = bvs*evs
      !j
      filter_weights(-4, 0) = evs
      filter_weights( 4, 0) = evs
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS
      !j-3
      filter_weights(-3,-3) = dvs*dvs
      filter_weights(-2,-3) = dvs*cvs
      filter_weights(-1,-3) = dvs*bvs
      filter_weights( 0,-3) = dvs
      filter_weights( 1,-3) = dvs*bvs
      filter_weights( 2,-3) = dvs*cvs
      filter_weights( 3,-3) = dvs*dvs
      !j-2
      filter_weights(-3,-2) = cvs*dvs
      filter_weights( 3,-2) = cvs*dvs
      !j-1
      filter_weights(-3,-1) = bvs*dvs
      filter_weights( 3,-1) = bvs*dvs
      !j
      filter_weights(-3, 0) = dvs
      filter_weights( 3, 0) = dvs
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS
      !j-2
      filter_weights(-2,-2) = cvs*cvs
      filter_weights(-1,-2) = cvs*bvs
      filter_weights( 0,-2) = cvs
      filter_weights( 1,-2) = cvs*bvs
      filter_weights( 2,-2) = cvs*cvs
      !j-1
      filter_weights(-2,-1) = bvs*cvs
      filter_weights( 2,-1) = bvs*cvs
      !j
      filter_weights(-2, 0) = cvs
      filter_weights( 2, 0) = cvs
# endif      
      !j-1
      filter_weights(-1,-1) = bvs*bvs
      filter_weights( 0,-1) = bvs
      filter_weights( 1,-1) = bvs*bvs
      !j
      filter_weights(-1, 0) = bvs
      filter_weights( 0, 0) = avs
      filter_weights( 1, 0) = bvs
      !j+1
      filter_weights(:, 1) = filter_weights(:, -1)
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS      
      !j+2
      filter_weights(:, 2) = filter_weights(:, -2)
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS      
      !j+3
      filter_weights(:, 3) = filter_weights(:, -3)
# endif
# if defined LERAY_FILTER_9PTS      
      !j+4
      filter_weights(:, 4) = filter_weights(:, -4)
# endif      

      weights_sum3 = sum( filter_weights(-1:1,-1:1) )
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS      
      weights_sum5 = sum( filter_weights(-2:2,-2:2) )
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS      
      weights_sum7 = sum( filter_weights(-3:3,-3:3) )
# endif
# if defined LERAY_FILTER_9PTS      
      weights_sum9 = sum( filter_weights(-4:4,-4:4) )
# endif      

      return
      end subroutine
!
!Count the number of available points around each point to determine the
!size of the filter window:
!        - no available points (domain boundary or mask points)
      subroutine fill_fwidth_array(tile)
      implicit none
      integer  tile, trd, omp_get_thread_num
# include "param.h"      
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call fill_fwidth_array_tile (istr,iend,jstr,jend)
      return
      end

      subroutine fill_fwidth_array_tile(istr,iend,jstr,jend)
      implicit none
# include "param.h"
# include "ocean3d.h"
# include "leray_smoothing.h" 

      integer istr,iend,jstr,jend,i,j,k
      integer imin,imax,jmin,jmax
# include "compute_auxiliary_bounds.h"
!compute imin,jmin,imax,jmax according to periodic boundaries, INTER or .not.INTER
# ifdef NS_PERIODIC
      jmin=Jstr-filter_width/2
      jmax=Jend+filter_width/2
# else
#  ifdef MPI
      if (SOUTH_INTER) then
        jmin=Jstr-filter_width/2
      else
        jmin=Jstr
      endif
      if (NORTH_INTER) then
        jmax=Jend+filter_width/2
      else
        jmax=Jend
      endif
#  else
      jmin=Jstr
      jmax=Jend
#  endif
# endif      
# ifdef EW_PERIODIC
      imin=IstrU-filter_width/2
      imax=Iend+filter_width/2
# else
#  ifdef MPI
      if (WEST_INTER) then
        imin=IstrU-filter_width/2
      else
        imin=IstrU
      endif
      if (EAST_INTER) then
        imax=Iend+filter_width/2
      else
        imax=Iend
      endif
#  else
      imin=IstrU
      imax=Iend
#  endif
# endif      
      u_fwidth_array(IstrU:Iend,Jstr:Jend)=filter_width
      if (WESTERN_EDGE  .or. EASTERN_EDGE .or.
     &    SOUTHERN_EDGE .or. NORTHERN_EDGE) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            do k=filter_width/2,1,-1
              if(i-k.ge.imin .and. i+k.le.imax
     &     .and. j-k.ge.jmin .and. j+k.le.jmax) then
                exit
              else
                u_fwidth_array(i,j)=2*(k-1)+1
              endif  
            enddo  
          enddo
        enddo
      endif  
            
# ifdef NS_PERIODIC
      jmin=JstrV-filter_width/2
      jmax=Jend+filter_width/2
# else
#  ifdef MPI
      if (SOUTH_INTER) then
        jmin=JstrV-filter_width/2
      else
        jmin=JstrV
      endif
      if (NORTH_INTER) then
        jmax=Jend+filter_width/2
      else
        jmax=Jend
      endif
#  else
      jmin=JstrV
      jmax=Jend
#  endif
# endif      
# ifdef EW_PERIODIC
      imin=Istr-filter_width/2
      imax=Iend+filter_width/2
# else
#  ifdef MPI
      if (WEST_INTER) then
        imin=Istr-filter_width/2
      else
        imin=Istr
      endif
      if (EAST_INTER) then
        imax=Iend+filter_width/2
      else
        imax=Iend
      endif
#  else
      imin=Istr
      imax=Iend
#  endif
# endif        
      v_fwidth_array(Istr:Iend,JstrV:Jend)=filter_width
      if (WESTERN_EDGE  .or. EASTERN_EDGE .or.
     &    SOUTHERN_EDGE .or. NORTHERN_EDGE) then
        do j=JstrV,Jend
          do i=Istr,Iend
            do k=filter_width/2,1,-1
              if(i-k.ge.imin .and. i+k.le.imax
     &     .and. j-k.ge.jmin .and. j+k.le.jmax) then
                exit
              else
                v_fwidth_array(i,j)=2*(k-1)+1
              endif  
            enddo  
          enddo
        enddo
      endif  

      return
      end subroutine


#else
      subroutine leray_smoothing_empty
      return
      end
#endif
