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
/* Weigths matrix for leray smoothing */
#ifdef LERAY_FILTER_3PTS
      real filter_weights(-1:1,-1:1)
#elif defined LERAY_FILTER_5PTS
      real filter_weights(-2:2,-2:2)
#elif defined LERAY_FILTER_7PTS
      real filter_weights(-3:3,-3:3)
#elif defined LERAY_FILTER_9PTS
      real filter_weights(-4:4,-4:4)
#endif    
      common/fweights/filter_weights
      real weights_sum3, weights_sum5, weights_sum7, weights_sum9
      common/fweights_sum/weights_sum3, weights_sum5, 
     &                    weights_sum7, weights_sum9
!
/* Array for filter window size reduction near boundaries*/
      integer u_fwidth_array(GLOBAL_2D_ARRAY), 
     &        v_fwidth_array(GLOBAL_2D_ARRAY)
      common/fwidth_array/u_fwidth_array,v_fwidth_array
