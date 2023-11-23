#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#    if defined LERAY_FILTER_9PTS
      call exchange_u2d_4pts_tile (Istr,Iend,Jstr,Jend,
     &                    Duon(START_2D_ARRAY))
      call exchange_v2d_4pts_tile (Istr,Iend,Jstr,Jend,
     &                    Dvom(START_2D_ARRAY))
#    elif defined THREE_GHOST_POINTS_UV | defined LERAY_FILTER_7PTS
      call exchange_u2d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                    Duon(START_2D_ARRAY))
      call exchange_v2d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                    Dvom(START_2D_ARRAY))
#    else
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   Duon(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   Dvom(START_2D_ARRAY))
#    endif      
#   endif      
      do j=Jstr,Jend
        do i=IstrU,Iend
#   ifdef LERAY_FILTER_9PTS
          if(u_fwidth_array(i,j).eq.9) then
            wrk3(i,j) = sum( DUon(i-4:i+4,j-4:j+4)
     &                     * filter_weights(-4:4,-4:4) )
     &                  / weights_sum9
          endif
#   endif          
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS
          if(u_fwidth_array(i,j).eq.7) then 
            wrk3(i,j) = sum( DUon(i-3:i+3,j-3:j+3)
     &                     * filter_weights(-3:3,-3:3) )
     &                  / weights_sum7
          endif  
#   endif
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS
          if(u_fwidth_array(i,j).eq.5) then 
            wrk3(i,j) = sum( DUon(i-2:i+2,j-2:j+2)
     &                     * filter_weights(-2:2,-2:2) )
     &                  / weights_sum5
          endif  
#   endif
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS | defined LERAY_FILTER_3PTS
          if(u_fwidth_array(i,j).eq.3) then 
            wrk3(i,j) = sum( DUon(i-1:i+1,j-1:j+1)
     &                     * filter_weights(-1:1,-1:1) )
     &                  / weights_sum3
          endif
#   endif
          if(u_fwidth_array(i,j).lt.3) then
              wrk3(i,j) = DUon(i,j)
          endif
        enddo
      enddo
      DUon(IstrU:Iend,Jstr:Jend) = wrk3(IstrU:Iend,Jstr:Jend)
!
      do j=JstrV,Jend
        do i=Istr,Iend
#   ifdef LERAY_FILTER_9PTS
          if(v_fwidth_array(i,j).eq.9) then
            wrk3(i,j) = sum( DVom(i-4:i+4,j-4:j+4)
     &                   * filter_weights(-4:4,-4:4) )
     &                  / weights_sum9
          endif
#   endif          
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS
          if(v_fwidth_array(i,j).eq.7) then 
            wrk3(i,j) = sum( DVom(i-3:i+3,j-3:j+3)
     &                     * filter_weights(-3:3,-3:3) )
     &                  / weights_sum7
          endif  
#   endif
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS
          if(v_fwidth_array(i,j).eq.5) then 
            wrk3(i,j) = sum( DVom(i-2:i+2,j-2:j+2)
     &                     * filter_weights(-2:2,-2:2) )
     &                  / weights_sum5
          endif  
#   endif
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS | defined LERAY_FILTER_3PTS
          if(v_fwidth_array(i,j).eq.3) then 
            wrk3(i,j) = sum( DVom(i-1:i+1,j-1:j+1)
     &                       * filter_weights(-1:1,-1:1) )
     &                  / weights_sum3
          endif
#   endif
          if(v_fwidth_array(i,j).lt.3) then
              wrk3(i,j) = DVom(i,j)
          endif
        enddo
      enddo
      DVom(Istr:Iend,JstrV:Jend) = wrk3(Istr:Iend,JstrV:Jend)
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   Duon(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   Dvom(START_2D_ARRAY))
#   endif 
