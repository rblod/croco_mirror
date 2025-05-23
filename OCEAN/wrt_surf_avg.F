!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
!
! ROMS_AGRIF website : http://www.romsagrif.org
!======================================================================
!
#include "cppdefs.h"
#if defined OUTPUTS_SURFACE && defined AVERAGES && ! defined XIOS

!
!---------------------------------------------------------------
! Write diagnostics fields at requested levels into diagnostics
! netCDF file.

!---------------------------------------------------------------
!
      subroutine wrt_surf_avg

      implicit none
      integer ierr,  record, lstr, lvar, lenstr
     &   , start(2), count(2), ibuff(4), nf_fwrite, type
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE), blank
# endif
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "forces.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "mixing.h"
# include "surf.h"
# include "mpi_cpl.h"
# ifdef SOLVE3D
      integer tile,itrc,i,j,k,ivar
#  include "work.h"
# endif
# include "netcdf.inc"
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, MPI_COMM_WORLD, status, ierr)
      endif
# endif
!
# undef DEBUG
!
!
! Create/open diagnostic file; write grid arrays, if so needed.
!
      call def_surf_avg (ncidsurf_avg, nrecsurf_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(surfname_avg)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been
!
                                           !!! started, it is assumed
      if (iic.eq.0) nrecsurf_avg=nrecsurf_avg+1  !!! that global history
      if (nrpfsurf_avg.eq.0) then            !!! record index "nrecsurf"
        record=nrecsurf_avg                  !!! is advanced by main.
      else
        record=1+mod(nrecsurf_avg-1, nrpfsurf_avg)
      endif


!
!---------------------------------------------------------------
! Write out evolving model variables:
!---------------------------------------------------------------
!
! Time step number and record numbers.
!
      type=filetype_surf_avg
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecsurf_avg
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1

      ierr=nf_put_vara_int (ncidsurf_avg, surfTstep_avg,
     &                       start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step_avg', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (ncidsurf_avg, surfTime_avg,
     &                       record, timesurf_avg
# ifdef USE_CALENDAR
     &  - origin_date_in_sec
# endif
     &  )
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time2
!
      ierr=nf_put_var1_FTYPE (ncidsurf_avg, surfTime2_avg,
     &                       record, time
# ifdef USE_CALENDAR
     &  - origin_date_in_sec
# endif
     &  )
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime2))
        write(stdout,1) vname(1,indxTime2)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

!---------------------------------------------------------------
! Momentum diagnostic variables.
!---------------------------------------------------------------
!
!
          itrc = 1
          if (wrtsurf_avg(itrc)) then



!
! indxSST
!
           work2d=surft_avg(:,:)
           call fillvalue2d(work2d,ncidsurf_avg,
     &                      surf_surft_avg(itrc),
     &                      indxsurft+itrc-1,
     &                      record,r2dvar,type)
# ifdef SALINITY
!
! indxSSS
!
           work2d=surfs_avg(:,:)
           call fillvalue2d(work2d,ncidsurf_avg,
     &                      surf_surfs_avg(itrc),
     &                      indxsurfs+itrc-1,
     &                      record,r2dvar,type)
# endif
!
! indxSSH
!
           work2d=surfz_avg(:,:)
           call fillvalue2d(work2d,ncidsurf_avg,
     &                      surf_surfz_avg(itrc),
     &                      indxsurfz+itrc-1,
     &                      record,r2dvar,type)
!
! indxU
!
           work2d=surfu_avg(:,:)
           call fillvalue2d(work2d,ncidsurf_avg,
     &                      surf_surfu_avg(itrc),
     &                      indxsurfu+itrc-1,
     &                      record,u2dvar,type)
!
! indxV
!
           work2d=surfv_avg(:,:)
           call fillvalue2d(work2d,ncidsurf_avg,
     &                      surf_surfv_avg(itrc),
     &                      indxsurfv+itrc-1,
     &                      record,v2dvar,type)
!

        endif
!


  1   format(/1x,'WRT_surf_AVG ERROR while writing variable ''', A,
     &          ''' into diag file.', /11x, 'Time record:',
     &                   I6,3x,'netCDF error code',i4,3x,a,i4)
      goto 100
  99  may_day_flag=3
 100  continue

!
!---------------------------------------------------------------
! Synchronize netCDF file to disk to allow other processes
! to access data immediately after it is written.
!---------------------------------------------------------------
!
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      ierr=nf_close (ncidsurf_avg)
      if (nrpfsurf_avg.gt.0 .and. record.ge.nrpfsurf_avg)
     &            ncidsurf_avg=-1
# else
      if (nrpfsurf_avg.gt.0 .and. record.ge.nrpfsurf_avg) then
        ierr=nf_close (ncidsurf_avg)
        ncidsurf_avg=-1
      else
        ierr=nf_sync(ncidsurf_avg)
      endif
# endif


      if (ierr .eq. nf_noerr) then
        MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_surf_AVG -- wrote',
     &            ' diag fields into time record =', record, '/',
     &             nrecsurf_avg  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)')
     &             'WRT_surf_AVG ERROR: Cannot ',
     &             'synchronize/close diag netCDF file.'
        may_day_flag=3
      endif

# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode .lt. NNODES-1) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, MPI_COMM_WORLD,  ierr)
      endif
# endif
      return
      end

#else
      subroutine wrt_surf_avg_empty
      end
#endif /* (DIAGNOSTICS_EK) */
