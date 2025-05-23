! $Id: wrt_bio_diags.F 1571 2014-07-01 12:38:05Z gcambon $
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
#ifdef DIAGNOSTICS_BIO
                                   ! Writes requested model
      subroutine wrt_bio_diags     ! fields at requested levels
                                   ! into diagnostics netCDF file.
      implicit none
      integer ierr,  record, lstr, lvar, lenstr, type
     &   , start(2), count(2), ibuff(4), nf_fwrite
# if defined MPI & !defined PARALLEL_FILES
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
# include "diagnostics.h"
# include "mpi_cpl.h"

# ifdef SOLVE3D
      integer tile, iflux,i,j,k
#  include "work.h"
# endif

# include "netcdf.inc"

# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &       1, MPI_COMM_WORLD, status, ierr)
      endif
# endif
!
!     Create/open diagnostic file; write grid arrays, if so needed.
!
      call def_bio_diags (nciddiabio, nrecdiabio, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(dianamebio)
!                                            !!! WARNING: Once time
! Set record within the file.                !!! stepping has been
                                             !!! started, it is assumed
      if (iic.eq.0) nrecdiabio=nrecdiabio+1  !!! that global history
      if (nrpfdiabio.eq.0) then              !!! record index is
        record=nrecdiabio                    !!! advanced by main.
      else
        record=1+mod(nrecdiabio-1, nrpfdiabio)
      endif
!
! Write out evolving model variables:
! ----- --- -------- ----- ----------
!
! Time step number and record numbers.
!
      type=filetype_diabio
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
!#ifdef AVERAGES
!      ibuff(4)=nrecavg
!#else
!      ibuff(4)=0
!#endif
      ibuff(4)=nrecdiabio
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiabio,diaTstepbio,start,count,ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
CDEBUG      write(*,*) 'Write Time'
      ierr=nf_put_var1_FTYPE (nciddiabio, diaTimebio, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
!     Time2
!
      ierr=nf_put_var1_FTYPE (nciddiabio, diaTime2bio, record, time)
      if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxTime2))
         write(stdout,1) vname(1,indxTime2)(1:lvar), record, ierr
     &        MYID
         goto 99                !--> ERROR
      endif

!
! Tracer diag variables.
!
!
! Flux terms

# ifdef PISCES

#  if defined key_trc_diaadd
! 3-D terms (it is called Flux but it means nothing) Christophe Menkes
       do iflux=1,NumFluxTerms
         if (wrtdiabioFlux(iflux)) then
           workr=bioFlux(:,:,:,iflux)
           call fillvalue3d(workr, nciddiabio, diabioFlux(iflux),
     &          indxbioFlux+iflux-1,record,r3dvar,type)
         endif
       end do
!     2-D terms (it is called Vsink but it means nothing) Christophe Menkes
       do iflux=1,NumVSinkTerms
         if (wrtdiabioVSink(iflux)) then
!     if (iflux.eq.11)  write(*,*) 'ZMEU WRT',bioVSink(30,30,iflux)
          work2d=bioVSink(:,:,iflux)
          call fillvalue2d(work2d, nciddiabio, diabioVsink(iflux),
     &         indxbioVSink+iflux-1,record,r2dvar,type)
        endif
      end do

#  endif

# else

!     Flux terms
      do iflux = 1, NumFluxTerms
        if (wrtdiabioFlux(iflux)) then
!     write(*,*)'write fluxes terms'
          workr=bioFlux(:,:,:,iflux)
          call fillvalue3d(workr, nciddiabio, diabioFlux(iflux),
     &         indxbioFlux+iflux-1,record,r3dvar,type)
        endif
      end do
!     vertical sinking fluxes
      do iflux = 1, NumVSinkTerms
        if (wrtdiabioVSink(iflux)) then
!     write(*,*)'write vertical sinking fluxes'
          work=bioVsink(:,:,:,iflux)
          call fillvalue3d(work, nciddiabio, diabioVsink(iflux),
     &         indxbioVSink+iflux-1,record,w3dvar,type)
        endif
      end do
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
!     gas exchange fluxes
      do iflux = 1, NumGasExcTerms
        if (wrtdiabioGasExc(iflux)) then
          work2d=GasExcFlux(:,:,iflux)
          call fillvalue2d(work2d, nciddiabio, diabioGasExc(iflux),
     &         indxGasExcFlux+iflux-1,record,r2dvar,type)
        endif
      end do
#  endif

# endif

 1    format(/1x,'WRT_BIO_DIAG ERROR while writing variable ''', A,
     &     ''' into diag file.', /11x, 'Time record:',
     &     I6,3x,'netCDF error code',i4,3x,a,i4)
      goto 100
 99   may_day_flag=3
 100  continue

!
!     Synchronize netCDF file to disk to allow other processes
!     to access data immediately after it is written.
!

# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      ierr=nf_close (nciddiabio)
      if (nrpfdiabio.gt.0 .and. record.ge.nrpfdiabio) nciddiabio=-1
# else
      if (nrpfdiabio.gt.0 .and. record.ge.nrpfdiabio) then
        ierr=nf_close (nciddiabio)
        nciddiabio=-1
      else
        ierr=nf_sync(nciddiabio)
      endif
# endif

      if (ierr .eq. nf_noerr) then
       MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &       'WRT_BIO_DIAG -- wrote ',
     &       'diag fields into time record =', record, '/',
     &       nrecdiabio  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)')
     &       'WRT_BIO_DIAG ERROR: Cannot ',
     &       'synchronize/close diag netCDF file.'
        may_day_flag=3
      endif

# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode .lt. NNODES-1) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &       1, MPI_COMM_WORLD,  ierr)
      endif
# endif
      return
      end

#else
      subroutine wrt_bio_diag_empty
      end
#endif /* DIAGNOSTICS_BIO*/

