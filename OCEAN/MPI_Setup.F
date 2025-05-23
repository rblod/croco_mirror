! $Id: MPI_Setup.F 1458 2014-02-03 15:01:25Z gcambon $
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
#ifdef MPI

      subroutine MPI_Setup (ierr)

!     implicit none
# include "param.h"
      integer ierr, nsize, ii_W, ii_E, jj_S, jj_N
      integer :: lbx,ubx,lby,uby
      integer chunk_size_X,margin_X,chunk_size_E,margin_E
      integer Istrmpi,Iendmpi,Jstrmpi,Jendmpi, i_X,j_E
      integer,dimension(0:NP_XI-1,0:NP_ETA-1)  ::   myproc
!      real,dimension(LLm+2,MMm+2)  ::   zmask
      real,dimension(0:LLm+1,0:MMm+1)  ::   zmask
      integer,dimension(0:NNODES-1) :: mymynode!,myistr,myjstr,
!    &    myiend,myjend
      integer,dimension(0:NNODES-1) :: mymynode2
      integer ncid,varid,nerr
      integer ip1,im1,jp1,jm1,inu2
      integer Npts
# ifndef MP_3PTS
      parameter (Npts=2)
# else
      parameter (Npts=3)
# endif


# include "scalars.h"
# include "ncscrum.h"
# include "mpi_cpl.h"
# include "netcdf.inc"
      include 'mpif.h'
      call MPI_Comm_size (MPI_COMM_WORLD, nsize,  ierr)
      call MPI_Comm_rank (MPI_COMM_WORLD, mynode, ierr)
      if (nsize.ne.NNODES) then
        MPI_master_only write(stdout,'(/1x,A,I5,1x,A,I5,A/)')
     &     'ERROR in MPI_Setup: number of MPI-nodes should be',
     &                         NNODES, 'instead of', nsize, '.'
        ierr=99
        return
      endif

!       write(200+mynode,*) mynode,nsize,NNODES
!       stop'toto'
!       write(200+mynode,*) "passe la 1"
!       print *,"passe la 1"
!       stop '1'
# ifndef MPI_NOLAND
      ii=mod(mynode,NP_XI)
      jj=mynode/NP_XI

      if (NP_XI.eq.1) then
        WEST_INTER=.false.
        EAST_INTER=.false.
      else
#  ifdef EW_PERIODIC
        WEST_INTER=.true.
        EAST_INTER=.true.
#  else
        if (ii.eq.0) then
          WEST_INTER=.false.
        else
          WEST_INTER=.true.
        endif
        if (ii.eq.NP_XI-1) then
          EAST_INTER=.false.
        else
          EAST_INTER=.true.
        endif
#  endif
      endif

      if (NP_ETA.eq.1) then
        SOUTH_INTER=.false.
        NORTH_INTER=.false.
      else
#  ifdef NS_PERIODIC
        SOUTH_INTER=.true.
        NORTH_INTER=.true.
#  else
        if (jj.eq.0) then
          SOUTH_INTER=.false.
        else
          SOUTH_INTER=.true.
        endif
        if (jj.eq.NP_ETA-1) then
          NORTH_INTER=.false.
        else
          NORTH_INTER=.true.
        endif
#  endif
      endif

      ii_W=mod(ii-1+NP_XI,NP_XI)
      ii_E=mod(ii+1       ,NP_XI)
      jj_S=mod(jj-1+NP_ETA,NP_ETA)
      jj_N=mod(jj+1       ,NP_ETA)

      p_W=ii_W +NP_XI*jj
      p_E=ii_E +NP_XI*jj
      p_S=ii   +NP_XI*jj_S
      p_N=ii   +NP_XI*jj_N

      p_NW=ii_W+NP_XI*jj_N
      p_SW=ii_W+NP_XI*jj_S
      p_NE=ii_E+NP_XI*jj_N
      p_SE=ii_E+NP_XI*jj_S


      j_E=mynode/NP_XI
      i_X=mynode-j_E*NP_XI


      WEST_INTER2=WEST_INTER
      EAST_INTER2=EAST_INTER
      SOUTH_INTER2=SOUTH_INTER
      NORTH_INTER2=NORTH_INTER

      CORNER_SW=SOUTH_INTER .and.  WEST_INTER
      CORNER_NW=NORTH_INTER .and.  WEST_INTER
      CORNER_SE=SOUTH_INTER .and.  EAST_INTER
      CORNER_NE=NORTH_INTER .and.  EAST_INTER

      NNODES2=NP_XI*NP_ETA
      mynode2=mynode
# else
  !   ! read mask
      nerr= nf_open('croco_grd.nc',
     &   NF_NOWRITE, ncid)
      nerr= nf_inq_varid(ncid, 'mask_rho', varid)
      nerr= nf_get_var_double(ncid, varid, zmask)
      if(nerr /= nf_noerr ) then
        write(*,*) 'Reading mask file failed'
        call mpi_abort(MPI_COMM_WORLD, ierr)
      endif
      nerr= nf_close(ncid)

      WEST_INTER=.true.
      EAST_INTER=.true.
      SOUTH_INTER=.true.
      NORTH_INTER=.true.
      WEST_INTER2=.true.
      EAST_INTER2=.true.
      SOUTH_INTER2=.true.
      NORTH_INTER2=.true.
      CORNER_SW= .true.
      CORNER_NW= .true.
      CORNER_SE=  .true.
      CORNER_NE= .true.

      inu=0 ; inu2=0
      myproc(:,:)=mpi_proc_null
      mymynode(:)=0
      mymynode2(:)=0
      ! loop on all sub-domain
      do j=0,NP_ETA-1
        do i=0,NP_XI-1
          j_E=inu2/NP_XI
          i_X=inu2-j_E*NP_XI
          !
          chunk_size_X=(LLm+NP_XI-1)/NP_XI
          margin_X=(NP_XI*chunk_size_X-LLm)/2
          chunk_size_E=(MMm+NP_ETA-1)/NP_ETA
          margin_E=(NP_ETA*chunk_size_E-MMm)/2
          !
          istrmpi=1+i_X*chunk_size_X-margin_X !-Npts
          iendmpi=istrmpi+chunk_size_X-1  ! +Npts
          istrmpi=max(istrmpi,1)
          iendmpi=min(iendmpi,LLm)
          !
          jstrmpi=1+j_E*chunk_size_E-margin_E !-NPTS
          jendmpi=jstrmpi+chunk_size_E-1 !  +Npts
          jstrmpi=max(jstrmpi,1)
          jendmpi=min(jendmpi,Mmm)

          lbx=max(istrmpi-Npts,1)
          ubx=min(iendmpi+Npts,LLm)
          lby=max(jstrmpi-Npts,1)
          uby=min(jendmpi+Npts,Mmm)

          if(sum(zmask(lbx:ubx,lby:uby))>0.)then
            myproc(i,j)=inu
            mymynode(inu)=inu
            mymynode2(inu)=inu2
            inu=inu+1
          endif
          inu2=inu2+1
        enddo
      enddo

      mynode2=mymynode2(mynode)
      mynode=mymynode(mynode)
      NNODES2=NP_XI*NP_ETA

      !
      do j=0,NP_ETA -1
        do i=0,NP_XI  -1
          if (myproc(i,j)==mynode)   then
            im1=mod(i-1+NP_XI,NP_XI)
            jm1=mod(j-1+NP_ETA,NP_ETA)
            ip1=mod(i+1,NP_XI)
            jp1=mod(j+1,NP_ETA)

            i_X=i
            ii=i_X
            j_E=j
            jj=j_e

            if(i==0 ) WEST_INTER2=.false.
            if(myproc(im1,j)==mpi_proc_null )WEST_INTER2=.false.
            p_W=myproc(im1,j)

            if(i==NP_XI-1) EAST_INTER2=.false.
            if(myproc(ip1,j)==mpi_proc_null ) EAST_INTER2=.false.
            p_E=myproc(ip1,j)

            if(j==0.) SOUTH_INTER2=.false.
            if(myproc(i,jm1) ==mpi_proc_null ) SOUTH_INTER2=.false.
            p_S=myproc(i,jm1)

            if(j==NP_ETA-1) NORTH_INTER2=.false.
            if(myproc(i,jp1)==mpi_proc_null ) NORTH_INTER2=.false.
            p_N=myproc(i,jp1)

            p_SW=myproc(im1,jm1)
            p_NW=myproc(im1,jp1)
            p_SE=myproc(ip1,jm1)
            p_NE=myproc(ip1,jp1)

            if(p_SW== mpi_proc_null) CORNER_SW=.false.
            if(p_NW== mpi_proc_null) CORNER_NW=.false.
            if(p_SE== mpi_proc_null) CORNER_SE=.false.
            if(p_NE== mpi_proc_null) CORNER_NE=.false.

          endif
        enddo
      enddo

      if (NP_XI.eq.1) then
        WEST_INTER=.false.
        EAST_INTER=.false.
      else
#  ifdef EW_PERIODIC
        WEST_INTER=.true.
        EAST_INTER=.true.
#  else
        if (ii.eq.0) then
          WEST_INTER=.false.
        else
          WEST_INTER=.true.
        endif
        if (ii.eq.NP_XI-1) then
          EAST_INTER=.false.
        else
          EAST_INTER=.true.
        endif
#  endif
      endif

      if (NP_ETA.eq.1) then
        SOUTH_INTER=.false.
        NORTH_INTER=.false.
      else
#  ifdef NS_PERIODIC
        SOUTH_INTER=.true.
        NORTH_INTER=.true.
#  else
        if (jj.eq.0) then
          SOUTH_INTER=.false.
        else
          SOUTH_INTER=.true.
        endif
        if (jj.eq.NP_ETA-1) then
          NORTH_INTER=.false.
        else
          NORTH_INTER=.true.
        endif
#  endif
      endif

# endif /* MPI_NOLAND */

       IF( .not. north_inter .or. .not.east_inter) CORNER_NE=.false.
       IF( .not. north_inter .or. .not.west_inter) CORNER_NW=.false.
       IF( .not. south_inter .or. .not.east_inter) CORNER_SE=.false.
       IF( .not. south_inter .or. .not.west_inter) CORNER_SW=.false.

      chunk_size_X=(LLm+NP_XI-1)/NP_XI
      margin_X=(NP_XI*chunk_size_X-LLm)/2
      chunk_size_E=(MMm+NP_ETA-1)/NP_ETA
      margin_E=(NP_ETA*chunk_size_E-MMm)/2

      istrmpi=1+i_X*chunk_size_X-margin_X
      iendmpi=istrmpi+chunk_size_X-1
      istrmpi=max(istrmpi,1)
      iendmpi=min(iendmpi,LLm)

      jstrmpi=1+j_E*chunk_size_E-margin_E
      jendmpi=jstrmpi+chunk_size_E-1
      jstrmpi=max(jstrmpi,1)
      jendmpi=min(jendmpi,MMm)

      Lmmpi=iendmpi-istrmpi+1
      Mmmpi=jendmpi-jstrmpi+1

      iminmpi=istrmpi
      imaxmpi=iendmpi

      jminmpi=jstrmpi
      jmaxmpi=jendmpi

     ! write(*,*) mynode,istrmpi,iendmpi,jstrmpi,jendmpi


!       if(mynode == 0) then
!       open (unit=100,file="results.txt",action="write",
!     &  status="replace")!

 !       write(100,*) 'DECOUPAGE : '
 !       write(100,*)
 !       do j=NP_ETA-1,0,-1
 !          write(100,*) (myproc(i,j),i=0,NP_XI-1)
 !       enddo
 !      endif
 !      call MPI_Barrier (MPI_COMM_WORLD, ierr)
 !      if(mynode == 0) close(100)

 !       if(jj==4)then
 !      write(*,*) 'diags : ', mynode, CORNER_NE,CORNER_NW,
 !    & CORNER_SE,CORNER_SW
 !      write(*,*) 'bounds : ',mynode,istrmpi,iendmpi,jstrmpi,jendmpi

 !      endif
 !      call MPI_Barrier (MPI_COMM_WORLD, ierr)

# ifdef PARALLEL_FILES
!# ifndef EW_PERIODIC
      xi_rho=Lmmpi
      xi_u=xi_rho
      if (ii.eq.0) xi_rho=xi_rho+1
      if (ii.eq.NP_XI-1) then
        xi_rho=xi_rho+1
        xi_u=xi_u+1
      endif
!# endif
!# ifndef NS_PERIODIC
      eta_rho=Mmmpi
      eta_v=eta_rho
      if (jj.eq.0) eta_rho=eta_rho+1
      if (jj.eq.NP_ETA-1) then
        eta_rho=eta_rho+1
        eta_v=eta_v+1
      endif
!# endif
# endif

      return
      end


# undef CHECK_MPI
# ifdef CHECK_MPI

      subroutine MPI_Test
      USE param
      implicit none
c#  include "param.h"
      integer tile
      do tile=0,NSUB_X*NSUB_E-1
        call MPI_Test1 (tile)
      enddo
      return
      end

      subroutine MPI_Test1 (tile)
      USE param
      implicit none
      integer tile
c#  include "param.h"
#  include "compute_tile_bounds.h"
      call MPI_Test1_tile  (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine MPI_Test1_tile (Istr,Iend,Jstr,Jend)
      implicit none
      integer Istr,Iend,Jstr,Jend
#  include "param.h"
#  include "scalars.h"
#  include "mpi_cpl.h"
      include 'mpif.h'

      integer i,j,k,ierr
      real temp2D(GLOBAL_2D_ARRAY),
     &     temp3D(GLOBAL_2D_ARRAY,0:N)
      common /MPI_Test_Arr/ temp2D,temp3D
      character string(128)
!
#  include "compute_extended_bounds.h"
!
      do j=JstrR,JendR
        do i=IstrR,IendR
c          temp2D(i,j)=1.*mynode+0.5
           temp2D(i,j)=-1.*(mynode+1)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          temp2D(i,j)=1.*mynode+1.
c          temp2D(i,j)=i+ii*Lm
c          temp2D(i,j)=j+jj*Mm
        enddo
      enddo


      call MessPass2D_tile (Istr,Iend,Jstr,Jend, temp2D)

      do i=mynode,NNODES
        call MPI_Barrier (MPI_COMM_WORLD, ierr)
      enddo
      write(*,*)
      write(*,*)
      write(*,'(A5,I3,2I4)') 'node=',mynode, ii,jj
      write(*,*)
      do j=JendR,JstrR,-1
       write(string,'(I3,3x,20F4.1)') j, (temp2D(i,j), i=IstrR,IendR)
        do i=7,6+4*(IendR-IstrR+1)
         if (string(i-1).eq.'.' .and. string(i).eq.'0') string(i)=' '
        enddo
       write(*,'(128A1)') (string(i), i=1,6+4*(IendR-IstrR+1))
      enddo
      write(*,*)
        write(*,'(5x,20I4)') (i,i=IstrR,IendR)
      write(*,*)
      write(*,*)
      do i=0,mynode
        call MPI_Barrier (MPI_COMM_WORLD, ierr)
      enddo


c      return


      do k=0,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            temp3D(i,j,k)=-1.*(mynode+1)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            temp3D(i,j,k)=1.*mynode+1.
             temp3D(i,j,k)=0.1*( float(j-1+(Jend-Jstr+1)*mynode) )
          enddo
        enddo
      enddo

      call MessPass3D_tile (Istr,Iend,Jstr,Jend, temp3D,N+1)

      do k=0,N
        write(*,*) 'k=',k
        do i=mynode,NNODES
          call MPI_Barrier (MPI_COMM_WORLD, ierr)
        enddo
        write(*,*)
        write(*,*)
        write(*,'(A5,I3,2I4)') 'node=',mynode, ii,jj
        write(*,*)
        do j=JendR,JstrR,-1
        write(string,'(I3,3x,20F4.1)')j,(temp3D(i,j,k),i=IstrR,IendR)
         do i=7,6+4*(IendR-IstrR+1)
          if (string(i-1).eq.'.'.and.string(i).eq.'0') string(i)=' '
         enddo
        write(*,'(128A1)') (string(i), i=1,6+4*(IendR-IstrR+1))
        enddo
        write(*,*)
        write(*,'(5x,20I4)') (i,i=IstrR,IendR)
        write(*,*)
        write(*,*)
        do i=0,mynode
          call MPI_Barrier (MPI_COMM_WORLD, ierr)
        enddo
      enddo
      return
      end
# endif        /* CHECK_MPI */
#else
      subroutine MPI_Setup_empty
      end
#endif    /* MPI */

