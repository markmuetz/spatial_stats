      PROGRAM READPP
!
! Example fortran code to read in pp data
!
      IMPLICIT NONE
!
      CHARACTER(80):: FILE
! FILE - filename
!
      INTEGER, PARAMETER:: MAXSIZE=100,MAXNCL=1000,nbin=20,ntrials=100,NUMLOC=500
! MAXSIZE - max number of grid points in a single cloud
! MAXNCL - maximum number of clouds
! nbin - number of bins in the CDF
! ntrials - number of Monte-carlo trials
! NUMLOC - number of locations for to-nearest-cloud distances

      INTEGER:: ERROR,NX,NY,n,m,p,ncl
! ERROR - used to detect errors
! NX - number of x grid points
! NY - number of y grid points
! n,m,p - loop variables
! ncl - number of cloud structures found
      INTEGER:: k,ksq
! k - number of locations on regular grid
! ksq - total number of locations considered
      INTEGER, DIMENSION(45):: IHEAD
! IHEAD - integer headers to pp file
      INTEGER, DIMENSION(MAXNCL):: sizes,nnid
! sizes - number of grid points constituting the cloud
! nnid - identifier for the nearest neighbour
      INTEGER, DIMENSION(MAXNCL,ntrials):: mcnnid
! mcnnid - identifier for the nearest neighbour for the Monte-Carlo trial given
!          by the second index
      INTEGER, DIMENSION(:,:), ALLOCATABLE:: BINDAT
! BINDAT - binary data 
      INTEGER, DIMENSION(MAXNCL,MAXSIZE):: loci,locj
! loci - x index of grid points constituting cloud
! locj - y index of grid points constituting cloud
!
      LOGICAL:: END1
! END - signals end of file 
      LOGICAL:: check
! check - flag to signal if column needs checking

      REAL, PARAMETER:: RLGE=1.0E30
! RLGE - large real number
      REAL:: PI
! pi - the 3.14 thing
      REAL, DIMENSION(19):: RHEAD
! RHEAD - real headers to pp file
      REAL, DIMENSION(:,:), ALLOCATABLE:: INDAT
! INDAT - data read in from FILE
      REAL, DIMENSION(MAXNCL):: locx,locy,nnsep
! locx - x location of cloud centre
! locy - y location of cloud centre
! nnsep - nearest neighbour separation
      REAL, DIMENSION(MAXNCL,ntrials):: mcnnsep
! mcnnsep - nearest neighbour separation for the Monte-Carlo trial given by the
!           second index
      REAL, DIMENSION(nbin):: valbin,cdf_nnsep,cdf_mcnnsep,cdfnnmax,cdfnnmin, &
     &                        G0,cdf_sep,cdf_mcsep,cdfmax,cdfmin,cdf_locsep, &
     &                        cdf_mclocsep,cdflocmax,cdflocmin
! valbin - upper range of distance bin
! cdf_nnsep - CDF for the nearest neighbour separations
! cdf_mcnnsep - CDF of nearest neighbour separations for Monte-Carlo trials
! cdfnnmax - max of nearest neighbour CDF over the Monte-Carlo trials
! cdfnnmin - min of nearest neighbour CDF over the Monte-Carlo trials
! cdf_sep - CDF for the separations
! cdf_mcsep - CDF of separations for Monte-Carlo trials
! cdfmax - max of CDF over the Monte-Carlo trials
! cdfmin - min of CDF over the Monte-Carlo trials
! G0 - reference CDF
! cdf_locsep - CDF for the point-to-cloud separations
! cdf_mclocsep - CDF of point-to-cloud separations over the Monte-Carlo trials
! cdflocmax - max of point-to-cloud CDF over the Monte-Carlo trials
! cdflocmin - min of point-to-cloud CDF over the Monte-Carlo trials
      REAL, DIMENSION(NUMLOC):: locsep
! locsep - distance of locations from nearest cloud
      REAL, DIMENSION(MAXNCL,MAXNCL):: sep
! sep - separation of the first cloud to its neighbours in the second index
      REAL, DIMENSION(NUMLOC,ntrials):: mclocsep
! mclocsep - distance of locations from nearest cloud, for the Monte Carlo 
!            trials in the 2nd index
      REAL,  DIMENSION(MAXNCL,MAXNCL,ntrials):: mcsep
! mcsep - separation of the first cloud to its neighbours in the second index
!         for the Monte-Carlo trial in the 3rd index
!
! Open files 
!  
      FILE='rates_fixed.pp' ! filename
      OPEN(UNIT=81,FILE=FILE,FORM='UNFORMATTED',STATUS='OLD', &
      & ACTION='READ',IOSTAT=ERROR)
      IF (ERROR.NE.0) GOTO 9000      
!
! Allocate arrays
!
      READ(81,ERR=9006)IHEAD,RHEAD
      NX=IHEAD(19)
      NY=IHEAD(18)
      ALLOCATE(INDAT(NX,NY),STAT=ERROR)
      IF (ERROR.NE.0) GOTO 9001
      REWIND(81)
!     
      END1=.FALSE.
!
! We are now going to read in data slices, one by one until the file ends. If
! there are multiple data slices then the INDAT array is going to get 
! overwritten each time some new data is found
!
      check=.true.
 100  CONTINUE ! Start of loop for differencing 
      IF (END1) THEN ! Get data from file  if there is any
         INDAT=0.0
      ELSE
         READ(81,ERR=9006,END=101)IHEAD,RHEAD
         READ(81,ERR=9006,END=101)INDAT
      END IF
      GOTO 111 ! If we are here then there may be more data to come
 101  CONTINUE ! If we are here then this is the last of the data
      END1=.TRUE.
 111  CONTINUE
      IF (.NOT.(END1)) THEN ! If we have some data, then use it 
!
! The data is now loaded up for a horizontal NX,NY section in the array INDAT
! You can output this to screen or file
!
! You should insert code to do the output or use it somehow
!
      ALLOCATE(BINDAT(NX,NY),STAT=ERROR)
      IF (ERROR.NE.0) GOTO 9001
!
! Now threshold the data into a binary array
!
      CALL IDENTIFY_CLOUDY_POINTS(NX,NY,INDAT,BINDAT)
!
! Just use the second lot of data!
         IF (check) THEN
            DEALLOCATE(BINDAT)
            check=.false.
            GOTO 100! this starts the loop again to look for more data in the file
         END IF
      END IF
!
      CLOSE(81)
!
! Now process the last set of binary data in the file.
! First collate the objects together and assign a central position for each
!
      CALL CONSTRUCT_CLOUDS(ncl,MAXNCL,loci,locj,MAXSIZE,sizes,NX,NY,BINDAT)
! 
! Now determine the central point location for each cloud
!
      CALL CENTRE_OF_CLOUDS(ncl,MAXNCL,locx,locy,sizes,loci,locj,MAXSIZE,NX,NY)
!
! Now calculate cloud separations
!
      CALL CLOUDSEP(ncl,MAXNCL,locx,locy,sep,RLGE,nnid,nnsep,NUMLOC,LOCSEP)
!
! Compute Monte-Carlo samples of completely random spatial distributions
!
      DO n=1,ntrials
         CALL RANDOM_CLOUDS(ncl,MAXNCL,locx,locy,n)
         CALL CLOUDSEP(ncl,MAXNCL,locx,locy,mcsep(:,:,n),RLGE,mcnnid(:,n), &
     &                 mcnnsep(:,n),NUMLOC,mclocsep(:,n))
      END DO
!
! Now compute CDF of nnsep
!
      do n=1,nbin
         valbin(n)=MAXVAL(nnsep(1:ncl))*n/nbin
         cdf_nnsep(n)=COUNT(nnsep(1:ncl).LT.valbin(n))/(ncl*1.0)
!         write(*,*) valbin(n),cdf_nnsep(n)
         cdf_mcnnsep(n)=0.0
         cdfnnmax(n)=0.0
         cdfnnmin(n)=RLGE
         do m=1,ntrials
            cdf_mcnnsep(n)=cdf_mcnnsep(n)+COUNT(mcnnsep(1:ncl,m).LT.valbin(n))
            cdfnnmax(n)=MAX(cdfnnmax(n),COUNT(mcnnsep(1:ncl,m).LT.valbin(n)))
            cdfnnmin(n)=MIN(cdfnnmin(n),COUNT(mcnnsep(1:ncl,m).LT.valbin(n)))
         end do
         cdf_mcnnsep(n)=cdf_mcnnsep(n)/(ncl*ntrials*1.0)
         cdfnnmax(n)=cdfnnmax(n)/(ncl*1.0)
         cdfnnmin(n)=cdfnnmin(n)/(ncl*1.0)
      end do      
!
! Evaluate G0
!
      PI=4.0*ATAN(1.0)
      do n=1,nbin
         G0(n)=1-EXP(-ncl*pi*(valbin(n)**2))
!         write(*,*) n,valbin(n),cdf_nnsep(n),cdf_mcnnsep(n),G0(n),cdfnnmax(n),cdfnnmin(n)
      end do
!
! Look at CDF of all inter-cloud separations
!
      do n=1,nbin
         valbin(n)=MAXVAL(sep(1:ncl,1:ncl))*n/nbin
         do m=1,ncl
            sep(m,m)=RLGE
            do p=1,ntrials
               mcsep(m,m,p)=RLGE
            end do
         end do
! Note that the separations appear twice so we normalize by twice n(n-1)/2
         cdf_sep(n)=COUNT(sep(1:ncl,1:ncl).LT.valbin(n))/(1.0*ncl*(ncl-1))
         cdf_mcsep(n)=0.0
         cdfmax(n)=0.0
         cdfmin(n)=RLGE
         do m=1,ntrials
            cdf_mcsep(n)=cdf_mcsep(n)+COUNT(mcsep(1:ncl,1:ncl,m).LT.valbin(n))
            cdfmax(n)=MAX(cdfmax(n),COUNT(mcsep(1:ncl,1:ncl,m).LT.valbin(n)))
            cdfmin(n)=MIN(cdfmin(n),COUNT(mcsep(1:ncl,1:ncl,m).LT.valbin(n)))
         end do
         cdf_mcsep(n)=cdf_mcsep(n)/(ncl*(ncl-1)*ntrials*1.0)
         cdfmax(n)=cdfmax(n)/(ncl*(ncl-1)*1.0)
         cdfmin(n)=cdfmin(n)/(ncl*(ncl-1)*1.0)
!         write(*,*) n,valbin(n),cdf_sep(n),cdf_mcsep(n),cdfmax(n),cdfmin(n)
         do m=1,ncl
            sep(m,m)=0.0
            do p=1,ntrials
               mcsep(m,m,p)=0.0
            end do
         end do
      end do
!
! Look at CDF of point-to-nearest-cloud separations
!
      k=SQRT(NUMLOC*1.0)
      ksq=k**2
      do n=1,nbin
         valbin(n)=MAXVAL(locsep(1:ksq))*n/nbin    
         cdf_locsep(n)=COUNT(locsep(1:ksq).LT.valbin(n))/(ksq*1.0)
!         write(*,*) valbin(n),cdf_locsep(n)
         cdf_mclocsep(n)=0.0
         cdflocmax(n)=0.0
         cdflocmin=RLGE
         do m=1,ntrials
            cdf_mclocsep(n)=cdf_mclocsep(n)+COUNT(mclocsep(1:ksq,m).LT.valbin(n))
            cdflocmax(n)=MAX(cdflocmax(n),COUNT(mclocsep(1:ksq,m).LT.valbin(n)))
            cdflocmin(n)=MIN(cdflocmin(n),COUNT(mclocsep(1:ksq,m).LT.valbin(n)))
         end do
         cdf_mclocsep(n)=cdf_mclocsep(n)/(ksq*ntrials*1.0)
         cdflocmax(n)=cdflocmax(n)/(ksq*1.0)
         cdflocmin(n)=cdflocmin(n)/(ksq*1.0)
         write(*,*) n,valbin(n),cdf_locsep(n),cdf_mclocsep(n),cdflocmax(n),cdflocmin(n)
      end do 
!
      STOP
!
 9000 CONTINUE
      WRITE(*,*) 'Could not open input file',FILE
      GOTO 10000
 9001 CONTINUE
      WRITE(*,*) 'Allocation error'
      GOTO 10000
 9006 CONTINUE
      WRITE(*,*) 'Error reading from input file'
      GOTO 10000
10000 CONTINUE
      WRITE(*,*) 'EXECUTION TERMINATED'
      STOP
!
      END PROGRAM
!
      SUBROUTINE RANDOM_CLOUDS(ncl,MAXNCL,locx,locy,SEED)
!
      INTEGER, INTENT(IN):: ncl,MAXNCL
! ncl - number of cloud structures 
! MAXNCL - maximum number of clouds
      INTEGER, INTENT(IN), DIMENSION(1):: SEED
! SEED - seed for random number generator

      REAL, INTENT(OUT), DIMENSION(MAXNCL):: locx,locy
! locx - x location of cloud centre, on a normalized grid
! locy - y location of cloud centre, on a normalized grid

      INTEGER:: n
! n - loop variable

      CALL RANDOM_SEED(PUT=SEED)
      DO n=1,ncl
         CALL RANDOM_NUMBER(locx(n))
         CALL RANDOM_NUMBER(locy(n))
      END DO
!
      RETURN
!
      END SUBROUTINE RANDOM_CLOUDS
!
      SUBROUTINE CONSTRUCT_CLOUDS(ncl,MAXNCL,loci,locj,MAXSIZE,sizes, &
     &                            NX,NY,BINDAT)
!
      INTEGER, PARAMETER:: ndir=8,minsize=1
! ndir - number of neighbouring points to consider
! minsize - minimum number of gridpoints needed to define a viable cloud
      INTEGER, DIMENSION(ndir):: ioff,joff
! ioff,joff - used to identify neighbouring points
      data ioff/0,0,1,-1,-1,-1,1,1/
      data joff/1,-1,0,0,-1,1,-1,1/
      LOGICAL, PARAMETER:: cyclicbc=.false.
! cyclicbc - true if the boundary conditions are bicyclic

      INTEGER, INTENT(IN):: NX,NY,MAXNCL,MAXSIZE
! NX - number of x grid points
! NY - number of y grid points
! MAXNCL - maximum number of clouds
! MAXSIZE - max number of grid points in a single cloud
      INTEGER, INTENT(OUT):: ncl
! ncl - number of cloud structures found
      INTEGER, INTENT(OUT), DIMENSION(MAXNCL):: sizes
! sizes - number of grid points constituting the cloud
      INTEGER, INTENT(IN), DIMENSION(NX,NY):: BINDAT
! BINDAT - binary data 
      INTEGER, INTENT(OUT), DIMENSION(MAXNCL,MAXSIZE):: loci,locj
! loci - x index of grid points constituting cloud
! locj - y index of grid points constituting cloud

      INTEGER:: I,J,residual,index,idir,newi,newj
! I,J,idir - loop variable
! residual - number of points that are cloudy but not part of recognized 
!            clouds
! index - temporary variable labelling grid points within cloud
! newi,newj - position of neighbouring point
      INTEGER, DIMENSION(MAXSIZE):: tempi,tempj
! tempi - temporary x index storing grid points contributing to cloud
! tempj - temporary y index storing grid points contributing to cloud
      INTEGER, DIMENSION(NX,NY):: mask
! mask - flag to signal if column has been tested
!
      LOGICAL:: check
! check - flag to signal if column needs checking
!
! Initialize arrays defining the clouds for this timestep
!
      ncl=0
      sizes=0
      loci=0
      locj=0
!
! Attempt to group cloudy grid points into clouds
!
      do i = 1,nx
         do j =1,ny
            mask(i,j) = 0
         enddo
      enddo
!
      residual=0
      do i = 1,nx
         do j = 1,ny
            if (mask(i,j).eq.0) then
               mask(i,j) = 1  
! signals that we have looked at this point
               check=.true.   
! signals that we need to look at this point
               if (BINDAT(i,j).eq.0) then
                  check=.false.
               endif
               if (check) then  ! have found a cloudy grid point
                  ncl=ncl+1
                  IF (ncl.GT.MAXNCL) GOTO 9002
                  index = 1
                  sizes(ncl) = 1
                  tempi(index)= i
                  tempj(index)= j
                  loci(ncl,sizes(ncl))= i
                  locj(ncl,sizes(ncl))= j
                  do while (index.gt.0)
                     do idir = 1,ndir 
! cycling through all neighbouring points
                        if (cyclicbc) then
                           newi = mod(tempi(index)+ioff(idir)+nx-1, &
                           &                       nx)+1
                           newj = mod(tempj(index)+joff(idir)+ny-1, &
                           &                       ny)+1
                        else
                           newi=tempi(index)+ioff(idir)
                           newj=tempj(index)+joff(idir)
                           newi=min(newi,nx)
                           newj=min(newj,ny)
                           newi=max(newi,1)
                           newj=max(newj,1)
                        end if
                        if (mask(newi,newj).eq.0) then 
! only look at points that we have not already dealt with
                           check=.true. 
! see if we may need to look at this point
                           if (BINDAT(newi,newj).eq.0) then
                              check = .false.
                           endif
                           if (check) then 
! found a cloud point next door to ours 
                              sizes(ncl) = sizes(ncl)+1
                              if (sizes(ncl).gt.maxsize) GOTO 9003
                              index = index+1
                              tempi(index)=newi
                              tempj(index)=newj
                              loci(ncl,sizes(ncl))=newi
                              locj(ncl,sizes(ncl))=newj
                              mask(newi,newj) = 1
                              goto 10
                           endif 
! condition for neighbour to be cloudy
                        endif 
! condition for neighbour to have been already tested			
                        mask(newi,newj)=1
                     enddo      ! loop over neighbours
                     index = index-1
 10                  continue
                  enddo
!                  write(*,*) 'Final size:',sizes(ncl),ncl
!                  write(*,*) 'We have a cloud at these points:'
!                  do idir=1,sizes(ncl)
!                     write(*,*) loci(ncl,idir),locj(ncl,idir)
!                  enddo
                  if (sizes(ncl).LT.minsize) then
! unset this "cloud"
                     residual=residual+sizes(ncl)
                     sizes(ncl)=0
                     ncl=ncl-1
                  endif
               endif            ! condition for cloudy grid point
            endif               ! mask condition
         enddo
      enddo
!
      RETURN
!
 9002 CONTINUE
      WRITE(*,*) 'Too many clouds!'
      GOTO 10000
 9003 CONTINUE
      WRITE(*,*) 'A cloud is too big'
      GOTO 10000
10000 CONTINUE
      WRITE(*,*) 'EXECUTION TERMINATED FROM CONSTRUCT_CLOUDS'
      STOP
!
      END SUBROUTINE CONSTRUCT_CLOUDS
!
      SUBROUTINE IDENTIFY_CLOUDY_POINTS(NX,NY,INDAT,BINDAT)
!
      INTEGER, INTENT(IN):: NX,NY
! NX - number of x grid points
! NY - number of y grid points
      INTEGER, INTENT(OUT), DIMENSION(NX,NY):: BINDAT
! BINDAT - binary data 

      REAL, PARAMETER:: THRESH=1.0E-5
! THRESH - threshold to identify cloud

      REAL, INTENT(IN), DIMENSION(NX,NY):: INDAT
! INDAT - data read in from FILE

      INTEGER:: I,J
! I,J - loop variables
!
      DO I=1,NX
         DO J=1,NY
            IF (INDAT(I,J).GE.THRESH) THEN
               BINDAT(I,J)=1
            ELSE
               BINDAT(I,J)=0
            END IF
         END DO
      END DO
!
      RETURN
!
      END SUBROUTINE IDENTIFY_CLOUDY_POINTS
!
      SUBROUTINE CENTRE_OF_CLOUDS(ncl,MAXNCL,locx,locy,sizes,loci,locj, &
     &                            MAXSIZE,NX,NY)
!
      INTEGER, INTENT(IN):: ncl,MAXNCL,MAXSIZE,NX,NY
! ncl - number of cloud structures 
! MAXNCL - maximum number of clouds
! MAXSIZE - max number of grid points in a single cloud
! NX - number of x grid points
! NY - number of y grid points
      INTEGER, INTENT(IN), DIMENSION(MAXNCL):: sizes
! sizes - number of grid points constituting the cloud
      INTEGER, DIMENSION(MAXNCL,MAXSIZE):: loci,locj
! loci - x index of grid points constituting cloud
! locj - y index of grid points constituting cloud

      REAL, INTENT(OUT), DIMENSION(MAXNCL):: locx,locy
! locx - x location of cloud centre, on a normalized grid
! locy - y location of cloud centre, on a normalized grid

      INTEGER:: n,m 
! n,m - loop variables

      do n=1,ncl
         locx(n)=0.0
         locy(n)=0.0
         do m=1,sizes(n)
! NB: this code would require modification for bicyclic bc's
            locx(n)=locx(n)+loci(n,m)
            locy(n)=locy(n)+locj(n,m)
         end do
         locx(n)=locx(n)/(sizes(n)*NX)
         locy(n)=locy(n)/(sizes(n)*NY)
      end do

      RETURN
! 
      END SUBROUTINE CENTRE_OF_CLOUDS
!
      SUBROUTINE CLOUDSEP(ncl,MAXNCL,locx,locy,sep,RLGE,nnid,nnsep,NUMLOC,LOCSEP)
!
      INTEGER, INTENT(IN):: ncl,MAXNCL,NUMLOC
! ncl - number of cloud structures 
! MAXNCL - maximum number of clouds
! NUMLOC - number of locations for which distance-to-cloud are to be computed
     INTEGER, INTENT(OUT), DIMENSION(MAXNCL):: nnid
! nnid - identifier for the nearest neighbour

      REAL, INTENT(IN):: RLGE
! RLGE - large real number
      REAL:: rtmp
! rtmp - temporary real number
      REAL, INTENT(IN), DIMENSION(MAXNCL):: locx,locy
! locx - x location of cloud centre
! locy - y location of cloud centre
      REAL, INTENT(OUT), DIMENSION(MAXNCL):: nnsep
! nnsep - nearest neighbour separation
      REAL, INTENT(OUT), DIMENSION(NUMLOC):: locsep
! locsep - separation from a grid of locations to the nearest cloud
      REAL, INTENT(OUT), DIMENSION(MAXNCL,MAXNCL):: sep
! sep - separation of the first cloud to its neighbours in the second index
      REAL, DIMENSION(NUMLOC):: locgridx,locgridy
! locgridx - x locations of grid of points chosen for distance-to-cloud calculations
! locgridy - y locations of grid of points chosen for distance-to-cloud calculations

      INTEGER:: n,m,p,q 
! n,m,p,q - loop variables
      INTEGER:: k,ksq
! k - number of locations on regular grid
! ksq - total number of locations considered
      INTEGER, DIMENSION(1):: itemp1
!
! Set up regular grid
!
      k=SQRT(NUMLOC*1.0)
      ksq=k**2
      do p=1,k
        do q=1,k
          locgridx((p-1)*k+q)=(p-1)*(1.0/k)
          locgridy((p-1)*k+q)=(q-1)*(1.0/k)
        end do
      end do
      locsep=RLGE
!
! Determine separations
!
      do n=1,ncl
         do m=1,ncl
            sep(n,m)=sqrt((locx(n)-locx(m))**2+(locy(n)-locy(m))**2)
         end do
         do m=ncl+1,MAXNCL
            sep(n,m)=RLGE
         end do
         sep(n,n)=RLGE
         itemp1=MINLOC(sep(n,:))
         nnid(n)=itemp1(1)
         sep(n,n)=0.0
         nnsep(n)=sep(n,nnid(n))
         do p=1,ksq
            rtmp=sqrt((locx(n)-locgridx(p))**2+(locy(n)-locgridy(p))**2)
            locsep(p)=MIN(locsep(p),rtmp)
         end do  
      end do
!
      RETURN
!
      END SUBROUTINE CLOUDSEP
