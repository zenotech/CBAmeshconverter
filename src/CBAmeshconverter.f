!
      PROGRAM CBAMESHCONVERTER
!
! Declare variables, memory to be assigned later
!
      IMPLICIT NONE
c
      INTEGER:: MI,MJ,MK,MB
c
      DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE:: x,y,z
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: xv,yv,zv
c
      INTEGER,DIMENSION(:,:,:,:),ALLOCATABLE:: numcell,numvert
      INTEGER,DIMENSION(:,:),ALLOCATABLE:: nfacecell,ncellfaces,nvf
      INTEGER,DIMENSION(:),ALLOCATABLE::
     &            i1tag,i2tag,j1tag,j2tag,k1tag,k2tag,nhalofaces,
     &            bound,ntag,cellnum,invec1,invec2,nptsi,nptsj,nptsk
!
      DOUBLE PRECISION:: dis,dx,dy,dz,xx,yy,zz,vnum
!
      INTEGER:: i,j,k,nc,nsym,nii,njj,nkk,nimax,njmax,nkmax,nblock,
     &          ii1,ii2,ii3,ic,nff,nf,i2,nn,i1,nboundtag,nvtot,nb,
     &          icommon,imatch,NCELLS,NFACES,NVERTS,nv,nv1,nv2,nv3,nv4,
     &          ncellshalo
!
       CHARACTER*80:: gridfile
!
       LOGICAL:: TWOD
!
      write(6,601)
601   format('Input mesh filename in CBA format:')
      read(5,888) gridfile
888   format((a))
!
      TWOD=.false.
!
      nimax=0
      njmax=0
      nkmax=0
!
! Scan CBA file to work out dimensionality
!
      open(193,file=gridfile)
      read(193,*) nblock,nsym,vnum
      do nb=1,nblock
        read(193,*) nii,njj,nkk
        if(nii.gt.nimax) nimax=nii
        if(njj.gt.njmax) njmax=njj
        if(nkk.gt.nkmax) nkmax=nkk
        do k=0,nkk-1
          do j=0,njj-1
            do i=0,nii-1
              read(193,*) xx,yy,zz
            enddo
          enddo
        enddo
        read(193,*) ii1,ii2,ii3
        read(193,*) ii1,ii2,ii3
        read(193,*) ii1,ii2,ii3
        read(193,*) ii1,ii2,ii3
        read(193,*) ii1,ii2,ii3
        read(193,*) ii1,ii2,ii3
      enddo
      close(193)
!
      IF(nkk.eq.1) then
        TWOD=.TRUE.
        nkmax=2
      ENDIF
!
      MI=nimax
      Mj=njmax
      Mk=nkmax
      MB=nblock
!
! Allocate memory based on CBA file size
!
      write(6,602)
602   format('ALLOCATING DATA')
!
      ALLOCATE(numcell(MI,MJ,MK,MB))
      ALLOCATE(numvert(0:MI,0:MJ,0:MK,MB))
      ALLOCATE(nfacecell(0:MI*MJ*MK*MB*6,2))
      ALLOCATE(ncellfaces(0:MI*MJ*MK*MB,6))
      ALLOCATE(nhalofaces(0:2*MI*MJ*MK*MB))
      ALLOCATE(nvf(0:MI*MJ*MK*MB*6,6))
      ALLOCATE(xv(0:(MI+1)*(MJ+1)*(MK+1)*MB))
      ALLOCATE(yv(0:(MI+1)*(MJ+1)*(MK+1)*MB))
      ALLOCATE(zv(0:(MI+1)*(MJ+1)*(MK+1)*MB))
      ALLOCATE(bound(0:MI*MJ*MK*MB))
      ALLOCATE(ntag(0:MI*MJ*MK*MB*6))
      ALLOCATE(cellnum(0:MI*MJ*MK*MB))
!
      bound=0
      nfacecell=-100
!
      ALLOCATE(x(0:MI,0:MJ,0:MK,MB))
      ALLOCATE(y(0:MI,0:MJ,0:MK,MB))
      ALLOCATE(z(0:MI,0:MJ,0:MK,MB))

      ALLOCATE(i1tag(MB))
      ALLOCATE(i2tag(MB))
      ALLOCATE(j1tag(MB))
      ALLOCATE(j2tag(MB))
      ALLOCATE(k1tag(MB))
      ALLOCATE(k2tag(MB))
      ALLOCATE(nptsi(MB))
      ALLOCATE(nptsj(MB))
      ALLOCATE(nptsk(MB))
!
! Read in data from CBA file
!
      open(193,file=gridfile)
      read(193,*) nblock,nsym,vnum
      do nb=1,nblock
        read(193,*) nptsi(nb),nptsj(nb),nptsk(nb)
        do k=0,nptsk(nb)-1
          do j=0,nptsj(nb)-1
            do i=0,nptsi(nb)-1
              read(193,*) x(i,j,k,nb),y(i,j,k,nb),z(i,j,k,nb)
            enddo
          enddo
        enddo
        read(193,*) i1tag(nb),ii2,ii3
        read(193,*) i2tag(nb),ii2,ii3
        read(193,*) j1tag(nb),ii2,ii3
        read(193,*) j2tag(nb),ii2,ii3
        read(193,*) k1tag(nb),ii2,ii3
        read(193,*) k2tag(nb),ii2,ii3
      enddo
      close(193)
!
! Check if case is 2D or 3D
!
      IF(TWOD) THEN
        do nb=1,nblock
          nptsk(nb)=2
          do j=0,nptsj(nb)-1
            do i=0,nptsi(nb)-1
              x(i,j,1,nb)=x(i,j,0,nb)
              y(i,j,1,nb)=y(i,j,0,nb) - 1.0d0
              z(i,j,1,nb)=z(i,j,0,nb)
            enddo
          enddo
        enddo
      ENDIF
!
!  SET CELL NUMBERS
!
      nc=0
      do nb=1,nblock
      do k=1,nptsk(nb)-1
        do j=1,nptsj(nb)-1
          do i=1,nptsi(nb)-1
            nc=nc+1
            numcell(i,j,k,nb)=nc
          enddo
        enddo
      enddo
      enddo
!
      NCELLS=nc
!
      write(6,604) ncells
!
!  SET VERTEX NUMBERS
!  CHECK FOR COMMON MESH POINTS:
!  ONLY CHECKS NODES THAT ARE ON BOUNDARIES
!
      nv=0
      nvtot=0
      do nb=1,nblock
        write(6,603) nb
603     format('SCANNING BLOCK',i9,' FOR COMMON MESH NODES')
      do k=0,nptsk(nb)-1
        do j=0,nptsj(nb)-1
          do i=0,nptsi(nb)-1
            nvtot=nvtot+1
            IF(nv.gt.1) THEN
              IF(i.eq.0.or. i.eq.nptsi(nb)-1.or.
     &           j.eq.0.or. j.eq.nptsj(nb)-1.or.
     &           k.eq.0.or. k.eq.nptsk(nb)-1) THEN
              icommon=0
              DO nv1=1,nv
                IF(bound(nv1).eq.1) THEN
                  dx=x(i,j,k,nb)-xv(nv1)
                  dy=y(i,j,k,nb)-yv(nv1)
                  dz=z(i,j,k,nb)-zv(nv1)
                  dis=dsqrt(dx*dx+dy*dy+dz*dz)
                  if(dis.lt.0.000000001d0) then
                    numvert(i,j,k,nb)=nv1
                    icommon=1
                  endif
                ENDIF
              ENDDO
              if(icommon.eq.0) then
                nv=nv+1
                xv(nv)=x(i,j,k,nb)
                yv(nv)=y(i,j,k,nb)
                zv(nv)=z(i,j,k,nb)
                numvert(i,j,k,nb)=nv
              IF(i.eq.0.or. i.eq.nptsi(nb)-1.or.
     &           j.eq.0.or. j.eq.nptsj(nb)-1.or.
     &           k.eq.0.or. k.eq.nptsk(nb)-1) THEN
                      bound(nv)=1
                   ENDIF
              endif
              ELSE
                nv=nv+1
                xv(nv)=x(i,j,k,nb)
                yv(nv)=y(i,j,k,nb)
                zv(nv)=z(i,j,k,nb)
                numvert(i,j,k,nb)=nv
              IF(i.eq.0.or. i.eq.nptsi(nb)-1.or.
     &           j.eq.0.or. j.eq.nptsj(nb)-1.or.
     &           k.eq.0.or. k.eq.nptsk(nb)-1) THEN
                      bound(nv)=1
                   ENDIF
              ENDIF
            ELSE
              nv=nv+1
              xv(nv)=x(i,j,k,nb)
              yv(nv)=y(i,j,k,nb)
              zv(nv)=z(i,j,k,nb)
              numvert(i,j,k,nb)=nv
              IF(i.eq.0.or. i.eq.nptsi(nb)-1.or.
     &           j.eq.0.or. j.eq.nptsj(nb)-1.or.
     &           k.eq.0.or. k.eq.nptsk(nb)-1) THEN
                      bound(nv)=1
                   ENDIF
            ENDIF
          enddo
        enddo
      enddo
      enddo
!
      NVERTS=nv
!
      write(6,606) nvtot-nv
      write(6,605) nverts
604   format('NCELLS (REAL) =',i8)
606   format('REMOVED',i16,' COMMON NODES' )
605   format('NNODES        =',i8)
!
259   format(3f18.12)
333   format(9i9)
!
      open(359,file='cellfacedata.dat')
      write(359,*) ncells
      do nb=1,nblock
      do k=1,nptsk(nb)-1
        do j=1,nptsj(nb)-1
          do i=1,nptsi(nb)-1
            write(359,333) numcell(i,j,k,nb),6,4
            write(359,333)
     &        numvert(i-1,j-1,k-1,nb),numvert(i-1,j-1,k,nb),
     &        numvert(i-1,j,k,nb),numvert(i-1,j,k-1,nb)
            if(i.eq.1.and.i1tag(nb).ne.2) then
              write(359,333) i1tag(nb)
            else
              write(359,333) 2
            endif
            write(359,333)
     &        numvert(i,j-1,k-1,nb),numvert(i,j-1,k,nb),
     &        numvert(i,j,k,nb),numvert(i,j,k-1,nb)
            if(i.eq.nptsi(nb)-1.and.i2tag(nb).ne.2) then
              write(359,333) i2tag(nb)
            else
              write(359,333) 2
            endif
            write(359,333)
     &        numvert(i-1,j-1,k-1,nb),numvert(i,j-1,k-1,nb),
     &        numvert(i,j-1,k,nb),numvert(i-1,j-1,k,nb)
            if(j.eq.1.and.j1tag(nb).ne.2) then
              write(359,333) j1tag(nb)
            else
              write(359,333) 2
            endif
            write(359,333)
     &        numvert(i-1,j,k-1,nb),numvert(i,j,k-1,nb),
     &        numvert(i,j,k,nb),numvert(i-1,j,k,nb)
              if(j.eq.nptsj(nb)-1.and.j2tag(nb).ne.2) then
                write(359,333) j2tag(nb)
              else
                write(359,333) 2
              endif
            write(359,333)
     &        numvert(i-1,j-1,k-1,nb),numvert(i,j-1,k-1,nb),
     &        numvert(i,j,k-1,nb),numvert(i-1,j,k-1,nb)
              if(k.eq.1.and.k1tag(nb).ne.2) then
                write(359,333) k1tag(nb)
              else
                write(359,333) 2
              endif
            write(359,333)
     &        numvert(i-1,j-1,k,nb),numvert(i,j-1,k,nb),
     &        numvert(i,j,k,nb),numvert(i-1,j,k,nb)
              if(k.eq.nptsk(nb)-1.and.k2tag(nb).ne.2) then
                write(359,333) k2tag(nb)
              else
                write(359,333) 2
              endif
          enddo
        enddo
      enddo
      enddo
      close(359)
!
! CHECK FACE NODES FOR COMMON FACES:
!
      ALLOCATE(invec1(4))
      ALLOCATE(invec2(4))
!
      nf=0

      open(359,file='cellfacedata.dat')
      read(359,*) ncells
      do nc=1,ncells
         read(359,*) cellnum(nc),i1,i2
!********************
         do nn=1,6
!********************
         read(359,*) nv1,nv2,nv3,nv4
         read(359,*) nboundtag
         if(nf.eq.0) then
           nf=nf+1
           nvf(nf,1)=nv1
           nvf(nf,2)=nv2
           nvf(nf,3)=nv3
           nvf(nf,4)=nv4
           nfacecell(nf,1)=cellnum(nc)
           ntag(nf)=nboundtag
         else
           imatch=0
           if(nboundtag.eq.2) then
           do nff=1,nf
              invec1(1)=nv1
              invec1(2)=nv2
              invec1(3)=nv3
              invec1(4)=nv4
              call isort(invec1,4)
              invec2(1)=nvf(nff,1)
              invec2(2)=nvf(nff,2)
              invec2(3)=nvf(nff,3)
              invec2(4)=nvf(nff,4)
              call isort(invec2,4)
!
              if(invec1(1).eq.invec2(1).and.invec1(2).eq.invec2(2).and.
     &           invec1(3).eq.invec2(3).and.invec1(4).eq.invec2(4)) then

                nfacecell(nff,2)=cellnum(nc)
                ntag(nff)=nboundtag
                imatch=1
              endif
           enddo
           endif
              if(imatch.eq.0) then
                nf=nf+1
                nvf(nf,1)=nv1
                nvf(nf,2)=nv2
                nvf(nf,3)=nv3
                nvf(nf,4)=nv4
                nfacecell(nf,1)=cellnum(nc)
                ntag(nf)=nboundtag
              endif
         endif
!********************
         enddo
!********************
       enddo
!
       DEALLOCATE(invec1,invec2)
!
       NFACES=nf
      write(6,607) nfaces
607   format('NFACES        =',i8)
!
       do nc=1,ncells
         ic=0
         do nf=1,nfaces
           if(nfacecell(nf,1).eq.nc) then
             ic=ic+1
             ncellfaces(nc,ic)=nf
           elseif(nfacecell(nf,2).eq.nc) then
             ic=ic+1
             ncellfaces(nc,ic)=nf
           endif
         enddo
       enddo
!
! CHECK FOR HALO CELLS TO ADD ON
!
       ncellshalo=0
       do nf=1,nfaces
         if(nfacecell(nf,2).eq.-100) then
           ncellshalo=ncellshalo+1
           nfacecell(nf,2)=ncells+ncellshalo
           nhalofaces(ncells+ncellshalo)=nf
         endif
       enddo
!
      write(6,608) ncellshalo
608   format('NCELLS (HALO) =',i8)
!
! OUTPUT GENERAL DATA
!
       open(100,file='facescells.dat')
       write(100,*) NFACES,NCELLS,ncellshalo
       close(100)
!
       write(6,609)
609    format('WRITING OUTPUT FILES')
!
! WRITE OUTPUT FILES
!
       open(100,file='cellFace.txt')
       write(100,*) ncells+ncellshalo,nfaces
       do nc=1,ncells
         do nf=1,6
           write(100,*) ncellfaces(nc,nf)-1
         enddo
       enddo
       do nc=ncells+1,ncells+ncellshalo
         write(100,*) nhalofaces(nc)-1
       enddo
       close(100)
!
       open(100,file='cellType.txt')
       write(100,*) ncells+ncellshalo
       do nc=1,ncells+ncellshalo
         write(100,*) 208
       enddo
       close(100)
!
       open(100,file='faceBC.txt')
       write(100,*) nfaces
       do nf=1,nfaces
         if(ntag(nf).eq.-2) then
           write(100,*) 7 !Symmetry
         elseif(ntag(nf).eq.-1.or.ntag(nf).eq.0) then
           write(100,*) 3 !Wall
         elseif(ntag(nf).eq.1) then
           write(100,*) 9 !Farfield
         elseif(ntag(nf).eq.2) then
           write(100,*) 2 !Interior
         endif
       enddo
       close(100)
!
      open(100,file='faceInfo.txt')
      write(100,*) nfaces
      do nf=1,nfaces
        if(ntag(nf).eq.-2) then
          write(100,*) 1,0
!          if(ntag(nf-1).eq.-2) then
!            write(100,*) -1,0
!          else
!            write(100,*) 1,0 !Symmetry
!          endif
        elseif(ntag(nf).eq.-1.or.ntag(nf).eq.0) then
          write(100,*) 4,0 !Wall
        elseif(ntag(nf).eq.1) then
          write(100,*) 2,0 !Farfield
        elseif(ntag(nf).eq.2) then
          write(100,*) 0,0 !Interior
        endif
      enddo
      close(100)
!
       open(100,file='faceCell.txt')
       write(100,*) nfaces
       do nf=1,nfaces
         write(100,*) nfacecell(nf,1)-1,nfacecell(nf,2)-1
       enddo
       close(100)
!
       open(100,file='faceNodes.txt')
       write(100,*) nfaces,4
       do nf=1,nfaces
           write(100,*) nvf(nf,1)-1
           write(100,*) nvf(nf,2)-1
           write(100,*) nvf(nf,3)-1
           write(100,*) nvf(nf,4)-1
       enddo
       close(100)
!
       open(100,file='faceType.txt')
       write(100,*) nfaces
       do nf=1,nfaces
           write(100,*) 4
       enddo
       close(100)
!
       open(100,file='nodeVertex.txt')
       write(100,*) nverts
       do nv=1,nverts
         write(100,259) xv(nv),yv(nv),zv(nv)
       enddo
       close(100)
!
       close(209)
       close(219)
       close(109)
!
      DEALLOCATE(numcell)
      DEALLOCATE(numvert)
      DEALLOCATE(nfacecell)
      DEALLOCATE(ncellfaces)
      DEALLOCATE(nhalofaces)
      DEALLOCATE(nvf)
      DEALLOCATE(bound)
      DEALLOCATE(ntag)
      DEALLOCATE(cellnum)
      DEALLOCATE(x)
      DEALLOCATE(y)
      DEALLOCATE(z)
      DEALLOCATE(xv)
      DEALLOCATE(yv)
      DEALLOCATE(zv)

      DEALLOCATE(i1tag)
      DEALLOCATE(i2tag)
      DEALLOCATE(j1tag)
      DEALLOCATE(j2tag)
      DEALLOCATE(k1tag)
      DEALLOCATE(k2tag)
      DEALLOCATE(nptsi)
      DEALLOCATE(nptsj)
      DEALLOCATE(nptsk)
!
      call system('rm cellfacedata.dat')
      stop
!
      end program cbameshconverter
!
!*********************************************************************
!**********************************************************************
!
       SUBROUTINE ISORT(invec,nsize)
!
       IMPLICIT NONE
       INTEGER:: nsize
       INTEGER,DIMENSION(nsize):: invec
       INTEGER:: i,j,itemp
       LOGICAL:: swapped

       do j=nsize-1,1,-1
         swapped=.FALSE.
         do i = 1,j
           if(invec(i).gt.invec(i+1)) then
             itemp=invec(i)
             invec(i)=invec(i+1)
             invec(i+1)=itemp
             swapped=.TRUE.
           endif
         enddo
         if(swapped.eqv..FALSE.) goto 999
       enddo
999    return
       end subroutine isort
