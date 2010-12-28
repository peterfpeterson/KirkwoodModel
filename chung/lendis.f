      program lendis
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C-
C-   Purpose and Methods : 
C-      Calculates length distribution
C-      using the final output of ALLSIG.F for sigma
C-      Output is written as a topdrawer file.
C-
C-      For now, we have in mind only the case of A_1-x B_x C
C-
C-
C-   Inputs  : 
C-   Outputs : 
C-   Controls: 
C-
C-   Created  10-MAR-1996   Jean Soo Chung
C-   Modified  2-AUG-1997
C-       nnn peak resolution added.  It uses the same numbering scheme
C-       as that of "allsig.f" and the same nn assignment, the 
C-       subroutine "findnn()"
C-
c-----------------------------------------------------------------------
      implicit none

      integer maxl,maxn,maxd
      parameter (maxl=4,maxn=maxl*maxl*maxl*4,maxd=maxn*2*3)
      real*8  unitv(4,2,3), phi(maxd)
      integer nn(maxn,2,4), ldnn

      integer maxfile,maxout,maxbl
      real*8  pi,pix2,pix4
      parameter (maxfile=20,maxout=2000,pi=3.14159265358979d0
     *          ,pix2=2*pi,pix4=4*pi)
      real*8  r,sigma,x,y,radius(2,2),mass(2,2), alpha,beta,lambda
     *       ,le,a0,del,del2,temp
     *       ,rmax
      integer nfile,unit,nout,lx,ly,lz,iseed,na,natom
     *       ,ibl,isl,iel,jbl,jsl,jel
c     *       ,ibl(maxn),isl(maxn),iel(maxn)
c     *       ,jbl(maxn),jsl(maxn),jel(maxn)
     *       ,nin,nerr,i,j,hmax,his(0:maxout)
     *       ,h11(0:maxout),h22(0:maxout),h12(0:maxout)
     *       ,i11(0:maxout),i22(0:maxout),i12(0:maxout)
     *       ,j11(0:maxout),j22(0:maxout),j12(0:maxout),j21(0:maxout)
      character elem(2,2)*2,infile(maxfile)*40,outfile*40,ans*20,job*3
c-----------------------------------------------------------------------
c     initial setup
c-----------------------------------------------------------------------
c     read input files
c     ------------------------------------------------------------------
110   write(6,*) 'enter number of input files'
      read(5,*,err=110) nfile
      if(nfile.gt.maxfile) stop '*** too many files, increase MAXFILE'
      do j=1,nfile
120      write(6,'('' enter input file name'',i3)') j
         read(5,'(a)') infile(j)
         unit=10+j
         open(unit=unit,file=infile(j),status='old',err=120)
         call readhead(unit,elem,x,y,radius,mass,lx,ly,lz,temp)
      end do

      na=lx*ly*lz*4
      natom=na*2
c     set parameters for output
c     -------------------------
c---  maximum length
      le= (x*radius(1,1)+(1-x)*radius(1,2))
     *   +(y*radius(2,1)+(1-y)*radius(2,2))
      a0=le/sqrt(3.d0)*4
      rmax=dfloat(int(a0/2*sqrt(3.d0)*10))/10
      if(nout.gt.maxout) nout=maxout
200   write(6,30) le,rmax
30    format(' nn distance =',f6.3/
     *       ' [NN] peak partial distribution'/
     *       ' [NNN] peak partial distribution'/
     *       ' Max length for distribution? [',f7.3,' ]')
      read(5,'(a)') ans
      if(ans(1:3).eq.'NNN' .or. ans(1:3).eq.'nnn') then
         job='nnn'
         rmax=int(le*1.4*10)/10
      else if(ans(1:2).eq.'NN' .or. ans(1:2).eq.'nnn') then
         job='nn '
      else if(ans(1:1).ne.' ') then
         job='tot'
         read(ans,*,err=200) rmax
      end if
c---  number of bins
      nout=rmax*100
250   write(6,'('' Number of bins? ['',i4,'']'')') nout
      read(5,'(a)') ans
      if(ans(1:1).ne.' ') then
         read(ans,*,err=250) nout
      end if
      if(nout.gt.maxout) then
         nout=maxout
         write(6,'(/'' *** nout is set to'',i6/)') maxout
      end if
c     ------------------------------------------------------------------
c     output setup
c     ------------------------------------------------------------------
350   write(6,*) 'enter output TopDrawer file name'
      read(5,'(a)') outfile
      open(unit=8,file=outfile,status='new',err=350)
      write(8,50) elem(1,1),1-x,elem(1,2),x,elem(2,1),1-y,elem(2,2),y,
     *             rmax
50    format(1x,'set font duplex'
     *      /1x,'title top ''length distribution from ',
     *      2(a2,'0',f4.2,'1',a2,'0',f4.2,'1'),'''',
     *      /1x,'case      ''',27x,
     *      'X    X  X    X  X    X  X    X'''
     *      /1x,'title left ''histogram''',
     *      /1x,'title bottom ''r''',
     *      /1x,'set limit x 0 ',f7.3)
c-----------------------------------------------------------------------
c     read and accumulate length distribution
c-----------------------------------------------------------------------
c     initialize histogram
c     --------------------
      do i=0,nout
         his(i)=0
         h11(i)=0
         h22(i)=0
         h12(i)=0
         i11(i)=0
         i22(i)=0
         i12(i)=0
         j11(i)=0
         j22(i)=0
         j12(i)=0
      end do
      hmax=0
c     read in and accumulate length distribution
c     ------------------------------------------
      do unit=11,10+nfile
         nerr=0
         nin=0
         do while(nerr.lt.100)
            read(unit,*,err=500,end=600) ibl,isl,iel,jbl,jsl,jel,r,sigma
            if(r.lt.rmax) then
               j=(r/rmax)*nout

c---           total histogram
               his(j)=his(j)+1
               if(his(j).gt.hmax) hmax=his(j)

c---           between atoms in sublattice 1
               if(isl.eq.1 .and. jsl.eq.1) then
                  if(iel.eq.1 .and. jel.eq.1) then
                     h11(j)=h11(j)+1
                  else if(iel.eq.2 .and. jel.eq.2) then
                     h22(j)=h22(j)+1
                  else
                     h12(j)=h12(j)+1
                  end if

c---           between atoms in sublattice 2
               else if(isl.eq.2 .and. jsl.eq.2) then
                  if(iel.eq.1 .and. jel.eq.1) then
                     i11(j)=i11(j)+1
                  else if(iel.eq.2 .and. jel.eq.2) then
                     i22(j)=i22(j)+1
                  else
                     i12(j)=i12(j)+1
                  end if

c---           between sublattice 1 and sublattice 2 atoms
               else if(isl.eq.1) then
                  if(iel.eq.1 .and. jel.eq.1) then
                     j11(j)=j11(j)+1
                  else if(iel.eq.2 .and. jel.eq.2) then
                     j22(j)=j22(j)+1
                  else if(iel.eq.1) then
                     j12(j)=j12(j)+1
                  else
                     j21(j)=j21(j)+1
                  end if
               else if(isl.eq.2) then
                  if(iel.eq.1 .and. jel.eq.1) then
                     j11(j)=j11(j)+1
                  else if(iel.eq.2 .and. jel.eq.2) then
                     j22(j)=j22(j)+1
                  else if(iel.eq.2) then
                     j12(j)=j12(j)+1
                  else
                     j21(j)=j21(j)+1
                  end if
               endif
            end if
            nin=nin+1
            go to 550
c---        error handling
500         nerr=nerr+1
550         continue
         end do
600      close(unit)
         write(6,60) nin,infile(unit-10)
60       format(i9,' lines read from file ',a30)
         if(nerr.ne.0) write(6,70) nerr
70       format(// ' *****',i8,' errors in reading, lines ignored'//)
      end do

c     write length distribution
c     -------------------------
      nin=hmax/15
      write(8,80) hmax+nin, hmax-nin,hmax-nin
c       A-A
     *            ,hmax-2*nin,hmax-2*nin,elem(1,1),elem(1,1)
c       B-B
     *            ,hmax-2*nin,hmax-2*nin,elem(1,2),elem(1,2)
c       A-B
     *            ,hmax-3*nin,hmax-3*nin,elem(1,1),elem(1,2)
c       C-C
     *            ,hmax-4*nin,hmax-4*nin,elem(2,1),elem(2,1)
c       A-C
     *            ,hmax-5*nin,hmax-5*nin,elem(1,1),elem(2,1)
c       B-C
     *            ,hmax-6*nin,hmax-6*nin,elem(1,2),elem(2,1)
80    format(' set limit y 0',i6/
     *       ' set symbol 1O; 0.1 ',i6,'; plot'/
     *       ' title 0.2 ',i6,' data ''all'''/
     *       ' set symbol 2O; 0.1 ',i6,'; plot'/
     *       ' title 0.2 ',i6,' data ''',a2,'-',a2,''''/
     *       ' set symbol 3O; 0.1 ',i6,'; plot'/
     *       ' title 0.2 ',i6,' data ''',a2,'-',a2,''''/
     *       ' set symbol 4O; 0.1 ',i6,'; plot'/
     *       ' title 0.2 ',i6,' data ''',a2,'-',a2,''''/
     *       ' set symbol 5O; 0.1 ',i6,'; plot'/
     *       ' title 0.2 ',i6,' data ''',a2,'-',a2,''''/
     *       ' set symbol 6O; 0.1 ',i6,'; plot'/
     *       ' title 0.2 ',i6,' data ''',a2,'-',a2,''''/
     *       ' set symbol 7O; 0.1 ',i6,'; plot'/
     *       ' title 0.2 ',i6,' data ''',a2,'-',a2,'''')
      del=rmax/nout
      del2=del/2
      write(8,*) 'set symbol 1O'
      do i=1,nout-1
         if(his(i-1).ne.0 .or. his(i).ne.0 .or. his(i+1).ne.0) 
     *      write(8,'(f8.4,i8)') del*i+del2,his(i)
      end do
      write(8,*) 'plot'
      write(8,*) 'join'
      write(8,*) 'set symbol 2O'
      do i=1,nout-1
         if(h11(i-1).ne.0 .or. h11(i).ne.0 .or. h11(i+1).ne.0) 
     *      write(8,'(f8.4,i8)') del*i+del2,h11(i)
      end do
      write(8,*) 'plot'
      write(8,*) 'join'
      write(8,*) 'set symbol 3O'
      do i=1,nout-1
         if(h22(i-1).ne.0 .or. h22(i).ne.0 .or. h22(i+1).ne.0) 
     *      write(8,'(f8.4,i8)') del*i+del2,h22(i)
      end do
      write(8,*) 'plot'
      write(8,*) 'join'
      write(8,*) 'set symbol 4O'
      do i=1,nout-1
         if(h12(i-1).ne.0 .or. h12(i).ne.0 .or. h12(i+1).ne.0) 
     *      write(8,'(f8.4,i8)') del*i+del2,h12(i)
      end do
      write(8,*) 'plot'
      write(8,*) 'join'
      write(8,*) 'set symbol 5O'
      do i=1,nout-1
         if(i11(i-1).ne.0 .or. i11(i).ne.0 .or. i11(i+1).ne.0) 
     *      write(8,'(f8.4,i8)') del*i+del2,i11(i)
      end do
      write(8,*) 'plot'
      write(8,*) 'join'
      write(8,*) 'set symbol 6O'
      do i=1,nout-1
         if(j11(i-1).ne.0 .or. j11(i).ne.0 .or. j11(i+1).ne.0) 
     *      write(8,'(f8.4,i8)') del*i+del2,j11(i)
      end do
      write(8,*) 'plot'
      write(8,*) 'join'
      write(8,*) 'set symbol 7O'
      do i=1,nout-1
         if(j21(i-1).ne.0 .or. j21(i).ne.0 .or. j21(i+1).ne.0) 
     *      write(8,'(f8.4,i8)') del*i+del2,j21(i)
      end do
      write(8,*) 'plot'
      write(8,*) 'join'
      stop
      end


      subroutine readhead(unit,elem,x,y,radius,mass,lx,ly,lz,temp)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read header of file with unit number 'unit'
c     compare parameters with existing values
c     and return necessary parameters
c-----------------------------------------------------------------------
      implicit none
      integer unit,lx,ly,lz,iseed, i,j,flag
      integer lx0,ly0,lz0
      real*8  radius(2,2),mass(2,2), alpha(2,4),beta(4,4),lambda, 
     *        x,y,temp,wmax,xmax,tmax
      real*8  radius0(2,2),mass0(2,2), alpha0(2,4),beta0(4,4),lambda0,
     *        x0,y0,temp0,wmax0,xmax0,tmax0
      character elem(2,2)*2,elem0(2,2)*2

      read(unit,'(4a2)') ((elem(i,j),j=1,2),i=1,2)
      read(unit,*) x,y
      read(unit,*) lx,ly,lz
      read(unit,*) ((radius(i,j),j=1,2),i=1,2)
      read(unit,*) ((mass(i,j),j=1,2),i=1,2)
      read(unit,*) (alpha(1,i),i=1,4)
      read(unit,*) (alpha(2,i),i=1,4)
      read(unit,*) (beta(1,i),i=1,4)
      read(unit,*) (beta(2,i),i=1,4)
      read(unit,*) (beta(3,i),i=1,4)
      read(unit,*) (beta(4,i),i=1,4)
      read(unit,*) lambda
      read(unit,*) iseed
      read(unit,*) temp
      read(unit,*) wmax,xmax,tmax
      read(unit,*)
      read(unit,*)

c---  if it is not a first call, compare with original values
      flag=0
10    format(' *** wrong elements : ',a2,1x,a2)
20    format(' *** wrong',a7,' : ',i9,1x,i9)
30    format(' *** wrong',a7,' : ',f9.4,1x,f9.4)
c      if(elem0(1,1).ne.' ') then
      if(unit.ne.11) then
         do i=1,2
         do j=1,2
            if(elem0(i,j).ne.elem(i,j)) then
               write(6,10) elem0(i,j),elem(i,j)
               flag=1
            end if
            if(radius0(i,j).ne.radius(i,j)) then
               write(6,30) 'radii',radius0(i,j),radius(i,j)
               flag=1
            end if
            if(mass0(i,j).ne.mass(i,j)) then
               write(6,30) 'masses',mass0(i,j),mass(i,j)
               flag=1
            end if
         end do
         end do
         do j=1,4
            do i=1,2
               if(alpha0(i,j).ne.alpha(i,j)) then
                  write(6,30) 'alpha',alpha0(i,j),alpha(i,j)
                  flag=1
               end if
            end do
            do i=1,4
               if(beta0(i,j).ne.beta(i,j)) then
                  write(6,30) 'beta',beta0(i,j),beta(i,j)
                  flag=1
               end if
            end do
         end do
         if(x0.ne.x) then
            write(6,30) 'x',x0,x
            flag=1
         end if
         if(y0.ne.y) then
            write(6,30) 'y',y0,y
            flag=1
         end if
         if(lx0.ne.lx) then
            write(6,20) 'lx',lx0,lx
            flag=1
         end if
         if(ly0.ne.ly) then
            write(6,20) 'ly',ly0,ly
            flag=1
         end if
         if(lz0.ne.lz) then
            write(6,20) 'lz',lz0,lz
            flag=1
         end if
         if(lambda0.ne.lambda) then
            write(6,30) 'lambda',lambda0,lambda
            flag=1
         end if
         if(wmax0.ne.wmax) then
            write(6,30) 'wmax',wmax0,wmax
            flag=1
         end if
         if(xmax0.ne.xmax) then
            write(6,30) 'xmax',xmax0,xmax
            flag=1
         end if
         if(tmax0.ne.tmax) then
            write(6,30) 'tmax',tmax0,tmax
            flag=1
         end if
      end if
      if(flag.ne.0) stop

c---  save current values as original values
      do i=1,2
      do j=1,2
         elem0(i,j)=elem(i,j)
         radius0(i,j)=radius(i,j)
         mass0(i,j)=mass(i,j)
      end do
      end do
      do j=1,4
         do i=1,2
            alpha0(i,j)=alpha(i,j)
         end do
         do i=1,4
            beta0(i,j)=beta(i,j)
         end do
      end do
      x0=x
      y0=y
      lx0=lx
      ly0=ly
      lz0=lz
      lambda0=lambda
      wmax0=wmax
      xmax0=xmax
      tmax0=tmax
      return
      end


      subroutine findnn(lx,ly,lz,nn,ldnn,unitv,phi)
cccc--------------------------------------------------------------------
c     assign nearest neighbor array
c     l   : # of cubic unit cells
c     nn  : indicate 4 neighbors of A & B sublattices
c     ldnn: leading dimension of nn in the calling routine
c     phi : work space for coordinates of atoms
c-----------------------------------------------------------------------
      implicit none
      integer lx,ly,lz,ldnn,nn(ldnn,2,4),
     *        na,ll,ix,iy,iz,ix1,iy1,iz1, ibl,ii,mx,my,mz
      real*8  unitv(4,2,3),phi(ldnn*6),sqrt3
      ll=lx*ly
      na=lx*ly*lz*4
c-----------------------------------------------------------------------
c     sweep through lattice
c-----------------------------------------------------------------------
      do iz=0,lz-1
         mz=iz*4
      do iy=0,ly-1
         my=iy*4
      do ix=0,lx-1
         mx=ix*4
c     ------------------------------------------------------------------
c     A sublattice
c     ------------------------------------------------------------------
         ix1=ix-1
         iy1=iy-1
         iz1=iz-1
c        atom #1
           ibl=(ix+iy*lx+iz*ll)*4+1
           ii=(ibl-1)*6
           phi(ii+1)=mx
           phi(ii+2)=my
           phi(ii+3)=mz
           nn(ibl,1,1)=ibl
           nn(ibl,1,2)=(ix1+iy1*lx+iz*ll)*4+2
           nn(ibl,1,3)=(ix1+iy*lx+iz1*ll)*4+3
           nn(ibl,1,4)=(ix+iy1*lx+iz1*ll)*4+4
c        atom #2
           ibl=ibl+1
           ii=(ibl-1)*6
           phi(ii+1)=mx+2
           phi(ii+2)=my+2
           phi(ii+3)=mz
           nn(ibl,1,1)=ibl
           nn(ibl,1,2)=ibl-1
           nn(ibl,1,3)=(ix+iy*lx+iz1*ll)*4+4
           nn(ibl,1,4)=(ix+iy*lx+iz1*ll)*4+3
c        atom #3
           ibl=ibl+1
           ii=(ibl-1)*6
           phi(ii+1)=mx+2
           phi(ii+2)=my
           phi(ii+3)=mz+2
           nn(ibl,1,1)=ibl
           nn(ibl,1,2)=(ix+iy1*lx+iz*ll)*4+4
           nn(ibl,1,3)=ibl-2
           nn(ibl,1,4)=(ix+iy1*lx+iz*ll)*4+2
c        atom #4
           ibl=ibl+1
           ii=(ibl-1)*6
           phi(ii+1)=mx
           phi(ii+2)=my+2
           phi(ii+3)=mz+2
           nn(ibl,1,1)=ibl
           nn(ibl,1,2)=(ix1+iy*lx+iz*ll)*4+3
           nn(ibl,1,3)=(ix1+iy*lx+iz*ll)*4+2
           nn(ibl,1,4)=ibl-3
c     ------------------------------------------------------------------
c     B sublattice
c     ------------------------------------------------------------------
         ix1=ix+1
         iy1=iy+1
         iz1=iz+1
c        atom #1
           ibl=(ix+iy*lx+iz*ll)*4+1
           ii=(ibl-1)*6
           phi(ii+4)=phi(ii+1)+1
           phi(ii+5)=phi(ii+2)+1
           phi(ii+6)=phi(ii+3)+1
           nn(ibl,2,1)=ibl
           nn(ibl,2,2)=ibl+1
           nn(ibl,2,3)=ibl+2
           nn(ibl,2,4)=ibl+3
c        atom #2
           ibl=ibl+1
           ii=(ibl-1)*6
           phi(ii+4)=phi(ii+1)+1
           phi(ii+5)=phi(ii+2)+1
           phi(ii+6)=phi(ii+3)+1
           nn(ibl,2,1)=ibl
           nn(ibl,2,2)=(ix1+iy1*lx+iz*ll)*4+1
           nn(ibl,2,3)=(ix1+iy*lx+iz*ll)*4+4
           nn(ibl,2,4)=(ix+iy1*lx+iz*ll)*4+3
c        atom #3
           ibl=ibl+1
           ii=(ibl-1)*6
           phi(ii+4)=phi(ii+1)+1
           phi(ii+5)=phi(ii+2)+1
           phi(ii+6)=phi(ii+3)+1
           nn(ibl,2,1)=ibl
           nn(ibl,2,2)=(ix1+iy*lx+iz*ll)*4+4
           nn(ibl,2,3)=(ix1+iy*lx+iz1*ll)*4+1
           nn(ibl,2,4)=(ix+iy*lx+iz1*ll)*4+2
c        atom #4
           ibl=ibl+1
           ii=(ibl-1)*6
           phi(ii+4)=phi(ii+1)+1
           phi(ii+5)=phi(ii+2)+1
           phi(ii+6)=phi(ii+3)+1
           nn(ibl,2,1)=ibl
           nn(ibl,2,2)=(ix+iy1*lx+iz*ll)*4+3
           nn(ibl,2,3)=(ix+iy*lx+iz1*ll)*4+2
           nn(ibl,2,4)=(ix+iy1*lx+iz1*ll)*4+1
      end do
      end do
      end do
c-----------------------------------------------------------------------
c     correction for planes parallel to xy-plane
c-----------------------------------------------------------------------
c     planes parallel to xy-plane
c     ------------------------------------------------------------------
      iz=lx*ly*lz*4
      do ix=0,lx-1
      do iy=0,ly-1
c        atom #1 on A sublattice
           ibl=(ix+iy*lx)*4+1
           nn(ibl,1,3)=nn(ibl,1,3)+iz
           nn(ibl,1,4)=nn(ibl,1,4)+iz
c        atom #2 on A sublattice
           ibl=ibl+1
           nn(ibl,1,3)=nn(ibl,1,3)+iz
           nn(ibl,1,4)=nn(ibl,1,4)+iz
c        atom #3 on B sublattice
           ibl=(ix+iy*lx+(lz-1)*ll)*4+3
           nn(ibl,2,3)=nn(ibl,2,3)-iz
           nn(ibl,2,4)=nn(ibl,2,4)-iz
c        atom #4 on B sublattice
           ibl=ibl+1
           nn(ibl,2,3)=nn(ibl,2,3)-iz
           nn(ibl,2,4)=nn(ibl,2,4)-iz
      end do
      end do
c     ------------------------------------------------------------------
c     correction for planes parallel to yz-plane
c     ------------------------------------------------------------------
      ix=lx*4
      do iy=0,ly-1
      do iz=0,lz-1
c        atom #1 on A sublattice
           ibl=(iy*lx+iz*ll)*4+1
           nn(ibl,1,2)=nn(ibl,1,2)+ix
           nn(ibl,1,3)=nn(ibl,1,3)+ix
c        atom #4 on A sublattice
           ibl=ibl+3
           nn(ibl,1,2)=nn(ibl,1,2)+ix
           nn(ibl,1,3)=nn(ibl,1,3)+ix
c        atom #2 on B sublattice
           ibl=((lx-1)+iy*lx+iz*ll)*4+2
           nn(ibl,2,2)=nn(ibl,2,2)-ix
           nn(ibl,2,3)=nn(ibl,2,3)-ix
c        atom #3 on B sublattice
           ibl=ibl+1
           nn(ibl,2,2)=nn(ibl,2,2)-ix
           nn(ibl,2,3)=nn(ibl,2,3)-ix
      end do
      end do
c     ------------------------------------------------------------------
c     correction for planes parallel to zx-plane
c     ------------------------------------------------------------------
      iy=lx*ly*4
      do iz=0,lz-1
      do ix=0,lx-1
c        atom #1 on A sublattice
           ibl=(ix+iz*ll)*4+1
           nn(ibl,1,2)=nn(ibl,1,2)+iy
           nn(ibl,1,4)=nn(ibl,1,4)+iy
c        atom #3 on A sublattice
           ibl=ibl+2
           nn(ibl,1,2)=nn(ibl,1,2)+iy
           nn(ibl,1,4)=nn(ibl,1,4)+iy
c        atom #2 on B sublattice
           ibl=(ix+(ly-1)*lx+iz*ll)*4+2
           nn(ibl,2,2)=nn(ibl,2,2)-iy
           nn(ibl,2,4)=nn(ibl,2,4)-iy
c        atom #4 on B sublattice
           ibl=ibl+2
           nn(ibl,2,2)=nn(ibl,2,2)-iy
           nn(ibl,2,4)=nn(ibl,2,4)-iy
      end do
      end do
c-----------------------------------------------------------------------
c     unit vectors to nn
c        For A sublattice:
c           nna=1 - ( 1/2  1/2  1/2)
c           nna=2 - (-1/2 -1/2  1/2)
c           nna=3 - (-1/2  1/2 -1/2)
c           nna=4 - ( 1/2 -1/2 -1/2)
c        For B sublattice:
c           nnb=1 - (-1/2 -1/2 -1/2)
c           nnb=2 - ( 1/2  1/2 -1/2)
c           nnb=3 - ( 1/2 -1/2  1/2)
c           nnb=4 - (-1/2  1/2  1/2)
c-----------------------------------------------------------------------
      sqrt3=sqrt(1.d0/3)
      unitv(1,1,1)=sqrt3
      unitv(1,1,2)=sqrt3
      unitv(1,1,3)=sqrt3
      unitv(2,1,1)=-sqrt3
      unitv(2,1,2)=-sqrt3
      unitv(2,1,3)=sqrt3
      unitv(3,1,1)=-sqrt3
      unitv(3,1,2)=sqrt3
      unitv(3,1,3)=-sqrt3
      unitv(4,1,1)=sqrt3
      unitv(4,1,2)=-sqrt3
      unitv(4,1,3)=-sqrt3
      do ix=1,4
      do iy=1,3
         unitv(ix,2,iy)=-unitv(ix,1,iy)
      end do
      end do
      return
      end
