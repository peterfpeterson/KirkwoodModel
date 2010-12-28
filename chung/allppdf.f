      program allppdf
C-----------------------------------------------------------------------
C-
C-   Purpose and Methods :
C-      Calculates total or partial Pair Distribution Function
C-      of a zns structure
C-      using the peak width (sigma) from the output of ZNSSIG.
C-      Output is written as a topdrawer file.
C-
C-   Inputs  :
C-   Outputs :
C-   Controls:
C-
C-   Created  8-MAR-1996   Jean Soo Chung
C-
C-   Modified jan.1999
C-    * According to the change in allsig which allows the relaxation 
C-      of the lattice size when force constants are different,
C-      le is read in instead of calculated.
C-    * TopDrawer format is changed to xmgr format
C-
C-----------------------------------------------------------------------
      implicit none
      integer maxfile,nin,nout,maxbl
      real*8  pi,pix2,pix4
      parameter (maxfile=20,nin=10000000,nout=2000,maxbl=500
     *          ,pi=3.14159265358979d0,pix2=2*pi,pix4=4*pi)
      real*8  r0,sigma,sigma2,sigma10,x,y,radius(2,2),mass(2,2)
     *       ,le,a0,del,temp, scfct(2,2),za,zb,zz
     *       ,tt,rmax,r(nout),dr, pdf(nout),rho
      integer ru,wu,nfile,lx,ly,lz,na,natom
     *       ,ibl,isl,iel,jbl,jsl,jel
     *       ,i,j,k,nerr,psl,pel,psl1,pel1,psl2,pel2
      character elem(2,2)*2,infile(maxfile)*30,outfile*30,ans*3,fcn*3

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
         read(5,'(a20)') infile(j)
         ru=10+j
         open(unit=ru,file=infile(j),status='old',err=120)
         call readhead(ru,elem,x,y,radius,mass,lx,ly,lz,temp,le,a0)
      end do

      na=lx*ly*lz*4
c      if(na.gt.maxbl) stop '*** too large lattice, increase maxbl'
      natom=na*2

c     ------------------------------------------------------------------
c     neutron? or x-ray? and proper factor
c     ------------------------------------------------------------------
150   write(6,*) '1) neutron   2) X-ray (Enter number/n/x)'
      read(5,'(a)') ans
      if(ans.eq.'1' .or. ans.eq.'N' .or. ans.eq.'n') then
         call nfact(elem(1,1),scfct(1,1))
         call nfact(elem(1,2),scfct(1,2))
         call nfact(elem(2,1),scfct(2,1))
         call nfact(elem(2,2),scfct(2,2))
      else if(ans.eq.'2' .or. ans.eq.'X' .or. ans.eq.'x') then
         call xfact(elem(1,1),scfct(1,1))
         call xfact(elem(1,2),scfct(1,2))
         call xfact(elem(2,1),scfct(2,1))
         call xfact(elem(2,2),scfct(2,2))
      else
         goto 150
      end if
c     renormalize the scattering strength
      za=(1-x)*scfct(1,1)+x*scfct(1,2)
      zb=(1-y)*scfct(2,1)+y*scfct(2,2)
      zz=(za+zb)/2
      do i=1,2
      do j=1,2
         scfct(i,j)=scfct(i,j)/zz
      end do
      end do
c     ------------------------------------------------------------------
c     rdf, pdf or partial pdf?
c     ------------------------------------------------------------------
250   write(6,*) '    [r] radial distribution function'
      write(6,*) '    [p] pair   distribution function or reduced rdf'
      write(6,*) '  Partial rdf''s'
      write(6,*) '    [1] partial rdf for '//elem(1,1)//
     *        '       [2] partial rdf for '//elem(1,2)
      write(6,*) '    [3] partial rdf for '//elem(2,1)//
     *        '       [4] partial rdf for '//elem(2,2)
      write(6,*) '   [11] partial rdf for '//elem(1,1)//'-'//elem(1,1)//
     *           '   [12] partial rdf for '//elem(1,1)//'-'//elem(1,2)
      write(6,*) '   [13] partial rdf for '//elem(1,1)//'-'//elem(2,1)//
     *           '   [14] partial rdf for '//elem(1,1)//'-'//elem(2,2)
      write(6,*) '   [22] partial rdf for '//elem(1,2)//'-'//elem(1,2)//
     *           '   [23] partial rdf for '//elem(1,2)//'-'//elem(2,1)
      write(6,*) '   [24] partial rdf for '//elem(1,2)//'-'//elem(2,2)//
     *           '   [33] partial rdf for '//elem(2,1)//'-'//elem(2,1)
      write(6,*) '   [34] partial rdf for '//elem(2,1)//'-'//elem(2,2)//
     *           '   [44] partial rdf for '//elem(2,2)//'-'//elem(2,2)
      write(6,*) '(Enter p/r/one or two digit number)'
      read(5,'(a)') ans
      if(ans.eq.'R' .or. ans.eq.'r') then
         fcn='rdf'
      else if(ans.eq.'P' .or. ans.eq.'p') then
         fcn='pdf'
      else if(ans.eq.'1' .or. ans.eq.'2' .or.
     *        ans.eq.'3' .or. ans.eq.'4' ) then
         fcn='par'
         read(ans,'(i1)') i
         psl=(i+1)/2
         pel=mod(i,2)
         pel=2-pel
      else if(ans.eq.'11' .or. ans.eq.'12' .or.
     *        ans.eq.'13' .or. ans.eq.'14' .or.
     *        ans.eq.'22' .or. ans.eq.'23' .or.
     *        ans.eq.'24' .or. ans.eq.'33' .or.
     *        ans.eq.'34' .or. ans.eq.'44'  ) then
         fcn='ppr'
         read(ans,'(2i1)') i,j
         psl1=(i+1)/2
         pel1=mod(i,2)
         pel1=2-pel1
         psl2=(j+1)/2
         pel2=mod(j,2)
         pel2=2-pel2
      else
         goto 250
      end if
c     ------------------------------------------------------------------
c     output setup
c     ------------------------------------------------------------------
c350   write(6,*) 'enter output TopDrawer file name'
350   write(6,*) 'enter output file name'
      read(5,'(a20)') outfile
      wu=51
      open(unit=wu,file=outfile,status='new',err=350)
      if(fcn.eq.'pdf' .or. fcn.eq.'rdf') then
         write(wu,20) fcn
     *               ,elem(1,1),1-x,elem(1,2),x
     *               ,elem(2,1),1-y,elem(2,2),y, int(temp)
         write(wu,21)
20       format('#S 1 "',a3,' for ',
     *          2(a2,f3.1,a2,f3.1),
     *          ' ( T=',i4,' )"')
21       format('#L r G(r) dr dG')
c20       format('@    subtitle "',a3,' for ',
c     *          2(a2,'\s',f3.1,'\N',a2,'\s',f3.1,'\N'),
c     *          ' ( T=',i4,' )"')
c21       format('@    yaxis  label "G(r)"'
c     *         /1x,'@    xaxis  label "r"')
c         write(wu,20) fcn
c     *               ,elem(1,1),1-x,elem(1,2),x
c     *               ,elem(2,1),1-y,elem(2,2),y, int(temp)
c         write(wu,21) fcn
c20       format(1x,'title top ''',a3,' for '
c     *         ,2(a2,'0',f3.1,'1',a2,'0',f3.1,'1')
c     *         ,' ( T=',i4,' )''',  /1x,'case      ''',10x
c     *         ,'X   X  X   X  X   X  X   X''')
c21       format(1x,'title left ''',a3,''''
c     *         /1x,'title bottom ''r'''
c     *         /1x,'set limit x 0 20' )
      else if(fcn.eq.'par') then
         write(wu,30) elem(psl,pel)
     *               ,elem(1,1),1-x,elem(1,2),x
     *               ,elem(2,1),1-y,elem(2,2),y, int(temp)
         write(wu,31) elem(psl,pel)
30       format('# title top ''partial pdf for ',a2,' in '
     *         ,2(a2,'0',f3.1,'1',a2,'0',f3.1,'1')
     *         ,' ( T=',i4,' )''',  /,'# case      ''',24x
     *         ,'X   X  X   X  X   X  X   X''')
31       format('# title left ''partial g(r) for ',a2,''''
     *         /,'# title bottom ''r'''
     *         /,'# set limit x 0 20' )
      else if(fcn.eq.'ppr') then
         write(wu,40) elem(psl1,pel1),elem(psl2,pel2)
     *               ,elem(1,1),1-x,elem(1,2),x
     *               ,elem(2,1),1-y,elem(2,2),y, int(temp)
         write(wu,41) elem(psl1,pel1),elem(psl2,pel2)
40       format('# title top ''partial pdf for ',a2,'-',a2,' in '
     *         ,2(a2,'0',f3.1,'1',a2,'0',f3.1,'1')
     *         ,' ( T=',i4,' )''',  /,'# case      ''',24x
     *         ,'X   X  X   X  X   X  X   X''')
41       format('# title left ''partial g(r) for ',a2,'-',a2,''''
     *         /,'# title bottom ''r'''
     *         /,'# set limit x 0 20' )
      end if

c-----------------------------------------------------------------------
c     calculate pdf
c-----------------------------------------------------------------------
c      le= ((1-x)*radius(1,1)+x*radius(1,2))
c     *   +((1-y)*radius(2,1)+y*radius(2,2))
cccccccccccccccccccccccccccccccccccccccccccccccc
c following line is temporary for Aug. 1997
c
c      le=2.54416
cccccccccccccccccccccccccccccccccccccccccccccccc
      a0=le/sqrt(3.d0)*4
c     maximum length for output is half the supercell size
c     ----------------------------------------------------
      rmax=lx
      if(ly.lt.lx) rmax=ly
      if(lz.lt.lz) rmax=lz
      rmax=rmax*(a0/2)
c      rmax=rmax*(a0/2)*sqrt(3.d0)
c      rmax=10

      del=rmax/nout
c     average number density rho=(N/V)
c     --------------------------------
      rho=8/(a0*a0*a0)
c     initialize pdf
c     --------------
      do i=1,nout
         pdf(i)=0
         r(i)=del*i
      end do
c     accumulate pdf by summing over many gaussians
c     ---------------------------------------------
      do k=1,nfile
         ru=10+k
         nerr=0
         do i=1,nin
            read(ru,*,err=600,end=800) ibl,isl,iel,jbl,jsl,jel,r0,sigma
c        ...following two if's test for partials specified by p
c        ...for 'par', if( not ( (i=p) or (j=p) ) ), then skip
            if(fcn.eq.'par') then
               if(.not. ( (isl.eq.psl .and. iel.eq.pel) .or.
     *                    (jsl.eq.psl .and. jel.eq.pel) ) ) goto 500
            end if
c        ...for 'ppr', if(not ( (i=p1 and j=p2) or (i=p2 and j=p1) )),
c                      then skip
            if(fcn.eq.'ppr') then
               if(.not.(   ( (isl.eq.psl1 .and. iel.eq.pel1) .and.
     *                       (jsl.eq.psl2 .and. jel.eq.pel2) )
     *                 .or.( (isl.eq.psl2 .and. iel.eq.pel2) .and.
     *                       (jsl.eq.psl1 .and. jel.eq.pel1) ) 
     *                 )
     *           ) goto 500
            end if
            sigma2=sigma*sigma
            sigma10=sigma*10
            zz=scfct(isl,iel)*scfct(jsl,jel)
            do j=1,nout
               dr=r(j)-r0
               if(dr.le.sigma10) then
                  tt=zz*exp( -dr*dr/(2*sigma2) )/( sqrt(pix2*sigma2) )
c                  tt=exp( -dr*dr/(2*sigma) )/( sqrt(pix2*sigma) )
                  pdf(j)=pdf(j)+tt
               end if
            end do
500         continue
         end do
         stop '***** too many input lines, increase ''nin'''
c        reading error handling
c        ----------------------
600      nerr=nerr+1
c         goto 500
800      close(ru)
         if(nerr.ne.0) write(6,70) nerr,infile(k)
      end do
70    format(/// ' *******',i8,' errors in reading, lines ignored'///)
c     rescale pdf or rdf and write
c     ----------------------------
      do i=1,nout
c     ...a specific type of neighbor is written for every atom
c     ...a pair is counted only once in allsig.f
c     ...therefore pdf should be divided by natom/2 to get right weight
         pdf(i)=pdf(i)/nfile/(natom/2)
c     ...for pdf or partials, rescale
c         if(fcn.eq.'pdf' .or. fcn.eq.'par' .or. fcn.eq.'ppr')
c     ...for pdf rescale
         if(fcn.eq.'pdf')
     *      pdf(i)=(pdf(i)-pix4*rho*r(i)*r(i))/r(i)
         write(wu,'(4f15.8)') r(i),pdf(i),0.0,1.0
      end do
c      write(wu,*) 'plot'
c      write(wu,*) 'join'
      stop
      end

      subroutine nfact(el,scfact)
C-----------------------------------------------------------------------
C     Neutron scattering length for elements
C-----------------------------------------------------------------------
      implicit none
      character el*2
      real*8 scfact
         if     (el.eq.'In') then
            scfact=4.065
         else if(el.eq.'Ga') then
            scfact=7.288
         else if(el.eq.'As') then
            scfact=6.58
         else if(el.eq.'Sb') then
            scfact=5.57
         else if(el.eq.'Zn') then
            scfact=5.680
         else if(el.eq.'Se') then
            scfact=7.970
         else if(el.eq.'Te') then
            scfact=5.80
         else if(el.eq.'Cd') then
            scfact=4.87
         else
160         write(6,'('' Neutron scattering length for '',a2)') el
            read(5,*,err=160) scfact
         end if
      return
      end

      subroutine xfact(el,scfact)
C-----------------------------------------------------------------------
C     X-ray form factor for elements
C     for now, it is just the charge
C-----------------------------------------------------------------------
      implicit none
      character el*2
      real*8 scfact
         if     (el.eq.'In') then
            scfact=49
         else if(el.eq.'Ga') then
            scfact=31
         else if(el.eq.'As') then
            scfact=33
         else if(el.eq.'Sb') then
            scfact=51
         else if(el.eq.'Si') then
            scfact=14
         else if(el.eq.'Ge') then
            scfact=32
         else
160         write(6,'('' X-ray form factor for '',a2)') el
            read(5,*,err=160) scfact
         end if
      return
      end
      subroutine readhead(unit,elem,x,y,radius,mass,lx,ly,lz,temp,le,a0)
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
     *        le,a0,x0,y0,temp0,wmax0,xmax0,tmax0
      character elem(2,2)*2,elem0(2,2)*2

      SAVE elem0,lx0,ly0,lz0
      SAVE radius0,mass0,alpha0,beta0,lambda0
      SAVE x0,y0,temp0,wmax0,xmax0,tmax0

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
      read(unit,*) le,a0
      read(unit,*)

      if(unit.eq.11) then
c        if it is the first call, save current values as original values
c-----------------------------------------------------------------------
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
         temp0=temp
         wmax0=wmax
         xmax0=xmax
         tmax0=tmax
      else
c        if it is not a first call, compare with original values
c-----------------------------------------------------------------------
10       format(' *** wrong elements : ',a2,1x,a2)
20       format(' *** wrong',a7,' : ',i9,1x,i9)
30       format(' *** wrong',a7,' : ',f9.4,1x,f9.4)
         flag=0
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
         if(temp0.ne.temp) then
            write(6,30) 'temp',temp0,temp
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
      if(flag.ne.0) stop
      end if

      return
      end
