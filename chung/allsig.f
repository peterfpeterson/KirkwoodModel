      program allsig
cccc--------------------------------------------------------------------
c 
c    Purpose and Methods :
c       Calculates sigma (width of peaks in pair distribution fcn)
c       for given neighbor of an alloy of A_{1-x}B_{x}C
c       in ZnS structure
c       Using supercell of NxNxN cubic cells and integrating over
c       k=0 point in reduced zone scheme
c       ( sigma = < [ (u_2 - u_1).R_12 ]^2 >^{1/2} )
c
c    modified jan 1999:
c       To accomodate non-Vegard behavior in the case of different
c       alpha's and beta's, the supercell size via le is varied.
c       Set control to 1 in input file to activate this feature,
c       and to 0 if not.
c
c    Sample input file: input file should be one of the following 4
c
c       1) check [nn|m|stress|eig] infile outfile
c          123     1.d-10  1.d-8     :iseed ftol gtol
c          0                         :emethod
c             nn  : check nearest neighbor assignments
c             m   : check the validity of the dynamical matrix 
c                   when no disorder (x=y=0)
c             str : stress field
c             read parameters from 'infile'
c             write results on 'outfile'
c
c       2) static  infile  ufile
c          123     1.d-10  1.d-8     :iseed ftol gtol
c             read parameters from 'infile'
c             write static displacements u on 'ufile'
c             iseed : random number seed to generate a random alloy
c             ftol,gtol
c
c       3) eigen   ufile  efile
c          0                         :emethod
c             read parameters from 'infile', solve eigenvalue problem
c             write e-values & e-vectors on 'efile'
c             emethod: method of solving eig prob, 0-num rec, 1-imsl
c
c       4) thermal ufile  efile  outfile
c          300                       :temp
c             read static displacements from 'ufile' 
c                  e-values & e-vectors from 'efile'
c             write sigma at temp on outfile
c             temp : temp in K
c
c    infile contains following parameters:
c
c       elem(i,j)  : name of the elements in the i-th sublattice 
c                    with j-th element
c       x,y        : concentration in A_1-x B_x C_1-y D_y
c       lx,ly,lz   : number of cubic unit cells in each direction
c       radius(i,j): radius of the elements in the i-th sublattice
c                    with j-th element
c       mass(i,j)  : masses of the elements in the i-th sublattice
c                    with j-th element
c       alpha(i,j,k): bond stretching force constant 
c       beta(i,j,k,l): bond bending force constant
c       lambda     : force parameters in lambda-model
c
c    Sample infile for 'check' and 'static':
c
c GaInAsSb                           elem * 4
c 0.5   1.                           x,y
c 4  4  4                            Lx Ly Lz
c 1.2124  1.3878  1.2355   1.4135    radius * 4 (in angstrom)
c 69.735  114.82  74.9216  121.75    mass *4    (in amu)
c 105.3  105.3  105.3 105.3       alpha from A sublat (11)(12)(21)(22)
c 105.3  105.3  105.3 105.3       alpha from B sublat (11)(12)(21)(22)
c 17.28  17.28  17.28  17.28      beta from A sublat (111)(112)(121)(122)
c 17.28  17.28  17.28  17.28      beta from A sublat (211)(212)(221)(222)
c 17.28  17.28  17.28  17.28      beta from B sublat (111)(112)(121)(122)
c 17.28  17.28  17.28  17.28      beta from b sublat (211)(212)(221)(222)
c 1.                              lambda (1 for Kirkwood, 0 for Keating)
c
c 
c    Outputs
c       xaa,yaa,zaa : coordinates of neighbors (in fcc)
c       sigma : width of peaks in pair distribution fcn
c
c
c-----------------------------------------------------------------------
c       label of atoms in each sublattice:
c          Cubic unit cells are swept in the order of x, y and z.
c          In a cubic unit cell following labeling scheme is used
c             1 - (0 0 0)  2 - (1 1 0)
c             3 - (1 0 1)  4 - (0 1 1)
c          for A sublattice.  B sublattice with the same label is 
c          shifted from A sublattice by (1/2,1/2,1/2).
c          For example, if l=2, the top view of A sublattice looks like;
c             1st layer           2nd layer
c               12 10   16 14       28 26   32 30
c                9 11   13 15       25 27   29 31
c                4  2    8  6       20 18   24 22
c                1  3    5  7       17 19   21 23
c 
c       label of nearest neighbors in the array nn:
c          The direction of nn's are as follows
c          For A sublattice:            For B sublattice:
c             nna=1 - ( 1/2  1/2  1/2)     nnb=1 - (-1/2 -1/2 -1/2)
c             nna=2 - (-1/2 -1/2  1/2)     nnb=2 - ( 1/2  1/2 -1/2)
c             nna=3 - (-1/2  1/2 -1/2)     nnb=3 - ( 1/2 -1/2  1/2)
c             nna=4 - ( 1/2 -1/2 -1/2)     nnb=4 - (-1/2  1/2  1/2)
c          So, the first atom in B sublattice has nn's having the atom #
c          in A sublattice shown above as its nn #.
c          Also, nna=i and nnb=i have opposite directions.
c 
c    Require
c       Numerical Recipes (tred2, tqli, pythag, ran1 or ran3)
c 
c 
c    Created  03-JAN-1996   Jean Soo Chung
c 
cccc--------------------------------------------------------------------
      implicit none
      include 'allsig.inc'
      integer i,j
c===================== variables defined in allsig.inc =================
c      integer maxl,maxn,maxd
c      parameter (maxl=4,maxn=maxl*maxl*maxl*4,maxd=maxn*2*3)
cc     define variables
cc----------------------------------------------------------------------
c      integer nn,ldnn,ldm,ldd,lx,ly,lz,na,nd,el,iseed,control
c      integer ru,wu,chknn,chkm,chkstr,chkeig,emethod
c      character job*80,elem*2,infile*40,outfile*40
c      real*8        m,d,unitv,phi,dd,sd
c      real*8        x,y,le,a0,u,radius,mass,mu, alpha,beta,lambda,cc
c      real*8        wmax,tmax,xmax,temp,sigma,ftol,gtol
c      common /ctrl/ job,chknn,chkm,chkstr,chkeig,emethod
c      common /file/ ru,wu,infile,outfile
c      common /parm/ ldnn,ldm,ldd,lx,ly,lz,na,nd,ftol,gtol
c      common /matr/ m(maxd,maxd),d(maxd,maxd)
c      common /vect/ dd(maxd),sd(maxd)
c      common /stru/ unitv(4,2,3),nn(maxn,2,4)
c      common /conf/ x,y,le,a0,phi(maxd),el(maxn,2),iseed,control
c      common /atom/ u(maxd),radius(2,2),mass(2,2),elem(2,2),mu
c      common /forc/ alpha(2,2,2),beta(2,2,2,2),lambda,cc
c      common /cnst/ wmax,tmax,xmax,temp
c      common /sigm/ sigma
c========================== end of allsig.inc ==========================
c-----------------------------------------------------------------------
c     initial setup
c-----------------------------------------------------------------------
c     ------------------------------------------------------------------
c     read in, setup output file, make everything dimensionless
c     and random assignment of atoms
c     ------------------------------------------------------------------
      call setup()
c     ------------------------------------------------------------------
c     find nearest neighbors and nearest neighbor vectors
c     ------------------------------------------------------------------
      ldm=maxd
      ldnn=maxn
      call findnn(lx,ly,lz,nn,ldnn,unitv,phi)
      if(chknn.ne.0) then
         call checknn(wu,nn,ldnn,unitv,lx,ly,lz,phi)
      end if
c-----------------------------------------------------------------------
c     check
c-----------------------------------------------------------------------
      if(job(1:5).eq.'check') then
c        check stress field
c        ------------------
         if(chkstr.ne.0) then
c----       construct stress field phi (in units of length)
            call makephi(na,el,nn,ldnn,le,radius,alpha,unitv,phi,ldm)
            call wrtfld(wu,phi,na,el,ldnn,elem,'stress field phi')
         end if
         if(chkm.ne.0 .or. chkeig.ne.0) then
c----       construct force matrix M (dimensionless)
            call makem(na,m,ldm,nn,ldnn,unitv,alpha,beta,lambda,el)
c           check m
c           -------
            if(chkm.ne.0) then
               call wrtmtx(wu,m,nd,ldm,'matrix M')
               call checkm(wu,na,m,ldm,phi)
c----          so far, m is checked for l=2,3,4
            end if
c           check eigenvalues of d
c           ----------------------
            if(chkeig.ne.0) then
               call maked(na,m,ldm,ldnn,el,mass,d)
               call findeig(emethod,m,d,ldm,nd,dd,sd)
               write(wu,'(/'' eigenvalues-phonon freq.'')')
               do i=1,nd
                  sd(i)=sqrt(dd(i))
                  sd(i)=dd(i)
               end do
               call sort(nd,sd)
               do i=1,nd-1,3
                  write(wu,'(3f9.5)') sd(i),sd(i+1),sd(i+2)
               end do
            end if
         end if
         stop
c-----------------------------------------------------------------------
c     find static strain field u (in units of length) of atoms
c     use conjugate gradient method
c     (num. rec. (2nd ed) routine and imsl routine failed)
c-----------------------------------------------------------------------
      else if(job(1:6).eq.'static') then
c----    construct stress field phi (in units of length)
         call makephi(na,el,nn,ldnn,le,radius,alpha,unitv,phi,ldm)
         call makem(na,m,ldm,nn,ldnn,unitv,alpha,beta,lambda,el)
         call findu()
c----    write in the form to be read when 'thermal'
ccc following line is added (jan99)
c         write(wu,'(1x,2f20.10,t45,''le  a0'')') le,a0
         write(wu,*) le,a0,'  le a0'

cccccc test for two energy routines
c         control=1
c         nd=nd+1
c         call enersubm(nd,u,x)
c         call enersubx(nd,u,y)
c         write(wu,*) x,y,x-y
cccccc

         do i=1,na
            j=(i-1)*6
            write(wu,'(i2,3f18.14)') el(i,1),u(j+1),u(j+2),u(j+3)
            write(wu,'(i2,3f18.14)') el(i,2),u(j+4),u(j+5),u(j+6)
         end do
c         call wrtfld(wu,u,na,el,ldnn,elem,'strain field u')
         stop
c-----------------------------------------------------------------------
c     solve the eigenvalue problems
c-----------------------------------------------------------------------
      else if(job(1:5).eq.'eigen') then
         call makem(na,m,ldm,nn,ldnn,unitv,alpha,beta,lambda,el)
         call maked(na,m,ldm,ldnn,el,mass,d)
         call findeig(emethod,m,d,ldm,nd,dd,sd)
         do i=1,nd
            write(wu) dd(i),(d(i,j),j=1,nd)
         end do
         stop
c-----------------------------------------------------------------------
c     from read-in static displacements, eigenvalues and eigenvetors
c     calculate thermal part, sigma
c-----------------------------------------------------------------------
      else if(job(1:7).eq.'thermal') then
ccc following line is added (jan99)
         write(wu,'(1x,2f20.16,t45,''le  a0'')') le,a0
         call findsig(wu)
         stop
      end if
      end

************************************************************************
************************** end of main routine *************************
************************************************************************


      subroutine findsig(iunit)
cccc--------------------------------------------------------------------
c     find sigma
c-----------------------------------------------------------------------
      implicit none
      include 'allsig.inc'
      integer iunit,ll,i,j,
     *        ia,ie,ix,iy,iz,ibl,isl,
     *        ja,je,jx,jy,jz,jbl,jsl,
     *        kx,ky,kz
      real*8  ra(4,3),rij(3),rmag,rmax,ws,ns,eri,erj,mi,mj,sig
     *       ,urij(3)
      ll=lx*ly
c     limit calculation to the diagonal of half minimum peroid box
c     a0 is the half lattice constant
      rmax=lx
      if(ly.lt.lx) rmax=ly
      if(lz.lt.lz) rmax=lz
      rmax=rmax*a0*sqrt(3.d0)
c      write(iunit,*)
      write(iunit,*) '    i       j      rij         sig'
c-----------------------------------------------------------------------
c     find vectors to 4 A sublattice sites (in unit of length)
c-----------------------------------------------------------------------
      do i=1,4
      do j=1,3
         ra(i,j)=(unitv(1,1,j)+unitv(i,2,j))*le
      end do
      end do
c-----------------------------------------------------------------------
c     sweep through the lattice for the 1st atom
c-----------------------------------------------------------------------
      do iz=0,lz-1
      do iy=0,ly-1
      do ix=0,lx-1
c        write to unit 6 how the program is doing
         write(6,'('' ix ='',i2,'' iy ='',i2,'' iz ='',i2)') ix,iy,iz
      do ibl=1,4
         ia=(ix+iy*lx+iz*ll)*4+ibl
      do isl=1,2
         ie=(ia-1)*6+(isl-1)*3
         mi=mass(isl,el(ia,isl))
c-----------------------------------------------------------------------
c     sweep through the lattice for the 2nd atom
c     while taking care of periodic boundary conditions
c-----------------------------------------------------------------------
      do jz=-lz/2,lz/2
         kz=iz+jz
         if(kz.ge.lz) kz=kz-lz
         if(kz.lt.0) kz=kz+lz
      do jy=-ly/2,ly/2
         ky=iy+jy
         if(ky.ge.ly) ky=ky-ly
         if(ky.lt.0) ky=ky+ly
      do jx=-lx/2,lx/2
         kx=ix+jx
         if(kx.ge.lx) kx=kx-lx
         if(kx.lt.0) kx=kx+lx
      do jbl=1,4
         ja=(kx+ky*lx+kz*ll)*4+jbl
      do jsl=1,2
         je=(ja-1)*6+(jsl-1)*3
         mj=mass(jsl,el(ja,jsl))
c        to count every pair only once, skip if i > j
c        --------------------------------------------
         if(ie.ge.je) go to 100
c        ---------------------------------------------------------------
c        find unit vector and magnitude (in units of length) 
c        of the displacement vector from i to j
c        ---------------------------------------------------------------
         i=1
         rij(1)=jx*2*a0 + ra(jbl,i)-ra(ibl,i)
     *          +(unitv(1,1,i)*(jsl-1)-unitv(1,1,i)*(isl-1))*le
c     *          +u(je+i)-u(ie+i)
         i=2
         rij(2)=jy*2*a0 + ra(jbl,i)-ra(ibl,i)
     *          +(unitv(1,1,i)*(jsl-1)-unitv(1,1,i)*(isl-1))*le
c     *          +u(je+i)-u(ie+i)
         i=3
         rij(3)=jz*2*a0 + ra(jbl,i)-ra(ibl,i)
     *          +(unitv(1,1,i)*(jsl-1)-unitv(1,1,i)*(isl-1))*le
c     *          +u(je+i)-u(ie+i)
         rmag=sqrt(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
c        if rij >~ l/2, skip
c        -------------------
         if(rmag.gt.rmax) go to 100
         if(rmag.eq.0.d0) go to 100
c         rij(1)=rij(1)/rmag
c         rij(2)=rij(2)/rmag
c         rij(3)=rij(3)/rmag
**********beg: to check rij projected along Rij
         urij(1)=rij(1)/rmag
         urij(2)=rij(2)/rmag
         urij(3)=rij(3)/rmag
         rij(1)=(rij(1)+u(je+1)-u(ie+1))
         rij(2)=(rij(2)+u(je+2)-u(ie+2))
         rij(3)=(rij(3)+u(je+3)-u(ie+3))
         rmag=sqrt(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
c         rmag=abs(rij(1)*urij(1)+rij(2)*urij(2)+rij(3)*urij(3))
c         rij(1)=urij(1)
c         rij(2)=urij(2)
c         rij(3)=urij(3)
**********end


c        ---------------------------------------------------------------
c        branch sum at k=0 to find sigma
c        ---------------------------------------------------------------

**********
         sig=0
         do i=1,nd
         if(dd(i).ne.0.d0) then
            ws=sqrt(dd(i))
            ns=1/(exp(ws/temp)-1)+.5d0
            eri=d(ie+1,i)*urij(1)+d(ie+2,i)*urij(2)+d(ie+3,i)*urij(3)
            erj=d(je+1,i)*urij(1)+d(je+2,i)*urij(2)+d(je+3,i)*urij(3)
            sig=sig+ns/ws*((eri*eri/mi+erj*erj/mj)*.5d0
     *                     -eri*erj/sqrt(mi*mj))
         end if
         end do
         sig=sqrt((nd/(nd-3))*sig)*xmax
***********

c        ---------------------------------------------------------------
c        write the result
c        ---------------------------------------------------------------
         write(iunit,'(2(i4,i2,i2),2f12.8)')
     *        ia,isl,el(ia,isl),ja,jsl,el(ja,jsl),rmag,sig
100      continue
      end do
      end do
      end do
      end do
      end do
      end do
      end do
      end do
      end do
      end do
      return
      end

      subroutine findeig(emethod,m,d,ldm,nd,dd,sd)
cccc--------------------------------------------------------------------
c     solve eigen value problem
c       emethod: 0-Num.Rec.routine, 1-IMSL routine
c       d : matrix (input) eigenvectors (output)
c       m : matrix (workspace for IMSL)
c       ldm: leading dimension of d and m
c       nd : dimension of the problem
c       dd : vector of length ldm containing eigenvalues (output)
c       sd : vector of length ldm (workspace)
c
c       call tred2(a,n,np,d,e)
c          Householder reduction of a real, symmetric, n by n matrix a,
c          stored in an np by np physical array. On output, a is 
c          replaced by the orthogonal matrix Q effecting the 
c          transformation.  d returns the diagonal elements of the
c          tridiagonal matrix, and e the off-diagonal elements with
c          e(1)=0.
c       call tqli(d,e,n,np,z)
c          QL algorithm with implicit shifts, to determine the 
c          eigenvalues and eigenvectors of a real symmetric tridiagonal matrix,
c          or real symmetric matrix previously reduced by tred2.
c          d is a vector of length np.  On input, its first n elements 
c          are the diagonal elements of the tridiagonal matrix.
c          On output, it returns the eigenvalues.  The vector e inputs
c          the subdiagonal elements of the tridiagonal matrix, with
c          e(1) arbitrary.  On output e is destroyed.
c          If the eigenvctors of a tridiagonal matrix are desired, the 
c          matrix z (n by n matrix strored in np by np array) is input
c          as the identity matrix.  If the eigenvectors of a matrix 
c          that has been reduced by tred2 are required, z is input as
c          the matrix output by tred2.  In either case, the kth column 
c          of z returns the normalized eigenvector corresponding to d(k)
c       IMSL
c       call devcsf(n,a,lda,eval,evec,ldevec)
c          n - order of matrix A (input)
c          a - real symmetric matrix of order n (input)
c          eval - real vector of length n of eigenvalues (output)
c          evec - real matrix of order n, j-th eigenvector of eval(j)
c                 is stored in the j-th column (output)
c          ldevec - leading dimension of evec in calling program (input)
c-----------------------------------------------------------------------
      implicit none
      integer emethod,ldm,nd, i,j
      real*8  m(ldm,ldm),d(ldm,ldm),dd(ldm),sd(ldm),eps

      if(emethod.eq.0) write(6,'(//'' *** using num rec routine''//)')

ccc      if(emethod.eq.0) then
c        --------------------------------------------------
c        num. rec routine
c        dd retruns eigenvalue, d returns eigenvector
c          m   : matrix declared as ldm by ldm
c          n   : dimension of the problem
c          ldm : dimension of m in declaration
c          dd  : diagonal element of tridiagonal matrix from tred2
c          sd  : subdiagonal element of tridiagonal matrix from tred2
c          m   : on output of tqli,
c                k-th column contains the k-th eigenvectors
c        --------------------------------------------------
         call tred2(d,nd,ldm,dd,sd)
         call tqli(dd,sd,nd,ldm,d)
ccc      else
c        --------------------------------------------------
c        imsl routine
c        dd retruns eigenvalue, m returns eigenvector
c        --------------------------------------------------
ccc         call devcsf(nd,d,ldm,dd,m,ldm)
ccc         do i=1,nd
ccc         do j=1,nd
ccc            d(i,j)=m(i,j)
ccc         end do
ccc         end do
ccc      end if
c     if eigenvalue is too small (numerical error), set it to zero
c     ------------------------------------------------------------------
      eps=1.d-12
      do i=1,nd
         if(abs(dd(i)).lt.eps) dd(i)=0.d0
      end do
      return
      end

      subroutine makephi(na,el,nn,ldnn,le,radius,alpha,unitv,phi,ldm)
cccc--------------------------------------------------------------------
c     construct stress field phi (in units of length)
c-----------------------------------------------------------------------
      implicit none
      integer na,ldnn,ldm, el(ldnn,2),nn(ldnn,2,4),
     *        ibl,isl,jsl,eli,elj,i,j,i1,j1
      real*8  le,radius(2,2),alpha(2,2,2),unitv(4,2,3),phi(ldm),lij
      do ibl=1,na
      do isl=1,2
         i=(ibl-1)*6+(isl-1)*3
         jsl=mod(isl,2)+1
         eli=el(ibl,isl)
         do i1=1,3
            j=i+i1
            phi(j)=0
            do j1=1,4
               elj=el(nn(ibl,isl,j1),jsl)
               lij=le-radius(isl,eli)-radius(jsl,elj)
               phi(j)=phi(j)+alpha(isl,eli,elj)
     *                      *lij*unitv(j1,isl,i1)
            end do
         end do
      end do
      end do
      return
      end

      subroutine checkm(wu,na,m,ldm,phi)
cccc--------------------------------------------------------------------
c     numerical differentiation of energy to check matrix m
c     ********* for now, no disorder only (Lij=Le) *********
c-----------------------------------------------------------------------
      implicit none
      integer wu,na,ldm,nd, i,j
      real*8 m(ldm,ldm),phi(na*6),mij,tol,eps,ep1,ep2,e1,e2,e3,e4,
     *       energy1
      nd=na*6
      eps=0.0001
      tol=0.00001
      ep1=eps*eps
      ep2=4*ep1
      do i=1,nd
        phi(i)=0
      end do
      do i=1,nd
      do j=1,nd
c        diagonl block
c        -------------
         if(i.eq.j) then
            phi(i)=eps
            e1=energy1()
            phi(i)=-eps
            e2=energy1()
            mij=(e1+e2)/ep1
c        off diagonl block
c        -----------------
         else
            phi(i)=eps
            phi(j)=eps
            e1=energy1()
            phi(i)=-eps
            e2=energy1()
            phi(j)=-eps
            e3=energy1()
            phi(j)=eps
            e4=energy1()
            phi(i)=0
            phi(j)=0
            mij=((e1-e2)-(e4-e3))/ep2
         end if
         if(abs(mij-m(i,j)).gt.eps)
     *      write(wu,10) i,j,m(i,j),mij
10    format(' *******',2i5, 2f10.5)
      end do
      end do
      end

      function energy1()
cccc--------------------------------------------------------------------
c     calculates energy associated with one atom due to u(0) and u(i)
c     for numerical construction of matrix m
c     ********* for now, no disorder only (Lij=Le) *********
c-----------------------------------------------------------------------
      implicit none
      include 'allsig.inc'
      integer ii,jj,ll,ibl,isl,jsl,jn,jnn,ln,lnn
      real*8 lam3,uij(3),uil(3),rijuij,rijuil,riluij,riluil,fijl,
     *       term1,term2,energy1
      lam3=lambda/3
      term1=0
      term2=0
      do ibl=1,na
      do isl=1,2
         jsl=mod(isl,2)+1
         ii=(ibl-1)*6+(isl-1)*3
         do jn=1,4
            jnn=nn(ibl,isl,jn)
            jj=(jnn-1)*6+(jsl-1)*3
            uij(1)=phi(jj+1)-phi(ii+1)
            uij(2)=phi(jj+2)-phi(ii+2)
            uij(3)=phi(jj+3)-phi(ii+3)
            rijuij=unitv(jn,isl,1)*uij(1)
     *            +unitv(jn,isl,2)*uij(2)
     *            +unitv(jn,isl,3)*uij(3)
            term1=term1+rijuij*rijuij
            do ln=jn+1,4
               lnn=nn(ibl,isl,ln)
               ll=(lnn-1)*6+(jsl-1)*3
               uil(1)=phi(ll+1)-phi(ii+1)
               uil(2)=phi(ll+2)-phi(ii+2)
               uil(3)=phi(ll+3)-phi(ii+3)
               rijuil=unitv(jn,isl,1)*uil(1)
     *               +unitv(jn,isl,2)*uil(2)
     *               +unitv(jn,isl,3)*uil(3)
               riluij=unitv(ln,isl,1)*uij(1)
     *               +unitv(ln,isl,2)*uij(2)
     *               +unitv(ln,isl,3)*uij(3)
               riluil=unitv(ln,isl,1)*uil(1)
     *               +unitv(ln,isl,2)*uil(2)
     *               +unitv(ln,isl,3)*uil(3)
               fijl=rijuil+riluij+lam3*(rijuij+riluil)
               term2=term2+fijl*fijl
            end do
         end do
      end do
      end do
      energy1=alpha(1,1,1)*term1/4+beta(1,1,1,1)*term2/8
      end

      subroutine makem(na,m,ldm,nn,ldnn,unitv,alpha,beta,lambda,el)
cccc--------------------------------------------------------------------
c     compose the force matrix M
c        we need 6 indices of M for one fcc unit cell
c        first 3 are used for 3 components of A sublattice
c        next  3 are used for 3 components of B sublattice
c-----------------------------------------------------------------------
      implicit none
      integer na,nd,ldm,ldnn,nn(ldnn,2,4),el(ldnn,2),eli,elj,ell,
     *        i,j,ibl,jbl,isl,jsl,jnn,lnn,mnn,al,be
      real*8 m(ldm,ldm),unitv(4,2,3),alpha(2,2,2),beta(2,2,2,2),
     *       lambda,lam3,
     *       term1,term2,term3
      lam3=lambda/3
      nd=na*6
      do i=1,nd
      do j=1,nd
         m(i,j)=0
      end do
      end do
c-----------------------------------------------------------------------
c     diagonal block
c-----------------------------------------------------------------------
c     sweep through Bravais Lattice and SubLattices
      do ibl=1,na
      do isl=1,2
c                             Bravais lattice no.-1
         i=(ibl-1)*6+(isl-1)*3
c                        the other sublattice
         jsl=mod(isl,2)+1
         eli=el(ibl,isl)
c        all alpha and beta components
         do al=1,3
         do be=1,3
            term1=0
            term2=0
            term3=0
            do jnn=1,4
ccc following line added (jan.99)
               jbl=nn(ibl,isl,jnn)
               elj=el(jbl,jsl)
               term1=term1+alpha(isl,eli,elj)
     *                    *unitv(jnn,isl,al)*unitv(jnn,isl,be)
               do mnn=jnn+1,jnn+3
                  lnn=mod(mnn-1,4)+1
                  ell=el(nn(ibl,isl,lnn),jsl)
                  term2=term2+beta(isl,eli,elj,ell)
     *                    *(unitv(jnn,isl,al)+unitv(lnn,isl,al))
     *                    *(unitv(jnn,isl,be)+unitv(lnn,isl,be))
ccc following line added (jan.99)
                  ell=el(nn(jbl,jsl,lnn),isl)
                  term3=term3+beta(jsl,elj,eli,ell)
     *                    *(unitv(lnn,isl,al)+lam3*unitv(jnn,isl,al))
     *                    *(unitv(lnn,isl,be)+lam3*unitv(jnn,isl,be))
               end do
            end do
            m(i+al,i+be)=m(i+al,i+be)+term1+(1+lam3)**2*term2/8+term3/4
         end do
         end do
      end do
      end do
c-----------------------------------------------------------------------
c     nn block
c-----------------------------------------------------------------------
      do ibl=1,na
      do isl=1,2
c                             Bravais lattice no.-1
         i=(ibl-1)*6+(isl-1)*3
c                        the other sublattice
         jsl=mod(isl,2)+1
         eli=el(ibl,isl)
         do jnn=1,4
            jbl=nn(ibl,isl,jnn)
            j=(jbl-1)*6+(jsl-1)*3
            elj=el(jbl,jsl)
         do al=1,3
         do be=1,3
            term1=alpha(isl,eli,elj)
     *            *unitv(jnn,isl,al)*unitv(jnn,isl,be)
            term2=0
            term3=0
            do mnn=jnn+1,jnn+3
               lnn=mod(mnn-1,4)+1
               ell=el(nn(ibl,isl,lnn),jsl)
               term2=term2+beta(isl,eli,elj,ell)
     *               *(unitv(jnn,isl,al)+unitv(lnn,isl,al))
     *               *(lam3*unitv(jnn,isl,be)+unitv(lnn,isl,be))
ccc following line modified (jan.99)
c               ell=el(nn(ibl,isl,lnn),jsl)
               ell=el(nn(jbl,jsl,lnn),isl)
               term3=term3+beta(jsl,elj,eli,ell)
     *               *(unitv(lnn,isl,al)+lam3*unitv(jnn,isl,al))
     *               *(unitv(lnn,isl,be)+unitv(jnn,isl,be))
            end do
            m(i+al,j+be)
     *         =m(i+al,j+be)-term1-(1+lam3)*(term2+term3)/4
         end do
         end do
         end do
      end do
      end do
c-----------------------------------------------------------------------
c     nnn block
c-----------------------------------------------------------------------
      do ibl=1,na
      do isl=1,2
c                             Bravais lattice no.-1
         i=(ibl-1)*6+(isl-1)*3
c                        the other sublattice
         jsl=mod(isl,2)+1
         eli=el(ibl,isl)
         do jnn=1,4
            jbl=nn(ibl,isl,jnn)
            elj=el(jbl,jsl)
         do mnn=jnn+1,jnn+3
            lnn=mod(mnn-1,4)+1
            j=(nn(jbl,jsl,lnn)-1)*6+(isl-1)*3
            ell=el(nn(jbl,jsl,lnn),isl)
         do al=1,3
         do be=1,3
            term3=beta(jsl,elj,eli,ell)
     *            *(unitv(lnn,isl,al)+lam3*unitv(jnn,isl,al))
     *            *(unitv(jnn,isl,be)+lam3*unitv(lnn,isl,be))
            m(i+al,j+be)=m(i+al,j+be)+term3/4
         end do
         end do
         end do
         end do
      end do
      end do

      return
      end

      subroutine maked(na,m,ldm,ldnn,el,mass,d)
cccc--------------------------------------------------------------------
c     compose the dynamical matrix D from the force matrix M
c
c                          -1/2                [-ik.(x(ir)-x(js)]
c        D(rs|k) = (m_r m_s)    sum M (ir;js) e
c                                j   ab
c                  i,j : bravais lattice index
c                  r,s : basis index
c                  a,b : axis index (x,y,z)
c
c     summation has already been done in getting M
c     so, we only need to multiply it by proper mass
c-----------------------------------------------------------------------
      implicit none
      integer na,nd,ldm,ldnn,el(ldnn,2)
      integer  i,j,isl,jsl,i1,j1,al,be
      real*8  m(ldm,ldm),d(ldm,ldm),mass(2,2),mm
c-----------------------------------------------------------------------
      nd=na*6
      do i=1,na
      do isl=1,2
         i1=(i-1)*6+(isl-1)*3
      do j=1,na
      do jsl=1,2
         j1=(j-1)*6+(jsl-1)*3
         mm=mass(isl,el(i,isl))*mass(jsl,el(j,jsl))
         do al=1,3
         do be=1,3
            d(i1+al,j1+be)=m(i1+al,j1+be)/sqrt(mm)
         end do
         end do
      end do
      end do
      end do
      end do
      return
      end

      subroutine findu()
cccc--------------------------------------------------------------------
c     solve linear system of equations and
c     find the equilibrium position of strained lattice
c       next two method have been tried but faild when L>2
c       num. rec. routines (singular value decomposition and back sub.)
c          dsvdcmp(a,m,n,mp,np,w,v)
c          dsvbksb(u,w,v,m,n,mp,np,b,x)
c       imsl routine
c          call dlslsf(n,a,lda,b,x)
c     Therefore, only conjugate gradient method is coded here.
c     call conjugate gradient method
c     ------------------------------
c        call cgmin (n,p,g,ftol,gtol,mxitr,iter,fret,fchg,res,
c                    enersub,gradsub,monsub, wk,nwk,irprt,lumsg)
c-----------------------------------------------------------------------
      implicit none
      include 'allsig.inc'
ccc following line is modified (jan99)
c      integer i,maxiter,iter,nwk,ireport,lumsg
      integer i,nd1,maxiter,iter,nwk,ireport,lumsg
      real*8  fret,fchg,res,p(maxd)
c      external enersub,gradsub,monsub
      external enersubm,gradsubm, enersubx,gradsubx,monsub

ccc start of modification (jan99)
      nd1=nd
      if(control.eq.1) then
         nd1=nd+1
         p(nd1)=le
      end if
ccc end of modification (jan99)
      ftol=na*ftol
ccc nd is changed to nd1 in following 2 lines (jan99)
      maxiter=nd1*1.5
      nwk=6*nd1
      ireport=2
      lumsg=6
      do i=1,nd
         p(i)=phi(i)
      end do
ccc in the following line, nd is changed to nd1 (jan99)
c      call cgmin (nd1,p,dd,ftol,gtol,maxiter,iter,fret,fchg,res,
c     &            enersub,gradsub,monsub, d,nwk,ireport,lumsg)


      if(elem(1,1).eq.'mm') then
      call cgmin (nd1,p,dd,ftol,gtol,maxiter,iter,fret,fchg,res,
     &            enersubm,gradsubm,monsub, d,nwk,ireport,lumsg)
      else if(elem(1,1).eq.'mx') then
      call cgmin (nd1,p,dd,ftol,gtol,maxiter,iter,fret,fchg,res,
     &            enersubm,gradsubx,monsub, d,nwk,ireport,lumsg)
      else if(elem(1,1).eq.'xm') then
      call cgmin (nd1,p,dd,ftol,gtol,maxiter,iter,fret,fchg,res,
     &            enersubx,gradsubm,monsub, d,nwk,ireport,lumsg)
      else
      call cgmin (nd1,p,dd,ftol,gtol,maxiter,iter,fret,fchg,res,
     &            enersubx,gradsubx,monsub, d,nwk,ireport,lumsg)
      end if


      write(6,'('' final energy ='',e12.4)') fret
      do i=1,nd
         u(i)=p(i)
      end do
ccc start of modification (jan99)
      if(control.eq.1) then
         le=p(nd1)
         a0=le/sqrt(3.d0)*2
      end if
ccc end of modification (jan99)
      return
      end

      subroutine enersubm (n,p,f)
cccc--------------------------------------------------------------------
c     returns energy for cgmin
c     uses matrix notation form of energy
c     only one of this and the following should be named 'gradsub'
c-----------------------------------------------------------------------
      implicit none
      include 'allsig.inc'
      integer n,i,j
      real*8  p(n),f
ccc start of modification (jan99)
      real*8  e0,dl
      integer eli,elj

      f=0
      if(control.eq.1) then
         le=p(n)
         call makephi(na,el,nn,ldnn,le,radius,alpha,unitv,phi,ldm)
c        add otherwise-constant energy associated with le
         do i=1,na
         do j=1,4
            eli=el(i,1)
            elj=el(nn(i,1,j),2)
            dl=le-(radius(1,eli)+radius(2,elj))
            f=f+alpha(1,eli,elj)*dl*dl
         end do
         end do
      end if
ccc end of modification (jan99)
ccc      f=0
      do i=1,nd
      do j=1,nd
         f=f+p(i)*m(i,j)*p(j)
      end do
      end do
      f=f*0.5
      do i=1,nd
         f=f-p(i)*phi(i)
      end do
      return
      end

      subroutine enersubx (n,p,f)
cccc--------------------------------------------------------------------
c     returns energy for cgmin
c     uses direct expression for energy
c     copied from relax.f for comparison
c-----------------------------------------------------------------------
      implicit none
      include 'allsig.inc'
      integer n
      real*8  p(n),f
      integer ii,jj,ll,eli,elj,ell,ibl,isl,jsl,jn,jnn,ln,lnn
      real*8 uij(3),uil(3),rijuij,rijuil,riluij,riluil,
     *       lij,fijl,term1,term2,lam3
      lam3=lambda/3
      term1=0
      term2=0
ccc following line is a modification (jan99)
      if(control.eq.1) le=p(n)
c  ...sweep through the bravais lattice
      do ibl=1,na
c  ...for both sublattices
      do isl=1,2
         ii=(ibl-1)*6+(isl-1)*3
         eli=el(ibl,isl)
c     ...the other sublattice where nn's reside
         jsl=mod(isl,2)+1
c     ...sum over 1st nn
         do jn=1,4
            jnn=nn(ibl,isl,jn)
            jj=(jnn-1)*6+(jsl-1)*3
            elj=el(jnn,jsl)
            uij(1)=p(jj+1)-p(ii+1)
            uij(2)=p(jj+2)-p(ii+2)
            uij(3)=p(jj+3)-p(ii+3)
            rijuij=unitv(jn,isl,1)*uij(1)
     *            +unitv(jn,isl,2)*uij(2)
     *            +unitv(jn,isl,3)*uij(3)
c        ...sum bond-stretching only over the sublattice A
            if(isl.eq.1) then
               lij=le-radius(isl,eli)-radius(jsl,elj)+rijuij
               term1=term1+alpha(isl,eli,elj)*lij*lij
            end if
c        ...sum over 2nd nn to make an angle
            do ln=jn+1,4
               lnn=nn(ibl,isl,ln)
               ll=(lnn-1)*6+(jsl-1)*3
               ell=el(lnn,jsl)
               uil(1)=p(ll+1)-p(ii+1)
               uil(2)=p(ll+2)-p(ii+2)
               uil(3)=p(ll+3)-p(ii+3)
               rijuil=unitv(jn,isl,1)*uil(1)
     *               +unitv(jn,isl,2)*uil(2)
     *               +unitv(jn,isl,3)*uil(3)
               riluij=unitv(ln,isl,1)*uij(1)
     *               +unitv(ln,isl,2)*uij(2)
     *               +unitv(ln,isl,3)*uij(3)
               riluil=unitv(ln,isl,1)*uil(1)
     *               +unitv(ln,isl,2)*uil(2)
     *               +unitv(ln,isl,3)*uil(3)
               fijl=rijuil+riluij+lam3*(rijuij+riluil)
               term2=term2+beta(isl,eli,elj,ell)*fijl*fijl
            end do
         end do
      end do
      end do
      f=term1/2+term2/8

      return
      end

      subroutine gradsubm (n,p,g,gmag)
cccc--------------------------------------------------------------------
c     returns gradient for cgmin
c     uses matrix notation form of energy
c     only one of this and the following should be named 'gradsub'
c-----------------------------------------------------------------------
      implicit none
      include 'allsig.inc'
ccc following line is modified (jan99)
c     integer n,i,j
      integer n,i,j, ii,jj,ibl,isl,jsl,eli,elj,jnn,jn
ccc following line is modified (jan99)
c     real*8  p(n),g(n),gmag
      real*8  p(n),g(n),gmag, uij(3),rijuij,lij
ccc start of modification (jan99)
c  ...the last component of gradient (dV/dL_e) when control=1
      if(control.eq.1) then
         le=p(n)
         call makephi(na,el,nn,ldnn,le,radius,alpha,unitv,phi,ldm)
         g(n)=0
c     ...sweep through the bravais lattice (A sublattice only)
         isl=1
         jsl=2
         do ibl=1,na
            ii=(ibl-1)*6+(isl-1)*3
            eli=el(ibl,isl)
c        ...sum over nn
            do jn=1,4
               jnn=nn(ibl,isl,jn)
               jj=(jnn-1)*6+(jsl-1)*3
               elj=el(jnn,jsl)
               uij(1)=p(jj+1)-p(ii+1)
               uij(2)=p(jj+2)-p(ii+2)
               uij(3)=p(jj+3)-p(ii+3)
               rijuij=unitv(jn,isl,1)*uij(1)
     *               +unitv(jn,isl,2)*uij(2)
     *               +unitv(jn,isl,3)*uij(3)
c           ...sum bond-stretching only over the sublattice A
               if(isl.eq.1) then
                  lij=le-radius(isl,eli)-radius(jsl,elj)+rijuij
                  g(n)=g(n)+alpha(isl,eli,elj)*lij
               end if
            end do
         end do
      end if
ccc end of modification (jan99)
      gmag=0
      do i=1,nd
         g(i)=0
         do j=1,nd
            g(i)=g(i)+m(i,j)*p(j)
         end do
         g(i)=g(i)-phi(i)
         gmag=gmag+g(i)*g(i)
      end do
ccc following if-block added (jan99)
      if(control.eq.1) then
         gmag=gmag+g(n)*g(n)
      end if
      return
      end

      subroutine gradsubx (n,p,g,gmag)
cccc--------------------------------------------------------------------
c     returns gradient for cgmin - analytic expression
c     uses direct expression for energy
c     copied from relax.f for comparison
c-----------------------------------------------------------------------
      implicit none
      include 'allsig.inc'
      integer n
      real*8  p(n),g(n),gmag
      integer ii,jj,ll,k,ik,eli,elj,ell,ibl,isl,jsl,jn,jnn,ln,lnn,mn
      real*8 uij(3),uil(3),rijuij,rijuil,riluij,riluil,
     *       lij,fijl,lam3
      lam3=lambda/3
      do k=1,n
         g(k)=0
      end do
ccc start of modification (jan99)
c  ...the last component of gradient (dV/dL_e) when control=1
      if(control.eq.1) then
         le=p(n)
         call makephi(na,el,nn,ldnn,le,radius,alpha,unitv,phi,ldm)
c     ...sweep through the bravais lattice (A sublattice only)
         isl=1
         jsl=2
         do ibl=1,na
            ii=(ibl-1)*6+(isl-1)*3
            eli=el(ibl,isl)
c        ...sum over nn
            do jn=1,4
               jnn=nn(ibl,isl,jn)
               jj=(jnn-1)*6+(jsl-1)*3
               elj=el(jnn,jsl)
               uij(1)=p(jj+1)-p(ii+1)
               uij(2)=p(jj+2)-p(ii+2)
               uij(3)=p(jj+3)-p(ii+3)
               rijuij=unitv(jn,isl,1)*uij(1)
     *               +unitv(jn,isl,2)*uij(2)
     *               +unitv(jn,isl,3)*uij(3)
c           ...sum bond-stretching only over the sublattice A
               if(isl.eq.1) then
                  lij=le-radius(isl,eli)-radius(jsl,elj)+rijuij
                  g(n)=g(n)+alpha(isl,eli,elj)*lij
               end if
            end do
         end do
      end if
ccc end of modification (jan99)
c  ...sweep through all the atoms for dV/du_i
      do ibl=1,na
      do isl=1,2
         ii=(ibl-1)*6+(isl-1)*3
         eli=el(ibl,isl)
c     ...the other sublattice where nn's reside
         jsl=mod(isl,2)+1
c     ...sum over 1st nn
c        ---------------------------------------------------------------
         do jn=1,4
            jnn=nn(ibl,isl,jn)
            jj=(jnn-1)*6+(jsl-1)*3
            elj=el(jnn,jsl)
            uij(1)=p(jj+1)-p(ii+1)
            uij(2)=p(jj+2)-p(ii+2)
            uij(3)=p(jj+3)-p(ii+3)
            rijuij=unitv(jn,isl,1)*uij(1)
     *            +unitv(jn,isl,2)*uij(2)
     *            +unitv(jn,isl,3)*uij(3)
            lij=le-radius(isl,eli)-radius(jsl,elj)+rijuij
            do ik=1,3
               k=ii+ik
               g(k)=g(k)-alpha(isl,eli,elj)*lij*unitv(jn,isl,ik)
            end do
c        ...sum over 2nd nn to make angles around (ibl,isl)
c           ------------------------------------------------------------
            do ln=jn+1,4
               lnn=nn(ibl,isl,ln)
               ll=(lnn-1)*6+(jsl-1)*3
               ell=el(lnn,jsl)
               uil(1)=p(ll+1)-p(ii+1)
               uil(2)=p(ll+2)-p(ii+2)
               uil(3)=p(ll+3)-p(ii+3)
               rijuil=unitv(jn,isl,1)*uil(1)
     *               +unitv(jn,isl,2)*uil(2)
     *               +unitv(jn,isl,3)*uil(3)
               riluij=unitv(ln,isl,1)*uij(1)
     *               +unitv(ln,isl,2)*uij(2)
     *               +unitv(ln,isl,3)*uij(3)
               riluil=unitv(ln,isl,1)*uil(1)
     *               +unitv(ln,isl,2)*uil(2)
     *               +unitv(ln,isl,3)*uil(3)
               fijl=(rijuil+riluij+lam3*(rijuij+riluil))*(1+lam3)*0.25
               do ik=1,3
                  k=ii+ik
                  g(k)=g(k)-beta(isl,eli,elj,ell)*fijl
     *                  *(unitv(jn,isl,ik)+unitv(ln,isl,ik))
               end do
            end do
c        ...sum over nnn to make angles to nnn
c           ------------------------------------------------------------
            do mn=jn+1,jn+3
               ln=mod(mn-1,4)+1
               lnn=nn(jnn,jsl,ln)
               ll=(lnn-1)*6+(isl-1)*3
               ell=el(lnn,isl)
c           ...uil <= -ujl = ulj
               uil(1)=p(jj+1)-p(ll+1)
               uil(2)=p(jj+2)-p(ll+2)
               uil(3)=p(jj+3)-p(ll+3)
c           ...rijulj=rjiujl
               rijuil=unitv(jn,isl,1)*uil(1)
     *               +unitv(jn,isl,2)*uil(2)
     *               +unitv(jn,isl,3)*uil(3)
c           ...rljuij=rjluji
               riluij=unitv(ln,isl,1)*uij(1)
     *               +unitv(ln,isl,2)*uij(2)
     *               +unitv(ln,isl,3)*uij(3)
c           ...rljulj=rjlujl
               riluil=unitv(ln,isl,1)*uil(1)
     *               +unitv(ln,isl,2)*uil(2)
     *               +unitv(ln,isl,3)*uil(3)
               fijl=(rijuil+riluij+lam3*(rijuij+riluil))*0.25
               do ik=1,3
                  k=ii+ik
                  g(k)=g(k)-beta(jsl,elj,eli,ell)*fijl
     *                    *(unitv(ln,isl,ik)+lam3*unitv(jn,isl,ik))
               end do
            end do
         end do
      end do
      end do
      gmag=0
      do k=1,n
         gmag=gmag+g(k)*g(k)
      end do
      return
      end

      subroutine monsub (iter,p,g,fp,fchg,res)
cccc--------------------------------------------------------------------
c     cgmin monitoring routine - does nothing
c-----------------------------------------------------------------------
      implicit none
      integer iter
      real*8  p(1),g(1),fp,fchg,res
      integer       ldnn,ldm,ldd,lx,ly,lz,na,nd
      common /parm/ ldnn,ldm,ldd,lx,ly,lz,na,nd
      if(mod(iter*5,na).eq.0) write(6,10) iter,fp,fchg,res,p(1),g(1)
10    format(' iteration #',i7,'    energy =',e12.4,4e12.4)
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
c           nna=1 - ( 1/2  1/2  1/2)     unitv(1,1,i)
c           nna=2 - (-1/2 -1/2  1/2)     unitv(2,1,i)
c           nna=3 - (-1/2  1/2 -1/2)     unitv(3,1,i)
c           nna=4 - ( 1/2 -1/2 -1/2)     unitv(4,1,i)
c        For B sublattice:
c           nnb=1 - (-1/2 -1/2 -1/2)     unitv(1,2,i)
c           nnb=2 - ( 1/2  1/2 -1/2)     unitv(2,2,i)
c           nnb=3 - ( 1/2 -1/2  1/2)     unitv(3,2,i)
c           nnb=4 - (-1/2  1/2  1/2)     unitv(4,2,i)
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

      subroutine checknn(wu,nn,ldnn,unitv,lx,ly,lz,phi)
cccc--------------------------------------------------------------------
c     check nn distance to see if nn assignment is right
c-----------------------------------------------------------------------
      implicit none
      integer wu,ldnn,nn(ldnn,2,4),lx,ly,lz,na,
     *        ix,iy,iz,ix1,iy1,iz1,ibl,isl,jsl,ii,ll,mx,iflag,i,j
      real*8  unitv(4,2,3),phi(ldnn*6),dx,dy,dz,dr
      na=lx*ly*lz*4
      ix=lx*4
      iy=ly*4
      iz=lz*4
      ix1=lx*2
      iy1=ly*2
      iz1=lz*2
      do ibl=1,na
      do isl=1,2
         jsl=mod(isl,2)+1
         ii=(ibl-1)*6+(isl-1)*3
         do mx=1,4
            ll=(nn(ibl,isl,mx)-1)*6+(jsl-1)*3
            dx=phi(ii+1)-phi(ll+1)
            dy=phi(ii+2)-phi(ll+2)
            dz=phi(ii+3)-phi(ll+3)
c           apply periodic boundary conditions
c           ----------------------------------
            if(dx.lt.-ix1) dx=dx+ix
            if(dx.gt. ix1) dx=dx-ix
            if(dy.lt.-iy1) dy=dy+iy
            if(dy.gt. iy1) dy=dy-iy
            if(dz.lt.-iz1) dz=dz+iz
            if(dz.gt. iz1) dz=dz-iz
            dr=dx*dx+dy*dy+dz*dz
            if(dr.lt.2.999d0 .or. dr.gt.3.001) then
               iflag=1
               write(wu,10) ibl,isl,mx,nn(ibl,isl,mx),dr,dx,dy,dz
     c                      ,ii-(isl-1)*3,ll-(jsl-1)*3
            end if
10          format(' ***** ibl=',i4,'  isl=',i1,'  inn=',i1,
     c             '  nn=',i4,'    dr=',f6.2,'  dx,dy,dz=',3f5.2,2i5) 
         end do
      end do
      end do
      if(iflag.eq.1) stop '***** wrong nn'
      write(wu,*)
      write(wu,*) 'nearest neighbor vectors'
      do ix=1,4
         write(wu,'(i4,2(a3,1x,3f6.3,2x))')
     *         ix,'A',(unitv(i,1,j),j=1,3),'B',(unitv(i,2,j),j=1,3)
      end do
      write(wu,*)
      write(wu,*) 'nearest neighbor definitions'
      write(wu,*) ' ibl sl1  nn1 nn2 nn3 nn3   sl2  nn1 nn2 nn3 nn3'
      do i=1,na
         write(wu,'(i4,2(a4,1x,4i4,2x))')
     *         i,'A',(nn(i,1,j),j=1,4),'B',(nn(i,2,j),j=1,4)
      end do
      return
      end

      subroutine wrtmtx(iunit,m,n,ldm,name)
cccc--------------------------------------------------------------------
c     write a matrix on the output file
c-----------------------------------------------------------------------
      implicit none
      integer iunit,n,ldm,i,j,k,l
      real*8 m(ldm,ldm),csum(96)
      character name*(*)
      do i=1,96
        csum(i)=0
      end do
      if(n.le.24) then
         write(iunit,*)
         write(iunit,*) name
         write(iunit,'(4x,24(1x,i4))') (j,j=1,n)
         do i=1,n
            write(iunit,'(i4,24(1x,f4.1))') i,(m(i,j),j=1,n)
         end do
         do j=1,n
         do i=1,n
            csum(j)=csum(j)+m(i,j)
         end do
         end do
         write(iunit,'('' csm'',24(1x,f4.1))') (csum(j),j=1,n)
      else
         write(iunit,*)
         write(iunit,*) name
         do k=0,n/96-1
            write(iunit,'(5x,96i4)') (j,j=k*96+1,k*96+96)
         do l=0,n/96-1
            do i=l*96+1,l*96+96
               write(iunit,'(i5,96f4.1)') i,(m(i,j),j=k*96+1,k*96+96)
            end do
         end do
            do j=1,96
            do i=1,n
               csum(j)=csum(j)+m(i,k*96+j)
            end do
            end do
            write(iunit,'('' csum'',96f4.1)') (csum(j),j=1,96)
         end do
      end if
      return
      end

      subroutine wrtfld(iunit,fld,na,el,ldnn,elem,name)
cccc--------------------------------------------------------------------
c     write a vector field on the output file
c-----------------------------------------------------------------------
      implicit none
      integer iunit,na,ldnn,ibl,i,j,el(ldnn,2)
      real*8 fld(1),cs1,cs2,cs3
      character elem(2,2)*2,name*(*)
         cs1=0
         cs2=0
         cs3=0
         write(iunit,*)
         write(iunit,*) name
         write(iunit,*) '               x        y        z'
         do ibl=1,na
            i=(ibl-1)*6
            write(iunit,'(i4,2a3,1x,3f9.6)')
     *           ibl,'A',elem(1,el(ibl,1)),(fld(j),j=i+1,i+3)
            write(iunit,'(4x,2a3,1x,3f9.6)')
     *               'B',elem(2,el(ibl,2)),(fld(j),j=i+4,i+6)
            cs1=cs1+fld(i+1)+fld(i+4)
            cs2=cs2+fld(i+2)+fld(i+5)
            cs3=cs3+fld(i+3)+fld(i+6)
         end do
         write(iunit,'('' column sum'',3f9.4)') cs1,cs2,cs3
      return
      end

      subroutine nxtwrd(string,offset,word,wlen)
cccc--------------------------------------------------------------------
c     gives next word in 'string'(input) delimited by ' ' or ','
c     to 'word' (output) and return.
c     offset (input) is the length of the first word to be deleted
c     on output, offset is the end of the next word.
c     wlen (output) is the length of word
c     if there is no more word, word=' ' and wlen=0 will be returned
c-----------------------------------------------------------------------
      implicit none
      integer offset,wlen,slen,i,j
      character string*(*),word*(*),del1,del2
      del1=' '
      del2=','
      slen=len(string)
      word=' '
c     find the starting position of the next word
c     -------------------------------------------
      i=offset+1
      do while (string(i:i).eq.del1 .or. string(i:i).eq.del2)
c        if there is no more word, return word=' ' and wlen=0
c        ----------------------------------------------------
         if(i.ge.slen) then
            wlen=0
            return
         end if
c        otherwise just increase the offset
c        ----------------------------------
         i=i+1
      end do
      offset=i-1
c     now offset=(beg of next word)-1
c     find the length of the next word
c     --------------------------------
      do i=offset+1,slen
         if(string(i:i).eq.del1 .or. string(i:i).eq.del2) go to 333
      end do
333   i=i-1
c     now i=(end of next word)
c     ------------------------
      wlen=i-offset
c     assign next word
c     ----------------
      do j=1,wlen
         word(j:j)=string(offset+j:offset+j)
      end do
      offset=i
      return
      end

      subroutine setup()
cccc--------------------------------------------------------------------
c     read in, setup output file, make everything dimensionless
c-----------------------------------------------------------------------
      implicit none
      include 'allsig.inc'
      integer i,j,k,n1,n2,jseed,off,wlen,def,new
      real*8 hbar,kb,amu,avga,avgb,ran1
      character word*40,infil2*40
      ru=11
      wu=21
c     ------------------------------------------------------------------
c     read in
c     ------------------------------------------------------------------
      read(*,'(a)') job
      write(6,*) job
      if(job(1:5).eq.'check') then
         off=5
         call nxtwrd(job,off,word,wlen)
         do while(wlen.ne.0)
            write(6,'(2x,a20)') word
            if(word(1:2).eq.'nn') then
               chknn=1
            else if(word(1:1).eq.'m') then
               chkm=1
            else if(word(1:3).eq.'str') then
               chkstr=1
            else if(word(1:3).eq.'eig') then
               chkeig=1
            else
               infile=word
               call nxtwrd(job,off,word,wlen)
               outfile=word
            end if
            call nxtwrd(job,off,word,wlen)
         end do
         read(*,*) iseed,ftol,gtol
         read(*,*) emethod
      else if(job(1:6).eq.'static') then
         off=6
         call nxtwrd(job,off,word,wlen)
         infile=word
         call nxtwrd(job,off,word,wlen)
         outfile=word
         read(*,*) iseed,ftol,gtol
      else if(job(1:5).eq.'eigen') then
         off=5
         call nxtwrd(job,off,word,wlen)
         infile=word
         call nxtwrd(job,off,word,wlen)
         outfile=word
         read(*,*) emethod
      else if(job(1:7).eq.'thermal') then
         off=7
         call nxtwrd(job,off,word,wlen)
         infile=word
         call nxtwrd(job,off,word,wlen)
         infil2=word
         call nxtwrd(job,off,word,wlen)
         outfile=word
         read(*,*) temp
      else
         write(6,*) '*****'
         write(6,*) '***** invalid input file ''allsig.in'''
         write(6,*) '*****'
         stop
      end if
c -------------------- old method --------------------
c      open(ru,file='allsig.in',status='old')
c      read(ru,'(a)') job
c      write(6,*) job
c      if(job(1:5).eq.'check') then
c         off=5
c         call nxtwrd(job,off,word,wlen)
c         do while(wlen.ne.0)
c            write(6,'(2x,a20)') word
c            if(word(1:2).eq.'nn') then
c               chknn=1
c            else if(word(1:1).eq.'m') then
c               chkm=1
c            else if(word(1:3).eq.'str') then
c               chkstr=1
c            else if(word(1:3).eq.'eig') then
c               chkeig=1
c            else
c               infile=word
c               call nxtwrd(job,off,word,wlen)
c               outfile=word
c            end if
c            call nxtwrd(job,off,word,wlen)
c         end do
c         read(ru,*) iseed,ftol,gtol
c         read(ru,*) emethod
c      else if(job(1:6).eq.'static') then
c         off=6
c         call nxtwrd(job,off,word,wlen)
c         infile=word
c         call nxtwrd(job,off,word,wlen)
c         outfile=word
c         read(ru,*) iseed,ftol,gtol
c      else if(job(1:5).eq.'eigen') then
c         off=5
c         call nxtwrd(job,off,word,wlen)
c         infile=word
c         call nxtwrd(job,off,word,wlen)
c         outfile=word
c         read(ru,*) emethod
c      else if(job(1:7).eq.'thermal') then
c         off=7
c         call nxtwrd(job,off,word,wlen)
c         infile=word
c         call nxtwrd(job,off,word,wlen)
c         infil2=word
c         call nxtwrd(job,off,word,wlen)
c         outfile=word
c         read(ru,*) temp
c      else
c         write(6,*) '*****'
c         write(6,*) '***** invalid input file ''allsig.in'''
c         write(6,*) '*****'
c         stop
c      end if
c      close(ru)
c     ------------------------------------------------------------------
c     read in data
c     ------------------------------------------------------------------
      open(ru,file=infile,status='old',err=990)
      read(ru,'(4a2)') ((elem(i,j),j=1,2),i=1,2)
      read(ru,*) x,y
      read(ru,*) lx,ly,lz
      read(ru,*) ((radius(i,j),j=1,2),i=1,2)
      read(ru,*) ((mass(i,j),j=1,2),i=1,2)
      read(ru,*) ((alpha(1,i,j),j=1,2),i=1,2)
      read(ru,*) ((alpha(2,i,j),j=1,2),i=1,2)
      read(ru,*) ((beta(1,1,i,j),j=1,2),i=1,2)
      read(ru,*) ((beta(1,2,i,j),j=1,2),i=1,2)
      read(ru,*) ((beta(2,1,i,j),j=1,2),i=1,2)
      read(ru,*) ((beta(2,2,i,j),j=1,2),i=1,2)
ccc following line is modified (jan99)
c      read(ru,*) lambda
      read(ru,*) lambda,control
      go to 995
990   write(6,*) '*** open error file ',infile
      stop
995   continue
c     ------------------------------------------------------------------
c     dimension of dynamical matrix
c     ------------------------------------------------------------------
      na=lx*ly*lz*4
      if(na.gt.maxn) stop '*** too large l ***'
      nd=na*6
c     ------------------------------------------------------------------
c     if 'eigen' or 'thermal',
c         read in random elements and static displacement
c     otherwise,
c         assign random elements of an alloy A_1-x B_x C_1-y D_y
c         and get actual fraction x and y
c     ------------------------------------------------------------------
      if(job(1:5).eq.'eigen' .or. job(1:7).eq.'thermal') then
c---     read in random elements and static displacements
         read(ru,*) iseed
         read(ru,*)
ccc following line is added
         read(ru,*) le,a0
         do i=1,na
            j=(i-1)*6
            read(ru,'(i2,3f18.14)') el(i,1),u(j+1),u(j+2),u(j+3)
            read(ru,'(i2,3f18.14)') el(i,2),u(j+4),u(j+5),u(j+6)
         end do
         if(job(1:7).eq.'thermal') then
c---        read in eigenvalues and eigenvectors
            infile=infil2
            open(14,file=infile,status='old',form='unformatted',err=990)
            do i=1,nd
               read(14) dd(i),(d(i,j),j=1,nd)
            end do
            close(14)
         end if
      else
c---     if the job is not 'eigen' or 'thermal'
c---     initialize random number generator
         jseed=iseed
         i=ran1(-jseed)
c---     assign random elements
c        A sublattice
         if(x.le.0.5d0) then
            def=1
            new=2
            n1=na*x+.000001
         else
            def=2
            new=1
            n1=na*(1-x)+.000001
         end if
         do i=1,na
            el(i,1)=def
         end do
         n2=0
         do while (n2.lt.n1)
            i=int(ran1(jseed)*na)+1
            if(el(i,1).eq.def) then
               el(i,1)=new
               n2=n2+1
            end if
         end do
         if(x.le.0.5d0) then
            x=dfloat(n2)/na
         else
            x=1-dfloat(n2)/na
         end if
c        B sublattice
         if(y.le.0.5d0) then
            def=1
            new=2
            n1=na*y+.000001
         else
            def=2
            new=1
            n1=na*(1-y)+.000001
         end if
         do i=1,na
            el(i,2)=def
         end do
         n2=0
         do while (n2.lt.n1)
            i=int(ran1(jseed)*na)+1
            if(el(i,2).eq.def) then
               el(i,2)=new
               n2=n2+1
            end if
         end do
         if(y.le.0.5d0) then
            y=dfloat(n2)/na
         else
            y=1-dfloat(n2)/na
         end if
ccc following lines are moved from further down (jan99)
c        le is the nn distance (in A) as given in the virtual crystal
         le=radius(1,1)*(1-x)+radius(1,2)*x
     *     +radius(2,1)*(1-y)+radius(2,2)*y
ccc test for EMT of Al P_1-x Sb_x
c         le=y*2.6559+(1-y)*2.3658+y*(1-y)*(-0.0537)
c        a0 is the half the lattice constant (in A) of the cubic unit cell
         a0=2*le/sqrt(3.d0)
ccc end of movement (jan99)
      end if
      close(ru)
c     ------------------------------------------------------------------
c     setup output file
c     ------------------------------------------------------------------
      if(job(1:5).eq.'eigen') then
         open(wu,file=outfile,status='new',form='unformatted')
      else
c         open(wu,file=outfile,status='new',reclen=400)
         open(wu,file=outfile,status='new')
         write(wu,'(4a2,     t45,''elements'')')
     *         ((elem(i,j),j=1,2),i=1,2)
         write(wu,'(1x,2f9.5,t45,''fraction'')')
     *         x,y
         write(wu,'(1x,3i9,  t45,''L'')')
     *         lx,ly,lz
         write(wu,'(1x,4f9.5,t45,''radii'')')
     *         ((radius(i,j),j=1,2),i=1,2)
         write(wu,'(1x,4f9.4,t45,''masses'')')
     *         ((mass(i,j),j=1,2),i=1,2)
         write(wu,'(1x,4f9.3,t45,''alpha'')')
     *         ((alpha(1,i,j), j=1,2),i=1,2)
         write(wu,'(1x,4f9.3,t45,''alpha'')')
     *         ((alpha(2,i,j), j=1,2),i=1,2)
         write(wu,'(1x,4f9.3,t45,''beta '')')
     *         ((beta(1,1,i,j),j=1,2),i=1,2)
         write(wu,'(1x,4f9.3,t45,''beta '')')
     *         ((beta(1,2,i,j),j=1,2),i=1,2)
         write(wu,'(1x,4f9.3,t45,''beta '')')
     *         ((beta(2,1,i,j),j=1,2),i=1,2)
         write(wu,'(1x,4f9.3,t45,''beta '')')
     *         ((beta(2,2,i,j),j=1,2),i=1,2)
ccc following line is modified (jan99)
c        write(wu,'(1x,f9.3, t45,''lambda'')')   lambda
         write(wu,'(1x,f9.3,i9,t45,''lambda  control'')') lambda,control
         write(wu,'(1x,i9,   t45,''iseed'')')    iseed
         if(job(1:7).eq.'thermal') then
            write(wu,'(1x,f9.2, t45,''temp'')')     temp
         end if
      end if

      write(6,'(1x,4a3,  t45,''elements'')') ((elem(i,j),j=1,2),i=1,2)
      write(6,'(1x,2f9.4,t45,''fraction'')') x,y
      write(6,'(1x,3i9,  t45,''L'')')        lx,ly,lz
      write(6,'(1x,4f9.3,t45,''alpha'')') ((alpha(1,i,j), j=1,2),i=1,2)
      write(6,'(1x,4f9.3,t45,''beta '')') ((beta(1,1,i,j),j=1,2),i=1,2)
ccc following line is modified (jan99)
c     write(6,'(1x,f9.3, t45,''lambda'')')   lambda
      write(6,'(1x,f9.3,i9,t45,''lambda  control'')') lambda,control
      write(6,'(1x,a20,t45,''outfile'')')    outfile
c     ------------------------------------------------------------------
c     make everything dimensionless
c     ------------------------------------------------------------------
c     masses
c     ------
      avga=mass(1,1)*(1-x)+mass(1,2)*x
      avgb=mass(2,1)*(1-y)+mass(2,2)*y
      mu=(avga*avgb)/(avga+avgb)
      do i=1,2
      do j=1,2
        mass(i,j)=mass(i,j)/mu
      end do
      end do
c     force constants
c     ---------------
      avga=0
      avgb=0
      do i=1,2
      do j=1,2
         avga=avga+alpha(1,i,j)+alpha(2,i,j)
      do k=1,2
         avgb=avgb+beta(1,i,j,k)+beta(2,i,j,k)
      end do
      end do
      end do
      avga=avga/8
      avgb=avgb/16
      cc=(4.d0/3)*(avga+avgb+7*lambda*avgb/9)
      do i=1,2
      do j=1,2
         alpha(1,i,j)=alpha(1,i,j)/cc
         alpha(2,i,j)=alpha(2,i,j)/cc
      do k=1,2
         beta(1,i,j,k)=beta(1,i,j,k)/cc
         beta(2,i,j,k)=beta(2,i,j,k)/cc
      end do
      end do
      end do
c     maximums
c     --------
      hbar=1.05457266d-34
      kb=1.380658d-23
      amu=1.6605402d-27
      wmax=sqrt(cc/(mu*amu))
      tmax=hbar*wmax/kb
      xmax=sqrt(2*hbar*wmax/cc)*1.d10
      temp=temp/tmax
      if(job(1:5).ne.'eigen')
     *   write(wu,'(1x,3e11.3, t45,''wmax xmax tmax'')') wmax,xmax,tmax
      return
      end
