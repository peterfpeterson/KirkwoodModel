      integer maxl,maxn,maxd
      parameter (maxl=4,maxn=maxl*maxl*maxl*4,maxd=maxn*2*3+1)
c     define variables
c-----------------------------------------------------------------------
      integer nn,ldnn,ldm,ldd,lx,ly,lz,na,nd,el,iseed,control
      integer ru,wu,chknn,chkm,chkstr,chkeig,emethod
      character job*80,elem*2,infile*40,outfile*40
      real*8        m,d,unitv,phi,dd,sd
      real*8        x,y,le,a0,u,radius,mass,mu, alpha,beta,lambda,cc
      real*8        wmax,tmax,xmax,temp,sigma,ftol,gtol
      common /ctrl/ job,chknn,chkm,chkstr,chkeig,emethod
      common /file/ ru,wu,infile,outfile
c      common /parm/ ldnn,ldm,ldd,lx,ly,lz,na,nd,ftol,gtol
      common /parm/ ldnn,ldm,ldd,lx,ly,lz,na,nd
      common /matr/ m(maxd,maxd),d(maxd,maxd)
      common /vect/ dd(maxd),sd(maxd)
      common /stru/ unitv(4,2,3),nn(maxn,2,4)
      common /conf/ x,y,le,a0,phi(maxd),el(maxn,2),iseed,control
      common /atom/ u(maxd),radius(2,2),mass(2,2),elem(2,2),mu
      common /forc/ alpha(2,2,2),beta(2,2,2,2),lambda,cc
      common /cnst/ wmax,tmax,xmax,temp
      common /sigm/ sigma
