c==------------------------------------------------------------------------
c Version: 14-DEC-1994 (v3.0)
c File : cgmin.f
c
c   Fortran subrutines to perform function minimization using the
c   conjugate gradient method in double precision.
c
c Authors:
c
c   Original version: FRPRMN.FOR
c         William H. Press, Saul A. Teukolsky, William T. Vetterling,
c         Brian P. Flannery,
c         "Numerical Recipes in FORTRAN: The Art of Scientific Computing",
c         2nd Edition, Cambridge (1992).
c
c   Version 1.0 :
c       JeanSoo Chung
c       Department of Physics
c       CHUNGBUK NATIONAL UNIVERSITY
c       Kae-Shin Dong, Cheong-Ju, Chung-Cheong-Buk-Do 360-763, KOREA
c       FAX : 011-82-431-63-0612
c       chung@cbucc.chungbuk.ac.kr
c
c   Verion 2.0 and higer:
c
c       Seong Gon Kim
c       Department of Physics and Astronomy
c       MICHIGAN STATE UNIVERSITY
c       East Lansing, Michigan 48824.
c       kimsg@msupa.pa.msu.edu
c
c Usage:
c
c      call cgmin (n,p,g,ftol,gtol,mxitr,
c     &  iter,fret,fchg,res,enersub,gradsub,monsub,
c     &  wk,nwk,irprt,lumsg)
c
c Arguments:
c
c     n      : dimension of the configuration space (input).
c     p(1:n) : initial guess of the configuration vector (input)
c              and calculated (or optimized) configuration vector
c              at which the enersub has minimum value (output).
c               -- real*8 array of size n.
c     g(1:n) : force, or "negative" gradient vector at minimum (output)
c               -- real*8 array of size n.
c     ftol   : tolerence for function value variation (input).
c     gtol   : tolerence for residual gradient (input).
c     mxitr  : maximum no. of iterations to perform (input).
c     iter   : no. of iteration performed (output).
c     fret   : minimized function value (output).
c     fchg   : last change in function value (output).
c     res    : residual force magnitude (output).
c     enersub: subroutine to evaluate the energy value (input).
c     gradsub: subroutine to evaluate the gradient vector (input).
c     monsub : subroutine to perform monitoring (input).
c     wk(nwk): work space.
c               -- real*8 array of size nwk at least.
c     nwk    : size of the work space (input).
c     irprt  : report index (input).
c              zero or negative value suppresses the message.
c     lumsg : logical unit for message printing.
c
c     "enersub","gradsub" and "monsub" should be declared as external
c       in the calling program and provided by the user with the
c       following structures :
c
c      subroutine enersub (n,p,f)
c      subroutine gradsub (n,p,g,gg0)
c      subroutine monsub (iter,p,g,fp,fchg,res)
c
c Algorithm:
c     Press et.al., Numerical Recipes, (Cambridge Univ. Press,1986).
c       See routine FRPRMN (Sec. 10.6).
c
c Routines referenced:
c
c      subroutine linmin (n,p,xi,pcom,xicom,wk1,wk2,
c     &  fret,enersub,gradsub,irprt,lumsg)
c      real*8 function f1dim (x,n,pcom,xicom,xt,enersub)
c      real*8 function df1dim (x,n,pcom,xicom,xt,df,gradsub)
c      subroutine mnbrak (ax,bx,cx,fa,fb,fc,n,pcom,xicom,
c     &  wk1,func,enersub)
c      real*8 function dbrent (ax,bx,cx,n,pcom,xicom,wk1,wk2,
c     &  func,dfunc,enersub,gradsub,tol,xmin,irprt,lumsg)
c***********************************************************************
c==
      subroutine cgmin (n,p,g,ftol,gtol,mxitr,
     &  iter,fret,fchg,res,enersub,gradsub,monsub,
     &  wk,nwk,irprt,lumsg)

c.... parameters.

      integer*4 n,mxitr,iter,nwk,irprt,lumsg
      real*8 p(*),g(*),ftol,gtol,fret,fchg,res,wk(*)
      external enersub,gradsub,monsub

c.... numerical constants.

      real*8 zero,one,two,three,eps,tiny
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter (eps=1.0d-10,tiny=1.0d-30)

c.... local variables.

      integer*4 j
      real*8 fp,dgg,gam,gg,gg0

c.... check work space.

      if (nwk .lt. n*6) stop '%cgmin: Not enough work space.'

c.... initialize.

      iter=0

      call enersub (n,p,fp)
      call gradsub (n,p,wk(n+1),gg0)

      gg=zero
      do j=1,n
        g(j)=-wk(n+j)
        wk(j)=g(j)
        wk(n+j)=wk(j)
        gg=gg+g(j)*g(j)
      enddo
      res=sqrt(gg0)

c.... initial guess is the solution.

      if (gg .le. tiny) then
        fret=fp
        fchg=0.0
        goto 100
      endif

c*******************************************************************
c Beginning of the main loop.
c*******************************************************************

      do iter=1,mxitr

c.... minimize along current directional vector.

        call linmin (n,p,wk(n+1),wk(n*2+1),wk(n*3+1),
     &    wk(n*4+1),wk(n*5+1),fret,enersub,gradsub,irprt,lumsg)

c.... function value change.

        fchg=fp-fret

c.... update function value and gradient vector.

        call enersub (n,p,fp)
        call gradsub (n,p,wk(n+1),gg0)

c.... update directional vector.

        dgg=zero
        do j=1,n
          dgg=dgg+(wk(n+j)+g(j))*wk(n+j)
        enddo

c.... new residual gradient sum.

        gam=dgg/gg
        gg=zero
        do j=1,n
          g(j)=-wk(n+j)
          gg=gg+g(j)*g(j)
          wk(j)=g(j)+gam*wk(j)
          wk(n+j)=wk(j)
        enddo
        res=sqrt(gg0)

c.... monitor progress of minimization.

        call monsub (iter,p,g,fp,fchg,res)

c.... check if converged.

        if ((abs(fchg) .le. ftol) .and. (res .le. gtol)) goto 120
c        if (fchg .eq. zero) goto 110
        if (gg .le. tiny) goto 100

      enddo

c*******************************************************************
c End of the main loop.
c*******************************************************************

c.... maximum iteration exceeds.

      iter=mxitr
      if (irprt .gt. 0) then
        write (lumsg,600)
        write (lumsg,610)
     &    '%cgmin: Failed to converge until ending iteration.'
      endif
      goto 200

c.... zero gradient vector.

  100 continue
      if (irprt .gt. 1) then
        write (lumsg,600)
        write (lumsg,610)
     &    '%cgmin: Converged with zero gradient vector.'
      endif
      goto 200

c.... zero energy change.

c  110 continue
c      if (irprt .gt. 1) then
c        write (lumsg,600)
c        write (lumsg,610)
c     &    '%cgmin: Converged with zero energy change.'
c      endif
c      goto 200

c.... converged successfully.

  120 continue
      if (irprt .gt. 1) then
        write (lumsg,600)
        write (lumsg,610) '%cgmin: Converged successfully.'
      endif

c.... write convergence report.

  200 if (irprt .gt. 1) then
        write (lumsg,610) '% ====  cgmin : CONEVRGENCE REPORT  ==== '
        write (lumsg,620) '% Iteration number   = ',iter
        write (lumsg,630) '% Last value change  = ',fchg
        write (lumsg,630) '% Residual force     = ',res
        write (lumsg,630) '% Function tolerance = ',ftol
        write (lumsg,630) '% Gradient tolerance = ',gtol
      endif

c.... done.

      return

c.... formats.

  600 format (60('='))
  610 format (a)
  620 format (a,i12)
  630 format (a,1pe22.15e2)

      end
c==
      subroutine linmin (n,p,xi,pcom,xicom,wk1,wk2,
     &  fret,enersub,gradsub,irprt,lumsg)
c***********************************************************************
c Line minimization subroutine.
c
c Input:
c   n : dimension of the vectors.
c   p : configuration vector.
c   xi : initial gradient vector.
c   pcom : working copy of p.
c   xicom : working copy of xi.
c   wk1,wk2 : work space.
c   enersub: subroutine to evaluate the energy value.
c   gradsub: subroutine to evaluate the gradient vector.
c   irprt  : report index (input).
c            zero or negative value suppresses the message.
c   lumsg : logical unit for message printing.
c
c Output:
c   fret : returned minimized function value.
c***********************************************************************

c.... parameters.

      integer*4 n,irprt,lumsg
      real*8 p(*),xi(*),pcom(*),xicom(*),fret,wk1(*),wk2(*)
      external enersub,gradsub

c.... numerical constants.

      real*8 tol
      parameter (tol=1.0d-6)

c.... local variables.

      integer j
      real*8 ax,bx,cx,fa,fb,fc,xmin

c.... external routines.

      external f1dim,df1dim
      real*8 dbrent

c.... make working copy.

      do j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
      enddo

c.... bracket the minimum.

      ax=0.0
      bx=0.5
      cx=1.0

      call mnbrak2 (ax,bx,cx,fa,fb,fc,n,pcom,xicom,
     &  wk1,f1dim,enersub)

c.... find minimum along this direction.

      fret=dbrent (ax,bx,cx,n,pcom,xicom,wk1,wk2,
     &  f1dim,df1dim,enersub,gradsub,tol,xmin,irprt,lumsg)

c.... new direction vector and optimized configuration vector.

      do j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
      enddo

      end
c==
      real*8 function f1dim (x,n,pcom,xicom,xt,enersub)

c.... parameters.

      integer*4 n
      real*8 x,pcom(*),xicom(*),xt(*)
      external enersub

c.... local variables

      integer*4 j
      real*8 fval

      do j=1,n
        xt(j)=pcom(j)+x*xicom(j)
      enddo

      call enersub (n,xt,fval)
      f1dim=fval

      end
c==
      real*8 function df1dim (x,n,pcom,xicom,xt,df,gradsub)

c.... parameters.

      integer*4 n
      real*8 x,pcom(*),xicom(*),xt(*),df(*),gg0
      external gradsub

c.... local variables.

      integer*4 j

      do j=1,n
        xt(j)=pcom(j)+x*xicom(j)
      enddo

      call gradsub (n,xt,df,gg0)

      df1dim=0.0
      do j=1,n
        df1dim=df1dim+df(j)*xicom(j)
      enddo

      end
c==
      subroutine mnbrak2 (ax,bx,cx,fa,fb,fc,n,pcom,xicom,
     &  wk1,func,enersub)

c.... parameters.

      real*8 ax,bx,cx,fa,fb,fc
      integer*4 n
      real*8 pcom(*),xicom(*),wk1(*)
      real*8 func
      external enersub

c.... numerical constants.

      real*8 gold,glimit,tiny
      parameter (gold=1.618034d00,glimit=100.0d00,tiny=1.0d-10)

c.... local variables.

      real*8 r,q,u,ulim,fu

      fa=func (ax,n,pcom,xicom,wk1,enersub)
      fb=func (bx,n,pcom,xicom,wk1,enersub)
      fc=func (cx,n,pcom,xicom,wk1,enersub)

c.... if f(b) > f(a), min. should be in (a,b).

    5 if (fb .gt. fa) then
        cx=bx
        fc=fb
        bx=0.5*(cx-ax)
        fb=func (bx,n,pcom,xicom,wk1,enersub)
        goto 5
      endif

c.... exit, if f(c) > f(b).

   10 if (fc .gt. fb) return

c.... compute u by parabolic extrapolation from a,b,c.
c.... "tiny" is used to prevent any possible division by zero.

      r=(bx-ax)*(fb-fc)
      q=(bx-cx)*(fb-fa)
      u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),tiny),q-r))

c.... we won't go farther than glimit.

      ulim=bx+glimit*(cx-bx)

      if ((bx-u)*(u-cx) .gt. 0.0) then

c.... if b<u<c, try it.

        fu=func (u,n,pcom,xicom,wk1,enersub)

        if (fu .lt. fc) then

c.... got a minimum between b and c.

          ax=bx
          fa=fb
          bx=u
          fb=fu
          goto 10

        elseif (fu .gt. fb) then

c.... got a minimum between a and u.

          cx=u
          fc=fu
          goto 10

        endif

c.... parabolic fit was no use. Use default magnification.

        u=cx+gold*(cx-bx)
        fu=func (u,n,pcom,xicom,wk1,enersub)

      elseif ((cx-u)*(u-ulim) .gt. 0.0) then

c.... if c<u<ulim,

        fu=func (u,n,pcom,xicom,wk1,enersub)
        if (fu .lt. fc) then
          bx=cx
          cx=u
          u=cx+gold*(cx-bx)
          fb=fc
          fc=fu
          fu=func (u,n,pcom,xicom,wk1,enersub)
        endif

      elseif ((u-ulim)*(ulim-cx) .ge. 0.0) then

c.... limit parabolic u to maximum allowed value.

        u=ulim
        fu=func (u,n,pcom,xicom,wk1,enersub)

      else

c.... reject parabolic u, use default magnification.

        u=cx+gold*(cx-bx)
        fu=func (u,n,pcom,xicom,wk1,enersub)

      endif

c.... move down hill.

      ax=bx
      bx=cx
      cx=u
      fa=fb
      fb=fc
      fc=fu

c.... continue until bracketed.

      goto 10

      end
c==
      real*8 function dbrent (ax,bx,cx,n,pcom,xicom,wk1,wk2,
     &  func,dfunc,enersub,gradsub,tol,xmin,irprt,lumsg)

c.... parameters.

      integer*4 n,irprt,lumsg
      real*8 ax,bx,cx
      real*8 pcom(*),xicom(*),wk1(*),wk2(*)
      real*8 func,dfunc
      external enersub,gradsub
      real*8 tol,xmin

c.... numerical constants.

      integer*4 iend
      real*8 zero,zeps
      parameter (iend=100,zero=0.0d00,zeps=1.0d-15)
      real*8 mzero
      parameter (mzero=1.0e-30)

c.... local variables.

      real*8 a,b,v,w,x,e,fx,fv,fw,dx,dv,dw,xm,tol1,tol2
      real*8 u1,u2,d1,d2,olde,d,u,fu,du
      integer*4 iter
      logical*4 ok1,ok2

      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=zero
      fx=func (x,n,pcom,xicom,wk1,enersub)
      fv=fx
      fw=fx
      dx=dfunc (x,n,pcom,xicom,wk1,wk2,gradsub)
      dv=dx
      dw=dx

      do iter=1,iend

        xm=0.5*(a+b)
        tol1=tol*abs(x)+zeps
        tol2=2.0d00*tol1

*********** exit test. *****************************************
        if (abs(x-xm) .le. (tol2-0.5d00*(b-a))) goto 30
****************************************************************

        if (abs(e) .gt. tol1) then
          d1=2.0d00*(b-a)
          d2=d1

c-----------------------------------------------------------------
          if (dw .ne. dx) d1=(w-x)*dx/(dx-dw)
          if (dv .ne. dx) d2=(v-x)*dx/(dx-dv)
c          if (abs(dw-dx) .gt. mzero) d1=(w-x)*dx/(dx-dw)
c          if (abs(dv-dx) .gt. mzero) d2=(v-x)*dx/(dx-dv)
c-----------------------------------------------------------------

          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b) .gt. zero) .and. (dx*d1 .le. zero)
          ok2=((a-u2)*(u2-b) .gt. zero) .and. (dx*d2 .le. zero)
          olde=e
          e=d
          if (.not. (ok1 .or. ok2)) then
            goto 10
          elseif (ok1 .and. ok2) then
            if (abs(d1) .lt. abs(d2)) then
              d=d1
            else
              d=d2
            endif
          elseif (ok1) then
            d=d1
          else
            d=d2
          endif
          if (abs(d) .gt. abs(0.5d00*olde)) goto 10
          u=x+d
          if (((u-a) .lt. tol2) .or.
     &        ((b-u) .lt. tol2)) d=sign(tol1,xm-x)
          goto 20
        endif

   10   if (dx .ge. zero) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5*e
   20   if (abs(d) .ge. tol1) then
          u=x+d
          fu=func (u,n,pcom,xicom,wk1,enersub)
        else
          u=x+sign(tol1,d)
          fu=func (u,n,pcom,xicom,wk1,enersub)
          if (fu .gt. fx) goto 30
        endif
        du=dfunc (u,n,pcom,xicom,wk1,wk2,gradsub)
        if (fu .le. fx) then
          if (u .ge .x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if (u .lt. x) then
            a=u
          else
            b=u
          endif
          if ((fu .le. fw) .or. (w .eq. x)) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          elseif ((fu .le. fv) .or. (v .eq. x) .or. (v .eq. w)) then
            v=u
            fv=fu
            dv=du
          endif
        endif
      enddo

c.... maximum iteration exceeds.

      if (irprt .gt. 0) then
        write (lumsg,*) '%%%%%%%  WARNING  %%%%%%%'
        write (lumsg,*) '%cgmin-dbrent: Maximum iteration exceeded.'
        write (lumsg,*) '%%%%%%%  WARNING  %%%%%%%'
      endif

c.... normal exit.

   30 xmin=x
      dbrent=fx

      end
