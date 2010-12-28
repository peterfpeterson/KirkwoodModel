      program convol
C----------------------------------------------------------------------
C-
C-   Purpose and Methods :
C-      Calculates experimental Pair Distribution Function (convoluted
C-      with the finite window of experiment) of a zns structure
C-      using the output of ZNSPDF (theoretical pdf).
C-      cf. Chung & Thorpe, PRB 55, 1545 (1997) Eq. B9
C-      Output is written as a topdrawer file.
C-
C-   Inputs  :
C-   Outputs :
C-   Controls:
C-
C-   Created  18-JUL-1997   Jean Soo Chung
C-   Modified    jan-1999
C-
C----------------------------------------------------------------------
      implicit none
      integer maxdat
      real*8  pi,pix2
      parameter (maxdat=9000,pi=3.14159265358979d0,pix2=2*pi)
      real*8  r(0:maxdat),gt(0:maxdat),ge(0:maxdat),del,qmax,t1,t2,
     *        rij(-maxdat:2*maxdat),sij(-maxdat:2*maxdat)
      integer ndat,nerr,i,j,skip
      character line*40,file*30,data(5)*40
c----------------------------------------------------------------------
c     initial setup
c----------------------------------------------------------------------
c     read in
c     -----------------------------------------------------------------
120   write(6,*) 'enter input file name'
      read(5,'(a)') file
      open(unit=11,file=file,status='old',err=120)
140   write(6,*) 'enter output TopDrawer file name'
      read(5,'(a)') file
      open(unit=21,file=file,status='new',err=140)
160   write(6,*) 'enter max q value (in A^-1)'
      read(5,'(a)') line
      read(line,*,err=160) qmax
      ndat=0
      nerr=0
      do i=1,maxdat
         read(11,'(a)',end=300) line
c         if(i.lt.15) write(6,*) line
         read(line,*,err=220) t1,t2
         if(t1.le.0) goto 220
         ndat=ndat+1
         r(ndat)=t1
         gt(ndat)=t2
         if(ndat.le.5) data(ndat)=line
200      continue
      end do
      goto 300
c-----begin  error handling
220      write(*,*) '*** error reading : '//line
         nerr=nerr+1
c         goto 200
c-----end of error handling
300   continue
      write(6,'(i5,'' lines read from input file, starting'')') ndat
      do i=1,5
         write(6,*) '   '//data(i)
c         write(6,*) r(i),gt(i)
      end do
      close(11)
c         nout=5000
c550      write(6,*) 'enter no. of output points [default:5000]'
c         read(5,'(a20)') file
c         if(file(1:1).ne.' ') read(file,*,err=550) nout
c     preparation
c     -----------------------------------------------------------------
      r(0)=0
      gt(0)=0
      rij(0)=0
      sij(0)=0
      do i=1,ndat
         rij(i)=r(i)
         rij(-i)=-r(i)
         rij(ndat+i)=r(ndat)+r(i)
         sij(i)=sin(qmax*rij(i))
         sij(-i)=-sij(i)
         sij(ndat+i)=sin(qmax*rij(ndat+i))
      end do
      del=r(2)-r(1)
c----------------------------------------------------------------------
c     convolution and write
c----------------------------------------------------------------------
      skip=ndat/2000+1
      write(6,'(i5,a23)') ndat/skip,' points will be written'
      do i=1,ndat,skip
         ge(i)=0
         do j=1,ndat
            if(i-j.eq.0) then
               t1=qmax
            else
               t1=sij(i-j)/rij(i-j)
            end if
            if(i+j.eq.0) then
               t2=qmax
            else
               t2=sij(i+j)/rij(i+j)
            end if
c        ...next line has two convolution terms
c            ge(i)=ge(i)+gt(j)*(t1-t2)
c        ...next line has one convolution term
            ge(i)=ge(i)+gt(j)*(t1)
         end do
         ge(i)=ge(i)*del/pi
         write(21,'(2f15.5)') r(i),ge(i)
      end do
         write(21,'('' join'')')
      stop
      end
