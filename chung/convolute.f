	program conv
c*****7*******************************************************************
c
c	Program to convolute PDF with SINC function 
c
c	Version:  1.0
c	Author:   Thomas Proffen, MSU
c
c*****7*******************************************************************
c
	parameter	(nmax=2500)
c
	character*80	line,ifile,ofile
	real		calc(nmax),r(nmax),sinc(2*nmax)
c
c------ first read input file
c
	write (*,*) 'Give PDF input file name : '
	read (*,9000) ifile
	write (*,*) 'Give PDF output file name : '
	read (*,9000) ofile
	write (*,*) 'Give Qmax : '
	read (*,*) qmax

c
c------ now read observed PDF
c
	open (11,file=ifile,status='old')
c
c------ Any other header lines starting with # ?
c
44      continue
          read(11,9000) line
        if (line(1:1).eq.'#') goto 44
        backspace(11)
c
c------ start reading data 
c
	ip = 1
10	continue
	  read (11,*,end=20) r(ip),calc(ip)
	  ip = ip+1
	  if (ip.gt.nmax) then
	    close (11)
	    write (*,'(a)') '*** ERROR: Too many points in PDF ***'
	    stop
	  endif
	goto 10
20	continue
	close (11)
c
c------ Now we do the convolution
c
	call setup(sinc,r,qmax,nmax)
	call convolute(r,calc,sinc,qmax,nmax,ip)
c
	open (11,file=ofile,status='unknown')
	do i=1,ip-1
	  write (11,*) r(i),calc(i)
	enddo
	close (11)
c
9000	format (a)
c
	end
c*****7*******************************************************************
	subroutine setup (sinc,r,qmax,nmax)
c+
c	Setting up SINC function etc. ...
c-
	real	sinc(2*nmax),r(nmax)
c
	rcut   = 0.0
	deltar = r(2)-r(1)
c
	if (qmax.gt.0.0) then
	  sincut = 0.025
	  rcut   = 1.0 / (qmax*sincut)
	  z      = deltar*qmax
	  nn     = int(rcut/deltar) + 100
	endif
c
	if (qmax.gt.0.0) then
	  do i=1,nn
	    sinc(i) = sin(z*float(i))/(deltar*float(i))
	  enddo
	  do i=nn+1,2*nmax
	    sinc(i) = 0.0
	  enddo
	endif
c
	end
c*****7*******************************************************************
	subroutine convolute(r,calc,sinc,qmax,nmax,ip)
c+
c	Do the convolution here ..
c-
	real		ppp(nmax),r(nmax),calc(nmax)
	real		sinc(2*nmax)
	logical		ldiff
c
	if (qmax.gt.0.0) then
	  do i=1,ip-1
	    ppp(i) = calc(i)*(qmax-sinc(2*i))
	    do k=1,i-1
	      ppp(i) = ppp(i) + calc(k)*(sinc(i-k) - sinc(i+k))
	    enddo
	    do k=i+1,ip
	      ppp(i) = ppp(i) + calc(k)*(sinc(k-i) - sinc(i+k))
	    enddo
	  enddo
	  do i=1,ip
	    calc(i) = ppp(i)*(r(2)-r(1))/3.141592654
	  enddo
	endif
c
	end
