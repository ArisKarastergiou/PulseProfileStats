c generic routine for creating a histogram
c running between xmin and xmax with nh bins
c the data is in x which is of dimension nx
c returns the histogram in hist and the x-axis in xhist
c    and ymin,ymax and tgood
	subroutine dohist(x,nx,xmin,xmax,ymin,ymax,tgood,hist,xhist,nh)

	implicit none
	integer nx,nh
	real x(nx),hist(nh),xhist(nh)
	integer i,hbin,tgood
	real xmin,xmax,ymin,ymax,step

c compute the step size and form the x array
	step=(xmax-xmin)/nh
	do i=1,nh
	  hist(i)=0.0
	  xhist(i)=xmin+(i-0.5)*step
	enddo

c decide which bin the data falls in
c and update the number of things in that histogram step
c check to make sure not out of bounds
	tgood=0
	do i=1,nx
	  hbin=int((x(i)-xmin)/step)+1
	  if(hbin.lt.1 .or. hbin.gt.nh)then
	    tgood=tgood+1
	  else
	    hist(hbin)=hist(hbin)+1
	  endif
	enddo
	tgood=nx-tgood

c compute the maximum and minimum of the histogram
	ymin=1.e9
	ymax=-1.e9
	do i=1,nh
	  if(hist(i).gt.ymax)ymax=hist(i)
	  if(hist(i).lt.ymin)ymin=hist(i)
	enddo
	ymax=1.05*ymax

c write out to junk file
	do i=1,nh
	  write(77,*)i,xhist(i),hist(i)
	enddo

	return
	end
