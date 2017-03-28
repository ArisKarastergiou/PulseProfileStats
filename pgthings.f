c -----------------------------------------------------------------------
	subroutine boxes(v1,v2,v3,v4,xmin,xmax,ymin,ymax,
     +                   xopt,yopt,xlab,ylab,tlab)
	implicit none
	real v1,v2,v3,v4,xmin,xmax,ymin,ymax
	character*(*) xopt,yopt,xlab,ylab,tlab

	call pgsci(1)
	call pgsvp(v1,v2,v3,v4)
	call pgswin(xmin,xmax,ymin,ymax)
	call pgbox(xopt,0.0,0,yopt,0.0,0)
	call pglabel(xlab,ylab,tlab)

	return
	end
c -----------------------------------------------------------------------
	subroutine ptline(c,n,x,y)
	implicit none
	integer c,n
	real x(*),y(*)

	call pgsci(c)
	call pgpt(n,x,y,1)
	call pgline(n,x,y)
	call pgsci(1)

	return
	end
c -----------------------------------------------------------------------
