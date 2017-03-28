	program pw 
	implicit none

	real d2r,alpha_decay
	parameter (d2r=3.141592/180.,alpha_decay=-1.e7)
	integer i,j,c,nargs,npsr,seed,tgood,nh,ninterpulse
	real patchy,this_patchy,obs_width(1000),period(1000)
	real logperiod(1000),obstmp_width(1000),simtmp_width(1000)
	real ran1,rho,alpha,beta,sim_width(1000),width,age(1000)
	real pmin,pmax,ymin,ymax,yymax,s,smin,smax
	real hmin,hmax,height,konst,law_width,diff
	real obshist(1000),simhist(1000),xhist(1000)
	real v1,v2,v3,v4
	character string*50,fname*50,device*50
	logical nosim,noobs,nodecay

c get the arguments
	do i=1,1000
	  age(i)=0.0
	enddo
	ninterpulse=0
	hmin=490.0
	hmax=490.0
	patchy=0.75
	smin=1.0
	smax=1.0
	device='/xs'
	nosim=.false.
	noobs=.false.
	nodecay=.false.
	nargs = IARGC()
	if(nargs.lt.2)then
	  write(*,'('' Specify at least file name and seed '')')
	  stop
	endif
	do i=1,nargs
	  call getarg(i,string)
	  if(index(string,'-fn').gt.0)then
	    read(string(4:),'(a)')fname
	  endif
	  if(index(string,'-dev').gt.0)then
	    read(string(5:),'(a)')device
	  endif
	  if(index(string,'-hmin').gt.0)then
	    read(string(6:),*)hmin
	  endif
	  if(index(string,'-hmax').gt.0)then
	    read(string(6:),*)hmax
	  endif
	  if(index(string,'-p').gt.0)then
	    read(string(3:),*)patchy
	  endif
	  if(index(string,'-smin').gt.0)then
	    read(string(6:),*)smin
	  endif
	  if(index(string,'-smax').gt.0)then
	    read(string(6:),*)smax
	  endif
	  if(index(string,'-ran').gt.0)then
	    read(string(5:),*)seed
	  endif
	  if(index(string,'-nosim').gt.0)then
	    nosim=.true.
	  endif
	  if(index(string,'-noobs').gt.0)then
	    noobs=.true.
	  endif
	  if(index(string,'-nodecay').gt.0)then
	    nodecay=.true.
	  endif
	enddo
	write(*,'('' filename is '',a)')fname
	write(*,'('' hmin is '',f12.2)')hmin
	write(*,'('' hmax is '',f12.2)')hmax
	write(*,'('' patchy is '',f12.2)')patchy
	write(*,'('' smin factor is '',f12.2)')smin
	write(*,'('' smax factor is '',f12.2)')smax
	konst = patchy * sqrt((9.0*3.1415*(hmin+hmax)/2.0)/(2.0*3.0e5))
	konst = konst/d2r
	write(*,'('' konst is '',f12.2)')konst
	

c read the observed pulsar values
	open(unit=10,file=fname,form='formatted',status='old')
	i=1
 10	read(10,*,end=20)width,period(i),age(i)
	obs_width(i) = log10(width)
	logperiod(i)=log10(period(i))
	i=i+1
	goto 10

 20	continue
	npsr=i-1
	do i=1,npsr
 25	  height = hmin + (hmax-hmin)*ran1(seed)
	  s = smin + (smax-smin)*ran1(seed)
	  rho = s * sqrt(height/period(i)) * sqrt(9.0*3.1415/2.0/3.0e5)
	  alpha = acos(ran1(seed))
c	  alpha = (ran1(seed)-0.5) * 3.14159
c	  alpha = 60.0 * d2r
	  beta = ran1(seed) * rho
c	  beta = 0.0

c add in alpha decay
	  if(.not. nodecay)
     +	    alpha = alpha * exp(age(i)/alpha_decay)

	  width=(cos(rho)-cos(alpha)*cos(alpha+beta)) / 
     +          (sin(alpha)*sin(alpha+beta))
	  if(width.gt.1.0 .or. width.lt.-1.0)then
	    sim_width(i) = log10(360.0)
	  else
	    sim_width(i)=log10(acos(width)/d2r*2.0)
	  endif

	  this_patchy = patchy + (1.0-patchy)*ran1(seed)
	  sim_width(i) = sim_width(i)*this_patchy

	  law_width = log10(konst / sqrt(period(i)))
	  diff = obs_width(i) - law_width
	  diff = obs_width(i) - sim_width(i)
	  write(33,*)i,logperiod(i),obs_width(i),sim_width(i),
     +               diff
	  if(abs(3.1415-2*alpha-beta) .lt. rho) then
	    write(77,*)period(i),alpha/d2r,beta/d2r
	    ninterpulse = ninterpulse + 1
	  endif
	enddo

	call pgbegin(0,device,1,1)
	call pgscf(2)
	call boxes(0.1,0.9,0.1,0.9,-1.5,1.0,0.3,2.3,'bcnst','bcnst',
     +             'log[period]','log[width]','')
c	call boxes(0.1,0.9,0.1,0.9,-3.0,1.0,0.0,2.5,'bcnst','bcnst',
c     +             'log[period]','log[width]','')
	if(.not. nosim)then
	  call pgsci(5)
	  call pgpt(npsr,logperiod,sim_width)
	endif
	if(.not. noobs)then
	  call pgsci(2)
	  call pgpt(npsr,logperiod,obs_width)
	endif

c	goto 100

	read(*,*)
	call pgpage
	v1 = 0.1
	v3 = 0.1
	do j=1,9
	  v2 = v1+0.27
	  v4 = v3+0.27
	  pmin=-1.50 + (j-1)*0.25
	  pmax=pmin+0.25
	  c=1
	  do i=1,npsr
	    if(logperiod(i).ge.pmin .and. logperiod(i).lt.pmax)then
	      obstmp_width(c) = obs_width(i)
	      simtmp_width(c) = sim_width(i)
	      c=c+1
	    endif
	  enddo
	  c=c-1
	  write(*,*)j,c
	  call dohist(obstmp_width,c,0.0,2.5,ymin,ymax,tgood,
     +                obshist,xhist,10)
	  yymax=ymax
	  call dohist(simtmp_width,c,0.0,2.5,ymin,ymax,tgood,
     +                simhist,xhist,10)
	  if(ymax.gt.yymax)yymax=ymax
	  call boxes(v1,v2,v3,v4,0.0,2.5,0.0,yymax,'bcst','bcst',
     +             '','','')
	  call pgsci(2)
	  call pgbin(10,xhist,obshist,.true.)
	  write(string,'(''N='',i3)')c
	  call pgtext(1.75,yymax*0.9,string)
	  call pgsci(4)
	  call pgbin(10,xhist,simhist,.true.)
	  v1=v1+0.27
	  if(v1.gt.0.75)then
	    v1=0.1
	    v3=v3+0.27
	  endif
	enddo

 100	continue
	call pgend
	write(*,'('' There were '',i8,'' interpulses'')') ninterpulse
	end
