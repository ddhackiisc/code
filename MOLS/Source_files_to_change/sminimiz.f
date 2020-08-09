c	main line program for subroutine frprmn(fletcher - Reeves
c	and Polak - Ribiere algorithm
c	program cg$main
	subroutine sminimiz(kk11)

	include 'mols.par'
cs      common /crda/x(maxatm,8)
        common /ligcrda/xlig(maxatm,8)
        common /ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,12)
	common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
cs      common /crdb/y(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /vectors/iv(maxpar,4)
        common /par/natom, ntor, nhb, ns, lres
        common /hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
        common /order/n
        common /inp/e_inp(maxstr),p1(maxstr,maxpar)
	common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
        common /pout/popmi(maxstr,maxpar),pemi(maxstr)
c	common /ctl/iopt,iff,icint,fseq(maxres)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /pctl/ifff,ioptt
	common /fnam/if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
        common /comment/icomment
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4

	real fret,ftol,p(maxpar)
        real y(maxatm,8)
	integer iter, n,nsm
cc	write(31,*) 'enter no. of atoms'
cc	read *, natom
c	n = 8
cc	write(31,*) 'enter no. of parameters'
c	natom=77
cc	read *, n
c	n = ntor
c	MODIFIED TO INCLUDE TRANS AND ROT PARAMETER
	n = ntor + 6 + np
cc	write(31,*) 'enter no. of hbond parameters'
cc	read *, nhb
cc	write(31,*) 'enter no. of structures'
cc	read *, ns

cs	do i=1,ns

        do i=kk11,kk11
	do j=1,n
	   p1(i,j)=0.0
	enddo
	enddo
cccccc	call input2

        if(icomment.eq.1)then
        do i=kk11,kk11
        do j=1,n
          print *,'smin-opmo',i,j,opmo(i,j)
        enddo
        enddo
        endif


	call sgetx(kk11)

c	print *,'getx completed'
	do i=1,natom
             y(i,1)=ylig(i,1)
             y(i,2)=ylig(i,2)
             y(i,3)=ylig(i,3)
cs      print *,'y-smin',ylig(i,1),ylig(i,2),ylig(i,3)   
        enddo

ccccc	call pdb2
cs	do i = 1,ns
        do i = kk11,kk11
c	print *,'stno & n = ',i,n
	n = ntor + 6 + np
	  do j=1,n
	   p(j) = p1(i,j)
c	  print *,'after-pgetx',p(j)
	  enddo
	
c	  print *,'entering function: ',p
c	  temp = func(p)
c	  stop
cs        print *,'before frprmn, stno & n = ',i,n
cs  	  print *,'before frprmn: ',(p(k),k=1,n)
        
          nsm = n
 
          call sfrprmn(p,nsm,ftol,iter,fret)
        
cs        print *,'after frprmn, stno & nsm = ',i,nsm
cs        print *,'after sfrprmn: ',(p(k),k=1,n)
               
	  call soutp(p,fret,i)

c	write(6,*) 'Struc. NO:	Energy value	No. fo iteration'
cs	write(6,*) 'minimized structure no: ',i,fret,iter
c*******Flushing standard output************
	call flush(6)
c*******End of Flushing ********************
cs	write(31,*) 'minimized structure no: ',i,fret,iter
cs      write(77,*)'MINI(LIGAND+PLP)',i,fret
	enddo
c	call pdbout
 	return
	end

c***********************************************************************
        subroutine sfrprmn(p,nsm,ftol,iter,fret)

        common /comment/icomment
        common /sminimize/itmaxs
        integer iter,nmax,nsm, itmaxs
        real fret, ftol, p(nsm), eps, sfunc
        external sfunc
cs      parameter (nmax = 150,itmax=1000,eps=1.e-10)
cs      parameter (nmax = 150,itmax=5,eps=1.e-10)
        parameter (nmax = 150,eps=1.e-5)
c       parameter (nmax = 150,itmax=1000,eps=1.e-5)
c uses dfunc, func, linmin
c Given a starting point p that is vector of length n. Fletcher-Reeves
c -Polak-Ribiere minimiZation is performed on a function func, using its
c gradient as calculated by a routine dfunc. The convergence tolerance on
c the function value is input as ftol. Returned quantities are p (the
c location of the minimum), iter (the number of iterations tha were
c performed), and fret (the minimum value of the function). The linmin
c is called to perform line minimization.
        integer its, j
        real dgg,sfp,gam,gg,g(nmax),h(nmax),xi(nmax)
        itmax = itmaxs !user-specified
	ftol = 1.e-5
c	ftol = 1.e-2
c	initializations
cs      print *,'before sfunc, nsm = ',nsm

        sfp = sfunc(p)

cs	print *,'after sfunc, nsm = ',nsm
c	print *,'inside frprmn: ',(p(l),l=1,n),sfp
c	print *,sfp
        call sdfunc(p,xi)
                
        do  j = 1,nsm
                g(j) = -xi(j)
                h(j) = g(j)
                xi(j) = h(j)
        enddo
c	loop over iterations
        do  its = 1, itmax
          iter = its
          if(icomment.eq.1) print *,'smin-iteration',its
c	  print *,'before linmin: ',(p(il),il=1,n)
	  call slinmin(p,xi,nsm,fret)
c	  print *,'after linmin: ',(p(il),il=1,n)
cccc	write(*,*) its,sfp,fret
c	next statement is the normal return
          if(2.*abs(fret-sfp).le.ftol*(abs(fret)+abs(sfp)+eps)) return
          sfp = sfunc(p)
          call sdfunc(p,xi)
ccc	write(31,*)(xi(i1),i1=7,n)
          gg = 0.
          dgg = 0.
          do  j = 1,nsm
            gg = gg + g(j)**2
c	this statement for Fletcher - Reeves
c           dgg = dgg + xi(j)**2
c	this statement for Polak - Ribiere
            dgg = dgg + (xi(j)+g(j))*xi(j)
          enddo
c	unlikely, if gradient is exactly zero then we are already done.
          if(gg.eq.0) return
          gam = dgg/gg
          do  j = 1,nsm
                g(j) = -xi(j)
                h(j) = g(j) + gam*h(j)
                xi(j) = h(j)
          enddo
        enddo
c       print *,'end-sfrprmn-n',nsm
c       pause 'frprmn maximum iterations exceeded'
        return
        end

c********************************************************************
        subroutine slinmin(p,xi,n,fret)
	include 'mols.par'
	common /par/ natom, ntor, nhb, ns, lres
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
        integer n, nmax
        real fret, p(n), xi(n), tol
        parameter (nmax =150, tol = 1.e-4)
c uses dbrent, f1dim, df1dim, mnbrak
c given an n-dimensional poin p(1:n) and an n-dimensional direction xi(1:n)
c moves and resets p to where the function func(p) takes on a minimum
c along the direction xi from p, and replaces xi by the actual vector
c displacement that p was moved. also returns as fret the value of func at
c the returned location p. this is actually all accomplished by calling
c the routines mnbrak and brent.

        integer j, ncom
        real ax,bx,fa,fb,fx,xmin,xx,pcom(nmax),xicom(nmax),brent
        common /f1com/ pcom,xicom,ncom
        external sf1dim,sdf1dim
c	set up the common block
	n = ntor + 6 + np
        ncom = n
c	print *,'line min n=',n
c	stop

        do  j = 1,n
                pcom(j) = p(j)
                xicom(j) = xi(j)
        enddo
c	intial guess for brackets.
        ax = 0.
        xx = 1.
        call smnbrak(ax,xx,bx,fa,fx,fb,sf1dim)
c	print *,'after mnbrak'

        fret = sdbrent(ax,xx,bx,sf1dim,sdf1dim,tol,xmin)
c	print *,'after dbrent'
c 	print *,(p(il),il=1,n)
c	construct the vector results to return
        do  j = 1,n
                xi(j) = xmin*xi(j)
                p(j) = p(j) + xi(j)
        enddo
c	print *,(xi(il),il=1,n)
        return
        end

c**********************************************************************
        function sdbrent(ax,bx,cx,f,df,tol,xmin)
        integer itmax
        real pdbrent,ax,bx,cx,tol,xmin,df,f,zeps
        external df,f
        parameter (itmax =150, zeps = 1.0e-10)
c   Given function f and its derivative function df, and given a bracketing
c triplet of abscissas ax, bx, cx (such that bx is between ax and cx, and
c f(bx) is less than both f(ax) and f(cx). this routine isolates the minimum
c to a fractional precision of about tol using a modification of Brents''s
c method that uses derivatives. the abscissa of the minimum is returned as xmin,
c and the minimum function value is returned as dbrent, the returned function
c value
        integer iter
        real a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,
     &u,u1,u2,v,w,x,xm
        logical ok1, ok2

c	a and b must be in ascending order, though the input abscissas
c	need not be.
        a = min(ax,cx)
        b = max(ax,cx)
c	intializations..
        v = bx
        w = v
        x = v
        e = 0.
        fx = f(x)
        fv = fx
        fw = fx
c	all our housekeeping chores are doubled by the necessity of moving
c	derivative values around as well as function values.
        dx = df(x)
        dv = dx
        dw = dx
c	main program loop
        do  iter = 1, itmax
           xm = 0.5*(a+b)
           tol1 = tol*abs(x) + zeps
           tol2 = 2.*tol1
c	test for done here
           if(abs(x-xm).le.(tol2-.5*(b-a))) go to 3
c	construct a trial parabolic fit.
           if(abs(e).gt.tol1) then
c	initialize these d''s to an out-of-bracket value.
             d1 = 2.*(b-a)
             d2 = d1
c	secant method with one point.
             if(dw.ne.dx) d1 = (w-x)*dx/(dx-dw)
c	and the other.
             if(dv.ne.dx) d2 = (v-x)*dx/(dx-dv)
c	which of these two estimates of d shall we take? we will insist that
c	they be within the bracket, and on the side pointed to by the
c	derivative at x
             u1 = x + d1
             u2 = x + d2
             ok1 = ((a-u1)*(u1-b).gt.0.).and.(dx*d1.le.0.)
             ok2 = ((a-u2)*(u2-b).gt.0.).and.(dx*d2.le.0.)
c	movement on the step before last
             olde = e
             e = d
c	take only an acceptable d, and if both are acceptable, then take
c	 the smallest one
             if(.not.(ok1.or.ok2)) then
                go to 1
             else if(ok1.and.ok2) then
                if(abs(d1).lt.abs(d2)) then
                  d = d1
                else
                  d = d2
                endif
            else if (ok1) then
                d = d1
            else
                d = d2
            endif
            if(abs(d).gt.abs(0.5*olde)) go to 1
            u = x + d
c	skip over the golden section step
            if(u-a.lt.tol2.or.b-u.lt.tol2) d = sign(tol1,xm-x)
                go to 2
            endif
c	decide which segment by the sign of the derivative
1           if(dx.ge.0.) then
                e = a - x
            else
                e = b - x
            endif
            d = 0.5*e
2           if(abs(d).ge.tol1) then
                u = x + d
                fu = f(u)
            else
                u = x + sign(tol1,d)
                fu = f(u)
c	if the minimum step in the downhill direction takes us uphill
c	then we are done.
                if(fu.gt.fx) go to 3
            endif
c	now all the housekeeping sigh.
            du = df(u)
            if(fu.le.fx) then
		if(u.ge.x) then
                a = x
            else
                b = x
            endif
            v = w
            fv = fw
            dv = dw
            w = x
            fw = fx
            dw = dx
            x = u
            fx = fu
            dx = du
         else
            if(u.lt.x) then
                a = u
            else
                b = u
            endif
            if(fu.le.fw.or.w.eq.x) then
                v = w
                fv = fw
                dv = dw
                w = u
                fw = fu
                dw = du
           else if(fu.le.fv.or.v.eq.x.or.v.eq.w) then
                v = u
                fv = fu
                dv = du
           endif
         endif
        enddo
c       pause 'dbrent exceeded maximum iteration'
c	arrive here ready to exit with best values
3       xmin = x
        sdbrent = fx
        return
        end
c********************************************************************
        subroutine smnbrak(ax,bx,cx,fa,fb,fc,sfunc)
        real ax, bx, cx, fa, fb, fc, func, gold, glimit, tiny
        external sfunc
        parameter (gold=1.618034,glimit=100.,tiny=1.e-20)
c        parameter (gold=1.618034,glimit=1000.,tiny=1.e-5)
c	Given a fucntion func, and given distinct intial points ax and bx,
c this routine searches in the downhill direction (defined by the function
c as evaluated at the initial points) and returns new points ax,bx,cx tha
c bracket a minimum of the function. also returned are the function values at
c three points, fa, fb and fc.
c gold is the default ratio by which successive intervals are magnified;
c glimit is the maximum magnification allowed for a parabolic-fit step.
        real dum,fu,q,r,u,ulim
        fa = sfunc(ax)
        fb = sfunc(bx)
c	switch roles of a and b so that we can go downhill in the direction
c from a to b
        if(fb.gt.fa) then
           dum = ax
           ax = bx
           bx = dum
           dum = fb
           fb = fa
           fa = dum
        endif
c	first guess for c
        cx = bx+gold*(bx-ax)
        fc = sfunc(cx)
c	do while keep returning here until we bracker.
1       if(fb.ge.fc) then
c	compute u by parabolic extrapolation from a,b,c
           r = (bx-ax)*(fb-fc)
           q = (bx-cx)*(fb-fa)
c	tiny is used to prevent any possible division by zero
          u = bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),tiny),q-r))
c	we won''t go farther than this.
           ulim = bx+glimit*(cx-bx)
c	test various possibilities
           if((bx-u)*(u-cx).gt.0.) then
c	parabolic u is between b and c try it
                fu = sfunc(u)
c	got a minimum between b and c
                if(fu.lt.fc) then
                   ax = bx
                   fa = fb
                   bx = u
                   fb = fu
                   return
c	got a minimum between a and u
                else if(fu.gt.fb) then
                   cx = u
                   fc = fu
                   return
                endif
c	parabolic fit was no use. use default magnification.
                u = cx + gold*(cx-bx)
                fu = sfunc(u)
c	parabolic fit is between c and its allowed limit
           else if((cx-u)*(u-ulim).gt.0.) then
                fu = sfunc(u)
                if(fu.lt.fc) then
                  bx = cx
                  cx = u
                  u = cx + gold*(cx-bx)
                  fb = fc
                  fc = fu
                  fu = sfunc(u)
                endif
c	limit parabolic u to maximum allowed value
          else if((u-ulim)*(ulim-cx).ge.0.) then
                u = ulim
                fu = sfunc(u)
c	reject parabolic u, use default magnification
          else
                u = cx + gold*(cx-bx)
                fu = sfunc(u)
          endif
c	eliminate oldest point and continue
          ax = bx
          bx = cx
          cx = u
          fa = fb
          fb = fc
          fc = fu
          go to 1
        endif
        return
        end
c****************************************************************************
	function sf1dim(x)
	integer nmax
	real sf1dim,sfunc,x
	parameter (nmax=150)
c	uses func
c	used by linimin as the function passed to mnbrak and dbrent
	external sfunc
	integer j, ncom
	real pcom(nmax),xicom(nmax),xt(nmax)
	common /f1com/ pcom, xicom, ncom
c       print *,'sf1dim-ncom',ncom
	do j = 1,ncom
	   xt(j) = pcom(j)+x*xicom(j)
	enddo
	sf1dim = sfunc(xt)
	return
	end
c****************************************************************************
        function sdf1dim(x)
        integer nmax
        real sdf1dim, x
        parameter (nmax =150)
c       uses dfunc
c	used by linimin as the function passed to mnbrak and dbrent
        integer j, ncom
        real df(nmax), pcom(nmax), xicom(nmax), xt(nmax)
        common /f1com/ pcom, xicom, ncom
        do  j = 1,ncom
                xt(j) = pcom(j)+x*xicom(j)
        enddo
c	print *,ncom
c	print *,'from dfldim ',(xt(j),j=1,ncom)
        call sdfunc(xt,df)
c	print *,'after dfunc ',(xt(j),j=1,ncom)
        sdf1dim=0.
        do  j = 1,ncom
                sdf1dim = sdf1dim + df(j)*xicom(j)
        enddo
        return
        end
c********************************************************************
c********************************************************************
        subroutine sgetx(kk11)
	include 'mols.par'
        common/inp/e_inp(maxstr),p1(maxstr,maxpar)
        common/par/natom, ntor, nhb, ns, lres
        common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
	common/mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
	common/recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
        common/order/n
        common/comment/icomment
c10	format(1x,8(1x,f5.1),5x,f20.4)
cc	open(unit=2,file='mols.out',status='unknown')
	if(icomment.eq.1) print *,'parameters: ',n,np
c	stop

        if(icomment.eq.1)then
        do i=kk11,kk11
        do j=1,n
          print *,'sgetx-opmo',i,j,opmo(i,j)
        enddo
        enddo
        endif

cs      IF(ifopt.eq.2)THEN
	DO i = kk11,kk11
	e_inp(i)=emo(i)
        if(icomment.eq.1)print *,'e_inp',i,e_inp(i)
	do j = 1,n
	p1(i,j) = opmo(i,j)
	enddo

css        do j = 1,n-np-6
css        p1(i,j) = opmo(i,j)
css        enddo

css        do j = n-np-5,n-np-3
css        tt1 = opmo(i,j)
css        if(tt1.eq.0.0) then
css        tt2 = 0.0
css        else
c       tt2 = (tt1/10.0)*.28
css        tt2 = (tt1/10.0)*0.13888889
c       tt2 = (tt1/10.0)*0.55555556
css        endif
c       tt2 = tt2-5.0
css        tt2 = tt2-2.5
c       tt2 = tt2-10.0
css        p1(i,j) = tt2
css        print *,'p1',p1(i,j)
css        enddo

css        do j = n-np-2,n
css        p1(i,j) = opmo(i,j)
css        print *,'p1',p1(i,j)
css        enddo
        
        ENDDO
cs      ENDIF
        
css        do i=kk11,kk11
css        do j=1,n
css          print *,'opmo,p1',i,j,opmo(i,j),p1(i,j)
css        enddo
css        enddo

cs	do j = n-np-5,n-np-3
cs	p1(i,j)=((opmo(i,j)*10.0)/0.13888889)
cs	enddo

cs	do j = n-np-2,n
cs	p1(i,j)=opmo(i,j)
cs	enddo

c	read(2,10)(p1(i,j),j=1,n),e_inp(i)
cc	read(2,*)(p1(i,j),j=1,n),e_inp(i)

cc	close(unit=2)
	return
	end
c********************************************************************
c*********************************************************************
	function sfunc(p2)
	include 'mols.par'
cs      common /order/n
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
c	common /ctl/iopt,iff,icint,fseq(maxres)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
        common /pctl/ifff,ioptt
	common /par/natom, ntor, nhb, ns, lres
	common /scang/frange,rang
	common /rctl/iscopt
        common /comment/icomment
	dimension p2(maxpar)
	integer n
	real inter_ene,intra_ene

        if(icomment.eq.1)print *,'ifff',ifff
c	n = 0
        n = ntor + 6 + np

	do j = 1,n-np
        p2(j) = mod(p2(j),360.0)
        if(p2(j).le.0) p2(j) = 360.0 + p2(j)
	phi(j) = p2(j)
c	print *,'phi-pfunc',j,phi(j)
	enddo
c
c	print *,'inside function: ',(p2(j),j=1,n)

	do j = n-np-5,n-np-3
cs	p2(j) = mod(p2(j),5.0)
cs	if(p2(j).le.0) p2(j) = 5.0 + p2(j) 
	tt1 = p2(j)
c-------sam added--------------------------------
cs	tt3 = ((tt1*10.0)/0.13888889)
cs	if (tt3.gt.360.0)then
cs	tt1 = mod(tt1,360.0)
cs	if(tt1.le.0) tt1 = 360.0 + tt1
cs	tt1 = ((tt1/10.0)*0.13888889)
cs	endif
c------------------------------------------------
	if(tt1.eq.0.0) then
	tt2 = 0.0
	else
c	tt2 = (tt1/10.0)*.28
	tt2 = (tt1/10.0)*0.13888889
cs	tt2 = tt1
c	tt2 = (tt1/10.0)*0.55555556
	endif
c	tt2 = tt2-5.0
	tt2 = tt2-2.5
c	tt2 = tt2-10.0
	phi(j) = tt2
        if(icomment.eq.1) print *,'phi-sfunc',j,phi(j)
	enddo
c--------------------------------------------------
cs	do j = n-np-2,n-np ! sam added for rotation
cs	p2(j) = mod(p2(j),360.0)
cs      if(p2(j).le.0) p2(j) = 360.0 + p2(j)
cs      phi(j) = p2(j)
cs      if(icomment.eq.1)print *,'phi-pfunc',j,phi(j)
cs	enddo
c--------------------------------------------------
        do j = n-np+1,n
        if(iscopt.eq.1)then
        p2(j) = mod(p2(j),360.0)
        if(p2(j).le.0) p2(j) = 360.0 + p2(j)
        phi(j) = p2(j)
        elseif(iscopt.eq.0)then
        p2(j) = mod(p2(j),frange)
c       p2(j) = mod(p2(j),15.0)
        if(p2(j).le.-(frange))p2(j) = p2(j)+frange
c       if(p2(j).le.-15.0)p2(j) = p2(j)+15.0
        endif
        phi(j) = p2(j)
c        print *,'func-phi',j,phi(j)
        enddo

c	print *,'after converting: ',(phi(j),j=1,n)
c	stop
        if(icomment.eq.1)then
        do ik=1,n
        print *,'phi-smin',ik,phi(ik)
        enddo
        endif
c	call molgen(p2)
c	MODIFIED P2->PHI
cs	call ppmolgen(phi,2)
	call smolgen(2,2,phi,2,ifopt,eflex)
c	call energee2(ej)
c	func = ej
c	func = ampene(1,1)
c		if(iff.eq.1) func=ampene(1,1)
c		if(iff.eq.2) func=ecpene(1,1)
c------------
		if(ifff.eq.1) then
c		 intra_ene = pampene(i,j)
                 intra_ene = rmmffene(i,j)
c		 call flexall(2,2,phi,e_flex,2)
		 call eplp(plpe,hb,steric,2)
		 inter_ene = plpe
		endif

		if(ifff.eq.2) then
		 intra_ene = rgaffene(i,j)
c		 call flexall(2,2,phi,e_flex,2)
		 call eplp(plpe,hb,steric,2)
		 inter_ene = plpe
		endif

		if(ifff.eq.3) then
c		 intra_ene = ecpene(i,j)
c		  print *,'Intra ene', intra_ene
c		 print *,'about to calculate'
c		 inter_ene = autointer()
c		  print *,'Inter ene', inter_ene
		endif

		sfunc=intra_ene+inter_ene+eflex
	if(icomment.eq.1) print *,'pep-plp-sfunc',intra_ene,inter_ene,
     &  eflex,sfunc

c------------

c		if(iff.eq.1) tfunc=ampene(1,1)
c		if(iff.eq.2) tfunc=ecpene(1,1)
c	call eplp(plpe,hb,steric,2)
c	print *,'function value: ecepp & plp ',tfunc,plpe
c	tfunc = tfunc + plpe
c	func = tfunc
	return
	end
c*******************************************************************
	subroutine sdfunc(p2,g1)
	include 'mols.par'
c       common/order/n
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
c	common /ctl/iopt,iff,icint,fseq(maxres)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /pctl/ifff,ioptt
	common /comment/icomment
	common /par/ natom, ntor, nhb, ns, lres
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
	common /scang/frange,rang
        common /rctl/iscopt

	dimension p2(maxpar),g1(maxpar)
	real deltaX,deltaE,f1,p2,temp,plp,inter_ene,intra_ene
	integer n
	n = 0
	deltaX = 0.1
	temp = 0.0
	plp = 0.0
	n = ntor + 6 + np

        do j = 1,n-np
        p2(j) = mod(p2(j),360.0)
        if(p2(j).le.0) p2(j) = 360.0 + p2(j)
	phi(j) = p2(j)
c	print *,'pdfunc-phi',j,phi(j)
	enddo

c	print *,'from dfunc ',(phi(j),j=1,n)
c	print *,'ntor=',ntor
c	stop
c
	do j = n-np-5,n-np-3
cs 	p2(j) = mod(p2(j),5.0)
cs      if(p2(j).le.0) p2(j) = 5.0 + p2(j)
	tt1 = p2(j)
c-------sam added--------------------------------
cs        tt3 = ((tt1*10.0)/0.13888889)
cs        if (tt3.gt.360.0)then
cs        tt1 = mod(tt1,360.0)
cs        if(tt1.le.0) tt1 = 360.0 + tt1
cs        tt1 = ((tt1/10.0)*0.13888889)
cs        endif
c------------------------------------------------
	if(tt1.eq.0.0) then
	tt2 = 0.0
	else
c	tt2 = (tt1/10.0)*.28
	tt2 = (tt1/10.0)*0.13888889
cs	tt2 = tt1
c	tt2 = (tt1/10.0)*0.55555556
	endif
c	tt2 = tt2-5.0
	tt2 = tt2-2.5
c	tt2 = tt2-10.0
	phi(j) = tt2
c	print *,'pdfunc-phi',j,phi(j)
	enddo

cs	do j = n-np-2,n-np
cs      p2(j) = mod(p2(j),360.0)
cs      if(p2(j).le.0) p2(j) = 360.0 + p2(j)
cs      phi(j) = p2(j)
c       print *,'pdfunc-phi',j,phi(j)
cs      enddo

	do j = n-np+1,n
        if(iscopt.eq.1)then
        p2(j) = mod(p2(j),360.0)
        if(p2(j).le.0) p2(j) = 360.0 + p2(j)
        phi(j) = p2(j)
        elseif(iscopt.eq.0)then
        p2(j) = mod(p2(j),frange)
c       p2(j) = mod(p2(j),15.0)
        if(p2(j).le.-(frange))p2(j) = p2(j)+frange
c       if(p2(j).le.-15.0)p2(j) = p2(j)+15.0
        phi(j) = p2(j)
c	print *,'dfunc-phi',phi(j)
        endif
        enddo


c	call molgen(p2)
c	MODIFIED P2->PHI
cs	call ppmolgen(phi,2)
	call smolgen(2,2,phi,2,ifopt,eflex)
c	print *,'from dfunc ',(phi(j),j=1,n)
c	print *,'initial p2 set: ',phi
c------------
		if(ifff.eq.1) then
cd		 intra_ene = pampene(i,j)
                 intra_ene = rmmffene(i,j)
cs		 call flexall(2,2,phi,e_flex,2)
		 call eplp(plpe,hb,steric,2)
		 inter_ene = plpe
		endif

		if(ifff.eq.2) then
		 intra_ene = rgaffene(i,j)
cs		 call flexall(2,2,phi,e_flex,2)
		 call eplp(plpe,hb,steric,2)
		 inter_ene = plpe
		endif

		if(ifff.eq.3) then
c		 intra_ene = ecpene(i,j)
c		  print *,'Intra ene', intra_ene
c		 print *,'about to calculate'
c		 inter_ene = autointer()

c		  print *,'Inter ene', inter_ene
		endif
c		func = intra_ene + inter_ene
c------------
c	if(iff.eq.1) f1=ampene(1,1)
c	if(iff.eq.2) f1=ecpene(1,1)
c	call eplp(plp,hb,steric,2)
c	print *,'initial f1 and plpe: ',f1,plp
	f1 = intra_ene + inter_ene + eflex
c	f1 = inter_ene
c	print *,'total: ',f1
c	stop

c	do j=1,n
	do j=1,n-np-6
	plp = 0.0
	g1(j) = 0.0
	deltaE = 0.0
	ej = 0.0
	p2(j) = p2(j) - deltaX
	phi(j) = p2(j)
c	print *,'altered p2 set: ',(phi(ij),ij=1,n)
c	stop
cs	call ppmolgen(phi,2)
	call smolgen(2,2,phi,2,ifopt,eflex)
c	print *,'altered p2 set: ',(phi(ij),ij=1,n)
c------------
		if(ifff.eq.1) then
cd		 intra_ene = pampene(i,j)
                 intra_ene = rmmffene(i,j)
cs		 call flexall(2,2,phi,e_flex,2)
		 call eplp(plpe,hb,steric,2)
		 inter_ene = plpe
		endif

		if(ifff.eq.2) then
		 intra_ene = rgaffene(i,j)
cs		 call flexall(2,2,phi,e_flex,2)
		 call eplp(plpe,hb,steric,2)
		 inter_ene = plpe
		endif

c		func = intra_ene + inter_ene
c------------
	ej = intra_ene + inter_ene + eflex

	deltaE = f1-ej
	g1(j) = deltaE/deltaX
	p2(j) = p2(j) + deltaX
	phi(j) = p2(j)
	enddo

c	deltaX = 0.3
	do j=n-np-5,n-np-3
	plp = 0.0
	g1(j) = 0.0
	deltaE = 0.0
	ej = 0.0
	temp = phi(j)
c---------------------------------------------
cs      tt1 = p2(j)
cs      tt3 = ((tt1*10.0)/0.13888889)
cs      if (tt3.gt.360.0)then
cs      tt1 = mod(tt1,360.0)
cs      if(tt1.le.0) tt1 = 360.0 + tt1
cs      tt1 = ((tt1/10.0)*0.13888889)
cs      endif
c--------------------------------------------
	p2(j) = p2(j) - deltaX
c--------------------------------------------
cs	p2(j) = mod(p2(j),5.0)
cs      if(p2(j).le.0) p2(j) = 5.0 + p2(j)
	tt1 = p2(j)
c--------------------------------------------
	if(tt1.eq.0.0) then
	tt2 = 0.0
	else
c	tt2 = (tt1/10.0)*.28
	tt2 = (tt1/10.0)*0.13888889
cs	tt2 = tt1
c	tt2 = (tt1/10.0)*0.55555556
	endif
c	tt2 = tt2-5.0
	tt2 = tt2-2.5
c	tt2 = tt2-10.0
	phi(j) = tt2

cs	call ppmolgen(phi,2)
	call smolgen(2,2,phi,2,ifopt,eflex)
c	print *,'altered phi set: ',phi
c------------
		if(ifff.eq.1) then
cd 		 intra_ene = pampene(i,j)
                 intra_ene = rmmffene(i,j)
cs		 call flexall(2,2,phi,e_flex,2)
		 call eplp(plpe,hb,steric,2)
		 inter_ene = plpe
		endif

		if(ifff.eq.2) then
		 intra_ene = rgaffene(i,j)
cs		 call flexall(2,2,phi,e_flex,2)
		 call eplp(plpe,hb,steric,2)
		 inter_ene = plpe
		endif

c		func = intra_ene + inter_ene
c------------
c	if(iff.eq.1) ej=ampene(1,1)
c	if(iff.eq.2) ej=ecpene(1,1)
c	call eplp(plp,hb,steric,2)
c	print *,'altered f1 and plpe: ',ej,plp
	ej = intra_ene + inter_ene + eflex
c	print *,'total: ',ej
	deltaE = f1-ej
	g1(j) = deltaE/deltaX
	p2(j) = p2(j) + deltaX
	phi(j) = temp
	enddo

c	deltaX = 0.1
	do j=n-np-2,n-np
	plp = 0.0
	g1(j) = 0.0
	deltaE = 0.0
	ej = 0.0
	p2(j) = p2(j) - deltaX
	phi(j) = p2(j)
c	print *,'altered phi set: ',phi
cs	call ppmolgen(phi,2)
	call smolgen(2,2,phi,2,ifopt,eflex)
c------------
		if(ifff.eq.1) then
cd		 intra_ene = pampene(i,j)
                 intra_ene = rmmffene(i,j)
cs		 call flexall(2,2,phi,e_flex,2)
		 call eplp(plpe,hb,steric,2)
		 inter_ene = plpe
		endif

		if(ifff.eq.2) then
		 intra_ene = rgaffene(i,j)
cs		 call flexall(2,2,phi,e_flex,2)
		 call eplp(plpe,hb,steric,2)
		 inter_ene = plpe
		endif

c		func = intra_ene + inter_ene
c------------
	ej = intra_ene + inter_ene + eflex
	deltaE = f1-ej
	g1(j) = deltaE/deltaX
	p2(j) = p2(j) + deltaX
	phi(j) = p2(j)
	enddo

        do j=n-np+1,n
        plp = 0.0
        g1(j) = 0.0
        deltaE = 0.0
        ej = 0.0
        p2(j) = p2(j) - deltaX
        phi(j) = p2(j)
c       print *,'altered phi set: ',phi
cs      call ppmolgen(phi,2)
	call smolgen(2,2,phi,2,ifopt,eflex)
c------------
                if(ifff.eq.1) then
cd               intra_ene = pampene(i,j)
                 intra_ene = rmmffene(i,j)
cs               call flexall(2,2,phi,e_flex,2)
                 call eplp(plpe,hb,steric,2)
                 inter_ene = plpe
                endif

                if(ifff.eq.2) then
                 intra_ene = rgaffene(i,j)
cs               call flexall(2,2,phi,e_flex,2)
                 call eplp(plpe,hb,steric,2)
                 inter_ene = plpe
                endif

c               func = intra_ene + inter_ene
c------------
        ej = intra_ene + inter_ene + eflex
        deltaE = f1-ej
        g1(j) = deltaE/deltaX
        p2(j) = p2(j) + deltaX
        phi(j) = p2(j)
	if(icomment.eq.1)print *,'phi-dfunc1',phi(j)
        enddo

c	print *,'pep-plp-pdfunc',intra_ene,inter_ene,eflex,ej
c214	format (/,4(2x,2hG(,i1,4h) = ,f10.4))
c	write (6,214) (j,g1(j), j=1,n)
	if(icomment.eq.1)print *,'g1: ',(g1(j),j=1,n)
c	stop
	return
	end
c*******************************************************************
	subroutine soutp(p3,ene,i)
	include 'mols.par'
c       common/order/n
        common /par/ natom, ntor, nhb, ns, lres,ie
	common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
        common /pout/popmi(maxstr,maxpar),pemi(maxstr)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
        common /comment/icomment        

	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	dimension p3(maxpar)
	integer n
	n = ntor + 6 + np

c	print *,'from mini output: '
	do k1=1,n
        opmi(i,k1) = p3(k1)
c	p3(k1) = p3(k1)-180.0
	enddo
	if(icomment.eq.1) print *,'opmi-pmin',i,(opmi(i,k1),k1=1,n)
        emi(i)=ene 
! 'pemi' added to avoid confusion with receptor 'emi'
	if(icomment.eq.1) print *,'pep-plp-pro-min-ene',emi(i)
	open(unit=14,file=of2,status='unknown')
c           write(14,*)(p3(j),j=1,n),ene
c           write(14,10)i,ene,(p3(j),j=1,n)
           write(14,*)i,ene,(p3(j),j=1,n)
c10	format(i4,1x,f15.2,<n>(1x,f6.1))
	if(i.eq.ie) close(unit=14)
        return
        end
c*******************************************************************
        subroutine smolgen(lx,ly,phi,tst,ifopt,eflex)

        include 'mols.par'
cd      common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
cs      common/crda/x(maxatm,8)
cs      common/crdb/y(maxatm,8)
cd      common/pepcrda/xpep(maxatm,8)
cd      common/pepcrdb/ypep(maxatm,8)
cd      common /par/ natom, ntor, nhb, ns, lres
cd      common/vectors/iv(maxpar,4)
cd      common /dcom/ ipatom,ax,ay,az,elen,gsize,gspace,pfile,lfile,
cd   &path,mname,interpol
cd      common/cen/bx,by,bz
cd      common/order/nn,mm

        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /par/ natom, ntor, nhb, ns, lres
        common /vectors/iv(maxpar,4)
        common /order/nn,mm
        common /getrot/inrot,irotb(maxi,2),el,ilsrot(maxi,maxi),
     &  iatrot(maxi),rx(100,3),ire(maxi,maxi),ind(maxi,maxi),big
        common /left/ le(maxi,maxi),innd(maxi),ielenum(maxi),
     &  bonum(maxi)
        common /cen/bx,by,bz
        common /recep/np,nres
        common /comment/icomment


cd        dimension x_one(maxatm,3),x_two(maxatm,3)
cd        dimension phi(maxpar)
cd        integer tst,nn
cd        real rotang,theta,psi
cd        real cx,cy,cz,bx,by,bz
c         real gsize,gspace,elen
cd        character pfile*128,lfile*128,path*128,mname*128

        dimension x_one(maxatm,3),x_two(maxatm,3)
        dimension phi(maxpar)
        integer tst,nn,ci
        real rotang,theta,psi
        real cx,cy,cz,bx,by,bz

        if(icomment.eq.1.and.lx.eq.3.and.ly.eq.3)then
        do ii=1,ntor+6+np
        print *,'phi-final',ii,phi(ii)
        enddo
        endif

c	print *,'ppmolgen',bx,by,bz
cd        nn = ntor
cd        do k=1,natom
cd           do ki=1,3
cd             x_one(k,ki)=xpep(k,ki)
cd           enddo
cd        enddo
cd        do if=1,nn

        nn=ntor
        do k=1,natom
           do ki=1,3
             x_one(k,ki)=rx(k,ki)
           enddo
        enddo
c------------------------------------------------------------------
        if(ntor.eq.0)then
        do k=1,natom
           do ki=1,3
             ylig(k,ki)=rx(k,ki)
           enddo
        enddo
        goto 334   
        endif     
c------------------------------------------------------------------
C###################### PHI ALL####################################

cd        call pelemen(x_one(iv(if,1),1),x_one(iv(if,1),2),
cd     &              x_one(iv(if,1),3),
cd     &              x_one(iv(if,2),1),x_one(iv(if,2),2),
cd     &              x_one(iv(if,2),3),
cd     &              el,em,en)

cd        do k=1,iv(if,3)-1
cd           do ki=1,3
cd             x_two(k,ki)=x_one(k,ki)
cd           enddo
cd        enddo

cd        do k=iv(if,3),iv(if,4)
cd           xin=x_one(k,1)-x_one(iv(if,1),1)
cd           yin=x_one(k,2)-x_one(iv(if,1),2)
cd           zin=x_one(k,3)-x_one(iv(if,1),3)
cd           call protor(el,em,en,phi(if),xin,yin,zin,
cd     &                xout,yout,zout)
cd           x_two(k,1)=xout+x_one(iv(if,1),1)
cd           x_two(k,2)=yout+x_one(iv(if,1),2)
cd           x_two(k,3)=zout+x_one(iv(if,1),3)
cd        enddo

cd        do k=iv(if,4)+1,natom
cd           do ki=1,3
cd             x_two(k,ki)=x_one(k,ki)
cd           enddo
cd        enddo

cd        do k=1,natom
cd           do ki=1,3
cd              x_one(k,ki)=x_two(k,ki)
cd           enddo
cd        enddo
C####################################################################

cd        enddo

cd        do k=1,natom
cd           do ki=1,3
cd           ypep(k,ki)=x_two(k,ki)
cd           enddo
cd        enddo

        do if=1,nn
C###################### PHI ALL ####################################

        call elemen(x_one(irotb(if,1),1),x_one(irotb(if,1),2),
     &              x_one(irotb(if,1),3),
     &              x_one(irotb(if,2),1),x_one(irotb(if,2),2),
     &              x_one(irotb(if,2),3),
     &              el,em,en)


50      format(i4,1x,3f8.3)
        do j=1,ielenum(if) !ielenum = total number of ligand atoms
        k=le(if,j)
           do ki=1,3
           x_two(k,ki)=x_one(k,ki)
           enddo
        enddo
c------------------------------------------------------
c------------------------------------------------------
         do j=1,iatrot(if)
          k=ilsrot(if,j)

           xin=x_one(k,1)-x_one(irotb(if,1),1)
           yin=x_one(k,2)-x_one(irotb(if,1),2)
           zin=x_one(k,3)-x_one(irotb(if,1),3)
           call rotor(el,em,en,phi(if),xin,yin,zin,
     &                xout,yout,zout)
           x_two(k,1)=xout+x_one(irotb(if,1),1)
           x_two(k,2)=yout+x_one(irotb(if,1),2)
           x_two(k,3)=zout+x_one(irotb(if,1),3)

        enddo

        do k=1,natom
           do ki=1,3
              x_one(k,ki)=x_two(k,ki)
           enddo
        enddo
C####################################################################

        enddo


c***********translation and rotation for ligand postioning***********************

cd        call protate(bx,by,bz,phi(nn+4),phi(nn+5),phi(nn+6),natom)
cd        call ptranslate(bx,by,bz,phi(nn+1),phi(nn+2),phi(nn+3),natom)

cd	if(ifopt.eq.2)then
cd        call flexall(2,2,phi,eflex,1)
cd        endif

        do k=1,natom
           do ki=1,3
             ylig(k,ki)=x_two(k,ki)
           enddo
        enddo

334     call rotate(bx,by,bz,phi(nn+4),phi(nn+5),phi(nn+6),natom)
        call translate(bx,by,bz,phi(nn+1),phi(nn+2),phi(nn+3),natom)
        if(ifopt.eq.2)then
        call flexall(lx,ly,phi,eflex,1)
        else
        eflex=0.0
        endif

        if(icomment.eq.1)then
        if(lx.eq.3.and.ly.eq.3)then
        do k=1,natom
        print *,'ylig',ylig(k,1),ylig(k,2),ylig(k,3)
        enddo
        endif
        endif

c*******************************************************************************
        return
        end
