c     library name : minimiz.f   

c     Copyright (c) 2013      K.Vengadesan and N.Gautham

c     This library is free software; you can redistribute it and/or
c     modify it under the terms of the GNU Lesser General Public
c     License as published by the Free Software Foundation; either
c     version 2.1 of the License, or (at your option) any later version.
c
c     This library is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     Lesser General Public License for more details.
c
c     You should have received a copy of the GNU Lesser General Public
c     License along with this library; if not, write to the Free
c     Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
c     02110-1301  USA
c
c
c     contact : n_gautham@hotmail.com  




c	main line program for subroutine frprmn(fletcher - Reeves
c	and Polak - Ribiere algorithm
c	program cg$main
	subroutine minimiz

	include 'mols.par'
        common/crda/x(maxatm,8) 
        common/ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,12)
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common/crdb/y(maxatm,8) 
        common/vectors/iv(maxpar,4)
        common /par/ natom, ntor, nhb, ns, lres
        common/hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
        common/order/n
        common/inp/e_inp(maxstr),p1(maxstr,maxpar)
	common/mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)

	common /ctl/iopt,iff,icint,fseq(maxres)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4

	real fret,ftol,p(maxpar)
	integer iter, n
cc	write(31,*) 'enter no. of atoms'
cc	read *, natom
c	n = 8
cc	write(31,*) 'enter no. of parameters'
c	natom=77
cc	read *, n
	n = ntor
cc	write(31,*) 'enter no. of hbond parameters'
cc	read *, nhb
cc	write(31,*) 'enter no. of structures'
cc	read *, ns

cccccc	call input2
	call getx
        do i=1,natom
          do j=4,8
             y(i,j)=x(i,j)
          enddo          
	enddo
ccccc	call pdb2
	do i = 1,ns
	  do j=1,n
	    p(j) = p1(i,j)
	  enddo
	  call frprmn(p,n,ftol,iter,fret)
	  call outp(p,fret,i)
c	write(6,*) 'Struc. NO:	Energy value	No. of iteration'	
	write(6,*) 'minimized structure no: ',i,fret,iter
c*******Flushing standard output************
	call flush(6)
c*******End of Flushing ********************
	write(31,*) 'minimized structure no: ',i,fret,iter
	enddo
c	call pdbout
	return
	end

c***********************************************************************
        subroutine frprmn(p,n,ftol,iter,fret)
        integer iter, n, nmax, itmax
        real fret, ftol, p(n), eps, func
        external func
        parameter (nmax = 150,itmax=1000,eps=1.e-10)
c uses dfunc, func, linmin
c Given a starting point p that is vector of length n. Fletcher-Reeves
c -Polak-Ribiere minimiZation is performed on a function func, using its
c gradient as calculated by a routine dfunc. The convergence tolerance on
c the function value is input as ftol. Returned quantities are p (the 
c location of the minimum), iter (the number of iterations that were
c performed), and fret (the minimum value of the function). The linmin
c is called to perform line minimization.
        integer its, j
        real dgg,fp,gam,gg,g(nmax),h(nmax),xi(nmax)
	ftol = 1.e-10
c	initializations	
        fp = func(p)
c	print *, n,(p(l),l=1,n),fp
c	print *,fp
        call dfunc(p,xi)
        do  j = 1,n
                g(j) = -xi(j)
                h(j) = g(j)
                xi(j) = h(j)
        enddo 
c	loop over iterations
        do  its = 1, itmax
          iter = its
          call linmin(p,xi,n,fret)
cccc	write(*,*) its,fp,fret
c	next statement is the normal return
          if(2.*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+eps)) return
          fp = func(p)
          call dfunc(p,xi)
ccc	write(31,*)(xi(i1),i1=7,n)
          gg = 0.
          dgg = 0.
          do  j = 1,n
            gg = gg + g(j)**2
c	this statement for Fletcher - Reeves
c           dgg = dgg + xi(j)**2
c	this statement for Polak - Ribiere
            dgg = dgg + (xi(j)+g(j))*xi(j)
          enddo 
c	unlikely, if gradient is exactly zero then we are already done.
          if(gg.eq.0) return
          gam = dgg/gg
          do  j = 1,n
                g(j) = -xi(j)
                h(j) = g(j) + gam*h(j)
                xi(j) = h(j)
          enddo 
        enddo 
        pause 'frprmn maximum iterations exceeded'
        return
        end

c********************************************************************
        subroutine linmin(p,xi,n,fret)
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
        external f1dim,df1dim
c	set up the common block
        ncom = n
        do  j = 1,n
                pcom(j) = p(j)
                xicom(j) = xi(j)
        enddo 
c	intial guess for brackets.
        ax = 0.
        xx = 1.
        call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)

        fret = dbrent(ax,xx,bx,f1dim,df1dim,tol,xmin)
c	construct the vector results to return
        do  j = 1,n
                xi(j) = xmin*xi(j)
                p(j) = p(j) + xi(j)
        enddo 
        return
        end

c**********************************************************************
        function dbrent(ax,bx,cx,f,df,tol,xmin)
        integer itmax
        real dbrent,ax,bx,cx,tol,xmin,df,f,zeps
        external df,f
        parameter (itmax =150, zeps = 1.0e-10)
c   Given function f and its derivative function df, and given a bracketing 
c triplet of abscissas ax, bx, cx (such that bx is between ax and cx, and
c f(bx) is less than both f(ax) and f(cx). this routine isolates the minimum
c to a fractional precision of about tol using a modification of Brents's
c method that uses derivatives. the abscissa of the minimum is returned as xmin,
c and the minimum function value is returned as dbrent, the returned function
c value	  
        integer iter
        real a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,
     &  u,u1,u2,v,w,x,xm
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
c	initialize these d's to an out-of-bracket value.
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
        pause 'dbrent exceeded maximum iteration'
c	arrive here ready to exit with best values
3       xmin = x
        dbrent = fx
        return
        end
c********************************************************************
        subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
        real ax, bx, cx, fa, fb, fc, func, gold, glimit, tiny
        external func
        parameter (gold=1.618034,glimit=100.,tiny=1.e-20)
c	Given a fucntion func, and given distinct intial points ax and bx,
c this routine searches in the downhill direction (defined by the function
c as evaluated at the initial points) and returns new points ax,bx,cx tha
c bracket a minimum of the function. also returned are the function values at
c three points, fa, fb and fc. 
c gold is the default ratio by which successive intervals are magnified;
c glimit is the maximum magnification allowed for a parabolic-fit step.
        real dum,fu,q,r,u,ulim
        fa = func(ax)
        fb = func(bx)
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
        fc = func(cx)
c	do while keep returning here until we bracker.
1       if(fb.ge.fc) then
c	compute u by parabolic extrapolation from a,b,c
           r = (bx-ax)*(fb-fc)
           q = (bx-cx)*(fb-fa)
c	tiny is used to prevent any possible division by zero
          u = bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),tiny),q-r))
c	we would not go farther than this. 
           ulim = bx+glimit*(cx-bx)
c	test various possibilities
           if((bx-u)*(u-cx).gt.0.) then
c	parabolic u is between b and c try it
                fu = func(u)
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
                fu = func(u)
c	parabolic fit is between c and its allowed limit
           else if((cx-u)*(u-ulim).gt.0.) then
                fu = func(u)
                if(fu.lt.fc) then
                  bx = cx
                  cx = u
                  u = cx + gold*(cx-bx)
                  fb = fc
                  fc = fu
                  fu = func(u)
                endif
c	limit parabolic u to maximum allowed value
          else if((u-ulim)*(ulim-cx).ge.0.) then
                u = ulim
                fu = func(u)
c	reject parabolic u, use default magnification
          else
                u = cx + gold*(cx-bx)
                fu = func(u)
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
	function f1dim(x)
	integer nmax
	real f1dim, func,x
	parameter (nmax=150)
c	uses func
c	used by linimin as the function passed to mnbrak and dbrent
	external func
	integer j, ncom
	real pcom(nmax),xicom(nmax),xt(nmax)
	common /f1com/ pcom, xicom, ncom
	do j = 1,ncom
	   xt(j) = pcom(j)+x*xicom(j)
	enddo
	f1dim = func(xt)
	return
	end
c****************************************************************************          
        function df1dim(x)
        integer nmax
        real df1dim, x
        parameter (nmax =150)
c       uses dfunc
c	used by linimin as the function passed to mnbrak and dbrent
        integer j, ncom
        real df(nmax), pcom(nmax), xicom(nmax), xt(nmax)
        common /f1com/ pcom, xicom, ncom
        do  j = 1,ncom
                xt(j) = pcom(j)+x*xicom(j)
        enddo 
        call dfunc(xt,df)
        df1dim=0.
        do  j = 1,ncom
                df1dim = df1dim + df(j)*xicom(j)
        enddo 
        return
        end
c********************************************************************
c********************************************************************
        subroutine getx
	include 'mols.par'
        common/inp/e_inp(maxstr),p1(maxstr,maxpar)
        common /par/ natom, ntor, nhb, ns, lres
	common/mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
        common/order/n
c10	format(1x,8(1x,f5.1),5x,f20.4)
cc	open(unit=2,file='mols.out',status='unknown')
	do i = 1,ns
	e_inp(i) = emo(i)
	do j = 1,n
	p1(i,j) = opmo(i,j)
	enddo
c	print *,(p1(i,j1),j1=1,n),e_inp(i) 
c	read(2,10)(p1(i,j),j=1,n),e_inp(i) 
cc	read(2,*)(p1(i,j),j=1,n),e_inp(i)
	enddo
cc	close(unit=2)
	return
	end
c********************************************************************
c*********************************************************************
	function func(p2)
	include 'mols.par'
        common/order/n
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /ctl/iopt,iff,icint,fseq(maxres)
	dimension p2(maxpar)
c
        do j = 1,n
        p2(j) = mod(p2(j),360.0)
        if(p2(j).le.0) p2(j) = 360.0 + p2(j)
	phi(j) = p2(j)
        enddo
c
	call molgen(p2)
c	call energee2(ej)
c	func = ej
c	func = ampene(1,1)
		if(iff.eq.1) func=ampene(1,1)
		if(iff.eq.2) func=ecpene(1,1)

	return
	end
c*******************************************************************
	subroutine dfunc(p2,g1)
	include 'mols.par'
        common/order/n
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /ctl/iopt,iff,icint,fseq(maxres)	
	dimension p2(maxpar),g1(maxpar)
	real deltaX,deltaE,f1,p2
	deltaX = 0.1
c
        do j = 1,n
        p2(j) = mod(p2(j),360.0)
        if(p2(j).le.0) p2(j) = 360.0 + p2(j)
	phi(j) = p2(j)
        enddo
c
	call molgen(p2)
c	call energee2(ej)
c	f1 = ej
		if(iff.eq.1) f1=ampene(1,1)
		if(iff.eq.2) f1=ecpene(1,1)
ccccc	f1 = ecpene(1,1)
	do j=1,n
	g1(j) = 0.0
	deltaE = 0.0
	ej = 0.0
	p2(j) = p2(j) - deltaX
	phi(j) = p2(j)
	call molgen(p2)
c	call energee2(ej)
c	ej = ampene(1,1)
		if(iff.eq.1) ej=ampene(1,1)
		if(iff.eq.2) ej=ecpene(1,1)
	deltaE = f1-ej
	g1(j) = deltaE/deltaX
	p2(j) = p2(j) + deltaX
	phi(j) = p2(j)
	enddo

c214	format (/,4(2x,2hG(,i1,4h) = ,f10.4))
c	write (6,214) (j,g1(j), j=1,n)

	return
	end
c*******************************************************************
	subroutine outp(p3,ene,i)
	include 'mols.par'
        common/order/n
        common /par/ natom, ntor, nhb, ns, lres
	common/mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	dimension p3(maxpar)
	do k1=1,n
	opmi(i,k1) = p3(k1)
	p3(k1) = p3(k1)-180.0
	enddo 
c	print *,i,(opmi(i,k1),k1=1,n)
	emi(i) =ene
	open(unit=14,file=of2,status='unknown')
c           write(14,*)(p3(j),j=1,n),ene 
c           write(14,10)i,ene,(p3(j),j=1,n)
           write(14,*)i,ene,(p3(j),j=1,n)
c10	format(i4,1x,f15.2,<n>(1x,f6.1))
	if(i.eq.ns) close(unit=14)
        return
        end
c*******************************************************************

