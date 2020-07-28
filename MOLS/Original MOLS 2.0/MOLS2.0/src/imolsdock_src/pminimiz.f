c-----------------------------------------------------------------------------
c     library name : pminimiz.f   

c     Copyright (c) 2020 P.Arun Prasad, S.Nehru Viji, D.Sam Paul and N.Gautham

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
c-----------------------------------------------------------------------------


c	main line program for subroutine frprmn(fletcher - Reeves
c	and Polak - Ribiere algorithm
c	program cg$main
	subroutine pminimiz(kk11)

	include 'mols.par'

        common /comment/icomment
	common /pepcrda/xpep(maxatm,8)
	common /pepcrdb/ypep(maxatm,8)
	common /vectors/iv(maxpar,4)
	common /par/natom,ntor,nhb,ns,lres
	common /inp/e_inp(maxstr),p1(maxstr,maxpar)
	common /pout/popmi(maxstr,maxpar),pemi(maxstr)
	common /hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
	common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,12)
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
	common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &  emo(maxstr),emi(maxstr)

	character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

	real fret,ftol,p(maxpar)
        real y(maxatm,8)
	integer iter, n

c	MODIFIED TO INCLUDE TRANS AND ROT PARAMETER
	n = ntor + 6 + np

        do i=kk11,kk11
	do j=1,n
		p1(i,j)=0.0
	enddo
	enddo

        if(icomment.eq.1)then
        do i=kk11,kk11
        do j=1,n
          print *,'pmin-opmo',i,j,opmo(i,j)
        enddo
        enddo
        endif

	call pgetx(kk11)

	do i=1,natom
          do j=4,8
             y(i,j)=ypep(i,j)
          enddo
	enddo

        do i = kk11,kk11
	n = ntor + 6 + np
	  do j=1,n
	   p(j) = p1(i,j)
	  enddo
	
        call pfrprmn(p,n,ftol,iter,fret)
        
	call poutp(p,fret,i)
        
c*******Flushing standard output************
	call flush(6)
c*******End of Flushing ********************
	enddo
 	return
	end
c***********************************************************************
        subroutine pfrprmn(p,n,ftol,iter,fret)
        integer iter, n, nmax, itmax
        real fret, ftol, p(n), eps, pfunc
        external pfunc
        parameter (nmax = 150,itmax=1000,eps=1.e-10)
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
        real dgg,fp,gam,gg,g(nmax),h(nmax),xi(nmax)
	ftol = 1.e-5

        fp = pfunc(p)

        call pdfunc(p,xi)
                
        do  j = 1,n
                g(j) = -xi(j)
                h(j) = g(j)
                xi(j) = h(j)
        enddo
c	loop over iterations
        do  its = 1, itmax
          iter = its
c	  print *,'iteration',its
	  call plinmin(p,xi,n,fret)

c	  next statement is the normal return
          if(2.*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+eps)) return
          fp = pfunc(p)
          call pdfunc(p,xi)

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
        return
        end
c********************************************************************
        subroutine plinmin(p,xi,n,fret)
	include 'mols.par'
	common /par/natom,ntor,nhb,ns,lres
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
        integer n,nmax
        real fret,p(n),xi(n),tol
        parameter(nmax =150,tol = 1.e-4)
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
        external pf1dim,pdf1dim

	n = ntor + 6 + np
        ncom = n

        do  j = 1,n
                pcom(j) = p(j)
                xicom(j) = xi(j)
        enddo
c	intial guess for brackets.
        ax = 0.
        xx = 1.
        call pmnbrak(ax,xx,bx,fa,fx,fb,pf1dim)

        fret = pdbrent(ax,xx,bx,pf1dim,pdf1dim,tol,xmin)

c	construct the vector results to return
        do  j = 1,n
                xi(j) = xmin*xi(j)
                p(j) = p(j) + xi(j)
        enddo

        return
        end
c**********************************************************************
        function pdbrent(ax,bx,cx,f,df,tol,xmin)
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
        pdbrent = fx
        return
        end
c********************************************************************
        subroutine pmnbrak(ax,bx,cx,fa,fb,fc,pfunc)
        real ax, bx, cx, fa, fb, fc, func, gold, glimit, tiny
        external pfunc
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
        fa = pfunc(ax)
        fb = pfunc(bx)
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
        fc = pfunc(cx)
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
                fu = pfunc(u)
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
                fu = pfunc(u)
c	parabolic fit is between c and its allowed limit
           else if((cx-u)*(u-ulim).gt.0.) then
                fu = pfunc(u)
                if(fu.lt.fc) then
                  bx = cx
                  cx = u
                  u = cx + gold*(cx-bx)
                  fb = fc
                  fc = fu
                  fu = pfunc(u)
                endif
c	limit parabolic u to maximum allowed value
          else if((u-ulim)*(ulim-cx).ge.0.) then
                u = ulim
                fu = pfunc(u)
c	reject parabolic u, use default magnification
          else
                u = cx + gold*(cx-bx)
                fu = pfunc(u)
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
	function pf1dim(x)
	integer nmax
	real pf1dim,pfunc,x
	parameter (nmax=150)
c	uses func
c	used by linimin as the function passed to mnbrak and dbrent
	external pfunc
	integer j, ncom
	real pcom(nmax),xicom(nmax),xt(nmax)
	common /f1com/pcom,xicom,ncom
	do j = 1,ncom
	   xt(j) = pcom(j)+x*xicom(j)
	enddo
	pf1dim = pfunc(xt)
	return
	end
c********************************************************************
        function pdf1dim(x)
        integer nmax
        real pdf1dim, x
        parameter (nmax =150)
c       uses dfunc
c	used by linimin as the function passed to mnbrak and dbrent
        integer j, ncom
        real df(nmax), pcom(nmax), xicom(nmax), xt(nmax)
        common /f1com/ pcom, xicom, ncom
        do  j = 1,ncom
                xt(j) = pcom(j)+x*xicom(j)
        enddo

        call pdfunc(xt,df)

        pdf1dim=0.
        do  j = 1,ncom
                pdf1dim = pdf1dim + df(j)*xicom(j)
        enddo
        return
        end
c********************************************************************
        subroutine pgetx(kk11)
	include 'mols.par'

        common /inp/e_inp(maxstr),p1(maxstr,maxpar)
        common /par/natom,ntor,nhb,ns,lres
	common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)

	n = ntor + 6 + np

	do i = kk11,kk11
	 e_inp(i) = emo(i)
	  do j = 1,n
	   p1(i,j) = opmo(i,j)
	  enddo
	enddo

	return
	end
c*********************************************************************
	function pfunc(p2)
	include 'mols.par'

	common /comment/icomment
	common /rctl/iscopt
	common /scang/frange,rang
	common /par/natom,ntor,nhb,ns,lres
	common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)

	dimension p2(maxpar)
	integer n
	real inter_ene,intra_ene

        if(icomment.eq.1)print *,'iff',iff
	n = 0
	n = ntor + 6 + np

	do j = 1,n-np
        p2(j) = mod(p2(j),360.0)
        if(p2(j).le.0) p2(j) = 360.0 + p2(j)
	phi(j) = p2(j)
	enddo

	do j = n-np-5,n-np-3
	tt1 = p2(j)
c---------------------------------------------------
	if(tt1.eq.0.0) then
	tt2 = 0.0
	else
	tt2 = (tt1/10.0)*0.13888889
	endif
	tt2 = tt2-2.5
	phi(j) = tt2
	enddo
c--------------------------------------------------
        do j = n-np+1,n
        if(iscopt.eq.1)then
        p2(j) = mod(p2(j),360.0)
        if(p2(j).le.0) p2(j) = 360.0 + p2(j)
        phi(j) = p2(j)
        elseif(iscopt.eq.0)then
        p2(j) = mod(p2(j),frange)
        if(p2(j).le.-(frange))p2(j) = p2(j)+frange
        endif
        phi(j) = p2(j)
        enddo
	
	call ppmolgen(2,2,phi,2,eflex)
c------------
		if(iff.eq.1) then
		 intra_ene = pampene(i,j)
		 call peplp(plpe,hb,steric,rep,2)
		 inter_ene = plpe
		endif

		if(iff.eq.2) then
		 intra_ene = pecpene(i,j)
		 call peplp(plpe,hb,steric,rep,2)
		 inter_ene = plpe
		endif

		pfunc = intra_ene + inter_ene + eflex
c------------
	return
	end
c*******************************************************************
	subroutine pdfunc(p2,g1)
	include 'mols.par'

	common /comment/icomment
	common /rctl/iscopt
	common /scang/frange,rang
	common /par/natom,ntor,nhb,ns,lres
	common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)

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
	enddo

	do j = n-np-5,n-np-3
	tt1 = p2(j)
c------------------------------------------------
	if(tt1.eq.0.0) then
	tt2 = 0.0
	else
	tt2 = (tt1/10.0)*0.13888889
	endif
	tt2 = tt2-2.5
	phi(j) = tt2
	enddo

	do j = n-np+1,n
        if(iscopt.eq.1)then
        p2(j) = mod(p2(j),360.0)
        if(p2(j).le.0) p2(j) = 360.0 + p2(j)
        phi(j) = p2(j)
        elseif(iscopt.eq.0)then
        p2(j) = mod(p2(j),frange)
        if(p2(j).le.-(frange))p2(j) = p2(j)+frange
        phi(j) = p2(j)
        endif
        enddo


c	MODIFIED P2->PHI
	call ppmolgen(2,2,phi,2,eflex)
c------------
		if(iff.eq.1) then
		 intra_ene = pampene(i,j)
		 call peplp(plpe,hb,steric,rep,2)
		 inter_ene = plpe
		endif

		if(iff.eq.2) then
		 intra_ene = pecpene(i,j)
		 call peplp(plpe,hb,steric,rep,2)
		 inter_ene = plpe
		endif
c------------
	f1 = intra_ene + inter_ene + eflex

	do j=1,n-np-6
	plp = 0.0
	g1(j) = 0.0
	deltaE = 0.0
	ej = 0.0
	p2(j) = p2(j) - deltaX
	phi(j) = p2(j)
	call ppmolgen(2,2,phi,2,eflex)
c------------
		if(iff.eq.1) then
		 intra_ene = pampene(i,j)
		 call peplp(plpe,hb,steric,rep,2)
		 inter_ene = plpe
		endif

		if(iff.eq.2) then
		 intra_ene = pecpene(i,j)
		 call peplp(plpe,hb,steric,rep,2)
		 inter_ene = plpe
		endif
c------------
	ej = intra_ene + inter_ene + eflex

	deltaE = f1-ej
	g1(j) = deltaE/deltaX
	p2(j) = p2(j) + deltaX
	phi(j) = p2(j)
	enddo

	do j=n-np-5,n-np-3
	plp = 0.0
	g1(j) = 0.0
	deltaE = 0.0
	ej = 0.0
	temp = phi(j)
c--------------------------------------------
	p2(j) = p2(j) - deltaX
c--------------------------------------------
	tt1 = p2(j)
c--------------------------------------------
	if(tt1.eq.0.0) then
	tt2 = 0.0
	else
	tt2 = (tt1/10.0)*0.13888889
	endif
	tt2 = tt2-2.5
	phi(j) = tt2

	call ppmolgen(2,2,phi,2,eflex)
c------------
		if(iff.eq.1) then
		 intra_ene = pampene(i,j)
		 call peplp(plpe,hb,steric,rep,2)
		 inter_ene = plpe
		endif

		if(iff.eq.2) then
		 intra_ene = pecpene(i,j)
		 call peplp(plpe,hb,steric,rep,2)
		 inter_ene = plpe
		endif
c------------
	ej = intra_ene + inter_ene + eflex
	deltaE = f1-ej
	g1(j) = deltaE/deltaX
	p2(j) = p2(j) + deltaX
	phi(j) = temp
	enddo

	do j=n-np-2,n-np
	plp = 0.0
	g1(j) = 0.0
	deltaE = 0.0
	ej = 0.0
	p2(j) = p2(j) - deltaX
	phi(j) = p2(j)
	call ppmolgen(2,2,phi,2,eflex)
c------------
		if(iff.eq.1) then
		 intra_ene = pampene(i,j)
		 call peplp(plpe,hb,steric,rep,2)
		 inter_ene = plpe
		endif

		if(iff.eq.2) then
		 intra_ene = pecpene(i,j)
		 call peplp(plpe,hb,steric,rep,2)
		 inter_ene = plpe
		endif
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
	call ppmolgen(2,2,phi,2,eflex)
c------------
                if(iff.eq.1) then
                 intra_ene = pampene(i,j)
                 call peplp(plpe,hb,steric,rep,2)
                 inter_ene = plpe
                endif

                if(iff.eq.2) then
                 intra_ene = pecpene(i,j)
                 call peplp(plpe,hb,steric,rep,2)
                 inter_ene = plpe
                endif
c------------
        ej = intra_ene + inter_ene + eflex
        deltaE = f1-ej
        g1(j) = deltaE/deltaX
        p2(j) = p2(j) + deltaX
        phi(j) = p2(j)
        enddo

	if(icomment.eq.1)print *,'g1: ',(g1(j),j=1,n)
	return
	end
c*******************************************************************
	subroutine poutp(p3,ene,i)
	include 'mols.par'

	common /comment/icomment
	common /par/natom,ntor,nhb,ns,lres,ie
	common /pout/popmi(maxstr,maxpar),pemi(maxstr)
	common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

	character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

	dimension p3(maxpar)
	integer n

	n = ntor + 6 + np

	do k1=1,n
        opmi(i,k1) = p3(k1)
	enddo
	if(icomment.eq.1) print *,'opmi-pmin',i,(opmi(i,k1),k1=1,n)

        emi(i) =ene 
	if(icomment.eq.1) print *,'pep-plp-pro-min-ene',emi(i)

	write(31,FMT='(5x,a23,2x,i4,2x,f10.3)')' Minimized 
     &Structure:',i,emi(i)
	write(*,*)'Minimized Structure No :',i,emi(i)

!-------'pemi' added to avoid confusion with receptor 'emi'--------
	open(unit=14,file=of2,status='unknown')
        write(14,*)i,ene,(p3(j),j=1,n)
	if(i.eq.ie) close(unit=14)
        return
        end
c*******************************************************************
        subroutine ppmolgen(lx,ly,phi,tst,eflex)

        include 'mols.par'

	common /pepcrda/xpep(maxatm,8)
        common /pepcrdb/ypep(maxatm,8)
	common /vectors/iv(maxpar,4)
	common /cen/bx,by,bz
	common /par/natom,ntor,nhb,ns,lres
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)

        dimension x_one(maxatm,3),x_two(maxatm,3)
        dimension phi(maxpar)
        integer tst,nn
        real rotang,theta,psi
        real cx,cy,cz,bx,by,bz

        nn = ntor
        do k=1,natom
           do ki=1,3
             x_one(k,ki)=xpep(k,ki)
           enddo
        enddo

        do if=1,nn
C###################### PHI ALL####################################
        call elemen(x_one(iv(if,1),1),x_one(iv(if,1),2),
     &              x_one(iv(if,1),3),
     &              x_one(iv(if,2),1),x_one(iv(if,2),2),
     &              x_one(iv(if,2),3),
     &              el,em,en)

        do k=1,iv(if,3)-1
           do ki=1,3
             x_two(k,ki)=x_one(k,ki)
           enddo
        enddo

        do k=iv(if,3),iv(if,4)
           xin=x_one(k,1)-x_one(iv(if,1),1)
           yin=x_one(k,2)-x_one(iv(if,1),2)
           zin=x_one(k,3)-x_one(iv(if,1),3)
           call rotor(el,em,en,phi(if),xin,yin,zin,
     &                xout,yout,zout)
           x_two(k,1)=xout+x_one(iv(if,1),1)
           x_two(k,2)=yout+x_one(iv(if,1),2)
           x_two(k,3)=zout+x_one(iv(if,1),3)
        enddo

        do k=iv(if,4)+1,natom
           do ki=1,3
             x_two(k,ki)=x_one(k,ki)
           enddo
        enddo

        do k=1,natom
           do ki=1,3
              x_one(k,ki)=x_two(k,ki)
           enddo
        enddo
C####################################################################

        enddo

        do k=1,natom
         do ki=1,3
          ypep(k,ki)=x_two(k,ki)
         enddo
        enddo

c***********translation and rotation for ligand postioning***********

        call protate(bx,by,bz,phi(nn+4),phi(nn+5),phi(nn+6),natom)
        call ptranslate(bx,by,bz,phi(nn+1),phi(nn+2),phi(nn+3),natom)

	if(ifopt.eq.2)then
        call flexall(2,2,phi,eflex,1)
        endif

c*********************************************************************
        return
        end
c---------------------------------------------------------------------
