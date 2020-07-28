c     library name : mols.f   

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



c       program mols
	subroutine mols(fseq,iopt,iff)

	include 'mols.par'
	parameter(mxtyat = 18)
	common /ecpp/ aij(mxtyat,mxtyat),cij(mxtyat,mxtyat),
     &  a14(mxtyat,mxtyat),ihb(mxtyat,mxtyat),ahb(mxtyat,mxtyat),
     &  chb(mxtyat,mxtyat)
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/crda/x(maxatm,8) 
        common/ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,12)
        common/crdb/y(maxatm,8) 
        common/energy/e_final(maxord,maxord) 
        common/mean/avrg1(maxpar,maxord)
        common/vectors/iv(maxpar,4)
        common /par/ natom, ntor, nhb, ns, lres
        common/hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
        common/order/nn,mm
        common/freq/ifreq(maxpar,maxord,maxord)
        common/out/e_out(maxatm),p_out(maxpar,maxatm)
	common/mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
        common/calls/ncalls
	common /scc/ sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
	common/pdbat/atom(maxatm),ele(maxatm)
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4,of0

         dimension angl(maxpar,maxord),angle(maxpar,maxord),jseed(5000),
     &ang(maxpar,maxord)
	character seq*(maxres),fseq(maxatm)*1

c        nn=19 ! No. of torsion angles ( or search dimensions)
        mm=37 ! Size of MOLS dimensions
c	mm = 13 !sam :: trial 
c	natom=77
cc	write(31,*) 'Enter the no. of parameters'
cc	read *, nn
	nn = ntor
	npar = ntor
cc	write(31,*) 'Enter the no. of atoms'
cc	read *,natom
cc	write(31,*) 'How many optimal structures to be generated (1 to 2000)'
cc	read *,ns
cc	write(31,*) 'No. of Hydrogen bond pairs'
cc	read *,nhb

ccc	call input
ccc	call pdb
	call anggen(ang,fseq,npar,iopt,aii) !here angles are generated
	if(iopt.eq.3) call scinp(fseq) !iopt = 3 = rotamer,fseq = sequence
        do i=1,natom
          do j=1,8
             y(i,j)=x(i,j)
          enddo          
	enddo
c	et=ecpene(1,1)
ccc	print *,'etot',et
ccc	stop
cccccccccc
ccc	return
ccccccccccc
c        open(unit=2,file='n37.inp',status='old')
c	print *,npar
        do i=1,npar
	  do j = 1,mm
          angl(i,j) =ang(i,j)
	  enddo
        enddo

c275      format(37f6.1)
c         close(unit=2)
c************************************************************
c        open(unit=3,file='seed.inp',status='unknown')
	call rand1(jseed)
         nt=1
         do ijx=1,ns
	iseed = jseed(ijx)
c          read(3,*)iseed
ccc	write(31,*)iseed
          call write_par(angl,iseed,angle,npar)
          call pargen1(angle,npar)
        do i=1,100
          e_out(i)=1.0e+25
        enddo

        do i=1,npar
          do j=1,mm
            do k=1,mm
              ifreq(i,j,k)=0
            enddo
          enddo
        enddo
c         ************ Main MOLS ***************
         do l1=1,npar
         do l2=1,mm
            p(l1,l2,2)=0.0
            p(l1,l2,3)=0.0
          enddo
         enddo

          do i=1,mm
           do j=1,mm
	     k = 1
	
	     call subpar(i,j,k,iopt,fseq,aii,phi)!output from subpar is 'phi'
	     call molgen(phi)!'phi' is taken as the input from subroutine 'subpar'
cc
		if(iff.eq.1) e_final(i,j)=ampene(i,j)
		if(iff.eq.2) e_final(i,j)=ecpene(i,j)
c		print *,e_final(i,j),i,j
cc
c	     call energee(i,j)
           enddo
          enddo
          call average(npar)
          do i=1,npar
            call sort_and_set_rank(i)
          enddo
         do i=1,npar 
           call rank_sort(i,nt)
         enddo
c         ************ LAST MOLS *************
	do m1=1, 15
          m2 = 2
	  kk = 1
	  call subpar(kk,m1,m2,iopt,fseq,aii,phi)
	  call molgen(phi)
cc
	if(iff.eq.1)   e_final(m1,m1) = ampene(m1,m1)
	if(iff.eq.2)  e_final(m1,m1) = ecpene(m1,m1)
c	print *,m1,(phi(ih),ih=1,nn),e_final(m1,m1)
cc
c	call energee(m1,m1)
	  call best(m1,npar)
	enddo
        call output(1,ijx,et,iopt,fseq,aii)
	write(*,*)'generated MOLS optimal Structure NO : ',ijx,et
c*******Flushing standard output************
	call flush(6)
c*******End of Flushing ********************
	write(31,*)'generated MOLS optimal Structure NO : ',ijx,et
333	enddo
c	stop
	return
        end

c********************************************************************
        subroutine write_par(angl,iseed,angle,npar)

	include 'mols.par'
        common/order/nn,mm
        dimension angll(maxpar,maxord),angle(maxpar,maxord),
     &  angl(maxpar,maxord)
       open(unit=4,file='molspargen.inp',status='unknown')

        do i=1,npar
          do j=1,mm
            angll(i,j)=angl(i,j)
          enddo
        enddo
        amm=float(mm)
        do j=1,npar
         ki=0
         do i=1,100000
c          y=ran(iseed)
          y=rand(iseed)
          iseed = y
          iy=nint(y*amm)
          if(iy.eq.0)iy=1
           if (angll(j,iy).lt.999) then
            ki=ki+1
            angle(j,ki)=angll(j,iy)
            angll(j,iy)=99999.0
            if(ki.ge.mm)go to 3 
           endif         
          enddo         
3       continue
        enddo 
     
        do i=1,npar
         write(4,10)(angle(i,j),j=1,mm)
10       format(37f6.1)
        enddo
    
        close(unit=4)
        return 
        end
c*********************************************************************
	subroutine pargen1(angle,npar)

c	pargen is for parameter generation
c	n  -   total number of parameters
c	m  -   maximum number of values for each parameter
      
	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/order/nn,mm
	dimension angle(maxpar,maxord)
	open (unit=5,file='molspargen.inp',status='unknown')

	do 150 i = 1, npar
	  do 150 j = 1,mm
	  p(i,j,1) = angle(i,j)
 150	continue
  10    format(37f6.1)
	do 500 l = 1, npar
	do 500 k = 1, mm
	do 500 j = 1, mm
	  i = mod( ((j-1) * (l-1) + (k-1)), mm) + 1
	  e(i,j,l)  = p(l,k,1)
 500	continue
        close (unit=5)
        
20	format(8f6.1)
	open (unit=6,file='pargen.out',status='unknown')
	 write (6,20) (((e(i,j,l),l=1,npar),j=1,mm),i=1,mm)
       close (unit=6)        

	return
	end

c************************************************************************
	subroutine molgen(phi)

	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/crda/x(maxatm,8) 
        common/crdb/y(maxatm,8) 
        common /par/ natom, ntor, nhb, ns, lres
        common/vectors/iv(maxpar,4)
        common/order/nn,mm
 
        dimension x_one(maxatm,3),x_two(maxatm,3)
        dimension phi(maxpar)

c	write(31,*)(phi(ii),ii=1,nn)
c	write(*,*)(phi(ii),ii=1,nn)

        do k=1,natom
           do ki=1,3
             x_one(k,ki)=x(k,ki) 
           enddo
        enddo
        do if=1,nn
 
C###################### PHI ALL ####################################      

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
c	   write(31,*)k,phi(if),xin,yin,zin,xout,yout,zout
c          write(*,*)'rotor-done'
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
             y(k,ki)=x_two(k,ki) 
           enddo
        enddo

c        do k=1,atom
c         write(*,*),(y(k,ki),ki=1,3)
c        enddo

        return
        end
c***********************************************************************
        function ampene(ie,je)
 
	include 'mols.par'
        common/crdb/y(maxatm,8) 
        common/ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,12)
        common/energy/e_final(maxord,maxord) 
        common /par/ natom, ntor, nhb, ns, lres
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common/hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
        common/calls/ncalls
c	open (unit=7,file='energy.out',status='unknown')
	
        ees=0.0
        enb=0.0
        enb1_4=0.0
        ees1_4=0.0
        ehb=0.0
        etot=0.0
	

	do i=1,natom
          n_range=ifix(y(i,7))
         do ij=1,n_range
          k1=jstart(i,ij)
         k2=jend(i,ij)
          do j=k1,k2
          dis = dist(y(i,1),y(i,2),y(i,3),y(j,1),y(j,2),y(j,3))
            if(dis.lt.0.01) then
!              e_final(ie,je)=1.0e+25
	       ampene = 1.0e+25
              return
            endif
            call force(i,j,dis,enbt,eest)
            ees=ees+eest
            enb=enb+enbt
          enddo
         enddo

         n1_4=ifix(y(i,8))
         do ij=1,n1_4
            j=j1_4(i,ij)
          dis = dist(y(i,1),y(i,2),y(i,3),y(j,1),y(j,2),y(j,3))
            if(dis.lt.0.5) then
              STOP 'Input coordinates are wrong !!'
            endif
            call force(i,j,dis,enbt,eest)
            ees1_4=ees1_4+0.5*eest
            enb1_4=enb1_4+0.5*enbt
         enddo           
       enddo

        do i=1,nhb
         dis = dist( y(ihb1(i),1),y(ihb1(i),2),y(ihb1(i),3),
     &                  y(ihb2(i),1),y(ihb2(i),2),y(ihb2(i),3))
c     	 print *,y(ihb1(i),1),y(ihb1(i),2),y(ihb1(i),3),
c    &                  y(ihb2(i),1),y(ihb2(i),2),y(ihb2(i),3)	
         ess=dis*dis
         artwo=1.0/ess
         arten=artwo**5
         artwelve=arten*artwo
         ehb=ehb+c(i)*artwelve-d(i)*arten
         arstar=y(ihb1(i),5)+y(ihb2(i),5)
         epsilon=sqrt(y(ihb1(i),6)*y(ihb2(i),6))
         aaa=epsilon*(arstar**12)
         ccc=2.0*epsilon*(arstar**6)
         arsix=artwo*artwo*artwo
         ef1=aaa*artwelve
         ef2=ccc*arsix
         enbt=ef1-ef2
         enb=enb-enbt
        enddo
 


        ampene=enb+ees+ehb+enb1_4+ees1_4
c        etot=enb+ees+ehb+enb1_4+ees1_4
c         e_final(ie,je)=etot
c	write(7,*) e_final(ie,je)
c	write(31,*)e_final(ie,je)
c	close(unit=7)
c	write(*,*) e_final(ie,je),ie,je
        return
        end
c*******************************************************************
        subroutine force(if,jf,diss,enbt,eest)

	include 'mols.par'
        common/crdb/y(maxatm,8) 

        arstar=y(if,5)+y(jf,5)
        epsilon=sqrt(y(if,6)*y(jf,6))
        aaa=epsilon*(arstar**12)
        ccc=2.0*epsilon*(arstar**6)
        ess=diss*diss
        artwo=1.0/ess
        arsix=artwo*artwo*artwo
        artwelve=arsix*arsix
        eest=332.0*y(if,4)*y(jf,4)*artwo/4.0
        ef1=aaa*artwelve
        ef2=ccc*arsix
        enbt=ef1-ef2        
        return                   	
        end

c***********************************************************************
        function ecpene(ie,je)
 
	include 'mols.par'
	parameter(mxtyat = 18)
        common/crdb/y(maxatm,8) 
        common/ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,12)
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
         common/energy/e_final(maxord,maxord) 
        common /par/ natom, ntor, nhb, ns, lres
        common/hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
        common /ecpp/ aij(mxtyat,mxtyat),cij(mxtyat,mxtyat),
     &  a14(mxtyat,mxtyat),ihb(mxtyat,mxtyat),ahb(mxtyat,mxtyat),
     &  chb(mxtyat,mxtyat)
        common/calls/ncalls
c	open (unit=7,file='energy.out',status='unknown')
	
	cdc=(22.0/(7.0*180))      !pi/180
        ees=0.0
        enb=0.0
        ehb=0.0
	etor=0.0
        enb14=0.0
        ees14=0.0
	ehb14=0.0
        etot=0.0
c	do i = 1,natom
c	print *,i,y(i,1),y(i,2),y(i,3)
c	enddo
c
 	do in = 1,ntor
	etor = etor + u0(in)*(1.0+sn(in)*cos(cdc*((phi(in)-180.0))*tn(in)))
c	print *,in, etor , tn(in),u0(in),sn(in),phi(in)
	enddo
c
	do i=1,natom
          n_range=ifix(y(i,7))
         do ij=1,n_range
          k1=jstart(i,ij)
         k2=jend(i,ij)
          do j=k1,k2
          dis = dist(y(i,1),y(i,2),y(i,3),y(j,1),y(j,2),y(j,3))
            if(dis.lt.0.01) then
!              e_final(ie,je)=1.0e+25
	       ecpene = 1.0e+25
              return
            endif
c
	ity = int(y(i,5))
	jty = int(y(j,5))
        ees = ees + (332.0*y(i,4)*y(j,4))/(dis*2.0)
	if(ihb(ity,jty).eq.0) then
	enb = enb + (aij(ity,jty)/(dis**12.0))-(cij(ity,jty)/(dis**6.0))
c       print *,i,j,dis,aij(ity,jty),cij(ity,jty),enb
	else
         ehb=ehb+(ahb(ity,jty)/(dis**12.0))-(chb(ity,jty)/(dis**10.0))
	endif
          enddo        
         enddo
c
         n1_4=ifix(y(i,8))
         do ij=1,n1_4
            j=j1_4(i,ij)
          dis = dist(y(i,1),y(i,2),y(i,3),y(j,1),y(j,2),y(j,3))
            if(dis.lt.0.5) then
              STOP 'Input coordinates are wrong !!'
            endif
c
	ity = int(y(i,5))
	jty = int(y(j,5))
        ees14 = ees14 + (332.0*y(i,4)*y(j,4))/(dis*2.0)
	if(ihb(ity,jty).eq.0) then
	enb14=enb14+(a14(ity,jty)/(dis**12.0))-(cij(ity,jty)/(dis**6.0))
	else
       ehb14=ehb14+(ahb(ity,jty)/(dis**12.0))-(chb(ity,jty)/(dis**10.0))
	endif
         enddo           
       enddo

c
	ees = ees + ees14
	enb = enb + enb14
	ehb = ehb + ehb14
        ecpene=enb+ees+ehb+etor
c	print *, ecpene,ie,je
c	ecpene = etot
c	print *,'elene1_4 =',ees14
c	print *,'vwene1_4 =',enb14
c	print *,'hbene1_4 =',ehb14
c	print *,'elene =',ees
c	print *,'vwene =',enb
c	print *,'hbene =',ehb
c	print *,ees,enb,ehb,etor
c	print *,'etot = ',etot
c         e_final(ie,je)=etot
c	close(unit=7)
        return
        end
c***********************************************************************
c********************************************************************
        subroutine average(npar)

	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/energy/e_final(maxord,maxord) 
        common/out/e_out(maxatm),p_out(maxpar,maxatm)
        common/mean/avrg1(maxpar,maxord)
        common/order/nn,mm

        data en_k_t/1.98578e+03/
c	open (unit=8,file='average.out',status='unknown')

c       en_k_t is boltzmannconst*avogadronumber*temp(=300 degree K)
c	average is for finding the average 
c	e_final   -  energy values of the conformation 
c	cutoff    -  for omitting impossible conformations 
c	avrg1     -  average without cutoff
c	n  -   total number of parameters = 9
c	m  -   maximum number of values for each parameter = 37

c	compute using mols

	do 400 l = 1, npar
	  do 500 k = 1, mm
	    avrg1(l,k) = 0.0 
            summ   = 0.0
            summ_mm= 0.0
	    do 600 j = 1, mm
	       i = mod( ( (j-1) * (l-1) + (k-1) ), mm) + 1
               weight=exp(-1.0*(e_final(i,j)/(en_k_t*0.01)))
  	       summ = summ + e_final(i,j)*weight
               summ_mm=summ_mm+weight
600	    continue
	    avrg1(l,k) = summ/summ_mm
c	write(8,*)avrg1(l,k)
c	close(unit=8) 
500	  continue
400    continue

        return
        end
c*********************************************************************
      
        subroutine sort_and_set_rank(is1)

	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/mean/avrg1(maxpar,maxord)
        common/order/nn,mm
        common/freq/ifreq(maxpar,maxord,maxord)

        dimension rank(maxord),icc(maxord)
c	open (unit=9,file='srank1.out',status='unknown')
c	open (unit=10,file='srank2.out',status='unknown')
c	open (unit=11,file='srank3.out',status='unknown')
c	open (unit=12,file='srank4.out',status='unknown')

c	do j=1,mm
c	write (9,*),avrg1(is1,j),p(is1,j,1)
c       enddo
c	close(unit=9)
        do j=2,mm
          aa=avrg1(is1,j)
          bb=p(is1,j,1)
          do i=j-1,1,-1
            if(avrg1(is1,i).le.aa) go to 10
            avrg1(is1,i+1)=avrg1(is1,i)
            p(is1,i+1,1)=p(is1,i,1)
          enddo
          i=0
10        avrg1(is1,i+1)=aa
          p(is1,i+1,1)=bb
        enddo
c	do j=1,mm
c	write (10,*),avrg1(is1,j),p(is1,j,1)
c        enddo
c	close(unit=10)

        do i=1,mm
          rank(i)=float(i)
        enddo

        do j=2,mm
          aa=rank(j)
          bb=p(is1,j,1)
          do i=j-1,1,-1
            if(p(is1,i,1).le.bb) go to 20
            p(is1,i+1,1)=p(is1,i,1)
            rank(i+1)=rank(i)
          enddo
          i=0
20        rank(i+1)=aa
          p(is1,i+1,1)=bb
        enddo

c	do j=1,mm
c	write (11,*)rank(j),p(is1,j,1)
c       enddo
c	close(unit=11)

        do i=1,mm
          jxx=ifix(rank(i))
          ifreq(is1,i,jxx)=ifreq(is1,i,jxx)+1
          p(is1,i,2)=p(is1,i,2)+rank(i)
          p(is1,i,3)=p(is1,i,3)+rank(i)*rank(i)
        enddo


c	do j=1,mm
c	write (12,*) p(is1,j,1),p(is1,j,2),p(is1,j,3)
c        enddo
c	close(unit=12)

        return
        end 
c********************************************************************
       
        subroutine rank_sort(ir1,nt)

	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/freq/ifreq(maxpar,maxord,maxord)
        common/order/nn,mm

        dimension ee(maxord)
c	open (unit=13,file='rank1.out',status='unknown')
c	open (unit=14,file='rank2.out',status='unknown')
c	open (unit=15,file='rank3.out',status='unknown')
c	open (unit=16,file='rank4.out',status='unknown')

c       calculate the average and sd
c	do j=1,mm
c	write (13,*) p(ir1,j,1),p(ir1,j,2),p(ir1,j,3)
c	enddo
c	close(unit =13)

        do i=1,mm
          sum_x=p(ir1,i,2)
          sum_x2=p(ir1,i,3)
          ave=sum_x/nt
          sd=sqrt((sum_x2/nt)-((sum_x/nt)*(sum_x/nt)))
          p(ir1,i,2)=ave
          p(ir1,i,3)=sd
        enddo       
c	do j=1,mm
c	write (14,*) p(ir1,j,1),p(ir1,j,2),p(ir1,j,3)
c	enddo
c	close(unit =14)

        do j=2,mm
          bb=p(ir1,j,1)
          cc=p(ir1,j,2)
          dd=p(ir1,j,3)
          do j1=1,mm
            ee(j1)=ifreq(ir1,j,j1)
          enddo
          do i=j-1,1,-1
            if(p(ir1,i,2).le.cc) go to 10
            p(ir1,i+1,3)=p(ir1,i,3)
            p(ir1,i+1,2)=p(ir1,i,2)
            p(ir1,i+1,1)=p(ir1,i,1)
            do j1=1,mm
               ifreq(ir1,i+1,j1)=ifreq(ir1,i,j1)
            enddo
          enddo
          i=0
10        p(ir1,i+1,1)=bb
          p(ir1,i+1,2)=cc
          p(ir1,i+1,3)=dd
          do j1=1,mm
            ifreq(ir1,i+1,j1)=ee(j1)
          enddo
        enddo

c	do j=1,mm
c	write (15,*) p(ir1,j,1),p(ir1,j,2),p(ir1,j,3)
c	enddo
c	close(unit =15)

	do j=2,mm
	 bb=p(ir1,j,1)
	 cc=p(ir1,j,2)
	 dd=p(ir1,j,3)
	 do j1=1,mm
         ee(j1)=ifreq(ir1,j,j1)
	 enddo
	  do i=j-1,1,-1
	   if (p(ir1,i,2).ne.cc) goto 40
	   if (p(ir1,i,3).le.dd) goto 40
	   p(ir1,i+1,3)=p(ir1,i,3)
	   p(ir1,i+1,2)=p(ir1,i,2)
	   p(ir1,i+1,1)=p(ir1,i,1)
	   do j1=1,mm
	   ifreq(ir1,i+1,j1)=ifreq(ir1,i,j1)
	   enddo
	  enddo
	  i=0
  40 	  p(ir1,i+1,1)=bb
	  p(ir1,i+1,2)=cc
	  p(ir1,i+1,3)=dd
	  do j1=1,mm
	  ifreq(ir1,i+1,j1)=ee(j1)
	  enddo
	enddo

c	do j=1,mm
c	write (16,*) p(ir1,j,1),p(ir1,j,2),p(ir1,j,3)
c	enddo
c	close(unit =16)

        return
        end

c*********************************************************************
       subroutine output(nou,ijx,et,iopt,fseq,aii)

	include 'mols.par'
        common/out/e_out(maxatm),p_out(maxpar,maxatm)
        common/order/nn,mm
        common /par/ natom, ntor, nhb, ns, lres
	common /scc/ sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
	common/mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4,of0
	character fseq(maxatm)*1 
	dimension pp(maxpar)

	i = nou
	if(iopt.eq.3) then
	j = 1
	k = 1
	jh = 1
172	if(fseq(k).eq.'P') go to 193
	pp(j) = p_out(jh,i)
	j = j + 1
	jh = jh + 1
193	pp(j) = p_out(jh,i)
	j = j + 1
	jh = jh + 1
	if(nrot(k).ne.0.) then
	j1 = p_out(jh,i)
c	j1 = int((j1/10.0) + 1.0)
	j1 = int((j1/aii) + 1.0)
c	print *,j1,aii
	j1 = mod(j1,nrot(k))
	if(j1.eq.0) j1 = nrot(k)
	do j2 = 1, nchi(k)
		pp(j) = sc(k,j1,j2)
		j = j + 1
	enddo
	jh = jh + 1
	endif
	k = k + 1
	if(k.le.lres) go to 172	
	endif

	open(unit=24,file=of1,status='unknown')
c          write(24,10)(p_out(j,i),j=1,npar),e_out(i)
c	   print *,'within output: ',(pp(j),j=1,nn)
c           write(24,*)(p_out(j,i),j=1,npar),e_out(i) 
	if(iopt.eq.3) then
	do k1=1,nn
	opmo(ijx,k1) = pp(k1)
	enddo 
	emo(ijx) = e_out(i)
	do k1 = 1,nn
	pp(k1) = pp(k1) - 180.0
	enddo
c           write(24,*)(pp(j),j=1,nn),e_out(i)
C           write(24,10) ijx,e_out(i),(pp(j),j=1,nn)
           write(24,*) ijx,e_out(i),(pp(j),j=1,nn)
c	   print *,'from output',ijx,e_out(i),(pp(j),j=1,nn)
c	write(*,*)(pp(j),j=1,nn),e_out(i)
	else
	do k1=1,nn
	opmo(ijx,k1) = p_out(k1,i)
	enddo
	emo(ijx) = e_out(i)
	do k1 = 1,nn
	p_out(k1,i) = p_out(k1,i) - 180.0
	enddo
c           write(24,*)(p_out(j,i),j=1,nn),e_out(i) 
C           write(24,10)ijx,e_out(i),(p_out(j,i),j=1,nn)
c	    print *,'from output',ijx,e_out(i),(p_out(j,i),j=1,nn)
           write(24,*)ijx,e_out(i),(p_out(j,i),j=1,nn)
	endif	
	et = e_out(i)
C10	format(i4,1x,f15.2,<nn>(1x,f6.1))
	if(ijx.eq.ns) close(unit=24)
        return
        end

c********************************************************************
c*********************************************************************
c***********************************************************************
	subroutine best(ib,npar)
	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/energy/e_final(maxord,maxord) 
        common/order/nn,mm
        common/out/e_out(maxatm),p_out(maxpar,maxatm)

        do i=1,25
          if(e_final(ib,ib).lt.e_out(i)) then            
            do j=25,i+1,-1
              e_out(j)=e_out(j-1)
              do k=1,npar
                p_out(k,j)=p_out(k,j-1)
              enddo
            enddo
            e_out(i)=e_final(ib,ib)
	    do l = 1,npar
	      p_out(l,i) = p(l,ib,1)
	    enddo
            go to 10
          endif
         enddo
10       continue
         return
         end
C###############################################################################
	subroutine rand1(jseed)
	integer a,b,c,d,e,r,jseed(5000)
c       open(unit=8,file='mseed.out',status='unknown')
	a=2346
	b=3465
	c=5421
	d=5323
	e=5000 
	  nss=0
           do i=1,10000
             r = mod(((a*b)+c),d)
	if (r.lt.0)then
	r = r * -1
	endif
             r=(r*1)*2+1
             if(r.le.999.or.r.gt.9999) go to 1
             nss=nss+1
             if(nss.gt.e)go to 2
c            write(8,*)r
	jseed(nss) = r	
1       b = r
	enddo

2       return
        end
c*************************************************************************
	subroutine anggen(ang,fseq,npar,iopt,aii)
	include 'mols.par'
        common/order/nn,mm
        common /par/ natom, ntor, nhb, ns, lres
	dimension ang(maxpar,maxord)
	character fseq(maxatm)*1 
	
c	open (unit=21,file='angles.inp',status='unknown')
c	type *,'enter nn (total no. of parameters)'
c	read *,nn
c	type *,'enter mm ((37)maximum no. of values for each parameters)'
c	read *,mm

  10    format(37f6.1)
cc        type *,'enter the interval to generate input angel values'
cc        read *,ii
cccc	ii = 10
cc        type *,'enter the starting value to generate input angl values'
cc	read *,ll

	if(iopt.eq.3) then
		npar = lres*3
		do i = 1,lres
		if(fseq(i).eq.'A'.or.fseq(i).eq.'G'.or.fseq(i).eq.'P') then
		npar = npar - 1
		endif
		enddo
	endif
cc-------------------------------------------------------------------
cc	check the mols order it should be greater than npar and prime
cc	number
	print *,mm,npar
	if(mm.lt.npar) then
	do i1 = 1,100
	k1 = 2
181	if(i1/k1*k1.eq.i1) go to 191
	k1 = k1 + 1
	if(i1.le.i1/2) go to 181
	if(i1.gt.npar) go to 201
c	print *,i1
191	enddo
201	mm = i1
	i1 = 0
	endif
c	print *,mm,npar
cc-------------------------------------------------------------------
	aii = 360.0/(mm-1)
	print *,'o.k',aii
	all = 0
	do i=1,npar
	ang(i,1)=all
        b=all
	do j=2,mm
c        write(21,10)ang(i,j)
	ang(i,j)=b+aii
        b=ang(i,j)
        enddo
	enddo
c        write(21,10)((ang(i,j),j=1,mm),i=1,npar)
c        write(6,10)((ang(i,j),j=1,mm),i=1,npar)
c	close(unit=21)
	return 
	end
c*************************************************************************
	subroutine scinp(fseq)
	include 'mols.par'
        common /par/ natom, ntor, nhb, ns, lres
	character fseq(maxatm)*1,ch*1, tcode*4, scode*1
	common /scc/ sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
cc	dimension sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
20	format(a4,a1,1x,i2,i2)
21	format(20x,5f7.1)
11	format(a18,i3)

c	write(6,*) 'The sequence   : ',(fseq(i),i=1,lres)

 	open(unit=33,file='SC.lib',status='old')
	do i = 1,lres
	rewind 33
	do i1 = 1,234
	read(33,20) tcode, scode, mrot, mchi
c	write(31,*)tcode, scode, mrot, mchi
	if(scode.eq.fseq(i).and.mrot.ne.0) then
	  do i2 = 1, mrot
	  read(33,21) (sc(i,i2,i3),i3=1,mchi)
c	write(31,*)(sc(i,i2,i3),i3=1,mchi)
	  enddo
	nrot(i) = mrot
	nchi(i) = mchi
	go to 212
	endif
	enddo
212	enddo
	do i = 1,lres
	do j = 1, nrot(i)
c	write(31,*)(sc(i,j,k),k=1,nchi(i))
	enddo
	enddo
	close(unit=33)
	return
	end
	
c************************************************************************	
	subroutine subpar(l1,l2,l3,iopt,fseq,aii,phi)
	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/crda/x(maxatm,8) 
        common/crdb/y(maxatm,8) 
        common /par/ natom, ntor, nhb, ns, lres
        common/vectors/iv(maxpar,4)
        common/order/nn,mm
	common /scc/ sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
 
	character fseq(maxatm)*1 
c	dimension sc(3000,50,5),nrot(50),nchi(5)
        dimension x_one(maxatm,3),x_two(maxatm,3)
        dimension phi(maxpar)
	if(iopt.eq.3) then !iopt = 3 = rotamer
	i = 1
	j = 1
	k = 1
161	if(l3.eq.1) then
	if(fseq(k).eq.'P') go to 191
	phi(i) = e(l1,l2,j)!e()=stored angles
	i = i + 1
	j = j + 1
191	phi(i) = e(l1,l2,j)
	i = i + 1
	j = j + 1
	if(nrot(k).ne.0.) then
	j1 = e(l1,l2,j)
c	j1 = int((j1/10.0) + 1.0)
	j1 = int((j1/aii) + 1.0)
	j1 = mod(j1,nrot(k))
	if(j1.eq.0) j1 = nrot(k)
c	write(31,*)j1,k,nrot(k),nchi(k)
	do j2 = 1, nchi(k)
		phi(i) = sc(k,j1,j2)!sc = stored rotamer values
		i = i + 1
	enddo
	j = j + 1
	endif
	k = k + 1
	endif

	if(l3.eq.2) then
	if(fseq(k).eq.'P') go to 192
	phi(i) = p(j,l2,1)
	i = i + 1
	j = j + 1
192	phi(i) = p(j,l2,1)
	i = i + 1
	j = j + 1
	if(nrot(k).ne.0.) then
	j1 = p(j,l2,1)
c	j1 = int((j1/10.0) + 1.0)
	j1 = int((j1/aii) + 1.0)
	j1 = mod(j1,nrot(k))
	if(j1.eq.0) j1 = nrot(k)
	do j2 = 1, nchi(k)
		phi(i) = sc(k,j1,j2)
		i = i + 1
	enddo
	j = j + 1
	endif
	k = k + 1
	endif

	if(k.le.lres) go to 161
	else
	do i = 1, nn
	if(l3.eq.1) phi(i) = e(l1,l2,i)
	if(l3.eq.2) phi(i) = p(i,l2,1)
	enddo
	endif

c	write(31,*)(phi(ii),ii=1,nn)
c	write(*,*)(phi(ii),ii=1,nn)
c	print *,nn
c	stop
	return
	
	end
	
c***********************************************************************
