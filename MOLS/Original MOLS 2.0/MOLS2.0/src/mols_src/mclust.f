c     library name : mclust.f   

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



        subroutine mclust(icnt,rk1,kr2)
	include 'mols.par'
        common /par/ natom, ntor, nhb, ns, lres
       common /roms/rmsd(maxstr,maxstr),ni1,ni2
       common /crd/xi(maxstr,maxatm,3),spdb(maxatm),ene(maxstr),
     &  stno(maxstr)
       common /scrd/x1(maxatm),y1(maxatm),z1(maxatm),
     &   x2(maxatm),y2(maxatm),z2(maxatm)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character spdb*30,epdb*25,mdl*5,endmdl*6,ter*3,
     &   mres1(maxstr)*2,mres2(maxstr)*4,atom1*4
	dimension epdb(maxatm),eene(maxstr),iclust(maxstr,maxstr),
     &   nm(maxstr),arms(maxstr),trms(maxstr,maxstr),rgy(maxstr),rr1(3),
     &   avene(maxstr)

	integer stno,resno2(maxstr)
	mdl = 'MODEL'
	endmdl = 'ENDMDL'
	ter = 'TER'
	
	ni1=ns
	ni2=ns
c	print *,'Enter the No of Atoms'
c	read(*,*)natom
c	print *,'Enter the value of rmscutoff'
c	read(*,*)rcutinp
c	rcutinp = 0.5



	call superimpose(kr2)
24	format(3000f6.3)

	do i = 1,ni1
	eene(i) = ene(i)
	do j = 1,ni1
	trms(i,j) = rmsd(i,j)
	enddo
	enddo

c	do j1=1,ni1
c	open (unit=7, file='crms.out',status='unknown')
c	write(7,24) (rmsd(j1,j2),j2=1,ni2)
c	enddo
c	close(unit=7)

c       go to 111

c	open(unit=7,file='crms.out'
c     & ,status='old')
c104	do j1=1,ni1
c	read(7,24) (rmsd(j1,j2),j2=1,ni2)
c	enddo
c	close(unit=7)

	go to 102
	temp1 = 10
	temp2 = 0.1
	do i = 1,ni1
	do j = i+1,ni2
	if(rmsd(i,j).lt.temp1) temp1 = rmsd(i,j)
	if(rmsd(i,j).gt.temp2) temp2 = rmsd(i,j)
	enddo
	enddo
	dmin = temp1
	dmax = temp2
c	print *,dmax,dmin
c    6.762000      8.2000002E-02
	icountt = 0
	rmin = 0.0
	rmax = 7.0
	rint = 0.2
	r1 = rmin
	do k = 1,maxstr
	r2 = r1 + rint
	if(r2.gt.rmax) go to 101
	icount = 0
	do i = 1,ni1
	do j = i+1,ni2
	if(rmsd(i,j).ge.r1.and.rmsd(i,j).lt.r2) icount = icount + 1
	if(rmsd(i,j).lt.0.0) print *,i,j,rmsd(i,j)
	enddo
	enddo
cc	print *, r1,r2,icount
	icountt = icountt + icount
	r1 = r2
	enddo
101	print *,icountt
c
102	rcut = rk1
c
20	format(a30,3f8.3)	
21	format(a5,4x,i4,3x,f10.2)
c21	format(a5,5x,i4)
22	format(a6)
23      format(a3)
31      format(a4,4x,i3,2x,a2,2x,a4,i5,4x,3f8.3)	
54	format(a50)
55	format(i3,a1,i4)
56	format(a1)

	temp1 = 0.0
	itot = 0
	do k = 1,maxstr
	icutt = 0
	temp  = 10000.0
		do i =1,ni1
		    if(ene(i).lt.temp) then
			temp = ene(i)
			ic1 = i
		    endif
		enddo
c	if(temp.eq.temp1) go to 133
	if(itot.eq.ni1) go to 133
ccc	print *,k,ic1,temp
		temp1 = temp
		i2 = ic1
		ene(i2) = 9999.0
		in = 0
		in = in + 1
		itot = itot + 1
		iclust(k,in) = i2
		nm(k) = in
cc	do j2 = i2+1,ni2
	do j2 = 1,ni2
	  if(rmsd(i2,j2).lt.rcut.and.i2.ne.j2) then
		in = in + 1
		itot = itot + 1
		iclust(k,in)=j2
		nm(k) = in
		ene(j2) = 9999.0
c		print *,i2,j2
		do i1 = 1,ni1
			rmsd(i1,j2) = 9999.0
		enddo
		do i1 = 1,ni1
			rmsd(j2,i1) = 9999.0
		enddo
	  endif
	enddo
	  do j1 = 1,ni1
		rmsd(i2,j1) = 9999.0
	 enddo
	  do j1 = 1,ni1
		rmsd(j1,i2) = 9999.0
	 enddo
	enddo
133     open(unit=4,file=of3,status='unknown')
40	format(500i5)
50	format(a6,i4,2x,i4,2x,f8.2,2x,f8.2,2x,f8.2,2x,f8.2,2x,f8.2)
	do i = 1,k-1
	write(4,40) i,nm(i),(iclust(i,j),j=1,nm(i))
c	write(*,40) i,nm(i),(iclust(i,j),j=1,nm(i))
	enddo
	nc = k-1
c
c	print the closest structure no. with centroid
c	
c	write(31,*) '    cluster #','   closest structure #'
	write(4,54) 'Reperesentative structure numbers to each cluster' 
	do i1 = 1,nc
	i = iclust(i1,1)
c	write(31,*)i1, i
	write (4,55) i1,':',i
	enddo
	close(unit = 4)
c
	call mclustplot(iclust,nm,nc,icnt)
c	
        open(unit=8,file=pf4,status='unknown')
	open (unit=1, file=pf1,status='old')
c          average rmsd
	do i = 1,nc
	srms=0.0
	do j = 2,nm(i)
	i1=iclust(i,1)
	j1 = iclust(i,j)
	srms = srms + trms(i1,j1)
c	print *,trms(i1,j1),j1,i1,i
	enddo
	if(nm(i).eq.1) then
	arms(i) = srms/1.0
	else
	arms(i) = srms/(nm(i)-1)
	endif
c	print *, arms(i)
	enddo
c
c
	do i = 1,nc
	rr2 = 0.0
	dot = 0.0
	do j3 = 1,3
	rr1(j3) = 0.0
	enddo
c
	i1=iclust(i,1)
	do i2 = 1,natom
	rr1(1) = rr1(1) + xi(i1,i2,1)
	rr1(2) = rr1(2) + xi(i1,i2,2)
	rr1(3) = rr1(3) + xi(i1,i2,3)
	rr2 = rr2 + xi(i1,i2,1)**2 + xi(i1,i2,2)**2 + xi(i1,i2,3)**2
	enddo
	dot = rr1(1)**2 + rr1(2)**2 + rr1(3)**2
	rgy(i) = (rr2 - dot/float(natom)) /float(natom)
	rgy(i) = sqrt (rgy(i))	
c
c	print *,rgy(i)
	enddo
c
34	format(a66)
35	format(a37,i5)
36	format(a37,f5.2)
	write (8,35),'REMARK  SEQUENCE LENGTH            = ',lres
	write (8,36),'REMARK  RMSD CUTOFF                = ',rcut
	write (8,35),'REMARK  TOTAL STRUCTURES           = ',ni1
	write (8,35),'REMARK  NUMBER OF CLUSTERS         = ',nc
	write (8,34),'REMARK  CN    NM    AVEENE    MINENE   
     &    MAXENE      AVERMSD   RGYR'
	do i = 1,nc
	aene=0.0
	dene=0.0
	dene=0.0
	a1 = 100.0
	a2 = -100.0

	do j = 1,nm(i)
	if(eene(iclust(i,j)).lt.a1) a1 = eene(iclust(i,j))
	if(eene(iclust(i,j)).gt.a2) a2 = eene(iclust(i,j))
	aene = aene + eene(iclust(i,j))
c	print *,eene(iclust(i,j)),iclust(i,j),aene,i
	enddo
	aene = aene/float(nm(i))
ccc	dene = (a2)-(a1)
	write(8,50) 'REMARK',i,nm(i),aene,a1,a2,arms(i),rgy(i)
	avene(i) = aene
	enddo
c
c	print the co-ordinates of centroids in the pdb format
c
	do i1 = 1,nc
	i = iclust(i1,1)
		write(8,21) mdl,i1,avene(i1)
	rewind (1)
	jc = 0
	  do j=1,natom
	j4 = 0
	read(1,31) atom1,iatmno1,mres1(j),mres2(j),resno2(j),
     &  xx1,xx2,xx3
	if(mres1(j).eq.'CA') then
	j4 = j
	jc = jc + 1
	write(8,31) atom1,jc,mres1(j),mres2(j),resno2(j),
     &  (xi(i,j4,k),k=1,3)
	endif
	  enddo
c		read(1,23) ter
		write(8,23) ter
c		read(1,22) endmdl
		write(8,22) endmdl
	enddo
	close(unit=1)
	close(unit=8)
c
	go to 111
cccccccccccccheck later cccccccccccccccccc

c111	open(unit=7,file='cprop.out',status='unknown')
c50	format(i4,2x,i3,2x,f8.2,2x,f6.2,2x,f8.2)
c	do i = 1,nc
c	write(7,50) i,nm(i),aene(i),arms(i),rgy(i)
c	enddo
c	close(unit=7)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	natom  = 75

	open (unit=1, file=pf3, status='old')
	do i=1,ni1
	read(1,21) mdl,st
c	print *,md1,st
		do j=1,natom
		  read(1,20) spdb(j),(xi(i,j,k),k=1,3)
		enddo
	read(1,23) tr
	read(1,22) emdl
	enddo
	close(unit=1)

	open(unit=3,file ='c1.pdb',status='unknown')
	do i = ic1,ic1
c	print *,ic1	
	write(3,21) mdl,ic1
		do j1=1,natom
		  write(3,20) spdb(j1),(xi(i,j1,k),k=1,3)
		enddo
	write(3,23) tr
	write(3,22) emdl

	do j = 1,ni2
	if(rmsd(i,j).lt.rcut) then
	print *,j	
	write(3,21) mdl,j
		do j1=1,natom
		  write(3,20) spdb(j1),(xi(j,j1,k),k=1,3)
		enddo
	write(3,23) tr
	write(3,22) emdl
	endif
	enddo
	enddo
	close(unit=3)

111	return
	end

c***********************************************************************
	subroutine superimpose(kr2)
	include 'mols.par'
        common /par/ natom, ntor, nhb, ns, lres
       common /roms/rmsd(maxstr,maxstr),ni1,ni2
       common /crd/xi(maxstr,maxatm,3),spdb(maxatm),ene(maxstr),
     &  stno(maxstr)
       common /scrd/x1(maxatm),y1(maxatm),z1(maxatm),
     & x2(maxatm),y2(maxatm),z2(maxatm)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	integer*2 jpick(maxatm),jpick2(maxatm)
	integer i2,iorie,i3,extra
cc	real ptms(maxatm),ptms2(maxatm)
cc	logical rmsmaton,ok
	dimension atom1(maxatm),epdb(maxatm)
cc	integer ovrlnpt,ovrlnpt2,molflaga,molflagb

c local
cc	double precision kabs(3,3),kabs2(3,3),a(3),e(3),b(3,3)
cc	double precision dkabs(3,3),rotat(3,3),check1(3,3),check2(3,3)
cc	double precision xcm,ycm,zcm,allmass,det,norm,allmass2
cc	double precision xcm1,ycm1,zcm1,tmpx,tmpy,tmpz
cc	double precision small
cc	real rms,mass
cccccccccc
        real dm(4),vm(4,4),cm(3),cf(3)                                  
        real tr(3),t(3,3),q(4,4),temp(3)                                
        real dxp(maxatm,3),dxm(maxatm,3),xf(maxatm,3),xm(maxatm,3)
        real*8 sumf(3),summ(3),xx(3)
ccccccccccc
	integer level,error,ismall
	integer i,j,k,namel
	character*7 name
ccc added decleration
	character*30 spdb
	character*25 epdb
	character*6 endmdl
	character*3 ter
	character*5 model
	character*4 atom1
	character*2 mres1(maxatm),tres
	character*10 atom2
	character*4 res3,mres2(maxatm)
	character*1 res1
	character*1 res
	character*1 atname(maxatm)
	character*2 an1,an2,an3,an4
	integer natom,stno,j1,j2,resno2(maxatm),ni1,ni2
ccc
c	rmsmaton=.false.
c	molflaga=1
c	molflagb=2
	extra=0
	
	small=1.d20
	
	namel = 7
	name  = 'overlap'
ccc
c20	format(a30,3f8.3,a23,a1)	
20	format(a30,3f8.3,a25)	
21      format(a5,3x,i5,3x,f10.2)
c21	format(a5,5x,i4,4x,f8.2)
22	format(a6)
23      format(a3)
24	format(3000f6.3)
c10	format(1x,6f10.4)	
31      format(a4,4x,i3,2x,a2,2x,a4,i5,4x,3f8.3,a25)	
c	natom=75
c	ni1=20
c	ni2=20
	do i = 1,ni1
	ene(i) = 0.0
	enddo
	open (unit=1, file=pf3, status='old')
	do i=1,ni1
	read(1,21) model,stno(i),ene(i)
c	write (6,21) model,stno(i),ene(i)
	do j=1,natom
	  read(1,31) atom1(j),atomno1,mres1(j),mres2(j),resno2(j),
     &    (xi(i,j,k),k=1,3),epdb(j)
c		  read(1,20) spdb(j),(xi(i,j,k),k=1,3),epdb(j)
		enddo
	read(1,23) ter
	read(1,22) endmdl
	enddo
	close(unit=1)
cccc
c	return
cccc
	do j1=1,ni1
	do j2=1,ni2
C	print *,j1,j2

	iorie=0
c
	do j=1,natom
	tres=mres1(j)
cccccccccccccccccccccccccccccccccccccccc
c	if(tres(1:1).ne.'H'.and.tres(1:1).ne.'L') then
c	if(tres(1:2).eq.'N '.or.tres(1:2).eq.'CA'.or.
c     &  tres(1:2).eq.'C '.or.tres(1:2).eq.'O ') then
	if(kr2.eq.3) then
	if(tres(1:2).eq.'CA') then
	iorie=iorie+1
		x1(iorie)=xi(j1,j,1)
		y1(iorie)=xi(j1,j,2)
		z1(iorie)=xi(j1,j,3)
		xf(iorie,1)=xi(j1,j,1)
		xf(iorie,2)=xi(j1,j,2)
		xf(iorie,3)=xi(j1,j,3)
	endif
	else if(kr2.eq.2) then
	if(tres(1:2).eq.'N '.or.tres(1:2).eq.'CA'.or.
     &  tres(1:2).eq.'C '.or.tres(1:2).eq.'O ') then
	iorie=iorie+1
		x1(iorie)=xi(j1,j,1)
		y1(iorie)=xi(j1,j,2)
		z1(iorie)=xi(j1,j,3)
		xf(iorie,1)=xi(j1,j,1)
		xf(iorie,2)=xi(j1,j,2)
		xf(iorie,3)=xi(j1,j,3)
	endif
	else if(kr2.eq.1) then
	if(tres(1:1).ne.'H'.and.tres(1:1).ne.'L') then
	iorie=iorie+1
		x1(iorie)=xi(j1,j,1)
		y1(iorie)=xi(j1,j,2)
		z1(iorie)=xi(j1,j,3)
		xf(iorie,1)=xi(j1,j,1)
		xf(iorie,2)=xi(j1,j,2)
		xf(iorie,3)=xi(j1,j,3)
	endif
	endif
cccccccccccccccccccccccccccccccccccccccccc
	enddo
c
	iorie=0
	do j=1,natom
	tres=mres1(j)
cccccccccccccccccccccccccccccccccccccccccc
c	if(tres(1:1).ne.'H'.and.tres(1:1).ne.'L') then
c	if(tres(1:2).eq.'N '.or.tres(1:2).eq.'CA'.or.
c     &  tres(1:2).eq.'C '.or.tres(1:2).eq.'O ') then
	if(kr2.eq.3) then
	if(tres(1:2).eq.'CA') then
	iorie=iorie+1
		x2(iorie)=xi(j2,j,1)
		y2(iorie)=xi(j2,j,2)
		z2(iorie)=xi(j2,j,3)
		xm(iorie,1)=xi(j2,j,1)
		xm(iorie,2)=xi(j2,j,2)
		xm(iorie,3)=xi(j2,j,3)
	endif
	else if(kr2.eq.2) then
	if(tres(1:2).eq.'N '.or.tres(1:2).eq.'CA'.or.
     &  tres(1:2).eq.'C '.or.tres(1:2).eq.'O ') then
	iorie=iorie+1
		x2(iorie)=xi(j2,j,1)
		y2(iorie)=xi(j2,j,2)
		z2(iorie)=xi(j2,j,3)
		xm(iorie,1)=xi(j2,j,1)
		xm(iorie,2)=xi(j2,j,2)
		xm(iorie,3)=xi(j2,j,3)
	endif
	else if(kr2.eq.1) then
	if(tres(1:1).ne.'H'.and.tres(1:1).ne.'L') then
	iorie=iorie+1
		x2(iorie)=xi(j2,j,1)
		y2(iorie)=xi(j2,j,2)
		z2(iorie)=xi(j2,j,3)
		xm(iorie,1)=xi(j2,j,1)
		xm(iorie,2)=xi(j2,j,2)
		xm(iorie,3)=xi(j2,j,3)
	endif
	endif
ccccccccccccccccccccccccccccccccccccccccccccc
	enddo
	ifix = iorie
	imov = ifix
c-----------------------------------------------------
c --- initialize all --                                                 
        do i=1,3                                                          
         sumf(i)=0                                                      
         summ(i)=0                                                      
        end do                                                            
        call filmat(4,4,q,0)                                              
                                                                        
c --- sum up all coordinates (in dble precision) to find centre ---     
       do k=1,imov                                                       
         do i=1,3                                                       
            sumf(i)=sumf(i)+xf(k,i)                                     
            summ(i)=summ(i)+xm(k,i)                                     
         end do                                                         
       end do                                                            
       do i=1,3                                                          
         cm(i)=sngl(summ(i)/imov)                                       
         cf(i)=sngl(sumf(i)/imov)                                       
         tr(i)=cf(i)-cm(i)                                              
       end do                                                            
                                                                        
c --- create coordinate differences delta x plus (dxp) and minus (dxm)  
       do k=1,imov                                                       
         do j=1,3                                                       
            dxm(k,j)=xm(k,j)-cm(j)-(xf(k,j)-cf(j))                      
            dxp(k,j)=xm(k,j)-cm(j)+(xf(k,j)-cf(j))                      
         end do                                                         
       end do                                                            
                                                                        
c --- fill upper triangle of (symmetric) quaternion matrix --           
       do k=1,imov                                                       
c ---    diags are sums of squared cyclic coordinate differences        
         q(1,1)=q(1,1)+dxm(k,1)**2+dxm(k,2)**2+dxm(k,3)**2              
         q(2,2)=q(2,2)+dxp(k,2)**2+dxp(k,3)**2+dxm(k,1)**2              
         q(3,3)=q(3,3)+dxp(k,1)**2+dxp(k,3)**2+dxm(k,2)**2              
         q(4,4)=q(4,4)+dxp(k,1)**2+dxp(k,2)**2+dxm(k,3)**2              
c ---    cross differences                                              
         q(1,2)=q(1,2)+dxp(k,2)*dxm(k,3)-dxm(k,2)*dxp(k,3)              
         q(1,3)=q(1,3)+dxm(k,1)*dxp(k,3)-dxp(k,1)*dxm(k,3)              
         q(1,4)=q(1,4)+dxp(k,1)*dxm(k,2)-dxm(k,1)*dxp(k,2)              
         q(2,3)=q(2,3)+dxm(k,1)*dxm(k,2)-dxp(k,1)*dxp(k,2)              
         q(2,4)=q(2,4)+dxm(k,1)*dxm(k,3)-dxp(k,1)*dxp(k,3)              
         q(3,4)=q(3,4)+dxm(k,2)*dxm(k,3)-dxp(k,2)*dxp(k,3)              
       end do                                                            
c --- fill the rest by transposing it onto itself                       
       call trpmat(4,q,q)                                                
c --- orthogonalization by jacobi rotation = solution of EV -problem -- 
       nx=4                                                               
       nxs=4                                                              
       call jacobi(q,nx,nxs,dm,vm,nmrot)                                   
c --- sort eigenvectors after eigenvalues, descending --                
       call eigsrt(dm,vm,nx,nxs)                                           
c --- the smallest eigenvector contains best fit srs                    
       rms=sqrt(abs(dm(4)/imov))                                        

	rmsd(j1,j2)=rms
cc	print *,j1,j2,rms,rmsd(j1,j2)
cc     	    print *,' rms is = ',j1,'-',j2,'=',rms
	rmsd(j1,j2)=rms
c
	enddo
	enddo
c
c	do j1=1,ni1
c	open (unit=7, file='rand.out',status='unknown')
c	write(7,24) (rmsd(j1,j2),j2=1,ni2)
c	enddo
c	close(unit=7)
	return
	end
c**********************************
	subroutine mclustplot(iclust,nm,kk,icnt)
	include 'mols.par'
        common /par/ natom, ntor, nhb, ns, lres
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	integer iclust(maxstr,maxstr),nm(maxstr),kk
	integer ic(maxatm),iz(maxatm)
	character aster*2,space*2,fp(maxatm,maxatm)*2,ffp(maxatm,maxatm)*2,
     &   hll(maxstr)*2,vl*1,hl*2	
	data aster/' *'/
	data space/'  '/
	data vl/'|'/
	data hl/'--'/
c
40	format(100(i5))
c	print *,kk
c	do i = 1,kk
c	write(*,40) i,nm(i),(iclust(i,j),j=1,nm(i))
c	enddo
c	close(unit=4)
c
ccc	write(31,*)'Enter the interval to draw cluster plot'
ccc	read *, in
c	read(30,*) in
c	in = ns/10
	in = icnt
	ni = ns/in
c	print *,ni,ns,in,kk
	
	do i = 1,kk
	ic(i) = 0
	iz(i) = 9999
	enddo
	do i = 1,maxstr
	hll(i) = hl
	enddo
	open (unit=22,file=of4,status='unknown')

	k = 1
c	kk = kk-1
1112	do i = 1,kk
	do j = 1,nm(i)
c	print *,i,j,iclust(i,j),k*in,(k-1)*in,kk,nm(i)
	if(iclust(i,j).le.k*in.and.iclust(i,j).gt.((k-1)*in).and.
     $  iz(i).ne.0) then
	ic(k) = ic(k) + 1
	iz(i) = 0
	go to 1111
	endif
	enddo
1111	enddo
	if(k.le.ni) then
	k = k + 1
	go to 1112
	endif
c
c	do i = 1,k
c	write(6,*) i*in, ic(i)
c	enddo
c
c	do i = 1,k
c	mp = ic(i)
c	write(*,'(100a2)') (aster, j = 1,mp)
c	enddo
	
	do i = 1,k
	mp = ic(i)
	do j = 1,mp
	fp(i,j) = aster
	enddo
c	write(*,'(100a2)') (fp(i,j), j = 1,mp)
	enddo
c
	kt = 0
	ktv = 0
	kth = 0
	do i = 1,k
	if(kt.lt.ic(i)) kt = ic(i)
	enddo
	kth = kt
	ktv = kt
	if(ktv.lt.k) ktv=k
	if(kt.lt.k) kt = k
c	ktv=ni
c	print *,kt

	do i = 1,kt
	do j = 1,kt
	if(fp(i,j).ne.aster) fp(i,j) = space
	enddo
c	write(*,'(100a2)') (fp(i,j), j = 1,mp)
	enddo

	do i = 1,kt
	do j = 1,kt
	ffp(j,i) = fp(i,j)
	enddo
	enddo
90	format(a60)
91	format(1x,a1,100a2)
92	format(1x,100a2)
	write(22,90) 'Structure Numbers Vs. Number of new Clusters Plot'
c	do i = fj(1),1,-1
	do i = kth,1,-1
	write(22,91) vl,(ffp(i,j),j=1,kt)
	enddo

c	write(22,92) (hll(j),j=1,nm(1))
	write(22,92) (hll(j),j=1,ktv)

93	format(a60)
94	format(5x,i4,a3,i4,25x,i4)
	write(22,93) 'Range of Struct. numbers		No. of new clusters'	
	do i = 1,ni
	write(22,94) (((i-1)*in)+1),' - ',i*in,ic(i)
	enddo

	if(mod(ns,in).ne.0)write(22,94) (((i-1)*in)+1),' - ',ns,ic(i)
	close(unit=22)
	return
	end

c********************************************************************

