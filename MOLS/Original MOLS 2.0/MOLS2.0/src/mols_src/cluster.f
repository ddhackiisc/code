c     library name : cluster.f   

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



c
c	Clustering of structures obtained by MOLS procedure
c
c	program cluster$main
	subroutine cluster(icnt)
	
	include 'mols.par'
        common /par/ natom, ntor, nhb, ns, lres
	common/crds/xi(maxstr,maxatm,3)
	common/crdr/yi(maxstr,maxatm,3)
	common/acrd/x1(maxatm),x2(maxatm),y1(maxatm),y2(maxatm),
     &   z1(maxatm),z2(maxatm)
	common/ccrd/x(maxatm,maxatm),y(maxatm,maxatm),z(maxatm,maxatm)
	common/dev/rmsd(maxstr,maxstr),rrmsd(maxstr,maxstr)
	common/parm/ni1,ni2,mres1
	common/rd/chainid(maxstr),pname(maxstr)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4

	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character endmdl*6,ter*3,model*5,atom1*4,atom*4,
     &  atom2*10,res3*4,res*1,chainn*1,resn*2,xpdb*42,xchain*5
	character spdb(maxstr)*30,epdb(maxstr)*23,atname(maxstr)*1,
     &  mres1(maxstr)*2,mres2(maxstr)*4, yn*1
	character*5 chainid
	character*42 pname

       integer natom,i,j,i2,j2,i3,j3,ni1,ni2,k,temp1
     &  ,clust(maxstr,maxstr)
       integer stno(maxstr),resno2(maxstr),flag,nspn(maxstr)
     &  ,nn(maxstr),clust1(maxstr,maxstr),nn1(maxstr),i4,j4,nm(maxstr),
     & dn,nclust(maxstr,maxstr),fclust(maxstr,maxstr),nj(maxstr),tc,
     & fj(maxstr),ffclust(maxstr,maxstr),nm1(maxstr),
     & nclust1(maxstr,maxstr),iseed,naene(maxstr),del(maxstr),
     & naa(maxstr),resno,iclose(maxatm),ifix1

	real xi,yi,ene(maxstr),r,rr,rmsd,rand(maxstr,maxstr),avrg,
     & sum,cgcut,rgcut,gcut,n,std,sigma,spn(maxstr),spn1(maxstr),temp,
     & lcut(maxatm),lcut1,lcut2,cpn(maxatm),x1,x2,y1,y2,z1,z2,
     & xsum(maxatm),ysum(maxatm),zsum(maxatm),nlcut(maxatm),nlcut1,
     & nlcut2,cpn1(maxatm),fcpn(maxatm),flcut(maxatm),srms,rg(maxatm),
     & sene,aene(maxstr),arms(maxatm),trms,rrms,arrms(maxatm),
     & aene1(maxstr),fcpd(maxatm),cc1,cc2,cc3,cpd(maxatm),
     & crms(maxatm), crrms(maxatm), pgcut(20)

	real*8 rr1(3),rr2,dot,rgy(maxatm)
     
cc	these are precalculated values for initial global cluster cutoff
	data (pgcut(i),i=5,15) /0.2352, 0.2773, 0.3091, 0.3267,
     &  0.3405, 0.3507, 0.3624, 0.3780, 0.3889, 0.3997, 0.3983/

cc	write(31,*)'Enter the no. of parameters'
cc	read *, ntor
cc	write(31,*)'Enter the no. of atoms'
cc	read *, natom
cc	write(31,*)'Enter the no. of optimal structures'
cc	read *, ns
cc	write(31,*)'Enter the no. residues'
cc	read *, lres

c	call corgen
ccccccccccccccccccc
c	return
cccccccccccccccccccccc
	
20	format(a30,3f8.3,a25)	
21	format(a5,4x,i4,3x,f10.2)
c21	format(a5,5x,i4,f8.2)
22	format(a6)
23      format(a3)
24	format(100f6.3)
31      format(a4,4x,i3,2x,a2,2x,a4,i5,4x,3f8.3,a25)	
32      format(30x,3f8.3)
33	format(a6,i4,i4,8f8.2)	
c34	format(a62)
34	format(a79)
35	format(a37,i5)
36	format(a37,f5.2)
c27	format(a1,i3,a1,500i5)
27	format(a1,i3,a1,/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5
     &  /15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5
     &  /15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5
     &  /15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5/15i5
     &  /15i5)

50	format(a5,1x,i4)
51	format(a42)
52 	format(a4)
53	format(a4,9x,a2,6x,a1,i4,4x,f8.3,f8.3,f8.3)
54	format(a42)
55	format(i3,a1,i4)
56	format(a1)
	
	sum=0.0
	avrg=0.0
	gcut=0.0
	ccgut=0.0
	rgcut=0.0
	n=0
	ni1= ns
	ni2 = ns

	gcut = pgcut(lres)

	if(lres.gt.15) then
	 write(31,*) 'you exceeded the upper limit of sequence
     & length, so you need to calculate the initial global cluter cutoff
     & value by uncommending the lines in the cluster routine or yourself
     & using non-homologus protein chains'
	write(31,*)'do you want to continue?(y or n)'
	write(31,*) 'if yes program will take the default value'
	write(31,*) 'if no you have to supply the global cutoff value'
	read(5,56) yn
	if(yn.eq.'n') stop
	if(yn.eq.'y') gcut = pgcut(15)
	endif

c
c	read co-ordinates of structures to be clustered
c
	open (unit=1, file=pf3,status='old')
	do i=1,ni1
		read(1,21) model,stno(i),ene(i)
	  do j=1,natom
		  read(1,31) atom1,atomno1,mres1(j),mres2(j),resno2(j),
     & (xi(i,j,k),k=1,3),epdb(j)
	  enddo
		read(1,23) ter
		read(1,22) endmdl
	enddo
	close(unit=1)
c
c	calculate average and sigma for rmsd distribution of structures to be clustered
c
	sum=0.0
	avrg=0.0
	flag=2
	std=0.0
	n=0
	do i = 1,ni1
	do j = 1,ni2
	call rootms(i,j,r,rr,flag)
	rmsd(i,j) = r
	rrmsd(i,j) = rr
c	print *,i,j,r,rr
	if(j.gt.i) then
	sum = sum + r
	n = n + 1
	endif
	enddo
	enddo
	
 
ccc	avrg = sum/n
ccc	do i = 1,ni1
ccc	do j = i+1,ni2
ccc	std = std + (rmsd(i,j) - avrg)**2
ccc	enddo
ccc	enddo
ccc	sigma = sqrt(std/n)
ccc	cgcut = avrg - sigma
c	write(31,*) avrg,'sigma =',sigma,cgcut,n
c
c	set default cluster cutoff: whichever is small cgcut or rgcut?
c
cc	if(rgcut.lt.cgcut) then
cc	gcut = rgcut
cc	else 
cc	gcut = cgcut
cc	endif

cc	gcut = rgcut
	write(31,*) 'global cluster cutoff = ',gcut

c
c	print the rmsd distribution of random structures and 
c	the structures being clustered
c
cc		open (unit=7, file='rand.out',status='unknown')
cc	do j1=1,nj1
cc		write(7,24) (rand(j1,j2),j2=1,nj2)
cc	enddo
cc	close(unit=7)

cc		open (unit=8, file='rms.out',status='unknown')
cc	do j1=1,ni1
cc		write(8,24) (rrmsd(j1,j2),j2=1,ni2)
cc	enddo
cc	close(unit=8)

cc	do j1=1,ni1
cc		open (unit=9, file='rrms.out',status='unknown')
cc		write(9,24) (rmsd(j1,j2),j2=1,ni2)
cc		enddo
cc	close(unit=9)

	do i = 1,ni1
		spn(i) = 0.0
		spn1(i) = 0
		nspn(i) = i
		naene(i) = i
		aene1(i) = 0.0
		aene(i) = 0.0
		clust(i,j) = 0
		clust1(i,j) = 0
		nclust(i,j) = 0
		nclust1(i,j) = 0
		fclust(i,j) = 0
		ffclust(i,j) = 0
		nn(i) = 0
		nn1(i) = 0
		nm(i) = 0
		nm1(i) = 0
	enddo
	do i = 1,100
		do j = 1,100
		cpn(i) = 0.0
		lcut(i) = 0.0
		nlcut(i) = 0.0
		cpn1(i) = 0.0
		cpd(i) = 0.0
		fcpd(i) = 0.0
		del(i) = 0
		arms(i) = 0.0
		iclose(i) = 0
		enddo
	enddo
	temp = 0
	temp1 = 0

c
c	Compute the structure packing number
c
	do i = 1,ni1
	do j = 1,ni2
	if(rmsd(i,j).le.gcut) then
		spn(i) = spn(i) + exp(-2.0*((rmsd(i,j)/gcut)**2))
		spn1(i) = spn(i)
	endif
	enddo
c	write(31,*)spn(i),i
	enddo

c
c	Arrange the structures in order of decreasing packing number
c
	do i = 1,ni1
	k = ni1 - i
	do j = 1,k
		if(spn1(j).ge.spn1(j+1)) go to 100
		temp=spn1(j)
		temp1=nspn(j)
		spn1(j)=spn1(j+1)
		nspn(j)=nspn(j+1)
		spn1(j+1)=temp
		nspn(j+1)=temp1
100	enddo
	enddo
	
	do i = 1,ni1
c	write(31,*) spn1(i),nspn(i)
	enddo


c
c	Select the center with highest structure packing number
c
	mm = 0
	do i = 1,ni1
	i1 = nspn(i)
c	check the center with previous cluster
	do l1 = 1,ni1
	do l2 = 1,ni2
	if(clust(l1,l2).eq.i1) go to 101
	enddo
	enddo

	mm =  mm + 1
	nn(mm) = 1
	clust(mm,nn(mm)) = i1
	do j = 1,ni2
	if(rmsd(i1,j).le.gcut.and.i1.ne.j) then
	nn(mm) =  nn(mm) + 1
	clust(mm,nn(mm)) = j
c	write(31,*)i1,j
	endif
	enddo
	if(nn(mm).lt.2) mm = mm - 1
101	enddo

	write(31,*)'No. of clusters = ',mm

cccc	write(31,*)'The initial cluster id and members ids using global cutoff'
cccc	do i = 1,mm
cccc	write(31,*)(clust(i,j),j=1,nn(i))
cccc	enddo
c
c	Compute the local cluster cutoff
c


c	do i=1,ni1
c	do j=i+1,ni2
	do i=1,mm
	lcut1=0.0
	lcut2=0.0
	do i2 = 1,nn(i)
	do j2 = i2+1,nn(i)
	i1 = clust(i,i2)
	j1 = clust(i,j2)
	lcut1= lcut1 + (rmsd(i1,j1)**2)*(exp(-2.0*((rmsd(i1,j1)/gcut)**2)))
	lcut2= lcut2 + (exp(-2.0*((rmsd(i1,j1)/gcut)**2)))
	enddo
	enddo
	lcut(i)=lcut1/lcut2
	lcut(i) = sqrt(lcut(i))
	enddo

cccc	do i = 1,mm
cccc	write(31,*) 'local cluster cutoff for cluster no: ',i, '= ',lcut(i)
cccc	enddo

c
c	Select the members to the centers using local cutoff
c
	mm1 = 0
	do i = 1,mm
	i1 = clust(i,1)
	mm1 =  mm1 + 1
	nn1(mm1) = 1
	clust1(mm1,nn1(mm1)) = i1
	do j = 1,ni2
	if(rmsd(i1,j).le.lcut(i).and.i1.ne.j) then
	nn1(mm1) =  nn1(mm1) + 1
	clust1(mm1,nn1(mm1)) = j
c	write(31,*)i1,j
	endif
	enddo
102	enddo

cccc	write(31,*)'The initial cluster id and members ids using local cutoff'
cccc	do i = 1,mm1
cccc	write(31,*)(clust1(i,j),j=1,nn1(i))
cccc	enddo

c
c	Compute the cluster packing number, which is used in the comparision
c	of next step
c

	do i = 1,mm1
	do i2 = 1,nn1(i)
	do j2 = i2+1,nn1(i)
	i1 = clust1(i,i2)
	j1 = clust1(i,j2)
	cpn(i) = cpn(i) + (exp(-2.0*((rmsd(i1,j1)/lcut(i))**2)))
	enddo
	enddo
	cpn(i) = (2.0/(nn1(i)-1))*cpn(i)
	enddo

cccc	write(31,*) 'cluster packing number for the initial clusters'
cccc	do i = 1,mm1
cccc	write(31,*) i, cpn(i)
cccc	enddo

c
c	Obtain the centroid by averaging the cluster structures
c
	flag = 2
	do i = 1,mm1
	i1 = clust1(i,1)
	do i4 = 1, lres
		xsum(i4) = 0.0
		ysum(i4) = 0.0
		zsum(i4) = 0.0
	enddo

	do j = 2,nn1(i)
	j1 = clust1(i,j)
	call rootms(i1,j1,rms,rr,flag)
	do j4 = 1, lres
		xsum(j4) = xsum(j4) + x2(j4)
		ysum(j4) = ysum(j4) + y2(j4)
		zsum(j4) = zsum(j4) + z2(j4)
	enddo
	enddo
	do j4 = 1, lres
		xsum(j4) = xsum(j4) + x1(j4)
		ysum(j4) = ysum(j4) + y1(j4)
		zsum(j4) = zsum(j4) + z1(j4)
		x(i,j4) = xsum(j4)/nn1(i)
		y(i,j4) = ysum(j4)/nn1(i)
		z(i,j4) = zsum(j4)/nn1(i)
	enddo
	enddo

c	do i = 1,mm1
c	do j = 1,lres
c	write(31,*)x(i,j),y(i,j),z(i,j)
c	enddo
c	enddo

c
c	Compute the rmsd between the cenroids and all structers and
c	pick up members with less than cluster cutoff
c
	flag = 4
	do i = 1,mm1
cc	temp = 10
	do j = 1,ni2
	call rootms(i,j,rms,rr,flag)
	if (rms.lt.lcut(i)) then
ccc	if (rms.lt.gcut) then

c	find the closest structure to the centroid
cc	if (rms.lt.temp) then
cc	temp = rms
cc	iclose(i) = j
cc	endif

	nm(i) = nm(i) + 1
	nclust(i,nm(i)) = j
	endif
	enddo
	enddo

cccc	write(31,*)'The clusters id and members ids using centroid'
cccc	do i = 1,mm1
cccc	write(31,*)(nclust(i,j),j=1,nm(i))
cccc	enddo

c
c	Comput the new local cluster cutoff
c
	do i = 1,mm1
	nlcut1 = 0.0
	nlcut2 = 0.0
	ifix = 0
	do i2 = 1,nm(i)
	do j2 = i2+1,nm(i)
	i1 = nclust(i,i2)
	j1 = nclust(i,j2)
	nlcut1 = nlcut1 + (rmsd(i1,j1)**2)*(exp(-2.0*((rmsd(i1,j1)/gcut)**2)))
	nlcut2 = nlcut2 + (exp(-2.0*((rmsd(i1,j1)/gcut)**2)))
	enddo
	enddo
	nlcut(i) = nlcut1/nlcut2
	nlcut(i) = sqrt(nlcut(i))
	enddo

cccc	do i = 1,mm1
cccc	write(31,*) 'new local cluster cutoff for cluster No: ',i,nlcut(i)
cccc	enddo


c
c	Select the members to the centers using new local cutoff
c

cc	flag = 4
cc	do i = 1,mm1
cc	do j = 1,ni2
cc	call rootms(i,j,rms,rr,flag)
check	if (rms.lt.nlcut(i)) then
cc	if (rms.lt.gcut) then
cc	nm1(i) = nm1(i) + 1
cc	nclust1(i,nm1(i)) = j
cc	endif
cc	enddo
cc	enddo

cc	write(31,*)'The clusters id and members ids using newlocal cutoff and centroid'
cc	do i = 1,mm1
cc	write(31,*)(nclust1(i,j),j=1,nm1(i))
cc	enddo

c
c	Compute the new cluster packing number
c

	do i = 1,mm1
cc	do i2 = 1,nm1(i)
cc	do j2 = i2+1,nm1(i)
cc	i1 = nclust1(i,i2)
cc	j1 = nclust1(i,j2)
	do i2 = 1,nm(i)
	do j2 = i2+1,nm(i)
	i1 = nclust(i,i2)
	j1 = nclust(i,j2)
	cpn1(i) = cpn1(i) + (exp(-2.0*((rmsd(i1,j1)/nlcut(i))**2)))
	enddo
	enddo
cc	cpn1(i) = (2.0/(nm1(i)-1))*cpn1(i)
	cpn1(i) = (2.0/(nm(i)-1))*cpn1(i)
	enddo
	
cccc	write(31,*) 'new cluster packing number '
cccc	do i = 1,mm1
cccc	write(31,*) i, cpn1(i)
cccc	enddo

c
c	Update the clusters
c

	do i = 1,mm1
	if (cpn(i).gt.cpn1(i)) then
		do j = 1,nn1(i)
		fclust(i,j) = clust1(i,j)
		enddo
		nj(i) = nn1(i)
		fcpn(i) = cpn(i)
		flcut(i) = lcut(i)
	else
		do j = 1,nm(i)
		fclust(i,j) = nclust(i,j)
		enddo
		nj(i) = nm(i)
		fcpn(i) = cpn1(i)
		flcut(i) = nlcut(i)
	endif
	enddo
cccc	write(31,*) 'cluster ids and members ids - after the refinement'
cccc	do i = 1,mm1
cccc	write(31,*)i,':',(fclust(i,j),j=1,nj(i))
cccc	enddo

c
c	Compute cluster packing density
c

	do i = 1,mm1
	cpd(i) = fcpn(i)/nj(i)
c	write(31,*)i,fcpn(i),nj(i),cpd(i)
	enddo

c
c	Compute the rmsd berween the centroids of clusters and eleminate
c	identical clusters
c

	flag = 5
	tc = 0
	dn = 0
	do i = 1,mm1
	do i4 = 1,ni1
		if(i.eq.del(i4)) go to 131
	enddo
	do j = i+1, mm1
	call rootms(i,j,rms,rr,flag)

	if(rms.lt.flcut(i).and.rms.lt.flcut(j)) then
c		write(31,*)'clusters',i,j,' overlap',flcut(i),flcut(j),rms
		if(cpd(i).lt.cpd(j)) then
			tc = tc + 1
			do j1 = 1,nj(j)
				ffclust(tc,j1) = fclust(j,j1)
			enddo
			fj(tc) = nj(j)
			fcpd(tc) = cpd(j)
			dn = dn +1
			del(dn) = i
		else
			tc = tc + 1
			do j1 = 1,nj(i)
				ffclust(tc,j1) = fclust(i,j1)
			enddo
			fj(tc) = nj(i)
			fcpd(tc) = cpd(i)
			dn = dn +1
			del(dn) = j
c	write(31,*)i,del(dn),dn
		endif	
	go to 131
c	write(31,*)i,j,flcut(i),flcut(j),rms,fcpn(i),fcpn(j)
	endif
	enddo
		tc = tc + 1
		do j1 = 1,nj(i)
			ffclust(tc,j1) = fclust(i,j1)
		enddo
		fj(tc) = nj(i)
		fcpd(tc) = cpd(i)
131	enddo

c	write(31,*) 'final cluster ids and members ids - after the refinement'
c	open(unit=10,file='cluster.out',status='unknown')
	do i = 1,tc
c	write(31,*)i,':',(ffclust(i,j),j=1,fj(i))
c	write(10,27) '#',i,':',(ffclust(i,j),j=1,fj(i))
	enddo

c	find the nearest structure to the cluster centroid	
ccccccc
	flag = 4
	do i = 1,tc
	temp = 10
	do j1 = 1,fj(i)
	j = ffclust(i,j1)
	call rootms(i,j,rms,rr,flag)
	if (rms.lt.temp) then
	temp = rms
	iclose(i) = j
	endif
	enddo
	enddo
ccccccc

c
c	Calculate the cluster properties
c

c
c	Compute the radius of gyration for each clusters
c
c	write(31,*)'  CN  NM   Energy  PACKD   CUTOF   RMSD    RRMSD   RGYR'
	flag = 4
	do i = 1,tc
	srms = 0.0
	sene = 0.0
	trms = 0.0
	rrms = 0.0
	rr2 = 0.0
	dot = 0.0
	do j4 = 1,3
	rr1(j4) = 0.0
	enddo

	do j = 1,fj(i)
	j1 = ffclust(i,j)
	call rootms(i,j1,rr,rms,flag)

	do i2 = 1,lres
	rr1(1) = rr1(1) + x2(i2)
	rr1(2) = rr1(2) + y2(i2)
	rr1(3) = rr1(3) + z2(i2)
	rr2 = rr2 + x2(i2)**2 + y2(i2)**2 + z2(i2)**2 
	enddo

	srms = srms + rms**2
	trms = trms + rms
	rrms = rrms + rr
	sene = sene + ene(j1)
	enddo
	dot = rr1(1)**2 + rr1(2)**2 + rr1(3)**2
	rgy(i) = (rr2 - dot/float(lres*fj(i))) /float(lres*fj(i))
	rgy(i) = dsqrt (rgy(i))	

	crms(i) = trms/fj(i)
	crrms(i) = rrms/fj(i)
	aene(i) = sene/fj(i)
	aene1(i) = aene(i)
cc	rg(i) = sqrt(srms/fj(i))

c	write (6,33), i,fj(i),aene(i),cpd(i),flcut(i),arms(i),arrms(i),rgy(i)
	enddo

	do i = 1,tc
	ifix1 = 0
	do i2 = 1,fj(i)
	i1 = ffclust(i,i2)
	do j = i2+1,fj(i)
	j1 = ffclust(i,j)
	ifix1 = ifix1 + 1
	arms(i) = arms(i) + rrmsd(i1,j1)
	arrms(i) = arrms(i) + rmsd(i1,j1)
c	write(31,*) i1,j1,rmsd(i1,j1),rrmsd(i1,j1)
	enddo
	enddo
	arms(i) = arms(i)/ifix1
	arrms(i) = arrms(i)/ifix1
	enddo

c	do i = 1,tc
c	write (6,33), i,fj(i),aene(i),cpd(i),flcut(i),arms(i),arrms(i),rgy(i)
c	enddo


c
c	Arrange the clusters in order of increasing energy number
c
	do i = 1,tc
	k = tc - i
	do j = 1,k
		if(aene1(j).le.aene1(j+1)) go to 200
		temp=aene1(j)
		temp1=naene(j)
		aene1(j)=aene1(j+1)
		naene(j)=naene(j+1)
		aene1(j+1)=temp
		naene(j+1)=temp1
200	enddo
	enddo

	open(unit=11,file=pf4,status='unknown')
	open(unit=10,file=of3,status='unknown')
	write (11,35),'REMARK  SEQUENCE LENGTH            = ',lres
	write (11,36),'REMARK  GLOBAL RRMSD CUTOFF        = ',gcut
	write (11,35),'REMARK  TOTAL STRUCTURES           = ',ni1
	write (11,35),'REMARK  NUMBER OF CLUSTERS         = ',tc


cc	write (11,34),'REMARK  CN  NM  Energy    PACKD   CUTOF   RMSD    RRMSD   RGYR'
	write (11,34),'REMARK  CN  NM  Energy    PACKD   CUTOF   RMSD    RRMSD   RGYR
     &    CRMS    CRRMS'
cc	write (6,34),'REMARK  CN  NM  Energy    PACKD   CUTOF   RMSD    RRMSD   RGYR'
	write (31,34),'REMARK  CN  NM  Energy    PACKD   CUTOF   RMSD    RRMSD   RGYR
     &    CRMS    CRRMS'
	do i1 = 1,tc
	i = naene(i1)
cc	write (11,33),'REMARK',i1,fj(i),aene(i),fcpd(i),flcut(i),arms(i),arrms(i),rgy(i)
	write (11,33),'REMARK',i1,fj(i),aene(i),fcpd(i),flcut(i),arms(i),
     &   arrms(i),rgy(i),crms(i),crrms(i)
c	write (6,33),'REMARK',i1,fj(i),aene(i),fcpd(i),flcut(i),arms(i),arrms(i),rgy(i)
	write (31,33),'REMARK',i1,fj(i),aene(i),fcpd(i),flcut(i),arms(i),
     &    arrms(i),rgy(i),crms(i),crrms(i)
40	format(500i5)
c	write(10,27) '#',i1,':',(ffclust(i,j),j=1,fj(i))
	write(10,40) i1,fj(i),(ffclust(i,j),j=1,fj(i))
cc	write(31,*) crms(i),crrms(i)
	enddo

c
c	print the co-ordinates of centroids in the pdb format
c
	open (unit=1, file=pf3,status='old')

	do i1 = 1,tc
	j4 = 0
	i = naene(i1)
		read(1,21) model,stno(i),ene(i)
		write(11,21) model,stno(i),aene(i)
	  do j=1,natom
		  read(1,31) atom1,atomno1,mres1(j),mres2(j),resno2(j),
     &  (xi(i,j,k),k=1,3),epdb(j)
	if(mres1(j).eq.'CA') then
	j4 = j4 + 1
	  write(11,31) atom1,j4,mres1(j),mres2(j),j4,
     &  x(i,j4),y(i,j4),z(i,j4),epdb(j)
	endif
	  enddo
		read(1,23) ter
		write(11,23) ter
		read(1,22) endmdl
		write(11,22) endmdl
	enddo
	close(unit=1)
	close(unit=11)

c
c	print the closest structure no. with centroid
c	
c	write(31,*) '    cluster #','   closest structure #'
	write(10,54) 'closest structure numbers to each clusters' 
	do i1 = 1,tc
	i = naene(i1)
c	write(31,*)i1, iclose(i)
	write (10,55) i1,':',iclose(i)
	enddo
	close(unit=10)

	call clusterplot(ffclust,fj,tc,icnt)	
2222	return
c2222	stop
	end


c********************************************************************
c********************************************************************
	subroutine rootms(s1,s2,rrms,rms,flag)

	include 'mols.par'
        common /par/ natom, ntor, nhb, ns, lres
	common/crds/xi(maxstr,maxatm,3)
	common/crdr/yi(maxstr,maxatm,3)
	common/acrd/x1(maxatm),x2(maxatm),y1(maxatm),y2(maxatm)
     &   ,z1(maxatm),z2(maxatm)
	common/ccrd/x(maxatm,maxatm),y(maxatm,maxatm),z(maxatm,maxatm)
	common/parm/ni1,ni2,mres1
cc	common/dev/rmsd(100,100),rrmsd(100,100)
cc	integer nsatm,jpick(maxatm),jpick2(maxatm),j1,j2,ifix
cc	real ptms(maxstr),ptms2(maxstr),dis,diss,rms,rrms
	real x1,y1,z1,x2,y2,z2,c
	real*8 rr1(3),rr21,rr22,dot1,dot2,rg1,rg2,rr2
	integer ovrlnpt,ovrlnpt2,molflaga,molflagb,extra,i2,i3,s1,s2
	character tres*2,mres1(maxstr)*2
cc	logical rmsmaton,ok

cc	double precision kabs(3,3),kabs2(3,3),a(3),e(3),b(3,3)
cc	double precision dkabs(3,3),rotat(3,3),check1(3,3),check2(3,3)
cc	double precision xcm,ycm,zcm,allmass,det,norm,allmass2
cc	double precision xcm1,ycm1,zcm1,tmpx,tmpy,tmpz
cc	double precision small
cc	real massdis

	integer level,error,ismall
	integer i,j,k,namel,l3,k3,flag
	character*7 name
cccccccccc
        real dm(4),vm(4,4),cm(3),cf(3)                                    
        real tr(3),t(3,3),q(4,4),temp(3)                                  
        real dxp(maxatm,3),dxm(maxatm,3),xf(maxatm,3),xm(maxatm,3)
        real*8 sumf(3),summ(3),xx(3)   
ccccccccccc

cc	rmsmaton=.false.
cc	molflaga=1
cc	molflagb=2
	extra=0
	c=0.24
	small=1.d20
	namel = 7
	name  = 'overlap'

24	format(100f6.3)

	j1 = s1
	j2 = s2
	ifix = 0
	if(flag.eq.1) then
	do j = 1,lres
	ifix = ifix + 1
		xf(ifix,1)=yi(j1,j,1)
		xf(ifix,2)=yi(j1,j,2)
		xf(ifix,3)=yi(j1,j,3)
		x1(ifix)=yi(j1,j,1)
		y1(ifix)=yi(j1,j,2)
		z1(ifix)=yi(j1,j,3)

c
		xm(ifix,1)=yi(j2,j,1)
		xm(ifix,2)=yi(j2,j,2)
		xm(ifix,3)=yi(j2,j,3)
		x2(ifix)=yi(j2,j,1)
		y2(ifix)=yi(j2,j,2)
		z2(ifix)=yi(j2,j,3)
	enddo
	elseif(flag.eq.2) then
	ifix = 0
	do j = 1,natom
	tres = mres1(j)
	if(tres(1:2).eq.'CA') then
	ifix = ifix + 1
		xf(ifix,1)=xi(j1,j,1)
		xf(ifix,2)=xi(j1,j,2)
		xf(ifix,3)=xi(j1,j,3)
		x1(ifix)=xi(j1,j,1)
		y1(ifix)=xi(j1,j,2)
		z1(ifix)=xi(j1,j,3)
c
		xm(ifix,1)=xi(j2,j,1)
		xm(ifix,2)=xi(j2,j,2)
		xm(ifix,3)=xi(j2,j,3)
		x2(ifix)=xi(j2,j,1)
		y2(ifix)=xi(j2,j,2)
		z2(ifix)=xi(j2,j,3)
	endif
	enddo
	elseif(flag.eq.4) then
	ifix = 0
	do j = 1,natom
	tres = mres1(j)
	if(tres(1:2).eq.'CA') then
	ifix = ifix + 1
		xf(ifix,1)=x(j1,ifix)
		xf(ifix,2)=y(j1,ifix)
		xf(ifix,3)=z(j1,ifix)
		x1(ifix)=x(j1,ifix)
		y1(ifix)=y(j1,ifix)
		z1(ifix)=z(j1,ifix)
c
		xm(ifix,1)=xi(j2,j,1)
		xm(ifix,2)=xi(j2,j,2)
		xm(ifix,3)=xi(j2,j,3)
		x2(ifix)=xi(j2,j,1)
		y2(ifix)=xi(j2,j,2)
		z2(ifix)=xi(j2,j,3)
	endif
	enddo
	elseif(flag.eq.5) then
	ifix = 0
	do j = 1,lres
	ifix = ifix + 1
		xf(ifix,1) = x(j1,j)
		xf(ifix,2) = y(j1,j)
		xf(ifix,3) = z(j1,j)
		x1(ifix) = x(j1,j)
		y1(ifix) = y(j1,j)
		z1(ifix) = z(j1,j)
c
		xm(ifix,1) = x(j2,j)
		xm(ifix,2) = y(j2,j)
		xm(ifix,3) = z(j2,j)
		x2(ifix) = x(j2,j)
		y2(ifix) = y(j2,j)
		z2(ifix) = z(j2,j)
	enddo
	endif
	imov = ifix
c--------------------------------------------------------------------
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
c --- write header to moved file                                        
cc      print *, 'r.m.s.d is : ',j1,j2,rmsd,' A'
cc	rms(j1,j2) = rmsd
cc	open (unit=7, file='brms.out',status='unknown')
cc	write(7,*) ni1
cc	do j1=1,ni1
cc	write(7,24) (rms(j1,j2),j2=1,ni2)
cc	enddo
cc	close(unit=7)

c--------------------------------------------------------------------
c************************************************
c	calculate rrmsd
	rr21 = 0.0
	dot1 = 0.0
	rg1 = 0.0
	do j3 = 1,3
	rr1(j3) = 0.0
	enddo

	do i2 = 1,lres
	rr1(1) = rr1(1) + x1(i2)
	rr1(2) = rr1(2) + y1(i2)
	rr1(3) = rr1(3) + z1(i2)
	rr21 = rr21 + x1(i2)**2 + y1(i2)**2 + z1(i2)**2 
	enddo

	dot1 = rr1(1)**2 + rr1(2)**2 + rr1(3)**2
	rg1 = (rr21 - dot1/float(lres)) /float(lres)
	rg1 = dsqrt (rg1)	

	rr22 = 0.0
	dot2 = 0.0
	rg2 = 0.0
	do j3 = 1,3
	rr1(j3) = 0.0
	enddo

	do i2 = 1,lres
	rr1(1) = rr1(1) + x2(i2)
	rr1(2) = rr1(2) + y2(i2)
	rr1(3) = rr1(3) + z2(i2)
	rr22 = rr22 + x2(i2)**2 + y2(i2)**2 + z2(i2)**2 
	enddo

	dot2 = rr1(1)**2 + rr1(2)**2 + rr1(3)**2
	rg2 = (rr22 - dot2/float(lres)) /float(lres)
	rg2 = dsqrt (rg2)
		
	rr2 = 0.0
	rr2 = ((rg1**2)+(rg2**2)) - 2.0*c*rg1*rg2
	rrms = rms/dsqrt(rr2)
c	rrms = rms / (dsqrt(((rg1**2)*(rg2**2))-(2.0*c*rg1*rg2)))

c	rmsd(j1,j2)=rms
c	rrmsd(j1,j2)=rrms

	return
	end
c********************************************************************
        subroutine inkey (answ)                                           
c ----------------------------------------------------------------------
c     reads a key as answer                                             
c ----------------------------------------------------------------------
       character answ                                                    
       read (*,'(a1)') answ                                              
       call upstrg (answ,1)                                              
       if ((answ.eq.' ').or.(answ.eq.char(13))) answ='Y'                 
       return                                                            
       end                                                               
                                                                        
       subroutine upstrg(strg,istrln)                                    
c ----------------------------------------------------------------------
c converts string str$ of lenght istrlen to upcase                      
c ----------------------------------------------------------------------
       character strg(istrln)                                            
       integer      iascii                                               
C-----change to upper case                                              
         do 3031 i=1,istrln                                             
            if (ichar(strg(i)).ge.95.and.ichar(strg(i)).le.122) then    
               iascii=ichar(strg(i))-32                                 
               strg(i)=char(iascii)                                     
            endif                                                       
3031    continue                                                       
       return                                                            
       end                                                               
                                                                        
       subroutine trpmat(n,t,tr)                                         
c --- transpose matrix -------------------------------------------------
       real t(n,n), tr(n,n)                                              
       do i=1,n                                                          
         do j=1,n                                                       
            tr(j,i)=t(i,j)                                              
         end do                                                         
       end do                                                            
       return                                                            
       end                                                               
                                                                        
       subroutine filmat(n,m,r,ifil)                                     
       real r(n,m)                                                       
c --- initialize matrix ------------------------------------------------
       do i=1,n                                                          
         do j=1,m                                                       
            r(i,j)=ifil                                                 
         end do                                                         
       end do                                                            
       return                                                            
       end                                                               
                                                                        
       subroutine rotvec (n,v,t)                                         
c --- multiply vector with matrix --------------------------------------
       real t(n,n), v(n),s(n)                                            
                                                                        
       do i=1,n                                                          
         s(i)=v(i)                                                      
         v(i)=0.0                                                       
       end do                                                            
       do i=1,n                                                          
         do j=1,n                                                       
            v(i)=v(i)+s(j)*t(i,j)                                       
         end do                                                         
       end do                                                            
       return                                                            
       end                                                               
                                                                        
       SUBROUTINE eigsrt(d,v,n,np)                                       
c ----------------------------------------------------------------------
       INTEGER n,np                                                      
       REAL d(np),v(np,np)                                               
       INTEGER i,j,k                                                     
       REAL p                                                            
       do 13 i=1,n-1                                                     
        k=i                                                             
        p=d(i)                                                          
        do 11 j=i+1,n                                                   
          if(d(j).ge.p)then                                             
            k=j                                                         
            p=d(j)                                                      
          endif                                                         
11      continue                                                        
        if(k.ne.i)then                                                  
          d(k)=d(i)                                                     
          d(i)=p                                                        
          do 12 j=1,n                                                   
            p=v(j,i)                                                    
            v(j,i)=v(j,k)                                               
            v(j,k)=p                                                    
12        continue                                                      
        endif                                                           
13     continue                                                          
       return                                                            
       END                                                               
c                                                                        
       SUBROUTINE jacobi(a,n,np,d,v,nrot)                                
c ----------------------------------------------------------------------
c     modified from numerical recipes book                              
c     one needs to set the threshold for sm from sm.eq.0 to sm.lt.10E-30
c     (anything in this range would be ok) due to underflow errors on   
c     some computers/compilers.                                         
c ----------------------------------------------------------------------
       PARAMETER (nmax=500)                                              
                                                                        
       INTEGER n,np,nrot                                                
       REAL a(np,np),d(np),v(np,np)                                      
       INTEGER i,ip,iq,j,maxrot                                          
       REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax),zero        
                                                                        
c --- zero set and iteration maximum                                 
       zero=10E-30                                                       
       maxrot=50                                                      
                                                                        
       do 12 ip=1,n                                                     
        do 11 iq=1,n                                                    
          v(ip,iq)=0.                                                 
11      continue                                                       
        v(ip,ip)=1.                                                  
12     continue                                                         
       do 13 ip=1,n                                                  
        b(ip)=a(ip,ip)                                                 
        d(ip)=b(ip)                                                   
        z(ip)=0.                                                       
13     continue                                                         
       nrot=0                                                           
       do 24 i=1,maxrot                                                 
        sm=0.                                                       
        do 15 ip=1,n-1                                                  
          do 14 iq=ip+1,n                                              
            sm=sm+abs(a(ip,iq))                                     
14        continue                                                  
15      continue                                                       
c ---   modified convergence threshold ---                            
        if(sm.lt.zero)return                                         
        if(i.lt.4)then                                              
          tresh=0.2*sm/n**2                                            
        else                                                         
          tresh=0.                                                     
        endif                                                        
        do 22 ip=1,n-1                                                 
          do 21 iq=ip+1,n                                            
            g=100.*abs(a(ip,iq))                                       
            if((i.gt.4).and.(abs(d(ip))+                               
     & g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then          
              a(ip,iq)=0.                                            
            else if(abs(a(ip,iq)).gt.tresh)then                       
              h=d(iq)-d(ip)                                             
              if(abs(h)+g.eq.abs(h))then                               
                t=a(ip,iq)/h                                            
              else                                                      
                theta=0.5*h/a(ip,iq)                                    
                t=1./(abs(theta)+sqrt(1.+theta**2))                     
                if(theta.lt.0.)t=-t                                     
              endif                                                     
              c=1./sqrt(1+t**2)                                         
              s=t*c                                                     
              tau=s/(1.+c)                                             
              h=t*a(ip,iq)                                             
              z(ip)=z(ip)-h                                             
              z(iq)=z(iq)+h                                             
              d(ip)=d(ip)-h                                             
              d(iq)=d(iq)+h                                             
              a(ip,iq)=0.                                               
              do 16 j=1,ip-1                                            
                g=a(j,ip)                                              
                h=a(j,iq)                                               
                a(j,ip)=g-s*(h+g*tau)                                   
                a(j,iq)=h+s*(g-h*tau)                                   
16            continue                                               
              do 17 j=ip+1,iq-1                                      
                g=a(ip,j)                                            
                h=a(j,iq)                                           
                a(ip,j)=g-s*(h+g*tau)                                   
                a(j,iq)=h+s*(g-h*tau)                                
17            continue                                                 
              do 18 j=iq+1,n                                         
                g=a(ip,j)                                           
                h=a(iq,j)                                          
                a(ip,j)=g-s*(h+g*tau)                             
                a(iq,j)=h+s*(g-h*tau)                             
18            continue                                              
              do 19 j=1,n                                             
                g=v(j,ip)                                            
                h=v(j,iq)                                             
                v(j,ip)=g-s*(h+g*tau)                                 
                v(j,iq)=h+s*(g-h*tau)                                 
19            continue                                              
              nrot=nrot+1                                             
            endif                                                   
21        continue                                                 
22      continue                                                    
        do 23 ip=1,n                                                
          b(ip)=b(ip)+z(ip)                                           
          d(ip)=b(ip)                                                
          z(ip)=0.                                                    
23      continue                                                      
24     continue                                                          
       pause 'too many iterations in jacobi'                           
       return                                                         
       END                                                           
c*********************************************************************
	subroutine clusterplot(ffclust,fj,tc,icnt)
	include 'mols.par'
        common /par/ natom, ntor, nhb, ns, lres
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	integer ffclust(maxstr,maxstr),fj(maxstr),tc, ic(maxstr),iz(maxstr)
	character aster*2, space*2, fp(maxstr,maxstr)*2, 
     &   ffp(maxstr,maxstr)*2,hll(maxstr)*2,vl*1,hl*2	
	data aster/' *'/
	data space/'  '/
	data vl/'|'/
	data hl/'--'/

ccc	write(31,*)'Enter the interval to draw cluster plot'
ccc	read *, in
c	read(30,*) in
c	in = ns/10
	in = icnt
	ni = ns/in
c	print *,in,ns,ni,tc

	
	do i = 1,tc
	ic(i)=0
	iz(i) = 9999
	enddo

c	do i=1,ni
c	ic(i)=0
c	end do

	do i = 1,maxstr
	hll(i) = hl
	enddo
	open (unit=22,file=of4,status='unknown')

	k = 1
1112	do i = 1,tc
	do j = 1,fj(i)
	if(ffclust(i,j).le.k*in.and.ffclust(i,j).gt.((k-1)*in).and.
     &  iz(i).ne.0) then
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

c	do i = 1,tc
c	write(31,*) i*in, ic(i)
c	enddo

c	do i = 1,tc
c	mp = ic(i)
c	write(*,'(100a2)') (aster, j = 1,mp)
c	enddo
	
c	do i = 1,tc
	do i = 1,k
	mp = ic(i)
	do j = 1,mp
	fp(i,j) = aster
	enddo
c	write(*,'(100a2)') (aster, j = 1,mp)
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

c	write(22,92) (hll(j),j=1,fj(1))
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
