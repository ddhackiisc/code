	program prelims
!	program prelims
!	program to do all the initial parameter calculations
!	like centroid calcluation, cutoff distance caln
!	and plp atpmtype assignment for ligand.
	include 'mols.par'
        common /prel/bx,by,bz
        common /cuto/big
        common /att/ hbtyp(maxi)
	common /atn/natom

	character*80 str*80,hbtyp*3
	call centroid
c	call cutoff
c	call plp_at
c	print*,'bx=',bx,'by=',by,'bz=',bz,big,natom
c	do m=1,natom
c	write(*,*)m,hbtyp(m)
c	enddo
	end
!#######################################################################

	subroutine centroid
!c	program to read the pdb file and extract the atom records of ligand
!	and calc centroid.
	include 'mols.par'
!	parameter (maxi=3000)
	common /prel/bx,by,bz
	character str*80,st1*54
	real xsum,ysum,zsum
	real bx,by,bz,rden
	dimension x(mnatp,3)
	
10	format(a80)
20	format(a54)
	open(unit=1,file='pocket.pdb',status='old')
	open(unit=2,file='part.pdb',status='unknown')

	do i=1,mnatp
	read(1,10,end=99),str
		if(str(1:6) .eq.'ATOM') then
		write(2,10),str
		if(str(18:20) .eq.'NIP')then
c		write(2,10),str
		endif
		endif
c		if(str(1:3) .eq.'END') then
c		goto 100
c		endif
	enddo
99	write(*,*),'line number',i
c100	pause
	close(unit=1)
	close(unit=2)

30	format(a30,3f8.3)
100	open(unit=2,file='part.pdb',status='old')
	do i=1,mnatp
	read(2,30,end=166)str,(x(i,j),j=1,3)
c	write(*,30)str,(x(i,j),j=1,3)
	enddo
166	 write(*,*),'line number',i

	k=i-1
	rden=k
	xsum=0
	ysum=0
	zsum=0

	do j=1,k
	xsum=(xsum)+(x(j,1))
	ysum=(ysum)+(x(j,2))
	zsum=(zsum)+(x(j,3))
c	print*,j
	enddo

	bx=xsum/rden
	by=ysum/rden
	bz=zsum/rden
40	format(3f8.3)
	print*,'xsum=',xsum,'ysum=',ysum,'zsum=',zsum
	print*,'bx=',bx,'by=',by,'bz=',bz

c	write(6,40)xsum,ysum,zsum
c	write(6,40)bx,by,bz

	CLOSE(UNIT=2)
	return

	end
!***********************************************************************
	subroutine cutoff
c	program to calculate the cut-off distance for grid box
	include 'mols.par'
!	parameter(maxi=1500)
	common /cuto/big
	character sdf_str*80
	integer atmno
	real rx(maxi,3)
	
10	format(a80)
20	format(i3,a50)
30	format(3f10.5,1x,a2,1x,a30)
	open(unit=1,file='rough.sdf',status='old') 

	do i=1,3
	read(1,10)sdf_str
c	write(*,10)sdf_str
	enddo

	read(1,20)atmno,sdf_str
c	write(*,20)atmno,sdf_str

	do i=1,atmno
	read(1,30)(rx(i,j),j=1,3),atomsym,sdf_str
c	write(*,30)(rx(i,j),j=1,3),atomsym,sdf_str
	enddo
	close(unit=1)

	big=0
	do j=1,atmno
	  do k=j+1,atmno
	   d=dist(rx(k,1),rx(k,2),rx(k,3),rx(j,1),rx(j,2),rx(j,3))
c	  print*,j,k,d
	    if(d .ge. big) then
	     big=d
	    endif
	 enddo
	enddo
c	print*,'big=',big
	return
	end


c*********************************************************************
        function dist(x1,y1,z1,x2,y2,z2)
          dist = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        return
        end
c*********************************************************************

!*********************************************************************
	subroutine plp_at
c	program to read atomtype from the file atomtype.dat
	include 'mols.par'
!	parameter (maxi=1500)
        common /string/ str1(maxi),str2(maxi),str3(maxi),res,ato(maxi)
        common /att/ hbtyp(maxi)
	common /atn/natom
  
        character*80 str1*30,str2*30,str3*30,res*10,ato*80,hbtyp*3
        integer natom
10      format(a80)
11      format(a25)
20      format(i3,i3,a8)
30      format(a30)

        open(unit=1,file='out.mol2',status='old')
        open(unit=2,file='atomtype.dat',status='old')

	do l=1,31
	read(2,30),str3(l)
c	write(*,30),str3(l)
	enddo

        do i=1,2
        read(1,11)str1(i)
c        write(*,11)str1(i)
        enddo

        read(1,20)natom,ibno,res
c        write(*,20)natom,ibno,res

        do j=1,5
        read(1,11)str2(j)
c        write(*,11)str2(j)
        enddo

        do  k=1,natom
	read(1,10),ato(k)
c	write(*,10),ato(k)
	do  n=1,31
	 if( ato(k)(48:55) .eq. str3(n)(6:13) ) then
	hbtyp(k)= str3(n)(23:24)
c	write(6,40),k,ato(k)(48:55),str3(n)(6:13),str3(n)(23:24)
	 endif
40	format(1x,i3,a10,1x,a10,1x,a3)
	enddo
        enddo
c	do m=1,natom
c	lat(m)=attyp(m)
c	write(*,*)m,attyp(m)
c	enddo

	close(1)
	close(2)
	end
!***********************************************************************




