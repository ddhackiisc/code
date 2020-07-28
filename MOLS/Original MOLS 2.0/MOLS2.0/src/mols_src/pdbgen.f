c     library name : pdbgen.f   

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

c	main program to generate pdb coordinates from aminio acid sequences
c	(single letter code(caps or small)-the molecule will be in the
c	extended conformation

c	program pdbgen
	subroutine pdbgen(seq,fseq)
	include 'mols.par'
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	common /par/ natom, ntor, nhb, ns, lres
	common /cord/ x1(maxatm), y1(maxatm), z1(maxatm)
	common /atm/ iatmno

	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character fseq(maxatm)*1, scode*1, ch*1, tcode*4, tres(maxatm)*4, 
     & atom(maxatm)*4, mscode(maxatm)*6, atomx(maxatm)*10, seq*(maxres),
     &fatom(maxatm)*4, fatomx(maxatm)*10, xtemp*10

	real x2(maxatm), y2(maxatm), z2(maxatm), xt, yt, zt, fx(maxatm), 
     &fy(maxatm), fz(maxatm)
	
	integer iresno2(maxatm), fatmno(maxatm), fresno(maxatm)
        CHARACTER*1 MINU(26),MAIU(26)
      DATA MINU/'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1 'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA MAIU/'A','B','C','D','E','F','G','H','I','J','K','L','M',
     1 'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/                 
10	format(a4,a1,1x,i2)
20      format(a4,4x,i3,a10,i5,4x,3f8.3)	
30      format(a4,4x,i3,a6,a4,i5,4x,3f8.3)
11	format(a18,i3)
c
c	initialize the integer variables	
c
1	iatmno=0
	iresno=0
	lres = 0
c
cc	write(*,*) 'Enter the sequence (single letter codes in caps or small)'
cc	read amino acid sequence single letter codes from terminal
cc	if the single letter code is in small letter convert it into caps
cc	read(*,'(A)') seq
c
12	format(i4)
	do i = 1,maxres
	  if(seq(i:i).ne.' ') then
	lres = lres + 1
	  do j = 1,26
      IF(seq(I:I).NE.MINU(J)) GO TO 181
      fseq(lres)=MAIU(J)
      GO TO 182
181     fseq(lres)=seq(i:i)
 	enddo
	endif
182	enddo  
	write(31,*) 'The sequence   : ',(fseq(i),i=1,lres)
	write(*,*) 'The sequence   : ',(fseq(i),i=1,lres)
	write(31,11)'sequence length: ',lres
	write(*,11)'sequence length: ',lres
	
c	if any wrong single letter codes are in the sequence
c	give caution and reread the sequece again
c
	do i=1,lres
	if(fseq(i).eq.' '.or.fseq(i).eq.'o'.or.fseq(i).eq.'O'.or.
     &fseq(i).eq.'j'.or.fseq(i).eq.'J'.or.fseq(i).eq.'u'.or.
     &fseq(i).eq.'U'.or.fseq(i).eq.'x'.or.fseq(i).eq.'X'.or.
     &fseq(i).eq.'z'.or.fseq(i).eq.'Z') then
	  write(31,*) 'Check the sequence and Enter the correct sequence'
	  write(*,*) 'Check the sequence and Enter the correct sequence'
	  goto 1
	endif
	enddo
c
c	pick up the coordinates for nterminal atoms from the lib.
c
	if(fseq(1).eq.'P') then
 	  open(unit=1,file='ALLAMINOMOLS_2.lib',status='old')
	else
	  open(unit=1,file='ALLAMINOMOLS_1.lib',status='old')
	endif
	do i=1,10
	 read(1,10)tcode,scode,natom
	 if(tcode.eq.'NH  ')then
	  do j=1,natom
	  read(1,20)atom(j),iatmno1,atomx(j),iresno1,x1(j),y1(j),z1(j)
	  if(atom(j).eq.'ATOM')then
		iatmno=iatmno+1
		iresno=iresno+1
		fx(iatmno) = x1(j)
		fy(iatmno) = y1(j)
		fz(iatmno) = z1(j)
		fatmno(iatmno) = iatmno
		fatomx(iatmno) = atomx(j)
		fatom(iatmno) = atom(j)
		fresno(iatmno) = iresno
	  endif
	 enddo
	endif
	close(unit=1)
	iresno=0
	goto 100
	enddo
c
c	pick up the coordinates for residues from lib.
c
100	do i=1,lres
	 close(unit=1)
	 xres=mod(i,2)
	 if(xres.ne.0)then
	  open(unit=1,file='ALLAMINOMOLS_1.lib',status='old')
	 else
	  open(unit=1,file='ALLAMINOMOLS_2.lib',status='old')
	 endif
	 do j=1,388
	  read(1,10)tcode,scode,natom
	  if(scode.eq.fseq(i))then
	   do k=1,natom
	    read(1,20)atom(k),iatmno1,atomx(k),iresno1,x1(k),y1(k),z1(k)
	   enddo
	   goto 200
	  endif
	 enddo
	write(31,*) 'Check your sequence'
	close(unit=1)
200	if(i.eq.1)then
	  iresno=iresno+1
	  do k = 1,natom
	    if(atom(k).eq.'ATOM')then
		iatmno=iatmno+1
		fx(iatmno) = x1(k)
		fy(iatmno) = y1(k)
		fz(iatmno) = z1(k)
		fatmno(iatmno) = iatmno
		fatomx(iatmno) = atomx(k)
		fatom(iatmno) = atom(k)
		fresno(iatmno) = iresno
 	    endif
	    xt=x1(natom)
	    yt=y1(natom)
	    zt=z1(natom)
	  enddo
	else
	  iresno=iresno+1
	  do k = 1,natom
	    x1(k)=x1(k)+xt
	    y1(k)=y1(k)+yt
	    z1(k)=z1(k)+zt
	    if(atom(k).eq.'ATOM')then
		iatmno=iatmno+1
		fx(iatmno) = x1(k)
		fy(iatmno) = y1(k)
		fz(iatmno) = z1(k)
		fatmno(iatmno) = iatmno
		fatomx(iatmno) = atomx(k)
		fatom(iatmno) = atom(k)
		fresno(iatmno) = iresno
	    endif
	 enddo
	 xt=x1(natom)
	 yt=y1(natom)
	 zt=z1(natom)
	endif
400	enddo
c
c	pick up the coordinates for c terminal atoms from lib.
c
	close(unit=1)
	xres=mod(lres,2)
	if(xres.ne.0)then
	  open(unit=1,file='ALLAMINOMOLS_1.lib',status='old')
	else
	  open(unit=1,file='ALLAMINOMOLS_2.lib',status='old')
	endif
	do j=1,388
	  read(1,10)tcode,scode,natom
	  if(tcode.eq.'COH ')then
	   do k=1,natom
	    read(1,20)atom(k),iatmno1,atomx(k),iresno1,x1(k),y1(k),z1(k)
	    x1(k)=x1(k)+xt
	    y1(k)=y1(k)+yt
	    z1(k)=z1(k)+zt
	    if(atom(k).eq.'ATOM')then
		iatmno=iatmno+1
		fx(iatmno) = x1(k)
		fy(iatmno) = y1(k)
		fz(iatmno) = z1(k)
		fatmno(iatmno) = iatmno
		fatomx(iatmno) = atomx(k)
		fatom(iatmno) = atom(k)
		fresno(iatmno) = iresno
	    endif
	   enddo
	   close(unit=1)
	 endif
	enddo
c
c	if proline residue is in the sequence, rotate it to fix its
c	phi as constant
c
	do k=1,iatmno
	   x1(k) = fx(k)
	   y1(k) = fy(k)
	   z1(k) = fz(k)
	   xtemp = fatomx(k)
	   atom(k) = fatom(k)
	   tres(k) = xtemp(7:10)
	   iatomno1 = fatmno(k)
	   iresno2(k) = fresno(k)
	   mscode(k) = xtemp(1:6)
	enddo
	k=1
500	if(tres(k).eq.'PRO ')then
	  call molgen1(k)
	  k=k+14
	else
	  k=k+1
	endif
	if(k.lt.iatmno) goto 500
	iatmno1=1
	tres(1)=tres(2)
c
c	write the final co-ordinates of molecule
c
	open(unit=4,file=pf1,status='unknown')
	print *,'iatmno :  ',iatmno
	do k=1,iatmno
	  write(4,30)atom(k),iatmno1,mscode(k),tres(k),iresno2(k),x1(k),
     &    y1(k),z1(k)
	  iatmno1=iatmno1+1
	enddo
	close(unit=4)
	write(31,11)'No. of atoms   : ',iatmno1-1
	write(*,11)'No. of atoms   : ',iatmno1-1
	natom = iatmno1-1

c	stop
	return
	end

**********************************************************
        subroutine molgen1(nn)
	include 'mols.par'
	common /par/ natom, ntor, nhb, ns, lres
	common /cord/ x1(maxatm), y1(maxatm), z1(maxatm)
	common/atm/iatmno
        dimension x_one(maxatm,3), x_two(maxatm,3)
        do k=1,iatmno
             x_one(k,1)=x1(k) 
             x_one(k,2)=y1(k) 
             x_one(k,3)=z1(k) 
        enddo
        phi=109.0
        do if=1,1
 
C###################### PHI ALL ####################################      

        call elemen(x_one(nn,1),x_one(nn,2),
     &              x_one(nn,3),
     &              x_one(nn+1,1),x_one(nn+1,2),
     &              x_one(nn+1,3),
     &              el,em,en)

        do k=1,nn+10
           do ki=1,3
             x_two(k,ki)=x_one(k,ki) 
           enddo
        enddo
        do k=nn+11,iatmno
           xin=x_one(k,1)-x_one(nn,1)
           yin=x_one(k,2)-x_one(nn,2)
           zin=x_one(k,3)-x_one(nn,3)
           call rotor(el,em,en,phi,xin,yin,zin,xout,yout,zout)
           x_two(k,1)=xout+x_one(nn,1)
           x_two(k,2)=yout+x_one(nn,2)
           x_two(k,3)=zout+x_one(nn,3)
        enddo

c        do k=nn+1,iatmno
c          do ki=1,3
c             x_two(k,ki)=x_one(k,ki) 
c          enddo
c        enddo

        do k=1,iatmno
           do ki=1,3
              x_one(k,ki)=x_two(k,ki) 
           enddo
        enddo
C####################################################################      
        enddo

         do k=1,iatmno
             x1(k)=x_two(k,1) 
             y1(k)=x_two(k,2) 
             z1(k)=x_two(k,3) 
        enddo
        return
        end
c***********************************************************************
	subroutine elemen(x1,y1,z1,x2,y2,z2,el,em,en)

c*****SUBROUTINE TO CALCULATE DIRECTION COSINES*****
	d = dist(x1,y1,z1,x2,y2,z2)
	   el=(x2-x1)/d
	   em=(y2-y1)/d
	   en=(z2-z1)/d
	return
	end

c***********************************************************************
	subroutine rotor(xr1,yr1,zr1,theta1,xr3,yr3,zr3,xr4,yr4,zr4)

c******SUBROUTINE TO ROTATE A POINT ABOUT A GENERAL AXIS BY THETA*****

            xp=xr3
	    yp=yr3
	    zp=zr3
	theta=theta1*3.1415927/180.0
c********* ROTATION MATRIX *********************************************
                    em11=cos(theta)+xr1*xr1*(1.0-cos(theta))
                    em12=xr1*yr1*(1.0-cos(theta))-zr1*sin(theta)
                    em13=xr1*zr1*(1.0-cos(theta))+yr1*sin(theta)
                    em21=xr1*yr1*(1.0-cos(theta))+zr1*sin(theta)
		    em22=cos(theta)+yr1*yr1*(1.0-cos(theta))
		    em23=yr1*zr1*(1.0-cos(theta))-xr1*sin(theta)
		    em31=xr1*zr1*(1.0-cos(theta))-yr1*sin(theta)
		    em32=yr1*zr1*(1.0-cos(theta))+xr1*sin(theta)
		    em33=cos(theta)+zr1*zr1*(1.0-cos(theta))
c************************************************************************
		xt=em11*xp+em12*yp+em13*zp
	        yt=em21*xp+em22*yp+em23*zp
		zt=em31*xp+em32*yp+em33*zp
	    xr4=xt
	    yr4=yt
	    zr4=zt
c	write(31,*)xp,yp,zp,' ',xr4,yr4,zr4
c	write(*,*)xp,yp,zp,' ',xr4,yr4,zr4
	return
	end

c*********************************************************************
	function dist(x1,y1,z1,x2,y2,z2)
	  dist = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
	return
	end
c*********************************************************************
