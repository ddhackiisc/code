c-----------------------------------------------------------------------------
c     library name : smpdbgen.f   

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

	subroutine pdbgen
c	program pdbgen
c	Program to call the smi23d program
        include 'mols.par'

        common /ligcrdb/ylig(maxatm,8)
	common /par/natom,ntor,nhb,ns,lres
	common /gen/ibno,iat1(maxatm),iat2(maxatm),isno(maxatm)
        common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &   iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),ind(maxatm,maxatm),
     &   big
	common /left/le(maxatm,maxatm),innd(maxatm),ielenum(maxatm),
     &	bonum(maxatm)
        common /string/atsym(maxatm),desc(maxatm),attyp(maxatm),
     &  str1(maxatm),str2(maxatm),str3,res,n1,n2,n3
    
	character*50 desc*30,attyp*10,atsym*6,str1*30,str2*30,
     &  str3*30,res*10
	integer el

	call prep_coord
c       write(31,*)'No.of rotatable bonds:',ntor
c       write(31,*)'No. of atoms         :',natom

	end

c	*******************************************************************
	subroutine prep_coord
	include 'mols.par'
        common /native/nat,pepfile
	character pdb*80,str*80,str1*80,str2*80,rpdb*80,pdb1*80


        if(nat.eq.0)then
	pdb ='sh script1'  ! this script calls the routines for 3d 
	                   ! coordinates and calcluates initial conformation 
                           ! energy
        call system(pdb)
c-------------------------
        call lchange
c-------------------------
        pdb1='sh script01'
	call system(pdb1)
	endif
c-------------------------
        if(nat.eq.1)then
        pdb ='sh script10'        
        call system(pdb)
        call lchange
        endif
c-------------------------
	call coorgen
c-------------------------
        if(nat.eq.0)then
	call store_index
        endif
c-------------------------
	end
c	--------------------------------------------------------------------
        subroutine lchange
        include 'mols.par'

        common /native/nat,pepfile

        character str*80,str3*80
        character desc*30,attyp*10,atsym*6,str1(maxatm)*30,
     &  str2(maxatm)*30
        character  res*10,str4(maxatm)*80,str5(maxatm)*80,
     &  str41(maxatm)*80,str42(maxatm)*80
        integer atno,ibno,natom,n1,n2,n3
        real r,big

10      format(a80)
11      format(a25)
21      format(i4,i4,i2,i2,i2)
50      format(a7,1x,a4,a64)
60      format(a7,1x,a2,2x,a64)

        open(unit=1,file='lout.mol2',status='old')
        open(unit=2,file='out.mol2',status='unknown')
        do i=1,2
        read(1,11)str1(i)
        write(2,11)str1(i)
        enddo

        read(1,*)atno,ibno,n1,n2,n3
        write(2,21)atno,ibno,n1,n2,n3

        natom=atno
        
        if(nat.eq.0)then
        do i=1,4
        read(1,11)str2(i)
        write(2,11)str2(i)
        enddo
        else if(nat.eq.1)then
        do i=1,5
        read(1,11)str2(i)
        write(2,11)str2(i)
        enddo        
        endif
        
        do i=1,atno

        read(1,50),str4(i),str41(i),str42(i)
        write(2,60),str4(i),str41(i),str42(i)
        enddo

        read(1,11)str3
        write(2,11)str3

        do j=1,ibno
        read(1,10)str5(j)
        write(2,10)str5(j)
        enddo

        return
        end

c       --------------------------------------------------------------------       
	subroutine store_index
        include 'mols.par'
c       program to store the indices of the rotatable bonds and its neighbours

       	common /ligcrdb/ylig(maxatm,8)
	common /par/natom,ntor,nhb,ns,lres
	common /left/le(maxatm,maxatm),innd(maxatm),ielenum(maxatm),
     &  bonum(maxatm)
        common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &   iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),
     &   ind(maxatm,maxatm),big

        character str*80
        integer el,tem(35)
c       integer numrot,bonum(50)   !declaration for remaining atoms list

10      format(a9,1x,i3)
20      format(a16,1x,i3)
30      format(i3,5x,i3)
40      format(i4)
50	format(a19,1x,i4)
60	format(a8,1x,i4)
70	format(a13,1x,i4)

        open(unit=4,file='out.out',status='old')
        open(unit=2,file='le.out',status='unknown')
        read(4,10)str,inrot !str=NROTBONDS,inrot=number of rotatable bonds
	 ntor=inrot !ntor = number of torsions = number of rotatable bonds
        do i=1,inrot
          read(4,30)(irotb(i,j),j=1,2)!irotb(i,1)=bonding atom no,irotb(i,2)=bond atom no
          read(4,20)str,el !str = Num_array_elment,el=0 
          iatrot(i) = el
          tem(i)=natom-iatrot(i)
              do k=1,iatrot(i)
                 read(4,40)ilsrot(i,k)
	         ind(i,ilsrot(i,k))=1
                enddo
        enddo
	write(2,50)'number of rot_bonds',inrot
	do i=1,inrot
	write(2,60)'Rot_Bond',i
	write(2,70)'No.of Element',tem(i)
	 do k=1,natom
	    if (ind(i,k).ne.1) then
	      write(2,40),k
	    endif
          enddo
        enddo
	close(unit=2)

c	This is to list the atoms that are not displaced from its position

	open(unit=3,file='le.out',status='old')

	read(3,50)str,inrot
	do i=1,inrot
	read(3,60)str,bonum(i)
	read(3,70)str,ielenum(i)
		do j=1,ielenum(i)
		read(3,40)le(i,j)
		enddo
	enddo
	close(unit=3)

	return
        end

c	*******************************************************************
	subroutine coorgen
	include 'mols.par'
c	program to extract the coordinates from the generated mol2. 

	common /comment/icomment
        common /native/nat,pepfile
        common /ligcrdb/ylig(maxatm,8)
        common /par/natom, ntor, nhb, ns, lres
	common /gen/ibno,iat1(maxatm),iat2(maxatm),isno(maxatm)
     	common /left/le(maxatm,maxatm),innd(maxatm),ielenum(maxatm),
     &	bonum(maxatm)
	common /string/atsym(maxatm),desc(maxatm),attyp(maxatm),
     &	str1(maxatm),str2(maxatm),str3,res,n1,n2,n3
        common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &   iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),
     &   ind(maxatm,maxatm),big

	character str*80
	character*50 desc*30,attyp*10,atsym*6,str1*30,str2*30,
     &	str3*30,res*10
	integer atno
	real r,big,d

11	format(a25)
21      format(i4,i4,i2,i2,i2)
30      format(3x,i4,1x,a6,2x,3f10.4,1x,a29)
50	format(1x,i5,1x,i5,1x,i5,a10)

	open(unit=1,file='out.mol2',status='old')
	
	do i=1,2
	read(1,11)str1(i)
	enddo

	read(1,21)atno,ibno,n1,n2,n3

	natom=atno
        if(nat.eq.0)then
	 do i=1,4
	 read(1,11)str2(i)
	 enddo
        else if(nat.eq.1)then
         do i=1,5
         read(1,11)str2(i)
         enddo
        endif

	do i=1,atno
        read(1,30),innd(i),atsym(i),(rx(i,j),j=1,3),desc(i)
	ylig(i,j)=rx(i,j)
	enddo

	read(1,11)str3

        do j=1,ibno
        read(1,50)isno(j),iat1(j),iat2(j),attyp(j)
        enddo
	big = 0.0
        do k=1,natom
          do l=k+1,natom
           d=(dist(rx(l,1),rx(l,2),rx(l,3),rx(k,1),rx(k,2),rx(k,3)))/2
            if(d.ge.big) then
             big=d
            endif
         enddo
        enddo
        if(icomment.eq.1)then
        print *,'big=',big !this is used to find the length of the ligand
        endif
	close(unit=1)
	return
	end
c*************************************************************************

