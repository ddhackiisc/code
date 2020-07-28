c	Peptide-Protein Docking
c	Generate MOLS parameter file for Peptide from inp.pdb
	subroutine pamppar(iopt)
	include 'mols.par'
	parameter (maxbnd = 819, jx = 25)

	common /comment/icomment
	common /hbpar/mnhb
	common /par/natom,ntor,nhb,ns,lres
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        
	character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
	character atname(maxatm)*4, rename(maxatm)*4, tcode*4, scode*1,
     &  an(maxatm,jx)*4, bn(maxatm,jx)*4, a(maxatm,jx)*4, 
     &  b(maxatm,jx)*4, c(maxatm,jx)*4, cn(maxatm,jx)*4, 
     &  dn(maxatm,jx)*4, m14(maxatm,jx)*4,atype*2, atnam*2, at1*4, 
     &  at2*1, at3*2, d1*4, d2*4, d3*4, d4*4
	integer resno(maxatm), atmno(maxatm), na(maxatm,jx), ka(maxatm),
     &  kc(maxatm),ma(maxatm,jx), kb(maxatm), mb(maxatm,jx), 
     &  nb(maxatm,jx), mc(maxatm,jx), nc(maxatm,jx), kd(maxatm), 
     &  nd(maxatm,maxatm), ke(maxatm), ne(maxatm,jx),k14(maxatm),
     &  k15(maxatm), n15(maxatm,jx), md(maxatm,jx), n14(maxatm,jx)

10	format(8x,i3,1x,a4,1x,a4,3x,i3)
20	format(a4,a1,1x,i2)
30	format(8x,2a4)
40	format(10i4)
60	format(2x,20i4)
70	format(2x,a2,6x,f6.3,5x,f6.3)
80	format(a4,7x,a1,2x,a2,f8.4)
90	format(a4,a1,1x,i2)
100	format(i4,1x,f8.4,f7.3,f7.3,f4.1,f4.1)
110	format(8x,4a4)
120	format(i4,i4,i6,i6)
130	format(5x,a27,2x,i4)
140	format(5x,a27,2x,i4)

	open(unit=1,file=pf1,status='old')
	do i=1,natom
		read(1,10) atmno(i),atname(i),rename(i),resno(i)
	enddo
	close(unit=1)
c****************************to find first neighbours*************************
      do i=1,natom
	ka(i)=1
	open(unit=2,file='PALLCONN.lib',status='old')
	rewind 2
	do i1=1,maxbnd
	  read(2,20)tcode,scode,iatom
	  if(tcode.eq.rename(i))then
	   do l=1,iatom
	    read(2,30) (a(i,j),j=1,2)
	    if(a(i,1).eq.atname(i))then
	      do j=1,natom
		if(atname(j).eq.a(i,2).and.a(i,2).ne.'    ')then
	 	  if(a(i,1).eq.' C  '.and.a(i,2).eq.' N  ')then
		     itemp=resno(i)+1
		  else
	 	     itemp=resno(i)
		  endif
	          if(resno(j).eq.itemp)then
			ma(i,1)=i
			ma(i,2)=j
		  endif
	        endif
	      enddo
	      if(ma(i,2).ne.na(i,1).and.ma(i,2).ne.na(i,2).and.
     &	      ma(i,2).ne.na(i,3).and.ma(i,2).ne.na(i,4))then
		na(i,ka(i))=ma(i,2)
		an(i,ka(i))=a(i,2)
		ka(i)=ka(i)+1
	      endif
	    endif
	  enddo
	 endif
	enddo
c*************************to find second neighbours**************************
	kb(i)=1
	do l2=1,ka(i)-1
	rewind 2
	  do i2=1,maxbnd
	    read(2,20)tcode,scode,iatom
	    if(tcode.eq.rename(i))then
	      do l1=1,iatom
		read(2,30) (b(i,j),j=1,2)
		if(an(i,l2).eq.b(i,1))then
		  do j=1,natom
	 	    if(atname(j).eq.b(i,2).and.b(i,2).ne.'    ')then
		      if(b(i,1).eq.' C  '.and.b(i,2).eq.' N  ')then
			itemp=resno(i)+1
		      else
			itemp=resno(i)
		      endif
		        if(b(i,1).eq.' N ') itemp = resno(na(i,l2))
			if(resno(j).eq.itemp)then
			  mb(i,1)=i
			  mb(i,2)=j
			endif
		    endif
		  enddo
		  if(mb(i,2).ne.i.and.mb(i,2).ne.nb(i,1).and.
     &		    mb(i,2).ne.nb(i,2).and.mb(i,2).ne.nb(i,3).and.
     &		    mb(i,2).ne.nb(i,4).and.mb(i,2).ne.nb(i,5).and.
     &		    mb(i,2).ne.nb(i,6))then
			nb(i,kb(i))=mb(i,2)
			bn(i,kb(i))=b(i,2)
			kb(i)=kb(i)+1
		  endif
	       endif
	     enddo
	   endif
	 enddo
	enddo
C************************to find third neighbours***************************
	kc(i)=1
	do l3=1,kb(i)-1
	rewind 2
	    do i3=1,maxbnd
	      read(2,20)tcode,scode,iatom
		if(tcode.eq.rename(nb(i,l3)))then
	  	  do l2=1,iatom
		    read(2,30) (c(i,j),j=1,2)
		    if(bn(i,l3).eq.c(i,1))then
		      do j=1,natom
			if(atname(j).eq.c(i,2).and.c(i,2).ne.'    ')then
			  if(c(i,1).eq.' C  '.and.c(i,2).eq.' N  ')then
				itemp=resno(i)+1
			  else
				itemp=resno(i)
			  endif
		          if(c(i,1).eq.' N  ') itemp = resno(nb(i,l3))
		          if(c(i,1).eq.' CA ') itemp = resno(nb(i,l3))
			  if(resno(j).eq.itemp)then
				mc(i,1)=i
				mc(i,2)=j
			  endif
			endif
		      enddo
	if(mc(i,2).ne.i.and.mc(i,2).ne.na(i,1).and.mc(i,2).ne.na(i,2)
     & .and.mc(i,2).ne.na(i,3).and.mc(i,2).ne.na(i,4).and.
     & mc(i,2).ne.na(i,5).and.mc(i,2).gt.i.and.mc(i,2).ne.nc(i,1).
     & and.mc(i,2).ne.nc(i,2).and.mc(i,2).ne.nc(i,3).and.mc(i,2).ne.
     & nc(i,4).and.mc(i,2).ne.nc(i,5))then
		nc(i,kc(i))=mc(i,2)
		cn(i,kc(i))=c(i,2)
		kc(i)=kc(i)+1
	endif
		endif
	      enddo
	    endif
	  enddo
	enddo
c****************************************************
      enddo
c
c*********************to find fourth and above neighbours**********************
	do i = 1,natom
	  kd(i)=1
	  do j = i+1,natom
	    do j1=1,ka(i)-1
	      if(na(i,j1).eq.j) goto 111
	    enddo
	    do j2=1,kb(i)-1
	      if(nb(i,j2).eq.j) goto 111
	    enddo
	    do j3=1,kc(i)-1
	      if(nc(i,j3).eq.j) goto 111
	    enddo
 	    nd(i,kd(i))=j
	    kd(i)=kd(i)+1
111	  enddo
	enddo
c*********************to count the atom pairs fourth & above neighbours********
	do i=1,natom
	  ke(i)=1
	  ne(i,ke(i))=nd(i,1)
	  do j=1,kd(i)-1
	    if(nd(i,j).ne.(nd(i,j-1)+1).and.j.gt.1)then
	      ke(i)=ke(i)+1
	      ne(i,ke(i))=nd(i,j-1)
	      ke(i)=ke(i)+1
	      ne(i,ke(i))=nd(i,j)
	    endif
	  enddo
	  ke(i)=ke(i)+1
	  ne(i,ke(i))=nd(i,kd(i)-1)
	  if(ne(i,1).eq.0) ke(i)=0
	enddo
c
c***************************************************************
c	to exclude the atompairs in around peptide bond
c	and in the rings(1-5 pairs)
C**************************************************************
	do i=1,natom
	k15(i)=0
	do j=1,ke(i)
	j1=ne(i,j)
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE1'.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE2'.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HH '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' OH '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' HH '.and.resno(i).eq.resno(j1)) goto 131

	if(rename(i).eq.'TYR '.and.atname(i).eq.' CD1'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TYR '.and.atname(i).eq.' HD1'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CD2'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' HD2'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CE1'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' HE1'.and.
     &resno(i).eq.resno(j1))  goto 151
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CE2'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' HE2'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CZ '.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' OH '.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' HH '.and.
     &resno(i).eq.resno(j1)) goto 151 

	if(rename(i).eq.'PHE '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE1'.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE2'.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HZ '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' HZ '.and.resno(i).eq.resno(j1)) goto 131

	if(rename(i).eq.'PHE '.and.atname(i).eq.' CD1'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'PHE '.and.atname(i).eq.' HD1'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CD2'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' HD2'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CE1'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' HE1'.and.
     &resno(i).eq.resno(j1))  goto 151
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CE2'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' HE2'.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CZ '.and.
     &resno(i).eq.resno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' HZ '.and.
     &resno(i).eq.resno(j1)) goto 151 

	if(rename(i).eq.'TRP '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE1'.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE3'.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HH2'.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' HZ2'.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' HZ3'.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' HH2'.and.resno(i).eq.resno(j1)) goto 131

	if(rename(i).eq.'TRP '.and.atname(i).eq.' CD1'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HD1'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CD2'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' NE1'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HE1'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CE2'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CE3'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HE3'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CZ2'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HZ2'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CZ3'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HZ3'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CH2'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HH2'.and.
     &resno(i).eq.resno(j1)) goto 151

	if(rename(i).eq.'HIS '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE1'.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'HIS '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE2'.and.resno(i).eq.resno(j1)) goto 131

	if(rename(i).eq.'HIS '.and.atname(i).eq.' HD2'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' ND1'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' CD2'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' HD2'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' CE1'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' HE1'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' NE2'.and.
     &resno(i).eq.resno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' HE2'.and.
     &resno(i).eq.resno(j1)) goto 151

	if(rename(j1).eq.'PRO '.and.atname(i).eq.' C  '
     &.and.(resno(i)+1).eq.resno(j1)) goto 151
	if(rename(j1).eq.'PRO '.and.atname(i).eq.' O  '
     &.and.(resno(i)+1).eq.resno(j1)) goto 151
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '
     &.and.atname(j1).eq.'1HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '
     &.and.atname(j1).eq.'2HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CA '
     &.and.atname(j1).eq.'1HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CA '
     &.and.atname(j1).eq.'2HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.'1HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.'2HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'1HB '
     &.and.atname(j1).eq.'2HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'1HB '
     &.and.atname(j1).eq.'1HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'1HB '
     &.and.atname(j1).eq.'2HD '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'1HB '
     &.and.atname(j1).eq.'1HD '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'2HB '
     &.and.atname(j1).eq.'2HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'2HB '
     &.and.atname(j1).eq.'1HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'2HB '
     &.and.atname(j1).eq.'2HD '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'2HB '
     &.and.atname(j1).eq.'1HD '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.'2HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.'1HG '.and.resno(i).eq.resno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'1HG '
     &.and.atname(j1).eq.'2HG '.and.resno(i).eq.resno(j1)) goto 131
c
	goto 161
c
151	do j2=i+1,natom
	  if(atname(j2).eq.' C  ') then
		k15(i)=k15(i)+1
		n15(i,k15(i))=j2
		k15(i)=k15(i)+1
		n15(i,k15(i))=natom
		goto 141
	  endif
	enddo
c
161	k15(i)=k15(i)+1
	n15(i,k15(i))=ne(i,j)
131	enddo
141	enddo
	
c******************************************************************
c	to exclude the atompairs in around peptide bond
c	and in the rings(in 1-4 pairs) and pickup needed 1-4 pairs
C******************************************************************
	do i = 1, natom
	k14(i)=1
	do l3=1,kb(i)-1
	rewind 2
	do i3=1,maxbnd
	read(2,20)tcode,scode,iatom
	if(tcode.eq.rename(nb(i,l3)))then
	do l2=1,iatom
	read(2,30) (dn(i,j),j=1,2)
	if(bn(i,l3).eq.dn(i,1))then
	do j=1,natom
	if(atname(i).eq.' CA '.and.atname(j).eq.' CA '.and.
     &(resno(i)+1).eq.resno(j)) goto 121
	if(atname(i).eq.' CA '.and.atname(j).eq.' H  '.and.
     &(resno(i)+1).eq.resno(j)) goto 121
	if(atname(i).eq.' O  '.and.atname(j).eq.' CA '.and.
     &(resno(i)+1).eq.resno(j)) goto 121
	if(atname(i).eq.' O  '.and.atname(j).eq.' H  '.and.
     &(resno(i)+1).eq.resno(j)) goto 121

	if(atname(i).eq.' CB '.and.atname(j).eq.' HD1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' HD2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CZ '.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CE1'.and.atname(j).eq.' HE2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CE1'.and.atname(j).eq.' HH '.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CE2'.and.atname(j).eq.' HH '.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CD1'.and.rename(i).eq.'TYR '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HD1'.and.rename(i).eq.'TYR '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CD2'.and.rename(i).eq.'TYR '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HD2'.and.rename(i).eq.'TYR '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HE1'.and.rename(i).eq.'TYR '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HE2'.and.rename(i).eq.'TYR '.and.
     &resno(i).eq.resno(j)) goto 121

	if(atname(i).eq.' CB '.and.atname(j).eq.' HD1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' HD2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CZ '.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'PHE ') goto 121

	if(atname(i).eq.' CD1'.and.rename(i).eq.'PHE '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HD1'.and.rename(i).eq.'PHE '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CD2'.and.rename(i).eq.'PHE '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HD2'.and.rename(i).eq.'PHE '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CE1'.and.rename(i).eq.'PHE '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CE2'.and.rename(i).eq.'PHE '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HE1'.and.rename(i).eq.'PHE '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HE2'.and.rename(i).eq.'PHE '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CZ '.and.rename(i).eq.'PHE '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HZ '.and.rename(i).eq.'PHE '.and.
     &resno(i).eq.resno(j)) goto 121

	if(atname(i).eq.' CB '.and.atname(j).eq.' HD1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' NE1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE3'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CE2 '.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CZ2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE3'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' NE1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CZ3'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HZ3'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'TRP ') goto 121

	if(atname(i).eq.' CD1'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HD1'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CD2'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' NE1'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HE1'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CE2'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CE3'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HE3'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CZ2'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HZ2'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CZ3'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HZ3'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CH2'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HH2'.and.rename(i).eq.'TRP '.and.
     &resno(i).eq.resno(j)) goto 121

	if(atname(i).eq.' CB '.and.atname(j).eq.' CE1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' NE2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' HD2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' NE2 '.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE2'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CE1'.and.
     &resno(i).eq.resno(j).and.rename(i).eq.'HIS ') goto 121

	if(atname(i).eq.' ND1'.and.rename(i).eq.'HIS '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CE1'.and.rename(i).eq.'HIS '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HE1'.and.rename(i).eq.'HIS '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' NE2'.and.rename(i).eq.'HIS '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HE2'.and.rename(i).eq.'HIS '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CD2'.and.rename(i).eq.'HIS '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' HD2'.and.rename(i).eq.'HIS '.and.
     &resno(i).eq.resno(j)) goto 121

	if(atname(i).eq.' C  '.and.rename(j).eq.'PRO ') goto 121

	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.'1HB '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.'2HB '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.' CG '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.' CB '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.'1HG '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.'2HG '.and.resno(i).eq.resno(j)) goto 121

	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.'1HB '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.'2HB '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.' CG '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.' CB '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.'1HG '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.'2HG '.and.resno(i).eq.resno(j)) goto 121

	if(rename(i).eq.'PRO '.and.atname(i).eq.' CB '.and.
     &atname(j).eq.'1HD '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CB '.and.
     &atname(j).eq.'2HD '.and.resno(i).eq.resno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CB '.and.
     &atname(j).eq.' CD '.and.resno(i).eq.resno(j)) goto 121

	if(atname(i).eq.' CA '.and.rename(i).eq.'PRO '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.'1HB '.and.rename(i).eq.'PRO '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.'2HB '.and.rename(i).eq.'PRO '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CG '.and.rename(i).eq.'PRO '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.'1HG '.and.rename(i).eq.'PRO '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.'2HG '.and.rename(i).eq.'PRO '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.' CD '.and.rename(i).eq.'PRO '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.'1HD '.and.rename(i).eq.'PRO '.and.
     &resno(i).eq.resno(j)) goto 121
	if(atname(i).eq.'2HD '.and.rename(i).eq.'PRO '.and.
     &resno(i).eq.resno(j)) goto 121
c
	if(atname(j).eq.dn(i,2).and.dn(i,2).ne.'    ')then
	  if(dn(i,1).eq.' C  '.and.dn(i,2).eq.' N  ')then
	    itemp=resno(i)+1
	  else
	    itemp=resno(i)
	  endif
          if(dn(i,1).eq.' N  ') itemp = resno(nb(i,l3))
          if(dn(i,1).eq.' CA ') itemp = resno(nb(i,l3))
	  if(resno(j).eq.itemp)then
	    md(i,1)=i
	    md(i,2)=j
	  endif
	endif
121	enddo
	if(md(i,2).ne.i.and.md(i,2).ne.na(i,1).and.md(i,2).ne.na(i,2)
     &.and.md(i,2).ne.na(i,3).and.md(i,2).ne.na(i,4).and.
     &md(i,2).ne.na(i,5).and.md(i,2).gt.i.and.md(i,2).ne.n14(i,1).
     &and.md(i,2).ne.n14(i,2).and.md(i,2).ne.n14(i,3).and.
     &md(i,2).ne.n14(i,4).and.md(i,2).ne.n14(i,5))then
	n14(i,k14(i))=md(i,2)
	m14(i,k14(i))=dn(i,2)
	k14(i)=k14(i)+1
	endif

	endif
	enddo
	endif
	enddo
	enddo
	enddo
	close(unit=2)
c
c*****write energy parameters and atom pairs for 1-4, 1-5 & above interactions**
	open(unit=4,file=if2,status='unknown')
	do i=1,natom
	if(i.eq.1)then
	e1 = 0.1260
	v1 = 1.000
	v2 = 0.020
	write(4,100) i,e1,v1,v2,float(k15(i)/2),float(k14(i)-1)
	write(4,60) (n15(i,j),j=1,k15(i)),(n14(i,j),j=1,k14(i)-1)
	endif
	close(unit=7)
	close(unit=8)
	open(unit=7,file='PENERGYPARAM.lib',status='old')
	open(unit=8,file='PATOMTYPE.lib',status='old')
	do j=1,368
	read(8,90) tcode,scode,natm
	if(tcode.eq.rename(i))then
	do j1=1,natm
	read(8,80) at1,at2,at3,e1
	if(at1.eq.atname(i))then
	if(resno(i).eq.1.and.atname(i).eq.' N  ') at3 = 'NT'
	do j2=1,42
	read(7,70) atype,v1,v2
	if(atype.eq.at3) then
	write(4,100) i,e1,v1,v2,float(k15(i)/2),float(k14(i)-1)
	write(4,60) (n15(i,j5),j5=1,k15(i)),(n14(i,j4),j4=1,k14(i)-1)
	endif
	enddo
	endif
	enddo
	endif
	enddo
	enddo
	close(unit=7)
	close(unit=8)
c*******************write dihedral angles*************************************
	ntor = 0
	do i=1,natom
	if(resno(i).eq.resno(i-1)) goto 333
	close(unit=10)
	open(unit=10,file='PDIHEDS.lib',status='old')
	do j=1,111
	read(10,90) tcode,scode,natm
	if(tcode.eq.rename(i))then
	if(iopt.eq.1) natm = 2
	do j1=1,natm
	read(10,110) d1,d2,d3,d4
	do j2=i,natom
	if(d1.eq.atname(j2).and.tcode.eq.rename(j2).and.
     &resno(i).eq.resno(j2)) id1=j2
	if(d2.eq.atname(j2).and.tcode.eq.rename(j2).and.
     &resno(i).eq.resno(j2)) id2=j2
	if(d3.eq.atname(j2).and.tcode.eq.rename(j2).and.
     &resno(i).eq.resno(j2)) id3=j2
	if(d4.eq.atname(j2).and.resno(i).eq.resno(j2)) id4=j2
	if(d4.eq.atname(j2).and.d4.eq.' HO ') id4=j2
	enddo
	write(4,40) id1, id2, id3, id4
	if(icomment.eq.1)write(*,40) id1, id2, id3, id4
	ntor = ntor + 1
	enddo
	endif
	enddo
333	enddo
	write(31,'(A)')
	write(31,130)'Dihedral angles in Peptide:',ntor
        if(icomment.eq.1)write(*,*)'No. of Parameters(dihedral angles)
     & = ',ntor

c******************write hydrogen bond parameters*****************************

	ihp1=7557
	ihp2=2385
        mnhb=0
	do i=1,natom
	if(atname(i).eq.' H  '.or.atname(i).eq.' HN ') then
	do j=i,natom
	if(resno(i).lt.resno(j).and.atname(j).eq.' O  ')then
	write(4,120) i,j,ihp1,ihp2
        mnhb=mnhb+1
	endif
	enddo
	elseif(atname(i).eq.' O  ') then
	do j=i,natom
	if(resno(i+1).lt.resno(j).and.atname(j).eq.' H  ')then
	write(4,120) i,j,ihp1,ihp2
        mnhb=mnhb+1
c	write(31,*)i,j
	endif
	enddo
	endif
	enddo
	write(31,'(A)')
        write(31,140) 'H-bond pairs in Peptide   :',mnhb
        if(icomment.eq.1)print *,'pamppar-mnhb',mnhb
	close(unit=4)
c***************************************************************
600	return
	end

