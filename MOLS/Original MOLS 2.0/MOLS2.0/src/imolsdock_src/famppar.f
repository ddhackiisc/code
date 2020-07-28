c-----------------------------------------------------------------------------
c     library name : famppar.f   

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



c       Program to generate parameter file for protein flexible residues

        subroutine famppar()
        include 'mols.par'
        parameter (maxbnd = 819, jx = 25, maxr = 10, maxrat = 400)

	common /comment/icomment
	common /part/ipart,iptor
	common /pcrdb/nl,nptor,nphb,yy(maxatm,8)
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

	character*7 atom(maxrat)
        character*4 atnam(maxrat),patnam(maxr,maxrat),a(maxr,maxrat),
     &  an(maxrat,jx),b(maxr,maxrat),bn(maxrat,jx),c(maxr,maxrat),
     &  cn(maxrat,jx),dn(maxrat,jx),m14(maxrat,jx),d1,d2,d3,d4,
     &  at1,atna(maxrat),renam(maxrat)
        character*3 rnam(maxrat),prnam(maxr,maxrat),tcode
        character*2 cid(maxr),pcid(maxr,maxr),at3,atype
        character*1 scode,at2,acid(maxrat),aacid,abcid

        integer nl,rnum(maxrat),prnum(maxr,maxrat),ifat(maxrat),
     &  atno(maxrat),patno(maxr,maxrat),ki(maxrat),ma(maxrat,jx),
     &  ka(maxrat),na(maxrat,jx),mb(maxrat,jx),kb(maxrat),nb(maxrat,jx),
     &  nc(maxrat,jx),mc(maxrat,jx),kc(maxrat),kd(maxrat),
     &  nd(maxrat,maxrat),ke(maxrat),ne(maxrat,jx),k15(maxrat),
     &  n15(maxrat,jx),k14(maxrat),n14(maxrat,jx),md(maxrat,jx),
     &  iv(maxrat,4),jstart(maxrat,10),jend(maxrat,10),j1_4(maxrat,25),
     &  ihb1(maxrat),ihb2(maxrat),ch(maxrat),dh(maxrat),resn(maxrat)

        real rx(maxrat),ry(maxrat),rz(maxrat),prx(maxr,maxrat),
     &  pry(maxr,maxrat),prz(maxr,maxrat),energy,fx1(maxrat),
     &  fx2(maxrat),fx3(maxrat)

10      format(7x,i4,1x,a4,1x,a3,1x,a2,i4,4x,3f8.3)
15      format(a7,i4,1x,a4,1x,a4,a1,i4,4x,3f8.3)
20      format(a3,1x,a1,1x,i2)
30      format(8x,2a4)
40      format(10i4)
50      format(i4,1x,f3.1,1x,f3.1)
60      format(2x,20i4)
70      format(2x,a2,6x,f6.3,5x,f6.3)
80      format(a4,7x,a1,2x,a2,f8.4)
90      format(a3,1x,a1,1x,i2)
100     format(i4,1x,f8.4,f7.3,f7.3,f4.1,f4.1)
110     format(8x,4a4)
120     format(i4,i4,i6,i6)
130     format(5x,a42,2x,i4)
c-------------------------------------------------------------------
        if(ipart.gt.maxrat)then
        print *,'Number of flexible prot. res. atoms exceeding limit'
	print *,'Increase the value of maxrat, hb(),and steric()'
        STOP
        endif
c------------------------------------------------------------------
        jk=1
        jcount=0
        open(unit=1,file='partner',status='old')
        open(unit=2,file='partners',status='unknown')
        do is=1,ipart
        read(1,15)atom(is),atno(is),atna(is),renam(is),acid(is),
     &  resn(is),fx1(is),fx2(is),fx3(is)
        aacid=acid(is)
        abcid=acid(is-1)
        if(is.eq.1)then
        write(2,15)atom(is),atno(is),atna(is),renam(is),acid(is),
     &  jk,fx1(is),fx2(is),fx3(is)
        else
        IF(is.gt.1.and.resn(is).ne.resn(is-1).or.aacid(1:1).ne.
     &  abcid(1:1))THEN
        jk=jk+1

        write(2,15)atom(is),atno(is),atna(is),renam(is),acid(is),
     &  jk,fx1(is),fx2(is),fx3(is)

        ELSE

        write(2,15)atom(is),atno(is),atna(is),renam(is),acid(is),
     &  jk,fx1(is),fx2(is),fx3(is)

        ENDIF
        endif
        enddo
        close(2)
        close(1)
	nl=ipart
c-------read flexible residue details from flexres-----------------
        open(unit=2,file='partners',status='old')
        do il=1,nl
        read(2,10)atno(il),atnam(il),rnam(il),cid(il),rnum(il),
     &  rx(il),ry(il),rz(il)
        enddo
        close(2)
c-------copy rx ry rz ---------------------------------------------
	do i=1,nl
        yy(i,1)=rx(i)
        yy(i,2)=ry(i)
        yy(i,3)=rz(i)
        enddo
c-------store flexible residues details residuewise--NOT USED HERE-
        ir=1
        is=0
        do il=1,nl        
        if((il.gt.1).and.rnum(il).ne.rnum(il-1))then
        ifat(ir)=is
        ir=ir+1
        is=0
        endif
        is=is+1
        patno(ir,is)=atno(il)
        patnam(ir,is)=atnam(il)
        prnam(ir,is)=rnam(il)
        pcid(ir,is)=cid(il)
        prnum(ir,is)=rnum(il)
        prx(ir,is)=rx(il)
        pry(ir,is)=ry(il)
        prz(ir,is)=rz(il)
        enddo

        ifat(ir)=is !for final update of ifat
        
c*******find first neighbours***************************************
	 DO i=1,nl
	 ka(i)=1
	 open(unit=2,file='PALLCONN.lib',status='old')
	 rewind 2
	  do i1=1,maxbnd
	  read(2,20)tcode,scode,iatom
	   if(tcode.eq.rnam(i))then
	    do l=1,iatom
	    read(2,30)(a(i,j),j=1,2)
	     if(a(i,1).eq.atnam(i))then
	      do j=1,nl
	       if(atnam(j).eq.a(i,2).and.a(i,2).ne.'    ')then
	        if(a(i,1).eq.' C  '.and.a(i,2).eq.' N  ')then
	        itemp=rnum(i)+1
	        else
	        itemp=rnum(i)
	        endif
	
		if(rnum(j).eq.itemp)then
	        ma(i,1)=i
	        ma(i,2)=j
	        endif
	       endif
	      enddo
		if(ma(i,2).ne.na(i,1).and.ma(i,2).ne.na(i,2).and.
     &		ma(i,2).ne.na(i,3).and.ma(i,2).ne.na(i,4))then
		 na(i,ka(i))=ma(i,2)
		 an(i,ka(i))=a(i,2)
		 ka(i)=ka(i)+1
		endif
	     endif
	    enddo
	   endif
	  enddo
c*******find second neighbours********************************
	kb(i)=1
	 do l2=1,ka(i)-1
	 rewind 2
	  do i2=1,maxbnd
	  read(2,20)tcode,scode,iatom
	   if(tcode.eq.rnam(i))then
	    do l1=1,iatom
	    read(2,30)(b(i,j),j=1,2)
	     if(an(i,l2).eq.b(i,1))then
	      do j=1,nl
	       if(atnam(j).eq.b(i,2).and.b(i,2).ne.'   ')then
		if(b(i,1).eq.' C  '.and.b(i,2).eq.' N  ')then
		 itemp=rnum(i)+1
		else
		 itemp=rnum(i)
		endif
		
		if(b(i,1).eq.' N  ') itemp=rnum(na(i,l2))
	        if(rnum(j).eq.itemp)then
	         mb(i,1)=i
		 mb(i,2)=j
		endif
	       endif
              enddo
		
		if(mb(i,2).ne.i.and.mb(i,2).ne.nb(i,1).and.
     &		 mb(i,2).ne.nb(i,2).and.mb(i,2).ne.nb(i,3).and.
     &		 mb(i,2).ne.nb(i,4).and.mb(i,2).ne.nb(i,5).and.
     &		 mb(i,2).ne.nb(i,6))then
		 nb(i,kb(i))=mb(i,2)
		 bn(i,kb(i))=b(i,2)
		 kb(i)=kb(i)+1
		endif
             endif
	    enddo
	   endif
	  enddo
	 enddo

c*******find third neighbours********************************
	kc(i)=1
	 do l3=1,kb(i)-1
	 rewind 2
	  do i3=1,maxbnd
	  read(2,20)tcode,scode,iatom
	   if(tcode.eq.rnam(nb(i,l3)))then
	    do l2=1,iatom
	     read(2,30)(c(i,j),j=1,2)
	     if(bn(i,l3).eq.c(i,1))then
              do j=1,nl
	       if(atnam(j).eq.c(i,2).and.c(i,2).ne.'    ')then
                if(c(i,1).eq.' C  '.and.c(i,2).eq.' N  ')then
	         itemp=rnum(i)+1
	        else
	         itemp=rnum(i)
	        endif
		if(c(i,1).eq.' N  ') itemp=rnum(nb(i,l3))              
		if(c(i,1).eq.' CA ') itemp=rnum(nb(i,l3))
 		 if(rnum(j).eq.itemp)then
	          mc(i,1)=i
	          mc(i,2)=j
		 endif
               endif
              enddo
     
	     if(mc(i,2).ne.i.and.mc(i,2).ne.na(i,1).and.
     & 	mc(i,2).ne.na(i,2).and.mc(i,2).ne.na(i,3).and.
     & 	mc(i,2).ne.na(i,4).and.mc(i,2).ne.na(i,5).and.
     & 	mc(i,2).gt.i.and.mc(i,2).ne.nc(i,1).and.
     & 	mc(i,2).ne.nc(i,2).and.mc(i,2).ne.nc(i,3).and.
     & 	mc(i,2).ne.nc(i,4).and.mc(i,2).ne.nc(i,5))then
		nc(i,kc(i))=mc(i,2)
		cn(i,kc(i))=c(i,2)
		kc(i)=kc(i)+1
	    endif
	     endif
            enddo
           endif
	  enddo
	 enddo

	ENDDO

c*******to find fourth and above neighbours*****************
	 DO i=1,nl
	 kd(i)=1
	  do j=i+1,nl
	   do j1=1,ka(i)-1
	    if(na(i,j1).eq.j)goto 111
	   enddo   
	   do j2=1,kb(i)-1
	    if(nb(i,j2).eq.j)goto 111
	   enddo
	   do j3=1,kc(i)-1
	    if(nc(i,j3).eq.j)goto 111
	   enddo
	   
	   nd(i,kd(i))=j
	   kd(i)=kd(i)+1
111	  enddo
	ENDDO
c*******to count the atom pairs fourth and above neighbours*
	  DO i=1,nl
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
	ENDDO

c***************************************************************
c       to exclude the atompairs in around peptide bond
c       and in the rings(1-5 pairs)
C**************************************************************
        do i=1,nl
        k15(i)=0
c       write(*,*) (ne(i,j),j=1,ke(i))
        do j=1,ke(i)
        j1=ne(i,j)
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.' HE1'.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.' HE2'.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.' HH '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' CG '
     &.and.atnam(j1).eq.' OH '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' CG '
     &.and.atnam(j1).eq.' HH '.and.rnum(i).eq.rnum(j1)) goto 131

        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' CD1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' HD1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' CD2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' HD2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' CE1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' HE1'.and.
     &rnum(i).eq.rnum(j1))  goto 151
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' CE2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' HE2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' CZ '.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' OH '.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TYR '.and.atnam(i).eq.' HH '.and.
     &rnum(i).eq.rnum(j1)) goto 151

        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.' HE1'.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.' HE2'.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.' HZ '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' CG '
     &.and.atnam(j1).eq.' HZ '.and.rnum(i).eq.rnum(j1)) goto 131

        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' CD1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' HD1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' CD2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' HD2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' CE1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' HE1'.and.
     &rnum(i).eq.rnum(j1))  goto 151
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' CE2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' HE2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' CZ '.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'PHE '.and.atnam(i).eq.' HZ '.and.
     &rnum(i).eq.rnum(j1)) goto 151


        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.' HE1'.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.' HE3'.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.' HH2'.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CG '
     &.and.atnam(j1).eq.' HZ2'.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CG '
     &.and.atnam(j1).eq.' HZ3'.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CG '
     &.and.atnam(j1).eq.' HH2'.and.rnum(i).eq.rnum(j1)) goto 131

        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CD1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' HD1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CD2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' NE1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' HE1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CE2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CE3'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' HE3'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CZ2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' HZ2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CZ3'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' HZ3'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' CH2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'TRP '.and.atnam(i).eq.' HH2'.and.
     &rnum(i).eq.rnum(j1)) goto 151

        if(rnam(i).eq.'HIS '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.' HE1'.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'HIS '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.' HE2'.and.rnum(i).eq.rnum(j1)) goto 131

        if(rnam(i).eq.'HIS '.and.atnam(i).eq.' HD2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'HIS '.and.atnam(i).eq.' ND1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'HIS '.and.atnam(i).eq.' CD2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'HIS '.and.atnam(i).eq.' HD2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'HIS '.and.atnam(i).eq.' CE1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'HIS '.and.atnam(i).eq.' HE1'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'HIS '.and.atnam(i).eq.' NE2'.and.
     &rnum(i).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'HIS '.and.atnam(i).eq.' HE2'.and.
     &rnum(i).eq.rnum(j1)) goto 151

        if(rnam(j1).eq.'PRO '.and.atnam(i).eq.' C  '
     &.and.(rnum(i)+1).eq.rnum(j1)) goto 151
        if(rnam(j1).eq.'PRO '.and.atnam(i).eq.' O  '
     &.and.(rnum(i)+1).eq.rnum(j1)) goto 151
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' N  '
     &.and.atnam(j1).eq.'1HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' N  '
     &.and.atnam(j1).eq.'2HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' CA '
     &.and.atnam(j1).eq.'1HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' CA '
     &.and.atnam(j1).eq.'2HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.'1HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' CB '
     &.and.atnam(j1).eq.'2HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.'1HB '
     &.and.atnam(j1).eq.'2HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.'1HB '
     &.and.atnam(j1).eq.'1HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.'1HB '
     &.and.atnam(j1).eq.'2HD '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.'1HB '
     &.and.atnam(j1).eq.'1HD '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.'2HB '
     &.and.atnam(j1).eq.'2HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.'2HB '
     &.and.atnam(j1).eq.'1HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.'2HB '
     &.and.atnam(j1).eq.'2HD '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.'2HB '
     &.and.atnam(j1).eq.'1HD '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' CG '
     &.and.atnam(j1).eq.'2HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' CG '
     &.and.atnam(j1).eq.'1HG '.and.rnum(i).eq.rnum(j1)) goto 131
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.'1HG '
     &.and.atnam(j1).eq.'2HG '.and.rnum(i).eq.rnum(j1)) goto 131

        goto 161

151     do j2=i+1,nl
          if(atnam(j2).eq.' C  ') then
                k15(i)=k15(i)+1
                n15(i,k15(i))=j2
                k15(i)=k15(i)+1
                n15(i,k15(i))=nl
                goto 141
          endif
        enddo

161     k15(i)=k15(i)+1
        n15(i,k15(i))=ne(i,j)
131     enddo
141     enddo

c******************************************************************
c       to exclude the atompairs in around peptide bond
c       and in the rings(in 1-4 pairs) and pickup needed 1-4 pairs
C******************************************************************
        do i = 1,nl
        k14(i)=1
        do l3=1,kb(i)-1
        rewind 2
        do i3=1,maxbnd
        read(2,20)tcode,scode,iatom
        if(tcode.eq.rnam(nb(i,l3)))then
        do l2=1,iatom
        read(2,30) (dn(i,j),j=1,2)
        if(bn(i,l3).eq.dn(i,1))then
        do j=1,nl
        if(atnam(i).eq.' CA '.and.atnam(j).eq.' CA '.and.
     &(rnum(i)+1).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CA '.and.atnam(j).eq.' H  '.and.
     &(rnum(i)+1).eq.rnum(j)) goto 121
        if(atnam(i).eq.' O  '.and.atnam(j).eq.' CA '.and.
     &(rnum(i)+1).eq.rnum(j)) goto 121
        if(atnam(i).eq.' O  '.and.atnam(j).eq.' H  '.and.
     &(rnum(i)+1).eq.rnum(j)) goto 121

        if(atnam(i).eq.' CB '.and.atnam(j).eq.' HD1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TYR ') goto 121
        if(atnam(i).eq.' CB '.and.atnam(j).eq.' HD2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TYR ') goto 121
        if(atnam(i).eq.' CB '.and.atnam(j).eq.' CE1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TYR ') goto 121
        if(atnam(i).eq.' CB '.and.atnam(j).eq.' CE2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TYR ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' CZ '.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TYR ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' HE1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TYR ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' HE2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TYR ') goto 121
        if(atnam(i).eq.' CE1'.and.atnam(j).eq.' HE2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TYR ') goto 121
        if(atnam(i).eq.' CE1'.and.atnam(j).eq.' HH '.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TYR ') goto 121
        if(atnam(i).eq.' CE2'.and.atnam(j).eq.' HH '.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TYR ') goto 121
        if(atnam(i).eq.' CD1'.and.rnam(i).eq.'TYR '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HD1'.and.rnam(i).eq.'TYR '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CD2'.and.rnam(i).eq.'TYR '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HD2'.and.rnam(i).eq.'TYR '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HE1'.and.rnam(i).eq.'TYR '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HE2'.and.rnam(i).eq.'TYR '.and.
     &rnum(i).eq.rnum(j)) goto 121

        if(atnam(i).eq.' CB '.and.atnam(j).eq.' HD1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'PHE ') goto 121
        if(atnam(i).eq.' CB '.and.atnam(j).eq.' HD2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'PHE ') goto 121
        if(atnam(i).eq.' CB '.and.atnam(j).eq.' CE1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'PHE ') goto 121
        if(atnam(i).eq.' CB '.and.atnam(j).eq.' CE2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'PHE ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' CZ '.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'PHE ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' HE1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'PHE ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' HE2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'PHE ') goto 121

        if(atnam(i).eq.' CD1'.and.rnam(i).eq.'PHE '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HD1'.and.rnam(i).eq.'PHE '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CD2'.and.rnam(i).eq.'PHE '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HD2'.and.rnam(i).eq.'PHE '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CE1'.and.rnam(i).eq.'PHE '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CE2'.and.rnam(i).eq.'PHE '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HE1'.and.rnam(i).eq.'PHE '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HE2'.and.rnam(i).eq.'PHE '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CZ '.and.rnam(i).eq.'PHE '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HZ '.and.rnam(i).eq.'PHE '.and.
     &rnum(i).eq.rnum(j)) goto 121

        if(atnam(i).eq.' CB '.and.atnam(j).eq.' HD1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TRP ') goto 121
        if(atnam(i).eq.' CB '.and.atnam(j).eq.' NE1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TRP ') goto 121
        if(atnam(i).eq.' CB '.and.atnam(j).eq.' CE2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TRP ') goto 121
        if(atnam(i).eq.' CB '.and.atnam(j).eq.' CE3'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TRP ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' CE2 '.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TRP ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' HE1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TRP ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' CZ2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TRP ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' HE3'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TRP ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' NE1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TRP ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' CZ3'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TRP ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' HZ3'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'TRP ') goto 121


        if(atnam(i).eq.' CD1'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HD1'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CD2'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' NE1'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HE1'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CE2'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CE3'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HE3'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CZ2'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HZ2'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CZ3'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HZ3'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CH2'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HH2'.and.rnam(i).eq.'TRP '.and.
     &rnum(i).eq.rnum(j)) goto 121

        if(atnam(i).eq.' CB '.and.atnam(j).eq.' CE1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'HIS ') goto 121
        if(atnam(i).eq.' CB '.and.atnam(j).eq.' NE2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'HIS ') goto 121
        if(atnam(i).eq.' CB '.and.atnam(j).eq.' HD2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'HIS ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' NE2 '.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'HIS ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' HE1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'HIS ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' HE2'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'HIS ') goto 121
        if(atnam(i).eq.' CG '.and.atnam(j).eq.' CE1'.and.
     &rnum(i).eq.rnum(j).and.rnam(i).eq.'HIS ') goto 121

        if(atnam(i).eq.' ND1'.and.rnam(i).eq.'HIS '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CE1'.and.rnam(i).eq.'HIS '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HE1'.and.rnam(i).eq.'HIS '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' NE2'.and.rnam(i).eq.'HIS '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HE2'.and.rnam(i).eq.'HIS '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CD2'.and.rnam(i).eq.'HIS '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' HD2'.and.rnam(i).eq.'HIS '.and.
     &rnum(i).eq.rnum(j)) goto 121

        if(atnam(i).eq.' C  '.and.rnam(j).eq.'PRO ') goto 121

        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' N  '.and.
     &atnam(j).eq.'1HB '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' N  '.and.
     &atnam(j).eq.'2HB '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' N  '.and.
     &atnam(j).eq.' CG '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' N  '.and.
     &atnam(j).eq.' CB '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' N  '.and.
     &atnam(j).eq.'1HG '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' N  '.and.
     &atnam(j).eq.'2HG '.and.rnum(i).eq.rnum(j)) goto 121

        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' HA '.and.
     &atnam(j).eq.'1HB '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' HA '.and.
     &atnam(j).eq.'2HB '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' HA '.and.
     &atnam(j).eq.' CG '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' HA '.and.
     &atnam(j).eq.' CB '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' HA '.and.
     &atnam(j).eq.'1HG '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' HA '.and.
     &atnam(j).eq.'2HG '.and.rnum(i).eq.rnum(j)) goto 121

        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' CB '.and.
     &atnam(j).eq.'1HD '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' CB '.and.
     &atnam(j).eq.'2HD '.and.rnum(i).eq.rnum(j)) goto 121
        if(rnam(i).eq.'PRO '.and.atnam(i).eq.' CB '.and.
     &atnam(j).eq.' CD '.and.rnum(i).eq.rnum(j)) goto 121

        if(atnam(i).eq.' CA '.and.rnam(i).eq.'PRO '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.'1HB '.and.rnam(i).eq.'PRO '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.'2HB '.and.rnam(i).eq.'PRO '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CG '.and.rnam(i).eq.'PRO '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.'1HG '.and.rnam(i).eq.'PRO '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.'2HG '.and.rnam(i).eq.'PRO '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.' CD '.and.rnam(i).eq.'PRO '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.'1HD '.and.rnam(i).eq.'PRO '.and.
     &rnum(i).eq.rnum(j)) goto 121
        if(atnam(i).eq.'2HD '.and.rnam(i).eq.'PRO '.and.
     &rnum(i).eq.rnum(j)) goto 121


        if(atnam(j).eq.dn(i,2).and.dn(i,2).ne.'    ')then
          if(dn(i,1).eq.' C  '.and.dn(i,2).eq.' N  ')then
            itemp=rnum(i)+1
          else
            itemp=rnum(i)
          endif
          if(dn(i,1).eq.' N  ') itemp = rnum(nb(i,l3))
          if(dn(i,1).eq.' CA ') itemp = rnum(nb(i,l3))
          if(rnum(j).eq.itemp)then
            md(i,1)=i
            md(i,2)=j
          endif
        endif
121     enddo

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
c*****write energy parameters and atom pairs for 1-4, 1-5 & above interactions**
        open(unit=4,file=if1,status='unknown')
        do i=1,nl
        close(unit=7)
        close(unit=8)
        open(unit=7,file='ENERGYPARAM.lib',status='old')
        open(unit=8,file='PATOMTYPE.lib',status='old')
        do j=1,366
        read(8,90) tcode,scode,natm
        if(tcode.eq.rnam(i))then
        do j1=1,natm
        read(8,80) at1,at2,at3,e1
        if(at1.eq.atnam(i))then
        if(rnum(i).eq.1.and.atnam(i).eq.' N  ') at3 = 'NT'
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
        nptor = 0
        do i=1,nl
        if(rnum(i).eq.rnum(i-1)) goto 333
        close(unit=10)
        open(unit=10,file='SCDIHEDS.lib',status='old')
        do j=1,111
        read(10,90) tcode,scode,natm
        if(tcode.eq.rnam(i))then
        if(iopt.eq.1) natm = 2
        do j1=1,natm
        read(10,110) d1,d2,d3,d4
        do j2=i,nl
        if(d1.eq.atnam(j2).and.tcode.eq.rnam(j2).and.
     &rnum(i).eq.rnum(j2)) id1=j2
        if(d2.eq.atnam(j2).and.tcode.eq.rnam(j2).and.
     &rnum(i).eq.rnum(j2)) id2=j2
        if(d3.eq.atnam(j2).and.tcode.eq.rnam(j2).and.
     &rnum(i).eq.rnum(j2)) id3=j2
        if(d4.eq.atnam(j2).and.rnum(i).eq.rnum(j2)) id4=j2
        if(d4.eq.atnam(j2).and.d4.eq.' HO ') id4=j2
        enddo
        write(4,40) id1, id2, id3, id4
        nptor = nptor + 1
        enddo
        endif
        enddo
333     enddo
	write(31,'(A)')
        write(31,130)'Dihedral angles in Protein Partners      :',nptor
c******************write hydrogen bond parameters*****************************
        ihp1=7557
        ihp2=2385
        nphb=0
        do i=1,nl
        if(atnam(i).eq.' H  '.or.atnam(i).eq.' HN ') then
        do j=i,nl
        if(rnum(i).lt.rnum(j).and.atnam(j).eq.' O  ')then
        write(4,120) i,j,ihp1,ihp2
        nphb=nphb+1
        endif
        enddo
        elseif(atnam(i).eq.' O  ') then
        do j=i,nl
        if(rnum(i+1).lt.rnum(j).and.atnam(j).eq.' H  ')then
        write(4,120) i,j,ihp1,ihp2
        nphb=nphb+1
        endif
        enddo
        endif
        enddo
	write(31,'(A)')
        write(31,130)'H-bond pairs in Protein Partners         :',nphb
        close(unit=4)
	return
	end
c*******read protein parameters*********************************
	subroutine finput()
	include 'mols.par'

	common /comment/icomment
	common /provectors/ivp(maxpar,4)
	common /pcrdb/nl,nptor,nphb,yy(maxatm,8)
	common /phb/iphb1(maxhb),iphb2(maxhb),ch(maxhb),dh(maxhb)
	common /pranges/jpstart(maxatm,10),jpend(maxatm,10),
     &  jp1_4(maxatm,25)
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

	character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

	open(unit=1,file=if1,status='old')
	do i=1,nl
	read(1,*)atn,(yy(i,j),j=4,8)
	i1=ifix(yy(i,7))
	i2=ifix(yy(i,8))
	read(1,*)(jpstart(i,j),jpend(i,j),j=1,i1),
     &  (jp1_4(i,k),k=1,i2)
	enddo
	
	do i=1,nptor
	read(1,*)(ivp(i,j),j=1,4)
	enddo

	do i=1,nphb
	read(1,*)iphb1(i),iphb2(i),ch(i),dh(i)
	enddo
	close(1)
	return
	end
c**************************************************************

