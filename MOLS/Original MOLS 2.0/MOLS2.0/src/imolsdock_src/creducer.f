c-----------------------------------------------------------------------------
c     library name : creducer.f   

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

c**********************************************************************
c-------copy protonated receptor file into rec.pdb in MOLS pdb format
c---------------------------------------------------------------------- 
        subroutine creducer(inat)
        include 'mols.par'
        parameter (maxn=20000)

	common /comment/icomment
	common /patom/ipatom
	common /native/nat,pepfile
	common /par/natom,ntor,nhb,ns,lres
        common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        character str*80,pepfile*128
        character*4 rATOM(maxn),mATOM(maxn),rresnam(maxn),
     &  mresnam(maxn)
        character*6 ratnam(maxn),matnam(maxn)
        character*1 rcid(maxn),mcid(maxn)
	character copy*80,rmcopy*80
        integer ratno(maxn),matno(maxn),rresno(maxn),
     &  mresno(maxn),inat
        real rx(maxn),ry(maxn),rz(maxn),mx(maxn),
     &  my(maxn),mz(maxn)
                
20      format(a80)
30      format(a4,3x,i4,a6,a4,a1,i4,4x,3f8.3,30x)

        ino=0
        idum=0
	IF(inat.eq.1) then
        open(unit=1,file='pfile_mod.pdb',status='old')
        open(unit=2,file='rec_temp.pdb',status='unknown')
        do i=1,maxn
        read(1,20,end=99)str
        if(str(1:4).eq.'ATOM')then
        write(2,20)str
        ino=ino+1
        endif
        enddo
99      if(icomment.eq.1)write(*,*)'protonated receptor :',ino
        ipatom=ino
        close(2)
        close(1)
	ELSE IF(inat.eq.2)then
	open(unit=1,file=pepfile,status='old')
        do i=1,maxn
        read(1,20,end=100)str
        ino=ino+1
        enddo
100	write(*,*)'native peptide atoms:',ino
        close(1)
	ENDIF
c---------------------------------------------------------
	if(inat.eq.1)then
	open(unit=1,file='rec_temp.pdb',status='old')
	open(unit=2,file='rec.pdb',status='unknown')

	else if(inat.eq.2)then
	open(unit=1,file=pepfile,status='old')
	open(unit=2,file='pepfile.pdb',status='unknown')
	endif
        do i=1,ino
        read(1,30)rATOM(i),ratno(i),ratnam(i),rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
cs      if(icomment.eq.1)then
cs      write(*,30)rATOM(i),ratno(i),ratnam(i),rresnam(i),rcid(i),
cs     &  rresno(i),rx(i),ry(i),rz(i)
cs        endif !icomment
        if(rresnam(i).eq.'ALA '.and.ratnam(i).eq.'  HB1 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ALA '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ALA '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 3HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'GLU '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'GLU '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'GLU '.and.ratnam(i).eq.'  HG2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'GLU '.and.ratnam(i).eq.'  HG3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'THR '.and.ratnam(i).eq.' HG21 ')then
        write(2,30)rATOM(i),ratno(i),' 1HG2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'THR '.and.ratnam(i).eq.' HG22 ')then
        write(2,30)rATOM(i),ratno(i),' 2HG2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'THR '.and.ratnam(i).eq.' HG23 ')then
        write(2,30)rATOM(i),ratno(i),' 3HG2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'PHE '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'PHE '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.'  HG2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.'  HG3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.'  HD2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HD  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.'  HD3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HD  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
c        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.' HG12 ')then
c        write(2,30)rATOM(i),ratno(i),' 2HH1 ',rresnam(i),rcid(i),
c     &  rresno(i),rx(i),ry(i),rz(i)        
	else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.' HH12 ')then
        write(2,30)rATOM(i),ratno(i),' 2HH1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
c        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.' HG11 ')then
c        write(2,30)rATOM(i),ratno(i),' 1HH1 ',rresnam(i),rcid(i),
c     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.' HH11 ')then
        write(2,30)rATOM(i),ratno(i),' 1HH1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i) 
c        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.' HG21 ')then
c        write(2,30)rATOM(i),ratno(i),' 1HH2 ',rresnam(i),rcid(i),
c     &  rresno(i),rx(i),ry(i),rz(i) 
        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.' HH21 ')then
        write(2,30)rATOM(i),ratno(i),' 1HH2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
c        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.' HG22 ')then
c        write(2,30)rATOM(i),ratno(i),' 2HH2 ',rresnam(i),rcid(i),
c     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.' HH22 ')then
        write(2,30)rATOM(i),ratno(i),' 2HH2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)

cs      else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.'  NH1 ')then
cs      write(2,30)rATOM(i),ratno(i),' 1NH1 ',rresnam(i),rcid(i),
cs   &  rresno(i),rx(i),ry(i),rz(i)
cs      else if(rresnam(i).eq.'ARG '.and.ratnam(i).eq.'  NH2 ')then
cs      write(2,30)rATOM(i),ratno(i),' 2NH  ',rresnam(i),rcid(i),
cs   &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ASN '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ASN '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ASN '.and.ratnam(i).eq.' HD21 ')then
        write(2,30)rATOM(i),ratno(i),' 1HD2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ASN '.and.ratnam(i).eq.' HD22 ')then
        write(2,30)rATOM(i),ratno(i),' 2HD2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ASP '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ASP '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'CYS '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'CYS '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'CYS '.and.ratnam(i).eq.'  HG  ')then
        write(2,30)rATOM(i),ratno(i),' 3HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'GLN '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'GLN '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'GLN '.and.ratnam(i).eq.'  HG2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'GLN '.and.ratnam(i).eq.'  HG3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'GLN '.and.ratnam(i).eq.' HE22 ')then
        write(2,30)rATOM(i),ratno(i),' 2HE2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i) 
        else if(rresnam(i).eq.'GLN '.and.ratnam(i).eq.' HE21 ')then
        write(2,30)rATOM(i),ratno(i),' 1HE2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'GLY '.and.ratnam(i).eq.'  HA2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HA  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i) 
        else if(rresnam(i).eq.'GLY '.and.ratnam(i).eq.'  HA3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HA  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'HIS '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'HIS '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ILE '.and.ratnam(i).eq.' HG12 ')then
        write(2,30)rATOM(i),ratno(i),' 1HG1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ILE '.and.ratnam(i).eq.' HG13 ')then
        write(2,30)rATOM(i),ratno(i),' 2HG1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ILE '.and.ratnam(i).eq.' HD13 ')then
        write(2,30)rATOM(i),ratno(i),' 3HD1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ILE '.and.ratnam(i).eq.' HD11 ')then
        write(2,30)rATOM(i),ratno(i),' 1HD1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ILE '.and.ratnam(i).eq.' HD12 ')then
        write(2,30)rATOM(i),ratno(i),' 2HD1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ILE '.and.ratnam(i).eq.' HG23 ')then
        write(2,30)rATOM(i),ratno(i),' 3HG2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ILE '.and.ratnam(i).eq.' HG21 ')then
        write(2,30)rATOM(i),ratno(i),' 1HG2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'ILE '.and.ratnam(i).eq.' HG22 ')then
        write(2,30)rATOM(i),ratno(i),' 2HG2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LEU '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LEU '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LEU '.and.ratnam(i).eq.' HD13 ')then
        write(2,30)rATOM(i),ratno(i),' 3HD1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LEU '.and.ratnam(i).eq.' HD11 ')then
        write(2,30)rATOM(i),ratno(i),' 1HD1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LEU '.and.ratnam(i).eq.' HD12 ')then
        write(2,30)rATOM(i),ratno(i),' 2HD1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LEU '.and.ratnam(i).eq.' HD23 ')then
        write(2,30)rATOM(i),ratno(i),' 3HD2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LEU '.and.ratnam(i).eq.' HD21 ')then
        write(2,30)rATOM(i),ratno(i),' 1HD2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LEU '.and.ratnam(i).eq.' HD22 ')then
        write(2,30)rATOM(i),ratno(i),' 2HD2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)        
        else if(rresnam(i).eq.'LYS '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LYS '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LYS '.and.ratnam(i).eq.'  HG2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LYS '.and.ratnam(i).eq.'  HG3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LYS '.and.ratnam(i).eq.'  HD2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HD  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LYS '.and.ratnam(i).eq.'  HD3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HD  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LYS '.and.ratnam(i).eq.'  HE2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HE  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LYS '.and.ratnam(i).eq.'  HE3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HE  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LYS '.and.ratnam(i).eq.'  HZ3 ')then
        write(2,30)rATOM(i),ratno(i),' 3HZ  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LYS '.and.ratnam(i).eq.'  HZ2 ')then
        write(2,30)rATOM(i),ratno(i),' 2HZ  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'LYS '.and.ratnam(i).eq.'  HZ1 ')then
        write(2,30)rATOM(i),ratno(i),' 1HZ  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
c       idum=idum+1
        else if(rresnam(i).eq.'MET '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'MET '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'MET '.and.ratnam(i).eq.'  HG2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'MET '.and.ratnam(i).eq.'  HG3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'MET '.and.ratnam(i).eq.'  HE1 ')then
        write(2,30)rATOM(i),ratno(i),' 1HE  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'MET '.and.ratnam(i).eq.'  HE2 ')then
        write(2,30)rATOM(i),ratno(i),' 2HE  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'MET '.and.ratnam(i).eq.'  HE3 ')then
        write(2,30)rATOM(i),ratno(i),' 3HE  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'PRO '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'PRO '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'PRO '.and.ratnam(i).eq.'  HG2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'PRO '.and.ratnam(i).eq.'  HG3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HG  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'PRO '.and.ratnam(i).eq.'  HD2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HD  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'PRO '.and.ratnam(i).eq.'  HD3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HD  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'SER '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'SER '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'TRP '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'TRP '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'TYR '.and.ratnam(i).eq.'  HB2 ')then
        write(2,30)rATOM(i),ratno(i),' 1HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'TYR '.and.ratnam(i).eq.'  HB3 ')then
        write(2,30)rATOM(i),ratno(i),' 2HB  ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'VAL '.and.ratnam(i).eq.' HG11 ')then
        write(2,30)rATOM(i),ratno(i),' 1HG1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'VAL '.and.ratnam(i).eq.' HG12 ')then
        write(2,30)rATOM(i),ratno(i),' 2HG1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'VAL '.and.ratnam(i).eq.' HG13 ')then
        write(2,30)rATOM(i),ratno(i),' 3HG1 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'VAL '.and.ratnam(i).eq.' HG21 ')then
        write(2,30)rATOM(i),ratno(i),' 1HG2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'VAL '.and.ratnam(i).eq.' HG22 ')then
        write(2,30)rATOM(i),ratno(i),' 2HG2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(rresnam(i).eq.'VAL '.and.ratnam(i).eq.' HG23 ')then
        write(2,30)rATOM(i),ratno(i),' 3HG2 ',rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(ratnam(i).eq.'  OXT ')then
        write(2,30)rATOM(i),ratno(i),'  HO  ','HO  ',rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        else if(ratnam(i).eq.'  H1  '.or.ratnam(i).eq.'  H2  '.or.
     &  ratnam(i).eq.'  H3  ')then
        idum=idum+1
        else if(ratnam(i).ne.'  H1  '.and.ratnam(i).ne.'  H2  '.and.
     &  ratnam(i).ne.'  H3  ')then 
        write(2,30)rATOM(i),ratno(i),ratnam(i),rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        endif        

        enddo
        close(2)
        close(1)
        if(icomment.eq.1)print *,'idum',idum,'ino',ino
        ino = ino-idum
        ipatom=ino
	if(inat.eq.1)goto 299        
c----------------------------------------------------------------------
        open(unit=1,file='pepfile.pdb',status='old')        
        do  i=1,ino
        read(1,30)rATOM(i),ratno(i),ratnam(i),rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        enddo
        close(1)        
c-------copying with new atom number-----------------------------------
	open(unit=2,file='pepfile.copy',status='unknown')
	do j=1,ino
	write(2,30)rATOM(j),j,ratnam(j),rresnam(j),rcid(j),
     &  rresno(j),rx(j),ry(j),rz(j)
	enddo
	close(2)
c------copy back to original file--------------------------------------
	open(unit=2,file='pepfile.pdb',status='unknown')
        do j=1,ino
        write(2,30)rATOM(j),ratno(j),ratnam(j),rresnam(j),rcid(j),
     &  rresno(j),rx(j),ry(j),rz(j)
        enddo
        close(2)

	copy='rm pepfile.copy'
	call system(copy)

c----------------------------------------------------------------------
        open(unit=1,file='pf1.temp',status='old')
        do i=1,natom
        read(1,30)mATOM(i),matno(i),matnam(i),mresnam(i),mcid(i),
     &  mresno(i),mx(i),my(i),mz(i)
        enddo
        close(1)

c-------changing the residue number ------------------------------------
        ir=1
        open(unit=2,file='pepfile.final2',status='unknown')
        do  i=1,ino
        if(i.eq.1)then
        write(2,30)rATOM(i),ratno(i),ratnam(i),rresnam(i),rcid(i),
     &  ir,rx(i),ry(i),rz(i)
        else
         if(rresno(i).ne.rresno(i-1))then
          ir=ir+1
          write(2,30)rATOM(i),ratno(i),ratnam(i),rresnam(i),rcid(i),
     &    ir,rx(i),ry(i),rz(i)
         else
         write(2,30)rATOM(i),ratno(i),ratnam(i),rresnam(i),rcid(i),
     &   ir,rx(i),ry(i),rz(i)
         endif
        endif
        enddo
        close(2)
c-----------------------------------------------------------------------       
 
        open(unit=3,file='pepfile.final2',status='old')
        do i=1,ino
        read(3,30)rATOM(i),ratno(i),ratnam(i),rresnam(i),rcid(i),
     &  rresno(i),rx(i),ry(i),rz(i)
        enddo
        close(3)
        
        iino=0
        iano=1
        open(unit=4,file=pf1,status='unknown')
        do i=1,natom 
        do j=1,ino
        if(mresnam(i).eq.rresnam(j).and.matnam(i).eq.ratnam(j)
     &  .and.mresno(i).eq.rresno(j))then        
        write(4,30)rATOM(j),iano,ratnam(j),rresnam(j),rcid(j),
     &  mresno(i),rx(j),ry(j),rz(j)
        iino=iino+1
        iano=iano+1
	else if(i.eq.natom.and.ratnam(j).eq.'  HO  '
     &  .and.rresnam(j).eq.'HO  ' )then
	write(4,30)rATOM(j),iano,ratnam(j),rresnam(j),rcid(j),
     &  mresno(i),rx(j),ry(j),rz(j)
	iino=iino+1
        iano=iano+1
        endif
        enddo !ino
        enddo !natom
        close(4)
        if(icomment.eq.1)print *,'pf1-new atoms:',iino
        natom=iino-1
        
299     return
        end
        
c----------------------------------------------------------------------

