c-----------------------------------------------------------------------------
c     library name : pvarinit.f   

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

c	Peptide-Protein Docking
c	Variable Initialization for Peptide and Protein
	subroutine pvarinit
	include 'mols.par'

	common /comment/icomment
	common /vectors/iv(maxpar,4)
	common /pepcrda/xpep(maxatm,8)
        common /pepcrdb/ypep(maxatm,8)
	common /native/nat,pepfile
	common /order/nn,mm
	common /pepdbat/pepatom(maxatm),pepele(maxatm)
	common /var/ivat(maxres*8,4),vty(maxres*8),ivres(maxres*8)
	common /ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,25)
	common /par/natom,ntor,nhb,ns,lres
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
	common /hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
	common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common /atm/iatmno(maxatm),atname(maxatm),rename(maxatm),
     &              iresno(maxatm)
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
	common /plp/pat(mnatp),lat(maxatm),patnam(mnatp),presnam(mnatp),
     &  pchid(mnatp),tpres(25),tatnam(25,50),natp,tatyp(25,50),
     &  ntatp(25),px(mnatp,3)
                
  	external pdihedr

	character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
	character*30 pepatom*30, pepele*24
	character dd1*4,dd2*4,dd3*4,dd4*4,dd5*4,vty*4,tcode*4,scode*1,
     &  atname*4, rename*4,fseq*1
     	character xatnam*4, xresnam*3,lat*2,pat*2,tatyp*2
	character patnam*4,presnam*3,tpres*3,tatnam*4

	dimension x_one(maxatm,3),x_two(maxatm,3)
	integer iff,sulphur,nat,natp

	nn = ntor

10      format(8x,i3,1x,a4,1x,a4,3x,i3)
41	format(i5,2x,4i5,2x,a4)
90	format(a4,a1,1x,i2)
112     format(8x,4a4,2x,a4,3f7.2)


	open(unit=1,file=pf1,status='old')
	do i=1,natom
		read(1,10) iatmno(i),atname(i),rename(i),iresno(i)
	enddo
	close(unit=1)
        
        if(icomment.eq.1)then
        print *,'iff-pvarinit',iff
        print *,'iopt',iopt
        print *,'lres',lres
        endif

	call ppdb

	call pinput
	if(nat.eq.1)goto 315
c-----------------define dihedral angles atoms-------------------------
	nvar = 0
	do i = 1,lres
	rewind 10
	open(unit=10,file='PVAR.lib',status='old')
	do j=1,111
	read(10,90) tcode,scode,natm
	if(scode.eq.fseq(i))then
        if(icomment.eq.1) print *,'scode',scode
	if(iopt.eq.1) natm = 2
	do j1=1,natm
	idd1= 0
	idd2= 0
	idd3= 0
	idd4= 0
	tu0 = 0.0
	tsn = 0.0
	ttn = 0.0
	read(10,112) dd1,dd2,dd3,dd4,dd5,tu0,tsn,ttn
	if(i.eq.1.and.dd5.eq.'phi ')goto 251
	if(i.eq.lres.and.dd5.eq.'psi ')goto 251
	do j2=1,natom
	if(dd5.eq.'phi ') then
	if(dd1.eq.atname(j2).and.iresno(j2).eq.(i-1)) idd1=j2
	if(dd2.eq.atname(j2).and.iresno(j2).eq.i) idd2=j2
	if(dd3.eq.atname(j2).and.iresno(j2).eq.i) idd3=j2
	if(dd4.eq.atname(j2).and.iresno(j2).eq.i) idd4=j2
	elseif(dd5.eq.'psi ')then
	if(dd1.eq.atname(j2).and.iresno(j2).eq.i) idd1=j2
	if(dd2.eq.atname(j2).and.iresno(j2).eq.i) idd2=j2
	if(dd3.eq.atname(j2).and.iresno(j2).eq.i) idd3=j2
	if(dd4.eq.atname(j2).and.iresno(j2).eq.(i+1)) idd4=j2
	else
	if(dd1.eq.atname(j2).and.iresno(j2).eq.i) idd1=j2
	if(dd2.eq.atname(j2).and.iresno(j2).eq.i) idd2=j2
	if(dd3.eq.atname(j2).and.iresno(j2).eq.i) idd3=j2
	if(dd4.eq.atname(j2).and.iresno(j2).eq.i) idd4=j2
	endif
	if(idd1.ne.0.and.idd2.ne.0.and.idd3.ne.0.and.idd4.ne.0)goto 251
	enddo
        
251	nvar = nvar + 1
	if(icomment.eq.1)write(*,41) i,idd1, idd2, idd3, idd4,dd5
	ivres(nvar) = i
	vty(nvar) = dd5
	ivat(nvar,1) = idd1
	ivat(nvar,2) = idd2
	ivat(nvar,3) = idd3
	ivat(nvar,4) = idd4
	u0(nvar) = tu0
	sn(nvar) = tsn
	tn(nvar) = ttn
	enddo
	endif
	enddo
	enddo
	close(unit=10)
c----------measure the dihedral angles of initial model -----------------
	do i=1,nn
	dih = 0.0
	  i1 = ivat(i,1)
	  i2 = ivat(i,2)
	  i3 = ivat(i,3)
	  i4 = ivat(i,4)
	if(i1.ne.0.and.i2.ne.0.and.i3.ne.0.and.i3.ne.0) then
	dih = pdihedr(i1,i2,i3,i4)
	endif
	phi(i) = dih
	enddo
	do ii = 1,nn
	phi(ii) = (-1.0*phi(ii))+180.0
	enddo

        do k=1,natom
           do ki=1,3
             x_one(k,ki)=xpep(k,ki)
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
c
        enddo

        do k=1,natom
           do ki=1,3
             xpep(k,ki)=x_two(k,ki)
             ypep(k,ki)=x_two(k,ki)
           enddo
        enddo

c
c*******INITIALIZATION FOR DOCKING CALCULATIONS**********************
315	if(icomment.eq.1)print *,'docking initialization'
	if(iff.eq.1) call pdockinit
 	if(iff.eq.2) call pdockinit
c********************************************************************
        return
        end

c********************************************************************
	subroutine ppdb
	include 'mols.par'

        common /pepcrda/xpep(maxatm,8)
        common /par/natom,ntor,nhb,ns,lres
	common /pepdbat/pepatom(maxatm),pepele(maxatm)
        common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

	character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
	character*30 pepatom
        character*24 pepele

15      format(a30,3f8.3,a24)
	open(unit=26,file=pf1,status='old')
	do i=1,natom
        read(26,15)pepatom(i),(xpep(i,j),j=1,3),pepele(i)
	enddo
	close(unit=26)
	return
	end
c***********************************************************************
        subroutine pinput
	include 'mols.par'

	common /hbpar/mnhb
	common /vectors/iv(maxpar,4)
	common /order/nn,mm
        common /pepcrda/xpep(maxatm,8)
	common /native/nat,pepfile
	common /par/natom,ntor,nhb,ns,lres
	common /hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,25)
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

	character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

        open(unit=1,file=if2,status='old')
            
        do i=1,natom
          read(1,*) atn,(xpep(i,j),j=4,8)
          i1=ifix(xpep(i,7))
          i2=ifix(xpep(i,8))
          read(1,*)(jstart(i,j),jend(i,j),j=1,i1),(j1_4(i,k),k=1,i2)
        enddo

        do i=1,nn
          read(1,*)(iv(i,j),j=1,4)
        enddo
	if(iff.eq.1) then
        do i=1,mnhb
          read(1,*)ihb1(i),ihb2(i),c(i),d(i)
        enddo
	endif
        close(unit=1)


       	return
        end
c************************************************************************
	function pdihedr(i1,i2,i3,i4)
	include 'mols.par'

        common /pepcrda/xpep(maxatm,8)
        common /par/natom,ntor,nhb,ns,lres

        real x(maxatm,8)
        acdc=((180.0*7.0)/22.0)
	one=1.d0

        do i=1,natom
        do j=1,8
        x(i,j)=xpep(i,j)        
        enddo
        enddo

        x1=x(i2,1)-x(i1,1)
        y1=x(i2,2)-x(i1,2)
        z1=x(i2,3)-x(i1,3)
        x2=x(i3,1)-x(i2,1)
        y2=x(i3,2)-x(i2,2)
        z2=x(i3,3)-x(i2,3)
        ux1=y1*z2-z1*y2
        uy1=z1*x2-x1*z2
        uz1=x1*y2-y1*x2
        x1=x(i4,1)-x(i3,1)
        y1=x(i4,2)-x(i3,2)
        z1=x(i4,3)-x(i3,3)
        ux2=z1*y2-y1*z2
        uy2=x1*z2-z1*x2
        uz2=y1*x2-x1*y2

        u1=ux1*ux1+uy1*uy1+uz1*uz1
        u2=ux2*ux2+uy2*uy2+uz2*uz2
        u=u1*u2

        if(u.ne.zero) then
        a=(ux1*ux2+uy1*uy2+uz1*uz2)/sqrt(u)
        a=max(a,-one)
        a=min(a,one)
        pdihedr=acos(a)*acdc
        if (ux1*(uy2*z2-uz2*y2)+uy1*(uz2*x2-ux2*z2)+
     #      uz1*(ux2*y2-uy2*x2).lt.zero) pdihedr =-pdihedr
        return
      else
        write (*,'(a,4i5)')' pdihedr> Error in coordinates of atoms #: '
     #                     ,i1,i2,i3,i4
      stop
      endif
      end
c************************************************************************
c***************   PLP INITIALIZATIONS   ********************************

	subroutine pdockinit()
	include 'mols.par'

	common /rad/pbig
	common /comment/icomment
	common /patom/ipatom
	common /cen/bx,by,bz
	common /pplp/prresid(mnatp),matp
	common /pepdbat/pepatom(maxatm),pepele(maxatm)
	common /par/natom,ntor,nhb,ns,lres
	common /plp/pat(mnatp),lat(maxatm),patnam(mnatp),
     &presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),
     &natp,tatyp(25,50),ntatp(25),px(mnatp,3)

	character yy4*4,zz3*3,z3*3,ltyp*4,ptyp*4,xx*1,xx4*4,xx3*3,
     & xatnam*4, xresnam*3,lat*2,pat*2,tatyp*2,yy5*4,xchid*1
	character patnam*4,presnam*3,tpres*3,tatnam*4,pepatom*30
	character resnum*4,pchid*1
        character xxatnam*7
	real x1,x2,x3,d,a,b,c,y,ene,cx,cy,cz,cutdist
        integer prresid,xresid,ic,natp

	if(icomment.eq.1)print *,'Read receptor file, assign atom type 
     & to receptor'
        if(icomment.eq.1)print *,'atoms and ligand atoms'

10      format(a7,i4,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
12      format(5x,a3,2x,i2)
13      format(5x,a4,17x,a2)
14      format(a1)

c***********read atom records of protein and ligand******************
        if(icomment.eq.1)print *,'pvarinit-pbig',pbig
        natp=ipatom
	cx = 0.0
	cy = 0.0
	cz = 0.0
	cutdist = 0.0

        open(unit=1,file='rec.pdb',status='unknown')
        open(unit=2,file='prot.pdb',status='unknown')
	ic = 0
	do i=1,natp
        read(1,10) xxatnam,ixnum,xatnam,xresnam,xchid,xresid,x1,x2,x3
	if(xatnam(1:1).eq.'H')goto 100
	if(xatnam(2:2).eq.'H')goto 100
	cutdist =dist(bx,by,bz,x1,x2,x3)
	if(cutdist.gt.(5.0+6.0+pbig))goto 100
c	5.0 = box
c       6.0 = max. range of the selected interaction energy
c       pbig=half length of extended peptide
	ic = ic+1
	patnam(ic)=xatnam
	prresid(ic)=xresid
        presnam(ic)=xresnam
	pchid(ic)=xchid
	px(ic,1) = x1
	px(ic,2) = x2
	px(ic,3) = x3
        write(2,10) xxatnam,ic,xatnam,xresnam,xchid,xresid,x1,x2,x3
100	enddo

	natp = ic
        matp = ic
	close(unit=1)
c***********Read template file for protein atom types*************
	open(unit=3,file='PATOMTYPE', status='unknown')
	do i = 1,25
	read(3,12) tpres(i),ntatp(i)
	do j = 1,ntatp(i)
	read(3,13) tatnam(i,j),tatyp(i,j)
	enddo
	read(3,14) xx
	enddo
c***********Assign atom type for protein ******************
	if(icomment.eq.1)print *,'Protein atom and its type'
	do i = 1,natp
	xx4 = patnam(i)
	xx3 = presnam(i)
	do j1 = 1,25
	zz3 = tpres(j1)
	if(zz3.eq.xx3) then
	do j2 = 1,ntatp(j1)
	yy4 = tatnam(j1,j2)
	yy5 = tpres(j1)
	if(yy4(1:4).eq.xx4(1:4).and.yy5(1:3).eq.xx3(1:3))then
	pat(i) = tatyp(j1,j2)
	endif
	enddo
	endif
	enddo
	enddo
c***********Assign atom type for ligand ******************
	if(icomment.eq.1)print *,'Ligand atom and its type'
	if(icomment.eq.1)print *,'no. of lig atoms: ',natom
	do i=1,natom
	xx4 = pepatom(i)(13:16)
	xx3 = pepatom(i)(18:20)
	k = 0
	do j1 = 1,25
	zz3 = tpres(j1)
	if(zz3.eq.xx3) then
	do j2 = 1,ntatp(j1)
	yy4 = tatnam(j1,j2)
	if(yy4(1:4).eq.xx4(1:4)) then
	lat(i) = tatyp(j1,j2)
	k = 1
	endif
	enddo
	endif
	enddo
	if(k.eq.0) lat(i) = 'X'
	enddo
	return
	end
c***************************************************************
