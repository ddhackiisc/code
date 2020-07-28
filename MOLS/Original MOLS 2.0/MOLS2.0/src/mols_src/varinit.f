c     library name : varinit.f   

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





	subroutine varinit
	include 'mols.par'
	common /par/ natom, ntor, nhb, ns, lres
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common /atm/ iatmno(maxatm),atname(maxatm),rename(maxatm),
     &              iresno(maxatm)
	common /var/ ivat(maxres*8,4),vty(maxres*8),ivres(maxres*8)
        common/order/nn,mm
        common/vectors/iv(maxpar,4)
        common/crda/x(maxatm,8) 
        common/crdb/y(maxatm,8) 
	common/pdbat/atom(maxatm),ele(maxatm)
        common/ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,25)
        common/hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
       	common/ctl/iopt,iff,icint,fseq(maxres)

        dimension x_one(maxatm,3),x_two(maxatm,3)
  	external dihedr
  	integer iff
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*30 atom*30, ele*24
	character dd1*4,dd2*4,dd3*4,dd4*4,dd5*4,vty*4,tcode*4,scode*1,
     &  atname*4, rename*4,fseq*1
	nn = ntor
41	format(i5,2x,4i5,2x,a4)
112	format(8x,4a4,2x,a4,3f7.2)
90	format(a4,a1,1x,i2)
10	format(8x,i3,1x,a4,1x,a4,3x,i3)
c
	
	open(unit=1,file=pf1,status='old')
	do i=1,natom
		read(1,10) iatmno(i),atname(i),rename(i),iresno(i)
c		write(*,*) iatmno(i),atname(i),rename(i),iresno(i)
	enddo
	close(unit=1)
	call pdb
	call input
c-----------------define dihedral angles atoms-------------------------
cc	open(unit=21,file='molsx.inp',status='unknown')
	nvar = 0
	do i = 1,lres
	rewind 10
	open(unit=10,file='VAR.lib',status='old')
	do j=1,111
	read(10,90) tcode,scode,natm
c	print *,tcode,scode,natm,iopt,fseq(i)
	if(scode.eq.fseq(i))then
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
c	print *,dd1,dd2,dd3,dd4,dd5,tu0,tsn,ttn
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
c	write(*,41) i,idd1, idd2, idd3, idd4,dd5
cc	write(21,41) i, idd1, idd2, idd3, idd4, dd5
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
cc	close(unit=21)
c----------measure the dihedral angles of initial model -----------------
	do i=1,nn
	dih = 0.0
	  i1 = ivat(i,1)
	  i2 = ivat(i,2)
	  i3 = ivat(i,3)
	  i4 = ivat(i,4)
	if(i1.ne.0.and.i2.ne.0.and.i3.ne.0.and.i3.ne.0) then
	dih = dihedr(i1,i2,i3,i4)
	endif
	phi(i) = dih
c	  print *,i,i1,i2,i3,i4,dih
	enddo
c
c	write(31,*)(phi(ii),ii=1,nn)
c	write(*,*)(phi(ii),ii=1,nn)
c
	do ii = 1,nn
	phi(ii) = (-1.0*phi(ii))+180.0
	enddo
ccccc
c	do ii = 1,nn
c	phi(ii) = 0.0
c	enddo
cccc
c
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
             x(k,ki)=x_two(k,ki) 
             y(k,ki)=x_two(k,ki) 
           enddo
        enddo
c
15	format(a30,3f8.3,a24)
	open(unit=50,file='molsx.pdb',status='unknown')
	do i=1,natom
	write (50,15) atom(i),y(i,1),y(i,2),y(i,3),ele(i)
	enddo
	close(unit=50)
c	stop
	
        return
        end

c********************************************************************
	subroutine pdb
	include 'mols.par'
        common/crda/x(maxatm,8) 
        common /par/ natom, ntor, nhb, ns, lres
	common/pdbat/atom(maxatm),ele(maxatm)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4,of0
	character*30 atom
        character*24 ele
	open(unit=26,file=pf1,status='old')
c	open(unit=26,file='molso.pdb',status='old')
c	open(unit=26,file='scheraga.pdb',status='old')
15	format(a30,3f8.3,a24)
	do i=1,natom
        read(26,15)atom(i),(x(i,j),j=1,3),ele(i)
c        read(26,15)atom(i),xt,yt,zt,ele(i)
	enddo
	close(unit=26)
	return
	end
c***********************************************************************
        subroutine input

	include 'mols.par'
        common/crda/x(maxatm,8) 
        common/ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,25)
        common /par/ natom, ntor, nhb, ns, lres
        common/vectors/iv(maxpar,4)
        common/order/nn,mm
        common/hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
       	common/ctl/iopt,iff,icint,fseq(maxres)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4,of0
      	character fseq*1


        open(unit=1,file=if1,status='old')
c        open(unit=1,file='gly.inp',status='old')

        do i=1,natom
          read(1,*) atn,(x(i,j),j=4,8)
          i1=ifix(x(i,7))
          i2=ifix(x(i,8))
          read(1,*)(jstart(i,j),jend(i,j),j=1,i1),(j1_4(i,k),k=1,i2)
        enddo

        do i=1,nn
          read(1,*)(iv(i,j),j=1,4)
          print *,(iv(i,j),j=1,4)
        enddo
	if(iff.eq.1) then
        do i=1,nhb
          read(1,*)ihb1(i),ihb2(i),c(i),d(i)
        enddo
	endif
        close(unit=1)
	
       	return
        end
c*******************************************************************
c************************************************************************	
       function dihedr(i1,i2,i3,i4)
	include 'mols.par'
        common/crda/x(maxatm,8)
        acdc=((180.0*7.0)/22.0)
	one=1.d0
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

      if (u.ne.zero) then
        a=(ux1*ux2+uy1*uy2+uz1*uz2)/sqrt(u)
        a=max(a,-one)
        a=min(a,one)
        dihedr=acos(a)*acdc
        if (ux1*(uy2*z2-uz2*y2)+uy1*(uz2*x2-ux2*z2)+
     #      uz1*(ux2*y2-uy2*x2).lt.zero) dihedr =-dihedr
        return
      else
        write (*,'(a,4i5)')' dihedr> Error in coordinates of atoms #: '
     #                     ,i1,i2,i3,i4
        stop
      endif
      end
c****************************************************************************
