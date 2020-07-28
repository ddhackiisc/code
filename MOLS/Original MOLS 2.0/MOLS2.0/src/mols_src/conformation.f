c     library name : conformation.f   

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




	subroutine conformation(kk)
c	program to write pdb file from cg output

	include 'mols.par'
	common/mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
        common/crda/x(maxatm,8) 
        common/crdb/y(maxatm,8) 
        common/vectors/iv(maxpar,4)
        common /par/ natom, ntor, nhb, ns, lres
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common/pdbat/atom(maxatm),ele(maxatm)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4,of0
	dimension pp(maxpar)
	character*10 pdbfile
	character*30 atom
        character*24 ele
15	format(a30,3f8.3,a24)
cc21	format(a5,4x,i4,3xf5.2)
21	format(a5,4x,i4,3x,f10.2)
22	format(a6)
23      format(a3)

	nn = ntor
c	call pdb
c        do i=1,natom
c          do j=1,8
c             y(i,j)=x(i,j)
c          enddo          
c	enddo
	if(kk.eq.1) then
	open(unit=50,file=pf2,status='unknown')
	elseif(kk.eq.2)then
	open(unit=50,file=pf3,status='unknown')
	elseif(kk.eq.0)then
	open(unit=50,file='molsx.pdb',status='unknown')
	goto 313
	endif
	do k11=1,ns
cc10	format(1x,19(1x,f5.1),2x,f15.4)
cc11	format(a10)
	if(kk.eq.1) then
	do j = 1,nn
	pp(j) = opmo(k11,j)
	enddo
	ent = emo(k11)
	elseif(kk.eq.2) then
	do j = 1,nn
	pp(j) = opmi(k11,j)
	enddo
	ent = emi(k11)
	endif
c	print *,'from conformation'
c	write(*,*) k11,ent,(opmo(k11,k),k=1,nn)

      call molgen3(pp)
c      print *,k11,ampene(0,0)

	write(50,21)'MODEL',k11,ent
313	do i=1,natom
	write (50,15) atom(i),y(i,1),y(i,2),y(i,3),ele(i)
	enddo
	write(50,23)'TER'
	write(50,22)'ENDMDL'
	if(kk.eq.0)return

	enddo
	close(unit=50)
	return
	end
c************************************************************************
       subroutine molgen3(phi)

	include 'mols.par'
        common/crda/x(maxatm,8) 
        common/crdb/y(maxatm,8) 
        common /par/ natom, ntor, nhb, ns, lres
        common/vectors/iv(maxpar,4)
        common/order/nn,mm
        dimension x_one(maxatm,3),x_two(maxatm,3)
        dimension phi(maxpar)
 	nn = ntor
 
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
C####################################################################      

        enddo

        do k=1,natom
           do ki=1,3
             y(k,ki)=x_two(k,ki) 
           enddo
        enddo

        return
        end

c********************************************************************
