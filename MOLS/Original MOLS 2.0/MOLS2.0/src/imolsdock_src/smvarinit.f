c-----------------------------------------------------------------------------
c     library name : smvarinit.f   

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

	subroutine varinit
	include 'mols.par'

	common /order/nn,mm
        common /vectors/iv(maxpar,4)
        common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
	common /cen/bx,by,bz
	common /par/natom,ntor,nhb,ns,lres
	common /pdbat/atom(maxatm),proele(maxatm)
	common /gen/ibno,iat1(maxatm),iat2(maxatm),isno(maxatm)
	common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /var/ivat(maxres*8,4),vty(maxres*8),ivres(maxres*8)
        common /hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        common /atm/iatmno(maxatm),atname(maxatm),rename(maxatm),
     &              iresno(maxatm)
	common /plp/pat(mnatp), attyp(maxatm),patnam(mnatp),
     &  presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),
     &  natp,tatyp(25,50),ntatp(25),px(mnatp,3)
        common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &   iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),
     &   ind(maxatm,maxatm),big
        common /left/le(maxatm,maxatm),innd(maxatm),ielenum(maxatm),
     &  bonum(maxatm)

        character*50 desc*30,attyp*10,atsym*6,str1*30,str2*30,
     &  str3*30,res*10

        dimension x_one(maxatm,3),x_two(maxatm,3)
  
  	external dihedr
  	integer iff,n1,n2,n3

        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
	character*30 atom*30, proele*24
	character dd1*4,dd2*4,dd3*4,dd4*4,dd5*4,vty*4,tcode*4,scode*1,
     &  atname*4, rename*4,fseq*1
     	character xatnam*4, xresnam*3,lat*2,pat*2,tatyp*2
	character patnam*4,presnam*3,tpres*3,tatnam*4
	nn = ntor

c-----------------define dihedral angles atoms-------------------------
c----------measure the dihedral angles of initial model ---------------

	do ii = 1,nn
	phi(ii) = (-1.0*phi(ii))+180.0
	enddo
        do k=1,natom
           do ki=1,3
             x_one(k,ki)=rx(k,ki) 
           enddo
        enddo
        do if=1,nn

C###################### PHI ALL ####################################      

        call elemen(x_one(irotb(if,1),1),x_one(irotb(if,1),2),
     &              x_one(irotb(if,1),3),
     &              x_one(irotb(if,2),1),x_one(irotb(if,2),2),
     &              x_one(irotb(if,2),3),
     &              el,em,en)
50      format(i4,1x,3f8.3)
        do j=1,ielenum(if)
        k=le(if,j)
           do ki=1,3
           x_two(k,ki)=x_one(k,ki)
           enddo
        enddo
c******************************************************************

         do j=1,iatrot(if)
          k=ilsrot(if,j)

           xin=x_one(k,1)-x_one(irotb(if,1),1)
           yin=x_one(k,2)-x_one(irotb(if,1),2)
           zin=x_one(k,3)-x_one(irotb(if,1),3)
           call rotor(el,em,en,phi(if),xin,yin,zin,
     &                xout,yout,zout)
           x_two(k,1)=xout+x_one(irotb(if,1),1)
           x_two(k,2)=yout+x_one(irotb(if,1),2)
           x_two(k,3)=zout+x_one(irotb(if,1),3)

        enddo

        do k=1,natom
           do ki=1,3
              x_one(k,ki)=x_two(k,ki)
           enddo
        enddo
C##################################################################
        enddo

        do k=1,natom
           do ki=1,3
             xlig(k,ki)=x_two(k,ki) 
             ylig(k,ki)=x_two(k,ki) 
           enddo
        enddo

	call dockinit
        return
        end
c********************************************************************
	subroutine pdb
	include 'mols.par'
        common /procrda/xpro(maxatm,8)
        common /part/ipart,iptor
        common /par/natom,ntor,nhb,ns,lres
	common /pdbat/atom(maxatm),proele(maxatm)
        common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
	character*30 atom
        character*24 proele

15      format(a30,3f8.3,a24)
	open(unit=26,file='partners',status='old')
	do i=1,ipart
        read(26,15)atom(i),(xpro(i,j),j=1,3),proele(i)
	enddo
	close(unit=26)
	return
	end
c***********************************************************************
        subroutine input()

	include 'mols.par'
        common /procrda/xpro(maxatm,8)
	common /provectors/ivp(maxpar,4)
        common /order/nn,mm
        common /part/ipart,iptor
        common /par/natom,ntor,nhb,ns,lres
	common /recsc/iline,idihed,ivx(maxatm,4)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /prohb/iprohb1(maxhb),iprohb2(maxhb),cpro(maxhb),
     &  dpro(maxhb)
	common /proranges/jstartpro(maxatm,10),jendpro(maxatm,10),
     &  j1_4pro(maxatm,25)
        common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

      	character fseq*1

        close(unit=2)


        open(unit=1,file=if1,status='old')

        do i=1,ipart
          read(1,*) atn,(xpro(i,j),j=4,8)
          i1=ifix(xpro(i,7))
          i2=ifix(xpro(i,8))
          read(1,*)(jstartpro(i,j),jendpro(i,j),j=1,i1),
     & (j1_4pro(i,k),k=1,i2)

        enddo
        
                
        do i=1,iptor
         read(1,*)(ivp(i,j),j=1,4)
        enddo

	if(iff.eq.1) then

        do i=1,nhb
          read(1,*)iprohb1(i),iprohb2(i),cpro(i),dpro(i)
        enddo
	endif
        close(unit=1)

       	return
        end
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
	subroutine dockinit()
	include 'mols.par'

	common /comment/icomment
	common /patom/ipatom
	common /cen/bx,by,bz
	common /pplp/prresid(mnatp),matp
	common /pdbat/atom(maxatm),proele(maxatm)
	common /par/natom,ntor,nhb,ns,lres
        common /strin/tem1(maxatm),tem2(maxatm),slig_str3(maxatm),res,
     &  ato(maxatm)
        common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &  iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),
     &  ind(maxatm,maxatm),big
	common /plp/pat(mnatp),lat(maxatm),patnam(mnatp),
     &  presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),
     &  natp,tatyp(25,50),ntatp(25),px(mnatp,3)
        
        character yy4*4,zz3*3,z3*3,ltyp*4,ptyp*4,xx*1,xx4*4,xx3*3,
     &  xatnam*4,xresnam*3,lat*2,pat*2,tatyp*2,yy5*4,xchid*1,
     &  patnam*4,presnam*3,tpres*3,tatnam*4,atom*30,resnum*4,pchid*1,
     &  tem1*30,tem2*30,slig_str3*30,res*10,ato*80,xxatnam*7
	real x1,x2,x3,d,a,b,c,y,ene,cx,cy,cz
	real cutdist
        integer xresid,prresid

10      format(a7,i4,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
12      format(5x,a3,2x,i2)
13      format(5x,a4,17x,a2)
14      format(a1)
22      format(a80)
33      format(a25)
55      format(a30)
441     format(i4,i4,i2,i2,i2)

c***********read atom records of protein and ligand******************
	
	natp=ipatom
        if(icomment.eq.1)print*,'smvarinit-big',big
	cx = 0.0
	cy = 0.0
	cz = 0.0
	cutdist = 0.0

        if(icomment.eq.1)print *,'scan-dis',5.0+6.0+big
	open(unit=1,file='rec.pdb',status='unknown')
        open(unit=2,file='prot.pdb',status='unknown')
	ic = 0
	do i=1,natp
        read(1,10) xxatnam,ixnum,xatnam,xresnam,xchid,xresid,x1,x2,x3
	if(xatnam(1:1).eq.'H')goto 100
	if(xatnam(2:2).eq.'H')goto 100
	cutdist =dist(bx,by,bz,x1,x2,x3)
	if(cutdist.gt.(5.0+6.0+big))goto 100
	ic = ic+1
	patnam(ic)=xatnam
	presnam(ic)=xresnam
        prresid(ic)=xresid
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
!*********************************************************************
        open(unit=7,file='tem.mol2',status='old')
        open(unit=8,file='atomtype.dat',status='old')

        do l=1,31
        read(8,55),slig_str3(l)
        enddo

        do i=1,2
        read(7,33)tem1(i)
        enddo

        read(7,441)nnatom,ibno,n1,n2,n3

        do j=1,4
        read(7,33)tem2(j)
        enddo

        do  k=1,natom
        read(7,22),ato(k)
         do  n=1,31
          if( ato(k)(48:55) .eq. slig_str3(n)(6:13) ) then
           lat(k)= slig_str3(n)(23:24)
          endif
         enddo
        enddo

        close(7)
        close(8)
	return
        end
!***********************************************************************
