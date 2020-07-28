c-----------------------------------------------------------------------------
c     library name : findflex.f   

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

	subroutine findflex()
	include 'mols.par'

	common /comment/icomment
	common /rad/pbig
	common /find/fdis
	common /patom/ipatom
	common /cen/bx,by,bz
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
	common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &  iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),
     &  ind(maxatm,maxatm),big

	character*6 atom(mnatp) 
	character*4 attyp(mnatp),atyp
	character*3 res(mnatp)
	character*1 cid(mnatp),csres(mnatp),cres(mnatp),chid
	real x1(mnatp),x2(mnatp),x3(mnatp),xx,yy,zz,cutd,
     &  rcut
	integer atno(mnatp),resno(mnatp),pres,sres(mnatp),
     &  ssres(mnatp),temp

	if(ilopt.eq.1)then !Peptide Docking
	if(icomment.eq.1)print *,'cen',bx,by,bz
	if(icomment.eq.1)print *,'pbig',pbig,'fdis',fdis
	rcut=pbig+fdis
	else if(ilopt.eq.2)then!small molecule
	if(icomment.eq.1)print *,'cen',bx,by,bz
        if(icomment.eq.1)print *,'big',big,'fdis',fdis
	rcut=big+fdis
	endif
c-------read rec.pdb-------------------------------------
10	format(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
20	format(a3,1x,a1,1x,i4)
	open(unit=1,file='rec.pdb',status='old')
	do i=1,ipatom
	read(1,10)atom(i),atno(i),attyp(i),res(i),cid(i),
     &  resno(i),x1(i),x2(i),x3(i)
	enddo
	close(1)
c--------------------------------------------------------
c------find distance between grid centre and protein-----
	jres=0
	do i=1,ipatom
	xx=x1(i)
	yy=x2(i)
	zz=x3(i)
	pres=resno(i)
        chid=cid(i)
	atyp=attyp(i)
	cutd=dist(bx,by,bz,xx,yy,zz)
	if(cutd.le.rcut.and.atyp.ne.' N  '.and.atyp.ne.' CA '.and.
     &  atyp.ne.' C  '.and.atyp.ne.' O  ')then
	jres=jres+1
	sres(jres)=pres
        cres(jres)=chid
	endif
	enddo
c------------------------------------------------------------
c-------picking id of partner residues begins here-----------
        l=0

        do i=1,jres !total number of atoms
        ik=0
        if(i.eq.1)then
         l=l+1
         ssres(l)=sres(i)!first residue id written into ssres(1)
         csres(l)=cres(i)
        else
          do j=1,l
           if(sres(i).ne.ssres(j).and.cres(i).ne.csres(i))then
            ik=ik+1
           endif
          enddo

          if(ik.eq.l)then
            l=l+1
            ssres(l)=sres(i)
            csres(l)=cres(i)
           if(icomment.eq.1) print *,'ssres',ssres(l),csres(l)
          endif
        endif
        enddo

       if(icomment.eq.1)print *,'flexible residues',l
c---------------------------------------------------------------
c------sort the partner residues in ascending order-------------
        do i=1,l-1
         do j=i+1,l
          if(ssres(i)>ssres(j))then
           temp=ssres(i)
           ssres(i)=ssres(j)
           ssres(j)=temp
          endif
         enddo
        enddo
        open(unit=110,file='sortfres',status='unknown')
        do i=1,l
         write(110,*)ssres(i)
        enddo
        close(110)
c----------------------------------------------------------------
	open(unit=2,file='residues.txt',status='unknown')
	do il=1,l
	im=ssres(il)
	do ik=1,ipatom
	imp=resno(ik)
	if(imp.eq.im.and.cid(ik).eq.csres(il))then
	write(2,20)res(ik),cid(ik),resno(ik)
	goto 100
	endif
	enddo
100	enddo
	close(2)
c----------------------------------------------------------------
	return
	end
