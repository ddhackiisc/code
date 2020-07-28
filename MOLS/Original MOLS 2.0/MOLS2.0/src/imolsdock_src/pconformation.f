c-----------------------------------------------------------------------------
c     library name : pconformation.f   

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
        subroutine pconformation(inp,kk)
c	program to write pdb file

	include 'mols.par'

	common /comment/icomment
	common /pepcrda/xpep(maxatm,8)
        common /pepcrdb/ypep(maxatm,8)
	common /vectors/iv(maxpar,4)
	common /conf/path,mname
	common /par/natom,ntor,nhb,ns,lres,is,ie
	common /pepdbat/pepatom(maxatm),pepele(maxatm)
	common /pout/popmi(maxstr,maxpar),pemi(maxstr)
	common /protene/penb,pees,penb1_4,pees1_4
	common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /ambecep/ees,enb,ees1_4,ehb,enb1_4,etor
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
	common /autoenrg/less,lnb,lhb,pess,pnb,phb,intra,inter
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
	common /pepcrdbf/ypepf(maxatm,8),atomf(maxatm),pepelef(maxatm)
	common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),emo(maxstr),
     &  emi(maxstr)
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
	common /plp/pat(mnatp),lat(maxatm),patnam(mnatp),
     &  presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),natp,
     &  tatyp(25,50),ntatp(25),px(mnatp,3)
        common /propdb/ihflex,patname(maxatm),px1(maxatm),px2(maxatm),
     &  px3(maxatm),presname(maxatm),pcid(maxatm),presnum(maxatm)

	character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6,temp,path,mname,pepmols,pepmini,plpmols,plpmini,
     &  promols,promini

	dimension pp(maxpar)
	character*30 atom,pepatom,atomf
	character*24 pepele,pepelef
	character*10 pdbfile

	real less,lnb,lhb,pess,pnb,phb,intra,inter

15	format(a30,3f8.3,a24)
21	format(a5,4x,i4,3x,f10.2)
22	format(a6)
23      format(a3)
40	format(a4,6a10)
41	format(a4,5a10)
42      format(a4,4a10)
43	format(i4,6f10.3)
44      format(i4,5f10.3)
45      format(i4,4f10.3)

	nn = ntor
	pepmols = path(:LNBLNK(path)) // '/' //
     &	mname(:LNBLNK(mname)) // '_' // 'pepmols.out'
	pepmini = path(:LNBLNK(path)) // '/' //
     &	mname(:LNBLNK(mname)) // '_' // 'pepmini.out'
        plpmols = path(:LNBLNK(path)) // '/' //
     &  mname(:LNBLNK(mname)) // '_' // 'plpmols.out'
	plpmini = path(:LNBLNK(path)) // '/' //
     &  mname(:LNBLNK(mname)) // '_' // 'plpmini.out'
        promols = path(:LNBLNK(path)) // '/' //
     &  mname(:LNBLNK(mname)) // '_' // 'promols.out'
        promini = path(:LNBLNK(path)) // '/' //
     &  mname(:LNBLNK(mname)) // '_' // 'promini.out'

        if(icomment.eq.1)print *,'pconf-natom',natom

	open(unit=80,file=pepmols,status='unknown',access='append')
	open(unit=81,file=pepmini,status='unknown',access='append')
	open(unit=82,file=plpmols,status='unknown',access='append')
	open(unit=83,file=plpmini,status='unknown',access='append')
	open(unit=84,file=promols,status='unknown',access='append')
	open(unit=85,file=promini,status='unknown',access='append')

	IF(inp.eq.is)THEN
	if(kk.eq.1) then
	if(iff.eq.1) write(80,40) 'SNo','vdw','Estatic','HB','NB1-4',
     &  'ES1-4','AMBER'
	if(iff.eq.2) write(80,41) 'SNo','Estatic','vdw','HB','Tor','ECEPP'
     	else if(kk.eq.2)then
	if(iff.eq.1) write(81,40) 'SNo','vdw','Estatic','HB','NB1-4',
     &  'ES1-4','AMBER'
	if(iff.eq.2) write(80,41) 'SNo','Estatic','vdw','HB','Tor','ECEPP'
     	endif

	if(kk.eq.1)then
	write(82,42)'SNo','HB','STERIC','REPULSIVE','PLP'
	else if(kk.eq.2)then
	write(83,42)'SNo','HB','STERIC','REPULSIVE','PLP'
	endif

        if(kk.eq.1) then
        write(84,41) 'SNo','vdw','Estatic','NB1-4','ES1-4','AMBER'
        else if(kk.eq.2)then
        write(85,41) 'SNo','vdw','Estatic','NB1-4','ES1-4','AMBER'
        endif
	ENDIF

	if(kk.eq.1) then
	open(unit=50,file=pf6,status='unknown',access='append')
	else if(kk.eq.2)then
	open(unit=50,file=pf7,status='unknown',access='append')
	endif
        
	do k11=inp,inp

	if(kk.eq.1) then

	do j = 1,nn+6+np
	pp(j) = opmo(k11,j)
	phi(j) = opmo(k11,j)
	enddo

	do j = nn+1,nn+3
	tt1 = opmo(k11,j)
	if(tt1.eq.0.0) then
	tt2 = 0.0
	else
	tt2 = (tt1/10.0)*0.13888889
	endif
	tt2 = tt2-2.5
	pp(j) = tt2
	enddo
	ent = emo(k11)

	else if(kk.eq.2) then

	do j = 1,nn+6+np
	pp(j) = opmi(k11,j)
	phi(j) = opmi(k11,j)
	enddo

	do j = nn+1,nn+3
	tt1 = opmi(k11,j)
	if(tt1.eq.0.0) then
	tt2 = 0.0
	else
	tt2 = (tt1/10.0)*0.13888889
	endif
	tt2 = tt2-2.5
	pp(j) = tt2
	enddo
	ent = emi(k11)
	endif
 
	if(icomment.eq.1)then
        do j=1,nn+6+np
        print *,'pconf-phi',pp(j)
        enddo
	endif

        call pmolgen3(pp,k11)

	if(ilopt.eq.1)then
        call peplp(plpe,hb,steric,rep,1)
        else if(ilopt.eq.2)then
        call eplp(plpe,hb,steric,rep,1)
        endif


	write(50,21)'MODEL',k11,ent
        if(icomment.eq.1)print *,'natom-pconf',natom
	do i=1,natom
	write (50,15) pepatom(i),ypep(i,1),ypep(i,2),
     &  ypep(i,3),pepele(i)
        enddo
c----------------------------------------------------------------

	if(iff.eq.1) energy = pampene(k11,k11)
	if(iff.eq.2) energy = pecpene(k11,k11)
	

	if(kk.eq.1) then

	if(iff.eq.1) then
	write (80,43)inp,enb,ees,ehb,enb1_4,ees1_4,energy
	write (82,45)inp,hb,steric,rep,plpe
	write (84,44)inp,penb,pees,penb1_4,pees1_4,penb+pees+penb1_4+pees1_4
	else if(iff.eq.2) then
	write (80,44)inp,ees,enb,ehb,etor,energy
	write (82,45)inp,hb,steric,rep,plpe
	write (84,44)inp,penb,pees,penb1_4,pees1_4,penb+pees+penb1_4+pees1_4
	endif

     	else if(kk.eq.2)then

        if(iff.eq.1) then
        write (81,43)inp,enb,ees,ehb,enb1_4,ees1_4,energy
        write (83,45)inp,hb,steric,rep,plpe
	write (85,44)inp,penb,pees,penb1_4,pees1_4,penb+pees+penb1_4+pees1_4
        else if(iff.eq.2) then
        write (81,44)inp,ees,enb,ehb,etor,energy
        write (82,45)inp,hb,steric,rep,plpe
	write (85,44)inp,penb,pees,penb1_4,pees1_4,penb+pees+penb1_4+pees1_4
        endif

     	endif

	write(50,23)'TER'
	write(50,22)'ENDMDL'
	if(kk.eq.0)return

	enddo
	close(unit=50)
	close(unit=80)
	close(unit=81)
	close(unit=82)
	close(unit=83)
	close(unit=84)
	close(unit=85)
	return
	end
c************************************************************************
       subroutine pmolgen3(phi,k11)

	include 'mols.par'
	common /comment/icomment
        common /pepcrda/xpep(maxatm,8)
        common /pepcrdb/ypep(maxatm,8)
	common /vectors/iv(maxpar,4)
	common /order/nn,mm
	common /cen/bx,by,bz
        common /par/natom,ntor,nhb,ns,lres
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt

        dimension x_one(maxatm,3),x_two(maxatm,3)
        dimension phi(maxpar)

	real bx,by,bz

 	nn = ntor

        if(icomment.eq.1)print *,'natom',natom
	do k=1,natom
           do ki=1,3
             x_one(k,ki)=xpep(k,ki)
           enddo
        enddo

        if(icomment.eq.1)then
	do k=1,nn+6
	print *,'phi-pmolgen',k,phi(k)
	enddo
        endif

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
             ypep(k,ki)=x_two(k,ki)
           enddo
        enddo

        if(icomment.eq.1)print *,'molgen-pconf-natom',natom

	call protate(bx,by,bz,phi(nn+4),phi(nn+5),phi(nn+6),natom)

	call ptranslate(bx,by,bz,phi(nn+1),phi(nn+2),phi(nn+3),natom)

        if(ifopt.eq.2)then

        call flexall(2,2,phi,e_flex,2)

        open(unit=121,file='rflex.txt',status='unknown')
        write(121,*)e_flex
        close(unit=121)
        endif

	return
        end
c----------------------------------------------------------------------
