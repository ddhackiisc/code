c-----------------------------------------------------------------------------
c     library name : conformation.f   

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
	subroutine conformation(kk,inp,kk1)

        include 'mols.par'

        common /comment/icomment
        common /patom/ipatom
        common /profile/pfile
        common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /pepcrdb/ypep(maxatm,8)
        common /vectors/iv(maxpar,4)
        common /conf/path,mname
        common /cen/bx,by,bz
        common /pdbat/atom(maxatm),proele(maxatm)
        common /pepdbat/pepatom(maxatm),pepele(maxatm)
        common /par/natom,ntor,nhb,ns,lres
        common /gen/ibno,iat1(maxatm),iat2(maxatm),isno(maxatm)
        common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &  emo(maxstr),emi(maxstr)
        common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        common /recep/np,nres,crnam(mnatp),crcid(mnatp),irrid(mnatp)
        common /pepcrdbf/ypepf(maxatm,8),atomf(maxatm),pepelef(maxatm)
        common /left/le(maxatm,maxatm),innd(maxatm),ielenum(maxatm),
     &  bonum(maxatm)
        common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &  iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),ind(maxatm,maxatm),
     &  big
        common /string/atsym(maxatm),desc(maxatm),attyp(maxatm),
     &  str1(maxatm),str2(maxatm),str3,res,n1,n2,n3
        common /plp/pat(mnatp),lat(maxatm),patnam(mnatp),
     &  presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),natp,
     &  tatyp(25,50),ntatp(25),px(mnatp,3)
        common /propdb/ihflex,patname(maxatm),px1(maxatm),px2(maxatm),
     &  px3(maxatm),presname(maxatm),pcid(maxatm),presnum(maxatm)

        dimension pp(maxpar)

        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        character*80 prt,prt1
        character*30 desc,str1,str2,str3,atom,pepatom
        character*24 proele,pepele
        character*12 prt2(mnatp)
        character*10 pdbfile,attyp,res
        character*6 atsym
        character*4 qatname(mnatp),xx3,yy3,patname
        character*3 qresname(mnatp),xx1,yy1,presname
        character*1 pcid,qcid(mnatp),xx2,yy2

        integer npar,nap,exresid,zero,presnum,qresnum(mnatp)
        real ex1,ex2,ex3,mine,ene_flex,qx1(mnatp),qx2(mnatp),
     &  qx3(mnatp),px1,px2,px3

11      format(a25)
12      format(a17)
13      format(a13)
15      format(a30,3f8.3,a24)
21      format(a5,4x,i4,3x,f10.2)
22      format(a6)
23      format(a3)
28      format(i4,i4,i2,i2,i2)
30      format(3x,i4,1x,a6,2x,3f10.4,1x,a29)
60      format(1x,i5,1x,i5,1x,i5,a10)
130     format(a12,a4,1x,a3,1x,a1,i4,4x,3f8.3,24x)
170     format(a23,i3,4x,3f8.3,a24)
180     format(a6,i5,2x,a4,a3,5x,i1,4x,3f8.3)

        nn = ntor
        npar=nn+6+np
        nap=ipatom
        zero=0

	IF(kk1.eq.1)THEN
c-------Peptide-Protein Docking---------------------------
	nap=ipatom
        if(icomment.eq.1)print *,'nap',nap,ipatom

        do k11=inp,inp!1,ns

         if(ifopt.eq.1)then

          do j = 1,nn+6
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
         endif

        if(ifopt.eq.2)then

        do j = nn+7,npar
        pp(j) = opmo(k11,j)
        phi(j) = opmo(k11,j)
        enddo

        ent = emo(k11)
        endif

        call molgen3(pp,ifopt,k11)

c-------writing MOLS complex file-----------------------------
	open(unit=11,file=pf6,status='old')
        open(unit=12,file='rec.pdb',status='old')
        open(unit=14,file=of6,access='append')

        if(ifopt.eq.1)then

         write(14,21)'MODEL',inp,ent
         do i=1,nap
         read(12,170)prt,exresid,ex1,ex2,ex3,prt1
         write(14,170)prt,exresid,ex1,ex2,ex3,prt1
         enddo
         write(14,23)'TER'

         do i=1,natom
         write(14,15)pepatom(i),ypep(i,1),ypep(i,2),
     &   ypep(i,3),pepele(i)
         enddo

         write(14,23)'TER'
         write(14,22)'ENDMDL'

        else if(ifopt.eq.2)then

         call flexall(2,2,phi,e_flex,2)
 
        write(14,21)'MODEL',inp,ent

 	do i=1,nap

         read(12,130)prt2(i),qatname(i),qresname(i),qcid(i),qresnum(i),
     &   qx1(i),qx2(i),qx3(i)

         do j=1,ihflex
         xx1=presname(j)
         xx2=pcid(j)
         xx3=patname(j)
         yy1=qresname(i)
         yy2=qcid(i)
         yy3=qatname(i)

         if(yy1(1:3).eq.xx1(1:3).and.yy2(1:1).eq.xx2(1:1)
     &   .and.yy3(1:4).eq.xx3(1:4)
     &   .and.presnum(j).eq.qresnum(i))then
         qx1(i)=px1(j)
         qx2(i)=px2(j)
         qx3(i)=px3(j)
         endif
         enddo
        enddo

	do i=1,nap
        write(14,130)prt2(i),qatname(i),qresname(i),qcid(i),qresnum(i),
     &   qx1(i),qx2(i),qx3(i)
        enddo
        write(14,23)'TER'

        do i=1,natom
         write(14,15)pepatom(i),ypep(i,1),ypep(i,2),
     &   ypep(i,3),pepele(i)
        enddo

        write(14,23)'TER'
        write(14,22)'ENDMDL'
        endif
        close(14)
        close(12)
        close(11)

        enddo!do k11=inp,inp for peptide-protein

c-------Protein-Ligand Docking-------------------------
	ELSE IF(kk1.eq.2)THEN

	if(kk.eq.1) then
         open(unit=50,file=pf2,access='append')
        elseif(kk.eq.2)then
         open(unit=50,file=pf3,access='append')
        endif

        do k11=inp,inp!1,ns
         if(ifopt.eq.1)then
          do j = 1,nn+6
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
        endif

        if(ifopt.eq.2)then
         do j = 1,nn+6
          pp(j) = opmo(k11,j)
          phi(j) = opmo(k11,j)
          if(icomment.eq.1)then
          print *,'fconformation,ifopt,pp',ifopt,j,pp(j)
          endif
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
         if(icomment.eq.1)print *,'fconformation,ifopt,pp',ifopt,j,pp(j)
         enddo

        do j = nn+7,npar
         pp(j) = opmo(k11,j)
         phi(j) = opmo(k11,j)
         if(icomment.eq.1)print *,'fconformation,ifopt,pp',ifopt,j,pp(j)
        enddo

        ent = emo(k11)

        endif

        call molgen3(pp,ifopt,k11)

        open(unit=551,file="mols.mol2",status="unknown")

        write(50,21)'MODEL',k11,ent
        write(551,21)'MODEL',k11,ent
        do i=1,2
         write(50,12)str1(i)
         write(551,12)str1(i)
        enddo
         write(50,28)natom,ibno,n1,n2,n3
         write(551,28)natom,ibno,n1,n2,n3
        do i=1,4
         write(50,13)str2(i)
         write(551,13)str2(i)
        enddo

313     do i=1,natom
         write(50,30),innd(i),atsym(i),ylig(i,1),ylig(i,2),
     &   ylig(i,3),desc(i)
         write(551,30),innd(i),atsym(i),ylig(i,1),ylig(i,2),
     &   ylig(i,3),desc(i)
        enddo
        write(50,11)str3
        write(551,11)str3
        do j=1,ibno
         write(50,60)isno(j),iat1(j),iat2(j),attyp(j)
         write(551,60)isno(j),iat1(j),iat2(j),attyp(j)
        enddo
        write(50,23)'TER'
        write(551,23)'TER'
        write(50,22)'ENDMDL'
        write(551,22)'ENDMDL'
        close(unit=551)

c--------sam added--writing MOLS complex file------------------------
        if(ifopt.eq.1)then
         open(unit=11,file='rec.pdb',status='old')
         open(unit=12,file=of6,access='append')
c--------------------------------------------------------------------
         write(12,21)'MODEL',k11,ent
         do i=1,nap
          read(11,170)prt,exresid,ex1,ex2,ex3,prt1
          write(12,170)prt,exresid,ex1,ex2,ex3,prt1
         enddo
         write(12,23)'TER'
         do i=1,natom
          write(12,180),'HETATM',nap+i+1,atsym(i),'LIG',zero,ylig(i,1),
     &    ylig(i,2),ylig(i,3)
         enddo
         write(12,23)'TER'
         write(12,22)'ENDMDL'
         close(12)
         close(11)

        elseif(ifopt.eq.2)then

         open(unit=11,file='rec.pdb',status='old')
         open(unit=12,file=of6,access='append')
         write(12,21)'MODEL',k11,ent

         call flexall(2,2,phi,e_flex,2)

         do i=1,nap

         read(11,130)prt2(i),qatname(i),qresname(i),qcid(i),qresnum(i),
     &   qx1(i),qx2(i),qx3(i)

         do j=1,ihflex
         xx1=presname(j)
         xx2=pcid(j)
         xx3=patname(j)
         yy1=qresname(i)
         yy2=qcid(i)
         yy3=qatname(i)
          if(yy1(1:3).eq.xx1(1:3).and.yy2(1:1).eq.xx2(1:1)
     &    .and.yy3(1:4).eq.xx3(1:4)
     &    .and.presnum(j).eq.qresnum(i))then
           qx1(i)=px1(j)
           qx2(i)=px2(j)
           qx3(i)=px3(j)
          endif
         enddo
        enddo

        do i=1,nap
         write(12,130)prt2(i),qatname(i),qresname(i),qcid(i),qresnum(i),
     &  qx1(i),qx2(i),qx3(i)
        enddo

         write(12,23)'TER'

         do i=1,natom
          write(12,180),'HETATM',nap+i+1,atsym(i),'LIG',zero,ylig(i,1),
     &    ylig(i,2),ylig(i,3)
         enddo
         write(12,23)'TER'
         write(12,22)'ENDMDL'
         close(12)
         close(11)

         endif

	enddo !do k11=inp,inp Protein-Ligand

	ENDIF

        return
	end
c----------------------------------------------------------------------
        subroutine molgen3(phi,ifopt,k11)

        include 'mols.par'

	common /comment/icomment
	common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
	common /vectors/iv(maxpar,4)
	common /order/nn,mm
	common /cen/bx,by,bz
	common /par/natom,ntor,nhb,ns,lres
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &   iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),
     &   ind(maxatm,maxatm),big
        common /left/le(maxatm,maxatm),innd(maxatm),ielenum(maxatm),
     &  bonum(maxatm)

        dimension x_one(maxatm,3),x_two(maxatm,3)
        dimension phi(maxpar)

        real bx,by,bz

        nn=ntor
        do k=1,natom
           do ki=1,3
             x_one(k,ki)=rx(k,ki)
           enddo
        enddo
c-------------------------------------------------------------------
        if(ntor.eq.0)then
        do k=1,natom
           do ki=1,3
             ylig(k,ki)=rx(k,ki)
           enddo
        enddo
        goto 333
        endif
c-------------------------------------------------------------------

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

c------------------------------------------------------

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
C####################################################################

        enddo

        do k=1,natom
           do ki=1,3
             ylig(k,ki)=x_two(k,ki)
           enddo
        enddo

333     call rotate(bx,by,bz,phi(nn+4),phi(nn+5),phi(nn+6),natom)
        call translate(bx,by,bz,phi(nn+1),phi(nn+2),phi(nn+3),natom)

        if(ifopt.eq.2)then
        call flexall(2,2,phi,e_flex,2)
        open(unit=121,file='rflex.txt',status='unknown')
        write(121,*)e_flex
        close(unit=121)
        endif


        return
        end
c***********************************************************************

