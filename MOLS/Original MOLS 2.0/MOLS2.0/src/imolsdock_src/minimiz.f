c-----------------------------------------------------------------------------
c     library name : minimiz.f   

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
        subroutine minimiz(kk,inp)
        include 'mols.par'

	common /comment/icomment
	common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /pepcrda/xpep(maxatm,8)
        common /pepcrdb/ypep(maxatm,8)
	common /patom/ipatom
	common /pout/popmi(maxstr,maxpar),pemi(maxstr)
	common /pepdbat/pepatom(maxatm),pepele(maxatm)
	common /par/natom,ntor,nhb,ns,lres
	common /gen/ibno,iat1(maxatm),iat2(maxatm),isno(maxatm)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
	common /pepcrdbf/ypepf(maxatm,8),atomf(maxatm),pepelef(maxatm)
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
        common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &         emo(maxstr),emi(maxstr)
        common /string1/str21(maxatm),str22(maxatm),str23,res2,
     &  atsym2(maxatm),adesc2(maxatm),attyp2(maxatm)
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        common /left/le(maxatm,maxatm),innd(maxatm),ielenum(maxatm),
     &  bonum(maxatm)
        common /string/atsym(maxatm),desc(maxatm),attyp(maxatm),
     &  str1(maxatm),str2(maxatm),str3,res,n1,n2,n3
        common /propdb/ihflex,patname(maxatm),px1(maxatm),px2(maxatm),
     &  px3(maxatm),presname(maxatm),pcid(maxatm),presnum(maxatm)
        common /plp/pat(mnatp),lat(maxatm),patnam(mnatp),
     &  presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),natp,
     &  tatyp(25,50),ntatp(25),px(mnatp,3)


        dimension phi(maxpar)

	character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
	character*80 str,prt,prt1,str11,str12,str13,res1,atsym1,adesc1,
     &  str21,str22,str23,res2,atsym2,adesc2,attyp1,attyp2,miniscr
	character*35 rmscr
	character*30 desc,str1,str2,str3,pepatom
	character*25 str4
	character*24 pepele
	character*17 str5
	character*13 str6
	character*12 prt2(mnatp)
	character*10 attyp,res
        character*6 atsym
	character*4 patname,qatname(mnatp),xx3,yy3
	character*3 presname,qresname(mnatp),xx1,yy1,matom(100)
	character*1 pcid,qcid(mnatp),xx2,yy2

        integer iatno1,ibno1,iatno2,ibno2,innd1(maxatm),isno1(maxatm),
     &  iat11(maxatm),iat12(maxatm),innd2(maxatm),iat21(maxatm),
     &  iat22(maxatm),isno2(maxatm),k11,nap,exresid,zero,fop,
     &  qresnum(mnatp),presnum,atno

        real molene,energy,rx1(100,3),rx2(100,3),fenergy,plpme,
     &  ex1,ex2,ex3,px1,px2,px3,qx1(mnatp),qx2(mnatp),qx3(mnatp),
     &  phi,rdihed,rvdw,reel,etot,mx1(100,3),yo(100,3),inter_ene,
     &  intra_ene,plpe,eflex

11      format(a25)
12      format(a17)
13      format(a13)
15      format(a30,3f8.3,a24)
21      format(a5,4x,i4,3x,f10.2)
22      format(a6)
23      format(a3)
63      format(a80)
67      format(1x,i5,1x,i5,1x,i5,a10)
68      format(3x,i4,1x,a6,2x,3f10.4,1x,a29)
130     format(a12,a4,1x,a3,1x,a1,i4,4x,3f8.3,24x)
170     format(a23,i3,4x,3f8.3,a24)
180     format(a6,i5,2x,a4,a3,5x,i1,4x,3f8.3)
641     format(i4,i4,i2,i2,i2)

	zero = 0
        k11=inp
        npar=ntor+6+np
        do i=1,npar
        phi(i)=opmi(k11,i)
        enddo

c=======PEPTIDE-PROTEIN=========================================
        IF(kk.eq.1)THEN

        nap=ipatom

        open(unit=11,file=pf7,status='old')
        open(unit=12,file='rec.pdb',status='old')
        open(unit=14,file=of5,access='append')

        if(ifopt.eq.1)then

         write(14,21)'MODEL',inp,emi(inp)
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

        write(14,21)'MODEL',inp,emi(inp)

        call flexall(100,100,phi,e_flex,1)

        do i=1,nap

        read(12,130)prt2(i),qatname(i),qresname(i),qcid(i),qresnum(i),
     &  qx1(i),qx2(i),qx3(i)

        do j=1,ihflex
        xx1=presname(j)
        xx2=pcid(j)
        xx3=patname(j)
        yy1=qresname(i)
        yy2=qcid(i)
        yy3=qatname(i)
        if(yy1(1:3).eq.xx1(1:3).and.yy2(1:1).eq.xx2(1:1).and.
     &  yy3(1:4).eq.xx3(1:4).and.presnum(j).eq.qresnum(i))then
        qx1(i)=px1(j)
        qx2(i)=px2(j)
        qx3(i)=px3(j)
        endif
        enddo
        enddo

        do i=1,nap
        write(14,130)prt2(i),qatname(i),qresname(i),qcid(i),qresnum(i),
     &  qx1(i),qx2(i),qx3(i)
        enddo

        write(14,23)'TER'
        do i=1,natom
        write(14,15)pepatom(i),ypep(i,1),ypep(i,2),
     &  ypep(i,3),pepele(i)
        enddo
        write(14,23)'TER'
        write(14,22)'ENDMDL'

        endif

        close(14)
        close(12)
        close(11)
c------------------------------------------------------------------------
	ELSE IF(kk.eq.2)THEN

        DO i = k11,k11

        do j = 1,npar-np-6
        phi(j) = opmi(i,j)
        enddo

        do j = npar-np-5,npar-np-3
        tt1 = opmi(i,j)
        if(tt1.eq.0.0) then
        tt2 = 0.0
        else
        tt2 = (tt1/10.0)*0.13888889
        endif
        tt2 = tt2-2.5
        phi(j) = tt2
        enddo

        do j = npar-np-2,npar
        phi(j) = opmi(i,j)
        enddo

        call smolgen(3,3,phi,2,ifopt,eflex)

        if(iff.eq.1)then
        intra_ene = rmmffene(i,j)
        else if(iff.eq.2)then
        intra_ene=rgaffene(i,j)
        endif
        call eplp(plpe,hb,steric,rep,2)
        inter_ene = plpe
        etot = intra_ene + inter_ene + eflex

        if(icomment.eq.1)print *,intra_ene,inter_ene,eflex,etot

        open(unit=1,file='mols.mol2',status='old')
        open(unit=2,file='min.mol2',status='unknown')
        open(unit=50,file=pf3,status='unknown',access='append')

        write(50,21)'MODEL',k11,etot
        read(1,63)str

        do j=1,2
        read (1,12)str5
        write(2,12)str5
        write(50,12)str5
        enddo

        read(1,641) inatom,ibno,n1,n2,n3
        write(2,641) inatom,ibno,n1,n2,n3
        write(50,641) inatom,ibno,n1,n2,n3

        do j=1,4
        read (1,13)str6
        write(2,13)str6
        write(50,13)str6
        enddo

        do j=1,inatom
        read(1,68)innd(j),atsym(j),yo(j,1),yo(j,2),yo(j,3),desc(j)
        write(2,68)innd(j),atsym(j),ylig(j,1),ylig(j,2),
     &  ylig(j,3),desc(j)
        write(50,68)innd(j),atsym(j),ylig(j,1),ylig(j,2),
     &  ylig(j,3),desc(j)
        enddo

        read (1,11)str4
        write(2,11)str4
        write(50,11)str4

        do k=1,ibno
        read(1,67)isno(k),iat1(k),iat2(k),attyp(k)
        write(2,67)isno(k),iat1(k),iat2(k),attyp(k)
        write(50,67)isno(k),iat1(k),iat2(k),attyp(k)
        enddo

        write(2,23)'TER'
        write(50,23)'TER'
        write(2,22)'ENDMDL'
        write(50,22)'ENDMDL'
        close(2)
        close(1)

c-------write ligand-protein complex--------------------------
        nap=ipatom

        if(ifopt.eq.1)then

        open(unit=1,file='rec.pdb',status='old')
        open(unit=2,file=of5,status='unknown',access='append')
        write(2,21)'MODEL',k11,etot
        do ij=1,nap
        read(1,170)prt,exresid,ex1,ex2,ex3,prt1
        write(2,170)prt,exresid,ex1,ex2,ex3,prt1
        enddo
        write(2,23)'TER'
        do ij=1,natom
        write(2,180)'HETATM',nap+ij+1,atsym(ij),'LIG',zero,ylig(ij,1),
     &  ylig(ij,2),ylig(ij,3)
        enddo
        write(2,23)'TER'
        write(2,22)'ENDMDL'
        close(2)
        close(1)

        else if(ifopt.eq.2)then

        open(unit=1,file='rec.pdb',status='old')
        open(unit=2,file=of5,status='unknown',access='Append')
        write(2,21)'MODEL',k11,etot

        do ij=1,nap
        read(1,130)prt2(ij),qatname(ij),qresname(ij),qcid(ij),
     &  qresnum(ij),qx1(ij),qx2(ij),qx3(ij)

        do j=1,ihflex
        xx1=presname(j)
        xx2=pcid(j)
        xx3=patname(j)
        yy1=qresname(ij)
        yy2=qcid(ij)
        yy3=qatname(ij)

        if(yy1(1:3).eq.xx1(1:3).and.yy2(1:1).eq.xx2(1:1)
     &  .and.yy3(1:4).eq.xx3(1:4)
     &  .and.presnum(j).eq.qresnum(ij))then
        qx1(ij)=px1(j)
        qx2(ij)=px2(j)
        qx3(ij)=px3(j)
        endif
        enddo
        enddo

        do ij=1,nap
         write(2,130)prt2(ij),qatname(ij),qresname(ij),qcid(ij),
     &  qresnum(ij),qx1(ij),qx2(ij),qx3(ij)
        enddo

        write(2,23)'TER'
        do ij=1,natom
        write(2,180)'HETATM',nap+ij+1,atsym(ij),'LIG',zero,ylig(ij,1),
     &  ylig(ij,2),ylig(ij,3)
        enddo
        write(2,23)'TER'
        write(2,22)'ENDMDL'
        close(2)
        close(1)

	endif

	ENDDO


        ENDIF
c=====================================================================
        return
        end
c*********************************************************************
        subroutine smolgen(lx,ly,phi,tst,ifopt,eflex)

        include 'mols.par'

	common /comment/icomment
	common /vectors/iv(maxpar,4)
	common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
	common /order/nn,mm
	common /cen/bx,by,bz
	common /par/natom,ntor,nhb,ns,lres
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &  iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),ind(maxatm,maxatm),
     &  big
        common /left/le(maxatm,maxatm),innd(maxatm),ielenum(maxatm),
     &  bonum(maxatm)


        dimension x_one(maxatm,3),x_two(maxatm,3),phi(maxpar)
	real rotang,theta,psi
        real cx,cy,cz,bx,by,bz
        integer tst,nn,ci

        if(icomment.eq.1.and.lx.eq.3.and.ly.eq.3)then
        do ii=1,ntor+6+np
        print *,'phi-final',ii,phi(ii)
        enddo
        endif


        nn=ntor
        do k=1,natom
           do ki=1,3
             x_one(k,ki)=rx(k,ki)
           enddo
        enddo
c------------------------------------------------------------------
        if(ntor.eq.0)then
        do k=1,natom
           do ki=1,3
             ylig(k,ki)=rx(k,ki)
           enddo
        enddo
        goto 334  
        endif    

        do if=1,nn
C###################### PHI ALL ####################################

        call elemen(x_one(irotb(if,1),1),x_one(irotb(if,1),2),
     &              x_one(irotb(if,1),3),
     &              x_one(irotb(if,2),1),x_one(irotb(if,2),2),
     &              x_one(irotb(if,2),3),
     &              el,em,en)


50      format(i4,1x,3f8.3)
        do j=1,ielenum(if) !ielenum = total number of ligand atoms
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

334     call rotate(bx,by,bz,phi(nn+4),phi(nn+5),phi(nn+6),natom)
        call translate(bx,by,bz,phi(nn+1),phi(nn+2),phi(nn+3),natom)
        if(ifopt.eq.2)then
        call flexall(lx,ly,phi,eflex,1)
        else
        eflex=0.0
        endif

        if(icomment.eq.1)then
        if(lx.eq.3.and.ly.eq.3)then
        do k=1,natom
        print *,'ylig',ylig(k,1),ylig(k,2),ylig(k,3)
        enddo
        endif
        endif

c*******************************************************************************
        return
        end

