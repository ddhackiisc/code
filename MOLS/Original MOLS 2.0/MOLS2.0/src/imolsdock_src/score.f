c-----------------------------------------------------------------------------
c     library name : score.f   

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

c-------Protein-Peptide Docking ---Peptide Energy--AMBER94---------------
        function pampene(ie,je)

        include 'mols.par'

	common /comment/icomment
	common /hbpar/mnhb
	common /calls/ncalls
        common /pepcrdb/ypep(maxatm,8)
	common /energy/e_final(maxord,maxord)
	common /tweight/xw1,xw2,xw3
	common /par/natom,ntor,nhb,ns,lres
	common /ambecep/ees,enb,ees1_4,ehb,enb1_4,etor
	common /hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
        common /ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,12)

        real y(maxatm,8),xw1,xw2,xw3

        ees=0.0
        enb=0.0
        enb1_4=0.0
        ees1_4=0.0
        ehb=0.0
        etot=0.0

        do i=1,natom
         do j=1,8
          y(i,j)=ypep(i,j)
         enddo
        enddo


        do i=1,natom
          n_range=ifix(y(i,7))
         do ij=1,n_range
          k1=jstart(i,ij)
         k2=jend(i,ij)
          do j=k1,k2
          dis = dist(y(i,1),y(i,2),y(i,3),y(j,1),y(j,2),y(j,3))

            if(dis.lt.0.01) then
               pampene = 1.0e+25
              return
            endif
            call pforce(i,j,dis,enbt,eest)
            ees=ees+eest
            enb=enb+enbt
          enddo
         enddo

         n1_4=ifix(y(i,8))
         do ij=1,n1_4
            j=j1_4(i,ij)
          dis = dist(y(i,1),y(i,2),y(i,3),y(j,1),y(j,2),y(j,3))
            if(dis.lt.0.5) then
              STOP 'Input coordinates are wrong !!'
            endif
            call pforce(i,j,dis,enbt,eest)
            ees1_4=ees1_4+0.5*eest
            enb1_4=enb1_4+0.5*enbt
         enddo
       enddo

        if(icomment.eq.1)print *,'mnhb',mnhb

        do i=1,mnhb
         dis = dist(y(ihb1(i),1),y(ihb1(i),2),y(ihb1(i),3),
     &                  y(ihb2(i),1),y(ihb2(i),2),y(ihb2(i),3))
         ess=dis*dis
         artwo=1.0/ess
         arten=artwo**5
         artwelve=arten*artwo
         ehb=ehb+c(i)*artwelve-d(i)*arten
         arstar=y(ihb1(i),5)+y(ihb2(i),5)
         epsilon=sqrt(y(ihb1(i),6)*y(ihb2(i),6))
         aaa=epsilon*(arstar**12)
         ccc=2.0*epsilon*(arstar**6)
         arsix=artwo*artwo*artwo
         ef1=aaa*artwelve
         ef2=ccc*arsix
         enbt=ef1-ef2
         enb=enb-enbt
        enddo
        
	pampene=(xw1*(enb+ees+ehb+enb1_4+ees1_4))

        if(icomment.eq.1)then
        print *,'PEPTIDE-AMBER'
        print *,'xw1',xw1
        print *,'enb',enb
        print *,'ees',ees
        print *,'ehb',ehb
        print *,'enb1_4',enb1_4
        print *,'ees1_4',ees1_4
        print *,'total',pampene
        endif
        return
        end
c-------------------------------------------------------------------
        subroutine pforce(if,jf,diss,enbt,eest)

        include 'mols.par'
	
	common /hbpar/mnhb
        common /pepcrdb/ypep(maxatm,8)
        common /par/natom,ntor,nhb,ns,lres

        real y(maxatm,8)

        do i=1,natom
         do j=1,8
          y(i,j)=ypep(i,j)
         enddo
        enddo

        arstar=y(if,5)+y(jf,5)
        epsilon=sqrt(y(if,6)*y(jf,6))
        aaa=epsilon*(arstar**12)
        ccc=2.0*epsilon*(arstar**6)
        ess=diss*diss
        artwo=1.0/ess
        arsix=artwo*artwo*artwo
        artwelve=arsix*arsix
        eest=332.0*y(if,4)*y(jf,4)*artwo/4.0
        ef1=aaa*artwelve
        ef2=ccc*arsix
        enbt=ef1-ef2
        return
        end
c-------Protein-Peptide Docking-----Peptide Energy--ECEPP/3----------
        function pecpene(ie,je)

        include 'mols.par'
        parameter(mxtyat = 18)

	common /comment/icomment
	common /calls/ncalls
        common /pepcrdb/ypep(maxatm,8)
	common /energy/e_final(maxord,maxord)
	common /tweight/xw1,xw2,xw3
	common /par/natom,ntor,nhb,ns,lres
	common /ambecep/ees,enb,ees1_4,ehb,enb1_4,etor
	common /hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
	common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common /ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,12)
        common /ecpp/aij(mxtyat,mxtyat),cij(mxtyat,mxtyat),
     &  a14(mxtyat,mxtyat),ihb(mxtyat,mxtyat),ahb(mxtyat,mxtyat),
     &  chb(mxtyat,mxtyat)

        real y(maxatm,8),xw1,xw2,xw3
        integer ie,je
        
        cdc=(22.0/(7.0*180))
        ees=0.0
        enb=0.0
        ehb=0.0
        etor=0.0
        enb14=0.0
        ees14=0.0
        ehb14=0.0
        etot=0.0

        do i=1,natom
         do j=1,8
          y(i,j)=ypep(i,j)
         enddo
        enddo

        print *,'ntor',ntor
        do in = 1,ntor
        etor = etor + 
     &         (u0(in)*(1.0+sn(in)*cos(cdc*((phi(in)-180.0))*tn(in))))
        enddo

        do i=1,natom
          n_range=ifix(y(i,7))

         do ij=1,n_range
          k1=jstart(i,ij)
         k2=jend(i,ij)
          do j=k1,k2
          dis = dist(y(i,1),y(i,2),y(i,3),y(j,1),y(j,2),y(j,3))
            if(dis.lt.0.01) then
               pecpene = 1.0e+25
              return
            endif
	ity = int(y(i,5))
        jty = int(y(j,5))
        ees = ees + (332.0*y(i,4)*y(j,4))/(dis*2.0)
        if(ihb(ity,jty).eq.0) then
        enb = enb + (aij(ity,jty)/(dis**12.0))-(cij(ity,jty)/(dis**6.0))
        else
         ehb=ehb+(ahb(ity,jty)/(dis**12.0))-(chb(ity,jty)/(dis**10.0))
        endif
          enddo
         enddo

         n1_4=ifix(y(i,8))
         do ij=1,n1_4
            j=j1_4(i,ij)
          dis = dist(y(i,1),y(i,2),y(i,3),y(j,1),y(j,2),y(j,3))
            if(dis.lt.0.5) then
              STOP 'Input coordinates are wrong !!'
            endif

        ity = int(y(i,5))
        jty = int(y(j,5))
        ees14 = ees14 + (332.0*y(i,4)*y(j,4))/(dis*2.0)
        if(ihb(ity,jty).eq.0) then
        enb14=enb14+(a14(ity,jty)/(dis**12.0))-(cij(ity,jty)/(dis**6.0))
        else
       ehb14=ehb14+(ahb(ity,jty)/(dis**12.0))-(chb(ity,jty)/(dis**10.0))
        endif
         enddo
       enddo

        ees = ees + ees14
        enb = enb + enb14
        ehb = ehb + ehb14

        pecpene=(xw1*(enb+ees+ehb+etor))

        if(icomment.eq.1)then
        print *,'ecp-elene1_4 =',ees14
        print *,'ecp-vwene1_4 =',enb14
        print *,'ecp-hbene1_4 =',ehb14
        print *,'ecp-elene =',ees
        print *,'ecp-vwene =',enb
        print *,'ecp-hbene =',ehb
        print *,ees,enb,ehb,etor
        print *,'total',pecpene
        endif
        return
        end
c-------Protein-Ligand Docking---Ligand Energy - MMFF94--------------
        function  rmmffene(ie,je) ! energy calculation using MMFF94

        include 'mols.par'


	common /comment/icomment
	common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /vectors/iv(maxpar,4)
	common /energy/e_final(maxord,maxord)
	common /native/nat,pepfile
	common /pdbat/atom(maxatm),ele(maxatm)
	common /tweight/xw1,xw2,xw3
	common /par/natom,ntor,nhb,ns,lres
	common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /gen/ibno,iat1(maxatm),iat2(maxatm),isno(maxatm)
        common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &  emo(maxstr),emi(maxstr)
        common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &  iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),ind(maxatm,maxatm),
     &  big
        common /left/le(maxatm,maxatm),innd(maxatm),ielenum(maxatm),
     &  bonum(maxatm)
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        common /string/atsym(maxatm),desc(maxatm),attyp(maxatm),
     &  str1(maxatm),str2(maxatm),str3,res,n1,n2,n3

        character desc*30,attyp*10,atsym*6,str1*30,str2*30,
     &  str3*30,res*10,rmscr*35,str*80,t*20
        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

        dimension pp(maxpar)

        real etot,rdihed,rvdw,reel,xw1,xw2,xw3

11      format(a25)
21      format(i4,i4,i2,i2,i2)
22      format(a6)
23      format(a3)
30      format(3x,i4,1x,a6,2x,3f10.4,1x,a29)
60      format(1x,i5,1x,i5,1x,i5,a10)

        open(unit=50,file='cal_energy.mol2',status='unknown')

        do i=1,2
        write(50,11)str1(i)
        enddo
        write(50,21)natom,ibno,n1,n2,n3

        if(nat.eq.0)then
        do i=1,4
        write(50,11)str2(i)
        enddo
        else if(nat.eq.1)then
        do i=1,5
        write(50,11)str2(i)
        enddo
        endif

        do i=1,natom
        write(50,30),innd(i),atsym(i),ylig(i,1),ylig(i,2),ylig(i,3),
     &  desc(i)
        enddo
        write(50,11)str3
        do j=1,ibno
        write(50,60)isno(j),iat1(j),iat2(j),attyp(j)
        enddo
        write(50,23)'TER'
        write(50,22)'ENDMDL'
        close(unit=50)

        str= 'sh ene_script.sh'
        call system(str)
        
        open(unit=113,file="output1",status= "old")

        read(113,*)rdihed
        read(113,*)rvdw
        read(113,*)reel

        etot= rdihed+rvdw+reel
        close(unit=113)

        rmmffene = etot
        rmscr='rm output1'
        call system(rmscr)
        if(icomment.eq.1)print *,'rmmffene-xw1',xw1
        rmmffene = (xw1*(etot))
        return
        end
c-------------------------------------------------------------------
        function rgaffene(ie,je) ! energy calculation using GAFF

        include 'mols.par'

	common /comment/icomment
	common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /vectors/iv(maxpar,4)
	common /energy/e_final(maxord,maxord)
	common /native/nat,pepfile
	common /pdbat/atom(maxatm),ele(maxatm)
	common /tweight/xw1,xw2,xw3
	common /par/natom,ntor,nhb,ns,lres
	common /gen/ibno,iat1(maxatm),iat2(maxatm),isno(maxatm)
	common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &  emo(maxstr),emi(maxstr)
        common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &  iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),
     &  ind(maxatm,maxatm),big
        common /left/le(maxatm,maxatm),innd(maxatm),ielenum(maxatm),
     &  bonum(maxatm)
        common /string/atsym(maxatm),desc(maxatm),attyp(maxatm),
     &  str1(maxatm),str2(maxatm),str3,res,n1,n2,n3


        character desc*30,attyp*10,atsym*6,str1*30,str2*30,
     &  str3*30,res*10,rmscr*35,str*80,t*20
        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

        dimension pp(maxpar)
      
        real etot,rdihed,rpdihed,rvdw,reel,xw1,xw2,xw3

11      format(a25)
21      format(i4,i4,i2,i2,i2)
22      format(a6)
23      format(a3)
30      format(3x,i4,1x,a6,2x,3f10.4,1x,a29)
60      format(1x,i5,1x,i5,1x,i5,a10)

        open(unit=50,file='gcal_energy.mol2',status='unknown')
        do i=1,2
        write(50,11)str1(i)
        enddo
        write(50,21)natom,ibno,n1,n2,n3

        if(nat.eq.0)then
        do i=1,4
        write(50,11)str2(i)
        enddo
        else if(nat.eq.1)then
        do i=1,5
        write(50,11)str2(i)
        enddo
        endif

        do i=1,natom
        write(50,30),innd(i),atsym(i),ylig(i,1),ylig(i,2),ylig(i,3),
     &  desc(i)
        enddo

        write(50,11)str3
        do j=1,ibno
        write(50,60)isno(j),iat1(j),iat2(j),attyp(j)
        enddo
        write(50,23)'TER'
        write(50,22)'ENDMDL'
        close(unit=50)

        str= 'sh gene_script.sh'
        call system(str)

        open(unit=113,file="output1",status= "old")

        read(113,*)rdihed
        read(113,*)rpdihed
        read(113,*)rvdw
        read(113,*)reel

        etot= rdihed+rpdihed+rvdw+reel
        close(unit=113)

        rgaffene = etot
        rmscr='rm output1'
        call system(rmscr)
        if(icomment.eq.1)print *,'rgaffene-xw1',xw1
        rgaffene = (xw1*(etot))
        if(icomment.eq.1)print *,'rgaffene',rgaffene
        return
        end
c--------------------------------------------------------------------
c-------Peptide-Protein PLP Scoring----------------------------------
c-------Precalculation of Interaction Pair---------------------------
        subroutine pprecal
        include 'mols.par'

        common /comment/icomment
        common /pplp/prresid(mnatp),matp
        common /par/natom,ntor,nhb,ns,lres
        common /preset/hbp(100000),hbl(100000),sp(100000),sl(100000),
     &          iscnt,ihcnt
        common /plp/pat(mnatp),lat(maxatm),patnam(mnatp),
     &  presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),
     &  natp,tatyp(25,50),ntatp(25),px(mnatp,3)

        character*4 xatnam,patnam,tatnam,ltyp,ptyp
        character*3 xresnam,presnam,tpres
        character*2 lat,pat,tatyp
        character*1 pchid

        ik = 1
        natp=matp
       if(icomment.eq.1)print *,'natp-precal',natp,'natom',natom
        do i = 1,natp
        ptyp = pat(i)
        do j = 1,natom
        ltyp = lat(j)
      if((ptyp.eq.'D'.and.ltyp.eq.'A').or.
     & (ptyp.eq.'D'.and.ltyp.eq.'DA').or.
     & (ptyp.eq.'A'.and.ltyp.eq.'D').or.
     & (ptyp.eq.'A'.and.ltyp.eq.'DA').or.
     & (ptyp.eq.'DA'.and.ltyp.eq.'D').or.
     & (ptyp.eq.'DA'.and.ltyp.eq.'A').or.
     & (ptyp.eq.'DA'.and.ltyp.eq.'DA')) then
        hbp(ik) = i
        hbl(ik) = j
        ik = ik + 1
        endif
        enddo
        enddo
        ihcnt = ik-1
        if(icomment.eq.1)print *,'ihcnt-precal',ihcnt
        ik = 1
        do i = 1,natp
        ptyp = pat(i)
        do j = 1,natom
        ltyp = lat(j)
      if((ptyp.eq.'D'.and.ltyp.eq.'D').or.
     & (ptyp.eq.'D'.and.ltyp.eq.'NP').or.
     & (ptyp.eq.'A'.and.ltyp.eq.'A').or.
     & (ptyp.eq.'A'.and.ltyp.eq.'NP').or.
     & (ptyp.eq.'DA'.and.ltyp.eq.'NP').or.
     & (ptyp.eq.'NP'.and.ltyp.eq.'A').or.
     & (ptyp.eq.'NP'.and.ltyp.eq.'D').or.
     & (ptyp.eq.'NP'.and.ltyp.eq.'DA').or.
     & (ptyp.eq.'NP'.and.ltyp.eq.'NP')) then
        sp(ik) = i
        sl(ik) = j
        ik = ik + 1
        endif
        enddo
        enddo
        iscnt = ik-1
        if(icomment.eq.1)print *,'iscnt-precal',iscnt

        return
        end
c-------End of Precalculation----------------------------------------
c**************PLP Energy calculation********************************
        subroutine peplp(plpe,hb,steric,rep,tst)
        include 'mols.par'

        common /comment/icomment
        common /pplp/prresid(mnatp)
        common /pepcrdb/ypep(maxatm,8)
        common /tweight/xw1,xw2,xw3
        common /par/natom,ntor,nhb,ns,lres
        common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /preset/hbp(100000),hbl(100000),sp(100000),sl(100000),
     &          iscnt,ihcnt
        common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        common /propdb/ihflex,patname(mnatp),px1(mnatp),px2(mnatp),
     &  px3(mnatp),presname(mnatp),pcid(mnatp),presnum(mnatp)
        common /plp/pat(mnatp),lat(maxatm),patnam(mnatp),
     &  presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),
     &  natp,tatyp(25,50),ntatp(25),px(mnatp,3)


        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6,pepfile
        character*50 strz(200)
        character*30 stry(200)
        character*4 xatnam,patnam,tatnam,ltyp,ptyp,patname,xx3,yy3,
     &  atname(200),rename(200)
        character*3 xresnam,presnam,tpres,presname,xx1,yy1
        character*2 lat,pat,tatyp
        character*1 pcid,pchid,xx2,yy2

        real d,a,b,c,y1,ene,plpe,hb,steric,xw1,xw2,xw3,
     &  rep,dr
        integer tst,presnum,ihflex,prresid,iresno(200),iatmno(200)

10      format(8x,i3,1x,a4,1x,a4,3x,i3)
12      format(12x,a4,1x,a3,1x,a1,i4,4x,3f8.3,24x)

        d = 0.0
        a = 0.0
        b = 0.0
        c = 0.0
        ene = 0.0
        hb = 0.0
        steric = 0.0
	rep = 0.0
	dr = 0.0

c*****************************************************************
        open(unit=1,file=pf1,status='old')
        do i=1,natom
            read(1,10) iatmno(i),atname(i),rename(i),iresno(i)
        enddo
        close(unit=1)

c*****************************************************************

c       -------for flexible receptor docking---------------------
        IF(ifopt.eq.2)then

        open(unit=9,file='hflexres_new.pdb',status='old')
        do i=1,ihflex
        read(9,12)patname(i),presname(i),pcid(i),presnum(i),
     &  px1(i),px2(i),px3(i)
        enddo
        close(9)

        if(icomment.eq.1)print *,'peplp-ihflex',ihflex,natp
        do i=1,ihflex
        do j=1,natp
        xx1=presname(i)
        xx2=pcid(i)
        xx3=patname(i)
        yy1=presnam(j)
        yy2=pchid(j)
        yy3=patnam(j)
        if(xx1(1:3).eq.yy1(1:3).and.xx2(1:1).eq.yy2(1:1)
     &  .and.xx3(1:4).eq.yy3(1:4)
     &  .and.presnum(i).eq.prresid(j))then
        px(j,1)=px1(i)
        px(j,2)=px2(i)
        px(j,3)=px3(i)
        endif
        enddo
        enddo    

        ENDIF
c-------PLP-HB---------------------------------------------------
        do i = 1,ihcnt
        ip = hbp(i)
        il = hbl(i)

        d = dist(px(ip,1),px(ip,2),px(ip,3),ypep(il,1),ypep(il,2),
     & ypep(il,3))

        if(d.le.3.4) then
         if(d.le.2.3.and.d.ge.0) then
         a = 20.0
         b = 2.3
         c = 46.0
         else
          if(d.le.2.6.and.d.ge.2.3) then
          a = 2.0
          b = 0.3
          c = 4.6
          else
           if(d.le.3.1.and.d.ge.2.6) then
           a = 0.0
           b = 1.0
           c = -2.0
         else
            if(d.le.3.4.and.d.ge.3.1) then
            a = 2.0
            b = -0.3
            c = 6.8
            endif
           endif
          endif
         endif
       	y1 = (c-(a*d))/b
        hb = hb + y1
        endif
        enddo
c-------PLP-STERIC---------------------------------------------
        do j = 1,iscnt
        ip = sp(j)
        il = sl(j)
        d = dist(px(ip,1),px(ip,2),px(ip,3),ypep(il,1),ypep(il,2),
     &  ypep(il,3))
        if(d.le.5.5) then
         if(d.le.3.4.and.d.ge.0) then
         a = 20.0
         b = 3.4
         c = 68.0
         else
          if(d.le.3.6.and.d.ge.3.4) then
          a = 0.4
          b = 0.2
          c = 1.36
          else
           if(d.le.4.5.and.d.ge.3.6) then
           a = 0.0
           b = 1.0
           c = -0.4
         else
            if(d.le.5.5.and.d.ge.4.5) then
            a = 0.4
            b = -1.0
            c = 2.2
            endif
           endif
          endif
         endif
        y1 = (c-(a*d))/b
        steric = steric + y1
        endif
        enddo
c*******PLP-REPULSIVE******************************************
        do j = 1,iscnt
        ip = sp(j)
        il = sl(j)
        dr = dist(px(ip,1),px(ip,2),px(ip,3),ypep(il,1),ypep(il,2),
     &   ypep(il,3))

        a = 3.2
        b = 5.0
        c = 0.1
        d = 20.0

        if(dr.lt.a)then
           y2 = ((dr*((c-d)/a))+d)
        else
         if(dr.le.b.and.dr.ge.a)then
          y2 = ((-c*((dr-a)/(b-a)))+c)
         else
          if(dr.gt.b)then
           y2 = 0.0
          endif
         endif
        endif
        rep = rep + y2
        enddo
c**************************************************************
        plpe = (xw2*(hb+steric+rep))
        if(icomment.eq.1)print *,'test-plp-xw2-plpe',xw2,plpe
        if(icomment.eq.1.and.tst.eq.2)write (*,*),'plp ene 
     &  from minimiz: ',plpe,hb,steric,rep
        return
        end
c**************End of Peptide-Protein PLP Score**********************

c*******Protein-Ligand Docking***************************************
c-------Precalculate the interaction pair----------------------------
        subroutine precal

        include 'mols.par'

	common /comment/icomment
	common /pplp/prresid(mnatp),matp
	common /par/natom,ntor,nhb,ns,lres
        common /preset/hbp(100000),hbl(100000),sp(100000),sl(100000),
     &          iscnt,ihcnt
        common /plp/pat(mnatp),lat(maxatm),patnam(mnatp),
     &          presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),
     &          natp,tatyp(25,50),ntatp(25),px(mnatp,3),hb,steric

	character*4 xatnam,patnam,tatnam,ltyp,ptyp
	character*3 xresnam,presnam,tpres
        character*2 lat,pat,tatyp

        if(icomment.eq.1)print *,'precal-matp',matp

        natp=matp
        ik = 1
        do i = 1,natp
        ptyp = pat(i)
        do j = 1,natom
        ltyp = lat(j)
      if((ptyp.eq.'D'.and.ltyp.eq.'A').or.
     & (ptyp.eq.'D'.and.ltyp.eq.'DA').or.
     & (ptyp.eq.'A'.and.ltyp.eq.'D').or.
     & (ptyp.eq.'A'.and.ltyp.eq.'DA').or.
     & (ptyp.eq.'DA'.and.ltyp.eq.'D').or.
     & (ptyp.eq.'DA'.and.ltyp.eq.'A').or.
     & (ptyp.eq.'DA'.and.ltyp.eq.'DA')) then
        hbp(ik) = i
        hbl(ik) = j

        ik = ik + 1
        endif
        enddo
        enddo

        ihcnt = ik-1

        ik = 1

        do i = 1,natp
        ptyp = pat(i)
        do j = 1,natom
        ltyp = lat(j)
        if((ptyp.eq.'D'.and.ltyp.eq.'D').or.
     & (ptyp.eq.'D'.and.ltyp.eq.'NP').or.
     & (ptyp.eq.'A'.and.ltyp.eq.'A').or.
     & (ptyp.eq.'A'.and.ltyp.eq.'NP').or.
     & (ptyp.eq.'DA'.and.ltyp.eq.'NP').or.
     & (ptyp.eq.'NP'.and.ltyp.eq.'A').or.
     & (ptyp.eq.'NP'.and.ltyp.eq.'D').or.
     & (ptyp.eq.'NP'.and.ltyp.eq.'DA').or.
     & (ptyp.eq.'NP'.and.ltyp.eq.'NP')) then
        sp(ik) = i
        sl(ik) = j
        ik = ik + 1
        endif
        enddo
        enddo

        iscnt = ik-1

        return
        end
c*******End of Precalculation*****************************************

c**************PLP Energy calculation*********************************
        subroutine eplp(plpe,hb,steric,rep,tst)

        include 'mols.par'

	common /comment/icomment
        common /ligcrdb/ylig(maxatm,8)
	common /pplp/prresid(mnatp),matp
	common /tweight/xw1,xw2,xw3
	common /par/natom,ntor,nhb,ns,lres
        common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /preset/ hbp(100000),hbl(100000),sp(100000),sl(100000),
     &          iscnt,ihcnt
        common /propdb/ihflex,patname(mnatp),px1(mnatp),px2(mnatp),
     &  px3(mnatp),presname(mnatp),pcid(mnatp),presnum(mnatp)
        common /plp/ pat(mnatp),lat(maxatm),patnam(mnatp),
     &  presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),natp,
     &  tatyp(25,50),ntatp(25),px(mnatp,3)

	character*4 xatnam,patnam,tatnam,ltyp,ptyp,patname,xx3,yy3
	character*3 xresnam,presnam,tpres,presname,xx1,yy1
        character*2 lat,pat,tatyp
	character*1 pcid,pchid,xx2,yy2

        real d,a,b,c,y1,y2,ene,plpe,hb,steric,rep,dr,xw1,xw2,xw3
        integer tst,presnum,ihflex,prresid

130     format(12x,a4,1x,a3,1x,a1,i4,4x,3f8.3,24x)

        d = 0.0
        a = 0.0
        b = 0.0
        c = 0.0
        ene = 0.0
        hb = 0.0
        steric = 0.0
        rep = 0.0
        dr = 0.0

c****************************************************************
c       -------for flexible receptor docking---------------------
        IF(ifopt.eq.2)then

        open(unit=9,file='hflexres_new.pdb',status='old')
        do i=1,ihflex
        read(9,130)patname(i),presname(i),pcid(i),presnum(i),
     &  px1(i),px2(i),px3(i)
        enddo
        close(9)

        if(icomment.eq.1)print *,'ihflex',ihflex,natp

        do i=1,ihflex
        do j=1,natp
        xx1=presname(i)
        xx2=pcid(i)
        xx3=patname(i)
        yy1=presnam(j)
        yy2=pchid(j)
        yy3=patnam(j)

        if(xx1(1:3).eq.yy1(1:3).and.xx2(1:1).eq.yy2(1:1)
     &  .and.xx3(1:4).eq.yy3(1:4)
     &  .and.presnum(i).eq.prresid(j))then
        px(j,1)=px1(i)
        px(j,2)=px2(i)
        px(j,3)=px3(i)

        endif
        enddo
        enddo

        ENDIF
c       --------------------------------------------------

        if(icomment.eq.1)print *,'ihcnt,iscnt',ihcnt,iscnt

c-----------------------------------------------------------------

        do i = 1,ihcnt
        ip = hbp(i)
        il = hbl(i)
        d = dist(px(ip,1),px(ip,2),px(ip,3),ylig(il,1),ylig(il,2),
     &  ylig(il,3))

        if(d.le.3.4) then
         if(d.le.2.3.and.d.ge.0) then
         a = 20.0
         b = 2.3
         c = 46.0
         else
          if(d.le.2.6.and.d.ge.2.3) then
          a = 2.0
          b = 0.3
          c = 4.6
          else
           if(d.le.3.1.and.d.ge.2.6) then
           a = 0.0
           b = 1.0
           c = -2.0
         else
         if(d.le.3.4.and.d.ge.3.1) then
            a = 2.0
            b = -0.3
            c = 6.8
            endif
           endif
        endif
       endif
       y1 = (c-(a*d))/b
        hb = hb + y1
       endif
        enddo
        do j = 1,iscnt
        ip = sp(j)
        il = sl(j)
        d = dist(px(ip,1),px(ip,2),px(ip,3),ylig(il,1),ylig(il,2),
     &   ylig(il,3))
c-----------------------------------------------------------------
        if(d.le.5.5) then
         if(d.le.3.4.and.d.ge.0) then
         a = 20.0
         b = 3.4
         c = 68.0
         else
          if(d.le.3.6.and.d.ge.3.4) then
          a = 0.4
          b = 0.2
          c = 1.36
          else
           if(d.le.4.5.and.d.ge.3.6) then
           a = 0.0
           b = 1.0
           c = -0.4
         else
          if(d.le.5.5.and.d.ge.4.5) then
            a = 0.4
            b = -1.0
            c = 2.2
            endif
           endif
        endif
       endif
       y1 = (c-(a*d))/b
        steric = steric + y1
      endif
      enddo
c*******PLP-REPULSIVE**********************************************
        do j = 1,iscnt
        ip = sp(j)
        il = sl(j)
        dr = dist(px(ip,1),px(ip,2),px(ip,3),ylig(il,1),ylig(il,2),
     &   ylig(il,3))

        a = 3.2
        b = 5.0
        c = 0.1
        d = 20.0

        if(dr.lt.a)then
           y2 = ((dr*((c-d)/a))+d)
        else
          if(dr.le.b.and.dr.ge.a)then
           y2 = ((-c*((dr-a)/(b-a)))+c)
          else
            if(dr.gt.b)then
            y2 = 0.0
            endif
          endif
        endif
           rep = rep + y2
        enddo
c*****************************************************************
        if(icomment.eq.1)print *,'plpe-smmols-xw2',xw2
        plpe = (xw2*(hb+steric+rep))
        if(icomment.eq.1)print *,'hb',hb,'steric',steric,'rep',rep

        return
        end
c-----------------------------------------------------------------
c*******Flexible Protein Energy Calculation***********************
        function afampene(hx,hy)
        include 'mols.par'
        
        common /comment/icomment
        common /flex/iponly
        common /provectors/ivp(maxpar,4)
        common /part/ipart,iptor
        common /tweight/xw1,xw2,xw3
        common /pcrdb/nl,nptor,nphb,yy(maxatm,8)
        common /protene/penb,pees,penb1_4,pees1_4
        common /phb/iphb1(maxhb),iphb2(maxhb),ch(maxhb),dh(maxhb)
        common /pranges/jpstart(maxatm,10),jpend(maxatm,10),
     &  jp1_4(maxatm,25)
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        common /partners/qatname(maxatm),qx1(maxatm),qx2(maxatm),
     &  qx3(maxatm),qresname(maxatm),qcid(maxatm),qresnum(maxatm)
        common /propdb/ihflex,patname(maxatm),px1(maxatm),px2(maxatm),
     &  px3(maxatm),presname(maxatm),pcid(maxatm),presnum(maxatm)

        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        character*4 patname,qatname,xx3,yy3
        character*3 qresname,presname,xx1,yy1
        character*1 pcid,qcid,xx2,yy2

        integer qresnum,presnum
        integer resno,jj,ji,hx,hy,jp(5000),jr(5000),iponly
        real enb,ees,enb1_4,ees1_4,penb,pees,penb1_4,pees1_4
c--------------------------------------------------------------------
        do i=1,ihflex
        do j=1,ipart
        if(icomment.eq.2)icount=icount+1
        xx1=presname(i)
        xx2=pcid(i)
        xx3=patname(i)
        yy1=qresname(j)
        yy2=qcid(j)
        yy3=qatname(j)
        if(xx1(1:3).eq.yy1(1:3).and.xx2(1:1).eq.yy2(1:1)
     &  .and.xx3(1:4).eq.yy3(1:4)
     &  .and.presnum(i).eq.qresnum(j))then
        yy(j,1)=px1(i)
        yy(j,2)=px2(i)
        yy(j,3)=px3(i)
        endif
        enddo
        enddo
c-------------------------------------------------------------------
        ees=0.0
        enb=0.0
        enb1_4=0.0
        ees1_4=0.0
        ehb=0.0
        etot=0.0
c--------------------------------------------------------------------
        do i=1,nl
          n_range=ifix(yy(i,7))
         do ij=1,n_range
          k1=jpstart(i,ij)
         k2=jpend(i,ij)
          do j=k1,k2
          dis = dist(yy(i,1),yy(i,2),yy(i,3),
     &    yy(j,1),yy(j,2),yy(j,3))
            if(dis.lt.0.01) then
               afampene = 1.0e+25
            return
            endif
            call fforce(i,j,dis,enbt,eest)
            ees=ees+eest
            enb=enb+enbt
          enddo
         enddo

         n1_4=ifix(yy(i,8))
         do ij=1,n1_4
            j=jp1_4(i,ij)
          dis = dist(yy(i,1),yy(i,2),yy(i,3),
     &   yy(j,1),yy(j,2),yy(j,3))

            call fforce(i,j,dis,enbt,eest)

            ees1_4=ees1_4+0.5*eest
            enb1_4=enb1_4+0.5*enbt
         enddo
        enddo


         do i=1,nphb
         dis = dist(yy(iphb1(i),1),yy(iphb1(i),2),yy(iphb1(i),3),
     &              yy(iphb2(i),1),yy(iphb2(i),2),yy(iphb2(i),3))
         ess=dis*dis
         artwo=1.0/ess
         arten=artwo**5
         artwelve=arten*artwo
         ehb=ehb+ch(i)*artwelve-dh(i)*arten
         arstar=yy(iphb1(i),5)+yy(iphb2(i),5)
         epsilon=sqrt(yy(iphb1(i),6)*yy(iphb2(i),6))
         aaa=epsilon*(arstar**12)
         ccc=2.0*epsilon*(arstar**6)
         arsix=artwo*artwo*artwo
         ef1=aaa*artwelve
         ef2=ccc*arsix
         enbt=ef1-ef2
         enb=enb-enbt
        enddo
        afampene=ees+enb+ees1_4+enb1_4

c-------Protein Energy for promols/promini energy log file------
        penb=enb
        pees=ees
        penb1_4=enb1_4
        pees1_4=ees1_4

        if(icomment.eq.1)then
        print *,'PROTEIN-AMBER'
        write(*,*)'ees',pees
        write(*,*)'enb-enbt',penb
        write(*,*)'ees1_4',pees1_4
        write(*,*)'enb1_4',penb1_4
        write(*,*)'Total Energy',afampene
        endif

        return
        end
c**************************************************************
        subroutine fforce(if,jf,diss,enbt,eest)

        include 'mols.par'
        common /pcrdb/nl,nptor,nphb,yy(maxatm,8)

        arstar=yy(if,5)+yy(jf,5)
        epsilon=sqrt(yy(if,6)*yy(jf,6))
        aaa=epsilon*(arstar**12)
        ccc=2.0*epsilon*(arstar**6)
        ess=diss*diss
        artwo=1.0/ess
        arsix=artwo*artwo*artwo
        artwelve=arsix*arsix
        eest=332.0*yy(if,4)*yy(jf,4)*artwo/4.0
        ef1=aaa*artwelve
        ef2=ccc*arsix
        enbt=ef1-ef2
        return
        end
c**************************************************************
