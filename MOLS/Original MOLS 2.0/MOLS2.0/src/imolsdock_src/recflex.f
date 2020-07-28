c-----------------------------------------------------------------------------
c     library name : recflex.f   

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

c-------Prepare Protein for Receptor Flexibility------------------
        subroutine recflex()
        include 'mols.par'
        parameter(maxn=2000)

	common /comment/icomment
	common /part/ipart
	common /rctl/iscopt
	common /patom/ipatom
	common /ddihed/ddih(maxatm)
	common /rcrda/xx_one(maxatm,3)
	common /pplp/prresid(mnatp),matp
	common /recsc/iline,idihed,ivx(maxatm,4)
        common /recep/np,nres,crnam(maxatm),crcid(maxatm),
     &  irrid(maxatm)
        common /sscc/ssc(maxatm,maxpar,5),nsrot(maxpar),
     &  nschi(maxpar)
        common /partners/qatname(maxatm),qx1(maxatm),qx2(maxatm),
     &  qx3(maxatm),qresname(maxatm),qcid(maxatm),qresnum(maxatm)

	character*80 rmf1,rcount,lcount,addh,pcount,rmpartner,
     &  frm,str,partexe,castr,cbstr,sstr(maxatm)
        character*4 pratnam(mnatp),rratnam(maxatm),qatname,
     &  tcode,xatnam,at1,at2,at3,at4,at5,at6,at7,attnam(maxatm),
     &  dd1,dd2,dd3,dd4
        character*3 presnam(mnatp),rresnam(maxatm),qresname,
     &  rname,rcode,ccode,crnam,xresnam
        character*1 pchid(mnatp),rchid(maxatm),qcid,chid(10),
     &  cid,rcode1,crcid,fseq(100),ch,scode,xchid
	character pl1*22,pl2*54

        integer presid(mnatp),rresid(maxatm),rres,pres,
     &  kres,sres(maxatm),ssres(maxatm),temp,ireno(maxatm),
     &  s,ivv(maxatm,4),rt,nrt(maxatm),dhd(maxatm,maxatm,4),
     &  nd,nres,np,nr,rid,ml,nl,qresnum,ncid

	real x1,x2,x3,cutdist,px(mnatp,3),rx(maxatm,3),
     &  ax,ay,az,cx,cy,cz,cutd,ddih,xx_one,qx1,qx2,qx3

10      format(a3,1x,a1,1x,i1)
20      format(a4,a1,1x,i2,i2)
21      format(20x,5f7.1)
22      format(8x,4a4)
30      format(i4)
40      format(a12,a4,6x,i4,a52)
50      format(30x,f8.3,f8.3,f8.3)
80      format(i4,2x,a4)
100     format(i4,2x,i4,2x,i4,2x,i4)
130     format(5x,a42,2x,i4)
190	format(a22,i4,a54)
191	format(12x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
192     format(a80)
200     format(12x,a4,6x,i4,52x)
210     format(a4,a1,1x,i2)
220     format(8x,4a4,26x)
221     format(12x,a4,1x,a3,1x,a1,i4,4x,3f8.3,24x)
223     format(5x,a27,2x,i4)
224	format(5x,i2,a1,2x,a3,2x,a1,2x,i4)

        np=0
        rcount='cat residues.txt|wc -l>rlines.txt'
        call system(rcount)
        open(unit=1,file='rlines.txt',status='unknown')
        read(1,*)nres !nres=number of flexible residues
	write(31,'(A)')
        write(31,223) 'Protein Flexible Residues :',nres
	write(31,'(/5X,20A/)')'--------------------'
        close(unit=1)

!-------store residue name and id in rnam & rrid -------------------- 
        open(unit=101,file='residues.txt',status='old')
        do ii=1,nres !number of protein residues for SC flexibility
        read(101,*),crnam(ii),crcid(ii),irrid(ii)
        write(31,224) ii,')',crnam(ii),crcid(ii),irrid(ii)
        enddo
	write(31,'(/5X,20A/)')'--------------------'
        close(101)
!.......rotamer......................................................
        if(iscopt.eq.1)then
        print *,'iscopt',iscopt
        open(unit=33,file='SCROT.lib',status='old')
        do i = 1,nres
        rewind 33
        do i1 = 1,234
        read(33,20) tcode, scode, mrot, mchi
c       write(31,*)tcode, scode, mrot, mchi
        ccode=crnam(i)
        if(tcode(1:3).eq.ccode(1:3).and.mrot.ne.0) then
          do i2 = 1, mrot
          read(33,21) (ssc(i,i2,i3),i3=1,mchi)

          if(icomment.eq.1)print *,'rotamer',(ssc(i,i2,i3),i3=1,mchi)
          enddo

        nsrot(i) = mrot
        nschi(i) = mchi

        go to 212
        endif
        enddo
212     enddo
        close(unit=33)
        endif

!------total number of side chain torsion angles ---------------------
        do ii=1,nres
        rname=crnam(ii)
        
          open(unit=3,file='SCDIHEDS.lib',status='unknown')
           do ll=1,111
            read(3,10)rcode,rcode1,nr
             if(rcode.eq.rname)then
              np=np+nr
            endif
           enddo
          close(3)
        enddo

        write(31,130)'Dihedral angles in Protein Flex. Residues:',np

        if(icomment.eq.1)print *,'total no. of side chain torsion 
     & angles = ',np

        frm='rm flexres'
        call system(frm)

        open(unit=116,file='flexres',access='append')
        DO  ij=1,nres
c#########################################################################

        rid=irrid(ij)
        cid=crcid(ij)

        open(unit=115,file='rec.pdb',status='old')
        do mi=1,ipatom
        read(115,190)pl1,ml,pl2
        if((pl1(22:22).eq.cid).and.(ml.eq.rid))then
        write(116,190)pl1,ml,pl2
        endif
        enddo
        close(115)
        ENDDO
        close(116)

c########################################################################
        lcount='cat flexres|wc -l>lines.txt'
        call system(lcount)

        open(unit=1,file='lines.txt',status='unknown')
        read(1,*)nl !nl=number of lines in flexres
        close(unit=1)
c-----------------------------------------------------------------------
        frm='rm rbno.txt'
        call system(frm)

        open(unit=5,file='rbno.txt',access='append',
     &             status='unknown')
        DO ij=1,nres
        open(unit=3,file='SCDIHEDS.lib',status='old')

        do j=1,111 !first do
         read(3,10)rcode,rcode1,nr
            if(rcode.eq.crnam(ij))then !first if
            if(icomment.eq.1)write(*,*)nr
            if(nr.ne.0)then !second if
             do k=1,nr ! second do
              read(3,22)at1,at2,at3,at4
              open(unit=4,file='flexres',status='old')
              if(icomment.eq.1)write(*,*)'nl',nl
               do s=1,nl
               read(4,40)at5,at6,irno,at7
               if(at6.eq.at1.and.irno.eq.irrid(ij))then
               write(5,80),s,at6
               endif
               if(at6.eq.at2.and.irno.eq.irrid(ij))then
               write(5,80),s,at6
               endif
               if(at6.eq.at3.and.irno.eq.irrid(ij))then
               write(5,80),s,at6
               endif
               if(at6.eq.at4.and.irno.eq.irrid(ij))then
               write(5,80),s,at6
               endif
              enddo
             close(4)

            enddo  !second do
           endif  !second if
        endif  !first if      
        enddo  !first do

         close(3)
        ENDDO
        close(5)

        frm='rm dihed.txt'
        call system(frm)

        idihed=0
        open(unit=6,file='dihed.txt',access='append')
        open(unit=5,file='rbno.txt',status='old')
        do ij=1,nres
        open(unit=3,file='SCDIHEDS.lib',status='unknown')
        do j=1,111!first do
          read(3,10)rcode,rcode1,nr
          if(rcode.eq.crnam(ij))then !first if
            if(nr.ne.0)then  !second if     
             do i=1,nr!second do
              idihed=idihed+1
              do k=1,4 !third do
               read(5,30)ivv(i,k)
                ivx(idihed,k)=ivv(i,k)
              enddo    !third do
              write(6,100),(ivv(i,l),l=1,4)
c             write(*,100),(ivv(i,l),l=1,4)
             enddo !second do
            endif !second if
          endif !first if
        enddo !first do
        close(3)

        ENDDO
        close(5)
        close(6)

        if(icomment.eq.1)then
        do i=1,idihed
        print *,'ivx-recflex',ivx(i,1),ivx(i,2),ivx(i,3),ivx(i,4)
        enddo
        endif
c--------------------------------------------------------------
        open(unit=1,file='flexres',status='old')
        do iline=0,maxn
        read(1,192,end=299),str
        enddo
299     if(icomment.eq.1)write(*,*),'line number',iline
        close(1)
c--------------------------------------------------------------
        open(unit=6,file='flexres',status='old')
        do k=1,iline
        read(6,50),(xx_one(k,ki),ki=1,3)
        enddo
        close(6)
c-------PROTEIN SIDE-CHAIN ROTAMER FLEXIBILITY-----------------
        if(iscopt.eq.1)then
        j3=0
cs      do ii=1,nres
        open(unit=61,file='flexres',status='old')
        do k=1,iline
        read(61,200),attnam(k),ireno(k)
        enddo
        close(61)

        do ii=1,nres
        open(unit=62,file='SCVAR.lib',status='old')
        do j=1,57
        read(62,210)tcode,scode,natm
         if(tcode.eq.crnam(ii))then
          do j1=1,natm
           idd1= 0
           idd2= 0
           idd3= 0
           idd4= 0
           read(62,220),dd1,dd2,dd3,dd4
            do j2=1,iline
             if(dd1.eq.attnam(j2).and.ireno(j2).eq.irrid(ii))idd1=j2
             if(dd2.eq.attnam(j2).and.ireno(j2).eq.irrid(ii))idd2=j2
             if(dd3.eq.attnam(j2).and.ireno(j2).eq.irrid(ii))idd3=j2
             if(dd4.eq.attnam(j2).and.ireno(j2).eq.irrid(ii))idd4=j2
            enddo
            if(icomment.eq.1)then
            print *,idd1,idd2,idd3,idd4
            endif

            i1=idd1
            i2=idd2
            i3=idd3
            i4=idd4
            j3=j3+1
            ddih(j3)=ddihedr(i1,i2,i3,i4)
          enddo
         endif
        enddo
        close(62)

        enddo
        endif

c-----------------------------------------------------------------------
c-------find partners for flexible residues-----------------------------        
        if(icomment.eq.1)print *,'finding partner residues',matp
c-------read protein----------------------------------------------------
        open(unit=1,file='prot.pdb',status='old')
        do i=1,matp
        read(1,191)pratnam(i),presnam(i),pchid(i),presid(i),
     &            (px(i,l),l=1,3)
        enddo
        close(1)

        open(unit=2,file='flexres',status='old')
        do j=1,iline
        read(2,191)rratnam(j),rresnam(j),rchid(j),rresid(j),
     &            (rx(j,l),l=1,3)
        enddo
        close(2)
c----------------------------------------------------------------------
        jres=0

        if(icomment.eq.1) print *,'iline-recflex',iline
        if(icomment.eq.1) print *,'matp-recflex',matp

        do i=1,iline
        do j=1,matp

        ax = rx(i,1)
        ay = rx(i,2)
        az = rx(i,3)
        bx = px(j,1)
        by = px(j,2)
        bz = px(j,3)
        rres=rresid(i)
        pres=presid(j)
        cutd=dist(ax,ay,az,bx,by,bz)
c--scan for residues between 2.2A and 6.0A around the selected residue --
        if(cutd.le.4.0.and.cutd.ge.0.5)then
        jres=jres+1
        sres(jres)=pres
        endif

        enddo
        enddo
c-picking id of partner residues begins here---------------------------
        l=0

        do i=1,jres !total number of atoms
        ik=0
        if(i.eq.1)then
         l=l+1
         ssres(l)=sres(i)!first residue id written into ssres(1)
        else
          do j=1,l
           if(sres(i).ne.ssres(j))then
            ik=ik+1
           endif
          enddo

          if(ik.eq.l)then
            l=l+1
            ssres(l)=sres(i)
           if(icomment.eq.1) print *,'ssres',ssres(l)
          endif
        endif
        enddo

       if(icomment.eq.1)print *,'partner residues',l
        if(l.eq.0)then
        print *,'flexible residue(s) not in the binding site'
        stop
        endif
c-------picking id of partner residues ends here------------------------
c------sort the partner residues in ascending order---------------------
        do i=1,l-1
         do j=i+1,l
          if(ssres(i)>ssres(j))then
           temp=ssres(i)
           ssres(i)=ssres(j)
           ssres(j)=temp
          endif
         enddo
        enddo
        open(unit=110,file='sortres',status='unknown')
        do i=1,l
         write(110,30)ssres(i)
        enddo
        close(110)
c------sorting ends here-----------------------------------------------
        rmpartner='rm apartner'
        call system(rmpartner)
        rmpartner='rm partner'
        call system(rmpartner)
        rmpartner='rm partners'
        call system(rmpartner)
c----------------------------------------------------------------------
        ipart=0
        open(unit=277,file='apartner',status='unknown')
        do il=1,l !l=number of partners

        open(unit=177,file='rec.pdb',status='unknown')
        iml=ssres(il)
        if(icomment.eq.1)print *,'apartner-iml',iml

        do imn=1,ipatom
        read(177,190)pl1,ml,pl2
        if(ml.eq.iml)then
        ipart=ipart+1
        write(277,190)pl1,ml,pl2
        endif
        enddo
        close(177)

        enddo
        close(277)
        if(icomment.eq.1) print *,'ipart-recflex',ipart
	if(ipart.ge.maxatm) STOP 'Partner atoms greater than maxatm'
c-------check number of chain IDs in partner-------------------------
        ncid=0
        open(unit=277,file='apartner',status='old')
        do i = 1,ipart
         read(277,192) sstr(i)
         ikc=0
         if(i.eq.1)then
         ncid=ncid+1
         castr=sstr(i)
         chid(ncid)=castr(22:22)
         else
         castr=sstr(i)
          do j=1,ncid
           cbstr=chid(j)
           if(castr(22:22).ne.cbstr(1:1))then
           ikc=ikc+1
           endif
          enddo
          if(ikc.eq.ncid)then
          ncid=ncid+1
          chid(ncid)=castr(22:22)
          endif
         endif
        enddo
        close(277)
        if(icomment.eq.1)print *,'chains : ',ncid
c-------write partner chain wise-------------------------------------
        open(unit=278,file='partner',status='unknown')
        do i=1,ncid
         cbstr=chid(i)
         open(unit=277,file='apartner',status='old')
         do j=1,ipart
         read(277,192)sstr(j)
         castr=sstr(j)
         if(castr(22:22).eq.cbstr(1:1))then
         write(278,192)sstr(j)
         endif
         enddo
         close(277)
        enddo
        close(278)
c------------------------------------------------------------------------
        rewind 278
        open(unit=278,file='partner',status='old')
        do i=1,ipart
        read(278,221)qatname(i),qresname(i),qcid(i),qresnum(i),
     &  qx1(i),qx2(i),qx3(i)
        enddo
c---------------------------------------------------------------------
c       'famppar' will parameterise 'partners' for AMBER.....
c...... energy calculation and writes parameters in 'inp.par'
        if(icomment.eq.1)print *,'parameterizing protein residues'
        call famppar()
c       'pdb' stores the x,y,z of 'partners' into (x(i,j),j=1,3)
        call pdb()
c       'input' takes parameters from 'inp.par' and stores......
c...... in (x(i,j),j=4,8)
        call finput()
        return
        end
c**********************************************************************
c----------------------------------------------------------------------
        function ddihedr(i1,i2,i3,i4)
        include 'mols.par'
        common /rcrda/xx_one(maxatm,3)
        common /comment/icomment
        real x(maxatm,3)
        acdc=((180.0*7.0)/22.0)
        one=1.d0
!---------------------------------------        
        x(i1,1)=xx_one(i1,1)
        x(i1,2)=xx_one(i1,2)
        x(i1,3)=xx_one(i1,3)
        x(i2,1)=xx_one(i2,1)
        x(i2,2)=xx_one(i2,2)
        x(i2,3)=xx_one(i2,3)
        x(i3,1)=xx_one(i3,1)
        x(i3,2)=xx_one(i3,2)
        x(i3,3)=xx_one(i3,3)
        x(i4,1)=xx_one(i4,1)
        x(i4,2)=xx_one(i4,2)
        x(i4,3)=xx_one(i4,3)
        if(icomment.eq.1)then
        print *,x(i1,1),x(i1,2),x(i1,3)
        print *,x(i2,1),x(i2,2),x(i2,3)
        print *,x(i3,1),x(i3,2),x(i3,3)
        print *,x(i4,1),x(i4,2),x(i4,3)
        endif
!--------------------------------------
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
          ddihedr=acos(a)*acdc
          if (ux1*(uy2*z2-uz2*y2)+uy1*(uz2*x2-ux2*z2)+
     #        uz1*(ux2*y2-uy2*x2).lt.zero) ddihedr =-ddihedr
          return
        else
        write (*,'(a,4i5)')' ddihedr> Error in coordinates of atoms #: '
     #                     ,i1,i2,i3,i4
        stop
      endif
      end
c------------------------------------------------------------------------- 

