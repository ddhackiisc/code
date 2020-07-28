c-----------------------------------------------------------------------------
c     library name : flexall.f   

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
        subroutine flexall(gx,gy,phi,e_flex,gz)
        include 'mols.par'
        parameter(dx=50,fx=2000)

	common /comment/icomment
	common /rctl/iscopt
	common /part/ipart
	common /ddihed/ddih(maxatm)
	common /native/nat,pepfile
	common /recsc/iline,idihed,ivx(maxatm,4)
	common /par/natom, ntor, nhb, ns, lres
        common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
        common /propdb/ihflex,patname(maxatm),px1(maxatm),px2(maxatm),
     &  px3(maxatm),presname(maxatm),pcid(maxatm),presnum(maxatm)

50      format(30x,f8.3,f8.3,f8.3)
120     format(a7,i4,1x,a4,a14,3f8.3,a24)
130     format(12x,a4,1x,a3,1x,a1,i4,4x,3f8.3,24x)

	character*80 str
	character*24 prt1,atm3
        character*4 latname(maxatm),patname,atm1
	character *3 lresname(maxatm),presname
        character*1 lcid(maxatm),pcid
	character prt*12,atm2*14,atm*7

        dimension phi(maxpar)
        real xx_two(maxatm,3),yy(maxatm,3),
     &       yy_old(maxatm,3),ex1,ex2,ex3,fene,ene,e_flex,
     &       en_flex(maxord,maxord),ddih,lx1(maxatm),
     &       lx2(maxatm),lx3(maxatm),px1,
     &       px2,px3,
     &       xx_one(maxatm,3)

        integer gz,exresid,atmno,lresnum(maxatm),presnum,ihflex,
     &  gx,gy,nat 
c------------------------------------------------------------------
        e_flex=0.0
        fene=0.0
        is=ntor+6
 	if(nat.eq.1)goto 314
 	if(gz.eq.2.and.icomment.eq.1)then
        do ik=is+1,is+np
        print *,'mini-phi-flexall',phi(ik)    
        enddo
        print *,'iline',iline
        endif
c------------------------------------------------------------------
	open(unit=1,file='flexres',status='old')
	do il=1,iline
	read(1,50),(xx_one(il,ki),ki=1,3)
	if(icomment.eq.1.and.gz.eq.2)print *,'xx_one',xx_one(il,1),
     &  xx_one(il,2),xx_one(il,3)
	enddo
	close(1)
c------------------------------------------------------------------
        do ip=1,idihed        
        is=is+1
c##################################################################
        call elemen(xx_one(ivx(ip,1),1),xx_one(ivx(ip,1),2),
     &              xx_one(ivx(ip,1),3),
     &              xx_one(ivx(ip,2),1),xx_one(ivx(ip,2),2),
     &              xx_one(ivx(ip,2),3),
     &              el,em,en)


        if(iscopt.eq.1)then
        if(icomment.eq.1)then
        print *,'rotamer-phi',is,phi(is)
        endif
        endif

        do k=1,ivx(ip,3)-1
           do ki=1,3
             xx_two(k,ki)=xx_one(k,ki)
           enddo
        enddo

        IF(iscopt.eq.1.and.gz.eq.1)then
        if(ddih(ip).lt.0.0)then
        ddih(ip)=ddih(ip)+360
        endif

        if(ddih(ip).lt.phi(is))then
        phi(is)=phi(is)-ddih(ip)
        if(icomment.eq.1)print *,'crotamer-phi',is,phi(is)
        else
        phi(is)=(360-ddih(ip))+phi(is)
        if(icomment.eq.1)print *,'crotamer-phi',is,phi(is)
        endif
 
        ENDIF

        do k=ivx(ip,3),ivx(ip,4)
           xxin=xx_one(k,1)-xx_one(ivx(ip,1),1)
           yyin=xx_one(k,2)-xx_one(ivx(ip,1),2)
           zzin=xx_one(k,3)-xx_one(ivx(ip,1),3)

           call rotor(el,em,en,phi(is),xxin,yyin,zzin,
     &                xxout,yyout,zzout)

           xx_two(k,1)=xxout+xx_one(ivx(ip,1),1)
           xx_two(k,2)=yyout+xx_one(ivx(ip,1),2)
           xx_two(k,3)=zzout+xx_one(ivx(ip,1),3)
        enddo

        do k=ivx(ip,4)+1,iline
           do ki=1,3
             xx_two(k,ki)=xx_one(k,ki)
           enddo
        enddo

        do k=1,iline
           do ki=1,3
              xx_one(k,ki)=xx_two(k,ki)
           enddo
        enddo

c#####################################################################
        enddo
        close(9)

         do k=1,iline
           do ki=1,3
             yy(k,ki)=xx_two(k,ki)
           enddo
        if(icomment.eq.1.and.gz.eq.2)print *,'rotated',yy(k,1),
     &  yy(k,2),yy(k,3)
         enddo
c-----------------------------------------------------------------------
314     ihflex=0
	if(nat.eq.1)then
	open(unit=8,file='flexres',status='old')
        open(unit=9,file='hflexres_new.pdb',!Access='append',
     &     status='unknown')
        do k=1,iline
         read(8,120)atm,atmno,atm1,atm2,(yy_old(k,ki),ki=1,3),atm3
           ihflex=ihflex+1 
             write(9,120)atm,atmno,atm1,atm2,(yy_old(k,ki),ki=1,3),atm3
        enddo
        close(9)
        close(8)
      else
      open(unit=8,file='flexres',status='old')
      open(unit=9,file='hflexres_new.pdb',!Access='append',
     &     status='unknown')
        do k=1,iline
         read(8,120)atm,atmno,atm1,atm2,(yy_old(k,ki),ki=1,3),atm3
          if(atmno.ne.0)then
          ihflex=ihflex+1 
            if(atm1.eq.' N  '.or.atm1.eq.' CA '.or.atm1.eq.' C  '.or.
     &         atm1.eq.' HA '.or.
     &         atm1.eq.' H  '.or.atm1.eq.' O  ')then
             write(9,120)atm,atmno,atm1,atm2,(yy_old(k,ki),ki=1,3),atm3
            else
                write(9,120)atm,atmno,atm1,atm2,(yy(k,ki),ki=1,3),atm3
            endif
          endif
        enddo
        close(9)
        close(8)
	endif
c-----------------------------------------------------------------------
        open(unit=9,file='hflexres_new.pdb',status='old')
        do i=1,ihflex
        read(9,130)patname(i),presname(i),pcid(i),presnum(i),
     &  px1(i),px2(i),px3(i)
        enddo
        close(9)
c-----------------------------------------------------------------------
        en_flex(gx,gy)=(afampene(gx,gy)*nres)
        e_flex=en_flex(gx,gy)
        if(icomment.eq.1)print *,'nres',nres,'e_flex',e_flex
c-----------------------------------------------------------------------
        return
        end

        
