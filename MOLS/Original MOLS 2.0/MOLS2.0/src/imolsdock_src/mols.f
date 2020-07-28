c-----------------------------------------------------------------------------
c     library name : mols.f   

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

c       program mols
	subroutine mols(inp)
	include 'mols.par'
	parameter(mxtyat = 18)

	common /comment/icomment
	common /calls/ncalls
	common /rctl/iscopt
	common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /pepcrda/xpep(maxatm,8)
        common /pepcrdb/ypep(maxatm,8)
	common /vectors/iv(maxpar,4)
	common /energy/e_final(maxord,maxord)
	common /optmo/opmophi(maxpar,maxpar)
	common /pdbat/atom(maxatm),ele(maxatm)
	common /mean/avrg1(maxpar,maxord)
	common /freq/ifreq(maxpar,maxord,maxord)
	common /order/nn,mm
	common /scang/frange,rang
        common /native/nat,pepfile
	common /cen/bx,by,bz
	common /samfnam/sf0,sf1,sf2
	common /par/natom,ntor,nhb,ns,lres
	common /out/e_out(maxatm),p_out(maxpar,maxatm)
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
 	common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
	common /emols/ligene(maxord),plpene(maxord),proene(maxord)
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
	common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
	common /scc/sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
        common /sscc/ssc(maxatm,maxpar,5),nsrot(maxpar),
     &  nschi(maxpar),mn(maxpar)       
        common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
	common /left/le(maxi,maxi),innd(maxi),ielenum(maxi),
     &  bonum(maxi)
	common /getrot/inrot,irotb(maxi,2),el,ilsrot(maxi,maxi),
     &  iatrot(maxi),rx(100,3),ire(maxi,maxi),ind(maxi,maxi),big,
     &  frotb(maxi,2)

        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

        dimension angl(maxpar,maxord),angle(maxpar,maxord),jseed(5000),
     &  ang(maxpar,maxord)

        integer fop
        real eflex,en_flex(maxord,maxord),enplp(maxord,maxord),
     &  en_final(maxord,maxord),een_final,plpe,rang,et
    
778     format(i4,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3)

        eflex=0.0
        fop=ifopt
        mm=37 ! Size of MOLS dimensions
	nn = ntor
	npar = ntor + 6  + np !np - protein SC torsion angles
        if(icomment.eq.1) print *,'npar-smmols',npar        
	if(icomment.eq.1)print *,'par-peptide',ntor,'par-protein',np,
     & 'total',npar
	if(npar.gt.200)then
	print *,'total parameters exceeding 200'
	STOP
	endif
c-------when 'npar' exceeds 36 MOLS squares will be increased ----------
        if(npar.gt.36)then
        icnt=0
        iy=npar+1

        DO
        do i=3,iy
          if(mod(iy,i)==0)then
          icnt=icnt+1
         endif
        enddo

         if(icnt.eq.1)exit
         iy=iy+1
         icnt=0
        ENDDO

        mm=iy
        endif
c------------------------------------------------------------------------
        rang=frange/mm !receptor side-chain fluctuation step size
	if(icomment.eq.1)print *,'mm',mm,'frange',frange,'rang',rang
c------------------------------------------------------------------------
        if(icomment.eq.1)print *,'isc-smmols',iscopt
c------------------------------------------------------------------------
        if(ilopt.eq.1)then
        call pprecal
        else
	call precal
        endif

	call anggen(ang,fseq,npar,iopt,aii)

	if(iopt.eq.3) call scinp(fseq)

	IF(ilopt.eq.1)THEN
        do i=1,natom
         do j=1,8
           ypep(i,j)=xpep(i,j)
         enddo
        enddo
        ELSE IF(ilopt.eq.2)THEN
        do i=1,natom
         do j=1,8
          ylig(i,j)=xlig(i,j)
         enddo          
	enddo   
        ENDIF
c------------------------------------------------------------
        do i=1,npar
	 do j = 1,mm
          angl(i,j) =ang(i,j)
	 enddo
        enddo
c************************************************************
	call rand1(jseed)
        nt=1
        do ijx=inp,inp!1,ns
	iseed = jseed(ijx)

        call write_par(angl,iseed,angle,npar)
    
        call pargen1(angle,npar)

        do i=1,100
          e_out(i)=1.0e+25
        enddo

        do i=1,npar
          do j=1,mm
            do k=1,mm
              ifreq(i,j,k)=0
            enddo
          enddo
        enddo
c      ************** Main MOLS ***************
         do l1=1,npar
         do l2=1,mm
            p(l1,l2,2)=0.0
            p(l1,l2,3)=0.0
         enddo
         enddo
        
        do i=1,nres
        mn(i)=1
        enddo
    
           do i=1,mm
           do j=1,mm
	     k = 1
	
	    call subpar(i,j,k,aii,phi,npar)

            if(ilopt.eq.1)then
            call pmolgen(i,j,phi,1,eflex)
            else if(ilopt.eq.2)then
	    call molgen(i,j,phi,1,eflex)
            endif
c---------------------------------------------------------------------
         if(ilopt.eq.1)then
          if(iff.eq.1)then
           en_final(i,j)=pampene(i,j)
          else
           en_final(i,j)=pecpene(i,j)
          endif
         else if(ilopt.eq.2)then
          if(iff.eq.1)then
           rmmff=rmmffene(i,j)     
           en_final(i,j)=rmmff
          else if(iff.eq.2)then
           rgaff=rgaffene(i,j)
           en_final(i,j)=rgaff
          endif
         endif
        
        if(ilopt.eq.1)then
        call peplp(plpe,hb,steric,rep,1)
        enplp(i,j)=plpe
        else if(ilopt.eq.2)then
        call eplp(plpe,hb,steric,rep,1)
        enplp(i,j)=plpe
        endif
        en_flex(i,j)=eflex
c----------------------------------------------------------------------
        print *,inp,i,j
        if(icomment.gt.0)then
        print *,'intra-ligand',en_final(i,j)
        print *,'inter-ene',enplp(i,j)
        print *,'intra-protein',en_flex(i,j)
        endif
        
        e_final(i,j) = enplp(i,j)+ en_final(i,j)+en_flex(i,j)!energy along with plp
        if(icomment.eq.1) print *,'total',e_final(i,j)
	if(nat.eq.1)then
	print *,'native-peptide',en_final(i,j)
	print *,'native-plp',enplp(i,j)
	print *,'native-protein',en_flex(i,j)
	print *,'native-total',e_final(i,j)
	stop
	endif
c----------------------------------------------------------------------
           enddo
          enddo          
       
          call average(npar)

          do i=1,npar
            call sort_and_set_rank(i)
          enddo
!--------------------------------------------
          do i=1,npar 
           call rank_sort(i,nt)
          enddo
!--------------------------------------------
c       ************ LAST MOLS **************
	do m1=1,15
          m2 = 2
	  kk = 1

	  call subpar(kk,m1,m2,aii,phi,npar)

          if(ilopt.eq.1)then
          call pmolgen(kk,m1,phi,1,eflex)
          else if(ilopt.eq.2)then
	  call molgen(kk,m1,phi,1,eflex)
          endif 
	  
cs------PHI of all structures of best 15------------
        do ik4=1,npar
        opmophi(m1,ik4)=phi(ik4)
        if(icomment.eq.1)print *,'opmophi',m1,opmophi(m1,ik4)
        enddo
cs--------------------------------------------------
        if(ilopt.eq.1)then
          if(iff.eq.1)then
           en_final(m1,m1)=pampene(m1,m1)
          else
           en_final(m1,m1)=pecpene(m1,m1)
          endif
        else if(ilopt.eq.2)then 
          if(iff.eq.1)then
            rmmff=rmmffene(m1,m1)
            en_final(m1,m1)=rmmff
          else if(iff.eq.2)then
            rgaff=rgaffene(m1,m1) 
            en_final(m1,m1)=rgaff
          endif
        endif

        if(ilopt.eq.1)then
        call peplp(plpe,hb,steric,rep,1)
        enplp(m1,m1)=plpe
        else
        call eplp(plpe,hb,steric,rep,1)
        enplp(m1,m1)=plpe
        endif
        
        en_flex(m1,m1)=eflex

	print *,'MOLS energy, PLP,AMBER, final energies: ',en_final(m1,m1),
     &	enplp(m1,m1),en_flex(m1,m1),
     &  en_final(m1,m1)+enplp(m1,m1)+en_flex(m1,m1)

        if(icomment.eq.1)then
        print *,'best_15'
        print *,kk,m1
        print *,'intra-lene',en_final(m1,m1)
        print *,'plp',enplp(m1,m1)
        print *,'intra-pene',en_flex(m1,m1),eflex
        endif
        
        e_final(m1,m1) = enplp(m1,m1)+en_final(m1,m1)
     &  + en_flex(m1,m1)!energy along with plp&flex

	call best(m1,npar)
	enddo

        call output(i,j,1,ijx,et,aii,npar) !included npar in arguments

	write(*,*)'MOLS Structure No      :',ijx,et
c*******Flushing standard output************
	call flush(6)
c*******End of Flushing ********************
	write(31,FMT='(5x,a23,2x,i4,2x,f10.3)')'MOLS Optimal Structure:',ijx,et
c-------mols.log------------------------------------------------------
	do iil=1,15
        if(e_final(iil,iil).eq.et)then
        ir=iil       
        endif
        enddo

        if(icomment.eq.1) print *,'ir-mols-best',ir
        if(icomment.eq.1)write(*,778)ijx,en_final(ir,ir),enplp(ir,ir),
     &   en_flex(ir,ir),e_final(ir,ir)
c--------------------------------------------------------------------
333	enddo

	return
        end
c********************************************************************
        subroutine write_par(angl,iseed,angle,npar)

	include 'mols.par'

        common/order/nn,mm

        dimension angll(maxpar,maxord),angle(maxpar,maxord),
     &  angl(maxpar,maxord)

c	open(unit=4,file='molspargen.inp',status='unknown')

        do i=1,npar
         do j=1,mm
          angll(i,j)=angl(i,j)
         enddo
        enddo
        amm=float(mm)
        do j=1,npar
         ki=0
         do i=1,100000
          y=rand(iseed)
          iseed = y
          iy=nint(y*amm)
          if(iy.eq.0)iy=1
          if (angll(j,iy).lt.999) then
           ki=ki+1
           angle(j,ki)=angll(j,iy)
           angll(j,iy)=99999.0
           if(ki.ge.mm)go to 3 
          endif         
         enddo         
3       continue
        enddo 

c       do i=1,npar
c       write(4,10)(angle(i,j),j=1,mm)
c	write(4,*)(angle(i,j),j=1,mm)
c10     format(37f6.1)
c       enddo
    
        close(unit=4)
        return 
        end
c*********************************************************************
	subroutine pargen1(angle,npar)

c	pargen is for parameter generation
c	n  -   total number of parameters
c	m  -   maximum number of values for each parameter
      
	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/order/nn,mm
	dimension angle(maxpar,maxord)
c	open (unit=5,file='molspargen.inp',status='unknown')

	do 150 i = 1, npar
	 do 150 j = 1,mm
	  p(i,j,1) = angle(i,j)
150     continue
	do 500 l = 1, npar
	do 500 k = 1, mm
	do 500 j = 1, mm
	 i = mod( ((j-1) * (l-1) + (k-1)), mm) + 1
	 e(i,j,l)  = p(l,k,1)
500	 continue
        close (unit=5)

c20	format(8f6.1)
c	open (unit=6,file='pargen.out',status='unknown')
c	write (6,20) (((e(i,j,l),l=1,npar),j=1,mm),i=1,mm)
c	close (unit=6)        

	return
	end
c*********************************************************************
c-------Peptide Docking-----------------------------------------------
        subroutine pmolgen(lx,ly,phi,tst,e_flex)

        include 'mols.par'

	common /pepcrda/xpep(maxatm,8)
        common /pepcrdb/ypep(maxatm,8)
	common /vectors/iv(maxpar,4)
	common /cen/bx,by,bz
	common /par/natom,ntor,nhb,ns,lres
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)

        dimension x_one(maxatm,3),x_two(maxatm,3),phi(maxpar)
	real cx,cy,cz
        integer tst,nn,nat

        IF(nat.eq.1)THEN !nat=1:Peptide-Protein Native Energy
	do k=1,natom
         do ki=1,3
          ypep(k,ki)=xpep(k,ki)
         enddo
        enddo     
        if(ifopt.eq.2) call flexall(lx,ly,phi,e_flex,1)
	goto 313
	ENDIF

        nn = ntor

        do k=1,natom
         do ki=1,3
          x_one(k,ki)=xpep(k,ki)
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
             ypep(k,ki)=x_two(k,ki)
           enddo
        enddo
c***********translation and rotation for ligand postioning***********
        call protate(bx,by,bz,phi(nn+4),phi(nn+5),phi(nn+6),natom)
        call ptranslate(bx,by,bz,phi(nn+1),phi(nn+2),phi(nn+3),natom)
        if(ifopt.eq.2) call flexall(lx,ly,phi,e_flex,1)
c---------------------------------------------------------------------
313     return
        end
c********************************************************************* 
c-------Protein-Ligand Docking----------------------------------------
        subroutine molgen(lx,ly,phi,tst,e_flex)

        include 'mols.par'

	common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
	common /vectors/iv(maxpar,4)
	common /order/nn,mm
	common /cen/bx,by,bz
	common /par/natom,ntor,nhb,ns,lres
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
	common /left/le(maxatm,maxatm),innd(maxatm),ielenum(maxatm),
     &  bonum(maxatm)
        common /getrot/inrot,irotb(maxatm,2),el,ilsrot(maxatm,maxatm),
     &  iatrot(maxatm),rx(100,3),ire(maxatm,maxatm),
     &  ind(maxatm,maxatm),big

        dimension x_one(maxatm,3),x_two(maxatm,3),phi(maxpar)
	real cx,cy,cz,bx,by,bz
        integer tst,nn,ci
 
        IF(nat.eq.1)THEN 
	 do k=1,natom
          do ki=1,3
           ylig(k,ki)=rx(k,ki)
          enddo
         enddo

         if(ifopt.eq.2) call flexall(lx,ly,phi,e_flex,1)
	goto 314
	ENDIF

	nn=ntor

        do k=1,natom
           do ki=1,3
             x_one(k,ki)=rx(k,ki)
           enddo
        enddo
c-------if rotatable bonds in ligand is zero------------------------
        IF(ntor.eq.0)THEN
        do k=1,natom
           do ki=1,3
             ylig(k,ki)=rx(k,ki)
           enddo
        enddo

        call rotate(bx,by,bz,phi(nn+4),phi(nn+5),phi(nn+6),natom)
        call translate(bx,by,bz,phi(nn+1),phi(nn+2),phi(nn+3),natom)
        if(ifopt.eq.2)then
        call flexall(lx,ly,phi,e_flex,1)
        endif
        goto 314
        ENDIF
c-------------------------------------------------------------------
        do if=1,nn
C###################### PHI ALL ####################################

        call elemen(x_one(irotb(if,1),1),x_one(irotb(if,1),2),
     &              x_one(irotb(if,1),3),
     &              x_one(irotb(if,2),1),x_one(irotb(if,2),2),
     &              x_one(irotb(if,2),3),
     &              el,em,en)

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
        
        call rotate(bx,by,bz,phi(nn+4),phi(nn+5),phi(nn+6),natom)
	call translate(bx,by,bz,phi(nn+1),phi(nn+2),phi(nn+3),natom)
        if(ifopt.eq.2) call flexall(lx,ly,phi,e_flex,1)

314     return
        end
c***********************************************************************
        subroutine average(npar)

	include 'mols.par'

	common /mean/avrg1(maxpar,maxord)
	common /energy/e_final(maxord,maxord)
	common /order/nn,mm
	common /out/e_out(maxatm),p_out(maxpar,maxatm)
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)

        data en_k_t/1.98578e+03/

c       en_k_t is boltzmannconst*avogadronumber*temp(=300 degree K)
c	average is for finding the average 
c	e_final   -  energy values of the conformation 
c	cutoff    -  for omitting impossible conformations 
c	avrg1     -  average without cutoff
c	n  -   total number of parameters = 9
c	m  -   maximum number of values for each parameter = 37

c	compute using mols

	do 400 l = 1, npar
	  do 500 k = 1, mm
	    avrg1(l,k) = 0.0 
            summ   = 0.0
            summ_mm= 0.0
	    do 600 j = 1, mm
	       i = mod( ( (j-1) * (l-1) + (k-1) ), mm) + 1
               weight=exp(-1.0*(e_final(i,j)/(en_k_t*0.01)))
  	       summ = summ + e_final(i,j)*weight
               summ_mm=summ_mm+weight
600	    continue
	    avrg1(l,k) = summ/summ_mm
500	continue
400    continue

        return
        end
c*********************************************************************
        subroutine sort_and_set_rank(is1)

	include 'mols.par'

	common /order/nn,mm
	common /mean/avrg1(maxpar,maxord)
	common /freq/ifreq(maxpar,maxord,maxord)
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)

        dimension rank(maxord),icc(maxord)

        do j=2,mm
          aa=avrg1(is1,j)
          bb=p(is1,j,1)
          do i=j-1,1,-1
            if(avrg1(is1,i).le.aa) go to 10
            avrg1(is1,i+1)=avrg1(is1,i)
            p(is1,i+1,1)=p(is1,i,1)
          enddo
          i=0
10        avrg1(is1,i+1)=aa
          p(is1,i+1,1)=bb
        enddo

        do i=1,mm
          rank(i)=float(i)
        enddo

        do j=2,mm
          aa=rank(j)
          bb=p(is1,j,1)
          do i=j-1,1,-1
            if(p(is1,i,1).le.bb) go to 20
            p(is1,i+1,1)=p(is1,i,1)
            rank(i+1)=rank(i)
          enddo
          i=0
20        rank(i+1)=aa
          p(is1,i+1,1)=bb
        enddo

        do i=1,mm
          jxx=ifix(rank(i))
          ifreq(is1,i,jxx)=ifreq(is1,i,jxx)+1
          p(is1,i,2)=p(is1,i,2)+rank(i)
          p(is1,i,3)=p(is1,i,3)+rank(i)*rank(i)
        enddo

        return
        end 
c********************************************************************
        subroutine rank_sort(ir1,nt)

	include 'mols.par'

	common/order/nn,mm
	common/freq/ifreq(maxpar,maxord,maxord)
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)

        dimension ee(maxord)

        do i=1,mm
          sum_x=p(ir1,i,2)
          sum_x2=p(ir1,i,3)
          ave=sum_x/nt
          sd=sqrt((sum_x2/nt)-((sum_x/nt)*(sum_x/nt)))
          p(ir1,i,2)=ave
          p(ir1,i,3)=sd
        enddo       

       do j=2,mm
          bb=p(ir1,j,1)
          cc=p(ir1,j,2)
          dd=p(ir1,j,3)
          do j1=1,mm
            ee(j1)=ifreq(ir1,j,j1)
          enddo
          do i=j-1,1,-1
            if(p(ir1,i,2).le.cc) go to 10
            p(ir1,i+1,3)=p(ir1,i,3)
            p(ir1,i+1,2)=p(ir1,i,2)
            p(ir1,i+1,1)=p(ir1,i,1)
            do j1=1,mm
               ifreq(ir1,i+1,j1)=ifreq(ir1,i,j1)
            enddo
          enddo
          i=0
10        p(ir1,i+1,1)=bb
          p(ir1,i+1,2)=cc
          p(ir1,i+1,3)=dd
          do j1=1,mm
            ifreq(ir1,i+1,j1)=ee(j1)
          enddo
        enddo

	do j=2,mm
	 bb=p(ir1,j,1)
	 cc=p(ir1,j,2)
	 dd=p(ir1,j,3)
	 do j1=1,mm
         ee(j1)=ifreq(ir1,j,j1)
	 enddo
	  do i=j-1,1,-1
	   if (p(ir1,i,2).ne.cc) goto 40
	   if (p(ir1,i,3).le.dd) goto 40
	   p(ir1,i+1,3)=p(ir1,i,3)
	   p(ir1,i+1,2)=p(ir1,i,2)
	   p(ir1,i+1,1)=p(ir1,i,1)
	   do j1=1,mm
	   ifreq(ir1,i+1,j1)=ifreq(ir1,i,j1)
	   enddo
	  enddo
	  i=0
40 	  p(ir1,i+1,1)=bb
	  p(ir1,i+1,2)=cc
	  p(ir1,i+1,3)=dd
	  do j1=1,mm
	  ifreq(ir1,i+1,j1)=ee(j1)
	  enddo
	enddo

        return
        end
c*********************************************************************
       subroutine output(cx,cy,nou,ijx,et,aii,npar)

	include 'mols.par'

	common /comment/icomment
	common /optmo/opmophi(maxpar,maxpar)
	common /energy/e_final(maxord,maxord)
	common /order/nn,mm
        common /scang/frange,rang
        common /out/e_out(maxatm),p_out(maxpar,maxatm)
        common /par/natom,ntor,nhb,ns,lres,is,ie
	common /scc/sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),
     &  irrid(maxatm)
        common /sscc/ssc(maxatm,maxpar,5),nsrot(maxpar),
     &  nschi(maxpar)
	common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
        common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt 

        character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

	character*80 prt,prt1 
        integer exresid,cx,cy,npar
        real ex1,ex2,ex3,et

	dimension pp(maxpar)

        npar = ntor+6+np

        i = nou

        if(icomment.eq.1)write(*,*),'output-ifopt',ifopt
        if(icomment.eq.1)write(*,*),'rang-output',rang        

c----------------------------------------------------------------
	open(unit=24,file=of1,status='unknown')
c----------------------------------------------------------------
	if(iopt.ne.3) then
	do k1=1,nn
	opmo(ijx,k1) = p_out(k1,i)
        if(icomment.eq.1)print *,'opmo',ijx,k1,opmo(ijx,k1)
	enddo 
	emo(ijx) = e_out(i)
	do k1 = 1,nn
        p_out(k1,i) = p_out(k1,i) - 180.0

	enddo
c************
        do k2 = nn+1,nn+3
        tt1 = p_out(k2,i)
        opmo(ijx,k2) = tt1
        if(icomment.eq.1)print *,'opmo-nn+1,nn+3',ijx,k2,opmo(ijx,k2)
        if(tt1.eq.0.0) then
        tt2 = 0.0
        else
        tt2 = (tt1/10.0)*0.13888889
        endif
        tt2 = tt2-2.5
        p_out(k2,i)=tt2
        enddo
        do k2=nn+4,nn+6
        opmo(ijx,k2) = p_out(k2,i)
        if(icomment.eq.1)print *,'opmo-nn+4,nn+6',ijx,k2,opmo(ijx,k2)
        enddo
        p_out(nn+6,i)=p_out(nn+6,i)/2
c-------opmo for receptor--sam added----------------------------------
	if(ifopt.eq.2)then
        do iil=1,15
        if(e_final(iil,iil).eq.e_out(i))then
        ir=iil
        endif
        enddo

        do ik4=nn+7,npar
        opmo(ijx,ik4)=opmophi(ir,ik4)
        if(icomment.eq.1)print *,'opmo-opmophi',ir,opmo(ijx,ik4)        
        enddo
	endif
c---------------------------------------------------------------------
        write(24,*) ijx,e_out(i),(p_out(j,i),j=1,npar)
        write(*,*)'one',ijx,e_out(i),(opmo(ijx,j),j=1,npar)
        else
        do k1=1,nn
 	opmo(ijx,k1) = pp(k1)
        if(icomment.eq.1)print *,'opmo-2',ijx,k1,opmo(ijx,k1)
	enddo
	emo(ijx) = e_out(i)
	do k1 = 1,nn
	pp(k1) = pp(k1) - 180.0
	enddo
        
C**************************************
c       STORING THE TRANS AND ROT PARAM
        ix = 5
        iy = 1

281     tt1 = p_out(npar-ix,i)

        opmo(ijx,nn+iy) = tt1
        if(icomment.eq.1)print *,'opmo-3',ijx,nn+iy,opmo(ijx,nn+iy)

        if(tt1.eq.0.0) then
        tt2 = 0.0
        else
        tt2 = (tt1/10.0)*0.13888889
        endif
        tt2 = tt2-2.5
        pp(nn+iy) = tt2
        ix = ix-1
        iy = iy+1
        if(iy.le.3) go to 281
        ix = 2
        iy = 4
282     pp(nn+iy) = p_out(npar-ix,i)
        opmo(ijx,nn+iy) = p_out(npar-ix,i)
        if(icomment.eq.1)print *,'opmo-4',ijx,nn+iy,opmo(ijx,nn+iy)
        ix = ix-1
        iy = iy+1
        if(iy.le.6) go to 282

        pp(nn+6) = pp(nn+6)/2
C***********************************************************************
        write(24,*) ijx,e_out(i),(pp(j),j=1,npar)
        write(*,*)'two',ijx,e_out(i),(pp(j),j=1,npar)
        endif 

	et = e_out(i)
        if(ijx.eq.ie) close(unit=24)
c-----------------------------------------------------------------------
        return
        end
c***********************************************************************
	subroutine best(ib,npar)
	include 'mols.par'

	common /comment/icomment
	common /order/nn,mm
	common /energy/e_final(maxord,maxord)
	common /out/e_out(maxatm),p_out(maxpar,maxatm)
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3) 

        do i=1,25
          if(e_final(ib,ib).lt.e_out(i)) then            
            do j=25,i+1,-1
              e_out(j)=e_out(j-1)
              do k=1,npar
                p_out(k,j)=p_out(k,j-1)
              enddo
            enddo
            e_out(i)=e_final(ib,ib)
   
	    do l = 1,npar
	      p_out(l,i) = p(l,ib,1)
            if(icomment.eq.1)then 
            print *,'best-p_out',p_out(l,i)
            endif
	    enddo
            go to 10
          endif
         enddo
10       continue
         ir=ib !sam added
         
         return
         end
C###############################################################################
	subroutine rand1(jseed)

	integer a,b,c,d,e,r,jseed(5000)

	a=2346
	b=3465
	c=5421
	d=5323
	e=5000 

	nss=0
        do i=1,10000
         r = mod(((a*b)+c),d)
	 if (r.lt.0)then
	  r = r * -1
	 endif
         r=(r*1)*2+1
         if(r.le.999.or.r.gt.9999) go to 1
         nss=nss+1
         if(nss.gt.e)go to 2
	 jseed(nss) = r
1        b = r
	enddo

2       return
        end
c*************************************************************************
	subroutine anggen(ang,fseq,npar,iopt,aii)
	include 'mols.par'
        common /order/nn,mm
        common /par/natom,ntor,nhb,ns,lres

	dimension ang(maxpar,maxord)
	character fseq(maxatm)*1 
	
	if(iopt.eq.3) then
		npar = lres*3
		do i = 1,lres
		if(fseq(i).eq.'A'.or.fseq(i).eq.'G'.or.fseq(i).eq.'P') then
		npar = npar - 1
		endif
		enddo
	 npar = npar + 6
	endif
cc-------------------------------------------------------------------
cc	check the mols order it should be greater than npar and prime
cc	number
	if(mm.lt.npar) then
	do i1 = 1,100
	k1 = 2
181	if(i1/k1*k1.eq.i1) go to 191
	k1 = k1 + 1
	if(i1.le.i1/2) go to 181
	if(i1.gt.npar) go to 201
191	enddo
201	mm = i1
	i1 = 0
	endif
cc-------------------------------------------------------------------
	aii = 360.0/(mm-1)
	all = 0
	do i=1,npar
	ang(i,1)=all
        b=all
	do j=2,mm
	ang(i,j)=b+aii
        b=ang(i,j)
        enddo
	enddo
	return 
	end
c*************************************************************************
	subroutine scinp(fseq)
	include 'mols.par'
        common /par/ natom, ntor, nhb, ns, lres
	character fseq(maxatm)*1,ch*1, tcode*4, scode*1
	common /scc/ sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)

20	format(a4,a1,1x,i2,i2)
21	format(20x,5f7.1)

 	open(unit=33,file='SC.lib',status='old')
	do i = 1,lres
	rewind 33
	 do i1 = 1,234
	 read(33,20) tcode, scode, mrot, mchi
	  if(scode.eq.fseq(i).and.mrot.ne.0) then
	   do i2 = 1, mrot
	   read(33,21) (sc(i,i2,i3),i3=1,mchi)
	   enddo
	   nrot(i) = mrot
	   nchi(i) = mchi
	   go to 212
	  endif
	 enddo
212	enddo
	do i = 1,lres
	do j = 1, nrot(i)
	write(31,*)(sc(i,j,k),k=1,nchi(i))
	enddo
	enddo
	close(unit=33)
	return
	end
c************************************************************************	
	subroutine subpar(l1,l2,l3,aii,phi,npar)
	include 'mols.par'

	common /comment/icomment
	common /rctl/iscopt
	common /vectors/iv(maxpar,4)
	common /order/nn,mm
	common /scang/frange,rang
	common /par/natom,ntor,nhb,ns,lres
	common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
	common /recep/np,nres,crnam(maxatm),crcid(maxatm),irrid(maxatm)
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
	common /scc/sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
        common /sscc/ssc(maxatm,maxpar,5),nsrot(maxpar),
     &  nschi(maxpar),mn(maxpar)

	dimension phi(maxpar)
        dimension x_one(maxatm,3),x_two(maxatm,3)

        if(icomment.eq.1) print *,'iscopt',iscopt
        if(icomment.eq.1) print *,'rang-subpar',rang

	if(iopt.eq.3) then   ! this portion is skipped as iopt !=3
	i = 1
	j = 1
	k = 1
161	if(l3.eq.1) then  
	if(fseq(k).eq.'P') go to 191 
	phi(i) = e(l1,l2,j)
	i = i + 1
	j = j + 1
191	phi(i) = e(l1,l2,j)
	i = i + 1
	j = j + 1
	if(nrot(k).ne.0.) then
	j1 = e(l1,l2,j)
	j1 = int((j1/aii) + 1.0)
	j1 = mod(j1,nrot(k))
	if(j1.eq.0) j1 = nrot(k)
	do j2 = 1, nchi(k)
	 phi(i) = sc(k,j1,j2)
	 i = i + 1
	enddo
	j = j + 1
	endif
	k = k + 1
	endif

	if(l3.eq.2) then
	if(fseq(k).eq.'P') go to 192
	phi(i) = p(j,l2,1)
	i = i + 1
	j = j + 1
192	phi(i) = p(j,l2,1)
	i = i + 1
	j = j + 1
	if(nrot(k).ne.0.) then
	j1 = p(j,l2,1)
	j1 = int((j1/aii) + 1.0)
	j1 = mod(j1,nrot(k))
	if(j1.eq.0) j1 = nrot(k)
	do j2 = 1, nchi(k)
	 phi(i) = sc(k,j1,j2)
	 i = i + 1
	enddo
	j = j + 1
	endif
	k = k + 1
	endif

	if(k.le.lres) go to 161
	else
	do i = 1, nn
	if(l3.eq.1) phi(i) = e(l1,l2,i)
	if(l3.eq.2) phi(i) = p(i,l2,1)
	enddo
	endif
        
c***modified to include rotation and translation parameters for ligand docking**************
        do i = nn+1, nn+3
        if(l3.eq.1) tt1 = e(l1,l2,i)
        if(l3.eq.2) tt1 = p(i,l2,1)
        if(tt1.eq.0.0) then
        tt2 = 0.0
        else
        tt2 = (tt1/10.0)*0.13888889
        endif
        tt2 = tt2-2.5
        phi(i) = tt2
        enddo
        do i = nn+4, nn+6
        if(l3.eq.1) phi(i) = e(l1,l2,i)
        if(l3.eq.2) phi(i) = p(i,l2,1)
        enddo

        if(iopt.eq.3) then
c       change ix according to the no. of parameters added for trans and rot
        ix = 5
        iy = 1

181     if(l3.eq.1) tt1 = e(l1,l2,npar-ix)
        if(l3.eq.2) tt1 = p(npar-ix,l2,1)
        if(tt1.eq.0.0) then
        tt2 = 0.0
        else
        tt2 = (tt1/10.0)*0.13888889
        endif
        tt2 = tt2-2.5
        phi(nn+iy) = tt2
        ix = ix-1
        iy = iy+1
        if(iy.le.3) go to 181
        ix = 2
        iy = 4
182     if(l3.eq.1) phi(nn+iy) = e(l1,l2,npar-ix)
        if(l3.eq.2) phi(nn+iy) = p(npar-ix,l2,1)
        ix = ix-1
        iy = iy+1
        if(iy.le.6) go to 182
        endif 
c*******step size for receptor SC flexibility*by*sam*******

        IF(iscopt.eq.1)then

        if(ifopt.eq.2)then
        mp=1
        do i=1,nres
        mchi=nschi(i)
        mrot=nsrot(i)        
        do j=1,mchi
        if(mn(i).gt.mrot)then
        mn(i)=1
        phi(nn+6+mp)=ssc(i,mn(i),j)
        if(icomment.eq.1) print *,'phi-rotamer',nn+6+mp,phi(nn+6+mp)
        mp=mp+1
        if(icomment.eq.1)then
        print *,'nres',i
        print *,'mchi',mchi
        print *,'mrot',mrot
        endif
        else

        phi(nn+6+mp)=ssc(i,mn(i),j) 
        if(icomment.eq.1) print *,'phi-rotamer',nn+6+mp,phi(nn+6+mp)
        mp=mp+1
        if(icomment.eq.1)then
        print *,'nres',i
        print *,'mchi',mchi
        print *,'mrot',mrot
        endif
        endif
        enddo

        if(icomment.eq.1)print *,l1,l2,'mn(i)',mn(i)
        mn(i)=mn(i)+1

        enddo

        endif

        ELSE IF(iscopt.eq.0)then

        if(ifopt.eq.2)then

        do i = nn+7,npar
        if(l3.eq.1) tt1 = e(l1,l2,i)
        if(l3.eq.2) tt1 = p(i,l2,1)

         if(tt1.eq.0.0)then
          tt2=0.0
         else
 
	 tt2=(tt1/(360.0/mm))*rang
         endif
     
c--------modulo has been added two times to shuffle '+' and '-' signs----
         if(mod(l2,2).ne.0.0)then
          if(mod(i,2).eq.0.0)then
           phi(i) = -tt2
          else 
           phi(i) = tt2
          endif
         else
          if(mod(i,2).eq.0.0)then
             phi(i) = tt2
          else
           phi(i) = -tt2       
          endif
         endif
cs      if(icomment.eq.1) print *,'tt1,phi-flex',i,tt1,phi(i)
        enddo
        
        endif !if(ifopt)

        ENDIF !if(isc)
c********************************************************************
 	return
	end

