c       program mols
	subroutine mols(inp,iopt,ifopt,ilopt)
	include 'mols.par'
	parameter(mxtyat = 18)
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /pepcrda/xpep(maxatm,8)
        common /pepcrdb/ypep(maxatm,8)
        common /energy/e_final(maxord,maxord) 
        common /mean/avrg1(maxpar,maxord)
        common /vectors/iv(maxpar,4)
        common /par/ natom, ntor, nhb, ns, lres
        common /order/nn,mm
        common /freq/ifreq(maxpar,maxord,maxord)
        common /out/e_out(maxatm),p_out(maxpar,maxatm)
	common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
        common /calls/ncalls
	common /scc/ sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
        common /sscc/ ssc(maxatm,maxpar,5),nsrot(maxpar),
     &  nschi(maxpar),mn(maxpar)
	common /pdbat/atom(maxatm),ele(maxatm)
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)       
        common /pctl/ifff,ioptt
        common /rctl/iscopt
        common /fnam/ if1,pf1,pf2,pf3,pf4,pf5,of1,of2,of3,of4,
     &  of5,of6
	common /samfnam/sf0,sf1,sf2
	common /getrot/inrot,irotb(maxi,2),el,ilsrot(maxi,maxi),
     &  iatrot(maxi),rx(100,3),ire(maxi,maxi),ind(maxi,maxi),big,
     &  frotb(maxi,2)   
        common /left/ le(maxi,maxi),innd(maxi),ielenum(maxi),
     &  bonum(maxi)

	common /cen/bx,by,bz
        common /recep/np,nres
        common /emols/ligene(maxord),plpene(maxord),proene(maxord) 
        common /comment/icomment 
        common /optmo/opmophi(maxpar,maxpar)
        common /scang/frange,rang
	common /native/nat,pepfile
        character*128 if1,pf1,pf2,pf3,pf4,pf5,pf6,of1,of2,of3,of4,
     &  of0,of5,of6,sf0,sf1,sf2
      

         
        dimension angl(maxpar,maxord),angle(maxpar,maxord),jseed(5000),
     &  ang(maxpar,maxord)
	character seq*(maxres),fseq(maxatm)*1
	real plpe  !added for dock
        integer fop
        real eflex,en_flex(maxord,maxord),enplp(maxord,maxord),
     &  en_final(maxord,maxord),een_final
        real rang
    
778     format(i4,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3)
788	format(i3,1x,i3,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2)
789	format(a9)

        eflex=0.0
        fop=ifopt
c       nn=19 ! No. of torsion angles ( or search dimensions)
        mm=37 ! Size of MOLS dimensions
cn	print*,mm
c	natom=77
cc	write(31,*) 'Enter the no. of parameters'
cc	read *, nn
	nn = ntor
	npar = ntor + 6  + np !np - protein SC torsion angles
        if(icomment.eq.1) print *,'npar-smmols',npar        
	if(icomment.eq.1)print *,'par-peptide',ntor,'par-protein',np,
     & 'total',npar
	if(npar.gt.200)then
	print *,'total parameters exceeding 200'
	stop
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
c       if(icomment.eq.1)print *,'mm',mm,'rang',rang
	if(icomment.eq.1)print *,'mm',mm,'frange',frange,'rang',rang
c------------------------------------------------------------------------
        if(icomment.eq.1)print *,'isc-smmols',iscopt

c------------------------------------------------------------------------
        if(ilopt.eq.1)then
        call pprecal
        else
	call precal
        endif

        if(ilopt.eq.1)then
        call panggen(ang,fseq,npar,iopt,aii)
        else
	call anggen(ang,fseq,npar,iopt,aii)
        endif

	if(iopt.eq.3) call pscinp(fseq)
        if(ilopt.eq.2)then
        do i=1,natom
          do j=1,8
             ylig(i,j)=xlig(i,j)
          enddo          
	enddo   
        endif
c-------sam added for pmols.f--------------------------
       if(ilopt.eq.1)then

       do i=1,natom
         do j=1,8
           ypep(i,j)=xpep(i,j)
         enddo
       enddo
       endif
c-----------------------------------------------------
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

          if(ilopt.eq.1)then
          call pwrite_par(angl,iseed,angle,npar)
          else
          call write_par(angl,iseed,angle,npar)
          endif


          if(ilopt.eq.1)then
          call ppargen1(angle,npar)
          else
          call pargen1(angle,npar)
          endif

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
c         ************ Main MOLS ***************
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
	
            if(ilopt.eq.1)then
            call psubpar(i,j,k,iopt,ifopt,fseq,aii,phi,npar)
            else
	    call subpar(i,j,k,iopt,ifopt,fseq,aii,phi,npar)
            endif

            if(ilopt.eq.1)then
            call pmolgen(i,j,phi,1,ifopt,eflex)
            else
	    call molgen(i,j,phi,1,ifopt,eflex)
            endif
c---------------------------------------------------------------------
c---------------------------------------------------------------------              
         if(ilopt.eq.1)then
          if(ifff.eq.1)then
           en_final(i,j)=pampene(i,j)
          else
           en_final(i,j)=pecpene(i,j)
          endif
         else if(ilopt.eq.2)then
          if(ifff.eq.1)then
           rmmff=rmmffene(i,j)     
           en_final(i,j)=rmmff
          else if(ifff.eq.2)then
           rgaff=rgaffene(i,j)
           en_final(i,j)=rgaff
          endif
         endif
        
        if(ilopt.eq.1)then
        call peplp(plpe,hb,steric,1)
        enplp(i,j)=plpe
        else
        call eplp(plpe,hb,steric,1)
        enplp(i,j)=plpe
        endif
        en_flex(i,j)=eflex
        write(82,788)i,j,hb,steric,plpe,en_final(i,j),eflex
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
c	print *,'total',e_final(i,j)
c	if(icomment.eq.1)stop        
c       stop
c	pause
c----------------------------------------------------------------------

           enddo
          enddo          
       
!-------------------------------------------
          if(ilopt.eq.1)then
          call paverage(npar)
          else
          call average(npar)
          endif
!--------------------------------------------
          if(ilopt.eq.1)then

          do i=1,npar
            call psort_and_set_rank(i)
          enddo

          else 

          do i=1,npar
            call sort_and_set_rank(i)
          enddo

          endif         
!--------------------------------------------
          if(ilopt.eq.1)then
          do i=1,npar
           call prank_sort(i,nt)
          enddo
          else               
          do i=1,npar 
           call rank_sort(i,nt)
          enddo
          endif
!--------------------------------------------
	write(82,789)'LAST MOLS'
c         ************ LAST MOLS *************
	do m1=1,15
          m2 = 2
	  kk = 1

          if(ilopt.eq.1)then
          call psubpar(kk,m1,m2,iopt,ifopt,fseq,aii,phi,npar)
          else
	  call subpar(kk,m1,m2,iopt,ifopt,fseq,aii,phi,npar)
          endif

          if(ilopt.eq.1)then
          call pmolgen(kk,m1,phi,1,ifopt,eflex)
          else
	  call molgen(kk,m1,phi,1,ifopt,eflex)
          endif 
        
	  
cs------PHI of all structures of best 15------------
        do ik4=1,npar
        opmophi(m1,ik4)=phi(ik4)
        if(icomment.eq.1)print *,'opmophi',m1,opmophi(m1,ik4)
        enddo
cs--------------------------------------------------
        if(ilopt.eq.1)then
          if(ifff.eq.1)then
           en_final(m1,m1)=pampene(m1,m1)
          else
           en_final(m1,m1)=pecpene(m1,m1)
          endif
        else if(ilopt.eq.2)then 
          if(ifff.eq.1)then
            rmmff=rmmffene(m1,m1)
            en_final(m1,m1)=rmmff
          else if(ifff.eq.2)then
            rgaff=rgaffene(m1,m1) 
            en_final(m1,m1)=rgaff
          endif
        endif

        if(ilopt.eq.1)then
        call peplp(plpe,hb,steric,1)
        enplp(m1,m1)=plpe
        else
        call eplp(plpe,hb,steric,1)
        enplp(m1,m1)=plpe
        endif
        
        en_flex(m1,m1)=eflex
	write(82,788)m1,m1,hb,steric,plpe,en_final(m1,m1),eflex

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

        call output(i,j,1,ijx,et,iopt,fseq,aii,npar,ifopt) !included npar in arguments

	write(*,*)'generated MOLS optimal Structure NO : ',ijx,et
c*******Flushing standard output************
	call flush(6)
c*******End of Flushing ********************
	write(31,*)'generated MOLS optimal Structure NO :',ijx,et
c-------mols.log------------------------------------------------------
      do iil=1,15
      if(e_final(iil,iil).eq.et)then
      ir=iil       
      endif
      enddo
        write(77,778)ijx,en_final(ir,ir),enplp(ir,ir),en_flex(ir,ir),
     &  e_final(ir,ir)

        if(icomment.eq.1) print *,'ir-mols-best',ir
        if(icomment.eq.1)write(*,778)ijx,en_final(ir,ir),enplp(ir,ir),
     &   en_flex(ir,ir),e_final(ir,ir)
c--------------------------------------------------------------------
333	enddo

cs------passing optimal parameters of receptor into 'opmo' ---------- 
cs	if(ifopt.eq.2)then
cs	do ik3=nn+7,npar
cs      opmo(ijx,ik3) = opmophi(ir,ik3)
cs      enddo
cs	endif

	return
        end

c********************************************************************
        subroutine write_par(angl,iseed,angle,npar)

	include 'mols.par'
        common/order/nn,mm
        dimension angll(maxpar,maxord),angle(maxpar,maxord),
     &  angl(maxpar,maxord)
       open(unit=4,file='molspargen.inp',status='unknown')

        do i=1,npar
          do j=1,mm
            angll(i,j)=angl(i,j)
          enddo
        enddo
        amm=float(mm)
        do j=1,npar
         ki=0
         do i=1,100000
c         y=ran(iseed)
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


        do i=1,npar
c       write(4,10)(angle(i,j),j=1,mm)
 	write(4,*)(angle(i,j),j=1,mm)
c10     format(37f6.1)
        enddo
    
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
10    format(37f6.1)
	do 500 l = 1, npar
	do 500 k = 1, mm
	do 500 j = 1, mm
	  i = mod( ((j-1) * (l-1) + (k-1)), mm) + 1
	  e(i,j,l)  = p(l,k,1)
c	write(31,*)i,'  ',j,'  ',l,'  ',e(i,j,l)

500	continue
        close (unit=5)
c20	format(8f6.1)
c	open (unit=6,file='pargen.out',status='unknown')
c	 write (6,20) (((e(i,j,l),l=1,npar),j=1,mm),i=1,mm)
c      close (unit=6)        

	return
	end

c*******************************************************************************************
 
        subroutine molgen(lx,ly,phi,tst,ifopt,e_flex)
c	program molgen(phi)
        include 'mols.par'
        common /parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /par/ natom, ntor, nhb, ns, lres
        common /vectors/iv(maxpar,4)
        common /order/nn,mm
        common /getrot/inrot,irotb(maxi,2),el,ilsrot(maxi,maxi),
     &   iatrot(maxi),rx(100,3),ire(maxi,maxi),ind(maxi,maxi),big
        common /left/ le(maxi,maxi),innd(maxi),ielenum(maxi),
     &  bonum(maxi)
	common /cen/bx,by,bz
        common /recep/np,nres
        common /native/nat,pepfile

c       integer le(maxi,maxi),elenum(maxi,maxi)

        dimension x_one(maxatm,3),x_two(maxatm,3)
        dimension phi(maxpar)
        integer tst,nn,ci
        real rotang,theta,psi

        real cx,cy,cz,bx,by,bz
 
        if(nat.eq.1) goto 314 

	nn=ntor
        do k=1,natom
           do ki=1,3
             x_one(k,ki)=rx(k,ki)
           enddo
        enddo

c-------if rotatable bonds in ligand is zero------------------------
        if(ntor.eq.0)then
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
        goto 315

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
        do j=1,ielenum(if) !ielenum = total number of ligand atoms
	k=le(if,j)
           do ki=1,3
           x_two(k,ki)=x_one(k,ki)
           enddo
        enddo
c------------------------------------------------------
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
        if(ifopt.eq.2)then 
        call flexall(lx,ly,phi,e_flex,1) 
        endif


csy      do k=1,natom
csy          print *,ylig(k,1),ylig(k,2),ylig(k,3)
csy      enddo
csy      stop


        goto 315
314     do k=1,natom
           do ki=1,3
             ylig(k,ki)=rx(k,ki)
           enddo
        enddo

        if(ifopt.eq.2)then
        call flexall(lx,ly,phi,e_flex,1)
        endif

315     return
        end
c***********************************************************************
c       function ampene(ie,je)
        function fampene(hx,hy,amberene)        
 
	include 'mols.par'
        common/procrda/xpro(maxatm,8)
        common/proranges/jstartpro(maxatm,10),jendpro(maxatm,10),
     & j1_4pro(maxatm,12)
        common/proenergy/proe_final(maxord,maxord),pro_final
        common/par/ natom, ntor, nhb, ns, lres
        common/prohb/iprohb1(maxhb),iprohb2(maxhb),cpro(maxhb),
     &  dpro(maxhb)
        common/calls/ncalls
        common /comment/icomment
        common /part/ipart,iptor
        common /partners/qatname(maxatm),qx1(maxatm),qx2(maxatm),
     &  qx3(maxatm),qresname(maxatm),qcid(maxatm),qresnum(maxatm)      
        common /propdb/ihflex,patname(maxatm),px1(maxatm),px2(maxatm),
     &  px3(maxatm),presname(maxatm),pcid(maxatm),presnum(maxatm)
	common /flex/iponly
	common /protres/resno(maxatm)
	common /tweight/xw1,xw2,xw3

        integer resno,jj,ji,hx,hy,jp(5000),jr(5000),iponly
        real y(maxatm,8),amberene
        real px1,px2,px3,qx1,qx2,qx3,xw1,xw2,xw3
        character*4 patname,qatname
	character*3 qresname,presname
        character*1 pcid,qcid
        integer qresnum,presnum
        character xx1*3,yy1*3,xx2*1,yy2*1,xx3*4,yy3*4
c----------------------------------------------------------------------        
	if(icomment.eq.2)icount=0
	if(hx.eq.1.and.hy.eq.1)then
	jq=0
        if(icomment.eq.1)print *,'ihflex-fampene',ihflex,ipart
c       print *,'fampene',qatname(1),qresname(1),qcid(1)
c       print *,'fampene',patname(1),presname(1),pcid(1)
        do i=1,ihflex
        do j=1,ipart
	if(icomment.eq.2)icount=icount+1
        xx1=presname(i)
        xx2=pcid(i)
        xx3=patname(i)
        yy1=qresname(j)
        yy2=qcid(j)
        yy3=qatname(j)
c        print *,'fampene',presnum(i),qresnum(j),
c     &  xx1(2:4),yy1(2:4),xx2(1:1),yy2(1:1),xx3(2:4),yy3(2:4)
        if(xx1(1:3).eq.yy1(1:3).and.xx2(1:1).eq.yy2(1:1)
     &  .and.xx3(1:4).eq.yy3(1:4)
     &  .and.presnum(i).eq.qresnum(j))then
cs        print *,'fampene-before',j,presnum(i),qresnum(j),xpro(j,1),
cs     &  xpro(j,2),xpro(j,3)
cs	print *,'qresnum',qresnum(j)
c	if(icomment.eq.2)print *,'i',i,'j',j
	jq=jq+1
	jp(jq)=j
	jr(jq)=i
        xpro(j,1)=px1(i)
        xpro(j,2)=px2(i)
        xpro(j,3)=px3(i)
cs        print *,'fampene-after',j,presnum(i),qresnum(j),xpro(j,1),
cs     &  xpro(j,2),xpro(j,3)
        endif
        enddo
        enddo
	iponly=jq
	if(icomment.eq.2)print *,'iponly',iponly
	ELSE
	do i=1,iponly
c	if(icomment.eq.2)print *,i,'jp',jp(i),'jr',jr(i)
	xpro(jp(i),1)=px1(jr(i))
        xpro(jp(i),2)=px2(jr(i))
	xpro(jp(i),3)=px3(jr(i))
	if(icomment.eq.2)icount=icount+1
	enddo
	ENDIF
	if(icomment.eq.2)print *,'fampene-coord',icount
c----------------------------------------------------------------------
c       open(unit=1,file='resnum.txt',status='old')
cs      do i=1,ipart
c	icount=icount+1
c       read(1,*)resno(i)
cs	print *,'resno',resno(i)
cs      enddo
c	pause
c       close(1)
	if(icomment.eq.2)print *,'after-resnum',icount

        ees=0.0
        enb=0.0
        enb1_4=0.0
        ees1_4=0.0
        ehb=0.0
        etot=0.0
        amberene=0.0

cs        do i=1,ipart
cs         do j=1,8
cs	  if(icomment.eq.2)icount=icount+1
cs          y(i,j)=xpro(i,j)
cs         enddo
cs        enddo
	if(icomment.eq.2)print *,'after-xpro',icount

	do i=1,ipart
          n_range=ifix(xpro(i,7))
         do ij=1,n_range
          k1=jstartpro(i,ij)
         k2=jendpro(i,ij)
          do j=k1,k2
          dis = dist(xpro(i,1),xpro(i,2),xpro(i,3),
     &    xpro(j,1),xpro(j,2),xpro(j,3))
            if(dis.lt.0.01) then
!              e_final(ie,je)=1.0e+25
               ampene = 1.0e+25
              return
            endif

cs         ji=ij-i
cs         if(ji.ge.3)then
c-------------------------------------------------------------
cs         if(dis.gt.3.0)then !added by sam
c..The above condition is given to remove atoms of the same
c..residue from van der Waals energy calculation failing which
c..caused high non-bonded energy..............................
c-------------------------------------------------------------
        if(resno(i).ne.resno(j))then!added by sam
	if(qatname(i).eq.' SG '.and.qatname(j).eq.' SG ')goto 212!added by sam
cs        'SG'added for relaxing disulphide bonds
cs        if(dis.gt.2.2)then         !added by sam  
	  if(icomment.eq.2)icount=icount+1
          call force(i,j,dis,enbt,eest)
           ees=ees+eest
           enb=enb+enbt
cs        print *,'force-enbt',i,j,qatname(i),qatname(j),dis,enb
cs        endif              !added by sam 
        endif               !added by sam
212    enddo
       enddo

c	print *,'after-force-nb',icount


         n1_4=ifix(xpro(i,8))
         do ij=1,n1_4
            j=j1_4pro(i,ij)
          dis = dist(xpro(i,1),xpro(i,2),xpro(i,3),xpro(j,1),
     &    xpro(j,2),xpro(j,3))
	  
          if(dis.lt.0.5) then
              STOP 'Input coordinates are wrong !!'
            endif
	    if(icomment.eq.2)icount=icount+1
c	    if(resno(i).eq.resno(j))then !calculating 1-4 interactions for each residues
            call force(i,j,dis,enbt,eest)
            ees1_4=ees1_4+0.5*eest
            enb1_4=enb1_4+0.5*enbt
c	    endif
c	print *,i,j,dis,'ees1_4',ees1_4,'enb1_4',enb1_4
        enddo    
c	print *,'after-force_1-4',icount       
         enddo
	if(icomment.eq.2)print *,'after-force-nb_1-4',icount        
        do i=1,nhb
	if(icomment.eq.2)icount=icount+1
         dis = dist(xpro(iprohb1(i),1),xpro(iprohb1(i),2),
     &   xpro(iprohb1(i),3),xpro(iprohb2(i),1),xpro(iprohb2(i),2),
     &   xpro(iprohb2(i),3))

         ess=dis*dis
         artwo=1.0/ess
         arten=artwo**5
         artwelve=arten*artwo
         ehb=ehb+cpro(i)*artwelve-dpro(i)*arten
         arstar=xpro(iprohb1(i),5)+xpro(iprohb2(i),5)
         epsilon=sqrt(xpro(iprohb1(i),6)*xpro(iprohb2(i),6))
         aaa=epsilon*(arstar**12)
         ccc=2.0*epsilon*(arstar**6)
         arsix=artwo*artwo*artwo
         ef1=aaa*artwelve
         ef2=ccc*arsix
         enbt=ef1-ef2
         enb=enb-enbt
        enddo
	if(icomment.eq.2)print *,'after-nhb',icount
c       xw3=1.0        
	amberene=(xw3*((enb)+(ees)+(ehb)+(enb1_4)+(ees1_4)))
c	amberene=(xw3*((enb)+(ees)+(ehb)))
        if(icomment.eq.1)then
        print *,'PROTEIN-AMBER'
	print *,'xw3',xw3
        print *,'enb',enb
        print *,'ees',ees
        print *,'ehb',ehb
        print *,'enb1_4',enb1_4
        print *,'ees1_4',ees1_4
        print *,'amberene',hx,hy,amberene
        endif
c------------------------------------------------------
        etot=amberene
        pro_final=etot
        return
        end
c*******************************************************************
        subroutine force(if,jf,diss,enbt,eest)

	include 'mols.par'
        common/procrda/xpro(maxatm,8)
        common /part/ipart
cs      real y(maxatm,8)
        

cs        do i=1,ipart
cs        do j=1,8
cs        y(i,j)=xpro(i,j)
cs        enddo
cs        enddo

        arstar=xpro(if,5)+xpro(jf,5)
        epsilon=sqrt(xpro(if,6)*xpro(jf,6))
        aaa=epsilon*(arstar**12)
        ccc=2.0*epsilon*(arstar**6)
        ess=diss*diss
        artwo=1.0/ess
        arsix=artwo*artwo*artwo
        artwelve=arsix*arsix
        eest=332.0*xpro(if,4)*xpro(jf,4)*artwo/4.0
        ef1=aaa*artwelve
        ef2=ccc*arsix
        enbt=ef1-ef2        
        return                  
        end

c***********************************************************************
c       function ecpene(ie,je)
        function fecpene(rt)                
 
	include 'mols.par'
	parameter(mxtyat = 18)
        common/pcrda/x(maxatm,8) 
        common/ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,12)
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common/energy/e_final(maxord,maxord) 
        common /par/ natom, ntor, nhb, ns, lres
        common/hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
        common /ecpp/ aij(mxtyat,mxtyat),cij(mxtyat,mxtyat),
     &  a14(mxtyat,mxtyat),ihb(mxtyat,mxtyat),ahb(mxtyat,mxtyat),
     &  chb(mxtyat,mxtyat)
        common/calls/ncalls
        common /part/ipart

        integer fatom,rt
        real y(maxatm,8)

        fatom=ipart
        
	cdc=(22.0/(7.0*180))
        ees=0.0
        enb=0.0
        ehb=0.0
	etor=0.0
        enb14=0.0
        ees14=0.0
	ehb14=0.0
        etot=0.0


        do i=1,fatom
          do j=1,8
            y(i,j)=x(i,j)
          enddo
        enddo


 	do in = 1,rt
	etor = etor + u0(in)*(1.0+sn(in)*cos(cdc*((phi(in)-180.0))*tn(in)))
	enddo

	do i=1,fatom
          n_range=ifix(y(i,7))
          write(*,*)'n_range',n_range,y(i,7)
         do ij=1,n_range
          k1=jstart(i,ij)
         k2=jend(i,ij)
          do j=k1,k2
          dis = dist(y(i,1),y(i,2),y(i,3),y(j,1),y(j,2),y(j,3))
        write(*,*)'dis',dis
            if(dis.lt.0.01) then
	       ecpene = 1.0e+25
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
        ecpene=enb+ees+ehb+etor
        print *,'ees',ees,'enb',enb,'ehb',ehb

        return
        end
c***********************************************************************
c********************************************************************
        subroutine average(npar)

	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/energy/e_final(maxord,maxord) 
        common/out/e_out(maxatm),p_out(maxpar,maxatm)
        common/mean/avrg1(maxpar,maxord)
        common/order/nn,mm

        data en_k_t/1.98578e+03/
	open (unit=8,file='average.out',status='unknown')

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
	write(8,*)l,k,avrg1(l,k)
	close(unit=8) 
500	continue
400    continue

        return
        end
c*********************************************************************
      
        subroutine sort_and_set_rank(is1)

	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/mean/avrg1(maxpar,maxord)
        common/order/nn,mm
        common/freq/ifreq(maxpar,maxord,maxord)

        dimension rank(maxord),icc(maxord)
c	open (unit=9,file='srank1.out',status='unknown')
c	open (unit=10,file='srank2.out',status='unknown')
c	open (unit=11,file='srank3.out',status='unknown')
c	open (unit=12,file='srank4.out',status='unknown')

c	do j=1,mm
c	write (9,*),avrg1(is1,j),p(is1,j,1)
c       enddo
c	close(unit=9)
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
c	do j=1,mm
c	write (10,*),avrg1(is1,j),p(is1,j,1)
c        enddo
c	close(unit=10)

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

c	do j=1,mm
c	write (11,*)rank(j),p(is1,j,1)
c       enddo
c	close(unit=11)

        do i=1,mm
          jxx=ifix(rank(i))
          ifreq(is1,i,jxx)=ifreq(is1,i,jxx)+1
          p(is1,i,2)=p(is1,i,2)+rank(i)
          p(is1,i,3)=p(is1,i,3)+rank(i)*rank(i)
        enddo


c	do j=1,mm
c	write (12,*) p(is1,j,1),p(is1,j,2),p(is1,j,3)
c        enddo
c	close(unit=12)

        return
        end 
c********************************************************************
       
        subroutine rank_sort(ir1,nt)

	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/freq/ifreq(maxpar,maxord,maxord)
        common/order/nn,mm

        dimension ee(maxord)
c	open (unit=13,file='rank1.out',status='unknown')
c	open (unit=14,file='rank2.out',status='unknown')
c	open (unit=15,file='rank3.out',status='unknown')
c	open (unit=16,file='rank4.out',status='unknown')

c       calculate the average and sd
c	do j=1,mm
c	write (13,*) p(ir1,j,1),p(ir1,j,2),p(ir1,j,3)
c	enddo
c	close(unit =13)

        do i=1,mm
          sum_x=p(ir1,i,2)
          sum_x2=p(ir1,i,3)
          ave=sum_x/nt
          sd=sqrt((sum_x2/nt)-((sum_x/nt)*(sum_x/nt)))
          p(ir1,i,2)=ave
          p(ir1,i,3)=sd
        enddo       
c	do j=1,mm
c	write (14,*) p(ir1,j,1),p(ir1,j,2),p(ir1,j,3)
c	enddo
c	close(unit =14)

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

c	do j=1,mm
c	write (15,*) p(ir1,j,1),p(ir1,j,2),p(ir1,j,3)
c	enddo
c	close(unit =15)

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

c	do j=1,mm
c	write (16,*) p(ir1,j,1),p(ir1,j,2),p(ir1,j,3)
c	enddo
c	close(unit =16)

        return
        end

c*********************************************************************
       subroutine output(cx,cy,nou,ijx,et,iopt,fseq,aii,npar,ifopt)

	include 'mols.par'
        common /out/e_out(maxatm),p_out(maxpar,maxatm)
        common /order/nn,mm
        common /par/ natom, ntor, nhb, ns, lres,is,ie
	common /scc/ sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
        common /sscc/ ssc(maxatm,maxpar,5),nsrot(maxpar),
     &  nschi(maxpar)
	common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
        common /fnam/if1,pf1,pf2,pf3,pf4,pf5,of1,of2,of3,of4,
     &  of5,of6

css	common /fnam/ if1,pf1,pf2,pf3,pf4,pf5,of1,of2,of3,of4
cs      common /fnam/ of0,pf2,pf3,pf4,pf5,of6,of1,of2,if1,of5       
        common /dcom/ipatom
        common /recep/np,nres
        common /scang/frange,rang
        common /comment/icomment
        common /rctl/iscopt
        common /optmo/opmophi(maxpar,maxpar)
        common /energy/e_final(maxord,maxord)
        

        character*128 if1,pf1,pf2,pf3,pf4,pf5,pf6,of1,of2,of3,of4,
     &  of0,of5,of6

css	character*128 if1,pf1,pf2,pf3,pf4,pf5,of1,of2,of3,of4,of0
cs      character*128 of0,pf2,pf3,pf4,pf5,of6,of1,of2,if1,of5

	character fseq(maxatm)*1,prt*80,prt1*80 
        integer exresid,cx,cy,npar
        real ex1,ex2,ex3

	dimension pp(maxpar)

c       print *,'npar=',npar,'nn=',nn
c       print *,'aii=',aii
c       print *,'from output routine'
c       print *,(p_out(j,1),j=1,npar)

170     format(a23,i3,4x,3f8.3,a24)
        npar = ntor+6+np
        i = nou

cs      if(icomment.eq.1)then
cs      do ik=1,npar        
cs      print *,'p_out',i,p_out(ik,i)
cs      enddo
cs      endif
        if(icomment.eq.1)write(*,*),'output-ifopt',ifopt
        if(icomment.eq.1)write(*,*),'rang-output',rang        
c        do ii=1,nres
c        write(*,*)rrid(ii)        
c        enddo
c        print *,'rid',rid
cs	i = nou
cs	if(iopt.eq.3) then  
cs	j = 1
cs	k = 1
cs	jh = 1
cs172	if(fseq(k).eq.'P') go to 193
cs	pp(j) = p_out(jh,i)
cs	j = j + 1
cs	jh = jh + 1
cs193	pp(j) = p_out(jh,i)
cs	j = j + 1
cs	jh = jh + 1
cs	if(nrot(k).ne.0.) then
cs	j1 = p_out(jh,i)
c	j1 = int((j1/10.0) + 1.0)
cs	j1 = int((j1/aii) + 1.0)
c	print *,j1,aii
cs	j1 = mod(j1,nrot(k))
cs	if(j1.eq.0) j1 = nrot(k)
cs	do j2 = 1, nchi(k)
cs		pp(j) = sc(k,j1,j2)
cs		j = j + 1
cs	enddo
cs	jh = jh + 1
cs	endif
cs	k = k + 1
cs	if(k.le.lres) go to 172	
cs	endif  
c------- write MOLS protein PDB-----------------------------------
        
c        natp=ipatom
c        open(unit=11,file='pfile1.pdb',status='old')
c        open(unit=12,file='pfile_mols.pdb',status='unknown')
c        do i=1,natp
c        read(11,170)prt,exresid,ex1,ex2,ex3,prt1
c        write(12,170)prt,exresid,ex1,ex2,ex3,prt1
c        enddo
c        close(12)
c        close(11)

c----------------------------------------------------------------
	open(unit=24,file=of1,status='unknown')
c          write(24,*)(p_out(j,i),ji=1,npar),e_out(i)
c           write(24,*)(p_out(j,i),j=1,npar),e_out(i)

css        IF(ifopt.eq.1)then
c----------------------------------------------------------------
	if(iopt.ne.3) then
	do k1=1,nn
	opmo(ijx,k1) = p_out(k1,i)
        if(icomment.eq.1)print *,'opmo',ijx,k1,opmo(ijx,k1)
	enddo 
	emo(ijx) = e_out(i)
	do k1 = 1,nn
c	pp(k1) = pp(k1) - 180.0
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
c       tt2 = (tt1/10.0)*.28
        tt2 = (tt1/10.0)*0.13888889
c       tt2 = (tt1/10.0)*0.55555556
        endif
c       tt2 = tt2-5.0
        tt2 = tt2-2.5
c       tt2 = tt2-10.0
        p_out(k2,i)=tt2
        enddo
        do k2=nn+4,nn+6
        opmo(ijx,k2) = p_out(k2,i)
        if(icomment.eq.1)print *,'opmo-nn+4,nn+6',ijx,k2,opmo(ijx,k2)
        enddo
        p_out(nn+6,i)=p_out(nn+6,i)/2
c-------opmo for receptor--sam added------------------------
	if(ifopt.eq.2)then
css      do k=nn+7,npar
css      opmo(ijx,k)=p_out(k,i)
css      enddo
css      endif
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
c----------------------------------------------------------
c************
C       write(24,10) ijx,e_out(i),(pp(j),j=1,nn)
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
        
C**********************
c       STORING THE TRANS AND ROT PARAM
        ix = 5
        iy = 1

281     tt1 = p_out(npar-ix,i)
c
        opmo(ijx,nn+iy) = tt1
        if(icomment.eq.1)print *,'opmo-3',ijx,nn+iy,opmo(ijx,nn+iy)
c
        if(tt1.eq.0.0) then
        tt2 = 0.0
        else
c       tt2 = (tt1/10.0)*.28
        tt2 = (tt1/10.0)*0.13888889
c       tt2 = (tt1/10.0)*0.55555556
        endif
c       tt2 = tt2-5.0
        tt2 = tt2-2.5
c       tt2 = tt2-10.0
        pp(nn+iy) = tt2
c       opmo(ijx,nn+iy) = tt2
c       print *,ix,iy,npar-ix,nn+iy
        ix = ix-1
        iy = iy+1
        if(iy.le.3) go to 281
        ix = 2
        iy = 4
282     pp(nn+iy) = p_out(npar-ix,i)
        opmo(ijx,nn+iy) = p_out(npar-ix,i)
        if(icomment.eq.1)print *,'opmo-4',ijx,nn+iy,opmo(ijx,nn+iy)
c       print *,ix,iy,npar-ix,nn+iy
        ix = ix-1
        iy = iy+1
        if(iy.le.6) go to 282


        pp(nn+6) = pp(nn+6)/2
c       opmo(ijx,nn+6) = opmo(ijx,nn+6)/2
C*******************************************************************
cs	if(ifopt.eq.2)then
cs	do k=nn+7,nn+7+np
cs	opmo(ijx,k)=p_out(k,i)
cs	enddo
cs	endif
C*******************************************************************
        write(24,*) ijx,e_out(i),(pp(j),j=1,npar)
        write(*,*)'two',ijx,e_out(i),(pp(j),j=1,npar)
        endif 

	et = e_out(i)
        if(ijx.eq.ie) close(unit=24)
22      format(a6)
23      format(a3)
c	stop
c---------------------------------------------------------------------
css        ELSE
css        et=e_out(i)
css        emo(ijx)=e_out(i)
css        do iil=1,15
css        if(e_final(iil,iil).eq.et)then
css        ir=iil
css        endif
css        enddo
    
css        do ik4=1,npar
css        opmo(ijx,ik4)=opmophi(ir,ik4)
css        if(icomment.eq.1)print *,'opmo-opmophi',ir,opmo(ijx,ik4)        
css        enddo
css        ENDIF
c-----------------------------------------------------------------------
css     if(ijx.eq.ie) close(unit=24)
cs----------------------------------------------------------------------
        return
        end

c********************************************************************
c*********************************************************************
c***********************************************************************
	subroutine best(ib,npar)
	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common/energy/e_final(maxord,maxord) 
        common/order/nn,mm
        common/out/e_out(maxatm),p_out(maxpar,maxatm)
        common /comment/icomment

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
c        open(unit=8,file='mseed.out',status='unknown')
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
c             write(8,*)r
	jseed(nss) = r	
1       b = r
	enddo

2       return
        end
c*************************************************************************
	subroutine anggen(ang,fseq,npar,iopt,aii)
	include 'mols.par'
        common/order/nn,mm
        common /par/ natom, ntor, nhb, ns, lres
	dimension ang(maxpar,maxord)
	character fseq(maxatm)*1 
	
c	open (unit=21,file='angles.inp',status='unknown')
c	type *,'enter nn (total no. of parameters)'
c	read *,nn
c	type *,'enter mm ((37)maximum no. of values for each parameters)'
c	read *,mm

  10    format(37f6.1)
cc        type *,'enter the interval to generate input angle values'
cc        read *,ii
cccc	ii = 10
cc        type *,'enter the starting value to generate input angle values'
cc	read *,ll

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
c	print *,mm,npar
	if(mm.lt.npar) then
	do i1 = 1,100
	k1 = 2
181	if(i1/k1*k1.eq.i1) go to 191
	k1 = k1 + 1
	if(i1.le.i1/2) go to 181
	if(i1.gt.npar) go to 201
c	print *,i1
191	enddo
201	mm = i1
	i1 = 0
	endif
c	print *,mm,npar
cc-------------------------------------------------------------------
	aii = 360.0/(mm-1)
	all = 0
	do i=1,npar
	ang(i,1)=all
        b=all
	do j=2,mm
c        write(21,10)ang(i,j)
	ang(i,j)=b+aii
        b=ang(i,j)
        enddo
	enddo
c        write(21,10)((ang(i,j),j=1,mm),i=1,npar)
c       write(6,10)((ang(i,j),j=1,mm),i=1,npar)
c	close(unit=21)
	return 
	end
c*************************************************************************
	subroutine scinp(fseq)
	include 'mols.par'
        common /par/ natom, ntor, nhb, ns, lres
	character fseq(maxatm)*1,ch*1, tcode*4, scode*1
	common /scc/ sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
cc	dimension sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
20	format(a4,a1,1x,i2,i2)
21	format(20x,5f7.1)
11	format(a18,i3)

c	write(6,*) 'The sequence   : ',(fseq(i),i=1,lres)

 	open(unit=33,file='SC.lib',status='old')
	do i = 1,lres
	rewind 33
	do i1 = 1,234
	read(33,20) tcode, scode, mrot, mchi
c	write(31,*)tcode, scode, mrot, mchi
	if(scode.eq.fseq(i).and.mrot.ne.0) then
	  do i2 = 1, mrot
	  read(33,21) (sc(i,i2,i3),i3=1,mchi)
c	write(31,*)(sc(i,i2,i3),i3=1,mchi) !this subroutine is not executed
	  enddo
	nrot(i) = mrot
	nchi(i) = mchi
	go to 212
	endif
	enddo
212	enddo
	do i = 1,lres
	do j = 1, nrot(i)
c	write(31,*)(sc(i,j,k),k=1,nchi(i))
	enddo
	enddo
	close(unit=33)
	return
	end
	
c************************************************************************	
	subroutine subpar(l1,l2,l3,iopt,ifopt,fseq,aii,phi,npar)
	include 'mols.par'
        common/parameters/e(maxord,maxord,maxpar),p(maxpar,maxord,3)
        common /par/ natom, ntor, nhb, ns, lres
        common/vectors/iv(maxpar,4)
        common/order/nn,mm
	common /scc/ sc(maxatm,maxpar,5),nrot(maxpar),nchi(maxpar)
        common /sscc/ ssc(maxatm,maxpar,5),nsrot(maxpar),
     &  nschi(maxpar),mn(maxpar)
        common /ctl/ipf,iff,icint
        common/pctl/ifff,ioptt
        common/rctl/iscopt
        common /recep/np,nres
        common/scang/frange,rang
        common/comment/icomment
	character fseq(maxatm)*1 
c	dimension sc(1500,50,5),nrot(50),nchi(5)
        dimension x_one(maxatm,3),x_two(maxatm,3)
        dimension phi(maxpar)

        if(icomment.eq.1) print *,'iscopt',iscopt
        if(icomment.eq.1) print *,'rang-subpar',rang

cs      IF(ifopt.eq.1)THEN
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
c	j1 = int((j1/10.0) + 1.0)
	j1 = int((j1/aii) + 1.0)
	j1 = mod(j1,nrot(k))
	if(j1.eq.0) j1 = nrot(k)
c	write(31,*)j1,k,nrot(k),nchi(k)
c	write(*,*)j1,k,nrot(k),nchi(k)
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
c	j1 = int((j1/10.0) + 1.0)
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
c       do i = nn+1, nn+3
c       if(l3.eq.1) tt1 = e(l1,l2,i)
c       if(l3.eq.2) tt1 = p(i,l2,1)
c       phi(i)=tt1
c       enddo
c       do i = nn+4, nn+6
c       if(l3.eq.1) phi(i) = e(l1,l2,i)
c       if(l3.eq.2) phi(i) = p(i,l2,1)
c       enddo
c       write(*,*)(phi(ii),ii=1,nn+6)

        do i = nn+1, nn+3
        if(l3.eq.1) tt1 = e(l1,l2,i)
        if(l3.eq.2) tt1 = p(i,l2,1)
        if(tt1.eq.0.0) then
        tt2 = 0.0
        else
c       tt2 = (tt1/10.0)*.28
        tt2 = (tt1/10.0)*0.13888889
c       tt2 = (tt1/10.0)*0.55555556
        endif
c       tt2 = tt2-5.0
        tt2 = tt2-2.5
c       tt2 = tt2-10.0
        phi(i) = tt2
        enddo
        do i = nn+4, nn+6
        if(l3.eq.1) phi(i) = e(l1,l2,i)
        if(l3.eq.2) phi(i) = p(i,l2,1)
        enddo

        if(iopt.eq.3) then
c       print *,nn,npar
c       change ix according to the no. of parameters added for trans and rot
c       ix = 5
c       iy = 1
        ix = 5
        iy = 1

c181    tt1 = e(l1,l2,npar-ix)
181     if(l3.eq.1) tt1 = e(l1,l2,npar-ix)
        if(l3.eq.2) tt1 = p(npar-ix,l2,1)
        if(tt1.eq.0.0) then
        tt2 = 0.0
        else
c       tt2 = (tt1/10.0)*.28
        tt2 = (tt1/10.0)*0.13888889
c       tt2 = (tt1/10.0)*0.55555556
        endif
c       tt2 = tt2-5.0
        tt2 = tt2-2.5
c       tt2 = tt2-10.0
        phi(nn+iy) = tt2
c       print *,ix,iy,npar-ix,nn+iy
        ix = ix-1
        iy = iy+1
        if(iy.le.3) go to 181
        ix = 2
        iy = 4
c182    phi(nn+iy) = e(l1,l2,npar-ix)
182     if(l3.eq.1) phi(nn+iy) = e(l1,l2,npar-ix)
        if(l3.eq.2) phi(nn+iy) = p(npar-ix,l2,1)
c       print *,ix,iy,npar-ix,nn+iy
        ix = ix-1
        iy = iy+1
c       if(iy.le.6) go to 182
        if(iy.le.6) go to 182
        endif 
cs      ENDIF
c       print *,'npar =',npar,'nn =',nn
c       write(*,*)(phi(ii),ii=1,nn+6)
c	stop

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

        ELSE

        if(ifopt.eq.2)then

        do i = nn+7,npar
        if(l3.eq.1) tt1 = e(l1,l2,i)
        if(l3.eq.2) tt1 = p(i,l2,1)

         if(tt1.eq.0.0)then
          tt2=0.0
         else
 
c        tt2=(tt1/10.0)*rang
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
c********************************************************

 	return
	
	end
c***********************************************************************

c**************Calculate ligand centroid*****************************
        subroutine calcent(natl,cx,cy,cz)
        include 'mols.par'
        common/ligcrdb/ylig(maxatm,8)
        integer natl
        real cx,cy,cz,cx1,cy1,cz1
c       print *,cx,' ',cy,' ',cz
c       print *,natl
c       print *,'START OF LIGATOM'
        cx=0.0
        cy=0.0
        cz=0.0

c       do i = 1,natl 
c          write (6,FMT='(f8.3,2x,f8.3,2x,f8.3)')y(i,1),y(i,2),y(i,3)
c       enddo

        do i = 1,natl

                cx = cx + ylig(i,1)
                cy = cy + ylig(i,2)
                cz = cz + ylig(i,3)

        enddo 
        cx1 = cx / natl
        cy1 = cy / natl
        cz1 = cz / natl
        cx = cx1
        cy = cy1
        cz = cz1
        return
        end
c**************End of ligand centroid******************************
c*****Translate ligand inside the grid as guided by MOLS******

        subroutine translate(bx,by,bz,tx,ty,tz,natl)
        include 'mols.par'
        common /ligcrdb/ylig(maxatm,8)
        integer natl
        real cx,cy,cz,dx,dy,dz
c       print *,'translation parameters'
c       print *,bx,'  ',by,'  ',bz,'  ',tx,'  ',ty,'  ',tz
c       print *,'centroid-translate',bx,by,bz
c       call calcent(natl,cx,cy,cz)
c       print *,'LIGAND CENTROID'
c       print *,cx,' ',cy,' ',cz

c       print *,'BEFORE TRANSLATION'
c       do i = 1,natl
c       write(6,FMT='(f8.3,2x,f8.3,2x,f8.3)')y(i,1),y(i,2),y(i,3)
c       enddo

c       print *,natl

c*******USED FOR MOVING CENTROID AGAIN TO ORGIN******
        call calcent(natl,cx,cy,cz)
        do i = 1,natl

                ylig(i,1) = cx-ylig(i,1)
                ylig(i,2) = cy-ylig(i,2)
                ylig(i,3) = cz-ylig(i,3)
        enddo
c       print *,'centroid after rotation: ',cx,cy,cz
c*******USED FOR MOVING CENTROID AGAIN TO ORGIN******


        do i = 1,natl
                
                ylig(i,1) = ylig(i,1)+(bx)+tx
                ylig(i,2) = ylig(i,2)+(by)+ty
                ylig(i,3) = ylig(i,3)+(bz)+tz
       enddo

       

c       print *,'AFTER TRANSLATION'
c       do i = 1,natl
c       write(6,FMT='(f8.3,2x,f8.3,2x,f8.3)')y(i,1),y(i,2),y(i,3)
c       enddo
c       stop
        return
        end
c*************End of Translation****************************

c*****Rotate ligand by angles Rx, Ry and Rz about****
c*****       x, y and z axiz respectively        ****

c       subroutine rotate(rx,ry,rz,natl)
        subroutine rotate(bx,by,bz,rotrang,theta,psi,natl)
        include 'mols.par'
        common/ligcrdb/ylig(maxatm,8)
        real cx,cy,cz,theta1,psi1,rotrang1
        integer natl,r
        real a,b,c,d,e,f,ad,bd,rmat(3,3)
        real w1,w2,qx,qy,qz,tw,tx,ty,tz,tw1,tx1
        real ax,ay,az,sx,sy,sz,s,t
10      format(a6,2x,i5,a10,i5,4x,3f8.3)
	
c	print*,'bx=',bx
c	print*,'by=',by
c	print*,'bz=',bz
        rotrang1 = rotrang
        theta1 = theta
        psi1 = psi
c	print*,psi1
c	print*,rotrang1
c	print*,theta1
c       print *, 'rotation angles:',rx,ry,rz
        theta1=theta1*3.1415927/180.0
c       psi angle range {0<=psi<=180}
        psi1=(psi1*3.1415927/180.0)/2
c       rz=rz*3.1415927/180.0
c       print *,'rotang:',rotrang
c       CONVERT DEGREES TO RADIAN
        rotrang1=rotrang1*3.1415927/180.0
c       print *,'rotrang(radian):',rotrang

        r = 1
        w2 = 0.0
        sx = 0.0
        sy = 0.0 
        sz = 0.0
c       CREATE POINTS ON THE UNIT SPHERE USING SPHERICAL COORD

        sx = r*cos(theta1)*sin(psi1)
        sy = r*sin(theta1)*sin(psi1)
        sz = r*cos(psi1) 

c       MOVE THE SPHERE CENTRE TO THE LIGAND CENTROID
c       print *,'sphere points',sx,sy,sz
c       sx = sx + bx
c       sy = sy + by
c       sz = sz + bz


c       vdis = dist(sx,sy,sz,bx,by,bz)
        vdis = dist(sx,sy,sz,0.0,0.0,0.0)
c        print *,'radius:',vdis
c       print *,'sphere points',sx,sy,sz
c       write(16,10) 'HETATM',r,"  O   HOH ",r,sx,sy,sz
c       print *,'ligand centroid',bx,by,bz

c       CREATE THE AXIS OF ROTATION

c       ax = (sx - bx)/vdis
c       ay = (sy - by)/vdis
c       az = (sz - bz)/vdis

        ax = sx/vdis
        ay = sy/vdis
        az = sz/vdis

c       print *,'axis of rot:',ax,ay,az

c***************QUATERNION ROTATION******************

c       CREATE UNIT QUATERNION FOR ROTATION

c       w1 = cos(rotrang/2)
c       qx = ax*sin(rotrang/2)
c       qy = ay*sin(rotrang/2)
c       qz = az*sin(rotrang/2)
c
c       print *,'quat mag:',sqrt(w1*w1+qx*qx+qy*qy+qz*qz)
cc      stop

c       do i = 1,natl
c
cc      QUATERNION MULTIPLICATION - TAKING ONLY VECTOR PART
cc      unit quat X coordinate points

c       tw = w1*w2-(qx*y(i,1)+qy*y(i,2)+qz*y(i,3))
c       tx = qx*w1+w2*y(i,1)+qy*y(i,3)-qz*y(i,2)
c       ty = w1*y(i,2)-qx*y(i,3)+qy*w2+qz*y(i,1)
c       tz = w1*y(i,3)+qx*y(i,2)-qy*y(i,1)+qz*w2
c
cc      create q inverse
c       qx = -qx
c       qy = -qy
c       qz = -qz
c
cc      unit quat X coordinate points X inv(unit_quat)
c       tw1 = tw*w1-(tx*qx+ty*qy+tz*qz)
c       tx1 = tx*w1+tw*qx+ty*qz-tz*qy
c       ty1 = tw*qy-tx*qz+ty*w1+tz*qx
c       tz1 = tw*qz+tx*qy-ty*qx+tz*w1
c
c       y(i,1) = tx1
c       y(i,2) = ty1
c       y(i,3) = tz1
cc      print *,tx1,ty1,tz1
c
cc      enddo
cc      stop
c*********END OF QUATERNION ROTATION**************

c*********ANGLE/AXIS ROTATION******************
c       print *,'Before rotation:'
c       do i=1,natl
c       print *,y(i,1),y(i,2),y(i,3)
c       enddo
c       print *,'After rotation:'
        c = cos(rotrang1)
        s = sin(rotrang1)
        t = 1-cos(rotrang1)
c       ROTATION MATRIX
        rmat(1,1) = t*ax*ax+c
        rmat(1,2) = t*ax*ay+s*az
        rmat(1,3) = t*ax*az-s*ay
        rmat(2,1) = t*ax*ay-s*az
        rmat(2,2) = t*ay*ay+c
        rmat(2,3) = t*ay*az+s*ax
        rmat(3,1) = t*ax*az+s*ay
        rmat(3,2) = t*ay*az-s*ax
        rmat(3,3) = t*az*az+c
c       print *,'rmat:',rmat(1,1),rmat(1,2),rmat(1,3)
c       print *,'rmat:',rmat(2,1),rmat(2,2),rmat(2,3)
c       print *,'rmat:',rmat(3,1),rmat(3,2),rmat(3,3)

        call calcent(natl,cx,cy,cz)
        do i = 1,natl

                ylig(i,1) = cx-ylig(i,1)
                ylig(i,2) = cy-ylig(i,2)
                ylig(i,3) = cz-ylig(i,3)
        enddo
c       print *,'centroid before moving: ',cx,cy,cz
c       call calcent(natl,cx,cy,cz)
c       print *,'centroid after moving: ',cx,cy,cz
        do k=1,natl


        tx=rmat(1,1)*ylig(k,1)+rmat(1,2)*ylig(k,2)+rmat(1,3)*ylig(k,3)
        ty=rmat(2,1)*ylig(k,1)+rmat(2,2)*ylig(k,2)+rmat(2,3)*ylig(k,3)
        tz=rmat(3,1)*ylig(k,1)+rmat(3,2)*ylig(k,2)+rmat(3,3)*ylig(k,3)


        ylig(k,1) = tx
        ylig(k,2) = ty
        ylig(k,3) = tz
c       print *,tx,ty,tz
        enddo


cs       do i=1,natl
cs       print *,ylig(i,1),ylig(i,2),ylig(i,3)
cs       enddo
cs       stop

c***********END OF ANGLE/AXIS REPRESENTATION**********

c************EULER ANGLE ROTATION*********************
c       print *, 'rotation angles:',rx,ry,rz
c       a = cos(rx)
c       b = sin(rx) 
c       c = cos(ry) 
c       d = sin(ry) 
c       e = cos(rz) 
c       f = sin(rz) 

ccc     Rotation matrix in zxz convention

c       rmat(1,1) = (e*a)-(c*b*f)
c       rmat(1,2) = (e*b)+(c*a*f)
c       rmat(1,3) = f*d
c       rmat(2,1) = -(f*a+c*b*e)
c       rmat(2,2) = -f*b+c*a*e
c       rmat(2,3) = e*d
c       rmat(3,1) = d*b
c       rmat(3,2) = -d*a 
c       rmat(3,3) = c
ccc     Transpose
cc      rmat(1,1) = (e*a)-(c*b*f)
cc      rmat(1,2) = -(f*a+c*b*e)
cc      rmat(1,3) = d*b
cc      rmat(2,1) = (e*b)+(c*a*f)
cc      rmat(2,2) = -f*b+c*a*e
cc      rmat(2,3) = -d*a
cc      rmat(3,1) = f*d
cc      rmat(3,2) = e*d
cc      rmat(3,3) = c


c       call calcent(natl,cx,cy,cz)
cc      write(6,FMT='(f8.3,2x,f8.3,2x,f8.3)')cx,cy,cz

cc      MOVE THE CENTRE OF MASS TO THE ORIGIN

cc      print *,'BEFORE'
cc      do i = 1,natl
cc      write(6,FMT='(f8.3,2x,f8.3,2x,f8.3)')y(i,1),y(i,2),y(i,3)
cc      enddo

cc      call conformation(0)

c       do i = 1,natl
c               y(i,1) = cx-y(i,1)
c               y(i,2) = cy-y(i,2)
c               y(i,3) = cz-y(i,3)
c       enddo

cc      APPLY ROTATION TRANSFORMATION ON COORDINATES

c       do i = 1,natl
c       tx=rmat(1,1)*y(i,1)+rmat(1,2)*y(i,2)+rmat(1,3)*y(i,3)
c       ty=rmat(2,1)*y(i,1)+rmat(2,2)*y(i,2)+rmat(2,3)*y(i,3)
c       tz=rmat(3,1)*y(i,1)+rmat(3,2)*y(i,2)+rmat(3,3)*y(i,3)
c       y(i,1) = tx
c       y(i,2) = ty
c       y(i,3) = tz
c       enddo
cc      call conformation(0)
c       do i = 1,natl
c               y(i,1) = cx + y(i,1)
c               y(i,2) = cy + y(i,2)
c               y(i,3) = cz + y(i,3)
c       enddo
cc      call conformation(0)
cc      print *,'AFTER ROTATION'
cc      do i = 1,natl
cc      write(6,FMT='(f8.3,2x,f8.3,2x,f8.3)')y(i,1),y(i,2),y(i,3)
cc      enddo
cc      stop
c************END OF EULER ANGLE ROTATION*********************



        return
        end


c*************Precalculate the interaction pair***********
        subroutine precal
        include 'mols.par'
        common /preset/ hbp(500000),hbl(500000),sp(500000),sl(500000),
     &          iscnt,ihcnt
        common /plp/ pat(mnatp), lat(maxatm),patnam(mnatp),
     &          presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),
     &          natp,tatyp(25,50),ntatp(25),px(mnatp,3),hb,steric
        common /par/ natom, ntor, nhb, ns, lres
        common /pplp/prresid(mnatp),matp
        common /comment/icomment
        character xatnam*4, xresnam*3,lat*2,pat*2,tatyp*2
        character patnam*4,presnam*3,tpres*3,tatnam*4
        character ltyp*4,ptyp*4
        
c       if(icomment.eq.1)print *,'precal-natp',natp
        if(icomment.eq.1)print *,'precal-matp',matp
c       print *,'HB SET'
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
c       print *,hbp(ik),hbl(ik),patnam(i)
        ik = ik + 1
        endif
        enddo
        enddo
c       stop
        ihcnt = ik-1

c       print *,'STERIC SET'
        ik = 1 
c       print *,natp,natom
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
c       print *,sp(ik),sl(ik),patnam(i),presnam(i)
c       print *,ptyp,ltyp
c	write(31,*)sp(ik),sl(ik),patnam(i),presnam(i)
c	write(31,*)ptyp,ltyp

        ik = ik + 1
        endif
        enddo
        enddo
        iscnt = ik-1

        return
c        stop
        end
c*************End pf Precalculation***********

c**************PLP Energy calculation**************
        subroutine eplp(plpe,hb,steric,tst)
        include 'mols.par'
cs      common/crdb/y(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /par/ natom, ntor, nhb, ns, lres
        common /preset/ hbp(500000),hbl(500000),sp(500000),sl(500000),
     &          iscnt,ihcnt
        common /propdb/ihflex,patname(mnatp),px1(mnatp),px2(mnatp),
     &  px3(mnatp),presname(mnatp),pcid(mnatp),presnum(mnatp)
        common /pplp/prresid(mnatp),matp

        common /plp/ pat(mnatp),lat(maxatm),patnam(mnatp),
     &  presnam(mnatp),pchid(mnatp),
     &  tpres(25),tatnam(25,50),natp,
     &  tatyp(25,50),ntatp(25),px(mnatp,3)
        common /comment/icomment
	common /tweight/xw1,xw2,xw3

        character xatnam*4, xresnam*3,lat*2,pat*2,tatyp*2
        character patnam*4,presnam*3,tpres*3,tatnam*4
        character ltyp*4,ptyp*4
c       integer hbp(90000),hbl(90000),sp(90000),sl(90000)
        real d,a,b,c,y1,ene,plpe,hb,steric
        integer tst

        character pcid*1,patname*4,presname*3,pchid*1
        integer presnum,ihflex,prresid
        character xx1*3,yy1*3,xx2*1,yy2*1,xx3*4,yy3*4
	real xw1,xw2,xw3

c       print *,'plp-natp',natp

        


        d = 0.0
        a = 0.0
        b = 0.0
        c = 0.0
        ene = 0.0
        hb = 0.0
        steric = 0.0

c****************************************************************

        if(icomment.eq.1)print *,'ihflex',ihflex,natp

        do i=1,ihflex
        do j=1,natp
        xx1=presname(i)
        xx2=pcid(i)
        xx3=patname(i)
        yy1=presnam(j)
        yy2=pchid(j)
        yy3=patnam(j)

c        if(tst.eq.2)then
c        print *,xx1(2:4),yy1(1:3),xx2(1:1),yy2(1:1),xx3(1:4),yy3(1:4)
c        endif

        if(xx1(1:3).eq.yy1(1:3).and.xx2(1:1).eq.yy2(1:1)
     &  .and.xx3(1:4).eq.yy3(1:4)
     &  .and.presnum(i).eq.prresid(j))then
        px(j,1)=px1(i)
        px(j,2)=px2(i)
        px(j,3)=px3(i)

c        if(tst.eq.2)then
c        print *,'pmols',j,presnum(i),prresid(j),px(j,1),px(j,2),px(j,3)
c        endif

        endif
        enddo
        enddo

c       do i=1,natp
c       print *,i,px(i,1),px(i,2),px(i,3)
c       enddo

        if(icomment.eq.1)print *,'ihcnt,iscnt',ihcnt,iscnt
        
cs        if(tst.eq.2)then
cs        do k=1,natom
cs        print *,'eplp',ylig(k,1),ylig(k,2),ylig(k,3)
cs        enddo
cs        endif
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
c       ene = ene + y1
        hb = hb + y1
       endif
        enddo
        do j = 1,iscnt
        ip = sp(j)
        il = sl(j)
        d = dist(px(ip,1),px(ip,2),px(ip,3),ylig(il,1),ylig(il,2),
     &   ylig(il,3))
c---------------------------------------------------------------
cs        if(icomment.eq.1.and.tst.eq.2.and.ip.eq.319)then
cs        print *,'clash-check',ip,il,d
cs        endif
c---------------------------------------------------------------
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
c       ene = ene + y1
        steric = steric + y1
      endif
      enddo
c	xw2=1.0
	if(icomment.eq.1)print *,'plpe-smmols-xw2',xw2
        plpe = (xw2*(hb+steric))

        if(icomment.eq.1)print *,'hb',hb,'steric',steric

cs      if(tst.eq.2) then
cs      endif
        return
        end
**************End of PLP****************************************
        function rgaffene(ie,je) ! energy calculation using GAFF

        include 'mols.par'
        common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &  emo(maxstr),emi(maxstr)
        common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /vectors/iv(maxpar,4)
        common /par/ natom, ntor, nhb, ns, lres
        common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common /pdbat/atom(maxatm),ele(maxatm)
        common /fnam/ if1,pf1,pf2,pf3,pf4,pf5,of1,of2,of3,of4,
     &  of5,of6
        common /energy/e_final(maxord,maxord)

        common /getrot/inrot,irotb(maxi,2),el,ilsrot(maxi,maxi),
     &  iatrot(maxi),rx(100,3),ire(maxi,maxi),ind(maxi,maxi),big
        common /left/ le(maxi,maxi),innd(maxi),ielenum(maxi),
     &  bonum(maxi)

        common /gen/ ibno,iat1(maxi),iat2(maxi),isno(maxi)
        common /string/ atsym(maxi),desc(maxi),attyp(maxi),
     &  str1(maxi),str2(maxi),str3,res,n1,n2,n3
        common /comment/icomment
        common /tweight/xw1,xw2,xw3
        common /native/nat,pepfile

        character*50 desc*30,attyp*10,atsym*6,str1*30,str2*30,
     &  str3*30,res*10,rmscr*35
        character*128 if1,pf1,pf2,pf3,pf4,pf5,pf6,of1,of2,of3,of4,
     &  of0,of5,of6
        
        dimension pp(maxpar)
        character str*80,t*20
        real etot,rdihed,rvdw,reel,xw1,xw2,xw3

10      format(a80)
11      format(a25)
20      format(i3,i3,a8)
21      format(i4,i4,i2,i2,i2)
22      format(a6)
23      format(a3)
30      format(3x,i4,1x,a6,2x,3f10.4,1x,a29)
40      format(3f8.3)
60      format(1x,i5,1x,i5,1x,i5,a10)
70      format(a10,7x,f10.5)

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

cd        if(nat.eq.0)then
        do i=1,natom
        write(50,30),innd(i),atsym(i),ylig(i,1),ylig(i,2),ylig(i,3),
     &  desc(i)
cs      write(*,30),innd(i),atsym(i),ylig(i,1),ylig(i,2),ylig(i,3),
cs   &  desc(i)
        enddo
cd        else if(nat.eq.1)then
cd         do i=1,natom
cd         write(50,30),innd(i),atsym(i),rx(i,1),rx(i,2),rx(i,3),
cd     &   desc(i)
cd         enddo
cd        endif

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

cs      read(113,*)rdihed
cs      read(113,*)rvdw
cs      read(113,*)reel

cs      etot= rdihed+rvdw+reel
        read(113,*)rene
        etot=rene
        close(unit=113)

        rgaffene = etot
        rmscr='rm output1'
        call system(rmscr)
c       xw1=1.0
        if(icomment.eq.1)print *,'rgaffene-xw1',xw1
        rgaffene = (xw1*(etot))
        if(icomment.eq.1)print *,'rgaffene',rgaffene
        return
        end
c******************************************************************
	function  rmmffene(ie,je) ! energy calculation using MMFF94

        include 'mols.par'
        common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &  emo(maxstr),emi(maxstr)
        common /ligcrda/xlig(maxatm,8)
        common /ligcrdb/ylig(maxatm,8)
        common /vectors/iv(maxpar,4)
        common /par/ natom, ntor, nhb, ns, lres
        common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common /pdbat/atom(maxatm),ele(maxatm)
        common /fnam/ if1,pf1,pf2,pf3,pf4,pf5,of1,of2,of3,of4,
     &  of5,of6
        common /energy/e_final(maxord,maxord)
 
        common /getrot/inrot,irotb(maxi,2),el,ilsrot(maxi,maxi),
     &  iatrot(maxi),rx(100,3),ire(maxi,maxi),ind(maxi,maxi),big
        common /left/ le(maxi,maxi),innd(maxi),ielenum(maxi),
     &  bonum(maxi)
        
        common /gen/ ibno,iat1(maxi),iat2(maxi),isno(maxi)
        common /string/ atsym(maxi),desc(maxi),attyp(maxi),
     &  str1(maxi),str2(maxi),str3,res,n1,n2,n3
	common /comment/icomment
	common /tweight/xw1,xw2,xw3
        common /native/nat,pepfile

        character*50 desc*30,attyp*10,atsym*6,str1*30,str2*30,
     &  str3*30,res*10,rmscr*35
        character*128 if1,pf1,pf2,pf3,pf4,pf5,pf6,of1,of2,of3,of4,
     &  of0,of5,of6

        dimension pp(maxpar)
        character str*80,t*20
        real etot,rdihed,rvdw,reel,xw1,xw2,xw3
 
10      format(a80)
11      format(a25)
20      format(i3,i3,a8)
21      format(i4,i4,i2,i2,i2)
22      format(a6)
23      format(a3)
30      format(3x,i4,1x,a6,2x,3f10.4,1x,a29)
40      format(3f8.3)
60      format(1x,i5,1x,i5,1x,i5,a10)
70	format(a10,7x,f10.5)
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

cd        if(nat.eq.0)then
	do i=1,natom
        write(50,30),innd(i),atsym(i),ylig(i,1),ylig(i,2),ylig(i,3),
     &  desc(i)
        enddo
cd        else if(nat.eq.1)then
cd        do i=1,natom
cd        write(50,30),innd(i),atsym(i),rx(i,1),rx(i,2),rx(i,3),
cd     &  desc(i)
cd        enddo
cd        endif

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
c	xw1=1.0
	if(icomment.eq.1)print *,'rmmffene-xw1',xw1
	rmmffene = (xw1*(etot))
	return
	end
c**********************************************************************
c-----------------------------------------------------------------------
        function ddihedr(i1,i2,i3,i4)
        include 'mols.par'
        common/rcrda/xx_one(maxatm,3)
        common /comment/icomment
        real x(maxatm,3)
        acdc=((180.0*7.0)/22.0)
        one=1.d0
!-----------------------------------        
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
!-----------------------------------
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
        write (*,'(a,4i5)')' dihedr> Error in coordinates of atoms #: '
     #                     ,i1,i2,i3,i4
        stop
      endif
      end
c------------------------------------------------------------------------- 


