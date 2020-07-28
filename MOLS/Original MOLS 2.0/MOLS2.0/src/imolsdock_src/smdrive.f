c-----------------------------------------------------------------------------
c     library name : smdrive.f   

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

        program MAIN
c	program smdrive
        include 'mols.par'
c-------------------------------------------------------------------
	common /comment/icomment
	common /ligconv/lfile
	common /prepare/receptor
	common /find/fdis
	common /rctl/iscopt
	common /sminimize/itmaxs
	common /scang/frange,rang
	common /native/nat,pepfile
	common /samfnam/sf0,sf1,sf2
	common /cen/bx,by,bz
	common /tweight/xw1,xw2,xw3
	common /conf/path,mname
        common /par/natom,ntor,nhb,ns,lres,is,ie
        common /tor/u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common /ctl/iopt,ipf,iff,icint,fseq(maxres),ifopt,ilopt
        common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6
        common /mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
        common /plp/pat(mnatp),lat(maxatm),patnam(mnatp),
     &  presnam(mnatp),pchid(mnatp),tpres(25),tatnam(25,50),
     &  natp,tatyp(25,50),ntatp(25),px(mnatp,3)
c--------------------------------------------------------------------
	integer is,ie,icomment,ilopt,iopt,irun,iff,ifopt,
     &  ipf,ifmin,iscopt,ifind,nat,itmaxs,inp
	real bx,by,bz,frange,fdis,xw1,xw2,xw3

	character seq*(maxres)
	character*128 path,mname,sf0,sf1,sf2,if1,if2,pf1,pf2,pf3,pf6,
     &  pf7,of0,of1,of2,of3,of5,of6,receptor,pepfile,lfile
	character*12 c1,c2
	character*24 date,job
c--------------------------------------------------------------------
10	format(5x,a27,1x,f9.3,2x,f9.3,2x,f9.3)
15	format(5x,a27,f6.1)
16	format(5x,a27,f4.1)
17	format(5x,a6,f4.1,1x,a4,f4.1,1x,a7,f4.1)
20	format(5x,a27,2x,i4,2x,i4)
c--------------------------------------------------------------------
        open(unit=18,file='input',status='old')
        read(18,*)is,ie
	if(is.gt.ie) STOP 'Error: Starting Structure Number higher'
	rewind(18)
        read(18,*)c1,c2
        close(18)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	open(unit=30,file='user.inp',status='old')
        read(30,*)icomment
        read(30,'(A)') path
        read(30,*) job
        mname=job(:LNBLNK(job))//'_'//
     &  c1(:LNBLNK(c1))//'-'//
     &  c2(:LNBLNK(c2))
c***********Output file names***************************************
        sf0 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'mols.log'
        sf1 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'mini.log'
        sf2 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'ene.log'
        if1 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'inp.par'
        if2 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'pinp.par'
        pf1 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'inp.pdb'
        pf2 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'mols.mol2'
        pf3 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'mini.mol2'
        pf6 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'pep-mols.pdb'
        pf7 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'pep-mini.pdb'
        of0 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'inf.log'
        of1 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'mols.out'
        of2 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'mini.out'
        of3 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'inp.info'
        of5 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'mini_complex.pdb'
        of6 = path(:LNBLNK(path)) // '/' //
     &        mname(:LNBLNK(mname)) // '_' // 'mols_complex.pdb'
c*********************************************************************

	read(30,*)ilopt
c---------------------------------------------------------------
	if(ilopt.eq.1)then	!1: Peptide Ligand		! 
	read(30,'(A)')seq					! 
	iopt=2		!Back-Bone + Side-Chain			!  
	irun=3 			!Peptide run upto minimiz	!	  
	elseif(ilopt.eq.2)then					!  
	read(30,'(A)')lfile					! 
	call ligandconv()					! 
	endif							! 
c---------------------------------------------------------------
	read(30,*)iff 	!Peptide force-field (1:AMBER)	  
				!Small-molecule (1:MMFF94;2:GAFF) 
	read(30,*)bx,by,bz 	!Grid Centre			  
	read(30,'(A)')receptor
c-------------------------------------------------------
	read(30,*)ifopt,frange				!
							!
        if(ifopt.eq.1)then !rigid-receptor docking	!
        ipf=0						!
        ifmin=0						!
        endif						!
							!
        if(ifopt.eq.2)then!flexible-receptor docking	!
        ipf=1						!
        ifmin=1						!
        endif 						!
c-------------------------------------------------------
	read(30,*)iscopt
        read(30,*)ifind,fdis
	read(30,*)xw1,xw2,xw3!weights
	read(30,*)nat!default-0;native energy-1		

	if(nat.eq.0)read(30,*) itmaxs
	if(nat.eq.1)read(30,'(A)') pepfile

	close(30)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c-----Details for inf.log--------------------------------------------
      open(unit=31,file=of0,status='unknown')
      WRITE (31,1)
      WRITE (*,1)

1     FORMAT (//
     *10X,'#########################################################'/
     *10X,'#                                                       #'/
     *10X,'#                     iMOLSDOCK                         #'/
     *10X,'#                                                       #'/
     *10X,'#          CAS in Crystallography and Biophysics,       #'/
     *10X,'#    University of Madras, Chennai - 600 025, INDIA     #'/
     *10X,'#                                                       #'/
     *10X,'#########################################################'/)
 
        call fdate(date)

        WRITE (*,'(/19X,3A/)') 'Job started on ',date
        WRITE (31,'(/19X,3A/)') 'Job started on ',date

	WRITE(31,'(/5X,3A/)') 'Project Name              :',job

        if(ilopt.eq.1)then
        write(31,'(/5X,3A/)') 'Docking                   :
     &Peptide-Protein'
        write(31,'(/5X,3A/)') 'Peptide Sequence          :',
     &seq(:LNBLNK(seq))
        else if(ilopt.eq.2) then
        write(31,'(/5X,3A/)') 'Docking                   :
     &Small molecule-Protein'
        write(31,'(/5X,3A/)') 'Ligand File               :',lfile
        endif

        write(31,'(/5X,3A/)') 'Protein File              :',
     &receptor(:LNBLNK(receptor))
	write(31,'(A)')
        write(31,10) 'Grid Centre               :',bx,by,bz

        if(ifopt.eq.1)then
        write(31,'(/5X,3A/)') 'Protein Flexibility       :Disabled'
        else if(ifopt.eq.2)then
        write(31,'(/5X,3A/)') 'Protein Flexibility       :Enabled'
        write(31,'(A)')
        write(31,15) 'Protein Flexibility Range :',frange
        endif

        if(ifopt.eq.2.and.ifind.eq.1)then
        write(31,'(/5X,3A/)') 'Flexible residues         :Auto-find 
     &using Fpocket2.0'
        write(31,16) 'Auto-find scan range      :',fdis
        else if(ifopt.eq.2.and.ifind.eq.0)then
        write(31,'(/5X,3A/)') 'Flexible residues         :Manually 
     &defined'
        endif

	if(ilopt.eq.1.and.iff.eq.1)then
        write(31,'(/5X,3A/)')'Intra-Pep:AMBER,Peptide-Protein:PLP,
     &Intra-Protein:AMBER'
        write(31,'(A)')
        write(31,17)'AMBER: ',xw1,'PLP:',xw2,'AMBER :',xw3
        else if(ilopt.eq.1.and.iff.eq.2)then
        write(31,'(/5X,3A/)')'Intra-Pep:ECEPP,Peptide-Protein:PLP,
     &Intra-Protein:AMBER'
        write(31,'(A)')
        write(31,17)'ECEPP: ',xw1,'PLP:',xw2,'AMBER :',xw3
        else if(ilopt.eq.2.and.iff.eq.1)then
        write(31,'(/5X,3A/)')'Intra-Lig:MMFF94,Peptide-Protein:PLP,
     &Intra-Protein:AMBER'
        write(31,'(A)')
        write(31,17)'MMFF94:',xw1,'PLP:',xw2,'AMBER :',xw3
        else if(ilopt.eq.2.and.iff.eq.2)then
        write(31,'(/5X,3A/)')'Intra-Lig:GAFF,Peptide-Protein:PLP,
     &Intra-Protein:AMBER'
        write(31,'(A)')
        write(31,17)'GAFF:  ',xw1,'PLP:',xw2,'AMBER :',xw3
        endif

        write(31,'(A)')
        write(31,20) 'Structures to be generated:',is,ie
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	call preparepro()!Receptor Protein Protonation
	call iread_pdb(ifopt)

	if(ilopt.eq.1)then
	call ppdbgen(seq,fseq)!Generate Peptide structure

	  if(nat.eq.1)then
	    call creducer(2)
	    call npamppar(iopt)
	  else if(nat.eq.0)then
           if(iff.eq.1)call pamppar(iopt)
	  endif
	call pvarinit
	else if(ilopt.eq.2)then
	  call pdbgen
	  call varinit
	endif

	if(ifind.eq.1)call findflex()
	if(ifopt.eq.2)call recflex()
	call flush(6)

	write(31,'(A)')
	write(31,'(/5X,3A/)')'RESULTS:'
	if(ilopt.eq.1)then
	 do inp = is,ie
	  call mols(inp)
	  call flush(6)
	  call pconformation(inp,1)
	  call flush(6)
          call conformation(1,inp,1)
 	  call flush(6)
	  call pminimiz(inp)
	  call flush(6)
	  call pconformation(inp,2)
	  call flush(6)
	  call minimiz(1,inp)
	  call flush(6)
	  call updateinp(inp,ie)
	 enddo
	else if(ilopt.eq.2)then
	 do inp = is,ie
	  call mols(inp)
	  call flush(6)
	  call conformation(1,inp,2)
	  call flush(6)
	  call sminimiz(inp)
	  call flush(6)
 	  call minimiz(2,inp)
	  call flush(6)
	  call updateinp(inp,ie)
	 enddo
	endif

	call fdate(date)
        WRITE (*,'(/19X,3A/)') 'Job ended on ', date
        WRITE (31,'(/19X,3A/)') 'Job ended on ', date

	end
c=======================================================
	subroutine updateinp(inp,ie)
	open(unit=108,file='input',status='unknown')
        write(108,*)inp+1,ie
        close(108)
	return
	end
c=======================================================
	subroutine ligandconv()
	parameter (maxn=2000)
	common /ligconv/lfile

10      format(a80)
11      format(a6)

	character*80 str,convert
        character*128 lfile

	integer i,maxn

        open(unit=1,file=lfile,status='unknown')
        do i=1,maxn
        read(1,10,end=99),str
        enddo
99      write(*,*),'line number',i
        close(1)
        latom=i-1

        if(icomment.eq.1) print *,latom

        open(unit=2,file=lfile,status='unknown')
        open(unit=3,file='ligand.mol',status='unknown')
        do k=1,latom
        read(2,10)str
        if(k.eq.1)then
        write(3,11)'ligand'
        else
        write(3,10)str
        endif
        enddo
        close(3)
        close(2)

        convert='sh fileform.sh'
        call system(convert)

	return
	end
c========================================================
	subroutine preparepro()
	include 'mols.par'
        parameter (maxm=20000)
	common /comment/icomment
	common /prepare/receptor
	common /profile/pfile
	common /fnam/if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6

	integer i,lpres
	character*128 if1,if2,pf1,pf2,pf3,pf6,pf7,of0,of1,of2,of3,
     &  of5,of6,receptor,protfile,prot,pfile
	character str*80,str1*16,str2*1,str3*63,prep*80

10      format(a80)
11      format(a16,a1,a63)
12      format(a16,1x,a63)

        protfile=receptor(:LNBLNK(receptor))
        lpres=0
        do i=1,128
        if(protfile(i:i).ne.' ')then
        lpres=lpres+1
        endif
        enddo
        prot=protfile(1:lpres-4)
        pfile=prot(:LNBLNK(prot))//'_'//'h.pdb'


	open(unit=1,file=protfile,status='unknown')
        do i=0,maxm
        read(1,10,end=19),str
        enddo
19      if(icomment.eq.1)write(*,*),'protein-atoms',i
        close(1)
c-------filter protein ATOMs----------------------------------------------      
        itemp=0
        open(unit=2,file=protfile,status='unknown')
        open(unit=3,file='temp1',status='unknown')
        do j=1,i
        read(2,10)str
        if(str(1:4).eq.'ATOM')then
        itemp=itemp+1
        write(3,10)str
        endif
        enddo
        close(3)
        close(2)
c------check isoforms-----------------------------------------------------      
        open(unit=4,file='temp1',status='unknown')
        open(unit=5,file='temp2',status='unknown')
        do k=1,itemp
        read(4,11)str1,str2,str3
        if(str2(1:1).eq.' '.or.str2(1:1).eq.'A')then
        write(5,12)str1,str3
        endif
        enddo
c------protonate protein-------------------------------------------------
cs      prep='chmod 777 reduce'
cs      call system(prep)
        prep='./reduce -DB red_dict.txt -Quiet -keep temp2>temp2_h.pdb'
        call system(prep)
        prep='grep ATOM temp2_h.pdb > temp3_h.pdb'
        call system(prep)
c-------copy protonated protein to pfile---------------------------------
        open(unit=6,file='temp3_h.pdb',status='unknown')
        do il=0,maxm
        read(6,10,end=29),str
        enddo
29	if(icomment.eq.1)write(*,*),'pfile-protein-atoms',il
        close(6)

        open(unit=7,file='temp3_h.pdb',status='unknown')
        open(unit=8,file=pfile,status='unknown')
        do im=1,il
        read(7,10)str
        write(8,10)str
        enddo
        close(8)
        close(7)

        return
        end
c*************************************************************************
!!Subroutine to find the number of atoms in the receptor (protein) molecule.
        subroutine iread_pdb(ifopt)
        include 'mols.par'
        parameter (maxm=20000)

	common /comment/icomment
	common /patom/ipatom
	common /profile/pfile

        character pfile*128,str*80,atm*4,str1*2,str2*67

        integer opatom,opatm,ll

10      format(a80)
30      format(a4,a2,i5,a67)

        open(unit=1,file=pfile,status='unknown')
        do i=1,maxm
        read(1,10,end=99),str
        enddo
99      if(icomment.eq.1)write(*,*),'line number',i
        close(unit=1)

        opatom=i-1 !no. of lines in original receptor PDB

        open(unit=1,file=pfile,status='unknown')
        open(unit=2,file='pfile_orig.pdb',status='unknown')
        do k=1,opatom
        read(1,10)str
        write(2,10)str
        enddo
        close(2)
        close(1)

        ll=0
        open(unit=1,file=pfile,status='unknown')
        open(unit=2,file='pfile_mod.pdb',status='unknown')
        do l=1,opatom
        read(1,30)atm,str1,is,str2
        if(atm.eq.'ATOM'.and.str2(3:5).ne.'OXT')then
        ll=ll+1
        write(2,30)atm,str1,ll,str2
        endif
        enddo
        close(2)
        close(1)

        ipatom=ll   !#1-- no. of lines in modified receptor PDB
                    !#2--modified receptor PDB[pfile_mod.pdb] has only 
                    !'ATOM' of the original receptor PDB[pfile_orig.pdb] 
        if(ipatom.ge.mnatp)then
        print *,'No. of Protein atoms exceeding mnatp in mols.par'
        STOP
        endif

        if(icomment.eq.1) write(*,*)'ipatom',ipatom

        call creducer(1)!convert receptor protein to MOLS atomtype format

!Q:why does the following loop copies the modified protein into pfile1.pdb?
!A:while writing complex result files,'pfile1.pdb' is the file that
!has the protein structure. In case of 'flexible
!docking','pfile1.pdb' is generated in subroutine 'flexall()'. For
!'rigid docking' protein structure remains static therefore 'pfile1.pdb'
!is generated here in the initial stage.        

        if(ifopt.eq.1)then
        if(icomment.eq.1) print *,"rigid docking"
        open(unit=1,file='rec.pdb',status='unknown')
        open(unit=2,file='pfile1.pdb',status='unknown')
        do m=1,ipatom
        read(1,10)str
        write(2,10)str
        enddo
        close(2)
        close(1)
        else
        if(icomment.eq.1) print *,"flexible docking"
        endif

        return
        end
c***********************************************************************
