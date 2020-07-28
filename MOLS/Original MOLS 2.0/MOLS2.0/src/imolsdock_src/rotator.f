c-----------------------------------------------------------------------------
c     library name : rotator.f   

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



c-------Peptide-Protein Docking-----------------------------
c*****Rotate ligand by angles Rx, Ry and Rz about***********
c*****       x, y and z axiz respectively        ***********
        subroutine protate(bx,by,bz,rotrang,theta,psi,natl)
        include 'mols.par'

        common /pepcrdb/ypep(maxatm,8)
        common /native/nat,pepfile

        real cx,cy,cz,theta1,psi1,rotrang1,a,b,c,d,e,f,
     &  ad,bd,rmat(3,3),w1,w2,qx,qy,qz,tw,tx,ty,tz,tw1,
     &  tx1,ax,ay,az,sx,sy,sz,s,t
        integer natl,r,nat

        rotrang1 = rotrang
        theta1 = theta
        psi1 = psi

        if(nat.eq.1)goto 199

        theta1=theta1*3.1415927/180.0
        psi1=(psi1*3.1415927/180.0)/2
        rotrang1=rotrang1*3.1415927/180.0

        r = 1
        w2 = 0.0
        sx = 0.0
        sy = 0.0
        sz = 0.0

c       CREATE POINTS ON THE UNIT SPHERE USING SPHERICAL COORD

        sx = r*cos(theta1)*sin(psi1)
        sy = r*sin(theta1)*sin(psi1)
        sz = r*cos(psi1)

        vdis = dist(sx,sy,sz,0.0,0.0,0.0)

c       CREATE THE AXIS OF ROTATION

        ax = sx/vdis
        ay = sy/vdis
        az = sz/vdis

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

        call pcalcent(natl,cx,cy,cz)
        do i = 1,natl
                ypep(i,1) = cx-ypep(i,1)
                ypep(i,2) = cy-ypep(i,2)
                ypep(i,3) = cz-ypep(i,3)
        enddo
        
        do k=1,natl
        tx=rmat(1,1)*ypep(k,1)+rmat(1,2)*ypep(k,2)+rmat(1,3)*ypep(k,3)
        ty=rmat(2,1)*ypep(k,1)+rmat(2,2)*ypep(k,2)+rmat(2,3)*ypep(k,3)
        tz=rmat(3,1)*ypep(k,1)+rmat(3,2)*ypep(k,2)+rmat(3,3)*ypep(k,3)
        ypep(k,1) = tx
        ypep(k,2) = ty
        ypep(k,3) = tz
        enddo

199     return
        end
c*******End of Rotate*************************************************

c**************Calculate ligand centroid******************************
        subroutine pcalcent(natl,cx,cy,cz)
        include 'mols.par'
        common/pepcrdb/ypep(maxatm,8)
        integer natl
        real cx,cy,cz,cx1,cy1,cz1
        cx=0.0
        cy=0.0
        cz=0.0

        do i = 1,natl
                cx = cx + ypep(i,1)
                cy = cy + ypep(i,2)
                cz = cz + ypep(i,3)
        enddo
        cx1 = cx / natl
        cy1 = cy / natl
        cz1 = cz / natl
        cx = cx1
        cy = cy1
        cz = cz1
        return
        end
c**************End of ligand centroid*********************************

c*******Translate ligand inside the grid as guided by MOLS************

        subroutine ptranslate(bx,by,bz,tx,ty,tz,natl)
        include 'mols.par'

        common /pepcrdb/ypep(maxatm,8)
        common /native/nat,pepfile

        real cx,cy,cz,dx,dy,dz
        integer natl,nat

        if(nat.eq.1.)goto 299
c*******USED FOR MOVING CENTROID AGAIN TO ORGIN***********************
        call pcalcent(natl,cx,cy,cz)

        do i = 1,natl
         ypep(i,1) = cx-ypep(i,1)
         ypep(i,2) = cy-ypep(i,2)
         ypep(i,3) = cz-ypep(i,3)
        enddo
c*******USED FOR MOVING CENTROID AGAIN TO ORGIN***********************
        do i = 1,natl
         ypep(i,1) = ypep(i,1)+(bx)+tx
         ypep(i,2) = ypep(i,2)+(by)+ty
         ypep(i,3) = ypep(i,3)+(bz)+tz
        enddo

299     return
        end
c*************End of Translation**************************************

c------Protein-Ligand Docking-----------------------------------------

c-------Rotate ligand by angles Rx, Ry and Rz about-------------------
c-------     x, y and z axiz respectively        ---------------------
        subroutine rotate(bx,by,bz,rotrang,theta,psi,natl)

        include 'mols.par'

        common /ligcrdb/ylig(maxatm,8)

        real cx,cy,cz,theta1,psi1,rotrang1,a,b,c,d,e,f,ad,bd,
     &  rmat(3,3),w1,w2,qx,qy,qz,tw,tx,ty,tz,tw1,tx1,ax,ay,az,
     &  sx,sy,sz,s,t

        integer natl,r

        rotrang1 = rotrang
        theta1 = theta
        psi1 = psi
        theta1=theta1*3.1415927/180.0
        psi1=(psi1*3.1415927/180.0)/2

c       CONVERT DEGREES TO RADIAN

        rotrang1=rotrang1*3.1415927/180.0

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

        vdis = dist(sx,sy,sz,0.0,0.0,0.0)

        ax = sx/vdis
        ay = sy/vdis
        az = sz/vdis

c*********ANGLE/AXIS ROTATION******************************************
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

        call calcent(natl,cx,cy,cz)
        do i = 1,natl

                ylig(i,1) = cx-ylig(i,1)
                ylig(i,2) = cy-ylig(i,2)
                ylig(i,3) = cz-ylig(i,3)
        enddo

        do k=1,natl

        tx=rmat(1,1)*ylig(k,1)+rmat(1,2)*ylig(k,2)+rmat(1,3)*ylig(k,3)
        ty=rmat(2,1)*ylig(k,1)+rmat(2,2)*ylig(k,2)+rmat(2,3)*ylig(k,3)
        tz=rmat(3,1)*ylig(k,1)+rmat(3,2)*ylig(k,2)+rmat(3,3)*ylig(k,3)

        ylig(k,1) = tx
        ylig(k,2) = ty
        ylig(k,3) = tz
        enddo

        return
        end
c--------------------------------------------------------------------
c**************Calculate ligand centroid*****************************
        subroutine calcent(natl,cx,cy,cz)

        include 'mols.par'

        common /ligcrdb/ylig(maxatm,8)

        integer natl
        real cx,cy,cz,cx1,cy1,cz1
        
	cx=0.0
        cy=0.0
        cz=0.0

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
c**************End of ligand centroid********************************

c*******Translate ligand inside the grid as guided by MOLS***********
        subroutine translate(bx,by,bz,tx,ty,tz,natl)

        include 'mols.par'

        common /ligcrdb/ylig(maxatm,8)

        real cx,cy,cz,dx,dy,dz
	integer natl

c-------USED FOR MOVING CENTROID AGAIN TO ORGIN----------------------

        call calcent(natl,cx,cy,cz)

        do i = 1,natl
         ylig(i,1) = cx-ylig(i,1)
         ylig(i,2) = cy-ylig(i,2)
         ylig(i,3) = cz-ylig(i,3)
        enddo

c-------USED FOR MOVING CENTROID AGAIN TO ORGIN---------------------

        do i = 1,natl
         ylig(i,1) = ylig(i,1)+(bx)+tx
         ylig(i,2) = ylig(i,2)+(by)+ty
         ylig(i,3) = ylig(i,3)+(bz)+tz
       enddo

        return
        end
c*************End of Translation************************************






