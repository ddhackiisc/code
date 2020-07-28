c      implicit integer*4 (i-n)
c      implicit real*8 (a-h,o-z)	
       
      parameter (izero=0, ione=1,mxtyat = 18,mxhbdo=4,mxhbac=6,
     #           zero=0.d0, one=1.d0)

c  Energy parametrization for SCHERAGA set
c
c  ehm   - 1.5*(elem. charge)/(Planck's const./2Pi * electr. mass)
c  atpl  - atomic polarizabilities
c  efel  - empirical numbers of effective electrons per atom
c  emin  - minima of pairwise Lennard-Jones terms of like atoms
c  rmin  - corresponding pair-distance of energy minima

      common /ecp/ ehm,atpl(mxtyat),efel(mxtyat),emin(mxtyat),
     #    rmin(mxtyat),chb_s(mxhbdo,mxhbac),ahb_s(mxhbdo,mxhbac)
      logical do_f(mxtyat),ac_f(mxtyat)
c ***************************************************************

c  Atom types ------------------------------------------------------------
c                                    Original types  -Scheraga:  -Flex:
c  H  1 - with aliphatic carbon                               1      12
c     2 - with aromatic carbon                                3      13
c     3 - with non-sp3 types of nitrogen                      2       1
c     4 - with sp3-hybr. nitrogen                             2       2
c     5 - with oxygen                                         4       1
c     6 - with sulfur                                         3(was 5)1
c  C  7 - sp3-hybr. carbon                                    6,9     3
c     8 - sp2-carbon (carbonyl,carboxyl,carboxylate)          7,11    4
c     9 - aromatic carbon                                     8,10    4
c  O 10 - hydroxyl, ester oxygen (inc. water)                 18,19   8
c    11 - carbonyl oxygen                                     17      9
c    12 - carboxylate oxygen                                  18,19  10
c  N 13 - aliph. nitrogen with 0/1 hydrogen & charged N       13-15   6
c    14 - nitrogen with two hydrogens                         13-15   5
c    15 - all other nitrogens (+ sp2-hybrid. in heteroc.)     13-15   7
c  S 16 - any sulfur                                          20,21  18,19
c  H 17 - H-delta of Pro, Hyp of ECEPP/3 dataset               5(new) -
c  C 18 - C-delta of Pro, Hyp of ECEPP/3 dataset              12(new) -
c  Classes for torsional potential ---------------------------------------
c
c   1 : 'Omega' = C'(pept.)-N(pept.)  [Cpept-Npept]
c   2 : 'Phi'   = N(pept.)-C(sp3)     [C4-Npept]
c   3 : 'Psi'   = C(sp3)-C'(pept.)    [C4-Cpept]
c   4 : 'Chi1'  = C(sp3)-C(sp3)       [C4-C4]
c   5 : C(sp3)-OH (Hydroxyl)          [C4-OH]
c   6 : C(sp3)-NH2                    [C4-NH2]
c   7 : C(sp3)-NH3+                   [C4-NH3+]
c   8 : C(sp3)-NH-(guanidyl)          [C4-NHX]
c   9 : C(sp3)-COOH(carboxyl)         [C4-COO]
c  10 : C(sp3)-COO-(carboxylate)      [C4-COO]
c  11 : C(sp3)-CO(sp2 of amide)       [C4-Cpept]
c  12 : C(sp3)-C(aromatic ring)       [C4-C3]
c  13 : C(sp3)-S                      [C4-SC4]
c  14 : C(sp3)-SH                     [C4-SH]
c  15 : C(aromatic ring)-OH           [C3-OH]
c ________________________________________________ "rigid" torsions:
c  16 : C(carboxyl)-OH                [C3-OH]
c  17 : -NH-C(sp2 of guanidyl)        [C3-NHX]
c  18 : -C(sp3)-NH2 (guanidyl)        [not in Flex]
c  19 : -C(sp3)-NH2 (amide)           [Cpept-Npept]

      data conv/332.d0/  ! to convert electrost. energy into [kcal/mole]

c ------------------------- ECEPP/3 potential --------------------------------
c 1) Momany F.A McGuire R.F Burgess A.W Scheraga H.A J Phys Chem v79 2361-2381
c    1975
c 2) Nemethy G Pottle M.S Scheraga H.A, J Phys Chem v87 1883-1887 1983
c 3) Sippl M.J Nemethy G Scheraga H.A J Phys Chem v88 6231-6233 1984
c 4) Nemethy G Gibson K.D Palmer K.A Yoon C.N Paterlini G Zagari A Rumsey S
c    Scheraga H.A J Phys Chem v96 6472-6484 1992
c ----------------------------------------------------------------------------

      data eps_s/2.d0/  ! Distance-INdependent diel. constant
c     data eps_s/6.d0/  ! Distance-INdependent diel. constant
      data plt/78.d0/,  slp/0.3d0/   ! Parameters for Epsilon(R)

      data ehm /362.55d0/  !  Angstrom**2/3 * kcal / mol  ! from KONF90
cc      data ehm /362.09561409d0/  !  Angstrom**2/3 * kcal / mol
c From:
c   1.5
c * elementary charge       = 4.80325   *e+2  Angstrom**3/2 * g**1/2 * s**(-1)
c * Planck's constant/2*Pi  = 1.0545887 *e-34 Joule * s
c * Avogadro's number       = 6.022045  *e+23 mol**(-1)
c / sqrt (mass of electron) = sqrt (9.109534  *e-28 g )
c / thermal equivalent      = 4.1868    *e+3  Joule * kcal**(-1)
cc      data ehm /362.36d0/         ! calculated using Tab II in ref. 2
cc      data 1/ehm /2.757670d-3/    ! 3*sqrt(m)/(2*e*h) taken from ICM

c ---------------------- atomic polarizabilties (*100,[Angstrom**3])
c                1   2   3   4   5   6   7    8    9  10  11  12
      data atpl/42.,42.,42.,42.,42.,42.,93.,151.,115.,59.,84.,59.,
c               13  14  15   16  17  18
     #          93.,93.,93.,220.,42.,93./
c ---------------------- effective numbers of electrons (*100,ref. 2)
c                1   2   3   4   5   6    7    8    9   10   11   12
      data efel/85.,85.,85.,85.,85.,85.,520.,520.,520.,700.,700.,700.,
c                13   14   15    16  17   18
     #          610.,610.,610.,1480.,85.,520./
c ------------------------- min. pairwise 6-12 energy (*1000,[kcal/mol])
c                1   2   3   4   5   6   7    8   9  10   11  12
      data emin/37.,36.,61.,61.,44.,36.,38.,140.,99.,94.,200.,94.,
c                13   14   15   16  17  18
     #          107.,107.,107.,223.,99.,38./
c ---------------------------- opt. pairwise distance (*100,[Angstrom])
c                 1    2    3    4    5    6    7    8    9   10   11
      data rmin/292.,293.,268.,268.,283.,293.,412.,374.,370.,324.,312.,
c                12   13   14   15   16   17   18
     #          324.,351.,351.,351.,415.,248.,412./
c --------------------------------- HB-parameters (/1000,attraction)
      data chb_s/2624.,2624.,4610.,.0,  ! given as:
     #           4014.,4014.,5783.,.0,  !  (ac_typ x do_typ)
     #           2624.,2624.,4610.,.0,  ! to be used:
     #           8244.,8244.,8244.,.0,  !  (DO_typ x AC_typ)
     #           8244.,8244.,8244.,.0,  ! i.e.:
     #           8244.,8244.,8244.,.0/  !  ( 3-5 x 10-15 )
c --------------------------------- HB-parameters (/1000,repulsion)
      data ahb_s/ 5890., 5890.,11220.,.0,
     #           12040.,12040.,16583.,.0,  ! 13344 -> 16583 = Ref. 3
     #            5890., 5890.,11220.,.0,
     #           32897.,32897.,32897.,.0,
     #           32897.,32897.,32897.,.0,
     #           32897.,32897.,32897.,.0/
c ---------------------------------------------- Hydrogen-bond donors
c                  1       2       3      4      5      6
      data do_f/.false.,.false.,.true.,.true.,.true.,.true.,
c                  7       8       9      10      11      12
     #          .false.,.false.,.false.,.false.,.false.,.false.,
c                 13      14      15      16      17      18
     #          .false.,.false.,.false.,.false.,.false.,.false./
c -------------------------------------------- Hydrogen-bond acceptors
c                  1       2       3       4       5       6
      data ac_f/.false.,.false.,.false.,.false.,.false.,.false.,
c                  7       8       9      10     11     12
     #          .false.,.false.,.false.,.true.,.true.,.true.,
c                 13      14      15     16     17      18
     #          .true.,.false.,.true.,.true.,.false.,.false./
c --------------------------------- HB-parameters (/1000,attraction)
