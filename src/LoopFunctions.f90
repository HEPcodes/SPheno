Module LoopFunctions
! comments
! This module collects the functions necessary for 1- and 2- loop
! calculations. 

! load modules
Use Control
Use Mathematics, Only: Li2, CLi2, Delt
! load modules

! interfaces
 Interface Log1minusX
  Module Procedure Log1minusX_r, Log1minusX_c
 End Interface

 Interface Log1minusXpXn
  Module Procedure Log1minusXpXn_r, Log1minusXpXn_c
 End Interface

 Interface roots
  Module Procedure roots_r, roots_r2, roots_c2
 End Interface
! interfaces

! data types
Type C_data
 Real(dp) :: mi(3)
 Real(dp) :: pi(3)
 Complex(dp) :: Ci(13)
End Type

! private variables
Real(dp), Private :: mudim2=1._dp, divergence=0._dp, lambda2 = 1._dp
!
!       xloss:  factor that the final result of a subtraction can be
!               smaller than the terms without warning (default 1/8)
!       precx:  precision of real numbers, determined at runtime by
!               ffinit (IEEE: 4.e-16)
!       precc:  same for complex numbers
!       xalogm: smallest real number of which a log can be taken,
!               determined at runtime by ffinit (IEEE: 2.e-308)
!       xclogm: same for complex.
!       xalog2: xalogm**2
!       xclog2: xclogm**2
!       reqprc: not used
!	   the precision wanted in the complex D0 (and hence E0) when
!	   nschem=7, these are calculated via Taylor exoansion in the real
!	   one and hence expensive.
!       pi:     pi
!       pi6:    pi**2/6
!       pi12:   pi**2/12
!       xlg2:   log(2)
!       bf:     factors in the expansion of dilog (~Bernouilli numbers)
!       xninv:  1/n
!       xn2inv: 1/n**2
!       xinfac: 1/n!
!       fpij2:  vi.vj for 2point function 1-2: si, 3-3:  pi
!       fpij3:  vi.vj for 3point function 1-3: si, 4-6:  pi
!       fpij4:  vi.vj for 4point function 1-4: si, 5-10: pi
!       fpij5:  vi.vj for 5point function 1-5: si, 6-15: pi
!       fpij6:  vi.vj for 6point function 1-6: si, 7-21: pi
!       fdel2:  del2 = delta_(p1,p2)^(p1,p2) = p1^2.p2^2 - p1.p2^2 in C0
!       fdel3:  del3 = delta_(p1,p2,p3)^(p1,p2,p3) in D0
!       fdel4s: del4s = delta_(s1,s2,s3,s4)^(s1,s2,s3,s4) in D0
!       fdel4:  del4 = delta_(p1,p2,p3,p4)^(p1,p2,p3,p4) in E0
!       fdl3i:  del3i = delta_(pj,pk,pl)^(pj,pk,pl) in E0, D0 without si
!       fdl4si: dl4si = del4s in E0, D0 without si
!       fdl3ij: same in F0 without si and sj.
!       fd4sij: dl4si = del4s in E0, D0 without si
!       fdl4i:  delta4 in F0 without si.
!       fodel2: same offshell (in case of complex or z-functions)
!       fodel3: -"-
!       cfdl4s: -"-
!       fodel4: -"-
!       fodl3i: -"-
!       fod3ij: -"-
!       fodl4i: -"-
!       fidel3: ier of del3 (is not included in D0)
!       fidel4: ier of del4 (is not included in E0)
!       fidl3i: ier of dl3i (is not included in E0)
!       fid3ij: ier of dl3ij (is not included in F0)
!       fidl4i: ier of dl4i (is not included in F0)
!
Real(dp), Private :: xloss = 0.125_dp, precx,precc,xalogm,xclogm,xalog2 &
 & ,xclog2,bf(20), xninv(30),xn2inv(30),xinfac(30), &
 &  fpij2(3,3),fpij3(6,6),fdel2, xlg2
! Real(dp), Private ::  reqprc=1.e-8_dp,fpij4(10,10),fpij5(15,15), &
! &  fpij6(21,21),fdel3,fdel4s,fdel4,fdl3i(5), &
! &  fdl4si(5),fdl3ij(6,6),fd4sij(6,6),fdl4i(6),fodel2, &
! &  fodel3,fodel4,fodl3i(5),fod3ij(6,6),fodl4i(6)
! Integer, Private :: fidel3,fidel4,fidl3i(5),fid3ij(6,6),fidl4i(6)

!
!	c[zero1]:0,1 complex
!	c2ipi:	2*i*pi
!	cipi2:	i*pi**2
!	cfp..:	complex version of fp..., only defined in ff[cz]*
!	cmipj:	(internal only) mi^2 - pj^2 in C0
!	c2sisj:	(internal only) 2*si.sj in D0
!	cfdl4s:	del4s in complex case (D0)
!	ca1:	(internal only) complex A1
!	csdl2p: (internal only) complex transformed sqrt(del2)
!
Complex(dp), Private :: czero,cone,c2ipi,cipi2,cmipj(3,3)
! Complex(dp), Private :: ,ca1,cfdl4s,c2sisj(4,4), cfpij6(21,21), cfpij2(3,3) &
!  & , cfpij3(6,6), cfpij4(10,10),cfpij5(15,15)

! error control
Integer, Private, Parameter :: N_err = 50
Logical, Private :: PrintToScreen = .False.
Character(len=80) :: St_Error(N_err)
Integer, Private :: ErrorOccur(N_err) = 0
!
!	        if .TRUE. then                        default (ffinit)
!	l4also: in C0 (and higher), also consider the algorithm with 16
!	        dilogs                                .TRUE.
!	ldc3c4: in D0 (and higher), also consider possible cancellations
!	        between the C0s                      .TRUE.
!	lmem:   before computing the C0 and higher, first check whether
!	        it has already been done recently     .FALSE.
!	ldot:   leave the dotproducts and some determinants in common
!	                                              .FALSE.
!	onshel: (in ffz?0 only): use onshell momenta  .TRUE.
!	lsmug:  internal use
!	  a flag to indicate the validity of d  Ifferences smuggled to the
!    	  IR routines in the C0 (ff internal only)
!	lnasty: internal use
!
Logical, Private :: l4also = .False., ldot = .False. ,lsmug = .False.
! Logical, Private ::  , lnasty, onshel = .True. ,lmem = .False., ldc3c4 = .False.
!
!       nwidth: number of widths within which the complex mass is used
!         in some schemes, for onshel=.FALSE.,
!         when |p^2-Re(m^2)| < nwidth*|Im(m^2)| special action is taken
!       nschem: scheme to handle the complex mass
!       nschem = 1:   Do not use the complex mass at all
!                2: only use the complex mass in linearly divergent terms
!                3: also use the complex mass in divergent logs UNDEFINED
!                4: use the complex mass in the C0   If there are
!                  divergent logs
!                5: include the almost-divergent threshold terms from
!                  (m,m,0) vertices
!                6: include the (s-m^2)*log(s-m^2) threshold terms from
!                  (m1+m2),m1,m2) vertices
!                7: full complex computation
!       idot:   internal flags to signal that some of the dotproducts
!               are input: 0: none; 1: external pi.pj, 2: external +
!               kinematical determinant, 3: all dotproducts + kindet.
!
Integer :: nwidth=5,nschem=7,idot=0
! check initialization
Logical, Private :: initialized = .False.
! coefficients for series
Real(dp), Private :: B0serie1(60), B0serie2(60), B0serie3(60), DB0serie1(60) &
   & , DB0serie2(60), DB0serie3(60)
! intermediate results for B and C functions
Complex(dp), Private, Save :: ca0i(2), cb0, cb1
!
!	nevent:	number in integration loop (to be updated by user)
!	ner:	can be used to signal numerical problems (see ffrcvr)
!	id:	identifier of scalar function (to be set by user)
!	idsub:	internal identifier to pinpoint errors
!	inx:	in D0: p(inx(i,j)) = isgn(i,j)*(s(i)-s(j))
!	inx5:	in E0: p(inx5(i,j)) = isgn5(i,j)*(s(i)-s(j))
!	inx6:	in F0: p(inx6(i,j)) = isgn6(i,j)*(s(i)-s(j))
!	isgn:	see inx
!	isgn5:	see inx5
!	isgn6:	see inx6
!	iold:	rotation matrix for 4point function
!	isgrot:	signs to iold
!	isgn34:	+1 or -1: which root to choose in the transformation (D0)
!	isgnal:	+1 or -1: which root to choose in the alpha-trick (C0)
!	irota3:	save the number of positions the C0 configuration has been 
!		rotated over
!	irota4:	same for the D0
!	irota5:	same for the E0
!	irota6:	same for the F0
!
Integer, Private :: id,idsub,inx(4,4),isgn(4,4),inx5(5,5)  &
        & , isgn5(5,5),inx6(6,6),isgn6(6,6),isgnal=1  &
        & , irota3 !,nevent = 0,ner,isgn34=1,irota4,irota5,irota6
!Integer, Private ::  iold(13,12)  = Reshape(  Source = (/ &
!  & 1,2,3,4, 5,6,7,8,9,10, 11,12,13, 4,1,2,3, 8,5,6,7,10,9, 11,13,12,  &
!  & 3,4,1,2, 7,8,5,6,9,10, 11,12,13, 2,3,4,1, 6,7,8,5,10,9, 11,13,12,  &
!  & 4,2,3,1, 10,6,9,8,7,5, 12,11,13, 1,3,2,4, 9,6,10,8,5,7, 12,11,13,  &
!  & 1,2,4,3, 5,10,7,9,8,6, 13,12,11, 1,4,3,2, 8,7,6,5,9,10, 11,13,12,  &
!  & 3,4,2,1, 7,10,5,9,6,8, 13,12,11, 2,3,1,4, 6,9,8,10,5,7, 12,13,11,  &
!  & 4,2,1,3, 10,5,9,7,8,6, 13,11,12, 1,3,4,2, 9,7,10,5,8,6, 13,11,12/) &
!  & , Shape = (/13, 12/) )
!Integer, Private :: isgrot(10,12) =  Reshape(  Source = (/ &
!  &  +1,+1,+1,+1, +1,+1,+1,+1, +1,+1, +1,+1,+1,+1, +1,+1,+1,+1, -1,+1, &
!  &  +1,+1,+1,+1, +1,+1,+1,+1, -1,-1, +1,+1,+1,+1, +1,+1,+1,+1, +1,-1, &
!  &  +1,+1,+1,+1, -1,+1,+1,-1, +1,-1, +1,+1,+1,+1, -1,-1,+1,+1, -1,+1, &
!  &  +1,+1,+1,+1, +1,+1,-1,+1, +1,+1, +1,+1,+1,+1, -1,-1,-1,-1, +1,-1, &
!  &  +1,+1,+1,+1, -1,+1,+1,+1, -1,-1, +1,+1,+1,+1, +1,+1,+1,-1, +1,-1, &
!  &  +1,+1,+1,+1, -1,+1,+1,-1, -1,-1, +1,+1,+1,+1, -1,-1,+1,+1, -1,-1/) &
!  & , Shape = (/10, 12/) )
Integer :: idum93(2)
! constants

!  functions and subroutines
Private :: AbsC, FindBound, ffRot3, WriteLFerror

! now my definitions
Logical, Private :: l_look_up_cache = .False. ! .True.
Integer, Private, Parameter :: n_c_max=5000
Integer, Private :: num_ci = 0
Type(C_data), Private :: C_cache(n_c_max)
! private variables

! private variables for FeynFunctionX (X=A,B,C,D,E)
Real(dp), Private, Parameter :: sa(7) = (/ -1._dp / 3._dp, 1._dp / 4._dp  &
    & ,  -1._dp / 5._dp,  1._dp / 6._dp,  -1._dp / 7._dp,   1._dp / 8._dp &
    & ,  -1._dp / 9._dp /)                                                &
    & , sb(7) = (/  1._dp / 6._dp, -1._dp /12._dp,  1._dp /20._dp         &
    &            , -1._dp /30._dp,  1._dp /42._dp, -1._dp /56._dp         &
    &            ,  1._dp /72._dp  /)                                     &
    & , sd(7) = (/ -1._dp / 20._dp,  1._dp / 15._dp, -1._dp / 14._dp      &
    &            ,  1._dp / 14._dp, -5._dp / 72._dp,  1._dp / 15._dp      &
    &            , -7._dp /110._dp /)                                     &
    & , se(7) = (/  1._dp / 30._dp, -1._dp / 30._dp,  1._dp / 35._dp      &
    &            , -1._dp / 42._dp,  5._dp /252._dp, -1._dp / 60._dp      &
    &            ,  7._dp /495._dp /)

Contains 


 Complex(dp) Function A_one(x)
 !----------------------------------------------------------------------------
 ! calculates the spin 1 loop function for the radiative decay of a scalar
 ! to two photons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! written by Werner Porod, 27.12.2008
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Real(dp) :: x2

  If (x.Eq.1._dp) Then
   A_one = -5._dp - 0.75_dp * Pi2 
  Else If (x.Lt.1._dp) Then
   x2 = x**2
   A_one = - ( 2._dp * x2 + 3._dp * x &
         &   + 3._dp * (2._dp * x-1._dp) * Asin(Sqrt(x))**2 ) / x2
  Else
   x2 = x**2
   A_one = - 0.25_dp * ( Log( (1._dp + Sqrt(1._dp - 1._dp/x))      &
       &                  / (1._dp - Sqrt(1._dp - 1._dp/x)) )    &
       &             - Cmplx(0._dp, Pi, dp) )**2
   A_one = - ( 2._dp * x2 + 3._dp * x &
         &   + 3._dp * (2._dp * x-1._dp) * A_one ) / x2
  End If

 End Function A_one


 Complex(dp) Function A_onehalf(x)
 !----------------------------------------------------------------------------
 ! calculates the spin 1/2 loop function for the radiative decay of a scalar
 ! to two photons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! written by Werner Porod, 27.12.2008
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  If (x.Eq.1._dp) Then
   A_onehalf = 2._dp 
  Else If (x.Lt.1._dp) Then
   A_onehalf = 2._dp * ( x + (x-1._dp) * Asin(Sqrt(x))**2 ) / x**2
  Else
   A_onehalf = - 0.25_dp * ( Log( (1._dp + Sqrt(1._dp - 1._dp/x))      &
       &                  / (1._dp - Sqrt(1._dp - 1._dp/x)) )    &
       &             - Cmplx(0._dp, Pi, dp) )**2
   A_onehalf = 2._dp * ( x + (x-1._dp) * A_onehalf ) / x**2
  End If

 End Function A_onehalf


 Complex(dp) Function AP_onehalf(x)
 !----------------------------------------------------------------------------
 ! calculates the spin 1/2 loop function for the radiative decay of a 
 ! pseudoscalar to two photons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! based on A_onehalf by Werner Porod
 ! written by Florian Staub, 18.06.2010
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

!  If (x.Eq.1._dp) Then
!   A_onehalf = 2._dp 
!  Else If (x.Lt.1._dp) Then
  If (x.LE.1._dp) Then
   AP_onehalf = Asin(Sqrt(x))**2  / x
  Else
   AP_onehalf = - 0.25_dp * ( Log( (1._dp + Sqrt(1._dp - 1._dp/x))      &
       &                  / (1._dp - Sqrt(1._dp - 1._dp/x)) )    &
       &             - Cmplx(0._dp, Pi, dp) )**2
   AP_onehalf = AP_onehalf  / x
  End If

 End Function AP_onehalf

 Complex(dp) Function A_zero(x)
 !----------------------------------------------------------------------------
 ! calculates the spin 0 loop function for the radiative decay of a scalar
 ! to two photons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! written by Werner Porod, 27.12.2008
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  If (x.Eq.1._dp) Then
   A_zero = -1._dp + 0.25_dp * Pi2 
  Else If (x.Lt.1._dp) Then
   A_zero = ( - x + Asin(Sqrt(x))**2 ) / x**2
  Else
   A_zero = - 0.25_dp * ( Log( (1._dp + Sqrt(1._dp - 1._dp/x))      &
       &                  / (1._dp - Sqrt(1._dp - 1._dp/x)) )    &
       &             - Cmplx(0._dp, Pi, dp) )**2
   A_zero = ( - x + A_zero ) / x**2
  End If

 End Function A_zero


 Complex(dp) Function A0(m)
 Implicit None
  Real(dp), Intent(in) :: m

  Real(dp) :: xmu, xlogm

  If ( mudim2 .Ne. 0 ) Then
   xmu = m/mudim2
  Else
   xmu = m
  End If
  If ( xmu .Gt. xalogm ) Then
   xlogm = Log(xmu)
  Else
   xlogm = 0
   If ( xmu .Ne. 0 ) Call WriteLFerror(1)
  End If 

  A0 = -(m*(xlogm - 1 - divergence))

 End Function A0


 Real(dp) Function AbsC(z)
 Implicit None
  Complex(dp), Intent(in) :: z

  AbsC = Abs( Real(z,dp) ) + Abs( Aimag(z) )

 End Function AbsC


 Complex(dp) Function B00(xp, xm1, xm2)
 Implicit None
  Real(dp), Intent(in) :: xp, xm1, xm2
  Complex(dp) :: cbii(2)

  cbii = Bii(xp, xm1, xm2)
  B00 = cbii(2)
 End Function B00


 Complex(dp) Function B0(xp, xma, xmb)
 Implicit None
  Real(dp), Intent(in) :: xp, xma, xmb

  Real(dp) :: dmamb, dmap, dmbp, xm, dmp, s, s1, sumI, y, xx, slam, xlam &
     & , s2, s2a, ax, x, s1a, s1b, s2b, s1p, s2p, xtel, xnoe, alpha, xm1 &
     & , alph1, dm1m2, dm1p, xlogmm, xm2, dm2p
  Complex(dp) :: cx, cs2a, cs2b, cs2p
  Integer :: i_max, i1, JSIGN 

  If (.Not.initialized) Then
   Write(ErrCan,*) "Sorry, you forgot to call InitialzeLoopFunctions"
   Call TerminateProgram
  End If

  dmamb = xma - xmb
  dmap = xma - xp
  dmbp = xmb - xp

  !-------------
  ! dot product
  !-------------
  If (ldot) Then
   fpij2(1,1) = xma
   fpij2(2,2) = xmb
   fpij2(3,3) = xp
   If ( Abs(dmap) .Lt. Abs(dmbp) ) Then
     fpij2(1,2) = (dmap + xmb)/2
   Else
     fpij2(1,2) = (dmbp + xma)/2
   End If
   fpij2(2,1) = fpij2(1,2)
   If ( Abs(dmamb) .Lt. Abs(dmbp) ) Then
     fpij2(1,3) = (-dmamb - xp)/2
   Else
     fpij2(1,3) = (dmbp - xma)/2
   End If
   fpij2(3,1) = fpij2(1,3)
   If ( Abs(dmamb) .Lt. Abs(dmap) ) Then
     fpij2(2,3) = (-dmamb + xp)/2
   Else
     fpij2(2,3) = (-dmap + xmb)/2
   End If
   fpij2(3,2) = fpij2(2,3)
  End If ! ldot

  !---------------
  ! various cases
  !---------------
  If ( (xma.Eq.0._dp) .And. (xmb.Eq.0._dp) ) Then      ! both masses are 0
   If ( xp .Lt. -xalogm ) Then
    B0 = Log(-xp) - 2._dp
   Else If ( xp .Gt. xalogm ) Then
    B0 = Cmplx( Log(xp) - 2._dp, -pi, dp )
   Else
    B0 = Huge(1._dp)
    Call WriteLFerror(2)
   End If

  Else If ( (xma.Eq.0._dp) .Or. (xmb.Eq.0._dp) ) Then  ! one mass is 0
   If (xma.Eq.0._dp) Then
    xm = xmb
    dmp = dmbp
   Else
    xm = xma
    dmp = dmap
   End If

   If ( xp .Eq. 0 ) Then        ! zero momentum
    B0 = - 1._dp
   Else If ( dmp.Eq.0 ) Then     ! p^2 = m^2
    B0 = - 2._dp
   Else                         ! p^2 != m^2
    s1 = xp/xm
    If ( Abs(s1) .Lt. xloss ) Then
     s = Log1minusX(s1)
    Else
      s = Log(Abs(dmp/xm))
    End If
    s = -s*dmp/xp
    B0 = s - 2._dp
    If ( xp .Gt. xm ) B0 = B0 - Cmplx(0.0_dp,-(dmp/xp)*pi, dp)
   End If

  Else If (dmamb.Eq.0._dp) Then                        ! both masses equal
   xm = xma
   dmp = dmap

   If ( Abs(xp).Lt.8._dp*xloss*xm) Then ! talor series, result goes like 1/p
    s = - xp / xm
    i_max = FindBound(s, 1, B0serie1)
    sumI = s * B0serie1(i_max)
    Do i1=i_max-1,1,-1
     sumI = s * (B0serie1(i1) + sumI)
    End Do
    B0 = Cmplx(sumI,0._dp,dp)
   Else

    xlam = Kappa2p(-xp,-xm,-xm,dmp,dmp,0._dp)
    If ( xlam .Ge. 0 ) Then
     slam = Sqrt(xlam)
     s2a = dmp + xm
     s2 = s2a + slam
     If ( Abs(s2) .Gt. xloss*slam ) Then
      jsign = 1
     Else
      s2 = s2a - slam
      jsign = -1
     End If
     ax = Abs(s2/(2*xm))
     If ( ax .Lt. xalogm ) Then
      s = 0._dp
     Else If ( ((ax-1._dp).Lt.(0.1_dp)).And.(s2.Gt.0._dp) ) Then
      s2 = (xp - slam)
      s = -slam/xp*Log1minusX(s2/(2._dp*xm))
     Else
      s = -slam/xp*Log(ax)
      If ( jsign .Eq. -1 ) s = -s
     End If
     If ( xp .Gt. 2*xm ) Then ! in this case these is an imaginary part
      y = -pi*slam/xp
     Else
      y = 0._dp
     End If
    Else  ! the root is complex (p**2 between 0 and (2*m1)^2)
     slam = Sqrt(-xlam)
     s = 2*slam/xp*Atan2(xp,slam)
     y = 0
    End If
    xx = s - 2
    B0 = Cmplx(xx,y,dp)
   End If

  Else                                                 ! general case

   If ( xma .Gt. xmb ) Then
    xm2 = xma
    xm1 = xmb
    dm1m2 = -dmamb
    dm1p = dmbp
    dm2p = dmap
   Else
    xm1 = xma
    xm2 = xmb
    dm1m2 = dmamb
    dm1p = dmap
    dm2p = dmbp
   Endif

   x = xm2/xm1

   If ( 1._dp .Lt. xalogm*x ) Then ! numerical problem
    Call WriteLFerror(3)
    xlogmm = 0._dp
   Elseif ( Abs(x-1) .Lt. xloss ) Then
    xlogmm = Log1minusX(dm1m2/xm1)
   Else
    xlogmm = Log(x)
   Endif

   If (xp.Eq.0) Then
    s2 = ((xm2+xm1) / dm1m2)*xlogmm
    s = - s2 - 2._dp
    sumI = s
    If ( Abs(s) .Lt. xloss*2 ) Then
     !  Taylor expansions: choose which one
     x = dm1m2/xm1
     ax = Abs(x)
     If ( (ax.Lt.0.15_dp) .Or. (precx.Gt.1.e-8_dp) .And. (ax.Lt.0.3_dp) ) Then
      i_max = FindBound(x, 2, B0serie2)
      sumI = x * B0serie2(i_max)
      Do i1=i_max-1,1,-1
       sumI = x * (B0serie2(i1) + sumI)
      End Do
      sumI = x*sumI

     Else   
      y = 2._dp*x/(2._dp-x)
      i_max = FindBound(y, 1, B0serie3)
      sumI = 0._dp
      Do i1=i_max,5,-1
       sumI = y * (B0serie3(i1) + sumI)
      End Do
      sumI = (1-x)*y**4*(B0serie3(4)+sumI)
      sumI = x*y**2*(1+y)/12._dp - sumI
      sumI = - 2._dp*Log1minusX(sumI)/y
     Endif
    End If
    B0 = 0.5_dp * sumI

   Else

    xlam = Kappa2p(-xp,-xm2,-xm1,dm2p,dm1p,dm1m2)
    If ( xlam .Gt. 0 ) Then
    ! reel cases k^2 < -(m2+m1)^2 or k^2 > -(m2-m1)^2:
    !       first try:
     slam = Sqrt(xlam)
     s2a = dm2p + xm1
     s2 = s2a + slam
     If ( Abs(s2) .Gt. xloss*slam ) Then
      jsign = 1
     Else
      s2 = s2a - slam
      jsign = -1
     Endif
     s2 = s2**2/(4._dp*xm1*xm2)
     If ( Abs(s2) .Lt. xalogm ) Then
      Call WriteLFerror(4)
      s2 = 0
     Elseif ( Abs(s2-1) .Lt. xloss ) Then
      If ( jsign.Eq.1 ) Then
       s2 = -slam*(s2a+slam)/(2*xm1*xm2)
       s2 = -slam/(2*xp)*Log1minusX(s2)
      Else
       s2 = +slam*(s2a-slam)/(2*xm1*xm2)
       s2 = +slam/(2*xp)*Log1minusX(s2)
      Endif
     Else
      s2 = -slam/(2*xp)*Log(s2)
      If ( jsign .Eq. -1 ) s2 = -s2
     Endif
     s1 = -dm1m2*xlogmm/(2._dp*xp)
     xx = s1+s2-2._dp

     If ( Abs(xx) .Lt. xloss*Max(Abs(s1),Abs(s2)) ) Then
     !  this is unacceptable, try a better solution
      s1a = dm1m2 + slam
      If ( Abs(s1a) .Gt. xloss*slam ) Then
       s1 = -s1a/(2*xp)
      Else
       s1 = ( -0.5_dp*xp + xm1 + xm2 ) / ( slam - dm1m2 )
      Endif
      s1 = s1*xlogmm
      If ( Abs(xp) .Lt. xm2 ) Then
       s2a = xp - dm1m2
      Else
       s2a = xm2 - dm1p
      Endif
      s2 = s2a - slam
      If ( Abs(s2) .Gt. xloss*slam ) Then
       s2 = s2 / (2._dp*xm2)
      Else
       s2 = (2*xp) / (s2a+slam)
      Endif
      If ( Abs(s2) .Lt. 0.1_dp ) Then
       s2 = Log1minusX(s2)
      Else
       s2a = Abs(1-s2)
       s2 = Log(s2a)
      Endif
      s2 = -(slam/xp)*s2
      xx = s1 + s2 - 2._dp

      If ( Abs(xx) .Lt. xloss**2*Max(Abs(s1),Abs(s2)) ) Then
      !        third try:
      ! (we accept two times xloss because that's the same as in this try)
      ! A Taylor expansion might work.  We expand
      ! inside the logs. Only do the necessary work.

       alpha = slam/(slam-dm1m2)
       alph1 = -dm1m2/(slam-dm1m2)

       s1p = s1 - 2*alph1
       If ( Abs(s1p) .Lt. xloss*Abs(s1) ) Then
        xnoe = -xp + 2._dp*xm1 + 2._dp*xm2
        x = 4._dp*dm1m2/xnoe
        i_max = FindBound(x, 1, B0serie3)
        If (i_max.Gt.17) i_max=17
        s1a = 0._dp 
        Do i1=i_max,4,-1
         s1a = x * (B0serie3(i1) + s1a)
        End Do
        s1a = x**3 *(B0serie3(3) + s1a) *xm2/xm1
        s1b = dm1m2*(4*dm1m2**2 - xp*(4*xm1-xp))/ (xm1*xnoe**2)
        s1p = s1b - s1a
        s1p = xnoe*Log1minusX(s1p)/(slam - dm1m2)/2
       Endif

       s2p = s2 - 2*alpha
       If ( Abs(s2p) .Lt. xloss*Abs(s2) ) Then
        xnoe = slam - dm1m2
        x = 2*xp/xnoe
        i_max = FindBound(x, 1, B0serie3)
!        if (i_max.lt.Size(B0serie3)) then
        If (i_max.Lt.19) Then
         s2a = 0._dp 
         Do i1=i_max,5,-1
          s2a = x * (B0serie3(i1) + s2a)
         End Do
         s2a = x**4*(B0serie3(4)+s2a)*(1._dp-2._dp*xp/(xnoe+xp))
         s2b = -2*xp**3 *(-2*xp - xnoe)/(3*(xnoe+xp)* xnoe**3)
         s2p = s2b - s2a
         s2p = -slam/xp*Log1minusX(s2p)
        Endif
       Endif

       xx = s1p + s2p
      Endif
     Endif ! second try

     If ( xp .Gt. xm1+xm2 ) Then
      y = -pi*slam/xp
     Else
     y = 0
     Endif

    Else ! the root is complex (k^2 between -(m1+m2)^2 and -(m2-m1)^2)

     slam = Sqrt(-xlam)
     xnoe = dm2p + xm1
     s1 = -(dm1m2/(2*xp))*xlogmm
     s2 = (slam/xp)*Atan2(slam,xnoe)
     xx = s1 + s2 - 2

! encountered special case where xnoe.eq.0
! in this case I accept the numerical worse case
     If ((Abs(xx) .Lt. xloss**2*Max(Abs(s1),Abs(s2))).and.(xnoe.ne.0._dp)) Then
!     If ( Abs(xx) .Lt. xloss**2*Max(Abs(s1),Abs(s2)) ) Then
     !      second try:
     !  Again two times xloss as we'll accept that in the next step as well.

      xtel = dm1m2**2 - xp**2
      alpha = -xlam/(2._dp*xp*xnoe)
      alph1 = xtel/(2._dp*xp*xnoe)

      s1p = s1 - 2*alph1
      If ( Abs(s1p) .Lt. xloss*Abs(s1) ) Then
       x = 2._dp*xtel/(dm1m2*xnoe)
       i_max = FindBound(x, 1, B0serie3)
!       if (i_max.lt.Size(B0serie3)) then
       If (i_max.Lt.18) Then
        s1a = 0._dp 
        Do i1=i_max,4,-1
         s1a = x * (B0serie3(i1) + s1a)
        End Do
        s1a = x**3 *(B0serie3(3) + s1a) *xm2/xm1
        s1b = (dm1m2**3*(dm1m2**2-2*xp*xm1)                       &
            & + xp**2*(4*dm1m2*xm1**2-dm1m2**2*(dm1m2+2*xm1))     &
            & - 2*xm2* xp**3*(dm1m2+xp))/(xm1*dm1m2**2*xnoe**2)
        s1p = s1b - s1a
        s1p = -dm1m2*Log1minusX(s1p)/(2._dp*xp)
       Endif
      Endif

      s2p = s2 - 2*alpha
      If ( Abs(s2p) .Lt. xloss*Abs(s2) ) Then
       cx = Cmplx(0._dp,-slam/xnoe,dp)
       ax = absc(cx)
       i_max = FindBound(ax, 1, B0serie3)
!       if (i_max.lt.Size(B0serie3)) then
       If (i_max.Lt.18) Then
        cs2a = 0._dp 
        Do i1=i_max,4,-1
         cs2a = cx * (B0serie3(i1) + cs2a)
        End Do
        cs2a = cx**3*(B0serie3(3)+cs2a)*Cmplx(xnoe,slam,dp)
        cs2b = Cmplx(xnoe-xlam/xnoe/2, -slam**3/xnoe**2/2,dp)
        cs2p = cs2b + cs2a
        s2p = slam*Atan2(Aimag(cs2p),Real(cs2p,dp))/xp
       Endif
      Endif

      xx = s1p + s2p

     Endif
     y = 0
    Endif
    B0 = Cmplx(xx,y,dp)
   End If

  End If  ! various cases for the masses

  If ((xma.Eq.0._dp).And.(xmb.Eq.0._dp)) Then
    xm = 1._dp
  Else If (xma.Eq.0._dp) Then
    xm = xmb**2
  Elseif (xmb.Eq.0._dp) Then
   xm = xma**2
  Else
   xm = xma*xmb
  Endif
  If (mudim2.Ne.0._dp) xm = xm/mudim2**2

  If (Abs(xm).Gt.xalogm) Then
   B0 = divergence - 0.5_dp * Log(xm) - B0
  Else
   Call WriteLFerror(5)
   B0 = divergence - B0
  Endif
  
 End Function B0


 Complex(dp) Function B11(xp, xm1, xm2)
 Implicit None
  Real(dp), Intent(in) :: xp, xm1, xm2
  Complex(dp) :: cbii(2)

  cbii = Bii(xp, xm1, xm2)
  B11 = cbii(1)
 End Function B11


 Complex(dp) Function B1(xp, xm1, xm2)
 Implicit None
  Real(dp), Intent(in) :: xp, xm1, xm2

  Integer :: i_max, i1
  Logical :: lneg
  Real(dp) :: dm1m2, xmax, s, s1, slam, h, xma, xmb, x, xlogm, small, dmbma &
     & , xlam, ts2Dp, xmxp, xlo3
  Complex(dp) :: cs(5), csom

  ca0i(1) = A0(xm1)
  ca0i(2) = A0(xm2)
  ldot = .True.
  cb0 = B0(xp, xm1, xm2)
  ldot = .False.

  dm1m2 = xm1 - xm2

  If ( xp .Ne. 0 ) Then   ! p^2 != 0
   If ( dm1m2 .Ne. 0 ) Then
    cs(1) = -ca0i(2)
    cs(2) = +ca0i(1)
   Else
    cs(1) = 0
    cs(2) = 0
   Endif
   cs(3) = Real(2._dp*fpij2(1,3),dp) * cb0
   B1 = cs(1) + cs(2) + cs(3)
   xmax = Max(absc(cs(2)),absc(cs(3)))

   If ( absc(B1) .Ge. xloss*xmax ) Then
    B1 = B1 *(1._dp/Real(2*xp,dp))

   Else

    If ( Abs(dm1m2) .Le. xloss*xm1 ) Then  ! almost equal masses
     cs(2) = Real(dm1m2/xm1,dp)*cs(2)
     cs(1) = -xm2*Log1minusX(-dm1m2/xm2)
     B1 = cs(1) + cs(2) + cs(3)
     xmax = Max(absc(cs(2)),absc(cs(3)))
     If ( absc(B1) .Ge. xloss*xmax ) Then
      B1 = B1 *(1._dp/Real(2*xp,dp))
      Return
     End If
    Endif

    If ( xloss**2*Max(xm1,xm2) .Gt. Abs(xp) ) Then ! p^2 -> 0

     If ( xm2.Gt.xm1 ) Then
      xma = xm1
      xmb = xm2
      ts2Dp = +2*fpij2(2,3)
      lneg = .False.
     Else
      xma = xm2
      xmb = xm1
      ts2Dp = -2*fpij2(1,3)
      lneg = .True.
     Endif

     dmbma = Abs(dm1m2)
     If ( xma.Eq.0 ) Then
      xlogm = 1
     Elseif ( dmbma .Gt. xloss*xmb ) Then
      xlogm = Log(xmb/xma)
     Else
      xlogm = Log1minusX(-dmbma/xma)
     Endif
     xlam =  (dmbma-xp)**2 - 4*xma*xp

     If ( xlam.Gt.0 ) Then ! real roots
      slam = Sqrt(xlam)
      small = xp*(-2*(xma+xmb) + xp)/(slam+dmbma)
      h = slam+2*fpij2(1,2)
      cs(1) = xlogm*xma*(4*xmb*(small-xp) + (small-xp)**2)/(2*(slam+dmbma)*h)
      x = xp/slam
      i_max = FindBound(x, 1, B0serie3)
      s = 0._dp
      Do i1=i_max,3,-1
       s = x * (B0serie3(i1) + s)
      End Do
      s = x**2*(0.5_dp + s)
      h = ts2Dp + slam
      s1 = 2*xp/h*(s + x)
      h = -4*xp**2*xmb/(slam*h**2) - s + s1
      If ( Abs(h) .Lt. .1 ) Then
       cs(2) = dmbma*slam/xp*Log1MinusX(h)
      Else
       B1 = B1 *(1._dp/Real(2*xp,dp))
       Return
      Endif
      If ( lneg ) Then
       cs(1) = -cs(1)
       cs(2) = -cs(2)
      Endif
      cs(3) = -Real(xp,dp)*cb0
      B1 = cs(1) + cs(2) + cs(3)
      B1 = B1 *(1._dp/Real(2*xp,dp))

     Else ! imaginary roots, not implemented correctly
      B1 = B1 *(1._dp/Real(2*xp,dp))
     Endif

    Else
     B1 = B1 *(1._dp/Real(2*xp,dp))
     Return
    Endif
   Endif

  Elseif ( dm1m2 .Ne. 0 ) Then ! p^2=0, m1 != m2
   cs(1) = +Real(xm2/(2*dm1m2**2),dp)*(ca0i(2)+xm2/2)
   cs(2) = -Real(xm1/(2*dm1m2**2),dp)*(ca0i(1)+xm1/2)
   cs(3) = +ca0i(2)*(1/dm1m2)
   B1 = cs(1) + cs(2) + cs(3)
   xmax = Max(absc(cs(1)),absc(cs(2)),absc(cs(3)))
   If ( absc(B1).Ge.xloss**2*xmax ) Return ! give up

   If ((xm1.Eq.0._dp).Or.(xm2.Eq.0._dp)) Return ! give up

   If ( Abs(dm1m2).Lt.xloss*xm1 ) Then
    xlogm = Log1minusX(dm1m2/xm1)
   Else
    xlogm = Log(xm2/xm1)
   Endif
   cs(1) = -(xm1/dm1m2)/2
   cs(2) = -xlogm/2*(xm1/dm1m2)**2
   cs(3) = 0.25_dp - ca0i(1)*Real(1._dp/(2._dp*xm1),dp)
   cs(4) = xlogm/2
   csom = cs(1) + cs(2) + cs(3) + cs(4)
   xmxp = Max(absc(cs(2)),absc(cs(3)),absc(cs(4)))
   If ( xmxp.Lt.xmax ) Then
    xmax = xmxp
    B1 = csom
    If ( absc(B1).Gt.xloss**2*xmax ) Return ! give up
   Endif

   xlo3 = Log1minusXpXn(dm1m2/xm1,2)
   cs(1) = -(dm1m2/xm1)**2/4
   cs(2) = -(dm1m2/xm1)/2
   cs(3) = -xlo3/(dm1m2/xm1)**2/2
   cs(4) = xlo3/2
   cs(5) = 0.5_dp - ca0i(1)*Real(1/(2*xm1),dp)
   csom = cs(1) + cs(2) + cs(3) + cs(4) + cs(5)
   xmxp = Max(absc(cs(2)),absc(cs(3)),absc(cs(4)),absc(cs(5)))
   If ( xmxp.Lt.xmax ) Then
    xmax = xmxp
    B1 = csom
    If ( absc(B1).Gt.xloss**2*xmax ) Return ! give up
   Endif

  Else ! p^2=0, m1 == m2
   B1 = -cb0/2
  Endif

 End Function B1


 Complex(dp) Function B22(p2,m12,m22)
 !-----------------------------------------------------------------------
 ! calculates the function \hat B_22 as defined in J. Bagger at al, Nucl.Phys.B
 ! written by Werner Porod, 4.8.1999
 ! 18.05.2001: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: p2, m12, m22

  B22 = B00(p2,m12,m22) - 0.25_dp * ( A0(m12) + A0(m22) )

 End Function B22


 Function Bii(xp, xm1, xm2) Result(cb2i)
  Implicit None
  Real(dp), Intent(in) :: xp, xm1, xm2
  Complex(dp) :: cb2i(2)

  Logical :: llogmm
  Integer :: i1
  Real(dp) :: dm1m2, xmax, slam, xlo3, xmxp, bet, xnoe, xnoe2, xlam &
     & , xlogmm
  Complex(dp) :: cs(16), csom, clo3

  cb1 = B1(xp, xm1, xm2)

  dm1m2= xm1 - xm2

  If ( xp .Ne. 0) Then
   cs(1) = ca0i(2)
   cs(2) = xm1*cb0
   cs(3) = Real(2*fpij2(1,3),dp)*cb1
   cs(4) = (xm1+xm2)/2
   cs(5) = -xp/6
   cb2i(1) = cs(1) - cs(2) + 2*cs(3) - cs(4) - cs(5)
   cb2i(2) = cs(1) + 2*cs(2) - cs(3) + 2*cs(4) + 2*cs(5)
   xmax = Max(absc(cs(2)),absc(cs(3)),absc(cs(4)),absc(cs(5)))

   If ( absc(cb2i(1)) .Ge. xloss*xmax ) Then
    cb2i(1) = Real(1/(3*xp),dp) * cb2i(1)
    cb2i(2) = 1/6._dp * cb2i(2)
    Return
   End If

   If ( dm1m2.Eq.0 .And. xm1.Ne.0 ) Then
    If ( xp.Lt.0 ) Then
     slam = Sqrt(xp**2-4*xm1*xp)
     xlo3 = Log1minusXpXn((xp-slam)/(2*xm1),2)
     cs(1) = xp*(-1/3._dp + slam/(4*xm1))
     cs(2) = xp**2*(-slam/(4*xm1**2) - 3/(4*xm1))
     cs(3) = xp**3/(4*xm1**2)
     cs(4) = Real(xp/xm1,dp)*ca0i(1)
     cs(5) = xlo3/xp*(-xm1*slam)
     cs(6) = xlo3*slam
    Else
     slam = isgnal*Sqrt(-xp**2+4*xm1*xp)
     clo3 = Log1minusXpXn(Cmplx(xp/(2*xm1),-slam/(2*xm1),dp),2)
     cs(1) = xp * Cmplx(-1/3._dp,slam/(4*xm1),dp)
     cs(2) = Real(xp**2,dp)*Cmplx(-3/(4*xm1),-slam/(4*xm1**2),dp)
     cs(3) = xp**3/(4*xm1**2)
     cs(4) = Real(xp/xm1,dp)*ca0i(1)
     cs(5) = clo3*Cmplx(0._dp,-xm1*slam/xp,dp)
     cs(6) = clo3*Cmplx(0._dp,slam,dp)
    Endif
    csom = cs(1) + cs(2) + cs(3) + cs(4) + cs(5) + cs(6)
    xmxp = Max(absc(cs(2)),absc(cs(3)),absc(cs(4)),absc(cs(5)),absc(cs(6)))
   
 !    get rid of noise in the imaginary part  
    If ( xloss*Abs(Aimag(csom)).Lt.precc*Abs(Real(csom,dp)) ) & 
     &      csom = Cmplx(Real(csom,dp),0._dp,dp)
    If ( xmxp.Lt.xmax ) Then
     cb2i(1) = csom
     xmax = xmxp
    Endif
    If ( absc(cb2i(1)).Ge.xloss**2*xmax ) Then
     cb2i(1) = Real(1/(3*xp),dp) * cb2i(1)
     cb2i(2) = 1/6._dp * cb2i(2)
     Return
    End If
   Endif

   xlam =  4*(fpij2(1,3)**2 - xm1*xp)
   If ( (xm1.Eq.0._dp).Or.(xm2.Eq.0._dp) ) Then
    xlogmm = 0
   Elseif ( Abs(dm1m2).Lt.xloss*xm1 ) Then
    xlogmm = Log1minusX(dm1m2/xm1)
   Else
    xlogmm = Log(xm2/xm1)
   Endif
   If ((xlam.Gt.0) .And. (Abs(xp).Lt.xloss*xm2) .And. (xm1.Lt.xm2) ) Then
    slam = Sqrt(xlam)
    bet = 4*xm1*xp/(2*fpij2(1,3)+slam)
    cs(1) = Real(xp/xm2,dp)*ca0i(2)
    cs(2) = xlogmm*bet*(-2*xm1**2*xm2 - 2*xm1**3)                    &
        & /((-dm1m2+slam)*(2*fpij2(1,2)+slam)*(2*fpij2(1,3)+slam))
    cs(3) = xlogmm*(-4*xp*xm1**3)                                    &
        & /((-dm1m2+slam)*(2*fpij2(1,2)+slam)*(2*fpij2(1,3)+slam))
    xnoe = 1/(2*fpij2(2,3)+slam)
    xnoe2 = xnoe**2
    cs(4) = xnoe2*xm1*bet*(xp-4*xm2)
    cs(5) = xnoe2*xm1*2*xp*xm2
    cs(6) = xnoe2*xm1**2*bet
    cs(7) = xnoe2*xm1**2*4*xp
    cs(8) = xnoe2*bet*(xp*xm2+3*xm2**2)
    cs(9) = xnoe2*(-6*xp*xm2**2)
    cs(10)= xp*(7/6._dp - 2*xm1*slam*xnoe2 + 4*xm2*slam*xnoe2 - 2*slam*xnoe)
    cs(11)= xp**2*( -2*slam*xnoe2 )
    xlo3 = Log1minusXpXn(2*xp*xnoe,2)
    cs(12) = xlo3*dm1m2**2*slam/xp**2
    cs(13) = xlo3*(xm1 - 2*xm2)*slam/xp
    cs(14) = xlo3*slam
    csom = 0
    xmxp = 0
    Do i1=1,14
     csom = csom + cs(i1)
     xmxp = Max(xmxp,absc(cs(i1)))
    End Do
    If ( xmxp.Lt.xmax ) Then
     cb2i(1) = csom
     xmax = xmxp
    Endif
    If ( absc(cb2i(1)).Ge.xloss**2*xmax ) Then
     cb2i(1) = Real(1/(3*xp),dp) * cb2i(1)
     cb2i(2) = 1/6._dp * cb2i(2)
     Return
    End If
   Endif

   If ( (xlam.Gt.0) .And. (Abs(xp).Lt.xloss*xm1) .And. (xm2.Lt.xm1) ) Then
    slam = Sqrt(xlam)
    bet = 4*xm2*xp/(-2*fpij2(2,3)+slam)
    xnoe = 1/(-2*fpij2(1,3)+slam)
    xnoe2 = xnoe**2
    cs(1) = Real(xp/xm1,dp)*ca0i(1)
    cs(2) = -xlogmm*bet*(12*xp*xm1*xm2+6*xp*xm2**2                   &
        &                - 6*xp**2*xm2-2*xm1*xm2**2-2*xm2**3)        &
        &    /((dm1m2+slam)*(2*fpij2(1,2)+slam)*(-2*fpij2(2,3)+slam))
    cs(3) = -xlogmm*(-24*xp*xm1**2*xm2-4*xp*xm2**3+36*               &
     &     xp**2*xm1*xm2+12*xp**2*xm2**2-12*xp**3*xm2)                &
     &    /((dm1m2+slam)*(2*fpij2(1,2)+slam)*(-2*fpij2(2,3)+slam))
    cs(4) = xnoe2*xm2*bet*(xp-4*xm1)
    cs(5) = xnoe2*xm2*(-10*xp*xm1)
    cs(6) = xnoe2*xm2**2*bet
    cs(7) = xnoe2*xm2**2*4*xp
    cs(8) = xnoe2*bet*(xp*xm1+3*xm1**2)
    cs(9) = xnoe2*6*xp*xm1**2
    cs(10)= xp*(7/6._dp - 2*xm1*slam*xnoe2 + 4*xm2*slam*xnoe2 - 2*slam*xnoe)
    cs(11)= xp**2*( -2*slam*xnoe2 )
    xlo3 = Log1minusXpXn(2*xp*xnoe,2)
    cs(12) = xlo3*dm1m2**2*slam/xp**2
    cs(13) = xlo3*(xm1 - 2*xm2)*slam/xp
    cs(14) = xlo3*slam
    csom = 0
    xmxp = 0
    Do i1=1,14
     csom = csom + cs(i1)
     xmxp = Max(xmxp,absc(cs(i1)))
    End Do
    If ( xmxp.Lt.xmax ) Then
     cb2i(1) = csom
     xmax = xmxp
    Endif
    If ( absc(cb2i(1)).Ge.xloss**2*xmax ) Then
     cb2i(1) = Real(1/(3*xp),dp) * cb2i(1)
     cb2i(2) = 1/6._dp * cb2i(2)
     Return
    End If
   Endif

   cb2i(1) = Real(1/(3*xp),dp) * cb2i(1)
   cb2i(2) = 1/6._dp * cb2i(2)

  Elseif (dm1m2 .Ne. 0) Then ! p^2=0, m1!=m2:

   llogmm = .False.
! first b11
   cs(1) = Real(xm1**2/3/dm1m2**3,dp)*ca0i(1)
   cs(2) = Real((-xm1**2 + xm1*xm2 - xm2**2/3)/dm1m2**3,dp)*ca0i(2)
   cs(3) = (5*xm1**3/18 - xm1*xm2**2/2 + 2*xm2**3/9)/dm1m2**3
   cb2i(1) = cs(1)+cs(2)+cs(3)
   xmax = Max(absc(cs(2)),absc(cs(3)))
  ! implementing the case that one of the masses is 0
   If ( ( absc(cb2i(1)).Le.xloss**2*xmax ) .And.  &
      & (xm1.Ne.0._dp) .And. (xm2.Ne.0._dp) ) Then
    If ( Abs(dm1m2).Lt.xloss*xm1 ) Then ! xma ~ xmb
     xlogmm = Log1minusX(dm1m2/xm1)
    Else
     xlogmm = Log(xm2/xm1)
    Endif
    llogmm = .True.
    cs(1) = (xm1/dm1m2)/6
    cs(2) = (xm1/dm1m2)**2/3
    cs(3) = (xm1/dm1m2)**3*xlogmm/3
    cs(4) = -2/9._dp+ ca0i(1)*Real(1/(3*xm1),dp)
    cs(5) = -xlogmm/3
    csom = cs(1)+cs(2)+cs(3)+cs(4)+cs(5)
    xmxp = Max(absc(cs(2)),absc(cs(3)),absc(cs(4)),absc(cs(5)))
    If ( xmxp.Lt.xmax ) Then
     xmax = xmxp
     cb2i(1) = csom
    Endif

    !    and last try
    If (.Not.((xmxp.Lt.xmax).And.(absc(cb2i(1)).Gt.xloss**2*xmax))) Then
     xlo3 = Log1minusXpXn(dm1m2/xm1,2)
     cs(1) = (dm1m2/xm1)**2/6
     cs(2) = (dm1m2/xm1)/3
     cs(3) = xlo3/(3*(dm1m2/xm1)**3)
    ! same cs(4) = -2/DBLE(9) + ca0i(1)*DBLE(1/(3*xm1))
     cs(5) = -xlo3/3
     csom = cs(1)+cs(2)+cs(3)+cs(4)+cs(5)
     xmxp = Max(absc(cs(2)),absc(cs(3)),absc(cs(4)),absc(cs(5)))
     If ( xmxp.Lt.xmax ) Then
      xmax = xmxp
      cb2i(1) = csom
     End If
    Endif
   Endif

   ! now B00 

   cs(1) = Real(xm1/(4*dm1m2),dp)*ca0i(1)
   cs(2) = -Real(xm2/(4*dm1m2),dp)*ca0i(2)
   cs(3) = (xm1+xm2)/8
   cb2i(2) = cs(1) + cs(2) + cs(3)
   xmax = Max(absc(cs(2)),absc(cs(3)))
   If ( ( absc(cb2i(2)).Le.xloss*xmax ) .And. (xm2.Ne.0._dp) &
     & .And. (xm1.Ne.0._dp) ) Then

    If ( .Not.llogmm ) Then
     If ( Abs(dm1m2).Lt.xloss*xm1 ) Then
      xlogmm = Log1minusX(dm1m2/xm1)
     Else
      xlogmm = Log(xm2/xm1)
     Endif
    Endif
    cs(1) = dm1m2*( -1/8._dp - ca0i(1)*Real(1/(4*xm1),dp) )
    cs(2) = dm1m2*xlogmm/4
    cs(3) = xm1*(xm1/dm1m2)/4*xlogmm
    cs(4) = xm1*( 0.25_dp + ca0i(1)*Real(1/(2*xm1),dp) )
    cs(5) = -xm1*xlogmm/2
    csom = cs(1) + cs(2) + cs(3) + cs(4) + cs(5)
    xmxp = Max(absc(cs(2)),absc(cs(3)),absc(cs(4)),absc(cs(5))) 
    If ( xmxp.Lt.xmax ) Then
     xmax = xmxp
     cb2i(2) = csom
    Endif
   End If

  Else ! p^2=0, m1==m2:

   cb2i(1) = cb0/3
   cb2i(2) = Real(xm1/2,dp)*(cb0 + 1._dp)
  Endif

 End Function Bii

 Real(dp) Function C0_3m(m1, m2, m3)
 Implicit None
  Real(dp), Intent(in) :: m1, m2, m3

  Integer :: i1, i2, seriesOrder
  Real(dp) :: m_j, xpi1, xpi2, xpi3, d12, d13, d23, sum, LogR &
         &  , coeff(0:7), r13, r12

  Iname = Iname + 1
  NameOfUnit(Iname) = "C0_3m"

  r12 = 0._dp
  r13 = 0._dp
  
  xpi1 = m1
  xpi2 = m2
  xpi3 = m3

!   sort the masses such that m1 >= m2 >= m3
!   this is important to avoid complex logs later
   If (xpi1 .Lt. xpi2) Then
    m_j = xpi2
    xpi2 = xpi1
    xpi1 = m_j
   Endif
   If (xpi2 .Lt. xpi3) Then
    m_j = xpi3
    xpi3 = xpi2
    xpi2 = m_j
   Endif
   If (xpi1 .Lt. xpi2) Then
    m_j = xpi2
    xpi2 = xpi1
    xpi1 = m_j
   Endif

   d12 = xpi1 - xpi2
   d13 = xpi1 - xpi3
   d23 = xpi2 - xpi3

   If (xpi3.Lt.0._dp) Then ! negative mass squares not included at all
    Write(ErrCan,*) "Negative mass square are not include in "//NameOfUnit(Iname)
    Write(ErrCan,*) "m^2_3",xpi3
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    C0_3m = - Huge(1._dp)
    Iname = Iname - 1
    Return

   Else If (xpi2.Eq.0._dp) Then ! at least two masses zero -> infrared singularity
    Write(ErrCan,*) "Infrared singularity: all masses are zero in"//NameOfUnit(Iname)
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    C0_3m = Huge(1._dp)
    Iname = Iname - 1
    Return

   End If
!-------------------------
!    m3 = 0
!-------------------------
   If (xpi3.Eq.0._dp) Then
    d12 = d12 / xpi1
    If (d12.Eq.0._dp) Then ! m1 = m2
     C0_3m = -1._dp/xpi1

    Else If (d12.Lt.1e-4_dp) Then ! m1 nearly m2
     C0_3m = d12 / 9._dp
     Do i1=7,1,-1
      C0_3m = (C0_3m + 1._dp/Real(i1+1,dp)) * d12
     End Do
     C0_3m = - (1._dp + C0_3m) / xpi1

    Else
     C0_3m = Log(xpi2/xpi1)/(xpi1 - xpi2)
    End If
    Iname = Iname - 1
    Return
   End If

   d12 = d12 / xpi1
   d13 = d13 / xpi1
   d23 = d23 / xpi2
   r13 = r13+ xpi3 / xpi1
   r12 = r12+ xpi2 / xpi1

   If ((d12.Eq.0._dp).And.(d23.Eq.0._dp)) Then ! all masses equal
    C0_3m = -0.5_dp / xpi1

   Else If ( (d12.Eq.0._dp).And.(d23.Lt.1.e-4_dp)) Then ! m1=m2, m2 nearly equal m3

    C0_3m = d13 / 90._dp

    Do i1=7,1,-1
     C0_3m = ( C0_3m + 1._dp / Real((i1+1)*(i1+2),dp) ) * d13
    End Do
    C0_3m = - (0.5_dp + C0_3m) / xpi1
    
   Else If ( (d12.Eq.0._dp).And.(r13.Lt.1.e-4_dp)) Then ! m1=m2, m2 >> m3

    LogR = log(r13)

    C0_3m = (1._dp + 8._dp * LogR) * r13
    Do i1=7,1,-1
     C0_3m = ( C0_3m + 1._dp + i1 * LogR ) * r13
    End Do
    C0_3m = - (1._dp + C0_3m) / xpi1
    

   Else If ( (d12.Lt.1.e-4_dp).And.(d23.Eq.0._dp)) Then ! m1 nearly m2, m2=m3

    C0_3m = d13 / 90._dp

    Do i1=7,1,-1
     C0_3m = ( C0_3m + 1._dp / Real((i1+1)*(i1+2),dp) ) * d13
    End Do
    C0_3m = - (0.5_dp + C0_3m) / xpi1
    
   Else If ( (d23.Eq.0._dp).And.(r13.Lt.1.e-4_dp)) Then ! m2=m3, m1 >> m3

    LogR = log(r13)

    C0_3m = (1._dp + 9._dp * LogR) * r13
    Do i1=7,1,-1
     C0_3m = ( C0_3m + 1._dp + (1+i1) * LogR ) * r13
    End Do
    C0_3m = (1._dp + LogR + C0_3m) / xpi1

   Else If (d13.Lt.1.e-3_dp) Then ! all masses nearly equal

    sum = d13**9
    Do i1=1,8
     sum = sum + d12**i1 * d13**(9-i1)
    End Do
    C0_3m = (d12**8 + sum)/110._dp

    Do i1=8,1,-1
     sum = d13**i1
     Do i2=1,i1-1
      sum = sum + d12**i2 * d13**(i1-i2)
     End Do
     C0_3m = C0_3m + (d12**i1 + sum) / Real((i1+1)*(i1+2),dp)
    End Do
    C0_3m = - (0.5_dp + C0_3m) / xpi1

   Else If (d12.Lt.1.e-4_dp) Then ! m1 nearly m2

    LogR = Log(xpi3/xpi1)
    C0_3m = 0._dp
    Coeff(0) = -xpi1 + xpi3 - LogR*xpi3    
    d13 = xpi1 - xpi3
    If (d12.Gt.0._dp) Then
     Coeff(1) = (-xpi1**2 - 2*LogR*xpi1*xpi3 + xpi3**2)/2._dp
     Coeff(2) = ( -2*xpi1**3 - 3*xpi1**2*xpi3 - 6*LogR*xpi1**2*xpi3 &
             & + 6*xpi1*xpi3**2 - xpi3**3)/6._dp
     Coeff(3) = (-3*xpi1**4 - 10*xpi1**3*xpi3 - 12*LogR*xpi1**3*xpi3 &
              & + 18*xpi1**2*xpi3**2 - 6*xpi1*xpi3**3 + xpi3**4)/12._dp
     Coeff(4) = -xpi1**5/5._dp - (13*xpi1**4*xpi3)/12._dp - LogR*xpi1**4*xpi3 &
              & + 2*xpi1**3*xpi3**2 - xpi1**2*xpi3**3 + (xpi1*xpi3**4)/3._dp  &
              & - xpi3**5/20._dp
     Coeff(5) = ( -10*xpi1**6 - 77*xpi1**5*xpi3 - 60*LogR*xpi1**5*xpi3          &
             & + 150*xpi1**4*xpi3**2 - 100*xpi1**3*xpi3**3 + 50*xpi1**2*xpi3**4 &
             & - 15*xpi1*xpi3**5 + 2*xpi3**6)/60._dp
     Coeff(6) = -xpi1**7/7._dp - (29*xpi1**6*xpi3)/20._dp - LogR*xpi1**6*xpi3 &
             & + 3*xpi1**5*xpi3**2 - (5*xpi1**4*xpi3**3)/2._dp                &
             & + (5*xpi1**3*xpi3**4)/3._dp - (3*xpi1**2*xpi3**5)/4._dp       &
             & + (xpi1*xpi3**6)/5._dp - xpi3**7/42._dp
     Coeff(7) = (223*xpi1**7*xpi3)/140._dp - LogR*xpi1**7*xpi3           &
             & + (7*xpi1**6*xpi3**2)/2._dp - (7*xpi1**5*xpi3**3)/2._dp  &
             & + (35*xpi1**4*xpi3**4)/12._dp - (7*xpi1**3*xpi3**5)/4._dp &
             & + (7*xpi1**2*xpi3**6)/10._dp - (xpi1*xpi3**7)/6._dp + xpi3**8/56._dp
     C0_3m = Coeff(7) * d12 / d13

     Do i1=6,1,-1

      C0_3m = (Coeff(i1) + C0_3m) * d12 / d13
     End Do
    End If ! d12 .gt. 0._dp

    C0_3m = (Coeff(0) + C0_3m) / d13**2

   Else If ((r12.Lt.1.e-4_dp).And.(r13.Lt.1.e-4_dp)) Then 
     LogR = Log(xpi2/xpi3)
     C0_3m = 0._dp
     seriesOrder = 10
     ! 13. Ordnung in RIJ
     If (13.Le.seriesOrder) Then
     Do i1=13,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(13-i1)
     End Do       
     End If
     ! 12. Ordnung in RIJ   
     If (12.Le.seriesOrder) Then     
     Do i1=12,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(12-i1)
     End Do       
     Do i1=13,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-12)
     End Do   
     End If 
     ! 11. Ordnung in RIJ
     If (11.Le.seriesOrder) Then
     Do i1=11,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(11-i1)
     End Do       
     Do i1=12,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-11)
     End Do  
     End If
     ! 10. Ordnung in RIJ
     If (10.Le.seriesOrder) Then     
     Do i1=10,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(10-i1)
     End Do       
     Do i1=11,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-10)
     End Do      
     End If
     ! 9. Ordnung in RIJ
     If (9.Le.seriesOrder) Then  
     Do i1=9,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(9-i1)
     End Do       
     Do i1=10,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-9)
     End Do       
     End If 
     ! 8. Ordnung in RIJ
     If (8.Le.seriesOrder) Then
     Do i1=8,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(8-i1)
     End Do       
     Do i1=9,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-8)
     End Do   
     End If
     ! 7. Ordnung in rIJ
     If (8.Le.seriesOrder) Then
     Do i1=7,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(7-i1)
     End Do   
     Do i1=8,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-7)
     End Do        
     End If
     ! 6. Ordnung in rIJ
     If (6.Le.seriesOrder) Then
     Do i1=6,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(6-i1)
     End Do     
     Do i1=7,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-6)
     End Do       
     End If
     ! 5. Ordnung in rIJ
     If (5.Le.seriesOrder) Then
     Do i1=5,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(5-i1)
     End Do
     Do i1=6,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-5)
     End Do      
     End If
     ! 4. Ordnung in rIJ
     If (4.Le.seriesOrder) Then
     Do i1=4,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(4-i1)
     End Do  
     Do i1=5,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-4)
     End Do     
     End If
     ! 3. Ordnung in rIJ
     If (3.Le.seriesOrder) Then
     Do i1=3,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(3-i1)
     End Do  
     Do i1=4,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-3)
     End Do  
     End If
     ! 2. Ordnung in rIJ
     If (2.Le.seriesOrder) Then
     Do i1=2,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(2-i1)
     End Do       
     Do i1=3,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-2)
     End Do
     End If
     ! 1. Ordnung in rIJ
     If (1.Le.seriesOrder) Then
     Do i1=1,0,-1
       C0_3m = C0_3m + Log(r12)*r12**i1*r13**(1-i1)
     End Do       
     Do i1=2,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**(i1-1)
     End Do   
     End If
     ! 0. Ordnung in rIJ
     C0_3m = C0_3m + Log(r12)
     Do i1=1,seriesOrder,1
       C0_3m = C0_3m + LogR*r13**i1/r12**i1
     End Do       
     
     C0_3m = C0_3m/xpi1

   
   Else If (d23.Lt.1.e-4_dp) Then ! m2 nearly m3
    LogR = Log(xpi2/xpi1)
    Coeff(0) = xpi1 + LogR*xpi1 - xpi2
    C0_3m = 0._dp
    d12 = xpi1 - xpi2
    If (d23.Gt.0._dp) Then
     Coeff(1) = (-xpi1**2 - 2*LogR*xpi1*xpi2 + xpi2**2)/2._dp
     Coeff(2) = ( -xpi1**3 + 6*xpi1**2*xpi2 - 3*xpi1*xpi2**2  &
             & + 6*LogR*xpi1*xpi2**2 - 2*xpi2**3)/6._dp
     Coeff(3) = (-xpi1**4 + 6*xpi1**3*xpi2 - 18*xpi1**2*xpi2**2 &
             & + 10*xpi1*xpi2**3 - 12*LogR*xpi1*xpi2**3 + 3*xpi2**4)/12._dp
     Coeff(4) = -xpi1**5/20._dp + (xpi1**4*xpi2)/3._dp - xpi1**3*xpi2**2 &
              & + 2*xpi1**2*xpi2**3 - (13*xpi1*xpi2**4)/12._dp           &
              & + LogR*xpi1*xpi2**4 - xpi2**5/5._dp
     Coeff(5) = (-2*xpi1**6 + 15*xpi1**5*xpi2 - 50*xpi1**4*xpi2**2           &
             & + 100*xpi1**3*xpi2**3 - 150*xpi1**2*xpi2**4 + 77*xpi1*xpi2**5 &
             & - 60*LogR*xpi1*xpi2**5 + 10*xpi2**6)/60._dp
     Coeff(6) = -xpi1**7/42._dp + (xpi1**6*xpi2)/5._dp                  &
             & - (3*xpi1**5*xpi2**2)/4._dp  + (5*xpi1**4*xpi2**3)/3._dp &
             & - (5*xpi1**3*xpi2**4)/2. + 3*xpi1**2*xpi2**5             &
             & - (29*xpi1*xpi2**6)/20._dp + LogR*xpi1*xpi2**6 - xpi2**7/7._dp
     Coeff(7) = -xpi1**8/56._dp + (xpi1**7*xpi2)/6._dp                   &
             & - (7*xpi1**6*xpi2**2)/10._dp + (7*xpi1**5*xpi2**3)/4._dp  &
             & - (35*xpi1**4*xpi2**4)/12._dp + (7*xpi1**3*xpi2**5)/2._dp &
             & - (7*xpi1**2*xpi2**6)/2._dp + (223*xpi1*xpi2**7)/140._dp  &
             & - LogR*xpi1*xpi2**7 + xpi2**8/8._dp

     C0_3m = Coeff(7) * d23 / d12

     Do i1=6,1,-1
      C0_3m = (Coeff(i1) + C0_3m) * d23 / d12
     End Do

    End If ! d23 .gt. 0._dp

    C0_3m = (Coeff(0) + C0_3m) / d12**2

   Else ! all masses are different 
    C0_3m = (Log(xpi3/xpi2) + xpi1/(xpi3-xpi1)*Log(xpi3/xpi1) &
           & - xpi1/(xpi2 - xpi1)*Log(xpi2/xpi1))/(xpi2 - xpi3)
   End If
   
   Iname = Iname - 1
  
 End Function C0_3m


 Complex(dp) Function C00_3m(a,b,c)
  Real(dp), Intent(in) :: a,b,c
  C00_3m = 0.125_dp + 0.25_dp*vertexC0tilde(a,b,c) 

 End function C00_3m


 Complex(dp) Function C0(p1, p2, p1p2, m1, m2, m3)
 !-----------------------------------------------------------------------
 !
 !	Calculates the threepoint function closely following
 !	recipe in 't Hooft & Veltman, NP B(183) 1979.
 !	Bjorken and Drell metric is used nowadays!
 !				
 !	    p2	| |		
 !		v |		
 !		 / \		
 !	      m2/   \m3 	
 !	p1     /     \	p3	
 !	->    /  m1   \ <-	
 !	------------------------
 !				
 !		1   /			     1			
 !	    = ----- \d^4Q----------------------------------------
 !	      ipi^2 /	 [Q^2-m1^2][(Q+p1)^2-m2^2][(Q-p3)^2-m3^2]
 !
 !	If the function is infra-red divergent (p1=m2,p3=m3,m1=0 or
 !	cyclic) the function is calculated with a user-supplied cutoff
 !	lambda2 in the common block /ffregul/.
 !
 !	Input:	xpi	(real)		i=1,3: mass^2, i=4,6: pi.pi
 !	Output: cc0	(complex)	C0, the threepoint function.
 !									*
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: p1, p2, p1p2, m1, m2, m3

  Integer :: i1, i2, ier
  Real(dp) :: xpi(6), dpipj(6,6), maxM

  xpi(1) = m1
  xpi(2) = m2
  xpi(3) = m3
  xpi(4) = p1
  xpi(5) = p2
  xpi(6) = p1p2
  maxm = Minval(Abs(xpi(1:3)))
  ! all momenta are 0 or very small
  If (Sum(Abs(xpi(4:6))/maxm) .Lt. 1.e-10_dp ) Then
   C0 = Cmplx( C0_3m(m1, m2, m3), 0._dp, dp)

  Else 

   Do i1=1,6
    dpipj(i1,i1) = 0._dp
    Do i2 = i1+1,6
     dpipj(i1,i2) = xpi(i1) - xpi(i2)
     dpipj(i2,i1) = - dpipj(i1,i2)
    End Do
   End Do

 ! infrared divergent diagrams
   If (     ((dpipj(2,4).Eq.0) .And. (dpipj(3,6).Eq.0) .And. (xpi(1).Eq.0))  &
     & .Or.((dpipj(3,5).Eq.0) .And. (dpipj(1,4).Eq.0) .And. (xpi(2).Eq.0))   &
     & .Or.((dpipj(1,6).Eq.0) .And. (dpipj(2,5).Eq.0) .And. (xpi(3).Eq.0)) ) &
  & Then

    Call ffxc0i(C0, xpi, dpipj, ier)
    If(ier > 9) Print *, "C0 lost ", ier, " digits (m1 = ", Sqrt(m1), ")"
     
!    end if
   Else ! non-divergent case
    C0 = C0_regular(xpi, dpipj)
   End If
  Endif

 End Function C0


 Complex(dp) Function C0_regular(xpi, dpipj)
 Implicit None
  Real(dp), Intent(inout) :: xpi(6), dpipj(6,6)

  Real(dp) :: xqi(6), dqiqj(6,6), qiDqj(6,6), dum66(6,6), del2, fdel2, del3 &
    & , del2s(3), delpsi(3), del3mi(3), xmax, alph(3),etalam,etami(6),sdel2 &
    & , blph(3)
  Complex(dp) :: cs3(80),cs,clogi(3),cslam,cetalm,cetami(6),cel2s(3)      &
    & ,calph(3),cblph(3),csdel2,cqi(6),cdqiqj(6,6),cqiDqj(6,6),celpsi(3)
  Integer :: i, j, k,isoort(8),ipi12(8),ipi12t,ilogi(3)
  Integer, Parameter :: inew(6,6) = Reshape( &
     & Source =  (/1,2,3,4,5,6,  &
     &             2,3,1,5,6,4,  &
     &             3,1,2,6,4,5,  &
     &             1,3,2,6,5,4,  &
     &             3,2,1,5,4,6,  &
     &             2,1,3,4,6,5/), Shape = (/6, 6/) )
  dum66 = 0._dp
  Call ffrot3(irota3,xqi,dqiqj,qiDqj,xpi,dpipj,dum66,2,3)
  Call ffdot3(qiDqj, xqi, dqiqj)
  ! save dotproducts for tensor functions if requested

  If ( ldot ) Then
   Do i=1,6
    Do j=1,6
     fpij3(j,i) = qiDqj(inew(i,irota3),inew(j,irota3))
    End Do
   End Do
   If ( irota3 .Gt. 3 ) Then ! the sign of the s's has been changed!
    Do i=1,3
     Do j=4,6
      fpij3(j,i) = -fpij3(j,i)
      fpij3(i,j) = -fpij3(i,j)
     End Do
    End Do
   Endif
  Endif

 ! some determinats
  del2 = FFdel2(qiDqj,4,5,6)

  If ( ldot ) fdel2 = del2

  If ( del2 .Gt. 0 ) Then
   If ( .Not.(xqi(4).Lt.0 .And. xqi(5).Lt.0 .And. xqi(6).Lt.0) ) Then
    Call WriteLFerror(15)
   Endif
  Elseif ( del2 .Eq. 0 ) Then
   Call WriteLFerror(16)
   C0_regular = 0._dp
   Return
  Endif

  del3 = ffdel3(qiDqj)

  Call ffdl3m(del3mi,.True.,del3,del2,xqi,dqiqj,qiDqj,6, 4,5,6,1,3)
  Do i=1,3
   j = i+1
   If ( j .Eq. 4 ) j = 1
   del2s(i) = ffdel2(qiDqj,i+3,i,j)
   k = i-1
   If ( k .Eq. 0 ) k = 3
   delpsi(i) = ffdl2p(xqi,dqiqj,qiDqj,i+3,j+3,k+3,i,j,k)
  End Do
  ! initialize cs3:
  cs3 = 0._dp
  ipi12 = 0
  clogi = 0._dp
  ilogi = 0 

  If ( del2 .Gt. 0 ) Then 
  ! complex case, in case of three spacelike momenta or unphysical real ones
   Do i=1,3
    cel2s(i) = del2s(i)
    celpsi(i) = delpsi(i)
    cetami(i) = del3mi(i)/del2
   End Do
   Do i=1,6
    cqi(i) = xqi(i)
    Do j=1,6
     cdqiqj(j,i) = dqiqj(j,i)
     cqiDqj(j,i) = qiDqj(j,i)
    End Do
   End Do
   cetalm = del3/del2
   csdel2 = isgnal*Cmplx(0._dp,Sqrt(del2), dp)

   !   get alpha,1-alpha
   Call roots(cqi(5),-cqiDqj(5,6),cqi(6),csdel2,cblph(1),calph(1))
   Call roots(cqi(5),-cqiDqj(5,4),cqi(4),csdel2,calph(3),cblph(3))
   cslam = 2*csdel2
   Call ffcc0p(cs3,ipi12,isoort,clogi,ilogi,cqi,cdqiqj,cqiDqj &
             &, csdel2,cel2s,cetalm,cetami,celpsi,calph,3)

  Else ! real case
   etalam = del3/del2
   Do i=1,3
    etami(i) = del3mi(i)/del2
   End Do
   If ( Abs(isgnal).Ne.1 ) Then
    Print *,'c0_regular: error: isgnal should be +/-1, not ',isgnal
    Print *,'        forgot to call InitializeLoopFunctions?'
    Call InitializeLoopFunctions
   Endif
   sdel2 = isgnal*Sqrt(-del2)

   !  get alpha,1-alpha
   Call roots(xqi(5),-qiDqj(5,6),xqi(6),sdel2,blph(1),alph(1))
   Call roots(xqi(5),-qiDqj(5,4),xqi(4),sdel2,alph(3),blph(3))

   If ( l4also .And. ( alph(1) > 1 .Or. alph(1) < 0 ) .And.&
    &  Abs(blph(1)-0.5_dp) < Abs(alph(1)-0.5_dp) ) Then
    alph(1) = blph(1)
    alph(3) = blph(3)
    sdel2 = -sdel2
    isgnal = -isgnal
   End If
   cslam = 2*sdel2
   Call ffxc0p(cs3,ipi12,isoort,clogi,ilogi,xqi,dqiqj,qiDqj &
              &, sdel2,del2s,etalam,etami,delpsi,alph,3)
  End If

  cs = 0
  xmax = 0

  Do i=1,80
   cs = cs + cs3(i)
   xmax = Max(xmax,absc(cs))
  End Do

  ipi12t = Sum(ipi12)

  cs = cs + ipi12t*Pi2o12
  ! A imaginary component less than precc times the real part is zero
  If ( Abs(Aimag(cs)) .Lt. precc*Abs(Real(cs,dp)) )  &
     & cs = Cmplx(Real(cs,dp),0._dp, dp)
  C0_regular = - cs/cslam

 End Function C0_regular


 Real(dp) Function C_2(x, y, z)
 !-------------------------------------------------------------------
 ! this is the function C_2 as defined in Buras et al.hep-ph/0210145
 ! written by Werner Porod, 22.03.03
 !-------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x, y, z

  Real(dp) :: m1, m2, m3, diff1, diff2, sum1
  Integer :: i1, i2, n1
  !----------------------
  ! first sorting masses
  !----------------------
  If ((x.Ge.y).And.(x.Ge.z)) Then
   m3 = x
   If (y.Ge.z) Then
    m2 = y
    m1 = z
   Else
    m2 = z
    m1 = y
   End If

  Else If ((y.Ge.x).And.(y.Ge.z)) Then
   m3 = y
   If (x.Ge.z) Then
    m2 = x
    m1 = z
   Else
    m2 = z
    m1 = x
   End If

  Else
   m3 = z
   If (x.Ge.y) Then
    m2 = x
    m1 = y
   Else
    m2 = y
    m1 = x
   End If

  End If

  !-----------------------------------
  ! m1=m2=0
  !-----------------------------------
  If ((m1.Eq.0._dp).And.(m2.Eq.0._dp)) Then
   C_2 = Log(m3/mudim2)

  !-----------------------------------
  ! m1=0
  !-----------------------------------
  Else If (m1.Eq.0._dp) Then
   diff1 = (m3 - m2) / m2
   If (diff1.Eq.0._dp) Then ! m2=m3
    C_2 = 1._dp + Log(m3/mudim2)

   Else If (diff1.Le.1.e-3_dp) Then ! m2 and m3 agree with per-mille
    n1 = Log10(diff1/Epsilon(1._dp)) + 2 ! maximal number of terms in the serie
    C_2 = (-1)**(n1+1) * diff1 / Real(n1*(n1+1),dp)
    Do i1= n1-1,1,-1
     C_2 = diff1 * (C_2 + (-1)**(i1+1) / Real(i1*(i1+1),dp) )
    End Do
    C_2 = C_2 + 1._dp + Log(m3/mudim2)

   Else
    C_2 = Log(m3/mudim2) + m2 * Log(m3/m2) / (m3-m2)
   End If

  !---------------------------------
  ! all masses are different from 0
  !---------------------------------
  Else ! general formula
   diff1 = (m3 - m2) / m2
   diff2 = (m2 - m1) / m1

   If ((diff1.Eq.0._dp).And.(diff2.Eq.0._dp)) Then ! m1=m2=m3
    C_2 = 1.5_dp + Log(m3/mudim2)

   Else If ((diff1.Eq.0._dp).And.(diff2.Lt.1.e-3_dp)) Then ! m2=m3 \simeq m1
    n1 = Log10(diff2/Epsilon(1._dp)) + 2 ! maximal number of terms in the serie
    C_2 = (-1)**(n1+1) * diff2 * 2._dp / Real(n1*(n1+2),dp)
    Do i1= n1-1,1,-1
     C_2 = diff2 * (C_2 + (-1)**(i1+1) * 2._dp / Real(i1*(i1+2),dp) )
    End Do
    C_2 = C_2 + 1.5_dp + Log(m1/mudim2)
    
   Else If (diff1.Eq.0._dp) Then ! m2=m3 \= m1
    C_2 = m3 / (m3-m1) + Log(m3/mudim2) + m1**2 * Log(m1/m3) / (m3-m1)**2

   Else If ((diff2.Eq.0._dp).And.(diff1.Lt.1.e-3_dp)) Then ! m2=m1 \simeq m3
    n1 = Log10(diff1/Epsilon(1._dp)) + 2 ! maximal number of terms in the serie
    C_2 = (-1)**(n1+1) * diff1 * 2._dp / Real(n1*(n1+1)*(n1+2),dp)
    Do i1= n1-1,1,-1
     C_2 = diff1 * (C_2 + (-1)**(i1+1) * 2._dp / Real(i1*(i1+1)*(i1+2),dp) )
    End Do
    C_2 = C_2 + 1.5_dp + Log(m2/mudim2)

   Else If (diff2.Eq.0._dp) Then ! m2=m1 \= m3
    C_2 = m1 / (m1-m3) + Log(m3/mudim2) &
        &              + m1*(m1-m3*2._dp) * Log(m1/m3) /(m3-m1)**2

    ! all masses similar
   Else If ((diff1.Lt.1.e-3_dp).And.(diff2.Lt.1.e-3_dp)) Then
    diff1 = (m3 - m1) / m1
    ! maximal number of terms in the serie
    n1 = Log10(Min(diff1,diff2)/Epsilon(1._dp)) + 2
    C_2 = 0
    Do i1=n1,1-1
     sum1 = diff1**i1
     Do i2=1,i1
      sum1 = sum1 + diff1**(i1-1) * diff2**i2
     End Do
     C_2 = C_2 + sum1 * (-1)**(i1+1) * 2._dp / Real(i1*(i1+1)*(i1+2),dp)
    End Do
    C_2 = C_2 + 1.5_dp + Log(m1/mudim2)

   Else   
    C_2 = Log(m3/mudim2) + m2**2 * Log(m3/m2) / ((m3-m2)*(m2-m1))  &
        &               - m1**2 * Log(m3/m1) / ((m3-m1)*(m2-m1))
   End If
  End If

 End Function C_2

 Real(dp) Function C11_Dedes(x, y)
 !----------------------------------------------------------------------------
 ! loop function taken from A.Dedes et al. PRD 79 (2009) 055006
 !----------------------------------------------------------------------------
 Implicit None 

  Real(dp), Intent(in) :: x, y

  If (x.eq.y) then

   C11_Dedes = - 1._dp / (6._dp * x)

  Else If (y.eq.0._dp) then

   C11_Dedes = - 0.25_dp / x

  else

   C11_Dedes = ( 0.25_dp * (3._dp*y-x ) + y**2*Log(y/x)/(2.*(x - y)) ) &
           & / (x - y)**2
   
  End If

 End Function C11_Dedes

 Real(dp) Function C12_Dedes(x, y)
 !----------------------------------------------------------------------------
 ! loop function taken from A.Dedes et al. PRD 79 (2009) 055006
 !----------------------------------------------------------------------------
 Implicit None 

  Real(dp), Intent(in) :: x, y

  If (x.eq.0._dp) then

   C12_Dedes = - 0.25_dp / y

  Else If (y.eq.0._dp) then

   C12_Dedes = - 0.25_dp / x

  Else

   C12_Dedes = - (0.5_dp * (x+y) + x*y*Log(y/x) / (3.*(x - y) ) ) / (x - y)**2

  end If

 End Function C12_Dedes

 Complex(dp) Function Cget(name, p1, p2, p1p2, m1, m2, m3)
 Implicit None
  Real(dp), Intent(in) :: p1, p2, p1p2, m1, m2, m3
  Character(len=4), Intent(in) :: name

  Complex(dp) :: C_0, C1, C2, C00, C11, C12, C22, C001, C002, C111, C112  &
    & , C122, C222

  Iname = Iname + 1
  NameOfUnit(Iname) = "Cget"

  Call Get_All_Ci(p1, p2, p1p2, m1, m2, m3, C_0, C1, C2, C00, C11, C12 &
    & , C22, C001, C002, C111, C112, C122, C222)

  If (name.Eq."C0") Then
   Cget = C_0

  Else If (name.Eq."C1") Then
   Cget = C1

  Else If (name.Eq."C2") Then
   Cget = C2

  Else If (name.Eq."C00") Then
   Cget = C00

  Else If (name.Eq."C11") Then
   Cget = C11

  Else If (name.Eq."C12") Then
   Cget = C12

  Else If (name.Eq."C22") Then
   Cget = C22

  Else If (name.Eq."C001") Then
   Cget = C001

  Else If (name.Eq."C002") Then
   Cget = C002

  Else If (name.Eq."C111") Then
   Cget = C111

  Else If (name.Eq."C112") Then
   Cget = C112

  Else If (name.Eq."C122") Then
   Cget = C122

  Else If (name.Eq."C222") Then
   Cget = C222

  Else
   Write (ErrCan,*) "Problem in function Cget, function "//name
   Write (ErrCan,*) "is not defined"
   If (Errorlevel.Ge.0) Call TerminateProgram
  End If

  Iname = Iname - 1

 End Function Cget

 Subroutine ClearCache()
 Implicit None 
  num_ci = 0
 End Subroutine ClearCache

 Real(dp) Function D0_Bagger(m12, m22, m32, m42)
 !-----------------------------------------------------------------------
 ! calculates the function D_0 as defined in J. Bagger at al, Nucl.Phys.B
 ! written by Werner Porod, 12.8.1999
 ! 18.05.2001: porting to f90
 ! 21.09.09: adding expansion for the case m12=m22, m32=m42
 !                 m12 ~= m32
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: m12, m22, m32, m42

  Integer, Parameter :: ord=9
  Real(dp) :: Dm, x, ci(ord), di(ord), diff31, diff41, diff32 &
    & , xl2, xh2, yl2, yh2
  !-------------------------------------------------------------------
  ! recursion s1(i1+1) = (-1)^i1 * s1(i1) * i1 (i1+3)/(i1+1)**2
  !-------------------------------------------------------------------
  Real(dp), Parameter :: s1(0:16) = (/ 1._dp, -1._dp, 0.9_dp, -0.8_dp           &
     & , 5._dp/7._dp, -9._dp/14._dp, 7._dp/12._dp, -8._dp/15._dp, 27._dp/55._dp &
     & , -5._dp/11._dp, 11._dp/26._dp, -36._dp/91._dp, 13._dp/35._dp            &
     & , -0.35_dp, 45._dp/136._dp, -16._dp/51._dp, 17._dp/57._dp /)
  Integer :: i1, i2

  
  If (     ((m12.eq.0._dp).and.(m22.eq.0._dp))  &
     & .or.((m12.eq.0._dp).and.(m32.eq.0._dp))  &
     & .or.((m12.eq.0._dp).and.(m42.eq.0._dp))  &
     & .or.((m22.eq.0._dp).and.(m32.eq.0._dp))  &
     & .or.((m22.eq.0._dp).and.(m42.eq.0._dp))  &
     & .or.((m32.eq.0._dp).and.(m42.eq.0._dp)) )  then
   D0_Bagger = Huge(1._dp) ! infinity
   return
  Else ! resorting according to the symmetries
       ! m12 <-> m22, m32 <-> m42


   If (m12.le.m22) then
    xl2 = m12
    xh2 = m22
   Else
    xl2 = m22
    xh2 = m12
   end if

   If (m32.le.m42) then
    yl2 = m32
    yh2 = m42
   Else
    yl2 = m42
    yh2 = m32
   end if

  End If

  If ((m12.Eq.m22).And.(m32.Eq.m42)) Then

   Dm = Abs(m12/m42 - 1._dp)

   If (m12.Eq.m42) Then
    D0_Bagger = 1._dp / (6._dp * m12**2)
   Else If (Dm.Lt.1.e-2_dp) Then
    Dm = (m42-m12) / m12
    D0_Bagger = s1(16) * Dm
    Do i1=15,1,-1
     D0_Bagger = (D0_Bagger + s1(i1)) * Dm
    End Do
    D0_Bagger = (D0_Bagger + s1(0)) / (6._dp * m12**2)
   Else
    D0_Bagger = (2._dp * (m42-m12) - (m12+m42) * Log(m42/m12) )  &
           & / (m12-m42)**3
   End If

  Else If ((m12.Eq.m22).And.(m12.Eq.m42)) Then

   If (Abs((m32-m42)/m42).Lt.1.e-3_dp) then
    x = 1._dp - m42/m32 
    D0_Bagger = 1._dp / Real((ord+1)*(ord+2),dp)
    Do i1=ord-1,1,-1
     D0_Bagger = D0_Bagger * x + 1._dp / Real((i1+1)*(i1+2),dp)
    End Do
    D0_Bagger = D0_Bagger / m32**2

   Else
    D0_Bagger =  (m32**2 - m42**2 - 2._dp*m32*m42*Log(m32/m42))  &
           & / (2._dp*(m32 - m42)**3*m42)
   End If

  Else If ((m12.Eq.m32).And.(m12.Eq.m42)) Then

   If (Abs((m12-m22)/m22).Lt.1.e-3_dp) then
    x = 1._dp - m22/m12 
    D0_Bagger = 1._dp / Real((ord+1)*(ord+2),dp)
    Do i1=ord-1,1,-1
     D0_Bagger = D0_Bagger * x + 1._dp / Real((i1+1)*(i1+2),dp)
    End Do
    D0_Bagger = D0_Bagger / m12**2

   Else
    D0_Bagger =  (m12**2 - m22**2 - 2._dp*m12*m22*Log(m12/m22))  &
            & / (2._dp*(m12 - m22)**3*m12)
   End If

  Else If ((m12.Eq.m22).And.(m12.Eq.m32)) Then

   If (Abs((m32-m42)/m42).Lt.1.e-3_dp) then
    x = 1._dp - m42/m32 
    D0_Bagger =  Real(ord,dp) / Real(2*(ord+2),dp)
    Do i1=ord-1,1,-1
     D0_Bagger = D0_Bagger * x + Real(i1,dp) / Real(2*(i1+2),dp)
    End Do
    D0_Bagger = D0_Bagger / m32**2

   Else
    D0_Bagger = (m32**2 - m42**2 + 2._dp*m32*m42*Log(m42/m32)) &
            & / (2._dp*m32*(m32 - m42)**3)
   End If

  Else If ((m22.Eq.m32).And.(m22.Eq.m42)) Then

   If (Abs((m12-m22)/m22).Lt.1.e-3_dp) then
    x = 1._dp - m22/m12 
    D0_Bagger =  Real(ord,dp) / Real(2*(ord+2),dp)
    Do i1=ord-1,1,-1
     D0_Bagger = D0_Bagger * x + Real(i1,dp) / Real(2*(i1+2),dp)
    End Do
    D0_Bagger = D0_Bagger / m12**2

   Else
    D0_Bagger =  (m12**2 - m22**2 + 2._dp*m12*m22*Log(m22/m12))  &
            & / (2._dp*(m12 - m22)**3*m22)
   End If

  Else If (m12.Eq.m22) Then

   If ( (Abs((m42-m32)/m42)).Lt.1.e-5_dp) then
    x = m42 / m32 - 1._dp
    diff31 = m32 - m12

    di(1) = m12+m32

    Do i1=1,ord-1
     di(i1+1) = (-1)**i1 * m12 * (m32/diff31)**i1 &
            & - (m32*di(i1))/diff31
    End Do
    di = di * Log(m32/m12) / diff31**3

    Do i1=1,ord
     ci(i1) = 0._dp
     Do i2=1,i1-1
      ci(i1) = ci(i1) + (-1)**i2 * i2 * (m32/diff31)**(i2-1) &
             &           / Real( (i1-i2+1)*(i1-i2),dp) 
     End Do
     ci(i1) = ci(i1) + (-1)**i1 * (i1+1) * (m32/diff31)**(i1-1)
    End Do

    ci = ci / diff31**2

    D0_Bagger = ci(ord) + di(ord)
    Do i1=ord-1,1,-1
     D0_Bagger = D0_Bagger * x + ci(i1) + di(i1)
    End Do

   Else
    D0_Bagger = ( m42 * Log(m42/m12) / (m12-m42)**2                 &
              & - m32 * Log(m32/m12) / (m12-m32)**2 ) / (m32-m42)   &
              & - 1._dp / ( (m12-m32) * (m12-m42) )
   end If

  Else If (m32.Eq.m42) Then

   If ((Abs(m12-m22)/m22).Lt.1.e-5_dp) Then
    x = m22 / m12 - 1._dp
    diff31 = m32 - m12 

    ci = 0._dp
    If (ord.gt.9) then
     Write(ErrCan,*) "Warning from D0_Bagger: for m32=m42, m22=m12(1+x)"
     Write(ErrCan,*) "the ci coefficients are only included up to order 9"
    end if
   
    ci(1) = - 2._dp
    ci(2) = 2.5_dp * m12 + 0.5_dp * m32
    ci(3) = (-17._dp*m12**2 - 8._dp*m12*m32 + m32**2 )/6._dp
    ci(4) = (37*m12**3 + 29*m12**2*m32 - 7*m12*m32**2 + m32**3)/12._dp
    ci(5) = (-197*m12**4)/60._dp - 3.7_dp*m12**3*m32 + 1.3_dp*m12**2*m32**2    &
          & - (11*m12*m32**3)/30._dp + m32**4/20._dp
    ci(6) = 3.45_dp*m12**5 + 5.15_dp*m12**4*m32 - 2.35_dp*m12**3*m32**2        &
          & + (59*m12**2*m32**3)/60._dp - (4*m12*m32**4)/15._dp + m32**5/30._dp
    ci(7) = (-503*m12**6)/140._dp - (236*m12**5*m32)/35._dp                    &
          & + (263*m12**4*m32**2)/70._dp - (218*m12**3*m32**3)/105._dp         &
          & + (353*m12**2*m32**4)/420._dp - (22*m12*m32**5)/105._dp            &
          & + m32**6/42._dp 
    ci(8) = (1041*m12**7)/280._dp + (2369*m12**6*m32)/280._dp                  &
          & - (1551*m12**5*m32**2)/280._dp + (3187*m12**4*m32**3)/840._dp      &
          & - (571*m12**3*m32**4)/280._dp + (213*m12**2*m32**5)/280._dp        &
          & - (29*m12*m32**6)/168._dp + m32**7/56._dp
    ci(9) = (-9649*m12**8)/2520._dp - (2593*m12**7*m32)/252._dp                &
          & + (1943*m12**6*m32**2)/252._dp - (1585*m12**5*m32**3)/252._dp      &
          & + (1061*m12**4*m32**4)/252._dp - (2633*m12**3*m32**5)/1260._dp     &
          & + (179*m12**2*m32**6)/252._dp - (37*m12*m32**7)/252._dp            &
          & + m32**8/72._dp
  
    Do i1=1,ord
     di(i1) = (m12/diff31)**(i1-1) * (m12 + i1*m32)/diff31**3
     ci(i1) = ci(i1) / (-diff31)**(i1+1)
    end do

    di = di * Log(m32/m12)

    D0_Bagger = ci(ord) + di(ord)
    Do i1=ord-1,1,-1
     D0_Bagger = D0_Bagger * x + ci(i1) + di(i1)
    End Do

   Else
    D0_Bagger = ( m12 * Log(m42/m12) / (m12-m42)**2                &
             & - m22 * Log(m42/m22) / (m22-m42)**2 ) / (m12-m22)   &
             & + 1._dp / ( (m12-m42) * (m42-m22) )
   End If

  Else If ((Abs(m12-m22)/m22).Lt.1.e-5_dp) Then

   x = m22 / m12 - 1._dp
   diff31 = m32 - m12 
   diff41 = m42 - m12 
   ci(1) = -1._dp/diff31 + (Log(m32/m12)*m12)/diff31**2
   di(1) = 1._dp/diff41 - (Log(m42/m12)*m12)/diff41**2
   Do i1=1,ord-1
    ci(i1+1) = (-1)**(1+i1) / (diff31 * (1+i1)) + (m12*ci(i1))/diff31
    di(i1+1) = (-1)**i1 / (diff41 * (1+i1)) + (m12*di(i1))/diff41
   End do
   ci = ci * m32
   di = di * m42
   D0_Bagger = ci(ord) + di(ord)
   Do i1=ord-1,1,-1
    D0_Bagger = D0_Bagger * x + ci(i1) + di(i1)
   End Do
   D0_Bagger = D0_Bagger / ( (m42 - m32) * m12)

  Else If ((Abs(m42-m32)/m42).Lt.1.e-5_dp) Then
   x = m42 / m32 - 1._dp

   diff31 = m32 - m12 
   diff32 = m32 - m22 

   ci(1) = (-1._dp - Log(m32/m12)*(1._dp - m32/diff31) ) / diff31
   di(1) = ( 1._dp + Log(m32/m22)*(1._dp - m32/diff32) )/diff32 

   Do i1=1,ord-1
    ci(i1+1) = (-1)**i1/(diff31*i1*(1 + i1)) - (m32*ci(i1))/diff31
    di(i1+1) = (-1)**(1 + i1)/(diff32*i1*(1 + i1)) - (m32*di(i1))/diff32
   End Do

   D0_Bagger = ci(ord) + di(ord)

   Do i1=ord-1,1,-1
    D0_Bagger = D0_Bagger * x + ci(i1) + di(i1)

   End Do

   D0_Bagger = D0_Bagger / (m12 - m22)

  Else
   D0_Bagger = (C0_3m(m12,m32,m42) - C0_3m(m22,m32,m42)) / (m12 - m22)
  End If

 End Function D0_Bagger


 Real(dp) Function D27_Bagger(m12,m22,m32,m42)
 !-----------------------------------------------------------------------
 ! calculates the function D_0 as defined in J. Bagger at al, Nucl.Phys.B
 ! written by Werner Porod, 12.8.1999
 ! 18.05.2001: porting to f90
 ! 21.09.09: adding expansion for the case m12=m22, m32=m42
 !                 m12 ~= m32
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: m12,m22,m32,m42
  Integer, Parameter :: ord=9
  Real(dp) :: Dm, x, ci(ord), di(ord), diff31, diff41, diff32
  !-------------------------------------------------------------------
  ! recursion s1(i1+1) = (-1)^i1 * s1(i1) * (i1+3)/(i1+1)
  !-------------------------------------------------------------------
  Real(dp), Parameter :: s1(0:16) = (/ - 1._dp, 0.5_dp, - 0.3_dp, 0.2_dp       &
    & , -1._dp/7._dp, 3._dp/28._dp, -1._dp/12._dp, 1._dp/15._dp, -3._dp/55._dp &
    & , 1._dp/22._dp, -1._dp/26._dp, 3._dp/91._dp, -1._dp/35._dp               &
    & , 0.025_dp, -3._dp/136._dp, 1._dp/51._dp, -1._dp/57._dp /)
  Integer :: i1


  If ((m32.eq.0._dp).and.(m42.eq.0._dp)) then

   If (m12.eq.m22) then
    D27_Bagger = - 1._dp / m12
   Else
    Dm = m12 / m22 - 1._dp
    If (Abs(Dm).Lt.1.e-4_dp) then
     D27_Bagger = 1._dp / 12._dp ! n=12 
     Do i1=11,1,-1
      D27_Bagger = D27_Bagger * Dm + (-1._dp)**i1 / Real(i1,dp)
     End do
     D27_Bagger = D27_Bagger / m22
    Else
     D27_Bagger = Log(m22/m12) / (m12 - m22)
    End If
   End If

  Else If (m32.eq.0._dp) then

   If ((m12.eq.m22).and.(m12.eq.m42)) then
    D27_Bagger = -1._dp/(2*m42)
   Else If (m12.eq.m22) then
    Dm = m22 / m42 - 1._dp
    If (Abs(Dm).Lt.1.e-4_dp) then
     D27_Bagger = 1._dp / 78._dp ! n=12 -> n(n+1)/2
     Do i1=11,1,-1
      D27_Bagger = D27_Bagger * Dm &
               & + (-1._dp)**i1 * 2._dp / Real(i1*(i1+1),dp)
     End do
     D27_Bagger = D27_Bagger / (2._dp * m22)
    Else
     D27_Bagger = (m42 - m22 - m42*Log(m42/m22))/ (m22 - m42)**2
    End If

   Else If (m12.eq.m42) then
    Dm = m42 / m22 - 1._dp
    If (Abs(Dm).Lt.1.e-4_dp) then
     D27_Bagger = 1._dp / 78._dp ! n=12 -> n(n+1)/2
     Do i1=11,1,-1
      D27_Bagger = D27_Bagger * Dm &
               & + (-1._dp)**i1 * 2._dp / Real(i1*(i1+1),dp)
     End do
     D27_Bagger = D27_Bagger / (2._dp * m42)
    Else
     D27_Bagger = (m22 - m42 + m22*Log(m42/m22))/ (m22 - m42)**2
    End If

   Else If (m22.eq.m42) then
    Dm = m42 / m12 - 1._dp
    If (Abs(Dm).Lt.1.e-4_dp) then
     D27_Bagger = 1._dp / 78._dp ! n=12 -> n(n+1)/2
     Do i1=11,1,-1
      D27_Bagger = D27_Bagger * Dm &
               & + (-1._dp)**i1 * 2._dp / Real(i1*(i1+1),dp)
     End do
     D27_Bagger = D27_Bagger / (2._dp * m42)
    Else
     D27_Bagger = (m12 - m42 + m12*Log(m42/m12))/ (m12 - m42)**2
    End If

   Else
    D27_Bagger = ( m12*(-m22 + m42)*Log(m42/m12)   &
            &    + m22*(m12 - m42)*Log(m42/m22)  ) &
            &  / ((m12 - m22)*(m12 - m42)*(m42-m22))
   End If

  Else If (m42.eq.0._dp) then


   If ((m12.eq.m22).and.(m12.eq.m32)) then
    D27_Bagger = -1._dp/(2*m32)
   Else If (m12.eq.m22) then
    Dm = m22 / m32 - 1._dp
    If (Abs(Dm).Lt.1.e-4_dp) then
     D27_Bagger = 1._dp / 78._dp ! n=12 -> n(n+1)/2
     Do i1=11,1,-1
      D27_Bagger = D27_Bagger * Dm &
               & + (-1._dp)**i1 * 2._dp / Real(i1*(i1+1),dp)
     End do
     D27_Bagger = D27_Bagger / (2._dp * m22)
    Else
     D27_Bagger = (m32 - m22 - m32*Log(m32/m22))/ (m22 - m32)**2
    End If

   Else If (m12.eq.m32) then
    Dm = m32 / m22 - 1._dp
    If (Abs(Dm).Lt.1.e-4_dp) then
     D27_Bagger = 1._dp / 78._dp ! n=12 -> n(n+1)/2
     Do i1=11,1,-1
      D27_Bagger = D27_Bagger * Dm &
               & + (-1._dp)**i1 * 2._dp / Real(i1*(i1+1),dp)
     End do
     D27_Bagger = D27_Bagger / (2._dp * m32)
    Else
     D27_Bagger = (m22 - m32 + m22*Log(m32/m22))/ (m22 - m32)**2
    End If

   Else If (m22.eq.m32) then
    Dm = m32 / m12 - 1._dp
    If (Abs(Dm).Lt.1.e-4_dp) then
     D27_Bagger = 1._dp / 78._dp ! n=12 -> n(n+1)/2
     Do i1=11,1,-1
      D27_Bagger = D27_Bagger * Dm &
               & + (-1._dp)**i1 * 2._dp / Real(i1*(i1+1),dp)
     End do
     D27_Bagger = D27_Bagger / (2._dp * m32)
    Else
     D27_Bagger = (m12 - m32 + m12*Log(m32/m12))/ (m12 - m32)**2
    End If

   Else
    D27_Bagger = ( m12*(-m22 + m32)*Log(m32/m12)   &
            &    + m22*(m12 - m32)*Log(m32/m22)  ) &
            &  / ((m12 - m22)*(m12 - m32)*(m32-m22))
   End If


  Else If ((m12.Eq.m22).And.(m32.Eq.m42)) Then
   
   Dm = Abs(m12/m42 - 1._dp)
   If (m12.Eq.m42) Then
    D27_Bagger = - 1._dp / (3._dp * m12)

   Else If (Dm.Lt.1.e-2_dp) Then
    Dm = (m42-m12) / m12
    D27_Bagger = s1(16) * Dm
    Do i1=15,1,-1
     D27_Bagger = (D27_Bagger + s1(i1)) * Dm
    End Do
    D27_Bagger = (D27_Bagger + s1(0)) / (3._dp * m12)
   Else
    D27_Bagger = ( m42**2 -m12**2 - 2._dp * m12 * m42 * Log(m42/m12) )  &
           & / (m12-m42)**3
   End If

  Else If ((m12.Eq.m22).And.(m12.Eq.m42)) Then

   If ( (Abs((m12-m32)/m32)).Lt.1.e-4_dp) Then
    x = m32/m12 - 1._dp
    If (ord.Gt.11) Then
     Write(ErrCan,*) "Warning from D27_Bagger: for m22=m42=m32, m12=m22(1+x)"
     Write(ErrCan,*) "the ci coefficients are only included up to order 10"
    End If
    ci(1) = 0.25_dp ! for x, not x^0
    ci(2) = - 0.1_dp
    ci(3) =  0.05_dp
    ci(4) = -1._dp/35._dp
    ci(5) =  1._dp/56._dp
    ci(6) = -1._dp/84._dp
    ci(7) =  1._dp/120._dp
    ci(8) = -1._dp/156._dp
    ci(9) =  1._dp/220._dp

    D27_Bagger = - 1._dp/286._dp

    Do i1=9,1,-1
     D27_Bagger = D27_Bagger * x + ci(i1)
    End Do
    D27_Bagger =  (D27_Bagger * x - 1._dp) / (3._dp*m12)
 
   Else
    D27_Bagger = (3._dp*m32**2 - 4._dp*m32*m42 + m42**2               &
             & - 2*m32**2*Log(m32/m42))/(2._dp*(m32 - m42)**3)
   End If

  Else If ((m22.Eq.m32).And.(m22.Eq.m42)) Then

   If ( (Abs((m12-m22)/m22)).Lt.1.e-4_dp) Then
    x = m12/m22 - 1._dp
    If (ord.Gt.11) Then
     Write(ErrCan,*) "Warning from D27_Bagger: for m22=m42=m32, m12=m22(1+x)"
     Write(ErrCan,*) "the ci coefficients are only included up to order 10"
    End If
    ci(1) = 0.25_dp ! for x, not x^0
    ci(2) = - 0.1_dp
    ci(3) =  0.05_dp
    ci(4) = -1._dp/35._dp
    ci(5) =  1._dp/56._dp
    ci(6) = -1._dp/84._dp
    ci(7) =  1._dp/120._dp
    ci(8) = -1._dp/156._dp
    ci(9) =  1._dp/220._dp

    D27_Bagger = - 1._dp/286._dp

    Do i1=9,1,-1
     D27_Bagger = D27_Bagger * x + ci(i1)
    End Do
    D27_Bagger =  (D27_Bagger * x - 1._dp) / (3._dp*m12)
 
   Else
    D27_Bagger = ( 3._dp*m12**2 - 4._dp*m12*m22 + m22**2               &
             &   + 2*m12**2*Log(m22/m12))/(2._dp*(m12 - m22)**3)
   End If


  Else If ((m12.Eq.m22).And.(m12.Eq.m32)) Then

   If ( (Abs((m32-m42)/m32)).Lt.1.e-4_dp) Then
    x = m42/m12 - 1._dp
    If (ord.Gt.11) Then
     Write(ErrCan,*) "Warning from D27_Bagger: for m12=m42=m32, m22=m12(1+x)"
     Write(ErrCan,*) "the ci coefficients are only included up to order 10"
    End If
    ci(1) = 0.25_dp ! for x, not x^0
    ci(2) = - 0.1_dp
    ci(3) =  0.05_dp
    ci(4) = -1._dp/35._dp
    ci(5) =  1._dp/56._dp
    ci(6) = -1._dp/84._dp
    ci(7) =  1._dp/120._dp
    ci(8) = -1._dp/156._dp
    ci(9) =  1._dp/220._dp

    D27_Bagger = - 1._dp/286._dp

    Do i1=9,1,-1
     D27_Bagger = D27_Bagger * x + ci(i1)
    End Do
    D27_Bagger =  (D27_Bagger * x - 1._dp) / (3._dp*m12)
 
   Else
    D27_Bagger = -(m32**2 - 4._dp*m32*m42 + 3._dp*m42**2           &
            & - 2._dp*m42**2*Log(m42/m32))/(2._dp*(m32 - m42)**3)
   End If


  Else If ((m12.Eq.m32).And.(m12.Eq.m42)) Then

   If ( (Abs((m12-m22)/m22)).Lt.1.e-4_dp) Then
    x = m22/m12 - 1._dp
    If (ord.Gt.11) Then
     Write(ErrCan,*) "Warning from D27_Bagger: for m12=m42=m32, m22=m12(1+x)"
     Write(ErrCan,*) "the ci coefficients are only included up to order 10"
    End If
    ci(1) = 0.25_dp ! for x, not x^0
    ci(2) = - 0.1_dp
    ci(3) =  0.05_dp
    ci(4) = -1._dp/35._dp
    ci(5) =  1._dp/56._dp
    ci(6) = -1._dp/84._dp
    ci(7) =  1._dp/120._dp
    ci(8) = -1._dp/156._dp
    ci(9) =  1._dp/220._dp

    D27_Bagger = - 1._dp/286._dp

    Do i1=9,1,-1
     D27_Bagger = D27_Bagger * x + ci(i1)
    End Do
    D27_Bagger =  (D27_Bagger * x - 1._dp) / (3._dp*m12)
 
   Else
    D27_Bagger = -(m12**2 - 4._dp*m12*m22 + 3._dp*m22**2           &
             &    + 2._dp*m22**2*Log(m12/m22))/(2._dp*(m12 - m22)**3)
   End If

  Else If (m12.Eq.m22) Then

   If ( (Abs((m42-m32)/m42)).Lt.1.e-5_dp) Then
    x = m42 / m32 - 1._dp
    diff31 = m32 - m12

    If (ord.Gt.9) Then
     Write(ErrCan,*) "Warning from D27_Bagger: for m12=m22, m42=m32(1+x)"
     Write(ErrCan,*) "the ci coefficients are only included up to order 9"
    End If

    di(1) = 1._dp
    Do i1=1,ord-1
     di(i1+1) = -0.5_dp * (-m32/diff31)**(i1-1) * (i1*m12 + 2*m32) /diff31
    End Do
    di = 2._dp * di * (m12*m32) * Log(m32/m12) / diff31**3
    ci = 0
    ci(1) = -((m12 + m32)/diff31**2)
    ci(2) = (m32*(5*m12 + m32))/(2.*diff31**3)
    ci(3) = -(m32*(m12**2 + 10*m12*m32 + m32**2))/(3.*diff31**4)
    ci(4) = -(m32*(m12**3 - 11*m12**2*m32 - 47*m12*m32**2 - 3*m32**3)) &
          &  / (12.*diff31**5)
    ci(5) = -(m32*(m12**4 - 9*m12**3*m32 + 51*m12**2*m32**2 + 131*m12*m32**3 &
          &  + 6*m32**4))/(30.*diff31**6)
    ci(6) = -(m32*(m12**5 - 9*m12**4*m32 + 41*m12**3*m32**2 - 159*m12**2*m32**3 &
          &   - 284*m12*m32**4 - 10*m32**5))/(60.*diff31**7)
    ci(7) = -(m32*(2*m12**6 - 19*m12**5*m32 + 86*m12**4*m32**2           &
          &   - 264*m12**3*m32**3 + 786*m12**2*m32**4 + 1059*m12*m32**5  &
          &   + 30*m32**6))/(210.*diff31**8)
    ci(8) = ( m32*(-5*m12**7 + 51*m12**6*m32 - 243*m12**5*m32**2          &
          & + 737*m12**4*m32**3 - 1713*m12**3*m32**4 + 4167*m12**2*m32**5 &
          & + 4461*m12*m32**6 + 105*m32**7))/(840.*diff31**9)
    ci(9) = -(m32*(5*m12**8 - 55*m12**7*m32 + 281*m12**6*m32**2            &
          &  - 895*m12**5*m32**3 + 2045*m12**4*m32**4 - 3835*m12**3*m32**5 &
          &  + 7925*m12**2*m32**6 + 6989*m12*m32**7 + 140*m32**8))         &
          & / (1260.*diff31**10)

    D27_Bagger = ci(ord) + di(ord)
    Do i1=ord-1,1,-1
     D27_Bagger = (ci(i1+1) + di(i1+1) ) * x + ci(i1) + di(i1)
    End Do

   Else
    D27_Bagger = ( m42**2 * Log(m42/m12) / (m12-m42)**2               &
               & - m32**2 * Log(m32/m12) / (m12-m32)**2 ) / (m32-m42) &
               &  - m12 / ( (m12-m32) * (m12-m42) )
   End If


  Else If (m32.Eq.m42) Then

   If ((Abs(m12-m22)/m22).Lt.1.e-5_dp) Then
    x = m22 / m12 - 1._dp
    diff31 = m32 - m12 

    ci = 0._dp
    If (ord.Gt.9) Then
     Write(ErrCan,*) "Warning from D27_Bagger: for m32=m42, m22=m12(1+x)"
     Write(ErrCan,*) "the ci coefficients are only included up to order 9"
    End If
   
    ci(1) = -m12 - m32
    ci(2) = -0.5_dp * (m12*(m12 + 5*m32))
    ci(3) = -(m12*(m12**2 + 10*m12*m32 + m32**2))/3._dp
    ci(4) = (m12*(-3*m12**3 - 47*m12**2*m32 - 11*m12*m32**2 + m32**3))/12._dp
    ci(5) = -(m12*(6*m12**4 + 131*m12**3*m32 + 51*m12**2*m32**2 - 9*m12*m32**3 &
          &  + m32**4))/30._dp
    ci(6) = ( m12*(-10*m12**5 - 284*m12**4*m32 - 159*m12**3*m32**2  &
          & + 41*m12**2*m32**3 - 9*m12*m32**4 + m32**5))/60._dp
    ci(7) = -(m12*(30*m12**6 + 1059*m12**5*m32 + 786*m12**4*m32**2   &
          &   - 264*m12**3*m32**3 + 86*m12**2*m32**4 - 19*m12*m32**5 &
          &   + 2*m32**6))/210._dp
    ci(8) = -(m12*(105*m12**7 + 4461*m12**6*m32 + 4167*m12**5*m32**2      &
          &  - 1713*m12**4*m32**3 + 737*m12**3*m32**4 - 243*m12**2*m32**5 &
          &  + 51*m12*m32**6 - 5*m32**7))/840._dp
    ci(9) = -(m12*(140*m12**8 + 6989*m12**7*m32 + 7925*m12**6*m32**2       &
          &  - 3835*m12**5*m32**3 + 2045*m12**4*m32**4 - 895*m12**3*m32**5 &
          &  + 281*m12**2*m32**6 - 55*m12*m32**7 + 5*m32**8))/1260._dp

    Do i1=1,ord
     di(i1) = (m12/diff31)**(i1-1) * m32 * (2*m12 + (i1-1)*m32) / diff31**3
     ci(i1) = ci(i1) / diff31**(i1+1)
    End Do

    di = di * Log(m32/m12)

    D27_Bagger = ci(ord) + di(ord)
    Do i1=ord-1,1,-1
     D27_Bagger = D27_Bagger * x + ci(i1) + di(i1)
    End Do

   Else
    D27_Bagger = ( m12**2 * Log(m42/m12) / (m12-m42)**2               &
               & - m22**2 * Log(m42/m22) / (m22-m42)**2 ) / (m12-m22) &
               & - m42 / ( (m12-m42) * (m22-m42) )

   End If


  Else If (m12.Eq.0._dp) Then

   If ((m22.Eq.m32).And.(m32.Eq.m42)) Then
    D27_Bagger = - 0.5_dp / m42
   Else If (m22.Eq.m32) Then
    D27_Bagger = (-m32 + m42 - m42*Log(m42/m32))/(m32 - m42)**2
   Else If (m22.Eq.m42) Then
    D27_Bagger = (m32 - m42 - m32*Log(m32/m42))/(m32 - m42)**2
   Else If (m32.Eq.m42) Then
    D27_Bagger = (m22 - m42 - m22*Log(m22/m42))/(m22 - m42)**2
   Else
    D27_Bagger =  (m32*(m22 - m42)*Log(m32/m22) +          &
              &    (-m22 + m32)*m42*Log(m42/m22))/        &
              &      ((m22 - m32)*(m22 - m42)*(m32 - m42))
   End If

  Else If ((Abs(m12-m22)/m12).Lt.1.e-5_dp) Then

   x = m22 / m12 - 1._dp
   diff31 = m32 - m12 
   diff41 = m42 - m12 
   ci(1) = -1._dp + Log(m32/m12) * (1._dp + m12/diff31)
   di(1) =  1._dp - Log(m42/m12) * (1._dp + m12/diff41)
   Do i1=1,ord-1
    ci(i1+1) = ci(i1) * m12/diff31 + (-1)**i1 / Real(i1 * (i1 + 1),dp)
    di(i1+1) = di(i1) * m12/diff41 - (-1)**i1 / Real(i1 * (i1 + 1),dp)
   End Do
   ci = ci * m32 / ( (m42-m32) * diff31)
   di = di * m42 / ( (m42-m32) * diff41)
   D27_Bagger = ci(ord) + di(ord)
   Do i1=ord-1,1,-1
    D27_Bagger = (ci(i1+1) + di(i1+1) ) * x + ci(i1) + di(i1)
   End Do


  Else If ((Abs(m42-m32)/m42).Lt.1.e-5_dp) Then
   x = m42 / m32 - 1._dp

   diff31 = m32 - m12 
   diff32 = m32 - m22 

   ci(1) = -1._dp - Log(m32/m12) * (1._dp - m32/diff31)
   di(1) = 1._dp + Log(m32/m22) * (1._dp - m32/diff32)

   Do i1=1,ord-1
    ci(i1+1) = (-1)**i1/Real(i1*(1 + i1),dp) - (m32*ci(i1))/diff31
    di(i1+1) = -(-1)**(i1)/Real(i1*(1 + i1),dp) - (m32*di(i1))/diff32
   End Do  
   ci = ci * m12 / ( (m12-m22) * diff31)
   di = di * m22 / ( (m12-m22) * diff32)

   D27_Bagger = ci(ord) + di(ord)
   Do i1=ord-1,1,-1
    D27_Bagger = (ci(i1+1) + di(i1+1) ) * x + ci(i1) + di(i1)
   End Do

  Else
   D27_Bagger = ( m12 * C0_3m(m12,m32,m42)     &
              & - m22 * C0_3m(m22,m32,m42) ) / (m12 - m22)
  End If

  D27_Bagger = 0.25_dp * D27_Bagger

 End Function D27_Bagger

 Complex(dp) Function DB0(xp, xma, xmb)
 Implicit None
  Real(dp), Intent(in) :: xp, xma, xmb

  Integer :: i_max, i1, jsign
  Real(dp) :: pdb0p, dm, dmamb, dmap, dmbp, xm, dmp, x, ax, s, sumI, slam, xx &
     & , y, s2, s2a, xlam, xm1, xm2, dm1m2, dm1p, dm2p, s1, s1a, s1b, s1c     &
     & , s1d, s1e, s1f, s1p, s2b, s2p, xmax, betm2n, alph1, alpha, d1, d2, a  &
     & , beta, xnoe, b, c, d, xlogmm, h, s3, diff
  Complex(dp) :: cdb0p

  dm = (Sqrt(xma) - Sqrt(xmb))**2
  dmamb = xma - xmb
  dmap = xma - xp
  dmbp = xmb - xp

  If((xp.Eq.dm).And.(xp.Ne.0).And.(xma.Ne.0).And.(xmb.Ne.0) ) Then
   pdb0p = (xmb - xma)/2._dp/dm*Log(xmb/xma) - 2._dp
   DB0 = pdb0p/dm

  Else If ((xma.Eq.0._dp).And.(xmb.Eq.0._dp)) Then ! both masses = 0
   DB0 = -1._dp / xp

  Else If ((xma.Eq.0._dp).Or.(xmb.Eq.0._dp)) Then ! one masses = 0
   If (xma.Eq.0._dp) Then
    xm = xmb
    dmp = dmbp
   Else
    xm = xma
    dmp = dmap
   End If

   If (xp.Eq.0 ) Then
    DB0 = 1._dp/(2._dp*xm)

   Elseif (dmp.Eq.0._dp) Then

    If ( lambda2 .Eq. 0._dp ) Then
     Call WriteLFerror(6)
     cdb0p = 0
    Else
     cdb0p = -1._dp + 0.5_dp * Log(xm/lambda2)
    Endif
    DB0 = cdb0p / xp

   Else

    x = xp/xm
    ax = Abs(x)

    If ( ax .Lt. xloss ) Then
     i_max = FindBound(x, 1, DB0serie1)
     sumI = 0._dp     
     Do i1=i_max,3,-1
      sumI = x * (DB0serie1(i1) + sumI)
     End Do
     cdb0p = x*(DB0serie1(2) + sumI)

    Else
     s = Log(Abs(dmp/xm))
     cdb0p = -(1 + s*xm/xp)
     If ( xp.Gt.xm ) cdb0p = cdb0p + Cmplx(0._dp,xm/xp*pi,dp)
    Endif
    DB0 = cdb0p / xp
   Endif

  Else If (dmamb.Eq.0._dp) Then ! both masses equal

   xm = xma
   dmp = dmap

   If ( Abs(xp) .Lt. 8*xloss*xm ) Then
    ! a Taylor expansion seems appropriate as the result will go
    ! as k^2 but seems to go as 1/k !!

    x = -xp/xm
    ax = Abs(x)
    i_max = FindBound(x, 1, DB0serie2)
    sumI = 0._dp     
    Do i1=i_max,1,-1
     sumI = x * (DB0serie2(i1) + sumI)
    End Do
    cdb0p = - sumI

    If ( xp.Ne.0 ) Then
     DB0 = cdb0p*(1._dp/xp)
    Else
     DB0 = DB0serie2(1)/xm
    Endif

   Else ! normal case

    xlam = Kappa2p(-xp,-xm,-xm,dmp,dmp,0._dp)
    If ( xlam .Eq. 0 ) Then
     Call WriteLFerror(7)
     Return
    Elseif ( xlam .Gt. 0 ) Then
     slam = Sqrt(xlam)
     s2a = dmp + xm
     s2 = s2a + slam
     If ( Abs(s2) .Gt. xloss*slam ) Then
      jsign = 1
     Else
      s2 = s2a - slam
      jsign = -1
     Endif
     ax = Abs(s2/(2*xm))
     If ( ax .Lt. xalogm ) Then
      s = 0._dp
     Elseif( (ax-1.Lt.0.1_dp).And.(s2.Gt.0._dp) ) Then
     !  In this case a quicker and more accurate way is to calculate log(1-x).
      s2 = (xp - slam)
      s = 2*xm/slam*Log1minusX(s2/(2*xm))
     Else
      s = 2._dp*xm/slam*Log(ax)
      If ( jsign .Eq. -1 ) s = -s
     Endif
     If (xp.Gt.(2._dp*xm) ) Then ! existence of a non-zero imaginary part
      y = pi*2*xm/slam
     Else
      y = 0
     Endif
    Else ! the root is complex (k^2 between 0 and (2*m1)^2)
     slam = Sqrt(-xlam)
     s = 4._dp*xm/slam*Atan2(xp,slam)
     y = 0._dp
    Endif
    xx = s - 1._dp
    cdb0p = Cmplx(xx,y,dp)
    DB0 = cdb0p*(1._dp/xp)
   Endif

  Else    ! general  case

   If (xma.Gt.xmb ) Then
    xm2 = xma
    xm1 = xmb
    dm1m2 = -dmamb
    dm1p = dmbp
    dm2p = dmap
   Else
    xm1 = xma
    xm2 = xmb
    dm1m2 = dmamb
    dm1p = dmap
    dm2p = dmbp
   Endif

   x = xm2/xm1
   If ( 1._dp .Lt. xalogm*x ) Then
    Call WriteLFerror(8)
    xlogmm = 0._dp
   Elseif ( Abs(x-1._dp) .Lt. xloss ) Then
    xlogmm = Log1minusX(dm1m2/xm1)
   Else
    xlogmm = Log(x)
   Endif

   If ( xp .Eq. 0 ) Then
    s1 = xm1*xm2*xlogmm/dm1m2**3
    s2 = (xm1+xm2)/(2*dm1m2**2)
    s = s1 + s2
    If ( Abs(s) .Lt. xloss**2*s2 ) Then !  second try
     h = Log1minusXpXn(dm1m2/xm1,2)
     s1 = -xm1*h/dm1m2**2
     s2 = 1/(2*xm1)
     s3 = xm1**2*h/dm1m2**3
     s = s1 + s2 + s3
     If ( Abs(s) .Lt. xloss*Max(Abs(s2),Abs(s3)) ) Call WriteLFerror(9)
    Endif
    DB0 = s
    cdb0p = 0

   Else

    xlam = Kappa2p(-xp,-xm2,-xm1,dm2p,dm1p,dm1m2)
    diff = xlam + xp*(dm2p+xm1)
    If ( Abs(diff) .Lt. xloss*xlam ) Then
     h = dm1m2**2 - xp*(xm1+xm2)
     If ((Abs(h).Lt.xloss*dm1m2**2).And.(dm1m2**2.Lt.Abs(xlam))) diff = h
    Endif

    If ( xlam .Eq. 0 ) Then
     Call WriteLFerror(7)
     Return

    Elseif ( xlam .Gt. 0 ) Then
     slam = Sqrt(xlam)
     s2a = dm2p + xm1
     s2 = s2a + slam
     If ( Abs(s2) .Gt. xloss*slam ) Then
      jsign = 1
     Else
      s2 = s2a - slam
      jsign = -1
     Endif
     s2 = s2**2/(4*xm1*xm2)
     If ( Abs(s2) .Lt. xalogm ) Then
      Call WriteLFerror(10)
      s2 = 0
     Elseif ( Abs(s2-1._dp).Lt.xloss) Then
      If ( jsign.Eq.1 ) Then
       s2 = -slam*(s2a+slam)/(2*xm1*xm2)
       s2 = -diff/(2*slam*xp)*Log1MinusX(s2)
      Else
       Call WriteLFerror(11)
       s2 = +slam*(s2a-slam)/(2*xm1*xm2)
       s2 = +diff/(2*slam*xp)*Log1minusX(s2)
      Endif
     Else
      s2 = -diff/(2*slam*xp)*Log(s2)
      If ( jsign .Eq. -1 ) s2 = -s2
     Endif
     s1 = -dm1m2*xlogmm/(2*xp)
     xx = s1+s2-1._dp
     If ( Abs(xx) .Lt. xloss**2*Max(Abs(s1),Abs(s2)) ) Then
     ! this is unacceptable, try a better solution
      s1a = diff + slam*dm1m2
      If ( Abs(s1a) .Gt. xloss*diff ) Then
       s1 = -s1a/(2*xp*slam)
      Else
       s1 = -2*xm1*xm2*xp/(slam*(diff - slam*dm1m2))
      Endif
      s = s1
      s1 = s1*xlogmm

      If ( Abs(xp) .Lt. xm2 ) Then
       s2a = xp - dm1m2
      Else
       s2a = xm2 - dm1p
      Endif
      s2 = s2a - slam
      If ( Abs(s2) .Gt. xloss*slam ) Then
       s2 = s2 / (2*xm2)
      Else
       s2 = (2*xp) / (s2a+slam)
      Endif
      If ( Abs(s2).Lt.0.1_dp) Then
       s2 = Log1minusX(s2)
      Elseif ( s2.Eq.1._dp ) Then
       Call WriteLFerror(12)
       If ( xp .Gt. xm1+xm2 ) Then
        y = -pi*diff/(slam*xp)
       Else
        y = 0
       Endif
       cdb0p = Cmplx(xx,y,dp)
       DB0 = cdb0p*(1._dp/xp)
       Return
      Else
       s2 = Log(Abs(1._dp - s2))
      Endif
      s2 = -diff/(slam*xp)*s2
      xx = s1 + s2 - 1

      If ( Abs(xx) .Lt. xloss**2*Max(Abs(s1),Abs(s2)) ) Then
       ! third try, we accept two times xloss because that's the same
       ! as in this try)
       ! A Taylor expansion might work.  We expand inside the logs
       xnoe = s2a+slam
       a = 1
       b = 2/xnoe-1/xp
       c = -4/(xp*xnoe)
       d = Sqrt((2/xnoe)**2 + 1/xp**2)
       Call roots(a,b,c,d,d1,d2)
       If ( xp.Gt.0 ) Then
        beta = d2
       Else
        beta = d1
       Endif
       alpha = beta*diff/slam
       alph1 = 1-alpha
       If ( alph1 .Lt. xloss ) Then
        s1a = 4*xp**2*xm1*xm2/(slam*dm1m2*(diff-slam*dm1m2))
        s1b = -diff/slam*4*xm1*xp/(dm1m2*xnoe*(2*xp-xnoe))
        b = -1/xp
        c = -(2/xnoe)**2
        Call roots(a,b,c,d,d1,d2)
        If ( xp.Gt.0 ) Then
         betm2n = d2
        Else
         betm2n = d1
        Endif
        d1 = s1a + s1b - diff/slam*betm2n
        xmax = Max(Abs(s1a),Abs(s1b))
        If ( xmax .Lt. 1 ) Then
         alph1 = d1
        Else
         xmax = 1
        Endif
       Else
        betm2n = beta - 2/xnoe
       Endif

       s2p = s2 - alpha
       If ( Abs(s2p) .Lt. xloss*Abs(s2) ) Then
        x = beta*xp
        i_max = FindBound(x, 1, DB0serie3)
        s2a = 0._dp     
        Do i1=i_max,4,-1
         s2a = x * (DB0serie3(i1) + s2a)
        End Do
        s2a = x**3*(DB0serie3(3)+s2a)
        s2b = 2*xp/xnoe*(s2a + x**2/2)
        s2p = s2b - s2a
        s2p = -diff/(xp*slam)*Log1minusX(s2p)
       Endif

       s1p = s1 - alph1
       If ( Abs(s1p) .Lt. xloss*Abs(s1) ) Then
        x = slam*(diff-slam*dm1m2)*alph1/(2*xp*xm1*xm2)
        h = (2*xp*(xm1+xm2) - xp**2)/(slam-dm1m2)
        ax = Abs(x)

        s1b = diff*(diff-slam*dm1m2)*betm2n/(2*xp*xm1*xm2)
        s1c = 1/(xm1*xnoe*(2*xp-xnoe))*(                      &
            &     xp*( 4*xp*xm2 + 2*dm1m2**2/xm2*(xp-h) +     &
            &     2*dm1m2*(3*xp-h) - 8*dm1m2**2 )             &
            &     - 2*dm1m2**3/xm2*(3*xp-h)                   &
            &     + 4*dm1m2**4/xm2 )
        s1d = x*dm1m2/xm1
        s1e = -x**2/2

        i_max = FindBound(x, 1, DB0serie3)
        s1a = 0._dp     
        Do i1=i_max,4,-1
         s1a = x * (DB0serie3(i1) + s1a)
        End Do
        s1a = -x**3 *(DB0serie3(3) + s1a)
        s1f = dm1m2/xm1*(x**2/2 - s1a)
        s1p = s1e + s1d + s1c + s1b + s1a + s1f
        xmax = Max(Abs(s1a),Abs(s1b),Abs(s1c),Abs(s1d),Abs(s1e))
        s1p = s*Log1MinusX(s1p)
       Endif

       xx = s1p + s2p
      Endif ! third try
     Endif  ! second try

     If ( xp .Gt. xm1+xm2 ) Then  ! imaginary part
      y = -pi*diff/(slam*xp)
     Else
      y = 0
     Endif

    Else ! the root is complex (k^2 between -(m1+m2)^2 and -(m2-m1)^2)
     slam = Sqrt(-xlam)
     xnoe = dm2p + xm1
     s1 = -(dm1m2/(2*xp))*xlogmm
     s2 = -diff/(slam*xp)*Atan2(slam,xnoe)
     xx = s1 + s2 - 1._dp

     y = 0
    Endif

    cdb0p = Cmplx(xx,y,dp)
    DB0 = cdb0p*(1._dp/xp)
   Endif
  End If

 End Function DB0

 Real(dp) Function E_t_0(x)
 !-----------------------------------------------------------
 ! function E^t_0(x) taken from C.Bobeth et al., NPB574 (2000) 291
 !                              ,  hep-ph/9904413
 ! written by Werner Porod, 02.12.2005
 !----------------------------------------------------------------------------
 Implicit None 
  Real(dp), Intent(in) :: x

   E_t_0 = (-9._dp * x**2 + 16._dp * x - 4._dp) * Log(x)         &
       &     / (6._dp * (1._dp-x)**4)                            &
       & + (-7._dp * x**3 - 21._dp * x**2 + 42._dp * x + 4._dp)  &
       &     / (36._dp * (1._dp-x)**3)

 End Function E_t_0



 Subroutine ff2d22(dl2d22,xpi,dpipj,piDpj, i, j,k,kj,iskj, m,n,nm,isnm)
 !--------------------------------------------------------------------
 !							
 !	Calculate					
 !							
 !	/ si mu	 mu nu \2				
 !	|d	d      |				
 !	\ sj sk	 sm sn /				
 !							
 !	=   si.sj^2*sk.sm^2*sn.sn			
 !	- 2*si.sj^2*sk.sm*sk.sn*sm.sn			
 !	+   si.sj^2*sk.sn^2*sm.sm			
 !	- 2*si.sj*si.sk*sj.sm*sk.sm*sn.sn		
 !	+ 2*si.sj*si.sk*sj.sm*sk.sn*sm.sn		
 !	+ 2*si.sj*si.sk*sj.sn*sk.sm*sm.sn		
 !	- 2*si.sj*si.sk*sj.sn*sk.sn*sm.sm		
 !	+   si.sk^2*sj.sm^2*sn.sn			
 !	- 2*si.sk^2*sj.sm*sj.sn*sm.sn			
 !	+   si.sk^2*sj.sn^2*sm.sm			
 !							
 !	Input:	xpi(ns)			as usual	
 !		dpipj(ns,ns)		  -"-		
 !		piDpj(ns,ns)		  -"-		
 !		i,j,k,kj,iskj		see above	
 !		m,n,nm,isnm		 -"-		
 !							
 !	Output:	dl2d22			see above	
 !							
 !-------------------------------------------------------------------- 
 Implicit None
  Integer :: i,j,k,kj,iskj,m,n,nm,isnm
  Real(dp) :: dl2d22,xpi(10),dpipj(10,10),piDpj(10,10)

  Integer :: ii,isii
  Real(dp) :: s(10),del2s,del23,del24,del27,som,sMax,xmax

 ! special cases:
  If ( i == n .Or. i == m ) Then
      Call ffdl2s(del2s,piDpj, j,k,kj,iskj, m,n,nm,isnm, 10)
    dl2d22 = xpi(i)*del2s**2
    Return
  End If

 ! calculations:
 !	We use the product form
  If ( i == 3 ) Then
    del23 = 0
  Else If ( i <= 4 ) Then
    ii = inx(3,i)
    isii = isgn(3,i)
      Call ffdl2s(del23,piDpj,i,3,ii,isii,j,k,kj,iskj,10)
  Else
      Call ffdl2t(del23,piDpj,i,3,j,k,kj,iskj,+1,10)
  End If
  If ( i == 4 ) Then
    del24 = 0
  Else If ( i <= 4 ) Then
    ii = inx(n,i)
    isii = isgn(n,i)
      Call ffdl2s(del24,piDpj,i,4,ii,isii,j,k,kj,iskj,10)
  Else
      Call ffdl2t(del24,piDpj,i,4,j,k,kj,iskj,+1,10)
  End If
  s(1) = xpi(4)*del23**2
  s(2) = -2*piDpj(3,4)*del23*del24
  s(3) = xpi(3)*del24**2
  Dl2d22 = s(1) + s(2) + s(3)
  sMax = Max(Abs(s(1)),Abs(s(2)),Abs(s(3)))
  If ( Abs(dl2d22) >= xloss*sMax ) Return
  som = dl2d22
  xMax = smax
 !	try the special case k=4 (for use in ee->mumu among others)
  If ( i < 4 .And. k == 4 .And. Abs(s(3)) < xloss*sMax&
    &.And. ( Abs(dpipj(i,inx(4,i))) < xloss*xpi(i) .Or.&
    &Abs(piDpj(j,inx(4,i))) < xloss*Abs(piDpj(j,4)) ) ) Then
    s(1) = -del23*piDpj(i,4)*piDpj(j,3)*xpi(4)
    s(2) =  del23*dpipj(i,inx(4,i))*piDpj(j,4)*piDpj(3,4)
    s(4) =  del23*piDpj(3,4)*xpi(4)*piDpj(j,inx(4,i))*isgn(4,i)
    dl2d22 = s(1) + s(2) + s(3) + s(4)
    sMax = Max(Abs(s(1)),Abs(s(2)),Abs(s(3)),Abs(s(4)))
      If ( Abs(dl2d22) >= xloss*sMax ) Return
      If ( sMax < xmax ) Then
      som = dl2d22
      xMax = smax
      End If
  End If
  Call ffdl2t(del27,piDpj,i,7,j,k,kj,iskj,+1,10)
  s(1) = xpi(7)*del24**2
  s(2) = -2*piDpj(4,7)*del24*del27
  s(3) = xpi(4)*del27**2
  Dl2d22 = s(1) + s(2) + s(3)
  sMax = Max(Abs(s(1)),Abs(s(2)),Abs(s(3)))
  If ( Abs(dl2d22) >= xloss*sMax ) Return
  If ( sMax < xmax ) Then
    som = dl2d22
    xMax = smax
  End If
  s(1) = xpi(7)*del23**2
  s(2) = -2*piDpj(3,7)*del23*del27
  s(3) = xpi(3)*del27**2
  Dl2d22 = s(1) + s(2) + s(3)
  sMax = Max(Abs(s(1)),Abs(s(2)),Abs(s(3)))
  If ( Abs(dl2d22) >= xloss*sMax ) Return
 !
 !	We'll have to think of something more intelligent ...
 !
  If ( sMax < xmax ) Then
    som = dl2d22
    xMax = smax
  End If
  Dl2d22 = som
 
 End Subroutine ff2d22

 Subroutine ff2dl2(del2d2,del2n,xpi,dpipj,piDpj, i, &
                 & j,k,kj,iskj,l, m,n,nm,isnm, ns)
 !--------------------------------------------------------------------
 !							
 !	Calculate					
 !							
 !	 si mu	 mu sl					
 !	d	d	= si.sj*sk.sm*sl.sn - si.sk*sj.sm*sl.sn
 !	 sj sk	 sm sn		- si.sj*sk.sn*sl.sm + si.sk*sj.sn*sl.sm
 !							
 !	with p(kj) = iskj*(sk-sj)			
 !	with p(nm) = isnm*(sn-sm)			
 !							
 !	Input:	xpi(ns)			as usual	
 !		dpipj(ns,ns)		  -"-		
 !		piDpj(ns,ns)		  -"-		
 !		i,j,k,kj,iskj		see above	
 !		l,m,n,nm,isnm		  -"-		
 !							
 !	Output:	del2d2			see above	
 !		del2n			it is needed in fftran anyway
 !							
 !-------------------------------------------------------------------- 
 Implicit None
  Integer :: i,j,k,kj,iskj,l,m,n,nm,isnm,ns
  Real(dp) :: del2d2,del2n,xpi(10),dpipj(10,10),piDpj(10,10)

  Integer :: isii,ii,ik,ij,im,in
  Real(dp) :: s(5),del2m,del2nm,som,xMax,smax

 ! get del2n:	we need this in any case !
  If ( i == n ) Then
    del2n = 0
  Else If ( i <= 4 ) Then
    ii = inx(n,i)
    isii = isgn(n,i)
      Call ffdl2s(del2n,piDpj,i,n,ii,isii,j,k,kj,iskj,10)
  Else
      Call ffdl2t(del2n,piDpj,i,n,j,k,kj,iskj,+1,10)
  End If

 !  #[ special cases:
  If ( i == l .And. j == m .And. k == n ) Then
      Call ffdl3m(s,.False.,0._dp,0._dp,xpi,dpipj,piDpj,ns,j,k,kj, i,1)
    del2d2 = -s(1)
    Return
  End If
  If ( k == l .And. j <= 4 ) Then
      Call ffdl2s(del2m,piDpj, j,l,inx(l,j),isgn(l,j), m,n,nm,isnm, 10)
    del2d2 = -piDpj(i,k)*del2m
    Return
  End If

 !  #[ calculations:
  If ( i == m ) Then
    del2m = 0
  Else If ( i <= 4 ) Then
    ii = inx(m,i)
    isii = isgn(m,i)
      Call ffdl2s(del2m,piDpj,i,m,ii,isii,j,k,kj,iskj,10)
  Else
      Call ffdl2t(del2m,piDpj,i,m,j,k,kj,iskj,+1,10)
  End If
  s(1) = del2m*piDpj(n,l)
  s(2) = del2n*piDpj(m,l)
  sMax = Abs(s(1))
  Del2d2 = s(1) - s(2)
  If ( Abs(del2d2) >= xloss*sMax ) Return
  som = del2d2
  xMax = smax
  Call ffdl2t(del2nm,piDpj,i,nm,j,k,kj,iskj,+1,10)
  s(1) = del2n*piDpj(nm,l)
  s(2) = del2nm*piDpj(n,l)
  Del2d2 = isnm*(s(1) - s(2))
  sMax = Abs(s(2))
  If ( Abs(del2d2) >= xloss*Abs(s(1)) )Return
  If ( sMax < xmax ) Then
    som = del2d2
    xMax = smax
  End If
  s(1) = del2m*piDpj(nm,l)
  s(2) = del2nm*piDpj(m,l)
  Del2d2 = isnm*(s(1) - s(2))
  sMax = Abs(s(2))
  If ( Abs(del2d2) >= xloss*Abs(s(1)) ) Return
  If ( sMax < xmax ) Then
    som = del2d2
    xMax = smax
  End If
 !	One more special case:
  If ( k == m ) Then
    isii = -1
    ik = j
    ij = k
    im = m
    in = n
  Else If ( j == m ) Then
    isii = +1
    ik = k
    ij = j
    im = m
    in = n
  Else If ( j == n ) Then
    isii = -1
    ik = k
    ij = j
    im = n
    in = m
  Else If ( k == n ) Then
    isii = +1
    ik = j
    ij = k
    im = n
    in = m
  Else
   Del2d2 = som
   Return
  End If
  If ( ij == im .And. i <= 4 .And. ij <= 4 .And. in <= 4 ) Then
      If ( inx(ij,i) > 0 .And. inx(im,l) > 0 ) Then
        If (  Abs(dpipj(i,inx(ij,i))) < xloss*Abs(xpi(ij)) &
        .And. Abs(dpipj(l,inx(im,l))) < xloss*Abs(xpi(im)) ) Then
        s(1) = piDpj(l,in)*piDpj(ik,ij)*dpipj(i,inx(ij,i))/2
        s(2) = isgn(ij,i)*piDpj(l,in)*xpi(ij)*piDpj(ik, inx(ij,i))/2
        s(3) = -piDpj(i,ij)*piDpj(ik,in)*piDpj(l,im)
        s(4) = piDpj(i,ik)*piDpj(im,in)*dpipj(l,inx(im,l))/2
        s(5) = isgn(im,l)*piDpj(i,ik)*xpi(im)*piDpj(in, inx(im,l))/2
        del2d2 = s(1) + s(2) + s(3) + s(4) + s(5)
          If ( isii < 0 ) del2d2 = -del2d2
        sMax = Max(Abs(s(1)),Abs(s(2)),Abs(s(3)),Abs(s(4)), Abs(s(5)))
          If ( Abs(del2d2) >= xloss**2*Abs(sMax) ) Return
          If ( sMax < xmax ) Then
          som = del2d2
          xMax = smax
          End If
        End If
      End If
  End If

 !	give up
  Del2d2 = som

 End Subroutine ff2dl2


 Subroutine ff3dl2(del3d2,xpi,dpipj,piDpj, i, &
                  & j,k,kj,iskj, l,m,ml,isml, n, o,p,po,ispo)
 !--------------------------------------------------------------------
 !							
 !	Calculate					
 !							
 !	 si mu	 mu nu	 mu sn				
 !	d	d	d	= ...			
 !	 sj sk	 sl sm	 so sp				
 !							
 !	with p(kj) = iskj*(sk-sj)			
 !	     p(ml) = isml*(sm-sl)			
 !	     p(po) = ispo*(sp-so)			
 !							
 !	Input:	xpi(ns)			as usual	
 !		dpipj(ns,ns)		  -"-		
 !		piDpj(ns,ns)		  -"-		
 !		i,j,k,kj,iskj		see above	
 !		l,m,ml,isml		  -"-		
 !		n,o,p,po,ispo		  -"-		
 !							
 !	Output:	del3d2			see above	
 !							
 !-------------------------------------------------------------------- 
 Implicit None
  Integer :: i,j,k,kj,iskj,l,m,ml,isml,n,o,p,po,ispo
  Real(dp) :: del3d2,xpi(10),dpipj(10,10),piDpj(10,10)

  Integer :: isii,ii
  Real(dp) :: s(2),dl2il,dl2im,dl2ln,dl2mn,dl2iml,dl2mln
  Real(dp) :: d2d2j,d2d2k,d2d2kj,dum,d2d2o,d2d2p,d2d2po
  Real(dp) :: som,xMax

 !  split up l,m:
  If ( i == l ) Then
    dl2il = 0
  Else If ( i <= 4 ) Then
    ii = inx(l,i)
    isii = isgn(l,i)
      Call ffdl2s(dl2il,piDpj,i,l,ii,isii,j,k,kj,iskj,10)
  Else
      Call ffdl2t(dl2il,piDpj,i,l,j,k,kj,iskj,+1,10)
  End If
  If ( m == n ) Then
    dl2mn = 0
  Else If ( i <= 4 ) Then
    ii = inx(n,m)
    isii = isgn(n,m)
      Call ffdl2s(dl2mn,piDpj,m,n,ii,isii,o,p,po,ispo,10)
  Else
      Call ffdl2t(dl2mn,piDpj,m,n,o,p,po,ispo,+1,10)
  End If
  s(1) = dl2il*dl2mn
  If ( i == m ) Then
    dl2im = 0
  Else If ( i <= 4 ) Then
    ii = inx(m,i)
    isii = isgn(m,i)
      Call ffdl2s(dl2im,piDpj,i,m,ii,isii,j,k,kj,iskj,10)
  Else
      Call ffdl2t(dl2im,piDpj,i,m,j,k,kj,iskj,+1,10)
  End If
  If ( l == n ) Then
    dl2ln = 0
  Else If ( i <= 4 ) Then
    ii = inx(n,l)
    isii = isgn(n,l)
      Call ffdl2s(dl2ln,piDpj,l,n,ii,isii,o,p,po,ispo,10)
  Else
      Call ffdl2t(dl2ln,piDpj,l,n,o,p,po,ispo,+1,10)
  End If
  s(2) = dl2im*dl2ln
  Del3d2 = s(1) - s(2)
  If ( Abs(del3d2) >= xloss*Abs(s(1)) ) Return
  som = del3d2
  xMax = Abs(s(1))
 !
 !	rotate l,m
 !
  Call ffdl2t(dl2mln,piDpj,ml,n,o,p,po,ispo,+1,10)
  Call ffdl2t(dl2iml,piDpj,i,ml,j,k,kj,iskj,+1,10)
  s(1) = dl2im*dl2mln
  s(2) = dl2iml*dl2mn
  Del3d2 = isml*(s(1) - s(2))
  If ( Abs(del3d2) >= xloss*Abs(s(1)) ) Return
  If ( Abs(s(1)) < xMax ) Then
    som = del3d2
    xMax = Abs(s(1))
  End If
  s(1) = dl2il*dl2mln
  s(2) = dl2iml*dl2ln
  Del3d2 = isml*(s(1) - s(2))
  If ( Abs(del3d2) >= xloss*Abs(s(1)) ) Return
  If ( Abs(s(1)) < xMax ) Then
    som = del3d2
    xMax = Abs(s(1))
  End If
 !  #] split up l,m: 
 !  #[ split up j,k:
  Call ff2dl2(d2d2k,dum,xpi,dpipj,piDpj, k, l,m,ml,isml, n, &
    o,p,po,ispo, 10)
  Call ff2dl2(d2d2j,dum,xpi,dpipj,piDpj, j, l,m,ml,isml, n, &
    o,p,po,ispo, 10)
  s(1) = piDpj(i,j)*d2d2k
  s(2) = piDpj(i,k)*d2d2j
  Del3d2 = s(1) - s(2)
  If ( Abs(del3d2) >= xloss*Abs(s(1)) ) Return
  If ( Abs(s(1)) < xMax ) Then
    som = del3d2
    xMax = Abs(s(1))
  End If
  Call ff2dl2(d2d2kj,dum,xpi,dpipj,piDpj, kj, l,m,ml,isml, n, &
    o,p,po,ispo, 10)
  s(1) = piDpj(i,k)*d2d2kj
  s(2) = piDpj(i,kj)*d2d2k
  Del3d2 = iskj*(s(1) - s(2))
  If ( Abs(del3d2) >= xloss*Abs(s(1)) ) Return
  If ( Abs(s(1)) < xMax ) Then
    som = del3d2
    xMax = Abs(s(1))
  End If
  s(1) = piDpj(i,j)*d2d2kj
  s(2) = piDpj(i,kj)*d2d2j
  Del3d2 = iskj*(s(1) - s(2))
  If ( Abs(del3d2) >= xloss*Abs(s(1)) ) Return
  If ( Abs(s(1)) < xMax ) Then
    som = del3d2
    xMax = Abs(s(1))
  End If
 !  #] split up j,k: 
 !  #[ split up o,p:
  Call ff2dl2(d2d2o,dum,xpi,dpipj,piDpj, i, j,k,kj,iskj, o, &
    l,m,ml,isml, 10)
  Call ff2dl2(d2d2p,dum,xpi,dpipj,piDpj, i, j,k,kj,iskj, p, &
    l,m,ml,isml, 10)
  s(1) = piDpj(p,n)*d2d2o
  s(2) = piDpj(o,n)*d2d2p
  Del3d2 = s(1) - s(2)
  If ( Abs(del3d2) >= xloss*Abs(s(1)) ) Return
  If ( Abs(s(1)) < xMax ) Then
    som = del3d2
    xMax = Abs(s(1))
  End If
  Call ff2dl2(d2d2po,dum,xpi,dpipj,piDpj, i, j,k,kj,iskj, po, &
    l,m,ml,isml, 10)
  s(1) = piDpj(po,n)*d2d2p
  s(2) = piDpj(p,n)*d2d2po
  Del3d2 = ispo*(s(1) - s(2))
  If ( Abs(del3d2) >= xloss*Abs(s(1)) ) Return
  If ( Abs(s(1)) < xMax ) Then
    som = del3d2
    xMax = Abs(s(1))
  End If
  s(1) = piDpj(po,n)*d2d2o
  s(2) = piDpj(o,n)*d2d2po
  Del3d2 = ispo*(s(1) - s(2))
  If ( Abs(del3d2) >= xloss*Abs(s(1)) ) Return
  If ( Abs(s(1)) < xMax ) Then
    som = del3d2
    xMax = Abs(s(1))
  End If
 !  #] split up o,p: 
 !  #[ give up:
  Del3d2 = som
 !  #] give up: 
 !###] ff3dl2: 
 End Subroutine ff3dl2

 Subroutine ffcc0p(cs3,ipi12,isoort,clogi,ilogi,cpi,cpipj, &
   &  cpiDpj,sdel2,cel2si,etalam,etami,delpsi,alpha,npoin)
 !------------------------------------------------------------------------
 !	Calculates the threepoint function closely following
 !	recipe in 't Hooft & Veltman, NP B(183) 1979.
 !	Bjorken and Drell metric is used nowadays!
 !
 !	    p2	^ |
 !		| |
 !		 / \
 !	      m2/   \m3
 !	p1     /     \	p3
 !	<-    /  m1   \ ->
 !	------------------------
 !
 !	Input:	cpi(1-3)   (complex)	pi squared (,2=untransformed
 !					when npoin=4)
 !		cpi(4-6)   (complex)	internal mass squared
 !		cpipj(6,6)   (complex)	cpi(i)-cpi(j)
 !		cpiDpj(6,6)   (complex)	pi(i).pi(j)
 !
 !	Output: cs3	 (complex)(48)	C0, not yet summed.
 !		ipi12	 (  Integer)(3)	factors pi^2/12, not yet summed
 !		cslam	 (complex)	lambda(p1,p2,p3).
 !		isoort	 (  Integer)(3)	indication of he method used
 !
 !	Calls:	ffcel2,ffcoot,ffccyz,ffcdwz,ffcs3,ffcs4
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ipi12(8),isoort(8),ilogi(3),npoin
  Complex(dp) :: cs3(80),clogi(3),cpi(6),cpipj(6,6), &
    cpiDpj(6,6),sdel2,cel2si(3),etalam,etami(6), delpsi(3),alpha(3)

  Integer :: i,j,k,ip,ierw,jsoort(8),iw,ismall(3)
  Logical :: l4,l4pos
  Complex(dp) :: c,cs,cs1,cs2,cs4,ci
  Complex(dp) :: cy(4,3),cz(4,3),cw(4,3),cdyz(2,2,3), &
    cdwy(2,2,3),cdwz(2,2,3),cd2yzz(3),cd2yww(3)
  Complex(dp) :: csdl2i(3)

 !  #[ get roots etc:
 !  #[   get z-roots:
  If ( npoin /= 3 ) Then
    l4pos = .False.
  Else
    l4pos = l4also
  End If
  Do i=1,3
 !
 !	    get roots (y,z)
 !
   ip = i+3
 !	    first get the roots
   j = i+1
   If ( j == 4 ) j = 1
   csdl2i(i) = Sqrt(-cel2si(i))
   If ( cpi(ip) == 0 ) Then
     If ( i == 1 .And. alpha(3) == 0 .Or.i == 3 .And. alpha(1) == 0 ) Then
      isoort(2*i-1) = 0
      isoort(2*i) = 0
      l4pos = .False.
      Cycle
     End If
   End If
   Call ffccyz(cy(1,i),cz(1,i),cdyz(1,1,i),cd2yzz(i),i, &
    & sdel2,csdl2i(i),etalam,etami,delpsi(i), cpi,cpiDpj,isoort(2*i-1))
  End Do
 !  #]   get z-roots:
 !  #[   get w-roots:
 !
 !	get w's:
 !
  ierw = 0
  l4 = .False.
  If ( isoort(4) == 0 ) Then
     Call WriteLFerror(19)

  Else
   Do iw = 1,3,2
    If ( .Not. l4pos .Or. alpha(4-iw) == 0 ) Then
     jsoort(2*iw-1) = 0
     jsoort(2*iw) = 0
     l4pos = .False.
    Else
     jsoort(2*iw-1) = -1
     jsoort(2*iw) = -1
     cd2yww(iw) = -cd2yzz(2)/alpha(4-iw)
      Do j=1,2
        cw(j+iw-1,iw) = cz(j+3-iw,2)/alpha(4-iw)
        cw(j+3-iw,iw) = 1 - cw(j+iw-1,iw)
        If ( absc(cw(j+3-iw,iw)) < xloss ) Then
        cs = cz(j+iw-1,2) - alpha(iw)
          If ( absc(cs) < xloss*absc(alpha(iw)) ) Then
          ierw = 1
          Cycle
          End If
        cw(j+3-iw,iw) = cs/alpha(4-iw)
        End If
        cdwy(j,2,iw) = cdyz(2,j,2)/alpha(4-iw)
        Do i=1,2
          cdwz(j,i,iw) = cw(j,iw) - cz(i,iw)
          If ( absc(cdwz(j,i,iw)) >= xloss*absc(cw(j,iw)) ) Cycle
          cdwz(j,i,iw) = cz(i+2,iw) - cw(j+2,iw)
          If ( absc(cdwz(j,i,iw)) >= xloss*absc(cw(j+2,iw)) ) Cycle
          cdwz(j,i,iw) = cdwy(j,2,iw) + cdyz(2,i,iw)
          If ( absc(cdwz(j,i,iw)) >= xloss*absc(cdwy(j,2,iw)) ) Cycle
          l4 = .True.
          Call ffcdwz(cdwz(1,1,iw),cz(1,iw),j,i,iw, &
                    & alpha(1),alpha(3),cpi,cpipj,cpiDpj,csdl2i, sdel2,6,ierw)
        End Do
      End Do
    End If
   End Do
  End If
 !  #]   get w-roots:
 !  #[   which case:
  If ( l4 ) Then
    If ( Aimag(alpha(1)) /= 0 ) Then
      l4pos = .False.
    Else If ( ierw >= 1 ) Then
      l4pos = .False.
    Else
     Write(ErrCan,*) "possible problem in ffcc0p",ierw
    End If
  End If
 !  #]   which case:
 !  #] get roots etc:
 !  #[ logarithms for 4point function:
  If ( npoin == 4 ) Then
    Do i = 1,3
      ismall(i) = 0
      If ( ilogi(i) /= -999 ) Cycle
      If ( isoort(2*i) /= 0 ) Then
 !		maybe add sophisticated factors i*pi later
        c = -cdyz(2,1,i)/cdyz(2,2,i)
        If ( absc(c-1) < xloss ) Then
         cs = cd2yzz(i)/cdyz(2,2,i)
         clogi(i) = Log1minusX(cs)
         ilogi(i) = 0
         ismall(i) = 1
        Else If ( Real(c,dp) > 0 ) Then
         clogi(i) = zfflog(c,0,czero)
         ilogi(i) = 0
        Else
          If ( absc(c+1) < xloss ) Then
           cs = -2*csdl2i(i)/cdyz(2,2,i)/ Real(cpi(i+3),dp)
           clogi(i) = Log1minusX(cs)
           ismall(i) = -1
          Else
           cs = 0
           clogi(i) = zfflog(-c,0,czero)
          End If
          If ( Aimag(c)<0 .Or. Aimag(cs)<0 ) Then
           ilogi(i) = -1
          Else If ( Aimag(c)>0 .Or. Aimag(cs)>0 ) Then
           ilogi(i) = +1
          Else If ( Real(cdyz(2,2,i),dp) == 0 ) Then
           ilogi(i)=-Nint(Sign(1._dp,Real(cpi(i+3),dp)))
           Print *,'  Doubtful imaginary part ',ilogi(i)
          End If
          If (Abs(Aimag(c))<precc*absc(c) .And.Abs(Aimag(cs))<precc*absc(cs)) &
          & Then
            Print *,'ffcc0p: error: imaginary part   Doubtful'
          End If
        End If
      End If
    End Do
    Do i=1,3
     j = i + 1
     If ( j == 4 ) j = 1
     If ( Abs(ismall(i)+ismall(j)) == 2 .And. absc(clogi(i)+&
      &clogi(j)) < xloss*absc(clogi(i)) ) Then
      Print *,'eerst: ',clogi(i)+clogi(j)
 !		assume that we got here because of complex sqrt(-delta)
      ci = Cmplx(0._dp,1._dp,dp)
      cs1=-2*ci*Aimag(cy(2,i))*csdl2i(j)/Real(cpi(j+3),dp) &
         & / (cdyz(2,2,i)*cdyz(2,2,j))
      cs2=-2*ci*Aimag(cy(2,j))*csdl2i(i)/Real(cpi(i+3),dp) &
         &  / (cdyz(2,2,i)*cdyz(2,2,j))
      cs = cs1 + cs2
      If ( absc(cs) < xloss*absc(cs1) ) Then
        k = j+1
          If ( k == 4 ) k = 1
        cs1 = cpipj(j+3,i+3)*cpi(j)
        cs2 = cpiDpj(k+3,j)*cpiDpj(j+3,j)
        cs4 = -cpiDpj(k+3,j)*cpiDpj(i+3,j)
        cs = cs1 + cs2 + cs4
        If ( absc(cs) < xloss*Max(absc(cs1),absc(cs2), absc(cs4)) ) Then
          Print *,'ffcc0p: cancellations in delj-deli'
          Cycle
        End If
        cs1=ci*Aimag(cy(2,j))*cs/(csdl2i(i)+csdl2i(j))
        Call ffcl2t(cs2,cpiDpj,k+3,j,4,5,6,+1,-1,6)
        cs2 = -cs2*csdl2i(j)/sdel2/Real(cpi(j+3),dp)
        cs = cs1 + cs2
        If ( absc(cs) < xloss*absc(cs1) ) Then
          Print *,'ffcc0p: cancellations in extra terms'
          Cycle
        End If
        cs = -2*cs/Real(cpi(i+3),dp)/(cdyz(2,2,i)*cdyz(2,2,j))
      End If
      clogi(i) = Log1minusX(cs)
      clogi(j) = 0
      Print *,'nu:    ',clogi(i)+clogi(j)
     End If
    End Do
  End If
 !  #] logarithms for 4point function:
 !  #[ integrals:
  If ( .Not. l4 .Or. .Not. l4pos ) Then
 !	    normal case
      Do i=1,3
       j = 2*i-1
        If ( isoort(2*i-1) /= 0 ) Then
          Call ffcs3(cs3(20*i-19),ipi12(2*i-1),cy(1,i), &
          cz(1,i),cdyz(1,1,i),cd2yzz(i),cpi,cpiDpj, i,6,isoort(j))
        End If
      End Do
    isoort(7) = 0
    isoort(8) = 0
  Else
    isoort(3) = jsoort(1)
    isoort(4) = jsoort(2)
    Call ffcs4(cs3(1),ipi12(1),cw(1,1),cy(1,1),cz(1,1),cdwy(1,1,1) &
              &   ,cdwz(1,1,1),cdyz(1,1,1), cd2yww(1),cd2yzz(1),cpi,cpiDpj, &
              & cpi(5)*alpha(3)**2,1,6,isoort(1))
    isoort(7) = jsoort(5)
    isoort(8) = jsoort(6)
    Call ffcs4(cs3(41),ipi12(1),cw(1,3),cy(1,3),cz(1,3),cdwy(1,1,3) &
              & ,cdwz(1,1,3),cdyz(1,1,3), cd2yww(3),cd2yzz(3),cpi,cpiDpj, &
              & cpi(5)*alpha(1)**2,3,6,isoort(5))
  End If
 !  #] integrals:
  End Subroutine ffcc0p

 Subroutine ffccyz(cy,cz,cdyz,cd2yzz,ivert,csdelp,csdels,etalam, &
            &      etami,delps,xpi,piDpj,isoort)
 !------------------------------------------------------------------------
 !	calculate in a numerically stable way
 !
 !	cz(1,2) = (-p(ip1).p(is2) +/- csdelp)/xpi(ip1)
 !	cy(1,2) = (-p(ip1).p(is2) +/- sdisc)/xpi(ip1)
 !			cdisc = csdels + etaslam*xpi(ip1)
 !
 !	cy(3,4) = 1-cy(1,2)
 !	cz(3,4) = 1-cz(1,2)
 !	cdyz(i,j) = cy(i) - cz(j)
 !
 !	Input:	ivert		(  Integer)	defines the vertex
 !		csdelp		(complex)	sqrt(lam(p1,p2,p3))/2
 !		csdels		(complex)	sqrt(lam(p,ma,mb))/2
 !		etalam		(complex)	det(si.sj)/det(pi.pj)
 !		etami(6)	(complex)	si.si - etalam
 !		xpi(ns)		(complex)	standard
 !		piDpj(ns,ns)	(complex)	standard
 !		ns		(  Integer)	dim of xpi,piDpj
 !
 !	Output:	cy(4),cz(4),cdyz(4,4)	(complex)	see above
 !		ier		(  Integer)	usual error flag
 !
 !	Calls:	fferr,ffroot
 !------------------------------------------------------------------------
  Implicit None
   Integer :: ivert,isoort(2)
   Complex(dp) :: cy(4),cz(4),cdyz(2,2),cd2yzz,csdelp,csdels
   Complex(dp) :: etalam,etami(6),delps,xpi(6),piDpj(6,6)

   Integer :: ip1,is1,is2,is3
   Complex(dp) :: cdisc

 !  #[ set up pointers:
   is1 = ivert
   is2 = ivert+1
   If ( is2 == 4 ) is2 = 1
   is3 = ivert-1
   If ( is3 == 0 ) is3 = 3
   ip1 = is1 + 3
 !  #] set up pointers:
 !  #[ xk = 0:
   If ( xpi(ip1) == 0 ) Then
     isoort(2) = 0
     If ( piDpj(is1,ip1) == 0 ) Then
      isoort(1) = 0
      Return
     End If
     If ( Aimag(etalam)/=0 ) Then
      isoort(1) = -1
     Else
      isoort(1) = -3
     End If
     cy(1) = etami(is2) / piDpj(is1,ip1) /2
     cy(2) = cy(1)
     cy(3) = - etami(is1) / piDpj(is1,ip1) /2
     cy(4) = cy(3)
     cz(1) = xpi(is2) / piDpj(is1,ip1) /2
     cz(2) = cz(1)
     cz(3) = - xpi(is1) / piDpj(is1,ip1) /2
     cz(4) = cz(3)
     cdyz(1,1) = - etalam / piDpj(is1,ip1) /2
     cdyz(1,2) = cdyz(1,1)
     cdyz(2,1) = cdyz(1,1)
     cdyz(2,2) = cdyz(1,1)
     Return
   End If
 !  #] xk = 0:
 !  #[ get cy(1,2),cz(1,2):
   If ( Aimag(etalam)/=0 ) Then
    isoort(1) = -1
    isoort(2) = -1
   Else
    isoort(1) = -3
    isoort(2) = -3
   End If
   Call Roots(xpi(ip1),piDpj(ip1,is2),xpi(is2), csdels,cz(1),cz(2))
   cdisc = delps/csdelp
   Call Roots(xpi(ip1),piDpj(ip1,is2),etami(is2), cdisc,cy(1),cy(2))
 !  #] get cy(1,2),cz(1,2):
 !  #[ get cy(3,4),cz(3,4):
   cz(4) = 1-cz(2)
   cz(3) = 1-cz(1)
   If ( absc(cz(3)) < xloss .Or. absc(cz(4)) < xloss ) Then
      Call Roots(xpi(ip1),-piDpj(ip1,is1), xpi(is1),csdels,cz(4),cz(3))
   End If
 !	the imaginary part may not be accurate in these cases, take
 !	some precautions:
   If ( cz(3) == 0 ) cz(1) = 1
   If ( cz(4) == 0 ) cz(2) = 1
   If ( Aimag(cz(1))==0 ) cz(1) = Cmplx(Real(cz(1),dp),-Aimag(cz(3)),dp)
   If ( Aimag(cz(2))==0 ) cz(2) = Cmplx(Real(cz(2),dp),-Aimag(cz(4)),dp)
   If ( Aimag(cz(1)) > 0 .Neqv. Aimag(cz(3)) < 0 ) Then
     If ( Abs(Real(cz(1),dp)) >= Abs(Real(cz(3),dp)) ) Then
      cz(1) = Cmplx(Real(cz(1),dp),-Aimag(cz(3)),dp)
     Else
      cz(3) = Cmplx(Real(cz(3),dp),-Aimag(cz(1)),dp)
     End If
   End If
   If ( Aimag(cz(2)) > 0 .Neqv. Aimag(cz(4)) < 0 ) Then
     If ( Abs(Real(cz(2),dp)) >= Abs(Real(cz(4),dp)) ) Then
      cz(2) = Cmplx(Real(cz(2),dp),-Aimag(cz(4)),dp)
     Else
      cz(4) = Cmplx(Real(cz(4),dp),-Aimag(cz(2)),dp)
     End If
   End If
   cy(4) = 1-cy(2)
   cy(3) = 1-cy(1)
   If ( absc(cy(3)) < xloss .Or. absc(cy(4)) < xloss ) Then
      Call Roots(xpi(ip1),-piDpj(ip1,is1), etami(is1),cdisc,cy(4),cy(3))
   End If
   If ( cy(3) == 0 ) cy(1) = 1
   If ( cy(4) == 0 ) cy(2) = 1
   If ( Aimag(cy(1))==0 ) cy(1) = Cmplx(Real(cy(1),dp),-Aimag(cy(3)),dp)
   If ( Aimag(cy(2))==0 ) cy(2) = Cmplx(Real(cy(2),dp),-Aimag(cy(4)),dp)
   If ( Aimag(cy(1)) > 0 .Neqv. Aimag(cy(3)) < 0 ) Then
     If ( Abs(Real(cy(1),dp)) >= Abs(Real(cy(3),dp)) ) Then
      cy(1) = Cmplx(Real(cy(1),dp),-Aimag(cy(3)),dp)
     Else
      cy(3) = Cmplx(Real(cy(3),dp),-Aimag(cy(1)),dp)
     End If
   End If
   If ( Aimag(cy(2)) > 0 .Neqv. Aimag(cy(4)) < 0 ) Then
     If ( Abs(Real(cy(2),dp)) >= Abs(Real(cy(4),dp)) ) Then
      cy(2) = Cmplx(Real(cy(2),dp),-Aimag(cy(4)),dp)
     Else
      cy(4) = Cmplx(Real(cy(4),dp),-Aimag(cy(2)),dp)
     End If
   End If
 !  #] get cy(3,4),cz(3,4):
 !  #[ get cdyz:
 !	Note that cdyz(i,j) only exists for i,j=1,2!
   If ( absc(cdisc+csdels) > xloss*absc(cdisc) ) Then
    cdyz(2,1) = ( cdisc + csdels )/xpi(ip1)
    cdyz(2,2) = etalam/(xpi(ip1)*cdyz(2,1))
   Else
    cdyz(2,2) = ( cdisc - csdels )/xpi(ip1)
    cdyz(2,1) = etalam/(xpi(ip1)*cdyz(2,2))
   End If
   cdyz(1,1) = -cdyz(2,2)
   cdyz(1,2) = -cdyz(2,1)
   cd2yzz = 2*cdisc/xpi(ip1)
 !  #] get cdyz:
  End Subroutine ffccyz


 Subroutine ffcdwz(cdwz,cz,i1,j1,l,calpha,calph1,cpi,cdpipj, &
                  & cpiDpj,csdeli,csdel2,ns,ier)
 !------------------------------------------------------------------------
 !	Recalculate cdwz(i1,j1) = cw(i1) - cz(j1)
 !------------------------------------------------------------------------
 Implicit None
  Integer :: i1,j1,l,ns,ier
  Complex(dp) :: cdwz(2,2),cz(4),calpha,calph1,cpi(ns)
  Complex(dp) :: cdpipj(ns,ns),cpiDpj(ns,ns),csdeli(3),csdel2

  Integer :: i,n
  Complex(dp) :: cs(8),csum,cfac,cddel
  Real(dp) :: xmax,afac

 !  #[ calculations:
  If ( l == 1 ) Then
      If ( j1 == 1 ) Then
        If ( absc(csdeli(1)+csdel2) < xloss*absc(csdel2) ) Then
 !		    for example in e-> e g* with eeg loop
 !		    first get the d  Ifference of csdeli(1) and csdel2:
        cs(1) = cpi(4)*cdpipj(2,5)
        cs(2) = -cpiDpj(4,3)*cpiDpj(4,2)
        cs(3) = cpiDpj(4,3)*cpiDpj(4,5)
        csum = cs(1)+cs(2)+cs(3)
        xmax = Max(absc(cs(1)),absc(cs(2)),absc(cs(3)))
        If ( absc(csum) < xloss*xmax ) Then
          ier = 1
          Return
        End If
        cddel = csum/(csdel2-csdeli(1))
        If ( i1 == 1 ) Then
          cs(1) = cpi(4)*csdeli(2)
        Else
          cs(1) = -cpi(4)*csdeli(2)
        End If
        cs(2) = cddel*cpiDpj(4,2)
        cs(3) = -cpiDpj(4,3)*csdeli(1)
        cs(4) = cpiDpj(4,3)*cpiDpj(4,5)
        cs(5) = -cpi(4)*cpiDpj(5,3)
        cs(6) = -cddel*csdel2
        n = 6
      Else
        ier = ier + 100
        Return
      End If
      csum = 0
      xmax = 0
      Do i=1,n
        csum = csum + cs(i)
        xmax = Max(xmax,absc(cs(i)))
      End Do
      If ( absc(csum) < xloss*xmax ) Then
        ier = ier + 1
      End If
      cdwz(i1,j1) = csum/calph1/cpi(4)/cpi(5)
      If ( cdwz(i1,j1) == 0 .And. csum /= 0 ) Then
        Print *,'?#$&!! cdwz = 0 but csum != 0, try again'
        afac = 1/absc(csum)
        csum = csum*Real(afac,dp)
        cdwz(i1,j1) = csum/calph1/cpi(4)/cpi(5)
        afac = 1/afac
        cdwz(i1,j1) = cdwz(i1,j1)*Real(afac,dp)
      End If
    Else
      ier = ier + 100
    End If

  Else If ( l == 3 ) Then
    If ( (i1==2 .And. j1==1) .Or. (i1==1 .And. j1==2 ) ) Then
      cfac = 1/(csdeli(2) + csdeli(3))
      cs(1) = cdpipj(6,5)*cz(j1)
      cs(2) = -calph1*cpi(5)*cz(j1+2)
      If ( Max(absc(cdpipj(2,1)),absc(cdpipj(5,6))) <&
        &  Max(absc(cdpipj(2,6)),absc(cdpipj(5,1))) ) Then
        cs(3) = cdpipj(2,1)/2
        cs(4) = cdpipj(5,6)/2
      Else
        cs(3) = cdpipj(2,6)/2
        cs(4) = cdpipj(5,1)/2
      End If
      cs(5) = cpiDpj(4,3)*cpiDpj(5,3)*cfac
      cs(6) = -cpiDpj(4,3)*cpiDpj(6,3)*cfac
      cs(7) = cpi(3)*cdpipj(5,6)*cfac
      If ( i1 == 1 ) Then
        csum = cs(1)+cs(2)+cs(3)+cs(4) - (cs(5)+cs(6)+cs(7))
      Else
        csum = cs(1)+cs(2)+cs(3)+cs(4) + cs(5)+cs(6)+cs(7)
      End If
      xmax = absc(cs(1))
      Do i=2,7
        xmax = Max(xmax,absc(cs(i)))
      End Do
      If ( absc(csum) < xloss*xmax ) Then
 !		    this result is not used   If it is not accurate (see
 !		    ffxc0p)
        ier = ier + 1
        xmax = xmax/absc(calpha*cpi(5))
        If ( xmax < Min(absc(cz(j1)),absc(cz(j1+2))) ) Then
          cdwz(i1,j1) = csum/(calpha*cpi(5))
        End If
      Else
        cdwz(i1,j1) = csum/(calpha*cpi(5))
      End If
    Else
      ier = ier + 100
    End If
  Else
    ier = ier + 100
  End If
 !  #] calculations:
  End Subroutine ffcdwz



 Subroutine ffcel2(del2,piDpj,ns,i1,i2,i3,lerr,ier)
 !--------------------------------------------------------------------
 !	calculate in a numerically stable way
 !	del2(piDpj(i1,i1),piDpj(i2,i2),piDpj(i3,i3)) =
 !		= piDpj(i1,i1)*piDpj(i2,i2) - piDpj(i1,i2)^2
 !		= piDpj(i1,i1)*piDpj(i3,i3) - piDpj(i1,i3)^2
 !		= piDpj(i2,i2)*piDpj(i3,i3) - piDpj(i2,i3)^2
 !	ier is the usual error flag.
 !--------------------------------------------------------------------
 Implicit None
  Integer :: ns,i1,i2,i3,lerr,ier
  Complex(dp) :: del2,piDpj(ns,ns)

  Complex(dp) :: s1,s2
 !
 !	calculations
 !
  If ( absc(piDpj(i1,i2)) < absc(piDpj(i1,i3)) .And.&
    &    absc(piDpj(i1,i2)) < absc(piDpj(i2,i3)) ) Then
    s1 = piDpj(i1,i1)*piDpj(i2,i2)
    s2 = piDpj(i1,i2)**2
  Else If ( absc(piDpj(i1,i3)) < absc(piDpj(i2,i3)) ) Then
    s1 = piDpj(i1,i1)*piDpj(i3,i3)
    s2 = piDpj(i1,i3)**2
  Else
    s1 = piDpj(i2,i2)*piDpj(i3,i3)
    s2 = piDpj(i2,i3)**2
  End If
  del2 = s1 - s2
  If ( absc(del2) < xloss*absc(s2) ) Then
      If ( lerr == 0 ) Then
 !		we know we have another chance
        If ( del2/=0 ) Then
         ier = ier + Int(Log10(xloss*absc(s2)/absc(del2)))
        Else
         ier = ier + Int(Log10(xloss*absc(s2)/xclogm))
        End If
      End If
  End If
  End Subroutine ffcel2


 Subroutine ffcl2p(delps1,xpi,dpipj,piDpj, ip1,ip2,ip3,is1,is2,is3,ns)
 !--------------------------------------------------------------------
 !	calculate in a numerically stable way
 !	delta_{ip1,is2}^{ip1,ip2}
 !	ier is the usual error flag.
 !--------------------------------------------------------------------
 Implicit None
  Integer :: ns,ip1,ip2,ip3,is1,is2,is3
  Complex(dp) :: delps1,xpi(ns),dpipj(ns,ns),piDpj(ns,ns)

  Complex(dp) :: s1,s2,s3,som
  Real(dp) :: xmax

 !  #[ stupid tree:
 !	1
  s1 = xpi(ip1)*piDpj(ip2,is2)
  s2 = piDpj(ip1,ip2)*piDpj(ip1,is2)
  delps1 = s1 - s2
  If ( absc(delps1) >= xloss*absc(s1) ) Return
  som = delps1
  xmax = absc(s1)
 !	2
  s1 = piDpj(ip1,ip2)*piDpj(ip3,is2)
  s2 = piDpj(ip1,ip3)*piDpj(ip2,is2)
  delps1 = s1 - s2
  If ( absc(delps1) >= xloss*absc(s1) ) Return
  If ( absc(s1) < xmax ) Then
    som = delps1
    xmax = absc(s1)
  End If
 !	3
  s1 = piDpj(ip1,ip3)*piDpj(ip1,is2)
  s2 = xpi(ip1)*piDpj(ip3,is2)
  delps1 = s1 - s2
  If ( absc(delps1) >= xloss*absc(s1) ) Return
  If ( absc(s1) < xmax ) Then
    som = delps1
    xmax = absc(s1)
  End If
 !	4
  s1 = xpi(ip1)*piDpj(ip2,is1)
  s2 = piDpj(ip1,is1)*piDpj(ip1,ip2)
  delps1 = s1 - s2
  If ( absc(delps1) >= xloss*absc(s1) ) Return
  If ( absc(s1) < xmax ) Then
    som = delps1
    xmax = absc(s1)
  End If
 !	5
  s1 = piDpj(ip1,is2)*piDpj(ip2,is1)
  s2 = piDpj(ip1,is1)*piDpj(ip2,is2)
  delps1 = s1 - s2
  If ( absc(delps1) >= xloss*absc(s1) ) Return
  If ( absc(s1) < xmax ) Then
    som = delps1
    xmax = absc(s1)
  End If
 !	6
  s1 = piDpj(ip1,ip2)*piDpj(ip3,is1)
  s2 = piDpj(ip1,ip3)*piDpj(ip2,is1)
  delps1 = s1 - s2
  If ( absc(delps1) >= xloss*absc(s1) ) Return
  If ( absc(s1) < xmax ) Then
    som = delps1
    xmax = absc(s1)
  End If
 !	7
  s1 = piDpj(ip2,is2)*piDpj(ip3,is1)
  s2 = piDpj(ip2,is1)*piDpj(ip3,is2)
  delps1 = s1 - s2
  If ( absc(delps1) >= xloss*absc(s1) ) Return
  If ( absc(s1) < xmax ) Then
    som = delps1
    xmax = absc(s1)
  End If
 !	8
  s1 = piDpj(ip1,ip3)*piDpj(ip1,is1)
  s2 = xpi(ip1)*piDpj(ip3,is1)
  delps1 = s1 - s2
  If ( absc(delps1) >= xloss*absc(s1) ) Return
  If ( absc(s1) < xmax ) Then
    som = delps1
    xmax = absc(s1)
  End If
 !	9
  s1 = piDpj(ip1,is1)*piDpj(ip3,is2)
  s2 = piDpj(ip1,is2)*piDpj(ip3,is1)
  delps1 = s1 - s2
  If ( absc(delps1) >= xloss*absc(s1) ) Return
  If ( absc(s1) < xmax ) Then
    som = delps1
    xmax = absc(s1)
  End If
 !10	22-nov-1993 yet another one
  If ( dpipj(1,1)==0 ) Then
    s1 = +xpi(ip1)*dpipj(is3,is2)/2
    s2 = -piDpj(ip1,ip2)*dpipj(is2,is1)/2
    s3 = +xpi(ip1)*piDpj(ip2,ip3)/2
    delps1 = s1+s2+s3
      If ( absc(delps1) >= xloss*Max(absc(s1),absc(s2)) ) Return
      If ( Max(absc(s1),absc(s2)) < xmax ) Then
      som = delps1
      xmax = absc(s1)
      End If
  End If
 !	NO possibility
  delps1 = som

 !  #] stupid tree:
  End Subroutine ffcl2p

 Subroutine ffcl2t(delps,piDpj,in,jn,kn,ln,lkn,islk,iss,ns)
 !--------------------------------------------------------------------
 !	calculate in a numerically stable way
 !
 !		\delta_{si,sj}^{sk,sl}
 !
 !	with p(lk) = islk*(iss*sl - sk)	(islk,iss = +/-1)
 !	and NO relationship between s1,s2 assumed (so 1/2 the
 !	possibilities of ffdl2s).
 !--------------------------------------------------------------------
  Implicit None
   Integer :: in,jn,kn,ln,lkn,islk,iss,ns
   Complex(dp) :: delps,piDpj(ns,ns)

   Complex(dp) :: s1,s2
 !  #[ calculations:
   If ( in == jn ) Then
    delps = 0._dp
    Return
   End If
   s1 = piDpj(kn,in)*piDpj(ln,jn)
   s2 = piDpj(ln,in)*piDpj(kn,jn)
   delps = s1 - s2
   If ( absc(delps) >= xloss*absc(s1) ) Return
    s1 = piDpj(kn,in)*piDpj(lkn,jn)
    s2 = piDpj(lkn,in)*piDpj(kn,jn)
    delps = iss*islk*(s1 - s2)
   If ( absc(delps) >= xloss*absc(s1) ) Return
    s1 = piDpj(lkn,in)*piDpj(ln,jn)
    s2 = piDpj(ln,in)*piDpj(lkn,jn)
    delps = islk*(- s1 + s2)
   If ( absc(delps) >= xloss*absc(s1) ) Return
 !  #] calculations:
  End Subroutine ffcl2t

  Subroutine ffcl3m(del3mi,ldel,del3,del2,xpi,dpipj,piDpj,ns,ip1n, &
    ip2n,ip3n,is,itime)
 !--------------------------------------------------------------------
 !	Calculate xpi(i)*del2 - del3(piDpj)
 !
 !	  /  si	mu \2		(This appears to be one of the harder
 !	= | d	   |		 determinants to calculate accurately.
 !	  \  p1	p2 /		 Note that we allow a loss of xloss^2)
 !
 !	Input:	ldel		  Iff .true. del2 and del3 exist
 !		del3		\delta^{s(1),p1,p2}_{s(1),p1,p2}
 !		del2		\delta^{p1,p2}_{p1,p2}
 !		xpi(ns)		standard
 !		dpipj(ns,ns)	standard
 !		piDpj(ns,ns)	standard
 !		ipi		pi = xpi(abs(ipi)) [p3=-p1 +/-p2]
 !		is		si = xpi(is,is+1,..,is+itime-1)
 !		itime		number of functions to calculate
 !
 !	Output:	del3mi(3)	(\delta^{s_i \mu}_{p_1 p_2})^2
 !--------------------------------------------------------------------
 Implicit None
  Integer :: ns,ip1n,ip2n,ip3n,is,itime
  Logical :: ldel
  Complex(dp) :: del3mi(itime),del3,del2,xpi(ns),dpipj(ns,ns), piDpj(ns,ns)

  Real(dp) :: smax,xmax
  Complex(dp) :: s(7),som,xsom
  Integer :: i,j,k,ip1,ip2,ip3,ipn,is1,is2,isi,is3,ihlp,iqn, &
    jsgn1,jsgn2,jsgn3,jsgnn,iadj(10,10,3:4),init,nm
  Save iadj,init
  Logical :: lmax,ltwist

 !
 !	data
 !
  Data iadj /200*0/
  Data init /0/
 !  #[ initialisations:
  If ( init == 0 ) Then
    init = 1
 !
 !	    Fill the array with adjacent values:   If
 !		x = iadj(i,j)
 !		k = abs(mod(k,100))
 !		jsgnk = sign(x)
 !		jsgnj = 1-2*theta(x-100)  (ie -1   Iff |x|>100)
 !	    then
 !		pi(k) = jsgnk*( p(i) - jsgnj*pi(j) )
 !
      Do nm=3,4
        Do i=1,nm
        is1 = i
        is2 = i+1
          If ( is2 > nm ) is2 = 1
        is3 = i-1
          If ( is3 == 0 ) is3 = nm
        ip1 = is1 + nm
        iadj(is1,is2,nm) = -ip1
        iadj(is2,is1,nm) = ip1
        iadj(ip1,is2,nm) = -is1
        iadj(is2,ip1,nm) = is1
        iadj(is1,ip1,nm) = 100+is2
        iadj(ip1,is1,nm) = 100+is2
          If ( nm == 3 ) Then
          iadj(ip1,is2+3,3) = -100-is3-3
          iadj(is2+3,ip1,3) = -100-is3-3
          End If
        End Do
      End Do
    iadj(3,1,4) = -9
    iadj(1,3,4) = 9
    iadj(9,1,4) = -3
    iadj(1,9,4) = 3
    iadj(3,9,4) = 100+1
    iadj(9,3,4) = 100+1
    iadj(2,4,4) = -10
    iadj(4,2,4) = 10
    iadj(10,4,4) = -2
    iadj(4,10,4) = 2
    iadj(2,10,4) = 100+4
    iadj(10,2,4) = 100+4
  End If
  If ( ns == 6 ) Then
    nm = 3
  Else
    nm = 4
  End If
 !  #] initialisations:
 !  #[ easy tries:
  Do i=1,itime
    isi = i+is-1
    lmax = .False.
 !
 !	    get xpi(isi)*del2 - del3 ...   If del3 and del2 are defined
 !
      If ( ldel ) Then
      s(1) = xpi(isi)*del2
      som = s(1) - del3
      smax = absc(s(1))
        If ( absc(som) >= xloss**2*smax ) Goto 35
      xsom = som
      xmax = smax
      lmax = .True.
      End If
    ip1 = ip1n
    ip2 = ip2n
    ip3 = ip3n
      Do j=1,3
 !
 !		otherwise use the simple threeterm formula
 !
      s(1) = xpi(ip2)*piDpj(ip1,isi)**2
      s(2) = xpi(ip1)*piDpj(ip2,isi)*piDpj(ip2,isi)
      s(3) = -2*piDpj(ip2,isi)*piDpj(ip2,ip1)*piDpj(ip1,isi)
      som = s(1) + s(2) + s(3)
      smax = Max(absc(s(1)),absc(s(2)),absc(s(3)))
        If ( absc(som) >= xloss**2*smax ) Goto 35
        If ( .Not. lmax .Or. smax < xmax ) Then
        xsom = som
        xmax = smax
        lmax = .True.
        End If
 !
 !		  If there are cancellations between two of the terms:
 !		we try mixing with isi.
 !
 !		First map cancellation to s(2)+s(3) (  Do not mess up
 !		rotations...)
 !
        If ( absc(s(1)+s(3)) < absc(s(3))/2 ) Then
        ihlp = ip1
        ip1 = ip2
        ip2 = ihlp
        som = s(1)
        s(1) = s(2)
        s(2) = som
        ltwist = .True.
        Else
        ltwist = .False.
        End If
        If ( absc(s(2)+s(3)) < absc(s(3))/2 ) Then
 !
 !		switch to the vector pn so that si = jsgn1*p1 + jsgnn*pn
 !
        k = iadj(isi,ip1,nm)
          If ( k /= 0 ) Then
          ipn = Abs(k)
          jsgnn = isign(1,k)
            If ( ipn > 100 ) Then
            ipn = ipn - 100
            jsgn1 = -1
            Else
            jsgn1 = +1
            End If
            If ( absc(dpipj(ipn,isi)) < xloss*absc(piDpj(ip1,isi)) .And.&
            &absc(piDpj(ipn,ip2)) < xloss*absc(piDpj(ip2,isi)) ) Then
 !		same:	s(1) = xpi(ip2)*piDpj(ip1,isi)**2
            s(2) = jsgnn*piDpj(isi,ip2)*piDpj(ipn,ip2)*xpi(ip1)
            s(3) = jsgn1*piDpj(isi,ip2)*piDpj(ip1,ip2)*dpipj(ipn,isi)
            som = s(1) + s(2) + s(3)
            smax = Max(absc(s(1)),absc(s(2)),absc(s(3)))
 !			print *,'    (isi+ip1) with isi,ip1,ip2,ipn: ',
 !     +				isi,ip1,ip2,ipn
 !			print *,'xpi(ip2),piDpj(ip1,isi)',xpi(ip2),
 !     +				piDpj(ip1,isi)
 !			print *,'piDpj(isi,ip2),piDpj(ipn,ip2),xpi(ip1)'
 !     +				,piDpj(isi,ip2),piDpj(ipn,ip2),xpi(ip1)
              If ( absc(som) >= xloss**2*smax ) Goto 35
              If ( smax < xmax ) Then
              xsom = som
              xmax = smax
              End If
 !
 !			there may be a cancellation between s(1) and
 !			s(2) left.  Introduce a vector q such that
 !			pn = jsgnq*q + jsgn2*p2.  We also need the sign
 !			jsgn3 in p3 = -p1 - jsgn3*p2
 !
            k = iadj(ipn,ip2,nm)
              If ( k /= 0 ) Then
              iqn = Abs(k)
 !not used		    jsgnq = isign(1,k)
                If ( iqn > 100 ) Then
                iqn = iqn - 100
                jsgn2 = -1
                Else
                jsgn2 = +1
                End If
              k = iadj(ip1,ip2,nm)
                If ( k == 0 .Or. k < 100 ) Then
 !				we have p1,p2,p3 all p's
                jsgn3 = +1
                Else If ( k < 0 ) Then
 !				ip1,ip2 are 2*s,1*p such that p2-p1=ip3
                jsgn3 = -1
                Else
                jsgn3 = 0
                End If
 !			    we need one condition on the signs for this
 !			    to work
                If ( ip3/=0 .And. jsgn1*jsgn2==jsgnn*&
                &jsgn3 .And. absc(s(3))<xloss*smax ) Then
                s(1) = piDpj(ip1,isi)**2*dpipj(iqn,ipn)
                s(2) = -jsgn2*jsgn1*piDpj(ipn,ip2)*piDpj(ip1,isi)*dpipj(ipn,isi)
 !				s(3) stays the same
                s(4) = -jsgn2*jsgn1*piDpj(ipn,ip2)*xpi(ip1)*piDpj(isi,ip3)
                som = s(1) + s(2) + s(3) + s(4)
                smax = Max(absc(s(1)),absc(s(2)), absc(s(3)),absc(s(4)))
                  If (absc(som)>=xloss**2*smax) Goto 35
                  If ( smax < xmax ) Then
                  xsom = som
                  xmax = smax
                  End If
                End If
              End If
            End If
          End If
        k = iadj(isi,ip2,nm)
          If ( k /= 0 ) Then
          ipn = Abs(k)
          jsgnn = isign(1,k)
            If ( ipn > 100 ) Then
            jsgn1 = -1
            ipn = ipn - 100
            Else
            jsgn1 = +1
            End If
            If ( absc(dpipj(ipn,isi)) < xloss*absc(piDpj(ip2,isi)) .And.&
            &absc(piDpj(ipn,ip1)) < xloss*absc(piDpj(ip1,isi)) ) Then
            s(1) = jsgnn*piDpj(isi,ip1)*piDpj(ipn,ip1)*xpi(ip2)
            s(2) = xpi(ip1)*piDpj(ip2,isi)**2
            s(3) = jsgn1*piDpj(isi,ip1)*piDpj(ip2,ip1)*dpipj(ipn,isi)
            som = s(1) + s(2) + s(3)
            smax = Max(absc(s(1)),absc(s(2)),absc(s(3)))
            Print *,'    (isi+ip2) with isi,ip1,ip2,ipn: ', isi,ip1,ip2,ipn
              If ( absc(som) >= xloss**2*smax ) Goto 35
              If ( smax < xmax ) Then
              xsom = som
              xmax = smax
              End If
            End If
          End If
        End If
 !
 !		rotate the ipi
 !
        If ( ip3 == 0 ) Goto 30
        If ( j /= 3 ) Then
          If ( .Not. ltwist ) Then
          ihlp = ip1
          ip1 = ip2
          ip2 = ip3
          ip3 = ihlp
          Else
          ihlp = ip2
          ip2 = ip3
          ip3 = ihlp
          End If
        End If
      End Do
 !  #] easy tries:
 !  #[ choose the best value:
 !
 !	    These values are the best found:
 !
 30   som = xsom
    smax = xmax
 35   del3mi(i) = som
  End Do
 !  #] choose the best value:
  End Subroutine ffcl3m

 Subroutine ffclgy(cs3,ipi12,ntot,cy,cz,cd2yzz,ier)
 !--------------------------------------------------------------------
 !	calculates the the difference of two S's with cy(3,4),cz(3,4),
 !	cy(4)cz(3)-cy(3)cz(4) given.  Note the difference with ffdcs4,
 !	in which the cy's are the same and only the cz's d  Ifferent.
 !	Here both can be d  Ifferent.	Also we skip an intermediat
 !	level.
 !
 !	input:	cy(4)	     (complex)	cy,1-cy in S with s3,s4
 !		cz(4)	     (complex)	cz,1-cz in S with s3,s4
 !		cdyz(2,2)    (complex)	cy - cz
 !		cd2yzz	     (complex)	2*cy - cz+ - cz-
 !		cdyzzy(4)    (complex)	cy(i,4)*cz(i,4)-cy(i,3)*cz(i,4)
 !		cpiDpj(6,6)  (complex)	usual
 !		cs3	     (complex)	assumed zero.
 !
 !	output: cs3	     (complex)	mod factors pi^2/12, in array
 !		ipi12	     (  Integer)	these factors
 !		isoort	     (  Integer)	returns kind of action taken
 !		ier	     (  Integer)	number of digits lost
 !--------------------------------------------------------------------
 Implicit None
  Complex(dp) :: cs3
  Complex(dp) :: cy(4),cz(4),cd2yzz
  Integer :: ipi12,ntot,ier

  Integer :: ipi
  Complex(dp) :: c,cc,clogy,c2y1,csum

 !  #[ calculations:
  ipi = 0
  If ( 1 < xloss*absc(cy(2)) ) Then
    clogy = Log1minusX(1._dp/cy(2))
  Else
    If ( absc(cy(2)) < xclogm .Or. absc(cy(4)) < xclogm ) Then
      If ( ntot /= 0 )   Call WriteLFerror(22)
      clogy = 0
    Else
      c = -cy(4)/cy(2)
      If ( Real(c,dp) > -Abs(Aimag(c)) ) Then
        clogy = zfflog(c,0,czero)
      Else
 !		    take out the factor 2*pi^2
        cc = c+1
        If ( absc(cc) < xloss ) Then
          c2y1 = -cd2yzz - cz(1) + cz(4)
          If ( absc(c2y1) < xloss*Max(absc(cz(1)), absc(cz(4))) ) Then
            c2y1 = -cd2yzz - cz(2) + cz(3)
          End If
          csum = -c2y1/cy(2)
          clogy = Log1minusX(csum)
        Else
          csum = 0
          clogy = zfflog(-c,0,czero)
        End If
        If (Aimag(c) < -precc*absc(c) .Or.Aimag(csum) < -precc*absc(csum)) Then
          ipi = -1
        Else If ( Aimag(c) > precc*absc(c) .Or.&
            &  Aimag(csum) > precc*absc(csum) ) Then
          ipi = +1
        Else
          Call WriteLFerror(18)
          ipi = 0
        End If
      End If
    End If
  End If
  cs3 = cs3 + ntot*c2ipi*clogy
  If ( ipi /= 0 ) Then
    ipi12 = ipi12 - 24*ntot*ipi
  End If
 !  #] calculations:
  End Subroutine ffclgy

 Subroutine ffcrr(crr,ipi12,cy,cy1,cz,cz1,cdyz,ld2yzz,cd2yzz,czz, &
                & czz1,isoort,ieps,ier)
 !------------------------------------------------------------------------
 !	calculates R as defined in appendix b:
 !
 !			/1   log(y-y1+ieps) - log(y0-y1+ieps)
 !	r(y0,y1,iesp) = \ dy --------------------------------
 !			/0		y-y0
 !
 !	    = li2(c1) - li2(c2)
 !		+ eta(-y1,1/(y0-y1))*log(c1)
 !		- eta(1-y1,1/(y0-y1))*log(c2)
 !	with
 !	    c1 = y0 / (y0-y1), c2 = (y0-1) / (y0-y1)
 !
 !	the factors pi^2/12 are passed separately in the integer ipi12
 !	ier is a status flag: 0=ok, 1=numerical problems, 2=error
 !
 !	Input:	cy	(complex)
 !		cy1	(complex)	1-y
 !		cz	(complex)
 !		cz1	(complex)	1-z
 !		cdyz	(complex)	y-z
 !		ieps	(  Integer)	denotes sign imaginary part of
 !					argument logs (0:   Don't care;
 !					+/-1: add -ieps to z; +/-2:
 !					direct in dilogs, no eta's)
 !
 !	Output	crr	(complex)	R modulo factors pi^2/12
 !		ipi12	(  Integer)	these factors
 !		ier	(  Integer)	lost ier digits, >100: error
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ipi12,isoort,ieps,ier
  Logical :: ld2yzz,lreal
  Complex(dp) :: crr(7),cy,cy1,cz,cz1,cdyz,cd2yzz,czz,czz1

  Complex(dp) :: cfact,cc1,cc2,cc1p,cc2p,carg1,carg2,carg3, &
    cli1,cli2,cli3,clo1,clo2,clo3,clog1p,clog2p,chill, &
    cd2,cd21,cd2n,cd21n1,cc1n,cterm,ctot,clog1,clog2, cli4,clo4
  Complex(dp) :: ctroep
  Real(dp), Save :: xprec,bndtay
  Real(dp) :: xa,xr, x00(3)
  Real(dp) :: y,y1,z,z1,dyz,d2yzz,zz,zz1
  Integer :: i,iclas1,iclas2,n1,n2,n3,ntot, i2pi,n3p

 !  #[ initialisations:
  Data xprec /-1._dp/
  If ( xprec /= precx ) Then
    xprec = precx
!    bndtay = ffbnd(2,18,xn2inv)
    bndtay = (precx*Abs(xn2inv(2)/xn2inv(20)))**(1./18._dp)
 !	    print *,'bndtay = ',bndtay
  End If
 !  #] initialisations:
 !  #[ real case:
  If ( Aimag(cy)==0 .And. Aimag(cy1)==0 .And. Aimag(cz)==0&
    & .And. Aimag(cz1)==0 ) Then
    y  = Real(cy,dp)
    y1 = Real(cy1,dp)
    z  = Real(cz,dp)
    z1 = Real(cz1,dp)
    dyz = Real(cdyz,dp)
    d2yzz = Real(cd2yzz,dp)
    zz = Real(czz,dp)
    zz1 = Real(czz1,dp)
    Call ffcxr(crr,ipi12,y,y1,z,z1,dyz,ld2yzz,d2yzz,zz,zz1 &
             & , .False.,x00,ieps,ier)
    Return
  End If
 !  #] real case:
 !
  xa = absc(cdyz)
  If ( xa == 0 ) Then
    Return
 !	This line is for 68000 compilers that have a limited range for
 !	complex division (Absoft, Apollo, Gould NP1):
  Else If ( Real(cdyz,dp) < xclogm .Or. Aimag(cdyz) < xclogm&
      & .Or. 1/xa < xclogm ) Then
    ctroep = cdyz*Real(1/xa,dp)
    cfact = 1/ctroep
    cfact = Real(1/xa,dp)*cfact
  Else
    cfact = 1/cdyz
  End If
  cc1 = cy * cfact
  cc2 = - cy1 * cfact
 !
 !	see   If we just need the real part
 !
  lreal = Mod(isoort,5) == 0
 !  #[ which area?:
 !
 !	determine the area:	1={|x|<=1,Re(x)<=1/2},
 !				2={|1-x|<=1,Re(x)>1/2}
 !				3={|x|>1,|1-x|>1}
 !
  xr = Real(cc1,dp)
  xa = absc(cc1)
  If ( xa > 1 .And. xa < 1+Sqrt(2.) ) Then
 !	    we need a more accurate estimate
    xa = xr**2 + Aimag(cc1)**2
  End If
  If ( ld2yzz .And. absc(cc1+1) < xloss/2 ) Then
    iclas1 = 4
    cc1p = cc1
  Else If ( xa <= 1 .And. xr <= 0.5 ) Then
    iclas1 = 1
    cc1p = cc1
  Else If ( xa < 1+Sqrt(2.) .And. xa < 2*xr ) Then
    iclas1 = 2
    cc1p = -cz * cfact
    If ( Abs(Aimag(cc1p)) < precc*Abs(Real(cc1p,dp)) ) cc1p = Real(cc1p,dp)
  Else
    iclas1 = 3
    If ( 1/xa < xclogm ) Then
      ctroep = cc1*Real(1/xa,dp)
      ctroep = 1/ctroep
      cc1p = ctroep*Real(1/xa,dp)
    Else
      cc1p = 1/cc1
    End If
  End If
  xr = Real(cc2,dp)
  xa = absc(cc2)
  If ( xa > 1 .And. xa < 1+Sqrt(2.) ) Then
    xa = xr**2 + Aimag(cc2)**2
  End If
  If ( ld2yzz .And. absc(cc2+1) < xloss ) Then
    iclas2 = 4
    cc2p = cc2
  Else If ( xa <= 1 .And. xr <= 0.5 ) Then
    iclas2 = 1
    cc2p = cc2
  Else If ( xa < 1+Sqrt(2._dp) .And. xa < 2*xr ) Then
    iclas2 = 2
    cc2p = cz1 * cfact
      If ( Abs(Aimag(cc2p)) < precc*Abs(Real(cc2p,dp)) ) cc2p = Real(cc2p,dp)
  Else
    iclas2 = 3
    If ( 1/xa < xclogm ) Then
      ctroep = cc2*Real(1/xa,dp)
      ctroep = 1/ctroep
      cc2p = ctroep*Real(1/xa,dp)
    Else
      cc2p = 1/cc2
    End If
  End If
 !
 !	throw together   If they are close
 !
  If ( iclas1 /= iclas2 .And. absc(cc1-cc2) < 2*xloss ) Then
 !	    we   Don't want trouble with iclasn = 4
    If ( iclas1 == 4 ) iclas1 = 1
    If ( iclas2 == 4 ) iclas2 = 1
    If ( iclas1 == iclas2 ) Goto 5
 !	    go on
    If ( iclas1 <= iclas2 ) Then
      iclas2 = iclas1
      If ( iclas1 == 1 ) Then
        cc2p = cc2
      Else
        cc2p = cz1*cfact
      End If
    Else
      iclas1 = iclas2
      If ( iclas1 == 1 ) Then
        cc1p = cc1
      Else
        cc1p = -cz*cfact
      End If
    End If
  End If

 !  #] which area?:
 !  #[ eta's:
 !
 !	get eta1 and eta2
 !
 5 If ( Abs(ieps) >= 2 .Or. isoort == -2 ) Then
    n1 = 0
    n2 = 0
  Else
      If ( Aimag(cz) == 0 .Or. Aimag(cz1) == 0 ) Then
        If ( Aimag(cz1) == 0 ) Then
          If ( Aimag(cz) == 0 ) Then
 !		    cz is really real, the hard case:
            If ( cz == 0 ) Then
 !			multiplied with log(1), so   Don't care:
              n1 = 0
 !			look at ieps for guidance
 !	n2 = nffet1(Cmplx(Real(0),Real(ieps)),cfact,cfact,ier) = 0
              n2 = 0
            Else If ( cz1 == 0 ) Then
              n1 = nffet1(Cmplx(0._dp,Real(ieps,dp),dp),cfact, -cfact,ier)
              n2 = 0
            Else
              n1 = nffet1(Cmplx(0._dp,Real(ieps,dp),dp),cfact, -cz*cfact,ier)
              n2 = nffet1(Cmplx(0._dp,Real(ieps,dp),dp),cfact, cz1*cfact,ier)
            End If
          Else
            n1 = nffet1(-cz,cfact,-cz*cfact,ier)
            n2 = nffet1(-cz,cfact,cz1*cfact,ier)
          End If
        Else
          n1 = nffet1(cz1,cfact,-cz*cfact,ier)
          n2 = nffet1(cz1,cfact,cz1*cfact,ier)
        End If
      Else
 !	    the imaginary part of cc1, cc1p is often very unstable.
 !	    make sure it agrees with the actual sign used.
        If ( iclas1 == 2 ) Then
          If ( Aimag(cc1p) == 0 ) Then
 !		      If y (or y1 further on) is purely imaginary
 !		    give a ran  Dom sh  Ift, this will also be used in
 !		    the transformation terms.  Checked 7-mar-94 that it
 !		    is independent of the sign used.
            If ( Real(cy,dp)==0 ) cy = cy +isgnal*Real(precc,dp)*Aimag(cy)
            n1 = nffet1(-cz,cfact,Cmplx(0._dp,ieps*Real(cy,dp),dp), ier)
          Else
            n1 = nffet1(-cz,cfact,cc1p,ier)
          End If
        Else
          If ( Aimag(cc1) == 0 ) Then
            If ( Real(cy1)==0 ) cy1 = cy1 +isgnal*Real(precc,dp)*Aimag(cy)
            n1 = nffet1(-cz,cfact,Cmplx(0._dp, -ieps*Real(cy1,dp),dp),ier)
          Else
            n1 = nffet1(-cz,cfact,-cc1,ier)
          End If
        End If
        If ( iclas2 == 2 ) Then
          If ( Aimag(cc2p) == 0 ) Then
            If ( Real(cy)==0 ) cy = cy +isgnal*Real(precc,dp)*Aimag(cy)
            n2 = nffet1(cz1,cfact,Cmplx(0._dp,ieps*Real(cy,dp),dp), ier)
          Else
            n2 = nffet1(cz1,cfact,cc2p,ier)
          End If
        Else
          If ( Aimag(cc2) == 0 ) Then
            If ( Real(cy1)==0 ) cy1 = cy1 +isgnal*Real(precc,dp)*Aimag(cy)
            n2 = nffet1(cz1,cfact,Cmplx(0._dp, -ieps*Real(cy1,dp),dp),ier)
          Else
            n2 = nffet1(cz1,cfact,-cc2,ier)
          End If
        End If
      End If
  End If
 !  #] eta's:
 !  #[ calculations:
 !	3-oct-1995 changed code to only use second criterium   If the
 !	Taylor expansion is used - otherwise the Hill identity will
 !	only make things worse
  If ( iclas1 == iclas2 .And. isoort /= -2 .And. &
          ( absc(cc1p-cc2p) < 2*xloss*absc(cc1p)  &
    .Or. lreal .And. Abs(Real(cc1p-cc2p,dp)) < 2*xloss*Abs(Real(cc1p,dp))&
    & .And. (Abs(Real(cc2p,dp)) +Aimag(cc2p)**2/4) < xloss .And.&
    & Abs(Aimag(cc2p)) < bndtay ) ) Then
 !	    Close together:
 ! -#[	    handle dilog's:
      If ( .Not. lreal .And. absc(cc2p) > xloss&
      & .Or. lreal .And. ( (Abs(Real(cc2p,dp)) + Aimag(cc2p)**2/4) &
      > xloss .Or. Abs(Aimag(cc2p)) > bndtay ) ) Then
 !--#[		Hill identity:
 !
 !		Use the Hill identity to get rid of the cancellations.
 !
 !
 !
        If ( iclas1 == 1 .Or. iclas1 == 4 ) Then
          carg1 = 1/cy
          carg2 = 1/cz1
          carg3 = carg2/cc1p
        Else If ( iclas1 == 2 ) Then
          carg1 = 1/cz
          carg2 = 1/cy1
          carg3 = carg2/cc1p
        Else If ( iclas1 == 3 ) Then
          carg1 = 1/cy1
          carg3 = 1/cz1
          carg2 = carg3*cc1p
        End If
        Call ffzli2(cli1,clo1,carg1)
        Call ffzli2(cli2,clo2,carg2)
        Call ffzli2(cli3,clo3,carg3)
        If ( absc(cc2p) < xloss ) Then
          clog2p = Log1minusX(cc2p)
        Else
          clog2p = zfflog(1-cc2p,0,czero)
        End If
        chill = clo1*clog2p
 !--#]		Hill identity:
      Else
 !--#[		Taylor expansion:
 !
 !		  If the points are close to zero   Do a Taylor
 !		expansion of the first and last dilogarithm
 !
 !			Li2(cc1p) - Li2(cc2p)
 !			  = sum cc1p^i ( 1-(1-cd2)^i ) /i^2
 !
 !		with cd2 = 1-cc2p/cc1p = ...
 !
        If ( iclas1 == 1 .Or. iclas1 == 4 ) Then
          cd2 = 1/cy
        Else If ( iclas1 == 2 ) Then
          cd2 = 1/cz
        Else If ( iclas1 == 3 ) Then
          cd2 = 1/cy1
        End If
        cd21 = 1-cd2
        cd21n1 = 1
        cc1n = cc1p
        cd2n = cd2
        ctot = cc1p*cd2
        Do i=2,20
         cc1n = cc1n*cc1p
         cd21n1 = cd21n1*cd21
         cd2n = cd2n + cd2*cd21n1
         cterm = cc1n*cd2n*Real(xn2inv(i),dp)
         ctot = ctot + cterm
          If ( absc(cterm) <= precc*absc(ctot) .Or.&
          &lreal .And. Abs(Real(cterm,dp)) <= precc*Abs(Real(ctot,dp)) ) Exit
        End Do

        cli1 = ctot
        cli2 = 0
        cli3 = 0
        chill = 0
 !		for the eta+transformation section we also need
        If ( iclas1/=1 .Or. n1/=0 .Or. n2/=0 ) clo1 = Log1minusX(cd2)
        If ( iclas1==2 ) clo2 = Log1minusX(1/cy1)
 !--#]		Taylor expansion:
      End If
 !
 ! -#]	    handle dilog's:
 ! -#[	    handle eta + transformation terms:
      If ( iclas1==1 .Or. iclas1==4 ) Then
 !--#[		no transformation:
 !
 !		no transformation was made.
 !
 !		crr(5) = 0
        If ( n1 /= n2 ) Then
          If ( absc(cc1) < xclogm ) Then
            Call WriteLFerror(17)
          Else
 !			imaginary part not checked
            ier = ier + 50
            crr(5) = (n1-n2)*c2ipi*zfflog(cc1,ieps,-cy)
          End If
        End If
 !		crr(6) = 0
 !		crr(7) = 0
        If ( n2/=0 ) Then
          crr(6) = - n2*c2ipi*clo1
          n3 = nffeta(cc2,1/cc1,ier)
          If ( n3 /= 0 ) Then
            crr(7) = n2*n3*c2ipi**2
 !		      Else
 !			crr(7) = 0
          End If
        End If
 !--#]		no transformation:
      Else If ( iclas1 == 2 ) Then
 !--#[		transform 1-x:
 !
 !		we tranformed to 1-x for both dilogs
 !
        If ( absc(cc1p) < xloss ) Then
          clog1 = Log1minusX(cc1p)
        Else
          clog1 = zfflog(cc1,ieps,-cy)
        End If
        If ( Aimag(cc2p)==0 ) Then
          If ( Aimag(cc1p)==0 ) Then
 !			use the ieps instead
            n3 = 0
          Else
            n3 = nffet1(Cmplx(0._dp,ieps*Real(cy,dp),dp), 1/cc1p,cc2p/cc1p,ier)
          End If
        Else
          If ( Aimag(cc1p)==0 ) Then
            n3 =nffet1(cc2p,Cmplx(0._dp,-ieps*Real(cy1,dp),dp), cc2p/cc1p,ier)
          Else
            n3 = nffet1(cc2p,1/cc1p,cz,ier)
          End If
        End If
        ntot = n1-n2-n3
        crr(5) = (ntot*c2ipi + clo1)*clog1
        clog2p = zfflog(cc2p,ieps,cy)
        crr(6) = clo2*(n2*c2ipi - clog2p)
 !--#]		transform 1-x:
      Else If ( iclas1 == 3 ) Then
 !--#[		transform 1/x:
 !
 !		we transformed to 1/x for both dilogs
 !
        clog2p = zfflog(-cc2p,ieps,cy1)
        If ( Aimag(cc2p)==0 .Or. Aimag(cc1)==0 ) Then
 !		    we chose the eta's already equal, no worry.
          n3 = 0
          n3p = 0
        Else
          n3 = nffet1(-cc2p,-cc1,-cy/cy1,ier)
          n3p = nffet1(cc2p,cc1,-cy/cy1,ier)
        End If
        If ( n3/=0 .Or. n3p/=0 .Or. n1/=n2 ) Then
 !		    for the time being the normal terms, I'll have to think of
 !		    something smarter one day
          clog1p = zfflog(-cc1p,ieps,-cy)
          crr(5) = -clog1p**2/2
          crr(6) = +clog2p**2/2
          crr(7) = (n1*zfflog(cc1,ieps,cy) -n2*zfflog(cc2,ieps,-cy1))*c2ipi
        Else
          crr(5) = clo1*(n2*c2ipi + clog2p - clo1/2)
        End If
 !--#]		transform 1/x:
      End If
 ! -#]	    handle eta + transformation terms:
 ! -#[	    add up:
      If ( iclas1 == 1 .Or. iclas1 == 4 ) Then
        crr(1) = cli1
        crr(2) = cli2
        crr(3) = - cli3
        crr(4) = chill
      Else
        crr(1) = - cli1
        crr(2) = - cli2
        crr(3) = cli3
        crr(4) = - chill
      End If
 ! -#]	    add up:
  Else
 !	    Normal case:
 ! -#[	    handle dilogs:
 !
 !	    the dilogs will not come close together so just go on
 !	    only the special case cc1p ~ (-1,0) needs special attention
 !
      If ( iclas1 /= 4 .Or. .Not. ld2yzz ) Then
        Call ffzli2(cli1,clo1,cc1p)
      Else
        cd2 = cd2yzz + czz
        If ( absc(cd2) < xloss*absc(cd2yzz) ) Then
         cd2 = cy + cdyz
        End If
        cd2 = cd2/cdyz
        cfact = 1/(2-cd2)
        Call ffzli2(cli1,clo1,cd2*cfact)
        Call ffzli2(cli3,clo3,-cd2*cfact)
        Call ffzli2(cli4,clo4,cd2)
      End If
      If ( iclas2 /= 4 .Or. .Not. ld2yzz ) Then
        Call ffzli2(cli2,clo2,cc2p)
      Else
        If ( iclas1 == 4 )   Call WriteLFerror(26)
         cd2 = cd2yzz - czz1
        If ( absc(cd2) < xloss*absc(cd2yzz) ) Then
          cd2 = cdyz - cy1
        End If
        cd2 = cd2/cdyz
        cfact = 1/(2-cd2)
        Call ffzli2(cli2,clo2,cd2*cfact)
        Call ffzli2(cli3,clo3,-cd2*cfact)
        Call ffzli2(cli4,clo4,cd2)
      End If
 ! -#]	    handle dilogs:
 ! -#[	    handle eta terms:
 !
 !	    the eta's
 !
      If ( n1 /= 0 ) Then
        If ( iclas1 /= 2 .Or. absc(cc1p) > xloss ) Then
          If ( Real(cc1,dp) > -Abs(Aimag(cc1)) ) Then
            clog1 = zfflog(cc1,ieps,cy)
          Else
 !			take apart the factor i*pi^2
            If ( iclas1 == 4 ) Then
              clog1 = Log1minusX(cd2)
            Else
              clog1 = zfflog(-cc1,0,cy)
            End If
            If ( Aimag(cc1) < 0 ) Then
              i2pi = -1
            Else If ( Aimag(cc1) > 0 ) Then
              i2pi = +1
            Else If ( Real(cy)*ieps < 0 ) Then
              i2pi = -1
            Else If ( Real(cy)*ieps > 0 ) Then
              i2pi = +1
            Else
              Call WriteLFerror(18)
              i2pi = 0
            End If
            ipi12 = ipi12 - n1*24*i2pi
          End If
        Else
          clog1 = Log1minusX(cc1p)
        End If
        crr(5) = n1*c2ipi*clog1
 !	      Else
 !		crr(5) = 0
      End If
      If ( n2 /= 0 ) Then
        If ( iclas2 /= 2 .Or. absc(cc2p) > xloss ) Then
          If ( Real(cc2,dp) > -Abs(Aimag(cc2)) ) Then
            clog2 = zfflog(cc2,ieps,cy)
          Else
 !			take apart the factor i*pi^2
            If ( iclas2 == 4 ) Then
              clog2 = Log1minusX(cd2)
            Else
              clog2 = zfflog(-cc2,0,czero)
            End If
            If ( Aimag(cc2) < 0 ) Then
              i2pi = -1
            Else If ( Aimag(cc2) > 0 ) Then
              i2pi = +1
            Else If ( Real(cy)*ieps < 0 ) Then
              i2pi = -1
            Else If ( Real(cy)*ieps > 0 ) Then
              i2pi = +1
            Else
              Call WriteLFerror(18)
              i2pi = 0
            End If
            ipi12 = ipi12 + n2*24*i2pi
          End If
        Else
          clog2 = Log1minusX(cc2p)
        End If
        crr(6) = n2*c2ipi*clog2
 !	      Else
 !		crr(6) = 0
      End If
 ! -#]	    handle eta terms:
 ! -#[	    handle transformation terms:
 !
 !	    transformation of cc1
 !
      If ( iclas1 == 1 ) Then
 !		crr(3) = 0
      Else If( iclas1 == 2 ) Then
        cli1 = -cli1
        ipi12 = ipi12 + 2
        crr(3) = - clo1*zfflog(cc1p,ieps,cy)
      Else If ( iclas1 == 3 ) Then
        cli1 = -cli1
        ipi12 = ipi12 - 2
        clog1p = zfflog(-cc1p,ieps,cy1)
        crr(3) = - clog1p**2/2
      Else If ( iclas1 == 4 ) Then
 !		Note that this sum   Does not cause problems as d2<<1
        crr(3) = -cli3 - cli4 + clo4*zfflog(cfact,0,czero)
        ipi12 = ipi12 - 1
      Else
        Call WriteLFerror(27)
      End If
 !
 !	    transformation of cc2
 !
      If ( iclas2 == 1 ) Then
      Else If( iclas2 == 2 ) Then
        cli2 = -cli2
        ipi12 = ipi12 - 2
        crr(4) = clo2*zfflog(cc2p,ieps,cy)
      Else If ( iclas2 == 3 ) Then
        cli2 = -cli2
        ipi12 = ipi12 + 2
        clog2p = zfflog(-cc2p,ieps,cy1)
        crr(4) = clog2p**2/2
      Else If ( iclas2 == 4 ) Then
 !		Note that this sum   Does not cause problems as d2<<1
        crr(4) = cli3 + cli4 - clo4*zfflog(cfact,0,czero)
        ipi12 = ipi12 + 1
      Else
        Call WriteLFerror(28)
      End If
 ! -#]	    handle transformation terms:
 ! -#[	    sum:
      crr(1) = cli1
      crr(2) = - cli2
      crr(6) = - crr(6)
 !	    crr(7) = 0
 ! -#]	    sum:
  End If
 !  #] calculations:
  End Subroutine ffcrr

 Subroutine ffcs3(cs3,ipi12,cy,cz,cdyz,cd2yzz,cpi,cpiDpj,ii,ns, isoort)
 !--------------------------------------------------------------------
 !	calculates the s3 as defined in appendix b.
 !
 !		  log( cpi(ii+3)*y^2 + (cpi(ii+3)+cpi(ii)-cpi(ii+1))*y
 !	     /1 		     +	cpi(ii+1))  - log( ... ) |y=cyi
 !	s3 = \ dy ----------------------------------------------------
 !	     /0 			y - cyi	
 !									*
 !	   = r(cyi,cy+) + r(cyi,cy-) +  ( eta(-cy-,-cy+) -
 !		eta(1-cy-,1-cy+) - eta(...) )*log(1-1/cyi)
 !
 !	with y+- the roots of the argument of the logarithm.
 !
 !	input:	cy(4)	 (complex)  cy(1)=y^-,cy(2)=y^+,cy(i+2)=1-cy(1)
 !		cz(4)	 (complex)  cz(1)=z^-,cz(2)=z^+,cz(i+2)=1-cz(1)
 !		cpi(6)   (complex)  masses & momenta (B&D)
 !		ii	 (  Integer)  position of cp,cma,cmb in cpi
 !		ns	 (  Integer)  size of arrays
 !		isoort(2)(  Integer)  returns the kind of action taken
 !		cs3	 (complex)(14)	assumed zero.
 !
 !	output: cs3	 (complex)  mod factors ipi12
 !		ipi12(2) (  Integer)  these factors
 !		ier	 (  Integer)  0=ok, 1=numerical problems, 2=error
 !
 !	  Calls:	ffcrr,Aimag,Real,zfflog
 !--------------------------------------------------------------------
 Implicit None
  Integer :: ipi12(2),ii,ns,isoort(2),ier
  Complex(dp) :: cs3(20),cpi(ns),cpiDpj(ns,ns), cy(4),cz(4),cdyz(2,2),cd2yzz

  Integer :: i,ip,ieps(2),ieps0,ni(4),ntot
  Logical :: ld2yzz
  Complex(dp) :: c,cyy,cyy1,czz,czz1,cdyyzz,zdilog
  Real(dp) :: y,y1,z,z1,dyz,d2yzz,zz,zz1, x00(3),sprec

 !  #[ get ieps:
  ip = ii+3
  Call ffieps(ieps,cz(1),cpi(ip),cpiDpj(ip,ii),isoort)
 !  #] get ieps:
 !  #[ special case |cz| >> |cy|:
  If ( isoort(2) /= 0 .And. Max(absc(cy(2)),absc(cy(4))) <&
    &    xloss*Min(absc(cz(1)),absc(cz(2)))/2 ) Then
 !
 !	    we will obtain cancellations of the type Li_2(x) + Li_2(-x)
 !	    with x small.
 !
    cyy = cdyz(2,1)/cd2yzz
    cyy1 = cdyz(2,2)/cd2yzz
    If ( absc(cy(2)) < xclogm ) Then
      If ( Aimag(cy(2)) == 0 .And. Abs(Real(cy(2),dp)) >xalogm ) Then
        czz = cz(2)*cyy*Cmplx(1/Real(cy(2),dp),0._dp,dp)
        cdyyzz = cyy*cdyz(2,2)*Cmplx(1/Real(cy(2),dp),0._dp,dp)
      Else If ( cy(2) == 0 .And. cz(2) /= 0 .And. cyy/= 0 ) Then
 !		    the answer IS zero
        Goto 30
      End If
    Else
      czz = cz(2)*cyy/cy(2)
      cdyyzz = cyy*cdyz(2,2)/cy(2)
    End If
    czz1 = 1-czz
    If ( isoort(1) == -10 ) Then
 !		no eta terms.
      ieps0 = 99
    Else
 !		  Do not know the im part
      ieps0 = 0
    End If
    Call ffcrr(cs3(1),ipi12(1),cyy,cyy1,czz,czz1,cdyyzz,.False., &
     &  czero,czero,czero,-1,ieps0,ier)

  30 If ( absc(cy(4)) < xclogm ) Then
      If ( Aimag(cy(4)) == 0 .And. Abs(Real(cy(4),dp)) >xalogm ) Then
        czz = cz(4)*cyy*Cmplx(1/Real(cy(4),dp),0._dp,dp)
        cdyyzz = -cyy*cdyz(2,2)*Cmplx(1/Real(cy(4),dp),0._dp,dp)
      Else If ( cy(4) == 0 .And. cz(4) /= 0 .And. cyy/= 0 ) Then
 !		    the answer IS zero
        Goto 50
      End If
     Else
      czz = cz(4)*cyy/cy(4)
      cdyyzz = -cyy*cdyz(2,2)/cy(4)
     End If
     czz1 = 1-czz
     Call ffcrr(cs3(8),ipi12(2),cyy,cyy1,czz,czz1,cdyyzz,.False., &
     czero,czero,czero,-1,ieps0,ier)
     Do i=8,14
      cs3(i) = -cs3(i)
     End Do
 !
 !	    And now the remaining Li_2(x^2) terms
 !	    stupid Gould NP1
 !
  50  c = cy(2)*cy(2)/(cdyz(2,1)*cdyz(2,1))
     zdilog = CLi2(c)
     cs3(15) = +zdilog/2
 !	    stupid Gould NP1
     c = cy(4)*cy(4)/(cdyz(2,1)*cdyz(2,1))
     zdilog = CLi2(c)
     cs3(16) = -zdilog/2
     Return
  End If
 !  #] special case |cz| >> |cy|:
 !  #[ normal:
  If ( isoort(2) == 0 ) Then
    ld2yzz = .False.
  Else
    ld2yzz = .True.
  End If
  If ( isoort(1) == 0 ) Then
 !	      Do nothing
  Else If ( Mod(isoort(1),10)==0 .Or. Mod(isoort(1),10)==-1&
      &  .Or. Mod(isoort(1),10)==-3 ) Then
   Call ffcrr(cs3(1),ipi12(1),cy(2),cy(4),cz(1),cz(3), &
             cdyz(2,1),ld2yzz,cd2yzz,cz(2),cz(4),isoort(1), ieps(1),ier)
  Else If ( Mod(isoort(1),10) == -5 .Or. Mod(isoort(1),10) == -6 ) Then
    y = Real(cy(2),dp)
    y1 = Real(cy(4),dp)
    z = Real(cz(1),dp)
    z1 = Real(cz(3),dp)
    dyz = Real(cdyz(2,1),dp)
    d2yzz = Real(cd2yzz,dp)
    zz = Real(cz(2),dp)
    zz1 = Real(cz(4),dp)
    sprec = precx
    precx = precc
      Call ffcxr(cs3(1),ipi12(1),y,y1,z,z1,dyz,ld2yzz,d2yzz,zz,zz1 &
      ,.False.,x00,ieps(1),ier)
    precx = sprec
  Else
      Call WriteLFerror(20)
  End If
  If ( isoort(2) == 0 ) Then
 !	      Do nothing
  Else If ( Mod(isoort(2),5) == 0 ) Then
    Do i=1,7
      cs3(i) = 2*Real(cs3(i),dp)
    End Do
    ipi12(1) = 2*ipi12(1)
  Else If ( Mod(isoort(2),10)==-1 .Or. Mod(isoort(1),10)==-3 ) Then
      Call ffcrr(cs3(8),ipi12(2),cy(2),cy(4),cz(2),cz(4), &
      cdyz(2,2),ld2yzz,cd2yzz,cz(1),cz(3),isoort(2), ieps(2),ier)
  Else If ( Mod(isoort(2),10) == -6 ) Then
    y = Real(cy(2),dp)
    y1 = Real(cy(4),dp)
    z = Real(cz(2),dp)
    z1 = Real(cz(4),dp)
    dyz = Real(cdyz(2,2),dp)
    d2yzz = Real(cd2yzz,dp)
    zz = Real(cz(1),dp)
    zz1 = Real(cz(3),dp)
    sprec = precx
    precx = precc
      Call ffcxr(cs3(8),ipi12(2),y,y1,z,z1,dyz,ld2yzz,d2yzz,zz,zz1 &
      ,.False.,x00,ieps(2),ier)
    precx = sprec
  Else
      Call WriteLFerror(21)
  End If
 !  #] normal:
 !  #[ eta's:
  If ( Mod(isoort(1),10)==-5 .Or. Mod(isoort(1),10)==-6 ) Then
      If ( Mod(isoort(2),10)/=-5 .And. Mod(isoort(1),10)/=-6 ) Then
      Print *,'ffcxs3: error: I assumed both would be real!'
      ier = ier + 50
      End If
 !	    we   Called ffcxr - no eta's
  Else If ( Aimag(cpi(ip))==0 ) Then
    Call ffgeta(ni,cz(1),cdyz(1,1), cpi(ip),cpiDpj(ii,ip),ieps,isoort,ier)
    ntot = ni(1) + ni(2) + ni(3) + ni(4)
    If ( ntot /= 0 ) Call ffclgy(cs3(15),ipi12(2),ntot, cy(1),cz(1),cd2yzz,ier)
  Else
 !
 !	    cpi(ip) is really complex (occurs in transformed
 !	    4pointfunction)
 !
    Print *,'THIS PART IS NOT READY ', 'and should not be reached'
    Call TerminateProgram()
  End If
 !  #] eta's:

 End Subroutine ffcs3

 Subroutine ffcs4(cs3,ipi12,cw,cy,cz,cdwy,cdwz,cdyz,cd2yww,cd2yzz &
                 &,cpi,cpiDpj,cp2p,ii,ns,isoort)
 !------------------------------------------------------------------------
 !	Calculate the 8 Spence functions = 4 R's = 2 dR's
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ipi12(4),ii,ns,isoort(4),ier
  Complex(dp) :: cs3(40)
  Complex(dp) :: cw(4),cy(4),cz(4),cdwy(2,2),cdwz(2,2),cdyz(2,2)
  Complex(dp) :: cd2yww,cd2yzz,cpi(ns),cp2p,cpiDpj(ns,ns)

  Logical :: ld2yzz,ld2yww
  Integer :: i,j,ip,iepz(2),iepw(2),nz(4),nw(4),ntot,i2pi
  Complex(dp) :: c,cc,clogy,c2y1,cdyw(2,2)

 !  #[ get counters:
  ip = ii+3
  If ( isoort(2) == 0 ) Then
    ld2yzz = .False.
  Else
    ld2yzz = .True.
  End If
  If ( isoort(4) == 0 ) Then
    ld2yww = .False.
  Else
    ld2yww = .True.
  End If
  Call ffieps(iepz,cz,cpi(ip),cpiDpj(ip,ii),isoort)
  Call ffieps(iepw,cw,cp2p,cpiDpj(ip,ii),isoort(3))
  If ( isoort(4) == 0 ) Then
    Print *,'ffcs4: error: case not implemented'
    ier = ier + 50
  End If
 !  #] get counters: 
 !  #[ R's:
  If ( isoort(4) == 0 ) Then
      Call ffcrr(cs3(1),ipi12(1),cy(2),cy(4),cz(1),cz(3),cdyz(2,1) &
                & ,ld2yzz,cd2yzz,cz(2),cz(4),isoort(4),iepz(1),ier)
  Else
      If ( .Not. ( cdwz(2,1)==0 .And. iepz(1)==iepw(2) ) ) &
        Call ffdcrr(cs3( 1),ipi12(1),cy(2),cy(4),cz(1),cz(3),cz(2), &
                   cz(4),cd2yzz,cw(2),cw(4),cw(1),cw(3),cd2yww,cdyz(2,1), &
                   cdwy(2,2),cdwz(2,1),isoort(4),iepz(1),iepw(2),ier)
  End If
  If ( isoort(2) == 0 ) Then
      Call ffcrr(cs3(1),ipi12(1),cy(2),cy(4),cw(1),cw(3),-cdwy(1,2) &
                ,ld2yww,cd2yww,cw(2),cw(4),isoort(2),iepw(1),ier)
  Else
      If ( .Not. ( cdwz(1,2)==0 .And. iepz(2)==iepw(1) ) ) &
        Call ffdcrr(cs3(21),ipi12(3),cy(2),cy(4),cz(2),cz(4),cz(1), &
                    cz(3),cd2yzz,cw(1),cw(3),cw(2),cw(4),cd2yww,cdyz(2,2), &
                     cdwy(1,2),cdwz(1,2),iepz(2),isoort(2),iepw(1),ier)
  End If
 !  #] R's: 
 !  #[ eta's:
  If ( Aimag(cpi(ip)) == 0 ) Then
      Call ffgeta(nz,cz,cdyz, cpi(ip),cpiDpj(ii,ip),iepz,isoort,ier)
      Do i=1,2
        Do j=1,2
          cdyw(i,j) = cdwy(j,i)
        End Do
      End Do
      Call ffgeta(nw,cw,cdyw, cp2p,cpiDpj(ii,ip),iepw,isoort(3),ier)
  Else
    Print *,'ffcs4: error: not ready for complex D0 yet'
  End If
  ntot = nz(1)+nz(2)+nz(3)+nz(4)-nw(1)-nw(2)-nw(3)-nw(4)
  If ( ntot /= 0 ) Then
    i2pi = 0
      If ( 1/absc(cy(2)) < xloss ) Then
        clogy = Log1MinusX(1/cy(2))
      Else
        c = -cy(4)/cy(2)
        If ( Real(c,dp) > -Abs(Aimag(c)) ) Then
          clogy = zfflog(c,0,czero)
        Else
 !		    take out the factor 2*pi^2
          cc = c+1
          If ( absc(cc) < xloss ) Then
            c2y1 = -cd2yzz - cz(1) + cz(4)
            If ( absc(c2y1) < xloss*Max(absc(cz(1)), absc(cz(4))) ) Then
              c2y1 = -cd2yzz - cz(2) + cz(3)
            End If
            clogy = Log1MinusX(-c2y1/cy(2))
          Else
            clogy = zfflog(-c,0,czero)
          End If
          If ( Aimag(c) < 0 ) Then
            i2pi = -1
          Else If ( Aimag(c) > 0 ) Then
            i2pi = +1
          Else
            Call WriteLFerror(18)
            i2pi = 0
          End If
          ipi12(2) = ipi12(2) - ntot*24*i2pi
        End If
      End If
      If ( cs3(40) /= 0 ) Print *,'ffcs4: error: cs3(40) != 0'
      cs3(40) = ntot*c2ipi*clogy
  End If
 !  #] eta's: 
  End Subroutine ffcs4



 Subroutine ffcxr(crr,ipi12,y,y1,z,z1,dyz,ld2yzz,d2yzz,zz,zz1, &
                & ldy2z,dy2z,ieps,ier)
 !------------------------------------------------------------------------
 !	calculates R as defined in appendix b:
 !
 !		   /1    log(x-z+i*eps) - log(y-z+i*eps)
 !	r(y,z)  =  \ dx  -----------------------------------
 !		   /0		      x-y
 !
 !	    = li2(y/(y-z)+i*eps') - li2((y-1)/(y-z)+i*eps')
 !
 !	y,z are real, ieps  integer denoting the sign of i*eps.
 !	factors pi^2/12 are passed in the   Integer ipi12.
 !
 !	Input:	y	(real)
 !		y1	(real)		1-y
 !		z	(real)
 !		z1	(real)		1-z
 !		dyz	(real)		y-z
 !
 !		ld2yzz	(logical)	  If .TRUE. also defined are:
 !		d2yzz	(real)		2*y - z^+ - z^-
 !		zz	(real)		the other z-root
 !		zz1	(real)		1 - zz
 !
 !		ieps	(  Integer)	if +/-1 denotes sign imaginary
 !					part of	argument logs
 !		ieps	(  Integer)	if +/-2 denotes sign imaginary
 !					part of	argument dilogs
 !
 !	Output	crr	(complex)	R modulo factors pi^2/12
 !		ipi12	(  Integer)	these factors
 !		ier	(intger)	0=ok, 1=num prob, 2=error
 !
 !	Calls:	ffxli2,(test: ffzxdl),dfflo1,zxfflg
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ipi12,ieps,ier
  Logical :: ld2yzz,ldy2z
  Real(dp) :: y,y1,z,z1,dyz,d2yzz,zz,zz1,dy2z(3)
  Complex(dp) :: crr(7)

  Integer :: i,iclas1,iclas2
  Real(dp) :: fact,xx1,xx2,xx1p,xx2p,arg2,arg3, &
    xli1,xli2,xli3,xlo1,xlo2,xlo3,xhill,xlog1, &
    xlog2p,xx1n,d2,d21,d2n,d21n1,term,tot,xtroep,xli4, xlo4,som,xmax
  Complex(dp) :: clog1p,clog2p
 !  #[ groundwork:

 !
  If ( dyz == 0 ) Return
  fact = 1/dyz
  xx1 = y * fact
  xx2 = - y1 * fact
 !
 !  #] groundwork:
 !  #[ which area?:
 !
 !	determine the area:	1 = [-1+xloss,1/2]
 !				2 = (1/2,2-xloss]
 !				3 = [2+xloss,->) U (<-,-1-xloss]
 !				4 = [-1-xloss,-1+xloss]
 !				5 = [2-xloss,2+xloss]
 !
  If ( xx1 < -1-xloss/2 ) Then
    iclas1 = 3
    xx1p = 1/xx1
  Else If( xx1 < -1+xloss/2 ) Then
    If ( ld2yzz ) Then
      iclas1 = 4
    Else
      iclas1 = 1
    End If
    xx1p = xx1
  Else If( xx1 <= 0.5_dp ) Then
    iclas1 = 1
    xx1p = xx1
  Else If ( xx1 < 2-xloss ) Then
    iclas1 = 2
    xx1p = -z*fact
  Else If ( ldy2z .And. xx1 < 2+xloss ) Then
    iclas1 = 5
    xx1p = dy2z(1)*fact
  Else
    iclas1 = 3
    xx1p = 1/xx1
  End If
  If ( xx2 < -1-xloss/2 ) Then
    iclas2 = 3
    xx2p = 1/xx2
  Else If( xx2 < -1+xloss/2 ) Then
    If ( ld2yzz ) Then
      iclas2 = 4
    Else
      iclas2 = 1
    End If
    xx2p = xx2
  Else If ( xx2 <= 0.5_dp ) Then
    iclas2 = 1
    xx2p = xx2
  Else If ( xx2 < 2-xloss ) Then
    iclas2 = 2
    xx2p = z1*fact
  Else If ( ldy2z .And. xx2 < 2+xloss ) Then
    iclas2 = 5
    xx2p = -dy2z(3)*fact
  Else
    iclas2 = 3
    xx2p = 1/xx2
  End If
 !
 !	throw together   If they are close
 !
  If ( iclas1 /= iclas2 .And. Abs(xx1-xx2) < 2*xloss ) Then
 !	    we   Don't want trouble with iclasn = 4,5
    If ( iclas1 == 4 ) Then
      iclas1 = 1
    Else If ( iclas1 == 5 ) Then
      iclas1 = 3
      xx1p = 1/xx1
    End If
    If ( iclas2 == 4 ) Then
      iclas2 = 1
    Else If ( iclas2 == 5 ) Then
      iclas2 = 3
      xx2p = 1/xx2
    End If
    If ( iclas1 == iclas2 ) Goto 5
 !	    go on
      If ( iclas1 <= iclas2 ) Then
      iclas2 = iclas1
      If ( iclas1 == 1 ) Then
        xx2p = xx2
      Else
        xx2p = z1*fact
      End If
      Else
      iclas1 = iclas2
      If ( iclas1 == 1 ) Then
        xx1p = xx1
      Else
        xx1p = -z*fact
      End If
    End If
  End If
 !  #] which area?:
 !  #[ calculations:
5   If ( iclas1 == iclas2 .And.Abs(xx1p-xx2p) < &
       &  2*xloss*Max(Abs(xx1p),Abs(xx2p)) .And. iclas1 /= 5 ) Then
 !		      |----->temporary!
 !	    Close together:
 ! -#[	    handle dilog's:
      If ( Abs(xx2p) > xloss ) Then
 !--#[		Hill identity:
 !
 !		Use the Hill identity to get rid of the cancellations.
 !
 !
 !
        If ( iclas1 == 1 .Or. iclas1 == 4 ) Then
          d2 = 1/y
          arg2 = 1/z1
          arg3 = arg2/xx1p
        Else If ( iclas1 == 2 ) Then
          d2 = 1/z
          arg2 = 1/y1
          arg3 = arg2/xx1p
        Else If ( iclas1 == 3 ) Then
          d2 = 1/y1
          arg3 = 1/z1
          arg2 = arg3*xx1p
        End If
        Call ffxli2(xli1,xlo1,d2)
        Call ffxli2(xli2,xlo2,arg2)
        Call ffxli2(xli3,xlo3,arg3)
        If ( Abs(xx2p) < xloss ) Then
          xlog2p = Log1minusX(xx2p)
        Else
          xlog2p = Real(zxfflg(1-xx2p,0,1._DP),dp)
        End If
        xhill = xlo1*xlog2p
 !--#]		Hill identity:
      Else
 !--#[		Taylor expansion:
 !
 !		  If the points are close to zero   Do a Taylor
 !		expansion of the first and last dilogarithm
 !
 !			Li2(xx1p) - Li2(xx2p)
 !			  = sum xx1p^i ( 1-(1-d2)^i ) /i^2
 !
 !		with d2 = 1-xx2p/xx1p = ...
 !
        If ( iclas1 == 1 .Or. iclas1 == 4 ) Then
          d2 = 1/y
        Else If ( iclas1 == 2 ) Then
          d2 = 1/z
        Else If ( iclas1 == 3 ) Then
          d2 = 1/y1
        End If
 !		flag to the print section that we did a Taylor expansion
        d21 = 1-d2
        d21n1 = 1
        xx1n = xx1p
        d2n = d2
        tot = xx1p*d2
 !		check for possible underflow on the next line
        If ( Abs(xx1p) .Ge. xalog2 ) Then
          Do i=2,20
           xx1n = xx1n*xx1p
           d21n1 = d21n1*d21
           d2n = d2n + d2*d21n1
           term = xx1n*d2n*xn2inv(i)
           tot = tot + term
           If ( Abs(term) <= precx*Abs(tot) ) Exit
          End Do
        End If
        xli1 = tot
        xli2 = 0
        xli3 = 0
        xhill = 0
 !		for the eta+transformation section we also need
        If ( iclas1 /= 1 ) Then
          If ( Abs(d2) < xloss ) Then
            xlo1 = Log1minusX(d2)
          Else
            xlo1 = Real(zxfflg(d21,0,1._DP),dp)
          End If
        End If
        If ( iclas1 == 2 ) xlo2 = Log1minusX(1/y1)
 !--#]		Taylor expansion:
      End If
 !
 ! -#]	    handle dilog's:
 ! -#[	    handle transformation terms:
      If ( iclas1 == 1 .Or. iclas1 == 4 ) Then
 !
 !		no transformation was made.
 !
 !		crr(5) = 0
 !		crr(6) = 0
      Else If ( iclas1 == 2 ) Then
 !
 !		we tranformed to 1-x for both dilogs
 !
        If ( Abs(xx1p) < xloss ) Then
          xlog1 = Log1minusX(xx1p)
        Else
          xlog1 = Real(zxfflg(xx1,0,1._DP),dp)
        End If
        crr(5) = xlo1*xlog1
        clog2p = zxfflg(xx2p,ieps,-y1)
        crr(6) = -Real(xlo2,dp)*clog2p
      Else If ( iclas1 == 3 ) Then
 !
 !		we transformed to 1/x for both dilogs
 !
        clog2p = zxfflg(-xx2p,-ieps,-y1)
        crr(5) = Real(xlo1,dp)*(clog2p - Real(xlo1,dp)/2)
      End If
 ! -#]	    handle transformation terms:
 ! -#[	    add up and print out:
      If ( iclas1 == 1 .Or. iclas1 == 4 ) Then
        crr(1) = xli1
        crr(2) = xli2
        crr(3) = - xli3
        crr(4) = xhill
      Else
        crr(1) = - xli1
        crr(2) = - xli2
        crr(3) = xli3
        crr(4) = - xhill
      End If
 ! -#]	    add up and print out:
  Else
 !	    Normal case:
 ! -#[	    handle dilogs:
 !
 !	    the dilogs will not come close together so just go on
 !	    only the special case xx1p ~ -1 needs special attention
 !	    - and the special case xx1 ~ 2 also needs special attention
 !
      If ( iclas1 == 4 ) Then
        d2 = d2yzz + zz
        xmax = Abs(d2yzz)
        If ( Abs(d2) < xloss*xmax ) Then
          som = y + dyz
          If ( Abs(y)<xmax ) Then
            d2 = som
            xmax = Abs(y)
          End If
        End If
        d2 = d2/dyz
        fact = 1/(2-d2)
        Call ffxli2(xli1,xlo1,d2*fact)
        Call ffxli2(xli3,xlo3,-d2*fact)
        Call ffxli2(xli4,xlo4,d2)
      Else If ( iclas1 == 5 ) Then
        Call ffxl22(xli1,xx1p)
        ipi12 = ipi12 + 3
      Else
        Call ffxli2(xli1,xlo1,xx1p)
      End If
      If ( iclas2 == 4 ) Then
        If ( iclas1 == 4 )   Call WriteLFerror(26)
        d2 = d2yzz - zz1
        xmax = Abs(d2yzz)
        If ( Abs(d2) < xloss*xmax ) Then
          som = dyz - y1
          If ( Abs(y1)<xmax ) Then
            d2 = som
            xmax = Abs(y1)
          End If
        End If
        d2 = d2/dyz
        fact = 1/(2-d2)
        Call ffxli2(xli2,xlo2,d2*fact)
        Call ffxli2(xli3,xlo3,-d2*fact)
        Call ffxli2(xli4,xlo4,d2)
      Else If ( iclas2 == 5 ) Then
        Call ffxl22(xli2,xx2p)
        ipi12 = ipi12 - 3
      Else
        Call ffxli2(xli2,xlo2,xx2p)
      End If
 ! -#]	    handle dilogs:
 ! -#[	    handle transformation terms xx1:
 !
 !	    transformation of c1
 !
      If ( iclas1 == 1 ) Then
        crr(1) = xli1
      Else If( iclas1 == 2 ) Then
        crr(1) = -xli1
        ipi12 = ipi12 + 2
        clog1p = zxfflg(xx1p,ieps,y)
        crr(3) = - Real(xlo1,dp)*clog1p
      Else If ( iclas1 == 3 ) Then
        crr(1) = -xli1
        ipi12 = ipi12 - 2
        clog1p = zxfflg(-xx1p,-ieps,y)
        crr(3) = - clog1p**2/2
      Else If ( iclas1 == 4 ) Then
        crr(1) = xli1
 !		Note that this sum   Does not cause problems as d2<<1
        crr(3) = Real(-xli3-xli4,dp) + Real(xlo4,dp)*zxfflg(fact,0,0._DP)
        ipi12 = ipi12 - 1
      Else If ( iclas1 == 5 ) Then
        crr(1) = xli1
 !		supply an imaginary part
        clog1p = zxfflg(-1/xx1,-ieps,y)
        xtroep = -Aimag(clog1p)*Real(clog1p,dp)
        crr(3) = Cmplx(0._DP,xtroep,dp)
      Else
        Call WriteLFerror(26)
      End If
 ! -#]	    handle transformation terms xx1:
 ! -#[	    handle transformation terms xx2:
 !
 !	    transformation of c2
 !
      If ( iclas2 == 1 ) Then
        crr(2) = -xli2
      Else If( iclas2 == 2 ) Then
        crr(2) = +xli2
        ipi12 = ipi12 - 2
        clog2p = zxfflg(xx2p,ieps,-y1)
        crr(4) = + Real(xlo2,dp)*clog2p
      Else If ( iclas2 == 3 ) Then
        crr(2) = +xli2
        ipi12 = ipi12 + 2
        clog2p = zxfflg(-xx2p,-ieps,-y1)
        crr(4) = clog2p**2/2
      Else If ( iclas2 == 4 ) Then
        crr(2) = -xli2
 !		Note that this sum   Does not cause problems as d2<<1
        crr(4) = Real(xli3+xli4,dp) - Real(xlo4,dp)*zxfflg(fact,0,0._DP)
        ipi12 = ipi12 + 1
      Else If ( iclas2 == 5 ) Then
        crr(2) = -xli2
 !		supply an imaginary part
        clog2p = zxfflg(-1/xx2,-ieps,-y1)
        xtroep = Aimag(clog2p)*Real(clog2p,dp)
        crr(4) = Cmplx(0._DP,xtroep,dp)
      Else
        Call  WriteLFerror(30)
      End If
 ! -#]	    handle transformation terms xx2:
  End If
 !  #] calculations:
  End Subroutine ffcxr



 Subroutine ffcxs3(cs3,ipi12,y,z,dyz,d2yzz,dy2z,xpi,piDpj,ii,ns, isoort,ier)
 !--------------------------------------------------------------------
 !	calculates the s3 as defined in appendix b.
 !		(ip = ii+3, is1 = ii, is2 = ii+1)
 !
 !		  log( xk*y^2 + (-xk+xm1-xm2)*y + xm2 - i*eps )
 !	     /1 				  - log( ... ) |y=yi
 !	s3 = \ dy --------------------------------------------------
 !	     /0 			y - yi	
 !
 !	    = r(yi,y-,+) + r(yi,y+,-)
 !
 !	with y+- the roots of the argument of the logarithm.
 !	the sign of the argument to the logarithms in r is passed
 !	in ieps
 !
 !	input:	y(4),z(4)	(real)	roots in form (z-,z+,1-z-,1-z+)
 !		dyz(2,2),d2yzz,	(real)	y() - z(), y+ - z- - z+
 !		dy2z(4)		(real)	y() - 2z()
 !		xpi	(real(ns))	p(i).p(i) (B&D metric)	i=1,3
 !					m(i)^2 = si.si		i=4,6
 !		ii	(  Integer)	xk = xpi(ii+3) etc
 !		ns	(  Integer)	size of arrays
 !		isoort	(  Integer)	returns kind of action taken
 !		cs3	(complex)(20)	assumed zero.
 !		ccy	(complex)(3)	  If i0 != 0: complex y
 !
 !	output: cs3	(complex)	mod factors pi^2/12, in array
 !		ipi12	(  Integer)	these factors
 !		ier	(  Integer)	0=ok 1=inaccurate 2=error
 !
 !	  Calls:	ffcrr,ffcxr,real/Real,Cmplx,log,ffadd1,ffadd2,ffadd3
 !--------------------------------------------------------------------
 Implicit None
  Integer :: ipi12(2),ii,ns,isoort(2),ier
  Complex(dp) :: cs3(20)
  Real(dp) :: y(4),z(4),dyz(2,2),d2yzz,dy2z(4), xpi(ns),piDpj(ns,ns)

  Integer :: i,ip,ieps(2)
  Real(dp) :: yy,yy1,zz,zz1,dyyzz,xdilog,x00(3)
  Logical :: ld2yzz
 !
 !  #[ get counters:
  ip = ii+3
  If ( isoort(2) /= 0 ) Then
   If ( (z(2)>z(1) .Or.  z(1)==z(2) .And. z(4)<z(3) ) .Eqv. xpi(ip) > 0 ) Then
      ieps(1) = +1
      ieps(2) = -1
   Else
      ieps(1) = -1
      ieps(2) = +1
   End If
  Else
   If ( piDpj(ip,ii) > 0 ) Then
      ieps(1) = +1
   Else
      ieps(1) = -1
   End If
  End If
 !  #] get counters:
 !  #[ special case |z| >> |y|:
  If ( xpi(ip)<0 .And. Max(Abs(y(2)),Abs(y(4))) < &
    &  xloss*Min(Abs(z(1)), Abs(z(2)))/2 ) Then
 !
 !	    we will obtain cancellations of the type Li_2(x) + Li_2(-x)
 !	    with x small.
 !
    yy = dyz(2,1)/d2yzz
    yy1 = dyz(2,2)/d2yzz
    If ( y(2) .Ne. 0 ) Then
     zz = z(2)*yy/y(2)
     zz1 = 1-zz
     dyyzz = dyz(2,2)*yy/y(2)
     Call ffcxr(cs3(1),ipi12(1),yy,yy1,zz,zz1,dyyzz,.False., &
      0._DP,0._DP,0._DP,.False.,x00,0,ier)
    End If
    If ( y(4) .Ne. 0 ) Then
      zz = yy*z(4)/y(4)
      zz1 = 1-zz
      dyyzz = -yy*dyz(2,2)/y(4)
      Call ffcxr(cs3(8),ipi12(2),yy,yy1,zz,zz1,dyyzz,.False., &
                0._DP,0._DP,0._DP,.False.,x00,0,ier)
      Do i=8,14
        cs3(i) = -cs3(i)
      End Do
    End If
 !	    And now the remaining Li_2(x^2) terms
    xdilog = Li2((y(2)/dyz(2,1))**2)
    cs3(15) = +xdilog/2
    xdilog = Li2((y(4)/dyz(2,1))**2)
    cs3(16) = -xdilog/2
    Return
  End If
 !  #] special case |z| >> |y|:
 !  #[ normal:
  If ( xpi(ip) == 0 ) Then
    ld2yzz = .False.
  Else
    ld2yzz = .True.
  End If
  If ( isoort(1) /= 0 )  &
     Call ffcxr(cs3(1),ipi12(1),y(2),y(4),z(1),z(3),dyz(2,1),ld2yzz &
               & ,d2yzz,z(2),z(4),.True.,dy2z(1), ieps(1),ier)
  If ( isoort(2) /= 0 ) Then
      If ( Mod(isoort(2),10) == 2 ) Then
 !		both roots are equal: multiply by 2
        Do i=1,7
        cs3(i) = 2*Real(cs3(i),dp)
        End Do
      ipi12(1) = 2*ipi12(1)
      Else
        Call ffcxr(cs3(8),ipi12(2),y(2),y(4),z(2),z(4),dyz(2,2), &
        ld2yzz,d2yzz,z(1),z(3),.True.,dy2z(2),ieps(2),ier)
      End If
  End If
 !
 !  #] normal:

 End Subroutine ffcxs3



 Subroutine ffcxs4(cs3,ipi12,w,y,z,dwy,dwz,dyz,d2yww,d2yzz, &
                 & xpi,piDpj,ii,ns,isoort,ier)
 !------------------------------------------------------------------------
 !	Calculate the 8 Spence functions = 4 R's = 2 dR's
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ipi12(4),ii,ns,isoort(4),ier
  Complex(dp) :: cs3(40)
  Real(dp) :: w(4),y(4),z(4),dwy(2,2),dwz(2,2),dyz(2,2), &
    d2yww,d2yzz,xpi(ns),piDpj(ns,ns),x00(3)

  Integer :: iepz(2),iepw(2)
  Logical :: ld2yzz,ld2yww

 !  #[ groundwork:
  If ( isoort(2) == 0 ) Then
    ld2yzz = .False.
  Else
    ld2yzz = .True.
  End If
  If ( isoort(4) == 0 ) Then
    ld2yww = .False.
  Else
    ld2yww = .True.
  End If
  If ( isoort(2) /= 0 ) Then
    If ( z(2) > z(1) .Eqv. xpi(ii+3) > 0 ) Then
      iepz(1) = +1
      iepz(2) = -1
    Else
      iepz(1) = -1
      iepz(2) = +1
    End If
  Else
    Print *,'ffcxs4: error: untested algorithm'
    If ( piDpj(ii,ii+3) > 0 ) Then
      iepz(1) = +1
    Else
      iepz(1) = -1
    End If
  End If
  If ( isoort(4) /= 0 ) Then
    If ( w(2) > w(1) .Eqv. xpi(5) > 0 ) Then
      iepw(1) = 1
      iepw(2) = -1
    Else
      iepw(1) = -1
      iepw(2) = 1
    End If
  Else
    Print *,'ffcxs4: error: untested algorithm'
    If ( piDpj(2,5) > 0 ) Then
      iepw(1) = +1
    Else
      iepw(1) = -1
    End If
  End If
 !  #] groundwork: 
 !  #[ zm and wp:
  If ( isoort(4) == 0 ) Then
   Call ffcxr(cs3(1),ipi12(1),y(2),y(4),z(1),z(3),dyz(2,1), &
         &  ld2yzz,d2yzz,z(2),z(4),.False.,x00,iepz(1),ier)
  Else
   If ( .Not. ( dwz(2,1)==0 .And. iepz(1)==iepw(2) ) ) &
        Call ffdcxr(cs3( 1),ipi12(1),y(2),y(4),z(1),z(3), &
               z(2),z(4),d2yzz,w(2),w(4),w(1),w(3),d2yww, &
               dyz(2,1),dwy(2,2),dwz(2,1),iepz(1),iepw(2),ier)
  End If
 !  #] zm and wp: 
 !  #[ zp and wm:
  If ( isoort(2) == 0 ) Then
      Call ffcxr(cs3(1),ipi12(1),y(2),y(4),w(1),w(3),-dwy(1,2), &
                 & ld2yww,d2yww,w(2),w(4),.False.,x00,iepw(1),ier)
  Else
      If ( .Not. ( dwz(1,2)==0 .And. iepz(2)==iepw(1) ) ) &
        Call ffdcxr(cs3(21),ipi12(3),y(2),y(4),z(2),z(4), &
                    z(1),z(3),d2yzz,w(1),w(3),w(2),w(4),d2yww, &
                    dyz(2,2),dwy(1,2),dwz(1,2),iepz(2),iepw(1),ier)
  End If
 !  #] zp and wm: 
  End Subroutine ffcxs4



 Subroutine ffcxyz(cy,cz,cdyz,cd2yzz,ivert,sdelpp,sdelps, &
                &    etami,delps,xpi,piDpj,isoort,ldel2s,ns,ier)
 !------------------------------------------------------------------------
 !	calculate in a numerically stable way
 !
 !	cz(1,2) = (-p(ip1).p(is2) +/- sdelpp)/xpi(ip1)
 !	cy(1,2) = (-p(ip1).p(is2) +/- sdisc)/xpi(ip1)
 !			disc = slam1 + 4*eta*xpi(ip)/slam
 !
 !	cy(3,4) = 1-cy(1,2)
 !	cz(3.4) = 1-cz(1,2)
 !	cdyz(i,j) = cy(i) - cz(j)
 !
 !	Input:	ivert		(  Integer)	1,2 of 3
 !		sdelpp		(real)		sqrt(lam(p1,p2,p3))/2
 !		sdelps		(real)		sqrt(-lam(p,ma,mb))/2
 !		etalam		(real)		det(si.sj)/det(pi.pj)
 !		etami(6)	(real)		si.si - etalam
 !		xpi(ns)		(real)		standard
 !		piDpj(ns,ns)	(real)		standard
 !		ns		(  Integer)	dim of xpi,piDpj
 !
 !	Output:	cy(4),cz(4),cdyz(4,4)	(complex)	see above
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ivert,isoort(2),ns,ier
  Logical :: ldel2s
  Complex(dp) :: cy(4),cz(4),cdyz(2,2),cd2yzz
  Real(dp) :: sdelpp,sdelps,etami(6),delps,xpi(ns), piDpj(ns,ns)

  Integer :: ip1,is1,is2,is3
  Real(dp) :: y(4)
  Real(dp) :: disc,hulp

 !  #[ set up pointers:
  If ( ldel2s .And. ivert /= 1 ) Then
 !  #[ special case, get indices:
    If ( ivert==2 ) Then
      is1 = 2
      ip1 = 5
    Else
      is1 = 1
      ip1 = 6
    End If
    isoort(1) = -100
    isoort(2) = -100
 !  #] special case, get indices:
 !  #[ get cypm,czpm:
 !
 !	special case del2s = 0, hence the roots are not the real roots
 !	but z_2'' = (z_2'-1)/delta, z''_3 = -z'_3/delta
 !
    hulp = sdelps/xpi(3)
    disc = delps/sdelpp
    If ( ivert == 3 ) Then
      hulp = -hulp
      disc = -disc
    End If
    cz(1) = Cmplx(piDpj(is1,3)/xpi(3),-hulp,dp)
    cz(2) = Cmplx(piDpj(is1,3)/xpi(3),+hulp,dp)
    Call Roots(xpi(3),piDpj(is1,3),etami(is1),disc,y(1),y(2))
    cy(1) = y(1)
    cy(2) = y(2)
 !  #] get cypm,czpm:
 !  #[ get cypm1,czpm1:
    cz(3) = 1 - cz(1)
    cz(4) = 1 - cz(2)
    If ( absc(cz(3))<xloss .Or. absc(cz(4))<xloss ) Then
      If ( ivert==2 ) Then
        cz(3) =Cmplx(piDpj(ip1,3)/xpi(3),+hulp,dp)
        cz(4) =Cmplx(piDpj(ip1,3)/xpi(3),-hulp,dp)
      Else
        cz(3) =Cmplx(-piDpj(ip1,3)/xpi(3),+hulp,dp)
        cz(4) =Cmplx(-piDpj(ip1,3)/xpi(3),-hulp,dp)
      End If
    End If
    y(3) = 1 - y(1)
    y(4) = 1 - y(2)
    If ( Abs(y(3)) < xloss .Or. Abs(y(4)) < xloss ) Then
      If ( ivert == 2 ) Then
        Call Roots(xpi(3),piDpj(ip1,3),etami(ip1), disc,y(4),y(3))
      Else
        Call Roots(xpi(3),-piDpj(ip1,3),etami(ip1), disc,y(4),y(3))
      End If
    End If
    cy(3) = y(3)
    cy(4) = y(4)
 !  #] get cypm1,czpm1:
 !  #[ get cdypzp, cdypzm:
    cdyz(2,1) = Cmplx(disc/xpi(3),+hulp,dp)
    cdyz(2,2) = Cmplx(disc/xpi(3),-hulp,dp)
    cdyz(1,1) = -cdyz(2,2)
    cdyz(1,2) = -cdyz(2,1)
    cd2yzz = 2*disc/xpi(3)
 !  #] get cdypzp, cdypzm:

  Else
    is1 = ivert
    is2 = ivert+1
    If ( is2 == 4 ) is2 = 1
    is3 = ivert-1
    If ( is3 == 0 ) is3 = 3
    ip1 = is1 + 3
 !	ip2 = is2 + 3
 !	ip3 = is3 + 3
    isoort(1) = -10
    isoort(2) = -10
 !  #] set up pointers:
 !  #[ get cypm,czpm:
    hulp = sdelps/xpi(ip1)
    cz(1) = Cmplx(piDpj(ip1,is2)/xpi(ip1),-hulp,dp)
    cz(2) = Cmplx(piDpj(ip1,is2)/xpi(ip1),+hulp,dp)
    disc = delps/sdelpp
    Call roots(xpi(ip1),piDpj(ip1,is2),etami(is2),disc,y(1),y(2))
    cy(1) = y(1)
    cy(2) = y(2)
 !  #] get cypm,czpm:
 !  #[ get cypm1,czpm1:
    If ( xpi(is1) == xpi(is2) ) Then
      cy(4) = cy(1)
      cy(3) = cy(2)
      cz(4) = cz(1)
      cz(3) = cz(2)
    Else
      cz(3) = 1 - cz(1)
      cz(4) = 1 - cz(2)
      If ( absc(cz(3))<xloss .Or. absc(cz(4))<xloss ) Then
        cz(3) =Cmplx(-piDpj(ip1,is1)/xpi(ip1),+hulp,dp)
        cz(4) =Cmplx(-piDpj(ip1,is1)/xpi(ip1),-hulp,dp)
      End If
      y(3) = 1 - y(1)
      y(4) = 1 - y(2)
      If ( Abs(y(3)) < xloss .Or. Abs(y(4)) < xloss ) Then
        Call Roots(xpi(ip1),-piDpj(ip1,is1), etami(is1),disc,y(4),y(3))
      End If
      cy(3) = y(3)
      cy(4) = y(4)
    End If
 !  #] get cypm1,czpm1:
 !  #[ get cdypzp, cdypzm:
    cdyz(2,1) = Cmplx(disc/xpi(ip1),+hulp,dp)
    cdyz(2,2) = Cmplx(disc/xpi(ip1),-hulp,dp)
    cdyz(1,1) = -cdyz(2,2)
    cdyz(1,2) = -cdyz(2,1)
    cd2yzz = 2*disc/xpi(ip1)
   End If
 !  #] get cdypzp, cdypzm:
  End Subroutine ffcxyz

 Subroutine ffdcrr(cs3,ipi12,cy,cy1,cz,cz1,czp,czp1,cd2yzz,cw,cw1 &
                 & ,cwp,cwp1,cd2yww,cdyz,cdwy,cdwz,isoort,iepsz,iepsw,ier)
 !------------------------------------------------------------------------
 !	Calculate
 !
 !		R(cy,cz,iepsz) - R(cy,cw,iepsw)
 !
 !	Input:
 !		a = [yzw]	(real)		see definition
 !		a1 = 1 - a	(real)
 !		dab = a - b	(real)
 !		ieps[zw]	(  Integer)	sign of imaginary part
 !						of argument logarithm
 !		cs3(20)		(complex)	assumed zero
 !
 !	Output:
 !		cs3(20)		(complex)	the results, not added
 !		ipi12(2)	(  Integer)	factors pi^2/12
 !
 !	Calls:	ffcrr
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ipi12(2),isoort,iepsz,iepsw,ier
  Complex(dp) :: cs3(20)
  Complex(dp) :: cy,cz,czp,cw,cwp,cy1,cz1,czp1,cw1,cwp1, &
    cdyz,cdwy,cdwz,cd2yzz,cd2yww

  Integer :: i,ieps,ieps1,ieps2, n1,n2,n3,n4,n5,n6
  Logical :: ld2yyz
  Complex(dp) :: cyy,cyy1,czz,czz1,cdyyzz,chulp, &
    cc1,cdw,cc1n,cterm,ctot,cd2,cd3, cd21,cd31,cd2n,cd3n,cd21n1,cd31n1, &
    cc2,cfactz,cfactw,czzp,czzp1,cd2yyz

 !  #[ groundwork:
  If ( cdwz == 0 ) Then
    If ( Abs(Aimag(cz)) > precc*Abs(Real(cz,dp)) .Or.iepsz == iepsw ) Return
    If ( Real(cz) >= 0 .And. Real(cz1) >= 0 ) Return
    Call WriteLfError(33)
    Return
  End If
  If ( cdyz == 0 ) Then
    Call WriteLfError(34)
    Return
  End If
  cc1 = cy/cdyz
  cdw = cdwz/cdyz
  If ( Real(cc1,dp) <= .5_DP .Or. Abs(cc1-1) > 1 ) Then
    cd2 = 1/cy
    cdw = cdw*cy/cw
  Else
    cd2 = 1/cz1
  End If
 !  #] groundwork:
 !  #[ trivial case:
  If ( absc(cdw) == 0 ) Then
 !  #] trivial case: 
 !  #[ normal case:
 !
 !	  If no cancellations are expected OR the imaginary signs differ
 !	and are sign  Ificant
 !
  Else If ( absc(cdw) > xloss .Or. (iepsz/=iepsw .And. &
      (Real(cy/cdyz,dp)>1._dp .Or. Real(-cy1/cdyz,dp)>1._dp) ) ) Then
 !	    nothing's the matter
 !	    special case to avoid bug found 15-oct=1995
      If ( iepsz==iepsw ) Then
        If ( Aimag(cz)==0 .And. Aimag(cz1)==0 ) Then
          Print *,'ffdcrr: flipping sign iepsz'
          iepsz = -iepsz
        Else If ( Aimag(cw)==0 .And. Aimag(cw1)==0 ) Then
          Print *,'ffdcrr: flipping sign iepsw'
          iepsw = -iepsw
        Else
          Print *,'ffdcrr: error: missing eta terms!'
          ier = ier + 100
        End If
      End If
      Call ffcrr(cs3(1),ipi12(1),cy,cy1,cz,cz1,cdyz,.True., &
                 cd2yzz,czp,czp1,isoort,iepsz,ier)
      Call ffcrr(cs3(8),ipi12(2),cy,cy1,cw,cw1,-cdwy,.True., &
                 cd2yww,cwp,cwp1,isoort,iepsw,ier)
      Do i=8,14
        cs3(i) = -cs3(i)
      End Do
      ipi12(2) = -ipi12(2)
 !  #] normal case:
 !  #[ only cancellations in cw, not in cy:
  Else If ( absc(cd2) > xloss ) Then
 !	    there are no cancellations the other way:
    cyy = cdwy/cdwz
    czz = cz*cyy/cy
    cyy1 = cdyz/cdwz
    czz1 = cyy1*cw/cy
    cdyyzz = cdyz*cyy/cy
    If ( Real(cy) > 0 ) Then
      ieps1 = -3*iepsz
    Else
      ieps1 = +3*iepsz
    End If
 !	    Often 2y-z-z is relevant, but 2*yy-zz-zz is not, solve by
 !	    introducing zzp.
    czzp = czp*cyy/cy
    cd2yyz = cd2yzz*cyy/cy
    czzp1 = 1 - czzp
    If ( absc(czzp1) < xloss ) Then
 !		later try more possibilities
      ld2yyz = .False.
    Else
      ld2yyz = .True.
    End If
    Call ffcrr(cs3(1),ipi12(1),cyy,cyy1,czz,czz1,cdyyzz, &
        ld2yyz,cd2yyz,czzp,czzp1,isoort,ieps1,ier)
    czz = cyy*cz1/cy1
    czz1 = cyy1*cw1/cy1
    If ( Real(-cy1) > 0 ) Then
      ieps2 = -3*iepsz
    Else
      ieps2 = +3*iepsz
    End If
    cdyyzz = -cyy*cdyz/cy1
    czzp = czp1*cyy/cy1
    cd2yyz = -cd2yzz*cyy/cy1
    czzp1 = 1 - czzp
    If ( absc(czzp1) < xloss ) Then
 !		later try more possibilities
      ld2yyz = .False.
    Else
      ld2yyz = .True.
    End If
    Call ffcrr(cs3(8),ipi12(2),cyy,cyy1,czz,czz1,cdyyzz, &
               .True.,cd2yyz,czzp,czzp1,isoort,ieps2,ier)
    Do i=8,14
      cs3(i) = -cs3(i)
    End Do
    ipi12(2) = -ipi12(2)
 !	    eta terms (are not calculated in ffcrr as ieps = 3)
    cfactz = 1/cdyz
    If ( Aimag(cz) == 0 ) Then
      If ( Aimag(cy) == 0 ) Then
        n1 = 0
        n2 = 0
      Else
        n1 = nffet1(Cmplx(0._dp,Real(iepsz,dp),dp),cfactz, -cz*cfactz,ier)
        n2 = nffet1(Cmplx(0._dp,Real(iepsz,dp),dp),cfactz, cz1*cfactz,ier)
      End If
    Else
      n1 = nffeta(-cz,cfactz,ier)
      n2 = nffeta(cz1,cfactz,ier)
    End If
    cfactw = -1/cdwy
    If ( Aimag(cw) == 0 ) Then
      If ( Aimag(cy) == 0 ) Then
        n4 = 0
        n5 = 0
      Else
        n4 = nffet1(Cmplx(0._dp,Real(iepsw,dp),dp),cfactw, -cw*cfactw,ier)
        n5 = nffet1(Cmplx(0._dp,Real(iepsw,dp),dp),cfactw, cw1*cfactw,ier)
      End If
    Else
      n4 = nffeta(-cw,cfactw,ier)
      n5 = nffeta(cw1,cfactw,ier)
    End If
 !
 !	    we assume that cs3(15-17) are not used, this is always true
 !
    n3 = 0
    n6 = 0
    If ( n1==n4 ) Then
      If ( n1==0 ) Then
 !		    nothing to   Do
      Else
        cc1 = cdwz/cdyz
        If ( absc(cc1) < xloss ) Then
          cs3(15) = n1*c2ipi*Log1MinusX(cc1)
        Else
          cc1 = -cdwy/cdyz
          cs3(15) = n1*c2ipi*zfflog(cc1,0,czero)
        End If
        cc1 = cy*cfactz
        cc2 = cy*cfactw
        If ( Aimag(cc1)==0 .Or. Aimag(cc2)==0 ) Then
          n3 = 0
        Else
          n3 = nffeta(cc1,1/cc2,ier)
        End If
        If ( n3/=0 ) Then
          Print *,'ffdcrr: error: untested algorithm'
          ier = ier + 50
          ipi12(1) = ipi12(1) + 4*12*n1*n3
        End If
      End If
    Else
      cc1 = cy*cfactz
      cc2 = cy*cfactw
      cs3(15) = (n1*zfflog(cc1,ieps1,czero) + n4*zfflog(cc2,ieps1,czero))*c2ipi
    End If
    If ( n2==n5 ) Then
      If ( n2==0 ) Then
 !		    nothing to   Do
      Else
        cc1 = cdwz/cdyz
        If ( absc(cc1) < xloss ) Then
          cs3(16) = n2*c2ipi*Log1MinusX(cc1)
        Else
          cc1 = -cdwy/cdyz
          cs3(16) = n2*c2ipi*zfflog(cc1,0,czero)
        End If
        cc1 = -cy1*cfactz
        cc2 = -cy1*cfactw
        If ( Aimag(cc1)==0 .Or. Aimag(cc2)==0 ) Then
          n6 = 0
        Else
          n6 = nffeta(cc1,1/cc2,ier)
        End If
        If ( n6/=0 ) Then
          Print *,'ffdcrr: error: untested algorithm'
          ier = ier + 50
          ipi12(2) = ipi12(2) + 4*12*n2*n6
        End If
      End If
    Else
      cc1 = -cy1*cfactz
      cc2 = -cy1*cfactw
      cs3(15) = (n2*zfflog(cc1,ieps2,czero) + n5*zfflog(cc2,ieps2,czero))*c2ipi
    End If
 !  #] only cancellations in cw, not in cy: 
 !  #[ Hill identity:
  Else If (  ( 1>xloss*absc(cy) .Or. absc(cc1)>xloss ) &
      .And. ( 1>xloss*absc(cz) .Or. absc(cz/cdyz)>xloss ) &
      .And. ( 1>xloss*absc(cy) .Or. absc(cdyz/cy)>xloss ) ) Then
 !	      Do a Hill identity on the cy,cy-1 direction
    cyy = -cy*cw1/cdwy
    cyy1 = cw*cy1/cdwy
    czz = -cz*cw1/cdwz
    czz1 = cw*cz1/cdwz
    cdyyzz = -cw*cw1*(cdyz/(cdwy*cdwz))
    ieps = -2*iepsz
    Call ffcrr(cs3(1),ipi12(1),cyy,cyy1,czz,czz1,cdyyzz, &
              .False.,czero,czero,czero,isoort,ieps,ier)
    cyy = cw1
    cyy1 = cw
    czz = -cw1*cz/cdwz
    czz1 = cw*cz1/cdwz
    cdyyzz = cw*cw1/cdwz
    Call ffcrr(cs3(8),ipi12(2),cyy,cyy1,czz,czz1,cdyyzz, &
               .False.,czero,czero,czero,isoort,0,ier)
    Do i=8,14
      cs3(i) = -cs3(i)
    End Do
    ipi12(2) = -ipi12(2)
 !	    the extra logarithms ...
    If ( 1 < xloss*absc(cw) ) Then
      chulp = Log1MinusX(1/cw)
    Else
      chulp = zfflog(-cw1/cw,0,czero)
    End If
    cs3(15) = -Log1MinusX(cdwz/cdwy)*chulp
 !  #] Hill identity: 
 !  #[ Taylor expansion:
  Else
 !	    Do a Taylor expansion
    If ( absc(cc1) < xloss ) Then
      cd3 = cdwz/cdwy
 !		isign = 1
      cc1n = cc1
      cd2n = cd2
      cd3n = cd3
      cd21 = 1-cd2
      cd21n1 = 1
      cd31 = 1-cd3
      cd31n1 = 1
      ctot = cc1*cd2*cd3
      Do i=2,20
        cc1n = cc1n*cc1
        cd21n1 = cd21n1*cd21
        cd31n1 = cd31n1*cd31
        cd2n = cd2n + cd2*cd21n1
        cd3n = cd3n + cd3*cd31n1
        cterm = cc1n*cd2n*cd3n*Real(xn2inv(i),dp)
        ctot = ctot + cterm
        If ( absc(cterm) < precc*absc(ctot) ) Exit
      End Do

      cs3(1) = ctot
    Else If ( absc(cz/cdyz) < xloss ) Then
      Call ffcrr(cs3(1),ipi12(1),cy,cy1,cz,cz1,cdyz,.True., &
          cd2yzz,czp,czp1,isoort,iepsz,ier)
      Call ffcrr(cs3(8),ipi12(2),cy,cy1,cw,cw1,-cdwy,.True., &
          cd2yww,cwp,cwp1,isoort,iepsw,ier)
      Do i=8,14
        cs3(i) = -cs3(i)
      End Do
      ipi12(2) = -ipi12(2)
    Else
        Call WriteLfError(35)
      Return
    End If
  End If
 !  #] Taylor expansion: 
  End Subroutine ffdcrr

 Subroutine ffdcxr(cs3,ipi12,y,y1,z,z1,zp,zp1,d2yzz, &
                 &  w,w1,wp,wp1,d2yww,dyz,dwy,dwz,iepsz,iepsw,ier)
 !------------------------------------------------------------------------
 !	Calculate
 !
 !		R(y,z,iepsz) - R(y,w,iepsw)
 !
 !	Input:
 !		a = [yzw]	(real)		see definition
 !		a1 = 1 - a	(real)
 !		dab = a - b	(real)
 !		ieps[zw]	(  Integer)	sign of imaginary part
 !						of argument logarithm
 !		cs3(20)		(complex)	assumed zero
 !
 !	Output:
 !		cs3(20)		(complex)	the results, not added
 !		ipi12(2)	(  Integer)	factors pi^2/12
 !
 !	Calls:	ffcxr
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ipi12(2),iepsz,iepsw,ier
  Complex(dp) :: cs3(20)
  Real(dp) :: y,z,w,y1,z1,w1,dyz,dwy,dwz,zp,zp1,d2yzz,wp,wp1, d2yww

  Integer :: i,ieps
  Logical :: again
  Real(dp) :: yy,yy1,zz,zz1,dyyzz,xx1,xx1n,term,tot,d2,d3, &
    d21,d31,d2n,d3n,d21n1,d31n1,dw,x00(3)
  Complex(dp) :: chulp

 !  #[ groundwork:
  If ( dwz == 0 .And. iepsz == iepsw ) Return
  If ( dyz == 0 ) Then
    Call WriteLfError(31)
    Return
  End If
  xx1 = y/dyz
  dw = dwz/dyz
  If ( xx1 <= .5_DP .Or. xx1 > 2 ) Then
    d2 = 1/y
    dw = dw*y/w
  Else
    d2 = 1/z1
  End If
  again = .False.
  123 Continue
 !  #] groundwork: 
 !  #[ trivial case:
  If ( dw == 0 ) Then
 !  #] trivial case: 
 !  #[ normal case:
  Else If ( Abs(dw) > xloss .Or. again ) Then
 !	    nothing's the matter
    Call ffcxr(cs3( 1),ipi12(1),y,y1,z,z1,dyz, &
      .True.,d2yzz,zp,zp1,.False.,x00,iepsz,ier)
    Call ffcxr(cs3(11),ipi12(2),y,y1,w,w1,-dwy, &
      .True.,d2yww,wp,wp1,.False.,x00,iepsw,ier)
    Do i=11,20
      cs3(i) = -cs3(i)
    End Do
    ipi12(2) = -ipi12(2)
 !  #] normal case: 
 !  #[ only cancellations in w, not in y:
  Else If ( Abs(d2) > xloss ) Then
 !	    there are no cancellations the other way:
    If ( iepsz /= iepsw .And. ( y/dyz > 1 .Or.-y/dwy >1 ) ) Then
      again = .True.
      Goto 123
    End If
    yy = dwy/dwz
    zz = yy*z/y
    yy1 = dyz/dwz
    zz1 = yy1*w/y
    dyyzz = yy*dyz/y
    If ( y < 0 ) Then
      ieps = iepsz
    Else
      ieps = -iepsz
    End If
    Call ffcxr(cs3( 1),ipi12(1),yy,yy1,zz,zz1,dyyzz,.False., &
      0._DP,0._DP,0._DP,.False.,x00,2*ieps,ier)
    zz = yy*z1/y1
    zz1 = yy1*w1/y1
    dyyzz = -yy*dyz/y1
    If ( y1 > 0 ) Then
      ieps = iepsz
    Else
      ieps = -iepsz
    End If
    Call ffcxr(cs3(11),ipi12(2),yy,yy1,zz,zz1,dyyzz,.False., &
      0._DP,0._DP,0._DP,.False.,x00,2*ieps,ier)
    Do i=11,20
      cs3(i) = -cs3(i)
    End Do
    ipi12(2) = -ipi12(2)
 !  #] only cancellations in w, not in y: 
 !  #[ Hill identity:
  Else If (  ( 1 > xloss*Abs(y) .Or. Abs(xx1) > xloss ) &
      .And. ( 1 > xloss*Abs(z) .Or. Abs(z/dyz) > xloss ) &
      .And. ( 1 > xloss*Abs(y) .Or. Abs(dyz/y) > xloss ) ) Then
 !	      Do a Hill identity on the y,y-1 direction
    yy = -y*w1/dwy
    yy1 = w*y1/dwy
    zz = -z*w1/dwz
    zz1 = w*z1/dwz
    dyyzz = -w*w1*(dyz/(dwy*dwz))
    If ( y*dwz > 0 .Eqv. (y+dwz) > 0 ) Then
      ieps = 2*iepsw
    Else
      ieps = -2*iepsw
    End If
    Call ffcxr(cs3( 1),ipi12(1),yy,yy1,zz,zz1,dyyzz,.False., &
      0._DP,0._DP,0._DP,.False.,x00,ieps,ier)
    yy = w1
    yy1 = w
    zz = -w1*z/dwz
    zz1 = w*z1/dwz
    dyyzz = w*w1/dwz
    Call ffcxr(cs3( 9),ipi12(2),yy,yy1,zz,zz1,dyyzz,.False., &
      0._DP,0._DP,0._DP,.False.,x00,ieps,ier)
    Do i=9,16
      cs3(i) = -cs3(i)
    End Do
    ipi12(2) = -ipi12(2)
 !	    the extra logarithms ...
    If ( 1 < xloss*Abs(w) ) Then
      chulp = Log1MinusX(1/w)
    Else If ( w1 < 0 .Or. w < 0 ) Then
      chulp = Log(-w1/w)
    Else
      chulp = Cmplx(Real(Log(w1/w),dp),Real(-iepsw*pi,dp),dp)
    End If
    cs3(20) = -Real(Log1MinusX(dwz/dwy),dp)*chulp
 !  #] Hill identity: 
 !  #[ Taylor expansion:
  Else If ( (w<0..Or.w1<0) .And. (z<0..Or.z1<0) ) Then
 !	      Do a Taylor expansion
    If ( Abs(xx1) < xloss ) Then
      d3 = dwz/dwy
      xx1n = xx1
      d2n = d2
      d3n = d3
      d21 = 1-d2
      d21n1 = 1
      d31 = 1-d3
      d31n1 = 1
      tot = xx1*d2*d3
      Do i=2,20
        xx1n = xx1n*xx1
        d21n1 = d21n1*d21
        d31n1 = d31n1*d31
        d2n = d2n + d2*d21n1
        d3n = d3n + d3*d31n1
        term = xx1n*d2n*d3n*xn2inv(i)
        tot = tot + term
        If ( Abs(term) <= precx*Abs(tot) ) Exit
      End Do

      cs3(1) = tot
    Else If ( Abs(z/dyz) < xloss ) Then
      Call ffcxr(cs3( 1),ipi12(1),y,y1,z,z1,dyz, &
        .True.,d2yzz,zp,zp1,.False.,x00,iepsz,ier)
      Call ffcxr(cs3(11),ipi12(2),y,y1,w,w1,-dwy, &
        .True.,d2yww,wp,wp1,.False.,x00,iepsw,ier)
      Do i=11,20
        cs3(i) = -cs3(i)
      End Do
    Else
      Call WriteLfError(32)
      Return
    End If
  Else
    Call ffcxr(cs3( 1),ipi12(1),y,y1,z,z1,dyz,.False., &
      0._DP,0._DP,0._DP,.False.,x00,iepsz,ier)
    Call ffcxr(cs3(11),ipi12(2),y,y1,w,w1,-dwy,.False., &
      0._DP,0._DP,0._DP,.False.,x00,iepsw,ier)
    Do i=11,20
      cs3(i) = -cs3(i)
    End Do
    ipi12(2) = -ipi12(2)
  End If
 !  #] Taylor expansion: 
  End Subroutine ffdcxr



 Real(dp) Function FFDel2(piDpj,i1,i2,i3)
 !------------------------------------------------------------
 ! calculates in a numerical stable way
 ! fdel2(piDpj(i1,i1),piDpj(i2,i2),piDpj(i3,i3)) =
 !              = piDpj(i1,i1)*piDpj(i2,i2) - piDpj(i1,i2)^2
 !              = piDpj(i1,i1)*piDpj(i3,i3) - piDpj(i1,i3)^2
 !              = piDpj(i2,i2)*piDpj(i3,i3) - piDpj(i2,i3)^2
 !------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: piDpj(:,:)
  Integer, Intent(in) :: i1, i2, i3

  Real(dp) :: s1, s2

  If ( Abs(piDpj(i1,i2)) .Lt. Abs(piDpj(i1,i3)) .And. &
     &       Abs(piDpj(i1,i2)) .Lt. Abs(piDpj(i2,i3)) ) Then
      s1 = piDpj(i1,i1)*piDpj(i2,i2)
      s2 = piDpj(i1,i2)**2
  Elseif ( Abs(piDpj(i1,i3)) .Lt. Abs(piDpj(i2,i3)) ) Then
      s1 = piDpj(i1,i1)*piDpj(i3,i3)
      s2 = piDpj(i1,i3)**2
  Else
      s1 = piDpj(i2,i2)*piDpj(i3,i3)
      s2 = piDpj(i2,i3)**2
  Endif

  ffdel2 = s1 - s2
  
 End Function FFDel2



 Real(dp) Function ffdel3(piDpj)
 !------------------------------------------------------------------------
 ! 	Calculate del3(piDpj) = det(si.sj)	with
 !	the momenta as follows:
 !	p(1-3) = s(i)
 !	p(4-6) = p(i)
 !
 !	Input:	xpi(ns)		(real)	m^2(i),i=1,3; p^2(i-3),i=4,10
 !		piDpj(ns,ns)	(real)
 !		ns		(  Integer)
 !		ier		(  Integer)
 !
 !	Output:	ffdel3		(real)	det(si.sj)
 !------------------------------------------------------------------------
 Implicit None
  Real(dp) :: piDpj(6,6)

  Integer, Parameter :: mem=10,nperm=16
  Integer :: i,jj(6),iperm(3,nperm),imem,memarr(mem,3),memind,inow
  Real(dp) :: s(6),xmax,del3p,xmaxp
  Save iperm,memind,memarr,inow

  Data memind /0/
  Data memarr /mem*0,mem*0,mem*1/
  Data inow /1/
 !
 !	these are all permutations that give a non-zero result with the
 !	correct sign.  This list was generated with getperm3.
 !
  Data iperm/ 1,2,3,  1,2,5,  1,6,2,  1,4,3, 1,3,5,  1,4,5,  1,6,4,  1,5,6, &
    2,4,3,  2,3,6,  2,4,5,  2,6,4, 2,5,6,  3,4,5,  3,6,4,  3,5,6/

  Do i=1,mem
     If ( id == memarr(i,1) .And. idsub == memarr(i,2) ) Then
      inow = memarr(i,3)
      Exit
     End If
  End Do

 !  #[ calculations:
  imem = inow
  ffdel3 = 0
  xmax = 0
10 Continue
  jj(1) = iperm(1,inow)
  jj(3) = iperm(2,inow)
  jj(5) = iperm(3,inow)
  jj(2) = iperm(1,inow)
  jj(4) = iperm(2,inow)
  jj(6) = iperm(3,inow)
  s(1) = +piDpj(jj(1),jj(2))*piDpj(jj(3),jj(4))*piDpj(jj(5),jj(6))
  s(2) = +piDpj(jj(1),jj(4))*piDpj(jj(3),jj(6))*piDpj(jj(5),jj(2))
  s(3) = +piDpj(jj(1),jj(6))*piDpj(jj(3),jj(2))*piDpj(jj(5),jj(4))
  s(4) = -piDpj(jj(1),jj(2))*piDpj(jj(3),jj(6))*piDpj(jj(5),jj(4))
  s(5) = -piDpj(jj(1),jj(6))*piDpj(jj(3),jj(4))*piDpj(jj(5),jj(2))
  s(6) = -piDpj(jj(1),jj(4))*piDpj(jj(3),jj(2))*piDpj(jj(5),jj(6))
  del3p = 0
  xmaxp = 0
  Do i=1,6
    del3p = del3p + s(i)
    xmaxp = Max(xmaxp,Abs(s(i)))
  End Do
  If ( Abs(del3p) < xloss*xmaxp ) Then
      If ( inow == imem .Or. xmaxp < xmax ) Then
      ffdel3 = del3p
      xmax = xmaxp
      End If
    inow = inow + 1
      If ( inow > nperm ) inow = 1
      If ( inow == imem ) Goto 800
    Goto 10
  End If
  ffdel3 = del3p
  xmax = xmaxp
 !  #] calculations:
 !  #[ into memory:
800 Continue
  memind = memind + 1
  If ( memind > mem ) memind = 1
  memarr(memind,1) = id
  memarr(memind,2) = idsub
  memarr(memind,3) = inow
 !  #] into memory:
 End Function ffdel3


 Real(dp) Function ffdl2p(xpi,dpipj,piDpj,ip1,ip2,ip3,is1,is2,is3)
 !---------------------------------------------------------------------
 ! calculate in a numerically stable way
 ! delta_{ip1,is2}^{ip1,ip2}
 !---------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: ip1,ip2,ip3,is1,is2,is3
  Real(dp), Intent(in) :: xpi(:), dpipj(:,:), piDpj(:,:)

  Real(dp) :: s1,s2,s3,xmax,som

  s1 = xpi(ip1)*piDpj(ip2,is2)
  s2 = piDpj(ip1,ip2)*piDpj(ip1,is2)
  ffdl2p = s1 - s2
  If ( Abs(ffdl2p) .Ge. xloss*Abs(s1) ) Return

  som = ffdl2p
  xmax = Abs(s1)
  s1 = piDpj(ip1,ip2)*piDpj(ip3,is2)
  s2 = piDpj(ip1,ip3)*piDpj(ip2,is2)
  ffdl2p = s1 - s2
  If ( Abs(ffdl2p) .Ge. xloss*Abs(s1) ) Return

  If ( Abs(s1) .Lt. xmax ) Then
      som = ffdl2p
      xmax = Abs(s1)
  Endif
  s1 = piDpj(ip1,ip3)*piDpj(ip1,is2)
  s2 = xpi(ip1)*piDpj(ip3,is2)
  ffdl2p = s1 - s2
  If ( Abs(ffdl2p) .Ge. xloss*Abs(s1) ) Return

  If ( Abs(s1) .Lt. xmax ) Then
      som = ffdl2p
      xmax = Abs(s1)
  Endif

  s1 = xpi(ip1)*piDpj(ip2,is1)
  s2 = piDpj(ip1,is1)*piDpj(ip1,ip2)
  ffdl2p = s1 - s2
  If ( Abs(ffdl2p) .Ge. xloss*Abs(s1) ) Return

  If ( Abs(s1) .Lt. xmax ) Then
      som = ffdl2p
      xmax = Abs(s1)
  Endif

  s1 = piDpj(ip1,is2)*piDpj(ip2,is1)
  s2 = piDpj(ip1,is1)*piDpj(ip2,is2)
  ffdl2p = s1 - s2
  If ( Abs(ffdl2p) .Ge. xloss*Abs(s1) ) Return

  If ( Abs(s1) .Lt. xmax ) Then
      som = ffdl2p
      xmax = Abs(s1)
  Endif

  s1 = piDpj(ip1,ip2)*piDpj(ip3,is1)
  s2 = piDpj(ip1,ip3)*piDpj(ip2,is1)
  ffdl2p = s1 - s2
  If ( Abs(ffdl2p) .Ge. xloss*Abs(s1) ) Return
  If ( Abs(s1) .Lt. xmax ) Then
      som = ffdl2p
      xmax = Abs(s1)
  Endif

  s1 = piDpj(ip2,is2)*piDpj(ip3,is1)
  s2 = piDpj(ip2,is1)*piDpj(ip3,is2)
  ffdl2p = s1 - s2
  If ( Abs(ffdl2p) .Ge. xloss*Abs(s1) ) Return
  If ( Abs(s1) .Lt. xmax ) Then
      som = ffdl2p
      xmax = Abs(s1)
  Endif

  s1 = piDpj(ip1,ip3)*piDpj(ip1,is1)
  s2 = xpi(ip1)*piDpj(ip3,is1)
  ffdl2p = s1 - s2
  If ( Abs(ffdl2p) .Ge. xloss*Abs(s1) ) Return

  If ( Abs(s1) .Lt. xmax ) Then
      som = ffdl2p
      xmax = Abs(s1)
  Endif

  s1 = piDpj(ip1,is1)*piDpj(ip3,is2)
  s2 = piDpj(ip1,is2)*piDpj(ip3,is1)
  ffdl2p = s1 - s2
  If ( Abs(ffdl2p) .Ge. xloss*Abs(s1) ) Return

  If ( Abs(s1) .Lt. xmax ) Then
      som = ffdl2p
      xmax = Abs(s1)
  Endif

  If ( dpipj(1,1).Eq.0 ) Then
      s1 = +xpi(ip1)*dpipj(is3,is2)/2
      s2 = -piDpj(ip1,ip2)*dpipj(is2,is1)/2
      s3 = +xpi(ip1)*piDpj(ip2,ip3)/2
      ffdl2p = s1+s2+s3
      If ( Abs(ffdl2p) .Ge. xloss*Max(Abs(s1),Abs(s2)) ) Return
      If ( Max(Abs(s1),Abs(s2)) .Lt. xmax ) Then
       som = ffdl2p
       xmax = Abs(s1)
      Endif
  Endif

  ffdl2p = som ! well, no other possibility

 End Function ffdl2p


 Subroutine ffdl2s(delps1,piDpj,in,jn,jin,isji, kn,ln,lkn,islk,ns)
 !---------------------------------------------------------------------
 !	calculate in a numerically stable way
 !
 !		\delta_{si,sj}^{sk,sl}
 !
 !	with p(ji) = isji*(sj-si)
 !	     p(lk) = islk*(sl-sk)
 !---------------------------------------------------------------------
 Implicit None
  Integer :: in,jn,jin,isji,kn,ln,lkn,islk,ns
  Real(dp) :: delps1,piDpj(ns,ns)

  Integer :: ii,jj,i,j,ji,k,l,lk,ihlp
  Real(dp) :: s1,s2,som,smax

  idsub = idsub + 1
  som = 0
  smax = 0
  i = in
  j = jn
  ji = jin
  k = kn
  l = ln
  lk = lkn
  Do ii=1,3
   Do jj=1,3
    s1 = piDpj(i,k)*piDpj(j,l)
    s2 = piDpj(i,l)*piDpj(j,k)
    delps1 = s1 - s2
      If ( ii > 1 ) delps1 = isji*delps1
      If ( jj > 1 ) delps1 = islk*delps1
      If ( ii == 3 .Neqv. jj == 3 ) delps1 = -delps1
      If ( Abs(delps1) >= xloss*Abs(s1) ) Exit
 !
 !		Save the most accurate estimate so far:
      If ( ii == 1 .And. jj == 1 .Or. Abs(s1) < smax ) Then
      som = delps1
      smax = Abs(s1)
      End If
 !
 !		rotate the jj's
      If ( lk == 0 ) Cycle
    ihlp = k
    k = l
    l = lk
    lk = ihlp
   End Do
 !
 !	    and the ii's
   If ( ji == 0 ) Exit
   ihlp = i
   i = j
   j = ji
   ji = ihlp
  End Do

  delps1 = som

 End Subroutine ffdl2s

 Subroutine ffdl2t(delps,piDpj,in,jn,kn,ln,lkn,islk,iss,ns)
 !---------------------------------------------------------------------
 !	calculate in a numerically stable way
 !
 !		\delta_{si,sj}^{sk,sl}
 !
 !	with p(lk) = islk*(iss*sl - sk)	(islk,iss = +/-1)
 !	and NO relationship between s1,s2 assumed (so 1/2 the
 !	possibilities of ffdl2s).
 !---------------------------------------------------------------------
 Implicit None
  Integer :: in,jn,kn,ln,lkn,islk,iss,ns
  Real(dp) :: delps,piDpj(ns,ns)

  Real(dp) :: s1,s2,som,smax

  If ( in == jn ) Then
    delps = 0
    Return
  End If
  s1 = piDpj(kn,in)*piDpj(ln,jn)
  s2 = piDpj(ln,in)*piDpj(kn,jn)
  delps = s1 - s2
  If ( Abs(delps) >= xloss*Abs(s1) ) Return
  som = delps
  smax = Abs(s1)
  s1 = piDpj(kn,in)*piDpj(lkn,jn)
  s2 = piDpj(lkn,in)*piDpj(kn,jn)
  delps = iss*islk*(s1 - s2)
  If ( Abs(delps) >= xloss*Abs(s1) ) Return
  If ( Abs(s1) < smax ) Then
    som = delps
    smax = Abs(s1)
  End If
  s1 = piDpj(lkn,in)*piDpj(ln,jn)
  s2 = piDpj(ln,in)*piDpj(lkn,jn)
  delps = islk*(- s1 + s2)
  If ( Abs(delps) >= xloss*Abs(s1) ) Return
  If ( Abs(s1) < smax ) Then
    som = delps
    smax = Abs(s1)
  End If
 !
 !	give up
 !
  delps = som

  End Subroutine ffdl2t

 Subroutine ffdl3m(del3mi,ldel,del3,del2,xpi,dpipj,piDpj,ns,ip1n &
                  &, ip2n,ip3n,is,itime)
 !-----------------------------------------------------------------------
 !  Calculate xpi(i)*del2 - del3(piDpj)
 !
 !    /  si  mu \2    (This appears to be one of the harder
 !  = | d     |        determinants to calculate accurately.
 !    \  p1  p2 /      Note that we allow a loss of xloss^2)
 !
 !  Input: ldel          iff .true. del2 and del3 exist
 !         del3          \delta^{s(1),p1,p2}_{s(1),p1,p2}
 !         del2          \delta^{p1,p2}_{p1,p2}
 !         xpi(ns)       standard
 !         dpipj(ns,ns)  standard
 !         piDpj(ns,ns)  standard
 !         ipi           pi = xpi(abs(ipi)) [p3=-p1 +/-p2]
 !         is            si = xpi(is,is+1,..,is+itime-1)
 !         itime         number of functions to calculate
 !   Output:  del3mi(3)  (\delta^{s_i \mu}_{p_1 p_2})^2
 !-----------------------------------------------------------------------
 Implicit None
  Integer :: ns,ip1n,ip2n,ip3n,is,itime
  Logical, Intent(in) :: ldel
  Real(dp) :: del3mi(itime),del3,del2,xpi(ns),dpipj(ns,ns),piDpj(ns,ns)

  Real(dp) :: s(7),som,smax,xsom,xmax
  Integer :: i,j,k,ip1,ip2,ip3,ipn,is1,is2,isi,is3,ihlp,iqn, &
       &    jsgn1,jsgn2,jsgn3,jsgnn,nm
  Integer, Save :: iadj(10,10,3:4),init
  Logical :: lmax,ltwist

  Data iadj /200*0/
  Data init /0/

  If ( init .Eq. 0 ) Then
   init = 1

 !      Fill the array with adjacent values: if
 !    x = iadj(i,j)
 !    k = abs(mod(k,100))
 !    jsgnk = sign(x)
 !    jsgnj = 1-2*theta(x-100)  (ie -1 iff |x|>100)
 !      then
 !    pi(k) = jsgnk*( p(i) - jsgnj*pi(j) )

   Do nm=3,4
    Do  i=1,nm
     is1 = i
     is2 = i+1
     If ( is2 .Gt. nm ) is2 = 1
     is3 = i-1
     If ( is3 .Eq. 0 ) is3 = nm
     ip1 = is1 + nm
     iadj(is1,is2,nm) = -ip1
     iadj(is2,is1,nm) = ip1
     iadj(ip1,is2,nm) = -is1
     iadj(is2,ip1,nm) = is1
     iadj(is1,ip1,nm) = 100+is2
     iadj(ip1,is1,nm) = 100+is2
     If ( nm .Eq. 3 ) Then
      iadj(ip1,is2+3,3) = -100-is3-3
      iadj(is2+3,ip1,3) = -100-is3-3
     Endif
    End Do
   End Do

   iadj(3,1,4) = -9
   iadj(1,3,4) = 9
   iadj(9,1,4) = -3
   iadj(1,9,4) = 3
   iadj(3,9,4) = 100+1
   iadj(9,3,4) = 100+1

   iadj(2,4,4) = -10
   iadj(4,2,4) = 10
   iadj(10,4,4) = -2
   iadj(4,10,4) = 2
   iadj(2,10,4) = 100+4
   iadj(10,2,4) = 100+4
  Endif ! init

  If ( ns .Eq. 6 ) Then
   nm = 3
  Else
   nm = 4
  Endif

outer:  Do i=1,itime
          isi = i+is-1
          lmax = .False.

  !      get xpi(isi)*del2 - del3 ... if del3 and del2 are defined

          If ( ldel ) Then
           s(1) = xpi(isi)*del2
           som = s(1) - del3
           smax = Abs(s(1))
           If ( Abs(som) .Ge. xloss**2*smax ) Then
            del3mi(i) = som
            Cycle outer
           End If
           xsom = som
           xmax = smax
           lmax = .True.
          Endif
          ip1 = ip1n
          ip2 = ip2n
          ip3 = ip3n
    inner: Do j=1,3 !    otherwise use the simple threeterm formula
             s(1) = xpi(ip2)*piDpj(ip1,isi)**2
             s(2) = xpi(ip1)*piDpj(ip2,isi)*piDpj(ip2,isi)
             s(3) = -2*piDpj(ip2,isi)*piDpj(ip2,ip1)*piDpj(ip1,isi)
             som = s(1) + s(2) + s(3)
             smax = Max(Abs(s(1)),Abs(s(2)),Abs(s(3)))
             If ( Abs(som) .Ge. xloss**2*smax ) Then
              del3mi(i) = som
              Cycle outer
             End If
             If ( .Not. lmax .Or. smax .Lt. xmax ) Then
              xsom = som
              xmax = smax
              lmax = .True.
             Endif

 !    if there are cancellations between two of the terms:
 !    we try mixing with isi.
 !
 !    First map cancellation to s(2)+s(3) (do not mess up rotations...)

             If ( Abs(s(1)+s(3)) .Lt. Abs(s(3))/2 ) Then
              ihlp = ip1
              ip1 = ip2
              ip2 = ihlp
              som = s(1)
              s(1) = s(2)
              s(2) = som
              ltwist = .True.
             Else
              ltwist = .False.
             Endif
             If ( Abs(s(2)+s(3)) .Lt. Abs(s(3))/2 ) Then
               !  switch to the vector pn so that si = jsgn1*p1 + jsgnn*pn

              k = iadj(isi,ip1,nm)
              If ( k .Ne. 0 ) Then
               ipn = Abs(k)
               jsgnn = isign(1,k)
               If ( ipn .Gt. 100 ) Then
                ipn = ipn - 100
                jsgn1 = -1
               Else
                jsgn1 = +1
               Endif
               If (Abs(dpipj(ipn,isi)).Lt.xloss*Abs(piDpj(ip1,isi)).And.  &
                  & Abs(piDpj(ipn,ip2)).Lt.xloss*Abs(piDpj(ip2,isi)) ) Then
                !    same:  s(1) = xpi(ip2)*piDpj(ip1,isi)**2
                s(2) = jsgnn*piDpj(isi,ip2)*piDpj(ipn,ip2)*xpi(ip1)
                s(3) = jsgn1*piDpj(isi,ip2)*piDpj(ip1,ip2)*dpipj(ipn,isi)
                som = s(1) + s(2) + s(3)
                smax = Max(Abs(s(1)),Abs(s(2)),Abs(s(3)))
                If ( Abs(som) .Ge. xloss**2*smax ) Then
                 del3mi(i) = som
                 Cycle outer
                End If
                If ( smax .Lt. xmax ) Then
                 xsom = som
                 xmax = smax
                Endif

              !  there may be a cancellation between s(1) and
              !  s(2) left.  Introduce a vector q such that
              !  pn = jsgnq*q + jsgn2*p2.  We also need the sign
              !  jsgn3 in p3 = -p1 - jsgn3*p2

                k = iadj(ipn,ip2,nm)
                If ( k .Ne. 0 ) Then
                 iqn = Abs(k)
                 If ( iqn .Gt. 100 ) Then
                  iqn = iqn - 100
                  jsgn2 = -1
                 Else
                  jsgn2 = +1
                 Endif
                 k = iadj(ip1,ip2,nm)
                 If ( k .Eq. 0 .Or. k .Lt. 100 ) Then !we have p1,p2,p3 all p's
                  jsgn3 = +1
                 Elseif ( k .Lt. 0 ) Then
                 ! ip1,ip2 are 2*s,1*p such that p2-p1=ip3
                  jsgn3 = -1
                 Else
                  jsgn3 = 0
                 Endif
                 ! we need one condition on the signs for this to work
                 If ( ip3.Ne.0 .And. jsgn1*jsgn2.Eq.jsgnn*jsgn3  &
                   & .And. Abs(s(3)).Lt.xloss*smax ) Then
                  s(1) = piDpj(ip1,isi)**2*dpipj(iqn,ipn)
                  s(2) = -jsgn2*jsgn1*piDpj(ipn,ip2) &
                      &   *piDpj(ip1,isi)*dpipj(ipn,isi)
                !  s(3) stays the same
                  s(4) = -jsgn2*jsgn1*piDpj(ipn,ip2)*xpi(ip1)*piDpj(isi,ip3)
                  som = s(1) + s(2) + s(3) + s(4)
                  smax =Max(Abs(s(1)),Abs(s(2)),Abs(s(3)),Abs(s(4)))
                  If ( Abs(som).Ge.xloss**2*smax ) Then
                   del3mi(i) = som
                   Cycle
                  End If
                  If ( smax .Lt. xmax ) Then
                     xsom = som
                     xmax = smax
                  Endif
                 Endif
                Endif
               Endif
              Endif
             k = iadj(isi,ip2,nm)
             If ( k .Ne. 0 ) Then
              ipn = Abs(k)
              jsgnn = isign(1,k)
              If ( ipn .Gt. 100 ) Then
               jsgn1 = -1
               ipn = ipn - 100
              Else
               jsgn1 = +1
              Endif
              If (Abs(dpipj(ipn,isi)).Lt.xloss*Abs(piDpj(ip2,isi)).And.  &
                 &  Abs(piDpj(ipn,ip1)).Lt.xloss*Abs(piDpj(ip1,isi)) ) Then
               s(1) = jsgnn*piDpj(isi,ip1)*piDpj(ipn,ip1)*xpi(ip2)
               s(2) = xpi(ip1)*piDpj(ip2,isi)**2
               s(3) = jsgn1*piDpj(isi,ip1)*piDpj(ip2,ip1)*dpipj(ipn,isi)
               som = s(1) + s(2) + s(3)
               smax = Max(Abs(s(1)),Abs(s(2)),Abs(s(3)))
               Print *,'    (isi+ip2) with isi,ip1,ip2,ipn: ',isi,ip1,ip2,ipn
               If ( Abs(som) .Ge. xloss**2*smax ) Then
                del3mi(i) = som
                Cycle outer
               End If
               If ( smax .Lt. xmax ) Then
                xsom = som
                xmax = smax
               Endif
              Endif
             Endif
            Endif

         !    rotate the ipi
            If ( ip3 .Eq. 0 ) Exit inner
            If ( j .Ne. 3 ) Then
             If ( .Not. ltwist ) Then
              ihlp = ip1
              ip1 = ip2
              ip2 = ip3
              ip3 = ihlp
             Else
              ihlp = ip2
              ip2 = ip3
              ip3 = ihlp
             Endif
            Endif
           End Do inner ! j
     ! These values are the best found:
          som = xsom
          smax = xmax
          del3mi(i) = som
      End Do outer ! i
 End Subroutine ffdl3m


 Subroutine ffdl3s(dl3s,piDpj,ii,ns)
 !------------------------------------------------------------------------
 !	Calculate dl3s(piDpj) = det(si.sj)	with
 !	the momenta indicated by the indices ii(1-6,1), ii(1-6,2)
 !	as follows:
 !	p(|ii(1,)|-|ii(3,)|) = s(i)
 !	p(|ii(4,)|-|ii(6,)|) = p(i) = sgn(ii())*(s(i+1) - s(i))
 !
 !	At this moment (26-apr-1990) only the diagonal is tried
 !
 !	Input:	xpi(ns)		(real)	m^2(i),i=1,3; p^2(i-3),i=4,10
 !		piDpj(ns,ns)	(real)
 !		ii(6,2)		(  Integer)	see above
 !		ns		(  Integer)
 !		ier		(  Integer)
 !
 !	Output:	dl3s		(real)	det(si.sj)
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ii(6,2),ns
  Real(dp) :: dl3s,piDpj(ns,ns)

  Integer,Parameter :: mem=10,nperm=16
  Integer :: i,jj(6),jsgn,iperm(3,nperm),imem,memarr(mem,3), memind,inow
  Real(dp) :: s(6),xmax,dl3sp,xmaxp
  Save iperm,memind,memarr,inow

  Data memind /0/
  Data memarr /mem*0,mem*0,mem*1/
  Data inow /1/
 !
 !	these are all permutations that give a non-zero result with the
 !	correct sign.  This list was generated with getperm3.
 !
  Data iperm/ 1,2,3,  1,2,5,  1,6,2,  1,4,3, 1,3,5,  1,4,5,  1,6,4,  1,5,6, &
    2,4,3,  2,3,6,  2,4,5,  2,6,4, 2,5,6,  3,4,5,  3,6,4,  3,5,6/

  Do i=1,mem
     If ( id == memarr(i,1) .And. idsub == memarr(i,2) ) Then
      inow = memarr(i,3)
      Exit
     End If
  End Do

 !  #] starting point in memory?:
 !  #[ calculations:
  imem = inow
  dl3s = 0
  xmax = 0
10 Continue
  jj(1) = Abs(ii(iperm(1,inow),1))
  jj(3) = Abs(ii(iperm(2,inow),1))
  jj(5) = Abs(ii(iperm(3,inow),1))
  jj(2) = Abs(ii(iperm(1,inow),2))
  jj(4) = Abs(ii(iperm(2,inow),2))
  jj(6) = Abs(ii(iperm(3,inow),2))
  jsgn =  Sign(1,ii(iperm(1,inow),1)) *Sign(1,ii(iperm(2,inow),1)) &
    *Sign(1,ii(iperm(3,inow),1)) *Sign(1,ii(iperm(1,inow),2)) &
    *Sign(1,ii(iperm(2,inow),2)) *Sign(1,ii(iperm(3,inow),2))
  s(1) = +piDpj(jj(1),jj(2))*piDpj(jj(3),jj(4))*piDpj(jj(5),jj(6))
  s(2) = +piDpj(jj(1),jj(4))*piDpj(jj(3),jj(6))*piDpj(jj(5),jj(2))
  s(3) = +piDpj(jj(1),jj(6))*piDpj(jj(3),jj(2))*piDpj(jj(5),jj(4))
  s(4) = -piDpj(jj(1),jj(2))*piDpj(jj(3),jj(6))*piDpj(jj(5),jj(4))
  s(5) = -piDpj(jj(1),jj(6))*piDpj(jj(3),jj(4))*piDpj(jj(5),jj(2))
  s(6) = -piDpj(jj(1),jj(4))*piDpj(jj(3),jj(2))*piDpj(jj(5),jj(6))
  dl3sp = 0
  xmaxp = 0
  Do i=1,6
    dl3sp = dl3sp + s(i)
    xmaxp = Max(xmaxp,Abs(s(i)))
  End Do
  If ( Abs(dl3sp) < xloss*xmaxp ) Then
      If ( inow == imem .Or. xmaxp < xmax ) Then
      dl3s = jsgn*dl3sp
      xmax = xmaxp
      End If
    inow = inow + 1
      If ( inow > nperm ) inow = 1
      If ( inow == imem ) Goto 800
    Goto 10
  End If
  dl3s = jsgn*dl3sp
  xmax = xmaxp
 !  #] calculations:
 !  #[ into memory:
800 Continue
  memind = memind + 1
  If ( memind > mem ) memind = 1
  memarr(memind,1) = id
  memarr(memind,2) = idsub
  memarr(memind,3) = inow
 !  #] into memory:

 End Subroutine ffdl3s


 Subroutine Ffdot3(piDpj,xpi,dpipj)
 Implicit None
  Real(dp), Intent(in) :: xpi(6), dpipj(6,6)
  Real(dp), Intent(out) :: piDpj(6,6)

  Integer :: is1, is2, is3, ip1, ip2, ip3, i, j
  Integer, Parameter :: inew(6,6) = Reshape( &
     & Source =  (/1,2,3,4,5,6,  &
     &             2,3,1,5,6,4,  &
     &             3,1,2,6,4,5,  &
     &             1,3,2,6,5,4,  &
     &             3,2,1,5,4,6,  &
     &             2,1,3,4,6,5/), Shape = (/6, 6/) )

   If ( idot.Ge.3 ) Then
    Do i=1,6
     Do j=1,6
      piDpj(inew(j,irota3),inew(i,irota3)) = fpij3(j,i)
     End Do
    End Do

    If ( irota3 .Gt. 3 ) Then ! the sign of the s's has been changed!
     Do i=1,3
      Do j=4,6
       piDpj(j,i) = -piDpj(j,i)
       piDpj(i,j) = -piDpj(i,j)
      End Do
     End Do
    End If
    Return
   Endif

   Do is1=1,3
    is2 = is1 + 1
    If ( is2 .Eq. 4 ) is2 = 1
    is3 = is2 + 1
    If ( is3 .Eq. 4 ) is3 = 1
    ip1 = is1 + 3
    ip2 = is2 + 3
    ip3 = is3 + 3

    !  pi.pj, si.sj

    piDpj(is1,is1) = xpi(is1)
    piDpj(ip1,ip1) = xpi(ip1)

    ! si.s(i+1)

    If ( xpi(is2) .Le. xpi(is1) ) Then
      piDpj(is1,is2) = (dpipj(is1,ip1) + xpi(is2))/2
    Else
      piDpj(is1,is2) = (dpipj(is2,ip1) + xpi(is1))/2
    Endif
    piDpj(is2,is1) = piDpj(is1,is2)

    !  pi.si

    If ( Abs(xpi(ip1)) .Le. xpi(is1) ) Then
      piDpj(ip1,is1) = (dpipj(is2,is1) - xpi(ip1))/2
    Else
      piDpj(ip1,is1) = (dpipj(is2,ip1) - xpi(is1))/2
    Endif
    piDpj(is1,ip1) = piDpj(ip1,is1)

    !  pi.s(i+1)

    If ( Abs(xpi(ip1)) .Le. xpi(is2) ) Then
      piDpj(ip1,is2) = (dpipj(is2,is1) + xpi(ip1))/2
    Else
      piDpj(ip1,is2) = (dpipj(ip1,is1) + xpi(is2))/2
    Endif
    piDpj(is2,ip1) = piDpj(ip1,is2)

    ! pi.s(i+2)

    If ( Min(Abs(dpipj(is2,is1)),Abs(dpipj(ip3,ip2))) .Le.    &
       &     Min(Abs(dpipj(ip3,is1)),Abs(dpipj(is2,ip2))) ) Then
     piDpj(ip1,is3) = (dpipj(ip3,ip2) + dpipj(is2,is1))/2
    Else
     piDpj(ip1,is3) = (dpipj(ip3,is1) + dpipj(is2,ip2))/2
    Endif
    piDpj(is3,ip1) = piDpj(ip1,is3)

    !  pi.p(i+1)

    If ( idot.Le.0 ) Then
     If ( Abs(xpi(ip2)) .Le. Abs(xpi(ip1)) ) Then
       piDpj(ip1,ip2) = (dpipj(ip3,ip1) - xpi(ip2))/2
     Else
       piDpj(ip1,ip2) = (dpipj(ip3,ip2) - xpi(ip1))/2
     Endif
     piDpj(ip2,ip1) = piDpj(ip1,ip2)
    Else
     piDpj(inew(ip2,irota3),inew(ip1,irota3)) = fpij3(ip1,ip2)
     piDpj(inew(ip1,irota3),inew(ip2,irota3)) =  &
          &      piDpj(inew(ip2,irota3),inew(ip1,irota3))
    Endif
   End Do
  
 End Subroutine Ffdot3

 Subroutine ffdwz(dwz,z,i1,j1,l,alpha,alph1,xpi,dpipj,piDpj, sdel2i,ns,ier)
 !----------------------------------------------------------------------
 !	Recalculate dwz(i1,j1) = w(i1) - z(j1)
 !----------------------------------------------------------------------
 Implicit None
  Integer :: i1,j1,l,ns,ier
  Real(dp) :: dwz(2,2),z(4)
  Real(dp) :: alpha,alph1,xpi(ns),dpipj(ns,ns),piDpj(ns,ns), sdel2i(3)

  Real(dp) :: s(8),sum,fac,xmax
  Integer :: i

 !  #[ calculations:
  If ( l == 1 ) Then
    ier = ier + 100
  Else If ( l == 3 ) Then
    If ( (i1==2 .And. j1==1) .Or. (i1==1 .And. j1==2) ) Then
      fac = 1._dp/(sdel2i(2) + sdel2i(3))
      s(1) = dpipj(6,5)*z(j1)
      s(2) = -alph1*xpi(5)*z(j1+2)
      If ( Max(Abs(dpipj(2,1)),Abs(dpipj(5,6))) <&
            &Max(Abs(dpipj(2,6)),Abs(dpipj(5,1))) ) Then
        s(3) = 0.5_dp*dpipj(2,1)
        s(4) = 0.5_dp*dpipj(5,6)
      Else
        s(3) = 0.5_dp*dpipj(2,6)
        s(4) = 0.5_dp*dpipj(5,1)
      End If
      s(5) = piDpj(4,3)*piDpj(5,3)*fac
      s(6) = -piDpj(4,3)*piDpj(6,3)*fac
      s(7) = xpi(3)*dpipj(5,6)*fac
      If ( i1 == 1 ) Then
        sum = s(1)+s(2)+s(3)+s(4) - (s(5)+s(6)+s(7))
      Else
        sum = s(1)+s(2)+s(3)+s(4) + s(5)+s(6)+s(7)
      End If
      xmax = Abs(s(1))
      Do i=2,7
        xmax = Max(xmax,Abs(s(i)))
      End Do
      If ( Abs(sum) < xloss*xmax ) Then
 !		    this result is not used   If it is not accurate (see
 !		    ffxc0p)
        ier = ier + 1
        xmax = xmax/Abs(alpha*xpi(5))
        dwz(i1,j1) = sum/(alpha*xpi(5))
      Else
        dwz(i1,j1) = sum/(alpha*xpi(5))
      End If
    Else
      ier = ier + 100
    End If
  End If
 !  #] calculations:
  End Subroutine ffdwz

 Subroutine ffgeta(ni,cz,cdyz,cp,cpDs,ieps,isoort,ier)
 !--------------------------------------------------------------------
 !	Get the eta terms which arise from splitting up
 !	log(p2(x-z-)(x-z+)) - log(p2(y-z-)(y-z+))
 !
 !	Input:	cz	complex(4)	the roots z-,z+,1-z-,1-z+
 !		cdyz	complex(2,2)	y-z
 !		cd2yzz	complex(2)	2y-(z-)-(z+)
 !		cp	complex		p^2
 !		cpDs	complex		p.s
 !		ieps	  Integer(2)	the assumed im part if Im(z)=0
 !		isoort	  Integer(2)	which type of Ri
 !
 !	Output:	ni	  Integer(4)	eta()/(2*pi*i)
 !--------------------------------------------------------------------
 Implicit None
  Integer :: ni(4),ieps(2),isoort(2),ier
  Complex(dp) :: cp,cpDs,cz(4),cdyz(2,2)

  Integer :: i
  Complex(dp) :: cmip
 !
 !  #[ complex  masses or imaginary roots:
 !
 !	only complex because of complex roots in y or z
 !	[checked and in agreement with ieps definition 23-sep-1991]
 !
 !	isoort = +1:        y is real, z is real
 !	isoort = -1-n*10:   y is complex, possibly z as well
 !	isoort = -3-n*10:   y,z complex, (y-z-)*(y-z+) real
 !	isoort = 0:         y is complex, one z root only
 !	isoort = -10-n*10:  y is real, z is complex
 !	isoort = -5,6-n*10: y,z real
 !
  If ( isoort(1) > 0 ) Then
 !
 !	    really a real case
 !
    ni(1) = 0
    ni(2) = 0
    ni(3) = 0
    ni(4) = 0
  Else If ( Mod(isoort(1),10) /= 0 .And. isoort(2) /= 0 ) Then
    cmip = Cmplx(0._DP,-Real(cp,dp),dp)
 !
 !	    ni(1) = eta(p2,(x-z-)(x-z+)) = 0 by definition (see ni(3))
 !	    ni(2) = eta(x-z-,x-z+)
 !
    ni(1) = 0
    If ( ieps(1) > 0 .Neqv. ieps(2) > 0 ) Then
      ni(2) = 0
    Else
      ni(2) = nffet1(-cz(1),-cz(2),cmip,ier)
      If ( cz(3)/=0 .And. cz(4)/=0 ) Then
        i = nffet1(cz(3),cz(4),cmip,ier)
        If ( i /= ni(2) )   Call WriteLFerror(23)
      End If
    End If
 !
 !	    ni(3) compensates for whatever convention we chose in ni(1)
 !	    ni(4) = -eta(y-z-,y-z+)
 !
    If ( Mod(isoort(1),10)==-3 ) Then
 !		follow the i*epsilon prescription as (y-z-)(y-z+) real
      ni(3) = 0
      ni(4) = -nffet1(cdyz(2,1),cdyz(2,2),cmip,ier)
    Else
      If ( Real(cp,dp) < 0 .And. Aimag(cdyz(2,1)*cdyz(2,2)) < 0 ) Then
        ni(3) = -1
      Else
        ni(3) = 0
      End If
      ni(4) = -nffeta(cdyz(2,1),cdyz(2,2),ier)
    End If
  Else If ( (Mod(isoort(1),10)==-1 .Or. Mod(isoort(1),10)==-3) &
      .And. isoort(2) == 0 ) Then
    ni(1) = 0
    If ( Aimag(cz(1)) /= 0 ) Then
      ni(2) = nffet1(-cpDs,-cz(1),Cmplx(0._dp,-1._dp,dp),ier)
    Else
      ni(2) = nffet1(-cpDs,Cmplx(0._dp,-1._dp,dp),Cmplx(0._dp,-1._dp,dp),ier)
    End If
    ni(3) = 0
    ni(4) = -nffeta(-cpDs,cdyz(2,1),ier)
  Else
    ni(1) = 0
    ni(2) = 0
    ni(3) = 0
    ni(4) = 0
  End If
 !  #] complex  masses or imaginary roots:
  End Subroutine ffgeta

 Subroutine ffieps(ieps,cz,cp,cpDs,isoort)
 !--------------------------------------------------------------------
 !	Get the ieps prescription in such a way that it is compatible
 !	with the imaginary part of cz   If non-zero, compatible with the
 !	real case   If zero.
 !
 !	Input:	cz	complex(4)	the roots z-,z+,1-z-,1-z+
 !		cp	complex		p^2
 !		cpDs	complex		p.s
 !		isoort	  Integer(2)	which type of Ri
 !
 !	Output:	ieps	  Integer(2)	z -> z-ieps*i*epsilon
 !					will give correct im part
 !--------------------------------------------------------------------
 Implicit None
  Integer :: ieps(2),isoort(2)
  Complex(dp) :: cp,cpDs,cz(4)
 !
 !  #[ work:
  If ( Aimag(cp) /= 0 ) Then
 !	      Do not calculate ANY eta terms, we'll do them ourselves.
    ieps(1) = 99
    ieps(2) = 99
  Else If ( isoort(2) /= 0 ) Then
    If ( Aimag(cz(1)) < 0 ) Then
      ieps(1) = +1
      If ( Aimag(cz(2)) < 0 ) Then
        ieps(2) = +1
      Else
        ieps(2) = -1
      End If
    Else If ( Aimag(cz(1)) > 0 ) Then
      ieps(1) = -1
      If ( Aimag(cz(2)) <= 0 ) Then
        ieps(2) = +1
      Else
        ieps(2) = -1
      End If
    Else
      If ( Aimag(cz(2)) < 0 ) Then
        ieps(1) = -1
        ieps(2) = +1
      Else If ( Aimag(cz(2)) > 0 ) Then
        ieps(1) = +1
        ieps(2) = -1
      Else
        If ( (Real(cz(2),dp)>Real(cz(1),dp) &
           &   .Or. (Real(cz(1),dp)==Real(cz(2),dp) &
          .And. Real(cz(4),dp)<Real(cz(3),dp)) ) .Eqv. Real(cp,dp)>0 ) Then
          ieps(1) = +1
          ieps(2) = -1
        Else
          ieps(1) = -1
          ieps(2) = +1
        End If
      End If
    End If
  Else
    If ( Aimag(cz(1)) < 0 ) Then
      ieps(1) = +1
    Else If ( Aimag(cz(1)) > 0 ) Then
      ieps(1) = -1
    Else If ( Real(cpDs,dp) > 0 ) Then
      ieps(1) = +1
    Else
      ieps(1) = -1
    End If
    ieps(2) = -9999
  End If
 !  #] work:
  End Subroutine ffieps

 Subroutine ffrot3(irota,xqi,dqiqj,qiDqj,xpi,dpipj,piDpj,iflag,npoin)
 !-----------------------------------------------------------------------
 ! rotates the arrays xpi, dpipj into xqi,dqiqj so that
 ! xpi(6),xpi(4) suffer the strongest outside cancellations and
 ! xpi(6) > xpi(4) if iflag = 1, so that xpi(5) largest and xpi(5)
 ! and xpi(6) suffer cancellations if iflag = 2.
 ! if iflag = 3 make xqi(3)=0.
 ! If npoin=4, rotate piDpj into qiDqj as well.
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: iflag, npoin
  Integer, Intent(inout) :: irota
  Real(dp), Intent(inout) :: xpi(6), dpipj(6,6), piDpj(6,6), xqi(6) &
     & , dqiqj(6,6), qiDqj(6,6)

  Real(dp) :: a1, a2, a3, xpimax
  Complex(dp) :: chulp(3,3)
  Integer :: i, j
  Integer, Parameter :: inew(6,6) = Reshape( &
     & Source =  (/1,2,3,4,5,6,  &
     &             2,3,1,5,6,4,  &
     &             3,1,2,6,4,5,  &
     &             1,3,2,6,5,4,  &
     &             3,2,1,5,4,6,  &
     &             2,1,3,4,6,5/), Shape = (/6, 6/) )

  If ( iflag .Eq. 1 ) Then
      a1 = Abs(dpipj(6,4))/Max(Abs(xpi(6)+xpi(4)),xalogm)
      a2 = Abs(dpipj(5,4))/Max(Abs(xpi(5)+xpi(4)),xalogm)
      a3 = Abs(dpipj(5,6))/Max(Abs(xpi(6)+xpi(5)),xalogm)
   If ( a1 .Le. a2 .And. a1 .Le. a3 ) Then
    irota = 1
    If ( Abs(xpi(6)) .Lt. Abs(xpi(4)) )  irota = 4
   Elseif ( a2 .Le. a3 ) Then
    irota = 3
    If ( Abs(xpi(4)) .Lt. Abs(xpi(5)) ) irota = 6
   Else
    irota = 2
    If ( Abs(xpi(5)) .Lt. Abs(xpi(6)) ) irota = 5
   Endif

  Elseif ( iflag .Eq. 2 ) Then
   xpimax = Max(xpi(4),xpi(5),xpi(6))
   If ( xpimax .Eq. 0 ) Then
    If ( xpi(5) .Ne. 0 ) Then
     irota = 1
    Elseif ( xpi(4) .Ne. 0 ) Then
     irota = 2
    Elseif ( xpi(6) .Ne. 0 ) Then
     irota = 3
    Else
     Call WriteLFerror(13)
     irota = 1
    Endif
   Elseif ( xpi(5) .Eq. xpimax ) Then
    If ( xpi(4) .Le. xpi(6) ) Then
     irota = 1
    Else
     irota = 4
    Endif
   Elseif ( xpi(4) .Eq. xpimax ) Then
    If ( xpi(5) .Ge. xpi(6) ) Then
     irota = 2
    Else
     irota = 5
    Endif
   Else
    If ( xpi(4) .Ge. xpi(6) ) Then
     irota = 3
    Else
     irota = 6
    Endif
   Endif
  Elseif ( iflag .Eq. 3 ) Then
   If ( dpipj(2,4).Eq.0 .And. dpipj(3,6).Eq.0 .And. xpi(1).Eq.0 ) Then
    irota = 3
   Elseif ( dpipj(1,6).Eq.0 .And. dpipj(2,5).Eq.0 .And. xpi(3).Eq.0 ) Then
    irota = 1
   Elseif ( dpipj(3,5).Eq.0 .And. dpipj(1,4).Eq.0 .And. xpi(2).Eq.0 ) Then
    irota = 2
   Else
    Call WriteLFerror(14)
    irota = 1
   Endif
  Else
   Call WriteLFerror(14)
   irota = 1
  Endif

  Do i=1,6
   xqi(inew(i,irota)) = xpi(i)
   Do j=1,6
    dqiqj(inew(i,irota),inew(j,irota)) = dpipj(i,j)
   End Do
  End Do

 ! when called in a 4pointfunction we already have the dotproducts

  If ( npoin .Eq. 4 ) Then
   Do j=1,6
    Do i=1,6
     qiDqj(inew(i,irota),inew(j,irota)) = piDpj(i,j)
    End Do
   End Do
  Endif

 !DEBUG	  If ( iflag .eq. 3 .and. lsmug ) then
  If ( lsmug ) Then
 !	    
 !	      Do not forget to rotate the smuggled d  Ifferences
 !	    
    Do j=1,3
      Do i=1,3
        chulp(i,j) = cmipj(i,j)
      End Do
    End Do
    Do j=1,3
      Do i=1,3
        cmipj(inew(i,irota),inew(j+3,irota)-3) = chulp(i,j)
      End Do
    End Do
  End If
 !  #] rotate:

 End Subroutine ffrot3

 Subroutine ffrt3p(clogip,ilogip,irota,clogi,ilogi,idir)
 !------------------------------------------------------------------------
 !	rotates the arrays clogi,ilogi also over irota (idir=+1) or
 !	back (-1)
 !
 !	Input:	irota	  (Integer)	index in rotation array
 !		clogi(3)  (complex)	only if idir=-1
 !		ilogi(3)  (Integer)	indicates which clogi are needed
 !					(idir=+1), i*pi terms (idir=-1)
 !		idir	  (Integer)	direction: forward (+1) or
 !					backward (-1)
 !	Output:	clogip(3)  (Integer)	clogi rotated
 !		ilogip(3)  (Integer)	ilogi rotated
 !------------------------------------------------------------------------
 Implicit None
  Integer :: irota,idir,ilogi(3),ilogip(3)
  Complex(dp) :: clogi(3),clogip(3)

  Integer :: i,inew(6,6)
  Save inew

 !	data
 !
  Data inew /1,2,3,4,5,6, 2,3,1,5,6,4, 3,1,2,6,4,5, 1,3,2,6,5,4, 3,2,1,5,4,6, &
    2,1,3,4,6,5/
 !  #[ rotate:
 !
 !	the clogi, ilogi are numbered according to the p_i
 !
  If ( idir == +1 ) Then
    Do i=1,3
      ilogip(inew(i+3,irota)-3) = ilogi(i)
      clogip(inew(i+3,irota)-3) = clogi(i)
    End Do
  Else
    Do i=1,3
      ilogip(i) = ilogi(inew(i+3,irota)-3)
      clogip(i) = clogi(inew(i+3,irota)-3)
    End Do
  End If
 !
 !  #] rotate:
  End Subroutine ffrt3p

 Subroutine ffxc0i(cc0,xpi,dpipj,ier)
 !------------------------------------------------------------------------
 !	Calculates the infrared finite part of a infrared divergent
 !	threepoint function with the factor ipi^2.  The cutoff
 !	parameter is assumed to be in /ffregul/.
 !
 !	Input:	xpi(6)		(real)	pi.pi (B&D)
 !		dpipj(6,6)	(real)	xpi(i)-xpi(j)
 !		lambda2		(real)	cutoff (either photon mass**2
 !						or radiation limit).
 !	Output: cc0	(complex)	C0, the threepoint function.
 !		ier	(  Integer)	usual error flag
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ier
  Complex(dp) :: cc0
  Real(dp) :: xpi(6),dpipj(6,6)

  Integer :: init,ipi12,i,ilogi(3),irota,n
  Integer :: j,inew(6,6)
  Complex(dp) :: cs(15),csum,clogi(3)
  Real(dp) :: xqi(6),dqiqj(6,6),qiDqj(6,6),sdel2,xmax, dum66(6,6),del2
 !
 !	data
 !
  Data init /0/
  Data inew /1,2,3,4,5,6, 2,3,1,5,6,4, 3,1,2,6,4,5, 1,3,2,6,5,4, 3,2,1,5,4,6, &
    2,1,3,4,6,5/
  Data ilogi /3*0/
 !
 !
 !	initialisations
 !
  Do i=1,15
    cs(i) = 0
  End Do
  ipi12 = 0
 !  #[ check input:
  If ( init == 0 .And. .Not.lsmug ) Then
    init = 1
    Print *,'ffxc0i: infra-red divergent threepoint function, ', &
            'working with a cutoff ',lambda2
  End If
  If ( .Not.lsmug .And. lambda2 == 0 ) Then
      Call WriteLfError(37)
    Return
  End If
 !  #] check input:
 !  #[ groundwork:
 !
 !	rotate to xpi(3)=0, xpi(1)=xpi(6), xpi(2)=xpi(5)
 !
  Call ffrot3(irota,xqi,dqiqj,qiDqj,xpi,dpipj,dum66,3,3)
 !
 !	get some   Dotproducts
 !
  If ( ldot ) Then
      Call ffDot3(qiDqj,xqi,dqiqj)
      Do i=1,6
        Do j=1,6
          fpij3(j,i) = qiDqj(inew(i,irota),inew(j,irota))
        End Do
      End Do
  Else
    If ( Abs(xqi(4)) < xqi(1) ) Then
      qiDqj(4,1) = dqiqj(2,1) - xqi(4)
      xmax = Abs(xqi(4))
    Else
      qiDqj(4,1) = dqiqj(2,4) - xqi(1)
      xmax = xqi(1)
    End If
    qiDqj(4,1) = qiDqj(4,1)/2
    qiDqj(1,4) = qiDqj(4,1)
    If ( Abs(xqi(4)) < xqi(2) ) Then
      qiDqj(4,2) = dqiqj(2,1) + xqi(4)
      xmax = Abs(xqi(4))
    Else
      qiDqj(4,2) = xqi(2) - dqiqj(1,4)
      xmax = xqi(2)
    End If
    qiDqj(4,2) = qiDqj(4,2)/2
    qiDqj(2,4) = qiDqj(4,2)
    If ( xqi(1) < xqi(2) ) Then
      qiDqj(1,2) = xqi(1) + dqiqj(2,4)
      xmax = xqi(1)
    Else
      qiDqj(1,2) = xqi(2) + dqiqj(1,4)
      xmax = xqi(2)
    End If
    qiDqj(1,2) = qiDqj(1,2)/2
    qiDqj(2,1) = qiDqj(1,2)
    qiDqj(1,1) = xqi(1)
    qiDqj(2,2) = xqi(2)
    qiDqj(4,4) = xqi(4)
  End If
 !  #] groundwork:
 !  #[ calculations:
 !
  del2 = ffdel2(qiDqj,1,2,4)
  If ( ldot ) fdel2 = del2
 !
 !	the case del2=0 is hopeless - this is really a two-point function
 !
  If ( del2 == 0 ) Then
    Call WriteLfError(38)
    Return
  End If
 !
 !	we cannot yet handle the complex case
 !
  If ( del2 > 0 ) Then
    Call WriteLfError(39)
    Return
  End If
 !
  sdel2 = isgnal*Sqrt(-del2)
 !
  Call ffxc0j(cs,ipi12,sdel2,clogi,ilogi,xqi,dqiqj,qiDqj, lambda2,3,ier)
 !  #] calculations:
 !  #[ sum:
 !
 !	Sum
 !
  xmax = 0
  csum = 0
  If ( .Not.lsmug ) Then
    n = 10
  Else
    n = 15
  End If
  Do i=1,n
    csum = csum + cs(i)
    xmax = Max(xmax,absc(csum))
  End Do
  csum = csum + ipi12*Pi2o12
  cc0 = -csum*Real(1/(2*sdel2),dp)
 !  #] sum:
  End Subroutine ffxc0i

 Subroutine ffxc0j(cs,ipi12,sdel2i,clogi,ilogi, &
                  xpi,dpipj,piDpj,delta,npoin,ier)
 !------------------------------------------------------------------------
 !	Calculates the infrared finite part of a infrared divergent
 !	threepoint function with the factor ipi^2.
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ipi12,ilogi(3),npoin,ier
  Complex(dp) :: cs(15),clogi(3)
  Real(dp) :: xpi(6),dpipj(6,6),piDpj(6,6),delta,sdel2i

  Integer :: i,ieps,ieps1,n
  Complex(dp) :: clog1,clog2,cdum(2),cel3,cdyzm,cdyzp,cli,chulp, &
    carg1,carg2,chulp1
  Real(dp) :: del2,zm,zp,zm1,zp1,sdel2,hulp,dum(3),dyzp,dyzm,wm,wp,arg1 &
       & ,arg2,del3
 !
 !  #[ get determinants, roots, ieps:
 !
  If ( lsmug ) Then
    del3 = (- Real(xpi(1),dp)*Real(cmipj(2,2),dp)**2  &
      &  - Real(xpi(2),dp)*Real(cmipj(1,3),dp)**2     &
      &  + 2*Real(piDpj(1,2),dp)*Real(cmipj(2,2),dp)*Real(cmipj(1,3),dp) )/4
      If ( nschem >= 3 ) Then
        cel3 = (- Real(xpi(1),dp)*cmipj(2,2)**2- Real(xpi(2),dp)*cmipj(1,3)**2&
           & + 2*Real(piDpj(1,2),dp)*cmipj(2,2)*cmipj(1,3) )/4
      Else
        cel3 = Real(del3,dp)
      End If
  End If
  del2 = -sdel2i**2
 !
 !	the routine as it stands can not handle sdel2<0.
 !	the simplest solution seems to be to switch to sdel2>0 for
 !	the time being - we calculate a complete 3point function so it
 !	should not be a problem (just a sign).  Of course this spoils a
 !	good check on the correctness.
 !
  sdel2 = Abs(sdel2i)
 !
  If ( xpi(4)==0 ) Then
    zm = xpi(2)/dpipj(2,1)
    zm1 = -xpi(1)/dpipj(2,1)
  Else
      Call Roots(xpi(4),piDpj(4,2),xpi(2),sdel2,zm,zp)
      If ( dpipj(1,2) /= 0 ) Then
        Call Roots(xpi(4),-piDpj(4,1),xpi(1),sdel2,zp1,zm1)
      Else
       zm1 = zp
       zp1 = zm
      End If
  End If
 !	imag sign ok 30-oct-1989
  ieps = -1
  If ( xpi(4)/=0 ) dyzp = -2*sdel2/xpi(4)
 !
 !  #] get determinants, roots, ieps:
 !  #[ the finite+divergent S1:
 !
  If ( xpi(4)/=0 ) Then
      Call ffcxr(cs(1),ipi12,zm,zm1,zp,zp1,dyzp, &
                .False.,0._DP,0._DP,0._DP,.False.,dum,ieps,ier)
  End If
 !
 !	Next the divergent piece
 !
  If ( .Not.lsmug ) Then
 !
 !	    Here we dropped the term log(lam/delta)*log(-zm/zm1)
 !
    If ( Abs(zm1) > 1/xloss ) Then
      clog1 = Log1MinusX(1/zm1)
    Else If ( zm/=0 ) Then
      clog1 = zxfflg(-zm/zm1,-2,0._DP)
    Else
      Call WriteLfError(40)
      Return
    End If
    hulp = zm*zm1*4*del2/delta**2
 !
 !	    14-jan-1994:   Do not count when this is small, this was 
 !	    meant to be so by the user carefully adjusting delta
 !
    If ( hulp==0 )   Call WriteLfError(40)
    clog2 = zxfflg(hulp,2,0._DP)
    cs(8) = -clog1*clog2/2
  Else
 !
 !	    checked 4-aug-1992, but found Yet Another Bug 30-sep-1992
 !
    cdyzm = cel3*Real(1/(-2*sdel2*del2),dp)
    dyzm = del3/(-2*sdel2*del2)
    carg1 = +cdyzm*Real(1/zm,dp)
    arg1 = +dyzm/zm
    clog1 = zfflog(-carg1,+ieps,Cmplx(Real(zm,dp),0._dp,dp))
    If (Aimag(cdyzm) < 0 .And. arg1 > 0 ) Then
      clog1 = clog1 - c2ipi
 !		ier = ier + 50
    End If
    cs(8) = -clog1**2/2
    carg2 = -cdyzm*Real(1/zm1,dp)
    arg2 = -dyzm/zm1
    clog2 = zfflog(-carg2,ieps,Cmplx(Real(-zm1,dp),0._dp,dp))
    If ( Aimag(cdyzm) < 0 .And. arg2 > 0 )  clog2 = clog2 + c2ipi
    cs(9) = +clog2**2/2
  End If
 !  #] the finite+divergent S1:
 !  #[ log(1) for npoin=4:
  If ( npoin == 4 ) Then
    If ( ilogi(1) == -999 ) Then
      If ( .Not.lsmug ) Then
        hulp = xpi(4)*delta/(4*del2)
        If ( hulp==0 )   Call WriteLfError(40)
        clogi(1) = -zxfflg(Abs(hulp),0,0._DP)
        If ( hulp < 0 ) Then
          If ( xpi(4) > 0 ) Then
            ilogi(1) = -1
          Else
            ilogi(1) = +1
          End If
        Else
          ilogi(1) = 0
        End If
      Else
        If ( xpi(4)==0 ) Then
          Print *,'ffxc0i: cannot handle t=0 yet, sorry'
          Print *,'Please regularize with a small mass'
          Call TerminateProgram()
        End If
        chulp = -cdyzm*Real(1/dyzp,dp)
        chulp1 = 1+chulp
        If ( absc(chulp1) < xloss ) Then
          Write(ErrCan,*) &
           & "zxfflg: warning: taking log of number close to 1, must be cured."
          Write(ErrCan,*) "Calling routine is ffxc0j."
        End If
        Call ffxclg(clogi(1),ilogi(1),chulp,chulp1,dyzp, ier)
      End If
    End If
  End If
 !  #] log(1) for npoin=4:
 !  #[ the log(lam) Si:
  If ( .Not.lsmug ) Then
 !
 !	    Next the divergent S_i (easy).
 !	    The term -2*log(lam/delta)*log(xpi(2)/xpi(1)) has been discarded
 !	    with lam the photon mass (regulator).
 !	    If delta = sqrt(xpi(1)*xpi(2)) the terms cancel as well
 !
    If ( dpipj(1,2)/=0 .And. xloss*Abs(xpi(1)*xpi(2)-delta**2) &
           >precx*delta**2 ) Then
      If ( xpi(1) /= delta ) Then
        If ( xpi(1)==0 )   Call WriteLfError(40)
        cs(9) = -zxfflg(xpi(1)/delta,0,0._DP)**2 /4
      End If
      If ( xpi(2) /= delta ) Then
        If ( xpi(2)==0 )   Call WriteLfError(40)
        cs(10) = zxfflg(xpi(2)/delta,0,0._DP)**2 /4
      End If
    End If
 !  #] the log(lam) Si:
 !  #[ the logs for A_i<0:
    If ( npoin==4 ) Then
      clogi(2) = 0
      ilogi(2) = 0
      clogi(3) = 0
      ilogi(3) = 0
    End If
 !  #] the logs for A_i<0:
 !  #[ the off-shell S3:
  Else
 !
 !	    the divergent terms in the offshell regulator scheme - not
 !	    quite as easy
 !	    wm = p3.p2/sqrtdel - 1 = -s1.s2/sqrtdel - 1
 !	    wp = p3.p2/sqrtdel + 1 = -s1.s2/sqrtdel + 1
 !	    Note that we took the choice sdel2<0 in S1 when
 !	    \delta^{p1 s2}_{p1 p2} < 0 by using yp=zm
 !
    wm = -1 - piDpj(1,2)/sdel2
    wp = wm + 2
    If ( Abs(wm) < Abs(wp) ) Then
      wm = -xpi(5)*xpi(6)/(del2*wp)
    Else
      wp = -xpi(5)*xpi(6)/(del2*wm)
    End If
 !
 !	    the im sign
 !
    If ( -Real(cmipj(1,3),dp) > 0 ) Then
      ieps = -1
    Else
      ieps = +1
    End If
 !
    If ( nschem < 3 .Or. Aimag(cmipj(1,3))==0 .And.Aimag(cmipj(2,2))==0 ) Then
 !  #[		real case:
 !
 !		first z-,z+
 !
      dyzp = -Real(cmipj(1,3),dp)*Real(wm,dp)/(2*Real(xpi(6),dp))  &
        &    -Real(cmipj(2,2),dp)/(2*Real(sdel2,dp))
      dyzm = -Real(cmipj(1,3),dp)*Real(wp,dp)/(2*Real(xpi(6),dp)) &
        &    -Real(cmipj(2,2),dp)/(2*Real(sdel2,dp))
 !
 !		the (di)logs
 !
      clog1 = zxfflg(-dyzp,-ieps,1._DP)
      cs(10) = -clog1**2/2
      ipi12 = ipi12 - 4
      clog2 = zxfflg(-dyzm,+ieps,1._DP)
      cs(11) = -clog2**2/2
      ipi12 = ipi12 - 2
      hulp = dyzp/dyzm
      If ( dyzp < 0 ) Then
        ieps1 = -ieps
      Else
        ieps1 = +ieps
      End If
      Call ffzxdl(cli,i,cdum(1),hulp,+ieps1,ier)
      cs(12) = -cli
      ipi12 = ipi12 - i
 !
 !		the log for npoin=4
 !
      If ( npoin==4 ) Then
        If ( ilogi(3) == -999 ) Then
          If ( Real(cmipj(1,3),dp) == 0 ) Then
            chulp = -1
            chulp1 = 0
          Else If ( dyzp < dyzm ) Then
            chulp = -dyzm/dyzp
            chulp1 = +Real(cmipj(1,3),dp)/Real(xpi(6)*dyzp,dp)
          Else
            chulp = -dyzp/dyzm
            chulp1 = -Real(cmipj(1,3),dp)/Real(xpi(6)*dyzm,dp)
          End If
          Call ffxclg(clogi(3),ilogi(3),chulp,chulp1,dyzp, ier)
        End If
      End If
 !  #]		real case:
    Else
 !  #[		complex case:
 !
 !		first z+
 !
      cdyzp = -cmipj(1,3)*Real(wm,dp)/(2*Real(xpi(6),dp)) &
          &   -cmipj(2,2)/(2*Real(sdel2,dp))
      clog1 = zfflog(-cdyzp,-ieps,cone)
      If ( ieps*Aimag(cdyzp)<0.And.Real(cdyzp,dp)>0) clog1 = clog1 - ieps*c2ipi
      cs(10) = -clog1**2/2
      ipi12 = ipi12 - 4
 !
 !		now z-
 !
      cdyzm = -cmipj(1,3)*Real(wp,dp)/(2*Real(xpi(6),dp)) &
          &   -cmipj(2,2)/(2*Real(sdel2,dp))
      clog2 = zfflog(-cdyzm,+ieps,cone)
      If ( ieps*Aimag(cdyzm)>0.And.Real(cdyzm)>0 ) Then
        clog2 = clog2 + ieps*c2ipi
      End If
      cs(11) = -clog2**2/2
      ipi12 = ipi12 - 2
 !
 !		the dilog
 !
      chulp = cdyzp/cdyzm
      hulp = Real(cdyzp,dp)/Real(cdyzm,dp)
      If ( Real(cdyzp) < 0 ) Then
        ieps1 = -ieps
      Else
        ieps1 = +ieps
      End If
      If ( Aimag(chulp) == 0 ) Then
        hulp = Real(chulp,dp)
        Call ffzxdl(cli,i,cdum(1),hulp,+ieps1,ier)
      Else
        Call ffzzdl(cli,i,cdum(1),chulp,ier)
        If ( hulp>1 .And. ieps1*Aimag(chulp)<0 ) Then
          cli = cli + ieps1*c2ipi*zfflog(chulp,0,czero)
        End If
      End If
      cs(12) = -cli
      ipi12 = ipi12 - i
 !
 !		the log for npoin=4
 !
      If ( npoin==4 ) Then
        If ( ilogi(3) == -999 ) Then
          If ( cmipj(1,3) == 0 ) Then
            chulp = -1
            chulp1 = 0
          Else If ( Real(cdyzp,dp) < Real(cdyzm,dp) ) Then
            chulp = -cdyzm/cdyzp
            chulp1 = +cmipj(1,3)/cdyzp*Real(1/xpi(6),dp)
          Else
            chulp = -cdyzp/cdyzm
            chulp1 = -cmipj(1,3)/cdyzm*Real(1/xpi(6),dp)
          End If
          dyzp = Real(cdyzp,dp)
          Call ffxclg(clogi(3),ilogi(3),chulp,chulp1,dyzp, ier)
        End If
      End If
 !  #]		complex case:
    End If
 !  #] the off-shell S3:
 !  #[ the off-shell S2:
 !
 !	    the im sign
 !
    If ( -Real(cmipj(2,2)) > 0 ) Then
      ieps = -1
    Else
      ieps = +1
    End If
 !
    If ( nschem < 3 ) Then
 !  #[		real case:
 !
 !		first z-
 !
      dyzm = -Real(cmipj(2,2),dp)*Real(wp,dp)/(2*Real(xpi(5),dp)) &
        &    -Real(cmipj(1,3),dp)/(2*Real(sdel2,dp))
      clog1 = zxfflg(+dyzm,-ieps,1._DP)
      cs(13) = +clog1**2/2
      ipi12 = ipi12 + 4
 !
 !		now z+
 !
      dyzp = -Real(cmipj(2,2),dp)*Real(wm,dp)/(2*Real(xpi(5),dp)) &
        &    -Real(cmipj(1,3),dp)/(2*Real(sdel2,dp))
      clog2 = zxfflg(+dyzp,+ieps,1._DP)
      cs(14) = +clog2**2/2
      ipi12 = ipi12 + 2
      hulp = dyzm/dyzp
      If ( dyzm < 0 ) Then
        ieps1 = -ieps
      Else
        ieps1 = +ieps
      End If
      Call ffzxdl(cli,i,cdum(1),hulp,-ieps1,ier)
      cs(15) = +cli
      ipi12 = ipi12 + i
 !
 !		the log for npoin=4
 !
      If ( npoin==4 ) Then
        If ( ilogi(2) == -999 ) Then
          If ( Real(cmipj(2,2)) == 0 ) Then
            chulp = -1
            chulp1 = 0
          Else If ( dyzp < dyzm ) Then
            chulp = -dyzm/dyzp
            chulp1 = +Real(cmipj(2,2),dp)/Real(xpi(5)*dyzp,dp)
          Else If ( dyzp > dyzm ) Then
            chulp = -dyzp/dyzm
            chulp1 = -Real(cmipj(2,2),dp)/Real(xpi(5)*dyzm,dp)
          End If
          Call ffxclg(clogi(2),ilogi(2),chulp,chulp1,dyzp, ier)
        End If
      End If
 !  #]		real case:
    Else
 !  #[		complex case:
 !
 !		first z-
 !
      cdyzm = -cmipj(2,2)*Real(wp,dp)/(2*Real(xpi(5),dp))  &
          &   -cmipj(1,3)/(2*Real(sdel2,dp))
      clog1 = zfflog(+cdyzm,-ieps,cone)
      If ( Real(cdyzm)<0.And.ieps*Aimag(cdyzm)>0 ) clog1 = clog1 - ieps*c2ipi
      cs(13) = +clog1**2/2
      ipi12 = ipi12 + 4
 !
 !		now z+
 !
      cdyzp = -cmipj(2,2)*Real(wm,dp)/(2*Real(xpi(5),dp)) &
          &   -cmipj(1,3)/(2*Real(sdel2,dp))
      clog2 = zfflog(+cdyzp,+ieps,cone)
      If ( Real(cdyzp)<0.And.ieps*Aimag(cdyzp)<0 )  clog2 = clog2 + ieps*c2ipi
      cs(14) = +clog2**2/2
      ipi12 = ipi12 + 2
 !		
 !		and ghe dilog
 !		
      chulp = cdyzm/cdyzp
      hulp = Real(dyzm,dp)/Real(dyzp,dp)
      If ( Real(cdyzm) < 0 ) Then
        ieps1 = -ieps
      Else
        ieps1 = +ieps
      End If
      If ( Aimag(chulp ) == 0 ) Then
        hulp = Real(chulp,dp)
        Call ffzxdl(cli,i,cdum(1),hulp,-ieps1,ier)
      Else
        Call ffzzdl(cli,i,cdum(1),chulp,ier)
        If ( hulp>1 .And. ieps1*Aimag(chulp)>0 ) Then
          cli = cli - ieps1*c2ipi*zfflog(chulp,0,czero)
        End If
      End If
      cs(15) = +cli
      ipi12 = ipi12 + i
 !
 !		the log for npoin=4
 !
      If ( npoin==4 ) Then
        If ( ilogi(2) == -999 ) Then
          If ( cmipj(2,2) == 0 ) Then
            chulp = -1
            chulp1 = 0
          Else If ( Real(cdyzp,dp) < Real(cdyzm,dp) ) Then
            chulp = -cdyzm/cdyzp
            chulp1 = +cmipj(2,2)/cdyzp*Real(1/xpi(5),dp)
          Else If ( Real(cdyzp,dp) > Real(cdyzm,dp) ) Then
            chulp = -cdyzp/cdyzm
            chulp1 = -cmipj(2,2)/cdyzm*Real(1/xpi(5),dp)
          End If
          dyzp = Real(cdyzp,dp)
          Call ffxclg(clogi(2),ilogi(2),chulp,chulp1,dyzp, ier)
        End If
      End If
 !  #]		complex case:
    End If
  End If
 !  #] the off-shell S2:
 !  #[ sdel2<0!:
  If ( sdel2i>0 .Neqv. xpi(4)==0.And.xpi(1)>xpi(2) ) Then
    If ( .Not.lsmug ) Then
      n = 10
    Else
      n = 15
    End If
    Do i=1,n
      cs(i) = -cs(i)
    End Do
    ipi12 = -ipi12
    If ( npoin==4 ) Then
      Do i=1,3
        ilogi(i) = -ilogi(i)
        clogi(i) = -clogi(i)
      End Do
    End If
  End If
 !  #] sdel2<0!:
 End Subroutine ffxc0j

 Subroutine ffxc0p(cs3,ipi12,isoort,clogi,ilogi,xpi,dpipj,piDpj, &
                  sdel2,del2s,etalam,etami,delpsi,alph,npoin)
 !------------------------------------------------------------------------
 !	calculates the threepoint function closely following
 !	recipe in 't Hooft & Veltman, NP B(183) 1979.
 !	Bjorken and Drell metric is used nowadays!
 !
 !	    p2	^ |
 !		| |
 !		 / \
 !	      m2/   \m3
 !	p1     /     \	p3
 !	<-    /  m1   \ ->
 !	------------------------
 !
 !	Input:	xpi(1-3)     (real)	pi squared
 !		xpi(4-6)     (real)	internal mass squared
 !		dpipj(6,6)   (real)	xpi(i)-xpi(j)
 !		piDpj(6,6)   (real)	pi(i).pi(j)
 !		sdel2	     (real)	sqrt(delta_{p_1 p_2}^{p_1 p_2})
 !		del2s(3)     (real)	delta_{p_i s_i}^{p_i s_i}
 !		etalam	     (real)	delta_{s_1 s_2 s_3}^{s_1 s_2 s_3}
 !					  /delta_{p_1 p_2}^{p_1 p_2}
 !		etami(6)     (real)	m_i^2 - etalam
 !		alph(3)	     (real)	alph(1)=alpha, alph(3)=1-alpha
 !
 !	Output: cs3(80)	     (complex)	C0, not yet summed.
 !		ipi12(8)     (Integer)  factors pi^2/12, not yet summed
 !		slam	     (complex)	lambda(p1,p2,p3).
 !		isoort(8)    (Integer)	indication of he method used
 !		clogi(3)     (complex)	log(-dyz(2,1,i)/dyz(2,2,i))
 !		ilogi(3)     (Integer)	factors i*pi in this
 !		ier	     (Integer)	number of digits inaccurate in
 !					answer
 !
 !	Calls:	ffroot,ffxxyz,ffcxyz,ffdwz,ffcdwz,
 !		ffcxs3,ffcs3,ffcxs4,ffcs4
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ipi12(8),isoort(8),ilogi(3),npoin,ier
  Complex(dp) :: cs3(80),clogi(3)
  Real(dp) :: xpi(6),dpipj(6,6),piDpj(6,6),sdel2,del2s(3), &
    etalam,etami(6),delpsi(3),alph(3)

  Integer :: i,j,k,m,ip,jsoort(8),ierw,iw,ier0,ier1,irota, ilogip(3)
  Logical :: l4,lcompl,lcpi,l4pos
  Complex(dp) :: cs,calph(3),csdl2i(3),csdel2
  Complex(dp) :: cy(4,3),cz(4,3),cw(4,3),cdyz(2,2,3),cdwy(2,2,3), &
    cdwz(2,2,3),cd2yzz(3),cd2yww(3)
  Complex(dp) :: cpi(6),cdpipj(6,6),cpiDpj(6,6),clogip(3)
  Real(dp) :: y(4,3),z(4,3),w(4,3),dyz(2,2,3),dwy(2,2,3), &
    dwz(2,2,3),d2yzz(3),d2yww(3),dy2z(4,3)
  Real(dp) :: sdel2i(3),s1,s2
  Real(dp) :: s,xqi(6),dqiqj(6,6),qiDqj(6,6)
 !
 !  #[ IR case:
 !
 !	but only the off-shell regulator case - the log(lam) has been
 !	caught before
 !
  ier = 0
  If ( lsmug ) Then
    Do i=1,3
      If ( xpi(i) == 0 ) Then
        j = Mod(i,3)+1
        k = Mod(j,3)+1
        If ( piDpj(i,j)==0 .And. piDpj(i,k)==0 ) Then
          Call ffrot3(irota,xqi,dqiqj,qiDqj, xpi,dpipj,piDpj,3,4)
          If ( npoin==4 ) Call ffrt3p(clogip,ilogip, irota,clogi,ilogi,+1)
          Call ffxc0j(cs3(1),ipi12(1),sdel2,clogip,ilogip, xqi,dqiqj,qiDqj &
                       ,0._DP,4,ier)
          If ( npoin==4 ) Call ffrt3p(clogi,ilogi,irota, clogip,ilogip,-1)
          Return
        End If
      End If
    End Do
  End If
 !  #] IR case:
 !  #[ get roots etc:
 !  #[   get z-roots:
 !	  If ( npoin .eq. 3 ) then
  l4pos = l4also
 !	  Else
 !	    l4pos = .FALSE.
 !	  End If
  lcompl = .False.
  ier1 = ier
  Do i=1,3
 !
 !	get roots (y,z,w) and flag what to   Do: 0=nothing, 1=normal,
 !	-1=complex
 !
   ip = i+3
 !	    first get the roots
   ier0 = ier
   If ( del2s(i) <= 0 ) Then
 !		real case
     sdel2i(i) = Sqrt(-del2s(i))
     csdl2i(i) = sdel2i(i)
 !		then handle the special case Si = 0
      If ( xpi(ip) == 0 ) Then
        If ( i == 1 .And. alph(3) == 0 .Or.i == 3 .And. alph(1) == 0 ) Then
          isoort(2*i-1) = 0
          isoort(2*i) = 0
          l4pos = .False.
          Cycle
        End If
      End If
      Call ffxxyz(y(1,i),z(1,i),dyz(1,1,i),d2yzz(i),dy2z(1,i), &
                 i,sdel2,sdel2i(i),etalam,etami,delpsi(i),xpi, &
                 dpipj,piDpj,isoort(2*i-1),.False.,6,ier0)
   Else
 !		complex case
    sdel2i(i) = Sqrt(del2s(i))
    csdl2i(i) = Cmplx(0._DP,sdel2i(i),dp)
    lcompl = .True.
    Call ffcxyz(cy(1,i),cz(1,i),cdyz(1,1,i),cd2yzz(i),i, &
      sdel2,sdel2i(i),etami,delpsi(i),xpi, piDpj,isoort(2*i-1),.False.,6,ier0)
   End If
   ier1 = Max(ier1,ier0)
  End Do

  ier = ier1
 !  #]   get z-roots:
 !  #[   get w-roots:
 !
 !	get w's:
 !
  ierw = ier
  l4 = .False.
  lcpi = .False.
  If ( isoort(4) == 0 ) Then
 !	    no error message; just bail out
    Write(ErrCan,*) "unreliable numerical result in ffxc0p"
    ierw = ierw + 100

  Else
  outer: Do iw = 1,3,2
    If ( .Not. l4pos .Or. alph(4-iw) == 0 ) Then
      jsoort(2*iw-1) = 0
      jsoort(2*iw) = 0
      l4pos = .False.
    Else
     If ( isoort(4) > 0 .And. isoort(2*iw) >= 0 ) Then
      jsoort(2*iw-1) = 1
      jsoort(2*iw) = 1
      d2yww(iw) = -d2yzz(2)/alph(4-iw)
      Do j=1,2
        w(j+iw-1,iw) = z(j+3-iw,2)/alph(4-iw)
        w(j+3-iw,iw) = 1 - w(j+iw-1,iw)
        If ( Abs(w(j+3-iw,iw)) < xloss ) Then
          s = z(j+iw-1,2) - alph(iw)
          If ( Abs(s) < xloss*alph(iw) ) Then
            ierw = ierw + 15
            Cycle outer
          End If
          w(j+3-iw,iw) = s/alph(4-iw)
        End If
        dwy(j,2,iw) = dyz(2,j,2)/alph(4-iw)
        Do i=1,2
          dwz(j,i,iw) = w(j,iw) - z(i,iw)
          If ( Abs(dwz(j,i,iw)) >= xloss*Abs(w(j,iw)) ) Cycle
          dwz(j,i,iw) = z(i+2,iw) - w(j+2,iw)
          If ( Abs(dwz(j,i,iw)) >= xloss*Abs(w(j+2,iw)) ) Cycle
          dwz(j,i,iw) = dwy(j,2,iw) + dyz(2,i,iw)
          If ( Abs(dwz(j,i,iw)) >= xloss*Abs(dwy(j,2,iw)) ) Cycle
          l4 = .True.
          Call ffdwz(dwz(1,1,iw),z(1,iw),j,i,iw, &
            alph(1),alph(3),xpi,dpipj,piDpj,sdel2i,6,ierw)
        End Do
      End Do
     Else
 !	    convert to complex ...
      jsoort(2*iw-1) = -10
      jsoort(2*iw) = -10
      If ( isoort(4)>=0 .And. (iw==1 .Or. isoort(2)>=0) ) Then
        cd2yzz(2) = d2yzz(2)
        Do i=1,4
          cy(i,2) = y(i,2)
          cz(i,2) = z(i,2)
        End Do
        Do i=1,2
          Do j=1,2
            cdyz(j,i,2) = dyz(j,i,2)
          End Do
        End Do
      End If
      If ( isoort(2*iw) >= 0 ) Then
        cd2yzz(iw) = d2yzz(iw)
        Do i=1,4
          cy(i,iw) = y(i,iw)
          cz(i,iw) = z(i,iw)
        End Do
        Do i=1,2
          Do j=1,2
            cdyz(j,i,iw) = dyz(j,i,iw)
          End Do
        End Do
      End If
      cd2yww(iw) = -cd2yzz(2)/Real(alph(4-iw),dp)
      Do j=1,2
        cw(j+iw-1,iw) = cz(j+3-iw,2)/Real(alph(4-iw),dp)
        cw(j+3-iw,iw) = 1 - cw(j+iw-1,iw)
        If ( absc(cw(j+3-iw,iw)) < xloss ) Then
          cs = cz(j+iw-1,2) - Real(alph(iw),dp)
          If ( absc(cs) < xloss*alph(iw) ) ierw = ierw + 15
          cw(j+3-iw,iw) = cs/Real(alph(4-iw),dp)
        End If
        cdwy(j,2,iw) = cdyz(2,j,2)/Real(alph(4-iw),dp)
        Do i=1,2
          cdwz(j,i,iw) = cw(j,iw) - cz(i,iw)
          If ( absc(cdwz(j,i,iw)) >= xloss*absc(cw(j,iw)) ) Cycle
          cdwz(j,i,iw) = cz(i+2,iw) - cw(j+2,iw)
          If ( absc(cdwz(j,i,iw)) >= xloss*absc(cw(j+2,iw))) Cycle
          cdwz(j,i,iw) = cdwy(j,2,iw) + cdyz(2,i,iw)
          If ( absc(cdwz(j,i,iw))>=xloss*absc(cdwy(j,2,iw))) Cycle
          l4 = .True.
          If ( .Not. lcpi ) Then
            lcpi = .True.
            calph(1) = alph(1)
            calph(3) = alph(3)
            csdel2 = sdel2
            Do k=1,6
              cpi(k) = xpi(k)
              Do m=1,6
                cdpipj(m,k) = dpipj(m,k)
                cpiDpj(m,k) = piDpj(m,k)
              End Do
            End Do
          End If
          Call ffcdwz(cdwz(1,1,iw),cz(1,iw),j,i,iw, &
            calph(1),calph(3),cpi,cdpipj,cpiDpj,csdl2i, csdel2,6,ierw)
        End Do
      End Do
     End If
    End If
   End Do outer
  End If

  ierw = ierw-ier
 !  #]   get w-roots:
 !  #[   which case:
  If ( l4 ) Then
 !	    21-aug-1995.  added check for isoort(2*i-1).eq.0 to avoid
 !	    undefined variables etc in ffdcs, ffdcrr.  They should be
 !	    able to handle this, but are not (yet?)
    If ( ierw >= 1 .Or. isoort(1)==0 .Or. isoort(3)==0.Or. isoort(5)==0 ) Then
      l4pos = .False.
    Else
      ier = ier + ierw
    End If
  End If
 !  #]   which case:
 !  #] get roots etc:
 !  #[ logarithms for 4point function:
  If ( npoin == 4 ) Then
    Do i = 1,3
      If ( ilogi(i) /= -999 ) Cycle
      If ( isoort(2*i) > 0 .And.isoort(2*i-1) >= 0 ) Then
        s1 = -dyz(2,1,i)/dyz(2,2,i)
        If ( Abs(s1-1) < xloss ) Then
          clogi(i) = Log1MinusX(d2yzz(i)/dyz(2,2,i))
          ilogi(i) = 0
        Else
          If ( Abs(s1+1) < xloss ) Then
            clogi(i) = Log1MinusX(-2*sdel2i(i)/(xpi(i+3)*dyz(2,2,i)))
          Else
            clogi(i) = zxfflg(Abs(s1),0,0._DP)
          End If
          If ( dyz(2,2,i)>0 .And. dyz(2,1,i)>0 ) Then
            ilogi(i) = -1
          Else If ( dyz(2,1,i)<0 .And. dyz(2,2,i)<0) Then
            ilogi(i) = +1
          Else
            ilogi(i) = 0
          End If
        End If
      Else If ( isoort(2*i-1) < 0 ) Then
 !		for stability split the unit circle up in 4*pi/2
 !		(this may have to be improved to 8*pi/4...)
        ier0 = 0
        If ( Real(cdyz(2,1,i),dp) > Aimag(cdyz(2,1,i)) ) Then
          s = 2*Atan2(Aimag(cdyz(2,1,i)),Real(cdyz(2,1,i),dp))
          clogi(i) = Cmplx(0._DP,s,dp)
          ilogi(i) = -1
        Else If ( Real(cdyz(2,1,i),dp) < -Aimag(cdyz(2,1,i))) Then
          If ( Aimag(cdyz(2,1,i)) == 0 )  Call WriteLfError(36)
          s = 2*Atan2(-Aimag(cdyz(2,1,i)),-Real(cdyz(2,1,i),dp))
          clogi(i) = Cmplx(0._DP,s,dp)
          ilogi(i) = 1
        Else
          s1 = -Real(cdyz(2,1,i),dp)
          s2 = Aimag(cdyz(2,1,i))
          s = 2*Atan2(s1,s2)
          clogi(i) = Cmplx(0._DP,s,dp)
          ilogi(i) = 0
        End If
      End If
    End Do
 !	An algorithm to obtain the sum of two small logarithms more
 !	accurately has been put in ffcc0p, not yet here
  End If
 !  #] logarithms for 4point function:
 !  #[ real case integrals:
  ier1 = ier
  If ( .Not. lcompl ) Then
    If ( .Not. l4 .Or. .Not. l4pos ) Then
 !		normal case
      Do i=1,3
        j = 2*i-1
        If ( isoort(j) /= 0 ) Then
          ier0 = ier
          Call ffcxs3(cs3(20*i-19),ipi12(j),y(1,i),z(1,i), &
              dyz(1,1,i),d2yzz(i),dy2z(1,i),xpi,piDpj, i,6,isoort(j),ier0)
          ier1 = Max(ier1,ier0)
        End If
      End Do
      isoort(7) = 0
      isoort(8) = 0
    Else
      Do i=1,3,2
        j = 2*i-1
        isoort(j+2) = jsoort(j)
        isoort(j+3) = jsoort(j+1)
        ier0 = ier
        Call ffcxs4(cs3(20*i-19),ipi12(j),w(1,i),y(1,i), &
            z(1,i),dwy(1,1,i),dwz(1,1,i),dyz(1,1,i), &
            d2yww(i),d2yzz(i),xpi,piDpj,i,6,isoort(j),ier0)
        ier1 = Max(ier1,ier0)
      End Do
    End If
 !  #] real case integrals:
 !  #[ complex case integrals:
  Else
 !	    convert xpi
    If ( .Not.lcpi ) Then
      Do i=1,6
        cpi(i) = xpi(i)
      End Do
      End If
      If ( .Not. l4 .Or. .Not. l4pos ) Then
 !		normal case
        Do i=1,3
          j = 2*i-1
          ier0 = ier
          If ( isoort(j) > 0 ) Then
            Call ffcxs3(cs3(20*i-19),ipi12(2*i-1),y(1,i),z(1,i),dyz(1,1,i) &
                   ,d2yzz(i),dy2z(1,i), xpi,piDpj,i,6,isoort(j),ier0)
          Else If( isoort(j) /= 0 ) Then
            Call ffcs3(cs3(20*i-19),ipi12(2*i-1),cy(1,i),cz(1,i),cdyz(1,1,i) &
                  ,cd2yzz(i),cpi, cpiDpj,i,6,isoort(j))
          End If
          ier1 = Max(ier1,ier0)
        End Do
        isoort(7) = 0
        isoort(8) = 0
      Else
        isoort(3) = jsoort(1)
        isoort(4) = jsoort(2)
        ier0 = ier
        If ( isoort(1) > 0 .And. isoort(3) > 0 ) Then
          Call ffcxs4(cs3(1),ipi12(1),w(1,1),y(1,1), &
             z(1,1),dwy(1,1,1),dwz(1,1,1),dyz(1,1,1), &
            d2yww(1),d2yzz(1),xpi,piDpj,1,6,isoort(1),ier0)
        Else
          Call ffcs4(cs3(1),ipi12(1),cw(1,1),cy(1,1),cz(1,1),cdwy(1,1,1) &
              ,cdwz(1,1,1),cdyz(1,1,1), cd2yww(1),cd2yzz(1),cpi,cpiDpj, &
               Cmplx(xpi(5)*alph(3)**2,0._dp,dp),1,6,isoort(1))
        End If
        ier1 = Max(ier1,ier0)
        isoort(7) = jsoort(5)
        isoort(8) = jsoort(6)
        ier0 = ier
        If ( isoort(5) > 0 .And. isoort(7) > 0 ) Then
          Call ffcxs4(cs3(41),ipi12(5),w(1,3),y(1,3), &
               z(1,3),dwy(1,1,3),dwz(1,1,3),dyz(1,1,3), &
               d2yww(3),d2yzz(3),xpi,piDpj,3,6,isoort(5),ier0)
        Else
          Call ffcs4(cs3(41),ipi12(5),cw(1,3),cy(1,3),cz(1,3),cdwy(1,1,3) &
                    ,cdwz(1,1,3),cdyz(1,1,3), cd2yww(3),cd2yzz(3),cpi,cpiDpj, &
                     Cmplx(xpi(5)*alph(1)**2,0._dp,dp),3,6,isoort(5))
        End If
        ier1 = Max(ier1,ier0)
      End If
  End If
  ier = ier1
 !  #] complex case integrals:
  End Subroutine ffxc0p

 Subroutine ffxclg(clg,ilg,chulp,chulp1,dyzp,ier)
 !------------------------------------------------------------------------
 !	compute the extra logs for npoin=4 given chulp=-cdyzm/cdyzp	*
 !	all flagchecking has already been   Done.				*
 !
 !	Input:	chulp	(complex)	see above			*
 !		chulp1	(complex)	1+chulp (in case chulp ~ -1)	*
 !		dyzp	(real)		(real part of) y-z+ for im part	*
 !	Output:	clg	(complex)	the log				*
 !		ilg	(  Integer)	factor i*pi split off clg	*
 !------------------------------------------------------------------------
 Implicit None
  Integer :: ilg,ier
  Real(dp) :: dyzp
  Complex(dp) :: clg,chulp,chulp1

  Real(dp) :: hulp,hulp1
 !
 !  #[ work:
 !
  If ( Aimag(chulp) == 0 ) Then
    hulp = Real(chulp,dp)
    hulp1 = Real(chulp1,dp)
    If ( Abs(hulp1) < xloss ) Then
      clg = Log1MinusX(hulp1)
    Else
      clg = zxfflg(Abs(hulp),0,0._DP)
    End If
    If ( hulp < 0 ) Then
      If ( dyzp<0 ) Then
        ilg = +1
      Else
        ilg = -1
      End If
    Else
      ilg = 0
    End If
  Else
 !
 !	    may have to be improved
 !
    If ( Abs(Real(chulp1,dp))+Abs(Aimag(chulp1)) < xloss ) Then
      clg = Log1MinusX(chulp1)
    Else
      clg = zfflog(chulp,0,czero)
    End If
    ilg = 0
    If ( Real(chulp) < 0 ) Then
      If ( dyzp<0 .And. Aimag(clg)<0 ) Then
        ilg = +2
      Else If ( dyzp>0 .And. Aimag(clg)>0 ) Then
        ilg = -2
      End If
    End If
  End If
 !  #] work:
  End Subroutine ffxclg


 Subroutine ffxl22(xl22,x)
 !---------------------------------------------------------------
 !      calculates Li2(2-x) for |x|<.14 in a faster way to ~15
 !      significant figures. Note, that series starts with x^2 
 !---------------------------------------------------------------
 Implicit None
  Real(dp) :: xl22, x

  Integer :: i1, nmax
  Real(dp) :: sumI
  Real(dp), Parameter :: dilog2(28) = (/  0.25_dp, 1._dp/6._dp              &
     & , 5._dp/48._dp, 1._dp / 15._dp, 2._dp / 45._dp, 13._dp / 420._dp     &
     & , 151._dp / 6720._dp, 16._dp / 945._dp, 83._dp / 6300._dp            &
     & , 73._dp / 6930._dp, 1433._dp / 166320._dp, 647._dp / 90090._dp      &
     & , 15341._dp / 2522520._dp, 28211._dp / 5405400._dp                   &
     & , 10447._dp / 2306304._dp, 608._dp / 153153._dp                      &
     & , 19345._dp / 5513508._dp, 18181._dp / 5819814._dp                   &
     & , 130349._dp / 46558512._dp, 771079._dp / 305540235._dp              &
     & , 731957._dp / 320089770._dp, 2786599._dp / 1338557220._dp           &
     & , 122289917._dp / 64250746560._dp, 14614772._dp / 8365982625._dp     &
     & , 140001721._dp / 87006219300._dp, 134354573._dp / 90352612350._dp   &
     & , 774885169._dp / 562194032400._dp, 745984697._dp / 582272390700._dp /)

   If (Abs(x).Gt.0.14_dp) Then
    xl22 = Li2(2._dp - x)
   Else
    nmax = FindBound(x, 2, dilog2)
    sumI = x * dilog2(nmax)
    Do i1=nmax-1,1,-1
     sumI = x * (dilog2(i1) + sumI)
    End Do
    xl22 = - x * sumI
   End If
 End Subroutine ffxl22
    

 Subroutine ffxli2(xli, xlog, x_in)
 !------------------------------------------------------
 ! substitutes the original one, makes portation easier
 !------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x_in
  Real(dp), Intent(out) :: xli, xlog
  xli = Li2(x_in)
  xlog = Log1minusX(x_in)
 End Subroutine ffxli2



 Subroutine ffxxyz(y,z,dyz,d2yzz,dy2z,ivert,sdel2p,sdel2s,etalam, &
                 &  etami,delps,xpi,dpipj,piDpj,isoort,ldel2s,ns,ier)
 !----------------------------------------------------------------------
 !	calculate in a numerically stable way
 !
 !	z(1,2) = (-p(ip1).p(is2) +/- sdel2s)/xpi(ip1)
 !	y(1,2) = (-p(ip1).p(is2) +/- sdisc)/xpi(ip1)
 !			disc = del2s + etaslam*xpi(ip1)
 !
 !	y(3,4) = 1-y(1,2)
 !	z(3,4) = 1-z(1,2)
 !	dyz(i,j) = y(i) - z(j)
 !	d2yzz = y(2) - z(1) - z(2)
 !	dy2z(j) = y(2) - 2*z(j)
 !
 !	Input:	ivert		(Integer)	defines the vertex
 !		sdel2p		(real)		sqrt(lam(p1,p2,p3))/2
 !		sdel2s		(real)		sqrt(lam(p,ma,mb))/2
 !		etalam		(real)		det(si.sj)/det(pi.pj)
 !		etami(6)	(real)		si.si - etalam
 !		xpi(ns)		(real)		standard
 !		piDpj(ns,ns)	(real)		standard
 !		ns		(Integer)	dim of xpi,piDpj
 !
 !	Output:	y(4),z(4),dyz(4,4)	(real)		see above
 !
 !	Calls:	fferr,ffroot
 !----------------------------------------------------------------------
 Implicit None
  Integer :: ivert,ns,ier,isoort(2)
  Logical :: ldel2s
  Real(dp) :: y(4),z(4),dyz(2,2),d2yzz,dy2z(4), &
    sdel2p,sdel2s,etalam,etami(6),delps,xpi(ns), dpipj(ns,ns),piDpj(ns,ns)

  Integer :: i,j,n,ip1,ip2,ip3,is1,is2,is3,iwarn,ier1
  Real(dp) :: disc,hulp,s,smax,som(51),xmax
  Real(dp) :: t1,t2,t4,t5,t8,t3,t7,t9,t12,t14,t21,t23,t24, &
    t28,t6,t35,t44,t42,t36,t55,t41,t19,t59,t25,t69,t82,t75,t84,t92, &
    t31,t98,t74,t101,t89,t106,t112,t113,t13,t117,t126,t127,t129, &
    t130,t133,t128,t132,t134,t137,t139,t146,t148,t149,t153,t131, &
    t160,t171,t169,t161,t182,t168,t144,t186,t150,t208,t201,t210, &
    t219,t156,t225,t200,t228,t215,t233,t239,t240,t138,t244
 !
 !

 !  #[ set up pointers:
  If ( ldel2s .And. ivert /= 1 ) Goto 100
  is1 = ivert
  is2 = ivert+1
  If ( is2 == 4 ) is2 = 1
  is3 = ivert-1
  If ( is3 == 0 ) is3 = 3
  ip1 = is1 + 3
  ip2 = is2 + 3
  ip3 = is3 + 3
 !  #] set up pointers:
 !  #[ xk = 0:
  If ( xpi(ip1) == 0 ) Then
    isoort(2) = 0
    If ( piDpj(is1,ip1) == 0 ) Then
      isoort(1) = 0
      Return
    End If
    isoort(1) = 1
    y(1) = etami(is2) / piDpj(is1,ip1) /2
    y(2) = y(1)
    y(3) = - etami(is1) / piDpj(is1,ip1) /2
    y(4) = y(3)
    z(1) = xpi(is2) / piDpj(is1,ip1) /2
    z(2) = z(1)
    z(3) = - xpi(is1) / piDpj(is1,ip1) /2
    z(4) = z(3)
    dyz(1,1) = - etalam / piDpj(is1,ip1) /2
    dyz(1,2) = dyz(1,1)
    dyz(2,1) = dyz(1,1)
    dyz(2,2) = dyz(1,1)
    ier1 = ier
    Do i=1,3,2
      dy2z(i) = y(i) - 2*z(i)
      smax = Abs(y(i))
      dy2z(i+1) = dy2z(i)
    End Do
    ier = ier1
    Return
  End If
 !  #] xk = 0:
 !  #[ get y(1,2),z(1,2):
  If ( sdel2s == 0 ) Then
    isoort(1) = 2
    isoort(2) = 2
    z(1) = piDpj(ip1,is2)/xpi(ip1)
    z(2) = z(1)
  Else
    isoort(1) = 1
    isoort(2) = 1
    Call Roots(xpi(ip1),piDpj(ip1,is2),xpi(is2), sdel2s,z(1),z(2))
  End If
  disc = delps/sdel2p
  Call Roots(xpi(ip1),piDpj(ip1,is2),etami(is2),disc, y(1),y(2))
 !  #] get y(1,2),z(1,2):
 !  #[ get y(3,4),z(3,4):
  If ( isoort(1) == 2 ) Then
    z(3) = -piDpj(ip1,is1)/xpi(ip1)
    z(4) = z(3)
  Else
    z(3) = 1-z(1)
    z(4) = 1-z(2)
      If ( Abs(z(3)) < xloss .Or. Abs(z(4)) < xloss ) &
      Call Roots(xpi(ip1),-piDpj(ip1,is1), xpi(is1),sdel2s,z(4),z(3))
  End If
  y(3) = 1-y(1)
  y(4) = 1-y(2)
  If ( Abs(y(3)) < xloss .Or. Abs(y(4)) < xloss ) Then
    Call Roots(xpi(ip1),-piDpj(ip1,is1), etami(is1),disc,y(4),y(3))
  End If
 !  #] get y(3,4),z(3,4):
 !  #[ get dyz:
 !	Note that dyz(i,j) only exists for i,j=1,2!
  If ( isoort(1) == 2 ) Then
    dyz(2,1) = disc/xpi(ip1)
    dyz(2,2) = dyz(2,1)
  Else If ( disc > 0 .Eqv. sdel2s > 0 ) Then
    dyz(2,1) = ( disc + sdel2s )/xpi(ip1)
    dyz(2,2) = etalam/(xpi(ip1)*dyz(2,1))
  Else
    dyz(2,2) = ( disc - sdel2s )/xpi(ip1)
    dyz(2,1) = etalam/(xpi(ip1)*dyz(2,2))
  End If
  dyz(1,1) = -dyz(2,2)
  dyz(1,2) = -dyz(2,1)
  d2yzz = 2*disc/xpi(ip1)
 !
 !	these are very rarely needed, but ...
 !
  iwarn = 0
  ier1 = ier
  Do i=1,4
    j = 2*((i+1)/2)
    dy2z(i) = y(j) - 2*z(i)
    smax = Abs(y(j))
    If ( Abs(dy2z(i)) < xloss*smax ) Then
      If ( i/2 == 1 ) Then
        s = -y(j-1) - 2*sdel2s/xpi(ip1)
      Else
        s = -y(j-1) + 2*sdel2s/xpi(ip1)
      End If
      If ( Abs(y(j-1)) < smax ) Then
        dy2z(i) = s
        smax = Abs(y(j-1))
      End If
      If ( Abs(dy2z(i)) < xloss*smax ) Then
        If ( iwarn /= 0 ) Then
        Else
          iwarn = i
          xmax = smax
        End If
      End If
    End If
  End Do
  If ( iwarn /= 0 ) Then
 !
 !	    we should import the d  Ifferences, but later...
 !
    If ( Abs(dpipj(is3,ip1)) < xloss*xpi(is3) &
      .And. Abs(dpipj(is1,is2)) < xloss*Abs(xpi(ip1))) Then
 !
 !		give it another try - multiply roots (see dy2z.frm)
 !
      If ( iwarn<3 ) Then
 !prod1=
 !	som(1)=+160*xpi(ip1)*xpi(ip2)*xpi(is2)*piDpj(ip1,ip2)**2*
 !     +	dpipj(is2,is1)**2
 !	som(2)=-40*xpi(ip1)*xpi(ip2)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is2)*dpipj(is2,is1)**3
 !	som(3)=-32*xpi(ip1)*xpi(ip2)*piDpj(ip1,ip2)**2*dpipj(is2,
 !     +	is1)**3
 !	som(4)=+9*xpi(ip1)*xpi(ip2)**2*dpipj(is2,is1)**4
 !	som(5)=-128*xpi(ip1)*xpi(is2)*piDpj(ip1,ip2)**3*piDpj(ip2,
 !     +	is2)*dpipj(is2,is1)
 !	som(6)=-128*xpi(ip1)*xpi(is2)*piDpj(ip1,ip2)**4*dpipj(is2,
 !     +	is1)
 !	som(7)=+256*xpi(ip1)*xpi(is2)**2*piDpj(ip1,ip2)**4
 !	som(8)=-16*xpi(ip1)*piDpj(ip1,ip2)**2*piDpj(ip2,is2)**2*
 !     +	dpipj(is2,is1)**2
 !	som(9)=+96*xpi(ip1)*piDpj(ip1,ip2)**3*piDpj(ip2,is2)*dpipj(is2,
 !     +	is1)**2
 !	som(10)=+128*xpi(ip1)**2*xpi(ip2)*xpi(is2)*piDpj(ip1,ip2)*piDpj(
 !     +	ip2,is2)*dpipj(is2,is1)
 !	som(11)=+320*xpi(ip1)**2*xpi(ip2)*xpi(is2)*piDpj(ip1,ip2)**2*
 !     +	dpipj(is2,is1)
 !	som(12)=-512*xpi(ip1)**2*xpi(ip2)*xpi(is2)**2*piDpj(ip1,ip2)**2
 !	som(13)=-120*xpi(ip1)**2*xpi(ip2)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is2)*dpipj(is2,is1)**2
 !	som(14)=-48*xpi(ip1)**2*xpi(ip2)*piDpj(ip1,ip2)**2*dpipj(is2,
 !     +	is1)**2
 !	som(15)=+40*xpi(ip1)**2*xpi(ip2)*piDpj(ip2,is2)**2*dpipj(is2,
 !     +	is1)**2
 !	som(16)=-96*xpi(ip1)**2*xpi(ip2)**2*xpi(is2)*dpipj(is2,is1)**2
 !	som(17)=+36*xpi(ip1)**2*xpi(ip2)**2*dpipj(is2,is1)**3
 !	som(18)=+128*xpi(ip1)**2*xpi(is2)*piDpj(ip1,ip2)**2*piDpj(ip2,
 !     +	is2)**2
 !	som(19)=-128*xpi(ip1)**2*xpi(is2)*piDpj(ip1,ip2)**3*piDpj(ip2,
 !     +	is2)
 !	som(20)=-64*xpi(ip1)**2*xpi(is2)*piDpj(ip1,ip2)**4
 !	som(21)=-32*xpi(ip1)**2*piDpj(ip1,ip2)*piDpj(ip2,is2)**3*
 !     +	dpipj(is2,is1)
 !	som(22)=-32*xpi(ip1)**2*piDpj(ip1,ip2)**2*piDpj(ip2,is2)**2*
 !     +	dpipj(is2,is1)
 !	som(23)=+96*xpi(ip1)**2*piDpj(ip1,ip2)**3*piDpj(ip2,is2)*
 !     +	dpipj(is2,is1)
 !	som(24)=+128*xpi(ip1)**3*xpi(ip2)*xpi(is2)*piDpj(ip1,ip2)*piDpj(
 !     +	ip2,is2)
 !	som(25)=+160*xpi(ip1)**3*xpi(ip2)*xpi(is2)*piDpj(ip1,ip2)**2
 !	som(26)=-128*xpi(ip1)**3*xpi(ip2)*xpi(is2)*piDpj(ip2,is2)**2
 !	som(27)=+32*xpi(ip1)**3*xpi(ip2)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is1)*piDpj(ip2,is2)
 !	som(28)=-120*xpi(ip1)**3*xpi(ip2)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is2)*dpipj(is2,is1)
 !	som(29)=-32*xpi(ip1)**3*xpi(ip2)*piDpj(ip1,ip2)**2*dpipj(is2,
 !     +	is1)
 !	som(30)=-16*xpi(ip1)**3*xpi(ip2)*piDpj(ip2,is1)*piDpj(ip2,
 !     +	is2)**2
 !	som(31)=+80*xpi(ip1)**3*xpi(ip2)*piDpj(ip2,is2)**2*dpipj(is2,
 !     +	is1)
 !	som(32)=-192*xpi(ip1)**3*xpi(ip2)**2*xpi(is2)*dpipj(is2,is1)
 !	som(33)=+256*xpi(ip1)**3*xpi(ip2)**2*xpi(is2)**2
 !	som(34)=+54*xpi(ip1)**3*xpi(ip2)**2*dpipj(is2,is1)**2
 !	som(35)=-16*xpi(ip1)**3*xpi(ip3)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is1)*piDpj(ip2,is2)
 !	som(36)=+8*xpi(ip1)**3*xpi(ip3)*piDpj(ip2,is1)*piDpj(ip2,is2)**2
 !	som(37)=+16*xpi(ip1)**3*xpi(is2)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is1)*piDpj(ip2,is2)
 !	som(38)=-8*xpi(ip1)**3*xpi(is2)*piDpj(ip2,is1)*piDpj(ip2,is2)**2
 !	som(39)=-16*xpi(ip1)**3*piDpj(ip1,ip2)*piDpj(ip2,is1)*piDpj(ip2,
 !     +	is2)*dpipj(is3,ip1)
 !	som(40)=+8*xpi(ip1)**3*piDpj(ip2,is1)*piDpj(ip2,is2)**2*
 !     +	dpipj(is3,ip1)
 !	som(41)=-40*xpi(ip1)**4*xpi(ip2)*piDpj(ip1,ip2)*piDpj(ip2,is2)
 !	som(42)=-8*xpi(ip1)**4*xpi(ip2)*piDpj(ip1,ip2)**2
 !	som(43)=+40*xpi(ip1)**4*xpi(ip2)*piDpj(ip2,is2)**2
 !	som(44)=-96*xpi(ip1)**4*xpi(ip2)**2*xpi(is2)
 !	som(45)=+36*xpi(ip1)**4*xpi(ip2)**2*dpipj(is2,is1)
 !	som(46)=+9*xpi(ip1)**5*xpi(ip2)**2
 !	som(47)=-8*xpi(ip2)*piDpj(ip1,ip2)**2*dpipj(is2,is1)**4
 !	som(48)=-64*xpi(is2)*piDpj(ip1,ip2)**4*dpipj(is2,is1)**2
 !	som(49)=+32*piDpj(ip1,ip2)**3*piDpj(ip2,is2)*dpipj(is2,is1)**3
 !	print '(7g20.12)',(som(i),i=1,49)
 !
 !	optimized by Maple (see ffxxyz.map)
 !
        t1 = xpi(ip1)
        t2 = xpi(ip2)
        t3 = t1*t2
        t4 = xpi(is2)
        t5 = piDpj(ip1,ip2)
        t6 = t5**2
        t7 = t4*t6
        t8 = dpipj(is2,is1)
        t9 = t8**2
        som(1) = 160*t3*t7*t9
        t12 = piDpj(ip2,is2)
        t13 = t5*t12
        t14 = t9*t8
        som(2) = -40*t3*t13*t14
        som(3) = -32*t3*t6*t14
        t19 = t2**2
        t21 = t9**2
        som(4) = 9*t1*t19*t21
        t23 = t1*t4
        t24 = t6*t5
        t25 = t24*t12
        som(5) = -128*t23*t25*t8
        t28 = t6**2
        som(6) = -128*t23*t28*t8
        t31 = t4**2
        som(7) = 256*t1*t31*t28
        t35 = t12**2
        t36 = t35*t9
        som(8) = -16*t1*t6*t36
        som(9) = 96*t1*t24*t12*t9
        t41 = t1**2
        t42 = t41*t2
        t44 = t13*t8
        som(10) = 128*t42*t4*t44
        som(11) = 320*t42*t7*t8
        som(12) = -512*t42*t31*t6
        som(13) = -120*t42*t13*t9
        som(14) = -48*t42*t6*t9
        som(15) = 40*t42*t36
        t55 = t41*t19
        som(16) = -96*t55*t4*t9
        som(17) = 36*t55*t14
        t59 = t41*t4
        som(18) = 128*t59*t6*t35
        som(19) = -128*t59*t25
        som(20) = -64*t59*t28
        som(21) = -32*t41*t5*t35*t12*t8
        t69 = t35*t8
        som(22) = -32*t41*t6*t69
        som(23) = 96*t41*t24*t12*t8
        t74 = t41*t1
        t75 = t74*t2
        som(24) = 128*t75*t4*t5*t12
        som(25) = 160*t75*t7
        som(26) = -128*t75*t4*t35
        t82 = piDpj(ip2,is1)
        t84 = t5*t82*t12
        som(27) = 32*t75*t84
        som(28) = -120*t75*t44
        som(29) = -32*t75*t6*t8
        t89 = t82*t35
        som(30) = -16*t75*t89
        som(31) = 80*t75*t69
        t92 = t74*t19
        som(32) = -192*t92*t4*t8
        som(33) = 256*t92*t31
        som(34) = 54*t92*t9
        t98 = t74*xpi(ip3)
        som(35) = -16*t98*t84
        som(36) = 8*t98*t89
        t101 = t74*t4
        som(37) = 16*t101*t84
        som(38) = -8*t101*t89
        t106 = dpipj(is3,ip1)
        som(39) = -16*t74*t5*t82*t12*t106
        som(40) = 8*t74*t82*t35*t106
        t112 = t41**2
        t113 = t112*t2
        som(41) = -40*t113*t13
        som(42) = -8*t113*t6
        som(43) = 40*t113*t35
        t117 = t112*t19
        som(44) = -96*t117*t4
        som(45) = 36*t117*t8
        som(46) = 9*t112*t1*t19
        som(47) = -8*t2*t6*t21
        som(48) = -64*t4*t28*t9
        som(49) = 32*t25*t14
 !	print '(7g20.12)',(som(i),i=1,49)
        n=49
      Else
 !prod3=
 !	som(1)=+160*xpi(ip1)*xpi(ip2)*xpi(is2)*piDpj(ip1,ip2)**2*
 !     +	dpipj(is2,is1)**2
 !	som(2)=-40*xpi(ip1)*xpi(ip2)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is2)*dpipj(is2,is1)**3
 !	som(3)=-88*xpi(ip1)*xpi(ip2)*piDpj(ip1,ip2)**2*dpipj(is2,
 !     +	is1)**3
 !	som(4)=+9*xpi(ip1)*xpi(ip2)**2*dpipj(is2,is1)**4
 !	som(5)=-128*xpi(ip1)*xpi(is2)*piDpj(ip1,ip2)**3*piDpj(ip2,
 !     +	is2)*dpipj(is2,is1)
 !	som(6)=-256*xpi(ip1)*xpi(is2)*piDpj(ip1,ip2)**4*dpipj(is2,is1)
 !	som(7)=+256*xpi(ip1)*xpi(is2)**2*piDpj(ip1,ip2)**4
 !	som(8)=-16*xpi(ip1)*piDpj(ip1,ip2)**2*piDpj(ip2,is2)**2*dpipj(
 !     +	is2,is1)**2
 !	som(9)=+64*xpi(ip1)*piDpj(ip1,ip2)**3*piDpj(ip2,is2)*dpipj(is2,
 !     +	is1)**2
 !	som(10)=+80*xpi(ip1)*piDpj(ip1,ip2)**4*dpipj(is2,is1)**2
 !	som(11)=+128*xpi(ip1)**2*xpi(ip2)*xpi(is2)*piDpj(ip1,ip2)*piDpj(
 !     +	ip2,is2)*dpipj(is2,is1)
 !	som(12)=+576*xpi(ip1)**2*xpi(ip2)*xpi(is2)*piDpj(ip1,ip2)**2*
 !     +	dpipj(is2,is1)
 !	som(13)=-512*xpi(ip1)**2*xpi(ip2)*xpi(is2)**2*piDpj(ip1,ip2)**2
 !	som(14)=-88*xpi(ip1)**2*xpi(ip2)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is2)*dpipj(is2,is1)**2
 !	som(15)=-192*xpi(ip1)**2*xpi(ip2)*piDpj(ip1,ip2)**2*dpipj(is2,
 !     +	is1)**2
 !	som(16)=+40*xpi(ip1)**2*xpi(ip2)*piDpj(ip2,is2)**2*dpipj(is2,
 !     +	is1)**2
 !	som(17)=-96*xpi(ip1)**2*xpi(ip2)**2*xpi(is2)*dpipj(is2,is1)**2
 !	som(18)=+60*xpi(ip1)**2*xpi(ip2)**2*dpipj(is2,is1)**3
 !	som(19)=+128*xpi(ip1)**2*xpi(is2)*piDpj(ip1,ip2)**2*piDpj(ip2,
 !     +	is2)**2
 !	som(20)=-128*xpi(ip1)**2*xpi(is2)*piDpj(ip1,ip2)**3*piDpj(ip2,
 !     +	is2)
 !	som(21)=-64*xpi(ip1)**2*xpi(is2)*piDpj(ip1,ip2)**4
 !	som(22)=-32*xpi(ip1)**2*piDpj(ip1,ip2)*piDpj(ip2,is2)**3*
 !     +	dpipj(is2,is1)
 !	som(23)=+64*xpi(ip1)**2*piDpj(ip1,ip2)**3*piDpj(ip2,is2)*
 !     +	dpipj(is2,is1)
 !	som(24)=+32*xpi(ip1)**2*piDpj(ip1,ip2)**4*dpipj(is2,is1)
 !	som(25)=+128*xpi(ip1)**3*xpi(ip2)*xpi(is2)*piDpj(ip1,ip2)*piDpj(
 !     +	ip2,is2)
 !	som(26)=+160*xpi(ip1)**3*xpi(ip2)*xpi(is2)*piDpj(ip1,ip2)**2
 !	som(27)=-128*xpi(ip1)**3*xpi(ip2)*xpi(is2)*piDpj(ip2,is2)**2
 !	som(28)=+32*xpi(ip1)**3*xpi(ip2)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is1)*piDpj(ip2,is2)
 !	som(29)=-88*xpi(ip1)**3*xpi(ip2)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is2)*dpipj(is2,is1)
 !	som(30)=-88*xpi(ip1)**3*xpi(ip2)*piDpj(ip1,ip2)**2*dpipj(is2,
 !     +	is1)
 !	som(31)=-16*xpi(ip1)**3*xpi(ip2)*piDpj(ip2,is1)*piDpj(ip2,
 !     +	is2)**2
 !	som(32)=+48*xpi(ip1)**3*xpi(ip2)*piDpj(ip2,is2)**2*dpipj(is2,
 !     +	is1)
 !	som(33)=-320*xpi(ip1)**3*xpi(ip2)**2*xpi(is2)*dpipj(is2,is1)
 !	som(34)=+256*xpi(ip1)**3*xpi(ip2)**2*xpi(is2)**2
 !	som(35)=+118*xpi(ip1)**3*xpi(ip2)**2*dpipj(is2,is1)**2
 !	som(36)=-16*xpi(ip1)**3*xpi(ip3)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is1)*piDpj(ip2,is2)
 !	som(37)=+8*xpi(ip1)**3*xpi(ip3)*piDpj(ip2,is1)*piDpj(ip2,is2)**2
 !	som(38)=+16*xpi(ip1)**3*xpi(is2)*piDpj(ip1,ip2)*piDpj(ip2,
 !     +	is1)*piDpj(ip2,is2)
 !	som(39)=-8*xpi(ip1)**3*xpi(is2)*piDpj(ip2,is1)*piDpj(ip2,is2)**2
 !	som(40)=-16*xpi(ip1)**3*piDpj(ip1,ip2)*piDpj(ip2,is1)*piDpj(ip2,
 !     +	is2)*dpipj(is3,ip1)
 !	som(41)=+8*xpi(ip1)**3*piDpj(ip2,is1)*piDpj(ip2,is2)**2*
 !     +	dpipj(is3,ip1)
 !	som(42)=-40*xpi(ip1)**4*xpi(ip2)*piDpj(ip1,ip2)*piDpj(ip2,is2)
 !	som(43)=-8*xpi(ip1)**4*xpi(ip2)*piDpj(ip1,ip2)**2
 !	som(44)=+40*xpi(ip1)**4*xpi(ip2)*piDpj(ip2,is2)**2
 !	som(45)=-96*xpi(ip1)**4*xpi(ip2)**2*xpi(is2)
 !	som(46)=+60*xpi(ip1)**4*xpi(ip2)**2*dpipj(is2,is1)
 !	som(47)=+9*xpi(ip1)**5*xpi(ip2)**2
 !	som(48)=-8*xpi(ip2)*piDpj(ip1,ip2)**2*dpipj(is2,is1)**4
 !	som(49)=-64*xpi(is2)*piDpj(ip1,ip2)**4*dpipj(is2,is1)**2
 !	som(50)=+32*piDpj(ip1,ip2)**3*piDpj(ip2,is2)*dpipj(is2,is1)**3
 !	som(51)=+32*piDpj(ip1,ip2)**4*dpipj(is2,is1)**3
 !	print '(7g20.12)',(som(i),i=1,51)
 !
 !	optimized by Maple (see ffxxyz.map)
 !
        t126 = xpi(ip1)
        t127 = xpi(ip2)
        t128 = t126*t127
        t129 = xpi(is2)
        t130 = piDpj(ip1,ip2)
        t131 = t130**2
        t132 = t129*t131
        t133 = dpipj(is2,is1)
        t134 = t133**2
        som(1) = 160*t128*t132*t134
        t137 = piDpj(ip2,is2)
        t138 = t130*t137
        t139 = t134*t133
        som(2) = -40*t128*t138*t139
        som(3) = -88*t128*t131*t139
        t144 = t127**2
        t146 = t134**2
        som(4) = 9*t126*t144*t146
        t148 = t126*t129
        t149 = t131*t130
        t150 = t149*t137
        som(5) = -128*t148*t150*t133
        t153 = t131**2
        som(6) = -256*t148*t153*t133
        t156 = t129**2
        som(7) = 256*t126*t156*t153
        t160 = t137**2
        t161 = t160*t134
        som(8) = -16*t126*t131*t161
        som(9) = 64*t126*t149*t137*t134
        som(10) = 80*t126*t153*t134
        t168 = t126**2
        t169 = t168*t127
        t171 = t138*t133
        som(11) = 128*t169*t129*t171
        som(12) = 576*t169*t132*t133
        som(13) = -512*t169*t156*t131
        som(14) = -88*t169*t138*t134
        som(15) = -192*t169*t131*t134
        som(16) = 40*t169*t161
        t182 = t168*t144
        som(17) = -96*t182*t129*t134
        som(18) = 60*t182*t139
        t186 = t168*t129
        som(19) = 128*t186*t131*t160
        som(20) = -128*t186*t150
        som(21) = -64*t186*t153
        som(22) = -32*t168*t130*t160*t137*t133
        som(23) = 64*t168*t149*t137*t133
        som(24) = 32*t168*t153*t133
        t200 = t168*t126
        t201 = t200*t127
        som(25) = 128*t201*t129*t130*t137
        som(26) = 160*t201*t132
        som(27) = -128*t201*t129*t160
        t208 = piDpj(ip2,is1)
        t210 = t130*t208*t137
        som(28) = 32*t201*t210
        som(29) = -88*t201*t171
        som(30) = -88*t201*t131*t133
        t215 = t208*t160
        som(31) = -16*t201*t215
        som(32) = 48*t201*t160*t133
        t219 = t200*t144
        som(33) = -320*t219*t129*t133
        som(34) = 256*t219*t156
        som(35) = 118*t219*t134
        t225 = t200*xpi(ip3)
        som(36) = -16*t225*t210
        som(37) = 8*t225*t215
        t228 = t200*t129
        som(38) = 16*t228*t210
        som(39) = -8*t228*t215
        t233 = dpipj(is3,ip1)
        som(40) = -16*t200*t130*t208*t137*t233
        som(41) = 8*t200*t208*t160*t233
        t239 = t168**2
        t240 = t239*t127
        som(42) = -40*t240*t138
        som(43) = -8*t240*t131
        som(44) = 40*t240*t160
        t244 = t239*t144
        som(45) = -96*t244*t129
        som(46) = 60*t244*t133
        som(47) = 9*t239*t126*t144
        som(48) = -8*t127*t131*t146
        som(49) = -64*t129*t153*t134
        som(50) = 32*t150*t139
        som(51) = 32*t153*t139
 !	print '(7g20.12)',(som(i),i=1,51)
        n=51
      End If
 !
      s = 0
      smax = 0
      Do j=1,n
        s = s + som(j)
        smax = Max(smax,som(j))
      End Do
      If ( iwarn < 3 ) Then
        hulp = 1/(16*xpi(ip1)**3*sdel2p**4*dy2z(3-iwarn)  &
           & * (y(1)-2*z(1))*(y(1)-2*z(2)))
      Else
        hulp = 1/(16*xpi(ip1)**3*sdel2p**4*dy2z(7-iwarn)  &
           & * (y(3)-2*z(3))*(y(3)-2*z(4)))
      End If
      s = s*hulp
      smax = smax*hulp
      If ( smax < xmax ) Then
        dy2z(iwarn) = s
        xmax = smax
      End If
    Else
      n=0
    End If
  End If
  ier = ier1
 !
Goto 200
 !  #] get dyz:
 !  #[ special case, get indices:
100 Continue
  If ( ivert==2 ) Then
    is1 = 2
    ip1 = 5
  Else
    is1 = 1
    ip1 = 6
  End If
 !  #] special case, get indices:
 !  #[ xk = 0:
  If ( xpi(ip1) == 0 )  Call WriteLfError(42)
 !  #] xk = 0:
 !  #[ get ypm,zpm:
 !
 !	special case del2s = 0, hence the roots are not the real roots
 !	but z_2'' = (z_2'-1)/delta, z''_3 = -z'_3/delta
 !
  hulp = sdel2s
  disc = delps/sdel2p
  If ( ivert == 3 ) Then
    hulp = -hulp
    disc = -disc
  End If
  If ( sdel2s == 0 ) Then
    isoort(1) = 102
    isoort(2) = 102
    z(1) = piDpj(is1,3)/xpi(3)
    z(2) = z(1)
  Else
    isoort(1) = 101
    isoort(2) = 101
    Call Roots(xpi(3),piDpj(is1,3),xpi(is1),hulp,z(1),z(2))
  End If
  Call Roots(xpi(3),piDpj(is1,3),etami(is1),disc,y(1),y(2))
 !  #] get ypm,zpm:
 !  #[ get ypm1,zpm1:
  z(3) = 1 - z(1)
  z(4) = 1 - z(2)
  If ( Abs(z(3))<xloss .Or. Abs(z(4))<xloss ) Then
    If ( ivert==2 ) Then
      Call Roots(xpi(3),piDpj(ip1,3),xpi(ip1),hulp, z(4),z(3))
    Else
      Call Roots(xpi(3),-piDpj(ip1,3),xpi(ip1),hulp ,z(4),z(3))
    End If
  End If
  y(3) = 1 - y(1)
  y(4) = 1 - y(2)
  If ( Abs(y(3)) < xloss .Or. Abs(y(4)) < xloss ) Then
    If ( ivert == 2 ) Then
      Call Roots(xpi(3),piDpj(ip1,3),etami(ip1), disc,y(4),y(3))
    Else
      Call Roots(xpi(3),-piDpj(ip1,3),etami(ip1), disc,y(4),y(3))
    End If
  End If
 !  #] get ypm1,zpm1:
 !  #[ get dypzp, dypzm:
  If ( isoort(1) == 2 ) Then
    dyz(2,1) = disc/xpi(3)
    dyz(2,2) = dyz(2,1)
  Else If ( disc > 0 .Eqv. sdel2s > 0 ) Then
    dyz(2,1) = ( disc + hulp )/xpi(3)
    dyz(2,2) = etalam/(xpi(3)*dyz(2,1))
  Else
    dyz(2,2) = ( disc - hulp )/xpi(3)
    dyz(2,1) = etalam/(xpi(3)*dyz(2,2))
  End If
  dyz(1,1) = -dyz(2,2)
  dyz(1,2) = -dyz(2,1)
  d2yzz = 2*disc/xpi(3)
 !
 !	these are very rarely needed, but ...
 !
  Do i=1,4
    j = 2*((i+1)/2)
    dy2z(i) = y(j) - 2*z(i)
    smax = Abs(y(j))
  End Do
 !  #] get dypzp, dypzm:
200 Continue
  End Subroutine ffxxyz

 Subroutine ffzli2(cli, clog, c_in)
 !------------------------------------------------------
 ! substitutes the original one, makes portation easier
 !------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: c_in
  Complex(dp), Intent(out) :: cli, clog
  cli = CLi2(c_in)
  clog = Log1minusX(c_in)
 End Subroutine ffzli2


 Subroutine ffzxdl(zxdilo,ipi12,zlog,x,ieps,ier)
 !---------------------------------------------------------------
 !	Computes the dilogarithm (Li2, Sp) for any (real) x
 !	to precision precx. If an error message is given add
 !	more bf's. For x > 1 the imaginary part is
 !	 -/+i*pi*log(x), corresponding to x+ieps.
 !	The number of factors pi^2/12 is passed separately in
 !	ipi12 for accuracy.  We also calculate log(1-x)
 !	which is likely to be needed.
 !
 !	Input:	x	(real)
 !		ieps	( Integer,+/-1)
 !
 !	Output: zxdilo	(complex) the dilog mod factors pi2/12
 !		ipi12	(Integer) these factors
 !		zlog	(complex) log(1-x)
 !
 !	Calls:	log,dfflo1
 !---------------------------------------------------------------
 Implicit None
  Integer :: ipi12,ieps,ier
  Real(dp) :: x
  Complex(dp) :: zxdilo,zlog

  Integer :: jsgn, nmax, i1
  Real(dp) :: fact,u,u2,a,xdilo
  Complex(dp) :: cy,cfact

 !  #] initialisations:
 !  #[ exceptional cases:
  If ( x == 1) Then
    zxdilo = 0
    zlog = -99999
    ipi12 = 2
    Return
  Else If (x == -1) Then
    zxdilo = 0
    zlog = xlg2
    ipi12 = -1
    Return
  Else If (x == 0.5_dp) Then
    zxdilo = - xlg2**2/2
    zlog = -xlg2
    ipi12 = 1
    Return
  Else If ( Abs(x) < precx ) Then
    zxdilo = x
    zlog = -x
    ipi12 = 0
    Return
  End If
 !  #] exceptional cases:
 !  #[ transform to (-1,.5):
  If (x < -1) Then
    fact = Log(-x)
    cy = - fact**2/2
    ipi12 = -2
    If ( -x*xloss > 1 ) Then
      u = -Log1MinusX(1/x)
    Else
      u = -Log(1-1/x)
    End If
    zlog = Log(1-x)
    jsgn = -1
  Else If ( x < 0.5_dp) Then
    cy = 0
    ipi12 = 0
    If ( Abs(x) < xloss ) Then
      zlog = Log1MinusX(x)
    Else
      zlog = Log(1-x)
    End If
    u = -Real(zlog,dp)
    jsgn = 1
  Else If ( x <= 2 ) Then
    u = -Log(x)
    If ( Abs(1-x) < xalogm ) Then
      cy = 0
    Else If ( x < 1 ) Then
      zlog = Log(1-x)
      cy = Real(u,dp)*zlog
    Else If ( ieps > 0 ) Then
      zlog = Cmplx(Log(x-1),-pi,dp)
      cy = Real(u,dp)*zlog
    Else
      zlog = Cmplx(Log(x-1),+pi,dp)
      cy = Real(u,dp)*zlog
    End If
    ipi12 = 2
    jsgn = -1
  Else
    If ( ieps > 0 ) Then
      cfact = Cmplx(Log(x),-pi,dp)
      zlog = Cmplx(Log(x-1),-pi,dp)
    Else
      cfact = Cmplx(Log(x),+pi,dp)
      zlog = Cmplx(Log(x-1),+pi,dp)
    End If
    cy = - cfact**2/2
    ipi12 = -2
    If ( x*xloss > 1 ) Then
      u = -Log1MinusX(1/x)
    Else
      u = -Log(1-1/x)
    End If
    jsgn = -1
  End If
 !  #] transform to (-1,.5):
 !  #[ calculate dilog:
  If ( Abs(u) < xalog2 ) Then
    xdilo = u
  Else
    u2 = u**2
    a = Abs(u2)
    nmax = Findbound(a, 1, bf)
    if (nmax.gt.2) then
     xdilo = u2 * bf(nmax)
     Do i1=nmax-1,3,-1
      xdilo = u2 * (bf(i1) + xdilo)
     End Do
    else
     xdilo = 0._dp
    end if
 !	watch the powers of u.
    xdilo = u + u2*(bf(1) + u*(bf(2) + xdilo))
  End If
  If(jsgn==1)Then
    zxdilo =  Real(xdilo,dp) + cy
  Else
    zxdilo = -Real(xdilo,dp) + cy
  End If
 !  #] calculate dilog:
 End Subroutine ffzxdl

 Subroutine ffzzdl(zdilog,ipi12,zlog,cx,ier)
 !-------------------------------------------------------------
 !	Computes the dilogarithm (Li2, Sp) for any (complex) cx
 !	to about 15 sign  Ificant figures. This can be improved
 !	by adding more of the bf's. For real cx > 1 an error is
 !	generated as the imaginary part is undefined then.
 !	For use in ffcdbd zlog = log(1-cx) is also calculated
 !
 !	Input:	cx	(complex)
 !
 !	Output:	zdilog	(complex) Li2(cx) mod factors pi^2/12
 !		ipi12	(Integer) these factors
 !		zlog	(complex) log(1-cx)
 !
 !	Calls:	log,zfflo1,(d/a)imag,real/Real
 !-------------------------------------------------------------
  Implicit None
   Integer :: ipi12,ier
   Complex(dp) :: zdilog,zlog,cx

   Integer :: jsgn, i1, nmax
   Real(dp) :: xi,xr,s1,s2,xa,a
   Complex(dp) :: cfact,cx1,cy,cz,cz2

    !  #[ exceptional cases:
    xi  = Aimag(cx)
    xr  = Real(cx,dp)
    If ( xi == 0 ) Then
      If ( xr > 1 )  Call WriteLfError(41)
      Call ffzxdl(zdilog,ipi12,zlog,xr,1,ier)
      Return
    End If
    If ( Abs(xi) < xalog2 ) Then
      s1 = 0
    Else
      s1 = xi**2
    End If
    If ( Abs(xr) < xalog2 ) Then
      s2 = 0
    Else
      s2 = xr**2
    End If
    xa = Sqrt(s1 + s2)
    If ( xa < precc ) Then
      zdilog = cx
      zlog = -cx
      ipi12 = 0
      Return
    End If
    !  #] exceptional cases: 
    !  #[ transform to |x|<1, Re(x) < 0.5:
    If ( xr <= 0.5_dp) Then
      If (xa > 1) Then
        If ( 1/xa < xalogm ) Then
          cfact = 0
        Else If ( 1/xa < xclogm ) Then
          cx1 = cx*Real(1/xa,dp)
          cfact = Log(-cx1) + Log(Real(xa,dp))
        Else
          cfact = Log(-cx)
        End If
        cy = - cfact**2/2
        ipi12 = -2
        If ( xa*xloss**2 > 1) Then
          If ( 1/xa < xclogm ) Then
            cx1 = cx*Real(1/xa,dp)
            cx1 = 1/cx1
            cx1 = cx1*Real(1/xa,dp)
          Else
            cx1 = 1/cx
          End If
          cz = -Log1MinusX(cx1)
        Else
          cz = -Log(1-1/cx)
        End If
        zlog = Log(1-cx)
        jsgn = -1
      Else
        cy = 0
        ipi12 = 0
        If ( xa < xloss**2 ) Then
          zlog = Log1MinusX(cx)
        Else
          zlog = Log(1-cx)
        End If
        cz = -zlog
        jsgn = 1
      End If
    Else
      If (xa <= Sqrt(2*xr)) Then
        cz = -Log(cx)
        If ( Abs(xr-1) + Abs(xi) < xclogm ) Then
          cy = 0
        Else
          zlog = Log(1-cx)
          cy = cz*zlog
        End If
        ipi12 = 2
        jsgn = -1
      Else
        If ( 1/xa < xalogm ) Then
          cfact = 0
        Else If ( 1/xa < xclogm ) Then
          cx1 = cx*Real(1/xa,dp)
          cfact = Log(-cx1) + Log(Real(xa,dp))
        Else
          cfact = Log(-cx)
        End If
        cy = - cfact**2/2
        ipi12 = -2
        If ( xa*xloss > 1) Then
          If ( 1/xa < xclogm ) Then
            cx1 = cx*Real(1/xa,dp)
            cx1 = 1/cx1
            cx1 = cx1*Real(1/xa,dp)
          Else
            cx1 = 1/cx
          End If
          cz = -Log1MinusX(cx1)
        Else
          cz = -Log(1-1/cx)
        End If
        zlog = Log(1-cx)
        jsgn = -1
      End If
    End If
    !  #] transform to |x|<1, Re(x) < 0.5: 
    !  #[ get dilog:
    If ( absc(cz) < xclogm ) Then
      zdilog = cz
    Else
      cz2 = cz*cz
      a = Real(cz,dp)**2 + Aimag(cz)**2
      nmax = Findbound(a, 1, bf)
      if (nmax.gt.2) then
       zdilog = cz2 * bf(nmax)
       Do i1=nmax-1,3,-1
        zdilog = cz2 * (bf(i1) + zdilog)
       End Do
      else
       zdilog = 0._dp
      end if
    !	watch the powers of z.
      zdilog = cz + cz2*(bf(1) + cz*(bf(2) + zdilog))
    End If
    If(jsgn==1)Then
      zdilog =  zdilog + cy
    Else
      zdilog = -zdilog + cy
    End If
    !  #] get dilog: 
  End Subroutine ffzzdl



 Integer Function FindBound(x, i1, array)
 !-----------------------------------------------------------------------
 ! finds the maximal power of a series which is numerical significant
 ! written by Werner Porod, 03.04.02
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i1
  Real(dp), Intent(in) :: array(:), x

  Integer :: nmax, i_u, i_l
  Real(dp) :: test1, test2, test3

  nmax = Size(array)

  If (i1.Eq.0) Then
   test1 = array(1)
  Else If (i1.Eq.1) Then
   test1 = array(1) * x
  Else
   test1 = array(1) * x**i1
  End If

  test2 = test1 + array(nmax) * x**(nmax+i1)
  If (test1.Ne.test2) Then
   FindBound = nmax
  Else
   i_l = 2
   i_u = nmax - 1
   Do
    test3 = test1 + array(i_u+1) * x**(i_u+1)
    test2 = test1 + array(i_u) * x**(i_u)
    If ( (test1.Ne.test2).And.(test1.Eq.test3)) Then
     Exit
    Else If (test1.Eq.test2) Then
     nmax = i_u
     i_u = (i_l + i_u)/2
    Else
     i_l = i_u
     i_u = (i_u + nmax) / 2
    End If
    If (i_u.Le.i_l) Exit
   End Do
   FindBound = i_u
  End If
 End Function FindBound


 Complex(dp) Function Floop(p2,m12,m22)
 !-----------------------------------------------------------------------
 ! calculates the function F as defined in J. Bagger at al, Nucl.Phys.B
 ! written by Werner Porod, 4.8.1999
 ! 18.05.2001: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: p2,m12,m22

  If ((m12.Eq.0._dp).and.(m22.eq.0._dp)) then
   Floop =  - 2._dp * p2 * B0(p2,m12,m22)
  
  Else If (m12.eq.0._dp) then
   Floop = - 2._dp * A0(m22) - (2._dp * p2 - m22) * B0(p2,m12,m22)

  Else If (m22.eq.0._dp) then
   Floop = A0(m12) - 2._dp * (p2 + m12) * B0(p2,m12,m22)

  Else If (m12.eq.m22) then
   Floop = - A0(m12) - (2._dp * p2 + m12) * B0(p2,m12,m12)

  Else
   Floop = A0(m12) - 2._dp * A0(m22) &
       & -  (2._dp * p2 + 2._dp * m12 - m22) * B0(p2,m12,m22)

  End If

 End Function Floop

 Subroutine Get_All_Ci(p1, p2, p1p2, m1, m2, m3, C_0, C1, C2, C00, C11, C12 &
    & , C22, C001, C002, C111, C112, C122, C222)
 Implicit None
  Real(dp), Intent(in) :: p1, p2, p1p2, m1, m2, m3
  Complex(dp), Intent(out) :: C_0, C1, C2, C00, C11, C12, C22, C001, C002 &
    & , C111, C112, C122, C222

  Real(dp) :: M11, M12, M22, det2, f1, f2
  Complex(dp) :: b1123, b023, b123, b1113, b113, b1112, b112, s1, s2
  Integer :: i1

  If (l_look_up_cache) Then
   Do i1=1,num_ci
    If (m1.Eq.C_cache(i1)%mi(1)) Then
     If (m2.Eq.C_cache(i1)%mi(2)) Then
      If (m3.Eq.C_cache(i1)%mi(3)) Then
       If (p1.Eq.C_cache(i1)%pi(1)) Then
        If (p2.Eq.C_cache(i1)%pi(2)) Then
         If (p1p2.Eq.C_cache(i1)%pi(3)) Then
          C_0 = C_cache(num_ci)%Ci(1)
          C1 = C_cache(num_ci)%Ci(2)
          C2 = C_cache(num_ci)%Ci(3)
          C00 = C_cache(num_ci)%Ci(4)
          C11 = C_cache(num_ci)%Ci(5)
          C12 = C_cache(num_ci)%Ci(6)
          C22 = C_cache(num_ci)%Ci(7)
          C001 = C_cache(num_ci)%Ci(8)
          C002 = C_cache(num_ci)%Ci(9)
          C111 = C_cache(num_ci)%Ci(10)
          C112 = C_cache(num_ci)%Ci(11)
          C122 = C_cache(num_ci)%Ci(12)
          C222 = C_cache(num_ci)%Ci(13)
          Return
         End If
        End If
       End If
      End If
     End If
    End If
   End Do
   If (num_ci .Lt. n_c_max) num_ci = num_ci + 1
  End If

  M12 = 0.5_dp * (p2 - p1p2 - p1)
  det2 = 2._dp * (p1*p1p2 - M12*M12)
  M12 = M12 / det2
  M11 = p1p2/det2
  M22 = p1/det2
  f1 = m1 - m2 + p1
  f2 = m1 - m3 + p1p2
  C_0 = C0(p1, p2, p1p2, m1, m2, m3)

  b1123 = B11(p2, m3, m2)
  b023 = cb0
  b123 = cb1
  b1113 = B11(p1p2, m1, m3)
  b113 = cb1
  s1 = cb0 - b023 - f1 * C_0

  b1112 = B11(p1, m1, m2)
  b112 = cb1
  s2 = cb0 - b023 - f2*C_0
  C1 = lc1()
  C2 = lc2()
  C00 = 0.5_DP*(m1*C_0 +0.5_DP*(f1*C1 + f2*C2 + b023)) + 0.25_DP

  s1 = -(f1*C1 + b123) - 2*C00
  s2 = -(f2*C1 + b123 - b112)
  C11 = lc1()
  C12 = lc2()

  b023 = b023 + b123
  s1 = b023 + b113 - f1*C2
  s2 = b023 - f2*C2 - 2*C00
  C22 = lc2()
  C001 = (m1*C1 + 0.5_DP*(f1*C11 +f2*C12 + b123))/3._DP - 1._dp/18._dp
  C002 = (m1*C2 + 0.5_DP*(f1*C12 +f2*C22 - b023))/3._DP - 1._dp/18._dp

  s1 = -(b1123 + f1*C11) - 4*C001
  s2 = -(b1123 + f2*C11 - b1112)
  C111 = lc1()
  C112 = lc2()

  b1123 = b1123 + b023 + b123
  s1 = -(b1123 + f1*C22 - b1113)
  s2 = -(b1123 + f2*C22) - 4*C002
  C122 = lc1()
  C222 = lc2()
  
  If (l_look_up_cache) Then  ! storing data in cache
   C_cache(num_ci)%mi(1) = m1
   C_cache(num_ci)%mi(2) = m2
   C_cache(num_ci)%mi(3) = m3
   C_cache(num_ci)%pi(1) = p1
   C_cache(num_ci)%pi(2) = p2
   C_cache(num_ci)%pi(3) = p1p2
   C_cache(num_ci)%Ci(1) = C_0
   C_cache(num_ci)%Ci(2) = C1
   C_cache(num_ci)%Ci(3) = C2
   C_cache(num_ci)%Ci(4) = C00
   C_cache(num_ci)%Ci(5) = C11
   C_cache(num_ci)%Ci(6) = C12
   C_cache(num_ci)%Ci(7) = C22
   C_cache(num_ci)%Ci(8) = C001
   C_cache(num_ci)%Ci(9) = C002
   C_cache(num_ci)%Ci(10) = C111
   C_cache(num_ci)%Ci(11) = C112
   C_cache(num_ci)%Ci(12) = C122
   C_cache(num_ci)%Ci(13) = C222
  End If

 Contains
  Complex(dp) Function lc1()
   lc1 = M11 * s1 + M12 * s2
  End Function lc1

  Complex(dp) Function lc2()
   lc2 = M12 * s1 + M22 * s2
  End Function lc2

 End Subroutine Get_All_Ci


 Real(dp) Function GetRenormalizationScale()
 !-----------------------------------------------------------------------
 ! gets renormalization scale, gives old scale as output for checking
 ! written by Werner Porod
 ! 15.11.01
 !-----------------------------------------------------------------------
 Implicit None

  GetRenormalizationScale = mudim2

 End Function GetRenormalizationScale


 Complex(dp) Function Gloop(p2,m12,m22)
 !-----------------------------------------------------------------------
 ! calculates the function G as defined in J. Bagger at al, Nucl.Phys.B
 ! written by Werner Porod, 4.8.1999
 ! 18.05.2001: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: p2,m12,m22

  If ((m12.eq.0._dp).and.(m22.eq.0._dp)) then
   Gloop = p2 * B0(p2,m12,m22)

  Else If (m12.eq.0._dp) then
   Gloop = - A0(m22) + (p2 - m22) * B0(p2,m12,m22)

  Else If (m22.eq.0._dp) then
   Gloop = - A0(m12) + (p2 - m12) * B0(p2,m12,m22)

  Else If (m12.eq.m22) then
   Gloop =  - 2._dp * A0(m12) + (p2 - m12 - m22) * B0(p2,m12,m22)

  Else 
   Gloop =  - A0(m12) - A0(m22) + (p2 - m12 - m22) * B0(p2,m12,m22)

  End If

 End Function Gloop


 Complex(dp) Function Hloop(p2,m12,m22)
 !-----------------------------------------------------------------------
 ! calculates the function H as defined in J. Bagger at al, Nucl.Phys.B
 ! written by Werner Porod, 4.8.1999
 ! 18.05.2001: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: p2,m12,m22

  Hloop = 4._dp * B00(p2,m12,m22) + GLoop(p2,m12,m22) 

 End Function Hloop

 Real(dp) Function h_0_1(x)
 !--------------------------------------------------------------
 ! loop function, PhD thesis, Ch. Bobeth, C.8
 !--------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x
!  Integer :: i1, sum

  If (x.Eq.1._dp) Then
   h_0_1 = 73._dp / 72._dp

  Else
   h_0_1 = ( 47._dp * x**3 - 237._dp * x**2 + 312._dp * x - 104._dp) / 108._dp &
         & - (3._dp * x**4 - 30._dp * x**3 + 54._dp * x**2 -32._dp * x + 8._dp) &
         &   * Log(x) / (18._dp * (x-1._dp) )
   h_0_1 = h_0_1 / (x-1._dp)**3
  End If

 End Function h_0_1


 Real(dp) Function h_0_5(x)
 !--------------------------------------------------------------
 ! loop function, Goto et al., PRF55 (1997) 4273, Gl. E5
 !--------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Logical, Save :: l_first=.True.
  Integer, Parameter :: n_max=16
  Real(dp), Save :: ci(0:n_max) 
  Integer :: i1, sum
  Real(dp) :: r

  If (l_first) Then
   l_first = .False.
   sum = 3
   Do i1=0,n_max
    sum = sum + 3 + i1
    ci(i1) = (-1)**i1 * (9._dp + 2._dp * i1) / (12._dp * sum)
   End Do
  End If

  If (x.Eq.1._dp) Then
   h_0_5 = 1._dp / 8._dp

  Else If (Abs(x-1._dp).Le.1.e-2_dp) Then
   r = x - 1._dp
   h_0_5 = ci(n_max)
   Do i1=n_max-1,0,-1
    h_0_5 = ci(i1) + h_0_5 * r
   End Do

  Else
   h_0_5 = ( -16._dp + 45._dp*x - 36._dp*x**2 + 7._dp*x**3  &
         & + (-12._dp + 18._dp*x) *Log(x) ) / (36._dp*(x-1._dp)**4)
  End If

 End Function h_0_5
 

 Real(dp) Function h_0_6(x)
 !--------------------------------------------------------------
 ! loop function, Goto et al., PRF55 (1997) 4273, Gl. E6
 !--------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Logical, Save :: l_first=.True.
  Integer, Parameter :: n_max=16
  Real(dp), Save :: ci(0:n_max) 
  Integer :: i1, sum
  Real(dp) :: r

  If (l_first) Then
   l_first = .False.
   sum = 24
   ci(0) = 1._dp / Real(sum,dp)
   Do i1=1,n_max
    sum = sum * (4+i1) / i1
    ci(i1) = (-1)**i1 / Real(sum,dp)
   End Do
  End If

  If (x.Eq.1._dp) Then
   h_0_6 = 1._dp/24._dp

  Else If (x.Eq.0._dp) Then
   h_0_6 = 1._dp/18._dp

  Else If (Abs(x-1._dp).Le.1.e-2_dp) Then
   r = x - 1._dp
   h_0_6 = ci(n_max)
   Do i1=n_max-1,0,-1
    h_0_6 = ci(i1) + h_0_6 * r
   End Do

  Else
   h_0_6 = ( 2._dp - 9._dp*x + 18._dp *x**2 - 11._dp*x**3 &
         & + 6._dp*x**3*Log(x) )/ (36._dp*(-1._dp + x)**4)
  End If

 End Function h_0_6
 

 Complex(dp) Function Igamma(mj2,mi2,mb2,mf2)
 !-----------------------------------------------------------------------
 ! calculates the funtion I(mj2,mi2,mb2,mf2) as given in Haber and Wyler
 ! written by Werner Porod, 25.12.99
 ! 18.05.2001: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: mj2,mi2,mb2,mf2

  Igamma = C0(mj2,mi2,0._dp,mf2,mb2,mf2)

 End Function Igamma


 Complex(dp) Function I2gamma(mj2,mi2,mb2,mf2)
 !-----------------------------------------------------------------------
 ! calculates the funtion I(mj2,mi2,mb2,mf2) as given in Haber and Wyler
 ! written by Werner Porod, 25.12.99
 ! 18.05.2001: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: mj2,mi2,mb2,mf2

  I2gamma = - Cget("C1  ",mj2,mi2,0._dp,mf2,mb2,mf2)
!  I2gamma = - Cget("C1  ",mj2,mi2,0._dp,mf2,mf2,mb2)
!  I2gamma = - Cget("C2  ",mj2,mi2,0._dp,mb2,mf2,mf2)
!  I2gamma = - Cget("C2  ",mj2,mi2,0._dp,mf2,mb2,mf2)
!  I2gamma = - Cget("C2  ",mj2,mi2,0._dp,mf2,mf2,mb2)

 End Function I2gamma


 Complex(dp) Function Jgamma(mj2,mi2,mb2,mf2)
 !-----------------------------------------------------------------------
 ! calculates the funtion I(mj2,mi2,mb2,mf2) as given in Haber and Wyler
 ! written by Werner Porod, 25.12.99
 ! 18.05.2001: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: mj2,mi2,mb2,mf2

  Jgamma = C0(mj2,mi2,0._dp,mb2,mf2,mb2)

 End Function Jgamma


  Complex(dp) Function Kgamma(mj2,mi2,mb2,mf2)
 !-----------------------------------------------------------------------
 ! calculates the funtion I(mj2,mi2,mb2,mf2) as given in Haber and Wyler
 ! written by Werner Porod, 25.12.99
 ! 18.05.2001: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: mj2,mi2,mb2,mf2

  Kgamma = ( 1._dp                                        &
         & + mf2 * C0(mj2,mi2,0._dp,mf2,mb2,mf2)           &
         & + mj2 * Cget("C1  ",mj2,mi2,0._dp,mf2,mb2,mf2)    &
         & + mb2 * C0(mj2,mi2,0._dp,mb2,mf2,mb2) ) / (mi2 - mj2)

 End Function Kgamma

!! Check if Boxes are called with one massless particle

Complex(dp) Function D00check(a,b,c,d)
Implicit None
Real(dp), Intent(in) :: a,b,c,d
D00check = KilianD00(a,b,c,d)
If ((Real(D00check,dp).ne.Real(D00check,dp)).or. &
  & (Abs(D00check).gt.1.0E+30_dp)) Then 

Write(ErrCan,*) "Numerical problem in D00check"
Write(ErrCan,*) "Involved masses: ",a,b,c,d

 If (ErrorLevel.gt.0) Then
   Call TerminateProgram
 Else 
   D00check = 0._dp
  End if
End if
End Function D00check

Complex(dp) Function C0check(a,b,c)
Implicit None
Real(dp), Intent(in) :: a,b,c
C0check = C0_3m(a,b,c)
If ((Real(C0check,dp).ne.Real(C0check,dp)).or. &
  & (Abs(C0check).gt.1.0E+30_dp))  Then 
    Write(ErrCan,*) "Numerical problem in C0check"
    Write(ErrCan,*) "Involved masses: ",a,b,c

     If (ErrorLevel.gt.0) Then
       Call TerminateProgram
     Else 
       C0check = 0._dp
     End if 
End if
End Function C0check

Complex(dp) Function D0check(a,b,c,d)
Implicit None
Real(dp), Intent(in) :: a,b,c,d
D0check = KilianD0(a,b,c,d)
If ((Real(D0check,dp).ne.Real(D0check,dp)).or.  & 
 & (Abs(D0check).gt.1.0E+30_dp))  Then 
   Write(ErrCan,*) "Numerical problem in D0check"
   Write(ErrCan,*) "Involved masses: ",a,b,c,d

    If (ErrorLevel.gt.0) Then
       Call TerminateProgram
    Else 
       D0check = 0._dp
    End if
End if
End Function D0check

Complex(dp) Function C0D0check(a,b,c,d)
Implicit None
Real(dp), Intent(in) :: a,b,c,d
C0D0check = KilianC0D0(a,b,c,d)
If ((Real(C0D0check,dp).ne.Real(C0D0check,dp)).or. &
  & (Abs(C0D0check).gt.1.0E+30_dp))  Then 
   Write(ErrCan,*) "Numerical problem in C0D0check"
   Write(ErrCan,*) "Involved masses: ",a,b,c,d

   If (ErrorLevel.gt.0) Then
     Call TerminateProgram
   Else 
     C0D0check = 0._dp
   End if 
End if
End Function C0D0check

Complex(dp) Function B0C0check(a,b,c)
Implicit None
Real(dp), Intent(in) :: a,b,c
B0C0check = KilianB0C0(a,b,c)
If ((Real(B0C0check,dp).ne.Real(B0C0check,dp)).or. &
 & (Abs(B0C0check).gt.1.0E+30_dp))  Then 

   Write(ErrCan,*) "Numerical problem in B0C0check"
   Write(ErrCan,*) "Involved masses: ",a,b,c

   If (ErrorLevel.gt.0) Then
     Call TerminateProgram
   Else 
     B0C0check = 0._dp
   End if
End if
End Function B0C0check

Real(dp) Function KilianD0(a,b,c,d)
Real(dp),Intent(in) :: a,b,c,d
Real(dp) :: eps = 1E-12_dp, max_mass, ra, rb, rc, rd

If ((a.ge.b).and.(a.ge.c).and.(a.ge.d)) Then
  max_mass = a
Else  If ((b.ge.c).and.(b.ge.d)) Then
  max_mass = b  
Else  If (c.ge.d) Then
  max_mass = c
Else
  max_mass = d
End if

ra = Abs(a/max_mass)
rb = Abs(b/max_mass)
rc = Abs(c/max_mass)
rd = Abs(d/max_mass)

If ((ra.gt.eps).and.(rb.gt.eps).and.(rc.gt.eps).and.(rd.gt.eps)) Then
 KilianD0 = D0_Bagger(a,b,c,d)
Else If ((ra.lt.eps).and.(rb.gt.eps).and.(rc.gt.eps).and.(rd.gt.eps)) Then
 KilianD0 = KilianD0_3m(b,c,d)
Else If ((rb.lt.eps).and.(ra.gt.eps).and.(rc.gt.eps).and.(rd.gt.eps)) Then
 KilianD0 = KilianD0_3m(a,c,d)
Else If ((rc.lt.eps).and.(rb.gt.eps).and.(ra.gt.eps).and.(rd.gt.eps)) Then
 KilianD0 = KilianD0_3m(a,b,d)
Else If ((rd.lt.eps).and.(rb.gt.eps).and.(rc.gt.eps).and.(ra.gt.eps)) Then
 KilianD0 = KilianD0_3m(a,b,c)
Else
   Write(*,*) "Numerical problem in KilianD0"
   Write(*,*) "Involved masses: ",a,b,c,d

  If (ErrorLevel.gt.0) Then 
    Call TerminateProgram 
  Else 
    KilianD0 = 0._dp
  End if
End If
End Function KilianD0

Real(dp) Function MMD0(m1,m2,a,b,c,d)
Real(dp),Intent(in) :: m1,m2,a,b,c,d
Real(dp) :: eps = 1E-12_dp, max_mass, r1, r2
! calculates m1*m2*D0(a,b,c,d)

If ((a.ge.b).and.(a.ge.c).and.(a.ge.d)) Then
  max_mass = a
Else  If ((b.ge.c).and.(b.ge.d)) Then
  max_mass = b  
Else  If (c.ge.d) Then
  max_mass = c
Else
  max_mass = d
End if

r1 = Abs(m1**2/max_mass)
r2 = Abs(m2**2/max_mass)


If ((r1.le.eps).or.(r2.le.eps)) Then
 MMD0 = 0._dp
Else if ((c.le.1E-24_dp).and.(d.le.1E-24_dp)) Then  ! Two-Photon boxes
 MMD0 = 0._dp
Else 
 MMD0=m1*m2*D0check(a,b,c,d)
End if

! Write(*,*) r1, r2, MMD0

End Function MMD0


Real(dp) Function KilianD00(a,b,c,d)
Real(dp),Intent(in) :: a,b,c,d
Real(dp) :: eps = 1E-12_dp, max_mass, ra, rb, rc, rd

If ((a.ge.b).and.(a.ge.c).and.(a.ge.d)) Then
  max_mass = a
Else  If ((b.ge.c).and.(b.ge.d)) Then
  max_mass = b  
Else  If (c.ge.d) Then
  max_mass = c
Else
  max_mass = d
End if

ra = Abs(a/max_mass)
rb = Abs(b/max_mass)
rc = Abs(c/max_mass)
rd = Abs(d/max_mass)

If ((ra.le.eps).and.(rb.le.eps)) Then
 If (Abs(c-d).le.eps) Then 
   KilianD00 = -1._dp/(4._dp*d)
 Else
   KilianD00 = Log(d/c)/(4._dp*(c-d))
 End if
Else If ((ra.le.eps).and.(rc.le.eps)) Then
 If (Abs(b-d).le.eps) Then 
   KilianD00 = -1._dp/(4._dp*d)
 Else
   KilianD00 = Log(d/b)/(4._dp*(b-d))
 End if
Else If ((ra.le.eps).and.(rd.le.eps)) Then
 If (Abs(b-c).le.eps) Then 
   KilianD00 = -1._dp/(4._dp*c)
 Else
   KilianD00 = Log(c/b)/(4._dp*(b-c))
 End if
Else If ((rb.le.eps).and.(rc.le.eps)) Then
 If (Abs(a-d).le.eps) Then 
   KilianD00 = -1._dp/(4._dp*d)
 Else
   KilianD00 = Log(d/a)/(4._dp*(a-d))
 End if
Else If ((rb.le.eps).and.(rd.le.eps)) Then
 If (Abs(a-c).le.eps) Then 
   KilianD00 = -1._dp/(4._dp*c)
 Else
   KilianD00 = Log(c/a)/(4._dp*(a-c))
 End if
Else If ((rc.le.eps).and.(rd.le.eps)) Then
 If (Abs(a-b).le.eps) Then 
   KilianD00 = -1._dp/(4._dp*b)
 Else
   KilianD00 = Log(b/a)/(4._dp*(a-b))
 End if
Else If ((rd.le.eps).and.(ra.gt.eps).and.(rb.gt.eps).and.(rc.gt.eps)) Then
 KilianD00 = KilianD00_3m(a,b,c)
Else If ((ra.le.eps).and.(rd.gt.eps).and.(rb.gt.eps).and.(rc.gt.eps)) Then
 KilianD00 = KilianD00_3m(b,c,d)
Else If ((rb.le.eps).and.(rd.gt.eps).and.(ra.gt.eps).and.(rc.gt.eps)) Then
 KilianD00 = KilianD00_3m(a,c,d)
Else If ((rc.le.eps).and.(rd.gt.eps).and.(ra.gt.eps).and.(rb.gt.eps)) Then
 KilianD00 = KilianD00_3m(a,b,d)
Else
 KilianD00 = D27_Bagger(a,b,c,d)
End If


End Function KilianD00


Real(dp) Function KilianD0_3m(a,b,c)
Real (dp), Intent(in) :: a,b,c
! calculates D0(0,a,b,c)
! all unequal 0
If ((a.ne.0._dp).and.(b.ne.0._dp).and.(c.ne.0._dp)) Then
 If ((a.eq.b)) Then
  KilianD0_3m = (-a+c-a*log(c/a))/(a*(a-c)**2)
 Else If (a.eq.c) Then 
  KilianD0_3m = (-a+b-a*log(b/a))/(a*(a-b)**2)
 Else If (b.eq.c) Then
  KilianD0_3m = (-b+a-b*log(a/b))/(b*(b-a)**2)
 Else
  KilianD0_3m = D0_Bagger(0._dp,a,b,c)
 End If
End If
End Function KilianD0_3m


Real(dp) Function KilianD00_3m(a,b,c)
Real (dp), Intent(in) :: a,b,c
! calculates D00(0,a,b,c)

! all unequal 0
If ((a.ne.0._dp).and.(b.ne.0._dp).and.(c.ne.0._dp)) Then
 If ((b.eq.c)) Then
  KilianD00_3m = (a-b+a*log(b/a))/(4._dp*(a-b)**2)
  Else If (a.eq.b) Then 
  KilianD00_3m = (c-b+c*log(b/c))/(4._dp*(c-b)**2)
  Else If (a.eq.c) Then
  KilianD00_3m = (b-a+b*log(a/b))/(4._dp*(b-a)**2)
  Else
  KilianD00_3m = D27_Bagger(0._dp,a,b,c)
 End If
Else If (a.eq.0._dp) Then
 If (b.ne.c) Then
  KilianD00_3m = log(c/b)/(4._dp*(b-c))
 Else
  KilianD00_3m = -1._dp/(8._dp*b)
 End If
Else If (b.eq.0._dp) Then
 If (a.ne.c) Then
  KilianD00_3m = log(c/a)/(4._dp*(a-c))
 Else
  KilianD00_3m = -1._dp/(8._dp*a)
 End If
Else If (c.eq.0._dp) Then
  If (a.ne.b) Then
  KilianD00_3m = log(b/a)/(4._dp*(a-b))
 Else
  KilianD00_3m = -1._dp/(8._dp*a)
 End If
End If
End Function KilianD00_3m

!! Analytical limit of sum of diagrams

! Real(dp) Function KilianC0D0(a,b,c,d)
!  ! This function is C0(a,b,c)+d*D0(a,b,c,d)
!  ! it catches the case a=b=0
!  ! cutoff for masses = sqrt(10^-10) GeV
!  Implicit None
!  Real(dp), Intent(in) :: a,b,c,d
!  Real(dp) :: D0user, C0user
! 
!  If ((a.le.1.E-10_dp).and.(b.le.1.E-10_dp)) Then
!   KilianC0D0 = 1./(c-d)*Log(d/c)
! Else If ((a.le.1.E-10_dp).and.(c.le.1.E-10_dp)) Then
!  KilianC0D0 = 1./(b-d)*Log(d/b)
! Else If ((c.le.1.E-10_dp).and.(b.le.1.E-10_dp)) Then
!  KilianC0D0 = 1./(a-d)*Log(d/a)
! Else
!   C0user = C0_3m(a,b,c)
!   D0user = D0_Bagger(a,b,c,d)
!   KilianC0D0 = C0user+d*D0user
!   End If
! End Function 

Real(dp) Function KilianC0D0(a,b,c,d)
 ! This function is C0(a,b,c)+d*D0(a,b,c,d)
 Implicit None
 Real(dp), Intent(in) :: a,b,c,d
 Real(dp) :: D0user, C0user
 Real(dp) :: eps = 1E-12_dp, max_mass, ra, rb, rc, rd

If ((a.ge.b).and.(a.ge.c).and.(a.ge.d)) Then
  max_mass = a
Else  If ((b.ge.c).and.(b.ge.d)) Then
  max_mass = b  
Else  If (c.ge.d) Then
  max_mass = c
Else
  max_mass = d
End if

ra = Abs(a/max_mass)
rb = Abs(b/max_mass)
rc = Abs(c/max_mass)
rd = Abs(d/max_mass)



 If ((ra.le.eps).and.(rb.le.eps)) Then
  If (Abs(c-d).le.eps) Then
   KilianC0D0 = 1/c
  Else 
   KilianC0D0 = 1./(c-d)*Log(d/c)
  End if
 Else If ((ra.le.eps).and.(rc.le.eps)) Then
  If (Abs(b-d).le.eps) Then
   KilianC0D0 = 1/b
  Else 
   KilianC0D0 = 1./(b-d)*Log(d/b)
  End if
 Else If ((ra.le.eps).and.(rd.le.eps)) Then
  If (Abs(b-c).le.eps) Then
   KilianC0D0 = 1/c
  Else 
   KilianC0D0 = 1./(b-c)*Log(c/b)
  End if
 Else If ((rb.le.eps).and.(rc.le.eps)) Then
  If (Abs(a-d).le.eps) Then
   KilianC0D0 = 1/a
  Else 
   KilianC0D0 = 1./(a-d)*Log(d/a)
  End if
 Else If ((rb.le.eps).and.(rd.le.eps)) Then
  If (Abs(a-c).le.eps) Then
   KilianC0D0 = 1/c
  Else 
   KilianC0D0 = 1./(a-c)*Log(c/a)
  End if
 Else If ((rc.le.eps).and.(rd.le.eps)) Then
  If (Abs(a-b).le.eps) Then
   KilianC0D0 = 1/a
  Else 
   KilianC0D0 = 1./(a-b)*Log(b/a)
  End if
 Else if ((Abs(a-b)/max_mass.lt.eps).and.(rc.lt.eps)) Then
   KilianC0D0 = (-b + d + d*Log(b/d))/(b - d)**2
 Else if ((Abs(a-b)/max_mass.lt.eps).and.(rd.lt.eps)) Then
   KilianC0D0 = (-a + c + c*Log(a/c))/(a - c)**2
 Else if ((Abs(a-c)/max_mass.lt.eps).and.(rd.lt.eps)) Then
   KilianC0D0 = (-a + b + b*Log(a/b))/(a - b)**2
 Else if ((Abs(a-c)/max_mass.lt.eps).and.(rb.lt.eps)) Then
   KilianC0D0 = (-a + d + d*Log(a/d))/(a - d)**2
 Else if ((Abs(a-d)/max_mass.lt.eps).and.(rb.lt.eps)) Then
   KilianC0D0 = (c - d + c*Log(d/c))/(c - d)**2
 Else if ((Abs(a-d)/max_mass.lt.eps).and.(rc.lt.eps)) Then
   KilianC0D0 = (b - d + c*Log(d/b))/(b - d)**2
 Else if ((Abs(b-c)/max_mass.lt.eps).and.(rd.lt.eps)) Then
   KilianC0D0 = (a - c + a*Log(c/a))/(a - c)**2
 Else if ((Abs(b-c)/max_mass.lt.eps).and.(ra.lt.eps)) Then
   KilianC0D0 = (d - c + d*Log(c/d))/(d - c)**2
 Else if ((Abs(b-d)/max_mass.lt.eps).and.(ra.lt.eps)) Then
   KilianC0D0 = (c - d + c*Log(d/c))/(c - d)**2
 Else if ((Abs(b-d)/max_mass.lt.eps).and.(rc.lt.eps)) Then
   KilianC0D0 = (a - d + a*Log(d/a))/(a - d)**2
 Else if ((Abs(c-d)/max_mass.lt.eps).and.(ra.lt.eps)) Then
   KilianC0D0 = (b - d + b*Log(b/d))/(b - d)**2
 Else if ((Abs(c-d)/max_mass.lt.eps).and.(rb.lt.eps)) Then
   KilianC0D0 = (a - d + a*Log(a/d))/(a - d)**2
 Else if ((Abs(a-b)/max_mass.lt.eps).and.(Abs(c-d)/max_mass.lt.eps)) Then
   KilianC0D0 = (-b**2 + c**2 + 2*b*c*Log(b/c))/(b - c)**3
 Else
  C0user = C0_3m(a,b,c)
  D0user = KilianD0(a,b,c,d)
  KilianC0D0 = C0user+d*D0user
 End If
! End if
End Function 

Real(dp) Function KilianB0C0(a,b,c)
 ! This function is B0(0,a,b)+c*C0_3m(a,b,c)
 ! it catches the case a=b=0
 ! cutoff for masses = sqrt(10^-10) GeV
 Implicit None
 Real(dp), Intent(in) :: a,b,c
 Real(dp) :: B0user, C0user

 If ((a.le.1.E-10_dp).and.(b.le.1.E-10_dp)) Then
  KilianB0C0 = 1._dp - Log(c/getRenormalizationScale())
 Else
  KilianB0C0 = B0(0._dp,a,b)+c*C0_3m(a,b,c)
 End If
End Function 


!! Arrange order masses in C-functions if two masses are the same

Complex(dp) Function vertexC11aux(a,b,c)
Implicit None
Real(dp), Intent(in) :: a,b,c
If (a.eq.b) Then
vertexC11aux = vertexC11(c,b,a)
Else 
vertexC11aux = vertexC11(a,b,c)
end if
End function vertexC11aux
Complex(dp) Function vertexC12aux(a,b,c)
Implicit None
Real(dp), Intent(in) :: a,b,c
If (a.eq.b) Then
vertexC12aux = vertexC12(c,b,a)
Else 
vertexC12aux = vertexC12(a,b,c)
end if
End function vertexC12aux
Complex(dp) Function vertexC0tildeaux(a,b,c)
Implicit None
Real(dp), Intent(in) :: a,b,c
If (a.eq.b) Then
vertexC0tildeaux = vertexC0tilde(c,b,a) + divergence
Else 
vertexC0tildeaux = vertexC0tilde(a,b,c) + divergence
end if
End function vertexC0tildeaux


 Subroutine InitializeLoopFunctions()
  Real(dp) :: s, cs, sold
  Integer :: i1

  s = 1._dp
  xalogm = 1._dp
  Do i1=1,10000
   s = 0.5_dp * s
   If ( (2._dp*s).Ne.xalogm ) Exit
   xalogm = Abs( s )
  End Do
  If (xalogm.Eq.0._dp) xalogm = Tiny(1._dp)
  s = 1._dp
  xclogm = (1._dp,0._dp)
  Do i1=1,10000
   s = 0.5_dp * s
   If ( (2._dp*Abs(Cmplx(s,0._dp,dp))).Ne.xclogm ) Exit
   xclogm = Abs( Cmplx(s,0._dp,dp) )
  End Do
  If (xclogm.Eq.0._dp) xclogm = Tiny(1._dp)
  xalog2 = Sqrt(xalogm)
  xclog2 = Sqrt(xclogm)

  ! the precision to which real calculations are done is
  precx = 1._dp
  sold = 0._dp
  Do i=1,1000
   precx = precx/2
   Call addone(s, precx)
   s = Exp(Log(s))
   If ( s .Eq. sold ) Exit
   sold = s
  End Do
  precx = precx*8
  ! the precision to which complex calculations are done is
  precc = 1._dp
  sold = 0._dp
  Do i=1,1000
   precc = precc/2
   cs = Exp(Log(Cmplx(1._dp+precc,0._dp,dp)))
   If ( Real(cs,dp) .Eq. sold ) Exit
   sold = Real(cs,dp)
  End Do
  precc = precc*8
 !
 !	for efficiency tke them equal   If they are not too different
 !
  If ( precx/precc < 4._dp .And. precx/precc > .25_dp ) Then
    precx = Max(precc,precx)
    precc = Max(precc,precx)
  End If

  ! some constants
  xlg2 = Log(2._dp)
  czero =(0._dp, 0._dp)
  cone=(1._dp, 0._dp)
  c2ipi = Cmplx(0._dp, 2._dp * Pi, dp)
  cipi2 = Cmplx(0._dp, Pi**2, dp)

 !
 !	inverses of   Integers:
 !
  Do i=1,30
    xninv(i) = 1._dp/i
    xn2inv(i) = 1._dp/(i*i)
  End Do
 !
 !	inverses of faculties of   Integers:
 !
  xinfac(1) = 1._dp
  Do i=2,30
    xinfac(i) = xinfac(i-1)/i
  End Do
 !
 !	inx: p(inx(i,j)) = isgn(i,j)*(s(i)-s(j))
 !
  inx(1,1) = -9999
  inx(2,1) = 5
  inx(3,1) = 9
  inx(4,1) = 8
  inx(1,2) = 5
  inx(2,2) = -9999
  inx(3,2) = 6
  inx(4,2) = 10
  inx(1,3) = 9
  inx(2,3) = 6
  inx(3,3) = -9999
  inx(4,3) = 7
  inx(1,4) = 8
  inx(2,4) = 10
  inx(3,4) = 7
  inx(4,4) = -9999
  isgn(1,1) = -9999
  isgn(2,1) = +1
  isgn(3,1) = -1
  isgn(4,1) = -1
  isgn(1,2) = -1
  isgn(2,2) = -9999
  isgn(3,2) = +1
  isgn(4,2) = +1
  isgn(1,3) = +1
  isgn(2,3) = -1
  isgn(3,3) = -9999
  isgn(4,3) = +1
  isgn(1,4) = +1
  isgn(2,4) = -1
  isgn(3,4) = -1
  isgn(4,4) = -9999
  inx5(1,1) = -9999
  inx5(1,2) =  6
  inx5(1,3) = 11
  inx5(1,4) = 14
  inx5(1,5) = 10
  inx5(2,1) =  6
  inx5(2,2) = -9999
  inx5(2,3) =  7
  inx5(2,4) = 12
  inx5(2,5) = 15
  inx5(3,1) = 11
  inx5(3,2) =  7
  inx5(3,3) = -9999
  inx5(3,4) =  8
  inx5(3,5) = 13
  inx5(4,1) = 14
  inx5(4,2) = 12
  inx5(4,3) =  8
  inx5(4,4) = -9999
  inx5(4,5) =  9
  inx5(5,1) = 10
  inx5(5,2) = 15
  inx5(5,3) = 13
  inx5(5,4) =  9
  inx5(5,5) = -9999
 !	isgn5 is not yet used.
  Do i=1,5
      Do j=1,5
      isgn5(i,j) = -9999
      End Do
  End Do
 !
  inx6(1,1) = -9999
  inx6(1,2) =  7
  inx6(1,3) = 13
  inx6(1,4) = 19
  inx6(1,5) = 17
  inx6(1,6) = 12
  inx6(2,1) =  7
  inx6(2,2) = -9999
  inx6(2,3) =  8
  inx6(2,4) = 14
  inx6(2,5) = 20
  inx6(2,6) = 18
  inx6(3,1) = 13
  inx6(3,2) =  8
  inx6(3,3) = -9999
  inx6(3,4) =  9
  inx6(3,5) = 15
  inx6(3,6) = 21
  inx6(4,1) = 19
  inx6(4,2) = 14
  inx6(4,3) =  9
  inx6(4,4) = -9999
  inx6(4,5) = 10
  inx6(4,6) = 16
  inx6(5,1) = 17
  inx6(5,2) = 20
  inx6(5,3) = 15
  inx6(5,4) = 10
  inx6(5,5) = -9999
  inx6(5,6) = 11
  inx6(6,1) = 12
  inx6(6,2) = 18
  inx6(6,3) = 21
  inx6(6,4) = 16
  inx6(6,5) = 11
  inx6(6,6) = -9999
 !	isgn6 is used.
  Do i=1,6
      Do j=1,6
      ji = j-i
        If ( ji>+3 ) ji = ji - 6
        If ( ji<-3 ) ji = ji + 6
        If ( ji==0 ) Then
        isgn6(j,i) = -9999
        Else If ( Abs(ji)==3 ) Then
          If ( i<0 ) Then
          isgn6(j,i) = -1
          Else
          isgn6(j,i) = +1
          End If
        Else If ( ji>0 ) Then
        isgn6(j,i) = +1
        Else If ( ji<0 ) Then
        isgn6(j,i) = -1
        Else
        Write(ErrCan,*) "InitialzeLoopFunctions: internal error in isgn6"
        Call TerminateProgram()
        End If
      End Do
  End Do
!
 !      calculate the coefficients of the series expansion
 !      li2(x) = sum bn*z^n/(n+1)!, z = -log(1-x), bn are the
 !      bernouilli numbers (zero for odd n>1).
 !
  bf(1) = - 1._dp/4._dp
  bf(2) = + 1._dp/36._dp
  bf(3) = - 1._dp/36.e2_dp
  bf(4) = + 1._dp/21168.e1_dp
  bf(5) = - 1._dp/108864.e2_dp
  bf(6) = + 1._dp/52690176.e1_dp
  bf(7) = - 691._dp/16999766784.e3_dp
  bf(8) = + 1._dp/1120863744.e3_dp
  bf(9) = - 3617._dp/18140058832896.e4_dp
  bf(10) = + 43867._dp/97072790126247936.e3_dp
  bf(11) = - 174611._dp/168600109166641152.e5_dp
  bf(12) = + 77683._dp/32432530090601152512.e4_dp
  bf(13) = - 236364091._dp/4234560341829359173632.e7_dp
  bf(14) = + 657931._dp/5025632054039239458816.e6_dp
  bf(15) = - 3392780147._dp/109890470493622010006470656.e7_dp
  bf(16)=+172.3168255201_dp/2355349904102724211909.3102313472e6_dp
  bf(17)=-770.9321041217_dp/4428491985594062112714.2791446528e8_dp
  bf(18)=  4.1576356446138997196178996207752e-29_dp
  bf(19)= -9.9621484882846221031940067024558e-31_dp
  bf(20)=  2.3940344248961653005211679878937e-32_dp

  ! coefficients for series in B0
  B0serie1(1) = 1._dp / 6._dp
  B0serie2(1) = B0serie1(1)
  B0serie3(1) = 1._dp
  DB0serie1(1) = 1._dp
  DB0serie2(1) = B0serie1(1)
  Do i1=2,60
   B0serie1(i1) = - B0serie1(i1-1) * Real(i1-1,dp) / Real(2*(2*i1+1),dp)
   B0serie2(i1) = Real(i1,dp) / Real((i1+1)*(i1+2), dp)
   B0serie3(i1) = B0serie3(i1-1) / Real(i1,dp)
   DB0serie1(i1) = 1._dp / Real(i1,dp)
   DB0serie2(i1) = - DB0serie2(i1-1) * Real(i1,dp) / Real(2*(2*i1+1),dp)
  End Do
  DB0serie3 = B0serie3

  ! error messages
  St_error(1) = "Numerical Problems in A0"
  St_error(2) = "Problem in Function B0, divergence for k->0, m1=m2=0"
  St_error(3) = "Numerical problem in Function B0, log for equal masses"
  St_error(4) = "Numerical problem in Function B0, log for unequal masses"
  St_error(5) = &
     & "Numerical problem in Function B0, min. value for log, change mu "
  St_error(6) = "Error in DB0, in case of IR divergence lambda must not be 0"
  St_error(7) = "Error in DB0, lam = 0, derivative cannot be computed"
  St_error(8) = &
     & "Numerical problem in Function B0, min. value for log, change mu "
  St_error(9) = "DB0, warning: cancellations in case p^2=0"
  St_error(10) = "Numerical problem in Function DB0, log for unequal masses"
  St_error(11) = "Numerical problem in Function DB0, cancellatios"
  St_error(12) = "DB0: error: arg log would be 0!"
  St_error(13) = "Rot3: all three external masses zero !"
  St_error(14) = "Rota: illegal flag (should not occur)"
  St_error(15) = "C0: lambda(p1,p2,p3) < 0, unphysical configuration"
  St_error(16) = &
     & "C0: cannot handle this case (p1,p2,p3 dependent, on threshold)"
  St_error(17) = "Minimum value complex log caues correction term to be wrong."
  St_error(18) = "zfflog: imaginary part undefined for real z < 0."
  St_error(19) = "ffcc0p: cannot handle two spacelike momenta and one zero."
  St_error(20) = "ffcs3: illegal code for isoort(1) (should not occur)"
  St_error(21) = "ffcs3: illegal code for isoort(2) (should not occur)"
  St_error(22) = "ffcs3: isoort = -1,0 not yet ready"
  St_error(23) = "ffcs3: eta changes within (0,1), add sophisticated terms..."
  St_error(24) = "nffet1: eta is not defined for real negative numbers a,b,ab."
  St_error(25) = "nffeta: eta is not defined for real negative numbers a,b,ab."
  St_error(26) = "ffcxr: illegal code for iclas1 (should not occur)"
  St_error(27) = "ffcrr: illegal code for iclas1 (should not occur)"
  St_error(28) = "ffcrr: illegal code for iclas2 (should not occur)"
  St_error(29) = "zxfflg: imaginary part undefined for x < 0."
  St_error(30) = "ffcxr: illegal code for iclas2 (should not occur)"
  St_error(31) = "ffdcxr: dyz=0, should not occur"
  St_error(32) = "ffdcxr: Taylor expansion in 1/x not yet ready"
  St_error(33) = "ffdcrr: cdwz=0, but iepsz!=iepsz and significant"
  St_error(34) = "ffdcrr: cdyz=0, should not occur"
  St_error(35) = "ffdcrr: Taylor expansion in 1/x not yet ready"
  St_error(36) = "fxc0p: imaginary part Ai < 0 terms uncertain"
  St_error(37) = "fxc0i: IR divergent C0 with lambda2=0."
  St_error(38) = "ffxc0j: IR divergent C0 with lambda(p1,p2,p3)=0."
  St_error(39) = "ffxc0j: orry, complex roots not yet supported here"
  St_error(40) = &
     & "ffxc0j: ill-defined IR-divergent C0 for massless charged particles."
  St_error(41) = "ffzzdl: imaginary part dilog is undefined for real x > 1."
  St_error(42) = "ffxxyz: xk = 0 not yet implemented"
  ! set check variable
  initialized = .True.
 Contains
  Subroutine addone(res, x)
  ! We need this stupid subroutine so that Fortran is FORCED to put
  ! it's result in memory and we don't end up with the precision of
  ! the FPU registers in precx.
  Implicit None
   Real(dp), Intent(in) :: x
   Real(dp), Intent(out) :: res
   res = x + 1._dp
  End Subroutine addone

 End Subroutine InitializeLoopFunctions
    
! Interface Log1minusXpXn

 Real(dp) Function Kappa2p(a1,a2,a3,a12,a13,a23)
 !------------------------------------------------------
 ! numerical more accurate calculation of
 ! kappa = (a1-a2-a3)**2 - 4*a2*a3
 !------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: a1,a2,a3,a12,a13,a23

  Real(dp) :: aa1, aa2, aa3, a, aff, asq

  aa1 = Abs(a1)
  aa2 = Abs(a2)
  aa3 = Abs(a3)

!  first see if there are input parameters with opposite sign:
  If ( ((a1.Lt.0).And.(a2.Gt.0)) .Or. ((a1.Gt.0).And.(a2.Lt.0)) ) Then
   If ( aa1 .Gt. aa2 ) Then
    a = a13 + a2
   Else
    a = a1 + a23
   End If
   aff = 4*a1*a2
  Else If ( ((a1.Lt.0).And.(a3.Gt.0)) .Or. ((a1.Gt.0).And.(a3.Lt.0)) ) Then
   If ( aa1 .Gt. aa3 ) Then
    a = a12 + a3
   Else
    a = a1 - a23
   End If
  aff = 4*a1*a3
!  all have the same sign, choose the smallest 4*ai*aj term
  Else If ( (aa1.Gt.aa2) .And. (aa1.Gt.aa3) ) Then
   If ( aa2 .Gt. aa3 ) Then
    a = a12 - a3
   Else
    a = a13 - a2
   End If
   aff = 4*a2*a3
  Else If ( aa2 .Gt. aa3 ) Then
   If ( aa1 .Gt. aa3 ) Then
    a = a12 + a3
   Else
    a = a1 - a23
   End If
   aff = 4*a1*a3
  Else
   If ( aa1 .Gt. aa2 ) Then
    a = a13 + a2
   Else
    a = a1 + a23
   End If
   aff = 4*a1*a2
  End If
  asq = a**2
  Kappa2p = asq - aff
  
 End Function Kappa2p

 Complex(dp) Function Log1minusX_c(x)
 !-------------------------------------------------
 ! calculates the series for for ln(1-x) 
 ! written by Werner Porod, 02.04.02
 !-------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: x
  
  Integer :: i1, nmax
  Complex(dp) :: test

  If (x.Eq.(1._dp,0._dp)) Then
   Log1minusX_c = Huge(1._dp)
  Else If (Abs(x).Gt.0.125_dp) Then
   Log1minusX_c = Log(1._dp - x )
  Else If (x.Eq.(0._dp,0._dp)) Then
   Log1minusX_c = 0._dp
  Else
   test = - x
   Log1minusX_c = - x
   i1 = 1
   ! find first maxmal power which is impartant
   Do 
    i1 = i1 + 1
    Log1minusX_c = Log1minusX_c - x**i1 / Real(i1,dp)
    If (Log1minusX_c.Eq.test) Exit
    test = Log1minusX_c
   End Do
   nmax = i1
   Log1minusX_c=0._dp
   Do i1 = nmax,1,-1
    Log1minusX_c = Log1minusX_c - x**i1 / Real(i1,dp)
   End Do

  End If

 End Function Log1minusX_c


! Interface Log1minusXpXn
 Complex(dp) Function Log1minusXpXn_c(x,n)
 !--------------------------------------------------------
 ! calculates the complex series for for ln(1-x) + sum_1^n x^n/n 
 ! written by Werner Porod, 08.04.02
 !--------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: n
  Complex(dp), Intent(in) :: x
  
  Integer :: i1, nmax
  Complex(dp) :: test

  If ((Real(x,dp).Eq.1._dp).And.(Aimag(x).Eq.0._dp)) Then
   Log1minusXpXn_c = Huge(1._dp)
  Else If (Abs(x).Gt.0.125_dp) Then
   Log1minusXpXn_c= Log(1._dp - x )
   Do i1=1,n
    Log1minusXpXn_c = Log1minusXpXn_c + x**i1 / Real(i1,dp)
   End Do
  Else If (Abs(x).Eq.0._dp) Then
   Log1minusXpXn_c = 0._dp
  Else
   i1 = n+1
   test = - (x)**(i1) / Real(i1, dp)
   Log1minusXpXn_c = test 
   ! find first maximal power which is important
   Do 
    i1 = i1 + 1
    Log1minusXpXn_c = Log1minusXpXn_c - x**i1 / Real(i1,dp)
    If (Log1minusXpXn_c.Eq.test) Exit
    test = Log1minusXpXn_c
   End Do
   nmax = i1
   Log1minusXpXn_c=0._dp
   Do i1 = nmax,n+1,-1
    Log1minusXpXn_c = Log1minusXpXn_c - x**i1 / Real(i1,dp)
   End Do

  End If

 End Function Log1minusXpXn_c


 Real(dp) Function Log1minusXpXn_r(x,n)
 !--------------------------------------------------------
 ! calculates the series for for ln(1-x) + sum_1^n x^n/n 
 ! written by Werner Porod, 05.04.02
 !--------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: n
  Real(dp), Intent(in) :: x
  
  Integer :: i1, nmax
  Real(dp) :: test

  If (x.Eq.1._dp) Then
   Log1minusXpXn_r = Huge(1._dp)
  Else If (Abs(x).Gt.0.125_dp) Then
   Log1minusXpXn_r= Log(1._dp - x )
   Do i1=1,n
    Log1minusXpXn_r = Log1minusXpXn_r + x**i1 / Real(i1,dp)
   End Do
  Else If (x.Eq.0._dp) Then
   Log1minusXpXn_r = 0._dp
  Else
   i1 = n+1
   test = - (x)**(i1) / Real(i1, dp)
   Log1minusXpXn_r = test 
   ! find first maximal power which is important
   Do 
    i1 = i1 + 1
    Log1minusXpXn_r = Log1minusXpXn_r - x**i1 / Real(i1,dp)
    If (Log1minusXpXn_r.Eq.test) Exit
    test = Log1minusXpXn_r
   End Do
   nmax = i1
   Log1minusXpXn_r=0._dp
   Do i1 = nmax,n+1,-1
    Log1minusXpXn_r = Log1minusXpXn_r - x**i1 / Real(i1,dp)
   End Do

  End If

 End Function Log1minusXpXn_r


 Real(dp) Function Log1minusX_r(x)
 !-------------------------------------------------
 ! calculates the series for for ln(1-x) 
 ! written by Werner Porod, 02.04.02
 !-------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x
  
  Integer :: i1, nmax
  Real(dp) :: test

  If (x.Eq.1._dp) Then
   Log1minusX_r = Huge(1._dp)
  Else If (Abs(x).Gt.0.125_dp) Then
   Log1minusX_r = Log(1._dp - x )
  Else If (x.Eq.0._dp) Then
   Log1minusX_r = 0._dp
  Else
   test = - x
   Log1minusX_r = - x
   i1 = 1
   ! find first maxmal power which is impartant
   Do 
    i1 = i1 + 1
    Log1minusX_r = Log1minusX_r - x**i1 / Real(i1,dp)
    If (Log1minusX_r.Eq.test) Exit
    test = Log1minusX_r
   End Do
   nmax = i1
   Log1minusX_r=0._dp
   Do i1 = nmax,1,-1
    Log1minusX_r = Log1minusX_r - x**i1 / Real(i1,dp)
   End Do

  End If

 End Function Log1minusX_r


 Subroutine LoopFunctionsErrorSummary
  Logical :: NoErrors = .True.
  Integer :: i1

  Write (ErrCan,*) "Report on Numerical Problems within the LoopFunctions:"
  If (PrintToScreen) &
      & Write (*,*) "Report on Numerical Problems within the LoopFunctions:"

  Do i1=1,N_err
   If (ErrorOccur(i1).Ne.0) Then
    Write(ErrCan,*) Trim(St_error(i1)), ErrorOccur(1)
    If (PrintToScreen) Write(*,*) Trim(St_error(i1)), ErrorOccur(1)
    NoErrors = .False.
   End If
  End Do

  If (NoErrors) Then
   Write(ErrCan,*) "No errors or warnings."
   If (PrintToScreen) Write(*,*) "No errors or warnings."
  End If

 End Subroutine LoopFunctionsErrorSummary

 Integer Function nffet1(ca,cb,cc,ier)
 !--------------------------------------------------------------------
 !	calculates the same eta with three input variables
 !
 !	et1(a,b)/(2*i*pi) = ( thIm(-a)*thIm(-b)*thIm(c)
 !				- thIm(a)*thIm(b)*thIm(-c) )
 !
 !	with thIm(a) = theta(Im(a))
 !--------------------------------------------------------------------
 Implicit None
  Integer :: ier
  Complex(dp) :: ca,cb,cc
  Real(dp) :: a,b,ab

 !  #[ calculations:
  a = Aimag(ca)
  b = Aimag(cb)
  If ( a > 0 .Neqv. b > 0 ) Then
    nffet1 = 0
    Return
  End If
  ab = Aimag(cc)
  If ( a < 0 .And. b < 0 .And. ab > 0 ) Then
    nffet1 = 1
  Else If ( a > 0 .And. b > 0 .And. ab < 0 ) Then
    nffet1 = -1
  Else If ( a == 0 .And. Real(ca) <= 0 .Or.b == 0 .And. Real(cb) <= 0 .Or.&
      &ab == 0 .And. Real(cc) <= 0 ) Then
      Call WriteLFerror(24)
    nffet1 = 1
  Else
    nffet1 = 0
  End If
 !  #] calculations:
 End Function nffet1

 Integer Function nffeta(ca,cb,ier)
 !--------------------------------------------------------------------
 !	calculates
 !
 !	eta(a,b)/(2*i*pi) = ( thIm(-a)*thIm(-b)*thIm(a*b)
 !				- thIm(a)*thIm(b)*thIm(-a*b) )
 !
 !	with thIm(a) = theta(Im(a))
 !--------------------------------------------------------------------
 Implicit None
  Integer :: ier
  Complex(dp) :: ca,cb
  Real(dp) :: a,b,ab,rab

 !  #[ calculations:
  a = Aimag(ca)
  b = Aimag(cb)
  If ( a*b < 0 ) Then
    nffeta = 0
    Return
  End If
  rab = Real(ca,dp)*Real(cb,dp) - a*b
  ab = Real(ca,dp)*b + a*Real(cb,dp)
  If ( Abs(ab) < precc*Abs(Real(ca,dp)*b) ) Then
      Call WriteLFerror(25)
  End If
  If ( a < 0 .And. b < 0 .And. ab > 0 ) Then
    nffeta = 1
  Else If ( a > 0 .And. b > 0 .And. ab < 0 ) Then
    nffeta = -1
  Else If ( a == 0 .And. Real(ca) <= 0 .Or.b == 0 .And. Real(cb) <= 0 .Or.&
      &ab == 0 .And. rab <= 0 ) Then
      Call WriteLFerror(25)
    nffeta = 0
  Else
    nffeta = 0
  End If
 !  #] calculations:
 End Function nffeta

 Real(dp) Function phi(x,y,z)
 !--------------------------------------------------------------
 ! from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23
 ! version by Pietro Slavich
 ! 12.03.02: portation to f90
 !--------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x, y, z
  
  If(x.Le.z.And.y.Le.z) Then
     phi = myphi(x,y,z)
  Elseif(z.Le.x.And.y.Le.x) Then
     phi = z/x*myphi(z,y,x)
  Elseif(z.Le.y.And.x.Le.y) Then
     phi = z/y*myphi(z,x,y)
  Endif

 Contains

  Real(dp) Function myphi(x,y,z)
   Implicit None
   Real(dp), Intent(in) :: x,y,z
   Real(dp) :: u,v
   Complex(dp) :: clam,cxp,cxm,ccphi

   !     auxiliary variables
   u = x/z
   v = y/z
   If(u<=1.e-8_dp) Then
     
     If(v/=1._dp) Then
       myphi = (Log(u)*Log(v)+2._dp*Li2(1._dp-v))/(1._dp-v)
     Else
       myphi = 2._dp-Log(u)
     Endif

   Elseif(v<=1.e-8_dp) Then
     If(u/=1._dp) Then
       myphi = (Log(v)*Log(u)+2._dp*Li2(1._dp-u))/(1._dp-u)
     Else
       myphi = 2._dp-Log(v)
     Endif
   Else
     
    If((1._dp-u-v)**2>=4._dp*u*v) Then         
      clam = Cmplx(Sqrt((1._dp-u-v)**2 - 4._dp*u*v),0._dp,dp)
    Else
      clam = Cmplx(0._dp,Sqrt(4._dp*u*v - (1._dp-u-v)**2),dp)
    Endif
    cxp = (1._dp+(u-v)-clam)/2._dp
    cxm = (1._dp-(u-v)-clam)/2._dp
       
   !     phi function from eq. (A4)
       
    ccphi = (2._dp*Log(cxp)*Log(cxm) - Log(u)*Log(v)  &
        & - 2._dp*(CLI2(cxp) + CLI2(cxm)) + Pi2o3)/clam
    myphi = Real(ccphi,dp)
       
   Endif
   Return

  End Function myphi
  
 End Function phi  

 Subroutine roots_c2(a,b,c,d,xm,xp)
 !--------------------------------------------------------------------
 ! calculates the roots for the equation
 !   a x^2 - 2 b x + c = 0
 ! xp = (b + d )/a and xm = (b - d )/ a
 ! with d = Sqrt( - 4*a*c + b**2)
 !--------------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: a, b, c, d
  Complex(dp), Intent(out) :: xp, xm

  Complex(dp) :: cc

  If (a.Eq.0._dp) Then
   If (ErrorLevel.Ge.-1) Then
    Write (ErrCan,*) "Error in roots_c2: a=0"
    If (ErrorLevel.Ge.1) Call TerminateProgram
   End If
   If ((Real(b).Gt.0._dp).Eqv.(Real(d).Gt.0._dp)) Then
    xp = Huge(1._dp)
    xm = - c / b
   Else
    xm = Huge(1._dp)
    xp = - c / b
   End If
  End If

  cc = b+d

  If (d.Eq.0._dp) Then
   xm = b / a
   xp = xm
  Else If (   (Abs(Real(cc,dp))+Abs(Aimag(cc))) .Gt. &
          &  xloss*(Abs(Real(d,dp))+Abs(Aimag(d))) ) Then
   xp = ( b + d ) / a
   xm = c / (a*xp)
  Else
   xm = ( b - d ) / a
   xp = c / (a*xm)
  End If

 End Subroutine roots_c2


 Subroutine roots_r2(a,b,c,d,xm,xp)
 !--------------------------------------------------------------------
 ! calculates the roots for the equation
 !   a x^2 - 2 b x + c = 0
 ! xp = (b + d )/a and xm = (b - d )/ a
 ! with d = Sqrt( - 4*a*c + b**2)
 !--------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: a, b, c, d
  Real(dp), Intent(out) :: xp, xm

  If (a.Eq.0._dp) Then
   If (ErrorLevel.Ge.-1) Then
    Write (ErrCan,*) "Error in roots_r: a=0"
    If (ErrorLevel.Ge.1) Call TerminateProgram
   End If
   If ((b.Gt.0._dp).Eqv.(d.Gt.0._dp)) Then
    xp = Huge(1._dp)
    xm = - c / b
   Else
    xm = Huge(1._dp)
    xp = - c / b
   End If
  End If

  If (d.Eq.0._dp) Then
   xm = b / a
   xp = xm
  Else If ((b.Gt.0._dp).Eqv.(d.Gt.0._dp)) Then
   xp = ( b + d ) / a
   xm = c / (a*xp)
  Else
   xm = ( b - d ) / a
   xp = c / (a*xm)
  End If

 End Subroutine roots_r2


! Interface roots
 Subroutine roots_r(a,b,c,xm,xp)
 !--------------------------------------------------------------------
 ! calculates the roots for the equation
 !   a x^2 - 2 b x + c = 0
 ! xp = (b + d )/a and xm = (b - d )/ a
 ! with d = Sqrt( - 4*a*c + b**2)
 !--------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: a, b, c
  Real(dp), Intent(out) :: xp, xm

  Real(dp) :: d
  
  If (a.Eq.0._dp) Then
   If (ErrorLevel.Ge.-1) Then
    Write (ErrCan,*) "Error in roots_r: a=0"
    If (ErrorLevel.Ge.1) Call TerminateProgram
   End If
   If (c/b.Gt.0._dp) Then
    xp = Huge(1._dp)
    xm = - c / b
   Else
    xm = Huge(1._dp)
    xp = - c / b
   End If
  End If

  d = b**2 - c * a
  If (d.Lt.0._dp) Then
   If (ErrorLevel.Ge.-1) Then
    Write (ErrCan,*) "Error in roots_r: d<0, d=",d
    If (ErrorLevel.Ge.1) Call TerminateProgram
   End If
   xp = 0._dp
   xm = 0._dp
  Else If (d.Eq.0._dp) Then
   xm = b / a
   xp = xm
  Else If ((b.Gt.0._dp).Eqv.(d.Gt.0._dp)) Then
   xp = ( b + d ) / a
   xm = c / (a*xp)
  Else
   xm = ( b - d ) / a
   xp = c / (a*xm)
  End If

 End Subroutine roots_r



 Subroutine SetLookUpCache(l_new)
 Implicit None
  Logical, Intent(in) :: l_new
  l_look_up_cache = l_new
  If (l_new) Call ClearCache()
 End Subroutine SetLookUpCache


 Logical Function SetPrintWarnings(lx)
 Implicit None
  Logical, Intent(in) :: lx
  SetPrintWarnings = PrintToScreen
  PrintToScreen = lx
 End Function SetPrintWarnings


 Real(dp) Function SetRenormalizationScale(mu2_in)
 !-----------------------------------------------------------------------
 ! resets renormalization scale, gives old scale as output for checking
 ! written by Werner Porod
 ! 15.11.01
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mu2_in

  SetRenormalizationScale = mudim2
  mudim2 = mu2_in

 End Function SetRenormalizationScale


 Subroutine WriteLFerror(n)
 Implicit None
  Integer, Intent(in) :: n
  ErrorOccur(n) = ErrorOccur(n) + 1
  If (ErrorLevel.Ge.-1) Then
   Write(ErrCan,*) Trim( St_error(n) )
   If (PrintToScreen) Write(*,*) Trim( St_error(n) )
   If (ErrorLevel.Ge.1) Call TerminateProgram
  End If
  
 End Subroutine WriteLFerror

  Complex(dp) Function zfflog(cx,ieps,cy)
  !----------------------------------------------------------------------
  !	Calculate the complex logarithm of cx.  The following cases
  !	are treted separately:
  !		|cx| too small:		give warning and return 0
  !					(for Absoft, Apollo DN300)
  !		Im(cx) = 0, Re(cx) < 0:	take sign according to ieps
  !----------------------------------------------------------------------
  Implicit None
    Integer :: ieps
    Complex(dp) :: cx,cy

    Complex(dp) :: ctroep
    Real(dp) :: xa, xlog1p
    !  #[ calculations:
    xa = absc(cx)
      If ( xa < xalogm ) Then
        If ( cx /= 0 )   Call WriteLFerror(17)
        zfflog = 0
      Else If ( Real(cx) < 0 .And. Aimag(cx) == 0 ) Then
    !     +		 abs(Aimag(cx)) .lt. precc*abs(Real(cx)) ) then
        xlog1p = Log(-Real(cx,dp))
    !	    checked imaginary parts 19-May-1988
        If ( Abs(ieps) == 1 ) Then
          If ( ieps*Real(cy) < 0 ) Then
           zfflog = Cmplx(xlog1p,-pi,dp)
          Else If ( ieps*Real(cy,dp) > 0 ) Then
           zfflog = Cmplx(xlog1p,pi,dp)
          Else
            Call WriteLFerror(18)
           zfflog = Cmplx(xlog1p,pi,dp)
          End If
        Else If ( ieps >= 2 .And. ieps <= 3 ) Then
         zfflog = Cmplx(xlog1p,-pi,dp)
        Else If ( ieps <= -2 .And. ieps >= -3 ) Then
         zfflog = Cmplx(xlog1p,pi,dp)
        Else
          Call WriteLFerror(18)
          zfflog = Cmplx(xlog1p,pi,dp)
        End If
      Else If ( xa < xclogm .Or. 1/xa < xclogm ) Then
       ctroep = cx*Real(1._dp/xa,dp)
       zfflog = Log(ctroep) + Real(Log(xa),dp)
      Else
    !	    print *,'zfflog: neem log van ',cx
       zfflog = Log(cx)
      End If
    !  #] calculations: 
    End Function zfflog

 Complex(dp) Function zxfflg(x,ieps,y)
 !----------------------------------------------------------------------
 !      Calculate the complex logarithm of x.  The following cases      *
 !      are treted separately:                                          *
 !              |x| too small:          give warning and return 0       *
 !                                      (for Absoft, Apollo DN300)      *
 !              |x| < 0:                take sign according to ieps     *
 !----------------------------------------------------------------------
 Implicit None
  Integer :: ieps
  Real(dp) :: x,y

  Real(dp) :: xlog

 !  #[ calculations:
  If ( Abs(x) < xalogm ) Then
    zxfflg = 0
  Else If ( x > 0 ) Then
    zxfflg = Log(x)
  Else
    xlog = Log(-x)
 !          checked imaginary parts 19-May-1988
      If ( Abs(ieps) == 1 ) Then
        If ( y*ieps < 0 ) Then
          zxfflg = Cmplx(xlog,-pi,dp)
        Else
          zxfflg = Cmplx(xlog,pi,dp)
        End If
      Else If ( ieps == 2 ) Then
        zxfflg = Cmplx(xlog,-pi,dp)
      Else If ( ieps == -2 ) Then
        zxfflg = Cmplx(xlog,+pi,dp)
      Else
        Call WriteLFerror(29)
        zxfflg = Cmplx(xlog,pi,dp)
      End If
  End If
 !  #] calculations:
 End Function zxfflg

!------------------------------------------------------------
! functions taken from C.Bobeth et al., NPB 630 (2002) 87 
! first index is index of functions whereas the second one determines if
! it is a lowest order function or a NLO function. Functions starting
! with D are derivatives of functions
!------------------------------------------------------------
 Real(dp) Function f_1_0(x)
!------------------------------------------------------------
! loop function taken from C.Bobeth et al., NPB 630 (2002) 87 
!------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: LogX, r
  Real(dp), Parameter :: ci(7) = (/ 0.75_dp, 5._dp/6._dp, -1._dp/24._dp &
                    &            , 0._dp, 1._dp/120._dp, -1./105._dp    &
                    &            , 1._dp/112._dp /)

  If (x.Eq.1._dp) Then
   f_1_0 = 0.75_dp
  Else If (Abs(x-1._dp).Lt.1.e-6_dp) Then
   r = x - 1._dp
   f_1_0 = ci(7) * r
   Do i1=5,1,-1
    f_1_0 = r * (f_1_0 + ci(i1) )
   End Do
   f_1_0 = f_1_0 + ci(1)

  Else If (x.Eq.0._dp) Then 
   f_1_0 = 0._dp

  Else If (Abs(x).Lt.1.e-6_dp) Then
   LogX = Log(x)
   f_1_0 = x * (3._dp + LogX                                        &
         &     +x * (2.5_dp + 3.5_dp * LogX                         &
         &          +x * (2.5_dp + 6._dp * LogX                     &
         &               +x * (2.5_dp + 8.5_dp * LogX               &
         &                    +x * (2.5_dp + 11._dp * LogX          &
         &                         + x* (2.5_dp + 13.5_dp * LogX ) ) ) ) ) )

  Else
   f_1_0 = ( (x-6._dp) + (2._dp + 3._dp * x) * Log(x) / (x-1._dp) ) &
         & * 0.5_dp * x / (x - 1._dp)
  End If
  
 End Function f_1_0


 Real(dp) Function Df_1_0(x)
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: r
  Real(dp), Parameter :: ci(7) = (/ 5._dp/6._dp, -1._dp/12._dp, 0._dp        &
                    &            , 1._dp/30._dp, -1./21._dp, 3._dp / 56._dp  &
                    &            , - 1._dp/18._dp /)

  If (x.Eq.1._dp) Then
   Df_1_0 = ci(1)
  Else If (Abs(x-1._dp).Lt.1.e-6_dp) Then
   r = x - 1._dp
   Df_1_0 = ci(7) * r
   Do i1=5,1,-1
    Df_1_0 = r * (Df_1_0 + ci(i1) )
   End Do
   Df_1_0 = Df_1_0 + ci(1)

  Else
   Df_1_0 = 0.5_dp * (-8._dp + 7._dp * x                    &
          &          - 2._dp * (1._dp + 4._dp * x) * Log(x) ) / (x - 1._dp)**3
  End If
  
 End Function Df_1_0

 Real(dp) Function f_1_1(x)
 Implicit none
  Real(dp), intent(in) :: x

   f_1_1 = 4._dp * x * (29._dp + 7._dp *x + 4._dp * x**2)             &
       &             / (3._dp * (x - 1._dp)**2 )                      &
       & - 4._dp * x * (23._dp + 14._dp *x + 3._dp * x**2) * Log(x)   &
       &             / (3._dp * (x - 1._dp)**3 )                      &
       & - 4._dp * x * (4._dp + x**2) * Li2(1._dp - 1._dp / x)        &
       &             / (x - 1._dp)**2

 End Function f_1_1

 Real(dp) Function f_2_0(x)
 !------------------------------------------------------------
 ! loop function taken from C.Bobeth et al., NPB 630 (2002) 87 
 !------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: LogX, r
  Real(dp), Parameter :: ci(7) = (/ -0.5_dp, -1._dp/6._dp, 1._dp/12._dp        &
                    &            , -1._dp/20._dp,  1._dp/30._dp, -1./42._dp    &
                    &            , 1._dp/56._dp /)

  If (x.Eq.1._dp) Then
   f_2_0 = ci(1)
  Else If (Abs(x-1._dp).Lt.1.e-6_dp) Then
   r = x - 1._dp
   f_2_0 = ci(7) * r
   Do i1=5,1,-1
    f_2_0 = r * (f_2_0 + ci(i1) )
   End Do
   f_2_0 = f_2_0 + ci(1)

  Else If (x.Eq.0._dp) Then 
   f_2_0 = 0._dp
  Else If (Abs(x).Lt.1.e-6_dp) Then
   LogX = Log(x)
   f_2_0 = 0._dp
   Do i1=6,1,-1
    f_2_0 = x * (f_2_0 + 1._dp + i1 * LogX )
   end do

  Else
   f_2_0 = ( -1._dp + Log(x) / (x-1._dp) ) * x / (x - 1._dp)
  End If
  
 End Function f_2_0


 Real(dp) Function Df_2_0(x)
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: r
  Real(dp), Parameter :: ci(7) = (/ -1._dp/6._dp, 1._dp/6._dp, -3._dp/20._dp  &
                    &            ,  2._dp/15._dp, -5./42._dp, 3._dp/28._dp    &
                    &            , - 7._dp / 72._dp /)

  If (x.Eq.1._dp) Then
   Df_2_0 = ci(1)
  Else If (Abs(x-1._dp).Lt.1.e-6_dp) Then
   r = x - 1._dp
   Df_2_0 = ci(7) * r
   Do i1=5,1,-1
    Df_2_0 = r * (Df_2_0 + ci(i1) )
   End Do
   Df_2_0 = Df_2_0 + ci(1)

  Else
   Df_2_0 = ( -2._dp + 2._dp * x - (1._dp + x ) * Log(x) ) / (x - 1._dp)**2
  End If
  
 End Function Df_2_0


 Real(dp) Function f_2_1(x)
 Implicit none
  Real(dp), intent(in) :: x

   f_2_1 = 32._dp * x * (3._dp - x)  / (3._dp * (x - 1._dp)**2 )      &
       & - 8._dp * x * (11._dp - 3._dp * x) * Log(x)                  &
       &             / (3._dp * (x - 1._dp)**3 )                      &
       & - 8._dp * x * (2._dp - x) * Li2(1._dp - 1._dp / x)           &
       &             / (x - 1._dp)**2

 End Function f_2_1

 Real(dp) Function f_3_0(x, y)
 !------------------------------------------------------------
 ! loop function taken from C.Bobeth et al., NPB 630 (2002) 87 
 !------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x, y

  Integer :: i1
  Real(dp) :: r

  If (x.Eq.y) Then
   If (y.Eq.1._dp) Then
    f_3_0 = 0.5_dp
   Else If (Abs(y-1._dp).Lt.1.e-6_dp) Then
    f_3_0 = 0._dp
    r = y - 1._dp
    Do i1=6,1,-1
     f_3_0 = r * (f_3_0 + (-1._dp)**i1 / Real(i1+2,dp))
    End Do
    f_3_0 = f_3_0 + 0.5_dp
   Else
    f_3_0 = (- 1._dp + y - Log(y) ) / (y-1._dp)**2
   End If

  Else If ((x.Eq.0._dp).And.(y.Eq.1._dp)) Then
   f_3_0 = 1._dp

  Else If ((x.Eq.1._dp).And.(y.Eq.0._dp)) Then
   f_3_0 = 1._dp

  Else If ((x.Eq.1._dp).And.(y.Eq.1._dp)) Then
   f_3_0 = 0.5_dp

  Else If (x.Eq.0._dp) Then
   f_3_0 = Log(y) / (y - 1._dp)

  Else If (y.Eq.0._dp) Then
   f_3_0 = Log(x) / (x - 1._dp)

  Else If (x.Eq.1._dp) Then
   f_3_0 = (1._dp - y + y * Log(y)) / (y - 1._dp)**2

  Else If (y.Eq.1._dp) Then
   f_3_0 = (1._dp - x + x * Log(x)) / (x - 1._dp)**2

  Else
   f_3_0 = (x * Log(x) / (x - 1._dp) -  y * Log(y) / (y - 1._dp)) / ( x- y)
  End If
  
 End Function f_3_0

 Real(dp) Function Dxf_3_0(x,y)
 Implicit None
  Real(dp), Intent(in) :: x, y

  Dxf_3_0 = (-((x**2 - y)*(-1 + y)*Log(x))                      &
        &    +(-1 + x)*((x - y)*(-1 + y) + (-1 + x)*y*Log(y)))  &
        &  /((-1 + x)**2*(x - y)**2*(-1 + y))

 End Function Dxf_3_0


 Real(dp) Function Dyf_3_0(x,y)
 Implicit None
  Real(dp), Intent(in) :: x, y

   Dyf_3_0 = ( x*(-1 + y)**2*Log(x)                                   & 
         &     -  (-1 + x)*((x - y)*(-1 + y) +  (-x + y**2)*Log(y)))  &
         &  /  ((-1 + x)*(x - y)**2*(-1 + y)**2)

 End Function Dyf_3_0

 Real(dp) Function F_3_1(x, y)
 Implicit None
  Real(dp), Intent(in) :: x, y

  ! some sign are absorbed by replacing e.g. 1-y by y-1

  F_3_1 = 28._dp * y / (3._dp * (x-y) * (1._dp - y) )                      &
      & + 2._dp * x * (11._dp * x + 3._dp * y) * Log(x)                    &
      &   / (3._dp * (x-y)**2 * (x-1._dp) )                                &
      & + 2._dp * y * (x*(25._dp-11._dp*y) - y*(11._dp+3._dp*y)) * Log(y)  &
      &   / (3._dp * (x-y)**2 * (y-1._dp)**2 )                             &
      & + 4._dp * (1._dp+y) * Li2(1._dp-1._dp/y) / ((x-1._dp)*(y-1._dp))   &
      & + 4._dp * (x+y) * Li2(1._dp-x/y) / ((x-1._dp)*(x-y))  
      
 End Function F_3_1


 Real(dp) Function f_4_0(x, y)
 !------------------------------------------------------------
 ! loop function taken from C.Bobeth et al., NPB 630 (2002) 87 
 !------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x, y

  Integer :: i1
  Real(dp) :: r

  If (x.Eq.y) Then
   If (y.Eq.0._dp) Then
    f_4_0 = 0._dp
   Else If (y.Eq.1._dp) Then
    f_4_0 = 1.5_dp
   Else If (Abs(y-1._dp).Lt.1.e-6_dp) Then
    f_4_0 = 0._dp
    r = y - 1._dp
    Do i1=6,1,-1
     f_4_0 = r * (f_4_0 + (-1._dp)**i1 / Real(i1+2,dp))
    End Do
    f_4_0 = f_4_0 + 0.5_dp
   Else
    f_4_0 = y * (-1._dp + y + (-2._dp + y) * Log(y) ) / (y-1._dp)**2
   End If

  Else If ((x.Eq.0._dp).And.(y.Eq.1._dp)) Then
   f_4_0 = 1._dp

  Else If ((x.Eq.1._dp).And.(y.Eq.0._dp)) Then
   f_4_0 = 1._dp

  Else If (x.Eq.0._dp) Then
   f_4_0 = Log(y) / (y - 1._dp)

  Else If (y.Eq.0._dp) Then
   f_4_0 = Log(x) / (x - 1._dp)

  Else If (x.Eq.1._dp) Then
   f_4_0 = (1._dp - y + y**2 * Log(y)) / (y - 1._dp)**2

  Else If (y.Eq.1._dp) Then
   f_4_0 = (1._dp - x + x**2 * Log(x)) / (x - 1._dp)**2

  Else
   f_4_0 = (x**2 * Log(x) / (x - 1._dp) -  y**2 * Log(y) / (y - 1._dp)) / ( x- y)
  End If
  
 End Function f_4_0


 Real(dp) Function Dxf_4_0(x, y)
 Implicit None
  Real(dp), Intent(in) :: x, y

  Dxf_4_0 = (-(x*(-1 + y)*(x - 2*y + x*y)*Log(x))                       &
        &    + (-1 + x)*(x*(x - y)*(-1 + y) + (-1 + x)*y**2*Log(y)))    &
        & /  ((-1 + x)**2*(x - y)**2*(-1 + y))

 End Function Dxf_4_0


 Real(dp) Function Dyf_4_0(x, y)
 Implicit None
  Real(dp), Intent(in) :: x, y

  Dyf_4_0 = (x**2*(-1 + y)**2*Log(x)                                     &
        &    - (-1 + x)*y*((x - y)*(-1 + y) + (x*(-2 + y) + y)*Log(y)))  &
        & /  ((-1 + x)*(x - y)**2*(-1 + y)**2)

 End Function Dyf_4_0

! dummy version
 Real(dp) Function F_4_1(x, y)
 Implicit None
  Real(dp), Intent(in) :: x, y

  F_4_1 = x+y

 End Function F_4_1


 Real(dp) Function f_5_0(x, y, z)
 !------------------------------------------------------------
 ! loop function taken from C.Bobeth et al., NPB 630 (2002) 87 
 !------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x, y, z

  Real(dp) :: yy, zz

   If ((x.Eq.y).And.(y.Eq.z)) Then
    If (z.Eq.1._dp) Then
     f_5_0 = 1._dp / 3._dp
    Else
     f_5_0 = (3._dp - 4._dp * z + z**2 + 2._dp *Log(z))/(2._dp*(z-1._dp)**3)
    End If

   Else If ((x.Eq.y).Or.(x.Eq.z).Or.(y.Eq.z)) Then
    If (x.Eq.y) Then
     yy = y
     zz = z
    Else If (x.Eq.z) Then
     yy = z
     zz = y
    Else If (y.Eq.z) Then
     yy = y
     zz = x
    End If

    If (yy.Eq.1._dp) Then
     If (zz.Eq.1._dp) Then
      F_5_0 = 1._dp / 3._dp
     Else
      F_5_0 = (-1._dp + 4._dp*zz - 3._dp*zz**2 + 2._dp*zz**2*Log(zz)) &
           & /(2._dp*(-1._dp + zz)**3)
     End If

    Else
     If (zz.eq.1._dp) then
      F_5_0 = (-1._dp + yy**2 - 2._dp*yy*Log(yy))/(-1._dp + yy)**3
     Else
      F_5_0 = (-(yy*(-1._dp + zz)*(yy - 2._dp*zz + yy*zz)*Log(yy)) +    &
         &    (-1._dp + yy)*(yy*(yy - zz)*(-1._dp + zz) + (-1._dp + yy)*zz**2*Log(zz)))/&
         &    ((-1._dp + yy)**2*(yy - zz)**2*(-1._dp + zz))
     End If
    End If

   Else If (x.Eq.1._dp) Then
    f_5_0 = 1._dp / ( 1._dp - y - z + y*z)              &
        & + y**2 * Log(y) / ( (y - 1._dp)**2 * (y-z) )   &
        & + z**2 * Log(z) / ( (z - 1._dp)**2 * ( z- y) )
   Else If (y.Eq.1._dp) Then
    f_5_0 = 1._dp / ( 1._dp - x - z + x*z)              &
        & + x**2 * Log(x) / ( (x - 1._dp)**2 * (x-z) )   &
        & + z**2 * Log(z) / ( (z - 1._dp)**2 * ( z- x) )
   Else If (z.Eq.1._dp) Then
    f_5_0 = 1._dp / ( 1._dp - y - x + y*x)              &
        & + y**2 * Log(y) / ( (y - 1._dp)**2 * (y-x) )   &
        & + x**2 * Log(x) / ( (x - 1._dp)**2 * ( x- y) )
   Else
    f_5_0 = x**2 * Log(x) / ( (x - 1._dp) * ( x- y) * (x-z))   &
        & + y**2 * Log(y) / ( (y - 1._dp) * ( y- x) * (y-z))   &
        & + z**2 * Log(z) / ( (z - 1._dp) * ( z- y) * (z-x))
   End If
  
 End Function f_5_0


 Real(dp) Function Dyf_5_0(x, y, z)
 Implicit None
  Real(dp), Intent(in) :: x, y, z

   Dyf_5_0 =  (x**2*Log(x))/((-1 + x)*(x - y)**2*(x - z)) +     &
     &  ((y*(-y + z))/((x - y)*(-1 + y)) +                      &
     &     (y*(-y**3 + y*z + x*(y - 2*z + y*z))*Log(y))/        &
     &      ((x - y)**2*(-1 + y)**2) +                          &
     &     (z**2*Log(z))/((-1 + z)*(-x + z)))/(y - z)**2


  
 End Function Dyf_5_0

! dummy version
 Real(dp) Function F_5_1(x, y)
 Implicit None
  Real(dp), Intent(in) :: x, y

  F_5_1 = x+y

 End Function F_5_1

 Real(dp) Function f_6_0(x, y, z)
 !------------------------------------------------------------
 ! loop function taken from C.Bobeth et al., NPB 630 (2002) 87 
 !------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x, y, z

  Real(dp) :: yy, zz

   If ((x.Eq.y).And.(y.Eq.z)) Then
    If (z.Eq.1._dp) Then
     f_6_0 = - 1._dp / 6._dp
    Else
     f_6_0 = (1._dp - z**2 + 2._dp*z*Log(z))/(2._dp *(z-1._dp )**3 *z)
    End If

   Else If ((x.Eq.y).Or.(x.Eq.z).Or.(y.Eq.z)) Then
    If (x.Eq.y) Then
     yy = y
     zz = z
    Else If (x.Eq.z) Then
     yy = z
     zz = y
    Else If (y.Eq.z) Then
     yy = y
     zz = x
    End If

    If (yy.Eq.1._dp) Then
     If (zz.Eq.1._dp) Then
      F_6_0 = -1._dp / 6._dp
     Else
      F_6_0 =(1._dp - zz**2 + 2._dp*zz*Log(zz))/(2._dp*(-1._dp + zz)**3)
     End If

    Else
     if (zz.eq.1._dp) then
      F_6_0 = (-2._dp + 2._dp*yy - (1._dp + yy)*Log(yy))/(-1._dp + yy)**3
     Else
      F_6_0 =(-((yy**2 - zz)*(-1._dp + zz)*Log(yy)) +  &
         &    (-1._dp + yy)*((yy - zz)*(-1._dp + zz) + (-1._dp + yy)*zz*Log(zz)))/ &
         &    ((-1._dp + yy)**2*(yy - zz)**2*(-1._dp + zz)) 
     End If
    End If

   Else If (x.Eq.1._dp) Then
    f_6_0 = 1._dp / ( 1._dp - y - z + y*z)            &
        & + y * Log(y) / ( (y - 1._dp)**2 * (y-z) )   &
        & + z * Log(z) / ( (z - 1._dp)**2 * ( z- y) )
   Else If (y.Eq.1._dp) Then
    f_6_0 = 1._dp / ( 1._dp - x - z + x*z)            &
        & + x * Log(x) / ( (x - 1._dp)**2 * (x-z) )   &
        & + z * Log(z) / ( (z - 1._dp)**2 * ( z- x) )
   Else If (z.Eq.1._dp) Then
    f_6_0 = 1._dp / ( 1._dp - y - x + y*x)            &
        & + y * Log(y) / ( (y - 1._dp)**2 * (y-x) )   &
        & + x * Log(x) / ( (x - 1._dp)**2 * ( x- y) )
   Else
    f_6_0 = x * Log(x) / ( (x - 1._dp) * ( x- y) * (x-z))   &
        & + y * Log(y) / ( (y - 1._dp) * ( y- x) * (y-z))   &
        & + z * Log(z) / ( (z - 1._dp) * ( z- y) * (z-x))
   End If

 End Function f_6_0


 Real(dp) Function Dyf_6_0(x, y, z)
 Implicit None                      ! only full formula is correct
  Real(dp), Intent(in) :: x, y, z

   Dyf_6_0 = (x*Log(x))/((-1 + x)*(x - y)**2*(x - z)) +   &
     &  ((-y + z)/((x - y)*(-1 + y)) +                    &
     &     ((-2*y**3 - x*z + y**2*(1 + x + z))*Log(y))/   &
     &      ((x - y)**2*(-1 + y)**2) +                    &
     &     (z*Log(z))/((-1 + z)*(-x + z)))/(y - z)**2

 End Function Dyf_6_0



 Real(dp) Function f_6_1(x)
 Implicit none
  Real(dp), intent(in) :: x

   f_6_1 = 2._dp * x * (29._dp +  3._dp * x)  / (3._dp * (x - 1._dp)**2 )  &
       & - 2._dp * x * (25._dp + 7._dp * x) * Log(x)                       &
       &             / (3._dp * (x - 1._dp)**3 )                           &
       & - 8._dp * x * Li2(1._dp - 1._dp / x) / (x - 1._dp)**2

 End Function f_6_1

 Real(dp) Function f_7_0(x, y)
 !------------------------------------------------------------
 ! loop function taken from C.Bobeth et al., NPB 630 (2002) 87 
 !------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x, y

!  Integer :: i1
!  Real(dp) :: r

  If (x.Eq.y) Then
   If (y.Eq.0._dp) Then
    f_7_0 = -1._dp
   Else If (y.Eq.1._dp) Then
    f_7_0 = -0.5_dp
   Else
    f_7_0 = (- 1._dp + y - Log(y) ) / (y-1._dp)**2
   End If

  Else If (x.Eq.0._dp) Then
   f_7_0 = 0._dp

  Else If (x.Eq.1._dp) Then
   f_7_0 = (1._dp - y + Log(y)) / (y - 1._dp)**2

  Else If (y.Eq.1._dp) Then
   f_7_0 = x * (1._dp - x + Log(x)) / (x - 1._dp)**2

  Else
   f_7_0 = x * (Log(x) / (x - 1._dp) - Log(y) / (y - 1._dp)) / ( x- y)
  End If
  
 End Function f_7_0


 Real(dp) Function Dxf_7_0(x, y)
 Implicit none
  Real(dp), Intent(in) :: x, y
  
  Dxf_7_0 = (-((x**2 - y)*(-1 + y)*Log(x))                      &
        &   + (-1 + x)*((x - y)*(-1 + y) + (-1 + x)*y*Log(y)))  &
        & /  ((-1 + x)**2*(x - y)**2*(-1 + y))

 End Function Dxf_7_0


 Real(dp) Function f_7_1(x)
 Implicit none
  Real(dp), intent(in) :: x

   f_7_1 = 4._dp * x * (27._dp - 11._dp * x)  / (3._dp * (x - 1._dp)**2 )  &
       & + 4._dp * x * Pi2 / 3._dp                                         &
       & - 4._dp * x * (37._dp - 33._dp * x + 12._dp * x**2) * Log(x)      &
       &             / (3._dp * (x - 1._dp)**3 )                           &
       & - 8._dp * x * (2._dp - 2._dp * x + x**2) * Li2(1._dp - 1._dp / x) &
       &             / (x - 1._dp)**2

 End Function f_7_1

 Real(dp) Function f_8_0(x)
 !------------------------------------------------------------
 ! loop function taken from C.Bobeth et al., NPB 630 (2002) 87 
 !------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  real(dp) :: r
  Real(dp), Parameter :: ci(7) = (/ 1._dp, 0.5_dp, -1./6._dp, 1._dp/12._dp &
                               &  , -1._dp/20._dp, 1._dp/30._dp, -1._dp/42._dp /)
  If (x.Eq.0._dp) Then
   f_8_0 = 0._dp
  Else If (x.Eq.1._dp) Then
   f_8_0 = 1._dp

  Else If (Abs(x).Lt.1.e-6_dp) Then
   f_8_0 = 0._dp
   Do i1=6,1,-1
    f_8_0 = x * (f_8_0 - 1._dp)
   End Do
   f_8_0 = f_8_0 * Log(x)
   
  Else If (Abs(x-1._dp).Lt.1.e-6_dp) Then
   f_8_0 = ci(7)
   r = x - 1._dp 
   Do i1=5,1,-1
    f_8_0 = r * (f_8_0 + ci(i1+1) )
   End Do
   f_8_0 = f_8_0 + ci(1)

  Else
   f_8_0 = x * Log(x) / (x-1._dp)
  End If

 End Function f_8_0


 Real(dp) Function Df_8_0(x)
 Implicit None
  Real(dp), Intent(in) :: x

  Df_8_0 = (-1 + x - Log(x))/(-1 + x)**2

 End Function Df_8_0

 Real(dp) Function f_9_0(w, x, y, z)
 !------------------------------------------------------------
 ! loop function taken from C.Bobeth et al., NPB 630 (2002) 87 
 !------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) ::  w, x, y, z

  If ((z.Eq.x).and.(z.eq.y).and.(z.eq.w)) Then
   If (z.eq.1._dp) then
    F_9_0 = -1._dp / 12._dp
   Else
    F_9_0 = -(2._dp + 3._dp*z - 6._dp *z**2 + z**3 + 6._dp*z*Log(z)) &
        &    / (6._dp*(-1._dp + z)**4*z)
   End If
  Else If (w.Eq.1) Then
   f_9_0 = 1._dp / ( (1._dp-x) * (1._dp-y) * (1._dp-z) )     &
       & + x**2 * Log(x) / ((x-1._dp)**2 * (x-y) * (x-z) )   &
       & + y**2 * Log(y) / ((y-1._dp)**2 * (y-x) * (y-z) )   &
       & + z**2 * Log(z) / ((z-1._dp)**2 * (z-x) * (z-y) )

  Else If (x.Eq.1) Then
   f_9_0 = w**2 * Log(w) / ((w-1._dp)**2 * (w-y) * (w-z) )   &
       & + 1._dp / ( (1._dp-w) * (1._dp-y) * (1._dp-z) )     &
       & + y**2 * Log(y) / ((y-1._dp)**2 * (y-w) * (y-z) )   &
       & + z**2 * Log(z) / ((z-1._dp)**2 * (z-y) * (z-w) )

  Else If (y.Eq.1) Then
   f_9_0 = w**2 * Log(w) / ((w-1._dp)**2 * (w-x) * (w-z) )   &
       & + x**2 * Log(x) / ((x-1._dp)**2 * (x-w) * (x-z) )   &
       & + 1._dp / ( (1._dp-w) * (1._dp-x) * (1._dp-z) )     &
       & + z**2 * Log(z) / ((z-1._dp)**2 * (z-x) * (z-w) )

  Else If (z.Eq.1) Then
   f_9_0 = w**2 * Log(w) / ((w-1._dp)**2 * (w-x) * (w-y) )   &
       & + x**2 * Log(x) / ((x-1._dp)**2 * (x-w) * (x-y) )   &
       & + y**2 * Log(y) / ((y-1._dp)**2 * (y-x) * (y-w) )   &
       & + 1._dp / ( (1._dp-w) * (1._dp-y) * (1._dp-x) )

  Else

   f_9_0 = w**2 * Log(w) / ((w-1._dp) * (w-x) * (w-y) * (w-z) )     &
       & + x**2 * Log(x) / ((x-1._dp) * (x-w) * (x-y) * (x-z) )     &
       & + y**2 * Log(y) / ((y-1._dp) * (y-x) * (y-w) * (y-z) )     &
       & + z**2 * Log(z) / ((z-1._dp) * (z-x) * (z-y) * (z-w) )
  End If

 End Function f_9_0

 Real(dp) Function f_10_0(w, x, y, z)
 !------------------------------------------------------------
 ! loop function taken from C.Bobeth et al., NPB 630 (2002) 87 
 !------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) ::  w, x, y, z


  If ((z.Eq.x).and.(z.eq.y).and.(z.eq.w)) Then
   If (z.eq.1._dp) then
    F_10_0 = 1._dp / 12._dp
   Else
    F_10_0 = (1._dp - 6._dp*z + 3._dp*z**2 + 2._dp*z**3 - 6._dp*z**2*Log(z)) &
         &  / (6._dp*(-1._dp + z)**4*z**2)
   End If
  Else If (w.Eq.1) Then
   f_10_0 = 1._dp / ( (1._dp-x) * (1._dp-y) * (1._dp-z) )  &
        & + x * Log(x) / ((x-1._dp)**2 * (x-y) * (x-z) )   &
        & + y * Log(y) / ((y-1._dp)**2 * (y-x) * (y-z) )   &
        & + z * Log(z) / ((z-1._dp)**2 * (z-x) * (z-y) )

  Else If (x.Eq.1) Then
   f_10_0 = w * Log(w) / ((w-1._dp)**2 * (w-y) * (w-z) )   &
        & + 1._dp / ( (1._dp-w) * (1._dp-y) * (1._dp-z) )  &
        & + y * Log(y) / ((y-1._dp)**2 * (y-w) * (y-z) )   &
        & + z * Log(z) / ((z-1._dp)**2 * (z-y) * (z-w) )

  Else If (y.Eq.1) Then
   f_10_0 = w * Log(w) / ((w-1._dp)**2 * (w-x) * (w-z) )   &
        & + x * Log(x) / ((x-1._dp)**2 * (x-w) * (x-z) )   &
        & + 1._dp / ( (1._dp-w) * (1._dp-x) * (1._dp-z) )  &
        & + z * Log(z) / ((z-1._dp)**2 * (z-x) * (z-w) )

  Else If (z.Eq.1) Then
   f_10_0 = w * Log(w) / ((w-1._dp)**2 * (w-x) * (w-y) )   &
        & + x * Log(x) / ((x-1._dp)**2 * (x-w) * (x-y) )   &
        & + y * Log(y) / ((y-1._dp)**2 * (y-x) * (y-w) )   &
        & + 1._dp / ( (1._dp-w) * (1._dp-y) * (1._dp-x) )

  Else
   f_10_0 = w * Log(w) / ((w-1._dp) * (w-x) * (w-y) * (w-z) )   &
        & + x * Log(x) / ((x-1._dp) * (x-w) * (x-y) * (x-z) )   &
        & + y * Log(y) / ((y-1._dp) * (y-x) * (y-w) * (y-z) )   &
        & + z * Log(z) / ((z-1._dp) * (z-x) * (z-y) * (z-w) )
  End If

 End Function f_10_0

 Real(dp) Function F_10_1(x) 
 Implicit None
  Real(dp), Intent(in) :: x

   F_10_1 = 4._dp * x * ( (19._dp - 3._dp * x) / 3._dp                &
          &             - (17._dp - x) * Log(x) / (3._dp * (x-1._dp)) &
          &             - 2._dp * Li2(1._dp - 1._dp/x)                &
          &             ) / (x-1._dp)**2

 End Function F_10_1


 Real(dp) Function f_11_0(x, y)
 !------------------------------------------------------------
 ! loop function taken from C.Bobeth et al., NPB 630 (2002) 87 
 !------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x, y

!  Integer :: i1
!  Real(dp) :: r
!  Real(dp), Parameter :: ci(7) = (/ 2._dp, 1._dp, -1./3._dp, 1._dp/6._dp &
!                               &  , -1._dp/10._dp, 1._dp/15._dp, -1._dp/21._dp /)

  If (x.Eq.y) Then

   f_11_0 = 1._dp

  Else If (x.Eq.0._dp) Then
   f_11_0 = 0._dp

  Else If (x.Eq.1._dp) Then
   f_11_0 = Log(y) / (1._dp - y )

  Else If (y.Eq.1._dp) Then
   f_11_0 = x * Log(x) / (x - 1._dp)

  Else
   f_11_0 = x * Log(x/y) / (x - y)  
  End If
  
 End Function f_11_0


 Real(dp) Function F1_0(x)
 !-----------------------------------------------------------
 ! loop function as given by Bobeth et al., NPB 630 (2002) 87
 !-----------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: r
  Real(dp), Parameter :: c1(7) = (/ 5._dp/6._dp, -1._dp/24._dp, 0._dp            &
          &                       , 1._dp/120._dp, -1._dp/105._dp, 1._dp/112._dp &
          &                       , -1._dp/126._dp /)

  If (x.Eq.0._dp) Then
   F1_0 = 0._dp

  Else If (x.Eq.1._dp) Then
   F1_0 = 0.75_dp

  Else If (Abs(x-1._dp).Lt.1.e-4_dp) Then
   r = x - 1._dp
   F1_0 = c1(7) * r
   Do i1=1,3
    F1_0 = (c1(7-i1) + F1_0)*r
   End Do
   Do i1=5,6
    F1_0 = (c1(7-i1) + F1_0)*r
   End Do
   F1_0 = 0.75_dp + F1_0

  Else
   F1_0 = - 0.5_dp * x * (6 - x)/(x - 1) &
      &  + 0.5_dp * x * (2 + 3 * x) * Log(x) / (x - 1)**2
  End If

 End Function F1_0


 Real(dp) Function F2_0(x)
 !-----------------------------------------------------------
 ! loop function as given by Bobeth et al., NPB 630 (2002) 87
 !-----------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: r
  Real(dp), Parameter :: c1(7) = (/ -1._dp/6._dp, 1._dp/12._dp, -0.05_dp         &
          &                       , 1._dp/30._dp, -1._dp/42._dp, 1._dp/56._dp &
          &                       , -1._dp/72._dp /)

  If (x.Eq.0._dp) Then
   F2_0 = 0._dp

  Else If (x.Eq.1._dp) Then
   F2_0 = - 0.5_dp

  Else If (Abs(x-1._dp).Lt.1.e-4_dp) Then
   r = x - 1._dp
   F2_0 = c1(7) * r
   Do i1=1,6
    F2_0 = (c1(7-i1) + F2_0)*r
   End Do
   F2_0 = -0.5_dp + F2_0

  Else
   F2_0 = - x / (x - 1) + x * Log(x) / (x - 1)**2
  End If

 End Function F2_0


 Real(dp) Function F1(x)
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp), Parameter :: c1(7) = (/ 1._dp / 24._dp, -0.025_dp, 1._dp/60._dp  &
          &                       , -1._dp / 84._dp , 1._dp/112._dp         &
          &                       , -1._dp / 144._dp , 1._dp / 180._dp /)
  Real(dp) :: r
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "F1"

  If (x.Lt.0) Then
   F1 = 0
   If (ErrorLevel.Gt.-2) Then
    Write(ErrCan,*) "Error in Function F1, x=",x
    If (ErrorLevel.Ge.0) Call TerminateProgram
   End If
  End If

  If (x.Eq.0._dp) Then
   F1 = 1._dp / 6._dp
  Else If (x.Eq.1._dp) Then
   F1 = 1._dp / 24._dp
  Else If (Abs(1._dp - x).Lt.1.e-2_dp) Then
   F1 = c1(1)
   r = x - 1._dp
   Do i1=1,6
    F1 = F1 + c1(i1+1) * r**i1
   End Do 
  Else
   F1 = ( 2._dp + 3._dp * x - 6._dp * x**2 + x**3 &
      & + 6._dp * x * Log(x) ) / (12._dp * (1._dp - x)**4 )
  End If

  Iname = Iname - 1
 End Function F1

 Real(dp) Function F2(x)
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp), Parameter :: c1(7) = (/ 1._dp/24._dp, -1._dp/ 60._dp         &
          &                       , 1._dp/120._dp, -1._dp / 210._dp      &
          &                       , 1._dp/336._dp, -1._dp / 504._dp       &
          &                       , 1._dp/720._dp  /)
  Real(dp) :: r

  If (x.Eq.0._dp) Then
   F2 = 1._dp / 12._dp
  Else If (x.Eq.1._dp) Then
   F2 = 1._dp / 24._dp
  Else If (Abs(1._dp - x).Lt.1.e-2_dp) Then
   F2 = c1(1)
   r = x - 1._dp
   Do i1=1,6
    F2 = F2 + c1(i1+1) * r**i1
   End Do 
  Else
   F2 = ( 1._dp - 6._dp * x + 3._dp * x**2 + 2._dp * x**3 &
      &  - 6._dp * x**2 * Log(x) ) / (12._dp * (1._dp - x)**4 )
  End If

 End Function F2


 Real(dp) Function F3(x)
 Implicit None
  Real(dp), Intent(in) :: x 

  Integer :: i1
  Real(dp) :: r

  If (x.Eq.0._dp) Then
   F3 = Huge(x) ! goes to infinity
  Else If (x.Eq.1._dp) Then
   F3 = 1._dp / 3._dp
  Else If (Abs(1._dp - x).Lt.1.e-2_dp) Then
   F3 = 1._dp / 3._dp
   r = x - 1._dp
   Do i1=1,10
    F3 = F3 + r**i1 * (-1)**i1 / (i1 + 3._dp)
   End Do 
  Else
   F3 = (-3._dp + 4._dp * x - x**2 - 2._dp * Log(x) ) &
    & / (2._dp * (1._dp - x)**3 ) 
  End If

 End Function F3

 Real(dp) Function F3gamma(x)
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: x2, r
  Real(dp), Parameter :: c1(10) = (/ -2._dp/3._dp, 1._dp/ 6._dp        &
          &                       , -1._dp/15._dp, 1._dp / 30._dp      &
          &                       , -2._dp/105._dp, 1._dp / 84._dp     &
          &                       , -1._dp/1260._dp, 1._dp/180._dp     &
          &                       , -2._dp/495._dp, 1._dp / 330._dp  /)
      
  If (x.Eq.1._dp) Then
   F3gamma = - 2._dp / 3._dp
  Else If (x.Eq.0._dp) Then
   F3gamma = - 1._dp
  Else If (Abs(1._dp-x).Lt.1.e-2_dp) Then
   F3gamma = c1(1)
   r = x - 1._dp
   Do i1=1,9
    F3gamma = F3gamma + c1(i1+1) * r**i1
   End Do 

  Else
   x2 = x*x
   F3gamma = (1 - 4*x + 3 * x2 - 2 * x2 * Log(x) ) / (x-1)**3
  End If

 End Function F3gamma

 Real(dp) Function F4(x)
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: r

  If (x.Eq.0._dp) Then
   F4 = 0.5_dp
  Else If (x.Eq.1._dp) Then
   F4 = 1._dp / 6._dp
  Else If (Abs(1._dp - x).Lt.1.e-2_dp) Then
   r = x - 1._dp
   F4 = 0._dp
   Do i1=10,1,-1
    F4 = F4 +  (-r)**i1 / (6._dp + 5._dp * i1 + i1**2)
   End Do 
   F4 = F4 + 1._dp / 6._dp
  Else
   F4 = 0.5_dp * (1+x) / (1-x)**2 + x * Log(x) /(1-x)**3
  End If

 End Function F4


!\section{Function FeynFunctionA}
!\begin{verbatim}
 Real(dp) Function FeynFunctionA(r)
 !----------------------------------------------------------------------
 ! calculation of the Function A(r) ... 0-order in the mass
 ! input:
 !   x ............ m_fermion**2 / m_boson**2
 ! written by Werner Porod, 31.1.01
 !----------------------------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: r

   Real(dp) :: x, x2
   Integer :: i

   If (r.Eq.1._dp) Then
    FeynFunctionA = - 1._dp / 3._dp
  
   Else

    x = r - 1._dp 
    x2 = x * x
    If (x2 .Lt. 1.e-3_dp) Then
     FeynFunctionA = sa(1)
     Do i=1,6
      FeynFunctionA = FeynFunctionA + sa(i+1)*x**i
     End Do
    Else
     FeynFunctionA = ( 1._dp - 0.5_dp * x - Log(r) / x ) / x2
    End If

   End If

 End Function FeynFunctionA
!\end{verbatim}

!\section{Function FeynFunctionB}
!\begin{verbatim}
 Real(dp) Function FeynFunctionB(r)
 !----------------------------------------------------------------------
 ! calculation of the Function B(r) ... 0-order in the mass
 ! input:
 !   x ............ m_fermion**2 / m_boson**2
 ! written by Werner Porod, 31.01.01
 !----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: r

  Real(dp) :: x, x2
  Integer :: i

  If (r.Eq.1._dp) Then
   FeynFunctionB = 1._dp / 6._dp

  Else
   x = r - 1._dp
   x2 = x * x
   If (x2 .Lt. 1.e-3_dp) Then
    FeynFunctionB  = sb(1)
    Do i=1,6
      FeynFunctionB  = FeynFunctionB  + sb(i+1)*x**i
    End Do
   Else
    FeynFunctionB = ( 1._dp + 0.5_dp * x - r * Log(r) / x ) / x2
   End If

  End If

 End Function FeynFunctionB
!\end{verbatim}

 Real(dp) Function FeynFunctionC(r)
 Implicit None
  Real(dp), Intent(in) :: r

  If (r.Eq.1._dp) Then
   FeynFunctionC = - 19._dp / 18._dp
  Else
   FeynFunctionC = 3._dp * FeynFunctionA(r) - FeynFunctionB(r) / 3._dp
  End If
 End Function FeynFunctionC

!\section{Function FeynFunctionD}
!\begin{verbatim}
 Real(dp) Function FeynFunctionD(r)
 !----------------------------------------------------------------------
 ! calculation of the Function D(r) ... 1-order in the mass
 ! input:
 !   x ............ m_fermion**2 / m_boson**2
 ! written by Werner Porod, 31.01.01
 !---------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: r

  Real(dp) :: x, x2
  Integer :: i

  If (r.Eq.1._dp) Then
   FeynFunctionD = - 0.05_dp
  
  Else
   x = r - 1._dp
   x2 = x*x
   If (x2 .Lt. 1.0d-3) Then
    FeynFunctionD = sd(1)
    Do i=1,6
     FeynFunctionD = FeynFunctionD + sd(i+1)*x**i
    End Do
   Else
    FeynFunctionD = ( 4._dp + x - x2/6._dp - (4._dp + 3._dp*x)*Log(r)/x ) &
                  & / (x2*x2)
   End If
  End If

 End Function FeynFunctionD
!\end{verbatim}

!\section{Function FeynFunctionB}
!\begin{verbatim}
 Real(dp) Function FeynFunctionE(r)
 !--------------------------------------------------------------------- 
 ! calculation of the Function E(r) ... 0-order in the mass
 ! input:
 !   x ............ m_fermion**2 / m_boson**2
 ! written by Werner Porod, 31.01.01
 !----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: r

  Real(dp) :: x, x2
  Integer :: i

  If (r.Eq.1._dp) Then
   FeynFunctionE = 1._dp / 30._dp
  
  Else
   x = r - 1._dp 
   x2 = x*x
   If (x2 .Lt. 1.e-3_dp) Then
    FeynFunctionE = se(1)
    Do i=1,6
     FeynFunctionE = FeynFunctionE + se(i+1)*x**i
    End Do
   Else
    FeynFunctionE = ( 4._dp * r + x2 / 3._dp   &
                  & - 2._dp * r * (1._dp + r) * Log(r) / x ) / (x2*x2)
   End If
  End If

 End Function FeynFunctionE

  Real(dp) Function S0low(x)
  Implicit None
   Real(dp), Intent(in) :: x

   S0low = 1._dp - 2.75_dp * x + 0.25_dp * x**2 &
    & - 1.5_dp * x**2 * Log(x) / (1-x)
   S0low = x *  S0low / (1 -x)**2

  End  Function S0low

  Real(dp) Function S0_2(xc, xt)
  !----------------------------------------------
  ! loopfunction for epsilon_K
  ! taken from Buras et al., arXiv:0909.1333
  !----------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: xc, xt

   S0_2 = Log(xt/xc) - 0.75_dp * xt /(1-xt) &
      & - 0.75_dp * xt**2 * Log(xt) / (1-xt)**2
   S0_2 = xc *  S0_2

  End  Function S0_2

 Complex(dp) Function vertexC0tilde(m12in, m22in, m32in)
   Implicit None

   Real(dp), Intent(in) :: m12in, m22in, m32in

   Real(dp) :: m12, m22, m32, r12, r13, r23, x, y, sum, coeff(9), log12, renScale2
   Logical :: check
   Integer :: i1, i2

   Iname = Iname + 1
   NameOfUnit(Iname) = "vertexC0tilde"

   vertexC0tilde = 0._dp
   renScale2 = mudim2

   !-------------------------------------------------------
   ! if all masses are zero, then we have a problem
   ! due to an IR singularity
   !-------------------------------------------------------
   If ((m12in.Eq.0._dp).And.(m22in.Eq.0._dp).And.(m32in.Eq.0._dp)) Then
    vertexC0tilde = Huge(1._dp)
    Iname = Iname - 1
    Return
   End If

   !-------------------------------------------------------
   ! two masses are zero
   !-------------------------------------------------------
   check = .False.
   If ((m12in.Eq.0._dp).And.(m22in.Eq.0._dp)) Then
    m12 = m32in
    check = .True.
   Else If ((m12in.Eq.0._dp).And.(m32in.Eq.0._dp)) Then
    m12 = m22in
    check = .True.
   Else If ((m22in.Eq.0._dp).And.(m32in.Eq.0._dp)) Then
    m12 = m12in
    check = .True.
   End If

   If (check) Then
    vertexC0tilde = 1._dp - Log(m12 /renScale2)
    Iname = Iname - 1
    Return
   End If

   !-------------------------------------------------------
   ! one mass is zero
   !-------------------------------------------------------
   If (m12in.Eq.0._dp) Then
    m12 = m22in
    m22 = m32in
    check = .True.
   Else If (m22in.Eq.0._dp) Then
    m12 = m12in
    m22 = m32in
    check = .True.
   Else If (m32in.Eq.0._dp) Then
    m12 = m12in
    m22 = m22in
    check = .True.
   End If

   If (check) Then
    r12 = (m12-m22) / m12
    
    If (r12.Eq.0._dp) Then
     vertexC0tilde = - Log(m12/renScale2) 
     Iname = Iname - 1
     Return

    Else If (Abs(r12).Lt.1.e-4_dp) Then
     vertexC0tilde = r12 / 56._dp
     Do i1=6,1,-1
      vertexC0tilde = ( vertexC0tilde + 1._dp / Real(i1*(i1+1),dp) ) * r12
     End Do
     vertexC0tilde = vertexC0tilde - Log(m12/renScale2)
     Iname = Iname - 1
     Return

    Else
     vertexC0tilde = 1._dp + (m22*Log(m22/renScale2) - m12*Log(m12/renScale2)) &
                    &         / (m12 - m22)
     Iname = Iname - 1
     Return
    End If
   End If

   !-------------------------------------------------------
   ! two masses are equal
   !-------------------------------------------------------
   r12 = Abs(m12in-m22in)/m12in
   r13 = Abs(m12in-m32in)/m12in
   r23 = Abs(m22in-m32in)/m22in
   
   If ((r12.Eq.0._dp).And.(r13.Eq.0._dp)) Then ! all masses equal
    vertexC0tilde = - Log(m12in/renScale2) 
    Iname = Iname - 1
    Return

   Else If (r12.Eq.0._dp) Then ! m1=m2
    m12 = m32in
    m22 = m12in
    check = .True.

   Else If (r13.Eq.0._dp) Then ! m1=m3
    m12 = m22in
    m22 = m12in
    check = .True.
  
   Else If (r23.Eq.0._dp) Then ! m2=m3
    m12 = m12in
    m22 = m22in
    check = .True.  
   End If 
! m22 die Masse, die 2x vorkommt
! m12 die Masse, die 1x vorkommt
   If (check) Then ! two masses are equal
    x = (m12-m22)/m12
    y = m22 / m12 
    If (Abs(x).Lt.1.e-4_dp) Then
     vertexC0tilde = 2._dp * x / 63._dp
     Do i1=6,1,-1
      vertexC0tilde = (vertexC0tilde + 2._dp/Real(i1*(i1+2),dp))*x
     End Do
     vertexC0tilde = vertexC0tilde - Log(m12/renScale2) - 0.5_dp

     ! MEK: corrected code 05.07.2013
    Else If ( y .Lt.1.e-4_dp) then ! third mass is much larger than the other two
     Log12 = Log(y)
     vertexC0tilde = (1._dp+9*Log12)*y
     Do i1=7,1,-1
       vertexC0tilde = (vertexC0tilde + (1._dp+(i1+1)*Log12) ) * y
     End Do
     vertexC0tilde=vertexC0tilde - Log(m12/renScale2) + 1._dp
     
    ! MEK: corrected code 05.07.2013
    Else If ( y .gt.1.e4_dp) then ! third mass is much smaller than the other two
     y = 1._dp / y  ! taking the inverse to expand
     Log12 = Log(y)
     vertexC0tilde = - (1._dp+7*Log12)*y
     Do i1=7,1,-1
      vertexC0tilde = (vertexC0tilde - (1._dp+(i1-1)*Log12) ) * y
     End Do 
     vertexC0tilde=vertexC0tilde-Log(m22/renScale2)    

    Else
     vertexC0tilde= (m12 * (m12 -m22 ) -m12**2 * Log(m12 /renScale2)      & 

        & +(2._dp*m12 -m22 )*m22 *Log(m22 /renScale2))/(m12 -m22 )**2
    End If

    Iname = Iname - 1
    Return
   End If  ! two masses are equal

   !-------------------------------------------------------
   ! all masses are nearly equal
   !-------------------------------------------------------
   If ((r12.Lt.1.e-4_dp).And.(r13.Lt.1.e-4_dp)) Then
    x = (m12in-m22in)/m12in
    y = (m12in-m32in)/m12in
    ! check which one has the larger modulus
    If (Abs(x).Lt.Abs(y)) Then ! swap order to simplify things below
     sum = x
     x = y
     y = sum
    End If
    sum = y**9
    Do i1=1,8
     sum = sum + x**i1 * y**(9-i1)
    End Do 
    sum = sum + x**9
    vertexC0tilde = sum / 495._dp
    Do i1=8,1,-1
     sum = y**i1
     Do i2=1,i1-1
      sum = sum + x**i2 * y**(i1-i2)
     End Do 
     sum = sum + x**i1
     vertexC0tilde = vertexC0tilde + 2._dp * sum / Real(i1*(i1+1)*(i1+2),dp)
    End Do
    vertexC0tilde = vertexC0tilde - Log(m12in/renScale2) - 0.5_dp
    Iname = Iname - 1
    Return
   End If

   !-------------------------------------------------------
   ! the remaining cases
   !-------------------------------------------------------
   If (r12.Lt.1.e-4_dp) Then      ! m1 nearly m2
    m12 = m32in
    m22 = m12in
    x = (m12in-m22in) / m12in
    check = .True.
   Else If (r13.Lt.1.e-4_dp) Then ! m1 nearly m3
    m12 = m22in
    m22 = m12in
    x = (m12in-m32in) / m12in
    check = .True.
   Else If (r23.Lt.1.e-4_dp) Then ! m2 nearly m3
    m12 = m12in
    m22 = m22in
    x = (m22in-m32in) / m22in
    check = .True.
   Else                           ! all masses are significantly different
    m12 = m12in
    m22 = m22in
    m32 = m32in
   End If

   If (check) Then ! two masses are nearly equal
    log12 = Log(m12/m22)
 
    coeff(1) = -(m22*((3 - 2*Log12)*m12**2 - 4*m12*m22 + m22**2))/2._dp
    coeff(2) = ( m22*(2*m12**3 + (3 - 6*Log12)*m12**2*m22 - 6*m12*m22**2 &
             & + m22**3) )/6._dp
    coeff(3) = ( m22*(m12**4 - 8*m12**3*m22 + 12*Log12*m12**2*m22**2     &
             & + 8*m12*m22**3 - m22**4))/12._dp
    coeff(4) = ( m22*(2*m12**5 - 15*m12**4*m22 + 60*m12**3*m22**2        & 
             & - 20*(1 + 3*Log12)*m12**2*m22**3 - 30*m12*m22**4          & 

             & + 3*m22**5) ) / 60._dp
    coeff(5) = ( m22*(m12**6 - 8*m12**5*m22 + 30*m12**4*m22**2           &
             & - 80*m12**3*m22**3 + 5*(7 + 12*Log12)*m12**2*m22**4       & 

             & + 24*m12*m22**5 - 2*m22**6) ) / 60._dp
    coeff(6) = ( m22*(4*m12**7 - 35*m12**6*m22 + 140*m12**5*m22**2       & 
             & - 350*m12**4*m22**3 + 700*m12**3*m22**4                   &
             & - 7*(47 + 60*Log12)*m12**2*m22**5 - 140*m12*m22**6        &
             & + 10*m22**7) ) / 420._dp
    coeff(7) = ( m22*(5*m12**8 - 48*m12**7*m22 + 210*m12**6*m22**2       &
             & - 560*m12**5*m22**3 + 1050*m12**4*m22**4                  &
             & - 1680*m12**3*m22**5                                      &
             & + 42*(19 + 20*Log12)*m12**2*m22**6 + 240*m12*m22**7       & 

             & - 15*m22**8) ) / 840._dp
    coeff(8) = ( m22*(10*m12**9 - 105*m12**8*m22 + 504*m12**7*m22**2     &
             & - 1470*m12**6*m22**3 + 2940*m12**5*m22**4                 &
             & - 4410*m12**4*m22**5 + 5880*m12**3*m22**6                 &
             & - 18*(153 + 140*Log12)*m12**2*m22**7 - 630*m12*m22**8     &
             & + 35*m22**9) ) / 2520._dp
    coeff(9) = ( m22*(7*m12**10 - 80*m12**9*m22 + 420*m12**8*m22**2      &
             & - 1344*m12**7*m22**3 + 2940*m12**6*m22**4                 &
             & - 4704*m12**5*m22**5 + 5880*m12**4*m22**6                 &
             & - 6720*m12**3*m22**7                                      &
             & + 9*(341 + 280*Log12)*m12**2*m22**8 + 560*m12*m22**9      &
             & - 28*m22**10) ) / 2520._dp
    x = x / (m12-m22)
    vertexC0tilde = coeff(9) * x
    Do i1=8,1,-1
     vertexC0tilde = (vertexC0tilde + coeff(i1) ) * x 
    End Do
    vertexC0tilde =                                                    &
        &  ( vertexC0tilde + m12 *(m12 -m22 )-m12**2 * Log(m12 /renScale2) & 
        &   +(2._dp*m12 -m22 )*m22 *Log(m22 /renScale2)) / (m12 -m22 )**2   


   Else

    vertexC0tilde = 1._dp-1._dp/(m22 -m32 )*((m12 **2*Log(m12 /renScale2) & 
        & -m22 **2*Log(m22 /renScale2))/(m12 -m22 )                    &
        & -(m12 **2*Log(m12 /renScale2) -m32 **2*Log(m32 /renScale2))     &
        &  /(m12 -m32 ))

   End If
   
!    vertexC0tilde = vertexC0tilde + divergence

   Iname = Iname - 1

 End Function vertexC0tilde

 Complex(dp) Function vertexC00(masa1c,masa2c,masa3c)
   Real(dp), Intent(in) :: masa1c,masa2c,masa3c

   vertexC00 = vertexC0tilde(masa1c,masa2c,masa3c) + 0.5_dp 

 End Function vertexC00


 Complex(dp) Function vertexC11(masa1c,masa2c,masa3c)
   Real(dp), Intent(in) :: masa1c,masa2c,masa3c

   vertexC11 = 0._dp

   If (masa2c.Eq.masa3c) Then

   vertexC11=-(3._dp*masa1c**2-4._dp*masa1c*masa2c+masa2c**2 & 
        & -2._dp*masa1c**2*Log(masa1c/masa2c)) & 
        & /(2._dp*(masa1c-masa2c)**3)   

   Else

   vertexC11=masa1c*(masa1c**2/masa2c**2-masa1c/masa2c & 
        & -masa1c**2/masa3c**2+masa1c**3/(masa2c*masa3c**2) & 
        & +masa1c/masa3c-masa1c**3/(masa2c**2*masa3c) & 
        & +masa1c**2/masa3c**2*Log(masa1c/masa2c) & 
        & -2._dp*masa1c**3/(masa2c*masa3c**2)*Log(masa1c/masa2c) & 
        & -masa1c**2/masa2c**2*Log(masa1c/masa3c) & 
        & +2._dp*masa1c**3/(masa2c**2*masa3c)*Log(masa1c/masa3c) & 
        & -Log(masa2c/masa3c)+2._dp*masa1c/masa2c*Log(masa2c/masa3c) & 
        & +2._dp*masa1c/masa3c*Log(masa2c/masa3c) & 
        & -4._dp*masa1c**2/(masa2c*masa3c)*Log(masa2c/masa3c)) & 
        & /(2._dp*(1._dp-masa1c/masa2c)**2*masa2c  & 
        & *(1._dp-masa1c/masa3c)**2*(masa1c/masa2c-masa1c/masa3c) & 
        & *masa3c)

   End If

 End Function vertexC11


 Complex(dp) Function vertexC12(masa1c,masa2c,masa3c)
   Real(dp), Intent(in) :: masa1c,masa2c,masa3c

   vertexC12 = 0._dp

   If (masa2c.Eq.masa3c) Then

   vertexC12=-(3._dp*masa1c**2-4._dp*masa1c*masa2c+masa2c**2 & 
        & -2._dp*masa1c**2*Log(masa1c/masa2c)) & 
        & /(4._dp*(masa1c-masa2c)**3)   

   Else

   vertexC12=((masa1c-masa2c)*(masa1c-masa3c)*(masa2c-masa3c)*masa3c      & 
        & +masa1c**2*(masa2c**2*Log(masa1c/masa2c)+masa3c                 & 
        & *(-2._dp*masa2c+masa3c)*Log(masa1c/masa3c))                     & 
        & +masa2c**2*(2._dp*masa1c-masa3c)*masa3c*Log(masa2c/masa3c))     &
        & /(2._dp*(masa1c-masa2c)*(masa1c-masa3c)**2*(masa2c-masa3c)**2)

   End If

 End Function vertexC12
Logical Function SmallDifference(m1,m2)
Implicit None
Real(dp), Intent(in) :: m1, m2
! avelino
! changed eps=1E-8_dp -> eps=1E-6_dp
Real(dp) :: eps=1E-6_dp

If (Abs(m1).lt.eps) Then
 If (Abs(m2).lt.eps) Then
   SmallDifference = .true.
 Else 
   SmallDifference = .False.
 End if 
Else if (Abs(m2).lt.eps) Then
  SmallDifference = .false.
Else 
 If ((Abs(m1-m2)/Max(m1,m2)).lt.eps) Then 
   SmallDifference = .true.
 Else
   SmallDifference = .False.
 End if
End if


End Function SmallDifference

Real(dp) Function C0m(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: eps=1E-10_dp, large = 1E+5_dp

! !    C0m  = -C0_3m(m1,m2,m3)  ! check sign

If ((SmallDifference(m1,m2)).And.(SmallDifference(m1,m3))) Then ! all masse equal
  If (Abs(m1).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C0m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C0m = large
  Else
   C0m = 1/(2.*m1)
  End if

Else If (SmallDifference(m1,m2)) Then ! m1 = m2
 If (Abs(m3).lt.eps) Then
   C0m = 1/m2
 Else If (Abs(m2).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C0m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C0m = large/m3
 Else
   C0m = (m2 - m3 + m3*Log(m3/m2))/(m2 - m3)**2
 End if

Else If (SmallDifference(m1,m3)) Then ! m1 = m3
 If (Abs(m2).lt.eps) Then
   C0m = 1/m3
 Else If (Abs(m3).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C0m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C0m = large/m2
 Else
   C0m =   (-m2 + m3 + m2*Log(m2/m3))/(m2 - m3)**2
 End if

Else If (SmallDifference(m3,m2)) Then ! m2 = m3
 If (Abs(m1).lt.eps) Then
   C0m = 1/m3
 Else If (Abs(m3).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C0m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C0m = large/m1
 Else
   C0m =  (-m1 + m3 - m1*Log(m3/m1))/(m1 - m3)**2
 End if

Else
 If (Abs(m1).lt.eps) Then
  C0m = -Log(m3/m2)/(m2 - m3)
 Else if (Abs(m2).lt.eps) Then
  C0m = -Log(m3/m1)/(m1 - m3)
 Else if (Abs(m3).lt.eps) Then
  C0m = -Log(m2/m1)/(m1 - m2)
 Else
  C0m =  (m2*(-m1 + m3)*Log(m2/m1) + (m1 - m2)*m3*Log(m3/m1))/ &
     &  ((m1 - m2)*(m1 - m3)*(m2 - m3))
 End if

End if

 C0m = - C0m

End Function C0m

Real(dp) Function C00m(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: eps=1E-10_dp, large = 1E+5_dp

If ((SmallDifference(m1,m2)).And.(SmallDifference(m1,m3))) Then ! all masse equal
 If (Abs(m3).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C00m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C00m = large
 Else
   C00m = -Log(m3)/4.
 End if

Else If (SmallDifference(m1,m2)) Then ! m1 = m2
  If (Abs(m3).lt.eps) Then
   C00m = 0.125_dp - Log(m2)/4._dp
  Else if (Abs(m2).lt.eps) Then
   C00m = 0.375_dp - Log(m3)/4._dp
  Else
   C00m = ((m2 - m3)*(m2 - 3*m3 + 2*(-m2 + m3)*Log(m2)) - &
   &  2*m3**2*Log(m3/m2))/(8.*(m2 - m3)**2)
  End if

Else If (SmallDifference(m1,m3)) Then ! m1 = m3
  If (Abs(m2).lt.eps) Then
   C00m = 0.125_dp - Log(m2)/4._dp
  Else if (Abs(m3).lt.eps) Then
   C00m = 0.375_dp - Log(m3)/4._dp
  Else
   C00m =    -(2*m2**2*Log(m2/m3) +  &
    &      (m2 - m3)*(-3*m2 + m3 + 2*(m2 - m3)*Log(m3)))/ &
    &   (8.*(m2 - m3)**2)
  End if

Else If (SmallDifference(m3,m2)) Then ! m2 = m3
  If (Abs(m1).lt.eps) Then
   C00m = 0.125_dp - Log(m2)/4._dp
  Else if (Abs(m3).lt.eps) Then
   C00m = 0.375_dp - Log(m3)/4._dp
  Else
   C00m =  (-((m1 - m3)*(-3*m1 + m3 + 2*(m1 - m3)*Log(m1))) + &
   &    2*(2*m1 - m3)*m3*Log(m3/m1))/(8.*(m1 - m3)**2)
  End if

Else
 If (Abs(m1).lt.eps) Then
  C00m = -(-3*m2 + 3*m3 + 2*m2*Log(m2) - 2*m3*Log(m3))/(8.*(m2 - m3))
 Else If (Abs(m2).lt.eps) Then
  C00m = -(-3*m1 + 3*m3 + 2*m1*Log(m1) - 2*m3*Log(m3))/(8.*(m1 - m3))
 Else If (Abs(m3).lt.eps) Then
  C00m = -(-3*m1 + 3*m2 + 2*m1*Log(m1) - 2*m2*Log(m2))/(8.*(m1 - m2))
 Else 
   C00m =  ((-m1 + m3)*((m1 - m2)*(m2 - m3)*(-3 + 2*Log(m1)) -  &
  &          2*m2**2*Log(m2/m1)) + 2*(-m1 + m2)*m3**2*Log(m3/m1))/ &
  &      (8.*(m1 - m2)*(m1 - m3)*(m2 - m3))
  End if

End if

  C00m = - C00m 

End Function C00m



Real(dp) Function C11m(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: eps=1E-10_dp, large = 1E+5_dp

If ((SmallDifference(m1,m2)).And.(SmallDifference(m1,m3))) Then ! all masse equal
  If (Abs(m1).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C11m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C11m = large
  Else
   C11m = 1/(12.*m1)
  End if

Else If (SmallDifference(m1,m2)) Then ! m1 = m2
  If (Abs(m1).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C11m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C11m = -large
  Else if (Abs(m3).lt.eps) Then
   C11m = -1/(9.*m2)
  Else
   C11m = ((m2 - m3)*(2*m2**2 - 7*m2*m3 + 11*m3**2) + &
   &      6*m3**3*Log(m3/m2))/(18.*(m2 - m3)**4)
  End if

Else If (SmallDifference(m1,m3)) Then ! m1 = m3
 If (Abs(m1).lt.eps) Then
   C11m = -1/(6.*m2)
  Else if (Abs(m2).lt.eps) Then
   C11m =-1/(3.*m3)
  Else
   C11m = (m2**3 - 6*m2**2*m3 + 3*m2*m3**2 + 2*m3**3 + &
   &      6*m2*m3**2*Log(m2/m3))/(6.*(m2 - m3)**4)
  End if

Else If (SmallDifference(m3,m2)) Then ! m2 = m3
 If (Abs(m1).lt.eps) Then
   C11m = -1/(9.*m3)
  Else if (Abs(m2).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C11m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C11m =-large
  Else
   C11m =  -((m1 - m3)*(11*m1**2 - 7*m1*m3 + 2*m3**2) + &
   &       6*m1**3*Log(m3/m1))/(18.*(m1 - m3)**4)
  End if

Else
 If (Abs(m1).lt.eps) Then
  C11m = (-((m2 - 3*m3)*(m2 - m3)) + 2*m3**2*Log(m3/m2))/(6.*(m2 - m3)**3)
 Else If (Abs(m2).lt.eps) Then
  C11m =Log(m3/m1)/(3*m1 - 3*m3)
 Else If (Abs(m3).lt.eps) Then
  C11m = (3*m1**2 - 4*m1*m2 + m2**2 + 2*m1**2*Log(m2/m1))/(6.*(m1 - m2)**3)
 Else 
   C11m =   (m2*(m1 - m3)*(-((-m1 + m2)*(m2 - m3)* &
     &          (-3*m1*m2 + m2**2 + 5*m1*m3 - 3*m2*m3)) -  &
     &       2*(m1*m2*(m2 - 3*m3)*m3 + m2**2*m3**2 +  &
     &          m1**2*(m2**2 - 3*m2*m3 + 3*m3**2))*Log(m2/m1)) +  &
     &    2*(m1 - m2)**3*m3**3*Log(m3/m1))/ &
     &  (6.*(m1 - m2)**3*(m1 - m3)*(m2 - m3)**3)
 End if
End if

 C11m = - C11m

End Function C11m


Real(dp) Function C12m(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: eps=1E-10_dp, large = 1E+5_dp

If ((SmallDifference(m1,m2)).And.(SmallDifference(m1,m3))) Then ! all masse equal
  If (Abs(m1).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C12m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C12m = large
  Else
   C12m = 1/(24.*m1)
  End if

Else If (SmallDifference(m1,m2)) Then ! m1 = m2
 If (Abs(m1).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C12m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C12m = -large
  Else if (Abs(m3).lt.eps) Then
   C12m = -1/(9.*m2)
  Else
   C12m = (m2**3 - 6*m2**2*m3 + 3*m2*m3**2 + 2*m3**3 &
    &  - 6*m2*m3**2*Log(m3/m2))/(12.*(m2 - m3)**4)
  End if

Else If (SmallDifference(m1,m3)) Then ! m1 = m3
 If (Abs(m2).lt.eps) Then
   C12m = -1/(3.*m3)
  Else if (Abs(m3).lt.eps) Then
   C12m = -1/(6.*m2)
  Else
   C12m = (2*m2**3 + 3*m2**2*m3 - 6*m2*m3**2 + m3**3 &
   & - 6*m2**2*m3*Log(m2/m3))/(12.*(m2 - m3)**4)
  End if

Else If (SmallDifference(m3,m2)) Then ! m2 = m3
 If (Abs(m1).lt.eps) Then
   C12m = -1/(9.*m3)
  Else if (Abs(m3).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C12m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C12m = -large
  Else
   C12m = -((m1 - m3)*(11*m1**2 - 7*m1*m3 + 2*m3**2) &
   &  + 6*m1**3*Log(m3/m1))/(36.*(m1 - m3)**4)
  End if

Else
 If (Abs(m1).lt.eps) Then
  C12m = (m2**2 - m3**2 + 2*m2*m3*Log(m3/m2))/(6.*(m2 - m3)**3)
 Else If (Abs(m2).lt.eps) Then
  C12m =(-m1 + m3 + m1*Log(m1/m3))/(6.*(m1 - m3)**2)
 Else If (Abs(m3).lt.eps) Then
  C12m = (-m1 + m2 + m1*Log(m1/m2))/(6.*(m1 - m2)**2)
 Else 
   C12m =   -(m2**2*(m1 - m3)**2*(m1*(m2 - 3*m3) + 2*m2*m3)* &
     &      Log(m2/m1) + (m1 - m2)* &
     &      ((m1 - m3)*(m2 - m3)*  &
     &         (-(m2*m3*(m2 + m3)) + m1*(m2**2 + m3**2)) +  &
     &        (m1 - m2)*m3**2*(3*m1*m2 - (m1 + 2*m2)*m3)*Log(m3/m1)) &
     &     )/(6.*(m1 - m2)**2*(m1 - m3)**2*(m2 - m3)**3)
 End if

End if

 C12m = - C12m

End Function C12m



Real(dp) Function C22m(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: r1, r2
Real(dp) :: eps=1E-10_dp, large = 1E+5_dp

r1 = m2/m1
r2 = m3/m1

If ((SmallDifference(m1,m2)).And.(SmallDifference(m1,m3))) Then ! all masse equal
  If (Abs(m1).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C22m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C22m = -large 
  Else
   C22m = -1/(12.*m1)
  End if

Else If (SmallDifference(m1,m2)) Then ! m1 = m2
 If (Abs(m1).lt.eps) Then
   C22m = 1/(6.*m3)
  Else if (Abs(m3).lt.eps) Then
   C22m = 1/(3.*m2)
  Else
    C22m=-(2*m2**3 + 3*m2**2*m3 - 6*m2*m3**2 + m3**3 -  &
    &       6*m2**2*m3*Log(m2/m3))/(6.*(m2 - m3)**4)
  Endif

Else If (SmallDifference(m1,m3)) Then ! m1 = m3
 If (Abs(m2).lt.eps) Then
   C22m = 1/(9.*m3)
  Else if (Abs(m3).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C22m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C22m = large
  Else
     C22m=((m2 - m3)*(11*m2**2 - 7*m2*m3 + 2*m3**2) -  &
    &      6*m2**3*Log(m2/m3))/(18.*(m2 - m3)**4)
  End if

Else If (SmallDifference(m3,m2)) Then ! m2 = m3
 If (Abs(m1).lt.eps) Then
   C22m = 1/(9.*m3)
  Else if (Abs(m3).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C22m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C22m = large
  Else
     C22m =  ((m1 - m3)*(11*m1**2 - 7*m1*m3 + 2*m3**2) +  &
    &     6*m1**3*Log(m3/m1))/(18.*(m1 - m3)**4)  
  End if

Else
 If (Abs(m1).lt.eps) Then
  C22m = -(3*m2**2 - 4*m2*m3 + m3**2 - 2*m2**2*Log(m2/m3))/(6.*(m2 - m3)**3)
 Else If (Abs(m2).lt.eps) Then
  C22m = -(3*m1**2 - 4*m1*m3 + m3**2 - 2*m1**2*Log(m1/m3))/(6.*(m1 - m3)**3)
 Else If (Abs(m3).lt.eps) Then
  C22m = Log(m1/m2)/(3*m1 - 3*m2)
 Else
   C22m =   (2*m1**3*(m2 - m3)**3*Log(m2/m1) + (m1 - m2)*m3* &
     &     (-((m1 - m3)*(m2 - m3)*(5*m1*m2 - 3*(m1 + m2)*m3 + m3**2)) + &
     &       2*(3*m1**2*m2**2 - 3*m1*m2*(m1 + m2)*m3 + &
     & (m1**2 + m1*m2 + m2**2)*m3**2)*Log(m2/m3)))/  &
     &  (6.*(m1 - m2)*(m1 - m3)**3*(m2 - m3)**3)
 End if

End if

! ! C22m = - C22m ! check sign!
End Function C22m



Real(dp) Function C2m(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: eps=1E-10_dp, large = 1E+5_dp

If ((SmallDifference(m1,m2)).And.(SmallDifference(m1,m3))) Then ! all masse equal
   C2m = -1/(6.*m1)

Else If (SmallDifference(m1,m2)) Then ! m1 = m2
 If (Abs(m1).lt.eps) Then
   C2m =-1/(2.*m3)
  Else if (Abs(m3).lt.eps) Then
   C2m = -1/(2.*m2)
  Else
   C2m = (-m2**2 + m3**2 + 2*m2*m3*Log(m2/m3))/(2.*(m2 - m3)**3)
  End if

Else If (SmallDifference(m1,m3)) Then ! m1 = m3
 If (Abs(m2).lt.eps) Then
   C2m = -1/(4.*m3)
  Else if (Abs(m3).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C2m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C2m = -large
  Else
   C2m = (3*m2**2 - 4*m2*m3 + m3**2 - 2*m2**2*Log(m2/m3))/(4.*(m2 - m3)**3)
  End if

Else If (SmallDifference(m3,m2)) Then ! m2 = m3
 If (Abs(m1).lt.eps) Then
   C2m = -1/(4.*m3)
  Else if (Abs(m3).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C2m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C2m = -large
  Else
   C2m = (3*m1**2 - 4*m1*m3 + m3**2 + 2*m1**2*Log(m3/m1))/(4.*(m1 - m3)**3)
  End if

Else

 If (Abs(m1).lt.eps) Then
  C2m = (m2 - m3 + m2*Log(m3/m2))/(2.*(m2 - m3)**2)
 Else If (Abs(m2).lt.eps) Then
  C2m = (m1 - m3 + m1*Log(m3/m1))/(2.*(m1 - m3)**2)
 Else If (Abs(m3).lt.eps) Then
  C2m = Log(m2/m1)/(2*m1 - 2*m2)
 Else
   C2m =  (m1**2*(m2 - m3)**2*Log(m2/m1) +   &
     &    (-m1 + m2)*m3*((m1 - m3)*(m2 - m3) +  & 
     &       (-2*m1*m2 + (m1 + m2)*m3)*Log(m2/m3)))/ &
     &  (2.*(m1 - m2)*(m1 - m3)**2*(m2 - m3)**2) 
  End if
End if

 C2m = - C2m

End Function C2m



Real(dp) Function C1m(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: eps=1E-10_dp, large = 1E+5_dp

If ((SmallDifference(m1,m2)).And.(SmallDifference(m1,m3))) Then ! all masse equal
  If (Abs(m1).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C12m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C1m = -large 
  Else
   C1m = -1/(6.*m1)
  End if

Else If (SmallDifference(m1,m2)) Then ! m1 = m2
 If (Abs(m2).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C12m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C1m = -large
  Else if (Abs(m3).lt.eps) Then
   C1m = -1/(4.*m2)
  Else
   C1m = -((m2 - 3*m3)*(m2 - m3) - 2*m3**2*Log(m3/m2))/(4.*(m2 - m3)**3)
  End if

Else If (SmallDifference(m1,m3)) Then ! m1 = m3
 If (Abs(m2).lt.eps) Then
   C1m = -1/(2.*m3)
  Else if (Abs(m3).lt.eps) Then
   C1m = -1/(2.*m2)
  Else
   C1m = (-m2**2 + m3**2 + 2*m2*m3*Log(m2/m3))/(2.*(m2 - m3)**3)
  End if

Else If (SmallDifference(m3,m2)) Then ! m2 = m3
 If (Abs(m2).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in C12m"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   C1m = -large
  Else if (Abs(m1).lt.eps) Then
   C1m = -1/(4.*m3)
  Else
   C1m = (3*m1**2 - 4*m1*m3 + m3**2 + 2*m1**2*Log(m3/m1))/(4.*(m1 - m3)**3)
  End if

Else
 If (Abs(m1).lt.eps) Then
  C1m = (-m2 + m3 + m3*Log(m2/m3))/(2.*(m2 - m3)**2)
 Else If (Abs(m2).lt.eps) Then
  C1m = Log(m3/m1)/(2*m1 - 2*m3)
 Else If (Abs(m3).lt.eps) Then
  C1m =(m1 - m2 + m1*Log(m2/m1))/(2.*(m1 - m2)**2)
 Else
   C1m =   (m2*(-m1 + m3)*((-m1 + m2)*(m2 - m3) -      &
     &       (m1*(m2 - 2*m3) + m2*m3)*Log(m2/m1)) +    &
     &    (m1 - m2)**2*m3**2*Log(m3/m1))/              &
     &  (2.*(m1 - m2)**2*(m1 - m3)*(m2 - m3)**2) 
 End if
End if

 C1m = - C1m

End Function C1m


Real(dp) Function MonopoleFSS(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: eps=1E-10_dp

If ((SmallDifference(m2,m3)).and.(SmallDifference(m1,m2))) Then
 If (Abs(m1).lt. eps) Then
    Write(ErrCan,*) "Numerical problem in MonopoleFSS"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
    MonopoleFSS = -1/eps
 Else
    MonopoleFSS = -1/(24.*m1)
 End if

Else If (SmallDifference(m2,m3)) Then
  If (Abs(m1).lt.eps) Then
    MonopoleFSS = -1/(18.*m3)
  Else If (Abs(m2).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in MonopoleFSS"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
    MonopoleFSS = -1/eps
  Else 
    MonopoleFSS = (11*m1**3 - 18*m1**2*m3 + 9*m1*m3**2 - 2*m3**3 & 
         & - 6*m1**3*Log(m1) + 6*m1**3*Log(m3))/(36.*(m1 - m3)**4)
  End if

Else If (SmallDifference(m1,m3)) Then
  If (Abs(m3).lt.eps) Then
    MonopoleFSS = -1/(6.*m2)
  Else If (Abs(m2).lt.eps) Then
    MonopoleFSS = -1/(12.*m3)
  Else 
    MonopoleFSS = -(2*m2**3 + 3*m2**2*m3 - 6*m2*m3**2 + m3**3 - 6*m2**2*m3*Log(m2) &
         &  + 6*m2**2*m3*Log(m3))/(12.*(m2 - m3)**4)
  End if

Else If (SmallDifference(m1,m2)) Then
  If (Abs(m3).lt.eps) Then
    MonopoleFSS = -1/(12.*m2)
  Else If (Abs(m2).lt.eps) Then
    MonopoleFSS = -1/(6.*m3)
  Else 
    MonopoleFSS = -(m2**3 - 6*m2**2*m3 + 3*m2*m3**2 + 2*m3**3 +  6*m2*m3**2*Log(m2) &
         & - 6*m2*m3**2*Log(m3))/(12.*(m2 - m3)**4)
  End if

Else 

 If (Abs(m1).lt.eps) Then
   MonopoleFSS = (-m2**2 + m3**2 + 2*m2*m3*Log(m2/m3))/(6.*(m2 - m3)**3)  
 Else If (Abs(m2).lt.eps) Then
   MonopoleFSS = (m1 - m3 + m1*Log(m3/m1))/(6.*(m1 - m3)**2)
 Else If (Abs(m3).lt.eps) Then
   MonopoleFSS = (m1 - m2 + m1*Log(m2/m1))/(6.*(m1 - m2)**2)
 Else
  MonopoleFSS = (-(m1**3*(m2 - m3)**3*Log(m1)) &
     &     + (m1 - m3)*((m1 - m2)*(m2 - m3)*(-(m2*m3*(m2 + m3)) &
     &  + m1*(m2**2 + m3**2)) +  &
     &       m2**2*(m1 - m3)*(m1*(m2 - 3*m3) + 2*m2*m3)*Log(m2)) &
     &  + (m1 - m2)**2*m3**2*(3*m1*m2 - (m1 + 2*m2)*m3)*Log(m3))/   &
     &  (6.*(m1 - m2)**2*(m1 - m3)**2*(m2 - m3)**3)
  End if
End If

End Function MonopoleFSS


Real(dp) Function MonopoleSFF(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: eps=1E-10_dp

If ((SmallDifference(m2,m3)).and.(SmallDifference(m1,m2))) Then
 If (Abs(m1).lt.eps) Then 
    Write(ErrCan,*) "Numerical problem in MonopoleSFF"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
  MonopoleSFF = -1/eps
 Else
  MonopoleSFF = -1/(8.*m1)
 End if

Else If (SmallDifference(m2,m3)) Then
  If (Abs(m2).lt.eps) Then 
    Write(ErrCan,*) "Numerical problem in MonopoleSFF"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   MonopoleSFF =1/eps
  Else If (Abs(m1).lt.eps) Then 
    Write(ErrCan,*) "Numerical problem in MonopoleSFF"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
   MonopoleSFF =1/eps
  Else 
   MonopoleSFF = ((m1 - m3)*(16*m1**2 - 29*m1*m3 + 7*m3**2) + 6*m1**2*(3*m3*Log(m1/m3) + 2*m1*Log(m3/m1)))/(36.*(m1 - m3)**4)
  End if

Else If (SmallDifference(m1,m3)) Then
    Write(ErrCan,*) "Numerical problem in MonopoleSFF"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
    MonopoleSFF = 1/eps

Else If (SmallDifference(m2,m3)) Then
    Write(ErrCan,*) "Numerical problem in MonopoleSFF"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
    MonopoleSFF = 1/eps

Else 
 If (Abs(m1).lt.eps) Then
   MonopoleSFF = (-m2**2 + m3**2 + 2*m2*m3*Log(m2/m3))/(6.*(m2 - m3)**3)  
 Else If (Abs(m2).lt.eps) Then
   MonopoleSFF = (m1 - m3 + m1*Log(m3/m1))/(6.*(m1 - m3)**2)
 Else If (Abs(m3).lt.eps) Then
   MonopoleSFF = (m1 - m2 + m1*Log(m2/m1))/(6.*(m1 - m2)**2)
 Else
  MonopoleSFF = (-2*m2**4*(m1 - m3)**2 - 3*m1**3*Sqrt(m2**5*m3) &
     & + 3*m1**2*Sqrt(m2**7*m3) + 3*m1**3*Sqrt(m2*m3**5) - 3*m1**2*Sqrt(m2*m3**7) + &
     &    3*m1**2*Sqrt(m2**7*m3)*Log(m1/m2) - 3*m2*(m2*m3)**2.5*Log(m2) &
     & - 3*m3*(m2*m3)**2.5*Log(m2) + 2*m1**2*m3**3*(-m1 + m3 + m1*Log(m1/m3)) + & 
     &    2*m1*m2*m3**2*((m1 - m3)*(m1 + 2*m3) - 3*m1**2*Log(m1/m3)) &
     & + 6*m1**3*(m2*m3)**1.5*Log(m2/m3) + 6*m1*(m2*m3)**2.5*Log(m2/m3) +  &
     &    3*m1*m2**1.5*m3**2.5*(-3*m1 + 3*m3 + 3*m1*Log(m1) &
     & + 2*(-2*m1 + m3)*Log(m2) + (m1 - 2*m3)*Log(m3)) +  &
     &    3*m2**3.5*m3**1.5*(-3*m1 + 2*m3 + 2*m1*Log(m2) + (-2*m1 + m3)*Log(m3)) +  &
     &    3*m2**2.5*m3**1.5*(3*m1**2 - 2*m3**2 - m1**2*Log(m1**3*m2)  &
     & + (4*m1**2 + m3**2)*Log(m3)) + 3*m1**2*Sqrt(m2*m3**7)*Log(m3/m1) -  &
     &    2*m2**3*(m1**2*(-m1 + m3) + m1**3*Log(m1/m2) + 3*m1*m3**2*Log(m2/m3) &
     & + 2*m3**3*Log(m3/m2)) +  &
     &    2*m2**2*m3*(-m1**3 + m3**3 + 3*m1**3*Log(m1/m2) &
     & + 3*m1*m3*(2*m1*Log(m2/m3) + m3*Log(m3/m2))))      &
     & /(6.*(m1 - m2)**2*(m1 - m3)**2*(m2 - m3)**3)
  End if
End If

End Function MonopoleSFF


Real(dp) Function MonopoleVFF(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: eps=1E-10_dp

If ((SmallDifference(m2,m3)).and.(SmallDifference(m1,m2))) Then
 If (Abs(m1).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in MonopoleVFF"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
    MonopoleVFF = 1/eps
 Else 
  MonopoleVFF = -7/(12.*m1)
 End if

Else If (SmallDifference(m2,m3)) Then
 If (Abs(m1).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in MonopoleVFF"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
    MonopoleVFF = 1/eps
 Else If (Abs(m2).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in MonopoleVFF"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
    MonopoleVFF = 1/eps
 Else
  MonopoleVFF = (-((m1 - m3)*(2*m1**2 + 29*m1*m3 - 25*m3**2)) + &
     & 6*m1*(2*m1**2 - 9*m1*m3 + 6*m3**2)*Log(m3/m1))/(18.*(m1 - m3)**4)
  End if

Else If (SmallDifference(m1,m3)) Then
    Write(ErrCan,*) "Numerical problem in MonopoleVFF"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
    MonopoleVFF = 1/eps

Else If (SmallDifference(m1,m2)) Then
    Write(ErrCan,*) "Numerical problem in MonopoleVFF"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
    MonopoleVFF = 1/eps

Else 
 If (Abs(m1).lt. eps) Then
  MonopoleVFF =  (2*(-m2**2 + m3**2 + 3*Sqrt(m2**3*m3) - 3*Sqrt(m2*m3**3)) &
     &  - 3*Sqrt(m2*m3)*(m2 + m3)*Log(m2/m3) +  &
     &    (3*m2**2 - 10*m2*m3 + 3*m3**2)*Log(m3/m2))/(3.*(m2 - m3)**3)
 Else If (Abs(m2).lt. eps) Then
  MonopoleVFF = (-m1 + m3 + (-2*m1 + 3*m3)*Log(m1/m3))/(3.*(m1 - m3)**2)
 Else If (Abs(m3).lt. eps) Then
  MonopoleVFF = (-m1 + m2 + (-2*m1 + 3*m2)*Log(m1/m2))/(3.*(m1 - m2)**2)
 Else 
  MonopoleVFF = (-3*m1**3*Sqrt(m2**5*m3) + 3*m1**2*Sqrt(m2**7*m3) &
     & + 3*m1**3*Sqrt(m2*m3**5) - 3*m1**2*Sqrt(m2*m3**7) &
     & + 9*m1**2*(m2*m3)**1.5*(-m2 + m3)*Log(m1) + &
     &    m1**2*(3*m2**4 + 6*m1*m2**2*m3 + 3*Sqrt(m2**7*m3) &
     & - 2*m2**3*(m1 + 3*m3))*Log(m1/m2) + m1**2*m3**3 &
     &       *(m1 - m3 + (2*m1 - 3*m3)*Log(m1/m3)) +  &
     &    6*m1**3*(m2*m3)**1.5*Log(m2/m3) + 6*m1*(m2*m3)**2.5*Log(m2/m3) -  &
     &    3*m2**2.5*m3**1.5*(-3*m1**2 + 2*m3**2 + m1**2*Log(m2/m3**4) &
     & + m3**2*Log(m2/m3)) - 3*m2**3.5*m3**1.5  &
     & *(3*m1 - 2*m3 + (-2*m1 + m3)*Log(m2/m3)) +  &
     &    3*m1*m2**1.5*m3**2.5*(-3*m1 + 3*m3 + 2*(-2*m1 + m3)*Log(m2) &
     & + (m1 - 2*m3)*Log(m3)) + 3*m1**2*Sqrt(m2*m3**7)*Log(m3/m1) +  &
     &    m1*m2*m3**2*(-7*m1**2 + 8*m1*m3 - m3**2 - 6*m1**2*Log(m1/m3) &
     & - 6*m3*(m1 + m3)*Log(m3/m1)) +  &
     &    m2**4*(m1**2 + m1*m3 - 2*m3**2 + 6*m1*m3*Log(m2/m1) &
     &  + 3*m3**2*Log(m3/m2)) +   &
     &    m2**2*m3*(7*m1**3 - 9*m1*m3**2 + 2*m3**3 &
     & + 9*m1*m3*(m1*Log(m2/m3) + 2*m3*Log(m3/m1)) + 3*m3**3*Log(m3/m2)) -   &
     &    m2**3*(m1*(m1 - m3)*(m1 + 9*m3) + 2*m3**2*(9*m1*Log(m2/m1) &
     &  + 5*m3*Log(m3/m2))))/(3.*(m1 - m2)**2*(m1 - m3)**2*(m2 - m3)**3)  
 End If
End if

End Function MonopoleVFF

Real(dp) Function MonopoleFVV(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: eps=1E-10_dp

If ((SmallDifference(m2,m3)).and.(SmallDifference(m1,m2))) Then
 If (Abs(m1).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in MonopoleFVV"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
  MonopoleFVV = 1/eps
 Else
  MonopoleFVV = 1/(3.*m1)
 End if

Else If (SmallDifference(m2,m3)) Then
 If (Abs(m1).lt.eps) Then
  MonopoleFVV = 5/(9.*m2)
 Else If (Abs(m2).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in MonopoleFVV"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
  MonopoleFVV = 1/eps
 Else
  MonopoleFVV = (-((m1 - m2)*(5*m1**2 - 22*m1*m2 + 5*m2**2)) + 6*m1**2*(m1 - 3*m2)*Log(m1/m2))/(9.*(m1 - m2)**4)
 End if
!! Assuming always identical vector bosons
Else 
 
 If (Abs(m1).lt.eps) Then
  MonopoleFVV = 5/(9.*m2)
 Else If (Abs(m2).lt.eps) Then
    Write(ErrCan,*) "Numerical problem in MonopoleFVV"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
  MonopoleFVV = 1/eps
 Else
  MonopoleFVV = (-((m1 - m2)*(5*m1**2 - 22*m1*m2 + 5*m2**2)) &
           & + 6*m1**2*(m1 - 3*m2)*Log(m1/m2))/(9.*(m1 - m2)**4)
 End if

End If

End Function MonopoleFVV

Real(dp) Function MonopoleFVS(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
Real(dp) :: eps=1E-10_dp

If ((SmallDifference(m2,m3)).and.(SmallDifference(m1,m2))) Then
 If (Abs(m1).lt.eps) Then 
    Write(ErrCan,*) "Numerical problem in MonopoleFVS"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
  MonopoleFVS = 1/eps
 Else
  MonopoleFVS = 1/(24.*m2**1.5)
 End if

Else If (SmallDifference(m2,m3)) Then
  If (Abs(m1).lt.eps) Then 
   MonopoleFVS = 0
 Else if (Abs(m2).lt.eps) Then 
    Write(ErrCan,*) "Numerical problem in MonopoleFVS"
    Write(ErrCan,*) "Involved masses: ",m1,m2,m3
    MonopoleFVS = 1/eps
 Else
   MonopoleFVS =(Sqrt(m1)*(2*m1**3 + 3*m1**2*m3 - 6*m1*m3**2 + m3**3 &
            & + 6*m1**2*m3*Log(m3/m1)))/(12.*(m1 - m3)**4*m3)
 End if

Else If (SmallDifference(m1,m3)) Then
 If (Abs(m1).lt.eps) Then 
    MonopoleFVS = 0
 Else if (Abs(m2).lt.eps) Then 
    MonopoleFVS = 1/(4.*m3**1.5)
 Else
    MonopoleFVS = (Sqrt(m3)*(-5*m2**2 + 4*m2*m3 + m3**2 &
       & + 2*m2*(m2 + 2*m3)*Log(m2/m3)))/(4.*(m2 - m3)**4)
 End if

Else If (SmallDifference(m1,m2)) Then
 If (Abs(m1).lt.eps) Then 
    MonopoleFVS = 0
  Else if (Abs(m3).lt.eps) Then 
    MonopoleFVS = 1/(4.*m2**1.5)
  Else
    MonopoleFVS = (Sqrt(m2)*((m2 - m3)*(m2 + 5*m3) &
             & + 2*m3*(2*m2 + m3)*Log(m3/m2)))/(4.*(m2 - m3)**4)
  End if

Else 
  If (Abs(m1).lt.eps) Then 
    MonopoleFVS = 0
  Else if (Abs(m2).lt.eps) Then 
    MonopoleFVS = (Sqrt(m1)*(m1 - m3 - m3*Log(m1) + m3*Log(m3)))/(2.*(m1 - m3)**2*m3)
  Else if (Abs(m3).lt.eps) Then 
    MonopoleFVS = (Sqrt(m1)*(m1 - m2 - m2*Log(m1) + m2*Log(m2)))/(2.*(m1 - m2)**2*m2)
  Else
  MonopoleFVS =  (Sqrt(m1)*((m1 - m2)*(m1 - m3)*(m2 - m3)*(-2*m2*m3 + m1*(m2 + m3))&
     &  - m1**2*(m2 - m3)**3*Log(m1) +  &
     &      m2*(m1 - m3)**2*(-2*m1*m3 + m2*(m2 + m3))*Log(m2) &
     & + (m1 - m2)**2*m3*(2*m1*m2 - m3*(m2 + m3))*Log(m3)))/ &
     &  (2.*(m1 - m2)**2*(m1 - m3)**2*(m2 - m3)**3)
  End if
End If

End Function MonopoleFVS

Real(dp) Function MonopoleFSV(m1, m2, m3)
Implicit None
Real(dp), Intent(in) :: m1, m2, m3
 MonopoleFSV = MonopoleFVS(m1,m3,m2)
End Function MonopoleFSV


End Module LoopFunctions
