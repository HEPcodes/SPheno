Module RGEs

! load modules
 Use Control
 Use Mathematics
! load modules

! global variables
  !-------------------------------------------------------------------
  ! coefficients for gauge contribution in RGEs for gauge couplings
  !-------------------------------------------------------------------
  Real(dp), Parameter :: b_1(3) = (/ 6.6_dp, 1._dp, -3._dp /)  &
    & , b_2(3,3) = Reshape( Source = (/ 199._dp / 25._dp, 1.8_dp, 2.2_dp    &
    &                                , 5.4_dp,           25._dp,  9._dp     &
    &                                , 17.6_dp,          24._dp, 14._dp /), &
    &                       Shape = (/3, 3/) )
  Real(dp) :: Delta_b_1(3) = 0._dp, Delta_b_2(3,3) = 0._dp
  !-------------------------------------------------------------------
  ! coefficients for Yukawa contribution in RGEs for gauge couplings
  ! the matrix is revised compared to usual notation: tau,bottom,top
  !-------------------------------------------------------------------
  Real(dp), Parameter :: a_2(3,3) =     &
    &       Reshape( Source = (/  3.6_dp,  2._dp,  0._dp       &
    &                          ,  2.8_dp,  6._dp,  4._dp       &
    &                          ,  5.2_dp,  6._dp,  4._dp  /),  &
    &                Shape = (/3, 3/) )
  ! including right handed neutrinos: tau,neutrino,bottom,top
  Real(dp), Parameter :: a_2a(3,4) =     &
    &       Reshape( Source = (/  3.6_dp,  2._dp,  0._dp       &
    &                          ,  1.2_dp,  2._dp,  0._dp       &
    &                          ,  2.8_dp,  6._dp,  4._dp       &
    &                          ,  5.2_dp,  6._dp,  4._dp  /),  &
    &                Shape = (/3, 4/) )
  ! including complete 15-plet: tau, T, bottom,top, S, Z, lam_1, lam_2
  Real(dp), Parameter :: a_2b(3,8) =     &
    &       Reshape( Source = (/  3.6_dp,  2._dp,  0._dp       &
    &                          ,  5.4_dp,  7._dp,  0._dp       &
    &                          ,  2.8_dp,  6._dp,  4._dp       &
    &                          ,  5.2_dp,  6._dp,  4._dp        &
    &                          ,  4.8_dp,  0._dp,  9._dp        &
    &                          ,  2.8_dp,  6._dp,  4._dp        &
    &                          ,  5.4_dp,  7._dp,  0._dp        &
    &                          ,  5.4_dp,  7._dp,  0._dp   /),  &
    &                Shape = (/3, 8/) )
  !-------------------------------------------------------------------
  ! coefficients for contribution in RGEs for Yukawa couplings
  ! the matrix is revised compared to usual notation: tau,bottom,top
  !-------------------------------------------------------------------
  ! gauge contributions
  !-------------------------
  Real(dp), Parameter :: c1_1(3,3) =     &
    &       Reshape( Source = (/-1.8_dp,  -7._dp /15._dp ,-13._dp /15._dp  &
    &                        , -3._dp  ,  -3._dp,           -3._dp         &
    &                        ,  0._dp  , -16._dp/3._dp,   -16._dp/3._dp /) &
    &              , Shape = (/3, 3/) )
  !-------------------------
  ! Yukawa contributions
  !-------------------------
  Real(dp), Parameter :: c2_1(3,3) =     &
    &       Reshape( Source = (/  4._dp,  1._dp,  0._dp    &
    &                          ,  3._dp,  6._dp,  1._dp    &
    &                          ,  0._dp,  1._dp,  6._dp /) &
    &              , Shape = (/3, 3/) )

 Logical, Save :: TwoLoopRGE=.True.
 Real(dp) :: M2S_GUT
 Complex(dp) :: Alam_GUT, Alamp_GUT
! global variables

! private variables
! simplifies matrix multiplication in case of diagonal entries only
Logical, Private, Save :: OnlyDiagonal
! fix to completely decouple heavy states
Logical, Private, Save :: decoupling_heavy_states = .False.
! private variables

!------------------------------------------------
! needed for programs by Florian Staub
!------------------------------------------------
# ifdef SEESAWIII
! Real(dp), parameter :: sqrt0d3 = Sqrt(0.3_dp), sqrt1d2 = Sqrt(1.2_dp) &
!        &  , sqrt7d5 = Sqrt(7.5_dp), sqrt30 = Sqrt(30._dp)             &
!        &  , sqrt4o3 = Sqrt(4._dp/3._dp)
 Real(dp), parameter :: & 
         &   sqrt0d3 = 0.54772255750516611345696978280080213395274469499798_dp &
         & , sqrt1d2 = 1.0954451150103322269139395656016042679054893899960_dp  &
         & , sqrt7d5 = 2.7386127875258305672848489140040106697637234749899_dp  &
         & , sqrt30  = 5.4772255750516611345696978280080213395274469499798_dp  &
         & , sqrt4o3 = 1.1547005383792515290182975610039149112952035025403_dp

 Real(dp), Save :: NGHb3,NGHg3, NGHw3, NGHx3,NGHxb3
 Integer, save :: ThresholdCrossed = 1
 Real(dp),Parameter::id3R(3,3)= Reshape(Source=(/ 1,0,0,&
                             &                    0,1,0,&
                             &                    0,0,1 &
                             & /),shape=(/3,3/))

# endif SEESAWIII
Contains


 Subroutine CouplingsToG(gauge, yuk_l, yuk_d, yuk_u, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings to a vector, splitting complex parameters into
 ! real and imaginary part
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)   ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_d(3,3)  & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)    ! u-quark Yukawa couplings
  Real(dp), Intent(out) :: g1(57)          ! vector containing these couplings

  Integer :: i1, i2, sumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CouplingsToG'

  g1(1:3) = gauge
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_d(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_d(i1,i2))
    g1(SumI+32) = Real(yuk_u(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_u(i1,i2))
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine CouplingsToG

 Subroutine CouplingsToG2(gauge, yuk_l, yuk_nu, yuk_d, yuk_u, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings to a vector, splitting complex parameters into
 ! real and imaginary part, including Y_nu for MSSM
 ! written by Werner Porod, 04.03.2001
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)   ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3) & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)  & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)    ! u-quark Yukawa couplings
  Real(dp), Intent(out) :: g1(75)          ! vector containing these couplings

  Integer :: i1, i2, sumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CouplingsToG2'

  g1(1:3) = gauge
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_nu(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_nu(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine CouplingsToG2

 Subroutine CouplingsToG3(gauge, yuk_l, yuk_nu, yuk_d, yuk_u, Mnu, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings to a vector, splitting complex parameters into
 ! real and imaginary part, including Y_nu for MSSM
 ! written by Werner Porod
 ! 30.08.02: taking CouplingsToG2 as basis, adding left neutrino mass
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3) & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)  & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)  & ! u-quark Yukawa couplings
                         & ,  Mnu(3,3)     ! dim-5 neutrino mass operator
  Real(dp), Intent(out) :: g1(93)          ! vector containing these couplings

  Integer :: i1, i2, sumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CouplingsToG3'

  g1(1:3) = gauge
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_nu(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_nu(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+68) = Real(Mnu(i1,i2),dp)
    g1(SumI+69) = Aimag(Mnu(i1,i2))
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine CouplingsToG3

 Subroutine CouplingsToG4(gauge, yuk_l, yuk_T, yuk_d, yuk_u, lam1, lam2, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings to a vector, splitting complex parameters into
 ! real and imaginary part, including Y_T for MSSM
 ! written by Werner Porod, 04.03.2001
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)   ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)  & ! triplet Yukawa couplings
                         & , yuk_d(3,3)  & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)  & ! u-quark Yukawa couplings
                         & , lam1        & ! coupling triplet H_d
                         & , lam2          ! coupling triplet H_u
  Real(dp), Intent(out) :: g1(79)          ! vector containing these couplings

  Integer :: i1, i2, sumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CouplingsToG4'

  g1(1:3) = gauge
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_T(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_T(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
   End Do
  End Do
  g1(76) = Real(lam1,dp)
  g1(77) = Aimag(lam1)
  g1(78) = Real(lam2,dp)
  g1(79) = Aimag(lam2)

  Iname = Iname - 1

 End Subroutine CouplingsToG4

 Subroutine CouplingsToG5(gauge, yuk_l, yuk_T, yuk_d, yuk_u, yuk_Z, yuk_S &
             & , lam1, lam2, M15, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings to a vector, splitting complex parameters into
 ! real and imaginary part, including Y_T for MSSM
 ! written by Werner Porod, 04.03.2001
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)   ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)  & ! triplet Yukawa couplings
                         & , yuk_d(3,3)  & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)  & ! u-quark Yukawa couplings
                         & , yuk_Z(3,3)   & ! triplet Yukawa couplings 
                         & , yuk_S(3,3)   & ! triplet Yukawa couplings
                         & , lam1        & ! coupling triplet H_d
                         & , lam2          ! coupling triplet H_u
  Real(dp), Intent(in) :: M15(3)            ! M_T, M_Z, M_S 15-plet masses
  Real(dp), Intent(out) :: g1(118)          ! vector containing these couplings

  Integer :: i1, i2, sumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CouplingsToG5'

  g1(1:3) = gauge
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_T(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_T(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+68) = Real(yuk_Z(i1,i2),dp)
    g1(SumI+69) = Aimag(yuk_Z(i1,i2))
    g1(SumI+86) = Real(yuk_S(i1,i2),dp)
    g1(SumI+87) = Aimag(yuk_S(i1,i2))
   End Do
  End Do
  g1(112) = Real(lam1,dp)
  g1(113) = Aimag(lam1)
  g1(114) = Real(lam2,dp)
  g1(115) = Aimag(lam2)
  g1(116:118) = M15
  Iname = Iname - 1

 End Subroutine CouplingsToG5

 Subroutine GToCouplings(g1, gauge, yuk_l, yuk_d, yuk_u)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings, putting real and imaginary part
 ! together
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(57)            ! vector containing these couplings
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
 
  Integer :: i1, i2, SumI

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToCouplings'

  gauge = g1(1:3)
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine GToCouplings

 Subroutine GToCouplings2(g1, gauge, yuk_l, yuk_nu, yuk_d, yuk_u)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings, putting real and imaginary part
 ! together
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(75)            ! vector containing these couplings
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3) &  ! neutrino Yukawa couplings
                         & , yuk_d(3,3)  &  ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
 
  Integer :: i1, i2, SumI

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToCouplings'

  gauge = g1(1:3)
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_nu(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2)  = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2)  = Cmplx( g1(SumI+50), g1(SumI+51),dp )
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine GToCouplings2

 Subroutine GToCouplings3(g1, gauge, yuk_l, yuk_nu, yuk_d, yuk_u, Mnu)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings, putting real and imaginary part
 ! together
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 ! 30.08.02: adding left neutrino mass matrix, taking GToCouplings2 as basis
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(93)            ! vector containing these couplings
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3)  & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & ,  Mnu(3,3)      ! dim-5 neutrino mass operator
 
  Integer :: i1, i2, SumI

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToCouplings3'

  gauge = g1(1:3)
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_nu(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2)  = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2)  = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    Mnu(i1,i2)  = Cmplx( g1(SumI+68), g1(SumI+69),dp )
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine GToCouplings3

 Subroutine GToCouplings4(g1, gauge, yuk_l, yuk_T, yuk_d, yuk_u, lam1, lam2)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings, putting real and imaginary part
 ! together
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(79)            ! vector containing these couplings
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
 
  Integer :: i1, i2, SumI

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToCouplings4'

  gauge = g1(1:3)
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_T(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2)  = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2)  = Cmplx( g1(SumI+50), g1(SumI+51),dp )
   End Do
  End Do
  lam1  = Cmplx( g1(76), g1(77),dp )
  lam2  = Cmplx( g1(78), g1(79),dp )

  Iname = Iname - 1

 End Subroutine GToCouplings4

 Subroutine GToCouplings5(g1, gauge, yuk_l, yuk_T, yuk_d, yuk_u, yuk_Z, yuk_S &
             & , lam1, lam2, M15)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings, putting real and imaginary part
 ! together
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(118)           ! vector containing these couplings
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , yuk_Z(3,3)   & ! triplet Yukawa couplings 
                         & , yuk_S(3,3)   & ! triplet Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
  Real(dp), Intent(out) :: M15(3)           ! M_T, M_Z, M_S 15-plet masses
 
  Integer :: i1, i2, SumI

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToCouplings5'

  gauge = g1(1:3)
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_T(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2)  = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2)  = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    yuk_Z(i1,i2)  = Cmplx( g1(SumI+68), g1(SumI+69),dp )
    yuk_S(i1,i2)  = Cmplx( g1(SumI+86), g1(SumI+87),dp )
   End Do
  End Do
  lam1  = Cmplx( g1(112), g1(113),dp )
  lam2  = Cmplx( g1(114), g1(115),dp )
  M15 = g1(116:118)

  Iname = Iname - 1

 End Subroutine GToCouplings5


 Subroutine GToParameters(g1, gauge, yuk_l, yuk_d, yuk_u  &
                         &, Mhlf, Ae, Ad, Au, Me, Ml, Md, Mq, Mu, Mh, mue, B)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings ans susy mass parameters, 
 ! putting real and imaginary part together
 ! written by Werner Porod, 20.8.99 
 ! 24.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(213)            ! vector containing the parameters
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(out) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(out) :: Ae(3,3)  & ! lepton A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(out) :: Me(3,3) & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3)   ! R u-squark mass parameters squared
  Real(dp), Intent(out) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(out) :: mue, B    ! mu, mu*B

  Integer :: i1, i2, SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToParameters'

  gauge = g1(1:3)
  Do i1=1,3
   Mhlf(i1) = Cmplx( g1(56 + 2*i1), g1(57 + 2*i1),dp )
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    Ae(i1,i2) = Cmplx( g1(SumI+56), g1(SumI+57),dp )
    Ad(i1,i2) = Cmplx( g1(SumI+74), g1(SumI+75),dp )
    Au(i1,i2) = Cmplx( g1(SumI+92), g1(SumI+93),dp )
    Me(i1,i2) = Cmplx( g1(SumI+110), g1(SumI+111),dp )
    Ml(i1,i2) = Cmplx( g1(SumI+128), g1(SumI+129),dp )
    Md(i1,i2) = Cmplx( g1(SumI+146), g1(SumI+147),dp )
    Mq(i1,i2) = Cmplx( g1(SumI+164), g1(SumI+165),dp )
    Mu(i1,i2) = Cmplx( g1(SumI+182), g1(SumI+183),dp )
   End Do
  End Do
  mH = g1(208:209)
  mue = Cmplx( g1(210), g1(211),dp )
  B = Cmplx( g1(212), g1(213),dp )

  Iname = Iname - 1

 End Subroutine GToParameters

 Subroutine GToParameters2(g1, gauge, yuk_l, yuk_nu, yuk_d, yuk_u  &
        &, Mhlf, Ae, Anu, Ad, Au, Me, Ml, Mr, Md, Mq, Mu, Mh, mue, B)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings ans susy mass parameters, 
 ! putting real and imaginary part together
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 24.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(267)            ! vector containing the parameters
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3)  & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(out) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(out) :: Ae(3,3)  & ! lepton A parameters
                         & , Anu(3,3)  & ! neutrino A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(out) :: Me(3,3) & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Mr(3,3) & ! R sneutrino mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3)   ! R u-squark mass parameters squared
  Real(dp), Intent(out) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(out) :: mue, B    ! mu, mu*B

  Integer :: i1, i2, SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToParameters2'

  gauge = g1(1:3)
  Do i1=1,3
   Mhlf(i1) = Cmplx( g1(74 + 2*i1), g1(75 + 2*i1),dp )
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_nu(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    Ae(i1,i2) = Cmplx( g1(SumI+74), g1(SumI+75),dp )
    Anu(i1,i2) = Cmplx( g1(SumI+92), g1(SumI+93),dp )
    Ad(i1,i2) = Cmplx( g1(SumI+110), g1(SumI+111),dp )
    Au(i1,i2) = Cmplx( g1(SumI+128), g1(SumI+129),dp )
    Me(i1,i2) = Cmplx( g1(SumI+146), g1(SumI+147),dp )
    Ml(i1,i2) = Cmplx( g1(SumI+164), g1(SumI+165),dp )
    Mr(i1,i2) = Cmplx( g1(SumI+182), g1(SumI+183),dp )
    Md(i1,i2) = Cmplx( g1(SumI+200), g1(SumI+201),dp )
    Mq(i1,i2) = Cmplx( g1(SumI+218), g1(SumI+219),dp )
    Mu(i1,i2) = Cmplx( g1(SumI+236), g1(SumI+237),dp )
   End Do
  End Do
  mH = g1(262:263)
  mue = Cmplx( g1(264), g1(265),dp )
  B = Cmplx( g1(266), g1(267),dp )

  Iname = Iname - 1

 End Subroutine GToParameters2

 Subroutine GToParameters3(g1, gauge, yuk_l, yuk_nu, yuk_d, yuk_u  &
        & , Mhlf, Ae, Anu, Ad, Au, Me, Ml, Mr, Md, Mq, Mu, Mh, mue &
        & , B, MnuL)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings ans susy mass parameters, 
 ! putting real and imaginary part together
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 24.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(285)            ! vector containing the parameters
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3)  & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(out) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(out) :: Ae(3,3)  & ! lepton A parameters
                         & , Anu(3,3)  & ! neutrino A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(out) :: Me(3,3) & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Mr(3,3) & ! R sneutrino mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(out) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(out) :: mue, B    ! mu, mu*B

  Integer :: i1, i2, SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToParameters3'

  gauge = g1(1:3)
  Do i1=1,3
   Mhlf(i1) = Cmplx( g1(74 + 2*i1), g1(75 + 2*i1),dp )
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_nu(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    Ae(i1,i2) = Cmplx( g1(SumI+74), g1(SumI+75),dp )
    Anu(i1,i2) = Cmplx( g1(SumI+92), g1(SumI+93),dp )
    Ad(i1,i2) = Cmplx( g1(SumI+110), g1(SumI+111),dp )
    Au(i1,i2) = Cmplx( g1(SumI+128), g1(SumI+129),dp )
    Me(i1,i2) = Cmplx( g1(SumI+146), g1(SumI+147),dp )
    Ml(i1,i2) = Cmplx( g1(SumI+164), g1(SumI+165),dp )
    Mr(i1,i2) = Cmplx( g1(SumI+182), g1(SumI+183),dp )
    Md(i1,i2) = Cmplx( g1(SumI+200), g1(SumI+201),dp )
    Mq(i1,i2) = Cmplx( g1(SumI+218), g1(SumI+219),dp )
    Mu(i1,i2) = Cmplx( g1(SumI+236), g1(SumI+237),dp )
    MnuL(i1,i2) = Cmplx( g1(SumI+260), g1(SumI+261),dp )
   End Do
  End Do
  mH = g1(262:263)
  mue = Cmplx( g1(264), g1(265),dp )
  B = Cmplx( g1(266), g1(267),dp )

  Iname = Iname - 1

 End Subroutine GToParameters3

 Subroutine GToParameters4(g1, gauge, yuk_l, yuk_T, yuk_d, yuk_u, lam1, lam2  &
         & ,  Mhlf, Ae, AT, Ad, Au, Alam1, Alam2, Me, Ml, Md, Mq, Mu, Mh  &
         & , MT, mue, B, MnuL)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings ans susy mass parameters, 
 ! putting real and imaginary part together
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 24.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(277)           ! vector containing the parameters
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
  Complex(dp), Intent(out) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(out) :: Ae(3,3)  & ! lepton A parameters
                         & , AT(3,3)   & ! triplet A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)   & ! u-quark A parameters
                         & , Alam1     & ! A-parameter triplet H_d
                         & , Alam2       ! A-parameter triplet H_u
  Complex(dp), Intent(out) :: Me(3,3) & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(out) :: Mh(2)        ! soft Higgs masses squared
  Real(dp), Intent(out) :: MT(2)        ! soft Higgs triplet masses squared
  Complex(dp), Intent(out) :: mue, B    ! mu, mu*B

  Integer :: i1, i2, SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToParameters4'

  gauge = g1(1:3)
  lam1 = Cmplx(g1(76),g1(77),dp)
  lam2 = Cmplx(g1(78),g1(79),dp)
  Alam1 = Cmplx(g1(158),g1(159),dp)
  Alam2 = Cmplx(g1(160),g1(161),dp)
  Do i1=1,3
   Mhlf(i1) = Cmplx( g1(78 + 2*i1), g1(79 + 2*i1),dp )
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_T(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    Ae(i1,i2) = Cmplx( g1(SumI+78), g1(SumI+79),dp )
    AT(i1,i2) = Cmplx( g1(SumI+96), g1(SumI+97),dp )
    Ad(i1,i2) = Cmplx( g1(SumI+114), g1(SumI+115),dp )
    Au(i1,i2) = Cmplx( g1(SumI+132), g1(SumI+133),dp )
    Me(i1,i2) = Cmplx( g1(SumI+154), g1(SumI+155),dp )
    Ml(i1,i2) = Cmplx( g1(SumI+172), g1(SumI+173),dp )
    Md(i1,i2) = Cmplx( g1(SumI+190), g1(SumI+191),dp )
    Mq(i1,i2) = Cmplx( g1(SumI+208), g1(SumI+209),dp )
    Mu(i1,i2) = Cmplx( g1(SumI+226), g1(SumI+227),dp )
    MnuL(i1,i2) = Cmplx( g1(SumI+252), g1(SumI+253),dp )
   End Do
  End Do
  mH = g1(252:253)
  mT = g1(254:255)
  mue = Cmplx( g1(256), g1(257),dp )
  B = Cmplx( g1(258), g1(259),dp )

  Iname = Iname - 1

 End Subroutine GToParameters4

 Subroutine GToParameters5(g1, gauge, yuk_l, yuk_T, yuk_d, yuk_u, yuk_Z, yuk_S &
         & , lam1, lam2, Mhlf, Ae, AT, Ad, Au, AZ, AS, Alam1, Alam2, Me, Ml    &
         & , Md, Mq, Mu, Mh, MT, mZ, mS, MT15, MZ15, MS15, mue, B, MnuL)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings ans susy mass parameters, 
 ! putting real and imaginary part together
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 24.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(356)           ! vector containing the parameters
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , yuk_Z(3,3)   & ! triplet Z Yukawa couplings
                         & , yuk_S(3,3)   & ! triplet S Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
  Complex(dp), Intent(out) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(out) :: Ae(3,3)  & ! lepton A parameters
                         & , AT(3,3)   & ! triplet A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)   & ! u-quark A parameters
                         & , AZ(3,3)   & ! triplet Z A parameters
                         & , AS(3,3)   & ! triplet S A parameters
                         & , Alam1     & ! A-parameter triplet H_d
                         & , Alam2       ! A-parameter triplet H_u
  Complex(dp), Intent(out) :: Me(3,3) & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(out) :: Mh(2)        ! soft Higgs masses squared
  Real(dp), Intent(out) :: MT(2)        ! soft Higgs triplet masses squared
  Real(dp), Intent(out) :: MZ(2)        ! soft Higgs Z triplet masses squared
  Real(dp), Intent(out) :: MS(2)        ! soft Higgs S triplet masses squared
  Real(dp), Intent(out) :: MT15         ! Higgs triplet masses squared
  Real(dp), Intent(out) :: MZ15         ! Higgs Z triplet masses squared
  Real(dp), Intent(out) :: MS15         ! Higgs S triplet masses squared
  Complex(dp), Intent(out) :: mue, B    ! mu, mu*B

  Integer :: i1, i2, SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToParameters5'

  gauge = g1(1:3)
  lam1 = Cmplx(g1(112),g1(113),dp)
  lam2 = Cmplx(g1(114),g1(115),dp)
  Alam1 = Cmplx(g1(230),g1(231),dp)
  Alam2 = Cmplx(g1(232),g1(233),dp)
  Do i1=1,3
   Mhlf(i1) = Cmplx( g1(114 + 2*i1), g1(115 + 2*i1),dp )
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_T(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    yuk_Z(i1,i2)  = Cmplx( g1(SumI+68), g1(SumI+69),dp )
    yuk_S(i1,i2)  = Cmplx( g1(SumI+86), g1(SumI+87),dp )

    Ae(i1,i2) = Cmplx( g1(SumI+114), g1(SumI+115),dp )
    AT(i1,i2) = Cmplx( g1(SumI+132), g1(SumI+133),dp )
    Ad(i1,i2) = Cmplx( g1(SumI+150), g1(SumI+151),dp )
    Au(i1,i2) = Cmplx( g1(SumI+168), g1(SumI+169),dp )
    AZ(i1,i2) = Cmplx( g1(SumI+186), g1(SumI+187),dp )
    AS(i1,i2) = Cmplx( g1(SumI+204), g1(SumI+205),dp )

    Me(i1,i2) = Cmplx( g1(SumI+226), g1(SumI+227),dp )
    Ml(i1,i2) = Cmplx( g1(SumI+244), g1(SumI+245),dp )
    Md(i1,i2) = Cmplx( g1(SumI+262), g1(SumI+263),dp )
    Mq(i1,i2) = Cmplx( g1(SumI+280), g1(SumI+281),dp )
    Mu(i1,i2) = Cmplx( g1(SumI+298), g1(SumI+299),dp )

    MnuL(i1,i2) = Cmplx( g1(SumI+331), g1(SumI+332),dp )
   End Do
  End Do

  mH = g1(324:325)
  mT = g1(326:327)
  mZ = g1(328:329)
  mS = g1(330:331)

  MT15 = g1(332)
  MZ15 = g1(333)
  MS15 = g1(334)
  
  mue = Cmplx( g1(335), g1(336),dp )
  B = Cmplx( g1(337), g1(338),dp )

  Iname = Iname - 1

 End Subroutine GToParameters5

 Subroutine ParametersToG(gauge, yuk_l, yuk_d, yuk_u                 &
               &,  Mhlf, Ae, Ad, Au, Me, Ml, Md, Mq, Mu, Mh, mue, B, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings and susy mass parameters to a vector, 
 ! splitting complex parameters into real and imaginary part
 ! written by Werner Porod, 20.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)          ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)   & ! lepton Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(in) :: Ae(3,3)   & ! lepton A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(in) :: Me(3,3)  & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3)   ! R u-squark mass parameters squared
  Real(dp), Intent(in) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(in) :: mue, B    ! mu, mu*B
  Real(dp), Intent(out) :: g1(213)            ! vector containing the parameters

  Integer i1,i2,SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ParametersToG'

  g1(1:3) = gauge
  Do i1=1,3
   g1(56 + 2*i1) = Real(Mhlf(i1),dp)
   g1(57 + 2*i1) = Aimag(Mhlf(i1))
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_d(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_d(i1,i2))
    g1(SumI+32) = Real(yuk_u(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_u(i1,i2))
    g1(SumI+56) = Real(Ae(i1,i2),dp)
    g1(SumI+57) = Aimag(Ae(i1,i2))
    g1(SumI+74) = Real(Ad(i1,i2),dp)
    g1(SumI+75) = Aimag(Ad(i1,i2))
    g1(SumI+92) = Real(Au(i1,i2),dp)
    g1(SumI+93) = Aimag(Au(i1,i2))
    g1(SumI+110) = Real(Me(i1,i2),dp)
    g1(SumI+111) = Aimag(Me(i1,i2))
    g1(SumI+128) = Real(Ml(i1,i2),dp)
    g1(SumI+129) = Aimag(Ml(i1,i2))
    g1(SumI+146) = Real(Md(i1,i2),dp)
    g1(SumI+147) = Aimag(Md(i1,i2))
    g1(SumI+164) = Real(Mq(i1,i2),dp)
    g1(SumI+165) = Aimag(Mq(i1,i2))
    g1(SumI+182) = Real(Mu(i1,i2),dp)
    g1(SumI+183) = Aimag(Mu(i1,i2))
   End Do
  End Do
  g1(208) = mH(1)
  g1(209) = mH(2)
  g1(210) = Real(mue,dp)
  g1(211) = Aimag(mue)
  g1(212) = Real(B,dp)
  g1(213) = Aimag(B)

  Iname = Iname - 1

 End Subroutine ParametersToG

 Subroutine ParametersToG2(gauge, yuk_l, yuk_nu, yuk_d, yuk_u             &
         &,  Mhlf, Ae, Anu, Ad, Au, Me, Ml, Mr, Md, Mq, Mu, Mh, mue, B, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings and susy mass parameters to a vector, 
 ! splitting complex parameters into real and imaginary part
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)          ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)   & ! lepton Yukawa couplings
                         & , yuk_nu(3,3)  & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(in) :: Ae(3,3)   & ! lepton A parameters
                         & , Anu(3,3)  & ! neutrino A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(in) :: Me(3,3)  & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Mr(3,3) & ! R sneutrino mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3)   ! R u-squark mass parameters squared
  Real(dp), Intent(in) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(in) :: mue, B    ! mu, mu*B
  Real(dp), Intent(out) :: g1(267)            ! vector containing the parameters

  Integer i1,i2,SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ParametersToG2'

  g1(1:3) = gauge
  Do i1=1,3
   g1(74 + 2*i1) = Real(Mhlf(i1),dp)
   g1(75 + 2*i1) = Aimag(Mhlf(i1))
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_nu(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_nu(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+74) = Real(Ae(i1,i2),dp)
    g1(SumI+75) = Aimag(Ae(i1,i2))
    g1(SumI+92) = Real(Anu(i1,i2),dp)
    g1(SumI+93) = Aimag(Anu(i1,i2))
    g1(SumI+110) = Real(Ad(i1,i2),dp)
    g1(SumI+111) = Aimag(Ad(i1,i2))
    g1(SumI+128) = Real(Au(i1,i2),dp)
    g1(SumI+129) = Aimag(Au(i1,i2))
    g1(SumI+146) = Real(Me(i1,i2),dp)
    g1(SumI+147) = Aimag(Me(i1,i2))
    g1(SumI+164) = Real(Ml(i1,i2),dp)
    g1(SumI+165) = Aimag(Ml(i1,i2))
    g1(SumI+182) = Real(Mr(i1,i2),dp)
    g1(SumI+183) = Aimag(Mr(i1,i2))
    g1(SumI+200) = Real(Md(i1,i2),dp)
    g1(SumI+201) = Aimag(Md(i1,i2))
    g1(SumI+218) = Real(Mq(i1,i2),dp)
    g1(SumI+219) = Aimag(Mq(i1,i2))
    g1(SumI+236) = Real(Mu(i1,i2),dp)
    g1(SumI+237) = Aimag(Mu(i1,i2))
   End Do
  End Do
  g1(262) = mH(1)
  g1(263) = mH(2)
  g1(264) = Real(mue,dp)
  g1(265) = Aimag(mue)
  g1(266) = Real(B,dp)
  g1(267) = Aimag(B)

  Iname = Iname - 1

 End Subroutine ParametersToG2

 Subroutine ParametersToG3(gauge, yuk_l, yuk_nu, yuk_d, yuk_u             &
         & ,  Mhlf, Ae, Anu, Ad, Au, Me, Ml, Mr, Md, Mq, Mu, Mh, mue, B   &
         & , MnuL5, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings and susy mass parameters to a vector, 
 ! splitting complex parameters into real and imaginary part
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 08.01.07: including dim5 operator for neutrinos
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)          ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)   & ! lepton Yukawa couplings
                         & , yuk_nu(3,3)  & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(in) :: Ae(3,3)   & ! lepton A parameters
                         & , Anu(3,3)  & ! neutrino A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(in) :: Me(3,3)  & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Mr(3,3) & ! R sneutrino mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL5(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(in) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(in) :: mue, B    ! mu, mu*B
  Real(dp), Intent(out) :: g1(285)            ! vector containing the parameters

  Integer i1,i2,SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ParametersToG3'

  g1(1:3) = gauge
  Do i1=1,3
   g1(74 + 2*i1) = Real(Mhlf(i1),dp)
   g1(75 + 2*i1) = Aimag(Mhlf(i1))
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_nu(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_nu(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+74) = Real(Ae(i1,i2),dp)
    g1(SumI+75) = Aimag(Ae(i1,i2))
    g1(SumI+92) = Real(Anu(i1,i2),dp)
    g1(SumI+93) = Aimag(Anu(i1,i2))
    g1(SumI+110) = Real(Ad(i1,i2),dp)
    g1(SumI+111) = Aimag(Ad(i1,i2))
    g1(SumI+128) = Real(Au(i1,i2),dp)
    g1(SumI+129) = Aimag(Au(i1,i2))
    g1(SumI+146) = Real(Me(i1,i2),dp)
    g1(SumI+147) = Aimag(Me(i1,i2))
    g1(SumI+164) = Real(Ml(i1,i2),dp)
    g1(SumI+165) = Aimag(Ml(i1,i2))
    g1(SumI+182) = Real(Mr(i1,i2),dp)
    g1(SumI+183) = Aimag(Mr(i1,i2))
    g1(SumI+200) = Real(Md(i1,i2),dp)
    g1(SumI+201) = Aimag(Md(i1,i2))
    g1(SumI+218) = Real(Mq(i1,i2),dp)
    g1(SumI+219) = Aimag(Mq(i1,i2))
    g1(SumI+236) = Real(Mu(i1,i2),dp)
    g1(SumI+237) = Aimag(Mu(i1,i2))
    g1(SumI+260) = Real(MnuL5(i1,i2),dp)
    g1(SumI+261) = Aimag(MnuL5(i1,i2))
   End Do
  End Do
  g1(262) = mH(1)
  g1(263) = mH(2)
  g1(264) = Real(mue,dp)
  g1(265) = Aimag(mue)
  g1(266) = Real(B,dp)
  g1(267) = Aimag(B)

  Iname = Iname - 1

 End Subroutine ParametersToG3

 Subroutine ParametersToG4(gauge, yuk_l, yuk_T, yuk_d, yuk_u, lam1, lam2  &
         & ,  Mhlf, Ae, AT, Ad, Au, Alam1, Alam2, Me, Ml, Md, Mq, Mu, Mh  &
         & , MT, mue, B, MnuL5, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings and susy mass parameters to a vector, 
 ! splitting complex parameters into real and imaginary part
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 08.01.07: including dim5 operator for neutrinos
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)          ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)   & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
  Complex(dp), Intent(in) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(in) :: Ae(3,3)  &  ! lepton A parameters
                         & , AT(3,3)   & ! triplet A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)   & ! u-quark A parameters
                         & , Alam1     & ! A-parameter triplet H_d
                         & , Alam2       ! A-parameter triplet H_u
  Complex(dp), Intent(in) :: Me(3,3) &  ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL5(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(in) :: Mh(2)         ! soft Higgs masses squared
  Real(dp), Intent(in) :: MT(2)         ! soft Higgs triplet masses squared
  Complex(dp), Intent(in) :: mue, B     ! mu, mu*B
  Real(dp), Intent(out) :: g1(277)           ! vector containing the parameters

  Integer i1,i2,SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ParametersToG4'

  g1(1:3) = gauge
  g1(76) = Real(lam1,dp)
  g1(77) = Aimag(lam1)
  g1(78) = Real(lam2,dp)
  g1(79) = Aimag(lam2)
  g1(158) = Real(Alam1,dp)
  g1(159) = Aimag(Alam1)
  g1(160) = Real(Alam2,dp)
  g1(161) = Aimag(Alam2)
  Do i1=1,3
   g1(78 + 2*i1) = Real(Mhlf(i1),dp)
   g1(79 + 2*i1) = Aimag(Mhlf(i1))
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_T(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_T(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+78) = Real(Ae(i1,i2),dp)
    g1(SumI+79) = Aimag(Ae(i1,i2))
    g1(SumI+96) = Real(AT(i1,i2),dp)
    g1(SumI+97) = Aimag(AT(i1,i2))
    g1(SumI+114) = Real(Ad(i1,i2),dp)
    g1(SumI+115) = Aimag(Ad(i1,i2))
    g1(SumI+132) = Real(Au(i1,i2),dp)
    g1(SumI+133) = Aimag(Au(i1,i2))
    g1(SumI+154) = Real(Me(i1,i2),dp)
    g1(SumI+155) = Aimag(Me(i1,i2))
    g1(SumI+172) = Real(Ml(i1,i2),dp)
    g1(SumI+173) = Aimag(Ml(i1,i2))
    g1(SumI+190) = Real(Md(i1,i2),dp)
    g1(SumI+191) = Aimag(Md(i1,i2))
    g1(SumI+208) = Real(Mq(i1,i2),dp)
    g1(SumI+209) = Aimag(Mq(i1,i2))
    g1(SumI+226) = Real(Mu(i1,i2),dp)
    g1(SumI+227) = Aimag(Mu(i1,i2))
    g1(SumI+252) = Real(MnuL5(i1,i2),dp)
    g1(SumI+253) = Aimag(MnuL5(i1,i2))
   End Do
  End Do
  g1(252) = mH(1)
  g1(253) = mH(2)
  g1(254) = mT(1)
  g1(255) = mT(2)
  g1(256) = Real(mue,dp)
  g1(257) = Aimag(mue)
  g1(258) = Real(B,dp)
  g1(259) = Aimag(B)

  Iname = Iname - 1

 End Subroutine ParametersToG4

 Subroutine ParametersToG5(gauge, yuk_l, yuk_T, yuk_d, yuk_u, yuk_Z, yuk_S  &
         & , lam1, lam2, Mhlf, Ae, AT, Ad, Au, AZ, AS, Alam1, Alam2, Me, Ml &
         & , Md, Mq, Mu, Mh, MT, MZ, MS, MT15, MZ15, MS15, mue, B, MnuL5, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings and susy mass parameters to a vector, 
 ! splitting complex parameters into real and imaginary part
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 08.01.07: including dim5 operator for neutrinos
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)          ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)   & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , yuk_Z(3,3)   & ! triplet Z Yukawa couplings
                         & , yuk_S(3,3)   & ! triplet S Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
  Complex(dp), Intent(in) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(in) :: Ae(3,3)  &  ! lepton A parameters
                         & , AT(3,3)   & ! triplet A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)   & ! u-quark A parameters
                         & , AZ(3,3)   & ! triplet Z A parameters
                         & , AS(3,3)   & ! triplet S A parameters
                         & , Alam1     & ! A-parameter triplet H_d
                         & , Alam2       ! A-parameter triplet H_u
  Complex(dp), Intent(in) :: Me(3,3) &  ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL5(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(in) :: Mh(2)         ! soft Higgs masses squared
  Real(dp), Intent(in) :: MT(2)         ! soft Higgs triplet masses squared
  Real(dp), Intent(in) :: MZ(2)         ! soft Higgs Z triplet masses squared
  Real(dp), Intent(in) :: MS(2)         ! soft Higgs S triplet masses squared
  Real(dp), Intent(in) :: MT15          ! Higgs triplet masses squared
  Real(dp), Intent(in) :: MZ15          ! Higgs Z triplet masses squared
  Real(dp), Intent(in) :: MS15          ! Higgs S triplet masses squared
  Complex(dp), Intent(in) :: mue, B     ! mu, mu*B
  Real(dp), Intent(out) :: g1(356)           ! vector containing the parameters

  Integer i1,i2,SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ParametersToG5'

  g1(1:3) = gauge
  g1(112) = Real(lam1,dp)
  g1(113) = Aimag(lam1)
  g1(114) = Real(lam2,dp)
  g1(115) = Aimag(lam2)

  g1(230) = Real(Alam1,dp)
  g1(231) = Aimag(Alam1)
  g1(232) = Real(Alam2,dp)
  g1(233) = Aimag(Alam2)
  Do i1=1,3
   g1(114 + 2*i1) = Real(Mhlf(i1),dp)
   g1(115 + 2*i1) = Aimag(Mhlf(i1))
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_T(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_T(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+68) = Real(yuk_Z(i1,i2),dp)
    g1(SumI+69) = Aimag(yuk_Z(i1,i2))
    g1(SumI+86) = Real(yuk_S(i1,i2),dp)
    g1(SumI+87) = Aimag(yuk_S(i1,i2))

    g1(SumI+114) = Real(Ae(i1,i2),dp)
    g1(SumI+115) = Aimag(Ae(i1,i2))
    g1(SumI+132) = Real(AT(i1,i2),dp)
    g1(SumI+133) = Aimag(AT(i1,i2))
    g1(SumI+150) = Real(Ad(i1,i2),dp)
    g1(SumI+151) = Aimag(Ad(i1,i2))
    g1(SumI+168) = Real(Au(i1,i2),dp)
    g1(SumI+169) = Aimag(Au(i1,i2))
    g1(SumI+186) = Real(AZ(i1,i2),dp)
    g1(SumI+187) = Aimag(AZ(i1,i2))
    g1(SumI+204) = Real(AS(i1,i2),dp)
    g1(SumI+205) = Aimag(AS(i1,i2))

    g1(SumI+226) = Real(Me(i1,i2),dp)
    g1(SumI+227) = Aimag(Me(i1,i2))
    g1(SumI+244) = Real(Ml(i1,i2),dp)
    g1(SumI+245) = Aimag(Ml(i1,i2))
    g1(SumI+262) = Real(Md(i1,i2),dp)
    g1(SumI+263) = Aimag(Md(i1,i2))
    g1(SumI+280) = Real(Mq(i1,i2),dp)
    g1(SumI+281) = Aimag(Mq(i1,i2))
    g1(SumI+298) = Real(Mu(i1,i2),dp)
    g1(SumI+299) = Aimag(Mu(i1,i2))

    g1(SumI+331) = Real(MnuL5(i1,i2),dp)
    g1(SumI+332) = Aimag(MnuL5(i1,i2))
   End Do
  End Do
  g1(324:325) = mH
  g1(326:327) = mT
  g1(328:329) = mZ
  g1(330:331) = mS
  g1(332) = MT15
  g1(333) = MZ15
  g1(334) = MS15
  g1(335) = Real(mue,dp)
  g1(336) = Aimag(mue)
  g1(337) = Real(B,dp)
  g1(338) = Aimag(B)

  Iname = Iname - 1

 End Subroutine ParametersToG5

 Subroutine rge6(len, T,GY,F)
 !-----------------------------------------------------------------
 !     Right hand side of truncated renormalization group equations
 !          dGY_i/dT = F_i(G)
 !     for the determination of M_GUT and the value of alpha_GUT
 !     and values of the Yukawas
 !  written by Werner Porod, 28.12.99
 !  6.10.2000: changing to f90
 !-----------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: len
  Real(Dp), Intent(in) :: T,GY(len)
  Real(Dp), Intent(out) :: F(len)

  Integer :: j,i
  Real(Dp) :: gy2(len),sumI,beta2(3), q

  q = t

  gy2 = gy**2

  Do i=1,3       ! gauge couplings two loop
   sumI = 0._dp
   Do j=1,3
    sumI = sumI + b_2(i,j) * gy2(j) - a_2(i,j) * gy2(3+j)
   Enddo
   f(i) = oo16pi2 * gy(i) * gy2(i) * ( b_1(i) + oo16pi2 * sumI )
  Enddo

  beta2(1) = ( 13.5_dp* gy2(1) + 1.8_dp * gy2(2) + 1.2_dp * gy2(4)   &
   &         - 0.4_dp * gy2(5) ) * gy2(1)                            &
   &       + (7.5_dp * gy2(2) + 6._dp * gy2(4) ) * gy2(2)            &
   &       + 16._dp * gy2(5) * gy2(3)                                &
   &       - (10._dp * gy2(4) + 9._dp * gy2(5) ) * gy2(4)            &
   &       - (9._dp * gy2(5) +  3._dp * gy2(6) ) * gy2(5)

  beta2(2) = ( 287._dp * gy2(1) / 90._dp + gy2(2) + 8._dp * gy2(3) / 9._dp    &
   &         + 1.2_dp * gy2(4) + 0.4_dp * gy2(5) + 0.8_dp * gy2(6) ) * gy2(1) &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(5) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(5) ) * gy2(3)                    &
   &       - 3._dp * (gy2(4) + gy2(5) ) * gy2(4)                              &
   &       - (22._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)                     &
   &       - 5._dp * gy2(6)**2

  beta2(3) = ( 2743._dp * gy2(1) / 450._dp + gy2(2) + 136._dp * gy2(3)/45._dp &
   &         + 0.4_dp * gy2(5) + 1.2_dp * gy2(6) ) * gy2(1)                   &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(6) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(6) ) * gy2(3)                    &
   &       - (gy2(4) + 5._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)             &
   &       - 22._dp * gy2(6)**2

   Do i=1,3       ! yukawa couplings two loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    f(i+3) = oo16pi2 * gy(i+3) * (sumI + oo16pi2 * beta2(i) )
   End Do

  End Subroutine rge6


 Subroutine rge7(len, T,GY,F)
 !-----------------------------------------------------------------
 !     Right hand side of truncated renormalization group equations
 !          dGY_i/dT = F_i(G)
 !     for the determination of M_GUT and the value of alpha_GUT
 !     and values of the Yukawas
 !  written by Werner Porod, 28.12.99
 !  10.01.00: including tan(beta)
 !  06.10.00: changing to f90
 !  27.07.13: adding gauge dependence as discussed in 1305.1548
 !            by Dominik Stoeckinger
 !-----------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: len
  Real(Dp), Intent(in) :: T, GY(len)
  Real(Dp), Intent(out) :: F(len)

  Integer :: j,i
  Real(Dp) :: gy2(6), sumI, beta2(3), gamma1, gamma2, q

  q = t

  gy2 = gy(1:6)**2

  If (TwoLoopRGE) Then
   Do i=1,3       ! gauge couplings two loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + b_2(i,j) * gy2(j) - a_2(i,j) * gy2(3+j)
    Enddo
    f(i) = oo16pi2 * gy(i) * gy2(i) * ( b_1(i) + oo16pi2 * sumI )
   Enddo

   beta2(1) = ( 13.5_dp* gy2(1) + 1.8_dp * gy2(2) + 1.2_dp * gy2(4)  &
   &         - 0.4_dp * gy2(5) ) * gy2(1)                            &
   &       + (7.5_dp * gy2(2) + 6._dp * gy2(4) ) * gy2(2 )           &
   &       + 16._dp * gy2(5) * gy2(3)                                &
   &       - (10._dp * gy2(4) + 9._dp * gy2(5) ) * gy2(4)            &
   &       - (9._dp * gy2(5) +  3._dp * gy2(6) ) * gy2(5)

   beta2(2) = ( 287._dp * gy2(1) / 90._dp + gy2(2) + 8._dp * gy2(3) / 9._dp   &
   &         + 1.2_dp * gy2(4) + 0.4_dp * gy2(5) + 0.8_dp * gy2(6) ) * gy2(1) &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(5) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(5) ) * gy2(3)                    &
   &       - 3._dp * (gy2(4) + gy2(5) ) * gy2(4)                              &
   &       - (22._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)                     &
   &       - 5._dp * gy2(6)**2

   beta2(3) = ( 2743._dp * gy2(1) / 450._dp + gy2(2) +136._dp * gy2(3)/45._dp &
   &         + 0.4_dp * gy2(5) + 1.2_dp * gy2(6) ) * gy2(1)                   &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(6) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(6) ) * gy2(3)                    &
   &       - (gy2(4) + 5._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)             &
   &       - 22._dp * gy2(6)**2

   Do i=1,3       ! yukawa couplings two loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    f(i+3) = oo16pi2 * gy(i+3) * (sumI + oo16pi2 * beta2(i) )
   End Do
  !---------------
  ! Ln(tan(beta))
  !---------------
   gamma1 = 3._dp * (gy2(5) - gy2(6)) + gy2(4)

   gamma2 = 0.75_dp * ( 3._dp * (gy2(6)**2 - gy2(5)**2) - gy2(4)**2)      &
   &    - (1.9_dp * gy2(1) + 4.5_dp * gy2(2) + 20._dp * gy2(3) ) * gy2(6) &
   &   + (0.4_dp * gy2(1) + 4.5_dp * gy2(2) + 20._dp * gy2(3) ) * gy2(5)  &
   &   + (1.8_dp * gy2(1) + 1.5_dp * gy2(2) ) * gy2(4)                    &
   &   + (0.3_dp * gy2(1) + 1.5_dp * gy2(2) ) * gamma1 ! gauge dependence
   f(7) = oo16pi2 * (gamma1 + oo16pi2 * gamma2 ) 

  Else ! Everything at 1-loop
   f(1:3) = oo16pi2 * gy(1:3) * gy2(1:3) * b_1   ! gauge couplings
   Do i=1,3       ! yukawa couplings one loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    f(i+3) = oo16pi2 * gy(i+3) * sumI 
   End Do
   f(7) = oo16pi2 * ( 3._dp * (gy2(5) - gy2(6)) + gy2(4) )
  End If

  End Subroutine rge7


 Subroutine rge8_NMSSM(len, T,GY,F)
 !-----------------------------------------------------------------
 !     Right hand side of truncated renormalization group equations
 !          dGY_i/dT = F_i(G)
 !     for the determination of M_GUT and the value of alpha_GUT
 !     and values of the Yukawas
 !  written by Werner Porod, 28.12.99
 !  6.10.2000: changing to f90
 !  4.12.2008: extension to include NMSSM couplings
 !-----------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: len
  Real(Dp), Intent(in) :: T,GY(len)
  Real(Dp), Intent(out) :: F(len)

  Integer :: j,i
  Real(Dp) :: gy2(len),sumI,beta2(3), q

  q = t

  gy2 = gy**2

  Do i=1,3       ! gauge couplings two loop
   sumI = 0._dp
   Do j=1,3
    sumI = sumI + b_2(i,j) * gy2(j) - a_2(i,j) * gy2(3+j)
   Enddo
   f(i) = oo16pi2 * gy(i) * gy2(i) * ( b_1(i) + oo16pi2 * sumI )
  Enddo

  beta2(1) = ( 13.5_dp* gy2(1) + 1.8_dp * gy2(2) + 1.2_dp * gy2(4)   &
   &         - 0.4_dp * gy2(5) ) * gy2(1)                            &
   &       + (7.5_dp * gy2(2) + 6._dp * gy2(4) ) * gy2(2)            &
   &       + 16._dp * gy2(5) * gy2(3)                                &
   &       - (10._dp * gy2(4) + 9._dp * gy2(5) ) * gy2(4)            &
   &       - (9._dp * gy2(5) +  3._dp * gy2(6) ) * gy2(5)

  beta2(2) = ( 287._dp * gy2(1) / 90._dp + gy2(2) + 8._dp * gy2(3) / 9._dp    &
   &         + 1.2_dp * gy2(4) + 0.4_dp * gy2(5) + 0.8_dp * gy2(6) ) * gy2(1) &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(5) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(5) ) * gy2(3)                    &
   &       - 3._dp * (gy2(4) + gy2(5) ) * gy2(4)                              &
   &       - (22._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)                     &
   &       - 5._dp * gy2(6)**2

  beta2(3) = ( 2743._dp * gy2(1) / 450._dp + gy2(2) + 136._dp * gy2(3)/45._dp &
   &         + 0.4_dp * gy2(5) + 1.2_dp * gy2(6) ) * gy2(1)                   &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(6) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(6) ) * gy2(3)                    &
   &       - (gy2(4) + 5._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)             &
   &       - 22._dp * gy2(6)**2

   Do i=1,3       ! yukawa couplings two loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    sumI = sumI + gy2(7) ! adding lambda**2
    f(i+3) = oo16pi2 * gy(i+3) * (sumI + oo16pi2 * beta2(i) )
   End Do
   f(7) = oo16pi2 * gy(7) * (-gy2(1) - 3._dp * gy2(2) + gy2(4)         &
        &                   + 3._dp * (gy2(5)+gy2(6)) + 4._dp * gy2(7) &
        &                   + 2._dp* gy2(8) )
   f(8) = oo16pi2 * gy(8) * (  gy2(7) +  gy2(8) )

  End Subroutine rge8_NMSSM


 Subroutine rge9_NMSSM(len, T,GY,F)
 !-----------------------------------------------------------------
 !     Right hand side of truncated renormalization group equations
 !          dGY_i/dT = F_i(G)
 !     for the determination of M_GUT and the value of alpha_GUT
 !     and values of the Yukawas
 !  written by Werner Porod, 28.12.99
 !  10.01.00: including tan(beta)
 !  6.10.2000: changing to f90
 !  4.12.2008: extension to include NMSSM couplings
 !  27.07.13: adding gauge dependence as discussed in 1305.1548
 !            by Dominik Stoeckinger
 !-----------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: len
  Real(Dp), Intent(in) :: T, GY(len)
  Real(Dp), Intent(out) :: F(len)

  Integer :: j,i
  Real(Dp) :: gy2(8), sumI, beta2(3), gamma1, gamma2, q

  q = t

  gy2 = gy(1:8)**2

  If (TwoLoopRGE) Then
   Do i=1,3       ! gauge couplings two loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + b_2(i,j) * gy2(j) - a_2(i,j) * gy2(3+j)
    Enddo
    f(i) = oo16pi2 * gy(i) * gy2(i) * ( b_1(i) + oo16pi2 * sumI )
   Enddo

   beta2(1) = ( 13.5_dp* gy2(1) + 1.8_dp * gy2(2) + 1.2_dp * gy2(4)  &
   &         - 0.4_dp * gy2(5) ) * gy2(1)                            &
   &       + (7.5_dp * gy2(2) + 6._dp * gy2(4) ) * gy2(2 )           &
   &       + 16._dp * gy2(5) * gy2(3)                                &
   &       - (10._dp * gy2(4) + 9._dp * gy2(5) ) * gy2(4)            &
   &       - (9._dp * gy2(5) +  3._dp * gy2(6) ) * gy2(5)

   beta2(2) = ( 287._dp * gy2(1) / 90._dp + gy2(2) + 8._dp * gy2(3) / 9._dp   &
   &         + 1.2_dp * gy2(4) + 0.4_dp * gy2(5) + 0.8_dp * gy2(6) ) * gy2(1) &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(5) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(5) ) * gy2(3)                    &
   &       - 3._dp * (gy2(4) + gy2(5) ) * gy2(4)                              &
   &       - (22._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)                     &
   &       - 5._dp * gy2(6)**2

   beta2(3) = ( 2743._dp * gy2(1) / 450._dp + gy2(2) +136._dp * gy2(3)/45._dp &
   &         + 0.4_dp * gy2(5) + 1.2_dp * gy2(6) ) * gy2(1)                   &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(6) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(6) ) * gy2(3)                    &
   &       - (gy2(4) + 5._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)             &
   &       - 22._dp * gy2(6)**2

   Do i=1,3       ! yukawa couplings two loop
    sumI = gy2(7)
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    f(i+3) = oo16pi2 * gy(i+3) * (sumI + oo16pi2 * beta2(i) )
   End Do
  !---------------
  ! Ln(tan(beta))
  !---------------
   gamma1 = 3._dp * (gy2(5) - gy2(6)) + gy2(4)

   gamma2 = 0.75_dp * ( 3._dp * (gy2(6)**2 - gy2(5)**2) - gy2(4)**2)      &
   &    - (1.9_dp * gy2(1) + 4.5_dp * gy2(2) + 20._dp * gy2(3) ) * gy2(6) &
   &   + (0.4_dp * gy2(1) + 4.5_dp * gy2(2) + 20._dp * gy2(3) ) * gy2(5)  &
   &   + (1.8_dp * gy2(1) + 1.5_dp * gy2(2) ) * gy2(4)                    &
   &   + (0.3_dp * gy2(1) + 1.5_dp * gy2(2) ) * gamma1 ! gauge dependence

   f(7) = oo16pi2 * gy(7) * (-gy2(1) - 3._dp * gy2(2) + gy2(4)         &
        &                   + 3._dp * (gy2(5)+gy2(6)) + 4._dp * gy2(7) &
        &                   + 2._dp* gy2(8) )
   f(8) = oo16pi2 * gy(8) * (  gy2(7) +  gy2(8) )
   f(9) = oo16pi2 * (gamma1 + oo16pi2 * gamma2 )

  Else ! Everything at 1-loop

   f(1:3) = oo16pi2 * gy(1:3) * gy2(1:3) * b_1   ! gauge couplings
   Do i=1,3       ! yukawa couplings one loop
    sumI = gy2(7)
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    f(i+3) = oo16pi2 * gy(i+3) * sumI 
   End Do
   f(7) = oo16pi2 * gy(7) * (-gy2(1) - 3._dp * gy2(2) + gy2(4)         &
        &                   + 3._dp * (gy2(5)+gy2(6)) + 4._dp * gy2(7) &
        &                   + 2._dp* gy2(8) )
   f(8) = oo16pi2 * gy(8) * (  gy2(7) +  gy2(8) )
   f(9) = oo16pi2 * ( 3._dp * (gy2(5) - gy2(6)) + gy2(4) )
  End If

  End Subroutine rge9_NMSSM


 Subroutine RGE10_SM(len,t,gy,f)
 !--------------------------------------------------------
 ! RGEs within the SM assuming the MSbar scheme
 ! 2-loop RGEs for e
 ! 4-loop RGEs for g_3
 ! 2-loop RGEs for lepton masses
 ! 4-loop QCD and 2-loop QED RGES for quark masses
 ! Assumption: the only threhold to be checked is m_b
 ! input: t = Log(Q^2)
 !        gy(i) ... i=1  -> e(Q)
 !                  i=2  -> g_3
 !                  i=3  -> m_e
 !                  i=4  -> m_mu
 !                  i=5  -> m_tau
 !                  i=6  -> m_u
 !                  i=7  -> m_c
 !                  i=8  -> m_d
 !                  i=9  -> m_s
 !                  i=10 -> m_b, is optional
 ! output:
 !   f = d(gy)/d(t)
 ! written by Werner Porod, 03.12.03
 !--------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: t, gy(len)
  Real(dp), Intent(out) :: f(len)

  Integer :: i1
  Real(dp) :: g32, g34, g36, g38, e2, e4, g32e2, q
  Real(dp), Parameter :: b_e1(2) = (/ 76._dp / 9._dp , 80._dp / 9._dp /)    &
       & , b_e2(2) = (/ 460._dp / 27._dp , 464._dp / 27._dp /)              & 
       & , b_e3(2) = (/ 160._dp / 9._dp , 176._dp / 9._dp  /)               & 
       & , b_g1(2) = (/ -25._dp / 3._dp, -23._dp/3._dp /)                   &
       & , b_g2(2) = (/ -154._dp / 3._dp, -116._dp/3._dp /)                 &
       & , b_g3(2) = (/ 20._dp / 3._dp, 22._dp/3._dp /)                     &
       & , b_g4(2) = (/ -21943._dp/54._dp, 9769._dp/54._dp /)               &
       & , b_g5(2) = (/ -4918247._dp/1458._dp-414140._dp*zeta3/81._dp       &
       &             , 598391._dp/1458._dp - 352864._dp*zeta3/81._dp /)     &
       & , g_el1(2) = (/ -6._dp, -6._dp /)                                  &
       & , g_el2(2) = (/ 353._dp / 9._dp,  373._dp / 9._dp /)               & 
       & , g_eu1(2) = (/ -8._dp/3._dp, -8._dp/3._dp /)                      &
       & , g_eu2(2) = (/ 1472._dp / 81._dp, 1552._dp / 81._dp/)             & 
       & , g_eu3(2) = (/ -32._dp / 9._dp,  -32._dp / 9._dp/)                & 
       & , g_ed1(2) = (/ -2._dp/3._dp, -2._dp/3._dp /)                      &
       & , g_ed2(2) = (/ 377._dp / 81._dp,  397._dp / 81._dp /)             & 
       & , g_ed3(2) = (/ -8._dp / 9._dp,  -8._dp / 9._dp /)                 & 
       & , g_q1(2) = (/ - 8._dp , -8._dp /)                                 &
       & , g_q2(2) = (/ -1052._dp / 9._dp ,  -1012._dp / 9._dp /)           &
       & , g_q3(2) = (/ -144674._dp/81._dp + 1280._dp * zeta3 / 3._dp       &
       &              , -128858._dp/81._dp + 1600._dp * zeta3 / 3._dp /)    &
       & , g_q4(2) = (/ -7330357._dp/243._dp + 51584._dp* zeta3/3._dp       &
       &                - 16000._dp*zeta4 / 3._dp + 11200._dp* zeta5 /9._dp &
       &             , -1911065._dp/81._dp + 618400._dp* zeta3/27._dp       &
       &                - 18400._dp*zeta4 / 3._dp - 25600._dp* zeta5 /9._dp  /)

       
  Iname = Iname + 1
  NameOfUnit(Iname) = 'RGE10_SM'

  q = t

  If (len.Eq.9) Then ! check which beta function (anomalous dimension) to use
   i1 = 1
  Else If (len.Eq.10) Then
   i1 = 2
  Else
   Write(ErrCan,*) "Error in routine "//Trim(NameOfUnit(Iname))
   Write(ErrCan,*) "Length of the vector gy = ",len
   Call TerminateProgram
  End If

  g32 = gy(1)**2
  g34 = gy(1)**4
  g36 = gy(1)**6
  g38 = gy(1)**8
  e2 = gy(2)**2
  e4 = gy(2)**4
  g32e2 = g32 * e2 
 !--------
 ! g_3
 !--------
  f(1) = oo16pi2 * gy(1) * ( b_g1(i1)*g32                                     &
       &                   + oo16pi2 * ( b_g2(i1)*g34 + b_g3(i1)*g32e2        &
       &                               + oo16pi2 * ( b_g4(i1)*g36             &
       &                                           + oo16pi2 * b_g5(i1)*g38 )))
 !--------
 ! e
 !--------
  f(2) = oo16pi2 * gy(2) * ( b_e1(i1) * e2                                &
       &                   + oo16pi2 * (b_e2(i1) * e4 + b_e3(i1) * g32e2 ))
 !-----------------
 ! m_l, l=e,mu,tau
 !-----------------
  f(3:5) =  oo16pi2 * gy(3:5) * (g_el1(i1) * e2 + oo16pi2 *g_el2(i1) * e4)
 !---------
 ! m_u, m_c
 !---------
  f(6:7) = oo16pi2 * gy(6:7) * (g_eu1(i1) * e2 + g_q1(i1) * g32              &
         &                     + oo16pi2 * (g_eu2(i1)*e4 + g_eu3(i1) * g32e2 &
         &                                 + g_q2(i1) * g34                  &
         &                                 + oo16pi2 * (g_q3(i1) * g36       &
         &                                       + oo16pi2 * g_q4(i1) * g38 )))
 !---------------
 ! m_d, m_s, m_b
 !---------------
  f(8:len) = oo16pi2 * gy(8:len) * (g_ed1(i1) * e2 + g_q1(i1) * g32          &
         &                     + oo16pi2 * (g_ed2(i1)*e4 + g_ed3(i1) * g32e2 &
         &                                 + g_q2(i1) * g34                  &
         &                                 + oo16pi2 * (g_q3(i1) * g36       &
         &                                       + oo16pi2 * g_q4(i1) * g38 )))

  Iname = Iname - 1

 End Subroutine RGE10_SM


 Subroutine rge57(len,T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(3), Dgauge(3), q, TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)  &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)        &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)            &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)          &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                 &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(5,2)       &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge57'

  q = t

  Call GToCouplings(gy,gauge,Ye,Yd,Yu)

  gauge2 = gauge**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYdYd = Matmul(aYd,Yd)
  aYeYe = Matmul(aYe,Ye)
  aYuYu = Matmul(aYu,Yu)
  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aYdYd(i1,i1) = Real(aYdYd(i1,i1), dp)
   aYeYe(i1,i1) = Real(aYeYe(i1,i1), dp)
   aYuYu(i1,i1) = Real(aYuYu(i1,i1), dp)
  End Do

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYdYd),dp )
  TraceY(3) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = (3._dp,0._dp) * TraceY(2) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = (3._dp,0._dp) * TraceY(2) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)

  diagonal(3,1) = (3._dp,0._dp) * TraceY(3)              &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(3,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
   aYdYdaYuYu = Matmul(aYdYd,aYuYu)

   !------------------------------------------------
   ! these are hermitian matrices, clean up to
   ! avoid numerical problems
   !------------------------------------------------
   Do i1=1,3
    aYdYdaYdYd(i1,i1) = Real(aYdYdaYdYd(i1,i1), dp)
    aYeYeaYeYe(i1,i1) = Real(aYeYeaYeYe(i1,i1), dp)
    aYuYuaYuYu(i1,i1) = Real(aYuYuaYuYu(i1,i1), dp)
   End Do

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(2)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(2) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(3)
   hd(2) = 9._dp * TraceY(2) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(3)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(3) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(2) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2                                  &
        & * ( b_1 + oo16pi2 * ( Matmul(b_2,gauge2) - Matmul(a_2,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )

 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
  End If

  !--------------------------------------------
  ! This helps avoiding numerical instabilities
  !--------------------------------------------
!  Call Chop(DYe)
!  Call Chop(DYd)
!  Call Chop(DYu)

  Call CouplingsToG(Dgauge,DYe,DYd,DYu,f)

  Iname = Iname - 1

 End Subroutine rge57


 Subroutine rge58(len,T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 25.09.01: portation to f90
 ! 01.06.11: include tan(beta) evolution
 !  27.07.13: adding gauge dependence as discussed in 1305.1548
 !            by Dominik Stoeckinger
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(3), Dgauge(3), q, TraceY2(4)    &
    & , gamma1, gamma2
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)  &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)        &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)            &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)          &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                 &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(5,2)       &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge58'

  q = t

  Call GToCouplings(gy(1:57),gauge,Ye,Yd,Yu)

  gauge2 = gauge**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYdYd = Matmul(aYd,Yd)
  aYeYe = Matmul(aYe,Ye)
  aYuYu = Matmul(aYu,Yu)
  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aYdYd(i1,i1) = Real(aYdYd(i1,i1), dp)
   aYeYe(i1,i1) = Real(aYeYe(i1,i1), dp)
   aYuYu(i1,i1) = Real(aYuYu(i1,i1), dp)
  End Do

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYdYd),dp )
  TraceY(3) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = (3._dp,0._dp) * TraceY(2) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = (3._dp,0._dp) * TraceY(2) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)

  diagonal(3,1) = (3._dp,0._dp) * TraceY(3)              &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(3,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
   aYdYdaYuYu = Matmul(aYdYd,aYuYu)

   !------------------------------------------------
   ! these are hermitian matrices, clean up to
   ! avoid numerical problems
   !------------------------------------------------
   Do i1=1,3
    aYdYdaYdYd(i1,i1) = Real(aYdYdaYdYd(i1,i1), dp)
    aYeYeaYeYe(i1,i1) = Real(aYeYeaYeYe(i1,i1), dp)
    aYuYuaYuYu(i1,i1) = Real(aYuYuaYuYu(i1,i1), dp)
   End Do

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(2)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(2) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(3)
   hd(2) = 9._dp * TraceY(2) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(3)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(3) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(2) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

  !---------------
  ! Ln(tan(beta))
  !---------------
  gamma1 = 3._dp * (TraceY(2) - TraceY(3)) + TraceY(1)

  If (TwoLoopRGE) Then 
   gamma2 = 0.75_dp * ( 3._dp * (TraceY(3)**2 - TraceY(2)**2) - TraceY(1)**2)  &
   &   - (1.9_dp * gauge2(1) + 4.5_dp * gauge2(2) + 20._dp * gauge2(3) )       &
   &                                                              * TraceY(3)  &
   &   + (0.4_dp * gauge2(1) + 4.5_dp * gauge2(2) + 20._dp * gauge2(3) )       &
   &                                                              * TraceY(2)  &
   &   + (1.8_dp * gauge2(1) + 1.5_dp * gauge2(2) ) * TraceY(1)                &
   &   + (0.3_dp * gauge2(1) + 1.5_dp * gauge2(2) ) * gamma1 ! gauge dependence
  End If

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2                                  &
        & * ( b_1 + oo16pi2 * ( Matmul(b_2,gauge2) - Matmul(a_2,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )

 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
  End If


  Call CouplingsToG(Dgauge,DYe,DYd,DYu,f(1:57))
  If (TwoLoopRGE) Then 
   f(58) = oo16pi2 * (gamma1 + oo16pi2 * gamma2 )
  Else
   f(58) = oo16pi2 * gamma1
  End If

  Iname = Iname - 1

 End Subroutine rge58


 Subroutine rge75(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 04.03.2001: including neutrino Yukawas as given in the MSSM
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(4), Dgauge(3), q
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)  &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)        &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)            &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)          &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                 &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), TraceY2(4), diagonal(5,2)       &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)             &
    & ,Ynu(3,3), aYnu(3,3), aYnuYnu(3,3), sumnu1(3,3), betaYnu1(3,3), DYnu(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge75'

  q = t

  Call GToCouplings2(gy,gauge,Ye,Ynu,Yd,Yu)

  gauge2 = gauge**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ynu,aYnu)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul2(aYe,Ye,OnlyDiagonal)
  aYnuYnu = Matmul2(aYnu,Ynu,OnlyDiagonal)
  aYdYd = Matmul2(aYd,Yd,OnlyDiagonal)
  aYuYu = Matmul2(aYu,Yu,OnlyDiagonal)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYnuYnu),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe + aYnuYnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = (3._dp,0._dp) * TraceY(4) + TraceY(2)     &
            &   - 0.6_dp * gauge2(1) - 3._dp * gauge2(2)
  sumnu1 = 3._dp * aYnuYnu + aYeYe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do

  betaYnu1 = Matmul2(Ynu,sumnu1,OnlyDiagonal)

  diagonal(3,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)

  diagonal(4,1) = (3._dp,0._dp) * TraceY(4)  + TraceY(2)             &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul2(Yu,sumu1,OnlyDiagonal)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul2(aYdYd,aYdYd,OnlyDiagonal)
   aYeYeaYeYe = Matmul2(aYeYe,aYeYe,OnlyDiagonal)
   aYuYuaYuYu = Matmul2(aYuYu,aYuYu,OnlyDiagonal)
   aYuYuaYdYd = Matmul2(aYuYu,aYdYd,OnlyDiagonal)
   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)

   TraceY2(1) = cTrace(aYeYeaYeYe)
   TraceY2(2) = cTrace(aYdYdaYdYd)
   TraceY2(3) = cTrace(aYuYuaYuYu)
   TraceY2(4) = cTrace(aYdYdaYuYu)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul2(Ye,sume2,OnlyDiagonal)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul2(Yd,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul2(Yu,sumu2,OnlyDiagonal)
    
  End If 

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2                                &
        & * ( b_1 + oo16pi2 * ( Matmul(b_2,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYnu = oo16pi2 * betaYnu1
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )

 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYnu = oo16pi2 * betaYnu1
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
  End If

  Call CouplingsToG2(Dgauge,DYe,DYnu,DYd,DYu,f)

  Iname = Iname - 1

 End Subroutine rge75


 Subroutine rge79(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.01.06: including neutrino dim 5 operator
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(4), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2)           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)

  Complex(dp), Dimension(3,3) :: YT, aYT, aYTYT, betaYT1, sumT1, DYT
  Complex(dp) :: lam1, lam2, Dlam1, Dlam2, betalam11, betalam21
  Real(dp) :: lam12, lam22, b_1a(3), Q, b_2a(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge79'

  q = t

  Call GToCouplings4(gy, gauge, Ye, YT, Yd, Yu, lam1, lam2)

  gauge2 = gauge**2
!-----------------
! beta functions
!-----------------
  lam12 = Abs(lam1)**2
  lam22 = Abs(lam2)**2

  Call Adjungate(Yd,aYd)
  Call Adjungate(YT,aYT)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul(aYe,Ye)
  aYTYT = Matmul(aYT,YT)
  aYdYd = Matmul(aYd,Yd)
  aYuYu = Matmul(aYu,Yu)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYTYT),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = 3._dp * (TraceY(3) + lam12) + TraceY(1)    &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * (aYeYe + aYTYT)
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
  sumT1 = aYeYe + 6._dp * aYTYT
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaYT1 = Matmul(YT,sumT1)  &
        & + Matmul(Transpose(aYeYe),YT)

  diagonal(3,1) = 3._dp * (TraceY(3)  + lam12) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)

  diagonal(4,1) = 3._dp * (TraceY(4) + lam22)                &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  betalam11 = lam1 * (7._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
            &        + 6._dp * TraceY(3)                            &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  betalam21 = lam2 * (7._dp * lam22  + 6._dp * TraceY(4)     &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
   aYdYdaYuYu = Matmul(aYdYd,aYuYu)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

 !---------------
 ! 2-loop RGEs
 !---------------
  b_1a = b_1 + Delta_b_1
  b_2a = b_2 + Delta_b_2

  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1a + oo16pi2 * (Matmul(b_2a,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYT = oo16pi2 * betaYT1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1a 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYT = oo16pi2 * betaYT1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
  End If

  Call CouplingsToG4(Dgauge, DYe, DYT, DYd, DYu, Dlam1, Dlam2, f)

  Iname = Iname - 1

 End Subroutine rge79


 Subroutine rge93(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 
 ! 30.08.2002: taking rge75 as basis and adding the RGEs for the
 !             left-neutrino mass matrix
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(4), Dgauge(3), q
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)  &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)        &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)            &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)          &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                 &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), TraceY2(4), diagonal(5,2)       &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)             &
    & ,Ynu(3,3), aYnu(3,3), aYnuYnu(3,3), sumnu1(3,3), betaYnu1(3,3), DYnu(3,3)
  Complex(dp) :: Mnu(3,3), DMnu(3,3), betamnu1(3,3), sumM1(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge93'

  Call GToCouplings3(gy,gauge,Ye,Ynu,Yd,Yu,Mnu)
  
  gauge2 = gauge**2
  q = t
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ynu,aYnu)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul(aYe,Ye)
  aYnuYnu = Matmul(aYnu,Ynu)
  aYdYd = Matmul(aYd,Yd)
  aYuYu = Matmul(aYu,Yu)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYnuYnu),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe + aYnuYnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = (3._dp,0._dp) * TraceY(4) + TraceY(2)     &
            &   - 0.6_dp * gauge2(1) - 3._dp * gauge2(2)
  sumnu1 = 3._dp * aYnuYnu + aYeYe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do

  betaYnu1 = Matmul(Ynu,sumnu1)

  diagonal(3,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)

  diagonal(4,1) = (3._dp,0._dp) * TraceY(4)  + TraceY(2)             &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  !--------------------------------
  ! neutrino dim. 5 operator
  !--------------------------------
  sumM1 = aYeYe + aYnuYnu
  diagonal(5,1) = 2._dp * TraceY(2) + 6._dp * TraceY(4)   &
              & - 1.2_dp * gauge2(1) - 6._dp * gauge2(2)
  betaMnu1 = Matmul( Transpose(sumM1), Mnu) + Matmul(Mnu, sumM1)  &
          & + diagonal(5,1) * Mnu
  

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
   aYdYdaYuYu = Matmul(aYdYd,aYuYu)

   TraceY2(1) = cTrace(aYeYeaYeYe)
   TraceY2(2) = cTrace(aYdYdaYdYd)
   TraceY2(3) = cTrace(aYuYuaYuYu)
   TraceY2(4) = cTrace(aYdYdaYuYu)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2                                &
        & * ( b_1 + oo16pi2 * ( Matmul(b_2,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYnu = oo16pi2 * betaYnu1
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
 !----------------------
 ! left neutrino mass
 !----------------------
   DMnu = oo16pi2 * betaMnu1
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYnu = oo16pi2 * betaYnu1
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
 !----------------------
 ! left neutrino mass
 !----------------------
   DMnu = oo16pi2 * betaMnu1
  End If

  Call CouplingsToG3(Dgauge,DYe,DYnu,DYd,DYu,DMnu,f)

  Iname = Iname - 1

 End Subroutine rge93


 Subroutine rge118(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.01.06: including neutrino dim 5 operator
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(6), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2)           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)

  Complex(dp), Dimension(3,3) :: YT, aYT, aYTYT, betaYT1, sumT1, DYT 
  Complex(dp), Dimension(3,3) :: YZ, aYZ, aYZYZ, YZaYZ, betaYZ1, sumZ1, DYZ  &
    & , YdaYd
  Complex(dp), Dimension(3,3) :: YS, aYS, aYSYS, YSaYS, betaYS1, sumS1, DYS 
  Complex(dp) :: lam1, lam2, Dlam1, Dlam2, betalam11, betalam21
  Real(dp) :: lam12, lam22, b_1a(3), Q, M15(3), betaM15(3), DM15(3), b_2a(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge118'

  q = t

  Call GToCouplings5(gy, gauge, Ye, YT, Yd, Yu, YZ, YS, lam1, lam2, M15)

  gauge2 = gauge**2
!-----------------
! beta functions
!-----------------
  lam12 = Abs(lam1)**2
  lam22 = Abs(lam2)**2

  Call Adjungate(Yd,aYd)
  Call Adjungate(YT,aYT)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)
  Call Adjungate(YZ,aYZ)
  Call Adjungate(YS,aYS)

  aYeYe = Matmul(aYe,Ye)
  aYTYT = Matmul(aYT,YT)
  aYdYd = Matmul(aYd,Yd)
  aYuYu = Matmul(aYu,Yu)
  aYZYZ = Matmul(aYZ,YZ)
  aYSYS = Matmul(aYS,YS)

  YdaYd = Matmul(Yd,aYd)
  YZaYZ = Matmul(YZ,aYZ)
  YSaYS = Matmul(YS,aYS)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYTYT),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )
  TraceY(5) = Real( cTrace(aYZYZ),dp )
  TraceY(6) = Real( cTrace(aYSYS),dp )

  diagonal(1,1) = 3._dp * (TraceY(3) + lam12) + TraceY(1)    &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * (aYeYe + aYTYT + aYZYZ)
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
  sumT1 = aYeYe + 6._dp * aYTYT + 3._dp * aYZYZ
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaYT1 = Matmul(YT,sumT1)                 &
        & + Matmul(Transpose(aYeYe+3._dp * aYZYZ),YT)

  diagonal(3,1) = 3._dp * (TraceY(3)  + lam12) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)  &
        & + 2._dp * Matmul(YZaYZ + 2._dp * YSaYS, Yd)

  diagonal(4,1) = 3._dp * (TraceY(4) + lam22)                &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  diagonal(5,1) = TraceY(5) + c1_1(2,1) * gauge2(1) &
            &   + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumZ1 = aYeYe + 3._dp * aYTYT + 5._dp * aYZYZ
  Do i1=1,3
   sumZ1(i1,i1) = sumZ1(i1,i1) + diagonal(5,1)
  End Do

  betaYZ1 = Matmul(YZ,sumZ1)                 &
        & + 2._dp * Matmul(YdaYd+2._dp * aYSYS,YZ)

  diagonal(6,1) = TraceY(6) - 0.8_dp * gauge2(1) - 12._dp * gauge2(3)
  sumS1 = 2._dp * Transpose(YdaYd + YZaYZ) + 8._dp * aYSYS
  Do i1=1,3
   sumS1(i1,i1) = sumS1(i1,i1) + diagonal(6,1)
  End Do

  betaYS1 = Matmul(YS,sumS1)                 &
        & + 2._dp * Matmul(YdaYd + YZaYZ,YS)


  betalam11 = lam1 * (7._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
            &        + 6._dp * TraceY(3)                            &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  betalam21 = lam2 * (7._dp * lam22  + 6._dp * TraceY(4)     &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  !---------------------------------------------
  ! beta function for the 15-plet, MT, M_Z, M_S
  !---------------------------------------------
   betaM15(1) = M15(1) * ( TraceY(2) + lam12 + lam22 &
            &        - 2.4_dp * gauge2(1) - 8._dp * gauge2(2)  )
   betaM15(2) = M15(2) * ( TraceY(5) - gauge2(1)/15._dp - 3._dp * gauge2(2) &
                     & - 16._dp * gauge2(3) / 3._dp  )
   betaM15(3) = M15(3) * ( TraceY(6)  &
                     & - (3.2_dp * gauge2(1) + 40._dp * gauge2(3))/3._dp  )

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
   aYdYdaYuYu = Matmul(aYdYd,aYuYu)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

 !---------------
 ! 2-loop RGEs
 !---------------
  b_1a = b_1 + Delta_b_1
  b_2a = b_2 + Delta_b_2

  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1a + oo16pi2 * (Matmul(b_2a,gauge2) - Matmul(a_2a,TraceY(1:4))))
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYT = oo16pi2 * betaYT1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
   DYZ = oo16pi2 * betaYZ1 
   DYS = oo16pi2 * betaYS1 
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
  !----------
  ! 15-plet
  !----------
   DM15 = oo16pi2 * betaM15
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1a 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYT = oo16pi2 * betaYT1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
   DYZ = oo16pi2 * betaYZ1 
   DYS = oo16pi2 * betaYS1 
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
  !----------
  ! 15-plet
  !----------
   DM15 = oo16pi2 * betaM15
  End If

  Call CouplingsToG5(Dgauge, DYe, DYT, DYd, DYu, DYZ, DYS, Dlam1, Dlam2 &
       & , DM15, f)

  Iname = Iname - 1

 End Subroutine rge118


 Subroutine rge213(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2
  Real(dp) :: gauge(3), gauge2(3), sumI, TraceY(3), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)  &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)        &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)            &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)          &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                 &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(5,2)       &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3), hc(4)

  Complex(dp) :: Mhlf(3),DMhlf(3)

  Complex(dp) :: Ae(3,3), Ad(3,3), Au(3,3), aAe(3,3), aAd(3,3), aAu(3,3)   &
     &  , DAe(3,3), DAd(3,3), DAu(3,3), aAdAd(3,3), aAeAe(3,3), aAuAu(3,3) &
     &  , aYdAd(3,3), aYeAe(3,3), aYuAu(3,3), TraceaYA(3), betaAd1(3,3)    &
     &  , betaAd2(3,3), betaAe1(3,3), betaAe2(3,3), betaAu1(3,3), betaAu2(3,3)
  Real(dp) :: TraceA(3)
  Complex(dp) :: aYdYdaYdAd(3,3), aYdAdaYdYd(3,3), TraceAY2(5)               &
     &  , aYeYeaYeAe(3,3), aYeAeaYeYe(3,3), aYuYuaYuAu(3,3), aYuAuaYuYu(3,3) &
     &  , aYuYuaYdAd(3,3), aYuAuaYdYd(3,3), aYdYdaYuAu(3,3), aYdAdaYuYu(3,3)
   
  Complex(dp) :: Me(3,3), Ml(3,3), Md(3,3), Mq(3,3), Mu(3,3), DMe(3,3)        &
     & , DMl(3,3), DMd(3,3), DMq(3,3), DMu(3,3), YdaYd(3,3), YeaYe(3,3)       &
     & , YuaYu(3,3), MdYdaYd(3,3), MeYeaYe(3,3), MuYuaYu(3,3), YdaYdMd(3,3)   &
     & , YeaYeMe(3,3), YuaYuMu(3,3), YdMqaYd(3,3), YeMlaYe(3,3)               &
     & ,  YuMqaYu(3,3), AdaAd(3,3), AeaAe(3,3), AuaAu(3,3), betaMd1(3,3)      &
     & , betaMd2(3,3), betaMe1(3,3), betaMe2(3,3), betaMl1(3,3), betaMl2(3,3) &
     & , betaMq1(3,3), betaMq2(3,3), betaMu1(3,3), betaMu2(3,3), MqaYdYd(3,3) &
     & , MqaYuYu(3,3), aYdYdMq(3,3), aYuYuMq(3,3), aYeYeMl(3,3), MlaYeYe(3,3) &
     & , aYeMeYe(3,3), aYdMdYd(3,3), aYuMuYu(3,3)                             &
     & , YdaYdYdaYd(3,3), YeaYeYeaYe(3,3), YuaYuYuaYu(3,3), MeYeaYeYeaYe(3,3) &
     & , YeaYeYeaYeMe(3,3), YeaYeMeYeaYe(3,3), YeMlaYeYeaYe(3,3), AeaYe(3,3)  &
     & , YeaYeYeMlaYe(3,3), AeaAeYeaYe(3,3), YeaYeAeaAe(3,3), YeaAe(3,3)      &
     & , AeaYeYeaAe(3,3), YeaAeAeaYe(3,3), Tr3aAdYdaAeYe

  Complex(dp) :: Tr3aYdAdaYeAe, AdaYd(3,3), YdaAd(3,3), MlaYeYeaYeYe(3,3)     &
     & , aYeYeaYeYeMl(3,3), aYeYeMlaYeYe(3,3), aYeYeaYeMeYe(3,3)              &
     & , aYeMeYeaYeYe(3,3), aAdYd(3,3),aAeYe(3,3), aAeAeaYeYe(3,3)            &
     & , aYeYeaAeAe(3,3), aAeYeaYeAe(3,3), aYeAeaAeYe(3,3), MdYdaYdYdaYd(3,3) &
     & , YdaYdYdaYdMd(3,3), YdMqaYdYdaYd(3,3),YdaYdMdYdaYd(3,3)               &
     & , YdaYdYdMqaYd(3,3), AdaAdYdaYd(3,3), YdaYDAdaAd(3,3), AdaYdYdaAd(3,3) &
     & , YdaAdAdaYd(3,3)
  Complex(dp) :: MdYdaYuYuaYd(3,3), YdaYuYuaYdMd(3,3), YdMqaYuYuaYd(3,3)      &
     & , YdaYuYuMqaYd(3,3), YdaYuMuYuaYd(3,3), AdaAuYuaYd(3,3)                &
     & , YdaYuAuaAd(3,3), AdaYuYuaAd(3,3), YdaAuAuaYd(3,3), YdaYuYuaYd(3,3)   &
     & , Tr3aYuAu, Tr3aAuYu, YuaAu(3,3)                                       &
     & , MqaYdYdaYdYd(3,3), aYdYdaYdYDMq(3,3), aYdMdYdaYdYd(3,3)              &
     & , aYdYdMqaYdYd(3,3), aYdYdaYdMdYd(3,3), aAdAdaYdYd(3,3)                &
     & , aYdYDaAdAd(3,3), aAdYdaYdAd(3,3), aYdAdaAdYd(3,3), MqaYuYuaYuYu(3,3) &
     & , aYuYuaYuYUMq(3,3), aYuMuYuaYuYu(3,3), aYuYuMqaYuYu(3,3)              &
     & , aYuYuaYuMuYu(3,3), aAuAuaYuYu(3,3), aYuYUaAuAu(3,3), aAuYuaYuAu(3,3) &
     & , aYuAuaAuYu(3,3), aAuYu(3,3), AuaYu(3,3), YuaYdYdaYu(3,3)             &
     & , AuaYdYdaAu(3,3), YuaAdAdaYu(3,3), AuaAdYdaYu(3,3), YuaYdAdaAu(3,3)   &
     & , YuMqaYuYuaYu(3,3), YuaYuYuMqaYu(3,3), MuYuaYuYuaYu(3,3)              &
     & , YuaYuYuaYuMu(3,3), YuaYuMuYuaYu(3,3), AuaAuYuaYu(3,3)                &
     & , YuaYuAuaAu(3,3), AuaYuYuaAu(3,3), YuaAuAuaYu(3,3), MuYuaYdYdaYu(3,3) &
     & , YuaYdYdaYuMu(3,3), YuMqaYdYdaYu(3,3), YuaYdYdMqaYu(3,3)              &
     & , YuaYdMdYdaYu(3,3)

  Real(dp) :: S1, S2, sig(3), Tr3aYdYdaYeYe, Tr3aAdAdaAeAe, AbsGM2(3)         &
     & , Tr3MqaYdYd3aYDMdYd , Tr3MqaYuYu3aYuMu, Tr3aAuAu

  Real(dp) :: Mh(2), DMh(2), TraceMH1(3), TraceMH2(2), betaMH11        &
     &  , betaMH12, betaMH21, betaMH22, q

  Complex(dp) :: mue, B, Dmue, DB, TraceMue(2), TraceB(2), betaMue1, betaMue2 &
     & , betaB1, betaB2, g2Mi(3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge213'

  OnlyDiagonal = .Not.GenerationMixing

  q = t

  Call GToParameters(gy, gauge, Ye, Yd, Yu                            &
                  & , Mhlf, Ae, Ad, Au, Me, Ml, Md, Mq, Mu, Mh, mue, B)

  gauge2 = gauge**2
  AbsGM2 = gauge2 * Abs( Mhlf )**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYdYd = Matmul(aYd,Yd)
  aYeYe = Matmul(aYe,Ye)
  aYuYu = Matmul(aYu,Yu)

  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aYdYd(i1,i1) = Real(aYdYd(i1,i1), dp)
   aYeYe(i1,i1) = Real(aYeYe(i1,i1), dp)
   aYuYu(i1,i1) = Real(aYuYu(i1,i1), dp)
  End Do
  
  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYdYd),dp )
  TraceY(3) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = 3._dp * TraceY(2) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = 3._dp * TraceY(2) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)

  diagonal(3,1) = 3._dp * TraceY(3)              &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(3,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
!   aYdYdaYuYu = Matmul(aYdYd,aYuYu)
   !------------------------------------------------
   ! this are hermitian matrices, clean up to
   ! avoid numerical problems
   !------------------------------------------------
   Do i1=1,3
    aYdYdaYdYd(i1,i1) = Real(aYdYdaYdYd(i1,i1), dp)
    aYeYeaYeYe(i1,i1) = Real(aYeYeaYeYe(i1,i1), dp)
    aYuYuaYuYu(i1,i1) = Real(aYuYuaYuYu(i1,i1), dp)
   End Do

   Call Adjungate(aYuYuaYdYd, aYdYdaYuYu)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(2)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(2) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(3)
   hd(2) = 9._dp * TraceY(2) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(3)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(3) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(2) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ad,aAd)
  Call Adjungate(Ae,aAe)
  Call Adjungate(Au,aAu)

  aAdAd = Matmul(aAd,Ad)
  aAeAe = Matmul(aAe,Ae)
  aAuAu = Matmul(aAu,Au)

  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aAdAd(i1,i1) = Real(aAdAd(i1,i1), dp)
   aAeAe(i1,i1) = Real(aAeAe(i1,i1), dp)
   aAuAu(i1,i1) = Real(aAuAu(i1,i1), dp)
  End Do

  TraceA(1) = Real( cTrace(aAeAe),dp )
  TraceA(2) = Real( cTrace(aAdAd),dp )
  TraceA(3) = Real( cTrace(aAuAu),dp )

  aYdAd = Matmul(aYd,Ad)
  aYeAe = Matmul(aYe,Ae)
  aYuAu = Matmul(aYu,Au)

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYdAd) 
  TraceaYA(3) = cTrace(aYuAu) 

  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = Matmul(Ae,sume1)

  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(2) + TraceaYA(1)  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    )
  sume1 = 4._dp * aYeAe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + Matmul(Ye,sume1)

  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = Matmul(Ad,sumd1)
  
  diagonal(2,1) = 2._dp * ( 3._dp * TraceaYA(2) + TraceaYA(1) &
                &         - c1_1(2,1) * g2Mi(1) - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do
  betaAd1 = betaAd1 + Matmul(Yd,sumd1)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = Matmul(Au,sumu1)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3)              &
                &         - c1_1(3,1) * g2Mi(1) - c1_1(3,2) * g2Mi(2)   &
                &         - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(3,1)
  End Do
  betaAu1 = betaAu1 + Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdAd = Matmul(aYdYd,aYdAd)
   aYdAdaYdYd = Matmul(aYdAd,aYdYd)
   aYeYeaYeAe = Matmul(aYeYe,aYeAe)
   aYeAeaYeYe = Matmul(aYeAe,aYeYe)
   aYuYuaYuAu = Matmul(aYuYu,aYuAu)
   aYuAuaYuYu = Matmul(aYuAu,aYuYu)
   aYuAuaYdYd = Matmul(aYuAu,aYdYd)
   aYuYuaYdAd = Matmul(aYuYu,aYdAd)
   aYdAdaYuYu = Matmul(aYdAd,aYuYu)
   aYdYdaYuAu = Matmul(aYdYd,aYuAu)
   TraceAY2(1) = cTrace(aYeYeaYeAe)
   TraceAY2(2) = cTrace(aYdYdaYdAd)
   TraceAY2(3) = cTrace(aYuYuaYuAu)
   TraceAY2(4) = cTrace(aYuYuaYdAd)
   TraceAY2(5) = cTrace(aYdYdaYuAu)

  !--------------
  ! A_e
  !--------------
   hd(1) = 6._dp * TraceY(2) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = Matmul(Ae,sume2)
    
   diagonal(1,2) = -6._dp * ( 6._dp * TraceAY2(2) + TraceAY2(4)       &
     &                      + TraceAY2(5) + 2._dp * TraceAY2(1)  )    &
     &  + ( 32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(2)   &
     &  - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(2)                &
     &  + 2.4_dp * gauge2(1) * TraceaYA(1) - 2.4_dp * g2Mi(1) * TraceY(1)    &
     &  - ( 30._dp * g2Mi(2)                              &
     &    + 3.6_dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)      &
     &  - 54._dp * gauge2(1)**2 * Mhlf(1)
   hd(1) = 12._dp * TraceY(2) + 4._dp * TraceY(1)     &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(1) = 18._dp * TraceaYA(2) + 6._dp * TraceaYA(1) &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2)
   sume2 = - 6._dp * aYeYeaYeAe - 8._dp * aYeAeaYeYe  &
         & - hd(1) * aYeAe - hc(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
   betaAe2 = betaAe2 + Matmul(Ye,sume2)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(2) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = Matmul(Ad,sumd2)
    
   diagonal(2,2) = diagonal(1,2)                                   &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                 &
     &              - gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp      &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)    &
     &  + 1.6_dp * gauge2(1) * gauge2(2) * (Mhlf(1)+Mhlf(2))       & 
     &  + 1.856e3_dp * gauge2(1)**2 * Mhlf(1) / 4.5e1_dp

   hd(1) = 12._dp * TraceY(2) + 4._dp * TraceY(1)  &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hd(2) = 6._dp * TraceY(3) - 1.6_dp * gauge2(1) 
   hc(1) = 18._dp * TraceaYA(2) + 6._dp * TraceaYA(1)                 &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2) + 1.6_dp * gauge2(1) *  Mhlf(1)
   hc(2) = 6._dp * TraceaYA(3) + 1.6_dp * gauge2(1) *  Mhlf(1)
   sumd2 = - 6._dp * aYdYdaYdAd - 8._dp * aYdAdaYdYd                  &
       &   - 4._dp * ( aYuAuaYuYu + aYuYuaYuAu + aYuAuaYdYd )         &
       &   - 2._dp * aYuYuaYdAd - hd(1) * aYdAd - hc(1) * aYdYd       &
       &  - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
   betaAd2 = betaAd2 + Matmul(Yd,sumd2)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(3) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = Matmul(Au,sumu2)
    
   diagonal(3,2) =  -6._dp * ( 6._dp * TraceAY2(3) + TraceAY2(4)        &
     &                        + TraceAY2(5)  )                          &
     &  + ( 32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(3)     &
     &  - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(3)           &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                      &
     &              - 3.4_dp * gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp  &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)         &
     &  - ( 30._dp * g2Mi(2)                                &
     &    + 2._dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)         &
     &  - 5486._dp * gauge2(1)**2 * Mhlf(1) / 225._dp
   hd(1) = 6._dp * TraceY(2) + 2._dp * TraceY(1) - 0.8_dp * gauge2(1)
   hc(1) = 6._dp * TraceaYA(2) + 2._dp * TraceaYA(1)   &
       &  + 0.8_dp * gauge2(1) *  Mhlf(1)
   hd(2) = 12._dp * TraceY(3) - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(2) = 18._dp * TraceaYA(3) + 1.2e1_dp * g2Mi(2) + 0.8_dp * g2Mi(1)
   sumu2 = - 6._dp * aYuYuaYuAu - 8._dp * aYuAuaYuYu                     &
       &   - 4._dp * ( aYdAdaYdYd + aYdYdaYdAd + aYdAdaYuYu )            &
       &   - 2._dp * aYdYdaYuAu - hd(1) * aYdAd - hc(1) * aYdYd          &
       &   - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
   betaAu2 = betaAu2 + Matmul(Yu,sumu2)

  End If 
!----------------------------------------------
! beta functions for Sfermion mass parameters
!----------------------------------------------
   S1 = mH(2) - mH(1)
   Do i1=1,3
    S1 = S1 + Real(Me(i1,i1),dp) - Real(Ml(i1,i1),dp) &
       &    + Real(Md(i1,i1),dp) + Real(Mq(i1,i1),dp) &
       &    - 2._dp * Real(Mu(i1,i1),dp)
   End Do
   S1 = S1 * gauge2(1)

   YdaYd = Matmul(Yd,aYd)
   YeaYe = Matmul(Ye,aYe)
   YuaYu = Matmul(Yu,aYu)

   MdYdaYd = Matmul(Md,YdaYd)
   MeYeaYe = Matmul(Me,YeaYe)
   MlaYeYe = Matmul(Ml,aYeYe)
   MqaYdYd = Matmul(Mq,aYdYd)
   MqaYuYu = Matmul(Mq,aYuYu)
   MuYuaYu = Matmul(Mu,YuaYu)

   Call Adjungate(MdYdaYd, YdaYdMd) ! YdaYdMd = Matmul(YdaYd,Md)
   Call Adjungate(MeYeaYe, YeaYeMe) ! YeaYeMe = Matmul(YeaYe,Me)
   Call Adjungate(MlaYeYe, aYeYeMl) ! aYeYeMl = Matmul(aYeYe,Ml)
   Call Adjungate(MqaYdYd, aYdYdMq) ! aYdYdMq = Matmul(aYdYd,Mq)
   Call Adjungate(MqaYuYu, aYuYuMq) ! aYuYuMq = Matmul(aYuYu,Mq)
   Call Adjungate(MuYuaYu, YuaYuMu) ! YuaYuMu = Matmul(YuaYu,Mu)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)

   AdaAd = Matmul(Ad,aAd)
   AeaAe = Matmul(Ae,aAe)
   AuaAu = Matmul(Au,aAu)

   !------------------------------------------------
   ! these are hermitian matrices, clean up to
   ! avoid numerical problems
   !------------------------------------------------
   Do i1=1,3
    aYdMdYd(i1,i1) = Real(aYdMdYd(i1,i1), dp)
    aYeMeYe(i1,i1) = Real(aYeMeYe(i1,i1), dp)
    aYuMuYu(i1,i1) = Real(aYuMuYu(i1,i1), dp)
    YdMqaYd(i1,i1) = Real(YdMqaYd(i1,i1), dp)
    YeMlaYe(i1,i1) = Real(YeMlaYe(i1,i1), dp)
    YuMqaYu(i1,i1) = Real(YuMqaYu(i1,i1), dp)
    AdaAd(i1,i1) = Real(AdaAd(i1,i1), dp)
    AeaAe(i1,i1) = Real(AeaAe(i1,i1), dp)
    AuaAu(i1,i1) = Real(AuaAu(i1,i1), dp)
   End Do

   diagonal(1,1) = - 4.8_dp * AbsGM2(1) + 1.2_dp * S1
   betaMe1 = 2._dp * (MeYeaYe + YeaYeMe)             &
         & + 4._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe )
   Do i1=1,3
    betaMe1(i1,i1) = betaMe1(i1,i1) + diagonal(1,1)
   End Do

   diagonal(2,1) = - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1
   betaMl1 = MlaYeYe + aYeYeMl + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe )
   Do i1=1,3
    betaMl1(i1,i1) = betaMl1(i1,i1) + diagonal(2,1)
   End Do

   diagonal(3,1) = - ( 32._dp * AbsGM2(3) + 1.6_dp * AbsGM2(1) ) / 3._dp &
               & + 0.4_dp * S1
   betaMd1 = 2._dp * (MdYdaYd + YdaYdMd)             &
         & + 4._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd )
   Do i1=1,3
    betaMd1(i1,i1) = betaMd1(i1,i1) + diagonal(3,1)
   End Do

   diagonal(4,1) = - ( 32._dp * AbsGM2(3) + 0.4_dp * AbsGM2(1) ) / 3._dp &
               & - 6._dp * AbsGM2(2) + 0.2_dp * S1
   betaMq1 = MqaYuYu + aYuYuMq + MqaYdYd + aYdYdMq             &
         & + 2._dp * ( mH(2) * aYuYu + mH(1) * aYdYd + aYuMuYu     &
         &           + aYdMdYd + aAuAu + aAdAd )
   Do i1=1,3
    betaMq1(i1,i1) = betaMq1(i1,i1) + diagonal(4,1)
   End Do

   diagonal(5,1) = - ( 32._dp * AbsGM2(3) + 6.4_dp * AbsGM2(1) ) / 3._dp &
               &   - 0.8_dp * S1
   betaMu1 = 2._dp * (MuYuaYu + YuaYuMu)             &
         & + 4._dp * ( mH(2) * YuaYu + YuMqaYu + AuaAu )
   Do i1=1,3
    betaMu1(i1,i1) = betaMu1(i1,i1) + diagonal(5,1)
   End Do

   If (TwoLoopRGE) Then
    YdaYdYdaYd = MatSquare(YdaYd,OnlyDiagonal)
    YeaYeYeaYe = MatSquare(YeaYe,OnlyDiagonal)
    YuaYuYuaYu = MatSquare(YuaYu,OnlyDiagonal)

    !------------------------------------------------
    ! these are hermitian matrices, clean up to
    ! avoid numerical problems
    !------------------------------------------------
    Do i1=1,3
     YdaYdYdaYd(i1,i1) = Real( YdaYdYdaYd(i1,i1), dp)
     YeaYeYeaYe(i1,i1) = Real( YeaYeYeaYe(i1,i1), dp)
     YuaYuYuaYu(i1,i1) = Real( YuaYuYuaYu(i1,i1), dp)
    End Do

    AdaYd = Matmul(Ad,aYd)
    AeaYe = Matmul(Ae,aYe)
    AuaYu = Matmul(Au,aYu)

    aAdYd = Matmul(aAd,Yd)
    aAeYe = Matmul(aAe,Ye)
    aAuYu = Matmul(aAu,Yu)

    Call Adjungate(AdaYd,YdaAd) ! YdaAd = Matmul(Yd,aAd)
    Call Adjungate(AeaYe,YeaAe) ! YeaAe = Matmul(Ye,aAe)
    Call Adjungate(AuaYu,YuaAu) ! YuaAu = Matmul(Yu,aAu)

    YdaYuYuaYd = MatMul3(Yd,aYuYu,aYd,OnlyDiagonal)
    AdaYuYuaAd = MatMul3(Ad,aYuYu,aAd,OnlyDiagonal)
    YdaAuAuaYd = MatMul3(Yd,aAuAu,aYd,OnlyDiagonal)
    AdaAuYuaYd = MatMul4(Ad,aAu,Yu,aYd,OnlyDiagonal)
    YdaYuAuaAd = MatMul3(Yd,aYuAu,aAd,OnlyDiagonal)

    YuaYdYdaYu = MatMul3(Yu,aYdYd,aYu,OnlyDiagonal)
    AuaYdYdaAu = MatMul3(Au,aYdYd,aAu,OnlyDiagonal)
    YuaAdAdaYu = MatMul3(Yu,aAdAd,aYu,OnlyDiagonal)
    AuaAdYdaYu = MatMul4(Au,aAd,Yd,aYu,OnlyDiagonal)
    YuaYdAdaAu = MatMul3(Yu,aYdAd,aAu,OnlyDiagonal)

    MdYdaYuYuaYd = Matmul(Md,YdaYuYuaYd)
    Call Adjungate(MdYdaYuYuaYd, YdaYuYuaYdMd)
    YdMqaYuYuaYd = MatMul3(Yd,MqaYuYu,aYd,OnlyDiagonal)
    Call Adjungate(YdMqaYuYuaYd, YdaYuYuMqaYd)
    YdaYuMuYuaYd = MatMul3(Yd,aYuMuYu,aYd,OnlyDiagonal)

    MuYuaYdYdaYu = Matmul(Mu,YuaYdYdaYu)
    Call Adjungate(MuYuaYdYdaYu, YuaYdYdaYuMu)
    YuMqaYdYdaYu = MatMul3(Yu,MqaYdYd,aYu,OnlyDiagonal)
    Call Adjungate(YuMqaYdYdaYu, YuaYdYdMqaYu)
    YuaYdMdYdaYu = MatMul3(Yu,aYdMdYd,aYu,OnlyDiagonal)

    MeYeaYeYeaYe = Matmul(MeYeaYe,YeaYe)
    Call Adjungate(MeYeaYeYeaYe,YeaYeYeaYeMe)
    aYeMeYeaYeYe = Matmul(aYeMeYe,aYeYe)
    Call Adjungate(aYeMeYeaYeYe,aYeYeaYeMeYe)
    YeaYeMeYeaYe = Matmul(YeaYeMe,YeaYe)

    MlaYeYeaYeYe = Matmul(MlaYeYe,aYeYe)
    Call Adjungate(MlaYeYeaYeYe, aYeYeaYeYeMl)
    YeMlaYeYeaYe = Matmul(YeMlaYe,YeaYe)
    Call Adjungate(YeMlaYeYeaYe, YeaYeYeMlaYe)
    aYeYeMlaYeYe = Matmul(aYeYeMl,aYeYe)

    MdYdaYdYdaYd = Matmul(MdYdaYd,YdaYd)
    Call Adjungate(MdYdaYdYdaYd, YdaYdYdaYdMd)
    aYdMdYdaYdYd = Matmul(aYdMdYd,aYdYd)
    Call Adjungate(aYdMdYdaYdYd, aYdYdaYdMdYd)
    YdaYdMdYdaYd = Matmul(YdaYdMd,YdaYd)

    MqaYdYdaYdYd = Matmul(MqaYdYd,aYdYd)
    Call Adjungate(MqaYdYdaYdYd, aYdYdaYdYdMq)
    YdMqaYdYdaYd = Matmul(YdMqaYd,YdaYd)
    Call Adjungate(YdMqaYdYdaYd, YdaYdYdMqaYd)
    aYdYdMqaYdYd = Matmul(aYdYdMq,aYdYd)

    MqaYuYuaYuYu = Matmul(MqaYuYu,aYuYu)
    Call Adjungate(MqaYuYuaYuYu, aYuYuaYuYuMq)
    YuMqaYuYuaYu = Matmul(YuMqaYu,YuaYu)
    Call Adjungate(YuMqaYuYuaYu, YuaYuYuMqaYu)
    aYuYuMqaYuYu = Matmul(aYuYuMq,aYuYu)

    MuYuaYuYuaYu = Matmul(MuYuaYu,YuaYu)
    Call Adjungate(MuYuaYuYuaYu, YuaYuYuaYuMu)
    aYuMuYuaYuYu = Matmul(aYuMuYu,aYuYu)
    Call Adjungate(aYuMuYuaYuYu, aYuYuaYuMuYu)
    YuaYuMuYuaYu = Matmul(YuaYuMu,YuaYu)

    AdaAdYdaYd = Matmul(AdaAd,YdaYd)
    Call Adjungate(AdaAdYdaYd, YdaYdAdaAd)
    AdaYdYdaAd = Matmul(AdaYd,YdaAd)
    YdaAdAdaYd = Matmul(YdaAd,AdaYd)

    aAdAdaYdYd = Matmul(aAdAd,aYdYd)
    Call Adjungate(aAdAdaYdYd, aYdYdaAdAd)
    aAdYdaYdAd = Matmul(aAdYd,aYdAd)
    aYdAdaAdYd = Matmul(aYdAd,aAdYd)

    AeaAeYeaYe = Matmul(AeaAe,YeaYe)
    Call Adjungate(AeaAeYeaYe, YeaYeAeaAe)
    AeaYeYeaAe = Matmul(AeaYe,YeaAe)
    YeaAeAeaYe = Matmul(YeaAe,AeaYe)

    aAeAeaYeYe = Matmul(aAeAe,aYeYe)
    Call Adjungate(aAeAeaYeYe, aYeYeaAeAe)
    aAeYeaYeAe = Matmul(aAeYe,aYeAe)
    aYeAeaAeYe = Matmul(aYeAe,aAeYe)

    AuaAuYuaYu = Matmul(AuaAu,YuaYu)
    Call Adjungate(AuaAuYuaYu, YuaYuAuaAu)
    AuaYuYuaAu = Matmul(AuaYu,YuaAu)
    YuaAuAuaYu = Matmul(YuaAu,AuaYu)

    aAuAuaYuYu = Matmul(aAuAu,aYuYu)
    Call Adjungate(aAuAuaYuYu, aYuYuaAuAu)
    aAuYuaYuAu = Matmul(aAuYu,aYuAu)
    aYuAuaAuYu = Matmul(aYuAu,aAuYu)

    S2 = (1.5_dp * gauge2(2) + 0.3_dp * gauge2(1) )            &
     &      * (MH(2) - MH(1) - Real(cTrace(ML),dp) )              &
     & + ( (8._dp * gauge2(3) + 0.1_dp*gauge2(1)) / 3._dp      &
     &   + 1.5_dp * gauge2(2) ) * Real( cTrace(Mq),dp   )         &
     & - (16._dp * gauge2(3) + 3.2_dp*gauge2(1) )              &
     &    * Real(cTrace(Mu),dp) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md),dp) / 3._dp                           &
     & + 1.2_dp*gauge2(1) * Real(cTrace(Me),dp)                   &
     & - 3._dp * (MH(2)*TraceY(3) - MH(1) * TraceY(2) )        &
     & + mH(1) * TraceY(1)
    
    Do i1=1,3
     S2 = S2 - Real(YuMqaYu(i1,i1),dp) + 4._dp * Real(aYuMuYu(i1,i1),dp)   &
        &    - Real(YdMqaYd(i1,i1),dp) - 2._dp * Real(aYdMdYd(i1,i1),dp)   &
        &    + Real(YeMlaYe(i1,i1),dp) - 2._dp * Real(aYeMeYe(i1,i1),dp)
    End Do

    sig(1) = 3._dp * (MH(1) + MH(2) + Real(cTrace(Ml),dp) )          &
         & + Real(cTrace(Mq),dp) + 8._dp * Real(cTrace(Mu),dp)          &
         & + 2._dp * Real(cTrace(Md),dp) + 6._dp * Real(cTrace(Me),dp) 
    sig(1) = 0.2_dp * gauge2(1) * sig(1)
    sig(2) = gauge2(2) * ( MH(1) + MH(2) + Real( cTrace(Ml),dp )    &
           &             + 3._dp * Real( cTrace(Mq),dp ) )
    sig(3) = gauge2(3) * ( 2._dp * Real( cTrace(Mq),dp )            &
           &             + Real(cTrace(Mu),dp) + Real(cTrace(Md),dp) )

    Tr3aYdYdaYeYe = 3._dp * TraceY(2) + TraceY(1)
    Tr3MqaYdYd3aYDMdYd = Real(cTrace(MlaYeYe),dp) + Real(cTrace(aYeMeYe),dp)  &
              & + 3._dp * ( Real(cTrace(MqaYdYd),dp)+Real(cTrace(aYdMdYd),dp) )
    Tr3aAdAdaAeAe = 3._dp * TraceA(2) + TraceA(1)
    Tr3aYdAdaYeAe = 3._dp * TraceaYA(2) + TraceaYA(1)
    Tr3aAdYdaAeYe = Conjg( Tr3aYdAdaYeAe )
    Tr3MqaYuYu3aYuMu = 3._dp * (Real(cTrace(MqaYuYu),dp) &
                     &         +Real(cTrace(aYuMuYu),dp) )
    Tr3aAuAu = 3._dp * TraceA(3)
    Tr3aYuAu = 3._dp * TraceaYA(3)
    Tr3aAuYu = Conjg( Tr3aYuAu )

    diagonal(1,2) = 2.4_dp * gauge2(1) * (S2 + sig(1) )   &
                & + 112.32_dp * gauge2(1) * AbsGM2(1)
    hd(1) = 6._dp*gauge2(2)-1.2_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 4.8_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 2.4_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMe2 = -2._dp * ( MeYeaYeYeaYe + YeaYeYeaYeMe )                       &
      & - 4._dp * ( YeMlaYeYeaYe + YeaYeMeYeaYe + YeaYeYeMlaYe )             &
      & - 8._dp * mH(1) * YeaYeYeaYe                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MeYeaYe + YeaYeMe + 4._dp * MH(1) * YeaYe &
      &                           + 2._dp * (YeMlaYe + AeaAE)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YeaYe + AeaAeYeaYe + YeaYEAeaAe     &
      &           + AeaYeYeaAe + YeaAeAeaYe + Tr3aAdAdaAeAe * YeaYe          &
      &           + Tr3aAdYdaAeYe * AeaYe + Tr3aYdAdaYeAe * YeaAe )          &
      & + hd(1) * ( MeYeaYe + YeaYeMe                                        &
      &           + 2._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe ) )            &
      & + hd(2) * YeaYe + hc(1) * AeaYe + hc(2) * YeaAe
    Do i1=1,3
     betaMe2(i1,i1) = betaMe2(i1,i1) + diagonal(1,2)
    End Do

    diagonal(2,2) = gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)        &
     &     + 3._dp * gauge2(2) * sig(2)                                 &
     &     + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)      &
     &                   + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2       &
     &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )     &
     &     + 24.84_dp * gauge2(1) * AbsGM2(1)
    betaMl2 = -2._dp * ( MlaYeYeaYeYe + aYeYeaYeYeMl )                   &
      & - 4._dp * ( aYeMeYeaYeYe + aYeYeMlaYeYe + aYeYeaYeMeYe )         &
      & - 8._dp * mH(1) * aYeYeaYeYe                                     &
      & - Tr3aYdYdaYeYe * ( MlaYeYe + aYeYeML + 4._dp * MH(1) * aYeYe    &
      &                   + 2._dp * (aYeMeYe + aAeAE)  )                 &
      & - 4._dp * (aAeAeaYeYe + aYeYeaAeAe + aAeYeaYeAe + aYeAeaAeYe )   &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYeYe +  Tr3aAdAdaAeAe * aYeYe  &
      &           + Tr3aAdYdaAeYe * aYeAe + Tr3aYdAdaYeAe * aAeYe     )  &
      & + 1.2_dp*gauge2(1) * ( MlaYeYe + aYeYeMl                         &
      &                      + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe &
      &                                - Mhlf(1) * aAeYe                 &
      &                                - Conjg(Mhlf(1)) * aYeAe ) )      &
      & + 4.8_dp * AbsGM2(1) * aYeYe
    Do i1=1,3
     betaMl2(i1,i1) = betaMl2(i1,i1) + diagonal(2,2)
    End Do

    diagonal(3,2) = 0.8_dp * gauge2(1) * (S2 + sig(1)/3._dp )             &
      &   + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp   &
      &   + 8.08e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                    &
      &   + 1.28e2_dp * ( gauge2(1) * AbsGM2(3)                           &
      &                 + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)     &
      &                   * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)+0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) + 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) - 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMd2 = -2._dp * ( MdYdaYdYdaYd + YdaYdYdaYdMd )                       &
      & - 4._dp * ( YdMqaYdYdaYd + YdaYdMdYdaYd + YdaYdYdMqaYd )             &
      & - 8._dp * mH(1) * YdaYdYdaYd                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MdYdaYd + YdaYdMd + 4._dp * MH(1) * YdaYd &
      &                           + 2._dp * (YdMqaYd + AdaAD)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YdaYd + AdaAdYdaYd + YdaYDAdaAd     &
      &           + AdaYdYdaAd + YdaAdAdaYd + Tr3aAdAdaAeAe * YdaYd          &
      &           + Tr3aAdYdaAeYe * AdaYd + Tr3aYdAdaYeAe * YdaAd )          &
      & + hd(1) * ( MdYdaYd + YdaYdMd                                        &
      &           + 2._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd ) )            &
      & + hd(2) * YdaYd + hc(1) * AdaYd + hc(2) * YdaAd                      &
      & - 2._dp * (MdYdaYuYuaYd + YdaYuYuaYdMd )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YdaYuYuaYd + YdMqaYuYuaYd+ YdaYuYuMqaYd  &
      &           + YdaYuMuYuaYd + AdaAuYuaYd + YdaYuAuaAd                   &
      &           + AdaYuYuaAd + YdaAuAuaYd )
    Do i1=1,3
     betaMd2(i1,i1) = betaMd2(i1,i1) + diagonal(3,2)
    End Do

    diagonal(4,2) = 0.2_dp * gauge2(1) * (2._dp * S2 + sig(1)/3._dp )     &
      &  + 3._dp * gauge2(2) * sig(2)                                     &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp    &
      &  + 1.99e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                     &
      &  + 32._dp * ( gauge2(1) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)         &
      &              * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp      &
      &  + 33._dp * gauge2(2) * AbsGM2(2)                                 &
      &  + 32._dp * ( gauge2(2) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(2)  + gauge2(3)*gauge2(2)        &
      &              * Real( Mhlf(2) * Conjg(Mhlf(3)),dp ) )                 &
      &  + 0.4_dp * ( gauge2(2) * AbsGM2(1)                               &
      &             + gauge2(1) * AbsGM2(2) + gauge2(1)*gauge2(2)         &
      &              * Real( Mhlf(2) * Conjg(Mhlf(1)),dp ) )
    hd(1) = 1.6_dp * AbsGM2(1)
    hc(1) = - 0.8_dp * g2Mi(1)
    hc(2) = Conjg( hc(1) )
    hd(2) = 3.2_dp * AbsGM2(1)
    hc(3) = - 1.6_dp * g2Mi(1)
    hc(4) = Conjg( hc(3) )
    betaMq2 = -2._dp * ( MqaYdYdaYdYd + aYdYdaYdYDMq )                       &
      & - 4._dp * ( aYdMdYdaYdYd + aYdYdMqaYdYd + aYdYdaYdMdYd )             &
      & - 8._dp * mH(1) * aYdYdaYdYd                                         &
      & - Tr3aYdYdaYeYe * ( MqaYdYd + aYdYdMq + 4._dp * MH(1) * aYdYd        &
      &                   + 2._dp * (aYdMdYd + aAdAd)  )                     &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYdYd                               &
      &           + 2._dp * ( aAdAdaYdYd + aYdYDaAdAd + aAdYdaYdAd           &
      &                     + aYdAdaAdYd )                                   &
      &           + Tr3aAdAdaAeAe * aYdYd + Tr3aAdYdaAeYe * aYdAd            &
      &           + Tr3aYdAdaYeAe * aAdYd )                                  &
      & + 0.4_dp*gauge2(1) * ( MqaYdYd + aYdYdMq                             &
      &                      + 2._dp * ( mH(1) * aYdYd + aYdMdYd + aAdAd ) ) &
      & + hd(1) * aYdYd + hc(1) * aAdYd + hc(2) * aYdAd                      &
      & - 2._dp * ( MqaYuYuaYuYu + aYuYuaYuYuMq )                            &
      & - 4._dp * ( aYuMuYuaYuYu + aYuYuMqaYuYu + aYuYuaYuMuYu )             &
      & - 8._dp * MH(2) * aYuYuaYuYu                                         &
      & - 3._dp * TraceY(3) * ( MqaYuYu + aYuYuMq + 4._dp * MH(2) * aYuYu    &
      &                       + 2._dp * (aYuMuYu + aAuAu) )                  &
      & - 2._dp * ( Tr3MqaYuYu3aYuMu * aYuYu                                 &
      &           + 2._dp * ( aAuAuaYuYu + aYuYuaAuAu                        &
      &                     + aAuYuaYuAu + aYuAuaAuYu )                      &
      &           + Tr3aAuAu * aYuYu + Tr3aAuYu * aYuAu + Tr3aYuAu * aAuYu ) &
      & + 0.8_dp*gauge2(1) * ( MqaYuYu + aYuYuMq                             &
      &                      + 2._dp * ( MH(2) * aYuYu + aYuMuYu + aAuAu ) ) &
      & + hd(2) * aYuYu + hc(3) * aAuYu + hc(4) * aYuAu
    Do i1=1,3
     betaMq2(i1,i1) = betaMq2(i1,i1) + diagonal(4,2)
    End Do

    diagonal(5,2) = 1.6_dp * gauge2(1) * (2._dp*sig(1)/3._dp - S2)           &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp       &
      &  + 3424._dp * gauge2(1) * AbsGM2(1) / 75._dp                         &
      &  + 512._dp * ( gauge2(1) * AbsGM2(3)                                 &
      &              + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)           &
      &                * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)-0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMu2 = -2._dp * ( MuYuaYuYuaYu + YuaYuYuaYuMu )                       &
      & - 4._dp * ( YuMqaYuYuaYu + YuaYuMuYuaYu + YuaYuYuMqaYu )             &
      & - 8._dp * MH(2) * YuaYuYuaYu                                         &
      & - 6._dp * TraceY(3) * ( MuYuaYu + YuaYuMu + 4._dp * MH(2) * YuaYu    &
      &                       + 2._dp * (YuMqaYu + AuaAu)  )                 &
      & - 4._dp * ( Tr3MqaYuYu3aYuMu * YuaYu + AuaAuYuaYu + YuaYuAuaAu       &
      &           + AuaYuYuaAu + YuaAuAuaYu + Tr3aAuAu * YuaYu               &
      &           + Tr3aAuYu * AuaYu + Tr3aYuAu * YuaAu )                    &
      & + hd(1) * ( MuYuaYu + YuaYuMu                                        &
      &           + 2._dp * ( MH(2) * YuaYu + YuMqaYu + AuaAu ) )            &
      & + hd(2) * YuaYu + hc(1) * AuaYu + hc(2) * YuaAu                      &
      & - 2._dp * (MuYuaYdYdaYu + YuaYdYdaYuMu )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YuaYdYdaYu + YuMqaYdYdaYu + YuaYdYdMqaYu &
      &           + YuaYdMdYdaYu + AuaAdYdaYu + YuaYdAdaAu                   &
      &           + AuaYdYdaAu + YuaAdAdaYu )
    Do i1=1,3
     betaMu2(i1,i1) = betaMu2(i1,i1) + diagonal(5,2)
    End Do

   End If 

  !------------------------------------------
  ! beta functions for Higgs mass parameters
  !------------------------------------------
   traceMH1(1) = mH(1) * TraceY(1) + Real( cTrace(YeMlaYe),dp ) &
             & + Real( cTrace(aYeMeYe),dp ) + TraceA(1)
   traceMH1(2) = mH(1) * TraceY(2) + Real( cTrace(YdMqaYd),dp ) &
             & + Real( cTrace(aYdMdYd),dp ) + TraceA(2)
   betamH11 = 6._dp * TraceMH1(2) + 2._dp * TraceMH1(1)      &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1

   traceMH2(1) = mH(2) * TraceY(3) + Real( cTrace(YuMqaYu),dp )  &
             & + Real( cTrace(aYuMuYu),dp ) + TraceA(3)

   betamH21 = 6._dp * TraceMH2(1) - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1)  &
          & + 0.6_dp * S1

   If (TwoLoopRGE) Then
    traceMH1(3) = MH(1) * (6._dp*TraceY2(2) + 2._dp*TraceY2(1) + TraceY2(4) ) &
              & + MH(2) * TraceY2(4)                                          &
              & + 6._dp * ( Real( cTrace(MqaYdYdaYdYd),dp )                   &
              &           + Real( cTrace(aYdMdYdaYdYd),dp )                   &
              &           + Real( cTrace(aAdAdaYdYd),dp )                     &
              &           + Real( cTrace(aAdYdaYdAd),dp )  )                  &
              & + 2._dp * ( Real( cTrace(MlaYeYeaYeYe),dp )                   &
              &           + Real( cTrace(aYeMeYeaYeYe),dp )                   &
              &           + Real( cTrace(aAeAeaYeYe),dp )                     &
              &           + Real( cTrace(aAeYeaYeAe),dp )  )                  &
           & + Real(cTrace(YdMqaYuYuaYd),dp) + Real(cTrace(YdaYuMuYuaYd),dp ) &
       & + Real( cTrace(YdaYuYuMqaYd),dp ) + Real( cTrace(YuaYdMdYdaYu),dp ) &
       & + Real( cTrace(YdaAuAuaYd),dp ) + Real( cTrace(AdaYuYuaAd),dp )     &
              & + Real( cTrace(AdaAuYuaYd),dp ) + Real( cTrace(YdaYuAuaAd),dp )
    betaMH12 = - 6._dp * traceMH1(3)                                       &
      &   + (32._dp*gauge2(3) - 0.8_dp*gauge2(1) ) * traceMH1(2)           &
      &   + 64._dp * ( AbsGM2(3) * TraceY(2)                               &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(2),dp ) )     &
      &   - 1.6_dp * ( AbsGM2(1) * TraceY(2)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(2),dp ) )     &
      &   + 2.4_dp*gauge2(1) * traceMH1(1)                                 &
      &   + 4.8_dp * ( AbsGM2(1) * TraceY(1)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(1),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)                   &
      &   + 3._dp * gauge2(2) * sig(2)                                     &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)          &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2           &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

    traceMH2(2) = MH(2) * ( 6._dp * TraceY2(3) + TraceY2(4) )                 &
      &  + MH(1) * TraceY2(4)                                                 &
      & + 6._dp*(Real(cTrace(MqaYuYuaYuYu),dp)+Real(cTrace(aYuMuYuaYuYu),dp)  &
      &        + Real(cTrace(aAuAuaYuYu),dp) + Real(cTrace(aAuYuaYuAu),dp))   &
      & + Real(cTrace(YuMqaYdYdaYu),dp)+Real(cTrace(YuaYdMdYdaYu),dp )        &
      & + Real(cTrace(YuaYdYdMqaYu),dp)+Real(cTrace(YdaYuMuYuaYd),dp )        &
      & + Real(cTrace(YuaAdAdaYu),dp) + Real( cTrace(AuaYdYdaAu),dp )         &
      &  + Real( cTrace(AuaAdYdaYu),dp ) + Real( cTrace(YuaYdAdaAu),dp )
    betaMH22 = - 6._dp * traceMH2(2)                                          &
      &   + (32._dp*gauge2(3) + 1.6_dp*gauge2(1) ) * traceMH2(1)              &
      &   + 64._dp * ( AbsGM2(3) * TraceY(3)                                  &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(3),dp ) )     &
      &   + 3.2_dp * ( AbsGM2(1) * TraceY(3)                                  &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(3),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) + 1.2_dp * S2)                      &
      &   + 3._dp * gauge2(2) * sig(2)                                        &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)             &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2              &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

   End If
!-----------------------------
! beta functions for mu and B
!-----------------------------
   TraceMue(1) = 3._dp * (TraceY(2)+TraceY(3)) + TraceY(1)  &
             & - 3._dp * gauge2(2) - 0.6_dp * gauge2(1)
   betaMue1 = mue * TraceMue(1)

   TraceB(1) = 6._dp * (TraceaYA(2)+TraceaYA(3)) + 2._dp * TraceaYA(1)    &
           & + 6._dp * g2Mi(2) + 1.2_dp * g2Mi(1)
   betaB1 = mue * TraceB(1) + B * TraceMue(1)

   If (TwoLoopRGE) Then
    TraceMue(2) = - 3._dp * ( 3._dp * (TraceY2(2) + TraceY2(3) )          &
      &                     + 2._dp * TraceY2(4) + TraceY2(1) )           &
      &         + (16._dp * gauge2(3) + 0.8_dp * gauge2(1)) * TraceY(3)   &
      &         + (16._dp * gauge2(3) - 0.4_dp * gauge2(1)) * TraceY(2)   &
      &         + 1.2_dp * gauge2(1) * TraceY(1)                          &
      &         + 7.5_dp * gauge2(2)**2                                   &
      &         + 1.8_dp * gauge2(2) * gauge2(1)                          &
      &         + 4.14_dp * gauge2(1)**2
    betaMue2 = mue * TraceMue(2)

    TraceB(1) = cTrace( 3._dp * (Matmul(AuaYu,YuaYu) + Matmul(AdaYd,YdaYd) )  &
              &       + Matmul(AeaYe,YeaYe) + Matmul(aYuAu,aYdYd)             &
              &       + Matmul(aYdAd,aYuYu) ) 
    TraceB(2) = -12._dp * TraceB(1)                                           &
      &   + (32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(3)          &
      &   + (32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(2)          &
      &   + 2.4_dp * gauge2(1) * TraceaYA(1)                                  &
      &   - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(3)               &
      &   - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(2)               &
      &   - 2.4_dp * g2Mi(1) * TraceY(1)                          &
      &   - 30._dp * gauge2(2)**2 * Mhlf(2)                                   &
      &   - 3.6_dp * gauge2(2) * gauge2(1) * (Mhlf(1) + Mhlf(2) )             &
      &   - 16.56_dp * gauge2(1)**2 * Mhlf(1) 
    betaB2 = mue * TraceB(2) + B * TraceMue(2)

   End If

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1 + oo16pi2 * (Matmul(b_2,gauge2) - Matmul(a_2,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   Do i1 = 1,3    
    sumI = 0._dp
    Do i2=1,3
     sumI = sumI + b_2(i1,i2) * gauge2(i2) * (Mhlf(i1) + Mhlf(i2) )     &
          &      + a_2(i1,i2) * ( TraceaYA(i2) - Mhlf(i1)*TraceY(i2) )
    End Do
    DMhlf(i1) = oo8pi2 * gauge2(i1) * ( b_1(i1) * Mhlf(i1) + oo16pi2 * sumI)
   End Do
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * ( betaAe1 + oo16pi2 * betaAe2 )
   DAd = oo16pi2 * ( betaAd1 + oo16pi2 * betaAd2 )
   DAu = oo16pi2 * ( betaAu1 + oo16pi2 * betaAu2 )
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * ( betaMe1 + oo16pi2 * betaMe2 )
   DMl = oo16pi2 * ( betaMl1 + oo16pi2 * betaMl2 )
   DMd = oo16pi2 * ( betaMd1 + oo16pi2 * betaMd2 )
   DMq = oo16pi2 * ( betaMq1 + oo16pi2 * betaMq2 )
   DMu = oo16pi2 * ( betaMu1 + oo16pi2 * betaMu2 )
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * ( betaMH11 + oo16pi2 * betaMH12 )
   DmH(2) = oo16pi2 * ( betaMH21 + oo16pi2 * betaMH22 )
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * ( betaMue1 + oo16pi2 * betaMue2 )
   DB = oo16pi2 * ( betaB1 + oo16pi2 * betaB2 )

 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   DMhlf = oo8pi2 * gauge2 * b_1 * Mhlf
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * betaAe1
   DAd = oo16pi2 * betaAd1
   DAu = oo16pi2 * betaAu1
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * betaMe1
   DMl = oo16pi2 * betaMl1
   DMd = oo16pi2 * betaMd1
   DMq = oo16pi2 * betaMq1
   DMu = oo16pi2 * betaMu1
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * betaMH11
   DmH(2) = oo16pi2 * betaMH21
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * betaMue1
   DB = oo16pi2 * betaB1
  End If

  !---------------------------------------
  ! to avoid numerical problems in odeint
  !---------------------------------------
  Do i1=1,3
   DMe(i1,i1) = Real(DMe(i1,i1),dp)
   DMl(i1,i1) = Real(DMl(i1,i1),dp)
   DMd(i1,i1) = Real(DMd(i1,i1),dp)
   DMu(i1,i1) = Real(DMu(i1,i1),dp)
   DMq(i1,i1) = Real(DMq(i1,i1),dp)
  End Do
  Dmd = 0.5_dp * ( Dmd + Transpose(Conjg(Dmd)) )
  Dme = 0.5_dp * ( Dme + Transpose(Conjg(Dme)) )
  Dml = 0.5_dp * ( Dml + Transpose(Conjg(Dml)) )
  Dmq = 0.5_dp * ( Dmq + Transpose(Conjg(Dmq)) )
  Dmu = 0.5_dp * ( Dmu + Transpose(Conjg(Dmu)) )

  !--------------------------------------------
  ! This helps avoiding numerical instabilities
  !--------------------------------------------
!   Call Chop(DYe)
!   Call Chop(DYd)
!   Call Chop(DYu)
!   Call Chop(DMhlf)
!   Call Chop(DAe)
!   Call Chop(DAd)
!   Call Chop(DAu)
!   Call Chop(DMe)
!   Call Chop(DMl)
!   Call Chop(DMd)
!   Call Chop(DMu)
!   Call Chop(DMq)
   Call Chop(Dmue)
   Call Chop(DB)


  Call ParametersToG(Dgauge, DYe, DYd, DYu, DMhlf, DAe, DAd, DAu &
                   &, DMe, DMl, DMd, DMq, DMu, DMh, Dmue, DB, f)

  Iname = Iname - 1

 End Subroutine rge213


 Subroutine rge214(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 25.09.01: portation to f90
 !  27.07.13: adding gauge dependence as discussed in 1305.1548
 !            by Dominik Stoeckinger
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2
  Real(dp) :: gauge(3), gauge2(3), sumI, TraceY(3), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)  &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)        &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)            &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)          &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                 &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(5,2)       &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3), hc(4)

  Complex(dp) :: Mhlf(3),DMhlf(3)

  Complex(dp) :: Ae(3,3), Ad(3,3), Au(3,3), aAe(3,3), aAd(3,3), aAu(3,3)   &
     &  , DAe(3,3), DAd(3,3), DAu(3,3), aAdAd(3,3), aAeAe(3,3), aAuAu(3,3) &
     &  , aYdAd(3,3), aYeAe(3,3), aYuAu(3,3), TraceaYA(3), betaAd1(3,3)    &
     &  , betaAd2(3,3), betaAe1(3,3), betaAe2(3,3), betaAu1(3,3), betaAu2(3,3)
  Real(dp) :: TraceA(3)
  Complex(dp) :: aYdYdaYdAd(3,3), aYdAdaYdYd(3,3), TraceAY2(5)               &
     &  , aYeYeaYeAe(3,3), aYeAeaYeYe(3,3), aYuYuaYuAu(3,3), aYuAuaYuYu(3,3) &
     &  , aYuYuaYdAd(3,3), aYuAuaYdYd(3,3), aYdYdaYuAu(3,3), aYdAdaYuYu(3,3)
   
  Complex(dp) :: Me(3,3), Ml(3,3), Md(3,3), Mq(3,3), Mu(3,3), DMe(3,3)        &
     & , DMl(3,3), DMd(3,3), DMq(3,3), DMu(3,3), YdaYd(3,3), YeaYe(3,3)       &
     & , YuaYu(3,3), MdYdaYd(3,3), MeYeaYe(3,3), MuYuaYu(3,3), YdaYdMd(3,3)   &
     & , YeaYeMe(3,3), YuaYuMu(3,3), YdMqaYd(3,3), YeMlaYe(3,3)               &
     & ,  YuMqaYu(3,3), AdaAd(3,3), AeaAe(3,3), AuaAu(3,3), betaMd1(3,3)      &
     & , betaMd2(3,3), betaMe1(3,3), betaMe2(3,3), betaMl1(3,3), betaMl2(3,3) &
     & , betaMq1(3,3), betaMq2(3,3), betaMu1(3,3), betaMu2(3,3), MqaYdYd(3,3) &
     & , MqaYuYu(3,3), aYdYdMq(3,3), aYuYuMq(3,3), aYeYeMl(3,3), MlaYeYe(3,3) &
     & , aYeMeYe(3,3), aYdMdYd(3,3), aYuMuYu(3,3)                             &
     & , YdaYdYdaYd(3,3), YeaYeYeaYe(3,3), YuaYuYuaYu(3,3), MeYeaYeYeaYe(3,3) &
     & , YeaYeYeaYeMe(3,3), YeaYeMeYeaYe(3,3), YeMlaYeYeaYe(3,3), AeaYe(3,3)  &
     & , YeaYeYeMlaYe(3,3), AeaAeYeaYe(3,3), YeaYeAeaAe(3,3), YeaAe(3,3)      &
     & , AeaYeYeaAe(3,3), YeaAeAeaYe(3,3), Tr3aAdYdaAeYe

  Complex(dp) :: Tr3aYdAdaYeAe, AdaYd(3,3), YdaAd(3,3), MlaYeYeaYeYe(3,3)     &
     & , aYeYeaYeYeMl(3,3), aYeYeMlaYeYe(3,3), aYeYeaYeMeYe(3,3)              &
     & , aYeMeYeaYeYe(3,3), aAdYd(3,3),aAeYe(3,3), aAeAeaYeYe(3,3)            &
     & , aYeYeaAeAe(3,3), aAeYeaYeAe(3,3), aYeAeaAeYe(3,3), MdYdaYdYdaYd(3,3) &
     & , YdaYdYdaYdMd(3,3), YdMqaYdYdaYd(3,3),YdaYdMdYdaYd(3,3)               &
     & , YdaYdYdMqaYd(3,3), AdaAdYdaYd(3,3), YdaYDAdaAd(3,3), AdaYdYdaAd(3,3) &
     & , YdaAdAdaYd(3,3)
  Complex(dp) :: MdYdaYuYuaYd(3,3), YdaYuYuaYdMd(3,3), YdMqaYuYuaYd(3,3)      &
     & , YdaYuYuMqaYd(3,3), YdaYuMuYuaYd(3,3), AdaAuYuaYd(3,3)                &
     & , YdaYuAuaAd(3,3), AdaYuYuaAd(3,3), YdaAuAuaYd(3,3), YdaYuYuaYd(3,3)   &
     & , Tr3aYuAu, Tr3aAuYu, YuaAu(3,3)                                       &
     & , MqaYdYdaYdYd(3,3), aYdYdaYdYDMq(3,3), aYdMdYdaYdYd(3,3)              &
     & , aYdYdMqaYdYd(3,3), aYdYdaYdMdYd(3,3), aAdAdaYdYd(3,3)                &
     & , aYdYDaAdAd(3,3), aAdYdaYdAd(3,3), aYdAdaAdYd(3,3), MqaYuYuaYuYu(3,3) &
     & , aYuYuaYuYUMq(3,3), aYuMuYuaYuYu(3,3), aYuYuMqaYuYu(3,3)              &
     & , aYuYuaYuMuYu(3,3), aAuAuaYuYu(3,3), aYuYUaAuAu(3,3), aAuYuaYuAu(3,3) &
     & , aYuAuaAuYu(3,3), aAuYu(3,3), AuaYu(3,3), YuaYdYdaYu(3,3)             &
     & , AuaYdYdaAu(3,3), YuaAdAdaYu(3,3), AuaAdYdaYu(3,3), YuaYdAdaAu(3,3)   &
     & , YuMqaYuYuaYu(3,3), YuaYuYuMqaYu(3,3), MuYuaYuYuaYu(3,3)              &
     & , YuaYuYuaYuMu(3,3), YuaYuMuYuaYu(3,3), AuaAuYuaYu(3,3)                &
     & , YuaYuAuaAu(3,3), AuaYuYuaAu(3,3), YuaAuAuaYu(3,3), MuYuaYdYdaYu(3,3) &
     & , YuaYdYdaYuMu(3,3), YuMqaYdYdaYu(3,3), YuaYdYdMqaYu(3,3)              &
     & , YuaYdMdYdaYu(3,3)

  Real(dp) :: S1, S2, sig(3), Tr3aYdYdaYeYe, Tr3aAdAdaAeAe, AbsGM2(3)         &
     & , Tr3MqaYdYd3aYDMdYd , Tr3MqaYuYu3aYuMu, Tr3aAuAu

  Real(dp) :: Mh(2), DMh(2), TraceMH1(3), TraceMH2(2), betaMH11        &
     &  , betaMH12, betaMH21, betaMH22, q, gamma1, gamma2

  Complex(dp) :: mue, B, Dmue, DB, TraceMue(2), TraceB(2), betaMue1, betaMue2 &
     & , betaB1, betaB2, g2Mi(3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge214'

  OnlyDiagonal = .Not.GenerationMixing

  q = t

  Call GToParameters(gy, gauge, Ye, Yd, Yu                            &
                  & , Mhlf, Ae, Ad, Au, Me, Ml, Md, Mq, Mu, Mh, mue, B)

  gauge2 = gauge**2
  AbsGM2 = gauge2 * Abs( Mhlf )**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYdYd = Matmul(aYd,Yd)
  aYeYe = Matmul(aYe,Ye)
  aYuYu = Matmul(aYu,Yu)

  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aYdYd(i1,i1) = Real(aYdYd(i1,i1), dp)
   aYeYe(i1,i1) = Real(aYeYe(i1,i1), dp)
   aYuYu(i1,i1) = Real(aYuYu(i1,i1), dp)
  End Do
  
  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYdYd),dp )
  TraceY(3) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = 3._dp * TraceY(2) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = 3._dp * TraceY(2) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)

  diagonal(3,1) = 3._dp * TraceY(3)              &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(3,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
!   aYdYdaYuYu = Matmul(aYdYd,aYuYu)
   !------------------------------------------------
   ! this are hermitian matrices, clean up to
   ! avoid numerical problems
   !------------------------------------------------
   Do i1=1,3
    aYdYdaYdYd(i1,i1) = Real(aYdYdaYdYd(i1,i1), dp)
    aYeYeaYeYe(i1,i1) = Real(aYeYeaYeYe(i1,i1), dp)
    aYuYuaYuYu(i1,i1) = Real(aYuYuaYuYu(i1,i1), dp)
   End Do

   Call Adjungate(aYuYuaYdYd, aYdYdaYuYu)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(2)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(2) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(3)
   hd(2) = 9._dp * TraceY(2) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(3)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(3) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(2) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ad,aAd)
  Call Adjungate(Ae,aAe)
  Call Adjungate(Au,aAu)

  aAdAd = Matmul(aAd,Ad)
  aAeAe = Matmul(aAe,Ae)
  aAuAu = Matmul(aAu,Au)

  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aAdAd(i1,i1) = Real(aAdAd(i1,i1), dp)
   aAeAe(i1,i1) = Real(aAeAe(i1,i1), dp)
   aAuAu(i1,i1) = Real(aAuAu(i1,i1), dp)
  End Do

  TraceA(1) = Real( cTrace(aAeAe),dp )
  TraceA(2) = Real( cTrace(aAdAd),dp )
  TraceA(3) = Real( cTrace(aAuAu),dp )

  aYdAd = Matmul(aYd,Ad)
  aYeAe = Matmul(aYe,Ae)
  aYuAu = Matmul(aYu,Au)

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYdAd) 
  TraceaYA(3) = cTrace(aYuAu) 

  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = Matmul(Ae,sume1)

  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(2) + TraceaYA(1)  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    )
  sume1 = 4._dp * aYeAe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + Matmul(Ye,sume1)

  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = Matmul(Ad,sumd1)
  
  diagonal(2,1) = 2._dp * ( 3._dp * TraceaYA(2) + TraceaYA(1) &
                &         - c1_1(2,1) * g2Mi(1) - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do
  betaAd1 = betaAd1 + Matmul(Yd,sumd1)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = Matmul(Au,sumu1)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3)              &
                &         - c1_1(3,1) * g2Mi(1) - c1_1(3,2) * g2Mi(2)   &
                &         - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(3,1)
  End Do
  betaAu1 = betaAu1 + Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdAd = Matmul(aYdYd,aYdAd)
   aYdAdaYdYd = Matmul(aYdAd,aYdYd)
   aYeYeaYeAe = Matmul(aYeYe,aYeAe)
   aYeAeaYeYe = Matmul(aYeAe,aYeYe)
   aYuYuaYuAu = Matmul(aYuYu,aYuAu)
   aYuAuaYuYu = Matmul(aYuAu,aYuYu)
   aYuAuaYdYd = Matmul(aYuAu,aYdYd)
   aYuYuaYdAd = Matmul(aYuYu,aYdAd)
   aYdAdaYuYu = Matmul(aYdAd,aYuYu)
   aYdYdaYuAu = Matmul(aYdYd,aYuAu)
   TraceAY2(1) = cTrace(aYeYeaYeAe)
   TraceAY2(2) = cTrace(aYdYdaYdAd)
   TraceAY2(3) = cTrace(aYuYuaYuAu)
   TraceAY2(4) = cTrace(aYuYuaYdAd)
   TraceAY2(5) = cTrace(aYdYdaYuAu)

  !--------------
  ! A_e
  !--------------
   hd(1) = 6._dp * TraceY(2) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = Matmul(Ae,sume2)
    
   diagonal(1,2) = -6._dp * ( 6._dp * TraceAY2(2) + TraceAY2(4)       &
     &                      + TraceAY2(5) + 2._dp * TraceAY2(1)  )    &
     &  + ( 32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(2)   &
     &  - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(2)                &
     &  + 2.4_dp * gauge2(1) * TraceaYA(1) - 2.4_dp * g2Mi(1) * TraceY(1)    &
     &  - ( 30._dp * g2Mi(2)                              &
     &    + 3.6_dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)      &
     &  - 54._dp * gauge2(1)**2 * Mhlf(1)
   hd(1) = 12._dp * TraceY(2) + 4._dp * TraceY(1)     &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(1) = 18._dp * TraceaYA(2) + 6._dp * TraceaYA(1) &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2)
   sume2 = - 6._dp * aYeYeaYeAe - 8._dp * aYeAeaYeYe  &
         & - hd(1) * aYeAe - hc(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
   betaAe2 = betaAe2 + Matmul(Ye,sume2)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(2) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = Matmul(Ad,sumd2)
    
   diagonal(2,2) = diagonal(1,2)                                   &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                 &
     &              - gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp      &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)    &
     &  + 1.6_dp * gauge2(1) * gauge2(2) * (Mhlf(1)+Mhlf(2))       & 
     &  + 1.856e3_dp * gauge2(1)**2 * Mhlf(1) / 4.5e1_dp

   hd(1) = 12._dp * TraceY(2) + 4._dp * TraceY(1)  &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hd(2) = 6._dp * TraceY(3) - 1.6_dp * gauge2(1) 
   hc(1) = 18._dp * TraceaYA(2) + 6._dp * TraceaYA(1)                 &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2) + 1.6_dp * gauge2(1) *  Mhlf(1)
   hc(2) = 6._dp * TraceaYA(3) + 1.6_dp * gauge2(1) *  Mhlf(1)
   sumd2 = - 6._dp * aYdYdaYdAd - 8._dp * aYdAdaYdYd                  &
       &   - 4._dp * ( aYuAuaYuYu + aYuYuaYuAu + aYuAuaYdYd )         &
       &   - 2._dp * aYuYuaYdAd - hd(1) * aYdAd - hc(1) * aYdYd       &
       &  - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
   betaAd2 = betaAd2 + Matmul(Yd,sumd2)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(3) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = Matmul(Au,sumu2)
    
   diagonal(3,2) =  -6._dp * ( 6._dp * TraceAY2(3) + TraceAY2(4)        &
     &                        + TraceAY2(5)  )                          &
     &  + ( 32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(3)     &
     &  - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(3)           &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                      &
     &              - 3.4_dp * gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp  &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)         &
     &  - ( 30._dp * g2Mi(2)                                &
     &    + 2._dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)         &
     &  - 5486._dp * gauge2(1)**2 * Mhlf(1) / 225._dp
   hd(1) = 6._dp * TraceY(2) + 2._dp * TraceY(1) - 0.8_dp * gauge2(1)
   hc(1) = 6._dp * TraceaYA(2) + 2._dp * TraceaYA(1)   &
       &  + 0.8_dp * gauge2(1) *  Mhlf(1)
   hd(2) = 12._dp * TraceY(3) - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(2) = 18._dp * TraceaYA(3) + 1.2e1_dp * g2Mi(2) + 0.8_dp * g2Mi(1)
   sumu2 = - 6._dp * aYuYuaYuAu - 8._dp * aYuAuaYuYu                     &
       &   - 4._dp * ( aYdAdaYdYd + aYdYdaYdAd + aYdAdaYuYu )            &
       &   - 2._dp * aYdYdaYuAu - hd(1) * aYdAd - hc(1) * aYdYd          &
       &   - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
   betaAu2 = betaAu2 + Matmul(Yu,sumu2)

  End If 
!----------------------------------------------
! beta functions for Sfermion mass parameters
!----------------------------------------------
   S1 = mH(2) - mH(1)
   Do i1=1,3
    S1 = S1 + Real(Me(i1,i1),dp) - Real(Ml(i1,i1),dp) &
       &    + Real(Md(i1,i1),dp) + Real(Mq(i1,i1),dp) &
       &    - 2._dp * Real(Mu(i1,i1),dp)
   End Do
   S1 = S1 * gauge2(1)

   YdaYd = Matmul(Yd,aYd)
   YeaYe = Matmul(Ye,aYe)
   YuaYu = Matmul(Yu,aYu)

   MdYdaYd = Matmul(Md,YdaYd)
   MeYeaYe = Matmul(Me,YeaYe)
   MlaYeYe = Matmul(Ml,aYeYe)
   MqaYdYd = Matmul(Mq,aYdYd)
   MqaYuYu = Matmul(Mq,aYuYu)
   MuYuaYu = Matmul(Mu,YuaYu)

   Call Adjungate(MdYdaYd, YdaYdMd) ! YdaYdMd = Matmul(YdaYd,Md)
   Call Adjungate(MeYeaYe, YeaYeMe) ! YeaYeMe = Matmul(YeaYe,Me)
   Call Adjungate(MlaYeYe, aYeYeMl) ! aYeYeMl = Matmul(aYeYe,Ml)
   Call Adjungate(MqaYdYd, aYdYdMq) ! aYdYdMq = Matmul(aYdYd,Mq)
   Call Adjungate(MqaYuYu, aYuYuMq) ! aYuYuMq = Matmul(aYuYu,Mq)
   Call Adjungate(MuYuaYu, YuaYuMu) ! YuaYuMu = Matmul(YuaYu,Mu)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)

   AdaAd = Matmul(Ad,aAd)
   AeaAe = Matmul(Ae,aAe)
   AuaAu = Matmul(Au,aAu)

   !------------------------------------------------
   ! these are hermitian matrices, clean up to
   ! avoid numerical problems
   !------------------------------------------------
   Do i1=1,3
    aYdMdYd(i1,i1) = Real(aYdMdYd(i1,i1), dp)
    aYeMeYe(i1,i1) = Real(aYeMeYe(i1,i1), dp)
    aYuMuYu(i1,i1) = Real(aYuMuYu(i1,i1), dp)
    YdMqaYd(i1,i1) = Real(YdMqaYd(i1,i1), dp)
    YeMlaYe(i1,i1) = Real(YeMlaYe(i1,i1), dp)
    YuMqaYu(i1,i1) = Real(YuMqaYu(i1,i1), dp)
    AdaAd(i1,i1) = Real(AdaAd(i1,i1), dp)
    AeaAe(i1,i1) = Real(AeaAe(i1,i1), dp)
    AuaAu(i1,i1) = Real(AuaAu(i1,i1), dp)
   End Do

   diagonal(1,1) = - 4.8_dp * AbsGM2(1) + 1.2_dp * S1
   betaMe1 = 2._dp * (MeYeaYe + YeaYeMe)             &
         & + 4._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe )
   Do i1=1,3
    betaMe1(i1,i1) = betaMe1(i1,i1) + diagonal(1,1)
   End Do

   diagonal(2,1) = - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1
   betaMl1 = MlaYeYe + aYeYeMl + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe )
   Do i1=1,3
    betaMl1(i1,i1) = betaMl1(i1,i1) + diagonal(2,1)
   End Do

   diagonal(3,1) = - ( 32._dp * AbsGM2(3) + 1.6_dp * AbsGM2(1) ) / 3._dp &
               & + 0.4_dp * S1
   betaMd1 = 2._dp * (MdYdaYd + YdaYdMd)             &
         & + 4._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd )
   Do i1=1,3
    betaMd1(i1,i1) = betaMd1(i1,i1) + diagonal(3,1)
   End Do

   diagonal(4,1) = - ( 32._dp * AbsGM2(3) + 0.4_dp * AbsGM2(1) ) / 3._dp &
               & - 6._dp * AbsGM2(2) + 0.2_dp * S1
   betaMq1 = MqaYuYu + aYuYuMq + MqaYdYd + aYdYdMq             &
         & + 2._dp * ( mH(2) * aYuYu + mH(1) * aYdYd + aYuMuYu     &
         &           + aYdMdYd + aAuAu + aAdAd )
   Do i1=1,3
    betaMq1(i1,i1) = betaMq1(i1,i1) + diagonal(4,1)
   End Do

   diagonal(5,1) = - ( 32._dp * AbsGM2(3) + 6.4_dp * AbsGM2(1) ) / 3._dp &
               &   - 0.8_dp * S1
   betaMu1 = 2._dp * (MuYuaYu + YuaYuMu)             &
         & + 4._dp * ( mH(2) * YuaYu + YuMqaYu + AuaAu )
   Do i1=1,3
    betaMu1(i1,i1) = betaMu1(i1,i1) + diagonal(5,1)
   End Do

   If (TwoLoopRGE) Then
    YdaYdYdaYd = MatSquare(YdaYd,OnlyDiagonal)
    YeaYeYeaYe = MatSquare(YeaYe,OnlyDiagonal)
    YuaYuYuaYu = MatSquare(YuaYu,OnlyDiagonal)

    !------------------------------------------------
    ! these are hermitian matrices, clean up to
    ! avoid numerical problems
    !------------------------------------------------
    Do i1=1,3
     YdaYdYdaYd(i1,i1) = Real( YdaYdYdaYd(i1,i1), dp)
     YeaYeYeaYe(i1,i1) = Real( YeaYeYeaYe(i1,i1), dp)
     YuaYuYuaYu(i1,i1) = Real( YuaYuYuaYu(i1,i1), dp)
    End Do

    AdaYd = Matmul(Ad,aYd)
    AeaYe = Matmul(Ae,aYe)
    AuaYu = Matmul(Au,aYu)

    aAdYd = Matmul(aAd,Yd)
    aAeYe = Matmul(aAe,Ye)
    aAuYu = Matmul(aAu,Yu)

    Call Adjungate(AdaYd,YdaAd) ! YdaAd = Matmul(Yd,aAd,OnlyDiagonal)
    Call Adjungate(AeaYe,YeaAe) ! YeaAe = Matmul(Ye,aAe)
    Call Adjungate(AuaYu,YuaAu) ! YuaAu = Matmul(Yu,aAu)

    YdaYuYuaYd = MatMul3(Yd,aYuYu,aYd,OnlyDiagonal)
    AdaYuYuaAd = MatMul3(Ad,aYuYu,aAd,OnlyDiagonal)
    YdaAuAuaYd = MatMul3(Yd,aAuAu,aYd,OnlyDiagonal)
    AdaAuYuaYd = MatMul4(Ad,aAu,Yu,aYd,OnlyDiagonal)
    YdaYuAuaAd = MatMul3(Yd,aYuAu,aAd,OnlyDiagonal)

    YuaYdYdaYu = MatMul3(Yu,aYdYd,aYu,OnlyDiagonal)
    AuaYdYdaAu = MatMul3(Au,aYdYd,aAu,OnlyDiagonal)
    YuaAdAdaYu = MatMul3(Yu,aAdAd,aYu,OnlyDiagonal)
    AuaAdYdaYu = MatMul4(Au,aAd,Yd,aYu,OnlyDiagonal)
    YuaYdAdaAu = MatMul3(Yu,aYdAd,aAu,OnlyDiagonal)

    MdYdaYuYuaYd = Matmul(Md,YdaYuYuaYd)
    Call Adjungate(MdYdaYuYuaYd, YdaYuYuaYdMd)
    YdMqaYuYuaYd = MatMul3(Yd,MqaYuYu,aYd,OnlyDiagonal)
    Call Adjungate(YdMqaYuYuaYd, YdaYuYuMqaYd)
    YdaYuMuYuaYd = MatMul3(Yd,aYuMuYu,aYd,OnlyDiagonal)

    MuYuaYdYdaYu = Matmul(Mu,YuaYdYdaYu)
    Call Adjungate(MuYuaYdYdaYu, YuaYdYdaYuMu)
    YuMqaYdYdaYu = MatMul3(Yu,MqaYdYd,aYu,OnlyDiagonal)
    Call Adjungate(YuMqaYdYdaYu, YuaYdYdMqaYu)
    YuaYdMdYdaYu = MatMul3(Yu,aYdMdYd,aYu,OnlyDiagonal)

    MeYeaYeYeaYe = Matmul(MeYeaYe,YeaYe)
    Call Adjungate(MeYeaYeYeaYe,YeaYeYeaYeMe)
    aYeMeYeaYeYe = Matmul(aYeMeYe,aYeYe)
    Call Adjungate(aYeMeYeaYeYe,aYeYeaYeMeYe)
    YeaYeMeYeaYe = Matmul(YeaYeMe,YeaYe)

    MlaYeYeaYeYe = Matmul(MlaYeYe,aYeYe)
    Call Adjungate(MlaYeYeaYeYe, aYeYeaYeYeMl)
    YeMlaYeYeaYe = Matmul(YeMlaYe,YeaYe)
    Call Adjungate(YeMlaYeYeaYe, YeaYeYeMlaYe)
    aYeYeMlaYeYe = Matmul(aYeYeMl,aYeYe)

    MdYdaYdYdaYd = Matmul(MdYdaYd,YdaYd)
    Call Adjungate(MdYdaYdYdaYd, YdaYdYdaYdMd)
    aYdMdYdaYdYd = Matmul(aYdMdYd,aYdYd)
    Call Adjungate(aYdMdYdaYdYd, aYdYdaYdMdYd)
    YdaYdMdYdaYd = Matmul(YdaYdMd,YdaYd)

    MqaYdYdaYdYd = Matmul(MqaYdYd,aYdYd)
    Call Adjungate(MqaYdYdaYdYd, aYdYdaYdYdMq)
    YdMqaYdYdaYd = Matmul(YdMqaYd,YdaYd)
    Call Adjungate(YdMqaYdYdaYd, YdaYdYdMqaYd)
    aYdYdMqaYdYd = Matmul(aYdYdMq,aYdYd)

    MqaYuYuaYuYu = Matmul(MqaYuYu,aYuYu)
    Call Adjungate(MqaYuYuaYuYu, aYuYuaYuYuMq)
    YuMqaYuYuaYu = Matmul(YuMqaYu,YuaYu)
    Call Adjungate(YuMqaYuYuaYu, YuaYuYuMqaYu)
    aYuYuMqaYuYu = Matmul(aYuYuMq,aYuYu)

    MuYuaYuYuaYu = Matmul(MuYuaYu,YuaYu)
    Call Adjungate(MuYuaYuYuaYu, YuaYuYuaYuMu)
    aYuMuYuaYuYu = Matmul(aYuMuYu,aYuYu)
    Call Adjungate(aYuMuYuaYuYu, aYuYuaYuMuYu)
    YuaYuMuYuaYu = Matmul(YuaYuMu,YuaYu)

    AdaAdYdaYd = Matmul(AdaAd,YdaYd)
    Call Adjungate(AdaAdYdaYd, YdaYdAdaAd)
    AdaYdYdaAd = Matmul(AdaYd,YdaAd)
    YdaAdAdaYd = Matmul(YdaAd,AdaYd)

    aAdAdaYdYd = Matmul(aAdAd,aYdYd)
    Call Adjungate(aAdAdaYdYd, aYdYdaAdAd)
    aAdYdaYdAd = Matmul(aAdYd,aYdAd)
    aYdAdaAdYd = Matmul(aYdAd,aAdYd)

    AeaAeYeaYe = Matmul(AeaAe,YeaYe)
    Call Adjungate(AeaAeYeaYe, YeaYeAeaAe)
    AeaYeYeaAe = Matmul(AeaYe,YeaAe)
    YeaAeAeaYe = Matmul(YeaAe,AeaYe)

    aAeAeaYeYe = Matmul(aAeAe,aYeYe)
    Call Adjungate(aAeAeaYeYe, aYeYeaAeAe)
    aAeYeaYeAe = Matmul(aAeYe,aYeAe)
    aYeAeaAeYe = Matmul(aYeAe,aAeYe)

    AuaAuYuaYu = Matmul(AuaAu,YuaYu)
    Call Adjungate(AuaAuYuaYu, YuaYuAuaAu)
    AuaYuYuaAu = Matmul(AuaYu,YuaAu)
    YuaAuAuaYu = Matmul(YuaAu,AuaYu)

    aAuAuaYuYu = Matmul(aAuAu,aYuYu)
    Call Adjungate(aAuAuaYuYu, aYuYuaAuAu)
    aAuYuaYuAu = Matmul(aAuYu,aYuAu)
    aYuAuaAuYu = Matmul(aYuAu,aAuYu)

    S2 = (1.5_dp * gauge2(2) + 0.3_dp * gauge2(1) )            &
     &      * (MH(2) - MH(1) - Real(cTrace(ML),dp) )              &
     & + ( (8._dp * gauge2(3) + 0.1_dp*gauge2(1)) / 3._dp      &
     &   + 1.5_dp * gauge2(2) ) * Real( cTrace(Mq),dp   )         &
     & - (16._dp * gauge2(3) + 3.2_dp*gauge2(1) )              &
     &    * Real(cTrace(Mu),dp) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md),dp) / 3._dp                           &
     & + 1.2_dp*gauge2(1) * Real(cTrace(Me),dp)                   &
     & - 3._dp * (MH(2)*TraceY(3) - MH(1) * TraceY(2) )        &
     & + mH(1) * TraceY(1)
    
    Do i1=1,3
     S2 = S2 - Real(YuMqaYu(i1,i1),dp) + 4._dp * Real(aYuMuYu(i1,i1),dp)   &
        &    - Real(YdMqaYd(i1,i1),dp) - 2._dp * Real(aYdMdYd(i1,i1),dp)   &
        &    + Real(YeMlaYe(i1,i1),dp) - 2._dp * Real(aYeMeYe(i1,i1),dp)
    End Do

    sig(1) = 3._dp * (MH(1) + MH(2) + Real(cTrace(Ml),dp) )          &
         & + Real(cTrace(Mq),dp) + 8._dp * Real(cTrace(Mu),dp)          &
         & + 2._dp * Real(cTrace(Md),dp) + 6._dp * Real(cTrace(Me),dp) 
    sig(1) = 0.2_dp * gauge2(1) * sig(1)
    sig(2) = gauge2(2) * ( MH(1) + MH(2) + Real( cTrace(Ml),dp )    &
           &             + 3._dp * Real( cTrace(Mq),dp ) )
    sig(3) = gauge2(3) * ( 2._dp * Real( cTrace(Mq),dp )            &
           &             + Real(cTrace(Mu),dp) + Real(cTrace(Md),dp) )

    Tr3aYdYdaYeYe = 3._dp * TraceY(2) + TraceY(1)
    Tr3MqaYdYd3aYDMdYd = Real(cTrace(MlaYeYe),dp) + Real(cTrace(aYeMeYe),dp)  &
              & + 3._dp * ( Real(cTrace(MqaYdYd),dp)+Real(cTrace(aYdMdYd),dp) )
    Tr3aAdAdaAeAe = 3._dp * TraceA(2) + TraceA(1)
    Tr3aYdAdaYeAe = 3._dp * TraceaYA(2) + TraceaYA(1)
    Tr3aAdYdaAeYe = Conjg( Tr3aYdAdaYeAe )
    Tr3MqaYuYu3aYuMu = 3._dp * (Real(cTrace(MqaYuYu),dp) &
                     &         +Real(cTrace(aYuMuYu),dp) )
    Tr3aAuAu = 3._dp * TraceA(3)
    Tr3aYuAu = 3._dp * TraceaYA(3)
    Tr3aAuYu = Conjg( Tr3aYuAu )

    diagonal(1,2) = 2.4_dp * gauge2(1) * (S2 + sig(1) )   &
                & + 112.32_dp * gauge2(1) * AbsGM2(1)
    hd(1) = 6._dp*gauge2(2)-1.2_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 4.8_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 2.4_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMe2 = -2._dp * ( MeYeaYeYeaYe + YeaYeYeaYeMe )                       &
      & - 4._dp * ( YeMlaYeYeaYe + YeaYeMeYeaYe + YeaYeYeMlaYe )             &
      & - 8._dp * mH(1) * YeaYeYeaYe                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MeYeaYe + YeaYeMe + 4._dp * MH(1) * YeaYe &
      &                           + 2._dp * (YeMlaYe + AeaAE)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YeaYe + AeaAeYeaYe + YeaYEAeaAe     &
      &           + AeaYeYeaAe + YeaAeAeaYe + Tr3aAdAdaAeAe * YeaYe          &
      &           + Tr3aAdYdaAeYe * AeaYe + Tr3aYdAdaYeAe * YeaAe )          &
      & + hd(1) * ( MeYeaYe + YeaYeMe                                        &
      &           + 2._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe ) )            &
      & + hd(2) * YeaYe + hc(1) * AeaYe + hc(2) * YeaAe
    Do i1=1,3
     betaMe2(i1,i1) = betaMe2(i1,i1) + diagonal(1,2)
    End Do

    diagonal(2,2) = gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)        &
     &     + 3._dp * gauge2(2) * sig(2)                                 &
     &     + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)      &
     &                   + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2       &
     &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )     &
     &     + 24.84_dp * gauge2(1) * AbsGM2(1)
    betaMl2 = -2._dp * ( MlaYeYeaYeYe + aYeYeaYeYeMl )                   &
      & - 4._dp * ( aYeMeYeaYeYe + aYeYeMlaYeYe + aYeYeaYeMeYe )         &
      & - 8._dp * mH(1) * aYeYeaYeYe                                     &
      & - Tr3aYdYdaYeYe * ( MlaYeYe + aYeYeML + 4._dp * MH(1) * aYeYe    &
      &                   + 2._dp * (aYeMeYe + aAeAE)  )                 &
      & - 4._dp * (aAeAeaYeYe + aYeYeaAeAe + aAeYeaYeAe + aYeAeaAeYe )   &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYeYe +  Tr3aAdAdaAeAe * aYeYe  &
      &           + Tr3aAdYdaAeYe * aYeAe + Tr3aYdAdaYeAe * aAeYe     )  &
      & + 1.2_dp*gauge2(1) * ( MlaYeYe + aYeYeMl                         &
      &                      + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe &
      &                                - Mhlf(1) * aAeYe                 &
      &                                - Conjg(Mhlf(1)) * aYeAe ) )      &
      & + 4.8_dp * AbsGM2(1) * aYeYe
    Do i1=1,3
     betaMl2(i1,i1) = betaMl2(i1,i1) + diagonal(2,2)
    End Do

    diagonal(3,2) = 0.8_dp * gauge2(1) * (S2 + sig(1)/3._dp )             &
      &   + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp   &
      &   + 8.08e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                    &
      &   + 1.28e2_dp * ( gauge2(1) * AbsGM2(3)                           &
      &                 + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)     &
      &                   * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)+0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) + 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) - 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMd2 = -2._dp * ( MdYdaYdYdaYd + YdaYdYdaYdMd )                       &
      & - 4._dp * ( YdMqaYdYdaYd + YdaYdMdYdaYd + YdaYdYdMqaYd )             &
      & - 8._dp * mH(1) * YdaYdYdaYd                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MdYdaYd + YdaYdMd + 4._dp * MH(1) * YdaYd &
      &                           + 2._dp * (YdMqaYd + AdaAD)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YdaYd + AdaAdYdaYd + YdaYDAdaAd     &
      &           + AdaYdYdaAd + YdaAdAdaYd + Tr3aAdAdaAeAe * YdaYd          &
      &           + Tr3aAdYdaAeYe * AdaYd + Tr3aYdAdaYeAe * YdaAd )          &
      & + hd(1) * ( MdYdaYd + YdaYdMd                                        &
      &           + 2._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd ) )            &
      & + hd(2) * YdaYd + hc(1) * AdaYd + hc(2) * YdaAd                      &
      & - 2._dp * (MdYdaYuYuaYd + YdaYuYuaYdMd )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YdaYuYuaYd + YdMqaYuYuaYd+ YdaYuYuMqaYd  &
      &           + YdaYuMuYuaYd + AdaAuYuaYd + YdaYuAuaAd                   &
      &           + AdaYuYuaAd + YdaAuAuaYd )
    Do i1=1,3
     betaMd2(i1,i1) = betaMd2(i1,i1) + diagonal(3,2)
    End Do

    diagonal(4,2) = 0.2_dp * gauge2(1) * (2._dp * S2 + sig(1)/3._dp )     &
      &  + 3._dp * gauge2(2) * sig(2)                                     &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp    &
      &  + 1.99e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                     &
      &  + 32._dp * ( gauge2(1) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)         &
      &              * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp      &
      &  + 33._dp * gauge2(2) * AbsGM2(2)                                 &
      &  + 32._dp * ( gauge2(2) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(2)  + gauge2(3)*gauge2(2)        &
      &              * Real( Mhlf(2) * Conjg(Mhlf(3)),dp ) )                 &
      &  + 0.4_dp * ( gauge2(2) * AbsGM2(1)                               &
      &             + gauge2(1) * AbsGM2(2) + gauge2(1)*gauge2(2)         &
      &              * Real( Mhlf(2) * Conjg(Mhlf(1)),dp ) )
    hd(1) = 1.6_dp * AbsGM2(1)
    hc(1) = - 0.8_dp * g2Mi(1)
    hc(2) = Conjg( hc(1) )
    hd(2) = 3.2_dp * AbsGM2(1)
    hc(3) = - 1.6_dp * g2Mi(1)
    hc(4) = Conjg( hc(3) )
    betaMq2 = -2._dp * ( MqaYdYdaYdYd + aYdYdaYdYDMq )                       &
      & - 4._dp * ( aYdMdYdaYdYd + aYdYdMqaYdYd + aYdYdaYdMdYd )             &
      & - 8._dp * mH(1) * aYdYdaYdYd                                         &
      & - Tr3aYdYdaYeYe * ( MqaYdYd + aYdYdMq + 4._dp * MH(1) * aYdYd        &
      &                   + 2._dp * (aYdMdYd + aAdAd)  )                     &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYdYd                               &
      &           + 2._dp * ( aAdAdaYdYd + aYdYDaAdAd + aAdYdaYdAd           &
      &                     + aYdAdaAdYd )                                   &
      &           + Tr3aAdAdaAeAe * aYdYd + Tr3aAdYdaAeYe * aYdAd            &
      &           + Tr3aYdAdaYeAe * aAdYd )                                  &
      & + 0.4_dp*gauge2(1) * ( MqaYdYd + aYdYdMq                             &
      &                      + 2._dp * ( mH(1) * aYdYd + aYdMdYd + aAdAd ) ) &
      & + hd(1) * aYdYd + hc(1) * aAdYd + hc(2) * aYdAd                      &
      & - 2._dp * ( MqaYuYuaYuYu + aYuYuaYuYuMq )                            &
      & - 4._dp * ( aYuMuYuaYuYu + aYuYuMqaYuYu + aYuYuaYuMuYu )             &
      & - 8._dp * MH(2) * aYuYuaYuYu                                         &
      & - 3._dp * TraceY(3) * ( MqaYuYu + aYuYuMq + 4._dp * MH(2) * aYuYu    &
      &                       + 2._dp * (aYuMuYu + aAuAu) )                  &
      & - 2._dp * ( Tr3MqaYuYu3aYuMu * aYuYu                                 &
      &           + 2._dp * ( aAuAuaYuYu + aYuYuaAuAu                        &
      &                     + aAuYuaYuAu + aYuAuaAuYu )                      &
      &           + Tr3aAuAu * aYuYu + Tr3aAuYu * aYuAu + Tr3aYuAu * aAuYu ) &
      & + 0.8_dp*gauge2(1) * ( MqaYuYu + aYuYuMq                             &
      &                      + 2._dp * ( MH(2) * aYuYu + aYuMuYu + aAuAu ) ) &
      & + hd(2) * aYuYu + hc(3) * aAuYu + hc(4) * aYuAu
    Do i1=1,3
     betaMq2(i1,i1) = betaMq2(i1,i1) + diagonal(4,2)
    End Do

    diagonal(5,2) = 1.6_dp * gauge2(1) * (2._dp*sig(1)/3._dp - S2)           &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp       &
      &  + 3424._dp * gauge2(1) * AbsGM2(1) / 75._dp                         &
      &  + 512._dp * ( gauge2(1) * AbsGM2(3)                                 &
      &              + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)           &
      &                * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)-0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMu2 = -2._dp * ( MuYuaYuYuaYu + YuaYuYuaYuMu )                       &
      & - 4._dp * ( YuMqaYuYuaYu + YuaYuMuYuaYu + YuaYuYuMqaYu )             &
      & - 8._dp * MH(2) * YuaYuYuaYu                                         &
      & - 6._dp * TraceY(3) * ( MuYuaYu + YuaYuMu + 4._dp * MH(2) * YuaYu    &
      &                       + 2._dp * (YuMqaYu + AuaAu)  )                 &
      & - 4._dp * ( Tr3MqaYuYu3aYuMu * YuaYu + AuaAuYuaYu + YuaYuAuaAu       &
      &           + AuaYuYuaAu + YuaAuAuaYu + Tr3aAuAu * YuaYu               &
      &           + Tr3aAuYu * AuaYu + Tr3aYuAu * YuaAu )                    &
      & + hd(1) * ( MuYuaYu + YuaYuMu                                        &
      &           + 2._dp * ( MH(2) * YuaYu + YuMqaYu + AuaAu ) )            &
      & + hd(2) * YuaYu + hc(1) * AuaYu + hc(2) * YuaAu                      &
      & - 2._dp * (MuYuaYdYdaYu + YuaYdYdaYuMu )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YuaYdYdaYu + YuMqaYdYdaYu + YuaYdYdMqaYu &
      &           + YuaYdMdYdaYu + AuaAdYdaYu + YuaYdAdaAu                   &
      &           + AuaYdYdaAu + YuaAdAdaYu )
    Do i1=1,3
     betaMu2(i1,i1) = betaMu2(i1,i1) + diagonal(5,2)
    End Do

   End If 

  !------------------------------------------
  ! beta functions for Higgs mass parameters
  !------------------------------------------
   traceMH1(1) = mH(1) * TraceY(1) + Real( cTrace(YeMlaYe),dp ) &
             & + Real( cTrace(aYeMeYe),dp ) + TraceA(1)
   traceMH1(2) = mH(1) * TraceY(2) + Real( cTrace(YdMqaYd),dp ) &
             & + Real( cTrace(aYdMdYd),dp ) + TraceA(2)
   betamH11 = 6._dp * TraceMH1(2) + 2._dp * TraceMH1(1)      &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1

   traceMH2(1) = mH(2) * TraceY(3) + Real( cTrace(YuMqaYu),dp )  &
             & + Real( cTrace(aYuMuYu),dp ) + TraceA(3)

   betamH21 = 6._dp * TraceMH2(1) - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1)  &
          & + 0.6_dp * S1

   If (TwoLoopRGE) Then
    traceMH1(3) = MH(1) * (6._dp*TraceY2(2) + 2._dp*TraceY2(1) + TraceY2(4) ) &
              & + MH(2) * TraceY2(4)                                          &
              & + 6._dp * ( Real( cTrace(MqaYdYdaYdYd),dp )                   &
              &           + Real( cTrace(aYdMdYdaYdYd),dp )                   &
              &           + Real( cTrace(aAdAdaYdYd),dp )                     &
              &           + Real( cTrace(aAdYdaYdAd),dp )  )                  &
              & + 2._dp * ( Real( cTrace(MlaYeYeaYeYe),dp )                   &
              &           + Real( cTrace(aYeMeYeaYeYe),dp )                   &
              &           + Real( cTrace(aAeAeaYeYe),dp )                     &
              &           + Real( cTrace(aAeYeaYeAe),dp )  )                  &
           & + Real(cTrace(YdMqaYuYuaYd),dp) + Real(cTrace(YdaYuMuYuaYd),dp ) &
       & + Real( cTrace(YdaYuYuMqaYd),dp ) + Real( cTrace(YuaYdMdYdaYu),dp ) &
       & + Real( cTrace(YdaAuAuaYd),dp ) + Real( cTrace(AdaYuYuaAd),dp )     &
              & + Real( cTrace(AdaAuYuaYd),dp ) + Real( cTrace(YdaYuAuaAd),dp )
    betaMH12 = - 6._dp * traceMH1(3)                                       &
      &   + (32._dp*gauge2(3) - 0.8_dp*gauge2(1) ) * traceMH1(2)           &
      &   + 64._dp * ( AbsGM2(3) * TraceY(2)                               &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(2),dp ) )     &
      &   - 1.6_dp * ( AbsGM2(1) * TraceY(2)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(2),dp ) )     &
      &   + 2.4_dp*gauge2(1) * traceMH1(1)                                 &
      &   + 4.8_dp * ( AbsGM2(1) * TraceY(1)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(1),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)                   &
      &   + 3._dp * gauge2(2) * sig(2)                                     &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)          &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2           &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

    traceMH2(2) = MH(2) * ( 6._dp * TraceY2(3) + TraceY2(4) )                 &
      &  + MH(1) * TraceY2(4)                                                 &
      & + 6._dp*(Real(cTrace(MqaYuYuaYuYu),dp)+Real(cTrace(aYuMuYuaYuYu),dp)  &
      &        + Real(cTrace(aAuAuaYuYu),dp) + Real(cTrace(aAuYuaYuAu),dp))   &
      & + Real(cTrace(YuMqaYdYdaYu),dp)+Real(cTrace(YuaYdMdYdaYu),dp )        &
      & + Real(cTrace(YuaYdYdMqaYu),dp)+Real(cTrace(YdaYuMuYuaYd),dp )        &
      & + Real(cTrace(YuaAdAdaYu),dp) + Real( cTrace(AuaYdYdaAu),dp )         &
      &  + Real( cTrace(AuaAdYdaYu),dp ) + Real( cTrace(YuaYdAdaAu),dp )
    betaMH22 = - 6._dp * traceMH2(2)                                          &
      &   + (32._dp*gauge2(3) + 1.6_dp*gauge2(1) ) * traceMH2(1)              &
      &   + 64._dp * ( AbsGM2(3) * TraceY(3)                                  &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(3),dp ) )     &
      &   + 3.2_dp * ( AbsGM2(1) * TraceY(3)                                  &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(3),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) + 1.2_dp * S2)                      &
      &   + 3._dp * gauge2(2) * sig(2)                                        &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)             &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2              &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

   End If
!-----------------------------
! beta functions for mu and B
!-----------------------------
   TraceMue(1) = 3._dp * (TraceY(2)+TraceY(3)) + TraceY(1)  &
             & - 3._dp * gauge2(2) - 0.6_dp * gauge2(1)
   betaMue1 = mue * TraceMue(1)

   TraceB(1) = 6._dp * (TraceaYA(2)+TraceaYA(3)) + 2._dp * TraceaYA(1)    &
           & + 6._dp * g2Mi(2) + 1.2_dp * g2Mi(1)
   betaB1 = mue * TraceB(1) + B * TraceMue(1)

   If (TwoLoopRGE) Then
    TraceMue(2) = - 3._dp * ( 3._dp * (TraceY2(2) + TraceY2(3) )          &
      &                     + 2._dp * TraceY2(4) + TraceY2(1) )           &
      &         + (16._dp * gauge2(3) + 0.8_dp * gauge2(1)) * TraceY(3)   &
      &         + (16._dp * gauge2(3) - 0.4_dp * gauge2(1)) * TraceY(2)   &
      &         + 1.2_dp * gauge2(1) * TraceY(1)                          &
      &         + 7.5_dp * gauge2(2)**2                                   &
      &         + 1.8_dp * gauge2(2) * gauge2(1)                          &
      &         + 4.14_dp * gauge2(1)**2
    betaMue2 = mue * TraceMue(2)

    TraceB(1) = cTrace( 3._dp * (Matmul(AuaYu,YuaYu) + Matmul(AdaYd,YdaYd) )  &
              &       + Matmul(AeaYe,YeaYe) + Matmul(aYuAu,aYdYd)             &
              &       + Matmul(aYdAd,aYuYu) ) 
    TraceB(2) = -12._dp * TraceB(1)                                           &
      &   + (32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(3)          &
      &   + (32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(2)          &
      &   + 2.4_dp * gauge2(1) * TraceaYA(1)                                  &
      &   - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(3)               &
      &   - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(2)               &
      &   - 2.4_dp * g2Mi(1) * TraceY(1)                          &
      &   - 30._dp * gauge2(2)**2 * Mhlf(2)                                   &
      &   - 3.6_dp * gauge2(2) * gauge2(1) * (Mhlf(1) + Mhlf(2) )             &
      &   - 16.56_dp * gauge2(1)**2 * Mhlf(1) 
    betaB2 = mue * TraceB(2) + B * TraceMue(2)

   End If

  !---------------
  ! Ln(tan(beta))
  !---------------
  gamma1 = 3._dp * (TraceY(2) - TraceY(3)) + TraceY(1)

  If (TwoLoopRGE) Then 
   gamma2 = 0.75_dp * ( 3._dp * (TraceY(3)**2 - TraceY(2)**2) - TraceY(1)**2)  &
   &   - (1.9_dp * gauge2(1) + 4.5_dp * gauge2(2) + 20._dp * gauge2(3) )       &
   &                                                              * TraceY(3)  &
   &   + (0.4_dp * gauge2(1) + 4.5_dp * gauge2(2) + 20._dp * gauge2(3) )       &
   &                                                              * TraceY(2)  &
   &   + (1.8_dp * gauge2(1) + 1.5_dp * gauge2(2) ) * TraceY(1)                &
   &   + (0.3_dp * gauge2(1) + 1.5_dp * gauge2(2) ) * gamma1 ! gauge dependence
  end If

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1 + oo16pi2 * (Matmul(b_2,gauge2) - Matmul(a_2,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   Do i1 = 1,3    
    sumI = 0._dp
    Do i2=1,3
     sumI = sumI + b_2(i1,i2) * gauge2(i2) * (Mhlf(i1) + Mhlf(i2) )     &
          &      + a_2(i1,i2) * ( TraceaYA(i2) - Mhlf(i1)*TraceY(i2) )
    End Do
    DMhlf(i1) = oo8pi2 * gauge2(i1) * ( b_1(i1) * Mhlf(i1) + oo16pi2 * sumI)
   End Do
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * ( betaAe1 + oo16pi2 * betaAe2 )
   DAd = oo16pi2 * ( betaAd1 + oo16pi2 * betaAd2 )
   DAu = oo16pi2 * ( betaAu1 + oo16pi2 * betaAu2 )
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * ( betaMe1 + oo16pi2 * betaMe2 )
   DMl = oo16pi2 * ( betaMl1 + oo16pi2 * betaMl2 )
   DMd = oo16pi2 * ( betaMd1 + oo16pi2 * betaMd2 )
   DMq = oo16pi2 * ( betaMq1 + oo16pi2 * betaMq2 )
   DMu = oo16pi2 * ( betaMu1 + oo16pi2 * betaMu2 )
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * ( betaMH11 + oo16pi2 * betaMH12 )
   DmH(2) = oo16pi2 * ( betaMH21 + oo16pi2 * betaMH22 )
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * ( betaMue1 + oo16pi2 * betaMue2 )
   DB = oo16pi2 * ( betaB1 + oo16pi2 * betaB2 )

 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   DMhlf = oo8pi2 * gauge2 * b_1 * Mhlf
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * betaAe1
   DAd = oo16pi2 * betaAd1
   DAu = oo16pi2 * betaAu1
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * betaMe1
   DMl = oo16pi2 * betaMl1
   DMd = oo16pi2 * betaMd1
   DMq = oo16pi2 * betaMq1
   DMu = oo16pi2 * betaMu1
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * betaMH11
   DmH(2) = oo16pi2 * betaMH21
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * betaMue1
   DB = oo16pi2 * betaB1
  End If

  !---------------------------------------
  ! to avoid numerical problems in odeint
  !---------------------------------------
  Do i1=1,3
   DMe(i1,i1) = Real(DMe(i1,i1),dp)
   DMl(i1,i1) = Real(DMl(i1,i1),dp)
   DMd(i1,i1) = Real(DMd(i1,i1),dp)
   DMu(i1,i1) = Real(DMu(i1,i1),dp)
   DMq(i1,i1) = Real(DMq(i1,i1),dp)
  End Do
  Dmd = 0.5_dp * ( Dmd + Transpose(Conjg(Dmd)) )
  Dme = 0.5_dp * ( Dme + Transpose(Conjg(Dme)) )
  Dml = 0.5_dp * ( Dml + Transpose(Conjg(Dml)) )
  Dmq = 0.5_dp * ( Dmq + Transpose(Conjg(Dmq)) )
  Dmu = 0.5_dp * ( Dmu + Transpose(Conjg(Dmu)) )

  !--------------------------------------------
  ! This helps avoiding numerical instabilities
  !--------------------------------------------
!   Call Chop(DYe)
!   Call Chop(DYd)
!   Call Chop(DYu)
!   Call Chop(DMhlf)
!   Call Chop(DAe)
!   Call Chop(DAd)
!   Call Chop(DAu)
!   Call Chop(DMe)
!   Call Chop(DMl)
!   Call Chop(DMd)
!   Call Chop(DMu)
!   Call Chop(DMq)
   Call Chop(Dmue)
   Call Chop(DB)


  Call ParametersToG(Dgauge, DYe, DYd, DYu, DMhlf, DAe, DAd, DAu &
                   &, DMe, DMl, DMd, DMq, DMu, DMh, Dmue, DB, f(1:213))
  If (TwoLoopRGE) Then 
   f(214) = oo16pi2 * (gamma1 + oo16pi2 * gamma2 )
  Else
   f(214) = oo16pi2 * gamma1
  End If


  Iname = Iname - 1

 End Subroutine rge214


 Subroutine rge267(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.03.2001: implementing right-handed neutrinos  at 1-loop
 !             up to now: Y_nu, A_nu,
 !             the parameters m_H1, M_H2, B and mu still need to be changed
 ! 07.10.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2
  Real(dp) :: gauge(3), gauge2(3), TraceY(4), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2)           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3), hc(4),Ynu(3,3) &
    & , aYnu(3,3), aYnuYnu(3,3), sumnu1(3,3), betaYnu1(3,3) ,DYnu(3,3), sumI

  Complex(dp) :: Mhlf(3),DMhlf(3)

  Complex(dp) :: Ae(3,3), Ad(3,3), Au(3,3), aAe(3,3), aAd(3,3), aAu(3,3)   &
     &  , DAe(3,3), DAd(3,3), DAu(3,3), aAdAd(3,3), aAeAe(3,3), aAuAu(3,3) &
     &  , aYdAd(3,3), aYeAe(3,3), aYuAu(3,3), TraceaYA(4), betaAd1(3,3)    &
     &  , betaAd2(3,3), betaAe1(3,3), betaAe2(3,3), betaAu1(3,3)           &
     &  , betaAu2(3,3), Anu(3,3), aAnu(3,3), DAnu(3,3) ,aAnuAnu(3,3)       &
     &  , betaAnu1(3,3), aYnuAnu(3,3)!, betaAnu2(3,3)
  Real(dp) :: TraceA(4)
  Complex(dp) :: aYdYdaYdAd(3,3), aYdAdaYdYd(3,3), TraceAY2(5)               &
     &  , aYeYeaYeAe(3,3), aYeAeaYeYe(3,3), aYuYuaYuAu(3,3), aYuAuaYuYu(3,3) &
     &  , aYuYuaYdAd(3,3), aYuAuaYdYd(3,3), aYdYdaYuAu(3,3), aYdAdaYuYu(3,3)
   
  Complex(dp) :: Me(3,3), Ml(3,3), Md(3,3), Mq(3,3), Mu(3,3), DMe(3,3)        &
     & , DMl(3,3), DMd(3,3), DMq(3,3), DMu(3,3), YdaYd(3,3), YeaYe(3,3)       &
     & , YuaYu(3,3), MdYdaYd(3,3), MeYeaYe(3,3), MuYuaYu(3,3), YdaYdMd(3,3)   &
     & , YeaYeMe(3,3), YuaYuMu(3,3), YdMqaYd(3,3), YeMlaYe(3,3)               &
     & ,  YuMqaYu(3,3), AdaAd(3,3), AeaAe(3,3), AuaAu(3,3), betaMd1(3,3)      &
     & , betaMd2(3,3), betaMe1(3,3), betaMe2(3,3), betaMl1(3,3), betaMl2(3,3) &
     & , betaMq1(3,3), betaMq2(3,3), betaMu1(3,3), betaMu2(3,3), MqaYdYd(3,3) &
     & , MqaYuYu(3,3), aYdYdMq(3,3), aYuYuMq(3,3), aYeYeMl(3,3), MlaYeYe(3,3) &
     & , aYeMeYe(3,3), aYdMdYd(3,3), aYuMuYu(3,3)                             &
     & , YdaYdYdaYd(3,3), YeaYeYeaYe(3,3), YuaYuYuaYu(3,3), MeYeaYeYeaYe(3,3) &
     & , YeaYeYeaYeMe(3,3), YeaYeMeYeaYe(3,3), YeMlaYeYeaYe(3,3), AeaYe(3,3)  &
     & , YeaYeYeMlaYe(3,3), AeaAeYeaYe(3,3), YeaYeAeaAe(3,3), YeaAe(3,3)      &
     & , AeaYeYeaAe(3,3), YeaAeAeaYe(3,3), Tr3aAdYdaAeYe

  Complex(dp) :: Tr3aYdAdaYeAe, AdaYd(3,3), YdaAd(3,3), MlaYeYeaYeYe(3,3)     &
     & , aYeYeaYeYeMl(3,3), aYeYeMlaYeYe(3,3), aYeYeaYeMeYe(3,3)              &
     & , aYeMeYeaYeYe(3,3), aAdYd(3,3),aAeYe(3,3), aAeAeaYeYe(3,3)            &
     & , aYeYeaAeAe(3,3), aAeYeaYeAe(3,3), aYeAeaAeYe(3,3), MdYdaYdYdaYd(3,3) &
     & , YdaYdYdaYdMd(3,3), YdMqaYdYdaYd(3,3),YdaYdMdYdaYd(3,3)               &
     & , YdaYdYdMqaYd(3,3), AdaAdYdaYd(3,3), YdaYDAdaAd(3,3), AdaYdYdaAd(3,3) &
     & , YdaAdAdaYd(3,3)
  Complex(dp) :: MdYdaYuYuaYd(3,3), YdaYuYuaYdMd(3,3), YdMqaYuYuaYd(3,3)      &
     & , YdaYuYuMqaYd(3,3), YdaYuMuYuaYd(3,3), AdaAuYuaYd(3,3)                &
     & , YdaYuAuaAd(3,3), AdaYuYuaAd(3,3), YdaAuAuaYd(3,3), YdaYuYuaYd(3,3)   &
     & , Tr3aYuAu, Tr3aAuYu, YuaAu(3,3)                                       &
     & , MqaYdYdaYdYd(3,3), aYdYdaYdYDMq(3,3), aYdMdYdaYdYd(3,3)              &
     & , aYdYdMqaYdYd(3,3), aYdYdaYdMdYd(3,3), aAdAdaYdYd(3,3)                &
     & , aYdYDaAdAd(3,3), aAdYdaYdAd(3,3), aYdAdaAdYd(3,3), MqaYuYuaYuYu(3,3) &
     & , aYuYuaYuYUMq(3,3), aYuMuYuaYuYu(3,3), aYuYuMqaYuYu(3,3)              &
     & , aYuYuaYuMuYu(3,3), aAuAuaYuYu(3,3), aYuYUaAuAu(3,3), aAuYuaYuAu(3,3) &
     & , aYuAuaAuYu(3,3), aAuYu(3,3), AuaYu(3,3), YuaYdYdaYu(3,3)             &
     & , AuaYdYdaAu(3,3), YuaAdAdaYu(3,3), AuaAdYdaYu(3,3), YuaYdAdaAu(3,3)   &
     & , YuMqaYuYuaYu(3,3), YuaYuYuMqaYu(3,3), MuYuaYuYuaYu(3,3)              &
     & , YuaYuYuaYuMu(3,3), YuaYuMuYuaYu(3,3), AuaAuYuaYu(3,3)                &
     & , YuaYuAuaAu(3,3), AuaYuYuaAu(3,3), YuaAuAuaYu(3,3), MuYuaYdYdaYu(3,3) &
     & , YuaYdYdaYuMu(3,3), YuMqaYdYdaYu(3,3), YuaYdYdMqaYu(3,3)              &
     & , YuaYdMdYdaYu(3,3)
  Complex(dp) :: Mr(3,3), YnuaYnu(3,3), MlaYnuYnu(3,3), MrYnuaYnu(3,3)        &
     & , aYnuYnuMl(3,3), YnuaYnuMr(3,3), YnuMlaYnu(3,3), aYnuMrYnu(3,3)       &
     & , AnuaAnu(3,3), DMr(3,3), betaMr1(3,3)

  Real(dp) :: S1, S2, sig(3), Tr3aYdYdaYeYe, Tr3aAdAdaAeAe, AbsGM2(3)         &
     & , Tr3MqaYdYd3aYDMdYd , Tr3MqaYuYu3aYuMu, Tr3aAuAu

  Real(dp) :: Mh(2), DMh(2), TraceMH1(3), TraceMH2(2), betaMH11        &
     &  , betaMH12, betaMH21, betaMH22, q

  Complex(dp) :: mue, B, Dmue, DB, TraceMue(2), TraceB(2), betaMue1, betaMue2 &
     & , betaB1, betaB2, g2Mi(3)


  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge267'

  OnlyDiagonal = .Not.GenerationMixing
  q = t

  Call GToParameters2(gy, gauge, Ye, Ynu, Yd, Yu                            &
                & , Mhlf, Ae, Anu, Ad, Au, Me, Ml, Mr, Md, Mq, Mu, Mh, mue, B)

  gauge2 = gauge**2
  AbsGM2 = gauge2 * Abs( Mhlf )**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ynu,aYnu)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul(aYe,Ye)
  aYnuYnu = Matmul(aYnu,Ynu)
  aYdYd = Matmul(aYd,Yd)
  aYuYu = Matmul(aYu,Yu)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYnuYnu),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe + aYnuYnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = 3._dp * TraceY(4) + TraceY(2)           &
            &   - 0.6_dp * gauge2(1) - 3._dp * gauge2(2)
  sumnu1 = 3._dp * aYnuYnu + aYeYe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do

  betaYnu1 = Matmul(Ynu,sumnu1)

  diagonal(3,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)

  diagonal(4,1) = (3._dp,0._dp) * TraceY(4)  + TraceY(2)             &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
   aYdYdaYuYu = Matmul(aYdYd,aYuYu)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ae,aAe)
  Call Adjungate(Anu,aAnu)
  Call Adjungate(Ad,aAd)
  Call Adjungate(Au,aAu)

  aAdAd = Matmul(aAd,Ad)
  aAnuAnu = Matmul(aAnu,Anu)
  aAeAe = Matmul(aAe,Ae)
  aAuAu = Matmul(aAu,Au)

  TraceA(1) = Real( cTrace(aAeAe),dp )
  TraceA(2) = Real( cTrace(aAnuAnu),dp )
  TraceA(3) = Real( cTrace(aAdAd),dp )
  TraceA(4) = Real( cTrace(aAuAu),dp )

  aYdAd = Matmul(aYd,Ad)
  aYnuAnu = Matmul(aYnu,Anu)
  aYeAe = Matmul(aYe,Ae)
  aYuAu = Matmul(aYu,Au)

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYnuAnu) 
  TraceaYA(3) = cTrace(aYdAd) 
  TraceaYA(4) = cTrace(aYuAu) 

  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = Matmul(Ae,sume1)
  
  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1)  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    )
  sume1 = 4._dp * aYeAe + 2._dp * aYnuAnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + Matmul(Ye,sume1)

  !--------------
  ! A_nu
  !--------------
  sumnu1 = sumnu1 + 2._dp * aYnuYnu
  betaAnu1 = Matmul(Anu,sumnu1)

  diagonal(2,1) = 2._dp * ( 3._dp * TraceaYA(4) + TraceaYA(2)     &
      &         + 0.6_dp * g2Mi(1) + 3._dp * g2Mi(2) )
  sumnu1 = 4._dp * aYnuAnu + 2._dp * aYeAe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do
  betaAnu1 = betaAnu1 + Matmul(Ynu,sumnu1)

  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = Matmul(Ad,sumd1)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1)           &
                &         - c1_1(2,1) * g2Mi(1) - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do
  betaAd1 = betaAd1 + Matmul(Yd,sumd1)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = Matmul(Au,sumu1)
  
  diagonal(4,1) = 2._dp * ( 3._dp * TraceaYA(4) + TraceaYA(2)           &
                &         - c1_1(3,1) * g2Mi(1) - c1_1(3,2) * g2Mi(2)   &
                &         - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do
  betaAu1 = betaAu1 + Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdAd = Matmul(aYdYd,aYdAd)
   aYdAdaYdYd = Matmul(aYdAd,aYdYd)
   aYeYeaYeAe = Matmul(aYeYe,aYeAe)
   aYeAeaYeYe = Matmul(aYeAe,aYeYe)
   aYuYuaYuAu = Matmul(aYuYu,aYuAu)
   aYuAuaYuYu = Matmul(aYuAu,aYuYu)
   aYuAuaYdYd = Matmul(aYuAu,aYdYd)
   aYuYuaYdAd = Matmul(aYuYu,aYdAd)
   aYdAdaYuYu = Matmul(aYdAd,aYuYu)
   aYdYdaYuAu = Matmul(aYdYd,aYuAu)
   TraceAY2(1) = cTrace(aYeYeaYeAe)
   TraceAY2(2) = cTrace(aYdYdaYdAd)
   TraceAY2(3) = cTrace(aYuYuaYuAu)
   TraceAY2(4) = cTrace(aYuYuaYdAd)
   TraceAY2(5) = cTrace(aYdYdaYuAu)

  !--------------
  ! A_e
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = Matmul(Ae,sume2)
    
   diagonal(1,2) = -6._dp * ( 6._dp * TraceAY2(2) + TraceAY2(4)       &
     &                      + TraceAY2(5) + 2._dp * TraceAY2(1)  )    &
     &  + ( 32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)   &
     &  - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)         &
     &  + 2.4_dp * gauge2(1) * TraceaYA(1)                            &
     &  - 2.4_dp * g2Mi(1) * TraceY(1)                    &
     &  - ( 30._dp * g2Mi(2)                              &
     &    + 3.6_dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)      &
     &  - 54._dp * gauge2(1)**2 * Mhlf(1)
   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)     &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1) &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2)
   sume2 = - 6._dp * aYeYeaYeAe - 8._dp * aYeAeaYeYe  &
         & - hd(1) * aYeAe - hc(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
   betaAe2 = betaAe2 + Matmul(Ye,sume2)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = Matmul(Ad,sumd2)
    
   diagonal(3,2) = diagonal(1,2)                                   &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                 &
     &              - gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp      &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)    &
     &  + 1.6_dp * gauge2(1) * gauge2(2) * (Mhlf(1)+Mhlf(2))       & 
     &  + 1.856e3_dp * gauge2(1)**2 * Mhlf(1) / 4.5e1_dp

   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)  &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hd(2) = 6._dp * TraceY(4) - 1.6_dp * gauge2(1) 
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1)                 &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2) + 1.6_dp * gauge2(1) *  Mhlf(1)
   hc(2) = 6._dp * TraceaYA(4) + 1.6_dp * gauge2(1) *  Mhlf(1)
   sumd2 = - 6._dp * aYdYdaYdAd - 8._dp * aYdAdaYdYd                  &
       &   - 4._dp * ( aYuAuaYuYu + aYuYuaYuAu + aYuAuaYdYd )         &
       &   - 2._dp * aYuYuaYdAd - hd(1) * aYdAd - hc(1) * aYdYd       &
       &  - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(3,2)
   End Do
   betaAd2 = betaAd2 + Matmul(Yd,sumd2)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(4) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = Matmul(Au,sumu2)
    
   diagonal(4,2) =  -6._dp * ( 6._dp * TraceAY2(3) + TraceAY2(4)        &
     &                        + TraceAY2(5)  )                          &
     &  + ( 32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)     &
     &  - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)           &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                      &
     &              - 3.4_dp * gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp  &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)         &
     &  - ( 30._dp * g2Mi(2)                                &
     &    + 2._dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)         &
     &  - 5486._dp * gauge2(1)**2 * Mhlf(1) / 225._dp
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1) - 0.8_dp * gauge2(1)
   hc(1) = 6._dp * TraceaYA(3) + 2._dp * TraceaYA(1)   &
       &  + 0.8_dp * gauge2(1) *  Mhlf(1)
   hd(2) = 12._dp * TraceY(4) - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(2) = 18._dp * TraceaYA(4) + 1.2e1_dp * g2Mi(2) + 0.8_dp * g2Mi(1)
   sumu2 = - 6._dp * aYuYuaYuAu - 8._dp * aYuAuaYuYu                     &
       &   - 4._dp * ( aYdAdaYdYd + aYdYdaYdAd + aYdAdaYuYu )            &
       &   - 2._dp * aYdYdaYuAu - hd(1) * aYdAd - hc(1) * aYdYd          &
       &   - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(4,2)
   End Do
   betaAu2 = betaAu2 + Matmul(Yu,sumu2)

  End If 
!----------------------------------------------
! beta functions for Sfermion mass parameters
!----------------------------------------------
   S1 = mH(2) - mH(1)
   Do i1=1,3
    S1 = S1 + Me(i1,i1) - Ml(i1,i1) &
       &    + Md(i1,i1) + Mq(i1,i1) - 2._dp * Mu(i1,i1)
   End Do
   S1 = S1 * gauge2(1)

   YdaYd = Matmul(Yd,aYd)
   YnuaYnu = Matmul(Ynu,aYnu)
   YeaYe = Matmul(Ye,aYe)
   YuaYu = Matmul(Yu,aYu)

   MeYeaYe = Matmul(Me,YeaYe)
   MlaYeYe = Matmul(Ml,aYeYe)
   MlaYnuYnu = Matmul(Ml,aYnuYnu)
   MrYnuaYnu = Matmul(Mr,YnuaYnu)

   MdYdaYd = Matmul(Md,YdaYd)
   MqaYdYd = Matmul(Mq,aYdYd)
   MqaYuYu = Matmul(Mq,aYuYu)
   MuYuaYu = Matmul(Mu,YuaYu)

   YeaYeMe = Matmul(YeaYe,Me)
   aYeYeMl = Matmul(aYeYe,Ml)
   aYnuYnuMl = Matmul(aYnuYnu,Ml)
   YnuaYnuMr = Matmul(YnuaYnu,Mr)

   YdaYdMd = Matmul(YdaYd,Md)
   aYdYdMq = Matmul(aYdYd,Mq)
   aYuYuMq = Matmul(aYuYu,Mq)
   YuaYuMu = Matmul(YuaYu,Mu)

   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   YnuMlaYnu = MatMul3(Ynu,Ml,aYnu,OnlyDiagonal)
   aYnuMrYnu = MatMul3(aYnu,Mr,Ynu,OnlyDiagonal)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)

   AeaAe = Matmul(Ae,aAe)
   AnuaAnu = Matmul(Anu,aAnu)
   AdaAd = Matmul(Ad,aAd)
   AuaAu = Matmul(Au,aAu)

   diagonal(1,1) = - 4.8_dp * AbsGM2(1) + 1.2_dp * S1
   betaMe1 = 2._dp * (MeYeaYe + YeaYeMe)             &
         & + 4._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe )
   Do i1=1,3
    betaMe1(i1,i1) = betaMe1(i1,i1) + diagonal(1,1)
   End Do

   diagonal(3,1) = - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1
   betaMl1 = MlaYeYe + aYeYeMl + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe )  &
         & + MlaYnuYnu + aYnuYnuMl                                            &
         & + 2._dp * ( mH(2) * aYnuYnu + aYnuMrYnu + aAnuAnu )
   Do i1=1,3
    betaMl1(i1,i1) = betaMl1(i1,i1) + diagonal(3,1)
   End Do

   betaMr1 = 2._dp * (MrYnuaYnu + YnuaYnuMr)                        &
         & + 4._dp * ( mH(2) * YnuaYnu + YnuMlaYnu + AnuaAnu )

   diagonal(4,1) = - ( 32._dp * AbsGM2(3) + 1.6_dp * AbsGM2(1) ) / 3._dp &
               & + 0.4_dp * S1
   betaMd1 = 2._dp * (MdYdaYd + YdaYdMd)             &
         & + 4._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd )
   Do i1=1,3
    betaMd1(i1,i1) = betaMd1(i1,i1) + diagonal(4,1)
   End Do

   diagonal(5,1) = - ( 32._dp * AbsGM2(3) + 0.4_dp * AbsGM2(1) ) / 3._dp &
               & - 6._dp * AbsGM2(2) + 0.2_dp * S1
   betaMq1 = MqaYuYu + aYuYuMq + MqaYdYd + aYdYdMq             &
         & + 2._dp * ( mH(2) * aYuYu + mH(1) * aYdYd + aYuMuYu     &
         &           + aYdMdYd + aAuAu + aAdAd )
   Do i1=1,3
    betaMq1(i1,i1) = betaMq1(i1,i1) + diagonal(5,1)
   End Do

   diagonal(6,1) = - ( 32._dp * AbsGM2(3) + 6.4_dp * AbsGM2(1) ) / 3._dp &
               &   - 0.8_dp * S1
   betaMu1 = 2._dp * (MuYuaYu + YuaYuMu)             &
         & + 4._dp * ( mH(2) * YuaYu + YuMqaYu + AuaAu )
   Do i1=1,3
    betaMu1(i1,i1) = betaMu1(i1,i1) + diagonal(6,1)
   End Do

   If (TwoLoopRGE) Then
    YdaYdYdaYd = MatSquare(YdaYd,OnlyDiagonal)
    YeaYeYeaYe = MatSquare(YeaYe,OnlyDiagonal)
    YuaYuYuaYu = MatSquare(YuaYu,OnlyDiagonal)

    AdaYd = Matmul(Ad,aYd)
    AeaYe = Matmul(Ae,aYe)
    AuaYu = Matmul(Au,aYu)

    aAdYd = Matmul(aAd,Yd)
    aAeYe = Matmul(aAe,Ye)
    aAuYu = Matmul(aAu,Yu)

    YdaAd = Matmul(Yd,aAd)
    YeaAe = Matmul(Ye,aAe)
    YuaAu = Matmul(Yu,aAu)

    YdaYuYuaYd = MatMul3(Yd,aYuYu,aYd,OnlyDiagonal)
    AdaYuYuaAd = MatMul3(Ad,aYuYu,aAd,OnlyDiagonal)
    YdaAuAuaYd = MatMul3(Yd,aAuAu,aYd,OnlyDiagonal)
    AdaAuYuaYd = MatMul4(Ad,aAu,Yu,aYd,OnlyDiagonal)
    YdaYuAuaAd = MatMul3(Yd,aYuAu,aAd,OnlyDiagonal)

    YuaYdYdaYu = MatMul3(Yu,aYdYd,aYu,OnlyDiagonal)
    AuaYdYdaAu = MatMul3(Au,aYdYd,aAu,OnlyDiagonal)
    YuaAdAdaYu = MatMul3(Yu,aAdAd,aYu,OnlyDiagonal)
    AuaAdYdaYu = MatMul4(Au,aAd,Yd,aYu,OnlyDiagonal)
    YuaYdAdaAu = MatMul3(Yu,aYdAd,aAu,OnlyDiagonal)

    MdYdaYuYuaYd = Matmul(Md,YdaYuYuaYd)
    Call Adjungate(MdYdaYuYuaYd, YdaYuYuaYdMd)
    YdMqaYuYuaYd = MatMul3(Yd,MqaYuYu,aYd,OnlyDiagonal)
    Call Adjungate(YdMqaYuYuaYd, YdaYuYuMqaYd)
    YdaYuMuYuaYd = MatMul3(Yd,aYuMuYu,aYd,OnlyDiagonal)

    MuYuaYdYdaYu = Matmul(Mu,YuaYdYdaYu)
    Call Adjungate(MuYuaYdYdaYu, YuaYdYdaYuMu)
    YuMqaYdYdaYu = MatMul3(Yu,MqaYdYd,aYu,OnlyDiagonal)
    Call Adjungate(YuMqaYdYdaYu, YuaYdYdMqaYu)
    YuaYdMdYdaYu = MatMul3(Yu,aYdMdYd,aYu,OnlyDiagonal)

    MeYeaYeYeaYe = Matmul(MeYeaYe,YeaYe)
    Call Adjungate(MeYeaYeYeaYe,YeaYeYeaYeMe)
    aYeMeYeaYeYe = Matmul(aYeMeYe,aYeYe)
    Call Adjungate(aYeMeYeaYeYe,aYeYeaYeMeYe)
    YeaYeMeYeaYe = Matmul(YeaYeMe,YeaYe)

    MlaYeYeaYeYe = Matmul(MlaYeYe,aYeYe)
    Call Adjungate(MlaYeYeaYeYe, aYeYeaYeYeMl)
    YeMlaYeYeaYe = Matmul(YeMlaYe,YeaYe)
    Call Adjungate(YeMlaYeYeaYe, YeaYeYeMlaYe)
    aYeYeMlaYeYe = Matmul(aYeYeMl,aYeYe)

    MdYdaYdYdaYd = Matmul(MdYdaYd,YdaYd)
    Call Adjungate(MdYdaYdYdaYd, YdaYdYdaYdMd)
    aYdMdYdaYdYd = Matmul(aYdMdYd,aYdYd)
    Call Adjungate(aYdMdYdaYdYd, aYdYdaYdMdYd)
    YdaYdMdYdaYd = Matmul(YdaYdMd,YdaYd)

    MqaYdYdaYdYd = Matmul(MqaYdYd,aYdYd)
    Call Adjungate(MqaYdYdaYdYd, aYdYdaYdYdMq)
    YdMqaYdYdaYd = Matmul(YdMqaYd,YdaYd)
    Call Adjungate(YdMqaYdYdaYd, YdaYdYdMqaYd)
    aYdYdMqaYdYd = Matmul(aYdYdMq,aYdYd)

    MqaYuYuaYuYu = Matmul(MqaYuYu,aYuYu)
    Call Adjungate(MqaYuYuaYuYu, aYuYuaYuYuMq)
    YuMqaYuYuaYu = Matmul(YuMqaYu,YuaYu)
    Call Adjungate(YuMqaYuYuaYu, YuaYuYuMqaYu)
    aYuYuMqaYuYu = Matmul(aYuYuMq,aYuYu)

    MuYuaYuYuaYu = Matmul(MuYuaYu,YuaYu)
    Call Adjungate(MuYuaYuYuaYu, YuaYuYuaYuMu)
    aYuMuYuaYuYu = Matmul(aYuMuYu,aYuYu)
    Call Adjungate(aYuMuYuaYuYu, aYuYuaYuMuYu)
    YuaYuMuYuaYu = Matmul(YuaYuMu,YuaYu)

    AdaAdYdaYd = Matmul(AdaAd,YdaYd)
    Call Adjungate(AdaAdYdaYd, YdaYdAdaAd)
    AdaYdYdaAd = Matmul(AdaYd,YdaAd)
    YdaAdAdaYd = Matmul(YdaAd,AdaYd)

    aAdAdaYdYd = Matmul(aAdAd,aYdYd)
    Call Adjungate(aAdAdaYdYd, aYdYdaAdAd)
    aAdYdaYdAd = Matmul(aAdYd,aYdAd)
    aYdAdaAdYd = Matmul(aYdAd,aAdYd)

    AeaAeYeaYe = Matmul(AeaAe,YeaYe)
    Call Adjungate(AeaAeYeaYe, YeaYeAeaAe)
    AeaYeYeaAe = Matmul(AeaYe,YeaAe)
    YeaAeAeaYe = Matmul(YeaAe,AeaYe)

    aAeAeaYeYe = Matmul(aAeAe,aYeYe)
    Call Adjungate(aAeAeaYeYe, aYeYeaAeAe)
    aAeYeaYeAe = Matmul(aAeYe,aYeAe)
    aYeAeaAeYe = Matmul(aYeAe,aAeYe)

    AuaAuYuaYu = Matmul(AuaAu,YuaYu)
    Call Adjungate(AuaAuYuaYu, YuaYuAuaAu)
    AuaYuYuaAu = Matmul(AuaYu,YuaAu)
    YuaAuAuaYu = Matmul(YuaAu,AuaYu)

    aAuAuaYuYu = Matmul(aAuAu,aYuYu)
    Call Adjungate(aAuAuaYuYu, aYuYuaAuAu)
    aAuYuaYuAu = Matmul(aAuYu,aYuAu)
    aYuAuaAuYu = Matmul(aYuAu,aAuYu)

    S2 = (1.5_dp * gauge2(2) + 0.3_dp * gauge2(1) )            &
     &      * (MH(2) - MH(1) - Real(cTrace(ML),dp) )              &
     & + ( (8._dp * gauge2(3) + 0.1_dp*gauge2(1)) / 3._dp      &
     &   + 1.5_dp * gauge2(2) ) * Real( cTrace(Mq),dp   )         &
     & - (16._dp * gauge2(3) + 3.2_dp*gauge2(1) )              &
     &    * Real(cTrace(Mu),dp) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md),dp) / 3._dp                           &
     & + 1.2_dp*gauge2(1) * Real(cTrace(Me),dp)                   &
     & - 3._dp * (MH(2)*TraceY(4) - MH(1) * TraceY(3) )        &
     & + mH(1) * TraceY(1)
    
    Do i1=1,3
     S2 = S2 - YuMqaYu(i1,i1) + 4._dp * aYuMuYu(i1,i1)   &
        &    - YdMqaYd(i1,i1) - 2._dp * aYdMdYd(i1,i1)   &
        &    + YeMlaYe(i1,i1) - 2._dp * aYeMeYe(i1,i1)
    End Do

    sig(1) = 3._dp * (MH(1) + MH(2) + Real(cTrace(Ml),dp) )          &
         & + Real(cTrace(Mq),dp) + 8._dp * Real(cTrace(Mu),dp)          &
         & + 2._dp * Real(cTrace(Md),dp) + 6._dp * Real(cTrace(Me),dp) 
    sig(1) = 0.2_dp * gauge2(1) * sig(1)
    sig(2) = gauge2(2) * ( MH(1) + MH(2) + Real( cTrace(Ml),dp )    &
           &             + 3._dp * Real( cTrace(Mq),dp ) )
    sig(3) = gauge2(3) * ( 2._dp * Real( cTrace(Mq),dp )            &
           &             + Real(cTrace(Mu),dp) + Real(cTrace(Md),dp) )

    Tr3aYdYdaYeYe = 3._dp * TraceY(3) + TraceY(1)
    Tr3MqaYdYd3aYDMdYd = Real(cTrace(MlaYeYe),dp) + Real(cTrace(aYeMeYe),dp)  &
            & + 3._dp * ( Real(cTrace(MqaYdYd),dp) + Real(cTrace(aYdMdYd),dp) )
    Tr3aAdAdaAeAe = 3._dp * TraceA(3) + TraceA(1)
    Tr3aYdAdaYeAe = 3._dp * TraceaYA(3) + TraceaYA(1)
    Tr3aAdYdaAeYe = Conjg( Tr3aYdAdaYeAe )
    Tr3MqaYuYu3aYuMu = 3._dp * (Real(cTrace(MqaYuYu),dp) &
                               + Real(cTrace(aYuMuYu),dp) )
    Tr3aAuAu = 3._dp * TraceA(4)
    Tr3aYuAu = 3._dp * TraceaYA(4)
    Tr3aAuYu = Conjg( Tr3aYuAu )

    diagonal(1,2) = 2.4_dp * gauge2(1) * (S2 + sig(1) )   &
                & + 112.32_dp * gauge2(1) * AbsGM2(1)
    hd(1) = 6._dp*gauge2(2)-1.2_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 4.8_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 2.4_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMe2 = -2._dp * ( MeYeaYeYeaYe + YeaYeYeaYeMe )                       &
      & - 4._dp * ( YeMlaYeYeaYe + YeaYeMeYeaYe + YeaYeYeMlaYe )             &
      & - 8._dp * mH(1) * YeaYeYeaYe                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MeYeaYe + YeaYeMe + 4._dp * MH(1) * YeaYe &
      &                           + 2._dp * (YeMlaYe + AeaAE)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YeaYe + AeaAeYeaYe + YeaYEAeaAe     &
      &           + AeaYeYeaAe + YeaAeAeaYe + Tr3aAdAdaAeAe * YeaYe          &
      &           + Tr3aAdYdaAeYe * AeaYe + Tr3aYdAdaYeAe * YeaAe )          &
      & + hd(1) * ( MeYeaYe + YeaYeMe                                        &
      &           + 2._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe ) )            &
      & + hd(2) * YeaYe + hc(1) * AeaYe + hc(2) * YeaAe
    Do i1=1,3
     betaMe2(i1,i1) = betaMe2(i1,i1) + diagonal(1,2)
    End Do

    diagonal(3,2) = gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)        &
     &     + 3._dp * gauge2(2) * sig(2)                                 &
     &     + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)      &
     &                   + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2       &
     &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )     &
     &     + 24.84_dp * gauge2(1) * AbsGM2(1)
    betaMl2 = -2._dp * ( MlaYeYeaYeYe + aYeYeaYeYeMl )                   &
      & - 4._dp * ( aYeMeYeaYeYe + aYeYeMlaYeYe + aYeYeaYeMeYe )         &
      & - 8._dp * mH(1) * aYeYeaYeYe                                     &
      & - Tr3aYdYdaYeYe * ( MlaYeYe + aYeYeML + 4._dp * MH(1) * aYeYe    &
      &                   + 2._dp * (aYeMeYe + aAeAE)  )                 &
      & - 4._dp * (aAeAeaYeYe + aYeYeaAeAe + aAeYeaYeAe + aYeAeaAeYe )   &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYeYe +  Tr3aAdAdaAeAe * aYeYe  &
      &           + Tr3aAdYdaAeYe * aYeAe + Tr3aYdAdaYeAe * aAeYe     )  &
      & + 1.2_dp*gauge2(1) * ( MlaYeYe + aYeYeMl                         &
      &                      + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe &
      &                                - Mhlf(1) * aAeYe                 &
      &                                - Conjg(Mhlf(1)) * aYeAe ) )      &
      & + 4.8_dp * AbsGM2(1) * aYeYe
    Do i1=1,3
     betaMl2(i1,i1) = betaMl2(i1,i1) + diagonal(3,2)
    End Do

    diagonal(4,2) = 0.8_dp * gauge2(1) * (S2 + sig(1)/3._dp )             &
      &   + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp   &
      &   + 8.08e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                    &
      &   + 1.28e2_dp * ( gauge2(1) * AbsGM2(3)                           &
      &                 + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)     &
      &                   * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)+0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) + 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) - 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMd2 = -2._dp * ( MdYdaYdYdaYd + YdaYdYdaYdMd )                       &
      & - 4._dp * ( YdMqaYdYdaYd + YdaYdMdYdaYd + YdaYdYdMqaYd )             &
      & - 8._dp * mH(1) * YdaYdYdaYd                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MdYdaYd + YdaYdMd + 4._dp * MH(1) * YdaYd &
      &                           + 2._dp * (YdMqaYd + AdaAD)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YdaYd + AdaAdYdaYd + YdaYDAdaAd     &
      &           + AdaYdYdaAd + YdaAdAdaYd + Tr3aAdAdaAeAe * YdaYd          &
      &           + Tr3aAdYdaAeYe * AdaYd + Tr3aYdAdaYeAe * YdaAd )          &
      & + hd(1) * ( MdYdaYd + YdaYdMd                                        &
      &           + 2._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd ) )            &
      & + hd(2) * YdaYd + hc(1) * AdaYd + hc(2) * YdaAd                      &
      & - 2._dp * (MdYdaYuYuaYd + YdaYuYuaYdMd )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YdaYuYuaYd + YdMqaYuYuaYd+ YdaYuYuMqaYd  &
      &           + YdaYuMuYuaYd + AdaAuYuaYd + YdaYuAuaAd                   &
      &           + AdaYuYuaAd + YdaAuAuaYd )
    Do i1=1,3
     betaMd2(i1,i1) = betaMd2(i1,i1) + diagonal(4,2)
    End Do

    diagonal(5,2) = 0.2_dp * gauge2(1) * (2._dp * S2 + sig(1)/3._dp )     &
      &  + 3._dp * gauge2(2) * sig(2)                                     &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp    &
      &  + 1.99e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                     &
      &  + 32._dp * ( gauge2(1) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)         &
      &              * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp      &
      &  + 33._dp * gauge2(2) * AbsGM2(2)                                 &
      &  + 32._dp * ( gauge2(2) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(2)  + gauge2(3)*gauge2(2)        &
      &              * Real( Mhlf(2) * Conjg(Mhlf(3)),dp ) )                 &
      &  + 0.4_dp * ( gauge2(2) * AbsGM2(1)                               &
      &             + gauge2(1) * AbsGM2(2) + gauge2(1)*gauge2(2)         &
      &              * Real( Mhlf(2) * Conjg(Mhlf(1)),dp ) )
    hd(1) = 1.6_dp * AbsGM2(1)
    hc(1) = - 0.8_dp * g2Mi(1)
    hc(2) = Conjg( hc(1) )
    hd(2) = 3.2_dp * AbsGM2(1)
    hc(3) = - 1.6_dp * g2Mi(1)
    hc(4) = Conjg( hc(3) )
    betaMq2 = -2._dp * ( MqaYdYdaYdYd + aYdYdaYdYDMq )                       &
      & - 4._dp * ( aYdMdYdaYdYd + aYdYdMqaYdYd + aYdYdaYdMdYd )             &
      & - 8._dp * mH(1) * aYdYdaYdYd                                         &
      & - Tr3aYdYdaYeYe * ( MqaYdYd + aYdYdMq + 4._dp * MH(1) * aYdYd        &
      &                   + 2._dp * (aYdMdYd + aAdAd)  )                     &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYdYd                               &
      &           + 2._dp * ( aAdAdaYdYd + aYdYDaAdAd + aAdYdaYdAd           &
      &                     + aYdAdaAdYd )                                   &
      &           + Tr3aAdAdaAeAe * aYdYd + Tr3aAdYdaAeYe * aYdAd            &
      &           + Tr3aYdAdaYeAe * aAdYd )                                  &
      & + 0.4_dp*gauge2(1) * ( MqaYdYd + aYdYdMq                             &
      &                      + 2._dp * ( mH(1) * aYdYd + aYdMdYd + aAdAd ) ) &
      & + hd(1) * aYdYd + hc(1) * aAdYd + hc(2) * aYdAd                      &
      & - 2._dp * ( MqaYuYuaYuYu + aYuYuaYuYuMq )                            &
      & - 4._dp * ( aYuMuYuaYuYu + aYuYuMqaYuYu + aYuYuaYuMuYu )             &
      & - 8._dp * MH(2) * aYuYuaYuYu                                         &
      & - 3._dp * TraceY(4) * ( MqaYuYu + aYuYuMq + 4._dp * MH(2) * aYuYu    &
      &                       + 2._dp * (aYuMuYu + aAuAu) )                  &
      & - 2._dp * ( Tr3MqaYuYu3aYuMu * aYuYu                                 &
      &           + 2._dp * ( aAuAuaYuYu + aYuYuaAuAu                        &
      &                     + aAuYuaYuAu + aYuAuaAuYu )                      &
      &           + Tr3aAuAu * aYuYu + Tr3aAuYu * aYuAu + Tr3aYuAu * aAuYu ) &
      & + 0.8_dp*gauge2(1) * ( MqaYuYu + aYuYuMq                             &
      &                      + 2._dp * ( MH(2) * aYuYu + aYuMuYu + aAuAu ) ) &
      & + hd(2) * aYuYu + hc(3) * aAuYu + hc(4) * aYuAu
    Do i1=1,3
     betaMq2(i1,i1) = betaMq2(i1,i1) + diagonal(5,2)
    End Do

    diagonal(6,2) = 1.6_dp * gauge2(1) * (2._dp*sig(1)/3._dp - S2)           &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp       &
      &  + 3424._dp * gauge2(1) * AbsGM2(1) / 75._dp                         &
      &  + 512._dp * ( gauge2(1) * AbsGM2(3)                                 &
      &              + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)           &
      &                * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)-0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMu2 = -2._dp * ( MuYuaYuYuaYu + YuaYuYuaYuMu )                       &
      & - 4._dp * ( YuMqaYuYuaYu + YuaYuMuYuaYu + YuaYuYuMqaYu )             &
      & - 8._dp * MH(2) * YuaYuYuaYu                                         &
      & - 6._dp * TraceY(4) * ( MuYuaYu + YuaYuMu + 4._dp * MH(2) * YuaYu    &
      &                       + 2._dp * (YuMqaYu + AuaAu)  )                 &
      & - 4._dp * ( Tr3MqaYuYu3aYuMu * YuaYu + AuaAuYuaYu + YuaYuAuaAu       &
      &           + AuaYuYuaAu + YuaAuAuaYu + Tr3aAuAu * YuaYu               &
      &           + Tr3aAuYu * AuaYu + Tr3aYuAu * YuaAu )                    &
      & + hd(1) * ( MuYuaYu + YuaYuMu                                        &
      &           + 2._dp * ( MH(2) * YuaYu + YuMqaYu + AuaAu ) )            &
      & + hd(2) * YuaYu + hc(1) * AuaYu + hc(2) * YuaAu                      &
      & - 2._dp * (MuYuaYdYdaYu + YuaYdYdaYuMu )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YuaYdYdaYu + YuMqaYdYdaYu + YuaYdYdMqaYu &
      &           + YuaYdMdYdaYu + AuaAdYdaYu + YuaYdAdaAu                   &
      &           + AuaYdYdaAu + YuaAdAdaYu )
    Do i1=1,3
     betaMu2(i1,i1) = betaMu2(i1,i1) + diagonal(6,2)
    End Do

   End If 

  !------------------------------------------
  ! beta functions for Higgs mass parameters
  !------------------------------------------
   traceMH1(1) = mH(1) * TraceY(1) + Real( cTrace(YeMlaYe),dp ) &
             & + Real( cTrace(aYeMeYe),dp ) + TraceA(1)
   traceMH1(2) = mH(1) * TraceY(3) + Real( cTrace(YdMqaYd),dp ) &
             & + Real( cTrace(aYdMdYd),dp ) + TraceA(3)
   betamH11 = 6._dp * TraceMH1(2) + 2._dp * TraceMH1(1)      &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1

   traceMH2(2) = mH(2) * TraceY(2) + Real( cTrace(YnuMlaYnu),dp )  &
             & + Real( cTrace(aYnuMrYnu),dp ) + TraceA(2)
   traceMH2(1) = mH(2) * TraceY(4) + Real( cTrace(YuMqaYu),dp )  &
             & + Real( cTrace(aYuMuYu),dp ) + TraceA(4)
   betamH21 = 2._dp * TraceMH2(2) + 6._dp * TraceMH2(1)       &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) + 0.6_dp * S1

   If (TwoLoopRGE) Then
    traceMH1(3) = MH(1) * (6._dp*TraceY2(2) + 2._dp*TraceY2(1) + TraceY2(4) ) &
              & + MH(2) * TraceY2(4)                                          &
              & + 6._dp * ( Real( cTrace(MqaYdYdaYdYd),dp )                   &
              &           + Real( cTrace(aYdMdYdaYdYd),dp )                   &
              &           + Real( cTrace(aAdAdaYdYd),dp )                     &
              &           + Real( cTrace(aAdYdaYdAd),dp )  )                  &
              & + 2._dp * ( Real( cTrace(MlaYeYeaYeYe),dp )                   &
              &           + Real( cTrace(aYeMeYeaYeYe),dp )                   &
              &           + Real( cTrace(aAeAeaYeYe),dp )                     &
              &           + Real( cTrace(aAeYeaYeAe),dp )  )                  &
       & + Real( cTrace(YdMqaYuYuaYd),dp ) + Real( cTrace(YdaYuMuYuaYd),dp ) &
       & + Real( cTrace(YdaYuYuMqaYd),dp ) + Real( cTrace(YuaYdMdYdaYu),dp ) &
       & + Real( cTrace(YdaAuAuaYd),dp ) + Real( cTrace(AdaYuYuaAd),dp )     &
              & + Real( cTrace(AdaAuYuaYd),dp ) + Real( cTrace(YdaYuAuaAd),dp )
    betaMH12 = - 6._dp * traceMH1(3)                                       &
      &   + (32._dp*gauge2(3) - 0.8_dp*gauge2(1) ) * traceMH1(2)           &
      &   + 64._dp * ( AbsGM2(3) * TraceY(3)                               &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(3),dp ) )     &
      &   - 1.6_dp * ( AbsGM2(1) * TraceY(3)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(3),dp ) )     &
      &   + 2.4_dp*gauge2(1) * traceMH1(1)                                 &
      &   + 4.8_dp * ( AbsGM2(1) * TraceY(1)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(1),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)                   &
      &   + 3._dp * gauge2(2) * sig(2)                                     &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)          &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2           &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

    traceMH2(2) = MH(2) * ( 6._dp * TraceY2(3) + TraceY2(4) )                 &
      &  + MH(1) * TraceY2(4)                                                 &
 &  + 6._dp * ( Real(cTrace(MqaYuYuaYuYu),dp) + Real(cTrace(aYuMuYuaYuYu),dp) &
 &            + Real(cTrace(aAuAuaYuYu),dp) + Real(cTrace(aAuYuaYuAu),dp)  )  &
 &  + Real( cTrace(YuMqaYdYdaYu),dp ) + Real( cTrace(YuaYdMdYdaYu),dp )       &
 &  + Real( cTrace(YuaYdYdMqaYu),dp ) + Real( cTrace(YdaYuMuYuaYd),dp )       &
 &  + Real( cTrace(YuaAdAdaYu),dp ) + Real( cTrace(AuaYdYdaAu),dp )           &
      &  + Real( cTrace(AuaAdYdaYu),dp ) + Real( cTrace(YuaYdAdaAu),dp )
    betaMH22 = - 6._dp * traceMH2(2)                                          &
      &   + (32._dp*gauge2(3) + 1.6_dp*gauge2(1) ) * traceMH2(1)              &
      &   + 64._dp * ( AbsGM2(3) * TraceY(4)                                  &
      &          - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(4),dp ) )        &
      &   + 3.2_dp * ( AbsGM2(1) * TraceY(4)                                  &
      &          - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(4),dp ) )        &
      &   + gauge2(1) * ( 0.6_dp * sig(1) + 1.2_dp * S2)                      &
      &   + 3._dp * gauge2(2) * sig(2)                                        &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)             &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2              &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

   End If
!-----------------------------
! beta functions for mu and B
!-----------------------------
   TraceMue(1) = 3._dp * (TraceY(3)+TraceY(4)) + TraceY(1) + TraceY(2)  &
             & - 3._dp * gauge2(2) - 0.6_dp * gauge2(1)
   betaMue1 = mue * TraceMue(1)

   TraceB(1) = 6._dp * (TraceaYA(3)+TraceaYA(4))                          &
           & + 2._dp * (TraceaYA(1)+ TraceaYA(2))                         &
           & + 6._dp * g2Mi(2) + 1.2_dp * g2Mi(1)
   betaB1 = mue * TraceB(1) + B * TraceMue(1)

   If (TwoLoopRGE) Then
    TraceMue(2) = - 3._dp * ( 3._dp * (TraceY2(2) + TraceY2(3) )          &
      &                     + 2._dp * TraceY2(4) + TraceY2(1) )           &
      &         + (16._dp * gauge2(3) + 0.8_dp * gauge2(1)) * TraceY(4)   &
      &         + (16._dp * gauge2(3) - 0.4_dp * gauge2(1)) * TraceY(3)   &
      &         + 1.2_dp * gauge2(1) * TraceY(1)                          &
      &         + 7.5_dp * gauge2(2)**2                                   &
      &         + 1.8_dp * gauge2(2) * gauge2(1)                          &
      &         + 4.14_dp * gauge2(1)**2
    betaMue2 = mue * TraceMue(2)

    TraceB(1) = cTrace( 3._dp * ( Matmul(AuaYu,YuaYu) + Matmul(AdaYd,YdaYd) ) &
              &       + Matmul(AeaYe,YeaYe) + Matmul(aYuAu,aYdYd)             &
              &       + Matmul(aYdAd,aYuYu) ) 
    TraceB(2) = -12._dp * TraceB(1)                                           &
      &   + (32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)          &
      &   + (32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)          &
      &   + 2.4_dp * gauge2(1) * TraceaYA(1)                                  &
      &   - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)               &
      &   - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)               &
      &   - 2.4_dp * g2Mi(1) * TraceY(1) - 30._dp * gauge2(2)**2 * Mhlf(2)    &
      &   - 3.6_dp * gauge2(2) * gauge2(1) * (Mhlf(1) + Mhlf(2) )             &
      &   - 16.56_dp * gauge2(1)**2 * Mhlf(1) 
    betaB2 = mue * TraceB(2) + B * TraceMue(2)

   End If

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1 + oo16pi2 * (Matmul(b_2,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYnu = oo16pi2 * betaYnu1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   Do i1 = 1,3    
    sumI = 0._dp
    Do i2=1,3
     sumI = sumI + b_2(i1,i2) * gauge2(i2) * (Mhlf(i1) + Mhlf(i2) ) 
    End Do
    Do i2=1,4
     sumI = sumI + a_2a(i1,i2) * ( TraceaYA(i2) - Mhlf(i1)*TraceY(i2) )
    End Do
    DMhlf(i1) = oo8pi2 * gauge2(i1) * ( b_1(i1) * Mhlf(i1) + oo16pi2 * sumI)
   End Do
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * ( betaAe1 + oo16pi2 * betaAe2 )
   DAnu = oo16pi2 * betaAnu1
   DAd = oo16pi2 * ( betaAd1 + oo16pi2 * betaAd2 )
   DAu = oo16pi2 * ( betaAu1 + oo16pi2 * betaAu2 )
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * ( betaMe1 + oo16pi2 * betaMe2 )
   DMl = oo16pi2 * ( betaMl1 + oo16pi2 * betaMl2 )
   DMr = oo16pi2 * betaMr1 ! + oo16pi2 * betaMl2 )
   DMd = oo16pi2 * ( betaMd1 + oo16pi2 * betaMd2 )
   DMq = oo16pi2 * ( betaMq1 + oo16pi2 * betaMq2 )
   DMu = oo16pi2 * ( betaMu1 + oo16pi2 * betaMu2 )
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * ( betaMH11 + oo16pi2 * betaMH12 )
   DmH(2) = oo16pi2 * ( betaMH21 + oo16pi2 * betaMH22 )
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * ( betaMue1 + oo16pi2 * betaMue2 )
   DB = oo16pi2 * ( betaB1 + oo16pi2 * betaB2 )

 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYnu = oo16pi2 * betaYnu1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   DMhlf = oo8pi2 * gauge2 * b_1 * Mhlf
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * betaAe1
   DAnu = oo16pi2 * betaAnu1
   DAd = oo16pi2 * betaAd1
   DAu = oo16pi2 * betaAu1
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * betaMe1
   DMl = oo16pi2 * betaMl1
   DMr = oo16pi2 * betaMr1
   DMd = oo16pi2 * betaMd1
   DMq = oo16pi2 * betaMq1
   DMu = oo16pi2 * betaMu1
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * betaMH11
   DmH(2) = oo16pi2 * betaMH21
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * betaMue1
   DB = oo16pi2 * betaB1
  End If


  !---------------------------------------
  ! to avoid numerical problems in odeint
  !---------------------------------------
  Do i1=1,3
   DMe(i1,i1) = Real(DMe(i1,i1),dp)
   DMr(i1,i1) = Real(DMr(i1,i1),dp)
   DMl(i1,i1) = Real(DMl(i1,i1),dp)
   DMd(i1,i1) = Real(DMd(i1,i1),dp)
   DMu(i1,i1) = Real(DMu(i1,i1),dp)
   DMq(i1,i1) = Real(DMq(i1,i1),dp)
  End Do
  Dmd = 0.5_dp * ( Dmd + Transpose(Conjg(Dmd)) )
  Dme = 0.5_dp * ( Dme + Transpose(Conjg(Dme)) )
  Dml = 0.5_dp * ( Dml + Transpose(Conjg(Dml)) )
  Dmq = 0.5_dp * ( Dmq + Transpose(Conjg(Dmq)) )
  Dmu = 0.5_dp * ( Dmu + Transpose(Conjg(Dmu)) )

  Call Chop(Dmue)
  Call Chop(DB)

  Call ParametersToG2(Dgauge, DYe, DYnu, DYd, DYu, DMhlf, DAe, DAnu, DAd, DAu &
                   &, DMe, DMl, DMr, DMd, DMq, DMu, DMh, Dmue, DB, f)

  Iname = Iname - 1

 End Subroutine rge267


 Subroutine rge277(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.01.06: including neutrino dim 5 operator
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2
  Real(dp) :: gauge(3), gauge2(3), TraceY(4), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2)           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3), hc(4), sumI

  Complex(dp) :: Mhlf(3),DMhlf(3)

  Complex(dp) :: Ae(3,3), Ad(3,3), Au(3,3), aAe(3,3), aAd(3,3), aAu(3,3)   &
     &  , DAe(3,3), DAd(3,3), DAu(3,3), aAdAd(3,3), aAeAe(3,3), aAuAu(3,3) &
     &  , aYdAd(3,3), aYeAe(3,3), aYuAu(3,3), TraceaYA(4), betaAd1(3,3)    &
     &  , betaAd2(3,3), betaAe1(3,3), betaAe2(3,3), betaAu1(3,3)           &
     &  , betaAu2(3,3)
  Real(dp) :: TraceA(4)
  Complex(dp) :: aYdYdaYdAd(3,3), aYdAdaYdYd(3,3), TraceAY2(5)               &
     &  , aYeYeaYeAe(3,3), aYeAeaYeYe(3,3), aYuYuaYuAu(3,3), aYuAuaYuYu(3,3) &
     &  , aYuYuaYdAd(3,3), aYuAuaYdYd(3,3), aYdYdaYuAu(3,3), aYdAdaYuYu(3,3)
   
  Complex(dp) :: Me(3,3), Ml(3,3), Md(3,3), Mq(3,3), Mu(3,3), DMe(3,3)        &
     & , DMl(3,3), DMd(3,3), DMq(3,3), DMu(3,3), YdaYd(3,3), YeaYe(3,3)       &
     & , YuaYu(3,3), MdYdaYd(3,3), MeYeaYe(3,3), MuYuaYu(3,3), YdaYdMd(3,3)   &
     & , YeaYeMe(3,3), YuaYuMu(3,3), YdMqaYd(3,3), YeMlaYe(3,3)               &
     & ,  YuMqaYu(3,3), AdaAd(3,3), AeaAe(3,3), AuaAu(3,3), betaMd1(3,3)      &
     & , betaMd2(3,3), betaMe1(3,3), betaMe2(3,3), betaMl1(3,3), betaMl2(3,3) &
     & , betaMq1(3,3), betaMq2(3,3), betaMu1(3,3), betaMu2(3,3), MqaYdYd(3,3) &
     & , MqaYuYu(3,3), aYdYdMq(3,3), aYuYuMq(3,3), aYeYeMl(3,3), MlaYeYe(3,3) &
     & , aYeMeYe(3,3), aYdMdYd(3,3), aYuMuYu(3,3)                             &
     & , YdaYdYdaYd(3,3), YeaYeYeaYe(3,3), YuaYuYuaYu(3,3), MeYeaYeYeaYe(3,3) &
     & , YeaYeYeaYeMe(3,3), YeaYeMeYeaYe(3,3), YeMlaYeYeaYe(3,3), AeaYe(3,3)  &
     & , YeaYeYeMlaYe(3,3), AeaAeYeaYe(3,3), YeaYeAeaAe(3,3), YeaAe(3,3)      &
     & , AeaYeYeaAe(3,3), YeaAeAeaYe(3,3), Tr3aAdYdaAeYe

  Complex(dp) :: Tr3aYdAdaYeAe, AdaYd(3,3), YdaAd(3,3), MlaYeYeaYeYe(3,3)     &
     & , aYeYeaYeYeMl(3,3), aYeYeMlaYeYe(3,3), aYeYeaYeMeYe(3,3)              &
     & , aYeMeYeaYeYe(3,3), aAdYd(3,3),aAeYe(3,3), aAeAeaYeYe(3,3)            &
     & , aYeYeaAeAe(3,3), aAeYeaYeAe(3,3), aYeAeaAeYe(3,3), MdYdaYdYdaYd(3,3) &
     & , YdaYdYdaYdMd(3,3), YdMqaYdYdaYd(3,3),YdaYdMdYdaYd(3,3)               &
     & , YdaYdYdMqaYd(3,3), AdaAdYdaYd(3,3), YdaYDAdaAd(3,3), AdaYdYdaAd(3,3) &
     & , YdaAdAdaYd(3,3)
  Complex(dp) :: MdYdaYuYuaYd(3,3), YdaYuYuaYdMd(3,3), YdMqaYuYuaYd(3,3)      &
     & , YdaYuYuMqaYd(3,3), YdaYuMuYuaYd(3,3), AdaAuYuaYd(3,3)                &
     & , YdaYuAuaAd(3,3), AdaYuYuaAd(3,3), YdaAuAuaYd(3,3), YdaYuYuaYd(3,3)   &
     & , Tr3aYuAu, Tr3aAuYu, YuaAu(3,3)                                       &
     & , MqaYdYdaYdYd(3,3), aYdYdaYdYDMq(3,3), aYdMdYdaYdYd(3,3)              &
     & , aYdYdMqaYdYd(3,3), aYdYdaYdMdYd(3,3), aAdAdaYdYd(3,3)                &
     & , aYdYDaAdAd(3,3), aAdYdaYdAd(3,3), aYdAdaAdYd(3,3), MqaYuYuaYuYu(3,3) &
     & , aYuYuaYuYUMq(3,3), aYuMuYuaYuYu(3,3), aYuYuMqaYuYu(3,3)              &
     & , aYuYuaYuMuYu(3,3), aAuAuaYuYu(3,3), aYuYUaAuAu(3,3), aAuYuaYuAu(3,3) &
     & , aYuAuaAuYu(3,3), aAuYu(3,3), AuaYu(3,3), YuaYdYdaYu(3,3)             &
     & , AuaYdYdaAu(3,3), YuaAdAdaYu(3,3), AuaAdYdaYu(3,3), YuaYdAdaAu(3,3)   &
     & , YuMqaYuYuaYu(3,3), YuaYuYuMqaYu(3,3), MuYuaYuYuaYu(3,3)              &
     & , YuaYuYuaYuMu(3,3), YuaYuMuYuaYu(3,3), AuaAuYuaYu(3,3)                &
     & , YuaYuAuaAu(3,3), AuaYuYuaAu(3,3), YuaAuAuaYu(3,3), MuYuaYdYdaYu(3,3) &
     & , YuaYdYdaYuMu(3,3), YuMqaYdYdaYu(3,3), YuaYdYdMqaYu(3,3)              &
     & , YuaYdMdYdaYu(3,3)

  Real(dp) :: S1, S2, sig(3), Tr3aYdYdaYeYe, Tr3aAdAdaAeAe, AbsGM2(3)         &
     & , Tr3MqaYdYd3aYDMdYd , Tr3MqaYuYu3aYuMu, Tr3aAuAu

  Real(dp) :: Mh(2), DMh(2), TraceMH1(3), TraceMH2(2), betaMH11        &
     &  , betaMH12, betaMH21, betaMH22, q

  Complex(dp) :: mue, B, Dmue, DB, TraceMue(2), TraceB(2), betaMue1, betaMue2 &
     & , betaB1, betaB2

  Complex(dp) :: Mnu(3,3), DMnu(3,3), betamnu1(3,3)

  Complex(dp), Dimension(3,3) :: YT, aYT, aYTYT, betaYT1, sumT1, AT, aAT     &
      & , aATAT, betaAT1, aYTAT, MlaYTYT, aYTMlYT, aYTYTMl , DYT, DAT
  Complex(dp) :: g2Mi(3)
  Complex(dp) :: lam1, lam2, Dlam1, Dlam2, betalam11, betalam21, lam1Alam1 &
      & , lam2Alam2, Alam1, Alam2, betaAlam11, betaAlam21, DAlam1, DAlam2
  Real(dp) :: MT(2), lam12, lam22, betaMT2(2), Alam12, Alam22, b_1a(3)     &
      & , DMT2(2), b_2a(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge277'

  OnlyDiagonal = .Not.GenerationMixing
  q = t

  Call GToParameters4(gy, gauge, Ye, YT, Yd, Yu, lam1, lam2, Mhlf, Ae, AT, Ad &
                & , Au, Alam1, Alam2, Me, Ml, Md, Mq, Mu, Mh, mT, mue, B, Mnu)

  gauge2 = gauge**2
  AbsGM2 = gauge2 * Abs( Mhlf )**2
!-----------------
! beta functions
!-----------------
  lam12 = Abs(lam1)**2
  lam22 = Abs(lam2)**2

  Call Adjungate(Yd,aYd)
  Call Adjungate(YT,aYT)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul(aYe,Ye)
  aYTYT = Matmul(aYT,YT)
  aYdYd = Matmul(aYd,Yd)
  aYuYu = Matmul(aYu,Yu)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYTYT),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = 3._dp * (TraceY(3) + lam12) + TraceY(1)    &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * (aYeYe + aYTYT)
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
  sumT1 = aYeYe + 6._dp * aYTYT
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaYT1 = Matmul(YT,sumT1)  &
        & + Matmul(Transpose(aYeYe),YT)

  diagonal(3,1) = 3._dp * (TraceY(3)  + lam12) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)

  diagonal(4,1) = 3._dp * (TraceY(4) + lam22)                &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  betalam11 = lam1 * (7._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
            &        + 6._dp * TraceY(3)                            &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  betalam21 = lam2 * (7._dp * lam22  + 6._dp * TraceY(4)     &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
   aYdYdaYuYu = Matmul(aYdYd,aYuYu)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ae,aAe)
  Call Adjungate(AT,aAT)
  Call Adjungate(Ad,aAd)
  Call Adjungate(Au,aAu)

  aAdAd = Matmul(aAd,Ad)
  aATAT = Matmul(aAT,AT)
  aAeAe = Matmul(aAe,Ae)
  aAuAu = Matmul(aAu,Au)

  TraceA(1) = Real( cTrace(aAeAe),dp )
  TraceA(2) = Real( cTrace(aATAT),dp )
  TraceA(3) = Real( cTrace(aAdAd),dp )
  TraceA(4) = Real( cTrace(aAuAu),dp )

  aYdAd = Matmul(aYd,Ad)
  aYTAT = Matmul(aYT,AT)
  aYeAe = Matmul(aYe,Ae)
  aYuAu = Matmul(aYu,Au)

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYTAT) 
  TraceaYA(3) = cTrace(aYdAd) 
  TraceaYA(4) = cTrace(aYuAu) 
  lam1Alam1 = Conjg(lam1) * Alam1
  lam2Alam2 = Conjg(lam2) * Alam2
  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = Matmul(Ae,sume1)
  
  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1)  &
                &         + 3._dp * lam1Alam1                  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    ) 
  sume1 = 4._dp * aYeAe + 6._dp * aYTAT
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + Matmul(Ye,sume1)

  !--------------
  ! A_T
  !--------------
   If (decoupling_heavy_states) Then
    betaAT1 = 0._dp
   Else  
    diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
    sumT1 = aYeYe + 9._dp * aYTYT
    betaAT1 = Matmul(AT,sumT1)
    betaAT1 = Transpose(betaAT1)
    Do i1=1,3
     sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
    End Do

    betaAT1 = betaAT1 + Matmul(sumT1,AT)
  
    diagonal(2,1) = 2._dp * ( TraceaYA(2) + lam1Alam1    &
      &         + 1.8_dp * g2Mi(1) + 7._dp * g2Mi(2) )
    betaAT1 = betaAT1 + diagonal(2,1) * YT
   End If
  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = Matmul(Ad,sumd1)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1) &
                &         + 3._dp * lam1Alam1                 &
                &         - c1_1(2,1) * g2Mi(1) - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do
  betaAd1 = betaAd1 + Matmul(Yd,sumd1)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = Matmul(Au,sumu1)
  
  diagonal(4,1) = 2._dp * ( 3._dp * TraceaYA(4)                       &
                &         + 3._dp * lam2Alam2 - c1_1(3,1) * g2Mi(1)   &
                &         - c1_1(3,2) * g2Mi(2) - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do
  betaAu1 = betaAu1 + Matmul(Yu,sumu1)

  !--------------
  ! A_1
  !--------------
  betaAlam11 = Alam1 * (21._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
          &           + 6._dp * TraceY(3)                              &
          &           - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) )       &
          & + 2._dp * lam1 * ( TraceaYA(2) + 2._dp * TraceaYA(1)       &
          &                  + 6._dp * TraceaYA(3)                     &
          &                  + 1.8_dp * g2Mi(1) + 7._dp * g2Mi(2) ) 

  !--------------
  ! A_2
  !--------------
  betaAlam21 = Alam2 * (21._dp * lam22 + 6._dp * TraceY(4)          &
          &           - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) )    &
          & + 2._dp * lam2 * ( 6._dp * TraceaYA(3)                  &
          &                  + 1.8_dp * g2Mi(1) + 7._dp * g2Mi(2) ) 

  If (TwoLoopRGE) Then
   aYdYdaYdAd = Matmul(aYdYd,aYdAd)
   aYdAdaYdYd = Matmul(aYdAd,aYdYd)
   aYeYeaYeAe = Matmul(aYeYe,aYeAe)
   aYeAeaYeYe = Matmul(aYeAe,aYeYe)
   aYuYuaYuAu = Matmul(aYuYu,aYuAu)
   aYuAuaYuYu = Matmul(aYuAu,aYuYu)
   aYuAuaYdYd = Matmul(aYuAu,aYdYd)
   aYuYuaYdAd = Matmul(aYuYu,aYdAd)
   aYdAdaYuYu = Matmul(aYdAd,aYuYu)
   aYdYdaYuAu = Matmul(aYdYd,aYuAu)
   TraceAY2(1) = cTrace(aYeYeaYeAe)
   TraceAY2(2) = cTrace(aYdYdaYdAd)
   TraceAY2(3) = cTrace(aYuYuaYuAu)
   TraceAY2(4) = cTrace(aYuYuaYdAd)
   TraceAY2(5) = cTrace(aYdYdaYuAu)

  !--------------
  ! A_e
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = Matmul(Ae,sume2)
    
   diagonal(1,2) = -6._dp * ( 6._dp * TraceAY2(2) + TraceAY2(4)       &
     &                      + TraceAY2(5) + 2._dp * TraceAY2(1)  )    &
     &  + ( 32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)   &
     &  - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)             &
     &  + 2.4_dp * gauge2(1) * TraceaYA(1) - 2.4_dp * g2Mi(1) * TraceY(1) &
     &  - ( 30._dp * g2Mi(2)                              &
     &    + 3.6_dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)      &
     &  - 54._dp * gauge2(1)**2 * Mhlf(1)
   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)     &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1) &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2)
   sume2 = - 6._dp * aYeYeaYeAe - 8._dp * aYeAeaYeYe  &
         & - hd(1) * aYeAe - hc(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
   betaAe2 = betaAe2 + Matmul(Ye,sume2)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = Matmul(Ad,sumd2)
    
   diagonal(3,2) = diagonal(1,2)                                   &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                 &
     &              - gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp      &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)    &
     &  + 1.6_dp * gauge2(1) * gauge2(2) * (Mhlf(1)+Mhlf(2))       & 
     &  + 1.856e3_dp * gauge2(1)**2 * Mhlf(1) / 4.5e1_dp

   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)  &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hd(2) = 6._dp * TraceY(4) - 1.6_dp * gauge2(1) 
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1)                 &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2) + 1.6_dp * gauge2(1) *  Mhlf(1)
   hc(2) = 6._dp * TraceaYA(4) + 1.6_dp * gauge2(1) *  Mhlf(1)
   sumd2 = - 6._dp * aYdYdaYdAd - 8._dp * aYdAdaYdYd                  &
       &   - 4._dp * ( aYuAuaYuYu + aYuYuaYuAu + aYuAuaYdYd )         &
       &   - 2._dp * aYuYuaYdAd - hd(1) * aYdAd - hc(1) * aYdYd       &
       &  - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(3,2)
   End Do
   betaAd2 = betaAd2 + Matmul(Yd,sumd2)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(4) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = Matmul(Au,sumu2)
    
   diagonal(4,2) =  -6._dp * ( 6._dp * TraceAY2(3) + TraceAY2(4)        &
     &                        + TraceAY2(5)  )                          &
     &  + ( 32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)     &
     &  - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)           &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                      &
     &              - 3.4_dp * gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp  &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)         &
     &  - ( 30._dp * g2Mi(2)                                &
     &    + 2._dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)         &
     &  - 5486._dp * gauge2(1)**2 * Mhlf(1) / 225._dp
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1) - 0.8_dp * gauge2(1)
   hc(1) = 6._dp * TraceaYA(3) + 2._dp * TraceaYA(1)   &
       &  + 0.8_dp * gauge2(1) *  Mhlf(1)
   hd(2) = 12._dp * TraceY(4) - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(2) = 18._dp * TraceaYA(4) + 1.2e1_dp * g2Mi(2) + 0.8_dp * g2Mi(1)
   sumu2 = - 6._dp * aYuYuaYuAu - 8._dp * aYuAuaYuYu                     &
       &   - 4._dp * ( aYdAdaYdYd + aYdYdaYdAd + aYdAdaYuYu )            &
       &   - 2._dp * aYdYdaYuAu - hd(1) * aYdAd - hc(1) * aYdYd          &
       &   - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(4,2)
   End Do
   betaAu2 = betaAu2 + Matmul(Yu,sumu2)

  End If 
!----------------------------------------------
! beta functions for Sfermion mass parameters
!----------------------------------------------
   S1 = mH(2) - mH(1)  + 3._dp * (mT(1) - mT(2) )
   Do i1=1,3
    S1 = S1 + Me(i1,i1) - Ml(i1,i1) &
       &    + Md(i1,i1) + Mq(i1,i1) - 2._dp * Mu(i1,i1)
   End Do
   S1 = S1 * gauge2(1)

   YdaYd = Matmul(Yd,aYd)
   YeaYe = Matmul(Ye,aYe)
   YuaYu = Matmul(Yu,aYu)

   MeYeaYe = Matmul(Me,YeaYe)
   MlaYeYe = Matmul(Ml,aYeYe)
   MlaYTYT = Matmul(Ml,aYTYT)

   MdYdaYd = Matmul(Md,YdaYd)
   MqaYdYd = Matmul(Mq,aYdYd)
   MqaYuYu = Matmul(Mq,aYuYu)
   MuYuaYu = Matmul(Mu,YuaYu)

   YeaYeMe = Matmul(YeaYe,Me)
   aYeYeMl = Matmul(aYeYe,Ml)
   aYTYTMl = Matmul(aYTYT,Ml)

   YdaYdMd = Matmul(YdaYd,Md)
   aYdYdMq = Matmul(aYdYd,Mq)
   aYuYuMq = Matmul(aYuYu,Mq)
   YuaYuMu = Matmul(YuaYu,Mu)

   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   aYTMlYT = MatMul3(aYT,Transpose(Ml),YT,OnlyDiagonal)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)

   AeaAe = Matmul(Ae,aAe)
   AdaAd = Matmul(Ad,aAd)
   AuaAu = Matmul(Au,aAu)
   Alam12 = Abs(Alam1)**2
   Alam22 = Abs(Alam2)**2

   diagonal(1,1) = - 4.8_dp * AbsGM2(1) + 1.2_dp * S1
   betaMe1 = 2._dp * (MeYeaYe + YeaYeMe)             &
         & + 4._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe )
   Do i1=1,3
    betaMe1(i1,i1) = betaMe1(i1,i1) + diagonal(1,1)
   End Do

   diagonal(3,1) = - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1
   betaMl1 = MlaYeYe + aYeYeMl + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe )  &
         & + 3._dp * ( MlaYTYT + aYTYTMl )                                    &
         & + 6._dp * ( aYTMlYT + aATAT + MT(1) * aYTYT)
   Do i1=1,3
    betaMl1(i1,i1) = betaMl1(i1,i1) + diagonal(3,1)
   End Do

   diagonal(4,1) = - ( 32._dp * AbsGM2(3) + 1.6_dp * AbsGM2(1) ) / 3._dp &
               & + 0.4_dp * S1
   betaMd1 = 2._dp * (MdYdaYd + YdaYdMd)             &
         & + 4._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd )
   Do i1=1,3
    betaMd1(i1,i1) = betaMd1(i1,i1) + diagonal(4,1)
   End Do

   diagonal(5,1) = - ( 32._dp * AbsGM2(3) + 0.4_dp * AbsGM2(1) ) / 3._dp &
               & - 6._dp * AbsGM2(2) + 0.2_dp * S1
   betaMq1 = MqaYuYu + aYuYuMq + MqaYdYd + aYdYdMq             &
         & + 2._dp * ( mH(2) * aYuYu + mH(1) * aYdYd + aYuMuYu     &
         &           + aYdMdYd + aAuAu + aAdAd )
   Do i1=1,3
    betaMq1(i1,i1) = betaMq1(i1,i1) + diagonal(5,1)
   End Do

   diagonal(6,1) = - ( 32._dp * AbsGM2(3) + 6.4_dp * AbsGM2(1) ) / 3._dp &
               &   - 0.8_dp * S1
   betaMu1 = 2._dp * (MuYuaYu + YuaYuMu)             &
         & + 4._dp * ( mH(2) * YuaYu + YuMqaYu + AuaAu )
   Do i1=1,3
    betaMu1(i1,i1) = betaMu1(i1,i1) + diagonal(6,1)
   End Do

   If (TwoLoopRGE) Then
    YdaYdYdaYd = MatSquare(YdaYd,OnlyDiagonal)
    YeaYeYeaYe = MatSquare(YeaYe,OnlyDiagonal)
    YuaYuYuaYu = MatSquare(YuaYu,OnlyDiagonal)

    AdaYd = Matmul(Ad,aYd)
    AeaYe = Matmul(Ae,aYe)
    AuaYu = Matmul(Au,aYu)

    aAdYd = Matmul(aAd,Yd)
    aAeYe = Matmul(aAe,Ye)
    aAuYu = Matmul(aAu,Yu)

    YdaAd = Matmul(Yd,aAd)
    YeaAe = Matmul(Ye,aAe)
    YuaAu = Matmul(Yu,aAu)

    YdaYuYuaYd = MatMul3(Yd,aYuYu,aYd,OnlyDiagonal)
    AdaYuYuaAd = MatMul3(Ad,aYuYu,aAd,OnlyDiagonal)
    YdaAuAuaYd = MatMul3(Yd,aAuAu,aYd,OnlyDiagonal)
    AdaAuYuaYd = MatMul4(Ad,aAu,Yu,aYd,OnlyDiagonal)
    YdaYuAuaAd = MatMul3(Yd,aYuAu,aAd,OnlyDiagonal)

    YuaYdYdaYu = MatMul3(Yu,aYdYd,aYu,OnlyDiagonal)
    AuaYdYdaAu = MatMul3(Au,aYdYd,aAu,OnlyDiagonal)
    YuaAdAdaYu = MatMul3(Yu,aAdAd,aYu,OnlyDiagonal)
    AuaAdYdaYu = MatMul4(Au,aAd,Yd,aYu,OnlyDiagonal)
    YuaYdAdaAu = MatMul3(Yu,aYdAd,aAu,OnlyDiagonal)

    MdYdaYuYuaYd = Matmul(Md,YdaYuYuaYd)
    Call Adjungate(MdYdaYuYuaYd, YdaYuYuaYdMd)
    YdMqaYuYuaYd = MatMul3(Yd,MqaYuYu,aYd,OnlyDiagonal)
    Call Adjungate(YdMqaYuYuaYd, YdaYuYuMqaYd)
    YdaYuMuYuaYd = MatMul3(Yd,aYuMuYu,aYd,OnlyDiagonal)

    MuYuaYdYdaYu = Matmul(Mu,YuaYdYdaYu)
    Call Adjungate(MuYuaYdYdaYu, YuaYdYdaYuMu)
    YuMqaYdYdaYu = MatMul3(Yu,MqaYdYd,aYu,OnlyDiagonal)
    Call Adjungate(YuMqaYdYdaYu, YuaYdYdMqaYu)
    YuaYdMdYdaYu = MatMul3(Yu,aYdMdYd,aYu,OnlyDiagonal)

    MeYeaYeYeaYe = Matmul(MeYeaYe,YeaYe)
    Call Adjungate(MeYeaYeYeaYe,YeaYeYeaYeMe)
    aYeMeYeaYeYe = Matmul(aYeMeYe,aYeYe)
    Call Adjungate(aYeMeYeaYeYe,aYeYeaYeMeYe)
    YeaYeMeYeaYe = Matmul(YeaYeMe,YeaYe)

    MlaYeYeaYeYe = Matmul(MlaYeYe,aYeYe)
    Call Adjungate(MlaYeYeaYeYe, aYeYeaYeYeMl)
    YeMlaYeYeaYe = Matmul(YeMlaYe,YeaYe)
    Call Adjungate(YeMlaYeYeaYe, YeaYeYeMlaYe)
    aYeYeMlaYeYe = Matmul(aYeYeMl,aYeYe)

    MdYdaYdYdaYd = Matmul(MdYdaYd,YdaYd)
    Call Adjungate(MdYdaYdYdaYd, YdaYdYdaYdMd)
    aYdMdYdaYdYd = Matmul(aYdMdYd,aYdYd)
    Call Adjungate(aYdMdYdaYdYd, aYdYdaYdMdYd)
    YdaYdMdYdaYd = Matmul(YdaYdMd,YdaYd)

    MqaYdYdaYdYd = Matmul(MqaYdYd,aYdYd)
    Call Adjungate(MqaYdYdaYdYd, aYdYdaYdYdMq)
    YdMqaYdYdaYd = Matmul(YdMqaYd,YdaYd)
    Call Adjungate(YdMqaYdYdaYd, YdaYdYdMqaYd)
    aYdYdMqaYdYd = Matmul(aYdYdMq,aYdYd)

    MqaYuYuaYuYu = Matmul(MqaYuYu,aYuYu)
    Call Adjungate(MqaYuYuaYuYu, aYuYuaYuYuMq)
    YuMqaYuYuaYu = Matmul(YuMqaYu,YuaYu)
    Call Adjungate(YuMqaYuYuaYu, YuaYuYuMqaYu)
    aYuYuMqaYuYu = Matmul(aYuYuMq,aYuYu)

    MuYuaYuYuaYu = Matmul(MuYuaYu,YuaYu)
    Call Adjungate(MuYuaYuYuaYu, YuaYuYuaYuMu)
    aYuMuYuaYuYu = Matmul(aYuMuYu,aYuYu)
    Call Adjungate(aYuMuYuaYuYu, aYuYuaYuMuYu)
    YuaYuMuYuaYu = Matmul(YuaYuMu,YuaYu)

    AdaAdYdaYd = Matmul(AdaAd,YdaYd)
    Call Adjungate(AdaAdYdaYd, YdaYdAdaAd)
    AdaYdYdaAd = Matmul(AdaYd,YdaAd)
    YdaAdAdaYd = Matmul(YdaAd,AdaYd)

    aAdAdaYdYd = Matmul(aAdAd,aYdYd)
    Call Adjungate(aAdAdaYdYd, aYdYdaAdAd)
    aAdYdaYdAd = Matmul(aAdYd,aYdAd)
    aYdAdaAdYd = Matmul(aYdAd,aAdYd)

    AeaAeYeaYe = Matmul(AeaAe,YeaYe)
    Call Adjungate(AeaAeYeaYe, YeaYeAeaAe)
    AeaYeYeaAe = Matmul(AeaYe,YeaAe)
    YeaAeAeaYe = Matmul(YeaAe,AeaYe)

    aAeAeaYeYe = Matmul(aAeAe,aYeYe)
    Call Adjungate(aAeAeaYeYe, aYeYeaAeAe)
    aAeYeaYeAe = Matmul(aAeYe,aYeAe)
    aYeAeaAeYe = Matmul(aYeAe,aAeYe)

    AuaAuYuaYu = Matmul(AuaAu,YuaYu)
    Call Adjungate(AuaAuYuaYu, YuaYuAuaAu)
    AuaYuYuaAu = Matmul(AuaYu,YuaAu)
    YuaAuAuaYu = Matmul(YuaAu,AuaYu)

    aAuAuaYuYu = Matmul(aAuAu,aYuYu)
    Call Adjungate(aAuAuaYuYu, aYuYuaAuAu)
    aAuYuaYuAu = Matmul(aAuYu,aYuAu)
    aYuAuaAuYu = Matmul(aYuAu,aAuYu)

    S2 = (1.5_dp * gauge2(2) + 0.3_dp * gauge2(1) )            &
     &      * (MH(2) - MH(1) - Real(cTrace(ML),dp) )              &
     & + ( (8._dp * gauge2(3) + 0.1_dp*gauge2(1)) / 3._dp      &
     &   + 1.5_dp * gauge2(2) ) * Real( cTrace(Mq),dp   )         &
     & - (16._dp * gauge2(3) + 3.2_dp*gauge2(1) )              &
     &    * Real(cTrace(Mu),dp) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md),dp) / 3._dp                           &
     & + 1.2_dp*gauge2(1) * Real(cTrace(Me),dp)                   &
     & - 3._dp * (MH(2)*TraceY(4) - MH(1) * TraceY(3) )        &
     & + mH(1) * TraceY(1)
    
    Do i1=1,3
     S2 = S2 - YuMqaYu(i1,i1) + 4._dp * aYuMuYu(i1,i1)   &
        &    - YdMqaYd(i1,i1) - 2._dp * aYdMdYd(i1,i1)   &
        &    + YeMlaYe(i1,i1) - 2._dp * aYeMeYe(i1,i1)
    End Do

    sig(1) = 3._dp * (MH(1) + MH(2) + Real(cTrace(Ml),dp) )          &
         & + Real(cTrace(Mq),dp) + 8._dp * Real(cTrace(Mu),dp)          &
         & + 2._dp * Real(cTrace(Md),dp) + 6._dp * Real(cTrace(Me),dp) 
    sig(1) = 0.2_dp * gauge2(1) * sig(1)
    sig(2) = gauge2(2) * ( MH(1) + MH(2) + Real( cTrace(Ml),dp )    &
           &             + 3._dp * Real( cTrace(Mq),dp ) )
    sig(3) = gauge2(3) * ( 2._dp * Real( cTrace(Mq),dp )            &
           &             + Real(cTrace(Mu),dp) + Real(cTrace(Md),dp) )

    Tr3aYdYdaYeYe = 3._dp * TraceY(3) + TraceY(1)
    Tr3MqaYdYd3aYDMdYd = Real(cTrace(MlaYeYe),dp) + Real(cTrace(aYeMeYe),dp)  &
            & + 3._dp * ( Real(cTrace(MqaYdYd),dp) + Real(cTrace(aYdMdYd),dp) )
    Tr3aAdAdaAeAe = 3._dp * TraceA(3) + TraceA(1)
    Tr3aYdAdaYeAe = 3._dp * TraceaYA(3) + TraceaYA(1)
    Tr3aAdYdaAeYe = Conjg( Tr3aYdAdaYeAe )
    Tr3MqaYuYu3aYuMu = 3._dp * (Real(cTrace(MqaYuYu),dp) &
                               + Real(cTrace(aYuMuYu),dp) )
    Tr3aAuAu = 3._dp * TraceA(4)
    Tr3aYuAu = 3._dp * TraceaYA(4)
    Tr3aAuYu = Conjg( Tr3aYuAu )

    diagonal(1,2) = 2.4_dp * gauge2(1) * (S2 + sig(1) )   &
                & + 112.32_dp * gauge2(1) * AbsGM2(1)
    hd(1) = 6._dp*gauge2(2)-1.2_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 4.8_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 2.4_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMe2 = -2._dp * ( MeYeaYeYeaYe + YeaYeYeaYeMe )                       &
      & - 4._dp * ( YeMlaYeYeaYe + YeaYeMeYeaYe + YeaYeYeMlaYe )             &
      & - 8._dp * mH(1) * YeaYeYeaYe                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MeYeaYe + YeaYeMe + 4._dp * MH(1) * YeaYe &
      &                           + 2._dp * (YeMlaYe + AeaAE)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YeaYe + AeaAeYeaYe + YeaYEAeaAe     &
      &           + AeaYeYeaAe + YeaAeAeaYe + Tr3aAdAdaAeAe * YeaYe          &
      &           + Tr3aAdYdaAeYe * AeaYe + Tr3aYdAdaYeAe * YeaAe )          &
      & + hd(1) * ( MeYeaYe + YeaYeMe                                        &
      &           + 2._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe ) )            &
      & + hd(2) * YeaYe + hc(1) * AeaYe + hc(2) * YeaAe
    Do i1=1,3
     betaMe2(i1,i1) = betaMe2(i1,i1) + diagonal(1,2)
    End Do

    diagonal(3,2) = gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)        &
     &     + 3._dp * gauge2(2) * sig(2)                                 &
     &     + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)      &
     &                   + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2       &
     &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )     &
     &     + 24.84_dp * gauge2(1) * AbsGM2(1)
    betaMl2 = -2._dp * ( MlaYeYeaYeYe + aYeYeaYeYeMl )                   &
      & - 4._dp * ( aYeMeYeaYeYe + aYeYeMlaYeYe + aYeYeaYeMeYe )         &
      & - 8._dp * mH(1) * aYeYeaYeYe                                     &
      & - Tr3aYdYdaYeYe * ( MlaYeYe + aYeYeML + 4._dp * MH(1) * aYeYe    &
      &                   + 2._dp * (aYeMeYe + aAeAE)  )                 &
      & - 4._dp * (aAeAeaYeYe + aYeYeaAeAe + aAeYeaYeAe + aYeAeaAeYe )   &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYeYe +  Tr3aAdAdaAeAe * aYeYe  &
      &           + Tr3aAdYdaAeYe * aYeAe + Tr3aYdAdaYeAe * aAeYe     )  &
      & + 1.2_dp*gauge2(1) * ( MlaYeYe + aYeYeMl                         &
      &                      + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe &
      &                                - Mhlf(1) * aAeYe                 &
      &                                - Conjg(Mhlf(1)) * aYeAe ) )      &
      & + 4.8_dp * AbsGM2(1) * aYeYe
    Do i1=1,3
     betaMl2(i1,i1) = betaMl2(i1,i1) + diagonal(3,2)
    End Do

    diagonal(4,2) = 0.8_dp * gauge2(1) * (S2 + sig(1)/3._dp )             &
      &   + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp   &
      &   + 8.08e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                    &
      &   + 1.28e2_dp * ( gauge2(1) * AbsGM2(3)                           &
      &                 + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)     &
      &                   * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)+0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) + 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) - 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMd2 = -2._dp * ( MdYdaYdYdaYd + YdaYdYdaYdMd )                       &
      & - 4._dp * ( YdMqaYdYdaYd + YdaYdMdYdaYd + YdaYdYdMqaYd )             &
      & - 8._dp * mH(1) * YdaYdYdaYd                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MdYdaYd + YdaYdMd + 4._dp * MH(1) * YdaYd &
      &                           + 2._dp * (YdMqaYd + AdaAD)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YdaYd + AdaAdYdaYd + YdaYDAdaAd     &
      &           + AdaYdYdaAd + YdaAdAdaYd + Tr3aAdAdaAeAe * YdaYd          &
      &           + Tr3aAdYdaAeYe * AdaYd + Tr3aYdAdaYeAe * YdaAd )          &
      & + hd(1) * ( MdYdaYd + YdaYdMd                                        &
      &           + 2._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd ) )            &
      & + hd(2) * YdaYd + hc(1) * AdaYd + hc(2) * YdaAd                      &
      & - 2._dp * (MdYdaYuYuaYd + YdaYuYuaYdMd )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YdaYuYuaYd + YdMqaYuYuaYd+ YdaYuYuMqaYd  &
      &           + YdaYuMuYuaYd + AdaAuYuaYd + YdaYuAuaAd                   &
      &           + AdaYuYuaAd + YdaAuAuaYd )
    Do i1=1,3
     betaMd2(i1,i1) = betaMd2(i1,i1) + diagonal(4,2)
    End Do

    diagonal(5,2) = 0.2_dp * gauge2(1) * (2._dp * S2 + sig(1)/3._dp )     &
      &  + 3._dp * gauge2(2) * sig(2)                                     &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp    &
      &  + 1.99e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                     &
      &  + 32._dp * ( gauge2(1) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)         &
      &              * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp      &
      &  + 33._dp * gauge2(2) * AbsGM2(2)                                 &
      &  + 32._dp * ( gauge2(2) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(2)  + gauge2(3)*gauge2(2)        &
      &              * Real( Mhlf(2) * Conjg(Mhlf(3)),dp ) )                 &
      &  + 0.4_dp * ( gauge2(2) * AbsGM2(1)                               &
      &             + gauge2(1) * AbsGM2(2) + gauge2(1)*gauge2(2)         &
      &              * Real( Mhlf(2) * Conjg(Mhlf(1)),dp ) )
    hd(1) = 1.6_dp * AbsGM2(1)
    hc(1) = - 0.8_dp * g2Mi(1)
    hc(2) = Conjg( hc(1) )
    hd(2) = 3.2_dp * AbsGM2(1)
    hc(3) = - 1.6_dp * g2Mi(1)
    hc(4) = Conjg( hc(3) )
    betaMq2 = -2._dp * ( MqaYdYdaYdYd + aYdYdaYdYDMq )                       &
      & - 4._dp * ( aYdMdYdaYdYd + aYdYdMqaYdYd + aYdYdaYdMdYd )             &
      & - 8._dp * mH(1) * aYdYdaYdYd                                         &
      & - Tr3aYdYdaYeYe * ( MqaYdYd + aYdYdMq + 4._dp * MH(1) * aYdYd        &
      &                   + 2._dp * (aYdMdYd + aAdAd)  )                     &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYdYd                               &
      &           + 2._dp * ( aAdAdaYdYd + aYdYDaAdAd + aAdYdaYdAd           &
      &                     + aYdAdaAdYd )                                   &
      &           + Tr3aAdAdaAeAe * aYdYd + Tr3aAdYdaAeYe * aYdAd            &
      &           + Tr3aYdAdaYeAe * aAdYd )                                  &
      & + 0.4_dp*gauge2(1) * ( MqaYdYd + aYdYdMq                             &
      &                      + 2._dp * ( mH(1) * aYdYd + aYdMdYd + aAdAd ) ) &
      & + hd(1) * aYdYd + hc(1) * aAdYd + hc(2) * aYdAd                      &
      & - 2._dp * ( MqaYuYuaYuYu + aYuYuaYuYuMq )                            &
      & - 4._dp * ( aYuMuYuaYuYu + aYuYuMqaYuYu + aYuYuaYuMuYu )             &
      & - 8._dp * MH(2) * aYuYuaYuYu                                         &
      & - 3._dp * TraceY(4) * ( MqaYuYu + aYuYuMq + 4._dp * MH(2) * aYuYu    &
      &                       + 2._dp * (aYuMuYu + aAuAu) )                  &
      & - 2._dp * ( Tr3MqaYuYu3aYuMu * aYuYu                                 &
      &           + 2._dp * ( aAuAuaYuYu + aYuYuaAuAu                        &
      &                     + aAuYuaYuAu + aYuAuaAuYu )                      &
      &           + Tr3aAuAu * aYuYu + Tr3aAuYu * aYuAu + Tr3aYuAu * aAuYu ) &
      & + 0.8_dp*gauge2(1) * ( MqaYuYu + aYuYuMq                             &
      &                      + 2._dp * ( MH(2) * aYuYu + aYuMuYu + aAuAu ) ) &
      & + hd(2) * aYuYu + hc(3) * aAuYu + hc(4) * aYuAu
    Do i1=1,3
     betaMq2(i1,i1) = betaMq2(i1,i1) + diagonal(5,2)
    End Do

    diagonal(6,2) = 1.6_dp * gauge2(1) * (2._dp*sig(1)/3._dp - S2)           &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp       &
      &  + 3424._dp * gauge2(1) * AbsGM2(1) / 75._dp                         &
      &  + 512._dp * ( gauge2(1) * AbsGM2(3)                                 &
      &              + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)           &
      &                * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)-0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMu2 = -2._dp * ( MuYuaYuYuaYu + YuaYuYuaYuMu )                       &
      & - 4._dp * ( YuMqaYuYuaYu + YuaYuMuYuaYu + YuaYuYuMqaYu )             &
      & - 8._dp * MH(2) * YuaYuYuaYu                                         &
      & - 6._dp * TraceY(4) * ( MuYuaYu + YuaYuMu + 4._dp * MH(2) * YuaYu    &
      &                       + 2._dp * (YuMqaYu + AuaAu)  )                 &
      & - 4._dp * ( Tr3MqaYuYu3aYuMu * YuaYu + AuaAuYuaYu + YuaYuAuaAu       &
      &           + AuaYuYuaAu + YuaAuAuaYu + Tr3aAuAu * YuaYu               &
      &           + Tr3aAuYu * AuaYu + Tr3aYuAu * YuaAu )                    &
      & + hd(1) * ( MuYuaYu + YuaYuMu                                        &
      &           + 2._dp * ( MH(2) * YuaYu + YuMqaYu + AuaAu ) )            &
      & + hd(2) * YuaYu + hc(1) * AuaYu + hc(2) * YuaAu                      &
      & - 2._dp * (MuYuaYdYdaYu + YuaYdYdaYuMu )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YuaYdYdaYu + YuMqaYdYdaYu + YuaYdYdMqaYu &
      &           + YuaYdMdYdaYu + AuaAdYdaYu + YuaYdAdaAu                   &
      &           + AuaYdYdaAu + YuaAdAdaYu )
    Do i1=1,3
     betaMu2(i1,i1) = betaMu2(i1,i1) + diagonal(6,2)
    End Do

   End If 

  !------------------------------------------
  ! beta functions for Higgs mass parameters
  !------------------------------------------
   traceMH1(1) = mH(1) * (TraceY(1) + 6._dp * lam12)            &
             & + Real( cTrace(YeMlaYe),dp ) &
             & + Real( cTrace(aYeMeYe),dp ) + TraceA(1)
   traceMH1(2) = mH(1) * TraceY(3) + Real( cTrace(YdMqaYd),dp )   &
             & + Real( cTrace(aYdMdYd),dp ) + TraceA(3) + Alam12  &
             & + lam12 * MT(1)
   betamH11 = 6._dp * TraceMH1(2) + 2._dp * TraceMH1(1)      &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1

   traceMH2(1) = 3._dp * (2._dp*mH(2)+ MT(2)) * lam22 + 3._dp * Alam22
   traceMH2(2) = mH(2) * TraceY(4) + Real( cTrace(YuMqaYu),dp )  &
             & + Real( cTrace(aYuMuYu),dp ) + TraceA(4)
   betamH21 = 2._dp * TraceMH2(1) + 6._dp * TraceMH2(2)       &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) + 0.6_dp * S1

   If (decoupling_heavy_states) Then
    betaMT2 = 0._dp 

   Else
    betaMT2(1) = MT(1) * (lam12 + TraceY(2)) + 2._dp * mH(1) * lam12   &
            & + 2._dp * Real( cTrace(aYTMlYT),dp ) + TraceA(2) + Alam12 &
            & - 2.4_dp * AbsGM2(1) - 8._dp * AbsGM2(2)

    betaMT2(2) = MT(2) * lam22 + 2._dp * mH(2) * lam22 + Alam22 &
            & - 2.4_dp * AbsGM2(1) - 8._dp * AbsGM2(2)

    betaMT2 = 2._dp * betaMT2
   End If

   If (TwoLoopRGE) Then
    traceMH1(3) = MH(1) * (6._dp*TraceY2(2) + 2._dp*TraceY2(1) + TraceY2(4) ) &
              & + MH(2) * TraceY2(4)                                          &
              & + 6._dp * ( Real( cTrace(MqaYdYdaYdYd),dp )                   &
              &           + Real( cTrace(aYdMdYdaYdYd),dp )                   &
              &           + Real( cTrace(aAdAdaYdYd),dp )                     &
              &           + Real( cTrace(aAdYdaYdAd),dp )  )                  &
              & + 2._dp * ( Real( cTrace(MlaYeYeaYeYe),dp )                   &
              &           + Real( cTrace(aYeMeYeaYeYe),dp )                   &
              &           + Real( cTrace(aAeAeaYeYe),dp )                     &
              &           + Real( cTrace(aAeYeaYeAe),dp )  )                  &
       & + Real( cTrace(YdMqaYuYuaYd),dp ) + Real( cTrace(YdaYuMuYuaYd),dp ) &
       & + Real( cTrace(YdaYuYuMqaYd),dp ) + Real( cTrace(YuaYdMdYdaYu),dp ) &
       & + Real( cTrace(YdaAuAuaYd),dp ) + Real( cTrace(AdaYuYuaAd),dp )     &
              & + Real( cTrace(AdaAuYuaYd),dp ) + Real( cTrace(YdaYuAuaAd),dp )
    betaMH12 = - 6._dp * traceMH1(3)                                       &
      &   + (32._dp*gauge2(3) - 0.8_dp*gauge2(1) ) * traceMH1(2)           &
      &   + 64._dp * ( AbsGM2(3) * TraceY(3)                               &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(3),dp ) )     &
      &   - 1.6_dp * ( AbsGM2(1) * TraceY(3)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(3),dp ) )     &
      &   + 2.4_dp*gauge2(1) * traceMH1(1)                                 &
      &   + 4.8_dp * ( AbsGM2(1) * TraceY(1)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(1),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)                   &
      &   + 3._dp * gauge2(2) * sig(2)                                     &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)          &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2           &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

    traceMH2(1) = traceMH2(2)

    traceMH2(2) = MH(2) * ( 6._dp * TraceY2(3) + TraceY2(4) )                 &
      &  + MH(1) * TraceY2(4)                                                 &
 &  + 6._dp * ( Real(cTrace(MqaYuYuaYuYu),dp) + Real(cTrace(aYuMuYuaYuYu),dp) &
 &            + Real(cTrace(aAuAuaYuYu),dp) + Real(cTrace(aAuYuaYuAu),dp)  )  &
 &  + Real( cTrace(YuMqaYdYdaYu),dp ) + Real( cTrace(YuaYdMdYdaYu),dp )       &
 &  + Real( cTrace(YuaYdYdMqaYu),dp ) + Real( cTrace(YdaYuMuYuaYd),dp )       &
 &  + Real( cTrace(YuaAdAdaYu),dp ) + Real( cTrace(AuaYdYdaAu),dp )           &
      &  + Real( cTrace(AuaAdYdaYu),dp ) + Real( cTrace(YuaYdAdaAu),dp )
    betaMH22 = - 6._dp * traceMH2(2)                                          &
      &   + (32._dp*gauge2(3) + 1.6_dp*gauge2(1) ) * traceMH2(1)              &
      &   + 64._dp * ( AbsGM2(3) * TraceY(4)                                  &
      &          - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(4),dp ) )        &
      &   + 3.2_dp * ( AbsGM2(1) * TraceY(4)                                  &
      &          - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(4),dp ) )        &
      &   + gauge2(1) * ( 0.6_dp * sig(1) + 1.2_dp * S2)                      &
      &   + 3._dp * gauge2(2) * sig(2)                                        &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)             &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2              &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

   End If
!-----------------------------
! beta functions for mu and B
!-----------------------------
   TraceMue(1) = 3._dp * (TraceY(3)+TraceY(4)+lam12+lam22) + TraceY(1)  &
             & - 3._dp * gauge2(2) - 0.6_dp * gauge2(1)
   betaMue1 = mue * TraceMue(1)

   TraceB(1) = 6._dp * (TraceaYA(3)+TraceaYA(4) + lam1Alam1 + lam2Alam2) &
           & + 2._dp * TraceaYA(1) + 6._dp * g2Mi(2) + 1.2_dp * g2Mi(1)
   betaB1 = mue * TraceB(1) + B * TraceMue(1)

   If (TwoLoopRGE) Then
    TraceMue(2) = - 3._dp * ( 3._dp * (TraceY2(2) + TraceY2(3) )          &
      &                     + 2._dp * TraceY2(4) + TraceY2(1) )           &
      &         + (16._dp * gauge2(3) + 0.8_dp * gauge2(1)) * TraceY(4)   &
      &         + (16._dp * gauge2(3) - 0.4_dp * gauge2(1)) * TraceY(3)   &
      &         + 1.2_dp * gauge2(1) * TraceY(1)                          &
      &         + 7.5_dp * gauge2(2)**2                                   &
      &         + 1.8_dp * gauge2(2) * gauge2(1)                          &
      &         + 4.14_dp * gauge2(1)**2
    betaMue2 = mue * TraceMue(2)

    TraceB(1) = cTrace( 3._dp * ( Matmul(AuaYu,YuaYu) + Matmul(AdaYd,YdaYd) ) &
              &       + Matmul(AeaYe,YeaYe) + Matmul(aYuAu,aYdYd)             &
              &       + Matmul(aYdAd,aYuYu) ) 
    TraceB(2) = -12._dp * TraceB(1)                                           &
      &   + (32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)          &
      &   + (32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)          &
      &   + 2.4_dp * gauge2(1) * TraceaYA(1)                                  &
      &   - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)               &
      &   - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)               &
      &   - 2.4_dp * g2Mi(1) * TraceY(1) - 30._dp * gauge2(2)**2 * Mhlf(2)    &
      &   - 3.6_dp * gauge2(2) * gauge2(1) * (Mhlf(1) + Mhlf(2) )             &
      &   - 16.56_dp * gauge2(1)**2 * Mhlf(1) 
    betaB2 = mue * TraceB(2) + B * TraceMue(2)

   End If

  !--------------------------------
  ! neutrino dim. 5 operator
  !--------------------------------
  diagonal(5,1) = 6._dp * TraceY(4) - 1.2_dp * gauge2(1) - 6._dp * gauge2(2)
  betaMnu1 = Matmul( Transpose(aYeYe), Mnu) + Matmul(Mnu, aYeYe)  &
          & + diagonal(5,1) * Mnu
  

 !---------------
 ! 2-loop RGEs
 !---------------
  b_1a = b_1 + Delta_b_1
  b_2a = b_2 + Delta_b_2

  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1a + oo16pi2 * (Matmul(b_2a,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYT = oo16pi2 * betaYT1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   Do i1 = 1,3    
    sumI = 0._dp
    Do i2=1,3
     sumI = sumI + b_2a(i1,i2) * gauge2(i2) * (Mhlf(i1) + Mhlf(i2) ) 
    End Do
    Do i2=1,4
     sumI = sumI + a_2a(i1,i2) * ( TraceaYA(i2) - Mhlf(i1)*TraceY(i2) )
    End Do
    DMhlf(i1) = oo8pi2 * gauge2(i1) * ( b_1a(i1) * Mhlf(i1) + oo16pi2 * sumI)
   End Do
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * ( betaAe1 + oo16pi2 * betaAe2 )
   DAd = oo16pi2 * ( betaAd1 + oo16pi2 * betaAd2 )
   DAu = oo16pi2 * ( betaAu1 + oo16pi2 * betaAu2 )
   DAT = oo16pi2 * betaAT1
   DAlam1 = oo16pi2 * betaAlam11
   DAlam2 = oo16pi2 * betaAlam21
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * ( betaMe1 + oo16pi2 * betaMe2 )
   DMl = oo16pi2 * ( betaMl1 + oo16pi2 * betaMl2 )
   DMd = oo16pi2 * ( betaMd1 + oo16pi2 * betaMd2 )
   DMq = oo16pi2 * ( betaMq1 + oo16pi2 * betaMq2 )
   DMu = oo16pi2 * ( betaMu1 + oo16pi2 * betaMu2 )
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * ( betaMH11 + oo16pi2 * betaMH12 )
   DmH(2) = oo16pi2 * ( betaMH21 + oo16pi2 * betaMH22 )
   DmT2 = oo16pi2 * betaMT2
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * ( betaMue1 + oo16pi2 * betaMue2 )
   DB = oo16pi2 * ( betaB1 + oo16pi2 * betaB2 )

 !----------------------
 ! left neutrino mass
 !----------------------
   DMnu = oo16pi2 * betaMnu1
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1a 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYT = oo16pi2 * betaYT1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   DMhlf = oo8pi2 * gauge2 * b_1a * Mhlf
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * betaAe1
   DAd = oo16pi2 * betaAd1
   DAu = oo16pi2 * betaAu1
   DAT = oo16pi2 * betaAT1
   DAlam1 = oo16pi2 * betaAlam11
   DAlam2 = oo16pi2 * betaAlam21
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * betaMe1
   DMl = oo16pi2 * betaMl1
   DMd = oo16pi2 * betaMd1
   DMq = oo16pi2 * betaMq1
   DMu = oo16pi2 * betaMu1
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * betaMH11
   DmH(2) = oo16pi2 * betaMH21
   DmT2 = oo16pi2 * betaMT2
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * betaMue1
   DB = oo16pi2 * betaB1
  !----------------------
  ! left neutrino mass
  !----------------------
   DMnu = oo16pi2 * betaMnu1
  End If


  !---------------------------------------
  ! to avoid numerical problems in odeint
  !---------------------------------------
  Do i1=1,3
   DMe(i1,i1) = Real(DMe(i1,i1),dp)
   DMl(i1,i1) = Real(DMl(i1,i1),dp)
   DMd(i1,i1) = Real(DMd(i1,i1),dp)
   DMu(i1,i1) = Real(DMu(i1,i1),dp)
   DMq(i1,i1) = Real(DMq(i1,i1),dp)
  End Do
  Dmd = 0.5_dp * ( Dmd + Transpose(Conjg(Dmd)) )
  Dme = 0.5_dp * ( Dme + Transpose(Conjg(Dme)) )
  Dml = 0.5_dp * ( Dml + Transpose(Conjg(Dml)) )
  Dmq = 0.5_dp * ( Dmq + Transpose(Conjg(Dmq)) )
  Dmu = 0.5_dp * ( Dmu + Transpose(Conjg(Dmu)) )

  Call Chop(Dmue)
  Call Chop(DB)

  Call ParametersToG4(Dgauge, DYe, DYT, DYd, DYu, Dlam1, Dlam2, DMhlf, DAe, DAT &
          & , DAd, DAu, DAlam1, DAlam2, DMe, DMl, DMd, DMq, DMu, DMh, DMT2      &
          & , Dmue, DB, DMnu, f)

  Iname = Iname - 1

 End Subroutine rge277


 Subroutine rge285(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.03.2001: implementing right-handed neutrinos  at 1-loop
 !             up to now: Y_nu, A_nu,
 !             the parameters m_H1, M_H2, B and mu still need to be changed
 ! 08.01.06: including neutrino dim 5 operator
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2
  Real(dp) :: gauge(3), gauge2(3), TraceY(4), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2), sumI           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3), hc(4),Ynu(3,3) &
    & , aYnu(3,3), aYnuYnu(3,3), sumnu1(3,3), betaYnu1(3,3) ,DYnu(3,3)

  Complex(dp) :: Mhlf(3),DMhlf(3)

  Complex(dp) :: Ae(3,3), Ad(3,3), Au(3,3), aAe(3,3), aAd(3,3), aAu(3,3)   &
     &  , DAe(3,3), DAd(3,3), DAu(3,3), aAdAd(3,3), aAeAe(3,3), aAuAu(3,3) &
     &  , aYdAd(3,3), aYeAe(3,3), aYuAu(3,3), TraceaYA(4), betaAd1(3,3)    &
     &  , betaAd2(3,3), betaAe1(3,3), betaAe2(3,3), betaAu1(3,3)           &
     &  , betaAu2(3,3), Anu(3,3), aAnu(3,3), DAnu(3,3) ,aAnuAnu(3,3)       &
     &  , betaAnu1(3,3), aYnuAnu(3,3)!, betaAnu2(3,3)
  Real(dp) :: TraceA(4)
  Complex(dp) :: aYdYdaYdAd(3,3), aYdAdaYdYd(3,3), TraceAY2(5)               &
     &  , aYeYeaYeAe(3,3), aYeAeaYeYe(3,3), aYuYuaYuAu(3,3), aYuAuaYuYu(3,3) &
     &  , aYuYuaYdAd(3,3), aYuAuaYdYd(3,3), aYdYdaYuAu(3,3), aYdAdaYuYu(3,3)
   
  Complex(dp) :: Me(3,3), Ml(3,3), Md(3,3), Mq(3,3), Mu(3,3), DMe(3,3)        &
     & , DMl(3,3), DMd(3,3), DMq(3,3), DMu(3,3), YdaYd(3,3), YeaYe(3,3)       &
     & , YuaYu(3,3), MdYdaYd(3,3), MeYeaYe(3,3), MuYuaYu(3,3), YdaYdMd(3,3)   &
     & , YeaYeMe(3,3), YuaYuMu(3,3), YdMqaYd(3,3), YeMlaYe(3,3)               &
     & ,  YuMqaYu(3,3), AdaAd(3,3), AeaAe(3,3), AuaAu(3,3), betaMd1(3,3)      &
     & , betaMd2(3,3), betaMe1(3,3), betaMe2(3,3), betaMl1(3,3), betaMl2(3,3) &
     & , betaMq1(3,3), betaMq2(3,3), betaMu1(3,3), betaMu2(3,3), MqaYdYd(3,3) &
     & , MqaYuYu(3,3), aYdYdMq(3,3), aYuYuMq(3,3), aYeYeMl(3,3), MlaYeYe(3,3) &
     & , aYeMeYe(3,3), aYdMdYd(3,3), aYuMuYu(3,3)                             &
     & , YdaYdYdaYd(3,3), YeaYeYeaYe(3,3), YuaYuYuaYu(3,3), MeYeaYeYeaYe(3,3) &
     & , YeaYeYeaYeMe(3,3), YeaYeMeYeaYe(3,3), YeMlaYeYeaYe(3,3), AeaYe(3,3)  &
     & , YeaYeYeMlaYe(3,3), AeaAeYeaYe(3,3), YeaYeAeaAe(3,3), YeaAe(3,3)      &
     & , AeaYeYeaAe(3,3), YeaAeAeaYe(3,3), Tr3aAdYdaAeYe

  Complex(dp) :: Tr3aYdAdaYeAe, AdaYd(3,3), YdaAd(3,3), MlaYeYeaYeYe(3,3)     &
     & , aYeYeaYeYeMl(3,3), aYeYeMlaYeYe(3,3), aYeYeaYeMeYe(3,3)              &
     & , aYeMeYeaYeYe(3,3), aAdYd(3,3),aAeYe(3,3), aAeAeaYeYe(3,3)            &
     & , aYeYeaAeAe(3,3), aAeYeaYeAe(3,3), aYeAeaAeYe(3,3), MdYdaYdYdaYd(3,3) &
     & , YdaYdYdaYdMd(3,3), YdMqaYdYdaYd(3,3),YdaYdMdYdaYd(3,3)               &
     & , YdaYdYdMqaYd(3,3), AdaAdYdaYd(3,3), YdaYDAdaAd(3,3), AdaYdYdaAd(3,3) &
     & , YdaAdAdaYd(3,3)
  Complex(dp) :: MdYdaYuYuaYd(3,3), YdaYuYuaYdMd(3,3), YdMqaYuYuaYd(3,3)      &
     & , YdaYuYuMqaYd(3,3), YdaYuMuYuaYd(3,3), AdaAuYuaYd(3,3)                &
     & , YdaYuAuaAd(3,3), AdaYuYuaAd(3,3), YdaAuAuaYd(3,3), YdaYuYuaYd(3,3)   &
     & , Tr3aYuAu, Tr3aAuYu, YuaAu(3,3)                                       &
     & , MqaYdYdaYdYd(3,3), aYdYdaYdYDMq(3,3), aYdMdYdaYdYd(3,3)              &
     & , aYdYdMqaYdYd(3,3), aYdYdaYdMdYd(3,3), aAdAdaYdYd(3,3)                &
     & , aYdYDaAdAd(3,3), aAdYdaYdAd(3,3), aYdAdaAdYd(3,3), MqaYuYuaYuYu(3,3) &
     & , aYuYuaYuYUMq(3,3), aYuMuYuaYuYu(3,3), aYuYuMqaYuYu(3,3)              &
     & , aYuYuaYuMuYu(3,3), aAuAuaYuYu(3,3), aYuYUaAuAu(3,3), aAuYuaYuAu(3,3) &
     & , aYuAuaAuYu(3,3), aAuYu(3,3), AuaYu(3,3), YuaYdYdaYu(3,3)             &
     & , AuaYdYdaAu(3,3), YuaAdAdaYu(3,3), AuaAdYdaYu(3,3), YuaYdAdaAu(3,3)   &
     & , YuMqaYuYuaYu(3,3), YuaYuYuMqaYu(3,3), MuYuaYuYuaYu(3,3)              &
     & , YuaYuYuaYuMu(3,3), YuaYuMuYuaYu(3,3), AuaAuYuaYu(3,3)                &
     & , YuaYuAuaAu(3,3), AuaYuYuaAu(3,3), YuaAuAuaYu(3,3), MuYuaYdYdaYu(3,3) &
     & , YuaYdYdaYuMu(3,3), YuMqaYdYdaYu(3,3), YuaYdYdMqaYu(3,3)              &
     & , YuaYdMdYdaYu(3,3)
  Complex(dp) :: Mr(3,3), YnuaYnu(3,3), MlaYnuYnu(3,3), MrYnuaYnu(3,3)        &
     & , aYnuYnuMl(3,3), YnuaYnuMr(3,3), YnuMlaYnu(3,3), aYnuMrYnu(3,3)       &
     & , AnuaAnu(3,3), DMr(3,3), betaMr1(3,3)

  Real(dp) :: S1, S2, sig(3), Tr3aYdYdaYeYe, Tr3aAdAdaAeAe, AbsGM2(3)         &
     & , Tr3MqaYdYd3aYDMdYd , Tr3MqaYuYu3aYuMu, Tr3aAuAu

  Real(dp) :: Mh(2), DMh(2), TraceMH1(3), TraceMH2(2), betaMH11        &
     &  , betaMH12, betaMH21, betaMH22, q

  Complex(dp) :: mue, B, Dmue, DB, TraceMue(2), TraceB(2), betaMue1, betaMue2 &
     & , betaB1, betaB2, g2Mi(3)

  Complex(dp) :: Mnu(3,3), DMnu(3,3), betamnu1(3,3), sumM1(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge285'

  OnlyDiagonal = .Not.GenerationMixing
  q = t

  Call GToParameters3(gy, gauge, Ye, Ynu, Yd, Yu, Mhlf, Ae, Anu, Ad, Au      &
                & , Me, Ml, Mr, Md, Mq, Mu, Mh, mue, B, Mnu)

  gauge2 = gauge**2
  AbsGM2 = gauge2 * Abs( Mhlf )**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ynu,aYnu)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul(aYe,Ye)
  aYnuYnu = Matmul(aYnu,Ynu)
  aYdYd = Matmul(aYd,Yd)
  aYuYu = Matmul(aYu,Yu)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYnuYnu),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = 3._dp * TraceY(3) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe + aYnuYnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = 3._dp * TraceY(4) + TraceY(2)           &
            &   - 0.6_dp * gauge2(1) - 3._dp * gauge2(2)
  sumnu1 = 3._dp * aYnuYnu + aYeYe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do

  betaYnu1 = Matmul(Ynu,sumnu1)

  diagonal(3,1) = 3._dp * TraceY(3) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)

  diagonal(4,1) = 3._dp * TraceY(4)  + TraceY(2)             &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
   aYdYdaYuYu = Matmul(aYdYd,aYuYu)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ae,aAe)
  Call Adjungate(Anu,aAnu)
  Call Adjungate(Ad,aAd)
  Call Adjungate(Au,aAu)

  aAdAd = Matmul(aAd,Ad)
  aAnuAnu = Matmul(aAnu,Anu)
  aAeAe = Matmul(aAe,Ae)
  aAuAu = Matmul(aAu,Au)

  TraceA(1) = Real( cTrace(aAeAe),dp )
  TraceA(2) = Real( cTrace(aAnuAnu),dp )
  TraceA(3) = Real( cTrace(aAdAd),dp )
  TraceA(4) = Real( cTrace(aAuAu),dp )

  aYdAd = Matmul(aYd,Ad)
  aYnuAnu = Matmul(aYnu,Anu)
  aYeAe = Matmul(aYe,Ae)
  aYuAu = Matmul(aYu,Au)

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYnuAnu) 
  TraceaYA(3) = cTrace(aYdAd) 
  TraceaYA(4) = cTrace(aYuAu) 

  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = Matmul(Ae,sume1)
  
  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1)  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    )
  sume1 = 4._dp * aYeAe + 2._dp * aYnuAnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + Matmul(Ye,sume1)

  !--------------
  ! A_nu
  !--------------
  sumnu1 = sumnu1 + 2._dp * aYnuYnu
  betaAnu1 = Matmul(Anu,sumnu1)

  diagonal(2,1) = 2._dp * ( 3._dp * TraceaYA(4) + TraceaYA(2)     &
      &         + 0.6_dp * g2Mi(1) + 3._dp * g2Mi(2) )
  sumnu1 = 4._dp * aYnuAnu + 2._dp * aYeAe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do
  betaAnu1 = betaAnu1 + Matmul(Ynu,sumnu1)

  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = Matmul(Ad,sumd1)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1) &
                &         - c1_1(2,1) * g2Mi(1) - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do
  betaAd1 = betaAd1 + Matmul(Yd,sumd1)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = Matmul(Au,sumu1)
  
  diagonal(4,1) = 2._dp * ( 3._dp * TraceaYA(4) + TraceaYA(2) &
                &         - c1_1(3,1) * g2Mi(1) - c1_1(3,2) * g2Mi(2)   &
                &         - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do
  betaAu1 = betaAu1 + Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdAd = Matmul(aYdYd,aYdAd)
   aYdAdaYdYd = Matmul(aYdAd,aYdYd)
   aYeYeaYeAe = Matmul(aYeYe,aYeAe)
   aYeAeaYeYe = Matmul(aYeAe,aYeYe)
   aYuYuaYuAu = Matmul(aYuYu,aYuAu)
   aYuAuaYuYu = Matmul(aYuAu,aYuYu)
   aYuAuaYdYd = Matmul(aYuAu,aYdYd)
   aYuYuaYdAd = Matmul(aYuYu,aYdAd)
   aYdAdaYuYu = Matmul(aYdAd,aYuYu)
   aYdYdaYuAu = Matmul(aYdYd,aYuAu)
   TraceAY2(1) = cTrace(aYeYeaYeAe)
   TraceAY2(2) = cTrace(aYdYdaYdAd)
   TraceAY2(3) = cTrace(aYuYuaYuAu)
   TraceAY2(4) = cTrace(aYuYuaYdAd)
   TraceAY2(5) = cTrace(aYdYdaYuAu)

  !--------------
  ! A_e
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = Matmul(Ae,sume2)
    
   diagonal(1,2) = -6._dp * ( 6._dp * TraceAY2(2) + TraceAY2(4)       &
     &                      + TraceAY2(5) + 2._dp * TraceAY2(1)  )    &
     &  + ( 32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)   &
     &  - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)                &
     &  + 2.4_dp * gauge2(1) * TraceaYA(1) - 2.4_dp * g2Mi(1) * TraceY(1)    &
     &  - ( 30._dp * g2Mi(2)                              &
     &    + 3.6_dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)      &
     &  - 54._dp * gauge2(1)**2 * Mhlf(1)
   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)     &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1) &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2)
   sume2 = - 6._dp * aYeYeaYeAe - 8._dp * aYeAeaYeYe  &
         & - hd(1) * aYeAe - hc(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
   betaAe2 = betaAe2 + Matmul(Ye,sume2)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = Matmul(Ad,sumd2)
    
   diagonal(3,2) = diagonal(1,2)                                   &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                 &
     &              - gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp      &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)    &
     &  + 1.6_dp * gauge2(1) * gauge2(2) * (Mhlf(1)+Mhlf(2))       & 
     &  + 1.856e3_dp * gauge2(1)**2 * Mhlf(1) / 4.5e1_dp

   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)  &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hd(2) = 6._dp * TraceY(4) - 1.6_dp * gauge2(1) 
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1)                 &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2) + 1.6_dp * gauge2(1) *  Mhlf(1)
   hc(2) = 6._dp * TraceaYA(4) + 1.6_dp * gauge2(1) *  Mhlf(1)
   sumd2 = - 6._dp * aYdYdaYdAd - 8._dp * aYdAdaYdYd                  &
       &   - 4._dp * ( aYuAuaYuYu + aYuYuaYuAu + aYuAuaYdYd )         &
       &   - 2._dp * aYuYuaYdAd - hd(1) * aYdAd - hc(1) * aYdYd       &
       &  - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(3,2)
   End Do
   betaAd2 = betaAd2 + Matmul(Yd,sumd2)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(4) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = Matmul(Au,sumu2)
    
   diagonal(4,2) =  -6._dp * ( 6._dp * TraceAY2(3) + TraceAY2(4)        &
     &                        + TraceAY2(5)  )                          &
     &  + ( 32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)     &
     &  - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)           &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                      &
     &              - 3.4_dp * gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp  &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)         &
     &  - ( 30._dp * g2Mi(2)                                &
     &    + 2._dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)         &
     &  - 5486._dp * gauge2(1)**2 * Mhlf(1) / 225._dp
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1) - 0.8_dp * gauge2(1)
   hc(1) = 6._dp * TraceaYA(3) + 2._dp * TraceaYA(1)   &
       &  + 0.8_dp * gauge2(1) *  Mhlf(1)
   hd(2) = 12._dp * TraceY(4) - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(2) = 18._dp * TraceaYA(4) + 1.2e1_dp * g2Mi(2) + 0.8_dp * g2Mi(1)
   sumu2 = - 6._dp * aYuYuaYuAu - 8._dp * aYuAuaYuYu                     &
       &   - 4._dp * ( aYdAdaYdYd + aYdYdaYdAd + aYdAdaYuYu )            &
       &   - 2._dp * aYdYdaYuAu - hd(1) * aYdAd - hc(1) * aYdYd          &
       &   - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(4,2)
   End Do
   betaAu2 = betaAu2 + Matmul(Yu,sumu2)

  End If 
!----------------------------------------------
! beta functions for Sfermion mass parameters
!----------------------------------------------
   S1 = mH(2) - mH(1)
   Do i1=1,3
    S1 = S1 + Me(i1,i1) - Ml(i1,i1) &
       &    + Md(i1,i1) + Mq(i1,i1) - 2._dp * Mu(i1,i1)
   End Do
   S1 = S1 * gauge2(1)

   YdaYd = Matmul(Yd,aYd)
   YnuaYnu = Matmul(Ynu,aYnu)
   YeaYe = Matmul(Ye,aYe)
   YuaYu = Matmul(Yu,aYu)

   MeYeaYe = Matmul(Me,YeaYe)
   MlaYeYe = Matmul(Ml,aYeYe)
   MlaYnuYnu = Matmul(Ml,aYnuYnu)
   MrYnuaYnu = Matmul(Mr,YnuaYnu)

   MdYdaYd = Matmul(Md,YdaYd)
   MqaYdYd = Matmul(Mq,aYdYd)
   MqaYuYu = Matmul(Mq,aYuYu)
   MuYuaYu = Matmul(Mu,YuaYu)

   YeaYeMe = Matmul(YeaYe,Me)
   aYeYeMl = Matmul(aYeYe,Ml)
   aYnuYnuMl = Matmul(aYnuYnu,Ml)
   YnuaYnuMr = Matmul(YnuaYnu,Mr)

   YdaYdMd = Matmul(YdaYd,Md)
   aYdYdMq = Matmul(aYdYd,Mq)
   aYuYuMq = Matmul(aYuYu,Mq)
   YuaYuMu = Matmul(YuaYu,Mu)

   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   YnuMlaYnu = MatMul3(Ynu,Ml,aYnu,OnlyDiagonal)
   aYnuMrYnu = MatMul3(aYnu,Mr,Ynu,OnlyDiagonal)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)

   AeaAe = Matmul(Ae,aAe)
   AnuaAnu = Matmul(Anu,aAnu)
   AdaAd = Matmul(Ad,aAd)
   AuaAu = Matmul(Au,aAu)

   diagonal(1,1) = - 4.8_dp * AbsGM2(1) + 1.2_dp * S1
   betaMe1 = 2._dp * (MeYeaYe + YeaYeMe)             &
         & + 4._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe )
   Do i1=1,3
    betaMe1(i1,i1) = betaMe1(i1,i1) + diagonal(1,1)
   End Do

   diagonal(3,1) = - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1
   betaMl1 = MlaYeYe + aYeYeMl + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe )  &
         & + MlaYnuYnu + aYnuYnuMl                                            &
         & + 2._dp * ( mH(2) * aYnuYnu + aYnuMrYnu + aAnuAnu )
   Do i1=1,3
    betaMl1(i1,i1) = betaMl1(i1,i1) + diagonal(3,1)
   End Do

   betaMr1 = 2._dp * (MrYnuaYnu + YnuaYnuMr)                        &
         & + 4._dp * ( mH(2) * YnuaYnu + YnuMlaYnu + AnuaAnu )

   diagonal(4,1) = - ( 32._dp * AbsGM2(3) + 1.6_dp * AbsGM2(1) ) / 3._dp &
               & + 0.4_dp * S1
   betaMd1 = 2._dp * (MdYdaYd + YdaYdMd)             &
         & + 4._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd )
   Do i1=1,3
    betaMd1(i1,i1) = betaMd1(i1,i1) + diagonal(4,1)
   End Do

   diagonal(5,1) = - ( 32._dp * AbsGM2(3) + 0.4_dp * AbsGM2(1) ) / 3._dp &
               & - 6._dp * AbsGM2(2) + 0.2_dp * S1
   betaMq1 = MqaYuYu + aYuYuMq + MqaYdYd + aYdYdMq             &
         & + 2._dp * ( mH(2) * aYuYu + mH(1) * aYdYd + aYuMuYu     &
         &           + aYdMdYd + aAuAu + aAdAd )
   Do i1=1,3
    betaMq1(i1,i1) = betaMq1(i1,i1) + diagonal(5,1)
   End Do

   diagonal(6,1) = - ( 32._dp * AbsGM2(3) + 6.4_dp * AbsGM2(1) ) / 3._dp &
               &   - 0.8_dp * S1
   betaMu1 = 2._dp * (MuYuaYu + YuaYuMu)             &
         & + 4._dp * ( mH(2) * YuaYu + YuMqaYu + AuaAu )
   Do i1=1,3
    betaMu1(i1,i1) = betaMu1(i1,i1) + diagonal(6,1)
   End Do

   If (TwoLoopRGE) Then
    YdaYdYdaYd = MatSquare(YdaYd,OnlyDiagonal)
    YeaYeYeaYe = MatSquare(YeaYe,OnlyDiagonal)
    YuaYuYuaYu = MatSquare(YuaYu,OnlyDiagonal)

    AdaYd = Matmul(Ad,aYd)
    AeaYe = Matmul(Ae,aYe)
    AuaYu = Matmul(Au,aYu)

    aAdYd = Matmul(aAd,Yd)
    aAeYe = Matmul(aAe,Ye)
    aAuYu = Matmul(aAu,Yu)

    YdaAd = Matmul(Yd,aAd)
    YeaAe = Matmul(Ye,aAe)
    YuaAu = Matmul(Yu,aAu)

    YdaYuYuaYd = MatMul3(Yd,aYuYu,aYd,OnlyDiagonal)
    AdaYuYuaAd = MatMul3(Ad,aYuYu,aAd,OnlyDiagonal)
    YdaAuAuaYd = MatMul3(Yd,aAuAu,aYd,OnlyDiagonal)
    AdaAuYuaYd = MatMul4(Ad,aAu,Yu,aYd,OnlyDiagonal)
    YdaYuAuaAd = MatMul3(Yd,aYuAu,aAd,OnlyDiagonal)

    YuaYdYdaYu = MatMul3(Yu,aYdYd,aYu,OnlyDiagonal)
    AuaYdYdaAu = MatMul3(Au,aYdYd,aAu,OnlyDiagonal)
    YuaAdAdaYu = MatMul3(Yu,aAdAd,aYu,OnlyDiagonal)
    AuaAdYdaYu = MatMul4(Au,aAd,Yd,aYu,OnlyDiagonal)
    YuaYdAdaAu = MatMul3(Yu,aYdAd,aAu,OnlyDiagonal)

    MdYdaYuYuaYd = Matmul(Md,YdaYuYuaYd)
    Call Adjungate(MdYdaYuYuaYd, YdaYuYuaYdMd)
    YdMqaYuYuaYd = MatMul3(Yd,MqaYuYu,aYd,OnlyDiagonal)
    Call Adjungate(YdMqaYuYuaYd, YdaYuYuMqaYd)
    YdaYuMuYuaYd = MatMul3(Yd,aYuMuYu,aYd,OnlyDiagonal)

    MuYuaYdYdaYu = Matmul(Mu,YuaYdYdaYu)
    Call Adjungate(MuYuaYdYdaYu, YuaYdYdaYuMu)
    YuMqaYdYdaYu = MatMul3(Yu,MqaYdYd,aYu,OnlyDiagonal)
    Call Adjungate(YuMqaYdYdaYu, YuaYdYdMqaYu)
    YuaYdMdYdaYu = MatMul3(Yu,aYdMdYd,aYu,OnlyDiagonal)

    MeYeaYeYeaYe = Matmul(MeYeaYe,YeaYe)
    Call Adjungate(MeYeaYeYeaYe,YeaYeYeaYeMe)
    aYeMeYeaYeYe = Matmul(aYeMeYe,aYeYe)
    Call Adjungate(aYeMeYeaYeYe,aYeYeaYeMeYe)
    YeaYeMeYeaYe = Matmul(YeaYeMe,YeaYe)

    MlaYeYeaYeYe = Matmul(MlaYeYe,aYeYe)
    Call Adjungate(MlaYeYeaYeYe, aYeYeaYeYeMl)
    YeMlaYeYeaYe = Matmul(YeMlaYe,YeaYe)
    Call Adjungate(YeMlaYeYeaYe, YeaYeYeMlaYe)
    aYeYeMlaYeYe = Matmul(aYeYeMl,aYeYe)

    MdYdaYdYdaYd = Matmul(MdYdaYd,YdaYd)
    Call Adjungate(MdYdaYdYdaYd, YdaYdYdaYdMd)
    aYdMdYdaYdYd = Matmul(aYdMdYd,aYdYd)
    Call Adjungate(aYdMdYdaYdYd, aYdYdaYdMdYd)
    YdaYdMdYdaYd = Matmul(YdaYdMd,YdaYd)

    MqaYdYdaYdYd = Matmul(MqaYdYd,aYdYd)
    Call Adjungate(MqaYdYdaYdYd, aYdYdaYdYdMq)
    YdMqaYdYdaYd = Matmul(YdMqaYd,YdaYd)
    Call Adjungate(YdMqaYdYdaYd, YdaYdYdMqaYd)
    aYdYdMqaYdYd = Matmul(aYdYdMq,aYdYd)

    MqaYuYuaYuYu = Matmul(MqaYuYu,aYuYu)
    Call Adjungate(MqaYuYuaYuYu, aYuYuaYuYuMq)
    YuMqaYuYuaYu = Matmul(YuMqaYu,YuaYu)
    Call Adjungate(YuMqaYuYuaYu, YuaYuYuMqaYu)
    aYuYuMqaYuYu = Matmul(aYuYuMq,aYuYu)

    MuYuaYuYuaYu = Matmul(MuYuaYu,YuaYu)
    Call Adjungate(MuYuaYuYuaYu, YuaYuYuaYuMu)
    aYuMuYuaYuYu = Matmul(aYuMuYu,aYuYu)
    Call Adjungate(aYuMuYuaYuYu, aYuYuaYuMuYu)
    YuaYuMuYuaYu = Matmul(YuaYuMu,YuaYu)

    AdaAdYdaYd = Matmul(AdaAd,YdaYd)
    Call Adjungate(AdaAdYdaYd, YdaYdAdaAd)
    AdaYdYdaAd = Matmul(AdaYd,YdaAd)
    YdaAdAdaYd = Matmul(YdaAd,AdaYd)

    aAdAdaYdYd = Matmul(aAdAd,aYdYd)
    Call Adjungate(aAdAdaYdYd, aYdYdaAdAd)
    aAdYdaYdAd = Matmul(aAdYd,aYdAd)
    aYdAdaAdYd = Matmul(aYdAd,aAdYd)

    AeaAeYeaYe = Matmul(AeaAe,YeaYe)
    Call Adjungate(AeaAeYeaYe, YeaYeAeaAe)
    AeaYeYeaAe = Matmul(AeaYe,YeaAe)
    YeaAeAeaYe = Matmul(YeaAe,AeaYe)

    aAeAeaYeYe = Matmul(aAeAe,aYeYe)
    Call Adjungate(aAeAeaYeYe, aYeYeaAeAe)
    aAeYeaYeAe = Matmul(aAeYe,aYeAe)
    aYeAeaAeYe = Matmul(aYeAe,aAeYe)

    AuaAuYuaYu = Matmul(AuaAu,YuaYu)
    Call Adjungate(AuaAuYuaYu, YuaYuAuaAu)
    AuaYuYuaAu = Matmul(AuaYu,YuaAu)
    YuaAuAuaYu = Matmul(YuaAu,AuaYu)

    aAuAuaYuYu = Matmul(aAuAu,aYuYu)
    Call Adjungate(aAuAuaYuYu, aYuYuaAuAu)
    aAuYuaYuAu = Matmul(aAuYu,aYuAu)
    aYuAuaAuYu = Matmul(aYuAu,aAuYu)

    S2 = (1.5_dp * gauge2(2) + 0.3_dp * gauge2(1) )            &
     &      * (MH(2) - MH(1) - Real(cTrace(ML),dp) )              &
     & + ( (8._dp * gauge2(3) + 0.1_dp*gauge2(1)) / 3._dp      &
     &   + 1.5_dp * gauge2(2) ) * Real( cTrace(Mq),dp   )         &
     & - (16._dp * gauge2(3) + 3.2_dp*gauge2(1) )              &
     &    * Real(cTrace(Mu),dp) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md),dp) / 3._dp                           &
     & + 1.2_dp*gauge2(1) * Real(cTrace(Me),dp)                   &
     & - 3._dp * (MH(2)*TraceY(4) - MH(1) * TraceY(3) )        &
     & + mH(1) * TraceY(1)
    
    Do i1=1,3
     S2 = S2 - YuMqaYu(i1,i1) + 4._dp * aYuMuYu(i1,i1)   &
        &    - YdMqaYd(i1,i1) - 2._dp * aYdMdYd(i1,i1)   &
        &    + YeMlaYe(i1,i1) - 2._dp * aYeMeYe(i1,i1)
    End Do

    sig(1) = 3._dp * (MH(1) + MH(2) + Real(cTrace(Ml),dp) )          &
         & + Real(cTrace(Mq),dp) + 8._dp * Real(cTrace(Mu),dp)          &
         & + 2._dp * Real(cTrace(Md),dp) + 6._dp * Real(cTrace(Me),dp) 
    sig(1) = 0.2_dp * gauge2(1) * sig(1)
    sig(2) = gauge2(2) * ( MH(1) + MH(2) + Real( cTrace(Ml),dp )    &
           &             + 3._dp * Real( cTrace(Mq),dp ) )
    sig(3) = gauge2(3) * ( 2._dp * Real( cTrace(Mq),dp )            &
           &             + Real(cTrace(Mu),dp) + Real(cTrace(Md),dp) )

    Tr3aYdYdaYeYe = 3._dp * TraceY(3) + TraceY(1)
    Tr3MqaYdYd3aYDMdYd = Real(cTrace(MlaYeYe),dp) + Real(cTrace(aYeMeYe),dp)  &
            & + 3._dp * ( Real(cTrace(MqaYdYd),dp) + Real(cTrace(aYdMdYd),dp) )
    Tr3aAdAdaAeAe = 3._dp * TraceA(3) + TraceA(1)
    Tr3aYdAdaYeAe = 3._dp * TraceaYA(3) + TraceaYA(1)
    Tr3aAdYdaAeYe = Conjg( Tr3aYdAdaYeAe )
    Tr3MqaYuYu3aYuMu = 3._dp * (Real(cTrace(MqaYuYu),dp) &
                               + Real(cTrace(aYuMuYu),dp) )
    Tr3aAuAu = 3._dp * TraceA(4)
    Tr3aYuAu = 3._dp * TraceaYA(4)
    Tr3aAuYu = Conjg( Tr3aYuAu )

    diagonal(1,2) = 2.4_dp * gauge2(1) * (S2 + sig(1) )   &
                & + 112.32_dp * gauge2(1) * AbsGM2(1)
    hd(1) = 6._dp*gauge2(2)-1.2_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 4.8_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 2.4_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMe2 = -2._dp * ( MeYeaYeYeaYe + YeaYeYeaYeMe )                       &
      & - 4._dp * ( YeMlaYeYeaYe + YeaYeMeYeaYe + YeaYeYeMlaYe )             &
      & - 8._dp * mH(1) * YeaYeYeaYe                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MeYeaYe + YeaYeMe + 4._dp * MH(1) * YeaYe &
      &                           + 2._dp * (YeMlaYe + AeaAE)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YeaYe + AeaAeYeaYe + YeaYEAeaAe     &
      &           + AeaYeYeaAe + YeaAeAeaYe + Tr3aAdAdaAeAe * YeaYe          &
      &           + Tr3aAdYdaAeYe * AeaYe + Tr3aYdAdaYeAe * YeaAe )          &
      & + hd(1) * ( MeYeaYe + YeaYeMe                                        &
      &           + 2._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe ) )            &
      & + hd(2) * YeaYe + hc(1) * AeaYe + hc(2) * YeaAe
    Do i1=1,3
     betaMe2(i1,i1) = betaMe2(i1,i1) + diagonal(1,2)
    End Do

    diagonal(3,2) = gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)        &
     &     + 3._dp * gauge2(2) * sig(2)                                 &
     &     + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)      &
     &                   + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2       &
     &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )     &
     &     + 24.84_dp * gauge2(1) * AbsGM2(1)
    betaMl2 = -2._dp * ( MlaYeYeaYeYe + aYeYeaYeYeMl )                   &
      & - 4._dp * ( aYeMeYeaYeYe + aYeYeMlaYeYe + aYeYeaYeMeYe )         &
      & - 8._dp * mH(1) * aYeYeaYeYe                                     &
      & - Tr3aYdYdaYeYe * ( MlaYeYe + aYeYeML + 4._dp * MH(1) * aYeYe    &
      &                   + 2._dp * (aYeMeYe + aAeAE)  )                 &
      & - 4._dp * (aAeAeaYeYe + aYeYeaAeAe + aAeYeaYeAe + aYeAeaAeYe )   &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYeYe +  Tr3aAdAdaAeAe * aYeYe  &
      &           + Tr3aAdYdaAeYe * aYeAe + Tr3aYdAdaYeAe * aAeYe     )  &
      & + 1.2_dp*gauge2(1) * ( MlaYeYe + aYeYeMl                         &
      &                      + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe &
      &                                - Mhlf(1) * aAeYe                 &
      &                                - Conjg(Mhlf(1)) * aYeAe ) )      &
      & + 4.8_dp * AbsGM2(1) * aYeYe
    Do i1=1,3
     betaMl2(i1,i1) = betaMl2(i1,i1) + diagonal(3,2)
    End Do

    diagonal(4,2) = 0.8_dp * gauge2(1) * (S2 + sig(1)/3._dp )             &
      &   + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp   &
      &   + 8.08e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                    &
      &   + 1.28e2_dp * ( gauge2(1) * AbsGM2(3)                           &
      &                 + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)     &
      &                   * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)+0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) + 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) - 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMd2 = -2._dp * ( MdYdaYdYdaYd + YdaYdYdaYdMd )                       &
      & - 4._dp * ( YdMqaYdYdaYd + YdaYdMdYdaYd + YdaYdYdMqaYd )             &
      & - 8._dp * mH(1) * YdaYdYdaYd                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MdYdaYd + YdaYdMd + 4._dp * MH(1) * YdaYd &
      &                           + 2._dp * (YdMqaYd + AdaAD)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YdaYd + AdaAdYdaYd + YdaYDAdaAd     &
      &           + AdaYdYdaAd + YdaAdAdaYd + Tr3aAdAdaAeAe * YdaYd          &
      &           + Tr3aAdYdaAeYe * AdaYd + Tr3aYdAdaYeAe * YdaAd )          &
      & + hd(1) * ( MdYdaYd + YdaYdMd                                        &
      &           + 2._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd ) )            &
      & + hd(2) * YdaYd + hc(1) * AdaYd + hc(2) * YdaAd                      &
      & - 2._dp * (MdYdaYuYuaYd + YdaYuYuaYdMd )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YdaYuYuaYd + YdMqaYuYuaYd+ YdaYuYuMqaYd  &
      &           + YdaYuMuYuaYd + AdaAuYuaYd + YdaYuAuaAd                   &
      &           + AdaYuYuaAd + YdaAuAuaYd )
    Do i1=1,3
     betaMd2(i1,i1) = betaMd2(i1,i1) + diagonal(4,2)
    End Do

    diagonal(5,2) = 0.2_dp * gauge2(1) * (2._dp * S2 + sig(1)/3._dp )     &
      &  + 3._dp * gauge2(2) * sig(2)                                     &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp    &
      &  + 1.99e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                     &
      &  + 32._dp * ( gauge2(1) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)         &
      &              * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp      &
      &  + 33._dp * gauge2(2) * AbsGM2(2)                                 &
      &  + 32._dp * ( gauge2(2) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(2)  + gauge2(3)*gauge2(2)        &
      &              * Real( Mhlf(2) * Conjg(Mhlf(3)),dp ) )                 &
      &  + 0.4_dp * ( gauge2(2) * AbsGM2(1)                               &
      &             + gauge2(1) * AbsGM2(2) + gauge2(1)*gauge2(2)         &
      &              * Real( Mhlf(2) * Conjg(Mhlf(1)),dp ) )
    hd(1) = 1.6_dp * AbsGM2(1)
    hc(1) = - 0.8_dp * g2Mi(1)
    hc(2) = Conjg( hc(1) )
    hd(2) = 3.2_dp * AbsGM2(1)
    hc(3) = - 1.6_dp * g2Mi(1)
    hc(4) = Conjg( hc(3) )
    betaMq2 = -2._dp * ( MqaYdYdaYdYd + aYdYdaYdYDMq )                       &
      & - 4._dp * ( aYdMdYdaYdYd + aYdYdMqaYdYd + aYdYdaYdMdYd )             &
      & - 8._dp * mH(1) * aYdYdaYdYd                                         &
      & - Tr3aYdYdaYeYe * ( MqaYdYd + aYdYdMq + 4._dp * MH(1) * aYdYd        &
      &                   + 2._dp * (aYdMdYd + aAdAd)  )                     &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYdYd                               &
      &           + 2._dp * ( aAdAdaYdYd + aYdYDaAdAd + aAdYdaYdAd           &
      &                     + aYdAdaAdYd )                                   &
      &           + Tr3aAdAdaAeAe * aYdYd + Tr3aAdYdaAeYe * aYdAd            &
      &           + Tr3aYdAdaYeAe * aAdYd )                                  &
      & + 0.4_dp*gauge2(1) * ( MqaYdYd + aYdYdMq                             &
      &                      + 2._dp * ( mH(1) * aYdYd + aYdMdYd + aAdAd ) ) &
      & + hd(1) * aYdYd + hc(1) * aAdYd + hc(2) * aYdAd                      &
      & - 2._dp * ( MqaYuYuaYuYu + aYuYuaYuYuMq )                            &
      & - 4._dp * ( aYuMuYuaYuYu + aYuYuMqaYuYu + aYuYuaYuMuYu )             &
      & - 8._dp * MH(2) * aYuYuaYuYu                                         &
      & - 3._dp * TraceY(4) * ( MqaYuYu + aYuYuMq + 4._dp * MH(2) * aYuYu    &
      &                       + 2._dp * (aYuMuYu + aAuAu) )                  &
      & - 2._dp * ( Tr3MqaYuYu3aYuMu * aYuYu                                 &
      &           + 2._dp * ( aAuAuaYuYu + aYuYuaAuAu                        &
      &                     + aAuYuaYuAu + aYuAuaAuYu )                      &
      &           + Tr3aAuAu * aYuYu + Tr3aAuYu * aYuAu + Tr3aYuAu * aAuYu ) &
      & + 0.8_dp*gauge2(1) * ( MqaYuYu + aYuYuMq                             &
      &                      + 2._dp * ( MH(2) * aYuYu + aYuMuYu + aAuAu ) ) &
      & + hd(2) * aYuYu + hc(3) * aAuYu + hc(4) * aYuAu
    Do i1=1,3
     betaMq2(i1,i1) = betaMq2(i1,i1) + diagonal(5,2)
    End Do

    diagonal(6,2) = 1.6_dp * gauge2(1) * (2._dp*sig(1)/3._dp - S2)           &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp       &
      &  + 3424._dp * gauge2(1) * AbsGM2(1) / 75._dp                         &
      &  + 512._dp * ( gauge2(1) * AbsGM2(3)                                 &
      &              + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)           &
      &                * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)-0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMu2 = -2._dp * ( MuYuaYuYuaYu + YuaYuYuaYuMu )                       &
      & - 4._dp * ( YuMqaYuYuaYu + YuaYuMuYuaYu + YuaYuYuMqaYu )             &
      & - 8._dp * MH(2) * YuaYuYuaYu                                         &
      & - 6._dp * TraceY(4) * ( MuYuaYu + YuaYuMu + 4._dp * MH(2) * YuaYu    &
      &                       + 2._dp * (YuMqaYu + AuaAu)  )                 &
      & - 4._dp * ( Tr3MqaYuYu3aYuMu * YuaYu + AuaAuYuaYu + YuaYuAuaAu       &
      &           + AuaYuYuaAu + YuaAuAuaYu + Tr3aAuAu * YuaYu               &
      &           + Tr3aAuYu * AuaYu + Tr3aYuAu * YuaAu )                    &
      & + hd(1) * ( MuYuaYu + YuaYuMu                                        &
      &           + 2._dp * ( MH(2) * YuaYu + YuMqaYu + AuaAu ) )            &
      & + hd(2) * YuaYu + hc(1) * AuaYu + hc(2) * YuaAu                      &
      & - 2._dp * (MuYuaYdYdaYu + YuaYdYdaYuMu )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YuaYdYdaYu + YuMqaYdYdaYu + YuaYdYdMqaYu &
      &           + YuaYdMdYdaYu + AuaAdYdaYu + YuaYdAdaAu                   &
      &           + AuaYdYdaAu + YuaAdAdaYu )
    Do i1=1,3
     betaMu2(i1,i1) = betaMu2(i1,i1) + diagonal(6,2)
    End Do

   End If 

  !------------------------------------------
  ! beta functions for Higgs mass parameters
  !------------------------------------------
   traceMH1(1) = mH(1) * TraceY(1) + Real( cTrace(YeMlaYe),dp ) &
             & + Real( cTrace(aYeMeYe),dp ) + TraceA(1)
   traceMH1(2) = mH(1) * TraceY(3) + Real( cTrace(YdMqaYd),dp ) &
             & + Real( cTrace(aYdMdYd),dp ) + TraceA(3)
   betamH11 = 6._dp * TraceMH1(2) + 2._dp * TraceMH1(1)      &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1

   traceMH2(2) = mH(2) * TraceY(2) + Real( cTrace(YnuMlaYnu),dp )  &
             & + Real( cTrace(aYnuMrYnu),dp ) + TraceA(2)
   traceMH2(1) = mH(2) * TraceY(4) + Real( cTrace(YuMqaYu),dp )  &
             & + Real( cTrace(aYuMuYu),dp ) + TraceA(4)
   betamH21 = 2._dp * TraceMH2(2) + 6._dp * TraceMH2(1)       &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) + 0.6_dp * S1

   If (TwoLoopRGE) Then
    traceMH1(3) = MH(1) * (6._dp*TraceY2(2) + 2._dp*TraceY2(1) + TraceY2(4) ) &
              & + MH(2) * TraceY2(4)                                          &
              & + 6._dp * ( Real( cTrace(MqaYdYdaYdYd),dp )                   &
              &           + Real( cTrace(aYdMdYdaYdYd),dp )                   &
              &           + Real( cTrace(aAdAdaYdYd),dp )                     &
              &           + Real( cTrace(aAdYdaYdAd),dp )  )                  &
              & + 2._dp * ( Real( cTrace(MlaYeYeaYeYe),dp )                   &
              &           + Real( cTrace(aYeMeYeaYeYe),dp )                   &
              &           + Real( cTrace(aAeAeaYeYe),dp )                     &
              &           + Real( cTrace(aAeYeaYeAe),dp )  )                  &
       & + Real( cTrace(YdMqaYuYuaYd),dp ) + Real( cTrace(YdaYuMuYuaYd),dp ) &
       & + Real( cTrace(YdaYuYuMqaYd),dp ) + Real( cTrace(YuaYdMdYdaYu),dp ) &
       & + Real( cTrace(YdaAuAuaYd),dp ) + Real( cTrace(AdaYuYuaAd),dp )     &
              & + Real( cTrace(AdaAuYuaYd),dp ) + Real( cTrace(YdaYuAuaAd),dp )
    betaMH12 = - 6._dp * traceMH1(3)                                       &
      &   + (32._dp*gauge2(3) - 0.8_dp*gauge2(1) ) * traceMH1(2)           &
      &   + 64._dp * ( AbsGM2(3) * TraceY(3)                               &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(3),dp ) )     &
      &   - 1.6_dp * ( AbsGM2(1) * TraceY(3)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(3),dp ) )     &
      &   + 2.4_dp*gauge2(1) * traceMH1(1)                                 &
      &   + 4.8_dp * ( AbsGM2(1) * TraceY(1)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(1),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)                   &
      &   + 3._dp * gauge2(2) * sig(2)                                     &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)          &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2           &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

    traceMH2(2) = MH(2) * ( 6._dp * TraceY2(3) + TraceY2(4) )                 &
      &  + MH(1) * TraceY2(4)                                                 &
 &  + 6._dp * ( Real(cTrace(MqaYuYuaYuYu),dp) + Real(cTrace(aYuMuYuaYuYu),dp) &
 &            + Real(cTrace(aAuAuaYuYu),dp) + Real(cTrace(aAuYuaYuAu),dp)  )  &
 &  + Real( cTrace(YuMqaYdYdaYu),dp ) + Real( cTrace(YuaYdMdYdaYu),dp )       &
 &  + Real( cTrace(YuaYdYdMqaYu),dp ) + Real( cTrace(YdaYuMuYuaYd),dp )       &
 &  + Real( cTrace(YuaAdAdaYu),dp ) + Real( cTrace(AuaYdYdaAu),dp )           &
      &  + Real( cTrace(AuaAdYdaYu),dp ) + Real( cTrace(YuaYdAdaAu),dp )
    betaMH22 = - 6._dp * traceMH2(2)                                          &
      &   + (32._dp*gauge2(3) + 1.6_dp*gauge2(1) ) * traceMH2(1)              &
      &   + 64._dp * ( AbsGM2(3) * TraceY(4)                                  &
      &          - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(4),dp ) )        &
      &   + 3.2_dp * ( AbsGM2(1) * TraceY(4)                                  &
      &          - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(4),dp ) )        &
      &   + gauge2(1) * ( 0.6_dp * sig(1) + 1.2_dp * S2)                      &
      &   + 3._dp * gauge2(2) * sig(2)                                        &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)             &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2              &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

   End If
!-----------------------------
! beta functions for mu and B
!-----------------------------
   TraceMue(1) = 3._dp * (TraceY(3)+TraceY(4)) + TraceY(1) + TraceY(2)  &
             & - 3._dp * gauge2(2) - 0.6_dp * gauge2(1)
   betaMue1 = mue * TraceMue(1)

   TraceB(1) = 6._dp * (TraceaYA(3)+TraceaYA(4))          &
           & + 2._dp * (TraceaYA(1)+ TraceaYA(2))         &
           & + 6._dp * g2Mi(2) + 1.2_dp * g2Mi(1)
   betaB1 = mue * TraceB(1) + B * TraceMue(1)

   If (TwoLoopRGE) Then
    TraceMue(2) = - 3._dp * ( 3._dp * (TraceY2(2) + TraceY2(3) )          &
      &                     + 2._dp * TraceY2(4) + TraceY2(1) )           &
      &         + (16._dp * gauge2(3) + 0.8_dp * gauge2(1)) * TraceY(4)   &
      &         + (16._dp * gauge2(3) - 0.4_dp * gauge2(1)) * TraceY(3)   &
      &         + 1.2_dp * gauge2(1) * TraceY(1)                          &
      &         + 7.5_dp * gauge2(2)**2                                   &
      &         + 1.8_dp * gauge2(2) * gauge2(1)                          &
      &         + 4.14_dp * gauge2(1)**2
    betaMue2 = mue * TraceMue(2)

    TraceB(1) = cTrace( 3._dp * ( Matmul(AuaYu,YuaYu) + Matmul(AdaYd,YdaYd) ) &
              &       + Matmul(AeaYe,YeaYe) + Matmul(aYuAu,aYdYd)             &
              &       + Matmul(aYdAd,aYuYu) ) 
    TraceB(2) = -12._dp * TraceB(1)                                           &
      &   + (32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)          &
      &   + (32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)          &
      &   + 2.4_dp * gauge2(1) * TraceaYA(1)                                  &
      &   - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)               &
      &   - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)               &
      &   - 2.4_dp * g2Mi(1) * TraceY(1) - 30._dp * gauge2(2)**2 * Mhlf(2)    &
      &   - 3.6_dp * gauge2(2) * gauge2(1) * (Mhlf(1) + Mhlf(2) )             &
      &   - 16.56_dp * gauge2(1)**2 * Mhlf(1) 
    betaB2 = mue * TraceB(2) + B * TraceMue(2)

   End If

  !--------------------------------
  ! neutrino dim. 5 operator
  !--------------------------------
  sumM1 = aYeYe + aYnuYnu
  diagonal(5,1) = 2._dp * TraceY(2) + 6._dp * TraceY(4)   &
              & - 1.2_dp * gauge2(1) - 6._dp * gauge2(2)
  betaMnu1 = Matmul( Transpose(sumM1), Mnu) + Matmul(Mnu, sumM1)  &
          & + diagonal(5,1) * Mnu
  

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1 + oo16pi2 * (Matmul(b_2,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYnu = oo16pi2 * betaYnu1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   Do i1 = 1,3    
    sumI = 0._dp
    Do i2=1,3
     sumI = sumI + b_2(i1,i2) * gauge2(i2) * (Mhlf(i1) + Mhlf(i2) ) 
    End Do
    Do i2=1,4
     sumI = sumI + a_2a(i1,i2) * ( TraceaYA(i2) - Mhlf(i1)*TraceY(i2) )
    End Do
    DMhlf(i1) = oo8pi2 * gauge2(i1) * ( b_1(i1) * Mhlf(i1) + oo16pi2 * sumI)
   End Do
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * ( betaAe1 + oo16pi2 * betaAe2 )
   DAnu = oo16pi2 * betaAnu1
   DAd = oo16pi2 * ( betaAd1 + oo16pi2 * betaAd2 )
   DAu = oo16pi2 * ( betaAu1 + oo16pi2 * betaAu2 )
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * ( betaMe1 + oo16pi2 * betaMe2 )
   DMl = oo16pi2 * ( betaMl1 + oo16pi2 * betaMl2 )
   DMr = oo16pi2 * betaMr1 ! + oo16pi2 * betaMl2 )
   DMd = oo16pi2 * ( betaMd1 + oo16pi2 * betaMd2 )
   DMq = oo16pi2 * ( betaMq1 + oo16pi2 * betaMq2 )
   DMu = oo16pi2 * ( betaMu1 + oo16pi2 * betaMu2 )
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * ( betaMH11 + oo16pi2 * betaMH12 )
   DmH(2) = oo16pi2 * ( betaMH21 + oo16pi2 * betaMH22 )
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * ( betaMue1 + oo16pi2 * betaMue2 )
   DB = oo16pi2 * ( betaB1 + oo16pi2 * betaB2 )

 !----------------------
 ! left neutrino mass
 !----------------------
   DMnu = oo16pi2 * betaMnu1
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYnu = oo16pi2 * betaYnu1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   DMhlf = oo8pi2 * gauge2 * b_1 * Mhlf
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * betaAe1
   DAnu = oo16pi2 * betaAnu1
   DAd = oo16pi2 * betaAd1
   DAu = oo16pi2 * betaAu1
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * betaMe1
   DMl = oo16pi2 * betaMl1
   DMr = oo16pi2 * betaMr1
   DMd = oo16pi2 * betaMd1
   DMq = oo16pi2 * betaMq1
   DMu = oo16pi2 * betaMu1
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * betaMH11
   DmH(2) = oo16pi2 * betaMH21
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * betaMue1
   DB = oo16pi2 * betaB1
  !----------------------
  ! left neutrino mass
  !----------------------
   DMnu = oo16pi2 * betaMnu1
  End If


  !---------------------------------------
  ! to avoid numerical problems in odeint
  !---------------------------------------
  Do i1=1,3
   DMe(i1,i1) = Real(DMe(i1,i1),dp)
   DMr(i1,i1) = Real(DMr(i1,i1),dp)
   DMl(i1,i1) = Real(DMl(i1,i1),dp)
   DMd(i1,i1) = Real(DMd(i1,i1),dp)
   DMu(i1,i1) = Real(DMu(i1,i1),dp)
   DMq(i1,i1) = Real(DMq(i1,i1),dp)
  End Do
  Dmd = 0.5_dp * ( Dmd + Transpose(Conjg(Dmd)) )
  Dme = 0.5_dp * ( Dme + Transpose(Conjg(Dme)) )
  Dml = 0.5_dp * ( Dml + Transpose(Conjg(Dml)) )
  Dmq = 0.5_dp * ( Dmq + Transpose(Conjg(Dmq)) )
  Dmu = 0.5_dp * ( Dmu + Transpose(Conjg(Dmu)) )

  Call Chop(Dmue)
  Call Chop(DB)


  Call ParametersToG3(Dgauge, DYe, DYnu, DYd, DYu, DMhlf, DAe, DAnu, DAd, DAu &
                   &, DMe, DMl, DMr, DMd, DMq, DMu, DMh, Dmue, DB, DMnu, f)

  Iname = Iname - 1

 End Subroutine rge285


 Subroutine rge356(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.01.06: including neutrino dim 5 operator
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2
  Real(dp) :: gauge(3), gauge2(3), TraceY(6), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2), sumI           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3), hc(4)

  Complex(dp) :: Mhlf(3),DMhlf(3)

  Complex(dp) :: Ae(3,3), Ad(3,3), Au(3,3), aAe(3,3), aAd(3,3), aAu(3,3)   &
     &  , DAe(3,3), DAd(3,3), DAu(3,3), aAdAd(3,3), aAeAe(3,3), aAuAu(3,3) &
     &  , aYdAd(3,3), aYeAe(3,3), aYuAu(3,3), TraceaYA(6), betaAd1(3,3)    &
     &  , betaAd2(3,3), betaAe1(3,3), betaAe2(3,3), betaAu1(3,3)           &
     &  , betaAu2(3,3)
  Real(dp) :: TraceA(6)
  Complex(dp) :: aYdYdaYdAd(3,3), aYdAdaYdYd(3,3), TraceAY2(5)               &
     &  , aYeYeaYeAe(3,3), aYeAeaYeYe(3,3), aYuYuaYuAu(3,3), aYuAuaYuYu(3,3) &
     &  , aYuYuaYdAd(3,3), aYuAuaYdYd(3,3), aYdYdaYuAu(3,3), aYdAdaYuYu(3,3)
   
  Complex(dp) :: Me(3,3), Ml(3,3), Md(3,3), Mq(3,3), Mu(3,3), DMe(3,3)        &
     & , DMl(3,3), DMd(3,3), DMq(3,3), DMu(3,3), YdaYd(3,3), YeaYe(3,3)       &
     & , YuaYu(3,3), MdYdaYd(3,3), MeYeaYe(3,3), MuYuaYu(3,3), YdaYdMd(3,3)   &
     & , YeaYeMe(3,3), YuaYuMu(3,3), YdMqaYd(3,3), YeMlaYe(3,3)               &
     & ,  YuMqaYu(3,3), AdaAd(3,3), AeaAe(3,3), AuaAu(3,3), betaMd1(3,3)      &
     & , betaMd2(3,3), betaMe1(3,3), betaMe2(3,3), betaMl1(3,3), betaMl2(3,3) &
     & , betaMq1(3,3), betaMq2(3,3), betaMu1(3,3), betaMu2(3,3), MqaYdYd(3,3) &
     & , MqaYuYu(3,3), aYdYdMq(3,3), aYuYuMq(3,3), aYeYeMl(3,3), MlaYeYe(3,3) &
     & , aYeMeYe(3,3), aYdMdYd(3,3), aYuMuYu(3,3)                             &
     & , YdaYdYdaYd(3,3), YeaYeYeaYe(3,3), YuaYuYuaYu(3,3), MeYeaYeYeaYe(3,3) &
     & , YeaYeYeaYeMe(3,3), YeaYeMeYeaYe(3,3), YeMlaYeYeaYe(3,3), AeaYe(3,3)  &
     & , YeaYeYeMlaYe(3,3), AeaAeYeaYe(3,3), YeaYeAeaAe(3,3), YeaAe(3,3)      &
     & , AeaYeYeaAe(3,3), YeaAeAeaYe(3,3), Tr3aAdYdaAeYe

  Complex(dp) :: Tr3aYdAdaYeAe, AdaYd(3,3), YdaAd(3,3), MlaYeYeaYeYe(3,3)     &
     & , aYeYeaYeYeMl(3,3), aYeYeMlaYeYe(3,3), aYeYeaYeMeYe(3,3)              &
     & , aYeMeYeaYeYe(3,3), aAdYd(3,3),aAeYe(3,3), aAeAeaYeYe(3,3)            &
     & , aYeYeaAeAe(3,3), aAeYeaYeAe(3,3), aYeAeaAeYe(3,3), MdYdaYdYdaYd(3,3) &
     & , YdaYdYdaYdMd(3,3), YdMqaYdYdaYd(3,3),YdaYdMdYdaYd(3,3)               &
     & , YdaYdYdMqaYd(3,3), AdaAdYdaYd(3,3), YdaYDAdaAd(3,3), AdaYdYdaAd(3,3) &
     & , YdaAdAdaYd(3,3)
  Complex(dp) :: MdYdaYuYuaYd(3,3), YdaYuYuaYdMd(3,3), YdMqaYuYuaYd(3,3)      &
     & , YdaYuYuMqaYd(3,3), YdaYuMuYuaYd(3,3), AdaAuYuaYd(3,3)                &
     & , YdaYuAuaAd(3,3), AdaYuYuaAd(3,3), YdaAuAuaYd(3,3), YdaYuYuaYd(3,3)   &
     & , Tr3aYuAu, Tr3aAuYu, YuaAu(3,3)                                       &
     & , MqaYdYdaYdYd(3,3), aYdYdaYdYDMq(3,3), aYdMdYdaYdYd(3,3)              &
     & , aYdYdMqaYdYd(3,3), aYdYdaYdMdYd(3,3), aAdAdaYdYd(3,3)                &
     & , aYdYDaAdAd(3,3), aAdYdaYdAd(3,3), aYdAdaAdYd(3,3), MqaYuYuaYuYu(3,3) &
     & , aYuYuaYuYUMq(3,3), aYuMuYuaYuYu(3,3), aYuYuMqaYuYu(3,3)              &
     & , aYuYuaYuMuYu(3,3), aAuAuaYuYu(3,3), aYuYUaAuAu(3,3), aAuYuaYuAu(3,3) &
     & , aYuAuaAuYu(3,3), aAuYu(3,3), AuaYu(3,3), YuaYdYdaYu(3,3)             &
     & , AuaYdYdaAu(3,3), YuaAdAdaYu(3,3), AuaAdYdaYu(3,3), YuaYdAdaAu(3,3)   &
     & , YuMqaYuYuaYu(3,3), YuaYuYuMqaYu(3,3), MuYuaYuYuaYu(3,3)              &
     & , YuaYuYuaYuMu(3,3), YuaYuMuYuaYu(3,3), AuaAuYuaYu(3,3)                &
     & , YuaYuAuaAu(3,3), AuaYuYuaAu(3,3), YuaAuAuaYu(3,3), MuYuaYdYdaYu(3,3) &
     & , YuaYdYdaYuMu(3,3), YuMqaYdYdaYu(3,3), YuaYdYdMqaYu(3,3)              &
     & , YuaYdMdYdaYu(3,3)

  Real(dp) :: S1, S2, sig(3), Tr3aYdYdaYeYe, Tr3aAdAdaAeAe, AbsGM2(3)         &
     & , Tr3MqaYdYd3aYDMdYd , Tr3MqaYuYu3aYuMu, Tr3aAuAu

  Real(dp) :: Mh(2), DMh(2), TraceMH1(3), TraceMH2(2), betaMH11        &
     &  , betaMH12, betaMH21, betaMH22, q

  Complex(dp) :: mue, B, Dmue, DB, TraceMue(2), TraceB(2), betaMue1, betaMue2 &
     & , betaB1, betaB2

  Complex(dp) :: Mnu(3,3), DMnu(3,3), betamnu1(3,3)

  Complex(dp), Dimension(3,3) :: YT, aYT, aYTYT, betaYT1, sumT1, AT, aAT     &
      & , aATAT, betaAT1, aYTAT, MlaYTYT, aYTMlYT, aYTYTMl, DYT, DAT
  Complex(dp), Dimension(3,3) :: YZ, aYZ, aYZYZ, YZaYZ, betaYZ1, sumZ1, DYZ  &
    & , AS, aAS, ASaAS, betaAS1, aYSAS, DAS                                  &
    & , AZ, aAZ, aAZAZ, AZaAZ, betaAZ1, aYZAZ, DAZ
  Complex(dp), Dimension(3,3) :: YS, aYS, aYSYS, YSaYS, betaYS1, sumS1, DYS &
    & , AZaYZ, ASaYS, aYSMdYS, MlaYZYZ, aYZYZMl, aYZMdYZ, MdYZaYZ, YZaYZMd  &
    & , MdYSaYS, YSaYSMd, YZMlaYZ, YSMdaYS
  Complex(dp) :: g2Mi(3)
  Complex(dp) :: lam1, lam2, Dlam1, Dlam2, betalam11, betalam21, lam1Alam1 &
      & , lam2Alam2, Alam1, Alam2, betaAlam11, betaAlam21, DAlam1, DAlam2
  Real(dp) :: MT(2), lam12, lam22, betaMT2(2), Alam12, Alam22, b_1a(3), DMT2(2) &
      & , MS(2), betaMS2(2), DMS2(2), MZ(2), betaMZ2(2), DMZ2(2), MS15, MZ15    &
      & , MT15, DMS15, DMZ15, DMT15, betaMS15, betaMZ15, betaMT15, b_2a(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge356'

  OnlyDiagonal = .Not.GenerationMixing
  q = t

  Call GToParameters5(gy, gauge, Ye, YT, Yd, Yu, YZ, YS, lam1, lam2, Mhlf    &
                & , Ae, AT, Ad, Au, AZ, AS, Alam1, Alam2, Me, Ml, Md, Mq, Mu &
                & , Mh, mT, mZ, mS, MT15, MZ15, MS15, mue, B, Mnu)

  gauge2 = gauge**2
  AbsGM2 = gauge2 * Abs( Mhlf )**2
!-----------------
! beta functions
!-----------------
  lam12 = Abs(lam1)**2
  lam22 = Abs(lam2)**2

  Call Adjungate(Yd,aYd)
  Call Adjungate(YT,aYT)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)
  Call Adjungate(YZ,aYZ)
  Call Adjungate(YS,aYS)

  aYeYe = Matmul(aYe,Ye)
  aYTYT = Matmul(aYT,YT)
  aYdYd = Matmul(aYd,Yd)
  aYuYu = Matmul(aYu,Yu)
  aYZYZ = Matmul(aYZ,YZ)
  aYSYS = Matmul(aYS,YS)

  YdaYd = Matmul(Yd,aYd)
  YZaYZ = Matmul(YZ,aYZ)
  YSaYS = Matmul(YS,aYS)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYTYT),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )
  TraceY(5) = Real( cTrace(aYZYZ),dp )
  TraceY(6) = Real( cTrace(aYSYS),dp )

  diagonal(1,1) = 3._dp * (TraceY(3) + lam12) + TraceY(1)    &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * (aYeYe + aYTYT + aYZYZ)
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
  sumT1 = aYeYe + 6._dp * aYTYT + 3._dp * aYZYZ
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaYT1 = Matmul(YT,sumT1)                 &
        & + Matmul(Transpose(aYeYe+3._dp * aYZYZ),YT)

  diagonal(3,1) = 3._dp * (TraceY(3)  + lam12) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)  &
        & + 2._dp * Matmul(YZaYZ + 2._dp * YSaYS, Yd)

  diagonal(4,1) = 3._dp * (TraceY(4) + lam22)                &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  diagonal(5,1) = TraceY(5) + c1_1(2,1) * gauge2(1) &
            &   + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumZ1 = aYeYe + 3._dp * aYTYT + 5._dp * aYZYZ
  Do i1=1,3
   sumZ1(i1,i1) = sumZ1(i1,i1) + diagonal(5,1)
  End Do

  betaYZ1 = Matmul(YZ,sumZ1)                 &
        & + 2._dp * Matmul(YdaYd+2._dp * aYSYS,YZ)

  diagonal(6,1) = TraceY(6) - 0.8_dp * gauge2(1) - 12._dp * gauge2(3)
  sumS1 = 2._dp * Transpose(YdaYd + YZaYZ) + 8._dp * aYSYS
  Do i1=1,3
   sumS1(i1,i1) = sumS1(i1,i1) + diagonal(6,1)
  End Do

  betaYS1 = Matmul(YS,sumS1)                 &
        & + 2._dp * Matmul(YdaYd + YZaYZ,YS)


  betalam11 = lam1 * (7._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
            &        + 6._dp * TraceY(3)                            &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  betalam21 = lam2 * (7._dp * lam22  + 6._dp * TraceY(4)     &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
   aYdYdaYuYu = Matmul(aYdYd,aYuYu)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ae,aAe)
  Call Adjungate(AT,aAT)
  Call Adjungate(Ad,aAd)
  Call Adjungate(Au,aAu)
  Call Adjungate(AZ,aAZ)
  Call Adjungate(AS,aAS)

  aAdAd = Matmul(aAd,Ad)
  aATAT = Matmul(aAT,AT)
  aAeAe = Matmul(aAe,Ae)
  aAuAu = Matmul(aAu,Au)
  aAZAZ = Matmul(aAZ,AZ)

  aYdAd = Matmul(aYd,Ad)
  aYTAT = Matmul(aYT,AT)
  aYeAe = Matmul(aYe,Ae)
  aYuAu = Matmul(aYu,Au)
  aYZAZ = Matmul(aYZ,AZ)
  aYSAS = Matmul(aYS,AS)

  AeaYe = Matmul(Ad,aYd)
  AdaYd = Matmul(Ad,aYd)
  AZaYZ = Matmul(AZ,aYZ)
  ASaYS = Matmul(AS,aYS)

  ASaAS = Matmul(AS,aAS)

  TraceA(1) = Real( cTrace(aAeAe),dp )
  TraceA(2) = Real( cTrace(aATAT),dp )
  TraceA(3) = Real( cTrace(aAdAd),dp )
  TraceA(4) = Real( cTrace(aAuAu),dp )
  TraceA(5) = Real( cTrace(aAZAZ),dp )
  TraceA(6) = Real( cTrace(ASaAS),dp )

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYTAT) 
  TraceaYA(3) = cTrace(aYdAd) 
  TraceaYA(4) = cTrace(aYuAu) 
  TraceaYA(5) = cTrace(aYZAZ) 
  TraceaYA(6) = cTrace(aYSAS) 
  lam1Alam1 = Conjg(lam1) * Alam1
  lam2Alam2 = Conjg(lam2) * Alam2
  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = Matmul(Ae,sume1)
  
  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1)  &
                &         + 3._dp * lam1Alam1                  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    ) 
  sume1 = 4._dp * aYeAe + 6._dp * (aYTAT + aYZAZ)
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + Matmul(Ye,sume1)

  !--------------
  ! A_T
  !--------------
  diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
  sumT1 = 2._dp * aYeYe + 9._dp * aYTYT + 3._dp * aYZYZ
  betaAT1 = Matmul(AT,sumT1)
  betaAT1 = Transpose(betaAT1)
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaAT1 = betaAT1 + Matmul(AT,sumT1)

  diagonal(2,1) = 2._dp * ( TraceaYA(2) + lam1Alam1    &
      &         + 1.8_dp * g2Mi(1) + 7._dp * g2Mi(2) )
  betaAT1 = betaAT1 + diagonal(2,1) * YT                         &
        & + 6._dp * ( Matmul(YT, aYZAZ)            &
        &           + Matmul(Transpose(aYZAZ), YT) )

  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = Matmul(Ad,sumd1)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1) &
                &         + 3._dp * lam1Alam1                 &
                &         - c1_1(2,1) * g2Mi(1)   &
                &         - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do
  betaAd1 = betaAd1 + Matmul(Yd,sumd1)            &
        & + 2._dp * Matmul(2._dp*YSaYS+YZaYZ, Ad) &
        & + 4._dp * Matmul(2._dp*ASaYS+AZaYZ, Yd)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = Matmul(Au,sumu1)
  
  diagonal(4,1) = 2._dp * ( 3._dp * TraceaYA(4) + 3._dp * lam2Alam2   &
                &         - c1_1(3,1) * g2Mi(1) - c1_1(3,2) * g2Mi(2) &
                &         - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do
  betaAu1 = betaAu1 + Matmul(Yu,sumu1)

  !--------------
  ! A_Z
  !--------------
  sumZ1 = sumZ1 + 3._dp * aYZYZ
  betaAZ1 = Matmul(AZ,sumZ1) + 6._dp * Matmul(YZ, aYTAT)      &
    & + 2._dp * Matmul(2._dp*YSaYs + YdaYd + 2._dp*YZaYZ, AZ)

  diagonal(5,1) = - 2._dp * ( c1_1(2,1) * g2Mi(1) + c1_1(2,2) * g2Mi(2) &
                &           + c1_1(2,3) * g2Mi(3)  ) + 2._dp * TraceaYA(5)

  betaAZ1 = betaAZ1 + diagonal(5,1) * YZ                                                &
        & +  2._dp * Matmul(YZ, AeaYe)                              &
        & + 4._dp * Matmul(AdaYd + 2._dp * ASaYS, YZ)

  !--------------
  ! A_S
  !--------------
  diagonal(6,1) = TraceY(6) - 0.8_dp * gauge2(1) - 12._dp * gauge2(3)
  sumS1 = 2._dp * Transpose(YdaYd + YZaYZ) + 12._dp * aYSYS
  betaAS1 = Matmul(AS,sumS1)
  betaAS1 = Transpose(betaAS1)
  Do i1=1,3
   sumS1(i1,i1) = sumS1(i1,i1) + diagonal(6,1)
  End Do

  betaAS1 = betaAS1 + Matmul(AS,sumS1)

  diagonal(6,1) = 2._dp * ( TraceaYA(6) + 0.8_dp * g2Mi(1) + 12._dp * g2Mi(3) )
  sumS1 = 4._dp * (AdaYd + AZaYZ) 
  betaAS1 = betaAS1 + diagonal(6,1) * YS               &
        & + Matmul(sumS1, YS) + Matmul(YS,Transpose(sumS1)) 

  !--------------
  ! A_1
  !--------------
  betaAlam11 = Alam1 * (21._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
          &           + 6._dp * TraceY(3)                              &
          &           - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) )       &
          & + 2._dp * lam1 * ( TraceaYA(2) + 2._dp * TraceaYA(1)       &
          &                  + 6._dp * TraceaYA(3)                     &
          &                  + 1.8_dp * g2Mi(1)  + 7._dp * g2Mi(2) ) 

  !--------------
  ! A_2
  !--------------
  betaAlam21 = Alam2 * (21._dp * lam22 + 6._dp * TraceY(4)          &
          &           - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) )    &
          & + 2._dp * lam2 * ( 6._dp * TraceaYA(4)                  &
          &                  + 1.8_dp * g2Mi(1) + 7._dp * g2Mi(2) ) 

  If (TwoLoopRGE) Then
   aYdYdaYdAd = Matmul(aYdYd,aYdAd)
   aYdAdaYdYd = Matmul(aYdAd,aYdYd)
   aYeYeaYeAe = Matmul(aYeYe,aYeAe)
   aYeAeaYeYe = Matmul(aYeAe,aYeYe)
   aYuYuaYuAu = Matmul(aYuYu,aYuAu)
   aYuAuaYuYu = Matmul(aYuAu,aYuYu)
   aYuAuaYdYd = Matmul(aYuAu,aYdYd)
   aYuYuaYdAd = Matmul(aYuYu,aYdAd)
   aYdAdaYuYu = Matmul(aYdAd,aYuYu)
   aYdYdaYuAu = Matmul(aYdYd,aYuAu)
   TraceAY2(1) = cTrace(aYeYeaYeAe)
   TraceAY2(2) = cTrace(aYdYdaYdAd)
   TraceAY2(3) = cTrace(aYuYuaYuAu)
   TraceAY2(4) = cTrace(aYuYuaYdAd)
   TraceAY2(5) = cTrace(aYdYdaYuAu)

  !--------------
  ! A_e
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = Matmul(Ae,sume2)
    
   diagonal(1,2) = -6._dp * ( 6._dp * TraceAY2(2) + TraceAY2(4)       &
     &                      + TraceAY2(5) + 2._dp * TraceAY2(1)  )    &
     &  + ( 32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)   &
     &  - ( 32._dp * g2Mi(3)                              &
     &    - 0.8_dp * g2Mi(1) ) * TraceY(3)                &
     &  + 2.4_dp * gauge2(1) * TraceaYA(1)                            &
     &  - 2.4_dp * g2Mi(1) * TraceY(1)                    &
     &  - ( 30._dp * g2Mi(2)                              &
     &    + 3.6_dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)      &
     &  - 54._dp * gauge2(1)**2 * Mhlf(1)
   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)     &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1) &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2)
   sume2 = - 6._dp * aYeYeaYeAe - 8._dp * aYeAeaYeYe  &
         & - hd(1) * aYeAe - hc(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
   betaAe2 = betaAe2 + Matmul(Ye,sume2)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = Matmul(Ad,sumd2)
    
   diagonal(3,2) = diagonal(1,2)                                   &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                 &
     &              - gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp      &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)    &
     &  + 1.6_dp * gauge2(1) * gauge2(2) * (Mhlf(1)+Mhlf(2))       & 
     &  + 1.856e3_dp * gauge2(1)**2 * Mhlf(1) / 4.5e1_dp

   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)  &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hd(2) = 6._dp * TraceY(4) - 1.6_dp * gauge2(1) 
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1)                 &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2) + 1.6_dp * gauge2(1) *  Mhlf(1)
   hc(2) = 6._dp * TraceaYA(4) + 1.6_dp * gauge2(1) *  Mhlf(1)
   sumd2 = - 6._dp * aYdYdaYdAd - 8._dp * aYdAdaYdYd                  &
       &   - 4._dp * ( aYuAuaYuYu + aYuYuaYuAu + aYuAuaYdYd )         &
       &   - 2._dp * aYuYuaYdAd - hd(1) * aYdAd - hc(1) * aYdYd       &
       &  - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(3,2)
   End Do
   betaAd2 = betaAd2 + Matmul(Yd,sumd2)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(4) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = Matmul(Au,sumu2)
    
   diagonal(4,2) =  -6._dp * ( 6._dp * TraceAY2(3) + TraceAY2(4)        &
     &                        + TraceAY2(5)  )                          &
     &  + ( 32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)     &
     &  - ( 32._dp * g2Mi(3)                                &
     &    + 1.6_dp * g2Mi(1) ) * TraceY(4)                  &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                      &
     &              - 3.4_dp * gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp  &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)         &
     &  - ( 30._dp * g2Mi(2)                                &
     &    + 2._dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)         &
     &  - 5486._dp * gauge2(1)**2 * Mhlf(1) / 225._dp
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1) - 0.8_dp * gauge2(1)
   hc(1) = 6._dp * TraceaYA(3) + 2._dp * TraceaYA(1)   &
       &  + 0.8_dp * gauge2(1) *  Mhlf(1)
   hd(2) = 12._dp * TraceY(4) - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(2) = 18._dp * TraceaYA(4) + 1.2e1_dp * g2Mi(2) & 
       & + 0.8_dp * g2Mi(1)
   sumu2 = - 6._dp * aYuYuaYuAu - 8._dp * aYuAuaYuYu                     &
       &   - 4._dp * ( aYdAdaYdYd + aYdYdaYdAd + aYdAdaYuYu )            &
       &   - 2._dp * aYdYdaYuAu - hd(1) * aYdAd - hc(1) * aYdYd          &
       &   - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(4,2)
   End Do
   betaAu2 = betaAu2 + Matmul(Yu,sumu2)

  End If 
!----------------------------------------------
! beta functions for Sfermion mass parameters
!----------------------------------------------
   S1 = mH(2) - mH(1) + 3._dp * (mT(1) - mT(2) ) &
    & + mZ(1) - mZ(2) + 4._dp * (mS(2) - mS(1))
   Do i1=1,3
    S1 = S1 + Me(i1,i1) - Ml(i1,i1) &
       &    + Md(i1,i1) + Mq(i1,i1) - 2._dp * Mu(i1,i1)
   End Do

   S1 = S1 * gauge2(1)

   YdaYd = Matmul(Yd,aYd)
   YeaYe = Matmul(Ye,aYe)
   YuaYu = Matmul(Yu,aYu)

   MeYeaYe = Matmul(Me,YeaYe)
   MlaYeYe = Matmul(Ml,aYeYe)
   MlaYTYT = Matmul(Ml,aYTYT)
   MlaYZYZ = Matmul(Ml,aYZYZ)

   MdYdaYd = Matmul(Md,YdaYd)
   MqaYdYd = Matmul(Mq,aYdYd)
   MqaYuYu = Matmul(Mq,aYuYu)
   MuYuaYu = Matmul(Mu,YuaYu)
   MdYZaYZ = Matmul(Md,YZaYZ)
   MdYSaYS = Matmul(Md,YSaYS)

   YeaYeMe = Matmul(YeaYe,Me)
   aYeYeMl = Matmul(aYeYe,Ml)
   aYTYTMl = Matmul(aYTYT,Ml)
   aYZYZMl = Matmul(aYZYZ,Ml)

   YdaYdMd = Matmul(YdaYd,Md)
   aYdYdMq = Matmul(aYdYd,Mq)
   aYuYuMq = Matmul(aYuYu,Mq)
   YuaYuMu = Matmul(YuaYu,Mu)
   YZaYZMd = Matmul(YZaYZ,Md)
   YSaYSMd = Matmul(YSaYS,Md)

   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   aYTMlYT = MatMul3(aYT,Transpose(Ml),YT,OnlyDiagonal)
   YZMlaYZ = MatMul3(YZ,Ml,YZ,OnlyDiagonal)
   YSMdaYS = MatMul3(YS,Transpose(Md),aYS,OnlyDiagonal)
   aYSMdYS = MatMul3(aYS,Md,YS,OnlyDiagonal)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)
   aYZMdYZ = MatMul3(aYZ,Md,YZ,OnlyDiagonal)

   AeaAe = Matmul(Ae,aAe)
   AdaAd = Matmul(Ad,aAd)
   AuaAu = Matmul(Au,aAu)
   AZaAZ = Matmul(AZ,aAZ)
   Alam12 = Abs(Alam1)**2
   Alam22 = Abs(Alam2)**2

   diagonal(1,1) = - 4.8_dp * AbsGM2(1) + 1.2_dp * S1
   betaMe1 = 2._dp * (MeYeaYe + YeaYeMe)             &
         & + 4._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe )
   Do i1=1,3
    betaMe1(i1,i1) = betaMe1(i1,i1) + diagonal(1,1)
   End Do

   diagonal(3,1) = - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1
   betaMl1 = MlaYeYe + aYeYeMl + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe )  &
         & + 3._dp * ( MlaYTYT + MlaYZYZ + aYTYTMl + aYZYZMl )                &
         & + 6._dp * ( aYTMlYT + aATAT + MT(1) * aYTYT)                       &
         & + 6._dp * ( aYZMdYZ + aAZAZ + MZ(1) * aYZYZ)

   Do i1=1,3
    betaMl1(i1,i1) = betaMl1(i1,i1) + diagonal(3,1)
   End Do

   diagonal(4,1) = - ( 32._dp * AbsGM2(3) + 1.6_dp * AbsGM2(1) ) / 3._dp &
               & + 0.4_dp * S1
   betaMd1 = 2._dp * (MdYdaYd + YdaYdMd + MdYZaYZ + YZaYZMd )   &
         & + 4._dp * (MdYSaYS + YSaYSMd)                        &
         & + 4._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd )        &
         & + 4._dp * ( mZ(1) * YZaYZ + YZMlaYZ + AZaAZ )        &
         & + 8._dp * ( mS(1) * YSaYS + YSMdaYS + ASaAS )
   Do i1=1,3
    betaMd1(i1,i1) = betaMd1(i1,i1) + diagonal(4,1)
   End Do

   diagonal(5,1) = - ( 32._dp * AbsGM2(3) + 0.4_dp * AbsGM2(1) ) / 3._dp &
               & - 6._dp * AbsGM2(2) + 0.2_dp * S1
   betaMq1 = MqaYuYu + aYuYuMq + MqaYdYd + aYdYdMq             &
         & + 2._dp * ( mH(2) * aYuYu + mH(1) * aYdYd + aYuMuYu     &
         &           + aYdMdYd + aAuAu + aAdAd )
   Do i1=1,3
    betaMq1(i1,i1) = betaMq1(i1,i1) + diagonal(5,1)
   End Do

   diagonal(6,1) = - ( 32._dp * AbsGM2(3) + 6.4_dp * AbsGM2(1) ) / 3._dp &
               &   - 0.8_dp * S1
   betaMu1 = 2._dp * (MuYuaYu + YuaYuMu)             &
         & + 4._dp * ( mH(2) * YuaYu + YuMqaYu + AuaAu )
   Do i1=1,3
    betaMu1(i1,i1) = betaMu1(i1,i1) + diagonal(6,1)
   End Do

   If (TwoLoopRGE) Then
    YdaYdYdaYd = MatSquare(YdaYd,OnlyDiagonal)
    YeaYeYeaYe = MatSquare(YeaYe,OnlyDiagonal)
    YuaYuYuaYu = MatSquare(YuaYu,OnlyDiagonal)

    AuaYu = Matmul(Au,aYu)

    aAdYd = Matmul(aAd,Yd)
    aAeYe = Matmul(aAe,Ye)
    aAuYu = Matmul(aAu,Yu)

    YdaAd = Matmul(Yd,aAd)
    YeaAe = Matmul(Ye,aAe)
    YuaAu = Matmul(Yu,aAu)

    YdaYuYuaYd = MatMul3(Yd,aYuYu,aYd,OnlyDiagonal)
    AdaYuYuaAd = MatMul3(Ad,aYuYu,aAd,OnlyDiagonal)
    YdaAuAuaYd = MatMul3(Yd,aAuAu,aYd,OnlyDiagonal)
    AdaAuYuaYd = MatMul4(Ad,aAu,Yu,aYd,OnlyDiagonal)
    YdaYuAuaAd = MatMul3(Yd,aYuAu,aAd,OnlyDiagonal)

    YuaYdYdaYu = MatMul3(Yu,aYdYd,aYu,OnlyDiagonal)
    AuaYdYdaAu = MatMul3(Au,aYdYd,aAu,OnlyDiagonal)
    YuaAdAdaYu = MatMul3(Yu,aAdAd,aYu,OnlyDiagonal)
    AuaAdYdaYu = MatMul4(Au,aAd,Yd,aYu,OnlyDiagonal)
    YuaYdAdaAu = MatMul3(Yu,aYdAd,aAu,OnlyDiagonal)

    MdYdaYuYuaYd = Matmul(Md,YdaYuYuaYd)
    Call Adjungate(MdYdaYuYuaYd, YdaYuYuaYdMd)
    YdMqaYuYuaYd = MatMul3(Yd,MqaYuYu,aYd,OnlyDiagonal)
    Call Adjungate(YdMqaYuYuaYd, YdaYuYuMqaYd)
    YdaYuMuYuaYd = MatMul3(Yd,aYuMuYu,aYd,OnlyDiagonal)

    MuYuaYdYdaYu = Matmul(Mu,YuaYdYdaYu)
    Call Adjungate(MuYuaYdYdaYu, YuaYdYdaYuMu)
    YuMqaYdYdaYu = MatMul3(Yu,MqaYdYd,aYu,OnlyDiagonal)
    Call Adjungate(YuMqaYdYdaYu, YuaYdYdMqaYu)
    YuaYdMdYdaYu = MatMul3(Yu,aYdMdYd,aYu,OnlyDiagonal)

    MeYeaYeYeaYe = Matmul(MeYeaYe,YeaYe)
    Call Adjungate(MeYeaYeYeaYe,YeaYeYeaYeMe)
    aYeMeYeaYeYe = Matmul(aYeMeYe,aYeYe)
    Call Adjungate(aYeMeYeaYeYe,aYeYeaYeMeYe)
    YeaYeMeYeaYe = Matmul(YeaYeMe,YeaYe)

    MlaYeYeaYeYe = Matmul(MlaYeYe,aYeYe)
    Call Adjungate(MlaYeYeaYeYe, aYeYeaYeYeMl)
    YeMlaYeYeaYe = Matmul(YeMlaYe,YeaYe)
    Call Adjungate(YeMlaYeYeaYe, YeaYeYeMlaYe)
    aYeYeMlaYeYe = Matmul(aYeYeMl,aYeYe)

    MdYdaYdYdaYd = Matmul(MdYdaYd,YdaYd)
    Call Adjungate(MdYdaYdYdaYd, YdaYdYdaYdMd)
    aYdMdYdaYdYd = Matmul(aYdMdYd,aYdYd)
    Call Adjungate(aYdMdYdaYdYd, aYdYdaYdMdYd)
    YdaYdMdYdaYd = Matmul(YdaYdMd,YdaYd)

    MqaYdYdaYdYd = Matmul(MqaYdYd,aYdYd)
    Call Adjungate(MqaYdYdaYdYd, aYdYdaYdYdMq)
    YdMqaYdYdaYd = Matmul(YdMqaYd,YdaYd)
    Call Adjungate(YdMqaYdYdaYd, YdaYdYdMqaYd)
    aYdYdMqaYdYd = Matmul(aYdYdMq,aYdYd)

    MqaYuYuaYuYu = Matmul(MqaYuYu,aYuYu)
    Call Adjungate(MqaYuYuaYuYu, aYuYuaYuYuMq)
    YuMqaYuYuaYu = Matmul(YuMqaYu,YuaYu)
    Call Adjungate(YuMqaYuYuaYu, YuaYuYuMqaYu)
    aYuYuMqaYuYu = Matmul(aYuYuMq,aYuYu)

    MuYuaYuYuaYu = Matmul(MuYuaYu,YuaYu)
    Call Adjungate(MuYuaYuYuaYu, YuaYuYuaYuMu)
    aYuMuYuaYuYu = Matmul(aYuMuYu,aYuYu)
    Call Adjungate(aYuMuYuaYuYu, aYuYuaYuMuYu)
    YuaYuMuYuaYu = Matmul(YuaYuMu,YuaYu)

    AdaAdYdaYd = Matmul(AdaAd,YdaYd)
    Call Adjungate(AdaAdYdaYd, YdaYdAdaAd)
    AdaYdYdaAd = Matmul(AdaYd,YdaAd)
    YdaAdAdaYd = Matmul(YdaAd,AdaYd)

    aAdAdaYdYd = Matmul(aAdAd,aYdYd)
    Call Adjungate(aAdAdaYdYd, aYdYdaAdAd)
    aAdYdaYdAd = Matmul(aAdYd,aYdAd)
    aYdAdaAdYd = Matmul(aYdAd,aAdYd)

    AeaAeYeaYe = Matmul(AeaAe,YeaYe)
    Call Adjungate(AeaAeYeaYe, YeaYeAeaAe)
    AeaYeYeaAe = Matmul(AeaYe,YeaAe)
    YeaAeAeaYe = Matmul(YeaAe,AeaYe)

    aAeAeaYeYe = Matmul(aAeAe,aYeYe)
    Call Adjungate(aAeAeaYeYe, aYeYeaAeAe)
    aAeYeaYeAe = Matmul(aAeYe,aYeAe)
    aYeAeaAeYe = Matmul(aYeAe,aAeYe)

    AuaAuYuaYu = Matmul(AuaAu,YuaYu)
    Call Adjungate(AuaAuYuaYu, YuaYuAuaAu)
    AuaYuYuaAu = Matmul(AuaYu,YuaAu)
    YuaAuAuaYu = Matmul(YuaAu,AuaYu)

    aAuAuaYuYu = Matmul(aAuAu,aYuYu)
    Call Adjungate(aAuAuaYuYu, aYuYuaAuAu)
    aAuYuaYuAu = Matmul(aAuYu,aYuAu)
    aYuAuaAuYu = Matmul(aYuAu,aAuYu)

    S2 = (1.5_dp * gauge2(2) + 0.3_dp * gauge2(1) )            &
     &      * (MH(2) - MH(1) - Real(cTrace(ML),dp) )              &
     & + ( (8._dp * gauge2(3) + 0.1_dp*gauge2(1)) / 3._dp      &
     &   + 1.5_dp * gauge2(2) ) * Real( cTrace(Mq),dp   )         &
     & - (16._dp * gauge2(3) + 3.2_dp*gauge2(1) )              &
     &    * Real(cTrace(Mu),dp) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md),dp) / 3._dp                           &
     & + 1.2_dp*gauge2(1) * Real(cTrace(Me),dp)                   &
     & - 3._dp * (MH(2)*TraceY(4) - MH(1) * TraceY(3) )        &
     & + mH(1) * TraceY(1)
    
    Do i1=1,3
     S2 = S2 - YuMqaYu(i1,i1) + 4._dp * aYuMuYu(i1,i1)   &
        &    - YdMqaYd(i1,i1) - 2._dp * aYdMdYd(i1,i1)   &
        &    + YeMlaYe(i1,i1) - 2._dp * aYeMeYe(i1,i1)
    End Do

    sig(1) = 3._dp * (MH(1) + MH(2) + Real(cTrace(Ml),dp) )          &
         & + Real(cTrace(Mq),dp) + 8._dp * Real(cTrace(Mu),dp)          &
         & + 2._dp * Real(cTrace(Md),dp) + 6._dp * Real(cTrace(Me),dp) 
    sig(1) = 0.2_dp * gauge2(1) * sig(1)
    sig(2) = gauge2(2) * ( MH(1) + MH(2) + Real( cTrace(Ml),dp )    &
           &             + 3._dp * Real( cTrace(Mq),dp ) )
    sig(3) = gauge2(3) * ( 2._dp * Real( cTrace(Mq),dp )            &
           &             + Real(cTrace(Mu),dp) + Real(cTrace(Md),dp) )

    Tr3aYdYdaYeYe = 3._dp * TraceY(3) + TraceY(1)
    Tr3MqaYdYd3aYDMdYd = Real(cTrace(MlaYeYe),dp) + Real(cTrace(aYeMeYe),dp)  &
            & + 3._dp * ( Real(cTrace(MqaYdYd),dp) + Real(cTrace(aYdMdYd),dp) )
    Tr3aAdAdaAeAe = 3._dp * TraceA(3) + TraceA(1)
    Tr3aYdAdaYeAe = 3._dp * TraceaYA(3) + TraceaYA(1)
    Tr3aAdYdaAeYe = Conjg( Tr3aYdAdaYeAe )
    Tr3MqaYuYu3aYuMu = 3._dp * (Real(cTrace(MqaYuYu),dp) &
                               + Real(cTrace(aYuMuYu),dp) )
    Tr3aAuAu = 3._dp * TraceA(4)
    Tr3aYuAu = 3._dp * TraceaYA(4)
    Tr3aAuYu = Conjg( Tr3aYuAu )

    diagonal(1,2) = 2.4_dp * gauge2(1) * (S2 + sig(1) )   &
                & + 112.32_dp * gauge2(1) * AbsGM2(1)
    hd(1) = 6._dp*gauge2(2)-1.2_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 4.8_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 2.4_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMe2 = -2._dp * ( MeYeaYeYeaYe + YeaYeYeaYeMe )                       &
      & - 4._dp * ( YeMlaYeYeaYe + YeaYeMeYeaYe + YeaYeYeMlaYe )             &
      & - 8._dp * mH(1) * YeaYeYeaYe                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MeYeaYe + YeaYeMe + 4._dp * MH(1) * YeaYe &
      &                           + 2._dp * (YeMlaYe + AeaAE)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YeaYe + AeaAeYeaYe + YeaYEAeaAe     &
      &           + AeaYeYeaAe + YeaAeAeaYe + Tr3aAdAdaAeAe * YeaYe          &
      &           + Tr3aAdYdaAeYe * AeaYe + Tr3aYdAdaYeAe * YeaAe )          &
      & + hd(1) * ( MeYeaYe + YeaYeMe                                        &
      &           + 2._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe ) )            &
      & + hd(2) * YeaYe + hc(1) * AeaYe + hc(2) * YeaAe
    Do i1=1,3
     betaMe2(i1,i1) = betaMe2(i1,i1) + diagonal(1,2)
    End Do

    diagonal(3,2) = gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)        &
     &     + 3._dp * gauge2(2) * sig(2)                                 &
     &     + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)      &
     &                   + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2       &
     &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )     &
     &     + 24.84_dp * gauge2(1) * AbsGM2(1)
    betaMl2 = -2._dp * ( MlaYeYeaYeYe + aYeYeaYeYeMl )                   &
      & - 4._dp * ( aYeMeYeaYeYe + aYeYeMlaYeYe + aYeYeaYeMeYe )         &
      & - 8._dp * mH(1) * aYeYeaYeYe                                     &
      & - Tr3aYdYdaYeYe * ( MlaYeYe + aYeYeML + 4._dp * MH(1) * aYeYe    &
      &                   + 2._dp * (aYeMeYe + aAeAE)  )                 &
      & - 4._dp * (aAeAeaYeYe + aYeYeaAeAe + aAeYeaYeAe + aYeAeaAeYe )   &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYeYe +  Tr3aAdAdaAeAe * aYeYe  &
      &           + Tr3aAdYdaAeYe * aYeAe + Tr3aYdAdaYeAe * aAeYe     )  &
      & + 1.2_dp*gauge2(1) * ( MlaYeYe + aYeYeMl                         &
      &                      + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe &
      &                                - Mhlf(1) * aAeYe                 &
      &                                - Conjg(Mhlf(1)) * aYeAe ) )      &
      & + 4.8_dp * AbsGM2(1) * aYeYe
    Do i1=1,3
     betaMl2(i1,i1) = betaMl2(i1,i1) + diagonal(3,2)
    End Do

    diagonal(4,2) = 0.8_dp * gauge2(1) * (S2 + sig(1)/3._dp )             &
      &   + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp   &
      &   + 8.08e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                    &
      &   + 1.28e2_dp * ( gauge2(1) * AbsGM2(3)                           &
      &                 + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)     &
      &                   * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)+0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) + 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) - 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMd2 = -2._dp * ( MdYdaYdYdaYd + YdaYdYdaYdMd )                       &
      & - 4._dp * ( YdMqaYdYdaYd + YdaYdMdYdaYd + YdaYdYdMqaYd )             &
      & - 8._dp * mH(1) * YdaYdYdaYd                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MdYdaYd + YdaYdMd + 4._dp * MH(1) * YdaYd &
      &                           + 2._dp * (YdMqaYd + AdaAD)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YdaYd + AdaAdYdaYd + YdaYDAdaAd     &
      &           + AdaYdYdaAd + YdaAdAdaYd + Tr3aAdAdaAeAe * YdaYd          &
      &           + Tr3aAdYdaAeYe * AdaYd + Tr3aYdAdaYeAe * YdaAd )          &
      & + hd(1) * ( MdYdaYd + YdaYdMd                                        &
      &           + 2._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd ) )            &
      & + hd(2) * YdaYd + hc(1) * AdaYd + hc(2) * YdaAd                      &
      & - 2._dp * (MdYdaYuYuaYd + YdaYuYuaYdMd )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YdaYuYuaYd + YdMqaYuYuaYd+ YdaYuYuMqaYd  &
      &           + YdaYuMuYuaYd + AdaAuYuaYd + YdaYuAuaAd                   &
      &           + AdaYuYuaAd + YdaAuAuaYd )
    Do i1=1,3
     betaMd2(i1,i1) = betaMd2(i1,i1) + diagonal(4,2)
    End Do

    diagonal(5,2) = 0.2_dp * gauge2(1) * (2._dp * S2 + sig(1)/3._dp )     &
      &  + 3._dp * gauge2(2) * sig(2)                                     &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp    &
      &  + 1.99e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                     &
      &  + 32._dp * ( gauge2(1) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)         &
      &              * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp      &
      &  + 33._dp * gauge2(2) * AbsGM2(2)                                 &
      &  + 32._dp * ( gauge2(2) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(2)  + gauge2(3)*gauge2(2)        &
      &              * Real( Mhlf(2) * Conjg(Mhlf(3)),dp ) )                 &
      &  + 0.4_dp * ( gauge2(2) * AbsGM2(1)                               &
      &             + gauge2(1) * AbsGM2(2) + gauge2(1)*gauge2(2)         &
      &              * Real( Mhlf(2) * Conjg(Mhlf(1)),dp ) )
    hd(1) = 1.6_dp * AbsGM2(1)
    hc(1) = - 0.8_dp * g2Mi(1)
    hc(2) = Conjg( hc(1) )
    hd(2) = 3.2_dp * AbsGM2(1)
    hc(3) = - 1.6_dp * g2Mi(1)
    hc(4) = Conjg( hc(3) )
    betaMq2 = -2._dp * ( MqaYdYdaYdYd + aYdYdaYdYDMq )                       &
      & - 4._dp * ( aYdMdYdaYdYd + aYdYdMqaYdYd + aYdYdaYdMdYd )             &
      & - 8._dp * mH(1) * aYdYdaYdYd                                         &
      & - Tr3aYdYdaYeYe * ( MqaYdYd + aYdYdMq + 4._dp * MH(1) * aYdYd        &
      &                   + 2._dp * (aYdMdYd + aAdAd)  )                     &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYdYd                               &
      &           + 2._dp * ( aAdAdaYdYd + aYdYDaAdAd + aAdYdaYdAd           &
      &                     + aYdAdaAdYd )                                   &
      &           + Tr3aAdAdaAeAe * aYdYd + Tr3aAdYdaAeYe * aYdAd            &
      &           + Tr3aYdAdaYeAe * aAdYd )                                  &
      & + 0.4_dp*gauge2(1) * ( MqaYdYd + aYdYdMq                             &
      &                      + 2._dp * ( mH(1) * aYdYd + aYdMdYd + aAdAd ) ) &
      & + hd(1) * aYdYd + hc(1) * aAdYd + hc(2) * aYdAd                      &
      & - 2._dp * ( MqaYuYuaYuYu + aYuYuaYuYuMq )                            &
      & - 4._dp * ( aYuMuYuaYuYu + aYuYuMqaYuYu + aYuYuaYuMuYu )             &
      & - 8._dp * MH(2) * aYuYuaYuYu                                         &
      & - 3._dp * TraceY(4) * ( MqaYuYu + aYuYuMq + 4._dp * MH(2) * aYuYu    &
      &                       + 2._dp * (aYuMuYu + aAuAu) )                  &
      & - 2._dp * ( Tr3MqaYuYu3aYuMu * aYuYu                                 &
      &           + 2._dp * ( aAuAuaYuYu + aYuYuaAuAu                        &
      &                     + aAuYuaYuAu + aYuAuaAuYu )                      &
      &           + Tr3aAuAu * aYuYu + Tr3aAuYu * aYuAu + Tr3aYuAu * aAuYu ) &
      & + 0.8_dp*gauge2(1) * ( MqaYuYu + aYuYuMq                             &
      &                      + 2._dp * ( MH(2) * aYuYu + aYuMuYu + aAuAu ) ) &
      & + hd(2) * aYuYu + hc(3) * aAuYu + hc(4) * aYuAu
    Do i1=1,3
     betaMq2(i1,i1) = betaMq2(i1,i1) + diagonal(5,2)
    End Do

    diagonal(6,2) = 1.6_dp * gauge2(1) * (2._dp*sig(1)/3._dp - S2)           &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp       &
      &  + 3424._dp * gauge2(1) * AbsGM2(1) / 75._dp                         &
      &  + 512._dp * ( gauge2(1) * AbsGM2(3)                                 &
      &              + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)           &
      &                * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)-0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMu2 = -2._dp * ( MuYuaYuYuaYu + YuaYuYuaYuMu )                       &
      & - 4._dp * ( YuMqaYuYuaYu + YuaYuMuYuaYu + YuaYuYuMqaYu )             &
      & - 8._dp * MH(2) * YuaYuYuaYu                                         &
      & - 6._dp * TraceY(4) * ( MuYuaYu + YuaYuMu + 4._dp * MH(2) * YuaYu    &
      &                       + 2._dp * (YuMqaYu + AuaAu)  )                 &
      & - 4._dp * ( Tr3MqaYuYu3aYuMu * YuaYu + AuaAuYuaYu + YuaYuAuaAu       &
      &           + AuaYuYuaAu + YuaAuAuaYu + Tr3aAuAu * YuaYu               &
      &           + Tr3aAuYu * AuaYu + Tr3aYuAu * YuaAu )                    &
      & + hd(1) * ( MuYuaYu + YuaYuMu                                        &
      &           + 2._dp * ( MH(2) * YuaYu + YuMqaYu + AuaAu ) )            &
      & + hd(2) * YuaYu + hc(1) * AuaYu + hc(2) * YuaAu                      &
      & - 2._dp * (MuYuaYdYdaYu + YuaYdYdaYuMu )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YuaYdYdaYu + YuMqaYdYdaYu + YuaYdYdMqaYu &
      &           + YuaYdMdYdaYu + AuaAdYdaYu + YuaYdAdaAu                   &
      &           + AuaYdYdaAu + YuaAdAdaYu )
    Do i1=1,3
     betaMu2(i1,i1) = betaMu2(i1,i1) + diagonal(6,2)
    End Do

   End If 

  !------------------------------------------
  ! beta functions for Higgs mass parameters
  !------------------------------------------
   traceMH1(1) = mH(1) * (TraceY(1) + 6._dp * lam12)            &
             & + Real( cTrace(YeMlaYe),dp ) &
             & + Real( cTrace(aYeMeYe),dp ) + TraceA(1)
   traceMH1(2) = mH(1) * TraceY(3) + Real( cTrace(YdMqaYd),dp )   &
             & + Real( cTrace(aYdMdYd),dp ) + TraceA(3) + Alam12  &
             & + lam12 * MT(1)
   betamH11 = 6._dp * TraceMH1(2) + 2._dp * TraceMH1(1)      &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1

   traceMH2(1) = 3._dp * (2._dp* mH(2)+ MT(2)) * lam22 + 3._dp * Alam22
   traceMH2(2) = mH(2) * TraceY(4) + Real( cTrace(YuMqaYu),dp )  &
             & + Real( cTrace(aYuMuYu),dp ) + TraceA(4)
   betamH21 = 2._dp * TraceMH2(1) + 6._dp * TraceMH2(2)       &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) + 0.6_dp * S1

   betaMT2(1) = MT(1) * (lam12 + TraceY(2)) + 2._dp * mH(1) * lam12   &
            & + 2._dp * Real( cTrace(aYTMlYT),dp ) + TraceA(2) + Alam12 &
            & - 2.4_dp * AbsGM2(1) - 8._dp * AbsGM2(2)

   betaMT2(2) = MT(2) * lam22 + 2._dp * mH(2) * lam22 + Alam22 &
            & - 2.4_dp * AbsGM2(1) - 8._dp * AbsGM2(2)

   betaMT2 = 2._dp * betaMT2
 
   betaMZ2(2) = - (0.4_dp * AbsGM2(1) + 32._dp * AbsGM2(3) ) / 3._dp  &
              & - 6._dp  * AbsGM2(2)

   betaMZ2(1) = 2._dp * MZ(1) * TraceY(5) + 2._dp * TraceA(5) + betaMZ2(2) &
            & + 2._dp * Real( cTrace(aYZMdYZ) + cTrace(YZMlaYZ) ,dp ) 

   betaMS2(2) = - (3.2_dp * AbsGM2(1) + 40._dp * AbsGM2(3) ) / 3._dp

   betaMS2(1) = MS(1) * TraceY(6) + 2._dp * Real( cTrace(aYSMdYS),dp ) &
            & + TraceA(6) + betaMS2(2)

   betaMS2 = 2._dp * betaMS2


   If (TwoLoopRGE) Then
    traceMH1(3) = MH(1) * (6._dp*TraceY2(2) + 2._dp*TraceY2(1) + TraceY2(4) ) &
              & + MH(2) * TraceY2(4)                                          &
              & + 6._dp * ( Real( cTrace(MqaYdYdaYdYd),dp )                   &
              &           + Real( cTrace(aYdMdYdaYdYd),dp )                   &
              &           + Real( cTrace(aAdAdaYdYd),dp )                     &
              &           + Real( cTrace(aAdYdaYdAd),dp )  )                  &
              & + 2._dp * ( Real( cTrace(MlaYeYeaYeYe),dp )                   &
              &           + Real( cTrace(aYeMeYeaYeYe),dp )                   &
              &           + Real( cTrace(aAeAeaYeYe),dp )                     &
              &           + Real( cTrace(aAeYeaYeAe),dp )  )                  &
       & + Real( cTrace(YdMqaYuYuaYd),dp ) + Real( cTrace(YdaYuMuYuaYd),dp ) &
       & + Real( cTrace(YdaYuYuMqaYd),dp ) + Real( cTrace(YuaYdMdYdaYu),dp ) &
       & + Real( cTrace(YdaAuAuaYd),dp ) + Real( cTrace(AdaYuYuaAd),dp )     &
              & + Real( cTrace(AdaAuYuaYd),dp ) + Real( cTrace(YdaYuAuaAd),dp )
    betaMH12 = - 6._dp * traceMH1(3)                                       &
      &   + (32._dp*gauge2(3) - 0.8_dp*gauge2(1) ) * traceMH1(2)           &
      &   + 64._dp * ( AbsGM2(3) * TraceY(3)                               &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(3),dp ) )     &
      &   - 1.6_dp * ( AbsGM2(1) * TraceY(3)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(3),dp ) )     &
      &   + 2.4_dp*gauge2(1) * traceMH1(1)                                 &
      &   + 4.8_dp * ( AbsGM2(1) * TraceY(1)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(1),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)                   &
      &   + 3._dp * gauge2(2) * sig(2)                                     &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)          &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2           &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

    traceMH2(2) = MH(2) * ( 6._dp * TraceY2(3) + TraceY2(4) )                 &
      &  + MH(1) * TraceY2(4)                                                 &
 &  + 6._dp * ( Real(cTrace(MqaYuYuaYuYu),dp) + Real(cTrace(aYuMuYuaYuYu),dp) &
 &            + Real(cTrace(aAuAuaYuYu),dp) + Real(cTrace(aAuYuaYuAu),dp)  )  &
 &  + Real( cTrace(YuMqaYdYdaYu),dp ) + Real( cTrace(YuaYdMdYdaYu),dp )       &
 &  + Real( cTrace(YuaYdYdMqaYu),dp ) + Real( cTrace(YdaYuMuYuaYd),dp )       &
 &  + Real( cTrace(YuaAdAdaYu),dp ) + Real( cTrace(AuaYdYdaAu),dp )           &
      &  + Real( cTrace(AuaAdYdaYu),dp ) + Real( cTrace(YuaYdAdaAu),dp )
    betaMH22 = - 6._dp * traceMH2(2)                                          &
      &   + (32._dp*gauge2(3) + 1.6_dp*gauge2(1) ) * traceMH2(1)              &
      &   + 64._dp * ( AbsGM2(3) * TraceY(4)                                  &
      &          - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(4),dp ) )        &
      &   + 3.2_dp * ( AbsGM2(1) * TraceY(4)                                  &
      &          - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(4),dp ) )        &
      &   + gauge2(1) * ( 0.6_dp * sig(1) + 1.2_dp * S2)                      &
      &   + 3._dp * gauge2(2) * sig(2)                                        &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)             &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2              &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

   End If
!------------------------------
! beta function for the 15-plet
!------------------------------
   betaMT15 = MT15 * ( TraceY(2) + lam12 + lam22 &
            &        - 2.4_dp * gauge2(1) - 8._dp * gauge2(2)  )
   betaMZ15 = MZ15 * ( TraceY(5) - gauge2(1)/15._dp - 3._dp * gauge2(2) &
                     & - 16._dp * gauge2(3) / 3._dp  )
   betaMS15 = MS15 * ( TraceY(6)  &
                     & - (3.2_dp * gauge2(1) + 40._dp * gauge2(3))/3._dp  )
!-----------------------------
! beta functions for mu and B
!-----------------------------
   TraceMue(1) = 3._dp * (TraceY(3)+TraceY(4)+lam12+lam22) + TraceY(1)  &
             & - 3._dp * gauge2(2) - 0.6_dp * gauge2(1)
   betaMue1 = mue * TraceMue(1)

   TraceB(1) = 6._dp * (TraceaYA(3)+TraceaYA(4) + lam1Alam1 + lam2Alam2) &
           & + 2._dp * TraceaYA(1)                         &
           & + 6._dp * g2Mi(2) + 1.2_dp * g2Mi(1)
   betaB1 = mue * TraceB(1) + B * TraceMue(1)

   If (TwoLoopRGE) Then
    TraceMue(2) = - 3._dp * ( 3._dp * (TraceY2(2) + TraceY2(3) )          &
      &                     + 2._dp * TraceY2(4) + TraceY2(1) )           &
      &         + (16._dp * gauge2(3) + 0.8_dp * gauge2(1)) * TraceY(4)   &
      &         + (16._dp * gauge2(3) - 0.4_dp * gauge2(1)) * TraceY(3)   &
      &         + 1.2_dp * gauge2(1) * TraceY(1)                          &
      &         + 7.5_dp * gauge2(2)**2                                   &
      &         + 1.8_dp * gauge2(2) * gauge2(1)                          &
      &         + 4.14_dp * gauge2(1)**2
    betaMue2 = mue * TraceMue(2)

    TraceB(1) = cTrace( 3._dp * ( Matmul(AuaYu,YuaYu) + Matmul(AdaYd,YdaYd) ) &
              &       + Matmul(AeaYe,YeaYe) + Matmul(aYuAu,aYdYd)             &
              &       + Matmul(aYdAd,aYuYu) ) 
    TraceB(2) = -12._dp * TraceB(1)                                           &
      &   + (32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)          &
      &   + (32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)          &
      &   + 2.4_dp * gauge2(1) * TraceaYA(1)                                  &
      &   - ( 32._dp * g2Mi(3)                                    &
      &     + 1.6_dp * g2Mi(1) ) * TraceY(4)                      &
      &   - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)   &
      &   - 2.4_dp * g2Mi(1) * TraceY(1)                          &
      &   - 30._dp * gauge2(2)**2 * Mhlf(2)                                   &
      &   - 3.6_dp * gauge2(2) * gauge2(1) * (Mhlf(1) + Mhlf(2) )             &
      &   - 16.56_dp * gauge2(1)**2 * Mhlf(1) 
    betaB2 = mue * TraceB(2) + B * TraceMue(2)

   End If

  !--------------------------------
  ! neutrino dim. 5 operator
  !--------------------------------
  diagonal(5,1) = 6._dp * TraceY(4) - 1.2_dp * gauge2(1) - 6._dp * gauge2(2)
  betaMnu1 = Matmul( Transpose(aYeYe), Mnu) + Matmul(Mnu, aYeYe)  &
          & + diagonal(5,1) * Mnu
  

 !---------------
 ! 2-loop RGEs
 !---------------
  b_1a = b_1 + Delta_b_1
  b_2a = b_2 + Delta_b_2

  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   If (Maxval(Delta_b_2).Eq.0._dp) Then
    Dgauge = oo16pi2 * gauge * gauge2  &
         & * ( b_1a + oo16pi2 * ( Matmul(b_2a,gauge2) - a_2(:,1) * TraceY(1) &
         &                      - Matmul(a_2(:,2:3),TraceY(3:4))  ))
   Else
    Dgauge = oo16pi2 * gauge * gauge2  &
         & * ( b_1a + oo16pi2 * (Matmul(b_2a,gauge2) - Matmul(a_2b(:,1:6),TraceY) )  &
         &                      - a_2b(:,7)*lam12 - a_2b(:,8)*lam22 ) 
   End If
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYT = oo16pi2 * betaYT1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
   DYZ = oo16pi2 * betaYZ1 
   DYS = oo16pi2 * betaYS1 
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   Do i1 = 1,3    
    sumI = 0._dp
    Do i2=1,3
     sumI = sumI + b_2a(i1,i2) * gauge2(i2) * (Mhlf(i1) + Mhlf(i2) ) 
    End Do
    Do i2=1,6
     sumI = sumI + a_2b(i1,i2) * ( TraceaYA(i2) - Mhlf(i1)*TraceY(i2) )
    End Do
    sumI = sumI + a_2b(i1,7) * lam1Alam1 + a_2b(i1,8) * lam2Alam2
    DMhlf(i1) = oo8pi2 * gauge2(i1) * ( b_1a(i1) * Mhlf(i1) + oo16pi2 * sumI)
   End Do
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * ( betaAe1 + oo16pi2 * betaAe2 )
   DAd = oo16pi2 * ( betaAd1 + oo16pi2 * betaAd2 )
   DAu = oo16pi2 * ( betaAu1 + oo16pi2 * betaAu2 )
   DAT = oo16pi2 * betaAT1
   DAZ = oo16pi2 * betaAZ1
   DAS = oo16pi2 * betaAS1
   DAlam1 = oo16pi2 * betaAlam11
   DAlam2 = oo16pi2 * betaAlam21
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * ( betaMe1 + oo16pi2 * betaMe2 )
   DMl = oo16pi2 * ( betaMl1 + oo16pi2 * betaMl2 )
   DMd = oo16pi2 * ( betaMd1 + oo16pi2 * betaMd2 )
   DMq = oo16pi2 * ( betaMq1 + oo16pi2 * betaMq2 )
   DMu = oo16pi2 * ( betaMu1 + oo16pi2 * betaMu2 )
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * ( betaMH11 + oo16pi2 * betaMH12 )
   DmH(2) = oo16pi2 * ( betaMH21 + oo16pi2 * betaMH22 )
   DmT2 = oo16pi2 * betaMT2
   DmZ2 = oo16pi2 * betaMZ2
   DmS2 = oo16pi2 * betaMS2
  !----------
  ! 15-plet
  !----------
  DMT15 = oo16pi2 * betaMT15
  DMZ15 = oo16pi2 * betaMZ15
  DMS15 = oo16pi2 * betaMS15
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * ( betaMue1 + oo16pi2 * betaMue2 )
   DB = oo16pi2 * ( betaB1 + oo16pi2 * betaB2 )

 !----------------------
 ! left neutrino mass
 !----------------------
   DMnu = oo16pi2 * betaMnu1
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1a 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYT = oo16pi2 * betaYT1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
   DYZ = oo16pi2 * betaYZ1 
   DYS = oo16pi2 * betaYS1 
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   DMhlf = oo8pi2 * gauge2 * b_1a * Mhlf
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * betaAe1
   DAd = oo16pi2 * betaAd1
   DAu = oo16pi2 * betaAu1
   DAT = oo16pi2 * betaAT1
   DAZ = oo16pi2 * betaAZ1
   DAS = oo16pi2 * betaAS1
   DAlam1 = oo16pi2 * betaAlam11
   DAlam2 = oo16pi2 * betaAlam21
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * betaMe1
   DMl = oo16pi2 * betaMl1
   DMd = oo16pi2 * betaMd1
   DMq = oo16pi2 * betaMq1
   DMu = oo16pi2 * betaMu1
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * betaMH11
   DmH(2) = oo16pi2 * betaMH21
   DmT2 = oo16pi2 * betaMT2
   DmZ2 = oo16pi2 * betaMZ2
   DmS2 = oo16pi2 * betaMS2
  !----------
  ! 15-plet
  !----------
  DMT15 = oo16pi2 * betaMT15
  DMZ15 = oo16pi2 * betaMZ15
  DMS15 = oo16pi2 * betaMS15
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * betaMue1
   DB = oo16pi2 * betaB1
  !----------------------
  ! left neutrino mass
  !----------------------
   DMnu = oo16pi2 * betaMnu1
  End If


  !---------------------------------------
  ! to avoid numerical problems in odeint
  !---------------------------------------
  Do i1=1,3
   DMe(i1,i1) = Real(DMe(i1,i1),dp)
   DMl(i1,i1) = Real(DMl(i1,i1),dp)
   DMd(i1,i1) = Real(DMd(i1,i1),dp)
   DMu(i1,i1) = Real(DMu(i1,i1),dp)
   DMq(i1,i1) = Real(DMq(i1,i1),dp)
  End Do

  Call ParametersToG5(Dgauge, DYe, DYT, DYd, DYu, DYZ, DYS, Dlam1, Dlam2    &
          & , DMhlf, DAe, DAT, DAd, DAu, DAZ, DAS, DAlam1, DAlam2, DMe, DMl &
          & , DMd, DMq, DMu, DMh, DMT2, DmZ2, DmS2, DMT15, DMZ15, DMS15     &
          & , Dmue, DB, DMnu, f)

  Iname = Iname - 1

 End Subroutine rge356


 Subroutine RGE_SU5(len,T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings, soft SUSY breaking parameters
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2, SumI
  Real(dp) :: g5, g52, b5, lam2, lamp2, GammaH, GammaHbar, GammaS, TraceY(3) &
    & , Q
  Complex(dp) :: lam, lamp
  Complex(dp), Dimension(3,3) :: Y_u, Y_d, Y_N, YdaYd, aYdYd, YnaYn, aYnYn &
    & , YuaYu, aYuYu, GammaT, GammaF, GammaN, aY_u, aY_d, aY_N 
  Real(dp) :: dG5
  Complex(dp) :: Dlam, Dlamp
  Complex(dp), Dimension(3,3) :: DYu1, DYn1, DYd1

  Complex(dp) :: MuH, DMuH, MuSig, DMuSig, M5, DM5, gM5, GammaAH, GammaAHbar &
    & , GammaAS
  Complex(dp), Dimension(3,3) :: A_u, A_d, A_N, AdaYd, aYdAd, aYnAn &
    & , AuaYu, aYuAu, GammaAT, GammaAF, GammaAN
  Complex(dp) :: Alam, Alamp, Alamlam, Alamplamp, DAlam, DAlamp, TraceAY(3)
  Complex(dp), Dimension(3,3) :: DAu1, DAn1, DAd1

  Complex(dp), Dimension(3,3) :: MT2, MF2, MN2, DMT2, DMF2, DMN2    &
    & , aA_u, aA_d, aA_N, AdaAd, aAdAd, AnaAn, aAnAn, AuaAu, aAuAu  &
    & , MT2YuaYu, YuaYuMT2, YuMT2aYu, MT2YdaYd, YdaYdMT2, aYdMT2Yd  &
    & , YdMF2aYd, MF2aYdYd, aYdYdMF2, MF2aYnYn, aYnYnMF2, YnMF2aYn  &
    & , aYnYnMN2, MN2aYnYn, aYnMN2Yn, aYnMF2Yn

  Complex(dp) :: TraceA(3)
  Real(dp) :: M52, Alam2, Alamp2, MH2, MHbar2, MS2, DMH2, DMHbar2, DMS2
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'RGE_SU5'

  q = t

  OnlyDiagonal = .Not.GenerationMixing

  !----------------------------------
  ! gauge couplings beta functions
  !----------------------------------
  g5 = gy(1)
  g52 = g5**2
  b5 = -3._dp

  !----------------------------
  ! Yukawa couplings
  !----------------------------
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    Y_u(i1,i2) = Cmplx( gy(SumI-6), gy(SumI-5),dp )
    Y_d(i1,i2) = Cmplx( gy(SumI+12), gy(SumI+13),dp )
    Y_n(i1,i2) = Cmplx( gy(SumI+30), gy(SumI+31),dp )
   End Do
  End Do

  lam = Cmplx( gy(56), gy(57),dp )
  lamp = Cmplx( gy(58), gy(59),dp )
  lam2 = Abs(lam)**2
  lamp2 = Abs(lamp)**2
  !-----------------
  ! beta functions
  !-----------------
  Call Adjungate(Y_d,aY_d)
  Call Adjungate(Y_n,aY_n)
  Call Adjungate(Y_u,aY_u)

  aYdYd = Matmul(aY_d, Y_d)
  aYnYn = Matmul(aY_n, Y_n)
  aYuYu = Matmul(aY_u, Y_u)
  YdaYd = Matmul(Y_d, aY_d)
  YnaYn = Matmul(Y_n, aY_n)
  YuaYu = Matmul(Y_u, aY_u)
  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aYdYd(i1,i1) = Real(aYdYd(i1,i1), dp)
   aYnYn(i1,i1) = Real(aYnYn(i1,i1), dp)
   aYuYu(i1,i1) = Real(aYuYu(i1,i1), dp)
   YdaYd(i1,i1) = Real(YdaYd(i1,i1), dp)
   YnaYn(i1,i1) = Real(YnaYn(i1,i1), dp)
   YuaYu(i1,i1) = Real(YuaYu(i1,i1), dp)
  End Do

  TraceY(1) = Real( cTrace(aYnYn),dp )
  TraceY(2) = Real( cTrace(aYdYd),dp )
  TraceY(3) = Real( cTrace(aYuYu),dp )

  GammaT = 2._dp * YdaYd + 3._dp * YuaYu
  GammaF = 4._dp * aYdYd + aYnYn
  GammaN = 5._dp * aYnYn

  Do i1=1,3
   GammaT(i1,i1) = GammaT(i1,i1) - 7.2_dp * g52
   GammaF(i1,i1) = GammaF(i1,i1) - 4.8_dp * g52
  End Do

  GammaHbar = 4._dp * TraceY(2) + 4.8_dp * (lam2 - g52)
  GammaH = 3._dp * TraceY(3) + TraceY(1) + 4.8_dp * (lam2 - g52)
  GammaS = 1.05_dp * lamp2 + lam2 - 10._dp * g52

  Dg5 = oo16pi2 * b5 * g5 * g52
 
  DYu1 = Y_u * GammaH + Matmul(GammaT,Y_u) &
       & + Matmul(Y_u,Transpose(GammaT))
  DYd1 = Y_d * GammaHbar + Matmul(GammaT,Y_d) &
       & + Matmul(Y_d,GammaF)
  DYn1 = Y_n * GammaH + Matmul(GammaT,Y_n) &
       & + Matmul(Y_n,GammaN)

  DYn1 = oo16pi2 * DYn1 
  DYd1 = oo16pi2 * DYd1 
  DYu1 = oo16pi2 * DYu1 

  Dlam = oo16pi2 * lam * (GammaH + GammaHbar + GammaS)
  Dlamp = oo16pi2 * lamp * 3._dp * GammaS

  f(1) = Dg5
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    f(SumI-6) = Real(DYu1(i1,i2),dp)
    f(SumI-5) = Aimag(DYu1(i1,i2))
    f(SumI+12) = Real(DYd1(i1,i2),dp)
    f(SumI+13) = Aimag(DYd1(i1,i2))
    f(SumI+30) = Real(DYn1(i1,i2),dp)
    f(SumI+31) = Aimag(DYn1(i1,i2))
   End Do
  End Do
  f(56) = Real(Dlam, dp)
  f(57) = Aimag(Dlam)
  f(58) = Real(Dlamp, dp)
  f(59) = Aimag(Dlamp)

  If (len.Eq.59) Then
   Iname = Iname - 1
   Return
  End If

  !----------------------------------
  ! bilinear Higgs parameters
  !----------------------------------
  MuH = Cmplx( gy(60), gy(61),dp )
  MuSig = Cmplx( gy(62), gy(63),dp )

  DMuH = oo16pi2 * MuH * (GammaH + GammaHbar)
  DMuSig = oo16pi2 * MuSig * 2._dp * GammaS

  f(60) = Real(DMuH, dp)
  f(61) = Aimag(DMuH)
  f(62) = Real(DMuSig, dp)
  f(63) = Aimag(DMuSig)

  If (len.Eq.63) Then
   Iname = Iname - 1
   Return
  End If

  !----------------------------------
  ! gaugino beta function
  !----------------------------------
  M5 = Cmplx( gy(64), gy(65) , dp)
  gM5 = g52 * M5

  !----------------------------
  ! A-parameters
  !----------------------------
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    A_u(i1,i2) = Cmplx( gy(SumI+58), gy(SumI+59),dp )
    A_d(i1,i2) = Cmplx( gy(SumI+76), gy(SumI+77),dp )
    A_n(i1,i2) = Cmplx( gy(SumI+94), gy(SumI+95),dp )
   End Do
  End Do

  Alam = Cmplx( gy(120), gy(121),dp )
  Alamp = Cmplx( gy(122), gy(123),dp )
  Alamlam = Alam * Conjg(lam)
  Alamplamp = Alamp * Conjg(lamp)
  !-----------------
  ! beta functions
  !-----------------

  aYdAd = Matmul(aY_d, A_d)
  aYnAn = Matmul(aY_n, A_n)
  aYuAu = Matmul(aY_u, A_u)
  AdaYd = Matmul(A_d, aY_d)
  AuaYu = Matmul(A_u, aY_u)

  TraceAY(1) = cTrace(aYnAn)
  TraceAY(2) = cTrace(aYdAd)
  TraceAY(3) = cTrace(aYuAu)

  GammaAT = -2._dp * AdaYd - 3._dp * AuaYu
  GammaAF = -4._dp * aYdAd - aYnAn
  GammaAN = -5._dp * aYnAn

  Do i1=1,3
   GammaAT(i1,i1) = GammaAT(i1,i1) - 7.2_dp * gM5
   GammaAF(i1,i1) = GammaAF(i1,i1) - 4.8_dp * gM5
  End Do

  GammaAHbar = - 4._dp * TraceAY(2) - 4.8_dp * (Alamlam + gM5)
  GammaAH = -3._dp * TraceAY(3) - TraceAY(1) - 4.8_dp * (Alamlam + gM5)
  GammaAS = - 1.05_dp * Alamplamp - Alamlam - 10._dp * gM5

  DM5 = oo8pi2 * b5 * gM5
 
  DAu1 = A_u * GammaH + Matmul(GammaT,A_u)               &
       & + Matmul(A_u,Transpose(GammaT))                 &
       & - 2._dp * ( Y_u * GammaAH + Matmul(GammaAT,Y_u) &
       &           + Matmul(Y_u,Transpose(GammaAT)) ) 
  DAd1 = A_d * GammaHbar + Matmul(GammaT,A_d)               &
       & + Matmul(A_d,GammaF)                               &
       & - 2._dp * ( Y_d * GammaAHbar + Matmul(GammaAT,Y_d) &
       &           + Matmul(Y_d,GammaAF) )
  DAn1 = A_n * GammaH + Matmul(GammaT,A_n)              &
       & + Matmul(A_n,GammaN)                           &
       & - 2._dp * (Y_n * GammaAH + Matmul(GammaAT,Y_n) &
       &           + Matmul(Y_n,GammaAN) )

  DAn1 = oo16pi2 * DAn1 
  DAd1 = oo16pi2 * DAd1 
  DAu1 = oo16pi2 * DAu1 

  DAlam = oo16pi2 * ( Alam * (GammaH + GammaHbar + GammaS)       &
        &           - 2._dp * lam * (GammaAH + GammaAHbar + GammaAS) )
  DAlamp = oo16pi2 * ( Alamp * 3._dp * GammaS - lamp * 6._dp * GammaAS )

  f(64) = Real(DM5, dp)
  f(65) = Aimag( DM5 )

  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    f(SumI+58) = Real(DAu1(i1,i2),dp)
    f(SumI+59) = Aimag(DAu1(i1,i2))
    f(SumI+76) = Real(DAd1(i1,i2),dp)
    f(SumI+77) = Aimag(DAd1(i1,i2))
    f(SumI+94) = Real(DAn1(i1,i2),dp)
    f(SumI+95) = Aimag(DAn1(i1,i2))
   End Do
  End Do
  f(120) = Real(DAlam, dp)
  f(121) = Aimag(DAlam)
  f(122) = Real(DAlamp, dp)
  f(123) = Aimag(DAlamp)

  !----------------------------
  ! scalar masses squared
  !----------------------------
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    MT2(i1,i2) = Cmplx( gy(SumI+116), gy(SumI+117),dp )
    MF2(i1,i2) = Cmplx( gy(SumI+134), gy(SumI+135),dp )
    MN2(i1,i2) = Cmplx( gy(SumI+152), gy(SumI+153),dp )
   End Do
  End Do

  MHbar2 = gy(178)
  MH2 = gy(179)
  mS2 = gy(180) 

  Call Adjungate(A_d,aA_d)
  Call Adjungate(A_n,aA_n)
  Call Adjungate(A_u,aA_u)

  aAdAd = Matmul(aA_d, A_d)
  aAnAn = Matmul(aA_n, A_n)
  aAuAu = Matmul(aA_u, A_u)
  AdaAd = Matmul(A_d, aA_d)
  AnaAn = Matmul(A_n, aA_n)
  AuaAu = Matmul(A_u, aA_u)
  
  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aAdAd(i1,i1) = Real(aAdAd(i1,i1), dp)
   aAnAn(i1,i1) = Real(aAnAn(i1,i1), dp)
   aAuAu(i1,i1) = Real(aAuAu(i1,i1), dp)
   AdaAd(i1,i1) = Real(AdaAd(i1,i1), dp)
   AnaAn(i1,i1) = Real(AnaAn(i1,i1), dp)
   AuaAu(i1,i1) = Real(AuaAu(i1,i1), dp)
  End Do

  Alam2 = Abs(Alam)**2
  Alamp2 = Abs(Alamp)**2

  TraceA(1) = Real( cTrace(aAnAn),dp )
  TraceA(2) = Real( cTrace(aAdAd),dp )
  TraceA(3) = Real( cTrace(aAuAu),dp )

  MT2YuaYu = Matmul( MT2, YuaYu )
  YuaYuMT2 = Matmul( YuaYu, MT2 )
  YuMT2aYu = MatMul3( Y_u, MT2, aY_u, OnlyDiagonal )
  MT2YdaYd = Matmul( MT2, YdaYd )
  YdaYdMT2 = Matmul( YdaYd, MT2 )
  aYdMT2Yd = MatMul3( aY_d, MT2, Y_d, OnlyDiagonal )

  YdMF2aYd = MatMul3( Y_d, MF2, aY_d, OnlyDiagonal )
  MF2aYdYd = Matmul( MF2, aYdYd )
  aYdYdMF2 = Matmul( aYdYd, MF2 )
  MF2aYnYn = Matmul( MF2, aYnYn )
  aYnYnMF2 = Matmul( aYnYn, MF2 )
  YnMF2aYn = MatMul3( Y_n, MF2, aY_n, OnlyDiagonal )
  aYnMF2Yn = MatMul3( aY_n, MF2, Y_n, OnlyDiagonal )

  aYnYnMN2 = Matmul( aYnYn, MN2 )
  MN2aYnYn = Matmul( MN2, aYnYn )
  aYnMN2Yn = MatMul3( aY_n, MN2, Y_n, OnlyDiagonal )

  DMT2 = 4._dp * AdaAd + 6._dp * AuaAu                      &
     & + 3._dp * ( MT2YuaYu + 2._dp * YuMT2aYu + YuaYuMT2)  &
     & + 2._dp * ( MT2YdaYd + 2._dp * YdMF2aYd + YdaYdMT2)  &
     & + 6._dp * MH2 * YuaYu + 4._dp * MHbar2 * YdaYd
  DMF2 = 8._dp * aAdAd + 2._dp * aAnAn                      &
     & + 4._dp * (MF2aYdYd + 2._dp * aYdMT2Yd + aYdYdMF2 )  &
     & + MF2aYnYn + aYnYnMF2 + 2._dp * aYnMN2Yn             &
     & + 2._dp * MHbar2 * aYdYd+ 2._dp * MH2 * aYnYn
  DMN2 = 10._dp * aAnAn                                    &
     & + 5._dp * ( MN2aYnYn + 2._dp * aYnMF2Yn + aYnYnMN2) &
     & + 10._dp * MH2 * aYnYn

  M52 = g52 * Abs(M5)**2

  Do i1=1,3
   DMT2(i1,i1) = DMT2(i1,i1) - 28.8_dp * M52
   DMF2(i1,i1) = DMF2(i1,i1) - 19.2_dp * M52
  End Do

  DMH2 = 6._dp * TraceA(3) + 2._dp * TraceA(1)  &
     & + 9.6_dp * (Alam2 - 2._dp * M52            &
     &            + lam2 * (MH2 + MHbar2 + MS2) )                &
     & + 2._dp * (3._dp * TraceY(3) + TraceY(1) ) * MH2       &
     & + 2._dp * Real( 6._dp * cTrace(YuMT2aYu)               &
     &               + cTrace(YnMF2aYn) + cTrace(aYnMN2Yn) ,dp)
  DMHbar2 = 9.6_dp * ( Alam2 - 2._dp * M52                           &
     &               + lam2 * (MH2 + MHbar2 + MS2) )                 &
     & + 8._dp * (TraceY(2) * MHbar2 + TraceA(2)                     &
     &           + Real( cTrace( aYdMT2Yd ) + cTrace(YdMF2aYd) , dp) ) 
  DMS2 = 2.1_dp * Alamp2 + 2._dp * Alamp2 - 40._dp * M52 &
     & + 6.3_dp * lamp2 * MS2 + 2._dp * lam2 * (MH2 + MHbar2 + MS2) 

  DMT2 = oo16pi2 * DMT2
  DMF2 = oo16pi2 * DMF2
  DMN2 = oo16pi2 * DMN2
  DMH2 = oo16pi2 * DMH2
  DMHbar2 = oo16pi2 * DMHbar2
  DMS2 = oo16pi2 * DMS2

  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    f(SumI+116) = Real(DMT2(i1,i2),dp)
    f(SumI+117) = Aimag(DMT2(i1,i2))
    f(SumI+134) = Real(DMF2(i1,i2),dp)
    f(SumI+135) = Aimag(DMF2(i1,i2))
    f(SumI+152) = Real(DMN2(i1,i2),dp)
    f(SumI+153) = Aimag(DMN2(i1,i2))
   End Do
  End Do

  f(178) = DMHbar2
  f(179) = DMH2
  f(180) = DMS2
  
  Iname = Iname - 1

 End Subroutine RGE_SU5


#ifdef SEESAWIII

 Subroutine GToParameters111(g,g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3)

 Implicit None 
  Real(dp), Intent(in) :: g(111) 
  Real(dp),Intent(out) :: g1,g2,g3

  Complex(dp), Intent(out), Dimension(3,3) :: Yu, Yd, Ye, Yb3, Yw3, Yx3

  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'GToParameters111' 
 
  g1= g(1) 
  g2= g(2) 
  g3= g(3) 
  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3 )
    Yu(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
    Yd(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp) 
    Ye(i1,i2) = Cmplx( g(SumI+40), g(SumI+41), dp) 
    Yb3(i1,i2) = Cmplx( g(SumI+58), g(SumI+59), dp) 
    Yw3(i1,i2) = Cmplx( g(SumI+76), g(SumI+77), dp) 
    Yx3(i1,i2) = Cmplx( g(SumI+94), g(SumI+95), dp) 
   End Do 
  End Do 
  
  Iname = Iname - 1 
 
 End Subroutine GToParameters111


 Subroutine ParametersToG111(g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3,g)

 Implicit None 
  Real(dp), Intent(out) :: g(111) 
  Real(dp), Intent(in) :: g1,g2,g3

  Complex(dp), Intent(in), Dimension(3,3) :: Yu, Yd, Ye, Yb3, Yw3, Yx3

  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'ParametersToG111' 
 
  g(1) = g1  
  g(2) = g2  
  g(3) = g3  

  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3 )
    g(SumI+4) = Real(Yu(i1,i2), dp) 
    g(SumI+5) = Aimag(Yu(i1,i2)) 
    g(SumI+22) = Real(Yd(i1,i2), dp) 
    g(SumI+23) = Aimag(Yd(i1,i2)) 
    g(SumI+40) = Real(Ye(i1,i2), dp) 
    g(SumI+41) = Aimag(Ye(i1,i2)) 
    g(SumI+58) = Real(Yb3(i1,i2), dp) 
    g(SumI+59) = Aimag(Yb3(i1,i2)) 
    g(SumI+76) = Real(Yw3(i1,i2), dp) 
    g(SumI+77) = Aimag(Yw3(i1,i2)) 
    g(SumI+94) = Real(Yx3(i1,i2), dp) 
    g(SumI+95) = Aimag(Yx3(i1,i2)) 
   End Do 
  End Do 

  Iname = Iname - 1 
 
 End Subroutine ParametersToG111


 Subroutine GToParameters555(g,g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3,mue,MXM3,MWM3,     &
 & MGM3,MBM3,TYu,TYd,TYe,TYb3,TYw3,TYx3,Bmue,BMXM3,BMWM3,BMGM3,BMBM3,mq2,ml2,   &
 & mHd2,mHu2,md2,mu2,me2,mHw32,mHg3p2,mHb32,mHx32,mHxb32,MassB,MassWB,MassG,MnuL)

 Implicit None 
  Real(dp), Intent(in) :: g(573) 
  Real(dp),Intent(out) :: g1,g2,g3,mHd2,mHu2

  Complex(dp),Intent(out) :: mue, Bmue, MassB, MassWB, MassG
 
  Complex(dp),Intent(out), Dimension(3,3) :: Yu,Yd,Ye,Yb3,Yw3,Yx3              &
     & ,MXM3,MWM3,MGM3,MBM3,TYu,TYd,TYe,TYb3,TYw3,TYx3,BMXM3,BMWM3,BMGM3,BMBM3 &
     & ,mq2,ml2,md2,mu2,me2,mHw32,mHg3p2,mHb32,mHx32,mHxb32,MnuL
 
  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'GToParameters555' 
 
  g1= g(1) 
  g2= g(2) 
  g3= g(3) 
  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3 )
    Yu(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
    Yd(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp)
    Ye(i1,i2) = Cmplx( g(SumI+40), g(SumI+41), dp) 
    Yb3(i1,i2) = Cmplx( g(SumI+58), g(SumI+59), dp) 
    Yw3(i1,i2) = Cmplx( g(SumI+76), g(SumI+77), dp) 
    Yx3(i1,i2) = Cmplx( g(SumI+94), g(SumI+95), dp) 
    MXM3(i1,i2) = Cmplx( g(SumI+114), g(SumI+115), dp) 
    MWM3(i1,i2) = Cmplx( g(SumI+132), g(SumI+133), dp) 
    MGM3(i1,i2) = Cmplx( g(SumI+150), g(SumI+151), dp) 
    MBM3(i1,i2) = Cmplx( g(SumI+168), g(SumI+169), dp) 
    TYu(i1,i2) = Cmplx( g(SumI+186), g(SumI+187), dp) 
    TYd(i1,i2) = Cmplx( g(SumI+204), g(SumI+205), dp) 
    TYe(i1,i2) = Cmplx( g(SumI+222), g(SumI+223), dp) 
    TYb3(i1,i2) = Cmplx( g(SumI+240), g(SumI+241), dp) 
    TYw3(i1,i2) = Cmplx( g(SumI+258), g(SumI+259), dp) 
    TYx3(i1,i2) = Cmplx( g(SumI+276), g(SumI+277), dp) 
    BMXM3(i1,i2) = Cmplx( g(SumI+296), g(SumI+297), dp) 
    BMWM3(i1,i2) = Cmplx( g(SumI+314), g(SumI+315), dp) 
    BMGM3(i1,i2) = Cmplx( g(SumI+332), g(SumI+333), dp) 
    BMBM3(i1,i2) = Cmplx( g(SumI+350), g(SumI+351), dp) 
    mq2(i1,i2) = Cmplx( g(SumI+368), g(SumI+369), dp) 
    ml2(i1,i2) = Cmplx( g(SumI+386), g(SumI+387), dp) 
    md2(i1,i2) = Cmplx( g(SumI+406), g(SumI+407), dp) 
    mu2(i1,i2) = Cmplx( g(SumI+424), g(SumI+425), dp) 
    me2(i1,i2) = Cmplx( g(SumI+442), g(SumI+443), dp) 
    mHw32(i1,i2) = Cmplx( g(SumI+460), g(SumI+461), dp) 
    mHg3p2(i1,i2) = Cmplx( g(SumI+478), g(SumI+479), dp) 
    mHb32(i1,i2) = Cmplx( g(SumI+496), g(SumI+497), dp) 
    mHx32(i1,i2) = Cmplx( g(SumI+514), g(SumI+515), dp) 
    mHxb32(i1,i2) = Cmplx( g(SumI+532), g(SumI+533), dp) 
    MnuL(i1,i2) = Cmplx( g(SumI+556), g(SumI+557), dp) 
   End Do 
  End Do 

  mue= Cmplx(g(112),g(113),dp) 
  Bmue= Cmplx(g(294),g(295),dp) 
  mHd2= g(404) 
  mHu2= g(405) 
  MassB= Cmplx(g(550),g(551),dp) 
  MassWB= Cmplx(g(552),g(553),dp) 
  MassG= Cmplx(g(554),g(555),dp) 

  Iname = Iname - 1 
 
 End Subroutine GToParameters555


 Subroutine ParametersToG555(g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3,mue,MXM3,MWM3,MGM3  &
  & ,MBM3,TYu,TYd,TYe,TYb3,TYw3,TYx3,Bmue,BMXM3,BMWM3,BMGM3,BMBM3,mq2,ml2,mHd2 &
  & ,mHu2,md2,mu2,me2,mHw32,mHg3p2,mHb32,mHx32,mHxb32,MassB,MassWB,MassG,MnuL,g)

 Implicit None 
  Real(dp), Intent(out) :: g(573) 
  Real(dp), Intent(in) :: g1,g2,g3,mHd2,mHu2

  Complex(dp),Intent(in) :: mue, Bmue, MassB, MassWB, MassG
 
  Complex(dp),Intent(in), Dimension(3,3) :: Yu,Yd,Ye,Yb3,Yw3,Yx3               &
     & ,MXM3,MWM3,MGM3,MBM3,TYu,TYd,TYe,TYb3,TYw3,TYx3,BMXM3,BMWM3,BMGM3,BMBM3 &
     & ,mq2,ml2,md2,mu2,me2,mHw32,mHg3p2,mHb32,mHx32,mHxb32,MnuL

  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'ParametersToG555' 
 
  g(1) = g1  
  g(2) = g2  
  g(3) = g3  

  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3 )
    g(SumI+4) = Real(Yu(i1,i2), dp) 
    g(SumI+5) = Aimag(Yu(i1,i2)) 
    g(SumI+22) = Real(Yd(i1,i2), dp) 
    g(SumI+23) = Aimag(Yd(i1,i2)) 
    g(SumI+40) = Real(Ye(i1,i2), dp) 
    g(SumI+41) = Aimag(Ye(i1,i2)) 
    g(SumI+58) = Real(Yb3(i1,i2), dp) 
    g(SumI+59) = Aimag(Yb3(i1,i2)) 
    g(SumI+76) = Real(Yw3(i1,i2), dp) 
    g(SumI+77) = Aimag(Yw3(i1,i2)) 
    g(SumI+94) = Real(Yx3(i1,i2), dp) 
    g(SumI+95) = Aimag(Yx3(i1,i2)) 
    g(SumI+114) = Real(MXM3(i1,i2), dp) 
    g(SumI+115) = Aimag(MXM3(i1,i2)) 
    g(SumI+132) = Real(MWM3(i1,i2), dp) 
    g(SumI+133) = Aimag(MWM3(i1,i2)) 
    g(SumI+150) = Real(MGM3(i1,i2), dp) 
    g(SumI+151) = Aimag(MGM3(i1,i2)) 
    g(SumI+168) = Real(MBM3(i1,i2), dp) 
    g(SumI+169) = Aimag(MBM3(i1,i2)) 
    g(SumI+186) = Real(TYu(i1,i2), dp) 
    g(SumI+187) = Aimag(TYu(i1,i2)) 
    g(SumI+204) = Real(TYd(i1,i2), dp) 
    g(SumI+205) = Aimag(TYd(i1,i2)) 
    g(SumI+222) = Real(TYe(i1,i2), dp) 
    g(SumI+223) = Aimag(TYe(i1,i2)) 
    g(SumI+240) = Real(TYb3(i1,i2), dp) 
    g(SumI+241) = Aimag(TYb3(i1,i2)) 
    g(SumI+258) = Real(TYw3(i1,i2), dp) 
    g(SumI+259) = Aimag(TYw3(i1,i2)) 
    g(SumI+276) = Real(TYx3(i1,i2), dp) 
    g(SumI+277) = Aimag(TYx3(i1,i2)) 
    g(SumI+296) = Real(BMXM3(i1,i2), dp) 
    g(SumI+297) = Aimag(BMXM3(i1,i2)) 
    g(SumI+314) = Real(BMWM3(i1,i2), dp) 
    g(SumI+315) = Aimag(BMWM3(i1,i2)) 
    g(SumI+332) = Real(BMGM3(i1,i2), dp) 
    g(SumI+333) = Aimag(BMGM3(i1,i2)) 
    g(SumI+350) = Real(BMBM3(i1,i2), dp) 
    g(SumI+351) = Aimag(BMBM3(i1,i2)) 
    g(SumI+368) = Real(mq2(i1,i2), dp) 
    g(SumI+369) = Aimag(mq2(i1,i2)) 
    g(SumI+386) = Real(ml2(i1,i2), dp) 
    g(SumI+387) = Aimag(ml2(i1,i2)) 
    g(SumI+406) = Real(md2(i1,i2), dp) 
    g(SumI+407) = Aimag(md2(i1,i2)) 
    g(SumI+424) = Real(mu2(i1,i2), dp) 
    g(SumI+425) = Aimag(mu2(i1,i2)) 
    g(SumI+442) = Real(me2(i1,i2), dp) 
    g(SumI+443) = Aimag(me2(i1,i2)) 
    g(SumI+460) = Real(mHw32(i1,i2), dp) 
    g(SumI+461) = Aimag(mHw32(i1,i2)) 
    g(SumI+478) = Real(mHg3p2(i1,i2), dp) 
    g(SumI+479) = Aimag(mHg3p2(i1,i2)) 
    g(SumI+496) = Real(mHb32(i1,i2), dp) 
    g(SumI+497) = Aimag(mHb32(i1,i2)) 
    g(SumI+514) = Real(mHx32(i1,i2), dp) 
    g(SumI+515) = Aimag(mHx32(i1,i2)) 
    g(SumI+532) = Real(mHxb32(i1,i2), dp) 
    g(SumI+533) = Aimag(mHxb32(i1,i2)) 
    g(SumI+556) = Real(MnuL(i1,i2), dp) 
    g(SumI+557) = Aimag(MnuL(i1,i2)) 
   End Do 
  End Do 

  g(112) = Real(mue,dp)  
  g(113) = Aimag(mue)
  g(294) = Real(Bmue,dp)  
  g(295) = Aimag(Bmue)  
  g(404) = mHd2  
  g(405) = mHu2
  g(550) = Real(MassB,dp)  
  g(551) = Aimag(MassB)  
  g(552) = Real(MassWB,dp)  
  g(553) = Aimag(MassWB)  
  g(554) = Real(MassG,dp)  
  g(555) = Aimag(MassG)  

  Iname = Iname - 1 
 
 End Subroutine ParametersToG555


Subroutine rge111(len, T, GY, F) 
Implicit None 
Integer, Intent(in) :: len 
Real(dp), Intent(in) :: T, GY(len) 
Real(dp), Intent(out) :: F(len) 
Integer :: i2
Real(dp) :: q 
Real(dp) :: g1,betag11,betag12,Dg1,g2,betag21,betag22,Dg2,g3,betag31,betag32,Dg3
Complex(dp) :: Yu(3,3),betaYu1(3,3),betaYu2(3,3),DYu(3,3),adjYu(3,3),Yd(3,3)          & 
& ,betaYd1(3,3),betaYd2(3,3),DYd(3,3),adjYd(3,3),Ye(3,3),betaYe1(3,3),betaYe2(3,3)       & 
& ,DYe(3,3),adjYe(3,3),Yb3(3,3),betaYb31(3,3),betaYb32(3,3),DYb3(3,3),adjYb3(3,3)        & 
& ,Yw3(3,3),betaYw31(3,3),betaYw32(3,3),DYw3(3,3),adjYw3(3,3),Yx3(3,3),betaYx31(3,3)     & 
& ,betaYx32(3,3),DYx3(3,3),adjYx3(3,3)
Complex(dp) :: Yb3adjYb3(3,3),YdadjYd(3,3),YeadjYe(3,3),YuadjYu(3,3),Yw3adjYw3(3,3),Yx3adjYx3(3,3),  & 
& adjYb3Yb3(3,3),adjYdYd(3,3),adjYeYe(3,3),adjYuYu(3,3),adjYw3Yw3(3,3),adjYx3Yx3(3,3),   & 
& CYdTpYd(3,3),CYx3Yd(3,3),Yb3adjYb3Yb3(3,3),Yb3adjYeYe(3,3),Yb3adjYw3Yw3(3,3),          & 
& YdadjYdYd(3,3),YdadjYuYu(3,3),YeadjYb3Yb3(3,3),YeadjYeYe(3,3),YeadjYw3Yw3(3,3),        & 
& YuadjYdYd(3,3),YuadjYuYu(3,3),Yw3adjYb3Yb3(3,3),Yw3adjYeYe(3,3),Yw3adjYw3Yw3(3,3),     & 
& Yx3adjYx3Yx3(3,3),Yx3CYdTpYd(3,3),TpYx3CYx3Yd(3,3)

Complex(dp) :: YeadjYb3(3,3),YuadjYd(3,3),Yw3adjYb3(3,3),Yw3adjYe(3,3),CYuTpYd(3,3),TpYx3CYx3(3,3),  & 
& adjYb3Yb3adjYb3(3,3),adjYdYdadjYd(3,3),adjYdTpYx3CYx3(3,3),adjYeYeadjYb3(3,3),         & 
& adjYeYeadjYe(3,3),adjYuYuadjYd(3,3),adjYuYuadjYu(3,3),adjYw3Yw3adjYb3(3,3),            & 
& adjYw3Yw3adjYe(3,3),adjYw3Yw3adjYw3(3,3),adjYx3Yx3adjYx3(3,3),TpYdadjYx3Yx3(3,3),      & 
& TpYdCYdTpYd(3,3),TpYuCYuTpYd(3,3),Yb3adjYb3Yb3adjYb3(3,3),Yb3adjYeYeadjYb3(3,3),       & 
& Yb3adjYw3Yw3adjYb3(3,3),YdadjYdYdadjYd(3,3),YdadjYdTpYx3CYx3(3,3),YdadjYuYuadjYd(3,3), & 
& YeadjYeYeadjYe(3,3),YeadjYw3Yw3adjYe(3,3),YuadjYuYuadjYu(3,3),Yw3adjYw3Yw3adjYw3(3,3), & 
& Yx3adjYx3Yx3adjYx3(3,3),adjYb3Yb3adjYb3Yb3(3,3),adjYb3Yb3adjYeYe(3,3),adjYb3Yb3adjYw3Yw3(3,3),& 
& adjYdYdadjYdYd(3,3),adjYdYdadjYuYu(3,3),adjYdTpYx3CYx3Yd(3,3),adjYeYeadjYb3Yb3(3,3),   & 
& adjYeYeadjYeYe(3,3),adjYeYeadjYw3Yw3(3,3),adjYuYuadjYdYd(3,3),adjYuYuadjYuYu(3,3),     & 
& adjYw3Yw3adjYb3Yb3(3,3),adjYw3Yw3adjYeYe(3,3),adjYw3Yw3adjYw3Yw3(3,3),adjYx3Yx3adjYx3Yx3(3,3),& 
& CYdTpYdadjYx3Yx3(3,3),CYdTpYdCYdTpYd(3,3),CYdTpYuCYuTpYd(3,3),CYx3TpYx3CYx3Yd(3,3),    & 
& Yb3adjYb3Yb3adjYb3Yb3(3,3),Yb3adjYb3Yb3adjYw3Yw3(3,3),Yb3adjYeYeadjYb3Yb3(3,3),        & 
& Yb3adjYeYeadjYeYe(3,3),Yb3adjYw3Yw3adjYb3Yb3(3,3),Yb3adjYw3Yw3adjYw3Yw3(3,3),          & 
& YdadjYdYdadjYdYd(3,3),YdadjYdTpYx3CYx3Yd(3,3),YdadjYuYuadjYdYd(3,3),YdadjYuYuadjYuYu(3,3),& 
& YeadjYb3Yb3adjYb3Yb3(3,3),YeadjYb3Yb3adjYeYe(3,3),YeadjYb3Yb3adjYw3Yw3(3,3),           & 
& YeadjYeYeadjYeYe(3,3),YeadjYw3Yw3adjYb3Yb3(3,3),YeadjYw3Yw3adjYeYe(3,3),               & 
& YeadjYw3Yw3adjYw3Yw3(3,3),YuadjYdYdadjYdYd(3,3),YuadjYdYdadjYuYu(3,3),YuadjYdTpYx3CYx3Yd(3,3),& 
& YuadjYuYuadjYuYu(3,3),Yw3adjYb3Yb3adjYb3Yb3(3,3),Yw3adjYb3Yb3adjYw3Yw3(3,3),           & 
& Yw3adjYeYeadjYeYe(3,3),Yw3adjYeYeadjYw3Yw3(3,3),Yw3adjYw3Yw3adjYb3Yb3(3,3),            & 
& Yw3adjYw3Yw3adjYw3Yw3(3,3),Yx3adjYx3Yx3adjYx3Yx3(3,3),Yx3CYdTpYdadjYx3Yx3(3,3),        & 
& Yx3CYdTpYdCYdTpYd(3,3),Yx3CYdTpYuCYuTpYd(3,3),TpYx3CYx3TpYx3CYx3Yd(3,3)

Complex(dp) :: TrYb3adjYb3,TrYdadjYd,TrYeadjYe,TrYuadjYu,TrYw3adjYw3,TrYx3adjYx3

Complex(dp) :: TrYb3adjYb3Yb3adjYb3,TrYb3adjYeYeadjYb3,TrYb3adjYw3Yw3adjYb3,TrYdadjYdYdadjYd,        & 
& TrYdadjYdTpYx3CYx3,TrYdadjYuYuadjYd,TrYeadjYeYeadjYe,TrYeadjYw3Yw3adjYe,               & 
& TrYuadjYuYuadjYu,TrYw3adjYw3Yw3adjYw3,TrYx3adjYx3Yx3adjYx3

Real(dp) :: g1p2,g1p3,g2p2,g2p3,g3p2,g3p3

Real(dp) :: g1p4,g2p4,g3p4

Iname = Iname +1 
NameOfUnit(Iname) = 'rge111' 
 
OnlyDiagonal = .Not.GenerationMixing 
q = t 
 
Call GToParameters111(gy,g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3)

If (ThresholdCrossed.Lt.1) Then 
Yx3(1,:) = 0._dp 
Yb3(1,:) = 0._dp 
Yw3(1,:) = 0._dp 
End If 

If (ThresholdCrossed.Lt.2) Then 
Yx3(2,:) = 0._dp 
Yb3(2,:) = 0._dp 
Yw3(2,:) = 0._dp 
End If 

If (ThresholdCrossed.Lt.3) Then 
Yx3(3,:) = 0._dp 
Yb3(3,:) = 0._dp 
Yw3(3,:) = 0._dp 
End If 

Call Adjungate(Yu,adjYu)
Call Adjungate(Yd,adjYd)
Call Adjungate(Ye,adjYe)
Call Adjungate(Yb3,adjYb3)
Call Adjungate(Yw3,adjYw3)
Call Adjungate(Yx3,adjYx3)
 Yb3adjYb3 = Matmul(Yb3,adjYb3) 
Forall(i2=1:3)  Yb3adjYb3(i2,i2) =  Real(Yb3adjYb3(i2,i2),dp) 
 YdadjYd = Matmul(Yd,adjYd) 
Forall(i2=1:3)  YdadjYd(i2,i2) =  Real(YdadjYd(i2,i2),dp) 
 YeadjYe = Matmul(Ye,adjYe) 
Forall(i2=1:3)  YeadjYe(i2,i2) =  Real(YeadjYe(i2,i2),dp) 
 YuadjYu = Matmul(Yu,adjYu) 
Forall(i2=1:3)  YuadjYu(i2,i2) =  Real(YuadjYu(i2,i2),dp) 
 Yw3adjYw3 = Matmul(Yw3,adjYw3) 
Forall(i2=1:3)  Yw3adjYw3(i2,i2) =  Real(Yw3adjYw3(i2,i2),dp) 
 Yx3adjYx3 = Matmul(Yx3,adjYx3) 
Forall(i2=1:3)  Yx3adjYx3(i2,i2) =  Real(Yx3adjYx3(i2,i2),dp) 
 adjYb3Yb3 = Matmul(adjYb3,Yb3) 
Forall(i2=1:3)  adjYb3Yb3(i2,i2) =  Real(adjYb3Yb3(i2,i2),dp) 
 adjYdYd = Matmul(adjYd,Yd) 
Forall(i2=1:3)  adjYdYd(i2,i2) =  Real(adjYdYd(i2,i2),dp) 
 adjYeYe = Matmul(adjYe,Ye) 
Forall(i2=1:3)  adjYeYe(i2,i2) =  Real(adjYeYe(i2,i2),dp) 
 adjYuYu = Matmul(adjYu,Yu) 
Forall(i2=1:3)  adjYuYu(i2,i2) =  Real(adjYuYu(i2,i2),dp) 
 adjYw3Yw3 = Matmul(adjYw3,Yw3) 
Forall(i2=1:3)  adjYw3Yw3(i2,i2) =  Real(adjYw3Yw3(i2,i2),dp) 
 adjYx3Yx3 = Matmul(adjYx3,Yx3) 
Forall(i2=1:3)  adjYx3Yx3(i2,i2) =  Real(adjYx3Yx3(i2,i2),dp) 
 CYdTpYd = Matmul(Conjg(Yd),Transpose(Yd)) 
Forall(i2=1:3)  CYdTpYd(i2,i2) =  Real(CYdTpYd(i2,i2),dp) 
 CYx3Yd = Matmul(Conjg(Yx3),Yd) 
 Yb3adjYb3Yb3 = Matmul(Yb3,adjYb3Yb3) 
 Yb3adjYeYe = Matmul(Yb3,adjYeYe) 
 Yb3adjYw3Yw3 = Matmul(Yb3,adjYw3Yw3) 
 YdadjYdYd = Matmul(Yd,adjYdYd) 
 YdadjYuYu = Matmul(Yd,adjYuYu) 
 YeadjYb3Yb3 = Matmul(Ye,adjYb3Yb3) 
 YeadjYeYe = Matmul(Ye,adjYeYe) 
 YeadjYw3Yw3 = Matmul(Ye,adjYw3Yw3) 
 YuadjYdYd = Matmul(Yu,adjYdYd) 
 YuadjYuYu = Matmul(Yu,adjYuYu) 
 Yw3adjYb3Yb3 = Matmul(Yw3,adjYb3Yb3) 
 Yw3adjYeYe = Matmul(Yw3,adjYeYe) 
 Yw3adjYw3Yw3 = Matmul(Yw3,adjYw3Yw3) 
 Yx3adjYx3Yx3 = Matmul(Yx3,adjYx3Yx3) 
 Yx3CYdTpYd = Matmul(Yx3,CYdTpYd) 
 TpYx3CYx3Yd = Matmul(Transpose(Yx3),CYx3Yd) 
 TrYb3adjYb3 = Real(cTrace(Yb3adjYb3),dp) 
 TrYdadjYd = Real(cTrace(YdadjYd),dp) 
 TrYeadjYe = Real(cTrace(YeadjYe),dp) 
 TrYuadjYu = Real(cTrace(YuadjYu),dp) 
 TrYw3adjYw3 = Real(cTrace(Yw3adjYw3),dp) 
 TrYx3adjYx3 = Real(cTrace(Yx3adjYx3),dp) 
 g1p2 =g1**2 
 g1p3 =g1**3 
 g2p2 =g2**2 
 g2p3 =g2**3 
 g3p2 =g3**2 
 g3p3 =g3**3 


If (TwoLoopRGE) Then 
 YeadjYb3 = Matmul(Ye,adjYb3) 
 YuadjYd = Matmul(Yu,adjYd) 
 Yw3adjYb3 = Matmul(Yw3,adjYb3) 
 Yw3adjYe = Matmul(Yw3,adjYe) 
 CYuTpYd = Matmul(Conjg(Yu),Transpose(Yd)) 
 TpYx3CYx3 = Matmul(Transpose(Yx3),Conjg(Yx3)) 
Forall(i2=1:3)  TpYx3CYx3(i2,i2) =  Real(TpYx3CYx3(i2,i2),dp) 
 adjYb3Yb3adjYb3 = Matmul(adjYb3,Yb3adjYb3) 
 adjYdYdadjYd = Matmul(adjYd,YdadjYd) 
 adjYdTpYx3CYx3 = Matmul(adjYd,TpYx3CYx3) 
 adjYeYeadjYb3 = Matmul(adjYe,YeadjYb3) 
 adjYeYeadjYe = Matmul(adjYe,YeadjYe) 
 adjYuYuadjYd = Matmul(adjYu,YuadjYd) 
 adjYuYuadjYu = Matmul(adjYu,YuadjYu) 
 adjYw3Yw3adjYb3 = Matmul(adjYw3,Yw3adjYb3) 
 adjYw3Yw3adjYe = Matmul(adjYw3,Yw3adjYe) 
 adjYw3Yw3adjYw3 = Matmul(adjYw3,Yw3adjYw3) 
 adjYx3Yx3adjYx3 = Matmul(adjYx3,Yx3adjYx3) 
 TpYdadjYx3Yx3 = Matmul(Transpose(Yd),adjYx3Yx3) 
 TpYdCYdTpYd = Matmul(Transpose(Yd),CYdTpYd) 
 TpYuCYuTpYd = Matmul(Transpose(Yu),CYuTpYd) 
 Yb3adjYb3Yb3adjYb3 = Matmul(Yb3,adjYb3Yb3adjYb3) 
Forall(i2=1:3)  Yb3adjYb3Yb3adjYb3(i2,i2) =  Real(Yb3adjYb3Yb3adjYb3(i2,i2),dp) 
 Yb3adjYeYeadjYb3 = Matmul(Yb3,adjYeYeadjYb3) 
Forall(i2=1:3)  Yb3adjYeYeadjYb3(i2,i2) =  Real(Yb3adjYeYeadjYb3(i2,i2),dp) 
 Yb3adjYw3Yw3adjYb3 = Matmul(Yb3,adjYw3Yw3adjYb3) 
Forall(i2=1:3)  Yb3adjYw3Yw3adjYb3(i2,i2) =  Real(Yb3adjYw3Yw3adjYb3(i2,i2),dp) 
 YdadjYdYdadjYd = Matmul(Yd,adjYdYdadjYd) 
Forall(i2=1:3)  YdadjYdYdadjYd(i2,i2) =  Real(YdadjYdYdadjYd(i2,i2),dp) 
 YdadjYdTpYx3CYx3 = Matmul(Yd,adjYdTpYx3CYx3) 
 YdadjYuYuadjYd = Matmul(Yd,adjYuYuadjYd) 
Forall(i2=1:3)  YdadjYuYuadjYd(i2,i2) =  Real(YdadjYuYuadjYd(i2,i2),dp) 
 YeadjYeYeadjYe = Matmul(Ye,adjYeYeadjYe) 
Forall(i2=1:3)  YeadjYeYeadjYe(i2,i2) =  Real(YeadjYeYeadjYe(i2,i2),dp) 
 YeadjYw3Yw3adjYe = Matmul(Ye,adjYw3Yw3adjYe) 
Forall(i2=1:3)  YeadjYw3Yw3adjYe(i2,i2) =  Real(YeadjYw3Yw3adjYe(i2,i2),dp) 
 YuadjYuYuadjYu = Matmul(Yu,adjYuYuadjYu) 
Forall(i2=1:3)  YuadjYuYuadjYu(i2,i2) =  Real(YuadjYuYuadjYu(i2,i2),dp) 
 Yw3adjYw3Yw3adjYw3 = Matmul(Yw3,adjYw3Yw3adjYw3) 
Forall(i2=1:3)  Yw3adjYw3Yw3adjYw3(i2,i2) =  Real(Yw3adjYw3Yw3adjYw3(i2,i2),dp) 
 Yx3adjYx3Yx3adjYx3 = Matmul(Yx3,adjYx3Yx3adjYx3) 
Forall(i2=1:3)  Yx3adjYx3Yx3adjYx3(i2,i2) =  Real(Yx3adjYx3Yx3adjYx3(i2,i2),dp) 
 adjYb3Yb3adjYb3Yb3 = Matmul(adjYb3,Yb3adjYb3Yb3) 
Forall(i2=1:3)  adjYb3Yb3adjYb3Yb3(i2,i2) =  Real(adjYb3Yb3adjYb3Yb3(i2,i2),dp) 
 adjYb3Yb3adjYeYe = Matmul(adjYb3,Yb3adjYeYe) 
 adjYb3Yb3adjYw3Yw3 = Matmul(adjYb3,Yb3adjYw3Yw3) 
 adjYdYdadjYdYd = Matmul(adjYd,YdadjYdYd) 
Forall(i2=1:3)  adjYdYdadjYdYd(i2,i2) =  Real(adjYdYdadjYdYd(i2,i2),dp) 
 adjYdYdadjYuYu = Matmul(adjYd,YdadjYuYu) 
 adjYdTpYx3CYx3Yd = Matmul(adjYd,TpYx3CYx3Yd) 
Forall(i2=1:3)  adjYdTpYx3CYx3Yd(i2,i2) =  Real(adjYdTpYx3CYx3Yd(i2,i2),dp) 
 adjYeYeadjYb3Yb3 = Matmul(adjYe,YeadjYb3Yb3) 
 adjYeYeadjYeYe = Matmul(adjYe,YeadjYeYe) 
Forall(i2=1:3)  adjYeYeadjYeYe(i2,i2) =  Real(adjYeYeadjYeYe(i2,i2),dp) 
 adjYeYeadjYw3Yw3 = Matmul(adjYe,YeadjYw3Yw3) 
 adjYuYuadjYdYd = Matmul(adjYu,YuadjYdYd) 
 adjYuYuadjYuYu = Matmul(adjYu,YuadjYuYu) 
Forall(i2=1:3)  adjYuYuadjYuYu(i2,i2) =  Real(adjYuYuadjYuYu(i2,i2),dp) 
 adjYw3Yw3adjYb3Yb3 = Matmul(adjYw3,Yw3adjYb3Yb3) 
 adjYw3Yw3adjYeYe = Matmul(adjYw3,Yw3adjYeYe) 
 adjYw3Yw3adjYw3Yw3 = Matmul(adjYw3,Yw3adjYw3Yw3) 
Forall(i2=1:3)  adjYw3Yw3adjYw3Yw3(i2,i2) =  Real(adjYw3Yw3adjYw3Yw3(i2,i2),dp) 
 adjYx3Yx3adjYx3Yx3 = Matmul(adjYx3,Yx3adjYx3Yx3) 
Forall(i2=1:3)  adjYx3Yx3adjYx3Yx3(i2,i2) =  Real(adjYx3Yx3adjYx3Yx3(i2,i2),dp) 
 CYdTpYdadjYx3Yx3 = Matmul(Conjg(Yd),TpYdadjYx3Yx3) 
 CYdTpYdCYdTpYd = Matmul(Conjg(Yd),TpYdCYdTpYd) 
Forall(i2=1:3)  CYdTpYdCYdTpYd(i2,i2) =  Real(CYdTpYdCYdTpYd(i2,i2),dp) 
 CYdTpYuCYuTpYd = Matmul(Conjg(Yd),TpYuCYuTpYd) 
Forall(i2=1:3)  CYdTpYuCYuTpYd(i2,i2) =  Real(CYdTpYuCYuTpYd(i2,i2),dp) 
 CYx3TpYx3CYx3Yd = Matmul(Conjg(Yx3),TpYx3CYx3Yd) 
 Yb3adjYb3Yb3adjYb3Yb3 = Matmul(Yb3,adjYb3Yb3adjYb3Yb3) 
 Yb3adjYb3Yb3adjYw3Yw3 = Matmul(Yb3,adjYb3Yb3adjYw3Yw3) 
 Yb3adjYeYeadjYb3Yb3 = Matmul(Yb3,adjYeYeadjYb3Yb3) 
 Yb3adjYeYeadjYeYe = Matmul(Yb3,adjYeYeadjYeYe) 
 Yb3adjYw3Yw3adjYb3Yb3 = Matmul(Yb3,adjYw3Yw3adjYb3Yb3) 
 Yb3adjYw3Yw3adjYw3Yw3 = Matmul(Yb3,adjYw3Yw3adjYw3Yw3) 
 YdadjYdYdadjYdYd = Matmul(Yd,adjYdYdadjYdYd) 
 YdadjYdTpYx3CYx3Yd = Matmul(Yd,adjYdTpYx3CYx3Yd) 
 YdadjYuYuadjYdYd = Matmul(Yd,adjYuYuadjYdYd) 
 YdadjYuYuadjYuYu = Matmul(Yd,adjYuYuadjYuYu) 
 YeadjYb3Yb3adjYb3Yb3 = Matmul(Ye,adjYb3Yb3adjYb3Yb3) 
 YeadjYb3Yb3adjYeYe = Matmul(Ye,adjYb3Yb3adjYeYe) 
 YeadjYb3Yb3adjYw3Yw3 = Matmul(Ye,adjYb3Yb3adjYw3Yw3) 
 YeadjYeYeadjYeYe = Matmul(Ye,adjYeYeadjYeYe) 
 YeadjYw3Yw3adjYb3Yb3 = Matmul(Ye,adjYw3Yw3adjYb3Yb3) 
 YeadjYw3Yw3adjYeYe = Matmul(Ye,adjYw3Yw3adjYeYe) 
 YeadjYw3Yw3adjYw3Yw3 = Matmul(Ye,adjYw3Yw3adjYw3Yw3) 
 YuadjYdYdadjYdYd = Matmul(Yu,adjYdYdadjYdYd) 
 YuadjYdYdadjYuYu = Matmul(Yu,adjYdYdadjYuYu) 
 YuadjYdTpYx3CYx3Yd = Matmul(Yu,adjYdTpYx3CYx3Yd) 
 YuadjYuYuadjYuYu = Matmul(Yu,adjYuYuadjYuYu) 
 Yw3adjYb3Yb3adjYb3Yb3 = Matmul(Yw3,adjYb3Yb3adjYb3Yb3) 
 Yw3adjYb3Yb3adjYw3Yw3 = Matmul(Yw3,adjYb3Yb3adjYw3Yw3) 
 Yw3adjYeYeadjYeYe = Matmul(Yw3,adjYeYeadjYeYe) 
 Yw3adjYeYeadjYw3Yw3 = Matmul(Yw3,adjYeYeadjYw3Yw3) 
 Yw3adjYw3Yw3adjYb3Yb3 = Matmul(Yw3,adjYw3Yw3adjYb3Yb3) 
 Yw3adjYw3Yw3adjYw3Yw3 = Matmul(Yw3,adjYw3Yw3adjYw3Yw3) 
 Yx3adjYx3Yx3adjYx3Yx3 = Matmul(Yx3,adjYx3Yx3adjYx3Yx3) 
 Yx3CYdTpYdadjYx3Yx3 = Matmul(Yx3,CYdTpYdadjYx3Yx3) 
 Yx3CYdTpYdCYdTpYd = Matmul(Yx3,CYdTpYdCYdTpYd) 
 Yx3CYdTpYuCYuTpYd = Matmul(Yx3,CYdTpYuCYuTpYd) 
 TpYx3CYx3TpYx3CYx3Yd = Matmul(Transpose(Yx3),CYx3TpYx3CYx3Yd) 
 TrYb3adjYb3Yb3adjYb3 = cTrace(Yb3adjYb3Yb3adjYb3) 
 TrYb3adjYeYeadjYb3 = cTrace(Yb3adjYeYeadjYb3) 
 TrYb3adjYw3Yw3adjYb3 = cTrace(Yb3adjYw3Yw3adjYb3) 
 TrYdadjYdYdadjYd = cTrace(YdadjYdYdadjYd) 
 TrYdadjYdTpYx3CYx3 = cTrace(YdadjYdTpYx3CYx3) 
 TrYdadjYuYuadjYd = cTrace(YdadjYuYuadjYd) 
 TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe) 
 TrYeadjYw3Yw3adjYe = cTrace(YeadjYw3Yw3adjYe) 
 TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu) 
 TrYw3adjYw3Yw3adjYw3 = cTrace(Yw3adjYw3Yw3adjYw3) 
 TrYx3adjYx3Yx3adjYx3 = cTrace(Yx3adjYx3Yx3adjYx3) 
 g1p4 =g1**4 
 g2p4 =g2**4 
 g3p4 =g3**4 
End If 
 
 
!-------------------- 
! g1 
!-------------------- 
 
betag11  = (g1p3*(66 + 25*NGHx3 + 25*NGHxb3))/10._dp

 
 
If (TwoLoopRGE) Then 
betag12 = (g1p3*(6*(199._dp*(g1p2) + 135._dp*(g2p2) + 440._dp*(g3p2) - 9._dp*(TrYb3adjYb3) -    & 
&  70._dp*(TrYdadjYd) - 90._dp*(TrYeadjYe) - 130._dp*(TrYuadjYu) - 60._dp*(TrYw3adjYw3) -& 
&  190._dp*(TrYx3adjYx3)) + 125*(5._dp*(g1p2) + 9._dp*(g2p2) + 16._dp*(g3p2))*(NGHx3 +& 
&  NGHxb3)))/150._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21  = (g2p3*(2 + 4*NGHw3 + 3*NGHx3 +               & 
&  3*NGHxb3))/2._dp

 
 
If (TwoLoopRGE) Then 
betag22 = (g2p3*(2*(27._dp*(g1p2) + 375._dp*(g2p2) + 360._dp*(g3p2) - 9._dp*(TrYb3adjYb3) -     & 
&  90._dp*(TrYdadjYd) - 30._dp*(TrYeadjYe) - 90._dp*(TrYuadjYu) - 140._dp*(TrYw3adjYw3) -& 
&  90._dp*(TrYx3adjYx3)) + 720*g2p2*NGHw3 + 15*(5._dp*(g1p2) +            & 
&  21._dp*(g2p2) + 16._dp*(g3p2))*(NGHx3 + NGHxb3)))/30._dp

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31  = g3p3*(-3 + 3*NGHg3 + NGHx3 + NGHxb3)

 
 
If (TwoLoopRGE) Then 
betag32 = (g3p3*(3*(11._dp*(g1p2) + 45._dp*(g2p2) + 70._dp*(g3p2) - 20._dp*(TrYdadjYd) -        & 
&  20._dp*(TrYuadjYu) - 20._dp*(TrYx3adjYx3)) + 810*g3p2*NGHg3 +          & 
&  5*(5._dp*(g1p2) + 9._dp*(g2p2) + 34._dp*(g3p2))*(NGHx3 +               & 
&  NGHxb3)))/15._dp

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yu 
!-------------------- 
 
betaYu1  = ((-26._dp*(g1p2) - 90._dp*(g2p2) - 160._dp*(g3p2) + 9._dp*(TrYb3adjYb3)    & 
&  + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3) + 90._dp*(TrYx3adjYx3))*Yu)               & 
& /30._dp + YuadjYdYd + 3._dp*(YuadjYuYu)

 
 
If (TwoLoopRGE) Then 
betaYu2 = (2._dp*(g1p2)/5._dp - 3._dp*(TrYdadjYd) - TrYeadjYe)*YuadjYdYd + ((4._dp*(g1p2) +     & 
&  60._dp*(g2p2) - 9._dp*(TrYb3adjYb3) - 90._dp*(TrYuadjYu) - 45._dp*(TrYw3adjYw3) -     & 
&  90._dp*(TrYx3adjYx3))*YuadjYuYu)/10._dp - 2*(YuadjYdTpYx3CYx3Yd + YuadjYdYdadjYdYd +  & 
&  YuadjYdYdadjYuYu + 2._dp*(YuadjYuYuadjYuYu)) + (Yu*(-486._dp*(TrYb3adjYb3Yb3adjYb3) - & 
&  540._dp*(TrYb3adjYeYeadjYb3) - 2565._dp*(TrYb3adjYw3Yw3adjYb3) + 2*(5486._dp*(g1p4) + & 
&  20*g1p2*(45._dp*(g2p2) + 136._dp*(g3p2)) + 50*(135._dp*(g2p4) + 144*g2p2*g3p2 -       & 
&  32._dp*(g3p4)) - 5400._dp*(TrYdadjYdTpYx3CYx3) - 2700._dp*(TrYdadjYuYuadjYd) -        & 
&  1350._dp*(TrYeadjYw3Yw3adjYe) - 8100._dp*(TrYuadjYuYuadjYu) - 3375._dp*(TrYw3adjYw3Yw3adjYw3) -& 
&  8100._dp*(TrYx3adjYx3Yx3adjYx3)) + 60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) +& 
&  65*g1p2*(NGHx3 + NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu +   & 
&  TrYx3adjYx3) + g3p2*(3*NGHg3 + NGHx3 + NGHxb3)) +& 
&  2700*g2p2*(4._dp*(TrYw3adjYw3) + g2p2*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))))/1800._dp

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
betaYd1  = 2._dp*(TpYx3CYx3Yd) + (-7._dp*(g1p2)/15._dp - 3._dp*(g2p2) -               & 
&  16._dp*(g3p2)/3._dp + 3._dp*(TrYdadjYd) + TrYeadjYe)*Yd + 3._dp*(YdadjYdYd)           & 
&  + YdadjYuYu

 
 
If (TwoLoopRGE) Then 
betaYd2 = (4._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 9._dp*(TrYdadjYd) - 3._dp*(TrYeadjYe))*YdadjYdYd +& 
&  (2*TpYx3CYx3Yd*(-3._dp*(TrYb3adjYb3) + 5*(2._dp*(g1p2) + 6._dp*(g2p2) -               & 
&  6._dp*(TrYuadjYu) - 3._dp*(TrYw3adjYw3) - 6._dp*(TrYx3adjYx3))) + (8._dp*(g1p2) -     & 
&  3._dp*(TrYb3adjYb3) - 30._dp*(TrYuadjYu) - 15._dp*(TrYw3adjYw3) - 30._dp*(TrYx3adjYx3))*YdadjYuYu)/10._dp -& 
&  2*(TpYx3CYx3TpYx3CYx3Yd + YdadjYdTpYx3CYx3Yd + 2._dp*(YdadjYdYdadjYdYd) +             & 
&  YdadjYuYuadjYdYd + YdadjYuYuadjYuYu) + (Yd*(287._dp*(g1p4) + 90*g1p2*g2p2 +           & 
&  675._dp*(g2p4) + 80*(g1p2 + 9._dp*(g2p2))*g3p2 - 160._dp*(g3p4) - 27._dp*(TrYb3adjYeYeadjYb3) -& 
&  540._dp*(TrYdadjYdTpYx3CYx3) - 810._dp*(TrYdadjYdYdadjYd) - 270._dp*(TrYdadjYuYuadjYd) -& 
&  270._dp*(TrYeadjYeYeadjYe) - 135._dp*(TrYeadjYw3Yw3adjYe) + 135*g2p4*(4*NGHw3 +& 
&  3*(NGHx3 + NGHxb3)) + 3*g1p2*(-12*(TrYdadjYd -          & 
&  3._dp*(TrYeadjYe)) + 35*g1p2*(NGHx3 + NGHxb3)) +        & 
&  480*g3p2*(3._dp*(TrYdadjYd) + g3p2*(3*NGHg3 + NGHx3 +   & 
&  NGHxb3))))/90._dp

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
!-------------------- 
! Ye 
!-------------------- 
 
betaYe1  = (-9._dp*(g1p2)/5._dp - 3._dp*(g2p2) + 3._dp*(TrYdadjYd) + TrYeadjYe)       & 
& *Ye + (3*(YeadjYb3Yb3 + 5*(2._dp*(YeadjYeYe) + YeadjYw3Yw3)))/10._dp

 
 
If (TwoLoopRGE) Then 
betaYe2 = (-36._dp*(YeadjYb3Yb3adjYb3Yb3) - 18*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) +            & 
&  TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))*(YeadjYb3Yb3 + 5._dp*(YeadjYw3Yw3)) -             & 
&  600*((-2._dp*(g2p2) + 3._dp*(TrYdadjYd) + TrYeadjYe)*YeadjYeYe - 2*g2p2*YeadjYw3Yw3) -& 
&  5*(24._dp*(YeadjYb3Yb3adjYeYe) + 12._dp*(YeadjYb3Yb3adjYw3Yw3) + 160._dp*(YeadjYeYeadjYeYe) +& 
&  9._dp*(YeadjYw3Yw3adjYb3Yb3) + 60*(2._dp*(YeadjYw3Yw3adjYeYe) + YeadjYw3Yw3adjYw3Yw3)) +& 
&  20*Ye*(-4*(g1p2 - 40._dp*(g3p2))*TrYdadjYd + 12*g1p2*TrYeadjYe + 3*(45._dp*(g1p4) +   & 
&  6*g1p2*g2p2 + 25._dp*(g2p4) - TrYb3adjYeYeadjYb3 - 20._dp*(TrYdadjYdTpYx3CYx3) -      & 
&  30._dp*(TrYdadjYdYdadjYd) - 10._dp*(TrYdadjYuYuadjYd) - 10._dp*(TrYeadjYeYeadjYe) -   & 
&  5._dp*(TrYeadjYw3Yw3adjYe)) + 45*g1p4*NGHx3 + 45*g1p4*NGHxb3 +& 
&  15*g2p4*(4*NGHw3 + 3*(NGHx3 + NGHxb3))))/200._dp

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
!-------------------- 
! Yb3 
!-------------------- 
 
betaYb31  = (3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) -              & 
&  2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))*Yb3)/10._dp + 9._dp*(Yb3adjYb3Yb3)     & 
& /10._dp + Yb3adjYeYe + 3._dp*(Yb3adjYw3Yw3)/2._dp

 
 
If (TwoLoopRGE) Then 
betaYb32 = (18*(4._dp*(g1p2) + 20._dp*(g2p2) - 3._dp*(TrYb3adjYb3) - 30._dp*(TrYuadjYu) -        & 
&  15._dp*(TrYw3adjYw3) - 30._dp*(TrYx3adjYx3))*Yb3adjYb3Yb3 - 72._dp*(Yb3adjYb3Yb3adjYb3Yb3) +& 
&  40*(6._dp*(g1p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*Yb3adjYeYe -               & 
&  30*(3._dp*(TrYb3adjYb3) + 5*(-8._dp*(g2p2) + 6._dp*(TrYuadjYu) + 3._dp*(TrYw3adjYw3) +& 
&  6._dp*(TrYx3adjYx3)))*Yb3adjYw3Yw3 - 5*(12._dp*(Yb3adjYb3Yb3adjYw3Yw3) +              & 
&  24._dp*(Yb3adjYeYeadjYb3Yb3) + 80._dp*(Yb3adjYeYeadjYeYe) + 45._dp*(Yb3adjYw3Yw3adjYb3Yb3) +& 
&  60._dp*(Yb3adjYw3Yw3adjYw3Yw3)) + Yb3*(3*(276._dp*(g1p4) + 120*g1p2*g2p2 +            & 
&  500._dp*(g2p4) - 18._dp*(TrYb3adjYb3Yb3adjYb3) - 20._dp*(TrYb3adjYeYeadjYb3) -        & 
&  95._dp*(TrYb3adjYw3Yw3adjYb3) - 400._dp*(TrYdadjYdTpYx3CYx3) - 200._dp*(TrYdadjYuYuadjYd) -& 
&  100._dp*(TrYeadjYw3Yw3adjYe) - 600._dp*(TrYuadjYuYuadjYu) - 250._dp*(TrYw3adjYw3Yw3adjYw3) -& 
&  600._dp*(TrYx3adjYx3Yx3adjYx3)) + 20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu +            & 
&  20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 + 15*g1p4*NGHx3 + 15*g1p4*NGHxb3) +& 
&  300*g2p2*(4._dp*(TrYw3adjYw3) + g2p2*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))))/200._dp

 
DYb3 = oo16pi2*( betaYb31 + oo16pi2 * betaYb32 ) 

 
Else 
DYb3 = oo16pi2* betaYb31 
End If 
 
 
!-------------------- 
! Yw3 
!-------------------- 
 
betaYw31  = (-3._dp*(g1p2)/5._dp - 7._dp*(g2p2) + 3._dp*(TrYb3adjYb3)/10._dp +        & 
&  3._dp*(TrYuadjYu) + 3._dp*(TrYw3adjYw3)/2._dp + 3._dp*(TrYx3adjYx3))*Yw3 +            & 
&  3._dp*(Yw3adjYb3Yb3)/10._dp + Yw3adjYeYe + 5._dp*(Yw3adjYw3Yw3)/2._dp

 
 
If (TwoLoopRGE) Then 
betaYw32 = (-36._dp*(Yw3adjYb3Yb3adjYb3Yb3) + 40*(6._dp*(g1p2) - 15._dp*(TrYdadjYd) -            & 
&  5._dp*(TrYeadjYe))*Yw3adjYeYe + 40*(3._dp*(g1p2) + 25._dp*(g2p2))*Yw3adjYw3Yw3 -      & 
&  6*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) + TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))*(3._dp*(Yw3adjYb3Yb3) +& 
&  25._dp*(Yw3adjYw3Yw3)) - 5*(24._dp*(Yw3adjYb3Yb3adjYw3Yw3) + 80._dp*(Yw3adjYeYeadjYeYe) +& 
&  40._dp*(Yw3adjYeYeadjYw3Yw3) + 9._dp*(Yw3adjYw3Yw3adjYb3Yb3) + 120._dp*(Yw3adjYw3Yw3adjYw3Yw3)) +& 
&  Yw3*(828._dp*(g1p4) + 360*g1p2*g2p2 + 5500._dp*(g2p4) - 54._dp*(TrYb3adjYb3Yb3adjYb3) -& 
&  60._dp*(TrYb3adjYeYeadjYb3) - 285._dp*(TrYb3adjYw3Yw3adjYb3) - 1200._dp*(TrYdadjYdTpYx3CYx3) -& 
&  600._dp*(TrYdadjYuYuadjYd) - 300._dp*(TrYeadjYw3Yw3adjYe) - 1800._dp*(TrYuadjYuYuadjYu) +& 
&  1200*g2p2*TrYw3adjYw3 - 750._dp*(TrYw3adjYw3Yw3adjYw3) - 1800._dp*(TrYx3adjYx3Yx3adjYx3) +& 
&  20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu + 20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 +       & 
&  15*g1p4*NGHx3 + 15*g1p4*NGHxb3) + 700*g2p4*(4*NGHw3 +& 
&  3*(NGHx3 + NGHxb3))))/200._dp

 
DYw3 = oo16pi2*( betaYw31 + oo16pi2 * betaYw32 ) 

 
Else 
DYw3 = oo16pi2* betaYw31 
End If 
 
 
!-------------------- 
! Yx3 
!-------------------- 
 
betaYx31  = ((-38._dp*(g1p2) - 90._dp*(g2p2) - 160._dp*(g3p2) + 9._dp*(TrYb3adjYb3)   & 
&  + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3) + 90._dp*(TrYx3adjYx3))*Yx3)              & 
& /30._dp + 3._dp*(Yx3adjYx3Yx3) + 2._dp*(Yx3CYdTpYd)

 
 
If (TwoLoopRGE) Then 
betaYx32 = (8._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 9._dp*(TrYb3adjYb3)/10._dp - 9._dp*(TrYuadjYu) - & 
&  9._dp*(TrYw3adjYw3)/2._dp - 9._dp*(TrYx3adjYx3))*Yx3adjYx3Yx3 + (2*(g1p2 +            & 
&  15._dp*(g2p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*Yx3CYdTpYd)/5._dp -           & 
&  2*(2._dp*(Yx3adjYx3Yx3adjYx3Yx3) + Yx3CYdTpYdadjYx3Yx3 + Yx3CYdTpYdCYdTpYd +          & 
&  Yx3CYdTpYuCYuTpYd) + (Yx3*(-486._dp*(TrYb3adjYb3Yb3adjYb3) - 540._dp*(TrYb3adjYeYeadjYb3) -& 
&  2565._dp*(TrYb3adjYw3Yw3adjYb3) + 2*(8246._dp*(g1p4) + 20*g1p2*(153._dp*(g2p2) +      & 
&  232._dp*(g3p2)) + 50*(135._dp*(g2p4) + 144*g2p2*g3p2 - 32._dp*(g3p4)) -               & 
&  5400._dp*(TrYdadjYdTpYx3CYx3) - 2700._dp*(TrYdadjYuYuadjYd) - 1350._dp*(TrYeadjYw3Yw3adjYe) -& 
&  8100._dp*(TrYuadjYuYuadjYu) - 3375._dp*(TrYw3adjYw3Yw3adjYw3) - 8100._dp*(TrYx3adjYx3Yx3adjYx3)) +& 
&  60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) + 95*g1p2*(NGHx3 +  & 
&  NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu + TrYx3adjYx3) + g3p2*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)) + 2700*g2p2*(4._dp*(TrYw3adjYw3) +     & 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))))/1800._dp

 
DYx3 = oo16pi2*( betaYx31 + oo16pi2 * betaYx32 ) 

 
Else 
DYx3 = oo16pi2* betaYx31 
End If 
 
 
If (ThresholdCrossed.Lt.1) Then 
DYx3(1,:) = 0._dp 
DYb3(1,:) = 0._dp 
DYw3(1,:) = 0._dp 
End If 

If (ThresholdCrossed.Lt.2) Then 
DYx3(2,:) = 0._dp 
DYb3(2,:) = 0._dp 
DYw3(2,:) = 0._dp 
End If 

If (ThresholdCrossed.Lt.3) Then 
DYx3(3,:) = 0._dp 
DYb3(3,:) = 0._dp 
DYw3(3,:) = 0._dp 
End If 

Call ParametersToG111(Dg1,Dg2,Dg3,DYu,DYd,DYe,DYb3,DYw3,DYx3,f)

Iname = Iname - 1 
 
End Subroutine rge111


Subroutine rge555(len, T, GY, F) 
Implicit None 
Integer, Intent(in) :: len 
Real(dp), Intent(in) :: T, GY(len) 
Real(dp), Intent(out) :: F(len) 
Integer :: i1,i2
Real(dp) :: q 
Real(dp) :: g1,betag11,betag12,Dg1,g2,betag21,betag22,Dg2,g3,betag31,betag32,         & 
& Dg3,mHd2,betamHd21,betamHd22,DmHd2,mHu2,betamHu21,betamHu22,DmHu2
Complex(dp) :: Yu(3,3),betaYu1(3,3),betaYu2(3,3),DYu(3,3),adjYu(3,3),Yd(3,3)          & 
& ,betaYd1(3,3),betaYd2(3,3),DYd(3,3),adjYd(3,3),Ye(3,3),betaYe1(3,3),betaYe2(3,3)       & 
& ,DYe(3,3),adjYe(3,3),Yb3(3,3),betaYb31(3,3),betaYb32(3,3),DYb3(3,3),adjYb3(3,3)        & 
& ,Yw3(3,3),betaYw31(3,3),betaYw32(3,3),DYw3(3,3),adjYw3(3,3),Yx3(3,3),betaYx31(3,3)     & 
& ,betaYx32(3,3),DYx3(3,3),adjYx3(3,3),mue,betamue1,betamue2,Dmue,MXM3(3,3)              & 
& ,betaMXM31(3,3),betaMXM32(3,3),DMXM3(3,3),adjMXM3(3,3),MWM3(3,3),betaMWM31(3,3)        & 
& ,betaMWM32(3,3),DMWM3(3,3),adjMWM3(3,3),MGM3(3,3),betaMGM31(3,3),betaMGM32(3,3)        & 
& ,DMGM3(3,3),adjMGM3(3,3),MBM3(3,3),betaMBM31(3,3),betaMBM32(3,3),DMBM3(3,3)            & 
& ,adjMBM3(3,3),TYu(3,3),betaTYu1(3,3),betaTYu2(3,3),DTYu(3,3),adjTYu(3,3)               & 
& ,TYd(3,3),betaTYd1(3,3),betaTYd2(3,3),DTYd(3,3),adjTYd(3,3),TYe(3,3),betaTYe1(3,3)     & 
& ,betaTYe2(3,3),DTYe(3,3),adjTYe(3,3),TYb3(3,3),betaTYb31(3,3),betaTYb32(3,3)           & 
& ,DTYb3(3,3),adjTYb3(3,3),TYw3(3,3),betaTYw31(3,3),betaTYw32(3,3),DTYw3(3,3)            & 
& ,adjTYw3(3,3),TYx3(3,3),betaTYx31(3,3),betaTYx32(3,3),DTYx3(3,3),adjTYx3(3,3)          & 
& ,Bmue,betaBmue1,betaBmue2,DBmue,BMXM3(3,3),DBMXM3(3,3),BMWM3(3,3),DBMWM3(3,3)      & 
& ,BMGM3(3,3),DBMGM3(3,3),BMBM3(3,3),DBMBM3(3,3),mq2(3,3)         & 
& ,betamq21(3,3),betamq22(3,3),Dmq2(3,3),adjmq2(3,3),ml2(3,3),betaml21(3,3)              & 
& ,betaml22(3,3),Dml2(3,3),adjml2(3,3),md2(3,3),betamd21(3,3),betamd22(3,3)              & 
& ,Dmd2(3,3),adjmd2(3,3),mu2(3,3),betamu21(3,3),betamu22(3,3),Dmu2(3,3),adjmu2(3,3)      & 
& ,me2(3,3),betame21(3,3),betame22(3,3),Dme2(3,3),adjme2(3,3),mHw32(3,3),betamHw321(3,3) & 
& ,betamHw322(3,3),DmHw32(3,3),adjmHw32(3,3),mHg32(3,3),betamHg321(3,3),betamHg322(3,3)  & 
& ,DmHg32(3,3),adjmHg32(3,3),mHb32(3,3),betamHb321(3,3),betamHb322(3,3),DmHb32(3,3)      & 
& ,adjmHb32(3,3),mHx32(3,3),betamHx321(3,3),betamHx322(3,3),DmHx32(3,3),adjmHx32(3,3)    & 
& ,mHxb32(3,3),betamHxb321(3,3),betamHxb322(3,3),DmHxb32(3,3),adjmHxb32(3,3)             & 
& ,MassB,betaMassB1,betaMassB2,DMassB,MassWB,betaMassWB1,betaMassWB2,DMassWB,MassG    &
& ,betaMassG1,betaMassG2,DMassG
Complex(dp) :: Tr1(3),Tr2(3),Tr3(3) 
Real(dp) :: AbsMassB,AbsMassWB,AbsMassG
Complex(dp) :: md2adjYx3(3,3),md2CYd(3,3),me2CYe(3,3),mHb32CYb3(3,3),mHw32CYw3(3,3),mHxb32CYx3(3,3), & 
& ml2adjYb3(3,3),ml2adjYe(3,3),ml2adjYw3(3,3),mq2adjYd(3,3),mq2adjYu(3,3),               & 
& mu2CYu(3,3),Yb3adjYb3(3,3),YdadjYd(3,3),YeadjYe(3,3),YuadjYu(3,3),Yw3adjYw3(3,3),      & 
& Yx3adjYx3(3,3),adjYb3MBM3(3,3),adjYb3mHb32(3,3),adjYb3Yb3(3,3),       & 
& adjYb3TYb3(3,3),adjYdmd2(3,3),adjYdYd(3,3),adjYdTYd(3,3),adjYeme2(3,3),adjYeYe(3,3),   & 
& adjYeTYe(3,3),adjYumu2(3,3),adjYuYu(3,3),adjYuTYu(3,3),adjYw3mHw32(3,3),               & 
& adjYw3MWM3(3,3),adjYw3Yw3(3,3),adjYw3TYw3(3,3),adjYx3mHxb32(3,3),     & 
& adjYx3Yx3(3,3),adjYx3TYx3(3,3),CYb3ml2(3,3),CYb3TpYb3(3,3),            & 
& CYdmq2(3,3),CYdTpYd(3,3),CYdTpTYd(3,3),CYeml2(3,3),CYumq2(3,3),CYw3ml2(3,3),           & 
& CYw3TpYw3(3,3),CYx3md2(3,3),CYx3Yd(3,3),CYx3TYd(3,3),CYx3TpYx3(3,3),   & 
& CTYb3TpTYb3(3,3),CTYdTpTYd(3,3),CTYeTpTYe(3,3),CTYuTpTYu(3,3),         & 
& CTYw3TpTYw3(3,3),CTYx3TpTYx3(3,3),TYb3adjTYb3(3,3),TYdadjTYd(3,3),TYeadjTYe(3,3),      & 
& TYuadjTYu(3,3),TYw3adjTYw3(3,3),TYx3adjTYx3(3,3),TpYb3CYb3(3,3),TpYdCYd(3,3),          & 
& TpYeCYe(3,3),TpYuCYu(3,3),TpYw3CYw3(3,3),TpYx3CYx3(3,3),TpTYb3CTYb3(3,3),              & 
& TpTYdCTYd(3,3),TpTYeCTYe(3,3),TpTYuCTYu(3,3),TpTYw3CTYw3(3,3),TpTYx3CTYx3(3,3),        & 
& MBM3CYb3TpYb3(3,3),md2YdadjYd(3,3),md2adjYx3Yx3(3,3),              & 
& md2TpYx3CYx3(3,3),me2YeadjYe(3,3),mHb32Yb3adjYb3(3,3),mHw32Yw3adjYw3(3,3),             & 
& mHxb32Yx3adjYx3(3,3),ml2adjYb3Yb3(3,3),ml2adjYeYe(3,3),ml2adjYw3Yw3(3,3),              & 
& ml2TpYb3CYb3(3,3),ml2TpYeCYe(3,3),ml2TpYw3CYw3(3,3),mq2adjYdYd(3,3),mq2adjYuYu(3,3),   & 
& mq2TpYdCYd(3,3),mq2TpYuCYu(3,3),mu2YuadjYu(3,3),MWM3CYw3TpYw3(3,3),& 
& MXM3CYx3TpYx3(3,3),Yb3ml2adjYb3(3,3),Yb3adjYb3MBM3(3,3),           & 
& Yb3adjYb3mHb32(3,3),Yb3adjYb3Yb3(3,3),Yb3adjYb3TYb3(3,3),          & 
& Yb3adjYeYe(3,3),Yb3adjYeTYe(3,3),Yb3adjYw3Yw3(3,3),Yb3adjYw3TYw3(3,3),Ydmq2adjYd(3,3), & 
& YdadjYdmd2(3,3),YdadjYdYd(3,3),YdadjYdTYd(3,3),YdadjYuYu(3,3),YdadjYuTYu(3,3),         & 
& Yeml2adjYe(3,3),YeadjYb3Yb3(3,3),YeadjYb3TYb3(3,3),YeadjYeme2(3,3),YeadjYeYe(3,3),     & 
& YeadjYeTYe(3,3),YeadjYw3Yw3(3,3),YeadjYw3TYw3(3,3),Yumq2adjYu(3,3),YuadjYdYd(3,3),     & 
& YuadjYdTYd(3,3),YuadjYumu2(3,3),YuadjYuYu(3,3),YuadjYuTYu(3,3),Yw3ml2adjYw3(3,3),      & 
& Yw3adjYb3Yb3(3,3),Yw3adjYb3TYb3(3,3),Yw3adjYeYe(3,3),Yw3adjYeTYe(3,3),Yw3adjYw3mHw32(3,3),& 
& Yw3adjYw3MWM3(3,3),Yw3adjYw3Yw3(3,3),Yw3adjYw3TYw3(3,3),           & 
& Yx3md2adjYx3(3,3),Yx3adjYx3mHxb32(3,3),Yx3adjYx3Yx3(3,3),Yx3adjYx3TYx3(3,3),           & 
& Yx3CYdTpYd(3,3),Yx3CYdTpTYd(3,3),              & 
& TYb3adjYb3Yb3(3,3),TYb3adjYeYe(3,3),           & 
& TYb3adjYw3Yw3(3,3),TYdadjYdYd(3,3),TYdadjYuYu(3,3),TYeadjYb3Yb3(3,3),TYeadjYeYe(3,3),  & 
& TYeadjYw3Yw3(3,3),TYuadjYdYd(3,3),TYuadjYuYu(3,3),TYw3adjYb3Yb3(3,3),TYw3adjYeYe(3,3), & 
& TYw3adjYw3Yw3(3,3),TYx3adjYx3Yx3(3,3),TYx3CYdTpYd(3,3),            & 
& TpYb3mHb32CYb3(3,3),TpYb3CYb3ml2(3,3),TpYdmd2CYd(3,3),TpYdCYdmq2(3,3),TpYeme2CYe(3,3), & 
& TpYeCYeml2(3,3),TpYumu2CYu(3,3),TpYuCYumq2(3,3),TpYw3mHw32CYw3(3,3),TpYw3CYw3ml2(3,3)

Complex(dp) :: TpYx3mHxb32CYx3(3,3),TpYx3CYx3md2(3,3),TpYx3CYx3Yd(3,3),TpYx3CYx3TYd(3,3),             & 
& TpTYx3CYx3Yd(3,3)

Complex(dp) :: Yb3adjYe(3,3),Yb3adjYw3(3,3),Yb3adjTYb3(3,3),Yb3adjTYe(3,3),Yb3adjTYw3(3,3),          & 
& YdadjYu(3,3),YdadjTYd(3,3),YdadjTYu(3,3),YeadjYb3(3,3),YeadjYw3(3,3),YeadjTYb3(3,3),   & 
& YeadjTYe(3,3),YeadjTYw3(3,3),YuadjYd(3,3),YuadjTYd(3,3),YuadjTYu(3,3),Yw3adjYb3(3,3),  & 
& Yw3adjYe(3,3),Yw3adjTYb3(3,3),Yw3adjTYe(3,3),Yw3adjTYw3(3,3),Yx3adjTYx3(3,3),          & 
& Yx3CYd(3,3),Yx3CTYd(3,3),adjYdadjTYx3(3,3),adjYdTpYx3(3,3),           & 
& adjTYdadjYx3(3,3),CYb3TpYw3(3,3),CYeTpYb3(3,3),CYeTpYw3(3,3),          & 
& CYuTpYd(3,3),CYuTpTYd(3,3),CYw3TpYb3(3,3),               & 
& CTYb3adjYb3(3,3),CTYb3adjYe(3,3),CTYb3adjYw3(3,3),CTYb3TpYb3(3,3),     & 
& CTYdadjYd(3,3),CTYdadjYu(3,3),CTYdTpYd(3,3),CTYeadjYb3(3,3),CTYeadjYe(3,3),            & 
& CTYeadjYw3(3,3),CTYeTpYe(3,3),CTYuadjYd(3,3),CTYuadjYu(3,3),CTYuTpYu(3,3),             & 
& CTYw3adjYb3(3,3),CTYw3adjYe(3,3),CTYw3adjYw3(3,3),CTYw3TpYw3(3,3),CTYx3adjYx3(3,3),    & 
& CTYx3TYd(3,3),CTYx3TpYd(3,3),CTYx3TpYx3(3,3),CTYx3TpTYd(3,3),TYb3adjYb3(3,3),          & 
& TYb3adjYe(3,3),TYb3adjYw3(3,3),TYb3adjTYe(3,3),TYb3adjTYw3(3,3),TYdadjYd(3,3),         & 
& TYdadjYu(3,3),TYdadjTYu(3,3),TYeadjYb3(3,3),TYeadjYe(3,3),TYeadjYw3(3,3),              & 
& TYeadjTYb3(3,3),TYeadjTYw3(3,3),TYuadjYd(3,3),TYuadjYu(3,3),TYuadjTYd(3,3),            & 
& TYw3adjYb3(3,3),TYw3adjYe(3,3),TYw3adjYw3(3,3),TYw3adjTYb3(3,3),TYw3adjTYe(3,3),       & 
& TYx3adjYx3(3,3),TYx3CTYd(3,3),TpYb3CTYb3(3,3),TpYdadjYx3(3,3),TpYdadjTYx3(3,3),        & 
& TpYdCTYd(3,3),TpYeCTYe(3,3),TpYuCTYu(3,3),TpYw3CTYw3(3,3),TpYx3CTYx3(3,3),             & 
& TpTYb3CYb3(3,3),TpTYdadjYx3(3,3),TpTYdCYd(3,3),TpTYeCYe(3,3),TpTYuCYu(3,3),            & 
& TpTYw3CYw3(3,3),TpTYx3CYx3(3,3),md2YdadjYu(3,3),me2YeadjYb3(3,3),me2YeadjYw3(3,3),     & 
& mHb32Yb3adjYe(3,3),mHb32Yb3adjYw3(3,3),mHw32Yw3adjYb3(3,3),mHw32Yw3adjYe(3,3),         & 
& mHxb32Yx3CYd(3,3),mq2TpYdadjYx3(3,3),mu2YuadjYd(3,3),Yb3ml2adjYe(3,3),Yb3ml2adjYw3(3,3),& 
& Yb3adjYeme2(3,3),Yb3adjYw3mHw32(3,3),Yb3adjYw3MWM3(3,3),           & 
& Ydmq2adjYu(3,3),YdadjYdTpYx3(3,3),YdadjYumu2(3,3),YdadjTYdadjYx3(3,3),& 
& Yeml2adjYb3(3,3),Yeml2adjYw3(3,3),YeadjYb3MBM3(3,3),YeadjYb3mHb32(3,3),& 
& YeadjYw3mHw32(3,3),YeadjYw3MWM3(3,3),Yumq2adjYd(3,3),               & 
& YuadjYdmd2(3,3),Yw3ml2adjYb3(3,3),Yw3ml2adjYe(3,3),Yw3adjYb3MBM3(3,3),Yw3adjYb3mHb32(3,3),& 
& Yw3adjYeme2(3,3),Yx3md2CYd(3,3),Yx3CYdmq2(3,3),adjYb3Yb3adjYb3(3,3),& 
& adjYb3Yb3adjYe(3,3),adjYb3Yb3adjYw3(3,3),adjYb3Yb3adjTYb3(3,3),adjYb3Yb3adjTYe(3,3),   & 
& adjYb3Yb3adjTYw3(3,3),adjYb3TYb3adjYb3(3,3),adjYb3TYb3adjYe(3,3),adjYb3TYb3adjYw3(3,3),& 
& adjYb3TYb3adjTYb3(3,3),adjYb3TYb3adjTYe(3,3),adjYb3TYb3adjTYw3(3,3),adjYdYdadjYd(3,3), & 
& adjYdYdadjYu(3,3),adjYdYdadjTYd(3,3),adjYdYdadjTYu(3,3),adjYdadjYx3Yx3(3,3),           & 
& adjYdTYdadjYd(3,3),adjYdTYdadjYu(3,3),adjYdTYdadjTYd(3,3),adjYdTYdadjTYu(3,3),         & 
& adjYdTpYx3CYx3(3,3),adjYdTpYx3CTYx3(3,3),adjYdTpTYx3CTYx3(3,3),adjYeYeadjYb3(3,3),     & 
& adjYeYeadjYe(3,3),adjYeYeadjYw3(3,3),adjYeYeadjTYb3(3,3),adjYeYeadjTYe(3,3),           & 
& adjYeYeadjTYw3(3,3),adjYeTYeadjYb3(3,3),adjYeTYeadjYe(3,3),adjYeTYeadjYw3(3,3),        & 
& adjYeTYeadjTYb3(3,3),adjYeTYeadjTYe(3,3),adjYeTYeadjTYw3(3,3),adjYuYuadjYd(3,3)

Complex(dp) :: adjYuYuadjYu(3,3),adjYuYuadjTYd(3,3),adjYuYuadjTYu(3,3),adjYuTYuadjYd(3,3),            & 
& adjYuTYuadjYu(3,3),adjYuTYuadjTYd(3,3),adjYuTYuadjTYu(3,3),adjYw3Yw3adjYb3(3,3),       & 
& adjYw3Yw3adjYe(3,3),adjYw3Yw3adjYw3(3,3),adjYw3Yw3adjTYb3(3,3),adjYw3Yw3adjTYe(3,3),   & 
& adjYw3Yw3adjTYw3(3,3),adjYw3TYw3adjYb3(3,3),adjYw3TYw3adjYe(3,3),adjYw3TYw3adjYw3(3,3),& 
& adjYw3TYw3adjTYb3(3,3),adjYw3TYw3adjTYe(3,3),adjYw3TYw3adjTYw3(3,3),adjYx3Yx3adjYx3(3,3),& 
& adjYx3Yx3adjTYx3(3,3),adjYx3Yx3CYd(3,3),adjYx3Yx3CTYd(3,3),adjYx3TYx3adjYx3(3,3),      & 
& adjYx3TYx3adjTYx3(3,3),adjYx3TYx3CTYd(3,3),adjTYb3TYb3adjYb3(3,3),adjTYb3TYb3adjYe(3,3),& 
& adjTYb3TYb3adjYw3(3,3),adjTYdTYdadjYd(3,3),adjTYdTYdadjYu(3,3),adjTYeTYeadjYb3(3,3),   & 
& adjTYeTYeadjYe(3,3),adjTYeTYeadjYw3(3,3),adjTYuTYuadjYd(3,3),adjTYuTYuadjYu(3,3),      & 
& adjTYw3TYw3adjYb3(3,3),adjTYw3TYw3adjYe(3,3),adjTYw3TYw3adjYw3(3,3),adjTYx3TYx3adjYx3(3,3),& 
& CYb3TpYb3CYb3(3,3),CYb3TpYb3CTYb3(3,3),CYb3TpYw3CYw3(3,3),CYb3TpYw3CTYw3(3,3),         & 
& CYb3TpTYb3CTYb3(3,3),CYb3TpTYw3CTYw3(3,3),CYdTpYdadjYx3(3,3),CYdTpYdadjTYx3(3,3),      & 
& CYdTpYdCYd(3,3),CYdTpYdCTYd(3,3),CYdTpTYdCTYd(3,3),CYeTpYeCYe(3,3),CYeTpYeCTYe(3,3),   & 
& CYeTpTYeCTYe(3,3),CYuTpYuCYu(3,3),CYuTpYuCTYu(3,3),CYuTpTYuCTYu(3,3),CYw3TpYb3CYb3(3,3),& 
& CYw3TpYb3CTYb3(3,3),CYw3TpYw3CYw3(3,3),CYw3TpYw3CTYw3(3,3),CYw3TpTYb3CTYb3(3,3),       & 
& CYw3TpTYw3CTYw3(3,3),CYx3YdadjYd(3,3),CYx3TpYx3CYx3(3,3),CYx3TpYx3CTYx3(3,3),          & 
& CYx3TpTYx3CTYx3(3,3),CTYb3TpTYb3CYb3(3,3),CTYb3TpTYw3CYw3(3,3),CTYdTpTYdadjYx3(3,3),   & 
& CTYdTpTYdCYd(3,3),CTYeTpTYeCYe(3,3),CTYuTpTYuCYu(3,3),CTYw3TpTYb3CYb3(3,3),            & 
& CTYw3TpTYw3CYw3(3,3),CTYx3TpTYx3CYx3(3,3),TYb3TpYb3CYb3(3,3),      & 
& TYb3TpYw3CYw3(3,3),TYdadjYdadjTYx3(3,3),TYdTpYdCYd(3,3),            & 
& TYeTpYeCYe(3,3),TYuTpYuCYu(3,3),& 
& TYw3TpYb3CYb3(3,3),TYw3TpYw3CYw3(3,3),TYx3CTYdTpYd(3,3),TYx3TpYx3CYx3(3,3),            & 
& TpYb3CYb3TpYb3(3,3),TpYb3CYb3TpYw3(3,3),     & 
& TpYb3CTYb3adjYb3(3,3),TpYb3CTYb3adjYe(3,3),TpYb3CTYb3adjYw3(3,3),TpYdmd2adjYx3(3,3),   & 
& TpYdadjYx3mHxb32(3,3),TpYdadjYx3Yx3(3,3),TpYdadjYx3TYx3(3,3),TpYdCYdTpYd(3,3),         & 
& TpYdCYdTpTYd(3,3),TpYdCTYdadjYd(3,3),TpYdCTYdadjYu(3,3),TpYeCYeTpYb3(3,3),             & 
& TpYeCYeTpYw3(3,3),TpYeCTYeadjYb3(3,3),           & 
& TpYeCTYeadjYe(3,3),TpYeCTYeadjYw3(3,3),TpYuCYuTpYd(3,3),TpYuCYuTpTYd(3,3),             & 
& TpYuCTYuadjYd(3,3),TpYuCTYuadjYu(3,3),TpYw3CYw3TpYb3(3,3),TpYw3CYw3TpYw3(3,3),         & 
& TpYw3CTYw3adjYb3(3,3),TpYw3CTYw3adjYe(3,3),  & 
& TpYw3CTYw3adjYw3(3,3),TpYx3CYx3TpYx3(3,3),TpYx3CTYx3adjYx3(3,3),  & 
& TpYx3CTYx3TYd(3,3),TpYx3CTYx3TpTYd(3,3),     & 
& TpTYb3CTYb3adjYb3(3,3),TpTYb3CTYb3adjYe(3,3),TpTYb3CTYb3adjYw3(3,3),TpTYdCYdTpYd(3,3), & 
& TpTYdCTYdadjYd(3,3),TpTYdCTYdadjYu(3,3),         & 
& TpTYeCTYeadjYb3(3,3),TpTYeCTYeadjYe(3,3),TpTYeCTYeadjYw3(3,3),TpTYuCYuTpYd(3,3),       & 
& TpTYuCTYuadjYd(3,3),TpTYuCTYuadjYu(3,3),     & 
& TpTYw3CTYw3adjYb3(3,3),TpTYw3CTYw3adjYe(3,3),TpTYw3CTYw3adjYw3(3,3)

Complex(dp) :: TpTYx3CTYx3adjYx3(3,3),TpTYx3CTYx3TpYd(3,3),md2adjYx3Yx3adjYx3(3,3),md2adjYx3Yx3CYd(3,3),& 
& md2CYdTpYdadjYx3(3,3),md2CYdTpYdCYd(3,3),me2CYeTpYeCYe(3,3),mHb32CYb3TpYb3CYb3(3,3),   & 
& mHb32CYb3TpYw3CYw3(3,3),mHw32CYw3TpYb3CYb3(3,3),mHw32CYw3TpYw3CYw3(3,3),               & 
& mHxb32CYx3TpYx3CYx3(3,3),ml2adjYb3Yb3adjYb3(3,3),ml2adjYb3Yb3adjYe(3,3),               & 
& ml2adjYb3Yb3adjYw3(3,3),ml2adjYeYeadjYb3(3,3),ml2adjYeYeadjYe(3,3),ml2adjYeYeadjYw3(3,3),& 
& ml2adjYw3Yw3adjYb3(3,3),ml2adjYw3Yw3adjYe(3,3),ml2adjYw3Yw3adjYw3(3,3),mq2adjYdYdadjYd(3,3),& 
& mq2adjYdYdadjYu(3,3),mq2adjYuYuadjYd(3,3),mq2adjYuYuadjYu(3,3),mq2adjYx3Yx3CYd(3,3),   & 
& mu2CYuTpYuCYu(3,3),Yb3adjYb3Yb3adjYb3(3,3),Yb3adjYb3TYb3adjYb3(3,3),Yb3adjYb3TYb3adjTYb3(3,3),& 
& Yb3adjYeYeadjYb3(3,3),Yb3adjYeTYeadjYb3(3,3),Yb3adjYeTYeadjTYb3(3,3),Yb3adjYw3Yw3adjYb3(3,3),& 
& Yb3adjYw3TYw3adjYb3(3,3),Yb3adjYw3TYw3adjTYb3(3,3),Yb3adjTYb3TYb3adjYb3(3,3),          & 
& Yb3adjTYeTYeadjYb3(3,3),Yb3adjTYw3TYw3adjYb3(3,3),Yb3TpTYb3CTYb3adjYb3(3,3),           & 
& Yb3TpTYeCTYeadjYb3(3,3),Yb3TpTYw3CTYw3adjYb3(3,3),YdadjYdYdadjYd(3,3),YdadjYdTYdadjYd(3,3),& 
& YdadjYdTYdadjTYd(3,3),YdadjYdTpYx3CYx3(3,3),YdadjYdTpTYx3CTYx3(3,3),YdadjYuYuadjYd(3,3),& 
& YdadjYuTYuadjYd(3,3),YdadjYuTYuadjTYd(3,3),YdadjTYdTYdadjYd(3,3),YdadjTYuTYuadjYd(3,3),& 
& YdTpTYdCTYdadjYd(3,3),YdTpTYuCTYuadjYd(3,3),YeadjYb3Yb3adjYe(3,3),YeadjYb3TYb3adjYe(3,3),& 
& YeadjYb3TYb3adjTYe(3,3),YeadjYeYeadjYe(3,3),YeadjYeTYeadjYe(3,3),YeadjYeTYeadjTYe(3,3),& 
& YeadjYw3Yw3adjYe(3,3),YeadjYw3TYw3adjYe(3,3),YeadjYw3TYw3adjTYe(3,3),YeadjTYb3TYb3adjYe(3,3),& 
& YeadjTYeTYeadjYe(3,3),YeadjTYw3TYw3adjYe(3,3),YeTpTYb3CTYb3adjYe(3,3),YeTpTYeCTYeadjYe(3,3),& 
& YeTpTYw3CTYw3adjYe(3,3),YuadjYdYdadjYu(3,3),YuadjYdTYdadjYu(3,3),YuadjYdTYdadjTYu(3,3),& 
& YuadjYuYuadjYu(3,3),YuadjYuTYuadjYu(3,3),YuadjYuTYuadjTYu(3,3),YuadjTYdTYdadjYu(3,3),  & 
& YuadjTYuTYuadjYu(3,3),YuTpTYdCTYdadjYu(3,3),YuTpTYuCTYuadjYu(3,3),Yw3adjYb3Yb3adjYw3(3,3),& 
& Yw3adjYb3TYb3adjYw3(3,3),Yw3adjYb3TYb3adjTYw3(3,3),Yw3adjYeYeadjYw3(3,3),              & 
& Yw3adjYeTYeadjYw3(3,3),Yw3adjYeTYeadjTYw3(3,3),Yw3adjYw3Yw3adjYw3(3,3),Yw3adjYw3TYw3adjYw3(3,3),& 
& Yw3adjYw3TYw3adjTYw3(3,3),Yw3adjTYb3TYb3adjYw3(3,3),Yw3adjTYeTYeadjYw3(3,3),           & 
& Yw3adjTYw3TYw3adjYw3(3,3),Yw3TpTYb3CTYb3adjYw3(3,3),Yw3TpTYeCTYeadjYw3(3,3),           & 
& Yw3TpTYw3CTYw3adjYw3(3,3),Yx3adjYx3Yx3adjYx3(3,3),Yx3adjYx3TYx3adjYx3(3,3),            & 
& Yx3adjYx3TYx3adjTYx3(3,3),Yx3adjTYx3TYx3adjYx3(3,3),Yx3CYdTpYdadjYx3(3,3),             & 
& Yx3CTYdTpTYdadjYx3(3,3),Yx3TYdadjYdadjTYx3(3,3),Yx3TpTYx3CTYx3adjYx3(3,3),             & 
& adjYb3mHb32Yb3adjYb3(3,3),adjYb3mHb32Yb3adjYe(3,3),adjYb3mHb32Yb3adjYw3(3,3),          & 
& adjYb3Yb3ml2adjYb3(3,3),adjYb3Yb3ml2adjYe(3,3),adjYb3Yb3ml2adjYw3(3,3),adjYb3Yb3adjYb3MBM3(3,3),& 
& adjYb3Yb3adjYb3mHb32(3,3),adjYb3Yb3adjYb3Yb3(3,3),           & 
& adjYb3Yb3adjYb3TYb3(3,3),adjYb3Yb3adjYeme2(3,3),adjYb3Yb3adjYeYe(3,3),adjYb3Yb3adjYeTYe(3,3),& 
& adjYb3Yb3adjYw3mHw32(3,3),adjYb3Yb3adjYw3MWM3(3,3),adjYb3Yb3adjYw3Yw3(3,3),            & 
& adjYb3Yb3adjYw3TYw3(3,3),          & 
& adjYb3TYb3adjYb3Yb3(3,3),adjYb3TYb3adjYeYe(3,3),             & 
& adjYb3TYb3adjYw3Yw3(3,3),adjYdmd2YdadjYd(3,3),adjYdmd2YdadjYu(3,3),adjYdYdmq2adjYd(3,3),& 
& adjYdYdmq2adjYu(3,3),adjYdYdadjYdmd2(3,3),adjYdYdadjYdYd(3,3),adjYdYdadjYdTYd(3,3)

Complex(dp) :: adjYdYdadjYumu2(3,3),adjYdYdadjYuYu(3,3),adjYdYdadjYuTYu(3,3),adjYdTYdadjYdYd(3,3),    & 
& adjYdTYdadjYuYu(3,3),adjYdTpYx3CYx3Yd(3,3),adjYdTpYx3CYx3TYd(3,3),adjYdTpYx3CTYx3TYd(3,3),& 
& adjYdTpYx3CTYx3TpTYd(3,3),adjYdTpTYx3CYx3Yd(3,3),adjYdTpTYx3CTYx3TpYd(3,3),            & 
& adjYeme2YeadjYb3(3,3),adjYeme2YeadjYe(3,3),adjYeme2YeadjYw3(3,3),adjYeYeml2adjYb3(3,3),& 
& adjYeYeml2adjYe(3,3),adjYeYeml2adjYw3(3,3),adjYeYeadjYb3MBM3(3,3),adjYeYeadjYb3mHb32(3,3),& 
& adjYeYeadjYb3Yb3(3,3),adjYeYeadjYb3TYb3(3,3),adjYeYeadjYeme2(3,3),& 
& adjYeYeadjYeYe(3,3),adjYeYeadjYeTYe(3,3),adjYeYeadjYw3mHw32(3,3),adjYeYeadjYw3MWM3(3,3),& 
& adjYeYeadjYw3Yw3(3,3),adjYeYeadjYw3TYw3(3,3), & 
& adjYeTYeadjYb3Yb3(3,3),adjYeTYeadjYeYe(3,3),adjYeTYeadjYw3Yw3(3,3),& 
& adjYumu2YuadjYd(3,3),adjYumu2YuadjYu(3,3),adjYuYumq2adjYd(3,3),adjYuYumq2adjYu(3,3),   & 
& adjYuYuadjYdmd2(3,3),adjYuYuadjYdYd(3,3),adjYuYuadjYdTYd(3,3),adjYuYuadjYumu2(3,3),    & 
& adjYuYuadjYuYu(3,3),adjYuYuadjYuTYu(3,3),adjYuTYuadjYdYd(3,3),adjYuTYuadjYuYu(3,3),    & 
& adjYw3mHw32Yw3adjYb3(3,3),adjYw3mHw32Yw3adjYe(3,3),adjYw3mHw32Yw3adjYw3(3,3),          & 
& adjYw3Yw3ml2adjYb3(3,3),adjYw3Yw3ml2adjYe(3,3),adjYw3Yw3ml2adjYw3(3,3),adjYw3Yw3adjYb3MBM3(3,3),& 
& adjYw3Yw3adjYb3mHb32(3,3),adjYw3Yw3adjYb3Yb3(3,3),           & 
& adjYw3Yw3adjYb3TYb3(3,3),adjYw3Yw3adjYeme2(3,3),adjYw3Yw3adjYeYe(3,3),adjYw3Yw3adjYeTYe(3,3),& 
& adjYw3Yw3adjYw3mHw32(3,3),adjYw3Yw3adjYw3MWM3(3,3),adjYw3Yw3adjYw3Yw3(3,3),            & 
& adjYw3Yw3adjYw3TYw3(3,3),          & 
& adjYw3TYw3adjYb3Yb3(3,3),adjYw3TYw3adjYeYe(3,3),             & 
& adjYw3TYw3adjYw3Yw3(3,3),adjYx3mHxb32Yx3adjYx3(3,3),adjYx3mHxb32Yx3CYd(3,3),           & 
& adjYx3Yx3md2adjYx3(3,3),adjYx3Yx3md2CYd(3,3),adjYx3Yx3adjYx3mHxb32(3,3),               & 
& adjYx3Yx3adjYx3Yx3(3,3),adjYx3Yx3adjYx3TYx3(3,3),adjYx3Yx3CYdmq2(3,3),adjYx3TYx3adjYx3Yx3(3,3),& 
& adjYx3TYx3CYdTpYd(3,3),adjYx3TYx3CTYdTpYd(3,3),adjTYb3TYb3TpYb3CYb3(3,3),              & 
& adjTYb3TYb3TpYw3CYw3(3,3),adjTYdTYdTpYdCYd(3,3),adjTYeTYeTpYeCYe(3,3),adjTYuTYuTpYuCYu(3,3),& 
& adjTYw3TYw3TpYb3CYb3(3,3),adjTYw3TYw3TpYw3CYw3(3,3),adjTYx3TYx3TpYx3CYx3(3,3),         & 
& CYb3ml2TpYb3CYb3(3,3),CYb3ml2TpYw3CYw3(3,3),CYb3TpYb3mHb32CYb3(3,3),CYb3TpYb3CYb3ml2(3,3),& 
& CYb3TpYb3CYb3TpYb3(3,3),CYb3TpYeCYeTpYb3(3,3),& 
& CYb3TpYw3mHw32CYw3(3,3),CYb3TpYw3CYw3ml2(3,3),CYb3TpYw3CYw3TpYb3(3,3),& 
& CYdmq2TpYdadjYx3(3,3),CYdmq2TpYdCYd(3,3),CYdTpYdmd2adjYx3(3,3),CYdTpYdmd2CYd(3,3),     & 
& CYdTpYdadjYx3mHxb32(3,3),CYdTpYdadjYx3Yx3(3,3),CYdTpYdadjYx3TYx3(3,3),CYdTpYdCYdmq2(3,3),& 
& CYdTpYdCYdTpYd(3,3),CYdTpYdCYdTpTYd(3,3),CYdTpYuCYuTpYd(3,3),CYdTpYuCYuTpTYd(3,3),     & 
& CYdTpTYdCYdTpYd(3,3),CYdTpTYuCYuTpYd(3,3),CYeml2TpYeCYe(3,3),CYeTpYeme2CYe(3,3),       & 
& CYeTpYeCYeml2(3,3),CYumq2TpYuCYu(3,3),CYuTpYumu2CYu(3,3),CYuTpYuCYumq2(3,3),           & 
& CYw3ml2TpYb3CYb3(3,3),CYw3ml2TpYw3CYw3(3,3),CYw3TpYb3mHb32CYb3(3,3),CYw3TpYb3CYb3ml2(3,3),& 
& CYw3TpYb3CYb3TpYw3(3,3),CYw3TpYeCYeTpYw3(3,3),& 
& CYw3TpYw3mHw32CYw3(3,3),CYw3TpYw3CYw3ml2(3,3),CYw3TpYw3CYw3TpYw3(3,3)

Complex(dp) ::               & 
& CYx3md2TpYx3CYx3(3,3),CYx3YdadjYdTpYx3(3,3),& 
& CYx3TpYx3mHxb32CYx3(3,3),CYx3TpYx3CYx3md2(3,3),CYx3TpYx3CYx3Yd(3,3),CYx3TpYx3CYx3TYd(3,3),& 
& CYx3TpYx3CYx3TpYx3(3,3),CYx3TpTYx3CYx3Yd(3,3),& 
& TYb3adjYb3Yb3adjTYb3(3,3),TYb3adjYeYeadjTYb3(3,3),TYb3adjYw3Yw3adjTYb3(3,3),           & 
& TYb3TpYb3CTYb3adjYb3(3,3),TYb3TpYeCTYeadjYb3(3,3),TYb3TpYw3CTYw3adjYb3(3,3),           & 
& TYdadjYdYdadjTYd(3,3),TYdadjYdadjYx3Yx3(3,3),TYdadjYuYuadjTYd(3,3),TYdTpYdCTYdadjYd(3,3),& 
& TYdTpYuCTYuadjYd(3,3),TYeadjYb3Yb3adjTYe(3,3),TYeadjYeYeadjTYe(3,3),TYeadjYw3Yw3adjTYe(3,3),& 
& TYeTpYb3CTYb3adjYe(3,3),TYeTpYeCTYeadjYe(3,3),TYeTpYw3CTYw3adjYe(3,3),TYuadjYdYdadjTYu(3,3),& 
& TYuadjYuYuadjTYu(3,3),TYuTpYdCTYdadjYu(3,3),TYuTpYuCTYuadjYu(3,3),TYw3adjYb3Yb3adjTYw3(3,3),& 
& TYw3adjYeYeadjTYw3(3,3),TYw3adjYw3Yw3adjTYw3(3,3),TYw3TpYb3CTYb3adjYw3(3,3),           & 
& TYw3TpYeCTYeadjYw3(3,3),TYw3TpYw3CTYw3adjYw3(3,3),TYx3YdadjTYdadjYx3(3,3),             & 
& TYx3adjYx3Yx3adjTYx3(3,3),TYx3CYdTpYdadjTYx3(3,3),TYx3TpYx3CTYx3adjYx3(3,3),           & 
& TpYb3CYb3TpYb3CYb3(3,3),TpYb3CYb3TpYw3CYw3(3,3),TpYb3CYb3TpTYb3CTYb3(3,3),             & 
& TpYb3CYb3TpTYw3CTYw3(3,3),TpYb3CTYb3TpTYb3CYb3(3,3),TpYb3CTYb3TpTYw3CYw3(3,3),         & 
& TpYdadjYdTpTYx3CTYx3(3,3),TpYdadjYx3Yx3CYd(3,3),TpYdadjYx3TYx3CTYd(3,3),               & 
& TpYdCYdTpYdCYd(3,3),TpYdCYdTpTYdCTYd(3,3),TpYdCTYdTpTYdCYd(3,3),TpYeCYeTpYeCYe(3,3),   & 
& TpYeCYeTpTYeCTYe(3,3),TpYeCTYeTpTYeCYe(3,3),TpYuCYuTpYuCYu(3,3),TpYuCYuTpTYuCTYu(3,3), & 
& TpYuCTYuTpTYuCYu(3,3),TpYw3CYw3TpYb3CYb3(3,3),TpYw3CYw3TpYw3CYw3(3,3),TpYw3CYw3TpTYb3CTYb3(3,3),& 
& TpYw3CYw3TpTYw3CTYw3(3,3),TpYw3CTYw3TpTYb3CYb3(3,3),TpYw3CTYw3TpTYw3CYw3(3,3),         & 
& TpYx3CYx3YdadjYd(3,3),TpYx3CYx3TpYx3CYx3(3,3),TpYx3CYx3TpTYx3CTYx3(3,3),               & 
& TpYx3CTYx3TpTYx3CYx3(3,3),TpTYb3CYb3TpYb3CTYb3(3,3),TpTYb3CYb3TpYw3CTYw3(3,3),         & 
& TpTYdadjYdTpYx3CTYx3(3,3),TpTYdadjYx3Yx3CTYd(3,3),TpTYdCYdTpYdCTYd(3,3),               & 
& TpTYeCYeTpYeCTYe(3,3),TpTYuCYuTpYuCTYu(3,3),TpTYw3CYw3TpYb3CTYb3(3,3),TpTYw3CYw3TpYw3CTYw3(3,3),& 
& TpTYx3CYx3TpYx3CTYx3(3,3),MBM3CYb3TpYb3CYb3TpYb3(3,3),    & 
& MBM3CYb3TpYeCYeTpYb3(3,3),MBM3CYb3TpYw3CYw3TpYb3(3,3),      & 
& md2YdadjYdYdadjYd(3,3),md2YdadjYdTpYx3CYx3(3,3),          & 
& md2YdadjYuYuadjYd(3,3),md2adjYx3Yx3adjYx3Yx3(3,3),md2TpYx3CYx3YdadjYd(3,3),            & 
& md2TpYx3CYx3TpYx3CYx3(3,3),me2YeadjYb3Yb3adjYe(3,3),me2YeadjYeYeadjYe(3,3),            & 
& me2YeadjYw3Yw3adjYe(3,3),mHb32Yb3adjYb3Yb3adjYb3(3,3),mHb32Yb3adjYeYeadjYb3(3,3),      & 
& mHb32Yb3adjYw3Yw3adjYb3(3,3),mHw32Yw3adjYb3Yb3adjYw3(3,3),mHw32Yw3adjYeYeadjYw3(3,3),  & 
& mHw32Yw3adjYw3Yw3adjYw3(3,3),mHxb32Yx3adjYx3Yx3adjYx3(3,3),mHxb32Yx3CYdTpYdadjYx3(3,3),& 
& mHxb32CYx3YdadjYdTpYx3(3,3),ml2adjYb3Yb3adjYb3Yb3(3,3),ml2adjYb3Yb3adjYeYe(3,3),       & 
& ml2adjYb3Yb3adjYw3Yw3(3,3),ml2adjYeYeadjYb3Yb3(3,3),ml2adjYeYeadjYeYe(3,3),            & 
& ml2adjYeYeadjYw3Yw3(3,3),ml2adjYw3Yw3adjYb3Yb3(3,3),ml2adjYw3Yw3adjYeYe(3,3),          & 
& ml2adjYw3Yw3adjYw3Yw3(3,3),ml2TpYb3CYb3TpYb3CYb3(3,3),ml2TpYb3CYb3TpYw3CYw3(3,3)

! the corresponding RGEs are commented out
! complex(dp) :: MBM3CYb3TpYw3CYw3TpTYb3(3,3),MBM3CYb3TpTYb3(3,3),CYb3TpYw3CYw3TpTYb3(3,3), &
!  & MBM3CYb3TpTYb3CYb3TpYb3(3,3),MBM3CYb3TpYb3CYb3TpTYb3(3,3),MBM3CYb3TpYeCYeTpTYb3(3,3), &
!  & MBM3CYb3TpTYeCYeTpYb3(3,3),MWM3CYw3TpTYw3CYw3TpYw3(3,3),MWM3CYw3TpTYw3(3,3),  & 
!  & MWM3CYw3TpTYeCYeTpYw3(3,3),MWM3CYw3TpTYb3CYb3TpYw3(3,3),MWM3CYw3TpYb3CYb3TpTYw3(3,3), &
!  & MWM3CYw3TpYeCYeTpTYw3(3,3),MXM3CYx3TpTYx3CYx3TpYx3(3,3),MXM3CYx3TpTYx3(3,3),
!  & MWM3CYw3TpYw3CYw3TpTYw3(3,3),MXM3CYx3TpYx3CYx3TpTYx3(3,3),TYb3adjYeYeadjYb3MBM3(3,3),  & 
!  & MXM3CYx3YdadjYdTpTYx3(3,3), MXM3CYx3TYdadjYdTpYx3(3,3),TYb3adjYb3Yb3adjYb3MBM3(3,3), &
!  & TYb3adjYw3Yw3adjYb3MBM3(3,3),TYw3adjYb3Yb3adjYw3MWM3(3,3),TYw3adjYeYeadjYw3MWM3(3,3), &
!  & TYw3adjYw3Yw3adjYw3MWM3(3,3),Yb3adjYb3TYb3adjYb3MBM3(3,3),Yb3adjYeTYeadjYb3MBM3(3,3), &
!  & Yb3adjYw3TYw3adjYb3MBM3(3,3),Yw3adjYb3TYb3adjYw3MWM3(3,3),Yw3adjYeTYeadjYw3MWM3(3,3), &
!  & Yw3adjYw3TYw3adjYw3MWM3(3,3),adjYb3TYb3adjYb3MBM3(3,3),adjYb3TYb3adjYw3MWM3(3,3),     &
!  & adjYeTYeadjYb3MBM3(3,3),adjYeTYeadjYw3MWM3(3,3),adjYw3TYw3adjYb3MBM3(3,3),            &
!  & adjYw3TYw3adjYw3MWM3(3,3),CYb3TpTYb3CYb3TpYb3(3,3),CYb3TpTYeCYeTpYb3(3,3),
!  & CYb3TpTYw3CYw3TpYb3(3,3),CYb3TpYb3CYb3TpTYb3(3,3),CYb3TpYeCYeTpTYb3(3,3),              &
!  & CYw3TpTYb3CYb3TpYw3(3,3),CYw3TpTYeCYeTpYw3(3,3),CYw3TpTYw3CYw3TpYw3(3,3),    &
!  & CYw3TpYb3CYb3TpTYw3(3,3),CYw3TpYeCYeTpTYw3(3,3),CYw3TpYw3CYw3TpTYw3(3,3), &
!  & CYx3TpTYx3CYx3TpYx3(3,3),CYx3TpYx3CYx3TpTYx3(3,3),CYx3TYdadjYdTpYx3(3,3), &
!  & CYx3YdadjYdTpTYx3(3,3),TpTYb3CYb3TpYb3(3,3),TpTYb3CYb3TpYw3(3,3),TpTYeCYeTpYb3(3,3), &
!  & TpTYeCYeTpYw3(3,3),TpTYw3CYw3TpYb3(3,3),TpTYw3CYw3TpYw3(3,3),TpTYx3CYx3TpYx3(3,3),   &
!  & TpYb3CYb3TpTYb3(3,3),TpYb3CYb3TpTYw3(3,3),TpYeCYeTpTYb3(3,3),TpYeCYeTpTYw3(3,3), &
!  & TpYw3CYw3TpTYb3(3,3),TpYw3CYw3TpTYw3(3,3),TpYx3CYx3TpTYx3(3,3),TYb3adjYb3MBM3(3,3), &
!  & TYb3adjYw3MWM3(3,3),TYdadjYdTpYx3(3,3),TYeadjYb3MBM3(3,3),TYeadjYw3MWM3(3,3), &
!  & TYw3adjYb3MBM3(3,3),TYw3adjYw3MWM3(3,3),YdadjYdTpTYx3(3,3),adjYdTpTYx3(3,3), &
!  & CYb3TpTYb3(3,3),CYb3TpTYw3(3,3),CYeTpTYb3(3,3),CYeTpTYw3(3,3),CYw3TpTYb3(3,3), &
!  & CYw3TpTYw3(3,3),CYx3TpTYx3(3,3),

Complex(dp) :: ml2TpYeCYeTpYeCYe(3,3),ml2TpYw3CYw3TpYb3CYb3(3,3),ml2TpYw3CYw3TpYw3CYw3(3,3),          & 
& mq2adjYdYdadjYdYd(3,3),mq2adjYdYdadjYuYu(3,3),mq2adjYdTpYx3CYx3Yd(3,3),mq2adjYuYuadjYdYd(3,3),& 
& mq2adjYuYuadjYuYu(3,3),mq2TpYdCYdTpYdCYd(3,3),mq2TpYuCYuTpYuCYu(3,3),mu2YuadjYdYdadjYu(3,3),& 
& mu2YuadjYuYuadjYu(3,3),MWM3CYw3TpYb3CYb3TpYw3(3,3),       & 
& MWM3CYw3TpYeCYeTpYw3(3,3),MWM3CYw3TpYw3CYw3TpYw3(3,3),      & 
& MXM3CYx3YdadjYdTpYx3(3,3),     & 
& MXM3CYx3TpYx3CYx3TpYx3(3,3),   & 
& Yb3ml2adjYb3Yb3adjYb3(3,3),Yb3ml2adjYeYeadjYb3(3,3),      & 
& Yb3ml2adjYw3Yw3adjYb3(3,3),Yb3adjYb3mHb32Yb3adjYb3(3,3),Yb3adjYb3Yb3ml2adjYb3(3,3),    & 
& Yb3adjYb3Yb3adjYb3MBM3(3,3),Yb3adjYb3Yb3adjYb3mHb32(3,3),Yb3adjYb3Yb3adjYb3Yb3(3,3),   & 
& Yb3adjYb3Yb3adjYb3TYb3(3,3),Yb3adjYb3Yb3adjYw3Yw3(3,3),   & 
& Yb3adjYb3Yb3adjYw3TYw3(3,3),Yb3adjYb3TYb3adjYb3Yb3(3,3),  & 
& Yb3adjYb3TYb3adjYw3Yw3(3,3),Yb3adjYeme2YeadjYb3(3,3),Yb3adjYeYeml2adjYb3(3,3),         & 
& Yb3adjYeYeadjYb3MBM3(3,3),Yb3adjYeYeadjYb3mHb32(3,3),Yb3adjYeYeadjYb3Yb3(3,3),         & 
& Yb3adjYeYeadjYb3TYb3(3,3),Yb3adjYeYeadjYeYe(3,3),           & 
& Yb3adjYeYeadjYeTYe(3,3),Yb3adjYeYeadjYw3TYw3(3,3),          & 
& Yb3adjYeTYeadjYb3Yb3(3,3),Yb3adjYeTYeadjYeYe(3,3),Yb3adjYeTYeadjYw3Yw3(3,3),           & 
& Yb3adjYw3mHw32Yw3adjYb3(3,3),Yb3adjYw3Yw3ml2adjYb3(3,3),Yb3adjYw3Yw3adjYb3MBM3(3,3),   & 
& Yb3adjYw3Yw3adjYb3mHb32(3,3),Yb3adjYw3Yw3adjYb3Yb3(3,3),  & 
& Yb3adjYw3Yw3adjYb3TYb3(3,3),Yb3adjYw3Yw3adjYw3Yw3(3,3),Yb3adjYw3Yw3adjYw3TYw3(3,3),    & 
& Yb3adjYw3TYw3adjYb3Yb3(3,3),Yb3adjYw3TYw3adjYw3Yw3(3,3),  & 
& Ydmq2adjYdYdadjYd(3,3),Ydmq2adjYuYuadjYd(3,3),Ydmq2adjYx3Yx3CYd(3,3),YdadjYdmd2YdadjYd(3,3),& 
& YdadjYdYdmq2adjYd(3,3),YdadjYdYdadjYdmd2(3,3),YdadjYdYdadjYdYd(3,3),YdadjYdYdadjYdTYd(3,3),& 
& YdadjYdTYdadjYdYd(3,3),YdadjYdTpYx3CYx3Yd(3,3),YdadjYdTpYx3CYx3TYd(3,3),               & 
& YdadjYdTpTYx3CYx3Yd(3,3),YdadjYumu2YuadjYd(3,3),YdadjYuYumq2adjYd(3,3),YdadjYuYuadjYdmd2(3,3),& 
& YdadjYuYuadjYdYd(3,3),YdadjYuYuadjYdTYd(3,3),YdadjYuYuadjYuYu(3,3),YdadjYuYuadjYuTYu(3,3),& 
& YdadjYuTYuadjYdYd(3,3),YdadjYuTYuadjYuYu(3,3),Yeml2adjYb3Yb3adjYe(3,3),Yeml2adjYeYeadjYe(3,3),& 
& Yeml2adjYw3Yw3adjYe(3,3),YeadjYb3mHb32Yb3adjYe(3,3),YeadjYb3Yb3ml2adjYe(3,3),          & 
& YeadjYb3Yb3adjYb3Yb3(3,3),YeadjYb3Yb3adjYb3TYb3(3,3),YeadjYb3Yb3adjYeme2(3,3),         & 
& YeadjYb3Yb3adjYeYe(3,3),YeadjYb3Yb3adjYeTYe(3,3),YeadjYb3Yb3adjYw3Yw3(3,3),            & 
& YeadjYb3Yb3adjYw3TYw3(3,3),YeadjYb3TYb3adjYb3Yb3(3,3),YeadjYb3TYb3adjYeYe(3,3),        & 
& YeadjYb3TYb3adjYw3Yw3(3,3),YeadjYeme2YeadjYe(3,3),YeadjYeYeml2adjYe(3,3),              & 
& YeadjYeYeadjYeme2(3,3),YeadjYeYeadjYeYe(3,3),YeadjYeYeadjYeTYe(3,3),YeadjYeTYeadjYeYe(3,3),& 
& YeadjYw3mHw32Yw3adjYe(3,3),YeadjYw3Yw3ml2adjYe(3,3),YeadjYw3Yw3adjYb3Yb3(3,3),         & 
& YeadjYw3Yw3adjYb3TYb3(3,3),YeadjYw3Yw3adjYeme2(3,3),YeadjYw3Yw3adjYeYe(3,3),           & 
& YeadjYw3Yw3adjYeTYe(3,3),YeadjYw3Yw3adjYw3Yw3(3,3),YeadjYw3Yw3adjYw3TYw3(3,3)

Complex(dp) :: YeadjYw3TYw3adjYb3Yb3(3,3),YeadjYw3TYw3adjYeYe(3,3),YeadjYw3TYw3adjYw3Yw3(3,3),        & 
& Yumq2adjYdYdadjYu(3,3),Yumq2adjYuYuadjYu(3,3),YuadjYdmd2YdadjYu(3,3),YuadjYdYdmq2adjYu(3,3),& 
& YuadjYdYdadjYdYd(3,3),YuadjYdYdadjYdTYd(3,3),YuadjYdYdadjYumu2(3,3),YuadjYdYdadjYuYu(3,3),& 
& YuadjYdYdadjYuTYu(3,3),YuadjYdTYdadjYdYd(3,3),YuadjYdTYdadjYuYu(3,3),YuadjYdTpYx3CYx3Yd(3,3),& 
& YuadjYdTpYx3CYx3TYd(3,3),YuadjYdTpTYx3CYx3Yd(3,3),YuadjYumu2YuadjYu(3,3),              & 
& YuadjYuYumq2adjYu(3,3),YuadjYuYuadjYumu2(3,3),YuadjYuYuadjYuYu(3,3),YuadjYuYuadjYuTYu(3,3),& 
& YuadjYuTYuadjYuYu(3,3),Yw3ml2adjYb3Yb3adjYw3(3,3),Yw3ml2adjYeYeadjYw3(3,3),            & 
& Yw3ml2adjYw3Yw3adjYw3(3,3),Yw3adjYb3mHb32Yb3adjYw3(3,3),Yw3adjYb3Yb3ml2adjYw3(3,3),    & 
& Yw3adjYb3Yb3adjYb3Yb3(3,3),Yw3adjYb3Yb3adjYb3TYb3(3,3),Yw3adjYb3Yb3adjYw3mHw32(3,3),   & 
& Yw3adjYb3Yb3adjYw3MWM3(3,3),Yw3adjYb3Yb3adjYw3Yw3(3,3),   & 
& Yw3adjYb3Yb3adjYw3TYw3(3,3),Yw3adjYb3TYb3adjYb3Yb3(3,3),  & 
& Yw3adjYb3TYb3adjYw3Yw3(3,3),Yw3adjYeme2YeadjYw3(3,3),Yw3adjYeYeml2adjYw3(3,3),         & 
& Yw3adjYeYeadjYb3TYb3(3,3),Yw3adjYeYeadjYeYe(3,3),Yw3adjYeYeadjYeTYe(3,3),              & 
& Yw3adjYeYeadjYw3mHw32(3,3),Yw3adjYeYeadjYw3MWM3(3,3),Yw3adjYeYeadjYw3Yw3(3,3),         & 
& Yw3adjYeYeadjYw3TYw3(3,3),Yw3adjYeTYeadjYb3Yb3(3,3),        & 
& Yw3adjYeTYeadjYeYe(3,3),Yw3adjYeTYeadjYw3Yw3(3,3),          & 
& Yw3adjYw3mHw32Yw3adjYw3(3,3),Yw3adjYw3Yw3ml2adjYw3(3,3),Yw3adjYw3Yw3adjYb3Yb3(3,3),    & 
& Yw3adjYw3Yw3adjYb3TYb3(3,3),Yw3adjYw3Yw3adjYw3mHw32(3,3),Yw3adjYw3Yw3adjYw3MWM3(3,3),  & 
& Yw3adjYw3Yw3adjYw3Yw3(3,3),Yw3adjYw3Yw3adjYw3TYw3(3,3),   & 
& Yw3adjYw3TYw3adjYb3Yb3(3,3),Yw3adjYw3TYw3adjYw3Yw3(3,3),  & 
& Yx3md2adjYx3Yx3adjYx3(3,3),Yx3md2CYdTpYdadjYx3(3,3),Yx3adjYx3mHxb32Yx3adjYx3(3,3),     & 
& Yx3adjYx3Yx3md2adjYx3(3,3),Yx3adjYx3Yx3adjYx3mHxb32(3,3),Yx3adjYx3Yx3adjYx3Yx3(3,3),   & 
& Yx3adjYx3Yx3adjYx3TYx3(3,3),Yx3adjYx3TYx3adjYx3Yx3(3,3),Yx3CYdmq2TpYdadjYx3(3,3),      & 
& Yx3CYdTpYdmd2adjYx3(3,3),Yx3CYdTpYdadjYx3mHxb32(3,3),Yx3CYdTpYdadjYx3Yx3(3,3),         & 
& Yx3CYdTpYdadjYx3TYx3(3,3),Yx3CYdTpYdCYdTpYd(3,3),Yx3CYdTpYdCYdTpTYd(3,3),              & 
& Yx3CYdTpYuCYuTpYd(3,3),Yx3CYdTpYuCYuTpTYd(3,3),Yx3CYdTpTYdCYdTpYd(3,3),Yx3CYdTpTYuCYuTpYd(3,3),& 
& Yx3TYdadjYdadjYx3Yx3(3,3),  & 
& TYb3adjYb3Yb3adjYb3Yb3(3,3),TYb3adjYb3Yb3adjYw3Yw3(3,3),  & 
& TYb3adjYeYeadjYb3Yb3(3,3),TYb3adjYeYeadjYeYe(3,3),          & 
& TYb3adjYeYeadjYw3Yw3(3,3),TYb3adjYw3Yw3adjYb3Yb3(3,3),    & 
& TYb3adjYw3Yw3adjYw3Yw3(3,3),TYdadjYdYdadjYdYd(3,3),TYdadjYdTpYx3CYx3Yd(3,3),           & 
& TYdadjYuYuadjYdYd(3,3),TYdadjYuYuadjYuYu(3,3),TYeadjYb3Yb3adjYb3Yb3(3,3),              & 
& TYeadjYb3Yb3adjYeYe(3,3),TYeadjYb3Yb3adjYw3Yw3(3,3),TYeadjYeYeadjYeYe(3,3),            & 
& TYeadjYw3Yw3adjYb3Yb3(3,3),TYeadjYw3Yw3adjYeYe(3,3),TYeadjYw3Yw3adjYw3Yw3(3,3),        & 
& TYuadjYdYdadjYdYd(3,3),TYuadjYdYdadjYuYu(3,3),TYuadjYdTpYx3CYx3Yd(3,3),TYuadjYuYuadjYuYu(3,3)

Complex(dp) :: TYw3adjYb3Yb3adjYb3Yb3(3,3),TYw3adjYb3Yb3adjYw3Yw3(3,3),  & 
& TYw3adjYeYeadjYb3Yb3(3,3),TYw3adjYeYeadjYeYe(3,3),          & 
& TYw3adjYeYeadjYw3Yw3(3,3),TYw3adjYw3Yw3adjYb3Yb3(3,3),    & 
& TYw3adjYw3Yw3adjYw3Yw3(3,3),TYx3adjYx3Yx3adjYx3Yx3(3,3),TYx3CYdTpYdadjYx3Yx3(3,3),     & 
& TYx3CYdTpYdCYdTpYd(3,3),TYx3CYdTpYuCYuTpYd(3,3),TpYb3mHb32CYb3TpYb3CYb3(3,3),          & 
& TpYb3mHb32CYb3TpYw3CYw3(3,3),TpYb3CYb3ml2TpYb3CYb3(3,3),TpYb3CYb3ml2TpYw3CYw3(3,3),    & 
& TpYb3CYb3TpYb3mHb32CYb3(3,3),TpYb3CYb3TpYb3CYb3ml2(3,3),TpYb3CYb3TpYw3mHw32CYw3(3,3),  & 
& TpYb3CYb3TpYw3CYw3ml2(3,3),TpYdmd2adjYx3Yx3CYd(3,3),TpYdmd2CYdTpYdCYd(3,3),            & 
& TpYdadjYx3mHxb32Yx3CYd(3,3),TpYdadjYx3Yx3md2CYd(3,3),TpYdadjYx3Yx3CYdmq2(3,3),         & 
& TpYdCYdmq2TpYdCYd(3,3),TpYdCYdTpYdmd2CYd(3,3),TpYdCYdTpYdCYdmq2(3,3),TpYeme2CYeTpYeCYe(3,3),& 
& TpYeCYeml2TpYeCYe(3,3),TpYeCYeTpYeme2CYe(3,3),TpYeCYeTpYeCYeml2(3,3),TpYumu2CYuTpYuCYu(3,3),& 
& TpYuCYumq2TpYuCYu(3,3),TpYuCYuTpYumu2CYu(3,3),TpYuCYuTpYuCYumq2(3,3),TpYw3mHw32CYw3TpYb3CYb3(3,3),& 
& TpYw3mHw32CYw3TpYw3CYw3(3,3),TpYw3CYw3ml2TpYb3CYb3(3,3),TpYw3CYw3ml2TpYw3CYw3(3,3),    & 
& TpYw3CYw3TpYb3mHb32CYb3(3,3),TpYw3CYw3TpYb3CYb3ml2(3,3),TpYw3CYw3TpYw3mHw32CYw3(3,3),  & 
& TpYw3CYw3TpYw3CYw3ml2(3,3),TpYx3mHxb32CYx3TpYx3CYx3(3,3),TpYx3CYx3md2TpYx3CYx3(3,3),   & 
& TpYx3CYx3TpYx3mHxb32CYx3(3,3),TpYx3CYx3TpYx3CYx3md2(3,3),TpYx3CYx3TpYx3CYx3Yd(3,3),    & 
& TpYx3CYx3TpYx3CYx3TYd(3,3),TpYx3CYx3TpTYx3CYx3Yd(3,3),TpTYx3CYx3TpYx3CYx3Yd(3,3)

Complex(dp) :: Trmd2,Trme2,TrmHx32,TrmHxb32,Trml2,Trmq2,Trmu2,TrYb3adjYb3,TrYdadjYd,TrYeadjYe,       & 
& TrYuadjYu,TrYw3adjYw3,TrYx3adjYx3,TradjYb3TYb3,TradjYdTYd,TradjYeTYe,TradjYuTYu,       & 
& TradjYw3TYw3,TradjYx3TYx3,TrCTYb3TpTYb3,TrCTYdTpTYd,TrCTYeTpTYe,TrCTYuTpTYu,           & 
& TrCTYw3TpTYw3,TrCTYx3TpTYx3,Trmd2YdadjYd,Trmd2adjYx3Yx3,Trme2YeadjYe,TrmHb32Yb3adjYb3, & 
& TrmHw32Yw3adjYw3,TrmHxb32Yx3adjYx3,Trml2adjYb3Yb3,Trml2adjYeYe,Trml2adjYw3Yw3,         & 
& Trmq2adjYdYd,Trmq2adjYuYu,Trmu2YuadjYu

Complex(dp) :: TrmHg32,TrmHw32,TrCTYb3TpYb3,TrCTYdTpYd,TrCTYeTpYe,TrCTYuTpYu,TrCTYw3TpYw3,           & 
& TrCTYx3TpYx3,TrYb3adjYb3Yb3adjYb3,TrYb3adjYb3TYb3adjYb3,TrYb3adjYb3TYb3adjTYb3,        & 
& TrYb3adjYeYeadjYb3,TrYb3adjYeTYeadjYb3,TrYb3adjYeTYeadjTYb3,TrYb3adjYw3Yw3adjYb3,      & 
& TrYb3adjYw3TYw3adjYb3,TrYb3adjYw3TYw3adjTYb3,TrYb3adjTYb3TYb3adjYb3,TrYb3adjTYeTYeadjYb3,& 
& TrYb3adjTYw3TYw3adjYb3,TrYb3TpTYb3CTYb3adjYb3,TrYb3TpTYeCTYeadjYb3,TrYb3TpTYw3CTYw3adjYb3,& 
& TrYdadjYdYdadjYd,TrYdadjYdTYdadjYd,TrYdadjYdTYdadjTYd,TrYdadjYdTpYx3CYx3,              & 
& TrYdadjYdTpTYx3CTYx3,TrYdadjYuYuadjYd,TrYdadjYuTYuadjYd,TrYdadjYuTYuadjTYd,            & 
& TrYdadjTYdTYdadjYd,TrYdadjTYuTYuadjYd,TrYdTpTYdCTYdadjYd,TrYdTpTYuCTYuadjYd,           & 
& TrYeadjYb3TYb3adjYe,TrYeadjYb3TYb3adjTYe,TrYeadjYeYeadjYe,TrYeadjYeTYeadjYe,           & 
& TrYeadjYeTYeadjTYe,TrYeadjYw3Yw3adjYe,TrYeadjYw3TYw3adjYe,TrYeadjYw3TYw3adjTYe,        & 
& TrYeadjTYb3TYb3adjYe,TrYeadjTYeTYeadjYe,TrYeadjTYw3TYw3adjYe,TrYeTpTYb3CTYb3adjYe,     & 
& TrYeTpTYeCTYeadjYe,TrYeTpTYw3CTYw3adjYe,TrYuadjYdTYdadjYu,TrYuadjYdTYdadjTYu,          & 
& TrYuadjYuYuadjYu,TrYuadjYuTYuadjYu,TrYuadjYuTYuadjTYu,TrYuadjTYdTYdadjYu,              & 
& TrYuadjTYuTYuadjYu,TrYuTpTYdCTYdadjYu,TrYuTpTYuCTYuadjYu,TrYw3adjYb3TYb3adjYw3,        & 
& TrYw3adjYb3TYb3adjTYw3,TrYw3adjYeTYeadjYw3,TrYw3adjYeTYeadjTYw3,TrYw3adjYw3Yw3adjYw3,  & 
& TrYw3adjYw3TYw3adjYw3,TrYw3adjYw3TYw3adjTYw3,TrYw3adjTYb3TYb3adjYw3,TrYw3adjTYeTYeadjYw3,& 
& TrYw3adjTYw3TYw3adjYw3,TrYw3TpTYb3CTYb3adjYw3,TrYw3TpTYeCTYeadjYw3,TrYw3TpTYw3CTYw3adjYw3,& 
& TrYx3adjYx3Yx3adjYx3,TrYx3adjYx3TYx3adjYx3,TrYx3adjYx3TYx3adjTYx3,TrYx3adjTYx3TYx3adjYx3,& 
& TrYx3CTYdTpTYdadjYx3,TrYx3TpTYx3CTYx3adjYx3,TradjYdTpYx3CYx3TYd,TradjYdTpYx3CTYx3TYd,  & 
& TradjYdTpYx3CTYx3TpTYd,TradjYdTpTYx3CTYx3TpYd,TradjYx3TYx3CYdTpYd,TradjYx3TYx3CTYdTpYd,& 
& Trmd2YdadjYdYdadjYd,Trmd2YdadjYdTpYx3CYx3,Trmd2YdadjYuYuadjYd,Trmd2adjYx3Yx3adjYx3Yx3, & 
& Trmd2TpYx3CYx3YdadjYd,Trme2YeadjYb3Yb3adjYe,Trme2YeadjYeYeadjYe,Trme2YeadjYw3Yw3adjYe, & 
& TrmHb32Yb3adjYb3Yb3adjYb3,TrmHb32Yb3adjYeYeadjYb3,TrmHb32Yb3adjYw3Yw3adjYb3,           & 
& TrmHw32Yw3adjYb3Yb3adjYw3,TrmHw32Yw3adjYeYeadjYw3,TrmHw32Yw3adjYw3Yw3adjYw3,           & 
& TrmHxb32Yx3adjYx3Yx3adjYx3,TrmHxb32Yx3CYdTpYdadjYx3,TrmHxb32CYx3YdadjYdTpYx3,          & 
& Trml2adjYb3Yb3adjYb3Yb3,Trml2adjYb3Yb3adjYeYe,Trml2adjYb3Yb3adjYw3Yw3,Trml2adjYeYeadjYb3Yb3,& 
& Trml2adjYeYeadjYeYe,Trml2adjYeYeadjYw3Yw3,Trml2adjYw3Yw3adjYb3Yb3,Trml2adjYw3Yw3adjYeYe,& 
& Trml2adjYw3Yw3adjYw3Yw3,Trmq2adjYdYdadjYdYd,Trmq2adjYdYdadjYuYu,Trmq2adjYdTpYx3CYx3Yd, & 
& Trmq2adjYuYuadjYdYd,Trmq2adjYuYuadjYuYu,Trmu2YuadjYdYdadjYu,Trmu2YuadjYuYuadjYu

Complex(dp) :: MnuL(3,3), DMnuL(3,3), betaMnuL1(3,3), betaMNuL2(3,3), gammaL1(3,3), &
 & gammaL2(3,3), gammaHd1, gammaHd2

Real(dp) :: g1p2,g1p3,g2p2,g2p3,g3p2,g3p3

Real(dp) :: g1p4,g2p4,g3p4

Iname = Iname +1 
NameOfUnit(Iname) = 'rge555' 
 
OnlyDiagonal = .Not.GenerationMixing 
q = t 
 
Call GToParameters555(gy,g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3,mue,MXM3,MWM3,MGM3,            & 
& MBM3,TYu,TYd,TYe,TYb3,TYw3,TYx3,Bmue,BMXM3,BMWM3,BMGM3,BMBM3,mq2,ml2,mHd2,             & 
& mHu2,md2,mu2,me2,mHw32,mHg32,mHb32,mHx32,mHxb32,MassB,MassWB,MassG,MnuL)

If (ThresholdCrossed.lt.1) Then 
MXM3(1,:) = 0._dp 
BMXM3(1,:) = 0._dp 
mHx32(1,:) = 0._dp 
mHx32(:,1) = 0._dp 
Yx3(1,:) = 0._dp 
TYx3(1,:) = 0._dp 
MXM3(:,1) = 0._dp 
BMXM3(:,1) = 0._dp 
mHxb32(1,:) = 0._dp 
mHxb32(:,1) = 0._dp 
MGM3(1,:) = 0._dp 
BMGM3(1,:) = 0._dp 
MGM3(:,1) = 0._dp 
BMGM3(:,1) = 0._dp 
mHg32(1,:) = 0._dp 
mHg32(:,1) = 0._dp 
Yb3(1,:) = 0._dp 
TYb3(1,:) = 0._dp 
MBM3(1,:) = 0._dp 
BMBM3(1,:) = 0._dp 
MBM3(:,1) = 0._dp 
BMBM3(:,1) = 0._dp 
mHb32(1,:) = 0._dp 
mHb32(:,1) = 0._dp 
Yw3(1,:) = 0._dp 
TYw3(1,:) = 0._dp 
MWM3(1,:) = 0._dp 
BMWM3(1,:) = 0._dp 
MWM3(:,1) = 0._dp 
BMWM3(:,1) = 0._dp 
mHw32(1,:) = 0._dp 
mHw32(:,1) = 0._dp 
End if 

If (ThresholdCrossed.lt.2) Then 
MXM3(2,:) = 0._dp 
BMXM3(2,:) = 0._dp 
mHx32(2,:) = 0._dp 
mHx32(:,2) = 0._dp 
Yx3(2,:) = 0._dp 
TYx3(2,:) = 0._dp 
MXM3(:,2) = 0._dp 
BMXM3(:,2) = 0._dp 
mHxb32(2,:) = 0._dp 
mHxb32(:,2) = 0._dp 
MGM3(2,:) = 0._dp 
BMGM3(2,:) = 0._dp 
MGM3(:,2) = 0._dp 
BMGM3(:,2) = 0._dp 
mHg32(2,:) = 0._dp 
mHg32(:,2) = 0._dp 
Yb3(2,:) = 0._dp 
TYb3(2,:) = 0._dp 
MBM3(2,:) = 0._dp 
BMBM3(2,:) = 0._dp 
MBM3(:,2) = 0._dp 
BMBM3(:,2) = 0._dp 
mHb32(2,:) = 0._dp 
mHb32(:,2) = 0._dp 
Yw3(2,:) = 0._dp 
TYw3(2,:) = 0._dp 
MWM3(2,:) = 0._dp 
BMWM3(2,:) = 0._dp 
MWM3(:,2) = 0._dp 
BMWM3(:,2) = 0._dp 
mHw32(2,:) = 0._dp 
mHw32(:,2) = 0._dp 
End if 

If (ThresholdCrossed.lt.3) Then 
MXM3(3,:) = 0._dp 
BMXM3(3,:) = 0._dp 
mHx32(3,:) = 0._dp 
mHx32(:,3) = 0._dp 
Yx3(3,:) = 0._dp 
TYx3(3,:) = 0._dp 
MXM3(:,3) = 0._dp 
BMXM3(:,3) = 0._dp 
mHxb32(3,:) = 0._dp 
mHxb32(:,3) = 0._dp 
MGM3(3,:) = 0._dp 
BMGM3(3,:) = 0._dp 
MGM3(:,3) = 0._dp 
BMGM3(:,3) = 0._dp 
mHg32(3,:) = 0._dp 
mHg32(:,3) = 0._dp 
Yb3(3,:) = 0._dp 
TYb3(3,:) = 0._dp 
MBM3(3,:) = 0._dp 
BMBM3(3,:) = 0._dp 
MBM3(:,3) = 0._dp 
BMBM3(:,3) = 0._dp 
mHb32(3,:) = 0._dp 
mHb32(:,3) = 0._dp 
Yw3(3,:) = 0._dp 
TYw3(3,:) = 0._dp 
MWM3(3,:) = 0._dp 
BMWM3(3,:) = 0._dp 
MWM3(:,3) = 0._dp 
BMWM3(:,3) = 0._dp 
mHw32(3,:) = 0._dp 
mHw32(:,3) = 0._dp 
End if 

AbsMassB = Abs(MassB)**2
AbsMassWB = Abs(MassWB)**2
AbsMassG = Abs(MassG)**2
Call Adjungate(Yu,adjYu)
Call Adjungate(Yd,adjYd)
Call Adjungate(Ye,adjYe)
Call Adjungate(Yb3,adjYb3)
Call Adjungate(Yw3,adjYw3)
Call Adjungate(Yx3,adjYx3)
Call Adjungate(MXM3,adjMXM3)
Call Adjungate(MWM3,adjMWM3)
Call Adjungate(MGM3,adjMGM3)
Call Adjungate(MBM3,adjMBM3)
Call Adjungate(TYu,adjTYu)
Call Adjungate(TYd,adjTYd)
Call Adjungate(TYe,adjTYe)
Call Adjungate(TYb3,adjTYb3)
Call Adjungate(TYw3,adjTYw3)
Call Adjungate(TYx3,adjTYx3)
!Call Adjungate(BMXM3,adjBMXM3)
!Call Adjungate(BMWM3,adjBMWM3)
!Call Adjungate(BMGM3,adjBMGM3)
!Call Adjungate(BMBM3,adjBMBM3)
Call Adjungate(mq2,adjmq2)
Call Adjungate(ml2,adjml2)
Call Adjungate(md2,adjmd2)
Call Adjungate(mu2,adjmu2)
Call Adjungate(me2,adjme2)
Call Adjungate(mHw32,adjmHw32)
Call Adjungate(mHg32,adjmHg32)
Call Adjungate(mHb32,adjmHb32)
Call Adjungate(mHx32,adjmHx32)
Call Adjungate(mHxb32,adjmHxb32)
 md2adjYx3 = Matmul(md2,adjYx3) 
 md2CYd = Matmul(md2,Conjg(Yd)) 
 me2CYe = Matmul(me2,Conjg(Ye)) 
 mHb32CYb3 = Matmul(mHb32,Conjg(Yb3)) 
 mHw32CYw3 = Matmul(mHw32,Conjg(Yw3)) 
 mHxb32CYx3 = Matmul(mHxb32,Conjg(Yx3)) 
 ml2adjYb3 = Matmul(ml2,adjYb3) 
 ml2adjYe = Matmul(ml2,adjYe) 
 ml2adjYw3 = Matmul(ml2,adjYw3) 
 mq2adjYd = Matmul(mq2,adjYd) 
 mq2adjYu = Matmul(mq2,adjYu) 
 mu2CYu = Matmul(mu2,Conjg(Yu)) 
 Yb3adjYb3 = Matmul(Yb3,adjYb3) 
Forall(i2=1:3)  Yb3adjYb3(i2,i2) =  Real(Yb3adjYb3(i2,i2),dp) 
 YdadjYd = Matmul(Yd,adjYd) 
Forall(i2=1:3)  YdadjYd(i2,i2) =  Real(YdadjYd(i2,i2),dp) 
 YeadjYe = Matmul(Ye,adjYe) 
Forall(i2=1:3)  YeadjYe(i2,i2) =  Real(YeadjYe(i2,i2),dp) 
 YuadjYu = Matmul(Yu,adjYu) 
Forall(i2=1:3)  YuadjYu(i2,i2) =  Real(YuadjYu(i2,i2),dp) 
 Yw3adjYw3 = Matmul(Yw3,adjYw3) 
Forall(i2=1:3)  Yw3adjYw3(i2,i2) =  Real(Yw3adjYw3(i2,i2),dp) 
 Yx3adjYx3 = Matmul(Yx3,adjYx3) 
Forall(i2=1:3)  Yx3adjYx3(i2,i2) =  Real(Yx3adjYx3(i2,i2),dp) 
 adjYb3MBM3 = Matmul(adjYb3,MBM3) 
 adjYb3mHb32 = Matmul(adjYb3,mHb32) 
 adjYb3Yb3 = Matmul(adjYb3,Yb3) 
Forall(i2=1:3)  adjYb3Yb3(i2,i2) =  Real(adjYb3Yb3(i2,i2),dp) 
! adjYb3BMBM3 = Matmul(adjYb3,BMBM3) 
 adjYb3TYb3 = Matmul(adjYb3,TYb3) 
 adjYdmd2 = Matmul(adjYd,md2) 
 adjYdYd = Matmul(adjYd,Yd) 
Forall(i2=1:3)  adjYdYd(i2,i2) =  Real(adjYdYd(i2,i2),dp) 
 adjYdTYd = Matmul(adjYd,TYd) 
 adjYeme2 = Matmul(adjYe,me2) 
 adjYeYe = Matmul(adjYe,Ye) 
Forall(i2=1:3)  adjYeYe(i2,i2) =  Real(adjYeYe(i2,i2),dp) 
 adjYeTYe = Matmul(adjYe,TYe) 
 adjYumu2 = Matmul(adjYu,mu2) 
 adjYuYu = Matmul(adjYu,Yu) 
Forall(i2=1:3)  adjYuYu(i2,i2) =  Real(adjYuYu(i2,i2),dp) 
 adjYuTYu = Matmul(adjYu,TYu) 
 adjYw3mHw32 = Matmul(adjYw3,mHw32) 
 adjYw3MWM3 = Matmul(adjYw3,MWM3) 
 adjYw3Yw3 = Matmul(adjYw3,Yw3) 
Forall(i2=1:3)  adjYw3Yw3(i2,i2) =  Real(adjYw3Yw3(i2,i2),dp) 
! adjYw3BMWM3 = Matmul(adjYw3,BMWM3) 
 adjYw3TYw3 = Matmul(adjYw3,TYw3) 
 adjYx3mHxb32 = Matmul(adjYx3,mHxb32) 
 adjYx3Yx3 = Matmul(adjYx3,Yx3) 
Forall(i2=1:3)  adjYx3Yx3(i2,i2) =  Real(adjYx3Yx3(i2,i2),dp) 
 adjYx3TYx3 = Matmul(adjYx3,TYx3) 
 CYb3ml2 = Matmul(Conjg(Yb3),ml2) 
 CYb3TpYb3 = Matmul(Conjg(Yb3),Transpose(Yb3)) 
Forall(i2=1:3)  CYb3TpYb3(i2,i2) =  Real(CYb3TpYb3(i2,i2),dp) 
! CYb3TpTYb3 = Matmul(Conjg(Yb3),Transpose(TYb3)) 
 CYdmq2 = Matmul(Conjg(Yd),mq2) 
 CYdTpYd = Matmul(Conjg(Yd),Transpose(Yd)) 
Forall(i2=1:3)  CYdTpYd(i2,i2) =  Real(CYdTpYd(i2,i2),dp) 
 CYdTpTYd = Matmul(Conjg(Yd),Transpose(TYd)) 
 CYeml2 = Matmul(Conjg(Ye),ml2) 
 CYumq2 = Matmul(Conjg(Yu),mq2) 
 CYw3ml2 = Matmul(Conjg(Yw3),ml2) 
 CYw3TpYw3 = Matmul(Conjg(Yw3),Transpose(Yw3)) 
Forall(i2=1:3)  CYw3TpYw3(i2,i2) =  Real(CYw3TpYw3(i2,i2),dp) 
! CYw3TpTYw3 = Matmul(Conjg(Yw3),Transpose(TYw3)) 
 CYx3md2 = Matmul(Conjg(Yx3),md2) 
 CYx3Yd = Matmul(Conjg(Yx3),Yd) 
 CYx3TYd = Matmul(Conjg(Yx3),TYd) 
 CYx3TpYx3 = Matmul(Conjg(Yx3),Transpose(Yx3)) 
Forall(i2=1:3)  CYx3TpYx3(i2,i2) =  Real(CYx3TpYx3(i2,i2),dp) 
! CYx3TpTYx3 = Matmul(Conjg(Yx3),Transpose(TYx3)) 
 CTYb3TpTYb3 = Matmul(Conjg(TYb3),Transpose(TYb3)) 
 CTYdTpTYd = Matmul(Conjg(TYd),Transpose(TYd)) 
 CTYeTpTYe = Matmul(Conjg(TYe),Transpose(TYe)) 
 CTYuTpTYu = Matmul(Conjg(TYu),Transpose(TYu)) 
 CTYw3TpTYw3 = Matmul(Conjg(TYw3),Transpose(TYw3)) 
 CTYx3TpTYx3 = Matmul(Conjg(TYx3),Transpose(TYx3)) 
 TYb3adjTYb3 = Matmul(TYb3,adjTYb3) 
 TYdadjTYd = Matmul(TYd,adjTYd) 
 TYeadjTYe = Matmul(TYe,adjTYe) 
 TYuadjTYu = Matmul(TYu,adjTYu) 
 TYw3adjTYw3 = Matmul(TYw3,adjTYw3) 
 TYx3adjTYx3 = Matmul(TYx3,adjTYx3) 
 TpYb3CYb3 = Matmul(Transpose(Yb3),Conjg(Yb3)) 
Forall(i2=1:3)  TpYb3CYb3(i2,i2) =  Real(TpYb3CYb3(i2,i2),dp) 
 TpYdCYd = Matmul(Transpose(Yd),Conjg(Yd)) 
Forall(i2=1:3)  TpYdCYd(i2,i2) =  Real(TpYdCYd(i2,i2),dp) 
 TpYeCYe = Matmul(Transpose(Ye),Conjg(Ye)) 
Forall(i2=1:3)  TpYeCYe(i2,i2) =  Real(TpYeCYe(i2,i2),dp) 
 TpYuCYu = Matmul(Transpose(Yu),Conjg(Yu)) 
Forall(i2=1:3)  TpYuCYu(i2,i2) =  Real(TpYuCYu(i2,i2),dp) 
 TpYw3CYw3 = Matmul(Transpose(Yw3),Conjg(Yw3)) 
Forall(i2=1:3)  TpYw3CYw3(i2,i2) =  Real(TpYw3CYw3(i2,i2),dp) 
 TpYx3CYx3 = Matmul(Transpose(Yx3),Conjg(Yx3)) 
Forall(i2=1:3)  TpYx3CYx3(i2,i2) =  Real(TpYx3CYx3(i2,i2),dp) 
 TpTYb3CTYb3 = Matmul(Transpose(TYb3),Conjg(TYb3)) 
 TpTYdCTYd = Matmul(Transpose(TYd),Conjg(TYd)) 
 TpTYeCTYe = Matmul(Transpose(TYe),Conjg(TYe)) 
 TpTYuCTYu = Matmul(Transpose(TYu),Conjg(TYu)) 
 TpTYw3CTYw3 = Matmul(Transpose(TYw3),Conjg(TYw3)) 
 TpTYx3CTYx3 = Matmul(Transpose(TYx3),Conjg(TYx3)) 
 MBM3CYb3TpYb3 = Matmul(MBM3,CYb3TpYb3) 
! MBM3CYb3TpTYb3 = Matmul(MBM3,CYb3TpTYb3) 
 md2YdadjYd = Matmul(md2,YdadjYd) 
 md2adjYx3Yx3 = Matmul(md2,adjYx3Yx3) 
 md2TpYx3CYx3 = Matmul(md2,TpYx3CYx3) 
 me2YeadjYe = Matmul(me2,YeadjYe) 
 mHb32Yb3adjYb3 = Matmul(mHb32,Yb3adjYb3) 
 mHw32Yw3adjYw3 = Matmul(mHw32,Yw3adjYw3) 
 mHxb32Yx3adjYx3 = Matmul(mHxb32,Yx3adjYx3) 
 ml2adjYb3Yb3 = Matmul(ml2,adjYb3Yb3) 
 ml2adjYeYe = Matmul(ml2,adjYeYe) 
 ml2adjYw3Yw3 = Matmul(ml2,adjYw3Yw3) 
 ml2TpYb3CYb3 = Matmul(ml2,TpYb3CYb3) 
 ml2TpYeCYe = Matmul(ml2,TpYeCYe) 
 ml2TpYw3CYw3 = Matmul(ml2,TpYw3CYw3) 
 mq2adjYdYd = Matmul(mq2,adjYdYd) 
 mq2adjYuYu = Matmul(mq2,adjYuYu) 
 mq2TpYdCYd = Matmul(mq2,TpYdCYd) 
 mq2TpYuCYu = Matmul(mq2,TpYuCYu) 
 mu2YuadjYu = Matmul(mu2,YuadjYu) 
 MWM3CYw3TpYw3 = Matmul(MWM3,CYw3TpYw3) 
! MWM3CYw3TpTYw3 = Matmul(MWM3,CYw3TpTYw3) 
 MXM3CYx3TpYx3 = Matmul(MXM3,CYx3TpYx3) 
! MXM3CYx3TpTYx3 = Matmul(MXM3,CYx3TpTYx3) 
 Yb3ml2adjYb3 = Matmul(Yb3,ml2adjYb3) 
 Yb3adjYb3MBM3 = Matmul(Yb3,adjYb3MBM3) 
 Yb3adjYb3mHb32 = Matmul(Yb3,adjYb3mHb32) 
 Yb3adjYb3Yb3 = Matmul(Yb3,adjYb3Yb3) 
! Yb3adjYb3BMBM3 = Matmul(Yb3,adjYb3BMBM3) 
 Yb3adjYb3TYb3 = Matmul(Yb3,adjYb3TYb3) 
 Yb3adjYeYe = Matmul(Yb3,adjYeYe) 
 Yb3adjYeTYe = Matmul(Yb3,adjYeTYe) 
 Yb3adjYw3Yw3 = Matmul(Yb3,adjYw3Yw3) 
 Yb3adjYw3TYw3 = Matmul(Yb3,adjYw3TYw3) 
 Ydmq2adjYd = Matmul(Yd,mq2adjYd) 
 YdadjYdmd2 = Matmul(Yd,adjYdmd2) 
 YdadjYdYd = Matmul(Yd,adjYdYd) 
 YdadjYdTYd = Matmul(Yd,adjYdTYd) 
 YdadjYuYu = Matmul(Yd,adjYuYu) 
 YdadjYuTYu = Matmul(Yd,adjYuTYu) 
 Yeml2adjYe = Matmul(Ye,ml2adjYe) 
 YeadjYb3Yb3 = Matmul(Ye,adjYb3Yb3) 
 YeadjYb3TYb3 = Matmul(Ye,adjYb3TYb3) 
 YeadjYeme2 = Matmul(Ye,adjYeme2) 
 YeadjYeYe = Matmul(Ye,adjYeYe) 
 YeadjYeTYe = Matmul(Ye,adjYeTYe) 
 YeadjYw3Yw3 = Matmul(Ye,adjYw3Yw3) 
 YeadjYw3TYw3 = Matmul(Ye,adjYw3TYw3) 
 Yumq2adjYu = Matmul(Yu,mq2adjYu) 
 YuadjYdYd = Matmul(Yu,adjYdYd) 
 YuadjYdTYd = Matmul(Yu,adjYdTYd) 
 YuadjYumu2 = Matmul(Yu,adjYumu2) 
 YuadjYuYu = Matmul(Yu,adjYuYu) 
 YuadjYuTYu = Matmul(Yu,adjYuTYu) 
 Yw3ml2adjYw3 = Matmul(Yw3,ml2adjYw3) 
 Yw3adjYb3Yb3 = Matmul(Yw3,adjYb3Yb3) 
 Yw3adjYb3TYb3 = Matmul(Yw3,adjYb3TYb3) 
 Yw3adjYeYe = Matmul(Yw3,adjYeYe) 
 Yw3adjYeTYe = Matmul(Yw3,adjYeTYe) 
 Yw3adjYw3mHw32 = Matmul(Yw3,adjYw3mHw32) 
 Yw3adjYw3MWM3 = Matmul(Yw3,adjYw3MWM3) 
 Yw3adjYw3Yw3 = Matmul(Yw3,adjYw3Yw3) 
! Yw3adjYw3BMWM3 = Matmul(Yw3,adjYw3BMWM3) 
 Yw3adjYw3TYw3 = Matmul(Yw3,adjYw3TYw3) 
 Yx3md2adjYx3 = Matmul(Yx3,md2adjYx3) 
 Yx3adjYx3mHxb32 = Matmul(Yx3,adjYx3mHxb32) 
 Yx3adjYx3Yx3 = Matmul(Yx3,adjYx3Yx3) 
 Yx3adjYx3TYx3 = Matmul(Yx3,adjYx3TYx3) 
 Yx3CYdTpYd = Matmul(Yx3,CYdTpYd) 
 Yx3CYdTpTYd = Matmul(Yx3,CYdTpTYd) 
! BMBM3CYb3TpYb3 = Matmul(BMBM3,CYb3TpYb3) 
! BMWM3CYw3TpYw3 = Matmul(BMWM3,CYw3TpYw3) 
! BMXM3CYx3TpYx3 = Matmul(BMXM3,CYx3TpYx3) 
! TYb3adjYb3MBM3 = Matmul(TYb3,adjYb3MBM3) 
 TYb3adjYb3Yb3 = Matmul(TYb3,adjYb3Yb3) 
 TYb3adjYeYe = Matmul(TYb3,adjYeYe) 
 TYb3adjYw3Yw3 = Matmul(TYb3,adjYw3Yw3) 
 TYdadjYdYd = Matmul(TYd,adjYdYd) 
 TYdadjYuYu = Matmul(TYd,adjYuYu) 
 TYeadjYb3Yb3 = Matmul(TYe,adjYb3Yb3) 
 TYeadjYeYe = Matmul(TYe,adjYeYe) 
 TYeadjYw3Yw3 = Matmul(TYe,adjYw3Yw3) 
 TYuadjYdYd = Matmul(TYu,adjYdYd) 
 TYuadjYuYu = Matmul(TYu,adjYuYu) 
 TYw3adjYb3Yb3 = Matmul(TYw3,adjYb3Yb3) 
 TYw3adjYeYe = Matmul(TYw3,adjYeYe) 
! TYw3adjYw3MWM3 = Matmul(TYw3,adjYw3MWM3) 
 TYw3adjYw3Yw3 = Matmul(TYw3,adjYw3Yw3) 
 TYx3adjYx3Yx3 = Matmul(TYx3,adjYx3Yx3) 
 TYx3CYdTpYd = Matmul(TYx3,CYdTpYd) 
 TpYb3mHb32CYb3 = Matmul(Transpose(Yb3),mHb32CYb3) 
 TpYb3CYb3ml2 = Matmul(Transpose(Yb3),CYb3ml2) 
 TpYdmd2CYd = Matmul(Transpose(Yd),md2CYd) 
 TpYdCYdmq2 = Matmul(Transpose(Yd),CYdmq2) 
 TpYeme2CYe = Matmul(Transpose(Ye),me2CYe) 
 TpYeCYeml2 = Matmul(Transpose(Ye),CYeml2) 
 TpYumu2CYu = Matmul(Transpose(Yu),mu2CYu) 
 TpYuCYumq2 = Matmul(Transpose(Yu),CYumq2) 
 TpYw3mHw32CYw3 = Matmul(Transpose(Yw3),mHw32CYw3) 
 TpYw3CYw3ml2 = Matmul(Transpose(Yw3),CYw3ml2) 
 TpYx3mHxb32CYx3 = Matmul(Transpose(Yx3),mHxb32CYx3) 
 TpYx3CYx3md2 = Matmul(Transpose(Yx3),CYx3md2) 
 TpYx3CYx3Yd = Matmul(Transpose(Yx3),CYx3Yd) 
 TpYx3CYx3TYd = Matmul(Transpose(Yx3),CYx3TYd) 
 TpTYx3CYx3Yd = Matmul(Transpose(TYx3),CYx3Yd) 
 Trmd2 = Real(cTrace(md2),dp) 
 Trme2 = Real(cTrace(me2),dp) 
 TrmHx32 = Real(cTrace(mHx32),dp) 
 TrmHxb32 = Real(cTrace(mHxb32),dp) 
 Trml2 = Real(cTrace(ml2),dp) 
 Trmq2 = Real(cTrace(mq2),dp) 
 Trmu2 = Real(cTrace(mu2),dp) 
 TrYb3adjYb3 = Real(cTrace(Yb3adjYb3),dp) 
 TrYdadjYd = Real(cTrace(YdadjYd),dp) 
 TrYeadjYe = Real(cTrace(YeadjYe),dp) 
 TrYuadjYu = Real(cTrace(YuadjYu),dp) 
 TrYw3adjYw3 = Real(cTrace(Yw3adjYw3),dp) 
 TrYx3adjYx3 = Real(cTrace(Yx3adjYx3),dp) 
 TradjYb3TYb3 = Real(cTrace(adjYb3TYb3),dp) 
 TradjYdTYd = Real(cTrace(adjYdTYd),dp) 
 TradjYeTYe = Real(cTrace(adjYeTYe),dp) 
 TradjYuTYu = Real(cTrace(adjYuTYu),dp) 
 TradjYw3TYw3 = Real(cTrace(adjYw3TYw3),dp) 
 TradjYx3TYx3 = Real(cTrace(adjYx3TYx3),dp) 
 TrCTYb3TpTYb3 = Real(cTrace(CTYb3TpTYb3),dp) 
 TrCTYdTpTYd = Real(cTrace(CTYdTpTYd),dp) 
 TrCTYeTpTYe = Real(cTrace(CTYeTpTYe),dp) 
 TrCTYuTpTYu = Real(cTrace(CTYuTpTYu),dp) 
 TrCTYw3TpTYw3 = Real(cTrace(CTYw3TpTYw3),dp) 
 TrCTYx3TpTYx3 = Real(cTrace(CTYx3TpTYx3),dp) 
 Trmd2YdadjYd = Real(cTrace(md2YdadjYd),dp) 
 Trmd2adjYx3Yx3 = Real(cTrace(md2adjYx3Yx3),dp) 
 Trme2YeadjYe = Real(cTrace(me2YeadjYe),dp) 
 TrmHb32Yb3adjYb3 = Real(cTrace(mHb32Yb3adjYb3),dp) 
 TrmHw32Yw3adjYw3 = Real(cTrace(mHw32Yw3adjYw3),dp) 
 TrmHxb32Yx3adjYx3 = Real(cTrace(mHxb32Yx3adjYx3),dp) 
 Trml2adjYb3Yb3 = Real(cTrace(ml2adjYb3Yb3),dp) 
 Trml2adjYeYe = Real(cTrace(ml2adjYeYe),dp) 
 Trml2adjYw3Yw3 = Real(cTrace(ml2adjYw3Yw3),dp) 
 Trmq2adjYdYd = Real(cTrace(mq2adjYdYd),dp) 
 Trmq2adjYuYu = Real(cTrace(mq2adjYuYu),dp) 
 Trmu2YuadjYu = Real(cTrace(mu2YuadjYu),dp) 
 g1p2 =g1**2 
 g1p3 =g1**3 
 g2p2 =g2**2 
 g2p3 =g2**3 
 g3p2 =g3**2 
 g3p3 =g3**3 


If (TwoLoopRGE) Then 
 Yb3adjYe = Matmul(Yb3,adjYe) 
 Yb3adjYw3 = Matmul(Yb3,adjYw3) 
 Yb3adjTYb3 = Matmul(Yb3,adjTYb3) 
 Yb3adjTYe = Matmul(Yb3,adjTYe) 
 Yb3adjTYw3 = Matmul(Yb3,adjTYw3) 
 YdadjYu = Matmul(Yd,adjYu) 
 YdadjTYd = Matmul(Yd,adjTYd) 
 YdadjTYu = Matmul(Yd,adjTYu) 
 YeadjYb3 = Matmul(Ye,adjYb3) 
 YeadjYw3 = Matmul(Ye,adjYw3) 
 YeadjTYb3 = Matmul(Ye,adjTYb3) 
 YeadjTYe = Matmul(Ye,adjTYe) 
 YeadjTYw3 = Matmul(Ye,adjTYw3) 
 YuadjYd = Matmul(Yu,adjYd) 
 YuadjTYd = Matmul(Yu,adjTYd) 
 YuadjTYu = Matmul(Yu,adjTYu) 
 Yw3adjYb3 = Matmul(Yw3,adjYb3) 
 Yw3adjYe = Matmul(Yw3,adjYe) 
 Yw3adjTYb3 = Matmul(Yw3,adjTYb3) 
 Yw3adjTYe = Matmul(Yw3,adjTYe) 
 Yw3adjTYw3 = Matmul(Yw3,adjTYw3) 
 Yx3adjTYx3 = Matmul(Yx3,adjTYx3) 
 Yx3CYd = Matmul(Yx3,Conjg(Yd)) 
 Yx3CTYd = Matmul(Yx3,Conjg(TYd)) 
 adjYdadjTYx3 = Matmul(adjYd,adjTYx3) 
 adjYdTpYx3 = Matmul(adjYd,Transpose(Yx3)) 
! adjYdTpTYx3 = Matmul(adjYd,Transpose(TYx3)) 
 adjTYdadjYx3 = Matmul(adjTYd,adjYx3) 
 CYb3TpYw3 = Matmul(Conjg(Yb3),Transpose(Yw3)) 
! CYb3TpTYw3 = Matmul(Conjg(Yb3),Transpose(TYw3)) 
 CYeTpYb3 = Matmul(Conjg(Ye),Transpose(Yb3)) 
 CYeTpYw3 = Matmul(Conjg(Ye),Transpose(Yw3)) 
! CYeTpTYb3 = Matmul(Conjg(Ye),Transpose(TYb3)) 
! CYeTpTYw3 = Matmul(Conjg(Ye),Transpose(TYw3)) 
 CYuTpYd = Matmul(Conjg(Yu),Transpose(Yd)) 
 CYuTpTYd = Matmul(Conjg(Yu),Transpose(TYd)) 
 CYw3TpYb3 = Matmul(Conjg(Yw3),Transpose(Yb3)) 
! CYw3TpTYb3 = Matmul(Conjg(Yw3),Transpose(TYb3)) 
 CTYb3adjYb3 = Matmul(Conjg(TYb3),adjYb3) 
 CTYb3adjYe = Matmul(Conjg(TYb3),adjYe) 
 CTYb3adjYw3 = Matmul(Conjg(TYb3),adjYw3) 
 CTYb3TpYb3 = Matmul(Conjg(TYb3),Transpose(Yb3)) 
 CTYdadjYd = Matmul(Conjg(TYd),adjYd) 
 CTYdadjYu = Matmul(Conjg(TYd),adjYu) 
 CTYdTpYd = Matmul(Conjg(TYd),Transpose(Yd)) 
 CTYeadjYb3 = Matmul(Conjg(TYe),adjYb3) 
 CTYeadjYe = Matmul(Conjg(TYe),adjYe) 
 CTYeadjYw3 = Matmul(Conjg(TYe),adjYw3) 
 CTYeTpYe = Matmul(Conjg(TYe),Transpose(Ye)) 
 CTYuadjYd = Matmul(Conjg(TYu),adjYd) 
 CTYuadjYu = Matmul(Conjg(TYu),adjYu) 
 CTYuTpYu = Matmul(Conjg(TYu),Transpose(Yu)) 
 CTYw3adjYb3 = Matmul(Conjg(TYw3),adjYb3) 
 CTYw3adjYe = Matmul(Conjg(TYw3),adjYe) 
 CTYw3adjYw3 = Matmul(Conjg(TYw3),adjYw3) 
 CTYw3TpYw3 = Matmul(Conjg(TYw3),Transpose(Yw3)) 
 CTYx3adjYx3 = Matmul(Conjg(TYx3),adjYx3) 
 CTYx3TYd = Matmul(Conjg(TYx3),TYd) 
 CTYx3TpYd = Matmul(Conjg(TYx3),Transpose(Yd)) 
 CTYx3TpYx3 = Matmul(Conjg(TYx3),Transpose(Yx3)) 
 CTYx3TpTYd = Matmul(Conjg(TYx3),Transpose(TYd)) 
 TYb3adjYb3 = Matmul(TYb3,adjYb3) 
 TYb3adjYe = Matmul(TYb3,adjYe) 
 TYb3adjYw3 = Matmul(TYb3,adjYw3) 
 TYb3adjTYe = Matmul(TYb3,adjTYe) 
 TYb3adjTYw3 = Matmul(TYb3,adjTYw3) 
 TYdadjYd = Matmul(TYd,adjYd) 
 TYdadjYu = Matmul(TYd,adjYu) 
 TYdadjTYu = Matmul(TYd,adjTYu) 
 TYeadjYb3 = Matmul(TYe,adjYb3) 
 TYeadjYe = Matmul(TYe,adjYe) 
 TYeadjYw3 = Matmul(TYe,adjYw3) 
 TYeadjTYb3 = Matmul(TYe,adjTYb3) 
 TYeadjTYw3 = Matmul(TYe,adjTYw3) 
 TYuadjYd = Matmul(TYu,adjYd) 
 TYuadjYu = Matmul(TYu,adjYu) 
 TYuadjTYd = Matmul(TYu,adjTYd) 
 TYw3adjYb3 = Matmul(TYw3,adjYb3) 
 TYw3adjYe = Matmul(TYw3,adjYe) 
 TYw3adjYw3 = Matmul(TYw3,adjYw3) 
 TYw3adjTYb3 = Matmul(TYw3,adjTYb3) 
 TYw3adjTYe = Matmul(TYw3,adjTYe) 
 TYx3adjYx3 = Matmul(TYx3,adjYx3) 
 TYx3CTYd = Matmul(TYx3,Conjg(TYd)) 
 TpYb3CTYb3 = Matmul(Transpose(Yb3),Conjg(TYb3)) 
 TpYdadjYx3 = Matmul(Transpose(Yd),adjYx3) 
 TpYdadjTYx3 = Matmul(Transpose(Yd),adjTYx3) 
 TpYdCTYd = Matmul(Transpose(Yd),Conjg(TYd)) 
 TpYeCTYe = Matmul(Transpose(Ye),Conjg(TYe)) 
 TpYuCTYu = Matmul(Transpose(Yu),Conjg(TYu)) 
 TpYw3CTYw3 = Matmul(Transpose(Yw3),Conjg(TYw3)) 
 TpYx3CTYx3 = Matmul(Transpose(Yx3),Conjg(TYx3)) 
 TpTYb3CYb3 = Matmul(Transpose(TYb3),Conjg(Yb3)) 
 TpTYdadjYx3 = Matmul(Transpose(TYd),adjYx3) 
 TpTYdCYd = Matmul(Transpose(TYd),Conjg(Yd)) 
 TpTYeCYe = Matmul(Transpose(TYe),Conjg(Ye)) 
 TpTYuCYu = Matmul(Transpose(TYu),Conjg(Yu)) 
 TpTYw3CYw3 = Matmul(Transpose(TYw3),Conjg(Yw3)) 
 TpTYx3CYx3 = Matmul(Transpose(TYx3),Conjg(Yx3)) 
 md2YdadjYu = Matmul(md2,YdadjYu) 
 me2YeadjYb3 = Matmul(me2,YeadjYb3) 
 me2YeadjYw3 = Matmul(me2,YeadjYw3) 
 mHb32Yb3adjYe = Matmul(mHb32,Yb3adjYe) 
 mHb32Yb3adjYw3 = Matmul(mHb32,Yb3adjYw3) 
 mHw32Yw3adjYb3 = Matmul(mHw32,Yw3adjYb3) 
 mHw32Yw3adjYe = Matmul(mHw32,Yw3adjYe) 
 mHxb32Yx3CYd = Matmul(mHxb32,Yx3CYd) 
 mq2TpYdadjYx3 = Matmul(mq2,TpYdadjYx3) 
 mu2YuadjYd = Matmul(mu2,YuadjYd) 
 Yb3ml2adjYe = Matmul(Yb3,ml2adjYe) 
 Yb3ml2adjYw3 = Matmul(Yb3,ml2adjYw3) 
 Yb3adjYeme2 = Matmul(Yb3,adjYeme2) 
 Yb3adjYw3mHw32 = Matmul(Yb3,adjYw3mHw32) 
 Yb3adjYw3MWM3 = Matmul(Yb3,adjYw3MWM3) 
! Yb3adjYw3BMWM3 = Matmul(Yb3,adjYw3BMWM3) 
 Ydmq2adjYu = Matmul(Yd,mq2adjYu) 
 YdadjYdTpYx3 = Matmul(Yd,adjYdTpYx3) 
! YdadjYdTpTYx3 = Matmul(Yd,adjYdTpTYx3) 
 YdadjYumu2 = Matmul(Yd,adjYumu2) 
 YdadjTYdadjYx3 = Matmul(Yd,adjTYdadjYx3) 
 Yeml2adjYb3 = Matmul(Ye,ml2adjYb3) 
 Yeml2adjYw3 = Matmul(Ye,ml2adjYw3) 
 YeadjYb3MBM3 = Matmul(Ye,adjYb3MBM3) 
 YeadjYb3mHb32 = Matmul(Ye,adjYb3mHb32) 
! YeadjYb3BMBM3 = Matmul(Ye,adjYb3BMBM3) 
 YeadjYw3mHw32 = Matmul(Ye,adjYw3mHw32) 
 YeadjYw3MWM3 = Matmul(Ye,adjYw3MWM3) 
! YeadjYw3BMWM3 = Matmul(Ye,adjYw3BMWM3) 
 Yumq2adjYd = Matmul(Yu,mq2adjYd) 
 YuadjYdmd2 = Matmul(Yu,adjYdmd2) 
 Yw3ml2adjYb3 = Matmul(Yw3,ml2adjYb3) 
 Yw3ml2adjYe = Matmul(Yw3,ml2adjYe) 
 Yw3adjYb3MBM3 = Matmul(Yw3,adjYb3MBM3) 
 Yw3adjYb3mHb32 = Matmul(Yw3,adjYb3mHb32) 
! Yw3adjYb3BMBM3 = Matmul(Yw3,adjYb3BMBM3) 
 Yw3adjYeme2 = Matmul(Yw3,adjYeme2) 
 Yx3md2CYd = Matmul(Yx3,md2CYd) 
 Yx3CYdmq2 = Matmul(Yx3,CYdmq2) 
 adjYb3Yb3adjYb3 = Matmul(adjYb3,Yb3adjYb3) 
 adjYb3Yb3adjYe = Matmul(adjYb3,Yb3adjYe) 
 adjYb3Yb3adjYw3 = Matmul(adjYb3,Yb3adjYw3) 
 adjYb3Yb3adjTYb3 = Matmul(adjYb3,Yb3adjTYb3) 
 adjYb3Yb3adjTYe = Matmul(adjYb3,Yb3adjTYe) 
 adjYb3Yb3adjTYw3 = Matmul(adjYb3,Yb3adjTYw3) 
 adjYb3TYb3adjYb3 = Matmul(adjYb3,TYb3adjYb3) 
 adjYb3TYb3adjYe = Matmul(adjYb3,TYb3adjYe) 
 adjYb3TYb3adjYw3 = Matmul(adjYb3,TYb3adjYw3) 
 adjYb3TYb3adjTYb3 = Matmul(adjYb3,TYb3adjTYb3) 
 adjYb3TYb3adjTYe = Matmul(adjYb3,TYb3adjTYe) 
 adjYb3TYb3adjTYw3 = Matmul(adjYb3,TYb3adjTYw3) 
 adjYdYdadjYd = Matmul(adjYd,YdadjYd) 
 adjYdYdadjYu = Matmul(adjYd,YdadjYu) 
 adjYdYdadjTYd = Matmul(adjYd,YdadjTYd) 
 adjYdYdadjTYu = Matmul(adjYd,YdadjTYu) 
 adjYdadjYx3Yx3 = Matmul(adjYd,adjYx3Yx3) 
 adjYdTYdadjYd = Matmul(adjYd,TYdadjYd) 
 adjYdTYdadjYu = Matmul(adjYd,TYdadjYu) 
 adjYdTYdadjTYd = Matmul(adjYd,TYdadjTYd) 
 adjYdTYdadjTYu = Matmul(adjYd,TYdadjTYu) 
 adjYdTpYx3CYx3 = Matmul(adjYd,TpYx3CYx3) 
 adjYdTpYx3CTYx3 = Matmul(adjYd,TpYx3CTYx3) 
 adjYdTpTYx3CTYx3 = Matmul(adjYd,TpTYx3CTYx3) 
 adjYeYeadjYb3 = Matmul(adjYe,YeadjYb3) 
 adjYeYeadjYe = Matmul(adjYe,YeadjYe) 
 adjYeYeadjYw3 = Matmul(adjYe,YeadjYw3) 
 adjYeYeadjTYb3 = Matmul(adjYe,YeadjTYb3) 
 adjYeYeadjTYe = Matmul(adjYe,YeadjTYe) 
 adjYeYeadjTYw3 = Matmul(adjYe,YeadjTYw3) 
 adjYeTYeadjYb3 = Matmul(adjYe,TYeadjYb3) 
 adjYeTYeadjYe = Matmul(adjYe,TYeadjYe) 
 adjYeTYeadjYw3 = Matmul(adjYe,TYeadjYw3) 
 adjYeTYeadjTYb3 = Matmul(adjYe,TYeadjTYb3) 
 adjYeTYeadjTYe = Matmul(adjYe,TYeadjTYe) 
 adjYeTYeadjTYw3 = Matmul(adjYe,TYeadjTYw3) 
 adjYuYuadjYd = Matmul(adjYu,YuadjYd) 
 adjYuYuadjYu = Matmul(adjYu,YuadjYu) 
 adjYuYuadjTYd = Matmul(adjYu,YuadjTYd) 
 adjYuYuadjTYu = Matmul(adjYu,YuadjTYu) 
 adjYuTYuadjYd = Matmul(adjYu,TYuadjYd) 
 adjYuTYuadjYu = Matmul(adjYu,TYuadjYu) 
 adjYuTYuadjTYd = Matmul(adjYu,TYuadjTYd) 
 adjYuTYuadjTYu = Matmul(adjYu,TYuadjTYu) 
 adjYw3Yw3adjYb3 = Matmul(adjYw3,Yw3adjYb3) 
 adjYw3Yw3adjYe = Matmul(adjYw3,Yw3adjYe) 
 adjYw3Yw3adjYw3 = Matmul(adjYw3,Yw3adjYw3) 
 adjYw3Yw3adjTYb3 = Matmul(adjYw3,Yw3adjTYb3) 
 adjYw3Yw3adjTYe = Matmul(adjYw3,Yw3adjTYe) 
 adjYw3Yw3adjTYw3 = Matmul(adjYw3,Yw3adjTYw3) 
 adjYw3TYw3adjYb3 = Matmul(adjYw3,TYw3adjYb3) 
 adjYw3TYw3adjYe = Matmul(adjYw3,TYw3adjYe) 
 adjYw3TYw3adjYw3 = Matmul(adjYw3,TYw3adjYw3) 
 adjYw3TYw3adjTYb3 = Matmul(adjYw3,TYw3adjTYb3) 
 adjYw3TYw3adjTYe = Matmul(adjYw3,TYw3adjTYe) 
 adjYw3TYw3adjTYw3 = Matmul(adjYw3,TYw3adjTYw3) 
 adjYx3Yx3adjYx3 = Matmul(adjYx3,Yx3adjYx3) 
 adjYx3Yx3adjTYx3 = Matmul(adjYx3,Yx3adjTYx3) 
 adjYx3Yx3CYd = Matmul(adjYx3,Yx3CYd) 
 adjYx3Yx3CTYd = Matmul(adjYx3,Yx3CTYd) 
 adjYx3TYx3adjYx3 = Matmul(adjYx3,TYx3adjYx3) 
 adjYx3TYx3adjTYx3 = Matmul(adjYx3,TYx3adjTYx3) 
 adjYx3TYx3CTYd = Matmul(adjYx3,TYx3CTYd) 
 adjTYb3TYb3adjYb3 = Matmul(adjTYb3,TYb3adjYb3) 
 adjTYb3TYb3adjYe = Matmul(adjTYb3,TYb3adjYe) 
 adjTYb3TYb3adjYw3 = Matmul(adjTYb3,TYb3adjYw3) 
 adjTYdTYdadjYd = Matmul(adjTYd,TYdadjYd) 
 adjTYdTYdadjYu = Matmul(adjTYd,TYdadjYu) 
 adjTYeTYeadjYb3 = Matmul(adjTYe,TYeadjYb3) 
 adjTYeTYeadjYe = Matmul(adjTYe,TYeadjYe) 
 adjTYeTYeadjYw3 = Matmul(adjTYe,TYeadjYw3) 
 adjTYuTYuadjYd = Matmul(adjTYu,TYuadjYd) 
 adjTYuTYuadjYu = Matmul(adjTYu,TYuadjYu) 
 adjTYw3TYw3adjYb3 = Matmul(adjTYw3,TYw3adjYb3) 
 adjTYw3TYw3adjYe = Matmul(adjTYw3,TYw3adjYe) 
 adjTYw3TYw3adjYw3 = Matmul(adjTYw3,TYw3adjYw3) 
 adjTYx3TYx3adjYx3 = Matmul(adjTYx3,TYx3adjYx3) 
 CYb3TpYb3CYb3 = Matmul(Conjg(Yb3),TpYb3CYb3) 
 CYb3TpYb3CTYb3 = Matmul(Conjg(Yb3),TpYb3CTYb3) 
 CYb3TpYw3CYw3 = Matmul(Conjg(Yb3),TpYw3CYw3) 
 CYb3TpYw3CTYw3 = Matmul(Conjg(Yb3),TpYw3CTYw3) 
 CYb3TpTYb3CTYb3 = Matmul(Conjg(Yb3),TpTYb3CTYb3) 
 CYb3TpTYw3CTYw3 = Matmul(Conjg(Yb3),TpTYw3CTYw3) 
 CYdTpYdadjYx3 = Matmul(Conjg(Yd),TpYdadjYx3) 
 CYdTpYdadjTYx3 = Matmul(Conjg(Yd),TpYdadjTYx3) 
 CYdTpYdCYd = Matmul(Conjg(Yd),TpYdCYd) 
 CYdTpYdCTYd = Matmul(Conjg(Yd),TpYdCTYd) 
 CYdTpTYdCTYd = Matmul(Conjg(Yd),TpTYdCTYd) 
 CYeTpYeCYe = Matmul(Conjg(Ye),TpYeCYe) 
 CYeTpYeCTYe = Matmul(Conjg(Ye),TpYeCTYe) 
 CYeTpTYeCTYe = Matmul(Conjg(Ye),TpTYeCTYe) 
 CYuTpYuCYu = Matmul(Conjg(Yu),TpYuCYu) 
 CYuTpYuCTYu = Matmul(Conjg(Yu),TpYuCTYu) 
 CYuTpTYuCTYu = Matmul(Conjg(Yu),TpTYuCTYu) 
 CYw3TpYb3CYb3 = Matmul(Conjg(Yw3),TpYb3CYb3) 
 CYw3TpYb3CTYb3 = Matmul(Conjg(Yw3),TpYb3CTYb3) 
 CYw3TpYw3CYw3 = Matmul(Conjg(Yw3),TpYw3CYw3) 
 CYw3TpYw3CTYw3 = Matmul(Conjg(Yw3),TpYw3CTYw3) 
 CYw3TpTYb3CTYb3 = Matmul(Conjg(Yw3),TpTYb3CTYb3) 
 CYw3TpTYw3CTYw3 = Matmul(Conjg(Yw3),TpTYw3CTYw3) 
 CYx3YdadjYd = Matmul(Conjg(Yx3),YdadjYd) 
 CYx3TpYx3CYx3 = Matmul(Conjg(Yx3),TpYx3CYx3) 
 CYx3TpYx3CTYx3 = Matmul(Conjg(Yx3),TpYx3CTYx3) 
 CYx3TpTYx3CTYx3 = Matmul(Conjg(Yx3),TpTYx3CTYx3) 
 CTYb3TpTYb3CYb3 = Matmul(Conjg(TYb3),TpTYb3CYb3) 
 CTYb3TpTYw3CYw3 = Matmul(Conjg(TYb3),TpTYw3CYw3) 
 CTYdTpTYdadjYx3 = Matmul(Conjg(TYd),TpTYdadjYx3) 
 CTYdTpTYdCYd = Matmul(Conjg(TYd),TpTYdCYd) 
 CTYeTpTYeCYe = Matmul(Conjg(TYe),TpTYeCYe) 
 CTYuTpTYuCYu = Matmul(Conjg(TYu),TpTYuCYu) 
 CTYw3TpTYb3CYb3 = Matmul(Conjg(TYw3),TpTYb3CYb3) 
 CTYw3TpTYw3CYw3 = Matmul(Conjg(TYw3),TpTYw3CYw3) 
 CTYx3TpTYx3CYx3 = Matmul(Conjg(TYx3),TpTYx3CYx3) 
! TYb3adjYw3MWM3 = Matmul(TYb3,adjYw3MWM3) 
 TYb3TpYb3CYb3 = Matmul(TYb3,TpYb3CYb3) 
 TYb3TpYw3CYw3 = Matmul(TYb3,TpYw3CYw3) 
 TYdadjYdadjTYx3 = Matmul(TYd,adjYdadjTYx3) 
! TYdadjYdTpYx3 = Matmul(TYd,adjYdTpYx3) 
 TYdTpYdCYd = Matmul(TYd,TpYdCYd) 
! TYeadjYb3MBM3 = Matmul(TYe,adjYb3MBM3) 
! TYeadjYw3MWM3 = Matmul(TYe,adjYw3MWM3) 
 TYeTpYeCYe = Matmul(TYe,TpYeCYe) 
 TYuTpYuCYu = Matmul(TYu,TpYuCYu) 
! TYw3adjYb3MBM3 = Matmul(TYw3,adjYb3MBM3) 
 TYw3TpYb3CYb3 = Matmul(TYw3,TpYb3CYb3) 
 TYw3TpYw3CYw3 = Matmul(TYw3,TpYw3CYw3) 
 TYx3CTYdTpYd = Matmul(TYx3,CTYdTpYd) 
 TYx3TpYx3CYx3 = Matmul(TYx3,TpYx3CYx3) 
 TpYb3CYb3TpYb3 = Matmul(Transpose(Yb3),CYb3TpYb3) 
 TpYb3CYb3TpYw3 = Matmul(Transpose(Yb3),CYb3TpYw3) 
! TpYb3CYb3TpTYb3 = Matmul(Transpose(Yb3),CYb3TpTYb3) 
! TpYb3CYb3TpTYw3 = Matmul(Transpose(Yb3),CYb3TpTYw3) 
 TpYb3CTYb3adjYb3 = Matmul(Transpose(Yb3),CTYb3adjYb3) 
 TpYb3CTYb3adjYe = Matmul(Transpose(Yb3),CTYb3adjYe) 
 TpYb3CTYb3adjYw3 = Matmul(Transpose(Yb3),CTYb3adjYw3) 
 TpYdmd2adjYx3 = Matmul(Transpose(Yd),md2adjYx3) 
 TpYdadjYx3mHxb32 = Matmul(Transpose(Yd),adjYx3mHxb32) 
 TpYdadjYx3Yx3 = Matmul(Transpose(Yd),adjYx3Yx3) 
 TpYdadjYx3TYx3 = Matmul(Transpose(Yd),adjYx3TYx3) 
 TpYdCYdTpYd = Matmul(Transpose(Yd),CYdTpYd) 
 TpYdCYdTpTYd = Matmul(Transpose(Yd),CYdTpTYd) 
 TpYdCTYdadjYd = Matmul(Transpose(Yd),CTYdadjYd) 
 TpYdCTYdadjYu = Matmul(Transpose(Yd),CTYdadjYu) 
 TpYeCYeTpYb3 = Matmul(Transpose(Ye),CYeTpYb3) 
 TpYeCYeTpYw3 = Matmul(Transpose(Ye),CYeTpYw3) 
! TpYeCYeTpTYb3 = Matmul(Transpose(Ye),CYeTpTYb3) 
! TpYeCYeTpTYw3 = Matmul(Transpose(Ye),CYeTpTYw3) 
 TpYeCTYeadjYb3 = Matmul(Transpose(Ye),CTYeadjYb3) 
 TpYeCTYeadjYe = Matmul(Transpose(Ye),CTYeadjYe) 
 TpYeCTYeadjYw3 = Matmul(Transpose(Ye),CTYeadjYw3) 
 TpYuCYuTpYd = Matmul(Transpose(Yu),CYuTpYd) 
 TpYuCYuTpTYd = Matmul(Transpose(Yu),CYuTpTYd) 
 TpYuCTYuadjYd = Matmul(Transpose(Yu),CTYuadjYd) 
 TpYuCTYuadjYu = Matmul(Transpose(Yu),CTYuadjYu) 
 TpYw3CYw3TpYb3 = Matmul(Transpose(Yw3),CYw3TpYb3) 
 TpYw3CYw3TpYw3 = Matmul(Transpose(Yw3),CYw3TpYw3) 
! TpYw3CYw3TpTYb3 = Matmul(Transpose(Yw3),CYw3TpTYb3) 
! TpYw3CYw3TpTYw3 = Matmul(Transpose(Yw3),CYw3TpTYw3) 
 TpYw3CTYw3adjYb3 = Matmul(Transpose(Yw3),CTYw3adjYb3) 
 TpYw3CTYw3adjYe = Matmul(Transpose(Yw3),CTYw3adjYe) 
 TpYw3CTYw3adjYw3 = Matmul(Transpose(Yw3),CTYw3adjYw3) 
 TpYx3CYx3TpYx3 = Matmul(Transpose(Yx3),CYx3TpYx3) 
! TpYx3CYx3TpTYx3 = Matmul(Transpose(Yx3),CYx3TpTYx3) 
 TpYx3CTYx3adjYx3 = Matmul(Transpose(Yx3),CTYx3adjYx3) 
 TpYx3CTYx3TYd = Matmul(Transpose(Yx3),CTYx3TYd) 
 TpYx3CTYx3TpTYd = Matmul(Transpose(Yx3),CTYx3TpTYd) 
! TpTYb3CYb3TpYb3 = Matmul(Transpose(TYb3),CYb3TpYb3) 
! TpTYb3CYb3TpYw3 = Matmul(Transpose(TYb3),CYb3TpYw3) 
 TpTYb3CTYb3adjYb3 = Matmul(Transpose(TYb3),CTYb3adjYb3) 
 TpTYb3CTYb3adjYe = Matmul(Transpose(TYb3),CTYb3adjYe) 
 TpTYb3CTYb3adjYw3 = Matmul(Transpose(TYb3),CTYb3adjYw3) 
 TpTYdCYdTpYd = Matmul(Transpose(TYd),CYdTpYd) 
 TpTYdCTYdadjYd = Matmul(Transpose(TYd),CTYdadjYd) 
 TpTYdCTYdadjYu = Matmul(Transpose(TYd),CTYdadjYu) 
! TpTYeCYeTpYb3 = Matmul(Transpose(TYe),CYeTpYb3) 
! TpTYeCYeTpYw3 = Matmul(Transpose(TYe),CYeTpYw3) 
 TpTYeCTYeadjYb3 = Matmul(Transpose(TYe),CTYeadjYb3) 
 TpTYeCTYeadjYe = Matmul(Transpose(TYe),CTYeadjYe) 
 TpTYeCTYeadjYw3 = Matmul(Transpose(TYe),CTYeadjYw3) 
 TpTYuCYuTpYd = Matmul(Transpose(TYu),CYuTpYd) 
 TpTYuCTYuadjYd = Matmul(Transpose(TYu),CTYuadjYd) 
 TpTYuCTYuadjYu = Matmul(Transpose(TYu),CTYuadjYu) 
! TpTYw3CYw3TpYb3 = Matmul(Transpose(TYw3),CYw3TpYb3) 
! TpTYw3CYw3TpYw3 = Matmul(Transpose(TYw3),CYw3TpYw3) 
 TpTYw3CTYw3adjYb3 = Matmul(Transpose(TYw3),CTYw3adjYb3) 
 TpTYw3CTYw3adjYe = Matmul(Transpose(TYw3),CTYw3adjYe) 
 TpTYw3CTYw3adjYw3 = Matmul(Transpose(TYw3),CTYw3adjYw3) 
! TpTYx3CYx3TpYx3 = Matmul(Transpose(TYx3),CYx3TpYx3) 
 TpTYx3CTYx3adjYx3 = Matmul(Transpose(TYx3),CTYx3adjYx3) 
 TpTYx3CTYx3TpYd = Matmul(Transpose(TYx3),CTYx3TpYd) 
 md2adjYx3Yx3adjYx3 = Matmul(md2,adjYx3Yx3adjYx3) 
 md2adjYx3Yx3CYd = Matmul(md2,adjYx3Yx3CYd) 
 md2CYdTpYdadjYx3 = Matmul(md2,CYdTpYdadjYx3) 
 md2CYdTpYdCYd = Matmul(md2,CYdTpYdCYd) 
 me2CYeTpYeCYe = Matmul(me2,CYeTpYeCYe) 
 mHb32CYb3TpYb3CYb3 = Matmul(mHb32,CYb3TpYb3CYb3) 
 mHb32CYb3TpYw3CYw3 = Matmul(mHb32,CYb3TpYw3CYw3) 
 mHw32CYw3TpYb3CYb3 = Matmul(mHw32,CYw3TpYb3CYb3) 
 mHw32CYw3TpYw3CYw3 = Matmul(mHw32,CYw3TpYw3CYw3) 
 mHxb32CYx3TpYx3CYx3 = Matmul(mHxb32,CYx3TpYx3CYx3) 
 ml2adjYb3Yb3adjYb3 = Matmul(ml2,adjYb3Yb3adjYb3) 
 ml2adjYb3Yb3adjYe = Matmul(ml2,adjYb3Yb3adjYe) 
 ml2adjYb3Yb3adjYw3 = Matmul(ml2,adjYb3Yb3adjYw3) 
 ml2adjYeYeadjYb3 = Matmul(ml2,adjYeYeadjYb3) 
 ml2adjYeYeadjYe = Matmul(ml2,adjYeYeadjYe) 
 ml2adjYeYeadjYw3 = Matmul(ml2,adjYeYeadjYw3) 
 ml2adjYw3Yw3adjYb3 = Matmul(ml2,adjYw3Yw3adjYb3) 
 ml2adjYw3Yw3adjYe = Matmul(ml2,adjYw3Yw3adjYe) 
 ml2adjYw3Yw3adjYw3 = Matmul(ml2,adjYw3Yw3adjYw3) 
 mq2adjYdYdadjYd = Matmul(mq2,adjYdYdadjYd) 
 mq2adjYdYdadjYu = Matmul(mq2,adjYdYdadjYu) 
 mq2adjYuYuadjYd = Matmul(mq2,adjYuYuadjYd) 
 mq2adjYuYuadjYu = Matmul(mq2,adjYuYuadjYu) 
 mq2adjYx3Yx3CYd = Matmul(mq2,adjYx3Yx3CYd) 
 mu2CYuTpYuCYu = Matmul(mu2,CYuTpYuCYu) 
 Yb3adjYb3Yb3adjYb3 = Matmul(Yb3,adjYb3Yb3adjYb3) 
Forall(i2=1:3)  Yb3adjYb3Yb3adjYb3(i2,i2) =  Real(Yb3adjYb3Yb3adjYb3(i2,i2),dp) 
 Yb3adjYb3TYb3adjYb3 = Matmul(Yb3,adjYb3TYb3adjYb3) 
 Yb3adjYb3TYb3adjTYb3 = Matmul(Yb3,adjYb3TYb3adjTYb3) 
 Yb3adjYeYeadjYb3 = Matmul(Yb3,adjYeYeadjYb3) 
Forall(i2=1:3)  Yb3adjYeYeadjYb3(i2,i2) =  Real(Yb3adjYeYeadjYb3(i2,i2),dp) 
 Yb3adjYeTYeadjYb3 = Matmul(Yb3,adjYeTYeadjYb3) 
 Yb3adjYeTYeadjTYb3 = Matmul(Yb3,adjYeTYeadjTYb3) 
 Yb3adjYw3Yw3adjYb3 = Matmul(Yb3,adjYw3Yw3adjYb3) 
Forall(i2=1:3)  Yb3adjYw3Yw3adjYb3(i2,i2) =  Real(Yb3adjYw3Yw3adjYb3(i2,i2),dp) 
 Yb3adjYw3TYw3adjYb3 = Matmul(Yb3,adjYw3TYw3adjYb3) 
 Yb3adjYw3TYw3adjTYb3 = Matmul(Yb3,adjYw3TYw3adjTYb3) 
 Yb3adjTYb3TYb3adjYb3 = Matmul(Yb3,adjTYb3TYb3adjYb3) 
 Yb3adjTYeTYeadjYb3 = Matmul(Yb3,adjTYeTYeadjYb3) 
 Yb3adjTYw3TYw3adjYb3 = Matmul(Yb3,adjTYw3TYw3adjYb3) 
 Yb3TpTYb3CTYb3adjYb3 = Matmul(Yb3,TpTYb3CTYb3adjYb3) 
 Yb3TpTYeCTYeadjYb3 = Matmul(Yb3,TpTYeCTYeadjYb3) 
 Yb3TpTYw3CTYw3adjYb3 = Matmul(Yb3,TpTYw3CTYw3adjYb3) 
 YdadjYdYdadjYd = Matmul(Yd,adjYdYdadjYd) 
Forall(i2=1:3)  YdadjYdYdadjYd(i2,i2) =  Real(YdadjYdYdadjYd(i2,i2),dp) 
 YdadjYdTYdadjYd = Matmul(Yd,adjYdTYdadjYd) 
 YdadjYdTYdadjTYd = Matmul(Yd,adjYdTYdadjTYd) 
 YdadjYdTpYx3CYx3 = Matmul(Yd,adjYdTpYx3CYx3) 
 YdadjYdTpTYx3CTYx3 = Matmul(Yd,adjYdTpTYx3CTYx3) 
 YdadjYuYuadjYd = Matmul(Yd,adjYuYuadjYd) 
Forall(i2=1:3)  YdadjYuYuadjYd(i2,i2) =  Real(YdadjYuYuadjYd(i2,i2),dp) 
 YdadjYuTYuadjYd = Matmul(Yd,adjYuTYuadjYd) 
 YdadjYuTYuadjTYd = Matmul(Yd,adjYuTYuadjTYd) 
 YdadjTYdTYdadjYd = Matmul(Yd,adjTYdTYdadjYd) 
 YdadjTYuTYuadjYd = Matmul(Yd,adjTYuTYuadjYd) 
 YdTpTYdCTYdadjYd = Matmul(Yd,TpTYdCTYdadjYd) 
 YdTpTYuCTYuadjYd = Matmul(Yd,TpTYuCTYuadjYd) 
 YeadjYb3Yb3adjYe = Matmul(Ye,adjYb3Yb3adjYe) 
Forall(i2=1:3)  YeadjYb3Yb3adjYe(i2,i2) =  Real(YeadjYb3Yb3adjYe(i2,i2),dp) 
 YeadjYb3TYb3adjYe = Matmul(Ye,adjYb3TYb3adjYe) 
 YeadjYb3TYb3adjTYe = Matmul(Ye,adjYb3TYb3adjTYe) 
 YeadjYeYeadjYe = Matmul(Ye,adjYeYeadjYe) 
Forall(i2=1:3)  YeadjYeYeadjYe(i2,i2) =  Real(YeadjYeYeadjYe(i2,i2),dp) 
 YeadjYeTYeadjYe = Matmul(Ye,adjYeTYeadjYe) 
 YeadjYeTYeadjTYe = Matmul(Ye,adjYeTYeadjTYe) 
 YeadjYw3Yw3adjYe = Matmul(Ye,adjYw3Yw3adjYe) 
Forall(i2=1:3)  YeadjYw3Yw3adjYe(i2,i2) =  Real(YeadjYw3Yw3adjYe(i2,i2),dp) 
 YeadjYw3TYw3adjYe = Matmul(Ye,adjYw3TYw3adjYe) 
 YeadjYw3TYw3adjTYe = Matmul(Ye,adjYw3TYw3adjTYe) 
 YeadjTYb3TYb3adjYe = Matmul(Ye,adjTYb3TYb3adjYe) 
 YeadjTYeTYeadjYe = Matmul(Ye,adjTYeTYeadjYe) 
 YeadjTYw3TYw3adjYe = Matmul(Ye,adjTYw3TYw3adjYe) 
 YeTpTYb3CTYb3adjYe = Matmul(Ye,TpTYb3CTYb3adjYe) 
 YeTpTYeCTYeadjYe = Matmul(Ye,TpTYeCTYeadjYe) 
 YeTpTYw3CTYw3adjYe = Matmul(Ye,TpTYw3CTYw3adjYe) 
 YuadjYdYdadjYu = Matmul(Yu,adjYdYdadjYu) 
Forall(i2=1:3)  YuadjYdYdadjYu(i2,i2) =  Real(YuadjYdYdadjYu(i2,i2),dp) 
 YuadjYdTYdadjYu = Matmul(Yu,adjYdTYdadjYu) 
 YuadjYdTYdadjTYu = Matmul(Yu,adjYdTYdadjTYu) 
 YuadjYuYuadjYu = Matmul(Yu,adjYuYuadjYu) 
Forall(i2=1:3)  YuadjYuYuadjYu(i2,i2) =  Real(YuadjYuYuadjYu(i2,i2),dp) 
 YuadjYuTYuadjYu = Matmul(Yu,adjYuTYuadjYu) 
 YuadjYuTYuadjTYu = Matmul(Yu,adjYuTYuadjTYu) 
 YuadjTYdTYdadjYu = Matmul(Yu,adjTYdTYdadjYu) 
 YuadjTYuTYuadjYu = Matmul(Yu,adjTYuTYuadjYu) 
 YuTpTYdCTYdadjYu = Matmul(Yu,TpTYdCTYdadjYu) 
 YuTpTYuCTYuadjYu = Matmul(Yu,TpTYuCTYuadjYu) 
 Yw3adjYb3Yb3adjYw3 = Matmul(Yw3,adjYb3Yb3adjYw3) 
Forall(i2=1:3)  Yw3adjYb3Yb3adjYw3(i2,i2) =  Real(Yw3adjYb3Yb3adjYw3(i2,i2),dp) 
 Yw3adjYb3TYb3adjYw3 = Matmul(Yw3,adjYb3TYb3adjYw3) 
 Yw3adjYb3TYb3adjTYw3 = Matmul(Yw3,adjYb3TYb3adjTYw3) 
 Yw3adjYeYeadjYw3 = Matmul(Yw3,adjYeYeadjYw3) 
Forall(i2=1:3)  Yw3adjYeYeadjYw3(i2,i2) =  Real(Yw3adjYeYeadjYw3(i2,i2),dp) 
 Yw3adjYeTYeadjYw3 = Matmul(Yw3,adjYeTYeadjYw3) 
 Yw3adjYeTYeadjTYw3 = Matmul(Yw3,adjYeTYeadjTYw3) 
 Yw3adjYw3Yw3adjYw3 = Matmul(Yw3,adjYw3Yw3adjYw3) 
Forall(i2=1:3)  Yw3adjYw3Yw3adjYw3(i2,i2) =  Real(Yw3adjYw3Yw3adjYw3(i2,i2),dp) 
 Yw3adjYw3TYw3adjYw3 = Matmul(Yw3,adjYw3TYw3adjYw3) 
 Yw3adjYw3TYw3adjTYw3 = Matmul(Yw3,adjYw3TYw3adjTYw3) 
 Yw3adjTYb3TYb3adjYw3 = Matmul(Yw3,adjTYb3TYb3adjYw3) 
 Yw3adjTYeTYeadjYw3 = Matmul(Yw3,adjTYeTYeadjYw3) 
 Yw3adjTYw3TYw3adjYw3 = Matmul(Yw3,adjTYw3TYw3adjYw3) 
 Yw3TpTYb3CTYb3adjYw3 = Matmul(Yw3,TpTYb3CTYb3adjYw3) 
 Yw3TpTYeCTYeadjYw3 = Matmul(Yw3,TpTYeCTYeadjYw3) 
 Yw3TpTYw3CTYw3adjYw3 = Matmul(Yw3,TpTYw3CTYw3adjYw3) 
 Yx3adjYx3Yx3adjYx3 = Matmul(Yx3,adjYx3Yx3adjYx3) 
Forall(i2=1:3)  Yx3adjYx3Yx3adjYx3(i2,i2) =  Real(Yx3adjYx3Yx3adjYx3(i2,i2),dp) 
 Yx3adjYx3TYx3adjYx3 = Matmul(Yx3,adjYx3TYx3adjYx3) 
 Yx3adjYx3TYx3adjTYx3 = Matmul(Yx3,adjYx3TYx3adjTYx3) 
 Yx3adjTYx3TYx3adjYx3 = Matmul(Yx3,adjTYx3TYx3adjYx3) 
 Yx3CYdTpYdadjYx3 = Matmul(Yx3,CYdTpYdadjYx3) 
Forall(i2=1:3)  Yx3CYdTpYdadjYx3(i2,i2) =  Real(Yx3CYdTpYdadjYx3(i2,i2),dp) 
 Yx3CTYdTpTYdadjYx3 = Matmul(Yx3,CTYdTpTYdadjYx3) 
 Yx3TYdadjYdadjTYx3 = Matmul(Yx3,TYdadjYdadjTYx3) 
 Yx3TpTYx3CTYx3adjYx3 = Matmul(Yx3,TpTYx3CTYx3adjYx3) 
 adjYb3mHb32Yb3adjYb3 = Matmul(adjYb3,mHb32Yb3adjYb3) 
 adjYb3mHb32Yb3adjYe = Matmul(adjYb3,mHb32Yb3adjYe) 
 adjYb3mHb32Yb3adjYw3 = Matmul(adjYb3,mHb32Yb3adjYw3) 
 adjYb3Yb3ml2adjYb3 = Matmul(adjYb3,Yb3ml2adjYb3) 
 adjYb3Yb3ml2adjYe = Matmul(adjYb3,Yb3ml2adjYe) 
 adjYb3Yb3ml2adjYw3 = Matmul(adjYb3,Yb3ml2adjYw3) 
 adjYb3Yb3adjYb3MBM3 = Matmul(adjYb3,Yb3adjYb3MBM3) 
 adjYb3Yb3adjYb3mHb32 = Matmul(adjYb3,Yb3adjYb3mHb32) 
 adjYb3Yb3adjYb3Yb3 = Matmul(adjYb3,Yb3adjYb3Yb3) 
Forall(i2=1:3)  adjYb3Yb3adjYb3Yb3(i2,i2) =  Real(adjYb3Yb3adjYb3Yb3(i2,i2),dp) 
! adjYb3Yb3adjYb3BMBM3 = Matmul(adjYb3,Yb3adjYb3BMBM3) 
 adjYb3Yb3adjYb3TYb3 = Matmul(adjYb3,Yb3adjYb3TYb3) 
 adjYb3Yb3adjYeme2 = Matmul(adjYb3,Yb3adjYeme2) 
 adjYb3Yb3adjYeYe = Matmul(adjYb3,Yb3adjYeYe) 
 adjYb3Yb3adjYeTYe = Matmul(adjYb3,Yb3adjYeTYe) 
 adjYb3Yb3adjYw3mHw32 = Matmul(adjYb3,Yb3adjYw3mHw32) 
 adjYb3Yb3adjYw3MWM3 = Matmul(adjYb3,Yb3adjYw3MWM3) 
 adjYb3Yb3adjYw3Yw3 = Matmul(adjYb3,Yb3adjYw3Yw3) 
! adjYb3Yb3adjYw3BMWM3 = Matmul(adjYb3,Yb3adjYw3BMWM3) 
 adjYb3Yb3adjYw3TYw3 = Matmul(adjYb3,Yb3adjYw3TYw3) 
! adjYb3TYb3adjYb3MBM3 = Matmul(adjYb3,TYb3adjYb3MBM3) 
 adjYb3TYb3adjYb3Yb3 = Matmul(adjYb3,TYb3adjYb3Yb3) 
 adjYb3TYb3adjYeYe = Matmul(adjYb3,TYb3adjYeYe) 
! adjYb3TYb3adjYw3MWM3 = Matmul(adjYb3,TYb3adjYw3MWM3) 
 adjYb3TYb3adjYw3Yw3 = Matmul(adjYb3,TYb3adjYw3Yw3) 
 adjYdmd2YdadjYd = Matmul(adjYd,md2YdadjYd) 
 adjYdmd2YdadjYu = Matmul(adjYd,md2YdadjYu) 
 adjYdYdmq2adjYd = Matmul(adjYd,Ydmq2adjYd) 
 adjYdYdmq2adjYu = Matmul(adjYd,Ydmq2adjYu) 
 adjYdYdadjYdmd2 = Matmul(adjYd,YdadjYdmd2) 
 adjYdYdadjYdYd = Matmul(adjYd,YdadjYdYd) 
Forall(i2=1:3)  adjYdYdadjYdYd(i2,i2) =  Real(adjYdYdadjYdYd(i2,i2),dp) 
 adjYdYdadjYdTYd = Matmul(adjYd,YdadjYdTYd) 
 adjYdYdadjYumu2 = Matmul(adjYd,YdadjYumu2) 
 adjYdYdadjYuYu = Matmul(adjYd,YdadjYuYu) 
 adjYdYdadjYuTYu = Matmul(adjYd,YdadjYuTYu) 
 adjYdTYdadjYdYd = Matmul(adjYd,TYdadjYdYd) 
 adjYdTYdadjYuYu = Matmul(adjYd,TYdadjYuYu) 
 adjYdTpYx3CYx3Yd = Matmul(adjYd,TpYx3CYx3Yd) 
Forall(i2=1:3)  adjYdTpYx3CYx3Yd(i2,i2) =  Real(adjYdTpYx3CYx3Yd(i2,i2),dp) 
 adjYdTpYx3CYx3TYd = Matmul(adjYd,TpYx3CYx3TYd) 
 adjYdTpYx3CTYx3TYd = Matmul(adjYd,TpYx3CTYx3TYd) 
 adjYdTpYx3CTYx3TpTYd = Matmul(adjYd,TpYx3CTYx3TpTYd) 
 adjYdTpTYx3CYx3Yd = Matmul(adjYd,TpTYx3CYx3Yd) 
 adjYdTpTYx3CTYx3TpYd = Matmul(adjYd,TpTYx3CTYx3TpYd) 
 adjYeme2YeadjYb3 = Matmul(adjYe,me2YeadjYb3) 
 adjYeme2YeadjYe = Matmul(adjYe,me2YeadjYe) 
 adjYeme2YeadjYw3 = Matmul(adjYe,me2YeadjYw3) 
 adjYeYeml2adjYb3 = Matmul(adjYe,Yeml2adjYb3) 
 adjYeYeml2adjYe = Matmul(adjYe,Yeml2adjYe) 
 adjYeYeml2adjYw3 = Matmul(adjYe,Yeml2adjYw3) 
 adjYeYeadjYb3MBM3 = Matmul(adjYe,YeadjYb3MBM3) 
 adjYeYeadjYb3mHb32 = Matmul(adjYe,YeadjYb3mHb32) 
 adjYeYeadjYb3Yb3 = Matmul(adjYe,YeadjYb3Yb3) 
! adjYeYeadjYb3BMBM3 = Matmul(adjYe,YeadjYb3BMBM3) 
 adjYeYeadjYb3TYb3 = Matmul(adjYe,YeadjYb3TYb3) 
 adjYeYeadjYeme2 = Matmul(adjYe,YeadjYeme2) 
 adjYeYeadjYeYe = Matmul(adjYe,YeadjYeYe) 
Forall(i2=1:3)  adjYeYeadjYeYe(i2,i2) =  Real(adjYeYeadjYeYe(i2,i2),dp) 
 adjYeYeadjYeTYe = Matmul(adjYe,YeadjYeTYe) 
 adjYeYeadjYw3mHw32 = Matmul(adjYe,YeadjYw3mHw32) 
 adjYeYeadjYw3MWM3 = Matmul(adjYe,YeadjYw3MWM3) 
 adjYeYeadjYw3Yw3 = Matmul(adjYe,YeadjYw3Yw3) 
! adjYeYeadjYw3BMWM3 = Matmul(adjYe,YeadjYw3BMWM3) 
 adjYeYeadjYw3TYw3 = Matmul(adjYe,YeadjYw3TYw3) 
! adjYeTYeadjYb3MBM3 = Matmul(adjYe,TYeadjYb3MBM3) 
 adjYeTYeadjYb3Yb3 = Matmul(adjYe,TYeadjYb3Yb3) 
 adjYeTYeadjYeYe = Matmul(adjYe,TYeadjYeYe) 
! adjYeTYeadjYw3MWM3 = Matmul(adjYe,TYeadjYw3MWM3) 
 adjYeTYeadjYw3Yw3 = Matmul(adjYe,TYeadjYw3Yw3) 
 adjYumu2YuadjYd = Matmul(adjYu,mu2YuadjYd) 
 adjYumu2YuadjYu = Matmul(adjYu,mu2YuadjYu) 
 adjYuYumq2adjYd = Matmul(adjYu,Yumq2adjYd) 
 adjYuYumq2adjYu = Matmul(adjYu,Yumq2adjYu) 
 adjYuYuadjYdmd2 = Matmul(adjYu,YuadjYdmd2) 
 adjYuYuadjYdYd = Matmul(adjYu,YuadjYdYd) 
 adjYuYuadjYdTYd = Matmul(adjYu,YuadjYdTYd) 
 adjYuYuadjYumu2 = Matmul(adjYu,YuadjYumu2) 
 adjYuYuadjYuYu = Matmul(adjYu,YuadjYuYu) 
Forall(i2=1:3)  adjYuYuadjYuYu(i2,i2) =  Real(adjYuYuadjYuYu(i2,i2),dp) 
 adjYuYuadjYuTYu = Matmul(adjYu,YuadjYuTYu) 
 adjYuTYuadjYdYd = Matmul(adjYu,TYuadjYdYd) 
 adjYuTYuadjYuYu = Matmul(adjYu,TYuadjYuYu) 
 adjYw3mHw32Yw3adjYb3 = Matmul(adjYw3,mHw32Yw3adjYb3) 
 adjYw3mHw32Yw3adjYe = Matmul(adjYw3,mHw32Yw3adjYe) 
 adjYw3mHw32Yw3adjYw3 = Matmul(adjYw3,mHw32Yw3adjYw3) 
 adjYw3Yw3ml2adjYb3 = Matmul(adjYw3,Yw3ml2adjYb3) 
 adjYw3Yw3ml2adjYe = Matmul(adjYw3,Yw3ml2adjYe) 
 adjYw3Yw3ml2adjYw3 = Matmul(adjYw3,Yw3ml2adjYw3) 
 adjYw3Yw3adjYb3MBM3 = Matmul(adjYw3,Yw3adjYb3MBM3) 
 adjYw3Yw3adjYb3mHb32 = Matmul(adjYw3,Yw3adjYb3mHb32) 
 adjYw3Yw3adjYb3Yb3 = Matmul(adjYw3,Yw3adjYb3Yb3) 
! adjYw3Yw3adjYb3BMBM3 = Matmul(adjYw3,Yw3adjYb3BMBM3) 
 adjYw3Yw3adjYb3TYb3 = Matmul(adjYw3,Yw3adjYb3TYb3) 
 adjYw3Yw3adjYeme2 = Matmul(adjYw3,Yw3adjYeme2) 
 adjYw3Yw3adjYeYe = Matmul(adjYw3,Yw3adjYeYe) 
 adjYw3Yw3adjYeTYe = Matmul(adjYw3,Yw3adjYeTYe) 
 adjYw3Yw3adjYw3mHw32 = Matmul(adjYw3,Yw3adjYw3mHw32) 
 adjYw3Yw3adjYw3MWM3 = Matmul(adjYw3,Yw3adjYw3MWM3) 
 adjYw3Yw3adjYw3Yw3 = Matmul(adjYw3,Yw3adjYw3Yw3) 
Forall(i2=1:3)  adjYw3Yw3adjYw3Yw3(i2,i2) =  Real(adjYw3Yw3adjYw3Yw3(i2,i2),dp) 
! adjYw3Yw3adjYw3BMWM3 = Matmul(adjYw3,Yw3adjYw3BMWM3) 
 adjYw3Yw3adjYw3TYw3 = Matmul(adjYw3,Yw3adjYw3TYw3) 
! adjYw3TYw3adjYb3MBM3 = Matmul(adjYw3,TYw3adjYb3MBM3) 
 adjYw3TYw3adjYb3Yb3 = Matmul(adjYw3,TYw3adjYb3Yb3) 
 adjYw3TYw3adjYeYe = Matmul(adjYw3,TYw3adjYeYe) 
! adjYw3TYw3adjYw3MWM3 = Matmul(adjYw3,TYw3adjYw3MWM3) 
 adjYw3TYw3adjYw3Yw3 = Matmul(adjYw3,TYw3adjYw3Yw3) 
 adjYx3mHxb32Yx3adjYx3 = Matmul(adjYx3,mHxb32Yx3adjYx3) 
 adjYx3mHxb32Yx3CYd = Matmul(adjYx3,mHxb32Yx3CYd) 
 adjYx3Yx3md2adjYx3 = Matmul(adjYx3,Yx3md2adjYx3) 
 adjYx3Yx3md2CYd = Matmul(adjYx3,Yx3md2CYd) 
 adjYx3Yx3adjYx3mHxb32 = Matmul(adjYx3,Yx3adjYx3mHxb32) 
 adjYx3Yx3adjYx3Yx3 = Matmul(adjYx3,Yx3adjYx3Yx3) 
Forall(i2=1:3)  adjYx3Yx3adjYx3Yx3(i2,i2) =  Real(adjYx3Yx3adjYx3Yx3(i2,i2),dp) 
 adjYx3Yx3adjYx3TYx3 = Matmul(adjYx3,Yx3adjYx3TYx3) 
 adjYx3Yx3CYdmq2 = Matmul(adjYx3,Yx3CYdmq2) 
 adjYx3TYx3adjYx3Yx3 = Matmul(adjYx3,TYx3adjYx3Yx3) 
 adjYx3TYx3CYdTpYd = Matmul(adjYx3,TYx3CYdTpYd) 
 adjYx3TYx3CTYdTpYd = Matmul(adjYx3,TYx3CTYdTpYd) 
 adjTYb3TYb3TpYb3CYb3 = Matmul(adjTYb3,TYb3TpYb3CYb3) 
 adjTYb3TYb3TpYw3CYw3 = Matmul(adjTYb3,TYb3TpYw3CYw3) 
 adjTYdTYdTpYdCYd = Matmul(adjTYd,TYdTpYdCYd) 
 adjTYeTYeTpYeCYe = Matmul(adjTYe,TYeTpYeCYe) 
 adjTYuTYuTpYuCYu = Matmul(adjTYu,TYuTpYuCYu) 
 adjTYw3TYw3TpYb3CYb3 = Matmul(adjTYw3,TYw3TpYb3CYb3) 
 adjTYw3TYw3TpYw3CYw3 = Matmul(adjTYw3,TYw3TpYw3CYw3) 
 adjTYx3TYx3TpYx3CYx3 = Matmul(adjTYx3,TYx3TpYx3CYx3) 
 CYb3ml2TpYb3CYb3 = Matmul(Conjg(Yb3),ml2TpYb3CYb3) 
 CYb3ml2TpYw3CYw3 = Matmul(Conjg(Yb3),ml2TpYw3CYw3) 
 CYb3TpYb3mHb32CYb3 = Matmul(Conjg(Yb3),TpYb3mHb32CYb3) 
 CYb3TpYb3CYb3ml2 = Matmul(Conjg(Yb3),TpYb3CYb3ml2) 
 CYb3TpYb3CYb3TpYb3 = Matmul(Conjg(Yb3),TpYb3CYb3TpYb3) 
Forall(i2=1:3)  CYb3TpYb3CYb3TpYb3(i2,i2) =  Real(CYb3TpYb3CYb3TpYb3(i2,i2),dp) 
! CYb3TpYb3CYb3TpTYb3 = Matmul(Conjg(Yb3),TpYb3CYb3TpTYb3) 
 CYb3TpYeCYeTpYb3 = Matmul(Conjg(Yb3),TpYeCYeTpYb3) 
Forall(i2=1:3)  CYb3TpYeCYeTpYb3(i2,i2) =  Real(CYb3TpYeCYeTpYb3(i2,i2),dp) 
! CYb3TpYeCYeTpTYb3 = Matmul(Conjg(Yb3),TpYeCYeTpTYb3) 
 CYb3TpYw3mHw32CYw3 = Matmul(Conjg(Yb3),TpYw3mHw32CYw3) 
 CYb3TpYw3CYw3ml2 = Matmul(Conjg(Yb3),TpYw3CYw3ml2) 
 CYb3TpYw3CYw3TpYb3 = Matmul(Conjg(Yb3),TpYw3CYw3TpYb3) 
Forall(i2=1:3)  CYb3TpYw3CYw3TpYb3(i2,i2) =  Real(CYb3TpYw3CYw3TpYb3(i2,i2),dp) 
! CYb3TpYw3CYw3TpTYb3 = Matmul(Conjg(Yb3),TpYw3CYw3TpTYb3) 
! CYb3TpTYb3CYb3TpYb3 = Matmul(Conjg(Yb3),TpTYb3CYb3TpYb3) 
! CYb3TpTYeCYeTpYb3 = Matmul(Conjg(Yb3),TpTYeCYeTpYb3) 
! CYb3TpTYw3CYw3TpYb3 = Matmul(Conjg(Yb3),TpTYw3CYw3TpYb3) 
 CYdmq2TpYdadjYx3 = Matmul(Conjg(Yd),mq2TpYdadjYx3) 
 CYdmq2TpYdCYd = Matmul(Conjg(Yd),mq2TpYdCYd) 
 CYdTpYdmd2adjYx3 = Matmul(Conjg(Yd),TpYdmd2adjYx3) 
 CYdTpYdmd2CYd = Matmul(Conjg(Yd),TpYdmd2CYd) 
 CYdTpYdadjYx3mHxb32 = Matmul(Conjg(Yd),TpYdadjYx3mHxb32) 
 CYdTpYdadjYx3Yx3 = Matmul(Conjg(Yd),TpYdadjYx3Yx3) 
 CYdTpYdadjYx3TYx3 = Matmul(Conjg(Yd),TpYdadjYx3TYx3) 
 CYdTpYdCYdmq2 = Matmul(Conjg(Yd),TpYdCYdmq2) 
 CYdTpYdCYdTpYd = Matmul(Conjg(Yd),TpYdCYdTpYd) 
Forall(i2=1:3)  CYdTpYdCYdTpYd(i2,i2) =  Real(CYdTpYdCYdTpYd(i2,i2),dp) 
 CYdTpYdCYdTpTYd = Matmul(Conjg(Yd),TpYdCYdTpTYd) 
 CYdTpYuCYuTpYd = Matmul(Conjg(Yd),TpYuCYuTpYd) 
Forall(i2=1:3)  CYdTpYuCYuTpYd(i2,i2) =  Real(CYdTpYuCYuTpYd(i2,i2),dp) 
 CYdTpYuCYuTpTYd = Matmul(Conjg(Yd),TpYuCYuTpTYd) 
 CYdTpTYdCYdTpYd = Matmul(Conjg(Yd),TpTYdCYdTpYd) 
 CYdTpTYuCYuTpYd = Matmul(Conjg(Yd),TpTYuCYuTpYd) 
 CYeml2TpYeCYe = Matmul(Conjg(Ye),ml2TpYeCYe) 
 CYeTpYeme2CYe = Matmul(Conjg(Ye),TpYeme2CYe) 
 CYeTpYeCYeml2 = Matmul(Conjg(Ye),TpYeCYeml2) 
 CYumq2TpYuCYu = Matmul(Conjg(Yu),mq2TpYuCYu) 
 CYuTpYumu2CYu = Matmul(Conjg(Yu),TpYumu2CYu) 
 CYuTpYuCYumq2 = Matmul(Conjg(Yu),TpYuCYumq2) 
 CYw3ml2TpYb3CYb3 = Matmul(Conjg(Yw3),ml2TpYb3CYb3) 
 CYw3ml2TpYw3CYw3 = Matmul(Conjg(Yw3),ml2TpYw3CYw3) 
 CYw3TpYb3mHb32CYb3 = Matmul(Conjg(Yw3),TpYb3mHb32CYb3) 
 CYw3TpYb3CYb3ml2 = Matmul(Conjg(Yw3),TpYb3CYb3ml2) 
 CYw3TpYb3CYb3TpYw3 = Matmul(Conjg(Yw3),TpYb3CYb3TpYw3) 
Forall(i2=1:3)  CYw3TpYb3CYb3TpYw3(i2,i2) =  Real(CYw3TpYb3CYb3TpYw3(i2,i2),dp) 
! CYw3TpYb3CYb3TpTYw3 = Matmul(Conjg(Yw3),TpYb3CYb3TpTYw3) 
 CYw3TpYeCYeTpYw3 = Matmul(Conjg(Yw3),TpYeCYeTpYw3) 
Forall(i2=1:3)  CYw3TpYeCYeTpYw3(i2,i2) =  Real(CYw3TpYeCYeTpYw3(i2,i2),dp) 
! CYw3TpYeCYeTpTYw3 = Matmul(Conjg(Yw3),TpYeCYeTpTYw3) 
 CYw3TpYw3mHw32CYw3 = Matmul(Conjg(Yw3),TpYw3mHw32CYw3) 
 CYw3TpYw3CYw3ml2 = Matmul(Conjg(Yw3),TpYw3CYw3ml2) 
 CYw3TpYw3CYw3TpYw3 = Matmul(Conjg(Yw3),TpYw3CYw3TpYw3) 
Forall(i2=1:3)  CYw3TpYw3CYw3TpYw3(i2,i2) =  Real(CYw3TpYw3CYw3TpYw3(i2,i2),dp) 
! CYw3TpYw3CYw3TpTYw3 = Matmul(Conjg(Yw3),TpYw3CYw3TpTYw3) 
! CYw3TpTYb3CYb3TpYw3 = Matmul(Conjg(Yw3),TpTYb3CYb3TpYw3) 
! CYw3TpTYeCYeTpYw3 = Matmul(Conjg(Yw3),TpTYeCYeTpYw3) 
! CYw3TpTYw3CYw3TpYw3 = Matmul(Conjg(Yw3),TpTYw3CYw3TpYw3) 
 CYx3md2TpYx3CYx3 = Matmul(Conjg(Yx3),md2TpYx3CYx3) 
 CYx3YdadjYdTpYx3 = Matmul(Conjg(Yx3),YdadjYdTpYx3) 
Forall(i2=1:3)  CYx3YdadjYdTpYx3(i2,i2) =  Real(CYx3YdadjYdTpYx3(i2,i2),dp) 
! CYx3YdadjYdTpTYx3 = Matmul(Conjg(Yx3),YdadjYdTpTYx3) 
! CYx3TYdadjYdTpYx3 = Matmul(Conjg(Yx3),TYdadjYdTpYx3) 
 CYx3TpYx3mHxb32CYx3 = Matmul(Conjg(Yx3),TpYx3mHxb32CYx3) 
 CYx3TpYx3CYx3md2 = Matmul(Conjg(Yx3),TpYx3CYx3md2) 
 CYx3TpYx3CYx3Yd = Matmul(Conjg(Yx3),TpYx3CYx3Yd) 
 CYx3TpYx3CYx3TYd = Matmul(Conjg(Yx3),TpYx3CYx3TYd) 
 CYx3TpYx3CYx3TpYx3 = Matmul(Conjg(Yx3),TpYx3CYx3TpYx3) 
Forall(i2=1:3)  CYx3TpYx3CYx3TpYx3(i2,i2) =  Real(CYx3TpYx3CYx3TpYx3(i2,i2),dp) 
! CYx3TpYx3CYx3TpTYx3 = Matmul(Conjg(Yx3),TpYx3CYx3TpTYx3) 
 CYx3TpTYx3CYx3Yd = Matmul(Conjg(Yx3),TpTYx3CYx3Yd) 
! CYx3TpTYx3CYx3TpYx3 = Matmul(Conjg(Yx3),TpTYx3CYx3TpYx3) 
 TYb3adjYb3Yb3adjTYb3 = Matmul(TYb3,adjYb3Yb3adjTYb3) 
 TYb3adjYeYeadjTYb3 = Matmul(TYb3,adjYeYeadjTYb3) 
 TYb3adjYw3Yw3adjTYb3 = Matmul(TYb3,adjYw3Yw3adjTYb3) 
 TYb3TpYb3CTYb3adjYb3 = Matmul(TYb3,TpYb3CTYb3adjYb3) 
 TYb3TpYeCTYeadjYb3 = Matmul(TYb3,TpYeCTYeadjYb3) 
 TYb3TpYw3CTYw3adjYb3 = Matmul(TYb3,TpYw3CTYw3adjYb3) 
 TYdadjYdYdadjTYd = Matmul(TYd,adjYdYdadjTYd) 
 TYdadjYdadjYx3Yx3 = Matmul(TYd,adjYdadjYx3Yx3) 
 TYdadjYuYuadjTYd = Matmul(TYd,adjYuYuadjTYd) 
 TYdTpYdCTYdadjYd = Matmul(TYd,TpYdCTYdadjYd) 
 TYdTpYuCTYuadjYd = Matmul(TYd,TpYuCTYuadjYd) 
 TYeadjYb3Yb3adjTYe = Matmul(TYe,adjYb3Yb3adjTYe) 
 TYeadjYeYeadjTYe = Matmul(TYe,adjYeYeadjTYe) 
 TYeadjYw3Yw3adjTYe = Matmul(TYe,adjYw3Yw3adjTYe) 
 TYeTpYb3CTYb3adjYe = Matmul(TYe,TpYb3CTYb3adjYe) 
 TYeTpYeCTYeadjYe = Matmul(TYe,TpYeCTYeadjYe) 
 TYeTpYw3CTYw3adjYe = Matmul(TYe,TpYw3CTYw3adjYe) 
 TYuadjYdYdadjTYu = Matmul(TYu,adjYdYdadjTYu) 
 TYuadjYuYuadjTYu = Matmul(TYu,adjYuYuadjTYu) 
 TYuTpYdCTYdadjYu = Matmul(TYu,TpYdCTYdadjYu) 
 TYuTpYuCTYuadjYu = Matmul(TYu,TpYuCTYuadjYu) 
 TYw3adjYb3Yb3adjTYw3 = Matmul(TYw3,adjYb3Yb3adjTYw3) 
 TYw3adjYeYeadjTYw3 = Matmul(TYw3,adjYeYeadjTYw3) 
 TYw3adjYw3Yw3adjTYw3 = Matmul(TYw3,adjYw3Yw3adjTYw3) 
 TYw3TpYb3CTYb3adjYw3 = Matmul(TYw3,TpYb3CTYb3adjYw3) 
 TYw3TpYeCTYeadjYw3 = Matmul(TYw3,TpYeCTYeadjYw3) 
 TYw3TpYw3CTYw3adjYw3 = Matmul(TYw3,TpYw3CTYw3adjYw3) 
 TYx3YdadjTYdadjYx3 = Matmul(TYx3,YdadjTYdadjYx3) 
 TYx3adjYx3Yx3adjTYx3 = Matmul(TYx3,adjYx3Yx3adjTYx3) 
 TYx3CYdTpYdadjTYx3 = Matmul(TYx3,CYdTpYdadjTYx3) 
 TYx3TpYx3CTYx3adjYx3 = Matmul(TYx3,TpYx3CTYx3adjYx3) 
 TpYb3CYb3TpYb3CYb3 = Matmul(Transpose(Yb3),CYb3TpYb3CYb3) 
Forall(i2=1:3)  TpYb3CYb3TpYb3CYb3(i2,i2) =  Real(TpYb3CYb3TpYb3CYb3(i2,i2),dp) 
 TpYb3CYb3TpYw3CYw3 = Matmul(Transpose(Yb3),CYb3TpYw3CYw3) 
 TpYb3CYb3TpTYb3CTYb3 = Matmul(Transpose(Yb3),CYb3TpTYb3CTYb3) 
 TpYb3CYb3TpTYw3CTYw3 = Matmul(Transpose(Yb3),CYb3TpTYw3CTYw3) 
 TpYb3CTYb3TpTYb3CYb3 = Matmul(Transpose(Yb3),CTYb3TpTYb3CYb3) 
 TpYb3CTYb3TpTYw3CYw3 = Matmul(Transpose(Yb3),CTYb3TpTYw3CYw3) 
 TpYdadjYdTpTYx3CTYx3 = Matmul(Transpose(Yd),adjYdTpTYx3CTYx3) 
 TpYdadjYx3Yx3CYd = Matmul(Transpose(Yd),adjYx3Yx3CYd) 
Forall(i2=1:3)  TpYdadjYx3Yx3CYd(i2,i2) =  Real(TpYdadjYx3Yx3CYd(i2,i2),dp) 
 TpYdadjYx3TYx3CTYd = Matmul(Transpose(Yd),adjYx3TYx3CTYd) 
 TpYdCYdTpYdCYd = Matmul(Transpose(Yd),CYdTpYdCYd) 
Forall(i2=1:3)  TpYdCYdTpYdCYd(i2,i2) =  Real(TpYdCYdTpYdCYd(i2,i2),dp) 
 TpYdCYdTpTYdCTYd = Matmul(Transpose(Yd),CYdTpTYdCTYd) 
 TpYdCTYdTpTYdCYd = Matmul(Transpose(Yd),CTYdTpTYdCYd) 
 TpYeCYeTpYeCYe = Matmul(Transpose(Ye),CYeTpYeCYe) 
Forall(i2=1:3)  TpYeCYeTpYeCYe(i2,i2) =  Real(TpYeCYeTpYeCYe(i2,i2),dp) 
 TpYeCYeTpTYeCTYe = Matmul(Transpose(Ye),CYeTpTYeCTYe) 
 TpYeCTYeTpTYeCYe = Matmul(Transpose(Ye),CTYeTpTYeCYe) 
 TpYuCYuTpYuCYu = Matmul(Transpose(Yu),CYuTpYuCYu) 
Forall(i2=1:3)  TpYuCYuTpYuCYu(i2,i2) =  Real(TpYuCYuTpYuCYu(i2,i2),dp) 
 TpYuCYuTpTYuCTYu = Matmul(Transpose(Yu),CYuTpTYuCTYu) 
 TpYuCTYuTpTYuCYu = Matmul(Transpose(Yu),CTYuTpTYuCYu) 
 TpYw3CYw3TpYb3CYb3 = Matmul(Transpose(Yw3),CYw3TpYb3CYb3) 
 TpYw3CYw3TpYw3CYw3 = Matmul(Transpose(Yw3),CYw3TpYw3CYw3) 
Forall(i2=1:3)  TpYw3CYw3TpYw3CYw3(i2,i2) =  Real(TpYw3CYw3TpYw3CYw3(i2,i2),dp) 
 TpYw3CYw3TpTYb3CTYb3 = Matmul(Transpose(Yw3),CYw3TpTYb3CTYb3) 
 TpYw3CYw3TpTYw3CTYw3 = Matmul(Transpose(Yw3),CYw3TpTYw3CTYw3) 
 TpYw3CTYw3TpTYb3CYb3 = Matmul(Transpose(Yw3),CTYw3TpTYb3CYb3) 
 TpYw3CTYw3TpTYw3CYw3 = Matmul(Transpose(Yw3),CTYw3TpTYw3CYw3) 
 TpYx3CYx3YdadjYd = Matmul(Transpose(Yx3),CYx3YdadjYd) 
 TpYx3CYx3TpYx3CYx3 = Matmul(Transpose(Yx3),CYx3TpYx3CYx3) 
Forall(i2=1:3)  TpYx3CYx3TpYx3CYx3(i2,i2) =  Real(TpYx3CYx3TpYx3CYx3(i2,i2),dp) 
 TpYx3CYx3TpTYx3CTYx3 = Matmul(Transpose(Yx3),CYx3TpTYx3CTYx3) 
 TpYx3CTYx3TpTYx3CYx3 = Matmul(Transpose(Yx3),CTYx3TpTYx3CYx3) 
 TpTYb3CYb3TpYb3CTYb3 = Matmul(Transpose(TYb3),CYb3TpYb3CTYb3) 
 TpTYb3CYb3TpYw3CTYw3 = Matmul(Transpose(TYb3),CYb3TpYw3CTYw3) 
 TpTYdadjYdTpYx3CTYx3 = Matmul(Transpose(TYd),adjYdTpYx3CTYx3) 
 TpTYdadjYx3Yx3CTYd = Matmul(Transpose(TYd),adjYx3Yx3CTYd) 
 TpTYdCYdTpYdCTYd = Matmul(Transpose(TYd),CYdTpYdCTYd) 
 TpTYeCYeTpYeCTYe = Matmul(Transpose(TYe),CYeTpYeCTYe) 
 TpTYuCYuTpYuCTYu = Matmul(Transpose(TYu),CYuTpYuCTYu) 
 TpTYw3CYw3TpYb3CTYb3 = Matmul(Transpose(TYw3),CYw3TpYb3CTYb3) 
 TpTYw3CYw3TpYw3CTYw3 = Matmul(Transpose(TYw3),CYw3TpYw3CTYw3) 
 TpTYx3CYx3TpYx3CTYx3 = Matmul(Transpose(TYx3),CYx3TpYx3CTYx3) 
 MBM3CYb3TpYb3CYb3TpYb3 = Matmul(MBM3,CYb3TpYb3CYb3TpYb3) 
! MBM3CYb3TpYb3CYb3TpTYb3 = Matmul(MBM3,CYb3TpYb3CYb3TpTYb3) 
 MBM3CYb3TpYeCYeTpYb3 = Matmul(MBM3,CYb3TpYeCYeTpYb3) 
! MBM3CYb3TpYeCYeTpTYb3 = Matmul(MBM3,CYb3TpYeCYeTpTYb3) 
 MBM3CYb3TpYw3CYw3TpYb3 = Matmul(MBM3,CYb3TpYw3CYw3TpYb3) 
! MBM3CYb3TpYw3CYw3TpTYb3 = Matmul(MBM3,CYb3TpYw3CYw3TpTYb3) 
! MBM3CYb3TpTYb3CYb3TpYb3 = Matmul(MBM3,CYb3TpTYb3CYb3TpYb3) 
! MBM3CYb3TpTYeCYeTpYb3 = Matmul(MBM3,CYb3TpTYeCYeTpYb3) 
! MBM3CYb3TpTYw3CYw3TpYb3 = Matmul(MBM3,CYb3TpTYw3CYw3TpYb3) 
 md2YdadjYdYdadjYd = Matmul(md2,YdadjYdYdadjYd) 
 md2YdadjYdTpYx3CYx3 = Matmul(md2,YdadjYdTpYx3CYx3) 
 md2YdadjYuYuadjYd = Matmul(md2,YdadjYuYuadjYd) 
 md2adjYx3Yx3adjYx3Yx3 = Matmul(md2,adjYx3Yx3adjYx3Yx3) 
 md2TpYx3CYx3YdadjYd = Matmul(md2,TpYx3CYx3YdadjYd) 
 md2TpYx3CYx3TpYx3CYx3 = Matmul(md2,TpYx3CYx3TpYx3CYx3) 
 me2YeadjYb3Yb3adjYe = Matmul(me2,YeadjYb3Yb3adjYe) 
 me2YeadjYeYeadjYe = Matmul(me2,YeadjYeYeadjYe) 
 me2YeadjYw3Yw3adjYe = Matmul(me2,YeadjYw3Yw3adjYe) 
 mHb32Yb3adjYb3Yb3adjYb3 = Matmul(mHb32,Yb3adjYb3Yb3adjYb3) 
 mHb32Yb3adjYeYeadjYb3 = Matmul(mHb32,Yb3adjYeYeadjYb3) 
 mHb32Yb3adjYw3Yw3adjYb3 = Matmul(mHb32,Yb3adjYw3Yw3adjYb3) 
 mHw32Yw3adjYb3Yb3adjYw3 = Matmul(mHw32,Yw3adjYb3Yb3adjYw3) 
 mHw32Yw3adjYeYeadjYw3 = Matmul(mHw32,Yw3adjYeYeadjYw3) 
 mHw32Yw3adjYw3Yw3adjYw3 = Matmul(mHw32,Yw3adjYw3Yw3adjYw3) 
 mHxb32Yx3adjYx3Yx3adjYx3 = Matmul(mHxb32,Yx3adjYx3Yx3adjYx3) 
 mHxb32Yx3CYdTpYdadjYx3 = Matmul(mHxb32,Yx3CYdTpYdadjYx3) 
 mHxb32CYx3YdadjYdTpYx3 = Matmul(mHxb32,CYx3YdadjYdTpYx3) 
 ml2adjYb3Yb3adjYb3Yb3 = Matmul(ml2,adjYb3Yb3adjYb3Yb3) 
 ml2adjYb3Yb3adjYeYe = Matmul(ml2,adjYb3Yb3adjYeYe) 
 ml2adjYb3Yb3adjYw3Yw3 = Matmul(ml2,adjYb3Yb3adjYw3Yw3) 
 ml2adjYeYeadjYb3Yb3 = Matmul(ml2,adjYeYeadjYb3Yb3) 
 ml2adjYeYeadjYeYe = Matmul(ml2,adjYeYeadjYeYe) 
 ml2adjYeYeadjYw3Yw3 = Matmul(ml2,adjYeYeadjYw3Yw3) 
 ml2adjYw3Yw3adjYb3Yb3 = Matmul(ml2,adjYw3Yw3adjYb3Yb3) 
 ml2adjYw3Yw3adjYeYe = Matmul(ml2,adjYw3Yw3adjYeYe) 
 ml2adjYw3Yw3adjYw3Yw3 = Matmul(ml2,adjYw3Yw3adjYw3Yw3) 
 ml2TpYb3CYb3TpYb3CYb3 = Matmul(ml2,TpYb3CYb3TpYb3CYb3) 
 ml2TpYb3CYb3TpYw3CYw3 = Matmul(ml2,TpYb3CYb3TpYw3CYw3) 
 ml2TpYeCYeTpYeCYe = Matmul(ml2,TpYeCYeTpYeCYe) 
 ml2TpYw3CYw3TpYb3CYb3 = Matmul(ml2,TpYw3CYw3TpYb3CYb3) 
 ml2TpYw3CYw3TpYw3CYw3 = Matmul(ml2,TpYw3CYw3TpYw3CYw3) 
 mq2adjYdYdadjYdYd = Matmul(mq2,adjYdYdadjYdYd) 
 mq2adjYdYdadjYuYu = Matmul(mq2,adjYdYdadjYuYu) 
 mq2adjYdTpYx3CYx3Yd = Matmul(mq2,adjYdTpYx3CYx3Yd) 
 mq2adjYuYuadjYdYd = Matmul(mq2,adjYuYuadjYdYd) 
 mq2adjYuYuadjYuYu = Matmul(mq2,adjYuYuadjYuYu) 
 mq2TpYdCYdTpYdCYd = Matmul(mq2,TpYdCYdTpYdCYd) 
 mq2TpYuCYuTpYuCYu = Matmul(mq2,TpYuCYuTpYuCYu) 
 mu2YuadjYdYdadjYu = Matmul(mu2,YuadjYdYdadjYu) 
 mu2YuadjYuYuadjYu = Matmul(mu2,YuadjYuYuadjYu) 
 MWM3CYw3TpYb3CYb3TpYw3 = Matmul(MWM3,CYw3TpYb3CYb3TpYw3) 
! MWM3CYw3TpYb3CYb3TpTYw3 = Matmul(MWM3,CYw3TpYb3CYb3TpTYw3) 
 MWM3CYw3TpYeCYeTpYw3 = Matmul(MWM3,CYw3TpYeCYeTpYw3) 
! MWM3CYw3TpYeCYeTpTYw3 = Matmul(MWM3,CYw3TpYeCYeTpTYw3) 
 MWM3CYw3TpYw3CYw3TpYw3 = Matmul(MWM3,CYw3TpYw3CYw3TpYw3) 
! MWM3CYw3TpYw3CYw3TpTYw3 = Matmul(MWM3,CYw3TpYw3CYw3TpTYw3) 
! MWM3CYw3TpTYb3CYb3TpYw3 = Matmul(MWM3,CYw3TpTYb3CYb3TpYw3) 
! MWM3CYw3TpTYeCYeTpYw3 = Matmul(MWM3,CYw3TpTYeCYeTpYw3) 
! MWM3CYw3TpTYw3CYw3TpYw3 = Matmul(MWM3,CYw3TpTYw3CYw3TpYw3) 
 MXM3CYx3YdadjYdTpYx3 = Matmul(MXM3,CYx3YdadjYdTpYx3) 
! MXM3CYx3YdadjYdTpTYx3 = Matmul(MXM3,CYx3YdadjYdTpTYx3) 
! MXM3CYx3TYdadjYdTpYx3 = Matmul(MXM3,CYx3TYdadjYdTpYx3) 
 MXM3CYx3TpYx3CYx3TpYx3 = Matmul(MXM3,CYx3TpYx3CYx3TpYx3) 
! MXM3CYx3TpYx3CYx3TpTYx3 = Matmul(MXM3,CYx3TpYx3CYx3TpTYx3) 
! MXM3CYx3TpTYx3CYx3TpYx3 = Matmul(MXM3,CYx3TpTYx3CYx3TpYx3) 
 Yb3ml2adjYb3Yb3adjYb3 = Matmul(Yb3,ml2adjYb3Yb3adjYb3) 
 Yb3ml2adjYeYeadjYb3 = Matmul(Yb3,ml2adjYeYeadjYb3) 
 Yb3ml2adjYw3Yw3adjYb3 = Matmul(Yb3,ml2adjYw3Yw3adjYb3) 
 Yb3adjYb3mHb32Yb3adjYb3 = Matmul(Yb3,adjYb3mHb32Yb3adjYb3) 
 Yb3adjYb3Yb3ml2adjYb3 = Matmul(Yb3,adjYb3Yb3ml2adjYb3) 
 Yb3adjYb3Yb3adjYb3MBM3 = Matmul(Yb3,adjYb3Yb3adjYb3MBM3) 
 Yb3adjYb3Yb3adjYb3mHb32 = Matmul(Yb3,adjYb3Yb3adjYb3mHb32) 
 Yb3adjYb3Yb3adjYb3Yb3 = Matmul(Yb3,adjYb3Yb3adjYb3Yb3) 
! Yb3adjYb3Yb3adjYb3BMBM3 = Matmul(Yb3,adjYb3Yb3adjYb3BMBM3) 
 Yb3adjYb3Yb3adjYb3TYb3 = Matmul(Yb3,adjYb3Yb3adjYb3TYb3) 
 Yb3adjYb3Yb3adjYw3Yw3 = Matmul(Yb3,adjYb3Yb3adjYw3Yw3) 
 Yb3adjYb3Yb3adjYw3TYw3 = Matmul(Yb3,adjYb3Yb3adjYw3TYw3) 
! Yb3adjYb3TYb3adjYb3MBM3 = Matmul(Yb3,adjYb3TYb3adjYb3MBM3) 
 Yb3adjYb3TYb3adjYb3Yb3 = Matmul(Yb3,adjYb3TYb3adjYb3Yb3) 
 Yb3adjYb3TYb3adjYw3Yw3 = Matmul(Yb3,adjYb3TYb3adjYw3Yw3) 
 Yb3adjYeme2YeadjYb3 = Matmul(Yb3,adjYeme2YeadjYb3) 
 Yb3adjYeYeml2adjYb3 = Matmul(Yb3,adjYeYeml2adjYb3) 
 Yb3adjYeYeadjYb3MBM3 = Matmul(Yb3,adjYeYeadjYb3MBM3) 
 Yb3adjYeYeadjYb3mHb32 = Matmul(Yb3,adjYeYeadjYb3mHb32) 
 Yb3adjYeYeadjYb3Yb3 = Matmul(Yb3,adjYeYeadjYb3Yb3) 
! Yb3adjYeYeadjYb3BMBM3 = Matmul(Yb3,adjYeYeadjYb3BMBM3) 
 Yb3adjYeYeadjYb3TYb3 = Matmul(Yb3,adjYeYeadjYb3TYb3) 
 Yb3adjYeYeadjYeYe = Matmul(Yb3,adjYeYeadjYeYe) 
 Yb3adjYeYeadjYeTYe = Matmul(Yb3,adjYeYeadjYeTYe) 
 Yb3adjYeYeadjYw3TYw3 = Matmul(Yb3,adjYeYeadjYw3TYw3) 
! Yb3adjYeTYeadjYb3MBM3 = Matmul(Yb3,adjYeTYeadjYb3MBM3) 
 Yb3adjYeTYeadjYb3Yb3 = Matmul(Yb3,adjYeTYeadjYb3Yb3) 
 Yb3adjYeTYeadjYeYe = Matmul(Yb3,adjYeTYeadjYeYe) 
 Yb3adjYeTYeadjYw3Yw3 = Matmul(Yb3,adjYeTYeadjYw3Yw3) 
 Yb3adjYw3mHw32Yw3adjYb3 = Matmul(Yb3,adjYw3mHw32Yw3adjYb3) 
 Yb3adjYw3Yw3ml2adjYb3 = Matmul(Yb3,adjYw3Yw3ml2adjYb3) 
 Yb3adjYw3Yw3adjYb3MBM3 = Matmul(Yb3,adjYw3Yw3adjYb3MBM3) 
 Yb3adjYw3Yw3adjYb3mHb32 = Matmul(Yb3,adjYw3Yw3adjYb3mHb32) 
 Yb3adjYw3Yw3adjYb3Yb3 = Matmul(Yb3,adjYw3Yw3adjYb3Yb3) 
! Yb3adjYw3Yw3adjYb3BMBM3 = Matmul(Yb3,adjYw3Yw3adjYb3BMBM3) 
 Yb3adjYw3Yw3adjYb3TYb3 = Matmul(Yb3,adjYw3Yw3adjYb3TYb3) 
 Yb3adjYw3Yw3adjYw3Yw3 = Matmul(Yb3,adjYw3Yw3adjYw3Yw3) 
 Yb3adjYw3Yw3adjYw3TYw3 = Matmul(Yb3,adjYw3Yw3adjYw3TYw3) 
! Yb3adjYw3TYw3adjYb3MBM3 = Matmul(Yb3,adjYw3TYw3adjYb3MBM3) 
 Yb3adjYw3TYw3adjYb3Yb3 = Matmul(Yb3,adjYw3TYw3adjYb3Yb3) 
 Yb3adjYw3TYw3adjYw3Yw3 = Matmul(Yb3,adjYw3TYw3adjYw3Yw3) 
 Ydmq2adjYdYdadjYd = Matmul(Yd,mq2adjYdYdadjYd) 
 Ydmq2adjYuYuadjYd = Matmul(Yd,mq2adjYuYuadjYd) 
 Ydmq2adjYx3Yx3CYd = Matmul(Yd,mq2adjYx3Yx3CYd) 
 YdadjYdmd2YdadjYd = Matmul(Yd,adjYdmd2YdadjYd) 
 YdadjYdYdmq2adjYd = Matmul(Yd,adjYdYdmq2adjYd) 
 YdadjYdYdadjYdmd2 = Matmul(Yd,adjYdYdadjYdmd2) 
 YdadjYdYdadjYdYd = Matmul(Yd,adjYdYdadjYdYd) 
 YdadjYdYdadjYdTYd = Matmul(Yd,adjYdYdadjYdTYd) 
 YdadjYdTYdadjYdYd = Matmul(Yd,adjYdTYdadjYdYd) 
 YdadjYdTpYx3CYx3Yd = Matmul(Yd,adjYdTpYx3CYx3Yd) 
 YdadjYdTpYx3CYx3TYd = Matmul(Yd,adjYdTpYx3CYx3TYd) 
 YdadjYdTpTYx3CYx3Yd = Matmul(Yd,adjYdTpTYx3CYx3Yd) 
 YdadjYumu2YuadjYd = Matmul(Yd,adjYumu2YuadjYd) 
 YdadjYuYumq2adjYd = Matmul(Yd,adjYuYumq2adjYd) 
 YdadjYuYuadjYdmd2 = Matmul(Yd,adjYuYuadjYdmd2) 
 YdadjYuYuadjYdYd = Matmul(Yd,adjYuYuadjYdYd) 
 YdadjYuYuadjYdTYd = Matmul(Yd,adjYuYuadjYdTYd) 
 YdadjYuYuadjYuYu = Matmul(Yd,adjYuYuadjYuYu) 
 YdadjYuYuadjYuTYu = Matmul(Yd,adjYuYuadjYuTYu) 
 YdadjYuTYuadjYdYd = Matmul(Yd,adjYuTYuadjYdYd) 
 YdadjYuTYuadjYuYu = Matmul(Yd,adjYuTYuadjYuYu) 
 Yeml2adjYb3Yb3adjYe = Matmul(Ye,ml2adjYb3Yb3adjYe) 
 Yeml2adjYeYeadjYe = Matmul(Ye,ml2adjYeYeadjYe) 
 Yeml2adjYw3Yw3adjYe = Matmul(Ye,ml2adjYw3Yw3adjYe) 
 YeadjYb3mHb32Yb3adjYe = Matmul(Ye,adjYb3mHb32Yb3adjYe) 
 YeadjYb3Yb3ml2adjYe = Matmul(Ye,adjYb3Yb3ml2adjYe) 
 YeadjYb3Yb3adjYb3Yb3 = Matmul(Ye,adjYb3Yb3adjYb3Yb3) 
 YeadjYb3Yb3adjYb3TYb3 = Matmul(Ye,adjYb3Yb3adjYb3TYb3) 
 YeadjYb3Yb3adjYeme2 = Matmul(Ye,adjYb3Yb3adjYeme2) 
 YeadjYb3Yb3adjYeYe = Matmul(Ye,adjYb3Yb3adjYeYe) 
 YeadjYb3Yb3adjYeTYe = Matmul(Ye,adjYb3Yb3adjYeTYe) 
 YeadjYb3Yb3adjYw3Yw3 = Matmul(Ye,adjYb3Yb3adjYw3Yw3) 
 YeadjYb3Yb3adjYw3TYw3 = Matmul(Ye,adjYb3Yb3adjYw3TYw3) 
 YeadjYb3TYb3adjYb3Yb3 = Matmul(Ye,adjYb3TYb3adjYb3Yb3) 
 YeadjYb3TYb3adjYeYe = Matmul(Ye,adjYb3TYb3adjYeYe) 
 YeadjYb3TYb3adjYw3Yw3 = Matmul(Ye,adjYb3TYb3adjYw3Yw3) 
 YeadjYeme2YeadjYe = Matmul(Ye,adjYeme2YeadjYe) 
 YeadjYeYeml2adjYe = Matmul(Ye,adjYeYeml2adjYe) 
 YeadjYeYeadjYeme2 = Matmul(Ye,adjYeYeadjYeme2) 
 YeadjYeYeadjYeYe = Matmul(Ye,adjYeYeadjYeYe) 
 YeadjYeYeadjYeTYe = Matmul(Ye,adjYeYeadjYeTYe) 
 YeadjYeTYeadjYeYe = Matmul(Ye,adjYeTYeadjYeYe) 
 YeadjYw3mHw32Yw3adjYe = Matmul(Ye,adjYw3mHw32Yw3adjYe) 
 YeadjYw3Yw3ml2adjYe = Matmul(Ye,adjYw3Yw3ml2adjYe) 
 YeadjYw3Yw3adjYb3Yb3 = Matmul(Ye,adjYw3Yw3adjYb3Yb3) 
 YeadjYw3Yw3adjYb3TYb3 = Matmul(Ye,adjYw3Yw3adjYb3TYb3) 
 YeadjYw3Yw3adjYeme2 = Matmul(Ye,adjYw3Yw3adjYeme2) 
 YeadjYw3Yw3adjYeYe = Matmul(Ye,adjYw3Yw3adjYeYe) 
 YeadjYw3Yw3adjYeTYe = Matmul(Ye,adjYw3Yw3adjYeTYe) 
 YeadjYw3Yw3adjYw3Yw3 = Matmul(Ye,adjYw3Yw3adjYw3Yw3) 
 YeadjYw3Yw3adjYw3TYw3 = Matmul(Ye,adjYw3Yw3adjYw3TYw3) 
 YeadjYw3TYw3adjYb3Yb3 = Matmul(Ye,adjYw3TYw3adjYb3Yb3) 
 YeadjYw3TYw3adjYeYe = Matmul(Ye,adjYw3TYw3adjYeYe) 
 YeadjYw3TYw3adjYw3Yw3 = Matmul(Ye,adjYw3TYw3adjYw3Yw3) 
 Yumq2adjYdYdadjYu = Matmul(Yu,mq2adjYdYdadjYu) 
 Yumq2adjYuYuadjYu = Matmul(Yu,mq2adjYuYuadjYu) 
 YuadjYdmd2YdadjYu = Matmul(Yu,adjYdmd2YdadjYu) 
 YuadjYdYdmq2adjYu = Matmul(Yu,adjYdYdmq2adjYu) 
 YuadjYdYdadjYdYd = Matmul(Yu,adjYdYdadjYdYd) 
 YuadjYdYdadjYdTYd = Matmul(Yu,adjYdYdadjYdTYd) 
 YuadjYdYdadjYumu2 = Matmul(Yu,adjYdYdadjYumu2) 
 YuadjYdYdadjYuYu = Matmul(Yu,adjYdYdadjYuYu) 
 YuadjYdYdadjYuTYu = Matmul(Yu,adjYdYdadjYuTYu) 
 YuadjYdTYdadjYdYd = Matmul(Yu,adjYdTYdadjYdYd) 
 YuadjYdTYdadjYuYu = Matmul(Yu,adjYdTYdadjYuYu) 
 YuadjYdTpYx3CYx3Yd = Matmul(Yu,adjYdTpYx3CYx3Yd) 
 YuadjYdTpYx3CYx3TYd = Matmul(Yu,adjYdTpYx3CYx3TYd) 
 YuadjYdTpTYx3CYx3Yd = Matmul(Yu,adjYdTpTYx3CYx3Yd) 
 YuadjYumu2YuadjYu = Matmul(Yu,adjYumu2YuadjYu) 
 YuadjYuYumq2adjYu = Matmul(Yu,adjYuYumq2adjYu) 
 YuadjYuYuadjYumu2 = Matmul(Yu,adjYuYuadjYumu2) 
 YuadjYuYuadjYuYu = Matmul(Yu,adjYuYuadjYuYu) 
 YuadjYuYuadjYuTYu = Matmul(Yu,adjYuYuadjYuTYu) 
 YuadjYuTYuadjYuYu = Matmul(Yu,adjYuTYuadjYuYu) 
 Yw3ml2adjYb3Yb3adjYw3 = Matmul(Yw3,ml2adjYb3Yb3adjYw3) 
 Yw3ml2adjYeYeadjYw3 = Matmul(Yw3,ml2adjYeYeadjYw3) 
 Yw3ml2adjYw3Yw3adjYw3 = Matmul(Yw3,ml2adjYw3Yw3adjYw3) 
 Yw3adjYb3mHb32Yb3adjYw3 = Matmul(Yw3,adjYb3mHb32Yb3adjYw3) 
 Yw3adjYb3Yb3ml2adjYw3 = Matmul(Yw3,adjYb3Yb3ml2adjYw3) 
 Yw3adjYb3Yb3adjYb3Yb3 = Matmul(Yw3,adjYb3Yb3adjYb3Yb3) 
 Yw3adjYb3Yb3adjYb3TYb3 = Matmul(Yw3,adjYb3Yb3adjYb3TYb3) 
 Yw3adjYb3Yb3adjYw3mHw32 = Matmul(Yw3,adjYb3Yb3adjYw3mHw32) 
 Yw3adjYb3Yb3adjYw3MWM3 = Matmul(Yw3,adjYb3Yb3adjYw3MWM3) 
 Yw3adjYb3Yb3adjYw3Yw3 = Matmul(Yw3,adjYb3Yb3adjYw3Yw3) 
! Yw3adjYb3Yb3adjYw3BMWM3 = Matmul(Yw3,adjYb3Yb3adjYw3BMWM3) 
 Yw3adjYb3Yb3adjYw3TYw3 = Matmul(Yw3,adjYb3Yb3adjYw3TYw3) 
 Yw3adjYb3TYb3adjYb3Yb3 = Matmul(Yw3,adjYb3TYb3adjYb3Yb3) 
! Yw3adjYb3TYb3adjYw3MWM3 = Matmul(Yw3,adjYb3TYb3adjYw3MWM3) 
 Yw3adjYb3TYb3adjYw3Yw3 = Matmul(Yw3,adjYb3TYb3adjYw3Yw3) 
 Yw3adjYeme2YeadjYw3 = Matmul(Yw3,adjYeme2YeadjYw3) 
 Yw3adjYeYeml2adjYw3 = Matmul(Yw3,adjYeYeml2adjYw3) 
 Yw3adjYeYeadjYb3TYb3 = Matmul(Yw3,adjYeYeadjYb3TYb3) 
 Yw3adjYeYeadjYeYe = Matmul(Yw3,adjYeYeadjYeYe) 
 Yw3adjYeYeadjYeTYe = Matmul(Yw3,adjYeYeadjYeTYe) 
 Yw3adjYeYeadjYw3mHw32 = Matmul(Yw3,adjYeYeadjYw3mHw32) 
 Yw3adjYeYeadjYw3MWM3 = Matmul(Yw3,adjYeYeadjYw3MWM3) 
 Yw3adjYeYeadjYw3Yw3 = Matmul(Yw3,adjYeYeadjYw3Yw3) 
! Yw3adjYeYeadjYw3BMWM3 = Matmul(Yw3,adjYeYeadjYw3BMWM3) 
 Yw3adjYeYeadjYw3TYw3 = Matmul(Yw3,adjYeYeadjYw3TYw3) 
 Yw3adjYeTYeadjYb3Yb3 = Matmul(Yw3,adjYeTYeadjYb3Yb3) 
 Yw3adjYeTYeadjYeYe = Matmul(Yw3,adjYeTYeadjYeYe) 
! Yw3adjYeTYeadjYw3MWM3 = Matmul(Yw3,adjYeTYeadjYw3MWM3) 
 Yw3adjYeTYeadjYw3Yw3 = Matmul(Yw3,adjYeTYeadjYw3Yw3) 
 Yw3adjYw3mHw32Yw3adjYw3 = Matmul(Yw3,adjYw3mHw32Yw3adjYw3) 
 Yw3adjYw3Yw3ml2adjYw3 = Matmul(Yw3,adjYw3Yw3ml2adjYw3) 
 Yw3adjYw3Yw3adjYb3Yb3 = Matmul(Yw3,adjYw3Yw3adjYb3Yb3) 
 Yw3adjYw3Yw3adjYb3TYb3 = Matmul(Yw3,adjYw3Yw3adjYb3TYb3) 
 Yw3adjYw3Yw3adjYw3mHw32 = Matmul(Yw3,adjYw3Yw3adjYw3mHw32) 
 Yw3adjYw3Yw3adjYw3MWM3 = Matmul(Yw3,adjYw3Yw3adjYw3MWM3) 
 Yw3adjYw3Yw3adjYw3Yw3 = Matmul(Yw3,adjYw3Yw3adjYw3Yw3) 
! Yw3adjYw3Yw3adjYw3BMWM3 = Matmul(Yw3,adjYw3Yw3adjYw3BMWM3) 
 Yw3adjYw3Yw3adjYw3TYw3 = Matmul(Yw3,adjYw3Yw3adjYw3TYw3) 
 Yw3adjYw3TYw3adjYb3Yb3 = Matmul(Yw3,adjYw3TYw3adjYb3Yb3) 
! Yw3adjYw3TYw3adjYw3MWM3 = Matmul(Yw3,adjYw3TYw3adjYw3MWM3) 
 Yw3adjYw3TYw3adjYw3Yw3 = Matmul(Yw3,adjYw3TYw3adjYw3Yw3) 
 Yx3md2adjYx3Yx3adjYx3 = Matmul(Yx3,md2adjYx3Yx3adjYx3) 
 Yx3md2CYdTpYdadjYx3 = Matmul(Yx3,md2CYdTpYdadjYx3) 
 Yx3adjYx3mHxb32Yx3adjYx3 = Matmul(Yx3,adjYx3mHxb32Yx3adjYx3) 
 Yx3adjYx3Yx3md2adjYx3 = Matmul(Yx3,adjYx3Yx3md2adjYx3) 
 Yx3adjYx3Yx3adjYx3mHxb32 = Matmul(Yx3,adjYx3Yx3adjYx3mHxb32) 
 Yx3adjYx3Yx3adjYx3Yx3 = Matmul(Yx3,adjYx3Yx3adjYx3Yx3) 
 Yx3adjYx3Yx3adjYx3TYx3 = Matmul(Yx3,adjYx3Yx3adjYx3TYx3) 
 Yx3adjYx3TYx3adjYx3Yx3 = Matmul(Yx3,adjYx3TYx3adjYx3Yx3) 
 Yx3CYdmq2TpYdadjYx3 = Matmul(Yx3,CYdmq2TpYdadjYx3) 
 Yx3CYdTpYdmd2adjYx3 = Matmul(Yx3,CYdTpYdmd2adjYx3) 
 Yx3CYdTpYdadjYx3mHxb32 = Matmul(Yx3,CYdTpYdadjYx3mHxb32) 
 Yx3CYdTpYdadjYx3Yx3 = Matmul(Yx3,CYdTpYdadjYx3Yx3) 
 Yx3CYdTpYdadjYx3TYx3 = Matmul(Yx3,CYdTpYdadjYx3TYx3) 
 Yx3CYdTpYdCYdTpYd = Matmul(Yx3,CYdTpYdCYdTpYd) 
 Yx3CYdTpYdCYdTpTYd = Matmul(Yx3,CYdTpYdCYdTpTYd) 
 Yx3CYdTpYuCYuTpYd = Matmul(Yx3,CYdTpYuCYuTpYd) 
 Yx3CYdTpYuCYuTpTYd = Matmul(Yx3,CYdTpYuCYuTpTYd) 
 Yx3CYdTpTYdCYdTpYd = Matmul(Yx3,CYdTpTYdCYdTpYd) 
 Yx3CYdTpTYuCYuTpYd = Matmul(Yx3,CYdTpTYuCYuTpYd) 
 Yx3TYdadjYdadjYx3Yx3 = Matmul(Yx3,TYdadjYdadjYx3Yx3) 
! BMBM3CYb3TpYb3CYb3TpYb3 = Matmul(BMBM3,CYb3TpYb3CYb3TpYb3) 
! BMBM3CYb3TpYeCYeTpYb3 = Matmul(BMBM3,CYb3TpYeCYeTpYb3) 
! BMBM3CYb3TpYw3CYw3TpYb3 = Matmul(BMBM3,CYb3TpYw3CYw3TpYb3) 
! BMWM3CYw3TpYb3CYb3TpYw3 = Matmul(BMWM3,CYw3TpYb3CYb3TpYw3) 
! BMWM3CYw3TpYeCYeTpYw3 = Matmul(BMWM3,CYw3TpYeCYeTpYw3) 
! BMWM3CYw3TpYw3CYw3TpYw3 = Matmul(BMWM3,CYw3TpYw3CYw3TpYw3) 
! BMXM3CYx3YdadjYdTpYx3 = Matmul(BMXM3,CYx3YdadjYdTpYx3) 
! BMXM3CYx3TpYx3CYx3TpYx3 = Matmul(BMXM3,CYx3TpYx3CYx3TpYx3) 
! TYb3adjYb3Yb3adjYb3MBM3 = Matmul(TYb3,adjYb3Yb3adjYb3MBM3) 
 TYb3adjYb3Yb3adjYb3Yb3 = Matmul(TYb3,adjYb3Yb3adjYb3Yb3) 
 TYb3adjYb3Yb3adjYw3Yw3 = Matmul(TYb3,adjYb3Yb3adjYw3Yw3) 
! TYb3adjYeYeadjYb3MBM3 = Matmul(TYb3,adjYeYeadjYb3MBM3) 
 TYb3adjYeYeadjYb3Yb3 = Matmul(TYb3,adjYeYeadjYb3Yb3) 
 TYb3adjYeYeadjYeYe = Matmul(TYb3,adjYeYeadjYeYe) 
 TYb3adjYeYeadjYw3Yw3 = Matmul(TYb3,adjYeYeadjYw3Yw3) 
! TYb3adjYw3Yw3adjYb3MBM3 = Matmul(TYb3,adjYw3Yw3adjYb3MBM3) 
 TYb3adjYw3Yw3adjYb3Yb3 = Matmul(TYb3,adjYw3Yw3adjYb3Yb3) 
 TYb3adjYw3Yw3adjYw3Yw3 = Matmul(TYb3,adjYw3Yw3adjYw3Yw3) 
 TYdadjYdYdadjYdYd = Matmul(TYd,adjYdYdadjYdYd) 
 TYdadjYdTpYx3CYx3Yd = Matmul(TYd,adjYdTpYx3CYx3Yd) 
 TYdadjYuYuadjYdYd = Matmul(TYd,adjYuYuadjYdYd) 
 TYdadjYuYuadjYuYu = Matmul(TYd,adjYuYuadjYuYu) 
 TYeadjYb3Yb3adjYb3Yb3 = Matmul(TYe,adjYb3Yb3adjYb3Yb3) 
 TYeadjYb3Yb3adjYeYe = Matmul(TYe,adjYb3Yb3adjYeYe) 
 TYeadjYb3Yb3adjYw3Yw3 = Matmul(TYe,adjYb3Yb3adjYw3Yw3) 
 TYeadjYeYeadjYeYe = Matmul(TYe,adjYeYeadjYeYe) 
 TYeadjYw3Yw3adjYb3Yb3 = Matmul(TYe,adjYw3Yw3adjYb3Yb3) 
 TYeadjYw3Yw3adjYeYe = Matmul(TYe,adjYw3Yw3adjYeYe) 
 TYeadjYw3Yw3adjYw3Yw3 = Matmul(TYe,adjYw3Yw3adjYw3Yw3) 
 TYuadjYdYdadjYdYd = Matmul(TYu,adjYdYdadjYdYd) 
 TYuadjYdYdadjYuYu = Matmul(TYu,adjYdYdadjYuYu) 
 TYuadjYdTpYx3CYx3Yd = Matmul(TYu,adjYdTpYx3CYx3Yd) 
 TYuadjYuYuadjYuYu = Matmul(TYu,adjYuYuadjYuYu) 
 TYw3adjYb3Yb3adjYb3Yb3 = Matmul(TYw3,adjYb3Yb3adjYb3Yb3) 
! TYw3adjYb3Yb3adjYw3MWM3 = Matmul(TYw3,adjYb3Yb3adjYw3MWM3) 
 TYw3adjYb3Yb3adjYw3Yw3 = Matmul(TYw3,adjYb3Yb3adjYw3Yw3) 
 TYw3adjYeYeadjYb3Yb3 = Matmul(TYw3,adjYeYeadjYb3Yb3) 
 TYw3adjYeYeadjYeYe = Matmul(TYw3,adjYeYeadjYeYe) 
! TYw3adjYeYeadjYw3MWM3 = Matmul(TYw3,adjYeYeadjYw3MWM3) 
 TYw3adjYeYeadjYw3Yw3 = Matmul(TYw3,adjYeYeadjYw3Yw3) 
 TYw3adjYw3Yw3adjYb3Yb3 = Matmul(TYw3,adjYw3Yw3adjYb3Yb3) 
! TYw3adjYw3Yw3adjYw3MWM3 = Matmul(TYw3,adjYw3Yw3adjYw3MWM3) 
 TYw3adjYw3Yw3adjYw3Yw3 = Matmul(TYw3,adjYw3Yw3adjYw3Yw3) 
 TYx3adjYx3Yx3adjYx3Yx3 = Matmul(TYx3,adjYx3Yx3adjYx3Yx3) 
 TYx3CYdTpYdadjYx3Yx3 = Matmul(TYx3,CYdTpYdadjYx3Yx3) 
 TYx3CYdTpYdCYdTpYd = Matmul(TYx3,CYdTpYdCYdTpYd) 
 TYx3CYdTpYuCYuTpYd = Matmul(TYx3,CYdTpYuCYuTpYd) 
 TpYb3mHb32CYb3TpYb3CYb3 = Matmul(Transpose(Yb3),mHb32CYb3TpYb3CYb3) 
 TpYb3mHb32CYb3TpYw3CYw3 = Matmul(Transpose(Yb3),mHb32CYb3TpYw3CYw3) 
 TpYb3CYb3ml2TpYb3CYb3 = Matmul(Transpose(Yb3),CYb3ml2TpYb3CYb3) 
 TpYb3CYb3ml2TpYw3CYw3 = Matmul(Transpose(Yb3),CYb3ml2TpYw3CYw3) 
 TpYb3CYb3TpYb3mHb32CYb3 = Matmul(Transpose(Yb3),CYb3TpYb3mHb32CYb3) 
 TpYb3CYb3TpYb3CYb3ml2 = Matmul(Transpose(Yb3),CYb3TpYb3CYb3ml2) 
 TpYb3CYb3TpYw3mHw32CYw3 = Matmul(Transpose(Yb3),CYb3TpYw3mHw32CYw3) 
 TpYb3CYb3TpYw3CYw3ml2 = Matmul(Transpose(Yb3),CYb3TpYw3CYw3ml2) 
 TpYdmd2adjYx3Yx3CYd = Matmul(Transpose(Yd),md2adjYx3Yx3CYd) 
 TpYdmd2CYdTpYdCYd = Matmul(Transpose(Yd),md2CYdTpYdCYd) 
 TpYdadjYx3mHxb32Yx3CYd = Matmul(Transpose(Yd),adjYx3mHxb32Yx3CYd) 
 TpYdadjYx3Yx3md2CYd = Matmul(Transpose(Yd),adjYx3Yx3md2CYd) 
 TpYdadjYx3Yx3CYdmq2 = Matmul(Transpose(Yd),adjYx3Yx3CYdmq2) 
 TpYdCYdmq2TpYdCYd = Matmul(Transpose(Yd),CYdmq2TpYdCYd) 
 TpYdCYdTpYdmd2CYd = Matmul(Transpose(Yd),CYdTpYdmd2CYd) 
 TpYdCYdTpYdCYdmq2 = Matmul(Transpose(Yd),CYdTpYdCYdmq2) 
 TpYeme2CYeTpYeCYe = Matmul(Transpose(Ye),me2CYeTpYeCYe) 
 TpYeCYeml2TpYeCYe = Matmul(Transpose(Ye),CYeml2TpYeCYe) 
 TpYeCYeTpYeme2CYe = Matmul(Transpose(Ye),CYeTpYeme2CYe) 
 TpYeCYeTpYeCYeml2 = Matmul(Transpose(Ye),CYeTpYeCYeml2) 
 TpYumu2CYuTpYuCYu = Matmul(Transpose(Yu),mu2CYuTpYuCYu) 
 TpYuCYumq2TpYuCYu = Matmul(Transpose(Yu),CYumq2TpYuCYu) 
 TpYuCYuTpYumu2CYu = Matmul(Transpose(Yu),CYuTpYumu2CYu) 
 TpYuCYuTpYuCYumq2 = Matmul(Transpose(Yu),CYuTpYuCYumq2) 
 TpYw3mHw32CYw3TpYb3CYb3 = Matmul(Transpose(Yw3),mHw32CYw3TpYb3CYb3) 
 TpYw3mHw32CYw3TpYw3CYw3 = Matmul(Transpose(Yw3),mHw32CYw3TpYw3CYw3) 
 TpYw3CYw3ml2TpYb3CYb3 = Matmul(Transpose(Yw3),CYw3ml2TpYb3CYb3) 
 TpYw3CYw3ml2TpYw3CYw3 = Matmul(Transpose(Yw3),CYw3ml2TpYw3CYw3) 
 TpYw3CYw3TpYb3mHb32CYb3 = Matmul(Transpose(Yw3),CYw3TpYb3mHb32CYb3) 
 TpYw3CYw3TpYb3CYb3ml2 = Matmul(Transpose(Yw3),CYw3TpYb3CYb3ml2) 
 TpYw3CYw3TpYw3mHw32CYw3 = Matmul(Transpose(Yw3),CYw3TpYw3mHw32CYw3) 
 TpYw3CYw3TpYw3CYw3ml2 = Matmul(Transpose(Yw3),CYw3TpYw3CYw3ml2) 
 TpYx3mHxb32CYx3TpYx3CYx3 = Matmul(Transpose(Yx3),mHxb32CYx3TpYx3CYx3) 
 TpYx3CYx3md2TpYx3CYx3 = Matmul(Transpose(Yx3),CYx3md2TpYx3CYx3) 
 TpYx3CYx3TpYx3mHxb32CYx3 = Matmul(Transpose(Yx3),CYx3TpYx3mHxb32CYx3) 
 TpYx3CYx3TpYx3CYx3md2 = Matmul(Transpose(Yx3),CYx3TpYx3CYx3md2) 
 TpYx3CYx3TpYx3CYx3Yd = Matmul(Transpose(Yx3),CYx3TpYx3CYx3Yd) 
 TpYx3CYx3TpYx3CYx3TYd = Matmul(Transpose(Yx3),CYx3TpYx3CYx3TYd) 
 TpYx3CYx3TpTYx3CYx3Yd = Matmul(Transpose(Yx3),CYx3TpTYx3CYx3Yd) 
 TpTYx3CYx3TpYx3CYx3Yd = Matmul(Transpose(TYx3),CYx3TpYx3CYx3Yd) 
 TrmHg32 = cTrace(mHg32) 
 TrmHw32 = cTrace(mHw32) 
 TrCTYb3TpYb3 = cTrace(CTYb3TpYb3) 
 TrCTYdTpYd = cTrace(CTYdTpYd) 
 TrCTYeTpYe = cTrace(CTYeTpYe) 
 TrCTYuTpYu = cTrace(CTYuTpYu) 
 TrCTYw3TpYw3 = cTrace(CTYw3TpYw3) 
 TrCTYx3TpYx3 = cTrace(CTYx3TpYx3) 
 TrYb3adjYb3Yb3adjYb3 = cTrace(Yb3adjYb3Yb3adjYb3) 
 TrYb3adjYb3TYb3adjYb3 = cTrace(Yb3adjYb3TYb3adjYb3) 
 TrYb3adjYb3TYb3adjTYb3 = cTrace(Yb3adjYb3TYb3adjTYb3) 
 TrYb3adjYeYeadjYb3 = cTrace(Yb3adjYeYeadjYb3) 
 TrYb3adjYeTYeadjYb3 = cTrace(Yb3adjYeTYeadjYb3) 
 TrYb3adjYeTYeadjTYb3 = cTrace(Yb3adjYeTYeadjTYb3) 
 TrYb3adjYw3Yw3adjYb3 = cTrace(Yb3adjYw3Yw3adjYb3) 
 TrYb3adjYw3TYw3adjYb3 = cTrace(Yb3adjYw3TYw3adjYb3) 
 TrYb3adjYw3TYw3adjTYb3 = cTrace(Yb3adjYw3TYw3adjTYb3) 
 TrYb3adjTYb3TYb3adjYb3 = cTrace(Yb3adjTYb3TYb3adjYb3) 
 TrYb3adjTYeTYeadjYb3 = cTrace(Yb3adjTYeTYeadjYb3) 
 TrYb3adjTYw3TYw3adjYb3 = cTrace(Yb3adjTYw3TYw3adjYb3) 
 TrYb3TpTYb3CTYb3adjYb3 = cTrace(Yb3TpTYb3CTYb3adjYb3) 
 TrYb3TpTYeCTYeadjYb3 = cTrace(Yb3TpTYeCTYeadjYb3) 
 TrYb3TpTYw3CTYw3adjYb3 = cTrace(Yb3TpTYw3CTYw3adjYb3) 
 TrYdadjYdYdadjYd = cTrace(YdadjYdYdadjYd) 
 TrYdadjYdTYdadjYd = cTrace(YdadjYdTYdadjYd) 
 TrYdadjYdTYdadjTYd = cTrace(YdadjYdTYdadjTYd) 
 TrYdadjYdTpYx3CYx3 = cTrace(YdadjYdTpYx3CYx3) 
 TrYdadjYdTpTYx3CTYx3 = cTrace(YdadjYdTpTYx3CTYx3) 
 TrYdadjYuYuadjYd = cTrace(YdadjYuYuadjYd) 
 TrYdadjYuTYuadjYd = cTrace(YdadjYuTYuadjYd) 
 TrYdadjYuTYuadjTYd = cTrace(YdadjYuTYuadjTYd) 
 TrYdadjTYdTYdadjYd = cTrace(YdadjTYdTYdadjYd) 
 TrYdadjTYuTYuadjYd = cTrace(YdadjTYuTYuadjYd) 
 TrYdTpTYdCTYdadjYd = cTrace(YdTpTYdCTYdadjYd) 
 TrYdTpTYuCTYuadjYd = cTrace(YdTpTYuCTYuadjYd) 
 TrYeadjYb3TYb3adjYe = cTrace(YeadjYb3TYb3adjYe) 
 TrYeadjYb3TYb3adjTYe = cTrace(YeadjYb3TYb3adjTYe) 
 TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe) 
 TrYeadjYeTYeadjYe = cTrace(YeadjYeTYeadjYe) 
 TrYeadjYeTYeadjTYe = cTrace(YeadjYeTYeadjTYe) 
 TrYeadjYw3Yw3adjYe = cTrace(YeadjYw3Yw3adjYe) 
 TrYeadjYw3TYw3adjYe = cTrace(YeadjYw3TYw3adjYe) 
 TrYeadjYw3TYw3adjTYe = cTrace(YeadjYw3TYw3adjTYe) 
 TrYeadjTYb3TYb3adjYe = cTrace(YeadjTYb3TYb3adjYe) 
 TrYeadjTYeTYeadjYe = cTrace(YeadjTYeTYeadjYe) 
 TrYeadjTYw3TYw3adjYe = cTrace(YeadjTYw3TYw3adjYe) 
 TrYeTpTYb3CTYb3adjYe = cTrace(YeTpTYb3CTYb3adjYe) 
 TrYeTpTYeCTYeadjYe = cTrace(YeTpTYeCTYeadjYe) 
 TrYeTpTYw3CTYw3adjYe = cTrace(YeTpTYw3CTYw3adjYe) 
 TrYuadjYdTYdadjYu = cTrace(YuadjYdTYdadjYu) 
 TrYuadjYdTYdadjTYu = cTrace(YuadjYdTYdadjTYu) 
 TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu) 
 TrYuadjYuTYuadjYu = cTrace(YuadjYuTYuadjYu) 
 TrYuadjYuTYuadjTYu = cTrace(YuadjYuTYuadjTYu) 
 TrYuadjTYdTYdadjYu = cTrace(YuadjTYdTYdadjYu) 
 TrYuadjTYuTYuadjYu = cTrace(YuadjTYuTYuadjYu) 
 TrYuTpTYdCTYdadjYu = cTrace(YuTpTYdCTYdadjYu) 
 TrYuTpTYuCTYuadjYu = cTrace(YuTpTYuCTYuadjYu) 
 TrYw3adjYb3TYb3adjYw3 = cTrace(Yw3adjYb3TYb3adjYw3) 
 TrYw3adjYb3TYb3adjTYw3 = cTrace(Yw3adjYb3TYb3adjTYw3) 
 TrYw3adjYeTYeadjYw3 = cTrace(Yw3adjYeTYeadjYw3) 
 TrYw3adjYeTYeadjTYw3 = cTrace(Yw3adjYeTYeadjTYw3) 
 TrYw3adjYw3Yw3adjYw3 = cTrace(Yw3adjYw3Yw3adjYw3) 
 TrYw3adjYw3TYw3adjYw3 = cTrace(Yw3adjYw3TYw3adjYw3) 
 TrYw3adjYw3TYw3adjTYw3 = cTrace(Yw3adjYw3TYw3adjTYw3) 
 TrYw3adjTYb3TYb3adjYw3 = cTrace(Yw3adjTYb3TYb3adjYw3) 
 TrYw3adjTYeTYeadjYw3 = cTrace(Yw3adjTYeTYeadjYw3) 
 TrYw3adjTYw3TYw3adjYw3 = cTrace(Yw3adjTYw3TYw3adjYw3) 
 TrYw3TpTYb3CTYb3adjYw3 = cTrace(Yw3TpTYb3CTYb3adjYw3) 
 TrYw3TpTYeCTYeadjYw3 = cTrace(Yw3TpTYeCTYeadjYw3) 
 TrYw3TpTYw3CTYw3adjYw3 = cTrace(Yw3TpTYw3CTYw3adjYw3) 
 TrYx3adjYx3Yx3adjYx3 = cTrace(Yx3adjYx3Yx3adjYx3) 
 TrYx3adjYx3TYx3adjYx3 = cTrace(Yx3adjYx3TYx3adjYx3) 
 TrYx3adjYx3TYx3adjTYx3 = cTrace(Yx3adjYx3TYx3adjTYx3) 
 TrYx3adjTYx3TYx3adjYx3 = cTrace(Yx3adjTYx3TYx3adjYx3) 
 TrYx3CTYdTpTYdadjYx3 = cTrace(Yx3CTYdTpTYdadjYx3) 
 TrYx3TpTYx3CTYx3adjYx3 = cTrace(Yx3TpTYx3CTYx3adjYx3) 
 TradjYdTpYx3CYx3TYd = cTrace(adjYdTpYx3CYx3TYd) 
 TradjYdTpYx3CTYx3TYd = cTrace(adjYdTpYx3CTYx3TYd) 
 TradjYdTpYx3CTYx3TpTYd = cTrace(adjYdTpYx3CTYx3TpTYd) 
 TradjYdTpTYx3CTYx3TpYd = cTrace(adjYdTpTYx3CTYx3TpYd) 
 TradjYx3TYx3CYdTpYd = cTrace(adjYx3TYx3CYdTpYd) 
 TradjYx3TYx3CTYdTpYd = cTrace(adjYx3TYx3CTYdTpYd) 
 Trmd2YdadjYdYdadjYd = cTrace(md2YdadjYdYdadjYd) 
 Trmd2YdadjYdTpYx3CYx3 = cTrace(md2YdadjYdTpYx3CYx3) 
 Trmd2YdadjYuYuadjYd = cTrace(md2YdadjYuYuadjYd) 
 Trmd2adjYx3Yx3adjYx3Yx3 = cTrace(md2adjYx3Yx3adjYx3Yx3) 
 Trmd2TpYx3CYx3YdadjYd = cTrace(md2TpYx3CYx3YdadjYd) 
 Trme2YeadjYb3Yb3adjYe = cTrace(me2YeadjYb3Yb3adjYe) 
 Trme2YeadjYeYeadjYe = cTrace(me2YeadjYeYeadjYe) 
 Trme2YeadjYw3Yw3adjYe = cTrace(me2YeadjYw3Yw3adjYe) 
 TrmHb32Yb3adjYb3Yb3adjYb3 = cTrace(mHb32Yb3adjYb3Yb3adjYb3) 
 TrmHb32Yb3adjYeYeadjYb3 = cTrace(mHb32Yb3adjYeYeadjYb3) 
 TrmHb32Yb3adjYw3Yw3adjYb3 = cTrace(mHb32Yb3adjYw3Yw3adjYb3) 
 TrmHw32Yw3adjYb3Yb3adjYw3 = cTrace(mHw32Yw3adjYb3Yb3adjYw3) 
 TrmHw32Yw3adjYeYeadjYw3 = cTrace(mHw32Yw3adjYeYeadjYw3) 
 TrmHw32Yw3adjYw3Yw3adjYw3 = cTrace(mHw32Yw3adjYw3Yw3adjYw3) 
 TrmHxb32Yx3adjYx3Yx3adjYx3 = cTrace(mHxb32Yx3adjYx3Yx3adjYx3) 
 TrmHxb32Yx3CYdTpYdadjYx3 = cTrace(mHxb32Yx3CYdTpYdadjYx3) 
 TrmHxb32CYx3YdadjYdTpYx3 = cTrace(mHxb32CYx3YdadjYdTpYx3) 
 Trml2adjYb3Yb3adjYb3Yb3 = cTrace(ml2adjYb3Yb3adjYb3Yb3) 
 Trml2adjYb3Yb3adjYeYe = cTrace(ml2adjYb3Yb3adjYeYe) 
 Trml2adjYb3Yb3adjYw3Yw3 = cTrace(ml2adjYb3Yb3adjYw3Yw3) 
 Trml2adjYeYeadjYb3Yb3 = cTrace(ml2adjYeYeadjYb3Yb3) 
 Trml2adjYeYeadjYeYe = cTrace(ml2adjYeYeadjYeYe) 
 Trml2adjYeYeadjYw3Yw3 = cTrace(ml2adjYeYeadjYw3Yw3) 
 Trml2adjYw3Yw3adjYb3Yb3 = cTrace(ml2adjYw3Yw3adjYb3Yb3) 
 Trml2adjYw3Yw3adjYeYe = cTrace(ml2adjYw3Yw3adjYeYe) 
 Trml2adjYw3Yw3adjYw3Yw3 = cTrace(ml2adjYw3Yw3adjYw3Yw3) 
 Trmq2adjYdYdadjYdYd = cTrace(mq2adjYdYdadjYdYd) 
 Trmq2adjYdYdadjYuYu = cTrace(mq2adjYdYdadjYuYu) 
 Trmq2adjYdTpYx3CYx3Yd = cTrace(mq2adjYdTpYx3CYx3Yd) 
 Trmq2adjYuYuadjYdYd = cTrace(mq2adjYuYuadjYdYd) 
 Trmq2adjYuYuadjYuYu = cTrace(mq2adjYuYuadjYuYu) 
 Trmu2YuadjYdYdadjYu = cTrace(mu2YuadjYdYdadjYu) 
 Trmu2YuadjYuYuadjYu = cTrace(mu2YuadjYuYuadjYu) 
 g1p4 =g1**4 
 g2p4 =g2**4 
 g3p4 =g3**4 
End If 
 
 
Tr1(1) = (3*(-1._dp*(mHd2) + mHu2 + Trmd2 + Trme2 + 5._dp*(TrmHx32) - 5._dp*(TrmHxb32)& 
&  - Trml2 + Trmq2 - 2._dp*(Trmu2)))/5._dp

If (TwoLoopRGE) Then 
Tr2(1) = (3._dp*(mHd2) + 3._dp*(mHu2) + 2._dp*(Trmd2) + 6._dp*(Trme2) +               & 
&  25._dp*(TrmHx32) + 25._dp*(TrmHxb32) + 3._dp*(Trml2) + Trmq2 + 8._dp*(Trmu2))/10._dp

Tr3(1) = (3*(-20._dp*(Trmd2adjYx3Yx3) - 20._dp*(Trmd2YdadjYd) - 20._dp*(Trme2YeadjYe) & 
&  + 50._dp*(TrmHxb32Yx3adjYx3) + 3._dp*(Trml2adjYb3Yb3) + 10._dp*(Trml2adjYeYe)         & 
&  + 15._dp*(Trml2adjYw3Yw3) - 15*g2p2*(mHd2 - mHu2 - 5._dp*(TrmHx32) + 5._dp*(TrmHxb32) & 
&  + Trml2 - Trmq2) - 10._dp*(Trmq2adjYdYd) - 10._dp*(Trmq2adjYuYu) + (80*g3p2*(Trmd2 +  & 
&  5._dp*(TrmHx32) - 5._dp*(TrmHxb32) + Trmq2 - 2._dp*(Trmu2)))/3._dp - (g1p2*(9._dp*(mHd2)& 
&  - 9._dp*(mHu2) - 4._dp*(Trmd2) - 36._dp*(Trme2) - 125._dp*(TrmHx32) + 125._dp*(TrmHxb32)& 
&  + 9._dp*(Trml2) - Trmq2 + 32._dp*(Trmu2)))/3._dp + 40._dp*(Trmu2YuadjYu)              & 
&  - 3*mHu2*TrYb3adjYb3 + 30*mHd2*TrYdadjYd + 10*mHd2*TrYeadjYe - 30*mHu2*TrYuadjYu -    & 
&  15*mHu2*TrYw3adjYw3 - 30*mHu2*TrYx3adjYx3))/100._dp

Tr2(2) = (mHd2 + mHu2 + 6._dp*(TrmHw32) + 3._dp*(TrmHx32) + 3._dp*(TrmHxb32)          & 
&  + Trml2 + 3._dp*(Trmq2))/2._dp

Tr2(3) = Trmd2/2._dp + 8._dp*(TrmHg32) + TrmHx32 + TrmHxb32 + Trmq2 + Trmu2/2._dp

End If 
 
 
!-------------------- 
! g1 
!-------------------- 
 
betag11  = (g1p3*(66 + 25*NGHx3 + 25*NGHxb3))/10._dp

 
 
If (TwoLoopRGE) Then 
betag12 = (g1p3*(6*(199._dp*(g1p2) + 135._dp*(g2p2) + 440._dp*(g3p2) - 9._dp*(TrYb3adjYb3) -    & 
&  70._dp*(TrYdadjYd) - 90._dp*(TrYeadjYe) - 130._dp*(TrYuadjYu) - 60._dp*(TrYw3adjYw3) -& 
&  190._dp*(TrYx3adjYx3)) + 125*(5._dp*(g1p2) + 9._dp*(g2p2) + 16._dp*(g3p2))*(NGHx3 +& 
&  NGHxb3)))/150._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21  = (g2p3*(2 + 4*NGHw3 + 3*NGHx3 +               & 
&  3*NGHxb3))/2._dp

 
 
If (TwoLoopRGE) Then 
betag22 = (g2p3*(2*(27._dp*(g1p2) + 375._dp*(g2p2) + 360._dp*(g3p2) - 9._dp*(TrYb3adjYb3) -     & 
&  90._dp*(TrYdadjYd) - 30._dp*(TrYeadjYe) - 90._dp*(TrYuadjYu) - 140._dp*(TrYw3adjYw3) -& 
&  90._dp*(TrYx3adjYx3)) + 720*g2p2*NGHw3 + 15*(5._dp*(g1p2) +            & 
&  21._dp*(g2p2) + 16._dp*(g3p2))*(NGHx3 + NGHxb3)))/30._dp

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31  = g3p3*(-3 + 3*NGHg3 + NGHx3 + NGHxb3)

 
 
If (TwoLoopRGE) Then 
betag32 = (g3p3*(3*(11._dp*(g1p2) + 45._dp*(g2p2) + 70._dp*(g3p2) - 20._dp*(TrYdadjYd) -        & 
&  20._dp*(TrYuadjYu) - 20._dp*(TrYx3adjYx3)) + 810*g3p2*NGHg3 +          & 
&  5*(5._dp*(g1p2) + 9._dp*(g2p2) + 34._dp*(g3p2))*(NGHx3 +               & 
&  NGHxb3)))/15._dp

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yu 
!-------------------- 
 
betaYu1  = ((-26._dp*(g1p2) - 90._dp*(g2p2) - 160._dp*(g3p2) + 9._dp*(TrYb3adjYb3)    & 
&  + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3) + 90._dp*(TrYx3adjYx3))*Yu)               & 
& /30._dp + YuadjYdYd + 3._dp*(YuadjYuYu)

 
 
If (TwoLoopRGE) Then 
betaYu2 = (2._dp*(g1p2)/5._dp - 3._dp*(TrYdadjYd) - TrYeadjYe)*YuadjYdYd + ((4._dp*(g1p2) +     & 
&  60._dp*(g2p2) - 9._dp*(TrYb3adjYb3) - 90._dp*(TrYuadjYu) - 45._dp*(TrYw3adjYw3) -     & 
&  90._dp*(TrYx3adjYx3))*YuadjYuYu)/10._dp - 2*(YuadjYdTpYx3CYx3Yd + YuadjYdYdadjYdYd +  & 
&  YuadjYdYdadjYuYu + 2._dp*(YuadjYuYuadjYuYu)) + (Yu*(-486._dp*(TrYb3adjYb3Yb3adjYb3) - & 
&  540._dp*(TrYb3adjYeYeadjYb3) - 2565._dp*(TrYb3adjYw3Yw3adjYb3) + 2*(5486._dp*(g1p4) + & 
&  20*g1p2*(45._dp*(g2p2) + 136._dp*(g3p2)) + 50*(135._dp*(g2p4) + 144*g2p2*g3p2 -       & 
&  32._dp*(g3p4)) - 5400._dp*(TrYdadjYdTpYx3CYx3) - 2700._dp*(TrYdadjYuYuadjYd) -        & 
&  1350._dp*(TrYeadjYw3Yw3adjYe) - 8100._dp*(TrYuadjYuYuadjYu) - 3375._dp*(TrYw3adjYw3Yw3adjYw3) -& 
&  8100._dp*(TrYx3adjYx3Yx3adjYx3)) + 60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) +& 
&  65*g1p2*(NGHx3 + NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu +   & 
&  TrYx3adjYx3) + g3p2*(3*NGHg3 + NGHx3 + NGHxb3)) +& 
&  2700*g2p2*(4._dp*(TrYw3adjYw3) + g2p2*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))))/1800._dp

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
betaYd1  = 2._dp*(TpYx3CYx3Yd) + (-7._dp*(g1p2)/15._dp - 3._dp*(g2p2) -               & 
&  16._dp*(g3p2)/3._dp + 3._dp*(TrYdadjYd) + TrYeadjYe)*Yd + 3._dp*(YdadjYdYd)           & 
&  + YdadjYuYu

 
 
If (TwoLoopRGE) Then 
betaYd2 = (4._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 9._dp*(TrYdadjYd) - 3._dp*(TrYeadjYe))*YdadjYdYd +& 
&  (2*TpYx3CYx3Yd*(-3._dp*(TrYb3adjYb3) + 5*(2._dp*(g1p2) + 6._dp*(g2p2) -               & 
&  6._dp*(TrYuadjYu) - 3._dp*(TrYw3adjYw3) - 6._dp*(TrYx3adjYx3))) + (8._dp*(g1p2) -     & 
&  3._dp*(TrYb3adjYb3) - 30._dp*(TrYuadjYu) - 15._dp*(TrYw3adjYw3) - 30._dp*(TrYx3adjYx3))*YdadjYuYu)/10._dp -& 
&  2*(TpYx3CYx3TpYx3CYx3Yd + YdadjYdTpYx3CYx3Yd + 2._dp*(YdadjYdYdadjYdYd) +             & 
&  YdadjYuYuadjYdYd + YdadjYuYuadjYuYu) + (Yd*(287._dp*(g1p4) + 90*g1p2*g2p2 +           & 
&  675._dp*(g2p4) + 80*(g1p2 + 9._dp*(g2p2))*g3p2 - 160._dp*(g3p4) - 27._dp*(TrYb3adjYeYeadjYb3) -& 
&  540._dp*(TrYdadjYdTpYx3CYx3) - 810._dp*(TrYdadjYdYdadjYd) - 270._dp*(TrYdadjYuYuadjYd) -& 
&  270._dp*(TrYeadjYeYeadjYe) - 135._dp*(TrYeadjYw3Yw3adjYe) + 135*g2p4*(4*NGHw3 +& 
&  3*(NGHx3 + NGHxb3)) + 3*g1p2*(-12*(TrYdadjYd -          & 
&  3._dp*(TrYeadjYe)) + 35*g1p2*(NGHx3 + NGHxb3)) +        & 
&  480*g3p2*(3._dp*(TrYdadjYd) + g3p2*(3*NGHg3 + NGHx3 +   & 
&  NGHxb3))))/90._dp

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
!-------------------- 
! Ye 
!-------------------- 
 
betaYe1  = (-9._dp*(g1p2)/5._dp - 3._dp*(g2p2) + 3._dp*(TrYdadjYd) + TrYeadjYe)       & 
& *Ye + (3*(YeadjYb3Yb3 + 5*(2._dp*(YeadjYeYe) + YeadjYw3Yw3)))/10._dp

 
 
If (TwoLoopRGE) Then 
betaYe2 = (-36._dp*(YeadjYb3Yb3adjYb3Yb3) - 18*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) +            & 
&  TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))*(YeadjYb3Yb3 + 5._dp*(YeadjYw3Yw3)) -             & 
&  600*((-2._dp*(g2p2) + 3._dp*(TrYdadjYd) + TrYeadjYe)*YeadjYeYe - 2*g2p2*YeadjYw3Yw3) -& 
&  5*(24._dp*(YeadjYb3Yb3adjYeYe) + 12._dp*(YeadjYb3Yb3adjYw3Yw3) + 160._dp*(YeadjYeYeadjYeYe) +& 
&  9._dp*(YeadjYw3Yw3adjYb3Yb3) + 60*(2._dp*(YeadjYw3Yw3adjYeYe) + YeadjYw3Yw3adjYw3Yw3)) +& 
&  20*Ye*(-4*(g1p2 - 40._dp*(g3p2))*TrYdadjYd + 12*g1p2*TrYeadjYe + 3*(45._dp*(g1p4) +   & 
&  6*g1p2*g2p2 + 25._dp*(g2p4) - TrYb3adjYeYeadjYb3 - 20._dp*(TrYdadjYdTpYx3CYx3) -      & 
&  30._dp*(TrYdadjYdYdadjYd) - 10._dp*(TrYdadjYuYuadjYd) - 10._dp*(TrYeadjYeYeadjYe) -   & 
&  5._dp*(TrYeadjYw3Yw3adjYe)) + 45*g1p4*NGHx3 + 45*g1p4*NGHxb3 +& 
&  15*g2p4*(4*NGHw3 + 3*(NGHx3 + NGHxb3))))/200._dp

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
!-------------------- 
! Yb3 
!-------------------- 
 
betaYb31  = (3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) -              & 
&  2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))*Yb3)/10._dp + 9._dp*(Yb3adjYb3Yb3)     & 
& /10._dp + Yb3adjYeYe + 3._dp*(Yb3adjYw3Yw3)/2._dp

 
 
If (TwoLoopRGE) Then 
betaYb32 = (18*(4._dp*(g1p2) + 20._dp*(g2p2) - 3._dp*(TrYb3adjYb3) - 30._dp*(TrYuadjYu) -        & 
&  15._dp*(TrYw3adjYw3) - 30._dp*(TrYx3adjYx3))*Yb3adjYb3Yb3 - 72._dp*(Yb3adjYb3Yb3adjYb3Yb3) +& 
&  40*(6._dp*(g1p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*Yb3adjYeYe -               & 
&  30*(3._dp*(TrYb3adjYb3) + 5*(-8._dp*(g2p2) + 6._dp*(TrYuadjYu) + 3._dp*(TrYw3adjYw3) +& 
&  6._dp*(TrYx3adjYx3)))*Yb3adjYw3Yw3 - 5*(12._dp*(Yb3adjYb3Yb3adjYw3Yw3) +              & 
&  24._dp*(Yb3adjYeYeadjYb3Yb3) + 80._dp*(Yb3adjYeYeadjYeYe) + 45._dp*(Yb3adjYw3Yw3adjYb3Yb3) +& 
&  60._dp*(Yb3adjYw3Yw3adjYw3Yw3)) + Yb3*(3*(276._dp*(g1p4) + 120*g1p2*g2p2 +            & 
&  500._dp*(g2p4) - 18._dp*(TrYb3adjYb3Yb3adjYb3) - 20._dp*(TrYb3adjYeYeadjYb3) -        & 
&  95._dp*(TrYb3adjYw3Yw3adjYb3) - 400._dp*(TrYdadjYdTpYx3CYx3) - 200._dp*(TrYdadjYuYuadjYd) -& 
&  100._dp*(TrYeadjYw3Yw3adjYe) - 600._dp*(TrYuadjYuYuadjYu) - 250._dp*(TrYw3adjYw3Yw3adjYw3) -& 
&  600._dp*(TrYx3adjYx3Yx3adjYx3)) + 20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu +            & 
&  20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 + 15*g1p4*NGHx3 + 15*g1p4*NGHxb3) +& 
&  300*g2p2*(4._dp*(TrYw3adjYw3) + g2p2*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))))/200._dp

 
DYb3 = oo16pi2*( betaYb31 + oo16pi2 * betaYb32 ) 

 
Else 
DYb3 = oo16pi2* betaYb31 
End If 
 
 
!-------------------- 
! Yw3 
!-------------------- 
 
betaYw31  = (-3._dp*(g1p2)/5._dp - 7._dp*(g2p2) + 3._dp*(TrYb3adjYb3)/10._dp +        & 
&  3._dp*(TrYuadjYu) + 3._dp*(TrYw3adjYw3)/2._dp + 3._dp*(TrYx3adjYx3))*Yw3 +            & 
&  3._dp*(Yw3adjYb3Yb3)/10._dp + Yw3adjYeYe + 5._dp*(Yw3adjYw3Yw3)/2._dp

 
 
If (TwoLoopRGE) Then 
betaYw32 = (-36._dp*(Yw3adjYb3Yb3adjYb3Yb3) + 40*(6._dp*(g1p2) - 15._dp*(TrYdadjYd) -            & 
&  5._dp*(TrYeadjYe))*Yw3adjYeYe + 40*(3._dp*(g1p2) + 25._dp*(g2p2))*Yw3adjYw3Yw3 -      & 
&  6*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) + TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))*(3._dp*(Yw3adjYb3Yb3) +& 
&  25._dp*(Yw3adjYw3Yw3)) - 5*(24._dp*(Yw3adjYb3Yb3adjYw3Yw3) + 80._dp*(Yw3adjYeYeadjYeYe) +& 
&  40._dp*(Yw3adjYeYeadjYw3Yw3) + 9._dp*(Yw3adjYw3Yw3adjYb3Yb3) + 120._dp*(Yw3adjYw3Yw3adjYw3Yw3)) +& 
&  Yw3*(828._dp*(g1p4) + 360*g1p2*g2p2 + 5500._dp*(g2p4) - 54._dp*(TrYb3adjYb3Yb3adjYb3) -& 
&  60._dp*(TrYb3adjYeYeadjYb3) - 285._dp*(TrYb3adjYw3Yw3adjYb3) - 1200._dp*(TrYdadjYdTpYx3CYx3) -& 
&  600._dp*(TrYdadjYuYuadjYd) - 300._dp*(TrYeadjYw3Yw3adjYe) - 1800._dp*(TrYuadjYuYuadjYu) +& 
&  1200*g2p2*TrYw3adjYw3 - 750._dp*(TrYw3adjYw3Yw3adjYw3) - 1800._dp*(TrYx3adjYx3Yx3adjYx3) +& 
&  20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu + 20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 +       & 
&  15*g1p4*NGHx3 + 15*g1p4*NGHxb3) + 700*g2p4*(4*NGHw3 +& 
&  3*(NGHx3 + NGHxb3))))/200._dp

 
DYw3 = oo16pi2*( betaYw31 + oo16pi2 * betaYw32 ) 

 
Else 
DYw3 = oo16pi2* betaYw31 
End If 
 
 
!-------------------- 
! Yx3 
!-------------------- 
 
betaYx31  = ((-38._dp*(g1p2) - 90._dp*(g2p2) - 160._dp*(g3p2) + 9._dp*(TrYb3adjYb3)   & 
&  + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3) + 90._dp*(TrYx3adjYx3))*Yx3)              & 
& /30._dp + 3._dp*(Yx3adjYx3Yx3) + 2._dp*(Yx3CYdTpYd)

 
 
If (TwoLoopRGE) Then 
betaYx32 = (8._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 9._dp*(TrYb3adjYb3)/10._dp - 9._dp*(TrYuadjYu) - & 
&  9._dp*(TrYw3adjYw3)/2._dp - 9._dp*(TrYx3adjYx3))*Yx3adjYx3Yx3 + (2*(g1p2 +            & 
&  15._dp*(g2p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*Yx3CYdTpYd)/5._dp -           & 
&  2*(2._dp*(Yx3adjYx3Yx3adjYx3Yx3) + Yx3CYdTpYdadjYx3Yx3 + Yx3CYdTpYdCYdTpYd +          & 
&  Yx3CYdTpYuCYuTpYd) + (Yx3*(-486._dp*(TrYb3adjYb3Yb3adjYb3) - 540._dp*(TrYb3adjYeYeadjYb3) -& 
&  2565._dp*(TrYb3adjYw3Yw3adjYb3) + 2*(8246._dp*(g1p4) + 20*g1p2*(153._dp*(g2p2) +      & 
&  232._dp*(g3p2)) + 50*(135._dp*(g2p4) + 144*g2p2*g3p2 - 32._dp*(g3p4)) -               & 
&  5400._dp*(TrYdadjYdTpYx3CYx3) - 2700._dp*(TrYdadjYuYuadjYd) - 1350._dp*(TrYeadjYw3Yw3adjYe) -& 
&  8100._dp*(TrYuadjYuYuadjYu) - 3375._dp*(TrYw3adjYw3Yw3adjYw3) - 8100._dp*(TrYx3adjYx3Yx3adjYx3)) +& 
&  60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) + 95*g1p2*(NGHx3 +  & 
&  NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu + TrYx3adjYx3) + g3p2*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)) + 2700*g2p2*(4._dp*(TrYw3adjYw3) +     & 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))))/1800._dp

 
DYx3 = oo16pi2*( betaYx31 + oo16pi2 * betaYx32 ) 

 
Else 
DYx3 = oo16pi2* betaYx31 
End If 
 
 
!-------------------- 
! mue 
!-------------------- 
 
betamue1  = ((3._dp*(TrYb3adjYb3) + 30._dp*(TrYdadjYd) + 10._dp*(TrYeadjYe)           & 
&  + 30._dp*(TrYuadjYu) + 15._dp*(TrYw3adjYw3) - 6*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))& 
& *mue)/10._dp

 
 
If (TwoLoopRGE) Then 
betamue2 = 6*g2p4*mue*NGHw3 + (mue*(80*(-((g1p2 - 40._dp*(g3p2))*TrYdadjYd) +     & 
&  3*g1p2*TrYeadjYe + 2*(g1p2 + 20._dp*(g3p2))*TrYuadjYu + 5*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3) +& 
&  3*(276._dp*(g1p4) + 120*g1p2*g2p2 + 500._dp*(g2p4) - 18._dp*(TrYb3adjYb3Yb3adjYb3) -  & 
&  40._dp*(TrYb3adjYeYeadjYb3) - 95._dp*(TrYb3adjYw3Yw3adjYb3) - 800._dp*(TrYdadjYdTpYx3CYx3) -& 
&  600._dp*(TrYdadjYdYdadjYd) - 400._dp*(TrYdadjYuYuadjYd) - 200._dp*(TrYeadjYeYeadjYe) -& 
&  200._dp*(TrYeadjYw3Yw3adjYe) - 600._dp*(TrYuadjYuYuadjYu) + 400*g2p2*TrYw3adjYw3 -    & 
&  250._dp*(TrYw3adjYw3Yw3adjYw3) - 600._dp*(TrYx3adjYx3Yx3adjYx3)) + 300*(g1p4 +        & 
&  3._dp*(g2p4))*NGHx3 + 300*(g1p4 + 3._dp*(g2p4))*NGHxb3))/200._dp

 
Dmue = oo16pi2*( betamue1 + oo16pi2 * betamue2 ) 

 
Else 
Dmue = oo16pi2* betamue1 
End If 
 
 
!-------------------- 
! MXM3 
!-------------------- 
 
betaMXM31  = -((5._dp*(g1p2) + 9._dp*(g2p2) + 16._dp*(g3p2))*MXM3)/3._dp +            & 
&  MXM3CYx3TpYx3

 
 
If (TwoLoopRGE) Then 
betaMXM32 = (-9*(20*(MXM3CYx3TpYx3CYx3TpYx3 + MXM3CYx3YdadjYdTpYx3) + MXM3CYx3TpYx3*(4._dp*(g1p2) +& 
&  3._dp*(TrYb3adjYb3) + 30._dp*(TrYuadjYu) + 15._dp*(TrYw3adjYw3) + 30._dp*(TrYx3adjYx3))) +& 
&  5*MXM3*(223._dp*(g1p4) + 90*g1p2*g2p2 + 135._dp*(g2p4) + 32*(5._dp*(g1p2) +           & 
&  9._dp*(g2p2))*g3p2 - 32._dp*(g3p4) + 108*g2p4*NGHw3 + 3*(25._dp*(g1p4) +& 
&  27._dp*(g2p4))*(NGHx3 + NGHxb3) + 96*g3p4*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)))/90._dp

 
DMXM3 = oo16pi2*( betaMXM31 + oo16pi2 * betaMXM32 ) 

 
Else 
DMXM3 = oo16pi2* betaMXM31 
End If 
 
 
!-------------------- 
! MWM3 
!-------------------- 
 
betaMWM31  = -8*g2p2*MWM3 + MWM3CYw3TpYw3 + Yw3adjYw3MWM3

 
 
If (TwoLoopRGE) Then 
betaMWM32 = (-3._dp*(MWM3CYw3TpYb3CYb3TpYw3) - 10._dp*(MWM3CYw3TpYeCYeTpYw3) - 15._dp*(MWM3CYw3TpYw3CYw3TpYw3) -& 
&  3._dp*(Yw3adjYb3Yb3adjYw3MWM3) - 10._dp*(Yw3adjYeYeadjYw3MWM3) + (6._dp*(g1p2) -      & 
&  10._dp*(g2p2) - 3._dp*(TrYb3adjYb3) - 30._dp*(TrYuadjYu) - 15._dp*(TrYw3adjYw3) -     & 
&  30._dp*(TrYx3adjYx3))*(MWM3CYw3TpYw3 + Yw3adjYw3MWM3) - 15._dp*(Yw3adjYw3Yw3adjYw3MWM3) +& 
&  40*g2p4*MWM3*(10 + 4*NGHw3 + 3*NGHx3 + 3*NGHxb3))/10._dp

 
DMWM3 = oo16pi2*( betaMWM31 + oo16pi2 * betaMWM32 ) 

 
Else 
DMWM3 = oo16pi2* betaMWM31 
End If 
 
 
!-------------------- 
! MGM3 
!-------------------- 
 
betaMGM31  = -12*g3p2*MGM3

 
 
If (TwoLoopRGE) Then 
betaMGM32 = 12*g3p4*MGM3*(3 + 3*NGHg3 + NGHx3 + NGHxb3)

 
DMGM3 = oo16pi2*( betaMGM31 + oo16pi2 * betaMGM32 ) 

 
Else 
DMGM3 = oo16pi2* betaMGM31 
End If 
 
 
!-------------------- 
! MBM3 
!-------------------- 
 
betaMBM31  = (3*(MBM3CYb3TpYb3 + Yb3adjYb3MBM3))/5._dp

 
 
If (TwoLoopRGE) Then 
betaMBM32 = (-3*(3._dp*(MBM3CYb3TpYb3CYb3TpYb3) + 10._dp*(MBM3CYb3TpYeCYeTpYb3) + 15._dp*(MBM3CYb3TpYw3CYw3TpYb3) +& 
&  3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) -  & 
&  5._dp*(TrYx3adjYx3)))*(MBM3CYb3TpYb3 + Yb3adjYb3MBM3) + 3._dp*(Yb3adjYb3Yb3adjYb3MBM3) +& 
&  10._dp*(Yb3adjYeYeadjYb3MBM3) + 15._dp*(Yb3adjYw3Yw3adjYb3MBM3)))/50._dp

 
DMBM3 = oo16pi2*( betaMBM31 + oo16pi2 * betaMBM32 ) 

 
Else 
DMBM3 = oo16pi2* betaMBM31 
End If 
 
 
!-------------------- 
! TYu 
!-------------------- 
 
betaTYu1  = TYuadjYdYd + 5._dp*(TYuadjYuYu) + ((26*g1p2*MassB)/15._dp +               & 
&  (32*g3p2*MassG)/3._dp + 6*g2p2*MassWB + 3._dp*(TradjYb3TYb3)/5._dp + 6._dp*(TradjYuTYu)& 
&  + 3._dp*(TradjYw3TYw3) + 6._dp*(TradjYx3TYx3))*Yu + 2._dp*(YuadjYdTYd) +              & 
&  4._dp*(YuadjYuTYu) + ((-26._dp*(g1p2) - 90._dp*(g2p2) - 160._dp*(g3p2) +              & 
&  9._dp*(TrYb3adjYb3) + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3) + 90._dp*(TrYx3adjYx3))& 
& *TYu)/30._dp

 
 
If (TwoLoopRGE) Then 
betaTYu2 = (2700*(8._dp*(g2p2) - TrYb3adjYb3 - 10._dp*(TrYuadjYu) - 5._dp*(TrYw3adjYw3) -        & 
&  10._dp*(TrYx3adjYx3))*TYuadjYuYu + 360*((2._dp*(g1p2) - 15._dp*(TrYdadjYd) -          & 
&  5._dp*(TrYeadjYe))*TYuadjYdYd + 2*(2._dp*(g1p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*YuadjYdTYd -& 
&  2*(2*g1p2*MassB + 15._dp*(TradjYdTYd) + 5._dp*(TradjYeTYe))*YuadjYdYd +               & 
&  6*(g1p2 + 5._dp*(g2p2) - TrYb3adjYb3 - 10._dp*(TrYuadjYu) - 5._dp*(TrYw3adjYw3) -     & 
&  10._dp*(TrYx3adjYx3))*YuadjYuTYu - (4*g1p2*MassB + 60*g2p2*MassWB + 9._dp*(TradjYb3TYb3) +& 
&  90._dp*(TradjYuTYu) + 45._dp*(TradjYw3TYw3) + 90._dp*(TradjYx3TYx3))*YuadjYuYu) -     & 
&  3600*(TYuadjYdTpYx3CYx3Yd + TYuadjYdYdadjYdYd + 2._dp*(TYuadjYdYdadjYuYu) +           & 
&  3._dp*(TYuadjYuYuadjYuYu) + 2._dp*(YuadjYdTpTYx3CYx3Yd) + 2._dp*(YuadjYdTpYx3CYx3TYd) +& 
&  2._dp*(YuadjYdTYdadjYdYd) + 2._dp*(YuadjYdTYdadjYuYu) + 2._dp*(YuadjYdYdadjYdTYd) +   & 
&  YuadjYdYdadjYuTYu + 4._dp*(YuadjYuTYuadjYuYu) + 3._dp*(YuadjYuYuadjYuTYu)) -          & 
&  4*Yu*(486._dp*(TrYb3adjYb3TYb3adjYb3) + 270._dp*(TrYb3adjYeTYeadjYb3) +               & 
&  1215._dp*(TrYb3adjYw3TYw3adjYb3) + 2700._dp*(TrYdadjYuTYuadjYd) + 270._dp*(TrYeadjYb3TYb3adjYe) +& 
&  1350._dp*(TrYeadjYw3TYw3adjYe) + 2700._dp*(TrYuadjYdTYdadjYu) + 16200._dp*(TrYuadjYuTYuadjYu) +& 
&  1215._dp*(TrYw3adjYb3TYb3adjYw3) + 1350._dp*(TrYw3adjYeTYeadjYw3) + 6075._dp*(TrYw3adjYw3TYw3adjYw3) +& 
&  4*(2743*g1p4*MassB + 5*g1p2*(136*g3p2*(MassB + MassG) + 45*g2p2*(MassB +              & 
&  MassWB)) + 25*(-32*g3p4*MassG + 135*g2p4*MassWB + 72*g2p2*g3p2*(MassG +               & 
&  MassWB)) + 1350._dp*(TradjYdTpYx3CYx3TYd) + 1350._dp*(TradjYx3TYx3CYdTpYd) +          & 
&  4050._dp*(TrYx3adjYx3TYx3adjYx3)) + 28800*g3p4*MassG*NGHg3 +           & 
&  60*(6*(-2*(g1p2 + 20._dp*(g3p2))*TradjYuTYu - 5*(g1p2 + 8._dp*(g3p2))*TradjYx3TYx3 +  & 
&  2*(g1p2*MassB + 20*g3p2*MassG)*TrYuadjYu + 5*(g1p2*MassB + 8*g3p2*MassG)*TrYx3adjYx3) +& 
&  5*(13*g1p4*MassB + 32*g3p4*MassG)*NGHx3 + 5*(13*g1p4*MassB +           & 
&  32*g3p4*MassG)*NGHxb3) + 2700*g2p2*(-2._dp*(TradjYw3TYw3) +            & 
&  2*MassWB*TrYw3adjYw3 + g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))) + (-486._dp*(TrYb3adjYb3Yb3adjYb3) - 540._dp*(TrYb3adjYeYeadjYb3) -& 
&  2430._dp*(TrYb3adjYw3Yw3adjYb3) - 10800._dp*(TrYdadjYdTpYx3CYx3) - 5400._dp*(TrYdadjYuYuadjYd) -& 
&  2700._dp*(TrYeadjYw3Yw3adjYe) - 16200._dp*(TrYuadjYuYuadjYu) - 6075._dp*(TrYw3adjYw3Yw3adjYw3) +& 
&  4*(2743._dp*(g1p4) + 450*g1p2*g2p2 + 3375._dp*(g2p4) + 80*(17._dp*(g1p2) +            & 
&  45._dp*(g2p2))*g3p2 - 800._dp*(g3p4) - 4050._dp*(TrYx3adjYx3Yx3adjYx3)) +             & 
&  60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) + 65*g1p2*(NGHx3 +  & 
&  NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu + TrYx3adjYx3) + g3p2*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)) + 2700*g2p2*(4._dp*(TrYw3adjYw3) +     & 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3))))*TYu)/1800._dp

 
DTYu = oo16pi2*( betaTYu1 + oo16pi2 * betaTYu2 ) 

 
Else 
DTYu = oo16pi2* betaTYu1 
End If 
 
 
!-------------------- 
! TYd 
!-------------------- 
 
betaTYd1  = 4._dp*(TpTYx3CYx3Yd) + 2._dp*(TpYx3CYx3TYd) + 5._dp*(TYdadjYdYd)          & 
&  + TYdadjYuYu + ((14*g1p2*MassB)/15._dp + (32*g3p2*MassG)/3._dp + 6*g2p2*MassWB +      & 
&  6._dp*(TradjYdTYd) + 2._dp*(TradjYeTYe))*Yd + 4._dp*(YdadjYdTYd) + 2._dp*(YdadjYuTYu) & 
&  - ((7._dp*(g1p2) + 45._dp*(g2p2) + 80._dp*(g3p2) - 45._dp*(TrYdadjYd) -               & 
&  15._dp*(TrYeadjYe))*TYd)/15._dp

 
 
If (TwoLoopRGE) Then 
betaTYd2 = 4*g1p2*TpTYx3CYx3Yd + 12*g2p2*TpTYx3CYx3Yd + 2*g1p2*TpYx3CYx3TYd + 6*g2p2*TpYx3CYx3TYd -& 
&  4*g1p2*MassB*TpYx3CYx3Yd - 12*g2p2*MassWB*TpYx3CYx3Yd - (6*TpYx3CYx3Yd*TradjYb3TYb3)/5._dp -& 
&  12*TpYx3CYx3Yd*TradjYuTYu - 6*TpYx3CYx3Yd*TradjYw3TYw3 - 12*TpYx3CYx3Yd*TradjYx3TYx3 -& 
&  (6*TpTYx3CYx3Yd*TrYb3adjYb3)/5._dp - (3*TpYx3CYx3TYd*TrYb3adjYb3)/5._dp -             & 
&  12*TpTYx3CYx3Yd*TrYuadjYu - 6*TpYx3CYx3TYd*TrYuadjYu - 6*TpTYx3CYx3Yd*TrYw3adjYw3 -   & 
&  3*TpYx3CYx3TYd*TrYw3adjYw3 - 12*TpTYx3CYx3Yd*TrYx3adjYx3 - 6*TpYx3CYx3TYd*TrYx3adjYx3 +& 
&  (6*g1p2*TYdadjYdYd)/5._dp + 12*g2p2*TYdadjYdYd - 15*TrYdadjYd*TYdadjYdYd -            & 
&  5*TrYeadjYe*TYdadjYdYd + (4*g1p2*TYdadjYuYu)/5._dp - (3*TrYb3adjYb3*TYdadjYuYu)/10._dp -& 
&  3*TrYuadjYu*TYdadjYuYu - (3*TrYw3adjYw3*TYdadjYuYu)/2._dp - 3*TrYx3adjYx3*TYdadjYuYu +& 
&  (6*g1p2*YdadjYdTYd)/5._dp + 6*g2p2*YdadjYdTYd - 12*TrYdadjYd*YdadjYdTYd -             & 
&  4*TrYeadjYe*YdadjYdTYd - (2*(4*g1p2*MassB + 30*g2p2*MassWB + 45._dp*(TradjYdTYd) +    & 
&  15._dp*(TradjYeTYe))*YdadjYdYd)/5._dp + (8*g1p2*YdadjYuTYu)/5._dp - (3*TrYb3adjYb3*YdadjYuTYu)/5._dp -& 
&  6*TrYuadjYu*YdadjYuTYu - 3*TrYw3adjYw3*YdadjYuTYu - 6*TrYx3adjYx3*YdadjYuTYu -        & 
&  (8*g1p2*MassB*YdadjYuYu)/5._dp - (3*TradjYb3TYb3*YdadjYuYu)/5._dp - 6*TradjYuTYu*YdadjYuYu -& 
&  3*TradjYw3TYw3*YdadjYuYu - 6*TradjYx3TYx3*YdadjYuYu - 2*(2*(TpTYx3CYx3TpYx3CYx3Yd +   & 
&  TpYx3CYx3TpTYx3CYx3Yd) + TpYx3CYx3TpYx3CYx3TYd + TYdadjYdTpYx3CYx3Yd + 3._dp*(TYdadjYdYdadjYdYd) +& 
&  2._dp*(TYdadjYuYuadjYdYd) + TYdadjYuYuadjYuYu + 2._dp*(YdadjYdTpTYx3CYx3Yd) +         & 
&  2._dp*(YdadjYdTpYx3CYx3TYd) + 4._dp*(YdadjYdTYdadjYdYd) + 3._dp*(YdadjYdYdadjYdTYd) + & 
&  2._dp*(YdadjYuTYuadjYdYd) + 2._dp*(YdadjYuTYuadjYuYu) + YdadjYuYuadjYdTYd +           & 
&  2._dp*(YdadjYuYuadjYuTYu)) - (Yd*(2*(287*g1p4*MassB + 5*g1p2*(8*g3p2*(MassB +         & 
&  MassG) + 9*g2p2*(MassB + MassWB)) + 5*(-32*g3p4*MassG + 135*g2p4*MassWB +             & 
&  72*g2p2*g3p2*(MassG + MassWB)) + 270._dp*(TradjYdTpYx3CYx3TYd) + 270._dp*(TradjYx3TYx3CYdTpYd)) +& 
&  27._dp*(TrYb3adjYeTYeadjYb3) + 1620._dp*(TrYdadjYdTYdadjYd) + 270._dp*(TrYdadjYuTYuadjYd) +& 
&  27._dp*(TrYeadjYb3TYb3adjYe) + 540._dp*(TrYeadjYeTYeadjYe) + 135._dp*(TrYeadjYw3TYw3adjYe) +& 
&  270._dp*(TrYuadjYdTYdadjYu) + 135._dp*(TrYw3adjYeTYeadjYw3) + 1080*g2p4*MassWB*NGHw3 +& 
&  30*(-48*g3p2*TradjYdTYd + 48*g3p2*MassG*TrYdadjYd + 96*g3p4*MassG*NGHg3 +& 
&  32*g3p4*MassG*NGHx3 + 27*g2p4*MassWB*NGHx3 +            & 
&  32*g3p4*MassG*NGHxb3 + 27*g2p4*MassWB*NGHxb3) +         & 
&  6*g1p2*(6*(TradjYdTYd - 3._dp*(TradjYeTYe)) + MassB*(-6*(TrYdadjYd - 3._dp*(TrYeadjYe)) +& 
&  35*g1p2*(NGHx3 + NGHxb3)))))/45._dp + ((287._dp*(g1p4) +& 
&  90*g1p2*g2p2 + 675._dp*(g2p4) + 80*(g1p2 + 9._dp*(g2p2))*g3p2 - 160._dp*(g3p4) -      & 
&  27._dp*(TrYb3adjYeYeadjYb3) - 540._dp*(TrYdadjYdTpYx3CYx3) - 810._dp*(TrYdadjYdYdadjYd) -& 
&  270._dp*(TrYdadjYuYuadjYd) - 270._dp*(TrYeadjYeYeadjYe) - 135._dp*(TrYeadjYw3Yw3adjYe) +& 
&  135*g2p4*(4*NGHw3 + 3*(NGHx3 + NGHxb3)) +& 
&  3*g1p2*(-12*(TrYdadjYd - 3._dp*(TrYeadjYe)) + 35*g1p2*(NGHx3 +         & 
&  NGHxb3)) + 480*g3p2*(3._dp*(TrYdadjYd) + g3p2*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)))*TYd)/90._dp

 
DTYd = oo16pi2*( betaTYd1 + oo16pi2 * betaTYd2 ) 

 
Else 
DTYd = oo16pi2* betaTYd1 
End If 
 
 
!-------------------- 
! TYe 
!-------------------- 
 
betaTYe1  = 3._dp*(TYeadjYb3Yb3)/10._dp + 5._dp*(TYeadjYeYe) + 3._dp*(TYeadjYw3Yw3)   & 
& /2._dp + ((18*g1p2*MassB)/5._dp + 6*g2p2*MassWB + 6._dp*(TradjYdTYd) + 2._dp*(TradjYeTYe))& 
& *Ye + 3._dp*(YeadjYb3TYb3)/5._dp + 4._dp*(YeadjYeTYe) + 3._dp*(YeadjYw3TYw3)           & 
&  + ((-9._dp*(g1p2) - 15._dp*(g2p2) + 15._dp*(TrYdadjYd) + 5._dp*(TrYeadjYe))           & 
& *TYe)/5._dp

 
 
If (TwoLoopRGE) Then 
betaTYe2 = (-36._dp*(TYeadjYb3Yb3adjYb3Yb3) - 240._dp*(TYeadjYb3Yb3adjYeYe) - 45._dp*(TYeadjYb3Yb3adjYw3Yw3) -& 
&  1200._dp*(TYeadjYeYeadjYeYe) - 45._dp*(TYeadjYw3Yw3adjYb3Yb3) - 1200._dp*(TYeadjYw3Yw3adjYeYe) -& 
&  225._dp*(TYeadjYw3Yw3adjYw3Yw3) - 72._dp*(YeadjYb3TYb3adjYb3Yb3) - 240._dp*(YeadjYb3TYb3adjYeYe) -& 
&  90._dp*(YeadjYb3TYb3adjYw3Yw3) - 72._dp*(YeadjYb3Yb3adjYb3TYb3) - 120._dp*(YeadjYb3Yb3adjYeTYe) -& 
&  90._dp*(YeadjYb3Yb3adjYw3TYw3) - 1600._dp*(YeadjYeTYeadjYeYe) - 1200*(2*g2p2*MassWB + & 
&  3._dp*(TradjYdTYd) + TradjYeTYe)*YeadjYeYe - 1200._dp*(YeadjYeYeadjYeTYe) -           & 
&  18*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) + TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))*(TYeadjYb3Yb3 +& 
&  5._dp*(TYeadjYw3Yw3) + 2._dp*(YeadjYb3TYb3) + 10._dp*(YeadjYw3TYw3)) - 90._dp*(YeadjYw3TYw3adjYb3Yb3) -& 
&  1200._dp*(YeadjYw3TYw3adjYeYe) - 450._dp*(YeadjYw3TYw3adjYw3Yw3) - 36*(TradjYb3TYb3 + & 
&  5*(2._dp*(TradjYuTYu) + TradjYw3TYw3 + 2._dp*(TradjYx3TYx3)))*(YeadjYb3Yb3 +          & 
&  5._dp*(YeadjYw3Yw3)) + 40*(-6*g1p2*TYeadjYeYe + 60*g2p2*TYeadjYeYe - 75*TrYdadjYd*TYeadjYeYe -& 
&  25*TrYeadjYe*TYeadjYeYe + 30*g2p2*TYeadjYw3Yw3 + (6._dp*(g1p2) + 30._dp*(g2p2) -      & 
&  60._dp*(TrYdadjYd) - 20._dp*(TrYeadjYe))*YeadjYeTYe + 60*g2p2*YeadjYw3TYw3 -          & 
&  60*g2p2*MassWB*YeadjYw3Yw3) - 90._dp*(YeadjYw3Yw3adjYb3TYb3) - 600._dp*(YeadjYw3Yw3adjYeTYe) -& 
&  450._dp*(YeadjYw3Yw3adjYw3TYw3) - 40*Ye*(-4*(-((g1p2 - 40._dp*(g3p2))*TradjYdTYd) +   & 
&  3*g1p2*TradjYeTYe + (g1p2*MassB - 40*g3p2*MassG)*TrYdadjYd - 3*g1p2*MassB*TrYeadjYe) +& 
&  3*(90*g1p4*MassB + 50*g2p4*MassWB + 6*g1p2*g2p2*(MassB + MassWB) + 20._dp*(TradjYdTpYx3CYx3TYd) +& 
&  20._dp*(TradjYx3TYx3CYdTpYd) + TrYb3adjYeTYeadjYb3 + 60._dp*(TrYdadjYdTYdadjYd) +     & 
&  10._dp*(TrYdadjYuTYuadjYd) + TrYeadjYb3TYb3adjYe + 20._dp*(TrYeadjYeTYeadjYe) +       & 
&  5._dp*(TrYeadjYw3TYw3adjYe) + 10._dp*(TrYuadjYdTYdadjYu) + 5._dp*(TrYw3adjYeTYeadjYw3)) +& 
&  90*g1p4*MassB*NGHx3 + 90*g1p4*MassB*NGHxb3 +            & 
&  30*g2p4*MassWB*(4*NGHw3 + 3*(NGHx3 + NGHxb3))) +& 
&  20*(-4*(g1p2 - 40._dp*(g3p2))*TrYdadjYd + 12*g1p2*TrYeadjYe + 3*(45._dp*(g1p4) +      & 
&  6*g1p2*g2p2 + 25._dp*(g2p4) - TrYb3adjYeYeadjYb3 - 20._dp*(TrYdadjYdTpYx3CYx3) -      & 
&  30._dp*(TrYdadjYdYdadjYd) - 10._dp*(TrYdadjYuYuadjYd) - 10._dp*(TrYeadjYeYeadjYe) -   & 
&  5._dp*(TrYeadjYw3Yw3adjYe)) + 45*g1p4*NGHx3 + 45*g1p4*NGHxb3 +& 
&  15*g2p4*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))*TYe)/200._dp

 
DTYe = oo16pi2*( betaTYe1 + oo16pi2 * betaTYe2 ) 

 
Else 
DTYe = oo16pi2* betaTYe1 
End If 
 
 
!-------------------- 
! TYb3 
!-------------------- 
 
betaTYb31  = (12*(2*g1p2*MassB + 10*g2p2*MassWB + TradjYb3TYb3 + 10._dp*(TradjYuTYu)  & 
&  + 5._dp*(TradjYw3TYw3) + 10._dp*(TradjYx3TYx3))*Yb3 + 24._dp*(Yb3adjYb3TYb3)          & 
&  + 5*(6._dp*(TYb3adjYb3Yb3) + 4._dp*(TYb3adjYeYe) + 12._dp*(TYb3adjYw3Yw3)             & 
&  + 8._dp*(Yb3adjYeTYe) + 15._dp*(Yb3adjYw3TYw3)) + 6*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) & 
&  + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))*TYb3)/20._dp

 
 
If (TwoLoopRGE) Then 
betaTYb32 = (144*g1p2*TYb3adjYb3Yb3 + 720*g2p2*TYb3adjYb3Yb3 - 90*TrYb3adjYb3*TYb3adjYb3Yb3 -     & 
&  900*TrYuadjYu*TYb3adjYb3Yb3 - 450*TrYw3adjYw3*TYb3adjYb3Yb3 - 900*TrYx3adjYx3*TYb3adjYb3Yb3 -& 
&  108._dp*(TYb3adjYb3Yb3adjYb3Yb3) - 135._dp*(TYb3adjYb3Yb3adjYw3Yw3) + 240*g1p2*TYb3adjYeYe -& 
&  600*TrYdadjYd*TYb3adjYeYe - 200*TrYeadjYe*TYb3adjYeYe - 240._dp*(TYb3adjYeYeadjYb3Yb3) -& 
&  400._dp*(TYb3adjYeYeadjYeYe) - 300._dp*(TYb3adjYeYeadjYw3Yw3) + 180*g1p2*TYb3adjYw3Yw3 +& 
&  2100*g2p2*TYb3adjYw3Yw3 - 180*TrYb3adjYb3*TYb3adjYw3Yw3 - 1800*TrYuadjYu*TYb3adjYw3Yw3 -& 
&  900*TrYw3adjYw3*TYb3adjYw3Yw3 - 1800*TrYx3adjYx3*TYb3adjYw3Yw3 - 405._dp*(TYb3adjYw3Yw3adjYb3Yb3) -& 
&  675._dp*(TYb3adjYw3Yw3adjYw3Yw3) + 72*g1p2*Yb3adjYb3TYb3 + 360*g2p2*Yb3adjYb3TYb3 -   & 
&  72*TrYb3adjYb3*Yb3adjYb3TYb3 - 720*TrYuadjYu*Yb3adjYb3TYb3 - 360*TrYw3adjYw3*Yb3adjYb3TYb3 -& 
&  720*TrYx3adjYx3*Yb3adjYb3TYb3 - 144._dp*(Yb3adjYb3TYb3adjYb3Yb3) - 180._dp*(Yb3adjYb3TYb3adjYw3Yw3) -& 
&  36*(4*g1p2*MassB + 20*g2p2*MassWB + 3._dp*(TradjYb3TYb3) + 30._dp*(TradjYuTYu) +      & 
&  15._dp*(TradjYw3TYw3) + 30._dp*(TradjYx3TYx3))*Yb3adjYb3Yb3 - 108._dp*(Yb3adjYb3Yb3adjYb3TYb3) -& 
&  135._dp*(Yb3adjYb3Yb3adjYw3TYw3) + 480*g1p2*Yb3adjYeTYe - 1200*TrYdadjYd*Yb3adjYeTYe -& 
&  400*TrYeadjYe*Yb3adjYeTYe - 240._dp*(Yb3adjYeTYeadjYb3Yb3) - 800._dp*(Yb3adjYeTYeadjYeYe) -& 
&  300._dp*(Yb3adjYeTYeadjYw3Yw3) - 480*g1p2*MassB*Yb3adjYeYe - 1200*TradjYdTYd*Yb3adjYeYe -& 
&  400*TradjYeTYe*Yb3adjYeYe - 120._dp*(Yb3adjYeYeadjYb3TYb3) - 800._dp*(Yb3adjYeYeadjYeTYe) -& 
&  150._dp*(Yb3adjYeYeadjYw3TYw3) + 90*g1p2*Yb3adjYw3TYw3 + 2850*g2p2*Yb3adjYw3TYw3 -    & 
&  225*TrYb3adjYb3*Yb3adjYw3TYw3 - 2250*TrYuadjYu*Yb3adjYw3TYw3 - 1125*TrYw3adjYw3*Yb3adjYw3TYw3 -& 
&  2250*TrYx3adjYx3*Yb3adjYw3TYw3 - 450._dp*(Yb3adjYw3TYw3adjYb3Yb3) - 900._dp*(Yb3adjYw3TYw3adjYw3Yw3) -& 
&  180*g1p2*MassB*Yb3adjYw3Yw3 - 3300*g2p2*MassWB*Yb3adjYw3Yw3 - 270*TradjYb3TYb3*Yb3adjYw3Yw3 -& 
&  2700*TradjYuTYu*Yb3adjYw3Yw3 - 1350*TradjYw3TYw3*Yb3adjYw3Yw3 - 2700*TradjYx3TYx3*Yb3adjYw3Yw3 -& 
&  270._dp*(Yb3adjYw3Yw3adjYb3TYb3) - 675._dp*(Yb3adjYw3Yw3adjYw3TYw3) - 4*Yb3*(40*(-    & 
& 2*(g1p2 + 20._dp*(g3p2))*TradjYuTYu - 5*(g1p2 + 8._dp*(g3p2))*TradjYx3TYx3 +           & 
&  2*(g1p2*MassB + 20*g3p2*MassG)*TrYuadjYu + 5*(g1p2*MassB + 8*g3p2*MassG)*TrYx3adjYx3) +& 
&  3*(18._dp*(TrYb3adjYb3TYb3adjYb3) + 10._dp*(TrYb3adjYeTYeadjYb3) + 45._dp*(TrYb3adjYw3TYw3adjYb3) +& 
&  100._dp*(TrYdadjYuTYuadjYd) + 10._dp*(TrYeadjYb3TYb3adjYe) + 50._dp*(TrYeadjYw3TYw3adjYe) +& 
&  100._dp*(TrYuadjYdTYdadjYu) + 600._dp*(TrYuadjYuTYuadjYu) + 45._dp*(TrYw3adjYb3TYb3adjYw3) +& 
&  50._dp*(TrYw3adjYeTYeadjYw3) + 225._dp*(TrYw3adjYw3TYw3adjYw3) + 4*(69*g1p4*MassB +   & 
&  125*g2p4*MassWB + 15*g1p2*g2p2*(MassB + MassWB) + 50._dp*(TradjYdTpYx3CYx3TYd) +      & 
&  50._dp*(TradjYx3TYx3CYdTpYd) + 150._dp*(TrYx3adjYx3TYx3adjYx3))) + 300*g1p4*MassB*NGHx3 +& 
&  300*g1p4*MassB*NGHxb3 + 300*g2p2*(-2._dp*(TradjYw3TYw3) +              & 
&  2*MassWB*TrYw3adjYw3 + g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))) + (3*(276._dp*(g1p4) + 120*g1p2*g2p2 + 500._dp*(g2p4) -     & 
&  18._dp*(TrYb3adjYb3Yb3adjYb3) - 20._dp*(TrYb3adjYeYeadjYb3) - 90._dp*(TrYb3adjYw3Yw3adjYb3) -& 
&  400._dp*(TrYdadjYdTpYx3CYx3) - 200._dp*(TrYdadjYuYuadjYd) - 100._dp*(TrYeadjYw3Yw3adjYe) -& 
&  600._dp*(TrYuadjYuYuadjYu) - 225._dp*(TrYw3adjYw3Yw3adjYw3) - 600._dp*(TrYx3adjYx3Yx3adjYx3)) +& 
&  20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu + 20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 +       & 
&  15*g1p4*NGHx3 + 15*g1p4*NGHxb3) + 300*g2p2*(4._dp*(TrYw3adjYw3) +& 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3))))*TYb3)/200._dp

 
DTYb3 = oo16pi2*( betaTYb31 + oo16pi2 * betaTYb32 ) 

 
Else 
DTYb3 = oo16pi2* betaTYb31 
End If 
 
 
!-------------------- 
! TYw3 
!-------------------- 
 
betaTYw31  = 9._dp*(TYw3adjYb3Yb3)/10._dp + TYw3adjYeYe + 3._dp*(TYw3adjYw3Yw3)       & 
&  + ((6*g1p2*MassB)/5._dp + 14*g2p2*MassWB + 3._dp*(TradjYb3TYb3)/5._dp +               & 
&  6._dp*(TradjYuTYu) + 3._dp*(TradjYw3TYw3) + 6._dp*(TradjYx3TYx3))*Yw3 +               & 
&  9._dp*(Yw3adjYb3TYb3)/10._dp + 2._dp*(Yw3adjYeTYe) + 15._dp*(Yw3adjYw3TYw3)           & 
& /4._dp + ((-6._dp*(g1p2) - 70._dp*(g2p2) + 3._dp*(TrYb3adjYb3) + 30._dp*(TrYuadjYu)    & 
&  + 15._dp*(TrYw3adjYw3) + 30._dp*(TrYx3adjYx3))*TYw3)/10._dp

 
 
If (TwoLoopRGE) Then 
betaTYw32 = (72*g1p2*TYw3adjYb3Yb3 - 120*g2p2*TYw3adjYb3Yb3 - 54*TrYb3adjYb3*TYw3adjYb3Yb3 -      & 
&  540*TrYuadjYu*TYw3adjYb3Yb3 - 270*TrYw3adjYw3*TYw3adjYb3Yb3 - 540*TrYx3adjYx3*TYw3adjYb3Yb3 -& 
&  72._dp*(TYw3adjYb3Yb3adjYb3Yb3) - 135._dp*(TYw3adjYb3Yb3adjYw3Yw3) + 240*g1p2*TYw3adjYeYe -& 
&  600*TrYdadjYd*TYw3adjYeYe - 200*TrYeadjYe*TYw3adjYeYe - 120._dp*(TYw3adjYeYeadjYb3Yb3) -& 
&  400._dp*(TYw3adjYeYeadjYeYe) - 300._dp*(TYw3adjYeYeadjYw3Yw3) + 180*g1p2*TYw3adjYw3Yw3 +& 
&  900*g2p2*TYw3adjYw3Yw3 - 180*TrYb3adjYb3*TYw3adjYw3Yw3 - 1800*TrYuadjYu*TYw3adjYw3Yw3 -& 
&  900*TrYw3adjYw3*TYw3adjYw3Yw3 - 1800*TrYx3adjYx3*TYw3adjYw3Yw3 - 225._dp*(TYw3adjYw3Yw3adjYb3Yb3) -& 
&  675._dp*(TYw3adjYw3Yw3adjYw3Yw3) + 36*g1p2*Yw3adjYb3TYb3 - 60*g2p2*Yw3adjYb3TYb3 -    & 
&  54*TrYb3adjYb3*Yw3adjYb3TYb3 - 540*TrYuadjYu*Yw3adjYb3TYb3 - 270*TrYw3adjYw3*Yw3adjYb3TYb3 -& 
&  540*TrYx3adjYx3*Yw3adjYb3TYb3 - 108._dp*(Yw3adjYb3TYb3adjYb3Yb3) - 180._dp*(Yw3adjYb3TYb3adjYw3Yw3) -& 
&  24*(3*g1p2*MassB - 5*g2p2*MassWB + 3._dp*(TradjYb3TYb3) + 30._dp*(TradjYuTYu) +       & 
&  15._dp*(TradjYw3TYw3) + 30._dp*(TradjYx3TYx3))*Yw3adjYb3Yb3 - 90._dp*(Yw3adjYb3Yb3adjYb3TYb3) -& 
&  135._dp*(Yw3adjYb3Yb3adjYw3TYw3) + 480*g1p2*Yw3adjYeTYe - 1200*TrYdadjYd*Yw3adjYeTYe -& 
&  400*TrYeadjYe*Yw3adjYeTYe - 120._dp*(Yw3adjYeTYeadjYb3Yb3) - 800._dp*(Yw3adjYeTYeadjYeYe) -& 
&  300._dp*(Yw3adjYeTYeadjYw3Yw3) - 480*g1p2*MassB*Yw3adjYeYe - 1200*TradjYdTYd*Yw3adjYeYe -& 
&  400*TradjYeTYe*Yw3adjYeYe - 60._dp*(Yw3adjYeYeadjYb3TYb3) - 800._dp*(Yw3adjYeYeadjYeTYe) -& 
&  150._dp*(Yw3adjYeYeadjYw3TYw3) + 90*g1p2*Yw3adjYw3TYw3 + 2250*g2p2*Yw3adjYw3TYw3 -    & 
&  225*TrYb3adjYb3*Yw3adjYw3TYw3 - 2250*TrYuadjYu*Yw3adjYw3TYw3 - 1125*TrYw3adjYw3*Yw3adjYw3TYw3 -& 
&  2250*TrYx3adjYx3*Yw3adjYw3TYw3 - 270._dp*(Yw3adjYw3TYw3adjYb3Yb3) - 900._dp*(Yw3adjYw3TYw3adjYw3Yw3) -& 
&  180*g1p2*MassB*Yw3adjYw3Yw3 - 2100*g2p2*MassWB*Yw3adjYw3Yw3 - 270*TradjYb3TYb3*Yw3adjYw3Yw3 -& 
&  2700*TradjYuTYu*Yw3adjYw3Yw3 - 1350*TradjYw3TYw3*Yw3adjYw3Yw3 - 2700*TradjYx3TYx3*Yw3adjYw3Yw3 -& 
&  180._dp*(Yw3adjYw3Yw3adjYb3TYb3) - 675._dp*(Yw3adjYw3Yw3adjYw3TYw3) - 4*Yw3*(54._dp*(TrYb3adjYb3TYb3adjYb3) +& 
&  30._dp*(TrYb3adjYeTYeadjYb3) + 135._dp*(TrYb3adjYw3TYw3adjYb3) + 300._dp*(TrYdadjYuTYuadjYd) +& 
&  30._dp*(TrYeadjYb3TYb3adjYe) + 150._dp*(TrYeadjYw3TYw3adjYe) + 300._dp*(TrYuadjYdTYdadjYu) +& 
&  1800._dp*(TrYuadjYuTYuadjYu) + 135._dp*(TrYw3adjYb3TYb3adjYw3) + 150._dp*(TrYw3adjYeTYeadjYw3) +& 
&  675._dp*(TrYw3adjYw3TYw3adjYw3) + 4*(207*g1p4*MassB + 1375*g2p4*MassWB +              & 
&  45*g1p2*g2p2*(MassB + MassWB) + 150._dp*(TradjYdTpYx3CYx3TYd) + 150._dp*(TradjYx3TYx3CYdTpYd) +& 
&  450._dp*(TrYx3adjYx3TYx3adjYx3)) + 300*g1p4*MassB*NGHx3 +              & 
&  20*(2*(-2*(g1p2 + 20._dp*(g3p2))*TradjYuTYu - 5*(g1p2 + 8._dp*(g3p2))*TradjYx3TYx3 +  & 
&  2*(g1p2*MassB + 20*g3p2*MassG)*TrYuadjYu + 5*(g1p2*MassB + 8*g3p2*MassG)*TrYx3adjYx3) +& 
&  15*g1p4*MassB*NGHxb3) + 100*g2p2*(-6._dp*(TradjYw3TYw3) +              & 
&  6*MassWB*TrYw3adjYw3 + 7*g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))) + (828._dp*(g1p4) + 360*g1p2*g2p2 + 5500._dp*(g2p4) -       & 
&  54._dp*(TrYb3adjYb3Yb3adjYb3) - 60._dp*(TrYb3adjYeYeadjYb3) - 270._dp*(TrYb3adjYw3Yw3adjYb3) -& 
&  1200._dp*(TrYdadjYdTpYx3CYx3) - 600._dp*(TrYdadjYuYuadjYd) - 300._dp*(TrYeadjYw3Yw3adjYe) -& 
&  1800._dp*(TrYuadjYuYuadjYu) + 1200*g2p2*TrYw3adjYw3 - 675._dp*(TrYw3adjYw3Yw3adjYw3) -& 
&  1800._dp*(TrYx3adjYx3Yx3adjYx3) + 20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu +            & 
&  20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 + 15*g1p4*NGHx3 + 15*g1p4*NGHxb3) +& 
&  700*g2p4*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))*TYw3)/200._dp

 
DTYw3 = oo16pi2*( betaTYw31 + oo16pi2 * betaTYw32 ) 

 
Else 
DTYw3 = oo16pi2* betaTYw31 
End If 
 
 
!-------------------- 
! TYx3 
!-------------------- 
 
betaTYx31  = 4._dp*(TYx3adjYx3Yx3) + 2._dp*(TYx3CYdTpYd) + ((38*g1p2*MassB)           & 
& /15._dp + (32*g3p2*MassG)/3._dp + 6*g2p2*MassWB + 3._dp*(TradjYb3TYb3)/5._dp +         & 
&  6._dp*(TradjYuTYu) + 3._dp*(TradjYw3TYw3) + 6._dp*(TradjYx3TYx3))*Yx3 +               & 
&  5._dp*(Yx3adjYx3TYx3) + 4._dp*(Yx3CYdTpTYd) + ((-38._dp*(g1p2) - 90._dp*(g2p2)        & 
&  - 160._dp*(g3p2) + 9._dp*(TrYb3adjYb3) + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3)    & 
&  + 90._dp*(TrYx3adjYx3))*TYx3)/30._dp

 
 
If (TwoLoopRGE) Then 
betaTYx32 = (180*(12*(g1p2 + 5._dp*(g2p2) - TrYb3adjYb3 - 10._dp*(TrYuadjYu) - 5._dp*(TrYw3adjYw3) -& 
&  10._dp*(TrYx3adjYx3))*TYx3adjYx3Yx3 + 4*(g1p2 + 15._dp*(g2p2) - 15._dp*(TrYdadjYd) -  & 
&  5._dp*(TrYeadjYe))*TYx3CYdTpYd + 3*(12._dp*(g1p2) + 40._dp*(g2p2) - 5._dp*(TrYb3adjYb3) -& 
&  50._dp*(TrYuadjYu) - 25._dp*(TrYw3adjYw3) - 50._dp*(TrYx3adjYx3))*Yx3adjYx3TYx3 -     & 
&  2*(16*g1p2*MassB + 60*g2p2*MassWB + 9._dp*(TradjYb3TYb3) + 90._dp*(TradjYuTYu) +      & 
&  45._dp*(TradjYw3TYw3) + 90._dp*(TradjYx3TYx3))*Yx3adjYx3Yx3 + 8*(g1p2 +               & 
&  15._dp*(g2p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*Yx3CYdTpTYd - 8*(g1p2*MassB + & 
&  15*g2p2*MassWB + 15._dp*(TradjYdTYd) + 5._dp*(TradjYeTYe))*Yx3CYdTpYd) -              & 
&  3600*(3._dp*(TYx3adjYx3Yx3adjYx3Yx3) + 2._dp*(TYx3CYdTpYdadjYx3Yx3) + TYx3CYdTpYdCYdTpYd +& 
&  TYx3CYdTpYuCYuTpYd + 4._dp*(Yx3adjYx3TYx3adjYx3Yx3) + 3._dp*(Yx3adjYx3Yx3adjYx3TYx3) +& 
&  2._dp*(Yx3CYdTpTYdCYdTpYd) + 2._dp*(Yx3CYdTpTYuCYuTpYd) + Yx3CYdTpYdadjYx3TYx3 +      & 
&  2._dp*(Yx3CYdTpYdCYdTpTYd) + 2._dp*(Yx3CYdTpYuCYuTpTYd) + 2._dp*(Yx3TYdadjYdadjYx3Yx3)) -& 
&  4*Yx3*(486._dp*(TrYb3adjYb3TYb3adjYb3) + 270._dp*(TrYb3adjYeTYeadjYb3) +              & 
&  1215._dp*(TrYb3adjYw3TYw3adjYb3) + 2700._dp*(TrYdadjYuTYuadjYd) + 270._dp*(TrYeadjYb3TYb3adjYe) +& 
&  1350._dp*(TrYeadjYw3TYw3adjYe) + 2700._dp*(TrYuadjYdTYdadjYu) + 16200._dp*(TrYuadjYuTYuadjYu) +& 
&  1215._dp*(TrYw3adjYb3TYb3adjYw3) + 1350._dp*(TrYw3adjYeTYeadjYw3) + 6075._dp*(TrYw3adjYw3TYw3adjYw3) +& 
&  4*(4123*g1p4*MassB + 5*g1p2*(232*g3p2*(MassB + MassG) + 153*g2p2*(MassB +             & 
&  MassWB)) + 25*(-32*g3p4*MassG + 135*g2p4*MassWB + 72*g2p2*g3p2*(MassG +               & 
&  MassWB)) + 1350._dp*(TradjYdTpYx3CYx3TYd) + 1350._dp*(TradjYx3TYx3CYdTpYd) +          & 
&  4050._dp*(TrYx3adjYx3TYx3adjYx3)) + 28800*g3p4*MassG*NGHg3 +           & 
&  60*(6*(-2*(g1p2 + 20._dp*(g3p2))*TradjYuTYu - 5*(g1p2 + 8._dp*(g3p2))*TradjYx3TYx3 +  & 
&  2*(g1p2*MassB + 20*g3p2*MassG)*TrYuadjYu + 5*(g1p2*MassB + 8*g3p2*MassG)*TrYx3adjYx3) +& 
&  5*(19*g1p4*MassB + 32*g3p4*MassG)*NGHx3 + 5*(19*g1p4*MassB +           & 
&  32*g3p4*MassG)*NGHxb3) + 2700*g2p2*(-2._dp*(TradjYw3TYw3) +            & 
&  2*MassWB*TrYw3adjYw3 + g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))) + (-486._dp*(TrYb3adjYb3Yb3adjYb3) - 540._dp*(TrYb3adjYeYeadjYb3) -& 
&  2430._dp*(TrYb3adjYw3Yw3adjYb3) - 10800._dp*(TrYdadjYdTpYx3CYx3) - 5400._dp*(TrYdadjYuYuadjYd) -& 
&  2700._dp*(TrYeadjYw3Yw3adjYe) - 16200._dp*(TrYuadjYuYuadjYu) - 6075._dp*(TrYw3adjYw3Yw3adjYw3) +& 
&  4*(4123._dp*(g1p4) + 1530*g1p2*g2p2 + 3375._dp*(g2p4) + 80*(29._dp*(g1p2) +           & 
&  45._dp*(g2p2))*g3p2 - 800._dp*(g3p4) - 4050._dp*(TrYx3adjYx3Yx3adjYx3)) +             & 
&  60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) + 95*g1p2*(NGHx3 +  & 
&  NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu + TrYx3adjYx3) + g3p2*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)) + 2700*g2p2*(4._dp*(TrYw3adjYw3) +     & 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3))))*TYx3)/1800._dp

 
DTYx3 = oo16pi2*( betaTYx31 + oo16pi2 * betaTYx32 ) 

 
Else 
DTYx3 = oo16pi2* betaTYx31 
End If 
 
 
!-------------------- 
! Bmue 
!-------------------- 
 
betaBmue1  = ((6*g1p2*MassB + 30*g2p2*MassWB + 3._dp*(TradjYb3TYb3) + 30._dp*(TradjYdTYd)& 
&  + 10._dp*(TradjYeTYe) + 30._dp*(TradjYuTYu) + 15._dp*(TradjYw3TYw3) + 30._dp*(TradjYx3TYx3))& 
& *mue)/5._dp + (-3._dp*(g1p2)/5._dp - 3._dp*(g2p2) + 3._dp*(TrYb3adjYb3)/10._dp +       & 
&  3._dp*(TrYdadjYd) + TrYeadjYe + 3._dp*(TrYuadjYu) + 3._dp*(TrYw3adjYw3)               & 
& /2._dp + 3._dp*(TrYx3adjYx3))*Bmue

 
 
If (TwoLoopRGE) Then 
betaBmue2 = (-48*(69*g1p4*MassB + 125*g2p4*MassWB + 15*g1p2*g2p2*(MassB + MassWB))*mue +          & 
&  Bmue*(80*(-((g1p2 - 40._dp*(g3p2))*TrYdadjYd) + 3*g1p2*TrYeadjYe + 2*(g1p2 +          & 
&  20._dp*(g3p2))*TrYuadjYu + 5*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3) + 3*(276._dp*(g1p4) + & 
&  120*g1p2*g2p2 + 500._dp*(g2p4) - 18._dp*(TrYb3adjYb3Yb3adjYb3) - 40._dp*(TrYb3adjYeYeadjYb3) -& 
&  90._dp*(TrYb3adjYw3Yw3adjYb3) - 800._dp*(TrYdadjYdTpYx3CYx3) - 600._dp*(TrYdadjYdYdadjYd) -& 
&  400._dp*(TrYdadjYuYuadjYd) - 200._dp*(TrYeadjYeYeadjYe) - 200._dp*(TrYeadjYw3Yw3adjYe) -& 
&  600._dp*(TrYuadjYuYuadjYu) - 225._dp*(TrYw3adjYw3Yw3adjYw3) - 600._dp*(TrYx3adjYx3Yx3adjYx3)) +& 
&  300*g1p4*NGHx3 + 300*g1p4*NGHxb3 + 300*g2p2*(4._dp*(TrYw3adjYw3) +& 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))) -  & 
&  4*mue*(40*(g1p2*TradjYdTYd - 40*g3p2*TradjYdTYd - 3*g1p2*TradjYeTYe - 2*g1p2*TradjYuTYu -& 
&  40*g3p2*TradjYuTYu - 5*(g1p2 + 8._dp*(g3p2))*TradjYx3TYx3 + (-(g1p2*MassB) +          & 
&  40*g3p2*MassG)*TrYdadjYd + 3*g1p2*MassB*TrYeadjYe + 2*g1p2*MassB*TrYuadjYu +          & 
&  40*g3p2*MassG*TrYuadjYu + 5*g1p2*MassB*TrYx3adjYx3 + 40*g3p2*MassG*TrYx3adjYx3) +     & 
&  3*(18._dp*(TrYb3adjYb3TYb3adjYb3) + 5*(4._dp*(TrYb3adjYeTYeadjYb3) + 9._dp*(TrYb3adjYw3TYw3adjYb3) +& 
&  120._dp*(TrYdadjYdTYdadjYd) + 40._dp*(TrYdadjYuTYuadjYd) + 4._dp*(TrYeadjYb3TYb3adjYe) +& 
&  40._dp*(TrYeadjYeTYeadjYe) + 20._dp*(TrYeadjYw3TYw3adjYe) + 40._dp*(TrYuadjYdTYdadjYu) +& 
&  120._dp*(TrYuadjYuTYuadjYu) + 9._dp*(TrYw3adjYb3TYb3adjYw3) + 5*(4._dp*(TrYw3adjYeTYeadjYw3) +& 
&  9._dp*(TrYw3adjYw3TYw3adjYw3) + 8*(2*(TradjYdTpYx3CYx3TYd + TradjYx3TYx3CYdTpYd) +    & 
&  3._dp*(TrYx3adjYx3TYx3adjYx3))))) + 300*g1p4*MassB*NGHx3 +             & 
&  300*g1p4*MassB*NGHxb3 + 300*g2p2*(-2._dp*(TradjYw3TYw3) +              & 
&  2*MassWB*TrYw3adjYw3 + g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))))/200._dp

 
DBmue = oo16pi2*( betaBmue1 + oo16pi2 * betaBmue2 ) 

 
Else 
DBmue = oo16pi2* betaBmue1 
End If 
 
 
!-------------------- 
! BMXM3 
!-------------------- 
 
!betaBMXM31  = (2*(5*g1p2*MassB + 16*g3p2*MassG + 9*g2p2*MassWB)*MXM3 + 3*(BMXM3CYx3TpYx3 +& 
!&  2._dp*(MXM3CYx3TpTYx3)) - (5._dp*(g1p2) + 9._dp*(g2p2) + 16._dp*(g3p2))               & 
!& *BMXM3)/3._dp

 
 
! If (TwoLoopRGE) Then 
! betaBMXM32 = -2*(BMXM3CYx3TpYx3CYx3TpYx3 + BMXM3CYx3YdadjYdTpYx3 + 2._dp*(MXM3CYx3TpTYx3CYx3TpYx3) +& 
! &  2._dp*(MXM3CYx3TpYx3CYx3TpTYx3) + 2._dp*(MXM3CYx3TYdadjYdTpYx3) + 2._dp*(MXM3CYx3YdadjYdTpTYx3)) -& 
! &  (2*MXM3*(223*g1p4*MassB - 32*g3p4*MassG + 135*g2p4*MassWB + 144*g2p2*g3p2*(MassG +    & 
! &  MassWB) + 5*g1p2*(16*g3p2*(MassB + MassG) + 9*g2p2*(MassB + MassWB)) + 288*g3p4*MassG*NGHg3 +& 
! &  108*g2p4*MassWB*NGHw3 + 75*g1p4*MassB*(NGHx3 +          & 
! &  NGHxb3) + 3*(32*g3p4*MassG + 27*g2p4*MassWB)*(NGHx3 +   & 
! &  NGHxb3)))/9._dp + (18*MXM3CYx3TpYx3*(4*g1p2*MassB - 3._dp*(TradjYb3TYb3) -& 
! &  30._dp*(TradjYuTYu) - 15._dp*(TradjYw3TYw3) - 30._dp*(TradjYx3TYx3)) - 9*(BMXM3CYx3TpYx3 +& 
! &  2._dp*(MXM3CYx3TpTYx3))*(4._dp*(g1p2) + 3._dp*(TrYb3adjYb3) + 30._dp*(TrYuadjYu) +    & 
! &  15._dp*(TrYw3adjYw3) + 30._dp*(TrYx3adjYx3)) + 5*BMXM3*(223._dp*(g1p4) +              & 
! &  90*g1p2*g2p2 + 135._dp*(g2p4) + 32*(5._dp*(g1p2) + 9._dp*(g2p2))*g3p2 -               & 
! &  32._dp*(g3p4) + 108*g2p4*NGHw3 + 3*(25._dp*(g1p4) + 27._dp*(g2p4))*(NGHx3 +& 
! &  NGHxb3) + 96*g3p4*(3*NGHg3 + NGHx3 +     & 
! &  NGHxb3)))/90._dp
! 
!  
! DBMXM3 = oo16pi2*( betaBMXM31 + oo16pi2 * betaBMXM32 ) 
! 
!  
! Else 
! DBMXM3 = oo16pi2* betaBMXM31 
! End If 
  DBMXM3 = 0._dp
 
!-------------------- 
! BMWM3 
!-------------------- 
 
! betaBMWM31  = BMWM3CYw3TpYw3 + 2._dp*(MWM3CYw3TpTYw3) + 2._dp*(TYw3adjYw3MWM3)        & 
! &  + Yw3adjYw3BMWM3 + 8*g2p2*(2*MassWB*MWM3 - BMWM3)
! 
!  
!  
! If (TwoLoopRGE) Then 
! betaBMWM32 = (-3._dp*(BMWM3CYw3TpYb3CYb3TpYw3) - 10._dp*(BMWM3CYw3TpYeCYeTpYw3) - 15._dp*(BMWM3CYw3TpYw3CYw3TpYw3) +& 
! &  6*BMWM3CYw3TpYw3*g1p2 - 10*BMWM3CYw3TpYw3*g2p2 - 6._dp*(MWM3CYw3TpTYb3CYb3TpYw3) -    & 
! &  20._dp*(MWM3CYw3TpTYeCYeTpYw3) - 30._dp*(MWM3CYw3TpTYw3CYw3TpYw3) - 6._dp*(MWM3CYw3TpYb3CYb3TpTYw3) -& 
! &  20._dp*(MWM3CYw3TpYeCYeTpTYw3) - 30._dp*(MWM3CYw3TpYw3CYw3TpTYw3) - 2*MWM3CYw3TpYw3*(6*g1p2*MassB -& 
! &  10*g2p2*MassWB + 3._dp*(TradjYb3TYb3) + 30._dp*(TradjYuTYu) + 15._dp*(TradjYw3TYw3) + & 
! &  30._dp*(TradjYx3TYx3)) - 3*BMWM3CYw3TpYw3*TrYb3adjYb3 - 30*BMWM3CYw3TpYw3*TrYuadjYu - & 
! &  15*BMWM3CYw3TpYw3*TrYw3adjYw3 - 30*BMWM3CYw3TpYw3*TrYx3adjYx3 - MWM3CYw3TpTYw3*(-     & 
! & 12._dp*(g1p2) + 20._dp*(g2p2) + 6._dp*(TrYb3adjYb3) + 60._dp*(TrYuadjYu) +             & 
! &  30._dp*(TrYw3adjYw3) + 60._dp*(TrYx3adjYx3)) - 6._dp*(TYw3adjYb3Yb3adjYw3MWM3) -      & 
! &  20._dp*(TYw3adjYeYeadjYw3MWM3) + 2*(6._dp*(g1p2) - 10._dp*(g2p2) - 3._dp*(TrYb3adjYb3) -& 
! &  30._dp*(TrYuadjYu) - 15._dp*(TrYw3adjYw3) - 30._dp*(TrYx3adjYx3))*TYw3adjYw3MWM3 -    & 
! &  30._dp*(TYw3adjYw3Yw3adjYw3MWM3) - 6._dp*(Yw3adjYb3TYb3adjYw3MWM3) - 3._dp*(Yw3adjYb3Yb3adjYw3BMWM3) -& 
! &  20._dp*(Yw3adjYeTYeadjYw3MWM3) - 10._dp*(Yw3adjYeYeadjYw3BMWM3) + 6*g1p2*Yw3adjYw3BMWM3 -& 
! &  10*g2p2*Yw3adjYw3BMWM3 - 3*TrYb3adjYb3*Yw3adjYw3BMWM3 - 30*TrYuadjYu*Yw3adjYw3BMWM3 - & 
! &  15*TrYw3adjYw3*Yw3adjYw3BMWM3 - 30*TrYx3adjYx3*Yw3adjYw3BMWM3 - 12*g1p2*MassB*Yw3adjYw3MWM3 +& 
! &  20*g2p2*MassWB*Yw3adjYw3MWM3 - 6*TradjYb3TYb3*Yw3adjYw3MWM3 - 60*TradjYuTYu*Yw3adjYw3MWM3 -& 
! &  30*TradjYw3TYw3*Yw3adjYw3MWM3 - 60*TradjYx3TYx3*Yw3adjYw3MWM3 - 30._dp*(Yw3adjYw3TYw3adjYw3MWM3) -& 
! &  15._dp*(Yw3adjYw3Yw3adjYw3BMWM3))/10._dp - 16*g2p4*MassWB*MWM3*(10 + 4*NGHw3 +& 
! &  3*NGHx3 + 3*NGHxb3) + 4*g2p4*BMWM3*(10 + 4*NGHw3 +& 
! &  3*NGHx3 + 3*NGHxb3)
! 
!  
! DBMWM3 = oo16pi2*( betaBMWM31 + oo16pi2 * betaBMWM32 ) 
! 
!  
! Else 
! DBMWM3 = oo16pi2* betaBMWM31 
! End If 
 DBMWM3 = 0._dp
 
!-------------------- 
! BMGM3 
!-------------------- 
 
! betaBMGM31  = 12*g3p2*(2*MassG*MGM3 - BMGM3)
! 
!  
!  
! If (TwoLoopRGE) Then 
! betaBMGM32 = -12*g3p4*(4*MassG*MGM3 - BMGM3)*(3 + 3*NGHg3 + NGHx3 +  & 
! &  NGHxb3)
! 
!  
! DBMGM3 = oo16pi2*( betaBMGM31 + oo16pi2 * betaBMGM32 ) 
! 
!  
! Else 
! DBMGM3 = oo16pi2* betaBMGM31 
! End If 
 DBMGM3 = 0._dp
 
!-------------------- 
! BMBM3 
!-------------------- 
 
! betaBMBM31  = (3*(BMBM3CYb3TpYb3 + 2._dp*(MBM3CYb3TpTYb3) + 2._dp*(TYb3adjYb3MBM3)    & 
! &  + Yb3adjYb3BMBM3))/5._dp
! 
!  
!  
! If (TwoLoopRGE) Then 
! betaBMBM32 = (-3*(3._dp*(BMBM3CYb3TpYb3CYb3TpYb3) + 10._dp*(BMBM3CYb3TpYeCYeTpYb3) +               & 
! &  15._dp*(BMBM3CYb3TpYw3CYw3TpYb3) + 6._dp*(MBM3CYb3TpTYb3CYb3TpYb3) + 20._dp*(MBM3CYb3TpTYeCYeTpYb3) +& 
! &  30._dp*(MBM3CYb3TpTYw3CYw3TpYb3) + 6._dp*(MBM3CYb3TpYb3CYb3TpTYb3) + 20._dp*(MBM3CYb3TpYeCYeTpTYb3) +& 
! &  30._dp*(MBM3CYb3TpYw3CYw3TpTYb3) + 6*MBM3CYb3TpYb3*(2*g1p2*MassB + 10*g2p2*MassWB +   & 
! &  TradjYb3TYb3 + 10._dp*(TradjYuTYu) + 5._dp*(TradjYw3TYw3) + 10._dp*(TradjYx3TYx3)) +  & 
! &  6*MBM3CYb3TpTYb3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) -            & 
! &  2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3))) + 6._dp*(TYb3adjYb3Yb3adjYb3MBM3) +    & 
! &  20._dp*(TYb3adjYeYeadjYb3MBM3) + 30._dp*(TYb3adjYw3Yw3adjYb3MBM3) + 3*((TrYb3adjYb3 + & 
! &  10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))*(BMBM3CYb3TpYb3 +& 
! &  2._dp*(TYb3adjYb3MBM3) + Yb3adjYb3BMBM3) + 2*(2*g1p2*MassB + 10*g2p2*MassWB +         & 
! &  TradjYb3TYb3 + 10._dp*(TradjYuTYu) + 5._dp*(TradjYw3TYw3) + 10._dp*(TradjYx3TYx3))*Yb3adjYb3MBM3) +& 
! &  6._dp*(Yb3adjYb3TYb3adjYb3MBM3) + 3._dp*(Yb3adjYb3Yb3adjYb3BMBM3) + 20._dp*(Yb3adjYeTYeadjYb3MBM3) +& 
! &  10._dp*(Yb3adjYeYeadjYb3BMBM3) + 30._dp*(Yb3adjYw3TYw3adjYb3MBM3) + 15._dp*(Yb3adjYw3Yw3adjYb3BMBM3)))/50._dp
! 
!  
! DBMBM3 = oo16pi2*( betaBMBM31 + oo16pi2 * betaBMBM32 ) 
! 
!  
! Else 
! DBMBM3 = oo16pi2* betaBMBM31 
! End If 
 DBMBM3 = 0._dp
 
!-------------------- 
! mq2 
!-------------------- 
 
betamq21  = mq2TpYdCYd + mq2TpYuCYu + 2._dp*(TpTYdCTYd) + 2._dp*(TpTYuCTYu)           & 
&  + 2*mHd2*TpYdCYd + TpYdCYdmq2 + 2._dp*(TpYdmd2CYd) + 2*mHu2*TpYuCYu + TpYuCYumq2 +    & 
&  2._dp*(TpYumu2CYu) - (id3R*(90*AbsMassWB*g2p2 + 160*AbsMassG*g3p2 + g1p2*(2._dp*(AbsMassB)& 
&  - 5*Tr1(1))))/15._dp

 
 
If (TwoLoopRGE) Then 
betamq22 = (2*AbsMassWB*g1p2*g2p2*id3R)/5._dp + 33*AbsMassWB*g2p4*id3R + 32*AbsMassWB*g2p2*g3p2*id3R +& 
&  (2*g1p2*mq2TpYdCYd)/5._dp + (4*g1p2*mq2TpYuCYu)/5._dp + (4*g1p2*TpTYdCTYd)/5._dp +    & 
&  (8*g1p2*TpTYuCTYu)/5._dp - 4*mHd2*TpYdadjYx3Yx3CYd - 4*mHu2*TpYdadjYx3Yx3CYd -        & 
&  (4*g1p2*MassB*TpYdCTYd)/5._dp + (4*g1p2*mHd2*TpYdCYd)/5._dp + (2*g1p2*TpYdCYdmq2)/5._dp -& 
&  8*mHd2*TpYdCYdTpYdCYd + (4*g1p2*TpYdmd2CYd)/5._dp - (8*g1p2*MassB*TpYuCTYu)/5._dp +   & 
&  (8*g1p2*mHu2*TpYuCYu)/5._dp + (4*g1p2*TpYuCYumq2)/5._dp - 8*mHu2*TpYuCYuTpYuCYu +     & 
&  (8*g1p2*TpYumu2CYu)/5._dp - (3*TpYuCTYu*TradjYb3TYb3)/5._dp - 6*TpYdCTYd*TradjYdTYd - & 
&  2*TpYdCTYd*TradjYeTYe - 6*TpYuCTYu*TradjYuTYu - 3*TpYuCTYu*TradjYw3TYw3 -             & 
&  6*TpYuCTYu*TradjYx3TYx3 - (3*TpYuCYu*TrCTYb3TpTYb3)/5._dp - 6*TpYdCYd*TrCTYdTpTYd -   & 
&  2*TpYdCYd*TrCTYeTpTYe - 2*TpTYdCYd*(3._dp*(TrCTYdTpYd) + TrCTYeTpYe) - 6*TpYuCYu*TrCTYuTpTYu -& 
&  3*TpYuCYu*TrCTYw3TpTYw3 - 6*TpYuCYu*TrCTYx3TpTYx3 - (3*TpTYuCYu*(TrCTYb3TpYb3 +       & 
&  5*(2._dp*(TrCTYuTpYu) + TrCTYw3TpYw3 + 2._dp*(TrCTYx3TpYx3))))/5._dp - 6*TpYuCYu*Trmd2adjYx3Yx3 -& 
&  6*TpYdCYd*Trmd2YdadjYd - 2*TpYdCYd*Trme2YeadjYe - (3*TpYuCYu*TrmHb32Yb3adjYb3)/5._dp -& 
&  3*TpYuCYu*TrmHw32Yw3adjYw3 - 6*TpYuCYu*TrmHxb32Yx3adjYx3 - (3*TpYuCYu*Trml2adjYb3Yb3)/5._dp -& 
&  2*TpYdCYd*Trml2adjYeYe - 3*TpYuCYu*Trml2adjYw3Yw3 - 6*TpYdCYd*Trmq2adjYdYd -          & 
&  6*TpYuCYu*Trmq2adjYuYu - 6*TpYuCYu*Trmu2YuadjYu - (3*mq2TpYuCYu*TrYb3adjYb3)/10._dp - & 
&  (3*TpTYuCTYu*TrYb3adjYb3)/5._dp - (6*mHu2*TpYuCYu*TrYb3adjYb3)/5._dp - (3*TpYuCYumq2*TrYb3adjYb3)/10._dp -& 
&  (3*TpYumu2CYu*TrYb3adjYb3)/5._dp - 3*mq2TpYdCYd*TrYdadjYd - 6*TpTYdCTYd*TrYdadjYd -   & 
&  12*mHd2*TpYdCYd*TrYdadjYd - 3*TpYdCYdmq2*TrYdadjYd - 6*TpYdmd2CYd*TrYdadjYd -         & 
&  mq2TpYdCYd*TrYeadjYe - 2*TpTYdCTYd*TrYeadjYe - 4*mHd2*TpYdCYd*TrYeadjYe -             & 
&  TpYdCYdmq2*TrYeadjYe - 2*TpYdmd2CYd*TrYeadjYe - 3*mq2TpYuCYu*TrYuadjYu  
betamq22 =  betamq22- 6*TpTYuCTYu*TrYuadjYu - 12*mHu2*TpYuCYu*TrYuadjYu - 3*TpYuCYumq2*TrYuadjYu -          & 
&  6*TpYumu2CYu*TrYuadjYu - (3*mq2TpYuCYu*TrYw3adjYw3)/2._dp - 3*TpTYuCTYu*TrYw3adjYw3 - & 
&  6*mHu2*TpYuCYu*TrYw3adjYw3 - (3*TpYuCYumq2*TrYw3adjYw3)/2._dp - 3*TpYumu2CYu*TrYw3adjYw3 -& 
&  3*mq2TpYuCYu*TrYx3adjYx3 - 6*TpTYuCTYu*TrYx3adjYx3 - 12*mHu2*TpYuCYu*TrYx3adjYx3 -    & 
&  3*TpYuCYumq2*TrYx3adjYx3 - 6*TpYumu2CYu*TrYx3adjYx3 - 2*(2._dp*(adjTYdTYdTpYdCYd) +   & 
&  2._dp*(adjTYuTYuTpYuCYu) + mq2TpYdCYdTpYdCYd + mq2TpYuCYuTpYuCYu + 2._dp*(TpTYdadjYdTpYx3CTYx3) +& 
&  2._dp*(TpTYdadjYx3Yx3CTYd) + 2._dp*(TpTYdCYdTpYdCTYd) + 2._dp*(TpTYuCYuTpYuCTYu) +    & 
&  2._dp*(TpYdadjYdTpTYx3CTYx3) + 2._dp*(TpYdadjYx3mHxb32Yx3CYd) + 2._dp*(TpYdadjYx3TYx3CTYd) +& 
&  TpYdadjYx3Yx3CYdmq2 + 2._dp*(TpYdadjYx3Yx3md2CYd) + 2._dp*(TpYdCTYdTpTYdCYd) +        & 
&  2._dp*(TpYdCYdmq2TpYdCYd) + 2._dp*(TpYdCYdTpTYdCTYd) + TpYdCYdTpYdCYdmq2 +            & 
&  2._dp*(TpYdCYdTpYdmd2CYd) + 2._dp*(TpYdmd2adjYx3Yx3CYd) + 2._dp*(TpYdmd2CYdTpYdCYd) + & 
&  2._dp*(TpYuCTYuTpTYuCYu) + 2._dp*(TpYuCYuTpTYuCTYu) + TpYuCYuTpYuCYumq2 +             & 
&  2*(TpYuCYumq2TpYuCYu + TpYuCYuTpYumu2CYu + TpYumu2CYuTpYuCYu) + Ydmq2adjYx3Yx3CYd) +  & 
&  (g1p2*g2p2*id3R*MassB*Conjg(MassWB))/5._dp + 16*g2p2*g3p2*id3R*MassG*Conjg(MassWB) +  & 
&  36*AbsMassWB*g2p4*id3R*NGHw3 + 27*AbsMassWB*g2p4*id3R*NGHx3 +& 
&  27*AbsMassWB*g2p4*id3R*NGHxb3 + (16*g3p2*id3R*Conjg(MassG)*(g1p2*(MassB +& 
&  2._dp*(MassG)) + 15*(-8*g3p2*MassG + 3*g2p2*(2._dp*(MassG) + MassWB)) +               & 
&  90*g3p2*MassG*(3*NGHg3 + NGHx3 + NGHxb3)))/45._dp +& 
&  (g1p2*Conjg(MassB)*(-180*(TpTYdCYd + 2._dp*(TpTYuCYu)) + 360*MassB*(TpYdCYd +         & 
&  2._dp*(TpYuCYu)) + id3R*(597*g1p2*MassB + 5*(16*g3p2*(2._dp*(MassB) + MassG) +        & 
&  9*g2p2*(2._dp*(MassB) + MassWB)) + 225*g1p2*MassB*(NGHx3 +             & 
&  NGHxb3))))/225._dp + (2*g1p4*id3R*Tr2(1))/15._dp + 6*g2p4*id3R*Tr2(2)  
betamq22 =  betamq22+ (32*g3p4*id3R*Tr2(3))/3._dp + (4*g1p2*id3R*Tr3(1))/3._dp

 
Dmq2 = oo16pi2*( betamq21 + oo16pi2 * betamq22 ) 

 
Else 
Dmq2 = oo16pi2* betamq21 
End If 
 
 
Forall(i1=1:3) Dmq2(i1,i1) =  Real(Dmq2(i1,i1),dp) 
!-------------------- 
! ml2 
!-------------------- 
 
betaml21  = 3._dp*(ml2TpYb3CYb3)/10._dp + ml2TpYeCYe + 3._dp*(ml2TpYw3CYw3)           & 
& /2._dp + 3._dp*(TpTYb3CTYb3)/5._dp + 2._dp*(TpTYeCTYe) + 3._dp*(TpTYw3CTYw3)           & 
&  + (3*mHu2*TpYb3CYb3)/5._dp + 3._dp*(TpYb3CYb3ml2)/10._dp + 3._dp*(TpYb3mHb32CYb3)     & 
& /5._dp + 2*mHd2*TpYeCYe + TpYeCYeml2 + 2._dp*(TpYeme2CYe) + 3*mHu2*TpYw3CYw3 +         & 
&  3._dp*(TpYw3CYw3ml2)/2._dp + 3._dp*(TpYw3mHw32CYw3) - (id3R*(30*AbsMassWB*g2p2 +      & 
&  g1p2*(6._dp*(AbsMassB) + 5*Tr1(1))))/5._dp

 
 
If (TwoLoopRGE) Then 
betaml22 = (6*g1p2*ml2TpYeCYe)/5._dp + 6*g2p2*ml2TpYw3CYw3 + (12*g1p2*TpTYeCTYe)/5._dp +         & 
&  12*g2p2*TpTYw3CTYw3 - (18*mHu2*TpYb3CYb3TpYb3CYb3)/25._dp - (9*mHu2*TpYb3CYb3TpYw3CYw3)/10._dp -& 
&  (12*g1p2*MassB*TpYeCTYe)/5._dp + (12*g1p2*mHd2*TpYeCYe)/5._dp + (6*g1p2*TpYeCYeml2)/5._dp -& 
&  8*mHd2*TpYeCYeTpYeCYe + (12*g1p2*TpYeme2CYe)/5._dp - 12*g2p2*MassWB*TpYw3CTYw3 +      & 
&  12*g2p2*mHu2*TpYw3CYw3 + 6*g2p2*TpYw3CYw3ml2 - (9*mHu2*TpYw3CYw3TpYb3CYb3)/10._dp -   & 
&  (9*mHu2*TpYw3CYw3TpYw3CYw3)/2._dp + 12*g2p2*TpYw3mHw32CYw3 - (9*TpYw3CTYw3*TradjYb3TYb3)/10._dp -& 
&  6*TpYeCTYe*TradjYdTYd - 2*TpYeCTYe*TradjYeTYe - 9*TpYw3CTYw3*TradjYuTYu -             & 
&  (9*TpYw3CTYw3*TradjYw3TYw3)/2._dp - 9*TpYw3CTYw3*TradjYx3TYx3 - (9*TpYw3CYw3*TrCTYb3TpTYb3)/10._dp -& 
&  6*TpYeCYe*TrCTYdTpTYd - 2*TpYeCYe*TrCTYeTpTYe - 9*TpYw3CYw3*TrCTYuTpTYu -             & 
&  (9*TpYw3CYw3*TrCTYw3TpTYw3)/2._dp - 9*TpYw3CYw3*TrCTYx3TpTYx3 - 9*TpYw3CYw3*Trmd2adjYx3Yx3 -& 
&  6*TpYeCYe*Trmd2YdadjYd - 2*TpYeCYe*Trme2YeadjYe - (9*TpYw3CYw3*TrmHb32Yb3adjYb3)/10._dp -& 
&  (9*TpYw3CYw3*TrmHw32Yw3adjYw3)/2._dp - 9*TpYw3CYw3*TrmHxb32Yx3adjYx3 - (9*TpYw3CYw3*Trml2adjYb3Yb3)/10._dp -& 
&  2*TpYeCYe*Trml2adjYeYe - (9*TpYw3CYw3*Trml2adjYw3Yw3)/2._dp - 6*TpYeCYe*Trmq2adjYdYd -& 
&  9*TpYw3CYw3*Trmq2adjYuYu - 9*TpYw3CYw3*Trmu2YuadjYu - (9*ml2TpYb3CYb3*TrYb3adjYb3)/100._dp -& 
&  (9*ml2TpYw3CYw3*TrYb3adjYb3)/20._dp - (9*TpTYb3CTYb3*TrYb3adjYb3)/50._dp -            & 
&  (9*TpTYw3CTYw3*TrYb3adjYb3)/10._dp - (9*TpYb3CYb3ml2*TrYb3adjYb3)/100._dp -           & 
&  (9*TpYb3mHb32CYb3*TrYb3adjYb3)/50._dp - (9*mHu2*TpYw3CYw3*TrYb3adjYb3)/5._dp -        & 
&  (9*TpYw3CYw3ml2*TrYb3adjYb3)/20._dp - (9*TpYw3mHw32CYw3*TrYb3adjYb3)/10._dp -         & 
&  3*ml2TpYeCYe*TrYdadjYd - 6*TpTYeCTYe*TrYdadjYd - 12*mHd2*TpYeCYe*TrYdadjYd -          & 
&  3*TpYeCYeml2*TrYdadjYd - 6*TpYeme2CYe*TrYdadjYd - ml2TpYeCYe*TrYeadjYe -              & 
&  2*TpTYeCTYe*TrYeadjYe - 4*mHd2*TpYeCYe*TrYeadjYe - TpYeCYeml2*TrYeadjYe  
betaml22 =  betaml22- 2*TpYeme2CYe*TrYeadjYe - (9*ml2TpYb3CYb3*TrYuadjYu)/10._dp - (9*ml2TpYw3CYw3*TrYuadjYu)/2._dp -& 
&  (9*TpTYb3CTYb3*TrYuadjYu)/5._dp - 9*TpTYw3CTYw3*TrYuadjYu - (9*TpYb3CYb3ml2*TrYuadjYu)/10._dp -& 
&  (9*TpYb3mHb32CYb3*TrYuadjYu)/5._dp - 18*mHu2*TpYw3CYw3*TrYuadjYu - (9*TpYw3CYw3ml2*TrYuadjYu)/2._dp -& 
&  9*TpYw3mHw32CYw3*TrYuadjYu - (9*ml2TpYb3CYb3*TrYw3adjYw3)/20._dp - (9*ml2TpYw3CYw3*TrYw3adjYw3)/4._dp -& 
&  (9*TpTYb3CTYb3*TrYw3adjYw3)/10._dp - (9*TpTYw3CTYw3*TrYw3adjYw3)/2._dp -              & 
&  (9*TpYb3CYb3ml2*TrYw3adjYw3)/20._dp - (9*TpYb3mHb32CYb3*TrYw3adjYw3)/10._dp -         & 
&  9*mHu2*TpYw3CYw3*TrYw3adjYw3 - (9*TpYw3CYw3ml2*TrYw3adjYw3)/4._dp - (9*TpYw3mHw32CYw3*TrYw3adjYw3)/2._dp -& 
&  (9*ml2TpYb3CYb3*TrYx3adjYx3)/10._dp - (9*ml2TpYw3CYw3*TrYx3adjYx3)/2._dp -            & 
&  (9*TpTYb3CTYb3*TrYx3adjYx3)/5._dp - 9*TpTYw3CTYw3*TrYx3adjYx3 - (9*TpYb3CYb3ml2*TrYx3adjYx3)/10._dp -& 
&  (9*TpYb3mHb32CYb3*TrYx3adjYx3)/5._dp - 18*mHu2*TpYw3CYw3*TrYx3adjYx3 - (9*TpYw3CYw3ml2*TrYx3adjYx3)/2._dp -& 
&  9*TpYw3mHw32CYw3*TrYx3adjYx3 + (-72._dp*(adjTYb3TYb3TpYb3CYb3) - 90._dp*(adjTYb3TYb3TpYw3CYw3) -& 
&  800._dp*(adjTYeTYeTpYeCYe) - 90._dp*(adjTYw3TYw3TpYb3CYb3) - 450._dp*(adjTYw3TYw3TpYw3CYw3) -& 
&  36._dp*(ml2TpYb3CYb3TpYb3CYb3) - 45._dp*(ml2TpYb3CYb3TpYw3CYw3) - 400._dp*(ml2TpYeCYeTpYeCYe) -& 
&  45._dp*(ml2TpYw3CYw3TpYb3CYb3) - 225._dp*(ml2TpYw3CYw3TpYw3CYw3) - 72._dp*(TpTYb3CYb3TpYb3CTYb3) -& 
&  90._dp*(TpTYb3CYb3TpYw3CTYw3) - 800._dp*(TpTYeCYeTpYeCTYe) - 90._dp*(TpTYw3CYw3TpYb3CTYb3) -& 
&  450._dp*(TpTYw3CYw3TpYw3CTYw3) - 72._dp*(TpYb3CTYb3TpTYb3CYb3) - 90._dp*(TpYb3CTYb3TpTYw3CYw3) -& 
&  72._dp*(TpYb3CYb3ml2TpYb3CYb3) - 90._dp*(TpYb3CYb3ml2TpYw3CYw3) - 72._dp*(TpYb3CYb3TpTYb3CTYb3) -& 
&  90._dp*(TpYb3CYb3TpTYw3CTYw3) - 36._dp*(TpYb3CYb3TpYb3CYb3ml2) - 72._dp*(TpYb3CYb3TpYb3mHb32CYb3) -& 
&  72._dp*(TpYb3mHb32CYb3TpYb3CYb3) - 90._dp*(TpYb3mHb32CYb3TpYw3CYw3) - 800._dp*(TpYeCTYeTpTYeCYe) -& 
&  800._dp*(TpYeCYeTpTYeCTYe) - 90._dp*(TpYw3CTYw3TpTYb3CYb3) - 450._dp*(TpYw3CTYw3TpTYw3CYw3) -& 
&  90._dp*(TpYw3CYw3TpTYb3CTYb3) - 450._dp*(TpYw3CYw3TpTYw3CTYw3) - 5*(9._dp*(TpYb3CYb3TpYw3CYw3ml2) +& 
&  18._dp*(TpYb3CYb3TpYw3mHw32CYw3) + 160._dp*(TpYeCYeml2TpYeCYe) + 80._dp*(TpYeCYeTpYeCYeml2) +& 
&  160._dp*(TpYeCYeTpYeme2CYe) + 160._dp*(TpYeme2CYeTpYeCYe) + 9*(2._dp*(TpYw3CYw3ml2TpYb3CYb3) +& 
&  10._dp*(TpYw3CYw3ml2TpYw3CYw3) + TpYw3CYw3TpYb3CYb3ml2 + 2._dp*(TpYw3CYw3TpYb3mHb32CYb3) +& 
&  5._dp*(TpYw3CYw3TpYw3CYw3ml2) + 10._dp*(TpYw3CYw3TpYw3mHw32CYw3) + 2._dp*(TpYw3mHw32CYw3TpYb3CYb3) +& 
&  10._dp*(TpYw3mHw32CYw3TpYw3CYw3))) - 36*TpYb3CTYb3*(TradjYb3TYb3 + 5*(2._dp*(TradjYuTYu) +& 
&  TradjYw3TYw3 + 2._dp*(TradjYx3TYx3))) - 400*TpTYeCYe*(3._dp*(TrCTYdTpYd) +            & 
&  TrCTYeTpYe) - 36*(TpTYb3CYb3 + 5._dp*(TpTYw3CYw3))*(TrCTYb3TpYb3 + 5*(2._dp*(TrCTYuTpYu) +& 
&  TrCTYw3TpYw3 + 2._dp*(TrCTYx3TpYx3))) - 36*TpYb3CYb3*(TrCTYb3TpTYb3 + 10._dp*(TrCTYuTpTYu) +& 
&  5._dp*(TrCTYw3TpTYw3) + 10._dp*(TrCTYx3TpTYx3) + 10._dp*(Trmd2adjYx3Yx3) +            & 
&  TrmHb32Yb3adjYb3 + 5._dp*(TrmHw32Yw3adjYw3) + 10._dp*(TrmHxb32Yx3adjYx3) +            & 
&  Trml2adjYb3Yb3 + 5*(Trml2adjYw3Yw3 + 2*(Trmq2adjYuYu + Trmu2YuadjYu)) +               & 
&  2*mHu2*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) + TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))))/200._dp  
betaml22 =  betaml22+ (3*g1p2*Conjg(MassB)*(-20._dp*(TpTYeCYe) + 40*MassB*TpYeCYe + 3*id3R*(69*g1p2*MassB + & 
&  5*g2p2*(2._dp*(MassB) + MassWB) + 25*g1p2*MassB*(NGHx3 +               & 
&  NGHxb3))))/25._dp + (3*g2p2*Conjg(MassWB)*(-20._dp*(TpTYw3CYw3) +      & 
&  40*MassWB*TpYw3CYw3 + id3R*(55*g2p2*MassWB + 3*g1p2*(MassB + 2._dp*(MassWB)) +        & 
&  15*g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))))/5._dp +& 
&  (6*g1p4*id3R*Tr2(1))/5._dp + 6*g2p4*id3R*Tr2(2) - 4*g1p2*id3R*Tr3(1)

 
Dml2 = oo16pi2*( betaml21 + oo16pi2 * betaml22 ) 

 
Else 
Dml2 = oo16pi2* betaml21 
End If 
 
 
Forall(i1=1:3) Dml2(i1,i1) =  Real(Dml2(i1,i1),dp)
!-------------------- 
! mHd2 
!-------------------- 
 
betamHd21  = 2*(-3*AbsMassWB*g2p2 + 3._dp*(TrCTYdTpTYd) + TrCTYeTpTYe +               & 
&  3._dp*(Trmd2YdadjYd) + Trme2YeadjYe + Trml2adjYeYe + 3._dp*(Trmq2adjYdYd)             & 
&  + 3*mHd2*TrYdadjYd + mHd2*TrYeadjYe) - (g1p2*(6._dp*(AbsMassB) + 5*Tr1(1)))/5._dp

 
 
If (TwoLoopRGE) Then 
betamHd22 = (g1p2*Conjg(MassB)*(9*(69*g1p2*MassB + 5*g2p2*(2._dp*(MassB) + MassWB)) +             & 
&  20._dp*(TradjYdTYd) - 60._dp*(TradjYeTYe) + 5*MassB*(-8*(TrYdadjYd - 3._dp*(TrYeadjYe)) +& 
&  45*g1p2*(NGHx3 + NGHxb3))) + 5*(-4*g1p2*TrCTYdTpTYd +   & 
&  160*g3p2*TrCTYdTpTYd + 4*g1p2*MassB*TrCTYdTpYd - 160*g3p2*MassG*TrCTYdTpYd +          & 
&  12*g1p2*TrCTYeTpTYe - 12*g1p2*MassB*TrCTYeTpYe - 4*g1p2*Trmd2YdadjYd + 160*g3p2*Trmd2YdadjYd +& 
&  12*g1p2*Trme2YeadjYe + 12*g1p2*Trml2adjYeYe - 4*g1p2*Trmq2adjYdYd + 160*g3p2*Trmq2adjYdYd -& 
&  3*mHd2*TrYb3adjYeYeadjYb3 - 3*mHu2*TrYb3adjYeYeadjYb3 - 4*(-80*AbsMassG*g3p2 +        & 
&  (g1p2 - 40._dp*(g3p2))*mHd2)*TrYdadjYd - 60*mHd2*TrYdadjYdTpYx3CYx3 - 60*mHu2*TrYdadjYdTpYx3CYx3 -& 
&  180*mHd2*TrYdadjYdYdadjYd - 30*mHd2*TrYdadjYuYuadjYd - 30*mHu2*TrYdadjYuYuadjYd +     & 
&  12*g1p2*mHd2*TrYeadjYe - 60*mHd2*TrYeadjYeYeadjYe - 15*(mHd2 + mHu2)*TrYeadjYw3Yw3adjYe -& 
&  3*(20._dp*(TradjYdTpYx3CTYx3TYd) + 20._dp*(TradjYx3TYx3CTYdTpYd) + 20._dp*(Trmd2TpYx3CYx3YdadjYd) +& 
&  20._dp*(Trmd2YdadjYdTpYx3CYx3) + 60._dp*(Trmd2YdadjYdYdadjYd) + 10._dp*(Trmd2YdadjYuYuadjYd) +& 
&  Trme2YeadjYb3Yb3adjYe + 20._dp*(Trme2YeadjYeYeadjYe) + 5._dp*(Trme2YeadjYw3Yw3adjYe) +& 
&  TrmHb32Yb3adjYeYeadjYb3 + 5._dp*(TrmHw32Yw3adjYeYeadjYw3) + 20._dp*(TrmHxb32Yx3CYdTpYdadjYx3) +& 
&  Trml2adjYb3Yb3adjYeYe + Trml2adjYeYeadjYb3Yb3 + 5*(4._dp*(Trml2adjYeYeadjYeYe) +      & 
&  Trml2adjYeYeadjYw3Yw3 + Trml2adjYw3Yw3adjYeYe + 2*(2._dp*(Trmq2adjYdTpYx3CYx3Yd) +    & 
&  6._dp*(Trmq2adjYdYdadjYdYd) + Trmq2adjYdYdadjYuYu + Trmq2adjYuYuadjYdYd +             & 
&  Trmu2YuadjYdYdadjYu)) + TrYb3adjYeTYeadjTYb3 + TrYb3TpTYeCTYeadjYb3 + 30._dp*(TrYdadjTYdTYdadjYd) +& 
&  10._dp*(TrYdadjTYuTYuadjYd) + 20._dp*(TrYdadjYdTpTYx3CTYx3) + 60._dp*(TrYdadjYdTYdadjTYd) +& 
&  10._dp*(TrYdadjYuTYuadjTYd) + 30._dp*(TrYdTpTYdCTYdadjYd) + TrYeadjTYb3TYb3adjYe +    & 
&  10._dp*(TrYeadjTYeTYeadjYe) + 5._dp*(TrYeadjTYw3TYw3adjYe) + TrYeadjYb3TYb3adjTYe +   & 
&  20._dp*(TrYeadjYeTYeadjTYe) + 5._dp*(TrYeadjYw3TYw3adjTYe) + 10._dp*(TrYeTpTYeCTYeadjYe) +& 
&  10._dp*(TrYuadjYdTYdadjTYu) + 10._dp*(TrYuTpTYdCTYdadjYu) + 5._dp*(TrYw3adjYeTYeadjTYw3) +& 
&  5._dp*(TrYw3TpTYeCTYeadjYw3) + 20._dp*(TrYx3CTYdTpTYdadjYx3)) - 160*g3p2*TradjYdTYd*Conjg(MassG) +& 
&  3*g2p2*Conjg(MassWB)*(55*g2p2*MassWB + 3*g1p2*(MassB + 2._dp*(MassWB)) +              & 
&  15*g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 + NGHxb3))) +& 
&  6*g1p4*Tr2(1) + 30*g2p4*Tr2(2) - 20*g1p2*Tr3(1)))/25._dp

 
DmHd2 = oo16pi2*( betamHd21 + oo16pi2 * betamHd22 ) 

 
Else 
DmHd2 = oo16pi2* betamHd21 
End If 
 
 
!-------------------- 
! mHu2 
!-------------------- 
 
betamHu21  = (-6*AbsMassB*g1p2)/5._dp + (3*(-10*AbsMassWB*g2p2 + TrCTYb3TpTYb3 +      & 
&  10._dp*(TrCTYuTpTYu) + 5._dp*(TrCTYw3TpTYw3) + 10._dp*(TrCTYx3TpTYx3) +               & 
&  10._dp*(Trmd2adjYx3Yx3) + TrmHb32Yb3adjYb3 + 5._dp*(TrmHw32Yw3adjYw3) +               & 
&  10._dp*(TrmHxb32Yx3adjYx3) + Trml2adjYb3Yb3 + 5*(Trml2adjYw3Yw3 + 2*(Trmq2adjYuYu +   & 
&  Trmu2YuadjYu)) + mHu2*TrYb3adjYb3 + 5*mHu2*(2._dp*(TrYuadjYu) + TrYw3adjYw3 +         & 
&  2._dp*(TrYx3adjYx3))))/5._dp + g1p2*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamHu22 = (8*g1p2*TrCTYuTpTYu)/5._dp + 32*g3p2*TrCTYuTpTYu - (8*g1p2*MassB*TrCTYuTpYu)/5._dp -  & 
&  32*g3p2*MassG*TrCTYuTpYu + 12*g2p2*TrCTYw3TpTYw3 - 12*g2p2*MassWB*TrCTYw3TpYw3 +      & 
&  4*g1p2*TrCTYx3TpTYx3 + 32*g3p2*TrCTYx3TpTYx3 - 4*g1p2*MassB*TrCTYx3TpYx3 -            & 
&  32*g3p2*MassG*TrCTYx3TpYx3 + 4*g1p2*Trmd2adjYx3Yx3 + 32*g3p2*Trmd2adjYx3Yx3 +         & 
&  12*g2p2*TrmHw32Yw3adjYw3 + 4*g1p2*TrmHxb32Yx3adjYx3 + 32*g3p2*TrmHxb32Yx3adjYx3 +     & 
&  12*g2p2*Trml2adjYw3Yw3 + (8*g1p2*Trmq2adjYuYu)/5._dp + 32*g3p2*Trmq2adjYuYu +         & 
&  (8*g1p2*Trmu2YuadjYu)/5._dp + 32*g3p2*Trmu2YuadjYu - (27*mHu2*TrYb3adjYb3Yb3adjYb3)/25._dp -& 
&  (3*mHu2*TrYb3adjYeYeadjYb3)/5._dp - (27*mHu2*TrYb3adjYw3Yw3adjYb3)/5._dp -            & 
&  12*mHu2*TrYdadjYdTpYx3CYx3 - 6*mHu2*TrYdadjYuYuadjYd - 3*mHu2*TrYeadjYw3Yw3adjYe +    & 
&  64*AbsMassG*g3p2*TrYuadjYu + (8*g1p2*mHu2*TrYuadjYu)/5._dp + 32*g3p2*mHu2*TrYuadjYu - & 
&  36*mHu2*TrYuadjYuYuadjYu + 12*g2p2*mHu2*TrYw3adjYw3 - (27*mHu2*TrYw3adjYw3Yw3adjYw3)/2._dp +& 
&  64*AbsMassG*g3p2*TrYx3adjYx3 + 4*g1p2*mHu2*TrYx3adjYx3 + 32*g3p2*mHu2*TrYx3adjYx3 -   & 
&  36*mHu2*TrYx3adjYx3Yx3adjYx3 - (3*(18._dp*(TrYb3adjTYb3TYb3adjYb3) + 20._dp*(TrYb3adjTYeTYeadjYb3) +& 
&  45._dp*(TrYb3adjTYw3TYw3adjYb3) + 36._dp*(TrYb3adjYb3TYb3adjTYb3) + 20._dp*(TrYb3adjYeTYeadjTYb3) +& 
&  90._dp*(TrYb3adjYw3TYw3adjTYb3) + 18._dp*(TrYb3TpTYb3CTYb3adjYb3) + 45._dp*(TrYb3TpTYw3CTYw3adjYb3) +& 
&  200._dp*(TrYdadjYuTYuadjTYd) + 200._dp*(TrYdTpTYuCTYuadjYd) + 20._dp*(TrYeadjYb3TYb3adjTYe) +& 
&  100._dp*(TrYeadjYw3TYw3adjTYe) + 20*mHd2*(TrYb3adjYeYeadjYb3 + 5*(4._dp*(TrYdadjYdTpYx3CYx3) +& 
&  2._dp*(TrYdadjYuYuadjYd) + TrYeadjYw3Yw3adjYe)) + 20._dp*(TrYeTpTYb3CTYb3adjYe) +     & 
&  100._dp*(TrYeTpTYw3CTYw3adjYe) + 200._dp*(TrYuadjTYdTYdadjYu) + 600._dp*(TrYuadjTYuTYuadjYu) +& 
&  200._dp*(TrYuadjYdTYdadjTYu) + 1200._dp*(TrYuadjYuTYuadjTYu) + 600._dp*(TrYuTpTYuCTYuadjYu) +& 
&  45._dp*(TrYw3adjTYb3TYb3adjYw3) + 100._dp*(TrYw3adjTYeTYeadjYw3) + 225._dp*(TrYw3adjTYw3TYw3adjYw3) +& 
&  90._dp*(TrYw3adjYb3TYb3adjTYw3) + 100._dp*(TrYw3adjYeTYeadjTYw3) + 450._dp*(TrYw3adjYw3TYw3adjTYw3) +& 
&  45._dp*(TrYw3TpTYb3CTYb3adjYw3) + 225._dp*(TrYw3TpTYw3CTYw3adjYw3) + 2*(200._dp*(TradjYdTpTYx3CTYx3TpYd) +& 
&  200._dp*(TradjYdTpYx3CTYx3TpTYd) + 200._dp*(TradjYx3TYx3CTYdTpYd) + 600._dp*(Trmd2adjYx3Yx3adjYx3Yx3) +& 
&  200._dp*(Trmd2TpYx3CYx3YdadjYd) + 200._dp*(Trmd2YdadjYdTpYx3CYx3) + 100._dp*(Trmd2YdadjYuYuadjYd) +& 
&  10._dp*(Trme2YeadjYb3Yb3adjYe) + 50._dp*(Trme2YeadjYw3Yw3adjYe) + 18._dp*(TrmHb32Yb3adjYb3Yb3adjYb3) +& 
&  10._dp*(TrmHb32Yb3adjYeYeadjYb3) + 45._dp*(TrmHb32Yb3adjYw3Yw3adjYb3) +               & 
&  45._dp*(TrmHw32Yw3adjYb3Yb3adjYw3) + 50._dp*(TrmHw32Yw3adjYeYeadjYw3) +               & 
&  225._dp*(TrmHw32Yw3adjYw3Yw3adjYw3) + 200._dp*(TrmHxb32CYx3YdadjYdTpYx3) +            & 
&  600._dp*(TrmHxb32Yx3adjYx3Yx3adjYx3) + 18._dp*(Trml2adjYb3Yb3adjYb3Yb3) +             & 
&  5*(2._dp*(Trml2adjYb3Yb3adjYeYe) + 9._dp*(Trml2adjYb3Yb3adjYw3Yw3) + 2._dp*(Trml2adjYeYeadjYb3Yb3) +& 
&  10._dp*(Trml2adjYeYeadjYw3Yw3) + 9._dp*(Trml2adjYw3Yw3adjYb3Yb3) + 5*(2._dp*(Trml2adjYw3Yw3adjYeYe) +& 
&  9._dp*(Trml2adjYw3Yw3adjYw3Yw3) + 4*(2._dp*(Trmq2adjYdTpYx3CYx3Yd) + Trmq2adjYdYdadjYuYu +& 
&  Trmq2adjYuYuadjYdYd + 6._dp*(Trmq2adjYuYuadjYuYu) + Trmu2YuadjYdYdadjYu +             & 
&  6._dp*(Trmu2YuadjYuYuadjYu)))) + 300._dp*(TrYx3adjTYx3TYx3adjYx3) + 600._dp*(TrYx3adjYx3TYx3adjTYx3) +& 
&  200._dp*(TrYx3CTYdTpTYdadjYx3) + 300._dp*(TrYx3TpTYx3CTYx3adjYx3))))/100._dp  
betamHu22 =  betamHu22- 32*g3p2*TradjYuTYu*Conjg(MassG) - 32*g3p2*TradjYx3TYx3*Conjg(MassG) &
& + (g1p2*Conjg(MassB)*(9*(69*g1p2*MassB +& 
&  5*g2p2*(2._dp*(MassB) + MassWB)) - 40._dp*(TradjYuTYu) - 100._dp*(TradjYx3TYx3) +     & 
&  5*MassB*(8*(2._dp*(TrYuadjYu) + 5._dp*(TrYx3adjYx3)) + 45*g1p2*(NGHx3 +& 
&  NGHxb3))))/25._dp + (3*g2p2*Conjg(MassWB)*(55*g2p2*MassWB +            & 
&  3*g1p2*(MassB + 2._dp*(MassWB)) - 20._dp*(TradjYw3TYw3) + 5*MassWB*(8._dp*(TrYw3adjYw3) +& 
&  3*g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))))/5._dp +& 
&  (6*g1p4*Tr2(1))/5._dp + 6*g2p4*Tr2(2) + 4*g1p2*Tr3(1)

 
DmHu2 = oo16pi2*( betamHu21 + oo16pi2 * betamHu22 ) 

 
Else 
DmHu2 = oo16pi2* betamHu21 
End If 
 
 
!-------------------- 
! md2 
!-------------------- 
 
betamd21  = 2*(md2TpYx3CYx3 + md2YdadjYd + 2._dp*(TpTYx3CTYx3) + 2*mHu2*TpYx3CYx3 +   & 
&  TpYx3CYx3md2 + 2._dp*(TpYx3mHxb32CYx3) + 2._dp*(TYdadjTYd) + 2*mHd2*YdadjYd +         & 
&  YdadjYdmd2 + 2._dp*(Ydmq2adjYd)) - (2*id3R*(80*AbsMassG*g3p2 + g1p2*(4._dp*(AbsMassB) & 
&  - 5*Tr1(1))))/15._dp

 
 
If (TwoLoopRGE) Then 
betamd22 = 2*g1p2*md2TpYx3CYx3 + 6*g2p2*md2TpYx3CYx3 + (2*g1p2*md2YdadjYd)/5._dp +               & 
&  6*g2p2*md2YdadjYd + 4*g1p2*TpTYx3CTYx3 + 12*g2p2*TpTYx3CTYx3 - 4*g1p2*MassB*TpYx3CTYx3 -& 
&  12*g2p2*MassWB*TpYx3CTYx3 + 24*AbsMassWB*g2p2*TpYx3CYx3 + 4*g1p2*mHu2*TpYx3CYx3 +     & 
&  12*g2p2*mHu2*TpYx3CYx3 + 2*g1p2*TpYx3CYx3md2 + 6*g2p2*TpYx3CYx3md2 - 8*mHu2*TpYx3CYx3TpYx3CYx3 +& 
&  4*g1p2*TpYx3mHxb32CYx3 + 12*g2p2*TpYx3mHxb32CYx3 - (6*TpYx3CTYx3*TradjYb3TYb3)/5._dp -& 
&  12*TpYx3CTYx3*TradjYuTYu - 6*TpYx3CTYx3*TradjYw3TYw3 - 12*TpYx3CTYx3*TradjYx3TYx3 -   & 
&  (6*TpYx3CYx3*TrCTYb3TpTYb3)/5._dp - (6*TpTYx3CYx3*TrCTYb3TpYb3)/5._dp -               & 
&  12*TpYx3CYx3*TrCTYuTpTYu - 12*TpTYx3CYx3*TrCTYuTpYu - 6*TpYx3CYx3*TrCTYw3TpTYw3 -     & 
&  6*TpTYx3CYx3*TrCTYw3TpYw3 - 12*TpYx3CYx3*TrCTYx3TpTYx3 - 12*TpTYx3CYx3*TrCTYx3TpYx3 - & 
&  12*TpYx3CYx3*Trmd2adjYx3Yx3 - (6*TpYx3CYx3*TrmHb32Yb3adjYb3)/5._dp - 6*TpYx3CYx3*TrmHw32Yw3adjYw3 -& 
&  12*TpYx3CYx3*TrmHxb32Yx3adjYx3 - (6*TpYx3CYx3*Trml2adjYb3Yb3)/5._dp - 6*TpYx3CYx3*Trml2adjYw3Yw3 -& 
&  12*TpYx3CYx3*Trmq2adjYuYu - 12*TpYx3CYx3*Trmu2YuadjYu - (3*md2TpYx3CYx3*TrYb3adjYb3)/5._dp -& 
&  (6*TpTYx3CTYx3*TrYb3adjYb3)/5._dp - (12*mHu2*TpYx3CYx3*TrYb3adjYb3)/5._dp -           & 
&  (3*TpYx3CYx3md2*TrYb3adjYb3)/5._dp - (6*TpYx3mHxb32CYx3*TrYb3adjYb3)/5._dp -          & 
&  6*md2YdadjYd*TrYdadjYd - 2*md2YdadjYd*TrYeadjYe - 6*md2TpYx3CYx3*TrYuadjYu -          & 
&  12*TpTYx3CTYx3*TrYuadjYu - 24*mHu2*TpYx3CYx3*TrYuadjYu - 6*TpYx3CYx3md2*TrYuadjYu -   & 
&  12*TpYx3mHxb32CYx3*TrYuadjYu - 3*md2TpYx3CYx3*TrYw3adjYw3 - 6*TpTYx3CTYx3*TrYw3adjYw3 -& 
&  12*mHu2*TpYx3CYx3*TrYw3adjYw3 - 3*TpYx3CYx3md2*TrYw3adjYw3 - 6*TpYx3mHxb32CYx3*TrYw3adjYw3 -& 
&  6*md2TpYx3CYx3*TrYx3adjYx3 - 12*TpTYx3CTYx3*TrYx3adjYx3 - 24*mHu2*TpYx3CYx3*TrYx3adjYx3 -& 
&  6*TpYx3CYx3md2*TrYx3adjYx3 - 12*TpYx3mHxb32CYx3*TrYx3adjYx3 + (4*g1p2*TYdadjTYd)/5._dp +& 
&  12*g2p2*TYdadjTYd - 12*TrYdadjYd*TYdadjTYd - 4*TrYeadjYe*TYdadjTYd - 12*TrCTYdTpYd*TYdadjYd  
betamd22 =  betamd22- 4*TrCTYeTpYe*TYdadjYd - (4*g1p2*MassB*YdadjTYd)/5._dp - 12*g2p2*MassWB*YdadjTYd -     & 
&  12*TradjYdTYd*YdadjTYd - 4*TradjYeTYe*YdadjTYd + 24*AbsMassWB*g2p2*YdadjYd +          & 
&  (4*g1p2*mHd2*YdadjYd)/5._dp + 12*g2p2*mHd2*YdadjYd - 12*TrCTYdTpTYd*YdadjYd -         & 
&  4*TrCTYeTpTYe*YdadjYd - 12*Trmd2YdadjYd*YdadjYd - 4*Trme2YeadjYe*YdadjYd -            & 
&  4*Trml2adjYeYe*YdadjYd - 12*Trmq2adjYdYd*YdadjYd - 24*mHd2*TrYdadjYd*YdadjYd -        & 
&  8*mHd2*TrYeadjYe*YdadjYd + (2*g1p2*YdadjYdmd2)/5._dp + 6*g2p2*YdadjYdmd2 -            & 
&  6*TrYdadjYd*YdadjYdmd2 - 2*TrYeadjYe*YdadjYdmd2 - 8*mHd2*YdadjYdYdadjYd -             & 
&  4*mHd2*YdadjYuYuadjYd - 4*mHu2*YdadjYuYuadjYd + (4*g1p2*Ydmq2adjYd)/5._dp +           & 
&  12*g2p2*Ydmq2adjYd - 12*TrYdadjYd*Ydmq2adjYd - 4*TrYeadjYe*Ydmq2adjYd -               & 
&  2*(2._dp*(adjTYx3TYx3TpYx3CYx3) + md2TpYx3CYx3TpYx3CYx3 + md2YdadjYdYdadjYd +         & 
&  md2YdadjYuYuadjYd + 2._dp*(TpTYx3CYx3TpYx3CTYx3) + 2._dp*(TpYx3CTYx3TpTYx3CYx3) +     & 
&  2._dp*(TpYx3CYx3TpTYx3CTYx3) + TpYx3CYx3TpYx3CYx3md2 + 2*(TpYx3CYx3md2TpYx3CYx3 +     & 
&  TpYx3CYx3TpYx3mHxb32CYx3 + TpYx3mHxb32CYx3TpYx3CYx3) + 2._dp*(TYdadjYdYdadjTYd) +     & 
&  2._dp*(TYdadjYuYuadjTYd) + 2._dp*(TYdTpYdCTYdadjYd) + 2._dp*(TYdTpYuCTYuadjYd) +      & 
&  2._dp*(YdadjTYdTYdadjYd) + 2._dp*(YdadjTYuTYuadjYd) + 2._dp*(YdadjYdmd2YdadjYd) +     & 
&  2._dp*(YdadjYdTYdadjTYd) + YdadjYdYdadjYdmd2 + 2._dp*(YdadjYdYdmq2adjYd) +            & 
&  2._dp*(YdadjYumu2YuadjYd) + 2._dp*(YdadjYuTYuadjTYd) + YdadjYuYuadjYdmd2 +            & 
&  2._dp*(YdadjYuYumq2adjYd) + 2._dp*(Ydmq2adjYdYdadjYd) + 2._dp*(Ydmq2adjYuYuadjYd)) -  & 
&  12*g2p2*TpTYx3CYx3*Conjg(MassWB) - 12*g2p2*TYdadjYd*Conjg(MassWB) + (32*g3p2*id3R*Conjg(MassG)*(-& 
& 60*g3p2*MassG + 2*g1p2*(MassB + 2._dp*(MassG)) + 45*g3p2*MassG*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)))/45._dp + (4*g1p2*Conjg(MassB)*(-      & 
& 45*(5._dp*(TpTYx3CYx3) + TYdadjYd - 2*MassB*(5._dp*(TpYx3CYx3) + YdadjYd)) +           & 
&  id3R*(606*g1p2*MassB + 80*g3p2*(2._dp*(MassB) + MassG) + 225*g1p2*MassB*(NGHx3 +& 
&  NGHxb3))))/225._dp + (8*g1p4*id3R*Tr2(1))/15._dp + (32*g3p4*id3R*Tr2(3))/3._dp  
betamd22 =  betamd22+ (8*g1p2*id3R*Tr3(1))/3._dp

 
Dmd2 = oo16pi2*( betamd21 + oo16pi2 * betamd22 ) 

 
Else 
Dmd2 = oo16pi2* betamd21 
End If 
 
 
Forall(i1=1:3) Dmd2(i1,i1) =  Real(Dmd2(i1,i1),dp) 
!-------------------- 
! mu2 
!-------------------- 
 
betamu21  = 2*(mu2YuadjYu + 2._dp*(TYuadjTYu) + 2*mHu2*YuadjYu + YuadjYumu2 +         & 
&  2._dp*(Yumq2adjYu)) - (4*id3R*(40*AbsMassG*g3p2 + g1p2*(8._dp*(AbsMassB)              & 
&  + 5*Tr1(1))))/15._dp

 
 
If (TwoLoopRGE) Then 
betamu22 = (-2*g1p2*mu2YuadjYu)/5._dp + 6*g2p2*mu2YuadjYu - (3*mu2YuadjYu*TrYb3adjYb3)/5._dp -   & 
&  6*mu2YuadjYu*TrYuadjYu - 3*mu2YuadjYu*TrYw3adjYw3 - 6*mu2YuadjYu*TrYx3adjYx3 -        & 
&  (4*g1p2*TYuadjTYu)/5._dp + 12*g2p2*TYuadjTYu - (6*TrYb3adjYb3*TYuadjTYu)/5._dp -      & 
&  12*TrYuadjYu*TYuadjTYu - 6*TrYw3adjYw3*TYuadjTYu - 12*TrYx3adjYx3*TYuadjTYu -         & 
&  (6*TrCTYb3TpYb3*TYuadjYu)/5._dp - 12*TrCTYuTpYu*TYuadjYu - 6*TrCTYw3TpYw3*TYuadjYu -  & 
&  12*TrCTYx3TpYx3*TYuadjYu + (4*g1p2*MassB*YuadjTYu)/5._dp - 12*g2p2*MassWB*YuadjTYu -  & 
&  (6*TradjYb3TYb3*YuadjTYu)/5._dp - 12*TradjYuTYu*YuadjTYu - 6*TradjYw3TYw3*YuadjTYu -  & 
&  12*TradjYx3TYx3*YuadjTYu - 4*mHu2*YuadjYdYdadjYu + 24*AbsMassWB*g2p2*YuadjYu -        & 
&  (4*g1p2*mHu2*YuadjYu)/5._dp + 12*g2p2*mHu2*YuadjYu - (6*TrCTYb3TpTYb3*YuadjYu)/5._dp -& 
&  12*TrCTYuTpTYu*YuadjYu - 6*TrCTYw3TpTYw3*YuadjYu - 12*TrCTYx3TpTYx3*YuadjYu -         & 
&  12*Trmd2adjYx3Yx3*YuadjYu - (6*TrmHb32Yb3adjYb3*YuadjYu)/5._dp - 6*TrmHw32Yw3adjYw3*YuadjYu -& 
&  12*TrmHxb32Yx3adjYx3*YuadjYu - (6*Trml2adjYb3Yb3*YuadjYu)/5._dp - 6*Trml2adjYw3Yw3*YuadjYu -& 
&  12*Trmq2adjYuYu*YuadjYu - 12*Trmu2YuadjYu*YuadjYu - (12*mHu2*TrYb3adjYb3*YuadjYu)/5._dp -& 
&  24*mHu2*TrYuadjYu*YuadjYu - 12*mHu2*TrYw3adjYw3*YuadjYu - 24*mHu2*TrYx3adjYx3*YuadjYu -& 
&  (2*g1p2*YuadjYumu2)/5._dp + 6*g2p2*YuadjYumu2 - (3*TrYb3adjYb3*YuadjYumu2)/5._dp -    & 
&  6*TrYuadjYu*YuadjYumu2 - 3*TrYw3adjYw3*YuadjYumu2 - 6*TrYx3adjYx3*YuadjYumu2 -        & 
&  8*mHu2*YuadjYuYuadjYu - (4*g1p2*Yumq2adjYu)/5._dp + 12*g2p2*Yumq2adjYu -              & 
&  (6*TrYb3adjYb3*Yumq2adjYu)/5._dp - 12*TrYuadjYu*Yumq2adjYu - 6*TrYw3adjYw3*Yumq2adjYu -& 
&  12*TrYx3adjYx3*Yumq2adjYu - 2*(mu2YuadjYdYdadjYu + mu2YuadjYuYuadjYu + 2._dp*(TYuadjYdYdadjTYu) +& 
&  2._dp*(TYuadjYuYuadjTYu) + 2._dp*(TYuTpYdCTYdadjYu) + 2._dp*(TYuTpYuCTYuadjYu) +      & 
&  2._dp*(YuadjTYdTYdadjYu) + 2._dp*(YuadjTYuTYuadjYu) + 2._dp*(YuadjYdmd2YdadjYu) +     & 
&  2._dp*(YuadjYdTYdadjTYu) + 2*mHd2*YuadjYdYdadjYu + YuadjYdYdadjYumu2 + 2._dp*(YuadjYdYdmq2adjYu) +& 
&  2._dp*(YuadjYuTYuadjTYu) + YuadjYuYuadjYumu2 + 2*(YuadjYumu2YuadjYu + YuadjYuYumq2adjYu) +& 
&  2._dp*(Yumq2adjYdYdadjYu) + 2._dp*(Yumq2adjYuYuadjYu)) - 12*g2p2*TYuadjYu*Conjg(MassWB)  
betamu22 =  betamu22+ (32*g3p2*id3R*Conjg(MassG)*(-60*g3p2*MassG + 8*g1p2*(MassB + 2._dp*(MassG)) +         & 
&  45*g3p2*MassG*(3*NGHg3 + NGHx3 + NGHxb3)))/45._dp +& 
&  (4*g1p2*Conjg(MassB)*(45*(TYuadjYu - 2*MassB*YuadjYu) + 4*id3R*(642*g1p2*MassB +      & 
&  80*g3p2*(2._dp*(MassB) + MassG) + 225*g1p2*MassB*(NGHx3 +              & 
&  NGHxb3))))/225._dp + (32*g1p4*id3R*Tr2(1))/15._dp + (32*g3p4*id3R*Tr2(3))/3._dp -& 
&  (16*g1p2*id3R*Tr3(1))/3._dp

 
Dmu2 = oo16pi2*( betamu21 + oo16pi2 * betamu22 ) 

 
Else 
Dmu2 = oo16pi2* betamu21 
End If 
 
 
Forall(i1=1:3) Dmu2(i1,i1) =  Real(Dmu2(i1,i1),dp) 
!-------------------- 
! me2 
!-------------------- 
 
betame21  = 2*(me2YeadjYe + 2._dp*(TYeadjTYe) + 2*mHd2*YeadjYe + YeadjYeme2 +         & 
&  2._dp*(Yeml2adjYe)) + (2*g1p2*id3R*(-12._dp*(AbsMassB) + 5*Tr1(1)))/5._dp

 
 
If (TwoLoopRGE) Then 
betame22 = (12*g1p2*Conjg(MassB)*(5*(TYeadjYe - 2*MassB*YeadjYe) + 3*g1p2*id3R*MassB*(78 +       & 
&  25*NGHx3 + 25*NGHxb3)))/25._dp + (-3._dp*(me2YeadjYb3Yb3adjYe) -& 
&  10._dp*(me2YeadjYeYeadjYe) - 15._dp*(me2YeadjYw3Yw3adjYe) - 6._dp*(TYeadjYb3Yb3adjTYe) -& 
&  20._dp*(TYeadjYeYeadjTYe) - 30._dp*(TYeadjYw3Yw3adjTYe) - 6._dp*(TYeTpYb3CTYb3adjYe) -& 
&  20._dp*(TYeTpYeCTYeadjYe) - 30._dp*(TYeTpYw3CTYw3adjYe) - 6._dp*(YeadjTYb3TYb3adjYe) -& 
&  20._dp*(YeadjTYeTYeadjYe) - 30._dp*(YeadjTYw3TYw3adjYe) - 6._dp*(YeadjYb3mHb32Yb3adjYe) -& 
&  6._dp*(YeadjYb3TYb3adjTYe) - 3._dp*(YeadjYb3Yb3adjYeme2) - 6._dp*(YeadjYb3Yb3ml2adjYe) -& 
&  4*(-30*AbsMassWB*g2p2 + 3*g1p2*mHd2 - 15*g2p2*mHd2 + 15._dp*(TrCTYdTpTYd) +           & 
&  5._dp*(TrCTYeTpTYe) + 15._dp*(Trmd2YdadjYd) + 5._dp*(Trme2YeadjYe) + 5._dp*(Trml2adjYeYe) +& 
&  15._dp*(Trmq2adjYdYd) + 10*mHd2*(3._dp*(TrYdadjYd) + TrYeadjYe))*YeadjYe -            & 
&  20._dp*(YeadjYeTYeadjTYe) - 30._dp*(YeadjYw3TYw3adjTYe) - 2*(3*(mHd2 + mHu2)*YeadjYb3Yb3adjYe +& 
&  20*mHd2*YeadjYeYeadjYe + 15*(mHd2 + mHu2)*YeadjYw3Yw3adjYe) - 5*(4._dp*(YeadjYeme2YeadjYe) +& 
&  2._dp*(YeadjYeYeadjYeme2) + 4._dp*(YeadjYeYeml2adjYe) + 3._dp*(YeadjYw3Yw3adjYeme2) + & 
&  6*(YeadjYw3mHw32Yw3adjYe + YeadjYw3Yw3ml2adjYe)) - 6._dp*(Yeml2adjYb3Yb3adjYe) -      & 
&  20._dp*(Yeml2adjYeYeadjYe) - 30._dp*(Yeml2adjYw3Yw3adjYe) - 2*(10*(3._dp*(TrCTYdTpYd) +& 
&  TrCTYeTpYe)*TYeadjYe + (-6*g1p2*MassB + 30*g2p2*MassWB + 30._dp*(TradjYdTYd) +        & 
&  10._dp*(TradjYeTYe))*YeadjTYe + (3*(g1p2 - 5._dp*(g2p2) + 5._dp*(TrYdadjYd)) +        & 
&  5._dp*(TrYeadjYe))*(me2YeadjYe + 2._dp*(TYeadjTYe) + YeadjYeme2 + 2._dp*(Yeml2adjYe)) +& 
&  30*g2p2*TYeadjYe*Conjg(MassWB)) + 8*id3R*(3*g1p4*Tr2(1) + 5*g1p2*Tr3(1)))/5._dp

 
Dme2 = oo16pi2*( betame21 + oo16pi2 * betame22 ) 

 
Else 
Dme2 = oo16pi2* betame21 
End If 
 
 
Forall(i1=1:3) Dme2(i1,i1) =  Real(Dme2(i1,i1),dp) 
!-------------------- 
! mHw32 
!-------------------- 
 
betamHw321  = -16*AbsMassWB*g2p2*id3R + mHw32Yw3adjYw3 + 2._dp*(TYw3adjTYw3)          & 
&  + 2*mHu2*Yw3adjYw3 + Yw3adjYw3mHw32 + 2._dp*(Yw3ml2adjYw3)

 
 
If (TwoLoopRGE) Then 
betamHw322 = -3._dp*(mHw32Yw3adjYb3Yb3adjYw3)/10._dp - mHw32Yw3adjYeYeadjYw3 + (3*g1p2*mHw32Yw3adjYw3)/5._dp -& 
&  g2p2*mHw32Yw3adjYw3 - 3._dp*(mHw32Yw3adjYw3Yw3adjYw3)/2._dp - (3*mHw32Yw3adjYw3*TrYb3adjYb3)/10._dp -& 
&  3*mHw32Yw3adjYw3*TrYuadjYu - (3*mHw32Yw3adjYw3*TrYw3adjYw3)/2._dp - 3*mHw32Yw3adjYw3*TrYx3adjYx3 +& 
&  (6*g1p2*TYw3adjTYw3)/5._dp - 2*g2p2*TYw3adjTYw3 - (3*TrYb3adjYb3*TYw3adjTYw3)/5._dp - & 
&  6*TrYuadjYu*TYw3adjTYw3 - 3*TrYw3adjYw3*TYw3adjTYw3 - 6*TrYx3adjYx3*TYw3adjTYw3 -     & 
&  3._dp*(TYw3adjYb3Yb3adjTYw3)/5._dp - 2._dp*(TYw3adjYeYeadjTYw3) - (3*TrCTYb3TpYb3*TYw3adjYw3)/5._dp -& 
&  6*TrCTYuTpYu*TYw3adjYw3 - 3*TrCTYw3TpYw3*TYw3adjYw3 - 6*TrCTYx3TpYx3*TYw3adjYw3 -     & 
&  3._dp*(TYw3adjYw3Yw3adjTYw3) - 3._dp*(TYw3TpYb3CTYb3adjYw3)/5._dp - 2._dp*(TYw3TpYeCTYeadjYw3) -& 
&  3._dp*(TYw3TpYw3CTYw3adjYw3) - 3._dp*(Yw3adjTYb3TYb3adjYw3)/5._dp - 2._dp*(Yw3adjTYeTYeadjYw3) -& 
&  (6*g1p2*MassB*Yw3adjTYw3)/5._dp + 2*g2p2*MassWB*Yw3adjTYw3 - (3*TradjYb3TYb3*Yw3adjTYw3)/5._dp -& 
&  6*TradjYuTYu*Yw3adjTYw3 - 3*TradjYw3TYw3*Yw3adjTYw3 - 6*TradjYx3TYx3*Yw3adjTYw3 -     & 
&  3._dp*(Yw3adjTYw3TYw3adjYw3) - 3._dp*(Yw3adjYb3mHb32Yb3adjYw3)/5._dp - 3._dp*(Yw3adjYb3TYb3adjTYw3)/5._dp -& 
&  (6*mHu2*Yw3adjYb3Yb3adjYw3)/5._dp - 3._dp*(Yw3adjYb3Yb3adjYw3mHw32)/10._dp -          & 
&  3._dp*(Yw3adjYb3Yb3ml2adjYw3)/5._dp - 2._dp*(Yw3adjYeme2YeadjYw3) - 2._dp*(Yw3adjYeTYeadjTYw3) -& 
&  2*mHd2*Yw3adjYeYeadjYw3 - 2*mHu2*Yw3adjYeYeadjYw3 - Yw3adjYeYeadjYw3mHw32 -           & 
&  2._dp*(Yw3adjYeYeml2adjYw3) + (12*AbsMassB*g1p2*Yw3adjYw3)/5._dp + (6*g1p2*mHu2*Yw3adjYw3)/5._dp -& 
&  2*g2p2*mHu2*Yw3adjYw3 - (3*TrCTYb3TpTYb3*Yw3adjYw3)/5._dp - 6*TrCTYuTpTYu*Yw3adjYw3 - & 
&  3*TrCTYw3TpTYw3*Yw3adjYw3 - 6*TrCTYx3TpTYx3*Yw3adjYw3 - 6*Trmd2adjYx3Yx3*Yw3adjYw3 -  & 
&  (3*TrmHb32Yb3adjYb3*Yw3adjYw3)/5._dp - 3*TrmHw32Yw3adjYw3*Yw3adjYw3 - 6*TrmHxb32Yx3adjYx3*Yw3adjYw3 -& 
&  (3*Trml2adjYb3Yb3*Yw3adjYw3)/5._dp - 3*Trml2adjYw3Yw3*Yw3adjYw3 - 6*Trmq2adjYuYu*Yw3adjYw3 -& 
&  6*Trmu2YuadjYu*Yw3adjYw3 - (6*mHu2*TrYb3adjYb3*Yw3adjYw3)/5._dp - 12*mHu2*TrYuadjYu*Yw3adjYw3  
betamHw322 =  betamHw322- 6*mHu2*TrYw3adjYw3*Yw3adjYw3 - 12*mHu2*TrYx3adjYx3*Yw3adjYw3 + (3*g1p2*Yw3adjYw3mHw32)/5._dp -& 
&  g2p2*Yw3adjYw3mHw32 - (3*TrYb3adjYb3*Yw3adjYw3mHw32)/10._dp - 3*TrYuadjYu*Yw3adjYw3mHw32 -& 
&  (3*TrYw3adjYw3*Yw3adjYw3mHw32)/2._dp - 3*TrYx3adjYx3*Yw3adjYw3mHw32 - 3._dp*(Yw3adjYw3mHw32Yw3adjYw3) -& 
&  3._dp*(Yw3adjYw3TYw3adjTYw3) - 6*mHu2*Yw3adjYw3Yw3adjYw3 - 3._dp*(Yw3adjYw3Yw3adjYw3mHw32)/2._dp -& 
&  3._dp*(Yw3adjYw3Yw3ml2adjYw3) - 3._dp*(Yw3ml2adjYb3Yb3adjYw3)/5._dp - 2._dp*(Yw3ml2adjYeYeadjYw3) +& 
&  (6*g1p2*Yw3ml2adjYw3)/5._dp - 2*g2p2*Yw3ml2adjYw3 - (3*TrYb3adjYb3*Yw3ml2adjYw3)/5._dp -& 
&  6*TrYuadjYu*Yw3ml2adjYw3 - 3*TrYw3adjYw3*Yw3ml2adjYw3 - 6*TrYx3adjYx3*Yw3ml2adjYw3 -  & 
&  3._dp*(Yw3ml2adjYw3Yw3adjYw3) - (6*g1p2*TYw3adjYw3*Conjg(MassB))/5._dp +              & 
&  2*g2p2*Conjg(MassWB)*(TYw3adjYw3 - 2*MassWB*Yw3adjYw3 + 4*g2p2*id3R*MassWB*(26 +      & 
&  12*NGHw3 + 9*NGHx3 + 9*NGHxb3)) +        & 
&  16*g2p4*id3R*Tr2(2)

 
DmHw32 = oo16pi2*( betamHw321 + oo16pi2 * betamHw322 ) 

 
Else 
DmHw32 = oo16pi2* betamHw321 
End If 
 
 
Forall(i1=1:3) DmHw32(i1,i1) =  Real(DmHw32(i1,i1),dp) 
!-------------------- 
! mHg32 
!-------------------- 
 
betamHg321  = -24*AbsMassG*g3p2*id3R

 
 
If (TwoLoopRGE) Then 
betamHg322 = 24*g3p4*id3R*(3*AbsMassG*(2 + 3*NGHg3 + NGHx3 +         & 
&  NGHxb3) + Tr2(3))

 
DmHg32 = oo16pi2*( betamHg321 + oo16pi2 * betamHg322 ) 

 
Else 
DmHg32 = oo16pi2* betamHg321 
End If 
 
 
Forall(i1=1:3) DmHg32(i1,i1) =  Real(DmHg32(i1,i1),dp) 
!-------------------- 
! mHb32 
!-------------------- 
 
betamHb321  = (3*(mHb32Yb3adjYb3 + 2._dp*(TYb3adjTYb3) + 2*mHu2*Yb3adjYb3 +           & 
&  Yb3adjYb3mHb32 + 2._dp*(Yb3ml2adjYb3)))/5._dp

 
 
If (TwoLoopRGE) Then 
betamHb322 = (3*(-3._dp*(mHb32Yb3adjYb3Yb3adjYb3) - 10._dp*(mHb32Yb3adjYeYeadjYb3) -               & 
&  15._dp*(mHb32Yb3adjYw3Yw3adjYb3) - 3*mHb32Yb3adjYb3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) +& 
&  5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3))) - 6*(TrYb3adjYb3 +& 
&  10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))*TYb3adjTYb3 -& 
&  6._dp*(TYb3adjYb3Yb3adjTYb3) - 20._dp*(TYb3adjYeYeadjTYb3) - 30._dp*(TYb3adjYw3Yw3adjTYb3) -& 
&  6._dp*(TYb3TpYb3CTYb3adjYb3) - 20._dp*(TYb3TpYeCTYeadjYb3) - 30._dp*(TYb3TpYw3CTYw3adjYb3) -& 
&  6*(2*g1p2*MassB + 10*g2p2*MassWB + TradjYb3TYb3 + 10._dp*(TradjYuTYu) +               & 
&  5._dp*(TradjYw3TYw3) + 10._dp*(TradjYx3TYx3))*Yb3adjTYb3 - 6._dp*(Yb3adjTYb3TYb3adjYb3) -& 
&  20._dp*(Yb3adjTYeTYeadjYb3) - 30._dp*(Yb3adjTYw3TYw3adjYb3) + 6*(4*AbsMassB*g1p2 +    & 
&  20*AbsMassWB*g2p2 + 2*g1p2*mHu2 + 10*g2p2*mHu2 - TrCTYb3TpTYb3 - 10._dp*(TrCTYuTpTYu) -& 
&  5._dp*(TrCTYw3TpTYw3) - 10._dp*(TrCTYx3TpTYx3) - 10._dp*(Trmd2adjYx3Yx3) -            & 
&  TrmHb32Yb3adjYb3 - 5._dp*(TrmHw32Yw3adjYw3) - 10._dp*(TrmHxb32Yx3adjYx3) -            & 
&  Trml2adjYb3Yb3 - 5._dp*(Trml2adjYw3Yw3) - 10._dp*(Trmq2adjYuYu) - 10._dp*(Trmu2YuadjYu) -& 
&  2*mHu2*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) + TrYw3adjYw3 + 2._dp*(TrYx3adjYx3))))*Yb3adjYb3 -& 
&  3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) -  & 
&  5._dp*(TrYx3adjYx3)))*Yb3adjYb3mHb32 - 6._dp*(Yb3adjYb3mHb32Yb3adjYb3) -              & 
&  6._dp*(Yb3adjYb3TYb3adjTYb3) - 3._dp*(Yb3adjYb3Yb3adjYb3mHb32) - 6._dp*(Yb3adjYb3Yb3ml2adjYb3) -& 
&  20._dp*(Yb3adjYeme2YeadjYb3) - 20._dp*(Yb3adjYeTYeadjTYb3) - 20*mHd2*Yb3adjYeYeadjYb3 -& 
&  10._dp*(Yb3adjYeYeadjYb3mHb32) - 20._dp*(Yb3adjYeYeml2adjYb3) - 30._dp*(Yb3adjYw3mHw32Yw3adjYb3) -& 
&  30._dp*(Yb3adjYw3TYw3adjTYb3) - 4*mHu2*(3._dp*(Yb3adjYb3Yb3adjYb3) + 5*(Yb3adjYeYeadjYb3 +& 
&  3._dp*(Yb3adjYw3Yw3adjYb3))) - 15._dp*(Yb3adjYw3Yw3adjYb3mHb32) - 30._dp*(Yb3adjYw3Yw3ml2adjYb3) -& 
&  6*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) -  & 
&  5._dp*(TrYx3adjYx3)))*Yb3ml2adjYb3 - 6._dp*(Yb3ml2adjYb3Yb3adjYb3) - 20._dp*(Yb3ml2adjYeYeadjYb3) -& 
&  30._dp*(Yb3ml2adjYw3Yw3adjYb3) - 6*TYb3adjYb3*(TrCTYb3TpYb3 + 5*(2._dp*(TrCTYuTpYu) + & 
&  TrCTYw3TpYw3 + 2._dp*(TrCTYx3TpYx3)) + 2*g1p2*Conjg(MassB) + 10*g2p2*Conjg(MassWB))))/50._dp

 
DmHb32 = oo16pi2*( betamHb321 + oo16pi2 * betamHb322 ) 

 
Else 
DmHb32 = oo16pi2* betamHb321 
End If 
 
 
Forall(i1=1:3) DmHb32(i1,i1) =  Real(DmHb32(i1,i1),dp) 
!-------------------- 
! mHx32 
!-------------------- 
 
betamHx321  = -(id3R*(18*AbsMassWB*g2p2 + 32*AbsMassG*g3p2 + 5*g1p2*(2._dp*(AbsMassB) & 
&  - Tr1(1))))/3._dp

 
 
If (TwoLoopRGE) Then 
betamHx322 = (id3R*(g1p2*Conjg(MassB)*(669*g1p2*MassB + 5*(16*g3p2*(2._dp*(MassB) + MassG) +       & 
&  9*g2p2*(2._dp*(MassB) + MassWB)) + 225*g1p2*MassB*(NGHx3 +             & 
&  NGHxb3)) + 16*g3p2*Conjg(MassG)*(5*g1p2*(MassB + 2._dp*(MassG)) +      & 
&  3*(-8*g3p2*MassG + 3*g2p2*(2._dp*(MassG) + MassWB)) + 18*g3p2*MassG*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)) + 3*(3*g2p2*Conjg(MassWB)*(33*g2p2*MassWB +& 
&  5*g1p2*(MassB + 2._dp*(MassWB)) + 16*g3p2*(MassG + 2._dp*(MassWB)) + 9*g2p2*MassWB*(4*NGHw3 +& 
&  3*(NGHx3 + NGHxb3))) + 2*(5*g1p4*Tr2(1) +               & 
&  9*g2p4*Tr2(2) + 16*g3p4*Tr2(3) + 10*g1p2*Tr3(1)))))/9._dp

 
DmHx32 = oo16pi2*( betamHx321 + oo16pi2 * betamHx322 ) 

 
Else 
DmHx32 = oo16pi2* betamHx321 
End If 
 
 
Forall(i1=1:3) DmHx32(i1,i1) =  Real(DmHx32(i1,i1),dp) 
!-------------------- 
! mHxb32 
!-------------------- 
 
betamHxb321  = mHxb32Yx3adjYx3 + 2._dp*(TYx3adjTYx3) + 2*mHu2*Yx3adjYx3 +             & 
&  Yx3adjYx3mHxb32 + 2._dp*(Yx3md2adjYx3) - (id3R*(18*AbsMassWB*g2p2 + 32*AbsMassG*g3p2 +& 
&  5*g1p2*(2._dp*(AbsMassB) + Tr1(1))))/3._dp

 
 
If (TwoLoopRGE) Then 
betamHxb322 = 10*AbsMassWB*g1p2*g2p2*id3R + 33*AbsMassWB*g2p4*id3R + 32*AbsMassWB*g2p2*g3p2*id3R -  & 
&  (2*g1p2*mHxb32Yx3adjYx3)/5._dp - (3*mHxb32Yx3adjYx3*TrYb3adjYb3)/10._dp -             & 
&  3*mHxb32Yx3adjYx3*TrYuadjYu - (3*mHxb32Yx3adjYx3*TrYw3adjYw3)/2._dp - 3*mHxb32Yx3adjYx3*TrYx3adjYx3 -& 
&  (4*g1p2*TYx3adjTYx3)/5._dp - (3*TrYb3adjYb3*TYx3adjTYx3)/5._dp - 6*TrYuadjYu*TYx3adjTYx3 -& 
&  3*TrYw3adjYw3*TYx3adjTYx3 - 6*TrYx3adjYx3*TYx3adjTYx3 - (3*(TrCTYb3TpYb3 +            & 
&  5*(2._dp*(TrCTYuTpYu) + TrCTYw3TpYw3 + 2._dp*(TrCTYx3TpYx3)))*TYx3adjYx3)/5._dp +     & 
&  (4*g1p2*MassB*Yx3adjTYx3)/5._dp - (3*TradjYb3TYb3*Yx3adjTYx3)/5._dp - 6*TradjYuTYu*Yx3adjTYx3 -& 
&  3*TradjYw3TYw3*Yx3adjTYx3 - 6*TradjYx3TYx3*Yx3adjTYx3 - (4*g1p2*mHu2*Yx3adjYx3)/5._dp -& 
&  (3*TrCTYb3TpTYb3*Yx3adjYx3)/5._dp - 6*TrCTYuTpTYu*Yx3adjYx3 - 3*TrCTYw3TpTYw3*Yx3adjYx3 -& 
&  6*TrCTYx3TpTYx3*Yx3adjYx3 - 6*Trmd2adjYx3Yx3*Yx3adjYx3 - (3*TrmHb32Yb3adjYb3*Yx3adjYx3)/5._dp -& 
&  3*TrmHw32Yw3adjYw3*Yx3adjYx3 - 6*TrmHxb32Yx3adjYx3*Yx3adjYx3 - (3*Trml2adjYb3Yb3*Yx3adjYx3)/5._dp -& 
&  3*Trml2adjYw3Yw3*Yx3adjYx3 - 6*Trmq2adjYuYu*Yx3adjYx3 - 6*Trmu2YuadjYu*Yx3adjYx3 -    & 
&  (6*mHu2*TrYb3adjYb3*Yx3adjYx3)/5._dp - 12*mHu2*TrYuadjYu*Yx3adjYx3 - 6*mHu2*TrYw3adjYw3*Yx3adjYx3 -& 
&  12*mHu2*TrYx3adjYx3*Yx3adjYx3 - (2*g1p2*Yx3adjYx3mHxb32)/5._dp - (3*TrYb3adjYb3*Yx3adjYx3mHxb32)/10._dp -& 
&  3*TrYuadjYu*Yx3adjYx3mHxb32 - (3*TrYw3adjYw3*Yx3adjYx3mHxb32)/2._dp - 3*TrYx3adjYx3*Yx3adjYx3mHxb32 -& 
&  8*mHu2*Yx3adjYx3Yx3adjYx3 - 4*mHu2*Yx3CYdTpYdadjYx3 - (4*g1p2*Yx3md2adjYx3)/5._dp -   & 
&  (3*TrYb3adjYb3*Yx3md2adjYx3)/5._dp - 6*TrYuadjYu*Yx3md2adjYx3 - 3*TrYw3adjYw3*Yx3md2adjYx3 -& 
&  6*TrYx3adjYx3*Yx3md2adjYx3 - 2*(mHxb32Yx3adjYx3Yx3adjYx3 + mHxb32Yx3CYdTpYdadjYx3 +   & 
&  2._dp*(TYx3adjYx3Yx3adjTYx3) + 2._dp*(TYx3CYdTpYdadjTYx3) + 2._dp*(TYx3TpYx3CTYx3adjYx3) +& 
&  2._dp*(TYx3YdadjTYdadjYx3) + 2._dp*(Yx3adjTYx3TYx3adjYx3) + 2._dp*(Yx3adjYx3mHxb32Yx3adjYx3) +& 
&  2._dp*(Yx3adjYx3TYx3adjTYx3) + Yx3adjYx3Yx3adjYx3mHxb32 + 2._dp*(Yx3adjYx3Yx3md2adjYx3) +& 
&  2._dp*(Yx3CTYdTpTYdadjYx3) + 2*mHd2*Yx3CYdTpYdadjYx3 + Yx3CYdTpYdadjYx3mHxb32 +       & 
&  2*(Yx3CYdmq2TpYdadjYx3 + Yx3CYdTpYdmd2adjYx3) + 2._dp*(Yx3md2adjYx3Yx3adjYx3) +       & 
&  2._dp*(Yx3md2CYdTpYdadjYx3) + 2._dp*(Yx3TYdadjYdadjTYx3)) + 5*g1p2*g2p2*id3R*MassB*Conjg(MassWB)  
betamHxb322 =  betamHxb322+ 16*g2p2*g3p2*id3R*MassG*Conjg(MassWB) + 36*AbsMassWB*g2p4*id3R*NGHw3 + & 
&  27*AbsMassWB*g2p4*id3R*NGHx3 + 27*AbsMassWB*g2p4*id3R*NGHxb3 +& 
&  (16*g3p2*id3R*Conjg(MassG)*(5*g1p2*(MassB + 2._dp*(MassG)) + 3*(-8*g3p2*MassG +       & 
&  3*g2p2*(2._dp*(MassG) + MassWB)) + 18*g3p2*MassG*(3*NGHg3 +            & 
&  NGHx3 + NGHxb3)))/9._dp + (g1p2*Conjg(MassB)*(36*(TYx3adjYx3 -& 
&  2*MassB*Yx3adjYx3) + 5*id3R*(669*g1p2*MassB + 5*(16*g3p2*(2._dp*(MassB) +             & 
&  MassG) + 9*g2p2*(2._dp*(MassB) + MassWB)) + 225*g1p2*MassB*(NGHx3 +    & 
&  NGHxb3))))/45._dp + (10*g1p4*id3R*Tr2(1))/3._dp + 6*g2p4*id3R*Tr2(2) + & 
&  (32*g3p4*id3R*Tr2(3))/3._dp - (20*g1p2*id3R*Tr3(1))/3._dp

 
DmHxb32 = oo16pi2*( betamHxb321 + oo16pi2 * betamHxb322 ) 

 
Else 
DmHxb32 = oo16pi2* betamHxb321 
End If 
 
 
Forall(i1=1:3) DmHxb32(i1,i1) =  Real(DmHxb32(i1,i1),dp) 
!-------------------- 
! MassB 
!-------------------- 
 
betaMassB1  = (g1p2*MassB*(66 + 25*NGHx3 + 25*NGHxb3))/5._dp

 
 
If (TwoLoopRGE) Then 
betaMassB2 = (g1p2*(6*(398*g1p2*MassB + 5*(88*g3p2*(MassB + MassG) + 27*g2p2*(MassB +              & 
&  MassWB)) + 9._dp*(TradjYb3TYb3) + 70._dp*(TradjYdTYd) + 90._dp*(TradjYeTYe) +         & 
&  130._dp*(TradjYuTYu) + 60._dp*(TradjYw3TYw3) + 190._dp*(TradjYx3TYx3) -               & 
&  9*MassB*TrYb3adjYb3 - 10*MassB*(7._dp*(TrYdadjYd) + 9._dp*(TrYeadjYe) +               & 
&  13._dp*(TrYuadjYu) + 6._dp*(TrYw3adjYw3) + 19._dp*(TrYx3adjYx3))) + 125*(10*g1p2*MassB +& 
&  16*g3p2*(MassB + MassG) + 9*g2p2*(MassB + MassWB))*NGHx3 +             & 
&  125*(10*g1p2*MassB + 16*g3p2*(MassB + MassG) + 9*g2p2*(MassB + MassWB))*NGHxb3))/75._dp

 
DMassB = oo16pi2*( betaMassB1 + oo16pi2 * betaMassB2 ) 

 
Else 
DMassB = oo16pi2* betaMassB1 
End If 
 
 
!-------------------- 
! MassWB 
!-------------------- 
 
betaMassWB1  = g2p2*MassWB*(2 + 4*NGHw3 + 3*NGHx3 +     & 
&  3*NGHxb3)

 
 
If (TwoLoopRGE) Then 
betaMassWB2 = (g2p2*(2*(9._dp*(TradjYb3TYb3) + 90._dp*(TradjYdTYd) + 30._dp*(TradjYeTYe) +          & 
&  90._dp*(TradjYuTYu) + 140._dp*(TradjYw3TYw3) + 3*(9*g1p2*(MassB + MassWB) +           & 
&  10*(25*g2p2*MassWB + 12*g3p2*(MassG + MassWB)) + 30._dp*(TradjYx3TYx3))) -            & 
&  18*MassWB*TrYb3adjYb3 - 180*MassWB*TrYdadjYd - 60*MassWB*TrYeadjYe - 180*MassWB*TrYuadjYu -& 
&  280*MassWB*TrYw3adjYw3 - 180*MassWB*TrYx3adjYx3 + 1440*g2p2*MassWB*NGHw3 +& 
&  15*(5*g1p2*(MassB + MassWB) + 2*(21*g2p2*MassWB + 8*g3p2*(MassG + MassWB)))*NGHx3 +& 
&  75*g1p2*MassB*NGHxb3 + 240*g3p2*MassG*NGHxb3 +          & 
&  75*g1p2*MassWB*NGHxb3 + 630*g2p2*MassWB*NGHxb3 +        & 
&  240*g3p2*MassWB*NGHxb3))/15._dp

 
DMassWB = oo16pi2*( betaMassWB1 + oo16pi2 * betaMassWB2 ) 

 
Else 
DMassWB = oo16pi2* betaMassWB1 
End If 
 
 
!-------------------- 
! MassG 
!-------------------- 
 
betaMassG1  = 2*g3p2*MassG*(-3 + 3*NGHg3 + NGHx3 + NGHxb3)

 
 
If (TwoLoopRGE) Then 
betaMassG2 = (2*g3p2*(3*(11*g1p2*(MassB + MassG) + 5*(28*g3p2*MassG + 9*g2p2*(MassG +              & 
&  MassWB)) + 20._dp*(TradjYdTYd) + 20._dp*(TradjYuTYu) + 20._dp*(TradjYx3TYx3)) +       & 
&  1620*g3p2*MassG*NGHg3 + 5*(68*g3p2*MassG + 5*g1p2*(MassB +             & 
&  MassG) + 9*g2p2*(MassG + MassWB))*NGHx3 + 5*(-12*MassG*(TrYdadjYd +    & 
&  TrYuadjYu + TrYx3adjYx3) + (68*g3p2*MassG + 5*g1p2*(MassB + MassG) + 9*g2p2*(MassG +  & 
&  MassWB))*NGHxb3)))/15._dp

 
DMassG = oo16pi2*( betaMassG1 + oo16pi2 * betaMassG2 ) 

 
Else 
DMassG = oo16pi2* betaMassG1 
End If 

!-------------------- 
! MnuL 
!-------------------- 
 

Do i1 = 1,3
Do i2 = 1,3
gammaL1(i1,i2) =  3._dp*adjYb3Yb3(i1,i2)/10._dp + adjYeYe(i1,i2) + 3._dp*adjYw3Yw3(i1,i2)/2._dp

If(i1.eq.i2) Then
 gammaL1(i1,i2) = gammaL1(i1,i2) -0.3_dp*g1p2 - 1.5_dp*g2p2
End If

End Do
End Do

gammaHd1 = -0.3_dp*g1p2 - 1.5_dp*g2p2 + 0.3_dp*TrYb3adjYb3 + &
& 3._dp*TrYuadjYu + 1.5_dp*TrYw3adjYw3 + 3._dp*TrYx3adjYx3

betaMnuL1 = MatMul(Transpose(gammaL1),MnuL) + MatMul(MnuL,gammaL1) + 2._dp*MnuL*gammaHd1

 
If (TwoLoopRGE) Then 

Do i1 = 1,3
Do i2 = 1,3
gammaL2(i1,i2) = -9._dp*TrYb3adjYb3*adjYb3Yb3(i1,i2)/100._dp - &
& 9._dp*TrYuadjYu*adjYb3Yb3(i1,i2)/10._dp - &
& 9._dp*TrYw3adjYw3*adjYb3Yb3(i1,i2)/20._dp - &
& 9._dp*TrYx3adjYx3*adjYb3Yb3(i1,i2)/10._dp - &
& 9._dp*adjYb3Yb3adjYb3Yb3(i1,i2)/50._dp - &
& 3._dp*adjYb3Yb3adjYw3Yw3(i1,i2)/10._dp + 6._dp*g1p2*adjYeYe(i1,i2)/5._dp - &
& 3._dp*TrYdadjYd*adjYeYe(i1,i2) - TrYeadjYe*adjYeYe(i1,i2) - &
& 2._dp*adjYeYeadjYeYe(i1,i2) + 6._dp*g2p2*adjYw3Yw3(i1,i2) - &
& 9._dp*TrYb3adjYb3*adjYw3Yw3(i1,i2)/20._dp - &
& 9._dp*TrYuadjYu*adjYw3Yw3(i1,i2)/2._dp - &
& 9._dp*TrYw3adjYw3*adjYw3Yw3(i1,i2)/4._dp - &
& 9._dp*TrYx3adjYx3*adjYw3Yw3(i1,i2)/2._dp - &
& 9._dp*adjYw3Yw3adjYb3Yb3(i1,i2)/40._dp - 3._dp*adjYw3Yw3adjYw3Yw3(i1,i2)/2._dp

If(i1.eq.i2) Then
 gammaL2(i1,i2) =  gammaL2(i1,i2) + 657._dp*g1**4/100._dp + 9._dp*g1p2*g2p2/10._dp &
	& + 105._dp*g2**4/4._dp
End If

End Do
End Do

gammaHd2 = 657._dp*g1**4/100._dp + 9._dp*g1p2*g2p2/10._dp + 105._dp*g2**4/4._dp - &
& 27._dp*TrYb3adjYb3Yb3adjYb3/100._dp - 3._dp*TrYb3adjYeYeadjYb3/10._dp - &
& 3._dp*TrYb3adjYw3Yw3adjYb3 /4._dp - 6._dp*TrYdadjYdTpYx3CYx3 - &
& 3._dp*TrYdadjYuYuadjYd + 4._dp*g1p2*TrYuadjYu/5._dp + 16._dp*g3p2*TrYuadjYu &
& - 9._dp*TrYuadjYuYuadjYu - 27._dp*TrYb3adjYw3Yw3adjYb3 /40._dp - &
& 3._dp*TrYeadjYw3Yw3adjYe/2._dp + 6._dp*g2p2*TrYw3adjYw3 - &
& 15._dp*TrYw3adjYw3Yw3adjYw3/4._dp + 2._dp*g1p2*TrYx3adjYx3 + &
& 16._dp*g3p2*TrYx3adjYx3 - 9._dp*TrYx3adjYx3Yx3adjYx3

betaMnuL2 = MatMul(Transpose(gammaL2),MnuL) + MatMul(MnuL,gammaL2) + 2._dp*MnuL*gammaHd2


 
DMnuL = oo16pi2*( betaMnuL1 + oo16pi2 * betaMnuL2 ) 

 
Else 
DMnuL = oo16pi2* betaMnuL1
End If 
 
 
If (ThresholdCrossed.lt.1) Then 
DMXM3(1,:) = 0._dp 
DBMXM3(1,:) = 0._dp 
DmHx32(1,:) = 0._dp 
DmHx32(:,1) = 0._dp 
DYx3(1,:) = 0._dp 
DTYx3(1,:) = 0._dp 
DMXM3(:,1) = 0._dp 
DBMXM3(:,1) = 0._dp 
DmHxb32(1,:) = 0._dp 
DmHxb32(:,1) = 0._dp 
DMGM3(1,:) = 0._dp 
DBMGM3(1,:) = 0._dp 
DMGM3(:,1) = 0._dp 
DBMGM3(:,1) = 0._dp 
DmHg32(1,:) = 0._dp 
DmHg32(:,1) = 0._dp 
DYb3(1,:) = 0._dp 
DTYb3(1,:) = 0._dp 
DMBM3(1,:) = 0._dp 
DBMBM3(1,:) = 0._dp 
DMBM3(:,1) = 0._dp 
DBMBM3(:,1) = 0._dp 
DmHb32(1,:) = 0._dp 
DmHb32(:,1) = 0._dp 
DYw3(1,:) = 0._dp 
DTYw3(1,:) = 0._dp 
DMWM3(1,:) = 0._dp 
DBMWM3(1,:) = 0._dp 
DMWM3(:,1) = 0._dp 
DBMWM3(:,1) = 0._dp 
DmHw32(1,:) = 0._dp 
DmHw32(:,1) = 0._dp 
End if 

If (ThresholdCrossed.lt.2) Then 
DMXM3(2,:) = 0._dp 
DBMXM3(2,:) = 0._dp 
DmHx32(2,:) = 0._dp 
DmHx32(:,2) = 0._dp 
DYx3(2,:) = 0._dp 
DTYx3(2,:) = 0._dp 
DMXM3(:,2) = 0._dp 
DBMXM3(:,2) = 0._dp 
DmHxb32(2,:) = 0._dp 
DmHxb32(:,2) = 0._dp 
DMGM3(2,:) = 0._dp 
DBMGM3(2,:) = 0._dp 
DMGM3(:,2) = 0._dp 
DBMGM3(:,2) = 0._dp 
DmHg32(2,:) = 0._dp 
DmHg32(:,2) = 0._dp 
DYb3(2,:) = 0._dp 
DTYb3(2,:) = 0._dp 
DMBM3(2,:) = 0._dp 
DBMBM3(2,:) = 0._dp 
DMBM3(:,2) = 0._dp 
DBMBM3(:,2) = 0._dp 
DmHb32(2,:) = 0._dp 
DmHb32(:,2) = 0._dp 
DYw3(2,:) = 0._dp 
DTYw3(2,:) = 0._dp 
DMWM3(2,:) = 0._dp 
DBMWM3(2,:) = 0._dp 
DMWM3(:,2) = 0._dp 
DBMWM3(:,2) = 0._dp 
DmHw32(2,:) = 0._dp 
DmHw32(:,2) = 0._dp 
End if 

If (ThresholdCrossed.lt.3) Then 
DMXM3(3,:) = 0._dp 
DBMXM3(3,:) = 0._dp 
DmHx32(3,:) = 0._dp 
DmHx32(:,3) = 0._dp 
DYx3(3,:) = 0._dp 
DTYx3(3,:) = 0._dp 
DMXM3(:,3) = 0._dp 
DBMXM3(:,3) = 0._dp 
DmHxb32(3,:) = 0._dp 
DmHxb32(:,3) = 0._dp 
DMGM3(3,:) = 0._dp 
DBMGM3(3,:) = 0._dp 
DMGM3(:,3) = 0._dp 
DBMGM3(:,3) = 0._dp 
DmHg32(3,:) = 0._dp 
DmHg32(:,3) = 0._dp 
DYb3(3,:) = 0._dp 
DTYb3(3,:) = 0._dp 
DMBM3(3,:) = 0._dp 
DBMBM3(3,:) = 0._dp 
DMBM3(:,3) = 0._dp 
DBMBM3(:,3) = 0._dp 
DmHb32(3,:) = 0._dp 
DmHb32(:,3) = 0._dp 
DYw3(3,:) = 0._dp 
DTYw3(3,:) = 0._dp 
DMWM3(3,:) = 0._dp 
DBMWM3(3,:) = 0._dp 
DMWM3(:,3) = 0._dp 
DBMWM3(:,3) = 0._dp 
DmHw32(3,:) = 0._dp 
DmHw32(:,3) = 0._dp 
End if 

!-------------------------------------------------------------------------------
! these matrices are hermitian but numerical effects induce non-hermiatian parts
! which need to be corrected
!-------------------------------------------------------------------------------
Dmd2 = 0.5_dp * ( Dmd2 + Transpose(Conjg(Dmd2)) )
Dme2 = 0.5_dp * ( Dme2 + Transpose(Conjg(Dme2)) )
Dml2 = 0.5_dp * ( Dml2 + Transpose(Conjg(Dml2)) )
Dmq2 = 0.5_dp * ( Dmq2 + Transpose(Conjg(Dmq2)) )
Dmu2 = 0.5_dp * ( Dmu2 + Transpose(Conjg(Dmu2)) )

DmHb32 = 0.5_dp * ( DmHb32 + Transpose(Conjg(DmHb32)) )
DmHg32 = 0.5_dp * ( DmHg32 + Transpose(Conjg(DmHg32)) )
DmHw32 = 0.5_dp * ( DmHw32 + Transpose(Conjg(DmHw32)) )
DmHx32 = 0.5_dp * ( DmHx32 + Transpose(Conjg(DmHx32)) )

Call Chop(Dmue)
Call Chop(DBmue)

Call ParametersToG555(Dg1,Dg2,Dg3,DYu,DYd,DYe,DYb3,DYw3,DYx3,Dmue,DMXM3,       & 
& DMWM3,DMGM3,DMBM3,DTYu,DTYd,DTYe,DTYb3,DTYw3,DTYx3,DBmue,DBMXM3,DBMWM3,      & 
& DBMGM3,DBMBM3,Dmq2,Dml2,DmHd2,DmHu2,Dmd2,Dmu2,Dme2,DmHw32,DmHg32,DmHb32,     & 
& DmHx32, DmHxb32,DMassB,DMassWB,DMassG,DMnuL,f)

Iname = Iname - 1 
 
End Subroutine rge555

#endif SEESAWIII


#ifdef SEESAWIII

! ----------------------------------------------------------------------  
! This file was automatically created by SARAH version SARAHVERSION 
! SARAH References: arXiv:0806.0538, arXiv:0909.2863, arXiv:1002.0840    
! (c) Florian Staub, 2011  
! ----------------------------------------------------------------------  
! File created at 13:17 on 29.2.2012   
! ----------------------------------------------------------------------  
 
 Subroutine GToParameters115(g,g1,g2,g3,Yu,Yd,Ye,YtII,YsII,YzII,L1II,L2II)

 Implicit None 
  Real(dp), Intent(in) :: g(115) 
  Real(dp),Intent(out) :: g1,g2,g3

  Complex(dp),Intent(out) :: Yu(3,3),Yd(3,3),Ye(3,3),YtII(3,3),YsII(3,3),YzII(3,3) &
     & ,L1II,L2II

  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'GToParameters115' 
 
  g1 = g(1) 
  g2 = g(2) 
  g3 = g(3)

  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ((i2-1) + (i1-1)*3)
    Yu(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
    Yd(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp) 
    Ye(i1,i2) = Cmplx( g(SumI+40), g(SumI+41), dp) 
    YtII(i1,i2) = Cmplx( g(SumI+58), g(SumI+59), dp) 
    YsII(i1,i2) = Cmplx( g(SumI+76), g(SumI+77), dp) 
    YzII(i1,i2) = Cmplx( g(SumI+94), g(SumI+95), dp) 
   End Do 
  End Do 
 
  L1II= Cmplx(g(112),g(113),dp) 
  L2II= Cmplx(g(114),g(115),dp) 
  Do i1=1,115 
   If (g(i1).ne.g(i1)) Then 
    Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
    Write(*,*) "At position ", i1 
    Call TerminateProgram 
   End if 
  End do 
  Iname = Iname - 1 
 
 End Subroutine GToParameters115


 Subroutine ParametersToG115(g1,g2,g3,Yu,Yd,Ye,YtII,YsII,YzII,L1II,L2II,g)

 Implicit None 
  Real(dp), Intent(out) :: g(115) 
  Real(dp), Intent(in) :: g1,g2,g3

  Complex(dp), Intent(in) :: Yu(3,3),Yd(3,3),Ye(3,3),YtII(3,3),YsII(3,3),YzII(3,3) &
    & ,L1II,L2II

  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'ParametersToG115' 
 
  g(1) = g1
  g(2) = g2
  g(3) = g3
  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3 )
    g(SumI+4) = Real(Yu(i1,i2), dp) 
    g(SumI+5) = Aimag(Yu(i1,i2)) 
    g(SumI+22) = Real(Yd(i1,i2), dp) 
    g(SumI+23) = Aimag(Yd(i1,i2)) 
    g(SumI+40) = Real(Ye(i1,i2), dp) 
    g(SumI+41) = Aimag(Ye(i1,i2)) 
    g(SumI+58) = Real(YtII(i1,i2), dp) 
    g(SumI+59) = Aimag(YtII(i1,i2)) 
    g(SumI+76) = Real(YsII(i1,i2), dp) 
    g(SumI+77) = Aimag(YsII(i1,i2)) 
    g(SumI+94) = Real(YzII(i1,i2), dp) 
    g(SumI+95) = Aimag(YzII(i1,i2)) 
   End Do 
  End Do 

  g(112) = Real(L1II,dp)  
  g(113) = Aimag(L1II)  
  g(114) = Real(L2II,dp)  
  g(115) = Aimag(L2II)  
  Iname = Iname - 1 
 
 End Subroutine ParametersToG115

Subroutine rge115(len, T, GY, F) 
Implicit None 
Integer, Intent(in) :: len 
Real(dp), Intent(in) :: T, GY(len) 
Real(dp), Intent(out) :: F(len) 
Integer :: i2
Real(dp) :: q 
Real(dp) :: g1,betag11,betag12,Dg1,g2,betag21,betag22,Dg2,g3,betag31,betag32,Dg3
Complex(dp) :: Yu(3,3),betaYu1(3,3),betaYu2(3,3),DYu(3,3),adjYu(3,3),Yd(3,3)          & 
& ,betaYd1(3,3),betaYd2(3,3),DYd(3,3),adjYd(3,3),Ye(3,3),betaYe1(3,3),betaYe2(3,3)       & 
& ,DYe(3,3),adjYe(3,3),YtII(3,3),betaYtII1(3,3),betaYtII2(3,3),DYtII(3,3),               & 
& adjYtII(3,3),YsII(3,3),betaYsII1(3,3),betaYsII2(3,3),DYsII(3,3),adjYsII(3,3)           & 
& ,YzII(3,3),betaYzII1(3,3),betaYzII2(3,3),DYzII(3,3),adjYzII(3,3),L1II,betaL1II1,       & 
& betaL1II2,DL1II,L2II,betaL2II1,betaL2II2,DL2II
Real(dp) :: AbsL1II,AbsL2II
Complex(dp) :: YdadjYd(3,3),YeadjYe(3,3),YsIICYsII(3,3),YtIICYtII(3,3),YuadjYu(3,3),YzIIadjYzII(3,3),& 
& adjYdYd(3,3),adjYdYsII(3,3),adjYdYzII(3,3),adjYeYe(3,3),adjYuYu(3,3),adjYzIIYd(3,3),   & 
& adjYzIIYsII(3,3),adjYzIIYzII(3,3),CYdTpYd(3,3),CYeYtII(3,3),CYsIIYd(3,3),              & 
& CYsIIYsII(3,3),CYsIIYzII(3,3),CYtIIYtII(3,3),CYzIIYtII(3,3),CYzIITpYzII(3,3),          & 
& YdadjYdYd(3,3),YdadjYdYsII(3,3),YdadjYdYzII(3,3),YdadjYuYu(3,3),YeadjYeYe(3,3),        & 
& YeadjYzIIYzII(3,3),YeCYtIIYtII(3,3),YsIICYdTpYd(3,3),YsIICYsIIYd(3,3),YsIICYsIIYsII(3,3),& 
& YsIICYsIIYzII(3,3),YsIICYzIITpYzII(3,3),YtIIadjYeYe(3,3),YtIIadjYzIIYzII(3,3),         & 
& YtIICYtIIYtII(3,3),YuadjYdYd(3,3),YuadjYuYu(3,3),YzIIadjYeYe(3,3),YzIIadjYzIIYd(3,3),  & 
& YzIIadjYzIIYsII(3,3),YzIIadjYzIIYzII(3,3),YzIICYtIIYtII(3,3),TpYeCYeYtII(3,3),         & 
& TpYzIICYzIIYtII(3,3)

Complex(dp) :: YtIIadjYe(3,3),YuadjYd(3,3),YzIIadjYe(3,3),YzIICYtII(3,3),CYeTpYzII(3,3),             & 
& CYtIITpYzII(3,3),CYuTpYd(3,3),YeadjYzIIYd(3,3),YeadjYzIIYsII(3,3),YsIICYzIIYtII(3,3),  & 
& YtIIadjYzIIYd(3,3),YtIIadjYzIIYsII(3,3),YtIICYtIITpYzII(3,3),YuadjYdYsII(3,3),         & 
& YuadjYdYzII(3,3),adjYdYdadjYd(3,3),adjYdYsIICYsII(3,3),adjYdYzIIadjYzII(3,3),          & 
& adjYeYeadjYe(3,3),adjYuYuadjYd(3,3),adjYuYuadjYu(3,3),adjYzIIYzIIadjYe(3,3),           & 
& adjYzIIYzIIadjYzII(3,3),adjYzIIYzIICYtII(3,3),CYsIIYsIICYsII(3,3),CYsIIYzIIadjYzII(3,3),& 
& CYtIIYtIIadjYe(3,3),CYtIIYtIICYtII(3,3),TpYdCYdTpYd(3,3),TpYdCYsIIYd(3,3),             & 
& TpYdCYsIIYsII(3,3),TpYdCYsIIYzII(3,3),TpYdCYzIIYtII(3,3),TpYeCYeTpYzII(3,3),           & 
& TpYuCYuTpYd(3,3),TpYzIICYsIIYd(3,3),TpYzIICYsIIYsII(3,3),TpYzIICYsIIYzII(3,3),         & 
& TpYzIICYzIITpYzII(3,3),YdadjYdYdadjYd(3,3),YdadjYdYsIICYsII(3,3),YdadjYdYzIIadjYzII(3,3),& 
& YdadjYuYuadjYd(3,3),YeadjYeYeadjYe(3,3),YeadjYzIIYzIIadjYe(3,3),YeCYtIIYtIIadjYe(3,3), & 
& YsIICYsIIYsIICYsII(3,3),YsIICYsIIYzIIadjYzII(3,3),YtIIadjYzIIYzIICYtII(3,3),           & 
& YtIICYtIIYtIICYtII(3,3),YuadjYuYuadjYu(3,3),YzIIadjYzIIYzIIadjYzII(3,3),               & 
& adjYdYdadjYdYd(3,3),adjYdYdadjYdYsII(3,3),adjYdYdadjYdYzII(3,3),adjYdYdadjYuYu(3,3),   & 
& adjYdYsIICYsIIYd(3,3),adjYdYzIIadjYzIIYd(3,3),adjYeYeadjYeYe(3,3),adjYeYeadjYzIIYd(3,3),& 
& adjYeYeadjYzIIYsII(3,3),adjYeYeadjYzIIYzII(3,3),adjYeYeCYtIIYtII(3,3),adjYuYuadjYdYd(3,3),& 
& adjYuYuadjYdYsII(3,3),adjYuYuadjYdYzII(3,3),adjYuYuadjYuYu(3,3),adjYzIIYdadjYdYzII(3,3),& 
& adjYzIIYsIICYsIIYzII(3,3),adjYzIIYzIIadjYeYe(3,3),adjYzIIYzIIadjYzIIYd(3,3),           & 
& adjYzIIYzIIadjYzIIYsII(3,3),adjYzIIYzIIadjYzIIYzII(3,3),adjYzIIYzIICYtIIYtII(3,3),     & 
& CYdTpYdCYdTpYd(3,3),CYdTpYdCYsIIYd(3,3),CYdTpYdCYsIIYsII(3,3),CYdTpYdCYsIIYzII(3,3),   & 
& CYdTpYdCYzIIYtII(3,3),CYdTpYuCYuTpYd(3,3),CYeTpYeCYeYtII(3,3),CYsIIYdadjYdYsII(3,3),   & 
& CYsIIYsIICYsIIYd(3,3),CYsIIYsIICYsIIYsII(3,3),CYsIIYsIICYsIIYzII(3,3),CYsIIYsIICYzIIYtII(3,3),& 
& CYsIIYzIIadjYzIIYsII(3,3),CYtIIYtIIadjYeYe(3,3),CYtIIYtIIadjYzIIYd(3,3),               & 
& CYtIIYtIIadjYzIIYsII(3,3),CYtIIYtIIadjYzIIYzII(3,3),CYtIIYtIICYtIIYtII(3,3),           & 
& CYtIITpYeCYeYtII(3,3),CYtIITpYzIICYzIIYtII(3,3),CYzIIYtIICYtIITpYzII(3,3),             & 
& CYzIITpYeCYeTpYzII(3,3),CYzIITpYzIICYsIIYd(3,3),CYzIITpYzIICYsIIYsII(3,3),             & 
& CYzIITpYzIICYsIIYzII(3,3),CYzIITpYzIICYzIIYtII(3,3),CYzIITpYzIICYzIITpYzII(3,3),       & 
& YdadjYdYdadjYdYd(3,3),YdadjYdYdadjYdYsII(3,3),YdadjYdYdadjYdYzII(3,3),YdadjYdYsIICYsIIYd(3,3),& 
& YdadjYdYzIIadjYzIIYd(3,3),YdadjYuYuadjYdYd(3,3),YdadjYuYuadjYdYsII(3,3),               & 
& YdadjYuYuadjYdYzII(3,3),YdadjYuYuadjYuYu(3,3),YeadjYeYeadjYeYe(3,3),YeadjYzIIYdadjYdYzII(3,3),& 
& YeadjYzIIYsIICYsIIYzII(3,3),YeadjYzIIYzIIadjYeYe(3,3),YeadjYzIIYzIIadjYzIIYzII(3,3),   & 
& YeCYtIIYtIIadjYeYe(3,3),YeCYtIIYtIICYtIIYtII(3,3),YeCYtIITpYeCYeYtII(3,3),             & 
& YeCYtIITpYzIICYzIIYtII(3,3),YsIICYdTpYdCYdTpYd(3,3),YsIICYdTpYdCYsIIYd(3,3),           & 
& YsIICYdTpYdCYsIIYsII(3,3),YsIICYdTpYdCYsIIYzII(3,3),YsIICYdTpYuCYuTpYd(3,3),           & 
& YsIICYsIIYdadjYdYsII(3,3),YsIICYsIIYsIICYsIIYd(3,3),YsIICYsIIYsIICYsIIYsII(3,3),       & 
& YsIICYsIIYsIICYsIIYzII(3,3),YsIICYsIIYzIIadjYzIIYsII(3,3),YsIICYzIIYtIICYtIITpYzII(3,3),& 
& YsIICYzIITpYeCYeTpYzII(3,3),YsIICYzIITpYzIICYsIIYd(3,3),YsIICYzIITpYzIICYsIIYsII(3,3)

Complex(dp) :: YsIICYzIITpYzIICYsIIYzII(3,3),YsIICYzIITpYzIICYzIITpYzII(3,3),YtIIadjYeYeadjYeYe(3,3), & 
& YtIIadjYeYeCYtIIYtII(3,3),YtIIadjYzIIYdadjYdYzII(3,3),YtIIadjYzIIYsIICYsIIYzII(3,3),   & 
& YtIIadjYzIIYzIIadjYzIIYzII(3,3),YtIIadjYzIIYzIICYtIIYtII(3,3),YtIICYtIIYtIICYtIIYtII(3,3),& 
& YtIICYtIITpYeCYeYtII(3,3),YtIICYtIITpYzIICYzIIYtII(3,3),YuadjYdYdadjYdYd(3,3),         & 
& YuadjYdYdadjYuYu(3,3),YuadjYdYsIICYsIIYd(3,3),YuadjYdYzIIadjYzIIYd(3,3),               & 
& YuadjYuYuadjYuYu(3,3),YzIIadjYeYeadjYeYe(3,3),YzIIadjYeYeadjYzIIYd(3,3),               & 
& YzIIadjYeYeadjYzIIYsII(3,3),YzIIadjYeYeadjYzIIYzII(3,3),YzIIadjYzIIYdadjYdYzII(3,3),   & 
& YzIIadjYzIIYsIICYsIIYzII(3,3),YzIIadjYzIIYzIIadjYzIIYd(3,3),YzIIadjYzIIYzIIadjYzIIYsII(3,3),& 
& YzIIadjYzIIYzIIadjYzIIYzII(3,3),YzIICYtIIYtIIadjYzIIYd(3,3),YzIICYtIIYtIIadjYzIIYsII(3,3),& 
& YzIICYtIIYtIIadjYzIIYzII(3,3),YzIICYtIIYtIICYtIIYtII(3,3),YzIICYtIITpYeCYeYtII(3,3),   & 
& YzIICYtIITpYzIICYzIIYtII(3,3),TpYeCYeTpYeCYeYtII(3,3),TpYzIICYdTpYdCYzIIYtII(3,3),     & 
& TpYzIICYsIIYsIICYzIIYtII(3,3),TpYzIICYzIITpYzIICYzIIYtII(3,3)

Complex(dp) :: TrYdadjYd,TrYeadjYe,TrYsIICYsII,TrYtIICYtII,TrYuadjYu,TrYzIIadjYzII

Complex(dp) :: TrYdadjYdYdadjYd,TrYdadjYdYsIICYsII,TrYdadjYdYzIIadjYzII,TrYdadjYuYuadjYd,            & 
& TrYeadjYeYeadjYe,TrYeadjYzIIYzIIadjYe,TrYeCYtIIYtIIadjYe,TrYsIICYsIIYsIICYsII,         & 
& TrYsIICYsIIYzIIadjYzII,TrYtIIadjYzIIYzIICYtII,TrYtIICYtIIYtIICYtII,TrYuadjYuYuadjYu,   & 
& TrYzIIadjYzIIYzIIadjYzII

Real(dp) :: g1p2,g1p3,g2p2,g2p3,g3p2,g3p3

Complex(dp) :: L1IIp2,L2IIp2

Real(dp) :: g1p4,g2p4,g3p4

Complex(dp) :: CL1IIp2,CL2IIp2

Iname = Iname +1 
NameOfUnit(Iname) = 'rge115' 
 
OnlyDiagonal = .Not.GenerationMixing 
q = t 
 
Call GToParameters115(gy,g1,g2,g3,Yu,Yd,Ye,YtII,YsII,YzII,L1II,L2II)

AbsL1II = Abs(L1II)**2
AbsL2II = Abs(L2II)**2
Call Adjungate(Yu,adjYu)
Call Adjungate(Yd,adjYd)
Call Adjungate(Ye,adjYe)
Call Adjungate(YtII,adjYtII)
Call Adjungate(YsII,adjYsII)
Call Adjungate(YzII,adjYzII)
 YdadjYd = Matmul2(Yd,adjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYd(i2,i2) =  Real(YdadjYd(i2,i2),dp) 
 YeadjYe = Matmul2(Ye,adjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYe(i2,i2) =  Real(YeadjYe(i2,i2),dp) 
 YsIICYsII = Matmul2(YsII,adjYsII,OnlyDiagonal) 
 YtIICYtII = Matmul2(YtII,adjYtII,OnlyDiagonal) 
 YuadjYu = Matmul2(Yu,adjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYu(i2,i2) =  Real(YuadjYu(i2,i2),dp) 
 YzIIadjYzII = Matmul2(YzII,adjYzII,OnlyDiagonal) 
Forall(i2=1:3)  YzIIadjYzII(i2,i2) =  Real(YzIIadjYzII(i2,i2),dp) 
 adjYdYd = Matmul2(adjYd,Yd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYd(i2,i2) =  Real(adjYdYd(i2,i2),dp) 
 adjYdYsII = Matmul2(adjYd,YsII,OnlyDiagonal) 
 adjYdYzII = Matmul2(adjYd,YzII,OnlyDiagonal) 
 adjYeYe = Matmul2(adjYe,Ye,OnlyDiagonal) 
Forall(i2=1:3)  adjYeYe(i2,i2) =  Real(adjYeYe(i2,i2),dp) 
 adjYuYu = Matmul2(adjYu,Yu,OnlyDiagonal) 
Forall(i2=1:3)  adjYuYu(i2,i2) =  Real(adjYuYu(i2,i2),dp) 
 adjYzIIYd = Matmul2(adjYzII,Yd,OnlyDiagonal) 
 adjYzIIYsII = Matmul2(adjYzII,YsII,OnlyDiagonal) 
 adjYzIIYzII = Matmul2(adjYzII,YzII,OnlyDiagonal) 
Forall(i2=1:3)  adjYzIIYzII(i2,i2) =  Real(adjYzIIYzII(i2,i2),dp) 
 CYdTpYd = Matmul2(Conjg(Yd),Transpose(Yd),OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYd(i2,i2) =  Real(CYdTpYd(i2,i2),dp) 
 CYeYtII = Matmul2(Conjg(Ye),YtII,OnlyDiagonal) 
 CYsIIYd = Matmul2(adjYsII,Yd,OnlyDiagonal) 
 CYsIIYsII = Matmul2(adjYsII,YsII,OnlyDiagonal) 
 CYsIIYzII = Matmul2(adjYsII,YzII,OnlyDiagonal) 
 CYtIIYtII = Matmul2(adjYtII,YtII,OnlyDiagonal) 
 CYzIIYtII = Matmul2(Conjg(YzII),YtII,OnlyDiagonal) 
 CYzIITpYzII = Matmul2(Conjg(YzII),Transpose(YzII),OnlyDiagonal) 
Forall(i2=1:3)  CYzIITpYzII(i2,i2) =  Real(CYzIITpYzII(i2,i2),dp) 
 YdadjYdYd = Matmul2(Yd,adjYdYd,OnlyDiagonal) 
 YdadjYdYsII = Matmul2(Yd,adjYdYsII,OnlyDiagonal) 
 YdadjYdYzII = Matmul2(Yd,adjYdYzII,OnlyDiagonal) 
 YdadjYuYu = Matmul2(Yd,adjYuYu,OnlyDiagonal) 
 YeadjYeYe = Matmul2(Ye,adjYeYe,OnlyDiagonal) 
 YeadjYzIIYzII = Matmul2(Ye,adjYzIIYzII,OnlyDiagonal) 
 YeCYtIIYtII = Matmul2(Ye,CYtIIYtII,OnlyDiagonal) 
 YsIICYdTpYd = Matmul2(YsII,CYdTpYd,OnlyDiagonal) 
 YsIICYsIIYd = Matmul2(YsII,CYsIIYd,OnlyDiagonal) 
 YsIICYsIIYsII = Matmul2(YsII,CYsIIYsII,OnlyDiagonal) 
 YsIICYsIIYzII = Matmul2(YsII,CYsIIYzII,OnlyDiagonal) 
 YsIICYzIITpYzII = Matmul2(YsII,CYzIITpYzII,OnlyDiagonal) 
 YtIIadjYeYe = Matmul2(YtII,adjYeYe,OnlyDiagonal) 
 YtIIadjYzIIYzII = Matmul2(YtII,adjYzIIYzII,OnlyDiagonal) 
 YtIICYtIIYtII = Matmul2(YtII,CYtIIYtII,OnlyDiagonal) 
 YuadjYdYd = Matmul2(Yu,adjYdYd,OnlyDiagonal) 
 YuadjYuYu = Matmul2(Yu,adjYuYu,OnlyDiagonal) 
 YzIIadjYeYe = Matmul2(YzII,adjYeYe,OnlyDiagonal) 
 YzIIadjYzIIYd = Matmul2(YzII,adjYzIIYd,OnlyDiagonal) 
 YzIIadjYzIIYsII = Matmul2(YzII,adjYzIIYsII,OnlyDiagonal) 
 YzIIadjYzIIYzII = Matmul2(YzII,adjYzIIYzII,OnlyDiagonal) 
 YzIICYtIIYtII = Matmul2(YzII,CYtIIYtII,OnlyDiagonal) 
 TpYeCYeYtII = Matmul2(Transpose(Ye),CYeYtII,OnlyDiagonal) 
 TpYzIICYzIIYtII = Matmul2(Transpose(YzII),CYzIIYtII,OnlyDiagonal) 
 TrYdadjYd = Real(cTrace(YdadjYd),dp) 
 TrYeadjYe = Real(cTrace(YeadjYe),dp) 
 TrYsIICYsII = Real(cTrace(YsIICYsII),dp) 
 TrYtIICYtII = Real(cTrace(YtIICYtII),dp) 
 TrYuadjYu = Real(cTrace(YuadjYu),dp) 
 TrYzIIadjYzII = Real(cTrace(YzIIadjYzII),dp) 
 g1p2 =g1**2 
 g1p3 =g1**3 
 g2p2 =g2**2 
 g2p3 =g2**3 
 g3p2 =g3**2 
 g3p3 =g3**3 
 L1IIp2 =L1II**2 
 L2IIp2 =L2II**2 
 g1p4 =g1**4 
 g2p4 =g2**4 
 g3p4 =g3**4 
 CL1IIp2 =Conjg(L1II)**2 
 CL2IIp2 =Conjg(L2II)**2 


If (TwoLoopRGE) Then 
 YtIIadjYe = Matmul2(YtII,adjYe,OnlyDiagonal) 
 YuadjYd = Matmul2(Yu,adjYd,OnlyDiagonal) 
 YzIIadjYe = Matmul2(YzII,adjYe,OnlyDiagonal) 
 YzIICYtII = Matmul2(YzII,adjYtII,OnlyDiagonal) 
 CYeTpYzII = Matmul2(Conjg(Ye),Transpose(YzII),OnlyDiagonal) 
 CYtIITpYzII = Matmul2(adjYtII,Transpose(YzII),OnlyDiagonal) 
 CYuTpYd = Matmul2(Conjg(Yu),Transpose(Yd),OnlyDiagonal) 
 YeadjYzIIYd = Matmul2(Ye,adjYzIIYd,OnlyDiagonal) 
 YeadjYzIIYsII = Matmul2(Ye,adjYzIIYsII,OnlyDiagonal) 
 YsIICYzIIYtII = Matmul2(YsII,CYzIIYtII,OnlyDiagonal) 
 YtIIadjYzIIYd = Matmul2(YtII,adjYzIIYd,OnlyDiagonal) 
 YtIIadjYzIIYsII = Matmul2(YtII,adjYzIIYsII,OnlyDiagonal) 
 YtIICYtIITpYzII = Matmul2(YtII,CYtIITpYzII,OnlyDiagonal) 
 YuadjYdYsII = Matmul2(Yu,adjYdYsII,OnlyDiagonal) 
 YuadjYdYzII = Matmul2(Yu,adjYdYzII,OnlyDiagonal) 
 adjYdYdadjYd = Matmul2(adjYd,YdadjYd,OnlyDiagonal) 
 adjYdYsIICYsII = Matmul2(adjYd,YsIICYsII,OnlyDiagonal) 
 adjYdYzIIadjYzII = Matmul2(adjYd,YzIIadjYzII,OnlyDiagonal) 
 adjYeYeadjYe = Matmul2(adjYe,YeadjYe,OnlyDiagonal) 
 adjYuYuadjYd = Matmul2(adjYu,YuadjYd,OnlyDiagonal) 
 adjYuYuadjYu = Matmul2(adjYu,YuadjYu,OnlyDiagonal) 
 adjYzIIYzIIadjYe = Matmul2(adjYzII,YzIIadjYe,OnlyDiagonal) 
 adjYzIIYzIIadjYzII = Matmul2(adjYzII,YzIIadjYzII,OnlyDiagonal) 
 adjYzIIYzIICYtII = Matmul2(adjYzII,YzIICYtII,OnlyDiagonal) 
 CYsIIYsIICYsII = Matmul2(adjYsII,YsIICYsII,OnlyDiagonal) 
 CYsIIYzIIadjYzII = Matmul2(adjYsII,YzIIadjYzII,OnlyDiagonal) 
 CYtIIYtIIadjYe = Matmul2(adjYtII,YtIIadjYe,OnlyDiagonal) 
 CYtIIYtIICYtII = Matmul2(adjYtII,YtIICYtII,OnlyDiagonal) 
 TpYdCYdTpYd = Matmul2(Transpose(Yd),CYdTpYd,OnlyDiagonal) 
 TpYdCYsIIYd = Matmul2(Transpose(Yd),CYsIIYd,OnlyDiagonal) 
 TpYdCYsIIYsII = Matmul2(Transpose(Yd),CYsIIYsII,OnlyDiagonal) 
 TpYdCYsIIYzII = Matmul2(Transpose(Yd),CYsIIYzII,OnlyDiagonal) 
 TpYdCYzIIYtII = Matmul2(Transpose(Yd),CYzIIYtII,OnlyDiagonal) 
 TpYeCYeTpYzII = Matmul2(Transpose(Ye),CYeTpYzII,OnlyDiagonal) 
 TpYuCYuTpYd = Matmul2(Transpose(Yu),CYuTpYd,OnlyDiagonal) 
 TpYzIICYsIIYd = Matmul2(Transpose(YzII),CYsIIYd,OnlyDiagonal) 
 TpYzIICYsIIYsII = Matmul2(Transpose(YzII),CYsIIYsII,OnlyDiagonal) 
 TpYzIICYsIIYzII = Matmul2(Transpose(YzII),CYsIIYzII,OnlyDiagonal) 
 TpYzIICYzIITpYzII = Matmul2(Transpose(YzII),CYzIITpYzII,OnlyDiagonal) 
 YdadjYdYdadjYd = Matmul2(Yd,adjYdYdadjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYdYdadjYd(i2,i2) =  Real(YdadjYdYdadjYd(i2,i2),dp) 
 YdadjYdYsIICYsII = Matmul2(Yd,adjYdYsIICYsII,OnlyDiagonal) 
 YdadjYdYzIIadjYzII = Matmul2(Yd,adjYdYzIIadjYzII,OnlyDiagonal) 
 YdadjYuYuadjYd = Matmul2(Yd,adjYuYuadjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYuYuadjYd(i2,i2) =  Real(YdadjYuYuadjYd(i2,i2),dp) 
 YeadjYeYeadjYe = Matmul2(Ye,adjYeYeadjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYeYeadjYe(i2,i2) =  Real(YeadjYeYeadjYe(i2,i2),dp) 
 YeadjYzIIYzIIadjYe = Matmul2(Ye,adjYzIIYzIIadjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYzIIYzIIadjYe(i2,i2) =  Real(YeadjYzIIYzIIadjYe(i2,i2),dp) 
 YeCYtIIYtIIadjYe = Matmul2(Ye,CYtIIYtIIadjYe,OnlyDiagonal) 
 YsIICYsIIYsIICYsII = Matmul2(YsII,CYsIIYsIICYsII,OnlyDiagonal) 
 YsIICYsIIYzIIadjYzII = Matmul2(YsII,CYsIIYzIIadjYzII,OnlyDiagonal) 
 YtIIadjYzIIYzIICYtII = Matmul2(YtII,adjYzIIYzIICYtII,OnlyDiagonal) 
 YtIICYtIIYtIICYtII = Matmul2(YtII,CYtIIYtIICYtII,OnlyDiagonal) 
 YuadjYuYuadjYu = Matmul2(Yu,adjYuYuadjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYuYuadjYu(i2,i2) =  Real(YuadjYuYuadjYu(i2,i2),dp) 
 YzIIadjYzIIYzIIadjYzII = Matmul2(YzII,adjYzIIYzIIadjYzII,OnlyDiagonal) 
Forall(i2=1:3)  YzIIadjYzIIYzIIadjYzII(i2,i2) =  Real(YzIIadjYzIIYzIIadjYzII(i2,i2),dp) 
 adjYdYdadjYdYd = Matmul2(adjYd,YdadjYdYd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYdadjYdYd(i2,i2) =  Real(adjYdYdadjYdYd(i2,i2),dp) 
 adjYdYdadjYdYsII = Matmul2(adjYd,YdadjYdYsII,OnlyDiagonal) 
 adjYdYdadjYdYzII = Matmul2(adjYd,YdadjYdYzII,OnlyDiagonal) 
 adjYdYdadjYuYu = Matmul2(adjYd,YdadjYuYu,OnlyDiagonal) 
 adjYdYsIICYsIIYd = Matmul2(adjYd,YsIICYsIIYd,OnlyDiagonal) 
 adjYdYzIIadjYzIIYd = Matmul2(adjYd,YzIIadjYzIIYd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYzIIadjYzIIYd(i2,i2) =  Real(adjYdYzIIadjYzIIYd(i2,i2),dp) 
 adjYeYeadjYeYe = Matmul2(adjYe,YeadjYeYe,OnlyDiagonal) 
Forall(i2=1:3)  adjYeYeadjYeYe(i2,i2) =  Real(adjYeYeadjYeYe(i2,i2),dp) 
 adjYeYeadjYzIIYd = Matmul2(adjYe,YeadjYzIIYd,OnlyDiagonal) 
 adjYeYeadjYzIIYsII = Matmul2(adjYe,YeadjYzIIYsII,OnlyDiagonal) 
 adjYeYeadjYzIIYzII = Matmul2(adjYe,YeadjYzIIYzII,OnlyDiagonal) 
 adjYeYeCYtIIYtII = Matmul2(adjYe,YeCYtIIYtII,OnlyDiagonal) 
 adjYuYuadjYdYd = Matmul2(adjYu,YuadjYdYd,OnlyDiagonal) 
 adjYuYuadjYdYsII = Matmul2(adjYu,YuadjYdYsII,OnlyDiagonal) 
 adjYuYuadjYdYzII = Matmul2(adjYu,YuadjYdYzII,OnlyDiagonal) 
 adjYuYuadjYuYu = Matmul2(adjYu,YuadjYuYu,OnlyDiagonal) 
Forall(i2=1:3)  adjYuYuadjYuYu(i2,i2) =  Real(adjYuYuadjYuYu(i2,i2),dp) 
 adjYzIIYdadjYdYzII = Matmul2(adjYzII,YdadjYdYzII,OnlyDiagonal) 
Forall(i2=1:3)  adjYzIIYdadjYdYzII(i2,i2) =  Real(adjYzIIYdadjYdYzII(i2,i2),dp) 
 adjYzIIYsIICYsIIYzII = Matmul2(adjYzII,YsIICYsIIYzII,OnlyDiagonal) 
 adjYzIIYzIIadjYeYe = Matmul2(adjYzII,YzIIadjYeYe,OnlyDiagonal) 
 adjYzIIYzIIadjYzIIYd = Matmul2(adjYzII,YzIIadjYzIIYd,OnlyDiagonal) 
 adjYzIIYzIIadjYzIIYsII = Matmul2(adjYzII,YzIIadjYzIIYsII,OnlyDiagonal) 
 adjYzIIYzIIadjYzIIYzII = Matmul2(adjYzII,YzIIadjYzIIYzII,OnlyDiagonal) 
Forall(i2=1:3)  adjYzIIYzIIadjYzIIYzII(i2,i2) =  Real(adjYzIIYzIIadjYzIIYzII(i2,i2),dp) 
 adjYzIIYzIICYtIIYtII = Matmul2(adjYzII,YzIICYtIIYtII,OnlyDiagonal) 
 CYdTpYdCYdTpYd = Matmul2(Conjg(Yd),TpYdCYdTpYd,OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYdCYdTpYd(i2,i2) =  Real(CYdTpYdCYdTpYd(i2,i2),dp) 
 CYdTpYdCYsIIYd = Matmul2(Conjg(Yd),TpYdCYsIIYd,OnlyDiagonal) 
 CYdTpYdCYsIIYsII = Matmul2(Conjg(Yd),TpYdCYsIIYsII,OnlyDiagonal) 
 CYdTpYdCYsIIYzII = Matmul2(Conjg(Yd),TpYdCYsIIYzII,OnlyDiagonal) 
 CYdTpYdCYzIIYtII = Matmul2(Conjg(Yd),TpYdCYzIIYtII,OnlyDiagonal) 
 CYdTpYuCYuTpYd = Matmul2(Conjg(Yd),TpYuCYuTpYd,OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYuCYuTpYd(i2,i2) =  Real(CYdTpYuCYuTpYd(i2,i2),dp) 
 CYeTpYeCYeYtII = Matmul2(Conjg(Ye),TpYeCYeYtII,OnlyDiagonal) 
 CYsIIYdadjYdYsII = Matmul2(adjYsII,YdadjYdYsII,OnlyDiagonal) 
 CYsIIYsIICYsIIYd = Matmul2(adjYsII,YsIICYsIIYd,OnlyDiagonal) 
 CYsIIYsIICYsIIYsII = Matmul2(adjYsII,YsIICYsIIYsII,OnlyDiagonal) 
 CYsIIYsIICYsIIYzII = Matmul2(adjYsII,YsIICYsIIYzII,OnlyDiagonal) 
 CYsIIYsIICYzIIYtII = Matmul2(adjYsII,YsIICYzIIYtII,OnlyDiagonal) 
 CYsIIYzIIadjYzIIYsII = Matmul2(adjYsII,YzIIadjYzIIYsII,OnlyDiagonal) 
 CYtIIYtIIadjYeYe = Matmul2(adjYtII,YtIIadjYeYe,OnlyDiagonal) 
 CYtIIYtIIadjYzIIYd = Matmul2(adjYtII,YtIIadjYzIIYd,OnlyDiagonal) 
 CYtIIYtIIadjYzIIYsII = Matmul2(adjYtII,YtIIadjYzIIYsII,OnlyDiagonal) 
 CYtIIYtIIadjYzIIYzII = Matmul2(adjYtII,YtIIadjYzIIYzII,OnlyDiagonal) 
 CYtIIYtIICYtIIYtII = Matmul2(adjYtII,YtIICYtIIYtII,OnlyDiagonal) 
 CYtIITpYeCYeYtII = Matmul2(adjYtII,TpYeCYeYtII,OnlyDiagonal) 
 CYtIITpYzIICYzIIYtII = Matmul2(adjYtII,TpYzIICYzIIYtII,OnlyDiagonal) 
 CYzIIYtIICYtIITpYzII = Matmul2(Conjg(YzII),YtIICYtIITpYzII,OnlyDiagonal) 
 CYzIITpYeCYeTpYzII = Matmul2(Conjg(YzII),TpYeCYeTpYzII,OnlyDiagonal) 
Forall(i2=1:3)  CYzIITpYeCYeTpYzII(i2,i2) =  Real(CYzIITpYeCYeTpYzII(i2,i2),dp) 
 CYzIITpYzIICYsIIYd = Matmul2(Conjg(YzII),TpYzIICYsIIYd,OnlyDiagonal) 
 CYzIITpYzIICYsIIYsII = Matmul2(Conjg(YzII),TpYzIICYsIIYsII,OnlyDiagonal) 
 CYzIITpYzIICYsIIYzII = Matmul2(Conjg(YzII),TpYzIICYsIIYzII,OnlyDiagonal) 
 CYzIITpYzIICYzIIYtII = Matmul2(Conjg(YzII),TpYzIICYzIIYtII,OnlyDiagonal) 
 CYzIITpYzIICYzIITpYzII = Matmul2(Conjg(YzII),TpYzIICYzIITpYzII,OnlyDiagonal) 
Forall(i2=1:3)  CYzIITpYzIICYzIITpYzII(i2,i2) =  Real(CYzIITpYzIICYzIITpYzII(i2,i2),dp) 
 YdadjYdYdadjYdYd = Matmul2(Yd,adjYdYdadjYdYd,OnlyDiagonal) 
 YdadjYdYdadjYdYsII = Matmul2(Yd,adjYdYdadjYdYsII,OnlyDiagonal) 
 YdadjYdYdadjYdYzII = Matmul2(Yd,adjYdYdadjYdYzII,OnlyDiagonal) 
 YdadjYdYsIICYsIIYd = Matmul2(Yd,adjYdYsIICYsIIYd,OnlyDiagonal) 
 YdadjYdYzIIadjYzIIYd = Matmul2(Yd,adjYdYzIIadjYzIIYd,OnlyDiagonal) 
 YdadjYuYuadjYdYd = Matmul2(Yd,adjYuYuadjYdYd,OnlyDiagonal) 
 YdadjYuYuadjYdYsII = Matmul2(Yd,adjYuYuadjYdYsII,OnlyDiagonal) 
 YdadjYuYuadjYdYzII = Matmul2(Yd,adjYuYuadjYdYzII,OnlyDiagonal) 
 YdadjYuYuadjYuYu = Matmul2(Yd,adjYuYuadjYuYu,OnlyDiagonal) 
 YeadjYeYeadjYeYe = Matmul2(Ye,adjYeYeadjYeYe,OnlyDiagonal) 
 YeadjYzIIYdadjYdYzII = Matmul2(Ye,adjYzIIYdadjYdYzII,OnlyDiagonal) 
 YeadjYzIIYsIICYsIIYzII = Matmul2(Ye,adjYzIIYsIICYsIIYzII,OnlyDiagonal) 
 YeadjYzIIYzIIadjYeYe = Matmul2(Ye,adjYzIIYzIIadjYeYe,OnlyDiagonal) 
 YeadjYzIIYzIIadjYzIIYzII = Matmul2(Ye,adjYzIIYzIIadjYzIIYzII,OnlyDiagonal) 
 YeCYtIIYtIIadjYeYe = Matmul2(Ye,CYtIIYtIIadjYeYe,OnlyDiagonal) 
 YeCYtIIYtIICYtIIYtII = Matmul2(Ye,CYtIIYtIICYtIIYtII,OnlyDiagonal) 
 YeCYtIITpYeCYeYtII = Matmul2(Ye,CYtIITpYeCYeYtII,OnlyDiagonal) 
 YeCYtIITpYzIICYzIIYtII = Matmul2(Ye,CYtIITpYzIICYzIIYtII,OnlyDiagonal) 
 YsIICYdTpYdCYdTpYd = Matmul2(YsII,CYdTpYdCYdTpYd,OnlyDiagonal) 
 YsIICYdTpYdCYsIIYd = Matmul2(YsII,CYdTpYdCYsIIYd,OnlyDiagonal) 
 YsIICYdTpYdCYsIIYsII = Matmul2(YsII,CYdTpYdCYsIIYsII,OnlyDiagonal) 
 YsIICYdTpYdCYsIIYzII = Matmul2(YsII,CYdTpYdCYsIIYzII,OnlyDiagonal) 
 YsIICYdTpYuCYuTpYd = Matmul2(YsII,CYdTpYuCYuTpYd,OnlyDiagonal) 
 YsIICYsIIYdadjYdYsII = Matmul2(YsII,CYsIIYdadjYdYsII,OnlyDiagonal) 
 YsIICYsIIYsIICYsIIYd = Matmul2(YsII,CYsIIYsIICYsIIYd,OnlyDiagonal) 
 YsIICYsIIYsIICYsIIYsII = Matmul2(YsII,CYsIIYsIICYsIIYsII,OnlyDiagonal) 
 YsIICYsIIYsIICYsIIYzII = Matmul2(YsII,CYsIIYsIICYsIIYzII,OnlyDiagonal) 
 YsIICYsIIYzIIadjYzIIYsII = Matmul2(YsII,CYsIIYzIIadjYzIIYsII,OnlyDiagonal) 
 YsIICYzIIYtIICYtIITpYzII = Matmul2(YsII,CYzIIYtIICYtIITpYzII,OnlyDiagonal) 
 YsIICYzIITpYeCYeTpYzII = Matmul2(YsII,CYzIITpYeCYeTpYzII,OnlyDiagonal) 
 YsIICYzIITpYzIICYsIIYd = Matmul2(YsII,CYzIITpYzIICYsIIYd,OnlyDiagonal) 
 YsIICYzIITpYzIICYsIIYsII = Matmul2(YsII,CYzIITpYzIICYsIIYsII,OnlyDiagonal) 
 YsIICYzIITpYzIICYsIIYzII = Matmul2(YsII,CYzIITpYzIICYsIIYzII,OnlyDiagonal) 
 YsIICYzIITpYzIICYzIITpYzII = Matmul2(YsII,CYzIITpYzIICYzIITpYzII,OnlyDiagonal) 
 YtIIadjYeYeadjYeYe = Matmul2(YtII,adjYeYeadjYeYe,OnlyDiagonal) 
 YtIIadjYeYeCYtIIYtII = Matmul2(YtII,adjYeYeCYtIIYtII,OnlyDiagonal) 
 YtIIadjYzIIYdadjYdYzII = Matmul2(YtII,adjYzIIYdadjYdYzII,OnlyDiagonal) 
 YtIIadjYzIIYsIICYsIIYzII = Matmul2(YtII,adjYzIIYsIICYsIIYzII,OnlyDiagonal) 
 YtIIadjYzIIYzIIadjYzIIYzII = Matmul2(YtII,adjYzIIYzIIadjYzIIYzII,OnlyDiagonal) 
 YtIIadjYzIIYzIICYtIIYtII = Matmul2(YtII,adjYzIIYzIICYtIIYtII,OnlyDiagonal) 
 YtIICYtIIYtIICYtIIYtII = Matmul2(YtII,CYtIIYtIICYtIIYtII,OnlyDiagonal) 
 YtIICYtIITpYeCYeYtII = Matmul2(YtII,CYtIITpYeCYeYtII,OnlyDiagonal) 
 YtIICYtIITpYzIICYzIIYtII = Matmul2(YtII,CYtIITpYzIICYzIIYtII,OnlyDiagonal) 
 YuadjYdYdadjYdYd = Matmul2(Yu,adjYdYdadjYdYd,OnlyDiagonal) 
 YuadjYdYdadjYuYu = Matmul2(Yu,adjYdYdadjYuYu,OnlyDiagonal) 
 YuadjYdYsIICYsIIYd = Matmul2(Yu,adjYdYsIICYsIIYd,OnlyDiagonal) 
 YuadjYdYzIIadjYzIIYd = Matmul2(Yu,adjYdYzIIadjYzIIYd,OnlyDiagonal) 
 YuadjYuYuadjYuYu = Matmul2(Yu,adjYuYuadjYuYu,OnlyDiagonal) 
 YzIIadjYeYeadjYeYe = Matmul2(YzII,adjYeYeadjYeYe,OnlyDiagonal) 
 YzIIadjYeYeadjYzIIYd = Matmul2(YzII,adjYeYeadjYzIIYd,OnlyDiagonal) 
 YzIIadjYeYeadjYzIIYsII = Matmul2(YzII,adjYeYeadjYzIIYsII,OnlyDiagonal) 
 YzIIadjYeYeadjYzIIYzII = Matmul2(YzII,adjYeYeadjYzIIYzII,OnlyDiagonal) 
 YzIIadjYzIIYdadjYdYzII = Matmul2(YzII,adjYzIIYdadjYdYzII,OnlyDiagonal) 
 YzIIadjYzIIYsIICYsIIYzII = Matmul2(YzII,adjYzIIYsIICYsIIYzII,OnlyDiagonal) 
 YzIIadjYzIIYzIIadjYzIIYd = Matmul2(YzII,adjYzIIYzIIadjYzIIYd,OnlyDiagonal) 
 YzIIadjYzIIYzIIadjYzIIYsII = Matmul2(YzII,adjYzIIYzIIadjYzIIYsII,OnlyDiagonal) 
 YzIIadjYzIIYzIIadjYzIIYzII = Matmul2(YzII,adjYzIIYzIIadjYzIIYzII,OnlyDiagonal) 
 YzIICYtIIYtIIadjYzIIYd = Matmul2(YzII,CYtIIYtIIadjYzIIYd,OnlyDiagonal) 
 YzIICYtIIYtIIadjYzIIYsII = Matmul2(YzII,CYtIIYtIIadjYzIIYsII,OnlyDiagonal) 
 YzIICYtIIYtIIadjYzIIYzII = Matmul2(YzII,CYtIIYtIIadjYzIIYzII,OnlyDiagonal) 
 YzIICYtIIYtIICYtIIYtII = Matmul2(YzII,CYtIIYtIICYtIIYtII,OnlyDiagonal) 
 YzIICYtIITpYeCYeYtII = Matmul2(YzII,CYtIITpYeCYeYtII,OnlyDiagonal) 
 YzIICYtIITpYzIICYzIIYtII = Matmul2(YzII,CYtIITpYzIICYzIIYtII,OnlyDiagonal) 
 TpYeCYeTpYeCYeYtII = Matmul2(Transpose(Ye),CYeTpYeCYeYtII,OnlyDiagonal) 
 TpYzIICYdTpYdCYzIIYtII = Matmul2(Transpose(YzII),CYdTpYdCYzIIYtII,OnlyDiagonal) 
 TpYzIICYsIIYsIICYzIIYtII = Matmul2(Transpose(YzII),CYsIIYsIICYzIIYtII,OnlyDiagonal) 
 TpYzIICYzIITpYzIICYzIIYtII = Matmul2(Transpose(YzII),CYzIITpYzIICYzIIYtII,OnlyDiagonal) 
 TrYdadjYdYdadjYd = cTrace(YdadjYdYdadjYd) 
 TrYdadjYdYsIICYsII = cTrace(YdadjYdYsIICYsII) 
 TrYdadjYdYzIIadjYzII = cTrace(YdadjYdYzIIadjYzII) 
 TrYdadjYuYuadjYd = cTrace(YdadjYuYuadjYd) 
 TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe) 
 TrYeadjYzIIYzIIadjYe = cTrace(YeadjYzIIYzIIadjYe) 
 TrYeCYtIIYtIIadjYe = cTrace(YeCYtIIYtIIadjYe) 
 TrYsIICYsIIYsIICYsII = cTrace(YsIICYsIIYsIICYsII) 
 TrYsIICYsIIYzIIadjYzII = cTrace(YsIICYsIIYzIIadjYzII) 
 TrYtIIadjYzIIYzIICYtII = cTrace(YtIIadjYzIIYzIICYtII) 
 TrYtIICYtIIYtIICYtII = cTrace(YtIICYtIIYtIICYtII) 
 TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu) 
 TrYzIIadjYzIIYzIIadjYzII = cTrace(YzIIadjYzIIYzIIadjYzII) 
 g1p4 =g1**4 
 g2p4 =g2**4 
 g3p4 =g3**4 
 CL1IIp2 =Conjg(L1II)**2 
 CL2IIp2 =Conjg(L2II)**2 
End If 
 
 
!-------------------- 
! g1 
!-------------------- 
 
betag11  = 68._dp*(g1p3)/5._dp

 
 
If (TwoLoopRGE) Then 
betag12 = (g1p3*(-405._dp*(AbsL1II) - 405._dp*(AbsL2II) + 1502._dp*(g1p2) + 2610._dp*(g2p2) +   & 
&  4600._dp*(g3p2) - 210._dp*(TrYdadjYd) - 270._dp*(TrYeadjYe) - 360._dp*(TrYsIICYsII) - & 
&  405._dp*(TrYtIICYtII) - 390._dp*(TrYuadjYu) - 210._dp*(TrYzIIadjYzII)))/75._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21  = 8._dp*(g2p3)

 
 
If (TwoLoopRGE) Then 
betag22 = (g2p3*(-35._dp*(AbsL1II) - 35._dp*(AbsL2II) + 58._dp*(g1p2) + 470._dp*(g2p2) +        & 
&  200._dp*(g3p2) - 30._dp*(TrYdadjYd) - 10._dp*(TrYeadjYe) - 35._dp*(TrYtIICYtII) -     & 
&  30._dp*(TrYuadjYu) - 30._dp*(TrYzIIadjYzII)))/5._dp

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31  = 4._dp*(g3p3)

 
 
If (TwoLoopRGE) Then 
betag32 = (g3p3*(23._dp*(g1p2) + 45._dp*(g2p2) + 400._dp*(g3p2) - 12._dp*(TrYdadjYd) -          & 
&  27._dp*(TrYsIICYsII) - 12._dp*(TrYuadjYu) - 12._dp*(TrYzIIadjYzII)))/3._dp

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yu 
!-------------------- 
 
betaYu1  = (3._dp*(AbsL2II) - 13._dp*(g1p2)/15._dp - 3._dp*(g2p2) - 16._dp*(g3p2)     & 
& /3._dp + 3._dp*(TrYuadjYu))*Yu + YuadjYdYd + 3._dp*(YuadjYuYu)

 
 
If (TwoLoopRGE) Then 
betaYu2 = (5473._dp*(g1p4)/450._dp + g1p2*g2p2 + 57._dp*(g2p4)/2._dp + (136*g1p2*g3p2)/45._dp + & 
&  8*g2p2*g3p2 + 320._dp*(g3p4)/9._dp - 12*CL2IIp2*L2IIp2 - 3._dp*(TrYdadjYuYuadjYd) +   & 
&  (3*AbsL2II*(6._dp*(g1p2) + 20._dp*(g2p2) - 15._dp*(TrYuadjYu)))/5._dp +               & 
&  (4*(g1p2 + 20._dp*(g3p2))*TrYuadjYu)/5._dp - 9._dp*(TrYuadjYuYuadjYu))*Yu +           & 
&  (-3._dp*(AbsL1II) + 2._dp*(g1p2)/5._dp - 3._dp*(TrYdadjYd) - TrYeadjYe)*YuadjYdYd -   & 
&  2._dp*(YuadjYdYdadjYdYd) - 2._dp*(YuadjYdYdadjYuYu) - 4._dp*(YuadjYdYsIICYsIIYd) -    & 
&  2._dp*(YuadjYdYzIIadjYzIIYd) - 9*AbsL2II*YuadjYuYu + (2*g1p2*YuadjYuYu)/5._dp +       & 
&  6*g2p2*YuadjYuYu - 9*TrYuadjYu*YuadjYuYu - 4._dp*(YuadjYuYuadjYuYu)

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
betaYd1  = (3._dp*(AbsL1II) - 7._dp*(g1p2)/15._dp - 3._dp*(g2p2) - 16._dp*(g3p2)      & 
& /3._dp + 3._dp*(TrYdadjYd) + TrYeadjYe)*Yd + 3._dp*(YdadjYdYd) + YdadjYuYu +           & 
&  4._dp*(YsIICYsIIYd) + 2._dp*(YzIIadjYzIIYd)

 
 
If (TwoLoopRGE) Then 
betaYd2 = (581._dp*(g1p4)/90._dp + g1p2*g2p2 + 57._dp*(g2p4)/2._dp + (8*g1p2*g3p2)/9._dp +      & 
&  8*g2p2*g3p2 + 320._dp*(g3p4)/9._dp - 12*CL1IIp2*L1IIp2 - (2*(g1p2 - 40._dp*(g3p2))*TrYdadjYd)/5._dp -& 
&  9._dp*(TrYdadjYdYdadjYd) - 12._dp*(TrYdadjYdYsIICYsII) - 6._dp*(TrYdadjYdYzIIadjYzII) -& 
&  3._dp*(TrYdadjYuYuadjYd) + (6*g1p2*TrYeadjYe)/5._dp - 3._dp*(TrYeadjYeYeadjYe) -      & 
&  3._dp*(TrYeadjYzIIYzIIadjYe) - 3._dp*(TrYeCYtIIYtIIadjYe) + (3*AbsL1II*(6._dp*(g1p2) +& 
&  20._dp*(g2p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe) - 5._dp*(TrYtIICYtII)))/5._dp)*Yd +& 
&  (-9._dp*(AbsL1II) + 4._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 9._dp*(TrYdadjYd) -           & 
&  3._dp*(TrYeadjYe))*YdadjYdYd - 4._dp*(YdadjYdYdadjYdYd) - 4._dp*(YdadjYdYsIICYsIIYd) -& 
&  2._dp*(YdadjYdYzIIadjYzIIYd) - 3*AbsL2II*YdadjYuYu + (4*g1p2*YdadjYuYu)/5._dp -       & 
&  3*TrYuadjYu*YdadjYuYu - 2._dp*(YdadjYuYuadjYdYd) - 2._dp*(YdadjYuYuadjYuYu) -         & 
&  8._dp*(YsIICYdTpYdCYsIIYd) + (32*g1p2*YsIICYsIIYd)/15._dp + (80*g3p2*YsIICYsIIYd)/3._dp -& 
&  4*TrYsIICYsII*YsIICYsIIYd - 16._dp*(YsIICYsIIYsIICYsIIYd) - 8._dp*(YsIICYzIITpYzIICYsIIYd) -& 
&  2._dp*(YzIIadjYeYeadjYzIIYd) + (2*g1p2*YzIIadjYzIIYd)/5._dp + 6*g2p2*YzIIadjYzIIYd -  & 
&  2*TrYzIIadjYzII*YzIIadjYzIIYd - 6._dp*(YzIIadjYzIIYzIIadjYzIIYd) - 6._dp*(YzIICYtIIYtIIadjYzIIYd)

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
!-------------------- 
! Ye 
!-------------------- 
 
betaYe1  = (3._dp*(AbsL1II) - 9._dp*(g1p2)/5._dp - 3._dp*(g2p2) + 3._dp*(TrYdadjYd)   & 
&  + TrYeadjYe)*Ye + 3*(YeadjYeYe + YeadjYzIIYzII + YeCYtIIYtII)

 
 
If (TwoLoopRGE) Then 
betaYe2 = -((120*CL1IIp2*L1IIp2 + 4*(g1p2 - 40._dp*(g3p2))*TrYdadjYd + 3*(-87._dp*(g1p4) -      & 
&  6*g1p2*g2p2 - 95._dp*(g2p4) + 30._dp*(TrYdadjYdYdadjYd) + 40._dp*(TrYdadjYdYsIICYsII) +& 
&  20._dp*(TrYdadjYdYzIIadjYzII) + 10._dp*(TrYdadjYuYuadjYd) - 4*g1p2*TrYeadjYe +        & 
&  10._dp*(TrYeadjYeYeadjYe) + 10._dp*(TrYeadjYzIIYzIIadjYe) + 10._dp*(TrYeCYtIIYtIIadjYe)) -& 
&  6*AbsL1II*(6._dp*(g1p2) + 20._dp*(g2p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe) -    & 
&  5._dp*(TrYtIICYtII)))*Ye)/10._dp + (-9._dp*(AbsL1II) + 6._dp*(g2p2) - 9._dp*(TrYdadjYd) -& 
&  3._dp*(TrYeadjYe))*YeadjYeYe - 4._dp*(YeadjYeYeadjYeYe) - 6._dp*(YeadjYzIIYdadjYdYzII) -& 
&  12._dp*(YeadjYzIIYsIICYsIIYzII) - (2*g1p2*YeadjYzIIYzII)/5._dp + 16*g3p2*YeadjYzIIYzII -& 
&  3*TrYzIIadjYzII*YeadjYzIIYzII - 6._dp*(YeadjYzIIYzIIadjYeYe) - 6._dp*(YeadjYzIIYzIIadjYzIIYzII) -& 
&  3._dp*(YeCYtIITpYeCYeYtII) - 9._dp*(YeCYtIITpYzIICYzIIYtII) - 3*AbsL1II*YeCYtIIYtII + & 
&  (18*g1p2*YeCYtIIYtII)/5._dp + 12*g2p2*YeCYtIIYtII - 3*TrYtIICYtII*YeCYtIIYtII -       & 
&  6._dp*(YeCYtIIYtIIadjYeYe) - 9._dp*(YeCYtIIYtIICYtIIYtII)

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
!-------------------- 
! YtII 
!-------------------- 
 
betaYtII1  = TpYeCYeYtII + 3._dp*(TpYzIICYzIIYtII) + (AbsL1II - 9._dp*(g1p2)          & 
& /5._dp - 7._dp*(g2p2) + TrYtIICYtII)*YtII + YtIIadjYeYe + 3._dp*(YtIIadjYzIIYzII)      & 
&  + 6._dp*(YtIICYtIIYtII)

 
 
If (TwoLoopRGE) Then 
betaYtII2 = ((261._dp*(g1p4) + 114*g1p2*g2p2 + 765._dp*(g2p4) - 60*CL1IIp2*L1IIp2 -               & 
&  2*AbsL1II*(3._dp*(g1p2) + 5._dp*(g2p2) + 30._dp*(TrYdadjYd) + 10._dp*(TrYeadjYe)) -   & 
&  20._dp*(TrYeCYtIIYtIIadjYe) - 60._dp*(TrYtIIadjYzIIYzIICYtII) - 2*(3._dp*(g1p2) +     & 
&  5._dp*(g2p2))*TrYtIICYtII - 60._dp*(TrYtIICYtIIYtIICYtII))*YtII)/10._dp +             & 
&  (-3._dp*(AbsL1II) + 6._dp*(g1p2)/5._dp - 3._dp*(TrYdadjYd) - TrYeadjYe)*YtIIadjYeYe + & 
&  (-10._dp*(TpYeCYeTpYeCYeYtII) - 15*AbsL1II*TpYeCYeYtII + 6*g1p2*TpYeCYeYtII -         & 
&  30._dp*(TpYzIICYdTpYdCYzIIYtII) - 60._dp*(TpYzIICYsIIYsIICYzIIYtII) - 30._dp*(TpYzIICYzIITpYzIICYzIIYtII) -& 
&  2*g1p2*TpYzIICYzIIYtII + 80*g3p2*TpYzIICYzIIYtII - 15*TpYeCYeYtII*TrYdadjYd -         & 
&  5*TpYeCYeYtII*TrYeadjYe - 15*TpYzIICYzIIYtII*TrYzIIadjYzII - 10._dp*(YtIIadjYeYeadjYeYe) -& 
&  15._dp*(YtIIadjYeYeCYtIIYtII) - 30._dp*(YtIIadjYzIIYdadjYdYzII) - 60._dp*(YtIIadjYzIIYsIICYsIIYzII) +& 
&  (-2._dp*(g1p2) + 80._dp*(g3p2) - 15._dp*(TrYzIIadjYzII))*YtIIadjYzIIYzII -            & 
&  30._dp*(YtIIadjYzIIYzIIadjYzIIYzII) - 45._dp*(YtIIadjYzIIYzIICYtIIYtII) -             & 
&  15._dp*(YtIICYtIITpYeCYeYtII) - 45._dp*(YtIICYtIITpYzIICYzIIYtII) + 6*(-              & 
& 5._dp*(AbsL1II) + 6._dp*(g1p2) + 20._dp*(g2p2) - 5._dp*(TrYtIICYtII))*YtIICYtIIYtII -  & 
&  90._dp*(YtIICYtIIYtIICYtIIYtII))/5._dp

 
DYtII = oo16pi2*( betaYtII1 + oo16pi2 * betaYtII2 ) 

 
Else 
DYtII = oo16pi2* betaYtII1 
End If 
 
 
!-------------------- 
! YsII 
!-------------------- 
 
betaYsII1  = ((-4*(g1p2 + 15._dp*(g3p2)))/5._dp + TrYsIICYsII)*YsII + 2*(YdadjYdYsII +& 
&  YsIICYdTpYd + 4._dp*(YsIICYsIIYsII) + YsIICYzIITpYzII + YzIIadjYzIIYsII)

 
 
If (TwoLoopRGE) Then 
betaYsII2 = -2._dp*(YdadjYdYdadjYdYsII) + (-6._dp*(AbsL1II) + 2._dp*(g1p2)/5._dp + 6._dp*(g2p2) - & 
&  6._dp*(TrYdadjYd) - 2._dp*(TrYeadjYe))*YdadjYdYsII - 2._dp*(YdadjYuYuadjYdYsII) +     & 
&  (4*(42._dp*(g1p4) + 32*g1p2*g3p2 + 400._dp*(g3p4) - 15._dp*(TrYdadjYdYsIICYsII) -     & 
&  (g1p2 + 5._dp*(g3p2))*TrYsIICYsII - 30._dp*(TrYsIICYsIIYsIICYsII) - 15._dp*(TrYsIICYsIIYzIIadjYzII))*YsII)/15._dp -& 
&  6*AbsL1II*YsIICYdTpYd + (2*g1p2*YsIICYdTpYd)/5._dp + 6*g2p2*YsIICYdTpYd -             & 
&  6*TrYdadjYd*YsIICYdTpYd - 2*TrYeadjYe*YsIICYdTpYd - 2._dp*(YsIICYdTpYdCYdTpYd) -      & 
&  8._dp*(YsIICYdTpYdCYsIIYsII) - 2._dp*(YsIICYdTpYuCYuTpYd) - 8._dp*(YsIICYsIIYdadjYdYsII) +& 
&  (64*g1p2*YsIICYsIIYsII)/15._dp + (160*g3p2*YsIICYsIIYsII)/3._dp - 8*TrYsIICYsII*YsIICYsIIYsII -& 
&  32._dp*(YsIICYsIIYsIICYsIIYsII) - 8._dp*(YsIICYsIIYzIIadjYzIIYsII) - 2._dp*(YsIICYzIITpYeCYeTpYzII) +& 
&  (2*g1p2*YsIICYzIITpYzII)/5._dp + 6*g2p2*YsIICYzIITpYzII - 2*TrYzIIadjYzII*YsIICYzIITpYzII -& 
&  8._dp*(YsIICYzIITpYzIICYsIIYsII) - 6._dp*(YsIICYzIITpYzIICYzIITpYzII) -               & 
&  6._dp*(YsIICYzIIYtIICYtIITpYzII) - 2._dp*(YzIIadjYeYeadjYzIIYsII) + (2*g1p2*YzIIadjYzIIYsII)/5._dp +& 
&  6*g2p2*YzIIadjYzIIYsII - 2*TrYzIIadjYzII*YzIIadjYzIIYsII - 6._dp*(YzIIadjYzIIYzIIadjYzIIYsII) -& 
&  6._dp*(YzIICYtIIYtIIadjYzIIYsII)

 
DYsII = oo16pi2*( betaYsII1 + oo16pi2 * betaYsII2 ) 

 
Else 
DYsII = oo16pi2* betaYsII1 
End If 
 
 
!-------------------- 
! YzII 
!-------------------- 
 
betaYzII1  = 2._dp*(YdadjYdYzII) + 4._dp*(YsIICYsIIYzII) + (-7._dp*(g1p2)             & 
& /15._dp - 3._dp*(g2p2) - 16._dp*(g3p2)/3._dp + TrYzIIadjYzII)*YzII + YzIIadjYeYe +     & 
&  5._dp*(YzIIadjYzIIYzII) + 3._dp*(YzIICYtIIYtII)

 
 
If (TwoLoopRGE) Then 
betaYzII2 = -2._dp*(YdadjYdYdadjYdYzII) + (-6._dp*(AbsL1II) + 2._dp*(g1p2)/5._dp + 6._dp*(g2p2) - & 
&  6._dp*(TrYdadjYd) - 2._dp*(TrYeadjYe))*YdadjYdYzII - 2._dp*(YdadjYuYuadjYdYzII) -     & 
&  8._dp*(YsIICYdTpYdCYsIIYzII) - 16._dp*(YsIICYsIIYsIICYsIIYzII) + (32*g1p2*YsIICYsIIYzII)/15._dp +& 
&  (80*g3p2*YsIICYsIIYzII)/3._dp - 4*TrYsIICYsII*YsIICYsIIYzII - 8._dp*(YsIICYzIITpYzIICYsIIYzII) +& 
&  (581._dp*(g1p4)/90._dp + g1p2*g2p2 + 57._dp*(g2p4)/2._dp + (8*g1p2*g3p2)/9._dp +      & 
&  8*g2p2*g3p2 + 320._dp*(g3p4)/9._dp - 2._dp*(TrYdadjYdYzIIadjYzII) - TrYeadjYzIIYzIIadjYe -& 
&  4._dp*(TrYsIICYsIIYzIIadjYzII) - 3._dp*(TrYtIIadjYzIIYzIICYtII) + (2*g1p2*TrYzIIadjYzII)/5._dp -& 
&  5._dp*(TrYzIIadjYzIIYzIIadjYzII))*YzII - 3*AbsL1II*YzIIadjYeYe + (6*g1p2*YzIIadjYeYe)/5._dp -& 
&  3*TrYdadjYd*YzIIadjYeYe - TrYeadjYe*YzIIadjYeYe - 2._dp*(YzIIadjYeYeadjYeYe) -        & 
&  2._dp*(YzIIadjYeYeadjYzIIYzII) - 6._dp*(YzIIadjYzIIYdadjYdYzII) - 12._dp*(YzIIadjYzIIYsIICYsIIYzII) +& 
&  6*g2p2*YzIIadjYzIIYzII + 16*g3p2*YzIIadjYzIIYzII - 5*TrYzIIadjYzII*YzIIadjYzIIYzII -  & 
&  12._dp*(YzIIadjYzIIYzIIadjYzIIYzII) - 3._dp*(YzIICYtIITpYeCYeYtII) - 9._dp*(YzIICYtIITpYzIICYzIIYtII) -& 
&  3*AbsL1II*YzIICYtIIYtII + (18*g1p2*YzIICYtIIYtII)/5._dp + 12*g2p2*YzIICYtIIYtII -     & 
&  3*TrYtIICYtII*YzIICYtIIYtII - 6._dp*(YzIICYtIIYtIIadjYzIIYzII) - 9._dp*(YzIICYtIIYtIICYtIIYtII)

 
DYzII = oo16pi2*( betaYzII1 + oo16pi2 * betaYzII2 ) 

 
Else 
DYzII = oo16pi2* betaYzII1 
End If 
 
 
!-------------------- 
! L1II 
!-------------------- 
 
betaL1II1  = (-9*g1p2*L1II)/5._dp - 7*g2p2*L1II + 6*L1II*TrYdadjYd + 2*L1II*TrYeadjYe +& 
&  L1II*TrYtIICYtII + 7*L1IIp2*Conjg(L1II)

 
 
If (TwoLoopRGE) Then 
betaL1II2 = -(L1II*(-261._dp*(g1p4) - 114*g1p2*g2p2 - 765._dp*(g2p4) + 300*CL1IIp2*L1IIp2 +       & 
&  8*(g1p2 - 40._dp*(g3p2))*TrYdadjYd + 180._dp*(TrYdadjYdYdadjYd) + 240._dp*(TrYdadjYdYsIICYsII) +& 
&  120._dp*(TrYdadjYdYzIIadjYzII) + 60._dp*(TrYdadjYuYuadjYd) - 24*g1p2*TrYeadjYe +      & 
&  60._dp*(TrYeadjYeYeadjYe) + 60._dp*(TrYeadjYzIIYzIIadjYe) + 80._dp*(TrYeCYtIIYtIIadjYe) +& 
&  60._dp*(TrYtIIadjYzIIYzIICYtII) + 6*g1p2*TrYtIICYtII + 10*g2p2*TrYtIICYtII +          & 
&  2*AbsL1II*(-33._dp*(g1p2) - 115._dp*(g2p2) + 120._dp*(TrYdadjYd) + 40._dp*(TrYeadjYe) +& 
&  30._dp*(TrYtIICYtII)) + 60._dp*(TrYtIICYtIIYtIICYtII)))/10._dp

 
DL1II = oo16pi2*( betaL1II1 + oo16pi2 * betaL1II2 ) 

 
Else 
DL1II = oo16pi2* betaL1II1 
End If 
 
 
!-------------------- 
! L2II 
!-------------------- 
 
betaL2II1  = (-9*g1p2*L2II)/5._dp - 7*g2p2*L2II + 6*L2II*TrYuadjYu + 7*L2IIp2*Conjg(L2II)

 
 
If (TwoLoopRGE) Then 
betaL2II2 = -(L2II*(300*CL2IIp2*L2IIp2 - 2*AbsL2II*(33._dp*(g1p2) + 115._dp*(g2p2) -              & 
&  120._dp*(TrYuadjYu)) - 16*(g1p2 + 20._dp*(g3p2))*TrYuadjYu - 3*(87._dp*(g1p4) +       & 
&  38*g1p2*g2p2 + 255._dp*(g2p4) - 20._dp*(TrYdadjYuYuadjYd) - 60._dp*(TrYuadjYuYuadjYu))))/10._dp

 
DL2II = oo16pi2*( betaL2II1 + oo16pi2 * betaL2II2 ) 

 
Else 
DL2II = oo16pi2* betaL2II1 
End If 
 
 
Call ParametersToG115(Dg1,Dg2,Dg3,DYu,DYd,DYe,DYtII,DYsII,DYzII,DL1II,DL2II,f)

Iname = Iname - 1 
 
End Subroutine rge115  

 Subroutine GToParameters365(g,g1,g2,g3,Yu,Yd,Ye,YtII,YsII,YzII,L1II,L2II,      & 
 & Mu,MTII,MZII,MSII,Tu,Td,Te,TtII,TsII,TzII,TL1II,TL2II,Bmu,BMTII,BMZII,BMSII, &
 & mq2,ml2,mHd2,mHu2,md2,mu2,me2,mt2,mtb2,ms2,msb2,mzz2,mzb2,M1,M2,M3,WOp)

 Implicit None 
  Real(dp), Intent(in) :: g(365) 
  Real(dp),Intent(out) :: g1,g2,g3,mHd2,mHu2,mt2,mtb2,ms2,msb2,mzz2,mzb2

  Complex(dp),Intent(out), Dimension(3,3) :: Yu,Yd,Ye,YtII,YsII,YzII,Tu,Td,Te &
      & ,TtII,TsII,TzII,mq2,ml2,md2,mu2,me2,WOp
  Complex(dp),Intent(out) :: L1II,L2II,Mu,MTII,MZII,MSII,TL1II,TL2II,Bmu,BMTII &
      & ,BMZII,BMSII,M1,M2,M3

  Integer :: i1, i2, SumI 

  Iname = Iname +1 
  NameOfUnit(Iname) = 'GToParameters365' 
 
  g1= g(1) 
  g2= g(2) 
  g3= g(3) 
  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ((i2-1) + (i1-1)*3)
    Yu(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
    Yd(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp) 
    Ye(i1,i2) = Cmplx( g(SumI+40), g(SumI+41), dp) 
    YtII(i1,i2) = Cmplx( g(SumI+58), g(SumI+59), dp) 
    YsII(i1,i2) = Cmplx( g(SumI+76), g(SumI+77), dp) 
    YzII(i1,i2) = Cmplx( g(SumI+94), g(SumI+95), dp) 
    Tu(i1,i2) = Cmplx( g(SumI+124), g(SumI+125), dp) 
    Td(i1,i2) = Cmplx( g(SumI+142), g(SumI+143), dp) 
    Te(i1,i2) = Cmplx( g(SumI+160), g(SumI+161), dp) 
    TtII(i1,i2) = Cmplx( g(SumI+178), g(SumI+179), dp) 
    TsII(i1,i2) = Cmplx( g(SumI+196), g(SumI+197), dp) 
    TzII(i1,i2) = Cmplx( g(SumI+214), g(SumI+215), dp) 
    mq2(i1,i2) = Cmplx( g(SumI+244), g(SumI+245), dp) 
    ml2(i1,i2) = Cmplx( g(SumI+262), g(SumI+263), dp) 
    md2(i1,i2) = Cmplx( g(SumI+282), g(SumI+283), dp) 
    mu2(i1,i2) = Cmplx( g(SumI+300), g(SumI+301), dp) 
    me2(i1,i2) = Cmplx( g(SumI+318), g(SumI+319), dp) 
    WOp(i1,i2) = Cmplx( g(SumI+348), g(SumI+349), dp) 
   End Do 
  End Do 
 
  L1II= Cmplx(g(112),g(113),dp) 
  L2II= Cmplx(g(114),g(115),dp) 
  Mu= Cmplx(g(116),g(117),dp) 
  MTII= Cmplx(g(118),g(119),dp) 
  MZII= Cmplx(g(120),g(121),dp) 
  MSII= Cmplx(g(122),g(123),dp)
 
  TL1II= Cmplx(g(232),g(233),dp) 
  TL2II= Cmplx(g(234),g(235),dp) 
  Bmu= Cmplx(g(236),g(237),dp) 
  BMTII= Cmplx(g(238),g(239),dp) 
  BMZII= Cmplx(g(240),g(241),dp) 
  BMSII= Cmplx(g(242),g(243),dp) 
 
  mHd2= g(280) 
  mHu2= g(281) 
 
  mt2= g(336) 
  mtb2= g(337) 
  ms2= g(338) 
  msb2= g(339) 
  mzz2= g(340) 
  mzb2= g(341) 
  M1= Cmplx(g(342),g(343),dp) 
  M2= Cmplx(g(344),g(345),dp) 
  M3= Cmplx(g(346),g(347),dp) 
 
  Do i1=1,365 
  If (g(i1).ne.g(i1)) Then 
   Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
   Write(*,*) "At position ", i1 
   Call TerminateProgram 
  End if 
  End do 

  Iname = Iname - 1 
 
 End Subroutine GToParameters365

 Subroutine ParametersToG365(g1,g2,g3,Yu,Yd,Ye,YtII,YsII,YzII,L1II,L2II,Mu,MTII &
   & ,MZII,MSII,Tu,Td,Te,TtII,TsII,TzII,TL1II,TL2II,Bmu,BMTII,BMZII,BMSII,mq2   &
   & ,ml2,mHd2,mHu2,md2,mu2,me2,mt2,mtb2,ms2,msb2,mzz2,mzb2,M1,M2,M3,WOp,g)

 Implicit None 
 
  Real(dp), Intent(out) :: g(365) 
  Real(dp), Intent(in) :: g1,g2,g3,mHd2,mHu2,mt2,mtb2,ms2,msb2,mzz2,mzb2

  Complex(dp), Intent(in),Dimension(3,3) :: Yu,Yd,Ye,YtII,YsII,YzII,Tu,Td,Te  & 
     & ,TtII,TsII,TzII,mq2,ml2,md2,mu2,me2,WOp
  Complex(dp), Intent(in) :: L1II,L2II,Mu,MTII,MZII,MSII,TL1II,TL2II,Bmu,BMTII &
     & ,BMZII,BMSII ,M1,M2,M3
  
  Integer i1, i2, SumI 

  Iname = Iname +1 
  NameOfUnit(Iname) = 'ParametersToG365' 

  g(1) = g1  
  g(2) = g2  
  g(3) = g3  

  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2*((i2-1) + (i1-1)*3)
    g(SumI+4) = Real(Yu(i1,i2), dp) 
    g(SumI+5) = Aimag(Yu(i1,i2)) 
    g(SumI+22) = Real(Yd(i1,i2), dp) 
    g(SumI+23) = Aimag(Yd(i1,i2)) 
    g(SumI+40) = Real(Ye(i1,i2), dp) 
    g(SumI+41) = Aimag(Ye(i1,i2)) 
    g(SumI+58) = Real(YtII(i1,i2), dp) 
    g(SumI+59) = Aimag(YtII(i1,i2)) 
    g(SumI+76) = Real(YsII(i1,i2), dp) 
    g(SumI+77) = Aimag(YsII(i1,i2)) 
    g(SumI+94) = Real(YzII(i1,i2), dp) 
    g(SumI+95) = Aimag(YzII(i1,i2)) 
    g(SumI+124) = Real(Tu(i1,i2), dp) 
    g(SumI+125) = Aimag(Tu(i1,i2)) 
    g(SumI+142) = Real(Td(i1,i2), dp) 
    g(SumI+143) = Aimag(Td(i1,i2)) 
    g(SumI+160) = Real(Te(i1,i2), dp) 
    g(SumI+161) = Aimag(Te(i1,i2)) 
    g(SumI+178) = Real(TtII(i1,i2), dp) 
    g(SumI+179) = Aimag(TtII(i1,i2)) 
    g(SumI+196) = Real(TsII(i1,i2), dp) 
    g(SumI+197) = Aimag(TsII(i1,i2)) 
    g(SumI+214) = Real(TzII(i1,i2), dp) 
    g(SumI+215) = Aimag(TzII(i1,i2)) 
    g(SumI+244) = Real(mq2(i1,i2), dp) 
    g(SumI+245) = Aimag(mq2(i1,i2)) 
    g(SumI+262) = Real(ml2(i1,i2), dp) 
    g(SumI+263) = Aimag(ml2(i1,i2)) 
    g(SumI+282) = Real(md2(i1,i2), dp) 
    g(SumI+283) = Aimag(md2(i1,i2)) 
    g(SumI+300) = Real(mu2(i1,i2), dp) 
    g(SumI+301) = Aimag(mu2(i1,i2)) 
    g(SumI+318) = Real(me2(i1,i2), dp) 
    g(SumI+319) = Aimag(me2(i1,i2)) 
    g(SumI+348) = Real(WOp(i1,i2), dp) 
    g(SumI+349) = Aimag(WOp(i1,i2)) 
   End Do 
  End Do 

  g(112) = Real(L1II,dp)  
  g(113) = Aimag(L1II)  
  g(114) = Real(L2II,dp)  
  g(115) = Aimag(L2II)  
  g(116) = Real(Mu,dp)  
  g(117) = Aimag(Mu)  
  g(118) = Real(MTII,dp)  
  g(119) = Aimag(MTII)  
  g(120) = Real(MZII,dp)  
  g(121) = Aimag(MZII)  
  g(122) = Real(MSII,dp)  
  g(123) = Aimag(MSII)  

  g(232) = Real(TL1II,dp)  
  g(233) = Aimag(TL1II)  
  g(234) = Real(TL2II,dp)  
  g(235) = Aimag(TL2II)  
  g(236) = Real(Bmu,dp)  
  g(237) = Aimag(Bmu)  
  g(238) = Real(BMTII,dp)  
  g(239) = Aimag(BMTII)  
  g(240) = Real(BMZII,dp)  
  g(241) = Aimag(BMZII)  
  g(242) = Real(BMSII,dp)  
  g(243) = Aimag(BMSII)

  g(280) = mHd2  
  g(281) = mHu2  

  g(336) = mt2  
  g(337) = mtb2  
  g(338) = ms2  
  g(339) = msb2  
  g(340) = mzz2  
  g(341) = mzb2  
  g(342) = Real(M1,dp)  
  g(343) = Aimag(M1)  
  g(344) = Real(M2,dp)  
  g(345) = Aimag(M2)  
  g(346) = Real(M3,dp)  
  g(347) = Aimag(M3)  

  Iname = Iname - 1 
 
 End Subroutine ParametersToG365

Subroutine rge365(len, T, GY, F) 
Implicit None 
Integer, Intent(in) :: len 
Real(dp), Intent(in) :: T, GY(len) 
Real(dp), Intent(out) :: F(len) 
Integer :: i1,i2
Real(dp) :: q 
Real(dp) :: g1,betag11,betag12,Dg1,g2,betag21,betag22,Dg2,g3,betag31,betag32,         & 
& Dg3,mHd2,betamHd21,betamHd22,DmHd2,mHu2,betamHu21,betamHu22,DmHu2,mt2,betamt21,        & 
& betamt22,Dmt2,mtb2,betamtb21,betamtb22,Dmtb2,ms2,betams21,betams22,Dms2,               & 
& msb2,betamsb21,betamsb22,Dmsb2,mzz2,betamzz21,betamzz22,Dmzz2,mzb2,betamzb21,          & 
& betamzb22,Dmzb2
Complex(dp) :: Yu(3,3),betaYu1(3,3),betaYu2(3,3),DYu(3,3),adjYu(3,3),Yd(3,3)          & 
& ,betaYd1(3,3),betaYd2(3,3),DYd(3,3),adjYd(3,3),Ye(3,3),betaYe1(3,3),betaYe2(3,3)       & 
& ,DYe(3,3),adjYe(3,3),YtII(3,3),betaYtII1(3,3),betaYtII2(3,3),DYtII(3,3),               & 
& adjYtII(3,3),YsII(3,3),betaYsII1(3,3),betaYsII2(3,3),DYsII(3,3),adjYsII(3,3)           & 
& ,YzII(3,3),betaYzII1(3,3),betaYzII2(3,3),DYzII(3,3),adjYzII(3,3),L1II,betaL1II1,       & 
& betaL1II2,DL1II,L2II,betaL2II1,betaL2II2,DL2II,Mu,betaMu1,betaMu2,DMu,MTII,            & 
& betaMTII1,betaMTII2,DMTII,MZII,betaMZII1,betaMZII2,DMZII,MSII,betaMSII1,               & 
& betaMSII2,DMSII,Tu(3,3),betaTu1(3,3),betaTu2(3,3),DTu(3,3),adjTu(3,3),Td(3,3)          & 
& ,betaTd1(3,3),betaTd2(3,3),DTd(3,3),adjTd(3,3),Te(3,3),betaTe1(3,3),betaTe2(3,3)       & 
& ,DTe(3,3),adjTe(3,3),TtII(3,3),betaTtII1(3,3),betaTtII2(3,3),DTtII(3,3),               & 
& adjTtII(3,3),TsII(3,3),betaTsII1(3,3),betaTsII2(3,3),DTsII(3,3),adjTsII(3,3)           & 
& ,TzII(3,3),betaTzII1(3,3),betaTzII2(3,3),DTzII(3,3),adjTzII(3,3),TL1II,betaTL1II1,     & 
& betaTL1II2,DTL1II,TL2II,betaTL2II1,betaTL2II2,DTL2II,Bmu,betaBmu1,betaBmu2,            & 
& DBmu,BMTII,betaBMTII1,betaBMTII2,DBMTII,BMZII,betaBMZII1,betaBMZII2,DBMZII,            & 
& BMSII,betaBMSII1,betaBMSII2,DBMSII,mq2(3,3),betamq21(3,3),betamq22(3,3),               & 
& Dmq2(3,3),ml2(3,3),betaml21(3,3),betaml22(3,3),Dml2(3,3),md2(3,3),betamd21(3,3)        & 
& ,betamd22(3,3),Dmd2(3,3),mu2(3,3),betamu21(3,3),betamu22(3,3),Dmu2(3,3),               & 
& me2(3,3),betame21(3,3),betame22(3,3),Dme2(3,3),M1,betaM11,betaM12,DM1,M2,              & 
& betaM21,betaM22,DM2,M3,betaM31,betaM32,DM3,WOp(3,3),betaWOp1(3,3),betaWOp2(3,3)        & 
& ,DWOp(3,3),adjWOp(3,3)
Real(dp) :: Tr1(3),Tr2(3),Tr3(3) 
Real(dp) :: Tr2U1(3,3) 
Real(dp) :: AbsL1II,AbsL2II,AbsTL1II,AbsTL2II,AbsM1,AbsM2,AbsM3
Complex(dp) :: md2Yd(3,3),md2YzII(3,3),me2Ye(3,3),ml2adjYe(3,3),ml2adjYzII(3,3),mq2adjYd(3,3),       & 
& mq2adjYu(3,3),mu2Yu(3,3),Ydmq2(3,3),YdadjYd(3,3),Yeml2(3,3),YeadjYe(3,3),              & 
& YsIICYsII(3,3),YtIIml2(3,3),YtIICYtII(3,3),Yumq2(3,3),YuadjYu(3,3),YzIIml2(3,3),       & 
& YzIIadjYzII(3,3),adjYdmd2(3,3),adjYdYd(3,3),adjYdYsII(3,3),adjYdYzII(3,3),             & 
& adjYdTd(3,3),adjYdTsII(3,3),adjYdTzII(3,3),adjYeme2(3,3),adjYeYe(3,3),adjYeTe(3,3),    & 
& adjYumu2(3,3),adjYuYu(3,3),adjYuTu(3,3),adjYzIImd2(3,3),adjYzIIYd(3,3),adjYzIIYsII(3,3),& 
& adjYzIIYzII(3,3),adjYzIITd(3,3),adjYzIITsII(3,3),adjYzIITzII(3,3),adjTdTd(3,3),        & 
& adjTeTe(3,3),adjTuTu(3,3),adjTzIITzII(3,3),Cmd2CYsII(3,3),Cml2YtII(3,3),               & 
& CYdTpYd(3,3),CYdTpTd(3,3),CYeWOp(3,3),CYeYtII(3,3),CYeTtII(3,3),CYsIImd2(3,3),         & 
& CYsIIYd(3,3),CYsIIYsII(3,3),CYsIIYzII(3,3),CYsIITd(3,3),CYsIITsII(3,3),CYsIITzII(3,3), & 
& CYtIIWOp(3,3),CYtIIYtII(3,3),CYtIITtII(3,3),CYzIIWOp(3,3),CYzIIYtII(3,3),              & 
& CYzIITtII(3,3),CYzIITpYzII(3,3),CYzIITpTzII(3,3),CTdTpTd(3,3),CTeTpTe(3,3),            & 
& CTsIITsII(3,3),CTtIITtII(3,3),CTuTpTu(3,3),CTzIITpTzII(3,3),TdadjTd(3,3),              & 
& TeadjTe(3,3),TsIICTsII(3,3),TuadjTu(3,3),TzIIadjTzII(3,3),md2YdadjYd(3,3),             & 
& md2YsIICYsII(3,3),md2YzIIadjYzII(3,3),me2YeadjYe(3,3),ml2adjYeYe(3,3),ml2adjYzIIYzII(3,3),& 
& ml2CYtIIYtII(3,3),mq2adjYdYd(3,3),mq2adjYuYu(3,3),mu2YuadjYu(3,3),WOpadjYeYe(3,3),     & 
& WOpadjYzIIYzII(3,3),WOpCYtIIYtII(3,3),Ydmq2adjYd(3,3),YdadjYdmd2(3,3),YdadjYdYd(3,3),  & 
& YdadjYdYsII(3,3),YdadjYdYzII(3,3),YdadjYdTd(3,3),YdadjYdTsII(3,3),YdadjYdTzII(3,3),    & 
& YdadjYuYu(3,3),YdadjYuTu(3,3),Yeml2adjYe(3,3),YeadjYeme2(3,3),YeadjYeYe(3,3),          & 
& YeadjYeTe(3,3),YeadjYzIIYzII(3,3),YeadjYzIITzII(3,3),YeCYtIIYtII(3,3),YeCYtIITtII(3,3),& 
& YsIICmd2CYsII(3,3),YsIICYdTpYd(3,3),YsIICYdTpTd(3,3),YsIICYsIImd2(3,3),YsIICYsIIYd(3,3),& 
& YsIICYsIIYsII(3,3),YsIICYsIIYzII(3,3),YsIICYsIITd(3,3),YsIICYsIITsII(3,3),             & 
& YsIICYsIITzII(3,3),YsIICYzIITpYzII(3,3),YsIICYzIITpTzII(3,3),YtIIadjYeYe(3,3),         & 
& YtIIadjYeTe(3,3),YtIIadjYzIIYzII(3,3),YtIIadjYzIITzII(3,3),YtIICYtIIWOp(3,3),          & 
& YtIICYtIIYtII(3,3),YtIICYtIITtII(3,3),Yumq2adjYu(3,3),YuadjYdYd(3,3),YuadjYdTd(3,3),   & 
& YuadjYumu2(3,3),YuadjYuYu(3,3),YuadjYuTu(3,3),YzIIml2adjYzII(3,3),YzIIadjYeYe(3,3),    & 
& YzIIadjYeTe(3,3),YzIIadjYzIImd2(3,3),YzIIadjYzIIYd(3,3),YzIIadjYzIIYsII(3,3),          & 
& YzIIadjYzIIYzII(3,3),YzIIadjYzIITd(3,3),YzIIadjYzIITsII(3,3),YzIIadjYzIITzII(3,3),     & 
& YzIICYtIIYtII(3,3),YzIICYtIITtII(3,3),adjYdmd2Yd(3,3),adjYdYdmq2(3,3),adjYeme2Ye(3,3), & 
& adjYeYeml2(3,3),adjYumu2Yu(3,3),adjYuYumq2(3,3),adjYzIImd2YzII(3,3),adjYzIIYzIIml2(3,3),& 
& CYtIIYtIIml2(3,3),CYtIICml2YtII(3,3),TdadjYdYd(3,3),TdadjYdYsII(3,3),TdadjYdYzII(3,3), & 
& TdadjYuYu(3,3),TeadjYeYe(3,3),TeadjYzIIYzII(3,3),TeCYtIIYtII(3,3),TsIICYdTpYd(3,3),    & 
& TsIICYsIIYd(3,3),TsIICYsIIYsII(3,3),TsIICYsIIYzII(3,3),TsIICYzIITpYzII(3,3),           & 
& TtIIadjYeYe(3,3),TtIIadjYzIIYzII(3,3),TtIICYtIIYtII(3,3),TuadjYdYd(3,3),               & 
& TuadjYuYu(3,3),TzIIadjYeYe(3,3),TzIIadjYzIIYd(3,3),TzIIadjYzIIYsII(3,3),               & 
& TzIIadjYzIIYzII(3,3),TzIICYtIIYtII(3,3),TpYeCYeWOp(3,3),TpYeCYeYtII(3,3),              & 
& TpYeCYeTtII(3,3),TpYzIICYzIIWOp(3,3),TpYzIICYzIIYtII(3,3),TpYzIICYzIITtII(3,3)

Complex(dp) :: TpTeCYeYtII(3,3),TpTzIICYzIIYtII(3,3)

Complex(dp) :: YdadjYu(3,3),YdadjTd(3,3),YdadjTu(3,3),YeadjYzII(3,3),YeadjTe(3,3),YeadjTzII(3,3),    & 
& YsIICTsII(3,3),YtIIadjYe(3,3),YtIIadjYzII(3,3),YtIIadjTe(3,3),YtIIadjTzII(3,3),        & 
& YtIICTtII(3,3),YuadjYd(3,3),YuadjTd(3,3),YuadjTu(3,3),YzIIadjYe(3,3),YzIIadjTe(3,3),   & 
& YzIIadjTzII(3,3),YzIICYtII(3,3),adjYdCmd2(3,3),adjYdCTsII(3,3),adjYeCme2(3,3),         & 
& adjYuCmu2(3,3),adjYzIICmd2(3,3),adjYzIICTsII(3,3),adjTdYd(3,3),adjTdYzII(3,3),         & 
& adjTdCYsII(3,3),adjTdTzII(3,3),adjTeYe(3,3),adjTuYu(3,3),adjTzIIYd(3,3),               & 
& adjTzIIYzII(3,3),adjTzIICYsII(3,3),adjTzIITd(3,3),Cml2adjYe(3,3),Cml2adjYzII(3,3),     & 
& Cmq2adjYd(3,3),Cmq2adjYu(3,3),CYeTpYzII(3,3),CYeTpTzII(3,3),CYtIICml2(3,3),            & 
& CYtIITpYzII(3,3),CYtIITpTzII(3,3),CYuTpYd(3,3),CYuTpTd(3,3),CTdTpYd(3,3),              & 
& CTeYtII(3,3),CTeTtII(3,3),CTeTpYe(3,3),CTsIIYd(3,3),CTsIIYzII(3,3),CTsIITd(3,3),       & 
& CTsIITzII(3,3),CTtIIYtII(3,3),CTuTpYu(3,3),CTzIIYtII(3,3),CTzIITtII(3,3),              & 
& CTzIITpYzII(3,3),TdadjYd(3,3),TdadjYu(3,3),TdadjTu(3,3),TeadjYe(3,3),TeadjYzII(3,3),   & 
& TeadjTzII(3,3),TeCYtII(3,3),TeCTtII(3,3),TsIICYsII(3,3),TtIIadjYe(3,3),TtIIadjYzII(3,3),& 
& TtIIadjTe(3,3),TtIIadjTzII(3,3),TtIICYtII(3,3),TtIICTtII(3,3),TuadjYd(3,3),            & 
& TuadjYu(3,3),TuadjTd(3,3),TzIIadjYe(3,3),TzIIadjYzII(3,3),TzIIadjTe(3,3),              & 
& TzIICYtII(3,3),TzIICTtII(3,3),TpYdCYsII(3,3),TpYdCTsII(3,3),TpYzIICYsII(3,3),          & 
& TpYzIICTsII(3,3),TpTdCYsII(3,3),TpTdCTsII(3,3),TpTeCTe(3,3),TpTzIICYsII(3,3),          & 
& TpTzIICTsII(3,3),TpTzIICTzII(3,3),md2YdadjYu(3,3),md2YzIIadjYe(3,3),md2CYsIIYsII(3,3), & 
& me2YeadjYzII(3,3),ml2YtIICYtII(3,3),ml2adjYzIIYd(3,3),mq2adjYdYzII(3,3),               & 
& mu2YuadjYd(3,3),Ydmq2adjYu(3,3),YdadjYdCmd2(3,3),YdadjYumu2(3,3),YdadjTdCYsII(3,3),    & 
& YdadjTdTd(3,3),YdadjTdTzII(3,3),YdCmq2adjYd(3,3),Yeml2adjYzII(3,3),YeadjYeCme2(3,3),   & 
& YeadjYzIImd2(3,3),YeadjYzIIYd(3,3),YeadjYzIIYsII(3,3),YeadjYzIITd(3,3),YeadjYzIITsII(3,3),& 
& YeadjTeTe(3,3),YeCml2adjYe(3,3),YeCYtIIWOp(3,3),YsIICYzIIWOp(3,3),YsIICYzIIYtII(3,3),  & 
& YsIICYzIITtII(3,3),YsIICTsIITd(3,3),YsIICTsIITzII(3,3),YtIIml2adjYe(3,3),              & 
& YtIIml2adjYzII(3,3),YtIIadjYeme2(3,3),YtIIadjYzIImd2(3,3),YtIIadjYzIIYd(3,3),          & 
& YtIIadjYzIIYsII(3,3),YtIIadjYzIITd(3,3),YtIIadjYzIITsII(3,3),YtIICYtIITpYzII(3,3),     & 
& YtIICYtIITpTzII(3,3),YtIICTtIITtII(3,3),Yumq2adjYd(3,3),YuadjYdmd2(3,3),               & 
& YuadjYdYsII(3,3),YuadjYdYzII(3,3),YuadjYdTsII(3,3),YuadjYdTzII(3,3),YuadjYuCmu2(3,3),  & 
& YuadjTuTu(3,3),YuCmq2adjYu(3,3),YzIIml2adjYe(3,3),YzIIadjYeme2(3,3),YzIIadjYzIICmd2(3,3),& 
& YzIIadjTzIICYsII(3,3),YzIIadjTzIITd(3,3),YzIIadjTzIITzII(3,3),YzIICml2adjYzII(3,3),    & 
& YzIICYtIIWOp(3,3),YzIICYtIICml2(3,3),adjYdmd2YzII(3,3),adjYdYdadjYd(3,3),              & 
& adjYdYdadjYu(3,3),adjYdYdadjTd(3,3),adjYdYdadjTu(3,3),adjYdYsIICYsII(3,3),             & 
& adjYdYsIICTsII(3,3),adjYdYzIIml2(3,3),adjYdYzIIadjYzII(3,3),adjYdCYsIIYd(3,3),         & 
& adjYdCYsIIYzII(3,3),adjYdTdadjYd(3,3),adjYdTdadjYu(3,3),adjYdTdadjTd(3,3),             & 
& adjYdTdadjTu(3,3),adjYdTsIICYsII(3,3),adjYdTsIICTsII(3,3),adjYdTzIIadjYzII(3,3),       & 
& adjYdTzIIadjTzII(3,3),adjYeYeadjYe(3,3),adjYeYeadjYzII(3,3),adjYeYeadjTe(3,3),         & 
& adjYeYeadjTzII(3,3),adjYeTeadjYe(3,3),adjYeTeadjYzII(3,3),adjYeTeadjTe(3,3)

Complex(dp) :: adjYeTeadjTzII(3,3),adjYeTeCYtII(3,3),adjYeTeCTtII(3,3),adjYuYuadjYd(3,3),             & 
& adjYuYuadjYu(3,3),adjYuYuadjTd(3,3),adjYuYuadjTu(3,3),adjYuTuadjYd(3,3),               & 
& adjYuTuadjYu(3,3),adjYuTuadjTd(3,3),adjYuTuadjTu(3,3),adjYzIImd2Yd(3,3),               & 
& adjYzIIYdmq2(3,3),adjYzIIYdadjYd(3,3),adjYzIIYsIICYsII(3,3),adjYzIIYsIICTsII(3,3),     & 
& adjYzIIYzIIadjYe(3,3),adjYzIIYzIIadjYzII(3,3),adjYzIIYzIIadjTe(3,3),adjYzIIYzIIadjTzII(3,3),& 
& adjYzIIYzIICYtII(3,3),adjYzIICYsIIYd(3,3),adjYzIICYsIIYzII(3,3),adjYzIITdadjYd(3,3),   & 
& adjYzIITdadjTd(3,3),adjYzIITsIICYsII(3,3),adjYzIITsIICTsII(3,3),adjYzIITzIIadjYe(3,3), & 
& adjYzIITzIIadjYzII(3,3),adjYzIITzIIadjTe(3,3),adjYzIITzIIadjTzII(3,3),adjYzIITzIICYtII(3,3),& 
& adjYzIITzIICTtII(3,3),adjTdYdadjYd(3,3),adjTdYdadjYu(3,3),adjTdTdadjYd(3,3),           & 
& adjTdTdadjYu(3,3),adjTdTsIICYsII(3,3),adjTdTzIIadjYzII(3,3),adjTeYeadjYe(3,3),         & 
& adjTeYeadjYzII(3,3),adjTeTeadjYe(3,3),adjTeTeadjYzII(3,3),adjTeTeCYtII(3,3),           & 
& adjTuYuadjYd(3,3),adjTuYuadjYu(3,3),adjTuTuadjYd(3,3),adjTuTuadjYu(3,3),               & 
& adjTzIIYzIIadjYe(3,3),adjTzIIYzIIadjYzII(3,3),adjTzIITdadjYd(3,3),adjTzIITsIICYsII(3,3),& 
& adjTzIITzIIadjYe(3,3),adjTzIITzIIadjYzII(3,3),adjTzIITzIICYtII(3,3),Cmd2CYsIIYd(3,3),  & 
& Cmd2CYsIIYzII(3,3),Cmd2CYzIIYtII(3,3),Cme2CYeYtII(3,3),Cml2YtIIadjYe(3,3),             & 
& Cml2YtIIadjYzII(3,3),Cml2TpYzIICYsII(3,3),Cmq2TpYdCYsII(3,3),CYdTpYdCYsII(3,3),        & 
& CYdTpYdCTsII(3,3),CYdTpTdCTsII(3,3),CYeYtIIml2(3,3),CYeCml2YtII(3,3),CYsIImd2Yd(3,3),  & 
& CYsIImd2YzII(3,3),CYsIIYdmq2(3,3),CYsIIYdadjYd(3,3),CYsIIYsIICYsII(3,3),               & 
& CYsIIYsIICTsII(3,3),CYsIIYzIIml2(3,3),CYsIIYzIIadjYzII(3,3),CYsIITdadjYd(3,3),         & 
& CYsIITdadjTd(3,3),CYsIITsIICYsII(3,3),CYsIITsIICTsII(3,3),CYsIITzIIadjYzII(3,3),       & 
& CYsIITzIIadjTzII(3,3),CYtIIYtIIadjYe(3,3),CYtIIYtIIadjYzII(3,3),CYtIIYtIIadjTe(3,3),   & 
& CYtIIYtIIadjTzII(3,3),CYtIIYtIICYtII(3,3),CYtIITtIIadjYe(3,3),CYtIITtIIadjYzII(3,3),   & 
& CYtIITtIIadjTe(3,3),CYtIITtIIadjTzII(3,3),CYtIITtIICYtII(3,3),CYtIITtIICTtII(3,3),     & 
& CYtIITpTeCTe(3,3),CYtIITpTzIICTzII(3,3),CYzIIYtIIml2(3,3),CYzIICml2YtII(3,3),          & 
& CYzIITpYzIICYsII(3,3),CYzIITpYzIICTsII(3,3),CYzIITpTzIICTsII(3,3),CTdTpYdCYsII(3,3),   & 
& CTdTpTdCYsII(3,3),CTsIIYsIICYsII(3,3),CTsIITdadjYd(3,3),CTsIITsIICYsII(3,3),           & 
& CTsIITzIIadjYzII(3,3),CTtIIYtIIadjYe(3,3),CTtIIYtIIadjYzII(3,3),CTtIITtIIadjYe(3,3),   & 
& CTtIITtIIadjYzII(3,3),CTtIITtIICYtII(3,3),CTzIITpYzIICYsII(3,3),CTzIITpTzIICYsII(3,3), & 
& TdadjYdCTsII(3,3),TdadjTdYd(3,3),TdadjTdYzII(3,3),TeadjYzIIYd(3,3),TeadjYzIIYsII(3,3), & 
& TeadjTeYe(3,3),TsIICYzIIYtII(3,3),TsIICTdTpYd(3,3),TsIICTsIIYd(3,3),TsIICTsIIYzII(3,3),& 
& TsIICTzIITpYzII(3,3),TtIIadjYzIIYd(3,3),TtIIadjYzIIYsII(3,3),TtIICYtIITpYzII(3,3),     & 
& TtIICTtIIYtII(3,3),TuadjYdYsII(3,3),TuadjYdYzII(3,3),TuadjTuYu(3,3),TzIIadjYzIICTsII(3,3),& 
& TzIIadjTzIIYd(3,3),TzIIadjTzIIYzII(3,3),TpYdCmd2CYsII(3,3),TpYdCYdTpYd(3,3),           & 
& TpYdCYdTpTd(3,3),TpYdCYsIImd2(3,3),TpYdCYsIIYd(3,3),TpYdCYsIIYsII(3,3),TpYdCYsIIYzII(3,3),& 
& TpYdCYsIITd(3,3),TpYdCYsIITsII(3,3),TpYdCYsIITzII(3,3),TpYdCYzIIWOp(3,3),              & 
& TpYdCYzIIYtII(3,3),TpYdCYzIITtII(3,3),TpYeCYeTpYzII(3,3),TpYeCYeTpTzII(3,3),           & 
& TpYeCTeTtII(3,3),TpYuCYuTpYd(3,3),TpYuCYuTpTd(3,3),TpYzIICmd2CYsII(3,3)

Complex(dp) :: TpYzIICYsIImd2(3,3),TpYzIICYsIIYd(3,3),TpYzIICYsIIYsII(3,3),TpYzIICYsIIYzII(3,3),      & 
& TpYzIICYsIITd(3,3),TpYzIICYsIITsII(3,3),TpYzIICYsIITzII(3,3),TpYzIICYzIITpYzII(3,3),   & 
& TpYzIICYzIITpTzII(3,3),TpYzIICTzIITtII(3,3),TpTdCYdTpYd(3,3),TpTdCYsIIYd(3,3),         & 
& TpTdCYsIIYsII(3,3),TpTdCYsIIYzII(3,3),TpTdCYzIIYtII(3,3),TpTeCYeTpYzII(3,3),           & 
& TpTeCTeYtII(3,3),TpTuCYuTpYd(3,3),TpTzIICYsIIYd(3,3),TpTzIICYsIIYsII(3,3),             & 
& TpTzIICYsIIYzII(3,3),TpTzIICYzIITpYzII(3,3),TpTzIICTzIIYtII(3,3),md2YdadjYdYd(3,3),    & 
& md2YdadjYdYzII(3,3),md2YsIICYsIIYd(3,3),md2YsIICYsIIYzII(3,3),md2YzIIadjYzIIYd(3,3),   & 
& md2YzIIadjYzIIYzII(3,3),me2YeadjYeYe(3,3),ml2adjYeYeadjYe(3,3),ml2adjYeYeadjYzII(3,3), & 
& ml2adjYzIIYzIIadjYe(3,3),ml2adjYzIIYzIIadjYzII(3,3),ml2CYtIIYtIIadjYe(3,3),            & 
& ml2CYtIIYtIIadjYzII(3,3),mq2adjYdYdadjYd(3,3),mq2adjYdYdadjYu(3,3),mq2adjYuYuadjYd(3,3),& 
& mq2adjYuYuadjYu(3,3),mu2YuadjYuYu(3,3),Ydmq2adjYdYd(3,3),Ydmq2adjYdYzII(3,3),          & 
& YdadjYdmd2Yd(3,3),YdadjYdmd2YzII(3,3),YdadjYdYdmq2(3,3),YdadjYdYdadjYd(3,3),           & 
& YdadjYdYsIICYsII(3,3),YdadjYdYzIIml2(3,3),YdadjYdYzIIadjYzII(3,3),YdadjYdTdadjYd(3,3), & 
& YdadjYdTdadjTd(3,3),YdadjYdTsIICYsII(3,3),YdadjYdTsIICTsII(3,3),YdadjYdTzIIadjYzII(3,3),& 
& YdadjYdTzIIadjTzII(3,3),YdadjYuYuadjYd(3,3),YdadjYuTuadjYd(3,3),YdadjYuTuadjTd(3,3),   & 
& YdadjTdTdadjYd(3,3),YdadjTdTsIICYsII(3,3),YdadjTdTzIIadjYzII(3,3),YdadjTuTuadjYd(3,3), & 
& Yeml2adjYeYe(3,3),YeadjYeme2Ye(3,3),YeadjYeYeml2(3,3),YeadjYeYeadjYe(3,3),             & 
& YeadjYeTeadjYe(3,3),YeadjYeTeadjTe(3,3),YeadjYzIIYzIIadjYe(3,3),YeadjYzIITzIIadjYe(3,3),& 
& YeadjYzIITzIIadjTe(3,3),YeadjTeTeadjYe(3,3),YeadjTzIITzIIadjYe(3,3),YeCYtIIYtIIadjYe(3,3),& 
& YeCYtIITtIIadjYe(3,3),YeCYtIITtIIadjTe(3,3),YeCTtIITtIIadjYe(3,3),YsIICmd2CYsIIYd(3,3),& 
& YsIICmd2CYsIIYzII(3,3),YsIICYdTpYdCYsII(3,3),YsIICYdTpTdCTsII(3,3),YsIICYsIImd2Yd(3,3),& 
& YsIICYsIImd2YzII(3,3),YsIICYsIIYdmq2(3,3),YsIICYsIIYdadjYd(3,3),YsIICYsIIYsIICYsII(3,3),& 
& YsIICYsIIYzIIml2(3,3),YsIICYsIIYzIIadjYzII(3,3),YsIICYsIITdadjYd(3,3),YsIICYsIITdadjTd(3,3),& 
& YsIICYsIITsIICYsII(3,3),YsIICYsIITsIICTsII(3,3),YsIICYsIITzIIadjYzII(3,3),             & 
& YsIICYsIITzIIadjTzII(3,3),YsIICYzIITpYzIICYsII(3,3),YsIICYzIITpTzIICTsII(3,3),         & 
& YsIICTdTpTdCYsII(3,3),YsIICTsIITdadjYd(3,3),YsIICTsIITsIICYsII(3,3),YsIICTsIITzIIadjYzII(3,3),& 
& YsIICTzIITpTzIICYsII(3,3),YsIITdadjYdCTsII(3,3),YsIITzIIadjYzIICTsII(3,3),             & 
& YtIIml2CYtIIYtII(3,3),YtIIadjYeTeCYtII(3,3),YtIIadjYeTeCTtII(3,3),YtIIadjYzIIYzIICYtII(3,3),& 
& YtIIadjYzIITzIICYtII(3,3),YtIIadjYzIITzIICTtII(3,3),YtIIadjTeTeCYtII(3,3),             & 
& YtIIadjTzIITzIICYtII(3,3),YtIICYtIIYtIIml2(3,3),YtIICYtIIYtIICYtII(3,3),               & 
& YtIICYtIICml2YtII(3,3),YtIICYtIITtIICYtII(3,3),YtIICYtIITtIICTtII(3,3),YtIICYtIITpTeCTe(3,3),& 
& YtIICYtIITpTzIICTzII(3,3),YtIICTtIITtIICYtII(3,3),Yumq2adjYuYu(3,3),YuadjYdYdadjYu(3,3),& 
& YuadjYdTdadjYu(3,3),YuadjYdTdadjTu(3,3),YuadjYumu2Yu(3,3),YuadjYuYumq2(3,3),           & 
& YuadjYuYuadjYu(3,3),YuadjYuTuadjYu(3,3),YuadjYuTuadjTu(3,3),YuadjTdTdadjYu(3,3),       & 
& YuadjTuTuadjYu(3,3),YzIIml2adjYzIIYd(3,3),YzIIml2adjYzIIYzII(3,3),YzIIadjYeYeadjYzII(3,3),& 
& YzIIadjYeTeadjYzII(3,3),YzIIadjYeTeadjTzII(3,3),YzIIadjYzIImd2Yd(3,3),YzIIadjYzIImd2YzII(3,3),& 
& YzIIadjYzIIYdmq2(3,3),YzIIadjYzIIYdadjYd(3,3),YzIIadjYzIIYsIICYsII(3,3)

Complex(dp) :: YzIIadjYzIIYzIIml2(3,3),YzIIadjYzIIYzIIadjYzII(3,3),YzIIadjYzIITdadjYd(3,3),           & 
& YzIIadjYzIITdadjTd(3,3),YzIIadjYzIITsIICYsII(3,3),YzIIadjYzIITsIICTsII(3,3),           & 
& YzIIadjYzIITzIIadjYzII(3,3),YzIIadjYzIITzIIadjTzII(3,3),YzIIadjTeTeadjYzII(3,3),       & 
& YzIIadjTzIITdadjYd(3,3),YzIIadjTzIITsIICYsII(3,3),YzIIadjTzIITzIIadjYzII(3,3),         & 
& YzIICYtIIYtIIadjYzII(3,3),YzIICYtIITtIIadjYzII(3,3),YzIICYtIITtIIadjTzII(3,3),         & 
& YzIICTtIITtIIadjYzII(3,3),adjYdmd2YdadjYd(3,3),adjYdmd2YdadjYu(3,3),adjYdYdmq2adjYd(3,3),& 
& adjYdYdmq2adjYu(3,3),adjYdYdadjYdmd2(3,3),adjYdYdadjYdYd(3,3),adjYdYdadjYdYsII(3,3),   & 
& adjYdYdadjYdYzII(3,3),adjYdYdadjYdTd(3,3),adjYdYdadjYdTsII(3,3),adjYdYdadjYdTzII(3,3), & 
& adjYdYdadjYumu2(3,3),adjYdYdadjYuYu(3,3),adjYdYdadjYuTu(3,3),adjYdYdadjTdTd(3,3),      & 
& adjYdYsIICmd2CYsII(3,3),adjYdYsIICYsIIYd(3,3),adjYdYsIICYsIITd(3,3),adjYdYsIICTsIITd(3,3),& 
& adjYdYzIIadjYzIIYd(3,3),adjYdYzIIadjYzIITd(3,3),adjYdYzIIadjTzIITd(3,3),               & 
& adjYdTdadjYdYd(3,3),adjYdTdadjYdYsII(3,3),adjYdTdadjYdYzII(3,3),adjYdTdadjYuYu(3,3),   & 
& adjYdTdadjTdYd(3,3),adjYdTsIICYsIIYd(3,3),adjYdTsIICTsIIYd(3,3),adjYdTzIIadjYzIIYd(3,3),& 
& adjYdTzIIadjTzIIYd(3,3),adjYeme2YeadjYe(3,3),adjYeme2YeadjYzII(3,3),adjYeYeml2adjYe(3,3),& 
& adjYeYeml2adjYzII(3,3),adjYeYeadjYeme2(3,3),adjYeYeadjYeYe(3,3),adjYeYeadjYeTe(3,3),   & 
& adjYeYeadjYzIImd2(3,3),adjYeYeadjYzIIYd(3,3),adjYeYeadjYzIIYsII(3,3),adjYeYeadjYzIIYzII(3,3),& 
& adjYeYeadjYzIITd(3,3),adjYeYeadjYzIITsII(3,3),adjYeYeadjYzIITzII(3,3),adjYeYeadjTeTe(3,3),& 
& adjYeYeCYtIIWOp(3,3),adjYeYeCYtIIYtII(3,3),adjYeYeCYtIITtII(3,3),adjYeTeadjYeYe(3,3),  & 
& adjYeTeadjYzIIYd(3,3),adjYeTeadjYzIIYsII(3,3),adjYeTeadjYzIIYzII(3,3),adjYeTeadjTeYe(3,3),& 
& adjYeTeCYtIIYtII(3,3),adjYumu2YuadjYd(3,3),adjYumu2YuadjYu(3,3),adjYuYumq2adjYd(3,3),  & 
& adjYuYumq2adjYu(3,3),adjYuYuadjYdmd2(3,3),adjYuYuadjYdYd(3,3),adjYuYuadjYdYsII(3,3),   & 
& adjYuYuadjYdYzII(3,3),adjYuYuadjYdTd(3,3),adjYuYuadjYdTsII(3,3),adjYuYuadjYdTzII(3,3), & 
& adjYuYuadjYumu2(3,3),adjYuYuadjYuYu(3,3),adjYuYuadjYuTu(3,3),adjYuYuadjTuTu(3,3),      & 
& adjYuTuadjYdYd(3,3),adjYuTuadjYdYsII(3,3),adjYuTuadjYdYzII(3,3),adjYuTuadjYuYu(3,3),   & 
& adjYuTuadjTuYu(3,3),adjYzIImd2YzIIadjYe(3,3),adjYzIImd2YzIIadjYzII(3,3),               & 
& adjYzIIYdadjYdYzII(3,3),adjYzIIYdadjYdTzII(3,3),adjYzIIYdadjTdTzII(3,3),               & 
& adjYzIIYsIICYsIIYzII(3,3),adjYzIIYsIICYsIITzII(3,3),adjYzIIYsIICTsIITzII(3,3),         & 
& adjYzIIYzIIml2adjYe(3,3),adjYzIIYzIIml2adjYzII(3,3),adjYzIIYzIIadjYeme2(3,3),          & 
& adjYzIIYzIIadjYeYe(3,3),adjYzIIYzIIadjYeTe(3,3),adjYzIIYzIIadjYzIImd2(3,3),            & 
& adjYzIIYzIIadjYzIIYd(3,3),adjYzIIYzIIadjYzIIYsII(3,3),adjYzIIYzIIadjYzIIYzII(3,3),     & 
& adjYzIIYzIIadjYzIITd(3,3),adjYzIIYzIIadjYzIITsII(3,3),adjYzIIYzIIadjYzIITzII(3,3),     & 
& adjYzIIYzIIadjTzIITzII(3,3),adjYzIIYzIICYtIIWOp(3,3),adjYzIIYzIICYtIIYtII(3,3),        & 
& adjYzIIYzIICYtIICml2(3,3),adjYzIIYzIICYtIITtII(3,3),adjYzIITdadjYdYzII(3,3),           & 
& adjYzIITdadjTdYzII(3,3),adjYzIITsIICYsIIYzII(3,3),adjYzIITsIICTsIIYzII(3,3),           & 
& adjYzIITzIIadjYeYe(3,3),adjYzIITzIIadjYzIIYd(3,3),adjYzIITzIIadjYzIIYsII(3,3),         & 
& adjYzIITzIIadjYzIIYzII(3,3),adjYzIITzIIadjTzIIYzII(3,3),adjYzIITzIICYtIIYtII(3,3),     & 
& adjTdYdadjYdTd(3,3),adjTdYsIICYsIITd(3,3),adjTdYzIIadjYzIITd(3,3),adjTdTdadjYdYd(3,3)

Complex(dp) :: adjTdTsIICYsIIYd(3,3),adjTdTzIIadjYzIIYd(3,3),adjTeYeadjYeTe(3,3),adjTeTeadjYeYe(3,3), & 
& adjTuYuadjYuTu(3,3),adjTuTuadjYuYu(3,3),adjTzIIYdadjYdTzII(3,3),adjTzIIYsIICYsIITzII(3,3),& 
& adjTzIIYzIIadjYzIITzII(3,3),adjTzIITdadjYdYzII(3,3),adjTzIITsIICYsIIYzII(3,3),         & 
& adjTzIITzIIadjYzIIYzII(3,3),Cmd2CYdTpYdCYsII(3,3),Cmd2CYsIIYsIICYsII(3,3),             & 
& Cmd2CYsIIYzIIadjYzII(3,3),Cmd2CYzIITpYzIICYsII(3,3),Cml2YtIICYtIIYtII(3,3),            & 
& Cml2TpYeCYeYtII(3,3),Cml2TpYzIICYzIIYtII(3,3),CYdCmq2TpYdCYsII(3,3),CYdTpYdCmd2CYsII(3,3),& 
& CYdTpYdCYdTpYd(3,3),CYdTpYdCYdTpTd(3,3),CYdTpYdCYsIImd2(3,3),CYdTpYdCYsIIYd(3,3),      & 
& CYdTpYdCYsIIYsII(3,3),CYdTpYdCYsIIYzII(3,3),CYdTpYdCYsIITd(3,3),CYdTpYdCYsIITsII(3,3), & 
& CYdTpYdCYsIITzII(3,3),CYdTpYdCYzIIWOp(3,3),CYdTpYdCYzIIYtII(3,3),CYdTpYdCYzIITtII(3,3),& 
& CYdTpYuCYuTpYd(3,3),CYdTpYuCYuTpTd(3,3),CYdTpTdCYdTpYd(3,3),CYdTpTdCYsIIYd(3,3),       & 
& CYdTpTdCYsIIYsII(3,3),CYdTpTdCYsIIYzII(3,3),CYdTpTdCYzIIYtII(3,3),CYdTpTuCYuTpYd(3,3), & 
& CYeTpYeCYeWOp(3,3),CYeTpYeCYeYtII(3,3),CYeTpYeCYeTtII(3,3),CYeTpTeCYeYtII(3,3),        & 
& CYsIImd2YsIICYsII(3,3),CYsIIYdadjYdYsII(3,3),CYsIIYdadjYdTsII(3,3),CYsIIYsIICmd2CYsII(3,3),& 
& CYsIIYsIICYsIImd2(3,3),CYsIIYsIICYsIIYd(3,3),CYsIIYsIICYsIIYsII(3,3),CYsIIYsIICYsIIYzII(3,3),& 
& CYsIIYsIICYsIITd(3,3),CYsIIYsIICYsIITsII(3,3),CYsIIYsIICYsIITzII(3,3),CYsIIYsIICYzIIWOp(3,3),& 
& CYsIIYsIICYzIIYtII(3,3),CYsIIYsIICYzIITtII(3,3),CYsIIYzIIadjYzIIYsII(3,3),             & 
& CYsIIYzIIadjYzIITsII(3,3),CYsIITdadjYdYsII(3,3),CYsIITsIICYsIIYd(3,3),CYsIITsIICYsIIYsII(3,3),& 
& CYsIITsIICYsIIYzII(3,3),CYsIITsIICYzIIYtII(3,3),CYsIITsIICTdTpYd(3,3),CYsIITsIICTzIITpYzII(3,3),& 
& CYsIITzIIadjYzIIYsII(3,3),CYtIIYtIIml2adjYe(3,3),CYtIIYtIIml2adjYzII(3,3),             & 
& CYtIIYtIIadjYeme2(3,3),CYtIIYtIIadjYeYe(3,3),CYtIIYtIIadjYeTe(3,3),CYtIIYtIIadjYzIImd2(3,3),& 
& CYtIIYtIIadjYzIIYd(3,3),CYtIIYtIIadjYzIIYsII(3,3),CYtIIYtIIadjYzIIYzII(3,3),           & 
& CYtIIYtIIadjYzIITd(3,3),CYtIIYtIIadjYzIITsII(3,3),CYtIIYtIIadjYzIITzII(3,3),           & 
& CYtIIYtIICYtIIWOp(3,3),CYtIIYtIICYtIIYtII(3,3),CYtIIYtIICYtIITtII(3,3),CYtIIYtIICTtIITtII(3,3),& 
& CYtIICml2YtIIadjYe(3,3),CYtIICml2YtIIadjYzII(3,3),CYtIITtIIadjYeYe(3,3),               & 
& CYtIITtIIadjYzIIYd(3,3),CYtIITtIIadjYzIIYsII(3,3),CYtIITtIIadjYzIIYzII(3,3),           & 
& CYtIITtIICYtIIYtII(3,3),CYtIITtIICTtIIYtII(3,3),CYtIITpYeCYeYtII(3,3),CYtIITpYeCYeTtII(3,3),& 
& CYtIITpYeCTeTtII(3,3),CYtIITpYzIICYzIIYtII(3,3),CYtIITpYzIICYzIITtII(3,3),             & 
& CYtIITpYzIICTzIITtII(3,3),CYtIITpTeCYeYtII(3,3),CYtIITpTeCTeYtII(3,3),CYtIITpTzIICYzIIYtII(3,3),& 
& CYtIITpTzIICTzIIYtII(3,3),CYzIIYtIICYtIITpYzII(3,3),CYzIIYtIICYtIITpTzII(3,3),         & 
& CYzIICml2TpYzIICYsII(3,3),CYzIITtIICYtIITpYzII(3,3),CYzIITpYeCYeTpYzII(3,3),           & 
& CYzIITpYeCYeTpTzII(3,3),CYzIITpYzIICmd2CYsII(3,3),CYzIITpYzIICYsIImd2(3,3),            & 
& CYzIITpYzIICYsIIYd(3,3),CYzIITpYzIICYsIIYsII(3,3),CYzIITpYzIICYsIIYzII(3,3),           & 
& CYzIITpYzIICYsIITd(3,3),CYzIITpYzIICYsIITsII(3,3),CYzIITpYzIICYsIITzII(3,3),           & 
& CYzIITpYzIICYzIIWOp(3,3),CYzIITpYzIICYzIIYtII(3,3),CYzIITpYzIICYzIITtII(3,3),          & 
& CYzIITpYzIICYzIITpYzII(3,3),CYzIITpYzIICYzIITpTzII(3,3),CYzIITpTeCYeTpYzII(3,3),       & 
& CYzIITpTzIICYsIIYd(3,3),CYzIITpTzIICYsIIYsII(3,3),CYzIITpTzIICYsIIYzII(3,3),           & 
& CYzIITpTzIICYzIIYtII(3,3),CYzIITpTzIICYzIITpYzII(3,3),CTtIIYtIICYtIITtII(3,3)

Complex(dp) :: CTtIITtIICYtIIYtII(3,3),CTtIITpYeCYeTtII(3,3),CTtIITpYzIICYzIITtII(3,3),               & 
& CTtIITpTeCYeYtII(3,3),CTtIITpTzIICYzIIYtII(3,3),TdadjYdYdadjTd(3,3),TdadjYdYsIICTsII(3,3),& 
& TdadjYdCYsIIYd(3,3),TdadjYdCYsIIYzII(3,3),TdadjYuYuadjTd(3,3),TdadjTdYdadjYd(3,3),     & 
& TdadjTuYuadjYd(3,3),TeadjYeYeadjTe(3,3),TeadjYzIIYzIIadjTe(3,3),TeadjTeYeadjYe(3,3),   & 
& TeadjTzIIYzIIadjYe(3,3),TeCYtIIYtIIadjTe(3,3),TeCTtIIYtIIadjYe(3,3),TsIIYdadjTdCYsII(3,3),& 
& TsIIYzIIadjTzIICYsII(3,3),TsIICYdTpYdCTsII(3,3),TsIICYsIIYsIICTsII(3,3),               & 
& TsIICYzIITpYzIICTsII(3,3),TsIICTdTpYdCYsII(3,3),TsIICTsIIYsIICYsII(3,3),               & 
& TsIICTzIITpYzIICYsII(3,3),TuadjYdYdadjTu(3,3),TuadjYuYuadjTu(3,3),TuadjTdYdadjYu(3,3), & 
& TuadjTuYuadjYu(3,3),TzIIadjYeYeadjTzII(3,3),TzIIadjYzIIYsIICTsII(3,3),TzIIadjYzIIYzIIadjTzII(3,3),& 
& TzIIadjYzIICYsIIYd(3,3),TzIIadjYzIICYsIIYzII(3,3),TzIIadjTeYeadjYzII(3,3),             & 
& TzIIadjTzIIYzIIadjYzII(3,3),TzIICYtIIYtIIadjTzII(3,3),TzIICTtIIYtIIadjYzII(3,3),       & 
& TpYeCme2CYeYtII(3,3),TpYeCYeYtIIml2(3,3),TpYeCYeCml2YtII(3,3),TpYzIICmd2CYzIIYtII(3,3),& 
& TpYzIICYzIIYtIIml2(3,3),TpYzIICYzIICml2YtII(3,3),md2YdadjYdYdadjYd(3,3),               & 
& md2YdadjYdYsIICYsII(3,3),md2YdadjYdYzIIadjYzII(3,3),md2YdadjYuYuadjYd(3,3),            & 
& md2YsIICYdTpYdCYsII(3,3),md2YsIICYsIIYdadjYd(3,3),md2YsIICYsIIYsIICYsII(3,3),          & 
& md2YsIICYsIIYzIIadjYzII(3,3),md2YsIICYzIITpYzIICYsII(3,3),md2YzIIadjYeYeadjYzII(3,3),  & 
& md2YzIIadjYzIIYdadjYd(3,3),md2YzIIadjYzIIYsIICYsII(3,3),md2YzIIadjYzIIYzIIadjYzII(3,3),& 
& md2YzIICYtIIYtIIadjYzII(3,3),me2YeadjYeYeadjYe(3,3),me2YeadjYzIIYzIIadjYe(3,3),        & 
& me2YeCYtIIYtIIadjYe(3,3),ml2adjYeYeadjYeYe(3,3),ml2adjYeYeadjYzIIYzII(3,3),            & 
& ml2adjYeYeCYtIIYtII(3,3),ml2adjYzIIYdadjYdYzII(3,3),ml2adjYzIIYsIICYsIIYzII(3,3),      & 
& ml2adjYzIIYzIIadjYeYe(3,3),ml2adjYzIIYzIIadjYzIIYzII(3,3),ml2adjYzIIYzIICYtIIYtII(3,3),& 
& ml2CYtIIYtIIadjYeYe(3,3),ml2CYtIIYtIIadjYzIIYzII(3,3),ml2CYtIIYtIICYtIIYtII(3,3),      & 
& ml2CYtIITpYeCYeYtII(3,3),ml2CYtIITpYzIICYzIIYtII(3,3),mq2adjYdYdadjYdYd(3,3),          & 
& mq2adjYdYdadjYuYu(3,3),mq2adjYdYsIICYsIIYd(3,3),mq2adjYdYzIIadjYzIIYd(3,3),            & 
& mq2adjYuYuadjYdYd(3,3),mq2adjYuYuadjYuYu(3,3),mu2YuadjYdYdadjYu(3,3),mu2YuadjYuYuadjYu(3,3),& 
& WOpadjYeYeadjYeYe(3,3),WOpadjYzIIYdadjYdYzII(3,3),WOpadjYzIIYsIICYsIIYzII(3,3),        & 
& WOpadjYzIIYzIIadjYzIIYzII(3,3),WOpCYtIIYtIICYtIIYtII(3,3),WOpCYtIITpYeCYeYtII(3,3),    & 
& WOpCYtIITpYzIICYzIIYtII(3,3),Ydmq2adjYdYdadjYd(3,3),Ydmq2adjYuYuadjYd(3,3),            & 
& YdadjYdmd2YdadjYd(3,3),YdadjYdYdmq2adjYd(3,3),YdadjYdYdadjYdmd2(3,3),YdadjYdYdadjYdYd(3,3),& 
& YdadjYdYdadjYdYsII(3,3),YdadjYdYdadjYdYzII(3,3),YdadjYdYdadjYdTd(3,3),YdadjYdYdadjYdTsII(3,3),& 
& YdadjYdYdadjYdTzII(3,3),YdadjYdYsIICmd2CYsII(3,3),YdadjYdYsIICYsIIYd(3,3),             & 
& YdadjYdYsIICYsIITd(3,3),YdadjYdYzIIadjYzIIYd(3,3),YdadjYdYzIIadjYzIITd(3,3),           & 
& YdadjYdTdadjYdYd(3,3),YdadjYdTdadjYdYsII(3,3),YdadjYdTdadjYdYzII(3,3),YdadjYdTsIICYsIIYd(3,3),& 
& YdadjYdTzIIadjYzIIYd(3,3),YdadjYumu2YuadjYd(3,3),YdadjYuYumq2adjYd(3,3),               & 
& YdadjYuYuadjYdmd2(3,3),YdadjYuYuadjYdYd(3,3),YdadjYuYuadjYdYsII(3,3),YdadjYuYuadjYdYzII(3,3),& 
& YdadjYuYuadjYdTd(3,3),YdadjYuYuadjYdTsII(3,3),YdadjYuYuadjYdTzII(3,3),YdadjYuYuadjYuYu(3,3),& 
& YdadjYuYuadjYuTu(3,3),YdadjYuTuadjYdYd(3,3),YdadjYuTuadjYdYsII(3,3),YdadjYuTuadjYdYzII(3,3)

Complex(dp) :: YdadjYuTuadjYuYu(3,3),Yeml2adjYeYeadjYe(3,3),Yeml2adjYzIIYzIIadjYe(3,3),               & 
& Yeml2CYtIIYtIIadjYe(3,3),YeadjYeme2YeadjYe(3,3),YeadjYeYeml2adjYe(3,3),YeadjYeYeadjYeme2(3,3),& 
& YeadjYeYeadjYeYe(3,3),YeadjYeYeadjYeTe(3,3),YeadjYeTeadjYeYe(3,3),YeadjYzIImd2YzIIadjYe(3,3),& 
& YeadjYzIIYdadjYdYzII(3,3),YeadjYzIIYdadjYdTzII(3,3),YeadjYzIIYsIICYsIIYzII(3,3),       & 
& YeadjYzIIYsIICYsIITzII(3,3),YeadjYzIIYzIIml2adjYe(3,3),YeadjYzIIYzIIadjYeme2(3,3),     & 
& YeadjYzIIYzIIadjYeYe(3,3),YeadjYzIIYzIIadjYeTe(3,3),YeadjYzIIYzIIadjYzIIYzII(3,3),     & 
& YeadjYzIIYzIIadjYzIITzII(3,3),YeadjYzIITdadjYdYzII(3,3),YeadjYzIITsIICYsIIYzII(3,3),   & 
& YeadjYzIITzIIadjYeYe(3,3),YeadjYzIITzIIadjYzIIYzII(3,3),YeCYtIIYtIIml2adjYe(3,3),      & 
& YeCYtIIYtIIadjYeme2(3,3),YeCYtIIYtIIadjYeYe(3,3),YeCYtIIYtIIadjYeTe(3,3),              & 
& YeCYtIIYtIICYtIIYtII(3,3),YeCYtIIYtIICYtIITtII(3,3),YeCYtIICml2YtIIadjYe(3,3),         & 
& YeCYtIITtIIadjYeYe(3,3),YeCYtIITtIICYtIIYtII(3,3),YeCYtIITpYeCYeYtII(3,3),             & 
& YeCYtIITpYeCYeTtII(3,3),YeCYtIITpYzIICYzIIYtII(3,3),YeCYtIITpYzIICYzIITtII(3,3),       & 
& YeCYtIITpTeCYeYtII(3,3),YeCYtIITpTzIICYzIIYtII(3,3),YsIICmd2CYdTpYdCYsII(3,3),         & 
& YsIICmd2CYsIIYsIICYsII(3,3),YsIICmd2CYsIIYzIIadjYzII(3,3),YsIICmd2CYzIITpYzIICYsII(3,3),& 
& YsIICYdCmq2TpYdCYsII(3,3),YsIICYdTpYdCmd2CYsII(3,3),YsIICYdTpYdCYdTpYd(3,3),           & 
& YsIICYdTpYdCYdTpTd(3,3),YsIICYdTpYdCYsIImd2(3,3),YsIICYdTpYdCYsIIYd(3,3),              & 
& YsIICYdTpYdCYsIIYsII(3,3),YsIICYdTpYdCYsIIYzII(3,3),YsIICYdTpYdCYsIITd(3,3),           & 
& YsIICYdTpYdCYsIITsII(3,3),YsIICYdTpYdCYsIITzII(3,3),YsIICYdTpYuCYuTpYd(3,3),           & 
& YsIICYdTpYuCYuTpTd(3,3),YsIICYdTpTdCYdTpYd(3,3),YsIICYdTpTdCYsIIYd(3,3),               & 
& YsIICYdTpTdCYsIIYsII(3,3),YsIICYdTpTdCYsIIYzII(3,3),YsIICYdTpTuCYuTpYd(3,3),           & 
& YsIICYsIImd2YsIICYsII(3,3),YsIICYsIIYdadjYdYsII(3,3),YsIICYsIIYdadjYdTsII(3,3),        & 
& YsIICYsIIYsIICmd2CYsII(3,3),YsIICYsIIYsIICYsIImd2(3,3),YsIICYsIIYsIICYsIIYd(3,3),      & 
& YsIICYsIIYsIICYsIIYsII(3,3),YsIICYsIIYsIICYsIIYzII(3,3),YsIICYsIIYsIICYsIITd(3,3),     & 
& YsIICYsIIYsIICYsIITsII(3,3),YsIICYsIIYsIICYsIITzII(3,3),YsIICYsIIYzIIadjYzIIYsII(3,3), & 
& YsIICYsIIYzIIadjYzIITsII(3,3),YsIICYsIITdadjYdYsII(3,3),YsIICYsIITsIICYsIIYd(3,3),     & 
& YsIICYsIITsIICYsIIYsII(3,3),YsIICYsIITsIICYsIIYzII(3,3),YsIICYsIITzIIadjYzIIYsII(3,3), & 
& YsIICYzIIYtIICYtIITpYzII(3,3),YsIICYzIIYtIICYtIITpTzII(3,3),YsIICYzIICml2TpYzIICYsII(3,3),& 
& YsIICYzIITtIICYtIITpYzII(3,3),YsIICYzIITpYeCYeTpYzII(3,3),YsIICYzIITpYeCYeTpTzII(3,3), & 
& YsIICYzIITpYzIICmd2CYsII(3,3),YsIICYzIITpYzIICYsIImd2(3,3),YsIICYzIITpYzIICYsIIYd(3,3),& 
& YsIICYzIITpYzIICYsIIYsII(3,3),YsIICYzIITpYzIICYsIIYzII(3,3),YsIICYzIITpYzIICYsIITd(3,3),& 
& YsIICYzIITpYzIICYsIITsII(3,3),YsIICYzIITpYzIICYsIITzII(3,3),YsIICYzIITpYzIICYzIITpYzII(3,3),& 
& YsIICYzIITpYzIICYzIITpTzII(3,3),YsIICYzIITpTeCYeTpYzII(3,3),YsIICYzIITpTzIICYsIIYd(3,3),& 
& YsIICYzIITpTzIICYsIIYsII(3,3),YsIICYzIITpTzIICYsIIYzII(3,3),YsIICYzIITpTzIICYzIITpYzII(3,3),& 
& YsIITdadjYdCYsIIYd(3,3),YsIITdadjYdCYsIIYzII(3,3),YsIITzIIadjYzIICYsIIYd(3,3),         & 
& YsIITzIIadjYzIICYsIIYzII(3,3),YtIIadjYeYeadjYeYe(3,3),YtIIadjYeYeadjYeTe(3,3),         & 
& YtIIadjYeYeCYtIIWOp(3,3),YtIIadjYeYeCYtIIYtII(3,3),YtIIadjYeYeCYtIITtII(3,3),          & 
& YtIIadjYeTeadjYeYe(3,3),YtIIadjYeTeCYtIIYtII(3,3),YtIIadjYzIIYdadjYdYzII(3,3)

Complex(dp) :: YtIIadjYzIIYdadjYdTzII(3,3),YtIIadjYzIIYsIICYsIIYzII(3,3),YtIIadjYzIIYsIICYsIITzII(3,3),& 
& YtIIadjYzIIYzIIadjYzIIYzII(3,3),YtIIadjYzIIYzIIadjYzIITzII(3,3),YtIIadjYzIIYzIICYtIIWOp(3,3),& 
& YtIIadjYzIIYzIICYtIIYtII(3,3),YtIIadjYzIIYzIICYtIICml2(3,3),YtIIadjYzIIYzIICYtIITtII(3,3),& 
& YtIIadjYzIITdadjYdYzII(3,3),YtIIadjYzIITsIICYsIIYzII(3,3),YtIIadjYzIITzIIadjYzIIYzII(3,3),& 
& YtIIadjYzIITzIICYtIIYtII(3,3),YtIICYtIIYtIICYtIIWOp(3,3),YtIICYtIIYtIICYtIIYtII(3,3),  & 
& YtIICYtIIYtIICYtIITtII(3,3),YtIICYtIITtIICYtIIYtII(3,3),YtIICYtIITpYeCYeYtII(3,3),     & 
& YtIICYtIITpYeCYeTtII(3,3),YtIICYtIITpYzIICYzIIYtII(3,3),YtIICYtIITpYzIICYzIITtII(3,3), & 
& YtIICYtIITpTeCYeYtII(3,3),YtIICYtIITpTzIICYzIIYtII(3,3),Yumq2adjYdYdadjYu(3,3),        & 
& Yumq2adjYuYuadjYu(3,3),YuadjYdmd2YdadjYu(3,3),YuadjYdYdmq2adjYu(3,3),YuadjYdYdadjYdYd(3,3),& 
& YuadjYdYdadjYdTd(3,3),YuadjYdYdadjYumu2(3,3),YuadjYdYdadjYuYu(3,3),YuadjYdYdadjYuTu(3,3),& 
& YuadjYdYsIICYsIIYd(3,3),YuadjYdYsIICYsIITd(3,3),YuadjYdYzIIadjYzIIYd(3,3),             & 
& YuadjYdYzIIadjYzIITd(3,3),YuadjYdTdadjYdYd(3,3),YuadjYdTdadjYuYu(3,3),YuadjYdTsIICYsIIYd(3,3),& 
& YuadjYdTzIIadjYzIIYd(3,3),YuadjYumu2YuadjYu(3,3),YuadjYuYumq2adjYu(3,3),               & 
& YuadjYuYuadjYumu2(3,3),YuadjYuYuadjYuYu(3,3),YuadjYuYuadjYuTu(3,3),YuadjYuTuadjYuYu(3,3),& 
& YzIIml2adjYeYeadjYzII(3,3),YzIIml2adjYzIIYzIIadjYzII(3,3),YzIIml2CYtIIYtIIadjYzII(3,3),& 
& YzIIadjYeme2YeadjYzII(3,3),YzIIadjYeYeml2adjYzII(3,3),YzIIadjYeYeadjYeYe(3,3),         & 
& YzIIadjYeYeadjYeTe(3,3),YzIIadjYeYeadjYzIImd2(3,3),YzIIadjYeYeadjYzIIYd(3,3),          & 
& YzIIadjYeYeadjYzIIYsII(3,3),YzIIadjYeYeadjYzIIYzII(3,3),YzIIadjYeYeadjYzIITd(3,3),     & 
& YzIIadjYeYeadjYzIITsII(3,3),YzIIadjYeYeadjYzIITzII(3,3),YzIIadjYeTeadjYeYe(3,3),       & 
& YzIIadjYeTeadjYzIIYd(3,3),YzIIadjYeTeadjYzIIYsII(3,3),YzIIadjYeTeadjYzIIYzII(3,3),     & 
& YzIIadjYzIImd2YzIIadjYzII(3,3),YzIIadjYzIIYdadjYdYzII(3,3),YzIIadjYzIIYdadjYdTzII(3,3),& 
& YzIIadjYzIIYsIICYsIIYzII(3,3),YzIIadjYzIIYsIICYsIITzII(3,3),YzIIadjYzIIYzIIml2adjYzII(3,3),& 
& YzIIadjYzIIYzIIadjYzIImd2(3,3),YzIIadjYzIIYzIIadjYzIIYd(3,3),YzIIadjYzIIYzIIadjYzIIYsII(3,3),& 
& YzIIadjYzIIYzIIadjYzIIYzII(3,3),YzIIadjYzIIYzIIadjYzIITd(3,3),YzIIadjYzIIYzIIadjYzIITsII(3,3),& 
& YzIIadjYzIIYzIIadjYzIITzII(3,3),YzIIadjYzIITdadjYdYzII(3,3),YzIIadjYzIITsIICYsIIYzII(3,3),& 
& YzIIadjYzIITzIIadjYzIIYd(3,3),YzIIadjYzIITzIIadjYzIIYsII(3,3),YzIIadjYzIITzIIadjYzIIYzII(3,3),& 
& YzIICYtIIYtIIml2adjYzII(3,3),YzIICYtIIYtIIadjYzIImd2(3,3),YzIICYtIIYtIIadjYzIIYd(3,3), & 
& YzIICYtIIYtIIadjYzIIYsII(3,3),YzIICYtIIYtIIadjYzIIYzII(3,3),YzIICYtIIYtIIadjYzIITd(3,3),& 
& YzIICYtIIYtIIadjYzIITsII(3,3),YzIICYtIIYtIIadjYzIITzII(3,3),YzIICYtIIYtIICYtIIYtII(3,3),& 
& YzIICYtIIYtIICYtIITtII(3,3),YzIICYtIICml2YtIIadjYzII(3,3),YzIICYtIITtIIadjYzIIYd(3,3), & 
& YzIICYtIITtIIadjYzIIYsII(3,3),YzIICYtIITtIIadjYzIIYzII(3,3),YzIICYtIITtIICYtIIYtII(3,3),& 
& YzIICYtIITpYeCYeYtII(3,3),YzIICYtIITpYeCYeTtII(3,3),YzIICYtIITpYzIICYzIIYtII(3,3),     & 
& YzIICYtIITpYzIICYzIITtII(3,3),YzIICYtIITpTeCYeYtII(3,3),YzIICYtIITpTzIICYzIIYtII(3,3), & 
& adjYdmd2YdadjYdYd(3,3),adjYdmd2YsIICYsIIYd(3,3),adjYdmd2YzIIadjYzIIYd(3,3),            & 
& adjYdYdmq2adjYdYd(3,3),adjYdYdadjYdmd2Yd(3,3),adjYdYdadjYdYdmq2(3,3),adjYdYsIICmd2CYsIIYd(3,3),& 
& adjYdYsIICYsIImd2Yd(3,3),adjYdYsIICYsIIYdmq2(3,3),adjYdYzIIml2adjYzIIYd(3,3),          & 
& adjYdYzIIadjYzIImd2Yd(3,3),adjYdYzIIadjYzIIYdmq2(3,3),adjYeme2YeadjYeYe(3,3)

Complex(dp) :: adjYeYeml2adjYeYe(3,3),adjYeYeadjYeme2Ye(3,3),adjYeYeadjYeYeml2(3,3),adjYumu2YuadjYuYu(3,3),& 
& adjYuYumq2adjYuYu(3,3),adjYuYuadjYumu2Yu(3,3),adjYuYuadjYuYumq2(3,3),adjYzIImd2YdadjYdYzII(3,3),& 
& adjYzIImd2YsIICYsIIYzII(3,3),adjYzIImd2YzIIadjYzIIYzII(3,3),adjYzIIYdmq2adjYdYzII(3,3),& 
& adjYzIIYdadjYdmd2YzII(3,3),adjYzIIYdadjYdYzIIml2(3,3),adjYzIIYsIICmd2CYsIIYzII(3,3),   & 
& adjYzIIYsIICYsIImd2YzII(3,3),adjYzIIYsIICYsIIYzIIml2(3,3),adjYzIIYzIIml2adjYzIIYzII(3,3),& 
& adjYzIIYzIIadjYzIImd2YzII(3,3),adjYzIIYzIIadjYzIIYzIIml2(3,3),CYtIIYtIIml2CYtIIYtII(3,3),& 
& CYtIIYtIICYtIIYtIIml2(3,3),CYtIIYtIICYtIICml2YtII(3,3),CYtIICml2YtIICYtIIYtII(3,3),    & 
& CYtIICml2TpYeCYeYtII(3,3),CYtIICml2TpYzIICYzIIYtII(3,3),CYtIITpYeCme2CYeYtII(3,3),     & 
& CYtIITpYeCYeYtIIml2(3,3),CYtIITpYeCYeCml2YtII(3,3),CYtIITpYzIICmd2CYzIIYtII(3,3),      & 
& CYtIITpYzIICYzIIYtIIml2(3,3),CYtIITpYzIICYzIICml2YtII(3,3),TdadjYdYdadjYdYd(3,3),      & 
& TdadjYdYdadjYdYsII(3,3),TdadjYdYdadjYdYzII(3,3),TdadjYdYsIICYsIIYd(3,3),               & 
& TdadjYdYzIIadjYzIIYd(3,3),TdadjYuYuadjYdYd(3,3),TdadjYuYuadjYdYsII(3,3),               & 
& TdadjYuYuadjYdYzII(3,3),TdadjYuYuadjYuYu(3,3),TeadjYeYeadjYeYe(3,3),TeadjYzIIYdadjYdYzII(3,3),& 
& TeadjYzIIYsIICYsIIYzII(3,3),TeadjYzIIYzIIadjYeYe(3,3),TeadjYzIIYzIIadjYzIIYzII(3,3),   & 
& TeCYtIIYtIIadjYeYe(3,3),TeCYtIIYtIICYtIIYtII(3,3),TeCYtIITpYeCYeYtII(3,3),             & 
& TeCYtIITpYzIICYzIIYtII(3,3),TsIICYdTpYdCYdTpYd(3,3),TsIICYdTpYdCYsIIYd(3,3),           & 
& TsIICYdTpYdCYsIIYsII(3,3),TsIICYdTpYdCYsIIYzII(3,3),TsIICYdTpYuCYuTpYd(3,3),           & 
& TsIICYsIIYdadjYdYsII(3,3),TsIICYsIIYsIICYsIIYd(3,3),TsIICYsIIYsIICYsIIYsII(3,3),       & 
& TsIICYsIIYsIICYsIIYzII(3,3),TsIICYsIIYzIIadjYzIIYsII(3,3),TsIICYzIIYtIICYtIITpYzII(3,3),& 
& TsIICYzIITpYeCYeTpYzII(3,3),TsIICYzIITpYzIICYsIIYd(3,3),TsIICYzIITpYzIICYsIIYsII(3,3), & 
& TsIICYzIITpYzIICYsIIYzII(3,3),TsIICYzIITpYzIICYzIITpYzII(3,3),TtIIadjYeYeadjYeYe(3,3), & 
& TtIIadjYeYeCYtIIYtII(3,3),TtIIadjYzIIYdadjYdYzII(3,3),TtIIadjYzIIYsIICYsIIYzII(3,3),   & 
& TtIIadjYzIIYzIIadjYzIIYzII(3,3),TtIIadjYzIIYzIICYtIIYtII(3,3),TtIICYtIIYtIICYtIIYtII(3,3),& 
& TtIICYtIITpYeCYeYtII(3,3),TtIICYtIITpYzIICYzIIYtII(3,3),TuadjYdYdadjYdYd(3,3),         & 
& TuadjYdYdadjYuYu(3,3),TuadjYdYsIICYsIIYd(3,3),TuadjYdYzIIadjYzIIYd(3,3),               & 
& TuadjYuYuadjYuYu(3,3),TzIIadjYeYeadjYeYe(3,3),TzIIadjYeYeadjYzIIYd(3,3),               & 
& TzIIadjYeYeadjYzIIYsII(3,3),TzIIadjYeYeadjYzIIYzII(3,3),TzIIadjYzIIYdadjYdYzII(3,3),   & 
& TzIIadjYzIIYsIICYsIIYzII(3,3),TzIIadjYzIIYzIIadjYzIIYd(3,3),TzIIadjYzIIYzIIadjYzIIYsII(3,3),& 
& TzIIadjYzIIYzIIadjYzIIYzII(3,3),TzIICYtIIYtIIadjYzIIYd(3,3),TzIICYtIIYtIIadjYzIIYsII(3,3),& 
& TzIICYtIIYtIIadjYzIIYzII(3,3),TzIICYtIIYtIICYtIIYtII(3,3),TzIICYtIITpYeCYeYtII(3,3),   & 
& TzIICYtIITpYzIICYzIIYtII(3,3),TpYeCYeTpYeCYeWOp(3,3),TpYeCYeTpYeCYeYtII(3,3),          & 
& TpYeCYeTpYeCYeTtII(3,3),TpYeCYeTpTeCYeYtII(3,3),TpYzIICYdTpYdCYzIIWOp(3,3),            & 
& TpYzIICYdTpYdCYzIIYtII(3,3),TpYzIICYdTpYdCYzIITtII(3,3),TpYzIICYdTpTdCYzIIYtII(3,3),   & 
& TpYzIICYsIIYsIICYzIIWOp(3,3),TpYzIICYsIIYsIICYzIIYtII(3,3),TpYzIICYsIIYsIICYzIITtII(3,3),& 
& TpYzIICYsIITsIICYzIIYtII(3,3),TpYzIICYzIITpYzIICYzIIWOp(3,3),TpYzIICYzIITpYzIICYzIIYtII(3,3),& 
& TpYzIICYzIITpYzIICYzIITtII(3,3),TpYzIICYzIITpTzIICYzIIYtII(3,3),TpTeCYeTpYeCYeYtII(3,3),& 
& TpTzIICYdTpYdCYzIIYtII(3,3),TpTzIICYsIIYsIICYzIIYtII(3,3),TpTzIICYzIITpYzIICYzIIYtII(3,3)

Complex(dp) :: Trmd2,Trme2,Trml2,Trmq2,Trmu2,TrYdadjYd,TrYeadjYe,TrYsIICYsII,TrYtIICYtII,            & 
& TrYuadjYu,TrYzIIadjYzII,TradjYdTd,TradjYeTe,TradjYuTu,TradjYzIITzII,TrCYsIITsII,       & 
& TrCYtIITtII,TrCTdTpTd,TrCTeTpTe,TrCTsIITsII,TrCTtIITtII,TrCTuTpTu,TrCTzIITpTzII,       & 
& Trmd2YdadjYd,Trmd2YsIICYsII,Trmd2YzIIadjYzII,Trme2YeadjYe,Trml2adjYeYe,Trml2adjYzIIYzII,& 
& Trml2CYtIIYtII,Trmq2adjYdYd,Trmq2adjYuYu,Trmu2YuadjYu

Complex(dp) :: TrYsIICTsII,TrYtIICTtII,TrCTdTpYd,TrCTeTpYe,TrCTuTpYu,TrCTzIITpYzII,Trmd2CYsIIYsII,   & 
& Trml2YtIICYtII,TrYdadjYdCmd2,TrYdCmq2adjYd,TrYeadjYeCme2,TrYeCml2adjYe,TrYuadjYuCmu2,  & 
& TrYuCmq2adjYu,TrYzIIadjYzIICmd2,TrYzIICml2adjYzII,TrYdadjYdYdadjYd,TrYdadjYdYsIICYsII, & 
& TrYdadjYdYzIIadjYzII,TrYdadjYdTdadjYd,TrYdadjYdTdadjTd,TrYdadjYdTsIICYsII,             & 
& TrYdadjYdTsIICTsII,TrYdadjYdTzIIadjYzII,TrYdadjYdTzIIadjTzII,TrYdadjYuYuadjYd,         & 
& TrYdadjYuTuadjYd,TrYdadjYuTuadjTd,TrYdadjTdTdadjYd,TrYdadjTdTsIICYsII,TrYdadjTdTzIIadjYzII,& 
& TrYdadjTuTuadjYd,TrYeadjYeYeadjYe,TrYeadjYeTeadjYe,TrYeadjYeTeadjTe,TrYeadjYzIIYzIIadjYe,& 
& TrYeadjYzIITzIIadjYe,TrYeadjYzIITzIIadjTe,TrYeadjTeTeadjYe,TrYeadjTzIITzIIadjYe,       & 
& TrYeCYtIIYtIIadjYe,TrYeCYtIITtIIadjYe,TrYeCYtIITtIIadjTe,TrYeCTtIITtIIadjYe,           & 
& TrYsIICYsIIYsIICYsII,TrYsIICYsIIYzIIadjYzII,TrYsIICYsIITdadjYd,TrYsIICYsIITdadjTd,     & 
& TrYsIICYsIITsIICYsII,TrYsIICYsIITsIICTsII,TrYsIICYsIITzIIadjYzII,TrYsIICYsIITzIIadjTzII,& 
& TrYsIICTdTpTdCYsII,TrYsIICTsIITdadjYd,TrYsIICTsIITsIICYsII,TrYsIICTsIITzIIadjYzII,     & 
& TrYsIICTzIITpTzIICYsII,TrYtIIadjYeTeCYtII,TrYtIIadjYeTeCTtII,TrYtIIadjYzIIYzIICYtII,   & 
& TrYtIIadjYzIITzIICYtII,TrYtIIadjYzIITzIICTtII,TrYtIIadjTeTeCYtII,TrYtIIadjTzIITzIICYtII,& 
& TrYtIICYtIIYtIICYtII,TrYtIICYtIITtIICYtII,TrYtIICYtIITtIICTtII,TrYtIICYtIITpTeCTe,     & 
& TrYtIICYtIITpTzIICTzII,TrYtIICTtIITtIICYtII,TrYuadjYdTdadjYu,TrYuadjYdTdadjTu,         & 
& TrYuadjYuYuadjYu,TrYuadjYuTuadjYu,TrYuadjYuTuadjTu,TrYuadjTdTdadjYu,TrYuadjTuTuadjYu,  & 
& TrYzIIadjYeTeadjYzII,TrYzIIadjYeTeadjTzII,TrYzIIadjYzIIYzIIadjYzII,TrYzIIadjYzIITdadjYd,& 
& TrYzIIadjYzIITdadjTd,TrYzIIadjYzIITsIICYsII,TrYzIIadjYzIITsIICTsII,TrYzIIadjYzIITzIIadjYzII,& 
& TrYzIIadjYzIITzIIadjTzII,TrYzIIadjTeTeadjYzII,TrYzIIadjTzIITdadjYd,TrYzIIadjTzIITsIICYsII,& 
& TrYzIIadjTzIITzIIadjYzII,TrYzIICYtIITtIIadjYzII,TrYzIICYtIITtIIadjTzII,TrYzIICTtIITtIIadjYzII,& 
& TrCYsIITsIICTdTpYd,TrCYsIITsIICTzIITpYzII,TrCYtIITpYeCTeTtII,TrCYtIITpYzIICTzIITtII,   & 
& Trmd2YdadjYdYdadjYd,Trmd2YdadjYdYsIICYsII,Trmd2YdadjYdYzIIadjYzII,Trmd2YdadjYuYuadjYd, & 
& Trmd2YsIICYsIIYdadjYd,Trmd2YsIICYsIIYsIICYsII,Trmd2YsIICYsIIYzIIadjYzII,               & 
& Trmd2YzIIadjYeYeadjYzII,Trmd2YzIIadjYzIIYdadjYd,Trmd2YzIIadjYzIIYsIICYsII,             & 
& Trmd2YzIIadjYzIIYzIIadjYzII,Trmd2YzIICYtIIYtIIadjYzII,Trme2YeadjYeYeadjYe,             & 
& Trme2YeadjYzIIYzIIadjYe,Trme2YeCYtIIYtIIadjYe,Trml2adjYeYeadjYeYe,Trml2adjYeYeadjYzIIYzII,& 
& Trml2adjYeYeCYtIIYtII,Trml2adjYzIIYdadjYdYzII,Trml2adjYzIIYsIICYsIIYzII,               & 
& Trml2adjYzIIYzIIadjYeYe,Trml2adjYzIIYzIIadjYzIIYzII,Trml2adjYzIIYzIICYtIIYtII,         & 
& Trml2CYtIIYtIIadjYeYe,Trml2CYtIIYtIIadjYzIIYzII,Trml2CYtIIYtIICYtIIYtII,               & 
& Trmq2adjYdYdadjYdYd,Trmq2adjYdYdadjYuYu,Trmq2adjYdYsIICYsIIYd,Trmq2adjYdYzIIadjYzIIYd, & 
& Trmq2adjYuYuadjYdYd,Trmq2adjYuYuadjYuYu,Trmu2YuadjYdYdadjYu,Trmu2YuadjYuYuadjYu,       & 
& TrYdadjYdYsIICmd2CYsII,TrYeCYtIICml2YtIIadjYe,TrYsIICmd2CYsIIYzIIadjYzII,              & 
& TrYtIIadjYzIIYzIICYtIICml2

Real(dp) :: g1p2,g1p3,g2p2,g2p3,g3p2,g3p3

Complex(dp) :: sqrt3ov5,ooSqrt15,sqrt15,L1IIp2,L2IIp2

Real(dp) :: g1p4,g2p4,g3p4

Complex(dp) :: CL1IIp2,CL2IIp2

Iname = Iname +1 
NameOfUnit(Iname) = 'rge365' 
 
OnlyDiagonal = .Not.GenerationMixing 
q = t 
 
Call GToParameters365(gy,g1,g2,g3,Yu,Yd,Ye,YtII,YsII,YzII,L1II,L2II,Mu,               & 
& MTII,MZII,MSII,Tu,Td,Te,TtII,TsII,TzII,TL1II,TL2II,Bmu,BMTII,BMZII,BMSII,              & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mt2,mtb2,ms2,msb2,mzz2,mzb2,M1,M2,M3,WOp)

AbsL1II = Abs(L1II)**2
AbsL2II = Abs(L2II)**2
AbsTL1II = Abs(TL1II)**2
AbsTL2II = Abs(TL2II)**2
AbsM1 = Abs(M1)**2
AbsM2 = Abs(M2)**2
AbsM3 = Abs(M3)**2
Call Adjungate(Yu,adjYu)
Call Adjungate(Yd,adjYd)
Call Adjungate(Ye,adjYe)
Call Adjungate(YtII,adjYtII)
Call Adjungate(YsII,adjYsII)
Call Adjungate(YzII,adjYzII)
Call Adjungate(Tu,adjTu)
Call Adjungate(Td,adjTd)
Call Adjungate(Te,adjTe)
Call Adjungate(TtII,adjTtII)
Call Adjungate(TsII,adjTsII)
Call Adjungate(TzII,adjTzII)
Call Adjungate(WOp,adjWOp)
 md2Yd = Matmul2(md2,Yd,OnlyDiagonal) 
 md2YzII = Matmul2(md2,YzII,OnlyDiagonal) 
 me2Ye = Matmul2(me2,Ye,OnlyDiagonal) 
 ml2adjYe = Matmul2(ml2,adjYe,OnlyDiagonal) 
 ml2adjYzII = Matmul2(ml2,adjYzII,OnlyDiagonal) 
 mq2adjYd = Matmul2(mq2,adjYd,OnlyDiagonal) 
 mq2adjYu = Matmul2(mq2,adjYu,OnlyDiagonal) 
 mu2Yu = Matmul2(mu2,Yu,OnlyDiagonal) 
 Ydmq2 = Matmul2(Yd,mq2,OnlyDiagonal) 
 YdadjYd = Matmul2(Yd,adjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYd(i2,i2) =  Real(YdadjYd(i2,i2),dp) 
 Yeml2 = Matmul2(Ye,ml2,OnlyDiagonal) 
 YeadjYe = Matmul2(Ye,adjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYe(i2,i2) =  Real(YeadjYe(i2,i2),dp) 
 YsIICYsII = Matmul2(YsII,adjYsII,OnlyDiagonal) 
 YtIIml2 = Matmul2(YtII,ml2,OnlyDiagonal) 
 YtIICYtII = Matmul2(YtII,adjYtII,OnlyDiagonal) 
 Yumq2 = Matmul2(Yu,mq2,OnlyDiagonal) 
 YuadjYu = Matmul2(Yu,adjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYu(i2,i2) =  Real(YuadjYu(i2,i2),dp) 
 YzIIml2 = Matmul2(YzII,ml2,OnlyDiagonal) 
 YzIIadjYzII = Matmul2(YzII,adjYzII,OnlyDiagonal) 
Forall(i2=1:3)  YzIIadjYzII(i2,i2) =  Real(YzIIadjYzII(i2,i2),dp) 
 adjYdmd2 = Matmul2(adjYd,md2,OnlyDiagonal) 
 adjYdYd = Matmul2(adjYd,Yd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYd(i2,i2) =  Real(adjYdYd(i2,i2),dp) 
 adjYdYsII = Matmul2(adjYd,YsII,OnlyDiagonal) 
 adjYdYzII = Matmul2(adjYd,YzII,OnlyDiagonal) 
 adjYdTd = Matmul2(adjYd,Td,OnlyDiagonal) 
 adjYdTsII = Matmul2(adjYd,TsII,OnlyDiagonal) 
 adjYdTzII = Matmul2(adjYd,TzII,OnlyDiagonal) 
 adjYeme2 = Matmul2(adjYe,me2,OnlyDiagonal) 
 adjYeYe = Matmul2(adjYe,Ye,OnlyDiagonal) 
Forall(i2=1:3)  adjYeYe(i2,i2) =  Real(adjYeYe(i2,i2),dp) 
 adjYeTe = Matmul2(adjYe,Te,OnlyDiagonal) 
 adjYumu2 = Matmul2(adjYu,mu2,OnlyDiagonal) 
 adjYuYu = Matmul2(adjYu,Yu,OnlyDiagonal) 
Forall(i2=1:3)  adjYuYu(i2,i2) =  Real(adjYuYu(i2,i2),dp) 
 adjYuTu = Matmul2(adjYu,Tu,OnlyDiagonal) 
 adjYzIImd2 = Matmul2(adjYzII,md2,OnlyDiagonal) 
 adjYzIIYd = Matmul2(adjYzII,Yd,OnlyDiagonal) 
 adjYzIIYsII = Matmul2(adjYzII,YsII,OnlyDiagonal) 
 adjYzIIYzII = Matmul2(adjYzII,YzII,OnlyDiagonal) 
Forall(i2=1:3)  adjYzIIYzII(i2,i2) =  Real(adjYzIIYzII(i2,i2),dp) 
 adjYzIITd = Matmul2(adjYzII,Td,OnlyDiagonal) 
 adjYzIITsII = Matmul2(adjYzII,TsII,OnlyDiagonal) 
 adjYzIITzII = Matmul2(adjYzII,TzII,OnlyDiagonal) 
 adjTdTd = Matmul2(adjTd,Td,OnlyDiagonal) 
 adjTeTe = Matmul2(adjTe,Te,OnlyDiagonal) 
 adjTuTu = Matmul2(adjTu,Tu,OnlyDiagonal) 
 adjTzIITzII = Matmul2(adjTzII,TzII,OnlyDiagonal) 
 Cmd2CYsII = Matmul2(Conjg(md2),adjYsII,OnlyDiagonal) 
 Cml2YtII = Matmul2(Conjg(ml2),YtII,OnlyDiagonal) 
 CYdTpYd = Matmul2(Conjg(Yd),Transpose(Yd),OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYd(i2,i2) =  Real(CYdTpYd(i2,i2),dp) 
 CYdTpTd = Matmul2(Conjg(Yd),Transpose(Td),OnlyDiagonal) 
 CYeWOp = Matmul2(Conjg(Ye),WOp,OnlyDiagonal) 
 CYeYtII = Matmul2(Conjg(Ye),YtII,OnlyDiagonal) 
 CYeTtII = Matmul2(Conjg(Ye),TtII,OnlyDiagonal) 
 CYsIImd2 = Matmul2(adjYsII,md2,OnlyDiagonal) 
 CYsIIYd = Matmul2(adjYsII,Yd,OnlyDiagonal) 
 CYsIIYsII = Matmul2(adjYsII,YsII,OnlyDiagonal) 
 CYsIIYzII = Matmul2(adjYsII,YzII,OnlyDiagonal) 
 CYsIITd = Matmul2(adjYsII,Td,OnlyDiagonal) 
 CYsIITsII = Matmul2(adjYsII,TsII,OnlyDiagonal) 
 CYsIITzII = Matmul2(adjYsII,TzII,OnlyDiagonal) 
 CYtIIWOp = Matmul2(adjYtII,WOp,OnlyDiagonal) 
 CYtIIYtII = Matmul2(adjYtII,YtII,OnlyDiagonal) 
 CYtIITtII = Matmul2(adjYtII,TtII,OnlyDiagonal) 
 CYzIIWOp = Matmul2(Conjg(YzII),WOp,OnlyDiagonal) 
 CYzIIYtII = Matmul2(Conjg(YzII),YtII,OnlyDiagonal) 
 CYzIITtII = Matmul2(Conjg(YzII),TtII,OnlyDiagonal) 
 CYzIITpYzII = Matmul2(Conjg(YzII),Transpose(YzII),OnlyDiagonal) 
Forall(i2=1:3)  CYzIITpYzII(i2,i2) =  Real(CYzIITpYzII(i2,i2),dp) 
 CYzIITpTzII = Matmul2(Conjg(YzII),Transpose(TzII),OnlyDiagonal) 
 CTdTpTd = Matmul2(Conjg(Td),Transpose(Td),OnlyDiagonal) 
 CTeTpTe = Matmul2(Conjg(Te),Transpose(Te),OnlyDiagonal) 
 CTsIITsII = Matmul2(adjTsII,TsII,OnlyDiagonal) 
 CTtIITtII = Matmul2(adjTtII,TtII,OnlyDiagonal) 
 CTuTpTu = Matmul2(Conjg(Tu),Transpose(Tu),OnlyDiagonal) 
 CTzIITpTzII = Matmul2(Conjg(TzII),Transpose(TzII),OnlyDiagonal) 
 TdadjTd = Matmul2(Td,adjTd,OnlyDiagonal) 
 TeadjTe = Matmul2(Te,adjTe,OnlyDiagonal) 
 TsIICTsII = Matmul2(TsII,adjTsII,OnlyDiagonal) 
 TuadjTu = Matmul2(Tu,adjTu,OnlyDiagonal) 
 TzIIadjTzII = Matmul2(TzII,adjTzII,OnlyDiagonal) 
 md2YdadjYd = Matmul2(md2,YdadjYd,OnlyDiagonal) 
 md2YsIICYsII = Matmul2(md2,YsIICYsII,OnlyDiagonal) 
 md2YzIIadjYzII = Matmul2(md2,YzIIadjYzII,OnlyDiagonal) 
 me2YeadjYe = Matmul2(me2,YeadjYe,OnlyDiagonal) 
 ml2adjYeYe = Matmul2(ml2,adjYeYe,OnlyDiagonal) 
 ml2adjYzIIYzII = Matmul2(ml2,adjYzIIYzII,OnlyDiagonal) 
 ml2CYtIIYtII = Matmul2(ml2,CYtIIYtII,OnlyDiagonal) 
 mq2adjYdYd = Matmul2(mq2,adjYdYd,OnlyDiagonal) 
 mq2adjYuYu = Matmul2(mq2,adjYuYu,OnlyDiagonal) 
 mu2YuadjYu = Matmul2(mu2,YuadjYu,OnlyDiagonal) 
 WOpadjYeYe = Matmul2(WOp,adjYeYe,OnlyDiagonal) 
 WOpadjYzIIYzII = Matmul2(WOp,adjYzIIYzII,OnlyDiagonal) 
 WOpCYtIIYtII = Matmul2(WOp,CYtIIYtII,OnlyDiagonal) 
 Ydmq2adjYd = Matmul2(Yd,mq2adjYd,OnlyDiagonal) 
Forall(i2=1:3)  Ydmq2adjYd(i2,i2) =  Real(Ydmq2adjYd(i2,i2),dp) 
 YdadjYdmd2 = Matmul2(Yd,adjYdmd2,OnlyDiagonal) 
 YdadjYdYd = Matmul2(Yd,adjYdYd,OnlyDiagonal) 
 YdadjYdYsII = Matmul2(Yd,adjYdYsII,OnlyDiagonal) 
 YdadjYdYzII = Matmul2(Yd,adjYdYzII,OnlyDiagonal) 
 YdadjYdTd = Matmul2(Yd,adjYdTd,OnlyDiagonal) 
 YdadjYdTsII = Matmul2(Yd,adjYdTsII,OnlyDiagonal) 
 YdadjYdTzII = Matmul2(Yd,adjYdTzII,OnlyDiagonal) 
 YdadjYuYu = Matmul2(Yd,adjYuYu,OnlyDiagonal) 
 YdadjYuTu = Matmul2(Yd,adjYuTu,OnlyDiagonal) 
 Yeml2adjYe = Matmul2(Ye,ml2adjYe,OnlyDiagonal) 
Forall(i2=1:3)  Yeml2adjYe(i2,i2) =  Real(Yeml2adjYe(i2,i2),dp) 
 YeadjYeme2 = Matmul2(Ye,adjYeme2,OnlyDiagonal) 
 YeadjYeYe = Matmul2(Ye,adjYeYe,OnlyDiagonal) 
 YeadjYeTe = Matmul2(Ye,adjYeTe,OnlyDiagonal) 
 YeadjYzIIYzII = Matmul2(Ye,adjYzIIYzII,OnlyDiagonal) 
 YeadjYzIITzII = Matmul2(Ye,adjYzIITzII,OnlyDiagonal) 
 YeCYtIIYtII = Matmul2(Ye,CYtIIYtII,OnlyDiagonal) 
 YeCYtIITtII = Matmul2(Ye,CYtIITtII,OnlyDiagonal) 
 YsIICmd2CYsII = Matmul2(YsII,Cmd2CYsII,OnlyDiagonal) 
 YsIICYdTpYd = Matmul2(YsII,CYdTpYd,OnlyDiagonal) 
 YsIICYdTpTd = Matmul2(YsII,CYdTpTd,OnlyDiagonal) 
 YsIICYsIImd2 = Matmul2(YsII,CYsIImd2,OnlyDiagonal) 
 YsIICYsIIYd = Matmul2(YsII,CYsIIYd,OnlyDiagonal) 
 YsIICYsIIYsII = Matmul2(YsII,CYsIIYsII,OnlyDiagonal) 
 YsIICYsIIYzII = Matmul2(YsII,CYsIIYzII,OnlyDiagonal) 
 YsIICYsIITd = Matmul2(YsII,CYsIITd,OnlyDiagonal) 
 YsIICYsIITsII = Matmul2(YsII,CYsIITsII,OnlyDiagonal) 
 YsIICYsIITzII = Matmul2(YsII,CYsIITzII,OnlyDiagonal) 
 YsIICYzIITpYzII = Matmul2(YsII,CYzIITpYzII,OnlyDiagonal) 
 YsIICYzIITpTzII = Matmul2(YsII,CYzIITpTzII,OnlyDiagonal) 
 YtIIadjYeYe = Matmul2(YtII,adjYeYe,OnlyDiagonal) 
 YtIIadjYeTe = Matmul2(YtII,adjYeTe,OnlyDiagonal) 
 YtIIadjYzIIYzII = Matmul2(YtII,adjYzIIYzII,OnlyDiagonal) 
 YtIIadjYzIITzII = Matmul2(YtII,adjYzIITzII,OnlyDiagonal) 
 YtIICYtIIWOp = Matmul2(YtII,CYtIIWOp,OnlyDiagonal) 
 YtIICYtIIYtII = Matmul2(YtII,CYtIIYtII,OnlyDiagonal) 
 YtIICYtIITtII = Matmul2(YtII,CYtIITtII,OnlyDiagonal) 
 Yumq2adjYu = Matmul2(Yu,mq2adjYu,OnlyDiagonal) 
Forall(i2=1:3)  Yumq2adjYu(i2,i2) =  Real(Yumq2adjYu(i2,i2),dp) 
 YuadjYdYd = Matmul2(Yu,adjYdYd,OnlyDiagonal) 
 YuadjYdTd = Matmul2(Yu,adjYdTd,OnlyDiagonal) 
 YuadjYumu2 = Matmul2(Yu,adjYumu2,OnlyDiagonal) 
 YuadjYuYu = Matmul2(Yu,adjYuYu,OnlyDiagonal) 
 YuadjYuTu = Matmul2(Yu,adjYuTu,OnlyDiagonal) 
 YzIIml2adjYzII = Matmul2(YzII,ml2adjYzII,OnlyDiagonal) 
Forall(i2=1:3)  YzIIml2adjYzII(i2,i2) =  Real(YzIIml2adjYzII(i2,i2),dp) 
 YzIIadjYeYe = Matmul2(YzII,adjYeYe,OnlyDiagonal) 
 YzIIadjYeTe = Matmul2(YzII,adjYeTe,OnlyDiagonal) 
 YzIIadjYzIImd2 = Matmul2(YzII,adjYzIImd2,OnlyDiagonal) 
 YzIIadjYzIIYd = Matmul2(YzII,adjYzIIYd,OnlyDiagonal) 
 YzIIadjYzIIYsII = Matmul2(YzII,adjYzIIYsII,OnlyDiagonal) 
 YzIIadjYzIIYzII = Matmul2(YzII,adjYzIIYzII,OnlyDiagonal) 
 YzIIadjYzIITd = Matmul2(YzII,adjYzIITd,OnlyDiagonal) 
 YzIIadjYzIITsII = Matmul2(YzII,adjYzIITsII,OnlyDiagonal) 
 YzIIadjYzIITzII = Matmul2(YzII,adjYzIITzII,OnlyDiagonal) 
 YzIICYtIIYtII = Matmul2(YzII,CYtIIYtII,OnlyDiagonal) 
 YzIICYtIITtII = Matmul2(YzII,CYtIITtII,OnlyDiagonal) 
 adjYdmd2Yd = Matmul2(adjYd,md2Yd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdmd2Yd(i2,i2) =  Real(adjYdmd2Yd(i2,i2),dp) 
 adjYdYdmq2 = Matmul2(adjYd,Ydmq2,OnlyDiagonal) 
 adjYeme2Ye = Matmul2(adjYe,me2Ye,OnlyDiagonal) 
Forall(i2=1:3)  adjYeme2Ye(i2,i2) =  Real(adjYeme2Ye(i2,i2),dp) 
 adjYeYeml2 = Matmul2(adjYe,Yeml2,OnlyDiagonal) 
 adjYumu2Yu = Matmul2(adjYu,mu2Yu,OnlyDiagonal) 
Forall(i2=1:3)  adjYumu2Yu(i2,i2) =  Real(adjYumu2Yu(i2,i2),dp) 
 adjYuYumq2 = Matmul2(adjYu,Yumq2,OnlyDiagonal) 
 adjYzIImd2YzII = Matmul2(adjYzII,md2YzII,OnlyDiagonal) 
Forall(i2=1:3)  adjYzIImd2YzII(i2,i2) =  Real(adjYzIImd2YzII(i2,i2),dp) 
 adjYzIIYzIIml2 = Matmul2(adjYzII,YzIIml2,OnlyDiagonal) 
 CYtIIYtIIml2 = Matmul2(adjYtII,YtIIml2,OnlyDiagonal) 
 CYtIICml2YtII = Matmul2(adjYtII,Cml2YtII,OnlyDiagonal) 
 TdadjYdYd = Matmul2(Td,adjYdYd,OnlyDiagonal) 
 TdadjYdYsII = Matmul2(Td,adjYdYsII,OnlyDiagonal) 
 TdadjYdYzII = Matmul2(Td,adjYdYzII,OnlyDiagonal) 
 TdadjYuYu = Matmul2(Td,adjYuYu,OnlyDiagonal) 
 TeadjYeYe = Matmul2(Te,adjYeYe,OnlyDiagonal) 
 TeadjYzIIYzII = Matmul2(Te,adjYzIIYzII,OnlyDiagonal) 
 TeCYtIIYtII = Matmul2(Te,CYtIIYtII,OnlyDiagonal) 
 TsIICYdTpYd = Matmul2(TsII,CYdTpYd,OnlyDiagonal) 
 TsIICYsIIYd = Matmul2(TsII,CYsIIYd,OnlyDiagonal) 
 TsIICYsIIYsII = Matmul2(TsII,CYsIIYsII,OnlyDiagonal) 
 TsIICYsIIYzII = Matmul2(TsII,CYsIIYzII,OnlyDiagonal) 
 TsIICYzIITpYzII = Matmul2(TsII,CYzIITpYzII,OnlyDiagonal) 
 TtIIadjYeYe = Matmul2(TtII,adjYeYe,OnlyDiagonal) 
 TtIIadjYzIIYzII = Matmul2(TtII,adjYzIIYzII,OnlyDiagonal) 
 TtIICYtIIYtII = Matmul2(TtII,CYtIIYtII,OnlyDiagonal) 
 TuadjYdYd = Matmul2(Tu,adjYdYd,OnlyDiagonal) 
 TuadjYuYu = Matmul2(Tu,adjYuYu,OnlyDiagonal) 
 TzIIadjYeYe = Matmul2(TzII,adjYeYe,OnlyDiagonal) 
 TzIIadjYzIIYd = Matmul2(TzII,adjYzIIYd,OnlyDiagonal) 
 TzIIadjYzIIYsII = Matmul2(TzII,adjYzIIYsII,OnlyDiagonal) 
 TzIIadjYzIIYzII = Matmul2(TzII,adjYzIIYzII,OnlyDiagonal) 
 TzIICYtIIYtII = Matmul2(TzII,CYtIIYtII,OnlyDiagonal) 
 TpYeCYeWOp = Matmul2(Transpose(Ye),CYeWOp,OnlyDiagonal) 
 TpYeCYeYtII = Matmul2(Transpose(Ye),CYeYtII,OnlyDiagonal) 
 TpYeCYeTtII = Matmul2(Transpose(Ye),CYeTtII,OnlyDiagonal) 
 TpYzIICYzIIWOp = Matmul2(Transpose(YzII),CYzIIWOp,OnlyDiagonal) 
 TpYzIICYzIIYtII = Matmul2(Transpose(YzII),CYzIIYtII,OnlyDiagonal) 
 TpYzIICYzIITtII = Matmul2(Transpose(YzII),CYzIITtII,OnlyDiagonal) 
 TpTeCYeYtII = Matmul2(Transpose(Te),CYeYtII,OnlyDiagonal) 
 TpTzIICYzIIYtII = Matmul2(Transpose(TzII),CYzIIYtII,OnlyDiagonal) 
 Trmd2 = Real(cTrace(md2),dp) 
 Trme2 = Real(cTrace(me2),dp) 
 Trml2 = Real(cTrace(ml2),dp) 
 Trmq2 = Real(cTrace(mq2),dp) 
 Trmu2 = Real(cTrace(mu2),dp) 
 TrYdadjYd = Real(cTrace(YdadjYd),dp) 
 TrYeadjYe = Real(cTrace(YeadjYe),dp) 
 TrYsIICYsII = Real(cTrace(YsIICYsII),dp) 
 TrYtIICYtII = Real(cTrace(YtIICYtII),dp) 
 TrYuadjYu = Real(cTrace(YuadjYu),dp) 
 TrYzIIadjYzII = Real(cTrace(YzIIadjYzII),dp) 
 TradjYdTd = Real(cTrace(adjYdTd),dp) 
 TradjYeTe = Real(cTrace(adjYeTe),dp) 
 TradjYuTu = Real(cTrace(adjYuTu),dp) 
 TradjYzIITzII = Real(cTrace(adjYzIITzII),dp) 
 TrCYsIITsII = Real(cTrace(CYsIITsII),dp) 
 TrCYtIITtII = Real(cTrace(CYtIITtII),dp) 
 TrCTdTpTd = Real(cTrace(CTdTpTd),dp) 
 TrCTeTpTe = Real(cTrace(CTeTpTe),dp) 
 TrCTsIITsII = Real(cTrace(CTsIITsII),dp) 
 TrCTtIITtII = Real(cTrace(CTtIITtII),dp) 
 TrCTuTpTu = Real(cTrace(CTuTpTu),dp) 
 TrCTzIITpTzII = Real(cTrace(CTzIITpTzII),dp) 
 Trmd2YdadjYd = Real(cTrace(md2YdadjYd),dp) 
 Trmd2YsIICYsII = Real(cTrace(md2YsIICYsII),dp) 
 Trmd2YzIIadjYzII = Real(cTrace(md2YzIIadjYzII),dp) 
 Trme2YeadjYe = Real(cTrace(me2YeadjYe),dp) 
 Trml2adjYeYe = Real(cTrace(ml2adjYeYe),dp) 
 Trml2adjYzIIYzII = Real(cTrace(ml2adjYzIIYzII),dp) 
 Trml2CYtIIYtII = Real(cTrace(ml2CYtIIYtII),dp) 
 Trmq2adjYdYd = Real(cTrace(mq2adjYdYd),dp) 
 Trmq2adjYuYu = Real(cTrace(mq2adjYuYu),dp) 
 Trmu2YuadjYu = Real(cTrace(mu2YuadjYu),dp) 
 sqrt3ov5 =Sqrt(3._dp/5._dp) 
 ooSqrt15 =1._dp/sqrt(15._dp) 
 sqrt15 =sqrt(15._dp) 
 g1p2 =g1**2 
 g1p3 =g1**3 
 g2p2 =g2**2 
 g2p3 =g2**3 
 g3p2 =g3**2 
 g3p3 =g3**3 
 L1IIp2 =L1II**2 
 L2IIp2 =L2II**2 
 g1p4 =g1**4 
 g2p4 =g2**4 
 g3p4 =g3**4 
 CL1IIp2 =Conjg(L1II)**2 
 CL2IIp2 =Conjg(L2II)**2 


If (TwoLoopRGE) Then 
 YdadjYu = Matmul2(Yd,adjYu,OnlyDiagonal) 
 YdadjTd = Matmul2(Yd,adjTd,OnlyDiagonal) 
 YdadjTu = Matmul2(Yd,adjTu,OnlyDiagonal) 
 YeadjYzII = Matmul2(Ye,adjYzII,OnlyDiagonal) 
 YeadjTe = Matmul2(Ye,adjTe,OnlyDiagonal) 
 YeadjTzII = Matmul2(Ye,adjTzII,OnlyDiagonal) 
 YsIICTsII = Matmul2(YsII,adjTsII,OnlyDiagonal) 
 YtIIadjYe = Matmul2(YtII,adjYe,OnlyDiagonal) 
 YtIIadjYzII = Matmul2(YtII,adjYzII,OnlyDiagonal) 
 YtIIadjTe = Matmul2(YtII,adjTe,OnlyDiagonal) 
 YtIIadjTzII = Matmul2(YtII,adjTzII,OnlyDiagonal) 
 YtIICTtII = Matmul2(YtII,adjTtII,OnlyDiagonal) 
 YuadjYd = Matmul2(Yu,adjYd,OnlyDiagonal) 
 YuadjTd = Matmul2(Yu,adjTd,OnlyDiagonal) 
 YuadjTu = Matmul2(Yu,adjTu,OnlyDiagonal) 
 YzIIadjYe = Matmul2(YzII,adjYe,OnlyDiagonal) 
 YzIIadjTe = Matmul2(YzII,adjTe,OnlyDiagonal) 
 YzIIadjTzII = Matmul2(YzII,adjTzII,OnlyDiagonal) 
 YzIICYtII = Matmul2(YzII,adjYtII,OnlyDiagonal) 
 adjYdCmd2 = Matmul2(adjYd,Conjg(md2),OnlyDiagonal) 
 adjYdCTsII = Matmul2(adjYd,adjTsII,OnlyDiagonal) 
 adjYeCme2 = Matmul2(adjYe,Conjg(me2),OnlyDiagonal) 
 adjYuCmu2 = Matmul2(adjYu,Conjg(mu2),OnlyDiagonal) 
 adjYzIICmd2 = Matmul2(adjYzII,Conjg(md2),OnlyDiagonal) 
 adjYzIICTsII = Matmul2(adjYzII,adjTsII,OnlyDiagonal) 
 adjTdYd = Matmul2(adjTd,Yd,OnlyDiagonal) 
 adjTdYzII = Matmul2(adjTd,YzII,OnlyDiagonal) 
 adjTdCYsII = Matmul2(adjTd,adjYsII,OnlyDiagonal) 
 adjTdTzII = Matmul2(adjTd,TzII,OnlyDiagonal) 
 adjTeYe = Matmul2(adjTe,Ye,OnlyDiagonal) 
 adjTuYu = Matmul2(adjTu,Yu,OnlyDiagonal) 
 adjTzIIYd = Matmul2(adjTzII,Yd,OnlyDiagonal) 
 adjTzIIYzII = Matmul2(adjTzII,YzII,OnlyDiagonal) 
 adjTzIICYsII = Matmul2(adjTzII,adjYsII,OnlyDiagonal) 
 adjTzIITd = Matmul2(adjTzII,Td,OnlyDiagonal) 
 Cml2adjYe = Matmul2(Conjg(ml2),adjYe,OnlyDiagonal) 
 Cml2adjYzII = Matmul2(Conjg(ml2),adjYzII,OnlyDiagonal) 
 Cmq2adjYd = Matmul2(Conjg(mq2),adjYd,OnlyDiagonal) 
 Cmq2adjYu = Matmul2(Conjg(mq2),adjYu,OnlyDiagonal) 
 CYeTpYzII = Matmul2(Conjg(Ye),Transpose(YzII),OnlyDiagonal) 
 CYeTpTzII = Matmul2(Conjg(Ye),Transpose(TzII),OnlyDiagonal) 
 CYtIICml2 = Matmul2(adjYtII,Conjg(ml2),OnlyDiagonal) 
 CYtIITpYzII = Matmul2(adjYtII,Transpose(YzII),OnlyDiagonal) 
 CYtIITpTzII = Matmul2(adjYtII,Transpose(TzII),OnlyDiagonal) 
 CYuTpYd = Matmul2(Conjg(Yu),Transpose(Yd),OnlyDiagonal) 
 CYuTpTd = Matmul2(Conjg(Yu),Transpose(Td),OnlyDiagonal) 
 CTdTpYd = Matmul2(Conjg(Td),Transpose(Yd),OnlyDiagonal) 
 CTeYtII = Matmul2(Conjg(Te),YtII,OnlyDiagonal) 
 CTeTtII = Matmul2(Conjg(Te),TtII,OnlyDiagonal) 
 CTeTpYe = Matmul2(Conjg(Te),Transpose(Ye),OnlyDiagonal) 
 CTsIIYd = Matmul2(adjTsII,Yd,OnlyDiagonal) 
 CTsIIYzII = Matmul2(adjTsII,YzII,OnlyDiagonal) 
 CTsIITd = Matmul2(adjTsII,Td,OnlyDiagonal) 
 CTsIITzII = Matmul2(adjTsII,TzII,OnlyDiagonal) 
 CTtIIYtII = Matmul2(adjTtII,YtII,OnlyDiagonal) 
 CTuTpYu = Matmul2(Conjg(Tu),Transpose(Yu),OnlyDiagonal) 
 CTzIIYtII = Matmul2(Conjg(TzII),YtII,OnlyDiagonal) 
 CTzIITtII = Matmul2(Conjg(TzII),TtII,OnlyDiagonal) 
 CTzIITpYzII = Matmul2(Conjg(TzII),Transpose(YzII),OnlyDiagonal) 
 TdadjYd = Matmul2(Td,adjYd,OnlyDiagonal) 
 TdadjYu = Matmul2(Td,adjYu,OnlyDiagonal) 
 TdadjTu = Matmul2(Td,adjTu,OnlyDiagonal) 
 TeadjYe = Matmul2(Te,adjYe,OnlyDiagonal) 
 TeadjYzII = Matmul2(Te,adjYzII,OnlyDiagonal) 
 TeadjTzII = Matmul2(Te,adjTzII,OnlyDiagonal) 
 TeCYtII = Matmul2(Te,adjYtII,OnlyDiagonal) 
 TeCTtII = Matmul2(Te,adjTtII,OnlyDiagonal) 
 TsIICYsII = Matmul2(TsII,adjYsII,OnlyDiagonal) 
 TtIIadjYe = Matmul2(TtII,adjYe,OnlyDiagonal) 
 TtIIadjYzII = Matmul2(TtII,adjYzII,OnlyDiagonal) 
 TtIIadjTe = Matmul2(TtII,adjTe,OnlyDiagonal) 
 TtIIadjTzII = Matmul2(TtII,adjTzII,OnlyDiagonal) 
 TtIICYtII = Matmul2(TtII,adjYtII,OnlyDiagonal) 
 TtIICTtII = Matmul2(TtII,adjTtII,OnlyDiagonal) 
 TuadjYd = Matmul2(Tu,adjYd,OnlyDiagonal) 
 TuadjYu = Matmul2(Tu,adjYu,OnlyDiagonal) 
 TuadjTd = Matmul2(Tu,adjTd,OnlyDiagonal) 
 TzIIadjYe = Matmul2(TzII,adjYe,OnlyDiagonal) 
 TzIIadjYzII = Matmul2(TzII,adjYzII,OnlyDiagonal) 
 TzIIadjTe = Matmul2(TzII,adjTe,OnlyDiagonal) 
 TzIICYtII = Matmul2(TzII,adjYtII,OnlyDiagonal) 
 TzIICTtII = Matmul2(TzII,adjTtII,OnlyDiagonal) 
 TpYdCYsII = Matmul2(Transpose(Yd),adjYsII,OnlyDiagonal) 
 TpYdCTsII = Matmul2(Transpose(Yd),adjTsII,OnlyDiagonal) 
 TpYzIICYsII = Matmul2(Transpose(YzII),adjYsII,OnlyDiagonal) 
 TpYzIICTsII = Matmul2(Transpose(YzII),adjTsII,OnlyDiagonal) 
 TpTdCYsII = Matmul2(Transpose(Td),adjYsII,OnlyDiagonal) 
 TpTdCTsII = Matmul2(Transpose(Td),adjTsII,OnlyDiagonal) 
 TpTeCTe = Matmul2(Transpose(Te),Conjg(Te),OnlyDiagonal) 
 TpTzIICYsII = Matmul2(Transpose(TzII),adjYsII,OnlyDiagonal) 
 TpTzIICTsII = Matmul2(Transpose(TzII),adjTsII,OnlyDiagonal) 
 TpTzIICTzII = Matmul2(Transpose(TzII),Conjg(TzII),OnlyDiagonal) 
 md2YdadjYu = Matmul2(md2,YdadjYu,OnlyDiagonal) 
 md2YzIIadjYe = Matmul2(md2,YzIIadjYe,OnlyDiagonal) 
 md2CYsIIYsII = Matmul2(md2,CYsIIYsII,OnlyDiagonal) 
 me2YeadjYzII = Matmul2(me2,YeadjYzII,OnlyDiagonal) 
 ml2YtIICYtII = Matmul2(ml2,YtIICYtII,OnlyDiagonal) 
 ml2adjYzIIYd = Matmul2(ml2,adjYzIIYd,OnlyDiagonal) 
 mq2adjYdYzII = Matmul2(mq2,adjYdYzII,OnlyDiagonal) 
 mu2YuadjYd = Matmul2(mu2,YuadjYd,OnlyDiagonal) 
 Ydmq2adjYu = Matmul2(Yd,mq2adjYu,OnlyDiagonal) 
 YdadjYdCmd2 = Matmul2(Yd,adjYdCmd2,OnlyDiagonal) 
 YdadjYumu2 = Matmul2(Yd,adjYumu2,OnlyDiagonal) 
 YdadjTdCYsII = Matmul2(Yd,adjTdCYsII,OnlyDiagonal) 
 YdadjTdTd = Matmul2(Yd,adjTdTd,OnlyDiagonal) 
 YdadjTdTzII = Matmul2(Yd,adjTdTzII,OnlyDiagonal) 
 YdCmq2adjYd = Matmul2(Yd,Cmq2adjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdCmq2adjYd(i2,i2) =  Real(YdCmq2adjYd(i2,i2),dp) 
 Yeml2adjYzII = Matmul2(Ye,ml2adjYzII,OnlyDiagonal) 
 YeadjYeCme2 = Matmul2(Ye,adjYeCme2,OnlyDiagonal) 
 YeadjYzIImd2 = Matmul2(Ye,adjYzIImd2,OnlyDiagonal) 
 YeadjYzIIYd = Matmul2(Ye,adjYzIIYd,OnlyDiagonal) 
 YeadjYzIIYsII = Matmul2(Ye,adjYzIIYsII,OnlyDiagonal) 
 YeadjYzIITd = Matmul2(Ye,adjYzIITd,OnlyDiagonal) 
 YeadjYzIITsII = Matmul2(Ye,adjYzIITsII,OnlyDiagonal) 
 YeadjTeTe = Matmul2(Ye,adjTeTe,OnlyDiagonal) 
 YeCml2adjYe = Matmul2(Ye,Cml2adjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeCml2adjYe(i2,i2) =  Real(YeCml2adjYe(i2,i2),dp) 
 YeCYtIIWOp = Matmul2(Ye,CYtIIWOp,OnlyDiagonal) 
 YsIICYzIIWOp = Matmul2(YsII,CYzIIWOp,OnlyDiagonal) 
 YsIICYzIIYtII = Matmul2(YsII,CYzIIYtII,OnlyDiagonal) 
 YsIICYzIITtII = Matmul2(YsII,CYzIITtII,OnlyDiagonal) 
 YsIICTsIITd = Matmul2(YsII,CTsIITd,OnlyDiagonal) 
 YsIICTsIITzII = Matmul2(YsII,CTsIITzII,OnlyDiagonal) 
 YtIIml2adjYe = Matmul2(YtII,ml2adjYe,OnlyDiagonal) 
 YtIIml2adjYzII = Matmul2(YtII,ml2adjYzII,OnlyDiagonal) 
 YtIIadjYeme2 = Matmul2(YtII,adjYeme2,OnlyDiagonal) 
 YtIIadjYzIImd2 = Matmul2(YtII,adjYzIImd2,OnlyDiagonal) 
 YtIIadjYzIIYd = Matmul2(YtII,adjYzIIYd,OnlyDiagonal) 
 YtIIadjYzIIYsII = Matmul2(YtII,adjYzIIYsII,OnlyDiagonal) 
 YtIIadjYzIITd = Matmul2(YtII,adjYzIITd,OnlyDiagonal) 
 YtIIadjYzIITsII = Matmul2(YtII,adjYzIITsII,OnlyDiagonal) 
 YtIICYtIITpYzII = Matmul2(YtII,CYtIITpYzII,OnlyDiagonal) 
 YtIICYtIITpTzII = Matmul2(YtII,CYtIITpTzII,OnlyDiagonal) 
 YtIICTtIITtII = Matmul2(YtII,CTtIITtII,OnlyDiagonal) 
 Yumq2adjYd = Matmul2(Yu,mq2adjYd,OnlyDiagonal) 
 YuadjYdmd2 = Matmul2(Yu,adjYdmd2,OnlyDiagonal) 
 YuadjYdYsII = Matmul2(Yu,adjYdYsII,OnlyDiagonal) 
 YuadjYdYzII = Matmul2(Yu,adjYdYzII,OnlyDiagonal) 
 YuadjYdTsII = Matmul2(Yu,adjYdTsII,OnlyDiagonal) 
 YuadjYdTzII = Matmul2(Yu,adjYdTzII,OnlyDiagonal) 
 YuadjYuCmu2 = Matmul2(Yu,adjYuCmu2,OnlyDiagonal) 
 YuadjTuTu = Matmul2(Yu,adjTuTu,OnlyDiagonal) 
 YuCmq2adjYu = Matmul2(Yu,Cmq2adjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuCmq2adjYu(i2,i2) =  Real(YuCmq2adjYu(i2,i2),dp) 
 YzIIml2adjYe = Matmul2(YzII,ml2adjYe,OnlyDiagonal) 
 YzIIadjYeme2 = Matmul2(YzII,adjYeme2,OnlyDiagonal) 
 YzIIadjYzIICmd2 = Matmul2(YzII,adjYzIICmd2,OnlyDiagonal) 
 YzIIadjTzIICYsII = Matmul2(YzII,adjTzIICYsII,OnlyDiagonal) 
 YzIIadjTzIITd = Matmul2(YzII,adjTzIITd,OnlyDiagonal) 
 YzIIadjTzIITzII = Matmul2(YzII,adjTzIITzII,OnlyDiagonal) 
 YzIICml2adjYzII = Matmul2(YzII,Cml2adjYzII,OnlyDiagonal) 
Forall(i2=1:3)  YzIICml2adjYzII(i2,i2) =  Real(YzIICml2adjYzII(i2,i2),dp) 
 YzIICYtIIWOp = Matmul2(YzII,CYtIIWOp,OnlyDiagonal) 
 YzIICYtIICml2 = Matmul2(YzII,CYtIICml2,OnlyDiagonal) 
 adjYdmd2YzII = Matmul2(adjYd,md2YzII,OnlyDiagonal) 
 adjYdYdadjYd = Matmul2(adjYd,YdadjYd,OnlyDiagonal) 
 adjYdYdadjYu = Matmul2(adjYd,YdadjYu,OnlyDiagonal) 
 adjYdYdadjTd = Matmul2(adjYd,YdadjTd,OnlyDiagonal) 
 adjYdYdadjTu = Matmul2(adjYd,YdadjTu,OnlyDiagonal) 
 adjYdYsIICYsII = Matmul2(adjYd,YsIICYsII,OnlyDiagonal) 
 adjYdYsIICTsII = Matmul2(adjYd,YsIICTsII,OnlyDiagonal) 
 adjYdYzIIml2 = Matmul2(adjYd,YzIIml2,OnlyDiagonal) 
 adjYdYzIIadjYzII = Matmul2(adjYd,YzIIadjYzII,OnlyDiagonal) 
 adjYdCYsIIYd = Matmul2(adjYd,CYsIIYd,OnlyDiagonal) 
 adjYdCYsIIYzII = Matmul2(adjYd,CYsIIYzII,OnlyDiagonal) 
 adjYdTdadjYd = Matmul2(adjYd,TdadjYd,OnlyDiagonal) 
 adjYdTdadjYu = Matmul2(adjYd,TdadjYu,OnlyDiagonal) 
 adjYdTdadjTd = Matmul2(adjYd,TdadjTd,OnlyDiagonal) 
 adjYdTdadjTu = Matmul2(adjYd,TdadjTu,OnlyDiagonal) 
 adjYdTsIICYsII = Matmul2(adjYd,TsIICYsII,OnlyDiagonal) 
 adjYdTsIICTsII = Matmul2(adjYd,TsIICTsII,OnlyDiagonal) 
 adjYdTzIIadjYzII = Matmul2(adjYd,TzIIadjYzII,OnlyDiagonal) 
 adjYdTzIIadjTzII = Matmul2(adjYd,TzIIadjTzII,OnlyDiagonal) 
 adjYeYeadjYe = Matmul2(adjYe,YeadjYe,OnlyDiagonal) 
 adjYeYeadjYzII = Matmul2(adjYe,YeadjYzII,OnlyDiagonal) 
 adjYeYeadjTe = Matmul2(adjYe,YeadjTe,OnlyDiagonal) 
 adjYeYeadjTzII = Matmul2(adjYe,YeadjTzII,OnlyDiagonal) 
 adjYeTeadjYe = Matmul2(adjYe,TeadjYe,OnlyDiagonal) 
 adjYeTeadjYzII = Matmul2(adjYe,TeadjYzII,OnlyDiagonal) 
 adjYeTeadjTe = Matmul2(adjYe,TeadjTe,OnlyDiagonal) 
 adjYeTeadjTzII = Matmul2(adjYe,TeadjTzII,OnlyDiagonal) 
 adjYeTeCYtII = Matmul2(adjYe,TeCYtII,OnlyDiagonal) 
 adjYeTeCTtII = Matmul2(adjYe,TeCTtII,OnlyDiagonal) 
 adjYuYuadjYd = Matmul2(adjYu,YuadjYd,OnlyDiagonal) 
 adjYuYuadjYu = Matmul2(adjYu,YuadjYu,OnlyDiagonal) 
 adjYuYuadjTd = Matmul2(adjYu,YuadjTd,OnlyDiagonal) 
 adjYuYuadjTu = Matmul2(adjYu,YuadjTu,OnlyDiagonal) 
 adjYuTuadjYd = Matmul2(adjYu,TuadjYd,OnlyDiagonal) 
 adjYuTuadjYu = Matmul2(adjYu,TuadjYu,OnlyDiagonal) 
 adjYuTuadjTd = Matmul2(adjYu,TuadjTd,OnlyDiagonal) 
 adjYuTuadjTu = Matmul2(adjYu,TuadjTu,OnlyDiagonal) 
 adjYzIImd2Yd = Matmul2(adjYzII,md2Yd,OnlyDiagonal) 
 adjYzIIYdmq2 = Matmul2(adjYzII,Ydmq2,OnlyDiagonal) 
 adjYzIIYdadjYd = Matmul2(adjYzII,YdadjYd,OnlyDiagonal) 
 adjYzIIYsIICYsII = Matmul2(adjYzII,YsIICYsII,OnlyDiagonal) 
 adjYzIIYsIICTsII = Matmul2(adjYzII,YsIICTsII,OnlyDiagonal) 
 adjYzIIYzIIadjYe = Matmul2(adjYzII,YzIIadjYe,OnlyDiagonal) 
 adjYzIIYzIIadjYzII = Matmul2(adjYzII,YzIIadjYzII,OnlyDiagonal) 
 adjYzIIYzIIadjTe = Matmul2(adjYzII,YzIIadjTe,OnlyDiagonal) 
 adjYzIIYzIIadjTzII = Matmul2(adjYzII,YzIIadjTzII,OnlyDiagonal) 
 adjYzIIYzIICYtII = Matmul2(adjYzII,YzIICYtII,OnlyDiagonal) 
 adjYzIICYsIIYd = Matmul2(adjYzII,CYsIIYd,OnlyDiagonal) 
 adjYzIICYsIIYzII = Matmul2(adjYzII,CYsIIYzII,OnlyDiagonal) 
 adjYzIITdadjYd = Matmul2(adjYzII,TdadjYd,OnlyDiagonal) 
 adjYzIITdadjTd = Matmul2(adjYzII,TdadjTd,OnlyDiagonal) 
 adjYzIITsIICYsII = Matmul2(adjYzII,TsIICYsII,OnlyDiagonal) 
 adjYzIITsIICTsII = Matmul2(adjYzII,TsIICTsII,OnlyDiagonal) 
 adjYzIITzIIadjYe = Matmul2(adjYzII,TzIIadjYe,OnlyDiagonal) 
 adjYzIITzIIadjYzII = Matmul2(adjYzII,TzIIadjYzII,OnlyDiagonal) 
 adjYzIITzIIadjTe = Matmul2(adjYzII,TzIIadjTe,OnlyDiagonal) 
 adjYzIITzIIadjTzII = Matmul2(adjYzII,TzIIadjTzII,OnlyDiagonal) 
 adjYzIITzIICYtII = Matmul2(adjYzII,TzIICYtII,OnlyDiagonal) 
 adjYzIITzIICTtII = Matmul2(adjYzII,TzIICTtII,OnlyDiagonal) 
 adjTdYdadjYd = Matmul2(adjTd,YdadjYd,OnlyDiagonal) 
 adjTdYdadjYu = Matmul2(adjTd,YdadjYu,OnlyDiagonal) 
 adjTdTdadjYd = Matmul2(adjTd,TdadjYd,OnlyDiagonal) 
 adjTdTdadjYu = Matmul2(adjTd,TdadjYu,OnlyDiagonal) 
 adjTdTsIICYsII = Matmul2(adjTd,TsIICYsII,OnlyDiagonal) 
 adjTdTzIIadjYzII = Matmul2(adjTd,TzIIadjYzII,OnlyDiagonal) 
 adjTeYeadjYe = Matmul2(adjTe,YeadjYe,OnlyDiagonal) 
 adjTeYeadjYzII = Matmul2(adjTe,YeadjYzII,OnlyDiagonal) 
 adjTeTeadjYe = Matmul2(adjTe,TeadjYe,OnlyDiagonal) 
 adjTeTeadjYzII = Matmul2(adjTe,TeadjYzII,OnlyDiagonal) 
 adjTeTeCYtII = Matmul2(adjTe,TeCYtII,OnlyDiagonal) 
 adjTuYuadjYd = Matmul2(adjTu,YuadjYd,OnlyDiagonal) 
 adjTuYuadjYu = Matmul2(adjTu,YuadjYu,OnlyDiagonal) 
 adjTuTuadjYd = Matmul2(adjTu,TuadjYd,OnlyDiagonal) 
 adjTuTuadjYu = Matmul2(adjTu,TuadjYu,OnlyDiagonal) 
 adjTzIIYzIIadjYe = Matmul2(adjTzII,YzIIadjYe,OnlyDiagonal) 
 adjTzIIYzIIadjYzII = Matmul2(adjTzII,YzIIadjYzII,OnlyDiagonal) 
 adjTzIITdadjYd = Matmul2(adjTzII,TdadjYd,OnlyDiagonal) 
 adjTzIITsIICYsII = Matmul2(adjTzII,TsIICYsII,OnlyDiagonal) 
 adjTzIITzIIadjYe = Matmul2(adjTzII,TzIIadjYe,OnlyDiagonal) 
 adjTzIITzIIadjYzII = Matmul2(adjTzII,TzIIadjYzII,OnlyDiagonal) 
 adjTzIITzIICYtII = Matmul2(adjTzII,TzIICYtII,OnlyDiagonal) 
 Cmd2CYsIIYd = Matmul2(Conjg(md2),CYsIIYd,OnlyDiagonal) 
 Cmd2CYsIIYzII = Matmul2(Conjg(md2),CYsIIYzII,OnlyDiagonal) 
 Cmd2CYzIIYtII = Matmul2(Conjg(md2),CYzIIYtII,OnlyDiagonal) 
 Cme2CYeYtII = Matmul2(Conjg(me2),CYeYtII,OnlyDiagonal) 
 Cml2YtIIadjYe = Matmul2(Conjg(ml2),YtIIadjYe,OnlyDiagonal) 
 Cml2YtIIadjYzII = Matmul2(Conjg(ml2),YtIIadjYzII,OnlyDiagonal) 
 Cml2TpYzIICYsII = Matmul2(Conjg(ml2),TpYzIICYsII,OnlyDiagonal) 
 Cmq2TpYdCYsII = Matmul2(Conjg(mq2),TpYdCYsII,OnlyDiagonal) 
 CYdTpYdCYsII = Matmul2(Conjg(Yd),TpYdCYsII,OnlyDiagonal) 
 CYdTpYdCTsII = Matmul2(Conjg(Yd),TpYdCTsII,OnlyDiagonal) 
 CYdTpTdCTsII = Matmul2(Conjg(Yd),TpTdCTsII,OnlyDiagonal) 
 CYeYtIIml2 = Matmul2(Conjg(Ye),YtIIml2,OnlyDiagonal) 
 CYeCml2YtII = Matmul2(Conjg(Ye),Cml2YtII,OnlyDiagonal) 
 CYsIImd2Yd = Matmul2(adjYsII,md2Yd,OnlyDiagonal) 
 CYsIImd2YzII = Matmul2(adjYsII,md2YzII,OnlyDiagonal) 
 CYsIIYdmq2 = Matmul2(adjYsII,Ydmq2,OnlyDiagonal) 
 CYsIIYdadjYd = Matmul2(adjYsII,YdadjYd,OnlyDiagonal) 
 CYsIIYsIICYsII = Matmul2(adjYsII,YsIICYsII,OnlyDiagonal) 
 CYsIIYsIICTsII = Matmul2(adjYsII,YsIICTsII,OnlyDiagonal) 
 CYsIIYzIIml2 = Matmul2(adjYsII,YzIIml2,OnlyDiagonal) 
 CYsIIYzIIadjYzII = Matmul2(adjYsII,YzIIadjYzII,OnlyDiagonal) 
 CYsIITdadjYd = Matmul2(adjYsII,TdadjYd,OnlyDiagonal) 
 CYsIITdadjTd = Matmul2(adjYsII,TdadjTd,OnlyDiagonal) 
 CYsIITsIICYsII = Matmul2(adjYsII,TsIICYsII,OnlyDiagonal) 
 CYsIITsIICTsII = Matmul2(adjYsII,TsIICTsII,OnlyDiagonal) 
 CYsIITzIIadjYzII = Matmul2(adjYsII,TzIIadjYzII,OnlyDiagonal) 
 CYsIITzIIadjTzII = Matmul2(adjYsII,TzIIadjTzII,OnlyDiagonal) 
 CYtIIYtIIadjYe = Matmul2(adjYtII,YtIIadjYe,OnlyDiagonal) 
 CYtIIYtIIadjYzII = Matmul2(adjYtII,YtIIadjYzII,OnlyDiagonal) 
 CYtIIYtIIadjTe = Matmul2(adjYtII,YtIIadjTe,OnlyDiagonal) 
 CYtIIYtIIadjTzII = Matmul2(adjYtII,YtIIadjTzII,OnlyDiagonal) 
 CYtIIYtIICYtII = Matmul2(adjYtII,YtIICYtII,OnlyDiagonal) 
 CYtIITtIIadjYe = Matmul2(adjYtII,TtIIadjYe,OnlyDiagonal) 
 CYtIITtIIadjYzII = Matmul2(adjYtII,TtIIadjYzII,OnlyDiagonal) 
 CYtIITtIIadjTe = Matmul2(adjYtII,TtIIadjTe,OnlyDiagonal) 
 CYtIITtIIadjTzII = Matmul2(adjYtII,TtIIadjTzII,OnlyDiagonal) 
 CYtIITtIICYtII = Matmul2(adjYtII,TtIICYtII,OnlyDiagonal) 
 CYtIITtIICTtII = Matmul2(adjYtII,TtIICTtII,OnlyDiagonal) 
 CYtIITpTeCTe = Matmul2(adjYtII,TpTeCTe,OnlyDiagonal) 
 CYtIITpTzIICTzII = Matmul2(adjYtII,TpTzIICTzII,OnlyDiagonal) 
 CYzIIYtIIml2 = Matmul2(Conjg(YzII),YtIIml2,OnlyDiagonal) 
 CYzIICml2YtII = Matmul2(Conjg(YzII),Cml2YtII,OnlyDiagonal) 
 CYzIITpYzIICYsII = Matmul2(Conjg(YzII),TpYzIICYsII,OnlyDiagonal) 
 CYzIITpYzIICTsII = Matmul2(Conjg(YzII),TpYzIICTsII,OnlyDiagonal) 
 CYzIITpTzIICTsII = Matmul2(Conjg(YzII),TpTzIICTsII,OnlyDiagonal) 
 CTdTpYdCYsII = Matmul2(Conjg(Td),TpYdCYsII,OnlyDiagonal) 
 CTdTpTdCYsII = Matmul2(Conjg(Td),TpTdCYsII,OnlyDiagonal) 
 CTsIIYsIICYsII = Matmul2(adjTsII,YsIICYsII,OnlyDiagonal) 
 CTsIITdadjYd = Matmul2(adjTsII,TdadjYd,OnlyDiagonal) 
 CTsIITsIICYsII = Matmul2(adjTsII,TsIICYsII,OnlyDiagonal) 
 CTsIITzIIadjYzII = Matmul2(adjTsII,TzIIadjYzII,OnlyDiagonal) 
 CTtIIYtIIadjYe = Matmul2(adjTtII,YtIIadjYe,OnlyDiagonal) 
 CTtIIYtIIadjYzII = Matmul2(adjTtII,YtIIadjYzII,OnlyDiagonal) 
 CTtIITtIIadjYe = Matmul2(adjTtII,TtIIadjYe,OnlyDiagonal) 
 CTtIITtIIadjYzII = Matmul2(adjTtII,TtIIadjYzII,OnlyDiagonal) 
 CTtIITtIICYtII = Matmul2(adjTtII,TtIICYtII,OnlyDiagonal) 
 CTzIITpYzIICYsII = Matmul2(Conjg(TzII),TpYzIICYsII,OnlyDiagonal) 
 CTzIITpTzIICYsII = Matmul2(Conjg(TzII),TpTzIICYsII,OnlyDiagonal) 
 TdadjYdCTsII = Matmul2(Td,adjYdCTsII,OnlyDiagonal) 
 TdadjTdYd = Matmul2(Td,adjTdYd,OnlyDiagonal) 
 TdadjTdYzII = Matmul2(Td,adjTdYzII,OnlyDiagonal) 
 TeadjYzIIYd = Matmul2(Te,adjYzIIYd,OnlyDiagonal) 
 TeadjYzIIYsII = Matmul2(Te,adjYzIIYsII,OnlyDiagonal) 
 TeadjTeYe = Matmul2(Te,adjTeYe,OnlyDiagonal) 
 TsIICYzIIYtII = Matmul2(TsII,CYzIIYtII,OnlyDiagonal) 
 TsIICTdTpYd = Matmul2(TsII,CTdTpYd,OnlyDiagonal) 
 TsIICTsIIYd = Matmul2(TsII,CTsIIYd,OnlyDiagonal) 
 TsIICTsIIYzII = Matmul2(TsII,CTsIIYzII,OnlyDiagonal) 
 TsIICTzIITpYzII = Matmul2(TsII,CTzIITpYzII,OnlyDiagonal) 
 TtIIadjYzIIYd = Matmul2(TtII,adjYzIIYd,OnlyDiagonal) 
 TtIIadjYzIIYsII = Matmul2(TtII,adjYzIIYsII,OnlyDiagonal) 
 TtIICYtIITpYzII = Matmul2(TtII,CYtIITpYzII,OnlyDiagonal) 
 TtIICTtIIYtII = Matmul2(TtII,CTtIIYtII,OnlyDiagonal) 
 TuadjYdYsII = Matmul2(Tu,adjYdYsII,OnlyDiagonal) 
 TuadjYdYzII = Matmul2(Tu,adjYdYzII,OnlyDiagonal) 
 TuadjTuYu = Matmul2(Tu,adjTuYu,OnlyDiagonal) 
 TzIIadjYzIICTsII = Matmul2(TzII,adjYzIICTsII,OnlyDiagonal) 
 TzIIadjTzIIYd = Matmul2(TzII,adjTzIIYd,OnlyDiagonal) 
 TzIIadjTzIIYzII = Matmul2(TzII,adjTzIIYzII,OnlyDiagonal) 
 TpYdCmd2CYsII = Matmul2(Transpose(Yd),Cmd2CYsII,OnlyDiagonal) 
 TpYdCYdTpYd = Matmul2(Transpose(Yd),CYdTpYd,OnlyDiagonal) 
 TpYdCYdTpTd = Matmul2(Transpose(Yd),CYdTpTd,OnlyDiagonal) 
 TpYdCYsIImd2 = Matmul2(Transpose(Yd),CYsIImd2,OnlyDiagonal) 
 TpYdCYsIIYd = Matmul2(Transpose(Yd),CYsIIYd,OnlyDiagonal) 
 TpYdCYsIIYsII = Matmul2(Transpose(Yd),CYsIIYsII,OnlyDiagonal) 
 TpYdCYsIIYzII = Matmul2(Transpose(Yd),CYsIIYzII,OnlyDiagonal) 
 TpYdCYsIITd = Matmul2(Transpose(Yd),CYsIITd,OnlyDiagonal) 
 TpYdCYsIITsII = Matmul2(Transpose(Yd),CYsIITsII,OnlyDiagonal) 
 TpYdCYsIITzII = Matmul2(Transpose(Yd),CYsIITzII,OnlyDiagonal) 
 TpYdCYzIIWOp = Matmul2(Transpose(Yd),CYzIIWOp,OnlyDiagonal) 
 TpYdCYzIIYtII = Matmul2(Transpose(Yd),CYzIIYtII,OnlyDiagonal) 
 TpYdCYzIITtII = Matmul2(Transpose(Yd),CYzIITtII,OnlyDiagonal) 
 TpYeCYeTpYzII = Matmul2(Transpose(Ye),CYeTpYzII,OnlyDiagonal) 
 TpYeCYeTpTzII = Matmul2(Transpose(Ye),CYeTpTzII,OnlyDiagonal) 
 TpYeCTeTtII = Matmul2(Transpose(Ye),CTeTtII,OnlyDiagonal) 
 TpYuCYuTpYd = Matmul2(Transpose(Yu),CYuTpYd,OnlyDiagonal) 
 TpYuCYuTpTd = Matmul2(Transpose(Yu),CYuTpTd,OnlyDiagonal) 
 TpYzIICmd2CYsII = Matmul2(Transpose(YzII),Cmd2CYsII,OnlyDiagonal) 
 TpYzIICYsIImd2 = Matmul2(Transpose(YzII),CYsIImd2,OnlyDiagonal) 
 TpYzIICYsIIYd = Matmul2(Transpose(YzII),CYsIIYd,OnlyDiagonal) 
 TpYzIICYsIIYsII = Matmul2(Transpose(YzII),CYsIIYsII,OnlyDiagonal) 
 TpYzIICYsIIYzII = Matmul2(Transpose(YzII),CYsIIYzII,OnlyDiagonal) 
 TpYzIICYsIITd = Matmul2(Transpose(YzII),CYsIITd,OnlyDiagonal) 
 TpYzIICYsIITsII = Matmul2(Transpose(YzII),CYsIITsII,OnlyDiagonal) 
 TpYzIICYsIITzII = Matmul2(Transpose(YzII),CYsIITzII,OnlyDiagonal) 
 TpYzIICYzIITpYzII = Matmul2(Transpose(YzII),CYzIITpYzII,OnlyDiagonal) 
 TpYzIICYzIITpTzII = Matmul2(Transpose(YzII),CYzIITpTzII,OnlyDiagonal) 
 TpYzIICTzIITtII = Matmul2(Transpose(YzII),CTzIITtII,OnlyDiagonal) 
 TpTdCYdTpYd = Matmul2(Transpose(Td),CYdTpYd,OnlyDiagonal) 
 TpTdCYsIIYd = Matmul2(Transpose(Td),CYsIIYd,OnlyDiagonal) 
 TpTdCYsIIYsII = Matmul2(Transpose(Td),CYsIIYsII,OnlyDiagonal) 
 TpTdCYsIIYzII = Matmul2(Transpose(Td),CYsIIYzII,OnlyDiagonal) 
 TpTdCYzIIYtII = Matmul2(Transpose(Td),CYzIIYtII,OnlyDiagonal) 
 TpTeCYeTpYzII = Matmul2(Transpose(Te),CYeTpYzII,OnlyDiagonal) 
 TpTeCTeYtII = Matmul2(Transpose(Te),CTeYtII,OnlyDiagonal) 
 TpTuCYuTpYd = Matmul2(Transpose(Tu),CYuTpYd,OnlyDiagonal) 
 TpTzIICYsIIYd = Matmul2(Transpose(TzII),CYsIIYd,OnlyDiagonal) 
 TpTzIICYsIIYsII = Matmul2(Transpose(TzII),CYsIIYsII,OnlyDiagonal) 
 TpTzIICYsIIYzII = Matmul2(Transpose(TzII),CYsIIYzII,OnlyDiagonal) 
 TpTzIICYzIITpYzII = Matmul2(Transpose(TzII),CYzIITpYzII,OnlyDiagonal) 
 TpTzIICTzIIYtII = Matmul2(Transpose(TzII),CTzIIYtII,OnlyDiagonal) 
 md2YdadjYdYd = Matmul2(md2,YdadjYdYd,OnlyDiagonal) 
 md2YdadjYdYzII = Matmul2(md2,YdadjYdYzII,OnlyDiagonal) 
 md2YsIICYsIIYd = Matmul2(md2,YsIICYsIIYd,OnlyDiagonal) 
 md2YsIICYsIIYzII = Matmul2(md2,YsIICYsIIYzII,OnlyDiagonal) 
 md2YzIIadjYzIIYd = Matmul2(md2,YzIIadjYzIIYd,OnlyDiagonal) 
 md2YzIIadjYzIIYzII = Matmul2(md2,YzIIadjYzIIYzII,OnlyDiagonal) 
 me2YeadjYeYe = Matmul2(me2,YeadjYeYe,OnlyDiagonal) 
 ml2adjYeYeadjYe = Matmul2(ml2,adjYeYeadjYe,OnlyDiagonal) 
 ml2adjYeYeadjYzII = Matmul2(ml2,adjYeYeadjYzII,OnlyDiagonal) 
 ml2adjYzIIYzIIadjYe = Matmul2(ml2,adjYzIIYzIIadjYe,OnlyDiagonal) 
 ml2adjYzIIYzIIadjYzII = Matmul2(ml2,adjYzIIYzIIadjYzII,OnlyDiagonal) 
 ml2CYtIIYtIIadjYe = Matmul2(ml2,CYtIIYtIIadjYe,OnlyDiagonal) 
 ml2CYtIIYtIIadjYzII = Matmul2(ml2,CYtIIYtIIadjYzII,OnlyDiagonal) 
 mq2adjYdYdadjYd = Matmul2(mq2,adjYdYdadjYd,OnlyDiagonal) 
 mq2adjYdYdadjYu = Matmul2(mq2,adjYdYdadjYu,OnlyDiagonal) 
 mq2adjYuYuadjYd = Matmul2(mq2,adjYuYuadjYd,OnlyDiagonal) 
 mq2adjYuYuadjYu = Matmul2(mq2,adjYuYuadjYu,OnlyDiagonal) 
 mu2YuadjYuYu = Matmul2(mu2,YuadjYuYu,OnlyDiagonal) 
 Ydmq2adjYdYd = Matmul2(Yd,mq2adjYdYd,OnlyDiagonal) 
 Ydmq2adjYdYzII = Matmul2(Yd,mq2adjYdYzII,OnlyDiagonal) 
 YdadjYdmd2Yd = Matmul2(Yd,adjYdmd2Yd,OnlyDiagonal) 
 YdadjYdmd2YzII = Matmul2(Yd,adjYdmd2YzII,OnlyDiagonal) 
 YdadjYdYdmq2 = Matmul2(Yd,adjYdYdmq2,OnlyDiagonal) 
 YdadjYdYdadjYd = Matmul2(Yd,adjYdYdadjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYdYdadjYd(i2,i2) =  Real(YdadjYdYdadjYd(i2,i2),dp) 
 YdadjYdYsIICYsII = Matmul2(Yd,adjYdYsIICYsII,OnlyDiagonal) 
 YdadjYdYzIIml2 = Matmul2(Yd,adjYdYzIIml2,OnlyDiagonal) 
 YdadjYdYzIIadjYzII = Matmul2(Yd,adjYdYzIIadjYzII,OnlyDiagonal) 
 YdadjYdTdadjYd = Matmul2(Yd,adjYdTdadjYd,OnlyDiagonal) 
 YdadjYdTdadjTd = Matmul2(Yd,adjYdTdadjTd,OnlyDiagonal) 
 YdadjYdTsIICYsII = Matmul2(Yd,adjYdTsIICYsII,OnlyDiagonal) 
 YdadjYdTsIICTsII = Matmul2(Yd,adjYdTsIICTsII,OnlyDiagonal) 
 YdadjYdTzIIadjYzII = Matmul2(Yd,adjYdTzIIadjYzII,OnlyDiagonal) 
 YdadjYdTzIIadjTzII = Matmul2(Yd,adjYdTzIIadjTzII,OnlyDiagonal) 
 YdadjYuYuadjYd = Matmul2(Yd,adjYuYuadjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYuYuadjYd(i2,i2) =  Real(YdadjYuYuadjYd(i2,i2),dp) 
 YdadjYuTuadjYd = Matmul2(Yd,adjYuTuadjYd,OnlyDiagonal) 
 YdadjYuTuadjTd = Matmul2(Yd,adjYuTuadjTd,OnlyDiagonal) 
 YdadjTdTdadjYd = Matmul2(Yd,adjTdTdadjYd,OnlyDiagonal) 
 YdadjTdTsIICYsII = Matmul2(Yd,adjTdTsIICYsII,OnlyDiagonal) 
 YdadjTdTzIIadjYzII = Matmul2(Yd,adjTdTzIIadjYzII,OnlyDiagonal) 
 YdadjTuTuadjYd = Matmul2(Yd,adjTuTuadjYd,OnlyDiagonal) 
 Yeml2adjYeYe = Matmul2(Ye,ml2adjYeYe,OnlyDiagonal) 
 YeadjYeme2Ye = Matmul2(Ye,adjYeme2Ye,OnlyDiagonal) 
 YeadjYeYeml2 = Matmul2(Ye,adjYeYeml2,OnlyDiagonal) 
 YeadjYeYeadjYe = Matmul2(Ye,adjYeYeadjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYeYeadjYe(i2,i2) =  Real(YeadjYeYeadjYe(i2,i2),dp) 
 YeadjYeTeadjYe = Matmul2(Ye,adjYeTeadjYe,OnlyDiagonal) 
 YeadjYeTeadjTe = Matmul2(Ye,adjYeTeadjTe,OnlyDiagonal) 
 YeadjYzIIYzIIadjYe = Matmul2(Ye,adjYzIIYzIIadjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYzIIYzIIadjYe(i2,i2) =  Real(YeadjYzIIYzIIadjYe(i2,i2),dp) 
 YeadjYzIITzIIadjYe = Matmul2(Ye,adjYzIITzIIadjYe,OnlyDiagonal) 
 YeadjYzIITzIIadjTe = Matmul2(Ye,adjYzIITzIIadjTe,OnlyDiagonal) 
 YeadjTeTeadjYe = Matmul2(Ye,adjTeTeadjYe,OnlyDiagonal) 
 YeadjTzIITzIIadjYe = Matmul2(Ye,adjTzIITzIIadjYe,OnlyDiagonal) 
 YeCYtIIYtIIadjYe = Matmul2(Ye,CYtIIYtIIadjYe,OnlyDiagonal) 
 YeCYtIITtIIadjYe = Matmul2(Ye,CYtIITtIIadjYe,OnlyDiagonal) 
 YeCYtIITtIIadjTe = Matmul2(Ye,CYtIITtIIadjTe,OnlyDiagonal) 
 YeCTtIITtIIadjYe = Matmul2(Ye,CTtIITtIIadjYe,OnlyDiagonal) 
 YsIICmd2CYsIIYd = Matmul2(YsII,Cmd2CYsIIYd,OnlyDiagonal) 
 YsIICmd2CYsIIYzII = Matmul2(YsII,Cmd2CYsIIYzII,OnlyDiagonal) 
 YsIICYdTpYdCYsII = Matmul2(YsII,CYdTpYdCYsII,OnlyDiagonal) 
 YsIICYdTpTdCTsII = Matmul2(YsII,CYdTpTdCTsII,OnlyDiagonal) 
 YsIICYsIImd2Yd = Matmul2(YsII,CYsIImd2Yd,OnlyDiagonal) 
 YsIICYsIImd2YzII = Matmul2(YsII,CYsIImd2YzII,OnlyDiagonal) 
 YsIICYsIIYdmq2 = Matmul2(YsII,CYsIIYdmq2,OnlyDiagonal) 
 YsIICYsIIYdadjYd = Matmul2(YsII,CYsIIYdadjYd,OnlyDiagonal) 
 YsIICYsIIYsIICYsII = Matmul2(YsII,CYsIIYsIICYsII,OnlyDiagonal) 
 YsIICYsIIYzIIml2 = Matmul2(YsII,CYsIIYzIIml2,OnlyDiagonal) 
 YsIICYsIIYzIIadjYzII = Matmul2(YsII,CYsIIYzIIadjYzII,OnlyDiagonal) 
 YsIICYsIITdadjYd = Matmul2(YsII,CYsIITdadjYd,OnlyDiagonal) 
 YsIICYsIITdadjTd = Matmul2(YsII,CYsIITdadjTd,OnlyDiagonal) 
 YsIICYsIITsIICYsII = Matmul2(YsII,CYsIITsIICYsII,OnlyDiagonal) 
 YsIICYsIITsIICTsII = Matmul2(YsII,CYsIITsIICTsII,OnlyDiagonal) 
 YsIICYsIITzIIadjYzII = Matmul2(YsII,CYsIITzIIadjYzII,OnlyDiagonal) 
 YsIICYsIITzIIadjTzII = Matmul2(YsII,CYsIITzIIadjTzII,OnlyDiagonal) 
 YsIICYzIITpYzIICYsII = Matmul2(YsII,CYzIITpYzIICYsII,OnlyDiagonal) 
 YsIICYzIITpTzIICTsII = Matmul2(YsII,CYzIITpTzIICTsII,OnlyDiagonal) 
 YsIICTdTpTdCYsII = Matmul2(YsII,CTdTpTdCYsII,OnlyDiagonal) 
 YsIICTsIITdadjYd = Matmul2(YsII,CTsIITdadjYd,OnlyDiagonal) 
 YsIICTsIITsIICYsII = Matmul2(YsII,CTsIITsIICYsII,OnlyDiagonal) 
 YsIICTsIITzIIadjYzII = Matmul2(YsII,CTsIITzIIadjYzII,OnlyDiagonal) 
 YsIICTzIITpTzIICYsII = Matmul2(YsII,CTzIITpTzIICYsII,OnlyDiagonal) 
 YsIITdadjYdCTsII = Matmul2(YsII,TdadjYdCTsII,OnlyDiagonal) 
 YsIITzIIadjYzIICTsII = Matmul2(YsII,TzIIadjYzIICTsII,OnlyDiagonal) 
 YtIIml2CYtIIYtII = Matmul2(YtII,ml2CYtIIYtII,OnlyDiagonal) 
 YtIIadjYeTeCYtII = Matmul2(YtII,adjYeTeCYtII,OnlyDiagonal) 
 YtIIadjYeTeCTtII = Matmul2(YtII,adjYeTeCTtII,OnlyDiagonal) 
 YtIIadjYzIIYzIICYtII = Matmul2(YtII,adjYzIIYzIICYtII,OnlyDiagonal) 
 YtIIadjYzIITzIICYtII = Matmul2(YtII,adjYzIITzIICYtII,OnlyDiagonal) 
 YtIIadjYzIITzIICTtII = Matmul2(YtII,adjYzIITzIICTtII,OnlyDiagonal) 
 YtIIadjTeTeCYtII = Matmul2(YtII,adjTeTeCYtII,OnlyDiagonal) 
 YtIIadjTzIITzIICYtII = Matmul2(YtII,adjTzIITzIICYtII,OnlyDiagonal) 
 YtIICYtIIYtIIml2 = Matmul2(YtII,CYtIIYtIIml2,OnlyDiagonal) 
 YtIICYtIIYtIICYtII = Matmul2(YtII,CYtIIYtIICYtII,OnlyDiagonal) 
 YtIICYtIICml2YtII = Matmul2(YtII,CYtIICml2YtII,OnlyDiagonal) 
 YtIICYtIITtIICYtII = Matmul2(YtII,CYtIITtIICYtII,OnlyDiagonal) 
 YtIICYtIITtIICTtII = Matmul2(YtII,CYtIITtIICTtII,OnlyDiagonal) 
 YtIICYtIITpTeCTe = Matmul2(YtII,CYtIITpTeCTe,OnlyDiagonal) 
 YtIICYtIITpTzIICTzII = Matmul2(YtII,CYtIITpTzIICTzII,OnlyDiagonal) 
 YtIICTtIITtIICYtII = Matmul2(YtII,CTtIITtIICYtII,OnlyDiagonal) 
 Yumq2adjYuYu = Matmul2(Yu,mq2adjYuYu,OnlyDiagonal) 
 YuadjYdYdadjYu = Matmul2(Yu,adjYdYdadjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYdYdadjYu(i2,i2) =  Real(YuadjYdYdadjYu(i2,i2),dp) 
 YuadjYdTdadjYu = Matmul2(Yu,adjYdTdadjYu,OnlyDiagonal) 
 YuadjYdTdadjTu = Matmul2(Yu,adjYdTdadjTu,OnlyDiagonal) 
 YuadjYumu2Yu = Matmul2(Yu,adjYumu2Yu,OnlyDiagonal) 
 YuadjYuYumq2 = Matmul2(Yu,adjYuYumq2,OnlyDiagonal) 
 YuadjYuYuadjYu = Matmul2(Yu,adjYuYuadjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYuYuadjYu(i2,i2) =  Real(YuadjYuYuadjYu(i2,i2),dp) 
 YuadjYuTuadjYu = Matmul2(Yu,adjYuTuadjYu,OnlyDiagonal) 
 YuadjYuTuadjTu = Matmul2(Yu,adjYuTuadjTu,OnlyDiagonal) 
 YuadjTdTdadjYu = Matmul2(Yu,adjTdTdadjYu,OnlyDiagonal) 
 YuadjTuTuadjYu = Matmul2(Yu,adjTuTuadjYu,OnlyDiagonal) 
 YzIIml2adjYzIIYd = Matmul2(YzII,ml2adjYzIIYd,OnlyDiagonal) 
 YzIIml2adjYzIIYzII = Matmul2(YzII,ml2adjYzIIYzII,OnlyDiagonal) 
 YzIIadjYeYeadjYzII = Matmul2(YzII,adjYeYeadjYzII,OnlyDiagonal) 
Forall(i2=1:3)  YzIIadjYeYeadjYzII(i2,i2) =  Real(YzIIadjYeYeadjYzII(i2,i2),dp) 
 YzIIadjYeTeadjYzII = Matmul2(YzII,adjYeTeadjYzII,OnlyDiagonal) 
 YzIIadjYeTeadjTzII = Matmul2(YzII,adjYeTeadjTzII,OnlyDiagonal) 
 YzIIadjYzIImd2Yd = Matmul2(YzII,adjYzIImd2Yd,OnlyDiagonal) 
 YzIIadjYzIImd2YzII = Matmul2(YzII,adjYzIImd2YzII,OnlyDiagonal) 
 YzIIadjYzIIYdmq2 = Matmul2(YzII,adjYzIIYdmq2,OnlyDiagonal) 
 YzIIadjYzIIYdadjYd = Matmul2(YzII,adjYzIIYdadjYd,OnlyDiagonal) 
 YzIIadjYzIIYsIICYsII = Matmul2(YzII,adjYzIIYsIICYsII,OnlyDiagonal) 
 YzIIadjYzIIYzIIml2 = Matmul2(YzII,adjYzIIYzIIml2,OnlyDiagonal) 
 YzIIadjYzIIYzIIadjYzII = Matmul2(YzII,adjYzIIYzIIadjYzII,OnlyDiagonal) 
Forall(i2=1:3)  YzIIadjYzIIYzIIadjYzII(i2,i2) =  Real(YzIIadjYzIIYzIIadjYzII(i2,i2),dp) 
 YzIIadjYzIITdadjYd = Matmul2(YzII,adjYzIITdadjYd,OnlyDiagonal) 
 YzIIadjYzIITdadjTd = Matmul2(YzII,adjYzIITdadjTd,OnlyDiagonal) 
 YzIIadjYzIITsIICYsII = Matmul2(YzII,adjYzIITsIICYsII,OnlyDiagonal) 
 YzIIadjYzIITsIICTsII = Matmul2(YzII,adjYzIITsIICTsII,OnlyDiagonal) 
 YzIIadjYzIITzIIadjYzII = Matmul2(YzII,adjYzIITzIIadjYzII,OnlyDiagonal) 
 YzIIadjYzIITzIIadjTzII = Matmul2(YzII,adjYzIITzIIadjTzII,OnlyDiagonal) 
 YzIIadjTeTeadjYzII = Matmul2(YzII,adjTeTeadjYzII,OnlyDiagonal) 
 YzIIadjTzIITdadjYd = Matmul2(YzII,adjTzIITdadjYd,OnlyDiagonal) 
 YzIIadjTzIITsIICYsII = Matmul2(YzII,adjTzIITsIICYsII,OnlyDiagonal) 
 YzIIadjTzIITzIIadjYzII = Matmul2(YzII,adjTzIITzIIadjYzII,OnlyDiagonal) 
 YzIICYtIIYtIIadjYzII = Matmul2(YzII,CYtIIYtIIadjYzII,OnlyDiagonal) 
 YzIICYtIITtIIadjYzII = Matmul2(YzII,CYtIITtIIadjYzII,OnlyDiagonal) 
 YzIICYtIITtIIadjTzII = Matmul2(YzII,CYtIITtIIadjTzII,OnlyDiagonal) 
 YzIICTtIITtIIadjYzII = Matmul2(YzII,CTtIITtIIadjYzII,OnlyDiagonal) 
 adjYdmd2YdadjYd = Matmul2(adjYd,md2YdadjYd,OnlyDiagonal) 
 adjYdmd2YdadjYu = Matmul2(adjYd,md2YdadjYu,OnlyDiagonal) 
 adjYdYdmq2adjYd = Matmul2(adjYd,Ydmq2adjYd,OnlyDiagonal) 
 adjYdYdmq2adjYu = Matmul2(adjYd,Ydmq2adjYu,OnlyDiagonal) 
 adjYdYdadjYdmd2 = Matmul2(adjYd,YdadjYdmd2,OnlyDiagonal) 
 adjYdYdadjYdYd = Matmul2(adjYd,YdadjYdYd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYdadjYdYd(i2,i2) =  Real(adjYdYdadjYdYd(i2,i2),dp) 
 adjYdYdadjYdYsII = Matmul2(adjYd,YdadjYdYsII,OnlyDiagonal) 
 adjYdYdadjYdYzII = Matmul2(adjYd,YdadjYdYzII,OnlyDiagonal) 
 adjYdYdadjYdTd = Matmul2(adjYd,YdadjYdTd,OnlyDiagonal) 
 adjYdYdadjYdTsII = Matmul2(adjYd,YdadjYdTsII,OnlyDiagonal) 
 adjYdYdadjYdTzII = Matmul2(adjYd,YdadjYdTzII,OnlyDiagonal) 
 adjYdYdadjYumu2 = Matmul2(adjYd,YdadjYumu2,OnlyDiagonal) 
 adjYdYdadjYuYu = Matmul2(adjYd,YdadjYuYu,OnlyDiagonal) 
 adjYdYdadjYuTu = Matmul2(adjYd,YdadjYuTu,OnlyDiagonal) 
 adjYdYdadjTdTd = Matmul2(adjYd,YdadjTdTd,OnlyDiagonal) 
 adjYdYsIICmd2CYsII = Matmul2(adjYd,YsIICmd2CYsII,OnlyDiagonal) 
 adjYdYsIICYsIIYd = Matmul2(adjYd,YsIICYsIIYd,OnlyDiagonal) 
 adjYdYsIICYsIITd = Matmul2(adjYd,YsIICYsIITd,OnlyDiagonal) 
 adjYdYsIICTsIITd = Matmul2(adjYd,YsIICTsIITd,OnlyDiagonal) 
 adjYdYzIIadjYzIIYd = Matmul2(adjYd,YzIIadjYzIIYd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYzIIadjYzIIYd(i2,i2) =  Real(adjYdYzIIadjYzIIYd(i2,i2),dp) 
 adjYdYzIIadjYzIITd = Matmul2(adjYd,YzIIadjYzIITd,OnlyDiagonal) 
 adjYdYzIIadjTzIITd = Matmul2(adjYd,YzIIadjTzIITd,OnlyDiagonal) 
 adjYdTdadjYdYd = Matmul2(adjYd,TdadjYdYd,OnlyDiagonal) 
 adjYdTdadjYdYsII = Matmul2(adjYd,TdadjYdYsII,OnlyDiagonal) 
 adjYdTdadjYdYzII = Matmul2(adjYd,TdadjYdYzII,OnlyDiagonal) 
 adjYdTdadjYuYu = Matmul2(adjYd,TdadjYuYu,OnlyDiagonal) 
 adjYdTdadjTdYd = Matmul2(adjYd,TdadjTdYd,OnlyDiagonal) 
 adjYdTsIICYsIIYd = Matmul2(adjYd,TsIICYsIIYd,OnlyDiagonal) 
 adjYdTsIICTsIIYd = Matmul2(adjYd,TsIICTsIIYd,OnlyDiagonal) 
 adjYdTzIIadjYzIIYd = Matmul2(adjYd,TzIIadjYzIIYd,OnlyDiagonal) 
 adjYdTzIIadjTzIIYd = Matmul2(adjYd,TzIIadjTzIIYd,OnlyDiagonal) 
 adjYeme2YeadjYe = Matmul2(adjYe,me2YeadjYe,OnlyDiagonal) 
 adjYeme2YeadjYzII = Matmul2(adjYe,me2YeadjYzII,OnlyDiagonal) 
 adjYeYeml2adjYe = Matmul2(adjYe,Yeml2adjYe,OnlyDiagonal) 
 adjYeYeml2adjYzII = Matmul2(adjYe,Yeml2adjYzII,OnlyDiagonal) 
 adjYeYeadjYeme2 = Matmul2(adjYe,YeadjYeme2,OnlyDiagonal) 
 adjYeYeadjYeYe = Matmul2(adjYe,YeadjYeYe,OnlyDiagonal) 
Forall(i2=1:3)  adjYeYeadjYeYe(i2,i2) =  Real(adjYeYeadjYeYe(i2,i2),dp) 
 adjYeYeadjYeTe = Matmul2(adjYe,YeadjYeTe,OnlyDiagonal) 
 adjYeYeadjYzIImd2 = Matmul2(adjYe,YeadjYzIImd2,OnlyDiagonal) 
 adjYeYeadjYzIIYd = Matmul2(adjYe,YeadjYzIIYd,OnlyDiagonal) 
 adjYeYeadjYzIIYsII = Matmul2(adjYe,YeadjYzIIYsII,OnlyDiagonal) 
 adjYeYeadjYzIIYzII = Matmul2(adjYe,YeadjYzIIYzII,OnlyDiagonal) 
 adjYeYeadjYzIITd = Matmul2(adjYe,YeadjYzIITd,OnlyDiagonal) 
 adjYeYeadjYzIITsII = Matmul2(adjYe,YeadjYzIITsII,OnlyDiagonal) 
 adjYeYeadjYzIITzII = Matmul2(adjYe,YeadjYzIITzII,OnlyDiagonal) 
 adjYeYeadjTeTe = Matmul2(adjYe,YeadjTeTe,OnlyDiagonal) 
 adjYeYeCYtIIWOp = Matmul2(adjYe,YeCYtIIWOp,OnlyDiagonal) 
 adjYeYeCYtIIYtII = Matmul2(adjYe,YeCYtIIYtII,OnlyDiagonal) 
 adjYeYeCYtIITtII = Matmul2(adjYe,YeCYtIITtII,OnlyDiagonal) 
 adjYeTeadjYeYe = Matmul2(adjYe,TeadjYeYe,OnlyDiagonal) 
 adjYeTeadjYzIIYd = Matmul2(adjYe,TeadjYzIIYd,OnlyDiagonal) 
 adjYeTeadjYzIIYsII = Matmul2(adjYe,TeadjYzIIYsII,OnlyDiagonal) 
 adjYeTeadjYzIIYzII = Matmul2(adjYe,TeadjYzIIYzII,OnlyDiagonal) 
 adjYeTeadjTeYe = Matmul2(adjYe,TeadjTeYe,OnlyDiagonal) 
 adjYeTeCYtIIYtII = Matmul2(adjYe,TeCYtIIYtII,OnlyDiagonal) 
 adjYumu2YuadjYd = Matmul2(adjYu,mu2YuadjYd,OnlyDiagonal) 
 adjYumu2YuadjYu = Matmul2(adjYu,mu2YuadjYu,OnlyDiagonal) 
 adjYuYumq2adjYd = Matmul2(adjYu,Yumq2adjYd,OnlyDiagonal) 
 adjYuYumq2adjYu = Matmul2(adjYu,Yumq2adjYu,OnlyDiagonal) 
 adjYuYuadjYdmd2 = Matmul2(adjYu,YuadjYdmd2,OnlyDiagonal) 
 adjYuYuadjYdYd = Matmul2(adjYu,YuadjYdYd,OnlyDiagonal) 
 adjYuYuadjYdYsII = Matmul2(adjYu,YuadjYdYsII,OnlyDiagonal) 
 adjYuYuadjYdYzII = Matmul2(adjYu,YuadjYdYzII,OnlyDiagonal) 
 adjYuYuadjYdTd = Matmul2(adjYu,YuadjYdTd,OnlyDiagonal) 
 adjYuYuadjYdTsII = Matmul2(adjYu,YuadjYdTsII,OnlyDiagonal) 
 adjYuYuadjYdTzII = Matmul2(adjYu,YuadjYdTzII,OnlyDiagonal) 
 adjYuYuadjYumu2 = Matmul2(adjYu,YuadjYumu2,OnlyDiagonal) 
 adjYuYuadjYuYu = Matmul2(adjYu,YuadjYuYu,OnlyDiagonal) 
Forall(i2=1:3)  adjYuYuadjYuYu(i2,i2) =  Real(adjYuYuadjYuYu(i2,i2),dp) 
 adjYuYuadjYuTu = Matmul2(adjYu,YuadjYuTu,OnlyDiagonal) 
 adjYuYuadjTuTu = Matmul2(adjYu,YuadjTuTu,OnlyDiagonal) 
 adjYuTuadjYdYd = Matmul2(adjYu,TuadjYdYd,OnlyDiagonal) 
 adjYuTuadjYdYsII = Matmul2(adjYu,TuadjYdYsII,OnlyDiagonal) 
 adjYuTuadjYdYzII = Matmul2(adjYu,TuadjYdYzII,OnlyDiagonal) 
 adjYuTuadjYuYu = Matmul2(adjYu,TuadjYuYu,OnlyDiagonal) 
 adjYuTuadjTuYu = Matmul2(adjYu,TuadjTuYu,OnlyDiagonal) 
 adjYzIImd2YzIIadjYe = Matmul2(adjYzII,md2YzIIadjYe,OnlyDiagonal) 
 adjYzIImd2YzIIadjYzII = Matmul2(adjYzII,md2YzIIadjYzII,OnlyDiagonal) 
 adjYzIIYdadjYdYzII = Matmul2(adjYzII,YdadjYdYzII,OnlyDiagonal) 
Forall(i2=1:3)  adjYzIIYdadjYdYzII(i2,i2) =  Real(adjYzIIYdadjYdYzII(i2,i2),dp) 
 adjYzIIYdadjYdTzII = Matmul2(adjYzII,YdadjYdTzII,OnlyDiagonal) 
 adjYzIIYdadjTdTzII = Matmul2(adjYzII,YdadjTdTzII,OnlyDiagonal) 
 adjYzIIYsIICYsIIYzII = Matmul2(adjYzII,YsIICYsIIYzII,OnlyDiagonal) 
 adjYzIIYsIICYsIITzII = Matmul2(adjYzII,YsIICYsIITzII,OnlyDiagonal) 
 adjYzIIYsIICTsIITzII = Matmul2(adjYzII,YsIICTsIITzII,OnlyDiagonal) 
 adjYzIIYzIIml2adjYe = Matmul2(adjYzII,YzIIml2adjYe,OnlyDiagonal) 
 adjYzIIYzIIml2adjYzII = Matmul2(adjYzII,YzIIml2adjYzII,OnlyDiagonal) 
 adjYzIIYzIIadjYeme2 = Matmul2(adjYzII,YzIIadjYeme2,OnlyDiagonal) 
 adjYzIIYzIIadjYeYe = Matmul2(adjYzII,YzIIadjYeYe,OnlyDiagonal) 
 adjYzIIYzIIadjYeTe = Matmul2(adjYzII,YzIIadjYeTe,OnlyDiagonal) 
 adjYzIIYzIIadjYzIImd2 = Matmul2(adjYzII,YzIIadjYzIImd2,OnlyDiagonal) 
 adjYzIIYzIIadjYzIIYd = Matmul2(adjYzII,YzIIadjYzIIYd,OnlyDiagonal) 
 adjYzIIYzIIadjYzIIYsII = Matmul2(adjYzII,YzIIadjYzIIYsII,OnlyDiagonal) 
 adjYzIIYzIIadjYzIIYzII = Matmul2(adjYzII,YzIIadjYzIIYzII,OnlyDiagonal) 
Forall(i2=1:3)  adjYzIIYzIIadjYzIIYzII(i2,i2) =  Real(adjYzIIYzIIadjYzIIYzII(i2,i2),dp) 
 adjYzIIYzIIadjYzIITd = Matmul2(adjYzII,YzIIadjYzIITd,OnlyDiagonal) 
 adjYzIIYzIIadjYzIITsII = Matmul2(adjYzII,YzIIadjYzIITsII,OnlyDiagonal) 
 adjYzIIYzIIadjYzIITzII = Matmul2(adjYzII,YzIIadjYzIITzII,OnlyDiagonal) 
 adjYzIIYzIIadjTzIITzII = Matmul2(adjYzII,YzIIadjTzIITzII,OnlyDiagonal) 
 adjYzIIYzIICYtIIWOp = Matmul2(adjYzII,YzIICYtIIWOp,OnlyDiagonal) 
 adjYzIIYzIICYtIIYtII = Matmul2(adjYzII,YzIICYtIIYtII,OnlyDiagonal) 
 adjYzIIYzIICYtIICml2 = Matmul2(adjYzII,YzIICYtIICml2,OnlyDiagonal) 
 adjYzIIYzIICYtIITtII = Matmul2(adjYzII,YzIICYtIITtII,OnlyDiagonal) 
 adjYzIITdadjYdYzII = Matmul2(adjYzII,TdadjYdYzII,OnlyDiagonal) 
 adjYzIITdadjTdYzII = Matmul2(adjYzII,TdadjTdYzII,OnlyDiagonal) 
 adjYzIITsIICYsIIYzII = Matmul2(adjYzII,TsIICYsIIYzII,OnlyDiagonal) 
 adjYzIITsIICTsIIYzII = Matmul2(adjYzII,TsIICTsIIYzII,OnlyDiagonal) 
 adjYzIITzIIadjYeYe = Matmul2(adjYzII,TzIIadjYeYe,OnlyDiagonal) 
 adjYzIITzIIadjYzIIYd = Matmul2(adjYzII,TzIIadjYzIIYd,OnlyDiagonal) 
 adjYzIITzIIadjYzIIYsII = Matmul2(adjYzII,TzIIadjYzIIYsII,OnlyDiagonal) 
 adjYzIITzIIadjYzIIYzII = Matmul2(adjYzII,TzIIadjYzIIYzII,OnlyDiagonal) 
 adjYzIITzIIadjTzIIYzII = Matmul2(adjYzII,TzIIadjTzIIYzII,OnlyDiagonal) 
 adjYzIITzIICYtIIYtII = Matmul2(adjYzII,TzIICYtIIYtII,OnlyDiagonal) 
 adjTdYdadjYdTd = Matmul2(adjTd,YdadjYdTd,OnlyDiagonal) 
 adjTdYsIICYsIITd = Matmul2(adjTd,YsIICYsIITd,OnlyDiagonal) 
 adjTdYzIIadjYzIITd = Matmul2(adjTd,YzIIadjYzIITd,OnlyDiagonal) 
 adjTdTdadjYdYd = Matmul2(adjTd,TdadjYdYd,OnlyDiagonal) 
 adjTdTsIICYsIIYd = Matmul2(adjTd,TsIICYsIIYd,OnlyDiagonal) 
 adjTdTzIIadjYzIIYd = Matmul2(adjTd,TzIIadjYzIIYd,OnlyDiagonal) 
 adjTeYeadjYeTe = Matmul2(adjTe,YeadjYeTe,OnlyDiagonal) 
 adjTeTeadjYeYe = Matmul2(adjTe,TeadjYeYe,OnlyDiagonal) 
 adjTuYuadjYuTu = Matmul2(adjTu,YuadjYuTu,OnlyDiagonal) 
 adjTuTuadjYuYu = Matmul2(adjTu,TuadjYuYu,OnlyDiagonal) 
 adjTzIIYdadjYdTzII = Matmul2(adjTzII,YdadjYdTzII,OnlyDiagonal) 
 adjTzIIYsIICYsIITzII = Matmul2(adjTzII,YsIICYsIITzII,OnlyDiagonal) 
 adjTzIIYzIIadjYzIITzII = Matmul2(adjTzII,YzIIadjYzIITzII,OnlyDiagonal) 
 adjTzIITdadjYdYzII = Matmul2(adjTzII,TdadjYdYzII,OnlyDiagonal) 
 adjTzIITsIICYsIIYzII = Matmul2(adjTzII,TsIICYsIIYzII,OnlyDiagonal) 
 adjTzIITzIIadjYzIIYzII = Matmul2(adjTzII,TzIIadjYzIIYzII,OnlyDiagonal) 
 Cmd2CYdTpYdCYsII = Matmul2(Conjg(md2),CYdTpYdCYsII,OnlyDiagonal) 
 Cmd2CYsIIYsIICYsII = Matmul2(Conjg(md2),CYsIIYsIICYsII,OnlyDiagonal) 
 Cmd2CYsIIYzIIadjYzII = Matmul2(Conjg(md2),CYsIIYzIIadjYzII,OnlyDiagonal) 
 Cmd2CYzIITpYzIICYsII = Matmul2(Conjg(md2),CYzIITpYzIICYsII,OnlyDiagonal) 
 Cml2YtIICYtIIYtII = Matmul2(Conjg(ml2),YtIICYtIIYtII,OnlyDiagonal) 
 Cml2TpYeCYeYtII = Matmul2(Conjg(ml2),TpYeCYeYtII,OnlyDiagonal) 
 Cml2TpYzIICYzIIYtII = Matmul2(Conjg(ml2),TpYzIICYzIIYtII,OnlyDiagonal) 
 CYdCmq2TpYdCYsII = Matmul2(Conjg(Yd),Cmq2TpYdCYsII,OnlyDiagonal) 
 CYdTpYdCmd2CYsII = Matmul2(Conjg(Yd),TpYdCmd2CYsII,OnlyDiagonal) 
 CYdTpYdCYdTpYd = Matmul2(Conjg(Yd),TpYdCYdTpYd,OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYdCYdTpYd(i2,i2) =  Real(CYdTpYdCYdTpYd(i2,i2),dp) 
 CYdTpYdCYdTpTd = Matmul2(Conjg(Yd),TpYdCYdTpTd,OnlyDiagonal) 
 CYdTpYdCYsIImd2 = Matmul2(Conjg(Yd),TpYdCYsIImd2,OnlyDiagonal) 
 CYdTpYdCYsIIYd = Matmul2(Conjg(Yd),TpYdCYsIIYd,OnlyDiagonal) 
 CYdTpYdCYsIIYsII = Matmul2(Conjg(Yd),TpYdCYsIIYsII,OnlyDiagonal) 
 CYdTpYdCYsIIYzII = Matmul2(Conjg(Yd),TpYdCYsIIYzII,OnlyDiagonal) 
 CYdTpYdCYsIITd = Matmul2(Conjg(Yd),TpYdCYsIITd,OnlyDiagonal) 
 CYdTpYdCYsIITsII = Matmul2(Conjg(Yd),TpYdCYsIITsII,OnlyDiagonal) 
 CYdTpYdCYsIITzII = Matmul2(Conjg(Yd),TpYdCYsIITzII,OnlyDiagonal) 
 CYdTpYdCYzIIWOp = Matmul2(Conjg(Yd),TpYdCYzIIWOp,OnlyDiagonal) 
 CYdTpYdCYzIIYtII = Matmul2(Conjg(Yd),TpYdCYzIIYtII,OnlyDiagonal) 
 CYdTpYdCYzIITtII = Matmul2(Conjg(Yd),TpYdCYzIITtII,OnlyDiagonal) 
 CYdTpYuCYuTpYd = Matmul2(Conjg(Yd),TpYuCYuTpYd,OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYuCYuTpYd(i2,i2) =  Real(CYdTpYuCYuTpYd(i2,i2),dp) 
 CYdTpYuCYuTpTd = Matmul2(Conjg(Yd),TpYuCYuTpTd,OnlyDiagonal) 
 CYdTpTdCYdTpYd = Matmul2(Conjg(Yd),TpTdCYdTpYd,OnlyDiagonal) 
 CYdTpTdCYsIIYd = Matmul2(Conjg(Yd),TpTdCYsIIYd,OnlyDiagonal) 
 CYdTpTdCYsIIYsII = Matmul2(Conjg(Yd),TpTdCYsIIYsII,OnlyDiagonal) 
 CYdTpTdCYsIIYzII = Matmul2(Conjg(Yd),TpTdCYsIIYzII,OnlyDiagonal) 
 CYdTpTdCYzIIYtII = Matmul2(Conjg(Yd),TpTdCYzIIYtII,OnlyDiagonal) 
 CYdTpTuCYuTpYd = Matmul2(Conjg(Yd),TpTuCYuTpYd,OnlyDiagonal) 
 CYeTpYeCYeWOp = Matmul2(Conjg(Ye),TpYeCYeWOp,OnlyDiagonal) 
 CYeTpYeCYeYtII = Matmul2(Conjg(Ye),TpYeCYeYtII,OnlyDiagonal) 
 CYeTpYeCYeTtII = Matmul2(Conjg(Ye),TpYeCYeTtII,OnlyDiagonal) 
 CYeTpTeCYeYtII = Matmul2(Conjg(Ye),TpTeCYeYtII,OnlyDiagonal) 
 CYsIImd2YsIICYsII = Matmul2(adjYsII,md2YsIICYsII,OnlyDiagonal) 
 CYsIIYdadjYdYsII = Matmul2(adjYsII,YdadjYdYsII,OnlyDiagonal) 
 CYsIIYdadjYdTsII = Matmul2(adjYsII,YdadjYdTsII,OnlyDiagonal) 
 CYsIIYsIICmd2CYsII = Matmul2(adjYsII,YsIICmd2CYsII,OnlyDiagonal) 
 CYsIIYsIICYsIImd2 = Matmul2(adjYsII,YsIICYsIImd2,OnlyDiagonal) 
 CYsIIYsIICYsIIYd = Matmul2(adjYsII,YsIICYsIIYd,OnlyDiagonal) 
 CYsIIYsIICYsIIYsII = Matmul2(adjYsII,YsIICYsIIYsII,OnlyDiagonal) 
 CYsIIYsIICYsIIYzII = Matmul2(adjYsII,YsIICYsIIYzII,OnlyDiagonal) 
 CYsIIYsIICYsIITd = Matmul2(adjYsII,YsIICYsIITd,OnlyDiagonal) 
 CYsIIYsIICYsIITsII = Matmul2(adjYsII,YsIICYsIITsII,OnlyDiagonal) 
 CYsIIYsIICYsIITzII = Matmul2(adjYsII,YsIICYsIITzII,OnlyDiagonal) 
 CYsIIYsIICYzIIWOp = Matmul2(adjYsII,YsIICYzIIWOp,OnlyDiagonal) 
 CYsIIYsIICYzIIYtII = Matmul2(adjYsII,YsIICYzIIYtII,OnlyDiagonal) 
 CYsIIYsIICYzIITtII = Matmul2(adjYsII,YsIICYzIITtII,OnlyDiagonal) 
 CYsIIYzIIadjYzIIYsII = Matmul2(adjYsII,YzIIadjYzIIYsII,OnlyDiagonal) 
 CYsIIYzIIadjYzIITsII = Matmul2(adjYsII,YzIIadjYzIITsII,OnlyDiagonal) 
 CYsIITdadjYdYsII = Matmul2(adjYsII,TdadjYdYsII,OnlyDiagonal) 
 CYsIITsIICYsIIYd = Matmul2(adjYsII,TsIICYsIIYd,OnlyDiagonal) 
 CYsIITsIICYsIIYsII = Matmul2(adjYsII,TsIICYsIIYsII,OnlyDiagonal) 
 CYsIITsIICYsIIYzII = Matmul2(adjYsII,TsIICYsIIYzII,OnlyDiagonal) 
 CYsIITsIICYzIIYtII = Matmul2(adjYsII,TsIICYzIIYtII,OnlyDiagonal) 
 CYsIITsIICTdTpYd = Matmul2(adjYsII,TsIICTdTpYd,OnlyDiagonal) 
 CYsIITsIICTzIITpYzII = Matmul2(adjYsII,TsIICTzIITpYzII,OnlyDiagonal) 
 CYsIITzIIadjYzIIYsII = Matmul2(adjYsII,TzIIadjYzIIYsII,OnlyDiagonal) 
 CYtIIYtIIml2adjYe = Matmul2(adjYtII,YtIIml2adjYe,OnlyDiagonal) 
 CYtIIYtIIml2adjYzII = Matmul2(adjYtII,YtIIml2adjYzII,OnlyDiagonal) 
 CYtIIYtIIadjYeme2 = Matmul2(adjYtII,YtIIadjYeme2,OnlyDiagonal) 
 CYtIIYtIIadjYeYe = Matmul2(adjYtII,YtIIadjYeYe,OnlyDiagonal) 
 CYtIIYtIIadjYeTe = Matmul2(adjYtII,YtIIadjYeTe,OnlyDiagonal) 
 CYtIIYtIIadjYzIImd2 = Matmul2(adjYtII,YtIIadjYzIImd2,OnlyDiagonal) 
 CYtIIYtIIadjYzIIYd = Matmul2(adjYtII,YtIIadjYzIIYd,OnlyDiagonal) 
 CYtIIYtIIadjYzIIYsII = Matmul2(adjYtII,YtIIadjYzIIYsII,OnlyDiagonal) 
 CYtIIYtIIadjYzIIYzII = Matmul2(adjYtII,YtIIadjYzIIYzII,OnlyDiagonal) 
 CYtIIYtIIadjYzIITd = Matmul2(adjYtII,YtIIadjYzIITd,OnlyDiagonal) 
 CYtIIYtIIadjYzIITsII = Matmul2(adjYtII,YtIIadjYzIITsII,OnlyDiagonal) 
 CYtIIYtIIadjYzIITzII = Matmul2(adjYtII,YtIIadjYzIITzII,OnlyDiagonal) 
 CYtIIYtIICYtIIWOp = Matmul2(adjYtII,YtIICYtIIWOp,OnlyDiagonal) 
 CYtIIYtIICYtIIYtII = Matmul2(adjYtII,YtIICYtIIYtII,OnlyDiagonal) 
 CYtIIYtIICYtIITtII = Matmul2(adjYtII,YtIICYtIITtII,OnlyDiagonal) 
 CYtIIYtIICTtIITtII = Matmul2(adjYtII,YtIICTtIITtII,OnlyDiagonal) 
 CYtIICml2YtIIadjYe = Matmul2(adjYtII,Cml2YtIIadjYe,OnlyDiagonal) 
 CYtIICml2YtIIadjYzII = Matmul2(adjYtII,Cml2YtIIadjYzII,OnlyDiagonal) 
 CYtIITtIIadjYeYe = Matmul2(adjYtII,TtIIadjYeYe,OnlyDiagonal) 
 CYtIITtIIadjYzIIYd = Matmul2(adjYtII,TtIIadjYzIIYd,OnlyDiagonal) 
 CYtIITtIIadjYzIIYsII = Matmul2(adjYtII,TtIIadjYzIIYsII,OnlyDiagonal) 
 CYtIITtIIadjYzIIYzII = Matmul2(adjYtII,TtIIadjYzIIYzII,OnlyDiagonal) 
 CYtIITtIICYtIIYtII = Matmul2(adjYtII,TtIICYtIIYtII,OnlyDiagonal) 
 CYtIITtIICTtIIYtII = Matmul2(adjYtII,TtIICTtIIYtII,OnlyDiagonal) 
 CYtIITpYeCYeYtII = Matmul2(adjYtII,TpYeCYeYtII,OnlyDiagonal) 
 CYtIITpYeCYeTtII = Matmul2(adjYtII,TpYeCYeTtII,OnlyDiagonal) 
 CYtIITpYeCTeTtII = Matmul2(adjYtII,TpYeCTeTtII,OnlyDiagonal) 
 CYtIITpYzIICYzIIYtII = Matmul2(adjYtII,TpYzIICYzIIYtII,OnlyDiagonal) 
 CYtIITpYzIICYzIITtII = Matmul2(adjYtII,TpYzIICYzIITtII,OnlyDiagonal) 
 CYtIITpYzIICTzIITtII = Matmul2(adjYtII,TpYzIICTzIITtII,OnlyDiagonal) 
 CYtIITpTeCYeYtII = Matmul2(adjYtII,TpTeCYeYtII,OnlyDiagonal) 
 CYtIITpTeCTeYtII = Matmul2(adjYtII,TpTeCTeYtII,OnlyDiagonal) 
 CYtIITpTzIICYzIIYtII = Matmul2(adjYtII,TpTzIICYzIIYtII,OnlyDiagonal) 
 CYtIITpTzIICTzIIYtII = Matmul2(adjYtII,TpTzIICTzIIYtII,OnlyDiagonal) 
 CYzIIYtIICYtIITpYzII = Matmul2(Conjg(YzII),YtIICYtIITpYzII,OnlyDiagonal) 
 CYzIIYtIICYtIITpTzII = Matmul2(Conjg(YzII),YtIICYtIITpTzII,OnlyDiagonal) 
 CYzIICml2TpYzIICYsII = Matmul2(Conjg(YzII),Cml2TpYzIICYsII,OnlyDiagonal) 
 CYzIITtIICYtIITpYzII = Matmul2(Conjg(YzII),TtIICYtIITpYzII,OnlyDiagonal) 
 CYzIITpYeCYeTpYzII = Matmul2(Conjg(YzII),TpYeCYeTpYzII,OnlyDiagonal) 
Forall(i2=1:3)  CYzIITpYeCYeTpYzII(i2,i2) =  Real(CYzIITpYeCYeTpYzII(i2,i2),dp) 
 CYzIITpYeCYeTpTzII = Matmul2(Conjg(YzII),TpYeCYeTpTzII,OnlyDiagonal) 
 CYzIITpYzIICmd2CYsII = Matmul2(Conjg(YzII),TpYzIICmd2CYsII,OnlyDiagonal) 
 CYzIITpYzIICYsIImd2 = Matmul2(Conjg(YzII),TpYzIICYsIImd2,OnlyDiagonal) 
 CYzIITpYzIICYsIIYd = Matmul2(Conjg(YzII),TpYzIICYsIIYd,OnlyDiagonal) 
 CYzIITpYzIICYsIIYsII = Matmul2(Conjg(YzII),TpYzIICYsIIYsII,OnlyDiagonal) 
 CYzIITpYzIICYsIIYzII = Matmul2(Conjg(YzII),TpYzIICYsIIYzII,OnlyDiagonal) 
 CYzIITpYzIICYsIITd = Matmul2(Conjg(YzII),TpYzIICYsIITd,OnlyDiagonal) 
 CYzIITpYzIICYsIITsII = Matmul2(Conjg(YzII),TpYzIICYsIITsII,OnlyDiagonal) 
 CYzIITpYzIICYsIITzII = Matmul2(Conjg(YzII),TpYzIICYsIITzII,OnlyDiagonal) 
 CYzIITpYzIICYzIIWOp = Matmul2(Conjg(YzII),TpYzIICYzIIWOp,OnlyDiagonal) 
 CYzIITpYzIICYzIIYtII = Matmul2(Conjg(YzII),TpYzIICYzIIYtII,OnlyDiagonal) 
 CYzIITpYzIICYzIITtII = Matmul2(Conjg(YzII),TpYzIICYzIITtII,OnlyDiagonal) 
 CYzIITpYzIICYzIITpYzII = Matmul2(Conjg(YzII),TpYzIICYzIITpYzII,OnlyDiagonal) 
Forall(i2=1:3)  CYzIITpYzIICYzIITpYzII(i2,i2) =  Real(CYzIITpYzIICYzIITpYzII(i2,i2),dp) 
 CYzIITpYzIICYzIITpTzII = Matmul2(Conjg(YzII),TpYzIICYzIITpTzII,OnlyDiagonal) 
 CYzIITpTeCYeTpYzII = Matmul2(Conjg(YzII),TpTeCYeTpYzII,OnlyDiagonal) 
 CYzIITpTzIICYsIIYd = Matmul2(Conjg(YzII),TpTzIICYsIIYd,OnlyDiagonal) 
 CYzIITpTzIICYsIIYsII = Matmul2(Conjg(YzII),TpTzIICYsIIYsII,OnlyDiagonal) 
 CYzIITpTzIICYsIIYzII = Matmul2(Conjg(YzII),TpTzIICYsIIYzII,OnlyDiagonal) 
 CYzIITpTzIICYzIIYtII = Matmul2(Conjg(YzII),TpTzIICYzIIYtII,OnlyDiagonal) 
 CYzIITpTzIICYzIITpYzII = Matmul2(Conjg(YzII),TpTzIICYzIITpYzII,OnlyDiagonal) 
 CTtIIYtIICYtIITtII = Matmul2(adjTtII,YtIICYtIITtII,OnlyDiagonal) 
 CTtIITtIICYtIIYtII = Matmul2(adjTtII,TtIICYtIIYtII,OnlyDiagonal) 
 CTtIITpYeCYeTtII = Matmul2(adjTtII,TpYeCYeTtII,OnlyDiagonal) 
 CTtIITpYzIICYzIITtII = Matmul2(adjTtII,TpYzIICYzIITtII,OnlyDiagonal) 
 CTtIITpTeCYeYtII = Matmul2(adjTtII,TpTeCYeYtII,OnlyDiagonal) 
 CTtIITpTzIICYzIIYtII = Matmul2(adjTtII,TpTzIICYzIIYtII,OnlyDiagonal) 
 TdadjYdYdadjTd = Matmul2(Td,adjYdYdadjTd,OnlyDiagonal) 
 TdadjYdYsIICTsII = Matmul2(Td,adjYdYsIICTsII,OnlyDiagonal) 
 TdadjYdCYsIIYd = Matmul2(Td,adjYdCYsIIYd,OnlyDiagonal) 
 TdadjYdCYsIIYzII = Matmul2(Td,adjYdCYsIIYzII,OnlyDiagonal) 
 TdadjYuYuadjTd = Matmul2(Td,adjYuYuadjTd,OnlyDiagonal) 
 TdadjTdYdadjYd = Matmul2(Td,adjTdYdadjYd,OnlyDiagonal) 
 TdadjTuYuadjYd = Matmul2(Td,adjTuYuadjYd,OnlyDiagonal) 
 TeadjYeYeadjTe = Matmul2(Te,adjYeYeadjTe,OnlyDiagonal) 
 TeadjYzIIYzIIadjTe = Matmul2(Te,adjYzIIYzIIadjTe,OnlyDiagonal) 
 TeadjTeYeadjYe = Matmul2(Te,adjTeYeadjYe,OnlyDiagonal) 
 TeadjTzIIYzIIadjYe = Matmul2(Te,adjTzIIYzIIadjYe,OnlyDiagonal) 
 TeCYtIIYtIIadjTe = Matmul2(Te,CYtIIYtIIadjTe,OnlyDiagonal) 
 TeCTtIIYtIIadjYe = Matmul2(Te,CTtIIYtIIadjYe,OnlyDiagonal) 
 TsIIYdadjTdCYsII = Matmul2(TsII,YdadjTdCYsII,OnlyDiagonal) 
 TsIIYzIIadjTzIICYsII = Matmul2(TsII,YzIIadjTzIICYsII,OnlyDiagonal) 
 TsIICYdTpYdCTsII = Matmul2(TsII,CYdTpYdCTsII,OnlyDiagonal) 
 TsIICYsIIYsIICTsII = Matmul2(TsII,CYsIIYsIICTsII,OnlyDiagonal) 
 TsIICYzIITpYzIICTsII = Matmul2(TsII,CYzIITpYzIICTsII,OnlyDiagonal) 
 TsIICTdTpYdCYsII = Matmul2(TsII,CTdTpYdCYsII,OnlyDiagonal) 
 TsIICTsIIYsIICYsII = Matmul2(TsII,CTsIIYsIICYsII,OnlyDiagonal) 
 TsIICTzIITpYzIICYsII = Matmul2(TsII,CTzIITpYzIICYsII,OnlyDiagonal) 
 TuadjYdYdadjTu = Matmul2(Tu,adjYdYdadjTu,OnlyDiagonal) 
 TuadjYuYuadjTu = Matmul2(Tu,adjYuYuadjTu,OnlyDiagonal) 
 TuadjTdYdadjYu = Matmul2(Tu,adjTdYdadjYu,OnlyDiagonal) 
 TuadjTuYuadjYu = Matmul2(Tu,adjTuYuadjYu,OnlyDiagonal) 
 TzIIadjYeYeadjTzII = Matmul2(TzII,adjYeYeadjTzII,OnlyDiagonal) 
 TzIIadjYzIIYsIICTsII = Matmul2(TzII,adjYzIIYsIICTsII,OnlyDiagonal) 
 TzIIadjYzIIYzIIadjTzII = Matmul2(TzII,adjYzIIYzIIadjTzII,OnlyDiagonal) 
 TzIIadjYzIICYsIIYd = Matmul2(TzII,adjYzIICYsIIYd,OnlyDiagonal) 
 TzIIadjYzIICYsIIYzII = Matmul2(TzII,adjYzIICYsIIYzII,OnlyDiagonal) 
 TzIIadjTeYeadjYzII = Matmul2(TzII,adjTeYeadjYzII,OnlyDiagonal) 
 TzIIadjTzIIYzIIadjYzII = Matmul2(TzII,adjTzIIYzIIadjYzII,OnlyDiagonal) 
 TzIICYtIIYtIIadjTzII = Matmul2(TzII,CYtIIYtIIadjTzII,OnlyDiagonal) 
 TzIICTtIIYtIIadjYzII = Matmul2(TzII,CTtIIYtIIadjYzII,OnlyDiagonal) 
 TpYeCme2CYeYtII = Matmul2(Transpose(Ye),Cme2CYeYtII,OnlyDiagonal) 
 TpYeCYeYtIIml2 = Matmul2(Transpose(Ye),CYeYtIIml2,OnlyDiagonal) 
 TpYeCYeCml2YtII = Matmul2(Transpose(Ye),CYeCml2YtII,OnlyDiagonal) 
 TpYzIICmd2CYzIIYtII = Matmul2(Transpose(YzII),Cmd2CYzIIYtII,OnlyDiagonal) 
 TpYzIICYzIIYtIIml2 = Matmul2(Transpose(YzII),CYzIIYtIIml2,OnlyDiagonal) 
 TpYzIICYzIICml2YtII = Matmul2(Transpose(YzII),CYzIICml2YtII,OnlyDiagonal) 
 md2YdadjYdYdadjYd = Matmul2(md2,YdadjYdYdadjYd,OnlyDiagonal) 
 md2YdadjYdYsIICYsII = Matmul2(md2,YdadjYdYsIICYsII,OnlyDiagonal) 
 md2YdadjYdYzIIadjYzII = Matmul2(md2,YdadjYdYzIIadjYzII,OnlyDiagonal) 
 md2YdadjYuYuadjYd = Matmul2(md2,YdadjYuYuadjYd,OnlyDiagonal) 
 md2YsIICYdTpYdCYsII = Matmul2(md2,YsIICYdTpYdCYsII,OnlyDiagonal) 
 md2YsIICYsIIYdadjYd = Matmul2(md2,YsIICYsIIYdadjYd,OnlyDiagonal) 
 md2YsIICYsIIYsIICYsII = Matmul2(md2,YsIICYsIIYsIICYsII,OnlyDiagonal) 
 md2YsIICYsIIYzIIadjYzII = Matmul2(md2,YsIICYsIIYzIIadjYzII,OnlyDiagonal) 
 md2YsIICYzIITpYzIICYsII = Matmul2(md2,YsIICYzIITpYzIICYsII,OnlyDiagonal) 
 md2YzIIadjYeYeadjYzII = Matmul2(md2,YzIIadjYeYeadjYzII,OnlyDiagonal) 
 md2YzIIadjYzIIYdadjYd = Matmul2(md2,YzIIadjYzIIYdadjYd,OnlyDiagonal) 
 md2YzIIadjYzIIYsIICYsII = Matmul2(md2,YzIIadjYzIIYsIICYsII,OnlyDiagonal) 
 md2YzIIadjYzIIYzIIadjYzII = Matmul2(md2,YzIIadjYzIIYzIIadjYzII,OnlyDiagonal) 
 md2YzIICYtIIYtIIadjYzII = Matmul2(md2,YzIICYtIIYtIIadjYzII,OnlyDiagonal) 
 me2YeadjYeYeadjYe = Matmul2(me2,YeadjYeYeadjYe,OnlyDiagonal) 
 me2YeadjYzIIYzIIadjYe = Matmul2(me2,YeadjYzIIYzIIadjYe,OnlyDiagonal) 
 me2YeCYtIIYtIIadjYe = Matmul2(me2,YeCYtIIYtIIadjYe,OnlyDiagonal) 
 ml2adjYeYeadjYeYe = Matmul2(ml2,adjYeYeadjYeYe,OnlyDiagonal) 
 ml2adjYeYeadjYzIIYzII = Matmul2(ml2,adjYeYeadjYzIIYzII,OnlyDiagonal) 
 ml2adjYeYeCYtIIYtII = Matmul2(ml2,adjYeYeCYtIIYtII,OnlyDiagonal) 
 ml2adjYzIIYdadjYdYzII = Matmul2(ml2,adjYzIIYdadjYdYzII,OnlyDiagonal) 
 ml2adjYzIIYsIICYsIIYzII = Matmul2(ml2,adjYzIIYsIICYsIIYzII,OnlyDiagonal) 
 ml2adjYzIIYzIIadjYeYe = Matmul2(ml2,adjYzIIYzIIadjYeYe,OnlyDiagonal) 
 ml2adjYzIIYzIIadjYzIIYzII = Matmul2(ml2,adjYzIIYzIIadjYzIIYzII,OnlyDiagonal) 
 ml2adjYzIIYzIICYtIIYtII = Matmul2(ml2,adjYzIIYzIICYtIIYtII,OnlyDiagonal) 
 ml2CYtIIYtIIadjYeYe = Matmul2(ml2,CYtIIYtIIadjYeYe,OnlyDiagonal) 
 ml2CYtIIYtIIadjYzIIYzII = Matmul2(ml2,CYtIIYtIIadjYzIIYzII,OnlyDiagonal) 
 ml2CYtIIYtIICYtIIYtII = Matmul2(ml2,CYtIIYtIICYtIIYtII,OnlyDiagonal) 
 ml2CYtIITpYeCYeYtII = Matmul2(ml2,CYtIITpYeCYeYtII,OnlyDiagonal) 
 ml2CYtIITpYzIICYzIIYtII = Matmul2(ml2,CYtIITpYzIICYzIIYtII,OnlyDiagonal) 
 mq2adjYdYdadjYdYd = Matmul2(mq2,adjYdYdadjYdYd,OnlyDiagonal) 
 mq2adjYdYdadjYuYu = Matmul2(mq2,adjYdYdadjYuYu,OnlyDiagonal) 
 mq2adjYdYsIICYsIIYd = Matmul2(mq2,adjYdYsIICYsIIYd,OnlyDiagonal) 
 mq2adjYdYzIIadjYzIIYd = Matmul2(mq2,adjYdYzIIadjYzIIYd,OnlyDiagonal) 
 mq2adjYuYuadjYdYd = Matmul2(mq2,adjYuYuadjYdYd,OnlyDiagonal) 
 mq2adjYuYuadjYuYu = Matmul2(mq2,adjYuYuadjYuYu,OnlyDiagonal) 
 mu2YuadjYdYdadjYu = Matmul2(mu2,YuadjYdYdadjYu,OnlyDiagonal) 
 mu2YuadjYuYuadjYu = Matmul2(mu2,YuadjYuYuadjYu,OnlyDiagonal) 
 WOpadjYeYeadjYeYe = Matmul2(WOp,adjYeYeadjYeYe,OnlyDiagonal) 
 WOpadjYzIIYdadjYdYzII = Matmul2(WOp,adjYzIIYdadjYdYzII,OnlyDiagonal) 
 WOpadjYzIIYsIICYsIIYzII = Matmul2(WOp,adjYzIIYsIICYsIIYzII,OnlyDiagonal) 
 WOpadjYzIIYzIIadjYzIIYzII = Matmul2(WOp,adjYzIIYzIIadjYzIIYzII,OnlyDiagonal) 
 WOpCYtIIYtIICYtIIYtII = Matmul2(WOp,CYtIIYtIICYtIIYtII,OnlyDiagonal) 
 WOpCYtIITpYeCYeYtII = Matmul2(WOp,CYtIITpYeCYeYtII,OnlyDiagonal) 
 WOpCYtIITpYzIICYzIIYtII = Matmul2(WOp,CYtIITpYzIICYzIIYtII,OnlyDiagonal) 
 Ydmq2adjYdYdadjYd = Matmul2(Yd,mq2adjYdYdadjYd,OnlyDiagonal) 
 Ydmq2adjYuYuadjYd = Matmul2(Yd,mq2adjYuYuadjYd,OnlyDiagonal) 
 YdadjYdmd2YdadjYd = Matmul2(Yd,adjYdmd2YdadjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYdmd2YdadjYd(i2,i2) =  Real(YdadjYdmd2YdadjYd(i2,i2),dp) 
 YdadjYdYdmq2adjYd = Matmul2(Yd,adjYdYdmq2adjYd,OnlyDiagonal) 
 YdadjYdYdadjYdmd2 = Matmul2(Yd,adjYdYdadjYdmd2,OnlyDiagonal) 
 YdadjYdYdadjYdYd = Matmul2(Yd,adjYdYdadjYdYd,OnlyDiagonal) 
 YdadjYdYdadjYdYsII = Matmul2(Yd,adjYdYdadjYdYsII,OnlyDiagonal) 
 YdadjYdYdadjYdYzII = Matmul2(Yd,adjYdYdadjYdYzII,OnlyDiagonal) 
 YdadjYdYdadjYdTd = Matmul2(Yd,adjYdYdadjYdTd,OnlyDiagonal) 
 YdadjYdYdadjYdTsII = Matmul2(Yd,adjYdYdadjYdTsII,OnlyDiagonal) 
 YdadjYdYdadjYdTzII = Matmul2(Yd,adjYdYdadjYdTzII,OnlyDiagonal) 
 YdadjYdYsIICmd2CYsII = Matmul2(Yd,adjYdYsIICmd2CYsII,OnlyDiagonal) 
 YdadjYdYsIICYsIIYd = Matmul2(Yd,adjYdYsIICYsIIYd,OnlyDiagonal) 
 YdadjYdYsIICYsIITd = Matmul2(Yd,adjYdYsIICYsIITd,OnlyDiagonal) 
 YdadjYdYzIIadjYzIIYd = Matmul2(Yd,adjYdYzIIadjYzIIYd,OnlyDiagonal) 
 YdadjYdYzIIadjYzIITd = Matmul2(Yd,adjYdYzIIadjYzIITd,OnlyDiagonal) 
 YdadjYdTdadjYdYd = Matmul2(Yd,adjYdTdadjYdYd,OnlyDiagonal) 
 YdadjYdTdadjYdYsII = Matmul2(Yd,adjYdTdadjYdYsII,OnlyDiagonal) 
 YdadjYdTdadjYdYzII = Matmul2(Yd,adjYdTdadjYdYzII,OnlyDiagonal) 
 YdadjYdTsIICYsIIYd = Matmul2(Yd,adjYdTsIICYsIIYd,OnlyDiagonal) 
 YdadjYdTzIIadjYzIIYd = Matmul2(Yd,adjYdTzIIadjYzIIYd,OnlyDiagonal) 
 YdadjYumu2YuadjYd = Matmul2(Yd,adjYumu2YuadjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYumu2YuadjYd(i2,i2) =  Real(YdadjYumu2YuadjYd(i2,i2),dp) 
 YdadjYuYumq2adjYd = Matmul2(Yd,adjYuYumq2adjYd,OnlyDiagonal) 
 YdadjYuYuadjYdmd2 = Matmul2(Yd,adjYuYuadjYdmd2,OnlyDiagonal) 
 YdadjYuYuadjYdYd = Matmul2(Yd,adjYuYuadjYdYd,OnlyDiagonal) 
 YdadjYuYuadjYdYsII = Matmul2(Yd,adjYuYuadjYdYsII,OnlyDiagonal) 
 YdadjYuYuadjYdYzII = Matmul2(Yd,adjYuYuadjYdYzII,OnlyDiagonal) 
 YdadjYuYuadjYdTd = Matmul2(Yd,adjYuYuadjYdTd,OnlyDiagonal) 
 YdadjYuYuadjYdTsII = Matmul2(Yd,adjYuYuadjYdTsII,OnlyDiagonal) 
 YdadjYuYuadjYdTzII = Matmul2(Yd,adjYuYuadjYdTzII,OnlyDiagonal) 
 YdadjYuYuadjYuYu = Matmul2(Yd,adjYuYuadjYuYu,OnlyDiagonal) 
 YdadjYuYuadjYuTu = Matmul2(Yd,adjYuYuadjYuTu,OnlyDiagonal) 
 YdadjYuTuadjYdYd = Matmul2(Yd,adjYuTuadjYdYd,OnlyDiagonal) 
 YdadjYuTuadjYdYsII = Matmul2(Yd,adjYuTuadjYdYsII,OnlyDiagonal) 
 YdadjYuTuadjYdYzII = Matmul2(Yd,adjYuTuadjYdYzII,OnlyDiagonal) 
 YdadjYuTuadjYuYu = Matmul2(Yd,adjYuTuadjYuYu,OnlyDiagonal) 
 Yeml2adjYeYeadjYe = Matmul2(Ye,ml2adjYeYeadjYe,OnlyDiagonal) 
 Yeml2adjYzIIYzIIadjYe = Matmul2(Ye,ml2adjYzIIYzIIadjYe,OnlyDiagonal) 
 Yeml2CYtIIYtIIadjYe = Matmul2(Ye,ml2CYtIIYtIIadjYe,OnlyDiagonal) 
 YeadjYeme2YeadjYe = Matmul2(Ye,adjYeme2YeadjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYeme2YeadjYe(i2,i2) =  Real(YeadjYeme2YeadjYe(i2,i2),dp) 
 YeadjYeYeml2adjYe = Matmul2(Ye,adjYeYeml2adjYe,OnlyDiagonal) 
 YeadjYeYeadjYeme2 = Matmul2(Ye,adjYeYeadjYeme2,OnlyDiagonal) 
 YeadjYeYeadjYeYe = Matmul2(Ye,adjYeYeadjYeYe,OnlyDiagonal) 
 YeadjYeYeadjYeTe = Matmul2(Ye,adjYeYeadjYeTe,OnlyDiagonal) 
 YeadjYeTeadjYeYe = Matmul2(Ye,adjYeTeadjYeYe,OnlyDiagonal) 
 YeadjYzIImd2YzIIadjYe = Matmul2(Ye,adjYzIImd2YzIIadjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYzIImd2YzIIadjYe(i2,i2) =  Real(YeadjYzIImd2YzIIadjYe(i2,i2),dp) 
 YeadjYzIIYdadjYdYzII = Matmul2(Ye,adjYzIIYdadjYdYzII,OnlyDiagonal) 
 YeadjYzIIYdadjYdTzII = Matmul2(Ye,adjYzIIYdadjYdTzII,OnlyDiagonal) 
 YeadjYzIIYsIICYsIIYzII = Matmul2(Ye,adjYzIIYsIICYsIIYzII,OnlyDiagonal) 
 YeadjYzIIYsIICYsIITzII = Matmul2(Ye,adjYzIIYsIICYsIITzII,OnlyDiagonal) 
 YeadjYzIIYzIIml2adjYe = Matmul2(Ye,adjYzIIYzIIml2adjYe,OnlyDiagonal) 
 YeadjYzIIYzIIadjYeme2 = Matmul2(Ye,adjYzIIYzIIadjYeme2,OnlyDiagonal) 
 YeadjYzIIYzIIadjYeYe = Matmul2(Ye,adjYzIIYzIIadjYeYe,OnlyDiagonal) 
 YeadjYzIIYzIIadjYeTe = Matmul2(Ye,adjYzIIYzIIadjYeTe,OnlyDiagonal) 
 YeadjYzIIYzIIadjYzIIYzII = Matmul2(Ye,adjYzIIYzIIadjYzIIYzII,OnlyDiagonal) 
 YeadjYzIIYzIIadjYzIITzII = Matmul2(Ye,adjYzIIYzIIadjYzIITzII,OnlyDiagonal) 
 YeadjYzIITdadjYdYzII = Matmul2(Ye,adjYzIITdadjYdYzII,OnlyDiagonal) 
 YeadjYzIITsIICYsIIYzII = Matmul2(Ye,adjYzIITsIICYsIIYzII,OnlyDiagonal) 
 YeadjYzIITzIIadjYeYe = Matmul2(Ye,adjYzIITzIIadjYeYe,OnlyDiagonal) 
 YeadjYzIITzIIadjYzIIYzII = Matmul2(Ye,adjYzIITzIIadjYzIIYzII,OnlyDiagonal) 
 YeCYtIIYtIIml2adjYe = Matmul2(Ye,CYtIIYtIIml2adjYe,OnlyDiagonal) 
 YeCYtIIYtIIadjYeme2 = Matmul2(Ye,CYtIIYtIIadjYeme2,OnlyDiagonal) 
 YeCYtIIYtIIadjYeYe = Matmul2(Ye,CYtIIYtIIadjYeYe,OnlyDiagonal) 
 YeCYtIIYtIIadjYeTe = Matmul2(Ye,CYtIIYtIIadjYeTe,OnlyDiagonal) 
 YeCYtIIYtIICYtIIYtII = Matmul2(Ye,CYtIIYtIICYtIIYtII,OnlyDiagonal) 
 YeCYtIIYtIICYtIITtII = Matmul2(Ye,CYtIIYtIICYtIITtII,OnlyDiagonal) 
 YeCYtIICml2YtIIadjYe = Matmul2(Ye,CYtIICml2YtIIadjYe,OnlyDiagonal) 
 YeCYtIITtIIadjYeYe = Matmul2(Ye,CYtIITtIIadjYeYe,OnlyDiagonal) 
 YeCYtIITtIICYtIIYtII = Matmul2(Ye,CYtIITtIICYtIIYtII,OnlyDiagonal) 
 YeCYtIITpYeCYeYtII = Matmul2(Ye,CYtIITpYeCYeYtII,OnlyDiagonal) 
 YeCYtIITpYeCYeTtII = Matmul2(Ye,CYtIITpYeCYeTtII,OnlyDiagonal) 
 YeCYtIITpYzIICYzIIYtII = Matmul2(Ye,CYtIITpYzIICYzIIYtII,OnlyDiagonal) 
 YeCYtIITpYzIICYzIITtII = Matmul2(Ye,CYtIITpYzIICYzIITtII,OnlyDiagonal) 
 YeCYtIITpTeCYeYtII = Matmul2(Ye,CYtIITpTeCYeYtII,OnlyDiagonal) 
 YeCYtIITpTzIICYzIIYtII = Matmul2(Ye,CYtIITpTzIICYzIIYtII,OnlyDiagonal) 
 YsIICmd2CYdTpYdCYsII = Matmul2(YsII,Cmd2CYdTpYdCYsII,OnlyDiagonal) 
 YsIICmd2CYsIIYsIICYsII = Matmul2(YsII,Cmd2CYsIIYsIICYsII,OnlyDiagonal) 
 YsIICmd2CYsIIYzIIadjYzII = Matmul2(YsII,Cmd2CYsIIYzIIadjYzII,OnlyDiagonal) 
 YsIICmd2CYzIITpYzIICYsII = Matmul2(YsII,Cmd2CYzIITpYzIICYsII,OnlyDiagonal) 
 YsIICYdCmq2TpYdCYsII = Matmul2(YsII,CYdCmq2TpYdCYsII,OnlyDiagonal) 
 YsIICYdTpYdCmd2CYsII = Matmul2(YsII,CYdTpYdCmd2CYsII,OnlyDiagonal) 
 YsIICYdTpYdCYdTpYd = Matmul2(YsII,CYdTpYdCYdTpYd,OnlyDiagonal) 
 YsIICYdTpYdCYdTpTd = Matmul2(YsII,CYdTpYdCYdTpTd,OnlyDiagonal) 
 YsIICYdTpYdCYsIImd2 = Matmul2(YsII,CYdTpYdCYsIImd2,OnlyDiagonal) 
 YsIICYdTpYdCYsIIYd = Matmul2(YsII,CYdTpYdCYsIIYd,OnlyDiagonal) 
 YsIICYdTpYdCYsIIYsII = Matmul2(YsII,CYdTpYdCYsIIYsII,OnlyDiagonal) 
 YsIICYdTpYdCYsIIYzII = Matmul2(YsII,CYdTpYdCYsIIYzII,OnlyDiagonal) 
 YsIICYdTpYdCYsIITd = Matmul2(YsII,CYdTpYdCYsIITd,OnlyDiagonal) 
 YsIICYdTpYdCYsIITsII = Matmul2(YsII,CYdTpYdCYsIITsII,OnlyDiagonal) 
 YsIICYdTpYdCYsIITzII = Matmul2(YsII,CYdTpYdCYsIITzII,OnlyDiagonal) 
 YsIICYdTpYuCYuTpYd = Matmul2(YsII,CYdTpYuCYuTpYd,OnlyDiagonal) 
 YsIICYdTpYuCYuTpTd = Matmul2(YsII,CYdTpYuCYuTpTd,OnlyDiagonal) 
 YsIICYdTpTdCYdTpYd = Matmul2(YsII,CYdTpTdCYdTpYd,OnlyDiagonal) 
 YsIICYdTpTdCYsIIYd = Matmul2(YsII,CYdTpTdCYsIIYd,OnlyDiagonal) 
 YsIICYdTpTdCYsIIYsII = Matmul2(YsII,CYdTpTdCYsIIYsII,OnlyDiagonal) 
 YsIICYdTpTdCYsIIYzII = Matmul2(YsII,CYdTpTdCYsIIYzII,OnlyDiagonal) 
 YsIICYdTpTuCYuTpYd = Matmul2(YsII,CYdTpTuCYuTpYd,OnlyDiagonal) 
 YsIICYsIImd2YsIICYsII = Matmul2(YsII,CYsIImd2YsIICYsII,OnlyDiagonal) 
 YsIICYsIIYdadjYdYsII = Matmul2(YsII,CYsIIYdadjYdYsII,OnlyDiagonal) 
 YsIICYsIIYdadjYdTsII = Matmul2(YsII,CYsIIYdadjYdTsII,OnlyDiagonal) 
 YsIICYsIIYsIICmd2CYsII = Matmul2(YsII,CYsIIYsIICmd2CYsII,OnlyDiagonal) 
 YsIICYsIIYsIICYsIImd2 = Matmul2(YsII,CYsIIYsIICYsIImd2,OnlyDiagonal) 
 YsIICYsIIYsIICYsIIYd = Matmul2(YsII,CYsIIYsIICYsIIYd,OnlyDiagonal) 
 YsIICYsIIYsIICYsIIYsII = Matmul2(YsII,CYsIIYsIICYsIIYsII,OnlyDiagonal) 
 YsIICYsIIYsIICYsIIYzII = Matmul2(YsII,CYsIIYsIICYsIIYzII,OnlyDiagonal) 
 YsIICYsIIYsIICYsIITd = Matmul2(YsII,CYsIIYsIICYsIITd,OnlyDiagonal) 
 YsIICYsIIYsIICYsIITsII = Matmul2(YsII,CYsIIYsIICYsIITsII,OnlyDiagonal) 
 YsIICYsIIYsIICYsIITzII = Matmul2(YsII,CYsIIYsIICYsIITzII,OnlyDiagonal) 
 YsIICYsIIYzIIadjYzIIYsII = Matmul2(YsII,CYsIIYzIIadjYzIIYsII,OnlyDiagonal) 
 YsIICYsIIYzIIadjYzIITsII = Matmul2(YsII,CYsIIYzIIadjYzIITsII,OnlyDiagonal) 
 YsIICYsIITdadjYdYsII = Matmul2(YsII,CYsIITdadjYdYsII,OnlyDiagonal) 
 YsIICYsIITsIICYsIIYd = Matmul2(YsII,CYsIITsIICYsIIYd,OnlyDiagonal) 
 YsIICYsIITsIICYsIIYsII = Matmul2(YsII,CYsIITsIICYsIIYsII,OnlyDiagonal) 
 YsIICYsIITsIICYsIIYzII = Matmul2(YsII,CYsIITsIICYsIIYzII,OnlyDiagonal) 
 YsIICYsIITzIIadjYzIIYsII = Matmul2(YsII,CYsIITzIIadjYzIIYsII,OnlyDiagonal) 
 YsIICYzIIYtIICYtIITpYzII = Matmul2(YsII,CYzIIYtIICYtIITpYzII,OnlyDiagonal) 
 YsIICYzIIYtIICYtIITpTzII = Matmul2(YsII,CYzIIYtIICYtIITpTzII,OnlyDiagonal) 
 YsIICYzIICml2TpYzIICYsII = Matmul2(YsII,CYzIICml2TpYzIICYsII,OnlyDiagonal) 
 YsIICYzIITtIICYtIITpYzII = Matmul2(YsII,CYzIITtIICYtIITpYzII,OnlyDiagonal) 
 YsIICYzIITpYeCYeTpYzII = Matmul2(YsII,CYzIITpYeCYeTpYzII,OnlyDiagonal) 
 YsIICYzIITpYeCYeTpTzII = Matmul2(YsII,CYzIITpYeCYeTpTzII,OnlyDiagonal) 
 YsIICYzIITpYzIICmd2CYsII = Matmul2(YsII,CYzIITpYzIICmd2CYsII,OnlyDiagonal) 
 YsIICYzIITpYzIICYsIImd2 = Matmul2(YsII,CYzIITpYzIICYsIImd2,OnlyDiagonal) 
 YsIICYzIITpYzIICYsIIYd = Matmul2(YsII,CYzIITpYzIICYsIIYd,OnlyDiagonal) 
 YsIICYzIITpYzIICYsIIYsII = Matmul2(YsII,CYzIITpYzIICYsIIYsII,OnlyDiagonal) 
 YsIICYzIITpYzIICYsIIYzII = Matmul2(YsII,CYzIITpYzIICYsIIYzII,OnlyDiagonal) 
 YsIICYzIITpYzIICYsIITd = Matmul2(YsII,CYzIITpYzIICYsIITd,OnlyDiagonal) 
 YsIICYzIITpYzIICYsIITsII = Matmul2(YsII,CYzIITpYzIICYsIITsII,OnlyDiagonal) 
 YsIICYzIITpYzIICYsIITzII = Matmul2(YsII,CYzIITpYzIICYsIITzII,OnlyDiagonal) 
 YsIICYzIITpYzIICYzIITpYzII = Matmul2(YsII,CYzIITpYzIICYzIITpYzII,OnlyDiagonal) 
 YsIICYzIITpYzIICYzIITpTzII = Matmul2(YsII,CYzIITpYzIICYzIITpTzII,OnlyDiagonal) 
 YsIICYzIITpTeCYeTpYzII = Matmul2(YsII,CYzIITpTeCYeTpYzII,OnlyDiagonal) 
 YsIICYzIITpTzIICYsIIYd = Matmul2(YsII,CYzIITpTzIICYsIIYd,OnlyDiagonal) 
 YsIICYzIITpTzIICYsIIYsII = Matmul2(YsII,CYzIITpTzIICYsIIYsII,OnlyDiagonal) 
 YsIICYzIITpTzIICYsIIYzII = Matmul2(YsII,CYzIITpTzIICYsIIYzII,OnlyDiagonal) 
 YsIICYzIITpTzIICYzIITpYzII = Matmul2(YsII,CYzIITpTzIICYzIITpYzII,OnlyDiagonal) 
 YsIITdadjYdCYsIIYd = Matmul2(YsII,TdadjYdCYsIIYd,OnlyDiagonal) 
 YsIITdadjYdCYsIIYzII = Matmul2(YsII,TdadjYdCYsIIYzII,OnlyDiagonal) 
 YsIITzIIadjYzIICYsIIYd = Matmul2(YsII,TzIIadjYzIICYsIIYd,OnlyDiagonal) 
 YsIITzIIadjYzIICYsIIYzII = Matmul2(YsII,TzIIadjYzIICYsIIYzII,OnlyDiagonal) 
 YtIIadjYeYeadjYeYe = Matmul2(YtII,adjYeYeadjYeYe,OnlyDiagonal) 
 YtIIadjYeYeadjYeTe = Matmul2(YtII,adjYeYeadjYeTe,OnlyDiagonal) 
 YtIIadjYeYeCYtIIWOp = Matmul2(YtII,adjYeYeCYtIIWOp,OnlyDiagonal) 
 YtIIadjYeYeCYtIIYtII = Matmul2(YtII,adjYeYeCYtIIYtII,OnlyDiagonal) 
 YtIIadjYeYeCYtIITtII = Matmul2(YtII,adjYeYeCYtIITtII,OnlyDiagonal) 
 YtIIadjYeTeadjYeYe = Matmul2(YtII,adjYeTeadjYeYe,OnlyDiagonal) 
 YtIIadjYeTeCYtIIYtII = Matmul2(YtII,adjYeTeCYtIIYtII,OnlyDiagonal) 
 YtIIadjYzIIYdadjYdYzII = Matmul2(YtII,adjYzIIYdadjYdYzII,OnlyDiagonal) 
 YtIIadjYzIIYdadjYdTzII = Matmul2(YtII,adjYzIIYdadjYdTzII,OnlyDiagonal) 
 YtIIadjYzIIYsIICYsIIYzII = Matmul2(YtII,adjYzIIYsIICYsIIYzII,OnlyDiagonal) 
 YtIIadjYzIIYsIICYsIITzII = Matmul2(YtII,adjYzIIYsIICYsIITzII,OnlyDiagonal) 
 YtIIadjYzIIYzIIadjYzIIYzII = Matmul2(YtII,adjYzIIYzIIadjYzIIYzII,OnlyDiagonal) 
 YtIIadjYzIIYzIIadjYzIITzII = Matmul2(YtII,adjYzIIYzIIadjYzIITzII,OnlyDiagonal) 
 YtIIadjYzIIYzIICYtIIWOp = Matmul2(YtII,adjYzIIYzIICYtIIWOp,OnlyDiagonal) 
 YtIIadjYzIIYzIICYtIIYtII = Matmul2(YtII,adjYzIIYzIICYtIIYtII,OnlyDiagonal) 
 YtIIadjYzIIYzIICYtIICml2 = Matmul2(YtII,adjYzIIYzIICYtIICml2,OnlyDiagonal) 
 YtIIadjYzIIYzIICYtIITtII = Matmul2(YtII,adjYzIIYzIICYtIITtII,OnlyDiagonal) 
 YtIIadjYzIITdadjYdYzII = Matmul2(YtII,adjYzIITdadjYdYzII,OnlyDiagonal) 
 YtIIadjYzIITsIICYsIIYzII = Matmul2(YtII,adjYzIITsIICYsIIYzII,OnlyDiagonal) 
 YtIIadjYzIITzIIadjYzIIYzII = Matmul2(YtII,adjYzIITzIIadjYzIIYzII,OnlyDiagonal) 
 YtIIadjYzIITzIICYtIIYtII = Matmul2(YtII,adjYzIITzIICYtIIYtII,OnlyDiagonal) 
 YtIICYtIIYtIICYtIIWOp = Matmul2(YtII,CYtIIYtIICYtIIWOp,OnlyDiagonal) 
 YtIICYtIIYtIICYtIIYtII = Matmul2(YtII,CYtIIYtIICYtIIYtII,OnlyDiagonal) 
 YtIICYtIIYtIICYtIITtII = Matmul2(YtII,CYtIIYtIICYtIITtII,OnlyDiagonal) 
 YtIICYtIITtIICYtIIYtII = Matmul2(YtII,CYtIITtIICYtIIYtII,OnlyDiagonal) 
 YtIICYtIITpYeCYeYtII = Matmul2(YtII,CYtIITpYeCYeYtII,OnlyDiagonal) 
 YtIICYtIITpYeCYeTtII = Matmul2(YtII,CYtIITpYeCYeTtII,OnlyDiagonal) 
 YtIICYtIITpYzIICYzIIYtII = Matmul2(YtII,CYtIITpYzIICYzIIYtII,OnlyDiagonal) 
 YtIICYtIITpYzIICYzIITtII = Matmul2(YtII,CYtIITpYzIICYzIITtII,OnlyDiagonal) 
 YtIICYtIITpTeCYeYtII = Matmul2(YtII,CYtIITpTeCYeYtII,OnlyDiagonal) 
 YtIICYtIITpTzIICYzIIYtII = Matmul2(YtII,CYtIITpTzIICYzIIYtII,OnlyDiagonal) 
 Yumq2adjYdYdadjYu = Matmul2(Yu,mq2adjYdYdadjYu,OnlyDiagonal) 
 Yumq2adjYuYuadjYu = Matmul2(Yu,mq2adjYuYuadjYu,OnlyDiagonal) 
 YuadjYdmd2YdadjYu = Matmul2(Yu,adjYdmd2YdadjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYdmd2YdadjYu(i2,i2) =  Real(YuadjYdmd2YdadjYu(i2,i2),dp) 
 YuadjYdYdmq2adjYu = Matmul2(Yu,adjYdYdmq2adjYu,OnlyDiagonal) 
 YuadjYdYdadjYdYd = Matmul2(Yu,adjYdYdadjYdYd,OnlyDiagonal) 
 YuadjYdYdadjYdTd = Matmul2(Yu,adjYdYdadjYdTd,OnlyDiagonal) 
 YuadjYdYdadjYumu2 = Matmul2(Yu,adjYdYdadjYumu2,OnlyDiagonal) 
 YuadjYdYdadjYuYu = Matmul2(Yu,adjYdYdadjYuYu,OnlyDiagonal) 
 YuadjYdYdadjYuTu = Matmul2(Yu,adjYdYdadjYuTu,OnlyDiagonal) 
 YuadjYdYsIICYsIIYd = Matmul2(Yu,adjYdYsIICYsIIYd,OnlyDiagonal) 
 YuadjYdYsIICYsIITd = Matmul2(Yu,adjYdYsIICYsIITd,OnlyDiagonal) 
 YuadjYdYzIIadjYzIIYd = Matmul2(Yu,adjYdYzIIadjYzIIYd,OnlyDiagonal) 
 YuadjYdYzIIadjYzIITd = Matmul2(Yu,adjYdYzIIadjYzIITd,OnlyDiagonal) 
 YuadjYdTdadjYdYd = Matmul2(Yu,adjYdTdadjYdYd,OnlyDiagonal) 
 YuadjYdTdadjYuYu = Matmul2(Yu,adjYdTdadjYuYu,OnlyDiagonal) 
 YuadjYdTsIICYsIIYd = Matmul2(Yu,adjYdTsIICYsIIYd,OnlyDiagonal) 
 YuadjYdTzIIadjYzIIYd = Matmul2(Yu,adjYdTzIIadjYzIIYd,OnlyDiagonal) 
 YuadjYumu2YuadjYu = Matmul2(Yu,adjYumu2YuadjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYumu2YuadjYu(i2,i2) =  Real(YuadjYumu2YuadjYu(i2,i2),dp) 
 YuadjYuYumq2adjYu = Matmul2(Yu,adjYuYumq2adjYu,OnlyDiagonal) 
 YuadjYuYuadjYumu2 = Matmul2(Yu,adjYuYuadjYumu2,OnlyDiagonal) 
 YuadjYuYuadjYuYu = Matmul2(Yu,adjYuYuadjYuYu,OnlyDiagonal) 
 YuadjYuYuadjYuTu = Matmul2(Yu,adjYuYuadjYuTu,OnlyDiagonal) 
 YuadjYuTuadjYuYu = Matmul2(Yu,adjYuTuadjYuYu,OnlyDiagonal) 
 YzIIml2adjYeYeadjYzII = Matmul2(YzII,ml2adjYeYeadjYzII,OnlyDiagonal) 
 YzIIml2adjYzIIYzIIadjYzII = Matmul2(YzII,ml2adjYzIIYzIIadjYzII,OnlyDiagonal) 
 YzIIml2CYtIIYtIIadjYzII = Matmul2(YzII,ml2CYtIIYtIIadjYzII,OnlyDiagonal) 
 YzIIadjYeme2YeadjYzII = Matmul2(YzII,adjYeme2YeadjYzII,OnlyDiagonal) 
Forall(i2=1:3)  YzIIadjYeme2YeadjYzII(i2,i2) =  Real(YzIIadjYeme2YeadjYzII(i2,i2),dp) 
 YzIIadjYeYeml2adjYzII = Matmul2(YzII,adjYeYeml2adjYzII,OnlyDiagonal) 
 YzIIadjYeYeadjYeYe = Matmul2(YzII,adjYeYeadjYeYe,OnlyDiagonal) 
 YzIIadjYeYeadjYeTe = Matmul2(YzII,adjYeYeadjYeTe,OnlyDiagonal) 
 YzIIadjYeYeadjYzIImd2 = Matmul2(YzII,adjYeYeadjYzIImd2,OnlyDiagonal) 
 YzIIadjYeYeadjYzIIYd = Matmul2(YzII,adjYeYeadjYzIIYd,OnlyDiagonal) 
 YzIIadjYeYeadjYzIIYsII = Matmul2(YzII,adjYeYeadjYzIIYsII,OnlyDiagonal) 
 YzIIadjYeYeadjYzIIYzII = Matmul2(YzII,adjYeYeadjYzIIYzII,OnlyDiagonal) 
 YzIIadjYeYeadjYzIITd = Matmul2(YzII,adjYeYeadjYzIITd,OnlyDiagonal) 
 YzIIadjYeYeadjYzIITsII = Matmul2(YzII,adjYeYeadjYzIITsII,OnlyDiagonal) 
 YzIIadjYeYeadjYzIITzII = Matmul2(YzII,adjYeYeadjYzIITzII,OnlyDiagonal) 
 YzIIadjYeTeadjYeYe = Matmul2(YzII,adjYeTeadjYeYe,OnlyDiagonal) 
 YzIIadjYeTeadjYzIIYd = Matmul2(YzII,adjYeTeadjYzIIYd,OnlyDiagonal) 
 YzIIadjYeTeadjYzIIYsII = Matmul2(YzII,adjYeTeadjYzIIYsII,OnlyDiagonal) 
 YzIIadjYeTeadjYzIIYzII = Matmul2(YzII,adjYeTeadjYzIIYzII,OnlyDiagonal) 
 YzIIadjYzIImd2YzIIadjYzII = Matmul2(YzII,adjYzIImd2YzIIadjYzII,OnlyDiagonal) 
Forall(i2=1:3)  YzIIadjYzIImd2YzIIadjYzII(i2,i2) =  Real(YzIIadjYzIImd2YzIIadjYzII(i2,i2),dp) 
 YzIIadjYzIIYdadjYdYzII = Matmul2(YzII,adjYzIIYdadjYdYzII,OnlyDiagonal) 
 YzIIadjYzIIYdadjYdTzII = Matmul2(YzII,adjYzIIYdadjYdTzII,OnlyDiagonal) 
 YzIIadjYzIIYsIICYsIIYzII = Matmul2(YzII,adjYzIIYsIICYsIIYzII,OnlyDiagonal) 
 YzIIadjYzIIYsIICYsIITzII = Matmul2(YzII,adjYzIIYsIICYsIITzII,OnlyDiagonal) 
 YzIIadjYzIIYzIIml2adjYzII = Matmul2(YzII,adjYzIIYzIIml2adjYzII,OnlyDiagonal) 
 YzIIadjYzIIYzIIadjYzIImd2 = Matmul2(YzII,adjYzIIYzIIadjYzIImd2,OnlyDiagonal) 
 YzIIadjYzIIYzIIadjYzIIYd = Matmul2(YzII,adjYzIIYzIIadjYzIIYd,OnlyDiagonal) 
 YzIIadjYzIIYzIIadjYzIIYsII = Matmul2(YzII,adjYzIIYzIIadjYzIIYsII,OnlyDiagonal) 
 YzIIadjYzIIYzIIadjYzIIYzII = Matmul2(YzII,adjYzIIYzIIadjYzIIYzII,OnlyDiagonal) 
 YzIIadjYzIIYzIIadjYzIITd = Matmul2(YzII,adjYzIIYzIIadjYzIITd,OnlyDiagonal) 
 YzIIadjYzIIYzIIadjYzIITsII = Matmul2(YzII,adjYzIIYzIIadjYzIITsII,OnlyDiagonal) 
 YzIIadjYzIIYzIIadjYzIITzII = Matmul2(YzII,adjYzIIYzIIadjYzIITzII,OnlyDiagonal) 
 YzIIadjYzIITdadjYdYzII = Matmul2(YzII,adjYzIITdadjYdYzII,OnlyDiagonal) 
 YzIIadjYzIITsIICYsIIYzII = Matmul2(YzII,adjYzIITsIICYsIIYzII,OnlyDiagonal) 
 YzIIadjYzIITzIIadjYzIIYd = Matmul2(YzII,adjYzIITzIIadjYzIIYd,OnlyDiagonal) 
 YzIIadjYzIITzIIadjYzIIYsII = Matmul2(YzII,adjYzIITzIIadjYzIIYsII,OnlyDiagonal) 
 YzIIadjYzIITzIIadjYzIIYzII = Matmul2(YzII,adjYzIITzIIadjYzIIYzII,OnlyDiagonal) 
 YzIICYtIIYtIIml2adjYzII = Matmul2(YzII,CYtIIYtIIml2adjYzII,OnlyDiagonal) 
 YzIICYtIIYtIIadjYzIImd2 = Matmul2(YzII,CYtIIYtIIadjYzIImd2,OnlyDiagonal) 
 YzIICYtIIYtIIadjYzIIYd = Matmul2(YzII,CYtIIYtIIadjYzIIYd,OnlyDiagonal) 
 YzIICYtIIYtIIadjYzIIYsII = Matmul2(YzII,CYtIIYtIIadjYzIIYsII,OnlyDiagonal) 
 YzIICYtIIYtIIadjYzIIYzII = Matmul2(YzII,CYtIIYtIIadjYzIIYzII,OnlyDiagonal) 
 YzIICYtIIYtIIadjYzIITd = Matmul2(YzII,CYtIIYtIIadjYzIITd,OnlyDiagonal) 
 YzIICYtIIYtIIadjYzIITsII = Matmul2(YzII,CYtIIYtIIadjYzIITsII,OnlyDiagonal) 
 YzIICYtIIYtIIadjYzIITzII = Matmul2(YzII,CYtIIYtIIadjYzIITzII,OnlyDiagonal) 
 YzIICYtIIYtIICYtIIYtII = Matmul2(YzII,CYtIIYtIICYtIIYtII,OnlyDiagonal) 
 YzIICYtIIYtIICYtIITtII = Matmul2(YzII,CYtIIYtIICYtIITtII,OnlyDiagonal) 
 YzIICYtIICml2YtIIadjYzII = Matmul2(YzII,CYtIICml2YtIIadjYzII,OnlyDiagonal) 
 YzIICYtIITtIIadjYzIIYd = Matmul2(YzII,CYtIITtIIadjYzIIYd,OnlyDiagonal) 
 YzIICYtIITtIIadjYzIIYsII = Matmul2(YzII,CYtIITtIIadjYzIIYsII,OnlyDiagonal) 
 YzIICYtIITtIIadjYzIIYzII = Matmul2(YzII,CYtIITtIIadjYzIIYzII,OnlyDiagonal) 
 YzIICYtIITtIICYtIIYtII = Matmul2(YzII,CYtIITtIICYtIIYtII,OnlyDiagonal) 
 YzIICYtIITpYeCYeYtII = Matmul2(YzII,CYtIITpYeCYeYtII,OnlyDiagonal) 
 YzIICYtIITpYeCYeTtII = Matmul2(YzII,CYtIITpYeCYeTtII,OnlyDiagonal) 
 YzIICYtIITpYzIICYzIIYtII = Matmul2(YzII,CYtIITpYzIICYzIIYtII,OnlyDiagonal) 
 YzIICYtIITpYzIICYzIITtII = Matmul2(YzII,CYtIITpYzIICYzIITtII,OnlyDiagonal) 
 YzIICYtIITpTeCYeYtII = Matmul2(YzII,CYtIITpTeCYeYtII,OnlyDiagonal) 
 YzIICYtIITpTzIICYzIIYtII = Matmul2(YzII,CYtIITpTzIICYzIIYtII,OnlyDiagonal) 
 adjYdmd2YdadjYdYd = Matmul2(adjYd,md2YdadjYdYd,OnlyDiagonal) 
 adjYdmd2YsIICYsIIYd = Matmul2(adjYd,md2YsIICYsIIYd,OnlyDiagonal) 
 adjYdmd2YzIIadjYzIIYd = Matmul2(adjYd,md2YzIIadjYzIIYd,OnlyDiagonal) 
 adjYdYdmq2adjYdYd = Matmul2(adjYd,Ydmq2adjYdYd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYdmq2adjYdYd(i2,i2) =  Real(adjYdYdmq2adjYdYd(i2,i2),dp) 
 adjYdYdadjYdmd2Yd = Matmul2(adjYd,YdadjYdmd2Yd,OnlyDiagonal) 
 adjYdYdadjYdYdmq2 = Matmul2(adjYd,YdadjYdYdmq2,OnlyDiagonal) 
 adjYdYsIICmd2CYsIIYd = Matmul2(adjYd,YsIICmd2CYsIIYd,OnlyDiagonal) 
 adjYdYsIICYsIImd2Yd = Matmul2(adjYd,YsIICYsIImd2Yd,OnlyDiagonal) 
 adjYdYsIICYsIIYdmq2 = Matmul2(adjYd,YsIICYsIIYdmq2,OnlyDiagonal) 
 adjYdYzIIml2adjYzIIYd = Matmul2(adjYd,YzIIml2adjYzIIYd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYzIIml2adjYzIIYd(i2,i2) =  Real(adjYdYzIIml2adjYzIIYd(i2,i2),dp) 
 adjYdYzIIadjYzIImd2Yd = Matmul2(adjYd,YzIIadjYzIImd2Yd,OnlyDiagonal) 
 adjYdYzIIadjYzIIYdmq2 = Matmul2(adjYd,YzIIadjYzIIYdmq2,OnlyDiagonal) 
 adjYeme2YeadjYeYe = Matmul2(adjYe,me2YeadjYeYe,OnlyDiagonal) 
 adjYeYeml2adjYeYe = Matmul2(adjYe,Yeml2adjYeYe,OnlyDiagonal) 
Forall(i2=1:3)  adjYeYeml2adjYeYe(i2,i2) =  Real(adjYeYeml2adjYeYe(i2,i2),dp) 
 adjYeYeadjYeme2Ye = Matmul2(adjYe,YeadjYeme2Ye,OnlyDiagonal) 
 adjYeYeadjYeYeml2 = Matmul2(adjYe,YeadjYeYeml2,OnlyDiagonal) 
 adjYumu2YuadjYuYu = Matmul2(adjYu,mu2YuadjYuYu,OnlyDiagonal) 
 adjYuYumq2adjYuYu = Matmul2(adjYu,Yumq2adjYuYu,OnlyDiagonal) 
Forall(i2=1:3)  adjYuYumq2adjYuYu(i2,i2) =  Real(adjYuYumq2adjYuYu(i2,i2),dp) 
 adjYuYuadjYumu2Yu = Matmul2(adjYu,YuadjYumu2Yu,OnlyDiagonal) 
 adjYuYuadjYuYumq2 = Matmul2(adjYu,YuadjYuYumq2,OnlyDiagonal) 
 adjYzIImd2YdadjYdYzII = Matmul2(adjYzII,md2YdadjYdYzII,OnlyDiagonal) 
 adjYzIImd2YsIICYsIIYzII = Matmul2(adjYzII,md2YsIICYsIIYzII,OnlyDiagonal) 
 adjYzIImd2YzIIadjYzIIYzII = Matmul2(adjYzII,md2YzIIadjYzIIYzII,OnlyDiagonal) 
 adjYzIIYdmq2adjYdYzII = Matmul2(adjYzII,Ydmq2adjYdYzII,OnlyDiagonal) 
Forall(i2=1:3)  adjYzIIYdmq2adjYdYzII(i2,i2) =  Real(adjYzIIYdmq2adjYdYzII(i2,i2),dp) 
 adjYzIIYdadjYdmd2YzII = Matmul2(adjYzII,YdadjYdmd2YzII,OnlyDiagonal) 
 adjYzIIYdadjYdYzIIml2 = Matmul2(adjYzII,YdadjYdYzIIml2,OnlyDiagonal) 
 adjYzIIYsIICmd2CYsIIYzII = Matmul2(adjYzII,YsIICmd2CYsIIYzII,OnlyDiagonal) 
 adjYzIIYsIICYsIImd2YzII = Matmul2(adjYzII,YsIICYsIImd2YzII,OnlyDiagonal) 
 adjYzIIYsIICYsIIYzIIml2 = Matmul2(adjYzII,YsIICYsIIYzIIml2,OnlyDiagonal) 
 adjYzIIYzIIml2adjYzIIYzII = Matmul2(adjYzII,YzIIml2adjYzIIYzII,OnlyDiagonal) 
Forall(i2=1:3)  adjYzIIYzIIml2adjYzIIYzII(i2,i2) =  Real(adjYzIIYzIIml2adjYzIIYzII(i2,i2),dp) 
 adjYzIIYzIIadjYzIImd2YzII = Matmul2(adjYzII,YzIIadjYzIImd2YzII,OnlyDiagonal) 
 adjYzIIYzIIadjYzIIYzIIml2 = Matmul2(adjYzII,YzIIadjYzIIYzIIml2,OnlyDiagonal) 
 CYtIIYtIIml2CYtIIYtII = Matmul2(adjYtII,YtIIml2CYtIIYtII,OnlyDiagonal) 
 CYtIIYtIICYtIIYtIIml2 = Matmul2(adjYtII,YtIICYtIIYtIIml2,OnlyDiagonal) 
 CYtIIYtIICYtIICml2YtII = Matmul2(adjYtII,YtIICYtIICml2YtII,OnlyDiagonal) 
 CYtIICml2YtIICYtIIYtII = Matmul2(adjYtII,Cml2YtIICYtIIYtII,OnlyDiagonal) 
 CYtIICml2TpYeCYeYtII = Matmul2(adjYtII,Cml2TpYeCYeYtII,OnlyDiagonal) 
 CYtIICml2TpYzIICYzIIYtII = Matmul2(adjYtII,Cml2TpYzIICYzIIYtII,OnlyDiagonal) 
 CYtIITpYeCme2CYeYtII = Matmul2(adjYtII,TpYeCme2CYeYtII,OnlyDiagonal) 
 CYtIITpYeCYeYtIIml2 = Matmul2(adjYtII,TpYeCYeYtIIml2,OnlyDiagonal) 
 CYtIITpYeCYeCml2YtII = Matmul2(adjYtII,TpYeCYeCml2YtII,OnlyDiagonal) 
 CYtIITpYzIICmd2CYzIIYtII = Matmul2(adjYtII,TpYzIICmd2CYzIIYtII,OnlyDiagonal) 
 CYtIITpYzIICYzIIYtIIml2 = Matmul2(adjYtII,TpYzIICYzIIYtIIml2,OnlyDiagonal) 
 CYtIITpYzIICYzIICml2YtII = Matmul2(adjYtII,TpYzIICYzIICml2YtII,OnlyDiagonal) 
 TdadjYdYdadjYdYd = Matmul2(Td,adjYdYdadjYdYd,OnlyDiagonal) 
 TdadjYdYdadjYdYsII = Matmul2(Td,adjYdYdadjYdYsII,OnlyDiagonal) 
 TdadjYdYdadjYdYzII = Matmul2(Td,adjYdYdadjYdYzII,OnlyDiagonal) 
 TdadjYdYsIICYsIIYd = Matmul2(Td,adjYdYsIICYsIIYd,OnlyDiagonal) 
 TdadjYdYzIIadjYzIIYd = Matmul2(Td,adjYdYzIIadjYzIIYd,OnlyDiagonal) 
 TdadjYuYuadjYdYd = Matmul2(Td,adjYuYuadjYdYd,OnlyDiagonal) 
 TdadjYuYuadjYdYsII = Matmul2(Td,adjYuYuadjYdYsII,OnlyDiagonal) 
 TdadjYuYuadjYdYzII = Matmul2(Td,adjYuYuadjYdYzII,OnlyDiagonal) 
 TdadjYuYuadjYuYu = Matmul2(Td,adjYuYuadjYuYu,OnlyDiagonal) 
 TeadjYeYeadjYeYe = Matmul2(Te,adjYeYeadjYeYe,OnlyDiagonal) 
 TeadjYzIIYdadjYdYzII = Matmul2(Te,adjYzIIYdadjYdYzII,OnlyDiagonal) 
 TeadjYzIIYsIICYsIIYzII = Matmul2(Te,adjYzIIYsIICYsIIYzII,OnlyDiagonal) 
 TeadjYzIIYzIIadjYeYe = Matmul2(Te,adjYzIIYzIIadjYeYe,OnlyDiagonal) 
 TeadjYzIIYzIIadjYzIIYzII = Matmul2(Te,adjYzIIYzIIadjYzIIYzII,OnlyDiagonal) 
 TeCYtIIYtIIadjYeYe = Matmul2(Te,CYtIIYtIIadjYeYe,OnlyDiagonal) 
 TeCYtIIYtIICYtIIYtII = Matmul2(Te,CYtIIYtIICYtIIYtII,OnlyDiagonal) 
 TeCYtIITpYeCYeYtII = Matmul2(Te,CYtIITpYeCYeYtII,OnlyDiagonal) 
 TeCYtIITpYzIICYzIIYtII = Matmul2(Te,CYtIITpYzIICYzIIYtII,OnlyDiagonal) 
 TsIICYdTpYdCYdTpYd = Matmul2(TsII,CYdTpYdCYdTpYd,OnlyDiagonal) 
 TsIICYdTpYdCYsIIYd = Matmul2(TsII,CYdTpYdCYsIIYd,OnlyDiagonal) 
 TsIICYdTpYdCYsIIYsII = Matmul2(TsII,CYdTpYdCYsIIYsII,OnlyDiagonal) 
 TsIICYdTpYdCYsIIYzII = Matmul2(TsII,CYdTpYdCYsIIYzII,OnlyDiagonal) 
 TsIICYdTpYuCYuTpYd = Matmul2(TsII,CYdTpYuCYuTpYd,OnlyDiagonal) 
 TsIICYsIIYdadjYdYsII = Matmul2(TsII,CYsIIYdadjYdYsII,OnlyDiagonal) 
 TsIICYsIIYsIICYsIIYd = Matmul2(TsII,CYsIIYsIICYsIIYd,OnlyDiagonal) 
 TsIICYsIIYsIICYsIIYsII = Matmul2(TsII,CYsIIYsIICYsIIYsII,OnlyDiagonal) 
 TsIICYsIIYsIICYsIIYzII = Matmul2(TsII,CYsIIYsIICYsIIYzII,OnlyDiagonal) 
 TsIICYsIIYzIIadjYzIIYsII = Matmul2(TsII,CYsIIYzIIadjYzIIYsII,OnlyDiagonal) 
 TsIICYzIIYtIICYtIITpYzII = Matmul2(TsII,CYzIIYtIICYtIITpYzII,OnlyDiagonal) 
 TsIICYzIITpYeCYeTpYzII = Matmul2(TsII,CYzIITpYeCYeTpYzII,OnlyDiagonal) 
 TsIICYzIITpYzIICYsIIYd = Matmul2(TsII,CYzIITpYzIICYsIIYd,OnlyDiagonal) 
 TsIICYzIITpYzIICYsIIYsII = Matmul2(TsII,CYzIITpYzIICYsIIYsII,OnlyDiagonal) 
 TsIICYzIITpYzIICYsIIYzII = Matmul2(TsII,CYzIITpYzIICYsIIYzII,OnlyDiagonal) 
 TsIICYzIITpYzIICYzIITpYzII = Matmul2(TsII,CYzIITpYzIICYzIITpYzII,OnlyDiagonal) 
 TtIIadjYeYeadjYeYe = Matmul2(TtII,adjYeYeadjYeYe,OnlyDiagonal) 
 TtIIadjYeYeCYtIIYtII = Matmul2(TtII,adjYeYeCYtIIYtII,OnlyDiagonal) 
 TtIIadjYzIIYdadjYdYzII = Matmul2(TtII,adjYzIIYdadjYdYzII,OnlyDiagonal) 
 TtIIadjYzIIYsIICYsIIYzII = Matmul2(TtII,adjYzIIYsIICYsIIYzII,OnlyDiagonal) 
 TtIIadjYzIIYzIIadjYzIIYzII = Matmul2(TtII,adjYzIIYzIIadjYzIIYzII,OnlyDiagonal) 
 TtIIadjYzIIYzIICYtIIYtII = Matmul2(TtII,adjYzIIYzIICYtIIYtII,OnlyDiagonal) 
 TtIICYtIIYtIICYtIIYtII = Matmul2(TtII,CYtIIYtIICYtIIYtII,OnlyDiagonal) 
 TtIICYtIITpYeCYeYtII = Matmul2(TtII,CYtIITpYeCYeYtII,OnlyDiagonal) 
 TtIICYtIITpYzIICYzIIYtII = Matmul2(TtII,CYtIITpYzIICYzIIYtII,OnlyDiagonal) 
 TuadjYdYdadjYdYd = Matmul2(Tu,adjYdYdadjYdYd,OnlyDiagonal) 
 TuadjYdYdadjYuYu = Matmul2(Tu,adjYdYdadjYuYu,OnlyDiagonal) 
 TuadjYdYsIICYsIIYd = Matmul2(Tu,adjYdYsIICYsIIYd,OnlyDiagonal) 
 TuadjYdYzIIadjYzIIYd = Matmul2(Tu,adjYdYzIIadjYzIIYd,OnlyDiagonal) 
 TuadjYuYuadjYuYu = Matmul2(Tu,adjYuYuadjYuYu,OnlyDiagonal) 
 TzIIadjYeYeadjYeYe = Matmul2(TzII,adjYeYeadjYeYe,OnlyDiagonal) 
 TzIIadjYeYeadjYzIIYd = Matmul2(TzII,adjYeYeadjYzIIYd,OnlyDiagonal) 
 TzIIadjYeYeadjYzIIYsII = Matmul2(TzII,adjYeYeadjYzIIYsII,OnlyDiagonal) 
 TzIIadjYeYeadjYzIIYzII = Matmul2(TzII,adjYeYeadjYzIIYzII,OnlyDiagonal) 
 TzIIadjYzIIYdadjYdYzII = Matmul2(TzII,adjYzIIYdadjYdYzII,OnlyDiagonal) 
 TzIIadjYzIIYsIICYsIIYzII = Matmul2(TzII,adjYzIIYsIICYsIIYzII,OnlyDiagonal) 
 TzIIadjYzIIYzIIadjYzIIYd = Matmul2(TzII,adjYzIIYzIIadjYzIIYd,OnlyDiagonal) 
 TzIIadjYzIIYzIIadjYzIIYsII = Matmul2(TzII,adjYzIIYzIIadjYzIIYsII,OnlyDiagonal) 
 TzIIadjYzIIYzIIadjYzIIYzII = Matmul2(TzII,adjYzIIYzIIadjYzIIYzII,OnlyDiagonal) 
 TzIICYtIIYtIIadjYzIIYd = Matmul2(TzII,CYtIIYtIIadjYzIIYd,OnlyDiagonal) 
 TzIICYtIIYtIIadjYzIIYsII = Matmul2(TzII,CYtIIYtIIadjYzIIYsII,OnlyDiagonal) 
 TzIICYtIIYtIIadjYzIIYzII = Matmul2(TzII,CYtIIYtIIadjYzIIYzII,OnlyDiagonal) 
 TzIICYtIIYtIICYtIIYtII = Matmul2(TzII,CYtIIYtIICYtIIYtII,OnlyDiagonal) 
 TzIICYtIITpYeCYeYtII = Matmul2(TzII,CYtIITpYeCYeYtII,OnlyDiagonal) 
 TzIICYtIITpYzIICYzIIYtII = Matmul2(TzII,CYtIITpYzIICYzIIYtII,OnlyDiagonal) 
 TpYeCYeTpYeCYeWOp = Matmul2(Transpose(Ye),CYeTpYeCYeWOp,OnlyDiagonal) 
 TpYeCYeTpYeCYeYtII = Matmul2(Transpose(Ye),CYeTpYeCYeYtII,OnlyDiagonal) 
 TpYeCYeTpYeCYeTtII = Matmul2(Transpose(Ye),CYeTpYeCYeTtII,OnlyDiagonal) 
 TpYeCYeTpTeCYeYtII = Matmul2(Transpose(Ye),CYeTpTeCYeYtII,OnlyDiagonal) 
 TpYzIICYdTpYdCYzIIWOp = Matmul2(Transpose(YzII),CYdTpYdCYzIIWOp,OnlyDiagonal) 
 TpYzIICYdTpYdCYzIIYtII = Matmul2(Transpose(YzII),CYdTpYdCYzIIYtII,OnlyDiagonal) 
 TpYzIICYdTpYdCYzIITtII = Matmul2(Transpose(YzII),CYdTpYdCYzIITtII,OnlyDiagonal) 
 TpYzIICYdTpTdCYzIIYtII = Matmul2(Transpose(YzII),CYdTpTdCYzIIYtII,OnlyDiagonal) 
 TpYzIICYsIIYsIICYzIIWOp = Matmul2(Transpose(YzII),CYsIIYsIICYzIIWOp,OnlyDiagonal) 
 TpYzIICYsIIYsIICYzIIYtII = Matmul2(Transpose(YzII),CYsIIYsIICYzIIYtII,OnlyDiagonal) 
 TpYzIICYsIIYsIICYzIITtII = Matmul2(Transpose(YzII),CYsIIYsIICYzIITtII,OnlyDiagonal) 
 TpYzIICYsIITsIICYzIIYtII = Matmul2(Transpose(YzII),CYsIITsIICYzIIYtII,OnlyDiagonal) 
 TpYzIICYzIITpYzIICYzIIWOp = Matmul2(Transpose(YzII),CYzIITpYzIICYzIIWOp,OnlyDiagonal) 
 TpYzIICYzIITpYzIICYzIIYtII = Matmul2(Transpose(YzII),CYzIITpYzIICYzIIYtII,OnlyDiagonal) 
 TpYzIICYzIITpYzIICYzIITtII = Matmul2(Transpose(YzII),CYzIITpYzIICYzIITtII,OnlyDiagonal) 
 TpYzIICYzIITpTzIICYzIIYtII = Matmul2(Transpose(YzII),CYzIITpTzIICYzIIYtII,OnlyDiagonal) 
 TpTeCYeTpYeCYeYtII = Matmul2(Transpose(Te),CYeTpYeCYeYtII,OnlyDiagonal) 
 TpTzIICYdTpYdCYzIIYtII = Matmul2(Transpose(TzII),CYdTpYdCYzIIYtII,OnlyDiagonal) 
 TpTzIICYsIIYsIICYzIIYtII = Matmul2(Transpose(TzII),CYsIIYsIICYzIIYtII,OnlyDiagonal) 
 TpTzIICYzIITpYzIICYzIIYtII = Matmul2(Transpose(TzII),CYzIITpYzIICYzIIYtII,OnlyDiagonal) 
 TrYsIICTsII = cTrace(YsIICTsII) 
 TrYtIICTtII = cTrace(YtIICTtII) 
 TrCTdTpYd = cTrace(CTdTpYd) 
 TrCTeTpYe = cTrace(CTeTpYe) 
 TrCTuTpYu = cTrace(CTuTpYu) 
 TrCTzIITpYzII = cTrace(CTzIITpYzII) 
 Trmd2CYsIIYsII = cTrace(md2CYsIIYsII) 
 Trml2YtIICYtII = cTrace(ml2YtIICYtII) 
 TrYdadjYdCmd2 = cTrace(YdadjYdCmd2) 
 TrYdCmq2adjYd = cTrace(YdCmq2adjYd) 
 TrYeadjYeCme2 = cTrace(YeadjYeCme2) 
 TrYeCml2adjYe = cTrace(YeCml2adjYe) 
 TrYuadjYuCmu2 = cTrace(YuadjYuCmu2) 
 TrYuCmq2adjYu = cTrace(YuCmq2adjYu) 
 TrYzIIadjYzIICmd2 = cTrace(YzIIadjYzIICmd2) 
 TrYzIICml2adjYzII = cTrace(YzIICml2adjYzII) 
 TrYdadjYdYdadjYd = cTrace(YdadjYdYdadjYd) 
 TrYdadjYdYsIICYsII = cTrace(YdadjYdYsIICYsII) 
 TrYdadjYdYzIIadjYzII = cTrace(YdadjYdYzIIadjYzII) 
 TrYdadjYdTdadjYd = cTrace(YdadjYdTdadjYd) 
 TrYdadjYdTdadjTd = cTrace(YdadjYdTdadjTd) 
 TrYdadjYdTsIICYsII = cTrace(YdadjYdTsIICYsII) 
 TrYdadjYdTsIICTsII = cTrace(YdadjYdTsIICTsII) 
 TrYdadjYdTzIIadjYzII = cTrace(YdadjYdTzIIadjYzII) 
 TrYdadjYdTzIIadjTzII = cTrace(YdadjYdTzIIadjTzII) 
 TrYdadjYuYuadjYd = cTrace(YdadjYuYuadjYd) 
 TrYdadjYuTuadjYd = cTrace(YdadjYuTuadjYd) 
 TrYdadjYuTuadjTd = cTrace(YdadjYuTuadjTd) 
 TrYdadjTdTdadjYd = cTrace(YdadjTdTdadjYd) 
 TrYdadjTdTsIICYsII = cTrace(YdadjTdTsIICYsII) 
 TrYdadjTdTzIIadjYzII = cTrace(YdadjTdTzIIadjYzII) 
 TrYdadjTuTuadjYd = cTrace(YdadjTuTuadjYd) 
 TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe) 
 TrYeadjYeTeadjYe = cTrace(YeadjYeTeadjYe) 
 TrYeadjYeTeadjTe = cTrace(YeadjYeTeadjTe) 
 TrYeadjYzIIYzIIadjYe = cTrace(YeadjYzIIYzIIadjYe) 
 TrYeadjYzIITzIIadjYe = cTrace(YeadjYzIITzIIadjYe) 
 TrYeadjYzIITzIIadjTe = cTrace(YeadjYzIITzIIadjTe) 
 TrYeadjTeTeadjYe = cTrace(YeadjTeTeadjYe) 
 TrYeadjTzIITzIIadjYe = cTrace(YeadjTzIITzIIadjYe) 
 TrYeCYtIIYtIIadjYe = cTrace(YeCYtIIYtIIadjYe) 
 TrYeCYtIITtIIadjYe = cTrace(YeCYtIITtIIadjYe) 
 TrYeCYtIITtIIadjTe = cTrace(YeCYtIITtIIadjTe) 
 TrYeCTtIITtIIadjYe = cTrace(YeCTtIITtIIadjYe) 
 TrYsIICYsIIYsIICYsII = cTrace(YsIICYsIIYsIICYsII) 
 TrYsIICYsIIYzIIadjYzII = cTrace(YsIICYsIIYzIIadjYzII) 
 TrYsIICYsIITdadjYd = cTrace(YsIICYsIITdadjYd) 
 TrYsIICYsIITdadjTd = cTrace(YsIICYsIITdadjTd) 
 TrYsIICYsIITsIICYsII = cTrace(YsIICYsIITsIICYsII) 
 TrYsIICYsIITsIICTsII = cTrace(YsIICYsIITsIICTsII) 
 TrYsIICYsIITzIIadjYzII = cTrace(YsIICYsIITzIIadjYzII) 
 TrYsIICYsIITzIIadjTzII = cTrace(YsIICYsIITzIIadjTzII) 
 TrYsIICTdTpTdCYsII = cTrace(YsIICTdTpTdCYsII) 
 TrYsIICTsIITdadjYd = cTrace(YsIICTsIITdadjYd) 
 TrYsIICTsIITsIICYsII = cTrace(YsIICTsIITsIICYsII) 
 TrYsIICTsIITzIIadjYzII = cTrace(YsIICTsIITzIIadjYzII) 
 TrYsIICTzIITpTzIICYsII = cTrace(YsIICTzIITpTzIICYsII) 
 TrYtIIadjYeTeCYtII = cTrace(YtIIadjYeTeCYtII) 
 TrYtIIadjYeTeCTtII = cTrace(YtIIadjYeTeCTtII) 
 TrYtIIadjYzIIYzIICYtII = cTrace(YtIIadjYzIIYzIICYtII) 
 TrYtIIadjYzIITzIICYtII = cTrace(YtIIadjYzIITzIICYtII) 
 TrYtIIadjYzIITzIICTtII = cTrace(YtIIadjYzIITzIICTtII) 
 TrYtIIadjTeTeCYtII = cTrace(YtIIadjTeTeCYtII) 
 TrYtIIadjTzIITzIICYtII = cTrace(YtIIadjTzIITzIICYtII) 
 TrYtIICYtIIYtIICYtII = cTrace(YtIICYtIIYtIICYtII) 
 TrYtIICYtIITtIICYtII = cTrace(YtIICYtIITtIICYtII) 
 TrYtIICYtIITtIICTtII = cTrace(YtIICYtIITtIICTtII) 
 TrYtIICYtIITpTeCTe = cTrace(YtIICYtIITpTeCTe) 
 TrYtIICYtIITpTzIICTzII = cTrace(YtIICYtIITpTzIICTzII) 
 TrYtIICTtIITtIICYtII = cTrace(YtIICTtIITtIICYtII) 
 TrYuadjYdTdadjYu = cTrace(YuadjYdTdadjYu) 
 TrYuadjYdTdadjTu = cTrace(YuadjYdTdadjTu) 
 TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu) 
 TrYuadjYuTuadjYu = cTrace(YuadjYuTuadjYu) 
 TrYuadjYuTuadjTu = cTrace(YuadjYuTuadjTu) 
 TrYuadjTdTdadjYu = cTrace(YuadjTdTdadjYu) 
 TrYuadjTuTuadjYu = cTrace(YuadjTuTuadjYu) 
 TrYzIIadjYeTeadjYzII = cTrace(YzIIadjYeTeadjYzII) 
 TrYzIIadjYeTeadjTzII = cTrace(YzIIadjYeTeadjTzII) 
 TrYzIIadjYzIIYzIIadjYzII = cTrace(YzIIadjYzIIYzIIadjYzII) 
 TrYzIIadjYzIITdadjYd = cTrace(YzIIadjYzIITdadjYd) 
 TrYzIIadjYzIITdadjTd = cTrace(YzIIadjYzIITdadjTd) 
 TrYzIIadjYzIITsIICYsII = cTrace(YzIIadjYzIITsIICYsII) 
 TrYzIIadjYzIITsIICTsII = cTrace(YzIIadjYzIITsIICTsII) 
 TrYzIIadjYzIITzIIadjYzII = cTrace(YzIIadjYzIITzIIadjYzII) 
 TrYzIIadjYzIITzIIadjTzII = cTrace(YzIIadjYzIITzIIadjTzII) 
 TrYzIIadjTeTeadjYzII = cTrace(YzIIadjTeTeadjYzII) 
 TrYzIIadjTzIITdadjYd = cTrace(YzIIadjTzIITdadjYd) 
 TrYzIIadjTzIITsIICYsII = cTrace(YzIIadjTzIITsIICYsII) 
 TrYzIIadjTzIITzIIadjYzII = cTrace(YzIIadjTzIITzIIadjYzII) 
 TrYzIICYtIITtIIadjYzII = cTrace(YzIICYtIITtIIadjYzII) 
 TrYzIICYtIITtIIadjTzII = cTrace(YzIICYtIITtIIadjTzII) 
 TrYzIICTtIITtIIadjYzII = cTrace(YzIICTtIITtIIadjYzII) 
 TrCYsIITsIICTdTpYd = cTrace(CYsIITsIICTdTpYd) 
 TrCYsIITsIICTzIITpYzII = cTrace(CYsIITsIICTzIITpYzII) 
 TrCYtIITpYeCTeTtII = cTrace(CYtIITpYeCTeTtII) 
 TrCYtIITpYzIICTzIITtII = cTrace(CYtIITpYzIICTzIITtII) 
 Trmd2YdadjYdYdadjYd = cTrace(md2YdadjYdYdadjYd) 
 Trmd2YdadjYdYsIICYsII = cTrace(md2YdadjYdYsIICYsII) 
 Trmd2YdadjYdYzIIadjYzII = cTrace(md2YdadjYdYzIIadjYzII) 
 Trmd2YdadjYuYuadjYd = cTrace(md2YdadjYuYuadjYd) 
 Trmd2YsIICYsIIYdadjYd = cTrace(md2YsIICYsIIYdadjYd) 
 Trmd2YsIICYsIIYsIICYsII = cTrace(md2YsIICYsIIYsIICYsII) 
 Trmd2YsIICYsIIYzIIadjYzII = cTrace(md2YsIICYsIIYzIIadjYzII) 
 Trmd2YzIIadjYeYeadjYzII = cTrace(md2YzIIadjYeYeadjYzII) 
 Trmd2YzIIadjYzIIYdadjYd = cTrace(md2YzIIadjYzIIYdadjYd) 
 Trmd2YzIIadjYzIIYsIICYsII = cTrace(md2YzIIadjYzIIYsIICYsII) 
 Trmd2YzIIadjYzIIYzIIadjYzII = cTrace(md2YzIIadjYzIIYzIIadjYzII) 
 Trmd2YzIICYtIIYtIIadjYzII = cTrace(md2YzIICYtIIYtIIadjYzII) 
 Trme2YeadjYeYeadjYe = cTrace(me2YeadjYeYeadjYe) 
 Trme2YeadjYzIIYzIIadjYe = cTrace(me2YeadjYzIIYzIIadjYe) 
 Trme2YeCYtIIYtIIadjYe = cTrace(me2YeCYtIIYtIIadjYe) 
 Trml2adjYeYeadjYeYe = cTrace(ml2adjYeYeadjYeYe) 
 Trml2adjYeYeadjYzIIYzII = cTrace(ml2adjYeYeadjYzIIYzII) 
 Trml2adjYeYeCYtIIYtII = cTrace(ml2adjYeYeCYtIIYtII) 
 Trml2adjYzIIYdadjYdYzII = cTrace(ml2adjYzIIYdadjYdYzII) 
 Trml2adjYzIIYsIICYsIIYzII = cTrace(ml2adjYzIIYsIICYsIIYzII) 
 Trml2adjYzIIYzIIadjYeYe = cTrace(ml2adjYzIIYzIIadjYeYe) 
 Trml2adjYzIIYzIIadjYzIIYzII = cTrace(ml2adjYzIIYzIIadjYzIIYzII) 
 Trml2adjYzIIYzIICYtIIYtII = cTrace(ml2adjYzIIYzIICYtIIYtII) 
 Trml2CYtIIYtIIadjYeYe = cTrace(ml2CYtIIYtIIadjYeYe) 
 Trml2CYtIIYtIIadjYzIIYzII = cTrace(ml2CYtIIYtIIadjYzIIYzII) 
 Trml2CYtIIYtIICYtIIYtII = cTrace(ml2CYtIIYtIICYtIIYtII) 
 Trmq2adjYdYdadjYdYd = cTrace(mq2adjYdYdadjYdYd) 
 Trmq2adjYdYdadjYuYu = cTrace(mq2adjYdYdadjYuYu) 
 Trmq2adjYdYsIICYsIIYd = cTrace(mq2adjYdYsIICYsIIYd) 
 Trmq2adjYdYzIIadjYzIIYd = cTrace(mq2adjYdYzIIadjYzIIYd) 
 Trmq2adjYuYuadjYdYd = cTrace(mq2adjYuYuadjYdYd) 
 Trmq2adjYuYuadjYuYu = cTrace(mq2adjYuYuadjYuYu) 
 Trmu2YuadjYdYdadjYu = cTrace(mu2YuadjYdYdadjYu) 
 Trmu2YuadjYuYuadjYu = cTrace(mu2YuadjYuYuadjYu) 
 TrYdadjYdYsIICmd2CYsII = cTrace(YdadjYdYsIICmd2CYsII) 
 TrYeCYtIICml2YtIIadjYe = cTrace(YeCYtIICml2YtIIadjYe) 
 TrYsIICmd2CYsIIYzIIadjYzII = cTrace(YsIICmd2CYsIIYzIIadjYzII) 
 TrYtIIadjYzIIYzIICYtIICml2 = cTrace(YtIIadjYzIIYzIICYtIICml2) 
 g1p4 =g1**4 
 g2p4 =g2**4 
 g3p4 =g3**4 
 CL1IIp2 =Conjg(L1II)**2 
 CL2IIp2 =Conjg(L2II)**2 
End If 
 
 
Tr1(1) = g1*sqrt3ov5*(-1._dp*(mHd2) + mHu2 - 4._dp*(ms2) + 4._dp*(msb2)               & 
&  + 3._dp*(mt2) - 3._dp*(mtb2) + mzz2 - mzb2 + Trmd2 + Trme2 - 2._dp*(Trmu2)            & 
&  - Conjg(Trml2) + Conjg(Trmq2))

If (TwoLoopRGE) Then 
Tr2U1(1, 1) = (g1p2*(3._dp*(mHd2) + 3._dp*(mHu2) + 16._dp*(ms2) + 16._dp*(msb2)       & 
&  + 18._dp*(mt2) + 18._dp*(mtb2) + mzz2 + mzb2 + 2._dp*(Trmd2) + 6._dp*(Trme2)          & 
&  + 8._dp*(Trmu2) + 3*Conjg(Trml2) + Conjg(Trmq2)))/10._dp

Tr3(1) = (g1*ooSqrt15*(-9*g1p2*mHd2 - 45*g2p2*mHd2 + 9*g1p2*mHu2 + 45*g2p2*mHu2 -     & 
&  64*g1p2*ms2 - 800*g3p2*ms2 + 64*g1p2*msb2 + 800*g3p2*msb2 + 90*AbsL1II*(mHd2 -        & 
&  mt2) + 108*g1p2*mt2 + 360*g2p2*mt2 - 90*AbsL2II*(mHu2 - mtb2) - 108*g1p2*mtb2 -       & 
&  360*g2p2*mtb2 + g1p2*mzz2 + 45*g2p2*mzz2 + 80*g3p2*mzz2 - g1p2*mzb2 - 45*g2p2*mzb2 -  & 
&  80*g3p2*mzb2 + 4*g1p2*Trmd2 + 80*g3p2*Trmd2 - 120._dp*(Trmd2CYsIIYsII) +              & 
&  36*g1p2*Trme2 + 90._dp*(Trml2YtIICYtII) - 32*g1p2*Trmu2 - 160*g3p2*Trmu2 +            & 
&  90*mHd2*TrYdadjYd - 60._dp*(TrYdadjYdCmd2) - 30._dp*(TrYdCmq2adjYd) + 30*mHd2*TrYeadjYe -& 
&  60._dp*(TrYeadjYeCme2) + 30._dp*(TrYeCml2adjYe) + 120*ms2*TrYsIICYsII -               & 
&  90*mt2*TrYtIICYtII - 90*mHu2*TrYuadjYu + 120._dp*(TrYuadjYuCmu2) - 30._dp*(TrYuCmq2adjYu)& 
&  - 30*mzz2*TrYzIIadjYzII - 60._dp*(TrYzIIadjYzIICmd2) + 90._dp*(TrYzIICml2adjYzII)     & 
&  - 9*g1p2*Conjg(Trml2) - 45*g2p2*Conjg(Trml2) + g1p2*Conjg(Trmq2) + 45*g2p2*Conjg(Trmq2)& 
&  + 80*g3p2*Conjg(Trmq2)))/20._dp

Tr2(2) = (mHd2 + mHu2 + 6._dp*(mt2) + 6._dp*(mtb2) + 3._dp*(mzz2) + 3._dp*(mzb2)      & 
&  + Conjg(Trml2) + 3*Conjg(Trmq2))/2._dp

Tr2(3) = 5._dp*(ms2) + 5._dp*(msb2) + mzz2 + mzb2 + Trmd2/2._dp + Trmu2/2._dp +       & 
&  Conjg(Trmq2)

End If 
 
 
!-------------------- 
! g1 
!-------------------- 
 
betag11  = 68._dp*(g1p3)/5._dp

 
 
If (TwoLoopRGE) Then 
betag12 = (g1p3*(-405._dp*(AbsL1II) - 405._dp*(AbsL2II) + 1502._dp*(g1p2) + 2610._dp*(g2p2) +   & 
&  4600._dp*(g3p2) - 210._dp*(TrYdadjYd) - 270._dp*(TrYeadjYe) - 360._dp*(TrYsIICYsII) - & 
&  405._dp*(TrYtIICYtII) - 390._dp*(TrYuadjYu) - 210._dp*(TrYzIIadjYzII)))/75._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21  = 8._dp*(g2p3)

 
 
If (TwoLoopRGE) Then 
betag22 = (g2p3*(-35._dp*(AbsL1II) - 35._dp*(AbsL2II) + 58._dp*(g1p2) + 470._dp*(g2p2) +        & 
&  200._dp*(g3p2) - 30._dp*(TrYdadjYd) - 10._dp*(TrYeadjYe) - 35._dp*(TrYtIICYtII) -     & 
&  30._dp*(TrYuadjYu) - 30._dp*(TrYzIIadjYzII)))/5._dp

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31  = 4._dp*(g3p3)

 
 
If (TwoLoopRGE) Then 
betag32 = (g3p3*(23._dp*(g1p2) + 45._dp*(g2p2) + 400._dp*(g3p2) - 12._dp*(TrYdadjYd) -          & 
&  27._dp*(TrYsIICYsII) - 12._dp*(TrYuadjYu) - 12._dp*(TrYzIIadjYzII)))/3._dp

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yu 
!-------------------- 
 
betaYu1  = (3._dp*(AbsL2II) - 13._dp*(g1p2)/15._dp - 3._dp*(g2p2) - 16._dp*(g3p2)     & 
& /3._dp + 3._dp*(TrYuadjYu))*Yu + YuadjYdYd + 3._dp*(YuadjYuYu)

 
 
If (TwoLoopRGE) Then 
betaYu2 = (5473._dp*(g1p4)/450._dp + g1p2*g2p2 + 57._dp*(g2p4)/2._dp + (136*g1p2*g3p2)/45._dp + & 
&  8*g2p2*g3p2 + 320._dp*(g3p4)/9._dp - 12*CL2IIp2*L2IIp2 - 3._dp*(TrYdadjYuYuadjYd) +   & 
&  (3*AbsL2II*(6._dp*(g1p2) + 20._dp*(g2p2) - 15._dp*(TrYuadjYu)))/5._dp +               & 
&  (4*(g1p2 + 20._dp*(g3p2))*TrYuadjYu)/5._dp - 9._dp*(TrYuadjYuYuadjYu))*Yu +           & 
&  (-3._dp*(AbsL1II) + 2._dp*(g1p2)/5._dp - 3._dp*(TrYdadjYd) - TrYeadjYe)*YuadjYdYd -   & 
&  2._dp*(YuadjYdYdadjYdYd) - 2._dp*(YuadjYdYdadjYuYu) - 4._dp*(YuadjYdYsIICYsIIYd) -    & 
&  2._dp*(YuadjYdYzIIadjYzIIYd) - 9*AbsL2II*YuadjYuYu + (2*g1p2*YuadjYuYu)/5._dp +       & 
&  6*g2p2*YuadjYuYu - 9*TrYuadjYu*YuadjYuYu - 4._dp*(YuadjYuYuadjYuYu)

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
betaYd1  = (3._dp*(AbsL1II) - 7._dp*(g1p2)/15._dp - 3._dp*(g2p2) - 16._dp*(g3p2)      & 
& /3._dp + 3._dp*(TrYdadjYd) + TrYeadjYe)*Yd + 3._dp*(YdadjYdYd) + YdadjYuYu +           & 
&  4._dp*(YsIICYsIIYd) + 2._dp*(YzIIadjYzIIYd)

 
 
If (TwoLoopRGE) Then 
betaYd2 = (581._dp*(g1p4)/90._dp + g1p2*g2p2 + 57._dp*(g2p4)/2._dp + (8*g1p2*g3p2)/9._dp +      & 
&  8*g2p2*g3p2 + 320._dp*(g3p4)/9._dp - 12*CL1IIp2*L1IIp2 - (2*(g1p2 - 40._dp*(g3p2))*TrYdadjYd)/5._dp -& 
&  9._dp*(TrYdadjYdYdadjYd) - 12._dp*(TrYdadjYdYsIICYsII) - 6._dp*(TrYdadjYdYzIIadjYzII) -& 
&  3._dp*(TrYdadjYuYuadjYd) + (6*g1p2*TrYeadjYe)/5._dp - 3._dp*(TrYeadjYeYeadjYe) -      & 
&  3._dp*(TrYeadjYzIIYzIIadjYe) - 3._dp*(TrYeCYtIIYtIIadjYe) + (3*AbsL1II*(6._dp*(g1p2) +& 
&  20._dp*(g2p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe) - 5._dp*(TrYtIICYtII)))/5._dp)*Yd +& 
&  (-9._dp*(AbsL1II) + 4._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 9._dp*(TrYdadjYd) -           & 
&  3._dp*(TrYeadjYe))*YdadjYdYd - 4._dp*(YdadjYdYdadjYdYd) - 4._dp*(YdadjYdYsIICYsIIYd) -& 
&  2._dp*(YdadjYdYzIIadjYzIIYd) - 3*AbsL2II*YdadjYuYu + (4*g1p2*YdadjYuYu)/5._dp -       & 
&  3*TrYuadjYu*YdadjYuYu - 2._dp*(YdadjYuYuadjYdYd) - 2._dp*(YdadjYuYuadjYuYu) -         & 
&  8._dp*(YsIICYdTpYdCYsIIYd) + (32*g1p2*YsIICYsIIYd)/15._dp + (80*g3p2*YsIICYsIIYd)/3._dp -& 
&  4*TrYsIICYsII*YsIICYsIIYd - 16._dp*(YsIICYsIIYsIICYsIIYd) - 8._dp*(YsIICYzIITpYzIICYsIIYd) -& 
&  2._dp*(YzIIadjYeYeadjYzIIYd) + (2*g1p2*YzIIadjYzIIYd)/5._dp + 6*g2p2*YzIIadjYzIIYd -  & 
&  2*TrYzIIadjYzII*YzIIadjYzIIYd - 6._dp*(YzIIadjYzIIYzIIadjYzIIYd) - 6._dp*(YzIICYtIIYtIIadjYzIIYd)

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
!-------------------- 
! Ye 
!-------------------- 
 
betaYe1  = (3._dp*(AbsL1II) - 9._dp*(g1p2)/5._dp - 3._dp*(g2p2) + 3._dp*(TrYdadjYd)   & 
&  + TrYeadjYe)*Ye + 3*(YeadjYeYe + YeadjYzIIYzII + YeCYtIIYtII)

 
 
If (TwoLoopRGE) Then 
betaYe2 = -((120*CL1IIp2*L1IIp2 + 4*(g1p2 - 40._dp*(g3p2))*TrYdadjYd + 3*(-87._dp*(g1p4) -      & 
&  6*g1p2*g2p2 - 95._dp*(g2p4) + 30._dp*(TrYdadjYdYdadjYd) + 40._dp*(TrYdadjYdYsIICYsII) +& 
&  20._dp*(TrYdadjYdYzIIadjYzII) + 10._dp*(TrYdadjYuYuadjYd) - 4*g1p2*TrYeadjYe +        & 
&  10._dp*(TrYeadjYeYeadjYe) + 10._dp*(TrYeadjYzIIYzIIadjYe) + 10._dp*(TrYeCYtIIYtIIadjYe)) -& 
&  6*AbsL1II*(6._dp*(g1p2) + 20._dp*(g2p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe) -    & 
&  5._dp*(TrYtIICYtII)))*Ye)/10._dp + (-9._dp*(AbsL1II) + 6._dp*(g2p2) - 9._dp*(TrYdadjYd) -& 
&  3._dp*(TrYeadjYe))*YeadjYeYe - 4._dp*(YeadjYeYeadjYeYe) - 6._dp*(YeadjYzIIYdadjYdYzII) -& 
&  12._dp*(YeadjYzIIYsIICYsIIYzII) - (2*g1p2*YeadjYzIIYzII)/5._dp + 16*g3p2*YeadjYzIIYzII -& 
&  3*TrYzIIadjYzII*YeadjYzIIYzII - 6._dp*(YeadjYzIIYzIIadjYeYe) - 6._dp*(YeadjYzIIYzIIadjYzIIYzII) -& 
&  3._dp*(YeCYtIITpYeCYeYtII) - 9._dp*(YeCYtIITpYzIICYzIIYtII) - 3*AbsL1II*YeCYtIIYtII + & 
&  (18*g1p2*YeCYtIIYtII)/5._dp + 12*g2p2*YeCYtIIYtII - 3*TrYtIICYtII*YeCYtIIYtII -       & 
&  6._dp*(YeCYtIIYtIIadjYeYe) - 9._dp*(YeCYtIIYtIICYtIIYtII)

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
!-------------------- 
! YtII 
!-------------------- 
 
betaYtII1  = TpYeCYeYtII + 3._dp*(TpYzIICYzIIYtII) + (AbsL1II - 9._dp*(g1p2)          & 
& /5._dp - 7._dp*(g2p2) + TrYtIICYtII)*YtII + YtIIadjYeYe + 3._dp*(YtIIadjYzIIYzII)      & 
&  + 6._dp*(YtIICYtIIYtII)

 
 
If (TwoLoopRGE) Then 
betaYtII2 = ((261._dp*(g1p4) + 114*g1p2*g2p2 + 765._dp*(g2p4) - 60*CL1IIp2*L1IIp2 -               & 
&  2*AbsL1II*(3._dp*(g1p2) + 5._dp*(g2p2) + 30._dp*(TrYdadjYd) + 10._dp*(TrYeadjYe)) -   & 
&  20._dp*(TrYeCYtIIYtIIadjYe) - 60._dp*(TrYtIIadjYzIIYzIICYtII) - 2*(3._dp*(g1p2) +     & 
&  5._dp*(g2p2))*TrYtIICYtII - 60._dp*(TrYtIICYtIIYtIICYtII))*YtII)/10._dp +             & 
&  (-3._dp*(AbsL1II) + 6._dp*(g1p2)/5._dp - 3._dp*(TrYdadjYd) - TrYeadjYe)*YtIIadjYeYe + & 
&  (-10._dp*(TpYeCYeTpYeCYeYtII) - 15*AbsL1II*TpYeCYeYtII + 6*g1p2*TpYeCYeYtII -         & 
&  30._dp*(TpYzIICYdTpYdCYzIIYtII) - 60._dp*(TpYzIICYsIIYsIICYzIIYtII) - 30._dp*(TpYzIICYzIITpYzIICYzIIYtII) -& 
&  2*g1p2*TpYzIICYzIIYtII + 80*g3p2*TpYzIICYzIIYtII - 15*TpYeCYeYtII*TrYdadjYd -         & 
&  5*TpYeCYeYtII*TrYeadjYe - 15*TpYzIICYzIIYtII*TrYzIIadjYzII - 10._dp*(YtIIadjYeYeadjYeYe) -& 
&  15._dp*(YtIIadjYeYeCYtIIYtII) - 30._dp*(YtIIadjYzIIYdadjYdYzII) - 60._dp*(YtIIadjYzIIYsIICYsIIYzII) +& 
&  (-2._dp*(g1p2) + 80._dp*(g3p2) - 15._dp*(TrYzIIadjYzII))*YtIIadjYzIIYzII -            & 
&  30._dp*(YtIIadjYzIIYzIIadjYzIIYzII) - 45._dp*(YtIIadjYzIIYzIICYtIIYtII) -             & 
&  15._dp*(YtIICYtIITpYeCYeYtII) - 45._dp*(YtIICYtIITpYzIICYzIIYtII) + 6*(-              & 
& 5._dp*(AbsL1II) + 6._dp*(g1p2) + 20._dp*(g2p2) - 5._dp*(TrYtIICYtII))*YtIICYtIIYtII -  & 
&  90._dp*(YtIICYtIIYtIICYtIIYtII))/5._dp

 
DYtII = oo16pi2*( betaYtII1 + oo16pi2 * betaYtII2 ) 

 
Else 
DYtII = oo16pi2* betaYtII1 
End If 
 
 
!-------------------- 
! YsII 
!-------------------- 
 
betaYsII1  = ((-4*(g1p2 + 15._dp*(g3p2)))/5._dp + TrYsIICYsII)*YsII + 2*(YdadjYdYsII +& 
&  YsIICYdTpYd + 4._dp*(YsIICYsIIYsII) + YsIICYzIITpYzII + YzIIadjYzIIYsII)

 
 
If (TwoLoopRGE) Then 
betaYsII2 = -2._dp*(YdadjYdYdadjYdYsII) + (-6._dp*(AbsL1II) + 2._dp*(g1p2)/5._dp + 6._dp*(g2p2) - & 
&  6._dp*(TrYdadjYd) - 2._dp*(TrYeadjYe))*YdadjYdYsII - 2._dp*(YdadjYuYuadjYdYsII) +     & 
&  (4*(42._dp*(g1p4) + 32*g1p2*g3p2 + 400._dp*(g3p4) - 15._dp*(TrYdadjYdYsIICYsII) -     & 
&  (g1p2 + 5._dp*(g3p2))*TrYsIICYsII - 30._dp*(TrYsIICYsIIYsIICYsII) - 15._dp*(TrYsIICYsIIYzIIadjYzII))*YsII)/15._dp -& 
&  6*AbsL1II*YsIICYdTpYd + (2*g1p2*YsIICYdTpYd)/5._dp + 6*g2p2*YsIICYdTpYd -             & 
&  6*TrYdadjYd*YsIICYdTpYd - 2*TrYeadjYe*YsIICYdTpYd - 2._dp*(YsIICYdTpYdCYdTpYd) -      & 
&  8._dp*(YsIICYdTpYdCYsIIYsII) - 2._dp*(YsIICYdTpYuCYuTpYd) - 8._dp*(YsIICYsIIYdadjYdYsII) +& 
&  (64*g1p2*YsIICYsIIYsII)/15._dp + (160*g3p2*YsIICYsIIYsII)/3._dp - 8*TrYsIICYsII*YsIICYsIIYsII -& 
&  32._dp*(YsIICYsIIYsIICYsIIYsII) - 8._dp*(YsIICYsIIYzIIadjYzIIYsII) - 2._dp*(YsIICYzIITpYeCYeTpYzII) +& 
&  (2*g1p2*YsIICYzIITpYzII)/5._dp + 6*g2p2*YsIICYzIITpYzII - 2*TrYzIIadjYzII*YsIICYzIITpYzII -& 
&  8._dp*(YsIICYzIITpYzIICYsIIYsII) - 6._dp*(YsIICYzIITpYzIICYzIITpYzII) -               & 
&  6._dp*(YsIICYzIIYtIICYtIITpYzII) - 2._dp*(YzIIadjYeYeadjYzIIYsII) + (2*g1p2*YzIIadjYzIIYsII)/5._dp +& 
&  6*g2p2*YzIIadjYzIIYsII - 2*TrYzIIadjYzII*YzIIadjYzIIYsII - 6._dp*(YzIIadjYzIIYzIIadjYzIIYsII) -& 
&  6._dp*(YzIICYtIIYtIIadjYzIIYsII)

 
DYsII = oo16pi2*( betaYsII1 + oo16pi2 * betaYsII2 ) 

 
Else 
DYsII = oo16pi2* betaYsII1 
End If 
 
 
!-------------------- 
! YzII 
!-------------------- 
 
betaYzII1  = 2._dp*(YdadjYdYzII) + 4._dp*(YsIICYsIIYzII) + (-7._dp*(g1p2)             & 
& /15._dp - 3._dp*(g2p2) - 16._dp*(g3p2)/3._dp + TrYzIIadjYzII)*YzII + YzIIadjYeYe +     & 
&  5._dp*(YzIIadjYzIIYzII) + 3._dp*(YzIICYtIIYtII)

 
 
If (TwoLoopRGE) Then 
betaYzII2 = -2._dp*(YdadjYdYdadjYdYzII) + (-6._dp*(AbsL1II) + 2._dp*(g1p2)/5._dp + 6._dp*(g2p2) - & 
&  6._dp*(TrYdadjYd) - 2._dp*(TrYeadjYe))*YdadjYdYzII - 2._dp*(YdadjYuYuadjYdYzII) -     & 
&  8._dp*(YsIICYdTpYdCYsIIYzII) - 16._dp*(YsIICYsIIYsIICYsIIYzII) + (32*g1p2*YsIICYsIIYzII)/15._dp +& 
&  (80*g3p2*YsIICYsIIYzII)/3._dp - 4*TrYsIICYsII*YsIICYsIIYzII - 8._dp*(YsIICYzIITpYzIICYsIIYzII) +& 
&  (581._dp*(g1p4)/90._dp + g1p2*g2p2 + 57._dp*(g2p4)/2._dp + (8*g1p2*g3p2)/9._dp +      & 
&  8*g2p2*g3p2 + 320._dp*(g3p4)/9._dp - 2._dp*(TrYdadjYdYzIIadjYzII) - TrYeadjYzIIYzIIadjYe -& 
&  4._dp*(TrYsIICYsIIYzIIadjYzII) - 3._dp*(TrYtIIadjYzIIYzIICYtII) + (2*g1p2*TrYzIIadjYzII)/5._dp -& 
&  5._dp*(TrYzIIadjYzIIYzIIadjYzII))*YzII - 3*AbsL1II*YzIIadjYeYe + (6*g1p2*YzIIadjYeYe)/5._dp -& 
&  3*TrYdadjYd*YzIIadjYeYe - TrYeadjYe*YzIIadjYeYe - 2._dp*(YzIIadjYeYeadjYeYe) -        & 
&  2._dp*(YzIIadjYeYeadjYzIIYzII) - 6._dp*(YzIIadjYzIIYdadjYdYzII) - 12._dp*(YzIIadjYzIIYsIICYsIIYzII) +& 
&  6*g2p2*YzIIadjYzIIYzII + 16*g3p2*YzIIadjYzIIYzII - 5*TrYzIIadjYzII*YzIIadjYzIIYzII -  & 
&  12._dp*(YzIIadjYzIIYzIIadjYzIIYzII) - 3._dp*(YzIICYtIITpYeCYeYtII) - 9._dp*(YzIICYtIITpYzIICYzIIYtII) -& 
&  3*AbsL1II*YzIICYtIIYtII + (18*g1p2*YzIICYtIIYtII)/5._dp + 12*g2p2*YzIICYtIIYtII -     & 
&  3*TrYtIICYtII*YzIICYtIIYtII - 6._dp*(YzIICYtIIYtIIadjYzIIYzII) - 9._dp*(YzIICYtIIYtIICYtIIYtII)

 
DYzII = oo16pi2*( betaYzII1 + oo16pi2 * betaYzII2 ) 

 
Else 
DYzII = oo16pi2* betaYzII1 
End If 
 
 
!-------------------- 
! L1II 
!-------------------- 
 
betaL1II1  = (-9*g1p2*L1II)/5._dp - 7*g2p2*L1II + 6*L1II*TrYdadjYd + 2*L1II*TrYeadjYe +& 
&  L1II*TrYtIICYtII + 7*L1IIp2*Conjg(L1II)

 
 
If (TwoLoopRGE) Then 
betaL1II2 = -(L1II*(-261._dp*(g1p4) - 114*g1p2*g2p2 - 765._dp*(g2p4) + 300*CL1IIp2*L1IIp2 +       & 
&  8*(g1p2 - 40._dp*(g3p2))*TrYdadjYd + 180._dp*(TrYdadjYdYdadjYd) + 240._dp*(TrYdadjYdYsIICYsII) +& 
&  120._dp*(TrYdadjYdYzIIadjYzII) + 60._dp*(TrYdadjYuYuadjYd) - 24*g1p2*TrYeadjYe +      & 
&  60._dp*(TrYeadjYeYeadjYe) + 60._dp*(TrYeadjYzIIYzIIadjYe) + 80._dp*(TrYeCYtIIYtIIadjYe) +& 
&  60._dp*(TrYtIIadjYzIIYzIICYtII) + 6*g1p2*TrYtIICYtII + 10*g2p2*TrYtIICYtII +          & 
&  2*AbsL1II*(-33._dp*(g1p2) - 115._dp*(g2p2) + 120._dp*(TrYdadjYd) + 40._dp*(TrYeadjYe) +& 
&  30._dp*(TrYtIICYtII)) + 60._dp*(TrYtIICYtIIYtIICYtII)))/10._dp

 
DL1II = oo16pi2*( betaL1II1 + oo16pi2 * betaL1II2 ) 

 
Else 
DL1II = oo16pi2* betaL1II1 
End If 
 
 
!-------------------- 
! L2II 
!-------------------- 
 
betaL2II1  = (-9*g1p2*L2II)/5._dp - 7*g2p2*L2II + 6*L2II*TrYuadjYu + 7*L2IIp2*Conjg(L2II)

 
 
If (TwoLoopRGE) Then 
betaL2II2 = -(L2II*(300*CL2IIp2*L2IIp2 - 2*AbsL2II*(33._dp*(g1p2) + 115._dp*(g2p2) -              & 
&  120._dp*(TrYuadjYu)) - 16*(g1p2 + 20._dp*(g3p2))*TrYuadjYu - 3*(87._dp*(g1p4) +       & 
&  38*g1p2*g2p2 + 255._dp*(g2p4) - 20._dp*(TrYdadjYuYuadjYd) - 60._dp*(TrYuadjYuYuadjYu))))/10._dp

 
DL2II = oo16pi2*( betaL2II1 + oo16pi2 * betaL2II2 ) 

 
Else 
DL2II = oo16pi2* betaL2II1 
End If 
 
 
!-------------------- 
! Mu 
!-------------------- 
 
betaMu1  = 3*AbsL1II*Mu + 3*AbsL2II*Mu - (3*g1p2*Mu)/5._dp - 3*g2p2*Mu +              & 
&  3*TrYdadjYd*Mu + TrYeadjYe*Mu + 3*TrYuadjYu*Mu

 
 
If (TwoLoopRGE) Then 
betaMu2 = ((417._dp*(g1p4) + 90*g1p2*g2p2 + 1425._dp*(g2p4) - 600*CL1IIp2*L1IIp2 -              & 
&  600*CL2IIp2*L2IIp2 - 20*g1p2*TrYdadjYd + 800*g3p2*TrYdadjYd - 450._dp*(TrYdadjYdYdadjYd) -& 
&  600._dp*(TrYdadjYdYsIICYsII) - 300._dp*(TrYdadjYdYzIIadjYzII) - 300._dp*(TrYdadjYuYuadjYd) +& 
&  60*g1p2*TrYeadjYe - 150._dp*(TrYeadjYeYeadjYe) - 150._dp*(TrYeadjYzIIYzIIadjYe) -     & 
&  150._dp*(TrYeCYtIIYtIIadjYe) + 30*AbsL1II*(6._dp*(g1p2) + 20._dp*(g2p2) -             & 
&  15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe) - 5._dp*(TrYtIICYtII)) + 30*AbsL2II*(6._dp*(g1p2) +& 
&  20._dp*(g2p2) - 15._dp*(TrYuadjYu)) + 40*g1p2*TrYuadjYu + 800*g3p2*TrYuadjYu -        & 
&  450._dp*(TrYuadjYuYuadjYu))*Mu)/50._dp

 
DMu = oo16pi2*( betaMu1 + oo16pi2 * betaMu2 ) 

 
Else 
DMu = oo16pi2* betaMu1 
End If 
 
 
!-------------------- 
! MTII 
!-------------------- 
 
betaMTII1  = AbsL1II*MTII + AbsL2II*MTII - (12*g1p2*MTII)/5._dp - 8*g2p2*MTII +       & 
&  MTII*TrYtIICYtII

 
 
If (TwoLoopRGE) Then 
betaMTII2 = (MTII*(888._dp*(g1p4) + 480*g1p2*g2p2 + 2400._dp*(g2p4) - 150*CL1IIp2*L1IIp2 -        & 
&  150*CL2IIp2*L2IIp2 - 5*AbsL1II*(3._dp*(g1p2) + 5._dp*(g2p2) + 30._dp*(TrYdadjYd) +    & 
&  10._dp*(TrYeadjYe)) - 50._dp*(TrYeCYtIIYtIIadjYe) - 150._dp*(TrYtIIadjYzIIYzIICYtII) -& 
&  15*g1p2*TrYtIICYtII - 25*g2p2*TrYtIICYtII - 150._dp*(TrYtIICYtIIYtIICYtII) -          & 
&  5*AbsL2II*(3._dp*(g1p2) + 5._dp*(g2p2) + 30._dp*(TrYuadjYu))))/25._dp

 
DMTII = oo16pi2*( betaMTII1 + oo16pi2 * betaMTII2 ) 

 
Else 
DMTII = oo16pi2* betaMTII1 
End If 
 
 
!-------------------- 
! MZII 
!-------------------- 
 
betaMZII1  = -(MZII*(g1p2 + 45._dp*(g2p2) + 80._dp*(g3p2) - 15._dp*(TrYzIIadjYzII)))/15._dp

 
 
If (TwoLoopRGE) Then 
betaMZII2 = (409*g1p4*MZII)/450._dp + (g1p2*g2p2*MZII)/5._dp + (57*g2p4*MZII)/2._dp +             & 
&  (16*g1p2*g3p2*MZII)/45._dp + 16*g2p2*g3p2*MZII + (320*g3p4*MZII)/9._dp -              & 
&  2*MZII*TrYdadjYdYzIIadjYzII - MZII*TrYeadjYzIIYzIIadjYe - 4*MZII*TrYsIICYsIIYzIIadjYzII -& 
&  3*MZII*TrYtIIadjYzIIYzIICYtII + (2*g1p2*MZII*TrYzIIadjYzII)/5._dp - 5*MZII*TrYzIIadjYzIIYzIIadjYzII

 
DMZII = oo16pi2*( betaMZII1 + oo16pi2 * betaMZII2 ) 

 
Else 
DMZII = oo16pi2* betaMZII1 
End If 
 
 
!-------------------- 
! MSII 
!-------------------- 
 
betaMSII1  = (-8*(2._dp*(g1p2) + 25._dp*(g3p2))*MSII)/15._dp + MSII*TrYsIICYsII

 
 
If (TwoLoopRGE) Then 
betaMSII2 = (4*MSII*(848._dp*(g1p4) + 800*g1p2*g3p2 + 8000._dp*(g3p4) - 225._dp*(TrYdadjYdYsIICYsII) -& 
&  15*(g1p2 + 5._dp*(g3p2))*TrYsIICYsII - 450._dp*(TrYsIICYsIIYsIICYsII) -               & 
&  225._dp*(TrYsIICYsIIYzIIadjYzII)))/225._dp

 
DMSII = oo16pi2*( betaMSII1 + oo16pi2 * betaMSII2 ) 

 
Else 
DMSII = oo16pi2* betaMSII1 
End If 
 
 
!-------------------- 
! Tu 
!-------------------- 
 
betaTu1  = TuadjYdYd + 5._dp*(TuadjYuYu) + 2._dp*(YuadjYdTd) + 4._dp*(YuadjYuTu)      & 
&  + Yu*((26*g1p2*M1)/15._dp + (32*g3p2*M3)/3._dp + 6*g2p2*M2 + 6._dp*(TradjYuTu)        & 
&  + 6*Conjg(L2II)*TL2II) + 3*AbsL2II*Tu - (13*g1p2*Tu)/15._dp - 3*g2p2*Tu -             & 
&  (16*g3p2*Tu)/3._dp + 3*TrYuadjYu*Tu

 
 
If (TwoLoopRGE) Then 
betaTu2 = -3*AbsL1II*TuadjYdYd + (2*g1p2*TuadjYdYd)/5._dp - 3*TrYdadjYd*TuadjYdYd -             & 
&  TrYeadjYe*TuadjYdYd - 2._dp*(TuadjYdYdadjYdYd) - 4._dp*(TuadjYdYdadjYuYu) -           & 
&  4._dp*(TuadjYdYsIICYsIIYd) - 2._dp*(TuadjYdYzIIadjYzIIYd) - 15*AbsL2II*TuadjYuYu +    & 
&  12*g2p2*TuadjYuYu - 15*TrYuadjYu*TuadjYuYu - 6._dp*(TuadjYuYuadjYuYu) -               & 
&  6*AbsL1II*YuadjYdTd + (4*g1p2*YuadjYdTd)/5._dp - 6*TrYdadjYd*YuadjYdTd -              & 
&  2*TrYeadjYe*YuadjYdTd - 4._dp*(YuadjYdTdadjYdYd) - 4._dp*(YuadjYdTdadjYuYu) -         & 
&  8._dp*(YuadjYdTsIICYsIIYd) - 4._dp*(YuadjYdTzIIadjYzIIYd) - 4._dp*(YuadjYdYdadjYdTd) -& 
&  2._dp*(YuadjYdYdadjYuTu) - 8._dp*(YuadjYdYsIICYsIITd) - 4._dp*(YuadjYdYzIIadjYzIITd) -& 
&  12*AbsL2II*YuadjYuTu + (6*g1p2*YuadjYuTu)/5._dp + 6*g2p2*YuadjYuTu - 12*TrYuadjYu*YuadjYuTu -& 
&  8._dp*(YuadjYuTuadjYuYu) - (4*g1p2*M1*YuadjYuYu)/5._dp - 12*g2p2*M2*YuadjYuYu -       & 
&  18*TradjYuTu*YuadjYuYu - 6._dp*(YuadjYuYuadjYuTu) - (2*YuadjYdYd*(2*g1p2*M1 +         & 
&  15._dp*(TradjYdTd) + 5._dp*(TradjYeTe) + 15*Conjg(L1II)*TL1II))/5._dp -               & 
&  18*YuadjYuYu*Conjg(L2II)*TL2II - (2*Yu*(5473*g1p4*M1 + 225*g1p2*g2p2*M1 +             & 
&  680*g1p2*g3p2*M1 + 680*g1p2*g3p2*M3 + 1800*g2p2*g3p2*M3 + 16000*g3p4*M3 +             & 
&  225*g1p2*g2p2*M2 + 12825*g2p4*M2 + 1800*g2p2*g3p2*M2 - 180*g1p2*TradjYuTu -           & 
&  3600*g3p2*TradjYuTu + 675._dp*(TrYdadjYuTuadjYd) + 675._dp*(TrYuadjYdTdadjYu) +       & 
&  180*(g1p2*M1 + 20*g3p2*M3)*TrYuadjYu + 4050._dp*(TrYuadjYuTuadjYu) + 5400*CL2IIp2*L2II*TL2II +& 
&  135*Conjg(L2II)*(L2II*(6*g1p2*M1 + 20*g2p2*M2 + 15._dp*(TradjYuTu)) + (-              & 
& 6._dp*(g1p2) - 20._dp*(g2p2) + 15._dp*(TrYuadjYu))*TL2II)))/225._dp + (18*AbsL2II*g1p2*Tu)/5._dp +& 
&  (5473*g1p4*Tu)/450._dp + 12*AbsL2II*g2p2*Tu + g1p2*g2p2*Tu + (57*g2p4*Tu)/2._dp +     & 
&  (136*g1p2*g3p2*Tu)/45._dp + 8*g2p2*g3p2*Tu + (320*g3p4*Tu)/9._dp - 12*CL2IIp2*L2IIp2*Tu  
betaTu2 =  betaTu2- 3*TrYdadjYuYuadjYd*Tu - 9*AbsL2II*TrYuadjYu*Tu + (4*g1p2*TrYuadjYu*Tu)/5._dp +        & 
&  16*g3p2*TrYuadjYu*Tu - 9*TrYuadjYuYuadjYu*Tu

 
DTu = oo16pi2*( betaTu1 + oo16pi2 * betaTu2 ) 

 
Else 
DTu = oo16pi2* betaTu1 
End If 
 
 
!-------------------- 
! Td 
!-------------------- 
 
betaTd1  = 5._dp*(TdadjYdYd) + TdadjYuYu + 8._dp*(TsIICYsIIYd) + 4._dp*(TzIIadjYzIIYd)& 
&  + 4._dp*(YdadjYdTd) + 2._dp*(YdadjYuTu) + 4._dp*(YsIICYsIITd) + 2._dp*(YzIIadjYzIITd) & 
&  + Yd*((14*g1p2*M1)/15._dp + (32*g3p2*M3)/3._dp + 6*g2p2*M2 + 6._dp*(TradjYdTd)        & 
&  + 2._dp*(TradjYeTe) + 6*Conjg(L1II)*TL1II) + 3*AbsL1II*Td - (7*g1p2*Td)               & 
& /15._dp - 3*g2p2*Td - (16*g3p2*Td)/3._dp + 3*TrYdadjYd*Td + TrYeadjYe*Td

 
 
If (TwoLoopRGE) Then 
betaTd2 = -15*AbsL1II*TdadjYdYd + (6*g1p2*TdadjYdYd)/5._dp + 12*g2p2*TdadjYdYd - 6._dp*(TdadjYdYdadjYdYd) -& 
&  4._dp*(TdadjYdYsIICYsIIYd) - 2._dp*(TdadjYdYzIIadjYzIIYd) - 3*AbsL2II*TdadjYuYu +     & 
&  (4*g1p2*TdadjYuYu)/5._dp - 4._dp*(TdadjYuYuadjYdYd) - 2._dp*(TdadjYuYuadjYuYu) -      & 
&  15*TdadjYdYd*TrYdadjYd - 5*TdadjYdYd*TrYeadjYe - 3*TdadjYuYu*TrYuadjYu -              & 
&  16._dp*(TsIICYdTpYdCYsIIYd) + (64*g1p2*TsIICYsIIYd)/15._dp + (160*g3p2*TsIICYsIIYd)/3._dp -& 
&  8*TrYsIICYsII*TsIICYsIIYd - 32._dp*(TsIICYsIIYsIICYsIIYd) - 16._dp*(TsIICYzIITpYzIICYsIIYd) -& 
&  4._dp*(TzIIadjYeYeadjYzIIYd) + (4*g1p2*TzIIadjYzIIYd)/5._dp + 12*g2p2*TzIIadjYzIIYd - & 
&  4*TrYzIIadjYzII*TzIIadjYzIIYd - 12._dp*(TzIIadjYzIIYzIIadjYzIIYd) - 12._dp*(TzIICYtIIYtIIadjYzIIYd) -& 
&  12*AbsL1II*YdadjYdTd + (6*g1p2*YdadjYdTd)/5._dp + 6*g2p2*YdadjYdTd - 12*TrYdadjYd*YdadjYdTd -& 
&  4*TrYeadjYe*YdadjYdTd - 8._dp*(YdadjYdTdadjYdYd) - 8._dp*(YdadjYdTsIICYsIIYd) -       & 
&  4._dp*(YdadjYdTzIIadjYzIIYd) - 6._dp*(YdadjYdYdadjYdTd) - 8._dp*(YdadjYdYsIICYsIITd) -& 
&  4._dp*(YdadjYdYzIIadjYzIITd) - 6*AbsL2II*YdadjYuTu + (8*g1p2*YdadjYuTu)/5._dp -       & 
&  6*TrYuadjYu*YdadjYuTu - 4._dp*(YdadjYuTuadjYdYd) - 4._dp*(YdadjYuTuadjYuYu) -         & 
&  (8*g1p2*M1*YdadjYuYu)/5._dp - 6*TradjYuTu*YdadjYuYu - 2._dp*(YdadjYuYuadjYdTd) -      & 
&  4._dp*(YdadjYuYuadjYuTu) - 12._dp*(YsIICYdTpTdCYsIIYd) - 8._dp*(YsIICYdTpYdCYsIITd) + & 
&  (32*g1p2*YsIICYsIITd)/15._dp + (80*g3p2*YsIICYsIITd)/3._dp - 4*TrYsIICYsII*YsIICYsIITd -& 
&  32._dp*(YsIICYsIITsIICYsIIYd) - (64*g1p2*M1*YsIICYsIIYd)/15._dp - (160*g3p2*M3*YsIICYsIIYd)/3._dp -& 
&  8*TrCYsIITsII*YsIICYsIIYd - 16._dp*(YsIICYsIIYsIICYsIITd) - 12._dp*(YsIICYzIITpTzIICYsIIYd) -& 
&  8._dp*(YsIICYzIITpYzIICYsIITd) - 4._dp*(YsIITdadjYdCYsIIYd) - 4._dp*(YsIITzIIadjYzIICYsIIYd) -& 
&  4._dp*(YzIIadjYeTeadjYzIIYd) - 2._dp*(YzIIadjYeYeadjYzIITd) + (2*g1p2*YzIIadjYzIITd)/5._dp +& 
&  6*g2p2*YzIIadjYzIITd - 2*TrYzIIadjYzII*YzIIadjYzIITd - 12._dp*(YzIIadjYzIITzIIadjYzIIYd)  
betaTd2 =  betaTd2- (4*g1p2*M1*YzIIadjYzIIYd)/5._dp - 12*g2p2*M2*YzIIadjYzIIYd - 4*TradjYzIITzII*YzIIadjYzIIYd -& 
&  6._dp*(YzIIadjYzIIYzIIadjYzIITd) - 12._dp*(YzIICYtIITtIIadjYzIIYd) - 6._dp*(YzIICYtIIYtIIadjYzIITd) -& 
&  (2*YdadjYdYd*(4*g1p2*M1 + 30*g2p2*M2 + 45._dp*(TradjYdTd) + 15._dp*(TradjYeTe) +      & 
&  45*Conjg(L1II)*TL1II))/5._dp - (2*Yd*(581*g1p4*M1 + 45*g1p2*g2p2*M1 + 40*g1p2*g3p2*M1 +& 
&  40*g1p2*g3p2*M3 + 360*g2p2*g3p2*M3 + 3200*g3p4*M3 + 45*g1p2*g2p2*M2 + 2565*g2p4*M2 +  & 
&  360*g2p2*g3p2*M2 + 18*g1p2*TradjYdTd - 720*g3p2*TradjYdTd - 54*g1p2*TradjYeTe -       & 
&  18*(g1p2*M1 - 40*g3p2*M3)*TrYdadjYd + 810._dp*(TrYdadjYdTdadjYd) + 540._dp*(TrYdadjYdTsIICYsII) +& 
&  270._dp*(TrYdadjYdTzIIadjYzII) + 135._dp*(TrYdadjYuTuadjYd) + 54*g1p2*M1*TrYeadjYe +  & 
&  270._dp*(TrYeadjYeTeadjYe) + 135._dp*(TrYeadjYzIITzIIadjYe) + 135._dp*(TrYeCYtIITtIIadjYe) +& 
&  540._dp*(TrYsIICYsIITdadjYd) + 135._dp*(TrYtIIadjYeTeCYtII) + 135._dp*(TrYuadjYdTdadjYu) +& 
&  135._dp*(TrYzIIadjYeTeadjYzII) + 270._dp*(TrYzIIadjYzIITdadjYd) + 1080*CL1IIp2*L1II*TL1II +& 
&  27*Conjg(L1II)*(L1II*(6*g1p2*M1 + 20*g2p2*M2 + 15._dp*(TradjYdTd) + 5._dp*(TradjYeTe) +& 
&  5._dp*(TrCYtIITtII)) + (-6._dp*(g1p2) - 20._dp*(g2p2) + 15._dp*(TrYdadjYd) +          & 
&  5._dp*(TrYeadjYe) + 5._dp*(TrYtIICYtII))*TL1II)))/45._dp - 6*YdadjYuYu*Conjg(L2II)*TL2II +& 
&  (18*AbsL1II*g1p2*Td)/5._dp + (581*g1p4*Td)/90._dp + 12*AbsL1II*g2p2*Td +              & 
&  g1p2*g2p2*Td + (57*g2p4*Td)/2._dp + (8*g1p2*g3p2*Td)/9._dp + 8*g2p2*g3p2*Td +         & 
&  (320*g3p4*Td)/9._dp - 12*CL1IIp2*L1IIp2*Td - 9*AbsL1II*TrYdadjYd*Td - (2*g1p2*TrYdadjYd*Td)/5._dp +& 
&  16*g3p2*TrYdadjYd*Td - 9*TrYdadjYdYdadjYd*Td - 12*TrYdadjYdYsIICYsII*Td -             & 
&  6*TrYdadjYdYzIIadjYzII*Td - 3*TrYdadjYuYuadjYd*Td - 3*AbsL1II*TrYeadjYe*Td +          & 
&  (6*g1p2*TrYeadjYe*Td)/5._dp - 3*TrYeadjYeYeadjYe*Td - 3*TrYeadjYzIIYzIIadjYe*Td -     & 
&  3*TrYeCYtIIYtIIadjYe*Td - 3*AbsL1II*TrYtIICYtII*Td

 
DTd = oo16pi2*( betaTd1 + oo16pi2 * betaTd2 ) 

 
Else 
DTd = oo16pi2* betaTd1 
End If 
 
 
!-------------------- 
! Te 
!-------------------- 
 
betaTe1  = 5._dp*(TeadjYeYe) + 3._dp*(TeadjYzIIYzII) + 3._dp*(TeCYtIIYtII)            & 
&  + 4._dp*(YeadjYeTe) + 6._dp*(YeadjYzIITzII) + 6._dp*(YeCYtIITtII) + Ye*((18*g1p2*M1)  & 
& /5._dp + 6*g2p2*M2 + 6._dp*(TradjYdTd) + 2._dp*(TradjYeTe) + 6*Conjg(L1II)             & 
& *TL1II) + 3*AbsL1II*Te - (9*g1p2*Te)/5._dp - 3*g2p2*Te + 3*TrYdadjYd*Te +              & 
&  TrYeadjYe*Te

 
 
If (TwoLoopRGE) Then 
betaTe2 = (-150*AbsL1II*TeadjYeYe - 12*g1p2*TeadjYeYe + 120*g2p2*TeadjYeYe - 60._dp*(TeadjYeYeadjYeYe) -& 
&  60._dp*(TeadjYzIIYdadjYdYzII) - 120._dp*(TeadjYzIIYsIICYsIIYzII) - 4*g1p2*TeadjYzIIYzII +& 
&  160*g3p2*TeadjYzIIYzII - 120._dp*(TeadjYzIIYzIIadjYeYe) - 60._dp*(TeadjYzIIYzIIadjYzIIYzII) -& 
&  30._dp*(TeCYtIITpYeCYeYtII) - 90._dp*(TeCYtIITpYzIICYzIIYtII) - 30*AbsL1II*TeCYtIIYtII +& 
&  36*g1p2*TeCYtIIYtII + 120*g2p2*TeCYtIIYtII - 120._dp*(TeCYtIIYtIIadjYeYe) -           & 
&  90._dp*(TeCYtIIYtIICYtIIYtII) - 150*TeadjYeYe*TrYdadjYd - 50*TeadjYeYe*TrYeadjYe -    & 
&  30*TeCYtIIYtII*TrYtIICYtII - 30*TeadjYzIIYzII*TrYzIIadjYzII - 120*AbsL1II*YeadjYeTe + & 
&  12*g1p2*YeadjYeTe + 60*g2p2*YeadjYeTe - 120*TrYdadjYd*YeadjYeTe - 40*TrYeadjYe*YeadjYeTe -& 
&  80._dp*(YeadjYeTeadjYeYe) - 60._dp*(YeadjYeYeadjYeTe) - 120._dp*(YeadjYzIITdadjYdYzII) -& 
&  240._dp*(YeadjYzIITsIICYsIIYzII) - 8*g1p2*YeadjYzIITzII + 320*g3p2*YeadjYzIITzII -    & 
&  60*TrYzIIadjYzII*YeadjYzIITzII - 120._dp*(YeadjYzIITzIIadjYeYe) - 120._dp*(YeadjYzIITzIIadjYzIIYzII) -& 
&  120._dp*(YeadjYzIIYdadjYdTzII) - 240._dp*(YeadjYzIIYsIICYsIITzII) + 8*g1p2*M1*YeadjYzIIYzII -& 
&  320*g3p2*M3*YeadjYzIIYzII - 60*TradjYzIITzII*YeadjYzIIYzII - 60._dp*(YeadjYzIIYzIIadjYeTe) -& 
&  120._dp*(YeadjYzIIYzIIadjYzIITzII) - 60._dp*(YeCYtIITpTeCYeYtII) - 180._dp*(YeCYtIITpTzIICYzIIYtII) -& 
&  60._dp*(YeCYtIITpYeCYeTtII) - 180._dp*(YeCYtIITpYzIICYzIITtII) - 60*AbsL1II*YeCYtIITtII +& 
&  72*g1p2*YeCYtIITtII + 240*g2p2*YeCYtIITtII - 60*TrYtIICYtII*YeCYtIITtII -             & 
&  120._dp*(YeCYtIITtIIadjYeYe) - 180._dp*(YeCYtIITtIICYtIIYtII) - 72*g1p2*M1*YeCYtIIYtII -& 
&  240*g2p2*M2*YeCYtIIYtII - 60*TrCYtIITtII*YeCYtIIYtII - 60._dp*(YeCYtIIYtIIadjYeTe) -  & 
&  180._dp*(YeCYtIIYtIICYtIITtII) - 60*YeCYtIIYtII*Conjg(L1II)*TL1II - 60*YeadjYeYe*(2*g2p2*M2 +& 
&  3._dp*(TradjYdTd) + TradjYeTe + 3*Conjg(L1II)*TL1II) - 4*Ye*(261*g1p4*M1 +            & 
&  9*g1p2*g2p2*M1 + 9*g1p2*g2p2*M2 + 285*g2p4*M2 + 2*g1p2*TradjYdTd - 80*g3p2*TradjYdTd -& 
&  6*g1p2*TradjYeTe + (-2*g1p2*M1 + 80*g3p2*M3)*TrYdadjYd + 90._dp*(TrYdadjYdTdadjYd) +  & 
&  60._dp*(TrYdadjYdTsIICYsII) + 30._dp*(TrYdadjYdTzIIadjYzII) + 15._dp*(TrYdadjYuTuadjYd) +& 
&  6*g1p2*M1*TrYeadjYe + 30._dp*(TrYeadjYeTeadjYe) + 15._dp*(TrYeadjYzIITzIIadjYe) +     & 
&  15._dp*(TrYeCYtIITtIIadjYe) + 60._dp*(TrYsIICYsIITdadjYd) + 15._dp*(TrYtIIadjYeTeCYtII) +& 
&  15._dp*(TrYuadjYdTdadjYu) + 15._dp*(TrYzIIadjYeTeadjYzII) + 30._dp*(TrYzIIadjYzIITdadjYd) +& 
&  120*CL1IIp2*L1II*TL1II + 3*Conjg(L1II)*(L1II*(6*g1p2*M1 + 20*g2p2*M2 + 15._dp*(TradjYdTd) +& 
&  5._dp*(TradjYeTe) + 5._dp*(TrCYtIITtII)) + (-6._dp*(g1p2) - 20._dp*(g2p2) +           & 
&  15._dp*(TrYdadjYd) + 5._dp*(TrYeadjYe) + 5._dp*(TrYtIICYtII))*TL1II)) +               & 
&  36*AbsL1II*g1p2*Te + 261*g1p4*Te + 120*AbsL1II*g2p2*Te + 18*g1p2*g2p2*Te +            & 
&  285*g2p4*Te - 120*CL1IIp2*L1IIp2*Te - 90*AbsL1II*TrYdadjYd*Te - 4*g1p2*TrYdadjYd*Te + & 
&  160*g3p2*TrYdadjYd*Te - 90*TrYdadjYdYdadjYd*Te - 120*TrYdadjYdYsIICYsII*Te -          & 
&  60*TrYdadjYdYzIIadjYzII*Te - 30*TrYdadjYuYuadjYd*Te - 30*AbsL1II*TrYeadjYe*Te +       & 
&  12*g1p2*TrYeadjYe*Te - 30*TrYeadjYeYeadjYe*Te - 30*TrYeadjYzIIYzIIadjYe*Te -          & 
&  30*TrYeCYtIIYtIIadjYe*Te - 30*AbsL1II*TrYtIICYtII*Te)/10._dp

 
DTe = oo16pi2*( betaTe1 + oo16pi2 * betaTe2 ) 

 
Else 
DTe = oo16pi2* betaTe1 
End If 
 
 
!-------------------- 
! TtII 
!-------------------- 
 
betaTtII1  = 2._dp*(TpTeCYeYtII) + 6._dp*(TpTzIICYzIIYtII) + TpYeCYeTtII +            & 
&  3._dp*(TpYzIICYzIITtII) + TtIIadjYeYe + 3._dp*(TtIIadjYzIIYzII) + 9._dp*(TtIICYtIIYtII)& 
&  + 2._dp*(YtIIadjYeTe) + 6._dp*(YtIIadjYzIITzII) + 9._dp*(YtIICYtIITtII)               & 
&  + (2*YtII*(9*g1p2*M1 + 35*g2p2*M2 + 5._dp*(TrCYtIITtII) + 5*Conjg(L1II)               & 
& *TL1II))/5._dp + AbsL1II*TtII - (9*g1p2*TtII)/5._dp - 7*g2p2*TtII + TrYtIICYtII*TtII

 
 
If (TwoLoopRGE) Then 
betaTtII2 = (-40._dp*(TpTeCYeTpYeCYeYtII) - 60*AbsL1II*TpTeCYeYtII + 24*g1p2*TpTeCYeYtII -        & 
&  120._dp*(TpTzIICYdTpYdCYzIIYtII) - 240._dp*(TpTzIICYsIIYsIICYzIIYtII) -               & 
&  120._dp*(TpTzIICYzIITpYzIICYzIIYtII) - 8*g1p2*TpTzIICYzIIYtII + 320*g3p2*TpTzIICYzIIYtII -& 
&  40._dp*(TpYeCYeTpTeCYeYtII) - 20._dp*(TpYeCYeTpYeCYeTtII) - 30*AbsL1II*TpYeCYeTtII +  & 
&  12*g1p2*TpYeCYeTtII - 24*g1p2*M1*TpYeCYeYtII - 120._dp*(TpYzIICYdTpTdCYzIIYtII) -     & 
&  60._dp*(TpYzIICYdTpYdCYzIITtII) - 240._dp*(TpYzIICYsIITsIICYzIIYtII) - 120._dp*(TpYzIICYsIIYsIICYzIITtII) -& 
&  120._dp*(TpYzIICYzIITpTzIICYzIIYtII) - 60._dp*(TpYzIICYzIITpYzIICYzIITtII) -          & 
&  4*g1p2*TpYzIICYzIITtII + 160*g3p2*TpYzIICYzIITtII + 8*g1p2*M1*TpYzIICYzIIYtII -       & 
&  320*g3p2*M3*TpYzIICYzIIYtII - 60*TpYeCYeYtII*TradjYdTd - 20*TpYeCYeYtII*TradjYeTe -   & 
&  60*TpYzIICYzIIYtII*TradjYzIITzII - 60*TpTeCYeYtII*TrYdadjYd - 30*TpYeCYeTtII*TrYdadjYd -& 
&  20*TpTeCYeYtII*TrYeadjYe - 10*TpYeCYeTtII*TrYeadjYe - 60*TpTzIICYzIIYtII*TrYzIIadjYzII -& 
&  30*TpYzIICYzIITtII*TrYzIIadjYzII - 30*AbsL1II*TtIIadjYeYe + 12*g1p2*TtIIadjYeYe -     & 
&  30*TrYdadjYd*TtIIadjYeYe - 10*TrYeadjYe*TtIIadjYeYe - 20._dp*(TtIIadjYeYeadjYeYe) -   & 
&  60._dp*(TtIIadjYeYeCYtIIYtII) - 60._dp*(TtIIadjYzIIYdadjYdYzII) - 120._dp*(TtIIadjYzIIYsIICYsIIYzII) -& 
&  4*g1p2*TtIIadjYzIIYzII + 160*g3p2*TtIIadjYzIIYzII - 30*TrYzIIadjYzII*TtIIadjYzIIYzII -& 
&  60._dp*(TtIIadjYzIIYzIIadjYzIIYzII) - 180._dp*(TtIIadjYzIIYzIICYtIIYtII) -            & 
&  30._dp*(TtIICYtIITpYeCYeYtII) - 90._dp*(TtIICYtIITpYzIICYzIIYtII) - 90*AbsL1II*TtIICYtIIYtII +& 
&  108*g1p2*TtIICYtIIYtII + 360*g2p2*TtIICYtIIYtII - 90*TrYtIICYtII*TtIICYtIIYtII -      & 
&  270._dp*(TtIICYtIIYtIICYtIIYtII) - 60*AbsL1II*YtIIadjYeTe + 24*g1p2*YtIIadjYeTe -     & 
&  60*TrYdadjYd*YtIIadjYeTe - 20*TrYeadjYe*YtIIadjYeTe - 40._dp*(YtIIadjYeTeadjYeYe) -   & 
&  60._dp*(YtIIadjYeTeCYtIIYtII) - 40._dp*(YtIIadjYeYeadjYeTe) - 30._dp*(YtIIadjYeYeCYtIITtII) -& 
&  120._dp*(YtIIadjYzIITdadjYdYzII) - 240._dp*(YtIIadjYzIITsIICYsIIYzII) -               & 
&  8*g1p2*YtIIadjYzIITzII + 320*g3p2*YtIIadjYzIITzII - 60*TrYzIIadjYzII*YtIIadjYzIITzII -& 
&  120._dp*(YtIIadjYzIITzIIadjYzIIYzII) - 180._dp*(YtIIadjYzIITzIICYtIIYtII) -           & 
&  120._dp*(YtIIadjYzIIYdadjYdTzII) - 240._dp*(YtIIadjYzIIYsIICYsIITzII) +               & 
&  8*g1p2*M1*YtIIadjYzIIYzII - 320*g3p2*M3*YtIIadjYzIIYzII - 60*TradjYzIITzII*YtIIadjYzIIYzII -& 
&  120._dp*(YtIIadjYzIIYzIIadjYzIITzII) - 90._dp*(YtIIadjYzIIYzIICYtIITtII) -            & 
&  60._dp*(YtIICYtIITpTeCYeYtII) - 180._dp*(YtIICYtIITpTzIICYzIIYtII) - 60._dp*(YtIICYtIITpYeCYeTtII) -& 
&  180._dp*(YtIICYtIITpYzIICYzIITtII) - 90*AbsL1II*YtIICYtIITtII + 108*g1p2*YtIICYtIITtII +& 
&  360*g2p2*YtIICYtIITtII - 90*TrYtIICYtII*YtIICYtIITtII - 360._dp*(YtIICYtIITtIICYtIIYtII) -& 
&  144*g1p2*M1*YtIICYtIIYtII - 480*g2p2*M2*YtIICYtIIYtII - 120*TrCYtIITtII*YtIICYtIIYtII -& 
&  270._dp*(YtIICYtIIYtIICYtIITtII) - 60*TpYeCYeYtII*Conjg(L1II)*TL1II - 120*YtIICYtIIYtII*Conjg(L1II)*TL1II -& 
&  4*YtIIadjYeYe*(6*g1p2*M1 + 15._dp*(TradjYdTd) + 5._dp*(TradjYeTe) + 15*Conjg(L1II)*TL1II) -& 
&  4*YtII*(261*g1p4*M1 + 57*g1p2*g2p2*M1 + 57*g1p2*g2p2*M2 + 765*g2p4*M2 +               & 
&  3*g1p2*TrCYtIITtII + 5*g2p2*TrCYtIITtII + 10._dp*(TrYeCYtIITtIIadjYe) +               & 
&  10._dp*(TrYtIIadjYeTeCYtII) + 30._dp*(TrYtIIadjYzIITzIICYtII) - (3*g1p2*M1 +          & 
&  5*g2p2*M2)*TrYtIICYtII + 60._dp*(TrYtIICYtIITtIICYtII) + 30._dp*(TrYzIICYtIITtIIadjYzII) +& 
&  60*CL1IIp2*L1II*TL1II + Conjg(L1II)*(L1II*(-3*g1p2*M1 - 5*g2p2*M2 + 30._dp*(TradjYdTd) +& 
&  10._dp*(TradjYeTe)) + (3._dp*(g1p2) + 5._dp*(g2p2) + 30._dp*(TrYdadjYd) +             & 
&  10._dp*(TrYeadjYe))*TL1II)) - 6*AbsL1II*g1p2*TtII + 261*g1p4*TtII - 10*AbsL1II*g2p2*TtII +& 
&  114*g1p2*g2p2*TtII + 765*g2p4*TtII - 60*CL1IIp2*L1IIp2*TtII - 60*AbsL1II*TrYdadjYd*TtII -& 
&  20*AbsL1II*TrYeadjYe*TtII - 20*TrYeCYtIIYtIIadjYe*TtII - 60*TrYtIIadjYzIIYzIICYtII*TtII -& 
&  6*g1p2*TrYtIICYtII*TtII - 10*g2p2*TrYtIICYtII*TtII - 60*TrYtIICYtIIYtIICYtII*TtII)/10._dp

 
DTtII = oo16pi2*( betaTtII1 + oo16pi2 * betaTtII2 ) 

 
Else 
DTtII = oo16pi2* betaTtII1 
End If 
 
 
!-------------------- 
! TsII 
!-------------------- 
 
betaTsII1  = 4._dp*(TdadjYdYsII) + 2._dp*(TsIICYdTpYd) + 12._dp*(TsIICYsIIYsII)       & 
&  + 2._dp*(TsIICYzIITpYzII) + 4._dp*(TzIIadjYzIIYsII) + 2._dp*(YdadjYdTsII)             & 
&  + ((8*g1p2*M1)/5._dp + 24*g3p2*M3 + 2._dp*(TrCYsIITsII))*YsII + 4._dp*(YsIICYdTpTd)   & 
&  + 12._dp*(YsIICYsIITsII) + 4._dp*(YsIICYzIITpTzII) + 2._dp*(YzIIadjYzIITsII)          & 
&  - (4*g1p2*TsII)/5._dp - 12*g3p2*TsII + TrYsIICYsII*TsII

 
 
If (TwoLoopRGE) Then 
betaTsII2 = (-2*(30._dp*(TdadjYdYdadjYdYsII) + 90*AbsL1II*TdadjYdYsII - 6*g1p2*TdadjYdYsII -      & 
&  90*g2p2*TdadjYdYsII + 30._dp*(TdadjYuYuadjYdYsII) + 90*TdadjYdYsII*TrYdadjYd +        & 
&  30*TdadjYdYsII*TrYeadjYe + 45*AbsL1II*TsIICYdTpYd - 3*g1p2*TsIICYdTpYd -              & 
&  45*g2p2*TsIICYdTpYd + 45*TrYdadjYd*TsIICYdTpYd + 15*TrYeadjYe*TsIICYdTpYd +           & 
&  15._dp*(TsIICYdTpYdCYdTpYd) + 120._dp*(TsIICYdTpYdCYsIIYsII) + 15._dp*(TsIICYdTpYuCYuTpYd) +& 
&  60._dp*(TsIICYsIIYdadjYdYsII) - 48*g1p2*TsIICYsIIYsII - 600*g3p2*TsIICYsIIYsII +      & 
&  90*TrYsIICYsII*TsIICYsIIYsII + 360._dp*(TsIICYsIIYsIICYsIIYsII) + 60._dp*(TsIICYsIIYzIIadjYzIIYsII) +& 
&  15._dp*(TsIICYzIITpYeCYeTpYzII) - 3*g1p2*TsIICYzIITpYzII - 45*g2p2*TsIICYzIITpYzII +  & 
&  15*TrYzIIadjYzII*TsIICYzIITpYzII + 120._dp*(TsIICYzIITpYzIICYsIIYsII) +               & 
&  45._dp*(TsIICYzIITpYzIICYzIITpYzII) + 45._dp*(TsIICYzIIYtIICYtIITpYzII) +             & 
&  30._dp*(TzIIadjYeYeadjYzIIYsII) - 6*g1p2*TzIIadjYzIIYsII - 90*g2p2*TzIIadjYzIIYsII +  & 
&  30*TrYzIIadjYzII*TzIIadjYzIIYsII + 90._dp*(TzIIadjYzIIYzIIadjYzIIYsII) +              & 
&  90._dp*(TzIICYtIIYtIIadjYzIIYsII) + 30._dp*(YdadjYdTdadjYdYsII) + 45*AbsL1II*YdadjYdTsII -& 
&  3*g1p2*YdadjYdTsII - 45*g2p2*YdadjYdTsII + 45*TrYdadjYd*YdadjYdTsII + 15*TrYeadjYe*YdadjYdTsII +& 
&  15._dp*(YdadjYdYdadjYdTsII) + 30._dp*(YdadjYuTuadjYdYsII) + 15._dp*(YdadjYuYuadjYdTsII) +& 
&  4*(84*g1p4*M1 + 32*g1p2*g3p2*M1 + 32*g1p2*g3p2*M3 + 800*g3p4*M3 + (g1p2 +             & 
&  5._dp*(g3p2))*TrCYsIITsII + 15._dp*(TrYdadjYdTsIICYsII) - (g1p2*M1 + 5*g3p2*M3)*TrYsIICYsII +& 
&  15._dp*(TrYsIICYsIITdadjYd) + 60._dp*(TrYsIICYsIITsIICYsII) + 15._dp*(TrYsIICYsIITzIIadjYzII) +& 
&  15._dp*(TrYzIIadjYzIITsIICYsII))*YsII + 90*AbsL1II*YsIICYdTpTd - 6*g1p2*YsIICYdTpTd - & 
&  90*g2p2*YsIICYdTpTd + 90*TrYdadjYd*YsIICYdTpTd + 30*TrYeadjYe*YsIICYdTpTd +           & 
&  30._dp*(YsIICYdTpTdCYdTpYd) + 120._dp*(YsIICYdTpTdCYsIIYsII) + 30._dp*(YsIICYdTpTuCYuTpYd) +& 
&  6*g1p2*M1*YsIICYdTpYd + 90*g2p2*M2*YsIICYdTpYd + 90*TradjYdTd*YsIICYdTpYd +           & 
&  30*TradjYeTe*YsIICYdTpYd + 30._dp*(YsIICYdTpYdCYdTpTd) + 60._dp*(YsIICYdTpYdCYsIITsII) +& 
&  30._dp*(YsIICYdTpYuCYuTpTd) + 120._dp*(YsIICYsIITdadjYdYsII) - 48*g1p2*YsIICYsIITsII -& 
&  600*g3p2*YsIICYsIITsII + 90*TrYsIICYsII*YsIICYsIITsII + 480._dp*(YsIICYsIITsIICYsIIYsII) +& 
&  120._dp*(YsIICYsIITzIIadjYzIIYsII) + 120._dp*(YsIICYsIIYdadjYdTsII) + 64*g1p2*M1*YsIICYsIIYsII +& 
&  800*g3p2*M3*YsIICYsIIYsII + 120*TrCYsIITsII*YsIICYsIIYsII + 360._dp*(YsIICYsIIYsIICYsIITsII) +& 
&  120._dp*(YsIICYsIIYzIIadjYzIITsII) + 30._dp*(YsIICYzIITpTeCYeTpYzII) - 6*g1p2*YsIICYzIITpTzII -& 
&  90*g2p2*YsIICYzIITpTzII + 30*TrYzIIadjYzII*YsIICYzIITpTzII + 120._dp*(YsIICYzIITpTzIICYsIIYsII) +& 
&  90._dp*(YsIICYzIITpTzIICYzIITpYzII) + 30._dp*(YsIICYzIITpYeCYeTpTzII) +               & 
&  6*g1p2*M1*YsIICYzIITpYzII + 90*g2p2*M2*YsIICYzIITpYzII + 30*TradjYzIITzII*YsIICYzIITpYzII +& 
&  60._dp*(YsIICYzIITpYzIICYsIITsII) + 90._dp*(YsIICYzIITpYzIICYzIITpTzII) +             & 
&  90._dp*(YsIICYzIITtIICYtIITpYzII) + 90._dp*(YsIICYzIIYtIICYtIITpTzII) +               & 
&  30._dp*(YzIIadjYeTeadjYzIIYsII) + 15._dp*(YzIIadjYeYeadjYzIITsII) - 3*g1p2*YzIIadjYzIITsII -& 
&  45*g2p2*YzIIadjYzIITsII + 15*TrYzIIadjYzII*YzIIadjYzIITsII + 90._dp*(YzIIadjYzIITzIIadjYzIIYsII) +& 
&  6*g1p2*M1*YzIIadjYzIIYsII + 90*g2p2*M2*YzIIadjYzIIYsII + 30*TradjYzIITzII*YzIIadjYzIIYsII +& 
&  45._dp*(YzIIadjYzIIYzIIadjYzIITsII) + 90._dp*(YzIICYtIITtIIadjYzIIYsII) +             & 
&  45._dp*(YzIICYtIIYtIIadjYzIITsII) + 90*YsIICYdTpYd*Conjg(L1II)*TL1II + 6*YdadjYdYsII*(g1p2*M1 +& 
&  15*g2p2*M2 + 15._dp*(TradjYdTd) + 5._dp*(TradjYeTe) + 15*Conjg(L1II)*TL1II) -         & 
&  84*g1p4*TsII - 64*g1p2*g3p2*TsII - 800*g3p4*TsII + 30*TrYdadjYdYsIICYsII*TsII +       & 
&  2*g1p2*TrYsIICYsII*TsII + 10*g3p2*TrYsIICYsII*TsII + 60*TrYsIICYsIIYsIICYsII*TsII +   & 
&  30*TrYsIICYsIIYzIIadjYzII*TsII))/15._dp

 
DTsII = oo16pi2*( betaTsII1 + oo16pi2 * betaTsII2 ) 

 
Else 
DTsII = oo16pi2* betaTsII1 
End If 
 
 
!-------------------- 
! TzII 
!-------------------- 
 
betaTzII1  = 4._dp*(TdadjYdYzII) + 8._dp*(TsIICYsIIYzII) + TzIIadjYeYe +              & 
&  7._dp*(TzIIadjYzIIYzII) + 3._dp*(TzIICYtIIYtII) + 2._dp*(YdadjYdTzII) +               & 
&  4._dp*(YsIICYsIITzII) + (2*(7*g1p2*M1 + 80*g3p2*M3 + 45*g2p2*M2 + 15._dp*(TradjYzIITzII))& 
& *YzII)/15._dp + 2._dp*(YzIIadjYeTe) + 8._dp*(YzIIadjYzIITzII) + 6._dp*(YzIICYtIITtII)  & 
&  - (7*g1p2*TzII)/15._dp - 3*g2p2*TzII - (16*g3p2*TzII)/3._dp + TrYzIIadjYzII*TzII

 
 
If (TwoLoopRGE) Then 
betaTzII2 = -4._dp*(TdadjYdYdadjYdYzII) - 12*AbsL1II*TdadjYdYzII + (4*g1p2*TdadjYdYzII)/5._dp +   & 
&  12*g2p2*TdadjYdYzII - 4._dp*(TdadjYuYuadjYdYzII) - 12*TdadjYdYzII*TrYdadjYd -         & 
&  4*TdadjYdYzII*TrYeadjYe - 16._dp*(TsIICYdTpYdCYsIIYzII) - 32._dp*(TsIICYsIIYsIICYsIIYzII) +& 
&  (64*g1p2*TsIICYsIIYzII)/15._dp + (160*g3p2*TsIICYsIIYzII)/3._dp - 8*TrYsIICYsII*TsIICYsIIYzII -& 
&  16._dp*(TsIICYzIITpYzIICYsIIYzII) - 3*AbsL1II*TzIIadjYeYe + (6*g1p2*TzIIadjYeYe)/5._dp -& 
&  3*TrYdadjYd*TzIIadjYeYe - TrYeadjYe*TzIIadjYeYe - 2._dp*(TzIIadjYeYeadjYeYe) -        & 
&  4._dp*(TzIIadjYeYeadjYzIIYzII) - 6._dp*(TzIIadjYzIIYdadjYdYzII) - 12._dp*(TzIIadjYzIIYsIICYsIIYzII) +& 
&  (2*g1p2*TzIIadjYzIIYzII)/5._dp + 12*g2p2*TzIIadjYzIIYzII + 16*g3p2*TzIIadjYzIIYzII -  & 
&  7*TrYzIIadjYzII*TzIIadjYzIIYzII - 18._dp*(TzIIadjYzIIYzIIadjYzIIYzII) -               & 
&  3._dp*(TzIICYtIITpYeCYeYtII) - 9._dp*(TzIICYtIITpYzIICYzIIYtII) - 3*AbsL1II*TzIICYtIIYtII +& 
&  (18*g1p2*TzIICYtIIYtII)/5._dp + 12*g2p2*TzIICYtIIYtII - 3*TrYtIICYtII*TzIICYtIIYtII - & 
&  12._dp*(TzIICYtIIYtIIadjYzIIYzII) - 9._dp*(TzIICYtIIYtIICYtIIYtII) - 4._dp*(YdadjYdTdadjYdYzII) -& 
&  6*AbsL1II*YdadjYdTzII + (2*g1p2*YdadjYdTzII)/5._dp + 6*g2p2*YdadjYdTzII -             & 
&  6*TrYdadjYd*YdadjYdTzII - 2*TrYeadjYe*YdadjYdTzII - 2._dp*(YdadjYdYdadjYdTzII) -      & 
&  4._dp*(YdadjYuTuadjYdYzII) - 2._dp*(YdadjYuYuadjYdTzII) - 12._dp*(YsIICYdTpTdCYsIIYzII) -& 
&  8._dp*(YsIICYdTpYdCYsIITzII) - 32._dp*(YsIICYsIITsIICYsIIYzII) + (32*g1p2*YsIICYsIITzII)/15._dp +& 
&  (80*g3p2*YsIICYsIITzII)/3._dp - 4*TrYsIICYsII*YsIICYsIITzII - 16._dp*(YsIICYsIIYsIICYsIITzII) -& 
&  (64*g1p2*M1*YsIICYsIIYzII)/15._dp - (160*g3p2*M3*YsIICYsIIYzII)/3._dp -               & 
&  8*TrCYsIITsII*YsIICYsIIYzII - 12._dp*(YsIICYzIITpTzIICYsIIYzII) - 8._dp*(YsIICYzIITpYzIICYsIITzII) -& 
&  4._dp*(YsIITdadjYdCYsIIYzII) - 4._dp*(YsIITzIIadjYzIICYsIIYzII) - (2*(581*g1p4*M1 +   & 
&  45*g1p2*g2p2*M1 + 40*g1p2*g3p2*M1 + 40*g1p2*g3p2*M3 + 360*g2p2*g3p2*M3 +              & 
&  3200*g3p4*M3 + 45*g1p2*g2p2*M2 + 2565*g2p4*M2 + 360*g2p2*g3p2*M2 - 18*g1p2*TradjYzIITzII +& 
&  90._dp*(TrYdadjYdTzIIadjYzII) + 45._dp*(TrYeadjYzIITzIIadjYe) + 180._dp*(TrYsIICYsIITzIIadjYzII) +& 
&  135._dp*(TrYtIIadjYzIITzIICYtII) + 45._dp*(TrYzIIadjYeTeadjYzII) + 18*g1p2*M1*TrYzIIadjYzII +& 
&  90._dp*(TrYzIIadjYzIITdadjYd) + 180._dp*(TrYzIIadjYzIITsIICYsII) + 450._dp*(TrYzIIadjYzIITzIIadjYzII) +& 
&  135._dp*(TrYzIICYtIITtIIadjYzII))*YzII)/45._dp - 6*AbsL1II*YzIIadjYeTe  
betaTzII2 =  betaTzII2+ (12*g1p2*YzIIadjYeTe)/5._dp - 6*TrYdadjYd*YzIIadjYeTe - 2*TrYeadjYe*YzIIadjYeTe -     & 
&  4._dp*(YzIIadjYeTeadjYeYe) - 4._dp*(YzIIadjYeTeadjYzIIYzII) - (12*g1p2*M1*YzIIadjYeYe)/5._dp -& 
&  6*TradjYdTd*YzIIadjYeYe - 2*TradjYeTe*YzIIadjYeYe - 4._dp*(YzIIadjYeYeadjYeTe) -      & 
&  2._dp*(YzIIadjYeYeadjYzIITzII) - 12._dp*(YzIIadjYzIITdadjYdYzII) - 24._dp*(YzIIadjYzIITsIICYsIIYzII) -& 
&  (2*g1p2*YzIIadjYzIITzII)/5._dp + 6*g2p2*YzIIadjYzIITzII + 32*g3p2*YzIIadjYzIITzII -   & 
&  8*TrYzIIadjYzII*YzIIadjYzIITzII - 24._dp*(YzIIadjYzIITzIIadjYzIIYzII) -               & 
&  12._dp*(YzIIadjYzIIYdadjYdTzII) - 24._dp*(YzIIadjYzIIYsIICYsIITzII) - 32*g3p2*M3*YzIIadjYzIIYzII -& 
&  12*g2p2*M2*YzIIadjYzIIYzII - 10*TradjYzIITzII*YzIIadjYzIIYzII - 18._dp*(YzIIadjYzIIYzIIadjYzIITzII) -& 
&  6._dp*(YzIICYtIITpTeCYeYtII) - 18._dp*(YzIICYtIITpTzIICYzIIYtII) - 6._dp*(YzIICYtIITpYeCYeTtII) -& 
&  18._dp*(YzIICYtIITpYzIICYzIITtII) - 6*AbsL1II*YzIICYtIITtII + (36*g1p2*YzIICYtIITtII)/5._dp +& 
&  24*g2p2*YzIICYtIITtII - 6*TrYtIICYtII*YzIICYtIITtII - 12._dp*(YzIICYtIITtIIadjYzIIYzII) -& 
&  18._dp*(YzIICYtIITtIICYtIIYtII) - (36*g1p2*M1*YzIICYtIIYtII)/5._dp - 24*g2p2*M2*YzIICYtIIYtII -& 
&  6*TrCYtIITtII*YzIICYtIIYtII - 6._dp*(YzIICYtIIYtIIadjYzIITzII) - 18._dp*(YzIICYtIIYtIICYtIITtII) -& 
&  6*YzIIadjYeYe*Conjg(L1II)*TL1II - 6*YzIICYtIIYtII*Conjg(L1II)*TL1II - (4*YdadjYdYzII*(g1p2*M1 +& 
&  15*g2p2*M2 + 15._dp*(TradjYdTd) + 5._dp*(TradjYeTe) + 15*Conjg(L1II)*TL1II))/5._dp +  & 
&  (581*g1p4*TzII)/90._dp + g1p2*g2p2*TzII + (57*g2p4*TzII)/2._dp + (8*g1p2*g3p2*TzII)/9._dp +& 
&  8*g2p2*g3p2*TzII + (320*g3p4*TzII)/9._dp - 2*TrYdadjYdYzIIadjYzII*TzII -              & 
&  TrYeadjYzIIYzIIadjYe*TzII - 4*TrYsIICYsIIYzIIadjYzII*TzII - 3*TrYtIIadjYzIIYzIICYtII*TzII +& 
&  (2*g1p2*TrYzIIadjYzII*TzII)/5._dp - 5*TrYzIIadjYzIIYzIIadjYzII*TzII

 
DTzII = oo16pi2*( betaTzII1 + oo16pi2 * betaTzII2 ) 

 
Else 
DTzII = oo16pi2* betaTzII1 
End If 
 
 
!-------------------- 
! TL1II 
!-------------------- 
 
betaTL1II1  = (2*L1II*(9*g1p2*M1 + 35*g2p2*M2 + 30._dp*(TradjYdTd) + 10._dp*(TradjYeTe)& 
&  + 5._dp*(TrCYtIITtII)))/5._dp + (21._dp*(AbsL1II) - 9._dp*(g1p2)/5._dp -              & 
&  7._dp*(g2p2) + 6._dp*(TrYdadjYd) + 2._dp*(TrYeadjYe) + TrYtIICYtII)*TL1II

 
 
If (TwoLoopRGE) Then 
betaTL1II2 = (-4*L1II*(261*g1p4*M1 + 57*g1p2*g2p2*M1 + 57*g1p2*g2p2*M2 + 765*g2p4*M2 +             & 
&  4*g1p2*TradjYdTd - 160*g3p2*TradjYdTd - 12*g1p2*TradjYeTe + 3*g1p2*TrCYtIITtII +      & 
&  5*g2p2*TrCYtIITtII - 4*(g1p2*M1 - 40*g3p2*M3)*TrYdadjYd + 180._dp*(TrYdadjYdTdadjYd) +& 
&  120._dp*(TrYdadjYdTsIICYsII) + 60._dp*(TrYdadjYdTzIIadjYzII) + 30._dp*(TrYdadjYuTuadjYd) +& 
&  12*g1p2*M1*TrYeadjYe + 60._dp*(TrYeadjYeTeadjYe) + 30._dp*(TrYeadjYzIITzIIadjYe) +    & 
&  40._dp*(TrYeCYtIITtIIadjYe) + 120._dp*(TrYsIICYsIITdadjYd) + 40._dp*(TrYtIIadjYeTeCYtII) +& 
&  30._dp*(TrYtIIadjYzIITzIICYtII) - 3*g1p2*M1*TrYtIICYtII - 5*g2p2*M2*TrYtIICYtII +     & 
&  60._dp*(TrYtIICYtIITtIICYtII) + 30._dp*(TrYuadjYdTdadjYu) + 30._dp*(TrYzIIadjYeTeadjYzII) +& 
&  60._dp*(TrYzIIadjYzIITdadjYd) + 30._dp*(TrYzIICYtIITtIIadjYzII)) - 1500*CL1IIp2*L1IIp2*TL1II -& 
&  (-261._dp*(g1p4) - 114*g1p2*g2p2 - 765._dp*(g2p4) + 8*(g1p2 - 40._dp*(g3p2))*TrYdadjYd +& 
&  180._dp*(TrYdadjYdYdadjYd) + 240._dp*(TrYdadjYdYsIICYsII) + 120._dp*(TrYdadjYdYzIIadjYzII) +& 
&  60._dp*(TrYdadjYuYuadjYd) - 24*g1p2*TrYeadjYe + 60._dp*(TrYeadjYeYeadjYe) +           & 
&  60._dp*(TrYeadjYzIIYzIIadjYe) + 80._dp*(TrYeCYtIIYtIIadjYe) + 60._dp*(TrYtIIadjYzIIYzIICYtII) +& 
&  6*g1p2*TrYtIICYtII + 10*g2p2*TrYtIICYtII + 60._dp*(TrYtIICYtIIYtIICYtII))*TL1II -     & 
&  2*AbsL1II*(2*L1II*(33*g1p2*M1 + 115*g2p2*M2 + 120._dp*(TradjYdTd) + 40._dp*(TradjYeTe) +& 
&  30._dp*(TrCYtIITtII)) + (-99._dp*(g1p2) - 345._dp*(g2p2) + 360._dp*(TrYdadjYd) +      & 
&  120._dp*(TrYeadjYe) + 90._dp*(TrYtIICYtII))*TL1II))/10._dp

 
DTL1II = oo16pi2*( betaTL1II1 + oo16pi2 * betaTL1II2 ) 

 
Else 
DTL1II = oo16pi2* betaTL1II1 
End If 
 
 
!-------------------- 
! TL2II 
!-------------------- 
 
betaTL2II1  = (2*L2II*(9*g1p2*M1 + 35*g2p2*M2 + 30._dp*(TradjYuTu)))/5._dp +          & 
&  (21._dp*(AbsL2II) - 9._dp*(g1p2)/5._dp - 7._dp*(g2p2) + 6._dp*(TrYuadjYu))*TL2II

 
 
If (TwoLoopRGE) Then 
betaTL2II2 = (-4*L2II*(-8*(g1p2 + 20._dp*(g3p2))*TradjYuTu + 8*(g1p2*M1 + 20*g3p2*M3)*TrYuadjYu +  & 
&  3*(87*g1p4*M1 + 19*g1p2*g2p2*M1 + 19*g1p2*g2p2*M2 + 255*g2p4*M2 + 10._dp*(TrYdadjYuTuadjYd) +& 
&  10._dp*(TrYuadjYdTdadjYu) + 60._dp*(TrYuadjYuTuadjYu))) - 1500*CL2IIp2*L2IIp2*TL2II + & 
&  (16*(g1p2 + 20._dp*(g3p2))*TrYuadjYu + 3*(87._dp*(g1p4) + 38*g1p2*g2p2 +              & 
&  255._dp*(g2p4) - 20._dp*(TrYdadjYuYuadjYd) - 60._dp*(TrYuadjYuYuadjYu)))*TL2II -      & 
&  2*AbsL2II*(2*L2II*(33*g1p2*M1 + 115*g2p2*M2 + 120._dp*(TradjYuTu)) + (-               & 
& 99._dp*(g1p2) - 345._dp*(g2p2) + 360._dp*(TrYuadjYu))*TL2II))/10._dp

 
DTL2II = oo16pi2*( betaTL2II1 + oo16pi2 * betaTL2II2 ) 

 
Else 
DTL2II = oo16pi2* betaTL2II1 
End If 
 
 
!-------------------- 
! Bmu 
!-------------------- 
 
betaBmu1  = (6*g1p2*M1*Mu)/5._dp + 6*g2p2*M2*Mu + 6*TradjYdTd*Mu + 2*TradjYeTe*Mu +   & 
&  6*TradjYuTu*Mu + (3._dp*(AbsL1II) + 3._dp*(AbsL2II) - 3._dp*(g1p2)/5._dp -            & 
&  3._dp*(g2p2) + 3._dp*(TrYdadjYd) + TrYeadjYe + 3._dp*(TrYuadjYu))*Bmu +               & 
&  6*Mu*Conjg(L1II)*TL1II + 6*Mu*Conjg(L2II)*TL2II

 
 
If (TwoLoopRGE) Then 
betaBmu2 = (417._dp*(g1p4)/50._dp + (9*g1p2*g2p2)/5._dp + 57._dp*(g2p4)/2._dp - 12*CL1IIp2*L1IIp2 -& 
&  12*CL2IIp2*L2IIp2 - (2*g1p2*TrYdadjYd)/5._dp + 16*g3p2*TrYdadjYd - 9._dp*(TrYdadjYdYdadjYd) -& 
&  12._dp*(TrYdadjYdYsIICYsII) - 6._dp*(TrYdadjYdYzIIadjYzII) - 6._dp*(TrYdadjYuYuadjYd) +& 
&  (6*g1p2*TrYeadjYe)/5._dp - 3._dp*(TrYeadjYeYeadjYe) - 3._dp*(TrYeadjYzIIYzIIadjYe) -  & 
&  3._dp*(TrYeCYtIIYtIIadjYe) + (3*AbsL1II*(6._dp*(g1p2) + 20._dp*(g2p2) -               & 
&  15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe) - 5._dp*(TrYtIICYtII)))/5._dp + (3*AbsL2II*(6._dp*(g1p2) +& 
&  20._dp*(g2p2) - 15._dp*(TrYuadjYu)))/5._dp + (4*g1p2*TrYuadjYu)/5._dp +               & 
&  16*g3p2*TrYuadjYu - 9._dp*(TrYuadjYuYuadjYu))*Bmu - (2*Mu*(417*g1p4*M1 +              & 
&  45*g1p2*g2p2*M1 + 45*g1p2*g2p2*M2 + 1425*g2p4*M2 + 10*g1p2*TradjYdTd - 400*g3p2*TradjYdTd -& 
&  30*g1p2*TradjYeTe - 20*g1p2*TradjYuTu - 400*g3p2*TradjYuTu - 10*g1p2*M1*TrYdadjYd +   & 
&  400*g3p2*M3*TrYdadjYd + 450._dp*(TrYdadjYdTdadjYd) + 300._dp*(TrYdadjYdTsIICYsII) +   & 
&  150._dp*(TrYdadjYdTzIIadjYzII) + 150._dp*(TrYdadjYuTuadjYd) + 30*g1p2*M1*TrYeadjYe +  & 
&  150._dp*(TrYeadjYeTeadjYe) + 75._dp*(TrYeadjYzIITzIIadjYe) + 75._dp*(TrYeCYtIITtIIadjYe) +& 
&  300._dp*(TrYsIICYsIITdadjYd) + 75._dp*(TrYtIIadjYeTeCYtII) + 150._dp*(TrYuadjYdTdadjYu) +& 
&  20*g1p2*M1*TrYuadjYu + 400*g3p2*M3*TrYuadjYu + 450._dp*(TrYuadjYuTuadjYu) +           & 
&  75._dp*(TrYzIIadjYeTeadjYzII) + 150._dp*(TrYzIIadjYzIITdadjYd) + 600*CL1IIp2*L1II*TL1II +& 
&  15*Conjg(L1II)*(L1II*(6*g1p2*M1 + 20*g2p2*M2 + 15._dp*(TradjYdTd) + 5._dp*(TradjYeTe) +& 
&  5._dp*(TrCYtIITtII)) + (-6._dp*(g1p2) - 20._dp*(g2p2) + 15._dp*(TrYdadjYd) +          & 
&  5._dp*(TrYeadjYe) + 5._dp*(TrYtIICYtII))*TL1II) + 600*CL2IIp2*L2II*TL2II +            & 
&  15*Conjg(L2II)*(L2II*(6*g1p2*M1 + 20*g2p2*M2 + 15._dp*(TradjYuTu)) + (-               & 
& 6._dp*(g1p2) - 20._dp*(g2p2) + 15._dp*(TrYuadjYu))*TL2II)))/25._dp

 
DBmu = oo16pi2*( betaBmu1 + oo16pi2 * betaBmu2 ) 

 
Else 
DBmu = oo16pi2* betaBmu1 
End If 
 
 
!-------------------- 
! BMTII 
!-------------------- 
!goto 123 
betaBMTII1  = (AbsL1II + AbsL2II - 12._dp*(g1p2)/5._dp - 8._dp*(g2p2) +               & 
&  TrYtIICYtII)*BMTII + (2*MTII*(12*g1p2*M1 + 40*g2p2*M2 + 5._dp*(TrCYtIITtII)           & 
&  + 5*Conjg(L1II)*TL1II + 5*Conjg(L2II)*TL2II))/5._dp

 
 
If (TwoLoopRGE) Then 
betaBMTII2 = (888._dp*(g1p4)/25._dp + (96*g1p2*g2p2)/5._dp + 96._dp*(g2p4) - 6*CL1IIp2*L1IIp2 -    & 
&  6*CL2IIp2*L2IIp2 - (AbsL1II*(3._dp*(g1p2) + 5._dp*(g2p2) + 30._dp*(TrYdadjYd) +       & 
&  10._dp*(TrYeadjYe)))/5._dp - 2._dp*(TrYeCYtIIYtIIadjYe) - 6._dp*(TrYtIIadjYzIIYzIICYtII) -& 
&  (3*g1p2*TrYtIICYtII)/5._dp - g2p2*TrYtIICYtII - 6._dp*(TrYtIICYtIIYtIICYtII) -        & 
&  (AbsL2II*(3._dp*(g1p2) + 5._dp*(g2p2) + 30._dp*(TrYuadjYu)))/5._dp)*BMTII -           & 
&  (2*MTII*(1776*g1p4*M1 + 480*g1p2*g2p2*M1 + 480*g1p2*g2p2*M2 + 4800*g2p4*M2 +          & 
&  15*g1p2*TrCYtIITtII + 25*g2p2*TrCYtIITtII + 50._dp*(TrYeCYtIITtIIadjYe) +             & 
&  50._dp*(TrYtIIadjYeTeCYtII) + 150._dp*(TrYtIIadjYzIITzIICYtII) - 15*g1p2*M1*TrYtIICYtII -& 
&  25*g2p2*M2*TrYtIICYtII + 300._dp*(TrYtIICYtIITtIICYtII) + 150._dp*(TrYzIICYtIITtIIadjYzII) +& 
&  300*CL1IIp2*L1II*TL1II - 5*Conjg(L1II)*(L1II*(3*g1p2*M1 + 5*g2p2*M2 - 30._dp*(TradjYdTd) -& 
&  10._dp*(TradjYeTe)) - (3._dp*(g1p2) + 5._dp*(g2p2) + 30._dp*(TrYdadjYd) +             & 
&  10._dp*(TrYeadjYe))*TL1II) + 300*CL2IIp2*L2II*TL2II + 5*Conjg(L2II)*(L2II*(-          & 
& 3*g1p2*M1 - 5*g2p2*M2 + 30._dp*(TradjYuTu)) + (3._dp*(g1p2) + 5._dp*(g2p2) +           & 
&  30._dp*(TrYuadjYu))*TL2II)))/25._dp

 
DBMTII = oo16pi2*( betaBMTII1 + oo16pi2 * betaBMTII2 ) 

 
Else 
DBMTII = oo16pi2* betaBMTII1 
End If 

 
!-------------------- 
! BMZII 
!-------------------- 
 
betaBMZII1  = (2*MZII*(g1p2*M1 + 80*g3p2*M3 + 45*g2p2*M2 + 15._dp*(TradjYzIITzII))    & 
&  - (g1p2 + 45._dp*(g2p2) + 80._dp*(g3p2) - 15._dp*(TrYzIIadjYzII))*BMZII)/15._dp

 
 
If (TwoLoopRGE) Then 
betaBMZII2 = (-2*MZII*(409*g1p4*M1 + 45*g1p2*g2p2*M1 + 80*g1p2*g3p2*M1 + 80*g1p2*g3p2*M3 +         & 
&  3600*g2p2*g3p2*M3 + 16000*g3p4*M3 + 45*g1p2*g2p2*M2 + 12825*g2p4*M2 + 3600*g2p2*g3p2*M2 -& 
&  90*g1p2*TradjYzIITzII + 450._dp*(TrYdadjYdTzIIadjYzII) + 225._dp*(TrYeadjYzIITzIIadjYe) +& 
&  900._dp*(TrYsIICYsIITzIIadjYzII) + 675._dp*(TrYtIIadjYzIITzIICYtII) + 225._dp*(TrYzIIadjYeTeadjYzII) +& 
&  90*g1p2*M1*TrYzIIadjYzII + 450._dp*(TrYzIIadjYzIITdadjYd) + 900._dp*(TrYzIIadjYzIITsIICYsII) +& 
&  2250._dp*(TrYzIIadjYzIITzIIadjYzII) + 675._dp*(TrYzIICYtIITtIIadjYzII)))/225._dp +    & 
&  (409._dp*(g1p4)/450._dp + (g1p2*g2p2)/5._dp + 57._dp*(g2p4)/2._dp + (16*g1p2*g3p2)/45._dp +& 
&  16*g2p2*g3p2 + 320._dp*(g3p4)/9._dp - 2._dp*(TrYdadjYdYzIIadjYzII) - TrYeadjYzIIYzIIadjYe -& 
&  4._dp*(TrYsIICYsIIYzIIadjYzII) - 3._dp*(TrYtIIadjYzIIYzIICYtII) + (2*g1p2*TrYzIIadjYzII)/5._dp -& 
&  5._dp*(TrYzIIadjYzIIYzIIadjYzII))*BMZII

 
DBMZII = oo16pi2*( betaBMZII1 + oo16pi2 * betaBMZII2 ) 

 
Else 
DBMZII = oo16pi2* betaBMZII1 
End If 
 
!-------------------- 
! BMSII 
!-------------------- 
 
betaBMSII1  = (2*MSII*(16*g1p2*M1 + 200*g3p2*M3 + 15._dp*(TrCYsIITsII))               & 
&  + (-8*(2._dp*(g1p2) + 25._dp*(g3p2)) + 15._dp*(TrYsIICYsII))*BMSII)/15._dp

 
 
If (TwoLoopRGE) Then 
betaBMSII2 = (-4*(2*MSII*(1696*g1p4*M1 + 800*g1p2*g3p2*M1 + 800*g1p2*g3p2*M3 + 16000*g3p4*M3 +     & 
&  15*(g1p2 + 5._dp*(g3p2))*TrCYsIITsII + 225._dp*(TrYdadjYdTsIICYsII) - 15*(g1p2*M1 +   & 
&  5*g3p2*M3)*TrYsIICYsII + 225._dp*(TrYsIICYsIITdadjYd) + 900._dp*(TrYsIICYsIITsIICYsII) +& 
&  225._dp*(TrYsIICYsIITzIIadjYzII) + 225._dp*(TrYzIIadjYzIITsIICYsII)) + (-             & 
& 848._dp*(g1p4) - 800*g1p2*g3p2 - 8000._dp*(g3p4) + 225._dp*(TrYdadjYdYsIICYsII) +      & 
&  15*(g1p2 + 5._dp*(g3p2))*TrYsIICYsII + 450._dp*(TrYsIICYsIIYsIICYsII) +               & 
&  225._dp*(TrYsIICYsIIYzIIadjYzII))*BMSII))/225._dp

 
DBMSII = oo16pi2*( betaBMSII1 + oo16pi2 * betaBMSII2 ) 

 
Else 
DBMSII = oo16pi2* betaBMSII1 
End If 
 
!123  DBMTII = 0._dp
! DBMZII =  0._dp
! DBMSII = 0._dp
!-------------------- 
! mq2 
!-------------------- 
 
betamq21  = 2._dp*(adjTdTd) + 2._dp*(adjTuTu) + 2._dp*(adjYdmd2Yd) + adjYdYdmq2 +     & 
&  2._dp*(adjYumu2Yu) + adjYuYumq2 - (2*AbsM1*g1p2*id3R)/15._dp - 6*AbsM2*g2p2*id3R -    & 
&  (32*AbsM3*g3p2*id3R)/3._dp + 2*adjYdYd*mHd2 + 2*adjYuYu*mHu2 + mq2adjYdYd +           & 
&  mq2adjYuYu + g1*id3R*ooSqrt15*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamq22 = -6*AbsL1II*adjTdTd - 4._dp*(adjTdTdadjYdYd) - 8._dp*(adjTdTsIICYsIIYd) -              & 
&  4._dp*(adjTdTzIIadjYzIIYd) - 4._dp*(adjTdYdadjYdTd) - 8._dp*(adjTdYsIICYsIITd) -      & 
&  4._dp*(adjTdYzIIadjYzIITd) - 6*AbsL2II*adjTuTu - 4._dp*(adjTuTuadjYuYu) -             & 
&  4._dp*(adjTuYuadjYuTu) - 6*AbsL1II*adjYdmd2Yd - 4._dp*(adjYdmd2YdadjYdYd) -           & 
&  8._dp*(adjYdmd2YsIICYsIIYd) - 4._dp*(adjYdmd2YzIIadjYzIIYd) - 4._dp*(adjYdTdadjTdYd) -& 
&  8._dp*(adjYdTsIICTsIIYd) - 4._dp*(adjYdTzIIadjTzIIYd) - 6*AbsTL1II*adjYdYd -          & 
&  4._dp*(adjYdYdadjTdTd) - 4._dp*(adjYdYdadjYdmd2Yd) - 2._dp*(adjYdYdadjYdYdmq2) -      & 
&  3*AbsL1II*adjYdYdmq2 - 4._dp*(adjYdYdmq2adjYdYd) - 8._dp*(adjYdYsIICmd2CYsIIYd) -     & 
&  8._dp*(adjYdYsIICTsIITd) - 8._dp*(adjYdYsIICYsIImd2Yd) - 4._dp*(adjYdYsIICYsIIYdmq2) -& 
&  4._dp*(adjYdYzIIadjTzIITd) - 4._dp*(adjYdYzIIadjYzIImd2Yd) - 2._dp*(adjYdYzIIadjYzIIYdmq2) -& 
&  4._dp*(adjYdYzIIml2adjYzIIYd) - 6*AbsL2II*adjYumu2Yu - 4._dp*(adjYumu2YuadjYuYu) -    & 
&  4._dp*(adjYuTuadjTuYu) - 6*AbsTL2II*adjYuYu - 4._dp*(adjYuYuadjTuTu) - 4._dp*(adjYuYuadjYumu2Yu) -& 
&  2._dp*(adjYuYuadjYuYumq2) - 3*AbsL2II*adjYuYumq2 - 4._dp*(adjYuYumq2adjYuYu) +        & 
&  (4*adjTdTd*g1p2)/5._dp + (8*adjTuTu*g1p2)/5._dp + (4*adjYdmd2Yd*g1p2)/5._dp +         & 
&  (2*adjYdYdmq2*g1p2)/5._dp + (8*adjYumu2Yu*g1p2)/5._dp + (4*adjYuYumq2*g1p2)/5._dp +   & 
&  (2*AbsM2*g1p2*g2p2*id3R)/5._dp + 159*AbsM2*g2p4*id3R + 32*AbsM2*g2p2*g3p2*id3R -      & 
&  (4*adjTdYd*g1p2*M1)/5._dp - (8*adjTuYu*g1p2*M1)/5._dp - 18*AbsL1II*adjYdYd*mHd2 -     & 
&  8*adjYdYdadjYdYd*mHd2 - 8*adjYdYsIICYsIIYd*mHd2 - 4*adjYdYzIIadjYzIIYd*mHd2 +         & 
&  (4*adjYdYd*g1p2*mHd2)/5._dp - 18*AbsL2II*adjYuYu*mHu2 - 8*adjYuYuadjYuYu*mHu2 +       & 
&  (8*adjYuYu*g1p2*mHu2)/5._dp - 3*AbsL1II*mq2adjYdYd + (2*g1p2*mq2adjYdYd)/5._dp -      & 
&  2._dp*(mq2adjYdYdadjYdYd) - 4._dp*(mq2adjYdYsIICYsIIYd) - 2._dp*(mq2adjYdYzIIadjYzIIYd)  
betamq22 =  betamq22- 3*AbsL2II*mq2adjYuYu + (4*g1p2*mq2adjYuYu)/5._dp - 2._dp*(mq2adjYuYuadjYuYu) -        & 
&  8*adjYdYsIICYsIIYd*ms2 - 6*AbsL1II*adjYdYd*mt2 - 6*AbsL2II*adjYuYu*mtb2 -             & 
&  4*adjYdYzIIadjYzIIYd*mzz2 - 6*adjTdYd*TradjYdTd - 2*adjTdYd*TradjYeTe -               & 
&  6*adjTuYu*TradjYuTu - 6*adjYdYd*TrCTdTpTd - 6*adjYdTd*TrCTdTpYd - 2*adjYdYd*TrCTeTpTe -& 
&  2*adjYdTd*TrCTeTpYe - 6*adjYuYu*TrCTuTpTu - 6*adjYuTu*TrCTuTpYu - 6*adjYdYd*Trmd2YdadjYd -& 
&  2*adjYdYd*Trme2YeadjYe - 2*adjYdYd*Trml2adjYeYe - 6*adjYdYd*Trmq2adjYdYd -            & 
&  6*adjYuYu*Trmq2adjYuYu - 6*adjYuYu*Trmu2YuadjYu - 6*adjTdTd*TrYdadjYd -               & 
&  6*adjYdmd2Yd*TrYdadjYd - 3*adjYdYdmq2*TrYdadjYd - 12*adjYdYd*mHd2*TrYdadjYd -         & 
&  3*mq2adjYdYd*TrYdadjYd - 2*adjTdTd*TrYeadjYe - 2*adjYdmd2Yd*TrYeadjYe -               & 
&  adjYdYdmq2*TrYeadjYe - 4*adjYdYd*mHd2*TrYeadjYe - mq2adjYdYd*TrYeadjYe -              & 
&  6*adjTuTu*TrYuadjYu - 6*adjYumu2Yu*TrYuadjYu - 3*adjYuYumq2*TrYuadjYu -               & 
&  12*adjYuYu*mHu2*TrYuadjYu - 3*mq2adjYuYu*TrYuadjYu + (g1p2*(180*(-1._dp*(adjYdTd) -   & 
&  2._dp*(adjYuTu) + 2*adjYdYd*M1 + 4*adjYuYu*M1) + id3R*(1227*g1p2*M1 + 5*(16*g3p2*(2._dp*(M1) +& 
&  M3) + 9*g2p2*(2._dp*(M1) + M2))))*Conjg(M1))/225._dp + (16*g3p2*id3R*(g1p2*(M1 +      & 
&  2._dp*(M3)) + 15*(34*g3p2*M3 + 3*g2p2*(2._dp*(M3) + M2)))*Conjg(M3))/45._dp +         & 
&  (g1p2*g2p2*id3R*M1*Conjg(M2))/5._dp + 16*g2p2*g3p2*id3R*M3*Conjg(M2) - 6*adjYdTd*L1II*Conjg(TL1II) -& 
&  6*adjYuTu*L2II*Conjg(TL2II) - 6*adjTdYd*Conjg(L1II)*TL1II - 6*adjTuYu*Conjg(L2II)*TL2II +& 
&  6*g2p4*id3R*Tr2(2) + (32*g3p4*id3R*Tr2(3))/3._dp + (2*g1p2*id3R*Tr2U1(1,              & 
& 1))/15._dp + 4*g1*id3R*ooSqrt15*Tr3(1)

 
Dmq2 = oo16pi2*( betamq21 + oo16pi2 * betamq22 ) 

 
Else 
Dmq2 = oo16pi2* betamq21 
End If 
 
 
Forall(i1=1:3) Dmq2(i1,i1) =  Real(Dmq2(i1,i1),dp) 
!-------------------- 
! ml2 
!-------------------- 
 
betaml21  = 2._dp*(adjTeTe) + 6._dp*(adjTzIITzII) + 2._dp*(adjYeme2Ye) +              & 
&  adjYeYeml2 + 6._dp*(adjYzIImd2YzII) + 3._dp*(adjYzIIYzIIml2) + 6._dp*(CTtIITtII)      & 
&  + 6._dp*(CYtIICml2YtII) + 3._dp*(CYtIIYtIIml2) - (6*AbsM1*g1p2*id3R)/5._dp -          & 
&  6*AbsM2*g2p2*id3R + 2*adjYeYe*mHd2 + ml2adjYeYe + 3._dp*(ml2adjYzIIYzII)              & 
&  + 3._dp*(ml2CYtIIYtII) + 6*CYtIIYtII*mt2 + 6*adjYzIIYzII*mzz2 - g1*id3R*sqrt3ov5*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betaml22 = -6*AbsL1II*adjTeTe - 4._dp*(adjTeTeadjYeYe) - 4._dp*(adjTeYeadjYeTe) - 12._dp*(adjTzIITdadjYdYzII) -& 
&  24._dp*(adjTzIITsIICYsIIYzII) - 12._dp*(adjTzIITzIIadjYzIIYzII) - 12._dp*(adjTzIIYdadjYdTzII) -& 
&  24._dp*(adjTzIIYsIICYsIITzII) - 12._dp*(adjTzIIYzIIadjYzIITzII) - 6*AbsL1II*adjYeme2Ye -& 
&  4._dp*(adjYeme2YeadjYeYe) - 4._dp*(adjYeTeadjTeYe) - 6*AbsTL1II*adjYeYe -             & 
&  4._dp*(adjYeYeadjTeTe) - 4._dp*(adjYeYeadjYeme2Ye) - 2._dp*(adjYeYeadjYeYeml2) -      & 
&  3*AbsL1II*adjYeYeml2 - 4._dp*(adjYeYeml2adjYeYe) - 12._dp*(adjYzIImd2YdadjYdYzII) -   & 
&  24._dp*(adjYzIImd2YsIICYsIIYzII) - 12._dp*(adjYzIImd2YzIIadjYzIIYzII) -               & 
&  12._dp*(adjYzIITdadjTdYzII) - 24._dp*(adjYzIITsIICTsIIYzII) - 12._dp*(adjYzIITzIIadjTzIIYzII) -& 
&  12._dp*(adjYzIIYdadjTdTzII) - 12._dp*(adjYzIIYdadjYdmd2YzII) - 6._dp*(adjYzIIYdadjYdYzIIml2) -& 
&  12._dp*(adjYzIIYdmq2adjYdYzII) - 24._dp*(adjYzIIYsIICmd2CYsIIYzII) - 24._dp*(adjYzIIYsIICTsIITzII) -& 
&  24._dp*(adjYzIIYsIICYsIImd2YzII) - 12._dp*(adjYzIIYsIICYsIIYzIIml2) - 12._dp*(adjYzIIYzIIadjTzIITzII) -& 
&  12._dp*(adjYzIIYzIIadjYzIImd2YzII) - 6._dp*(adjYzIIYzIIadjYzIIYzIIml2) -              & 
&  12._dp*(adjYzIIYzIIml2adjYzIIYzII) - 6._dp*(CTtIITpTeCYeYtII) - 18._dp*(CTtIITpTzIICYzIIYtII) -& 
&  6._dp*(CTtIITpYeCYeTtII) - 18._dp*(CTtIITpYzIICYzIITtII) - 6*AbsL1II*CTtIITtII -      & 
&  18._dp*(CTtIITtIICYtIIYtII) - 18._dp*(CTtIIYtIICYtIITtII) - 6._dp*(CYtIICml2TpYeCYeYtII) -& 
&  18._dp*(CYtIICml2TpYzIICYzIIYtII) - 6*AbsL1II*CYtIICml2YtII - 18._dp*(CYtIICml2YtIICYtIIYtII) -& 
&  6._dp*(CYtIITpTeCTeYtII) - 18._dp*(CYtIITpTzIICTzIIYtII) - 6._dp*(CYtIITpYeCme2CYeYtII) -& 
&  6._dp*(CYtIITpYeCTeTtII) - 6._dp*(CYtIITpYeCYeCml2YtII) - 3._dp*(CYtIITpYeCYeYtIIml2) -& 
&  18._dp*(CYtIITpYzIICmd2CYzIIYtII) - 18._dp*(CYtIITpYzIICTzIITtII) - 18._dp*(CYtIITpYzIICYzIICml2YtII) -& 
&  9._dp*(CYtIITpYzIICYzIIYtIIml2) - 18._dp*(CYtIITtIICTtIIYtII) - 6*AbsTL1II*CYtIIYtII -& 
&  18._dp*(CYtIIYtIICTtIITtII) - 18._dp*(CYtIIYtIICYtIICml2YtII) - 9._dp*(CYtIIYtIICYtIIYtIIml2)  
betaml22 =  betaml22- 3*AbsL1II*CYtIIYtIIml2 - 18._dp*(CYtIIYtIIml2CYtIIYtII) + (12*adjTeTe*g1p2)/5._dp -   & 
&  (4*adjTzIITzII*g1p2)/5._dp + (12*adjYeme2Ye*g1p2)/5._dp + (6*adjYeYeml2*g1p2)/5._dp - & 
&  (4*adjYzIImd2YzII*g1p2)/5._dp - (2*adjYzIIYzIIml2*g1p2)/5._dp + (36*CTtIITtII*g1p2)/5._dp +& 
&  (36*CYtIICml2YtII*g1p2)/5._dp + (18*CYtIIYtIIml2*g1p2)/5._dp + 24*CTtIITtII*g2p2 +    & 
&  24*CYtIICml2YtII*g2p2 + 12*CYtIIYtIIml2*g2p2 + 32*adjTzIITzII*g3p2 + 32*adjYzIImd2YzII*g3p2 +& 
&  64*AbsM3*adjYzIIYzII*g3p2 + 16*adjYzIIYzIIml2*g3p2 - (12*adjTeYe*g1p2*M1)/5._dp +     & 
&  (4*adjTzIIYzII*g1p2*M1)/5._dp - (36*CTtIIYtII*g1p2*M1)/5._dp - 32*adjTzIIYzII*g3p2*M3 -& 
&  24*CTtIIYtII*g2p2*M2 - 18*AbsL1II*adjYeYe*mHd2 - 8*adjYeYeadjYeYe*mHd2 -              & 
&  12*adjYzIIYdadjYdYzII*mHd2 - 6*CYtIITpYeCYeYtII*mHd2 - 12*AbsL1II*CYtIIYtII*mHd2 +    & 
&  (12*adjYeYe*g1p2*mHd2)/5._dp - 3*AbsL1II*ml2adjYeYe + (6*g1p2*ml2adjYeYe)/5._dp -     & 
&  2._dp*(ml2adjYeYeadjYeYe) - 6._dp*(ml2adjYzIIYdadjYdYzII) - 12._dp*(ml2adjYzIIYsIICYsIIYzII) -& 
&  (2*g1p2*ml2adjYzIIYzII)/5._dp + 16*g3p2*ml2adjYzIIYzII - 6._dp*(ml2adjYzIIYzIIadjYzIIYzII) -& 
&  3._dp*(ml2CYtIITpYeCYeYtII) - 9._dp*(ml2CYtIITpYzIICYzIIYtII) - 3*AbsL1II*ml2CYtIIYtII +& 
&  (18*g1p2*ml2CYtIIYtII)/5._dp + 12*g2p2*ml2CYtIIYtII - 9._dp*(ml2CYtIIYtIICYtIIYtII) - & 
&  24*adjYzIIYsIICYsIIYzII*ms2 - 6*AbsL1II*adjYeYe*mt2 - 6*CYtIITpYeCYeYtII*mt2 -        & 
&  18*CYtIITpYzIICYzIIYtII*mt2 - 12*AbsL1II*CYtIIYtII*mt2 - 36*CYtIIYtIICYtIIYtII*mt2 +  & 
&  (36*CYtIIYtII*g1p2*mt2)/5._dp + 24*CYtIIYtII*g2p2*mt2 - 12*adjYzIIYdadjYdYzII*mzz2 -  & 
&  24*adjYzIIYsIICYsIIYzII*mzz2 - 24*adjYzIIYzIIadjYzIIYzII*mzz2 - 18*CYtIITpYzIICYzIIYtII*mzz2 -& 
&  (4*adjYzIIYzII*g1p2*mzz2)/5._dp + 32*adjYzIIYzII*g3p2*mzz2 - 6*adjTeYe*TradjYdTd -    & 
&  2*adjTeYe*TradjYeTe - 6*adjTzIIYzII*TradjYzIITzII - 6*adjYeYe*TrCTdTpTd -             & 
&  6*adjYeTe*TrCTdTpYd - 2*adjYeYe*TrCTeTpTe - 2*adjYeTe*TrCTeTpYe - 6*CYtIIYtII*TrCTtIITtII  
betaml22 =  betaml22- 6*adjYzIIYzII*TrCTzIITpTzII - 6*adjYzIITzII*TrCTzIITpYzII - 6*CTtIIYtII*TrCYtIITtII - & 
&  6*adjYeYe*Trmd2YdadjYd - 6*adjYzIIYzII*Trmd2YzIIadjYzII - 2*adjYeYe*Trme2YeadjYe -    & 
&  2*adjYeYe*Trml2adjYeYe - 6*adjYzIIYzII*Trml2adjYzIIYzII - 12*CYtIIYtII*Trml2CYtIIYtII -& 
&  6*adjYeYe*Trmq2adjYdYd - 6*adjTeTe*TrYdadjYd - 6*adjYeme2Ye*TrYdadjYd -               & 
&  3*adjYeYeml2*TrYdadjYd - 12*adjYeYe*mHd2*TrYdadjYd - 3*ml2adjYeYe*TrYdadjYd -         & 
&  2*adjTeTe*TrYeadjYe - 2*adjYeme2Ye*TrYeadjYe - adjYeYeml2*TrYeadjYe - 4*adjYeYe*mHd2*TrYeadjYe -& 
&  ml2adjYeYe*TrYeadjYe - 6*CYtIITtII*TrYtIICTtII - 6*CTtIITtII*TrYtIICYtII -            & 
&  6*CYtIICml2YtII*TrYtIICYtII - 3*CYtIIYtIIml2*TrYtIICYtII - 3*ml2CYtIIYtII*TrYtIICYtII -& 
&  12*CYtIIYtII*mt2*TrYtIICYtII - 6*adjTzIITzII*TrYzIIadjYzII - 6*adjYzIImd2YzII*TrYzIIadjYzII -& 
&  3*adjYzIIYzIIml2*TrYzIIadjYzII - 3*ml2adjYzIIYzII*TrYzIIadjYzII - 12*adjYzIIYzII*mzz2*TrYzIIadjYzII +& 
&  (g1p2*(20*(-3._dp*(adjYeTe) + adjYzIITzII - 9._dp*(CYtIITtII) + 6*adjYeYe*M1 -        & 
&  2*adjYzIIYzII*M1 + 18*CYtIIYtII*M1) + 9*id3R*(139*g1p2*M1 + 5*g2p2*(2._dp*(M1) +      & 
&  M2)))*Conjg(M1))/25._dp - 32*adjYzIITzII*g3p2*Conjg(M3) + (3*g2p2*(-40._dp*(CYtIITtII) +& 
&  80*CYtIIYtII*M2 + id3R*(265*g2p2*M2 + 3*g1p2*(M1 + 2._dp*(M2))))*Conjg(M2))/5._dp -   & 
&  6*adjYeTe*L1II*Conjg(TL1II) - 6*CYtIITtII*L1II*Conjg(TL1II) - 6*adjTeYe*Conjg(L1II)*TL1II -& 
&  6*CTtIIYtII*Conjg(L1II)*TL1II + 6*g2p4*id3R*Tr2(2) + (6*g1p2*id3R*Tr2U1(1,            & 
& 1))/5._dp - 4*g1*id3R*sqrt3ov5*Tr3(1)

 
Dml2 = oo16pi2*( betaml21 + oo16pi2 * betaml22 ) 

 
Else 
Dml2 = oo16pi2* betaml21 
End If 
 
 
Forall(i1=1:3) Dml2(i1,i1) =  Real(Dml2(i1,i1),dp) 
!-------------------- 
! mHd2 
!-------------------- 
 
betamHd21  = 6._dp*(AbsTL1II) - (6*AbsM1*g1p2)/5._dp - 6*AbsM2*g2p2 + 6*AbsL1II*(2._dp*(mHd2)& 
&  + mt2) + 6._dp*(TrCTdTpTd) + 2._dp*(TrCTeTpTe) + 6._dp*(Trmd2YdadjYd) +               & 
&  2._dp*(Trme2YeadjYe) + 2._dp*(Trml2adjYeYe) + 6._dp*(Trmq2adjYdYd) + 6*mHd2*TrYdadjYd +& 
&  2*mHd2*TrYeadjYe - g1*sqrt3ov5*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamHd22 = (36*AbsTL1II*g1p2)/5._dp + 24*AbsTL1II*g2p2 + (18*AbsM2*g1p2*g2p2)/5._dp +            & 
&  159*AbsM2*g2p4 - 48*CL1IIp2*L1IIp2*(2._dp*(mHd2) + mt2) - (4*g1p2*TrCTdTpTd)/5._dp +  & 
&  32*g3p2*TrCTdTpTd + (4*g1p2*M1*TrCTdTpYd)/5._dp - 32*g3p2*M3*TrCTdTpYd +              & 
&  (12*g1p2*TrCTeTpTe)/5._dp - (12*g1p2*M1*TrCTeTpYe)/5._dp - 12._dp*(TrCYsIITsIICTdTpYd) -& 
&  3._dp*(TrCYtIITpYeCTeTtII) - (4*g1p2*Trmd2YdadjYd)/5._dp + 32*g3p2*Trmd2YdadjYd -     & 
&  36._dp*(Trmd2YdadjYdYdadjYd) - 24._dp*(Trmd2YdadjYdYsIICYsII) - 12._dp*(Trmd2YdadjYdYzIIadjYzII) -& 
&  6._dp*(Trmd2YdadjYuYuadjYd) - 24._dp*(Trmd2YsIICYsIIYdadjYd) - 6._dp*(Trmd2YzIIadjYeYeadjYzII) -& 
&  12._dp*(Trmd2YzIIadjYzIIYdadjYd) + (12*g1p2*Trme2YeadjYe)/5._dp - 12._dp*(Trme2YeadjYeYeadjYe) -& 
&  6._dp*(Trme2YeadjYzIIYzIIadjYe) - 6._dp*(Trme2YeCYtIIYtIIadjYe) + (12*g1p2*Trml2adjYeYe)/5._dp -& 
&  12._dp*(Trml2adjYeYeadjYeYe) - 6._dp*(Trml2adjYeYeadjYzIIYzII) - 6._dp*(Trml2adjYeYeCYtIIYtII) -& 
&  12._dp*(Trml2adjYzIIYdadjYdYzII) - 6._dp*(Trml2adjYzIIYzIIadjYeYe) - 6._dp*(Trml2CYtIIYtIIadjYeYe) -& 
&  (4*g1p2*Trmq2adjYdYd)/5._dp + 32*g3p2*Trmq2adjYdYd - 36._dp*(Trmq2adjYdYdadjYdYd) -   & 
&  6._dp*(Trmq2adjYdYdadjYuYu) - 24._dp*(Trmq2adjYdYsIICYsIIYd) - 12._dp*(Trmq2adjYdYzIIadjYzIIYd) -& 
&  6._dp*(Trmq2adjYuYuadjYdYd) - 6._dp*(Trmu2YuadjYdYdadjYu) - 36._dp*(TrYdadjTdTdadjYd) -& 
&  12._dp*(TrYdadjTdTsIICYsII) - 12._dp*(TrYdadjTdTzIIadjYzII) - 6._dp*(TrYdadjTuTuadjYd) -& 
&  18*AbsTL1II*TrYdadjYd + 64*AbsM3*g3p2*TrYdadjYd - (4*g1p2*mHd2*TrYdadjYd)/5._dp +     & 
&  32*g3p2*mHd2*TrYdadjYd - 36._dp*(TrYdadjYdTdadjTd) - 24._dp*(TrYdadjYdTsIICTsII) -    & 
&  12._dp*(TrYdadjYdTzIIadjTzII) - 36*mHd2*TrYdadjYdYdadjYd - 24._dp*(TrYdadjYdYsIICmd2CYsII) -& 
&  24*mHd2*TrYdadjYdYsIICYsII - 24*ms2*TrYdadjYdYsIICYsII - 12*mHd2*TrYdadjYdYzIIadjYzII -& 
&  12*mzz2*TrYdadjYdYzIIadjYzII - 6._dp*(TrYdadjYuTuadjTd) - 6*mHd2*TrYdadjYuYuadjYd -   & 
&  6*mHu2*TrYdadjYuYuadjYd - 12._dp*(TrYeadjTeTeadjYe) - 6._dp*(TrYeadjTzIITzIIadjYe)  
betamHd22 =  betamHd22- 6*AbsTL1II*TrYeadjYe + (12*g1p2*mHd2*TrYeadjYe)/5._dp - 12._dp*(TrYeadjYeTeadjTe) -   & 
&  12*mHd2*TrYeadjYeYeadjYe - 6._dp*(TrYeadjYzIITzIIadjTe) - 6*mHd2*TrYeadjYzIIYzIIadjYe -& 
&  6*mzz2*TrYeadjYzIIYzIIadjYe - 6._dp*(TrYeCTtIITtIIadjYe) - 6._dp*(TrYeCYtIICml2YtIIadjYe) -& 
&  3._dp*(TrYeCYtIITtIIadjTe) - 6*mHd2*TrYeCYtIIYtIIadjYe - 6*mt2*TrYeCYtIIYtIIadjYe -   & 
&  6._dp*(TrYsIICTdTpTdCYsII) - 24._dp*(TrYsIICTsIITdadjYd) - 18._dp*(TrYsIICYsIITdadjTd) -& 
&  3._dp*(TrYtIIadjTeTeCYtII) - 6._dp*(TrYtIIadjYeTeCTtII) - 6*AbsTL1II*TrYtIICYtII -    & 
&  3._dp*(TrYtIICYtIITpTeCTe) - 6._dp*(TrYuadjTdTdadjYu) - 6._dp*(TrYuadjYdTdadjTu) -    & 
&  6._dp*(TrYzIIadjTeTeadjYzII) - 12._dp*(TrYzIIadjTzIITdadjYd) - 6._dp*(TrYzIIadjYeTeadjTzII) -& 
&  12._dp*(TrYzIIadjYzIITdadjTd) + (g1p2*(1251*g1p2*M1 + 90*g2p2*M1 + 45*g2p2*M2 +       & 
&  20._dp*(TradjYdTd) - 60._dp*(TradjYeTe) - 40*M1*TrYdadjYd + 120*M1*TrYeadjYe)*Conjg(M1))/25._dp -& 
&  32*g3p2*TradjYdTd*Conjg(M3) + (9*g1p2*g2p2*M1*Conjg(M2))/5._dp - (36*g1p2*L1II*M1*Conjg(TL1II))/5._dp -& 
&  24*g2p2*L1II*M2*Conjg(TL1II) - 18*L1II*TradjYdTd*Conjg(TL1II) - 6*L1II*TradjYeTe*Conjg(TL1II) -& 
&  6*L1II*TrCYtIITtII*Conjg(TL1II) + (6*Conjg(L1II)*(-80*AbsTL1II*L1II + 12*g1p2*L1II*mHd2 +& 
&  40*g2p2*L1II*mHd2 + 6*g1p2*L1II*mt2 + 20*g2p2*L1II*mt2 - 15*L1II*TrCTdTpTd -          & 
&  5*L1II*TrCTeTpTe - 5*L1II*TrCTtIITtII - 15*L1II*Trmd2YdadjYd - 5*L1II*Trme2YeadjYe -  & 
&  5*L1II*Trml2adjYeYe - 10*L1II*Trml2CYtIIYtII - 15*L1II*Trmq2adjYdYd - 45*L1II*mHd2*TrYdadjYd -& 
&  15*L1II*mt2*TrYdadjYd - 15*L1II*mHd2*TrYeadjYe - 5*L1II*mt2*TrYeadjYe -               & 
&  10*L1II*mHd2*TrYtIICYtII - 10*L1II*mt2*TrYtIICYtII + 6*g1p2*Conjg(M1)*(2*L1II*M1 -    & 
&  TL1II) + 20*g2p2*Conjg(M2)*(2*L1II*M2 - TL1II) - 15*TrCTdTpYd*TL1II - 5*TrCTeTpYe*TL1II -& 
&  5*TrYtIICTtII*TL1II))/5._dp + 6*g2p4*Tr2(2) + (6*g1p2*Tr2U1(1,1))/5._dp -             & 
&  4*g1*sqrt3ov5*Tr3(1)

 
DmHd2 = oo16pi2*( betamHd21 + oo16pi2 * betamHd22 ) 

 
Else 
DmHd2 = oo16pi2* betamHd21 
End If 
 
 
!-------------------- 
! mHu2 
!-------------------- 
 
betamHu21  = 6._dp*(AbsTL2II) - (6*AbsM1*g1p2)/5._dp - 6*AbsM2*g2p2 + 6*AbsL2II*(2._dp*(mHu2)& 
&  + mtb2) + 6._dp*(TrCTuTpTu) + 6._dp*(Trmq2adjYuYu) + 6._dp*(Trmu2YuadjYu)             & 
&  + 6*mHu2*TrYuadjYu + g1*sqrt3ov5*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamHu22 = (36*AbsTL2II*g1p2)/5._dp + 24*AbsTL2II*g2p2 + (18*AbsM2*g1p2*g2p2)/5._dp +            & 
&  159*AbsM2*g2p4 - 48*CL2IIp2*L2IIp2*(2._dp*(mHu2) + mtb2) + (8*g1p2*TrCTuTpTu)/5._dp + & 
&  32*g3p2*TrCTuTpTu - (8*g1p2*M1*TrCTuTpYu)/5._dp - 32*g3p2*M3*TrCTuTpYu -              & 
&  6._dp*(Trmd2YdadjYuYuadjYd) - 6._dp*(Trmq2adjYdYdadjYuYu) + (8*g1p2*Trmq2adjYuYu)/5._dp +& 
&  32*g3p2*Trmq2adjYuYu - 6._dp*(Trmq2adjYuYuadjYdYd) - 36._dp*(Trmq2adjYuYuadjYuYu) -   & 
&  6._dp*(Trmu2YuadjYdYdadjYu) + (8*g1p2*Trmu2YuadjYu)/5._dp + 32*g3p2*Trmu2YuadjYu -    & 
&  36._dp*(Trmu2YuadjYuYuadjYu) - 6._dp*(TrYdadjTuTuadjYd) - 6._dp*(TrYdadjYuTuadjTd) -  & 
&  6*mHd2*TrYdadjYuYuadjYd - 6*mHu2*TrYdadjYuYuadjYd - 6._dp*(TrYuadjTdTdadjYu) -        & 
&  36._dp*(TrYuadjTuTuadjYu) - 6._dp*(TrYuadjYdTdadjTu) - 18*AbsTL2II*TrYuadjYu +        & 
&  64*AbsM3*g3p2*TrYuadjYu + (8*g1p2*mHu2*TrYuadjYu)/5._dp + 32*g3p2*mHu2*TrYuadjYu -    & 
&  36._dp*(TrYuadjYuTuadjTu) - 36*mHu2*TrYuadjYuYuadjYu + (g1p2*(1251*g1p2*M1 +          & 
&  90*g2p2*M1 + 45*g2p2*M2 - 40._dp*(TradjYuTu) + 80*M1*TrYuadjYu)*Conjg(M1))/25._dp -   & 
&  32*g3p2*TradjYuTu*Conjg(M3) + (9*g1p2*g2p2*M1*Conjg(M2))/5._dp - (36*g1p2*L2II*M1*Conjg(TL2II))/5._dp -& 
&  24*g2p2*L2II*M2*Conjg(TL2II) - 18*L2II*TradjYuTu*Conjg(TL2II) + (6*Conjg(L2II)*(-     & 
& 80*AbsTL2II*L2II + 12*g1p2*L2II*mHu2 + 40*g2p2*L2II*mHu2 + 6*g1p2*L2II*mtb2 +          & 
&  20*g2p2*L2II*mtb2 - 15*L2II*TrCTuTpTu - 15*L2II*Trmq2adjYuYu - 15*L2II*Trmu2YuadjYu - & 
&  45*L2II*mHu2*TrYuadjYu - 15*L2II*mtb2*TrYuadjYu + 6*g1p2*Conjg(M1)*(2*L2II*M1 -       & 
&  TL2II) + 20*g2p2*Conjg(M2)*(2*L2II*M2 - TL2II) - 15*TrCTuTpYu*TL2II))/5._dp +         & 
&  6*g2p4*Tr2(2) + (6*g1p2*Tr2U1(1,1))/5._dp + 4*g1*sqrt3ov5*Tr3(1)

 
DmHu2 = oo16pi2*( betamHu21 + oo16pi2 * betamHu22 ) 

 
Else 
DmHu2 = oo16pi2* betamHu21 
End If 
 
 
!-------------------- 
! md2 
!-------------------- 
 
betamd21  = (-8*AbsM1*g1p2*id3R)/15._dp - (32*AbsM3*g3p2*id3R)/3._dp + 2._dp*(md2YdadjYd)& 
&  + 4._dp*(md2YsIICYsII) + 2._dp*(md2YzIIadjYzII) + 4._dp*(TdadjTd) + 8._dp*(TsIICTsII) & 
&  + 4._dp*(TzIIadjTzII) + 4*mHd2*YdadjYd + 2._dp*(YdadjYdmd2) + 4._dp*(Ydmq2adjYd)      & 
&  + 8._dp*(YsIICmd2CYsII) + 8*ms2*YsIICYsII + 4._dp*(YsIICYsIImd2) + 4*mzz2*YzIIadjYzII +& 
&  2._dp*(YzIIadjYzIImd2) + 4._dp*(YzIIml2adjYzII) + 2*g1*id3R*ooSqrt15*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamd22 = -6*AbsL1II*md2YdadjYd + (2*g1p2*md2YdadjYd)/5._dp + 6*g2p2*md2YdadjYd -               & 
&  2._dp*(md2YdadjYdYdadjYd) - 2._dp*(md2YdadjYuYuadjYd) - 8._dp*(md2YsIICYdTpYdCYsII) + & 
&  (32*g1p2*md2YsIICYsII)/15._dp + (80*g3p2*md2YsIICYsII)/3._dp - 16._dp*(md2YsIICYsIIYsIICYsII) -& 
&  8._dp*(md2YsIICYzIITpYzIICYsII) - 2._dp*(md2YzIIadjYeYeadjYzII) + (2*g1p2*md2YzIIadjYzII)/5._dp +& 
&  6*g2p2*md2YzIIadjYzII - 6._dp*(md2YzIIadjYzIIYzIIadjYzII) - 6._dp*(md2YzIICYtIIYtIIadjYzII) -& 
&  12*AbsL1II*TdadjTd + (4*g1p2*TdadjTd)/5._dp + 12*g2p2*TdadjTd - 4._dp*(TdadjTdYdadjYd) -& 
&  4._dp*(TdadjTuYuadjYd) - 4._dp*(TdadjYdYdadjTd) - 4._dp*(TdadjYdYsIICTsII) -          & 
&  4._dp*(TdadjYuYuadjTd) - 12*TdadjYd*TrCTdTpYd - 4*TdadjYd*TrCTeTpYe - 6*md2YdadjYd*TrYdadjYd -& 
&  12*TdadjTd*TrYdadjYd - 2*md2YdadjYd*TrYeadjYe - 4*TdadjTd*TrYeadjYe - 4*md2YsIICYsII*TrYsIICYsII -& 
&  2*md2YzIIadjYzII*TrYzIIadjYzII - 12._dp*(TsIICTdTpYdCYsII) + (64*g1p2*TsIICTsII)/15._dp +& 
&  (160*g3p2*TsIICTsII)/3._dp - 8*TrYsIICYsII*TsIICTsII - 32._dp*(TsIICTsIIYsIICYsII) -  & 
&  12._dp*(TsIICTzIITpYzIICYsII) - 16._dp*(TsIICYdTpYdCTsII) - 8*TrYsIICTsII*TsIICYsII - & 
&  32._dp*(TsIICYsIIYsIICTsII) - 16._dp*(TsIICYzIITpYzIICTsII) - 4._dp*(TsIIYdadjTdCYsII) -& 
&  4._dp*(TsIIYzIIadjTzIICYsII) - 4._dp*(TzIIadjTeYeadjYzII) + (4*g1p2*TzIIadjTzII)/5._dp +& 
&  12*g2p2*TzIIadjTzII - 4*TrYzIIadjYzII*TzIIadjTzII - 12._dp*(TzIIadjTzIIYzIIadjYzII) - & 
&  4._dp*(TzIIadjYeYeadjTzII) - 4*TrCTzIITpYzII*TzIIadjYzII - 4._dp*(TzIIadjYzIIYsIICTsII) -& 
&  12._dp*(TzIIadjYzIIYzIIadjTzII) - 12._dp*(TzIICTtIIYtIIadjYzII) - 12._dp*(TzIICYtIIYtIIadjTzII) -& 
&  (4*g1p2*M1*YdadjTd)/5._dp - 12*g2p2*M2*YdadjTd - 12*TradjYdTd*YdadjTd -               & 
&  4*TradjYeTe*YdadjTd - 4._dp*(YdadjTdTdadjYd) - 4._dp*(YdadjTuTuadjYd) -               & 
&  12*AbsTL1II*YdadjYd + 24*AbsM2*g2p2*YdadjYd - 36*AbsL1II*mHd2*YdadjYd +               & 
&  (4*g1p2*mHd2*YdadjYd)/5._dp + 12*g2p2*mHd2*YdadjYd - 12*AbsL1II*mt2*YdadjYd  
betamd22 =  betamd22- 12*TrCTdTpTd*YdadjYd - 4*TrCTeTpTe*YdadjYd - 12*Trmd2YdadjYd*YdadjYd - 4*Trme2YeadjYe*YdadjYd -& 
&  4*Trml2adjYeYe*YdadjYd - 12*Trmq2adjYdYd*YdadjYd - 24*mHd2*TrYdadjYd*YdadjYd -        & 
&  8*mHd2*TrYeadjYe*YdadjYd - 6*AbsL1II*YdadjYdmd2 + (2*g1p2*YdadjYdmd2)/5._dp +         & 
&  6*g2p2*YdadjYdmd2 - 6*TrYdadjYd*YdadjYdmd2 - 2*TrYeadjYe*YdadjYdmd2 - 4._dp*(YdadjYdmd2YdadjYd) -& 
&  4._dp*(YdadjYdTdadjTd) - 8*mHd2*YdadjYdYdadjYd - 2._dp*(YdadjYdYdadjYdmd2) -          & 
&  4._dp*(YdadjYdYdmq2adjYd) - 4._dp*(YdadjYumu2YuadjYd) - 4._dp*(YdadjYuTuadjTd) -      & 
&  4*mHd2*YdadjYuYuadjYd - 4*mHu2*YdadjYuYuadjYd - 2._dp*(YdadjYuYuadjYdmd2) -           & 
&  4._dp*(YdadjYuYumq2adjYd) - 12*AbsL1II*Ydmq2adjYd + (4*g1p2*Ydmq2adjYd)/5._dp +       & 
&  12*g2p2*Ydmq2adjYd - 12*TrYdadjYd*Ydmq2adjYd - 4*TrYeadjYe*Ydmq2adjYd -               & 
&  4._dp*(Ydmq2adjYdYdadjYd) - 4._dp*(Ydmq2adjYuYuadjYd) - 16._dp*(YsIICmd2CYdTpYdCYsII) +& 
&  (64*g1p2*YsIICmd2CYsII)/15._dp + (160*g3p2*YsIICmd2CYsII)/3._dp - 8*TrYsIICYsII*YsIICmd2CYsII -& 
&  32._dp*(YsIICmd2CYsIIYsIICYsII) - 16._dp*(YsIICmd2CYzIITpYzIICYsII) - 16._dp*(YsIICTdTpTdCYsII) -& 
&  (64*g1p2*M1*YsIICTsII)/15._dp - (160*g3p2*M3*YsIICTsII)/3._dp - 8*TrCYsIITsII*YsIICTsII -& 
&  32._dp*(YsIICTsIITsIICYsII) - 16._dp*(YsIICTzIITpTzIICYsII) - 16._dp*(YsIICYdCmq2TpYdCYsII) -& 
&  8._dp*(YsIICYdTpTdCTsII) - 16._dp*(YsIICYdTpYdCmd2CYsII) - 16*mHd2*YsIICYdTpYdCYsII - & 
&  16*ms2*YsIICYdTpYdCYsII - 8._dp*(YsIICYdTpYdCYsIImd2) + (64*g1p2*ms2*YsIICYsII)/15._dp +& 
&  (160*g3p2*ms2*YsIICYsII)/3._dp - 8*TrCTsIITsII*YsIICYsII - 16*Trmd2YsIICYsII*YsIICYsII -& 
&  16*ms2*TrYsIICYsII*YsIICYsII + (32*g1p2*YsIICYsIImd2)/15._dp + (80*g3p2*YsIICYsIImd2)/3._dp -& 
&  4*TrYsIICYsII*YsIICYsIImd2 - 32._dp*(YsIICYsIImd2YsIICYsII) - 32._dp*(YsIICYsIITsIICTsII) -& 
&  32._dp*(YsIICYsIIYsIICmd2CYsII) - 64*ms2*YsIICYsIIYsIICYsII - 16._dp*(YsIICYsIIYsIICYsIImd2) -& 
&  16._dp*(YsIICYzIICml2TpYzIICYsII) - 8._dp*(YsIICYzIITpTzIICTsII) - 16._dp*(YsIICYzIITpYzIICmd2CYsII)  
betamd22 =  betamd22- 16*ms2*YsIICYzIITpYzIICYsII - 16*mzz2*YsIICYzIITpYzIICYsII - 8._dp*(YsIICYzIITpYzIICYsIImd2) -& 
&  4._dp*(YsIITdadjYdCTsII) - 4._dp*(YsIITzIIadjYzIICTsII) - 4._dp*(YzIIadjTeTeadjYzII) -& 
&  (4*g1p2*M1*YzIIadjTzII)/5._dp - 12*g2p2*M2*YzIIadjTzII - 4*TradjYzIITzII*YzIIadjTzII -& 
&  12._dp*(YzIIadjTzIITzIIadjYzII) - 4._dp*(YzIIadjYeme2YeadjYzII) - 4._dp*(YzIIadjYeTeadjTzII) -& 
&  4*mHd2*YzIIadjYeYeadjYzII - 4*mzz2*YzIIadjYeYeadjYzII - 2._dp*(YzIIadjYeYeadjYzIImd2) -& 
&  4._dp*(YzIIadjYeYeml2adjYzII) + 24*AbsM2*g2p2*YzIIadjYzII + (4*g1p2*mzz2*YzIIadjYzII)/5._dp +& 
&  12*g2p2*mzz2*YzIIadjYzII - 4*TrCTzIITpTzII*YzIIadjYzII - 4*Trmd2YzIIadjYzII*YzIIadjYzII -& 
&  4*Trml2adjYzIIYzII*YzIIadjYzII - 8*mzz2*TrYzIIadjYzII*YzIIadjYzII + (2*g1p2*YzIIadjYzIImd2)/5._dp +& 
&  6*g2p2*YzIIadjYzIImd2 - 2*TrYzIIadjYzII*YzIIadjYzIImd2 - 12._dp*(YzIIadjYzIImd2YzIIadjYzII) -& 
&  12._dp*(YzIIadjYzIITzIIadjTzII) - 24*mzz2*YzIIadjYzIIYzIIadjYzII - 6._dp*(YzIIadjYzIIYzIIadjYzIImd2) -& 
&  12._dp*(YzIIadjYzIIYzIIml2adjYzII) - 12._dp*(YzIICTtIITtIIadjYzII) - 12._dp*(YzIICYtIICml2YtIIadjYzII) -& 
&  12._dp*(YzIICYtIITtIIadjTzII) - 12*mt2*YzIICYtIIYtIIadjYzII - 12*mzz2*YzIICYtIIYtIIadjYzII -& 
&  6._dp*(YzIICYtIIYtIIadjYzIImd2) - 12._dp*(YzIICYtIIYtIIml2adjYzII) - 4._dp*(YzIIml2adjYeYeadjYzII) +& 
&  (4*g1p2*YzIIml2adjYzII)/5._dp + 12*g2p2*YzIIml2adjYzII - 4*TrYzIIadjYzII*YzIIml2adjYzII -& 
&  12._dp*(YzIIml2adjYzIIYzIIadjYzII) - 12._dp*(YzIIml2CYtIIYtIIadjYzII) +               & 
&  (4*g1p2*(4*id3R*(309*g1p2*M1 + 20*g3p2*(2._dp*(M1) + M3)) + 15*(-3._dp*(TdadjYd) -    & 
&  16._dp*(TsIICYsII) - 3._dp*(TzIIadjYzII) + 6*M1*YdadjYd + 32*M1*YsIICYsII +           & 
&  6*M1*YzIIadjYzII))*Conjg(M1))/225._dp + (32*g3p2*(id3R*(255*g3p2*M3 + 2*g1p2*(M1 +    & 
&  2._dp*(M3))) + 75*(-1._dp*(TsIICYsII) + 2*M3*YsIICYsII))*Conjg(M3))/45._dp -          & 
&  12*g2p2*TdadjYd*Conjg(M2) - 12*g2p2*TzIIadjYzII*Conjg(M2) - 12*L1II*TdadjYd*Conjg(TL1II) -& 
&  12*YdadjTd*Conjg(L1II)*TL1II + (32*g3p4*id3R*Tr2(3))/3._dp + (8*g1p2*id3R*Tr2U1(1,    & 
& 1))/15._dp + 8*g1*id3R*ooSqrt15*Tr3(1)

 
Dmd2 = oo16pi2*( betamd21 + oo16pi2 * betamd22 ) 

 
Else 
Dmd2 = oo16pi2* betamd21 
End If 
 
 
Forall(i1=1:3) Dmd2(i1,i1) =  Real(Dmd2(i1,i1),dp) 
!-------------------- 
! mu2 
!-------------------- 
 
betamu21  = (-32*AbsM1*g1p2*id3R)/15._dp - (32*AbsM3*g3p2*id3R)/3._dp +               & 
&  2._dp*(mu2YuadjYu) + 4._dp*(TuadjTu) + 4*mHu2*YuadjYu + 2._dp*(YuadjYumu2)            & 
&  + 4._dp*(Yumq2adjYu) - 4*g1*id3R*ooSqrt15*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamu22 = -2._dp*(mu2YuadjYdYdadjYu) - 6*AbsL2II*mu2YuadjYu - (2*g1p2*mu2YuadjYu)/5._dp +       & 
&  6*g2p2*mu2YuadjYu - 2._dp*(mu2YuadjYuYuadjYu) - 6*mu2YuadjYu*TrYuadjYu -              & 
&  4._dp*(TuadjTdYdadjYu) - 12*AbsL2II*TuadjTu - (4*g1p2*TuadjTu)/5._dp + 12*g2p2*TuadjTu -& 
&  12*TrYuadjYu*TuadjTu - 4._dp*(TuadjTuYuadjYu) - 4._dp*(TuadjYdYdadjTu) -              & 
&  12*TrCTuTpYu*TuadjYu - 4._dp*(TuadjYuYuadjTu) - 4._dp*(YuadjTdTdadjYu) +              & 
&  (4*g1p2*M1*YuadjTu)/5._dp - 12*g2p2*M2*YuadjTu - 12*TradjYuTu*YuadjTu -               & 
&  4._dp*(YuadjTuTuadjYu) - 4._dp*(YuadjYdmd2YdadjYu) - 4._dp*(YuadjYdTdadjTu) -         & 
&  4*mHd2*YuadjYdYdadjYu - 4*mHu2*YuadjYdYdadjYu - 2._dp*(YuadjYdYdadjYumu2) -           & 
&  4._dp*(YuadjYdYdmq2adjYu) - 12*AbsTL2II*YuadjYu + 24*AbsM2*g2p2*YuadjYu -             & 
&  36*AbsL2II*mHu2*YuadjYu - (4*g1p2*mHu2*YuadjYu)/5._dp + 12*g2p2*mHu2*YuadjYu -        & 
&  12*AbsL2II*mtb2*YuadjYu - 12*TrCTuTpTu*YuadjYu - 12*Trmq2adjYuYu*YuadjYu -            & 
&  12*Trmu2YuadjYu*YuadjYu - 24*mHu2*TrYuadjYu*YuadjYu - 6*AbsL2II*YuadjYumu2 -          & 
&  (2*g1p2*YuadjYumu2)/5._dp + 6*g2p2*YuadjYumu2 - 6*TrYuadjYu*YuadjYumu2 -              & 
&  4._dp*(YuadjYumu2YuadjYu) - 4._dp*(YuadjYuTuadjTu) - 8*mHu2*YuadjYuYuadjYu -          & 
&  2._dp*(YuadjYuYuadjYumu2) - 4._dp*(YuadjYuYumq2adjYu) - 4._dp*(Yumq2adjYdYdadjYu) -   & 
&  12*AbsL2II*Yumq2adjYu - (4*g1p2*Yumq2adjYu)/5._dp + 12*g2p2*Yumq2adjYu -              & 
&  12*TrYuadjYu*Yumq2adjYu - 4._dp*(Yumq2adjYuYuadjYu) + (4*g1p2*(32*id3R*(159*g1p2*M1 + & 
&  10*g3p2*(2._dp*(M1) + M3)) + 45*(TuadjYu - 2*M1*YuadjYu))*Conjg(M1))/225._dp +        & 
&  (32*g3p2*id3R*(255*g3p2*M3 + 8*g1p2*(M1 + 2._dp*(M3)))*Conjg(M3))/45._dp -            & 
&  12*g2p2*TuadjYu*Conjg(M2) - 12*L2II*TuadjYu*Conjg(TL2II) - 12*YuadjTu*Conjg(L2II)*TL2II +& 
&  (32*g3p4*id3R*Tr2(3))/3._dp + (32*g1p2*id3R*Tr2U1(1,1))/15._dp - 16*g1*id3R*ooSqrt15*Tr3(1)

 
Dmu2 = oo16pi2*( betamu21 + oo16pi2 * betamu22 ) 

 
Else 
Dmu2 = oo16pi2* betamu21 
End If 
 
 
Forall(i1=1:3) Dmu2(i1,i1) =  Real(Dmu2(i1,i1),dp) 
!-------------------- 
! me2 
!-------------------- 
 
betame21  = (-24*AbsM1*g1p2*id3R)/5._dp + 2*(me2YeadjYe + 2._dp*(TeadjTe)             & 
&  + 2*mHd2*YeadjYe + YeadjYeme2 + 2._dp*(Yeml2adjYe)) + 2*g1*id3R*sqrt3ov5*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betame22 = (2*(6*g1p2*(444*g1p2*id3R*M1 + 5*(TeadjYe - 2*M1*YeadjYe))*Conjg(M1) - 5*(15*AbsL1II*me2YeadjYe +& 
&  3*g1p2*me2YeadjYe - 15*g2p2*me2YeadjYe + 5._dp*(me2YeadjYeYeadjYe) + 15._dp*(me2YeadjYzIIYzIIadjYe) +& 
&  15._dp*(me2YeCYtIIYtIIadjYe) + 30*AbsL1II*TeadjTe + 6*g1p2*TeadjTe - 30*g2p2*TeadjTe +& 
&  10._dp*(TeadjTeYeadjYe) + 30._dp*(TeadjTzIIYzIIadjYe) + 10._dp*(TeadjYeYeadjTe) +     & 
&  30._dp*(TeadjYzIIYzIIadjTe) + 30._dp*(TeCTtIIYtIIadjYe) + 30._dp*(TeCYtIIYtIIadjTe) + & 
&  30*TeadjYe*TrCTdTpYd + 10*TeadjYe*TrCTeTpYe + 15*me2YeadjYe*TrYdadjYd +               & 
&  30*TeadjTe*TrYdadjYd + 5*me2YeadjYe*TrYeadjYe + 10*TeadjTe*TrYeadjYe + 10._dp*(YeadjTeTeadjYe) +& 
&  30._dp*(YeadjTzIITzIIadjYe) + 2*(15._dp*(AbsTL1II) - 30*AbsM2*g2p2 + 3*g1p2*mHd2 -    & 
&  15*g2p2*mHd2 + 15*AbsL1II*(3._dp*(mHd2) + mt2) + 15._dp*(TrCTdTpTd) + 5._dp*(TrCTeTpTe) +& 
&  15._dp*(Trmd2YdadjYd) + 5._dp*(Trme2YeadjYe) + 5._dp*(Trml2adjYeYe) + 15._dp*(Trmq2adjYdYd) +& 
&  30*mHd2*TrYdadjYd + 10*mHd2*TrYeadjYe)*YeadjYe + 15*AbsL1II*YeadjYeme2 +              & 
&  3*g1p2*YeadjYeme2 - 15*g2p2*YeadjYeme2 + 15*TrYdadjYd*YeadjYeme2 + 5*TrYeadjYe*YeadjYeme2 +& 
&  10._dp*(YeadjYeme2YeadjYe) + 10._dp*(YeadjYeTeadjTe) + 20*mHd2*YeadjYeYeadjYe +       & 
&  5._dp*(YeadjYeYeadjYeme2) + 10._dp*(YeadjYeYeml2adjYe) + 30._dp*(YeadjYzIImd2YzIIadjYe) +& 
&  30._dp*(YeadjYzIITzIIadjTe) + 30*mHd2*YeadjYzIIYzIIadjYe + 30*mzz2*YeadjYzIIYzIIadjYe +& 
&  15._dp*(YeadjYzIIYzIIadjYeme2) + 30._dp*(YeadjYzIIYzIIml2adjYe) + 30._dp*(YeCTtIITtIIadjYe) +& 
&  30._dp*(YeCYtIICml2YtIIadjYe) + 30._dp*(YeCYtIITtIIadjTe) + 30*mHd2*YeCYtIIYtIIadjYe +& 
&  30*mt2*YeCYtIIYtIIadjYe + 15._dp*(YeCYtIIYtIIadjYeme2) + 30._dp*(YeCYtIIYtIIml2adjYe) +& 
&  30*AbsL1II*Yeml2adjYe + 6*g1p2*Yeml2adjYe - 30*g2p2*Yeml2adjYe + 30*TrYdadjYd*Yeml2adjYe +& 
&  10*TrYeadjYe*Yeml2adjYe + 10._dp*(Yeml2adjYeYeadjYe) + 30._dp*(Yeml2adjYzIIYzIIadjYe) +& 
&  30._dp*(Yeml2CYtIIYtIIadjYe) + 30*g2p2*TeadjYe*Conjg(M2) + 30*L1II*TeadjYe*Conjg(TL1II) +& 
&  YeadjTe*(-6*g1p2*M1 + 30*g2p2*M2 + 30._dp*(TradjYdTd) + 10._dp*(TradjYeTe) +          & 
&  30*Conjg(L1II)*TL1II)) + 20*g1*id3R*(3*g1*Tr2U1(1,1) + sqrt15*Tr3(1))))/25._dp

 
Dme2 = oo16pi2*( betame21 + oo16pi2 * betame22 ) 

 
Else 
Dme2 = oo16pi2* betame21 
End If 
 
 
Forall(i1=1:3) Dme2(i1,i1) =  Real(Dme2(i1,i1),dp) 
!-------------------- 
! mt2 
!-------------------- 
 
betamt21  = 2._dp*(AbsTL1II) - (24*AbsM1*g1p2)/5._dp - 16*AbsM2*g2p2 + 2*AbsL1II*(2._dp*(mHd2)& 
&  + mt2) + 2._dp*(TrCTtIITtII) + 4._dp*(Trml2CYtIIYtII) + 2*mt2*TrYtIICYtII +           & 
&  2*g1*sqrt3ov5*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamt22 = (-6*AbsTL1II*g1p2)/5._dp - 2*AbsTL1II*g2p2 + (192*AbsM2*g1p2*g2p2)/5._dp +            & 
&  544*AbsM2*g2p4 - 24*CL1IIp2*L1IIp2*(2._dp*(mHd2) + mt2) - (6*g1p2*TrCTtIITtII)/5._dp -& 
&  2*g2p2*TrCTtIITtII - 2._dp*(TrCYtIITpYeCTeTtII) - 6._dp*(TrCYtIITpYzIICTzIITtII) -    & 
&  12._dp*(Trmd2YzIICYtIIYtIIadjYzII) - 4._dp*(Trme2YeCYtIIYtIIadjYe) - 4._dp*(Trml2adjYeYeCYtIIYtII) -& 
&  12._dp*(Trml2adjYzIIYzIICYtIIYtII) - (12*g1p2*Trml2CYtIIYtII)/5._dp - 4*g2p2*Trml2CYtIIYtII -& 
&  4._dp*(Trml2CYtIIYtIIadjYeYe) - 12._dp*(Trml2CYtIIYtIIadjYzIIYzII) - 48._dp*(Trml2CYtIIYtIICYtIIYtII) -& 
&  12*AbsTL1II*TrYdadjYd - 4*AbsTL1II*TrYeadjYe - 4._dp*(TrYeCTtIITtIIadjYe) -           & 
&  4._dp*(TrYeCYtIICml2YtIIadjYe) - 2._dp*(TrYeCYtIITtIIadjTe) - 4*mHd2*TrYeCYtIIYtIIadjYe -& 
&  4*mt2*TrYeCYtIIYtIIadjYe - 2._dp*(TrYtIIadjTeTeCYtII) - 6._dp*(TrYtIIadjTzIITzIICYtII) -& 
&  4._dp*(TrYtIIadjYeTeCTtII) - 12._dp*(TrYtIIadjYzIITzIICTtII) - 12*mt2*TrYtIIadjYzIIYzIICYtII -& 
&  12*mzz2*TrYtIIadjYzIIYzIICYtII - 12._dp*(TrYtIIadjYzIIYzIICYtIICml2) + (6*g1p2*M1*TrYtIICTtII)/5._dp +& 
&  2*g2p2*M2*TrYtIICTtII - 21._dp*(TrYtIICTtIITtIICYtII) - 4*AbsM2*g2p2*TrYtIICYtII -    & 
&  (6*g1p2*mt2*TrYtIICYtII)/5._dp - 2*g2p2*mt2*TrYtIICYtII - 2._dp*(TrYtIICYtIITpTeCTe) -& 
&  6._dp*(TrYtIICYtIITpTzIICTzII) - 27._dp*(TrYtIICYtIITtIICTtII) - 24*mt2*TrYtIICYtIIYtIICYtII -& 
&  12._dp*(TrYzIICTtIITtIIadjYzII) - 6._dp*(TrYzIICYtIITtIIadjTzII) + (6*g1p2*(888*g1p2*M1 +& 
&  160*g2p2*M1 + 80*g2p2*M2 + 5._dp*(TrCYtIITtII) - 10*M1*TrYtIICYtII)*Conjg(M1))/25._dp +& 
&  (96*g1p2*g2p2*M1*Conjg(M2))/5._dp + 2*g2p2*TrCYtIITtII*Conjg(M2) + (6*g1p2*L1II*M1*Conjg(TL1II))/5._dp +& 
&  2*g2p2*L1II*M2*Conjg(TL1II) - 12*L1II*TradjYdTd*Conjg(TL1II) - 4*L1II*TradjYeTe*Conjg(TL1II) -& 
&  (2*Conjg(L1II)*(120*AbsTL1II*L1II + 6*g1p2*L1II*mHd2 + 10*g2p2*L1II*mHd2 +            & 
&  3*g1p2*L1II*mt2 + 5*g2p2*L1II*mt2 + 30*L1II*TrCTdTpTd + 10*L1II*TrCTeTpTe +           & 
&  30*L1II*Trmd2YdadjYd + 10*L1II*Trme2YeadjYe + 10*L1II*Trml2adjYeYe + 30*L1II*Trmq2adjYdYd +& 
&  90*L1II*mHd2*TrYdadjYd + 30*L1II*mt2*TrYdadjYd + 30*L1II*mHd2*TrYeadjYe +             & 
&  10*L1II*mt2*TrYeadjYe + 3*g1p2*Conjg(M1)*(2*L1II*M1 - TL1II) + 5*g2p2*Conjg(M2)*(2*L1II*M2 -& 
&  TL1II) + 30*TrCTdTpYd*TL1II + 10*TrCTeTpYe*TL1II))/5._dp + 16*g2p4*Tr2(2)  
betamt22 =  betamt22+ (24*g1p2*Tr2U1(1,1))/5._dp + 8*g1*sqrt3ov5*Tr3(1)

 
Dmt2 = oo16pi2*( betamt21 + oo16pi2 * betamt22 ) 

 
Else 
Dmt2 = oo16pi2* betamt21 
End If 
 
 
!-------------------- 
! mtb2 
!-------------------- 
 
betamtb21  = 2*AbsL2II*(2._dp*(mHu2) + mtb2) - (2*(-5._dp*(AbsTL2II) + 12*AbsM1*g1p2 +& 
&  40*AbsM2*g2p2 + g1*sqrt15*Tr1(1)))/5._dp

 
 
If (TwoLoopRGE) Then 
betamtb22 = (-2*(300*CL2IIp2*L2IIp2*(2._dp*(mHu2) + mtb2) - 24*g1p2*(111*g1p2*M1 + 10*g2p2*(2._dp*(M1) +& 
&  M2))*Conjg(M1) + 5*Conjg(L2II)*(120*AbsTL2II*L2II + 6*g1p2*L2II*mHu2 + 10*g2p2*L2II*mHu2 +& 
&  3*g1p2*L2II*mtb2 + 5*g2p2*L2II*mtb2 + 30*L2II*TrCTuTpTu + 30*L2II*Trmq2adjYuYu +      & 
&  30*L2II*Trmu2YuadjYu + 90*L2II*mHu2*TrYuadjYu + 30*L2II*mtb2*TrYuadjYu +              & 
&  3*g1p2*Conjg(M1)*(2*L2II*M1 - TL2II) + 5*g2p2*Conjg(M2)*(2*L2II*M2 - TL2II) +         & 
&  30*TrCTuTpYu*TL2II) - 5*(16*g2p2*(85*g2p2*M2 + 3*g1p2*(M1 + 2._dp*(M2)))*Conjg(M2) +  & 
&  Conjg(TL2II)*(L2II*(3*g1p2*M1 + 5*g2p2*M2 - 30._dp*(TradjYuTu)) - (3._dp*(g1p2) +     & 
&  5._dp*(g2p2) + 30._dp*(TrYuadjYu))*TL2II) + 40*g2p4*Tr2(2) + 12*g1p2*Tr2U1(1,         & 
& 1) - 4*g1*sqrt15*Tr3(1))))/25._dp

 
Dmtb2 = oo16pi2*( betamtb21 + oo16pi2 * betamtb22 ) 

 
Else 
Dmtb2 = oo16pi2* betamtb21 
End If 
 
 
!-------------------- 
! ms2 
!-------------------- 
 
betams21  = (-32*AbsM1*g1p2)/15._dp - (80*AbsM3*g3p2)/3._dp + 2._dp*(TrCTsIITsII)     & 
&  + 4._dp*(Trmd2YsIICYsII) + 2*ms2*TrYsIICYsII - 4*g1*ooSqrt15*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betams22 = (-8*g1p2*TrCTsIITsII)/15._dp - (8*g3p2*TrCTsIITsII)/3._dp - 4._dp*(TrCYsIITsIICTdTpYd) -& 
&  4._dp*(TrCYsIITsIICTzIITpYzII) - 8._dp*(Trmd2YdadjYdYsIICYsII) - (16*g1p2*Trmd2YsIICYsII)/15._dp -& 
&  (16*g3p2*Trmd2YsIICYsII)/3._dp - 8._dp*(Trmd2YsIICYsIIYdadjYd) - 64._dp*(Trmd2YsIICYsIIYsIICYsII) -& 
&  8._dp*(Trmd2YsIICYsIIYzIIadjYzII) - 8._dp*(Trmd2YzIIadjYzIIYsIICYsII) -               & 
&  8._dp*(Trml2adjYzIIYsIICYsIIYzII) - 8._dp*(Trmq2adjYdYsIICYsIIYd) - 4._dp*(TrYdadjTdTsIICYsII) -& 
&  8._dp*(TrYdadjYdTsIICTsII) - 8._dp*(TrYdadjYdYsIICmd2CYsII) - 8*mHd2*TrYdadjYdYsIICYsII -& 
&  8*ms2*TrYdadjYdYsIICYsII - 8._dp*(TrYsIICmd2CYsIIYzIIadjYzII) - 2._dp*(TrYsIICTdTpTdCYsII) +& 
&  (8*g1p2*M1*TrYsIICTsII)/15._dp + (8*g3p2*M3*TrYsIICTsII)/3._dp - 8._dp*(TrYsIICTsIITdadjYd) -& 
&  28._dp*(TrYsIICTsIITsIICYsII) - 8._dp*(TrYsIICTsIITzIIadjYzII) - 2._dp*(TrYsIICTzIITpTzIICYsII) -& 
&  (8*g1p2*ms2*TrYsIICYsII)/15._dp - (8*g3p2*ms2*TrYsIICYsII)/3._dp - 6._dp*(TrYsIICYsIITdadjTd) -& 
&  36._dp*(TrYsIICYsIITsIICTsII) - 6._dp*(TrYsIICYsIITzIIadjTzII) - 32*ms2*TrYsIICYsIIYsIICYsII -& 
&  8*ms2*TrYsIICYsIIYzIIadjYzII - 8*mzz2*TrYsIICYsIIYzIIadjYzII - 4._dp*(TrYzIIadjTzIITsIICYsII) -& 
&  8._dp*(TrYzIIadjYzIITsIICTsII) + (8*g1p2*(2544*g1p2*M1 + 800*g3p2*M1 + 400*g3p2*M3 +  & 
&  15._dp*(TrCYsIITsII) - 30*M1*TrYsIICYsII)*Conjg(M1))/225._dp + (8*g3p2*(16*g1p2*M1 +  & 
&  32*g1p2*M3 + 870*g3p2*M3 + 3._dp*(TrCYsIITsII) - 6*M3*TrYsIICYsII)*Conjg(M3))/9._dp + & 
&  (80*g3p4*Tr2(3))/3._dp + (32*g1p2*Tr2U1(1,1))/15._dp - 16*g1*ooSqrt15*Tr3(1)

 
Dms2 = oo16pi2*( betams21 + oo16pi2 * betams22 ) 

 
Else 
Dms2 = oo16pi2* betams21 
End If 
 
 
!-------------------- 
! msb2 
!-------------------- 
 
betamsb21  = (-4*(8*AbsM1*g1p2 + 100*AbsM3*g3p2 - g1*sqrt15*Tr1(1)))/15._dp

 
 
If (TwoLoopRGE) Then 
betamsb22 = (16*(8*g1p2*(159*g1p2*M1 + 25*g3p2*(2._dp*(M1) + M3))*Conjg(M1) + 5*(5*g3p2*(435*g3p2*M3 +& 
&  8*g1p2*(M1 + 2._dp*(M3)))*Conjg(M3) + 75*g3p4*Tr2(3) + 6*g1p2*Tr2U1(1,1) +            & 
&  3*g1*sqrt15*Tr3(1))))/225._dp

 
Dmsb2 = oo16pi2*( betamsb21 + oo16pi2 * betamsb22 ) 

 
Else 
Dmsb2 = oo16pi2* betamsb21 
End If 
 
 
!-------------------- 
! mzz2 
!-------------------- 
 
betamzz21  = (-2*AbsM1*g1p2)/15._dp - 6*AbsM2*g2p2 - (32*AbsM3*g3p2)/3._dp +          & 
&  2._dp*(TrCTzIITpTzII) + 2._dp*(Trmd2YzIIadjYzII) + 2._dp*(Trml2adjYzIIYzII)           & 
&  + 2*mzz2*TrYzIIadjYzII + g1*ooSqrt15*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamzz22 = (2*AbsM2*g1p2*g2p2)/5._dp + 159*AbsM2*g2p4 + 32*AbsM2*g2p2*g3p2 + (4*g1p2*TrCTzIITpTzII)/5._dp -& 
&  (4*g1p2*M1*TrCTzIITpYzII)/5._dp - 4._dp*(TrCYsIITsIICTzIITpYzII) - 3._dp*(TrCYtIITpYzIICTzIITtII) -& 
&  4._dp*(Trmd2YdadjYdYzIIadjYzII) - 8._dp*(Trmd2YsIICYsIIYzIIadjYzII) - 2._dp*(Trmd2YzIIadjYeYeadjYzII) +& 
&  (4*g1p2*Trmd2YzIIadjYzII)/5._dp - 4._dp*(Trmd2YzIIadjYzIIYdadjYd) - 8._dp*(Trmd2YzIIadjYzIIYsIICYsII) -& 
&  20._dp*(Trmd2YzIIadjYzIIYzIIadjYzII) - 6._dp*(Trmd2YzIICYtIIYtIIadjYzII) -            & 
&  2._dp*(Trme2YeadjYzIIYzIIadjYe) - 2._dp*(Trml2adjYeYeadjYzIIYzII) - 4._dp*(Trml2adjYzIIYdadjYdYzII) -& 
&  8._dp*(Trml2adjYzIIYsIICYsIIYzII) + (4*g1p2*Trml2adjYzIIYzII)/5._dp - 2._dp*(Trml2adjYzIIYzIIadjYeYe) -& 
&  20._dp*(Trml2adjYzIIYzIIadjYzIIYzII) - 6._dp*(Trml2adjYzIIYzIICYtIIYtII) -            & 
&  6._dp*(Trml2CYtIIYtIIadjYzIIYzII) - 4._dp*(Trmq2adjYdYzIIadjYzIIYd) - 4._dp*(TrYdadjTdTzIIadjYzII) -& 
&  4._dp*(TrYdadjYdTzIIadjTzII) - 4*mHd2*TrYdadjYdYzIIadjYzII - 4*mzz2*TrYdadjYdYzIIadjYzII -& 
&  2._dp*(TrYeadjTzIITzIIadjYe) - 2._dp*(TrYeadjYzIITzIIadjTe) - 2*mHd2*TrYeadjYzIIYzIIadjYe -& 
&  2*mzz2*TrYeadjYzIIYzIIadjYe - 8._dp*(TrYsIICmd2CYsIIYzIIadjYzII) - 8._dp*(TrYsIICTsIITzIIadjYzII) -& 
&  2._dp*(TrYsIICTzIITpTzIICYsII) - 6._dp*(TrYsIICYsIITzIIadjTzII) - 8*ms2*TrYsIICYsIIYzIIadjYzII -& 
&  8*mzz2*TrYsIICYsIIYzIIadjYzII - 3._dp*(TrYtIIadjTzIITzIICYtII) - 6._dp*(TrYtIIadjYzIITzIICTtII) -& 
&  6*mt2*TrYtIIadjYzIIYzIICYtII - 6*mzz2*TrYtIIadjYzIIYzIICYtII - 6._dp*(TrYtIIadjYzIIYzIICYtIICml2) -& 
&  3._dp*(TrYtIICYtIITpTzIICTzII) - 2._dp*(TrYzIIadjTeTeadjYzII) - 4._dp*(TrYzIIadjTzIITdadjYd) -& 
&  4._dp*(TrYzIIadjTzIITsIICYsII) - 20._dp*(TrYzIIadjTzIITzIIadjYzII) - 2._dp*(TrYzIIadjYeTeadjTzII) +& 
&  (4*g1p2*mzz2*TrYzIIadjYzII)/5._dp - 4._dp*(TrYzIIadjYzIITdadjTd) - 8._dp*(TrYzIIadjYzIITsIICTsII) -& 
&  20._dp*(TrYzIIadjYzIITzIIadjTzII) - 20*mzz2*TrYzIIadjYzIIYzIIadjYzII - 6._dp*(TrYzIICTtIITtIIadjYzII) -& 
&  3._dp*(TrYzIICYtIITtIIadjTzII) + (g1p2*(1227*g1p2*M1 + 90*g2p2*M1 + 160*g3p2*M1 +     & 
&  80*g3p2*M3 + 45*g2p2*M2 - 180._dp*(TradjYzIITzII) + 360*M1*TrYzIIadjYzII)*Conjg(M1))/225._dp  
betamzz22 =  betamzz22+ (16*g3p2*(g1p2*(M1 + 2._dp*(M3)) + 15*(34*g3p2*M3 + 3*g2p2*(2._dp*(M3) +              & 
&  M2)))*Conjg(M3))/45._dp + (g1p2*g2p2*M1*Conjg(M2))/5._dp + 16*g2p2*g3p2*M3*Conjg(M2) +& 
&  6*g2p4*Tr2(2) + (32*g3p4*Tr2(3))/3._dp + (2*g1p2*Tr2U1(1,1))/15._dp + 4*g1*ooSqrt15*Tr3(1)

 
Dmzz2 = oo16pi2*( betamzz21 + oo16pi2 * betamzz22 ) 

 
Else 
Dmzz2 = oo16pi2* betamzz21 
End If 
 
 
!-------------------- 
! mzb2 
!-------------------- 
 
betamzb21  = (-2*AbsM1*g1p2)/15._dp - 6*AbsM2*g2p2 - (32*AbsM3*g3p2)/3._dp -          & 
&  g1*ooSqrt15*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamzb22 = (g1p2*(1227*g1p2*M1 + 5*(16*g3p2*(2._dp*(M1) + M3) + 9*g2p2*(2._dp*(M1) +             & 
&  M2)))*Conjg(M1) + 5*(16*g3p2*(g1p2*(M1 + 2._dp*(M3)) + 15*(34*g3p2*M3 +               & 
&  3*g2p2*(2._dp*(M3) + M2)))*Conjg(M3) + 3*(3*g2p2*(795*g2p2*M2 + g1p2*(M1 +            & 
&  2._dp*(M2)) + 80*g3p2*(M3 + 2._dp*(M2)))*Conjg(M2) + 2*(45*g2p4*Tr2(2) +              & 
&  80*g3p4*Tr2(3) + g1*(g1*Tr2U1(1,1) - 2*sqrt15*Tr3(1))))))/225._dp

 
Dmzb2 = oo16pi2*( betamzb21 + oo16pi2 * betamzb22 ) 

 
Else 
Dmzb2 = oo16pi2* betamzb21 
End If 
 
 
!-------------------- 
! M1 
!-------------------- 
 
betaM11  = (136*g1p2*M1)/5._dp

 
 
If (TwoLoopRGE) Then 
betaM12 = (2*g1p2*(3004*g1p2*M1 + 2610*g2p2*M1 + 4600*g3p2*M1 + 4600*g3p2*M3 + 2610*g2p2*M2 +   & 
&  210._dp*(TradjYdTd) + 270._dp*(TradjYeTe) + 390._dp*(TradjYuTu) + 210._dp*(TradjYzIITzII) +& 
&  360._dp*(TrCYsIITsII) + 405._dp*(TrCYtIITtII) - 210*M1*TrYdadjYd - 270*M1*TrYeadjYe - & 
&  360*M1*TrYsIICYsII - 405*M1*TrYtIICYtII - 390*M1*TrYuadjYu - 210*M1*TrYzIIadjYzII -   & 
&  405*Conjg(L1II)*(L1II*M1 - TL1II) - 405*Conjg(L2II)*(L2II*M1 - TL2II)))/75._dp

 
DM1 = oo16pi2*( betaM11 + oo16pi2 * betaM12 ) 

 
Else 
DM1 = oo16pi2* betaM11 
End If 
 
 
!-------------------- 
! M2 
!-------------------- 
 
betaM21  = 16*g2p2*M2

 
 
If (TwoLoopRGE) Then 
betaM22 = (2*g2p2*(58*g1p2*M1 + 200*g3p2*M3 + 58*g1p2*M2 + 940*g2p2*M2 + 200*g3p2*M2 +          & 
&  30._dp*(TradjYdTd) + 10._dp*(TradjYeTe) + 30._dp*(TradjYuTu) + 30._dp*(TradjYzIITzII) +& 
&  35._dp*(TrCYtIITtII) - 30*M2*TrYdadjYd - 10*M2*TrYeadjYe - 35*M2*TrYtIICYtII -        & 
&  30*M2*TrYuadjYu - 30*M2*TrYzIIadjYzII - 35*Conjg(L1II)*(L1II*M2 - TL1II) -            & 
&  35*Conjg(L2II)*(L2II*M2 - TL2II)))/5._dp

 
DM2 = oo16pi2*( betaM21 + oo16pi2 * betaM22 ) 

 
Else 
DM2 = oo16pi2* betaM21 
End If 
 
 
!-------------------- 
! M3 
!-------------------- 
 
betaM31  = 8*g3p2*M3

 
 
If (TwoLoopRGE) Then 
betaM32 = (2*g3p2*(23*g1p2*M1 + 23*g1p2*M3 + 45*g2p2*M3 + 800*g3p2*M3 + 45*g2p2*M2 +            & 
&  12._dp*(TradjYdTd) + 12._dp*(TradjYuTu) + 12._dp*(TradjYzIITzII) + 27._dp*(TrCYsIITsII) -& 
&  12*M3*TrYdadjYd - 27*M3*TrYsIICYsII - 12*M3*TrYuadjYu - 12*M3*TrYzIIadjYzII))/3._dp

 
DM3 = oo16pi2*( betaM31 + oo16pi2 * betaM32 ) 

 
Else 
DM3 = oo16pi2* betaM31 
End If 
 
 
!-------------------- 
! WOp 
!-------------------- 
 
betaWOp1  = TpYeCYeWOp + 3._dp*(TpYzIICYzIIWOp) - (6*(-5._dp*(AbsL2II) +              & 
&  g1p2 + 5._dp*(g2p2) - 5._dp*(TrYuadjYu))*WOp)/5._dp + WOpadjYeYe + 3._dp*(WOpadjYzIIYzII)& 
&  + 3._dp*(WOpCYtIIYtII) + 3._dp*(YtIICYtIIWOp)

 
 
If (TwoLoopRGE) Then 
betaWOp2 = (417._dp*(g1p4)/25._dp + (18*g1p2*g2p2)/5._dp + 57._dp*(g2p4) - 24*CL2IIp2*L2IIp2 -   & 
&  6._dp*(TrYdadjYuYuadjYd) + (6*AbsL2II*(6._dp*(g1p2) + 20._dp*(g2p2) - 15._dp*(TrYuadjYu)))/5._dp +& 
&  (8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu)/5._dp - 18._dp*(TrYuadjYuYuadjYu))*WOp +         & 
&  (-3._dp*(AbsL1II) + 6._dp*(g1p2)/5._dp - 3._dp*(TrYdadjYd) - TrYeadjYe)*WOpadjYeYe +  & 
&  (-10._dp*(TpYeCYeTpYeCYeWOp) - 15*AbsL1II*TpYeCYeWOp + 6*g1p2*TpYeCYeWOp -            & 
&  30._dp*(TpYzIICYdTpYdCYzIIWOp) - 60._dp*(TpYzIICYsIIYsIICYzIIWOp) - 30._dp*(TpYzIICYzIITpYzIICYzIIWOp) -& 
&  2*g1p2*TpYzIICYzIIWOp + 80*g3p2*TpYzIICYzIIWOp - 15*TpYeCYeWOp*TrYdadjYd -            & 
&  5*TpYeCYeWOp*TrYeadjYe - 15*TpYzIICYzIIWOp*TrYzIIadjYzII - 10._dp*(WOpadjYeYeadjYeYe) -& 
&  30._dp*(WOpadjYzIIYdadjYdYzII) - 60._dp*(WOpadjYzIIYsIICYsIIYzII) + (-2._dp*(g1p2) +  & 
&  80._dp*(g3p2) - 15._dp*(TrYzIIadjYzII))*WOpadjYzIIYzII - 30._dp*(WOpadjYzIIYzIIadjYzIIYzII) -& 
&  15._dp*(WOpCYtIITpYeCYeYtII) - 45._dp*(WOpCYtIITpYzIICYzIIYtII) + 3*(-5._dp*(AbsL1II) +& 
&  6._dp*(g1p2) + 20._dp*(g2p2) - 5._dp*(TrYtIICYtII))*WOpCYtIIYtII - 45._dp*(WOpCYtIIYtIICYtIIYtII) -& 
&  15._dp*(YtIIadjYeYeCYtIIWOp) - 45._dp*(YtIIadjYzIIYzIICYtIIWOp) - 15*AbsL1II*YtIICYtIIWOp +& 
&  18*g1p2*YtIICYtIIWOp + 60*g2p2*YtIICYtIIWOp - 15*TrYtIICYtII*YtIICYtIIWOp -           & 
&  45._dp*(YtIICYtIIYtIICYtIIWOp))/5._dp

 
DWOp = oo16pi2*( betaWOp1 + oo16pi2 * betaWOp2 ) 

 
Else 
DWOp = oo16pi2* betaWOp1 
End If 
 
!-------------------------------------------------------------------------------
! these matrices are hermitian but numerical effects induce non-hermiatian parts
! which need to be corrected
!-------------------------------------------------------------------------------
Dmd2 = 0.5_dp * ( Dmd2 + Transpose(Conjg(Dmd2)) )
Dme2 = 0.5_dp * ( Dme2 + Transpose(Conjg(Dme2)) )
Dml2 = 0.5_dp * ( Dml2 + Transpose(Conjg(Dml2)) )
Dmq2 = 0.5_dp * ( Dmq2 + Transpose(Conjg(Dmq2)) )
Dmu2 = 0.5_dp * ( Dmu2 + Transpose(Conjg(Dmu2)) )

   Call Chop(Dmu)
   Call Chop(DBmu)
   Call Chop(DMTII)
   Call Chop(DBMTII)
   Call Chop(DMSII)
   Call Chop(DBMSII)
   Call Chop(DMZII)
   Call Chop(DBMZII)

Call ParametersToG365(Dg1,Dg2,Dg3,DYu,DYd,DYe,DYtII,DYsII,DYzII,DL1II,DL2II,          & 
& DMu,DMTII,DMZII,DMSII,DTu,DTd,DTe,DTtII,DTsII,DTzII,DTL1II,DTL2II,DBmu,DBMTII,         & 
& DBMZII,DBMSII,Dmq2,Dml2,DmHd2,DmHu2,Dmd2,Dmu2,Dme2,Dmt2,Dmtb2,Dms2,Dmsb2,              & 
& Dmzz2,Dmzb2,DM1,DM2,DM3,DWOp,f)

Iname = Iname - 1 
 
End Subroutine rge365  

#endif SEESAWIII


 Subroutine Set_Decoupling_Heavy_States(set)
  Implicit none
  Logical, Intent(in) :: set

   decoupling_heavy_states = set

 End Subroutine Set_Decoupling_Heavy_States


End Module RGEs

