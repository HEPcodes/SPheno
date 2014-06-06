Module LoopCouplings

Use Control
Use Mathematics
Use LoopFunctions
Use RGEs
Use StandardModel
Use Couplings

Logical, Save :: MZ_input = .False.
 Real(dp), Private :: vevSM(2)

Contains

 Real(dp) Function AlphaEwDR(Q, mHpm, mSqU, mSqD, mSl, mChar)
 !-----------------------------------------------------------------------
 ! This function calculates the electromagnetig alpha(Q) in the DR scheme.
 ! The formula is taken from J. Bagger et al., Nucl.Phys.B 491, 3, (1997)
 ! The input is:
 !  - Alpha at Q=0, something like 1/137
 !  - Q the energy
 !  - mW, the W-boson mass
 !  - mT, the top mass
 !  - mHpm, the mass of the charged Higgs boson
 !  - mSqU(i), i=1,..,6 the masses of the u-squarks
 !  - mSqD(i), i=1,..,6 the masses of the d-squarks
 !  - mSl(i), i=1,..,6 the masses of the sleptons
 !  - mChar(i), i=1,2 the masses of the charginos
 ! written by Werner Porod, 19.7.1999
 ! 16.11.01: - writting it such that R-parity violation is included
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: Q, mHpm(:), mSqU(:), mSqD(:), mSl(:), mChar(:)

  Integer :: i, n_C, n_H
  Real(dp) :: DeltaAlpha,sumI(4)

  n_H = Size( mHpm )
  n_C = Size( mChar )

  sumI = 0._dp
  Do i=1,6
   sumI(1) = sumI(1) + Log( mSqU(i) / Q )
   sumI(2) = sumI(2) + Log( mSqD(i) / Q )
  End Do

  Do i=2,n_H
   sumI(3) = sumI(3) + Log( mHpm(i) / Q )
  End Do
  Do i=1,(5-n_C)*2
   sumI(3) = sumI(3) + Log( mSl(i) / Q )
  End Do
  Do i=n_C-1,n_C
   sumI(4) = sumI(4) + Log( mChar(i) / Q )
  End Do

  If (MZ_input) then
   DeltaAlpha = 1._dp - Alpha / Alpha_MZ_MS  ! MSbar value ^= mW + light fermions
   DeltaAlpha = DeltaAlpha + alpha / (6._dp * Pi) ! conversion to DRbar
   DeltaAlpha = DeltaAlpha - ( 16._dp * Log(mf_u(3) / Q ) / 9._dp &
            & + sumI(3) / 3._dp + (4._dp * sumI(1) + sumI(2) ) / 9._dp    &
            & + 4._dp * sumI(4) / 3._dp ) * Alpha / (2._dp * pi) ! SUSY contr.
  else
   DeltaAlpha = 7._dp * Log(Q / mW) + 16._dp * Log(mf_u(3) / Q ) / 9._dp &
            & + sumI(3) / 3._dp + (4._dp * sumI(1) + sumI(2) ) / 9._dp    &
            & + 4._dp * sumI(4) / 3._dp
   DeltaAlpha = Delta_Alpha_Lepton + Delta_Alpha_Hadron &
            & - alpha * DeltaAlpha / ( 2._dp * Pi)
  End If

  AlphaEwDR = Alpha / (1._dp - DeltaAlpha)

 End Function AlphaEwDR


 Real(dp) Function AlphaEwMS(Q, alphaDR, mHpm, mSqU, mSqD, mSl, mChar, mt)
 !-----------------------------------------------------------------------
 ! This function calculates the electromagnetig alpha(Q) in the MS scheme
 ! starting from the DRbar value
 ! The formula is taken from J. Bagger et al., Nucl.Phys.B 491, 3, (1997)
 ! The input is:
 !  - Alpha at Q=0, something like 1/137
 !  - Q the energy
 !  - mW, the W-boson mass
 !  - mT, the top mass
 !  - mHpm, the mass of the charged Higgs boson
 !  - mSqU(i), i=1,..,6 the masses of the u-squarks
 !  - mSqD(i), i=1,..,6 the masses of the d-squarks
 !  - mSl(i), i=1,..,6 the masses of the sleptons
 !  - mChar(i), i=1,2 the masses of the charginos
 !  - mt, top-quark mass
 ! written by Werner Porod, 05.02.2010
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: Q, alphaDR, mHpm(:), mSqU(:), mSqD(:), mSl(:) &
                      & , mChar(:), mt

  Integer :: i, n_C, n_H
  Real(dp) :: DeltaAlpha,sumI(4)

  n_H = Size( mHpm )
  n_C = Size( mChar )

  sumI = 0._dp
  Do i=1,6
   sumI(1) = sumI(1) + Log( mSqU(i) / Q )
   sumI(2) = sumI(2) + Log( mSqD(i) / Q )
  End Do

  Do i=2,n_H
   sumI(3) = sumI(3) + Log( mHpm(i) / Q )
  End Do
  Do i=1,(5-n_C)*2
   sumI(3) = sumI(3) + Log( mSl(i) / Q )
  End Do
  Do i=n_C-1,n_C
   sumI(4) = sumI(4) + Log( mChar(i) / Q )
  End Do

  DeltaAlpha = 1._dp / (6._dp * Pi) ! conversion from DRbar
  DeltaAlpha = DeltaAlpha - ( 16._dp * Log(mt / Q ) / 9._dp &
            & + sumI(3) / 3._dp + (4._dp * sumI(1) + sumI(2) ) / 9._dp    &
            & + 4._dp * sumI(4) / 3._dp )  / (2._dp * pi) ! SUSY contr.

  AlphaEwMS = AlphaDR / (1._dp + AlphaDR * DeltaAlpha)

 End Function AlphaEwMS


 Real(dp) Function Alpha_MSbar(Q, mW, mt)
 Implicit None
  Real(dp), Intent(in) :: Q, mW
  Real(dp), Intent(in), Optional :: mt

  Real(dp) :: DeltaAlpha

  If (MZ_input) Then
   Alpha_MSbar = Alpha_mZ_MS
   If (Present(mt)) Then
    DeltaAlpha = - 8._dp * Log(Q / mt) / (9._dp * Pi)
    Alpha_MSbar = Alpha_MSbar / ( 1._dp + DeltaAlpha * alpha)
   End If 
  Else
   DeltaAlpha = 3.5_dp * Log(Q / mW) / Pi + 0.5_dp * oo3pi
   If (Present(mt)) DeltaAlpha = DeltaAlpha - 8._dp * Log(Q / mt) / (9._dp * Pi)
   Alpha_MSbar = Alpha / ( 1._dp - Delta_Alpha_Lepton - Delta_Alpha_Hadron  &
              &         + DeltaAlpha * alpha)
   Alpha_MZ_MS = Alpha_MSbar
  End If

 End Function Alpha_MSbar


 Real(dp) Function AlphaSDR(Q, mG, mSqU, mSqD, mt_in)
 !-----------------------------------------------------------------------
 ! This function calculates the strong coupling alpha_s(Q) in the DR scheme.
 ! The formula is taken from J. Bagger et al., Nucl.Phys.B
 ! The input is:
 !  - Q the energy
 !  - mG, the gluino mass
 !  - mT_in, the top mass
 !  - mSqU(i), i=1,..,6 the masses of the u-squarks
 !  - mSqD(i), i=1,..,6 the masses of the d-squarks
 ! written by Werner Porod, 19.7.1999
 ! 16.11.01: portation to f90 
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: Q,mG,mSqU(6),mSqD(6)
  Real(dp), Intent(in), Optional :: mt_in

  Integer :: i
  Real(dp) :: DeltaAlpha, sumI, mt

  sumI = 0._dp
  Do i=1,6
   sumI = sumI + Log( mSqU(i) / Q ) + Log( mSqD(i) / Q )
  Enddo

  If (Present(mt_in)) Then
   mt = mt_in
  Else
   mt = mf_u(3)
  End If

  DeltaAlpha = 0.5_dp - 2._dp * Log(mG /Q) - 2._dp * Log(mt/Q ) / 3._dp &
           & - sumI / 6._dp
  DeltaAlpha = AlphaS_mZ * DeltaAlpha / ( 2._dp * Pi)
  
  AlphaSDR = AlphaS_mZ / (1._dp - DeltaAlpha)

 End Function AlphaSDR

 Subroutine CoupPseudoScalarGluon(m_H2, mf_u2, g_u, mf_d2, g_d, coup, coupSM)
 !----------------------------------------------------------------------------
 ! calculates the lowest order coupling for the radiative decay of a pseudoscalar
 ! to two gluons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! based on  CoupScalarGluon by Werner Porod
 ! written by Florian Staub, 18.06.2010
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: m_H2, mf_u2(3), mf_d2(3)
  Complex(dp), Intent(in) :: g_u(3), g_d(3)
  Complex(dp), Intent(out) :: coup
  Complex(dp), Optional, Intent(out) :: coupSM
 
  Integer :: i1
  Real(dp) :: Mh2p
  !--------------------------------------
  ! W-boson contribution
  !--------------------------------------
  mH2p = 0.25_dp * m_H2

  If (Present(coupSM)) Then
   coupSM = 0._dp
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coupSM = coupSM + AP_onehalf(mH2p/ mf_d2(i1))+ AP_onehalf(mH2p/ mf_u2(i1))
   End Do
   coupSM = 0.75_dp * coupSM
  End If

  coup = 0._dp
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coup = coup + g_d(i1) * AP_onehalf(mH2p/ mf_d2(i1))  &
         &      + g_u(i1) * AP_onehalf(mH2p/ mf_u2(i1))
   End Do
  
  coup = 0.75_dp * coup
 
 End Subroutine CoupPseudoScalarGluon


 Subroutine CoupPseudoScalarPhoton(m_H2,mf_u2, g_u, mf_d2, g_d, mf_l2    &
    & , g_l, mC2, g_C, coup, coupSM)
 !----------------------------------------------------------------------------
 ! calculates the lowest order coupling for the radiative decay of a pseudoscalar
 ! to two photons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! based on  CoupScalarPhoton by Werner Porod
 ! written by Florian Staub, 18.06.2010
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: m_H2,  mf_u2(3), mf_d2(3), mf_l2(3), mC2(2)
  Complex(dp), Intent(in) :: g_u(3), g_d(3), g_l(3), g_C(2)
  Complex(dp), Intent(out) :: coup
  Complex(dp), Optional, Intent(out) :: coupSM
 
  Integer :: i1
  Real(dp) :: Mh2p
  !--------------------------------------
  ! W-boson contribution
  !--------------------------------------
  mH2p = 0.25_dp * m_H2

  If (Present(coupSM)) Then
   coupSM = A_one(mH2p / mW2)
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coupSM = coupSM + AP_onehalf(mH2p/ mf_l2(i1)) + (AP_onehalf(mH2p/ mf_d2(i1)) &
                                + 4._dp * AP_onehalf(mH2p/ mf_u2(i1)) ) / 3._dp
   End Do
  End If

  coup = 0._dp
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coup = coup + g_l(i1) * AP_onehalf(mH2p/ mf_l2(i1))    &
         &      + ( g_d(i1) * AP_onehalf(mH2p/ mf_d2(i1))  &
         &        + 4._dp * g_u(i1) * AP_onehalf(mH2p/ mf_u2(i1)) ) / 3._dp
   End Do
  !--------------------------------------
  ! chargino contributions
  !--------------------------------------
  Do i1=1,2
   coup = coup + g_C(i1) * AP_onehalf(mH2p/ mC2(i1))
  End Do
  
  
 End Subroutine CoupPseudoScalarPhoton


 Subroutine CoupScalarGluon(m_H2, mf_u2, g_u, mf_d2, g_d, mSup2, g_su &
                          & , mSdown2, g_sd, coup, coupSM)
 !----------------------------------------------------------------------------
 ! calculates the lowest order coupling for the radiative decay of a scalar
 ! to two gluons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! written by Werner Porod, 3.1.2009
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: m_H2, mf_u2(3), mf_d2(3), mSup2(6), mSdown2(6)
  Complex(dp), Intent(in) :: g_u(3), g_d(3), g_su(6), g_sd(6)
  Complex(dp), Intent(out) :: coup
  Complex(dp), Optional, Intent(out) :: coupSM
 
  Integer :: i1
  Real(dp) :: Mh2p
  !--------------------------------------
  ! W-boson contribution
  !--------------------------------------
  mH2p = 0.25_dp * m_H2

  If (Present(coupSM)) Then
   coupSM = 0._dp
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coupSM = coupSM + A_onehalf(mH2p/ mf_d2(i1))+ A_onehalf(mH2p/ mf_u2(i1))
   End Do
   coupSM = 0.75_dp * coupSM
  End If

  coup = 0._dp
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coup = coup + g_d(i1) * A_onehalf(mH2p/ mf_d2(i1))  &
         &      + g_u(i1) * A_onehalf(mH2p/ mf_u2(i1))
   End Do
  !--------------------------------------
  ! Squark contributions
  !--------------------------------------
  Do i1=1,6
   coup = coup + g_sd(i1) * A_zero(mH2p/ mSdown2(i1)) &
        &      + g_su(i1) * A_zero(mH2p/ mSup2(i1))
  End Do
  
  coup = 0.75_dp * coup
 
 End Subroutine CoupScalarGluon

 Subroutine CoupScalarPhoton(m_H2, mW2, g_W, mf_u2, g_u, r_T, mf_d2, g_d, mf_l2 &
    & , g_l, mC2, g_C, m_Hp2, g_Hp, mSup2, g_su, mSdown2, g_sd, r_sq            &
    & , mSlept2, g_sl, coup, coupSM)
 !----------------------------------------------------------------------------
 ! calculates the lowest order coupling for the radiative decay of a scalar
 ! to two photons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! written by Werner Porod, 27.12.2008
 ! 8.6.2012: adding approximate QCD corrections based on above paper,
 !           using form factors which are however only useful for light Higgs
 !           which should be set to 1 for heavier Higgs bosons
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: m_H2, mW2, mf_u2(3), mf_d2(3), mf_l2(3), m_Hp2 &
    & , mSup2(6), mSdown2(6), mSlept2(6), mC2(2)
  Real(dp), Intent(in) :: g_W, g_Hp, r_T, r_sq 
  Complex(dp), Intent(in) :: g_u(3), g_d(3), g_l(3), g_C(2), g_su(6), g_sd(6) &
    & , g_sl(6)
  Complex(dp), Intent(out) :: coup
  Complex(dp), Optional, Intent(out) :: coupSM
 
  Integer :: i1
  Real(dp) :: Mh2p
  !--------------------------------------
  ! W-boson contribution
  !--------------------------------------
  mH2p = 0.25_dp * m_H2

  If (Present(coupSM)) Then
   coupSM = A_one(mH2p / mW2)
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coupSM = coupSM + A_onehalf(mH2p/ mf_l2(i1)) &
         & + ( A_onehalf(mH2p/ mf_d2(i1))  &
         &   + r_T * 4._dp * A_onehalf(mH2p/ mf_u2(i1)) ) / 3._dp
   End Do
  End If

  coup = g_W * A_one(mH2p / mW2)
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coup = coup + g_l(i1) * A_onehalf(mH2p/ mf_l2(i1))   &
        &      + ( g_d(i1) * A_onehalf(mH2p/ mf_d2(i1))  &
        &        + r_T * 4._dp * g_u(i1) * A_onehalf(mH2p/ mf_u2(i1)) ) / 3._dp
   End Do
  !--------------------------------------
  ! chargino contributions
  !--------------------------------------
  Do i1=1,2
   coup = coup + g_C(i1) * A_onehalf(mH2p/ mC2(i1))
  End Do
  !--------------------------------------
  ! charged Higgs contributions
  !--------------------------------------
  coup = coup + g_Hp * A_zero(mH2p/ m_Hp2)
  !--------------------------------------
  ! Sfermion contributions
  !--------------------------------------
  Do i1=1,6
   coup = coup + g_sl(i1) * A_zero(mH2p/ mSlept2(i1))          &
        &      + r_sq * ( g_sd(i1) * A_zero(mH2p/ mSdown2(i1)) &
        &               + 4._dp * g_su(i1) * A_zero(mH2p/ mSup2(i1)) ) / 3._dp
  End Do
  
 End Subroutine CoupScalarPhoton

 Subroutine InitializeLoopCouplings(vevs)
 implicit none
  Real(dp), intent(in) :: vevs(2)

  vevSM = vevs
 end Subroutine InitializeLoopCouplings

 Subroutine RunningCouplings(Q,g,yuk)
 !-----------------------------------------------------------------------
 ! calculates running gauge couplings and the Yukawa couplings of the
 ! third generation to two-loop order
 ! written by Werner Porod, 6.10.2000
 ! - correct change to DRbar scheme is missing, using On-shell and MSbar
 !   values as input
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: Q
  Real(dp), Intent(out) :: g(3), yuk(3)

  Integer :: kont
  Real(dp) :: g1(6), tz, dt, cosW2, sinW2

  g1 = 0._dp
  cosW2 = mW2 / mZ2
  sinW2 = 1._dp - cosW2
  !----------------------------------
  ! running from Z-scale to scale Q
  !----------------------------------
  g1(1) = Sqrt( 20._dp*pi*alpha_mZ/(3._dp*cosW2) )
  g1(2) = Sqrt( 4._dp*pi*alpha_mZ/sinW2)
  g1(3) = Sqrt( 4._dp * pi * alphaS_mZ )

  g1(4) = sqrt2 * mf_l(3) / vevSM(1)
  g1(5) = sqrt2 * mf_d(3) / vevSM(1)
  g1(6) = sqrt2 * mf_u(3) / vevSM(2)

  tz = Log(mZ/Q)
  dt = - tz / 50._dp

  Call odeint(g1, 6, tz, 0._dp, 1.e-5_dp, dt, 0._dp, rge6, kont)

  g1(1) = Sqrt(3._dp/5._dp) * g1(1)

  g = g1(1:3)
  yuk = g1(4:6)

 End Subroutine RunningCouplings


 Real(Dp) Function SUALFE(QS)
 !----------------------------------------------------------------------
 !     Returns the running EM coupling alpha_em(q**2)
 !     taken from ISASUSY
 !  is used to calculate the coupling at low energies 
 !
 ! SEE BARGER/PHILLIPS, P. 202 
 !----------------------------------------------------------------------
 Implicit None

  Real(Dp), Intent(in) :: qs
  Integer :: i1
  Real(Dp) :: sum, qu2, qd2, qs4

  SUM=0._dp
  QU2=4._dp/3._dp
  QD2=1._dp/3._dp
  qs4 = 0.25_dp * QS 
  Do i1=1,3
   If (qs4.Gt.mf_l2(i1)) SUM=SUM+Log(qs4/mf_l2(i1) )
   If (qs4.Gt.mf_u2(i1)) SUM=SUM+QU2*Log(qs4/mf_u2(i1) )
   If (qs4.Gt.mf_d2(i1)) SUM=SUM+QD2*Log(qs4/mf_d2(i1) )
  End Do

  SUALFE = Alpha/(1._dp- Alpha * sum /(3._dp*PI) )

 End Function SUALFE

End Module LoopCouplings
