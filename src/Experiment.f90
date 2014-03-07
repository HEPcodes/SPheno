Module Experiment

Use Control
Use StandardModel
Use Couplings
Use Model_Data
Use MSSM_Data

Integer, Public :: check_exp(11)      

Real(dp), Public, Dimension(2) ::  m_sq_atm, tan_sq_atm, m_sq_sol, tan_sq_sol &
 & , Ue3_sq

! neutrino data taken from T.~Schwetz, M.~Tortola and J.~W.~F.~Valle
!                          arXiv:1103.0734
! correspond to 1 sigma but for U^2_e3 which is at 90% CL
! Min/Max of atm. mass squared
Real(dp), Public :: m2_atm_min = 2.36e-21_dp, m2_atm_max = 2.54e-21_dp
! best fit for atm. mass squared
Real(dp), Public :: m2_atm = 2.45e-21_dp
! Min/Max of tan^2(theta_atm)
Real(dp), Public :: tan2_atm_min = 0.8182_dp , tan2_atm_max = 1.3256_dp
! Min/Max of solar mass squared
Real(dp), Public :: m2_sol_min = 7.46e-23_dp, m2_sol_max = 7.83e-23_dp
! best fit for solar mass squared    
Real(dp), Public :: m2_sol = 7.64e-23_dp 
! Min/Max of tan^2(theta_sol)
Real(dp), Public :: tan2_sol_min = 0.4286_dp, tan2_sol_max = 0.4970_dp
! Min/Max (U_e3)^2
Real(dp), Public :: Ue32_min = 0._dp, Ue32_max = 0.035_dp
Real(dp), Public :: mC_min, mGlu_min           ! lower bound on chargino, gluinos
Real(dp), Public :: mSquark_min, mStop_min     ! lower bound on 1st gen. ~q, ~t_1

Contains


 Subroutine CheckBounds(i_mod)
 !-----------------------------------------------------------------
 ! This subroutine checks experimental bounds on SUSY particles
 ! 14.03.04: -implementing Higgs bound, ALEPH, PLB xx.xx.xx
 !-----------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_mod    ! defines model: 1 -> MSSM
                                  !                2 -> NMSSM
                                  !                3 -> RP explicit
!                                  !                4 -> RP spontaneous

  Integer :: i1, i_h
  Real(dp) :: sin2_ab(4), mhi(4), vev, mSel, mChar, m_bound, dm, tan2
  !----------------------------------
  ! initialisation
  !----------------------------------
  check_exp = 0 

  If (i_mod.Eq.1) Then
   sin2_ab(1) = (RS0(1,1) + RS0(1,2) * tanb)**2 / (1._dp + tanb**2)
   mhi(1) = mS0(1)
   i_h = 1
   mChar = mC(1)
   mSel = mSlepton(1)

  Else If (i_mod.Eq.3) Then
   mChar = mC5(4)

! neutrino physics
   Dm = mN7(3)**2 - mN7(2)**2
   If ((Dm.Lt.m2_atm_min).Or.(Dm.Gt.m2_atm_max)) check_exp(7) = 1
   Dm = mN7(2)**2 - mN7(1)**2
   If ((Dm.Lt.m2_sol_min).Or.(Dm.Gt.m2_sol_max)) check_exp(8) = 1

   If (Abs(N7(3,7)).Eq.0._dp) Then
    tan2 = 1.e99_dp
   Else 
    tan2 = Abs(N7(3,6)/N7(3,7))**2
   End If
   If ((tan2.Lt.tan2_atm_min).Or.(tan2.Gt.tan2_atm_max)) check_exp(9) = 1

   If (Abs(N7(1,5)).Eq.0._dp) Then
    tan2 = 1.e99_dp
   Else 
    tan2 = Abs(N7(2,5)/N7(1,5))**2
   End If
   If ((tan2.Lt.tan2_sol_min).Or.(tan2.Gt.tan2_sol_max)) check_exp(10) = 1

   tan2 = Abs(N7(3,5))**2
   If (tan2.Gt.Ue32_max) check_exp(11) = 1

!  Else If (i_mod.Eq.4) Then
!   i_h = 4
!   mhi = mS08(1:4)
!   !--------------------------------------------------------
!   ! vev = vev_SM / 2
!   ! the input in CoupScalarZ gives the reduced coupling
!   !--------------------------------------------------------
!   vev = mW / gauge_mZ(2)
!   Do i1=1,4
!    Call CoupScalarZ(i1, 1._dp, vev, vevSM, vevL, RS08, sin2_ab(i1) )
!   End Do
!   mChar = mC5(4)
!   mSel = 150._dp ! to be improved

! neutrino physics
!   Dm = mN10(3)**2 - mN10(2)**2
!   If ((Dm.Lt.m2_atm_min).Or.(Dm.Gt.m2_atm_max)) check_exp(7) = 1
!   Dm = mN10(2)**2 - mN10(1)**2
!   If ((Dm.Lt.m2_sol_min).Or.(Dm.Gt.m2_sol_max)) check_exp(8) = 1
!
!   If (Abs(N10(3,7)).Eq.0._dp) Then
!    tan2 = 1.e99_dp
!   Else 
!    tan2 = Abs(N10(3,6)/N10(3,7))**2
!   End If
!   If ((tan2.Lt.tan2_atm_min).Or.(tan2.Gt.tan2_atm_max)) check_exp(9) = 1
!
!   If (Abs(N10(1,5)).Eq.0._dp) Then
!    tan2 = 1.e99_dp
!   Else 
!    tan2 = Abs(N10(2,5)/N10(1,5))**2
!   End If
!   If ((tan2.Lt.tan2_sol_min).Or.(tan2.Gt.tan2_sol_max)) check_exp(10) = 1
!
!   tan2 = Abs(N10(3,5))**2
!   If (tan2.Gt.Ue32_max) check_exp(11) = 1
!
  End If

  !----------------------------------
  ! Higgs mass bound
  !----------------------------------
  Do i1=1,i_h
   If (sin2_ab(i1).Lt.0.01_dp) Then
    m_bound = 2500._dp * sin2_ab(i1)
   Else If  (sin2_ab(i1).Lt.0.05_dp) Then
    m_bound = 25._dp + 1100._dp * sin2_ab(i1)
   Else If (sin2_ab(i1).Lt.0.1_dp) Then
    m_bound = 80._dp + 138._dp  * sin2_ab(i1)
   Else If (sin2_ab(i1).Lt.0.2_dp) Then
    m_bound = 86.4_dp + 74._dp * sin2_ab(i1)
   Else If (sin2_ab(i1).Lt.0.3_dp) Then
    m_bound = 89.2_dp + 60._dp * sin2_ab(i1)
   Else If (sin2_ab(i1).Lt.0.4_dp) Then
    m_bound = 98.8_dp + 28._dp * sin2_ab(i1)
   Else If (sin2_ab(i1).Lt.0.5_dp) Then
    m_bound = 104._dp + 15._dp * sin2_ab(i1)
   Else If (sin2_ab(i1).Lt.0.5_dp) Then
    m_bound = 104._dp + 15._dp * sin2_ab(i1)
   Else If (sin2_ab(i1).Lt.0.6_dp) Then
    m_bound = 108._dp + 7._dp * sin2_ab(i1)
   Else If (sin2_ab(i1).Lt.0.8_dp) Then
    m_bound = 109._dp + 5._dp * sin2_ab(i1)
   Else 
    m_bound = 110._dp + 4._dp * sin2_ab(i1)
   End If
 ! including theoretical uncertainty
   If (m_bound.Gt.(mhi(i1)+3._dp)) check_exp(1) = 1
  End Do

! SUSY searches
  If (mChar.Lt.101._dp) check_exp(2) = 1    ! LEP
  If (mSel.Lt.95._dp) check_exp(3) = 1   ! LEP
  If (mSup(5).Lt.100._dp) check_exp(4) = 1 ! Tevatron
  If (mGlu.Lt.200._dp) check_exp(5) = 1  ! Tevatron


 End Subroutine CheckBounds


End Module Experiment

