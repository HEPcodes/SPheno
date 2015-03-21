
Module Model_Data
!-------------------------------------------------------------
! This module contains global definition of MSSM parameters,
! SUSY and Higgs masses, mixing matrices, decay widths and
! production cross sections.
! The SUSY parameters (denoted generically xxx) are given at four scales:
! mZ ............... xxx_mZ
! Q_EWSB ........... xxx
! M_GUT, M_GMSB .... xxx_0
! M_nu_R ........... xxx_mR
! Some of them may not be given at all scales
!-------------------------------------------------------------
Use Control
Use StandardModel
Use RGEs
Use LoopFunctions
!------------------------------
! parameters of the Lagrangian
!-----------------------------
 !------------------
 ! gauge couplings
 !------------------
 Real(dp), Dimension(3) :: gauge, gauge_mZ, gauge_0, gauge_mR, gauge_mH3
 !----------------------------
 ! Yukawa couplings
 !----------------------------
 Complex(dp), Dimension(3,3) :: Y_l, Y_d, Y_u, Y_l_mZ, Y_d_mZ, Y_u_mZ &
     & , Y_l_0, Y_nu_0, Y_d_0, Y_u_0, Y_nu, Y_T_0, Y_T_MH3, Y_l_mH3   &
     & , Y_d_mH3, Y_u_mH3, Y_T, Y_S_MH3, Y_Z_MH3
 ! there are 3 R-neutrinos, first index is nu_R index
 Complex(dp), Dimension(3,3,3) :: Y_l_mR, Y_nu_mR, Y_d_mR, Y_u_mR
 Complex(dp) :: lam_0, lamp_0, lam12_0(2), lam12_mH3(2), lam12(2)
 !----------------------------
 ! bilinear Higgs parameters
 !----------------------------
 Complex(dp) :: phase_mu, mu, B, mu_mZ, B_mZ, mu_mR, B_mR, mu_0, B_0, B_MH3 &
     & , mu_MH3, M_RS, M_phi, BM_rs, BM_phi, mu_loop, B_loop
 !----------------------------
 ! bilinear R-parity parameters
 !----------------------------
 Complex(dp), Dimension(3) :: eps, Beps
 Real(dp), Dimension(3) :: vevL, lam_ex
 !----------------------------
 ! gaugino mass parameters
 !----------------------------
 Complex(dp), Dimension(3) :: Mi, Mi_mZ, Mi_mR, Mi_0, Mi_mH3
 !----------------------------
 ! trilinear couplings
 !----------------------------
 Complex(dp) :: At_save, Ab_save, Atau_save
 Complex(dp), Dimension(3,3) :: A_l, A_d, A_u, A_l_mZ, A_d_mZ, A_u_mZ, A_l_0   &
     & , A_nu_0, A_d_0, A_u_0, A_nu, A_T_0, A_l_MH3, A_T_MH3, A_d_MH3, A_u_MH3 &
     & , A_T, A_S_MH3, A_Z_MH3
 Complex(dp), Dimension(3,3,3) :: A_l_mR, A_nu_mR, A_d_mR, A_u_mR
 Complex(dp) :: Alam_0, Alamp_0, Alam12_0(2), Alam12_mH3(2), Alam12(2)
 !------------------------------------------------------
 ! trilinear couplings divided by the Yukawa couplings
 !------------------------------------------------------
 Complex(dp), Dimension(3,3) :: AoY_l_0, AoY_nu_0, AoY_d_0, AoY_u_0, AoY_q_0   &
       & , AoY_l, AoY_nu, AoY_d, AoY_u, AoY_l_mZ, AoY_nu_mZ, AoY_d_mZ          &
       & , AoY_u_mZ, AoT_0, AoT_MH3
 Complex(dp) :: Aolam_0, Apolamp_0, Aolam12_0(2), Aolam12_mH3(2)
 !------------------------------------
 ! trilinear couplings, NMSSM
 ! Ao_h0 = A_h0/h0, Ao_lam=A_lam/lam
 ! + spontaneous R-parity violation
 !------------------------------------
 Complex(dp) :: h0, lam, A_h0, A_lam, Ao_h0, Ao_lam, h_pns, A_pns, Ao_hpns
 !-------------------------------------
 ! singlet vev, NMSSM + spontaneous R-parity violation
 !-------------------------------------
 Real(dp) :: vP, vR, vS
 !----------------------------
 ! trilinear R-parity parameters
 !----------------------------
 Complex(dp), Dimension(3,3,3) :: Rp_lam, Rp_lamp, Rp_Alam, Rp_Alamp &
     & , Rp_AoYlam, Rp_AoYlamp
 !----------------------------------
 ! sfermion mass parameters squared
 !----------------------------------
 Complex(dp), Dimension(3,3) :: M2_E, M2_L, M2_D, M2_U, M2_Q, M2_E_mZ, M2_L_mZ &
     & , M2_D_mZ, M2_U_mZ, M2_Q_mZ, M2_E_0, M2_L_0, M2_R_0, M2_D_0, M2_U_0     &
     & , M2_Q_0, M2_R, M2_E_MH3, M2_L_MH3, M2_D_MH3, M2_U_MH3, M2_Q_MH3
 Complex(dp), Dimension(3,3,3) :: M2_E_mR, M2_R_mR, M2_L_mR, M2_D_mR, M2_U_mR &
     & , M2_Q_mR
 !-----------------------------------------------------------
 ! sfermion mass parameters in the super-CKM basis
 !-----------------------------------------------------------
 Complex(dp), Dimension(3,3) :: M2Q_sckm , M2D_sckm , M2U_sckm, Au_sckm, Ad_sckm &
     & , M2Q_0_sckm , M2D_0_sckm , M2U_0_sckm, Au_0_sckm, Ad_0_sckm
 !-----------------------------------------------------------
 ! sfermion mass parameters in the super-PMNS basis
 !-----------------------------------------------------------
 Complex(dp), Dimension(3,3) :: M2L_pmns , M2E_pmns , M2R_pmns, Anu_pmns, Al_pmns &
     & , M2L_0_pmns , M2E_0_pmns , M2R_0_pmns, Anu_0_pmns, Al_0_pmns
 !----------------------------------------------------------------
 ! Higgs mass parameters, tan(beta) and vacuum expectation values
 !----------------------------------------------------------------
 Real(dp) :: M2_H(2), tanb, vevSM(2), M2_H_mZ(2), tanb_mZ, vevSM_mZ(2)      &
     & , M2_H_mR(2), M2_H_0(2), M2_S_0, M2_T_0(2), M2_T_MH3(2), M2_H_MH3(2) &
     & , M_H3(2), M2_T(2), M2_S_MH3(2), M2_Z_MH3(2), MT15_mH3, MZ15_mH3     &
     & , MS15_mH3, M2_P, M2_S
 Logical :: Fifteen_plet = .True., tanb_in_at_Q = .False.
 !--------------------------------------------
 ! neutrino dim. 5 operator + nu mass-matrix
 !--------------------------------------------
 Complex(dp), Dimension(3,3) :: MnuL5, MatNu
 !----------------------------
 ! mass of L- and R-neutrinos 
 !----------------------------
 Real(dp), Dimension(3) :: MnuR
 !-------------------------------------------------------------------
 ! SO(10) models, SO(10) scale and D-term of additional gauge boson
 !-------------------------------------------------------------------
 Real(dp), Save :: M_SO_10=0._dp, D_SO_10=0._dp
 !------------------------------------------------------
 ! Munoz model with 2 NuR, additional parameters
 ! Ao_h02 = A_h02/h02, Ao_lam112=A_lam112/lam112
 ! Ao_lam122=A_lam122/lam122, Ao_lam222=A_lam222/lam2
 !------------------------------------------------------
 Real(dp) :: vP2
 Complex(dp) :: h02, lam2, lam112, lam122, A_h02, A_lam112, A_lam122 &
      & , A_lam222, Ao_h02, Ao_lam112, Ao_lam122, Ao_lam222
!------------------------------
! masses and mixing angles
!------------------------------
! MSSM
!----------------------------- 
 !------------------------------------------------------------------
 ! scalar masses (h,H), masses squared, corresponding mixing matrix, 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------
 Real(dp) :: mS0(2), mS02(2), gT_S0(2), gP_S0(2,200), BR_S0(2,200)
 !------------------------------------------------------------------------
 ! pseudoscalar masses (G,A), masses squared, corresponding mixing matrix 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------------
 Real(dp) :: mP0(2), mP02(2), RP0(2,2), gT_P0(2), gP_P0(2,200), BR_P0(2,200)
 !---------------------------------------------------------------------------
 ! charged scalar masses (G+,H+), masses squared, corresponding mixing matrix 
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSpm(2), mSpm2(2)
 Complex(dp) :: RSpm(2,2)
 Real(dp) :: gT_Spm(2), gP_Spm(2,200), BR_Spm(2,200)
 !---------------------------------------------------------------------------
 ! chargino masses, masses squared, corresponding mixing matrices
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mC(2), mC2(2)
 Complex(dp) :: U(2,2), V(2,2)
 Real(dp) :: gT_C(2), gP_C(2,267), BR_C(2,267), gP_C2(2,120), BR_C2(2,120) &
     & , BR_C3(2,300), gP_C3(2,300)
 !---------------------------------------------------------------------------
 ! neutralino masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) ::  mN(4), mN2(4)
 Complex(dp) :: N(4,4)
 Real(dp) :: gT_N(4), gP_N4_2(4,200), BR_N4_2(4,200) &
           &        , gP_N4_3(4,400), BR_N4_3(4,400)
 Real(dp) :: gT_N5(5), gP_N5(5,350), gP_N(4,350), BR_N(4,350)  &
   & , BR_N5(5,350), gP_N2(5,200), BR_N2(5,200), gP_N3(5,400), BR_N3(5,400)
 !---------------------------------------------------------------------------
 ! gluino mass, phase of the parameter M_3 (=Mi(3))
 ! total decay width, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mGlu
 Complex(dp) :: PhaseGlu
 Real(dp) :: gT_Glu, BR_Glu(230), gP_Glu(230), BR_Glu2(82), gP_Glu2(82) &
      & , BR_Glu3(151), gP_Glu3(151)
 !---------------------------------------------------------------------------
 ! sneutrino masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSneut(3), mSneut2(3)
 Complex(dp) :: Rsneut(3,3)
 Real(dp) :: gT_Sn(3), gP_Sn(3,30), BR_Sn(3,30)
 !---------------------------------------------------------------------------
 ! charged slepton masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSlepton(6), mSlepton2(6)
 Complex(dp) :: RSlepton(6,6)
 Real(dp) :: gT_Sl(6), gP_Sl(6,45), BR_Sl(6,45)
 !---------------------------------------------------------------------------
 ! d-squark masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSdown(6), mSdown2(6)
 Complex(dp) :: RSdown(6,6)
 Real(dp) :: gT_Sd(6), gP_Sd(6,54), BR_Sd(6,54)
 !---------------------------------------------------------------------------
 ! u-squark masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSup(6), mSup2(6)
 Complex(dp) :: RSup(6,6)
 Real(dp) :: gT_Su(6), gP_Su(6,66), BR_Su(6,66), gP_Su3(6,66), BR_Su3(6,66)
!------------------------------
! NMSSM
!-----------------------------
 !------------------------------------------------------------------
 ! scalar masses (h,H,Re(snu)), masses squared, corresponding mixing matrix, 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------
 Real(dp) :: mS03(3), mS032(3), RS03(3,3), gT_S03(3)  &
    & , gP_S03(3,200), BR_S03(3,200)
 !------------------------------------------------------------------------
 ! pseudoscalar masses (G,A,Im(snu)), masses squared, corresponding mixing
 !  matrix, total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------------
 Real(dp) :: mP03(3), mP032(3), RP03(3,3), gT_P03(3)  &
    & , gP_P03(3,200), BR_P03(3,200)
 !---------------------------------------------------------------------------
 ! neutralino masses, masses squared, corresponding mixing matrix
 !---------------------------------------------------------------------------
 Real(dp) ::  mN5(5), mN52(5)
 Complex(dp) :: N5(5,5)
  
!------------------------------
! explicit R-parity violation
!----------------------------- 
 !------------------------------------------------------------------
 ! scalar masses (h,H,Re(snu)), masses squared, corresponding mixing matrix, 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------
 Real(dp) :: mS05(5), mS052(5), RS05(5,5), gT_S05(5)  &
    & , gP_S05(5,200), BR_S05(5,200)
 !------------------------------------------------------------------------
 ! pseudoscalar masses (G,A,Im(snu)), masses squared, corresponding mixing
 !  matrix, total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------------
 Real(dp) :: mP05(5), mP052(5), RP05(5,5), gT_P05(5)  &
    & , gP_P05(5,200), BR_P05(5,200)
 !---------------------------------------------------------------------------
 ! charged scalar masses (G+,H+,sleptons), masses squared, corresponding mixing
 !  matrix, total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSpm8(8), mSpm82(8)
 Complex(dp) :: RSpm8(8,8)
 Real(dp) :: gT_Spm8(8), gP_Spm8(8,200), BR_Spm8(8,200)
 !---------------------------------------------------------------------------
 ! (lepton,chargino) masses, masses squared, corresponding mixing matrices
 !---------------------------------------------------------------------------
 Real(dp) :: mC5(5), mC52(5)
 Complex(dp) :: U5(5,5), V5(5,5)
 !---------------------------------------------------------------------------
 ! neutralino masses, masses squared, corresponding mixing matrix
 !---------------------------------------------------------------------------
 Real(dp) ::  mN7(7), mN72(7)
 Complex(dp) :: N7(7,7)
 !---------------------------------------------------------------------------
 ! mixing matrices for fermions, mainly needed for later extensions
 !---------------------------------------------------------------------------
 Complex(dp), Dimension(3,3) :: uL_L, uL_R, uD_L, uD_R, uU_L, uU_R
 !-----------------------------------------------------------------------
 ! masses of right handed neutrinos
 !-----------------------------------------------------------------------
 Real(dp) :: m_nu_R(3)

!-------------------------------------
! Seesaw II + III a la Florian Staub
!-------------------------------------
 Complex(dp) :: Amue, MassB_H24, MassWB_H24, MassG_H24, mu_H24, mu_h15
 Complex(dp), Dimension(3) ::  mi_h24
 Complex(dp), Dimension(3,3) :: Ye_H24, Yd_H24,Yu_H24, &
  & AYe_H24,AYd_H24,AYu_H24,Ye0_H24, Yd0_H24,Yu0_H24 
 Complex(dp), Dimension(3,3) :: mq2_h24, md2_h24,mu2_h24,ml2_h24,me2_h24
 Complex(dp), Dimension(3) :: Yb0_H24, Yw0_H24, Yx0_H24
 Real(dp) :: g1_H24, g2_H24, g3_H24,gauge_H24(3), Mnu1_H24, Mnu2_H24, &
  & g10_H24, g20_H24, g30_H24,gauge0_H24(3), mHd2_h24, mHu2_h24, mh2_H24(2)
 Complex(dp), Dimension(3,3) :: mHw32, mHg32, mHb32, mHx32, mHxb32, MXM3 &
  & , MWM3, MWM30,MWM3_gut, MGM3, MBM3, AMXM3, AMWM3, AMGM3, AMBM3
 Complex(dp), Dimension(3,3) :: Yb3_H24, Yb3_H24_gut,Yw3_H24, Yx3_H24 &
  & , AYb3_H24, AYw3_H24, AYx3_H24, mi3_h24 
 Complex(dp), Dimension(3,3,3) :: Yb30_H24, Yw30_H24, Yx30_H24
 Complex(dp) :: MBM3running(3,3,3),MGM3running(3,3,3),MWM3running(3,3,3) &
  &  ,MXM3running(3,3,3)
 Real(dp) :: MassMXM3(3), MassMWM3(3), MassMGM3(3), MassMBM3(3)
 ! Seesaw II 
 Complex(dp), Dimension(3,3) :: YT0_H15, YZ0_H15, YS0_15, YT_H15_GUT, AYT_h15  &
  & , AYZ_H15, AYS_H15, Ye_h15, Yd_h15, Yu_h15, Yd0_h15, Yu0_h15, Ye0_h15, AYz &
  & , Ys0_h15, YT_h15, YZ_h15, YS_h15
 Complex(dp) :: Lambda1, Lambda2, MTM, MSM, MZM, AMTM, AMSM, AMZM, ALambda1   &
  & , ALambda2, MTM0 &
  & , Lambda1_gut, Lambda2_gut, MassB_h15, MassG_h15, MassWB_h15, Mi_H15(3)   &
  & , MTM_gut, Lambda10, Lambda20
 Real(dp) ::  gauge_H15(3), g1_h15, g2_h15, g3_h15, g10_h15, g20_h15 &
  & , g30_h15, gauge0_H15(3), mHd2_h15, mHu2_h15, mh2_h15(2),mt2_H15,mtb2_H15 &
  & , mz2_H15, mzb2_H15, ms2_H15, msb2_H15
 Complex(dp), Dimension(3,3) :: mq2_h15, md2_h15,mu2_h15,ml2_h15,me2_h15      &
  & , AYe_h15, AYd_h15, AYu_h15
!-----------------------------------------------------------
! checks if extpar is used for input and collects data
!-----------------------------------------------------------
 Integer :: in_extpar(0:203,2)=0
 Real(dp) :: r_extpar(0:203)=0._dp, i_extpar(203)=0._dp
 Character(len=26) :: n_extpar(0:203)=" "
Contains


 Real(dp) Function CalcScale_from_Stops(M2U, M2Q, Y, T, vevs, mu, g, gp)
 !-----------------------------------------------------------------
 ! using tree-level mass formula to obtain the square of the
 ! renormalisation scale related to the tree-level stop masses
 ! [ CalcScale_from_Stops ] = m^2
 ! written by Werner Porod, 11.09.2014
 !-----------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: Y, T, mu
  Real(dp), Intent(in) :: M2U, M2Q, vevs(2), g, gp
  
  Real(dp) :: diag(2), det, vev2
  Complex(dp) :: offdiag

  Real(dp), Parameter :: T3 = 0.5_dp, Yl = 1._dp / 3._dp  &
   &                     , Yr = -4._dp / 3._dp


  vev2 = 0.25_dp * (vevs(1)**2 - vevs(2)**2)
  diag(1) = M2Q + 0.5_dp * Abs(Y)**2 * vevs(2)**2        &
     &     + (T3 * g**2 - 0.5_dp * Yl * gp**2) * vev2
  diag(2) = M2U + 0.5_dp * Abs(Y)**2 * vevs(2)**2        &
     &     - 0.5_dp * Yr * gp**2 * vev2

  offdiag = (Conjg(T) * vevs(2) - mu * vevs(1) * Conjg(Y)) * oosqrt2

  det = diag(1)*diag(2) - Abs(offdiag)**2

  CalcScale_from_Stops = Sqrt(det)

 End Function CalcScale_from_Stops

 Subroutine GetParameters_at_Q(Q_out, g_i_Q, Y_l_Q, Y_d_Q, Y_u_Q,  Mi_Q, A_l_Q &
              & , A_d_Q, A_u_Q, M2_E_Q, M2_L_Q, M2_D_Q, M2_Q_Q, M2_U_Q, M2_H_Q &
              & , mu_Q, B_Q )
 !---------------------------------------------------------------------
 ! uses RGE running to recalculate the parameters at a different scale
 !---------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: Q_out
  Real(dp), Intent(out) :: g_i_Q(3), M2_H_Q(2)
  Complex(dp), Intent(out) :: Y_l_Q(3,3), Y_d_Q(3,3), Y_u_Q(3,3),  Mi_Q(3)   &
            & , A_l_Q(3,3), A_d_Q(3,3), A_u_Q(3,3), M2_E_Q(3,3), M2_L_Q(3,3) &
            & , M2_D_Q(3,3), M2_Q_Q(3,3), M2_U_Q(3,3), mu_Q, B_Q

  Real(dp) :: g1(213), tz, dt, mudim
  Integer :: kont

  Iname = Iname + 1
  NameOfUnit(Iname) = "GetParameters_at_Q"

  Call ParametersToG(gauge, Y_l, Y_d, Y_u,  Mi, A_l, A_d, A_u, M2_E, M2_L      &
              & , M2_D, M2_Q, M2_U, M2_H, mu, B, g1)

  mudim = Sqrt( GetRenormalizationScale() )

  If (mudim.Ne.Q_out) Then
   tz = Log(Q_out/mudim)
   dt = tz / 50._dp
   kont = 0

   g1(1) = Sqrt(5._dp/3._dp) * g1(1)

   Call odeint(g1, 213, 0._dp, tz, 1.0e-5_dp, dt, 0._dp, rge213, kont)

   g1(1) = Sqrt(3._dp/5._dp) * g1(1)

   If (kont.Ne.0) Then
    Write(ErrCan,*) "Problem in GetParameters_at_Q, RGE running gives",kont
    Call TerminateProgram
   End If
  End If

  Call GToParameters(g1, g_i_Q, Y_l_Q, Y_d_Q, Y_u_Q,  Mi_Q, A_l_Q, A_d_Q    &
       & , A_u_Q, M2_E_Q, M2_L_Q, M2_D_Q, M2_Q_Q, M2_U_Q, M2_H_Q, mu_Q, B_Q )

  Iname = Iname - 1

 End Subroutine GetParameters_at_Q

 Subroutine Set_All_Parameters_0()
 !----------------------------------------------------------------------
 ! this routines sets all model parameters defined above to zero
 !----------------------------------------------------------------------
 Implicit None
  gauge = 0._dp
  gauge_mZ = 0._dp
  gauge_0 = 0._dp
  gauge_mR = 0._dp
  
  Y_l = 0._dp
  Y_d = 0._dp
  Y_u = 0._dp
  Y_l_mZ = 0._dp
  Y_d_mZ = 0._dp
  Y_u_mZ = 0._dp
  Y_l_mR = 0._dp
  Y_nu_mR = 0._dp
  Y_d_mR = 0._dp
  Y_u_mR = 0._dp
  Y_l_0 = 0._dp
  Y_nu_0 = 0._dp
  Y_d_0 = 0._dp
  Y_u_0 = 0._dp
  Y_T_0 = 0._dp
  lam12_0 = 0._dp
  Y_T_MH3 = 0._dp
  Y_Z_MH3 = 0._dp
  Y_S_MH3 = 0._dp
  lam12_MH3 = 0._dp

  phase_mu = 0._dp
  mu = 0._dp
  B = 0._dp
  mu_mZ = 0._dp
  B_mZ = 0._dp
  mu_mR = 0._dp
  B_mR = 0._dp
  mu_0 = 0._dp
  B_0 = 0._dp
  
  eps = 0._dp
  Beps = 0._dp
  vevL = 0._dp
  lam_ex = 0._dp
  
  Rp_lam = 0._dp
  Rp_lamp = 0._dp
  Rp_Alam = 0._dp
  Rp_Alamp = 0._dp
  Rp_AoYlam = 0._dp
  Rp_AoYlamp = 0._dp

  Mi = 0._dp
  Mi_mZ = 0._dp
  Mi_mR = 0._dp
  Mi_0 = 0._dp
  
  A_l = 0._dp
  A_d = 0._dp
  A_u = 0._dp
  A_nu = 0._dp
  A_l_mZ = 0._dp
  A_d_mZ = 0._dp
  A_u_mZ = 0._dp
  A_l_mR = 0._dp
  A_nu_mR = 0._dp
  A_d_mR = 0._dp
  A_u_mR = 0._dp
  A_l_0 = 0._dp
  A_nu_0 = 0._dp
  A_d_0 = 0._dp
  A_u_0 = 0._dp
  AoY_l_0 = 0._dp
  AoY_nu_0 = 0._dp
  AoY_d_0 = 0._dp
  AoY_u_0 = 0._dp
  AoY_q_0 = 0._dp
  AoY_l = 0._dp
  AoY_nu = 0._dp
  AoY_d = 0._dp
  AoY_u = 0._dp
  AoY_l_mZ = 0._dp
  AoY_nu_mZ = 0._dp
  AoY_d_mZ = 0._dp
  AoY_u_mZ = 0._dp
  A_T_0 = 0._dp
  Alam12_0 = 0._dp
  A_T_MH3 = 0._dp
  A_S_MH3 = 0._dp
  A_Z_MH3 = 0._dp
  Alam12_MH3 = 0._dp
  AoT_0 = 0._dp
  Aolam12_0 = 0._dp
  AoT_MH3 = 0._dp
  Aolam12_MH3 = 0._dp
  At_save = 0._dp
  Ab_save = 0._dp
  Atau_save = 0._dp

  Ad_sckm = ZeroC
  Au_sckm = ZeroC
  Al_pmns = ZeroC
  Ad_0_sckm = ZeroC
  Au_0_sckm = ZeroC
  Al_0_pmns = ZeroC

  h0 = 0._dp
  lam = 0._dp
  A_h0 = 0._dp
  A_lam = 0._dp
  Ao_h0 = 0._dp
  Ao_lam = 0._dp

  vP = 0._dp
  
  h02 = 0._dp
  lam2 = 0._dp
  lam112 = 0._dp
  lam122 = 0._dp
  A_h02 = 0._dp
  A_lam112 = 0._dp
  A_lam122 = 0._dp
  A_lam222 = 0._dp
  Ao_h02 = 0._dp
  Ao_lam112 = 0._dp
  Ao_lam122 = 0._dp
  Ao_lam222 = 0._dp

  vP2 = 0._dp

  M2_E = 0._dp
  M2_L = 0._dp
  M2_R = 0._dp
  M2_D = 0._dp
  M2_U = 0._dp
  M2_Q = 0._dp
  M2_E_mZ = 0._dp
  M2_L_mZ = 0._dp
  M2_D_mZ = 0._dp
  M2_U_mZ = 0._dp
  M2_Q_mZ = 0._dp
  M2_E_mR = 0._dp
  M2_R_mR = 0._dp
  M2_L_mR = 0._dp
  M2_D_mR = 0._dp
  M2_U_mR = 0._dp
  M2_Q_mR = 0._dp
  M2_E_0 = 0._dp
  M2_L_0 = 0._dp
  M2_R_0 = 0._dp
  M2_D_0 = 0._dp
  M2_U_0 = 0._dp
  M2_Q_0 = 0._dp

  M2Q_sckm = 0._dp
  M2D_sckm = 0._dp
  M2U_sckm = 0._dp
  M2Q_0_sckm = 0._dp
  M2D_0_sckm = 0._dp
  M2U_0_sckm = 0._dp

  M2E_pmns = 0._dp
  M2L_pmns = 0._dp
  M2E_0_pmns = 0._dp
  M2L_0_pmns = 0._dp

  M2_S_MH3 = 0._dp
  M2_Z_MH3 = 0._dp
  MT15_mH3 = 0._dp
  MZ15_mH3 = 0._dp
  MS15_mH3= 0._dp

  M2_H = 0._dp
  tanb = 0._dp
  vevSM = 0._dp
  M2_H_mZ = 0._dp
  tanb_mZ = 0._dp
  vevSM_mZ = 0._dp
  M2_H_mR = 0._dp
  M2_H_0 = 0._dp
  
  MnuL5 = 0._dp
  MnuR = 0._dp

  lam_0 = 0._dp
  Alam_0 = 0._dp
  Aolam_0 = 0._dp
  lamp_0 = 0._dp
  Alamp_0 = 0._dp
  Apolamp_0 = 0._dp

  M_SO_10=0._dp
  D_SO_10=0._dp

  U = 0._dp
  U5 = 0._dp
  V = 0._dp
  V5 = 0._dp
  N = 0._dp
  N5 = 0._dp
  N7 = 0._dp
  RSup = 0._dp
  RSdown = 0._dp
  RSneut = 0._dp
  RSlepton = 0._dp
  
 !------------------------------------------------
 ! Seesaw II + III a la Florian Staub
 !------------------------------------------------
  Amue = 0._dp
  MassB_H24 = 0._dp
  MassWB_H24 = 0._dp
  MassG_H24 = 0._dp
  mu_H24 = 0._dp
  mu_h15 = 0._dp
  mi_h24 = 0._dp
  Ye_H24 = 0._dp
  Yd_H24 = 0._dp
  Yu_H24 = 0._dp 
  AYe_H24 = 0._dp
  AYd_H24 = 0._dp
  AYu_H24 = 0._dp
  Ye0_H24 = 0._dp
  Yd0_H24 = 0._dp
  Yu0_H24 = 0._dp
  mq2_h24 = 0._dp
  md2_h24 = 0._dp
  mu2_h24 = 0._dp
  ml2_h24 = 0._dp
  me2_h24 = 0._dp
  Yb0_H24 = 0._dp
  Yw0_H24 = 0._dp
  Yx0_H24 = 0._dp
  g1_H24 = 0._dp
  g2_H24 = 0._dp
  g3_H24 = 0._dp
  gauge_H24 = 0._dp
  Mnu1_H24 = 0._dp
  Mnu2_H24 = 0._dp
  g10_H24 = 0._dp
  g20_H24 = 0._dp
  g30_H24 = 0._dp
  gauge0_H24 = 0._dp
  mHd2_h24 = 0._dp
  mHu2_h24 = 0._dp
  mh2_H24 = 0._dp
  mHw32 = 0._dp 
  mHg32 = 0._dp
  mHb32 = 0._dp
  mHx32 = 0._dp
  mHxb32 = 0._dp
  MXM3 = 0._dp
  MWM3 = 0._dp
  MWM30 = 0._dp
  MWM3_gut = 0._dp
  MGM3 = 0._dp
  MBM3 = 0._dp
  AMXM3 = 0._dp
  AMWM3 = 0._dp
  AMGM3 = 0._dp
  AMBM3 = 0._dp
  Yb3_H24 = 0._dp
  Yb3_H24_gut = 0._dp
  Yw3_H24 = 0._dp
  Yx3_H24 = 0._dp
  AYb3_H24 = 0._dp
  AYw3_H24 = 0._dp
  AYx3_H24 = 0._dp
  mi3_h24 = 0._dp  
  Yb30_H24 = 0._dp
  Yw30_H24 = 0._dp
  Yx30_H24 = 0._dp
  YT0_H15 = 0._dp
  YZ0_H15 = 0._dp
  YS0_15 = 0._dp
  YT_H15_GUT = 0._dp
  AYT_h15 = 0._dp
  AYZ_H15 = 0._dp
  AYS_H15 = 0._dp
  Ye_h15 = 0._dp
  Yd_h15 = 0._dp
  Yu_h15 = 0._dp
  Yd0_h15 = 0._dp
  Yu0_h15 = 0._dp
  Ye0_h15 = 0._dp
  AYz = 0._dp
  Ys0_h15 = 0._dp
  YT_h15 = 0._dp
  YZ_h15 = 0._dp 
  YS_h15 = 0._dp
  Lambda1 = 0._dp
  Lambda2 = 0._dp
  MTM = 0._dp 
  MSM = 0._dp
  MZM = 0._dp
  AMTM = 0._dp
  AMSM = 0._dp
  AMZM = 0._dp
  ALambda1 = 0._dp
  ALambda2 = 0._dp
  MTM0 = 0._dp
  mt2_H15 = 0._dp
  mtb2_H15 = 0._dp
  mz2_H15 = 0._dp
  mzb2_H15 = 0._dp
  ms2_H15 = 0._dp
  msb2_H15 = 0._dp
  Lambda1_gut = 0._dp
  Lambda2_gut = 0._dp
  MassB_h15 = 0._dp
  MassG_h15 = 0._dp
  MassWB_h15 = 0._dp
  Mi_H15 = 0._dp
  MTM_gut = 0._dp
  Lambda10 = 0._dp
  Lambda20 = 0._dp
  gauge_H15 = 0._dp
  g1_h15 = 0._dp
  g2_h15 = 0._dp
  g3_h15 = 0._dp
  g10_h15 = 0._dp
  g20_h15 = 0._dp
  g30_h15 = 0._dp
  gauge0_H15 = 0._dp
  mHd2_h15 = 0._dp
  mHu2_h15 = 0._dp
  mh2_h15 = 0._dp
  mq2_h15 = 0._dp
  md2_h15 = 0._dp
  mu2_h15 = 0._dp
  ml2_h15 = 0._dp
  me2_h15 = 0._dp  
  AYe_h15 = 0._dp
  AYd_h15 = 0._dp
  AYu_h15 = 0._dp

 End Subroutine Set_All_Parameters_0


 Real(dp) Function product_stop_masses(mSup, RSup, GenerationMixing)
 !--------------------------------------------------------------
 ! is needed to determine the renormalisation scale
 ! written by Werner Porod, 09.12.2013
 !--------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mSup(6)
  Complex(dp), Intent(in) :: RSup(6,6)
  Logical, Intent(in) :: GenerationMixing

  Integer :: i1, count
  Real(dp) :: v(6), max_stop

  Iname = Iname + 1
  NameOfUnit(Iname) = "product_stop_masses"

  If (GenerationMixing) Then
   v = Abs(RSup(:,3))**2 + Abs(RSup(:,6))**2
   max_stop = Maxval(v)
   count = 0
   product_stop_masses = 1._dp
 
   Do i1=1,6
    If (max_stop.Eq.v(i1)) Then
     product_stop_masses = mSup(i1)
     v(i1) = 0._dp
     count = count + 1
     Exit
    End If
   End Do

   max_stop = Maxval(v)

   Do i1=1,6
    If (max_stop.Eq.v(i1)) Then
     product_stop_masses = product_stop_masses * mSup(i1)
     v(i1) = 0._dp
     count = count + 1
     Exit
    End If
   End Do

   If (count.Ne.2) Then
    Write(ErrCan,*) "Problem in routine "//Trim(NameOfUnit(Iname))
    Write(ErrCan,*) "count =", count
    Write(ErrCan,*) "Setting scale to 1 TeV"
    product_stop_masses = 1.6_dp
    If (ErrorLevel.Eq.2) Call TerminateProgram
   End If

  Else ! .not.GenerationMixing
   product_stop_masses = mSup(5)*mSup(6)
  End If

  Iname = Iname - 1

 End Function product_stop_masses

 Subroutine Switch_from_superCKM(Y_d, Y_u, Ad_in, Au_in, MD_in, MQ_in, MU_in &
                      &, Ad_out, Au_out, MD_out, MQ_out, MU_out, tr        &
                      &, RSd_in, RSu_in, RSd_out, RSu_out, CKM_out, Yd, Yu )
 !---------------------------------------------------------------------------
 ! shifts the parameter from the  super CKM basis to the electroweak basis
 ! written by werner Porod, 12.03.08
 !---------------------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in), Dimension(3,3) :: Y_d, Y_u, Au_in, Ad_in, MD_in &
        & , MQ_in, MU_in
  Complex(dp), Optional, Intent(in), Dimension(6,6) :: RSu_in, RSd_in
  Logical, Intent(in) :: tr  ! if true, then the matrices are transposed 
                             ! compared to low energy definition
  Complex(dp), Intent(out), Dimension(3,3) :: Au_out, Ad_out, MD_out, MQ_out &
        & , MU_out
  Complex(dp), Optional, Intent(out), Dimension(6,6) :: RSu_out, RSd_out
  Complex(dp), Optional, Intent(out) :: CKM_out(3,3)
  Real(dp), Optional, Intent(out) :: Yd(3), Yu(3)

  Complex(dp), Dimension(3,3) :: uU_L, uU_R, uD_L, uD_R, CKM_Q
  Complex(dp) :: rot(6,6), Ephi

  Real(dp) :: mf(3), s12, s23, aR, aI, s13, c13
  Integer :: ierr

  !------------------------------------------
  ! diagonalizing d- and u-Yukawa couplings
  ! I am only interested in the mixing matrices
  !------------------------------------------
  If (tr) Then
   Call FermionMass(Transpose(Y_u), 1._dp, mf, uU_L, uU_R, ierr)
   If (Present(Yu)) Yu = sqrt2 * mf
   Call FermionMass(Transpose(Y_d), 1._dp, mf, uD_L, uD_R, ierr)
   If (Present(Yd)) Yd = sqrt2 * mf
  Else
   Call FermionMass(Y_u, 1._dp, mf, uU_L, uU_R, ierr)
   If (Present(Yu)) Yu = sqrt2 * mf
   Call FermionMass(Y_d, 1._dp, mf, uD_L, uD_R, ierr)
   If (Present(Yd)) Yd = sqrt2 * mf
  End If
  !---------------------------------------------------------
  ! CKM matrix at Q, shifting phases according to PDG form
  !---------------------------------------------------------
  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )
  uD_L(1,:) = uD_L(1,:) / Conjg(CKM_Q(1,1)) * Abs(CKM_Q(1,1))
  uD_L(2,:) = uD_L(2,:) / Conjg(CKM_Q(1,2)) * Abs(CKM_Q(1,2))
  uU_L(2,:) = uU_L(2,:) / CKM_Q(2,3) * Abs(CKM_Q(2,3))
  uU_L(3,:) = uU_L(3,:) / CKM_Q(3,3) * Abs(CKM_Q(3,3))
  !-------------------------------------------------------------------
  ! also the right quark must be multiplied with the conjugate phase
  ! as otherwise the masses get complex
  !-------------------------------------------------------------------
  uD_R(1,:) = uD_R(1,:) / CKM_Q(1,1) * Abs(CKM_Q(1,1))
  uD_R(2,:) = uD_R(2,:) / CKM_Q(1,2) * Abs(CKM_Q(1,2))
  uU_R(2,:) = uU_R(2,:) / Conjg(CKM_Q(2,3)) * Abs(CKM_Q(2,3))
  uU_R(3,:) = uU_R(3,:) / Conjg(CKM_Q(3,3)) * Abs(CKM_Q(3,3))
  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )

  !--------------------------------------------------------------
  ! one more freedom left
  !--------------------------------------------------------------
  s13 = Abs(CKM_Q(1,3))
  c13 = sqrt(1._dp - s13**2)
  s23 = Abs(CKM_Q(2,3))/c13
  s12 = Abs(CKM_Q(1,2))/c13

  aR = Real(CKM_Q(2,2),dp) + s12 * s23 * Real(CKM_Q(1,3),dp)
  aI =  s12 * s23 * Aimag(CKM_Q(1,3)) - Aimag(CKM_Q(2,2))
  Ephi = Cmplx(aR/Sqrt(aR**2+aI**2),aI/Sqrt(aR**2+aI**2),dp)

  uU_L(2:3,:) = Ephi * uU_L(2:3,:)
  uD_L(3,:) = Ephi * uD_L(3,:)
  Ephi = Conjg(Ephi)
  uU_R(2:3,:) = Ephi * uU_R(2:3,:)
  uD_R(3,:) = Ephi * uD_R(3,:)

  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )

  If (Present(CKM_out)) CKM_out = CKM_Q
  !-------------------------------------------------------------------
  ! shifting the parameters from the super CKM basis
  !-------------------------------------------------------------------
  If (tr) Then
   Au_out = Matmul( Matmul(Transpose(uU_R), Au_in), uU_L)
   Ad_out = Matmul( Matmul(Transpose(uD_R), Ad_in), uD_L)

   MD_out = Matmul( Matmul( Transpose(Conjg(uD_R)), MD_in), uD_R)
   MU_out = Matmul( Matmul( Transpose(Conjg(uU_R)), MU_in), uU_R)
   MQ_out = Matmul( Matmul( Transpose(uD_L), MQ_in), Conjg(uD_L) )

  Else
   Au_out = Matmul( Matmul(Transpose(uU_L), Au_in), uU_R)
   Ad_out = Matmul( Matmul(Transpose(uD_L), Ad_in), uD_R)

   MD_out = Matmul( Matmul( Transpose(uD_R), MD_in), Conjg(uD_R))
   MU_out = Matmul( Matmul( Transpose(uU_R), MU_in), Conjg(uU_R))
   MQ_out = Matmul( Matmul( Transpose(Conjg(uD_L)), MQ_in), uD_L )

  End If
  !------------------------------------------------------------------
  ! to avoid numerical problems ensure that matrices are hermitian
  !-----------------------------------------------------------------
  MD_out = 0.5_dp * ( MD_out + Conjg(Transpose(MD_out)) )
  MU_out = 0.5_dp * ( MU_out + Conjg(Transpose(MU_out)) )
  MQ_out = 0.5_dp * ( MQ_out + Conjg(Transpose(MQ_out)) )

   If (Present(RSu_in).And.Present(RSu_out)) Then
    rot = 0._dp
    rot(1:3,1:3) = Conjg(uU_L)
    rot(4:6,4:6) = uU_R
    RSu_out = Matmul(RSu_in, rot)
   End If
   If (Present(RSd_in).And.Present(RSd_out)) Then
    rot = 0._dp
    rot(1:3,1:3) = Conjg(uD_L)
    rot(4:6,4:6) = uD_R
    RSd_out = Matmul(RSd_in, rot)
   End If

 End Subroutine Switch_from_superCKM

 Subroutine Switch_to_superCKM(Y_d, Y_u, Ad_in, Au_in, MD_in, MQ_in, MU_in &
                      &, Ad_out, Au_out, MD_out, MQ_out, MU_out, tr        &
                      &, RSd_in, RSu_in, RSd_out, RSu_out, CKM_out, Yd, Yu )
 !---------------------------------------------------------------------------
 ! shifts the parameter from the electroweak basis to the super CKM basis
 ! written by werner Porod, 12.03.08
 !---------------------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in), Dimension(3,3) :: Y_d, Y_u, Au_in, Ad_in, MD_in &
        & , MQ_in, MU_in
  Complex(dp), Optional, Intent(in), Dimension(6,6) :: RSu_in, RSd_in
  Logical, Intent(in) :: tr  ! if true, then the matrices are transposed 
                             ! compared to low energy definition
  Complex(dp), Intent(out), Dimension(3,3) :: Au_out, Ad_out, MD_out, MQ_out &
        & , MU_out
  Complex(dp), Optional, Intent(out), Dimension(6,6) :: RSu_out, RSd_out
  Complex(dp), Optional, Intent(out) :: CKM_out(3,3)
  Real(dp), Optional, Intent(out) :: Yd(3), Yu(3)

  Complex(dp), Dimension(3,3) :: uU_L, uU_R, uD_L, uD_R, CKM_Q
  Complex(dp) :: rot(6,6), Ephi

  Real(dp) :: mf(3), s12, s23, aR, aI, s13, c13
  Integer :: ierr

  !------------------------------------------
  ! diagonalizing d- and u-Yukawa couplings
  ! I am only interested in the mixing matrices
  !------------------------------------------
  If (tr) Then
   Call FermionMass(Transpose(Y_u), 1._dp, mf, uU_L, uU_R, ierr)
   If (Present(Yu)) Yu = sqrt2 * mf
   Call FermionMass(Transpose(Y_d), 1._dp, mf, uD_L, uD_R, ierr)
   If (Present(Yd)) Yd = sqrt2 * mf
  Else
   Call FermionMass(Y_u, 1._dp, mf, uU_L, uU_R, ierr)
   If (Present(Yu)) Yu = sqrt2 * mf
   Call FermionMass(Y_d, 1._dp, mf, uD_L, uD_R, ierr)
   If (Present(Yd)) Yd = sqrt2 * mf
  End If
  !---------------------------------------------------------
  ! CKM matrix at Q, shifting phases according to PDG form
  !---------------------------------------------------------
  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )
  uD_L(1,:) = uD_L(1,:) / Conjg(CKM_Q(1,1)) * Abs(CKM_Q(1,1))
  uD_L(2,:) = uD_L(2,:) / Conjg(CKM_Q(1,2)) * Abs(CKM_Q(1,2))
  uU_L(2,:) = uU_L(2,:) / CKM_Q(2,3) * Abs(CKM_Q(2,3))
  uU_L(3,:) = uU_L(3,:) / CKM_Q(3,3) * Abs(CKM_Q(3,3))
  !-------------------------------------------------------------------
  ! also the right quark must be multiplied with the conjugate phase
  ! as otherwise the masses get complex
  !-------------------------------------------------------------------
  uD_R(1,:) = uD_R(1,:) / CKM_Q(1,1) * Abs(CKM_Q(1,1))
  uD_R(2,:) = uD_R(2,:) / CKM_Q(1,2) * Abs(CKM_Q(1,2))
  uU_R(2,:) = uU_R(2,:) * Abs(CKM_Q(2,3)) / Conjg(CKM_Q(2,3))
  uU_R(3,:) = uU_R(3,:) * Abs(CKM_Q(3,3)) / Conjg(CKM_Q(3,3))
  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )

  !--------------------------------------------------------------
  ! one more freedom left
  !--------------------------------------------------------------
  s13 = Abs(CKM_Q(1,3))
  c13 = sqrt(1._dp - s13**2)
  s23 = Abs(CKM_Q(2,3))/c13
  s12 = Abs(CKM_Q(1,2))/c13

  aR = Real(CKM_Q(2,2),dp) + s12 * s23 * Real(CKM_Q(1,3),dp)
  aI =  s12 * s23 * Aimag(CKM_Q(1,3)) - Aimag(CKM_Q(2,2))
  Ephi = Cmplx(aR/Sqrt(aR**2+aI**2),aI/Sqrt(aR**2+aI**2),dp)

  uU_L(2:3,:) = Ephi * uU_L(2:3,:)
  uD_L(3,:) = Ephi * uD_L(3,:)
  Ephi = Conjg(Ephi)
  uU_R(2:3,:) = Ephi * uU_R(2:3,:)
  uD_R(3,:) = Ephi * uD_R(3,:)

  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )

  If (Present(CKM_out)) CKM_out = CKM_Q
  !-------------------------------------------------------------------
  ! shifting the parameters to the super CKM basis
  !-------------------------------------------------------------------
  If (tr) Then
   Au_out = Matmul( Matmul(Conjg(uU_R), Au_in), Transpose(Conjg(uU_L)))
   Ad_out = Matmul( Matmul(Conjg(uD_R), Ad_in), Transpose(Conjg(uD_L)))

   MD_out = Matmul( Matmul( uD_R, MD_in), Transpose(Conjg(uD_R)))
   MU_out = Matmul( Matmul( uU_R, MU_in), Transpose(Conjg(uU_R)))
   MQ_out = Matmul( Matmul( Conjg(uD_L), MQ_in), Transpose(uD_L) )
  Else
   Au_out = Matmul( Matmul(Conjg(uU_L), Au_in), Transpose(Conjg(uU_R)))
   Ad_out = Matmul( Matmul(Conjg(uD_L), Ad_in), Transpose(Conjg(uD_R)))

   MD_out = Matmul( Matmul( Conjg(uD_R), MD_in), Transpose(uD_R))
   MU_out = Matmul( Matmul( Conjg(uU_R), MU_in), Transpose(uU_R))
   MQ_out = Matmul( Matmul( uD_L, MQ_in), Transpose(Conjg(uD_L)) )

  End If
  !------------------------------------------------------------------
  ! to avoid numerical problems ensure that matrices are hermitian
  !-----------------------------------------------------------------
  MD_out = 0.5_dp * ( MD_out + Conjg(Transpose(MD_out)) )
  MU_out = 0.5_dp * ( MU_out + Conjg(Transpose(MU_out)) )
  MQ_out = 0.5_dp * ( MQ_out + Conjg(Transpose(MQ_out)) )

   If (Present(RSu_in).And.Present(RSu_out)) Then
    rot = 0._dp
    rot(1:3,1:3) = Transpose(uU_L)
    rot(4:6,4:6) = Transpose(Conjg(uU_R))
    RSu_out = Matmul(RSu_in, rot)
   End If
   If (Present(RSd_in).And.Present(RSd_out)) Then
    rot = 0._dp
    rot(1:3,1:3) = Transpose(uD_L)
    rot(4:6,4:6) = Transpose(Conjg(uD_R))
    RSd_out = Matmul(RSd_in, rot)
   End If

 End Subroutine Switch_to_superCKM

 Subroutine Switch_from_superPMNS(Y_l, M5_nu, Al_in, ME_in, ML_in &
                      &, Al_out, ME_out, ML_out, tr            &
                      &, RSl_in, RSn_in, RSl_out, RSn_out, PMNS_out, Yl )
 !---------------------------------------------------------------------------
 ! shifts the parameter from the electroweak basis to the super PMNS basis
 ! written by werner Porod, 12.03.08
 !---------------------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in), Dimension(3,3) :: Y_l, M5_nu, Al_in, ME_in, ML_in
  Complex(dp), Optional, Intent(in) :: RSl_in(6,6), RSn_in(3,3)
  Logical, Intent(in) :: tr  ! if true, then the matrices are transposed 
                             ! compared to low energy definition
  Complex(dp), Intent(out), Dimension(3,3) :: Al_out, ME_out, ML_out
  Complex(dp), Optional, Intent(out) :: RSl_out(6,6)
  Complex(dp), Optional, Intent(out) :: PMNS_out(3,3), RSn_out(3,3)
  Real(dp), Optional, Intent(out) :: Yl(3)

  Complex(dp), Dimension(3,3) :: uN_L, uL_L, uL_R, PMNS_Q
  Complex(dp) :: rot(6,6)

  Real(dp) :: mf(3)
  Integer :: ierr

  !------------------------------------------
  ! diagonalizing d- and u-Yukawa couplings
  ! I am only interested in the mixing matrices
  !------------------------------------------
  If (tr) Then
   Call FermionMass(Transpose(Y_l), 1._dp, mf, uL_L, uL_R, ierr)
   If (Present(Yl)) Yl = sqrt2 * mf
  Else
   Call FermionMass(Y_l, 1._dp, mf, uL_L, uL_R, ierr)
   If (Present(Yl)) Yl = sqrt2 * mf
  End If

  Call Neutrinomasses(M5_nu, mf, uN_L, ierr)
  !---------------------------------------------------------
  ! PMNS matrix at Q, shifting phases according to PDG form
  !---------------------------------------------------------
  PMNS_Q =  Matmul(uL_L, Transpose(Conjg(uN_L)))

  If (Present(PMNS_out)) PMNS_out = PMNS_Q
  !-------------------------------------------------------------------
  ! shifting the parameters to the super PMNS basis
  !-------------------------------------------------------------------
  If (tr) Then
   Al_out = Matmul( Matmul( Transpose(uL_R), Al_in), uL_L)

   ME_out = Matmul( Matmul( Transpose(Conjg(uL_R)), ME_in), uL_R)
   ML_out = Matmul( Matmul( Transpose(uL_L), ML_in), Conjg(uL_L) )

  Else
   Al_out = Matmul( Matmul( Transpose(uL_L), Al_in), uL_R)

   ME_out = Matmul( Matmul( Transpose(uL_R), ME_in), Conjg(uL_R))
   ML_out = Matmul( Matmul( Transpose(Conjg(uL_L)), ML_in), uL_L)

  End If
  !------------------------------------------------------------------
  ! to avoid numerical problems ensure that matrices are hermitian
  !-----------------------------------------------------------------
  ME_out = 0.5_dp * ( ME_out + Conjg(Transpose(ME_out)) )
  ML_out = 0.5_dp * ( ML_out + Conjg(Transpose(ML_out)) )

   If (Present(RSn_in).And.Present(RSn_out)) Then
    RSn_out = Matmul(RSn_in, Conjg(uN_L) )
   End If
   If (Present(RSl_in).And.Present(RSl_out)) Then
    rot = 0._dp
    rot(1:3,1:3) = Conjg(uL_L)
    rot(4:6,4:6) = uL_R
    RSl_out = Matmul(RSl_in, rot)
   End If

 End Subroutine Switch_from_superPMNS

 Subroutine Switch_to_superPMNS(Y_l, M5_nu, Al_in, ME_in, ML_in &
                      &, Al_out, ME_out, ML_out, tr            &
                      &, RSl_in, RSn_in, RSl_out, RSn_out, PMNS_out, Yl )
 !---------------------------------------------------------------------------
 ! shifts the parameter from the electroweak basis to the super PMNS basis
 ! written by werner Porod, 12.03.08
 !---------------------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in), Dimension(3,3) :: Y_l, M5_nu, Al_in, ME_in, ML_in
  Complex(dp), Optional, Intent(in) :: RSl_in(6,6), RSn_in(3,3)
  Logical, Intent(in) :: tr  ! if true, then the matrices are transposed 
                             ! compared to low energy definition
  Complex(dp), Intent(out), Dimension(3,3) :: Al_out, ME_out, ML_out
  Complex(dp), Optional, Intent(out) :: RSl_out(6,6)
  Complex(dp), Optional, Intent(out) :: PMNS_out(3,3), RSn_out(3,3)
  Real(dp), Optional, Intent(out) :: Yl(3)

  Complex(dp), Dimension(3,3) :: uN_L, uL_L, uL_R, PMNS_Q
  Complex(dp) :: rot(6,6)

  Real(dp) :: mf(3)
  Integer :: ierr

  !------------------------------------------
  ! diagonalizing d- and u-Yukawa couplings
  ! I am only interested in the mixing matrices
  !------------------------------------------
  If (tr) Then
   Call FermionMass(Transpose(Y_l), 1._dp, mf, uL_L, uL_R, ierr)
   If (Present(Yl)) Yl = sqrt2 * mf
  Else
   Call FermionMass(Y_l, 1._dp, mf, uL_L, uL_R, ierr)
   If (Present(Yl)) Yl = sqrt2 * mf
  End If

  Call Neutrinomasses(M5_nu, mf, uN_L, ierr)
  !---------------------------------------------------------
  ! PMNS matrix at Q, shifting phases according to PDG form
  !---------------------------------------------------------
  PMNS_Q =  Matmul(uL_L, Transpose(Conjg(uN_L)))

  If (Present(PMNS_out)) PMNS_out = PMNS_Q
  !-------------------------------------------------------------------
  ! shifting the parameters to the super PMNS basis
  !-------------------------------------------------------------------
  If (tr) Then
   Al_out = Matmul( Matmul(Conjg(uL_R), Al_in), Transpose(Conjg(uL_L)))

   ME_out = Matmul( Matmul( uL_R, ME_in), Transpose(Conjg(uL_R)))
   ML_out = Matmul( Matmul( Conjg(uL_L), ML_in), Transpose(uL_L) )

  Else
   Al_out = Matmul( Matmul(Conjg(uL_L), Al_in), Transpose(Conjg(uL_R)))

   ME_out = Matmul( Matmul( Conjg(uL_R), ME_in), Transpose(uL_R))
   ML_out = Matmul( Matmul( uL_L, ML_in), Transpose(Conjg(uL_L)) )

  End If
  !------------------------------------------------------------------
  ! to avoid numerical problems ensure that matrices are hermitian
  !-----------------------------------------------------------------
  ME_out = 0.5_dp * ( ME_out + Conjg(Transpose(ME_out)) )
  ML_out = 0.5_dp * ( ML_out + Conjg(Transpose(ML_out)) )

   If (Present(RSn_in).And.Present(RSn_out)) Then
    RSn_out = Matmul(RSn_in, Transpose(uN_L) )
   End If
   If (Present(RSl_in).And.Present(RSl_out)) Then
    rot = 0._dp
    rot(1:3,1:3) = Transpose(uL_L)
    rot(4:6,4:6) = Transpose(Conjg(uL_R))
    RSl_out = Matmul(RSl_in, rot)
   End If

 End Subroutine Switch_to_superPMNS

End Module Model_Data

Module MSSM_data
!-----------------------------------------------------------
! in this modul all basic information concerning the MSSM are 
! stored
! by Werner Porod, 15.09.2010
!-----------------------------------------------------------
Use Control

 !-------------------------------------------------------------
 ! number of different particle species, can be changed using
 ! the initialisation routine for other SUSY models
 !-------------------------------------------------------------
 Integer :: n_nu=3, n_l=3, n_d=3, n_u=3 ! generation of SM fermions
 Integer :: n_Z=1, n_W=1                ! number of vector bosons
 Integer :: n_snu=3, n_sle=6, n_sd=6, n_su = 6 ! number of sfermions
 Integer :: n_n=4, n_c=2, n_g=1         ! number of neutralinos, charginos
                                        ! and gluinos
 Integer :: n_s0=2, n_P0=2, n_Spm=2     ! number of scalar Higgs bosons, 
                                        ! pseudoscalars and charged scalars
                                        ! including Goldstone bosons
 !------------------------------------------------------------------
 ! scalar masses (h,H), masses squared, corresponding mixing matrix, 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------
 Real(dp) :: RS0(2,2)
 Type(particle23) :: S0(2)
 Real(dp) :: r_GlGlS0(2) ! ratio of h-g-g over SM values squared
 !------------------------------------------------------------------------
 ! pseudoscalar masses (G,A), masses squared, corresponding mixing matrix 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------------
! Real(dp) :: RP0(2,2)
 Type(particle2) :: P0(2)
 Real(dp) :: r_GlGlP0 ! ratio of A-g-g over SM values squared
 !---------------------------------------------------------------------------
 ! charged scalar masses (G+,H+), masses squared, corresponding mixing matrix 
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: RSpm(2,2)
 Type(particle2) :: Spm(2)
 !---------------------------------------------------------------------------
 ! chargino masses, masses squared, corresponding mixing matrices
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: U(2,2), V(2,2)
 Type(particle23) :: ChiPM(2)
 !---------------------------------------------------------------------------
 ! neutralino masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: N(4,4)
 Type(particle23) :: Chi0(4)
 !---------------------------------------------------------------------------
 ! gluino mass, phase of the parameter M_3 (=Mi(3))
 ! total decay width, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: PhaseGlu
 Type(particle23) :: Glu
 !---------------------------------------------------------------------------
 ! gluino mass, phase of the parameter M_3 (=Mi(3))
 ! total decay width, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: PhaseGlu
 Type(particle2) :: Grav
 !---------------------------------------------------------------------------
 ! sneutrino masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: Rsneut(3,3)
 Type(particle23) :: Sneut(3)
 !---------------------------------------------------------------------------
 ! charged slepton masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: RSlepton(6,6)
 Type(particle23) :: Slepton(6)
 !---------------------------------------------------------------------------
 ! d-squark masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: RSdown(6,6)
 Type(particle2) :: Sdown(6)
 !---------------------------------------------------------------------------
 ! u-squark masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: RSup(6,6)
 Type(particle23) :: Sup(6)
 !------------
 ! PDG codes
 !------------
 Integer, Parameter ::id_A0 = 36, id_Hp = 37, id_h0(2) = (/25, 35 /)          &
      & , id_Wboson = 24, id_Zboson = 23, id_photon = 22, id_gluon = 21       &
      & , id_lept(3) = (/ 11, 13, 15 /), id_neut(3) = (/ 12, 14, 16 /)        &
      & , id_u_quark(3) = (/ 2, 4, 6 /), id_d_quark(3) = (/ 1, 3, 5 /)        &
      & , id_gravitino = 1000039, id_gluino = 1000021                         &
      & , id_sneut(3) = (/ 1000012, 1000014, 1000016 /)                       &
      & , id_slept(6) = (/1000011,2000011, 1000013,2000013, 1000015,2000015/) &
      & , id_u_sq(6) = (/1000002,2000002, 1000004,2000004, 1000006,2000006/)  &
      & , id_d_sq(6) = (/1000001,2000001, 1000003,2000003, 1000005,2000005/)  &
      & , id_neutralino(4) = (/1000022, 1000023, 1000025, 1000035 /)          &
      & , id_chargino(2) = (/1000024, 1000037 /)

 !---------------------------------------------
 ! internal particle codes + name for output
 !---------------------------------------------
 Integer, Save :: id_p(92)
 Character(len=12) :: c_name(92)
!-----------------------------------------------
! for decays into pions
!-----------------------------------------------
 Integer, public :: id_pi0=90, id_piP=91, id_piM=92

Contains

 Subroutine Check_Charge(id_m, id_d)
 Implicit None
  Integer, Intent(in) :: id_m, id_d(:,:)

  Integer :: e_m, e_d1, e_d2, len1, i1, sum_e, len2, e_d3

  Iname = Iname + 1
  NameOfUnit(Iname) = "Check_Charge"

  len1 = Size(id_d, dim=1)
  len2 = Size(id_d, dim=2)

  e_m = find_charge(id_p(id_m))
  If (len2.Eq.2) Then
   Do i1=1,len1
    If ((id_d(i1,1).Eq.0).And.(id_d(i1,2).Eq.0)) Exit
    e_d1 = find_charge(id_p(id_d(i1,1)))
    e_d2 = find_charge(id_p(id_d(i1,2)))
    sum_e = e_m-e_d1-e_d2
    If (sum_e.Ne.0) Then
     Write(ErrCan,*) "charge not conserved",id_p(id_m),id_p(id_d(i1,1)) &
                   &                       ,id_p(id_d(i1,2))
     Write(ErrCan,*) e_m,e_d1,e_d2,sum_e
     Call TerminateProgram()
    End If
   End Do

  Else If (len2.Eq.3) Then
   Do i1=1,len1
    If ((id_d(i1,1).Eq.0).And.(id_d(i1,2).Eq.0).And.(id_d(i1,3).Eq.0)) Exit
    e_d1 = find_charge(id_p(id_d(i1,1)))
    e_d2 = find_charge(id_p(id_d(i1,2)))
    e_d3 = find_charge(id_p(id_d(i1,3)))
    sum_e = e_m-e_d1-e_d2-e_d3
    If (sum_e.Ne.0) Then
     Write(ErrCan,*) "charge not conserved",id_p(id_m),id_p(id_d(i1,1)) &
                   &                       ,id_p(id_d(i1,2)),id_p(id_d(i1,3))
     Write(ErrCan,*) e_m,e_d1,e_d2,e_d3,sum_e
     Call TerminateProgram()
    End If
   End Do

  End If
  
  Iname = Iname - 1

 End Subroutine Check_Charge

 Integer Function Find_Charge(id)
 Implicit None
  Integer, Intent(in) :: id

  Iname = Iname + 1
  NameOfUnit(Iname) = "Find_Charge"

  Select Case(id)
   ! W-, H-,~e_L,~e_R,~mu_L,~mu_R,~tau_i,~chi+_i, e, mu, tau
   Case(24,37,-1000011,-2000011,-1000013,-2000013,-1000015,-2000015,1000024,1000037,-11,-13,-15) 
    find_charge = 3
   Case(-24,-37,1000011,2000011,1000013,2000013,1000015,2000015,-1000024,-1000037,11,13,15) 
    find_charge = -3

   ! ~d_L,~d_R,~s_L,~s_R,~b_1,~b_2, d,s,b
   Case(1000001,2000001,1000003,2000003,1000005,2000005,1,3,5)
    find_charge = -1
   Case(-1000001,-2000001,-1000003,-2000003,-1000005,-2000005,-1,-3,-5)
    find_charge = 1

   ! ~u_L,~u_R,~c_L,~c_R,~t_1,~t_2, u, c, t
   Case(1000002,2000002,1000004,2000004,1000006,2000006,2,4,6)
    find_charge = 2
   Case(-1000002,-2000002,-1000004,-2000004,-1000006,-2000006,-2,-4,-6)
    find_charge = -2

   ! nu_i, ~nu_i, g, ~g, ~chi^0_1, gamma,Z, h0, H0, A0
   Case(12,-12,14,-14,16,-16,1000012,-1000012,1000014,-1000014,1000016 &
       & ,-1000016,21,1000021,1000022,1000023,1000025,1000035,22,23,25,35,36)
    find_charge = 0

   ! scalar/pseudoscalars for MSSM extensions
   Case(1000017, 1000018, 1000019, 1000039, 45, 46)
    find_charge = 0

  Case default
   Write(Errcan,*) "Routine find_charge: unknown id",id
   Call TerminateProgram
  End Select

  Iname = Iname - 1

 End Function Find_Charge


 Subroutine Initialize_MSSM(GenerationMixing, id_gl, id_ph, id_Z, id_Wp &
               & , id_nu, id_l, id_d, id_u, id_grav)
 !------------------------------------------------------------------------------
 ! intializes particle content and internal particle numbers
 ! gives relationship to PDG and output information
 ! Note: in case of an anti-particle one has to add 1 to the particle id
 !------------------------------------------------------------------------------
 Implicit None
  Logical, Intent(in) :: GenerationMixing
  Integer, Intent(out) :: id_gl, id_ph, id_grav
  Integer, Intent(out), Dimension(1) :: id_Z, id_Wp
  Integer, Intent(out), Dimension(3) :: id_nu, id_l, id_d, id_u

  Integer :: i1
  !---------------------------------------
  ! gauge and vector bosons
  !---------------------------------------
  c_name(1) = "g"
  id_gl = 1 
  id_p(1) = id_gluon
  c_name(2) = "photon"
  id_ph = 2
  id_p(2) = id_photon
  c_name(3) = "Z"
  id_Z = 3
  id_p(3) = id_Zboson
  c_name(4) = "W^+"
  id_Wp = 4
  id_p(4) = id_Wboson
  c_name(5) = "W^-"
  id_p(5) = -id_Wboson
  !--------------------------
  ! neutrinos 
  !--------------------------
  id_nu = (/6,8,10/)
  Do i1=1,3
   id_p(4+2*i1) = id_neut(i1)
   id_p(5+2*i1) = -id_neut(i1)
  End Do
  If (GenerationMixing) Then
   c_name(6) = "nu_1"
   c_name(7) = "nu_bar_1"
   c_name(8) = "nu_2"
   c_name(9) = "nu_bar_2"
   c_name(10) = "nu_3"
   c_name(11) = "nu_bar_3"
  Else
   c_name(6) = "nu_e"
   c_name(7) = "nu_bar_e"
   c_name(8) = "nu_mu"
   c_name(9) = "nu_bar_mu"
   c_name(10) = "nu_tau"
   c_name(11) = "nu_bar_tau"
  End If
  !--------------------------
  ! leptons
  !--------------------------
  id_l = (/12,14,16/)
  Do i1=1,3
   id_p(10+2*i1) = id_lept(i1)
   id_p(11+2*i1) = -id_lept(i1)
  End Do
  c_name(12) = "e^-"
  c_name(13) = "e^+"
  c_name(14) = "mu^-"
  c_name(15) = "mu^+"
  c_name(16) = "tau^-"
  c_name(17) = "tau^+"
  !--------------------------
  ! u-quarks
  !--------------------------
  id_u = (/18,20,22/)
  Do i1=1,3
   id_p(16+2*i1) = id_u_quark(i1)
   id_p(17+2*i1) = -id_u_quark(i1)
  End Do
  c_name(18) = "u"
  c_name(19) = "u_bar"
  c_name(20) = "c"
  c_name(21) = "c_bar"
  c_name(22) = "t"
  c_name(23) = "t_bar"
  !--------------------------
  ! d-quarks
  !--------------------------
  id_d = (/24,26,28/)
  Do i1=1,3
   id_p(22+2*i1) = id_d_quark(i1)
   id_p(23+2*i1) = -id_d_quark(i1)
  End Do
  c_name(24) = "d"
  c_name(25) = "d_bar"
  c_name(26) = "s"
  c_name(27) = "s_bar"
  c_name(28) = "b"
  c_name(29) = "b_bar"
  !--------------------------
  ! sneutrinos 
  !--------------------------
  Sneut%id = (/30,32,34/)
  Do i1=1,3
   id_p(28+2*i1) = id_Sneut(i1)
   id_p(29+2*i1) = -id_Sneut(i1)
  End Do
  If (GenerationMixing) Then
   c_name(30) = "~nu_1"
   c_name(31) = "~nu^*_1"
   c_name(32) = "~nu_2"
   c_name(33) = "~nu^*_2"
   c_name(34) = "~nu_3"
   c_name(35) = "~nu^*_3"
  Else
   c_name(30) = "~nu_e"
   c_name(31) = "~nu^*_e"
   c_name(32) = "~nu_mu"
   c_name(33) = "~nu^*_mu"
   c_name(34) = "~nu_tau"
   c_name(35) = "~nu^*_tau"
  End If
  !--------------------------
  ! sleptons
  !--------------------------
  Slepton%id = (/36,38,40,42,44,46/)
  If (GenerationMixing) Then
   Do i1=1,3
    id_p(34+2*i1) = id_Slept(2*i1-1)
    id_p(40+2*i1) = id_Slept(2*i1)
   End Do
   Do i1=1,6
    id_p(35+2*i1) = -id_p(34+2*i1)
    c_name(34+2*i1) = "~l^-_"//Bu(i1)
    c_name(35+2*i1) = "~l^+_"//Bu(i1)
   End Do
  Else
   Do i1=1,6
    id_p(34+2*i1) = id_Slept(i1)
    id_p(35+2*i1) = -id_Slept(i1)
   End Do
   c_name(36) = "~e^-_L"
   c_name(37) = "~e^+_L"
   c_name(38) = "~e^-_R"
   c_name(39) = "~e^+_R"
   c_name(40) = "~mu^-_L"
   c_name(41) = "~mu^+_L"
   c_name(42) = "~mu^-_R"
   c_name(43) = "~mu^+_R"
   c_name(44) = "~tau^-_1"
   c_name(45) = "~tau^+_1"
   c_name(46) = "~tau^-_2"
   c_name(47) = "~tau^+_2"
  End If
  !--------------------------
  ! u-squarks
  !--------------------------
  Sup%id = (/48,50,52,54,56,58/)
  If (GenerationMixing) Then
   Do i1=1,3
    id_p(46+2*i1) = id_u_sq(2*i1-1)
    id_p(52+2*i1) = id_u_sq(2*i1)
   End Do
   Do i1=1,6
    id_p(47+2*i1) = -id_p(46+2*i1)
    c_name(46+2*i1) = "~u_"//Bu(i1)
    c_name(47+2*i1) = "~u^*_"//Bu(i1)
   End Do
  Else
   Do i1=1,6
    id_p(46+2*i1) = id_u_sq(i1)
    id_p(47+2*i1) = -id_u_sq(i1)
   End Do
   c_name(48) = "~u_L"
   c_name(49) = "~u^*_L"
   c_name(50) = "~u_R"
   c_name(51) = "~u^*_R"
   c_name(52) = "~c_L"
   c_name(53) = "~c^*_L"
   c_name(54) = "~c_R"
   c_name(55) = "~c^*_R"
   c_name(56) = "~t_1"
   c_name(57) = "~t^*_1"
   c_name(58) = "~t_2"
   c_name(59) = "~t^*_2"
  End If
  !--------------------------
  ! d-squarks
  !--------------------------
  Sdown%id = (/60,62,64,66,68,70/)
  If (GenerationMixing) Then
   Do i1=1,3
    id_p(58+2*i1) = id_d_sq(2*i1-1)
    id_p(64+2*i1) = id_d_sq(2*i1)
   End Do
   Do i1=1,6
    id_p(59+2*i1) = -id_p(58+2*i1)
    c_name(58+2*i1) = "~d_"//Bu(i1)
    c_name(59+2*i1) = "~d^*_"//Bu(i1)
   End Do
  Else
   Do i1=1,6
    id_p(58+2*i1) = id_d_sq(i1)
    id_p(59+2*i1) = -id_d_sq(i1)
   End Do
   c_name(60) = "~d_L"
   c_name(61) = "~d^*_L"
   c_name(62) = "~d_R"
   c_name(63) = "~d^*_R"
   c_name(64) = "~s_L"
   c_name(65) = "~s^*_L"
   c_name(66) = "~s_R"
   c_name(67) = "~s^*_R"
   c_name(68) = "~b_1"
   c_name(69) = "~b^*_1"
   c_name(70) = "~b_2"
   c_name(71) = "~b^*_2"
  End If
  !-----------------------------
  ! neutral Higgs bosons
  !-----------------------------
  S0%id = (/ 72,73 /)
  id_p(72:73) = id_h0
  c_name(72) = "h^0"
  c_name(73) = "H^0"
  !-----------------------------
  ! pseudoscalar Higgs boson
  !-----------------------------
  P0%id = 74
  id_p(74) = id_A0
  c_name(74) = "A^0"
  !-----------------------------
  ! charged Higgs boson
  !-----------------------------
  Spm%id = 75
  id_p(75:76) = (/ id_Hp, -id_Hp /)
  c_name(75) = "H^+"
  c_name(76) = "H^-"
  !-----------------------------
  ! neutralinos
  !-----------------------------
  Chi0%id = (/ 77, 78, 79, 80 /)
  id_p(77:80) = id_neutralino
  Do i1=1,4
   c_name(76+i1) = "chi^0_"//Bu(i1)
  End Do
  !-----------------------------
  ! charginos
  !-----------------------------
  ChiPm%id = (/ 81, 83 /) 
  Do i1=1,2
   id_p(79+2*i1) = id_chargino(i1)
   id_p(80+2*i1) = - id_chargino(i1)
   c_name(79+2*i1) = "chi^+_"//Bu(i1)
   c_name(80+2*i1) = "chi^-_"//Bu(i1)
  End Do
  !-----------------------------
  ! gluino
  !-----------------------------
  Glu%id = 85
  id_p(85) = id_gluino
  c_name(85) = "~g"
  !-----------------------------
  ! Gravitino
  !-----------------------------
  id_grav = 86
  Grav%id = 86
  id_p(id_grav) = id_gravitino
  c_name(id_grav) = "~G"

  !----------------------------------
  ! pions
  !----------------------------------
  id_p(id_pi0) = 111
  c_name(id_pi0) = "pi^0"
  id_p(id_piP) = 211
  c_name(id_piP) = "pi^+"
  id_p(id_piM) = -211
  c_name(id_piM) = "pi^-"

 End Subroutine Initialize_MSSM

End Module MSSM_data

Module NMSSM_data
!---------------------------------------------------------------------
! in this modul all basic information concerning the NMSSM are 
! stored, using MSSM_data and give only the additional information
! by Werner Porod, 23.09.2010
!---------------------------------------------------------------------
Use Control
Use MSSM_data ! most of the paramters are the same
 !------------------------------------------------------------------
 ! scalar masses (h,H,Re(S)), masses squared, corresponding mixing matrix, 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------
! Real(dp) :: RS03(3,3)
 Type(particle23) :: S03(3)
 !------------------------------------------------------------------------
 ! pseudoscalar masses (G,A,Im(S)), masses squared, corresponding mixing matrix 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------------
! Real(dp) :: RP03(3,3)
 Type(particle2) :: P03(3)
 !---------------------------------------------------------------------------
 ! neutralino masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: N5(5,5)
 Type(particle23) :: Chi05(5)
 !------------
 ! PDG codes
 !------------
 Integer, Parameter ::id_A02(2) = (/ id_A0, 46 /), id_h03(3) = (/ id_h0, 45 /) &
    & , id_neutralino5(5) = (/id_neutralino, 1000045 /) 
Contains

 Subroutine Initialize_NMSSM(GenerationMixing, id_gl, id_ph, id_Z, id_Wp &
               & , id_nu, id_l, id_d, id_u, id_grav)
 !------------------------------------------------------------------------------
 ! intializes particle content and internal particle numbers
 ! gives relationship to PDG and output information
 ! Note: in case of an anti-particle one has to add 1 to the particle id
 !------------------------------------------------------------------------------
 Implicit None
  Logical, Intent(in) :: GenerationMixing
  Integer, Intent(out) :: id_gl, id_ph, id_grav
  Integer, Intent(out), Dimension(1) :: id_Z, id_Wp
  Integer, Intent(out), Dimension(3) :: id_nu, id_l, id_d, id_u

  Integer :: i1
  !------------------------------------------------------------------------
  ! re-setting the number of neutralinos and neutral scalars/pseudoscalars
  !------------------------------------------------------------------------
  n_n = 5
  n_S0 = 3
  n_P0 = 3
  !---------------------------------------
  ! gauge and vector bosons
  !---------------------------------------
  c_name(1) = "g"
  id_gl = 1 
  id_p(1) = id_gluon
  c_name(2) = "photon"
  id_ph = 2
  id_p(2) = id_photon
  c_name(3) = "Z"
  id_Z = 3
  id_p(3) = id_Zboson
  c_name(4) = "W^+"
  id_Wp = 4
  id_p(4) = id_Wboson
  c_name(5) = "W^-"
  id_p(5) = -id_Wboson
  !--------------------------
  ! neutrinos 
  !--------------------------
  id_nu = (/6,8,10/)
  Do i1=1,3
   id_p(4+2*i1) = id_neut(i1)
   id_p(5+2*i1) = -id_neut(i1)
  End Do
  If (GenerationMixing) Then
   c_name(6) = "nu_1"
   c_name(7) = "nu_bar_1"
   c_name(8) = "nu_2"
   c_name(9) = "nu_bar_2"
   c_name(10) = "nu_3"
   c_name(11) = "nu_bar_3"
  Else
   c_name(6) = "nu_e"
   c_name(7) = "nu_bar_e"
   c_name(8) = "nu_mu"
   c_name(9) = "nu_bar_mu"
   c_name(10) = "nu_tau"
   c_name(11) = "nu_bar_tau"
  End If
  !--------------------------
  ! leptons
  !--------------------------
  id_l = (/12,14,16/)
  Do i1=1,3
   id_p(10+2*i1) = id_lept(i1)
   id_p(11+2*i1) = -id_lept(i1)
  End Do
  c_name(12) = "e^-"
  c_name(13) = "e^+"
  c_name(14) = "mu^-"
  c_name(15) = "mu^+"
  c_name(16) = "tau^-"
  c_name(17) = "tau^+"
  !--------------------------
  ! u-quarks
  !--------------------------
  id_u = (/18,20,22/)
  Do i1=1,3
   id_p(16+2*i1) = id_u_quark(i1)
   id_p(17+2*i1) = -id_u_quark(i1)
  End Do
  c_name(18) = "u"
  c_name(19) = "u_bar"
  c_name(20) = "c"
  c_name(21) = "c_bar"
  c_name(22) = "t"
  c_name(23) = "t_bar"
  !--------------------------
  ! d-quarks
  !--------------------------
  id_d = (/24,26,28/)
  Do i1=1,3
   id_p(22+2*i1) = id_d_quark(i1)
   id_p(23+2*i1) = -id_d_quark(i1)
  End Do
  c_name(24) = "d"
  c_name(25) = "d_bar"
  c_name(26) = "s"
  c_name(27) = "s_bar"
  c_name(28) = "b"
  c_name(29) = "b_bar"
  !--------------------------
  ! sneutrinos 
  !--------------------------
  Sneut%id = (/30,32,34/)
  Do i1=1,3
   id_p(28+2*i1) = id_Sneut(i1)
   id_p(29+2*i1) = -id_Sneut(i1)
  End Do
  If (GenerationMixing) Then
   c_name(30) = "~nu_1"
   c_name(31) = "~nu^*_1"
   c_name(32) = "~nu_2"
   c_name(33) = "~nu^*_2"
   c_name(34) = "~nu_3"
   c_name(35) = "~nu^*_3"
  Else
   c_name(30) = "~nu_e"
   c_name(31) = "~nu^*_e"
   c_name(32) = "~nu_mu"
   c_name(33) = "~nu^*_mu"
   c_name(34) = "~nu_tau"
   c_name(35) = "~nu^*_tau"
  End If
  !--------------------------
  ! sleptons
  !--------------------------
  Slepton%id = (/36,38,40,42,44,46/)
  If (GenerationMixing) Then
   Do i1=1,3
    id_p(34+2*i1) = id_Slept(2*i1-1)
    id_p(40+2*i1) = id_Slept(2*i1)
   End Do
   Do i1=1,6
    id_p(35+2*i1) = -id_p(34+2*i1)
    c_name(34+2*i1) = "~l^-_"//Bu(i1)
    c_name(35+2*i1) = "~l^+_"//Bu(i1)
   End Do
  Else
   Do i1=1,6
    id_p(34+2*i1) = id_Slept(i1)
    id_p(35+2*i1) = -id_Slept(i1)
   End Do
   c_name(36) = "~e^-_L"
   c_name(37) = "~e^+_L"
   c_name(38) = "~e^-_R"
   c_name(39) = "~e^+_R"
   c_name(40) = "~mu^-_L"
   c_name(41) = "~mu^+_L"
   c_name(42) = "~mu^-_R"
   c_name(43) = "~mu^+_R"
   c_name(44) = "~tau^-_1"
   c_name(45) = "~tau^+_1"
   c_name(46) = "~tau^-_2"
   c_name(47) = "~tau^+_2"
  End If
  !--------------------------
  ! u-squarks
  !--------------------------
  Sup%id = (/48,50,52,54,56,58/)
  If (GenerationMixing) Then
   Do i1=1,3
    id_p(46+2*i1) = id_u_sq(2*i1-1)
    id_p(52+2*i1) = id_u_sq(2*i1)
   End Do
   Do i1=1,6
    id_p(47+2*i1) = -id_p(46+2*i1)
    c_name(46+2*i1) = "~u_"//Bu(i1)
    c_name(47+2*i1) = "~u^*_"//Bu(i1)
   End Do
  Else
   Do i1=1,6
    id_p(46+2*i1) = id_u_sq(i1)
    id_p(47+2*i1) = -id_u_sq(i1)
   End Do
   c_name(48) = "~u_L"
   c_name(49) = "~u^*_L"
   c_name(50) = "~u_R"
   c_name(51) = "~u^*_R"
   c_name(52) = "~c_L"
   c_name(53) = "~c^*_L"
   c_name(54) = "~c_R"
   c_name(55) = "~c^*_R"
   c_name(56) = "~t_1"
   c_name(57) = "~t^*_1"
   c_name(58) = "~t_2"
   c_name(59) = "~t^*_2"
  End If
  !--------------------------
  ! d-squarks
  !--------------------------
  Sdown%id = (/60,62,64,66,68,70/)
  If (GenerationMixing) Then
   Do i1=1,3
    id_p(58+2*i1) = id_d_sq(2*i1-1)
    id_p(64+2*i1) = id_d_sq(2*i1)
   End Do
   Do i1=1,6
    id_p(59+2*i1) = -id_p(58+2*i1)
    c_name(58+2*i1) = "~d_"//Bu(i1)
    c_name(59+2*i1) = "~d^*_"//Bu(i1)
   End Do
  Else
   Do i1=1,6
    id_p(58+2*i1) = id_d_sq(i1)
    id_p(59+2*i1) = -id_d_sq(i1)
   End Do
   c_name(60) = "~d_L"
   c_name(61) = "~d^*_L"
   c_name(62) = "~d_R"
   c_name(63) = "~d^*_R"
   c_name(64) = "~s_L"
   c_name(65) = "~s^*_L"
   c_name(66) = "~s_R"
   c_name(67) = "~s^*_R"
   c_name(68) = "~b_1"
   c_name(69) = "~b^*_1"
   c_name(70) = "~b_2"
   c_name(71) = "~b^*_2"
  End If
  !-----------------------------
  ! neutral Higgs bosons
  !-----------------------------
  S03%id = (/ 72,73, 74 /)
  id_p(72:74) = id_h03
  c_name(72) = "H^0_1"
  c_name(73) = "H^0_2"
  c_name(74) = "H^0_2"
  !-----------------------------
  ! pseudoscalar Higgs boson
  !-----------------------------
  P03%id = (/ 75, 75, 76 /)
  id_p(75:76) = id_A02
  c_name(75) = "A^0_1"
  c_name(76) = "A^0_2"
  !-----------------------------
  ! charged Higgs boson
  !-----------------------------
  Spm%id = 77
  id_p(77:78) = (/ id_Hp, -id_Hp /)
  c_name(77) = "H^+"
  c_name(78) = "H^-"
  !-----------------------------
  ! neutralinos
  !-----------------------------
  Chi05%id = (/ 79, 80, 81, 82, 83 /)
  id_p(79:83) = id_neutralino5
  Do i1=1,5
   c_name(78+i1) = "chi^0_"//Bu(i1)
  End Do
  !-----------------------------
  ! charginos
  !-----------------------------
  ChiPm%id = (/ 84, 86 /) 
  Do i1=1,2
   id_p(82+2*i1) = id_chargino(i1)
   id_p(83+2*i1) = - id_chargino(i1)
   c_name(82+2*i1) = "chi^+_"//Bu(i1)
   c_name(83+2*i1) = "chi^-_"//Bu(i1)
  End Do
  !-----------------------------
  ! gluino
  !-----------------------------
  Glu%id = 88
  id_p(88) = id_gluino
  c_name(88) = "~g"
  !-----------------------------
  ! Gravitino
  !-----------------------------
  id_grav = 89
  id_p(id_grav) = id_gravitino
  c_name(id_grav) = "~G"

 End Subroutine Initialize_NMSSM

End Module NMSSM_data
Module RP_data
!---------------------------------------------------------------------
! in this modul all basic information concerning the NMSSM are 
! stored, using MSSM_data and give only the additional information
! by Werner Porod, 23.09.2010
!---------------------------------------------------------------------
Use Control
Use MSSM_data ! most of the paramters are the same
 !------------------------------------------------------------------
 ! scalar masses (h,H), masses squared, corresponding mixing matrix, 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------
! Real(dp) :: RS05(5,5)
 Type(particle23) :: S05(5)
 !------------------------------------------------------------------------
 ! pseudoscalar masses (G,A), masses squared, corresponding mixing matrix 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------------
! Real(dp) :: RP05(5,5)
 Type(particle2) :: P05(5)
 !---------------------------------------------------------------------------
 ! charged scalar masses (G+,H+.~l_L,~l_R), masses squared,
 ! corresponding mixing matrix 
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: RSpm8(8,8)
 Type(particle2) :: Spm8(8)
 !---------------------------------------------------------------------------
 ! neutralino masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: N7(7,7)
 Type(particle23) :: Chi07(7)
 !---------------------------------------------------------------------------
 ! lepton/chargino masses, masses squared, corresponding mixing matrices
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
! Complex(dp) :: U5(5,5), V5(5,5)
 Type(particle23) :: ChiPM5(5)
 !---------------------------------------------------------------------------
 ! check if trilinear couplings are present
 !---------------------------------------------------------------------------
 Logical, Save :: RP_trilinear = .False.
 !------------
 ! PDG codes
 !------------
 Integer, Parameter :: id_sleptp(6) = -id_slept , id_leptp(3) = -id_lept
 Integer, Parameter :: id_P05(4) = (/ id_A0, 1000017, 1000018, 1000019 /)      &
    & , id_S05(5) = (/ id_h0, id_sneut /), id_Spm8(7) = (/ id_hp, id_sleptp /) &
    & , id_neutralino7(7) = (/id_neut, id_neutralino /)                        &
    & , id_chargino5(5) = (/id_leptp, id_chargino /)
Contains

 Subroutine Initialize_RPexplicit(GenerationMixing, id_gl, id_ph, id_Z, id_Wp &
               & , id_nu, id_l, id_d, id_u, id_grav)
 !------------------------------------------------------------------------------
 ! intializes particle content and internal particle numbers
 ! gives relationship to PDG and output information
 ! Note: in case of an anti-particle one has to add 1 to the particle id
 !------------------------------------------------------------------------------
 Implicit None
  Logical, Intent(in) :: GenerationMixing
  Integer, Intent(out) :: id_gl, id_ph, id_grav
  Integer, Intent(out), Dimension(1) :: id_Z, id_Wp
  Integer, Intent(out), Dimension(3) :: id_nu, id_l, id_d, id_u

  Integer :: i1
  id_p = 0
  Do i1=1,6
   Sup(i1)%id2 = 0
   Sup(i1)%id3 = 0
  End Do
  !---------------------------------------
  ! gauge and vector bosons
  !---------------------------------------
  c_name(1) = "g"
  id_gl = 1 
  id_p(1) = id_gluon
  c_name(2) = "photon"
  id_ph = 2
  id_p(2) = id_photon
  c_name(3) = "Z"
  id_Z = 3
  id_p(3) = id_Zboson
  c_name(4) = "W^+"
  id_Wp = 4
  id_p(4) = id_Wboson
  c_name(5) = "W^-"
  id_p(5) = -id_Wboson
  !--------------------------
  ! u-quarks
  !--------------------------
  id_u = (/6,8,10/)
  Do i1=1,3
   id_p(4+2*i1) = id_u_quark(i1)
   id_p(5+2*i1) = -id_u_quark(i1)
  End Do
  c_name(6) = "u"
  c_name(7) = "u_bar"
  c_name(8) = "c"
  c_name(9) = "c_bar"
  c_name(10) = "t"
  c_name(11) = "t_bar"
  !--------------------------
  ! d-quarks
  !--------------------------
  id_d = (/12,14,16/)
  Do i1=1,3
   id_p(10+2*i1) = id_d_quark(i1)
   id_p(11+2*i1) = -id_d_quark(i1)
  End Do
  c_name(12) = "d"
  c_name(13) = "d_bar"
  c_name(14) = "s"
  c_name(15) = "s_bar"
  c_name(16) = "b"
  c_name(17) = "b_bar"
  !--------------------------
  ! u-squarks
  !--------------------------
  Sup%id = (/18,20,22,24,26,28/)
  If (GenerationMixing) Then
   Do i1=1,3
    id_p(16+2*i1) = id_u_sq(2*i1-1)
    id_p(522+2*i1) = id_u_sq(2*i1)
   End Do
   Do i1=1,6
    id_p(17+2*i1) = -id_p(16+2*i1)
    c_name(16+2*i1) = "~u_"//Bu(i1)
    c_name(17+2*i1) = "~u^*_"//Bu(i1)
   End Do
  Else
   Do i1=1,6
    id_p(16+2*i1) = id_u_sq(i1)
    id_p(17+2*i1) = -id_u_sq(i1)
   End Do
   c_name(18) = "~u_L"
   c_name(19) = "~u^*_L"
   c_name(20) = "~u_R"
   c_name(21) = "~u^*_R"
   c_name(22) = "~c_L"
   c_name(23) = "~c^*_L"
   c_name(24) = "~c_R"
   c_name(25) = "~c^*_R"
   c_name(26) = "~t_1"
   c_name(27) = "~t^*_1"
   c_name(28) = "~t_2"
   c_name(29) = "~t^*_2"
  End If
  !--------------------------
  ! d-squarks
  !--------------------------
  Sdown%id = (/30,32,34,36,38,40/)
  If (GenerationMixing) Then
   Do i1=1,3
    id_p(28+2*i1) = id_d_sq(2*i1-1)
    id_p(34+2*i1) = id_d_sq(2*i1)
   End Do
   Do i1=1,6
    id_p(29+2*i1) = -id_p(28+2*i1)
    c_name(28+2*i1) = "~d_"//Bu(i1)
    c_name(29+2*i1) = "~d^*_"//Bu(i1)
   End Do
  Else
   Do i1=1,6
    id_p(28+2*i1) = id_d_sq(i1)
    id_p(29+2*i1) = -id_d_sq(i1)
   End Do
   c_name(30) = "~d_L"
   c_name(31) = "~d^*_L"
   c_name(32) = "~d_R"
   c_name(33) = "~d^*_R"
   c_name(34) = "~s_L"
   c_name(35) = "~s^*_L"
   c_name(36) = "~s_R"
   c_name(37) = "~s^*_R"
   c_name(38) = "~b_1"
   c_name(39) = "~b^*_1"
   c_name(40) = "~b_2"
   c_name(41) = "~b^*_2"
  End If
  !----------------------------------------
  ! neutral Higgs bosons/ Real(sneutrinos)
  !----------------------------------------
  S05%id = (/ 42, 43, 44, 45, 46 /)
  id_p(42:46) = id_S05
  Do i1=1,5
   c_name(41+i1) = "S^0_"//Bu(i1)
  End Do
  !----------------------------------------
  ! pseudoscalar Higgs boson/ Im(sneutrinos)
  !----------------------------------------
  P05%id = (/ 47, 47, 48, 49, 50 /)
  id_p(47:50) = id_P05
  Do i1=1,4
   c_name(46+i1) = "P^0_"//Bu(i1)
  End Do
  !-----------------------------
  ! charged Higgs boson/sleptons
  !-----------------------------
  Spm8%id = (/51, 51, 53, 55, 57, 59, 61, 63/)
  Do i1=1,7
   id_p(49+2*i1) = id_Spm8(i1)
   id_p(50+2*i1) = -id_Spm8(i1)
   c_name(49+2*i1) = "S^+_"//Bu(i1)
   c_name(50+2*i1) = "S^-_"//Bu(i1)
  End Do
  !-----------------------------
  ! neutrinos/neutralinos
  !-----------------------------
  Chi07%id = (/ 65, 66, 67, 68, 69, 70, 71 /)
  id_nu = Chi07(1:3)%id
  id_p(65:71) = id_neutralino7
  Do i1=1,3
   c_name(64+i1) = "nu_"//Bu(i1)
  End Do
  Do i1=1,4
   c_name(67+i1) = "chi^0_"//Bu(i1)
  End Do
  !-----------------------------
  ! leptons/charginos
  !-----------------------------
  ChiPm5%id = (/ 72, 74, 76, 78, 80 /) 
  id_l = ChiPm5(1:3)%id
  c_name(73) = "e^-"
  c_name(72) = "e^+"
  c_name(75) = "mu^-"
  c_name(74) = "mu^+"
  c_name(77) = "tau^-"
  c_name(76) = "tau^+"
  Do i1=1,5
   id_p(70+2*i1) = id_chargino5(i1)
   id_p(71+2*i1) = - id_chargino5(i1)
   If (i1.Gt.3) Then
    c_name(70+2*i1) = "chi^+_"//Bu(i1-3)
    c_name(71+2*i1) = "chi^-_"//Bu(i1-3)
   End If
  End Do
  !-----------------------------
  ! gluino
  !-----------------------------
  Glu%id = 88
  id_p(88) = id_gluino
  c_name(88) = "~g"
  !-----------------------------
  ! Gravitino
  !-----------------------------
  id_grav = 89
  id_p(id_grav) = id_gravitino
  c_name(id_grav) = "~G"

 End Subroutine Initialize_RPexplicit

End Module RP_data
