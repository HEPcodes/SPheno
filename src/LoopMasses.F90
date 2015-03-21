Module LoopMasses
! In this Module all 1-loop self energies are collected including the
! tad-poles. Starting point is the paper by Bagger et al., NPB 491, 3, (1997). 
! porting routines for two-loop contributions to Higgs masses provided 
! by Pietro Slavich

! load modules
Use Control
Use Mathematics, Only: Li2, Delt
Use MathematicsQP
Use LoopFunctions, Only: A0, B0, B1, B22, Floop, Gloop, HLoop, C0_3m &
    &        , D0_Bagger, D27_Bagger, GetRenormalizationScale, Phi
Use RGEs
Use StandardModel
Use Couplings
Use SusyMasses

#ifdef THREELOOP
Use ThreeLoopMasses
#endif
Use TwoLoopHiggsMass
Use Model_Data
! load modules

! public variables
 Real(dp) :: tanb_Q, vev_Q, mA2_Q
 Logical :: SPA_Convention
 Logical, Save :: LoopContributions(8) = .True.

Logical :: test_write
! private variables
 Integer, Private :: WriteOneLoopContributions = 0 &
        & , n_char, n_neut, n_P0, n_S0, n_Spm, n_Slept, n_Sneut
 Real(dp), Private :: mC_1L(2), mN_1L(4), mS0_1L(2), mP0_1L(2)           &
    & , mUsquark_1L(6), mDsquark_1L(6), mSlepton_1L(6), mSneutrino_1L(3)      &
    & , mUsquark2_1L(6), mDsquark2_1L(6), mSlepton2_1L(6), mSneutrino2_1L(3)  &
    & , mS02_1L(2), mP02_1L(2), mSpm2_1L(2), mglu_1L, mSpm_1L(2), RS0_1L(2,2) &
    & , mC2_1L(2), mN2_1L(4)
 Real(dp), Private :: mudim, WriteWert=0
 Real(dp), Private, Dimension(3) :: mf_d_save, mf_u_save, mf_l_save 
 Complex(dp), Private :: U_1L(2,2), V_1L(2,2), N_1L(4,4), RUSquark_1L(6,6)    &
    &  , RDSquark_1L(6,6), RSlepton_1L(6,6), RSneutrino_1L(3,3)

 Logical, Save :: Only_1loop_Higgsmass=.False., ThreeLoopHiggsMass=.False.

! numerical constants
 Integer, Private, Parameter :: Nc = 3
! private variables

Contains


 Subroutine LoopMassesMSSM(delta0, tanb_mZ, tanb_Q, g, Y_l, Y_d, Y_u, Mi      &
    & , A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, phase_mu, mu, B    &
    & , ic, uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                                &
    & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP0, mP02, RP0             &
    & , mSpm, mSpm2, RSpm, mDsquark, mDsquark2, RDsquark, mUsquark, mUsquark2 &
    & , RUsquark, mSlepton, mSlepton2, RSlepton, mSneutrino, mSneutrino2      &
    & , RSneutrino, mGlu, phase_glu, kont)
 !-----------------------------------------------------------------------
 ! In this Subroutine the 1-loop masses within the MSSM are calculated
 ! It is assumed that the parameters are provided at the same renormalisation
 ! scale mudim
 ! The corresponding Subroutine from LoopMasses.F is taken as basis
 ! written by Werner Porod
 ! 08.10.01: - defining Interface
 !           - first Call of TreeMassesMSSM and AllCouplings implemented
 !             here it is still assumed that there is no Generation mixing
 !           - implementing the running of Tan(beta)
 !           - implementing PiZZT1
 ! 09.10.01: - implementing One_Loop_Tadpoles_MSSM
 !           - implementing PiPseudoscalar
 ! 10.10.01: - implementing PiChargedScalar
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: ic
  Integer, Intent(inout) :: kont
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3), Mi(3), A_l(3,3)   &
    & , A_d(3,3), A_u(3,3), M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_Q(3,3)       &
    & , M2_U(3,3), phase_mu
  Real(dp), Intent(in) :: delta0, g(3), M2_H(2)
  Real(dp), Intent(inout) :: tanb_mZ, tanb_Q

  Real(dp), Intent(out) :: mC(2), mN(4), mS0(2), mP0(2), mSpm(2)   &
    & , mUsquark(6), mDsquark(6), mSlepton(6), mSneutrino(3)       &
    & , mUsquark2(6), mDsquark2(6), mSlepton2(6), mSneutrino2(3)   &
    & , mS02(2), mP02(2), mSpm2(2), RP0(2,2), mglu, mC2(2), mN2(4), RS0(2,2)
  Complex(dp), Intent(out) :: U(2,2), V(2,2), N(4,4), RSpm(2,2)           &
    & , RDsquark(6,6), RUsquark(6,6), RSlepton(6,6), RSneutrino(3,3)      &
    & , phase_Glu, mu, B, uU_L(3,3), uU_R(3,3) ,uD_L(3,3), uD_R(3,3)      &
    & , uL_L(3,3), uL_R(3,3)

  Real(dp) :: vevSM(2), tanbQ, cosb2, cos2b, sinb2, Abs_Mu2               &
    & , Abs_Mu, g57(57), g58(58), Q, dt, tz, vev2, sinW2_DR               &
    & , mZ2_mZ, vevs_DR(2), tadpoles_1L(2), M2_H_1L(2), p2, b_A, mW2_mW   &
    & , mA2_mA, eq, T3, Q2, Pi2A0, tadpoles_2L(2), mW2_run, mglu_sign     &
    & , comp2(2), wert
  Real(dp) :: M2Ein, M2Lin, M2Din, M2Qin, M2Uin, delta1
  Complex(dp) :: dmZ2, PiA0, dmW2, PiSpm, M1, M2, dmglu, B_T, mu_T
  Complex(dp), Dimension(3,3) :: CKM_Q
  Integer :: i1, i2, i_count
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'LoopMassesMSSM'

  kont = 0
  mudim = GetRenormalizationScale() ! from LoopFunctions
  !---------------------------------
  ! definning model
  !---------------------------------
  n_S0 = 2
  n_P0 = 2
  n_Spm = 2
  n_char = 2
  n_neut = 4
  n_Sneut = 3
  n_Slept = 6

  M1 = Mi(1)
  M2 = Mi(2)
  !----------------------------------------------------------------------------
  ! require internally a factor 10 smaller delta to be sure, except of smaller
  ! than 100*Epsilon(1._dp) to avoid numerical problems
  !----------------------------------------------------------------------------
  delta1 = delta0
  If (delta1.Gt.Epsilon(1._dp)*100._dp) delta1 = 0.1_dp * delta1
  !---------------------------
  ! running of Tan(beta)
  !---------------------------
  If (tanb_in_at_Q) Then
   Q = Sqrt(mudim)
   Q2 = mudim
   tanbQ = tanb_Q
   Call CouplingsToG(g,Y_l,Y_d,Y_u,g58(1:57))
   g58(1) = Sqrt(5._dp / 3._dp ) * g58(1)  ! rescaling
   g58(58) = Log(tanb_Q)
   tz = Log(Q/mZ)
   dt = - tz / 50._dp
   Call odeint(g58, 58, tz, 0._dp, delta1, dt, 0._dp, rge58, kont)

   tanb_mZ = Exp(g58(58))
  Else If (mudim.Ne.mZ2) Then
   Q = Sqrt(mudim)
   Q2 = mudim
   Call CouplingsToG(g,Y_l,Y_d,Y_u,g57)
   g57(1) = Sqrt(5._dp / 3._dp ) * g57(1)  ! rescaling

   !------------------------------------
   ! calculate first the couplings at mZ
   !------------------------------------
   tz = Log(Q/mZ)
   dt = - tz / 50._dp
   
   Call odeint(g57, 57, tz, 0._dp,  delta1, dt, 0._dp, rge57, kont)
   !---------------------------------------
   ! and now back to Q including Tan(beta)
   !---------------------------------------
   g58(1:57) = g57
   g58(58) = Log( tanb_mZ )
   tz = Log(Q/mZ)
   dt = tz / 50._dp

   Call odeint(g58, 58, 0._dp, tz,  delta1 , dt, 0._dp, rge58, kont)

   tanbQ = Exp( g58(58) )
   tanb_Q = tanbQ
  Else
   tanbQ = tanb_mZ
   tanb_Q = tanbQ
   Q2 = mZ2
  End If

  cosb2 = 1._dp / (1._dp + tanbQ**2)
  sinb2 = 1._dp - cosb2
  cos2b = cosb2 - sinb2

  vev2 =  4._dp * mZ2 / (g(1)**2 + g(2)**2)
  vevSM(1) = Sqrt(vev2 * cosb2 )
  vevSM(2) = tanbQ * vevSM(1)
  !------------------------------------------------
  ! replacing fermion pole masses by running masses
  !------------------------------------------------
  mf_d_save = mf_d
  mf_l_save = mf_l
  mf_u_save = mf_u
  If (GenerationMixing) Then
   Call FermionMass(Y_l,vevSM(1),mf_l,uL_L, uL_R,kont)
   Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevSM, mf_d, uD_L, uD_R &
                                     & , mf_u, uU_L, uU_R, CKM_Q)
  Else
   Do i1=1,3
    mf_d(i1) = oosqrt2 * Abs(Y_d(i1,i1)) * vevSM(1)
    mf_l(i1) = oosqrt2 * Abs(Y_l(i1,i1)) * vevSM(1)
    mf_u(i1) = oosqrt2 * Abs(Y_u(i1,i1)) * vevSM(2)
   End Do
   CKM_Q = id3C
  End If
  mf_d2 = mf_d**2
  mf_l2 = mf_l**2
  mf_u2 = mf_u**2
  !-----------------------------------
  ! mu and B parameters at tree-level
  !-----------------------------------
  Abs_Mu2 = (M2_H(2) * sinb2 - M2_H(1) * cosb2) / cos2b  - 0.5_dp * mZ2
  If (Abs_Mu2.Le.0) Then
   WriteWert = Abs_Mu2
   Call WriteLoopMassesError(1, "LoopMassesMSSM", kont)
   Abs_Mu2 = Abs(Abs_Mu2)
  Else If (Abs_Mu2.Ge.1.e20_dp) Then
   WriteWert = Abs_Mu2
   Call WriteLoopMassesError(2, "LoopMassesMSSM", kont)
   Abs_Mu2 = 1.e4_dp
  End If

  Abs_Mu = Sqrt( Abs_Mu2 )
  mu = Abs_mu * phase_mu
  B = (M2_H(1) + M2_H(2) + 2._dp *  Abs_Mu2) * tanbQ / (1+tanbQ**2)
  If (Real(B).Le.0._dp) B = -B+1 ! problem with large tan(beta)
  mu_T = mu
  B_T = B
  !------------------------------------------------
  ! now some iterations to get all values precisely
  !------------------------------------------------
  i_count = 0
  mZ2_mZ = mZ2
  mW2_run = 0.25_dp * g(2)**2 * vev2

  comp2(1) = Real(mu,dp)
  comp2(2) = mZ2_mZ
  vevs_DR = vevSM

  !-------------------------------------------------
  ! tree level masses and mixings as starting point
  !-------------------------------------------------
   kont = 0
   Call TreeMassesMSSM2(g(1), g(2), vevs_DR, Mi(1), Mi(2), Mi(3), mu, B      &
      & , tanbQ, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u  &
      & , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                                 &
      & , mGlu, Phase_Glu, mC, mC2, U, V, mN, mN2, N                         &
      & , mSneutrino, mSneutrino2, Rsneutrino, mSlepton, mSlepton2, RSlepton &
      & , mDSquark, mDSquark2, RDSquark, mUSquark, mUSquark2, RUSquark       &
      & , mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm, mZ2_mZ, mW2_run &
      & , GenerationMixing, kont, .True., .False.)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If
   mglu_sign = phase_glu * mglu
   kont = 0
   sinW2_DR = g(1)**2 / (g(1)**2+g(2)**2)

  !----------------------
  ! mZ(mZ)
  !----------------------
   Call PiZZT1(mZ2, g(2), sinW2_DR, vevs_DR, mZ2, mW2, mS02, RS0, mP02, RP0 &
      &  , mSpm2, RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton        &
      &  , mUSquark2, RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2    &
      &  , mC, mC2, U, V, mN, mN2, N ,dmZ2)
   vev2 =  4._dp * Real(mZ2+dmZ2,dp) / (g(1)**2 + g(2)**2)
   vevs_DR(1) = Sqrt(vev2 / (1._dp+tanbQ**2) )
   vevs_DR(2) = tanbQ * vevs_DR(1)
   mZ2_mZ = Real(dmZ2+mZ2,dp)
   mW2_run = mZ2_mZ * (1._dp - sinW2_DR)

   If (mZ2_mZ.Lt.0._dp) Then
    WriteWert = mZ2_mZ
    Call WriteLoopMassesError(4, "LoopMassesMSSM", kont)
    mZ2_mZ = mZ2
   End If
   !------------------------------------------------
   ! replacing fermion pole masses by running masses
   !------------------------------------------------
   If (GenerationMixing) Then
    Call FermionMass(Y_l,vevs_DR(1),mf_l,uL_L, uL_R,kont)
    Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevs_DR, mf_d, uD_L, uD_R &
                                  & , mf_u, uU_L, uU_R, CKM_Q)
   Else
    Do i1=1,3
     mf_d(i1) = oosqrt2 * Abs(Y_d(i1,i1)) * vevs_DR(1)
     mf_l(i1) = oosqrt2 * Abs(Y_l(i1,i1)) * vevs_DR(1)
     mf_u(i1) = oosqrt2 * Abs(Y_u(i1,i1)) * vevs_DR(2)
    End Do
    CKM_Q = id3C
   End If
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
  !-----------------------------------
  ! mu and B parameters at tree-level
  !-----------------------------------
  Abs_Mu2 = (M2_H(2) * sinb2 - M2_H(1) * cosb2) / cos2b  - 0.5_dp * mZ2_mZ
  If (Abs_Mu2.Le.0) Then
   WriteWert = Abs_Mu2
   Call WriteLoopMassesError(1, "LoopMassesMSSM", kont)
   Abs_Mu2 = Abs(Abs_Mu2)
  Else If (Abs_Mu2.Ge.1.e20_dp) Then
   WriteWert = Abs_Mu2
   Call WriteLoopMassesError(2, "LoopMassesMSSM", kont)
   Abs_Mu2 = 1.e4_dp
  End If

  Abs_Mu = Sqrt( Abs_Mu2 )
  mu = Abs_mu * phase_mu
  B = (M2_H(1) + M2_H(2) + 2._dp *  Abs_Mu2) * tanbQ / (1+tanbQ**2)
  If (Real(B).Le.0._dp) B = -B+1 ! problem with large tan(beta)
  mu_T = mu
  B_T = B
  !------------------------------------------------
  ! tree level masses and mixings as starting point
  !-------------------------------------------------
   kont = 0
   Call TreeMassesMSSM2(g(1), g(2), vevs_DR, Mi(1), Mi(2), Mi(3), mu, B      &
      & , tanbQ, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u  &
      & , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                                 &
      & , mGlu, Phase_Glu, mC, mC2, U, V, mN, mN2, N                         &
      & , mSneutrino, mSneutrino2, Rsneutrino, mSlepton, mSlepton2, RSlepton &
      & , mDSquark, mDSquark2, RDSquark, mUSquark, mUSquark2, RUSquark       &
      & , mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm, mZ2_mZ, mW2_run &
      & , GenerationMixing, kont, .True., .False.)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If

  !----------------------
  ! mZ(mZ)
  !----------------------
   Call PiZZT1(mZ2, g(2), sinW2_DR, vevs_DR, mZ2, mW2, mS02, RS0, mP02, RP0 &
      &  , mSpm2, RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton        &
      &  , mUSquark2, RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2    &
      &  , mC, mC2, U, V, mN, mN2, N ,dmZ2)
   vev2 =  4._dp * Real(mZ2+dmZ2,dp) / (g(1)**2 + g(2)**2)
   vevs_DR(1) = Sqrt(vev2 / (1._dp+tanbQ**2) )
   vevs_DR(2) = tanbQ * vevs_DR(1)
   mZ2_mZ = Real(dmZ2+mZ2,dp)
   mW2_run = mZ2_mZ * (1._dp - sinW2_DR)

   If (mZ2_mZ.Lt.0._dp) Then
    WriteWert = mZ2_mZ
    Call WriteLoopMassesError(4, "LoopMassesMSSM", kont)
    mZ2_mZ = mZ2
   End If
   !------------------------------------------------
   ! replacing fermion pole masses by running masses
   !------------------------------------------------
   If (GenerationMixing) Then
    Call FermionMass(Y_l,vevs_DR(1),mf_l,uL_L, uL_R,kont)
    Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevs_DR, mf_d, uD_L, uD_R &
                                  & , mf_u, uU_L, uU_R, CKM_Q)
   Else
    Do i1=1,3
     mf_d(i1) = oosqrt2 * Abs(Y_d(i1,i1)) * vevs_DR(1)
     mf_l(i1) = oosqrt2 * Abs(Y_l(i1,i1)) * vevs_DR(1)
     mf_u(i1) = oosqrt2 * Abs(Y_u(i1,i1)) * vevs_DR(2)
    End Do
    CKM_Q = id3C
   End If
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
  !-----------------------------------
  ! mu and B parameters at tree-level
  !-----------------------------------
  Abs_Mu2 = (M2_H(2) * sinb2 - M2_H(1) * cosb2) / cos2b  - 0.5_dp * mZ2_mZ
  If (Abs_Mu2.Le.0) Then
   WriteWert = Abs_Mu2
   Call WriteLoopMassesError(1, "LoopMassesMSSM", kont)
   Abs_Mu2 = Abs(Abs_Mu2)
  Else If (Abs_Mu2.Ge.1.e20_dp) Then
   WriteWert = Abs_Mu2
   Call WriteLoopMassesError(2, "LoopMassesMSSM", kont)
   Abs_Mu2 = 1.e4_dp
  End If

  Abs_Mu = Sqrt( Abs_Mu2 )
  mu = Abs_mu * phase_mu
  B = (M2_H(1) + M2_H(2) + 2._dp *  Abs_Mu2) * tanbQ / (1+tanbQ**2)
  If (Real(B).Le.0._dp) B = -B+1 ! problem with large tan(beta)
  mu_T = mu
  B_T = B
  !------------------------------------------------
  ! tree level masses and mixings as starting point
  !-------------------------------------------------
   kont = 0
   Call TreeMassesMSSM2(g(1), g(2), vevs_DR, Mi(1), Mi(2), Mi(3), mu, B      &
      & , tanbQ, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u  &
      & , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                                 &
      & , mGlu, Phase_Glu, mC, mC2, U, V, mN, mN2, N                         &
      & , mSneutrino, mSneutrino2, Rsneutrino, mSlepton, mSlepton2, RSlepton &
      & , mDSquark, mDSquark2, RDSquark, mUSquark, mUSquark2, RUSquark       &
      & , mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm, mZ2_mZ, mW2_run &
      & , GenerationMixing, kont, .True., .False.)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If

  !---------------------------
  ! tadpoles at 1-loop
  !---------------------------
#ifdef GENERATIONMIXING
   Call One_Loop_Tadpoles_MSSM(g(1), g(2), vevs_DR, Y_l, Y_d, Y_u, uU_L, uU_R &
    & , uD_L, uD_R, uL_L, uL_R, mu_T, A_l, A_d, A_u, mSneutrino2, Rsneutrino  &
    & , mSlepton2, Rslepton, mDSquark2, RDSquark, mUSquark2, RUSquark         &
    & , mSpm2, RSpm, mC, mC2, U, V, mP02, RP0, mS02, RS0, mN, mN2, N          &
    & , mZ2_mZ, mW2_run, tadpoles_1L)
#else
   Call One_Loop_Tadpoles_MSSM(g(1), g(2), vevs_DR, Y_l, Y_d, Y_u, mu_T, A_l  &
    & , A_d, A_u, mSneutrino2, mSlepton2, Rslepton, mDSquark2, RDSquark       &
    & , mUSquark2, RUSquark, mSpm2, RSpm, mC, mC2, U, V, mP02, RP0, mS02, RS0 &
    & , mN, mN2, N, mZ2_mZ, mW2_run, tadpoles_1L)
#endif
  !---------------------------
  ! tadpoles at 2-loop
  !---------------------------------------------------------------------
  ! Define first input, numerical problems can occur if left and right
  ! sfermion mass parameters are identical
  !---------------------------------------------------------------------
   M2Din = Real(M2_D(3,3),dp)
   M2Ein = Real(M2_E(3,3),dp)
   M2Lin = Real(M2_L(3,3),dp)
   M2Qin = Real(M2_Q(3,3),dp)
   M2Uin = Real(M2_U(3,3),dp)
   If (M2Ein.Eq.M2Lin) M2Ein = (1._dp + Max(1.0e-8_dp,delta1)) * M2Ein
   If (M2Din.Eq.M2Qin) M2Din = (1._dp + Max(1.0e-8_dp,delta1)) * M2Din
   If (M2Uin.Eq.M2Qin) M2Uin = (1._dp + Max(1.0e-8_dp,delta1)) * M2Uin

   If (Only_1loop_Higgsmass) Then
    tadpoles_2L = 0._dp
   Else
    Call Two_Loop_Tadpoles_MSSM(g(3), mglu_sign, mP02(2), vevs_DR          &
          & , M2Din, M2Uin, M2Qin, M2Ein, M2Lin, A_d(3,3), A_u(3,3)        &
          & , A_l(3,3), Y_d(3,3), Y_u(3,3), Y_l(3,3), mu_T, tadpoles_2L, kont)
   End If

   If ((ic.Eq.1).And.(kont.Ne.0)) Then
    kont = 0
    tadpoles_2L = 0._dp
   End If

   If (kont.Ne.0) Then
    Iname = Iname - 1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If
  !----------------------
  ! mu-Parameter
  !----------------------
   M2_H_1L = M2_H  - tadpoles_1L + tadpoles_2L

   Abs_Mu2 = (M2_H_1L(2) * sinb2 - M2_H_1L(1) * cosb2 )/ cos2b - 0.5_dp * mZ2_mZ
   If (Abs_Mu2.Le.0) Then
    If (ic.Eq.1) Then
     abs_mu2=1.e4_dp
    Else
     WriteWert = Abs_Mu2
     Call WriteLoopMassesError(-5, "LoopMassesMSSM", kont)
     Abs_Mu2 = Abs(Abs_Mu2)
    End If
   Else If (Abs_Mu2.Ge.1.e20_dp) Then
    WriteWert = Abs_Mu2
    Call WriteLoopMassesError(6, "LoopMassesMSSM", kont)
    Abs_Mu2 = 1.e4_dp
   End If

   Abs_mu = Sqrt( Abs_Mu2 )
   mu = Abs_mu * phase_mu
   B = (M2_H_1L(1) + M2_H_1L(2) + 2._dp *  Abs_Mu2) * tanbQ / (1+tanbQ**2)

  vev_Q = Sqrt(vev2)
  mA2_Q = mP02(2)
  mglu_sign = phase_glu * mglu

  !--------------------
  ! pseudoscalar Higgs
  !--------------------
  p2 = mP02(2)
  b_a = tadpoles_1L(1) * sinb2 + tadpoles_1L(2) * cosb2

  mP0_1L(1) = mZ
  mP02_1L(1) = mZ2
#ifdef GENERATIONMIXING
  Call PiPseudoScalar(p2, g(1), g(2), vevs_DR, mP02, RP0                   &
     & , Y_d, uD_L, uD_R, Y_u, uU_L, uU_R, Y_l, uL_L, uL_R                 &
     & , A_d, A_u, A_l, mu_T                                               & 
     & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSlepton2, RSlepton     &
     & , mSneutrino2, RSneutrino, mSpm2, RSpm, mC, mC2, U, V, mS02, RS0    &
     & , mN, mN2, N, mZ2_mZ, mW2_run, PiA0)
#else
  Call PiPseudoScalar(p2, g(1), g(2), vevs_DR, mP02, RP0, Y_d, Y_u, Y_l     &
     & , A_d, A_u, A_l, mu_T, mUSquark2, RUSquark, mDSquark2, RDSquark      & 
     & , mSlepton2, RSlepton, mSneutrino2, mSpm2, RSpm, mC, mC2, U, V, mS02 &
     & , RS0, mN, mN2, N, mZ2_mZ, mW2_run, PiA0)
#endif

  If (Only_1loop_Higgsmass) Then
   Pi2A0 = 0._dp
  Else
   Call PiPseudoScalar2(g(3), mglu_sign, mP02(2), vevs_DR, M2Din, M2Uin, M2Qin &
      & , M2Ein, M2Lin, A_d(3,3), A_u(3,3), A_l(3,3), Y_d(3,3), Y_u(3,3)       &
      & , Y_l(3,3), mu_T, Pi2A0, kont)
  End If

  If ((ic.Eq.1).And.(kont.Ne.0)) Then
   kont = 0
   Pi2A0 = 0._dp
  End If
  
  mP02_1L(2) = (M2_H_1L(2) - M2_H_1L(1)) / cos2b - mZ2_mZ - Real(PiA0,dp) &
      & + Pi2A0 + b_a
  !--------------------------------------------
  ! iteration using on-shell mass for p^2
  !--------------------------------------------
  i_count = 0
  Do
   i_count = i_count + 1
   p2 = mP02_1L(2)

#ifdef GENERATIONMIXING
   Call PiPseudoScalar(p2, g(1), g(2), vevs_DR, mP02, RP0, Y_d, uD_L, uD_R &
     & , Y_u, uU_L, uU_R, Y_l, uL_L, uL_R, A_d, A_u, A_l, mu_T             &
     & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSlepton2, RSlepton     &
     & , mSneutrino2, RSneutrino, mSpm2, RSpm, mC, mC2, U, V, mS02, RS0    &
     & , mN, mN2, N, mZ2_mZ, mW2_run, PiA0)
#else
   Call PiPseudoScalar(p2, g(1), g(2), vevs_DR, mP02, RP0, Y_d, Y_u, Y_l    &
     & , A_d, A_u, A_l, mu_T , mUSquark2, RUSquark, mDSquark2, RDSquark     & 
     & , mSlepton2, RSlepton, mSneutrino2, mSpm2, RSpm, mC, mC2, U, V, mS02 &
     & , RS0, mN, mN2, N, mZ2_mZ, mW2_run, PiA0)
#endif

    mP02_1L(2) = (M2_H_1L(2) - M2_H_1L(1)) / cos2b - mZ2_mZ - Real(PiA0,dp) &
      & + Pi2A0 + b_a

   If (p2.Ne.0._dp) Then
    wert = Abs(mP02_1L(2) - p2) / p2
   Else
    wert = Abs(mP02_1L(2))
   End If
   If (wert.Lt. delta1) Exit
   If (i_count.Gt.30) Then
    Write(ErrCan,*) "Problem in loop for A0",wert,mP02_1L(2)
    Exit
   End If
  End Do

  If (mP02_1L(2).Le.0._dp) Then
   WriteWert = mP02_1L(2)
   Call WriteLoopMassesError(-7, "LoopMasses", kont)
   Write(ErrCan,*) "M^2_A:",mP02_1L(2)
   mP02_1L(2) = 1.e4_dp
   mP0_1L(2) = 1.e2_dp
  Else
   mP0_1L(2) = Sqrt( mP02_1L(2) )
  End If
  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !----------------------
  ! mW(mW)
  !----------------------
  p2 = mW2
  Call PiWWT1(p2, g(2), sinW2_DR, mS02, RS0, mSpm2, RSpm, vevs_DR          &
         & , mP02, RP0, mSneutrino2, RSneutrino, mSlepton2, RSlepton       &
         & , mUSquark2, RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2 &
         & , CKM_Q , mN, mN2, N, mC, mC2, U, V, mZ2_mZ, mW2_run, dmW2)
  mW2_mW = Real(dmW2+mW2,dp)
  !--------------------
  ! charged Higgs
  !--------------------
  mA2_mA = mP02_1L(2) + Real( PiA0,dp )

  mSpm_1L = 0._dp
  mSpm2_1L = 0._dp
  mSpm_1L(1) = mW
  mSpm2_1L(1) = mW2
  p2 = mSpm2(2)
#ifdef GENERATIONMIXING
  Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm  &
   & , A_d, A_u, A_l, mu_T, uU_L, uU_R, uD_L, uD_R, uL_L, uL_R              &
   & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino    &
   & , mSlepton2, RSlepton, mP02, RP0, mS02, RS0, mC, mC2, U, V, mN, mN2, N &
   & , mZ2_mZ, mW2_run, PiSpm)
#else
  Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm   &
   & , A_d, A_u, A_l, mu_T, mUSquark2, RUSquark, mDSquark2, RDSquark         &
   & , mSneutrino2, mSlepton2, RSlepton, mP02, RP0, mS02, RS0, mC, mC2, U, V &
   & , mN, mN2, N, mZ2_mZ, mW2_run, PiSpm)
#endif

  mSpm2_1L(2) = mW2_mW + mA2_mA - Real( PiSpm,dp )
  !--------------------------------------------
  !  iteration using on-shell mass for p^2
  !--------------------------------------------
  i_count = 0
  Do 
   i_count = i_count + 1
   p2 = mSpm2_1L(2)
#ifdef GENERATIONMIXING
   Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm &
   & , A_d, A_u, A_l, mu_T, uU_L, uU_R, uD_L, uD_R, uL_L, uL_R              &
   & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino    &
   & , mSlepton2, RSlepton, mP02, RP0, mS02, RS0, mC, mC2, U, V, mN, mN2, N &
   & , mZ2_mZ, mW2_run, PiSpm)
#else
   Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm  &
   & , A_d, A_u, A_l, mu_T, mUSquark2, RUSquark, mDSquark2, RDSquark         &
   & , mSneutrino2, mSlepton2, RSlepton, mP02, RP0, mS02, RS0, mC, mC2, U, V &
   & , mN, mN2, N, mZ2_mZ, mW2_run, PiSpm)
#endif

   mSpm2_1L(2) = mW2_mW + mA2_mA - Real( PiSpm,dp )

   If (p2.Ne.0._dp) Then
    wert = Abs(mSpm2_1L(2) - p2) / p2
   Else
    wert = Abs(mSpm2_1L(2))
   End If
   If (wert.Lt.delta1) Exit
   If (i_count.Gt.30) Then
    Write(ErrCan,*) "Problem in loop for H+",wert,mSpm2_1L(2)
    Exit
   End If
  End Do


  If (mSpm2_1L(2).Le.0._dp) Then
   WriteWert = mSpm2_1L(2)
   Call WriteLoopMassesError(-8, "LoopMasses", kont)
   mSpm2_1L(2) = 1.e4_dp
   mSpm_1L(2) = 1.e2_dp
  Else
   mSpm_1L(2) = Sqrt( mSpm2_1L(2) )
  End If
  mA2_mA = mP02_1L(2) + Real( PiA0,dp ) 
  !---------------
  ! neutral Higgs
  !---------------
#ifdef GENERATIONMIXING
  Call ScalarMass_Loop_MSSM(mS02, mS02, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ    &
       & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u               &
       & , uD_L, uD_R, uL_L, uL_R, uU_L, uU_R                                  &
       & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino   &
       & , mSlepton2, RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02, RP0, mSpm2 &
       & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein, M2Lin  &
       & , mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, ic, kont)
#else
  Call ScalarMass_Loop_MSSM(mS02, mS02, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ    &
       & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u               &
       & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2               &
       & , mSlepton2, RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02, RP0, mSpm2 &
       & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein, M2Lin  &
       & , mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, ic, kont)
#endif
  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !--------------------------------------------
  ! iteration using on-shell mass for p^2
  !--------------------------------------------
  i_count = 0
  Do 
   i_count = i_count + 1
   comp2 = mS02_1L
#ifdef GENERATIONMIXING
   Call ScalarMass_Loop_MSSM(mS02, mS02_1L, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ&
       & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u               &
       & , uD_L, uD_R, uL_L, uL_R, uU_L, uU_R                                  &
       & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino   &
       & , mSlepton2, RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02, RP0, mSpm2 &
       & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein, M2Lin  &
       & , mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, ic, kont)
#else
   Call ScalarMass_Loop_MSSM(mS02, mS02_1L, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ&
       & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u               &
       & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2               &
       & , mSlepton2, RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02, RP0, mSpm2 &
       & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein, M2Lin  &
       & , mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, ic, kont)
#endif
   If (kont.Ne.0) Then
    Iname = Iname -1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If

   Do i2=1,2
    If (comp2(i2).Ne.0) Then
     comp2(i2) = Abs(comp2(i2) -  mS02_1L(i2)) / comp2(i2)
    Else
     comp2(i2) = mS02_1L(i2)
    End If
   End Do
   If (Maxval(comp2).Lt.delta1) Exit
   If (i_count.Gt.30) Then
    Write(ErrCan,*) "Problem in loop for h0/H0",comp2,mS02_1L
    Exit
   End If

  End Do

  !-----------
  ! Charginos 
  !-----------
#ifdef GENERATIONMIXING
  Call CharginoMass_Loop(mC, mC2, U, V, g(1), g(2), Y_D, Y_L, Y_U, M2, mu &
      & , vevs_DR, mN, mN2, N, mSpm2, RSpm, mS02, RS0, mP02, RP0          &
      & , mSneutrino2, RSneutrino, uD_L, uD_R, uU_L, uU_R, uL_L, uL_R     &
      & , mSlepton2, RSlepton, mDSquark2, RDSquark, mUSquark2, RUSquark   &
      & , mZ2_mZ, mW2_run, delta1, mC_1L, mC2_1L, U_1L, V_1L, kont)
#else
  Call CharginoMass_Loop(mC, mC2, U, V, g(1), g(2), Y_D, Y_L, Y_U, M2, mu     &
      & , vevs_DR, mN, mN2, N, mSpm2, RSpm, mS02, RS0, mP02, RP0, mSneutrino2 &
      & , mSlepton2, RSlepton, mDSquark2, RDSquark, mUSquark2, RUSquark       &
      & , mZ2_mZ, mW2_run, delta1, mC_1L, mC2_1L, U_1L, V_1L, kont)
#endif
  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !------------- 
  ! Neutralinos
  !-------------
#ifdef GENERATIONMIXING
  Call NeutralinoMass_Loop(g(1), g(2), Y_d, Y_l, Y_u, vevs_DR, M1, M2, mu   &
    & , mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark &
    & , mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02, RS0             &
    & , mP02, RP0, mSpm2, RSpm, uL_L, uL_R, uD_L, uD_R, uU_L, uU_R          &
    & , mZ2_mZ, mW2_run, delta1, mN_1L, mN2_1L, N_1L, kont)
#else
  Call NeutralinoMass_Loop(g(1), g(2), Y_d, Y_l, Y_u, vevs_DR, M1, M2, mu   &
    & , mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark &
    & , mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02, RP0, mSpm2, RSpm &
    & , mZ2_mZ, mW2_run, delta1, mN_1L, mN2_1L, N_1L, kont)
#endif
  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
 !-------------
 ! Gluino
 !-------------
#ifdef GENERATIONMIXING
  Call Sigma_Gluino(mglu**2, mudim, g(3), mglu,phase_glu, mUSquark2, RUSquark &
                  &, uU_L, uU_R, mDSquark2, RDSquark, uD_L, uD_R, dmglu)
#else
  Call Sigma_Gluino(mglu**2, mudim, g(3), mglu,phase_glu, mUSquark2, RUSquark &
                  &, mDSquark2, RDSquark, dmglu)
#endif
  mglu_1L = Abs( mglu - dmglu )
  phase_glu = phase_glu * ( mglu - dmglu ) / mglu_1L

  ! recalculation with improved mass
  i_count = 0
  Do 
   i_count = i_count + 1
   comp2(1) = mglu_1L**2
#ifdef GENERATIONMIXING
   Call Sigma_Gluino(comp2(1), mudim, g(3), mglu,phase_glu, mUSquark2 &
            &, RUSquark, uU_L, uU_R, mDSquark2, RDSquark, uD_L, uD_R, dmglu)
#else
   Call Sigma_Gluino(comp2(1), mudim, g(3), mglu,phase_glu, mUSquark2   &
                  &, RUSquark, mDSquark2, RDSquark, dmglu)
#endif
   mglu_1L = Abs( mglu - dmglu )
   phase_glu = phase_glu * ( mglu - dmglu ) / mglu_1L
   If (comp2(1).Ne. 0._dp) Then
    comp2(1) = Abs( Sqrt(comp2(1)) -  mglu_1L) / Sqrt(comp2(1))
   Else
    comp2(1) = mglu_1L 
   End If 
   If (comp2(1).Lt. delta1) Exit

   If (i_count.Gt.30) Then
    Write(ErrCan,*) "Problem in loop for ~g",comp2(1),mglu_1L
    Exit
   End If

  End Do
  !--------------------------
  ! Up Squarks
  !--------------------------
  T3 = 0.5_dp
  eq = 2._dp / 3._dp
#ifdef GENERATIONMIXING
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_U, A_u, A_d   &
       & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
       & , RDSquark, mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02        &
       & , RS0, mP02, RP0, mSpm2, RSpm, mglu, phase_glu, uD_L, uD_R, uU_L      &
       & , uU_R, mZ2_mZ, mW2_run, delta1                                       &
       & , mUSquark_1L, mUSquark2_1L, RUSquark_1L, kont)
#else
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_U, A_u, A_d   &
       & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
       & , RDSquark, mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02, RP0    &
       & , mSpm2, RSpm, mglu, phase_glu , mZ2_mZ, mW2_run, delta1              &
       & , mUSquark_1L, mUSquark2_1L, RUSquark_1L, kont)
#endif
  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !--------------------------
  ! Down Squarks
  !--------------------------
  T3 = -0.5_dp
  eq = -1._dp / 3._dp
#ifdef GENERATIONMIXING
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_D, A_d, A_u   &
       & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
       & , RDSquark, mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02        &
       & , RS0, mP02, RP0, mSpm2, RSpm, mglu, phase_glu, uD_L, uD_R, uU_L      &
       & , uU_R, mZ2_mZ, mW2_run, delta1                                       &
       & , mDSquark_1L, mDSquark2_1L, RDSquark_1L, kont)
#else
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_D, A_d, A_u   &
       & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
       & , RDSquark, mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02, RP0    &
       & , mSpm2, RSpm, mglu, phase_glu , mZ2_mZ, mW2_run, delta1              &
       & , mDSquark_1L, mDSquark2_1L, RDSquark_1L, kont)
#endif
  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !--------------------------
  ! Sleptons
  !--------------------------
#ifdef GENERATIONMIXING
  Call SleptonMass_1L(g(1), g(2), Y_l, Y_d, vevs_DR, M2_L, M2_E, A_l, mu_T, mu &
        & , mN, mN2, N, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark  &
        & , mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02, RS0, mP02    &
        & , RP0, mSpm2, RSpm, uL_L, uL_R, mZ2_mZ, mW2_run, delta1            &
        & , mSlepton_1L, mSlepton2_1L , RSlepton_1L, kont)
#else
  Call SleptonMass_1L(g(1), g(2), Y_l, Y_d, vevs_DR, M2_L, M2_E, A_l, mu_T, mu     &
        & , mN, mN2, N, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark  &
        & , mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02, RP0           &
        & , mSpm2, RSpm, mZ2_mZ, mW2_run, delta1                             &
        & , mSlepton_1L, mSlepton2_1L, RSlepton_1L, kont)
#endif
  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !--------------------------
  ! Sneutrinos
  !--------------------------
#ifdef GENERATIONMIXING
  Call SneutrinoMass_1L(g(1), g(2), Y_l, vevs_DR, M2_L, A_l, mu_T, mN2, N, mC &
        & , mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark, mSlepton2    &
        & , RSlepton, mSneutrino2, RSneutrino, mS02, RS0, mP02, RP0, mSpm2    &
        & , RSpm, uL_L, uL_R, mZ2_mZ, mW2_run, delta1                         &
        & , mSneutrino_1L, mSneutrino2_1L, RSneutrino_1L, kont)
#else
  Call SneutrinoMass_1L(g(1), g(2), Y_l, vevs_DR, M2_L, A_l, mu_T, mN2, N, mC &
       & , mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark, mSlepton2     &
       & , RSlepton, mSneutrino2, mS02, RS0, mP02, RP0, mSpm2, RSpm           &
       & , mZ2_mZ, mW2_run, delta1                                            &
       & , mSneutrino_1L, mSneutrino2_1L, RSneutrino_1L, kont)
#endif
  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !------------------------------------------------------------
  ! shifting masses and mixing angles from tree-level to 1-loop
  !------------------------------------------------------------
  mP0 = mP0_1L
  mP02 = mP02_1L

  mS0 = mS0_1L
  mS02 = mS02_1L
  RS0 = RS0_1L

  mSpm = mSpm_1L
  mSpm2 = mSpm2_1L

  mC = mC_1L
  mC2 = mC2_1L

  U = U_1L
  V = V_1L

  mN = mN_1L
  mN2 = mN2_1L
  N = N_1L

  mglu = mglu_1L

  mSneutrino = mSneutrino_1L
  mSneutrino2 = mSneutrino2_1L
  mUSquark = mUSquark_1L
  mUSquark2 = mUSquark2_1L
  mDSquark = mDSquark_1L
  mDSquark2 = mDSquark2_1L
  mSlepton = mSlepton_1L
  mSlepton2 = mSlepton2_1L

  RUsquark = RUsquark_1L
  RDSquark = RDSquark_1L
  RSlepton = RSlepton_1L
  RSneutrino = RSneutrino_1L

 !---------------------------
 ! resetting RP0 and RSpm
 !---------------------------
 If (mudim.Ne.mZ2) Then
   cosb2 = 1._dp / (1._dp + tanbQ**2)
   sinb2 = 1._dp - cosb2
   RP0(1,1) = - Sqrt(cosb2)
   RP0(1,2) = Sqrt(sinb2)
   RP0(2,2) = - RP0(1,1)
   RP0(2,1) = RP0(1,2)
   RSpm(1,1) = - Sqrt(cosb2)
   RSpm(1,2) = Sqrt(sinb2)
   RSpm(2,2) = - RP0(1,1)
   RSpm(2,1) = RP0(1,2)
  End If
  B_loop = B
  mu_loop = mu
  B = B_T
  mu = mu_T
  !------------------------------------------------
  ! replacing running fermion masses by pole masses
  !------------------------------------------------
  mf_d = mf_d_save
  mf_l = mf_l_save
  mf_u = mf_u_save
  mf_d2 = mf_d**2
  mf_l2 = mf_l**2
  mf_u2 = mf_u**2

  Iname = Iname -1


 End Subroutine LoopMassesMSSM


 Subroutine LoopMassesMSSM_2(delta, tanb_mZ, tanb_Q, g, Y_l, Y_d, Y_u, Mi     &
   & , A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, mu                        &
   & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP02, RP0                   &
   & , mSpm, mSpm2, RSpm, mDsquark, mDsquark2, RDsquark, mUsquark, mUsquark2  &
   & , RUsquark, mSlepton, mSlepton2, RSlepton, mSneutrino, mSneutrino2       &
   & , RSneutrino, mGlu, phase_glu, M2_H, B, kont)
!-----------------------------------------------------------------------
! In this Subroutine the 1-loop masses within the MSSM are calculated
! It is assumed that the parameters are provided at the same renormalisation
! scale mudim
! The corresponding Subroutine from LoopMasses.F is taken as basis
! written by Werner Porod
! 08.10.01: - defining Interface
!           - first Call of TreeMassesMSSM and AllCouplings implemented
!             here it is still assumed that there is no Generation mixing
!           - implementing the running of Tan(beta)
!           - implementing PiZZT1
! 09.10.01: - implementing One_Loop_Tadpoles_MSSM
!           - implementing PiPseudoscalar
! 10.10.01: - implementing PiChargedScalar
! 10.11.03: taking routine LoopMassesMSSM as basis to get a routine
!           where mu and m_A are input and not an output.
! 15.04.03: tan(beta) given at Q
!-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3), A_l(3,3)   &
      , A_d(3,3), A_u(3,3), M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_Q(3,3)       &
      , M2_U(3,3), mu, Mi(3)
  Real(dp), Intent(in) :: delta, g(3), mP02(2) 
  Real(dp), Intent(inout) :: tanb_mZ, tanb_Q

  Real(dp), Intent(out) :: mC(2), mN(4), mS0(2), mSpm(2), mglu   &
    , mUsquark(6), mDsquark(6), mSlepton(6), mSneutrino(3)       &
    , mUsquark2(6), mDsquark2(6), mSlepton2(6), mSneutrino2(3)   &
    , mS02(2), mSpm2(2), RP0(2,2), mC2(2), mN2(4), RS0(2,2), M2_H(2)
  Complex(dp), Intent(out) :: U(2,2), V(2,2), N(4,4), RSpm(2,2), phase_Glu  &
        , RDsquark(6,6), RUsquark(6,6), RSlepton(6,6), RSneutrino(3,3), B 
 

  Real(dp) :: vevSM(2), tanbQ, cosb2, cos2b, sinb2, g57(57), g58(58) , Q, dt &
    & , tz, vev2, sinW2_DR, mZ2_mZ, vevs_DR(2), tadpoles_1L(2), p2, b_A      &
    & , mW2_mW, mA2_mA, eq, T3, Q2, Pi2A0, tadpoles_2L(2), mW2_run, mP0_T(2) &
    & , mP02_T(2), mglu_sign, M2_H_1L(2), RP0_T(2,2), comp2(2), wert
  Real(dp) :: M2Ein, M2Lin, M2Din, M2Qin, M2Uin, delta1, Abs_mu, abs_mu2
  Complex(dp) :: dmZ2, PiA0, dmW2, PiSpm, dmglu, M1, M2, mu_T, B_T
  Complex(dp), Dimension(3,3) :: uL_L, uL_R, uD_L, uD_R, uU_L, uU_R, CKM_Q
  Integer :: i1, i_loop, i2, i_count
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'LoopMassesMSSM_2'

  kont = 0
  mudim = GetRenormalizationScale() ! from LoopFunctions
!---------------------------------
! definning model
!---------------------------------
  n_S0 = 2
  n_P0 = 2
  n_Spm = 2
  n_char = 2
  n_neut = 4
  n_Sneut = 3
  n_Slept = 6

  uL_L = id3C
  uL_R = id3C
  uD_L = id3C
  uD_R = id3C
  uU_L = id3C
  uU_R = id3C

  M1 = Mi(1)
  M2 = Mi(2)
  !----------------------------------------------------------------------------
  ! require internally a factor 10 smaller delta to be sure, except of smaller
  ! than 100*Epsilon(1._dp) to avoid numerical problems
  !----------------------------------------------------------------------------
  delta1 = delta
  If (delta1.Gt.Epsilon(1._dp)*100._dp) delta1 = 0.1_dp * delta1
  !---------------------------
  ! running of Tan(beta)
  !---------------------------
  If (tanb_in_at_Q) Then
   Q = Sqrt(mudim)
   Q2 = mudim
   Call CouplingsToG(g,Y_l,Y_d,Y_u,g58(1:57))
   g58(1) = Sqrt(5._dp / 3._dp ) * g58(1)  ! rescaling
   g58(58) = Log(tanb_Q)
   tz = Log(Q/mZ)
   dt = - tz / 50._dp
   Call odeint(g58, 58, tz, 0._dp, delta1, dt, 0._dp, rge58, kont)

   tanb_mZ = Exp(g58(58))
   tanbQ = tanb_Q
  Else If (mudim.Ne.mZ2) Then
   Q = Sqrt(mudim)
   Q2 = mudim
   Call CouplingsToG(g,Y_l,Y_d,Y_u,g57)
   g57(1) = Sqrt(5._dp / 3._dp ) * g57(1)  ! rescaling

   !------------------------------------
   ! calculate first the couplings at mZ
   !------------------------------------
   tz = Log(Q/mZ)
   dt = - tz / 50._dp
   
   Call odeint(g57, 57, tz, 0._dp,  delta1, dt, 0._dp, rge57, kont)
   !---------------------------------------
   ! and now back to Q including Tan(beta)
   !---------------------------------------
   g58(1:57) = g57
   g58(58) = Log( tanb_mZ )
   tz = Log(Q/mZ)
   dt = tz / 50._dp

   Call odeint(g58, 58, 0._dp, tz,  delta1 , dt, 0._dp, rge58, kont)

   tanbQ = Exp( g58(58) )
   tanb_Q = tanbQ
  Else
   tanbQ = tanb_mZ
   tanb_Q = tanbQ
   Q2 = mZ2
  End If

  cosb2 = 1._dp / (1._dp + tanbQ**2)
  sinb2 = 1._dp - cosb2
  cos2b = cosb2 - sinb2

  vev2 =  4._dp * mZ2 / (g(1)**2 + g(2)**2)
  vevSM(1) = Sqrt(vev2 * cosb2 )
  vevSM(2) = tanbQ * vevSM(1)
!------------------------------------------------
! replacing fermion pole masses by running masses
!------------------------------------------------
  mf_d_save = mf_d
  mf_l_save = mf_l
  mf_u_save = mf_u
  If (GenerationMixing) Then
   Call FermionMass(Y_l,vevSM(1),mf_l,uL_L, uL_R,kont)
   Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevSM, mf_d, uD_L, uD_R &
                                     & , mf_u, uU_L, uU_R, CKM_Q)
  Else
   Do i1=1,3
    mf_d(i1) = oosqrt2 * Abs(Y_d(i1,i1)) * vevSM(1)
    mf_l(i1) = oosqrt2 * Abs(Y_l(i1,i1)) * vevSM(1)
    mf_u(i1) = oosqrt2 * Abs(Y_u(i1,i1)) * vevSM(2)
   End Do
   CKM_Q = id3C
  End If
  mf_d2 = mf_d**2
  mf_l2 = mf_l**2
  mf_u2 = mf_u**2
!-------------------------------------------------
! tree level masses and mixings as starting point
! first guess of B and M2_H
!-------------------------------------------------
  B =  mP02(2) * tanbQ / (1._dp + tanbQ**2)
  Call TreeMassesMSSM2(g(1), g(2), vevSM, Mi(1), Mi(2), Mi(3), mu, B, tanbQ   &
        , M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u          &
        , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                                  &
        , mGlu, Phase_Glu, mC, mC2, U, V, mN, mN2, N                          &
        , mSneutrino, mSneutrino2, Rsneutrino, mSlepton, mSlepton2, RSlepton  &
        , mDSquark, mDSquark2, RDSquark, mUSquark, mUSquark2, RUSquark        &
        , mP0_T, mP02_T, RP0_T, mS0, mS02, RS0, mSpm, mSpm2, RSpm, mZ2, mW2   &
        , GenerationMixing, kont, .True., .False.)

  kont = 0
  sinW2_DR = g(1)**2 / (g(1)**2+g(2)**2)

!----------------------
! mZ(mZ)
!----------------------
  Call PiZZT1(mZ2, g(2), sinW2_DR, vevSM, mZ2, mW2, mS02, RS0, mP02_T, RP0_T &
   , mSpm2 , RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton, mUSquark2   &
   , RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2, mC, mC2, U, V       &
   , mN, mN2, N ,dmZ2)
  vev2 =  4._dp * Real(mZ2+dmZ2,dp) / (g(1)**2 + g(2)**2)
  vevs_DR(1) = Sqrt(vev2 / (1._dp+tanbQ**2) )
  vevs_DR(2) = tanbQ * vevs_DR(1)
  mZ2_mZ = Real(dmZ2+mZ2,dp)
  mW2_run = mZ2_mZ * (1._dp - sinW2_DR)

  If (mZ2_mZ.Lt.0._dp) Then
   WriteWert = mZ2_mZ
   Call WriteLoopMassesError(4, "LoopMassesMSSM_2", kont)
   mZ2_mZ = mZ2
  End If
!------------------------------------------------
! replacing fermion pole masses by running masses
!------------------------------------------------
   If (GenerationMixing) Then
    Call FermionMass(Y_l,vevs_DR(1),mf_l,uL_L, uL_R,kont)
    Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevs_DR, mf_d, uD_L, uD_R &
                                  & , mf_u, uU_L, uU_R, CKM_Q)
   Else
    Do i1=1,3
     mf_d(i1) = oosqrt2 * Abs(Y_d(i1,i1)) * vevs_DR(1)
     mf_l(i1) = oosqrt2 * Abs(Y_l(i1,i1)) * vevs_DR(1)
     mf_u(i1) = oosqrt2 * Abs(Y_u(i1,i1)) * vevs_DR(2)
    End Do
    CKM_Q = id3C
   End If
  mf_d2 = mf_d**2
  mf_l2 = mf_l**2
  mf_u2 = mf_u**2

  mglu_sign = mGlu * phase_Glu
  mu_T = mu
  B_T = B

  Do i_loop=1,2
  !-------------------------------------------------------------
  ! tree level masses and mixings using the loop corrected vevs
  !-------------------------------------------------------------

   Call PiZZT1(mZ2, g(2), sinW2_DR, vevs_DR, mZ2_mZ, mW2_run, mS02, RS0, mP02_T&
           , RP0_T, mSpm2, RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton  &
           , mUSquark2 , RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2    &
           , mC, mC2, U, V, mN, mN2, N ,dmZ2)
   vev2 =  4._dp * Real(mZ2+dmZ2,dp) / (g(1)**2 + g(2)**2)
   vevs_DR(1) = Sqrt(vev2 / (1._dp+tanbQ**2) )
   vevs_DR(2) = tanbQ * vevs_DR(1)
   mZ2_mZ = Real(dmZ2+mZ2,dp)
   mW2_run = mZ2_mZ * (1._dp - sinW2_DR)

   vev_Q = Sqrt(vev2) ! for output

   If (mZ2_mZ.Lt.0._dp) Then
    WriteWert = mZ2_mZ
    Call WriteLoopMassesError(-4, "LoopMassesMSSM_2", kont)
    Iname = Iname - 1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If
   !------------------------------------------------
   ! replacing fermion pole masses by running masses
   !------------------------------------------------
   If (GenerationMixing) Then
    Call FermionMass(Y_l,vevs_DR(1),mf_l,uL_L, uL_R,kont)
    Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevs_DR, mf_d, uD_L, uD_R &
                                  & , mf_u, uU_L, uU_R, CKM_Q)
   Else
    Do i1=1,3
     mf_d(i1) = oosqrt2 * Abs(Y_d(i1,i1)) * vevs_DR(1)
     mf_l(i1) = oosqrt2 * Abs(Y_l(i1,i1)) * vevs_DR(1)
     mf_u(i1) = oosqrt2 * Abs(Y_u(i1,i1)) * vevs_DR(2)
    End Do
    CKM_Q = id3C
   End If
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2

   Call TreeMassesMSSM2(g(1), g(2), vevs_DR, Mi(1), Mi(2), Mi(3), mu_T, B_T  &
        , tanbQ, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u  &
        , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                                 &
        , mGlu, Phase_Glu, mC, mC2, U, V, mN, mN2, N                         &
        , mSneutrino, mSneutrino2, Rsneutrino, mSlepton, mSlepton2, RSlepton &
        , mDSquark, mDSquark2, RDSquark, mUSquark, mUSquark2, RUSquark       &
        , mP0_T, mP02_T, RP0_T, mS0, mS02, RS0, mSpm, mSpm2, RSpm, mZ2, mW2  &
        , GenerationMixing, kont, .True., .False.)
   RP0 = RP0_T
   !---------------------------
   ! tadpoles at 1-loop
   !---------------------------
#ifdef GENERATIONMIXING
   Call One_Loop_Tadpoles_MSSM(g(1), g(2), vevs_DR, Y_l, Y_d, Y_u, uU_L, uU_R &
    & , uD_L, uD_R, uL_L, uL_R, mu_T, A_l, A_d, A_u, mSneutrino2, Rsneutrino  &
    & , mSlepton2, Rslepton, mDSquark2, RDSquark, mUSquark2, RUSquark         &
    & , mSpm2, RSpm, mC, mC2, U, V, mP02_T, RP0, mS02, RS0, mN, mN2, N        &
    & , mZ2_mZ, mW2_run, tadpoles_1L)
#else
   Call One_Loop_Tadpoles_MSSM(g(1), g(2), vevs_DR, Y_l, Y_d, Y_u, mu_T, A_l &
    & , A_d, A_u, mSneutrino2, mSlepton2, Rslepton, mDSquark2, RDSquark      &
    & , mUSquark2, RUSquark, mSpm2, RSpm, mC, mC2, U, V, mP02_T, RP0, mS02   &
    & , RS0, mN, mN2, N, mZ2_mZ, mW2_run, tadpoles_1L)
#endif
   !---------------------------
   ! tadpoles at 2-loop
  !---------------------------------------------------------------------
  ! Define first input, numerical problems can occur if left and right
  ! sfermion mass parameters are identical
  !---------------------------------------------------------------------
   M2Din = Real(M2_D(3,3),dp)
   M2Ein = Real(M2_E(3,3),dp)
   M2Lin = Real(M2_L(3,3),dp)
   M2Qin = Real(M2_Q(3,3),dp)
   M2Uin = Real(M2_U(3,3),dp)
   If (M2Ein.Eq.M2Lin) M2Ein = (1._dp + Max(1.0e-8_dp,delta1)) * M2Ein
   If (M2Din.Eq.M2Qin) M2Din = (1._dp + Max(1.0e-8_dp,delta1)) * M2Din
   If (M2Uin.Eq.M2Qin) M2Uin = (1._dp + Max(1.0e-8_dp,delta1)) * M2Uin

   If (Only_1loop_Higgsmass) Then
    tadpoles_2L = 0._dp
   Else
    Call Two_Loop_Tadpoles_MSSM(g(3), mglu_sign, mP02_T(2), vevs_DR  &
       & , M2Din, M2Uin, M2Qin, M2Ein, M2Lin, A_d(3,3), A_u(3,3)     &
       & , A_l(3,3), Y_d(3,3), Y_u(3,3), Y_l(3,3), mu_T, tadpoles_2L, kont)
   End If

   If (kont.Ne.0) Then
    Iname = Iname - 1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If

   !--------------------
   ! pseudoscalar Higgs 
   !--------------------
   b_a = tadpoles_1L(1) * sinb2 + tadpoles_1L(2) * cosb2

   p2 =  mP02(2)

#ifdef GENERATIONMIXING
  Call PiPseudoScalar(p2, g(1), g(2), vevs_DR, mP02_T, RP0              &
     & , Y_d, uD_L, uD_R, Y_u, uU_L, uU_R, Y_l, uL_L, uL_R              &
     & , A_d, A_u, A_l, mu_T                                            & 
     & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSlepton2, RSlepton  &
     & , mSneutrino2, RSneutrino, mSpm2, RSpm, mC, mC2, U, V, mS02, RS0 &
     & , mN, mN2, N, mZ2_mZ, mW2_run, PiA0)
#else
  Call PiPseudoScalar(p2, g(1), g(2), vevs_DR, mP02_T, RP0, Y_d, Y_u, Y_l   &
     & , A_d, A_u, A_l, mu_T , mUSquark2, RUSquark, mDSquark2, RDSquark     & 
     & , mSlepton2, RSlepton, mSneutrino2, mSpm2, RSpm, mC, mC2, U, V, mS02 &
     & , RS0, mN, mN2, N, mZ2_mZ, mW2_run, PiA0)
#endif

  If (Only_1loop_Higgsmass) Then
   Pi2A0 = 0._dp
  Else
   Call PiPseudoScalar2(g(3), mglu_sign, mP02_T(2), vevs_DR, M2Din, M2Uin &
     & , M2Qin, M2Ein, M2Lin, A_d(3,3), A_u(3,3), A_l(3,3), Y_d(3,3)      &
     & , Y_u(3,3), Y_l(3,3), mu_T, Pi2A0, kont)
  End If

   mA2_mA = -Pi2A0+ PiA0 + mP02(2) - b_a 
   M2_H_1L(1) = 0.5_dp * ( mA2_mA - 2._dp * Abs(mu)**2    &
              &          - (mZ2_mZ + mA2_mA) * cos2b )
   M2_H_1L(2) = 0.5_dp * ( mA2_mA - 2._dp * Abs(mu)**2    &
              &          + (mZ2_mZ + mA2_mA) * cos2b )
   B = (M2_H_1L(1) + M2_H_1L(2) + 2._dp *  Abs(Mu)**2) * tanbQ / (1+tanbQ**2)

   M2_H = M2_H_1L + tadpoles_1L - tadpoles_2L

   Abs_Mu2 = (M2_H(2) * sinb2 - M2_H(1) * cosb2) / cos2b  - 0.5_dp * mZ2_mZ
   If (Abs_Mu2.Le.0) Then
    WriteWert = Abs_Mu2
    Call WriteLoopMassesError(1, "LoopMassesMSSM_2", kont)
    Abs_Mu2 = Abs(Abs_Mu2)
   Else If (Abs_Mu2.Ge.1.e20_dp) Then
    WriteWert = Abs_Mu2
    Call WriteLoopMassesError(2, "LoopMassesMSSM_2", kont)
    Abs_Mu2 = 1.e4_dp
   End If

   Abs_Mu = Sqrt( Abs_Mu2 )
   mu_T = Abs_mu * phase_mu
   B_T = (M2_H(1) + M2_H(2) + 2._dp *  Abs_Mu2) * tanbQ / (1+tanbQ**2)
  End Do ! i_loop


  mA2_Q = -Pi2A0+ PiA0 + mP02(2)  ! for output
  !------------------------
  ! mW(mW)
  !------------------------
  p2 = mW2
  Call PiWWT1(p2, g(2), sinW2_DR, mS02, RS0, mSpm2, RSpm, vevs_DR          &
         & , mP02, RP0, mSneutrino2, RSneutrino, mSlepton2, RSlepton       &
         & , mUSquark2, RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2 &
         & , CKM_Q , mN, mN2, N, mC, mC2, U, V, mZ2_mZ, mW2_run, dmW2)
  mW2_mW = Real(dmW2+mW2,dp)
  !--------------------
  ! charged Higgs
  !--------------------
  mA2_mA = mP02(2) + Real( PiA0,dp ) 

  mSpm_1L = 0._dp
  mSpm2_1L = 0._dp
  mSpm_1L(1) = mW
  mSpm2_1L(1) = mW2
  p2 = mSpm2(2)

#ifdef GENERATIONMIXING
  Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm    &
   & , A_d, A_u, A_l, mu_T, uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                &
   & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino      &
   & , mSlepton2, RSlepton, mP02_T, RP0, mS02, RS0, mC, mC2, U, V, mN, mN2, N &
   & , mZ2_mZ, mW2_run, PiSpm)
#else
  Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm &
   & , A_d, A_u, A_l, mu_T, mUSquark2, RUSquark, mDSquark2, RDSquark       &
   & , mSneutrino2, mSlepton2, RSlepton, mP02_T, RP0, mS02, RS0, mC, mC2   &
   & , U, V, mN, mN2, N, mZ2_mZ, mW2_run, PiSpm)
#endif


  mSpm2_1L(2) = mW2_mW + mA2_mA - Real( PiSpm,dp )
  

  If (mSpm2_1L(2).Le.0._dp) Then
   Call WriteLoopMassesError(-8, "LoopMasses_2", kont)
   mSpm2_1L(2) = 1.e4_dp
   mSpm_1L(2) = 1.e2_dp
  Else
   mSpm_1L(2) = Sqrt( mSpm2_1L(2) )
  End If

  !--------------------------------------------
  !  iteration using on-shell mass for p^2
  !--------------------------------------------
  i_count = 0
  Do 
   i_count = i_count + 1
   p2 = mSpm2_1L(2)
#ifdef GENERATIONMIXING
   Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm   &
   & , A_d, A_u, A_l, mu_T, uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                &
   & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino      &
   & , mSlepton2, RSlepton, mP02_T, RP0, mS02, RS0, mC, mC2, U, V, mN, mN2, N &
   & , mZ2_mZ, mW2_run, PiSpm)
#else
   Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm &
   & , A_d, A_u, A_l, mu_T, mUSquark2, RUSquark, mDSquark2, RDSquark        &
   & , mSneutrino2, mSlepton2, RSlepton, mP02_T, RP0, mS02, RS0, mC, mC2    &
   & , U, V, mN, mN2, N, mZ2_mZ, mW2_run, PiSpm)
#endif

   mSpm2_1L(2) = mW2_mW + mA2_mA - Real( PiSpm,dp )
   If (p2.Ne.0._dp) Then
    wert = Abs(mSpm2_1L(2) - p2) / p2
   Else
    wert = Abs(mSpm2_1L(2))
   End If
   If (wert.Lt. delta1) Exit
   If (i_count.Gt.30) Then
    Write(ErrCan,*) "Problem in loop for H+",wert,mSpm2_1L(2)
    Exit
   End If
  End Do

  If (mSpm2_1L(2).Le.0._dp) Then
   Call WriteLoopMassesError(-8, "LoopMasses_2", kont)
   mSpm2_1L(2) = 1.e4_dp
   mSpm_1L(2) = 1.e2_dp
  Else
   mSpm_1L(2) = Sqrt( mSpm2_1L(2) )
  End If
  !---------------
  ! neutral Higgs
  !---------------
  mA2_mA = mP02(2) + Real( PiA0,dp ) 

#ifdef GENERATIONMIXING
  Call ScalarMass_Loop_MSSM(mS02, mS02, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ      &
       & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u                 &
       & , uD_L, uD_R, uL_L, uL_R, uU_L, uU_R                                    &
       & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino     &
       & , mSlepton2, RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02_T, RP0, mSpm2 &
       & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein, M2Lin    &
       & , mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, 2, kont)
#else
  Call ScalarMass_Loop_MSSM(mS02, mS02, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ      &
       & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u                 &
       & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2                 &
       & , mSlepton2, RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02_T, RP0, mSpm2 &
       & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein, M2Lin    &
       & , mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, 2, kont)
#endif
  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !--------------------------------------------
  ! iteration using on-shell mass for p^2
  !--------------------------------------------
  i_count = 0
  Do 
   i_count = i_count + 1
   comp2 = mS02_1L
#ifdef GENERATIONMIXING
   Call ScalarMass_Loop_MSSM(mS02, mS02_1L, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ&
     & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u                 &
     & , uD_L, uD_R, uL_L, uL_R, uU_L, uU_R                                    &
     & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino     &
     & , mSlepton2, RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02_T, RP0, mSpm2 &
     & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein, M2Lin    &
     & , mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, 2, kont)
#else
   Call ScalarMass_Loop_MSSM(mS02, mS02_1L, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ&
       & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u               &
       & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, mSlepton2    &
       & , RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02_T, RP0, mSpm2          &
       & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein, M2Lin  &
       & , mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, 2, kont)
#endif

   If (kont.Ne.0) Then
    Iname = Iname -1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If

   Do i2=1,2
    If (comp2(i2).Ne.0) Then
     comp2(i2) = Abs(comp2(i2) -  mS02_1L(i2)) / comp2(i2)
    Else
     comp2(i2) = mS02_1L(i2)
    End If
   End Do
   If (Maxval(comp2).Lt.delta1) Exit
   If (i_count.Gt.30) Then
    Write(ErrCan,*) "Problem in loop for h0,H0",comp2,mS02_1L
    Exit
   End If

  End Do

  !-----------
  ! Charginos 
  !-----------
#ifdef GENERATIONMIXING
  Call CharginoMass_Loop(mC, mC2, U, V, g(1), g(2), Y_D, Y_L, Y_U, M2, mu &
      & , vevs_DR, mN, mN2, N, mSpm2, RSpm, mS02, RS0, mP02_T, RP0        &
      & , mSneutrino2, RSneutrino, uD_L, uD_R, uU_L, uU_R, uL_L, uL_R     &
      & , mSlepton2, RSlepton, mDSquark2, RDSquark, mUSquark2, RUSquark   &
      & , mZ2_mZ, mW2_run, delta1, mC_1L, mC2_1L, U_1L, V_1L, kont)
#else
  Call CharginoMass_Loop(mC, mC2, U, V, g(1), g(2), Y_D, Y_L, Y_U, M2, mu   &
      & , vevs_DR, mN, mN2, N, mSpm2, RSpm, mS02, RS0, mP02_T, RP0          &
      & , mSneutrino2, mSlepton2, RSlepton, mDSquark2, RDSquark, mUSquark2  &
      & , RUSquark, mZ2_mZ, mW2_run, delta1, mC_1L, mC2_1L, U_1L, V_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !------------- 
  ! Neutralinos
  !-------------
#ifdef GENERATIONMIXING
  Call NeutralinoMass_Loop(g(1), g(2), Y_d, Y_l, Y_u, vevs_DR, M1, M2, mu    &
    & , mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark  &
    & , mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02, RS0              &
    & , mP02_T, RP0, mSpm2, RSpm, uL_L, uL_R, uD_L, uD_R, uU_L, uU_R         &
    & , mZ2_mZ, mW2_run, delta1, mN_1L, mN2_1L, N_1L, kont)
#else
  Call NeutralinoMass_Loop(g(1), g(2), Y_d, Y_l, Y_u, vevs_DR, M1, M2, mu      &
    & , mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark    &
    & , mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02_T, RP0, mSpm2, RSpm  &
    & , mZ2_mZ, mW2_run, delta1, mN_1L, mN2_1L, N_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !----------------------
  ! Gluino
  !----------------------
#ifdef GENERATIONMIXING
  Call Sigma_Gluino(mglu**2, mudim, g(3), mglu,phase_glu, mUSquark2, RUSquark &
                  &, uU_L, uU_R, mDSquark2, RDSquark, uD_L, uD_R, dmglu)
#else
  Call Sigma_Gluino(mglu**2, mudim, g(3), mglu,phase_glu, mUSquark2, RUSquark &
                  &, mDSquark2, RDSquark, dmglu)
#endif
  mglu_1L = Abs( mglu - dmglu )
  phase_glu = phase_glu * ( mglu - dmglu ) / mglu_1L

  ! recalculation with improved mass
  i_count = 0
  Do 
   i_count = i_count + 1
   comp2(1) = mglu_1L**2
#ifdef GENERATIONMIXING
   Call Sigma_Gluino(comp2(1), mudim, g(3), mglu,phase_glu, mUSquark2 &
              &, RUSquark, uU_L, uU_R, mDSquark2, RDSquark, uD_L, uD_R, dmglu)
#else
   Call Sigma_Gluino(comp2(1), mudim, g(3), mglu,phase_glu, mUSquark2   &
                  &, RUSquark, mDSquark2, RDSquark, dmglu)
#endif
   mglu_1L = Abs( mglu - dmglu )
   phase_glu = phase_glu * ( mglu - dmglu ) / mglu_1L
   If (comp2(1).Ne. 0._dp) Then
    comp2(1) = Abs( Sqrt(comp2(1)) -  mglu_1L) / Sqrt(comp2(1))
   Else
    comp2(1) = mglu_1L 
   End If 
   If (comp2(1).Lt. delta1) Exit

   If (i_count.Gt.30) Then
    Write(ErrCan,*) "Problem in loop for ~g",comp2(1),mglu_1L
    Exit
   End If

  End Do

  !--------------------------
  ! Up Squarks
  !--------------------------
  T3 = 0.5_dp
  eq = 2._dp / 3._dp

#ifdef GENERATIONMIXING
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_U, A_u, A_d   &
       & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
       & , RDSquark, mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02        &
       & , RS0, mP02_T, RP0, mSpm2, RSpm, mglu, phase_glu, uD_L, uD_R, uU_L    &
       & , uU_R, mZ2_mZ, mW2_run, delta1                                       &
       & , mUSquark_1L, mUSquark2_1L, RUSquark_1L, kont)
#else
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_U, A_u, A_d    &
       & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
       & , RDSquark, mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02_t, RP0  &
       & , mSpm2, RSpm, mglu, phase_glu , mZ2_mZ, mW2_run, delta1              &
       & , mUSquark_1L, mUSquark2_1L, RUSquark_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !--------------------------
  ! Down Squarks
  !--------------------------
  T3 = -0.5_dp
  eq = -1._dp / 3._dp

#ifdef GENERATIONMIXING
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_D, A_d, A_u   &
       & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
       & , RDSquark, mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02        &
       & , RS0, mP02_T, RP0, mSpm2, RSpm, mglu, phase_glu, uD_L, uD_R, uU_L    &
       & , uU_R, mZ2_mZ, mW2_run, delta1                                       &
       & , mDSquark_1L, mDSquark2_1L, RDSquark_1L, kont)
#else
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_D, A_d, A_u   &
       & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
       & , RDSquark, mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02_T, RP0  &
       & , mSpm2, RSpm, mglu, phase_glu , mZ2_mZ, mW2_run, delta1              &
       & , mDSquark_1L, mDSquark2_1L, RDSquark_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !--------------------------
  ! Sleptons
  !--------------------------
#ifdef GENERATIONMIXING
  Call SleptonMass_1L(g(1), g(2), Y_l, Y_d, vevs_DR, M2_L, M2_E, A_l, mu_T, mu &
        & , mN, mN2, N, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark    &
        & , mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02, RS0, mP02_T    &
        & , RP0, mSpm2, RSpm, uL_L, uL_R, mZ2_mZ, mW2_run, delta1              &
        & , mSlepton_1L, mSlepton2_1L , RSlepton_1L, kont)
#else
  Call SleptonMass_1L(g(1), g(2), Y_l, Y_d, vevs_DR, M2_L, M2_E, A_l, mu_T, mu &
        & , mN, mN2, N, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark    &
        & , mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02_T, RP0           &
        & , mSpm2, RSpm, mZ2_mZ, mW2_run, delta1                               &
        & , mSlepton_1L, mSlepton2_1L, RSlepton_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !--------------------------
  ! Sneutrinos
  !--------------------------
#ifdef GENERATIONMIXING
  Call SneutrinoMass_1L(g(1), g(2), Y_l, vevs_DR, M2_L, A_l, mu, mN2, N, mC  &
        & , mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark, mSlepton2   &
        & , RSlepton, mSneutrino2, RSneutrino, mS02, RS0, mP02_T, RP0, mSpm2 &
        & , RSpm, uL_L, uL_R, mZ2_mZ, mW2_run, delta1                        &
        & , mSneutrino_1L, mSneutrino2_1L, RSneutrino_1L, kont)
#else
  Call SneutrinoMass_1L(g(1), g(2), Y_l, vevs_DR, M2_L, A_l, mu, mN2, N, mC  &
       & , mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark, mSlepton2    &
       & , RSlepton, mSneutrino2, mS02, RS0, mP02_T, RP0, mSpm2, RSpm        &
       & , mZ2_mZ, mW2_run, delta1                                           &
       & , mSneutrino_1L, mSneutrino2_1L, RSneutrino_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
!------------------------------------------------------------
! shifting masses and mixing angles from tree-level to 1-loop
!------------------------------------------------------------
  mglu = mglu_1L

  mS0 = mS0_1L
  mS02 = mS02_1L
  RS0 = RS0_1L

  mSpm = mSpm_1L
  mSpm2 = mSpm2_1L

  mC = mC_1L
  mC2 = mC2_1L

  U = U_1L
  V = V_1L

  mN = mN_1L
  mN2 = mN2_1L
  N = N_1L

  mSneutrino = mSneutrino_1L
  mSneutrino2 = mSneutrino2_1L
  mUSquark = mUSquark_1L
  mUSquark2 = mUSquark2_1L
  mDSquark = mDSquark_1L
  mDSquark2 = mDSquark2_1L
  mSlepton = mSlepton_1L
  mSlepton2 = mSlepton2_1L
  RUsquark = RUsquark_1L
  RDSquark = RDSquark_1L
  RSlepton = RSlepton_1L
  RSneutrino = RSneutrino_1L
!---------------------------
! resetting RP0 and RSpm
!---------------------------
 If (mudim.Ne.mZ2) Then
   cosb2 = 1._dp / (1._dp + tanbQ**2)
   sinb2 = 1._dp - cosb2
   RP0(1,1) = - Sqrt(cosb2)
   RP0(1,2) = Sqrt(sinb2)
   RP0(2,2) = - RP0(1,1)
   RP0(2,1) = RP0(1,2)
   RSpm = RP0
  End If
!------------------------------------------------
! replacing running fermion masses by pole masses
!------------------------------------------------
  mf_d = mf_d_save
  mf_l = mf_l_save
  mf_u = mf_u_save
  mf_d2 = mf_d**2
  mf_l2 = mf_l**2
  mf_u2 = mf_u**2

  Iname = Iname -1

 End Subroutine LoopMassesMSSM_2


 Subroutine LoopMassesMSSM_3(tanb_mZ, tanb_Q, g, Y_l, Y_d, Y_u, Mi           &
   & , A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, mu, B, delta             &
   & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP0, mP02, RP0             &
   & , mSpm, mSpm2, RSpm, mDsquark, mDsquark2, RDsquark, mUsquark, mUsquark2 &
   & , RUsquark, mSlepton, mSlepton2, RSlepton, mSneutrino, mSneutrino2      &
   & , RSneutrino, mGlu, phase_glu, M2_H, kont)
!-----------------------------------------------------------------------
! In this Subroutine the 1-loop masses within the MSSM are calculated
! It is assumed that the parameters are provided at the same renormalisation
! scale mudim
! The corresponding Subroutine from LoopMasses.F is taken as basis
! written by Werner Porod
! 08.10.01: - defining Interface
!           - first Call of TreeMassesMSSM and AllCouplings implemented
!             here it is still assumed that there is no Generation mixing
!           - implementing the running of Tan(beta)
!           - implementing PiZZT1
! 09.10.01: - implementing One_Loop_Tadpoles_MSSM
!           - implementing PiPseudoscalar
! 10.10.01: - implementing PiChargedScalar
! 10.11.03: taking routine LoopMassesMSSM as basis to get a routine
!           where mu and m_A are input and not an output.
! 15.04.03: tan(beta) given at Q
!-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3), A_l(3,3)   &
      , A_d(3,3), A_u(3,3), M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_Q(3,3)       &
      , M2_U(3,3), mu, B, Mi(3)
  Real(dp), Intent(in) :: g(3), delta
  Real(dp), Intent(inout) :: tanb_mZ, tanb_Q

  Real(dp), Intent(out) :: mC(2), mN(4), mS0(2), mSpm(2), mP0(2), mP02(2)   &
    , mUsquark(6), mDsquark(6), mSlepton(6), mSneutrino(3)       &
    , mUsquark2(6), mDsquark2(6), mSlepton2(6), mSneutrino2(3)   &
    , mS02(2), mSpm2(2), RP0(2,2), mC2(2), mN2(4), RS0(2,2), M2_H(2), mglu
  Complex(dp), Intent(out) :: U(2,2), V(2,2), N(4,4), RSpm(2,2), phase_Glu   &
        , RDsquark(6,6), RUsquark(6,6), RSlepton(6,6), RSneutrino(3,3)
 

  Real(dp) :: vevSM(2), tanbQ, cosb2, cos2b, sinb2, Q, vev2, sinW2_DR    &
    & , mZ2_mZ, vevs_DR(2), tadpoles_1L(2), p2, b_A, mW2_mW, mA2_mA      &
    & , eq, T3, Q2, Pi2A0, tadpoles_2L(2), mW2_run, mglu_sign    &
    & , M2_H_1L(2), comp2(2), wert, g57(57), g58(58), tz, dt
  Real(dp) :: M2Ein, M2Lin, M2Din, M2Qin, M2Uin, delta1, Abs_mu, Abs_mu2
  Complex(dp) :: dmZ2, PiA0, dmW2, PiSpm, M1, M2, dmglu, mu_T, B_T
  Complex(dp), Dimension(3,3) :: uL_L, uL_R, uD_L, uD_R, uU_L, uU_R, CKM_Q
  Integer :: i1, i2, i_loop, i_count
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'LoopMassesMSSM_3'

  kont = 0
  mudim = GetRenormalizationScale() ! from LoopFunctions
  Q = Sqrt(mudim)
  Q2 = mudim
  !---------------------------------
  ! definning model
  !---------------------------------
  n_S0 = 2
  n_P0 = 2
  n_Spm = 2
  n_char = 2
  n_neut = 4
  n_Sneut = 3
  n_Slept = 6

  uL_L = id3C
  uL_R = id3C
  uD_L = id3C
  uD_R = id3C
  uU_L = id3C
  uU_R = id3C

  M1 = Mi(1)
  M2 = Mi(2)
  delta1 = delta
  !---------------------------
  ! running of Tan(beta)
  !---------------------------
  If (tanb_in_at_Q) Then
   Q = Sqrt(mudim)
   Q2 = mudim
   tanbQ = tanb_Q
   Call CouplingsToG(g,Y_l,Y_d,Y_u,g58(1:57))
   g58(1) = Sqrt(5._dp / 3._dp ) * g58(1)  ! rescaling
   g58(58) = Log(tanb_Q)
   tz = Log(Q/mZ)
   dt = - tz / 50._dp
   Call odeint(g58, 58, tz, 0._dp, delta1, dt, 0._dp, rge58, kont)

   tanb_mZ = Exp(g58(58))
  Else If (mudim.Ne.mZ2) Then
   Q = Sqrt(mudim)
   Q2 = mudim
   Call CouplingsToG(g,Y_l,Y_d,Y_u,g57)
   g57(1) = Sqrt(5._dp / 3._dp ) * g57(1)  ! rescaling

   !------------------------------------
   ! calculate first the couplings at mZ
   !------------------------------------
   tz = Log(Q/mZ)
   dt = - tz / 50._dp
   
   Call odeint(g57, 57, tz, 0._dp,  delta1, dt, 0._dp, rge57, kont)
   !---------------------------------------
   ! and now back to Q including Tan(beta)
   !---------------------------------------
   g58(1:57) = g57
   g58(58) = Log( tanb_mZ )
   tz = Log(Q/mZ)
   dt = tz / 50._dp

   Call odeint(g58, 58, 0._dp, tz,  delta1 , dt, 0._dp, rge58, kont)

   tanbQ = Exp( g58(58) )
   tanb_Q = tanbQ
  Else
   tanbQ = tanb_mZ
   tanb_Q = tanbQ
   Q2 = mZ2
  End If

  cosb2 = 1._dp / (1._dp + tanbQ**2)
  sinb2 = 1._dp - cosb2
  cos2b = cosb2 - sinb2

  vev2 =  4._dp * mZ2 / (g(1)**2 + g(2)**2)
  vevSM(1) = Sqrt(vev2 * cosb2 )
  vevSM(2) = tanbQ * vevSM(1)

  !----------------------------------------------------------------------------
  ! require internally a factor 10 smaller delta to be sure, except of smaller
  ! than 100*Epsilon(1._dp) to avoid numerical problems
  !----------------------------------------------------------------------------
  If (delta1.Gt.Epsilon(1._dp)*100._dp) delta1 = 0.1_dp * delta1
  !------------------------------------------------
  ! replacing fermion pole masses by running masses
  !------------------------------------------------
  mf_d_save = mf_d
  mf_l_save = mf_l
  mf_u_save = mf_u
  If (GenerationMixing) Then
   Call FermionMass(Y_l,vevSM(1),mf_l,uL_L, uL_R,kont)
   Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevSM, mf_d, uD_L, uD_R &
                                     & , mf_u, uU_L, uU_R, CKM_Q)
  Else
   Do i1=1,3
    mf_d(i1) = oosqrt2 * Abs(Y_d(i1,i1)) * vevSM(1)
    mf_l(i1) = oosqrt2 * Abs(Y_l(i1,i1)) * vevSM(1)
    mf_u(i1) = oosqrt2 * Abs(Y_u(i1,i1)) * vevSM(2)
   End Do
   CKM_Q = id3C
  End If
  mf_d2 = mf_d**2
  mf_l2 = mf_l**2
  mf_u2 = mf_u**2
  !-------------------------------------------------
  ! tree level masses and mixings as starting point
  !-------------------------------------------------
   Call TreeMassesMSSM2(g(1), g(2), vevSM, Mi(1), Mi(2), Mi(3), mu, B, tanbQ  &
      & , M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u          &
      & , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                                  &
      & , mGlu, Phase_Glu, mC, mC2, U, V, mN, mN2, N                          &
      & , mSneutrino, mSneutrino2, Rsneutrino, mSlepton, mSlepton2, RSlepton  &
      & , mDSquark, mDSquark2, RDSquark, mUSquark, mUSquark2, RUSquark        &
      & , mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm, mZ2, mW2         &
      & , GenerationMixing, kont, .True., .False.)
  kont = 0
  sinW2_DR = g(1)**2 / (g(1)**2+g(2)**2)

!----------------------
! mZ(mZ)
!----------------------
  Call PiZZT1(mZ2, g(2), sinW2_DR, vevSM, mZ2, mW2, mS02, RS0, mP02, RP0      &
      , mSpm2 , RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton, mUSquark2 &
      , RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2, mC, mC2, U, V     &
      , mN, mN2, N ,dmZ2)
  vev2 =  4._dp * Real(mZ2+dmZ2,dp) / (g(1)**2 + g(2)**2)
  vevs_DR(1) = Sqrt(vev2 / (1._dp+tanbQ**2) )
  vevs_DR(2) = tanbQ * vevs_DR(1)
  mZ2_mZ = Real(dmZ2+mZ2,dp)
  mW2_run = mZ2_mZ * (1._dp - sinW2_DR)

  If (mZ2_mZ.Lt.0._dp) Then
   WriteWert = mZ2_mZ
   Call WriteLoopMassesError(4, "LoopMassesMSSM_3", kont)
   mZ2_mZ = mZ2
  End If
  !------------------------------------------------
  ! replacing fermion pole masses by running masses
  !------------------------------------------------
   If (GenerationMixing) Then
    Call FermionMass(Y_l,vevs_DR(1),mf_l,uL_L, uL_R,kont)
    Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevs_DR, mf_d, uD_L, uD_R &
                                  & , mf_u, uU_L, uU_R, CKM_Q)
   Else
    Do i1=1,3
     mf_d(i1) = oosqrt2 * Abs(Y_d(i1,i1)) * vevs_DR(1)
     mf_l(i1) = oosqrt2 * Abs(Y_l(i1,i1)) * vevs_DR(1)
     mf_u(i1) = oosqrt2 * Abs(Y_u(i1,i1)) * vevs_DR(2)
    End Do
    CKM_Q = id3C
   End If
  mf_d2 = mf_d**2
  mf_l2 = mf_l**2
  mf_u2 = mf_u**2

  !-----------------------------------------
  ! intialization for comparision
  !-----------------------------------------
  comp2(1) = mP02(2)
  comp2(2) = mZ2_mZ
  mP02_1L(2) = mP02(2)
  i_count = 1
  mglu_sign = mGlu * phase_Glu
  mu_T = mu
  B_T = B

  Do i_loop=1,31
  !-------------------------------------------------------------
  ! tree level masses and mixings using the loop corrected vevs
  !-------------------------------------------------------------
   Call PiZZT1(mZ2, g(2), sinW2_DR, vevs_DR, mZ2_mZ, mW2_run, mS02, RS0, mP02 &
            , RP0, mSpm2, RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton  &
            , mUSquark2 , RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2  &
            , mC, mC2, U, V, mN, mN2, N ,dmZ2)
   vev2 =  4._dp * Real(mZ2+dmZ2,dp) / (g(1)**2 + g(2)**2)
   vevs_DR(1) = Sqrt(vev2 / (1._dp+tanbQ**2) )
   vevs_DR(2) = tanbQ * vevs_DR(1)
   mZ2_mZ = Real(dmZ2+mZ2,dp)
   mW2_run = mZ2_mZ * (1._dp - sinW2_DR)

   If (mZ2_mZ.Lt.0._dp) Then
    WriteWert = mZ2_mZ
    Call WriteLoopMassesError(-4, "LoopMassesMSSM_3", kont)
    Iname = Iname - 1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If

   Call TreeMassesMSSM2(g(1), g(2), vevs_DR, Mi(1), Mi(2), Mi(3), mu_T, B_T  &
      & , tanbQ, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u  &
      & , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                                 &
      & , mGlu, Phase_Glu, mC, mC2, U, V, mN, mN2, N                         &
      & , mSneutrino, mSneutrino2, Rsneutrino, mSlepton, mSlepton2, RSlepton &
      & , mDSquark, mDSquark2, RDSquark, mUSquark, mUSquark2, RUSquark       &
      & , mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm, mZ2, mW2        &
      & , GenerationMixing, kont, .True., .False.)
   !---------------------------
   ! tadpoles at 1-loop
   !---------------------------
#ifdef GENERATIONMIXING
   Call One_Loop_Tadpoles_MSSM(g(1), g(2), vevs_DR, Y_l, Y_d, Y_u, uU_L, uU_R &
    & , uD_L, uD_R, uL_L, uL_R, mu_T, A_l, A_d, A_u, mSneutrino2, Rsneutrino  &
    & , mSlepton2, Rslepton, mDSquark2, RDSquark, mUSquark2, RUSquark         &
    & , mSpm2, RSpm, mC, mC2, U, V, mP02, RP0, mS02, RS0, mN, mN2, N          &
    & , mZ2_mZ, mW2_run, tadpoles_1L)
#else
   Call One_Loop_Tadpoles_MSSM(g(1), g(2), vevs_DR, Y_l, Y_d, Y_u, mu_T, A_l  &
    & , A_d, A_u, mSneutrino2, mSlepton2, Rslepton, mDSquark2, RDSquark       &
    & , mUSquark2, RUSquark, mSpm2, RSpm, mC, mC2, U, V, mP02, RP0, mS02, RS0 &
    & , mN, mN2, N, mZ2_mZ, mW2_run, tadpoles_1L)
#endif
   !---------------------------
   ! tadpoles at 2-loop
  !---------------------------------------------------------------------
  ! Define first input, numerical problems can occur if left and right
  ! sfermion mass parameters are identical
  !---------------------------------------------------------------------
   M2Din = Real(M2_D(3,3),dp)
   M2Ein = Real(M2_E(3,3),dp)
   M2Lin = Real(M2_L(3,3),dp)
   M2Qin = Real(M2_Q(3,3),dp)
   M2Uin = Real(M2_U(3,3),dp)
   If (M2Ein.Eq.M2Lin) M2Ein = (1._dp + Max(1.0e-8_dp,delta1)) * M2Ein
   If (M2Din.Eq.M2Qin) M2Din = (1._dp + Max(1.0e-8_dp,delta1)) * M2Din
   If (M2Uin.Eq.M2Qin) M2Uin = (1._dp + Max(1.0e-8_dp,delta1)) * M2Uin

   If (Only_1loop_Higgsmass) Then
    tadpoles_2L = 0._dp
   Else
    Call Two_Loop_Tadpoles_MSSM(g(3), mglu_sign, mP02(2), vevs_DR   &
       & , M2Din, M2Uin, M2Qin, M2Ein, M2Lin, A_d(3,3), A_u(3,3)    &
       & , A_l(3,3), Y_d(3,3), Y_u(3,3), Y_l(3,3), mu_T, tadpoles_2L, kont)
   End If

   If (kont.Ne.0) Then
    Iname = Iname - 1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If

   M2_H_1L(1) = - Abs(mu)**2 - 0.5_dp * mZ2_mZ * cos2b + B * tanbQ
   M2_H_1L(2) = - Abs(mu)**2 + 0.5_dp * mZ2_mZ * cos2b + B / tanbQ

   M2_H = M2_H_1L + tadpoles_1L - tadpoles_2L
   Abs_Mu2 = (M2_H(2) * sinb2 - M2_H(1) * cosb2) / cos2b  - 0.5_dp * mZ2
   If (Abs_Mu2.Le.0) Then
    WriteWert = Abs_Mu2
    Call WriteLoopMassesError(1, "LoopMassesMSSM_3", kont)
    Abs_Mu2 = Abs(Abs_Mu2)
   Else If (Abs_Mu2.Ge.1.e20_dp) Then
    WriteWert = Abs_Mu2
    Call WriteLoopMassesError(2, "LoopMassesMSSM_3", kont)
    Abs_Mu2 = 1.e4_dp
   End If

   Abs_Mu = Sqrt( Abs_Mu2 )
   mu_T = Abs_mu * phase_mu
   B_T = (M2_H(1) + M2_H(2) + 2._dp *  Abs(Mu_T)**2) * tanbQ / (1+tanbQ**2)

   !--------------------
   ! pseudoscalar Higgs 
   !--------------------
   b_a = tadpoles_1L(1) * sinb2 + tadpoles_1L(2) * cosb2

   p2 =  mP02_1L(2)
   mP0_1L(1) = mZ
   mP02_1L(1) = mZ2

#ifdef GENERATIONMIXING
   Call PiPseudoScalar(p2, g(1), g(2), vevs_DR, mP02, RP0                &
     & , Y_d, uD_L, uD_R, Y_u, uU_L, uU_R, Y_l, uL_L, uL_R               &
     & , A_d, A_u, A_l, mu_T                                             & 
     & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSlepton2, RSlepton   &
     & , mSneutrino2, RSneutrino, mSpm2, RSpm, mC, mC2, U, V, mS02, RS0  &
     & , mN, mN2, N, mZ2_mZ, mW2_run, PiA0)
#else
   Call PiPseudoScalar(p2, g(1), g(2), vevs_DR, mP02, RP0, Y_d, Y_u, Y_l    &
     & , A_d, A_u, A_l, mu_T, mUSquark2, RUSquark, mDSquark2, RDSquark      & 
     & , mSlepton2, RSlepton, mSneutrino2, mSpm2, RSpm, mC, mC2, U, V, mS02 &
     & , RS0, mN, mN2, N, mZ2_mZ, mW2_run, PiA0)
#endif

  If (Only_1loop_Higgsmass) Then
   Pi2A0 = 0._dp
  Else
   Call PiPseudoScalar2(g(3), mglu_sign, mP02(2), vevs_DR, M2Din, M2Uin, M2Qin &
     & , M2Ein, M2Lin, A_d(3,3), A_u(3,3), A_l(3,3), Y_d(3,3), Y_u(3,3)        &
     & , Y_l(3,3), mu_T, Pi2A0, kont)
  End If

   mP02_1L(2) = (M2_H_1L(2) - M2_H_1L(1)) / cos2b - mZ2_mZ - Real(PiA0,dp) &
              + Pi2A0 + b_a

   If (mP02_1L(2).Le.0._dp) Then
    Call WriteLoopMassesError(-7, "LoopMasses", kont)
    Write(ErrCan,*) "M^2_A:",mP02_1L(2)
    mP02_1L(2) = 1.e4_dp
    mP0_1L(2) = 1.e2_dp
   Else
    mP0_1L(2) = Sqrt( mP02_1L(2) )
   End If

   comp2(1) = Abs(comp2(1) - mP02_1L(2)) / comp2(1)
   comp2(2) = Abs(comp2(2) - mZ2_mZ) / comp2(2)
   If (Maxval(Abs(comp2)).Lt.delta1) Then
    Exit
   Else If (i_count.Gt.30) Then
    Write(ErrCan,*) "Problem for A0, mZ loop in LoopMassesMSSM_3"
    Write(ErrCan,123)   mP02_1L(2),mZ2_mZ,comp2
123 Format(4f18.10)
   Else
    comp2(1) = mP02_1L(2)
    comp2(2) = mZ2_mZ
    i_count = i_count + 1
   End If 
  End Do ! i_loop
  vev_Q = Sqrt(vev2)
  
!------------------------
! mW(mW)
!------------------------
  p2 = mW2
  Call PiWWT1(p2, g(2), sinW2_DR, mS02, RS0, mSpm2, RSpm, vevs_DR          &
         & , mP02, RP0, mSneutrino2, RSneutrino, mSlepton2, RSlepton       &
         & , mUSquark2, RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2 &
         & , CKM_Q , mN, mN2, N, mC, mC2, U, V, mZ2_mZ, mW2_run, dmW2)
  mW2_mW = Real(dmW2+mW2,dp)
  !--------------------
  ! charged Higgs
  !--------------------
  mA2_mA = mP02_1L(2) + Real( PiA0,dp )

  mSpm_1L = 0._dp
  mSpm2_1L = 0._dp
  mSpm_1L(1) = mW
  mSpm2_1L(1) = mW2
  p2 = mSpm2(2)

#ifdef GENERATIONMIXING
  Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm   &
   & , A_d, A_u, A_l, mu_T, uU_L, uU_R, uD_L, uD_R, uL_L, uL_R               &
   & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino     &
   & , mSlepton2, RSlepton, mP02, RP0, mS02, RS0, mC, mC2, U, V, mN, mN2, N  &
   & , mZ2_mZ, mW2_run, PiSpm)
#else
  Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm    &
   & , A_d, A_u, A_l, mu_T, mUSquark2, RUSquark, mDSquark2, RDSquark          &
   & , mSneutrino2, mSlepton2, RSlepton, mP02, RP0, mS02, RS0, mC, mC2, U, V  &
   & , mN, mN2, N, mZ2_mZ, mW2_run, PiSpm)
#endif
  mSpm2_1L(2) = mW2_mW + mA2_mA - Real( PiSpm,dp )
  !--------------------------------------------
  !  iteration using on-shell mass for p^2
  !--------------------------------------------
  i_count = 0
  Do 
   i_count = i_count + 1
   p2 = mSpm2_1L(2)
#ifdef GENERATIONMIXING
   Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm  &
   & , A_d, A_u, A_l, mu_T, uU_L, uU_R, uD_L, uD_R, uL_L, uL_R               &
   & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino     &
   & , mSlepton2, RSlepton, mP02, RP0, mS02, RS0, mC, mC2, U, V, mN, mN2, N  &
   & , mZ2_mZ, mW2_run, PiSpm)
#else
   Call PiChargedScalar(p2, g(1), g(2), Y_d, Y_u, Y_l, vevs_DR, mSpm2, RSpm  &
   & , A_d, A_u, A_l, mu_T, mUSquark2, RUSquark, mDSquark2, RDSquark         &
   & , mSneutrino2, mSlepton2, RSlepton, mP02, RP0, mS02, RS0, mC, mC2, U, V &
   & , mN, mN2, N, mZ2_mZ, mW2_run, PiSpm)
#endif

   mSpm2_1L(2) = mW2_mW + mA2_mA - Real( PiSpm,dp )
   If (p2.Ne.0._dp) Then
    wert = Abs(mSpm2_1L(2) - p2) / p2
   Else
    wert = Abs(mSpm2_1L(2))
   End If
   If (wert.Lt.delta1) Exit
   If (i_count.Gt.30) Then
    Write(ErrCan,*) "Problem in loop for H+",wert,mSpm2_1L(2)
    Exit
   End If
  End Do

  mSpm2_1L(2) = mW2_mW + mA2_mA - Real( PiSpm,dp )
  If (mSpm2_1L(2).Le.0._dp) Then
   Call WriteLoopMassesError(-8, "LoopMasses", kont)
   mSpm2_1L(2) = 1.e4_dp
   mSpm_1L(2) = 1.e2_dp
  Else
   mSpm_1L(2) = Sqrt( mSpm2_1L(2) )
  End If

!---------------
! neutral Higgs
!---------------
  mA2_mA = mP02_1L(2) + Real( PiA0,dp )
#ifdef GENERATIONMIXING
  Call ScalarMass_Loop_MSSM(mS02, mS02, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ   &
      & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u               &
      & , uD_L, uD_R, uL_L, uL_R, uU_L, uU_R                                  &
      & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino   &
      & , mSlepton2, RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02, RP0, mSpm2 &
      & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein         &
      & , M2Lin, mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, 2, kont)
#else
  Call ScalarMass_Loop_MSSM(mS02, mS02, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ  &
       & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u             &
       & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, mSlepton2  &
       & , RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02, RP0, mSpm2          &
       & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein       &
       & , M2Lin, mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, 2, kont)
#endif
  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !--------------------------------------------
  ! iteration using on-shell mass for p^2
  !--------------------------------------------
  i_count = 0
  Do 
   i_count = i_count + 1
   comp2 = mS02_1L
#ifdef GENERATIONMIXING
   Call ScalarMass_Loop_MSSM(mS02, mS02_1L, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ&
       & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u               &
       & , uD_L, uD_R, uL_L, uL_R, uU_L, uU_R                                  &
       & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino   &
       & , mSlepton2, RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02, RP0, mSpm2 &
       & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein         &
       & , M2Lin, mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, 2, kont)
#else
   Call ScalarMass_Loop_MSSM(mS02, mS02_1L, RS0, Q2, tanbQ, tadpoles_1L, mZ2_mZ&
       & , mW2_run, mA2_mA, b_a, g(1), g(2), g(3), Y_l, Y_d, Y_u               &
       & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, mSlepton2    &
       & , RSlepton, A_d, A_u, A_l, mu_T, vevs_DR, mP02, RP0, mSpm2            &
       & , RSpm, mC, mC2, U, V, mN, mN2, N, M2Qin, M2Uin, M2Din, M2Ein, M2Lin  &
       & , mglu_sign, 0, mS0_1L, mS02_1L, RS0_1L, 2, kont)
#endif
   If (kont.Ne.0) Then
    Iname = Iname -1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If

   Do i2=1,2
    If (comp2(i2).Ne.0) Then
     comp2(i2) = Abs(comp2(i2) -  mS02_1L(i2)) / comp2(i2)
    Else
     comp2(i2) = mS02_1L(i2)
    End If
   End Do
   If (Maxval(comp2).Lt.delta1) Exit
   If (i_count.Gt.30) Then
    Write(ErrCan,*) "Problem in loop for h0/H0",comp2,mS02_1L
    Exit
   End If

  End Do

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !-----------
  ! Charginos 
  !-----------
#ifdef GENERATIONMIXING
  Call CharginoMass_Loop(mC, mC2, U, V, g(1), g(2), Y_D, Y_L, Y_U, M2, mu &
      & , vevs_DR, mN, mN2, N, mSpm2, RSpm, mS02, RS0, mP02, RP0          &
      & , mSneutrino2, RSneutrino, uD_L, uD_R, uU_L, uU_R, uL_L, uL_R     &
      & , mSlepton2, RSlepton, mDSquark2, RDSquark, mUSquark2, RUSquark   &
      & , mZ2_mZ, mW2_run, delta1, mC_1L, mC2_1L, U_1L, V_1L, kont)
#else
  Call CharginoMass_Loop(mC, mC2, U, V, g(1), g(2), Y_D, Y_L, Y_U, M2, mu     &
      & , vevs_DR, mN, mN2, N, mSpm2, RSpm, mS02, RS0, mP02, RP0, mSneutrino2 &
      & , mSlepton2, RSlepton, mDSquark2, RDSquark, mUSquark2, RUSquark       &
      & , mZ2_mZ, mW2_run, delta1, mC_1L, mC2_1L, U_1L, V_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !------------- 
  ! Neutralinos
  !-------------
#ifdef GENERATIONMIXING
  Call NeutralinoMass_Loop(g(1), g(2), Y_d, Y_l, Y_u, vevs_DR, M1, M2, mu    &
    & , mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark  &
    & , mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02, RS0              &
    & , mP02, RP0, mSpm2, RSpm, uL_L, uL_R, uD_L, uD_R, uU_L, uU_R           &
    & , mZ2_mZ, mW2_run, delta1, mN_1L, mN2_1L, N_1L, kont)
#else
  Call NeutralinoMass_Loop(g(1), g(2), Y_d, Y_l, Y_u, vevs_DR, M1, M2, mu    &
    & , mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark  &
    & , mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02, RP0, mSpm2, RSpm  &
    & , mZ2_mZ, mW2_run, delta1, mN_1L, mN2_1L, N_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !----------------------
  ! Gluino
  !----------------------
#ifdef GENERATIONMIXING
  Call Sigma_Gluino(mglu**2, mudim, g(3), mglu,phase_glu, mUSquark2, RUSquark &
                  &, uU_L, uU_R, mDSquark2, RDSquark, uD_L, uD_R, dmglu)
#else
  Call Sigma_Gluino(mglu**2, mudim, g(3), mglu,phase_glu, mUSquark2, RUSquark &
                  &, mDSquark2, RDSquark, dmglu)
#endif
  mglu_1L = Abs( mglu - dmglu )
  phase_glu = phase_glu * ( mglu - dmglu ) / mglu_1L

  ! recalculation with improved mass
  i_count = 0
  Do 
   i_count = i_count + 1
   comp2(1) = mglu_1L**2
#ifdef GENERATIONMIXING
   Call Sigma_Gluino(comp2(1), mudim, g(3), mglu,phase_glu, mUSquark2 &
              &, RUSquark, uU_L, uU_R, mDSquark2, RDSquark, uD_L, uD_R, dmglu)
#else
   Call Sigma_Gluino(comp2(1), mudim, g(3), mglu,phase_glu, mUSquark2   &
                  &, RUSquark, mDSquark2, RDSquark, dmglu)
#endif
   mglu_1L = Abs( mglu - dmglu )
   phase_glu = phase_glu * ( mglu - dmglu ) / mglu_1L
   If (comp2(1).Ne. 0._dp) Then
    comp2(1) = Abs( Sqrt(comp2(1)) -  mglu_1L) / Sqrt(comp2(1))
   Else
    comp2(1) = mglu_1L 
   End If 
   If (comp2(1).Lt. delta1) Exit

   If (i_count.Gt.30) Then
    Write(ErrCan,*) "Problem in loop for ~g",comp2(1),mglu_1L
    Exit
   End If

  End Do
  !--------------------------
  ! Up Squarks
  !--------------------------
  T3 = 0.5_dp
  eq = 2._dp / 3._dp

#ifdef GENERATIONMIXING
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_U, A_u, A_d  &
      & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
      & , RDSquark, mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02        &
      & , RS0, mP02, RP0, mSpm2, RSpm, mglu, phase_glu, uD_L, uD_R, uU_L      &
      & , uU_R, mZ2_mZ, mW2_run, delta1                                       &
      & , mUSquark_1L, mUSquark2_1L, RUSquark_1L, kont)
#else
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_U, A_u, A_d  &
      & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
      & , RDSquark, mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02, RP0    &
      & , mSpm2, RSpm, mglu, phase_glu , mZ2_mZ, mW2_run, delta1              &
      & , mUSquark_1L, mUSquark2_1L, RUSquark_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
!--------------------------
! Down Squarks
!--------------------------
  T3 = -0.5_dp
  eq = -1._dp / 3._dp

#ifdef GENERATIONMIXING
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_D, A_d, A_u &
     & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
     & , RDSquark, mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02        &
     & , RS0, mP02, RP0, mSpm2, RSpm, mglu, phase_glu, uD_L, uD_R, uU_L      &
     & , uU_R, mZ2_mZ, mW2_run, delta1                                       &
     & , mDSquark_1L, mDSquark2_1L, RDSquark_1L, kont)
#else
  Call SquarkMass_1L(g, T3, eq, Y_l, Y_d, Y_u, vevs_DR, M2_Q, M2_D, A_d, A_u &
     & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mUSquark2, RUSquark, mDSquark2 &
     & , RDSquark, mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02, RP0    &
     & , mSpm2, RSpm, mglu, phase_glu , mZ2_mZ, mW2_run, delta1              &
     & , mDSquark_1L, mDSquark2_1L, RDSquark_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !--------------------------
  ! Sleptons
  !--------------------------
#ifdef GENERATIONMIXING
  Call SleptonMass_1L(g(1), g(2), Y_l, Y_d, vevs_DR, M2_L, M2_E, A_l, mu_T, mu &
        & , mN, mN2, N, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark    &
        & , mSlepton2, RSlepton, mSneutrino2, RSneutrino, mS02, RS0, mP02      &
        & , RP0, mSpm2, RSpm, uL_L, uL_R, mZ2_mZ, mW2_run, delta1              &
        & , mSlepton_1L, mSlepton2_1L , RSlepton_1L, kont)
#else
  Call SleptonMass_1L(g(1), g(2), Y_l, Y_d, vevs_DR, M2_L, M2_E, A_l, mu_T, mu &
        & , mN, mN2, N, mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark    &
        & , mSlepton2, RSlepton, mSneutrino2, mS02, RS0, mP02, RP0             &
        & , mSpm2, RSpm, mZ2_mZ, mW2_run, delta1                               &
        & , mSlepton_1L, mSlepton2_1L, RSlepton_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
  !--------------------------
  ! Sneutrinos
  !--------------------------
#ifdef GENERATIONMIXING
  Call SneutrinoMass_1L(g(1), g(2), Y_l, vevs_DR, M2_L, A_l, mu, mN2, N, mC  &
        & , mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark, mSlepton2   &
        & , RSlepton, mSneutrino2, RSneutrino, mS02, RS0, mP02, RP0, mSpm2   &
        & , RSpm, uL_L, uL_R, mZ2_mZ, mW2_run, delta1                        &
        & , mSneutrino_1L, mSneutrino2_1L, RSneutrino_1L, kont)
#else
  Call SneutrinoMass_1L(g(1), g(2), Y_l, vevs_DR, M2_L, A_l, mu, mN2, N, mC  &
       & , mC2, U, V, mUSquark2, RUSquark, mDSquark2, RDSquark, mSlepton2    &
       & , RSlepton, mSneutrino2, mS02, RS0, mP02, RP0, mSpm2, RSpm          &
       & , mZ2_mZ, mW2_run, delta1                                           &
       & , mSneutrino_1L, mSneutrino2_1L, RSneutrino_1L, kont)
#endif

  If (kont.Ne.0) Then
   Iname = Iname -1
   mf_d = mf_d_save
   mf_l = mf_l_save
   mf_u = mf_u_save
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2
   Return
  End If
!------------------------------------------------------------
! shifting masses and mixing angles from tree-level to 1-loop
!------------------------------------------------------------
  mP0 = mP0_1L
  mP02 = mP02_1L

  mglu = mglu_1L

  mS0 = mS0_1L
  mS02 = mS02_1L
  RS0 = RS0_1L

  mSpm = mSpm_1L
  mSpm2 = mSpm2_1L

  mC = mC_1L
  mC2 = mC2_1L

  U = U_1L
  V = V_1L

  mN = mN_1L
  mN2 = mN2_1L
  N = N_1L

  mSneutrino = mSneutrino_1L
  mSneutrino2 = mSneutrino2_1L
  mUSquark = mUSquark_1L
  mUSquark2 = mUSquark2_1L
  mDSquark = mDSquark_1L
  mDSquark2 = mDSquark2_1L
  mSlepton = mSlepton_1L
  mSlepton2 = mSlepton2_1L
  RUsquark = RUsquark_1L
  RDSquark = RDSquark_1L
  RSlepton = RSlepton_1L
  RSneutrino = RSneutrino_1L
!------------------------------------------------
! replacing running fermion masses by pole masses
!------------------------------------------------
  mf_d = mf_d_save
  mf_l = mf_l_save
  mf_u = mf_u_save
  mf_d2 = mf_d**2
  mf_l2 = mf_l**2
  mf_u2 = mf_u**2

  Iname = Iname -1

 End Subroutine LoopMassesMSSM_3


 Subroutine Calcualte_MH2(delta, tanbQ, Q, g, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
   & , M2_E, M2_L, M2_D, M2_Q, M2_U, mu, mA2, IsOnShellMass, M2_H, B, kont)
!-----------------------------------------------------------------------
! In this subroutine the parameters M2_H and B  are calculated taking
! either the on-shell mass squared of the pseudoscalar Higgs boson as
! input or its running mass squared at the scale Q
! 07.03.12: written by Werner Porod
!-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3), A_l(3,3), A_d(3,3) &
    &  , A_u(3,3), M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_Q(3,3), M2_U(3,3), mu  &
    &  , Mi(3)
  Real(dp), Intent(in) :: delta, tanbQ, g(3), mA2, Q
  Logical, Intent(in) :: IsOnShellMass
  Real(dp), Intent(out) :: M2_H(2)
  Complex(dp), Intent(out) :: B 
 

  Real(dp) :: mC(2), mN(4), mS0(2), mSpm(2), mglu, mUsquark(6), mDsquark(6)  &
    & , mSlepton(6), mSneutrino(3), mUsquark2(6), mDsquark2(6), mSlepton2(6) &
    & , mSneutrino2(3), mS02(2), mSpm2(2), RP0(2,2), mC2(2), mN2(4), RS0(2,2)
  Complex(dp) :: U(2,2), V(2,2), N(4,4), RSpm(2,2), phase_Glu  &
    &    , RDsquark(6,6), RUsquark(6,6), RSlepton(6,6), RSneutrino(3,3)

  Real(dp) :: vevSM(2), cosb2, cos2b, sinb2, vev2, sinW2_DR, mZ2_mZ, vevs_DR(2) &
    & , tadpoles_1L(2), p2, b_A, mA2_mA, Pi2A0, tadpoles_2L(2), mW2_run         &
    & , mP0_T(2), mP02_T(2), mglu_sign, M2_H_1L(2), RP0_T(2,2)
  Real(dp) :: M2Ein, M2Lin, M2Din, M2Qin, M2Uin, delta1
  Complex(dp) :: dmZ2, PiA0
  Complex(dp), Dimension(3,3) :: uL_L, uL_R, uD_L, uD_R, uU_L, uU_R, CKM_Q
  Integer :: i1, i_loop
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'Calcualte_MH2'

  kont = 0
  mudim = GetRenormalizationScale() ! from LoopFunctions
!---------------------------------
! definning model
!---------------------------------
  n_S0 = 2
  n_P0 = 2
  n_Spm = 2
  n_char = 2
  n_neut = 4
  n_Sneut = 3
  n_Slept = 6

  uL_L = id3C
  uL_R = id3C
  uD_L = id3C
  uD_R = id3C
  uU_L = id3C
  uU_R = id3C
  !----------------------------------------------------------------------------
  ! require internally a factor 10 smaller delta to be sure, except of smaller
  ! than 100*Epsilon(1._dp) to avoid numerical problems
  !----------------------------------------------------------------------------
  delta1 = delta
  If (delta1.Gt.Epsilon(1._dp)*100._dp) delta1 = 0.1_dp * delta1

  cosb2 = 1._dp / (1._dp + tanbQ**2)
  sinb2 = 1._dp - cosb2
  cos2b = cosb2 - sinb2

  vev2 =  4._dp * mZ2 / (g(1)**2 + g(2)**2)
  vevSM(1) = Sqrt(vev2 * cosb2 )
  vevSM(2) = tanbQ * vevSM(1)
!------------------------------------------------
! replacing fermion pole masses by running masses
!------------------------------------------------
  mf_d_save = mf_d
  mf_l_save = mf_l
  mf_u_save = mf_u
  If (GenerationMixing) Then
   Call FermionMass(Y_l,vevSM(1),mf_l,uL_L, uL_R,kont)
   Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevSM, mf_d, uD_L, uD_R &
                                     & , mf_u, uU_L, uU_R, CKM_Q)
  Else
   Do i1=1,3
    mf_d(i1) = oosqrt2 * Abs(Y_d(i1,i1)) * vevSM(1)
    mf_l(i1) = oosqrt2 * Abs(Y_l(i1,i1)) * vevSM(1)
    mf_u(i1) = oosqrt2 * Abs(Y_u(i1,i1)) * vevSM(2)
   End Do
   CKM_Q = id3C
  End If
  mf_d2 = mf_d**2
  mf_l2 = mf_l**2
  mf_u2 = mf_u**2
!-------------------------------------------------
! saving old renormalisation scale and putting the
! new one
!-------------------------------------------------
  mudim = SetRenormaliZationScale(Q**2)
!-------------------------------------------------
! tree level masses and mixings as starting point
! first guess of B and M2_H
!-------------------------------------------------
  B =  mA2 * tanbQ / (1._dp + tanbQ**2)

  Call TreeMassesMSSM2(g(1), g(2), vevSM, Mi(1), Mi(2), Mi(3), mu, B, tanbQ   &
        , M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u          &
        , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                                  &
        , mGlu, Phase_Glu, mC, mC2, U, V, mN, mN2, N                          &
        , mSneutrino, mSneutrino2, Rsneutrino, mSlepton, mSlepton2, RSlepton  &
        , mDSquark, mDSquark2, RDSquark, mUSquark, mUSquark2, RUSquark        &
        , mP0_T, mP02_T, RP0_T, mS0, mS02, RS0, mSpm, mSpm2, RSpm, mZ2, mW2   &
        , GenerationMixing, kont, .True., .False.)

  kont = 0
  sinW2_DR = g(1)**2 / (g(1)**2+g(2)**2)

!----------------------
! mZ(mZ)
!----------------------
  Call PiZZT1(mZ2, g(2), sinW2_DR, vevSM, mZ2, mW2, mS02, RS0, mP02_T, RP0_T &
   , mSpm2 , RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton, mUSquark2   &
   , RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2, mC, mC2, U, V       &
   , mN, mN2, N ,dmZ2)
  vev2 =  4._dp * Real(mZ2+dmZ2,dp) / (g(1)**2 + g(2)**2)
  vevs_DR(1) = Sqrt(vev2 / (1._dp+tanbQ**2) )
  vevs_DR(2) = tanbQ * vevs_DR(1)
  mZ2_mZ = Real(dmZ2+mZ2,dp)
  mW2_run = mZ2_mZ * (1._dp - sinW2_DR)

  If (mZ2_mZ.Lt.0._dp) Then
   WriteWert = mZ2_mZ
   Call WriteLoopMassesError(4, "Calcualte_MH2", kont)
   mZ2_mZ = mZ2
  End If
!------------------------------------------------
! replacing fermion pole masses by running masses
!------------------------------------------------
   If (GenerationMixing) Then
    Call FermionMass(Y_l,vevs_DR(1),mf_l,uL_L, uL_R,kont)
    Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevs_DR, mf_d, uD_L, uD_R &
                                  & , mf_u, uU_L, uU_R, CKM_Q)
   Else
    Do i1=1,3
     mf_d(i1) = oosqrt2 * Abs(Y_d(i1,i1)) * vevs_DR(1)
     mf_l(i1) = oosqrt2 * Abs(Y_l(i1,i1)) * vevs_DR(1)
     mf_u(i1) = oosqrt2 * Abs(Y_u(i1,i1)) * vevs_DR(2)
    End Do
    CKM_Q = id3C
   End If
  mf_d2 = mf_d**2
  mf_l2 = mf_l**2
  mf_u2 = mf_u**2

  mglu_sign = mGlu * phase_Glu

  Do i_loop=1,2
  !-------------------------------------------------------------
  ! tree level masses and mixings using the loop corrected vevs
  !-------------------------------------------------------------

   Call PiZZT1(mZ2, g(2), sinW2_DR, vevs_DR, mZ2_mZ, mW2_run, mS02, RS0, mP02_T&
           , RP0_T, mSpm2, RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton  &
           , mUSquark2 , RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2    &
           , mC, mC2, U, V, mN, mN2, N ,dmZ2)
   vev2 =  4._dp * Real(mZ2+dmZ2,dp) / (g(1)**2 + g(2)**2)
   vevs_DR(1) = Sqrt(vev2 / (1._dp+tanbQ**2) )
   vevs_DR(2) = tanbQ * vevs_DR(1)
   mZ2_mZ = Real(dmZ2+mZ2,dp)
   mW2_run = mZ2_mZ * (1._dp - sinW2_DR)

   vev_Q = Sqrt(vev2) ! for output

   If (mZ2_mZ.Lt.0._dp) Then
    WriteWert = mZ2_mZ
    Call WriteLoopMassesError(-4, "Calcualte_MH2", kont)
    Iname = Iname - 1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If
   !------------------------------------------------
   ! replacing fermion pole masses by running masses
   !------------------------------------------------
   If (GenerationMixing) Then
    Call FermionMass(Y_l,vevs_DR(1),mf_l,uL_L, uL_R,kont)
    Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevs_DR, mf_d, uD_L, uD_R &
                                  & , mf_u, uU_L, uU_R, CKM_Q)
   Else
    Do i1=1,3
     mf_d(i1) = oosqrt2 * Abs(Y_d(i1,i1)) * vevs_DR(1)
     mf_l(i1) = oosqrt2 * Abs(Y_l(i1,i1)) * vevs_DR(1)
     mf_u(i1) = oosqrt2 * Abs(Y_u(i1,i1)) * vevs_DR(2)
    End Do
    CKM_Q = id3C
   End If
   mf_d2 = mf_d**2
   mf_l2 = mf_l**2
   mf_u2 = mf_u**2

   RP0 = RP0_T
   !---------------------------
   ! tadpoles at 1-loop
   !---------------------------
#ifdef GENERATIONMIXING
   Call One_Loop_Tadpoles_MSSM(g(1), g(2), vevs_DR, Y_l, Y_d, Y_u, uU_L, uU_R &
    & , uD_L, uD_R, uL_L, uL_R, mu, A_l, A_d, A_u, mSneutrino2, Rsneutrino    &
    & , mSlepton2, Rslepton, mDSquark2, RDSquark, mUSquark2, RUSquark         &
    & , mSpm2, RSpm, mC, mC2, U, V, mP02_T, RP0, mS02, RS0, mN, mN2, N        &
    & , mZ2_mZ, mW2_run, tadpoles_1L)
#else
   Call One_Loop_Tadpoles_MSSM(g(1), g(2), vevs_DR, Y_l, Y_d, Y_u, mu, A_l   &
    & , A_d, A_u, mSneutrino2, mSlepton2, Rslepton, mDSquark2, RDSquark      &
    & , mUSquark2, RUSquark, mSpm2, RSpm, mC, mC2, U, V, mP02_T, RP0, mS02   &
    & , RS0, mN, mN2, N, mZ2_mZ, mW2_run, tadpoles_1L)
#endif
   !---------------------------
   ! tadpoles at 2-loop
  !---------------------------------------------------------------------
  ! Define first input, numerical problems can occur if left and right
  ! sfermion mass parameters are identical
  !---------------------------------------------------------------------
   M2Din = Real(M2_D(3,3),dp)
   M2Ein = Real(M2_E(3,3),dp)
   M2Lin = Real(M2_L(3,3),dp)
   M2Qin = Real(M2_Q(3,3),dp)
   M2Uin = Real(M2_U(3,3),dp)
   If (M2Ein.Eq.M2Lin) M2Ein = (1._dp + Max(1.0e-8_dp,delta1)) * M2Ein
   If (M2Din.Eq.M2Qin) M2Din = (1._dp + Max(1.0e-8_dp,delta1)) * M2Din
   If (M2Uin.Eq.M2Qin) M2Uin = (1._dp + Max(1.0e-8_dp,delta1)) * M2Uin

   If (Only_1loop_Higgsmass) Then
    tadpoles_2L = 0._dp
   Else
    Call Two_Loop_Tadpoles_MSSM(g(3), mglu_sign, mP02_T(2), vevs_DR  &
        & , M2Din, M2Uin, M2Qin, M2Ein, M2Lin, A_d(3,3), A_u(3,3)    &
       & , A_l(3,3), Y_d(3,3), Y_u(3,3), Y_l(3,3), mu, tadpoles_2L, kont)
   End If

   If (kont.Ne.0) Then
    Iname = Iname - 1
    mf_d = mf_d_save
    mf_l = mf_l_save
    mf_u = mf_u_save
    mf_d2 = mf_d**2
    mf_l2 = mf_l**2
    mf_u2 = mf_u**2
    Return
   End If

   b_a = tadpoles_1L(1) * sinb2 + tadpoles_1L(2) * cosb2

   If (IsOnShellMass) Then  
   !--------------------
   ! pseudoscalar Higgs 
   !--------------------
    p2 =  mP02(2)

#ifdef GENERATIONMIXING
    Call PiPseudoScalar(p2, g(1), g(2), vevs_DR, mP02_T, RP0               &
       & , Y_d, uD_L, uD_R, Y_u, uU_L, uU_R, Y_l, uL_L, uL_R               &
       & , A_d, A_u, A_l, mu                                               & 
       & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSlepton2, RSlepton   &
       & , mSneutrino2, RSneutrino, mSpm2, RSpm, mC, mC2, U, V, mS02, RS0  &
       & , mN, mN2, N, mZ2_mZ, mW2_run, PiA0)
#else
    Call PiPseudoScalar(p2, g(1), g(2), vevs_DR, mP02_T, RP0, Y_d, Y_u, Y_l   &
       & , A_d, A_u, A_l, mu , mUSquark2, RUSquark, mDSquark2, RDSquark       & 
       & , mSlepton2, RSlepton, mSneutrino2, mSpm2, RSpm, mC, mC2, U, V, mS02 &
       & , RS0, mN, mN2, N, mZ2_mZ, mW2_run, PiA0)
#endif

    If (Only_1loop_Higgsmass) Then
     Pi2A0 = 0._dp
    Else
     Call PiPseudoScalar2(g(3), mglu_sign, mP02_T(2), vevs_DR, M2Din, M2Uin &
      & , M2Qin, M2Ein, M2Lin, A_d(3,3), A_u(3,3), A_l(3,3), Y_d(3,3)       &
      & , Y_u(3,3), Y_l(3,3), mu, Pi2A0, kont)
    End If

    mA2_mA = -Pi2A0+ PiA0 + mA2 - b_a 

    mA2_Q = -Pi2A0+ PiA0 + mA2  ! for output
   Else  ! I have already the running mass
    mA2_mA = mA2 - b_a
   End If

   M2_H_1L(1) = 0.5_dp * ( mA2_mA - 2._dp * Abs(mu)**2    &
              &          - (mZ2_mZ + mA2_mA) * cos2b )
   M2_H_1L(2) = 0.5_dp * ( mA2_mA - 2._dp * Abs(mu)**2    &
              &          + (mZ2_mZ + mA2_mA) * cos2b )

   M2_H = M2_H_1L + tadpoles_1L - tadpoles_2L
   B = (M2_H_1L(1) + M2_H_1L(2) + 2._dp *  Abs(Mu)**2) * tanbQ / (1+tanbQ**2)
  End Do ! i_loop

!------------------------------------------------
! replacing running fermion masses by pole masses
!------------------------------------------------
  mf_d = mf_d_save
  mf_l = mf_l_save
  mf_u = mf_u_save
  mf_d2 = mf_d**2
  mf_l2 = mf_l**2
  mf_u2 = mf_u**2

!------------------------------------------------
! resetting renormalisation scale
!------------------------------------------------
  mudim = SetRenormaliZationScale(mudim)

  Iname = Iname -1

 End Subroutine Calcualte_MH2



#ifdef GENERATIONMIXING
 Subroutine CharginoMass_Loop(mC, mC2, U, V, gU1, gSU2, Y_D, Y_L, Y_U, M2, mu &
      & , vevs_in , mN, mN2, N, mSpm2, RSpm, mS02, RS0, mP02, RP0, mSneut2    &
      & , RSneut, uD_L, uD_R, uU_L, uU_R, uL_L, uL_R                          &
      & , mSlepton2, RSlepton, mSdown2, RSdown, mSup2, RSup, mZ2, mW2, delta  &
      & , mC1L, mC21L, U1L, V1L, kont)
#else
 Subroutine CharginoMass_Loop(mC, mC2, U, V, gU1, gSU2, Y_D, Y_L, Y_U, M2, mu &
      & , vevs_in , mN, mN2, N, mSpm2, RSpm, mS02, RS0, mP02, RP0, mSneut2    &
      & , mSlepton2, RSlepton, mSdown2, RSdown, mSup2, RSup, mZ2, mW2, delta  &
      & , mC1L, mC21L, U1L, V1L, kont)
#endif
 !-----------------------------------------------------------------
 ! calculates chargino masses + mixing matrices U,V in the MSSM,
 ! written by Werner Porod, 1.8.99
 ! 03.11.01: portation to f90
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: gU1, gSU2, vevs_in(2), mC(2), mC2(2), mN(4), mN2(4) &
      & , mSneut2(3), mSlepton2(6), mSdown2(6), mSup2(6), mS02(2), RS0(2,2)   &
      & , mP02(2), RP0(2,2), mSpm2(2), mZ2, mW2, delta
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: Y_D(3,3), Y_L(3,3), Y_U(3,3), U(2,2), N(4,4)  &
      & , RSpm(2,2), RSlepton(6,6), RSdown(6,6), RSup(6,6)                 &
      & , uD_L(3,3), uD_R(3,3), uU_L(3,3), uU_R(3,3), uL_L(3,3), uL_R(3,3) &
      & , RSneut(3,3), M2, mu, V(2,2)
#else
  Complex(dp), Intent(in) :: Y_D(3,3), Y_L(3,3), Y_U(3,3), U(2,2), N(4,4)  &
      & , RSpm(2,2), RSlepton(6,6), RSdown(6,6), RSup(6,6), M2, mu, V(2,2)
#endif
  Real(dp), Intent(out) :: mC1L(2), mC21L(2)
  Complex(dp), Intent(out) :: U1L(2,2), V1L(2,2)

  Integer :: i1, i2, i3, ierr, i_l, i_count
  Real(dp) :: p2, sinW2_Q, v1(2,2), u1(2,2), test(2), cosW, mC2_L(2)
  Complex(dp) :: mat2(2,2), phaseM, mat2a(2,2), SigL(2,2), SigR(2,2) &
      & , SigS(2,2), mat22(2,2), U2(2,2), V2(2,2), coupLC, coupRC
  Complex(dp) :: c_SmpCN_L(2,2,4), c_SmpCN_R(2,2,4), c_CNuSl_L(2,3,6)         &
      & , c_CNuSl_R(2,3,6), c_CLSn_L(2,3,3), c_CLSn_R(2,3,3), c_CUSd_L(2,3,6) &
      & , c_CUSd_R(2,3,6), c_CDSu_L(2,3,6), c_CDSu_R(2,3,6), Rsd(2,2)         &
      & , Rsu(2,2), Rsl(2,2), yuk, YukP
#ifdef GENERATIONMIXING
  Complex(dp) :: mat3(3,3), mat6(6,6)
#endif
  Complex(dp) :: c_CCP0_L(2,2,2), c_CCP0_R(2,2,2), c_CCS0_L(2,2,2)   &
      & , c_CCS0_R(2,2,2), c_CCZ_L(2,2), c_CCZ_R(2,2), c_CNW_L(2,4)  &
      & , c_CNW_R(2,4), c_CCG_L(2,2), c_CCG_R(2,2)
  Logical :: WriteOut

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CharginoMass_Loop'

 !----------------
 ! Initialization
 !----------------
  If ((WriteOneLoopContributions.Eq.7).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions in CharginoMass_Loop:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If

  mC1L = 0._dp
  mC21L = 0._dp
  U1L = (0._dp,0._dp)
  V1L = (0._dp,0._dp)
  !---------------------
  ! tree level
  !---------------------
  mat2a(1,1) = M2
  mat2a(1,2) = gSU2 * vevs_in(2) * oosqrt2
  mat2a(2,1) = gSU2 * vevs_in(1) * oosqrt2
  mat2a(2,2) = mu
  If (WriteOut) Then
   Write(ErrCan,*) "Tree level mass matrix"
   Write(ErrCan,*) mat2a(1,:)
   Write(ErrCan,*) mat2a(2,:)
  End If
  !----------------------------------
  ! couplings for 1-loop calculation
  !----------------------------------
  c_CNuSl_L = 0._dp
  c_CNuSl_R = 0._dp
  c_CLSn_L = 0._dp
  c_CLSn_R = 0._dp
  c_CUSd_L = 0._dp
  c_CUSd_R = 0._dp
  c_CDSu_L = 0._dp
  c_CDSu_R = 0._dp


#ifdef GENERATIONMIXING
  If (GenerationMixing) Then

   mat6 = 0._dp
   mat6(1:3,1:3) = Rsneut
   mat3 = 0._dp
   Do i1=1,2
    Do i2=1,3
     Do i3=1,6
      Call CoupCharginoSfermion(i1, i2, i3, gSU2, 0.5_dp, RSlepton, Y_l &
                              &, mat3, id3c, id3c, id2c, id2c, coupLC, coupRC)
      c_CNuSl_L(i1, i2, i3) = coupLC
      c_CNuSl_R(i1, i2, i3) = coupRC
      Call CoupCharginoSfermion(i1, i2, i3, gSU2, 0.5_dp, RSdown, Y_d &
                              &, Y_u, uU_L, uU_R, id2c, id2c, coupLC, coupRC)
      c_CUSd_L(i1, i2, i3) = coupLC
      c_CUSd_R(i1, i2, i3) = coupRC
      Call CoupCharginoSfermion(i1, i2, i3, gSU2, -0.5_dp, RSup, Y_d &
                              &, Y_u, uD_L, uD_R, id2c, id2c, coupLC, coupRC)
      c_CDSu_L(i1, i2, i3) = coupLC
      c_CDSu_R(i1, i2, i3) = coupRC
     End Do
     Do i3=1,3
      Call CoupCharginoSfermion(i1, i2, i3, gSU2, -0.5_dp, mat6, Y_l, mat3 &
                              &, uL_L, uL_R, id2c, id2c, coupLC, coupRC)
      c_CLSn_L(i1, i2, i3) = coupLC
      c_CLSn_R(i1, i2, i3) = coupRC
     End Do
    End Do
   End Do

  Else
#endif

   Do i1=1,3
    Rsd = RSdown(2*i1-1:2*i1, 2*i1-1:2*i1)
    Rsu = RSup(2*i1-1:2*i1, 2*i1-1:2*i1)
    Rsl = RSlepton(2*i1-1:2*i1, 2*i1-1:2*i1)

    Yuk = Y_l(i1,i1)
    YukP = 0._dp
    Do i2=1,2
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSl, Yuk, Yukp, id2c, &
                              & id2c, coupLC, coupRC)
      c_CNuSl_L(i2, i1, (i1-1)*2 + i3) = coupLC
      c_CNuSl_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
     Call CoupCharginoSfermion(i2, 1, gSU2, -0.5_dp, id2C, Yuk, Yukp, id2c, &
                              & id2c, coupLC, coupRC)
     c_CLSn_L(i2, i1, i1) = coupLC
     c_CLSn_R(i2, i1, i1) = coupRC
    End Do

    Yuk = Y_d(i1,i1)
    YukP = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSd, Yuk, Yukp, id2c, &
                              & id2c, coupLC, coupRC)
      c_CUSd_L(i2, i1, (i1-1)*2 + i3) = coupLC
      c_CUSd_R(i2, i1, (i1-1)*2 + i3) = coupRC
      Call CoupCharginoSfermion(i2, i3, gSU2, -0.5_dp, RSu, Yuk, Yukp, id2c, &
                              & id2c, coupLC, coupRC)
      c_CDSu_L(i2, i1, (i1-1)*2 + i3) = coupLC
      c_CDSu_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
    End Do

   End Do

#ifdef GENERATIONMIXING
  End If
#endif

  c_SmpCN_L = 0._dp
  c_SmpCN_R = 0._dp
  Do i1 = 1,2
   Do i2 = 1,2
    Do i3 = 1,4
     Call CoupCSCharginoNeutralino(i1, i2, i3, N, id2c, id2c, RSpm     &
               &, gU1, gSU2, c_SmpCN_L(i1,i2,i3), c_SmpCN_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  c_CCP0_L = 0.0_dp
  c_CCP0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_P0
    Call CoupCharginoPseudoScalara(i1, i2, i3, U, V, RP0, gSU2  &
                     &, c_CCP0_L(i1,i2,i3), c_CCP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  c_CCS0_L = 0.0_dp
  c_CCS0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_S0
    Call CoupCharginoScalara(i1, i2, i3, U, V, RS0, gSU2  &
                     &, c_CCS0_L(i1,i2,i3), c_CCS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  c_CCZ_L = 0.0_dp
  c_CCZ_R = 0.0_dp
  cosw = gSU2 / Sqrt(gSU2**2 + gU1**2)
  Do i1=1,n_char
   Do i2=1,n_char
    Call CoupCharginoZa(i1, i2, U, V, gSU2, cosW  &
                     &, c_CCZ_L(i1,i2), c_CCZ_R(i1,i2) )
   End Do
  End Do

  c_CNW_L = 0._dp
  c_CNW_R = 0._dp

  Do i1 = 1,n_char
   Do i2= 1,n_neut
    Call CoupCharginoNeutralinoW(i1, i2, N, id2C, id2C, gSU2 &
                               &, c_CNW_L(i1,i2), c_CNW_R(i1,i2) )
   End Do
  End Do
  c_CCG_L = 0._dp
  c_CCG_R = 0._dp
  sinW2_Q = gSU2 * gU1 / Sqrt(gSU2**2 + gU1**2)
  Do i1 = 1,n_char
   Do i2= 1,n_char
    c_CCG_L(i1,i2) = sinW2_Q * Conjg( U(i2,i1) )
    c_CCG_R(i1,i2) = sinW2_Q * V(i2,i1)
   End Do
  End Do
  !--------------------
  ! 1-loop calculation
  !--------------------
  Do i_l=n_char,1,-1
   p2 = mC2(i_l)
   Call Sigma_Chargino(p2, mC, mC2, mSdown2, c_CUSd_L, c_CUSd_R             &
          , mSup2, c_CDSu_L, c_CDSu_R, mSlepton2, c_CNuSl_L, c_CNuSl_R      &
          , mSneut2, c_CLSn_L, c_CLSn_R, mN, mN2, c_CNW_L, c_CNW_R          &
          , mSpm2, c_SmpCN_L, c_SmpCN_R, c_CCZ_L, c_CCZ_R, c_CCG_L, c_CCG_R &
          , mP02, c_CCP0_L, c_CCP0_R, mS02, c_CCS0_L, c_CCS0_R, mZ2, mW2    &
          , WriteOut, SigL, SigR, SigS)

   mat2 = mat2a - SigS - Matmul(SigR,mat2a) - Matmul(mat2a,SigL)

   ierr = 0
   mat22 = Matmul( Transpose( Conjg( mat2 ) ), mat2 )
   If ( Maxval( Abs( Aimag(mat22) ) ).Eq.0._dp) Then ! reel matrix
    Call EigenSystem(Real(mat22,dp), mC2_L, v1, ierr, test)
    v2 = v1
   Else
    Call EigenSystem(mat22, mC2_L, v2, ierr, test)
   End If

   If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(ErrCan,*) "test =",test
     Write(ErrCan,*) " "
     If (ErrorLevel.Eq.2) Call TerminateProgram
     ierr = 0
   Else If (ierr.Ne.0) Then
    Write(ErrCan,*) 'Warning in subroutine CharginoMass_Loop, ierr = ',ierr
    Write(ErrCan,*) 'M2, mu ',M2,mu
    Write(ErrCan,*) 'gU1, gSU2',gU1, gSU2
    Write(ErrCan,*) 'vevs ',vevs_in
    Write(ErrCan,*) ' '
    kont = ierr
    mC1L(1) = Abs(M2)
    mC1L(2) = Abs(mu)
    mC21L = mC1L**2
    U1L = id2C
    V1L = U1L
    Iname = Iname - 1
    Return
   End If
 
   !---------------------------------------
   ! redoing calculation using 1-loop mass 
   !---------------------------------------
   i_count = 0
   p2_loop: Do
    i_count = i_count + 1 
    p2 = mC2_L(i_l)
    Call Sigma_Chargino(p2, mC, mC2, mSdown2, c_CUSd_L, c_CUSd_R            &
          , mSup2, c_CDSu_L, c_CDSu_R, mSlepton2, c_CNuSl_L, c_CNuSl_R      &
          , mSneut2, c_CLSn_L, c_CLSn_R, mN, mN2, c_CNW_L, c_CNW_R          &
          , mSpm2, c_SmpCN_L, c_SmpCN_R, c_CCZ_L, c_CCZ_R, c_CCG_L, c_CCG_R &
          , mP02, c_CCP0_L, c_CCP0_R, mS02, c_CCS0_L, c_CCS0_R, mZ2, mW2    &
          , WriteOut, SigL, SigR, SigS)

    mat2 = mat2a - SigS - Matmul(SigR,mat2a) - Matmul(mat2a,SigL)

    ierr = 0
    mat22 = Matmul( Transpose( Conjg( mat2 ) ), mat2 )
    If ( Maxval( Abs( Aimag(mat22) ) ).Eq.0._dp) Then ! reel matrix
     Call EigenSystem(Real(mat22,dp), mC2_L, v1, ierr, test)
     v2 = v1
    Else
     Call EigenSystem(mat22, mC2_L, v2, ierr, test)
    End If

    mC21L(i_l) = mC2_L(i_l)
    mC1L(i_l) = Sqrt( mC21L(i_l) )
    If (p2.Ne.0._dp) Then
     test(1) = Abs(mC2_L(i_l) - p2) / p2
    Else
     test(1) = Abs(mC2_L(i_l))
    End If
    If (test(1).Lt.0.1_dp*delta) Exit p2_loop
    If (i_count.Gt.30) Then
     Write(*,*) "Problem in loop for mC",i_l,test(1),mC2_L(i_l)
     Exit p2_loop
    End If
   End Do p2_loop

   If (i_l.Gt.1)  Cycle ! Currently I take only the mixing matrices of chi_1


   mat22 = Matmul( mat2, Transpose( Conjg( mat2 ) ) )
   If ( Maxval( Abs( Aimag(mat22) ) ).Eq.0._dp) Then ! reel matrix
    Call EigenSystem(Real(mat22,dp), mC2_L, u1, ierr, test)
    u2 = u1
   Else
    Call EigenSystem(mat22, mC2_L, u2, ierr, test)
   End If
   u2 = Conjg(u2)

   If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(ErrCan,*) "test =",test
     Write(ErrCan,*) " "
     If (ErrorLevel.Eq.2) Call TerminateProgram
     ierr = 0
   Else If (ierr.Ne.0) Then
    Write(ErrCan,*) 'Warning in subroutine CharginoMass_Loop, ierr = ',ierr
    Write(ErrCan,*) 'M2, mu ',M2,mu
    Write(ErrCan,*) 'gU1, gSU2',gU1, gSU2
    Write(ErrCan,*) 'vevs ',vevs_in
    Write(ErrCan,*) ' '
    kont = ierr
    mC1L(1) = Abs(M2)
    mC1L(2) = Abs(mu)
    mC21L = mC1L**2
    U1L = id2C
    V1L = U1L
    Iname = Iname - 1
    Return
   Endif

   mat22 = Matmul( Matmul( Conjg(u2), mat2), Transpose( Conjg(v2) ) )
   Do i1=1,2
    phaseM =   mat22(i1,i1)   / Abs( mat22(i1,i1) )
    v2(i1,:) = phaseM * v2(i1,:)
   End Do


   If (Abs(v2(1,2)).Gt.0._dp) Then
    phaseM = v2(1,2) / Abs(v2(1,2))
    v2(1,2) = Abs(v2(1,2))
    v2(1,1) = v2(1,1) * Conjg(phaseM)
    u2(1,1) = u2(1,1) * phaseM
    u2(1,2) = u2(1,2) * phaseM
   End If

   If (Abs(v2(2,1)).Gt.0._dp) Then
    phaseM = v2(2,1) / Abs(v2(2,1))
    v2(2,1) = Abs(v2(2,1))
    v2(2,2) = v2(2,2) * Conjg(phaseM)
    u2(2,2) = u2(2,2) * phaseM
    u2(2,1) = u2(2,1) * phaseM
   End If

   If (WriteOut) Then
    mat22 = Matmul( Matmul( Conjg(u2), mat2), Transpose( Conjg(v2) ) )
    Write(ErrCan,*) "test",Sqrt(mc2)
    Write(ErrCan,*) "    ",mat22(1,:)
    Write(ErrCan,*) "    ",mat22(2,:)
    Write(ErrCan,*) " "
   End If

   U1L = u2
   V1L = v2

  End Do ! i_l

  If (WriteOut) Then
   Write(ErrCan,*) "mC_1L",mC1L
   Write(ErrCan,*) "U_1L_1i",U1L(1,:)
   Write(ErrCan,*) "U_1L_2i",U1L(2,:)
   Write(ErrCan,*) "V_1L_1i",V1L(1,:)
   Write(ErrCan,*) "V_1L_2i",V1L(2,:)
   Write(ErrCan,*) " "
  End If

  Iname = Iname - 1

 End Subroutine CharginoMass_Loop


 Subroutine delta_VB(gSU2, sinW2, sinW2_DR, rho, mC, mC2, U, V, mN, mN2, N &
     & , Y_l, UL_L, UL_R, mSlept2, Rslept, mSneut2, RSneut, UNu_Lin, res)
 !-----------------------------------------------------------------------
 ! Calculates the the nonuniversal corrections to delta_r
 ! The formula of J. Bagger et al, Nucl.Phys.B is used. The renormalization
 ! scale has to be set in the main program due to the structure of LoopTools. 
 ! the input is:
 ! - gSU2 ....... the SU(2) gauge coupling at p2
 ! - sinW2 ...... sin(theta_W) squared
 ! - sinW2_DR ... sin(theta_W) squared, DR-scheme
 ! - rho ........ Rho parameter
 ! - mC ......... masses of the charginos
 ! - mC2 ........ masses of the charginos squared
 ! - U,V ........ mixing matrices of the charginos
 ! - mN ......... masses of the neutralinos
 ! - mN2 ........ masses of the neutralinos squared
 ! - N .......... mixing matrix of the neutrlinos
 ! - Y_l ........ lepton yukawas
 ! - UL_L, UL_R . lepton mixing matrices
 ! - mSlept2 .... masses of the sleptons squared
 ! - RSlept ..... mixing matrix of the sleptons
 ! - mSneut2 .... masses of the sneutrinos squared
 ! - RSneut ..... mixing matrix of the sneutrinos
 ! - UNu_L ...... neutrino mixing matrix
 ! output
 !  res    
 ! written by Werner Porod, 30.11.01
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: sinW2, sinW2_DR, gSU2, rho, mC(2), mC2(2), mN(4) &
     & , mN2(4), mSneut2(3), mSlept2(6)
  Complex(dp), Intent(in) :: RSlept(6,6), U(2,2), V(2,2), N(4,4)
  Complex(dp), Intent(in), Dimension(3,3) :: RSneut, Y_l, UL_L, UL_R, UNu_Lin
  Real(dp), Intent(out) :: res

  Integer :: i1, i2, i3, i4, i5, i6
  Real(dp) :: cosW2, cosW2_DR, g2, gp
  Complex(dp) :: sumI, c_CNuSl_L(2,3,6), c_CNuSl_R(2,3,6), c_LNSl_L(3,4,6) &
     & , c_LNSl_R(3,4,6), B1_C_Sl(2,6), B1_N_Sl(4,6), B1_C_Sn(2,3)         &
     & , B1_N_Sn(4,3), mat3(3,3), c_CLSn_L(2,3,3), c_CLSn_R(2,3,3)         &
     & , c_NuNSn_R(3,4,3), c_CNW_L(2,4), c_CNW_R(2,4), B0_CN(2,4)          &
     & , C0_SnCN(3,2,4), C0_SlCN(6,2,4), teil, D0_SlSnCN(6,3,2,4)          &
     & , D27_SlSlCN(6,6,2,4), D27_SnSnCN(3,3,2,4), D0_SlSnCC(6,3,2,2)      &
     & , D0_SlSnNN(6,3,4,4), D27_SlSnNN(6,3,4,4), Rsel(2,2), Rsmu(2,2)     &
     & , C0_NSlSn(4,6,3), B0_SlSn(6,3), RSlSn, UNu_L(3,3), phase, sumIJ(3)
  Logical :: WriteOut

  Iname = Iname + 1
  NameOfUnit(Iname) = "delta_VB"

  WriteOut = .False.
!  WriteOut = .True.
  If ( (WriteOneLoopContributions.Eq.9).Or.(WriteOneLoopContributions.Lt.0)) &
   &  WriteOut = .True.
  
  !----------------------------------------------------------------------------
  ! phase conventions, only needed for neutrinos and charged leptons
  ! is not needed for sfermions as an overall phase change cancels in the loops
  !----------------------------------------------------------------------------
  Do i1=1,3
   phase = (1._dp,0._dp) ! in case none of the conditions below apply, do nothing
   
   If (Abs(UNu_Lin(i1,i1)).Ne.0._dp) Then ! first diagonal, then potentially the others
    phase = Conjg(UNu_Lin(i1,i1))/Abs(UNu_Lin(i1,i1))
   Else If (Abs(UNu_Lin(i1,1)).Ne.0._dp) Then
    phase = Conjg(UNu_Lin(i1,1))/Abs(UNu_Lin(i1,1))
   Else If (Abs(UNu_Lin(i1,2)).Ne.0._dp) Then
    phase = Conjg(UNu_Lin(i1,2))/Abs(UNu_Lin(i1,2))
   Else If (Abs(UNu_Lin(i1,3)).Ne.0._dp) Then
    phase = Conjg(UNu_Lin(i1,3))/Abs(UNu_Lin(i1,3))
   End If

   UNu_L(i1,:) = UNu_Lin(i1,:) * phase
  End Do

  If (WriteOut) Write(ErrCan,*) &
   & "Contributions to the non-universal corrections of the rho parameter."
  !--------------------------
  ! SM contribution
  !--------------------------
  cosW2 = 1._dp - sinW2
  cosW2_DR = 1._dp - sinW2_DR
  g2 = gSU2**2

  sumI = 6._dp                                                       &
     & + Log(cosW2) * (3.5_dp - 2.5_dp * sinW2                       &
     &                - sinW2_DR * (5._dp - 1.5_dp*cosW2/cosW2_DR)   &
     &                ) / sinW2

  sumI = sumI * g2 * rho
  If (WriteOut) Write(ErrCan,*) "SM :",sumI
  res = sumI
  !------------------------------------------------------------------------
  !  some preperations, attention, Pierce et al. define B1 with
  !  a relative minus sign 
  !-----------------------------------------------------------------------
  mat3 = 0._dp ! zero matrix, needed in case of sneutrinos
  gp = gSU2 * Sqrt(sinW2_DR/cosW2_DR)
  sumI = 0._dp

  c_CNuSl_L = 0._dp
  c_CNuSl_R = 0._dp
  c_LNSl_L = 0._dp
  c_LNSl_R = 0._dp
  c_CNW_L = 0._dp
  c_CNW_R = 0._dp
  B1_C_Sl = 0._dp
  B1_N_Sl = 0._dp
  B1_C_Sn = 0._dp
  B1_N_Sn = 0._dp
  Do i1=1,2
   Do i2=1,4
    Call CoupCharginoNeutralinoW(i1, i2, N, U, V, gSU2 &
                               &, c_CNW_L(i1,i2), c_CNW_R(i1,i2))
    B0_CN(i1,i2) = B0(0._dp, mC2(i1), mN2(i2) )
   End Do
  End Do

#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i3=1,3
    Do i1=1,6
     Do i2=1,2
      Call CoupCharginoSfermion(i2, i3, i1, gSU2, 0.5_dp, RSlept, Y_l, mat3  &
           &, UNu_L, mat3, U, V, c_CNuSl_L(i2, i3, i1), c_CNuSl_R(i2, i3, i1) )
     End Do
     Do i2=1,4
      Call CoupNeutralinoSlepton(i3, i2, i1, gp, gSU2, RSlept, UL_L, UL_R  &
                  &, Y_l, N, c_LNSl_L(i3, i2, i1), c_LNSl_R(i3, i2, i1) )
     End Do    
    End Do 
    Do i1=1,3
     Do i2=1,2
      Call CoupCharginoSfermion(i2, i3, i1, gSU2, -0.5_dp, RSneut, Y_l, mat3 &
           & , UL_L, UL_R, U, V, c_CLSn_L(i2,i3,i1), c_CLSn_R(i2,i3,i1) )
     End Do
     Do i2=1,4
      Call CoupNeutralinoSneutrino(i3, i2, i1, gp, gSU2, N, RSneut, UNu_L &
                                 &, c_NuNSn_R(i3,i2,i1) )
     End Do
    End Do
   End Do 
   Do i1=1,6
    Do i2=1,2
     B1_C_Sl(i2,i1) = B1(0._dp, mC2(i2), mSlept2(i1) )
    End Do
    Do i2=1,4
     B1_N_Sl(i2,i1) = B1(0._dp, mN2(i2), mSlept2(i1) )
     Do i3=1,2
      C0_SlCN(i1,i3,i2) = C0_3m(mSlept2(i1), mC2(i3), mN2(i2))
     End Do
     Do i3=1,3
      C0_NSlSn(i2,i1,i3) = C0_3m(mN2(i2), mSlept2(i1), mSneut2(i3) )
     End Do
    End Do
   End Do
   Do i1=1,3
    Do i2=1,2
     B1_C_Sn(i2,i1) = B1(0._dp, mC2(i2), mSneut2(i1) )
    End Do
    Do i2=1,4
     B1_N_Sn(i2,i1) = B1(0._dp, mN2(i2), mSneut2(i1) )
     Do i3=1,2
      C0_SnCN(i1,i3,i2) = C0_3m(mSneut2(i1), mC2(i3), mN2(i2))
     End Do
    End Do
   End Do

   Do i1=1,2
    Do i2=1,4
     Do i3=1,6
      Do i4=1,3
       D0_SlSnCN(i3,i4,i1,i2) =  &
                 & D0_Bagger(mSlept2(i3), mSneut2(i4), mC2(i1), mN2(i2) )
      End Do
      Do i4=1,6
       D27_SlSlCN(i3,i4,i1,i2) =  &
                 & D27_Bagger(mSlept2(i3), mSlept2(i4), mC2(i1), mN2(i2) )
      End Do
     End Do
     Do i3=1,3
      Do i4=1,3
        D27_SnSnCN(i3,i4,i1,i2) =  &
                 & D27_Bagger(mSneut2(i3), mSneut2(i4), mC2(i1), mN2(i2) )
      End Do
     End Do
    End Do
   End Do
   Do i1=1,6
    Do i2=1,3
     B0_SlSn(i1,i2) = B0(0._dp, mSlept2(i1), mSneut2(i2) )
     Do i3=1,2
      Do i4=1,2
       D0_SlSnCC(i1,i2,i3,i4) =   &
                 & D0_Bagger(mSlept2(i1), mSneut2(i2), mC2(i3), mC2(i4) )
      End Do
     End Do
     Do i3=1,4
      Do i4=1,4
       D0_SlSnNN(i1,i2,i3,i4) =   &
                 & D0_Bagger(mSlept2(i1), mSneut2(i2), mN2(i3), mN2(i4) )
       D27_SlSnNN(i1,i2,i3,i4) =   &
                 & D27_Bagger(mSlept2(i1), mSneut2(i2), mN2(i3), mN2(i4) )
      End Do
     End Do
    End Do
   End Do

  Else ! .not.GenerationMixing
#endif

   RSel = RSlept(1:2,1:2)
   RSmu = RSlept(3:4,3:4)

   Do i1=1,2
    Do i2=1,2
     Call CoupCharginoSfermion(i2, i1, gSU2, 0.5_dp, RSel, Y_l(1,1), ZeroC &
                           &, U, V, c_CNuSl_L(i2, 1, i1), c_CNuSl_R(i2, 1, i1))
     Call CoupCharginoSfermion(i2, i1, gSU2, 0.5_dp, RSmu, Y_l(2,2), ZeroC &
                     &, U, V, c_CNuSl_L(i2, 2, i1+2), c_CNuSl_R(i2, 2, i1+2))
     Call CoupCharginoSfermion(i2, 1, gSU2, -0.5_dp, id2c, Y_l(i1,i1), ZeroC &
                             &, U, V, c_CLSn_L(i2,i1,i1), c_CLSn_R(i2,i1,i1) )
    End Do
    Do i2=1,4
     Call CoupNeutralinoSlepton(i2, i1, gp, gSU2, RSel, Y_l(1,1), N &
                               &, c_LNSl_L(1, i2, i1), c_LNSl_R(1, i2, i1) )
     Call CoupNeutralinoSlepton(i2, i1, gp, gSU2, RSel, Y_l(2,2), N &
                             &, c_LNSl_L(2, i2, i1+2), c_LNSl_R(2, i2, i1+2) )
     Call CoupNeutralinoSneutrino(i2, gp, gSU2, N, c_NuNSn_R(i1,i2,i1) )
    End Do
   End Do

   Do i1=1,4
    Do i2=1,2
     B1_C_Sl(i2,i1) = B1(0._dp, mC2(i2), mSlept2(i1) )
    End Do
    Do i2=1,4
     B1_N_Sl(i2,i1) = B1(0._dp, mN2(i2), mSlept2(i1) )
     Do i3=1,2
      C0_SlCN(i1,i3,i2) = C0_3m(mSlept2(i1), mC2(i3), mN2(i2))
     End Do
     Do i3=1,2
      C0_NSlSn(i2,i1,i3) = C0_3m(mN2(i2), mSlept2(i1), mSneut2(i3) )
     End Do
    End Do
   End Do
   Do i1=1,2
    Do i2=1,2
     B1_C_Sn(i2,i1) = B1(0._dp, mC2(i2), mSneut2(i1) )
    End Do
    Do i2=1,4
     B1_N_Sn(i2,i1) = B1(0._dp, mN2(i2), mSneut2(i1) )
     Do i3=1,2
      C0_SnCN(i1,i3,i2) = C0_3m(mSneut2(i1), mC2(i3), mN2(i2))
     End Do
    End Do
   End Do
   D27_SnSnCN = 0._dp
   Do i1=1,2
    Do i2=1,4
     Do i3=1,4
      Do i4=1,2
       D0_SlSnCN(i3,i4,i1,i2) =  &
                 & D0_Bagger(mSlept2(i3), mSneut2(i4), mC2(i1), mN2(i2) )
      End Do
      Do i4=1,4
       D27_SlSlCN(i3,i4,i1,i2) =  &
                 & D27_Bagger(mSlept2(i3), mSlept2(i4), mC2(i1), mN2(i2) )
      End Do
     End Do
     Do i3=1,2
      Do i4=1,2
        D27_SnSnCN(i3,i4,i1,i2) =  &
                 & D27_Bagger(mSneut2(i3), mSneut2(i4), mC2(i1), mN2(i2) )
      End Do
     End Do
    End Do
   End Do

   B0_SlSn = 0._dp
   D0_SlSnNN = 0._dp
   D27_SlSnNN = 0._dp
   Do i2=1,2
    Do i1=2*i2-1,2*i2
     B0_SlSn(i1,i2) = B0(0._dp, mSlept2(i1), mSneut2(i2) )
     Do i3=1,2
      Do i4=1,2
       D0_SlSnCC(i1,i2,i3,i4) =   &
                 & D0_Bagger(mSlept2(i1), mSneut2(i2), mC2(i3), mC2(i4) )
      End Do
     End Do
     Do i3=1,4
      Do i4=1,4
       D0_SlSnNN(i1,i2,i3,i4) =   &
                 & D0_Bagger(mSlept2(i1), mSneut2(i2), mN2(i3), mN2(i4) )
       D27_SlSnNN(i1,i2,i3,i4) =   &
                 & D27_Bagger(mSlept2(i1), mSneut2(i2), mN2(i3), mN2(i4) )
      End Do
     End Do
    End Do
   End Do

#ifdef GENERATIONMIXING
  End If
#endif

  !------------------------------------------------------------------------
  !  SUSY wave function, sleptons
  !-----------------------------------------------------------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then

   Do i3=1,2
    sumI = 0._dp
    Do i1=1,6
     Do i2=1,4
      sumI = sumI + 0.5_dp * Abs(c_LNSl_R(i3,i2,i1))**2 * B1_N_Sl(i2,i1)
     End Do
    End Do
    If (WriteOut) Write(ErrCan,*) "Z_sle l",i3,sumI
    res = res + SumI
   End Do

   Do i3=1,3
    sumI = 0._dp
    Do i5=1,3
     sumIJ(i5) = 0._dp
     Do i1=1,6
      Do i2=1,2
       sumIJ(i5) = sumIJ(i5) + 0.5_dp * c_CNuSl_R(i2,i3,i1) &
           &                * Conjg(c_CNuSl_R(i2,i5,i1)) * B1_C_Sl(i2,i1)
      End Do
     End Do
     sumI = sumI + sumIJ(i5) * (Unu_L(i5,1)+Unu_L(i5,2))
    End Do ! i5
    If (WriteOut) Write(ErrCan,*) "Z_sle nu",i3,sumI
    res = res + SumI
   End Do

  Else ! .not.GenerationMixing
#endif

   Do i3=1,2
    sumI = 0._dp
    Do i1=2*(i3-1)+1,2*(i3-1)+2
     Do i2=1,4
      sumI = sumI + 0.5_dp * Abs(c_LNSl_R(i3, i2, i1))**2 * B1_N_Sl(i2,i1)
     End Do
    End Do
    If (WriteOut) Write(ErrCan,*) "Z_sle l",i3,sumI
    res = res + sumI
    sumI = 0._dp
    Do i1=2*(i3-1)+1,2*(i3-1)+2
     Do i2=1,2
      sumI = sumI + 0.5_dp * Abs(c_CNuSl_R(i2, i3, i1))**2 * B1_C_Sl(i2,i1)
     End Do
    End Do
    If (WriteOut) Write(ErrCan,*) "Z_sle nu",i3,sumI
    res = res + SumI
   End Do

#ifdef GENERATIONMIXING
  End If
#endif
  !-------------------------------
  !  SUSY wave function, sneutrino
  !-------------------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then

   Do i1=1,2
    sumI = 0._dp
    Do i3=1,3
     Do i2=1,2
      sumI = sumI + 0.5_dp * Abs(c_CLSn_R(i2,i1,i3))**2 * B1_C_Sn(i2,i3)
     End Do
    End Do
    If (WriteOut) Write(ErrCan,*) "Z_snu l",i1,sumI
    res = res + sumI
   End Do

   Do i1=1,3
    sumI = 0._dp
    Do i5=1,3
     sumIJ(i5) = 0._dp
     Do i3=1,3
      Do i2=1,4
       sumIJ(i5) = sumIJ(i5) + 0.5_dp * c_NuNSn_R(i1,i2,i3) &
                               * Conjg(c_NuNSn_R(i5,i2,i3)) * B1_N_Sn(i2,i3)
      End Do
     End Do
     sumI = sumI + sumIJ(i5) * (Unu_L(i5,1)+Unu_L(i5,2))
    End Do
    If (WriteOut) Write(ErrCan,*) "Z_snu nu",i1,sumI 
    res = res + sumI
   End Do

  Else ! .not.GenerationMixing
#endif

   Do i1=1,2
    sumI = 0._dp
    Do i2=1,2
     sumI = sumI + 0.5_dp * Abs(c_CLSn_R(i2,i1,i1))**2 * B1_C_Sn(i2,i1)
    End Do
    If (WriteOut) Write(ErrCan,*) "Z_snu l",i1,sumI
    res = res + sumI

    sumI = 0._dp
    Do i2=1,4
     sumI = sumI + 0.5_dp * Abs(c_NuNSn_R(i1,i2,i1) )**2 * B1_N_Sn(i2,i1)
    End Do
    If (WriteOut) Write(ErrCan,*) "Z_snu nu",i1,sumI
    res = res + sumI
   End Do

#ifdef GENERATIONMIXING
  End If
#endif
  !------------------------------------------------------
  !  vertex, see also 
  ! P.H.Chankowski et al, Nucl. Phys. B417, 101 (1994)
  !------------------------------------------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then

   Do i1=1,2 ! leptons
    Do i5=1,3 ! neutrinos
    sumI = 0._dp
    Do i2=1,4 ! neutralinos
     Do i3=1,2 ! charginos
      Do i4=1,6 ! sleptons
       sumI = sumI                                                            &
          & + c_CNuSl_R(i3,i5,i4) * Conjg( c_LNSl_R(i1,i2,i4) )               &
          &    * ( - sqrt2 * Conjg(c_CNW_R(i3,i2)) * mC(i3) * mN(i2)          &
          &                * C0_SlCN(i4,i3,i2)                                &
          &      + oosqrt2 * Conjg(c_CNW_L(i3,i2))                            &
          &        * (B0_CN(i3,i2) - 0.5_dp + mSlept2(i4) *C0_SlCN(i4,i3,i2)) &
          &      ) / gSU2
      End Do
      Do i4=1,3 ! sneutrinos
       sumI = sumI                                                            &
          & - Conjg(c_CLSn_R(i3,i1,i4)) * c_NuNSn_R(i5,i2,i4)                 &
          &    * ( - sqrt2 * Conjg(c_CNW_L(i3,i2)) * mC(i3) * mN(i2)          &
          &                                         * C0_SnCN(i4,i3,i2)       &
          &      + oosqrt2 * Conjg(c_CNW_R(i3,i2))                            &
          &        * (B0_CN(i3,i2) - 0.5_dp + mSneut2(i4) *C0_SnCN(i4,i3,i2)) &
          &      ) / gSU2
      End Do
     End Do
     Do i3=1,6   ! sleptons
      Do i4=1,3  ! sneutrinos
       RSlSn = 0._dp
       Do i6=1,3
        RSlSn = RSlSn + Conjg(RSlept(i3,i6)) * RSneut(i4,i6)
       End Do
       sumI = sumI                                                            &
          & + 0.5_dp * c_NuNSn_R(i5,i2,i4) * Conjg( c_LNSl_R(i1,i2,i3) )      &
          &  *RSlSn * (B0_SlSn(i3,i4) + 0.5_dp + mN2(i2) * C0_NSlSn(i2,i3,i4) )
      End Do
     End Do
    End Do

    If (WriteOut) Write(ErrCan,*) "v_sl",i1,i5,sumI
    res = res + sumI
    End Do ! i5
   End Do

  Else ! .not.GenerationMixing
#endif

   Do i1=1,2
    sumI = 0._dp
    Do i2=1,4 ! neutralinos
     Do i3=1,2 ! charginos
      Do i4=2*i1-1,2*i1  ! sleptons
       sumI = sumI                                                            &
          & + c_CNuSl_R(i3,i1,i4) * Conjg( c_LNSl_R(i1,i2,i4) )               &
          &    * ( - sqrt2 * Conjg(c_CNW_R(i3,i2)) * mC(i3) * mN(i2)          &
          &                * C0_SlCN(i4,i3,i2)                                &
          &      + oosqrt2 * Conjg(c_CNW_L(i3,i2))                            &
          &        * (B0_CN(i3,i2) - 0.5_dp + mSlept2(i4) *C0_SlCN(i4,i3,i2)) &
          &      ) / gSU2
      End Do
      ! sneutrinos
      sumI = sumI                                                             &
          & - Conjg(c_CLSn_R(i3,i1,i1)) * c_NuNSn_R(i1,i2,i1)                 &
          &    * ( - sqrt2 * Conjg(c_CNW_L(i3,i2)) * mC(i3) * mN(i2)          &
          &                                        * C0_SnCN(i1,i3,i2)        &
          &      + oosqrt2 * Conjg(c_CNW_R(i3,i2))                            &
          &        * (B0_CN(i3,i2) - 0.5_dp + mSneut2(i1) *C0_SnCN(i1,i3,i2)) &
          &      ) / gSU2
     End Do
     Do i3=2*i1-1,2*i1             ! sleptons, sneutrinos
      sumI = sumI                                                            &
         & + 0.5_dp * c_NuNSn_R(i1,i2,i1) * Conjg( c_LNSl_R(i1,i2,i3) )      &
         &   * RSlept(i3,2*i1-1)                                             &
         &   * (B0_SlSn(i3,i1) + 0.5_dp + mN2(i2) * C0_NSlSn(i2,i3,i1) )
     End Do
    End Do
    If (WriteOut) Write(ErrCan,*) "v_sl",i1,sumI
    res = res + sumI
   End Do

#ifdef GENERATIONMIXING
  End If ! GenerationMixing
#endif

  !------------------------------------------------------
  ! Box-diagrams, see also 
  ! P.H.Chankowski et al, Nucl. Phys. B417, 101 (1994)
  !------------------------------------------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   teil = 0._dp 
   Do i1=1,2
    Do i2=1,4
     Do i3=1,6     ! sleptons
      Do i5=1,3   ! summing over neutrinos 1 -> i5, 2 -> i6
       Do i6=1,3
        Do i4=1,3    ! sneutrinos
         teil = teil                   &
            & + 0.5_dp * ( Conjg(c_CLSn_R(i1,2,i4) * c_CNuSl_R(i1,i5,i3) )    &
            &              * c_NuNSn_R(i6,i2,i4) * c_LNSl_R(1,i2,i3)          &
            &            + Conjg(c_LNSl_R(2,i2,i3) * c_NuNSn_R(i5,i2,i4) )    &
            &              * c_CNuSl_R(i1,i6,i3) * c_CLSn_R(i1,1,i4)       )  &
            &   * mC(i1) * mN(i2) * D0_SlSnCN(i3,i4,i1,i2)
        End Do
        Do i4=1,6    ! sleptons
         teil = teil                   &
            & + Conjg(c_LNSl_R(2,i2,i3) * c_CNuSl_R(i1,i5,i4) )    &
            &   * c_CNuSl_R(i1,i6,i3) * c_LNSl_R(1,i2,i4)          &
            &   * D27_SlSlCN(i3,i4,i1,i2)
        End Do ! i4
       End Do ! i6
      End Do ! i5
      
     End Do ! i3
     Do i3=1,3     ! sneutrinos
      Do i4=1,3    ! sneutrinos
       Do i5=1,3   ! summing over neutrinos 1 -> i5, 2 -> i6
        Do i6=1,3
         teil = teil                   &
            & + Conjg(c_CLSn_R(i1,2,i3) * c_NuNSn_R(i5,i2,i4) )    &
            &   * c_NuNSn_R(i6,i2,i3) * c_CLSn_R(i1,1,i4)          &
            &   * D27_SnSnCN(i3,i4,i1,i2)
        End Do ! i6
       End Do ! i5
      End Do
     End Do
    End Do
   End Do

   Do i1=1,6
    Do i2=1,3
     Do i3=1,2
      Do i4=1,2
       Do i5=1,3   ! summing over neutrinos 1 -> i5, 2 -> i6
        Do i6=1,3
         teil = teil &
            & + 0.5_dp * Conjg( c_CLSn_R(i3,2,i2) * c_CNuSl_R(i3,i5,i1) )    &
            &   * c_CNuSl_R(i4,i6,i1) * c_CLSn_R(i4,1,i2) * mC(i3) * mC(i4)  &
            &   * D0_SlSnCC(i1,i2,i3,i4)
        End Do ! i6
       End Do ! i5
      End Do
     End Do
     Do i3=1,4
      Do i4=1,4
       Do i5=1,3   ! summing over neutrinos 1 -> i5, 2 -> i6
        Do i6=1,3
         teil = teil &
            & + 0.5_dp * Conjg( c_LNSl_R(2,i3,i1) * c_NuNSn_R(i5,i3,i2) )    &
            &   * c_NuNSn_R(i6,i4,i2) * c_LNSl_R(1,i4,i1) * mN(i3) * mN(i4)  &
            &   * D0_SlSnNN(i1,i2,i3,i4)                                     &
            & + Conjg( c_LNSl_R(2,i3,i1) * c_NuNSn_R(i5,i4,i2) )             &
            &    * c_NuNSn_R(i6,i3,i2) * c_LNSl_R(1,i4,i1)                   &
            &    * D27_SlSnNN(i1,i2,i3,i4)
        End Do ! i6
       End Do ! i5
      End Do
     End Do
    End Do
   End Do

  Else ! .not.GenerationMixing
#endif

   teil = 0._dp 
   Do i1=1,2
    Do i2=1,4
     Do i3=1,2     ! sleptons
       teil = teil                   &
          & + 0.5_dp * ( Conjg(c_CLSn_R(i1,2,2) * c_CNuSl_R(i1,1,i3) )      &
          &              * c_NuNSn_R(2,i2,2) * c_LNSl_R(1,i2,i3)            &
          &              * D0_SlSnCN(i3,2,i1,i2)                            &
          &            + Conjg(c_LNSl_R(2,i2,2+i3) * c_NuNSn_R(1,i2,1) )    &
          &              * c_CNuSl_R(i1,2,2+i3) * c_CLSn_R(i1,1,1)          &
          &              * D0_SlSnCN(2+i3,1,i1,i2)  ) * mC(i1) * mN(i2)
      Do i4=3,4    ! sleptons
       teil = teil                   &
          & + Conjg(c_LNSl_R(2,i2,i4) * c_CNuSl_R(i1,1,i3) )    &
          &   * c_CNuSl_R(i1,2,i4) * c_LNSl_R(1,i2,i3)          &
          &   * D27_SlSlCN(i4,i3,i1,i2)
      End Do
     End Do
     teil = teil                   &  ! sneutrinos! sneutrinos
        & + Conjg(c_CLSn_R(i1,2,2) * c_NuNSn_R(1,i2,1) )    &
        &   * c_NuNSn_R(2,i2,2) * c_CLSn_R(i1,1,1) * D27_SnSnCN(2,1,i1,i2)
    End Do
   End Do

#ifdef GENERATIONMIXING
  End If ! GenerationMixing
#endif

  sumI =  - 2._dp * cosW2_DR * mZ2 * Real( teil,dp ) / g2
  If (WriteOut) Write(ErrCan,*) "Box contribution",sumI
  res = res + sumI

  res = res  * oo16pi2

  WriteOut = .False.
  Iname = Iname - 1

 End Subroutine delta_VB


 Subroutine delta_VB_SM(gSU2, sinW2, sinW2_DR, rho, res)
 !-----------------------------------------------------------------------
 ! Calculates the the nonuniversal corrections to delta_r in the SM
 ! The formula of J. Bagger et al, Nucl.Phys.B is used. The renormalization
 ! scale has to be set in the main program due to the structure of LoopTools. 
 ! the input is:
 ! - gSU2 ....... the SU(2) gauge coupling at p2
 ! - sinW2 ...... sin(theta_W) squared
 ! - sinW2_DR ... sin(theta_W) squared, DR-scheme
 ! - rho ........ Rho parameter
 ! output
 !  res    
 ! written by Werner Porod, 07.01.09
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: sinW2, sinW2_DR, gSU2, rho
  Real(dp), Intent(out) :: res

  Real(dp) :: cosW2, cosW2_DR

  Iname = Iname + 1
  NameOfUnit(Iname) = "delta_VB_SM"

  cosW2 = 1._dp - sinW2
  cosW2_DR = 1._dp - sinW2_DR

  res = 6._dp                                                       &
     & + Log(cosW2) * (3.5_dp - 2.5_dp * sinW2                       &
     &                - sinW2_DR * (5._dp - 1.5_dp*cosW2/cosW2_DR)   &
     &                ) / sinW2

  res = gSU2**2 * rho * res

  Iname = Iname - 1

 End Subroutine delta_VB_SM



#ifdef GENERATIONMIXING
 Subroutine NeutralinoMass_Loop(gU1, gSU2, Y_d, Y_l, Y_u, vevSM, M1, M2, mu  &
          & , mN, mN2, N, mC, mC2, U, V, mSup2, RSup, mSdown2, RSdown        &
          & , mSlepton2, RSlepton, mSneut2, RSneut, mS02, RS0, mP02, RP0     &
          & , mSpm2, RSpm, uL_L, uL_R, uD_L, uD_R, uU_L, uU_R, mZ2, mW2      &
          & , delta, mN1L, mN1L2, N1L, kont)
#else
 Subroutine NeutralinoMass_Loop(gU1, gSU2, Y_d, Y_l, Y_u, vevSM, M1, M2, mu   &
          & , mN, mN2, N, mC, mC2, U, V, mSup2, RSup, mSdown2, RSdown         &
          & , mSlepton2, RSlepton, mSneut2, mS02, RS0, mP02, RP0, mSpm2, RSpm &
          & , mZ2, mW2, delta, mN1L, mN1L2, N1L, kont)
#endif
 !-----------------------------------------------------------------
 ! calculates neutralino masses + mixing matrix N in the MSSM,
 ! input:
 !  M1 ......... U(1) gaugino mass
 !  M2 ......... SU(2) gaugino mass
 !  mu ......... mu-parameter
 !  vevSM(i) ... i=1 v_d
 !               i=2 v_u
 !  gSU2 ....... SU(2) gauge coupling
 !  gU1 ........ U(1) gauge coupling
 ! written by Werner Porod, 27.12.01
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: vevSM(2), gU1, gSU2, mN(4), mN2(4), mC(2), mC2(2)  &
      &  , mSup2(6), mSdown2(6), mSlepton2(6), mSneut2(3), mS02(2), RS0(2,2) &
      &  , mP02(2), RP0(2,2), mSpm2(2), mZ2, mW2, delta
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: M1, M2, mu, Y_d(3,3), Y_l(3,3), Y_u(3,3)    &
      &  , N(4,4), U(2,2), V(2,2), RSup(6,6), RSdown(6,6), RSlepton(6,6) &
      &  , RSneut(3,3), RSpm(2,2), uL_L(3,3), uL_R(3,3), uD_L(3,3)       &
      &  , uD_R(3,3), uU_L(3,3), uU_R(3,3)
#else
  Complex(dp), Intent(in) :: M1, M2, mu, Y_d(3,3), Y_l(3,3), Y_u(3,3)    &
      &  , N(4,4), U(2,2), V(2,2), RSup(6,6), RSdown(6,6), RSlepton(6,6) &
      &  , RSpm(2,2)
#endif
  Real(dp), Intent(out) :: mN1L(4), mN1L2(4)
  Complex(dp), Intent(out) :: N1L(4,4)

  Integer :: i1, i2, i3, i_l, i_count 
  Complex(dp) :: mat4(4,4), mat42(4,4), phaseM, mat4a(4,4), E4(4) &
     & , SigL(4,4), SigR(4,4), SigS(4,4), c_NNZ_L(4,4), c_NNZ_R(4,4)        &
     & , c_NNS0_L(4,4,2), c_NNS0_R(4,4,2), c_NNP0_L(4,4,2), c_NNP0_R(4,4,2) &
     & , c_UNSu_L(3,4,6), c_UNSu_R(3,4,6), c_LNSl_L(3,4,6), c_LNSl_R(3,4,6) &
     & , c_NuNSn_L(3,4,3), c_NuNSn_R(3,4,3), c_CNW_L(2,4), c_CNW_R(2,4)     &
     & , c_SmpCN_L(2,2,4), c_SmpCN_R(2,2,4), c_DNSd_L(3,4,6), c_DNSd_R(3,4,6)
  Complex(dp) :: Rsl(2,2), Rsu(2,2), Rsd(2,2), coupLC, coupRC, Yuk
  Real(dp) :: g1, gp1, p2, N4a(4,4), work, test(2), cosW, mN_L(4)
  Logical :: WriteOut

  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutralinoMass_Loop'

  !---------------
  ! Initialzation
  !---------------
  g1 = 0.5_dp * gSU2 
  gp1 = 0.5_dp * gU1
  mN1L = 0._dp
  mN1L2 = 0._dp
  N1L = ZeroC
  If ((WriteOneLoopContributions.Eq.11).Or.(WriteOneLoopContributions.Lt.0)) &
   & Then
     Write(ErrCan,*) "Contributions in NeutralinoMass_Loop:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If
  
  !-----------
  ! Tree-level
  !-----------
  mat4a(1,1) = M1
  mat4a(1,2) = 0._dp
  mat4a(1,3) = - gp1 * vevSM(1)
  mat4a(1,4) = gp1 * vevSM(2)
  mat4a(2,2) =  M2
  mat4a(2,3) = g1 * vevSM(1)
  mat4a(2,4) = - g1 * vevSM(2)
  mat4a(3,3) = 0._dp
  mat4a(3,4) = - mu
  mat4a(4,4) = 0._dp

  Do i1=2,4
   Do i2=1,i1-1
    mat4a(i1,i2) = mat4a(i2,i1)
   End Do
  End Do
  If (WriteOut) Then
   Write(ErrCan,*) "Tree level mass matrix:"
   Do i1=1,4
    Write(ErrCan,*) mat4a(i1,:)
   End Do
  End If

  !-----------
  ! Couplings
  !-----------
  c_LNSl_L = ZeroC
  c_LNSl_R = ZeroC
  c_NuNSn_L = ZeroC
  c_NuNSn_R = ZeroC
  c_DNSd_L = ZeroC
  c_DNSd_R = ZeroC
  c_UNSu_L = ZeroC
  c_UNSu_R = ZeroC

#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,4
     Do i3 = 1,6
      Call CoupNeutralinoSlepton(i1, i2, i3, gU1, gSU2, RSlepton, uL_L, uL_R &
                               &, Y_l, id4C, coupLC, coupRC)
      c_LNSl_L(i1, i2, i3 ) = coupLC
      c_LNSl_R(i1, i2, i3 ) = coupRC
      Call CoupNeutralinoSdown(i1, i2, i3, gU1, gSU2, RSdown, uD_L, uD_R &
                               &, Y_d, id4C, coupLC, coupRC)
      c_DNSd_L(i1, i2, i3 ) = coupLC
      c_DNSd_R(i1, i2, i3 ) = coupRC
      Call CoupNeutralinoSup(i1, i2, i3, gU1, gSU2, RSup, uU_L, uU_R &
                               &, Y_u, id4C, coupLC, coupRC)
      c_UNSu_L(i1, i2, i3 ) = coupLC
      c_UNSu_R(i1, i2, i3 ) = coupRC
     End Do
     Do i3 = 1,3
      Call CoupNeutralinoSneutrino(i1, i2, i3, gU1, gSU2, id4C, RSneut, id3C &
                               &, coupRC)
      c_NuNSn_R(i1, i2, i3 ) = coupRC
     End Do
    End Do
   End Do

  Else
#endif
   Do i1 = 1,4
    Call CoupNeutralinoSneutrino(i1, gU1, gSU2, id4C, coupRC)
    Do i2 = 1,3
     c_NuNSn_R(i2, i1, i2 ) = coupRC
    End Do
   End Do

   Do i1 = 1,3
    Rsl = RSlepton(2*i1-1:2*i1, 2*i1-1:2*i1)
    Rsd = RSdown(2*i1-1:2*i1, 2*i1-1:2*i1)
    Rsu = RSup(2*i1-1:2*i1, 2*i1-1:2*i1)

    Yuk = Y_l(i1,i1)
    Do i2 = 1,2
     Do i3 = 1,4
      Call CoupNeutralinoSlepton(i3, i2, gU1, gSU2, RSl, Yuk, id4c &
                               &, coupLC, coupRC)
      c_LNSl_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      c_LNSl_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_d(i1,i1)
    Do i2=1,2
     Do i3=1,4
      Call CoupNeutralinoSdown(i3, i2, gU1, gSU2, RSd, Yuk, id4C, &
                             & coupLC, coupRC)
      c_DNSd_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      c_DNSd_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,4
      Call CoupNeutralinoSup(i3, i2, gU1, gSU2, RSu, Yuk, id4C, &
                           & coupLC, coupRC)
      c_UNSu_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      c_UNSu_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

   End Do ! i1

#ifdef GENERATIONMIXING
  End If
#endif

  c_CNW_L = 0._dp
  c_CNW_R = 0._dp

  Do i1 = 1,n_char
   Do i2= 1,n_neut
    Call CoupCharginoNeutralinoW(i1, i2, id4C, U, V, gSU2               &
                               &, c_CNW_L(i1,i2), c_CNW_R(i1,i2) )
   End Do
  End Do

  c_SmpCN_L = 0._dp
  c_SmpCN_R = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_char
    Do i3 = 1,n_neut
     Call CoupCSCharginoNeutralino(i1, i2, i3, id4C, U, V, RSpm, gU1, gSU2    &
                                 &, c_SmpCN_L(i1,i2,i3), c_SmpCN_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  c_NNZ_L = 0.0_dp
  c_NNZ_R = 0.0_dp
  cosW = gSU2 / Sqrt(gU1**2 + gSU2**2)
  Do i1=1,n_neut
   Do i2=1,n_neut
    Call CoupNeutralinoZa(i1, i2, N, gSU2, cosW  &
                       & , c_NNZ_L(i1,i2), c_NNZ_R(i1,i2) )
   End Do
  End Do

  c_NNP0_L = 0.0_dp
  c_NNP0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_P0
    Call CoupNeutralinoPseudoscalara(i1, i2, i3, N, RP0, gU1, gSU2, &
                       & c_NNP0_L(i1,i2,i3), c_NNP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  c_NNS0_L = 0.0_dp
  c_NNS0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_S0
    Call CoupNeutralinoScalara(i1, i2, i3, N, RS0, gU1, gSU2, &
                       & c_NNS0_L(i1,i2,i3), c_NNS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !--------------------
  ! Loop Contributions
  !--------------------
  Do i_l=n_neut,1,-1
   p2 = mN2(i_l)
   sigL = 0._dp
   sigR = 0._dp
   sigS = 0._dp

   Call Sigma_Neutralino(p2, mN, mN2, c_NNZ_L, c_NNZ_R, mS02, c_NNS0_L    &
     , c_NNS0_R, mP02, c_NNP0_L, c_NNP0_R, mC, mC2, c_CNW_L, c_CNW_R      &
     , mSpm2, c_SmpCN_L, c_SmpCN_R, mSup2, c_UNSu_L, c_UNSu_R             &
     , mSdown2, c_DNSd_L, c_DNSd_R, mSlepton2, c_LNSl_L, c_LNSl_R         &
     , mSneut2, c_NuNSn_L, c_NuNSn_R, mZ2, mW2, WriteOut, SigL, SigR, SigS)

   If (WriteOut) Then
       Write(ErrCan,*) "Neut, sigL",i_l
       Do i1=1,4
        Write(ErrCan,*) (sigL(i1,i2),i2=1,4)
       End Do
       Write(ErrCan,*) "Neut, sigR",i_l
       Do i1=1,4
        Write(ErrCan,*) (sigr(i1,i2),i2=1,4)
       End Do
       Write(ErrCan,*) "Neut, sigS",i_l
       Do i1=1,4
        Write(ErrCan,*) (sigs(i1,i2),i2=1,4)
       End Do
   End If

   mat4 = mat4a - 0.5_dp * ( SigS + Transpose(SigS)          &
                           + Matmul(Transpose(SigL), mat4a)  &
                           + Matmul(SigR, mat4a)             &
                           + Matmul(mat4a,Transpose(SigR))   &
                           + Matmul(mat4a, SigL) )

   If (Maxval(Abs(Aimag(mat4))).Eq.0._dp) Then ! matrix is reel
    Call EigenSystem(Real(mat4,dp), mN_L, N4a, kont, test)

    Do i1=1,4
     If (mN_L(i1).Lt.0._dp) Then
      mN_L(i1) = - mN_L(i1)
      N1L(i1,:) = (0._dp,1._dp) * N4a(i1,:)
     Else
      N1L(i1,:) =N4a(i1,:)
     End If
    End Do

    Do i1=1,3
     Do i2=i1+1,4
      If (Abs(mN_L(i1)).Gt.Abs(mN_L(i2))) Then
       work = mN_L(i1)
       mN_L(i1) = mN_L(i2)
       mN_L(i2) = work
       E4 = N1L(i1,:)
       N1L(i1,:) = N1L(i2,:)
       N1L(i2,:) = E4
      End If
     End Do
    End Do
    mN1L(i_L) = mN_L(i_L)
    mN1L2(i_L) = mN1L(i_L)**2

   Else

    mat42 = Matmul( Transpose(Conjg( mat4 ) ), mat4 )
    Call EigenSystem(mat42, mN_L, N1L, kont, test)

    mat42 = Matmul(Conjg(N1L), Matmul( mat4, Transpose( Conjg( N1L ) ) ) )
    Do i1=1,4
     phaseM =   Sqrt( mat42(i1,i1)   / Abs( mat42(i1,i1) ) )
     N1L(i1,:) = phaseM * N1L(i1,:)
    End Do
    mN1L2(i_L) = mN_L(i_L)
    mN1L(i_L) = Sqrt( mN1L2(i_L) )

   End If

   If ((kont.Eq.-14).Or.(kont.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(ErrCan,*) "test =",test
     Write(ErrCan,*) " "
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = 0
   End If


   !------------------------------------------
   ! redoing calculation using refined p2
   !------------------------------------------
   i_count = 0
   p2_loop: Do 
    i_count = i_count + 1 
    p2 = mN1L2(i_l)
    sigL = 0._dp
    sigR = 0._dp
    sigS = 0._dp

    Call Sigma_Neutralino(p2, mN, mN2, c_NNZ_L, c_NNZ_R, mS02, c_NNS0_L   &
     , c_NNS0_R, mP02, c_NNP0_L, c_NNP0_R, mC, mC2, c_CNW_L, c_CNW_R      &
     , mSpm2, c_SmpCN_L, c_SmpCN_R, mSup2, c_UNSu_L, c_UNSu_R             &
     , mSdown2, c_DNSd_L, c_DNSd_R, mSlepton2, c_LNSl_L, c_LNSl_R         &
     , mSneut2, c_NuNSn_L, c_NuNSn_R, mZ2, mW2, WriteOut, SigL, SigR, SigS)

    If (WriteOut) Then
       Write(ErrCan,*) "Neut, sigL",i_l
       Do i1=1,4
        Write(ErrCan,*) (sigL(i1,i2),i2=1,4)
       End Do
       Write(ErrCan,*) "Neut, sigR",i_l
       Do i1=1,4
        Write(ErrCan,*) (sigr(i1,i2),i2=1,4)
       End Do
       Write(ErrCan,*) "Neut, sigS",i_l
       Do i1=1,4
        Write(ErrCan,*) (sigs(i1,i2),i2=1,4)
       End Do
    End If

    mat4 = mat4a - 0.5_dp * ( SigS + Transpose(SigS)         &
                           + Matmul(Transpose(SigL), mat4a)  &
                           + Matmul(SigR, mat4a)             &
                           + Matmul(mat4a,Transpose(SigR))   &
                           + Matmul(mat4a, SigL) )
    Call chop(mat4) ! to avoid problems with tiny numbers

    If (Maxval(Abs(Aimag(mat4))).Eq.0._dp) Then ! matrix is reel
     Call EigenSystem(Real(mat4,dp), mN_L, N4a, kont, test)

     Do i1=1,4
      If (mN_L(i1).Lt.0._dp) Then
       mN_L(i1) = - mN_L(i1)
       N1L(i1,:) = (0._dp,1._dp) * N4a(i1,:)
      Else
       N1L(i1,:) =N4a(i1,:)
      End If
     End Do

     Do i1=1,3
      Do i2=i1+1,4
       If (Abs(mN_L(i1)).Gt.Abs(mN_L(i2))) Then
        work = mN_L(i1)
        mN_L(i1) = mN_L(i2)
        mN_L(i2) = work
        E4 = N1L(i1,:)
        N1L(i1,:) = N1L(i2,:)
        N1L(i2,:) = E4
       End If
      End Do
     End Do
     mN1L(i_L) = mN_L(i_L)
     mN1L2(i_L) = mN1L(i_L)**2

    Else

     mat42 = Matmul( Transpose(Conjg( mat4 ) ), mat4 )
     Call EigenSystem(mat42, mN_L, N1L, kont, test)

     mat42 = Matmul(Conjg(N1L), Matmul( mat4, Transpose( Conjg( N1L ) ) ) )
     Do i1=1,4
      phaseM =   Sqrt( mat42(i1,i1)   / Abs( mat42(i1,i1) ) )
      N1L(i1,:) = phaseM * N1L(i1,:)
     End Do
     mN1L2(i_L) = mN_L(i_L)
     mN1L(i_L) = Sqrt( mN1L2(i_L) )

    End If

    If ((kont.Ne.0).And.(ErrorLevel.Gt.-1)) Then
     Write(ErrCan,*) 'Warning in subroutine NeutralinoMass_Loop, ierr =',kont
     Write(ErrCan,*) 'M1,M2 ',M1,M2
     Write(ErrCan,*) 'gp,g ',gU1, gSU2
     Write(ErrCan,*) 'mu ',mu
     Write(ErrCan,*) 'vevs ',vevSM
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If

    If (p2.Ne.0._dp) Then
     test(1) = Abs(mN1L2(i_l) - p2) / p2
    Else
     test(1) = Abs(mN1L2(i_l))
    End If
    If (test(1).Lt.0.1_dp*delta) Exit p2_loop
    If (i_count.Gt.30) Then
     Write(*,*) "Problem in loop for mN",i_l,test(1),mN1L2(i_l)
     Exit p2_loop
    End If

   End Do p2_loop
   
  End Do !  i_l

  Iname = Iname - 1

 End Subroutine NeutralinoMass_Loop


 Subroutine NeutralinoMass_Loop_RP(gU1, gSU2, Y_d, Y_l, Y_u, vevSM, vevL     &
          & , M1, M2, mu, eps, mC, mC2, U, V, mSup2, RSup, mSdown2, RSdown   &
          & , mS02, RS0, mP02, RP0, mSpm2, RSpm, uD_L, uD_R, uU_L, uU_R      &
          & , mN, mN2, N, mN1L, mN1L2, N1L, kont, WriteOut)
 !-----------------------------------------------------------------
 ! calculates neutralino masses + mixing matrix N in the MSSM,
 ! 1-generation epsilon model and 3-generation epsilon model. In
 ! addition the lepton yukawas are calculated.
 ! input:
 !  n_neut .... specifies the model: n_neut = 2 -> MSSM
 !                                   n_neut = 3 -> 1-generation epsilon model
 !                                   n_neut = 5 -> 3-generation epsilon model
 !  M1 ........ U(1) gaugino mass
 !  M2 ........ SU(2) gaugino mass
 !  bi(i) ..... i=1 mu-parameter
 !              i=2 epsilon_1
 !              i=3 epsilon_2
 !              i=4 epsilon_3
 !  vevs(i) ... i=1 v_d
 !              i=2 v_u
 !              i=3 v_1
 !              i=4 v_2
 !              i=5 v_3
 !  g ......... SU(2) gauge coupling
 !  gp ........ U(1) gauge coupling
 ! written by Werner Porod, 27.12.01
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: vevSM(2), gU1, gSU2, mC(5), mC2(5)  &
      &  , mSup2(6), mSdown2(6), mS02(5), RS0(5,5) &
      &  , mP02(5), RP0(5,5), mSpm2(8), vevL(3)
  Complex(dp), Intent(in) :: M1, M2, mu, Y_d(3,3), Y_l(3,3), Y_u(3,3)     &
      &  , U(5,5), V(5,5), RSup(6,6), RSdown(6,6), RSpm(8,8)              &
      &  , uD_L(3,3), uD_R(3,3), uU_L(3,3), uU_R(3,3), eps(3)
  Logical, Intent(in) :: WriteOut
  Real(dp), Intent(out) :: mN1L(7), mN1L2(7), mN(7), mN2(7)
  Complex(dp), Intent(out) :: N1L(7,7), N(7,7)

  Integer :: i1, i2, i3, n_neut, n_char, n_S0, n_P0, n_Spm
  Complex(dp) :: mat7(7,7), mat72(7,7), test7(7,7), phaseM, mat7a(7,7), E7(7) &
     & , SigS(7,7), c_NNZ_L(7,7), c_NNZ_R(7,7)                                &
     & , c_NNS0_L(7,7,5), c_NNS0_R(7,7,5), c_NNP0_L(7,7,5), c_NNP0_R(7,7,5)   &
     & , c_UNSu_L(3,7,6), c_UNSu_R(3,7,6), c_CNW_L(5,7), c_CNW_R(5,7)         &
     & , c_SmpCN_L(8,5,7), c_SmpCN_R(8,5,7), c_DNSd_L(3,7,6), c_DNSd_R(3,7,6)
  Complex(dp) :: Rsu(2,2), Rsd(2,2), coupLC, coupRC, Yuk, Ylep(3)
  Real(dp) :: g1, gp1, N7a(7,7), work, test(2), cosW
  ! for approximate diagonalization
  Complex(dp) :: mat4(4,4), N4(4,4), Xi(3,4), mat42(4,4), mat3(3,3), mat32(3,3) &
    & , Vnu(3,3), Minv(4,4), R1(7,7), R2(7,7), XiR(3,4)
  Real(dp) :: mN4(4),  N4R(4,4), MinvR(4,4), E4R(4), mat3R(3,3), mnu(3)  &
    & , VnuR(3,3) , mat4R(4,4)
  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutralinoMass_Loop_RP'

  !---------------
  ! Initialzation
  !---------------
  g1 = 0.5_dp * gSU2 
  gp1 = 0.5_dp * gU1
  mN1L = 0._dp
  mN1L2 = 0._dp
  N1L = ZeroC
  If (WriteOut)  Write(ErrCan,*) "Contributions in NeutralinoMass_Loop_RP:"

  n_neut = 7
  n_char = 5
  n_P0 = 5
  n_S0 = 5
  n_Spm = 8
  Do i1=1,3
   Ylep(i1) = Y_l(i1,i1)
  End Do 

  !-----------
  ! Tree-level
  !-----------
  mat7a  = 0._dp
  mat7a(1,1) = M1
  mat7a(1,3) = - gp1 * vevSM(1)
  mat7a(1,4) = gp1 * vevSM(2)
  mat7a(1,5) = - gp1 * vevL(1)
  mat7a(1,6) = - gp1 * vevL(2)
  mat7a(1,7) = - gp1 * vevL(3)
  mat7a(2,2) =  M2
  mat7a(2,3) = g1 * vevSM(1)
  mat7a(2,4) = - g1 * vevSM(2)
  mat7a(2,5) = g1 * vevL(1)
  mat7a(2,6) = g1 * vevL(2)
  mat7a(2,7) = g1 * vevL(3)
  mat7a(3,4) = - mu
  mat7a(4,5) = eps(1)
  mat7a(4,6) = eps(2)
  mat7a(4,7) = eps(3)

  Do i1=2,7
   Do i2=1,i1-1
    mat7a(i1,i2) = mat7a(i2,i1)
   Enddo
  Enddo

  If (WriteOut) Then
   Write(ErrCan,*) "Tree level mass matrix:"
   Do i1=1,7
    Write(ErrCan,*) mat7a(i1,:)
   End Do
  End If

  If (Maxval(Abs(Aimag(mat7a))).Eq.0._dp) Then ! matrix is reel
   Call EigenSystemQP(Real(mat7a,dp), mN, N7a, kont, test)

   N = N7a
   Do i1=1,6
    Do i2=i1+1,7
     If (Abs(mN(i1)).Gt.Abs(mN(i2))) Then
      work = mN(i1)
      mN(i1) = mN(i2)
      mN(i2) = work
      E7 = N(i1,:)
      N(i1,:) = N(i2,:)
      N(i2,:) = E7
     End If
    End Do
   End Do
   mN2 = mN**2

  Else

   Call EigenSystemQP(mat7a, mN, N, kont, test)
   mat72 = Matmul(Conjg(N), Matmul( mat7a, Transpose( Conjg( N ) ) ) )
   Do i1=1,7
    phaseM =   Sqrt( mat72(i1,i1)   / Abs( mat72(i1,i1) ) )
    N(i1,:) = phaseM * N(i1,:)
   End Do
   mN2 = mN**2 

  End If

  !-----------
  ! Couplings
  !-----------
  c_DNSd_L = ZeroC
  c_DNSd_R = ZeroC
  c_UNSu_L = ZeroC
  c_UNSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,7
     Do i3 = 1,6
      Call CoupNeutralinoSdown(i1, i2, i3, gU1, gSU2, RSdown, uD_L, uD_R &
                               &, Y_d, N, coupLC, coupRC)
      c_DNSd_L(i1, i2, i3 ) = coupLC
      c_DNSd_R(i1, i2, i3 ) = coupRC
      Call CoupNeutralinoSup(i1, i2, i3, gU1, gSU2, RSup, uU_L, uU_R &
                               &, Y_u, N, coupLC, coupRC)
      c_UNSu_L(i1, i2, i3 ) = coupLC
      c_UNSu_R(i1, i2, i3 ) = coupRC
     End Do
    End Do
   End Do

  Else
   Do i1 = 1,3
    Rsd = RSdown(2*i1-1:2*i1, 2*i1-1:2*i1)
    Rsu = RSup(2*i1-1:2*i1, 2*i1-1:2*i1)

    Yuk = Y_d(i1,i1)
    Do i2=1,2
     Do i3=1,7
      Call CoupNeutralinoSdown(i3, i2, gU1, gSU2, RSd, Yuk, N,coupLC,coupRC)
      c_DNSd_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      c_DNSd_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,7
      Call CoupNeutralinoSup(i3, i2, gU1, gSU2, RSu, Yuk, N, coupLC, coupRC)
      c_UNSu_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      c_UNSu_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

   End Do ! i1

  End If

  c_CNW_L = 0._dp
  c_CNW_R = 0._dp

  Do i1 = 1,n_char
   Do i2= 1,n_neut
    Call CoupCharginoNeutralinoW(i1, i2, N, U, V, gSU2               &
                               &, c_CNW_L(i1,i2), c_CNW_R(i1,i2) )
   End Do
  End Do

  c_SmpCN_L = 0._dp
  c_SmpCN_R = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_char
    Do i3 = 1,n_neut
     Call CoupCSCharginoNeutralino(i1, i2, i3, N, U, V, RSpm, Ylep     &
                  &, gU1, gSU2, c_SmpCN_L(i1,i2,i3), c_SmpCN_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  c_NNZ_L = 0.0_dp
  c_NNZ_R = 0.0_dp
  cosW = gSU2 / Sqrt(gU1**2 + gSU2**2)
  Do i1=1,n_neut
   Do i2=1,n_neut
    Call CoupNeutralinoZ(i1, i2, N, gSU2, cosW, c_NNZ_L(i1,i2), c_NNZ_R(i1,i2) )
   End Do
  End Do

  c_NNP0_L = 0.0_dp
  c_NNP0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_P0
     Call CoupNeutralinoPseudoscalar(i1, i2, i3, N, RP0, gU1, gSU2, &
                                  & c_NNP0_L(i1,i2,i3), c_NNP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  c_NNS0_L = 0.0_dp
  c_NNS0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_S0
     Call CoupNeutralinoScalar(i1, i2, i3, N, RS0, gU1, gSU2, &
                             & c_NNS0_L(i1,i2,i3), c_NNS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !--------------------
  ! Loop Contributions
  !--------------------
! LoopContributions(6:7) = .false.

  Call Sigma_Neutralino_2(mN, mN2, c_NNZ_L, c_NNZ_R, mS02, c_NNS0_L       &
        & , c_NNS0_R, mP02, c_NNP0_L, c_NNP0_R, mC, mC2, c_CNW_L, c_CNW_R &
        & , mSpm2, c_SmpCN_L, c_SmpCN_R, mSup2, c_UNSu_L, c_UNSu_R        &
        & , mSdown2, c_DNSd_L, c_DNSd_R, WriteOut, SigS)

  If (WriteOut) Then
   Write(ErrCan,*) "Neut, sigS"
   Do i1=1,7
    If (Maxval(Abs(Aimag(sigs(i1,:)))).Le.Epsilon(1._dp)) Then
     Write(ErrCan,*) (Real(sigs(i1,i2),dp),i2=1,7)
    Else
     Write(ErrCan,*) (sigs(i1,i2),i2=1,7)
    End If
   End Do
  End If
  mat7 = SigS
  Do i1=1,7
   mat7(i1,i1) = mN(i1) + SigS(i1,i1)
  End Do

  If (Maxval(Abs(Aimag(mat7))).Eq.0._dp) Then ! matrix is reel
   Call EigenSystemQP(Real(mat7,dp), mN1L, N7a, kont, test)

   N1L = N7a

   Do i1=1,6
    Do i2=i1+1,7
     If (Abs(mN1L(i1)).Gt.Abs(mN1L(i2))) Then
      work = mN1L(i1)
      mN1L(i1) = mN1L(i2)
      mN1L(i2) = work
      E7 = N1L(i1,:)
      N1L(i1,:) = N1L(i2,:)
      N1L(i2,:) = E7
     End If
    End Do
   End Do
   mN1L2 = mN1L**2

!-------------------------------------------------------------------------
! there is a huge hierarchy between neutrinos and neutralinos
! therefore lets try seesaw approximation if the first attempt fails
!-------------------------------------------------------------------------
   If (kont.Eq.-1006) Then
    kont = 0
    mat4R = Real(mat7(4:7,4:7) ,dp)
    Call EigenSystem(mat4R, mN4, N4R, kont, test)
    Do i1=1,4
     Do i2=i1+1,4
      If (Abs(mN4(i1)).Gt.Abs(mN4(i2))) Then
       work = mN4(i1)
       mN4(i1) = mN4(i2)
       mN4(i2) = work
       E4R = N4R(i1,:)
       N4R(i1,:) = N4R(i2,:)
       N4R(i2,:) = E4R
      End If
     End Do
    End Do

    XiR = Real(mat7(1:3,4:7),dp)
    MinvR = 0._dp
    Do i1=1,4
     Do i2=1,4
      Do i3=1,4
       MinvR(i1,i2) = MinvR(i1,i2) + N4R(i3,i1) * N4R(i3,i2)  / mN4(i3)
      End Do
     End Do
    End Do  
    mat3R = mat7(1:3,1:3) - Matmul( XiR, Matmul(MinvR, Transpose(XiR) ) )
    Call EigenSystem(mat3R, mnu, VnuR, kont, test)
    Do i1=1,3
     Do i2=i1+1,3
      If (Abs(mNu(i1)).Gt.Abs(mNu(i2))) Then
       work = mNu(i1)
       mNu(i1) = mNu(i2)
       mNu(i2) = work
       E4R(1:3) = VNUR(i1,:)
       VNUR(i1,:) = VNUR(i2,:)
       VNUR(i2,:) = E4R(1:3)
      End If
     End Do
    End Do

    R1 = 0._dp
    R1(1:3,1:3) = VnuR
    R1(4:7,4:7) = N4R 
    XiR = Matmul(XiR,MinvR)
    R2(1:3,1:3) = - 0.5_dp * Matmul(XiR, Transpose(XiR))
    R2(4:7,4:7) = - 0.5_dp * Matmul(Transpose(XiR), XiR)
    Do i1=1,7
     R2(i1,i1) = 1._dp + R2(i1,i1)
    End Do
    R2(1:3,4:7) = - XiR
    R2(4:7,1:3) = Transpose(XiR)
    
    N1L = Matmul(R1,R2)
    mn1l(1:3) = mnu
    mn1l(4:7) = mn4
    mN1L2 = mN1L**2
   End If  ! seesaw approximation
 
  Else
   mat72 = Matmul( Transpose(Conjg( mat7)),  mat7)
   Call EigenSystemQP(mat72, mN1L2, N1L, kont, test)
   mat72 = Matmul(Conjg(N1L), Matmul( mat7, Transpose( Conjg( N1L ) ) ) )
   Do i1=1,7
    phaseM =   Sqrt( mat72(i1,i1)   / Abs( mat72(i1,i1) ) )
    N1L(i1,:) = phaseM * N1L(i1,:)
   End Do
   mN1L = Sqrt( mN1L2 )
!-------------------------------------------------------------------------
! there is a huge hierarchy between neutrinos and neutralinos
! therefore lets try seesaw approximation if the first attempt fails
!-------------------------------------------------------------------------
   If (kont.Eq.-1006) Then
    kont = 0
    mat4 = mat7(4:7,4:7) 
    mat42 = Matmul( Transpose(Conjg( mat4)),  mat4)
    Call EigenSystem(mat42, mN4, N4, kont, test)
    mat42 = Matmul(Conjg(N4), Matmul( mat4, Transpose( Conjg( N4 ) ) ) )
    Do i1=1,4
     phaseM =   Sqrt( mat42(i1,i1)   / Abs( mat42(i1,i1) ) )
     N4(i1,:) = phaseM * N4(i1,:)
    End Do
    mN4 = Sqrt( mN4)

    Xi = mat7(1:3,4:7)
    Minv = 0._dp
    Do i1=1,4
     Do i2=1,4
      Do i3=1,4
       Minv(i1,i2) = Minv(i1,i2) + Conjg(N4(i3,i1) * N4(i3,i2) ) / mN4(i3)
      End Do
     End Do
    End Do  
    mat3 = mat7(1:3,1:3) - Matmul( Xi, Matmul(Minv, Transpose(Xi) ) )

    mat32 = Matmul( Transpose(Conjg( mat3)),  mat3)
    Call EigenSystem(mat32, mnu, Vnu, kont, test)
    mat32 = Matmul(Conjg(Vnu), Matmul( mat3, Transpose( Conjg( Vnu) ) ) )
    Do i1=1,3
     phaseM =   Sqrt( mat32(i1,i1)   / Abs( mat32(i1,i1) ) )
     Vnu(i1,:) = phaseM * Vnu(i1,:)
    End Do
    mNu = Sqrt( mNu)

    R1 = 0._dp
    R1(1:3,1:3) = Conjg(Vnu)
    R1(4:7,4:7) = Conjg(N4) 
    Xi = Matmul(Xi,Minv)
    R2(1:3,1:3) = - 0.5_dp * Matmul(Xi, Transpose(Conjg(Xi)))
    R2(4:7,4:7) = - 0.5_dp * Matmul(Transpose(Conjg(Xi)), Xi)
    Do i1=1,7
     R2(i1,i1) = 1._dp + R2(i1,i1)
    End Do
    R2(1:3,4:7) = - Xi
    R2(4:7,1:3) = Transpose(Conjg(Xi))
    
    N1L = Matmul(R1,R2)
    mn1l(1:3) = mnu
    mn1l(4:7) = mn4
    mN1l2 = mN1L**2
   End If  ! seesaw approximation
 
  End If

  N1L = Matmul(N1L, N)

  Do i1=1,7
   If (mN1L(i1).Lt.0._dp) Then
    mN1L(i1) = -mN1L(i1)
    N1L(i1,:) = N1L(i1,:) * (0._dp,1._dp)
   End If
  End Do

  If ((kont.Ne.0).And.(ErrorLevel.Gt.-1)) Then
   Write (ErrCan,*) 'Warning in subroutine NeutralinoMass_Loop_RP, ierr =',kont
   Write(Errcan,*) "test",test
   Write (ErrCan,*) 'M1,M2 ',M1,M2
   Write (ErrCan,*) 'gp,g ',gU1, gSU2
   Write (ErrCan,*) 'mu ',mu
   Write (ErrCan,*) 'vevs ',vevSM
   If (ErrorLevel.Eq.2) Call TerminateProgram
  End If

  If (WriteOut) Then
   Write(ErrCan,*) " "
   Write(ErrCan,*) "Test of Diagonalization, resulting matrix is:"
   Write(ErrCan,*) "mN",mN1L
   Call adjungate(N1L,test7)
   test7 = Matmul( Transpose(test7), Matmul(mat7, test7) )
   Do i1=1,7
    Write(ErrCan,*) test7(i1,:)
   End Do
   Write(ErrCan,*) " "
  End If
  Iname = Iname - 1

 End Subroutine NeutralinoMass_Loop_RP


 Subroutine NeutralinoMass_Loop_RPtri(gU1, gSU2, Y_d, Y_l, Y_u, RP_lam          &
          & , RP_lamp, vevSM, vevL, M1, M2, mu, eps, mC, mC2, U, V, mSup2, RSup &
          & , mSdown2, RSdown, mS02, RS0, mP02, RP0, mSpm2, RSpm, uD_L, uD_R    &
          & , uU_L, uU_R, mN, mN2, N, mN1L, mN1L2, N1L, kont, WriteOut)
 !-----------------------------------------------------------------
 ! calculates neutralino masses + mixing matrix N in the MSSM,
 ! 1-generation epsilon model and 3-generation epsilon model. In
 ! addition the lepton yukawas are calculated.
 ! input:
 !  n_neut .... specifies the model: n_neut = 2 -> MSSM
 !                                   n_neut = 3 -> 1-generation epsilon model
 !                                   n_neut = 5 -> 3-generation epsilon model
 !  M1 ........ U(1) gaugino mass
 !  M2 ........ SU(2) gaugino mass
 !  bi(i) ..... i=1 mu-parameter
 !              i=2 epsilon_1
 !              i=3 epsilon_2
 !              i=4 epsilon_3
 !  vevs(i) ... i=1 v_d
 !              i=2 v_u
 !              i=3 v_1
 !              i=4 v_2
 !              i=5 v_3
 !  g ......... SU(2) gauge coupling
 !  gp ........ U(1) gauge coupling
 ! written by Werner Porod, 27.12.01
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: vevSM(2), gU1, gSU2, mC(5), mC2(5)  &
      &  , mSup2(6), mSdown2(6), mS02(5), RS0(5,5) &
      &  , mP02(5), RP0(5,5), mSpm2(8), vevL(3)
  Complex(dp), Intent(in) :: M1, M2, mu, Y_d(3,3), Y_l(3,3), Y_u(3,3)  &
      &  , U(5,5), V(5,5), RSup(6,6), RSdown(6,6), RSpm(8,8)           &
      &  , uD_L(3,3), uD_R(3,3), uU_L(3,3), uU_R(3,3), eps(3)          &
      &  , RP_lam(3,3,3), RP_lamp(3,3,3)
  Logical, Intent(in) :: WriteOut
  Real(dp), Intent(out) :: mN1L(7), mN1L2(7), mN(7), mN2(7)
  Complex(dp), Intent(out) :: N1L(7,7), N(7,7)

  Integer :: i1, i2, i3, n_neut, n_char, n_S0, n_P0, n_Spm,ii
  Complex(dp) :: mat7(7,7), mat72(7,7), test7(7,7), phaseM, mat7a(7,7), E7(7) &
     & , SigS(7,7), c_NNZ_L(7,7), c_NNZ_R(7,7)                                &
     & , c_NNS0_L(7,7,5), c_NNS0_R(7,7,5), c_NNP0_L(7,7,5), c_NNP0_R(7,7,5)   &
     & , c_UNSu_L(3,7,6), c_UNSu_R(3,7,6), c_CNW_L(5,7), c_CNW_R(5,7)         &
     & , c_SmpCN_L(8,5,7), c_SmpCN_R(8,5,7), c_DNSd_L(3,7,6), c_DNSd_R(3,7,6)
  Complex(dp) :: Rsu(2,2), coupLC, coupRC, Yuk
  Real(dp) :: g1, gp1, N7a(7,7), work, test(2), cosW
  ! for approximate diagonalization
  Complex(dp) :: mat4(4,4), N4(4,4), Xi(3,4), mat42(4,4), mat3(3,3), mat32(3,3) &
    & , Vnu(3,3), Minv(4,4), RSd6(6,6)
  Real(dp) :: mN4(4),  N4R(4,4), XiR(3,4), MinvR(4,4), E4R(4), mat3R(3,3), mnu(3)  &
    & , VnuR(3,3) , mat4R(4,4), R1(7,7), R2(7,7)
  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutralinoMass_Loop_RPtri'

  !---------------
  ! Initialzation
  !---------------
  g1 = 0.5_dp * gSU2 
  gp1 = 0.5_dp * gU1
  mN1L = 0._dp
  mN1L2 = 0._dp
  N1L = ZeroC
  If (WriteOut)  Write(ErrCan,*) "Contributions in NeutralinoMass_Loop_RPtri:"

  n_neut = 7
  n_char = 5
  n_P0 = 5
  n_S0 = 5
  n_Spm = 8

  !-----------
  ! Tree-level
  !-----------
  mat7a  = 0._dp
  mat7a(1,1) = M1
  mat7a(1,3) = - gp1 * vevSM(1)
  mat7a(1,4) = gp1 * vevSM(2)
  mat7a(1,5) = - gp1 * vevL(1)
  mat7a(1,6) = - gp1 * vevL(2)
  mat7a(1,7) = - gp1 * vevL(3)
  mat7a(2,2) =  M2
  mat7a(2,3) = g1 * vevSM(1)
  mat7a(2,4) = - g1 * vevSM(2)
  mat7a(2,5) = g1 * vevL(1)
  mat7a(2,6) = g1 * vevL(2)
  mat7a(2,7) = g1 * vevL(3)
  mat7a(3,4) = - mu
  mat7a(4,5) = eps(1)
  mat7a(4,6) = eps(2)
  mat7a(4,7) = eps(3)

  Do i1=2,7
   Do i2=1,i1-1
    mat7a(i1,i2) = mat7a(i2,i1)
   Enddo
  Enddo

  If (WriteOut) Then
   Write(ErrCan,*) "Tree level mass matrix:"
   Do i1=1,7
    Write(ErrCan,*) mat7a(i1,:)
   End Do
  End If

  If (Maxval(Abs(Aimag(mat7a))).Eq.0._dp) Then ! matrix is reel
   Call EigenSystemQP(Real(mat7a,dp), mN, N7a, kont, test)

   N = N7a
   Do i1=1,6
    Do i2=i1+1,7
     If (Abs(mN(i1)).Gt.Abs(mN(i2))) Then
      work = mN(i1)
      mN(i1) = mN(i2)
      mN(i2) = work
      E7 = N(i1,:)
      N(i1,:) = N(i2,:)
      N(i2,:) = E7
     End If
    End Do
   End Do
   mN2 = mN**2

  Else

   Call EigenSystemQP(mat7a, mN, N, kont, test)
   mat72 = Matmul(Conjg(N), Matmul( mat7a, Transpose( Conjg( N ) ) ) )
   Do i1=1,7
    phaseM =   Sqrt( mat72(i1,i1)   / Abs( mat72(i1,i1) ) )
    N(i1,:) = phaseM * N(i1,:)
   End Do
   mN2 = mN**2 

  End If

  !-----------
  ! Couplings
  !-----------
  c_DNSd_L = ZeroC
  c_DNSd_R = ZeroC
  c_UNSu_L = ZeroC
  c_UNSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,7
     Do i3 = 1,6
      Call CoupNeutralinoSdown(i1, i2, i3, gU1, gSU2, RSdown, uD_L, uD_R &
                               &, Y_d, N, RP_lamp, coupLC, coupRC)
      c_DNSd_L(i1, i2, i3 ) = coupLC
      c_DNSd_R(i1, i2, i3 ) = coupRC
      Call CoupNeutralinoSup(i1, i2, i3, gU1, gSU2, RSup, uU_L, uU_R &
                               &, Y_u, N, coupLC, coupRC)
      c_UNSu_L(i1, i2, i3 ) = coupLC
      c_UNSu_R(i1, i2, i3 ) = coupRC
     End Do
    End Do
   End Do

  Else
   !-----------------------------------------------------------
   ! re-shuffle d-squark mixing matrix
   !-----------------------------------------------------------
   Rsd6 = 0
   Do i1=1,3
    Rsd6(i1,i1) = RSdown(1+(i1-1)*2,1+(i1-1)*2)
    Rsd6(i1,i1+3) = RSdown(1+(i1-1)*2,2+(i1-1)*2)
    Rsd6(i1+3,i1) = RSdown(2+(i1-1)*2,1+(i1-1)*2)
    Rsd6(i1+3,i1+3) = RSdown(2+(i1-1)*2,2+(i1-1)*2)
   End Do
   Do i1 = 1,3
    Do i2 = 1,7
     Do i3 = 1,6
      If ((i1.Eq.i3).Or.(i1+3.Eq.i3)) Then
       Call CoupNeutralinoSdown(i1, i2, i3, gU1, gSU2, RSd6, id3C, id3C &
                                &, Y_d, N, RP_lamp, coupLC, coupRC)
       !-------------------------------------------------
       ! and now the couplings back to the orginal basis
       !-------------------------------------------------
       Select Case(i3)
       Case(2)
        ii = 3
       Case(3)
        ii = 5
       Case(4)
        ii = 2
       Case(5)
        ii = 4
       Case Default
        ii = i3
       End Select
       c_DNSd_L(i1, i2, ii ) = coupLC
       c_DNSd_R(i1, i2, ii ) = coupRC
      End If
     End Do
    End Do
   End Do
   !--------------------------
   ! u-squarks are easier
   !--------------------------
   Do i1 = 1,3
    Rsu = RSup(2*i1-1:2*i1, 2*i1-1:2*i1)

    Yuk = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,7
      Call CoupNeutralinoSup(i3, i2, gU1, gSU2, RSu, Yuk, N, coupLC, coupRC)
      c_UNSu_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      c_UNSu_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

   End Do ! i1

  End If

  c_CNW_L = 0._dp
  c_CNW_R = 0._dp

  Do i1 = 1,n_char
   Do i2= 1,n_neut
    Call CoupCharginoNeutralinoW(i1, i2, N, U, V, gSU2               &
                               &, c_CNW_L(i1,i2), c_CNW_R(i1,i2) )
   End Do
  End Do

  c_SmpCN_L = 0._dp
  c_SmpCN_R = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_char
    Do i3 = 1,n_neut
     Call CoupCSCharginoNeutralino(i1, i2, i3, N, U, V, RSpm, Y_l     &
                  &, gU1, gSU2, RP_lam, c_SmpCN_L(i1,i2,i3), c_SmpCN_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  c_NNZ_L = 0.0_dp
  c_NNZ_R = 0.0_dp
  cosW = gSU2 / Sqrt(gU1**2 + gSU2**2)
  Do i1=1,n_neut
   Do i2=1,n_neut
    Call CoupNeutralinoZ(i1, i2, N, gSU2, cosW, c_NNZ_L(i1,i2), c_NNZ_R(i1,i2) )
   End Do
  End Do

  c_NNP0_L = 0.0_dp
  c_NNP0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_P0
     Call CoupNeutralinoPseudoscalar(i1, i2, i3, N, RP0, gU1, gSU2, &
                                  & c_NNP0_L(i1,i2,i3), c_NNP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  c_NNS0_L = 0.0_dp
  c_NNS0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_S0
     Call CoupNeutralinoScalar(i1, i2, i3, N, RS0, gU1, gSU2, &
                             & c_NNS0_L(i1,i2,i3), c_NNS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !--------------------
  ! Loop Contributions
  !--------------------
! LoopContributions(6:7) = .false.

  Call Sigma_Neutralino_2(mN, mN2, c_NNZ_L, c_NNZ_R, mS02, c_NNS0_L       &
        & , c_NNS0_R, mP02, c_NNP0_L, c_NNP0_R, mC, mC2, c_CNW_L, c_CNW_R &
        & , mSpm2, c_SmpCN_L, c_SmpCN_R, mSup2, c_UNSu_L, c_UNSu_R        &
        & , mSdown2, c_DNSd_L, c_DNSd_R, WriteOut, SigS)

  If (WriteOut) Then
   Write(ErrCan,*) "Neut, sigS"
   Do i1=1,7
    If (Maxval(Abs(Aimag(sigs(i1,:)))).Le.Epsilon(1._dp)) Then
     Write(ErrCan,*) (Real(sigs(i1,i2),dp),i2=1,7)
    Else
     Write(ErrCan,*) (sigs(i1,i2),i2=1,7)
    End If
   End Do
  End If
  mat7 = SigS
  Do i1=1,7
   mat7(i1,i1) = mN(i1) + SigS(i1,i1)
  End Do

  If (Maxval(Abs(Aimag(mat7))).Eq.0._dp) Then ! matrix is reel
   Call EigenSystemQP(Real(mat7,dp), mN1L, N7a, kont, test)

   N1L = N7a

   Do i1=1,6
    Do i2=i1+1,7
     If (Abs(mN1L(i1)).Gt.Abs(mN1L(i2))) Then
      work = mN1L(i1)
      mN1L(i1) = mN1L(i2)
      mN1L(i2) = work
      E7 = N1L(i1,:)
      N1L(i1,:) = N1L(i2,:)
      N1L(i2,:) = E7
     End If
    End Do
   End Do
   mN1L2 = mN1L**2

!-------------------------------------------------------------------------
! there is a huge hierarchy between neutrinos and neutralinos
! therefore lets try seesaw approximation if the first attempt fails
!-------------------------------------------------------------------------
   If (kont.Eq.-1006) Then
    kont = 0
    mat4R = Real(mat7(4:7,4:7) ,dp)
    Call EigenSystem(mat4R, mN4, N4R, kont, test)
    Do i1=1,4
     Do i2=i1+1,4
      If (Abs(mN4(i1)).Gt.Abs(mN4(i2))) Then
       work = mN4(i1)
       mN4(i1) = mN4(i2)
       mN4(i2) = work
       E4R = N4R(i1,:)
       N4R(i1,:) = N4R(i2,:)
       N4R(i2,:) = E4R
      End If
     End Do
    End Do

    XiR = Real(mat7(1:3,4:7),dp)
    MinvR = 0._dp
    Do i1=1,4
     Do i2=1,4
      Do i3=1,4
       MinvR(i1,i2) = MinvR(i1,i2) + N4R(i3,i1) * N4R(i3,i2)  / mN4(i3)
      End Do
     End Do
    End Do  
    mat3R = mat7(1:3,1:3) - Matmul( XiR, Matmul(MinvR, Transpose(XiR) ) )
    Call EigenSystem(mat3R, mnu, VnuR, kont, test)
    Do i1=1,3
     Do i2=i1+1,3
      If (Abs(mNu(i1)).Gt.Abs(mNu(i2))) Then
       work = mNu(i1)
       mNu(i1) = mNu(i2)
       mNu(i2) = work
       E4R(1:3) = VNUR(i1,:)
       VNUR(i1,:) = VNUR(i2,:)
       VNUR(i2,:) = E4R(1:3)
      End If
     End Do
    End Do

    R1 = 0._dp
    R1(1:3,1:3) = VnuR
    R1(4:7,4:7) = N4R 
    XiR = Matmul(XiR,MinvR)
    R2(1:3,1:3) = - 0.5_dp * Matmul(XiR, Transpose(XiR))
    R2(4:7,4:7) = - 0.5_dp * Matmul(Transpose(XiR), XiR)
    Do i1=1,7
     R2(i1,i1) = 1._dp + R2(i1,i1)
    End Do
    R2(1:3,4:7) = - XiR
    R2(4:7,1:3) = Transpose(XiR)
    
    N1L = Matmul(R1,R2)
    mn1l(1:3) = mnu
    mn1l(4:7) = mn4
    mN1L2 = mN1L**2
   End If  ! seesaw approximation
 
  Else
   mat72 = Matmul( Transpose(Conjg( mat7)),  mat7)
   Call EigenSystemQP(mat72, mN1L2, N1L, kont, test)
   mat72 = Matmul(Conjg(N1L), Matmul( mat7, Transpose( Conjg( N1L ) ) ) )
   Do i1=1,7
    phaseM =   Sqrt( mat72(i1,i1)   / Abs( mat72(i1,i1) ) )
    N1L(i1,:) = phaseM * N1L(i1,:)
   End Do
   mN1L = Sqrt( mN1L2 )
!-------------------------------------------------------------------------
! there is a huge hierarchy between neutrinos and neutralinos
! therefore lets try seesaw approximation if the first attempt fails
!-------------------------------------------------------------------------
   If (kont.Eq.-1006) Then
    kont = 0
    mat4 = mat7(4:7,4:7) 
    mat42 = Matmul( Transpose(Conjg( mat4)),  mat4)
    Call EigenSystem(mat42, mN4, N4, kont, test)
    mat42 = Matmul(Conjg(N4), Matmul( mat4, Transpose( Conjg( N4 ) ) ) )
    Do i1=1,4
     phaseM =   Sqrt( mat42(i1,i1)   / Abs( mat42(i1,i1) ) )
     N4(i1,:) = phaseM * N4(i1,:)
    End Do
    mN4 = Sqrt( mN4)

    Xi = mat7(1:3,4:7)
    Minv = 0._dp
    Do i1=1,4
     Do i2=1,4
      Do i3=1,4
       Minv(i1,i2) = Minv(i1,i2) + Conjg(N4(i3,i1) * N4(i3,i2) ) / mN4(i3)
      End Do
     End Do
    End Do  
    mat3 = mat7(1:3,1:3) - Matmul( Xi, Matmul(Minv, Transpose(Xi) ) )

    mat32 = Matmul( Transpose(Conjg( mat3)),  mat3)
    Call EigenSystem(mat32, mnu, Vnu, kont, test)
    mat32 = Matmul(Conjg(Vnu), Matmul( mat3, Transpose( Conjg( Vnu) ) ) )
    Do i1=1,3
     phaseM =   Sqrt( mat32(i1,i1)   / Abs( mat32(i1,i1) ) )
     Vnu(i1,:) = phaseM * Vnu(i1,:)
    End Do
    mNu = Sqrt( mNu)

    R1 = 0._dp
    R1(1:3,1:3) = Conjg(Vnu)
    R1(4:7,4:7) = Conjg(N4) 
    Xi = Matmul(Xi,Minv)
    R2(1:3,1:3) = - 0.5_dp * Matmul(Xi, Transpose(Conjg(Xi)))
    R2(4:7,4:7) = - 0.5_dp * Matmul(Transpose(Conjg(Xi)), Xi)
    Do i1=1,7
     R2(i1,i1) = 1._dp + R2(i1,i1)
    End Do
    R2(1:3,4:7) = - Xi
    R2(4:7,1:3) = Transpose(Conjg(Xi))
    
    N1L = Matmul(R1,R2)
    mn1l(1:3) = mnu
    mn1l(4:7) = mn4
    mN1l2 = mN1L**2
   End If  ! seesaw approximation
 
  End If

  N1L = Matmul(N1L, N)

  Do i1=1,7
   If (mN1L(i1).Lt.0._dp) Then
    mN1L(i1) = -mN1L(i1)
    N1L(i1,:) = N1L(i1,:) * (0._dp,1._dp)
   End If
  End Do

  If ((kont.Ne.0).And.(ErrorLevel.Gt.-1)) Then
   Write (ErrCan,*) 'Warning in subroutine NeutralinoMass_Loop_RPtri, ierr =',kont
   Write(Errcan,*) "test",test
   Write (ErrCan,*) 'M1,M2 ',M1,M2
   Write (ErrCan,*) 'gp,g ',gU1, gSU2
   Write (ErrCan,*) 'mu ',mu
   Write (ErrCan,*) 'vevs ',vevSM
   If (ErrorLevel.Eq.2) Call TerminateProgram
  End If

  If (WriteOut) Then
   Write(ErrCan,*) " "
   Write(ErrCan,*) "Test of Diagonalization, resulting matrix is:"
   Write(ErrCan,*) "mN",mN1L
   Call adjungate(N1L,test7)
   test7 = Matmul( Transpose(test7), Matmul(mat7, test7) )
   Do i1=1,7
    Write(ErrCan,*) test7(i1,:)
   End Do
   Write(ErrCan,*) " "
  End If
  Iname = Iname - 1

 End Subroutine NeutralinoMass_Loop_RPtri


 Subroutine NeutrinoMass_1L(MnuL5, gp, g, Y_l, mC2, U, V, mN2, N, mSle2, Rsl &
                         & , mSnu2, Rsn, mNu, Nnu, kont)
 !----------------------------------------------------------------------
 ! calculates neutrino masses + mixing matrix Nnu at the 1-loop level
 ! using the formulas of
 ! M.~Hirsch, E.~Ma, J.~C.~Romao, J.~W.~F.~Valle and A.~Villanova del Moral,
 ! Phys. Rev. D75 (2007) 053006, hep-ph/0606082
 ! input:
 !  MnuL5 .......... dim 5 operator
 !  gp, g .......... U(1) and SU(2) gauge couplings
 !  Y_l ............ lepton Yukawa coupling
 !  mC2 ............ chargino masses squared
 !  U, V ........... chargino mixing matrices
 !  mN2 ............ neutralino masses squared
 !  N .............. neutralino mixing matrix
 !  mSle2 .......... slepton masses squared
 !  Rsl ............ slepton mixing matrix
 !  mSnu2 .......... sneutrino masses squared
 !  Rsn ............ sneutrino mixing matrix
 ! output 
 !  mNu(i) ........... neutrino mass_i
 !  Nnu(i,j) ......... neutrino mixing matrix
 ! written by Werner Porod, 09.01.2006
 ! Note the factor of 1/2 appears because in the RGE running the Higgs field
 ! is still the complex field
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont 
  Complex(dp), Intent(in) :: MnuL5(3,3), Y_l(3,3), U(2,2), V(2,2), N(4,4) &
                          & , Rsl(6,6), Rsn(3,3)
  Real(dp), Intent(in) :: gp, g, mC2(2), mN2(4), mSle2(6), mSnu2(3)
  Real(dp), Intent(out) :: mNu(3)
  Complex(dp), Intent(out) :: Nnu(3,3)

  Integer :: i1,i2,ierr,i3,i4,i5
  Complex(dp) :: mat32(3,3), E3(3), phaseM, mat3(3,3), del(3,3) &
        &  , c_CNuSl(2,3,6), c_N(4), w1
  Real(dp) :: N3a(3,3), eig(3), test(2), c, b, mu2, c_N2(4), N42(4)


  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutrinoMass_1L'

  !------------------------
  ! effective couplings
  !------------------------
  c_N = g * N(:,2) - gp * N(:,1)
  c_N2 = Abs(c_N)**2
  N42 = Abs(N(:,4))**2

  c_CNuSl = 0._dp
   
  Do i2 = 1,3
   Do i3 = 1,6
    w1 = 0._dp
    Do i1 = 1,3
     w1 = w1 -  Conjg( Y_l(i2,i1) * RSl(i3,i1+3) )
    End Do
    Do i1=1,2
     c_CNuSl(i1,i2,i3) = g * U(i1,1) * Conjg(RSl(i3,i2)) + w1 * U(i1,2)
    End Do
   End Do
  End Do
  !------------------------------------
  ! loop contributions
  !------------------------------------
  mu2 = GetRenormalizationScale()

  del = 0._dp

  Do i1=1,6
   Do i2=1,2
    b = B1(0._dp, mC2(i2), mSle2(i1) )
    Do i3=1,3
     Do i4=1,3
      del(i3,i4) = del(i3,i4) + c_CNuSl(i2,i4,i1) * Conjg(c_CNuSl(i2,i3,i1)) * b
     End Do
    End Do
    Do i3=1,2
     c = 4._dp *g *U(i2,1) *Abs(V(i3,2))**2 * C00(mC2(i2), mC2(i3), mSle2(i1) )
     Do i4=1,3
      Do i5=1,3
       del(i4,i5) = del(i4,i5) + Conjg(c_CNuSl(i2,i3,i1)*Rsl(i1,i5) ) * c
      End Do
     End Do
    End Do
   End Do
  End Do

  Do i1=1,3
   Do i2=1,4
    b = 0.5_dp * c_N2(i2) * B1(0._dp, mN2(i2), mSnu2(i1) )
    Do i3=1,3
     Do i4=1,3
      del(i3,i4) = del(i3,i4) + Rsn(i1,i3) * Conjg(Rsn(i1,i4) ) * b
     End Do
    End Do
    Do i3=1,4
     c = 2._dp * c_N2(i2) * N42(i3) * C00(mN2(i2), mN2(i3), mSnu2(i1) )
     Do i4=1,3
      Do i5=1,3
       del(i4,i5) = del(i4,i5) + Rsn(i1,i4) * Conjg(Rsn(i1,i5) ) * c
      End Do
     End Do
    End Do
   End Do
  End Do

  del = oo16pi2 * del 
  mat3 = 0.5_dp * (MnuL5 + Matmul(Transpose(del),MnuL5) + Matmul(MnuL5,del))


  If (Maxval(Abs(Aimag(Mat3))).Lt.  &                           ! matrix is reel
     & (Maxval(Abs(Real(Mat3,dp)))*1.e-3_dp * Epsilon(1._dp))) Then

   Call EigenSystem(Real(Mat3,dp), Eig, N3a, ierr, test)

   Do i1=1,3
    If (Eig(i1).Lt.0._dp) Then
     mNu(i1) = - Eig(i1)
     Nnu(i1,:) = (0._dp,1._dp) * N3a(i1,:)
    Else
     mNu(i1) = Eig(i1)
     Nnu(i1,:) =N3a(i1,:)
    End If
   End Do

   Do i1=1,2
    Do i2=i1+1,3
     If (mNu(i1).Gt.mNu(i2)) Then
      Eig(1) = mNu(i1)
      mNu(i1) = mNu(i2)
      mNu(i2) = Eig(1)
      E3 = Nnu(i1,:)
      Nnu(i1,:) = Nnu(i2,:)
      Nnu(i2,:) = E3
     End If
    End Do
   End Do

  Else

   mat32 = Matmul( Transpose(Conjg( Mat3 ) ), Mat3 )
   Call EigenSystem(mat32, Eig, Nnu, ierr, test)
   mat32 = Matmul(Conjg(Nnu), Matmul( Mat3, Transpose( Conjg( Nnu ) ) ) )
   Do i1=1,3
    phaseM =   Sqrt( mat32(i1,i1)   / Abs( mat32(i1,i1) ) )
    Nnu(i1,:) = phaseM * Nnu(i1,:)
   End Do
   mNu = Sqrt( Eig )

  End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If (ierr.Ne.0) Then
   Write (Errcan,*) 'Warning in subroutine NeutrinoMasses, ierr =',ierr
   Write (Errcan,*) 'MnuL5:',Cmplx(MnuL5(1,:))
   Write (Errcan,*) '      ',Cmplx(MnuL5(2,:))
   Write (Errcan,*) '      ',Cmplx(MnuL5(3,:))
   kont = ierr
   Iname = Iname - 1
   Return
  Endif
  !-------------------------------------------
  ! note, that my mixing matrix is the
  !  Transpose( Conjg( Nnu ) )
  ! compared to standard convention in neutrino physics
  !-------------------------------------------

  Iname = Iname - 1

 Contains 

  Real(dp) Function C00(m02, m12, m22)
  Implicit None
   Real(dp), Intent(in) :: m02, m12, m22
  
   Real(dp), Parameter :: eps=1.e-6_dp, r38 = 3._dp / 8._dp
   Real(dp) :: r1, r2, dr1, dr2

   r1 = m12 / m02 
   r2 = m22 / m02

   dr1 = r1 - 1._dp
   dr2 = r2 - 1._dp

   If ((dr1.Eq.0._dp).And.(dr2.Eq.0._dp)) Then
    C00 = - 0.25_dp * Log(m02/mu2)

   Else If ((dr1.Eq.0._dp).And.(dr2.Le.eps)) Then
    C00 = dr2 * (1._dp + 0.25_dp * dr2)/ 12._dp - 0.25_dp * Log(m02/mu2)

   Else If ((dr1.Le.eps).And.(dr2.Eq.0._dp)) Then
    C00 = dr1 * (1._dp + 0.25_dp * dr1)/ 12._dp - 0.25_dp * Log(m02/mu2)

   Else If (dr1.Eq.0._dp) Then
    C00 = 0.25_dp * ( dr2 - r2**2*Log(r2)) / dr2**2 &
      & + r38 - 0.25_dp * Log(m02/mu2)

   Else If (dr2.Eq.0._dp) Then
    C00 = 0.25_dp * ( dr1 - r1**2*Log(r1)) / dr1**2 &
      & + r38 - 0.25_dp * Log(m02/mu2)

   Else If ((dr1.Le.eps).And.(dr2.Le.eps)) Then
    C00 = ( dr1 * (1._dp + 0.25_dp * dr1)                     &
      &   + dr2 * (1._dp + 0.25_dp * (dr2 + dr1) ) ) / 12._dp &
      & - 0.25_dp * Log(m02/mu2)

   Else If (r1.Eq.r2) Then
    C00 = - r2 * ( dr2 + ( dr2 - 1._dp ) * Log(r2) ) / (4._dp * dr2**2) &
      & + r38 - 0.25_dp * Log(m02/mu2)
   Else
    C00 = (-r1**2 *dr2 * Log(r1) + dr1 * r2**2 * Log(r2) ) &
      &   / (4._dp * dr1 * dr2 * (r1 - r2) )               &
      & + r38 - 0.25_dp * Log(m02/mu2)

   End If

  End Function C00

 End Subroutine NeutrinoMass_1L



#ifdef GENERATIONMIXING
 Subroutine One_Loop_Tadpoles_MSSM(gU1, gSU2, vevs_DR, Y_l, Y_d, Y_u       &
         & , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R, mu, Ae, Ad, Au            &
         & , mSneutrino2, Rsneut, mSlepton2, Rslepton, mSdown2, Rsdown     &
         & , mSup2, RSup, mSpm2, RSpm, mC, mC2, U, V, mP02, RP0, mS02, RS0 &
         & , mN, mN2, N, mZ2, mW2, tadpole)
#else
 Subroutine One_Loop_Tadpoles_MSSM(gU1, gSU2, vevs_DR, Y_l, Y_d, Y_u, mu    &
         & , Ae, Ad, Au, mSneutrino2, mSlepton2, Rslepton, mSdown2, Rsdown  &
         & , mSup2, RSup, mSpm2, RSpm, mC, mC2, U, V, mP02, RP0, mS02, RS0  &
         & , mN, mN2, N, mZ2, mW2, tadpole)
#endif
 !----------------------------------------------------------------
 ! In this subroutine the 1-loop tadpole contributions are calculated
 ! The formulae of J. Bagger et al, Nucl. Phys. B491 are used.
 ! input:
 !  - gU1 ........ the U(1) gauge coupling
 !  - gSU2 ....... the SU(2) gauge coupling 
 !  - vevs_DR .... MSSM vevs
 !  - Y_l ........ lepton Yukawa couplings
 !  - Y_d ........ d-quark Yukawa couplings
 !  - Y_u ........ u-quark Yukawa couplings
 !  - mu ......... higgsino parameter
 !  - Ae ......... slepton trilinear couplings
 !  - Ad ......... d-squark trilinear couplings
 !  - Au ......... u-squark trilinear couplings
 !  - mSneutrino2  masses of the sneutrinos squared
 !  - RSneut ..... mixing matrix of the sneutrinos
 !  - mSlepton2 .. masses of the sleptons squared
 !  - RSlepton ... mixing matrix of the sleptons  
 !  - mSdown2 .... masses of the s-downs squared
 !  - RSdown ..... mixing matrix of the s-downs
 !  - mSup2 ...... masses of the s-ups squared
 !  - RSup ....... mixing matrix of the s-ups
 !  - mSpm2 ...... masses of the charged Higgs squared
 !  - RSpm ....... mixing matrix of the charged scalar Higgs
 !  - mC ......... masses of the charginos
 !  - mC2 ........ masses of the charginos squared
 !  - U, V ....... mixing matrices of charginos
 !  - mP02 ....... masses of the pseudoscalar Higgs squared
 !  - RP0 ........ mixing matrix of the pseudoscalar Higgs
 !  - mS02 ....... masses of the scalar Higgs squared
 !  - RS0 ........ mixing matrix of the scalar Higgs
 !  - mN ......... masses of the neutralinos
 !  - mN2 ........ masses of the neutralinos squared
 !  - N .......... mixing matrix of neutralinos
 ! output:
 !  - tadpole ..... the MSSM 1-loop tadpole contribution divided by the vevs
 ! written by Werner Porod
 ! 09.10.01: portation to f90
 ! 10.10.01: including flavour mixing in quark sector
 ! 12.12.03: - taking mW and mZ as input instead of assuming that
 !             they are taking from the module StandardModel
 !           - scalar Higgs contribution has been a factor 2 too large 
 !----------------------------------------------------------------
 Implicit None
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3), Rsneut(3,3)    &
     & , Rslepton(6,6), Rsdown(6,6), RSup(6,6), mu, Ae(3,3), Ad(3,3)      &
     & , Au(3,3), RSpm(2,2), U(2,2), V(2,2), N(4,4), uU_L(3,3), uU_R(3,3) &
     & , uD_L(3,3), uD_R(3,3), uL_L(3,3), uL_R(3,3)
#else
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3), Rslepton(6,6)    &
     & , Rsdown(6,6), RSup(6,6), mu, Ae(3,3), Ad(3,3), Au(3,3), RSpm(2,2)   &
     & , U(2,2), V(2,2), N(4,4)
#endif
  Real(dp), Intent(in) :: gU1, gSU2, mSneutrino2(3), mSlepton2(6)           &
     & , mSdown2(6), mSup2(6), vevs_DR(2), mSpm2(2), mC(2), mC2(2), mP02(2) &
     & , RP0(2,2), mS02(2), RS0(2,2), mN(4), mN2(4), mZ2, mW2
  Real(dp), Intent(out) :: tadpole(2)

  Integer :: i1, i2
  Real(dp) :: id2R(2,2), coup, e_d, e_u
  Complex(dp) :: sumI(2), A0m, coupC, coupLC, coupRC, Rsf2(2,2)
#ifdef GENERATIONMIXING
  Complex(dp) :: mat3(3,3)
#endif
  Logical :: WriteOut

  Iname = Iname + 1
  NameOfUnit(Iname) = "One_Loop_Tadpoles_MSSM"

  !----------------------------------
  ! In case I need the contributions
  !----------------------------------
  If ((WriteOneLoopContributions.Eq.2).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions in One_Loop_Tadpoles_MSSM:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If

  sumI = 0._dp
  tadpole = 0._dp  
  id2R = 0._dp
  id2R(1,1) = 1._dp
  id2R(2,2) = 1._dp
  !-----------------
  ! SM contribtions
  !-----------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,3
    Call CoupFermionScalar(i1, i1, 1, -0.5_dp, Y_d, uD_L, uD_R, id2R &
                         &, coupC, coupRC)
    sumI(1) = 3._dp * mf_d(i1) * Real( CoupC,dp ) * A0(mf_d2(i1))
    Call CoupFermionScalar(i1, i1, 1, -0.5_dp, Y_l, uL_L, uL_R, id2R &
                         &, coupC, coupRC)
    sumI(1) = sumI(1) +  mf_l(i1) * Real( coupC,dp ) * A0(mf_l2(i1))
    Call CoupFermionScalar(i1, i1, 2, 0.5_dp, Y_u, uU_L, uU_R, id2R &
                         &, coupC, coupRC)
    sumI(2) = 3._dp * mf_u(i1) * Real( coupC,dp ) * A0(mf_u2(i1))
    sumI = 4.d0 * sumI
    If (WriteOut) Write(ErrCan,*) "fermions",i1,sumI
    tadpole = tadpole + sumI 
   End Do
  Else ! .not. GenerationMixing
#endif
   Do i1=1,3
    sumI(1) = - 3._dp * mf_d(i1) * Y_d(i1,i1)  * A0(mf_d2(i1))
    sumI(1) = sumI(1) - mf_l(i1) * Y_l(i1,i1) * A0(mf_l2(i1))
    sumI(2) = - 3._dp * mf_u(i1) * Y_u(i1,i1) * A0(mf_u2(i1))
    If (WriteOut) Write(ErrCan,*) "fermions",i1,sumI
    sumI = 2._dp * sqrt2 * sumI
    tadpole = tadpole + sumI
   End Do
#ifdef GENERATIONMIXING
  End If
#endif
  !-------------
  ! sneutrinos
  !-------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   mat3 = zeroC
   Do i1=1,3
    A0m = A0( mSneutrino2(i1) )
    Call CoupScalarSfermion3(1, i1, i1, id2R, 0.5_dp, 0._dp, mat3, Rsneut &
                           &, mat3, mu, vevs_DR, gU1, gSU2, coupC )
    sumI(1) =  - coupC * A0m
    Call CoupScalarSfermion3(2, i1, i1, id2R, 0.5_dp, 0._dp, mat3, Rsneut &
                           &, mat3, mu, vevs_DR, gU1, gSU2, coupC )
    sumI(2) =  - coupC * A0m
    If (WriteOut) Write(ErrCan,*) "sneutrinos",i1,sumI
    tadpole = tadpole + sumI
   End Do
  Else ! .not. GenerationMixing
#endif
   Call CoupScalarSfermion3(1, 1, 1, id2R, 0.5_dp, 0._dp, ZeroC, id2C &
                          &, ZeroC, mu, vevs_DR, gU1, gSU2, coupC)
   Call CoupScalarSfermion3(2, 1, 1, id2R, 0.5_dp, 0._dp, ZeroC, id2C &
                          &, ZeroC, mu, vevs_DR, gU1, gSU2, coupRC)
   
   Do i1=1,3
    A0m = A0( mSneutrino2(i1) )
    sumI(1) =  - coupC * A0m
    sumI(2) =  - coupRC * A0m
    If (WriteOut) Write(ErrCan,*) "sneutrinos",i1,sumI
    tadpole = tadpole + sumI
   End Do
#ifdef GENERATIONMIXING
  End If
#endif
  !---------------------
  ! remaining sfermions
  !--------------------
  e_d = -1._dp / 3._dp
  e_u = 2._dp / 3._dp
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,6
    A0m = A0( mSlepton2(i1) )
    Call CoupScalarSfermion3(1, i1, i1, id2R, -0.5_dp, -1._dp, Y_l, Rslepton  &
                            &, Ae, mu, vevs_DR, gU1, gSU2, coupC )
    sumI(1) =  - coupC * A0m
    Call CoupScalarSfermion3(2, i1, i1, id2R, -0.5_dp, -1._dp, Y_l, Rslepton  &
                            &, Ae, mu, vevs_DR, gU1, gSU2, coupC )
    sumI(2) =  - coupC * A0m
    If (WriteOut) Write(ErrCan,*) "sleptons",i1,sumI
    tadpole = tadpole + sumI
   End Do
   Do i1=1,6
    A0m = A0( mSup2(i1) )
    Call CoupScalarSfermion3(1, i1, i1, id2R, 0.5_dp, e_u, Y_u, Rsup   &
                            &, Au, mu, vevs_DR, gU1, gSU2, coupC )
    sumI(1) =  - 3._dp * coupC * A0m
    Call CoupScalarSfermion3(2, i1, i1, id2R, 0.5_dp, e_u, Y_u, Rsup   &
                            &, Au, mu, vevs_DR, gU1, gSU2, coupC )
    sumI(2) =  - 3._dp * coupC * A0m
    If (WriteOut) Write(ErrCan,*) "sups",i1,sumI
    tadpole = tadpole + sumI
   End Do
   Do i1=1,6
    A0m = A0( mSdown2(i1) )
    Call CoupScalarSfermion3(1, i1, i1, id2R, -0.5_dp, e_d, Y_d, Rsdown   &
                            &, Ad, mu, vevs_DR, gU1, gSU2, coupC )
    sumI(1) =  - 3._dp * coupC * A0m
    Call CoupScalarSfermion3(2, i1, i1, id2R, -0.5_dp, e_d, Y_d, Rsdown   &
                            &, Ad, mu, vevs_DR, gU1, gSU2, coupC )
    sumI(2) =  - 3._dp * coupC * A0m
    If (WriteOut) Write(ErrCan,*) "sdowns",i1,sumI
    tadpole = tadpole + sumI
   End Do
  Else ! .not. GenerationMixing
#endif
   Do i1=1,3
    Rsf2 = RSlepton(2*i1-1:2*i1,2*i1-1:2*i1)
    Do i2 =1,2
     A0m = A0( mSlepton2(2*i1-2+i2) )
     Call CoupScalarSfermion3(1, i2, i2, id2R, -0.5_dp, -1._dp, Y_l(i1,i1) &
                            &, Rsf2, Ae(i1,i1), mu, vevs_DR, gU1, gSU2, coupC )
     sumI(1) =  - coupC * A0m
     Call CoupScalarSfermion3(2, i2, i2, id2R, -0.5_dp, -1._dp, Y_l(i1,i1) &
                            &, Rsf2, Ae(i1,i1), mu, vevs_DR, gU1, gSU2, coupC )
     sumI(2) =  - coupC * A0m
     If (WriteOut) Write(ErrCan,*) "sleptons",i1,i2,sumI
     tadpole = tadpole + sumI
    End Do
   End Do
   Do i1=1,3
    Rsf2 = RSup(2*i1-1:2*i1,2*i1-1:2*i1)
    Do i2 =1,2
     A0m = A0( mSup2(2*i1-2+i2) )
     Call CoupScalarSfermion3(1, i2, i2, id2R, 0.5_dp, e_u, Y_u(i1,i1) &
                            &, Rsf2, Au(i1,i1), mu, vevs_DR, gU1, gSU2, coupC )
     sumI(1) =  - 3._dp * coupC * A0m
     Call CoupScalarSfermion3(2, i2, i2, id2R, 0.5_dp, e_u, Y_u(i1,i1) &
                            &, Rsf2, Au(i1,i1), mu, vevs_DR, gU1, gSU2, coupC )
     sumI(2) =  - 3._dp * coupC * A0m
     If (WriteOut) Write(ErrCan,*) "sups",i1,i2,sumI
     tadpole = tadpole + sumI
    End Do
   End Do
   Do i1=1,3
    Rsf2 = RSdown(2*i1-1:2*i1,2*i1-1:2*i1)
    Do i2 =1,2
     A0m = A0( mSdown2(2*i1-2+i2) )
     Call CoupScalarSfermion3(1, i2, i2, id2R, -0.5_dp, e_d, Y_d(i1,i1) &
                            &, Rsf2, Ad(i1,i1), mu, vevs_DR, gU1, gSU2, coupC )
     sumI(1) =  - 3._dp * coupC * A0m
     Call CoupScalarSfermion3(2, i2, i2, id2R, -0.5_dp, e_d, Y_d(i1,i1) &
                            &, Rsf2, Ad(i1,i1), mu, vevs_DR, gU1, gSU2, coupC )
     sumI(2) =  - 3._dp * coupC * A0m
     If (WriteOut) Write(ErrCan,*) "sdowns",i1,i2,sumI
     tadpole = tadpole + sumI
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  !----------------------
  ! gauge bosons
  !----------------------
  sumI = 1.5_dp *  gSU2**2 * A0(mW2) + 0.75_dp * (gSU2**2 + gU1**2) * A0(mZ2)
  sumI = sumI * vevs_DR 
  If (WriteOut) Write(ErrCan,*) "gauge bosons",sumI
  tadpole = tadpole + sumI
  !---------------------------------------------------------
  ! charged Higgs bosons
  ! I have negative sign in the definition of the couplings
  !---------------------------------------------------------
  Do i1 = 1,n_Spm
   A0m = A0(mSpm2(i1))
   Call CoupChargedScalarScalar3(i1, i1, 1, RSpm, id2R, vevs_DR, gU1, gSU2  &
                               &, coupC)
   sumI(1) = - coupC * A0m
   Call CoupChargedScalarScalar3(i1, i1, 2, RSpm, id2R, vevs_DR, gU1, gSU2  &
                               &, coupC)
   sumI(2) = - coupC * A0m
   If (WriteOut) Write(ErrCan,*) "S+ ",i1,sumI
   tadpole = tadpole + sumI
  End Do
  !-----------------------
  ! Charginos
  !-----------------------
  Do i1=1,n_char
   A0m = 4._dp *  mC(i1) * A0(mC2(i1))
   Call CoupCharginoScalar(i1, i1, 1, U, V, id2R, gSU2, coupLC, coupRC)
   sumI(1) = coupLC * A0m
   Call CoupCharginoScalar(i1, i1, 2, U, V, id2R, gSU2, coupLC, coupRC)
   sumI(2) = coupLC * A0m
   If (WriteOut) Write(ErrCan,*) "Charginos ",i1,sumI
   tadpole = tadpole + sumI
  End Do
 !-----------------------------------------------------
 ! pseudoscalar Higgs bosons, including Goldstone boson
 ! I have negative sign in the definition of the couplings
 !-----------------------------------------------------
  Do i1 = 1,n_P0
   A0m = A0(mP02(i1))
   Call CoupPseudoScalarScalar3(i1, i1, 1, RP0, id2R, gU1, gSU2, vevs_DR &
                              &, coup )
   sumI(1) = - coup * A0m
   Call CoupPseudoScalarScalar3(i1, i1, 2, RP0, id2R, gU1, gSU2, vevs_DR &
                              &, coup )
   sumI(2) = - coup * A0m
   If (WriteOut) Write(ErrCan,*) "P0 ",i1,sumI
   tadpole = tadpole + sumI
  End Do
 !-----------------------------------------------------
 ! scalar Higgs bosons
 ! I have negative sign in the definition of the couplings
 !-----------------------------------------------------
  Do i1 = 1,n_S0
   A0m = 0.5_dp * A0(mS02(i1))
   Call CoupScalar3a(1, i1, i1, RS0, gU1, gSU2, vevs_DR, coup )
   sumI(1) = - coup * A0m
   Call CoupScalar3a(2, i1, i1, RS0, gU1, gSU2, vevs_DR, coup )
   sumI(2) = - coup * A0m
   If (WriteOut) Write(ErrCan,*) "S0 ",i1,sumI
   tadpole = tadpole + sumI
  End Do
 !-----------------------
 ! Neutralinos
 !-----------------------
  Do i1=1,n_neut
   A0m = 2._dp *  mN(i1) * A0(mN2(i1))
   Call CoupNeutralinoScalar(i1, i1, 1, N, id2R, gU1, gSU2, coupLC, coupRC )
   sumI(1) = coupLC * A0m
   Call CoupNeutralinoScalar(i1, i1, 2, N, id2R, gU1, gSU2, coupLC, coupRC )
   sumI(2) = coupLC * A0m
   If (WriteOut) Write(ErrCan,*) "Neutralinos ",i1,sumI
   tadpole = tadpole + sumI
  End Do
  If (WriteOut) Write(ErrCan,*) "sum",tadpole
  tadpole =  oo16pi2 * tadpole / vevs_DR

  If (WriteOut) Write(ErrCan,*) "tadpoles",tadpole

  Iname = Iname - 1

 End Subroutine One_Loop_Tadpoles_MSSM


#ifdef GENERATIONMIXING
 Subroutine PiChargedScalar(p2, gU1, gSU2, Y_d, Y_u, Y_l, vevs_DR, mSpm2,RSpm &
   & , A_d, A_u, A_l, mu, uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                  &
   & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino      &
   & , mSlepton2, RSlepton, mP02, RP0, mS02, RS0, mC, mC2, U, V, mN, mN2, N   &
   & , mZ2, mW2, res)
#else
 Subroutine PiChargedScalar(p2, gU1, gSU2, Y_d, Y_u, Y_l, vevs_DR, mSpm2     &
   & ,RSpm , A_d, A_u, A_l, mu, mUSquark2, RUSquark, mDSquark2, RDSquark     &
   & , mSneutrino2, mSlepton2, RSlepton, mP02, RP0, mS02, RS0, mC, mC2, U, V &
   & , mN, mN2, N, mZ2, mW2, res)
#endif
 !-------------------------------------------------------------------------
 ! calculates the 1-loop self energies of the pseudoscalar higgs
 ! The formulae of J. Bagger et al, Nucl. Phys. B491 are used.
 ! Note, that gSU3 has to be set to 0 for sleptons
 ! the input is:
 !  - p2 ......... the outer momentum squared 
 !  - gU1 ........ the U(1) gauge coupling
 !  - gSU2 ....... the SU(2) gauge coupling
 !  - Y_d ........ d-quark Yukawa couplings
 !  - Y_u ........ u-quark Yukawa couplings
 !  - Y_l ........ lepton Yukawa couplings
 !  - vevs_DR .... MSSM vevs
 !  - mSpm2 ...... masses of the charged Higgs squared
 !  - RSpm ....... mixing matrix of the charged scalar Higgs
 !  - A_d ........ d-squark trilinear couplings
 !  - A_u ........ u-squark trilinear couplings
 !  - A_l ........ slepton trilinear couplings
 !  - mu ......... higgsino parameter
 ! 
 ! the output is:
 !  res  
 ! written by Werner Porod, 30.11.99
 ! 10.10.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp) :: p2, vevs_DR(2), gU1, gSU2, mUSquark2(6), mDSquark2(6)         &
    & , mSneutrino2(3), mSlepton2(6), mSpm2(:), mP02(:), RP0(:,:), mS02(:)  &
    & , RS0(:,:), mC(:), mC2(:), mN(:), mN2(:), mZ2, mW2
  Complex(dp), Intent(in) :: Y_d(3,3), Y_u(3,3), Y_l(3,3), RUSquark(6,6)    &
    & , RDSquark(6,6), RSlepton(6,6), RSpm(:,:), A_d(3,3)  &
    & , A_u(3,3), A_l(3,3), mu, U(:,:), V(:,:), N(:,:)
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: RSneutrino(3,3), uU_L(3,3), uU_R(3,3)  &
    & , uD_L(3,3), uD_R(3,3), uL_L(3,3), uL_R(3,3)
#endif
  Complex(dp), Intent(out) :: res

  Integer :: i1, i2, i3, i_gen
  Real(dp) :: sinW2_DR, cosW2_DR, e_u, e_d, B0m2, G0m2
  Complex(dp) :: sumI, coupC, coupCL, coupCR, A0m2, Rsf(2,2), Rsfp(2,2)     &
    & , yuku, yukd, Ad, Au, Zero3(3,3)
  Logical :: WriteOut

  Iname = Iname + 1
  NameOfUnit(Iname) = "PiChargedScalar"

  !----------------------------------
  ! In case I need the contributions
  !----------------------------------
  If ((WriteOneLoopContributions.Eq.4).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions in PiChargedScalar:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If
 !-----------------------
 ! intialisation
 !-----------------------
  res = (0._dp, 0._dp)

  e_u = 2._dp / 3._dp
  e_d = -1._dp / 3._dp

  sinW2_DR = gU1**2 / (gSU2**2 + gU1**2)
  cosW2_DR = 1._dp - sinW2_DR

  Zero3 = 0._dp

#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   !---------------
   ! quarks
   !---------------
   Do i1=1,3
    Do i2=1,3
     G0m2 = Real( Gloop(p2, mf_u2(i1), mf_d2(i2)),dp )
     B0m2 = 4._dp * mf_u(i1) * mf_d(i1) * Real( B0(p2,mf_u2(i1), mf_d2(i1)),dp)
     Call CoupChargedScalarFermion(2, i2, i1, RSpm, Y_d, uD_L, uD_R, Y_u &
                                  &, uU_L, uU_R, coupCL, coupCR)
     sumI = 3._dp * ( (Abs( coupCL )**2 + Abs( coupCR )**2) * G0m2  &
         &         - Conjg( coupCL ) * coupCR * B0m2 )
     If (WriteOut) Write(ErrCan,*) "quarks ",i1,i2,sumI
     res = res + sumI
    End Do
   End Do
   !-------------
   ! squarks
   !-------------
   Do i1=1,6
    Call CoupChargedScalarSfermion4(2, 2, i1, i1, RSpm, 0.5_dp, e_u  &
                                  & , Y_u, Y_d, RUSquark, gU1, gSU2, coupC)
    SumI = -3._dp * coupC *  A0( mUSquark2(i1) )
    If (WriteOut) Write(ErrCan,*) "u-squarks ",i1,sumI
    res = res + SumI

    Call CoupChargedScalarSfermion4(2, 2, i1, i1, RSpm, -0.5_dp, e_d  &
                                  & , Y_d, Y_u, RDSquark, gU1, gSU2, coupC)
    SumI = -3._dp * coupC * A0( mDSquark2(i1) )
    If (WriteOut) Write(ErrCan,*) "d-squarks ",i1,sumI
    res = res + SumI

    Do i2=1,6
     Call CoupChargedScalarSfermion3(2, i2, i1, RSpm, gSU2, vevs_DR, mu  &
                            &, Y_d, A_d, RDSquark, Y_u, A_u, RUSquark, coupC)
     SumI = 3._dp *  Abs(coupC)**2*Real(B0(p2,mUSquark2(i1), mDSquark2(i1)),dp)
     If (WriteOut) Write(ErrCan,*) "d-, u-squarks ",i1,i2,i3,sumI
     res = res + SumI
    End Do 
   End Do 
   !---------------
   ! leptons
   !---------------
   Do i1=1,3
    Do i2=1,3
     Call CoupChargedScalarFermion(2, i2, i1, RSpm, Y_l, uL_L, uL_R, Zero3 &
                                  &, id3C, id3C, coupCL, coupCR)
     sumI = Abs( coupCL )**2 * Real( Gloop(p2, 0._dp, mf_l2(i2)),dp )
     If (WriteOut) Write(ErrCan,*) "leptons ",i1,i2,sumI
     res = res + sumI
    End Do
   End Do
   !-------------
   ! Sleptons
   !-------------
   Do i1=1,3
    Call CoupChargedScalarSfermion4(2, 2, 1, 1, RSpm, 0.5_dp, 0._dp      &
                                  &, Zero3, Y_l, RSneutrino, gU1, gSU2, coupC)
    SumI = - coupC * A0(  mSneutrino2(i1) )
    If (WriteOut) Write(ErrCan,*) "sneutrinos ",i1,sumI
    res = res + SumI

    Call CoupChargedScalarSfermion4(2, 2, i1, i1, RSpm, -0.5_dp, -1._dp  &
                                  & , Y_l, Zero3, RSlepton, gU1, gSU2, coupC)
    SumI = -coupC * A0( mSlepton2(i1) )
    Call CoupChargedScalarSfermion4(2, 2, i1+3, i1+3, RSpm, -0.5_dp, -1._dp  &
                                  & , Y_l, Zero3, RSlepton, gU1, gSU2, coupC)
    SumI = sumI - coupC * A0( mSlepton2(i1+3) )
    If (WriteOut) Write(ErrCan,*) "sleptons ",i1,sumI
    res = res + SumI

    Do i2=1,6
     B0m2 = Real( B0(p2, mSneutrino2(i1), mSlepton2( i2)),dp )
     Call CoupChargedScalarSfermion3(2, i2, i1, RSpm, gSU2, vevs_DR , mu &
                       &, Y_l, A_l, RSlepton, Zero3, Zero3, RSneutrino, coupC)
     SumI = Abs( coupC )**2 * B0m2
     If (WriteOut) Write(ErrCan,*) "sneutrinos, sleptons ",i1,i2,sumI
     res = res + SumI
    End Do 
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,3
   !---------------
   ! quarks
   !---------------
    yuku = Y_u(i1,i1)
    yukd = Y_d(i1,i1)
    G0m2 = Real( Gloop(p2, mf_u2(i1), mf_d2(i1)),dp )
    B0m2 = 4._dp * mf_u(i1) * mf_d(i1) * Real( B0(p2,mf_u2(i1), mf_d2(i1)),dp )
    Call CoupChargedScalarFermion(2,RSpm,yukD,yukU,coupCL,coupCR)
    sumI = 3._dp * ( (Abs( coupCL )**2 + Abs( coupCR )**2) * G0m2  &
         &         - Conjg( coupCL ) * coupCR * B0m2 )
    If (WriteOut) Write(ErrCan,*) "quarks ",i1,sumI
    res = res + sumI
   End Do
   !-------------
   ! squarks
   !-------------
   Do i1=1,3
    i_gen = 2 * i1 - 1
    Rsf = RUSquark( i_gen:i_gen+1, i_gen:i_gen+1)
    Rsfp = RDSquark( i_gen:i_gen+1, i_gen:i_gen+1)
    Au = A_u(i1,i1)
    Ad = A_d(i1,i1)
    i_gen = i_gen - 1
    yuku = Y_u(i1,i1)
    yukd = Y_d(i1,i1)
    Do i2=1,2
     A0m2 = - A0( mUSquark2(i_gen + i2) ) 
     Call CoupChargedScalarSfermion4(2, 2, i2, i2, RSpm, 0.5_dp, e_u  &
                                   & , yuku, yukd, Rsf, gU1, gSU2, coupC)
     SumI =  3._dp * coupC * A0m2
     If (WriteOut) Write(ErrCan,*) "u-squarks ",i1,i2,sumI
     res = res + SumI

     A0m2 = - A0( mDSquark2(i_gen + i2) ) 
     Call CoupChargedScalarSfermion4(2, 2, i2, i2, RSpm, -0.5_dp, e_d  &
                                   & , yukd, yuku, Rsfp, gU1, gSU2, coupC)
     SumI =  3._dp * coupC * A0m2
     If (WriteOut) Write(ErrCan,*) "d-squarks ",i1,i2,sumI
     res = res + SumI

     Do i3=1,2
      B0m2 = Real( B0(p2,mUSquark2(i_gen + i2), mDSquark2(i_gen + i3) ),dp )
      Call CoupChargedScalarSfermion3(2, i3, i2, RSpm, gSU2, vevs_DR  &
                                 &, mu, yukd, Ad, Rsfp, yuku, Au, Rsf, coupC)
      SumI = 3._dp *  Abs( coupC )**2 * B0m2
      If (WriteOut) Write(ErrCan,*) "d-, u-squarks ",i1,i2,i3,sumI,coupc
      res = res + SumI
     End Do
    End Do 
   End Do 
   !---------------
   ! leptons
   !---------------
   yuku = (0._dp,0._dp)
   Do i1=1,3
    yukd = Y_l(i1,i1)
    G0m2 = Real( Gloop(p2, 0._dp,mf_L2(i1)),dp )
    Call CoupChargedScalarFermion(2, RSpm, yukD, yukU, coupCL, coupCR)
    sumI = Abs( coupCL )**2 * G0m2 
    If (WriteOut) Write(ErrCan,*) "leptons ",i1,sumI
    res = res + sumI
   End Do
   !-------------
   ! Sleptons
   !-------------
   Rsfp = id2C
   Au = 0._dp
   Do i1=1,3
    yukd = Y_l(i1,i1)
    i_gen = 2 * i1 - 1
    Rsf = RSlepton( i_gen:i_gen+1, i_gen:i_gen+1)
    Ad = A_l(i1,i1)
    
    A0m2 = - A0(  mSneutrino2(i1) )
    Call CoupChargedScalarSfermion4(2, 2, 1, 1, RSpm, 0.5_dp, 0._dp  &
                                  &, yuku, yukd, Rsfp, gU1, gSU2, coupC)
    SumI = coupC * A0m2
    If (WriteOut) Write(ErrCan,*) "sneutrinos ",i1,sumI
    res = res + SumI

    i_gen = i_gen - 1
    Do i2=1,2
     A0m2 = - A0(mSlepton2(i_gen + i2))
     Call CoupChargedScalarSfermion4(2, 2, i2, i2, RSpm, -0.5_dp,-1._dp  &
                                    & ,yukd, yuku, Rsf, gU1, gSU2, coupC)
     SumI = coupC * A0m2
     If (WriteOut) Write(ErrCan,*) "sleptons ",i1,i2,sumI
     res = res + SumI
  
     B0m2 = Real( B0(p2, mSneutrino2(i1), mSlepton2(i_gen + i2)),dp )
     Call CoupChargedScalarSfermion3(2, i2, 1, RSpm, gSU2, vevs_DR  &
                                 &, mu, yukd, Ad, Rsf, yuku, Au, Rsfp, coupC)
     SumI = Abs( coupC )**2 * B0m2
     If (WriteOut) Write(ErrCan,*) "sneutrinos, sleptons ",i1,i2,sumI
     res = res + SumI
    End Do 
   End Do
#ifdef GENERATIONMIXING
  End If
#endif
 
 !-------------------
 ! gauge bosons
 !-------------------
  sumI =  Real( Floop(p2, mSpm2(2), 0._dp),dp ) * gSU2**2 * sinW2_DR
  If (WriteOut) Write(ErrCan,*) "photon ",sumI
  res = res + sumI

  sumI = gSU2**2 * (2._dp*A0(mW2) + (cosW2_DR-sinW2_DR)**2 * A0(mZ2)/cosW2_DR )
  If (WriteOut) Write(ErrCan,*) "W, Z ",sumI
  res = res + sumI
 !--------------------
 ! charged Higgs
 !--------------------
  sumI = (gSU2**2 + gU1**2) * ( 2._dp * (Rspm(1,1)*RSpm(1,2))**2 - 0.25_dp ) &
       &                    * A0(mSpm2(1) )
  If (WriteOut) Write(ErrCan,*) "G+ ",sumI
  res = res + sumI

  sumI = 0.5_dp * (gSU2**2 + gU1**2)  &
   &            * (Rspm(1,1)**2 - RSpm(1,2)**2 )**2 * A0(mSpm2(2) )
  If (WriteOut) Write(ErrCan,*) "H+ ",sumI
  res = res + sumI

  Call CoupChargedScalarZ(2, 2, gSU2, sinW2_DR, RSpm, coupC)
  sumI = Abs( coupC )**2 * Real( Floop(p2,mSpm2(2) ,mZ2),dp )
  If (WriteOut) Write(ErrCan,*) "Z H+ ",sumI
  res = res + sumI

 !-------------------------
 ! pseudoscalar Higgs
 !-------------------------
  sumI = 0.25_dp * gSU2**2 * ( Real( Floop(p2,mP02(2),mW2),dp )   &
                             & + mW2 *Real( B0(p2,mP02(2),mW2),dp ) )
  If (WriteOut) Write(ErrCan,*) "A0 W+ ",sumI
  res = res + sumI

  Do i1=1,n_P0
   Call CoupChargedScalarPseudoScalar4(2, 2, i1, i1, RSpm, RP0, gU1, gSU2 &
                                     &, coupC)
   sumI = 0.5_dp * coupC * A0( mP02(i1) )
   If (WriteOut) Write(ErrCan,*) "P0 ",i1,sumI
   res = res + sumI
  End Do      
  !-------------------------
  ! neutral Higgs
  !-------------------------
  Do i1=1,n_S0
   Call CoupChargedScalarScalarW(2, i1, gSU2, RSpm, RS0, coupC)
   sumI = Abs( coupC )**2 * Real( Floop(p2, mS02(i1), mW2),dp )
   If (WriteOut) Write(ErrCan,*) "S0 W+ ",i1,sumI,coupc
   res = res + sumI

   Call CoupChargedScalarScalar4(2, 2, i1, i1, RSpm, RS0, gU1, gSU2, coupC)
   sumI = 0.5_dp * coupC * A0( mS02(i1) )
   If (WriteOut) Write(ErrCan,*) "S0 ",i1,sumI
   res = res + sumI

   Do i2=1,n_Spm
    Call CoupChargedScalarScalar3(2, i2, i1, RSpm, RS0, vevs_DR, gU1, gSU2 &
                                 &, coupC)
    sumI =  Abs( coupC )**2 * Real( B0(p2, mSpm2(i2), mS02(i1)),dp )
    If (WriteOut) Write(ErrCan,*) "S0 S+ ",i1,i2,sumI
    res = res + sumI
   End Do
  End Do
  !-------------------------
  ! charginos + neutralinos
  !-------------------------
  Do i1=1,n_char
   Do i2=1,n_neut
    Call CoupCSCharginoNeutralino(2, i1, i2, N, U, V, RSpm, gU1, gSU2  &
                                &, coupCL, coupCR)
    sumI = (Abs(coupCL)**2+Abs(coupCR)**2)*Real(Gloop(p2,mC2(i1),mN2(i2)),dp) &
       & - 4._dp * Real( Conjg(coupCR) * coupCL,dp ) * mC(i1) * mN(i2)        &
       &         * Real( B0(p2, mC2(i1), mN2(i2)),dp)
    If (WriteOut) Write(ErrCan,*) "C N",i1,i2,sumI
    res = res + sumI
   End Do
  End Do

  res = oo16pi2 * res 

  Iname = Iname - 1

 End Subroutine PiChargedScalar



#ifdef GENERATIONMIXING
 Subroutine PiPseudoScalar(p2, gU1, gSU2, vevs_DR, mP02, RP0             &
     & , Y_d, uD_L, uD_R, Y_u, uU_L, uU_R, Y_l, uL_L, uL_R               &
     & , A_d, A_u, A_l, mu, mUSquark2, RUSquark, mDSquark2, RDSquark     &
     & , mSlepton2, RSlepton, mSneutrino2, Rsneutrino, mSpm2, RSpm       &
     & , mC, mC2, U, V, mS02, RS0, mN, mN2, N, mZ2, mW2, res)
#else
 Subroutine PiPseudoScalar(p2, gU1, gSU2, vevs_DR, mP02, RP0, Y_d, Y_u     &
     & , Y_l, A_d, A_u, A_l, mu, mUSquark2, RUSquark, mDSquark2, RDSquark  &
     & , mSlepton2, RSlepton, mSneutrino2, mSpm2, RSpm, mC, mC2, U, V      &
     & , mS02, RS0, mN, mN2, N, mZ2, mW2, res)
#endif
 !-------------------------------------------------------------------------
 ! calculates the 1-loop self energies of the pseudoscalar higgs
 ! The formulae of J. Bagger et al, Nucl. Phys. B491 are used.
 ! Note, that gSU3 has to be set to 0 for sleptons
 ! the input is:
 ! - p2 ......... the outer momentum squared 
 ! - gU1 ........ the U(1) gauge coupling
 ! - gSU2 ....... the SU(2) gauge coupling
 ! 
 ! the output is:
 !  res  
 ! written by Werner Porod, 30.11.99
 ! last change: 30.11.99
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: p2, vevs_DR(2), gU1, gSU2, mUSquark2(6)       &
     & , mDSquark2(6), mSlepton2(:), mSneutrino2(:), mP02(:), RP0(:,:)  &
     & , mSpm2(:), mC(:), mC2(:), mS02(:), RS0(:,:), mN(:), mN2(:), mZ2, mW2
  Complex(dp), Intent(in) :: Y_d(3,3), Y_u(3,3), Y_l(3,3), RUSquark(6,6)  &
     & , RDSquark(6,6), RSlepton(:,:), A_d(3,3), A_u(3,3), A_l(3,3), mu   &
     & , RSpm(:,:), U(:,:), V(:,:), N(:,:)
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: RSneutrino(3,3), uU_L(3,3), uU_R(3,3)  &
    & , uD_L(3,3), uD_R(3,3), uL_L(3,3), uL_R(3,3)
#endif
  Complex(dp), Intent(out) :: res

  Integer :: i1, i2, i3, i_gen
  Complex(dp) :: sumI, coupC, coupC2, Rsf(2,2), Af, A0m2, bi(1), Zero3(3,3) &
     & , yuk
  Real(dp) :: B0m2, e_u, e_d, G0m2, coup, F0m2, cosW_DR
  Logical :: WriteOut

  Iname = Iname + 1
  NameOfUnit(Iname) = "PiPseudoScalar"

  !----------------------------------
  ! In case I need the contributions
  !----------------------------------
  If ((WriteOneLoopContributions.Eq.3).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions in PiPseudoScalar:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If
 !-----------------------
 ! intialisation
 !-----------------------
  res = (0._dp, 0._dp) 

  e_u = 2._dp / 3._dp
  e_d = -1._dp / 3._dp

  cosW_DR = gSU2 / Sqrt(gSU2**2 + gU1**2)

  bi(1) = mu

  Zero3 = 0._dp

  Do i1=1,3
  !---------------
  ! u-quarks
  !---------------
#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Call CoupFermionPseudoScalar(i1, i1, 2, 0.5_dp, Y_u, uU_L, uU_R, RP0 &
                                 &, coupC, coupC2)
   Else
#endif
    Call CoupFermionPseudoScalar(2, 0.5_dp, Y_u(i1,i1), RP0, coupC, coupC2)
#ifdef GENERATIONMIXING
   End If
#endif
   B0m2 = Real( B0(p2, mf_u2(i1),mf_u2(i1) ),dp )
   A0m2 = A0( mf_u2(i1) )
   sumI = 6._dp * Abs( coupC )**2 * (p2 * B0m2 - 2._dp * A0m2)
   If (WriteOut) Write(ErrCan,*) "u-quarks ",i1,sumI,coupc
   res = res + sumI
  !---------------
  ! d-quarks
  !---------------
#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Call CoupFermionPseudoScalar(i1, i1, 2, -0.5_dp, Y_d, uD_L, uD_R, RP0 &
                                 &, coupC, coupC2)
   Else
#endif
    Call CoupFermionPseudoScalar(2, -0.5_dp, Y_d(i1,i1), RP0, coupC, coupC2)
#ifdef GENERATIONMIXING
   End If
#endif
   B0m2 = Real( B0(p2, mf_d2(i1),mf_d2(i1) ),dp )
   A0m2 = A0( mf_d2(i1) )
   sumI = 6._dp * Abs( coupC )**2 * (p2 * B0m2 - 2._dp * A0m2)
   If (WriteOut) Write(ErrCan,*) "d-quarks ",i1,sumI,coupc
   res = res + sumI
  !---------------
  ! leptons
  !---------------
#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Call CoupFermionPseudoScalar(i1, i1, 2, -0.5_dp, Y_l, uL_L, uL_R, RP0 &
                                 &, coupC, coupC2)
   Else
#endif
    Call CoupFermionPseudoScalar(2, -0.5_dp, Y_l(i1,i1), RP0, coupC, coupC2)
#ifdef GENERATIONMIXING
   End If
#endif
   B0m2 = Real( B0(p2, mf_l2(i1),mf_l2(i1) ),dp )
   A0m2 = A0( mf_l2(i1) )
   sumI = 2._dp * Abs( coupC )**2 * (p2 * B0m2 - 2._dp * A0m2)
   If (WriteOut) Write(ErrCan,*) "leptons ",i1,sumI,coupc
   res = res + sumI
  End Do

#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,6
   !-------------
   ! u-squarks
   !-------------
    A0m2 = A0( mUSquark2(i1) )
    Call CoupPseudoscalarSfermion4(2, 2, i1, i1, RP0, 0.5_dp, e_u  &
                                  & , Y_u, RUSquark, gU1, gSU2, coupC)
    sumI = - 6._dp * coupC * A0m2
    If (WriteOut) Write(ErrCan,*) "u-squarks ",i1,sumI
    res = res + sumI

    Do i2=i1,6
     B0m2 = Real( B0(p2, mUSquark2(i1), mUSquark2(i2)),dp )
     Call CoupPseudoScalarSfermion3(2, i1, i2, RP0, 0.5_dp, Y_u, RUSquark &
                                    &, A_u, bi, coupC )
     sumI = 3._dp * Abs( coupC )**2 * B0m2
     If (i1.Ne.i2) sumI = 2._dp * sumI
     If (WriteOut) Write(ErrCan,*) "u-squarks ",i1,i2,sumI
     res = res + sumI
    End Do 
   !-------------
   ! d-squarks
   !-------------
    A0m2 = A0( mDSquark2(i1) )
    Call CoupPseudoscalarSfermion4(2, 2, i1, i1, RP0, -0.5_dp, e_d  &
                                  & , Y_d, RDSquark, gU1, gSU2, coupC)
    sumI = - 6._dp * coupC * A0m2
    If (WriteOut) Write(ErrCan,*) "d-squarks ",i1,sumI
    res = res + sumI

    Do i2=i1,6
     B0m2 = Real( B0(p2, mDSquark2(i1), mDSquark2(i2)),dp )
     Call CoupPseudoScalarSfermion3(2, i1, i2, RP0, -0.5_dp, Y_d, RDSquark  &
                                   &, A_d, bi, coupC )
     sumI = 3._dp * Abs( coupC )**2 * B0m2
     If (i1.Ne.i2) sumI = 2._dp * sumI
     If (WriteOut) Write(ErrCan,*) "d-squarks ",i1,i2,sumI
     res = res + sumI
    End Do 
   !-------------
   ! sleptons
   !-------------
    If (i1.Le.n_Slept) Then
     A0m2 = A0( mSlepton2(i1) )
     Call CoupPseudoscalarSfermion4(2, 2, i1, i1, RP0, -0.5_dp, -1._dp  &
                                  & , Y_l, RSlepton, gU1, gSU2, coupC)
     sumI = - 2._dp * coupC * A0m2
     If (WriteOut) Write(ErrCan,*) "sleptons ",i1,sumI
     res = res + sumI

     Do i3=i2,n_Slept
      B0m2 = Real( B0(p2, mSlepton2(i1), mSlepton2(i2)),dp )
      Call CoupPseudoScalarSfermion3(2, i1, i2, RP0, -0.5_dp, Y_l, RSlepton  &
                                    &, A_l, bi, coupC )
      sumI = Abs( coupC )**2 * B0m2
      If (i1.Ne.i2) sumI = 2._dp * sumI
      If (WriteOut) Write(ErrCan,*) "sleptons ",i1,i2,sumI
      res = res + sumI
     End Do 
    End If
   !-------------
   ! sneutrino
   !-------------
    If (i1.Le.n_Sneut) Then
     A0m2 = A0(mSneutrino2(i1))
     Call CoupPseudoscalarSfermion4(2, 2, i1, i1, RP0, 0.5_dp, 0._dp  &
                                  & , Zero3, RSneutrino, gU1, gSU2, coupC)
     sumI = - 2._dp * coupC * A0m2
     If (WriteOut) Write(ErrCan,*) "sneutrinos ",i1,sumI
     res = res + sumI
    End If
   End Do     ! i1

  Else
#endif

   Do i1=1,5,2
   !-------------
   ! u-squarks
   !-------------
    Rsf = RUSquark(i1:i1+1,i1:i1+1)
    i_gen = (i1+1) / 2
    Af = A_u(i_gen,i_gen)
    yuk = Y_u(i_gen,i_gen)
    Do i2=1,2
     A0m2 = A0( mUSquark2(i1-1+i2) )
     Call CoupPseudoscalarSfermion4(2, 2, i2, i2, RP0, 0.5_dp, e_u  &
                                  & , yuk, Rsf, gU1, gSU2, coupC)
     sumI = - 6._dp * coupC * A0m2
     If (WriteOut) Write(ErrCan,*) "u-squarks ",i1,i2,sumI
     res = res + sumI

     Do i3=i2,2
      Call CoupPseudoScalarSfermion3(2, i2, i3, RP0, 0.5_dp, Yuk, Rsf   &
                                    &, Af, bi, coupC )
      If (Abs(coupC).Gt.1.e-20_dp) Then
       B0m2 = Real( B0(p2, mUSquark2(i1-1+i2), mUSquark2(i1-1+i3)) ,dp)
       sumI = 3._dp * Abs( coupC )**2 * B0m2
       If (i2.Ne.i3) sumI = 2._dp * sumI
      Else
       sumI = 0._dp
      End If
      If (WriteOut) Write(ErrCan,*) "u-squarks ",i1,i2,i3,sumI
      res = res + sumI
     End Do
    End Do 
   !-------------
   ! d-squarks
   !-------------
    Rsf = RDSquark(i1:i1+1,i1:i1+1)
    Af = A_d(i_gen,i_gen)
    yuk = Y_d(i_gen,i_gen)
    Do i2=1,2
     A0m2 = A0( mDSquark2(i1-1+i2) )
     Call CoupPseudoscalarSfermion4(2, 2, i2, i2, RP0, -0.5_dp, e_d  &
                                  & , yuk, Rsf, gU1, gSU2, coupC)
     sumI = - 6._dp * coupC * A0m2
     If (WriteOut) Write(ErrCan,*) "d-squarks ",i1,i2,sumI
     res = res + sumI

     Do i3=i2,2
      Call CoupPseudoScalarSfermion3(2, i2, i3, RP0, -0.5_dp, Yuk, Rsf   &
                                    &, Af, bi, coupC )
      If (Abs(coupC).Gt.1.e-20_dp) Then
       B0m2 = Real( B0(p2, mDSquark2(i1-1+i2), mDSquark2(i1-1+i3)),dp )
       sumI = 3._dp * Abs( coupC )**2 * B0m2
       If (i2.Ne.i3) sumI = 2._dp * sumI
      Else
       sumI = 0._dp
      End If
      If (WriteOut) Write(ErrCan,*) "d-squarks ",i1,i2,i3,sumI
      res = res + sumI
     End Do
    End Do 
   !-------------
   ! sleptons
   !-------------
    If (i1.Lt.n_Slept) Then
     Rsf = RSlepton(i1:i1+1,i1:i1+1)
     Af = A_l(i_gen,i_gen)
     yuk = Y_l(i_gen,i_gen)
     Do i2=1,2
      A0m2 = A0( mSlepton2(i1-1+i2) )
      Call CoupPseudoscalarSfermion4(2, 2, i2, i2, RP0, -0.5_dp, -1._dp  &
                                   & , yuk, Rsf, gU1, gSU2, coupC)
      sumI = - 2._dp * coupC * A0m2
      If (WriteOut) Write(ErrCan,*) "sleptons ",i1,i2,sumI
      res = res + sumI

      Do i3=i2,2
       Call CoupPseudoScalarSfermion3(2, i2, i3, RP0, -0.5_dp, Yuk, Rsf   &
                                     &, Af, bi, coupC )
       If (Abs(coupC).Gt.1.e-20_dp) Then
        B0m2 = Real( B0(p2, mSlepton2(i1-1+i2), mSlepton2(i1-1+i3)),dp )
        sumI = Abs( coupC )**2 * B0m2
        If (i2.Ne.i3) sumI = 2._dp * sumI
       Else
        sumI = 0._dp
       End If
       If (WriteOut) Write(ErrCan,*) "sleptons ",i1,i2,i3,sumI
       res = res + sumI
      End Do
     End Do 
    End If
   !-------------
   ! sneutrino
   !-------------
    If (i_gen.Le.n_Sneut) Then
     Af = (0._dp,0._dp)
     yuk = Af
     A0m2 = A0(mSneutrino2(i_gen))
     Call CoupPseudoscalarSfermion4(2, 2, 1, 1, RP0, 0.5_dp, 0._dp  &
                                  & , yuk, id2C, gU1, gSU2, coupC)
     sumI = - 2._dp * coupC * A0m2
     If (WriteOut) Write(ErrCan,*) "sneutrinos ",i_gen,sumI
     res = res + sumI
    End If
   End Do     ! i1
#ifdef GENERATIONMIXING
  End If ! GenerationMixing
#endif
 !-------------------
 ! gauge bosons
 !-------------------
  SumI = gSU2**2 * ( 2._dp * A0(mW2) + A0(mZ2) / cosW_DR**2 )
  If (WriteOut) Write(ErrCan,*) "W, Z ",sumI
  res = res + sumI

 !-----------------
 ! charged Higgs
 !-----------------
  sumI = 0.5_dp * gSU2**2 * ( Real( Floop(p2,mSpm2(2),mW2),dp )  &
       &                    + mW2 * Real( B0(p2,mW2,mSpm2(2)), dp ) )
  If (WriteOut) Write(ErrCan,*) "W + S+ ",sumI
  res = res + sumI

  Do i1=1,n_Spm
   A0m2 = A0( mSpm2(i1) )
   Call CoupChargedScalarPseudoScalar4(i1, i1, 2, 2, RSpm, RP0, gU1, gSU2  &
                                     &, coupC)
   sumI = coupC * A0m2
   If (WriteOut) Write(ErrCan,*) "S+ ",i1,sumI
   res = res + sumI
  End Do
 !------------------
 ! charginos
 !------------------
  Do i1=1,n_char
   Call CoupCharginoPseudoScalar(i1, i1, 2, U, V, RP0, gSU2, coupC, coupC2)
   G0m2 = Real( Gloop(p2,mC2(i1),mC2(i1)),dp )
   B0m2 = Real( B0(p2,mC2(i1),mC2(i1)),dp )
   sumI = ( Abs(coupC)**2 + Abs(coupC2)**2 ) * G0m2               &
      & - 4._dp * Real( Conjg(coupC2) * coupC,dp ) * mC2(i1) * B0m2
   If (WriteOut) Write(ErrCan,*) "Charginos ",i1,i1,sumI
   res = res + sumI
   Do i2=i1+1,n_char
    Call CoupCharginoPseudoScalar(i1, i2, 2, U, V, RP0, gSU2, coupC, coupC2)
    G0m2 = Real( Gloop(p2,mC2(i1),mC2(i2)),dp )
    B0m2 = Real( B0(p2,mC2(i1),mC2(i2)),dp )
    sumI = 2._dp * ( Abs(coupC)**2 + Abs(coupC2)**2 ) * G0m2       &
       & - 8._dp * Real( Conjg(coupC2) * coupC,dp ) * mC(i1) * mC(i2) * B0m2
    If (WriteOut) Write(ErrCan,*) "Charginos ",i1,i2,sumI
    res = res + sumI
   End Do
  End Do
 !--------------------
 ! pseudoscalar Higgs
 !--------------------
  A0m2 = A0( mP02(1) ) 
  coup = (gSU2**2 + gU1**2) * ( 1.5_dp * (RP0(1,1)*RP0(1,2))**2 - 0.125_dp )
  sumI = coup * A0m2
  If (WriteOut) Write(ErrCan,*) "G0 ",sumI
  res = res + sumI

  A0m2 = A0( mP02(2) ) 
  coup = 0.375_dp * (gSU2**2 + gU1**2) * (RP0(1,1)**2 - RP0(1,2)**2 )**2 
  sumI = coup * A0m2
  If (WriteOut) Write(ErrCan,*) "A0 ",sumI
  res = res + sumI
 !--------------------
 ! neutral Higgs
 !--------------------
  Do i1=1,n_S0
   F0m2 = Real( Floop(p2, mS02(i1), mZ2),dp )
   Call CoupPseudoscalarScalarZ(2, i1, gSU2, cosW_DR, RP0, RS0, coupC)
   sumI = Abs( coupC )**2 * F0m2
   If (WriteOut) Write(ErrCan,*) "S0, Z ",i1,sumI
   res = res + sumI 

   Do i2=1,n_P0
    B0m2 = 2._dp * Real( B0(p2, mP02(i2), mS02(i1)),dp )
    Call CoupPseudoScalarScalar3(i2, 2, i1, RP0, RS0, gU1, gSU2, vevs_DR, coup)
    sumI = coup**2 * B0m2
    If (WriteOut) Write(ErrCan,*) "S0, P0 ",i1,i2,sumI
    res = res + sumI 
   End Do

   A0m2 = A0( mS02(i1) ) 
   Call CoupPseudoScalarScalar4(2, 2, i1, i1, RP0, RS0, gU1, gSU2, coup)
   sumI = -  coup * A0m2
   If (WriteOut) Write(ErrCan,*) "S0, P0 ",i1,sumI
   res = res + sumI 
  End Do
 !------------------
 ! neutralinos
 !------------------
  Do i1=1,n_neut
   Call CoupNeutralinoPseudoscalar(i1, i1, 2, N, RP0, gU1, gSU2, coupC,coupC2)
   G0m2 = Real( Gloop(p2, mN2(i1), mN2(i1)),dp )
   B0m2 = Real( B0(p2, mN2(i1), mN2(i1)),dp )
   sumI = 0.5_dp * ( Abs(coupC)**2 + Abs(coupC2)**2 ) * G0m2    &
      & - 2._dp * Real( Conjg(coupC2) * coupC,dp ) * mN2(i1) * B0m2
   If (WriteOut) Write(ErrCan,*) "Neutralinos ",i1,i1,sumI
   res = res + sumI
   Do i2=i1+1,n_neut
    Call CoupNeutralinoPseudoscalar(i1, i2, 2, N, RP0, gU1, gSU2, coupC,coupC2)
    G0m2 = Real( Gloop(p2, mN2(i1), mN2(i2)),dp )
    B0m2 = Real( B0(p2, mN2(i1), mN2(i2)),dp )
    sumI = ( Abs(coupC)**2 + Abs(coupC2)**2 ) * G0m2      &
       & - 4._dp * Real( Conjg(coupC2) * coupC,dp ) * mN(i1) * mN(i2)  * B0m2
    If (WriteOut) Write(ErrCan,*) "Neutralinos ",i1,i2,sumI
    res = res + sumI
   End Do
  End Do

  res = oo16pi2 * res 

  Iname = Iname - 1

 End Subroutine PiPseudoScalar


 Subroutine PiPseudoScalar_NMSSM(p2, c_DDP0_L, c_UUP0_L, c_LLP0_L             &
     & , mDSquark2, c_SdSdP0P0, mUSquark2, c_SuSuP0P0, mSlepton2, c_SlSlP0P0  &
     & , mSneutrino2, c_SnSnP0P0, c_SdSdP0, c_SuSuP0, c_SlSlP0, c_SnSnP0      &
     & , mS02, mP02, mSpm2, c_SpmP0W, c_P0S0Z, c_P0P0WW, c_P0P0ZZ, c_SpSmP0P0 &
     & , mC, mC2, c_CCP0_L, c_CCP0_R, c_P0P0S0S0, c_S0P0P0, c_P04, mN, mN2    &
     & , c_NNP0_L, c_NNP0_R, mZ2, mW2, WriteOut, res)
 !-------------------------------------------------------------------------
 ! calculates the 1-loop self energies of the pseudoscalar higgs
 ! The formulae of J. Bagger et al, Nucl. Phys. B491 are used.
 ! Note, that gSU3 has to be set to 0 for sleptons
 ! the input is:
 ! - p2 ......... the outer momentum squared 
 ! - gU1 ........ the U(1) gauge coupling
 ! - gSU2 ....... the SU(2) gauge coupling
 ! the output is:
 ! res
 ! written 03.11.2008
 !-----------------------------------------------------------------------
 Implicit None

 Real(dp), Intent(in) :: p2, mUSquark2(6), mDSquark2(6), mSneutrino2(:)    &
     &, mS02(:), mP02(:), mSpm2(:), c_P0S0Z(:,:), c_P0P0WW(:), c_P0P0ZZ(:) &
     &, mSlepton2(:), mC(:), mC2(:), c_P0P0S0S0(:,:,:), c_S0P0P0(:,:,:)    &
     &, c_P04(:,:,:), mN(:), mN2(:), mZ2, mW2
  Complex(dp), Intent(in) :: c_DDP0_L(:), c_UUP0_L(:), c_LLP0_L(:)           &
     &, c_SdSdP0P0(:,:,:), c_SuSuP0P0(:,:,:), c_SlSlP0P0(:,:,:)              &
     &, c_SnSnP0P0(:,:,:), c_SdSdP0(:,:,:), c_SuSuP0(:,:,:), c_SlSlP0(:,:,:) &
     &, c_SnSnP0(:,:,:), c_SpmP0W(:,:), c_SpSmP0P0(:,:,:), c_CCP0_L(:,:,:)   &
     &, c_CCP0_R(:,:,:), c_NNP0_L(:,:,:), c_NNP0_R(:,:,:)
  Logical, Intent(in) :: WriteOut
  Complex(dp), Intent(inout) :: res(:,:)

  Integer :: i1, i2, i3, i4, i5, i_sf1, i_sf2
  Complex(dp) :: sumI(n_P0,n_P0)

  Real(dp) :: A0m2, B0m2, F0m2, G0m2

  Iname = Iname + 1
  NameOfUnit(Iname) = "PiPseudoScalar"

  !-----------------------
  ! initialisation
  !-----------------------
  res = 0._dp
  !-------------------
  ! fermions
  !-------------------

  sumI = 0._dp
  Do i1=1,3
   sumI(1,1) = 6._dp * Abs( c_DDP0_L(i1) )**2                              &
           &   * ((p2 * Real(B0(p2, mf_d2(i1), mf_d2(i1)),dp)              &
           &     - 2._dp * Real( A0(mf_d2(i1)),dp )))                      &
           & +  2._dp * Abs( c_LLP0_L(i1) )**2                             &
           &   * ((p2 * Real(B0(p2, mf_l2(i1), mf_l2(i1)),dp)              &
           &     - 2._dp * Real( A0( mf_l2(i1) ),dp )))
   sumI(2,2) = 6._dp * Abs( c_UUP0_L(i1) )**2                              &
           &   * (p2 * Real(B0(p2, mf_u2(i1), mf_u2(i1)),dp)               &
           &     - 2._dp * Real( A0(mf_u2(i1)),dp) )
   If (WriteOut) Write(ErrCan,*) "fermions",i1,sumI(1,1),sumI(2,2)
   res = res + sumI
  End Do

  !--------------------
  ! sfermions
  !--------------------
  sumI = 0._dp
  Do i1=1,6
   A0m2 = -6._dp * A0( mDsquark2(i1) )
   Do i2=1,n_P0
    Do i3=i2,n_P0
       sumI(i2,i3) = c_SdSdP0P0(i1,i2,i3) * A0m2
    End Do
   End Do 
   res = res + sumI
If (WriteOut) Write(ErrCan,*) "Sd Sd P0 P0",i1,sumI(1,1),sumI(1,2)&
      &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
  End Do

  sumI = 0._dp
  Do i1=1,6
   A0m2 = -6._dp * A0( mUsquark2(i1) )
   Do i2=1,n_P0 
    Do i3=i2,n_P0
      sumI(i2,i3) = c_SuSuP0P0(i1,i2,i3) * A0m2
    End Do
   End Do
   res = res + sumI
If (WriteOut) Write(ErrCan,*) "Su Su P0 P0",i1,sumI(1,1),sumI(1,2)&
      &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
  End Do

  sumI = 0._dp
  Do i1=1,n_Slept
   A0m2 = -2._dp * A0( mSlepton2(i1) )
   Do i2=1,n_P0 
    Do i3=i2,n_P0
     sumI(i2,i3) = c_SlSlP0P0(i1,i2,i3) * A0m2
    End Do
   End Do
   res = res + sumI
If (WriteOut) Write(ErrCan,*) "Sl Sl P0 P0",i1,sumI(1,1),sumI(1,2)&
      &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
  End Do

  sumI = 0._dp
  Do i1=1,n_Sneut
   A0m2 = -2._dp * A0( mSneutrino2(i1) )
   Do i2=1,n_P0
    Do i3=i2,n_P0
     sumI(i2,i3) = c_SnSnP0P0((i1-1)*2+1,i2,i3) * A0m2
    End Do
   End Do
   res = res + sumI
 If (WriteOut) Write(ErrCan,*) "Sn Sn P0 P0",i1,sumI(1,1),sumI(1,2)&
      &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
  End Do

 ! .not.GenerationMixing
   sumI = 0._dp
   Do i1=1,3
     i_sf1 = (i1-1)*2+1
      i_sf2 = (i1-1)*2+2
      B0m2 = 3._dp * Real( B0(p2, mDSquark2(i_sf1), mDSquark2(i_sf2)),dp )
      Do i4=1,n_P0
       Do i5=i4,n_P0
        sumI(i4,i5) = Conjg( c_SdSdP0(i_sf1,i_sf2,i4) ) &
                    &  * c_SdSdP0(i_sf1,i_sf2,i5) * B0m2
       End Do ! i5
      End Do ! i4
      res = res + sumI
      If (WriteOut) Write(ErrCan,*) "Sd Sd P0",i_sf1,i_sf2 &
          & ,sumI(1,1),sumI(1,2),sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)

      B0m2 = 3._dp * Real( B0(p2, mUSquark2(i_sf1), mUSquark2( i_sf2) ),dp )
      Do i4=1,n_P0
       Do i5=i4,n_P0
        sumI(i4,i5) = Conjg( c_SuSuP0(i_sf1,i_sf2,i4) ) &
                    &  * c_SuSuP0(i_sf1,i_sf2,i5) * B0m2
       End Do ! i5
      End Do ! i4
      res = res + sumI
      If (WriteOut) Write(ErrCan,*) "Su Su P0",i_sf1,i_sf2 &
          & ,sumI(1,1),sumI(1,2),sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)  
   End Do

   sumI = 0._dp
   Do i1=1,3
     i_sf1 = (i1-1)*2+1
      i_sf2 = (i1-1)*2+2
      B0m2 = Real( B0(p2, mSlepton2(i_sf1), mSlepton2( i_sf2) ),dp )
      Do i4=1,n_P0
       Do i5=i4,n_P0
        sumI(i4,i5) = Conjg( c_SlSlP0(i_sf1,i_sf2,i4) ) &
                    &  * c_SlSlP0(i_sf1,i_sf2,i5) * B0m2
       End Do ! i5
      End Do ! i4
      If (WriteOut) Write(ErrCan,*) "Sl Sl P0",i_sf1,i_sf2 &
          & ,sumI(1,1),sumI(1,2),sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
      res = res + sumI
   End Do

 !-------------------
 ! gauge bosons
 !-------------------

 sumI = 0._dp
  Do i2=2,n_Spm
   F0m2 = Real( Floop(p2, mSpm2(i2), mW2),dp )
   sumI(1,1) = Abs( c_SpmP0W(i2,1) )**2 * F0m2
   sumI(1,2) = Conjg( c_SpmP0W(i2,1) ) * c_SpmP0W(i2,2) * F0m2
   sumI(2,2) = Abs( c_SpmP0W(i2,2) )**2 * F0m2
   If (WriteOut) Write(ErrCan,*) "S+ P0 W",i2 ,sumI(1,1),sumI(1,2),sumI(2,2)
   res = res + sumI
  End Do
  
  sumI = 0._dp
  F0m2 = Real( Floop(p2,mW2,mW2),dp ) + 7._dp * mW2 * Real( B0(p2,mW2,mW2),dp )
  If (WriteOut) Write(ErrCan,*) "P0 W+ W-",sumI(1,1),sumI(1,2),sumI(2,2)
  res = res + sumI

  sumI = 0._dp
  Do i2=1,n_S0
   F0m2 = Real( Floop(p2, mS02(i2), mZ2),dp )
   sumI(1,1) = Abs( c_P0S0Z(i2,1) )**2 * F0m2
   sumI(1,2) = c_P0S0Z(i2,1) * c_P0S0Z(i2,2) * F0m2
   sumI(1,3) = c_P0S0Z(i2,1) * c_P0S0Z(i2,3) * F0m2
   sumI(2,2) = Abs( c_P0S0Z(i2,2) )**2 * F0m2
   sumI(2,3) = c_P0S0Z(i2,2) * c_P0S0Z(i2,3) * F0m2
   If (WriteOut) Write(ErrCan,*) "P0 S0 Z",i2 ,sumI(1,1),sumI(1,2),sumI(2,2)
   res = res + sumI
  End Do

  sumI = 0._dp
  Do i1=1,n_P0
   SumI(i1,i1) = c_P0P0WW(i1) * A0(mW2) +  c_P0P0ZZ(i1) * A0(mZ2)
   If (WriteOut) Write(ErrCan,*) "P0 P0 WW, ZZ",i1,sumI(i1,i1)
   res(i1,i1) = res(i1,i1) + sumI(i1,i1)
  End Do

 !-----------------
 ! charged Higgs
 !-----------------

  Do i1=1,n_Spm
   A0m2 = A0( mSpm2(i1) )
   Do i3=1,n_P0
    Do i4=i3,n_P0
     sumI(i3,i4) = c_SpSmP0P0(i1,i3,i4) * A0m2
    End Do
   End Do
   If (WriteOut) Write(ErrCan,*) &
           & "Sp Sm P0 P0",i1,sumI(1,1),sumI(1,2),sumI(2,2)
   res = res + sumI
  End Do

 !------------------
 ! charginos
 !------------------

  sumI = 0._dp
  Do i1=1,n_char
   G0m2 = Real( Gloop(p2, mC2(i1), mC2(i1)),dp )
   B0m2 = 4._dp * mC2(i1) * Real( B0(p2, mC2(i1), mC2(i1)), dp )
   Do i3=1,n_P0
    sumI(i3,i3) =                                                           &
       &  (Abs(c_CCP0_L(i1,i1,i3))**2 + Abs(c_CCP0_R(i1,i1,i3))**2) * G0m2  &
       &  - Real( Conjg(c_CCP0_R(i1,i1,i3)) * c_CCP0_L(i1,i1,i3),dp ) * B0m2
    Do i4=i3+1,n_P0
     sumI(i3,i4) = ( Conjg(c_CCP0_L(i1,i1,i3))*c_CCP0_L(i1,i1,i4)           &
            &   + Conjg(c_CCP0_R(i1,i1,i3))* c_CCP0_R(i1,i1,i4) ) * G0m2    &
            & - 0.5_dp *( Conjg(c_CCP0_R(i1,i1,i3)) * c_CCP0_L(i1,i1,i4)    &
            &           + Conjg(c_CCP0_L(i1,i1,i3)) * c_CCP0_R(i1,i1,i4))*B0m2
    End Do
   End Do
   If (WriteOut) Write(ErrCan,*) "C C P0",i1,i1,sumI(1,1),sumI(1,2),sumI(1,3) &
       &,sumI(2,2),sumI(2,3),sumI(3,3) 
   res = res + sumI

   sumI = 0._dp
   Do i2=i1+1,n_char
    G0m2 = 2._dp * Real( Gloop(p2, mC2(i1), mC2(i2)),dp )
    B0m2 = 8._dp * mC(i1) * mC(i2) * Real( B0(p2, mC2(i1), mC2(i2)),dp )

    Do i3=1,n_s0
     sumI(i3,i3) =                                                           &
          &  (Abs(c_CCP0_L(i1,i2,i3))**2 +Abs(c_CCP0_R(i1,i2,i3))**2) * G0m2 &
          & - Real( Conjg(c_CCP0_R(i1,i2,i3)) * c_CCP0_L(i1,i2,i3),dp ) * B0m2
    Do i4=i3+1,n_P0     
     sumI(i3,i4) = ( Conjg(c_CCP0_L(i1,i2,i3))*c_CCP0_L(i1,i2,i4)          &
            &   + Conjg(c_CCP0_R(i1,i2,i3))* c_CCP0_R(i1,i2,i4) ) * G0m2   &
            & - 0.5_dp *( Conjg(c_CCP0_R(i1,i2,i3)) * c_CCP0_L(i1,i2,i4)   &
            &         + Conjg(c_CCP0_L(i1,i2,i3)) * c_CCP0_R(i1,i2,i4)) * B0m2
     End Do ! i4
    End Do ! i3
    If (WriteOut) Write(ErrCan,*) "C C P0",i1,i2,sumI(1,1),sumI(1,2),sumI(1,3)&
       &,sumI(2,2),sumI(2,3),sumI(3,3) 
    res = res + sumI
   End Do
  End Do

 !--------------------
 ! pseudoscalar Higgs
 !--------------------

 sumI = 0._dp
  Do i1=1,n_P0
   A0m2 = A0( mP02(i1) )
   Do i2=1,n_P0
    Do i3=i2,n_p0
    sumI(i2,i3) = c_P04(i1,i2,i3) * A0m2
    End Do
   End Do
   res = res + sumI
 If (WriteOut) Write(ErrCan,*) "P0 P0 P0 P0",i1,sumI(1,1),sumI(1,2),sumI(1,3)&
       &,sumI(2,2),sumI(2,3),sumI(3,3) 
  End Do

 !--------------------
 ! neutral Higgs
 !--------------------

 sumI = 0._dp
  Do i1=1,n_S0
   A0m2 = - 0.5_dp * A0( mS02(i1) )
   Do i2=1,n_P0
    Do i3=i2,n_P0
    sumI(i2,i3) = c_P0P0S0S0(i2,i3,i1) * A0m2
    If (WriteOut) Write(ErrCan,*) "P0 P0 S0 S0",i1,i2,sumI(1,1),sumI(1,2)&
           &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
    End Do ! i3
   End Do ! i2
   res = res + sumI
If (WriteOut) Write(ErrCan,*) "P0 P0 S0 S0",i1,i2,sumI(1,1),sumI(1,2)&
           &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
  End Do


  sumI = 0._dp
  Do i1=1,n_S0
   Do i2=1,n_P0
    B0m2 = 0.5_dp * Real( B0(p2,mS02(i1), mP02(i2)),dp )
    Do i3=1,n_S0
     sumI(i3,i3) = c_S0P0P0(i3,i1,i2)**2 * B0m2
     Do i4=i3+1,n_S0
      sumI(i3,i4) = c_S0P0P0(i3,i1,i2) * c_S0P0P0(i4,i1,i2) * B0m2
     End Do ! i4 
    End Do ! i3
    res = res + sumI
If (WriteOut) Write(ErrCan,*) "P0 P0 S0",i1,i2,sumI(1,1),sumI(1,2)&
           &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
   End Do
  End Do

 !------------------
 ! neutralinos
 !------------------

   sumI = 0._dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    G0m2 = 0.5_dp * Real( Gloop(p2,mN2(i1),mN2(i2)),dp )
    B0m2 = 2._dp * mN(i1) * mN(i2) * Real( B0(p2,mN2(i1),mN2(i2)),dp )
    Do i3=1,n_P0
     sumI(i3,i3) = (Abs(c_NNP0_L(i1,i2,i3))**2 + Abs(c_NNP0_R(i1,i2,i3))**2) &
              &   * G0m2                                                     &
              & - Real( Conjg(c_NNP0_R(i1,i2,i3)) * c_NNP0_L(i1,i2,i3),dp )   &
              &   * B0m2
     Do i4=i3+1,n_P0
      sumI(i3,i4) = ( Conjg(c_NNP0_L(i1,i2,i3)) * c_NNP0_L(i1,i2,i4)     &
               &    + Conjg(c_NNP0_R(i1,i2,i3)) * c_NNP0_R(i1,i2,i4) )   &
               &   * G0m2                                                &
               &  - ( Conjg(c_NNP0_R(i1,i2,i3)) * c_NNP0_L(i1,i2,i4)     &
               &    + Conjg(c_NNP0_L(i1,i2,i3)) * c_NNP0_R(i1,i2,i4) )   &
               &   * 0.5_dp *  B0m2
     End Do !i4
     res(i3,i3:n_P0) = res(i3,i3:n_P0) + sumI(i3,i3:n_P0)
    End Do ! i3
If (WriteOut) Write(ErrCan,*) "N N P0",i1,i2,sumI(1,1),sumI(1,2)&
           &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
   End Do ! i2
  End Do ! i1 

  Do i2=1,n_S0
   Do i1=2,n_S0
     res(i1,i2) = Conjg( res(i2,i1) )
   End Do
  End Do
  res = oo16pi2 * res

 End Subroutine PiPseudoScalar_NMSSM


 Subroutine PiScalar(p2, c_DDS0_L, c_UUS0_L, c_LLS0_L                         &
     & , mDSquark2, c_SdSdS0S0, mUSquark2, c_SuSuS0S0, mSlepton2, c_SlSlS0S0  &
     & , mSneutrino2, c_SnSnS0S0, c_SdSdS0, c_SuSuS0, c_SlSlS0, c_SnSnS0      &
     & , mS02, mP02, mSpm2, c_SpmS0W, c_S0WW, c_P0S0Z, c_S0ZZ, c_S0S0WW       &
     & , c_S0S0ZZ, c_SpSmS0S0, c_SpSmS0, mC, mC2, c_CCS0_L, c_CCS0_R          &
     & , c_P0P0S0S0, c_P0S0S0, c_S03, c_S04, mN, mN2, c_NNS0_L, c_NNS0_R      &
     & , mZ2, mW2, WriteOut, res)
 !-------------------------------------------------------------------------
 ! calculates the 1-loop self energies of higgs bosons
 ! The formulae of J. Bagger et al, Nucl. Phys. B491 are used.
 ! Note, that gSU3 has to be set to 0 for sleptons
 ! the input is:
 ! - p2 ......... the outer momentum squared 
 ! 
 ! the output is:
 !  res(i,j) 
 ! written by Werner Porod, 1.12.99
 ! 10.10.01: portation to f90
 ! 11.06.03: implementing alpha_b alpha_t + alpha^2_b corrections
 !           based on A. Dedes, G. Degrassi and P. Slavich, hep-ph/0305127
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: p2, mUSquark2(6), mDSquark2(6), mSneutrino2(:)    &
     & , mS02(:), mP02(:), mSpm2(:), c_S0WW(:), c_P0S0Z(:,:), c_S0ZZ(:)     &
     & , c_S0S0WW(:), c_S0S0ZZ(:), mSlepton2(:), mC(:), mC2(:)              &
     & , c_P0P0S0S0(:,:,:), c_P0S0S0(:,:,:), c_S03(:,:,:), c_S04(:,:,:)     &
     & , mN(:), mN2(:), mZ2, mW2
  Complex(dp), Intent(in) :: c_DDS0_L(:), c_UUS0_L(:), c_LLS0_L(:)          &
     & , c_SdSdS0S0(:,:), c_SuSuS0S0(:,:), c_SlSlS0S0(:,:), c_SnSnS0S0(:,:) &
     & , c_SdSdS0(:,:,:), c_SuSuS0(:,:,:), c_SlSlS0(:,:,:), c_SnSnS0(:,:,:) &
     & , c_SpmS0W(:,:), c_SpSmS0S0(:,:,:), c_SpSmS0(:,:,:), c_CCS0_L(:,:,:) &
     & , c_CCS0_R(:,:,:), c_NNS0_L(:,:,:), c_NNS0_R(:,:,:)
  Logical, Intent(in) :: WriteOut
  Complex(dp), Intent(inout) :: res(:,:)

  Integer :: i1, i2, i3, i4, i5, i_sf1, i_sf2
  Complex(dp) :: sumI(n_S0,n_S0 )

  Real(dp) :: A0m2, B0m2, F0m2, G0m2

  !-----------------------
  ! initialisation
  !-----------------------
  res = 0._dp
  !-------------------
  ! fermions
  !-------------------
  sumI = 0._dp
  Do i1=1,3
   sumI(1,1) = 6._dp * Abs( c_DDS0_L(i1) )**2                              &
           &   * ( (p2 - 4._dp * mf_d2(i1) )                               &
           &           * Real(B0(p2, mf_d2(i1), mf_d2(i1)),dp)             &
           &     - 2._dp * Real( A0(mf_d2(i1)),dp ) )                      &
           & +  2._dp * Abs( c_LLS0_L(i1) )**2                             &
           &   * ( (p2 - 4._dp * mf_l2(i1) )                               &
           &            * Real( B0(p2, mf_l2(i1), mf_l2(i1)),dp)           &
           &     - 2._dp * Real( A0( mf_l2(i1) ),dp ) )         
   sumI(2,2) = 6._dp * Abs( c_UUS0_L(i1) )**2                              &
           &   * ( (p2 - 4._dp * mf_u2(i1) )                               &
           &         * Real(B0(p2, mf_u2(i1), mf_u2(i1)),dp)               &
           &     - 2._dp * Real(A0(mf_u2(i1)),dp) )
   If (WriteOut) Write(ErrCan,*) "fermions",i1,sumI(1,1),sumI(2,2)
   res = res + sumI
  End Do

  !--------------------
  ! sfermions
  !--------------------
  Do i1=1,6
   A0m2 = -6._dp * A0( mDsquark2(i1) )
   Forall(i2=1:n_S0) sumI(i2,i2) = c_SdSdS0S0(i1,i2) * A0m2
   If (WriteOut) Write(ErrCan,*) "Sd Sd S0 S0",i1,sumI(1,1),sumI(2,2)
   res = res + sumI

   A0m2 = -6._dp * A0( mUsquark2(i1) )
   Forall(i2=1:n_S0) sumI(i2,i2) = c_SuSuS0S0(i1,i2) * A0m2
   If (WriteOut) Write(ErrCan,*) "Su Su S0 S0",i1,sumI(1,1),sumI(2,2)
   res = res + sumI
  End Do

  Do i1=1,n_Slept
   A0m2 = -2._dp * A0( mSlepton2(i1) )
   Forall(i2=1:n_S0) sumI(i2,i2) = c_SlSlS0S0(i1,i2) * A0m2
   If (WriteOut) Write(ErrCan,*) "Sl Sl S0 S0",i1,sumI(1,1),sumI(2,2)
   res = res + sumI
  End Do

  Do i1=1,n_Sneut
   A0m2 = -2._dp * A0( mSneutrino2(i1) )
   Forall(i2=1:n_S0) sumI(i2,i2) = c_SnSnS0S0(i1,i2) * A0m2
   If (WriteOut) Write(ErrCan,*) "Sn Sn S0 S0",i1,sumI(1,1),sumI(2,2)
   res = res + sumI
  End Do

#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,6
    Do i2=1,6
     B0m2 = 3._dp * Real( B0(p2, mDsquark2(i1), mDsquark2(i2)),dp )
     Do i3=1,n_s0
      Do i4=i3,n_s0
       sumI(i3,i4) = Conjg( c_SdSdS0(i1,i2,i3) ) * c_SdSdS0(i1,i2,i4) * B0m2
      End Do
     End Do
     If (WriteOut) Write(ErrCan,*) &
              & "Sd Sd S0",i1,i2,sumI(1,1),sumI(1,2),sumI(2,2)
     res = res + sumI

     B0m2 = 3._dp * Real( B0(p2, mUsquark2(i1), mUsquark2(i2)),dp )
     Do i3=1,n_s0
      Do i4=i3,n_s0
       sumI(i3,i4) = Conjg( c_SuSuS0(i1,i2,i3) ) * c_SuSuS0(i1,i2,i4) * B0m2
      End Do
     End Do
     If (WriteOut) Write(ErrCan,*) &
             &  "Su Su S0",i1,i2,sumI(1,1),sumI(1,2),sumI(2,2)
     res = res + sumI
    End Do
   End Do  

   Do i1=1,n_Slept
    Do i2=1,n_Slept
     B0m2 = Real( B0(p2, mSlepton2(i1), mSlepton2(i2)),dp )
     Do i3=1,n_s0
      Do i4=i3,n_s0
       sumI(i3,i4) = Conjg( c_SlSlS0(i1,i2,i3) ) * c_SlSlS0(i1,i2,i4) * B0m2
      End Do
     End Do
     If (WriteOut) Write(ErrCan,*) &
           & "Sl Sl S0",i1,i2,sumI(1,1),sumI(1,2),sumI(2,2)
     res = res + sumI
    End Do
   End Do  

   Do i1=1,n_Sneut
    Do i2=1,n_Sneut
     B0m2 = Real( B0(p2, mSneutrino2(i1), mSneutrino2(i2)),dp )
     Do i3=1,n_s0
      Do i4=i3,n_s0
       sumI(i3,i4) = Conjg( c_SnSnS0(i1,i2,i3) ) * c_SnSnS0(i1,i2,i4) * B0m2
      End Do
     End Do
     If (WriteOut) Write(ErrCan,*) &
            & "Sn Sn S0",i1,i2,sumI(1,1),sumI(1,2),sumI(2,2)
     res = res + sumI
    End Do
   End Do  

  Else ! .not.GenerationMixing
#endif

   Do i1=1,3
    Do i2=1,2
     i_sf1 = (i1-1)*2+i2
     Do i3=1,2
      i_sf2 = (i1-1)*2+i3
      B0m2 = 3._dp * Real( B0(p2, mDSquark2(i_sf1), mDSquark2( i_sf2) ),dp )
      Do i4=1,n_s0
       Do i5=i4,n_s0
        sumI(i4,i5) = Conjg( c_SdSdS0(i_sf1,i_sf2,i4) ) &
                    &  * c_SdSdS0(i_sf1,i_sf2,i5) * B0m2
       End Do
      End Do
      If (WriteOut) Write(ErrCan,*) "Sd Sd S0",i_sf1,i_sf2 &
                              & ,sumi(1,1),sumi(1,2),sumi(2,2)
      res = res + sumI

      B0m2 = 3._dp * Real( B0(p2, mUSquark2(i_sf1), mUSquark2( i_sf2) ),dp )
      Do i4=1,n_s0
       Do i5=i4,n_s0
        sumI(i4,i5) = Conjg( c_SuSuS0(i_sf1,i_sf2,i4) ) &
                    &  * c_SuSuS0(i_sf1,i_sf2,i5) * B0m2
       End Do
      End Do
      If (WriteOut) Write(ErrCan,*) "Su Su S0",i_sf1,i_sf2 &
                              & ,sumI(1,1),sumI(1,2),sumI(2,2)
      res = res + sumI
     End Do
    End Do  
   End Do

   Do i1=1,n_sneut
    Do i2=1,2
     i_sf1 = (i1-1)*2+i2
     Do i3=1,2
      i_sf2 = (i1-1)*2+i3
      B0m2 = Real( B0(p2, mSlepton2(i_sf1), mSlepton2( i_sf2) ),dp )
      Do i4=1,n_s0
       Do i5=i4,n_s0
        sumI(i4,i5) = Conjg( c_SlSlS0(i_sf1,i_sf2,i4) ) &
                    &  * c_SlSlS0(i_sf1,i_sf2,i5) * B0m2
       End Do
      End Do
      If (WriteOut) Write(ErrCan,*) "Sl Sl S0",i_sf1,i_sf2 &
                              & ,sumI(1,1),sumI(1,2),sumI(2,2)
      res = res + sumI
     End Do
    End Do

    B0m2 = Real( B0(p2, mSneutrino2(i1), mSneutrino2( i1) ),dp )
    Do i4=1,n_s0
     Do i5=i4,n_s0
      sumI(i4,i5) = Conjg( c_SnSnS0(i1,i1,i4) ) * c_SnSnS0(i1,i1,i5) * B0m2
     End Do
    End Do
    If (WriteOut) Write(ErrCan,*) &
          & "Sn Sn S0",i1,i1 ,sumI(1,1),sumI(1,2),sumI(2,2)
    res = res + sumI
   End Do
#ifdef GENERATIONMIXING
  End If
#endif
 !-----------------
 ! gauge bosons
 !-----------------
  Do i2=2,n_Spm
   F0m2 = Real( Floop(p2, mSpm2(i2), mW2),dp )
   sumI(1,1) = Abs( c_SpmS0W(i2,1) )**2 * F0m2
   sumI(1,2) = Conjg( c_SpmS0W(i2,1) ) * c_SpmS0W(i2,2) * F0m2
   sumI(2,2) = Abs( c_SpmS0W(i2,2) )**2 * F0m2
   If (WriteOut) Write(ErrCan,*) "S+ S0 W",i1 ,sumI(1,1),sumI(1,2),sumI(2,2)
   res = res + sumI
  End Do

  F0m2 = Real( Floop(p2,mW2,mW2),dp ) + 7._dp * mW2 * Real( B0(p2,mW2,mW2),dp )
  sumI(1,1) = Abs( c_S0WW(1) )**2 * F0m2
  sumI(1,2) = c_S0WW(1) * c_S0WW(2) * F0m2
  sumI(2,2) = Abs( c_S0WW(2) )**2 * F0m2
  If (WriteOut) Write(ErrCan,*) "S0 W+ W-",sumI(1,1),sumI(1,2),sumI(2,2)
  res = res + sumI

  Do i2=2,n_P0
   F0m2 = Real( Floop(p2, mP02(i2), mZ2),dp )
   sumI(1,1) = Abs( c_P0S0Z(i2,1) )**2 * F0m2
   sumI(1,2) = c_P0S0Z(i2,1) * c_P0S0Z(i2,2) * F0m2
   sumI(2,2) = Abs( c_P0S0Z(i2,2) )**2 * F0m2
   If (WriteOut) Write(ErrCan,*) "P0 S0 Z",i1 ,sumI(1,1),sumI(1,2),sumI(2,2)
   res = res + sumI
  End Do

  F0m2 = Real( Floop(p2,mZ2,mZ2),dp ) + 7._dp * mZ2 * Real( B0(p2,mZ2,mZ2),dp )
  sumI(1,1) = Abs( c_S0ZZ(1) )**2 * F0m2
  sumI(1,2) = c_S0ZZ(1) * c_S0ZZ(2) * F0m2
  sumI(2,2) = Abs( c_S0ZZ(2) )**2 * F0m2
  If (WriteOut) Write(ErrCan,*) "S0 Z Z",sumI(1,1),sumI(1,2),sumI(2,2)
  res = res + sumI

  Do i1=1,n_S0
   SumI(i1,i1) = c_S0S0WW(i1) * A0(mW2) +  c_S0S0ZZ(i1) * A0(mZ2)
   If (WriteOut) Write(ErrCan,*) "S0 S0 WW, ZZ",i1,sumI(i1,i1)
   res(i1,i1) = res(i1,i1) + sumI(i1,i1)
  End Do
 !-----------------
 ! charged Higgs
 !-----------------
  Do i1=1,n_Spm
   A0m2 = A0( mSpm2(i1) )
   Do i3=1,n_s0
    Do i4=i3,n_s0
     sumI(i3,i4) = c_SpSmS0S0(i1,i3,i4) * A0m2
    End Do
   End Do
   If (WriteOut) Write(ErrCan,*) &
           & "Sp Sm S0 S0",i1,sumI(1,1),sumI(1,2),sumI(2,2)
   res = res + sumI

   Do i2=1,n_Spm
    B0m2 = Real( B0(p2, mSpm2(i1), mSpm2(i2)),dp )
    Do i3=1,n_s0
     Do i4=i3,n_s0
      sumI(i3,i4) = Conjg( c_SpSmS0(i1,i2,i3) ) * c_SpSmS0(i1,i2,i4) * B0m2
     End Do
    End Do
    If (WriteOut) Write(ErrCan,*) &
            & "Sp Sm S0",i1,i2,sumI(1,1),sumI(1,2),sumI(2,2)
    res = res + sumI
   End Do
  End Do
 !------------------------------
 ! Charginos
 !------------------------------
  Do i1=1,n_char
   G0m2 = Real( Gloop(p2, mC2(i1), mC2(i1)),dp )
   B0m2 = 4._dp * mC2(i1) * Real( B0(p2, mC2(i1), mC2(i1)), dp )
   Do i3=1,n_s0
    sumI(i3,i3) =                                                           &
        &  (Abs(c_CCS0_L(i1,i1,i3))**2 + Abs(c_CCS0_R(i1,i1,i3))**2) * G0m2 &
        &  - Real( Conjg(c_CCS0_R(i1,i1,i3)) * c_CCS0_L(i1,i1,i3),dp ) * B0m2
    Do i4=i3+1,n_S0
     sumI(i3,i4) = ( Conjg(c_CCS0_L(i1,i1,i3))*c_CCS0_L(i1,i1,i4)              &
            &   + Conjg(c_CCS0_R(i1,i1,i3))* c_CCS0_R(i1,i1,i4) ) * G0m2       &
            & - 0.5_dp *( Conjg(c_CCS0_R(i1,i1,i3)) * c_CCS0_L(i1,i1,i4)       &
            &           + Conjg(c_CCS0_L(i1,i1,i3)) * c_CCS0_R(i1,i1,i4)) * B0m2
    End Do
   End Do
   If (WriteOut) Write(ErrCan,*) "C C S0",i1,i1,sumI(1,1),sumI(1,2),sumI(2,2)
   res = res + sumI

   Do i2=i1+1,n_char
    G0m2 = 2._dp * Real( Gloop(p2, mC2(i1), mC2(i2)),dp )
    B0m2 = 8._dp * mC(i1) * mC(i2) * Real( B0(p2, mC2(i1), mC2(i2)),dp )

    Do i3=1,n_s0
     sumI(i3,i3) =                                                           &
          &  (Abs(c_CCS0_L(i1,i2,i3))**2 +Abs(c_CCS0_R(i1,i2,i3))**2) * G0m2 &
          & - Real( Conjg(c_CCS0_R(i1,i2,i3)) * c_CCS0_L(i1,i2,i3), dp ) * B0m2
    Do i4=i3+1,n_S0     
     sumI(i3,i4) = ( Conjg(c_CCS0_L(i1,i2,i3))*c_CCS0_L(i1,i2,i4)          &
            &   + Conjg(c_CCS0_R(i1,i2,i3))* c_CCS0_R(i1,i2,i4) ) * G0m2   &
            & - 0.5_dp *( Conjg(c_CCS0_R(i1,i2,i3)) * c_CCS0_L(i1,i2,i4)   &
            &           + Conjg(c_CCS0_L(i1,i2,i3)) * c_CCS0_R(i1,i2,i4)) * B0m2
     End Do
    End Do
    If (WriteOut) Write(ErrCan,*) "C C S0",i1,i2,sumI(1,1),sumI(1,2),sumI(2,2)
    res = res + sumI
   End Do
  End Do
 !---------------------
 ! pseudoscalar Higgs
 !---------------------
  sumI = 0._dp
  Do i1=1,n_P0
   A0m2 = - A0( mP02(i1) )
   Do i2=1,n_S0
    sumI(i2,i2) = c_P0P0S0S0(i1,i2,i2) * A0m2
    If (WriteOut) Write(ErrCan,*) "P0 P0 S0 S0",i1,i2,sumI(i2,i2)
    res(i2,i2) = res(i2,i2) + sumI(i2,i2)
   End Do

   Do i2=1,n_P0
    B0m2 = 2._dp * Real( B0(p2, mP02(i1), mP02(i2)),dp )
    Do i3=1,n_S0
     sumI(i3,i3) = c_P0S0S0(i1,i2,i3)**2 * B0m2
     Do i4=i3+1,n_S0
      sumI(i3,i4) = c_P0S0S0(i1,i2,i3) * c_P0S0S0(i1,i2,i4) * B0m2
     End Do
     If (WriteOut) Write(ErrCan,*) "P0 S0 S0",i1,i2,i3,sumI(i3,i3:n_S0)
     res(i3,i3:n_S0) = res(i3,i3:n_S0) + sumI(i3,i3:n_S0)
    End Do
   End Do
  End Do
 !-----------------
 ! neutral Higgs
 !-----------------
  sumI = 0._dp
  Do i1=1,n_S0
   A0m2 = - 0.5_dp * A0( mS02(i1) )
   Do i2=1,n_S0
    sumI(i2,i2:n_S0) = c_S04(i1,i2,i2:n_S0) * A0m2
    If (WriteOut) Write(ErrCan,*) "S0 S0 S0 S0",i1,i2,sumI(i2,i2:n_S0)
    res(i2,i2:n_S0) = res(i2,i2:n_S0) + sumI(i2,i2:n_S0)

    B0m2 = 0.5_dp * Real( B0(p2,mS02(i1), mS02(i2)),dp )
    Do i3=1,n_S0
     sumI(i3,i3) = c_S03(i3,i1,i2)**2 * B0m2
     Do i4=i3+1,n_S0
      sumI(i3,i4) = c_S03(i3,i1,i2) * c_S03(i4,i1,i2) * B0m2
     End Do 
     If (WriteOut) Write(ErrCan,*) "S0 S0 S0",i1,i2,i3,sumI(i3,i3:n_S0)
     res(i3,i3:n_S0) = res(i3,i3:n_S0) + sumI(i3,i3:n_S0)
    End Do
   End Do
  End Do
 !------------------
 ! neutralinos
 !------------------
  sumI = 0._dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    G0m2 = 0.5_dp * Real( Gloop(p2,mN2(i1),mN2(i2)),dp )
    B0m2 = 2._dp * mN(i1) * mN(i2) * Real( B0(p2,mN2(i1),mN2(i2)),dp )
    Do i3=1,n_S0
     sumI(i3,i3) = (Abs(c_NNS0_L(i1,i2,i3))**2 + Abs(c_NNS0_R(i1,i2,i3))**2) &
              &   * G0m2                                                     &
              & - Real( Conjg(c_NNS0_R(i1,i2,i3)) * c_NNS0_L(i1,i2,i3),dp )   &
              &   * B0m2
     Do i4=i3+1,n_S0
      sumI(i3,i4) = ( Conjg(c_NNS0_L(i1,i2,i3)) * c_NNS0_L(i1,i2,i4)     &
               &    + Conjg(c_NNS0_R(i1,i2,i3)) * c_NNS0_R(i1,i2,i4) )   &
               &   * G0m2                                                &
               &  - ( Conjg(c_NNS0_R(i1,i2,i3)) * c_NNS0_L(i1,i2,i4)     &
               &    + Conjg(c_NNS0_L(i1,i2,i3)) * c_NNS0_R(i1,i2,i4) )   &
               &   * 0.5_dp *  B0m2
     End Do
     If (WriteOut) Write(ErrCan,*) "N N S0",i1,i2,i3,sumI(i3,i3:n_S0)
     res(i3,i3:n_S0) = res(i3,i3:n_S0) + sumI(i3,i3:n_S0)
    End Do
   End Do
  End Do

  Do i2=1,n_S0
   Do i1=2,n_S0
     res(i1,i2) = Conjg( res(i2,i1) )
   End Do
  End Do
  res = oo16pi2 * res

 End Subroutine PiScalar


 Subroutine PiScalar_NMSSM(p2, c_DDS0_L, c_UUS0_L, c_LLS0_L                   &
     & , mDSquark2, c_SdSdS0S0, mUSquark2, c_SuSuS0S0, mSlepton2, c_SlSlS0S0  &
     & , mSneutrino2, c_SnSnS0S0, c_SdSdS0, c_SuSuS0, c_SlSlS0, c_SnSnS0      &
     & , mS02, mP02, mSpm2, c_SpmS0W, c_S0WW, c_P0S0Z, c_S0ZZ, c_S0S0WW       &
     & , c_S0S0ZZ, c_SpSmS0S0, c_SpSmS0, mC, mC2, c_CCS0_L, c_CCS0_R          &
     & , c_P0P0S0S0, c_S0P0P0, c_S03, c_S04, mN, mN2, c_NNS0_L, c_NNS0_R      &
     & , mZ2, mW2, WriteOut, res)
 !-------------------------------------------------------------------------
 ! calculates the 1-loop self energies of higgs bosons
 ! The formulae of J. Bagger et al, Nucl. Phys. B491 are used.
 ! Note, that gSU3 has to be set to 0 for sleptons
 ! the input is:
 ! - p2 ......... the outer momentum squared 
 ! 
 ! the output is:
 !  res(i,j) 
 ! written 03.11.2008
 ! 10.10.01: portation to f90
 ! 11.06.03: implementing alpha_b alpha_t + alpha^2_b corrections
 !           based on A. Dedes, G. Degrassi and P. Slavich, hep-ph/0305127
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: p2, mUSquark2(6), mDSquark2(6), mSneutrino2(:)    &
     & , mS02(:), mP02(:), mSpm2(:), c_S0WW(:), c_P0S0Z(:,:), c_S0ZZ(:)     &
     & , c_S0S0WW(:), c_S0S0ZZ(:), mSlepton2(:), mC(:), mC2(:)              &
     & , c_P0P0S0S0(:,:,:), c_S0P0P0(:,:,:), c_S03(:,:,:), c_S04(:,:,:) &
     & , mN(:), mN2(:), mZ2, mW2
  Complex(dp), Intent(in) :: c_DDS0_L(:), c_UUS0_L(:), c_LLS0_L(:)           &
     & , c_SdSdS0S0(:,:,:), c_SuSuS0S0(:,:,:), c_SlSlS0S0(:,:,:)             & 
     & , c_SnSnS0S0(:,:,:)                                                   &
     & , c_SdSdS0(:,:,:), c_SuSuS0(:,:,:), c_SlSlS0(:,:,:), c_SnSnS0(:,:,:)  &
     & , c_SpmS0W(:,:), c_SpSmS0S0(:,:,:), c_SpSmS0(:,:,:), c_CCS0_L(:,:,:)&
     & , c_CCS0_R(:,:,:), c_NNS0_L(:,:,:), c_NNS0_R(:,:,:)
  Logical, Intent(in) :: WriteOut
  Complex(dp), Intent(inout) :: res(:,:)

  Integer :: i1, i2, i3, i4, i5, i_sf1, i_sf2
  Complex(dp) :: sumI(n_S0,n_S0)

  Real(dp) :: A0m2, B0m2, F0m2, G0m2

  !-----------------------
  ! initialisation
  !-----------------------
  res = 0._dp
  !-------------------
  ! fermions
  !-------------------
 Print*, "PiScalar_NMSSM"

  sumI = 0._dp
  Do i1=1,3
   sumI(1,1) = 6._dp * Abs( c_DDS0_L(i1) )**2                              &
           &   * ( (p2 - 4._dp * mf_d2(i1) )                               &
           &           * Real(B0(p2, mf_d2(i1), mf_d2(i1)),dp)             &
           &     - 2._dp * Real( A0(mf_d2(i1)),dp ) )                      &
           & +  2._dp * Abs( c_LLS0_L(i1) )**2                             &
           &   * ( (p2 - 4._dp * mf_l2(i1) )                               &
           &            * Real( B0(p2, mf_l2(i1), mf_l2(i1)),dp)           &
           &     - 2._dp * Real( A0( mf_l2(i1) ),dp ) )         
   sumI(2,2) = 6._dp * Abs( c_UUS0_L(i1) )**2                              &
           &   * ( (p2 - 4._dp * mf_u2(i1) )                               &
           &         * Real(B0(p2, mf_u2(i1), mf_u2(i1)),dp)               &
           &     - 2._dp * Real(A0(mf_u2(i1)),dp) )
   If (WriteOut) Write(ErrCan,*) "fermions",i1,sumI(1,1),sumI(2,2)
! print*, "fermions",i1,sumI(1,1),sumI(2,2)
     res = res + sumI
  End Do

  !--------------------
  ! sfermions
  !--------------------
  sumI = 0._dp
  Do i1=1,6
   A0m2 = -6._dp * A0( mDsquark2(i1) )
   Do i2=1,n_S0
    Do i3=i2,n_S0
       sumI(i2,i3) = c_SdSdS0S0(i1,i2,i3) * A0m2
    End Do
   End Do
    res = res + sumI 
If (WriteOut) Write(ErrCan,*) "Sd Sd S0 S0",i1,sumI(1,1),sumI(1,2)&
      &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
!  print*, "Sd Sd S0 S0",i1,sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
  End Do

  sumI = 0._dp
  Do i1=1,6
   A0m2 = -6._dp * A0( mUsquark2(i1) )
   Do i2=1,n_S0 
    Do i3=i2,n_S0
      sumI(i2,i3) = c_SuSuS0S0(i1,i2,i3) * A0m2
    End Do
   End Do
    res = res + sumI
If (WriteOut) Write(ErrCan,*) "Su Su S0 S0",i1,sumI(1,1),sumI(1,2)&
      &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "Su Su S0 S0",i1,sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
  End Do

  sumI = 0._dp
  Do i1=1,n_Slept
   A0m2 = -2._dp * A0( mSlepton2(i1) )
   Do i2=1,n_S0 
    Do i3=i2,n_S0
     sumI(i2,i3) = c_SlSlS0S0(i1,i2,i3) * A0m2
    End Do
   End Do
    res = res + sumI
If (WriteOut) Write(ErrCan,*) "Sl Sl S0 S0",i1,sumI(1,1),sumI(1,2)&
      &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "Sl Sl S0 S0",i1,sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
  End Do

  sumI = 0._dp
  Do i1=1,n_Sneut
   A0m2 = -2._dp * A0( mSneutrino2(i1) )
   Do i2=1,n_S0
    Do i3=i2,n_S0
     sumI(i2,i3) = c_SnSnS0S0((i1-1)*2+1,i2,i3) * A0m2
    End Do
   End Do
    res = res + sumI
 If (WriteOut) Write(ErrCan,*) "Sn Sn S0 S0",i1,sumI(1,1),sumI(1,2)&
      &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "Sn Sn S0 S0",i1,sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
  End Do

 ! .not.GenerationMixing
   sumI = 0._dp
   Do i1=1,3
    Do i2=1,2
     i_sf1 = (i1-1)*2+i2
     Do i3=1,2
      i_sf2=(i1-1)*2+i3
     B0m2 = 3._dp * Real( B0(p2, mDSquark2(i_sf1),mDSquark2(i_sf2)),dp )
      Do i4=1,n_S0
       Do i5=i4,n_S0
        sumI(i4,i5) = Conjg( c_SdSdS0(i_sf1,i_sf2,i4) ) &
                    &  * c_SdSdS0(i_sf1,i_sf2,i5) * B0m2
       End Do ! i5
      End Do ! i4
! print*, "Sd Sd S0 ",i_sf1,i_sf2,sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
      res = res + sumI
     End Do ! i3
    End Do ! i2
   End Do ! i1
   
   sumI = 0._dp
   Do i1=1,3
    Do i2=1,2
     i_sf1 = (i1-1)*2+i2
     Do i3=1,2
      i_sf2=(i1-1)*2+i3
     B0m2 = 3._dp * Real( B0(p2, mUSquark2(i_sf1), mUSquark2(i_sf2)),dp )
      Do i4=1,n_S0
       Do i5=i4,n_S0
        sumI(i4,i5) = Conjg( c_SuSuS0(i_sf1,i_sf2,i4) ) &
                    &  * c_SuSuS0(i_sf1,i_sf2,i5) * B0m2
       End Do ! i5
      End Do ! i4
      If (WriteOut) Write(ErrCan,*) "Su Su S0",i_sf1,i_sf2 &
          & ,sumI(1,1),sumI(1,2),sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "Su Su S0 ",i_sf1,i_sf2,sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
     res = res + sumI
     End Do ! i3
    End Do ! i2
   End Do ! i1

   sumI = 0._dp
   Do i1=1,3
    Do i2=1,2
     i_sf1 = (i1-1)*2+i2
     Do i3=1,2
      i_sf2 = (i1-1)*2+i3
      B0m2 = Real( B0(p2, mSlepton2(i_sf1), mSlepton2(i_sf2)),dp )
      Do i4=1,n_s0
       Do i5=i4,n_s0
        sumI(i4,i5) = Conjg( c_SlSlS0(i_sf1,i_sf2,i4) ) &
                    &  * c_SlSlS0(i_sf1,i_sf2,i5) * B0m2
       End Do ! i5
      End Do ! i4
      If (WriteOut) Write(ErrCan,*) "Sl Sl S0",i_sf1,i_sf2 &
          & ,sumI(1,1),sumI(1,2),sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "Sl Sl S0 ",i_sf1,i_sf2,sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
      res = res + sumI
     End Do
    End Do
   End Do

   sumI = 0._dp
   Do i1=1,n_sneut
    Do i2=1,2
     i_sf1 = (i1-1)*2+i2
     Do i3=1,2
      i_sf2 = (i1-1)*2+i3
    B0m2 = Real( B0(p2, mSneutrino2(i1), mSneutrino2(i1)),dp )
     Do i4=1,n_s0
      Do i5=i4,n_s0
      sumI(i4,i5) = Conjg( c_SnSnS0(i_sf1,i_sf2,i4) )&
                  & * c_SnSnS0(i_sf1,i_sf2,i5) * B0m2
      End Do ! i5
     End Do ! i4
     If (WriteOut) Write(ErrCan,*) "Sn Sn S0",i_sf1,i_sf2&
          &,sumI(1,1),sumI(1,2),sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "Sn Sn S0 ",i_sf1,i_sf2,sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
     res = res + sumI
    End Do
   End Do
  End Do

 !-----------------
 ! gauge bosons
 !-----------------
  sumI = 0._dp
  Do i2=2,n_Spm
   F0m2 = Real(Floop(p2, mSpm2(i2), mW2),dp)
   sumI(1,1) = Abs(c_SpmS0W(i2,1) )**2 * F0m2
   sumI(1,2) = Conjg(c_SpmS0W(i2,1)) * c_SpmS0W(i2,2) * F0m2
   sumI(2,2) = Abs(c_SpmS0W(i2,2))**2 * F0m2
   If (WriteOut) Write(ErrCan,*) "S+ S0 W",i2 ,sumI(1,1),sumI(1,2),sumI(2,2)
! print*, "Spm S0 W ",i2,sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
  res = res + sumI
  End Do

  sumI = 0._dp
  F0m2 = Real(Floop(p2, mSpm2(1),mW2),dp) + 7._dp * mW2 * Real(B0(p2,mW2,mW2),dp)
  sumI(1,1) = Abs(c_S0WW(1))**2 * F0m2
  sumI(1,2) = c_S0WW(1) * c_S0WW(2) * F0m2
  sumI(2,2) = Abs(c_S0WW(2))**2 * F0m2
  If (WriteOut) Write(ErrCan,*) "S0 W+ W-",sumI(1,1),sumI(1,2),sumI(2,2)
! print*, "S0 W W ",sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
  res = res + sumI

  sumI = 0._dp
  Do i2=2,n_P0
   F0m2 = Real(Floop(p2, mP02(i2), mZ2),dp)
   sumI(1,1) = Abs(c_P0S0Z(i2,1))**2 * F0m2
   sumI(1,2) = c_P0S0Z(i2,1) * c_P0S0Z(i2,2) * F0m2
   sumI(1,3) = c_P0S0Z(i2,1) * c_P0S0Z(i2,3) * F0m2
   sumI(2,2) = Abs(c_P0S0Z(i2,2))**2 * F0m2
   sumI(2,3) = c_P0S0Z(i2,2) * c_P0S0Z(i2,3) * F0m2
   If (WriteOut) Write(ErrCan,*) "P0 S0 Z",i2 ,sumI(1,1),sumI(1,2),sumI(2,2)
! print*, "P0 S0 Z ",i2, sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
   res = res + sumI
  End Do
  
  sumI = 0._dp
  F0m2 = Real(Floop(p2, mP02(1), mZ2),dp) + 7._dp * mZ2 * Real(B0(p2,mZ2,mZ2),dp)
  sumI(1,1) = Abs(c_S0ZZ(1))**2 * F0m2
  sumI(1,2) = c_S0ZZ(1) * c_S0ZZ(2) * F0m2
  sumI(2,2) = Abs(c_S0ZZ(2))**2 * F0m2
  If (WriteOut) Write(ErrCan,*) "S0 Z Z",sumI(1,1),sumI(1,2),sumI(2,2)
! print*, "S0 Z Z", sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
  res = res + sumI

  sumI = 0._dp
  Do i1=1,n_S0
   SumI(i1,i1) = c_S0S0WW(i1) * A0(mW2) +  c_S0S0ZZ(i1) * A0(mZ2)
   If (WriteOut) Write(ErrCan,*) "S0 S0 WW, ZZ",i1,sumI(i1,i1)
! print*, "S0 S0 WW, ZZ", i1, sumI(i1,i1)
   res(i1,i1) = res(i1,i1) + sumI(i1,i1)
  End Do

 !-----------------
 ! charged Higgs
 !-----------------
  sumI = 0._dp
  Do i1=1,n_Spm
   A0m2 = A0(mSpm2(i1))
   Do i3=1,n_s0
    Do i4=i3,n_s0
     sumI(i3,i4) = c_SpSmS0S0(i1,i3,i4) * A0m2
    End Do ! i4
   End Do ! i3
   If (WriteOut) Write(ErrCan,*) &
           & "Sp Sm S0 S0",i1,sumI(1,1),sumI(1,2)&
           &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "Sp Sm S0 S0", i1, sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
    res = res + sumI
  
   Do i2=1,n_Spm
    B0m2 = Real(B0(p2, mSpm2(i1), mSpm2(i2)),dp)
    Do i3=1,n_s0
     Do i4=i3,n_s0
      sumI(i3,i4) = Conjg(c_SpSmS0(i1,i2,i3)) * c_SpSmS0(i1,i2,i4) * B0m2
     End Do ! i4
    End Do ! i3
    If (WriteOut) Write(ErrCan,*) &
            & "Sp Sm S0",i1,i2,sumI(1,1),sumI(1,2),sumI(1,3)&
            &,sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "Sp Sm S0", i1, i2, sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
    res = res + sumI
   End Do
  End Do

 !------------------------------
 ! Charginos
 !------------------------------
  sumI = 0._dp
  Do i1=1,n_char
   G0m2 = Real(Gloop(p2, mC2(i1), mC2(i1)),dp)
   B0m2 = 4._dp * mC2(i1) * Real(B0(p2, mC2(i1), mC2(i1)), dp)
   Do i3=1,n_s0
    sumI(i3,i3) =                                                           &
       &  (Abs(c_CCS0_L(i1,i1,i3))**2 + Abs(c_CCS0_R(i1,i1,i3))**2) * G0m2  &
       &  - Real(Conjg(c_CCS0_R(i1,i1,i3)) * c_CCS0_L(i1,i1,i3),dp) * B0m2
    Do i4=i3+1,n_S0
     sumI(i3,i4) = (Conjg(c_CCS0_L(i1,i1,i3))*c_CCS0_L(i1,i1,i4)           &
            &   + Conjg(c_CCS0_R(i1,i1,i3))* c_CCS0_R(i1,i1,i4)) * G0m2    &
            & - 0.5_dp *(Conjg(c_CCS0_R(i1,i1,i3)) * c_CCS0_L(i1,i1,i4)    &
            &           + Conjg(c_CCS0_L(i1,i1,i3)) * c_CCS0_R(i1,i1,i4))*B0m2
    End Do
   End Do
   If (WriteOut) Write(ErrCan,*) "C C S0",i1,i1,sumI(1,1),sumI(1,2),sumI(1,3) &
       &,sumI(2,2),sumI(2,3),sumI(3,3) 
! print*, "C C S0", i1,i1, sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
    res = res + sumI

   sumI = 0._dp
   Do i2=i1+1,n_char
    G0m2 = 2._dp * Real(Gloop(p2, mC2(i1), mC2(i2)),dp)
    B0m2 = 8._dp * mC(i1) * mC(i2) * Real(B0(p2, mC2(i1), mC2(i2)),dp)

    Do i3=1,n_s0
     sumI(i3,i3) =                                                           &
          &  (Abs(c_CCS0_L(i1,i2,i3))**2 +Abs(c_CCS0_R(i1,i2,i3))**2) * G0m2 &
          & - Real(Conjg(c_CCS0_R(i1,i2,i3)) * c_CCS0_L(i1,i2,i3), dp) * B0m2
    Do i4=i3+1,n_S0
     sumI(i3,i4) = (Conjg(c_CCS0_L(i1,i2,i3))*c_CCS0_L(i1,i2,i4)          &
            &   + Conjg(c_CCS0_R(i1,i2,i3))* c_CCS0_R(i1,i2,i4)) * G0m2   &
            & - 0.5_dp * (Conjg(c_CCS0_R(i1,i2,i3)) * c_CCS0_L(i1,i2,i4)  &
            &         + Conjg(c_CCS0_L(i1,i2,i3)) * c_CCS0_R(i1,i2,i4)) * B0m2
     End Do ! i4
    End Do ! i3
    If (WriteOut) Write(ErrCan,*) "C C S0",i1,i2,sumI(1,1),sumI(1,2),sumI(1,3)&
       &,sumI(2,2),sumI(2,3),sumI(3,3) 
! print*, "C C S0", i1,i2, sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
    res = res + sumI
   End Do
  End Do

 !---------------------
 ! pseudoscalar Higgs
 !---------------------
  sumI = 0._dp
  Do i1=1,n_P0
   A0m2 = - A0( mP02(i1) )
   Do i2=1,n_S0
    Do i3=i2,n_S0
    sumI(i2,i3) = c_P0P0S0S0(i1,i2,i3) * A0m2
    End Do
   End Do
   res = res + sumI
 If (WriteOut) Write(ErrCan,*) "P0 P0 S0 S0",i1,sumI(1,1),sumI(1,2),sumI(1,3)&
       &,sumI(2,2),sumI(2,3),sumI(3,3) 
! print*, "P0 P0 S0 S0",i1, sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
  End Do

  sumI = 0._dp
   Do i1=1,n_P0
    Do i2=1,n_P0
     B0m2 = 2._dp * Real( B0(p2, mP02(i1), mP02(i2)),dp )
     Do i3=1,n_S0
      sumI(i3,i3) = c_S0P0P0(i1,i2,i3)**2 * B0m2
      Do i4=i3+1,n_S0
       sumI(i3,i4) = c_S0P0P0(i1,i2,i3) * c_S0P0P0(i1,i2,i4) * B0m2
      End Do ! i4
    End Do ! i3
    res = res + sumI
 If (WriteOut) Write(ErrCan,*) "S0 P0 P0",i1,i2,sumI(1,1),sumI(1,2),sumI(1,3)&
       &,sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "S0 P0 P0",i1, i2, sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
   End Do 
  End Do

 !-----------------
 ! neutral Higgs
 !-----------------
  sumI = 0._dp
  Do i1=1,n_S0
   A0m2 = - 0.5_dp * A0( mS02(i1) )
   Do i2=1,n_S0
    Do i3=i2,n_S0
    sumI(i2,i3) = c_S04(i2,i3,i1) * A0m2
    End Do ! i3
   End Do ! i2
  res = res + sumI
If (WriteOut) Write(ErrCan,*) "S0 S0 S0 S0",i1,i2,sumI(1,1),sumI(1,2)&
           &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "S0 S0 S0 S0",i1, sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
  End Do


  sumI = 0._dp
  Do i1=1,n_S0
   Do i2=1,n_S0
    B0m2 = 0.5_dp * Real( B0(p2,mS02(i1), mS02(i2)),dp )
    Do i3=1,n_S0
     sumI(i3,i3) = c_S03(i3,i1,i2)**2 * B0m2
     Do i4=i3+1,n_S0
      sumI(i3,i4) = c_S03(i3,i1,i2) * c_S03(i4,i1,i2) * B0m2
     End Do ! i4 
    End Do ! i3
    res = res + sumI
If (WriteOut) Write(ErrCan,*) "S0 S0 S0",i1,i2,sumI(1,1),sumI(1,2)&
           &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "S0 S0 S0",i1,i2, sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
   End Do
  End Do

 !------------------
 ! neutralinos
 !------------------
  sumI = 0._dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    G0m2 = 0.5_dp * Real( Gloop(p2,mN2(i1),mN2(i2)),dp )
    B0m2 = 2._dp * mN(i1) * mN(i2) * Real( B0(p2,mN2(i1),mN2(i2)),dp )
    Do i3=1,n_S0
     sumI(i3,i3) = (Abs(c_NNS0_L(i1,i2,i3))**2 + Abs(c_NNS0_R(i1,i2,i3))**2) &
              &   * G0m2                                                     &
              & - Real( Conjg(c_NNS0_R(i1,i2,i3)) * c_NNS0_L(i1,i2,i3),dp )   &
              &   * B0m2
     Do i4=i3+1,n_S0
      sumI(i3,i4) = ( Conjg(c_NNS0_L(i1,i2,i3)) * c_NNS0_L(i1,i2,i4)     &
               &    + Conjg(c_NNS0_R(i1,i2,i3)) * c_NNS0_R(i1,i2,i4) )   &
               &   * G0m2                                                &
               &  - ( Conjg(c_NNS0_R(i1,i2,i3)) * c_NNS0_L(i1,i2,i4)     &
               &    + Conjg(c_NNS0_L(i1,i2,i3)) * c_NNS0_R(i1,i2,i4) )   &
               &   * 0.5_dp *  B0m2
     End Do !i4
    res(i3,i3:n_S0) = res(i3,i3:n_S0) + sumI(i3,i3:n_S0)
    End Do ! i3
If (WriteOut) Write(ErrCan,*) "N N S0",i1,i2,sumI(1,1),sumI(1,2)&
           &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3)
! print*, "N N S0",i1, i2, sumI(1,1),sumI(1,2)&
!        &,sumI(1,3),sumI(2,2),sumI(2,3),sumI(3,3) 
   End Do ! i2
  End Do ! i1 

  Do i2=1,n_S0
   Do i1=2,n_S0
     res(i1,i2) = Conjg( res(i2,i1) )
   End Do
  End Do
  res = oo16pi2 * res

 End Subroutine PiScalar_NMSSM

 Subroutine PiSlepton(p2, mSlepton2, c_Sq4e, c_Sq4Z, c_Sq2Z, c_Sq2W, mSneut2  &
        & , c_Sq4W, mN, mN2, c_FNSf_L, c_FNSf_R, mC2, c_CFpSf_L, c_CFpSf_R    &
        & , mP02, c_P0SfSf, c_P0P0SfSf, mS02, c_S0SfSf, c_S0S0SfSf            &
        & , mSpm2, c_SpmSfSfp, c_SpmSpmSfSfp, mSup2, c_SuSf4                  &
        & , mSdown2, c_SdSf4, c_SlSf4, c_SnSf4, mZ2, mW2, WriteOut, res)
 !-------------------------------------------------------------------------
 ! calculates the 1-loop self energies of sleptons
 ! The formulae of J. Bagger et al, Nucl. Phys. B491 are used.
 ! the input is:
 ! - p2 ......... the outer momentum squared 
 ! 
 ! the output is:
 !  res(i,j)
 ! written by Werner Porod, 31.12.01
 ! 05.02.02: adding the 4-vertices Z-Z-Sl-Sl and  W-W-Sl-Sl
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: p2, mN(:), mN2(:), mC2(:), mP02(:), mS02(:)        &
     & , mSpm2(:), mSup2(6), mSdown2(6), mSlepton2(:), mSneut2(:), c_Sq2Z(:) &
     & , c_Sq2W(:), mZ2, mW2
  Complex(dp), Intent(in) :: c_FNSf_L(:,:,:), c_FNSf_R(:,:,:)                 &
     & , c_CFpSf_L(:,:,:), c_CFpSf_R(:,:,:), c_Sq4e(:,:,:), c_Sq4Z(:,:,:)     &
     & , c_Sq4W(:,:,:), c_P0P0SfSf(:,:,:), c_P0SfSf(:,:,:), c_S0S0SfSf(:,:,:) &
     & , c_S0SfSf(:,:,:), c_SpmSfSfp(:,:,:) , c_SpmSpmSfSfp(:,:,:)            &
     & , c_SuSf4(:,:,:), c_SdSf4(:,:,:), c_SlSf4(:,:,:), c_SnSf4(:,:,:)
  Logical, Intent(in) :: WriteOut
  Complex(dp), Intent(out) :: res(:,:)

  Integer :: i1, i2, i3, i2_min, i2_max, i3_min, i3_max, i4, i4_min, i4_max &
     & , ind_1, ind_2
  Complex(dp) :: sumI(6,6)
  Real(dp) :: A0m2, F0m2, B0m2, G0m2

 !-----------------------
 ! initialisation
 !-----------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "PiSlepton"

  res = ZeroC
 !--------------------------------
 ! electroweak gauge bosons
 !--------------------------------
  A0m2 = A0(mZ2)
  sumI = ZeroC
  Do i1=1,n_slept
   sumI(i1,i1) = c_Sq2Z(i1) * A0m2
  End Do
  If (WriteOut) Write(ErrCan,*) "Sl Sl Z Z",(sumI(i1,i1),i1=1,n_slept)
  res = res + sumI

  Do i1=1,n_slept
   F0m2 = Floop(p2, mSlepton2(i1), 0._dp )
   If (GenerationMixing) Then
    i2_min = 1
    i2_max = n_slept
   Else
    If (i1.Le.2) Then
     i2_min = 1
    Elseif (i1.Le.4) Then
     i2_min = 3
    Else
     i2_min = 5
    End If
    i2_max = i2_min + 1
   End If
   i3_min = i2_min
   i3_max = i2_max
   sumI = ZeroC
   Do i2=i2_min,i2_max
    Do i3=i3_min,i3_max
     sumI(i2,i3) = c_Sq4e(i1,i2,i3) * F0m2
    End Do
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "slepton photon",i1
    Do i2=i2_min,i2_max
     Write(ErrCan,*) "SumI(i2,:),",i2,sumI(i2,i3_min:i3_max)
    End Do
   End If
   res = res + sumI

   F0m2 = Floop(p2, mSlepton2(i1), mZ2 )
   sumI = ZeroC
   Do i2=i2_min,i2_max
    Do i3=i3_min,i3_max
     sumI(i2,i3) = c_Sq4Z(i1,i2,i3) * F0m2
    End Do
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "slepton Z",i1
    Do i2=i2_min,i2_max
     Write(ErrCan,*) "SumI(i2,:),",i2,sumI(i2,i3_min:i3_max)
    End Do
   End If
   res = res + sumI
  End Do

  A0m2 = A0(mW2)
  sumI = ZeroC
  Do i1=1,n_slept
   sumI(i1,i1) = c_Sq2W(i1) * A0m2
  End Do
  If (WriteOut) Write(ErrCan,*) "Sl Sl W W",(sumI(i1,i1),i1=1,n_slept)
  res = res + sumI

  sumI = ZeroC
  Do i1=1,n_sneut
   F0m2 = Floop(p2, mSneut2(i1), mW2 )
   If (GenerationMixing) Then
    i2_min = 1
    i2_max = n_slept
   Else
    i2_max = 2*i1
    i2_min = i2_max - 1
   End If
   i3_min = i2_min
   i3_max = i2_max
   sumI = ZeroC
   Do i2=i2_min,i2_max
    Do i3=i3_min,i3_max
     sumI(i2,i3) = c_Sq4W(i1,i2,i3) * F0m2
    End Do
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "sneutrino W",i1
    Do i2=i2_min,i2_max
     Write(ErrCan,*) "SumI(i2,:),",i2,sumI(i2,i3_min:i3_max)
    End Do
   End If
   res = res + sumI
  End Do

  !-----------------------
  ! neutralinos
  !-----------------------
  Do i1=1,n_neut
   Do i2=1,n_sneut ! number of effective slepton generations
    G0m2 = Gloop(p2, mN2(i1), mf_l2(i2) )
    B0m2 = - 2._dp * mN(i1) * mf_l(i2) * B0(p2, mN2(i1), mf_l2(i2) )
    If (GenerationMixing) Then
     i3_min = 1
     i3_max = n_slept
    Else
     i3_min = 2*i2-1
     i3_max = 2*i2
    End If
    i4_min = i3_min
    i4_max = i3_max
    SumI = ZeroC
    Do i3=i3_min,i3_max
     Do i4=i4_min,i4_max
      SumI(i3,i4) = ( Conjg(c_FNSf_L(i2,i1,i3))* c_FNSf_L(i2,i1,i4)           &
                &   + Conjg(c_FNSf_R(i2,i1,i3))* c_FNSf_R(i2,i1,i4) ) * G0m2  &
                & + ( c_FNSf_L(i2,i1,i4) * Conjg(c_FNSf_R(i2,i1,i3))          &
                &   + c_FNSf_R(i2,i1,i4) * Conjg(c_FNSf_L(i2,i1,i3)) ) * B0m2
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "chi0 f",i1,i2
     Do i3=i3_min,i3_max
      Write(ErrCan,*) "SumI(i3,:)",SumI(i3,i4_min:i4_max)
     End Do
    End If
    res = res + SumI
   End Do
  End Do

 !-----------------------
 ! charginos
 !-----------------------
  Do i1=1,n_char
   G0m2 = Gloop(p2, mC2(i1), 0._dp )
   Do i2=1,n_sneut ! number of effective slepton generations
    If (GenerationMixing) Then
     i3_min = 1
     i3_max = n_slept
    Else
     i3_min = 2*i2-1
     i3_max = 2*i2
    End If
    i4_min = i3_min
    i4_max = i3_max
    SumI = ZeroC
    Do i3=i3_min,i3_max
     Do i4=i4_min,i4_max
      SumI(i3,i4) = ( Conjg(c_CFpSf_L(i1,i2,i3))* c_CFpSf_L(i1,i2,i4)         &
                &  + Conjg(c_CFpSf_R(i1,i2,i3))* c_CFpSf_R(i1,i2,i4) ) * G0m2
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "chi+ f'",i1,i2
     Do i3=i3_min,i3_max
      Write(ErrCan,*) "SumI(i3,:)",SumI(i3,i4_min:i4_max)
     End Do
    End If
    res = res + SumI
   End Do
  End Do

  !-------------------------------
  ! sfermion quartic interactions
  !-------------------------------
  ! up squarks
  !------------
  Do i1=1,6
   A0m2 = A0( mSup2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_slept
     Do i3=1,n_slept
      SumI(i2,i3) = - c_SuSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Su Su Sf Sf",i1
     Do i2=1,6
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    res = res + sumI
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,n_sneut ! number of effective slepton generations
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_SuSf4(i1,ind_1,ind_2) * A0m2
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Su Su Sf Sf",i1
     Do i2=1,n_slept
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
   End If !
  End Do

  !--------------
  ! down squarks
  !--------------
  Do i1=1,6
   A0m2 = A0( mSdown2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_slept
     Do i3=1,n_slept
      SumI(i2,i3) = - c_SdSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sd Sd Sf Sf",i1
     Do i2=1,n_slept
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    res = res + sumI
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,n_sneut ! number of effective slepton generations
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_SdSf4(i1,ind_1,ind_2) * A0m2
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sd Sd Sf Sf",i1
     Do i2=1,n_slept
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
   End If !
  End Do

  !------------------
  ! charged sleptons
  !------------------
  Do i1=1,n_slept
   A0m2 = A0( mSlepton2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_slept
     Do i3=1,n_slept
      SumI(i2,i3) = - c_SlSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sl Sl Sf Sf",i1
     Do i2=1,n_slept
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    res = res + sumI
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,n_sneut ! number of effective slepton generations
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_SlSf4(i1,ind_1,ind_2) * A0m2
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sl Sl Sf Sf",i1
     Do i2=1,n_slept
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
   End If !
  End Do

  !----------------
  ! sneutrinos
  !----------------
  Do i1=1,n_sneut
   A0m2 = A0( mSneut2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_slept
     Do i3=1,n_slept
      SumI(i2,i3) = - c_SnSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sn Sn Sf Sf",i1
     Do i2=1,n_slept
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    res = res + sumI
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,n_sneut ! number of effective slepton generations
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_SnSf4(i1,ind_1,ind_2) * A0m2
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sn Sn Sf Sf",i1
     Do i2=1,n_slept
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
   End If !
  End Do

  !---------------------
  ! pseudoscalar Higgs
  !---------------------
  Do i1=1,n_P0
   A0m2 = A0( mP02(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_slept
     Do i3=1,n_slept
      SumI(i2,i3) = - c_P0P0SfSf(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "P0 P0 Sf Sf",i1
     Do i2=1,n_slept
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    res = res + sumI
 
    Do i2=1,n_slept
     SumI = ZeroC
     B0m2 = B0(p2, mP02(i1), mSlepton2(i2) )
     Do i3=1,n_slept
      Do i4=1,n_slept
       SumI(i3,i4) = Conjg(c_P0SfSf(i1,i3,i2)) * c_P0SfSf(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "P0 Sf Sf",i1,i2
      Do i3=1,n_slept
       Write(ErrCan,*) "SumI(i1,:)",i3,sumI(i3,:)
      End Do
     End If
     res = res + sumI
    End Do

   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,n_sneut ! number of effective slepton generations
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_P0P0SfSf(i1,ind_1,ind_2) * A0m2 
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "P0 P0 Sf Sf",i1
     Do i2=1,n_slept
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
 
    Do i2=1,n_slept
     SumI = ZeroC
     B0m2 = B0(p2, mP02(i1), mSlepton2(i2) )
     If (i2.Le.2) Then
      i3_min=1
     Else If (i2.Le.4) Then
      i3_min=3
     Else 
      i3_min=5
     End If
     i3_max = i3_min + 1
     i4_min = i3_min
     i4_max = i3_max
     Do i3=i3_min, i3_max
      Do i4=i4_min, i4_max
       SumI(i3,i4) = Conjg(c_P0SfSf(i1,i3,i2)) * c_P0SfSf(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "P0 Sf Sf",i1,i2
      Do i3=i3_min, i3_max
       Write(ErrCan,*) "SumI(i1,:)",sumI(i3,i4_min:i4_max)
      End Do
     End If
     res = res + sumI
    End Do
   End If !
  End Do

  !---------------------
  ! neutral Higgs
  !---------------------
  Do i1=1,n_S0
   A0m2 = A0( mS02(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_slept
     Do i3=1,n_slept
      SumI(i2,i3) = - c_S0S0SfSf(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "S0 S0 Sf Sf",i1
     Do i2=1,n_slept
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    res = res + sumI
 
    Do i2=1,n_slept
     SumI = ZeroC
     B0m2 = B0(p2, mS02(i1), mSlepton2(i2) )
     Do i3=1,n_slept
      Do i4=1,n_slept
       SumI(i3,i4) = Conjg(c_S0SfSf(i1,i3,i2)) * c_S0SfSf(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "S0 Sf Sf",i1,i2
      Do i3=1,n_slept
       Write(ErrCan,*) "SumI(i1,:)",i3,sumI(i3,:)
      End Do
     End If
     res = res + sumI
    End Do

   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,n_sneut ! number of effective slepton generations
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_S0S0SfSf(i1,ind_1,ind_2) * A0m2 
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "S0 S0 Sf Sf",i1
     Do i2=1,n_slept
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
 
    Do i2=1,n_slept
     SumI = ZeroC
     B0m2 = B0(p2, mS02(i1), mSlepton2(i2) )
     If (i2.Le.2) Then
      i3_min=1
     Else If (i2.Le.4) Then
      i3_min=3
     Else 
      i3_min=5
     End If
     i3_max = i3_min + 1
     i4_min = i3_min
     i4_max = i3_max
     Do i3=i3_min, i3_max
      Do i4=i4_min, i4_max
       SumI(i3,i4) = Conjg(c_S0SfSf(i1,i3,i2)) * c_S0SfSf(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "S0 Sf Sf",i1,i2
      Do i3=i3_min, i3_max
       Write(ErrCan,*) "SumI(i1,:)",sumI(i3,i4_min:i4_max)
      End Do
     End If
     res = res + sumI
    End Do
   End If !
  End Do

  !---------------------
  ! charged Higgs
  !---------------------
  Do i1=1,n_Spm
   A0m2 = A0( mSpm2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_slept
     Do i3=1,n_slept
      SumI(i2,i3) = - c_SpmSpmSfSfp(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Spm Spm Sf Sf",i1
     Do i2=1,n_slept
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    res = res + sumI
 
    Do i2=1,n_sneut
     SumI = ZeroC
     B0m2 = B0(p2, mSpm2(i1), mSneut2(i2) )
     Do i3=1,n_slept
      Do i4=1,n_slept
       SumI(i3,i4) = Conjg(c_SpmSfSfp(i1,i3,i2)) * c_SpmSfSfp(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "Spm Sf Sfp",i1,i2
      Do i3=1,n_slept
       Write(ErrCan,*) "SumI(i1,:)",i3,sumI(i3,:)
      End Do
     End If
     res = res + sumI
    End Do

   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,3
     Do i3=1,2
      Do i4=1,2
      ind_1 = 2*(i2-1)+i3
      ind_2 = 2*(i2-1)+i4
      SumI(ind_1,ind_2) = - c_SpmSpmSfSfp(i1,ind_1,ind_2) * A0m2 
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Spm Spm Sf Sf",i1
     Do i2=1,6
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
 
    Do i2=1,n_sneut
     SumI = ZeroC
     B0m2 = B0(p2, mSpm2(i1), mSneut2(i2) )
     i3_min=2*i2-1
     i3_max = i3_min + 1
     i4_min = i3_min
     i4_max = i3_max
     Do i3=i3_min, i3_max
      Do i4=i4_min, i4_max
       SumI(i3,i4) = Conjg(c_SpmSfSfp(i1,i3,i2)) * c_SpmSfSfp(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "Spm Sf Sfp",i1,i2
      Do i3=i3_min, i3_max
       Write(ErrCan,*) "SumI(i1,:)",sumI(i3,i4_min:i4_max)
      End Do
     End If
     res = res + sumI
    End Do
   End If !
  End Do

  res = oo16pi2 * res  

  Iname = Iname - 1
   
 End Subroutine PiSlepton 


 Subroutine PiSneutrino(p2, mSneut2, c_Sq4Z, mSlepton2, c_Sq4W, c_Sq2Z        &
   & , c_Sq2W, mN2, c_FNSf_R, mC, mC2, c_CFpSf_L, c_CFpSf_R, mP02, c_P0P0SfSf &
   & , mS02, c_S0SfSf, c_S0S0SfSf, mSpm2, c_SpmSfSfp, c_SpmSpmSfSfp           &
   & , mSup2, c_SuSf4, mSdown2, c_SdSf4, c_SlSf4, c_SnSf4, mZ2, mW2           &
   & , WriteOut, res)
 !-------------------------------------------------------------------------
 ! calculates the 1-loop self energies of sneutrinos
 ! The formulae of J. Bagger et al, Nucl. Phys. B491 are used.
 ! the input is:
 ! - p2 ......... the outer momentum squared 
 ! 
 ! the output is:
 !  res(i,j)
 ! written by Werner Porod, 05.01.02
 ! 05.02.02: adding contributions of 4-vertices Sn-Sn-Z-Z and  Sn-Sn-W-W
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: p2, mN2(:), mC(:), mC2(:), mP02(:), mS02(:)        &
     & , mSpm2(:), mSup2(6), mSdown2(6), mSlepton2(:), mSneut2(:), c_Sq2Z(:) &
     & , c_Sq2W(:), mZ2, mW2
  Complex(dp), Intent(in) :: c_FNSf_R(:,:,:), c_CFpSf_L(:,:,:)                &
     & , c_CFpSf_R(:,:,:), c_Sq4Z(:,:,:), c_Sq4W(:,:,:), c_P0P0SfSf(:,:,:)    &
     & , c_S0S0SfSf(:,:,:), c_S0SfSf(:,:,:), c_SpmSfSfp(:,:,:)                &
     & , c_SpmSpmSfSfp(:,:,:), c_SuSf4(:,:,:), c_SdSf4(:,:,:), c_SlSf4(:,:,:) &
     & , c_SnSf4(:,:,:)
  Logical, Intent(in) :: WriteOut
  Complex(dp), Intent(out) :: res(:,:)

  Integer :: i1, i2, i3, i4, i2_min, i2_max
  Complex(dp) :: sumI(3,3)
  Real(dp) :: A0m2, F0m2, B0m2, G0m2

 !-----------------------
 ! initialisation
 !-----------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "PiSneutrino"

  res = ZeroC
 !--------------------------------
 ! electroweak gauge bosons
 !--------------------------------
  A0m2 = A0(mZ2)
  sumI = ZeroC
  Do i1=1,n_sneut
   sumI(i1,i1) = c_Sq2Z(i1) * A0m2
  End Do
  If (WriteOut) Write(ErrCan,*) "Sn Sn Z Z",(sumI(i1,i1),i1=1,n_sneut)
  res = res + sumI

  Do i1=1,n_Sneut
   If (GenerationMixing) Then
    i2_min = 1
    i2_max = n_Sneut
   Else
    i2_min = i1
    i2_max = i1
   End If
   sumI = ZeroC

   F0m2 = Floop(p2, mSneut2(i1), mZ2 )
   sumI = ZeroC
   Do i2=i2_min,i2_max
    Do i3=i2_min,i2_max
     sumI(i2,i3) = c_Sq4Z(i1,i2,i3) * F0m2
    End Do
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "Sn Z",i1
    Do i2=i2_min,i2_max
     Write(ErrCan,*) "SumI(i2,:),",i2,sumI(i2,i2_min:i2_max)
    End Do
   End If
   res = res + sumI
  End Do

  A0m2 = A0(mW2)
  sumI = ZeroC
  Do i1=1,n_sneut
   sumI(i1,i1) = c_Sq2W(i1) * A0m2
  End Do
  If (WriteOut) Write(ErrCan,*) "Sn Sn W W",(sumI(i1,i1),i1=1,n_sneut)
  res = res + sumI

  Do i1=1,n_slept
   F0m2 = Floop(p2, mSlepton2(i1), mW2 )
   If (GenerationMixing) Then
    i2_min = 1
    i2_max = n_sneut
   Else
    i2_min = (i1+1)/2
    i2_max = i2_min
   End If
   sumI = ZeroC
   Do i2=i2_min,i2_max
    Do i3=i2_min,i2_max
     sumI(i2,i3) = c_Sq4W(i1,i2,i3) * F0m2
    End Do
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "slepton W",i1
    Do i2=i2_min,i2_max
     Write(ErrCan,*) "SumI(i2,:),",i2,sumI(i2,i2_min:i2_max)
    End Do
   End If
   res = res + sumI
  End Do

  !-----------------------
  ! neutralinos
  !-----------------------
  Do i1=1,n_neut
   G0m2 = Gloop(p2, mN2(i1), 0._dp )
   Do i2=1,n_sneut
    If (GenerationMixing) Then
     i2_min = 1
     i2_max = n_sneut
    Else
     i2_min = i2
     i2_max = i2
    End If
    SumI = ZeroC
    Do i3=i2_min,i2_max
     Do i4=i2_min,i2_max
      SumI(i3,i4) = Conjg(c_FNSf_R(i2,i1,i3))* c_FNSf_R(i2,i1,i4) * G0m2
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "chi0 f",i1,i2
     Do i3=i2_min,i2_max
      Write(ErrCan,*) "SumI(i3,:)",SumI(i3,i2_min:i2_max)
     End Do
    End If
    res = res + SumI
   End Do
  End Do
 !-----------------------
 ! charginos
 !-----------------------
  Do i1=1,n_char
   Do i2=1,n_sneut
    G0m2 = Gloop(p2, mC2(i1), mf_l2(i2) )
    B0m2 = - 2._dp * mC(i1) * mf_l(i2) * B0(p2, mC2(i1), mf_l2(i2) )
    If (GenerationMixing) Then
     i2_min = 1
     i2_max = n_sneut
    Else
     i2_min = i2
     i2_max = i2
    End If
    SumI = ZeroC
    Do i3=i2_min,i2_max
     Do i4=i2_min,i2_max
      SumI(i3,i4) = ( Conjg(c_CFpSf_L(i1,i2,i3))* c_CFpSf_L(i1,i2,i4)         &
                &  + Conjg(c_CFpSf_R(i1,i2,i3))* c_CFpSf_R(i1,i2,i4) ) * G0m2 &
                & + ( c_CFpSf_L(i1,i2,i4) * Conjg(c_CFpSf_R(i1,i2,i3))        &
                &   + c_CFpSf_R(i1,i2,i4) * Conjg(c_CFpSf_L(i1,i2,i3)) ) * B0m2
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "chi+ f'",i1,i2
     Do i3=i2_min,i2_max
      Write(ErrCan,*) "SumI(i3,:)",SumI(i3,i2_min:i2_max)
     End Do
    End If
    res = res + SumI
   End Do
  End Do
  !-------------------------------
  ! sfermion quartic interactions
  !-------------------------------
  ! up squarks
  !------------
  Do i1=1,6
   A0m2 = A0( mSup2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_sneut
     Do i3=1,n_sneut
      SumI(i2,i3) = - c_SuSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Su Su Sf Sf",i1
     Do i2=1,n_sneut
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:n_sneut)
     End Do
    End If
    res = res + sumI
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,3
     SumI(i2,i2) = - c_SuSf4(i1,i2,i2) * A0m2
     If (WriteOut) Write(ErrCan,*) "Su Su Sf Sf",i1,i2,sumI(i2,i2)
    End Do
    res = res + sumI
   End If !
  End Do
  !--------------
  ! down squarks
  !--------------
  Do i1=1,6
   A0m2 = A0( mSdown2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_sneut
     Do i3=1,n_sneut
      SumI(i2,i3) = - c_SdSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sd Sd Sf Sf",i1
     Do i2=1,n_sneut
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:n_sneut)
     End Do
    End If
    res = res + sumI
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,n_sneut
     SumI(i2,i2) = - c_SdSf4(i1,i2,i2) * A0m2
     If (WriteOut) Write(ErrCan,*) "Sd Sd Sf Sf",i1,i2,sumI(i2,i2)
    End Do
    res = res + sumI
   End If !
  End Do
  !------------------
  ! charged sleptons
  !------------------
  Do i1=1,n_slept
   A0m2 = A0( mSlepton2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_sneut
     Do i3=1,n_sneut
      SumI(i2,i3) = - c_SlSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sl Sl Sf Sf",i1
     Do i2=1,n_sneut
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:n_sneut)
     End Do
    End If
    res = res + sumI
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,n_sneut
     SumI(i2,i2) = - c_SlSf4(i1,i2,i2) * A0m2
     If (WriteOut) Write(ErrCan,*) "Sl Sl Sf Sf",i1,i2,sumI(i2,i2)
    End Do
    res = res + sumI
   End If !
  End Do
  !----------------
  ! sneutrinos
  !----------------
  Do i1=1,n_sneut
   A0m2 = A0( mSneut2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_sneut
     Do i3=1,n_sneut
      SumI(i2,i3) = - c_SnSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sn Sn Sf Sf",i1
     Do i2=1,n_sneut
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:n_sneut)
     End Do
    End If
    res = res + sumI
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,n_sneut
     SumI(i2,i2) = - c_SnSf4(i2,i1,i1) * A0m2
    End Do
    If (WriteOut) Write(ErrCan,*) "Sn Sn Sf Sf" &
              & ,i1,sumI(1,1),sumi(2,2),sumi(3,3)
    res = res + sumI
   End If !
  End Do
  !---------------------
  ! pseudoscalar Higgs
  !---------------------
  Do i1=1,n_P0
   A0m2 = A0( mP02(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_sneut
     Do i3=1,n_sneut
      SumI(i2,i3) = - c_P0P0SfSf(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "P0 P0 Sf Sf",i1
     Do i2=1,n_sneut
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:n_sneut)
     End Do
    End If
    res = res + sumI
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,n_sneut
     SumI(i2,i2) = - c_P0P0SfSf(i1,i2,i2) * A0m2 
     If (WriteOut) Write(ErrCan,*) "P0 P0 Sf Sf",i1,i2,sumI(i2,i2)
    End Do
    res = res + sumI
 
   End If !
  End Do
  !---------------------
  ! neutral Higgs
  !---------------------
  Do i1=1,n_S0
   A0m2 = A0( mS02(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_sneut
     Do i3=1,n_sneut
      SumI(i2,i3) = - c_S0S0SfSf(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "S0 S0 Sf Sf",i1
     Do i2=1,n_sneut
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:n_sneut)
     End Do
    End If
    res = res + sumI
 
    Do i2=1,n_sneut
     SumI = ZeroC
     B0m2 = B0(p2, mS02(i1), mSneut2(i2) )
     Do i3=1,n_sneut
      Do i4=1,n_sneut
       SumI(i3,i4) = Conjg(c_S0SfSf(i1,i3,i2)) * c_S0SfSf(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "S0 Sf Sf",i1,i2
      Do i3=1,n_sneut
       Write(ErrCan,*) "SumI(i1,:)",i3,sumI(i3,1:n_sneut)
      End Do
     End If
     res = res + sumI
    End Do

   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,3
     SumI(i2,i2) = - c_S0S0SfSf(i1,i2,i2) * A0m2 
     If (WriteOut) Write(ErrCan,*) "S0 S0 Sf Sf",i1,i2,sumI(i2,i2)
    End Do
    res = res + sumI
 
    Do i2=1,n_sneut
     SumI = ZeroC
     B0m2 = B0(p2, mS02(i1), mSneut2(i2) )
     SumI(i2,i2) = Conjg(c_S0SfSf(i1,i2,i2)) * c_S0SfSf(i1,i2,i2) * B0m2 
     If (WriteOut) Write(ErrCan,*) "S0 Sf Sf",i1,i2,sumI(i2,i2)
     res = res + sumI
    End Do
   End If !
  End Do
  !---------------------
  ! charged Higgs
  !---------------------
  Do i1=1,n_Spm
   A0m2 = A0( mSpm2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,n_sneut
     Do i3=1,n_sneut
      SumI(i2,i3) = - c_SpmSpmSfSfp(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Spm Spm Sf Sf",i1
     Do i2=1,n_sneut
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:n_sneut)
     End Do
    End If
    res = res + sumI
 
    Do i2=1,n_slept
     SumI = ZeroC
     B0m2 = B0(p2, mSpm2(i1), mSlepton2(i2) )
     Do i3=1,n_sneut
      Do i4=1,n_sneut
       SumI(i3,i4) = Conjg(c_SpmSfSfp(i1,i3,i2)) * c_SpmSfSfp(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "Spm Sf Sfp",i1,i2
      Do i3=1,n_sneut
       Write(ErrCan,*) "SumI(i1,:)",i3,sumI(i3,1:n_sneut)
      End Do
     End If
     res = res + sumI
    End Do

   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,3
     SumI(i2,i2) = - c_SpmSpmSfSfp(i1,i2,i2) * A0m2 
     If (WriteOut) Write(ErrCan,*) "Spm Spm Sf Sf",i1,i2,sumI(i2,i2)
    End Do
    res = res + sumI
 
    Do i2=1,6
     SumI = ZeroC
     B0m2 = B0(p2, mSpm2(i1), mSlepton2(i2) )
     i3 = (i2+1)/2
     SumI(i3,i3) = Abs(c_SpmSfSfp(i1,i3,i2))**2 * B0m2 
     If (WriteOut) Write(ErrCan,*) "Spm Sf Sfp",i1,i3,i2,sumI(i3,i3)
     res = res + sumI
    End Do
   End If !
  End Do

  res = oo16pi2 * res  

  Iname = Iname - 1
   
 End Subroutine PiSneutrino 


 Subroutine PiSquark(p2, mSf2, mf, mf2, mglu, c_GQSq, c_Sq4g3, c_Sq4e, c_Sq4Z &
         & , c_Sq2Z, c_Sq2W, c_Sq2g3, mSfp2, c_Sq4W, mfp, mfp2                &
         & , mN, mN2, c_FNSf_L, c_FNSf_R, mC, mC2, c_CFpSf_L, c_CFpSf_R       &
         & , mP02, c_P0SqSq, c_P0P0SqSq, mS02, c_S0SqSq, c_S0S0SqSq           &
         & , mSpm2, c_SpmSqSqp, c_SpmSpmSqSqp, mSup2, c_SuSf4                 &
         & , mSdown2, c_SdSf4, mSlepton2, c_SlSf4, mSneut2, c_SnSf4, mZ2, mW2 &
          & , WriteOut, res)
 !-------------------------------------------------------------------------
 ! calculates the 1-loop self energies of squarks
 ! The formulae of J. Bagger et al, Nucl. Phys. B491 are used.
 ! Note, that gSU3 has to be set to 0 for sleptons
 ! the input is:
 ! - p2 ......... the outer momentum squared 
 ! 
 ! the output is:
 !  res(i,j)
 ! written by Werner Porod, 31.12.01
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: p2, mSf2(6), mf(3), mf2(3), mSfp2(6), mfp(3)   &
     & , mfp2(3), mglu, mN(:), mN2(:), mC(:), mC2(:), mP02(:), mS02(:)   &
     & , mSpm2(:), mSup2(6), mSdown2(6), mSlepton2(:), mSneut2(:)        &
     & , c_Sq2Z(6), c_Sq2W(6), mZ2, mW2
  Complex(dp), Intent(in) :: c_Sq4g3(6,6,6), c_FNSf_L(:,:,:), c_FNSf_R(:,:,:) &
     & , c_CFpSf_L(:,:,:), c_CFpSf_R(:,:,:), c_GQSq(3,6,6), c_Sq4e(6,6,6)     &
     & , c_Sq4Z(6,6,6), c_Sq4W(6,6,6), c_P0P0SqSq(:,:,:), c_P0SqSq(:,:,:)     &
     & , c_S0S0SqSq(:,:,:), c_S0SqSq(:,:,:), c_SpmSqSqp(:,:,:)                &
     & , c_SpmSpmSqSqp(:,:,:), c_SuSf4(6,6,6), c_SdSf4(6,6,6), c_SlSf4(6,6,6) &
     & , c_SnSf4(3,6,6), c_Sq2g3(6,6,6)
  Logical, Intent(in) :: WriteOut
  Complex(dp), Intent(out) :: res(6,6)

  Integer :: i1, i2, i3, i2_min, i2_max, i3_min, i3_max, i4, i4_min, i4_max &
     & , ind_1, ind_2
  Complex(dp) :: sumI(6,6)
  Real(dp) :: mglu2, A0m2, F0m2, B0m2, G0m2

 !-----------------------
 ! initialisation
 !-----------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "PiSquark"

  mglu2 = mglu**2 

  res = ZeroC
  !------------------------
  ! strong interactions
  !------------------------
  Do i1=1,3 ! Generations of fermions
   G0m2 = Gloop(p2, mglu2, mf2(i1) )
   B0m2 = 2._dp * mf(i1) * mglu * Real( B0(p2, mglu2, mf2(i1) ), dp)
   sumI = ZeroC
   If (GenerationMixing) Then
    Do i2=1,3
     Do i3=1,3
      sumI(i2,i3) = c_GQSq(i1,i2,i3) * G0m2
      sumI(i2+3,i3+3) = c_GQSq(i1,i2+3,i3+3) * G0m2
     End Do
     Do i3=1,3
      sumI(i2,i3+3) = c_GQSq(i1,i2,i3+3) * B0m2
      sumI(i2+3,i3) = c_GQSq(i1,i2+3,i3) * B0m2
     End Do
    End Do
    i2_min = 1 
    i2_max = 6
    i3_min = i2_min
    i3_max = i2_max

   Else
    i2_min = 2*i1-1
    i2_max = 2*i1
    i3_min = i2_min
    i3_max = i2_max
    Do i2=i2_min,i2_max
     Do i3=i3_min,i3_max
      If (i2.Eq.i3) Then
       sumI(i2,i3) = c_GQSq(i1,i2,i3) * G0m2
      Else
       sumI(i2,i3) = c_GQSq(i1,i2,i3) * B0m2
      End If
     End Do
    End Do
   End If

   If (WriteOut) Then
    Write(ErrCan,*) "gluino quark",i1
    Do i2=i2_min,i2_max
     Write(ErrCan,*) "SumI(i2,:),",i2,sumI(i2,i3_min:i3_max)
    End Do
   End If
   Call CheckHermitian(sumI,"gluino quark "//Bu(i1) )
   res = res + sumI
  End Do

  Do i1=1,6
   A0m2 = A0(msf2(i1))
   F0m2 = Floop(p2, msf2(i1), 0._dp )
   If (GenerationMixing) Then
    i2_min = 1
    i2_max = 6
   Else
    If (i1.Le.2) Then
     i2_min = 1
    Elseif (i1.Le.4) Then
     i2_min = 3
    Else
     i2_min = 5
    End If
    i2_max = i2_min + 1
   End If
   i3_min = i2_min
   i3_max = i2_max
   sumI = ZeroC
   Do i2=i2_min,i2_max
    Do i3=i3_min,i3_max
     sumI(i2,i3) = c_Sq2g3(i1,i2,i3) * F0m2 + c_Sq4g3(i1,i2,i3) * A0m2
    End Do
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "squark gluon",i1
    Do i2=i2_min,i2_max
     Write(ErrCan,*) "SumI(i2,:),",i2,sumI(i2,i3_min:i3_max)
    End Do
   End If
   res = res + sumI
   Call CheckHermitian(sumI,"squark gluon "//Bu(i1) )
  End Do
 !--------------------------------
 ! electroweak gauge bosons
 !--------------------------------
  A0m2 = A0(mZ2)
  sumI = ZeroC
  Do i1=1,6
   sumI(i1,i1) = c_Sq2Z(i1) * A0m2
  End Do
  If (WriteOut) Write(ErrCan,*) "Sq Sq Z Z",(sumI(i1,i1),i1=1,6)
  res = res + sumI

  Do i1=1,6
   F0m2 = Floop(p2, msf2(i1), 0._dp )
   If (GenerationMixing) Then
    i2_min = 1
    i2_max = 6
   Else
    If (i1.Le.2) Then
     i2_min = 1
    Elseif (i1.Le.4) Then
     i2_min = 3
    Else
     i2_min = 5
    End If
    i2_max = i2_min + 1
   End If
   i3_min = i2_min
   i3_max = i2_max
   sumI = ZeroC
   Do i2=i2_min,i2_max
    Do i3=i3_min,i3_max
     sumI(i2,i3) = c_Sq4e(i1,i2,i3) * F0m2
    End Do
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "squark photon",i1
    Do i2=i2_min,i2_max
     Write(ErrCan,*) "SumI(i2,:),",i2,sumI(i2,i3_min:i3_max)
    End Do
   End If
   Call CheckHermitian(sumI,"squark photon "//Bu(i1) )
   res = res + sumI

   F0m2 = Floop(p2, msf2(i1), mZ2 )
   sumI = ZeroC
   Do i2=i2_min,i2_max
    Do i3=i3_min,i3_max
     sumI(i2,i3) = c_Sq4Z(i1,i2,i3) * F0m2
    End Do
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "squark Z",i1
    Do i2=i2_min,i2_max
     Write(ErrCan,*) "SumI(i2,:),",i2,sumI(i2,i3_min:i3_max)
    End Do
   End If
   Call CheckHermitian(sumI,"squark Z "//Bu(i1) )
   res = res + sumI
  End Do

  A0m2 = A0(mW2)
  sumI = ZeroC
  Do i1=1,6
   sumI(i1,i1) = c_Sq2W(i1) * A0m2
  End Do
  If (WriteOut) Write(ErrCan,*) "Sq Sq W W",(sumI(i1,i1),i1=1,6)
  res = res + sumI

  Do i1=1,6
   F0m2 = Floop(p2, msfp2(i1), mW2 )
   If (GenerationMixing) Then
    i2_min = 1
    i2_max = 3
   Else
    If (i1.Le.2) Then
     i2_min = 1
    Elseif (i1.Le.4) Then
     i2_min = 3
    Else
     i2_min = 5
    End If
    i2_max = i2_min
   End If
   i3_min = i2_min
   i3_max = i2_max
   sumI = ZeroC
   Do i2=i2_min,i2_max
    Do i3=i3_min,i3_max
     sumI(i2,i3) = c_Sq4W(i1,i2,i3) * F0m2
    End Do
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "squark' W",i1
    Do i2=i2_min,i2_max
     Write(ErrCan,*) "SumI(i2,:),",i2,sumI(i2,i3_min:i3_max)
    End Do
   End If
   Call CheckHermitian(sumI,"squark W "//Bu(i1) )
   res = res + sumI
     If (Aimag(res(6,6)).Gt.0.1) Then
      Write(errcan,*) "W+ sf'",i1
      Call terminateprogram
     End If
  End Do

  !-----------------------
  ! neutralinos
  !-----------------------
  Do i1=1,n_neut
   Do i2=1,3
    G0m2 = Gloop(p2, mN2(i1), mf2(i2) )
    B0m2 = - 2._dp * mN(i1) * mf(i2) * B0(p2, mN2(i1), mf2(i2) )
    If (GenerationMixing) Then
     i3_min = 1
     i3_max = 6
    Else
     i3_min = 2*i2-1
     i3_max = 2*i2
    End If
    i4_min = i3_min
    i4_max = i3_max
    SumI = ZeroC
    Do i3=i3_min,i3_max
     Do i4=i4_min,i4_max
      SumI(i3,i4) = ( Conjg(c_FNSf_L(i2,i1,i3))* c_FNSf_L(i2,i1,i4)           &
                &   + Conjg(c_FNSf_R(i2,i1,i3))* c_FNSf_R(i2,i1,i4) ) * G0m2  &
                & + ( c_FNSf_L(i2,i1,i4) * Conjg(c_FNSf_R(i2,i1,i3))          &
                &   + c_FNSf_R(i2,i1,i4) * Conjg(c_FNSf_L(i2,i1,i3)) ) * B0m2
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "chi0 f",i1,i2
     Do i3=i3_min,i3_max
      Write(ErrCan,*) "SumI(i3,:)",SumI(i3,i4_min:i4_max)
     End Do
    End If
    Call CheckHermitian(sumI,"chi0 f "//Bu(i1)//" "//Bu(i2) )
   res = res + SumI
     If (Aimag(res(6,6)).Gt.0.1) Then
      Write(errcan,*) "chi0 f",i1,i2
      Call terminateprogram
     End If
   End Do
  End Do

 !-----------------------
 ! charginos
 !-----------------------
  Do i1=1,n_char
   Do i2=1,3
    G0m2 = Gloop(p2, mC2(i1), mfp2(i2) )
    B0m2 = - 2._dp * mC(i1) * mfp(i2) * B0(p2, mC2(i1), mfp2(i2) )
    If (GenerationMixing) Then
     i3_min = 1
     i3_max = 6
    Else
     i3_min = 2*i2-1
     i3_max = 2*i2
    End If
    i4_min = i3_min
    i4_max = i3_max
    SumI = ZeroC
    Do i3=i3_min,i3_max
     Do i4=i4_min,i4_max
      SumI(i3,i4) = ( Conjg(c_CFpSf_L(i1,i2,i3))* c_CFpSf_L(i1,i2,i4)         &
                &  + Conjg(c_CFpSf_R(i1,i2,i3))* c_CFpSf_R(i1,i2,i4) ) * G0m2 &
                & + ( c_CFpSf_L(i1,i2,i4) * Conjg(c_CFpSf_R(i1,i2,i3))        &
                &   + c_CFpSf_R(i1,i2,i4) * Conjg(c_CFpSf_L(i1,i2,i3)) ) * B0m2
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "chi+ f'",i1,i2
     Do i3=i3_min,i3_max
      Write(ErrCan,*) "SumI(i3,:)",SumI(i3,i4_min:i4_max)
     End Do
    End If
    Call CheckHermitian(sumI,"chi+ f "//Bu(i1)//" "//Bu(i2) )
    res = res + SumI
     If (Aimag(res(6,6)).Gt.0.1) Then
      Write(errcan,*) "chi+ f'",i1,i2
      Call terminateprogram
     End If
   End Do
  End Do

  !-------------------------------
  ! sfermion quartic interactions
  !-------------------------------
  ! up squarks
  !------------
  Do i1=1,6
   A0m2 = A0( mSup2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,6
     Do i3=1,6
      SumI(i2,i3) = - c_SuSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Su Su Sf Sf",i1
     Do i2=1,6
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    Call CheckHermitian(sumI,"Su Su Sf Sf "//Bu(i1) )
    res = res + sumI
      If (Aimag(res(6,6)).Gt.0.1) Then
      Write(errcan,*) "Su Su Sf Sf",i1,i2
      Call terminateprogram
     End If

   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,3
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_SuSf4(i1,ind_1,ind_2) * A0m2
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Su Su Sf Sf",i1
     Do i2=1,6
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
   End If !
  End Do

  !--------------
  ! down squarks
  !--------------
  Do i1=1,6
   A0m2 = A0( mSdown2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,6
     Do i3=1,6
      SumI(i2,i3) = - c_SdSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sd Sd Sf Sf",i1
     Do i2=1,6
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
     If (Aimag(res(6,6)).Gt.0.1) Then
      Write(errcan,*) "Sd Sd Sf Sf",i1,i2
      Call terminateprogram
     End If
    res = res + sumI
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,3
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_SdSf4(i1,ind_1,ind_2) * A0m2
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sd Sd Sf Sf",i1
     Do i2=1,6
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
   End If !
  End Do

  !------------------
  ! charged sleptons
  !------------------
  Do i1=1,n_slept
   A0m2 = A0( mSlepton2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,6
     Do i3=1,6
      SumI(i2,i3) = - c_SlSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sl Sl Sf Sf",i1
     Do i2=1,6
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    res = res + sumI
     If (Aimag(res(6,6)).Gt.0.1) Then
      Write(errcan,*) "Sl Sl Sf Sf",i1,i2
      Call terminateprogram
     End If
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,3
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_SlSf4(i1,ind_1,ind_2) * A0m2
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sl Sl Sf Sf",i1
     Do i2=1,6
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
   End If !
  End Do

  !----------------
  ! sneutrinos
  !----------------
  Do i1=1,n_sneut
   A0m2 = A0( mSneut2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,6
     Do i3=1,6
      SumI(i2,i3) = - c_SnSf4(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sn Sn Sf Sf",i1
     Do i2=1,6
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    res = res + sumI
     If (Aimag(res(6,6)).Gt.0.1) Then
      Write(errcan,*) "Sn sn Sf Sf",i1,i2
      Call terminateprogram
     End If
 
   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,3
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_SnSf4(i1,ind_1,ind_2) * A0m2
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Sn Sn Sf Sf",i1
     Do i2=1,6
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
   End If !
  End Do

  !---------------------
  ! pseudoscalar Higgs
  !---------------------
  Do i1=1,n_P0
   A0m2 = A0( mP02(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,6
     Do i3=1,6
      SumI(i2,i3) = - c_P0P0SqSq(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "P0 P0 Sf Sf",i1
     Do i2=1,6
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    Call CheckHermitian(sumI,"P0 P0 Sf Sf "//Bu(i1) )
    res = res + sumI
 
    Do i2=1,6
     SumI = ZeroC
     B0m2 = B0(p2, mP02(i1), mSf2(i2) )
     Do i3=1,6
      Do i4=1,6
       SumI(i3,i4) = Conjg(c_P0SqSq(i1,i3,i2)) * c_P0SqSq(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "P0 Sf Sf",i1,i2
      Do i3=1,6
       Write(ErrCan,*) "SumI(i1,:)",i3,sumI(i3,:)
      End Do
     End If
     Call CheckHermitian(sumI,"P0 Sf Sf "//Bu(i1)//" "//Bu(i2) )
     res = res + sumI
     If (Aimag(res(6,6)).Gt.0.1) Then
      Write(errcan,*) "P0 Sf Sf",i1,i2
      Call terminateprogram
     End If
    End Do

   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,3
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_P0P0SqSq(i1,ind_1,ind_2) * A0m2 
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "P0 P0 Sf Sf",i1
     Do i2=1,6
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
 
    Do i2=1,6
     SumI = ZeroC
     B0m2 = B0(p2, mP02(i1), mSf2(i2) )
     If (i2.Le.2) Then
      i3_min=1
     Else If (i2.Le.4) Then
      i3_min=3
     Else 
      i3_min=5
     End If
     i3_max = i3_min + 1
     i4_min = i3_min
     i4_max = i3_max
     Do i3=i3_min, i3_max
      Do i4=i4_min, i4_max
       SumI(i3,i4) = Conjg(c_P0SqSq(i1,i3,i2)) * c_P0SqSq(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "P0 Sf Sf",i1,i2
      Do i3=i3_min, i3_max
       Write(ErrCan,*) "SumI(i1,:)",sumI(i3,i4_min:i4_max)
      End Do
     End If
     res = res + sumI
    End Do
   End If !
  End Do
  !---------------------
  ! neutral Higgs
  !---------------------
  Do i1=1,n_S0
   A0m2 = A0( mS02(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,6
     Do i3=1,6
      SumI(i2,i3) = - c_S0S0SqSq(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "S0 S0 Sf Sf",i1
     Do i2=1,6
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    Call CheckHermitian(sumI,"S0 S0 Sf Sf "//Bu(i1) )
    res = res + sumI
 
    Do i2=1,6
     SumI = ZeroC
     B0m2 = B0(p2, mS02(i1), mSf2(i2) )
     Do i3=1,6
      Do i4=1,6
       SumI(i3,i4) = Conjg(c_S0SqSq(i1,i3,i2)) * c_S0SqSq(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "S0 Sf Sf",i1,i2
      Do i3=1,6
       Write(ErrCan,*) "SumI(i1,:)",i3,sumI(i3,:)
      End Do
     End If
     Call CheckHermitian(sumI,"S0 Sf Sf "//Bu(i1)//" "//Bu(i2) )
     res = res + sumI
     If (Aimag(res(6,6)).Gt.0.1) Then
      Write(errcan,*) "S0 Sf Sf",i1,i2
      Call terminateprogram
     End If
    End Do

   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,3
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_S0S0SqSq(i1,ind_1,ind_2) * A0m2 
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "S0 S0 Sf Sf",i1
     Do i2=1,6
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
 
    Do i2=1,6
     SumI = ZeroC
     B0m2 = B0(p2, mS02(i1), mSf2(i2) )
     If (i2.Le.2) Then
      i3_min=1
     Else If (i2.Le.4) Then
      i3_min=3
     Else 
      i3_min=5
     End If
     i3_max = i3_min + 1
     i4_min = i3_min
     i4_max = i3_max
     Do i3=i3_min, i3_max
      Do i4=i4_min, i4_max
       SumI(i3,i4) = Conjg(c_S0SqSq(i1,i3,i2)) * c_S0SqSq(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "S0 Sf Sf",i1,i2
      Do i3=i3_min, i3_max
       Write(ErrCan,*) "SumI(i1,:)",sumI(i3,i4_min:i4_max)
      End Do
     End If
     res = res + sumI
    End Do
   End If !
  End Do
  !---------------------
  ! charged Higgs
  !---------------------
  Do i1=1,n_Spm
   A0m2 = A0( mSpm2(i1) )
   If (GenerationMixing) Then
    SumI = ZeroC
    Do i2=1,6
     Do i3=1,6
      SumI(i2,i3) = - c_SpmSpmSqSqp(i1,i2,i3) * A0m2 
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Spm Spm Sf Sf",i1
     Do i2=1,6
      Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,:)
     End Do
    End If
    Call CheckHermitian(sumI,"Spm Spm Sf Sf "//Bu(i1) )
    res = res + sumI
 
    Do i2=1,6
     SumI = ZeroC
     B0m2 = B0(p2, mSpm2(i1), mSfp2(i2) )
     Do i3=1,6
      Do i4=1,6
       SumI(i3,i4) = Conjg(c_SpmSqSqp(i1,i3,i2)) * c_SpmSqSqp(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "Spm Sf Sfp",i1,i2
      Do i3=1,6
       Write(ErrCan,*) "SumI(i1,:)",i3,sumI(i3,:)
      End Do
     End If
     Call CheckHermitian(sumI,"Spm Sf Sfp "//Bu(i1)//" "//Bu(i2) )
     res = res + sumI
     If (Aimag(res(6,6)).Gt.0.1) Then
      Write(errcan,*) "Spm Sf Sfp",i1,i2
      Call terminateprogram
     End If
    End Do

   Else ! .not. GenerationMixing
    SumI = ZeroC
    Do i2=1,3
     Do i3=1,2
      ind_1 = 2*(i2-1)+i3
      Do i4=1,2
       ind_2 = 2*(i2-1)+i4
       SumI(ind_1,ind_2) = - c_SpmSpmSqSqp(i1,ind_1,ind_2) * A0m2 
      End Do
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "Spm Spm Sf Sf",i1
     Do i2=1,6
      If (i2.Le.2) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,1:2)
      Else If (i2.Le.4) Then
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,3:4)
      Else
       Write(ErrCan,*) "SumI(i1,:)",i2,sumI(i2,5:6)
      End If
     End Do
    End If
    res = res + sumI
 
    Do i2=1,6
     SumI = ZeroC
     B0m2 = B0(p2, mSpm2(i1), mSfp2(i2) )
     If (i2.Le.2) Then
      i3_min=1
     Else If (i2.Le.4) Then
      i3_min=3
     Else 
      i3_min=5
     End If
     i3_max = i3_min + 1
     i4_min = i3_min
     i4_max = i3_max
     Do i3=i3_min, i3_max
      Do i4=i4_min, i4_max
       SumI(i3,i4) = Conjg(c_SpmSqSqp(i1,i3,i2)) * c_SpmSqSqp(i1,i4,i2) * B0m2 
      End Do
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "Spm Sf Sfp",i1,i2
      Do i3=i3_min, i3_max
       Write(ErrCan,*) "SumI(i1,:)",sumI(i3,i4_min:i4_max)
      End Do
     End If
     res = res + sumI
    End Do
   End If !
  End Do
  res = oo16pi2 * res  

  Iname = Iname - 1
   
 End Subroutine PiSquark 


 Subroutine PiWWT1(p2, gSU2, sinW2, mS02, RS0, mSpm2, RSpm, vevs            &
         & , mP02, RP0, mSneutrino2, RSneutrino, mSlepton2, RSlepton        &
         & , mUSquark2, RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2  &
         & , CKM, mN, mN2, N, mC, mC2, U, V, mZ2, mW2, res)
 !-----------------------------------------------------------------------
 ! Calculates the 1-loop self energy of Z-boson. 
 ! The formula of J. Bagger et al, Nucl.Phys.B is used. The renormalization
 ! scale has to be set in the main program due to the structure of LoopTools. 
 ! the input is:
 ! - p2 ...... the outer momentum squared 
 ! - gSU2 .... the SU(2) gauge coupling
 ! - sinW2 ... sin(theta_W) squared
 ! - mZ2 ..... Z-boson mass
 ! - mW2 ..... W-boson mass
 ! - mS02 .... masses of the scalar Higgs
 ! - RS0 ..... mixing matrix of the scalar Higgs
 ! - mP02 .... masses of the pseudoscalar Higgs
 ! - RP0 ..... mixing matrix of the pseudoscalar Higgs
 ! - mSpm2 ... masses of the charged Higgs
 ! - RSpm .... mixing matrix of the charged Higgs
 ! - mN ...... masses of the neutralinos
 ! - N ....... mixing matrix of the neutrlinos
 ! - mC ...... masses of the charginos
 ! - U,V ..... mixing matrices of the charginos
 ! written by Werner Porod, 19.7.99
 ! 10.10.01: portation to f90
 ! 14.01.02: - adding factor 3 in case of squarks
 !           - decoupling neutralino and chargino indices
 ! 14.12.03: correcting incorrect sneutrino index
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: p2, gSU2, sinW2, mS02(:), RS0(:,:), mSpm2(:) &
    & , mP02(:), RP0(:,:), mSneutrino2(:), mSlepton2(:), mUSquark2(6)  &
    & , mDSquark2(6), mN(:), mC(:), mN2(:), mC2(:), vevs(:), mZ2, mW2  &
    & , mf_l2(3), mf_d2(3), mf_u2(3)
  Complex(dp), Intent(in) :: N(:,:), U(:,:), V(:,:), RSpm(:,:)         &
    &  , RUSquark(6,6), RDSquark(6,6), RSlepton(:,:), RSneutrino(:,:)  &
    &  , CKM(3,3)
  Complex(dp), Intent(out) :: res 

  Integer :: i1, i2, i3
  Real(dp) :: cosW2, g2, coup
  Complex(dp) :: sumI, coupC, coupLC, coupRC

  Logical :: WriteOut

  Iname = Iname + 1
  NameOfUnit(Iname) = "PiWWT1"

  !----------------------------------
  ! In case I need the contributions
  !----------------------------------
  If ((WriteOneLoopContributions.Eq.5).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions in PiWWT1:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If
 !-----------------------
 ! intialisation
 !-----------------------
  res = (0._dp,0._dp)
  cosW2 = 1._dp - sinW2
  g2 = gSU2**2
 !--------------------------
 ! gauge bosons
 !--------------------------
  sumI = - g2 * ( ( (4._dp * p2 + mW2 + mZ2) * cosW2 - mZ2 * sinW2**2 )  &
     &            * B0(p2,mZ2,mW2)  &
     &          + ( 1._dp + 8._dp * cosW2 ) * B22(p2,mZ2,mW2) )
  If (WriteOut) Write(ErrCan,*) "Z",sumI
  res = res + sumI

  sumI = - g2 * sinW2 * ( 8._dp * B22(p2,mW2,0._dp)     &
       &                +  4._dp * p2 * B0(p2,mW2,0._dp) )
  If (WriteOut) Write(ErrCan,*) "photon",sumI
  res = res + sumI
  !--------------------------
  !  Higgs bosons
  !--------------------------
  Do i1=1,n_S0
   Do i2=1,n_Spm
    Call CoupChargedScalarScalarW(i2, i1, gSU2, RSpm, RS0, coupC)
    sumI = - 4._dp * Abs(coupC)**2 * B22(p2,mS02(i1), mSpm2(i2) )
    If (WriteOut) Write(ErrCan,*) "S0 S+",i1,i2,sumI
    res = res + sumI
   End Do
   Call CoupScalarW(i1, gSU2, vevs, RS0, coup)
   sumI = Abs(coup)**2 * B0(p2,mW2,mS02(i1))
   If (WriteOut) Write(ErrCan,*) "S0",i1,sumI
   res = res + sumI
  End Do
  !--------------------------
  ! pseudoscalar Higgs
  !--------------------------
  Do i1=2,n_P0
   Do i2=2,n_Spm
    Call CoupChargedScalarPseudoScalarW(i2, i1, gSU2, RSpm, RP0, coupC)
    sumI = - 4._dp * Abs(coupC)**2 * B22(p2,mP02(i1), mSpm2(i2) )
    If (WriteOut) Write(ErrCan,*) "P0 S+",i1,i2,sumI
    res = res + sumI
   End Do
  End Do
  !-----------------------
  ! Sneutrinos + Sleptons
  !-----------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1 = 1, n_Sneut ! number of sneutrinos/sleptons depend on model
    Do i2 =1, n_Slept
     Call CoupSfermionW3(i2, i1, gSU2, RSlepton, RSneutrino, coupC )
     sumI = - 4._dp * Abs(coupC)**2 * B22(p2, mSneutrino2(i1), mSlepton2(i2) )
     If (WriteOut) Write(ErrCan,*) "Snu Sle",i1,i2,sumI
     res = res + sumI
    End Do
   End Do

  Else 
#endif
   Do i1 = 1, n_Sneut ! number of sneutrinos/sleptons depend on model
    Do i2 =1,2
     sumI = - 2._dp * g2 * Abs(RSlepton((i1-1)*2+i2,2*i1-1))**2   &
          &         * B22(p2, mSneutrino2(i1), mSlepton2((i1-1)*2+i2) )
     If (WriteOut) Write(ErrCan,*) "Snu Sle",i1,i2,sumI
     res = res + sumI
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif
  !--------------------------
  ! Scalar up /down
  !--------------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1 = 1, 6
    Do i2 =1, 6
     Call CoupSfermionW3(i2, i1, gSU2, RDSquark, RUSquark, coupC )
     sumI = - 12._dp * Abs(coupC)**2 * B22(p2, mDSquark2(i2), mUSquark2(i1) )
     If (WriteOut) Write(ErrCan,*) "Sd Su",i1,i2,sumI
     res = res + sumI
    End Do
   End Do

  Else 
#endif
   Do i1 = 1,3
    Do i2 = 1,2
     Do i3 = 1,2
      sumI = -6._dp * g2 * Abs(RUSquark((i1-1)*2+i2,2*i1-1)               &
           &                   * Conjg(RDSquark((i1-1)*2+i3,2*i1-1)))**2  &
           &        * B22(p2, mUSquark2((i1-1)*2+i2), mDSquark2((i1-1)*2+i3))
      If (WriteOut) Write(ErrCan,*) "Sd Su",i1,i2,i3,sumI
      res = res + sumI
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
 End If
#endif
  !--------------------------
  ! neutrinos/leptons
  !--------------------------
  Do i1 = 1,5-n_char ! number of leptons depend on model
   sumI = 0.5_dp * g2 * HLoop(p2,mf_l2(i1), 0._dp)
   If (WriteOut) Write(ErrCan,*) "Nu L",i1,sumI
   res = res + sumI
  End Do
  !--------------------------
  ! quarks
  !--------------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     sumI = 1.5_dp * g2 * Abs(CKM(i1,i2))**2 * HLoop(p2, mf_u2(i1), mf_d2(i2)) 
     If (WriteOut) Write(ErrCan,*) "U D",i1,i2,sumI
     res = res + sumI
    End Do
   End Do
  Else 
#endif
   Do i1 = 1,3
    sumI = 1.5_dp * g2 * HLoop(p2, mf_u2(i1), mf_d2(i1)) 
    If (WriteOut) Write(ErrCan,*) "U D",i1,sumI
    res = res + sumI
   End Do
#ifdef GENERATIONMIXING
  End If
#endif
  !--------------------------
  ! charginos/neutralinos
  !--------------------------
  Do i1 = 1,n_char
   Do i2 = 1,n_neut
    Call CoupCharginoNeutralinoW(i1, i2, N, U, V, gSU2, coupLC, coupRC)
    sumI = ( Abs(coupLC)**2 + Abs(coupRC)**2 ) * HLoop(p2,mC2(i1),mN2(i2))  &
       & + 4._dp * Real( Conjg(coupLC) * coupRC,dp ) * mC(i1) * mN(i2)      &
       &         * B0(p2, mC2(i1), mN2(i2))
     If (WriteOut) Write(ErrCan,*) "Char Neut",i1,i2,sumI
     res = res + sumI
   End Do
  End Do

  res = oo16pi2 * res
  Iname = Iname - 1

 End Subroutine PiWWT1


 Subroutine PiWWT1_SM(p2, gSU2, sinW2, mS02, vev, mf_l2, mf_u2, mf_d2  &
         & , CKM, mZ2, mW2, res)
 !-----------------------------------------------------------------------
 ! Calculates the 1-loop self energy of W-boson. 
 ! The formula of J. Bagger et al, Nucl.Phys.B is used. The renormalization
 ! scale has to be set in the main program due to the structure of LoopTools. 
 ! the input is:
 ! - p2 ...... the outer momentum squared 
 ! - gSU2 .... the SU(2) gauge coupling
 ! - sinW2 ... sin(theta_W) squared
 ! - mZ2 ..... Z-boson mass
 ! - mW2 ..... W-boson mass
 ! - mS02 .... masses of the scalar Higgs
 ! written by Werner Porod, 19.7.99
 ! 10.10.01: portation to f90
 ! 14.01.02: - adding factor 3 in case of squarks
 !           - decoupling neutralino and chargino indices
 ! 14.12.03: correcting incorrect sneutrino index
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: p2, gSU2, sinW2, mS02, vev, mZ2, mW2, mf_l2(3) &
                      & , mf_d2(3), mf_u2(3)
  Complex(dp), Intent(in) :: CKM(3,3)
  Complex(dp), Intent(out) :: res 

  Integer :: i1, i2, i3
  Real(dp) :: cosW2, g2, coup
  Complex(dp) :: sumI, coupC, coupLC, coupRC

  Logical :: WriteOut

  Iname = Iname + 1
  NameOfUnit(Iname) = "PiWWT1_SM"

  !----------------------------------
  ! In case I need the contributions
  !----------------------------------
  If ((WriteOneLoopContributions.Eq.5).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions in PiWWT1_SM:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If
 !-----------------------
 ! intialisation
 !-----------------------
  res = (0._dp,0._dp)
  cosW2 = 1._dp - sinW2
  g2 = gSU2**2
 !--------------------------
 ! gauge bosons
 !--------------------------
  sumI = - g2 * ( ( (4._dp * p2 + mW2 + mZ2) * cosW2 - mZ2 * sinW2**2 )  &
     &            * B0(p2,mZ2,mW2)  &
     &          + ( 1._dp + 8._dp * cosW2 ) * B22(p2,mZ2,mW2) )
  If (WriteOut) Write(ErrCan,*) "Z",sumI
  res = res + sumI

  sumI = - g2 * sinW2 * ( 8._dp * B22(p2,mW2,0._dp)     &
       &                +  4._dp * p2 * B0(p2,mW2,0._dp) )
  If (WriteOut) Write(ErrCan,*) "photon",sumI
  res = res + sumI
  !--------------------------
  !  Higgs bosons
  !--------------------------
  coup = 0.5_dp * g2 * vev
  sumI = Abs(coup)**2 * B0(p2,mW2,mS02)
  If (WriteOut) Write(ErrCan,*) "S0",sumI
  res = res + sumI
  !--------------------------
  ! neutrinos/leptons
  !--------------------------
  Do i1 = 1,3 ! number of leptons depend on model
   sumI = 0.5_dp * g2 * HLoop(p2,mf_l2(i1), 0._dp)
   If (WriteOut) Write(ErrCan,*) "Nu L",i1,sumI
   res = res + sumI
  End Do
  !--------------------------
  ! quarks
  !--------------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     sumI = 1.5_dp * g2 * Abs(CKM(i1,i2))**2 * HLoop(p2, mf_u2(i1), mf_d2(i2)) 
     If (WriteOut) Write(ErrCan,*) "U D",i1,i2,sumI
     res = res + sumI
    End Do
   End Do
  Else 
#endif
   Do i1 = 1,3
    sumI = 1.5_dp * g2 * HLoop(p2, mf_u2(i1), mf_d2(i1)) 
    If (WriteOut) Write(ErrCan,*) "U D",i1,sumI
    res = res + sumI
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  res = oo16pi2 * res
  Iname = Iname - 1

 End Subroutine PiWWT1_SM


 Subroutine PiZZT1(p2, gSU2, sinW2, vevSM, mZ2, mW2, mS02, RS0, mP02, RP0  &
   &    , mSpm2, RSpm, mSneut2, RSneut, mSlepton2, RSlepton, mSup2, RSup   &
   &    , mSdown2, RSdown, mf_l2, mf_u2, mf_d2, mC, mC2, U, V, mN, mN2, N  &
   &    , res, NMSSM)
 !-----------------------------------------------------------------------
 ! Calculates the 1-loop self energy of Z-boson. 
 ! The formula of J. Bagger et al, Nucl.Phys.B is used. The renormalization
 ! scale has to be set in the calling routine due to the structure of
 ! LoopTools. 
 ! the input is:
 !  - p2 ......... the outer momentum squared 
 !  - gSU2 ....... the SU(2) gauge coupling at p2
 !  - sinW2 ...... sin(theta_W) squared
 !  - vevSM(2) ... MSSM vevs
 !  - mZ2 ........ Z-boson mass squared
 !  - mW2 ........ W-boson mass squared
 !  - mS02 ....... masses of the scalar Higgs squared
 !  - RS0 ........ mixing matrix of the scalar Higgs
 !  - mP02 ....... masses of the pseudoscalar Higgs squared
 !  - RP0 ........ mixing matrix of the pseudoscalar Higgs
 !  - mSpm2 ...... masses of the charged Higgs squared
 !  - RSpm ....... mixing matrix of the charged scalar Higgs
 !  - mSneut2 .... masses of the sneutrinos squared
 !  - RSneut ..... mixing matrix of the sneutrinos
 !  - mSlepton2 .. masses of the sleptons squared
 !  - RSlepton ... mixing matrix of the sleptons  
 !  - mSup2 ...... masses of the s-ups squared
 !  - RSup ....... mixing matrix of the s-ups
 !  - mSdown2 .... masses of the s-downs squared
 !  - RSdown ..... mixing matrix of the s-downs
 !  - mC ......... masses of the charginos
 !  - mC2 ........ masses of the charginos squared
 !  - U, V ....... mixing matrices of charginos
 !  - mN ......... masses of the neutralinos
 !  - mN2 ........ masses of the neutralinos squared
 !  - N .......... mixing matrix of neutralinos
 !  - NMSSM ...... logical variable, has to be set true in case of the NMSSM
 !                 so that the correct routine for the couplings of the 
 !                 Z-bosons to neutralinos is used.
 ! The output:
 !  - res    
 ! written by Werner Porod, 19.7.99
 ! 08.10.01: portation to f90, changing interface
 ! 14.01.02: adding a factor 3 in case of squarks
 ! 12.12.03: adding mW2 and mZ2 to interface
 ! 14.11.08: adding NMSSM
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: p2, gSU2, sinW2, vevSM(2), mS02(:), RS0(:,:)        &
        & , mP02(:), RP0(:,:), mSpm2(:), mSneut2(:), mSlepton2(:), mSup2(6)   &
        & , mSdown2(6), mC(:), mC2(:), mN(:), mN2(:), mW2, mZ2, mf_l2(3)      &
        & , mf_u2(3), mf_d2(3)
  Complex(dp), Intent(in) :: RSpm(:,:), RSlepton(:,:)  &
        & , RSup(6,6), RSdown(6,6), U(:,:), V(:,:), N(:,:)
  Complex(dp), Intent(in) :: RSneut(:,:)
  Complex(dp), Intent(out) :: res
  Logical, Intent(in), Optional :: NMSSM

  Integer :: i1, i2, i3
  Real(dp) :: cosW, cosW2, coup, coup2
  Complex(dp) :: sumI, coupC, coupC2, RSf(2,2)
  Logical :: WriteOut

  Iname = Iname + 1
  NameOfUnit(Iname) = "PiZZT1"

  res = ZeroC
  cosW2 = 1._dp - sinW2
  cosW = Sqrt( cosW2 )
  If (n_Char.Lt.2) n_char = Size( mC )
  !----------------------------------
  ! In case I need the contributions
  !----------------------------------
  If ((WriteOneLoopContributions.Eq.1).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions to PiZZT1:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If

  !--------------------------
  ! neutral Higgs bosons
  !--------------------------
  Do i1=1,n_S0
   Do i2=1,n_P0
    Call CoupPseudoscalarScalarZ(i2, i1, gSU2, cosW, RP0, RS0, coupC)
    sumI = - 4._dp * Abs(coupC)**2 * B22(p2, mP02(i2), mS02(i1) )
    If (WriteOut) Write(ErrCan,*) "S0, P0", i1,i2,sumI
    res = res + sumI
   End Do
   Call CoupScalarZ(i1, gSU2, cosW2, vevSM, RS0, coup )
   sumI = coup**2 * B0(p2, mZ2, mS02(i1))
   If (WriteOut) Write(ErrCan,*) "S0", i1,sumI
   res = res + sumI
  End Do
  !--------------------------
  ! W boson
  !--------------------------
  cosW2 = 1._dp - sinW2
  sumI = - gSU2**2 * ( 2._dp * ( (2._dp * p2 + mW2)*cosW2 - mZ2 * sinW2**2 ) &
   &                         * B0(p2,mW2,mW2)                                &
   &                 + ( 8._dp * cosW2 + (cosW2-sinW2)**2/cosW2 )            &
   &                    * B22(p2,mW2,mW2) )
  If (WriteOut) Write(ErrCan,*) "W", sumI
  res = res + sumI
  !--------------------------
  ! S+ bosons
  !--------------------------
  Do i1=2,n_Spm
   Do i2=2,n_Spm
    Call CoupChargedScalarZ(i1, i2, gSU2, sinW2, RSpm, coupC )
    sumI = - 4._dp * Abs(coupC)**2 * B22(p2, mSpm2(i1), mSpm2(i2)) 
    If (WriteOut) Write(ErrCan,*) "Spm", i1,i2,sumI
    res = res + sumI
   End Do
  End Do
  !--------------------------
  ! Sneutrinos
  !--------------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1 = 1, n_Sneut
    Do i2 = 1, n_Sneut
     Call CoupSneutrinoZ(i1, i2, gSU2, sinW2, Rsneut, coupC )
     sumI =  - 4._dp * Abs(coupC)**2  * B22(p2, mSneut2(i1), mSneut2(i2) )
     If (WriteOut) Write(ErrCan,*) "Sneutrino", i1,i2,sumI
     res = res + sumI
    End Do
   End Do
  Else
#endif
   Do i1 = 1, n_Sneut
    Call CoupSneutrinoZ(gSU2, sinW2, coupC )
    sumI =  - 4._dp * Abs(coupC)**2  * B22(p2, mSneut2(i1), mSneut2(i1) )
    If (WriteOut) Write(ErrCan,*) "Sneutrino", i1,sumI
    res = res + sumI
   End Do
#ifdef GENERATIONMIXING
  End If 
#endif
 !--------------------------
 ! Sleptons
 !--------------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1 = 1, n_Slept
    Do i2 = 1, n_Slept
     Call CoupSleptonZ(i1, i2, gSU2, sinW2, Rslepton, coupC )
     sumI =  - 4._dp * Abs(coupC)**2 * B22(p2, mSlepton2(i1), mSlepton2(i2) )
     If (WriteOut) Write(ErrCan,*) "Sleptons", i1,i2,sumI
     res = res + sumI
    End Do
   End Do
  Else
#endif
   Do i1 = 1, n_Slept, 2
    Rsf = RSlepton(i1:i1+1, i1:i1+1)
    Do i2=1,2
     Do i3=i2,2
      Call CoupSleptonZ(i2, i3, gSU2, sinW2, Rsf, coupC )
      sumI = - 4._dp * Abs(coupC)**2       &
           &          * B22(p2, mSlepton2(i1 - 1 + i2), mSlepton2(i1 - 1 + i3))
      If (i2.Ne.i3) sumI = 2._dp * sumI
      If (WriteOut) Write(ErrCan,*) "Sleptons", i1, i2,i3, sumI
      res = res + sumI
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If 
#endif
  !--------------------------
  ! Scalar up
  !--------------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1 = 1, 6
    Do i2 = 1, 6
     Call CoupSupZ(i1, i2, gSU2, sinW2, Rsup, coupC )
     sumI =  - 12._dp * Abs(coupC)**2 * B22(p2, mSup2(i1), mSup2(i2) )
     If (WriteOut) Write(ErrCan,*) "S-up", i1,i2,sumI
     res = res + sumI
    End Do
   End Do
  Else
#endif
   Do i1 = 1, 5, 2
    Rsf = RSup(i1:i1+1, i1:i1+1)
    Do i2=1,2
     Do i3=i2,2
      Call CoupSupZ(i2, i3, gSU2, sinW2, Rsf, coupC )
      sumI =  - 12._dp * Abs(coupC)**2 &
           &          * B22(p2, mSup2(i1 - 1 + i2), mSup2(i1 - 1 + i3) )
      If (i2.Ne.i3) sumI = 2._dp * sumI
      If (WriteOut) Write(ErrCan,*) "S-up", i1, i2, i3, sumI
      res = res + sumI
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If 
#endif
  !--------------------------
  ! Scalar down
  !--------------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1 = 1, 6
    Do i2 = 1, 6
     Call CoupSdownZ(i1, i2, gSU2, sinW2, Rsdown, coupC )
     sumI =  - 12._dp * Abs(coupC)**2 * B22(p2, mSdown2(i1), mSdown2(i2) )
     If (WriteOut) Write(ErrCan,*) "S-down", i1,i2,sumI
     res = res + sumI
    End Do
   End Do
  Else
#endif
   Do i1 = 1, 5, 2
    Rsf = RSdown(i1:i1+1, i1:i1+1)
    Do i2=1,2
     Do i3=i2,2
      Call CoupSdownZ(i2, i3, gSU2, sinW2, Rsf, coupC )
      sumI =  - 12._dp * Abs(coupC)**2 &
          &           * B22(p2, mSdown2(i1 - 1 + i2), mSdown2(i1 - 1 + i3) )
      If (i2.Ne.i3) sumI = 2._dp * sumI
      If (WriteOut) Write(ErrCan,*) "S-down", i1, i2, i3, sumI
      res = res + sumI
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If 
#endif
 !--------------------------
 ! neutrinos
 !--------------------------
  If ((n_neut.Lt.7).And.(p2.Gt.0._dp)) Then
   Call CoupFermionZ(0.5_dp, 0._dp, gSU2, sinW2, coup, coup2)
   sumI = Real(7-n_neut,dp) * coup**2  * HLoop(p2, 0._dp, 0._dp)
   If (WriteOut) Write(ErrCan,*) "Neutrinos", sumI
   res = res + sumI
  End If
  !--------------------------
  ! leptons
  !--------------------------
  Call CoupFermionZ(-0.5_dp, -1._dp, gSU2, sinW2, coup, coup2)
  Do i1 = 1,5-n_char ! number of leptons depend on model
   sumI = ( coup**2 + coup2**2 ) * HLoop(p2, mf_l2(i1), mf_l2(i1))   &
      & + 4._dp * coup * coup2 * mf_l2(i1) * B0(p2,mf_l2(i1),mf_l2(i1))
   If (WriteOut) Write(ErrCan,*) "Leptons",i1, sumI
   res = res + sumI
  End Do
  !--------------------------
  ! u-quarks
  !--------------------------
  Call CoupFermionZ(0.5_dp, 2._dp/3._dp, gSU2, sinW2, coup, coup2)
  Do i1 = 1,3
   sumI = ( coup**2 + coup2**2 ) * HLoop(p2, mf_u2(i1), mf_u2(i1))   &
      & + 4._dp * coup * coup2 * mf_u2(i1) * B0(p2,mf_u2(i1),mf_u2(i1))
   sumI = 3._dp * sumI
   If (WriteOut) Write(ErrCan,*) "u-quarks",i1, sumI
   res = res + sumI
  End Do
 !--------------------------
 ! d-quarks
 !--------------------------
  Call CoupFermionZ(-0.5_dp, -1._dp/3._dp, gSU2, sinW2, coup, coup2)
  Do i1 = 1,3
   sumI = ( coup**2 + coup2**2 ) * HLoop(p2, mf_d2(i1), mf_d2(i1))   &
      & + 4._dp * coup * coup2 * mf_d2(i1) * B0(p2,mf_d2(i1),mf_d2(i1))
   sumI = 3._dp * sumI
   If (WriteOut) Write(ErrCan,*) "d-quarks",i1, sumI
   res = res + sumI
  End Do
 !--------------------------
 ! charginos
 !--------------------------
  Do i1 = 1,n_char
   Call CoupCharginoZ(i1, i1, U, V, gSU2, cosW, coupC, coupC2)
   sumI = ( Abs(coupC)**2 + Abs(coupC2)**2 ) * HLoop(p2, mC2(i1), mC2(i1))  &
      & + 4._dp * Real( Conjg(coupC) * coupC2,dp )      &
       &        * mC2(i1) * B0(p2, mC2(i1), mC2(i1)) 
   If (WriteOut) Write(ErrCan,*) "charginos",i1,i1, sumI
   res = res + sumI
   Do i2 = i1+1,n_char
    Call CoupCharginoZ(i1, i2, U, V, gSU2, cosW, coupC, coupC2)
    sumI = 2._dp * ( Abs(coupC)**2 + Abs(coupC2)**2 )  &
       &         * HLoop(p2, mC2(i1), mC2(i2))                              &
       & + 8._dp * Real( Conjg(coupC) * coupC2,dp )      &
       &     * mC(i1) * mC(i2) * B0(p2, mC2(i1), mC2(i2)) 
     If (WriteOut) Write(ErrCan,*) "charginos",i1,i2, sumI
    res = res + sumI
   End Do
  End Do
 !--------------------------
 ! neutralino
 !--------------------------
  Do i1 = 1,n_neut
   If (Present(NMSSM)) Then
    Call CoupNeutralinoZ(i1, i1, N, gSU2, cosW, coupC, coupC2,NMSSM)
   Else
    Call CoupNeutralinoZ(i1, i1, N, gSU2, cosW, coupC, coupC2)
   End If
    sumI = 0.5_dp * ( Abs(coupC)**2 + Abs(coupC2)**2 )               &
       &          * HLoop(p2, mN2(i1), mN2(i1))                      &
       & + 2._dp * Real( Conjg(coupC) * coupC2,dp )                     &
       &         * mN2(i1) * B0(p2, mN2(i1), mN2(i1)) 
   If (WriteOut) Write(ErrCan,*) "neutralinos",i1,i1, sumI
   res = res + sumI
   Do i2 = i1+1,n_neut
    If (Present(NMSSM)) Then
     Call CoupNeutralinoZ(i1, i2, N, gSU2, cosW, coupC, coupC2, NMSSM)
    Else
     Call CoupNeutralinoZ(i1, i2, N, gSU2, cosW, coupC, coupC2)
    End If
    sumI = ( Abs(coupC)**2 + Abs(coupC2)**2 ) * HLoop(p2, mN2(i1), mN2(i2))  &
       & + 4._dp * Real( Conjg(coupC) * coupC2,dp )   &
       &         * mN(i1) * mN(i2) * B0(p2, mN2(i1), mN2(i2)) 
    If (WriteOut) Write(ErrCan,*) "neutralinos",i1,i2, sumI
    res = res + sumI
   End Do
  End Do

  res = oo16pi2 * res

  Iname = Iname - 1

 End Subroutine PiZZT1


 Subroutine PiZZT1_SM(p2, gSU2, sinW2, vev, mZ2, mW2, mS02, mf_l2, mf_u2, mf_d2 &
   &    , res)
 !-----------------------------------------------------------------------
 ! Calculates the 1-loop self energy of Z-boson. 
 ! The formula of J. Bagger et al, Nucl.Phys.B is used. The renormalization
 ! scale has to be set in the calling routine due to the structure of
 ! LoopTools. 
 ! the input is:
 !  - p2 ......... the outer momentum squared 
 !  - gSU2 ....... the SU(2) gauge coupling at p2
 !  - sinW2 ...... sin(theta_W) squared
 !  - vev ........ SM vev
 !  - mZ2 ........ Z-boson mass squared
 !  - mW2 ........ W-boson mass squared
 !  - mS02 ....... mass of the scalar Higgs squared
 ! The output:
 !  - res    
 ! written by Werner Porod, 4.1.09
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: p2, gSU2, sinW2, vev, mS02, mW2, mZ2, mf_l2(3) &
        & , mf_u2(3), mf_d2(3)
  Complex(dp), Intent(out) :: res

  Integer :: i1
  Real(dp) :: cosW2, coup, coup2
  Complex(dp) :: sumI
  Logical :: WriteOut

  Iname = Iname + 1
  NameOfUnit(Iname) = "PiZZT1_SM"

  res = ZeroC
  cosW2 = 1._dp - sinW2
  !----------------------------------
  ! In case I need the contributions
  !----------------------------------
  If ((WriteOneLoopContributions.Eq.1).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions to PiZZT1_SM:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If

  !--------------------------
  ! neutral Higgs bosons
  !--------------------------
  coup = 0.5_dp * gSU2**2 * vev / cosW2
  sumI = coup**2 * B0(p2, mZ2, mS02)
  If (WriteOut) Write(ErrCan,*) "S0", sumI
  res = res + sumI
  !--------------------------
  ! W boson
  !--------------------------
  sumI = - gSU2**2 * ( 2._dp * ( (2._dp * p2 + mW2)*cosW2 - mZ2 * sinW2**2 ) &
   &                         * B0(p2,mW2,mW2)                                &
   &                 + ( 8._dp * cosW2 + (cosW2-sinW2)**2/cosW2 )            &
   &                    * B22(p2,mW2,mW2) )
  If (WriteOut) Write(ErrCan,*) "W", sumI
  res = res + sumI
 !--------------------------
 ! neutrinos
 !--------------------------
  If ((n_neut.Lt.7).And.(p2.Gt.0._dp)) Then
   Call CoupFermionZ(0.5_dp, 0._dp, gSU2, sinW2, coup, coup2)
   sumI = Real(7-n_neut,dp) * coup**2  * HLoop(p2, 0._dp, 0._dp)
   If (WriteOut) Write(ErrCan,*) "Neutrinos", sumI
   res = res + sumI
  End If
  !--------------------------
  ! leptons
  !--------------------------
  Call CoupFermionZ(-0.5_dp, -1._dp, gSU2, sinW2, coup, coup2)
  Do i1 = 1,5-n_char ! number of leptons depend on model
   sumI = ( coup**2 + coup2**2 ) * HLoop(p2, mf_l2(i1), mf_l2(i1))   &
      & + 4._dp * coup * coup2 * mf_l2(i1) * B0(p2,mf_l2(i1),mf_l2(i1))
   If (WriteOut) Write(ErrCan,*) "Leptons",i1, sumI
   res = res + sumI
  End Do
  !--------------------------
  ! u-quarks
  !--------------------------
  Call CoupFermionZ(0.5_dp, 2._dp/3._dp, gSU2, sinW2, coup, coup2)
  Do i1 = 1,3
   sumI = ( coup**2 + coup2**2 ) * HLoop(p2, mf_u2(i1), mf_u2(i1))   &
      & + 4._dp * coup * coup2 * mf_u2(i1) * B0(p2,mf_u2(i1),mf_u2(i1))
   sumI = 3._dp * sumI
   If (WriteOut) Write(ErrCan,*) "u-quarks",i1, sumI
   res = res + sumI
  End Do
 !--------------------------
 ! d-quarks
 !--------------------------
  Call CoupFermionZ(-0.5_dp, -1._dp/3._dp, gSU2, sinW2, coup, coup2)
  Do i1 = 1,3
   sumI = ( coup**2 + coup2**2 ) * HLoop(p2, mf_d2(i1), mf_d2(i1))   &
      & + 4._dp * coup * coup2 * mf_d2(i1) * B0(p2,mf_d2(i1),mf_d2(i1))
   sumI = 3._dp * sumI
   If (WriteOut) Write(ErrCan,*) "d-quarks",i1, sumI
   res = res + sumI
  End Do

  res = oo16pi2 * res

  Iname = Iname - 1

 End Subroutine PiZZT1_SM


#ifdef GENERATIONMIXING
 Subroutine ScalarMass_Loop_MSSM(mS02, mpole2, RS0                            &
      & , Q2, tanbQ, vevs_1L, mZ2_mZ, mW2, mA2_mA, b_a                        &
      & , gU1, gSU2, gSU3, Y_l, Y_d, Y_u, uD_L, uD_R, uL_L, uL_R, uU_L, uU_R  &
      & , mUSquark2, RUSquark, mDSquark2, RDSquark, mSneutrino2, RSneutrino   &
      & , mSlepton2, RSlepton, A_d, A_u, A_l, mu, vevs_DR, mP02, RP0, mSpm2   &
      & , RSpm, mC, mC2, U, V, mN, mN2, N, M_Q2, M_U2, M_D2, M_E2, M_L2, mglu &
      & , i_os, mass, mass2, R, i_run, kont)
#else
 Subroutine ScalarMass_Loop_MSSM(mS02, RS0, Q2, tanbQ, vevs_1L, mZ2_mZ        &
      & , mA2_mA, b_a, gU1, gSU2, gSU3, Y_l, Y_d, Y_u, mUSquark2, RUSquark    &
      & , mDSquark2, RDSquark, mSneutrino2, mSlepton2, RSlepton, A_d, A_u     &
      & , A_l, mu, vevs_DR, mP02, RP0, mSpm2, RSpm, mC, mC2, U, V, mN, mN2, N &
      & , M_Q2, M_U2, M_D2, M_E2, M_L2, mglu, i_os, mass, mass2, R, i_run, kont)
#endif
 !-----------------------------------------------------------------
 ! calculates scalar masses + mixing matrix R in the MSSM,
 ! 1-generation epsilon model and 3-generation epsilon model. 
 ! input:
 ! output:
 !  mass(i) .... pseudoscalar masses
 !  mass2(i) ... = mass(i)**2
 !  R(i,j) .. mixing matrix of Scalars 
 ! written by Werner Porod, 1.12.99
 ! 10.10.01: portation to f90, new main strategy: calculate all couplings
 !           in this routine and pass to PiScalar only general couplings
 !           In this way one can use PiScalar for different models 
 ! 20.02.02: - correcting a bug in the definition of p2
 !           - implementing 2-loop corrections using the modules
 !             BDSZHiggs.f  DSZHiggs.f  funcs.f provided by Pietro Slavich
 ! 09.07.02: - implementing 2-loop corrections using the module BDSZodd.f 
 !             provided by Pietro Slavich
 ! 11.06.03: implementing alpha_b alpha_t + alpha^2_b corrections
 !           based on A. Dedes, G. Degrassi and P. Slavich, hep-ph/0305127
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) ::  i_os, i_run
  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: gU1, gSU2, vevs_DR(2), mA2_mA, b_a, vevs_1L(2) &
     & , mS02(2), RS0(2,2), mUSquark2(6), mDSquark2(6), mSneutrino2(3)      &
     & , mSlepton2(6), mP02(2), RP0(2,2), mSpm2(2), mC(2), mC2(2), mN(4)    &
     & , mN2(4), M_Q2, M_U2, M_D2, mglu, Q2, gSU3, mW2, mZ2_mZ, tanbQ      &
     & , mpole2(2), M_E2, M_L2
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3), A_u(3,3)         &
     & , RUSquark(6,6), RDSquark(6,6), RSlepton(6,6), A_l(3,3), A_d(3,3)    &
     & , mu, RSpm(2,2), U(2,2), V(2,2), N(4,4)
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: uD_L(3,3), uD_R(3,3), uL_L(3,3), uL_R(3,3)    &
     & , uU_L(3,3), uU_R(3,3), RSneutrino(3,3)
#endif
  Real(dp), Intent(inout) :: mass(2), mass2(2), R(2,2)

  Integer :: i1, i2, i3, i4
  Real(dp) :: mat2(2,2), cosb2, sinb2, sinbcosb, p2, e_d, e_u, sinW2_DR     &
     & , cosW_DR, sinb, cosb, c_S0WW(2), c_P0S0Z(2,2), c_S0ZZ(2)            &
     & , c_S0S0WW(2), c_S0S0ZZ(2), c_P0P0S0S0(2,2,2), c_P0S0S0(2,2,2)       &
     & , c_S03(2,2,2), c_S04(2,2,2), Mi2(2), mat2a(2,2), Pi2S(2,2), nen     &
     & , det, sumI, nen2
  Logical :: WriteOut
  Complex(dp) :: PiS1(2,2), c_DDS0_L(3), c_DDS0_R(3), c_LLS0_L(3), RSd(2,2) &
     & , RSu(2,2), RSl(2,2), c_LLS0_R(3), c_UUS0_L(3), c_UUS0_R(3)          &
     & , c_SdSdS0S0(6,2), c_SuSuS0S0(6,2), c_SlSlS0S0(6,2), c_SnSnS0S0(3,2) &
     & , c_SdSdS0(6,6,2), c_SuSuS0(6,6,2), c_SlSlS0(6,6,2), c_SnSnS0(3,3,2) &
     & , Zero3(3,3), c_SpmS0W(2,2), c_SpSmS0S0(2,2,2), c_SpSmS0(2,2,2)      &
     & , c_CCS0_L(2,2,2), c_CCS0_R(2,2,2), c_NNS0_L(4,4,2), c_NNS0_R(4,4,2)

  Integer :: i_lh
#ifdef THREELOOP
  !----------------------------
  ! for 3-loop part
  !----------------------------
  Real(dp) :: mT, mSt(2), mSq, Pi3S(2,2), Rst2(2)
  Complex(dp) :: Rstop(2,2)
#endif

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ScalarMass_Loop_MSSM'

  !----------------------------------
  ! In case I need the contributions
  !----------------------------------
  If ((WriteOneLoopContributions.Eq.6).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions to ScalarMass_Loop_MSSM:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If

 !-----------------------
 ! initialisation
 !-----------------------
  cosb2 = 1._dp / (1._dp + tanbQ**2) 
  sinb2 = 1._dp - cosb2
  sinbcosb = Sqrt(cosb2*sinb2)
  sinb = Sqrt( sinb2 )
  cosb = Sqrt( cosb2 )

  mat2a(1,1) = mZ2_mZ * cosb2 + (mA2_mA  - b_A) * sinb2 + vevs_1L(1)
  mat2a(1,2) = -(mZ2_mZ + mA2_mA - b_A) * sinbcosb
  mat2a(2,1) = mat2a(1,2)
  mat2a(2,2) = mZ2_mZ * sinb2 +( mA2_mA - b_A) * cosb2 + vevs_1L(2)

  If (WriteOut) Then 
   Write(ErrCan,*) "Tree level + pseudoscalar part:", mat2a(1,1), mat2a(1,2)
   Write(ErrCan,*) "                               ", mat2a(2,1), mat2a(2,2)
   Write(ErrCan,*) " The contribution von PiScalar:"
  End If

  !------------------------------------------------------------
  ! Calculate first the couplings needed for the loop  
  !------------------------------------------------------------
  e_d = -1._dp / 3._dp
  e_u = 2._dp / 3._dp

  sinW2_DR = gU1**2 / (gSU2**2 + gU1**2)
  cosW_DR = Sqrt( 1._dp - sinW2_DR )

  Zero3 = 0._dp

#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,3
    Call CoupFermionScalar3(i1, i1, 1, -0.5_dp, Y_d, uD_L, uD_R, id2R &
                                   &, c_DDS0_L(i1), c_DDS0_R(i1) )
    Call CoupFermionScalar3(i1, i1, 1, -0.5_dp, Y_l, uL_L, uL_R, id2R &
                                   &, c_LLS0_L(i1), c_LLS0_R(i1) )
    Call CoupFermionScalar3(i1, i1, 2, 0.5_dp, Y_u, uU_L, uU_R, id2R &
                                   &, c_UUS0_L(i1), c_UUS0_R(i1) )
    Do i2=1,2
     Call CoupScalarSfermion4(i2, i2, i1, i1, id2R, 0.5_dp, 0._dp, Zero3 &
                            &, RSneutrino, gU1, gSU2, c_SnSnS0S0(i1,i2) )
    End Do
    Do i2=1,3
     Do i3=1,2
      Call CoupScalarSfermion3(i3, i1, i2, id2R, 0.5_dp , 0._dp, Zero3  &
         & , Rsneutrino, Zero3, mu, vevs_DR, gU1, gSU2, c_SnSnS0(i1,i2,i3) )
     End Do
    End Do
   End Do

   Do i1=1,6
    Do i2=1,2
     Call CoupScalarSfermion4(i2, i2, i1, i1, id2R, -0.5_dp, e_d, Y_d &
                             &, RDSquark, gU1, gSU2, c_SdSdS0S0(i1,i2) )
     Call CoupScalarSfermion4(i2, i2, i1, i1, id2R, -0.5_dp, -1._dp &
                   &, Y_l, RSlepton, gU1, gSU2, c_SlSlS0S0(i1,i2) )
     Call CoupScalarSfermion4(i2, i2, i1, i1, id2R, 0.5_dp, e_u, Y_u &
                             &, RUSquark, gU1, gSU2, c_SuSuS0S0(i1,i2) )
    End Do
    Do i2=1,6
     Do i3=1,2
      Call CoupScalarSfermion3(i3, i1, i2, id2R, -0.5_dp , e_d, Y_d  &
         & , RDSquark, A_d, mu, vevs_DR, gU1, gSU2, c_SdSdS0(i1,i2,i3) )
      Call CoupScalarSfermion3(i3, i1, i2, id2R, -0.5_dp , -1._dp, Y_l  &
         & , RSlepton, A_l, mu, vevs_DR, gU1, gSU2, c_SlSlS0(i1,i2,i3) )
      Call CoupScalarSfermion3(i3, i1, i2, id2R, 0.5_dp , e_u, Y_u  &
         & , RUSquark, A_u, mu, vevs_DR, gU1, gSU2, c_SuSuS0(i1,i2,i3) )
     End Do
    End Do
   End Do
    
  Else ! .not.GenerationMixing
#endif
   Do i1=1,3
    Call CoupFermionScalar(1, -0.5_dp, Y_d(i1,i1), id2R &
                         &, c_DDS0_L(i1), c_DDS0_R(i1))
    Call CoupFermionScalar(1, -0.5_dp, Y_l(i1,i1), id2R &
                         &, c_LLS0_L(i1), c_LLS0_R(i1))
    Call CoupFermionScalar(2, 0.5_dp, Y_u(i1,i1), id2R &
                         &, c_UUS0_L(i1), c_UUS0_R(i1))
    i2 = 2*i1 - 1
    Rsd = RDSquark(i2:i2+1,i2:i2+1)
    Rsu = RUSquark(i2:i2+1,i2:i2+1)
    Rsl = Rslepton(i2:i2+1,i2:i2+1)
    Do i2=1,2
     Do i3=1,2
      Call CoupScalarSfermion4(i3, i3, i2, i2, id2R, -0.5_dp, e_d, Y_d(i1,i1) &
                             &, RSd, gU1, gSU2, c_SdSdS0S0((i1-1)*2+i2,i3) )
      Call CoupScalarSfermion4(i3, i3, i2, i2, id2R, -0.5_dp, -1._dp &
                   &, Y_l(i1,i1), RSl, gU1, gSU2, c_SlSlS0S0((i1-1)*2+i2,i3) )
      Call CoupScalarSfermion4(i3, i3, i2, i2, id2R, 0.5_dp, e_u, Y_u(i1,i1) &
                             &, RSu, gU1, gSU2, c_SuSuS0S0((i1-1)*2+i2,i3) )
      Do i4=1,2
       Call CoupScalarSfermion3(i2, i3, i4, id2R, -0.5_dp, e_d, Y_d(i1,i1)   &
                               &, Rsd, A_d(i1,i1), mu, vevs_DR, gU1, gSU2    &
                               &, c_SdSdS0((i1-1)*2+i3,(i1-1)*2+i4,i2) )
       Call CoupScalarSfermion3(i2, i3, i4, id2R, -0.5_dp, -1._dp, Y_l(i1,i1) &
                               &, Rsl, A_l(i1,i1), mu, vevs_DR, gU1, gSU2     &
                               &, c_SlSlS0((i1-1)*2+i3,(i1-1)*2+i4,i2) )
       Call CoupScalarSfermion3(i2, i3, i4, id2R, 0.5_dp, e_u, Y_u(i1,i1)    &
                               &, Rsu, A_u(i1,i1), mu, vevs_DR, gU1, gSU2    &
                               &, c_SuSuS0((i1-1)*2+i3,(i1-1)*2+i4,i2) )
      End Do
     End Do
     Call CoupScalarSfermion4(i2, i2, 1, 1, id2R, 0.5_dp, 0._dp, ZeroC &
                            &, id2C, gU1, gSU2, c_SnSnS0S0(i1,i2) )
     Call CoupScalarSfermion3(i2, 1, 1, id2R, 0.5_dp, 0._dp, ZeroC    &
                             &, id2C, ZeroC, mu, vevs_DR, gU1, gSU2     &
                             &, c_SnSnS0(i1,i1,i2) )
    End Do

   End Do ! i1

#ifdef GENERATIONMIXING
  End If
#endif

  c_SpmS0W(2,1) = gSU2 * oosqrt2 * sinb
  c_SpmS0W(2,2) = - gSU2 * oosqrt2 * cosb
  c_S0WW(1) = c_SpmS0W(2,2)
  c_S0WW(2) = - c_SpmS0W(2,1)

  c_P0S0Z(2,1) = gSU2 * 0.5_dp * sinb / cosW_DR
  c_P0S0Z(2,2) = - gSU2 * 0.5_dp * cosb / cosW_DR
  c_S0ZZ(1) = c_P0S0Z(2,2)
  c_S0ZZ(2) = - c_P0S0Z(2,1)

  c_S0S0WW = 2._dp * gSU2**2
  c_S0S0ZZ = gSU2**2 + gU1**2

  Do i1=1,2
   Do i2=1,2
    Do i3=1,2
     Call CoupChargedScalarScalar4(i1,i1, i2, i3, RSpm, id2R, gU1, gSU2  &
                                 &, c_SpSmS0S0(i1,i2,i3) ) 
     Call CoupChargedScalarScalar3(i1, i2, i3, RSpm, id2R, vevs_DR, gU1 &
                                &, gSU2, c_SpSmS0(i1,i2,i3) )

     Call CoupPseudoScalarScalar3(i1, i2, i3, RP0, id2R, gU1, gSU2, vevs_DR &
                                &, c_P0S0S0(i1,i2,i3))
     Call CoupPseudoScalarScalar4(i1, i1, i2, i3, RP0, id2R, gU1, gSU2 &
                                 &, c_P0P0S0S0(i1,i2,i3) )

     Call CoupScalar4a(i2,i3,i1,i1, RS0, gU1, gSU2, c_S04(i1,i2,i3) )
     Call CoupScalar3a(i1,i2,i3, RS0, gU1, gSU2, vevs_DR, c_S03(i1,i2,i3) )

     Call CoupCharginoScalar(i1, i2, i3, U, V, id2R, gSU2  &
                            &, c_CCS0_L(i1,i2,i3), c_CCS0_R(i1,i2,i3))
    End Do
   End Do
   Do i2=1,4
    Do i3=1,4
     Call CoupNeutralinoScalar(i2, i3, i1, N, id2R, gU1, gSU2  &
                            &, c_NNS0_L(i2,i3,i1), c_NNS0_R(i2,i3,i1))
    End Do
   End Do
  End Do

  !---------------------------------------------------------------
  ! 2-loop part does not depend on momentum, therefore done first
  !---------------------------------------------------------------
  If (Only_1loop_Higgsmass) Then
   Pi2s = 0._dp
  Else
   Call PiScalar2(Q2, gSU3, mglu, mP02(2), vevs_DR, m_D2, m_U2, m_Q2, m_E2  &
      & , m_L2, A_d(3,3), A_u(3,3), A_l(3,3), Y_d(3,3), Y_u(3,3), Y_l(3,3)  &
      & , mu, i_os, Pi2s, kont)
#ifdef THREELOOP
   If (ThreeLoopHiggsMass) Then
    mT = oosqrt2 * Real(Y_u(3,3)) * vevs_DR(2)
    If (GenerationMixing) Then
     mSq = Sum(Sqrt(mDsquark2))
     i_lh = 1
     Do i1=1,6
      Rst2(i_lh) = Abs(RUsquark(i1,3))**2+Abs(RUsquark(i1,6))**2
      If (RSt2(i_lh).Gt.0.7_dp) Then
       mSt(i_lh) = Sqrt(mUsquark2(i1))
       If (i_lh.Eq.1) Then
        RStop(1,1) = RUsquark(i1,3)
        RStop(1,2) = RUsquark(i1,6)
        ! to ensure unitarity
        RStop(1,:) = RStop(1,:) / Rst2(i_lh)
        RStop(2,2) = RStop(1,1)
        RStop(2,1) = -RStop(1,2)
       ! in case the stop content of the 2nd state is larger:
       Else If ((i_lh.Eq.2).And.(Rst2(2).Gt.Rst2(1))) Then
        RStop(2,1) = RUsquark(i1,3)
        RStop(2,2) = RUsquark(i1,6)
        RStop(2,:) = RStop(2,:) / Rst2(i_lh)
        RStop(1,1) = RStop(2,2)
        RStop(1,2) = -RStop(2,1)
       End If
       i_lh = i_lh + 1
      Else
       mSq = mSq + Sqrt(mUsquark2(i1))
      End If
     End Do
     If (i_lh.Ne.3) Then
      Write(ErrCan,*) "Problem for the calculation of 3-loop contribution to m_h"
      Call TerminateProgram()
     End If
     mSq = mSq / 10._dp
    Else
     mSt = Sqrt(mUsquark2(5:6))
     Rstop = RUsquark(5:6,5:6)
     mSq = (Sum(Sqrt(mUsquark2(1:4)))+Sum(Sqrt(mDsquark2))) / 10._dp
    End If
!msq=Sqrt(mUsquark2(5))
! Write(*,*) mSt,mSq,mglu,Sqrt(q2),gsu3
    Call PiScalar3(Q2, gSU3, mT, mSt, Rstop, mSq, mglu, mu &
                & , vevs_DR(2)/vevs_DR(1), Pi3s)
!Write(*,*) "Pi2",pi2s(1,:)
!Write(*,*) "   ",pi2s(2,:)
!Write(*,*) "Pi3",pi3s(1,:)
!Write(*,*) "   ",pi3s(2,:)
    Pi2s=Pi2s+Pi3s
   End If
#endif
  End If

  If ((i_run.Eq.1).And.(kont.Ne.0)) Then
   kont = 0
   Pi2s = 0._dp
  End If
  !------------------------------------------------------------
  ! now the loop calculation + masses and mixing angle
  !------------------------------------------------------------
  Do i_lh=2,1,-1
   p2 = mpole2(i_lh)
   Call PiScalar(p2, c_DDS0_L, c_UUS0_L, c_LLS0_L                             &
     & , mDSquark2, c_SdSdS0S0, mUSquark2, c_SuSuS0S0, mSlepton2, c_SlSlS0S0  &
     & , mSneutrino2, c_SnSnS0S0, c_SdSdS0, c_SuSuS0, c_SlSlS0, c_SnSnS0      &
     & , mS02, mP02, mSpm2, c_SpmS0W, c_S0WW, c_P0S0Z, c_S0ZZ, c_S0S0WW       &
     & , c_S0S0ZZ, c_SpSmS0S0, c_SpSmS0, mC, mC2, c_CCS0_L, c_CCS0_R          &
     & , c_P0P0S0S0, c_P0S0S0, c_S03, c_S04, mN, mN2, c_NNS0_L, c_NNS0_R      &
     & , mZ2_mZ, mW2, WriteOut, PiS1)

   mat2 = mat2a - Real(PiS1,dp) + Pi2s

   det = 0.5_dp * Sqrt((mat2(1,1)-mat2(2,2))**2+4*Abs(mat2(1,2))**2)
   sumI = 0.5*(mat2(1,1)+mat2(2,2))
   mi2(1) = sumI - det
   mi2(2) = sumI + det
   nen2 = (mat2(1,1) - mi2(1) )**2 + mat2(1,2)**2
   nen = Sqrt(nen2)
   If (nen.Lt.1.e-7_dp) Then
    R(1,1) = 1._dp
    R(1,2) = 0._dp
   Else
    R(1,1) = - mat2(1,2) / nen
    R(1,2) = ( mat2(1,1) - mi2(1) ) / nen
   End If
   R(2,1) = - R(1,2)
   R(2,2) = R(1,1)

    If (Mi2(i_lh).Lt.0._dp) Then
     WriteWert = Mi2(i_lh)
     Call WriteLoopMassesError(-9, "ScalarMass_Loop_MSSM", kont)
     Iname = Iname - 1
     Return
    Else
     mass2(i_lh) = Mi2(i_lh)
     mass(i_lh) = Sqrt( mass2(i_lh) )
    End If
  End Do ! i_lh


  Iname = Iname - 1

 End Subroutine ScalarMass_Loop_MSSM


 Subroutine SetLoopMassModel(i_c,i_n,i_p,i_s,i_spm,i_sl,i_sn)
 !----------------------------------------------------------------------
 ! this is needed if one accesses the various subroutines from outside
 ! written by Werner Porod, 14.01.02
 !----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_c, i_n, i_p, i_s, i_spm, i_sl, i_sn
  n_S0 = i_s
  n_P0 = i_p
  n_Spm = i_spm
  n_char = i_c
  n_neut = i_n
  n_Sneut = i_sn
  n_Slept = i_sl  
 End Subroutine SetLoopMassModel


 Subroutine SetWriteOneLoopContributions(n)
 Implicit None
  !-------------------------------------------------------------------
  ! This subroutine sets the flag controlling whether the contributions
  ! to the different 1-loop quantities are written out and if yes
  ! which ones
  ! written by Werner Porod, 09.10.01
  ! The various numbers are:
  !   n <  0 ................. everything is written
  !   n =  1 ................. PiZZT1
  !   n =  2 ................. One_Loop_Tadpoles_MSSM
  !   n =  3 ................. PiPseudoScalar
  !   n =  4 ................. PiChargedScalar
  !   n =  5 ................. PiWWT1
  !   n =  6 ................. PiScalar + Scalar_Loop
  !   n =  7 ................. CharginoMass_Loop + Sigma_Chargino
  !   n =  8 ................. Sigma_Gluino
  !   n =  9 ................. delta_VB
  !   n = 10 ................. Sigma_Fermion
  !   n = 11 ................. NeutralinoMass_Loop + Sigma_Neutralino
  !   n = 12 ................. SquarkMass_1L + PiSquark
  !   n = 13 ................. SleptonMass_1L + PiSlepton
  !   n = 13 ................. SneutrinoMass_1L + PiSneutrino
  !-------------------------------------------------------------------
  Integer, Intent(in) :: n

  WriteOneLoopContributions = n
 End Subroutine SetWriteOneLoopContributions


 Subroutine Sigma_Chargino(p2, mC, mC2, mSdown2, c_CUSd_L, c_CUSd_R        &
       & , mSup2, c_CDSu_L, c_CDSu_R, mSlepton2, c_CNuSl_L, c_CNuSl_R      &
       & , mSneut2, c_CLSn_L, c_CLSn_R, mN, mN2, c_CNW_L, c_CNW_R          &
       & , mSpm2, c_SmpCN_L, c_SmpCN_R, c_CCZ_L, c_CCZ_R, c_CCG_L, c_CCG_R &
       & , mP02, c_CCP0_L, c_CCP0_R, mS02, c_CCS0_L, c_CCS0_R, mZ2, mW2    &
       & , WriteOut, SigL, SigR, SigS)
 !-----------------------------------------------------------------------
 ! written by Werner Porod, 14.10.1999
 ! 07.11.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: p2, mC(:), mC2(:), mSdown2(:), mSup2(:), mZ2, mW2   &
       & , mSlepton2(:), mSneut2(:), mN(:), mN2(:), mSpm2(:), mP02(:), mS02(:)
  Complex(dp), Intent(in) :: c_CUSd_L(:,:,:), c_CUSd_R(:,:,:)                 &
     & , c_CDSu_L(:,:,:), c_CDSu_R(:,:,:), c_CNuSl_L(:,:,:), c_CNuSl_R(:,:,:) &
     & , c_CLSn_L(:,:,:), c_CLSn_R(:,:,:), c_CNW_L(:,:), c_CNW_R(:,:)         &
     & , c_SmpCN_L(:,:,:), c_SmpCN_R(:,:,:), c_CCZ_L(:,:), c_CCZ_R(:,:)       &
     & , c_CCG_L(:,:), c_CCG_R(:,:), c_CCP0_L(:,:,:), c_CCP0_R(:,:,:)         &
     & , c_CCS0_L(:,:,:), c_CCS0_R(:,:,:)
  Logical, Intent(in) :: WriteOut
  Complex(dp), Intent(out) :: SigL(:,:), SigR(:,:), sigS(:,:)

  Real(dp) :: B0m2, B1m2, B0m2a, B1m2a
  Integer :: i1, i2, i3, i4, i_sf
  Complex(dp) :: sumL(n_char,n_char), sumR(n_char,n_char), sumS(n_char,n_char)

 !------------------
 ! Inititalisation
 !------------------
  SigL = ZeroC
  SigR = ZeroC
  SigS = ZeroC

 !---------------------
 ! quark - squark
 !---------------------
  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,6
     sumL = 0._dp
     sumR = 0._dp
     sumS = 0._dp
     B1m2 = - 1.5_dp * Real( B1(p2, mf_U2(i1), mSdown2(i2) ),dp )
     B0m2 = 3._dp * mf_u(i1) * Real( B0(p2, mf_U2(i1), mSdown2(i2)),dp )
     Do i3=1,n_char
      Do i4=1,n_char
       sumL(i3,i4) = Conjg( c_CUSd_L(i4,i1,i2) ) * c_CUSd_L(i3,i1,i2) * B1m2
       sumR(i3,i4) = Conjg( c_CUSd_R(i4,i1,i2) ) * c_CUSd_R(i3,i1,i2) * B1m2
       sumS(i3,i4) = Conjg( c_CUSd_L(i4,i1,i2) ) * c_CUSd_R(i3,i1,i2)* B0m2
      End Do
      If (WriteOut) Then
       Write (ErrCan,*) "C U Sd",i3,i1,i2
       Write (ErrCan,*) "  sig_L",sumL(i3,:)
       Write (ErrCan,*) "  sig_R",sumR(i3,:)
       Write (ErrCan,*) "  sig_S",sumS(i3,:)
      End If
     End Do
     sigL = sigL + sumL
     sigR = sigR + sumR
     sigS = sigS + sumS

     sumL = 0._dp
     sumR = 0._dp
     sumS = 0._dp
     B1m2 = - 1.5_dp * Real( B1(p2, mf_D2(i1), mSup2(i2) ),dp )
     B0m2 = 3._dp * mf_d(i1) * Real( B0(p2, mf_D2(i1), mSup2(i2)),dp )
     Do i3=1,n_char
      Do i4=1,n_char
       sumL(i3,i4) = Conjg( c_CDSu_R(i4,i1,i2) ) * c_CDSu_R(i3,i1,i2) * B1m2
       sumR(i3,i4) = Conjg( c_CDSu_L(i4,i1,i2) ) * c_CDSu_L(i3,i1,i2) * B1m2
       sumS(i3,i4) = Conjg( c_CDSu_R(i4,i1,i2) ) * c_CDSu_L(i3,i1,i2)* B0m2
      End Do
      If (WriteOut) Then
       Write (ErrCan,*) "C D Su",i3,i1,i2
       Write (ErrCan,*) "  sig_L",sumL(i3,:)
       Write (ErrCan,*) "  sig_R",sumR(i3,:)
       Write (ErrCan,*) "  sig_S",sumS(i3,:)
      End If
     End Do
     sigL = sigL + sumL
     sigR = sigR + sumR
     sigS = sigS + sumS

    End Do
   End Do

  Else

   Do i1=1,3
    Do i_sf=1,2
     i2 = (i1-1)*2 + i_sf
     sumL = 0._dp
     sumR = 0._dp
     sumS = 0._dp
     B1m2 = - 1.5_dp * Real( B1(p2, mf_U2(i1), mSdown2(i2) ),dp )
     B0m2 = 3._dp * mf_u(i1) * Real( B0(p2, mf_U2(i1), mSdown2(i2)),dp )
     Do i3=1,n_char
      Do i4=1,n_char
       sumL(i3,i4) = Conjg( c_CUSd_L(i4,i1,i2) ) * c_CUSd_L(i3,i1,i2) * B1m2
       sumR(i3,i4) = Conjg( c_CUSd_R(i4,i1,i2) ) * c_CUSd_R(i3,i1,i2) * B1m2
       sumS(i3,i4) = Conjg( c_CUSd_L(i4,i1,i2) ) * c_CUSd_R(i3,i1,i2)* B0m2
      End Do
      If (WriteOut) Then
       Write (ErrCan,*) "C U Sd",i3,i1,i2
       Write (ErrCan,*) "  sig_L",sumL(i3,:)
       Write (ErrCan,*) "  sig_R",sumR(i3,:)
       Write (ErrCan,*) "  sig_S",sumS(i3,:)
      End If
     End Do
     sigL = sigL + sumL
     sigR = sigR + sumR
     sigS = sigS + sumS

     sumL = 0._dp
     sumR = 0._dp
     sumS = 0._dp
     B1m2 = - 1.5_dp * Real( B1(p2, mf_D2(i1), mSup2(i2) ),dp )
     B0m2 = 3._dp * mf_d(i1) * Real( B0(p2, mf_D2(i1), mSup2(i2)),dp )
     Do i3=1,n_char
      Do i4=1,n_char
       sumL(i3,i4) = Conjg( c_CDSu_R(i4,i1,i2) ) * c_CDSu_R(i3,i1,i2) * B1m2
       sumR(i3,i4) = Conjg( c_CDSu_L(i4,i1,i2) ) * c_CDSu_L(i3,i1,i2) * B1m2
       sumS(i3,i4) = Conjg( c_CDSu_R(i4,i1,i2) ) * c_CDSu_L(i3,i1,i2) * B0m2
      End Do
      If (WriteOut) Then
       Write (ErrCan,*) "C D Su",i3,i1,i2
       Write (ErrCan,*) "  sig_L",sumL(i3,:)
       Write (ErrCan,*) "  sig_R",sumR(i3,:)
       Write (ErrCan,*) "  sig_S",sumS(i3,:)
      End If
     End Do
     sigL = sigL + sumL
     sigR = sigR + sumR
     sigS = sigS + sumS

    End Do
   End Do

  End If
 !---------------------
 ! lepton - slepton
 !---------------------
  If (GenerationMixing) Then
   Do i1=1,5-n_char
    Do i2=1,2*(5-n_char)
     sumL = 0._dp
     sumR = 0._dp
     B1m2 = - 0.5_dp * Real( B1(p2, 0._dp, mSlepton2(i2) ) ,dp)
     Do i3=1,n_char
      Do i4=1,n_char
       sumR(i3,i4) = Conjg( c_CNuSl_R(i4,i1,i2) ) * c_CNuSl_R(i3,i1,i2) * B1m2
      End Do
      If (WriteOut) Then
       Write (ErrCan,*) "C Nu Sl",i3,i1,i2
       Write (ErrCan,*) "  sig_R",sumR(i3,:)
      End If
     End Do
     sigR = sigR + sumR
    End Do

    Do i2=1,5-n_char
     sumL = 0._dp
     sumR = 0._dp
     sumS = 0._dp
     B1m2 = - 0.5_dp * Real( B1(p2, mf_L2(i1), mSneut2(i2) ),dp )
     B0m2 = mf_l(i1) * Real( B0(p2, mf_L2(i1), mSneut2(i2)),dp )
     Do i3=1,n_char
      Do i4=1,n_char
       sumL(i3,i4) = Conjg( c_CLSn_R(i4,i1,i2) ) * c_CLSn_R(i3,i1,i2) * B1m2
       sumR(i3,i4) = Conjg( c_CLSn_L(i4,i1,i2) ) * c_CLSn_L(i3,i1,i2) * B1m2
       sumS(i3,i4) = Conjg( c_CLSn_R(i4,i1,i2) ) * c_CLSn_L(i3,i1,i2)* B0m2
      End Do
      If (WriteOut) Then
       Write (ErrCan,*) "C L Sn",i3,i1,i2
       Write (ErrCan,*) "  sig_L",sumL(i3,:)
       Write (ErrCan,*) "  sig_R",sumR(i3,:)
       Write (ErrCan,*) "  sig_S",sumS(i3,:)
      End If
     End Do
     sigL = sigL + sumL
     sigR = sigR + sumR
     sigS = sigS + sumS
    End Do
   End Do

  Else

   Do i1=1,5 - n_char
    Do i_sf=1,2
     i2 = (i1-1)*2 + i_sf

     sumL = 0._dp
     sumR = 0._dp
     B1m2 = - 0.5_dp * Real( B1(p2, 0._dp, mSlepton2(i2) ),dp )
     Do i3=1,n_char
      Do i4=1,n_char
       sumR(i3,i4) = Conjg( c_CNuSl_R(i4,i1,i2) ) * c_CNuSl_R(i3,i1,i2) * B1m2
      End Do
      If (WriteOut) Then
       Write (ErrCan,*) "C Nu Sl",i3,i1,i2
       Write (ErrCan,*) "  sig_L",sumL(i3,:)
       Write (ErrCan,*) "  sig_R",sumR(i3,:)
      End If
     End Do
     sigR = sigR + sumR
    End Do

    i2=i1
    sumL = 0._dp
    sumR = 0._dp
    sumS = 0._dp
    B1m2 = - 0.5_dp * Real( B1(p2, mf_L2(i1), mSneut2(i2) ),dp )
    B0m2 = mf_l(i1) * Real( B0(p2, mf_L2(i1), mSneut2(i2)),dp )
    Do i3=1,n_char
     Do i4=1,n_char
      sumL(i3,i4) = Conjg( c_CLSn_R(i4,i1,i2) ) * c_CLSn_R(i3,i1,i2) * B1m2
      sumR(i3,i4) = Conjg( c_CLSn_L(i4,i1,i2) ) * c_CLSn_L(i3,i1,i2) * B1m2
      sumS(i3,i4) = Conjg( c_CLSn_R(i4,i1,i2) ) * c_CLSn_L(i3,i1,i2) * B0m2
     End Do
     If (WriteOut) Then
      Write (ErrCan,*) "C L Sn",i3,i1,i2
      Write (ErrCan,*) "  sig_L",sumL(i3,:)
      Write (ErrCan,*) "  sig_R",sumR(i3,:)
      Write (ErrCan,*) "  sig_S",sumS(i3,:)
     End If
    End Do
    sigL = sigL + sumL
    sigR = sigR + sumR
    sigS = sigS + sumS
   End Do

  End If

 !-------------------------------------
 ! neutralino - W
 !-------------------------------------
  Do i1 = 1,n_neut
   B1m2 = - Real( B1(p2, mN2(i1), mW2),dp )
   B0m2 = - 4._dp * mN(i1) * Real( B0(p2, mN2(i1), mW2),dp )
   sumL = 0._dp
   sumR = 0._dp
   sumS = 0._dp
   Do i2 = 1,n_char
    Do i3 = 1,n_char
     sumL(i2,i3) = Conjg( c_CNW_R(i3,i1) ) * c_CNW_R(i2,i1) * B1m2
     sumR(i2,i3) = Conjg( c_CNW_L(i3,i1) ) * c_CNW_L(i2,i1) * B1m2
     sumS(i2,i3) = Conjg( c_CNW_R(i3,i1) ) * c_CNW_L(i2,i1) * B0m2     
    End Do 
    If (WriteOut) Then
     Write (ErrCan,*) "C N W",i2,i1
     Write (ErrCan,*) "  sig_L",sumL(i2,:)
     Write (ErrCan,*) "  sig_R",sumR(i2,:)
     Write (ErrCan,*) "  sig_S",sumS(i2,:)
    End If
   End Do
   sigL = sigL + sumL
   sigR = sigR + sumR
   sigS = sigS + sumS
  End Do
 !-------------------------------------
 ! neutralino - S+
 !-------------------------------------
  Do i1 = 1,n_neut
   Do i2 = 1,n_Spm
    B1m2 = - 0.5_dp * Real( B1(p2, mN2(i1),  mSpm2(i2)),dp )
    B0m2 = mN(i1) * Real( B0(p2, mN2(i1),  mSpm2(i2)),dp )
    sumL = 0._dp
    sumR = 0._dp
    sumS = 0._dp
    Do i3 = 1,n_char
     Do i4 = 1,n_char
      sumL(i3,i4) = Conjg( c_SmpCN_L(i2,i4,i1) ) * c_SmpCN_L(i2,i3,i1) * B1m2
      sumR(i3,i4) = Conjg( c_SmpCN_R(i2,i4,i1) ) * c_SmpCN_R(i2,i3,i1) * B1m2
      sumS(i3,i4) = Conjg( c_SmpCN_L(i2,i4,i1) ) * c_SmpCN_R(i2,i3,i1) * B0m2
     End Do 
     If (WriteOut) Then
      Write (ErrCan,*) "S^+ C N",i2,i3,i1
      Write (ErrCan,*) "  sig_L",sumL(i3,:)
      Write (ErrCan,*) "  sig_R",sumR(i3,:)
      Write (ErrCan,*) "  sig_S",sumS(i3,:)
     End If
    End Do 
    sigL = sigL + sumL
    sigR = sigR + sumR
    sigS = sigS + sumS
   End Do
  End Do
 !-------------------------------------
 ! neutral gauge bosons
 !-------------------------------------
  Do i1 = 1,n_char
   B1m2 = - Real( B1(p2, mC2(i1), mZ2),dp )
   B0m2 = - 4._dp * mC(i1) * Real( B0(p2, mC2(i1), mZ2),dp )
   B1m2a = - Real( B1(p2, mC2(i1), 0._dp),dp )
   B0m2a = - 4._dp * mC(i1) * Real( B0(p2, mC2(i1), 0._dp),dp )
   sumL = 0._dp
   sumR = 0._dp
   sumS = 0._dp
   Do i2 = 1,n_char
    Do i3 = 1,n_char
     sumL(i2,i3) = Conjg( c_CCZ_R(i2,i1) ) * c_CCZ_R(i3,i1) * B1m2   &
               & + Conjg( c_CCG_R(i2,i1) ) * c_CCG_R(i3,i1) * B1m2a
     sumR(i2,i3) = Conjg( c_CCZ_L(i2,i1) ) * c_CCZ_L(i3,i1) * B1m2   &
               & + Conjg( c_CCG_L(i2,i1) ) * c_CCG_L(i3,i1) * B1m2a
     sumS(i2,i3) = Conjg( c_CCZ_L(i2,i1) ) * c_CCZ_R(i3,i1) * B0m2   &
               & + Conjg( c_CCG_L(i2,i1) ) * c_CCG_R(i3,i1) * B0m2a
    End Do 
    If (WriteOut) Then
     Write (ErrCan,*) "C C Z/G",i2,i1
     Write (ErrCan,*) "  sig_L",sumL(i2,:)
     Write (ErrCan,*) "  sig_R",sumR(i2,:)
     Write (ErrCan,*) "  sig_S",sumS(i2,:)
    End If
   End Do
   sigL = sigL + sumL
   sigR = sigR + sumR
   sigS = sigS + sumS
  End Do
 !-------------------------------------
 ! chargino - P0
 !-------------------------------------
  Do i1 = 1,n_char
   Do i2 = 1,n_P0
    B1m2 = - 0.5_dp * Real( B1(p2,  mC2(i1), mP02(i2)),dp )
    B0m2 = mC(i1) * Real( B0(p2, mC2(i1), mP02(i2)),dp )
    sumL = 0._dp
    sumR = 0._dp
    sumS = 0._dp
    Do i3 = 1,n_char
     Do i4 = 1,n_char
      sumL(i3,i4) = Conjg( c_CCP0_L(i3,i1,i2) ) * c_CCP0_L(i4,i1,i2) * B1m2
      sumR(i3,i4) = Conjg( c_CCP0_R(i3,i1,i2) ) * c_CCP0_R(i4,i1,i2) * B1m2
      sumS(i3,i4) = Conjg( c_CCP0_R(i3,i1,i2) ) * c_CCP0_L(i4,i1,i2) * B0m2
     End Do 
     If (WriteOut) Then
      Write (ErrCan,*) "C C P0",i1,i3,i2
      Write (ErrCan,*) "  sig_L",sumL(i3,:)
      Write (ErrCan,*) "  sig_R",sumR(i3,:)
      Write (ErrCan,*) "  sig_S",sumS(i3,:)
     End If
    End Do 
    sigL = sigL + sumL
    sigR = sigR + sumR
    sigS = sigS + sumS
   End Do
  End Do

 !-------------------------------------
 ! chargino - S0
 !-------------------------------------
  Do i1 = 1,n_char
   Do i2 = 1,n_S0
    B1m2 = - 0.5_dp * Real( B1(p2,  mC2(i1), mS02(i2)),dp )
    B0m2 = mC(i1) * Real( B0(p2, mC2(i1), mS02(i2)),dp )
    sumL = 0._dp
    sumR = 0._dp
    sumS = 0._dp
    Do i3 = 1,n_char
     Do i4 = 1,n_char
      sumL(i3,i4) = Conjg( c_CCS0_L(i3,i1,i2) ) * c_CCS0_L(i4,i1,i2) * B1m2
      sumR(i3,i4) = Conjg( c_CCS0_R(i3,i1,i2) ) * c_CCS0_R(i4,i1,i2) * B1m2
      sumS(i3,i4) = Conjg( c_CCS0_R(i3,i1,i2) ) * c_CCS0_L(i4,i1,i2) * B0m2
     End Do 
     If (WriteOut) Then
      Write (ErrCan,*) "C C S0",i1,i3,i2
      Write (ErrCan,*) "  sig_L",sumL(i3,:)
      Write (ErrCan,*) "  sig_R",sumR(i3,:)
      Write (ErrCan,*) "  sig_S",sumS(i3,:)
     End If
    End Do 
    sigL = sigL + sumL
    sigR = sigR + sumR
    sigS = sigS + sumS
   End Do
  End Do

  sigL = oo16pi2 * SigL
  sigR = oo16pi2 * SigR
  sigs = oo16pi2 * SigS

 End Subroutine Sigma_Chargino


 Subroutine Sigma_SM_chirally_enhanced(gi, vevS, mu, V0, Yl, Yd, Yu    &
      & , Mi, Tl, Td, Tu, M2E, M2L, M2D, M2Q, M2U, epsD, epsL, epsD_FC &
      & , kont, RS0, RP0, m_d, C_ddH, C_ddA, m_l, C_llH, C_llA )
 !--------------------------------------------------------------------
 ! chirally enhanced terms of the d-quark selfenergies, based on
 ! A. Crivellin et al., hep-ph/1103.4272
 ! written by Werner Porod, 20.03.2014
 !--------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: gi(3), vevS(2)
  Complex(dp), Intent(in) :: Mi(3), mu
  Complex(dp), Intent(in), Dimension(3,3) :: V0, Yl, Yd, Yu, Tl, Td, Tu  &
    & , M2Q, M2D, M2U, M2E, M2L
  Real(dp), Intent(in), Optional :: RS0(2,2), RP0(2,2), m_d(3), m_l(3)

  Complex(dp), Intent(out) :: epsD(3), epsL(3), epsD_FC
  Complex(dp), Intent(out), Optional, Dimension(3,3,2) :: C_ddH, C_ddA  &
    & , C_llH, C_llA

  Complex(dp) :: sigma(3,3), alphaM(3)
  Complex(dp), Dimension(3,3) :: WdR, WdL, WuL, WuR, WlL, WlR, SigDLR  &
     & , SigULR, SigLLR, DelU, DelDnoY, DelDY, DelLnoY, DelLY, SigU    &
     & , SigDnoY, SigDY, SigLnoY, SigLY, E_ij, Etil_ij
  Complex(dp), Dimension(3,3,3) :: LamDLL, LamDRR, LamULL, LamURR, LamLLL &
     & , LamLRR
  Integer :: ierr, i1, i2, i3, i4, i5, i6
  Real(dp) :: mSdR2(3), mSqL2(3), mSuR2(3), mSlL2(3), mSlR2(3), test(2) &
     & , m12, m22, m32, alpha(3), C0m, vevSM(2), tanb

  Iname = Iname + 1
  NameOfUnit(Iname) = "Sigma_SM_chirally_enhanced"
  !------------------------------------
  ! some initialisaitons
  !------------------------------------
  alpha = oo4pi * gi**2

  vevSM = oosqrt2 * vevS ! they use a different normalisation

  kont = 0
  !------------------------------------
  ! diagonalisation of LL and RR parts
  !------------------------------------
  Call EigenSystem(M2D,mSdR2,WdR,ierr, test)
  Call Adjungate(WdR)
  Call EigenSystem(M2Q,mSqL2,WdL,ierr, test)
  Call Adjungate(WdL)
  WuL = Matmul(V0,WdL)
  Call EigenSystem(M2U,mSuR2,WuR,ierr, test)
  Call Adjungate(WuR)
  Call EigenSystem(M2E,mSlR2,WlR,ierr, test)
  Call Adjungate(WlR)
  Call EigenSystem(M2L,mSlL2,WlL,ierr, test)
  Call Adjungate(WlL)

  !-----------------------------------
  ! checking for negative mass squares
  !-----------------------------------
  If ( (mSdR2(1).Lt.0._dp).Or.(mSqL2(1).Lt.0._dp).Or.(mSuR2(1).Lt.0._dp) &
     & .Or.(mSlR2(1).Lt.0._dp).Or.(mSlL2(1).Lt.0._dp)) Then
   kont = -525
   Call AddError(525)
   Iname = Iname -1
   Return
  End If
  !-----------------------------------------------------------
  ! initialisation of basic quantities
  !-----------------------------------------------------------
  Do i1=1,3
   Do i2=1,3
    Do i3=1,3
     LamDLL(i1,i2,i3) = WdL(i2,i1) * Conjg(WdL(i3,i1))
     LamDRR(i1,i2,i3) = WdR(i2,i1) * Conjg(WdR(i3,i1))
     LamULL(i1,i2,i3) = WuL(i2,i1) * Conjg(WuL(i3,i1))
     LamURR(i1,i2,i3) = WuR(i2,i1) * Conjg(WuR(i3,i1))
     LamLLL(i1,i2,i3) = WlL(i2,i1) * Conjg(WlL(i3,i1))
     LamLRR(i1,i2,i3) = WlR(i2,i1) * Conjg(WlR(i3,i1))
    End Do
   End Do
  End Do
  DelU = vevSM(2) * Conjg(Tu) - vevSM(1) * mu * Conjg(Yu)
  DelDnoY = vevSM(1) * Conjg(Td)
  DelDY = - vevSM(2) * mu * Conjg(Yd)
  DelLnoY = vevSM(1) * Conjg(Tl)
  DelLY = - vevSM(2) * mu * Conjg(Yl)
 !-----------------------------------------------------------
  ! the different parts of the self energies
  !-----------------------------------------------------------
  SigDY = 0._dp
  SigDnoY = 0._dp
  SigLY = 0._dp
  SigLnoY = 0._dp
  SigU = 0._dp
  !--------------------------------------
  ! to avoid some unnessary calculations
  !--------------------------------------
  alphaM(1:2) = oo4pi * alpha(1:2) * mi(1:2)
  alphaM(3) = oo3pi * alpha(3) * mi(3)
  !------------------------
  ! double sfermion loops
  !------------------------
  Do i1=1,3
   Do i2=1,3
   !--------------------------
   ! gluino
   !--------------------------
    m12 = Abs(mi(3))**2
    m22 = mSqL2(i1)
    m32 = mSdR2(i2)
    C0m = 2._dp * alphaM(3) * C0_3m(m12, m22, m32)
    Do i3=1,3
     Do i4=1,3
      Do i5=1,3
       SigDY(i3,i4) = SigDY(i3,i4) + LamDLL(i1,i3,i5) * DelDY(i5,i5)   &
                    &                  * LamDRR(i2,i5,i4) * C0m
       Do i6=1,3
        SigDnoY(i3,i4) = SigDnoY(i3,i4) + LamDLL(i1,i3,i5) * DelDnoY(i5,i6) &
                       &                  * LamDRR(i2,i6,i4) * C0m
       End Do
      End Do
     End Do
    End Do
    m32 = mSuR2(i2)
    C0m = 2._dp * alphaM(3) * C0_3m(m12, m22, m32)
    Do i3=1,3
     Do i4=1,3
      Do i5=1,3
       Do i6=1,3
        SigU(i3,i4) = SigU(i3,i4) + LamULL(i1,i3,i5) * DelU(i5,i6) &
                       &                  * LamURR(i2,i6,i4) * C0m
       End Do
      End Do
     End Do
    End Do

   !--------------------------
   ! bino
   !--------------------------
    m12 = Abs(mi(1))**2
    m22 = mSqL2(i1)
    m32 = mSdR2(i2)
    C0m = - alphaM(1) * C0_3m(m12, m22, m32) / 9._dp
    Do i3=1,3
     Do i4=1,3
      Do i5=1,3
       SigDY(i3,i4) = SigDY(i3,i4) + LamDLL(i1,i3,i5) * DelDY(i5,i5)   &
                    &                  * LamDRR(i2,i5,i4) * C0m
       Do i6=1,3
        SigDnoY(i3,i4) = SigDnoY(i3,i4) + LamDLL(i1,i3,i5) * DelDnoY(i5,i6) &
                       &                  * LamDRR(i2,i6,i4) * C0m
       End Do
      End Do
     End Do
    End Do
    m32 = mSuR2(i2)
    C0m = 2._dp * alphaM(1) * C0_3m(m12, m22, m32) / 9._dp
    Do i3=1,3
     Do i4=1,3
      Do i5=1,3
       Do i6=1,3
        SigU(i3,i4) = SigU(i3,i4) + LamULL(i1,i3,i5) * DelU(i5,i6) &
                       &                  * LamURR(i2,i6,i4) * C0m
       End Do
      End Do
     End Do
    End Do
    m22 = mSlL2(i1)
    m32 = mSlR2(i2)
    C0m = alphaM(1) * C0_3m(m12, m22, m32) 
    Do i3=1,3
     Do i4=1,3
      Do i5=1,3
       SigLY(i3,i4) = SigLY(i3,i4) + LamLLL(i1,i3,i5) * DelLY(i5,i5)   &
                    &                  * LamLRR(i2,i5,i4) * C0m
       Do i6=1,3
        SigLnoY(i3,i4) = SigLnoY(i3,i4) + LamLLL(i1,i3,i5) * DelLnoY(i5,i6) &
                       &                  * LamLRR(i2,i6,i4) * C0m
       End Do
      End Do
     End Do
    End Do

   !--------------------------
   ! charged higgsino
   !--------------------------
    m12 = Abs(mu)**2
    m22 = mSqL2(i1)
    m32 = mSuR2(i2)
    C0m = - mu * Yu(3,3) * Yd(3,3) * C0_3m(m12, m22, m32) * V0(3,3) &
        &      * oo16pi2 * LamDLL(i1,3,3) * Conjg(DelU(3,3)) * LamuRR(i2,3,3)
    Do i3=1,3
     SigDY(i3,3) = SigDY(i3,3) + Conjg(V0(3,i3)) * C0m
    End Do

   End Do
  End Do
  !------------------------
  ! single sfermion loop
  !------------------------
  Do i1=1,3
   !-------------------------------------
   ! charged + neutral wino and sleptons
   !-------------------------------------
   m12 = Abs(mi(2))**2
   m22 = Abs(mu)**2
   m32 = mSlL2(i1)
   C0m = mu * alphaM(2) * C0_3m(m12, m22, m32) * vevSM(2) * 1.5_dp
   Do i2=1,3
    Do i3=1,3
     SigLY(i2,i3) = SigLY(i2,i3) + LamLLL(i1,i2,i3) * Yl(i3,i3) * C0m
    End Do
   End Do

   !-------------------------------------
   ! bino and L-sleptons
   !-------------------------------------
   m12 = Abs(mi(1))**2
   m22 = Abs(mu)**2
   m32 = mSlL2(i1)
   C0m = - mu * alphaM(1) * C0_3m(m12, m22, m32) * vevSM(2) * 0.5_dp
   Do i2=1,3
    Do i3=1,3
     SigLY(i2,i3) = SigLY(i2,i3) + LamLLL(i1,i2,i3) * Yl(i3,i3) * C0m
    End Do
   End Do

   !-------------------------------------
   ! bino and R-sleptons
   !-------------------------------------
   m12 = Abs(mi(1))**2
   m22 = Abs(mu)**2
   m32 = mSlR2(i1)
   C0m = mu * alphaM(1) * C0_3m(m12, m22, m32) * vevSM(2) 
   Do i2=1,3
    Do i3=1,3
     SigLY(i2,i3) = SigLY(i2,i3) + LamLRR(i1,i2,i3) * Yl(i2,i2) * C0m
    End Do
   End Do

   !---------------------------------------------
   ! charged + neutral wino and left D/U-squark
   !---------------------------------------------
   m12 = Abs(mi(2))**2
   m22 = Abs(mu)**2
   m32 = mSqL2(i1)
   C0m = mu * alphaM(2) * C0_3m(m12, m22, m32) * vevSM(2) * 1.5_dp
   Do i2=1,3
    Do i3=1,3
     SigDY(i2,i3) = SigDY(i2,i3) + LamDLL(i1,i2,i3) * Yd(i3,i3) * C0m
    End Do
   End Do
   !-------------------------------------
   ! bino and left D-squark
   !-------------------------------------
   m12 = Abs(mi(1))**2
   m22 = Abs(mu)**2
   m32 = mSqL2(i1)
   C0m = mu * alphaM(1) * C0_3m(m12, m22, m32) * vevSM(2) / 6._dp
   Do i2=1,3
    Do i3=1,3
     SigDY(i2,i3) = SigDY(i2,i3) + LamDLL(i1,i2,i3) * Yd(i3,i3) * C0m
    End Do
   End Do

   !-------------------------------------
   ! bino and right D-squarks
   !-------------------------------------
   m12 = Abs(mi(1))**2
   m22 = Abs(mu)**2
   m32 = mSdR2(i1)
   C0m = mu * alphaM(1) * C0_3m(m12, m22, m32) * vevSM(2) / 3._dp
   Do i2=1,3
    Do i3=1,3
     SigDY(i2,i3) = SigDY(i2,i3) + LamDRR(i1,i2,i3) * Yd(i2,i2) * C0m
    End Do
   End Do

  End Do

  Do i1=1,3
   epsD(i1) = SigDY(i1,i1) / (vevSM(2)*Yd(i1,i1) )
   epsL(i1) = SigLY(i1,i1) / (vevSM(2)*Yl(i1,i1) )
  End Do

  epsd_FC = 0._dp
  m12 = Abs(mu)**2
  Do i1=1,3
   Do i2=1,3
   !--------------------------
   ! charged higgsino
   !--------------------------
    m22 = mSqL2(i1)
    m32 = mSuR2(i2)
    C0m = - mu * Yu(3,3) * Yd(3,3) * C0_3m(m12, m22, m32) &
        &      * oo16pi2
    epsd_FC = epsd_FC + LamdLL(i1,3,3) * Conjg(DelU(3,3)) * LamuRR(i2,3,3) * C0m
   End Do
  End Do
  epsd_FC = epsd_FC /(Yd(3,3)*vevSM(1))

  If (Present(C_ddH).Or.Present(C_ddA)) Then
   tanb = vevSM(2)/vevSM(1)
   E_ij = SigDY / vevSM(2)

   sigma = 0._dp
   sigma(1,2) = (SigDY(1,2) + SigDnoY(1,2)) / m_d(2) 
   sigma(1,3) = (SigDY(1,3) + SigDnoY(1,3)) / m_d(3) 
   sigma(2,3) = (SigDY(2,3) + SigDnoY(2,3)) / m_d(3) 
   sigma(2,1) = (SigDY(2,1) + SigDnoY(2,1)) / m_d(2) 
   sigma(3,1) = (SigDY(3,1) + SigDnoY(3,1)) / m_d(3) 
   sigma(3,2) = (SigDY(3,2) + SigDnoY(3,2)) / m_d(3)

   Etil_ij = 0._dp

   Etil_ij(1,2) = E_ij(1,2) - sigma(1,2) * E_ij(2,2)
   Etil_ij(1,3) = E_ij(1,3)                                   &
              & - (sigma(1,3)-sigma(1,2)*sigma(2,3)) * E_ij(3,3) &
              & - sigma(1,2) * E_ij(2,3)
   Etil_ij(2,3) = E_ij(2,3) - sigma(2,3) * E_ij(3,3)
   Etil_ij(2,1) = E_ij(2,1) - sigma(2,1) * E_ij(2,2)
   Etil_ij(3,1) = E_ij(3,1)                                   &
              & - (sigma(3,1)-sigma(2,1)*sigma(3,2)) * E_ij(3,3) &
              & - sigma(2,1) * E_ij(3,2)
   Etil_ij(3,2) = E_ij(3,2) - sigma(3,2) * E_ij(3,3)

   If (Present(C_ddH)) Then
    C_ddH = 0._dp
    Do i3=1,2
     Do i1=1,3
      C_ddH(i1,i1,i3) = - m_d(i1) * RS0(i3,1) * oosqrt2 / vevSM(1)
      Do i2=1,3
       If (i1.Ne.i2) Then
        C_ddH(i1,i2,i3) = ( Etil_ij(i1,i2)*tanb * RS0(i3,1)  &
                     &    - Etil_ij(i1,i2)* RS0(i3,2) ) * oosqrt2
       End If
      End Do
     End Do
    End Do
   End If

   If (Present(C_ddA)) Then
    C_ddA = 0._dp
    Do i3=1,2
     Do i1=1,3
      C_ddA(i1,i1,i3) = - m_d(i1)  * Cmplx(0._dp,RP0(i3,1),dp) * oosqrt2 / vevSM(1)
      Do i2=1,3
       If (i1.Ne.i2) Then
        C_ddA(i1,i2,i3) = ( Etil_ij(i1,i2)*tanb * Cmplx(0._dp,RP0(i3,1),dp) &
                       &  - Etil_ij(i1,i2)* Cmplx(0._dp,RP0(i3,2),dp) ) * oosqrt2
       End If
      End Do
     End Do
    End Do

   End If
  End If

  If (Present(C_llH).Or.Present(C_llA)) Then
   tanb = vevSM(2)/vevSM(1)
   E_ij = SigLY / vevSM(2)
   sigma(1,2) = (SigLY(1,2) + SigLnoY(1,2)) / m_l(2) 
   sigma(1,3) = (SigLY(1,3) + SigLnoY(1,3)) / m_l(3) 
   sigma(2,3) = (SigLY(2,3) + SigLnoY(2,3)) / m_l(3) 
   sigma(2,1) = (SigLY(2,1) + SigLnoY(2,1)) / m_l(2) 
   sigma(3,1) = (SigLY(3,1) + SigLnoY(3,1)) / m_l(3) 
   sigma(3,2) = (SigLY(3,2) + SigLnoY(3,2)) / m_l(3) 
   Etil_ij = E_ij
   Etil_ij(1,2) = Etil_ij(1,2) + sigma(1,2) * E_ij(2,2)
   Etil_ij(1,3) = Etil_ij(1,3)                                   &
              & + (sigma(1,3)-sigma(1,2)*sigma(2,3)) * E_ij(3,3) &
              & + sigma(1,2) * E_ij(2,3)
   Etil_ij(2,3) = Etil_ij(2,3) + sigma(2,3) * E_ij(3,3)
   Etil_ij(2,1) = Etil_ij(2,1) + sigma(2,1) * E_ij(2,2)
   Etil_ij(3,1) = Etil_ij(3,1)                                   &
              & + (sigma(3,1)-sigma(2,1)*sigma(3,2)) * E_ij(3,3) &
              & + sigma(2,1) * E_ij(3,2)
   Etil_ij(3,2) = Etil_ij(3,2) + sigma(3,2) * E_ij(3,3)
   If (Present(C_llH)) Then
    C_llH = 0._dp
    Do i3=1,2
     Do i1=1,3
      C_llH(i1,i1,i3) = - m_l(i1) * RS0(i3,1) * oosqrt2 / vevSM(1)
      Do i2=1,3
       If (i1.Ne.i2) Then
        C_llH(i1,i2,i3) = ( Etil_ij(i1,i2)*tanb * RS0(i3,1)  &
                     &    - Etil_ij(i1,i2)* RS0(i3,2) ) * oosqrt2
       End If
      End Do
     End Do
    End Do
   End If
   If (Present(C_llA)) Then
    C_llA = 0._dp
    Do i3=1,2
     Do i1=1,3
      C_llA(i1,i1,i3) = - m_l(i1)  * Cmplx(0._dp,RP0(i3,1),dp) * oosqrt2/ vevSM(1)
      Do i2=1,3
       If (i1.Ne.i2) Then
        C_llA(i1,i2,i3) = ( Etil_ij(i1,i2)*tanb * Cmplx(0._dp,RP0(i3,1),dp) &
                       &  - Etil_ij(i1,i2)* Cmplx(0._dp,RP0(i3,2),dp) ) * oosqrt2
       End If
      End Do
     End Do
    End Do
   End If
  End If

  Iname = Iname - 1

 End Subroutine Sigma_SM_chirally_enhanced


 Subroutine Sigma_Fermion3(p2, mf, yuk_f, RfL, RfR, gSU2, gSU3, sinW2 ,T3, e  &
        & , mfp, yuk_fp, RfpL, RfpR, mSf2, RSf, mSf2p, Rsfp, mglu , phase_glu &
        & , mN, mN2, N, mC, mC2, U, V, mS02, RS0, mP02, RP0,mSpm2 , RSpm      &
        & , mZ2, mW2, l_SU3, SigS, SigL, SigR, SigQCD)
 !-----------------------------------------------------------------------
 ! Calculates the 1-loop self energy of SM-fermions in the three generation
 ! case. 
 ! The formulas of J. Bagger et al, Nucl.Phys.B are usedas basis.
 ! The renormalization scale has to be set in the main program
 !  due to the structure of Loopfuntcions. 
 ! the input is:
 ! - p2 ......... the outer momentum squared 
 ! - gSU3 ....... the SU(3) gauge coupling at p2
 ! - gSU2 ....... the SU(2) gauge coupling at p2
 ! - sinW2 ...... sin(theta_W) squared
 ! - T3 ......... weak isospin of the fermion
 ! - e .......... charge of the fermion
 ! - mf ......... fermion mass
 ! - mSf(i) ..... Sfermion masses
 ! - RSf(i,j) ... mixing matrix of the sfermions
 ! - mfp ........ fermion mass'
 ! - mSfp(i) .... Sfermion masses'
 ! - RSfp(i,j) .. mixing matrix of the sfermions'
 ! - mglu ....... gluino mass
 ! - phase_glu .. phase of the gluino 
 ! - mZ2 ........ Z-boson mass
 ! - mW2 ........ W-boson mass
 ! - mS02 ....... masses of the scalar Higgs
 ! - RS0 ........ mixing matrix of the scalar Higgs
 ! - mP02 ....... masses of the pseudoscalar Higgs
 ! - RP0 ........ mixing matrix of the pseudoscalar Higgs
 ! - mSpm2 ...... masses of the charged Higgs
 ! - RSpm ....... mixing matrix of the charged Higgs
 ! - mN ......... masses of the neutralinos
 ! - N .......... mixing matrix of the neutrlinos
 ! - mC ......... masses of the charginos
 ! - U,V ........ mixing matrices of the charginos
 !    
 ! written by Werner Porod, 04.02.03
 ! taking the routine Sigma_Fermion as basis. The pure gluino QCD part is
 ! given extra, because heavy and light quarks have to be treated different
 ! due potentially large logs, also the photon contributions are added there
 ! 05.01.04: taking gauge boson masses as input to allow for running masses
 ! 11.06.04: be aware that for 2-component spinors
 !             Y_ij f_L_i f_R_j
 !           gives a transpose for 4-compent spinors (f_L, f^*_R)
 !             Y_ij bar(f)_j P_L f_i
 !           and thus a somewhat unusual index structure
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: p2, gSU3, gSU2, sinW2, mS02(:), RS0(:,:), mP02(:) &
     & ,RP0(:,:), mSpm2(:), mN(:), mN2(:), mC(:), mC2(:), mf(3), mSf2(6)    &
     & ,mfp(3), mSf2p(:), e, T3, mglu, mZ2, mW2
  Complex(dp), Intent(in) :: N(:,:), U(:,:), V(:,:), RSpm(:,:), Rsf(6,6)  &
     & ,phase_glu, Rsfp(:,:), yuk_f(3,3), yuk_fp(3,3), RfL(3,3), RfR(3,3) &
     & , RfpL(3,3), RfpR(3,3)
  Logical, Intent(in) :: l_SU3
  Complex(dp), Dimension(3,3), Intent(out) :: SigS, SigL, SigR
  Real(dp), Optional, Intent(out) :: SigQCD

  Integer :: i1, i2, i3, i4, n_fp
  Real(dp) :: cosW2, coupL, coupR, mf2(3), mfp2(3), mglu2, gp, f_SU3 &
     & , Q2, logQ, B1m2, B0m2, B1m2a, B0m2a
  Complex(dp) :: sumS(3,3), sumL(3,3), sumR(3,3), c_L(3), c_R(3)
  Logical :: WriteOut 
  Real(dp), Parameter :: &
    & as2loop = 1._dp / 24._dp + 2011._dp * oo32Pi2 / 12._dp           &
    &         + Log2 / 12._dp - oo8Pi2 * Zeta3                        &
    & , log2loop_a = 123._dp * oo32Pi2, log2loop_b = 33._dp * oo32Pi2

 !----------------
 ! Initialization
 !----------------
  Iname = Iname + 1
  NameOfUnit(Iname) ="Sigma_Fermion3"
  If ((WriteOneLoopContributions.Eq.10).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions in Sigma_Fermion3:",mf
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If

  Q2 = GetRenormalizationScale() ! from LoopTools
  cosW2 = 1._dp - sinW2
  mf2 = mf**2
  mfp2 = mfp**2
  gp = gSU2 * Sqrt(sinW2/cosW2)
  SigS = ZeroC
  SigL = ZeroC
  SigR = ZeroC
  !-------------------------------
  ! strong interactions + photon
  !-------------------------------
  If (l_SU3) Then
   mglu2 = mglu**2
   f_SU3 = 4._dp * gSU3**2 / 3._dp
   If (Present(SigQCD)) Then ! need 2-loop for top-quark
    logQ = Log(Q2/mf2(3))
    SigQCD = - f_SU3 * mf(3) * (5._dp + 3._dp * LogQ                   &
         &                         + (as2loop + log2loop_a * logQ      &
         &                         + log2loop_b * logQ**2) * gSU3**2 ) &
         & - mf(3) * (e*gSU2)**2 * sinW2 * (5._dp + 3._dp * LogQ )
    If (WriteOut) Write(ErrCan,*) "gluon :",SigQCD,gsu3
    SigQCD = oo16pi2 * SigQCD
   End If

   Do i1=1,6
    B0m2 = 16._dp * mglu * B0(p2,mglu2, mSf2(i1)) / 3._dp
    B1m2 = - 8._dp * B1(p2,mglu2, mSf2(i1)) / 3._dp
    Do i2=1,3
     Call CoupGluinoSquark(gSU3, phase_glu, i2, i1, RSf, id3C, id3C &
                          &, c_L(i2), c_R(i2))
    End Do    
    Do i2=1,3
     B0m2a = 16._dp * mglu * B0(mf2(i2),mglu2, mSf2(i1)) / 3._dp
     B1m2a = - 8._dp * B1(mf2(i2),mglu2, mSf2(i1)) / 3._dp
  !---------------------------------------------------------------
  ! note the somewhat misleading notation, because in coupL, coupR
  ! the chirality of the gluino is meant with L, R
  !---------------------------------------------------------------
     SumR(i2,i2) = Abs( c_L(i2) )**2 * B1m2a
     SumL(i2,i2) = Abs( c_R(i2) )**2 * B1m2a
     SumS(i2,i2) = Conjg( c_R(i2) ) * c_L(i2) * B0m2a
     Do i3=1,3
      If (i2.Ne.i3) Then
       SumR(i2,i3) = Conjg( c_L(i2) ) * c_L(i3) * B1m2
       SumL(i2,i3) = Conjg( c_R(i2) ) * c_R(i3) * B1m2
       SumS(i2,i3) = Conjg( c_R(i2) ) * c_L(i3) * B0m2
      End If
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "gluino squark quark",i1,i2
      Write(ErrCan,*) " sig_L ",sumL(i2,:)
      Write(ErrCan,*) " sig_R ",sumR(i2,:)
      Write(ErrCan,*) " sig_S ",sumS(i2,:)
     End If
    End Do
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End If
  !--------------
  ! gauge bosons
  !--------------
   Do i1=1,3
    B1m2 = - 0.5_dp * gSU2**2 * B1(p2,mfp2(i1), mW2)
    c_L = RfpL(i1,:) 
    Do i2=1,3
     B1m2a = - 0.5_dp * gSU2**2 * B1(mf2(i2),mfp2(i1), mW2)
     Do i3=1,3
      If (i2.Eq.i3) Then
       sumL(i2,i3) = c_L(i3) * Conjg(c_L(i2) ) * B1m2a
      Else
       sumL(i2,i3) = c_L(i3) * Conjg(c_L(i2) ) * B1m2
      End If
     End Do
     If (WriteOut) Then
      Write (ErrCan,*) "W f f'",i1,i2
      Write(ErrCan,*) " sig_L ",sumL(i2,:)
     End If
    End Do
    SigL = SigL + sumL
   End Do

  Call CoupFermionZ(T3,e,gSU2,sinW2,coupL,coupR)
  sumL =0._dp
  sumR =0._dp
  sumS =0._dp
  Do i1=1,3
   C_L = coupL * RfL(i1,:) 
   C_R = coupR * RfR(i1,:)
   B1m2 = - B1(p2,mf2(i1), mZ2)
   B0m2 = - 4._dp * mf(i1) * B0(p2,mf2(i1), mZ2)
   Do i2=1,3
    B1m2a = - B1(mf2(i2),mf2(i1), mZ2)
    B0m2a = - 4._dp * mf(i1) * B0(mf2(i2),mf2(i1), mZ2)
    sumL(i2,i2) = Abs(c_L(i2))**2 * B1m2a
    sumR(i2,i2) = Abs(c_R(i2))**2 * B1m2a
    sumS(i2,i2) = c_L(i2) * Conjg(c_R(i2)) * B0m2a
    Do i3=1,3
     If (i2.Ne.i3) Then
      sumL(i2,i3) = c_L(i3) * Conjg(c_L(i2) ) * B1m2
      sumR(i2,i3) = c_R(i3) * Conjg(c_R(i2) ) * B1m2
      sumS(i2,i3) = c_L(i3) * Conjg(c_R(i2)) * B0m2 
     End If
    End Do
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "Z-boson",i1
     Write(ErrCan,*) " sig_L ",sumL(i1,:)
     Write(ErrCan,*) " sig_R ",sumR(i1,:)
     Write(ErrCan,*) " sig_S ",sumS(i1,:)
   End If
   SigL = SigL + sumL
   SigR = SigR + sumR
   SigS = SigS + sumS
  End Do
  !----------------------
  ! charged Higgs bosons
  !---------------------
  Do i1=1,n_Spm
   Do i2=1,3
    B1m2 = - 0.5_dp * B1(p2,mfp2(i2), mSpm2(i1))
    B0m2 = mfp(i2) * B0(p2,mfp2(i2), mSpm2(i1))
    Do i3=1,3
     If (T3.Gt.0._dp) Then
      Call CoupChargedScalarFermion(i1, i2, i3, RSpm, yuk_fp, RfpL, RfpR &
                                  & , yuk_f, id3C, id3C, c_L(i3), c_R(i3))
     Else
      Call CoupChargedScalarFermion(i1, i3, i2, RSpm, yuk_f, id3C, id3C &
                                  & , yuk_fp, RfpL, RfpR, c_R(i3), c_L(i3))
      c_R(i3) = Conjg(c_R(i3))
      c_L(i3) = Conjg(c_L(i3))
     End If
    End Do
    Do i3=1,3
     B1m2a = - 0.5_dp * B1(mf2(i3),mfp2(i2), mSpm2(i1))
     B0m2a = mfp(i2) * B0(mf2(i3),mfp2(i2), mSpm2(i1))
     SumL(i3,i3) = Abs(c_L(i3))**2 * B1m2a
     SumR(i3,i3) = Abs(c_R(i3))**2 * B1m2a
     SumS(i3,i3) = c_L(i3)*Conjg(c_R(i3)) * B0m2a
     Do i4=1,3
      If (i3.Ne.i4) Then
       SumL(i3,i4) = c_L(i4) * Conjg( c_L(i3) ) * B1m2
       SumR(i3,i4) = c_R(i4) * Conjg( c_R(i3) ) * B1m2
       SumS(i3,i4) = c_L(i4) * Conjg( c_R(i3) ) * B0m2
      End If
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "S+ f'",i1,i2,i3
      Write(ErrCan,*) "sumL",sumL(i3,:)
      Write(ErrCan,*) "sumR",sumR(i3,:)
      Write(ErrCan,*) "sumS",sumS(i3,:)
     End If
    End Do
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End Do 

  !-------------
  !  charginos
  !-------------
  n_fp = 6 ! number of sfermions'
  If (e.Eq.-1._dp) n_fp = 3  ! sneutrinos
  Do i1=1,n_char
   Do i2=1,n_fp
    B1m2 = - 0.5_dp * B1(p2,mC2(i1), mSf2p(i2))
    B0m2 = mC(i1) * B0(p2,mC2(i1), mSf2p(i2))
    Do i3=1,3
     If (T3.Gt.0._dp) Then
      Call CoupCharginoSfermion(i1, i3, i2, gSU2, T3, RSfp, Yuk_fp, Yuk_f &
                                &, id3C, id3C, U, V, c_L(i3), c_R(i3))
     Else
      Call CoupCharginoSfermion(i1, i3, i2, gSU2, T3, RSfp, Yuk_f, Yuk_fp &
                                &, id3C, id3C, U, V, c_L(i3), c_R(i3))
     End If
    End Do
    Do i3=1,3
     B1m2a = - 0.5_dp * B1(mf2(i3),mC2(i1), mSf2p(i2))
     B0m2a = mC(i1) * B0(mf2(i3),mC2(i1), mSf2p(i2))
  !---------------------------------------------------------------
  ! note the somewhat misleading notation, because in coupL, coupR
  ! the chirality of the chargino is meant with L, R
  !---------------------------------------------------------------
     SumR(i3,i3) = Abs(c_L(i3))**2 * B1m2a
     SumL(i3,i3) = Abs(c_R(i3))**2 * B1m2a
     SumS(i3,i3) = c_L(i3)*Conjg(c_R(i3)) * B0m2a
     Do i4=1,3
      If (i3.Ne.i4) Then
       SumR(i3,i4) = c_L(i4) * Conjg( c_L(i3) ) * B1m2
       SumL(i3,i4) = c_R(i4) * Conjg( c_R(i3) ) * B1m2
       SumS(i3,i4) = c_L(i4) * Conjg( c_R(i3) ) * B0m2
      End If
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "Chi Sf'",i1,i2,i3
      Write(ErrCan,*) "sumL",sumL(i3,:)
      Write(ErrCan,*) "sumR",sumR(i3,:)
      Write(ErrCan,*) "sumS",sumS(i3,:)
     End If
    End Do
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End Do

  !------------------------------
  ! neutral Higgs contributions
  !------------------------------
  Do i1=1,n_P0
   Do i2=1,3
    B1m2 = - 0.5_dp * B1(p2,mf2(i2), mP02(i1))
    B0m2 = mf(i2) * B0(p2,mf2(i2), mP02(i1))
    Do i3=1,3
     Call CoupFermionPseudoScalar3a(i3, i2, i1, T3, yuk_f, Rfl, RfR, RP0 &
                                &, c_L(i3), c_R(i3))
    End Do
    Do i3=1,3
     B1m2a = - 0.5_dp * B1(mf2(i3),mf2(i2), mP02(i1))
     B0m2a = mf(i2) * B0(mf2(i3),mf2(i2), mP02(i1))
     SumR(i3,i3) = Abs(c_L(i3))**2 * B1m2a
     SumL(i3,i3) = Abs(c_R(i3))**2 * B1m2a
     SumS(i3,i3) = c_L(i3)*Conjg(c_R(i3)) * B0m2a
     Do i4=1,3
      If (i3.Ne.i4) Then
       SumR(i3,i4) = c_L(i4) * Conjg( c_L(i3) ) * B1m2
       SumL(i3,i4) = c_R(i4) * Conjg( c_R(i3) ) * B1m2
       SumS(i3,i4) = c_L(i4) * Conjg( c_R(i3) ) * B0m2
      End If
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "P0 f",i1,i2,i3
      Write(ErrCan,*) "sumL",sumL(i3,:)
      Write(ErrCan,*) "sumR",sumR(i3,:)
      Write(ErrCan,*) "sumS",sumS(i3,:)
     End If
    End Do
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End Do

  Do i1=1,n_S0
   Do i2=1,3
    B1m2 = - 0.5_dp * B1(p2,mf2(i2), mS02(i1))
    B0m2 = mf(i2) * B0(p2,mf2(i2), mS02(i1))
    Do i3=1,3
     Call CoupFermionScalar3a(i3, i2, i1, T3, yuk_f, Rfl, RfR, RS0 &
                          &, c_L(i3), c_R(i3))
    End Do
    Do i3=1,3
     B1m2a = - 0.5_dp * B1(mf2(i3),mf2(i2), mS02(i1))
     B0m2a = mf(i2) * B0(mf2(i3),mf2(i2), mS02(i1))
     SumR(i3,i3) = Abs(c_L(i3))**2 * B1m2a
     SumL(i3,i3) = Abs(c_R(i3))**2 * B1m2a
     SumS(i3,i3) = c_L(i3)*Conjg(c_R(i3)) * B0m2a
     Do i4=1,3
      If (i3.Ne.i4) Then
       SumR(i3,i4) = c_L(i4) * Conjg( c_L(i3) ) * B1m2
       SumL(i3,i4) = c_R(i4) * Conjg( c_R(i3) ) * B1m2
       SumS(i3,i4) = c_L(i4) * Conjg( c_R(i3) ) * B0m2
      End If
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "S0 f",i1,i2,i3
      Write(ErrCan,*) "sumL",sumL(i3,:)
      Write(ErrCan,*) "sumR",sumR(i3,:)
      Write(ErrCan,*) "sumS",sumS(i3,:)
     End If
    End Do
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End Do

  !-------------
  ! neutralinos
  !-------------
  Do i1=1,n_neut
   Do i2=1,6
    B1m2 = - 0.5_dp * B1(p2,mN2(i1), mSf2(i2))
    B0m2 = mN(i1) * B0(p2,mN2(i1), mSf2(i2))
    Do i3=1,3
     Call CoupNeutralinoSfermion(i3, i1, i2, gp, gSU2, T3, e, RSf, id3C &
                               &, id3C, Yuk_f, N, c_L(i3), c_R(i3))
    End Do
    Do i3=1,3
     B1m2a = - 0.5_dp * B1(mf2(i3),mN2(i1), mSf2(i2))
     B0m2a = mN(i1) * B0(mf2(i3),mN2(i1), mSf2(i2))
  !---------------------------------------------------------------
  ! note the somewhat misleading notation, because in coupL, coupR
  ! the chirality of the neutralino is meant with L, R
  !---------------------------------------------------------------
     SumR(i3,i3) = Abs(c_L(i3))**2 * B1m2a 
     SumL(i3,i3) = Abs(c_R(i3))**2 * B1m2a
     SumS(i3,i3) = c_L(i3)*Conjg(c_R(i3)) * B0m2a
     Do i4=1,3
      If (i3.Ne.i4) Then
       SumR(i3,i4) = c_L(i4) * Conjg( c_L(i3) ) * B1m2
       SumL(i3,i4) = c_R(i4) * Conjg( c_R(i3) ) * B1m2
       SumS(i3,i4) = c_L(i4) * Conjg( c_R(i3) ) * B0m2
      End If
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "Chi0 Sf",i1,i2,i3
      Write(ErrCan,*) "sumL",sumL(i3,:)
      Write(ErrCan,*) "sumR",sumR(i3,:)
      Write(ErrCan,*) "sumS",sumS(i3,:)
     End If
    End Do
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End Do

  SigL = oo16pi2 * SigL
  SigR = oo16pi2 * SigR
  SigS = oo16pi2 * SigS

  Iname = Iname - 1

 End Subroutine Sigma_Fermion3


 Subroutine Sigma_Fermion3sckm(p2, gSU2, gSU3, sinW2 ,T3, e, CKM, mf, Y_f_in  &
       & , mfp, Y_fp_in, mSf2, RSf, mSf2p, Rsfp, mglu , phase_glu, mN, mN2, N &
       & , mC, mC2, U, V, mS02, RS0, mP02, RP0 ,mSpm2, RSpm, mZ2, mW2, l_SU3  &
       & , SigS, SigL, SigR, SigQCD)
 !-----------------------------------------------------------------------
 ! Calculates the 1-loop self energy of SM-fermions in the three generation
 ! case in the superCKM basis, taking Sigma_Fermion3 as starting point
 ! The formulas are based on J. Bagger et al, NPB491 (1997) 3.
 ! The renormalization scale has to be set in the main program
 !  due to the structure of Loopfuntcions. 
 ! the input is:
 ! - p2 ......... the outer momentum squared 
 ! - gSU3 ....... the SU(3) gauge coupling at p2
 ! - gSU2 ....... the SU(2) gauge coupling at p2
 ! - sinW2 ...... sin(theta_W) squared
 ! - T3 ......... weak isospin of the fermion
 ! - e .......... charge of the fermion
 ! - CKM ........ CKM/PMNS matrix depending on whether quarks or leptons
 ! - mf ......... fermion mass
 ! - mSf(i) ..... Sfermion masses
 ! - RSf(i,j) ... mixing matrix of the sfermions
 ! - mfp ........ fermion mass'
 ! - mSfp(i) .... Sfermion masses'
 ! - RSfp(i,j) .. mixing matrix of the sfermions'
 ! - mglu ....... gluino mass
 ! - phase_glu .. phase of the gluino 
 ! - mZ2 ........ Z-boson mass
 ! - mW2 ........ W-boson mass
 ! - mS02 ....... masses of the scalar Higgs
 ! - RS0 ........ mixing matrix of the scalar Higgs
 ! - mP02 ....... masses of the pseudoscalar Higgs
 ! - RP0 ........ mixing matrix of the pseudoscalar Higgs
 ! - mSpm2 ...... masses of the charged Higgs
 ! - RSpm ....... mixing matrix of the charged Higgs
 ! - mN ......... masses of the neutralinos
 ! - N .......... mixing matrix of the neutrlinos
 ! - mC ......... masses of the charginos
 ! - U,V ........ mixing matrices of the charginos
 !    
 ! written by Werner Porod, 03.04.14
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: p2, gSU3, gSU2, sinW2, mS02(:), RS0(:,:), mP02(:) &
     & ,RP0(:,:), mSpm2(:), mN(:), mN2(:), mC(:), mC2(:), mf(3), mSf2(6)    &
     & ,mfp(3), mSf2p(:), e, T3, mglu, mZ2, mW2
  Complex(dp), Intent(in) :: N(:,:), U(:,:), V(:,:), RSpm(:,:), Rsf(6,6)  &
     & ,phase_glu, Rsfp(:,:), Y_f_in(3,3), Y_fp_in(3,3), CKM(3,3)
  Logical, Intent(in) :: l_SU3
  Complex(dp), Dimension(3,3), Intent(out) :: SigS, SigL, SigR
  Real(dp), Optional, Intent(out) :: SigQCD

  Integer :: i1, i2, i3, i4, n_fp
  Real(dp) :: cosW2, coupL, coupR, mf2(3), mfp2(3), mglu2, gp, f_SU3 &
     & , Q2, logQ, B1m2, B0m2, B1m2a, B0m2a
  Complex(dp) :: sumS(3,3), sumL(3,3), sumR(3,3), c_L(3), c_R(3), yuk_f(3,3) &
     & , yuk_fp(3,3)
  Logical :: WriteOut 
  Real(dp), Parameter :: &
    & as2loop = 1._dp / 24._dp + 2011._dp * oo32Pi2 / 12._dp           &
    &         + Log2 / 12._dp - oo8Pi2 * Zeta3                        &
    & , log2loop_a = 123._dp * oo32Pi2, log2loop_b = 33._dp * oo32Pi2

 !----------------
 ! Initialization
 !----------------
  Iname = Iname + 1
  NameOfUnit(Iname) ="Sigma_Fermion3sckm"
  If ((WriteOneLoopContributions.Eq.10).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions in Sigma_Fermion3:",mf
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If

  Q2 = GetRenormalizationScale() ! from LoopTools
  cosW2 = 1._dp - sinW2
  mf2 = mf**2
  mfp2 = mfp**2
  gp = gSU2 * Sqrt(sinW2/cosW2)
  SigS = ZeroC
  SigL = ZeroC
  SigR = ZeroC
  !----------------------------------------------------------------------
  ! to be able to re-use the orginal coupling routines for charged Higgs
  ! bosons and charginos
  !----------------------------------------------------------------------
  If (T3.Gt.0._dp) Then
   yuk_f = Matmul(Transpose(CKM),Y_f_in)
   yuk_fp = Y_fp_in
  Else
   yuk_fp = Matmul(Transpose(CKM),Y_fp_in)
   yuk_f = Y_f_in
  End If
  !-------------------------------
  ! strong interactions + photon
  !-------------------------------
  If (l_SU3) Then
   mglu2 = mglu**2
   f_SU3 = 4._dp * gSU3**2 / 3._dp
   If (Present(SigQCD)) Then ! need 2-loop for top-quark
    logQ = Log(Q2/mf2(3))
    SigQCD = - f_SU3 * mf(3) * (5._dp + 3._dp * LogQ                   &
         &                         + (as2loop + log2loop_a * logQ      &
         &                         + log2loop_b * logQ**2) * gSU3**2 ) &
         & - mf(3) * (e*gSU2)**2 * sinW2 * (5._dp + 3._dp * LogQ )
    If (WriteOut) Write(ErrCan,*) "gluon :",SigQCD,gsu3
    SigQCD = oo16pi2 * SigQCD
   End If

   Do i1=1,6
    B0m2 = 16._dp * mglu * B0(p2,mglu2, mSf2(i1)) / 3._dp
    B1m2 = - 8._dp * B1(p2,mglu2, mSf2(i1)) / 3._dp
    Do i2=1,3
     Call CoupGluinoSquark(gSU3, phase_glu, i2, i1, RSf, id3C, id3C &
                          &, c_L(i2), c_R(i2))
    End Do    
    Do i2=1,3
     B0m2a = 16._dp * mglu * B0(mf2(i2),mglu2, mSf2(i1)) / 3._dp
     B1m2a = - 8._dp * B1(mf2(i2),mglu2, mSf2(i1)) / 3._dp
  !---------------------------------------------------------------
  ! note the somewhat misleading notation, because in coupL, coupR
  ! the chirality of the gluino is meant with L, R
  !---------------------------------------------------------------
     SumR(i2,i2) = Abs( c_L(i2) )**2 * B1m2a
     SumL(i2,i2) = Abs( c_R(i2) )**2 * B1m2a
     SumS(i2,i2) = Conjg( c_R(i2) ) * c_L(i2) * B0m2a
     Do i3=1,3
      If (i2.Ne.i3) Then
       SumR(i2,i3) = Conjg( c_L(i2) ) * c_L(i3) * B1m2
       SumL(i2,i3) = Conjg( c_R(i2) ) * c_R(i3) * B1m2
       SumS(i2,i3) = Conjg( c_R(i2) ) * c_L(i3) * B0m2
      End If
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "gluino squark quark",i1,i2
      Write(ErrCan,*) " sig_L ",sumL(i2,:)
      Write(ErrCan,*) " sig_R ",sumR(i2,:)
      Write(ErrCan,*) " sig_S ",sumS(i2,:)
     End If
    End Do
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End If
  !--------------
  ! gauge bosons
  !--------------
   Do i1=1,3
    B1m2 = - 0.5_dp * gSU2**2 * B1(p2,mfp2(i1), mW2)
    If (T3.Gt.0._dp) Then
     c_L = CKM(i1,:) 
    Else
     c_L = Conjg(CKM(:,i1)) 
    End If
    Do i2=1,3
     B1m2a = - 0.5_dp * gSU2**2 * B1(mf2(i2),mfp2(i1), mW2)
     Do i3=1,3
      If (i2.Eq.i3) Then
       sumL(i2,i3) = c_L(i3) * Conjg(c_L(i2) ) * B1m2a
      Else
       sumL(i2,i3) = c_L(i3) * Conjg(c_L(i2) ) * B1m2
      End If
     End Do
     If (WriteOut) Then
      Write (ErrCan,*) "W f f'",i1,i2
      Write(ErrCan,*) " sig_L ",sumL(i2,:)
     End If
    End Do
    SigL = SigL + sumL
   End Do

  Call CoupFermionZ(T3,e,gSU2,sinW2,coupL,coupR)
  sumL =0._dp
  sumR =0._dp
  sumS =0._dp
  Do i1=1,3
   B1m2 = - B1(mf2(i1),mf2(i1), mZ2)
   B0m2 = - 4._dp * mf(i1) * B0(mf2(i1),mf2(i1), mZ2)
   sumL(i1,i1) = coupL**2 * B1m2
   sumR(i1,i1) = coupR**2 * B1m2
   sumS(i1,i1) = coupL*coupR * B0m2
   If (WriteOut) Then
    Write(ErrCan,*) "Z-boson",i1
    Write(ErrCan,*) " sig_L ",sumL(i1,i1)
    Write(ErrCan,*) " sig_R ",sumR(i1,i1)
    Write(ErrCan,*) " sig_S ",sumS(i1,i1)
   End If
   SigL = SigL + sumL
   SigR = SigR + sumR
   SigS = SigS + sumS
  End Do
  !----------------------
  ! charged Higgs bosons
  !---------------------
  Do i1=1,n_Spm
   Do i2=1,3
    B1m2 = - 0.5_dp * B1(p2,mfp2(i2), mSpm2(i1))
    B0m2 = mfp(i2) * B0(p2,mfp2(i2), mSpm2(i1))
    Do i3=1,3
     If (T3.Gt.0._dp) Then
      Call CoupChargedScalarFermion(i1, i2, i3, RSpm, yuk_fp, id3C, id3C &
                                  & , yuk_f, id3C, id3C, c_L(i3), c_R(i3))
     Else
      Call CoupChargedScalarFermion(i1, i3, i2, RSpm, yuk_f, id3C, id3C &
                                  & , yuk_fp, id3C, id3C, c_R(i3), c_L(i3))
      c_R(i3) = Conjg(c_R(i3))
      c_L(i3) = Conjg(c_L(i3))
     End If
    End Do
    Do i3=1,3
     B1m2a = - 0.5_dp * B1(mf2(i3),mfp2(i2), mSpm2(i1))
     B0m2a = mfp(i2) * B0(mf2(i3),mfp2(i2), mSpm2(i1))
     SumL(i3,i3) = Abs(c_L(i3))**2 * B1m2a
     SumR(i3,i3) = Abs(c_R(i3))**2 * B1m2a
     SumS(i3,i3) = c_L(i3)*Conjg(c_R(i3)) * B0m2a
     Do i4=1,3
      If (i3.Ne.i4) Then
       SumL(i3,i4) = c_L(i4) * Conjg( c_L(i3) ) * B1m2
       SumR(i3,i4) = c_R(i4) * Conjg( c_R(i3) ) * B1m2
       SumS(i3,i4) = c_L(i4) * Conjg( c_R(i3) ) * B0m2
      End If
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "S+ f'",i1,i2,i3
      Write(ErrCan,*) "sumL",sumL(i3,:)
      Write(ErrCan,*) "sumR",sumR(i3,:)
      Write(ErrCan,*) "sumS",sumS(i3,:)
     End If
    End Do
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End Do 

  !-------------
  !  charginos
  !-------------
  n_fp = 6 ! number of sfermions'
  If (e.Eq.-1._dp) n_fp = 3  ! sneutrinos
  Do i1=1,n_char
   Do i2=1,n_fp
    B1m2 = - 0.5_dp * B1(p2,mC2(i1), mSf2p(i2))
    B0m2 = mC(i1) * B0(p2,mC2(i1), mSf2p(i2))
    Do i3=1,3
     If (T3.Gt.0._dp) Then
      Call CoupCharginoSfermion(i1, i3, i2, gSU2, T3, RSfp, Yuk_fp, Yuk_f &
                                &, id3C, id3C, U, V, c_L(i3), c_R(i3))
     Else
      Call CoupCharginoSfermion(i1, i3, i2, gSU2, T3, RSfp, Yuk_f, Yuk_fp &
                                &, id3C, id3C, U, V, c_L(i3), c_R(i3))
     End If
    End Do
    Do i3=1,3
     B1m2a = - 0.5_dp * B1(mf2(i3),mC2(i1), mSf2p(i2))
     B0m2a = mC(i1) * B0(mf2(i3),mC2(i1), mSf2p(i2))
  !---------------------------------------------------------------
  ! note the somewhat misleading notation, because in coupL, coupR
  ! the chirality of the chargino is meant with L, R
  !---------------------------------------------------------------
     SumR(i3,i3) = Abs(c_L(i3))**2 * B1m2a
     SumL(i3,i3) = Abs(c_R(i3))**2 * B1m2a
     SumS(i3,i3) = c_L(i3)*Conjg(c_R(i3)) * B0m2a
     Do i4=1,3
      If (i3.Ne.i4) Then
       SumR(i3,i4) = c_L(i4) * Conjg( c_L(i3) ) * B1m2
       SumL(i3,i4) = c_R(i4) * Conjg( c_R(i3) ) * B1m2
       SumS(i3,i4) = c_L(i4) * Conjg( c_R(i3) ) * B0m2
      End If
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "Chi Sf'",i1,i2,i3
      Write(ErrCan,*) "sumL",sumL(i3,:)
      Write(ErrCan,*) "sumR",sumR(i3,:)
      Write(ErrCan,*) "sumS",sumS(i3,:)
     End If
    End Do
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End Do

  !------------------------------
  ! neutral Higgs contributions
  !------------------------------
  sumL =0._dp
  sumR =0._dp
  sumS =0._dp
  Do i1=1,n_P0
   Do i2=1,3
    Call CoupFermionPseudoScalar(i1, T3, Y_f_in(i2, i2), RP0, c_L(i2), c_R(i2))

    B1m2 = - 0.5_dp * B1(mf2(i2),mf2(i2), mP02(i1))
    B0m2 = mf(i2) * B0(mf2(i2),mf2(i2), mP02(i1))
    SumR(i2,i2) = Abs(c_L(i2))**2 * B1m2
    SumL(i2,i2) = Abs(c_R(i2))**2 * B1m2
    SumS(i2,i2) = c_L(i2) * Conjg(c_R(i2)) * B0m2

    If (WriteOut) Then
     Write(ErrCan,*) "P0 f",i1,i2
     Write(ErrCan,*) "sumL",sumL(i2,i2)
     Write(ErrCan,*) "sumR",sumR(i2,i2)
     Write(ErrCan,*) "sumS",sumS(i2,i2)
    End If
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End Do

  sumL =0._dp
  sumR =0._dp
  sumS =0._dp
  Do i1=1,n_S0
   Do i2=1,3
    Call CoupFermionScalar(i1, T3, Y_f_in(i2, i2), RS0, c_L(i2), c_R(i2))

    B1m2 = - 0.5_dp * B1(mf2(i2),mf2(i2), mS02(i1))
    B0m2 = mf(i2) * B0(mf2(i2),mf2(i2), mS02(i1))

    SumR(i2,i2) = Abs(c_L(i2))**2 * B1m2
    SumL(i2,i2) = Abs(c_R(i2))**2 * B1m2
    SumS(i2,i2) = c_L(i2) * Conjg(c_R(i2)) * B0m2

    If (WriteOut) Then
     Write(ErrCan,*) "S0 f",i1,i2
     Write(ErrCan,*) "sumL",sumL(i2,i2)
     Write(ErrCan,*) "sumR",sumR(i2,i2)
     Write(ErrCan,*) "sumS",sumS(i2,i2)
    End If
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End Do

  !-------------
  ! neutralinos
  !-------------
  Do i1=1,n_neut
   Do i2=1,6
    B1m2 = - 0.5_dp * B1(p2,mN2(i1), mSf2(i2))
    B0m2 = mN(i1) * B0(p2,mN2(i1), mSf2(i2))
    Do i3=1,3
     Call CoupNeutralinoSfermion(i3, i1, i2, gp, gSU2, T3, e, RSf, id3C &
                               &, id3C, Y_f_in, N, c_L(i3), c_R(i3))
    End Do
    Do i3=1,3
     B1m2a = - 0.5_dp * B1(mf2(i3),mN2(i1), mSf2(i2))
     B0m2a = mN(i1) * B0(mf2(i3),mN2(i1), mSf2(i2))
  !---------------------------------------------------------------
  ! note the somewhat misleading notation, because in coupL, coupR
  ! the chirality of the neutralino is meant with L, R
  !---------------------------------------------------------------
     SumR(i3,i3) = Abs(c_L(i3))**2 * B1m2a 
     SumL(i3,i3) = Abs(c_R(i3))**2 * B1m2a
     SumS(i3,i3) = c_L(i3)*Conjg(c_R(i3)) * B0m2a
     Do i4=1,3
      If (i3.Ne.i4) Then
       SumR(i3,i4) = c_L(i4) * Conjg( c_L(i3) ) * B1m2
       SumL(i3,i4) = c_R(i4) * Conjg( c_R(i3) ) * B1m2
       SumS(i3,i4) = c_L(i4) * Conjg( c_R(i3) ) * B0m2
      End If
     End Do
     If (WriteOut) Then
      Write(ErrCan,*) "Chi0 Sf",i1,i2,i3
      Write(ErrCan,*) "sumL",sumL(i3,:)
      Write(ErrCan,*) "sumR",sumR(i3,:)
      Write(ErrCan,*) "sumS",sumS(i3,:)
     End If
    End Do
    SigL = SigL + sumL
    SigR = SigR + sumR
    SigS = SigS + sumS
   End Do
  End Do

  SigL = oo16pi2 * SigL
  SigR = oo16pi2 * SigR
  SigS = oo16pi2 * SigS

  Iname = Iname - 1

 End Subroutine Sigma_Fermion3sckm


 Subroutine Sigma_Fermion(p2, i_gen, mf, yuk_f, RfL, RfR, gSU2, gSU3          &
          & , sinW2 ,T3, e, mfp, yuk_fp, RfpL, RfpR, mSf2, RSf, mSf2p, Rsfp   &
          & , mglu , phase_glu, mN, mN2, N, mC, mC2, U, V, mS02, RS0, mP02    &
          & , RP0, mSpm2 , RSpm, mZ2, mW2, l_SU3, Resummed, res)
 !-----------------------------------------------------------------------
 ! Calculates the 1-loop self energy of SM-fermion. 
 ! The formula of J. Bagger et al, Nucl.Phys.B is used. The renormalization
 ! scale has to be set in the main program due to the structure of LoopTools. 
 ! the input is:
 ! - p2 ......... the outer momentum squared 
 ! - gSU3 ....... the SU(3) gauge coupling at p2
 ! - gSU2 ....... the SU(2) gauge coupling at p2
 ! - sinW2 ...... sin(theta_W) squared
 ! - T3 ......... weak isospin of the fermion
 ! - e .......... charge of the fermion
 ! - mf ......... fermion mass
 ! - mSf(i) ..... Sfermion masses
 ! - RSf(i,j) ... mixing matrix of the sfermions
 ! - mfp ........ fermion mass'
 ! - mSfp(i) .... Sfermion masses'
 ! - RSfp(i,j) .. mixing matrix of the sfermions'
 ! - mglu ....... gluino mass
 ! - phase_glu .. phase of the gluino 
 ! - mZ2 ........ Z-boson mass
 ! - mW2 ........ W-boson mass
 ! - mS02 ....... masses of the scalar Higgs
 ! - RS0 ........ mixing matrix of the scalar Higgs
 ! - mP02 ....... masses of the pseudoscalar Higgs
 ! - RP0 ........ mixing matrix of the pseudoscalar Higgs
 ! - mSpm2 ...... masses of the charged Higgs
 ! - RSpm ....... mixing matrix of the charged Higgs
 ! - mN ......... masses of the neutralinos
 ! - N .......... mixing matrix of the neutrlinos
 ! - mC ......... masses of the charginos
 ! - U,V ........ mixing matrices of the charginos
 !    
 ! written by Werner Porod, 25.12.01
 ! 30.01.03: correcting bug in 2-loop gluonic QCD corrections, forgot
 !           the 2-loop logs
 ! 05.01.04: taking gauge boson masses as input to allow for running masses
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_gen
  Real(dp), Intent(in) :: p2, gSU3, gSU2, sinW2, mS02(:), RS0(:,:), mP02(:) &
     & ,RP0(:,:), mSpm2(:), mN(:), mN2(:), mC(:), mC2(:), mf(3), mSf2(6)    &
     & ,mfp(3), mSf2p(:), e, T3, mglu, mZ2, mW2
  Complex(dp), Intent(in) :: N(:,:), U(:,:), V(:,:), RSpm(:,:), Rsf(6,6)  &
     & ,phase_glu, Rsfp(:,:), yuk_f(3,3), yuk_fp(3,3), RfL(3,3), RfR(3,3) &
     & , RfpL(3,3), RfpR(3,3)
  Logical, Intent(in) :: Resummed, l_SU3
  Complex(dp), Intent(out) :: res

  Integer :: i1, i2
  Real(dp) :: cosW2, coupL, coupR, mf2, mfp2, mglu2, gp, m22, f_SU3, Q2, logQ
  Complex(dp) :: sumI, coupLC, coupRC, yukD, yukU, Rsf2(2,2), Rsfp2(2,2)
  Logical :: WriteOut
  Real(dp), Parameter :: &
    & as2loop = 1._dp / 24._dp + 2011._dp * oo32Pi2 / 12._dp           &
    &         + Log2 / 12._dp - oo8Pi2 * Zeta3                        &
    & , log2loop_a = 123._dp * oo32Pi2, log2loop_b = 33._dp * oo32Pi2
 !----------------
 ! Initialization
 !----------------
  Iname = Iname + 1
  NameOfUnit(Iname) ="Sigma_Fermion"
  If ((WriteOneLoopContributions.Eq.10).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions in Sigma_Fermion:",mf(i_gen),mfp(i_gen)
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If

  Q2 = GetRenormalizationScale() ! from LoopTools
  cosW2 = 1._dp - sinW2
  mf2 = mf(i_gen)**2
  mfp2 = mfp(i_gen)**2
  mglu2 = mglu**2
  gp = gSU2 * Sqrt(sinW2/cosW2)
  res = ZeroC
  If (.Not.GenerationMixing) Then
   If (e.Eq.-1._dp) Then
    Rsfp2 = id2C
   Else
    RSfp2 = Rsfp(2*i_gen-1:2*i_gen, 2*i_gen-1:2*i_gen)
   End If
   RSf2 = Rsf(2*i_gen-1:2*i_gen, 2*i_gen-1:2*i_gen)
  End If
  !------------------------
  ! strong interactions
  !------------------------
  If (l_SU3) Then
   f_SU3 = 4._dp * gSU3**2 / 3._dp
   If (Resummed) Then ! for small masses a resummation has be done before
    sumI = 0._dp
   Else
    logQ = Log(Q2/mf2)
    sumI = - f_SU3 * mf(i_gen) * (5._dp + 3._dp * LogQ               &
         &                       + (as2loop + log2loop_a * logQ      &
         &                         + log2loop_b * logQ**2) * gSU3**2 )
   End If
   If (WriteOut) Write(ErrCan,*) "gluon :",sumI,gsu3
   res = res + sumI

#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Do i1=1,6
     Call CoupGluinoSquark(gSU3, phase_glu, i_gen, i1, RSf, RfL, RfR &
                         &, coupLC, coupRC)
     SumI = - mf(i_gen) * (Abs(coupLC)**2 + Abs(coupRC)**2)                  &
          &             * B1(p2,mglu2, mSf2(i1))                             &
          & + 2._dp * mglu * coupLC * Conjg( coupRC ) * B0(p2,mglu2, mSf2(i1))
     sumI = 8._dp * sumI / 3._dp 
     If (WriteOut) Write(ErrCan,*) "gluino squark",i1,sumI
     res = res + sumI
    End Do
   Else
#endif
    Do i1=2*i_gen-1,2*i_gen
     SumI = - mf(i_gen) * B1(p2,mglu2, mSf2(i1))                             &
          & - 2._dp * mglu * Real( phase_glu * Conjg( Rsf(i1,2*i_gen-1) )    &
          &        * Rsf(i1,2*i_gen),dp ) * B0(p2,mglu2, mSf2(i1))
     sumI = f_SU3 * sumI 
     If (WriteOut) Write(ErrCan,*) "gluino squark",i1,i_gen,sumI
    res = res + sumI
    End Do
   End If
#ifdef GENERATIONMIXING
  End If
#endif
  !--------------
  ! gauge bosons
  !--------------
  sumI = - 0.5_dp * gSU2**2 * mf(i_gen) * B1(p2,mfp2,mW2)
  If (WriteOut) Write(ErrCan,*) "W-boson",sumI
  res = res + sumI
  Call CoupFermionZ(T3,e,gSU2,sinW2,coupL,coupR)
  sumI = - mf(i_gen) * ( (coupL**2 + coupR**2) * B1(p2,mf2,mZ2)  &
       &               + 4._dp * coupL * coupR * B0(p2,mf2,mZ2) )
  If (WriteOut) Write(ErrCan,*) "Z-boson",sumI
  res = res + sumI
  If (ReSummed) Then ! for smaller masses a resummation should be done before
   sumI = 0._dp
  Else
   sumI = - mf(i_gen) * (e*gSU2)**2 * sinW2 * (5._dp + 3._dp * Log(Q2/mf2) )
  End If
  If (WriteOut) Write(ErrCan,*) "Photon",sumI
  res = res + sumI
  !----------------------
  ! charged Higgs bosons
  !---------------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_Spm
    Do i2=1,3
     If (T3.Gt.0._dp) Then
      Call CoupChargedScalarFermion(i1, i2, i_gen, RSpm, yuk_fp, RfpL, RfpR &
                                  & , yuk_f, RfL, RfR, coupLC, coupRC)
     Else
      Call CoupChargedScalarFermion(i1, i_gen, i2, RSpm, yuk_f, RfL, RfR &
                                  & , yuk_fp, RfpL, RfpR, coupLC, coupRC)
     End If
     sumI = - 0.5_dp * mf(i_gen) * ( Abs(coupLC)**2 + Abs(coupRC)**2 )    &
        &                        * B1(p2,mfp2, mSpm2(i1))                 &
        & +  mfp(i2) * Conjg(coupLC) * coupRC * B0(p2, mfp(i2)**2, mSpm2(i1))
     If (WriteOut) Write(ErrCan,*) "H+ F'",i1,i2,sumI
     res = res + sumI
    End Do
   End Do
 
  Else
#endif
   If (T3.Gt.0._dp) Then
    yukD = yuk_fp(i_gen,i_gen)
    yukU = yuk_f(i_gen,i_gen)
   Else
    yukD = yuk_f(i_gen,i_gen)
    yukU = yuk_fp(i_gen,i_gen)
   End If
   Do i1=1,n_Spm
    Call CoupChargedScalarFermion(i1, RSpm, yukD, yukU, coupLC, coupRC)
    sumI = - 0.5_dp * mf(i_gen) * ( Abs(coupLC)**2 + Abs(coupRC)**2 )       &
       &                        * B1(p2,mfp2, mSpm2(i1))                    &
       &   +  mfp(i_gen) * Conjg(coupLC) * coupRC * B0(p2,mfp2, mSpm2(i1))
    If (WriteOut) Write(ErrCan,*) "H+",i1,sumI
    res = res + sumI
   End Do
#ifdef GENERATIONMIXING
 End If
#endif
  !-------------
  !  charginos
  !-------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_char
    If (e.Eq.-1._dp) Then
     Do i2=1,3
      Call CoupCharginoSfermion(i1, i_gen, i2, gSU2, T3, RSfp, Yuk_f, Yuk_fp &
                               &, RfL, RfR, U, V, coupLC, coupRC)
      sumI = - 0.5_dp* mf(i_gen) * (Abs(coupLC)**2 + Abs(coupRC)**2)         &
           &                      * B1(p2,mC2(i1), mSf2p(i2))                &
           & + Conjg(coupLC) * coupRC * mC(i1) * B0(p2, mC2(i1), mSf2p(i2))
      If (WriteOut) Write(ErrCan,*) "Chi+ Sf'",i1,i2,sumI
      res = res + sumI
     End Do
    Else 
     Do i2=1,6
      If (T3.Gt.0._dp) Then
       Call CoupCharginoSfermion(i1, i_gen, i2, gSU2, T3, RSfp, Yuk_fp, Yuk_f &
                                &, RfL, RfR, U, V, coupLC, coupRC)
      Else
       Call CoupCharginoSfermion(i1, i_gen, i2, gSU2, T3, RSfp, Yuk_f, Yuk_fp &
                                &, RfL, RfR, U, V, coupLC, coupRC)
      End If
      sumI = - 0.5_dp * mf(i_gen)* (Abs(coupLC)**2 + Abs(coupRC)**2)         &
           &                     * B1(p2,  mC2(i1), mSf2p(i2) )              &
           & + Conjg(coupRC) * coupLC * mC(i1) * B0(p2, mC2(i1), mSf2p(i2) )
      If (WriteOut) Write(ErrCan,*) "Chi+ Sf'",i1,i2,sumI
      res = res + sumI
     End Do
    End If
   End Do

  Else ! .not.GenerationMixing
#endif
   If (T3.Gt.0._dp) Then
    yukD = yuk_fp(i_gen,i_gen)
    yukU = yuk_f(i_gen,i_gen)
   Else
    yukD = yuk_f(i_gen,i_gen)
    yukU = yuk_fp(i_gen,i_gen)
   End If
   Do i1=1,n_char
    If (e.Eq.-1._dp) Then
     Call CoupCharginoSfermion(i1, 1, gSU2, T3,RSfp2, yukD, yukU, U, V  &
                              &, coupLC, coupRC)
     sumI = - 0.5_dp* mf(i_gen) * (Abs(coupLC)**2 + Abs(coupRC)**2)         &
          &                     * B1(p2,mC2(i1), mSf2p(i_gen))              &
          & + Conjg(coupRC) * coupLC * mC(i1) * B0(p2, mC2(i1), mSf2p(i_gen))
     If (WriteOut) Write(ErrCan,*) "Chi+ Sf'",i1,1,sumI
     res = res + sumI
    Else 
     Do i2=1,2
      m22 = mSf2p(2*(i_gen-1)+i2)
      Call CoupCharginoSfermion(i1, i2, gSU2, T3, RSfp2, yukD, yukU, U, V  &
                               &, coupLC, coupRC)
      sumI = - 0.5_dp * mf(i_gen)* (Abs(coupLC)**2 + Abs(coupRC)**2)  &
           &                     * B1(p2,mC2(i1),m22)                 &
           & + Conjg(coupRC) * coupLC * mC(i1) * B0(p2,mC2(i1), m22)
      If (WriteOut) Write(ErrCan,*) "Chi+ Sf'",i1,i2,sumI
      res = res + sumI
     End Do
    End If
   End Do
#ifdef GENERATIONMIXING
  End If
#endif
  !------------------------------
  ! neutral Higgs contributions
  !------------------------------
  Do i1=1,n_P0
#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Call CoupFermionPseudoScalar(i_gen, i_gen, i1, T3, yuk_f, Rfl, RfR, RP0 &
                                &, coupLC, coupRC)
   Else
#endif
    Call CoupFermionPseudoScalar(i1, T3, yuk_f(i_gen, i_gen), RP0 &
                                &, coupLC, coupRC)
#ifdef GENERATIONMIXING
   End If
#endif
   sumI = - mf(i_gen) * Abs(coupLC)**2 &
        &             * ( B1(p2,mf2,mP02(i1)) + B0(p2,mf2,mP02(i1)) )
   If (WriteOut) Write(ErrCan,*) "P0'",i1,sumI
   res = res + sumI
  End Do

  Do i1=1,n_S0
#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Call CoupFermionScalar(i_gen, i_gen, i1, T3, yuk_f, Rfl, RfR, RS0 &
                          &, coupLC, coupRC)
   Else
#endif
    Call CoupFermionScalar(i1, T3, yuk_f(i_gen, i_gen), RS0, coupLC, coupRC)
#ifdef GENERATIONMIXING
   End If
#endif
   sumI = - mf(i_gen) * Abs(coupLC)**2  &
        &             * (B1(p2,mf2,mS02(i1)) - B0(p2,mf2,mS02(i1)) )
   If (WriteOut) Write(ErrCan,*) "S0'",i1,sumI
   res = res + sumI
  End Do
  !-------------
  ! neutralinos
  !-------------
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_neut
    Do i2=1,6
     Call CoupNeutralinoSfermion(i_gen, i1, i2, gp, gSU2, T3, e, RSf, RfL &
                               &, RfR, Yuk_f, N, coupLC, coupRC)
     sumI = - 0.5_dp * mf(i_gen) * (Abs(coupLC)**2 + Abs(coupRC)**2)  &
          &                      * B1(p2, mN2(i1), mSf2(i2))          &
          & + Conjg(coupRC) * coupLC * mN(i1) * B0(p2, mN2(i1), mSf2(i2))
     If (WriteOut) Write(ErrCan,*) "Chi0 Sf",i1,i2,sumI
     res = res + sumI
    End Do
   End Do

  Else  ! .not.GenerationMixing
#endif
   Do i1=1,n_neut
    Do i2=1,2
     m22 = mSf2(2*(i_gen-1)+i2)
     Call CoupNeutralinoSfermion(i1, i2, gp, gSU2, T3, e, RSf2 &
                                &, yuk_f(i_gen,i_gen), N, coupLC, coupRC)
     sumI = - 0.5_dp * mf(i_gen) * (Abs(coupLC)**2 + Abs(coupRC)**2)  &
          &                      * B1(p2, mN2(i1), m22)               &
          & + Conjg(coupRC) * coupLC * mN(i1) * B0(p2, mN2(i1), m22)
     If (WriteOut) Write(ErrCan,*) "Chi0 Sf",i1,2*(i_gen-1)+i2,sumI
     res = res + sumI
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  res = oo16pi2 * res

  Iname = Iname - 1
 End Subroutine Sigma_Fermion


#ifdef GENERATIONMIXING
 Subroutine Sigma_Gluino(p2, Q2, gSU3, mglu, phase_glu, mSup2, RSup, uU_L  &
                     & , uU_R, mSdown2, RSdown, uD_L, uD_R, res)
#else
 Subroutine Sigma_Gluino(p2, Q2, gSU3, mglu, phase_glu, mSup2, RSup  &
                     & , mSdown2, RSdown, res)
#endif
 !-----------------------------------------------------------------------
 ! Calculates the 1-loop self energy of SM-fermion. 
 ! The formula of J. Bagger et al, Nucl.Phys.B is used. The renormalization
 ! scale has to be set in the main program due to the structure of LoopTools. 
 ! the input is:
 ! - p2 ........... the outer momentum squared 
 ! - Q2 ........... renormalization scale
 ! - mglu ......... gluino mass
 ! - phase_glu .... phase of the gluino 
 ! - gSU3 ......... the SU(3) gauge coupling at p2
 ! - mSup2(i) ..... u-squark masses squared
 ! - RSup(i,j) .... u-squark mixing matrix
 ! - mSdown(i) .... d-squark masses squared
 ! - RSdown(i,j) .. d-squark mixing matrix
 !    
 ! written by Werner Porod, 13.8.99
 ! 27.11.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: p2, Q2, gSU3, mglu, mSup2(6), mSdown2(6)
  Complex(dp), Intent(in) :: RSup(6,6), RSdown(6,6), phase_glu
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: uU_L(3,3), uU_R(3,3), uD_L(3,3), uD_R(3,3)
#endif
  Complex(dp), Intent(out) :: res

  Integer :: i1, i2, i_sq
  Real(dp) :: mglu2
  Complex(dp) :: sumI
  Logical :: WriteOut
#ifdef GENERATIONMIXING
  Complex(dp) :: coupL, coupR
#endif

  Iname = Iname + 1
  NameOfUnit(Iname) = "Sigma_Gluino"

  WriteOut = .False.
  If ( (WriteOneLoopContributions.Eq.8).Or.(WriteOneLoopContributions.Lt.0)) &
   &  WriteOut = .True.

  If (WriteOut) Write(ErrCan,*) "Contributions to Sigma_Gluino"

  mglu2 = mglu**2

  If (mglu2.Eq.p2) Then
   sumI = - mglu * (15._dp + 9._dp * Log(Q2/mglu2) )
  Else
   sumI = - 6._dp * mglu   &
        &         * Real(B1(p2, mglu2, 0._dp) + 2._dp * B0(p2, mglu2, 0._dp), dp )
  End If

  If (WriteOut) Write(ErrCan,*) "Gluon ",sumI
  res = sumI

#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,6
    sumI = 0._dp
    Do i2=1,3
     Call CoupGluinoSquark(1._dp, phase_glu, i2, i1, RSup, uU_L, uU_R &
                         &, coupL, coupR)
     sumI = sumI                                             &
        & - 2._dp * mglu * (Abs(coupL)**2 + Abs(coupR)**2)   &
        &        * Real(B1(p2, mf_u2(i2), mSup2(i1)),dp)     &
        & - 4._dp * mf_u(i2) * coupL * Conjg(coupR)          &
        &         * Real(B0(p2, mf_u2(i2), mSup2(i1)),dp)
    End Do
    If (WriteOut) Write(ErrCan,*) "U-Squark",i1,sumI
    res = res + sumI

    sumI = 0._dp
    Do i2=1,3
     Call CoupGluinoSquark(1._dp, phase_glu, i2, i1, RSdown, uD_L, uD_R &
                         &, coupL, coupR)
     sumI = sumI                                            &
        & - 2._dp * mglu * (Abs(coupL)**2 + Abs(coupR)**2)  &
        &        * Real(B1(p2, mf_d2(i2), mSdown2(i1)),dp)  &
        & - 4._dp * mf_d(i2) * coupL * Conjg(coupR)         &
        &         * Real(B0(p2, mf_d2(i2), mSdown2(i1)),dp) 
    End Do
    If (WriteOut) Write(ErrCan,*) "D-Squark",i1,sumI
    res = res + sumI

   End Do
  Else ! .not.GenerationMixing
#endif
   sumI = 0._dp
   Do i1=1,3
    i_sq = (i1-1)*2
    Do i2=1,2
     sumI = sumI  - Real(B1(p2, mf_u2(i1), mSup2( i_sq+i2)),dp)
     sumI = sumI  - Real(B1(p2, mf_d2(i1), mSdown2( i_sq+i2)),dp)
    End Do
   End Do
   If (WriteOut) Write(ErrCan,*) "Squark, diagonal:",mglu*sumI
   res = res + sumI * mglu

   Do i1=1,3
    sumI = 0._dp
    i_sq = (i1-1)*2
    Do i2=1,2
     sumI = sumI + 2._dp * Real( phase_glu * Conjg( RSup(i_sq+i2,i_sq+1) )    &
          &                    * RSup(i_sq+i2,i_sq+2),dp )                    &
          &              * Real(B0(p2, mf_u2(i1), mSup2(i_sq+i2)),dp)
    End Do
    If (WriteOut) Write(ErrCan,*) "u-Squark, off-diagonal:",i1,mf_u(i1)*sumI
    res = res + mf_u(i1)*sumI
  
    sumI = 0._dp
    Do i2=1,2
     sumI = sumI + 2._dp * Real( phase_glu * Conjg( RSdown(i_sq+i2,i_sq+1) )  &
       &                       * RSdown(i_sq+i2,i_sq+2),dp )                  &
       &                 * Real(B0(p2, mf_d2(i1), mSdown2(i_sq+i2)),dp)
    End Do
    If (WriteOut) Write(ErrCan,*) "d-Squark, off-diagonal:",i1,mf_d(i1)*sumI
    res = res + mf_d(i1)*sumI
   End Do
#ifdef GENERATIONMIXING
  End If ! GenerationMixing
#endif

  If (WriteOut) Write(ErrCan,*) "total ",res

  res = oo16pi2 * gSU3**2 * res
  If (WriteOut) Write(ErrCan,*) "m_g - dm_g ",mglu,Real(mglu-res)

  Iname = Iname - 1

 End Subroutine Sigma_Gluino

 Subroutine Sigma_Neutralino(p2, mN, mN2, c_NNZ_L, c_NNZ_R, mS02, c_NNS0_L    &
            & , c_NNS0_R, mP02, c_NNP0_L, c_NNP0_R, mC, mC2, c_CNW_L, c_CNW_R &
            & , mSpm2, c_SmpCN_L, c_SmpCN_R, mSup2, c_UNSu_L, c_UNSu_R        &
            & , mSdown2, c_DNSd_L, c_DNSd_R, mSlepton2, c_LNSl_L, c_LNSl_R    &
            & , mSneut2, c_NuNSn_L, c_NuNSn_R, mZ2, mW2                       &
            & , WriteOut, SigL, SigR, SigS)
 !-----------------------------------------------------------------------
 ! written by Werner Porod, 27.12.01
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: p2, mN(:), mN2(:), mC(:), mC2(:), mSup2(:)        &
      &  , mSdown2(:), mSlepton2(:), mSneut2(:), mS02(:), mP02(:), mSpm2(:) &
      &  , mZ2, mW2
  Complex(dp), Intent(in) :: c_NNZ_L(:,:), c_NNZ_R(:,:)                     &
     & , c_NNS0_L(:,:,:), c_NNS0_R(:,:,:), c_NNP0_L(:,:,:), c_NNP0_R(:,:,:) &
     & , c_UNSu_L(:,:,:), c_UNSu_R(:,:,:), c_LNSl_L(:,:,:), c_LNSl_R(:,:,:) &
     & , c_NuNSn_L(:,:,:), c_NuNSn_R(:,:,:), c_CNW_L(:,:), c_CNW_R(:,:)     &
     & , c_SmpCN_L(:,:,:), c_SmpCN_R(:,:,:), c_DNSd_L(:,:,:), c_DNSd_R(:,:,:)
  Logical, Intent(in) :: WriteOut
  Complex(dp), Intent(out) :: SigL(:,:), SigR(:,:), SigS(:,:)

  Real(dp) :: B1m2, B0m2
  Integer :: i1, i2, i3, i4, i2_min, i2_max
  Complex(dp) :: sumL(n_neut,n_neut), sumR(n_neut,n_neut), sumS(n_neut,n_neut)

  !------------------
  ! Inititalisation
  !------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "Sigma_Neutralino"
  SigL = ZeroC
  SigR = ZeroC
  SigS = ZeroC
  !---------------------
  ! quark - squark
  !---------------------
  Do i1=1,3
   If (GenerationMixing) Then
    i2_min = 1
    i2_max = 6
   Else
    i2_min = 2*i1-1
    i2_max = 2*i1
   End If
   Do i2=i2_min,i2_max
    B1m2 = - 3._dp * B1(p2, mf_u2(i1), mSup2(i2))
    B0m2 = 6._dp * mf_u(i1) * B0(p2, mf_u2(i1), mSup2(i2))
    Do i3=1,n_neut
     Do i4=1,n_neut
      SumL(i3,i4) = Conjg( c_UNSu_L(i1,i3,i2) ) * c_UNSu_L(i1,i4,i2) * B1m2
      SumR(i3,i4) = Conjg( c_UNSu_R(i1,i3,i2) ) * c_UNSu_R(i1,i4,i2) * B1m2
      SumS(i3,i4) = Conjg( c_UNSu_R(i1,i3,i2) ) * c_UNSu_L(i1,i4,i2) * B0m2
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "u-quark, u-squark",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigR(i,:)",i3,SumR(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
     End Do 
    End If
    SigL = SigL + SumL
    SigR = SigR + SumR
    SigS = SigS + SumS

    B1m2 = - 3._dp * B1(p2, mf_d2(i1), mSdown2(i2))
    B0m2 = 6._dp * mf_d(i1) * B0(p2, mf_d2(i1), mSdown2(i2))
    Do i3=1,n_neut
     Do i4=1,n_neut
      SumL(i3,i4) = Conjg( c_DNSd_L(i1,i3,i2) ) * c_DNSd_L(i1,i4,i2) * B1m2
      SumR(i3,i4) = Conjg( c_DNSd_R(i1,i3,i2) ) * c_DNSd_R(i1,i4,i2) * B1m2
      SumS(i3,i4) = Conjg( c_DNSd_R(i1,i3,i2) ) * c_DNSd_L(i1,i4,i2) * B0m2
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "d-quark, d-squark",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigR(i,:)",i3,SumR(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
     End Do 
    End If
    SigL = SigL + SumL
    SigR = SigR + SumR
    SigS = SigS + SumS
   End Do
  End Do

 !-----------------
 ! lepton-sleptons 
 !-----------------
  Do i1=1,5-n_char ! needed for R-parity violation
   If (GenerationMixing) Then
    i2_min = 1
    i2_max = 2*(5-n_char)
   Else
    i2_min = 2*i1-1
    i2_max = 2*i1
   End If
   Do i2=i2_min,i2_max
    B1m2 = - B1(p2, mf_l2(i1), mSlepton2(i2))
    B0m2 = 2._dp * mf_l(i1) * B0(p2, mf_l2(i1), mSlepton2(i2))
    Do i3=1,n_neut
     Do i4=1,n_neut
      SumL(i3,i4) = Conjg( c_LNSl_L(i1,i3,i2) ) * c_LNSl_L(i1,i4,i2) * B1m2
      SumR(i3,i4) = Conjg( c_LNSl_R(i1,i3,i2) ) * c_LNSl_R(i1,i4,i2) * B1m2
      SumS(i3,i4) = Conjg( c_LNSl_R(i1,i3,i2) ) * c_LNSl_L(i1,i4,i2) * B0m2
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "lepton, slepton",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigR(i,:)",i3,SumR(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
     End Do 
    End If
    SigL = SigL + SumL
    SigR = SigR + SumR
    SigS = SigS + SumS
   End Do
  End Do

 !--------------------
 ! neutrino-sneutrino
 !--------------------
  Do i1=1,5-n_char ! needed for R-parity violation
   If (GenerationMixing) Then
    i2_min = 1
    i2_max = 5-n_char
   Else
    i2_min = i1
    i2_max = i1
   End If
   Do i2=i2_min,i2_max
    B1m2 = - B1(p2, 0._dp, mSneut2(i2))
    Do i3=1,n_neut
     Do i4=1,n_neut
      SumL(i3,i4) = Conjg( c_NuNSn_L(i1,i3,i2) ) * c_NuNSn_L(i1,i4,i2) * B1m2
      SumR(i3,i4) = Conjg( c_NuNSn_R(i1,i3,i2) ) * c_NuNSn_R(i1,i4,i2) * B1m2
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "neutrino, sneutrino",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigR(i,:)",i3,SumR(i3,:)
     End Do 
    End If
    SigL = SigL + SumL
    SigR = SigR + SumR
   End Do
  End Do


 !-------------------------------------
 ! chargino - W
 !-------------------------------------
  Do i1 = 1,n_char
   B1m2 = - 2._dp * B1(p2,mC2(i1),mW2)
   B0m2 = - 8._dp * mC(i1) * B0(p2,mC2(i1),mW2)
   Do i2 = 1, n_neut
    Do i3 = 1,n_neut
     SumL(i2,i3) = Conjg( c_CNW_R(i1,i2) ) * c_CNW_R(i1,i3) * B1m2
     SumR(i2,i3) = Conjg( c_CNW_L(i1,i2) ) * c_CNW_L(i1,i3) * B1m2
     SumS(i2,i3) = Conjg( c_CNW_L(i1,i2) ) * c_CNW_R(i1,i3) * B0m2
    End Do 
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "chargino, W",i1
    Do i3=1,n_neut
     Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
    End Do 
    Do i3=1,n_neut
     Write(ErrCan,*) "SigR(i,:)",i3,SumR(i3,:)
    End Do 
    Do i3=1,n_neut
     Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
    End Do 
   End If
   SigL = SigL + SumL
   SigR = SigR + SumR
   SigS = SigS + SumS
  End Do
 !-------------------------------------
 ! chargino - S+
 !-------------------------------------
  Do i1 = 1,n_char
   Do i2 = 1,n_Spm
    B1m2 = - B1(p2,mC2(i1),mSpm2(i2))
    B0m2 = 2._dp * mC(i1) * B0(p2,mC2(i1),mSpm2(i2))
    Do i3 = 1,n_neut
     Do i4 = 1,n_neut
      SumL(i3,i4) = Conjg( c_SmpCN_R(i2,i1,i3) ) * c_SmpCN_R(i2,i1,i4) * B1m2
      SumR(i3,i4) = Conjg( c_SmpCN_L(i2,i1,i3) ) * c_SmpCN_L(i2,i1,i4) * B1m2
      SumS(i3,i4) = Conjg( c_SmpCN_L(i2,i1,i3) ) * c_SmpCN_R(i2,i1,i4) * B0m2
     End Do 
    End Do 
    If (WriteOut) Then
     Write(ErrCan,*) "Chargino, S+",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigR(i,:)",i3,SumR(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
     End Do 
    End If
    SigL = SigL + SumL
    SigR = SigR + SumR
    SigS = SigS + SumS
   End Do
  End Do
 !-------------------------------------
 ! neutral gauge bosons
 !-------------------------------------
  Do i1 = 1,n_neut
   B1m2 = - B1(p2,mN2(i1),mZ2) 
   B0m2 = - 4._dp * mN(i1) * B0(p2,mN2(i1),mZ2)
   Do i2 = 1,n_neut
    Do i3 = 1,n_neut
     SumL(i2,i3) = Conjg( c_NNZ_R(i2,i1) ) * c_NNZ_R(i3,i1) * B1m2
     SumR(i2,i3) = Conjg( c_NNZ_L(i2,i1) ) * c_NNZ_L(i3,i1) * B1m2
     SumS(i2,i3) = Conjg( c_NNZ_L(i2,i1) ) * c_NNZ_R(i3,i1) * B0m2
    End Do 
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "neutralino, Z",i1
    Do i3=1,n_neut
     Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
    End Do 
    Do i3=1,n_neut
     Write(ErrCan,*) "SigR(i,:)",i3,SumR(i3,:)
    End Do 
    Do i3=1,n_neut
     Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
    End Do 
   End If
   SigL = SigL + SumL
   SigR = SigR + SumR
   SigS = SigS + SumS
  End Do

 !-------------------------------------
 ! neutralino - P0
 !-------------------------------------
  Do i1 = 1,n_neut
   Do i2 = 1,n_P0
    B1m2 = - 0.5_dp * B1(p2,mN2(i1),mP02(i2))
    B0m2 = mN(i1) * B0(p2,mN2(i1),mP02(i2))
    Do i3 = 1,n_neut
     Do i4 = 1,n_neut
      SumL(i3,i4) = Conjg( c_NNP0_L(i3,i1,i2) ) * c_NNP0_L(i4,i1,i2) * B1m2
      SumR(i3,i4) = Conjg( c_NNP0_R(i3,i1,i2) ) * c_NNP0_R(i4,i1,i2) * B1m2
      SumS(i3,i4) = Conjg( c_NNP0_R(i3,i1,i2) ) * c_NNP0_L(i4,i1,i2) * B0m2
     End Do 
    End Do 
    If (WriteOut) Then
     Write(ErrCan,*) "Neutralino, P0",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigR(i,:)",i3,SumR(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
     End Do 
    End If
    SigL = SigL + SumL
    SigR = SigR + SumR
    SigS = SigS + SumS
   End Do
  End Do

 !-------------------------------------
 ! neutralino - S0
 !-------------------------------------
  Do i1 = 1,n_neut
   Do i2 = 1,n_S0
    B1m2 = - 0.5_dp * B1(p2,mN2(i1),mS02(i2))
    B0m2 = mN(i1) * B0(p2,mN2(i1),mS02(i2))
    Do i3 = 1,n_neut
     Do i4 = 1,n_neut
      SumL(i3,i4) = Conjg( c_NNS0_L(i3,i1,i2) ) * c_NNS0_L(i4,i1,i2) * B1m2
      SumR(i3,i4) = Conjg( c_NNS0_R(i3,i1,i2) ) * c_NNS0_R(i4,i1,i2) * B1m2
      SumS(i3,i4) = Conjg( c_NNS0_R(i3,i1,i2) ) * c_NNS0_L(i4,i1,i2) * B0m2
     End Do 
    End Do 
    If (WriteOut) Then
     Write(ErrCan,*) "Neutralino, S0",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigR(i,:)",i3,SumR(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
     End Do 
    End If
    SigL = SigL + SumL
    SigR = SigR + SumR
    SigS = SigS + SumS
   End Do
  End Do

  SigL = oo16pi2 * SigL
  SigR = oo16pi2 * SigR
  SigS = oo16pi2 * SigS

  Iname = Iname - 1

 End Subroutine Sigma_Neutralino


 Subroutine Sigma_Neutralino_2( mN, mN2, c_NNZ_L, c_NNZ_R, mS02, c_NNS0_L    &
            & , c_NNS0_R, mP02, c_NNP0_L, c_NNP0_R, mC, mC2, c_CNW_L, c_CNW_R &
            & , mSpm2, c_SmpCN_L, c_SmpCN_R, mSup2, c_UNSu_L, c_UNSu_R        &
            & , mSdown2, c_DNSd_L, c_DNSd_R, WriteOut, SigS)
 !-----------------------------------------------------------------------
 ! the following data are taken from the common /SM/
 ! the following data are taken from the common /numerical/
 !     - sqrt2
 ! the following data are taken from the common /EWparameters/
 !     - vevs, mu, M2
 ! the following data are taken from the common /Masses/
 ! the following data are taken from the common /Sugra/
 ! the following data are taken from the common /Sugra1/
 ! written by Werner Porod, 27.12.01
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mN(:), mN2(:), mC(:), mC2(:), mSup2(:)        &
      &  , mSdown2(:), mS02(:), mP02(:), mSpm2(:)
  Complex(dp), Intent(in) :: c_NNZ_L(:,:), c_NNZ_R(:,:)                     &
     & , c_NNS0_L(:,:,:), c_NNS0_R(:,:,:), c_NNP0_L(:,:,:), c_NNP0_R(:,:,:) &
     & , c_UNSu_L(:,:,:), c_UNSu_R(:,:,:), c_CNW_L(:,:), c_CNW_R(:,:)     &
     & , c_SmpCN_L(:,:,:), c_SmpCN_R(:,:,:), c_DNSd_L(:,:,:), c_DNSd_R(:,:,:)
  Logical, Intent(in) :: WriteOut
  Complex(dp), Intent(out) :: SigS(:,:)

  Integer :: n_neut, n_char, n_S0, n_P0, n_Spm
  Integer :: i1, i2, i3, i4, i2_min, i2_max
  Real(dp), Allocatable :: B1m2_i(:), B0m2_i(:)
  Complex(dp), Allocatable :: SumL(:,:), SumR(:,:), SumS(:,:)

  !------------------
  ! Inititalisation
  !------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "Sigma_Neutralino"
  SigS = ZeroC
  n_neut = Size( mN )
  n_char = Size( mC )
  n_S0 = Size( mS02 )
  n_P0 = Size( mP02 )
  n_Spm = Size( mSpm2 )

  Allocate(SumL(n_neut,n_neut))
  Allocate(SumR(n_neut,n_neut))
  Allocate(SumS(n_neut,n_neut))
  Allocate(B0m2_i(n_neut))
  Allocate(B1m2_i(n_neut))
  !---------------------
  ! quark - squark
  !---------------------
  Do i1=1,3
   If (GenerationMixing) Then
    i2_min = 1
    i2_max = 6
   Else
    i2_min = 2*i1-1
    i2_max = 2*i1
   End If
   Do i2=i2_min,i2_max
    If (LoopContributions(1)) Then
     Do i3=1,n_neut
      B1m2_i(i3) = - 3._dp * mN(i3) * B1(mN2(i3), mf_u2(i1), mSup2(i2))
      B0m2_i(i3) = 6._dp * mf_u(i1) * B0(mN2(i3), mf_u2(i1), mSup2(i2))
     End Do
    Else
     B1m2_i = 0._dp
     B0m2_i = 0._dp
    End If
    Do i3=1,n_neut
     SumL(i3,i3) = ( c_UNSu_L(i1,i3,i2) * c_UNSu_L(i1,i3,i2)     &
               &   + c_UNSu_R(i1,i3,i2) * c_UNSu_R(i1,i3,i2)  )  &
               & * B1m2_i(i3)
     SumS(i3,i3) = c_UNSu_R(i1,i3,i2) * c_UNSu_L(i1,i3,i2) * B0m2_i(i3)
     Do i4=i3+1,n_neut
      SumL(i3,i4) = ( c_UNSu_L(i1,i4,i2) * c_UNSu_L(i1,i3,i2)       &
                &   + c_UNSu_R(i1,i4,i2) * c_UNSu_R(i1,i3,i2)  )  &
                & * 0.5_dp * (B1m2_i(i3) + B1m2_i(i4) )
      SumS(i3,i4) = ( c_UNSu_R(i1,i3,i2) * c_UNSu_L(i1,i4,i2)    &
                &   + c_UNSu_R(i1,i4,i2) * c_UNSu_L(i1,i3,i2) )  &
                & * 0.25_dp * ( B0m2_i(i3) + B0m2_i(i4) )
      SumL(i4,i3) = SumL(i3,i4)
      SumS(i4,i3) = SumS(i3,i4)
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "u-quark, u-squark",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
     End Do 
    End If
    SigS = SigS - SumS - SumL

    If (LoopContributions(2)) Then
     Do i3=1,n_neut
      B1m2_i(i3) = - 3._dp * mN(i3) * B1(mN2(i3), mf_d2(i1), mSdown2(i2))
      B0m2_i(i3) = 6._dp * mf_d(i1) * B0(mN2(i3), mf_d2(i1), mSdown2(i2))
     End Do
    Else
     B1m2_i = 0._dp
     B0m2_i = 0._dp
    End If
    Do i3=1,n_neut
     SumL(i3,i3) = ( c_DNSd_L(i1,i3,i2) * c_DNSd_L(i1,i3,i2)     &
               &   + c_DNSd_R(i1,i3,i2) * c_DNSd_R(i1,i3,i2)  )  &
               & * B1m2_i(i3)
     SumS(i3,i3) = c_DNSd_R(i1,i3,i2) * c_DNSd_L(i1,i3,i2) * B0m2_i(i3)
!     If (i3.Eq.n_neut) Write(*,*) i1,i2,oo16pi2*sums(i3,i3),b0m2_i(i3),b1m2_i(i3)
     Do i4=i3+1,n_neut
      SumL(i3,i4) = ( c_DNSd_L(i1,i4,i2) * c_DNSd_L(i1,i3,i2)       &
                &   + c_DNSd_R(i1,i4,i2) * c_DNSd_R(i1,i3,i2)  )  &
                & * 0.5_dp * (B1m2_i(i3) + B1m2_i(i4) )
      SumS(i3,i4) = ( c_DNSd_R(i1,i3,i2) * c_DNSd_L(i1,i4,i2)    &
                &   + c_DNSd_R(i1,i4,i2) * c_DNSd_L(i1,i3,i2) )  &
                & * 0.25_dp * ( B0m2_i(i3) + B0m2_i(i4) )
      SumL(i4,i3) = SumL(i3,i4)
      SumS(i4,i3) = SumS(i3,i4)
     End Do
    End Do
    If (WriteOut) Then
     Write(ErrCan,*) "d-quark, d-squark",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
     End Do 
    End If
    SigS = SigS - SumS - SumL
   End Do
  End Do


 !-------------------------------------
 ! chargino - W
 !-------------------------------------
  Do i1 = 1,n_char
   If (LoopContributions(3)) Then
     Do i3=1,n_neut
      B1m2_i(i3) =  - 2._dp * mN(i3) *B1(mN2(i3),mC2(i1),mW2)
      B0m2_i(i3) =  - 8._dp * mC(i1) * B0(mN2(i3),mC2(i1),mW2)
     End Do
    Else
     B1m2_i = 0._dp
     B0m2_i = 0._dp
    End If

   Do i2 = 1, n_neut
     SumL(i2,i2) = ( c_CNW_R(i1,i2) * c_CNW_R(i1,i2)     &
               &   + c_CNW_L(i1,i2) * c_CNW_L(i1,i2)  )  &
               & * B1m2_i(i2)
     SumS(i2,i2) = c_CNW_L(i1,i2) * c_CNW_R(i1,i2) * B0m2_i(i2)
    Do i3 = i2+1,n_neut
      SumL(i2,i3) = ( c_CNW_R(i1,i3) * c_CNW_R(i1,i2)       &
                &   + c_CNW_L(i1,i3) * c_CNW_L(i1,i2)  )  &
                & * 0.5_dp * (B1m2_i(i2) + B1m2_i(i3) )
      SumS(i2,i3) = ( c_CNW_L(i1,i2) * c_CNW_R(i1,i3)    &
                &   + c_CNW_L(i1,i3) * c_CNW_R(i1,i2) )  &
                & * 0.25_dp * ( B0m2_i(i2) + B0m2_i(i3) )
      SumL(i3,i2) = SumL(i2,i3)
      SumS(i3,i2) = SumS(i2,i3)
    End Do 
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "chargino, W",i1
    Do i3=1,n_neut
     Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
    End Do 
    Do i3=1,n_neut
     Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
    End Do 
   End If
   SigS = SigS - SumS - sumL
  End Do
 !-------------------------------------
 ! chargino - S+
 !-------------------------------------
  Do i1 = 1,n_char
   Do i2 = 1,n_Spm
    If (LoopContributions(4)) Then
     Do i3=1,n_neut
      B1m2_i(i3) =  - mN(i3) * B1(mN2(i3),mC2(i1),mSpm2(i2))
      B0m2_i(i3) =  2._dp * mC(i1) * B0(mN2(i3),mC2(i1),mSpm2(i2))
     End Do
    Else
     B1m2_i = 0._dp
     B0m2_i = 0._dp
    End If
    Do i3 = 1,n_neut
     SumL(i3,i3) = ( c_SmpCN_L(i2,i1,i3) * c_SmpCN_L(i2,i1,i3)     &
               &   + c_SmpCN_R(i2,i1,i3) * c_SmpCN_R(i2,i1,i3)  )  &
               & * B1m2_i(i3)
     SumS(i3,i3) = c_SmpCN_R(i2,i1,i3) * c_SmpCN_L(i2,i1,i3) * B0m2_i(i3)
     Do i4 = i3+1,n_neut
      SumL(i3,i4) = ( c_SmpCN_L(i2,i1,i4) * c_SmpCN_L(i2,i1,i3)       &
                &   + c_SmpCN_R(i2,i1,i4) * c_SmpCN_R(i2,i1,i3)  )  &
                & * 0.5_dp * (B1m2_i(i3) + B1m2_i(i4) )
      SumS(i3,i4) = ( c_SmpCN_R(i2,i1,i3) * c_SmpCN_L(i2,i1,i4)    &
                &   + c_SmpCN_R(i2,i1,i4) * c_SmpCN_L(i2,i1,i3) )  &
                & * 0.25_dp * ( B0m2_i(i3) + B0m2_i(i4) )
      SumL(i4,i3) = SumL(i3,i4)
      SumS(i4,i3) = SumS(i3,i4)
     End Do 
    End Do 
    If (WriteOut) Then
     Write(ErrCan,*) "Chargino, S+",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
     End Do 
    End If
    SigS = SigS - SumS - sumL
   End Do
  End Do
 !-------------------------------------
 ! neutral gauge bosons
 !-------------------------------------
  Do i1 = 1,n_neut
   If (LoopContributions(5)) Then
     Do i3=1,n_neut
      B1m2_i(i3) =  - mN(i3) * B1(mN2(i3),mN2(i1),mZ2) 
      B0m2_i(i3) =  - 4._dp * mN(i1) * B0(mN2(i3),mN2(i1),mZ2)
     End Do
   Else
    B1m2_i = 0._dp
    B0m2_i = 0._dp
   End If
   Do i2 = 1,n_neut
     SumL(i2,i2) = ( c_NNZ_R(i1,i2) * c_NNZ_R(i1,i2)     &
               &   + c_NNZ_L(i1,i2) * c_NNZ_L(i1,i2)  )  &
               & * B1m2_i(i2)
     SumS(i2,i2) = c_NNZ_L(i1,i2) * c_NNZ_R(i1,i2) * B0m2_i(i2)
    Do i3 = i2+1,n_neut
      SumL(i2,i3) = ( c_NNZ_R(i1,i3) * c_NNZ_R(i1,i2)       &
                &   + c_NNZ_L(i1,i3) * c_NNZ_L(i1,i2)  )  &
                & * 0.5_dp * (B1m2_i(i2) + B1m2_i(i3) )
      SumS(i2,i3) = ( c_NNZ_L(i1,i2) * c_NNZ_R(i1,i3)    &
                &   + c_NNZ_L(i1,i3) * c_NNZ_R(i1,i2) )  &
                & * 0.25_dp * ( B0m2_i(i2) + B0m2_i(i3) )
      SumL(i3,i2) = SumL(i2,i3)
      SumS(i3,i2) = SumS(i2,i3)
    End Do 
   End Do
   If (WriteOut) Then
    Write(ErrCan,*) "neutralino, Z",i1
    Do i3=1,n_neut
     Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
    End Do 
    Do i3=1,n_neut
     Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
    End Do 
   End If
   SigS = SigS - SumS - sumL
  End Do

 !-------------------------------------
 ! neutralino - P0 
 !-------------------------------------
  Do i1 = 1,n_neut
   Do i2 = 1,n_P0
    If (LoopContributions(6)) Then
     Do i3=1,n_neut
      B1m2_i(i3) =  - 0.5_dp * mN(i3) * B1(mN2(i3),mN2(i1),mP02(i2))
      B0m2_i(i3) =  mN(i1) * B0(mN2(i3),mN2(i1),mP02(i2))
     End Do
    Else
     B1m2_i = 0._dp
     B0m2_i = 0._dp
    End If
    Do i3 = 1,n_neut
     SumL(i3,i3) = ( c_NNP0_L(i3,i1,i2) * c_NNP0_L(i3,i1,i2)     &
               &   + c_NNP0_R(i3,i1,i2) * c_NNP0_R(i3,i1,i2)  )  &
               & * B1m2_i(i3)
     SumS(i3,i3) = c_NNP0_R(i3,i1,i2) * c_NNP0_L(i3,i1,i2) * B0m2_i(i3)
     Do i4 = i3+1,n_neut
      SumL(i3,i4) = ( c_NNP0_L(i4,i1,i2) * c_NNP0_L(i3,i1,i2)       &
                &   + c_NNP0_R(i4,i1,i2) * c_NNP0_R(i3,i1,i2)  )  &
                & * 0.5_dp * (B1m2_i(i3) + B1m2_i(i4) )
      SumS(i3,i4) = ( c_NNP0_R(i3,i1,i2) * c_NNP0_L(i4,i1,i2)    &
                &   + c_NNP0_R(i4,i1,i2) * c_NNP0_L(i3,i1,i2) )  &
                & * 0.25_dp * ( B0m2_i(i3) + B0m2_i(i4) )
      SumL(i4,i3) = SumL(i3,i4)
      SumS(i4,i3) = SumS(i3,i4)
     End Do 
    End Do 
    If (WriteOut) Then
     Write(ErrCan,*) "Neutralino, P0",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
     End Do 
    End If
     ! the + is here because my couplings are multiplied by i
    SigS = SigS + 0.5_dp * SumS + 0.5_dp * SumL
   End Do
  End Do

 !-------------------------------------
 ! neutralino - S0
 !-------------------------------------
  Do i1 = 1,n_neut
   Do i2 = 1,n_S0
    If (LoopContributions(7)) Then
     Do i3=1,n_neut
      B1m2_i(i3) =  - 0.5_dp * mN(i3) * B1(mN2(i3),mN2(i1),mS02(i2))
      B0m2_i(i3) =  mN(i1) * B0(mN2(i3),mN2(i1),mS02(i2))
     End Do
    Else
     B1m2_i = 0._dp
     B0m2_i = 0._dp
    End If
    Do i3 = 1,n_neut
     SumL(i3,i3) = ( c_NNS0_L(i3,i1,i2) * c_NNS0_L(i3,i1,i2)     &
               &   + c_NNS0_R(i3,i1,i2) * c_NNS0_R(i3,i1,i2)  )  &
               & * B1m2_i(i3)
     SumS(i3,i3) = c_NNS0_R(i3,i1,i2) * c_NNS0_L(i3,i1,i2) * B0m2_i(i3)
     Do i4 = i3+1,n_neut
      SumL(i3,i4) = ( c_NNS0_L(i4,i1,i2) * c_NNS0_L(i3,i1,i2)       &
                &   + c_NNS0_R(i4,i1,i2) * c_NNS0_R(i3,i1,i2)  )  &
                & * 0.5_dp * (B1m2_i(i3) + B1m2_i(i4) )
      SumS(i3,i4) = ( c_NNS0_R(i3,i1,i2) * c_NNS0_L(i4,i1,i2)    &
                &   + c_NNS0_R(i4,i1,i2) * c_NNS0_L(i3,i1,i2) )  &
                & * 0.25_dp * ( B0m2_i(i3) + B0m2_i(i4) )
      SumL(i4,i3) = SumL(i3,i4)
      SumS(i4,i3) = SumS(i3,i4)
     End Do 
    End Do 
    If (WriteOut) Then
     Write(ErrCan,*) "Neutralino, S0",i1,i2
     Do i3=1,n_neut
      Write(ErrCan,*) "SigL(i,:)",i3,SumL(i3,:)
     End Do 
     Do i3=1,n_neut
      Write(ErrCan,*) "SigS(i,:)",i3,SumS(i3,:)
     End Do 
    End If
    SigS = SigS -  0.5_dp * SumS - 0.5_dp * SumL
   End Do
  End Do

  SigS = oo16pi2 * SigS

  Deallocate(SumL, SumR, SumS, B0m2_i, B1m2_i)

  Iname = Iname - 1

 End Subroutine Sigma_Neutralino_2


#ifdef GENERATIONMIXING
 Subroutine SleptonMass_1L(gU1, gSU2, Y_l, Y_d, vevSM, M_L2, M_R2, Af     &
        & , mu_T, mu, mN, mN2, N, mC2, U, V, mSup2, RSup, mSdown2, RSdown &
        & , mSlepton2, RSlepton, mSneut2, RSneut, mS02, RS0, mP02, RP0    &
        & , mSpm2, RSpm, uL_L, uL_R, mZ2, mW2, delta, mass, mass2, Rsf, kont)
#else
 Subroutine SleptonMass_1L(gU1, gSU2, Y_l, Y_d, vevSM, M_L2, M_R2, Af       &
        & , mu_T, mu, mN, mN2, N, mC2, U, V, mSup2, RSup, mSdown2, RSdown   &
        & , mSlepton2, RSlepton, mSneut2, mS02, RS0, mP02, RP0, mSpm2, RSpm &
        & , mZ2, mW2, delta, mass, mass2, Rsf, kont)
#endif
 !-------------------------------------------------------------------------
 ! calculates the 1-loop corrected slepton masses + mixing angles for
 ! three generations but without generation mixing.
 ! input:
 ! output:
 !  - mass(i) ...... masses
 !  - mass2(i) ..... masses squared
 !  - Rsf(i,j) ..... mixing matrix
 ! written by Werner Porod, 28.12.01
 !  28.12.01: write first tree-level mass matrix as in case of generation
 !            mixing; check during the calculation of the couplings and
 !            diagonalisation whether generation mixing is included or not
 !-------------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: gU1, gSU2, vevSM(2), mN(4), mN2(4), mC2(2) &
      &  , mSup2(6), mSdown2(6), mSlepton2(6), mSneut2(3), mS02(2), RS0(2,2)  &
      &  , mP02(2), RP0(2,2), mSpm2(2), mZ2, mW2, delta
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), M_L2(3,3), mu, N(4,4)     &
     & , U(2,2), V(2,2), RSup(6,6), M_R2(3,3), Af(3,3), RSdown(6,6)        &
     & , RSlepton(6,6), RSpm(2,2), mu_T
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: RSneut(3,3), uL_L(3,3), uL_R(3,3)
#endif
  Complex(dp), Intent(out) :: Rsf(6,6)
  Real(dp), Intent(out) :: mass(6), mass2(6)

  Integer :: i1, i2, i3, i4, i5, i_min, i_max, i_gen1, i_gen2, i_count
  Real(dp) :: vev2, p2, Yl, Yr, sinW2, ed, eu, e, T3, c_Sq2Z(6), c_Sq2W(6) &
     & , mi2(2), test_m2(6)
#ifdef GENERATIONMIXING
  Real(dp) :: mi6(6), test(2)
  Complex(dp) :: mat6(6,6), Rsfp3(3,3)
#endif
  Complex(dp) :: off(3,3), YukC(3,3), YukT(3,3), mat6a(6,6)           &
     & , mat2(2,2), PiSf(6,6,6), Rsf2(2,2), Rsf6(6,6)                  &
     & , c_Sq4e(6,6,6), c_Sq4Z(6,6,6), c_Sq4W(3,6,6), coupC1, coupC2          &
     & , c_FNSf_L(3,4,6), c_FNSf_R(3,4,6), c_CFpSf_L(2,3,6), c_CFpSf_R(2,3,6) &
     & , c_P0P0SfSf(2,6,6), c_P0SfSf(2,6,6), bi(1), c_S0S0SfSf(2,6,6)         &
     & , c_S0SfSf(2,6,6), c_SpmSpmSfSfp(2,6,6), c_SpmSfSfp(2,6,3), mat3(3,3)  &
     & , c_SuSf4(6,6,6), c_SdSf4(6,6,6), c_SlSf4(6,6,6), c_SnSf4(3,6,6)       &
     & , Rsfp2(2,2), coupC
  Logical :: WriteOut

  !------------------------------------
  ! initialisation
  !------------------------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'SleptonMass_1L'

  If ((WriteOneLoopContributions.Eq.13).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions to SleptonMass_1L:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If

  PiSf = ZeroC
  mass = 0._dp
  mass2 = 0._dp
  RSf = ZeroC

  sinW2 = gU1**2 / (gSU2**2 + gU1**2)
  T3 = - 0.5_dp
  e = -1._dp
  mat3 = ZeroC ! needed instead of A_nu and Y_nu
  Rsf6 = RSlepton
#ifdef GENERATIONMIXING
  Rsfp3 = Rsneut
#endif
  !-------------------------
  ! tree level mass matrix
  !-------------------------
  Yl = -1._dp
  Yr = 2._dp
  vev2 = 0.25_dp * (vevSM(1)**2 - vevSM(2)**2)

  YukT = Transpose(Y_l)
  YukC = Conjg(Y_l)

  mat6a(1:3,1:3) = M_L2 + 0.5_dp * vevSM(1)**2 * Matmul(YukC,YukT) &
      &         + (T3 * gSU2**2 - 0.5_dp * Yl * gU1**2) * vev2 * id3C
  mat6a(4:6,4:6) = M_R2 + 0.5_dp * vevSM(1)**2 * Matmul(YukT,YukC) &
      &         - 0.5_dp * Yr * gU1**2 * vev2 * id3C
  off = (vevSM(1) * Af - Conjg(mu) * vevSM(2) * Y_l ) * oosqrt2

  off = Conjg(off)
  mat6a(1:3,4:6) = off
  Call Adjungate(off)
  mat6a(4:6,1:3) = off

  !-----------
  ! couplings
  !-----------
   c_Sq4e = ZeroC
   c_Sq4Z = ZeroC
   c_Sq4W = ZeroC
   c_Sq2Z = ZeroC
   c_Sq2W = ZeroC
#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Do i1=1,6
     Call CoupSfermionZ(i1, i1, gSU2, sinW2, e, T3, id6c, coupC1)
     c_Sq2Z(i1) = 4._dp * Abs( coupC1 )**2
     Do i2=1,6
      Call CoupSfermionZa(i2, i1, gSU2, sinW2, e, T3, RSf6, coupC1)
      Do i3=1,6
       c_Sq4e(i1,i2,i3) = Rsf6(i1,i2) * Conjg( Rsf6(i1,i3) )
       Call CoupSfermionZa(i3, i1, gSU2, sinW2, e, T3, RSf6, coupC2)
       c_Sq4Z(i1,i2,i3) = coupC1 * Conjg(coupC2)
      End Do
     End Do
    End Do  
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       c_Sq4W(i1,i2,i3) = Rsfp3(i1,i2) * Conjg( Rsfp3(i1,i3) ) 
      End Do
     End Do
    End Do  
    c_Sq2W(1) = 2._dp * gSU2**2
    c_Sq2W(2) = c_Sq2W(1)
    c_Sq2W(3) = c_Sq2W(1)

   Else ! .not. GenerationMixing
#endif
    Rsfp2 = id2C
    Do i1=1,3
     Rsf2 = Rsf6(2*i1-1:2*i1,2*i1-1:2*i1)
     Do i2=1,2
      Do i3=1,2
       i_gen1 = 2*(i1-1) + i3
       Call CoupSfermionZa(i3, i2, gSU2, sinW2, e, T3, RSf2, coupC1)
       Do i4=1,2
        i_gen2 = 2*(i1-1) + i4
        c_Sq4e(2*(i1-1)+i2,i_gen1,i_gen2) = Rsf2(i2,i3) * Conjg( Rsf2(i2,i4) )
        Call CoupSfermionZa(i4, i2, gSU2, sinW2, e, T3, RSf2, coupC2)
        c_Sq4Z(2*(i1-1)+i2,i_gen1,i_gen2) = coupC1 * Conjg(coupC2)
       End Do
      End Do
     End Do
     i_gen1 = 2*i1-1
     c_Sq4W(i1,i_gen1,i_gen1) = 1._dp
    End Do
    Call CoupSfermionZ(1, 1, gSU2, sinW2, e, T3, id2c, coupC1)
    c_Sq2Z(1) = 4._dp * Abs( coupC1 )**2
    c_Sq2Z(3) = c_Sq2Z(1)
    c_Sq2Z(5) = c_Sq2Z(1)
    Call CoupSfermionZ(2, 2, gSU2, sinW2, e, T3, id2c, coupC1)
    c_Sq2Z(2) = 4._dp * Abs( coupC1 )**2
    c_Sq2Z(4) = c_Sq2Z(2)
    c_Sq2Z(6) = c_Sq2Z(2)    
    c_Sq2W(1) = 2._dp * gSU2**2
    c_Sq2W(3) = c_Sq2W(1)
    c_Sq2W(5) = c_Sq2W(1)
#ifdef GENERATIONMIXING
   End If ! GenerationMixing
#endif
   c_Sq4e = c_Sq4e * gSU2**2 * sinW2
   c_Sq4W = 0.5_dp * gSU2**2 * c_Sq4W
   

  c_FNSf_L = ZeroC
  c_FNSf_R = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,4
     Do i3=1,6
      Call CoupNeutralinoSfermion(i1, i2, i3, gU1, gSU2, T3, e, id6C, uL_L  &
                       &, uL_R, Y_l, N, c_FNSf_L(i1,i2,i3), c_FNSf_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else ! .not. GenerationMixing
#endif
   Do i1=1,3
    Do i2=1,4
     Do i3=1,2
      Call CoupNeutralinoSfermion(i2, i3, gU1, gSU2, T3, e, id2C, Y_l(i1,i1) &
               & , N, c_FNSf_L(i1,i2,2*(i1-1)+i3), c_FNSf_R(i1,i2,2*(i1-1)+i3))
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_CFpSf_L = ZeroC
  c_CFpSf_R = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,2
    Do i2=1,3
     Do i3=1,6
      Call CoupCharginoSfermion(i1, i2, i3, gSU2, -T3, id6C, Y_l, mat3, id3C  &
                      &, id3C, U, V, c_CFpSf_L(i1,i2,i3), c_CFpSf_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else ! .not. GenerationMixing
#endif
   Do i1=1,2
    Do i2=1,3
     Do i3=1,2
      Call CoupCharginoSfermion(i1,i3,gSU2, -T3, id2C, Y_l(i2,i2), ZeroC &
          & , U, V, c_CFpSf_L(i1,i2,2*(i2-1)+i3), c_CFpSf_R(i1,i2,2*(i2-1)+i3))
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  bi(1) = mu_T
  c_P0SfSf = ZeroC
  c_P0P0SfSf = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_P0
    Do i2=1,6
     Do i3=1,6
      Call CoupPseudoscalarSfermion4(i1, i1, i2, i3, RP0, T3, e, Y_l  &
                                    &, id6C, gU1, gSU2, c_P0P0SfSf(i1,i2,i3) )
      Call CoupPseudoScalarSfermion3a(i1, i2, i3, RP0, T3, Y_l, Rsf6, Af, bi &
                                    &, c_P0SfSf(i1,i2,i3) )
     End Do   
    End Do   
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,n_P0
    Do i2=1,3
     Rsf2 = Rsf6(2*i2-1:2*i2,2*i2-1:2*i2)      
     Do i3=1,2
      Do i4=1,2
       Call CoupPseudoscalarSfermion4(i1, i1, i3, i4, RP0, T3, e, Y_l(i2,i2)  &
              &, id2C, gU1, gSU2, c_P0P0SfSf(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
       Call CoupPseudoScalarSfermion3a(i1, i3, i4, RP0, T3, Y_l(i2,i2), Rsf2 &
                     &, Af(i2,i2), bi, c_P0SfSf(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
      End Do   
     End Do   
    End Do   
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_S0SfSf = ZeroC
  c_S0S0SfSf = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_S0
    Do i2=1,6
     Do i3=1,6
      Call CoupScalarSfermion4(i1, i1, i2, i3, RS0, T3, e, Y_l  &
                             &, id6C, gU1, gSU2, c_S0S0SfSf(i1,i2,i3) )
      Call CoupScalarSfermion3a(i1, i2, i3, RS0, T3, e, Y_l, Rsf6, Af , mu_T &
                              &, vevSM, gU1, gSU2, c_S0SfSf(i1,i2,i3) )
     End Do   
    End Do   
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,n_S0
    Do i2=1,3
     Rsf2 = Rsf6(2*i2-1:2*i2,2*i2-1:2*i2)      
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion4(i1, i1, i3, i4, RS0, T3, e, Y_l(i2,i2)  &
              &, id2C, gU1, gSU2, c_S0S0SfSf(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
       Call CoupScalarSfermion3a(i1, i3, i4, RS0, T3, e, y_l(i2,i2), Rsf2 &
                                &, Af(i2,i2), mu_T, vevSM, gU1, gSU2        &
                                &, c_S0SfSf(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
      End Do   
     End Do   
    End Do   
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_SpmSfSfp = ZeroC
  c_SpmSpmSfSfp = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_Spm
    Do i2=1,6
     Do i3=1,6
      Call CoupChargedScalarSfermion4(i1, i1, i2, i3, RSpm, T3, e, Y_l  &
                           &, mat3, id6C, gU1, gSU2, c_SpmSpmSfSfp(i1,i2,i3) )
     End Do   
    End Do   
    Do i2=1,6
     Do i3=1,3
      Call CoupChargedScalarSfermion3(i1, i2, i3, RSpm, gSU2, vevSM, mu_T  &
                 & , Y_l, Af, id6c, mat3, mat3, RSneut,  c_SpmSfSfp(i1,i2,i3) )
     End Do   
    End Do   
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,n_Spm
    Do i2=1,3
     Do i3=1,2
      Call CoupChargedScalarSfermion3(i1, i3, 1, RSpm, gSU2, vevSM, mu_T &
               & , Y_l(i2,i2), Af(i2,i2), id2c, ZeroC, ZeroC, id2C &
               & , c_SpmSfSfp(i1,2*(i2-1)+i3,i2) )

      Do i4=1,2
       Call CoupChargedScalarSfermion4(i1, i1, i3, i4, RSpm, T3, e  &
                                &, Y_l(i2,i2), ZeroC, id2C, gU1, gSU2   &
                                &, c_SpmSpmSfSfp(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
      End Do   
     End Do   
    End Do   
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_SuSf4 = ZeroC
  c_SdSf4 = ZeroC
  c_SlSf4 = ZeroC
  c_SnSf4 = ZeroC
  ed =  -1._dp / 3._dp
  eu =  2._dp / 3._dp
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,6
    i_gen1 = (i1+1)/2
    Do i2=1,6
     i_gen2 = (i2+1)/2
     Do i3=1,6
      Call CoupSfermion4Y(i2, i3, i1, i1, Y_l, RSlepton, 1, c_SlSf4(i1,i2,i3), 1)
      Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, e, T3, RSf6, 1, coupC, 1)
      c_SlSf4(i1,i2,i3) = c_SlSf4(i1,i2,i3) + coupC
      Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, T3, e, id6C, -T3, eu &
                         &, RSup, 3, .False., c_SuSf4(i1,i2,i3) )
      Call CoupSfermion4Y(i1, i1, i2, i3, T3, Y_d, RSdown, T3, Y_l, id6c, 3 &
                          &,c_SdSf4(i1,i2,i3) )
      Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, T3, e, id6C, T3, ed &
                         &, RSdown, 3, .False., coupC )
      ! be aware of colour
      c_SdSf4(i1,i2,i3) = 3._dp * c_SdSf4(i1,i2,i3) + coupC
     End Do
    End Do
   End Do
  Else ! .not. GenerationMixing
#endif

   Do i1=1,3
    i_min = 2*i1-1
    i_max = i_min + 1
    Rsf2 = Rsf6(i_min:i_max,i_min:i_max)
    Rsfp2 = RSup(i_min:i_max,i_min:i_max)
    Do i2=1,2
     i_min = 2*(i1-1)+i2
     Do i3=1,3
      Do i4=1,2
       i_gen1 = 2*(i3-1)+i4
       Do i5=1,2
        i_gen2 = 2*(i3-1)+i5
        If (i1.Eq.i3) Then
         Call CoupSfermion4Y(i4, i5, i2, i2, Y_l(i1,i1), Rsf2, 1 &
                           &, c_SlSf4(i_min,i_gen1,i_gen2), 1)
         Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, e, T3, RSf2,1,coupC,1)
         c_SlSf4(i_min,i_gen1,i_gen2) = c_SlSf4(i_min,i_gen1,i_gen2) + coupC
        Else
         Call CoupSfermion4Y(i2, i2, i4, i5, T3, Y_l(i1,i1), RSf2, T3 &
                        &, Y_l(i3,i3), id2c, 1, c_SlSf4(i_min,i_gen1,i_gen2) )
         Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, T3, e, id2C &
                         &, T3, e, RSf2, 1, .False., coupC)
         c_SlSf4(i_min,i_gen1,i_gen2) = c_SlSf4(i_min,i_gen1,i_gen2) + coupC
        End If
        Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, T3, e, id2C, -T3, eu &
                         &, RSfp2, 3, .False., c_SuSf4(i_min,i_gen1,i_gen2) )
       End Do
      End Do
     End Do
    End Do
    i_min = 2*i1-1
    i_max = i_min + 1
    Rsfp2 = RSdown(i_min:i_max,i_min:i_max)
    Do i2=1,2
     i_min = 2*(i1-1)+i2
     Do i3=1,3
       Do i4=1,2
        i_gen1 = 2*(i3-1)+i4
        Do i5=1,2
         i_gen2 = 2*(i3-1)+i5
         Call CoupSfermion4Y(i2, i2, i4, i5, T3, Y_d(i1,i1), RSfp2, T3 &
                           &, Y_l(i3,i3), id2c, 3, c_SdSf4(i_min,i_gen1,i_gen2) )
         Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, T3, e, id2C, T3    &
                           & , ed, RSfp2, 3, .False., coupC )
      ! be aware of colour
         c_SdSf4(i_min,i_gen1,i_gen2) = 3._dp*c_SdSf4(i_min,i_gen1,i_gen2) &
               & + coupC
        End Do
      End Do
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  eu = 0._dp
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,6
    i_gen1 = (i2+1)/2
     Do i3=1,6
      Call CoupSfermion4Y(i1, i1, i2, i3, -T3, mat3, RSneut, T3, Y_l, id6c, 1 &
                         &,c_SnSf4(i1,i2,i3) )
      Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, T3, e, id6C, -T3, eu &
                         &, RSneut, 1, .True., coupC )
      c_SnSf4(i1,i2,i3) = c_SnSf4(i1,i2,i3) + coupC
     End Do
    End Do
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,3
    Do i3=1,2
     i_gen1 = 2*(i1-1)+i3
     Do i4=1,2
      i_gen2 = 2*(i1-1)+i4
      Call CoupSfermion4Y(1, 1, i3, i4, -T3, ZeroC, id2C, T3 &
                        &, Y_l(i1,i1), id2c, 1, c_SnSf4(i1,i_gen1,i_gen2) )
     End Do
    End Do
    Do i2=1,3
     Do i3=1,2
      i_gen1 = 2*(i2-1)+i3
      Do i4=1,2
       i_gen2 = 2*(i2-1)+i4
       If (i2.Eq.i1) Then
        Call CoupSfermion4G(i3, i4, 1, 1, gU1, gSU2, T3, e, id2C, 0.5_dp, eu &
                          &, id2c, 1, .True., coupC)
       Else 
        Call CoupSfermion4G(i3, i4, 1, 1, gU1, gSU2, T3, e, id2C, 0.5_dp, eu &
                          &, id2c, 1, .False., coupC)
       End If
       c_SnSf4(i1,i_gen1,i_gen2) = c_SnSf4(i1,i_gen1,i_gen2) + coupC
      End Do
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif
  !-------------------
  ! loop contribution
  !-------------------
  Do i1=1,6
   p2 = mSlepton2(i1) ! Sum(mSlepton2) / 6._dp

   Call PiSlepton(p2, mSlepton2, c_Sq4e, c_Sq4Z, c_Sq2Z, c_Sq2W, mSneut2      &
        & , c_Sq4W, mN, mN2, c_FNSf_L, c_FNSf_R, mC2, c_CFpSf_L, c_CFpSf_R    &
        & , mP02, c_P0SfSf, c_P0P0SfSf, mS02, c_S0SfSf, c_S0S0SfSf            &
        & , mSpm2, c_SpmSfSfp, c_SpmSpmSfSfp, mSup2, c_SuSf4                  &
        & , mSdown2, c_SdSf4, c_SlSf4, c_SnSf4, mZ2, mW2, WriteOut            &
        & , PiSf(i1,:,:) )
  End Do
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=6,1,-1
    mat6 = mat6a - PiSf(i1,:,:)
    Call Chop(mat6) ! to avoid numerical problems 
    Call EigenSystem(mat6,mi6,Rsf,kont, test)

    If ((kont.Eq.-14).Or.(kont.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(ErrCan,*) "test =",test
     Write(ErrCan,*) " "
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = 0
    End If

    If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
     Write(ErrCan,*) 'Diagonalization did not work in routine SleptonMass_1L!'
     Write(ErrCan,*) 'msf2 ', mass2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'M_R2, ',M_R2
     Write(ErrCan,*) 'A_f, ',Af
     Write(ErrCan,*) 'Y_f, ',Y_l
     Write(ErrCan,*) 'mu, vevs, T3, Yl, Yr, g, gU1'
     Write(ErrCan,*) mu, vevSM, T3, Yl, Yr, gSU2, gU1
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    mass2(i1) = mi6(i1)
   End Do

  Else ! .not.GenerationMixing
#endif
   Do i1=1,3
    Do i2=2,1,-1
     mat2(1,1) = mat6a(i1,i1) - PiSf(2*(i1-1)+i2,2*i1-1,2*i1-1)
     mat2(2,2) = mat6a(3+i1,3+i1) - PiSf(2*(i1-1)+i2,2*i1,2*i1)
     mat2(1,2) = mat6a(i1,3+i1) - PiSf(2*(i1-1)+i2,2*i1-1,2*i1)
     mat2(2,1) = mat6a(3+i1,i1) - PiSf(2*(i1-1)+i2,2*i1,2*i1-1)
     Call Diag2(mat2,mi2,Rsf2)
     If (WriteOut) Write(ErrCan,*) "tree(11,12,22)",i1,mat6a(i1,i1) &
   &        , mat6a(i1,3+i1), mat6a(3+i1,3+i1)
     If (WriteOut) Write(ErrCan,*) "pi(11,12,22)",i1    &
        & , PiSf(2*(i1-1)+i2,2*i1-1,2*i1-1) &
        & , PiSf(2*(i1-1)+i2,2*i1-1,2*i1), PiSf(2*(i1-1)+i2,2*i1,2*i1)
     If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
      Write(ErrCan,*) 'msf2 ', mass2
      Write(ErrCan,*) 'M_L2, ',M_L2
      Write(ErrCan,*) 'M_R2, ',M_R2
      Write(ErrCan,*) 'A_f, ',Af
      Write(ErrCan,*) 'Y_f, ',Y_l
      Write(ErrCan,*) 'mu, vevs, g, gU1'
      Write(ErrCan,*) mu, vevSM, gSU2, gU1
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
     mass2(2*(i1-1)+i2) = mi2(i2)
     Rsf(2*i1-1:2*i1,2*i1-1:2*i1) = Rsf2
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If ! GenerationMixing 
#endif

  Do i1=1,6
   If (mass2(i1).Gt.0._dp) Then
    mass(i1) = Sqrt( mass2(i1) )
   Else 
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from routine SleptonMass_1L!'
     Write(ErrCan,*) 'in the calculation of the masses'
     Write(ErrCan,*) 'occurred a negative mass squared!!!'
     Write(ErrCan,*) 'Diagonalization did not work in routine SleptonMass_1L'
     Write(ErrCan,*) "i1 = ",i1
     Do i2=1,6
      Write(ErrCan,*) "mat6a",mat6a(i2,:)
     End Do
     Write(ErrCan,*) 'msf2 ', mass2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'M_R2, ',M_R2
     Write(ErrCan,*) 'A_f, ',Af
     Write(ErrCan,*) 'Y_f, ',Y_l
     Write(ErrCan,*) 'mu, vevs, g, gU1'
     Write(ErrCan,*) mu, vevSM, gSU2, gU1

     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    kont = -501
    Call AddError(501)
    mass(i1) = 0._dp
   End If
  End Do

  !-------------------------------------
  ! redoing calculation using refined p2
  !-------------------------------------
  i_count = 0
  Do 
   i_count = i_count + 1
   test_m2 = mass2
  Do i1=1,6
   p2 = mass2(i1) ! Sum(mSlepton2) / 6._dp

   Call PiSlepton(p2, mSlepton2, c_Sq4e, c_Sq4Z, c_Sq2Z, c_Sq2W, mSneut2      &
        & , c_Sq4W, mN, mN2, c_FNSf_L, c_FNSf_R, mC2, c_CFpSf_L, c_CFpSf_R    &
        & , mP02, c_P0SfSf, c_P0P0SfSf, mS02, c_S0SfSf, c_S0S0SfSf            &
        & , mSpm2, c_SpmSfSfp, c_SpmSpmSfSfp, mSup2, c_SuSf4                  &
        & , mSdown2, c_SdSf4, c_SlSf4, c_SnSf4, mZ2, mW2, WriteOut            &
        & , PiSf(i1,:,:) )
  End Do
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=6,1,-1
    mat6 = mat6a - PiSf(i1,:,:)
    Call Chop(mat6) ! to avoid numerical problems 
    Call EigenSystem(mat6,mi6,Rsf,kont, test)

    If ((kont.Eq.-14).Or.(kont.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(ErrCan,*) "test =",test
     Write(ErrCan,*) " "
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = 0
    End If

    If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
     Write(ErrCan,*) 'Diagonalization did not work in routine SleptonMass_1L!'
     Write(ErrCan,*) 'msf2 ', mass2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'M_R2, ',M_R2
     Write(ErrCan,*) 'A_f, ',Af
     Write(ErrCan,*) 'Y_f, ',Y_l
     Write(ErrCan,*) 'mu, vevs, T3, Yl, Yr, g, gU1'
     Write(ErrCan,*) mu, vevSM, T3, Yl, Yr, gSU2, gU1
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    mass2(i1) = mi6(i1)
   End Do

  Else ! .not.GenerationMixing
#endif
   Do i1=1,3
    Do i2=2,1,-1
     mat2(1,1) = mat6a(i1,i1) - PiSf(2*(i1-1)+i2,2*i1-1,2*i1-1)
     mat2(2,2) = mat6a(3+i1,3+i1) - PiSf(2*(i1-1)+i2,2*i1,2*i1)
     mat2(1,2) = mat6a(i1,3+i1) - PiSf(2*(i1-1)+i2,2*i1-1,2*i1)
     mat2(2,1) = mat6a(3+i1,i1) - PiSf(2*(i1-1)+i2,2*i1,2*i1-1)
     Call Diag2(mat2,mi2,Rsf2)
     If (WriteOut) Write(ErrCan,*) "tree(11,12,22)",i1,mat6a(i1,i1) &
   &        , mat6a(i1,3+i1), mat6a(3+i1,3+i1)
     If (WriteOut) Write(ErrCan,*) "pi(11,12,22)",i1    &
        & , PiSf(2*(i1-1)+i2,2*i1-1,2*i1-1) &
        & , PiSf(2*(i1-1)+i2,2*i1-1,2*i1), PiSf(2*(i1-1)+i2,2*i1,2*i1)
     If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
      Write(ErrCan,*) 'msf2 ', mass2
      Write(ErrCan,*) 'M_L2, ',M_L2
      Write(ErrCan,*) 'M_R2, ',M_R2
      Write(ErrCan,*) 'A_f, ',Af
      Write(ErrCan,*) 'Y_f, ',Y_l
      Write(ErrCan,*) 'mu, vevs, g, gU1'
      Write(ErrCan,*) mu, vevSM, gSU2, gU1
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
     mass2(2*(i1-1)+i2) = mi2(i2)
     Rsf(2*i1-1:2*i1,2*i1-1:2*i1) = Rsf2
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If ! GenerationMixing 
#endif

  Do i1=1,6
    If (test_m2(i1).Ne.0._dp) Then
     test_m2(i1) = Abs(test_m2(i1) - mass2(i1)) / test_m2(i1)
    Else
     test_m2(i1) = Abs(mass2(i1))
    End If
   If (mass2(i1).Gt.0._dp) Then
    mass(i1) = Sqrt( mass2(i1) )
   Else 
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from routine SleptonMass_1L!'
     Write(ErrCan,*) 'in the calculation of the masses'
     Write(ErrCan,*) 'occurred a negative mass squared!!!'
     Write(ErrCan,*) 'Diagonalization did not work in routine SleptonMass_1L'
     Write(ErrCan,*) "i1 = ",i1
     Do i2=1,6
      Write(ErrCan,*) "mat6a",mat6a(i2,:)
     End Do
     Write(ErrCan,*) 'msf2 ', mass2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'M_R2, ',M_R2
     Write(ErrCan,*) 'A_f, ',Af
     Write(ErrCan,*) 'Y_f, ',Y_l
     Write(ErrCan,*) 'mu, vevs, g, gU1'
     Write(ErrCan,*) mu, vevSM, gSU2, gU1

     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    kont = -501
    Call AddError(501)
    mass(i1) = 0._dp
   End If
  End Do
   If (Maxval(test_m2).Lt. 0.1_dp*delta) Exit
   If (i_count.Gt.30) Then
    Write(*,*) "Problem in SleptonMassLoop",test_m2,mass2
    kont = - 502
     Call AddError(502)
    Exit
   End If
  End Do ! i_count


  Iname = Iname - 1

 End Subroutine SleptonMass_1L


#ifdef GENERATIONMIXING
 Subroutine SneutrinoMass_1L(gU1, gSU2, Y_l, vevSM, M_L2, A_l, mu, mN2, N     &
        & , mC, mC2, U, V, mSup2, RSup, mSdown2, RSdown, mSlepton2, RSlepton  &
        & , mSneut2, RSneut, mS02, RS0, mP02, RP0, mSpm2, RSpm, uL_L, uL_R    &
        & , mZ2, mW2, delta, mass, mass2, Rsf, kont)
#else
 Subroutine SneutrinoMass_1L(gU1, gSU2, Y_l, vevSM, M_L2, A_l, mu, mN2, N     &
        & , mC, mC2, U, V, mSup2, RSup, mSdown2, RSdown, mSlepton2, RSlepton  &
        & , mSneut2, mS02, RS0, mP02, RP0, mSpm2, RSpm, mZ2, mW2, delta       &
        & , mass, mass2, Rsf, kont)
#endif
 !-------------------------------------------------------------------------
 ! calculates the 1-loop corrected sneutrino masses + mixing angles for
 ! three generations but without generation mixing.
 ! input:
 ! output:
 !  - mass(i) ...... masses
 !  - mass2(i) ..... masses squared
 !  - Rsf(i,j) ..... mixing matrix
 ! written by Werner Porod, 28.12.01
 !  28.12.01: write first tree-level mass matrix as in case of generation
 !            mixing; check during the calculation of the couplings and
 !            diagonalisation whether generation mixing is included or not
 !-------------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: gU1, vevSM(2), mN2(4), mC(2), mC2(2) &
      &  , mSup2(6), mSdown2(6), mSlepton2(6), mSneut2(3), mS02(2), RS0(2,2)  &
      &  , mP02(2), RP0(2,2), mSpm2(2), gSU2, mZ2, mW2, delta
  Complex(dp), Intent(in) :: Y_l(3,3), M_L2(3,3), A_l(3,3), mu, RSpm(2,2)  &
     & , N(4,4), U(2,2), V(2,2), RSup(6,6), RSdown(6,6), RSlepton(6,6)
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: uL_L(3,3), uL_R(3,3), RSneut(3,3)
#endif
  Complex(dp), Intent(out) :: Rsf(3,3)
  Real(dp), Intent(out) :: mass(3), mass2(3)

  Integer :: i1, i2, i3, i4, i_min, i_max, i_count
  Real(dp) :: vev2, p2, sinW2, e2, T3, eu, ed, e, c_Sq2Z(3), c_Sq2W(3) &
     & , mi2(3), test_m2(3)
#ifdef GENERATIONMIXING
  Real(dp) :: test(2)
  Complex(dp) :: Rsf3(3,3), coupC2
#endif
  Complex(dp) :: mat3(3,3), mat3a(3,3), PiSf(3,3,3), Rsf2(2,2)                &
     & , c_Sq4Z(3,3,3), c_Sq4W(6,3,3), Rsfp6(6,6), c_FNSf_R(3,4,3)            &
     & , c_CFpSf_L(2,3,3), c_CFpSf_R(2,3,3), c_P0P0SfSf(2,3,3)                &
     & , c_S0S0SfSf(2,3,3), c_S0SfSf(2,3,3), c_SpmSpmSfSfp(2,3,3)             &
     & , c_SpmSfSfp(2,3,6), c_SuSf4(6,3,3), c_SdSf4(6,3,3), c_SlSf4(6,3,3)    &
     & , c_SnSf4(3,3,3), Rsfp2(2,2), coupC, coupC1
  Logical :: WriteOut

  !------------------------------------
  ! initialisation
  !------------------------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'SneutrinoMass_1L'

  If ((WriteOneLoopContributions.Eq.14).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions to SneutrinoMass_1L:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If

  PiSf = ZeroC
  mass = 0._dp
  mass2 = 0._dp
  RSf = ZeroC

  sinW2 = gU1**2 / (gSU2**2 + gU1**2)
  T3 = 0.5_dp
  e = 0._dp
  Rsfp6 = RSlepton
  mat3 = ZeroC ! needed instead of Y_nu and A_nu

  !-------------------------
  ! tree level mass matrix
  !-------------------------
  vev2 = 0.25_dp * (vevSM(1)**2 - vevSM(2)**2)

  mat3a = M_L2 + 0.5_dp * (gSU2**2 + gU1**2) * vev2 * id3C

  !-----------
  ! couplings
  !-----------
   c_Sq4Z = ZeroC
   c_Sq4W = ZeroC
   c_Sq2Z = ZeroC
   c_Sq2W = ZeroC
#ifdef GENERATIONMIXING
   Rsf3 = Rsneut
   If (GenerationMixing) Then
    Do i1=1,3
     Do i2=1,3
      Call CoupSfermionZa(i2, i1, gSU2, sinW2, e, T3, RSf3, coupC1)
      Do i3=1,3
       Call CoupSfermionZa(i3, i1, gSU2, sinW2, e, T3, RSf3, coupC2)
       c_Sq4Z(i1,i2,i3) = coupC1 * Conjg(coupC2)
      End Do
     End Do
    End Do  
    Do i1=1,6
     Do i2=1,3
      Do i3=1,3
       c_Sq4W(i1,i2,i3) = Rsfp6(i1,i2) * Conjg( Rsfp6(i1,i3) ) 
      End Do
     End Do
    End Do  
    c_Sq4W = 0.5_dp * gSU2**2 * c_Sq4W

   Else ! .not. GenerationMixing
#endif
    Do i1=1,3
     c_Sq4Z(i1,i1,i1) = 0.25_dp * gSU2**2 / (1._dp - sinW2)
     Rsfp2 = Rsfp6(2*i1-1:2*i1,2*i1-1:2*i1)
     Do i2=1,2
      c_Sq4W(2*(i1-1)+i2,i1,i1) = Rsfp2(i2,1) * Conjg( Rsfp2(i2,1) ) 
     End Do
    End Do
    c_Sq4W = 0.5_dp * gSU2**2 * c_Sq4W
#ifdef GENERATIONMIXING
   End If ! GenerationMixing
#endif
   Call CoupSfermionZ(1, 1, gSU2, sinW2, 0._dp, T3, id2c, coupC1)
   c_Sq2Z = 4._dp * Abs( coupC1 )**2
   c_Sq2W = 2._dp * gSU2**2
   

  c_FNSf_R = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,4
     Do i3=1,3
      Call CoupNeutralinoSneutrino(i1, i2, i3, gU1, gSU2, N, id3C, id3C  &
                                 &, c_FNSf_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else ! .not. GenerationMixing
#endif
  Do i2=1,4
    Call CoupNeutralinoSneutrino(i2, gU1, gSU2, N, c_FNSf_R(1,i2,1) )
    c_FNSf_R(2,i2,2) = c_FNSf_R(1,i2,1)
    c_FNSf_R(3,i2,3) = c_FNSf_R(1,i2,1)
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_CFpSf_L = ZeroC
  c_CFpSf_R = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,2
    Do i2=1,3
     Do i3=1,3
      Call CoupCharginoSfermion(i1, i2, i3, gSU2, -T3, id3C, Y_l, mat3, uL_L  &
                      &, uL_R, U, V, c_CFpSf_L(i1,i2,i3), c_CFpSf_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else ! .not. GenerationMixing
#endif
   Do i1=1,2
    Do i2=1,3
      Call CoupCharginoSfermion(i1,1,gSU2, -T3, id2C, Y_L(i2,i2), ZeroC, U, V &
                             & , c_CFpSf_L(i1,i2,i2), c_CFpSf_R(i1,i2,i2) )
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_P0P0SfSf = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_P0
    Do i2=1,3
     Do i3=1,3
      Call CoupPseudoscalarSfermion4(i1, i1, i2, i3, RP0, T3, e, mat3  &
                                    &, id3C, gU1, gSU2, c_P0P0SfSf(i1,i2,i3) )
     End Do   
    End Do   
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,n_P0
    Call CoupPseudoscalarSfermion4(i1, i1, 1, 1, RP0, T3, e, ZeroC, id2C  &
              &, gU1, gSU2, c_P0P0SfSf(i1,1,1) )
    c_P0P0SfSf(i1,2,2) = c_P0P0SfSf(i1,1,1)
    c_P0P0SfSf(i1,3,3) = c_P0P0SfSf(i1,1,1)
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_S0SfSf = ZeroC
  c_S0S0SfSf = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_S0
    Do i2=1,3
     Do i3=1,3
      Call CoupScalarSfermion4(i1, i1, i2, i3, RS0, T3, e, mat3  &
                             &, id3C, gU1, gSU2, c_S0S0SfSf(i1,i2,i3) )
      Call CoupScalarSfermion3a(i1, i2, i3, RS0, T3, e, mat3, Rsf3, mat3, mu &
                              &, vevSM, gU1, gSU2, c_S0SfSf(i1,i2,i3) )
     End Do   
    End Do   
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,n_S0
    Call CoupScalarSfermion4(i1, i1, 1, 1, RS0, T3, e, ZeroC, id2C  &
                            &, gU1, gSU2, c_S0S0SfSf(i1,1,1) )
    c_S0S0SfSf(i1,2,2) = c_S0S0SfSf(i1,1,1)
    c_S0S0SfSf(i1,3,3) = c_S0S0SfSf(i1,1,1)
    Call CoupScalarSfermion3a(i1, 1, 1, RS0, T3, e, ZeroC, id2c, ZeroC &
                             &, mu, vevSM, gU1, gSU2, c_S0SfSf(i1,1,1) )
    c_S0SfSf(i1,2,2) = c_S0SfSf(i1,1,1)
    c_S0SfSf(i1,3,3) = c_S0SfSf(i1,1,1)
   End Do
#ifdef GENERATIONMIXING
  End If
#endif


  c_SpmSfSfp = ZeroC
  c_SpmSpmSfSfp = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_Spm
    Do i2=1,3
     Do i3=1,3
      Call CoupChargedScalarSfermion4(i1, i1, i2, i3, RSpm, T3, e, mat3, Y_l  &
                             &, id3C, gU1, gSU2, c_SpmSpmSfSfp(i1,i2,i3) )
     End Do   
     Do i3=1,6
      Call CoupChargedScalarSfermion3(i1, i3, i2, RSpm, gSU2, vevSM, mu, Y_l &
                 & , A_l, RSlepton, mat3, mat3, id3C,  c_SpmSfSfp(i1,i2,i3) )
     End Do   
    End Do   
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,n_Spm
    Do i2=1,3
     Call CoupChargedScalarSfermion4(i1, i1, 1, 1, RSpm, T3, e, ZeroC  &
                     &, Y_L(i2,i2), id2c, gU1, gSU2, c_SpmSpmSfSfp(i1,i2,i2) )
     Rsf2 = RSlepton(2*i2-1:2*i2,2*i2-1:2*i2)
     Do i4=1,2
      Call CoupChargedScalarSfermion3(i1, i4, 1, RSpm, gSU2, vevSM, mu     &
               & , Y_l(i2,i2), A_l(i2,i2), Rsf2, ZeroC, ZeroC, id2c &
               & , c_SpmSfSfp(i1,i2,2*(i2-1)+i4) )
     End Do   
    End Do   
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_SuSf4 = ZeroC
  c_SdSf4 = ZeroC
  c_SlSf4 = ZeroC
  c_SnSf4 = ZeroC

  ed = -1._dp / 3._dp
  eu = 2._dp / 3._dp
  e2 = -1._dp
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,6
    Do i2=1,3
     Do i3=1,3
      Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, T3, e, id3C, T3, eu &
                         &, RSup, 3, .False., c_SuSf4(i1,i2,i3) )
      Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, T3, e, id3C, -T3, ed &
                         &, RSdown, 3, .False., c_SdSf4(i1,i2,i3) )
      Call CoupSfermion4Y(i1, i1, i2, i3, -T3, Y_l, RSlepton, T3, mat3, id3c, 1 &
                         &, c_SlSf4(i1,i2,i3) )
      Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, T3, e, id6C, -T3, e2 &
                        &, RSlepton, 1, .True., coupC)
      c_SlSf4(i1,i2,i3) = c_SlSf4(i1,i2,i3) + coupC
     End Do
    End Do
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,3
    Do i2=1,3
     i_min = 2*i2-1
     i_max = i_min + 1
     Rsf2 = RSup(i_min:i_max,i_min:i_max)
     Rsfp2 = RSdown(i_min:i_max,i_min:i_max)
     Do i3=1,2
      i_min = 2*(i2-1)+i3
      Call CoupSfermion4G(1, 1, i3, i3, gU1, gSU2, T3, e, id2C &
                         &, T3, eu, RSf2, 3, .False., c_SuSf4(i_min, i1, i1))
      Call CoupSfermion4G(1, 1, i3, i3, gU1, gSU2, T3, e, id2C, -T3, ed &
                         &, RSfp2, 3, .False., c_SdSf4(i_min, i1, i1))
     End Do
    End Do
    i_min = 2*i1-1
    i_max = i_min + 1
    Call CoupSfermion4Y(1, 1, 1, 1, -T3, Y_l(i1,i1), RSfp2, T3 &
                      &, ZeroC, id2c, 1, c_SlSf4(i_min,i1,i1) )
    Call CoupSfermion4Y(2, 2, 1, 1, -T3, Y_l(i1,i1), RSfp2, T3 &
                      &, ZeroC, id2c, 1, c_SlSf4(i_max,i1,i1) )
    Do i2=1,3
     i_min = 2*i2-1
     i_max = i_min + 1
     Rsfp2 = RSlepton(i_min:i_max,i_min:i_max)
     Do i3=1,2
      i_min = 2*(i2-1)+i3
      If (i1.Eq.i2) Then
       Call CoupSfermion4G(1, 1, i3, i3, gU1, gSU2, T3, e, id2C, -T3, e2    &
                          & , RSfp2, 1, .True., coupC)
      Else
       Call CoupSfermion4G(1, 1, i3, i3, gU1, gSU2, T3, e, id2C, -T3, e2    &
                          & , RSfp2, 1, .False., coupC)
      End If
      c_SlSf4(i_min,i1,i1) = c_SlSf4(i_min,i1,i1) + coupC
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  e2 = 0._dp
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,3
     Call CoupSfermion4G(i2, i2, i1, i1, gU1, gSU2, e2, T3, Rsneut, 1 &
                       &, c_SnSf4(i1,i2,i2), 1 )
    End Do
   End Do

  Else ! .not. GenerationMixing
#endif
   Call CoupSfermion4G(1, 1, 1, 1, gU1, gSU2, e2, T3, id2c, 1&
                      &, c_SnSf4(1,1,1), 1 )
   c_SnSf4(2,2,2) = c_SnSf4(1,1,1)
   c_SnSf4(3,3,3) = c_SnSf4(1,1,1)
   Call CoupSfermion4G(1, 1, 1, 1, gU1, gSU2, T3, e2, id2c, T3, e2, id2c, 1&
                      &, .False., c_SnSf4(2,1,1) )
   c_SnSf4(1,2,2) = c_SnSf4(2,1,1)
   c_SnSf4(1,3,3) = c_SnSf4(2,1,1)
   c_SnSf4(2,3,3) = c_SnSf4(2,1,1)
   c_SnSf4(3,1,1) = c_SnSf4(2,1,1)
   c_SnSf4(3,2,2) = c_SnSf4(2,1,1)
#ifdef GENERATIONMIXING
  End If
#endif
  !-------------------
  ! loop contribution
  !-------------------
  Do i1=1,3
   PiSf(i1,:,:) = ZeroC
   p2 = mSneut2(i1)
   Call PiSneutrino(p2, mSneut2, c_Sq4Z, mSlepton2, c_Sq4W, c_Sq2Z, c_Sq2W    &
          & , mN2, c_FNSf_R, mC, mC2, c_CFpSf_L, c_CFpSf_R, mP02, c_P0P0SfSf  &
          & , mS02, c_S0SfSf, c_S0S0SfSf, mSpm2, c_SpmSfSfp, c_SpmSpmSfSfp    &
          & , mSup2, c_SuSf4, mSdown2, c_SdSf4, c_SlSf4, c_SnSf4, mZ2, mW2    &
          & , WriteOut, PiSf(i1,:,:))
  End Do

  Do i1=3,1,-1
   mat3 = mat3a - PiSf(i1,:,:)
   Call Chop(mat3) ! to avoid numerical problems with tiny numbers
#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Call EigenSystem(mat3,mi2,Rsf,kont, test)
    If ((kont.Eq.-14).Or.(kont.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(ErrCan,*) "test =",test
     Write(ErrCan,*) " "
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = 0
    End If
    If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
    Write(ErrCan,*) 'Diagonalization did not work in routine SneutrinoMass_1L!'
     Write(ErrCan,*) 'msf2 ', mi2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'mu, vevs, g, gU1'
     Write(ErrCan,*) mu, vevSM, gSU2, gU1
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    mass2(i1) = mi2(i1)
   Else ! .not.GenerationMixing
#endif
    mass2(i1) = Real( mat3(i1,i1),dp )
    Rsf = id3C
#ifdef GENERATIONMIXING
   End If ! GenerationMixing 
#endif
  End Do

  Do i1=1,3
   If (mass2(i1).Gt.0._dp) Then
    mass(i1) = Sqrt( mass2(i1) )
   Else 
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from routine SneutrinoMass_1L!'
     Write(ErrCan,*) 'in the calculation of the masses'
     Write(ErrCan,*) 'occurred a negative mass squared!!!'
     Write(ErrCan,*) 'msf2 ', mass2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'mu, vevs, g, gU1'
     Write(ErrCan,*) mu, vevSM, gSU2, gU1

     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    kont = -503
    Call AddError(503)
    mass(i1) = 0._dp
   End If
  End Do

  !-------------------------------------
  ! redoing calculation using refined p2
  !-------------------------------------
  i_count = 0
  Do 
   i_count = i_count + 1
   test_m2 = mass2

  Do i1=1,3
   PiSf(i1,:,:) = ZeroC
   p2 = mass2(i1)
   Call PiSneutrino(p2, mSneut2, c_Sq4Z, mSlepton2, c_Sq4W, c_Sq2Z, c_Sq2W    &
          & , mN2, c_FNSf_R, mC, mC2, c_CFpSf_L, c_CFpSf_R, mP02, c_P0P0SfSf  &
          & , mS02, c_S0SfSf, c_S0S0SfSf, mSpm2, c_SpmSfSfp, c_SpmSpmSfSfp    &
          & , mSup2, c_SuSf4, mSdown2, c_SdSf4, c_SlSf4, c_SnSf4, mZ2, mW2    &
          & , WriteOut, PiSf(i1,:,:))
  End Do

  Do i1=3,1,-1
   mat3 = mat3a - PiSf(i1,:,:)

#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Call Chop(mat3) ! to avoid numerical problems with tiny numbers
    Call EigenSystem(mat3,mi2,Rsf,kont, test)
    If ((kont.Eq.-14).Or.(kont.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(ErrCan,*) "test =",test
     Write(ErrCan,*) " "
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = 0
    End If
    If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
    Write(ErrCan,*) 'Diagonalization did not work in routine SneutrinoMass_1L!'
     Write(ErrCan,*) 'msf2 ', mi2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'mu, vevs, g, gU1'
     Write(ErrCan,*) mu, vevSM, gSU2, gU1
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    mass2(i1) = mi2(i1)
   Else ! .not.GenerationMixing
#endif
    mass2(i1) = Real( mat3(i1,i1),dp )
    Rsf = id3C
#ifdef GENERATIONMIXING
   End If ! GenerationMixing 
#endif
  End Do

  Do i1=1,3
    If (test_m2(i1).Ne.0._dp) Then
     test_m2(i1) = Abs(test_m2(i1) - mass2(i1)) / test_m2(i1)
    Else
     test_m2(i1) = Abs(mass2(i1))
    End If
   If (mass2(i1).Gt.0._dp) Then
    mass(i1) = Sqrt( mass2(i1) )
   Else 
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from routine SneutrinoMass_1L!'
     Write(ErrCan,*) 'in the calculation of the masses'
     Write(ErrCan,*) 'occurred a negative mass squared!!!'
     Write(ErrCan,*) 'msf2 ', mass2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'mu, vevs, g, gU1'
     Write(ErrCan,*) mu, vevSM, gSU2, gU1

     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    kont = -503
    Call AddError(503)
    mass(i1) = 0._dp
   End If
  End Do
   If (Maxval(test_m2).Lt. 0.1_dp*delta) Exit
   If (i_count.Gt.30) Then
    Write(*,*) "Problem in SneutrinoMassLoop",test_m2,mass2
    kont = - 504
     Call AddError(504)
    Exit
   End If
  End Do ! i_count


  Iname = Iname - 1

 End Subroutine SneutrinoMass_1L


#ifdef GENERATIONMIXING
 Subroutine SquarkMass_1L(g, T3, e, Y_l, Y_d, Y_u, vevSM, M_L2, M_R2, Af, Afp &
        & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mSup2, RSup, mSdown2, RSdown &
        & , mSlepton2, RSlepton, mSneut2, RSneut, mS02, RS0, mP02, RP0        &
        & , mSpm2, RSpm, mglu, phase_glu, uD_L, uD_R, uU_L, uU_R, mZ2, mW2    &
        & , delta, mass, mass2, Rsf, kont)
#else
 Subroutine SquarkMass_1L(g, T3, e, Y_l, Y_d, Y_u, vevSM, M_L2, M_R2, Af, Afp &
        & , mu_T, mu, mN, mN2, N, mC, mC2, U, V, mSup2, RSup, mSdown2, RSdown &
        & , mSlepton2, RSlepton, mSneut2, mS02, RS0, mP02, RP0, mSpm2, RSpm   &
        & , mglu, phase_glu, mZ2, mW2, delta, mass, mass2, Rsf, kont)
#endif
 !-------------------------------------------------------------------------
 ! calculates the 1-loop corrected squark masses + mixing angles for
 ! three generations but without generation mixing.
 ! input:
 ! output:
 !  - mass(i) ...... masses
 !  - mass2(i) ..... masses squared
 !  - Rsf(i,j) ..... mixing matrix
 ! written by Werner Porod, 28.12.01
 !  28.12.01: write first tree-level mass matrix as in case of generation
 !            mixing; check during the calculation of the couplings and
 !            diagonalisation whether generation mixing is included or not
 !  30.01.02: fixing bug in quatric QCD squark coupling 
 !-------------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: g(3), vevSM(2), T3, e, mN(4), mN2(4), mC(2), mC2(2) &
      &  , mSup2(6), mSdown2(6), mSlepton2(6), mSneut2(3), mS02(2), RS0(2,2)  &
      &  , mP02(2), RP0(2,2), mSpm2(2), mglu, mZ2, mW2, delta
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3), M_L2(3,3)         &
     & , M_R2(3,3), Af(3,3), Afp(3,3), mu, N(4,4), U(2,2), V(2,2), RSup(6,6) &
     & , RSdown(6,6), RSlepton(6,6), RSpm(2,2), phase_glu, mu_T
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: uD_L(3,3), uD_R(3,3), uU_L(3,3), uU_R(3,3)   &
     & , RSneut(3,3)
#endif
  Complex(dp), Intent(out) :: Rsf(6,6)
  Real(dp), Intent(out) :: mass(6), mass2(6)

  Integer :: i1, i2, i3, i4, i5, i_min, i_max, i_gen1, i_gen2, i_count
#ifdef GENERATIONMIXING
  Real(dp) :: mi6(6), test(2)
  Complex(dp) :: mat6(6,6), uF_L(3,3), uF_R(3,3), uFp_L(3,3), uFp_R(3,3)
#endif
  Real(dp) :: gU1, gSU2, gSU3, vev2, p2, Yl, Yr, msf2(6), mf2(3), mf(3)     &
     & , msfp2(6), mfp2(3), mfp(3), sinW2, e2, c_Sq2Z(6), c_Sq2W(6), mi2(2) &
     & , test_m2(6)
  Complex(dp) :: off(3,3), Yuk(3,3), YukC(3,3), YukT(3,3), mat6a(6,6)         &
     & , mat2(2,2), PiSf(6,6,6), Rsf2(2,2), Rsf6(6,6), Rsfp6(6,6), coupC1     &
     & , c_Sq4g3(6,6,6), c_Sq4e(6,6,6), c_Sq4Z(6,6,6), c_Sq4W(6,6,6), coupC2  &
     & , c_FNSf_L(3,4,6), c_FNSf_R(3,4,6), c_CFpSf_L(2,3,6), c_CFpSf_R(2,3,6) &
     & , c_GQSq(3,6,6), c_P0P0SqSq(2,6,6), c_P0SqSq(2,6,6), bi(1)             &
     & , c_S0S0SqSq(2,6,6), c_S0SqSq(2,6,6), c_SpmSpmSqSqp(2,6,6)             &
     & , c_SpmSqSqp(2,6,6), c_SuSf4(6,6,6), c_SdSf4(6,6,6), c_SlSf4(6,6,6)    &
     & , c_SnSf4(3,6,6), Rsfp2(2,2), coupC, c_Sq2g3(6,6,6)
  Logical :: WriteOut

  !------------------------------------
  ! initialisation
  !------------------------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'SquarkMass_1L'

  If ((WriteOneLoopContributions.Eq.12).Or.(WriteOneLoopContributions.Lt.0) ) &
   & Then
     Write(ErrCan,*) "Contributions to SquarkMass_1L:"
     WriteOut = .True.
  Else
   WriteOut = .False.
  End If

  PiSf = ZeroC
  mass = 0._dp
  mass2 = 0._dp
  RSf = ZeroC

  gU1 = g(1)
  gSU2 = g(2)
  gSU3 = g(3)
  sinW2 = gU1**2 / (gSU2**2 + gU1**2)
  If (T3.Gt.0._dp) Then
   Yuk = Y_u
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   msf2 = mSup2
   Rsf6 = RSup
   mf = mf_u
   mf2 = mf_u2
   msfp2 = mSdown2
   Rsfp6 = RSdown
   mfp = mf_d
   mfp2 = mf_d2
#ifdef GENERATIONMIXING
   uF_L = uU_L
   uF_R = uU_R
   uFp_L = uD_L
   uFp_R = uD_R
#endif

  Else ! T3 < 0
   Yuk = Y_d
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   msf2 = mSdown2
   Rsf6 = RSdown
   mf = mf_d
   mf2 = mf_d2
   msfp2 = mSup2
   Rsfp6 = RSup
   mfp = mf_u
   mfp2 = mf_u2
#ifdef GENERATIONMIXING
   uF_L = uD_L
   uF_R = uD_R
   uFp_L = uU_L
   uFp_R = uU_R
#endif
  End If

  !-------------------------
  ! tree level mass matrix
  !-------------------------
  vev2 = 0.25_dp * (vevSM(1)**2 - vevSM(2)**2)

  YukC = Conjg(Yuk)
  YukT = Transpose(Yuk)

  If (T3.Gt.0) Then
   mat6a(1:3,1:3) = M_L2 + 0.5_dp * vevSM(2)**2 * Matmul(YukC,YukT) &
        &         + (T3 * gSU2**2 - 0.5_dp * Yl * gU1**2) * vev2 * id3C
   mat6a(4:6,4:6) = M_R2 + 0.5_dp * vevSM(2)**2 * Matmul(YukT,YukC) &
        &         - 0.5_dp * Yr * gU1**2 * vev2 * id3C
    off = (vevSM(2) * Af - Conjg(mu) * vevSM(1) * Yuk ) * oosqrt2
  Else
   mat6a(1:3,1:3) = M_L2 + 0.5_dp * vevSM(1)**2 * Matmul(YukC,YukT) &
        &         + (T3 * gSU2**2 - 0.5_dp * Yl * gU1**2) * vev2 * id3C
   mat6a(4:6,4:6) = M_R2 + 0.5_dp * vevSM(1)**2 * Matmul(YukT,YukC) &
        &         - 0.5_dp * Yr * gU1**2 * vev2 * id3C
   off = (vevSM(1) * Af - Conjg(mu) * vevSM(2) * Yuk ) * oosqrt2
  End If

  off = Conjg(off)
  mat6a(1:3,4:6) = off
  Call Adjungate(off)
  mat6a(4:6,1:3) = off
    
  !-----------
  ! couplings
  !-----------
   c_GQSq = 0._dp
#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       c_GQSq(i1,i2,i3) = Conjg(uF_L(i1,i2)) * uF_L(i1,i3)
       c_GQSq(i1,i2+3,i3+3) = Conjg(uF_R(i1,i2)) * uF_R(i1,i3)
       c_GQSq(i1,i2,i3+3) = Conjg(uF_R(i1,i3))*uF_L(i1,i2) * phase_glu
       c_GQSq(i1,i2+3,i3) = Conjg(uF_L(i1,i3))*uF_R(i1,i2) * Conjg(phase_glu)
      End Do
     End Do
    End Do
    c_GQSq = c_GQSq * 8._dp * gSU3**2 / 3._dp

   Else ! .not. GenerationMixing
#endif
    c_GQSq(1,1,1) = 8._dp * gSU3**2 / 3._dp
    Do i1=2,6
      c_GQSq(1,i1,i1) = c_GQSq(1,1,1)
      c_GQSq(1,i1-1,i1) = c_GQSq(1,1,1) * phase_glu
      c_GQSq(1,i1,i1-1) = c_GQSq(1,1,1) * Conjg(phase_glu)
    End Do
    c_GQSq(2,:,:) = c_GQSq(1,:,:)
    c_GQSq(3,:,:) = c_GQSq(1,:,:)
#ifdef GENERATIONMIXING
   End If 
#endif

   c_Sq2g3 = ZeroC
   c_Sq4g3 = ZeroC
   c_Sq4e = ZeroC
   c_Sq4Z = ZeroC
   c_Sq2Z = ZeroC
   c_Sq4W = ZeroC
   c_Sq2W = ZeroC
#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Do i1=1,6
     Call CoupSfermionZ(i1, i1, gSU2, sinW2, e, T3, id6c, coupC1)
     c_Sq2Z(i1) = 4._dp * Abs( coupC1 )**2
     Do i2=1,6
      Call CoupSfermionZa(i2, i1, gSU2, sinW2, e, T3, RSf6, coupC1)
      Do i3=1,6
       c_Sq4g3(i1,i2,i3) = Rsf6(i1,i2) * Conjg( Rsf6(i1,i3) )
       c_Sq2g3(i1,i2,i3) = c_Sq4g3(i1,i2,i3)
       ! right squark belongs to SU(3)bar representation
       ! this gives in this case a sign
       If (i2.Gt.3) c_Sq4g3(i1,i2,i3) = - c_Sq4g3(i1,i2,i3)
       If (i3.Gt.3) c_Sq4g3(i1,i2,i3) = - c_Sq4g3(i1,i2,i3)
       Call CoupSfermionZa(i3, i1, gSU2, sinW2, e, T3, RSf6, coupC2)
       c_Sq4Z(i1,i2,i3) = coupC1 * Conjg(coupC2)
      End Do
     End Do
     Do i2=1,3
      Do i3=1,3
       c_Sq4W(i1,i2,i3) = Rsfp6(i1,i2) * Conjg( Rsfp6(i1,i3) ) 
      End Do
     End Do
    End Do  
    c_Sq2W(1) = 2._dp * gSU2**2
    c_Sq2W(2) = c_Sq2W(1)
    c_Sq2W(3) = c_Sq2W(1)

   Else ! .not. GenerationMixing
#endif
    Do i1=1,3
     Rsf2 = Rsf6(2*i1-1:2*i1,2*i1-1:2*i1)
     Rsfp2 = Rsfp6(2*i1-1:2*i1,2*i1-1:2*i1)
     Do i2=1,2
      Do i3=1,2
       i_gen1 = 2*(i1-1) + i3
       Call CoupSfermionZa(i3, i2, gSU2, sinW2, e, T3, RSf2, coupC1)
       Do i4=1,2
        i_gen2 = 2*(i1-1) + i4
        c_Sq4g3(2*(i1-1)+i2,i_gen1,i_gen2) = Rsf2(i2,i3) * Conjg( Rsf2(i2,i4) )
        c_Sq2g3(2*(i1-1)+i2,i_gen1,i_gen2) = Rsf2(i2,i3) * Conjg( Rsf2(i2,i4) )
        Call CoupSfermionZa(i4, i2, gSU2, sinW2, e, T3, RSf2, coupC2)
        c_Sq4Z(2*(i1-1)+i2,i_gen1,i_gen2) = coupC1 * Conjg(coupC2)
        ! right squark belongs to SU(3)bar representation
        ! this gives in this case a sign
        If (i3.Ne.i4) c_Sq4g3(2*(i1-1)+i2,i_gen1,i_gen2) =  &
                   & - c_Sq4g3(2*(i1-1)+i2,i_gen1,i_gen2)
       End Do
      End Do
     End Do
     i_gen1 = 2*(i1-1)
     i_gen2 = i_gen1 + 1
     Do i2=1,2
      c_Sq4W(i_gen1+i2,i_gen2,i_gen2) = Rsfp2(i2,1) * Conjg( Rsfp2(i2,1) )
     End Do
    End Do
    Call CoupSfermionZ(1, 1, gSU2, sinW2, e, T3, id2c, coupC1)
    c_Sq2Z(1) = 4._dp * Abs( coupC1 )**2
    c_Sq2Z(3) = c_Sq2Z(1)
    c_Sq2Z(5) = c_Sq2Z(1)
    Call CoupSfermionZ(2, 2, gSU2, sinW2, e, T3, id2c, coupC1)
    c_Sq2Z(2) = 4._dp * Abs( coupC1 )**2
    c_Sq2Z(4) = c_Sq2Z(2)
    c_Sq2Z(6) = c_Sq2Z(2)    
    c_Sq2W(1) = 2._dp * gSU2**2
    c_Sq2W(3) = c_Sq2W(1)
    c_Sq2W(5) = c_Sq2W(1)
#ifdef GENERATIONMIXING
   End If ! GenerationMixing
#endif
   c_Sq4e = c_Sq2g3 * gSU2**2 * e**2 * sinW2
   c_Sq2g3 = c_Sq2g3 * 4._dp * gSU3**2 / 3._dp
   c_Sq4g3 = c_Sq4g3 * 4._dp * gSU3**2 / 3._dp
   c_Sq4W = 0.5_dp * gSU2**2 * c_Sq4W

  c_FNSf_L = ZeroC
  c_FNSf_R = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,4
     Do i3=1,6
      Call CoupNeutralinoSfermion(i1, i2, i3, gU1, gSU2, T3, e, id6C, uf_L  &
                       &, uf_R, Yuk, N, c_FNSf_L(i1,i2,i3), c_FNSf_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else ! .not. GenerationMixing
#endif
   Do i1=1,3
    Do i2=1,4
     Do i3=1,2
      Call CoupNeutralinoSfermion(i2, i3, gU1, gSU2, T3, e, id2C, Yuk(i1,i1) &
               & , N, c_FNSf_L(i1,i2,2*(i1-1)+i3), c_FNSf_R(i1,i2,2*(i1-1)+i3))
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_CFpSf_L = ZeroC
  c_CFpSf_R = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,2
    Do i2=1,3
     Do i3=1,6
      Call CoupCharginoSfermion(i1, i2, i3, gSU2, -T3, id6C, Y_D, Y_U, ufp_L  &
                      &, ufP_R, U, V, c_CFpSf_L(i1,i2,i3), c_CFpSf_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else ! .not. GenerationMixing
#endif
   Do i1=1,2
    Do i2=1,3
     Do i3=1,2
      Call CoupCharginoSfermion(i1,i3,gSU2, -T3, id2C, Y_D(i2,i2), Y_U(i2,i2) &
          & , U, V, c_CFpSf_L(i1,i2,2*(i2-1)+i3), c_CFpSf_R(i1,i2,2*(i2-1)+i3))
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  bi(1) = mu_T
  c_P0SqSq = ZeroC
  c_P0P0SqSq = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_P0
    Do i2=1,6
     Do i3=1,6
      Call CoupPseudoscalarSfermion4(i1, i1, i2, i3, RP0, T3, e, yuk  &
                                    &, id6C, gU1, gSU2, c_P0P0SqSq(i1,i2,i3) )
      Call CoupPseudoScalarSfermion3a(i1, i2, i3, RP0, T3, yuk, RSf6, Af, bi &
                                    &, c_P0SqSq(i1,i2,i3) )
     End Do   
    End Do   
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,n_P0
    Do i2=1,3
     Rsf2 = Rsf6(2*i2-1:2*i2,2*i2-1:2*i2)      
     Do i3=1,2
      Do i4=1,2
       Call CoupPseudoscalarSfermion4(i1, i1, i3, i4, RP0, T3, e, yuk(i2,i2)  &
              &, id2C, gU1, gSU2, c_P0P0SqSq(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
       Call CoupPseudoScalarSfermion3a(i1, i3, i4, RP0, T3, yuk(i2,i2), Rsf2 &
                     &, Af(i2,i2), bi, c_P0SqSq(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
      End Do   
     End Do   
    End Do   
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_S0SqSq = ZeroC
  c_S0S0SqSq = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_S0
    Do i2=1,6
     Do i3=1,6
      Call CoupScalarSfermion4(i1, i1, i2, i3, RS0, T3, e, yuk  &
                             &, id6C, gU1, gSU2, c_S0S0SqSq(i1,i2,i3) )
      Call CoupScalarSfermion3a(i1, i2, i3, RS0, T3, e, yuk, RSf6, Af, mu_T &
                              &, vevSM, gU1, gSU2, c_S0SqSq(i1,i2,i3) )
     End Do   
    End Do   
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,n_S0
    Do i2=1,3
     Rsf2 = Rsf6(2*i2-1:2*i2,2*i2-1:2*i2)      
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion4(i1, i1, i3, i4, RS0, T3, e, yuk(i2,i2)  &
              &, id2C, gU1, gSU2, c_S0S0SqSq(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
       Call CoupScalarSfermion3a(i1, i3, i4, RS0, T3, e, yuk(i2,i2), Rsf2 &
                                &, Af(i2,i2), mu_T, vevSM, gU1, gSU2      &
                                &, c_S0SqSq(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
      End Do   
     End Do   
    End Do   
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_SpmSqSqp = ZeroC
  c_SpmSpmSqSqp = ZeroC
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,n_Spm
    Do i2=1,6
     Do i3=1,6
      If (T3.Gt.0._dp) Then
       Call CoupChargedScalarSfermion4(i1, i1, i2, i3, RSpm, T3, e, Y_u  &
                            &, Y_d, id6C, gU1, gSU2, c_SpmSpmSqSqp(i1,i2,i3) )
       Call CoupChargedScalarSfermion3(i1, i3, i2, RSpm, gSU2, vevSM, mu_T  &
                 & , Y_d, Afp, RSdown, Y_u, Af, id6C,  c_SpmSqSqp(i1,i2,i3) )
      Else
       Call CoupChargedScalarSfermion4(i1, i1, i2, i3, RSpm, T3, e, Y_d  &
                            &, Y_u, id6C, gU1, gSU2, c_SpmSpmSqSqp(i1,i2,i3) )
       Call CoupChargedScalarSfermion3(i1, i2, i3, RSpm, gSU2, vevSM, mu_T  &
                   & , Y_d, Af, id6c, Y_u, Afp, RSup,  c_SpmSqSqp(i1,i2,i3) )
      End If
     End Do   
    End Do   
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,n_Spm
    Do i2=1,3
     Rsf2 = RSfp6(2*i2-1:2*i2,2*i2-1:2*i2)
     Do i3=1,2
      Do i4=1,2
       If (T3.Gt.0._dp) Then
        Call CoupChargedScalarSfermion4(i1, i1, i3, i4, RSpm, T3, e  &
              &, Y_u(i2,i2), Y_d(i2,i2), id2C, gU1, gSU2 &
              &, c_SpmSpmSqSqp(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
        Call CoupChargedScalarSfermion3(i1, i4, i3, RSpm, gSU2, vevSM, mu_T  &
               & , Y_d(i2,i2), Afp(i2,i2), Rsf2, Y_u(i2,i2), Af(i2,i2), id2c &
               & , c_SpmSqSqp(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
       Else
        Call CoupChargedScalarSfermion4(i1, i1, i3, i4, RSpm, T3, e  &
              &, Y_d(i2,i2), Y_u(i2,i2), id2C, gU1, gSU2 &
              &, c_SpmSpmSqSqp(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
        Call CoupChargedScalarSfermion3(i1, i3, i4, RSpm, gSU2, vevSM, mu_T  &
               & , Y_d(i2,i2), Af(i2,i2), id2c, Y_u(i2,i2), Afp(i2,i2), Rsf2 &
               & , c_SpmSqSqp(i1,2*(i2-1)+i3,2*(i2-1)+i4) )
       End If
      End Do   
     End Do   
    End Do   
   End Do
#ifdef GENERATIONMIXING
  End If
#endif

  c_SuSf4 = ZeroC
  c_SdSf4 = ZeroC
  c_SlSf4 = ZeroC
  c_SnSf4 = ZeroC
  If (T3.Gt.0._dp) Then
   e2 = - 1._dp / 3._dp
#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Do i1=1,6
     Do i2=1,6
      Do i3=1,6
       Call CoupSfermion4Y(i2, i3, i1, i1, Yuk, Rsup, 3, c_SuSf4(i1,i2,i3), 1)
       Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, e, T3, RSup, 3, coupC, 1)
       c_SuSf4(i1,i2,i3) = c_SuSf4(i1,i2,i3) + coupC
       Call CoupSfermion4Y(i1, i1, i2, i3, -T3, Y_d, RSdown, T3, Y_u, id6c, 3 &
                          &,c_SdSf4(i1,i2,i3) )
       Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, T3, e, id6C, -T3, e2 &
                         &, RSdown, 3, .True., coupC)
       c_SdSf4(i1,i2,i3) = c_SdSf4(i1,i2,i3) + coupC
       Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, T3, e, id6C, -T3,-1._dp &
                            &, RSlepton, 1, .False., c_SlSf4(i1,i2,i3) )
      End Do
     End Do
    End Do
   Else ! .not. GenerationMixing
#endif
    Do i1=1,3
     i_min = 2*i1-1
     i_max = i_min + 1
     Rsf2 = Rsf6(i_min:i_max,i_min:i_max)
     Rsfp2 = Rsfp6(i_min:i_max,i_min:i_max)
     Do i2=1,2
      i_min = 2*(i1-1)+i2
      Do i3=1,3
       Do i4=1,2
        i_gen1 = 2*(i3-1)+i4
        Do i5=1,2
         i_gen2 = 2*(i3-1)+i5
         If (i1.Eq.i3) Then
          Call CoupSfermion4Y(i4, i5, i2, i2, Yuk(i1,i1), Rsf2, 3  &
                           &, c_SuSf4(i_min,i_gen1,i_gen2), 1)
          Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, e, T3, RSf2,3,coupC,1)
          c_SuSf4(i_min,i_gen1,i_gen2) = c_SuSf4(i_min,i_gen1,i_gen2) + coupC
          Call CoupSfermion4Y(i2, i2, i4, i5, -T3, Y_d(i1,i1), RSfp2, T3 &
                           &, Y_u(i3,i3), id2c, 3, c_SdSf4(i_min,i_gen1,i_gen2) )
          Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, T3, e, id2C &
                            &, -T3, e2, RSfp2, 3, .True., coupC)
          c_SdSf4(i_min,i_gen1,i_gen2) = c_SdSf4(i_min,i_gen1,i_gen2) + coupC
         Else
          Call CoupSfermion4Y(i2, i2, i4, i5, T3, Y_u(i1,i1), RSf2, T3 &
                         &, Y_u(i3,i3), id2c, 3, c_SuSf4(i_min,i_gen1,i_gen2) )
          Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, T3, e, id2C &
                          &, T3, e, RSf2, 3, .False., coupC)
          c_SuSf4(i_min,i_gen1,i_gen2) = c_SuSf4(i_min,i_gen1,i_gen2) + coupC
          Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, T3, e, id2C &
                  &, -T3, e2, RSfp2, 3, .False., c_SdSf4(i_min,i_gen1,i_gen2) )
         End If
        End Do
       End Do
      End Do
     End Do
     i_min = 2*i1-1
     i_max = i_min + 1
     Rsfp2 = RSlepton(i_min:i_max,i_min:i_max)
     Do i2=1,2
      i_min = 2*(i1-1)+i2
      Do i3=1,3
       Do i4=1,2
        i_gen1 = 2*(i3-1)+i4
        Do i5=1,2
         i_gen2 = 2*(i3-1)+i5
         Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, T3, e, id2C, -T3    &
                 & , -1._dp, RSfp2, 1, .False.,c_SlSf4(i_min,i_gen1,i_gen2) )
        End Do
       End Do
      End Do
     End Do
    End Do
#ifdef GENERATIONMIXING
   End If
#endif

  Else ! T3 < 0
   e2 = 2._dp / 3._dp
#ifdef GENERATIONMIXING
   If (GenerationMixing) Then
    Do i1=1,6
     i_gen1 = (i1+1)/2
     Do i2=1,6
      i_gen2 = (i2+1)/2
      Do i3=1,6
       Call CoupSfermion4Y(i2, i3, i1, i1, Yuk, Rsdown, 3, c_SdSf4(i1,i2,i3), 1)
       Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, e, T3, RSf6, 3, coupC, 1)
       c_SdSf4(i1,i2,i3) = c_SdSf4(i1,i2,i3) + coupC
       Call CoupSfermion4Y(i1, i1, i2, i3, -T3, Y_u, RSup, T3, Y_d, id6c, 3 &
                          &,c_SuSf4(i1,i2,i3) )
       Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, T3, e, id6C, -T3, e2 &
                         &, RSup, 3, .True., coupC)
       c_SuSf4(i1,i2,i3) = c_SuSf4(i1,i2,i3) + coupC
       Call CoupSfermion4Y(i1, i1, i2, i3, T3, Y_l, RSlepton, T3, Y_d, id6c, 1 &
                          &,c_SlSf4(i1,i2,i3) )
       Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, T3, e, id6C, T3, -1._dp &
                            &, RSlepton, 1, .False., coupC )
       c_SlSf4(i1,i2,i3) = c_SlSf4(i1,i2,i3) + coupC
      End Do
     End Do
    End Do
   Else ! .not. GenerationMixing
#endif
    Do i1=1,3
     i_min = 2*i1-1
     i_max = i_min + 1
     Rsf2 = Rsf6(i_min:i_max,i_min:i_max)
     Rsfp2 = Rsfp6(i_min:i_max,i_min:i_max)
     Do i2=1,2
      i_min = 2*(i1-1)+i2
      Do i3=1,3
       Do i4=1,2
        i_gen1 = 2*(i3-1)+i4
        Do i5=1,2
         i_gen2 = 2*(i3-1)+i5
         If (i1.Eq.i3) Then
          Call CoupSfermion4Y(i4, i5, i2, i2, Yuk(i1,i1), Rsf2, 3 &
                            &, c_SdSf4(i_min,i_gen1,i_gen2), 1)
          Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, e, T3, RSf2,3,coupC,1)
          c_SdSf4(i_min,i_gen1,i_gen2) = c_SdSf4(i_min,i_gen1,i_gen2) + coupC
          Call CoupSfermion4Y(i2, i2, i4, i5, -T3, Y_u(i1,i1), RSfp2, T3 &
                            &, Y_d(i3,i3), id2c, 3, c_SuSf4(i_min,i_gen1,i_gen2) )
          Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, T3, e, id2C &
                            &, -T3, e2, RSfp2, 3, .True., coupC)
          c_SuSf4(i_min,i_gen1,i_gen2) = c_SuSf4(i_min,i_gen1,i_gen2) + coupC
         Else
          Call CoupSfermion4Y(i2, i2, i4, i5, T3, Y_d(i1,i1), RSf2, T3 &
                         &, Y_d(i3,i3), id2c, 3, c_SdSf4(i_min,i_gen1,i_gen2) )
          Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, T3, e, id2C &
                          &, T3, e, RSf2, 3, .False., coupC)
          c_SdSf4(i_min,i_gen1,i_gen2) = c_SdSf4(i_min,i_gen1,i_gen2) + coupC
          Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, T3, e, id2C &
                 &, -T3, e2, RSfp2, 3, .False., c_SuSf4(i_min,i_gen1,i_gen2) )
         End If
        End Do
       End Do
      End Do
     End Do
     i_min = 2*i1-1
     i_max = i_min + 1
     Rsfp2 = RSlepton(i_min:i_max,i_min:i_max)
     Do i2=1,2
      i_min = 2*(i1-1)+i2
      Do i3=1,3
       Do i4=1,2
        i_gen1 = 2*(i3-1)+i4
        Do i5=1,2
         i_gen2 = 2*(i3-1)+i5
         Call CoupSfermion4Y(i2, i2, i4, i5, T3, Y_l(i1,i1), RSfp2, T3 &
                           &, Y_d(i3,i3), id2c, 1, c_SlSf4(i_min,i_gen1,i_gen2) )
         Call CoupSfermion4G(i4, i5, i2, i2, gU1, gSU2, T3, e, id2C, T3    &
                 & , -1._dp, RSfp2, 1, .False., coupC )
         c_SlSf4(i_min,i_gen1,i_gen2) = c_SlSf4(i_min,i_gen1,i_gen2) + coupC
        End Do
       End Do
      End Do
     End Do
    End Do
   End If
#ifdef GENERATIONMIXING
  End If
#endif

  e2 = 0._dp
#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,6
     Do i3=1,6
      Call CoupSfermion4G(i2, i3, i1, i1, gU1, gSU2, T3, e, id6C, 0.5_dp, e2 &
                         &, RSneut, 1, .False., c_SnSf4(i1,i2,i3) )
     End Do
    End Do
   End Do

  Else ! .not. GenerationMixing
#endif
   Do i1=1,3
    i_min = i1
    i_max = i_min+1
    Rsf2 = Rsf6(i_min:i_max,i_min:i_max)
    Do i2=1,3
     Do i3=1,2
      i_gen1 = 2*(i2-1)+i3
      Do i4=1,2
       i_gen2 = 2*(i2-1)+i4
       Call CoupSfermion4G(i3, i4, 1, 1, gU1, gSU2, T3, e, id2C, 0.5_dp, e2 &
                         &, id2c, 1, .False., c_SnSf4(i1,i_gen1,i_gen2))
      End Do
     End Do
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If
#endif
  !-------------------
  ! loop contribution
  !-------------------
  Do i1=1,6
   If (T3.Gt.0._dp) Then
    p2 = msup2(i1)
   Else
    p2 = mSdown2(i1)
   End If 
   Call PiSquark(p2, mSf2, mf, mf2, mglu, c_GQSq, c_Sq4g3, c_Sq4e, c_Sq4Z     &
          & , c_Sq2Z, c_Sq2W, c_Sq2g3, mSfp2, c_Sq4W, mfp, mfp2               &
          & , mN, mN2, c_FNSf_L, c_FNSf_R, mC, mC2, c_CFpSf_L, c_CFpSf_R      &
          & , mP02, c_P0SqSq, c_P0P0SqSq, mS02, c_S0SqSq, c_S0S0SqSq          &
          & , mSpm2, c_SpmSqSqp, c_SpmSpmSqSqp, mSup2, c_SuSf4                &
          & , mSdown2, c_SdSf4, mSlepton2, c_SlSf4, mSneut2, c_SnSf4, mZ2, mW2&
          & , WriteOut, PiSf(i1,:,:))   
  End Do

#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=6,1,-1
    mat6 = mat6a - PiSf(i1,:,:)
    Call Chop(mat6) ! to avoid numerical problems
    Call EigenSystem(mat6,mi6,Rsf,kont,test)
    If ((kont.Eq.-14).Or.(kont.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(errcan,*) "t3",t3 
     Write(ErrCan,*) "test =",test
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = 0
    End If
    If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
     Write(ErrCan,*) 'Diagonalization did not work in routine SquarkMass_1L!'
     Write(ErrCan,*) 'msf2 ', mass2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'M_R2, ',M_R2
     Write(ErrCan,*) 'A_f, ',Af
     Write(ErrCan,*) 'Y_f, ',Yuk
     Write(ErrCan,*) 'mu, vevs, T3, Yl, Yr, g, gU1'
     Write(ErrCan,*) mu, vevSM, T3, Yl, Yr, gSU2, gU1
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    mass2(i1) = mi6(i1)
   End Do

  Else ! .not.GenerationMixing
#endif
   Do i1=1,3
    Do i2=2,1,-1
     mat2(1,1) = mat6a(i1,i1) - PiSf(2*(i1-1)+i2,2*i1-1,2*i1-1)
     mat2(2,2) = mat6a(3+i1,3+i1) - PiSf(2*(i1-1)+i2,2*i1,2*i1)
     mat2(1,2) = mat6a(i1,3+i1) - PiSf(2*(i1-1)+i2,2*i1-1,2*i1)
     mat2(2,1) = mat6a(3+i1,i1) - PiSf(2*(i1-1)+i2,2*i1,2*i1-1)
     If (WriteOut) Write(ErrCan,*) "tree(11,12,22)",i1,mat6a(i1,i1) &
   &         , mat6a(i1,3+i1), mat6a(3+i1,3+i1)
     If (WriteOut) Write(ErrCan,*) "pi(11,12,22)"   &
   &         , i1,PiSf(2*(i1-1)+i2,2*i1-1,2*i1-1) &
   &         , PiSf(2*(i1-1)+i2,2*i1-1,2*i1), PiSf(2*(i1-1)+i2,2*i1,2*i1)
     Call Diag2(mat2,mi2,Rsf2)
     If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
      Write(ErrCan,*) 'Diagonalization did not work in routine SquarkMass_1L'
      Write(ErrCan,*) 'msf2 ', mass2
      Write(ErrCan,*) 'M_L2, ',M_L2
      Write(ErrCan,*) 'M_R2, ',M_R2
      Write(ErrCan,*) 'A_f, ',Af
      Write(ErrCan,*) 'Y_f, ',Yuk
      Write(ErrCan,*) 'mu, vevs, T3, Yl, Yr, g, gU1'
      Write(ErrCan,*) mu, vevSM, T3, Yl, Yr, gSU2, gU1
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
     mass2(2*(i1-1)+i2) = mi2(i2)
     Rsf(2*i1-1:2*i1,2*i1-1:2*i1) = Rsf2
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If ! GenerationMixing 
#endif

  Do i1=1,6
   If (mass2(i1).Gt.0._dp) Then
    mass(i1) = Sqrt( mass2(i1) )
   Else 
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from routine SquarkMass_1L!'
     Write(ErrCan,*) 'in the calculation of the masses'
     Write(ErrCan,*) 'occurred a negative mass squared!!!'
     Write(ErrCan,*) 'msf2 ', mass2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'M_R2, ',M_R2
     Write(ErrCan,*) 'A_f, ',Af
     Write(ErrCan,*) 'Y_f, ',Yuk
     Write(ErrCan,*) 'mu, vevs, T3, Yl, Yr, g, gU1'
     Write(ErrCan,*) mu, vevSM, T3, Yl, Yr, gSU2, gU1

     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    kont = - 505
    Call AddError(505)
    mass(i1) = 0._dp
   End If
  End Do

  !-------------------------------------
  ! redoing calculation using refined p2
  !-------------------------------------
  i_count = 0
  Do 
   i_count = i_count + 1
   test_m2 = mass2
   Do i1=1,6
    p2 = mass2(i1)
    Call PiSquark(p2, mSf2, mf, mf2, mglu, c_GQSq, c_Sq4g3, c_Sq4e, c_Sq4Z    &
          & , c_Sq2Z, c_Sq2W, c_Sq2g3, mSfp2, c_Sq4W, mfp, mfp2               &
          & , mN, mN2, c_FNSf_L, c_FNSf_R, mC, mC2, c_CFpSf_L, c_CFpSf_R      &
          & , mP02, c_P0SqSq, c_P0P0SqSq, mS02, c_S0SqSq, c_S0S0SqSq          &
          & , mSpm2, c_SpmSqSqp, c_SpmSpmSqSqp, mSup2, c_SuSf4                &
          & , mSdown2, c_SdSf4, mSlepton2, c_SlSf4, mSneut2, c_SnSf4, mZ2, mW2&
          & , WriteOut, PiSf(i1,:,:))   
   End Do

#ifdef GENERATIONMIXING
  If (GenerationMixing) Then
   Do i1=6,1,-1
    mat6 = mat6a - PiSf(i1,:,:)
    Call Chop(mat6) ! to avoid numerical problems
    Call EigenSystem(mat6,mi6,Rsf,kont,test)
    If ((kont.Eq.-14).Or.(kont.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(errcan,*) "t3",t3 
     Write(ErrCan,*) "test =",test
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = 0
    End If
    If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
     Write(ErrCan,*) 'Diagonalization did not work in routine SquarkMass_1L!'
     Write(ErrCan,*) 'msf2 ', mass2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'M_R2, ',M_R2
     Write(ErrCan,*) 'A_f, ',Af
     Write(ErrCan,*) 'Y_f, ',Yuk
     Write(ErrCan,*) 'mu, vevs, T3, Yl, Yr, g, gU1'
     Write(ErrCan,*) mu, vevSM, T3, Yl, Yr, gSU2, gU1
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    mass2(i1) = mi6(i1)
   End Do

  Else ! .not.GenerationMixing
#endif
   Do i1=1,3
    Do i2=2,1,-1
     mat2(1,1) = mat6a(i1,i1) - PiSf(2*(i1-1)+i2,2*i1-1,2*i1-1)
     mat2(2,2) = mat6a(3+i1,3+i1) - PiSf(2*(i1-1)+i2,2*i1,2*i1)
     mat2(1,2) = mat6a(i1,3+i1) - PiSf(2*(i1-1)+i2,2*i1-1,2*i1)
     mat2(2,1) = mat6a(3+i1,i1) - PiSf(2*(i1-1)+i2,2*i1,2*i1-1)
     If (WriteOut) Write(ErrCan,*) "tree(11,12,22)",i1,mat6a(i1,i1) &
   &         , mat6a(i1,3+i1), mat6a(3+i1,3+i1)
     If (WriteOut) Write(ErrCan,*) "pi(11,12,22)"   &
   &         , i1,PiSf(2*(i1-1)+i2,2*i1-1,2*i1-1) &
   &         , PiSf(2*(i1-1)+i2,2*i1-1,2*i1), PiSf(2*(i1-1)+i2,2*i1,2*i1)
     Call Diag2(mat2,mi2,Rsf2)
     If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
      Write(ErrCan,*) 'Diagonalization did not work in routine SquarkMass_1L'
      Write(ErrCan,*) 'msf2 ', mass2
      Write(ErrCan,*) 'M_L2, ',M_L2
      Write(ErrCan,*) 'M_R2, ',M_R2
      Write(ErrCan,*) 'A_f, ',Af
      Write(ErrCan,*) 'Y_f, ',Yuk
      Write(ErrCan,*) 'mu, vevs, T3, Yl, Yr, g, gU1'
      Write(ErrCan,*) mu, vevSM, T3, Yl, Yr, gSU2, gU1
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
     mass2(2*(i1-1)+i2) = mi2(i2)
     Rsf(2*i1-1:2*i1,2*i1-1:2*i1) = Rsf2
    End Do
   End Do
#ifdef GENERATIONMIXING
  End If ! GenerationMixing 
#endif

  Do i1=1,6
    If (test_m2(i1).Ne.0._dp) Then
     test_m2(i1) = Abs(test_m2(i1) - mass2(i1)) / test_m2(i1)
    Else
     test_m2(i1) = Abs(mass2(i1))
    End If
   If (mass2(i1).Gt.0._dp) Then
    mass(i1) = Sqrt( mass2(i1) )
   Else 
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from routine SquarkMass_1L!'
     Write(ErrCan,*) 'in the calculation of the masses'
     Write(ErrCan,*) 'occurred a negative mass squared!!!'
     Write(ErrCan,*) 'msf2 ', mass2
     Write(ErrCan,*) 'M_L2, ',M_L2
     Write(ErrCan,*) 'M_R2, ',M_R2
     Write(ErrCan,*) 'A_f, ',Af
     Write(ErrCan,*) 'Y_f, ',Yuk
     Write(ErrCan,*) 'mu, vevs, T3, Yl, Yr, g, gU1'
     Write(ErrCan,*) mu, vevSM, T3, Yl, Yr, gSU2, gU1

     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    kont = - 505
    Call AddError(505)
    mass(i1) = 0._dp
   End If
  End Do
   If (Maxval(test_m2).Lt. 0.1_dp*delta) Exit
   If (i_count.Gt.30) Then
    Write(*,*) "Problem in SquarkMassLoop",test_m2,mass2
     kont = - 506
     Call AddError(506)
    Exit
   End If
  End Do ! i_count


  Iname = Iname - 1

 End Subroutine SquarkMass_1L


 Subroutine WriteLoopMassesError(n, name, kont)
 !------------------------------------------------------------------------
 ! This Subroutine Contains and writes the error messages for this Module
 ! In addition it sets the variable kont and in Case it is necessary
 ! it terminates the Program
 ! written by Werner Porod
 ! 08.10.02: - implementation for mu at tree level
 !------------------------------------------------------------------------
 Implicit None
 Integer, Intent(in) :: n
 Character(len=*), Intent(in) :: name
 Integer, Intent(inout) :: kont
 !------------------------------------
 ! in these cases Do nothing
 !------------------------------------
 If (ErrorLevel.Eq.-2) Return
 
 Select Case(n)
 Case(-9)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "m_h0^2 <0 in the 1-loop calculation", WriteWert
  kont = -507
  Call AddError(507)

 Case(-8)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "m_H+^2 <0 in the 1-loop calculation", WriteWert
  kont = -508
  Call AddError(507)

 Case(-7)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "m_A0^2 <0 in the 1-loop calculation", WriteWert
  kont = -509
  Call AddError(509)

 Case(-6)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "|mu|^2 > 10^20 in the 1-loop calculation", WriteWert
  kont = -510
  Call AddError(510)

 Case(-5)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "|mu|^2 < 0 in the 1-loop calculation", WriteWert
  kont = -511
  Call AddError(511)

 Case(-4)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "mZ^2(mZ) < 0",WriteWert
  kont = -512
  Call AddError(512)

 Case(1)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "|mu|^2 < 0 at tree level,  reversing sign",WriteWert

 Case(2)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "|mu|^2 > 10^20, setting it to 10^4",WriteWert

 Case(4)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "mZ^2(mZ) < 0 in the first run, setting it to mZ^2",WriteWert

 Case(5)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "|mu|^2 < 0 in the first 1-loop calc., reversing sign" &
              &, WriteWert

 Case(6)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) &
    &"|mu|^2  10^20 in the first 1-loop calc., setting it to 10^4", WriteWert

 Case(7)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "m_A0^2 <0 in the first 1-loop calculation", WriteWert
  Write(ErrCan,*) "Setting it to 10^4" 

 Case(8)
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "m_H+^2 <0 in the first 1-loop calculation", WriteWert
  Write(ErrCan,*) "Setting it to 10^4" 

 Case Default
  Write(ErrCan,*) "Problem in unit "//name
  Write(ErrCan,*) "An unknown error has occured :",n
  If  (ErrorLevel.Ge.0) Call TerminateProgram
 End Select
 If ((n.Lt.0).And.(Index(name,"_2").Gt.0)) Then
  kont = kont - 6
  Call AddError(-kont)
 End If
 If ((n.Lt.0).And.(Index(name,"_3").Gt.0)) Then
  kont = kont - 12
  Call AddError(-kont)
 End If
 !---------------------------------------
 ! in these cases the Program terminates
 !---------------------------------------
 If  ( (ErrorLevel.Eq.1).And.(n.Le.0) ) Call TerminateProgram
 If  (ErrorLevel.Eq.2) Call TerminateProgram

 End Subroutine WriteLoopMassesError



End Module LoopMasses
