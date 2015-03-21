Module SugraRuns

! load modules
 Use Control
 Use Mathematics, Only: Li2, odeint, odeintB, odeintC
 Use LoopFunctions
 Use RGEs
 Use StandardModel
 Use LoopCouplings
 Use Model_Data, Only: MnuR, in_extpar, r_extpar, i_extpar &
     & , product_stop_masses, CalcScale_from_Stops
 Use MSSM_Data, Only: P0
 Use LoopMasses
! load modules

! interfaces
 Interface BoundaryEW
  Module Procedure BoundaryEWnew, BoundaryEWold
 End Interface
! end interfaces

! global variables
Real(dp), save :: n5plets
Integer, save :: n10plets, YukScen
Real(dp), Save :: Lambda, MlambdaS, F_GMSB
! string scenarios
Integer, save :: num_t
Integer, save :: nH1_ai(36), nH2_ai(36), nE_ai(36,3), nL_ai(36,3)          &
     &   , nD_ai(36,3), nU_ai(36,3), nQ_ai(36,3), mH1_ai(36), mH2_ai(36)   &
     &   , mE_ai(36,3), mL_ai(36,3), mD_ai(36,3), mU_ai(36,3)              &
     &   , mQ_ai(36,3), mtH1_ai(36), mtH2_ai(36), mtE_ai(36,3)             &
     &   , mtL_ai(36,3), mtD_ai(36,3), mtU_ai(36,3), mtQ_ai(36,3)          &
     &   , n_ai(3,17), pLE(36,3,3), pEH1(36,3), pLH1(36,3), pQD(36,3,3)    &
     &   , pDH1(36,3), pQH1(36,3), pQU(36,3,3), pUH2(36,3), pQH2(36,3)     &
     &   , mG(36,3), pGE(36,3,3), pGL(36,3,3), pGD(36,3,3), pGU(36,3,3)    &
     &   , pGQ(36,3,3), pGH1(36,3), pGH2(36,3)
Real(dp), Save :: C_a(3), C_ai(3,20),  SumC_O1(3)
Real(dp), save :: m32, g_s, k_s, k_sb, k_ss, ThetaA(36)                     &
     &   , cosT, sumC1(3), sumC2(3), LnReDedekind(36), ReT(36)        &
     &   , oosqrt_k_ss, sinT, g_s2, delta_GS, sinT2, cosT2
Complex(dp), save :: phase_s, t(36), phase_t(36), G2t(36),ReG2ThetaT(36)
Real(dp), save :: mGUT_save, sinW2_DR_mZ      &
    & , mf_l_DR_SM(3), mf_d_DR_SM(3), mf_u_DR_SM(3)
Real(dp) :: m_32, M0_amsb
Complex(dp), save :: Yl_mZ(3,3), Yu_mZ(3,3), Yd_mZ(3,3)
Real(dp), Save :: vevs_DR_save(2)
! flag to decide if nu-Yukawas are given from outside or not
Logical, Save :: Fixed_Nu_Yukawas = .False., Ynu_eq_Yu = .False.
Logical, Save :: Off_GMSB=.False.
! flag to decide that nu-Yukawas are given at m_R_3
Logical, Save :: Ynu_at_MR3 = .False.
! check if in the super CKM basis are given
Logical, Save :: l_Au=.False., l_Ad=.False., l_Al=.False., l_MD = .False. &
    & , l_MQ = .False., l_MU = .False., l_ME = .False. , l_ML = .False. 
! quick and dirty way to implement model by Suchita Kulkarni
Logical, Save :: Model_Suchita = .False.      

! global variables

! private variables
Integer, Private, save :: YukawaScheme=1
Logical, Private, save :: CheckSugraDetails(10) = .False. &
                & , SugraErrors(10) = .False.        &
                & , StrictUnification = .False.      &
                & , UseFixedScale = .False.          &
                & , UseFixedGUTScale = .False.
Real(dp), Private, Save :: GUT_scale
Character(len=15), Private, save :: HighScaleModel
! private variables

Contains


 Subroutine BoundaryEWold(i_run, vevSM, mC, U, V, mN, N, mS02, RS0, mP02, RP0 &
    & , mSpm, mSpm2, RSpm, mDsquark, mDsquark2, RDsquark, mUsquark, mUsquark2 &
    & , RUsquark, mSlepton, mSlepton2, RSlepton, mSneutrino2, RSneutrino      &
    & , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R, Unu_L, mGlu, phase_glu, mZ2_run   &
    & , mW2_run, delta0, g1, kont)
 !-----------------------------------------------------------------------
 ! Calculates gauge and yukawa couplings at m_Z in the DRbar scheme
 ! written by Werner Porod
 ! 15.11.01: - dummy version to get consistent program structure
 ! 16.11.01: - implementing gauge couplings
 !             please note that the call to PiZZ1 has to be changed for
 !             complex neutral scalar
 !           - implementing first steps for Yukawa couplings including
 !             an iterative process due to bottom Yukawa at large tan(beta)
 ! 25.12.01: - adding Sigma_Fermion
 ! 10.03.03: - adding call to Sigma_Fermion3
 !           - adding new routine Yukawas
 ! 26.03.04: - changing resummation in case of generation mixing
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_run
  Real(dp), Intent(inout) :: mSpm2(:), mP02(:)
  Real(dp), Intent(in) :: mC(:), mN(:), mSpm(:), mUsquark(:), mDsquark(:)     &
    & , mSlepton(:), mUsquark2(:), mDsquark2(:), mSlepton2(:), mSneutrino2(:) &
    & , mS02(:), RP0(:,:), mglu, RS0(:,:), vevSM(2), delta0
  Complex(dp), Intent(in) :: U(:,:), V(:,:), N(:,:), RSpm(:,:)           &
    & , RDsquark(:,:), RUsquark(:,:), RSlepton(:,:), RSneutrino(:,:)     &
    & , phase_glu
  Complex(dp), Intent(in) :: Unu_L(3,3)
  Complex(dp), Intent(inout), Dimension(3,3) :: uU_L, uU_R, uD_L, uD_R   &
    & , uL_L, uL_R
  Real(dp), Intent(out) :: g1(57), mW2_run, mZ2_run
  Integer, Intent(inout) :: kont

  Real(dp), Save :: vevs_DR(2)

  Integer :: n_S0, i1, i_loop, i_loop_max, i2
  Real(dp) :: mC2( Size(mC) ), mN2(Size(mN) ), D_mat(3,3)
  Real(dp) :: test, alphaMZ, alpha3, gSU2, rho, delta_rho, sinW2_DR, vev2    &
    & , mZ2_mZ, CosW2SinW2, gauge(3), delta, sinW2_old, delta_r     &
    & , p2, gSU3, tanb, xt2, fac(2), SigQCD, delta_rw, sinW2, cosW2 
  Real(dp), Dimension(3) :: mf_d_DR, mf_l_DR, mf_u_DR, mf_d2, mf_l2, mf_u2
  Complex(dp) :: dmZ2, dmW2, dmW2_0, yuk_tau, yuk_t, yuk_b, SigLep, Sigdown  &
    & , SigUp
  Complex(dp), Dimension(3,3) :: SigS_u, sigR_u, SigL_u, SigS_d, SigR_d    &
    & , SigL_d, SigS_l, sigR_l, SigL_l, Y_u, Y_d, Y_l, adCKM, uU_L_T, uU_R_T &
    & , uD_L_T, uD_R_T, uL_L_T, uL_R_T, Y_l_old, Y_d_old, Y_u_old
  Logical :: converge
  Real(dp), Parameter :: e_d=-1._dp/3._dp, e_u=2._dp/3._dp, e_e=-1._dp &
    & , T3_d=-0.5_dp, T3_u=0.5_dp, mf_nu(3) = (/0._dp, 0._dp, 0._dp /)
  Complex(dp), Parameter :: Y_nu(3,3) = ZeroC
  Complex(dp), Dimension(6,6) :: rot, RUsq_ckm, RDsq_ckm, RUsq_in, RDsq_in

  Iname = Iname + 1
  nameOfUnit(Iname) = "BoundaryEWold"
  !----------------------------------------
  ! checking if masses squared are positiv
  !----------------------------------------
  If (Min(Minval(mUsquark2), Minval(mDSquark2), Minval(mSlepton2)           &
     &    ,Minval(mSneutrino2), Minval(mS02), Minval(mP02), Minval(mSpm2))  &
     & .Lt. 0._dp ) Then
   kont = -401
   Call AddError(401)
   Iname = Iname - 1
   Return
  End If

  sinW2 = 1._dp - mW2/mZ2
  mC2 = mC**2
  mN2 = mN**2
  n_s0 = Size(mS02)
  !-------------------------------------------------------------------
  ! setting renormalisation scale to m_Z, because the RGEs start there
  !-------------------------------------------------------------------
  test = SetRenormalizationScale(mZ2)
  tanb = vevSM(2) / vevSM(1)
  !-------------------------------------------------------------------
  ! initialization of LoopMasses
  !-------------------------------------------------------------------
  Call  SetLoopMassModel(Size(mC), Size(mN), n_s0, n_s0, Size(mSpm) &
                       & , Size(mSlepton), Size(msneutrino2))
  !-----------------
  ! sin(theta_W)^2
  !-----------------
  If (i_run.Eq.1) Then
   sinW2_DR = sinW2
   sinW2_old = sinW2_DR
   Y_l = 0._dp
   Do i1=1,3
    y_l(i1,i1) = sqrt2 * mf_l_mZ(i1) / vevSM(1)
   End Do
   mf_l2 = mf_l_mZ**2
   mf_d2 = mf_d_mZ**2
   mf_u2 = mf_u_mZ**2
   uL_L_T = id3C
   uL_R_T = id3C
  Else
   sinW2_DR = sinW2_DR_mZ
   sinW2_old = sinW2_DR
   Y_l = Yl_mZ

   Call FermionMass(Yl_mZ,vevs_DR_save(1),mf_l2,uL_L_T,uL_R_T,kont)
   Call QuarkMasses_and_PhaseShifts(Yd_mZ, Yu_mZ, vevs_DR, mf_d2, uD_L_T, uD_R_T &
                                     & , mf_u2, uU_L_T, uU_R_T)

   mf_l2 = mf_l2**2
   mf_d2 = mf_d2**2
   mf_u2 = mf_u2**2
  End If
  !-----------
  ! alpha(mZ)
  !-----------
  alphaMZ = AlphaEwDR(mZ, mSpm, mUsquark, mDSquark, mSlepton, mC)
  !-----------
  ! alpha_s(mZ)
  !-----------
  alpha3 = AlphaSDR(mZ, mglu, mUSquark, mDSquark, Sqrt(mf_u2(3)) )
  gSU3 = Sqrt( 4._dp*pi*alpha3)
  !--------------------
  ! for 2-loop parts
  !--------------------
   xt2 = 3._dp * (G_F * mf_u2(3) * oo8pi2 * oosqrt2)**2        &
      & * Abs(RS0(1,2))**2 * rho_2(Sqrt(mS02(1))/mf_U(3))     &
      & * ((1._dp+tanb**2)/tanb**2) 
   fac(1) = alphaMZ * alphaS_mZ * oo4pi                                    &
        & * (2.145_dp * mf_u2(3)/mZ2 + 0.575 * Log(mf_u(3)/mZ) - 0.224_dp  &
        &   - 0.144_dp * mZ2 / mf_u2(3)) / Pi
   fac(2) = alphamZ * alphaS_mZ * oo4pi                                    &
       & * (-2.145_dp * mf_u2(3)/mW2 + 1.262 * Log(mf_u(3)/mZ) - 2.24_dp  &
       &   - 0.85_dp * mZ2 / mf_u2(3)) / Pi 

  Do i1=1,100
   !------------------------------
   ! for gauge invariance
   !------------------------------
   mP02(1) = mZ2
   mSpm2(1) = mW2
   gSU2 = Sqrt( 4._dp*pi*alphamZ/sinW2_DR)
   Call PiZZT1(mZ2, gSU2, sinW2_DR, vevSM, mZ2, mW2, mS02, RS0, mP02, RP0   &
   & , mSpm2, RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton, mUsquark2 &
   & , RUsquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2, mC, mC2, U, V    &
   & , mN, mN2, N, dmZ2)

   mZ2_mZ = Real(dmZ2+mZ2,dp)
   If (mZ2_mZ.Lt.0._dp) Then
    Iname = Iname - 1
    kont = -402
    Call AddError(402) 
    Return
   End If
   mZ2_run = mZ2_mZ
   mW2_run = mZ2_mZ * (1._dp - sinW2_DR)
   !------------------------------
   ! for gauge invariance
   !------------------------------
   mP02(1) = mZ2_run
   mSpm2(1) = mW2_run

   CosW2SinW2 = (1._dp - sinW2_DR) * sinW2_DR
   vev2 =  mZ2_mZ * CosW2SinW2 / (pi * alphamZ)
   vevs_DR(1) = Sqrt(vev2 / (1._dp+tanb**2) )
   vevs_DR(2) = tanb * vevs_DR(1)
   !---------------------------------------
   ! recalculation, using running masses
   !---------------------------------------
   Call PiZZT1(mZ2, gSU2, sinW2_DR, vevs_DR, mZ2_mZ, mW2_run, mS02, RS0       &
   & , mP02, RP0, mSpm2, RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton   &
   & , mUsquark2, RUsquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2, mC, mC2 &
   & , U, V, mN, mN2, N  &
   & , dmZ2)
   mZ2_mZ = Real(dmZ2+mZ2,dp)
   If (mZ2_mZ.Lt.0._dp) Then
    Iname = Iname - 1
    kont = -402
    Call AddError(402) 
    Return
   End If
   mZ2_run = mZ2_mZ
   mW2_run = mZ2_mZ * (1._dp - sinW2_DR)
   mP02(1) = mZ2_run
   mSpm2(1) = mW2_run

   CosW2SinW2 = (1._dp - sinW2_DR) * sinW2_DR
   vev2 =  mZ2_mZ * CosW2SinW2 / (pi * alphamZ)
   vevs_DR(1) = Sqrt(vev2 / (1._dp+tanb**2) )
   vevs_DR(2) = tanb * vevs_DR(1)

   Call PiWWT1(mW2, gSU2, sinW2_DR, mS02, RS0, mSpm2, RSpm, vevs_DR        &
         & , mP02, RP0, mSneutrino2, RSneutrino, mSlepton2, RSlepton       &
         & , mUSquark2, RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2 &
         & , CKM, mN, mN2, N, mC, mC2, U, V, mZ2_mZ, mW2_run, dmW2)

   Call PiWWT1(0._dp, gSU2, sinW2_DR, mS02, RS0, mSpm2, RSpm , vevs_DR     &
         & , mP02, RP0, mSneutrino2, RSneutrino, mSlepton2, RSlepton       &
         & , mUSquark2, RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2 &
         & , CKM, mN, mN2, N, mC, mC2, U, V, mZ2_mZ, mW2_run, dmW2_0)

   rho = (1._dp + Real(dmZ2,dp)/mZ2) / (1._dp + Real(dmW2,dp) / mW2)
   delta_rho = 1._dp - 1._dp / rho 

   CosW2SinW2 = (1._dp - sinW2_DR) * sinW2_DR
   Call delta_VB(gSU2, sinW2, sinW2_DR, rho, mC, mC2, U, V, mN, mN2, N &
        & , Y_l, UL_L_T, UL_R_T, mSlepton2, Rslepton, mSneutrino2      &
        & , RSneutrino, UNu_L, delta)
   delta_r = rho*Real(dmW2_0,dp)/mW2 - Real(dmZ2,dp) / mZ2 + delta
   rho = 1._dp /  (1._dp - delta_rho - fac(2) / sinW2_DR - xt2)
   delta_r = rho*Real(dmW2_0,dp)/mW2 - Real(dmZ2,dp) / mZ2 + delta   &
         & + fac(1) / CosW2SinW2 - xt2 * (1-delta_r) * rho
   CosW2SinW2 = pi * alphamZ / (sqrt2 * mZ2 * G_F * (1-delta_r) )
   sinW2_DR = 0.5_dp - Sqrt(0.25_dp - CosW2SinW2)
   If (sinW2_DR.Lt.0._dp) Then
    kont = -403
    Call AddError(403)
    Iname = Iname -1
    Return
   End If
   !-----------------------------------------------------------------
   ! have at least two runs before leaving the loop
   !-----------------------------------------------------------------
   If ((Abs(sinW2_DR-sinW2_old).Lt. 0.1_dp*delta0).And.(i1.Gt.1)) Exit
   sinW2_old = sinW2_DR
   delta_rw = delta_rho*(1._dp-delta_r) + delta_r
   If ((0.25_dp-alphamz*pi/(sqrt2*G_F*mz2*rho*(1._dp-delta_rw))).Lt.0._dp) Then
    If (Errorlevel.Ge.0) Then
     Write(Errcan,*) "Problem in subroutine "//NameofUnit(Iname)
     Write(Errcan,*) "In the calculation of mW", &
          & 0.25_dp-alphamz*pi/(sqrt2*G_F*mz2*rho*(1._dp-delta_rw))
     Write(errcan,*) "mf_l_mZ, vevSM"
     Write(errcan,*) mf_l_mZ,vevSM
     Write(errcan,*) "mZ2_mZ, dmW2, dmW2_0, rho, delta_rho, sinW2_DR"
     Write(errcan,*) mZ2_mZ, dmW2, dmW2_0, rho, delta_rho, sinW2_DR
     Write(ErrCan,*) "delta_r, cosW2SinW2, delta_rw"
     Write(ErrCan,*) delta_r, cosW2SinW2, delta_rw
    End If
    kont = -404
    Call AddError(404)
    Iname = Iname - 1
    Return
   End If
   mW2 = mZ2 * rho * ( 0.5_dp &
      &         +Sqrt(0.25_dp-alphamz*pi/(sqrt2*G_F*mz2*rho*(1._dp-delta_rw))))

   mSpm2(1) = mW2_run ! for this loop
   cosW2 = mW2 / mZ2
   sinW2 = 1._dp - cosW2
  End Do
  mSpm2(1) = mW2
  !--------------------------------------------------------------------------
  ! recalcuating m_W and sin^2_W, the formula for m_W is base on Eq.25 of
  ! G.Degrassi et al., NPB351, 49 (1991)  
  !--------------------------------------------------------------------------
  delta_rw = delta_rho*(1._dp-delta_r) + delta_r
  mW2 = mZ2 * rho * ( 0.5_dp &
      &         +Sqrt(0.25_dp-alphamz*pi/(sqrt2*G_F*mz2*rho*(1._dp-delta_rw))))
  mW = Sqrt(mW2)
  cosW2 = mW2 / mZ2
  sinW2 = 1._dp - cosW2
  !---------------------------
  ! gauge couplings and vevs
  !---------------------------
  gauge(1) = Sqrt( 4._dp*pi*alphamZ/(1._dp-sinW2_DR) )
  gauge(2) = Sqrt( 4._dp*pi*alphamZ/sinW2_DR)
  gauge(3) = Sqrt( 4._dp*pi*alpha3)
  vev2 =  mZ2_mZ * CosW2SinW2 / (pi * alphamZ)
  vevs_DR(1) = Sqrt(vev2 / (1._dp+tanb**2) )
  vevs_DR(2) = tanb * vevs_DR(1)

  !-------------------------------------
  ! Initialize fermion mixing matrices
  !-------------------------------------
  uU_L = id3C
  uU_R = id3C
  uD_L = id3C
  uD_R = id3C
  uL_L = id3C
  uL_R = id3C
  If (GenerationMixing) Then
   Call Adjungate(CKM, adCKM)
   If (YukawaScheme.Eq.1) Then
    uU_L = CKM
   Else
    uD_L = adCKM
   End If
  End If

  If (i_run.Eq.1) Then
   !--------------------------------------------------------------------------
   ! shifting light fermion masses to DR-scheme, only gluon and photon part
   ! except for m_t
   !--------------------------------------------------------------------------
   mf_l_DR_SM = &
    &     mf_l_mZ * (1._dp - oo8pi2 *3._dp *(gauge(1)**2-gauge(2)**2)/16._dp)
   mf_d_DR_SM = mf_d_mZ * (1._dp - alpha3 / (3._dp*pi)                  &
         &               - 23._dp * alpha3**2 / (72._dp * Pi2 )         &
         &               + oo8pi2 * 3._dp * gauge(2)**2 / 16._dp        &
         &               - oo8pi2 * 13._dp * gauge(1)**2 / 144._dp  )
   mf_u_DR_SM(1:2) = mf_u_mZ(1:2)  * (1._dp - alpha3 / (3._dp*pi)       &
         &               - 23._dp * alpha3**2 / (72._dp * Pi2 )         &
         &               + oo8pi2 * 3._dp * gauge(2)**2 / 16._dp        &
         &               - oo8pi2 * 7._dp * gauge(1)**2 / 144._dp  )
   mf_u_DR_SM(3) = mf_u(3) ! QCD + QED shift will be added later
   mf_l_DR = mf_l_DR_SM
   mf_d_DR = mf_d_DR_SM
   mf_u_DR = mf_u_DR_SM
   !---------------------------------------------------------------------
   ! Yukawa couplings
   !--------------------------------------------------------------------
   Y_d = 0._dp
   Y_u = 0._dp 
   Y_l = 0._dp
   Do i1=1,3
    Y_u(i1,i1) = sqrt2 * mf_u_DR_SM(i1) / vevs_DR(2)
    Y_l(i1,i1) = sqrt2 * mf_l_DR_SM(i1) / vevs_DR(1)
    Y_d(i1,i1) = sqrt2 * mf_d_DR_SM(i1) / vevs_DR(1)
   End Do
   If (GenerationMixing) Then
    If (YukawaScheme.Eq.1) Then
     Y_u = Matmul(Transpose(uU_L),Y_u) 
    Else
     Y_d = Matmul(Transpose(uD_L),Y_d) 
    End If
   End If
   !--------------------------------------------
   ! the starting point of the tree-level mixing
   !--------------------------------------------
   uU_L_T = uU_L
   uU_R_T = uU_R
   uD_L_T = uD_L
   uD_R_T = uD_R
   uL_L_T = uL_L
   uL_R_T = uL_R
  Else 
   !--------------------------------------------------------------------------
   ! take Yukawas from previous run
   !--------------------------------------------------------------------------
   Y_l = Yl_mZ
   Y_d = Yd_mZ
   Y_u = Yu_mZ
   Call FermionMass(Y_l,vevs_DR(1),mf_l_DR,uL_L_T,uL_R_T,kont)
   Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevs_DR, mf_d_DR, uD_L_T, uD_R_T &
                                     & , mf_u_DR, uU_L_T, uU_R_T)
  End If ! i_run.eq.1

  !---------------------------------------------
  ! shifting mixing matrices to superCKM basis 
  !---------------------------------------------
  rot = 0._dp
  rot(1:3,1:3) = Conjg(uU_L_T)
  rot(4:6,4:6) = uU_R_T
  RUsq_ckm = Matmul(RUSquark, Transpose(rot))

  rot = 0._dp
  rot(1:3,1:3) = Conjg(uD_L_T)
  rot(4:6,4:6) = uD_R_T
  RDsq_ckm = Matmul(RDSquark, Transpose(rot))

  converge = .False.

  Y_l_old = Y_l
  Y_d_old = Y_d
  Y_u_old = Y_u

  !------------------------------
  ! now the iteration
  !------------------------------
  If (FermionMassResummation) Then
   i_loop_max = 100 ! this should be sufficient
  Else
   i_loop_max = 1
  End If
  Do i_loop =1,i_loop_max
   yuk_b = Y_d(3,3)! for checking of convergence
   yuk_t = Y_u(3,3)
   yuk_tau = Y_l(3,3)

   If (GenerationMixing) Then
    !---------------------------------------------------------------
    ! rotate squarks from superCKM basis to new electroweak basis
    !---------------------------------------------------------------
    rot = 0._dp
    rot(1:3,1:3) = uU_L_T
    rot(4:6,4:6) = Conjg(uU_R_T)
    RUsq_in = Matmul(RUsq_ckm, rot)
    rot = 0._dp
    rot(1:3,1:3) = uD_L_T
    rot(4:6,4:6) = Conjg(uD_R_T)
    RDsq_in = Matmul(RDsq_ckm, rot)

    p2 = 0._dp ! for off-diagonal elements
    ! u-quarks
    Call Sigma_Fermion3(p2, mf_u_DR, Y_u, uU_L_T, uU_R_T, gSU2, gSU3, sinW2_DR &
        & , T3_u, e_u, mf_d_DR, Y_d, uD_L_T, uD_R_T, mUSquark2,RUsq_in         &
        & , mDSquark2, RDsq_in, mglu , phase_glu, mN, mN2, N, mC, mC2, U, V   &
        & , mS02, RS0, mP02, RP0, mSpm2 , RSpm, mZ2_run, mW2_run, .True.       &
        & , SigS_u, SigL_u, SigR_u, SigQCD)
    ! d-quarks
    Call Sigma_Fermion3(p2, mf_d_DR, Y_d, uD_L_T, uD_R_T, gSU2, gSU3, sinW2_DR &
        & , T3_d, e_d, mf_u_DR, Y_u, uU_L_T, uU_R_T, mDSquark2, RDsq_in        &
        & , mUSquark2, RUsq_in,  mglu , phase_glu, mN, mN2, N, mC, mC2, U, V   &
        & , mS02, RS0, mP02, RP0, mSpm2 , RSpm, mZ2_run, mW2_run , .True.      &
        & , SigS_d, SigL_d, SigR_d)
    ! leptons
    Call Sigma_Fermion3(p2, mf_l_DR, Y_l, uL_L, uL_R, gSU2, gSU3, sinW2_DR    &
        & , T3_d, e_e, mf_nu, Y_nu, id3C, id3C, mSlepton2, RSlepton           &
        & , mSneutrino2, RSneutrino, mglu , phase_glu, mN, mN2, N, mC, mC2, U &
        & , V, mS02, RS0, mP02, RP0, mSpm2 , RSpm, mZ2_run, mW2_run, .False.  &
        & , SigS_l, SigL_l, SigR_l)

    mf_u_DR_SM(3) = mf_u(3) + SigQCD

    Call Yukawas(mf_u_DR_SM, vevs_DR(2), uU_L, uU_R, SigS_u, SigL_u, SigR_u &
          & , Y_u, .False., kont)

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
    Call Yukawas(mf_d_DR_SM, vevs_DR(1), uD_L, uD_R, SigS_d, SigL_d, SigR_d &
          & , Y_d, FermionMassResummation, kont)

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
    Call Yukawas(mf_l_DR_SM, vevs_DR(1), id3C, id3C, SigS_l, SigL_l, SigR_l &
          & , Y_l, FermionMassResummation, kont)

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If

    !----------------------------------------------------------------
    ! I am only interested in the mixing matrices and, thus, it does
    ! not matter which vev I am using
    !----------------------------------------------------------------
    Call FermionMass(Y_l,vevs_DR(1),mf_l_DR,uL_L_T,uL_R_T,kont)
    Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevs_DR, mf_d_DR, uD_L_T, uD_R_T &
                                     & , mf_u_DR, uU_L_T, uU_R_T)

    converge = .True.
    D_mat = Abs(Abs(Y_l) - Abs(Y_l_old))
    Where (Abs(Y_l).Ne.0._dp) D_mat = D_mat / Abs(Y_l)
    Do i1=1,3
     If (D_mat(i1,i1).Gt.0.1_dp*delta0) converge = .False.
     Do i2=i1+1,3
      If (D_mat(i1,i2).Gt.delta0) converge = .False.
      If (D_mat(i2,i1).Gt.delta0) converge = .False.
     End Do
    End Do
    D_mat = Abs(Abs(Y_d) - Abs(Y_d_old))
    Where (Abs(Y_d).Ne.0._dp) D_mat = D_mat / Abs(Y_d)
    Do i1=1,3
     If (D_mat(i1,i1).Gt.0.1_dp*delta0) converge = .False.
     Do i2=i1+1,3
      If (D_mat(i1,i2).Gt.10._dp*delta0) converge = .False.
      If (D_mat(i2,i1).Gt.10._dp*delta0) converge = .False.
     End Do
    End Do
    D_mat = Abs(Abs(Y_u) - Abs(Y_u_old))
    Where (Abs(Y_u).Ne.0._dp) D_mat = D_mat / Abs(Y_u)
    Do i1=1,3
     If (D_mat(i1,i1).Gt.0.1_dp*delta0) converge = .False.
     Do i2=i1+1,3
      If (D_mat(i1,i2).Gt.10._dp*delta0) converge = .False.
      If (D_mat(i2,i1).Gt.10._dp*delta0) converge = .False.
     End Do
    End Do

    If (converge) Exit

    Y_l_old = Y_l
    Y_u_old = Y_u
    Y_d_old = Y_d

  Else ! .not.GenerationMixing
   Do i1=1,3

     p2 = mf_d_DR(i1)**2
     Call Sigma_Fermion(p2, i1, mf_d_DR, Y_d, id3C, id3C,gSU2,gSU3,sinW2_DR   &
      & ,T3_d, e_d, mf_u_DR, Y_u, id3C, id3C, mDSquark2, RDSquark, mUSquark2  &
      & ,RUSquark, mglu, phase_glu, mN, mN2, N, mC, mC2, U, V, mS02, RS0      &
      & ,mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run, .True., .True., SigDown)
     If (FermionMassResummation) Then
      mf_d_DR(i1) = mf_d_DR_SM(i1) / (1- Real(SigDown,dp) / mf_d_DR(i1) )
     Else
      mf_d_DR(i1) = mf_d_DR_SM(i1) + Real(SigDown,dp)
     End If
     Y_d(i1,i1) = sqrt2 * mf_d_DR(i1) / vevs_DR(1)

     p2 = mf_u_DR(i1)**2
     If (i1.Lt.3) Then
      Call Sigma_Fermion(p2, i1, mf_u_dR, Y_u, id3C, id3C, gSU2, gSU3,sinW2_DR&
        & ,T3_u, e_u, mf_d_DR, Y_d, id3C, id3C,mUSquark2,RUSquark,mDSquark2   &
        & ,RDSquark, mglu, phase_glu, mN, mN2, N, mC, mC2, U, V, mS02, RS0    &
        & ,mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run, .True., .True., SigUp)
     Else
      Call Sigma_Fermion(p2, i1, mf_u_DR, Y_u, id3C, id3C, gSU2, gSU3,sinW2_DR&
        & ,T3_u, e_u, mf_d_DR, Y_d, id3C, id3C,mUSquark2,RUSquark,mDSquark2   &
        & ,RDSquark, mglu, phase_glu, mN, mN2, N, mC, mC2, U, V, mS02, RS0    &
        & ,mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run, .True., .False., SigUp)
     End If
     mf_u_DR(i1) = mf_u_DR_SM(i1) + Real(SigUp,dp)
     Y_u(i1,i1) = sqrt2 * mf_u_DR(i1) / vevs_DR(2)

     p2 = mf_l_DR(i1)**2
     Call Sigma_Fermion(p2, i1, mf_l, Y_l, id3C, id3C, gSU2, gSU3, sinW2_DR  &
        & , T3_d, e_e, mf_nu, Y_nu, id3C, id3C, mSlepton2, RSlepton          &
        & , mSneutrino2, RSneutrino, mglu, phase_glu, mN, mN2, N, mC, mC2, U &
        & , V, mS02, RS0, mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run           &
        & , .False., .True., SigLep)
     If (FermionMassResummation) Then
      mf_l_DR(i1) = mf_l_DR_SM(i1) / (1- Real(SigLep,dp) / mf_l_DR(i1) )
     Else
      mf_l_DR(i1) = mf_l_DR_SM(i1) + Real(SigLep,dp)
     End If
     Y_l(i1,i1)  = sqrt2 * mf_l_DR(i1) / vevs_DR(1)
    
    End Do

    If (    (      Abs((yuk_tau-y_l(3,3))/y_l(3,3)).Lt. 0.1_dp*delta0) &
        &    .And.(Abs((yuk_t-y_u(3,3))  /y_u(3,3)).Lt. 0.1_dp*delta0) &
        &    .And.(Abs((yuk_b-y_d(3,3))  /y_d(3,3)).Lt. 0.1_dp*delta0) ) Then
     converge = .True.
     Exit
    End If
   End If  ! GenerationMixing

   !--------------------------------------------------
   ! Either we have run into a numerical problem or
   ! perturbation theory breaks down
   !--------------------------------------------------

   If (    (Abs(mf_l_DR(3)/mf_l_mZ(3)).Lt.0.1_dp)  &
     & .Or.(Abs(mf_l_DR(3)/mf_l_mZ(3)).Gt.10._dp) ) Then
    Iname = Iname - 1
    kont = -405
    Call AddError(405)
    Return
   Else If (    (Abs(mf_d_DR(3)/mf_d_mZ(3)).Lt.0.1_dp)  &
          & .Or.(Abs(mf_d_DR(3)/mf_d_mZ(3)).Gt.10._dp) ) Then
    Iname = Iname - 1
    kont = -406
    Call AddError(406)
    Return
   Else If (    (Abs(mf_u_DR(3)/mf_u_mZ(3)).Lt.0.1_dp)  &
          & .Or.(Abs(mf_u_DR(3)/mf_u_mZ(3)).Gt.10._dp) ) Then
    Iname = Iname - 1
    kont = -407
    Call AddError(407)
    Return
   End If

  End Do ! i_loop

  If ((.Not.converge).And.FermionMassResummation) Then
   Write (ErrCan,*) 'Problem in subroutine BoundaryEWold!!'
   Write (ErrCan,*) "After",i_loop-1,"iterations no convergence of Yukawas"
   Write (ErrCan,*) 'yuk_tau,yuk_l(3,3)',yuk_tau,y_l(3,3)
   Write (ErrCan,*) 'yuk_b,yuk_d(3,3)',yuk_b,y_d(3,3)
   Write (ErrCan,*) 'yuk_t,yuk_u(3,3)',yuk_t,y_u(3,3)
  End If
  !----------------------------------------------------------------
  ! the RGE paper defines the Yukawas transposed to my conventions
  !----------------------------------------------------------------

  Yl_mZ = Y_l
  Yd_mZ = Y_d
  Yu_mZ = Y_u
  vevs_DR_save = vevs_DR
  Y_u = Transpose(Y_u)
  Y_d = Transpose(Y_d)
  Y_l = Transpose(Y_l)
  sinW2_DR_mZ = sinW2_DR
  gauge(1) = Sqrt( 5._dp/3._dp) * gauge(1)
  gauge_mZ = gauge

  Call  CouplingsToG(gauge, y_l, y_d, y_u, g1)

  !----------------------------------------------
  ! resetting scale
  !----------------------------------------------
  test = SetRenormalizationScale(test)

  Iname = Iname - 1

 Contains

  Real(dp) Function rho_2(r)
  Implicit None
   Real(dp), Intent(in) :: r
   Real(dp) :: r2, r3
   r2 = r*r
   r3 = r2*r
   rho_2 = 19._dp - 16.5_dp * r + 43._dp * r2 / 12._dp             &
       & + 7._dp * r3 / 120._dp                                    &
       & - Pi * Sqrt(r) * (4._dp - 1.5_dp * r + 3._dp * r2/32._dp  &
       &                  + r3/256._dp)                            &
       & - Pi2 * (2._dp - 2._dp * r + 0.5_dp * r2)                 &
       & - Log(r) * (3._dp * r - 0.5_dp * r2) 
  End  Function rho_2

  Subroutine Yukawas(mf, vev, uL, uR, SigS, SigL, SigR, Y, ReSum, kont)
  !--------------------------------------------------------
  ! solves the matrix equation for Y by a transformation to
  ! a linear system of 9 equations in 9 unknowns
  ! written by Werner Porod, 19.03.03
  !--------------------------------------------------------
  Implicit None
   Integer, Intent(inout) :: kont
   Real(dp), Intent(in) :: mf(3), vev
   Complex(dp), Dimension(3,3), Intent(in) :: uL, uR, SigS, SigL, SigR
   Logical, Intent(in) :: ReSum
   Complex(dp), Intent(inout) :: Y(3,3)

   Integer :: i1
   Complex(dp), Dimension(3,3) :: mass, uLa, uRa, f, invf, invY

   !-------------------------------------
   ! first the mass matrix in DR scheme
   !-------------------------------------
   Call Adjungate(uL, uLa)
   Call Adjungate(uR, uRa)
   mass = ZeroC
   Do i1=1,3
    mass(i1,i1) = mf(i1)
   End Do
   mass = Matmul( Transpose(uL), Matmul(mass, uR) )
   !----------------------------------------
   ! setting up the equations
   !----------------------------------------
   Y = Y * vev * oosqrt2
   If (ReSum) Then
    kont = 0
    Call chop(Y)
    invY = Y
    Call gaussj(kont,invY,3,3)
    If (kont.Ne.0) Return

    f = id3C - Matmul(SigS,invY) - Transpose(SigL) - Matmul(Y,Matmul(SigR,invY))
    invf = f
    Call gaussj(kont,invf,3,3)
    If (kont.Ne.0) Return

    Y = Matmul(invf,mass)

   Else

    Y = mass + SigS + Matmul(Transpose(SigL),Y) + Matmul(Y,SigR)

   End If

   Y = sqrt2 * Y / vev

   Call chop(y)

  End Subroutine Yukawas

 End Subroutine BoundaryEWold

 Subroutine BoundaryEWnew(i_run, Q, tanb, Mi, T_l, T_d, T_u, M2_E, M2_L  &
    & , M2_D, M2_Q, M2_U, mu, mA, delta0, GenerationMixing, resum        &
    & , mZ2_run, mW2_run, g1, kont)
 !-----------------------------------------------------------------------
 ! Calculates gauge and yukawa couplings at m_Z in the DRbar scheme
 ! written by Werner Porod
 ! 26.03.04: new version to improve the resummation of chirally enhance terms
 !           is a combination of my previous routine and A.Crivellin et al
 !           arXiv:1103.4272
 !           All sfermion mass parameters are understood to be in the super-CKM
 !           basis which diagonalizes the tree-level Yukawa couplings
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_run
  Logical, Intent(in) :: GenerationMixing ! if generation mixing has to 
                                          ! be included
  Logical, Intent(in) :: resum ! if the resummation of chirally enhanced 
                               ! correction should be done
  Real(dp), Intent(in) :: Q    ! scale where calculation should be done
  Real(dp), Intent(in) :: tanb ! tan(beta)
  Complex(dp), Intent(in) :: Mi(3) ! gaugino mass parameters
  Complex(dp), Dimension(3,3), Intent(in) :: T_l, T_u, T_d & ! trilinears
    & , M2_E, M2_L, M2_D, M2_Q, M2_U                  ! soft masses squared 
  Complex(dp), Intent(in) :: mu  ! mu and the corresponding soft-parameter
  Real(dp), Intent(in) :: mA     ! mass of pseudoscalar boson
  Real(dp), Intent(in) :: delta0 ! relative precision which should be achived

  Real(dp), Intent(out) :: mZ2_run, mW2_run ! running vector boson masses
  Real(dp), Intent(out) :: g1(57)  ! vector containg gauge and Yukawa couplings
  Integer, Intent(inout) :: kont   ! for error message system
                                          

  Real(dp) :: mSpm2(2), mP02(2), mC(2), mN(4), mSpm(2), mSup(6), mSdown(6) &
    & , mSlept(6), mSup2(6), mSdown2(6), mSlept2(6), mSneut2(3), mS02(2)   &
    & , RP0(2,2), mglu, RS0(2,2), mC2(2), mN2(4), mat2R(2,2), vev2     &
    & , D_sneut, mSf(2), mSf2(2), cosb, sinb, test(2), mudim2, ML, MR

  !--------------------------
  ! save these for later runs
  !--------------------------
  Real(dp), Save :: gp, gSU2, vevSM(2), sinW2_DR
  Real(dp), Dimension(3), Save :: mf_d_DR, mf_l_DR, mf_u_DR, mf_d2, mf_l2, mf_u2

  Complex(dp), Dimension(3,3), Save :: Y_u, Y_d, Y_l, V0ckm


  Complex(dp) :: U(2,2), V(2,2), N(4,4), RSpm(2,2), RSdown(6,6), RSup(6,6) &
    & , RSlept(6,6), RSneut(3,3), phase_glu, mat3(3,3), A, Y       &
    & , Rsf(2,2), epsD(3), epsL(3), epsD_FC

  Complex(dp), Dimension(3,3) :: uU_L, uD_L, Del_ZuL(3,3), Del_ZdL(3,3)

  Real(dp), Save :: vevs_DR(2)

  Integer :: n_S0, i1, i_loop, i_loop_max, i2
  Real(dp) :: D_mat(3,3), mf(3), mf2(3)
  Real(dp) :: alphaMZ, alpha3, rho, delta_rho    &
    & , mZ2_mZ, CosW2SinW2, gauge(3), delta, sinW2_old, delta_r     &
    & , p2, gSU3, xt2, fac(2), SigQCD, delta_rw, sinW2, cosW2 
  Complex(dp) :: dmZ2, dmW2, dmW2_0, yuk_tau, yuk_t, yuk_b, SigLep, Sigdown  &
    & , SigUp, M2Q_ckm(3,3), V0ckmAd(3,3), Tu_ckm(3,3), TestC(3,3),ephi
  Complex(dp), Dimension(3,3) :: SigS_u, sigR_u, SigL_u, SigS_d, SigR_d    &
    & , SigL_d, SigS_l, sigR_l, SigL_l, adCKM, Y_l_old, Y_d_old, Y_u_old
  Complex(dp), Save :: Y0_d(3,3), Y0_u(3,3), Y0_l(3,3),CKM0(3,3)
  Real(dp) :: msu(6), msu2(6), msd(6), msd2(6), s13,c13,s12,s23,ar,ai
  Logical :: converge
  Real(dp), Parameter :: e_d=-1._dp/3._dp, e_u=2._dp/3._dp, e_e=-1._dp &
    & , T3_d=-0.5_dp, T3_u=0.5_dp, mf_nu(3) = (/0._dp, 0._dp, 0._dp /) &
    & , YL_q = 1._dp/3._dp, YR_d = 2._dp/3._dp, YR_u = -4._dp/3._dp    &
    & , YL_l = -1._dp, YR_l = 2._dp
  Complex(dp), Parameter :: Y_nu(3,3) = ZeroC, Unu_L(3,3) = id3C
  Complex(dp), Dimension(6,6) :: rot

  Iname = Iname + 1
  nameOfUnit(Iname) = "BoundaryEWnew"

  sinW2 = 1._dp - mW2/mZ2  ! on-shell weak mixing angle
  kont = 0 ! no problem so far :-)

  If (i_run.Eq.1) Then ! starting values
   !-------------------------------------------------------------------
   ! initialization of LoopMasses
   !-------------------------------------------------------------------
   Call  SetLoopMassModel(2, 4, 2, 2, 2, 6, 3) ! MSSM particle content

   vev2 = 1._dp / Sqrt(sqrt2*G_F)         ! this actually v and not v^2
   vevSM(1) = vev2 / Sqrt(1._dp + tanb**2)
   vevSM(2) = vevsM(1) * tanb

   gp = Sqrt(alpha_mZ * 4._dp * Pi / (1._dp-sinW2) )
   gSU2 = Sqrt(alpha_mZ * 4._dp * Pi / sinW2 )

   sinW2_DR = sinW2
   sinW2_old = sinW2_DR
   Y_l = 0._dp
   Y_d = 0._dp
   Y_u = 0._dp
   Do i1=1,3
    y_l(i1,i1) = sqrt2 * mf_l_mZ(i1) / vevSM(1)
    y_d(i1,i1) = sqrt2 * mf_d_mZ(i1) / vevSM(1)
    y_u(i1,i1) = sqrt2 * mf_u_mZ(i1) / vevSM(2)
   End Do
   mf_l2 = mf_l_mZ**2
   mf_d2 = mf_d_mZ**2
   mf_u2 = mf_u_mZ**2
   V0ckm = CKM
  Else ! take information from previous run
   sinW2_old = sinW2_DR
  End If

  V0ckmAd = Conjg(Transpose(V0ckm))

  If (GenerationMixing) Then
   M2Q_ckm = Matmul(Matmul(V0ckm,M2_Q),V0ckmAd)
   Call Chop(M2Q_ckm)
   Tu_ckm = Matmul(Conjg(V0ckm), T_u)
   Call Chop(Tu_ckm)
  Else
   M2Q_ckm = M2_Q
   Tu_ckm = T_u
  End If
  !-------------------------
  ! SUSY masses and mixings
  !-------------------------------
  ! gluino
  !-------------------------------
  mGlu = Abs( Mi(3) )
  Phase_Glu = Exp( (0._dp,1._dp) * Arg(Mi(3)) )
  !-------------------------------
  ! charginos + neutralinos
  !-------------------------------
  Call CharginoMass(Mi(2), mu, vevSM, gSU2, mC, U, V, kont)
  If (kont.Ne.0) Then
   Iname = Iname - 1
   Return
  End If
  mC2 = mC**2
  Call NeutralinoMass(Mi(1), Mi(2), mu, vevSM, gp, gSU2, mN, N, kont)
  If (kont.Ne.0) Then
   Iname = Iname - 1
   Return
  End If
  mN2 = mN**2
  !-----------
  ! sfermions
  !-----------
  If (GenerationMixing) Then
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (gSU2**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    mat3 =  M2_L + D_sneut * id3C
    Call ComplexEigenSystem(mat3,mSneut2,Rsneut,kont, test)

    If (kont.Eq.-14) Then
      Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
      Write(ErrCan,*) "test =",test
      Write(ErrCan,*) " "
      If (ErrorLevel.Eq.2) Call TerminateProgram
      kont = 0
    End If
  
    If (kont.Ne.0) Then
     Write(ErrCan,*) 'Problems with the diagonalization of sneutrinos'
     Write(ErrCan,*) 'in routine ',NameOfUnit(Iname),'. kont = ',kont
     Iname = Iname - 1
     Return
    Else
     Do i1=1,3
      If (mSneut2(i1).Ge.0._dp) Then
       mSneut(i1) = Sqrt( mSneut2(i1) )
      Else
       kont = -229
       Call AddError(229)
       If (ErrorLevel.Ge.0) Then
        Write(ErrCan,*) 'Warning from ',NameOfUnit(Iname),' mSneut2 ',i1
        Write(ErrCan,*) '< 0, : ',mSneut2(i1),'is set to 0 '
       End If
       mSneut2(i1) = 0
       mSneut(i1) = 0
       If (ErrorLevel.Eq.2) Call TerminateProgram
      End If
     Enddo
    End If
    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   !--------------
   ! Sleptons
   !--------------
    Call SfermionMass(M2_L, M2_E, T_l, mu, vevSM, Y_l, T3_d, YL_l, YR_l &
                    & , gSU2, gp, kont, mSlept, mSlept2, Rslept) 
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   !--------------
   ! D-squarks
   !--------------
   Call SfermionMass(M2_Q, M2_D, T_d, mu, vevSM, Y_d, T3_d, YL_q, YR_d &
                    & , gSU2, gp, kont, mSdown, mSdown2, Rsdown) 
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   !--------------
   ! U-squarks
   !--------------
   Call SfermionMass(M2Q_ckm, M2_U, Tu_ckm, mu, vevSM, Y_u, T3_u, YL_q, YR_u &
                    & , gSU2, gp, kont, mSup, mSup2, Rsup) 
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
  Else
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (gSU2**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    Do i1=1,3
     mSneut2(i1) = Real(M2_L(i1,i1),dp)  + D_sneut
     If (mSneut2(i1).Lt.0._dp) Then
      kont = -229
      Call AddError(229)
      If (ErrorLevel.Ge.0) Then
        Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
        Write(ErrCan,*) 'mSneutrino^2 ',i1,' <= 0 :',mSneut2(i1)
        Write(ErrCan,*) 'setting it to 10.'
        mSneut2(i1) = 10._dp
      End If
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
     mSneut(i1) = Sqrt(mSneut2(i1))
    End Do
    RSneut = id3C
    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   !----------------
   ! sleptons
   !---------------
   RSlepton = 0
   Do i1 = 1,3
    Ml = M2_L(i1,i1)
    Mr = M2_E(i1,i1)
    A = T_l(i1,i1)
    Y = Y_l(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3_d, YL_l, YR_l, gSU2, gp &
                    &, kont, msf, msf2, Rsf)
    mSlept(2*(i1-1)+1:2*i1) = msf
    mSlept2(2*(i1-1)+1:2*i1) = msf2
    RSlept(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do
   !---------------
   ! down-squarks
   !---------------
   RSdown = 0
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_D(i1,i1)
    A = T_d(i1,i1)
    Y = Y_d(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3_d, YL_q, YR_d, gSU2, gp &
                    &, kont, msf, msf2, Rsf)
    mSdown(2*(i1-1)+1:2*i1) = msf
    mSdown2(2*(i1-1)+1:2*i1) = msf2
    RSdown(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do
   !---------------
   ! up-squarks
   !---------------
   RSup = 0
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_U(i1,i1)
    A = T_u(i1,i1)
    Y = Y_u(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3_u, YL_q, YR_u, gSU2, gp &
                    &, kont, msf, msf2, Rsf)
    mSup(2*(i1-1)+1:2*i1) = msf
    mSup2(2*(i1-1)+1:2*i1) = msf2
    RSup(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do
  End If
  !---------------
  ! Higgs bosons
  !---------------
  ! pseudoscalar
  !---------------
   cosb = Sqrt(1._dp / (1._dp + tanb**2))
   sinb = cosb * tanb
   mP02(1) = mZ2
   mP02(2) = mA**2
   RP0(1,1) = - cosb
   RP0(1,2) = sinb 
   RP0(2,1) = sinb 
   RP0(2,2) = cosb
  !---------------
  ! charged
  !---------------
   mSpm2(1) = mW2
   mSpm2(2) = mSpm2(1) + mP02(2)
   mSpm = Sqrt(mSpm2)
   RSpm = RP0
  !---------------
  ! scalar
  !---------------
  mat2R(1,1) = mP02(1) * cosb**2 + mP02(2) * sinb**2
  mat2R(1,2) = -(mP02(1) + mP02(2)) * cosb * sinb
  mat2R(2,1) = mat2R(1,2)
  mat2R(2,2) = mP02(2) * cosb**2 + mP02(1) * sinb**2

  Call EigenSystem(mat2R,mS02,RS0,kont, test)
  !-------------------------------------------------------------------
  ! setting renormalisation scale to Q, because the RGEs start there
  !-------------------------------------------------------------------
  mudim2 = SetRenormalizationScale(Q**2)
  !-----------
  ! alpha(mZ)
  !-----------
  alphaMZ = AlphaEwDR(mZ, mSpm, mSup, mSdown, mSlept, mC)
  !-----------
  ! alpha_s(mZ)
  !-----------
  alpha3 = AlphaSDR(mZ, mglu, mSup, mSdown, Sqrt(mf_u2(3)) )
  gSU3 = Sqrt( 4._dp*pi*alpha3)
  !-----------------
  ! sin(theta_W)^2
  !-----------------
  !--------------------
  ! for 2-loop parts
  !--------------------
  xt2 = 3._dp * (G_F * mf_u2(3) * oo8pi2 * oosqrt2)**2        &
      & * Abs(RS0(1,2))**2 * rho_2(Sqrt(mS02(1))/mf_U(3))     &
      & * ((1._dp+tanb**2)/tanb**2) 
  fac(1) = alphaMZ * alphaS_mZ * oo4pi                                    &
        & * (2.145_dp * mf_u2(3)/mZ2 + 0.575 * Log(mf_u(3)/mZ) - 0.224_dp &
        &   - 0.144_dp * mZ2 / mf_u2(3)) / Pi
  fac(2) = alphamZ * alphaS_mZ * oo4pi                                    &
       & * (-2.145_dp * mf_u2(3)/mW2 + 1.262 * Log(mf_u(3)/mZ) - 2.24_dp  &
       &   - 0.85_dp * mZ2 / mf_u2(3)) / Pi 

  Call PiZZT1(mZ2, gSU2, sinW2_DR, vevSM, mZ2, mW2, mS02, RS0, mP02, RP0    &
    & , mSpm2, RSpm, mSneut2, RSneut, mSlept2, RSlept, mSup2, RSup, mSdown2 &
    & , RSdown, mf_l2, mf_u2, mf_d2, mC, mC2, U, V, mN, mN2, N, dmZ2)
  mZ2_mZ = Real(dmZ2+mZ2,dp)

  If (mZ2_mZ.Lt.0._dp) Then
   Iname = Iname - 1
   kont = -402
   Call AddError(402) 
   Return
  End If
  mZ2_run = mZ2_mZ
  mW2_run = mZ2_mZ * (1._dp - sinW2_DR)

  Do i1=1,100
   !------------------------------
   ! for gauge invariance
   !------------------------------
   ! pseudoscalar
   !---------------
   mP02(1) = mZ2_run
   !---------------
   ! charged
   !---------------
   mSpm2(1) = mW2_run
   mSpm2(2) = mSpm2(1) + mP02(2)
   !---------------
   ! scalar
   !---------------
   mat2R(1,1) = mP02(1) * cosb**2 + mP02(2) * sinb**2
   mat2R(1,2) = -(mP02(1) + mP02(2)) * cosb * sinb
   mat2R(2,1) = mat2R(1,2)
   mat2R(2,2) = mP02(2) * cosb**2 + mP02(1) * sinb**2

   Call EigenSystem(mat2R,mS02,RS0,kont, test)

   gSU2 = Sqrt( 4._dp*pi*alphamZ/sinW2_DR)

   CosW2SinW2 = (1._dp - sinW2_DR) * sinW2_DR
   vev2 =  mZ2_mZ * CosW2SinW2 / (pi * alphamZ)
   vevSM(1) = Sqrt(vev2 / (1._dp+tanb**2) )
   vevSM(2) = tanb * vevSM(1)
   !---------------------------------------
   ! recalculation, using running masses
   !---------------------------------------
   Call PiZZT1(mZ2, gSU2, sinW2_DR, vevSM, mZ2_mZ, mW2_run, mS02, RS0 &
   & , mP02, RP0, mSpm2, RSpm, mSneut2, RSneut, mSlept2, RSlept       &
   & , mSup2, RSup, mSdown2, RSdown, mf_l2, mf_u2, mf_d2, mC, mC2     &
   & , U, V, mN, mN2, N, dmZ2)
   mZ2_mZ = Real(dmZ2+mZ2,dp)

   If (mZ2_mZ.Lt.0._dp) Then
    Iname = Iname - 1
    kont = -402
    Call AddError(402) 
    Return
   End If
   mZ2_run = mZ2_mZ
   mW2_run = mZ2_mZ * (1._dp - sinW2_DR)
   mP02(1) = mZ2_run
   mSpm2(1) = mW2_run

   CosW2SinW2 = (1._dp - sinW2_DR) * sinW2_DR
   vev2 =  mZ2_mZ * CosW2SinW2 / (pi * alphamZ)
   vevSM(1) = Sqrt(vev2 / (1._dp+tanb**2) )
   vevSM(2) = tanb * vevSM(1)

   Call PiWWT1(mW2, gSU2, sinW2_DR, mS02, RS0, mSpm2, RSpm, vevSM  &
         & , mP02, RP0, mSneut2, RSneut, mSlept2, RSlept           &
         & , mSup2, RSup, mSdown2, RSdown, mf_l2, mf_u2, mf_d2     &
         & , CKM, mN, mN2, N, mC, mC2, U, V, mZ2_mZ, mW2_run, dmW2)

   Call PiWWT1(0._dp, gSU2, sinW2_DR, mS02, RS0, mSpm2, RSpm , vevSM &
         & , mP02, RP0, mSneut2, RSneut, mSlept2, RSlept             &
         & , mSup2, RSup, mSdown2, RSdown, mf_l2, mf_u2, mf_d2       &
         & , CKM, mN, mN2, N, mC, mC2, U, V, mZ2_mZ, mW2_run, dmW2_0)

   rho = (1._dp + Real(dmZ2,dp)/mZ2) / (1._dp + Real(dmW2,dp) / mW2)
   delta_rho = 1._dp - 1._dp / rho 

   CosW2SinW2 = (1._dp - sinW2_DR) * sinW2_DR
   Call delta_VB(gSU2, sinW2, sinW2_DR, rho, mC, mC2, U, V, mN, mN2, N     &
        & , Y_l, id3C, id3C, mSlept2, RSlept, mSneut2, RSneut, UNu_L, delta)

   delta_r = rho*Real(dmW2_0,dp)/mW2 - Real(dmZ2,dp) / mZ2 + delta
   rho = 1._dp /  (1._dp - delta_rho - fac(2) / sinW2_DR - xt2)
   delta_r = rho*Real(dmW2_0,dp)/mW2 - Real(dmZ2,dp) / mZ2 + delta   &
         & + fac(1) / CosW2SinW2 - xt2 * (1-delta_r) * rho
   CosW2SinW2 = pi * alphamZ / (sqrt2 * mZ2 * G_F * (1-delta_r) )
   sinW2_DR = 0.5_dp - Sqrt(0.25_dp - CosW2SinW2)
   If (sinW2_DR.Lt.0._dp) Then
    kont = -403
    Call AddError(403)
    Iname = Iname -1
    Return
   End If
   !-----------------------------------------------------------------
   ! have at least two runs before leaving the loop
   !-----------------------------------------------------------------
   If ((Abs(sinW2_DR-sinW2_old).Lt. 0.1_dp*delta0).And.(i1.Gt.1)) Exit

   sinW2_old = sinW2_DR
   delta_rw = delta_rho*(1._dp-delta_r) + delta_r
   If ((0.25_dp-alphamz*pi/(sqrt2*G_F*mz2*rho*(1._dp-delta_rw))).Lt.0._dp) Then
    If (Errorlevel.Ge.0) Then
     Write(Errcan,*) "Problem in subroutine "//NameofUnit(Iname)
     Write(Errcan,*) "In the calculation of mW", &
          & 0.25_dp-alphamz*pi/(sqrt2*G_F*mz2*rho*(1._dp-delta_rw))
     Write(errcan,*) "mf_l_mZ, vevSM"
     Write(errcan,*) mf_l_mZ,vevSM
     Write(errcan,*) "mZ2_mZ, dmW2, dmW2_0, rho, delta_rho, sinW2_DR"
     Write(errcan,*) mZ2_mZ, dmW2, dmW2_0, rho, delta_rho, sinW2_DR
     Write(ErrCan,*) "delta_r, cosW2SinW2, delta_rw"
     Write(ErrCan,*) delta_r, cosW2SinW2, delta_rw
    End If
    kont = -404
    Call AddError(404)
    Iname = Iname - 1
    Return
   End If
   mW2 = mZ2 * rho * ( 0.5_dp &
      &     + Sqrt(0.25_dp-alphamz*pi/(sqrt2*G_F*mz2*rho*(1._dp-delta_rw))))

   cosW2 = mW2 / mZ2
   sinW2 = 1._dp - cosW2
  End Do
  mSpm2(1) = mW2
  !--------------------------------------------------------------------------
  ! recalcuating m_W and sin^2_W, the formula for m_W is base on Eq.25 of
  ! G.Degrassi et al., NPB351, 49 (1991)  
  !--------------------------------------------------------------------------
  delta_rw = delta_rho*(1._dp-delta_r) + delta_r
  mW2 = mZ2 * rho * ( 0.5_dp &
      &         +Sqrt(0.25_dp-alphamz*pi/(sqrt2*G_F*mz2*rho*(1._dp-delta_rw))))
  mW = Sqrt(mW2)
  cosW2 = mW2 / mZ2
  sinW2 = 1._dp - cosW2
  !---------------------------
  ! gauge couplings and vevs
  !---------------------------
  gauge(1) = Sqrt( 4._dp*pi*alphamZ/(1._dp-sinW2_DR) )
  gauge(2) = Sqrt( 4._dp*pi*alphamZ/sinW2_DR)
  gauge(3) = gSU3
  gp = gauge(1)
  gSU2 = gauge(2)

  vev2 =  mZ2_mZ * CosW2SinW2 / (pi * alphamZ)
  vevSM(1) = Sqrt(vev2 / (1._dp+tanb**2) )
  vevSM(2) = tanb * vevSM(1)
  !-------------------------------------
  ! Initialize fermion mixing matrices
  !-------------------------------------
  uU_L = id3C
  uD_L = id3C
  If (GenerationMixing) Then
   Call Adjungate(CKM, adCKM)
   If (YukawaScheme.Eq.1) Then
    uU_L = CKM
   Else
    uD_L = adCKM
   End If
  End If

  If (i_run.Eq.1) Then
   !--------------------------------------------------------------------------
   ! shifting light fermion masses to DR-scheme, only gluon and photon part
   ! except for m_t
   !--------------------------------------------------------------------------
   mf_l_DR_SM = &
    &     mf_l_mZ * (1._dp - oo8pi2 *3._dp *(gauge(1)**2-gauge(2)**2)/16._dp)
   mf_d_DR_SM = mf_d_mZ * (1._dp - alpha3 / (3._dp*pi)                  &
         &               - 23._dp * alpha3**2 / (72._dp * Pi2 )         &
         &               + oo8pi2 * 3._dp * gauge(2)**2 / 16._dp        &
         &               - oo8pi2 * 13._dp * gauge(1)**2 / 144._dp  )
   mf_u_DR_SM(1:2) = mf_u_mZ(1:2)  * (1._dp - alpha3 / (3._dp*pi)       &
         &               - 23._dp * alpha3**2 / (72._dp * Pi2 )         &
         &               + oo8pi2 * 3._dp * gauge(2)**2 / 16._dp        &
         &               - oo8pi2 * 7._dp * gauge(1)**2 / 144._dp  )
   mf_u_DR_SM(3) = mf_u(3) ! QCD + QED shift will be added later

   mf_l_DR = mf_l_DR_SM
   mf_d_DR = mf_d_DR_SM
   mf_u_DR = mf_u_DR_SM
   !---------------------------------------------------------------------
   ! Yukawa couplings
   !--------------------------------------------------------------------
   Y_d = 0._dp
   Y_u = 0._dp 
   Y_l = 0._dp
   Do i1=1,3
    Y_u(i1,i1) = sqrt2 * mf_u_DR_SM(i1) / vevSM(2)
    Y_l(i1,i1) = sqrt2 * mf_l_DR_SM(i1) / vevSM(1)
    Y_d(i1,i1) = sqrt2 * mf_d_DR_SM(i1) / vevSM(1)
   End Do

  End If ! i_run.eq.1

  converge = .False.

  Y_l_old = Y_l
  Y_d_old = Y_d
  Y_u_old = Y_u
  !------------------------------
  ! now the iteration
  !------------------------------
  i_loop_max = 100 ! this should be sufficient
  Do i_loop =1,i_loop_max
   yuk_b = Y_d(3,3)! for checking of convergence
   yuk_t = Y_u(3,3)
   yuk_tau = Y_l(3,3)

   If (GenerationMixing) Then
    M2Q_ckm = Matmul(Matmul(V0ckm,M2_Q),V0ckmAd)
    Call Chop(M2Q_ckm)
    Tu_ckm = Matmul(Conjg(V0ckm), T_u)
    Call Chop(Tu_ckm)
    !--------------------------------------------------------------
    ! recalculate sfermion masses and mixing in the super-CKm basis
    !--------------------------------------------------------------
    D_sneut = 0.125_dp * (gSU2**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    mat3 =  M2_L + D_sneut * id3C
    Call ComplexEigenSystem(mat3,mSneut2,Rsneut,kont, test)

    Call SfermionMass3mssm(M2_L, M2_E, T_l, mu, vevSM, Y_l, T3_d, YL_l, YR_l &
                      &  ,gSU2, gp, kont, mSlept, mSlept2, RSlept)

    Call SfermionMass3mssm(M2_Q, M2_D, T_d, mu, vevSM, Y_d, T3_d, YL_q, YR_d &
                      &  ,gSU2, gp, kont, mSdown, mSdown2, RSdown)

    Call SfermionMass3mssm(M2Q_ckm, M2_U, Tu_ckm, mu, vevSM, Y_u, T3_u, YL_q &
                     &  , YR_u,gSU2, gp, kont, mSup, mSup2, RSup)

    p2 = 0._dp ! for off-diagonal elements
    ! u-quarks
    Call Sigma_Fermion3sckm(p2, gSU2, gSU3, sinW2_DR, T3_u, e_u, V0CKM, mf_u_DR &
        & , Y_u, mf_d_DR, Y_d, mSup2, RSup, mSdown2, RSdown, mglu, phase_glu  &
        & , mN, mN2, N, mC, mC2, U, V, mS02, RS0, mP02, RP0, mSpm2, RSpm      &
        & , mZ2_run, mW2_run, .True., SigS_u, SigL_u, SigR_u, SigQCD)
    mf_u_DR_SM(3) = mf_u(3) + SigQCD

    ! d-quarks, total self-energy
    Call Sigma_Fermion3sckm(p2, gSU2, gSU3, sinW2_DR, T3_d, e_d, V0CKM, mf_d_DR &
        & , Y_d, mf_u_DR, Y_u, mSdown2, RSdown, mSup2, RSup,  mglu, phase_glu &
        & , mN, mN2, N, mC, mC2, U, V, mS02, RS0, mP02, RP0, mSpm2, RSpm      &
        & , mZ2_run, mW2_run , .True., SigS_d, SigL_d, SigR_d)

    ! leptons
    Call Sigma_Fermion3sckm(p2, gSU2, gSU3, sinW2_DR, T3_d, e_e, Unu_L        &
        & , mf_l_DR, Y_l, mf_nu, Y_nu, mSlept2, RSlept, mSneut2, RSneut       &
        & , mglu, phase_glu, mN, mN2, N, mC, mC2, U, V, mS02, RS0 , mP02, RP0 &
        & , mSpm2 , RSpm, mZ2_run, mW2_run, .False., SigS_l, SigL_l, SigR_l)
    !-------------------------------------------------------
    ! now the chirally enhanced terms
    !-------------------------------------------------------
    epsD = 0._dp
    epsL = 0._dp
    If (Resum) Then
     Call Sigma_SM_chirally_enhanced(gauge, vevSM, mu, V0ckm, Y_l, Y_d, Y_u &
       & , Mi, T_l, T_d, T_u, M2_E, M2_L, M2_D, M2_Q, M2_U &
       & , epsD, epsL, epsD_FC, kont )
     If (kont.ne.0) then
      Iname = Iname - 1
      return
     End If 
    End If

    !---------------------------------------------------------
    !  Yukawa couplings
    !---------------------------------------------------------
    Call Yukawas(mf_u_DR_SM, vevSM(2), tanb, SigS_u, SigL_u, SigR_u, Y_u)
    If (Resum) Then
     Call Yukawas(mf_d_DR_SM, vevSM(1), tanb, SigS_d, SigL_d, SigR_d, Y_d, epsD)
     Call Yukawas(mf_l_DR_SM, vevSM(1), tanb, SigS_l, SigL_l, SigR_l, Y_l, epsL)
    Else
     Call Yukawas(mf_d_DR_SM, vevSM(1), tanb, SigS_d, SigL_d, SigR_d, Y_d)
     Call Yukawas(mf_l_DR_SM, vevSM(1), tanb, SigS_l, SigL_l, SigR_l, Y_l)
    End If

    !---------------------------------------------------------
    ! tree-level CKM (PMNS will be included at a later stage)
    !---------------------------------------------------------
    Del_ZuL = 0._dp
    mf = mf_u_DR_SM
    mf2 = mf**2
    Do i1=1,3
     Del_ZuL(i1,i1) = 1._dp
     Do i2=1,3
      If (i1.Ne.i2) Then
       Del_ZuL(i1,i2) = 2._dp * ( mf(i2)**2 * SigL_u(i1,i2)       &
                      &         + mf(i1)*mf(i2) * SigR_u(i1,i2)   &
                      &         + mf(i1) * SigS_u(i1,i2)          &
                      &         + mf(i2) * Conjg(SigR_u(i2,i1)) ) &
                      &       / (mf2(i1) - mf2(i2) )
      End If
     End Do
    End Do
    Del_ZdL = 0._dp
    mf = mf_d_DR_SM
    mf2 = mf**2
    Do i1=1,3
     Del_ZdL(i1,i1) = 1._dp
     Do i2=1,3
      If (i1.Ne.i2) Then
       Del_ZdL(i1,i2) = 2._dp * ( mf(i2)**2 * SigL_d(i1,i2)       &
                      &         + mf(i1)*mf(i2) * SigR_d(i1,i2)   &
                      &         + mf(i1) * SigS_d(i1,i2)          &
                      &         + mf(i2) * Conjg(SigR_d(i2,i1)) ) &
                      &       / (mf2(i1) - mf2(i2) )
      End If
     End Do
    End Do

    ckm0 = Matmul( Del_ZuL, CKM)
    ckm0 = Matmul( CKM0, Transpose(Conjg(Del_ZdL)))
    !---------------------------------------------
    ! the enhanced terms
    !---------------------------------------------
    Del_ZdL = id3c
    Del_ZdL(1,3) = Conjg(CKM0(3,1)) * CKM0(3,3) * epsD_FC
    Del_ZdL(2,3) = Conjg(CKM0(3,2)) * CKM0(3,3) * epsD_FC
    Del_ZdL(1,3) = - CKM0(3,1) * Conjg(CKM0(3,3) * epsD_FC)
    Del_ZdL(2,3) = - CKM0(3,2) * Conjg(CKM0(3,3) * epsD_FC)
    !------------------------------------------------------
    ! Trick to get the tree-level CKM unitary
    ! is needed as otherwise there are numerical problems
    !------------------------------------------------------
    TestC = 0._dp
    Do i1=1,3
     TestC(i1,i1) = Abs(Y_u(i1,i1))**2
    end do
    TestC = Matmul(Transpose(CKM0),TestC)
    TestC = Matmul(TestC, Conjg(CKM0) )
    Call Eigensystem(TestC, mf2, CKM0, kont, test)
    CKM0 = Conjg(CKM0)
    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
    !--------------------------------------------
    ! use phase freedom
    !--------------------------------------------
    CKM0(1,:) = CKM0(1,:) * Conjg(CKM0(1,1))/Abs(CKM0(1,1))
    CKM0(:,2) = CKM0(:,2) * Conjg(CKM0(1,2))/Abs(CKM(1,2))
    s13 = Abs(CKM0(1,3))   
    c13 = sqrt(1._dp - s13**2)
    s23 = Abs(CKM0(2,3))/c13
    s12 = Abs(CKM0(1,2))/c13
    aR = Real(CKM0(2,2),dp) + s12 * s23 * Real(CKM0(1,3),dp)
    aI =  s12 * s23 * Aimag(CKM0(1,3)) - Aimag(CKM0(2,2))
    Ephi = Cmplx(aR/Sqrt(aR**2+aI**2),aI/Sqrt(aR**2+aI**2),dp)

    CKM0(2:3,:) = Ephi * CKM0(2:3,:)
    CKM0(2,:) = CKM0(2,:) * Conjg(CKM0(2,3))/Abs(CKM0(2,3))
    CKM0(3,:) = CKM0(3,:) * Conjg(CKM0(3,3))/Abs(CKM0(3,3))

    !---------------------
    ! check convergence
    !---------------------
    converge = .True.
    Call Chop(Y_l)
    D_mat = Abs(Abs(Y_l) - Abs(Y_l_old))
    Where (Abs(Y_l).Ne.0._dp) D_mat = D_mat / Abs(Y_l)

    Do i1=1,3
     If (D_mat(i1,i1).Gt.delta0) converge = .False.
     Do i2=i1+1,3
      If (D_mat(i1,i2).Gt.delta0) converge = .False.
      If (D_mat(i2,i1).Gt.delta0) converge = .False.
     End Do
    End Do
    Call Chop(Y_u)
    D_mat = Abs(Abs(Y_u) - Abs(Y_u_old))
    Where (Abs(Y_u).Ne.0._dp) D_mat = D_mat / Abs(Y_u)

    Do i1=1,3
     If (D_mat(i1,i1).Gt.delta0) converge = .False.
     Do i2=i1+1,3
      If (D_mat(i1,i2).Gt.delta0) converge = .False.
      If (D_mat(i2,i1).Gt.delta0) converge = .False.
     End Do
    End Do
    Call Chop(Y_d)
    D_mat = Abs(Abs(Y_d) - Abs(Y_d_old))
    Where (Abs(Y_d).Ne.0._dp) D_mat = D_mat / Abs(Y_d)
    Do i1=1,3
     If (D_mat(i1,i1).Gt.delta0) converge = .False.
     Do i2=i1+1,3
      If (D_mat(i1,i2).Gt.delta0) converge = .False.
      If (D_mat(i2,i1).Gt.delta0) converge = .False.
     End Do
    End Do

    D_mat = Abs(Abs(CKM0) - Abs(V0CKM))
    Where (Abs(V0CKM).Ne.0._dp) D_mat = D_mat / Abs(V0CKM)
    Do i1=1,3
     If (D_mat(i1,i1).Gt.delta0) converge = .False.
     Do i2=i1+1,3
      If (D_mat(i1,i2).Gt.delta0) converge = .False.
      If (D_mat(i2,i1).Gt.delta0) converge = .False.
     End Do
    End Do

    If (converge) Exit

    V0CKM = CKM0
    V0CKMad = Transpose(Conjg(V0CKM))
    Y_l_old = Y_l
    Y_u_old = Y_u
    Y_d_old = Y_d

  Else ! .not.GenerationMixing

   !-------------------------------------------------------
   ! recalculate sfermion masses
   !-------------------------------------------------------
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (gSU2**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    Do i1=1,3
     mSneut2(i1) = Real(M2_L(i1,i1),dp)  + D_sneut
     If (mSneut2(i1).Lt.0._dp) Then
      kont = -229
      Call AddError(229)
      If (ErrorLevel.Ge.0) Then
        Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
        Write(ErrCan,*) 'mSneutrino^2 ',i1,' <= 0 :',mSneut2(i1)
        Write(ErrCan,*) 'setting it to 10.'
        mSneut2(i1) = 10._dp
      End If
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
     mSneut(i1) = Sqrt(mSneut2(i1))
    End Do
    RSneut = id3C
    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   !----------------
   ! sleptons
   !---------------
   RSlepton = 0
   Do i1 = 1,3
    Ml = M2_L(i1,i1)
    Mr = M2_E(i1,i1)
    A = T_l(i1,i1)
    Y = Y_l(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3_d, YL_l, YR_l, gSU2, gp &
                    &, kont, msf, msf2, Rsf)
    mSlept(2*(i1-1)+1:2*i1) = msf
    mSlept2(2*(i1-1)+1:2*i1) = msf2
    RSlept(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do
   !---------------
   ! down-squarks
   !---------------
   RSdown = 0
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_D(i1,i1)
    A = T_d(i1,i1)
    Y = Y_d(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3_d, YL_q, YR_d, gSU2, gp &
                    &, kont, msf, msf2, Rsf)
    mSdown(2*(i1-1)+1:2*i1) = msf
    mSdown2(2*(i1-1)+1:2*i1) = msf2
    RSdown(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do
   !---------------
   ! up-squarks
   !---------------
   RSup = 0
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_U(i1,i1)
    A = T_u(i1,i1)
    Y = Y_u(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3_u, YL_q, YR_u, gSU2, gp &
                    &, kont, msf, msf2, Rsf)
    mSup(2*(i1-1)+1:2*i1) = msf
    mSup2(2*(i1-1)+1:2*i1) = msf2
    RSup(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do

   !-------------------------------------------------------
   ! the chirally enhanced terms
   !-------------------------------------------------------
   epsD = 0._dp
   epsL = 0._dp
   If (Resum) Then
    Call Sigma_SM_chirally_enhanced(gauge, vevSM, mu, V0ckm, Y_l, Y_d, Y_u &
       & , Mi, T_l, T_d, T_u, M2_E, M2_L, M2_D, M2_Q, M2_U &
       & , epsD, epsL, epsD_FC, kont )
    If (kont.ne.0) then
     Iname = Iname - 1
     return
    End If 
   End If

   Do i1=1,3
     p2 = mf_d_DR(i1)**2
     Call Sigma_Fermion(p2, i1, mf_d_DR, Y_d, id3C, id3C,gSU2,gSU3,sinW2_DR  &
      & ,T3_d, e_d, mf_u_DR, Y_u, id3C, id3C, mSdown2, RSdown, mSup2         &
      & ,RSup, mglu, phase_glu, mN, mN2, N, mC, mC2, U, V, mS02, RS0         &
      & ,mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run, .True., .True., SigDown)
     If (Resum) then
      mf_d_DR(i1) = mf_d_DR_SM(i1) + Real(SigDown,dp) &
                  &                + oosqrt2 * epsD(i1) * vevSM(2) * Y_d(i1,i1)
      mf_d_DR(i1) = mf_d_DR(i1) / (1 + tanb * epsD(i1) )
     Else
      mf_d_DR(i1) = mf_d_DR_SM(i1) + Real(SigDown,dp)
     End If
     Y_d(i1,i1) = sqrt2 * mf_d_DR(i1) / vevSM(1)

     p2 = mf_u_DR(i1)**2
     If (i1.Lt.3) Then
      Call Sigma_Fermion(p2, i1, mf_u_dR, Y_u, id3C, id3C, gSU2, gSU3,sinW2_DR&
        & ,T3_u, e_u, mf_d_DR, Y_d, id3C, id3C,mSup2,RSup,mSdown2   &
        & ,RSdown, mglu, phase_glu, mN, mN2, N, mC, mC2, U, V, mS02, RS0    &
        & ,mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run, .True., .True., SigUp)
     Else
      Call Sigma_Fermion(p2, i1, mf_u_DR, Y_u, id3C, id3C, gSU2, gSU3,sinW2_DR&
        & ,T3_u, e_u, mf_d_DR, Y_d, id3C, id3C,mSup2,RSup,mSdown2   &
        & ,RSdown, mglu, phase_glu, mN, mN2, N, mC, mC2, U, V, mS02, RS0    &
        & ,mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run, .True., .False., SigUp)
     End If
     mf_u_DR(i1) = mf_u_DR_SM(i1) + Real(SigUp,dp)
     Y_u(i1,i1) = sqrt2 * mf_u_DR(i1) / vevSM(2)

     p2 = mf_l_DR(i1)**2
     Call Sigma_Fermion(p2, i1, mf_l, Y_l, id3C, id3C, gSU2, gSU3, sinW2_DR  &
        & , T3_d, e_e, mf_nu, Y_nu, id3C, id3C, mSlept2, RSlept          &
        & , mSneut2, RSneut, mglu, phase_glu, mN, mN2, N, mC, mC2, U &
        & , V, mS02, RS0, mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run           &
        & , .False., .True., SigLep)
     If (Resum) then
      mf_l_DR(i1) = mf_l_DR_SM(i1) + Real(SigLep,dp) &
                  &                + oosqrt2 * epsL(i1) * vevSM(2) * Y_l(i1,i1)
      mf_l_DR(i1) = mf_l_DR(i1) / (1 + tanb * epsL(i1) )
     Else
      mf_l_DR(i1) = mf_l_DR_SM(i1) + Real(SigLep,dp)
     End If
     Y_l(i1,i1)  = sqrt2 * mf_l_DR(i1) / vevSM(1)
    
    End Do

    If (    (      Abs((yuk_tau-y_l(3,3))/y_l(3,3)).Lt. delta0) &
        &    .And.(Abs((yuk_t-y_u(3,3))  /y_u(3,3)).Lt. delta0) &
        &    .And.(Abs((yuk_b-y_d(3,3))  /y_d(3,3)).Lt. delta0) ) Then
     converge = .True.
     Exit
    End If
   End If  ! GenerationMixing

   !--------------------------------------------------
   ! Either we have run into a numerical problem or
   ! perturbation theory breaks down
   !--------------------------------------------------

   If (    (Abs(mf_l_DR(3)/mf_l_mZ(3)).Lt.0.1_dp)  &
     & .Or.(Abs(mf_l_DR(3)/mf_l_mZ(3)).Gt.10._dp) ) Then
    Iname = Iname - 1
    kont = -405
    Call AddError(405)
    Return
   Else If (    (Abs(mf_d_DR(3)/mf_d_mZ(3)).Lt.0.1_dp)  &
          & .Or.(Abs(mf_d_DR(3)/mf_d_mZ(3)).Gt.10._dp) ) Then
    Iname = Iname - 1
    kont = -406
    Call AddError(406)
    Return
   Else If (    (Abs(mf_u_DR(3)/mf_u_mZ(3)).Lt.0.1_dp)  &
          & .Or.(Abs(mf_u_DR(3)/mf_u_mZ(3)).Gt.10._dp) ) Then
    Iname = Iname - 1
    kont = -407
    Call AddError(407)
    Return
   End If

  End Do ! i_loop

  If ((.Not.converge).And.FermionMassResummation) Then
   Write (ErrCan,*) 'Problem in subroutine BoundaryEWnew!!'
   Write (ErrCan,*) "After",i_loop-1,"iterations no convergence of Yukawas"
   Write (ErrCan,*) 'yuk_tau,yuk_l(3,3)',yuk_tau,y_l(3,3)
   Write (ErrCan,*) 'yuk_b,yuk_d(3,3)',yuk_b,y_d(3,3)
   Write (ErrCan,*) 'yuk_t,yuk_u(3,3)',yuk_t,y_u(3,3)
  End If
  !----------------------------------------------------------------
  ! the RGE paper defines the Yukawas transposed to my conventions
  !----------------------------------------------------------------
   If (GenerationMixing) Then
    If (YukawaScheme.Eq.1) Then
     Y_u = Matmul(Transpose(CKM0),Y_u) 
    Else
     Y_d = Matmul(Conjg(CKM0),Y_d) 
    End If
   End If

  Yl_mZ = Y_l
  Yd_mZ = Y_d
  Yu_mZ = Y_u
  vevs_DR_save = vevSM
  Y_u = Transpose(Y_u)
  Y_d = Transpose(Y_d)
  Y_l = Transpose(Y_l)
  sinW2_DR_mZ = sinW2_DR
  gauge(1) = Sqrt( 5._dp/3._dp) * gauge(1)
  gauge_mZ = gauge

  Call  CouplingsToG(gauge, y_l, y_d, y_u, g1)

  !----------------------------------------------
  ! resetting scale
  !----------------------------------------------
  mudim2 = SetRenormalizationScale(mudim2)

  Iname = Iname - 1

 Contains

  Real(dp) Function rho_2(r)
  Implicit None
   Real(dp), Intent(in) :: r
   Real(dp) :: r2, r3
   r2 = r*r
   r3 = r2*r
   rho_2 = 19._dp - 16.5_dp * r + 43._dp * r2 / 12._dp             &
       & + 7._dp * r3 / 120._dp                                    &
       & - Pi * Sqrt(r) * (4._dp - 1.5_dp * r + 3._dp * r2/32._dp  &
       &                  + r3/256._dp)                            &
       & - Pi2 * (2._dp - 2._dp * r + 0.5_dp * r2)                 &
       & - Log(r) * (3._dp * r - 0.5_dp * r2) 
  End  Function rho_2

  Subroutine Yukawas(mf, vev, tanb, SigS, SigL, SigR, Y, eps)
  !--------------------------------------------------------
  ! solves the matrix equation for Y by a transformation to
  ! a linear system of 9 equations in 9 unknowns
  ! written by Werner Porod, 19.03.03
  !--------------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: mf(3), vev, tanb
   Complex(dp), Dimension(3,3), Intent(in) :: SigS, SigL, SigR
   Complex(dp), Intent(in), Optional :: eps(3)
   Complex(dp), Intent(inout) :: Y(3,3)

   Integer :: i1
   Complex(dp), Dimension(3,3) :: mass, Y0, Yold

   !-------------------------------------
   ! first the mass matrix in DR scheme
   !-------------------------------------
   mass = ZeroC
   Do i1=1,3
    mass(i1,i1) = mf(i1)
   End Do
   !----------------------------------------
   ! setting up the equations
   !----------------------------------------
   Yold = Y
   Y0 = Y * vev * oosqrt2
   Y = mass + SigS 
   If (Present(eps)) Then ! resummation included
    Do i1=1,3
     Y(i1,i1) = ( Y(i1,i1) + oosqrt2 * eps(i1) * vev * tanb  * Yold(i1,i1)) &
            &   / (1 + tanb * eps(i1) )
    End Do
   End If
   Y = Y + Matmul(Transpose(SigL),Y0) + Matmul(Y0,SigR)

   Y(1,2) = 0._dp
   Y(1,3) = 0._dp
   Y(2,3) = 0._dp
   Y(2,1) = 0._dp
   Y(3,1) = 0._dp
   Y(3,2) = 0._dp
   ! this is possible as I can always rotate the R-fields
   Y = sqrt2 * Abs(Y) / vev 

  End Subroutine Yukawas

 End Subroutine BoundaryEWnew

 Subroutine BoundaryHS(g1,g2)
 !-----------------------------------------------------------------------
 ! calculates the  boundary at th high scale
 ! written by Werner Porod, 28.8.99
 ! last change: 28.8.99
 !     and back, putting in correct thresholds. For the first iteration
 !     only the first 6 couplings are included and a common threshold
 !     is used.
 ! 25.09.01: Portation to f90
 !  - the string scenario A: eqs. (3.11), (3.15) and (3.19)
 !    the string scenario B: eqs. (3.11), (3.16) and (3.20)
 !    from P.Binetruy at al., NPB 604, 32 (2001), hep-ph/0011081
 ! 31.10.02: including AMSB
 ! 06.12.02: including OI model of 
 !           P.Binetruy at al., NPB 604, 32 (2001), hep-ph/0011081
 !       -  Note, that I do have the opposite sign convention concerning
 !          the anomalous dimensions 
 !       - It is assumed that all terms of the form Ln(mu_R) vanish
 !-----------------------------------------------------------------------
  Implicit None

  Real(dp), Intent(in) :: g1(:)
  Real(dp), Intent(out) :: g2(:)

  Real(dp) :: gGMSB, fGMSB, ratio, Mhlf2(3), Mhlf1(3), gauge2(3), gauge4(3), M02
  Integer :: i1, i2, i3, ierr

  Real(dp) :: GammaH1, GammaH2, GammaGE(2,3), GammaGL(2,3), GammaGD(3,3)    &
     & , GammaGU(3,3), GammaGQ(3,3), GammaGH1(3), GammaGH2(3), GammaYH1     &
     & , GammaYH2, LnG2(3), fac, m15(3)
  Complex(dp), Dimension(3,3) :: GammaE, GammaL, GammaD, GammaU, GammaQ     &
     & , GammaYE, GammaYL, GammaYD, GammaYU, GammaYQ, GammaYQu, GammaYQd    &
     & , Ynu, d3, d3a, d3b, Yeff, UL, UR, MnuL5a, A_S_0, A_Z_0
  Complex(dp) :: d1, d2, wert
  Real(dp), Dimension(3) :: mf, YeGUT, YdGUT, YuGUT

  Complex(dp), Parameter :: MMnu(3,3) = Reshape( &
    &            Source = (/ (1._dp,0._dp), ZeroC, ZeroC, ZeroC, (2._dp,0._dp)  &
    &                        , ZeroC, ZeroC, ZeroC, (3._dp,0._dp) /), &
    &                     Shape = (/3, 3/)  )
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'BoundaryHS'

  If (Size(g1).Eq.57) Then
   Call GToCouplings(g1, gauge_0, Y_l_0, Y_d_0, Y_u_0)
   MnuL5a = id3c
  Else If (Size(g1).Eq.75) Then
    Call GToCouplings2(g1, gauge_0, Y_l_0, Y_nu_0, Y_d_0, Y_u_0)
    MnuL5a = Matmul( Matmul(Transpose(Y_nu_0),MMnu), Y_nu_0)

  Else If (Size(g1).Eq.79) Then
   Call GToCouplings4(g1, gauge_0, Y_l_0, d3, Y_d_0, Y_u_0, d1, d2)
   MnuL5a = id3c

  Else If (Size(g1).Eq.93) Then
   If ((Fixed_Nu_Yukawas).Or.(Ynu_at_MR3)) Then ! do not use Y_nu from running
                  !  but re-use the ones given from outside -> dummy argument Ynu
                  ! or in case of Ynu_at_MR3
    Call GToCouplings3(g1, gauge_0, Y_l_0, Ynu, Y_d_0, Y_u_0, mNuL5)
   Else
    Call GToCouplings3(g1, gauge_0, Y_l_0, Y_nu_0, Y_d_0, Y_u_0, mNuL5)
   End If
   If (Maxval(Abs(MnuL5)).Gt.0._dp) Then
    MnuL5a = MnuL5
   Else
    MnuL5a = id3c
   End If

  Else If (Size(g1).Eq.118) Then
   Call GToCouplings5(g1, gauge_0, Y_l_0, d3, Y_d_0, Y_u_0, d3a, d3b &
                    & , d1, d2, M15)
   MnuL5a = Matmul( Matmul(d3,MMnu), d3)

! Florian Staub Seesaw II+III
# ifdef SEESAWIII
  Else If (Size(g1).Eq.115) Then ! Seesaw II (SARAH)
   Call GToParameters115(g1, g10_H15,g20_H15,g30_H15, Yu0_H15, Yd0_H15, Ye0_H15 &
     & , Yt_H15,Ys_H15,Yz_H15,Lambda1,Lambda2)
   gauge_0(1) = g10_H15   
   gauge_0(2) = g20_H15
   gauge_0(3) = g30_H15
   Y_l_0 = Ye0_H15
   Y_d_0 = Yd0_H15
   Y_u_0 = Yu0_H15
   MnuL5a = Matmul( Matmul(Yt_H15,MMnu), Yt_H15)

  Else If (Size(g1).Eq.111) Then ! Seesaw III (SARAH)
   Call GToParameters111(g1, g10_H24,g20_H24,g30_H24, Yu0_H24, Yd0_H24, Ye0_H24 &
    & , Yb3_H24,Yw3_H24,Yx3_H24)

   gauge_0(1) = g10_H24   
   gauge_0(2) = g20_H24
   gauge_0(3) = g30_H24
   Y_l_0 = Ye0_H24
   Y_d_0 = Yd0_H24
   Y_u_0 = Yu0_H24
   MnuL5a = Matmul( Matmul(Transpose(Yw3_H24),MMnu), Yw3_H24)

# endif SEESAWIII
! Florian Staub Seesaw II+III

  Else
   Write(ErrCan,*) "Error in routine BoundaryHS"
   Write(ErrCan,*) "Size of g1",Size(g1)
   Call TerminateProgram
  End If  

  !--------------------------------------------------
  ! the following models need anomalous dimensions
  !--------------------------------------------------
  If (    (HighScaleModel.Eq.'AMSB').Or.(HighScaleModel.Eq.'Str_A')  &
     &.Or.(HighScaleModel.Eq.'Str_B').Or.(HighScaleModel.Eq.'Str_C') ) Then
   !-------------------------------------
   ! some General Definitions
   !-------------------------------------
   gauge2 = gauge_0**2
   gauge4 = gauge2**2
   GammaYE = oo16pi2 * 2._dp * Matmul( Y_l_0, Conjg(Transpose(Y_l_0)) )
   GammaYL = oo16pi2 * Matmul( Conjg(Transpose(Y_l_0)), Y_l_0 )
   GammaYD = oo16pi2 * 2._dp * Matmul( Y_d_0, Conjg(Transpose(Y_d_0)) )
   GammaYU = oo16pi2 * 2._dp * Matmul( Y_u_0, Conjg(Transpose(Y_u_0)) )
   GammaYQd = oo16pi2 * Matmul( Conjg(Transpose(Y_d_0)), Y_d_0 )
   GammaYQu = oo16pi2 * Matmul( Conjg(Transpose(Y_u_0)), Y_u_0 )
   GammaYQ = GammaYQd + GammaYQu
   !-----------------------------------------------------------------------
   ! anomalous dimension, attention, there are different sign conventions
   ! in the literature
   !------------------------------------------------------------------------
   GammaGH1(1) = - oo16pi2 * 0.3_dp * gauge2(1)
   GammaGH1(2) = - oo16pi2 * 1.5_dp * gauge2(2)
   GammaGH1(3) = 0._dp
   GammaGH2 = GammaGH1
   GammaGE(1,:) = - oo16pi2 * 1.2_dp * gauge2(1) 
   GammaGE(2,:) = 0._dp
   GammaGL(1,:) = GammaGH1(1) 
   GammaGL(2,:) = GammaGH1(2)
   GammaGD(1,:) = - oo16pi2 * 0.4_dp * gauge2(1) / 3._dp
   GammaGD(2,:) = 0._dp
   GammaGD(3,:) = - oo16pi2 * 8._dp * gauge2(3) / 3._dp 
   GammaGU(1,:) = - oo16pi2 * 1.6_dp * gauge2(1) / 3._dp
   GammaGU(2,:) = 0._dp
   GammaGU(3,:) = GammaGD(3,:)
   GammaGQ(1,:) = - oo16pi2 * 0.1_dp * gauge2(1) / 3._dp
   GammaGQ(2,:) = GammaGL(2,:)
   GammaGQ(3,:) = GammaGD(3,:)
 
   GammaH1 =  Sum(GammaGH1)
   GammaH2 = GammaH1
   GammaE = GammaYE 
   GammaL = GammaYL 
   GammaD = GammaYD 
   GammaU = GammaYU 
   GammaQ = GammaYQ 
   Do i1=1,3
    GammaH1 = GammaH1 + 3._dp * GammaYQd(i1,i1) + GammaYL(i1,i1)
    GammaH2 = GammaH2 + 3._dp * GammaYQu(i1,i1) 
    GammaE(i1,i1) = GammaE(i1,i1) + Sum(GammaGE(:,i1))
    GammaL(i1,i1) = GammaL(i1,i1) + Sum(GammaGL(:,i1))
    GammaD(i1,i1) = GammaD(i1,i1) + Sum(GammaGD(:,i1))
    GammaU(i1,i1) = GammaU(i1,i1) + Sum(GammaGU(:,i1))
    GammaQ(i1,i1) = GammaQ(i1,i1) + Sum(GammaGQ(:,i1))
   End Do
   GammaYH1 = GammaH1 - Sum(GammaGH1)
   GammaYH2 = GammaH2 - Sum(GammaGH2)
  End If

  !--------------------------------------------------
  ! now the boundary conditions
  !--------------------------------------------------
  If (HighScaleModel.Eq.'GMSB') Then
   M2_E_0 = ZeroC
   M2_L_0 = ZeroC
   M2_D_0 = ZeroC
   M2_Q_0 = ZeroC
   M2_U_0 = ZeroC

   ratio = Lambda / MlambdaS

   gGMSB = (1._dp+ratio)/ratio**2 * Log(1._dp + ratio)        &
       & + (1._dp-ratio)/ratio**2 * Log(1._dp - ratio)
   fGMSB = (1._dp+ratio)/ratio**2 * ( Log(1._dp + ratio)                &
       &                          - 2._dp * Li2(ratio/(1._dp+ratio) )   &
       &               + 0.5_dp * Li2(2._dp*ratio/(1._dp+ratio) ) )     &
       & + (1._dp-ratio)/ratio**2 * ( Log(1._dp - ratio)                &
       &               - 2._dp * Li2(ratio/(ratio-1._dp) )              &
       &               + 0.5_dp * Li2(2._dp*ratio/(ratio-1._dp) ) )

   Mhlf1 = gauge_0**2 * Lambda * oo16pi2
   Mhlf2 = Mhlf1**2
   Mi_0 = gGMSB * (n5plets + 3*n10plets) * Mhlf1

   If (Off_GMSB) Then
    Mhlf2(1) = 1.000581_dp * Mhlf2(1)
    Mhlf2(2) = 1.000937_dp * Mhlf2(2)
    Mhlf2(3) = 1.00154_dp * Mhlf2(3)
   End If
   M2_H_0 = fGMSB * (n5plets + 3*n10plets)                  &
         &       * ( 1.5_dp * Mhlf2(2) + 0.3_dp * Mhlf2(1) )
   M2_E_0(1,1) = fGMSB * (n5plets + 3*n10plets) * 1.2_dp * Mhlf2(1)
   M2_L_0(1,1) = M2_H_0(1)
   M2_D_0(1,1) = fGMSB * (n5plets + 3*n10plets)              &
         &       * (2._dp / 15._dp * Mhlf2(1) + 8._dp / 3._dp * Mhlf2(3) )
   M2_Q_0(1,1) = fGMSB * (n5plets + 3*n10plets)              &
      &  * (Mhlf2(1) / 30._dp + 1.5_dp * Mhlf2(2) + 8._dp / 3._dp * Mhlf2(3) )
   M2_U_0(1,1) = fGMSB * (n5plets + 3*n10plets)              &
         &       * (8._dp / 15._dp * Mhlf2(1) + 8._dp / 3._dp * Mhlf2(3) )
   Do i1=2,3
    M2_E_0(i1,i1) = M2_E_0(1,1)
    M2_L_0(i1,i1) = M2_L_0(1,1)
    M2_D_0(i1,i1) = M2_D_0(1,1)
    M2_Q_0(i1,i1) = M2_Q_0(1,1)
    M2_U_0(i1,i1) = M2_U_0(1,1)
   End Do

  Else If (HighScaleModel.Eq.'AMSB') Then
   M2_E_0 = ZeroC
   M2_L_0 = ZeroC
   M2_D_0 = ZeroC
   M2_Q_0 = ZeroC
   M2_U_0 = ZeroC
   Mi_0 = ZeroC
   A_l_0 = ZeroC
   A_d_0 = ZeroC
   A_u_0 = ZeroC
   !----------------------------------
   ! gaugino mass parameters 
   !----------------------------------
   Mi_0 = 0._dp
   Do i1=1,3
    Mi_0(i1) = m_32 * b_1(i1) * gauge2(i1) * oo16pi2
   End Do
   !----------------------------------------
   ! A parameter
   !----------------------------------------
   A_l_0 = ZeroC
   A_d_0 = ZeroC
   A_u_0 = ZeroC

   If (GenerationMixing) Then
    A_l_0 = Y_l_0 * GammaH1 + Matmul(GammaL, Conjg(Transpose(Y_l_0)) ) &
         &               + Matmul(Y_l_0, GammaE)
    A_d_0 = Y_d_0 * GammaH1 + Matmul(GammaQ, Conjg(Transpose(Y_d_0)) ) &
         &               + Matmul(Y_d_0, GammaD)
    A_u_0 = Y_u_0 * GammaH2 + Matmul(GammaQ, Conjg(Transpose(Y_u_0)) ) &
         &               + Matmul(Y_u_0, GammaU)
    A_l_0 = - m_32 * A_l_0
    A_d_0 = - m_32 * A_d_0
    A_u_0 = - m_32 * A_u_0
   Else ! .not. GenerationMixing
    Do i1=1,3
     A_l_0(i1,i1) = - m_32 * Y_l_0(i1,i1) * (GammaL(i1,i1)+GammaE(i1,i1)+GammaH1)
     A_d_0(i1,i1) = - m_32 * Y_d_0(i1,i1) * (GammaQ(i1,i1)+GammaD(i1,i1)+GammaH1)
     A_u_0(i1,i1) = - m_32 * Y_u_0(i1,i1) * (GammaQ(i1,i1)+GammaU(i1,i1)+GammaH2)
    End Do
   End If
   !----------------------------------------
   ! scalar parameters
   !----------------------------------------
   M02 = M0_amsb**2
   If (GenerationMixing) Then
    M2_E_0 = - oo16pi2 * m_32 * ( Matmul( Transpose(Conjg(Y_l_0)), A_l_0)   &
          &                    + Matmul( Transpose(Conjg(A_l_0)), Y_l_0)   )
    M2_D_0 = - oo16pi2 * m_32 * ( Matmul( Transpose(Conjg(Y_d_0)), A_d_0)   &
          &                    + Matmul( Transpose(Conjg(A_d_0)), Y_d_0)   )
    M2_U_0 = - oo16pi2 * m_32 * ( Matmul( Transpose(Conjg(Y_u_0)), A_u_0)   &
          &                    + Matmul( Transpose(Conjg(A_u_0)), Y_u_0)   )
    M2_L_0 = 0.5_dp * Transpose( M2_E_0 )
    M2_Q_0 = 0.5_dp * Transpose( M2_D_0 + M2_U_0 )
    Do i1=1,3
     M2_E_0(i1,i1) = M2_E_0(i1,i1) + M02                                      &
         & - (oo16pi2 * m_32)**2 * 198._dp / 25._dp * gauge4(1)
     M2_L_0(i1,i1) = M2_L_0(i1,i1) + M02                                      &
         & - (oo16pi2 * m_32)**2                                              &
         &               * (99._dp / 50._dp * gauge4(1) + 1.5_dp * gauge4(2))
     M2_D_0(i1,i1) = M2_D_0(i1,i1) + M02                                     &
         & - (oo16pi2 * m_32)**2                                             &
         &               * (22._dp / 25._dp * gauge4(1) - 8._dp * gauge4(3))
     M2_Q_0(i1,i1) = M2_Q_0(i1,i1) + M02                                     &
         & - (oo16pi2 * m_32)**2                                             &
         &          * (11._dp / 50._dp * gauge4(1)  + 1.5_dp * gauge4(2)     &
         &             - 8._dp * gauge4(3))
     M2_U_0(i1,i1) = M2_U_0(i1,i1) + M02                                     &
         & - (oo16pi2 * m_32)**2                                             &
         &               * (88._dp / 25._dp * gauge4(1) - 8._dp * gauge4(3))
    End Do

   Else ! .not. GenerationMixing
    Do i1=1,3
     M2_E_0(i1,i1) = M02                                         &
         & - (oo16pi2 * m_32)**2 * 198._dp / 25._dp * gauge4(1) &
         & - oo16pi2 * m_32 * 2._dp * Y_l_0(i1,i1) * A_l_0(i1,i1)
     M2_L_0(i1,i1) = M02                                                      &
         & - (oo16pi2 * m_32)**2                                              &
         &               * (99._dp / 50._dp * gauge4(1) + 1.5_dp * gauge4(2)) &
         & - oo16pi2 * m_32 * Y_l_0(i1,i1) * A_l_0(i1,i1)
     M2_D_0(i1,i1) = M02                                                      &
         & - (oo16pi2 * m_32)**2                                              &
         &               * (22._dp / 25._dp * gauge4(1) - 8._dp * gauge4(3))  &
         & - oo16pi2 * m_32 * 2._dp * Y_d_0(i1,i1) * A_d_0(i1,i1)
     M2_Q_0(i1,i1) = M02                                                     &
         & - (oo16pi2 * m_32)**2                                             &
         &          * (11._dp / 50._dp * gauge4(1)  + 1.5_dp * gauge4(2)     &
         &             - 8._dp * gauge4(3))                                  &
         & - oo16pi2 * m_32                                                  &
         &       * (Y_d_0(i1,i1) * A_d_0(i1,i1) + Y_u_0(i1,i1) * A_u_0(i1,i1) )
     M2_U_0(i1,i1) = M02                                                      &
         & - (oo16pi2 * m_32)**2                                              &
         &               * (88._dp / 25._dp * gauge4(1) - 8._dp * gauge4(3)) &
         & - oo16pi2 * m_32 * 2._dp * Y_u_0(i1,i1) * A_u_0(i1,i1)
    End Do
   End If

   M2_H_0 = M02                                                             &
       & - (oo16pi2 * m_32)**2                                              &
       &               * (99._dp / 50._dp * gauge4(1) + 1.5_dp * gauge4(2))
   Do i1=1,3
    M2_H_0(1) = M2_H_0(1) - oo16pi2 * m_32 * (3._dp *Y_d_0(i1,i1) *A_d_0(i1,i1) &
             &                             +  Y_l_0(i1,i1) * A_l_0(i1,i1)  )
    M2_H_0(2) = M2_H_0(2) - oo16pi2 * m_32 * 3._dp * Y_u_0(i1,i1) * A_u_0(i1,i1)
   End Do  

  Else If (HighScaleModel.Eq.'Str_A') Then
   !----------------------------------------
   ! Gaugino mass parameter
   !----------------------------------------
   Do i1=1,3
    Mi_0(i1) = sinT * (1 + oo8pi2 * g_s2 * sumC2(i1)) * phase_s * oosqrt_k_ss &
         &   + 2._dp * oo8pi2 * b_1(i1) * oosqrt3
    If (cosT.Ne.0._dp) Then
     Do i2=1,num_t
      Mi_0(i1) = Mi_0(i1)                                   &
          &    + oo8pi2 * cosT * ReG2ThetaT(i2) * phase_t(i2)   &
          &             * (delta_GS + 2._dp*b_1(i1) ) * oosqrt3
     End Do
    End If
    Mi_0(i1) = - 0.5_dp * sqrt3 * gauge2(i1) * m32 * Mi_0(i1)
   End Do
  !----------------------------------------
   ! A parameter
   !----------------------------------------
   AoY_l_0 = ZeroC
   AoY_d_0 = ZeroC
   AoY_u_0 = ZeroC

   wert = - oosqrt3 * k_s * sinT * oosqrt_k_ss 
  
   Do i1=1,3
    AoY_l_0(i1,i1) = -m32 * (GammaL(i1,i1) + GammaE(i1,i1) + GammaH1 + wert)
    AoY_d_0(i1,i1) = -m32 * (GammaQ(i1,i1) + GammaD(i1,i1) + GammaH1 + wert)
    AoY_u_0(i1,i1) = -m32 * (GammaQ(i1,i1) + GammaU(i1,i1) + GammaH2 + wert)
   End Do

   !----------------------------------------
   ! scalar mass parameters
   !----------------------------------------
   M2_E_0 = ZeroC
   M2_L_0 = ZeroC
   M2_D_0 = ZeroC
   M2_Q_0 = ZeroC
   M2_U_0 = ZeroC

   fac = sqrt3 * sinT * oosqrt_k_ss
   Do i1=1,3
    M2_E_0(i1,i1) = sinT2 - GammaE(i1,i1) 
    M2_L_0(i1,i1) = sinT2 - GammaL(i1,i1)
    M2_D_0(i1,i1) = sinT2 - GammaD(i1,i1)
    M2_Q_0(i1,i1) = sinT2 - GammaQ(i1,i1)
    M2_U_0(i1,i1) = sinT2 - GammaU(i1,i1)

   ! sin(Theta) part
    wert = - gauge2(1) * GammaGE(1,i1) * Real( phase_s,dp )            &
     &   + 0.5_dp * GammaYE(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
    M2_E_0(i1,i1) = M2_E_0(i1,i1) + fac * wert

    wert = - ( gauge2(1) * GammaGL(1,i1)                             &
     &       + gauge2(2) * GammaGL(2,i1) ) * Real( phase_s,dp )         &
     &   + 0.5_dp * GammaYL(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
    M2_L_0(i1,i1) = M2_L_0(i1,i1) + fac * wert

    wert = 0._dp
    Do i3=1,3
     wert = wert - gauge2(i3) * GammaGD(i3,i1)
    End Do
    wert = wert * Real( phase_s, dp )                                  &
     &   + 0.5_dp * GammaYD(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
    M2_D_0(i1,i1) = M2_D_0(i1,i1) + fac * wert

    wert = 0._dp
    Do i3=1,3
     wert = wert - gauge2(i3) * GammaGQ(i3,i1)
    End Do
    wert = wert * Real( phase_s, dp )                                    &
     &   + 0.5_dp * GammaYQ(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
    M2_Q_0(i1,i1) = M2_Q_0(i1,i1) + fac * wert

    wert = 0._dp
    Do i3=1,3
     wert = wert - gauge2(i3) * GammaGU(i3,i1)
    End Do
    wert = wert * Real( phase_s, dp )                                     &
     &   + 0.5_dp * GammaYU(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
    M2_U_0(i1,i1) = M2_U_0(i1,i1) + fac * wert

    ! overall part
    M2_E_0(i1,i1) = M2_E_0(i1,i1) * m32**2
    M2_L_0(i1,i1) = M2_L_0(i1,i1) * m32**2
    M2_D_0(i1,i1) = M2_D_0(i1,i1) * m32**2
    M2_Q_0(i1,i1) = M2_Q_0(i1,i1) * m32**2
    M2_U_0(i1,i1) = M2_U_0(i1,i1) * m32**2
   End Do

   M2_H_0(1) = SinT2 - GammaH1
   M2_H_0(2) = SinT2 - GammaH2
   Do i2=1,num_t
   ! sin(Theta) part
    wert = - Real( phase_s, dp )                                                &
     &            * ( gauge2(1) * GammaGH1(1) + gauge2(2) * GammaGH1(2) )   &
     &   + 0.5_dp * GammaYH1 * ( k_s * Conjg(phase_s) + k_sb * phase_s )
    M2_H_0(1) = M2_H_0(1) + fac * wert

    wert = - Real( phase_s,dp )                                         & 
     &   * ( gauge2(1) * GammaGH2(1) + gauge2(2) * GammaGH2(2) )   &
     &   + 0.5_dp * GammaYH2 * ( k_s * Conjg(phase_s) + k_sb * phase_s )
    M2_H_0(2) = M2_H_0(2) + fac * wert

   End Do

   M2_H_0 = M2_H_0 * m32**2

  Else If (HighScaleModel.Eq.'Str_B') Then
   !----------------------------------------
   ! Gaugino mass parameter
   !----------------------------------------
   Do i1=1,3
    Mi_0(i1) = sinT * (1 + oo8pi2 * g_s2 * sumC2(i1)) * phase_s * oosqrt_k_ss &
      &   + 2._dp * oo8pi2 * b_1(i1) * oosqrt3
    If (cosT.Ne.0._dp) Then
     Do i2=1,num_t
      Mi_0(i1) = Mi_0(i1)                                   &
      &    + oo8pi2 * cosT * ReG2ThetaT(i2) * phase_t(i2)   &
      &         * (delta_GS + 2._dp*b_1(i1) ) * oosqrt3
     End Do
    End If
    Mi_0(i1) = - 0.5_dp * sqrt3 * gauge2(i1) * m32 * Mi_0(i1)
   End Do
   !----------------------------------------
   ! A parameter
   !----------------------------------------
   AoY_l_0 = ZeroC
   AoY_d_0 = ZeroC
   AoY_u_0 = ZeroC

   wert = - k_s * sinT * oosqrt_k_ss 
   LnG2 = 0.5_dp * (Log( Gauge2 ) - 1._dp)
   i2 = 1
   Do i1=1,3
    AoY_l_0(i1,i1) = (GammaL(i1,i1) + GammaE(i1,i1) + GammaH1)            &
     &             * (1._dp + cosT * ReG2ThetaT(i2) ) / Sqrt3          &
     &    - sinT * ( (GammaGE(1,i1) + GammaGL(1,i1) +  GammaGH1(1))    &
     &                  * gauge2(1) * ( LnG2(1)  -  LnReDedekind(i2) ) &
     &             + (GammaGL(2,i1) +  GammaGH1(2))                    &
     &                  * gauge2(2) * ( LnG2(2) -  LnReDedekind(i2) )  &
     &             + (GammaYE(i1,i1) + GammaYL(i1,i1) + GammaH1) * k_s &
     &                *  LnReDedekind(i2)                              &
     &             ) * oosqrt_k_ss
    AoY_d_0(i1,i1) = - (GammaQ(i1,i1) + GammaD(i1,i1) + GammaH1)          &
     &             * (1._dp + cosT * ReG2ThetaT(i2) ) / Sqrt3          &
     &    - sinT * ( (GammaGD(1,i1) + GammaGQ(1,i1) +  GammaGH1(1))    &
     &                  * gauge2(1) * ( LnG2(1)  -  LnReDedekind(i2) ) &
     &             + (GammaGQ(2,i1) +  GammaGH1(2))                    &
     &                  * gauge2(2) * ( LnG2(2) -  LnReDedekind(i2) )  &
     &             + (GammaGD(3,i1) +  GammaGQ(3,i1))                  &
     &                  * gauge2(3) * ( LnG2(3) -  LnReDedekind(i2) )  &
     &             + (GammaYD(i1,i1) + GammaYQ(i1,i1) + GammaH1) * k_s &
     &                *  LnReDedekind(i2)                              &
     &             ) * oosqrt_k_ss
    AoY_u_0(i1,i1) = (GammaQ(i1,i1) + GammaU(i1,i1) + GammaH2)         & 
     &             * (1._dp + cosT * ReG2ThetaT(i2) ) / Sqrt3         &
     &    - sinT * ( (GammaGU(1,i1) + GammaGQ(1,i1) +  GammaGH2(1))    &
     &                  * gauge2(1) * ( LnG2(1)  -  LnReDedekind(i2) ) &
     &             + (GammaGQ(2,i1) +  GammaGH2(2))                    &
     &                  * gauge2(2) * ( LnG2(2) -  LnReDedekind(i2) )  &
     &             + (GammaGU(3,i1) +  GammaGQ(3,i1))                  &
     &                  * gauge2(3) * ( LnG2(3) -  LnReDedekind(i2) )  &
     &             + (GammaYU(i1,i1) + GammaYQ(i1,i1) + GammaH2) * k_s &
     &                *  LnReDedekind(i2)                              &
     &             ) * oosqrt_k_ss

    AoY_l_0(i1,i1) = - sqrt3 * m32 * (AoY_l_0(i1,i1) + wert)
    AoY_d_0(i1,i1) = - sqrt3 * m32 * (AoY_d_0(i1,i1) + wert)
    AoY_u_0(i1,i1) = - sqrt3 * m32 * (AoY_u_0(i1,i1) + wert)
   End Do

   !----------------------------------------
   ! scalar mass parameters
   !----------------------------------------
   M2_E_0 = ZeroC
   M2_L_0 = ZeroC
   M2_D_0 = ZeroC
   M2_Q_0 = ZeroC
   M2_U_0 = ZeroC
   M2_H_0 = 0._dp

   If (SinT.Ne.0._dp) Then
    i2=1
    LnG2 = Log( Gauge2 )
    Do i1=1,3
     ! sin(Theta) part
     fac = sqrt3 * sinT *  oosqrt_k_ss * (1._dp + cosT * ReG2ThetaT(i2) )
     wert = - gauge2(1) * GammaGE(1,i1) * Real( phase_s, dp )            &
      &   + 0.5_dp * GammaYE(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_E_0(i1,i1) = fac * wert

     wert = - ( gauge2(1) * GammaGL(1,i1)                            &
      &        + gauge2(2) * GammaGL(2,i1) ) * Real( phase_s,dp )         &
      &   + 0.5_dp * GammaYL(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_L_0(i1,i1) = fac * wert

     wert = 0._dp
     Do i3=1,3
      wert = wert - gauge2(i3) * GammaGD(i3,i1)
     End Do
     wert = wert * Real( phase_s, dp )                                   &
      &   + 0.5_dp * GammaYD(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_D_0(i1,i1) = M2_D_0(i1,i1) + fac * wert

     wert = 0._dp
     Do i3=1,3
      wert = wert - gauge2(i3) * GammaGQ(i3,i1)
     End Do
     wert = wert * Real( phase_s, dp )                                    &
      &   + 0.5_dp * GammaYQ(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_Q_0(i1,i1) = M2_Q_0(i1,i1) + fac * wert

     wert = 0._dp
     Do i3=1,3
      wert = wert - gauge2(i3) * GammaGU(i3,i1)
     End Do
     wert = wert * Real( phase_s, dp )                                     &
      &   + 0.5_dp * GammaYU(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_U_0(i1,i1) = M2_U_0(i1,i1) + fac * wert

     ! sin(Theta)^2  part
     wert = 1._dp - GammaE(i1,i1)                                      &
       & + GammaGE(1,i1) * (LnG2(1) -  LnReDedekind(i2))               &
       & + 2._dp * GammaYE(i1,i1) * LnReDedekind(i2)                   &
       & + ( 2.25_dp * gauge4(1) * GammaGE(1,i1)                       &
       &             * ( LnG2(1) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 3._dp * GammaYE(i1,i1) * k_s * k_sb * LnReDedekind(i2)    &
       &   )   / k_ss
     M2_E_0(i1,i1) = M2_E_0(i1,i1) + SinT2 * wert

     wert = 1._dp - GammaL(i1,i1)                                      &
       & + GammaGL(1,i1) * (LnG2(1) -  LnReDedekind(i2))               &
       & + GammaGL(2,i1) * (LnG2(2) -  LnReDedekind(i2))               &
       & + 2._dp * GammaYL(i1,i1) * LnReDedekind(i2)                   &
       & + ( 2.25_dp * gauge4(1) * GammaGL(1,i1)                       &
       &             * ( LnG2(1) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 2.25_dp * gauge4(2) * GammaGL(2,i1)                       &
       &             * ( LnG2(2) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 3._dp * GammaYL(i1,i1) * k_s * k_sb * LnReDedekind(i2)    &
       &   )   / k_ss
     M2_L_0(i1,i1) = M2_L_0(i1,i1) + SinT2 * wert

     wert = 1._dp - GammaD(i1,i1)                                      &
       & + GammaGD(1,i1) * (LnG2(1) -  LnReDedekind(i2))               &
       & + GammaGD(3,i1) * (LnG2(3) -  LnReDedekind(i2))               &
       & + 2._dp * GammaYD(i1,i1) * LnReDedekind(i2)                   &
       & + ( 2.25_dp * gauge4(1) * GammaGD(1,i1)                       &
       &             * ( LnG2(1) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 2.25_dp * gauge4(3) * GammaGD(3,i1)                       &
       &             * ( LnG2(3) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 3._dp * GammaYD(i1,i1) * k_s * k_sb * LnReDedekind(i2)    &
       &   )   / k_ss
     M2_D_0(i1,i1) = M2_D_0(i1,i1) + SinT2 * wert

     wert = 1._dp - GammaU(i1,i1)                                      &
       & + GammaGU(1,i1) * (LnG2(1) -  LnReDedekind(i2))               &
       & + GammaGU(3,i1) * (LnG2(3) -  LnReDedekind(i2))               &
       & + 2._dp * GammaYU(i1,i1) * LnReDedekind(i2)                   &
       & + ( 2.25_dp * gauge4(1) * GammaGU(1,i1)                       &
       &             * ( LnG2(1) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 2.25_dp * gauge4(3) * GammaGU(3,i1)                       &
       &             * ( LnG2(3) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 3._dp * GammaYU(i1,i1) * k_s * k_sb * LnReDedekind(i2)    &
       &   )   / k_ss
     M2_U_0(i1,i1) = M2_U_0(i1,i1) + SinT2 * wert

     wert = 1._dp - GammaQ(i1,i1)                                      &
       & + GammaGQ(1,i1) * (LnG2(1) -  LnReDedekind(i2))               &
       & + GammaGQ(2,i1) * (LnG2(2) -  LnReDedekind(i2))               &
       & + GammaGQ(3,i1) * (LnG2(3) -  LnReDedekind(i2))               &
       & + 2._dp * GammaYQ(i1,i1) * LnReDedekind(i2)                   &
       & + ( 2.25_dp * gauge4(1) * GammaGQ(1,i1)                       &
       &             * ( LnG2(1) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 2.25_dp * gauge4(2) * GammaGQ(2,i1)                       &
       &             * ( LnG2(2) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 2.25_dp * gauge4(3) * GammaGQ(3,i1)                       &
       &             * ( LnG2(3) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 3._dp * GammaYQ(i1,i1) * k_s * k_sb * LnReDedekind(i2)    &
       &   )   / k_ss
     M2_Q_0(i1,i1) = M2_Q_0(i1,i1) + SinT2 * wert

     ! overall part
     M2_E_0(i1,i1) = M2_E_0(i1,i1) * m32**2
     M2_L_0(i1,i1) = M2_L_0(i1,i1) * m32**2
     M2_D_0(i1,i1) = M2_D_0(i1,i1) * m32**2
     M2_Q_0(i1,i1) = M2_Q_0(i1,i1) * m32**2
     M2_U_0(i1,i1) = M2_U_0(i1,i1) * m32**2
    End Do

    Do i2=1,num_t
     fac = sqrt3 * sinT *  oosqrt_k_ss * (1._dp + cosT * ReG2ThetaT(i2) )
    ! sin(Theta) part
     wert = - Real( phase_s, dp )                                            &
      &            * ( gauge2(1) * GammaGH1(1) + gauge2(2) * GammaGH1(2) )   &
      &   + 0.5_dp * GammaYH1 * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_H_0(1) = M2_H_0(1) + fac * wert

     wert = - Real( phase_s, dp )                                         & 
      &   * ( gauge2(1) * GammaGH2(1) + gauge2(2) * GammaGH2(2) )   &
      &   + 0.5_dp * GammaYH2 * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_H_0(2) = M2_H_0(2) + fac * wert

     ! sin(Theta)^2  part
     wert = 1._dp - GammaH1                                            &
       & + GammaGH1(1) * (LnG2(1) -  LnReDedekind(i2))                 &
       & + GammaGH1(2) * (LnG2(2) -  LnReDedekind(i2))                 &
       & + 2._dp * GammaYH1 * LnReDedekind(i2)                         &
       & + ( 2.25_dp * gauge4(1) * GammaGH1(1)                         &
       &             * ( LnG2(1) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 2.25_dp * gauge4(2) * GammaGH1(2)                         &
       &             * ( LnG2(2) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 3._dp * GammaYH1 * k_s * k_sb * LnReDedekind(i2)          &
       &   )   / k_ss
     M2_H_0(1) = M2_H_0(1) + SinT2 * wert

     wert = 1._dp - GammaH2                                            &
       & + GammaGH2(1) * (LnG2(1) -  LnReDedekind(i2))                 &
       & + GammaGH2(2) * (LnG2(2) -  LnReDedekind(i2))                 &
       & + 2._dp * GammaYH2 * LnReDedekind(i2)                         &
       & + ( 2.25_dp * gauge4(1) * GammaGH2(1)                         &
       &             * ( LnG2(1) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 2.25_dp * gauge4(2) * GammaGH2(2)                         &
       &             * ( LnG2(2) + 5._dp / 3._dp + LnReDedekind(i2) )  &
       &   + 3._dp * GammaYH2 * k_s * k_sb * LnReDedekind(i2)          &
       &   )   / k_ss
     M2_H_0(2) = M2_H_0(2) + SinT2 * wert
    End Do

    M2_H_0 = M2_H_0 * m32**2
   End If ! SinT.ne.0._dp

  Else If (HighScaleModel.Eq.'Str_C') Then
   !----------------------------------------
   ! Gaugino mass parameter
   !----------------------------------------
   fac = oosqrt_k_ss * sqrt3
   Do i1=1,3
    Mi_0(i1) = sinT * (1._dp + oo8pi2 * g_s2 * sumC2(i1)) * phase_s * fac  &
           & + 2._dp * oo8pi2 * b_1(i1)
    If (cosT.Ne.0._dp) Then
     Do i2=1,num_t
      Mi_0(i1) = Mi_0(i1) + oo8pi2 * cosT * ReG2ThetaT(i2) * phase_t(i2)   &
            &         * (delta_GS + 2._dp * (b_1(i1) - SumC_O1(i1)) )
     End Do
    End If
    Mi_0(i1) = - 0.5_dp * gauge2(i1) * m32 * Mi_0(i1)
   End Do
   !----------------------------------------
   ! A parameter
   !----------------------------------------
   AoY_l_0 = ZeroC
   AoY_d_0 = ZeroC
   AoY_u_0 = ZeroC

   wert = - k_s * sinT * oosqrt_k_ss * sqrt3
   i2 = 1
   fac = cosT * ReG2ThetaT(i2)

   Do i1=1,3
    AoY_l_0(i1,i1) = GammaL(i1,i1) + GammaE(i1,i1) + GammaH1
    AoY_d_0(i1,i1) = GammaQ(i1,i1) + GammaD(i1,i1) + GammaH1
    AoY_u_0(i1,i1) = GammaQ(i1,i1) + GammaU(i1,i1) + GammaH2
    If (cosT.Ne.0._dp) Then
     AoY_l_0(i1,i1) = AoY_l_0(i1,i1)                                    &
               & + fac * (nE_ai(i2,i1) + nL_ai(i2,i1) + nH1_ai(i2) +3)
     AoY_d_0(i1,i1) = AoY_d_0(i1,i1)                                    &
               & + fac * (nD_ai(i2,i1) + nQ_ai(i2,i1) + nH1_ai(i2) +3)
     AoY_u_0(i1,i1) = AoY_u_0(i1,i1)                                    &
               & + fac * (nU_ai(i2,i1) + nQ_ai(i2,i1) + nH2_ai(i2) +3)
    End If

    AoY_l_0(i1,i1) = - m32 * (AoY_l_0(i1,i1) + wert)
    AoY_d_0(i1,i1) = - m32 * (AoY_d_0(i1,i1) + wert)
    AoY_u_0(i1,i1) = - m32 * (AoY_u_0(i1,i1) + wert)
   End Do

   !----------------------------------------
   ! scalar mass parameters
   !----------------------------------------
   M2_E_0 = ZeroC
   M2_L_0 = ZeroC
   M2_D_0 = ZeroC
   M2_Q_0 = ZeroC
   M2_U_0 = ZeroC
   M2_H_0 = 0._dp

   i2=1
   fac = sqrt3 * sinT * oosqrt_k_ss
   Do i1=1,3
     wert = - gauge2(1) * GammaGE(1,i1) * Real( phase_s, dp)    &
        & + 0.5_dp * GammaYE(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_E_0(i1,i1) = fac * wert                                         &
             &    + cosT * ReG2ThetaT(1) *  GammaYE(i1,i1)             &
             &           * (nE_ai(1,i1) + nL_ai(1,i1) + nH1_ai(1) +3)  &
             &    + CosT2 * nE_ai(1,i1)

     wert = - gauge2(1) * GammaGL(1,i1) * Real( phase_s, dp)    &
        & - gauge2(2) * GammaGL(2,i1) * Real( phase_s, dp)      &
        & + 0.5_dp * GammaYL(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_L_0(i1,i1) = fac * wert                                         &
             &    + cosT * ReG2ThetaT(1) *  GammaYL(i1,i1)             &
             &           * (nE_ai(1,i1) + nL_ai(1,i1) + nH1_ai(1) +3)  &
             &    + CosT2 * nL_ai(1,i1)

     wert = - gauge2(1) * GammaGD(1,i1) * Real( phase_s, dp)    &
        & - gauge2(3) * GammaGD(3,i1) * Real( phase_s, dp)      &
        & + 0.5_dp * GammaYD(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_D_0(i1,i1) = fac * wert                                         &
             &    + cosT * ReG2ThetaT(1) *  GammaYD(i1,i1)             &
             &           * (nD_ai(1,i1) + nQ_ai(1,i1) + nH1_ai(1) +3)  &
             &    + CosT2 * nD_ai(1,i1)

     wert = - gauge2(1) * GammaGQ(1,i1) * Real( phase_s, dp)    &
        & - gauge2(2) * GammaGQ(2,i1) * Real( phase_s, dp)      &
        & - gauge2(3) * GammaGQ(3,i1) * Real( phase_s, dp)      &
        & + 0.5_dp * GammaYQ(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_Q_0(i1,i1) = fac * wert                                         &
             &    + cosT * ReG2ThetaT(1) *  GammaYQd(i1,i1)             &
             &           * (nD_ai(1,i1) + nQ_ai(1,i1) + nH1_ai(1) +3)  &
             &    + cosT * ReG2ThetaT(1) *  GammaYQu(i1,i1)             &
             &           * (nU_ai(1,i1) + nQ_ai(1,i1) + nH2_ai(1) +3)  &
             &    + CosT2 * nQ_ai(1,i1)

     wert = - gauge2(1) * GammaGU(1,i1) * Real( phase_s, dp)    &
        & - gauge2(3) * GammaGU(3,i1) * Real( phase_s, dp)      &
        & + 0.5_dp * GammaYU(i1,i1) * ( k_s * Conjg(phase_s) + k_sb * phase_s )
     M2_U_0(i1,i1) = fac * wert                                         &
             &    + cosT * ReG2ThetaT(1) *  GammaYU(i1,i1)             &
             &           * (nU_ai(1,i1) + nQ_ai(1,i1) + nH2_ai(1) +3)  &
             &    + CosT2 * nU_ai(1,i1)

    ! overall part
    M2_E_0(i1,i1) = ( 1._dp + GammaE(i1,i1) + M2_E_0(i1,i1)) * m32**2
    M2_L_0(i1,i1) = ( 1._dp + GammaL(i1,i1) + M2_L_0(i1,i1)) * m32**2
    M2_D_0(i1,i1) = ( 1._dp + GammaD(i1,i1) + M2_D_0(i1,i1)) * m32**2
    M2_Q_0(i1,i1) = ( 1._dp + GammaQ(i1,i1) + M2_Q_0(i1,i1)) * m32**2
    M2_U_0(i1,i1) = ( 1._dp + GammaU(i1,i1) + M2_U_0(i1,i1)) * m32**2
   End Do

   Do i2=1,num_t
     wert = - gauge2(1) * GammaGH1(1) * Real( phase_s, dp)    &
        & - gauge2(2) * GammaGH1(2) * Real( phase_s, dp)      &
        & + 0.5_dp * GammaYH1 * ( k_s * Conjg(phase_s) + k_sb * phase_s )

     M2_H_0(1) = sqrt3 * sinT * wert * oosqrt_k_ss
     wert = 0._dp
     Do i1=1,3
      wert = wert + GammaYL(i1,i1) * (nE_ai(1,i1)+nL_ai(1,i1)+nH1_ai(1)+3) &
        &  + 3._dp * GammaYD(i1,i1) * (nD_ai(1,i1)+nQ_ai(1,i1)+nH1_ai(1)+3)
     End Do
     M2_H_0(1) = M2_H_0(1) + cosT * ReG2ThetaT(1) * wert + CosT2 * nH1_ai(1) 

     wert = - gauge2(1) * GammaGH2(1) * Real( phase_s, dp)    &
        & - gauge2(2) * GammaGH2(2) * Real( phase_s, dp)      &
        & + 0.5_dp * GammaYH2 * ( k_s * Conjg(phase_s) + k_sb * phase_s )

     M2_H_0(2) = sqrt3 * sinT * wert * oosqrt_k_ss
     wert = 0._dp
     Do i1=1,3
      wert = wert &
        &  + 3._dp * GammaYU(i1,i1) * (nU_ai(1,i1)+nQ_ai(1,i1)+nH2_ai(1)+3)
     End Do
     M2_H_0(2) = M2_H_0(2) + cosT * ReG2ThetaT(1) * wert + CosT2 * nH2_ai(1) 

   End Do

   M2_H_0(1) = (1._dp + GammaH1 + M2_H_0(1)) * m32**2
   M2_H_0(2) = (1._dp + GammaH2 + M2_H_0(2)) * m32**2

  End If  ! end of different models

  If (HighScaleModel.Eq.'Oscar') Then
   A_l_0 = AoY_l_0 * Y_l_0
   A_d_0 = Matmul(AoY_q_0, Y_d_0) + Matmul(Y_d_0, AoY_d_0)
   A_u_0 = Matmul(AoY_q_0, Y_u_0) + Matmul(Y_u_0, AoY_u_0)
  Else If (HighScaleModel.Eq.'AMSB') Then 
   AoY_l_0 = 0._dp
   AoY_d_0 = 0._dp
   AoY_u_0 = 0._dp
   Where(Abs(Y_l_0).Ne.0._dp) AoY_l_0 = A_l_0 / Y_l_0
   Where(Abs(Y_d_0).Ne.0._dp) AoY_d_0 = A_d_0 / Y_d_0
   Where(Abs(Y_u_0).Ne.0._dp) AoY_u_0 = A_u_0 / Y_u_0
  Else  If ((HighScaleModel.Ne.'GMSB').And.(HighScaleModel(1:3).Ne.'Str')) Then 
   !--------------------------------------------------------------
   ! check if parameters in the super CKM/PMNS basis are given 
   !--------------------------------------------------------------
   If (GenerationMixing) Then
    If (.Not.l_Al) Al_0_pmns = 0
    If (.Not.l_Ad) Ad_0_sckm = 0
    If (.Not.l_Au) Au_0_sckm = 0

    If (Model_Suchita) Then
     Call FermionMass(Y_l_0,1._dp,YeGUT,UR,UL,ierr)
     Call FermionMass(Y_d_0,1._dp,YdGUT,UR,UL,ierr)
     Call FermionMass(Y_u_0,1._dp,YuGUT,UR,UL,ierr)
     Call fv_boundary(YuGUT, YdGUT, YeGUT, M2_H_0(1), AoY_d_0(1,1) &
      & , M2Q_0_sckm, M2U_0_sckm, M2D_0_sckm, Au_0_sckm, Ad_0_sckm, Al_0_pmns)
    End If

    Call Switch_from_superCKM(Y_d_0, Y_u_0, Ad_0_sckm, Au_0_sckm, M2D_0_sckm  &
      &, M2Q_0_sckm, M2U_0_sckm, A_d_0, A_u_0, M2_D_0, M2_Q_0, M2_U_0, .True. )

    Call Switch_from_superPMNS(Y_l_0, MnuL5a, Al_0_pmns, M2E_0_pmns, M2L_0_pmns &
      & , A_l_0, M2_E_0, M2_L_0, .True. )

    If (.Not.l_Al) A_l_0 = AoY_l_0 * Y_l_0 
    If (.Not.l_Ad) A_d_0 = AoY_d_0 * Y_d_0 
    If (.Not.l_Au) A_u_0 = AoY_u_0 * Y_u_0 

    Call Chop(M2_D_0)
    Call Chop(M2_E_0)
    Call Chop(M2_L_0)
    Call Chop(M2_Q_0)
    Call Chop(M2_U_0)
    Call Chop(A_D_0)
    Call Chop(A_L_0)
    Call Chop(A_U_0)

   Else
    M2_D_0 = M2D_0_sckm
    M2_E_0 = M2E_0_pmns
    M2_L_0 = M2L_0_pmns
    M2_Q_0 = M2Q_0_sckm
    M2_U_0 = M2U_0_sckm

    A_l_0 = AoY_l_0 * Y_l_0
    A_d_0 = AoY_d_0 * Y_d_0
    A_u_0 = AoY_u_0 * Y_u_0

   End If

  End If

  mu_0 = 0._dp
  B_0 = 0._dp

  If (Size(g2).Eq.213) Then
   Call ParametersToG(gauge_0, Y_l_0, Y_d_0, Y_u_0, Mi_0, A_l_0, A_d_0, A_u_0 &
          & , M2_E_0, M2_L_0, M2_D_0, M2_Q_0, M2_U_0, M2_H_0, mu_0, B_0, g2)
  Else If (Size(g2).Eq.267) Then
   If (Ynu_eq_Yu)  Y_nu_0 =  Y_u_0
   If (Ynu_at_MR3) Then
    A_nu_0 = AoY_nu_0 * Ynu
    Call ParametersToG2(gauge_0, y_l_0, Ynu, y_d_0, y_u_0, Mi_0, A_l_0  &
       & , A_nu_0, A_d_0, A_u_0, M2_E_0, M2_L_0, M2_R_0, M2_D_0, M2_Q_0 &
       & , M2_U_0, M2_H_0, mu_0, B_0, g2)
   Else
    Call FermionMass(Y_l_0,1._dp,mf,UR,UL,ierr)
    Yeff = Matmul( Y_nu_0, Transpose(Conjg(UL)))
    A_nu_0 = AoY_nu_0 * Yeff

    Call ParametersToG2(gauge_0, y_l_0, Y_nu_0, y_d_0, y_u_0, Mi_0, A_l_0 &
       & , A_nu_0, A_d_0, A_u_0, M2_E_0, M2_L_0, M2_R_0, M2_D_0, M2_Q_0   &
       & , M2_U_0, M2_H_0, mu_0, B_0, g2)
   End If

  Else If (Size(g2).Eq.277) Then
   If (Ynu_at_MR3) Then
    Yeff = d3
   Else
    Call FermionMass(Y_l_0, 1._dp, mf, UR, UL, ierr)
    UL = Conjg(UL)
    Yeff = Matmul(UL,Matmul(Y_T_0,Transpose(UL)))
   End If
   A_l_0 = AoY_L_0 * Y_l_0
   A_T_0 = AoT_0 * Yeff
   Alam12_0 = Aolam12_0 * lam12_0
   Call ParametersToG4(gauge_0, y_l_0, Yeff, y_d_0, y_u_0, lam12_0(1)  &
      & , lam12_0(2), Mi_0, A_l_0, A_T_0, A_d_0, A_u_0, Alam12_0(1)     &
      & , Alam12_0(2), M2_E_0, M2_L_0, M2_D_0, M2_Q_0   &
      & , M2_U_0, M2_H_0, M2_T_0, mu_0, B_0, MnuL5, g2)

  Else If (Size(g2).Eq.285) Then
   If (Ynu_eq_Yu)  Y_nu_0 =  Y_u_0
   If (Ynu_at_MR3) Then
    A_nu_0 = AoY_nu_0 * Ynu
   Else
    Call FermionMass(Y_l_0, 1._dp, mf, UR, UL, ierr)
    Ynu = Matmul(Y_nu_0,Transpose(Conjg(UL)))
    A_nu_0 = AoY_nu_0 * Ynu
   End If
   Call ParametersToG3(gauge_0, y_l_0, Ynu, y_d_0, y_u_0, Mi_0, A_l_0   &
       & , A_nu_0, A_d_0, A_u_0, M2_E_0, M2_L_0, M2_R_0, M2_D_0, M2_Q_0 &
       & , M2_U_0, M2_H_0, mu_0, B_0, MnuL5, g2)

  Else If (Size(g2).Eq.356) Then
   If (Ynu_at_mR3) Then
    A_T_0 = AoT_0 * d3
    A_S_0 = AoT_0 * d3a
    A_Z_0 = AoT_0 * d3b
    
    Call ParametersToG5(gauge_0, y_l_0, d3, y_d_0, y_u_0, d3a, d3b     &
       & , lam12_0(1), lam12_0(2), Mi_0, A_l_0, A_T_0, A_d_0, A_u_0, A_S_0 &
       & , A_Z_0, Alam12_0(1), Alam12_0(2), M2_E_0, M2_L_0, M2_D_0, M2_Q_0 &
       & , M2_U_0, M2_H_0, M2_T_0, M2_T_0, M2_T_0, M15(1), M15(1), M15(1)  &
       & , mu_0, B_0, MnuL5, g2)
   Else
    Call FermionMass(Y_l_0, 1._dp, mf, UR, UL, ierr)
    UL = Conjg(UL)
    Yeff = Matmul(UL,Matmul(Y_T_0,Transpose(UL)))
    A_T_0 = AoT_0 * Yeff
    Alam12_0 = Aolam12_0 * lam12_0
    Call ParametersToG5(gauge_0, y_l_0, Yeff, y_d_0, y_u_0, Yeff, Yeff     &
       & , lam12_0(1), lam12_0(2), Mi_0, A_l_0, A_T_0, A_d_0, A_u_0, A_T_0 &
       & , A_T_0, Alam12_0(1), Alam12_0(2), M2_E_0, M2_L_0, M2_D_0, M2_Q_0 &
       & , M2_U_0, M2_H_0, M2_T_0, M2_T_0, M2_T_0, M15(1), M15(1), M15(1)  &   
       & , mu_0, B_0, MnuL5, g2)
   End If

! Florian Staub Seesaw II+III
# ifdef SEESAWIII
  Else If (Size(g2).Eq.365) Then
   If (.Not.Ynu_at_MR3) Then
    Call FermionMass(Y_l_0, 1._dp, mf, UR, UL, ierr)
    UL = Conjg(UL)
    Yeff = Matmul(UL,Matmul(Y_T_0,Transpose(UL)))
    Yt_h15 = Yeff
    Yz_h15 = Yeff
    Ys_h15 = Yeff
   Else
    Yeff = YT_h15
   End If

   MTM = MTM_GUT
   MSM = MTM
   MZM = MTM
   AMZM = AoY_l_0(1,1)*MTM  ! nachfragen
   AMTM = AMZM
   AMSM = AMZM
   AYt_h15 = AoY_l_0*Yeff
   AYz_h15 = AYt_h15
   AYs_h15 = AYt_h15
   ALambda1 = AoY_l_0(1,1)*Lambda1_gut
   ALambda2 = AoY_l_0(1,1)*Lambda2_gut
   mt2_H15 = m2d_0_sckm(1,1)
   mtb2_H15 = m2d_0_sckm(1,1)
   ms2_H15 = m2d_0_sckm(1,1)
   msb2_H15 = m2d_0_sckm(1,1)
   mz2_H15 = m2d_0_sckm(1,1)
   mzb2_H15 = m2d_0_sckm(1,1)

   Call ParametersToG365(g10_h15,g20_h15,g30_h15,Y_u_0,Y_d_0,Y_l_0,Yt_h15       &
     & ,Ys_h15,Yz_h15,Lambda1_gut,Lambda2_gut,mu_0,MTM,MZM,MSM,A_u_0,A_d_0      &
     & ,A_l_0,AYt_h15,AYs_h15,AYz_h15,ALambda1,ALambda2,B_0,AMTM,AMZM,AMSM      &
     & ,M2_Q_0,M2_L_0,M2_H_0(1),M2_H_0(2),M2_D_0,M2_U_0,M2_E_0,mt2_H15,mtb2_H15 &
     & ,ms2_H15,msb2_H15,mz2_H15,mzb2_H15,Mi_0(1),Mi_0(2),Mi_0(3),MnuL5a,g2)

  Else If (Size(g2).Eq.573) Then
   MWM3 = MWM3_gut
   MXM3 = MWM3
   MBM3 = MWM3
   MGM3 = MWM3
   Yb3_h24 = Yb3_h24_gut
   Yw3_h24 = Yb3_h24
   Yx3_h24 = Yb3_h24
   mHw32 = m2d_0_sckm
   mHb32 = mHw32
   mHx32 = mHw32
   mHxb32 = mHw32
   mHg32 = mHw32
   AYw3_H24 = AoY_l_0*Yw3_h24
   AYb3_H24 = AYw3_H24
   AYx3_H24 = AYw3_H24
   AMWM3 = AoY_l_0 * MWM3 ! nachfragen
   AMXM3 = AMWM3
   AMGM3 = AMWM3
   AMBM3 = AMWM3

   Call ParametersToG555(g10_H24,g20_H24,g30_H24, Yu0_H24,Yd0_H24,Ye0_H24      &
      & ,Yb3_H24,Yw3_H24,Yx3_H24,mu_0,MXM3, MWM3,MGM3, MBM3,A_u_0,A_d_0,A_l_0  & 
      & ,AYb3_H24,AYw3_H24,AYx3_H24, B_0,AMXM3,AMWM3,AMGM3,AMBM3,M2_Q_0,M2_L_0 &
      & ,M2_H_0(1),M2_H_0(2),M2_D_0, M2_U_0,M2_E_0,mHw32,mHg32,mHb32,mHx32     & 
      & ,mHxb32,Mi_0(1),Mi_0(2),Mi_0(3),zero33c,g2)
# endif SEESAWIII
! Florian Staub Seesaw II+III

  Else 
   Write(ErrCan,*) "Error in routine BoundaryHS"
   Write(ErrCan,*) "Size of g2",Size(g2)
   Call TerminateProgram
  End If

  Iname = Iname - 1

 End Subroutine BoundaryHS

 Subroutine BoundaryHS2(g1,g2)
 !-----------------------------------------------------------------------
 ! calculates the  boundary at th high scale
 ! written by Werner Porod, 28.8.99
 ! last change: 28.8.99
 !     and back, putting in correct thresholds. For the first iteration
 !     only the first 6 couplings are included and a ommon threshold
 !     is used.
 ! 25.09.01: Portation to f90
 !  - the string scenario A: eqs. (3.11), (3.15) and (3.19)
 !    the string scenario B: eqs. (3.11), (3.16) and (3.20)
 !    from P.Binetruy at al., NPB 604, 32 (2001), hep-ph/0011081
 ! 31.10.02: including AMSB
 ! 06.12.02: including OI model of 
 !           P.Binetruy at al., NPB 604, 32 (2001), hep-ph/0011081
 !       -  Note, that I do have the opposite sign convention concerning
 !          the anomalous dimensions 
 !       - It is assumed that all terms of the form Ln(mu_R) vanish
 !-----------------------------------------------------------------------
  Implicit None

  Real(dp), Intent(in) :: g1(:)
  Real(dp), Intent(out) :: g2(:)

  Integer :: ierr

  Real(dp) :: m15(3)
  Complex(dp), Dimension(3,3) :: Ynu, d3, d3a, d3b, Yeff, UL, UR, MnuL5a, A_S_0 &
      & , A_Z_0, A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U
  Complex(dp) :: Mi(3), mu, B
  Real(dp) :: mf(3), M2_H(2), M2_S_0(3,3), M2_Z_0(3,3)

  Complex(dp), Parameter :: MMnu(3,3) = Reshape( &
    &            Source = (/ (1._dp,0._dp), ZeroC, ZeroC, ZeroC, (2._dp,0._dp)  &
    &                        , ZeroC, ZeroC, ZeroC, (3._dp,0._dp) /), &
    &                     Shape = (/3, 3/)  )
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'BoundaryHS2'

  If (Size(g1).Ne.Size(g2)) Then
   Write(ErrCan,*) "Problem in BoundaryHS2, the sizes of the vectors do not"
   Write(ErrCan,*) "agree. size(g1):",Size(g1),"size(g2)",Size(g2)
   Call TerminateProgram
  End If

  If (Size(g1).Eq.213) Then
   Call GToParameters(g1, gauge_0, Y_l_0, Y_d_0, Y_u_0, Mi, A_l, A_d, A_u &
          & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B)
   MnuL5a = id3c

  Else If (Size(g1).Eq.267) Then
   Call GToParameters2(g1, gauge_0, y_l_0, Ynu, y_d_0, y_u_0, Mi, A_l, A_nu  &
               & , A_d, A_u, M2_E, M2_L, M2_R, M2_D, M2_Q, M2_U, M2_H, mu, B)
   MnuL5a = Matmul( Matmul(Transpose(Y_nu_0),MMnu), Y_nu_0)

  Else If (Size(g1).Eq.285) Then

   If (Ynu_at_MR3) Then 
    Call GToParameters3(g1, gauge_0, y_l_0, Y_nu_0, y_d_0, y_u_0, Mi, A_l &
       & , A_nu_0, A_d, A_u, M2_E, M2_L, M2_R, M2_D, M2_Q, M2_U, M2_H, mu &
       & , B, MnuL5)
   Else
    Call GToParameters3(g1, gauge_0, y_l_0, Ynu, y_d_0, y_u_0, Mi, A_l &
       & , A_nu_0, A_d, A_u, M2_E, M2_L, M2_R, M2_D, M2_Q, M2_U, M2_H, mu &
       & , B, MnuL5)
   End If
   If (Maxval(Abs(MnuL5)).Gt.0._dp) Then
    MnuL5a = MnuL5
   Else
    MnuL5a = id3c
   End If

  Else If (Size(g1).Eq.356) Then

   Call GToParameters5(g1, gauge_0, y_l_0, d3, y_d_0, y_u_0, d3a, d3b     &
       & , lam12_0(1), lam12_0(2), Mi, A_l, A_T_0, A_d, A_u, A_S_0, A_Z_0 &
       & , Alam12_0(1), Alam12_0(2), M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H   &
       & , M2_T_0, M2_S_0, M2_Z_0, M15(1), M15(2), M15(3), mu, B, MnuL5a)

! Florian Staub Seesaw II+III
# ifdef SEESAWIII
  Else If (Size(g1).Eq.365) Then ! Seesaw II (SARAH)
   Call GToParameters365(g1, g10_H15,g20_H15,g30_H15, Yu0_H15, Yd0_H15, Ye0_H15 &
     & , Yt_H15,Ys_H15,Yz_H15,Lambda1,Lambda2,mu,MTM,MZM,MSM,A_u,A_d      &
     & , A_l,AYt_h15,AYs_h15,AYz_h15,ALambda1,ALambda2,B,AMTM,AMZM,AMSM     &
     & ,M2_Q,M2_L,M2_H(1),M2_H(2),M2_D,M2_U,M2_E,mt2_H15,mtb2_H15 &
     & ,ms2_H15,msb2_H15,mz2_H15,mzb2_H15,Mi(1),Mi(2),Mi(3),MnuL5a)
   gauge_0(1) = g10_H15   
   gauge_0(2) = g20_H15
   gauge_0(3) = g30_H15
   Y_l_0 = Ye0_H15
   Y_d_0 = Yd0_H15
   Y_u_0 = Yu0_H15

  Else If (Size(g1).Eq.555) Then ! Seesaw III (SARAH)
   Call GToParameters555(g1, g10_H24,g20_H24,g30_H24, Yu0_H24, Yd0_H24, Ye0_H24 &
    & , Yb3_H24,Yw3_H24,Yx3_H24,mu,MXM3, MWM3,MGM3, MBM3,A_u,A_d,A_l  & 
      & ,AYb3_H24,AYw3_H24,AYx3_H24, B,AMXM3,AMWM3,AMGM3,AMBM3,M2_Q,M2_L &
      & ,M2_H(1),M2_H(2),M2_D, M2_U,M2_E,mHw32,mHg32,mHb32,mHx32     & 
      & ,mHxb32,Mi(1),Mi(2),Mi(3),MnuL5a)

   gauge_0(1) = g10_H24   
   gauge_0(2) = g20_H24
   gauge_0(3) = g30_H24
   Y_l_0 = Ye0_H24
   Y_d_0 = Yd0_H24
   Y_u_0 = Yu0_H24

# endif SEESAWIII
! Florian Staub Seesaw II+III

  Else
   Write(ErrCan,*) "Error in routine BoundaryHS2"
   Write(ErrCan,*) "Size of g1",Size(g1)
   Call TerminateProgram
  End If  

  !--------------------------------------------------
  ! now the boundary conditions
  !--------------------------------------------------------------
  ! check if parameters in the super CKM/PMNS basis are given 
  !--------------------------------------------------------------
  If (GenerationMixing) Then
   If (.Not.l_Al) Al_0_pmns = 0
   If (.Not.l_Ad) Ad_0_sckm = 0
   If (.Not.l_Au) Au_0_sckm = 0

   Call Switch_from_superCKM(Y_d_0, Y_u_0, Ad_0_sckm, Au_0_sckm, M2D_0_sckm  &
      &, M2Q_0_sckm, M2U_0_sckm, A_d_0, A_u_0, M2_D_0, M2_Q_0, M2_U_0, .True. )

   Call Switch_from_superPMNS(Y_l_0, MnuL5a, Al_0_pmns, M2E_0_pmns, M2L_0_pmns &
     & , A_l_0, M2_E_0, M2_L_0, .True. )

   If (.Not.l_Al) A_l_0 = AoY_l_0 * Y_l_0 
   If (.Not.l_Ad) A_d_0 = AoY_d_0 * Y_d_0 
   If (.Not.l_Au) A_u_0 = AoY_u_0 * Y_u_0 

  Else
   M2_D_0 = M2D_0_sckm
   M2_E_0 = M2E_0_pmns
   M2_L_0 = M2L_0_pmns
   M2_Q_0 = M2Q_0_sckm
   M2_U_0 = M2U_0_sckm

   A_l_0 = AoY_l_0 * Y_l_0
   A_d_0 = AoY_d_0 * Y_d_0
   A_u_0 = AoY_u_0 * Y_u_0

  End If

  If (Sum(in_extpar(1,:)).Gt.0) Mi_0(1) = Mi(1)
  If (Sum(in_extpar(2,:)).Gt.0) Mi_0(2) = Mi(2)
  If (Sum(in_extpar(3,:)).Gt.0) Mi_0(3) = Mi(3)

  If (Sum(in_extpar(21,:)).Gt.0) M2_H_0(1) = M2_H(1)
  If (Sum(in_extpar(22,:)).Gt.0) M2_H_0(2) = M2_H(2)

  If ( ((Sum(in_extpar(23,:)).Gt.0).And.(Sum(in_extpar(24,:)).Gt.0) )         &
   & .Or. ((Sum(in_extpar(23,:)).Gt.0).And.(Sum(in_extpar(26,:)).Gt.0) ) ) Then
   M2_H_0 = M2_H
   mu_0 = mu
   B_0 = B
  Else
   mu_0 = 0._dp
   B_0 = 0._dp
  End If

  If (Size(g2).Eq.213) Then
   Call ParametersToG(gauge_0, Y_l_0, Y_d_0, Y_u_0, Mi_0, A_l_0, A_d_0, A_u_0 &
          & , M2_E_0, M2_L_0, M2_D_0, M2_Q_0, M2_U_0, M2_H_0, mu_0, B_0, g2)

  Else If (Size(g2).Eq.267) Then
   If (Ynu_eq_Yu)  Y_nu_0 =  Y_u_0
   If (Ynu_at_MR3) Then
    A_nu_0 = AoY_nu_0 * Ynu
    Call ParametersToG2(gauge_0, y_l_0, Ynu, y_d_0, y_u_0, Mi_0, A_l_0  &
       & , A_nu_0, A_d_0, A_u_0, M2_E_0, M2_L_0, M2_R_0, M2_D_0, M2_Q_0 &
       & , M2_U_0, M2_H_0, mu_0, B_0, g2)
   Else
    Call FermionMass(Y_l_0,1._dp,mf,UR,UL,ierr)
    Yeff = Matmul( Y_nu_0, Transpose(Conjg(UL)))
    A_nu_0 = AoY_nu_0 * Yeff

    Call ParametersToG2(gauge_0, y_l_0, Y_nu_0, y_d_0, y_u_0, Mi_0, A_l_0 &
       & , A_nu_0, A_d_0, A_u_0, M2_E_0, M2_L_0, M2_R_0, M2_D_0, M2_Q_0   &
       & , M2_U_0, M2_H_0, mu_0, B_0, g2)
   End If

  Else If (Size(g2).Eq.277) Then
   If (Ynu_at_MR3) Then
    Yeff = d3
   Else
    Call FermionMass(Y_l_0, 1._dp, mf, UR, UL, ierr)
    UL = Conjg(UL)
    Yeff = Matmul(UL,Matmul(Y_T_0,Transpose(UL)))
   End If
   A_l_0 = AoY_L_0 * Y_l_0
   A_T_0 = AoT_0 * Yeff
   Alam12_0 = Aolam12_0 * lam12_0
   Call ParametersToG4(gauge_0, y_l_0, Yeff, y_d_0, y_u_0, lam12_0(1)  &
      & , lam12_0(2), Mi_0, A_l_0, A_T_0, A_d_0, A_u_0, Alam12_0(1)     &
      & , Alam12_0(2), M2_E_0, M2_L_0, M2_D_0, M2_Q_0   &
      & , M2_U_0, M2_H_0, M2_T_0, mu_0, B_0, MnuL5, g2)

  Else If (Size(g2).Eq.285) Then
   If (Ynu_at_MR3) Then
    A_nu_0 = AoY_nu_0 * Ynu
   Else
    Call FermionMass(Y_l_0, 1._dp, mf, UR, UL, ierr)
    Ynu = Matmul(Y_nu_0,Transpose(Conjg(UL)))
    A_nu_0 = AoY_nu_0 * Ynu
   End If
   Call ParametersToG3(gauge_0, y_l_0, Ynu, y_d_0, y_u_0, Mi_0, A_l_0   &
       & , A_nu_0, A_d_0, A_u_0, M2_E_0, M2_L_0, M2_R_0, M2_D_0, M2_Q_0 &
       & , M2_U_0, M2_H_0, mu_0, B_0, MnuL5, g2)

  Else If (Size(g2).Eq.356) Then
   If (Ynu_at_mR3) Then
    A_T_0 = AoT_0 * d3
    A_S_0 = AoT_0 * d3a
    A_Z_0 = AoT_0 * d3b
    
    Call ParametersToG5(gauge_0, y_l_0, d3, y_d_0, y_u_0, d3a, d3b     &
       & , lam12_0(1), lam12_0(2), Mi_0, A_l_0, A_T_0, A_d_0, A_u_0, A_S_0 &
       & , A_Z_0, Alam12_0(1), Alam12_0(2), M2_E_0, M2_L_0, M2_D_0, M2_Q_0 &
       & , M2_U_0, M2_H_0, M2_T_0, M2_T_0, M2_T_0, M15(1), M15(1), M15(1)  &
       & , mu_0, B_0, MnuL5, g2)
   Else
    Call FermionMass(Y_l_0, 1._dp, mf, UR, UL, ierr)
    UL = Conjg(UL)
    Yeff = Matmul(UL,Matmul(Y_T_0,Transpose(UL)))
    A_T_0 = AoT_0 * Yeff
    Alam12_0 = Aolam12_0 * lam12_0
    Call ParametersToG5(gauge_0, y_l_0, Yeff, y_d_0, y_u_0, Yeff, Yeff     &
       & , lam12_0(1), lam12_0(2), Mi_0, A_l_0, A_T_0, A_d_0, A_u_0, A_T_0 &
       & , A_T_0, Alam12_0(1), Alam12_0(2), M2_E_0, M2_L_0, M2_D_0, M2_Q_0 &
       & , M2_U_0, M2_H_0, M2_T_0, M2_T_0, M2_T_0, M15(1), M15(1), M15(1)  &   
       & , mu_0, B_0, MnuL5, g2)
   End If

! Florian Staub Seesaw II+III
# ifdef SEESAWIII
  Else If (Size(g2).Eq.365) Then
   If (.Not.Ynu_at_MR3) Then
    Call FermionMass(Y_l_0, 1._dp, mf, UR, UL, ierr)
    UL = Conjg(UL)
    Yeff = Matmul(UL,Matmul(Y_T_0,Transpose(UL)))
    Yt_h15 = Yeff
    Yz_h15 = Yeff
    Ys_h15 = Yeff
   End If

   MTM = MTM_GUT
   MSM = MTM
   MZM = MTM
   AMZM = AoY_l_0(1,1)*MTM  ! nachfragen
   AMTM = AMZM
   AMSM = AMZM
   AYt_h15 = AoY_l_0*Yeff
   AYz_h15 = AYt_h15
   AYs_h15 = AYt_h15
   ALambda1 = AoY_l_0(1,1)*Lambda1_gut
   ALambda2 = AoY_l_0(1,1)*Lambda2_gut
   mt2_H15 = m2d_0_sckm(1,1)
   mtb2_H15 = m2d_0_sckm(1,1)
   ms2_H15 = m2d_0_sckm(1,1)
   msb2_H15 = m2d_0_sckm(1,1)
   mz2_H15 = m2d_0_sckm(1,1)
   mzb2_H15 = m2d_0_sckm(1,1)

   Call ParametersToG365(g10_h15,g20_h15,g30_h15,Y_u_0,Y_d_0,Y_l_0,Yt_h15       &
     & ,Ys_h15,Yz_h15,Lambda1_gut,Lambda2_gut,mu_0,MTM,MZM,MSM,A_u_0,A_d_0      &
     & ,A_l_0,AYt_h15,AYs_h15,AYz_h15,ALambda1,ALambda2,B_0,AMTM,AMZM,AMSM      &
     & ,M2_Q_0,M2_L_0,M2_H_0(1),M2_H_0(2),M2_D_0,M2_U_0,M2_E_0,mt2_H15,mtb2_H15 &
     & ,ms2_H15,msb2_H15,mz2_H15,mzb2_H15,Mi_0(1),Mi_0(2),Mi_0(3),MnuL5a,g2)

  Else If (Size(g2).Eq.573) Then
   MWM3 = MWM3_gut
   MXM3 = MWM3
   MBM3 = MWM3
   MGM3 = MWM3
   Yb3_h24 = Yb3_h24_gut
   Yw3_h24 = Yb3_h24
   Yx3_h24 = Yb3_h24
   mHw32 = m2d_0_sckm
   mHb32 = mHw32
   mHx32 = mHw32
   mHxb32 = mHw32
   mHg32 = mHw32
   AYw3_H24 = AoY_l_0*Yw3_h24
   AYb3_H24 = AYw3_H24
   AYx3_H24 = AYw3_H24
   AMWM3 = AoY_l_0 * MWM3 ! nachfragen
   AMXM3 = AMWM3
   AMGM3 = AMWM3
   AMBM3 = AMWM3

   Call ParametersToG555(g10_H24,g20_H24,g30_H24, Yu0_H24,Yd0_H24,Ye0_H24      &
      & ,Yb3_H24,Yw3_H24,Yx3_H24,mu_0,MXM3, MWM3,MGM3, MBM3,A_u_0,A_d_0,A_l_0  & 
      & ,AYb3_H24,AYw3_H24,AYx3_H24, B_0,AMXM3,AMWM3,AMGM3,AMBM3,M2_Q_0,M2_L_0 &
      & ,M2_H_0(1),M2_H_0(2),M2_D_0, M2_U_0,M2_E_0,mHw32,mHg32,mHb32,mHx32     & 
      & ,mHxb32,Mi_0(1),Mi_0(2),Mi_0(3),zero33c,g2)
# endif SEESAWIII
! Florian Staub Seesaw II+III

  Else 
   Write(ErrCan,*) "Error in routine BoundaryHS2"
   Write(ErrCan,*) "Size of g2",Size(g2)
   Call TerminateProgram
  End If

  Iname = Iname - 1

 End Subroutine BoundaryHS2

 Subroutine Calculate_gi_Yi(Q, gauge, tanb, Y_l, Y_d, Y_u, mC, U, V, mN, N     &
    & , mS02, RS0, mP02, RP0, mSpm, mSpm2, RSpm, mDsquark, mDsquark2, RDsquark &
    & , mUsquark, mUsquark2, RUsquark, mSlepton, mSlepton2, RSlepton           &
    & , mSneutrino2 , RSneutrino, mGlu, phase_glu                              &
    & , alphaMZ, AlphaS_5, mf_d_SM, mf_u_SM, CKM, mf_l_SM, kont         &
    & , Mnu5, mf_nu, PMNS)
 !-----------------------------------------------------------------------
 ! Calculates gauge and yukawa couplings from the DRbar scheme to
 ! pole and MSbar quantities. The SUSY mass vectors and mixing matrices
 ! are not fixed so that the model can easily be extended
 ! input:
 !  Q ................ scale, where all the parameters and running masses
 !                     are given
 !  gauge ............ 3-vector for gauge couplings: i=1 -> U(1)
 !                                                   i=2 -> SU(2)
 !                                                   i=3 -> SU(3)
 !  tanb ............. tan(beta)
 !  Y_l .............. lepton Yukawa couplings
 !  Y_d .............. d-quarks Yukawa couplings
 !  Y_u .............. u-quarks Yukawa couplings
 !  mC ............... chargino masses
 !  U,V .............. chargino mixing matrices
 !  mN ............... neutralino masses
 !  N ................ neutralino mixing matrices
 !  mS02 ............. masses squared of the neutral scalar Higgs bosons
 !  RS0 .............. mixing matrix of the neutral scalar Higgs bosons
 !  mP02 ............. masses squared of the neutral pseudo-scalar Higgs bosons
 !  RP0 .............. mixing matrix of the neutral pseudo-scalar Higgs bosons
 !  mSpm2 ............ masses squared of the charged scalar Higgs bosons
 !  RSpm ............. mixing matrix of the charged scalar Higgs bosons
 !  mDsquark ......... masses of d-squarks
 !  mDsquark2 ........ masses of d-squarks squared
 !  RDsquark ......... mixing matrix of d-squarks
 !  mUsquark ......... masses of u-squarks
 !  mUsquark2 ........ masses of u-squarks squared
 !  RUsquark ......... mixing matrix of u-squarks
 !  mSlepton ......... masses of sleptons
 !  mSlepton2 ........ masses of sleptons squared
 !  RSlepton ......... mixing matrix of sleptons
 !  mSneutrino2 ...... masses of sneutrinos squared
 !  RSneutrino ....... mixing matrix of sneutrinos
 !  mGlu ............. gluino mass
 !  phase_glu ........ phase of gluino mass parameter
 ! input, optional
 !  mNu5 ............. dim-5 operator for neutrino masses
 ! output:
 !  AlphaMZ .......... electromagnetic coupling alpha(Q) at the scale Q
 !                     without t-quark but including W-boson
 !                     (corresponding to the usual alpha(mZ) definition
 !  AlphaS_5 ......... strong coupling alpha_S(Q) at the scale Q with
 !                     5 flavours of quarks
 !  mf_d_SM .......... MSbar d-quark masses: i=1 -> m_d(2 GeV)
 !                                           i=2 -> m_s(2 GeV)
 !                                           i=3 -> mb(mb)
 !  mf_u_SM .......... MSbar u-quark masses: i=1 -> m_u(2 GeV)
 !                                           i=2 -> m_c(2 GeV)
 !                                           i=3 -> pole mass of t-quark
 !  CKM .............. CKM matrix at scale Q
 !  mf_l_SM .......... lepton pole masses: i=1 -> m_e
 !                                         i=2 -> m_mu
 !                                         i=3 -> m_tau
 !  kont ............. if 0 everything is fine, if =!0 -> problem
 ! output, optional but must be given once mNu5 is given!
 !  mf_nu ............ neutrino masses
 !  PMNS ............. PMNS matrix
 ! written by Werner Porod, 10.02.2010
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: Q, gauge(3), tanb  ! scale, DRbar gauge coupling
  Complex(dp), Dimension(3,3), Intent(in) :: Y_u, Y_d, Y_l ! Yukawas
  Real(dp), Intent(inout) :: mSpm2(:)
  Real(dp), Intent(in) :: mC(:), mN(:), mSpm(:), mUsquark(:), mDsquark(:)     &
    & , mSlepton(:), mUsquark2(:), mDsquark2(:), mSlepton2(:), mSneutrino2(:) &
    & , mS02(:), mP02(:), RP0(:,:), mglu, RS0(:,:)
  Complex(dp), Intent(in) :: U(:,:), V(:,:), N(:,:), RSpm(:,:)           &
    & , RDsquark(:,:), RUsquark(:,:), RSlepton(:,:), RSneutrino(:,:)     &
    & , phase_glu
  Complex(dp), Optional, Intent(in) :: Mnu5(3,3)
  Real(dp), Intent(out) :: AlphaMZ  ! alpha MSbar
  Real(dp), Intent(out) :: AlphaS_5 ! alpha_s MSbar, 5 flavours
  Complex(dp), Intent(out) :: CKM(3,3)
  Integer, Intent(inout) :: kont
  Real(dp), Intent(out), Dimension(3) :: mf_l_SM, mf_d_SM, mf_u_SM
  Real(dp), Intent(out), Optional :: mf_nu(3)
  Complex(dp), Optional, Intent(out) :: PMNS(3,3)

  Integer :: n_S0, i1
  Real(dp) :: mC2( Size(mC) ), mN2(Size(mN) )
  Real(dp) :: test, gU1, gU1sq, gSU2, gSU2sq, sinW2, vev2    &
    & , vevs_DR(2), mZ2_mZ, CosW2SinW2, p2, gSU3, gSU3sq, SigQCD
  Complex(dp) :: dmZ2, SigLep, Sigdown, SigUp
  Complex(dp), Dimension(3,3) :: SigS_u, sigR_u, SigL_u, SigS_d, SigR_d, SigL_d &
    & , SigS_l, sigR_l, SigL_l, uU_L_T, uU_R_T, uD_L_T, uD_R_T, uL_L_T, uL_R_T
  Real(dp), Parameter :: e_d=-1._dp/3._dp, e_u=2._dp/3._dp, e_e=-1._dp &
    & , T3_d=-0.5_dp, T3_u=0.5_dp
  Complex(dp), Parameter :: Y_nu(3,3) = ZeroC

  Real(dp) :: sumI, AlphaSDR , DeltaAlpha, alpha_in, mZ2_run, mW2_run
  Real(dp), Dimension(3) :: mf_u_DR, mf_u_DR_SM, mf_d_DR &
                        &  , mf_d_DR_SM, mf_l_DR, mf_l_DR_SM, mf_l2 &
                        &  , mf_d2, mf_u2
  Complex(dp), Dimension(3,3) :: uL_u, uR_u, uL_d, uR_d, uL_l, uR_l, Unu, MnuL5
  !---------------------------------------------------
  ! needed for threshold corrections and rge running
  !---------------------------------------------------
  Real(dp) :: g9(9), g10(10), as_nf, as_nf_minus_1, as_nf_d_pi, aem , tz, dt  &
     & , mb_guess, mb_check
  Logical :: converge
  Real(dp), Parameter :: zeta3 = 1.202056903159594285399738161511449990765_dp &
     & , zeta4 = 1.082323233711138191516003696541167902775_dp                 &
     & , B4 = -1.762800087073770864061897634679818807215_dp                   &
     & , c_as(2) = (/ 11._dp / 72._dp                                         &
     &             , 58067._dp / 13824._dp - 82043._dp * zeta3 / 27648._dp /) &
     & , c_m(2) = (/ 89._dp / 432._dp                                         &
     &            , 713._dp/486._dp - B4 / 36._dp - 221._dp * zeta3 / 288._dp &
     &          + 1.25_dp * zeta4  /)

  Iname = Iname + 1
  nameOfUnit(Iname) = "Calculate_gi_Yi"
  !----------------------------------------
  ! checking if masses squared are positiv
  !----------------------------------------
  If (Min(Minval(mUsquark2), Minval(mDSquark2), Minval(mSlepton2)           &
     &    ,Minval(mSneutrino2), Minval(mS02), Minval(mP02), Minval(mSpm2))  &
     & .Lt. 0._dp ) Then
   kont = -401
   Call AddError(401)
   Iname = Iname - 1
   Return
  End If

  mf_nu = 0._dp ! only for initialisation
  n_s0 = Size(mS02)
  !-------------------------------------------------------------------
  ! setting renormalisation scale to m_Z, because the RGEs start there
  !-------------------------------------------------------------------
  test = SetRenormalizationScale(Q**2)
  !-------------------------------------------------------------------
  ! initialization of LoopMasses
  !-------------------------------------------------------------------
  Call  SetLoopMassModel(Size(mC), Size(mN), n_s0, n_s0, Size(mSpm) &
                       & , Size(mSlepton), Size(msneutrino2))
  !------------
  ! couplings
  !------------
  gU1 = gauge(1)
  gU1sq = gU1**2
  gSU2 = gauge(2)
  gSU2sq = gSU2**2
  gSU3 = gauge(3)
  gSU3sq = gSU3**2
  sinW2 = gU1sq / (gU1sq + gSu2sq)
  CosW2SinW2 = (1._dp - sinW2) * sinW2
  alpha_in = oo4pi * sinW2 * gSu2sq
  !------------
  ! vevs
  !------------
  vev2 =  mZ2 * CosW2SinW2 / (pi * alpha_in)
  vevs_DR(1) = Sqrt(vev2 / (1._dp+tanb**2) )
  vevs_DR(2) = tanb * vevs_DR(1)

  If (GenerationMixing) Then
   !---------------------
   ! running masses at Q
   !---------------------
   Call FermionMass(Y_l,vevs_DR(1),mf_l_DR,uL_L_T,uL_R_T,kont)
   Call FermionMass(Y_d,vevs_DR(1),mf_d_DR,uD_L_T,uD_R_T,kont)
   Call FermionMass(Y_u,vevs_DR(2),mf_u_DR,uU_L_T,uU_R_T,kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
  Else
   Do i1=1,3
    mf_l_DR(i1) = oosqrt2 * vevs_DR(1) * Abs(Y_l(i1,i1))
    mf_d_DR(i1) = oosqrt2 * vevs_DR(1) * Abs(Y_d(i1,i1))
    mf_u_DR(i1) = oosqrt2 * vevs_DR(2) * Abs(Y_u(i1,i1))
   End Do
  End If
  mf_l2 =  mf_l_DR**2
  mf_d2 =  mf_d_DR**2
  mf_u2 =  mf_u_DR**2
  !-------------
  ! alpha(mZ)
  !-------------
  alphaMZ = AlphaEwMS(Q, alpha_in, mSpm, mUsquark, mDSquark, mSlepton &
          &          , mC, mf_u_DR(3))
  !---------------
  ! alpha_s(mZ)
  !---------------
  sumI = 0
  Do i1=1,6
   sumI = sumI + Log( mDSquark(i1) / Q ) + Log( mUSquark(i1) / Q )
  End Do

  AlphaSDR = oo4pi * gSu3sq

  DeltaAlpha = 0.5_dp - 2._dp * Log(mf_u_DR(3)/Q ) / 3._dp &
           & - sumI / 6._dp  - 2._dp * Log(mGlu /Q)
  DeltaAlpha = AlphaSDR * DeltaAlpha / ( 2._dp * Pi)
  
  AlphaS_5 = AlphaSDR / (1._dp + DeltaAlpha)
  !-----------------
  ! vevs
  !-----------------
  vev2 =  mZ2 * CosW2SinW2 / (pi * alpha_in)
  vevs_DR(1) = Sqrt(vev2 / (1._dp+tanb**2) )
  vevs_DR(2) = tanb * vevs_DR(1)
  mC2 = mC**2
  mN2 = mN**2
  Call PiZZT1(mZ2, gSU2, sinW2, vevs_DR, mZ2, mW2, mS02, RS0, mP02, RP0   &
   & , mSpm2, RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton, mUsquark2 &
   & , RUsquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2, mC, mC2, U, V    &
   & , mN, mN2, N, dmZ2)
  mZ2_mZ = Real(dmZ2+mZ2,dp)
  If (mZ2_mZ.Lt.0._dp) Then
   Iname = Iname - 1
   kont = -413
   Call AddError(413) 
   Return
  End If
  mZ2_run = mZ2_mZ
  mW2_run = mZ2_mZ * (1._dp - sinW2)
  If (GenerationMixing) Then
   !---------------------
   ! running masses at Q
   !---------------------
   Call FermionMass(Y_l,vevs_DR(1),mf_l_DR,uL_L_T,uL_R_T,kont)
   Call FermionMass(Y_d,vevs_DR(1),mf_d_DR,uD_L_T,uD_R_T,kont)
   Call FermionMass(Y_u,vevs_DR(2),mf_u_DR,uU_L_T,uU_R_T,kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
  Else
   Do i1=1,3
    mf_l_DR(i1) = oosqrt2 * vevs_DR(1) * Abs(Y_l(i1,i1))
    mf_d_DR(i1) = oosqrt2 * vevs_DR(1) * Abs(Y_d(i1,i1))
    mf_u_DR(i1) = oosqrt2 * vevs_DR(2) * Abs(Y_u(i1,i1))
   End Do
  End If
  mf_l2 =  mf_l_DR**2
  mf_d2 =  mf_d_DR**2
  mf_u2 =  mf_u_DR**2
  !---------------------------------------
  ! recalculation, using running masses
  !---------------------------------------
  Call PiZZT1(mZ2, gSU2, sinW2, vevs_DR, mZ2_mZ, mW2_run, mS02, RS0        &
   & , mP02, RP0, mSpm2, RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton   &
   & , mUsquark2, RUsquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2, mC, mC2 &
   & , U, V, mN, mN2, N  &
   & , dmZ2)
  mZ2_mZ = Real(dmZ2+mZ2,dp)
  If (mZ2_mZ.Lt.0._dp) Then
    Iname = Iname - 1
    kont = -413
    Call AddError(413) 
    Return
  End If
  mZ2_run = mZ2_mZ
  mW2_run = mZ2_mZ * (1._dp - sinW2)

  vev2 =  mZ2_mZ * CosW2SinW2 / (pi * alpha_in)
  vevs_DR(1) = Sqrt(vev2 / (1._dp+tanb**2) )
  vevs_DR(2) = tanb * vevs_DR(1)

  !----------------------------
  ! fermion masses, CKM, PMNS
  !----------------------------

  If (GenerationMixing) Then
   !---------------------
   ! running masses at Q
   !---------------------
   Call FermionMass(Y_l,vevs_DR(1),mf_l_DR,uL_L_T,uL_R_T,kont)
   Call FermionMass(Y_d,vevs_DR(1),mf_d_DR,uD_L_T,uD_R_T,kont)
   Call FermionMass(Y_u,vevs_DR(2),mf_u_DR,uU_L_T,uU_R_T,kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   p2 = 0._dp ! for off-diagonal elements
   !---------------------------------------------------
   ! 1-loop SUSY contributions + QCD & QED for t-quark
   !---------------------------------------------------
   ! u-quarks
   Call Sigma_Fermion3(p2, mf_u_DR, Y_u, uU_L_T, uU_R_T, gSU2, gSU3, sinW2  &
       & , T3_u, e_u, mf_d_DR, Y_d, uD_L_T, uD_R_T, mUSquark2, RUsquark     &
       & , mDSquark2, RDsquark, mglu , phase_glu, mN, mN2, N, mC, mC2, U, V &
       & , mS02, RS0, mP02, RP0, mSpm2 , RSpm, mZ2_run, mW2_run, .True.     &
       & , SigS_u, SigL_u, SigR_u, SigQCD)
   Call Masses(Y_u, vevs_DR(2), SigS_u, SigL_u, SigR_u, .False.    &
            & , uL_u, uR_u, mf_u_DR_SM, kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   
   mf_u_SM(3) = mf_u_DR_SM(3)- SigQCD
   ! d-quarks
   Call Sigma_Fermion3(p2, mf_d_DR, Y_d, uD_L_T, uD_R_T, gSU2, gSU3, sinW2 &
       & , T3_d, e_d, mf_u_DR, Y_u, uU_L_T, uU_R_T, mDSquark2, RDsquark       &
       & , mUSquark2, RUsquark,  mglu , phase_glu, mN, mN2, N, mC, mC2, U, V  &
       & , mS02, RS0, mP02, RP0, mSpm2 , RSpm, mZ2_run, mW2_run , .True.      &
       & , SigS_d, SigL_d, SigR_d)
   Call Masses(Y_d, vevs_DR(1), SigS_d, SigL_d, SigR_d, .True.    &
            & , uL_d, uR_d, mf_d_DR_SM, kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If

   ! leptons
   Call Sigma_Fermion3(p2, mf_l_DR, Y_l, uL_L_T, uL_R_T, gSU2, gSU3, sinW2    &
       & , T3_d, e_e, mf_nu, Y_nu, id3C, id3C, mSlepton2, RSlepton           &
       & , mSneutrino2, RSneutrino, mglu , phase_glu, mN, mN2, N, mC, mC2, U &
       & , V, mS02, RS0, mP02, RP0, mSpm2 , RSpm, mZ2_run, mW2_run, .False.  &
       & , SigS_l, SigL_l, SigR_l)
   Call Masses(Y_l, vevs_DR(1), SigS_l, SigL_l, SigR_l, .True.    &
            & , uL_l, uR_l, mf_l_DR_SM, kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If

   CKM = Matmul(uL_u, Transpose(Conjg(uL_D)) )
   ! neutrinos
   If (Present(mNu5)) Then
    MnuL5 = Mnu5 * vevs_DR(2)**2
    Call NeutrinoMass_1L(MnuL5, gU1, gSu2, Y_l, mC2, U, V, mN2, N &
           & , mSlepton2, Rslepton, mSneutrino2, Rsneut, mf_nu, Unu, kont)
    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
    PMNS = Matmul(uL_L, Transpose(Conjg(Unu)))
   End If

  Else ! GenerationMixing
   Do i1=1,3
    mf_u_DR(i1) = Y_u(i1,i1) * vevs_DR(2) * oosqrt2
    mf_d_DR(i1) = Y_d(i1,i1) * vevs_DR(1) * oosqrt2
    mf_l_DR(i1) = Y_l(i1,i1) * vevs_DR(1) * oosqrt2

    p2 = mf_d_DR(i1)**2
    Call Sigma_Fermion(p2, i1, mf_d_DR, Y_d, id3C, id3C,gSU2,gSU3,sinW2   &
      & ,T3_d, e_d, mf_u_DR, Y_u, id3C, id3C, mDSquark2, RDSquark, mUSquark2  &
      & ,RUSquark, mglu, phase_glu, mN, mN2, N, mC, mC2, U, V, mS02, RS0      &
      & ,mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run, .True., .True., SigDown)
    mf_d_DR_SM(i1) = mf_d_DR(i1) * (1- Real(SigDown,dp) / mf_d_DR(i1) )

    p2 = mf_u_DR(i1)**2
    If (i1.Lt.3) Then
     Call Sigma_Fermion(p2, i1, mf_u_dR, Y_u, id3C, id3C, gSU2, gSU3,sinW2&
       & ,T3_u, e_u, mf_d_DR, Y_d, id3C, id3C,mUSquark2,RUSquark,mDSquark2   &
       & ,RDSquark, mglu, phase_glu, mN, mN2, N, mC, mC2, U, V, mS02, RS0    &
       & ,mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run, .True., .True., SigUp)
    Else
     Call Sigma_Fermion(p2, i1, mf_u_DR, Y_u, id3C, id3C, gSU2, gSU3,sinW2&
       & ,T3_u, e_u, mf_d_DR, Y_d, id3C, id3C,mUSquark2,RUSquark,mDSquark2   &
       & ,RDSquark, mglu, phase_glu, mN, mN2, N, mC, mC2, U, V, mS02, RS0    &
       & ,mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run, .True., .False., SigUp)
    End If
    mf_u_DR_SM(i1) = mf_u_DR(i1) - Real(SigUp,dp)

    p2 = mf_l_DR(i1)**2
    Call Sigma_Fermion(p2, i1, mf_l, Y_l, id3C, id3C, gSU2, gSU3, sinW2  &
        & , T3_d, e_e, mf_nu, Y_nu, id3C, id3C, mSlepton2, RSlepton          &
        & , mSneutrino2, RSneutrino, mglu, phase_glu, mN, mN2, N, mC, mC2, U &
        & , V, mS02, RS0, mP02, RP0, mSpm2, RSpm, mZ2_run, mW2_run           &
        & , .False., .True., SigLep)
    mf_l_DR_SM(i1) = mf_l_DR(i1) * (1- Real(SigLep,dp) / mf_l_DR(i1) )
    
   End Do
  End If ! GenerationMixing

  !--------------------------------------------------------------------------
  ! shifting light fermion masses to MS-scheme, only gluon and photon part
  ! except for m_t
  !--------------------------------------------------------------------------
  mf_l_SM = &
    &    mf_l_DR_SM / (1._dp - oo8pi2 *3._dp *(gU1sq-gSu2sq)/16._dp)
  mf_d_SM = mf_d_DR_SM / (1._dp - AlphaSDR / (3._dp*pi)                  &
         &               - 23._dp * AlphaSDR**2 / (72._dp * Pi2 )         &
         &               + oo8pi2 * 3._dp * gSu2sq / 16._dp        &
         &               - oo8pi2 * 13._dp * gU1sq / 144._dp  )
  mf_u_SM(1:2) = mf_u_DR_SM(1:2)  / (1._dp - AlphaSDR / (3._dp*pi)       &
         &               - 23._dp * AlphaSDR**2 / (72._dp * Pi2 )         &
         &               + oo8pi2 * 3._dp * gSu2sq / 16._dp        &
         &               - oo8pi2 * 7._dp * gU1sq / 144._dp  )


  !---------------------------------------
  ! now we have to find mb(mb)
  !---------------------------------------
  mb_guess = mf_d_SM(3)
  converge = .False.
  Do i1=1,40
   g10(1) = Sqrt(4._dp * Pi * alphaS_5)
   g10(2) = Sqrt(4._dp * Pi * alphaMZ)
   g10(3:5) = mf_l_SM
   g10(6:7) = mf_u_SM(1:2)
   g10(8:10) = mf_d_SM
   tz = Log( mb_guess / Q)
   dt = tz / 50._dp
   Call odeint(g10, 10, 0._dp, tz, 1.e-7_dp, dt, 0._dp, RGE10_SM, kont)
   as_nf_d_Pi = oo4pi * g10(1)**2 / Pi
   mb_check = g10(10) * (1._dp + as_nf_d_Pi**2 * (c_m(1) + c_m(2)*as_nf_d_Pi)) 

   If ((Abs(mb_guess-mb_check) / mb_guess).Lt.1.e-5_dp) Then
    converge = .True.
    Exit
   End If
   mb_guess = mb_check
  End Do

  If (.Not.converge) Then
    kont = -413
    Call AddError(413) 
  End If

  g9 = g10(1:9)
  mf_d_SM(3) = mb_check

  !---------------------------
  ! and now from m_b to 2 GeV
  !---------------------------
  ! thresholds
  !------------
  g9(6:9) = g9(6:9) * (1._dp + as_nf_d_Pi**2 * (c_m(1) + c_m(2)*as_nf_d_Pi)) 
  as_nf = oo4pi * g10(1)**2
  as_nf_minus_1 = as_nf *(1._dp+ as_nf_d_pi**2 *(c_as(1) +c_as(2)*as_nf_d_pi))
  g9(1) = Sqrt(4._dp * Pi * as_nf_minus_1)
  tz = Log(2._dp/mb_check) 
  dt = tz / 50._dp

  Call odeint(g9, 9, 0._dp, tz, 1.e-7_dp, dt, 0._dp, RGE10_SM, kont)

  aem = g9(2)**2 / (4._dp * Pi**2)
  mf_l_SM = g9(3:5) / (1._dp - aem)
  mf_u_SM(1:2) = g9(6:7)
  mf_d_SM(1:2) = g9(8:9)
  !----------------------------------------------
  ! resetting scale
  !----------------------------------------------
  test = SetRenormalizationScale(test)

  Iname = Iname - 1

 Contains

  Subroutine Masses(Yin, vev, SigS, SigL, SigR, ReSum, uL, uR, mf, kont)
  !--------------------------------------------------------
  ! solves the matrix equation for Y by a transformation to
  ! a linear system of 9 equations in 9 unknowns
  ! written by Werner Porod, 19.03.03
  !--------------------------------------------------------
  Implicit None
   Integer, Intent(inout) :: kont
   Real(dp), Intent(in) :: vev
   Complex(dp), Dimension(3,3), Intent(in) :: Yin, SigS, SigL, SigR
   Logical, Intent(in) :: ReSum
   Real(dp), Intent(out) :: mf(3)
   Complex(dp), Intent(out) :: uL(3,3), uR(3,3)

   Complex(dp), Dimension(3,3) :: Y, mass, f, invY

   !-------------------------------------
   ! first the mass matrix in DR scheme
   !-------------------------------------
   Y = Yin * vev * oosqrt2

   If (ReSum) Then
    kont = 0
    Call chop(Y)
    invY = Y
    Call gaussj(kont,invY,3,3)
    If (kont.Ne.0) Return

    f = id3C - Matmul(SigS,invY) - Transpose(SigL) - Matmul(Y,Matmul(SigR,invY))

     mass = Matmul(f,Y)

   Else

    mass = Y - SigS - Matmul(Transpose(SigL),Y) - Matmul(Y,SigR)

   End If

   Call FermionMass(mass, sqrt2, mF, uL, uR,kont)

  End Subroutine Masses

 End Subroutine Calculate_gi_Yi

 Subroutine FirstGuess(phase_mu, tanb, Mi, M_E2, M_L2, A_e, M_D2 &
                    & , M_Q2, M_U2, A_d, A_u, mu, BImu, M_H2, gU1, gSU2 &
                    & , Y_l, Y_d, Y_u, vevs, mP02, mP0, kont, Take_old_data, delta, name)
 !-----------------------------------------------------------------------
 ! calculates approximate values of the electroweak MSSM parameters,
 ! saving one run of subroutine sugra by running 1-loop RGEs
 ! written by Werner Porod,  08.10.01
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: phase_mu
  Real(dp), Intent(in) :: tanb
  Complex(dp), Intent(out) :: Mi(3), M_E2(3,3), M_L2(3,3), A_e(3,3)           &
     & , M_D2(3,3), M_Q2(3,3), M_U2(3,3), A_d(3,3), A_u(3,3)                  &
     & , mu, BImu, Y_l(3,3), Y_d(3,3), Y_u(3,3)
  Real(dp), Intent(out) :: M_H2(2), mP02(2), mP0(2), gSU2, gU1, vevS(2)
  Logical, Intent(in), Optional :: Take_old_data
  Character(len=*) , Intent(in), Optional :: name
  Real(dp), Intent(in), Optional :: delta

  Real(dp) :: sinb2, cosb2, abs_mu2, abs_mu
  Integer :: i1, kont

  Real(dp) :: gauge(3), vev, g1(57), g0(213), mgut, mudim, cosW2, sinW2
  Real(dp), Parameter ::    e_d = -1._dp / 3._dp,  e_u = -2._dp * e_d   &
    & , oo2pi = 1._dp/(2._dp*pi), oo6pi = oo2pi / 3._dp 

  Integer :: set_mod_par(45)=0, i_par, id
  Logical :: TwoLoopRGE_save, read_data
  Character(len=80) :: read_line
  Character(len=10) :: c1,c2,c3,c4
   Real(dp) :: wert, g2(213), tz, dt, k_fac

  Iname = Iname + 1
  NameOfUnit(Iname) = "FirstGuess"

 read_data = .False.
 If (Present(Take_old_data)) read_data = Take_old_data

 If (read_data) Then ! reading the data from an external file
  If (Present(name)) Then
   Open(99,file=Trim(name),err=200,status="old")
  Else ! use old data from exisiting SPheno.spc
   Open(99,file="Spheno.spc",err=200,status="old")
  End If
   Mi = 0._dp
   M_E2 = 0._dp
   M_L2 = 0._dp
   A_e = 0._dp
   M_D2 = 0._dp
   M_Q2 = 0._dp
   M_U2 = 0._dp
   A_d = 0._dp
   A_u = 0._dp
   mu = 0._dp
   BImu = 0._dp 
   Y_l = 0._dp 
   Y_d = 0._dp
   Y_u = 0._dp

   set_mod_par = 0

   Do ! reading file
    Read(99,"(a80)",End=200,err=200) read_line
! Write(*,*) trim(read_line)
    If (read_line(1:1).Eq."#") Cycle ! ignore comments for the moment
    If (read_line.Eq." ") Cycle ! ignore empty lines for the moment
    Call PutUpperCase(read_line)
    If (read_line(1:5).Eq."BLOCK") Then ! assigning values for the select case
     If (read_line(7:11).Eq."MSOFT") Then
     Backspace(99)
     Read(99,*) c1, c2, c3,mudim,c4
     Do 
      Read(99,*,End=200) read_line
!     Write(*,*) trim(read_line)
      If (read_line(1:1).Eq."#") Cycle ! this loop
      Backspace(99) ! resetting to the beginning of the line
      If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

      Read(99,*) i_par,wert,read_line

     If ((i_par.Ge.1).And.(i_par.Le.3)) Then 
      Mi(i_par) = wert
      set_mod_par(i_par) = 1
     Else If (i_par.Eq.21) Then 
      M2_H(1) = wert
      set_mod_par(4) = 1

     Else If (i_par.Eq.22) Then 
      M2_H(2) = wert
      set_mod_par(5) = 1

     Else If (i_par.Eq.31) Then 
      M2_L(1,1) = wert**2
      set_mod_par(6) = 1
     Else If (i_par.Eq.32) Then 
      M2_L(2,2) = wert**2
      set_mod_par(7) = 1
     Else If (i_par.Eq.33) Then 
      M2_L(3,3) = wert**2
      set_mod_par(8) = 1
     Else If (i_par.Eq.34) Then 
      M2_E(1,1) = wert**2
      set_mod_par(9) = 1
     Else If (i_par.Eq.35) Then 
      M2_E(2,2) = wert**2
      set_mod_par(10) = 1
     Else If (i_par.Eq.36) Then 
      M2_E(3,3) = wert**2
      set_mod_par(11) = 1
     Else If (i_par.Eq.41) Then 
      M2_Q(1,1) = wert**2
      set_mod_par(12) = 1
     Else If (i_par.Eq.42) Then 
      M2_Q(2,2) = wert**2
      set_mod_par(13) = 1
     Else If (i_par.Eq.43) Then 
      M2_Q(3,3) = wert**2
      set_mod_par(14) = 1
     Else If (i_par.Eq.44) Then 
      M2_U(1,1) = wert**2
      set_mod_par(15) = 1
     Else If (i_par.Eq.45) Then 
      M2_U(2,2) = wert**2
      set_mod_par(16) = 1
     Else If (i_par.Eq.46) Then 
      M2_U(3,3) = wert**2
      set_mod_par(17) = 1
     Else If (i_par.Eq.47) Then 
      M2_D(1,1) = wert**2
      set_mod_par(18) = 1
     Else If (i_par.Eq.48) Then 
      M2_D(2,2) = wert**2
      set_mod_par(19) = 1
     Else If (i_par.Eq.49) Then 
      M2_D(3,3) = wert**2
      set_mod_par(20) = 1

       End If
!     Write(errcan,*) i_par,wert
      End Do  ! i_par 

     Else If  (read_line(7:10).Eq."HMIX") Then

     Do 
      Read(99,*,End=200) read_line
!     Write(*,*) trim(read_line)
      If (read_line(1:1).Eq."#") Cycle ! this loop
      Backspace(99) ! resetting to the beginning of the line
      If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

      Read(99,*) i_par,wert,read_line
      If (i_par.Eq.1) Then
      mu = wert
      set_mod_par(21) = 1

     Else If (i_par.Eq.4) Then  ! check with next item
      BImu = wert
      BImu = BImu * tanb /(1._dp+tanb**2)
      set_mod_par(22) = 1
     End If
     End Do

     Else If  (read_line(7:11).Eq."GAUGE") Then
     Do 
      Read(99,*,End=200) read_line
!     Write(*,*) trim(read_line)
      If (read_line(1:1).Eq."#") Cycle ! this loop
      Backspace(99) ! resetting to the beginning of the line
      If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

      Read(99,*) i_par,wert,read_line
      If ((i_par.Ge.1).And.(i_par.Le.3)) Then
       gauge(i_par) = wert
        set_mod_par(i_par+22) = 1
      End If
      End Do

     Else If  (read_line(7:8).Eq."YU") Then
     Do 
      Read(99,*,End=200) read_line
!     Write(*,*) trim(read_line)
      If (read_line(1:1).Eq."#") Cycle ! this loop
      Backspace(99) ! resetting to the beginning of the line
      If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

      Read(99,*) i_par,id,wert,read_line
      If ((i_par.Ge.1).And.(i_par.Le.3)) Then
       Y_u(i_par,id) = wert
        set_mod_par(i_par+25) = 1
      End If
      End Do

     Else If  (read_line(7:8).Eq."YD") Then
     Do 
      Read(99,*,End=200) read_line
!     Write(*,*) trim(read_line)
      If (read_line(1:1).Eq."#") Cycle ! this loop
      Backspace(99) ! resetting to the beginning of the line
      If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

      Read(99,*) i_par,id,wert,read_line
      If ((i_par.Ge.1).And.(i_par.Le.3)) Then
       Y_d(i_par,id) = wert
        set_mod_par(i_par+28) = 1
      End If
      End Do

     Else If  (read_line(7:8).Eq."YE") Then
     Do 
      Read(99,*,End=200) read_line
!     Write(*,*) trim(read_line)
      If (read_line(1:1).Eq."#") Cycle ! this loop
      Backspace(99) ! resetting to the beginning of the line
      If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

      Read(99,*) i_par,id,wert,read_line
      If ((i_par.Ge.1).And.(i_par.Le.3)) Then
       Y_l(i_par,id) = wert
        set_mod_par(i_par+31) = 1
      End If
      End Do

     Else If  (read_line(7:8).Eq."AU") Then
     Do 
      Read(99,*,End=200) read_line
!     Write(*,*) trim(read_line)
      If (read_line(1:1).Eq."#") Cycle ! this loop
      Backspace(99) ! resetting to the beginning of the line
      If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

      Read(99,*) i_par,id,wert,read_line
      If ((i_par.Ge.1).And.(i_par.Le.3)) Then
       A_u(i_par,id) = wert
        set_mod_par(i_par+34) = 1
      End If
      End Do

     Else If  (read_line(7:8).Eq."AD") Then
     Do 
      Read(99,*,End=200) read_line
!     Write(*,*) trim(read_line)
      If (read_line(1:1).Eq."#") Cycle ! this loop
      Backspace(99) ! resetting to the beginning of the line
      If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

      Read(99,*) i_par,id,wert,read_line
      If ((i_par.Ge.1).And.(i_par.Le.3)) Then
       A_d(i_par,id) = wert
        set_mod_par(i_par+37) = 1
      End If
      End Do

     Else If  (read_line(7:8).Eq."AE") Then
     Do 
      Read(99,*,End=200) read_line
!     Write(*,*) trim(read_line)
      If (read_line(1:1).Eq."#") Cycle ! this loop
      Backspace(99) ! resetting to the beginning of the line
      If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

      Read(99,*) i_par,id,wert,read_line
      If ((i_par.Ge.1).And.(i_par.Le.3)) Then
       A_e(i_par,id) = wert
        set_mod_par(i_par+40) = 1
      End If
      End Do

     Else If  (read_line(7:11).Eq."MASS") Then
     Do 
      Read(99,*,End=200) read_line
!     Write(*,*) trim(read_line)
      If (read_line(1:1).Eq."#") Cycle ! this loop
      Backspace(99) ! resetting to the beginning of the line
      If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

      Read(99,*) i_par,wert,read_line

     If (i_par.Eq.24) Then
       mW = wert
       mW2 = mW**2
       set_mod_par(44) = 1
      Else If (i_par.Eq.36) Then
       mP0(2) = wert
       mP02(2) = wert**2
       set_mod_par(45) = 1
      End If
      End Do

     End If ! ignore all other blocks without warning
    End If
   End Do
   ! check if all data have been given

200 If (Sum(set_mod_par).Eq.45) Then
   A_e = A_e * Y_l
   A_d = A_d * Y_d
   A_u = A_u * Y_u
   If (mudim.Gt.mZ) Then
    Call ParametersToG(gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
          & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, BImu, g2)
    tz = Log(mZ/mudim)
    dt = tz / 100._dp
    Call odeint(g2, 213, 0._dp, tz, delta, dt, 0._dp, rge213, kont)
    Call GToParameters(g2,gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
          & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, BImu)
   End If
   vev =  2._dp * mW / gauge(2)
   vevs(1) = vev / Sqrt(1._dp + tanb**2)
   vevs(2) = tanb * vevs(1)
   gU1 = gauge(1)
   gSU2 = gauge(2)
   End If
  Close(99)

 End If

 If ((.Not.read_data).Or.(Sum(set_mod_par).Ne.45) ) Then ! calculate data for the first time
 !---------
 ! W-boson, first rough estimate
 !---------
   mW2 = mZ2 * (0.5_dp + Sqrt(0.25_dp-Alpha_Mz*pi / (sqrt2*G_F*mZ2)))

   mW = Sqrt(mW2)      ! mass
   cosW2 = mw2 / mZ2
   sinW2 = 1._dp - cosW2

   k_fac = 1._dp - alpha * ( oo6pi &
              &  - oo2pi * (57._dp*Log(10._dp)+16._dp*Log(mf_u(3)/mZ))/9._dp)

   gauge(1) = Sqrt( 20._dp*pi*alpha_mZ/(k_fac*3._dp*(1._dp-sinW2)) )
   gauge(2) = Sqrt( 4._dp*pi*alpha_mZ/(k_fac*sinW2))
   k_fac = 1 - AlphaS_mZ * oo2pi * ( 0.5_dp - 4._dp * Log(10._dp) &
         &                         - 2._dp * Log(mf_u(3)/mZ) / 3._dp)
   gauge(3) = Sqrt( 4._dp*pi*alphas_mZ)
   gauge(3) = Sqrt( 4._dp*pi*alphas_mZ / k_fac)

   vev =  2._dp * mW / gauge(2)
   vevs(1) = vev / Sqrt(1._dp + tanb**2)
   vevs(2) = tanb * vevs(1)

   Y_l = 0._dp
   Y_d = 0._dp
   Y_u = 0._dp
   Do i1=1,3
    y_l(i1,i1) = sqrt2 * mf_L_mZ(i1) / vevS(1)
    If (i1.Eq.3) Then ! top and bottom are special
     y_u(i1,i1) = sqrt2 * mf_U(i1)  / vevS(2) &
               & * (1._dp - oo3pi *alphas_mZ *(5._dp +3._dp*Log(mZ2/mf_u2(3))))
     y_d(i1,i1) = sqrt2 * mf_D_mZ(i1) / ( vevS(1) * (1._dp+0.015*tanb*phase_mu))
    Else
     y_u(i1,i1) = sqrt2 * mf_U_mZ(i1) / vevS(2)
     y_d(i1,i1) = sqrt2 * mf_D_mZ(i1) / vevS(1)
    End If
   End Do

   If (GenerationMixing) Then
    If (YukawaScheme.Eq.1) Then
     Y_u = Matmul(Transpose(CKM),Y_u)
     Y_u = Transpose(Y_u)
    Else
     Y_d = Matmul(Conjg(CKM),Y_d)
     Y_d = Transpose(Y_d)
    End If
   End If
   
   Call  CouplingsToG(gauge, y_l, y_d, y_u, g1)
   TwoLoopRGE_save = TwoLoopRGE

   If (.Not.UseFixedScale) Then
    mudim =  0.5_dp*Abs(M2_U_0(3,3)+M2_Q_0(3,3))+4._dp*Abs(Mi_0(3))**2
    mudim = Max(mf_u2(3),mudim) 
    Call SetRGEScale(mudim)
    UseFixedScale = .False.
   Else
    mudim = GetRenormalizationScale() ! from LoopFunctions
   End If

   TwoLoopRGE = .False.

   kont = 0
   If (UseFixedGUTScale) Then
    Call RunRGE(kont, 0.001_dp, vevS, g1, g0, mGUT)

   Else 
    GUT_scale = mZ * Sqrt( Exp( 20._dp * Pi2  &
              &          * (1._dp/gauge(1)**2 - 1._dp/gauge(2)**2) / 7._dp) )
    UseFixedGUTScale = .True.
    Call RunRGE(kont, 0.001_dp, vevS, g1, g0, mGUT)
    UseFixedGUTScale = .False.
   End If

   TwoLoopRGE = TwoLoopRGE_save

   If (kont.Ne.0) Then
    Write(*,*) "Initialization failed, please send the input files used to"
    Write(*,*) "porod@physik.uni-wuerzburg.de so that the problem can be analyzed"
    Write(*,*) "kont",kont
    Iname = Iname - 1
    Return
   End If

   Call GToParameters(g0, gauge, Y_l, Y_d, Y_u, Mi, A_e, A_d, A_u &
                  & , M_E2, M_L2, M_D2, M_Q2, M_U2, M_H2, mu, BiMu)
   Y_u = Transpose(Y_u)
   Y_d = Transpose(Y_d)
   Y_l = Transpose(Y_l)
   A_u = Transpose(A_u)
   A_d = Transpose(A_d)
   A_e = Transpose(A_e)

   cosb2 = 1._dp / (1._dp + tanb**2)
   sinb2 = cosb2 * tanb**2

   abs_mu2 = ( M_H2(2) * sinb2 - M_H2(1) *cosb2 ) /(cosb2-sinb2) - 0.5_dp * mZ2
   If (abs_mu2.Lt.0._dp) Then
    Write (ErrCan,*) 'Warning, in subroutine FirstGuess abs(mu)^2'
    Write (ErrCan,*) 'is smaller 0 :',abs_mu2
    Write (ErrCan,*) 'Setting it to 10^4.'
    Write(ErrCan,*) "Y_t, Y_b, Y_tau, tanb"
    Write(ErrCan,*) y_u(3,3), Y_d(3,3), y_l(3,3), vevS(2)/vevS(1),tanb
    Write(ErrCan,*) "m^2_H"
    Write(ErrCan,*) M_H2
    abs_mu2 = 1.e4_dp
   End If

   abs_mu = Sqrt(abs_mu2)
   mu = abs_mu * phase_mu
   mP02(2) = M_H2(2) + M_H2(1) + 2._dp * abs_mu2
   If (mP02(2).Lt.0._dp) Then
    Write (ErrCan,*) 'Warning, in subroutine FirstGuess'
    Write (ErrCan,*) 'mP02(2) is smaller 0 :',mP02(2)
    Write (ErrCan,*) 'Setting it to its modulus'
    mP02(2) = Abs(mP02(2))
   End If
   mP0(2) = Sqrt(mP02(2))
   Bimu = mP02(2) * tanb / (1 + tanb*tanb)

   gU1 =  Sqrt(3._dp / 5._dp ) * gauge(1)    
   gSU2 = gauge(2)    

  End If ! Take_old_data

  Iname = Iname - 1

 Contains
  Real(dp) Function my_a0(m2)
  Implicit None
   Real(dp), Intent(in) :: m2
   If (m2.Le.0._dp) Then
    my_a0 = 0._dp
   Else
    my_a0 = m2 * (1._dp - Log(m2/mudim)) 
   End If
  End  Function my_a0

  Subroutine PutUpperCase(name)
  Implicit None
   Character(len=80), Intent(inout) :: name
   Integer :: len=80, i1
   Do i1=1,len
    If (name(i1:i1).Eq."a") name(i1:i1) = "A"
    If (name(i1:i1).Eq."b") name(i1:i1) = "B"
    If (name(i1:i1).Eq."c") name(i1:i1) = "C"
    If (name(i1:i1).Eq."d") name(i1:i1) = "D"
    If (name(i1:i1).Eq."e") name(i1:i1) = "E"
    If (name(i1:i1).Eq."f") name(i1:i1) = "F"
    If (name(i1:i1).Eq."g") name(i1:i1) = "G"
    If (name(i1:i1).Eq."h") name(i1:i1) = "H"
    If (name(i1:i1).Eq."i") name(i1:i1) = "I"
    If (name(i1:i1).Eq."j") name(i1:i1) = "J"
    If (name(i1:i1).Eq."k") name(i1:i1) = "K"
    If (name(i1:i1).Eq."l") name(i1:i1) = "L"
    If (name(i1:i1).Eq."m") name(i1:i1) = "M"
    If (name(i1:i1).Eq."n") name(i1:i1) = "N"
    If (name(i1:i1).Eq."o") name(i1:i1) = "O"
    If (name(i1:i1).Eq."p") name(i1:i1) = "P"
    If (name(i1:i1).Eq."q") name(i1:i1) = "Q"
    If (name(i1:i1).Eq."r") name(i1:i1) = "R"
    If (name(i1:i1).Eq."s") name(i1:i1) = "S"
    If (name(i1:i1).Eq."t") name(i1:i1) = "T"
    If (name(i1:i1).Eq."u") name(i1:i1) = "U"
    If (name(i1:i1).Eq."v") name(i1:i1) = "V"
    If (name(i1:i1).Eq."w") name(i1:i1) = "W"
    If (name(i1:i1).Eq."x") name(i1:i1) = "X"
    If (name(i1:i1).Eq."y") name(i1:i1) = "Y"
    If (name(i1:i1).Eq."z") name(i1:i1) = "Z"
   End Do
  End Subroutine PutUpperCase


 End Subroutine FirstGuess
!----------------------
! by Suchita Kulkarni
!----------------------
  Subroutine fv_boundary(YuGUT, YdGUT, YeGUT, M02, A0 &
           & , MSQ2IN, MSU2IN, MSD2IN, TUIN, TDIN, TEIN)
  Implicit None
  
  Real(dp), Intent(in) :: YuGUT(3), YdGUT(3), YeGUT(3), M02
  Complex(dp), Intent(in) :: A0
  Complex(dp), Intent(out), Dimension(3,3) ::  MSQ2IN, MSU2IN, MSD2IN &
      &      , TUIN, TDIN, TEIN
!  Real :: A0, m0  ! I assume they are global variables
    
  Integer :: i,j
    
! initiate the trilinears and soft terms to zero
   TUIN = ZeroC
   TDIN = ZeroC
   TEIN = ZeroC
   MSQ2IN = ZeroC
   MSU2IN = ZeroC
   MSD2IN = ZeroC

! set the soft mass terms  
  Do  i = 1,3
!   Do j = 1,3
     MSQ2IN(i,i) = (1 - 0.9/(YuGUT(3)*YuGUT(3))*(YuGUT(i)*YuGUT(i)))*m02
     MSU2IN(i,i) = (1 - 0.9/(YdGUT(3)*YdGUT(3))*(YuGUT(i)*YuGUT(i)))*m02
     MSD2IN(i,i) = (1 - 0.9/(YdGUT(3)*YdGUT(3))*(YdGUT(i)*YdGUT(i)))*m02
!     MSQ2IN(i,j) = (1 - 0.9/(YuGUT(3)*YuGUT(3))*(YuGUT(i)*YuGUT(i)))*m02
!     MSU2IN(i,j) = (1 - 0.9/(YdGUT(3)*YdGUT(3))*(YuGUT(i)*YuGUT(i)))*m02
!     MSD2IN(i,j) = (1 - 0.9/(YdGUT(3)*YdGUT(3))*(YdGUT(i)*YdGUT(i)))*m02
!   End Do
  End Do

 
! set the diagonal parts of trilinears
  Do  i = 1,3
   TUIN(i,i) = A0*YuGUT(i)
   TDIN(i,i) = A0*YdGUT(i)
   TEIN(i,i) = A0*YeGUT(i)
  End Do

 End Subroutine fv_boundary
!----------------------
! by Suchita Kulkarni
!----------------------
 

 Integer Function GetYukawaScheme()
 !-----------------------------------------------------------------------
 ! Sets the parameter YukawaScheme, which controls wheter the top (=1) or the
 ! down (=2) Yukawa couplings stay diagonal at the low scale 
 ! written by Werner Porod, 20.11.01
 !-----------------------------------------------------------------------
 Implicit None

  GetYukawaScheme = YukawaScheme

 End Function GetYukawaScheme

 Subroutine RunRGE(kont, delta, vevSM, g1, g2, mGUT)
 !-----------------------------------------------------------------------
 ! Uses Runge-Kutta method to integrate RGE's from M_Z to M_GUT
 ! and back, putting in correct thresholds. For the first iteration
 ! only the first 6 couplings are included and a common threshold is used.
 ! Written by Werner Porod, 10.07.99
 ! 07.03.2001: including right handed neutrinos
 ! 24.09.01: portation to f90
 ! 27.03.02: including a check if perturbation theory is valid
 ! 16.09.02: the electroweak scale is now set entirely outside, either 
 !           in the main program or in the routine sugra(...)
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: delta, vevSM(2)
  Real(dp), Intent(inout) :: g1(57)
  Real(dp), Intent(out) :: g2(213), mGUT
  
  Integer:: i1, i2, SumI
  Real(dp) :: g1a(93), g2a(285), g5_a(59), g5_b(180), g1b(75), g2b(267) &
      & , g1c(79), g2c(277), g1d(118), g2d(356), g1f(111), g2f(573)     &
      & , g1g(115), g2g(365)
  Real(dp) :: tz, dt, t_out
  Real(dp) :: mudim, gGUT, g1_h(57), m_hi, m_lo, M15(3)
  Logical :: FoundUnification

  Real(dp), Parameter :: Umns(3,3) = Reshape(   Source = (/  &
      &    Sqrt2/Sqrt3, -ooSqrt2*ooSqrt3, -ooSqrt2*ooSqrt3   &
      &  , ooSqrt3,      ooSqrt3,          ooSqrt3           &
      &  , 0._dp ,       ooSqrt2,         -ooSqrt2 /), shape = (/3, 3/) )
  Real(dp), Parameter :: ZeroR2(2) = 0._dp
  Complex(dp), Dimension(3,3) :: mat3, UnuR, Ynu, Anu, Mr2

  Complex(dp), Dimension(3,3) :: RotXl, RotXr, RotGl, RotGr, RotBl, RotBr &
      & , RotWl, RotWr, ZER, ZEL
  Real(dp) :: EigMWM3(3),EigMBM3(3),EigMXM3(3), EigMGM3(3), test2(2)
  Integer :: ierr

  Iname = Iname + 1
  NameOfUnit(Iname) = 'runRGE'

  !-------------------------------------
  ! running to the high scale
  !-------------------------------------
  g1_h = g1
  If ((HighScaleModel(1:9).Eq.'SUGRA_NuR').Or.(HighScaleModel.Eq.'SUGRA_SU5')) &
  Then
   m_lo = MnuR(1)
   FoundUnification = .False.
  Else If (HighScaleModel.Eq.'SEESAW_II') Then
   m_lo = M_H3(1)
   FoundUnification = .False.
  Else If (HighScaleModel.Eq.'GMSB') Then
   m_lo = MlambdaS
   FoundUnification = .True.

! Florian Staub Seesaw II+III
# ifdef SEESAWIII
  Else If (HighScaleModel.Eq.'SEESAW_III_3G') Then
   m_lo = Abs(MWM30(1,1))
   FoundUnification = .False.
  Else If (HighScaleModel.Eq.'SEESAW_II_SARAH') Then
   m_lo = Abs(MTM_GUT)
   FoundUnification = .False.
# endif SEESAWIII
! Florian Staub Seesaw II+III

  Else If (UseFixedGUTScale) Then ! GUT scale is fixed
   m_lo = GUT_scale
   mGUT = GUT_scale
   FoundUnification = .True.
  Else ! Sugra, strings with minimal particle content
   m_lo = 5.e14_dp
   FoundUnification = .False.
  End If

  tz = Log(mZ/m_lo)
  dt = - tz / 50._dp

  Call odeint(g1, 57, tz, 0._dp, delta, dt, 0._dp, rge57, kont)
  If (kont.Ne.0) Then
   Iname = Iname - 1
   Return
  End If
  !--------------------------
  ! check for perturbativity
  !--------------------------
  If ( (oo4pi*Maxval(g1**2)).Gt.0.5_dp) Then
   Write(ErrCan,*) "Non perturbative regime at high scale"
   If (ErrorLevel.Ge.2) Call TerminateProgram
!   Do i1=1,57
!    If (Abs(g1(i1)).Gt.1.e-12) Write(errcan,*)  i1,g1(i1),oo4pi*g1(i1)**2
!   end do
   Write(errcan,*) " "
   kont = -408
   Call AddError(408)
   Iname = Iname - 1
   Return
  End If
    
  !---------------------------
  ! looking for the GUT scale
  !---------------------------
  If (   (HighScaleModel.Eq.'SUGRA').Or.(HighScaleModel.Eq.'Oscar')         &
   & .Or.(HighScaleModel.Eq.'Str_A').Or.(HighScaleModel.Eq.'Str_B')         &
   & .Or.(HighScaleModel.Eq.'Str_C').Or.(HighScaleModel.Eq.'AMSB')  ) Then 

   If (.Not.UseFixedGUTScale) Then
    tz = Log(m_lo/1.e18_dp)
    dt = - tz / 50._dp

    Call odeintB(g1, 57, tz, 0._dp, delta, dt, 0._dp, rge57, t_out, kont)
    If (kont.Eq.0) Then
     FoundUnification = .True.
     mGUT = 1.e18_dp * Exp(t_out)
     gGUT = Sqrt( 0.5_dp * (g1(1)**2+g1(2)**2) )
     g1(1) = gGUT
     g1(2) = gGUT
     If (StrictUnification) g1(3) = gGUT
    Else
     Write(ErrCan,*) "kont",kont,delta,tz,dt
     Write(ErrCan,*) "m_t",mf_u(3)
     Write (ErrCan,*) "t_out",t_out,1.e18_dp * Exp(t_out)
     Do i1=1,57
      If ((g1(i1).Ne.0._dp).Or.(g1_h(i1).Ne.0._dp)) &
                 & Write(ErrCan,*) i1,g1_h(i1),g1(i1)
     End Do
     Write(ErrCan,*) " " 
     Iname = Iname - 1
     Return
    End If
   End If ! .not.UseFixedGUTScale

  Else If ((HighScaleModel.Eq.'SUGRA_NuR').Or.(HighScaleModel.Eq.'SUGRA_SU5')) &
   &  Then

   Call GToCouplings(g1,gauge_mR,Y_l_mR(1,:,:),Y_d_mR(1,:,:),Y_u_mR(1,:,:))

   Call CouplingsToG3(gauge_mR, Y_l_mR(1,:,:), Y_nu_mR(1,:,:), Y_d_mR(1,:,:) &
          & , Y_u_mR(1,:,:), MnuL5, g1a)
   !---------------------------
   ! running from m_R1 -> m_R2
   !---------------------------
   m_lo = MnuR(1)
   m_hi = MnuR(2)
   If (MnuR(1).Ne.MnuR(2)) Then
    tz = Log(m_lo / m_hi)
    dt = - tz / 50._dp  
    Call odeint(g1a, 93, tz, 0._dp, delta, dt, 0._dp, rge93, kont)
   End If

   Call GToCouplings3(g1a, gauge_mR, Y_l_mR(1,:,:), Y_nu_mR(1,:,:) &
           & , Y_d_mR(1,:,:), Y_u_mR(1,:,:), MnuL5 )

   Call CouplingsToG3(gauge_mR, Y_l_mR(1,:,:), Y_nu_mR(2,:,:), Y_d_mR(1,:,:) &
           & , Y_u_mR(1,:,:), MnuL5, g1a)
   m_lo = m_hi
   !---------------------------
   ! running from m_R2 -> m_R3
   !---------------------------
   m_hi = MnuR(3)
   If (m_lo.Ne.MnuR(3)) Then
    tz = Log(m_lo / m_hi)
    dt = - tz / 50._dp  

    Call odeint(g1a, 93, tz, 0._dp, delta, dt, 0._dp, rge93, kont)
   End If

   Call GToCouplings3(g1a, gauge_mR, Y_l_mR(1,:,:), Y_nu_mR(2,:,:) &
           & , Y_d_mR(1,:,:), Y_u_mR(1,:,:), MnuL5 )

   If (Ynu_at_MR3) Then
    Call FermionMass(Y_l_mR(2,:,:),1._dp,EigMGM3,ZER,ZEL,kont)
    Y_l_mR(3,:,:) = Matmul( Matmul(Conjg(ZER), Y_l_mR(1,:,:)) &
                  &       , Transpose(Conjg(ZEL)))
    Y_nu_mR(3,:,:) = Y_nu_0
   Else 
    Y_l_mR(3,:,:) = Y_l_mR(1,:,:)
   End If

   Call CouplingsToG3(gauge_mR, Y_l_mR(3,:,:), Y_nu_mR(3,:,:), Y_d_mR(1,:,:) &
          & , Y_u_mR(1,:,:), MnuL5, g1a)
   m_lo = m_hi
   !---------------------------
   ! running from m_R3 -> m_GUT
   !---------------------------
   If (UseFixedGUTScale) Then
    tz = Log(MnuR(3)/GUT_scale)
    mGUT = GUT_scale
    dt = - tz / 50._dp
    Call odeint(g1a, 93, tz, 0._dp, delta, dt, 0._dp, rge93, kont)
    If (kont.Ne.0) Then
     Iname = Iname -1
     Return
    End If

   Else
    If (g1a(1).Lt.g1a(2)) Then ! I am still below GUT scale
     tz = Log(MnuR(3)/1.e18_dp)
     dt = - tz / 50._dp
     Call odeintB(g1a, 93, tz, 0._dp, delta, dt, 0._dp, rge93, t_out, kont) 
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e18_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1a(1)**2+g1a(2)**2) )
      g1a(1) = gGUT
      g1a(2) = gGUT
      If (StrictUnification) g1a(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If

    Else If (g1a(1).Eq.g1a(2)) Then ! I am at the GUT scale, very unlikely
                                    ! but possible
     FoundUnification = .True.
     mGUT = 1.e15_dp * Exp(t_out)
     gGUT = g1a(1)
     If (StrictUnification) g1a(3) = gGUT

    Else ! I have already crossed the GUT scale
     tz = Log(MnuR(3)/1.e15_dp)
     dt = - tz / 50._dp
     Call odeintC(g1a, 93, tz, 0._dp, delta, dt, 0._dp, rge93, t_out, kont)
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1a(1)**2+g1a(2)**2) )
      g1a(1) = gGUT
      g1a(2) = gGUT
      If (StrictUnification) g1a(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If
    End If

   End If


   !----------------------------------------------------
   ! run up to SO(10) scale in case of the SU(5) model
   !----------------------------------------------------
   If (HighScaleModel.Eq.'SUGRA_SU5') Then
    tz = Log(mGUT/M_SO_10)
    dt = - tz / 50._dp

    g5_a = 0._dp
    g5_a(1) = g1a(1)
    g5_a(2:19) = g1a(58:75)    ! u-quark Yukawa couplines

    If (g1a(56).Gt.g1a(20)) Then
     g5_a(20:37) = g1a(40:57)   ! d-quark Yukawa couplines
    Else
     g5_a(20:37) = g1a(4:21)   ! lepton Yukawa couplines
    End If
    g5_a(38:55) = g1a(22:39)   ! neutrino Yukawa couplines
!    Write(*,*) "SU(5) Y_tau , Y_b, %",g1a(20),g1a(56),(g1a(20)-g1a(56))/g1a(20)
    g5_a(56) = Real(lam_0, dp)
    g5_a(57) = Aimag(lam_0)
    g5_a(58) = Real(lamp_0, dp)
    g5_a(59) = Aimag(lamp_0)

    Call odeint(g5_a, 59, tz, 0._dp, delta, dt, 0._dp, rge_SU5, kont)
   End If

  Else If (HighScaleModel.Eq.'SUGRA_NuR1') Then

    Call GToCouplings(g1,gauge_mR,Y_l_mR(1,:,:),Y_d_mR(1,:,:),Y_u_mR(1,:,:))


    If (Ynu_at_MR3) Then
     Call FermionMass(Y_l_mR(2,:,:),1._dp,EigMGM3,ZER,ZEL,kont)
     Y_l_mR(1,:,:) = Matmul( Matmul(Conjg(ZER), Y_l_mR(1,:,:)) &
                   &       , Transpose(Conjg(ZEL)))
     Y_nu_mR(1,:,:) = Y_nu_0
    End If

    Call CouplingsToG2(gauge_mR, Y_l_mR(1,:,:), Y_nu_mR(1,:,:), Y_d_mR(1,:,:) &
             & , Y_u_mR(1,:,:), g1a)

    tz = Log(mGUT/GUT_scale)
    dt = tz / 50._dp  

   If (UseFixedGUTScale) Then
    tz = Log(mGUT/GUT_scale)
    mGUT = GUT_scale
    dt = - tz / 50._dp
    Call odeint(g1a, 75, tz, 0._dp, delta, dt, 0._dp, rge75, kont)
    If (kont.Ne.0) Then
     Iname = Iname -1
     Return
    End If

   Else
    tz = Log(mGUT/1.e18_dp)
    dt = - tz / 50._dp
    Call odeintB(g1a, 75, tz, 0._dp, delta, dt, 0._dp, rge75, t_out, kont)
    If (kont.Eq.0) Then
     FoundUnification = .True.
     mGUT = 1.e18_dp * Exp(t_out)
     gGUT = Sqrt( 0.5_dp * (g1a(1)**2+g1a(2)**2) )
     g1a(1) = gGUT
     g1a(2) = gGUT
     If (StrictUnification) g1a(3) = gGUT
    Else
     Iname = Iname - 1
     Return
    End If
   End If

  Else If ((HighScaleModel.Eq.'SEESAW_II').And.Fifteen_plet) Then

   If (g1(1).Gt.g1(2)) Then ! already above the GUT scale
    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
    Call Set_Decoupling_Heavy_States(.True.)

    If (UseFixedGUTScale) Then
     tz = Log(m_lo/GUT_scale)
     mGUT = GUT_scale
     dt = - tz / 50._dp
     Call odeint(g1, 57, tz, 0._dp, delta, dt, 0._dp, rge57, kont)
     If (kont.Ne.0) Then
      Iname = Iname -1
      Return
     End If

    Else ! .not.UseFixedGUTScale
     tz = Log(m_lo/1.e15_dp)
     dt = - tz / 50._dp
     Call odeintC(g1, 57, tz, 0._dp, delta, dt, 0._dp, rge57, t_out, kont)
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1(1)**2+g1(2)**2) )
      g1(1) = gGUT
      g1(2) = gGUT
      If (StrictUnification) g1(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If

    End If ! UseFixedGUTScale

    Call GToCouplings(g1, gauge_mH3, Y_l_mH3, Y_d_mH3, Y_u_mH3)

    M15 = M_H3(1)
    Call CouplingsToG5(gauge_mH3, Y_l_mH3, Y_T_0, Y_d_mH3, Y_u_mH3 &
                 & , Y_T_0, Y_T_0, lam12_0(1), lam12_0(2), M15, g1d )

   Else ! still below the GUT scale

    !---------------------------------------------------
    ! heavy states are to be included
    !---------------------------------------------------
    Call Set_Decoupling_Heavy_States(.False.)

    Call GToCouplings(g1, gauge_mH3, Y_l_mH3, Y_d_mH3, Y_u_mH3)

    M15 = M_H3(1)

    Delta_b_1 = 7._dp
    If (TwoLoopRGE) Then
     Delta_b_2(1,1) = 181._dp/15._dp
     Delta_b_2(1,2) = 29.4_dp 
     Delta_b_2(1,3) = 656._dp/15._dp
     Delta_b_2(2,1) = 9.8_dp
     Delta_b_2(2,2) = 69._dp
     Delta_b_2(2,3) = 16._dp
     Delta_b_2(3,1) = 82._dp/15._dp
     Delta_b_2(3,2) = 6._dp
     Delta_b_2(3,3) = 358._dp/3._dp

    !-----------------------------------------------------
    ! adding shifts to gauge couplings
    !-----------------------------------------------------
     gauge_mH3(1) = gauge_mH3(1) * (1._dp - oo16pi2 * gauge_mH3(1)**2          &
                 &                       * (8._dp/3._dp*Log(MS15_mH3/MT15_mH3) &
                 &                         + Log(MZ15_mH3/MT15_mH3) /6._dp ) )
     gauge_mH3(2) = gauge_mH3(2) * (1._dp - oo16pi2 * gauge_mH3(2)**2          &
                 &                       * 1.5_dp *Log(MZ15_mH3/MT15_mH3) )
     gauge_mH3(3) = gauge_mH3(3) * (1._dp - oo16pi2 * gauge_mH3(3)**2          &
                 &                       * (2.5_dp*Log(MS15_mH3/MT15_mH3) &
                 &                         + Log(MZ15_mH3/MT15_mH3) ) )
    End If


    If (Ynu_at_MR3) Then
     Call FermionMass(Y_l_mH3,1._dp,EigMGM3,ZER,ZEL,kont)
     Y_l_mH3 = Matmul( Matmul(Conjg(ZER), Y_l_mH3), Transpose(Conjg(ZEL)))
     Call CouplingsToG5(gauge_mH3, Y_l_mH3, Y_T_0, Y_d_mH3, Y_u_mH3 &
                 & , Y_T_0, Y_T_0, lam12_MH3(1), lam12_MH3(2), M15, g1d )
    Else 
     Call CouplingsToG5(gauge_mH3, Y_l_mH3, Y_T_mH3, Y_d_mH3, Y_u_mH3 &
                 & , Y_Z_mH3, Y_S_mH3, lam12_MH3(1), lam12_MH3(2), M15, g1d )
    End If

    If (UseFixedGUTScale) Then
     tz = Log(m_lo/GUT_scale)
     mGUT = GUT_scale
     dt = - tz / 50._dp
     Call odeint(g1d, 118, tz, 0._dp, delta, dt, 0._dp, rge118, kont)
     If (kont.Ne.0) Then
      Iname = Iname -1
      Return
     End If

    Else

     If (g1d(1).Lt.g1d(2)) Then ! I am still below GUT scale
      tz = Log(m_lo/1.e18_dp)
      dt = - tz / 50._dp
      Call odeintB(g1d, 118, tz, 0._dp, delta, dt, 0._dp, rge118, t_out, kont) 

      If (kont.Eq.0) Then
       FoundUnification = .True.
       mGUT = 1.e18_dp * Exp(t_out)
       gGUT = Sqrt( 0.5_dp * (g1d(1)**2+g1d(2)**2) )
       g1d(1) = gGUT
       g1d(2) = gGUT
       If (StrictUnification) g1d(3) = gGUT
      Else
       Iname = Iname - 1
       Return
      End If

     Else If (g1d(1).Eq.g1d(2)) Then ! I am at the GUT scale, very unlikely
                                    ! but possible
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = g1d(1)
      If (StrictUnification) g1d(3) = gGUT

     End If
 
    End If

   End If ! check if below or above GUT scale

  Else If (HighScaleModel.Eq.'SEESAW_II') Then

    Call GToCouplings(g1, gauge_mH3, Y_l_mH3, Y_d_mH3, Y_u_mH3)

    If (Ynu_at_MR3) Then
     Call FermionMass(Y_l_mH3,1._dp,EigMGM3,ZER,ZEL,kont)
     Y_l_mH3 = Matmul( Matmul(Conjg(ZER), Y_l_mH3), Transpose(Conjg(ZEL)))
     Y_T_mH3 = Y_T_0
    End If

    Call CouplingsToG4(gauge_mH3, Y_l_mH3, Y_T_mH3, Y_d_mH3, Y_u_mH3 &
            & , lam12_MH3(1), lam12_MH3(2), g1c)

    Delta_b_1(1) = 3.6_dp
    Delta_b_1(2) = 4._dp

    Delta_b_2(1,1) = 8.64_dp
    Delta_b_2(1,2) = 28.8_dp 
    Delta_b_2(2,1) = 9.6_dp
    Delta_b_2(2,2) = 48._dp

   If (UseFixedGUTScale) Then
    tz = Log(m_lo/GUT_scale)
    mGUT = GUT_scale
    dt = - tz / 50._dp
    Call odeint(g1c, 79, tz, 0._dp, delta, dt, 0._dp, rge79, kont)
    If (kont.Ne.0) Then
     Iname = Iname -1
     Return
    End If

   Else

    If (g1c(1).Lt.g1c(2)) Then ! I am still below GUT scale
     tz = Log(m_lo/1.e18_dp)
     dt = - tz / 50._dp
     Call odeintB(g1c, 79, tz, 0._dp, delta, dt, 0._dp, rge79, t_out, kont) 

     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e18_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1c(1)**2+g1c(2)**2) )
      g1c(1) = gGUT
      g1c(2) = gGUT
      If (StrictUnification) g1c(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If

    Else If (g1c(1).Eq.g1c(2)) Then ! I am at the GUT scale, very unlikely
                                    ! but possible
     FoundUnification = .True.
     mGUT = 1.e15_dp * Exp(t_out)
     gGUT = g1c(1)
     If (StrictUnification) g1c(3) = gGUT

    Else ! I have already crossed the GUT scale, therefore I have to use
            ! the MSSM RGEs without Higgs triplets

     delta_b_1 = 0.0_dp
     delta_b_2 = 0.0_dp

     tz = Log(m_lo/1.e15_dp)
     dt = - tz / 50._dp
     Call odeintC(g1, 75, tz, 0._dp, delta, dt, 0._dp, rge75, t_out, kont)
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1(1)**2+g1(2)**2) )
      g1(1) = gGUT
      g1(2) = gGUT
      If (StrictUnification) g1(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If
    End If

   End If


! Florian Staub Seesaw II+III
# ifdef SEESAWIII
  Else If (HighScaleModel.Eq.'SEESAW_III_3G') Then

   If (g1(1).Gt.g1(2)) Then ! already above the GUT scale
    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
    Call Set_Decoupling_Heavy_States(.True.)

    If (UseFixedGUTScale) Then
     tz = Log(m_lo/GUT_scale)
     mGUT = GUT_scale
     dt = - tz / 50._dp
     Call odeint(g1, 57, tz, 0._dp, delta, dt, 0._dp, rge57, kont)
     If (kont.Ne.0) Then
      Iname = Iname -1
      Return
     End If

    Else ! .not.UseFixedGUTScale
     tz = Log(m_lo/1.e15_dp)
     dt = - tz / 50._dp
     Call odeintC(g1, 57, tz, 0._dp, delta, dt, 0._dp, rge57, t_out, kont)
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1(1)**2+g1(2)**2) )
      g1(1) = gGUT
      g1(2) = gGUT
      If (StrictUnification) g1(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If

    End If ! UseFixedGUTScale

    Call GToCouplings(g1, gauge_H24, Ye_H24, Yd_H24, Yu_H24)
    Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24 &
       & ,Ye_H24,Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:),g1f)

    NGHb3 = 3._dp
    NGHg3 = 3._dp 
    NGHw3 = 3._dp 
    NGHx3 = 3._dp 
    NGHxb3 = 3._dp 
    ThresholdCrossed = 3

   Else

    Call GToCouplings(g1, gauge_H24, Ye_H24, Yd_H24, Yu_H24)

   !-----------------------------------------------------
   ! 1 -> 2
   ! adding shifts to gauge couplings
   !-----------------------------------------------------
    If (TwoLoopRGE) Then

     gauge_h24(1) = gauge_h24(1) * (1._dp - oo16pi2 * gauge_h24(1)**2          &
                  &                  * 5._dp/2._dp*Log(MassMXM3(1)/MWM30(1,1)) )
     gauge_h24(2) = gauge_h24(2) * (1._dp - oo16pi2 * gauge_h24(2)**2          &
                &                     *( 1.5_dp *Log(MassMXM3(1)/MWM30(1,1))   &
                &                      + 2._dp *Log(MassMWM3(1)/MWM30(1,1))) )
     gauge_h24(3) = gauge_h24(3) * (1._dp - oo16pi2 * gauge_h24(3)**2          &
                &                       * (Log(MassMXM3(1)/MWM30(1,1)) &
                &                         + 3._dp*Log(MassMGM3(1)/MWM30(1,1)) ))
    End If

    MWM3 = MWM3running(1,:,:)

    NGHb3 = 1._dp
    NGHg3 = 1._dp 
    NGHw3 = 1._dp 
    NGHx3 = 1._dp 
    NGHxb3 = 1._dp 
    ThresholdCrossed = 1 

    Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24 &
      &  ,Ye_H24,Yb30_H24(1,:,:),Yw30_H24(1,:,:),Yx30_H24(1,:,:),g1f)

    m_lo = MWM30(1,1)
    If (MWM30(1,1).Ne.MWM30(2,2)) Then
     m_hi = MWM30(2,2)
     tz = Log(m_lo / m_hi)
     dt = - tz / 50._dp  
     Call odeint(g1f, 111, tz, 0._dp, delta, dt, 0._dp, rge111, kont)

     Call GToParameters111(g1f,gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24 &
       &  ,Yd_H24,Ye_H24,Yb30_H24(1,:,:),Yw30_H24(1,:,:),Yx30_H24(1,:,:))

     m_lo = m_hi
    End If

    Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24 &
       &  ,Ye_H24,Yb30_H24(2,:,:),Yw30_H24(2,:,:),Yx30_H24(2,:,:),g1f)

    !---------
    ! 2 -> 3
    !---------
    If (g1f(1).Gt.g1f(2)) Then ! already above the GUT scale
    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
     Call Set_Decoupling_Heavy_States(.True.)

     tz = Log(m_lo/1.e15_dp)
     dt = - tz / 50._dp
     Call odeintC(g1f, 111, tz, 0._dp, delta, dt, 0._dp, rge111, t_out, kont)
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1f(1)**2+g1f(2)**2) )
      g1f(1) = gGUT
      g1f(2) = gGUT
      If (StrictUnification) g1f(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If

     Call GToParameters111(g1f,gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24   &
       & ,Yd_H24,Ye_H24,Yb30_H24(2,:,:),Yw30_H24(2,:,:),Yx30_H24(2,:,:))
     Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24&
       & ,Ye_H24,Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:),g1f)

     NGHb3 = 3._dp
     NGHg3 = 3._dp 
     NGHw3 = 3._dp 
     NGHx3 = 3._dp 
     NGHxb3 = 3._dp 
     ThresholdCrossed = 3 

    Else   
     MWM3 = MWM3running(2,:,:)
    !-----------------------------------------------------
    ! adding shifts to gauge couplings
    !-----------------------------------------------------
     If (TwoLoopRGE) Then
      gauge_h24(1) = gauge_h24(1) * (1._dp - oo16pi2 * gauge_h24(1)**2         &
                &                    * 5._dp/2._dp*Log(MassMXM3(2)/MWM30(2,2)) )
      gauge_h24(2) = gauge_h24(2) * (1._dp - oo16pi2 * gauge_h24(2)**2         &
                &                     *( 1.5_dp *Log(MassMXM3(2)/MWM30(2,2)) + &
                &                    2._dp *Log(MassMWM3(2)/MWM30(2,2))) )
      gauge_h24(3) = gauge_h24(3) * (1._dp - oo16pi2 * gauge_h24(3)**2         &
                &                       * (Log(MassMXM3(2)/MWM30(2,2)) &
                &                        + 3._dp*Log(MassMGM3(2)/MWM30(2,2)) ) )

     End If

     Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24&
        &  ,Ye_H24,Yb30_H24(2,:,:),Yw30_H24(2,:,:),Yx30_H24(2,:,:),g1f)

     NGHb3 = 2._dp
     NGHg3 = 2._dp 
     NGHw3 = 2._dp 
     NGHx3 = 2._dp 
     NGHxb3 = 2._dp 
     ThresholdCrossed = 2 

     m_lo = MWM30(2,2)
     If (MWM30(2,2).Ne.MWM30(3,3)) Then
      m_hi = MWM30(3,3)
      tz = Log(m_lo / m_hi)
      dt = - tz / 50._dp  
      Call odeint(g1f, 111, tz, 0._dp, delta, dt, 0._dp, rge111, kont)

      Call GToParameters111(g1f,gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24  &
         &  ,Yd_H24,Ye_H24,Yb30_H24(2,:,:),Yw30_H24(2,:,:),Yx30_H24(2,:,:))

      m_lo = m_hi
     End If

     Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24&
         &  ,Ye_H24,Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:),g1f)

     If (g1f(1).Gt.g1f(2)) Then ! already above the GUT scale
    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
      Call Set_Decoupling_Heavy_States(.True.)

      tz = Log(m_lo/1.e15_dp)
      dt = - tz / 50._dp
      Call odeintC(g1f, 111, tz, 0._dp, delta, dt, 0._dp, rge111, t_out, kont)
      If (kont.Eq.0) Then
       FoundUnification = .True.
       mGUT = 1.e15_dp * Exp(t_out)
       gGUT = Sqrt( 0.5_dp * (g1f(1)**2+g1f(2)**2) )
       g1f(1) = gGUT
       g1f(2) = gGUT
       If (StrictUnification) g1f(3) = gGUT
      Else
       Iname = Iname - 1
       Return
      End If

      Call GToParameters111(g1f,gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24  &
        &  ,Yd_H24,Ye_H24,Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:))

      Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24 &
       & ,Yd_H24,Ye_H24,Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:),g1f)


      NGHb3 = 3._dp
      NGHg3 = 3._dp 
      NGHw3 = 3._dp 
      NGHx3 = 3._dp 
      NGHxb3 = 3._dp 
      ThresholdCrossed = 3 

     Else

     !-----------------------------------------------------
     ! 3 -> GUT
     !-----------------------------------------------------

      MWM3 = MWM3running(3,:,:)

      NGHb3 = 3._dp
      NGHg3 = 3._dp 
      NGHw3 = 3._dp 
      NGHx3 = 3._dp 
      NGHxb3 = 3._dp 
      ThresholdCrossed = 3 

     !-----------------------------------------------------
     ! adding shifts to gauge couplings
     !-----------------------------------------------------
      If (TwoLoopRGE) Then
       gauge_h24(1) = gauge_h24(1) * (1._dp - oo16pi2 * gauge_h24(1)**2        &
                &                    * 5._dp/2._dp*Log(MassMXM3(3)/MWM30(3,3)) )
       gauge_h24(2) = gauge_h24(2) * (1._dp - oo16pi2 * gauge_h24(2)**2        &
                &                     *( 1.5_dp *Log(MassMXM3(3)/MWM30(3,3)) + &
                &                     2._dp *Log(MassMWM3(3)/MWM30(3,3))) )
       gauge_h24(3) = gauge_h24(3) * (1._dp - oo16pi2 * gauge_h24(3)**2        &
                &                       * (Log(MassMXM3(3)/MWM30(3,3)) &
                &                        + 3._dp*Log(MassMGM3(3)/MWM30(3,3)) ) )
      
      End If

      If (Ynu_at_MR3) Then
       Call FermionMass(Ye_H24,1._dp,EigMGM3,ZER,ZEL,kont)
       Ye_H24 = Matmul( Matmul(Conjg(ZER), Ye_H24), Transpose(Conjg(ZEL)))

       Yb30_H24(3,:,:) = Yb3_H24_GUT
       Yw30_H24(3,:,:) = Yb3_H24_GUT
       Yx30_H24(3,:,:) = Yb3_H24_GUT
       AYe_H24 = Matmul( Matmul(Conjg(ZER), AYe_H24), Transpose(Conjg(ZEL)))
       ME2_H24 = Matmul( Matmul(ZER, ME2_H24), Transpose(Conjg(ZER)))
       ML2_H24 = Matmul( Matmul(ZEL, ML2_H24), Transpose(Conjg(ZEL)))
      End If

      Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24     &
      & ,Yd_H24,Ye_H24,Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:),g1f)

      If (UseFixedGUTScale) Then
       tz = Log(m_lo/GUT_scale)
       mGUT = GUT_scale
       dt = - tz / 50._dp
       Call odeint(g1f, 111, tz, 0._dp, delta, dt, 0._dp, rge111, kont)
       If (kont.Ne.0) Then
        Iname = Iname -1
        Return
       End If

      Else

       If (g1f(1).Lt.g1f(2)) Then ! I am still below GUT scale
        tz = Log(m_lo/1.e18_dp)
        dt = - tz / 50._dp
        Call odeintB(g1f, 111, tz, 0._dp, delta, dt, 0._dp, rge111, t_out, kont)
        If (kont.Eq.0) Then
         FoundUnification = .True.
         mGUT = 1.e18_dp * Exp(t_out)
         gGUT = Sqrt( 0.5_dp * (g1f(1)**2+g1f(2)**2) )
         g1f(1) = gGUT
         g1f(2) = gGUT
         If (StrictUnification) g1f(3) = gGUT
        Else
         Iname = Iname - 1
         Return
        End If

       Else If (g1f(1).Eq.g1f(2)) Then ! I am at the GUT scale, very unlikely
                                    ! but possible
        FoundUnification = .True.
        mGUT = 1.e15_dp * Exp(t_out)
        gGUT = g1f(1)
        If (StrictUnification) g1f(3) = gGUT

       Else ! I have already crossed the GUT scale
        tz = Log(m_lo/1.e15_dp)
        dt = - tz / 50._dp
        Call odeintC(g1f, 111, tz, 0._dp, delta, dt, 0._dp, rge111, t_out, kont)
        If (kont.Eq.0) Then
         FoundUnification = .True.
         mGUT = 1.e15_dp * Exp(t_out)
         gGUT = Sqrt( 0.5_dp * (g1f(1)**2+g1f(2)**2) )
         g1f(1:2) = gGUT
         If (StrictUnification) g1f(3) = gGUT
        Else
         Iname = Iname - 1
         Return
        End If
       End If

      End If

     End If
    End If


   End If


  Else If (HighScaleModel.Eq.'SEESAW_II_SARAH') Then

   If (g1(1).Gt.g1(2)) Then ! already above the GUT scale
    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
    If (UseFixedGUTScale) Then
     tz = Log(m_lo/GUT_scale)
     mGUT = GUT_scale
     dt = - tz / 50._dp
     Call odeint(g1, 57, tz, 0._dp, delta, dt, 0._dp, rge57, kont)
     If (kont.Ne.0) Then
      Iname = Iname -1
      Return
     End If

    Else ! .not.UseFixedGUTScale
     tz = Log(m_lo/1.e15_dp)
     dt = - tz / 50._dp
     Call odeintC(g1, 57, tz, 0._dp, delta, dt, 0._dp, rge57, t_out, kont)
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1(1)**2+g1(2)**2) )
      g1(1) = gGUT
      g1(2) = gGUT
      If (StrictUnification) g1(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If

    End If ! UseFixedGUTScale

    Call GToCouplings(g1, gauge_h15, Ye_h15, Yd_h15, Yu_h15)
    M15 = M_H3(1)

    Call ParametersToG115(gauge_h15(1),gauge_h15(2),gauge_h15(3), &
    & Yu_h15,Yd_h15,Ye_h15,Yt0_h15,Ys0_h15,Yz0_h15,Lambda10,Lambda20,g1g)

   Else ! still below the GUT scale


    !---------------------------------------------------
    ! heavy states are to be included
    !---------------------------------------------------
    Call GToCouplings(g1, gauge_h15, Ye_h15, Yd_h15, Yu_h15)


    If (TwoLoopRGE) Then
     !-----------------------------------------------------
     ! adding shifts to gauge couplings
     !-----------------------------------------------------
     gauge_h15(1) = gauge_h15(1) * (1._dp - 0.2_dp * oo16pi2 * gauge_h15(1)**2 &
                 &                       * (8._dp * Log(MSM/MTM_GUT)           &
                 &                         + 9._dp * Log(MTM0/MTM_GUT)         &
                 &                         + 0.5_dp*Log(MZM/MTM_GUT) ) )
     gauge_h15(2) = gauge_h15(2) * (1._dp - oo16pi2 * gauge_h15(2)**2          &
                 &    * (2._dp *Log(MTM0/MTM_GUT) + 1.5_dp *Log(MZM/MTM_GUT) ))
     gauge_h15(3) = gauge_h15(3) * (1._dp - oo16pi2 * gauge_h15(3)**2          &
                 &    * (2.5_dp*Log(MSM/MTM_GUT) + Log(MZM/MTM_GUT) ) )
    End If


    If (Ynu_at_MR3) Then
     Call FermionMass(Ye_H15,1._dp,EigMGM3,ZER,ZEL,kont)
     Ye_H15 = Matmul( Matmul(Conjg(ZER), Ye_H15), Transpose(Conjg(ZEL)))
     YT0_H15 = Y_T_0
     Ys0_H15 = Y_T_0
     Yz0_H15 = Y_T_0
    End If

    Call ParametersToG115(gauge_h15(1),gauge_h15(2),gauge_h15(3),Yu_h15,Yd_h15 &
       & ,Ye_h15,Yt0_h15,Ys0_h15,Yz0_h15,Lambda10,Lambda20,g1g)

    If (UseFixedGUTScale) Then
     tz = Log(m_lo/GUT_scale)
     mGUT = GUT_scale
     dt = - tz / 50._dp
     Call odeint(g1g, 115, tz, 0._dp, delta, dt, 0._dp, rge115, kont)
     If (kont.Ne.0) Then
      Iname = Iname -1
      Return
     End If

    Else

     If (g1g(1).Lt.g1g(2)) Then ! I am still below GUT scale
      tz = Log(m_lo/1.e18_dp)
      dt = - tz / 50._dp
      Call odeintB(g1g, 115, tz, 0._dp, delta, dt, 0._dp, rge115, t_out, kont) 

      If (kont.Eq.0) Then
       FoundUnification = .True.
       mGUT = 1.e18_dp * Exp(t_out)
       gGUT = Sqrt( 0.5_dp * (g1g(1)**2+g1g(2)**2) )
       g1g(1) = gGUT
       g1g(2) = gGUT
       If (StrictUnification) g1g(3) = gGUT
      Else
       Iname = Iname - 1
       Return
      End If

     Else If (g1g(1).Eq.g1g(2)) Then ! I am at the GUT scale, very unlikely
                                    ! but possible
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = g1g(1)
      If (StrictUnification) g1g(3) = gGUT

     End If
 
    End If


   End If ! check if below or above GUT scale
# endif SEESAWIII
! Florian Staub Seesaw II+III

  End If

  If ((.Not.UseFixedGUTScale).And.(.Not.FoundUnification)) Then
   Write (ErrCan,*) 'SUGRA: no unification found'
   SugraErrors(1) = .True.
   kont = -409
   Call AddError(409)
   Iname = Iname - 1
   Return
  End If

  !------------------------------------
  ! Saving parameters at high scale
  !------------------------------------
  If (HighScaleModel.Eq."GMSB") Then
   mGUT_Save = MlambdaS
  Else
   mGUT_Save = mGUT
  End If
  !---------------------------------------
  ! boundary condition at the high scale
  !---------------------------------------
  If (HighScaleModel.Eq.'SUGRA_NuR') Then
   Call BoundaryHS(g1a,g2a)

  Else If (HighScaleModel.Eq.'SUGRA_NuR1') Then
   Call BoundaryHS(g1b,g2b)

  Else If ((HighScaleModel.Eq.'SEESAW_II').And.Fifteen_plet) Then
    Call BoundaryHS(g1d,g2d)

  Else If (HighScaleModel.Eq.'SEESAW_II') Then
    Call BoundaryHS(g1c,g2c)

  Else If (HighScaleModel.Eq.'SUGRA_SU5') Then
   g5_b = 0._dp
   g5_b(1:59) = g5_a   ! gauge and Yukawa coupling
   g5_b(38:55) = g5_b(2:19) ! setting Y_nu = Y_u
   g5_b(60:63) = 0._dp ! bilinear parameters
   g5_b(64) = Real( Mi_0(1), dp )
   g5_b(65) = Aimag( Mi_0(1) )
   g5_b(66:123) = AoY_u_0(1,1) * g5_a(2:59)
   Do i1=1,3
    g5_b(8*i1+116) = M2_H_0(1) + D_SO_10          ! 10-plets
    g5_b(8*i1+134) = M2_H_0(1) - 3._dp * D_SO_10  ! 5-plets, matter
    g5_b(8*i1+152) = M2_H_0(1) + 5._dp * D_SO_10  ! singlet, nu_R
   End Do
   g5_b(178) = M2_H_0(1) - 2._dp * D_SO_10  ! 5-plets, u-type Higgs
   g5_b(179) = M2_H_0(1) + 2._dp * D_SO_10  ! 5-plets, d-type Higgs
   g5_b(180) = M2_H_0(1)                    ! Higgs 24-plet

! Florian Staub Seesaw II+III
# ifdef SEESAWIII
 Else If (HighScaleModel.Eq.'SEESAW_III_3G') Then
   Call BoundaryHS(g1f,g2f)

  Else If (HighScaleModel.Eq.'SEESAW_II_SARAH') Then

   Call BoundaryHS(g1g,g2g)
# endif SEESAWIII
! Florian Staub Seesaw II+III
  Else
   Call BoundaryHS(g1,g2)
  End If
  
  !--------------------------------------
  ! running down to the electroweak scale
  !--------------------------------------
  If ((HighScaleModel.Eq.'SUGRA_NuR').Or.(HighScaleModel.Eq.'SUGRA_SU5')) Then
   !-------------------------------------
   ! m_SO(10) -> m_GUT in case of SU(5)
   !-------------------------------------
   If (HighScaleModel.Eq.'SUGRA_SU5') Then
    tz = Log(mGUT/M_SO_10)
    dt = tz / 50._dp
    Call odeint(g5_b, 180, 0._dp, tz, delta, dt, 0._dp, rge_SU5, kont)

    g2a(1:75) = g1a(1:75)  ! using orginal boundary conditions
    g2a(22:39) = g5_b(38:55) ! neutrino Yukawa couplings
    g2a(76:77) = g5_b(64:65)  ! gaugino mass parameters
    g2a(78:79) = g5_b(64:65)
    g2a(80:81) = g5_b(64:65)
    g2a(82:99) = g5_b(84:101)    ! A_l = A_d
    g2a(100:117) = g5_b(102:119) ! A_nu
    g2a(118:135) = g5_b(84:101)  ! A_d = A_l
    g2a(136:153) = g5_b(66:83)   ! A_u
    g2a(154:171) = g5_b(124:141) ! M_E = M_Q = M_U
    g2a(172:189) = g5_b(142:159) ! M_L = M_D
    g2a(190:207) = g5_b(160:177) ! M_R
    g2a(208:225) = g5_b(142:159) ! M_L = M_D
    g2a(226:243) = g5_b(124:141) ! M_Q = M_U = M_E
    g2a(244:261) = g5_b(124:141) ! M_U = M_E = M_Q
    g2a(262:263) = g5_b(178:179) ! M_H
    g2a(264:285) = 0._dp         ! mu, B, MnuL5
    M2S_GUT = g5_b(180)
    Alam_GUT = Cmplx(g5_b(120),g5_b(121),dp)
    Alamp_GUT = Cmplx(g5_b(122),g5_b(123),dp)
    Do i1=1,3
     Do i2=1,3
      SumI = 6*i1+2*i2
      Ynu(i1,i2) = Cmplx( g2a(SumI+14), g2a(SumI+15),dp )
      Anu(i1,i2) = Cmplx( g2a(SumI+92), g2a(SumI+93),dp )
      Mr2(i1,i2) = Cmplx( g2a(SumI+182), g2a(SumI+183),dp )
     End Do
    End Do

   Else

    Ynu = Y_nu_0
    Anu = A_nu_0
    Mr2 = M2_R_0
   End If

    !--------------------------------------
    ! recalculating m_Nu_R, if necessary
    !--------------------------------------
    If ((.Not.Fixed_Nu_Yukawas).And.(.Not. Ynu_at_MR3)) Then
     mat3 = 0._dp
     mat3(1,1) = 1._dp / mf_nu(1)
     mat3(2,2) = 1._dp / mf_nu(2)
     mat3(3,3) = 1._dp / mf_nu(3)

     mat3 = Matmul(Umns,Matmul(mat3,Transpose(Umns)))
     mat3 = vevSM(2)**2 * Matmul(Ynu,Matmul(mat3,Transpose(Ynu)))

     Call Neutrinomasses(mat3, mNuR, UnuR, kont)

!    Else ! for the rotation below
!     UnuR = id3C
    !--------------------------------------------
    ! rotating R-neutrinos to mass eigenbasis
    !--------------------------------------------
     Ynu = Matmul(UnuR,Ynu)
     Anu = Matmul(UnuR,Anu)
     MR2 = Matmul(Transpose(Conjg(UnuR)),Matmul(MR2,UnuR))

     Do i1=1,3
      Do i2=1,3
       SumI = 6*i1+2*i2
       g2a(SumI+14) = Real(Ynu(i1,i2),dp)
       g2a(SumI+15) = Aimag(Ynu(i1,i2))
       g2a(SumI+92) = Real(Anu(i1,i2),dp)
       g2a(SumI+93) = Aimag(Anu(i1,i2))
       g2a(SumI+182) = Real(MR2(i1,i2),dp)
       g2a(SumI+183) = Aimag(MR2(i1,i2))
      End Do
     End Do
    End If

    
   !------------------------
   ! m_GUT -> m_nuR_3
   !------------------------
   m_hi = mGUT
   m_lo = MNuR(3)
   If (Abs(m_lo).Lt.Abs(m_hi)) Then
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)
    m_hi = m_lo
   Endif

   Call GToParameters3(g2a, gauge_mR, y_l_mR(3,:,:), y_nu_mR(3,:,:)          &
      & , y_d_mR(3,:,:), y_u_mR(3,:,:), Mi_mR, A_l_mR(3,:,:), A_nu_mR(3,:,:) &
      & , A_d_mR(3,:,:), A_u_mR(3,:,:), M2_E_mR(3,:,:), M2_L_mR(3,:,:)       &
      & , M2_R_mR(3,:,:), M2_D_mR(3,:,:), M2_Q_mR(3,:,:), M2_U_mR(3,:,:)     &
      & , M2_H_mR, mu_mR, B_mR, MnuL5)

   If (Ynu_at_MR3) Then  ! to enhance numerical stability and speed
    Call FermionMass(Y_l_mR(3,:,:),1._dp,EigMGM3,ZER,ZEL,kont)
    Y_l_mR(3,:,:) = Matmul( Matmul(Conjg(ZER), Y_l_mR(3,:,:)) &
                  &       , Transpose(Conjg(ZEL)))
    Y_nu_mR(3,:,:) = Y_nu_0
    A_l_mR(3,:,:) = Matmul( Matmul(Conjg(ZER), A_l_mR(3,:,:)) &
                  &       , Transpose(Conjg(ZEL)))
    M2_E_mR(3,:,:) = Matmul( Matmul(ZER, M2_E_mR(3,:,:)) &
                   &       , Transpose(Conjg(ZER)))
    M2_L_mR(3,:,:) = Matmul( Matmul(ZEL, M2_L_mR(3,:,:)) &
                   &       , Transpose(Conjg(ZEL)))
   End If

   Do i1=1,3
    Do i2=1,3
     MnuL5(i1,i2) = - Y_nu_mR(3,3,i1) * Y_nu_mR(3,3,i2) / MNuR(3)
    End Do
   End Do
   Y_nu_mR(2,:,:) = Y_nu_mR(3,:,:)
   Y_nu_mR(2,3,:) = 0._dp
   A_nu_mR(2,:,:) = A_nu_mR(3,:,:)
   A_nu_mR(2,3,:) = 0._dp
   M2_R_mR(2,:,:) = M2_R_mR(3,:,:)
   M2_R_mR(2,3,:) = 0._dp
   M2_R_mR(2,:,3) = 0._dp
   !------------------------
   ! m_nuR_3 -> m_nuR_2
   !------------------------
   Call ParametersToG3(gauge_mR, y_l_mR(3,:,:), y_nu_mR(2,:,:), y_d_mR(3,:,:) &
      & , y_u_mR(3,:,:), Mi_mR, A_l_mR(3,:,:), A_nu_mR(2,:,:), A_d_mR(3,:,:)  &
      & , A_u_mR(3,:,:), M2_E_mR(3,:,:), M2_L_mR(3,:,:), M2_R_mR(2,:,:)       &
      & , M2_D_mR(3,:,:), M2_Q_mR(3,:,:), M2_U_mR(3,:,:), M2_H_mR, mu_mR      &
      & , B_mR, MnuL5, g2a)

   m_lo = MNuR(2)
   If (Abs(m_lo).Lt.Abs(m_hi)) Then
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)
    m_hi = m_lo
   Endif

   Call GToParameters3(g2a, gauge_mR, y_l_mR(2,:,:), y_nu_mR(2,:,:)          &
      & , y_d_mR(2,:,:), y_u_mR(2,:,:), Mi_mR, A_l_mR(2,:,:), A_nu_mR(2,:,:) &
      & , A_d_mR(2,:,:), A_u_mR(2,:,:), M2_E_mR(2,:,:), M2_L_mR(2,:,:)       &
      & , M2_R_mR(2,:,:), M2_D_mR(2,:,:), M2_Q_mR(2,:,:), M2_U_mR(2,:,:)     &
      & , M2_H_mR, mu_mR, B_mR, MnuL5)

   Do i1=1,3
    Do i2=1,3
     MnuL5(i1,i2) = MnuL5(i1,i2) - Y_nu_mR(2,2,i1) * Y_nu_mR(2,2,i2) / MNuR(2)
    End Do
   End Do
   Y_nu_mR(1,:,:) = Y_nu_mR(2,:,:)
   Y_nu_mR(1,2,:) = 0._dp
   A_nu_mR(1,:,:) = A_nu_mR(2,:,:)
   A_nu_mR(1,2,:) = 0._dp
   M2_R_mR(1,:,:) = M2_R_mR(2,:,:)
   M2_R_mR(1,2,:) = 0._dp
   M2_R_mR(1,:,2) = 0._dp
   !------------------------
   ! m_nuR_2 -> m_nuR_1
   !------------------------
   Call ParametersToG3(gauge_mR, y_l_mR(2,:,:), y_nu_mR(1,:,:), y_d_mR(2,:,:) &
      & , y_u_mR(2,:,:), Mi_mR, A_l_mR(2,:,:), A_nu_mR(1,:,:), A_d_mR(2,:,:)  &
      & , A_u_mR(2,:,:), M2_E_mR(2,:,:), M2_L_mR(2,:,:), M2_R_mR(1,:,:)       &
      & , M2_D_mR(2,:,:), M2_Q_mR(2,:,:), M2_U_mR(2,:,:), M2_H_mR, mu_mR      &
      & , B_mR, MnuL5, g2a)

   m_lo = MNuR(1)
   If (Abs(m_lo).Lt.Abs(m_hi)) Then
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)
    m_hi = m_lo
   Endif

   Call GToParameters3(g2a, gauge_mR, y_l_mR(1,:,:), y_nu_mR(1,:,:)          &
      & , y_d_mR(1,:,:), y_u_mR(1,:,:), Mi_mR, A_l_mR(1,:,:), A_nu_mR(1,:,:) &
      & , A_d_mR(1,:,:), A_u_mR(1,:,:), M2_E_mR(1,:,:), M2_L_mR(1,:,:)       &
      & , M2_R_mR(1,:,:), M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:)     &
      & , M2_H_mR, mu_mR, B_mR, MnuL5)

   Do i1=1,3
    Do i2=1,3
     MnuL5(i1,i2) = MnuL5(i1,i2) - Y_nu_mR(1,1,i1) * Y_nu_mR(1,1,i2) / MNuR(1)
    End Do
   End Do
   !------------------------
   ! m_nuR_1 -> Q_EWSB
   !------------------------
   Call ParametersToG3(gauge_mR, y_l_mR(1,:,:), Zero33C, y_d_mR(1,:,:)     &
        & , y_u_mR(1,:,:), Mi_mR, A_l_mR(1,:,:), Zero33C, A_d_mR(1,:,:)    &
        & , A_u_mR(1,:,:), M2_E_mR(1,:,:), M2_L_mR(1,:,:), Zero33C         &
        & , M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:), M2_H_mR, mu_mR &
        & , B_mR, MnuL5, g2a)

   mudim = GetRenormalizationScale()
   mudim = Max(mudim, mZ2)
   tz = 0.5_dp * Log(m_hi**2/mudim)
   dt = - tz / 100._dp

   Call odeint(g2a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)

   Call GToParameters3(g2a, gauge, y_l, Y_nu, y_d, y_u, Mi, A_l, A_nu, A_d &
      & , A_u, M2_E, M2_L, M2_R, M2_D, M2_Q, M2_U, M2_H, mu, B, MnuL5)

   Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u, M2_E, M2_L  &
      & , M2_D, M2_Q, M2_U, M2_H, mu, B, g2)

  Else If (HighScaleModel.Eq.'SUGRA_NuR1') Then
   mudim = MNuR(1)**2
   tz = 0.5_dp * Log(mudim/mGUT**2)
   dt = tz / 50._dp
   Call odeint(g2a, 267, 0._dp, tz, delta, dt, 0._dp, rge267, kont)


   Call GToParameters2(g2a, gauge_mR, y_l_mR(1,:,:), y_nu_mR(1,:,:)          &
      & , y_d_mR(1,:,:), y_u_mR(1,:,:), Mi_mR, A_l_mR(1,:,:), A_nu_mR(1,:,:) &
      & , A_d_mR(1,:,:), A_u_mR(1,:,:), M2_E_mR(1,:,:), M2_L_mR(1,:,:)       &
      & , M2_R_mR(1,:,:), M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:)     &
      & , M2_H_mR, mu_mR, B_mR)


   Call ParametersToG(gauge_mR, y_l_mR(1,:,:), y_d_mR(1,:,:), y_u_mR(1,:,:)  &
      & , Mi_mR, A_l_mR(1,:,:), A_d_mR(1,:,:), A_u_mR(1,:,:), M2_E_mR(1,:,:) &
      & , M2_L_mR(1,:,:), M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:)     &
      & , M2_H_mR, mu_mR, B_mR, g2)

   mudim = GetRenormalizationScale()
   mudim = Max(mudim, mZ2)

   tz = 0.5_dp * Log(mudim/mNuR(1)**2)
   dt = tz / 100._dp
   Call odeint(g2, 213, 0._dp, tz, delta, dt, 0._dp, rge213, kont)

  Else If ((HighScaleModel.Eq.'SEESAW_II').And.Fifteen_plet) Then

   If ( (oo4pi*Maxval(g2d(1:115)**2)).Gt.0.5_dp) Then
    Write(ErrCan,*) "Non perturbative regime at M_GUT"
    If (ErrorLevel.Ge.2) Call TerminateProgram
    Write(errcan,*) " "
    kont = -410
    Call AddError(410)
    Iname = Iname - 1
    Return
   End If
   
   !-------------------------
   ! run only if m_H < m_GUT
   !-------------------------
   If (m_H3(1).Lt.mGUT) Then
    !---------------------------------------------------
    ! heavy states are to be included
    !---------------------------------------------------
    Call Set_Decoupling_Heavy_States(.False.)

    mudim = M_H3(1)**2
    tz = 0.5_dp * Log(mudim/mGUT**2)
    dt = tz / 50._dp
    Call odeint(g2d, 356 , 0._dp, tz, delta, dt, 0._dp, rge356, kont)
    m_lo = M_H3(1)
    If ( (oo4pi*Maxval(g2d(1:115)**2)).Gt.0.5_dp) Then
     Write(ErrCan,*) "Non perturbative regime at M_H3"
     If (ErrorLevel.Ge.2) Call TerminateProgram
     Write(errcan,*) " "
     kont = -411
     Call AddError(411)
     Iname = Iname - 1
     Return
    End If

   Else
     m_lo = mGUT
   End If

   Call GToParameters5(g2d, gauge_mH3, y_l_mH3, y_T_mH3, y_d_mH3, y_u_mH3     &
          & , y_Z_mH3, y_S_mH3, lam12_mH3(1), lam12_mH3(2), Mi_mH3, A_l_mH3   &
          & , A_T_mH3, A_d_mH3, A_u_mH3, A_Z_mH3, A_S_mH3, Alam12_MH3(1)      &
          & , Alam12_MH3(2), M2_E_mH3, M2_L_mH3, M2_D_mH3, M2_Q_mH3, M2_U_mH3 &
          & , M2_H_mH3, M2_T_mH3, M2_Z_mH3, M2_S_mH3, MT15_mH3, MZ15_mH3      &
          & , MS15_mH3, mu_mH3, B_mH3, MnuL5)

   MnuL5 = - lam12_MH3(2) * y_T_mH3 / M_H3(1)

   Delta_b_1 = 0._dp ! decoupling the Higgs triplets
   Delta_b_2 = 0._dp ! decoupling the Higgs triplets

   If ((m_H3(1).Lt.mGUT).And.TwoLoopRGE) Then
   !-----------------------------------------------------
   ! adding shifts to gauge couplings
   !-----------------------------------------------------
    gauge_mH3(1) = gauge_mH3(1) / (1._dp - oo16pi2 * gauge_mH3(1)**2           &
                &                       * (8._dp/3._dp*Log(MS15_mH3/MT15_mH3) &
                &                         + Log(MZ15_mH3/MT15_mH3) /6._dp ) )
    gauge_mH3(2) = gauge_mH3(2) / (1._dp - oo16pi2 * gauge_mH3(2)**2           &
                &                       * 1.5_dp *Log(MZ15_mH3/MT15_mH3) )
    gauge_mH3(3) = gauge_mH3(3) / (1._dp - oo16pi2 * gauge_mH3(3)**2           &
                &                       * (2.5_dp*Log(MS15_mH3/MT15_mH3) &
                &                         + Log(MZ15_mH3/MT15_mH3) ) )
    Mi_mH3(1) = Mi_mH3(1) / (1._dp - oo16pi2 * gauge_mH3(1)**2           &
                &                       * (8._dp/3._dp*Log(MS15_mH3/MT15_mH3) &
                &                         + Log(MZ15_mH3/MT15_mH3) /6._dp ) )
    Mi_mH3(2) = Mi_mH3(2) / (1._dp - oo16pi2 * gauge_mH3(2)**2           &
                &                       * 1.5_dp *Log(MZ15_mH3/MT15_mH3) )
    Mi_mH3(3) = Mi_mH3(3) / (1._dp - oo16pi2 * gauge_mH3(3)**2           &
                &                       * (2.5_dp*Log(MS15_mH3/MT15_mH3) &
                &                         + Log(MZ15_mH3/MT15_mH3) ) )
   End If

    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
    Call Set_Decoupling_Heavy_States(.True.)

   Call ParametersToG4(gauge_mH3, y_l_mH3, Zero33C, y_d_mH3, y_u_mH3   &
          & , ZeroC, ZeroC, Mi_mH3, A_l_mH3, Zero33C, A_d_mH3 &
          & , A_u_mH3, ZeroC, ZeroC, M2_E_mH3, M2_L_mH3     &
          & , M2_D_mH3, M2_Q_mH3, M2_U_mH3, M2_H_mH3, ZeroR2, mu_mH3      &
          & , B_mH3, MnuL5, g2c)

   mudim = GetRenormalizationScale()
   mudim = Max(mudim, mZ2)

   tz = 0.5_dp * Log(mudim/m_lo**2)
   dt = tz / 100._dp
   Call odeint(g2c, 277, 0._dp, tz, delta, dt, 0._dp, rge277, kont)

   Call GToParameters4(g2c, gauge, y_l, y_T, y_d, y_u, lam12(1), lam12(2), Mi &
          & , A_l, A_T, A_d, A_u, Alam12(1), Alam12(2), M2_E, M2_L, M2_D      &
          & , M2_Q, M2_U, M2_H, M2_T, mu, B, MnuL5)
   Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u, M2_E, M2_L  &
      & , M2_D, M2_Q, M2_U, M2_H, mu, B, g2)

  Else If (HighScaleModel.Eq.'SEESAW_II') Then

   If ( (oo4pi*Maxval(g2c(1:79)**2)).Gt.0.5_dp) Then
    Write(ErrCan,*) "Non perturbative regime at M_GUT"
    If (ErrorLevel.Ge.2) Call TerminateProgram
    Write(errcan,*) " "
    kont = -410
    Call AddError(410)
    Iname = Iname - 1
    Return
   End If
   
   !-------------------------
   ! run only if m_H < m_GUT
   !-------------------------
   If (m_H3(1).Lt.mGUT) Then
    mudim = M_H3(1)**2
    tz = 0.5_dp * Log(mudim/mGUT**2)
    dt = tz / 50._dp
    Call odeint(g2c, 277, 0._dp, tz, delta, dt, 0._dp, rge277, kont)
    m_lo = M_H3(1)
   Else
     m_lo = mGUT
   End If

   If ( (oo4pi*Maxval(g2c(1:79)**2)).Gt.0.5_dp) Then
    Write(ErrCan,*) "Non perturbative regime at M_H3"
    If (ErrorLevel.Ge.2) Call TerminateProgram
    Write(errcan,*) " "
    kont = -411
    Call AddError(411)
    Iname = Iname - 1
    Return
   End If

   Call GToParameters4(g2c, gauge_mH3, y_l_mH3, y_T_mH3, y_d_mH3, y_u_mH3   &
          & , lam12_mH3(1), lam12_mH3(2), Mi_mH3, A_l_mH3, A_T_mH3, A_d_mH3 &
          & , A_u_mH3, Alam12_MH3(1), Alam12_MH3(2), M2_E_mH3, M2_L_mH3     &
          & , M2_D_mH3, M2_Q_mH3, M2_U_mH3, M2_H_mH3, M2_T_mH3, mu_mH3      &
          & , B_mH3, MnuL5)

   MnuL5 = - lam12_MH3(2) * y_T_mH3 / M_H3(1)

   Delta_b_1 = 0._dp ! decoupling the Higgs triplets
   Delta_b_2 = 0._dp ! decoupling the Higgs triplets
   Call ParametersToG4(gauge_mH3, y_l_mH3, Zero33C, y_d_mH3, y_u_mH3   &
          & , ZeroC, ZeroC, Mi_mH3, A_l_mH3, Zero33C, A_d_mH3 &
          & , A_u_mH3, ZeroC, ZeroC, M2_E_mH3, M2_L_mH3     &
          & , M2_D_mH3, M2_Q_mH3, M2_U_mH3, M2_H_mH3, ZeroR2, mu_mH3      &
          & , B_mH3, MnuL5, g2c)

   mudim = GetRenormalizationScale()
   mudim = Max(mudim, mZ2)

   tz = 0.5_dp * Log(mudim/m_lo**2)
   dt = tz / 100._dp
   Call odeint(g2c, 277, 0._dp, tz, delta, dt, 0._dp, rge277, kont)

   Call GToParameters4(g2c, gauge, y_l, y_T, y_d, y_u, lam12(1), lam12(2), Mi &
          & , A_l, A_T, A_d, A_u, Alam12(1), Alam12(2), M2_E, M2_L, M2_D      &
          & , M2_Q, M2_U, M2_H, M2_T, mu, B, MnuL5)

   Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u, M2_E, M2_L  &
      & , M2_D, M2_Q, M2_U, M2_H, mu, B, g2)

! Florian Staub Seesaw II+III
# ifdef SEESAWIII
  Else If (HighScaleModel.Eq.'SEESAW_III_3G') Then 
   !---------------
   ! GUT -> 3
   !---------------
   m_hi = mGUT
   If (Abs(MWM30(3,3)).Lt.Abs(m_hi)) Then
    m_lo = Abs(MWM30(3,3))
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2f, 573, tz, 0._dp, delta, dt, 0._dp, rge555, kont)
    m_hi = m_lo
   Endif

   Call GToParameters555(g2f,g1_H24,g2_H24,g3_H24,Yu_H24,Yd_H24,Ye_H24,Yb3_H24 &
      & ,Yw3_H24,Yx3_H24,mu_H24,MXM3,MWM3,MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24    &
      & ,AYb3_H24,AYw3_H24,AYx3_H24,Amue,AMXM3,AMWM3,AMGM3,AMBM3, mq2_H24      &
      & ,ml2_H24,mHd2_H24,mHu2_H24,md2_H24,mu2_H24,me2_H24,mHw32,mHg32,mHb32   &
      & ,mHx32,mHxb32,MassB_H24,MassWB_H24,MassG_H24,MnuL5)

   MBM3Running(3,:,:)=MBM3
   MGM3Running(3,:,:)=MGM3
   MWM3Running(3,:,:)=MWM3
   MXM3Running(3,:,:)=MXM3

   NGHb3 = 2._dp
   NGHg3 = 2._dp 
   NGHw3 = 2._dp 
   NGHx3 = 2._dp 
   NGHxb3 = 2._dp 
   ThresholdCrossed = 2

   Call FermionMass(MWM3,sqrt2,EigMWM3,RotWL, RotWR,kont)  
   Call FermionMass(MXM3,sqrt2,EigMXM3,RotXL, RotXR,kont)
   Call FermionMass(MGM3,sqrt2,EigMGM3,RotGL, RotGR,kont)
   Call FermionMass(MBM3,sqrt2,EigMBM3,RotBL, RotBR,kont)

   Yb3_H24 = Matmul(RotBL,Yb3_H24)
   Yw3_H24 = Matmul(RotWL,Yw3_H24)
   Yx3_H24 = Matmul(RotXL,Yx3_H24)
   AYb3_H24 = Matmul(RotBL,AYb3_H24)
   AYw3_H24 = Matmul(RotwL,AYb3_H24)
   AYx3_H24 = Matmul(RotxL,AYb3_H24)
   MWM3 = Matmul(Conjg(Transpose(RotWR)),Matmul(MWM3,RotWL))
   MGM3 = Matmul(Conjg(Transpose(RotGR)),Matmul(MGM3,RotGL))
   MXM3 = Matmul(Conjg(Transpose(RotXR)),Matmul(MXM3,RotXL))
   MBM3 = Matmul(Conjg(Transpose(RotBR)),Matmul(MBM3,RotBL))
   mHx32 = Matmul(Conjg(Transpose(RotXL)),Matmul(mHx32,RotXL)) 
   mHxb32 = Matmul(Conjg(Transpose(RotXR)),Matmul(mHxb32,RotXR))  
   mHg32 = Matmul(Conjg(Transpose(RotGL)),Matmul(mHg32,RotGL)) 
   mHb32 = Matmul(Conjg(Transpose(RotBL)),Matmul(mHb32,RotBL))  
   mHw32 = Matmul(Conjg(Transpose(RotWL)),Matmul(mHw32,RotWL))

   MassMWM3(3) = Maxval(EigMWM3)
   MassMXM3(3) = Maxval(EigMXM3)
   MassMGM3(3) = Maxval(EigMGM3)
   MassMBM3(3) = Maxval(EigMBM3)

   If (MassMWM3(3).Lt.Abs(mGUT)) Then
    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = 0.5_dp*Yw3_H24(3,i1)*Yw3_H24(3,i2)/MassMWM3(3) + &
                   & 0.3_dp*Yb3_H24(3,i1)*Yb3_H24(3,i2)/MassMBM3(3) 
     End Do
    End Do
    If (TwoLoopRGE) Then
     g1_h24 = g1_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                &                    * 5._dp/2._dp*Log(MassMXM3(3)/MWM30(3,3)) )
     g2_h24 = g2_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                &                     *( 1.5_dp *Log(MassMXM3(3)/MWM30(3,3)) + &
                &                      2._dp *Log(MassMWM3(3)/MWM30(3,3))) )
     g3_h24 = g3_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                       * (Log(MassMXM3(3)/MWM30(3,3)) &
                &                        + 3._dp*Log(MassMGM3(3)/MWM30(3,3)) ) )

     MassB_h24 = MassB_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                &                    * 5._dp/2._dp*Log(MassMXM3(3)/MWM30(3,3)) )
     MassWB_h24 = MassWB_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                &                     *( 1.5_dp *Log(MassMXM3(3)/MWM30(3,3)) + &
                &                     2._dp *Log(MassMWM3(3)/MWM30(3,3))) )
     MassG_h24 = MassG_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                      * (Log(MassMXM3(3)/MWM30(3,3)) &
                &                        + 3._dp*Log(MassMGM3(3)/MWM30(3,3)) ) )
    End If

   Else 

    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = 0.5_dp * Yw3_H24(3,i1)*Yw3_H24(3,i2)/MWM30(3,3) + &
                   & 0.3_dp * Yb3_H24(3,i1)*Yb3_H24(3,i2)/MWM30(3,3) 
     End Do
    End Do

   End If

   Yb30_H24(3,:,:) = Yb3_H24
   Yw30_H24(3,:,:) = Yw3_H24
   Yx30_H24(3,:,:) = Yx3_H24
   Yb3_H24(3,:) = 0._dp
   Yw3_H24(3,:) = 0._dp
   Yx3_H24(3,:) = 0._dp
   AYb3_H24(3,:) = 0._dp
   AYw3_H24(3,:) = 0._dp
   AYx3_H24(3,:) = 0._dp
   MWM3(:,3) = 0._dp
   MWM3(3,:) = 0._dp
   MGM3(:,3) = 0._dp
   MGM3(3,:) = 0._dp
   MXM3(:,3) = 0._dp
   MXM3(3,:) = 0._dp
   MBM3(:,3) = 0._dp
   MBM3(3,:) = 0._dp

   mHx32(3,:) = 0._dp 
   mHx32(:,3) = 0._dp 
   mHxb32(3,:) = 0._dp 
   mHxb32(:,3) = 0._dp 
   mHg32(3,:) = 0._dp 
   mHg32(:,3) = 0._dp 
   mHb32(3,:) = 0._dp 
   mHb32(:,3) = 0._dp 
   mHw32(3,:) = 0._dp 
   mHw32(:,3) = 0._dp 

   Call ParametersToG555(g1_H24,g2_H24,g3_H24,Yu_H24,Yd_H24,Ye_H24,Yb3_H24   &
      & ,Yw3_H24,Yx3_H24,mu_H24,MXM3,MWM3,MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24  &
      & ,AYb3_H24,AYw3_H24,AYx3_H24,Amue,AMXM3,AMWM3,AMGM3,AMBM3, mq2_H24    &
      & ,ml2_H24,mHd2_H24,mHu2_H24,md2_H24,mu2_H24,me2_H24,mHw32,mHg32,mHb32 &
      & ,mHx32,mHxb32,MassB_H24,MassWB_H24,MassG_H24,MnuL5,g2f)

   !------------
   ! 3 -> 2
   !------------
   If (Abs(MWM30(2,2)).Lt.Abs(m_hi)) Then
    m_lo = Abs(MWM30(2,2))
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2f, 573, tz, 0._dp, delta, dt, 0._dp, rge555, kont)
    m_hi = m_lo
   Endif

   Call GToParameters555(g2f,g1_H24,g2_H24,g3_H24,Yu_H24,Yd_H24,Ye_H24,Yb3_H24 &
      & ,Yw3_H24,Yx3_H24,mu_H24,MXM3,MWM3,MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24    &
      & ,AYb3_H24,AYw3_H24,AYx3_H24,Amue,AMXM3,AMWM3,AMGM3,AMBM3, mq2_H24      &
      & ,ml2_H24,mHd2_H24,mHu2_H24,md2_H24,mu2_H24,me2_H24,mHw32,mHg32,mHb32   &
      & ,mHx32,mHxb32,MassB_H24,MassWB_H24,MassG_H24,MnuL5)

   MBM3Running(2,:,:)=MBM3
   MGM3Running(2,:,:)=MGM3
   MWM3Running(2,:,:)=MWM3
   MXM3Running(2,:,:)=MXM3

   NGHb3 = 1._dp
   NGHg3 = 1._dp 
   NGHw3 = 1._dp 
   NGHx3 = 1._dp 
   NGHxb3 = 1._dp 
   ThresholdCrossed = 1

   RotWR = 0._dp
   RotWL = 0._dp
   RotXR = 0._dp
   RotXL = 0._dp
   RotGR = 0._dp
   RotGL = 0._dp
   RotBR = 0._dp
   RotBL = 0._dp
   EigMWM3 = 0._dp
   EigMGM3 = 0._dp
   EigMXM3 = 0._dp
   EigMBM3 = 0._dp

   Call EigenSystem(Matmul(Transpose(Conjg(MWM3(1:2,1:2))),MWM3(1:2,1:2)) &
    & , EigMWM3(1:2),RotWR(1:2,1:2), ierr, test2)
   Call EigenSystem(Matmul(MWM3(1:2,1:2),Transpose(Conjg(MWM3(1:2,1:2)))) &
    & , EigMWM3(1:2),RotWL(1:2,1:2), ierr, test2)

   Call EigenSystem(Matmul(Transpose(Conjg(MBM3(1:2,1:2))),MBM3(1:2,1:2)) &
    & , EigMBM3(1:2),RotBR(1:2,1:2), ierr, test2)
   Call EigenSystem(Matmul(MBM3(1:2,1:2),Transpose(Conjg(MBM3(1:2,1:2)))) &
    & , EigMBM3(1:2),RotBL(1:2,1:2), ierr, test2)

   Call EigenSystem(Matmul(Transpose(Conjg(MGM3(1:2,1:2))),MGM3(1:2,1:2)) &
    & , EigMGM3(1:2),RotGR(1:2,1:2), ierr, test2)
   Call EigenSystem(Matmul(MGM3(1:2,1:2),Transpose(Conjg(MGM3(1:2,1:2)))) &
    & , EigMGM3(1:2),RotGL(1:2,1:2), ierr, test2)

   Call EigenSystem(Matmul(Transpose(Conjg(MXM3(1:2,1:2))),MXM3(1:2,1:2)) &
    & , EigMXM3(1:2),RotXR(1:2,1:2), ierr, test2)
   Call EigenSystem(Matmul(MXM3(1:2,1:2),Transpose(Conjg(MXM3(1:2,1:2)))) &
    & , EigMXM3(1:2),RotXL(1:2,1:2), ierr, test2)


   Yb3_H24 = Matmul(RotBL,Yb3_H24)
   Yw3_H24 = Matmul(RotWL,Yw3_H24)
   Yx3_H24 = Matmul(RotXL,Yx3_H24)
   AYb3_H24 = Matmul(RotBL,AYb3_H24)
   AYw3_H24 = Matmul(RotwL,AYb3_H24)
   AYx3_H24 = Matmul(RotxL,AYb3_H24)
   MWM3 = Matmul(Conjg(Transpose(RotWR)),Matmul(MWM3,RotWL))
   MGM3 = Matmul(Conjg(Transpose(RotGR)),Matmul(MGM3,RotGL))
   MXM3 = Matmul(Conjg(Transpose(RotXR)),Matmul(MXM3,RotXL))
   MBM3 = Matmul(Conjg(Transpose(RotBR)),Matmul(MBM3,RotBL))
   mHx32 = Matmul(Conjg(Transpose(RotXL)),Matmul(mHx32,RotXL)) 
   mHxb32 = Matmul(Conjg(Transpose(RotXR)),Matmul(mHxb32,RotXR))  
   mHg32 = Matmul(Conjg(Transpose(RotGL)),Matmul(mHg32,RotGL)) 
   mHb32 = Matmul(Conjg(Transpose(RotBL)),Matmul(mHb32,RotBL))  
   mHw32 = Matmul(Conjg(Transpose(RotWL)),Matmul(mHw32,RotWL))
  
   MassMWM3(2) = Sqrt(EigMWM3(2))
   MassMXM3(2) = Sqrt(EigMXM3(2))
   MassMGM3(2) = Sqrt(EigMGM3(2))
   MassMBM3(2) = Sqrt(EigMBM3(2))

   If (MassMWM3(2).Lt.Abs(mGUT)) Then

    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = MnuL5(i1,i2) &
                & + 0.5_dp * Yw3_H24(2,i1)*Yw3_H24(2,i2)/MassMWM3(2)  &
                & + 0.3_dp * Yb3_H24(2,i1)*Yb3_H24(2,i2)/MassMBM3(2) 
     End Do
    End Do

    If (TwoLoopRGE) Then
     g1_h24 = g1_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                &                   * 5._dp/2._dp*Log( MassMXM3(2)/MWM30(2,2)) )
     g2_h24 = g2_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                &                    *( 1.5_dp *Log(MassMXM3(2)/MWM30(2,2)) + & 
                &                      2._dp *Log(MassMWM3(2)/MWM30(2,2))) )
     g3_h24 = g3_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                     * (Log( MassMXM3(2)/MWM30(2,2)) &
                &                       + 3._dp*Log( MassMGM3(2)/MWM30(2,2)) ) )

     MassB_h24 = MassB_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                 &                  * 5._dp/2._dp*Log( MassMXM3(2)/MWM30(2,2)) )
     MassWB_h24 = MassWB_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                 &                    *( 1.5_dp *Log(MassMXM3(2)/MWM30(2,2)) + &
                 &                       2._dp *Log(MassMWM3(2)/MWM30(2,2))) )
     MassG_h24 = MassG_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                     * (Log( MassMXM3(2)/MWM30(2,2)) &
                &                       + 3._dp*Log( MassMGM3(2)/MWM30(2,2)) ) )
    End If

   Else 

    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = MnuL5(i1,i2)   &
                 & + 0.5_dp * Yw3_H24(2,i1)*Yw3_H24(2,i2)/mGut  &
                 & + 0.3_dp * Yb3_H24(2,i1)*Yb3_H24(2,i2)/mGut 
     End Do
    End Do

   End If

   Yb30_H24(2,:,:) = Yb3_H24
   Yw30_H24(2,:,:) = Yw3_H24
   Yx30_H24(2,:,:) = Yx3_H24
   Yb3_H24(2,:) = 0._dp
   Yw3_H24(2,:) = 0._dp
   Yx3_H24(2,:) = 0._dp
   AYb3_H24(2,:) = 0._dp
   AYw3_H24(2,:) = 0._dp
   AYx3_H24(2,:) = 0._dp
   MWM3(:,2) = 0._dp
   MWM3(2,:) = 0._dp
   MGM3(:,2) = 0._dp
   MGM3(2,:) = 0._dp
   MXM3(:,2) = 0._dp
   MXM3(2,:) = 0._dp
   MBM3(:,2) = 0._dp
   MBM3(2,:) = 0._dp

   mHx32(2,:) = 0._dp 
   mHx32(:,2) = 0._dp 
   mHxb32(2,:) = 0._dp 
   mHxb32(:,2) = 0._dp 
   mHg32(2,:) = 0._dp 
   mHg32(:,2) = 0._dp 
   mHb32(2,:) = 0._dp 
   mHb32(:,2) = 0._dp 
   mHw32(2,:) = 0._dp 
   mHw32(:,2) = 0._dp 

   Call ParametersToG555(g1_H24,g2_H24,g3_H24,Yu_H24,Yd_H24,Ye_H24,Yb3_H24 &
      & ,Yw3_H24,Yx3_H24,mu_H24,MXM3,MWM3,MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24    &
      & ,AYb3_H24,AYw3_H24,AYx3_H24,Amue,AMXM3,AMWM3,AMGM3,AMBM3, mq2_H24      &
      & ,ml2_H24,mHd2_H24,mHu2_H24,md2_H24,mu2_H24,me2_H24,mHw32,mHg32,mHb32   &
      & ,mHx32,mHxb32,MassB_H24,MassWB_H24,MassG_H24,MnuL5,g2f)

   !------------
   ! 2 -> 1
   !------------

   If (Abs(MWM30(1,1)).Lt.Abs(m_hi)) Then
    m_lo = MWM30(1,1)
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2f, 573, tz, 0._dp, delta, dt, 0._dp, rge555, kont)
    m_hi = m_lo
   Endif

   Call GToParameters555(g2f,g1_H24,g2_H24,g3_H24,Yu_H24,Yd_H24,Ye_H24,Yb3_H24 &
      & ,Yw3_H24,Yx3_H24,mu_H24,MXM3,MWM3,MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24    &
      & ,AYb3_H24,AYw3_H24,AYx3_H24,Amue,AMXM3,AMWM3,AMGM3,AMBM3, mq2_H24      &
      & ,ml2_H24,mHd2_H24,mHu2_H24,md2_H24,mu2_H24,me2_H24,mHw32,mHg32,mHb32   &
      & ,mHx32,mHxb32,MassB_H24,MassWB_H24,MassG_H24,MnuL5)

   MBM3Running(1,:,:)=MBM3
   MGM3Running(1,:,:)=MGM3
   MWM3Running(1,:,:)=MWM3
   MXM3Running(1,:,:)=MXM3

   NGHb3 = 0._dp
   NGHg3 = 0._dp 
   NGHw3 = 0._dp 
   NGHx3 = 0._dp 
   NGHxb3 = 0._dp 
   ThresholdCrossed = 0

   MassMWM3(1) = Abs(MWM3(1,1))
   MassMXM3(1) = Abs(MXM3(1,1))
   MassMGM3(1) = Abs(MGM3(1,1))
   MassMBM3(1) = Abs(MBM3(1,1))

   Yb30_H24(1,:,:) = Yb3_H24
   Yw30_H24(1,:,:) = Yw3_H24
   Yx30_H24(1,:,:) = Yx3_H24
 
   gauge_h24(1) = g1_h24
   gauge_h24(2) = g2_h24
   gauge_h24(3) = g3_h24
   Mi_h24(1) = MassB_h24
   Mi_h24(2) = MassWB_h24
   Mi_h24(3) = MassG_h24
   mh2_h24(1) = mHd2_h24
   mh2_h24(2) = mHu2_h24

   If (MassMWM3(1).Lt.Abs(mGUT)) Then

    If (TwoLoopRGE) Then
     g1_h24 = g1_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                &                  * 5._dp/12._dp*Log( MassMXM3(1)/MWM30(1,1)) )
     g2_h24 = g2_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                &                     *( 0.5_dp *Log(MassMXM3(1)/MWM30(1,1)) + &
                &                       2._dp *Log(MassMWM3(1)/MWM30(1,1))) )
     g3_h24 = g3_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                     * (0.5_dp*Log( MassMXM3(1)/MWM30(1,1)) &
                &                       + 3._dp*Log( MassMGM3(1)/MWM30(1,1)) ) )

     MassB_h24 = MassB_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                &                  * 5._dp/12._dp*Log( MassMXM3(1)/MWM30(1,1)) )
     MassWB_h24 = MassWB_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                 &                    *( 0.5_dp *Log(MassMXM3(1)/MWM30(1,1)) + &
                 &                       2._dp *Log(MassMWM3(1)/MWM30(1,1))) )
     MassG_h24 = MassG_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                     * (0.5_dp*Log( MassMXM3(1)/MWM30(1,1)) &
                &                       + 3._dp*Log( MassMGM3(1)/MWM30(1,1)) ) )
    End If
 
    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = MnuL5(i1,i2) &
                & + 0.5_dp * Yw3_H24(1,i1)*Yw3_H24(1,i2)/ MassMWM3(1)  &
                & + 0.3_dp * Yb3_H24(1,i1)*Yb3_H24(1,i2)/ MassMBM3(1) 
     End Do
    End Do

   Else 

    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = MnuL5(i1,i2)   &
                & + 0.5_dp * Yw3_H24(1,i1)*Yw3_H24(1,i2)/Abs(MWM30(1,1))  &
                & + 0.3_dp * Yb3_H24(1,i1)*Yb3_H24(1,i2)/Abs(MWM30(1,1)) 
     End Do
    End Do

   End If


  Call ParametersToG4(gauge_H24, ye_H24, Zero33C, yd_H24, yu_H24   &
          & , ZeroC, ZeroC, Mi_H24, AYe_h24, Zero33C, AYd_h24 &
          & , AYu_h24, ZeroC, ZeroC, Me2_h24, ML2_h24     &
          & , MD2_h24, MQ2_h24, MU2_h24, MH2_h24, ZeroR2, mu_h24      &
          & , Amue, MnuL5, g2c)

  mudim = GetRenormalizationScale()
  mudim = Max(mudim, mZ2)

  tz = 0.5_dp * Log(mudim/m_hi**2)
  dt = tz / 100._dp
  Call odeint(g2c, 277, 0._dp, tz, delta, dt, 0._dp, rge277, kont)

  Call GToParameters4(g2c, gauge, y_l, y_T, y_d, y_u, lam12(1), lam12(2), Mi &
          & , A_l, A_T, A_d, A_u, Alam12(1), Alam12(2), M2_E, M2_L, M2_D      &
          & , M2_Q, M2_U, M2_H, M2_T, mu, B, MnuL5)

  Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u, M2_E, M2_L  &
      & , M2_D, M2_Q, M2_U, M2_H, mu, B, g2)


 Else If (HighScaleModel.Eq.'SEESAW_II_SARAH') Then

   If ( (oo4pi*Maxval(g2g(1:115)**2)).Gt.0.5_dp) Then
    Write(ErrCan,*) "Non perturbative regime at M_GUT"
    If (ErrorLevel.Ge.2) Call TerminateProgram
    Write(errcan,*) " "
    kont = -410
    Call AddError(410)
    Iname = Iname - 1
    Return
   End If

   !-------------------------
   ! run only if m_H < m_GUT
   !-------------------------

   If (Abs(MTM_GUT).Lt.mGUT) Then
    !---------------------------------------------------
    ! heavy states are to be included
    !---------------------------------------------------
    

    mudim = Abs(MTM_GUT)**2
    tz = 0.5_dp * Log(mudim/mGUT**2)
    dt = tz / 50._dp
    Call odeint(g2g, 365 , 0._dp, tz, delta, dt, 0._dp, rge365, kont)
    m_lo = Abs(MTM_GUT)

    If ( (oo4pi*Maxval(g2g(1:115)**2)).Gt.0.5_dp) Then
     Write(ErrCan,*) "Non perturbative regime at M_H3"
     If (ErrorLevel.Ge.2) Call TerminateProgram
     Write(errcan,*) " "
     kont = -411
     Call AddError(411)
     Iname = Iname - 1
     Return
    End If

   Else
     m_lo = mGUT
   End If


   Call GToParameters365(g2g,g1_h15,g2_h15,g3_h15,Yu_h15,Yd_h15,Ye_h15,Yt0_h15 &
      & ,Ys0_h15,Yz0_h15, Lambda10,Lambda20,mu_h15,MTM,MZM,MSM,AYu_h15,AYd_h15 &
      & ,AYe_h15,AYt_h15,AYs_h15,AYz_h15,ALambda1,ALambda2,Amue,AMTM,AMZM,AMSM &
      & ,mq2_h15,ml2_h15,mHd2_h15, mHu2_h15,md2_h15,mu2_h15,me2_h15,mt2_H15    &
      & ,mtb2_H15,ms2_H15,msb2_H15,mz2_H15,mzb2_H15,MassB_h15,MassWB_h15       &
      & ,MassG_h15,MnuL5)

   gauge_h15(1) = g1_h15
   gauge_h15(2) = g2_h15
   gauge_h15(3) = g3_h15
   Mi_h15(1) = MassB_h15
   Mi_h15(2) = MassWB_h15
   Mi_h15(3) = MassG_h15
   mh2_h15(1) = mHd2_h15
   mh2_h15(2) = mHu2_h15

   MnuL5 = - Lambda20 * YT0_h15 / MTM

   If ((Abs(MTM_GUT).Lt.mGUT).And.TwoLoopRGE) Then
   !-----------------------------------------------------
   ! adding shifts to gauge couplings and gaugino masses
   !-----------------------------------------------------
    gauge_h15(1) = gauge_h15(1) / (1._dp - 0.2_dp * oo16pi2 * gauge_h15(1)**2 &
                &                       * ( 8._dp * Log(Abs(MSM/MTM_GUT))     &
                &                         + 9._dp * Log(Abs(MTM/MTM_GUT))     &
                &                         + 0.5_dp * Log(Abs(MZM/MTM_GUT)) ) )
    gauge_h15(2) = gauge_h15(2) / (1._dp - oo16pi2 * gauge_h15(2)**2          &
              & *(2._dp*Log(Abs(MTM/MTM_GUT))+ 1.5_dp *Log(Abs(MZM/MTM_GUT)) ))
    gauge_h15(3) = gauge_h15(3) / (1._dp - oo16pi2 * gauge_h15(3)**2          &
                &   * (2.5_dp*Log(Abs(MSM/MTM_GUT)) + Log(Abs(MZM/MTM_GUT)) ) )
    Mi_h15(1) = Mi_h15(1) / (1._dp - 0.2_dp * oo16pi2 * gauge_h15(1)**2       &
                &                       * ( 8._dp * Log(Abs(MSM/MTM_GUT))     &
                &                         + 9._dp * Log(Abs(MTM/MTM_GUT))     &
                &                         + 0.5_dp * Log(Abs(MZM/MTM_GUT)) ) )
    Mi_h15(2) = Mi_h15(2) / (1._dp - oo16pi2 * gauge_h15(2)**2                &
             & *(2._dp *Log(Abs(MTM/MTM_GUT))+ 1.5_dp *Log(Abs(MZM/MTM_GUT)) ))
    Mi_h15(3) = Mi_h15(3) / (1._dp - oo16pi2 * gauge_h15(3)**2           &
                &   * (2.5_dp*Log(Abs(MSM/MTM_GUT)) + Log(Abs(MZM/MTM_GUT)) ) )
   End If

   !---------------------------------------------------
   ! heavy states are not to be included
   !---------------------------------------------------

   Call ParametersToG4(gauge_h15, ye_h15, Zero33C, yd_h15, yu_h15   &
          & , ZeroC, ZeroC, Mi_h15, AYe_h15, Zero33C, AYd_h15 &
          & , AYu_h15, ZeroC, ZeroC, ME2_h15, ML2_h15     &
          & , MD2_h15, MQ2_h15, MU2_h15, MH2_h15, ZeroR2, mu_h15      &
          & , Amue, MnuL5, g2c)

   mudim = GetRenormalizationScale()
   mudim = Max(mudim, mZ2)

   tz = 0.5_dp * Log(mudim/m_lo**2)
   dt = tz / 100._dp
   Call odeint(g2c, 277, 0._dp, tz, delta, dt, 0._dp, rge277, kont)



   Call GToParameters4(g2c, gauge, y_l, y_T, y_d, y_u, lam12(1), lam12(2), Mi &
          & , A_l, A_T, A_d, A_u, Alam12(1), Alam12(2), M2_E, M2_L, M2_D      &
          & , M2_Q, M2_U, M2_H, M2_T, mu, B, MnuL5)
   Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u, M2_E, M2_L  &
      & , M2_D, M2_Q, M2_U, M2_H, mu, B, g2)

# endif SEESAWIII
! Florian Staub Seesaw II+III

  Else

   mudim = GetRenormalizationScale()
   mudim = Max(mudim, mZ2)

   tz = 0.5_dp * Log(mudim/mGUT_save**2)
   dt = tz / 100._dp

   Call odeint(g2, 213, 0._dp, tz, delta, dt, 0._dp, rge213, kont)

  End If

  Iname = Iname - 1


  900 Format(a20,e15.6)
  910 Format(a15,3e15.6)

 End Subroutine RunRGE

 Subroutine RunRGE_2_boundaries(kont, i_in, delta, Qin, g1, g2, mGUT)
 !-----------------------------------------------------------------------
 ! Uses Runge-Kutta method to integrate RGE's from M_Z to M_GUT
 ! and back, putting in correct thresholds. For the first iteration
 ! only the first 6 couplings are included and a common threshold is used.
 ! Written by Werner Porod, 10.07.99
 ! 07.03.2001: including right handed neutrinos
 ! 24.09.01: portation to f90
 ! 27.03.02: including a check if perturbation theory is valid
 ! 16.09.02: the electroweak scale is now set entirely outside, either 
 !           in the main program or in the routine sugra(...)
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in
  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: delta, Qin
  Real(dp), Intent(inout) :: g1(213)
  Real(dp), Intent(out) :: g2(213), mGUT
  
  Integer:: i1, i2, SumI
  Real(dp) :: g1a(285), g2a(285), g1b(75), g2b(267) &
      & , g2c(277), g1d(356), g2d(356), g1f(573), g2f(573)     &
      & , g1g(365), g2g(365)
  Real(dp) :: tz, dt, t_out
  Real(dp) :: mudim, gGUT, g1_h(213), m_hi, m_lo, M15
  Logical :: FoundUnification

  Real(dp), Parameter :: Umns(3,3) = Reshape(   Source = (/  &
      &    Sqrt2/Sqrt3, -ooSqrt2*ooSqrt3, -ooSqrt2*ooSqrt3   &
      &  , ooSqrt3,      ooSqrt3,          ooSqrt3           &
      &  , 0._dp ,       ooSqrt2,         -ooSqrt2 /), shape = (/3, 3/) )
  Real(dp), Parameter :: ZeroR2(2) = 0._dp
  Complex(dp), Dimension(3,3) :: mat3, UnuR, Ynu, Anu, Mr2

  Complex(dp), Dimension(3,3) :: RotXl, RotXr, RotGl, RotGr, RotBl, RotBr &
      & , RotWl, RotWr, ZER, ZEL
  Real(dp) :: EigMWM3(3),EigMBM3(3),EigMXM3(3), EigMGM3(3), test2(2)
  Integer :: ierr

  Complex(dp), Dimension(3,3,3) :: mHw32_s, mHg32_s, mHb32_s, mHx32_s, mHxb32_s &
      & , AYb3_s, AYw3_s, AYx3_s

  Iname = Iname + 1
  NameOfUnit(Iname) = 'RunRGE_2_boundaries'

  !-------------------------------------
  ! running to the high scale
  !-------------------------------------
  g1_h = g1
  If (HighScaleModel(1:9).Eq.'SUGRA_NuR') &
  Then
   m_lo = MnuR(1)
   FoundUnification = .False.
  Else If (HighScaleModel.Eq.'SEESAW_II') Then
   m_lo = M_H3(1)
   FoundUnification = .False.

! Florian Staub Seesaw II+III
# ifdef SEESAWIII
  Else If (HighScaleModel.Eq.'SEESAW_III_3G') Then
   m_lo = Abs(MWM30(1,1))
   FoundUnification = .False.
  Else If (HighScaleModel.Eq.'SEESAW_II_SARAH') Then
   m_lo = Abs(MTM_GUT)
   FoundUnification = .False.
# endif SEESAWIII
! Florian Staub Seesaw II+III

  Else If (UseFixedGUTScale) Then ! GUT scale is fixed
   m_lo = GUT_scale
   mGUT = GUT_scale
   FoundUnification = .True.
  Else ! Sugra, strings with minimal particle content
   m_lo = 5.e14_dp
   FoundUnification = .False.
  End If

  tz = Log(Qin/m_lo)
  dt = - tz / 50._dp

  Call odeint(g1, 213, tz, 0._dp, delta, dt, 0._dp, rge213, kont)
  If (kont.Ne.0) Then
   Iname = Iname - 1
   Return
  End If
  !--------------------------
  ! check for perturbativity
  !--------------------------
  If ( (oo4pi*Maxval(g1(1:57)**2)).Gt.0.5_dp) Then
   Write(ErrCan,*) "Non perturbative regime at high scale"
   If (ErrorLevel.Ge.2) Call TerminateProgram
   Write(errcan,*) " "
   kont = -416
   Call AddError(416)
   Iname = Iname - 1
   Return
  End If
    
  !---------------------------
  ! looking for the GUT scale
  !---------------------------
  If (HighScaleModel.Eq.'SUGRA') Then 

   If (.Not.UseFixedGUTScale) Then
    tz = Log(m_lo/1.e18_dp)
    dt = - tz / 50._dp

    Call odeintB(g1, 213, tz, 0._dp, delta, dt, 0._dp, rge213, t_out, kont)
    If (kont.Eq.0) Then
     FoundUnification = .True.
     mGUT = 1.e18_dp * Exp(t_out)
     gGUT = Sqrt( 0.5_dp * (g1(1)**2+g1(2)**2) )
     g1(1) = gGUT
     g1(2) = gGUT
     If (StrictUnification) g1(3) = gGUT
    Else
     Write(ErrCan,*) "kont",kont,delta,tz,dt
     Write(ErrCan,*) "m_t",mf_u(3)
     Write (ErrCan,*) "t_out",t_out,1.e18_dp * Exp(t_out)
     Do i1=1,213
      If ((g1(i1).Ne.0._dp).Or.(g1_h(i1).Ne.0._dp)) &
                 & Write(ErrCan,*) i1,g1_h(i1),g1(i1)
     End Do
     Write(ErrCan,*) " " 
     Iname = Iname - 1
     Return
    End If
   End If ! .not.UseFixedGUTScale

  Else If (HighScaleModel.Eq.'SUGRA_NuR')  Then

   Call GToParameters(g1, gauge_mR,Y_l_mR(1,:,:),Y_d_mR(1,:,:),Y_u_mR(1,:,:)   &
      & , Mi_mR, A_l_mR(1,:,:), A_d_mR(1,:,:), A_u_mR(1,:,:), M2_E_mR(1,:,:)   &
      & , M2_L_mR(1,:,:), M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:)       &
      & , M2_H_mR, mu_mR, B_mR)

   Call ParametersToG3(gauge_mR, Y_l_mR(1,:,:), Y_nu_mR(1,:,:), Y_d_mR(1,:,:)  &
      & , Y_u_mR(1,:,:), Mi_mR, A_l_mR(1,:,:), A_nu_mR(1,:,:), A_d_mR(1,:,:)   &
      & , A_u_mR(1,:,:), M2_E_mR(1,:,:), M2_L_mR(1,:,:), M2_R_mR(1,:,:)        &
      & , M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:), M2_H_mR, mu_mR, B_mR &
      & , MnuL5, g1a)
   !---------------------------
   ! running from m_R1 -> m_R2
   !---------------------------
   m_lo = MnuR(1)
   m_hi = MnuR(2)
   If (MnuR(1).Ne.MnuR(2)) Then
    tz = Log(m_lo / m_hi)
    dt = - tz / 50._dp  
    Call odeint(g1a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)
   End If

   Call GToParameters3(g2a, gauge_mR, y_l_mR(1,:,:), y_nu_mR(1,:,:)          &
      & , y_d_mR(1,:,:), y_u_mR(1,:,:), Mi_mR, A_l_mR(1,:,:), A_nu_mR(1,:,:) &
      & , A_d_mR(1,:,:), A_u_mR(1,:,:), M2_E_mR(1,:,:), M2_L_mR(1,:,:)       &
      & , M2_R_mR(1,:,:), M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:)     &
      & , M2_H_mR, mu_mR, B_mR, MnuL5)

   Call ParametersToG3(gauge_mR, Y_l_mR(2,:,:), Y_nu_mR(2,:,:), Y_d_mR(2,:,:)  &
      & , Y_u_mR(2,:,:), Mi_mR, A_l_mR(2,:,:), A_nu_mR(2,:,:), A_d_mR(2,:,:)   &
      & , A_u_mR(2,:,:), M2_E_mR(2,:,:), M2_L_mR(2,:,:), M2_R_mR(2,:,:)        &
      & , M2_D_mR(2,:,:), M2_Q_mR(2,:,:), M2_U_mR(2,:,:), M2_H_mR, mu_mR, B_mR &
      & , MnuL5, g1a)

   m_lo = m_hi
   !---------------------------
   ! running from m_R2 -> m_R3
   !---------------------------
   m_hi = MnuR(3)
   If (m_lo.Ne.MnuR(3)) Then
    tz = Log(m_lo / m_hi)
    dt = - tz / 50._dp  

    Call odeint(g1a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)
   End If

   Call GToParameters3(g2a, gauge_mR, y_l_mR(2,:,:), y_nu_mR(2,:,:)          &
      & , y_d_mR(2,:,:), y_u_mR(2,:,:), Mi_mR, A_l_mR(2,:,:), A_nu_mR(2,:,:) &
      & , A_d_mR(2,:,:), A_u_mR(2,:,:), M2_E_mR(2,:,:), M2_L_mR(2,:,:)       &
      & , M2_R_mR(2,:,:), M2_D_mR(2,:,:), M2_Q_mR(2,:,:), M2_U_mR(2,:,:)     &
      & , M2_H_mR, mu_mR, B_mR, MnuL5)

   If (Ynu_at_MR3) Then
    Call FermionMass(Y_l_mR(2,:,:),1._dp,EigMGM3,ZER,ZEL,kont)
    Y_l_mR(3,:,:) = Matmul( Matmul(Conjg(ZER), Y_l_mR(2,:,:)) &
                  &       , Transpose(Conjg(ZEL)))
    Y_nu_mR(3,:,:) = Y_nu_0
    A_l_mR(3,:,:) = Matmul( Matmul(Conjg(ZER), A_l_mR(2,:,:)) &
                  &       , Transpose(Conjg(ZEL)))
    M2_E_mR(3,:,:) = Matmul( Matmul(ZER, M2_E_mR(2,:,:)) &
                   &       , Transpose(Conjg(ZER)))
    M2_L_mR(3,:,:) = Matmul( Matmul(ZEL, M2_L_mR(2,:,:)) &
                   &       , Transpose(Conjg(ZEL)))
   Else
    Y_l_mR(3,:,:) = Y_l_mR(2,:,:)
    Y_nu_mR(3,:,:) = y_nu_mR(2,:,:)
    A_l_mR(3,:,:) = A_l_mR(2,:,:)
    M2_E_mR(3,:,:) = M2_E_mR(2,:,:)
    M2_L_mR(3,:,:) = M2_L_mR(2,:,:)

   End If

   Call ParametersToG3(gauge_mR, Y_l_mR(3,:,:), Y_nu_mR(3,:,:), Y_d_mR(3,:,:)  &
      & , Y_u_mR(3,:,:), Mi_mR, A_l_mR(3,:,:), A_nu_mR(3,:,:), A_d_mR(3,:,:)   &
      & , A_u_mR(3,:,:), M2_E_mR(3,:,:), M2_L_mR(3,:,:), M2_R_mR(3,:,:)        &
      & , M2_D_mR(3,:,:), M2_Q_mR(3,:,:), M2_U_mR(3,:,:), M2_H_mR, mu_mR, B_mR &
      & , MnuL5, g1a)

   m_lo = m_hi
   !---------------------------
   ! running from m_R3 -> m_GUT
   !---------------------------
   If (UseFixedGUTScale) Then
    tz = Log(MnuR(3)/GUT_scale)
    mGUT = GUT_scale
    dt = - tz / 50._dp
    Call odeint(g1a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)
    If (kont.Ne.0) Then
     Iname = Iname -1
     Return
    End If

   Else
    If (g1a(1).Lt.g1a(2)) Then ! I am still below GUT scale
     tz = Log(MnuR(3)/1.e18_dp)
     dt = - tz / 50._dp
     Call odeintB(g1a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, t_out, kont) 
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e18_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1a(1)**2+g1a(2)**2) )
      g1a(1) = gGUT
      g1a(2) = gGUT
      If (StrictUnification) g1a(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If

    Else If (g1a(1).Eq.g1a(2)) Then ! I am at the GUT scale, very unlikely
                                    ! but possible
     FoundUnification = .True.
     mGUT = 1.e15_dp * Exp(t_out)
     gGUT = g1a(1)
     If (StrictUnification) g1a(3) = gGUT

    Else ! I have already crossed the GUT scale
     tz = Log(MnuR(3)/1.e15_dp)
     dt = - tz / 50._dp
     Call odeintC(g1a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, t_out, kont)
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1a(1)**2+g1a(2)**2) )
      g1a(1) = gGUT
      g1a(2) = gGUT
      If (StrictUnification) g1a(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If
    End If

   End If

  Else If (HighScaleModel.Eq.'SUGRA_NuR1') Then

   Call GToParameters(g1, gauge_mR,Y_l_mR(1,:,:),Y_d_mR(1,:,:),Y_u_mR(1,:,:)   &
      & , Mi_mR, A_l_mR(1,:,:), A_d_mR(1,:,:), A_u_mR(1,:,:), M2_E_mR(1,:,:)   &
      & , M2_L_mR(1,:,:), M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:)       &
      & , M2_H_mR, mu_mR, B_mR)

   If (Ynu_at_MR3) Then
    Call FermionMass(Y_l_mR(1,:,:),1._dp,EigMGM3,ZER,ZEL,kont)
    Y_l_mR(1,:,:) = Matmul( Matmul(Conjg(ZER), Y_l_mR(1,:,:)) &
                  &       , Transpose(Conjg(ZEL)))
    Y_nu_mR(1,:,:) = Y_nu_0
    A_l_mR(1,:,:) = Matmul( Matmul(Conjg(ZER), A_l_mR(1,:,:)) &
                  &       , Transpose(Conjg(ZEL)))
    M2_E_mR(1,:,:) = Matmul( Matmul(ZER, M2_E_mR(1,:,:)) &
                   &       , Transpose(Conjg(ZER)))
    M2_L_mR(1,:,:) = Matmul( Matmul(ZEL, M2_L_mR(1,:,:)) &
                   &       , Transpose(Conjg(ZEL)))
   End If

   Call ParametersToG3(gauge_mR, Y_l_mR(1,:,:), Y_nu_mR(1,:,:), Y_d_mR(1,:,:)  &
      & , Y_u_mR(1,:,:), Mi_mR, A_l_mR(1,:,:), A_nu_mR(1,:,:), A_d_mR(1,:,:)   &
      & , A_u_mR(1,:,:), M2_E_mR(1,:,:), M2_L_mR(1,:,:), M2_R_mR(1,:,:)        &
      & , M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:), M2_H_mR, mu_mR, B_mR &
      & , MnuL5, g1a)

   tz = Log(mGUT/GUT_scale)
   dt = tz / 50._dp  

   If (UseFixedGUTScale) Then
    tz = Log(mGUT/GUT_scale)
    mGUT = GUT_scale
    dt = - tz / 50._dp
    Call odeint(g1a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)
    If (kont.Ne.0) Then
     Iname = Iname -1
     Return
    End If

   Else
    tz = Log(mGUT/1.e18_dp)
    dt = - tz / 50._dp
    Call odeintB(g1a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, t_out, kont)
    If (kont.Eq.0) Then
     FoundUnification = .True.
     mGUT = 1.e18_dp * Exp(t_out)
     gGUT = Sqrt( 0.5_dp * (g1a(1)**2+g1a(2)**2) )
     g1a(1) = gGUT
     g1a(2) = gGUT
     If (StrictUnification) g1a(3) = gGUT
    Else
     Iname = Iname - 1
     Return
    End If
   End If

  Else If ((HighScaleModel.Eq.'SEESAW_II').And.Fifteen_plet) Then

   If (g1(1).Gt.g1(2)) Then ! already above the GUT scale
    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
    Call Set_Decoupling_Heavy_States(.True.)

    If (UseFixedGUTScale) Then
     tz = Log(m_lo/GUT_scale)
     mGUT = GUT_scale
     dt = - tz / 50._dp
     Call odeint(g1, 213, tz, 0._dp, delta, dt, 0._dp, rge213, kont)
     If (kont.Ne.0) Then
      Iname = Iname -1
      Return
     End If

    Else ! .not.UseFixedGUTScale
     tz = Log(m_lo/1.e15_dp)
     dt = - tz / 50._dp
     Call odeintC(g1, 213, tz, 0._dp, delta, dt, 0._dp, rge213, t_out, kont)
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1(1)**2+g1(2)**2) )
      g1(1) = gGUT
      g1(2) = gGUT
      If (StrictUnification) g1(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If

    End If ! UseFixedGUTScale

    Call GToParameters(g1, gauge_mH3, Y_l_mH3, Y_d_mH3, Y_u_mH3, Mi_mH3  &
      & , A_l_mH3, A_d_mH3, A_u_mH3, M2_E_mH3   &
      & , M2_L_mH3, M2_D_mH3, M2_Q_mH3, M2_U_mH3       &
      & , M2_H_mH3, mu_mH3, B_mH3)

    MT15_mH3 = M_H3(1)
    Call ParametersToG5(gauge_mH3, Y_l_mH3, Y_T_0, Y_d_mH3, Y_u_mH3 &
          & , Y_T_0, Y_T_0, lam12_0(1), lam12_0(2), Mi_mH3, A_l_mH3   &
          & , A_T_mH3, A_d_mH3, A_u_mH3, A_Z_mH3, A_S_mH3, Alam12_MH3(1)      &
          & , Alam12_MH3(2), M2_E_mH3, M2_L_mH3, M2_D_mH3, M2_Q_mH3, M2_U_mH3 &
          & , M2_H_mH3, M2_T_mH3, M2_Z_mH3, M2_S_mH3, MT15_mH3, MT15_mH3      &
          & , MT15_mH3, mu_mH3, B_mH3, MnuL5, g1d )

   Else ! still below the GUT scale

    !---------------------------------------------------
    ! heavy states are to be included
    !---------------------------------------------------
    Call Set_Decoupling_Heavy_States(.False.)

    Call GToParameters(g1, gauge_mR,Y_l_mR(1,:,:),Y_d_mR(1,:,:),Y_u_mR(1,:,:)  &
      & , Mi_mR, A_l_mR(1,:,:), A_d_mR(1,:,:), A_u_mR(1,:,:), M2_E_mR(1,:,:)   &
      & , M2_L_mR(1,:,:), M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:)       &
      & , M2_H_mR, mu_mR, B_mR)


    M15 = M_H3(1)

    Delta_b_1 = 7._dp
    If (TwoLoopRGE) Then
     Delta_b_2(1,1) = 181._dp/15._dp
     Delta_b_2(1,2) = 29.4_dp 
     Delta_b_2(1,3) = 656._dp/15._dp
     Delta_b_2(2,1) = 9.8_dp
     Delta_b_2(2,2) = 69._dp
     Delta_b_2(2,3) = 16._dp
     Delta_b_2(3,1) = 82._dp/15._dp
     Delta_b_2(3,2) = 6._dp
     Delta_b_2(3,3) = 358._dp/3._dp

    !-----------------------------------------------------
    ! adding shifts to gauge couplings
    !-----------------------------------------------------
     gauge_mH3(1) = gauge_mH3(1) * (1._dp - oo16pi2 * gauge_mH3(1)**2          &
                 &                       * (8._dp/3._dp*Log(MS15_mH3/MT15_mH3) &
                 &                         + Log(MZ15_mH3/MT15_mH3) /6._dp ) )
     gauge_mH3(2) = gauge_mH3(2) * (1._dp - oo16pi2 * gauge_mH3(2)**2          &
                 &                       * 1.5_dp *Log(MZ15_mH3/MT15_mH3) )
     gauge_mH3(3) = gauge_mH3(3) * (1._dp - oo16pi2 * gauge_mH3(3)**2          &
                 &                       * (2.5_dp*Log(MS15_mH3/MT15_mH3) &
                 &                         + Log(MZ15_mH3/MT15_mH3) ) )
    End If

    If (Ynu_at_MR3) Then
     Call FermionMass(Y_l_mH3,1._dp,EigMGM3,ZER,ZEL,kont)
     Y_l_mH3 = Matmul( Matmul(Conjg(ZER), Y_l_mH3), Transpose(Conjg(ZEL)))

     Y_T_mH3 = Y_T_0
     Y_Z_mH3 = Y_T_0
     Y_S_mH3 = Y_T_0
     A_l_mH3 = Matmul( Matmul(Conjg(ZER), A_l_mH3), Transpose(Conjg(ZEL)))
     M2_E_mH3 = Matmul( Matmul(ZER, M2_E_mH3), Transpose(Conjg(ZER)))
     M2_L_mH3 = Matmul( Matmul(ZEL, M2_L_mH3), Transpose(Conjg(ZEL)))
    End If

    Call ParametersToG5(gauge_mH3, Y_l_mH3, Y_T_mH3, Y_d_mH3, Y_u_mH3, Y_Z_mH3 &
        & , Y_S_mH3, lam12_0(1), lam12_0(2), Mi_mH3, A_l_mH3, A_T_mH3, A_d_mH3 &
        & , A_u_mH3, A_Z_mH3, A_S_mH3, Alam12_MH3(1), Alam12_MH3(2), M2_E_mH3  &
        & , M2_L_mH3, M2_D_mH3, M2_Q_mH3, M2_U_mH3, M2_H_mH3, M2_T_mH3         &
        & , M2_Z_mH3, M2_S_mH3, M15, M15, M15, mu_mH3, B_mH3, MnuL5, g1d )

    If (UseFixedGUTScale) Then
     tz = Log(m_lo/GUT_scale)
     mGUT = GUT_scale
     dt = - tz / 50._dp
     Call odeint(g1d, 356, tz, 0._dp, delta, dt, 0._dp, rge356, kont)
     If (kont.Ne.0) Then
      Iname = Iname -1
      Return
     End If

    Else

     If (g1d(1).Lt.g1d(2)) Then ! I am still below GUT scale
      tz = Log(m_lo/1.e18_dp)
      dt = - tz / 50._dp
      Call odeintB(g1d, 356, tz, 0._dp, delta, dt, 0._dp, rge356, t_out, kont) 

      If (kont.Eq.0) Then
       FoundUnification = .True.
       mGUT = 1.e18_dp * Exp(t_out)
       gGUT = Sqrt( 0.5_dp * (g1d(1)**2+g1d(2)**2) )
       g1d(1) = gGUT
       g1d(2) = gGUT
       If (StrictUnification) g1d(3) = gGUT
      Else
       Iname = Iname - 1
       Return
      End If

     Else If (g1d(1).Eq.g1d(2)) Then ! I am at the GUT scale, very unlikely
                                    ! but possible
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = g1d(1)
      If (StrictUnification) g1d(3) = gGUT

     End If
 
    End If

   End If ! check if below or above GUT scale

! Florian Staub Seesaw II+III
# ifdef SEESAWIII
  Else If (HighScaleModel.Eq.'SEESAW_III_3G') Then

   If (g1(1).Gt.g1(2)) Then ! already above the GUT scale
    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
    Call Set_Decoupling_Heavy_States(.True.)

    If (UseFixedGUTScale) Then
     tz = Log(m_lo/GUT_scale)
     mGUT = GUT_scale
     dt = - tz / 50._dp
     Call odeint(g1, 213, tz, 0._dp, delta, dt, 0._dp, rge213, kont)
     If (kont.Ne.0) Then
      Iname = Iname -1
      Return
     End If

    Else ! .not.UseFixedGUTScale
     tz = Log(m_lo/1.e15_dp)
     dt = - tz / 50._dp
     Call odeintC(g1, 213, tz, 0._dp, delta, dt, 0._dp, rge213, t_out, kont)
     If (kont.Eq.0) Then
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1(1)**2+g1(2)**2) )
      g1(1) = gGUT
      g1(2) = gGUT
      If (StrictUnification) g1(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If

    End If ! UseFixedGUTScale

    Call GToParameters(g1, gauge_H24, Ye_H24, Yd_H24, Yu_H24, Mi_H24, AYe_H24 &
      & , AYd_H24, AYu_H24, Me2_H24, Ml2_H24, Md2_H24, Mq2_H24, Mu2_H24       &
      & , Mh2_H24, mu_H24, Amue)

    If (i_in.Eq.1) Then
     Yb30_H24(3,:,:) = ZeroC
     Yw30_H24(3,:,:) = ZeroC
     Yx30_H24(3,:,:) = ZeroC
     AYx3_H24 = ZeroC
     AYb3_H24 = ZeroC
     AYw3_H24 = ZeroC
     AMBM3 = ZeroC
     AMGM3 = ZeroC
     AMWM3 = ZeroC
     AMXM3 = ZeroC
     mHw32 = ZeroC
     mHg32 = ZeroC
     mHb32 = ZeroC
     mHx32 = ZeroC
     mHxb32 = ZeroC
     MnuL5 = ZeroC
    End If

    Call ParametersToG555(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24  &
      & , Ye_H24, Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:), mu_H24, MXM3 &
      & , MWM3, MGM3, MBM3, AYu_H24, AYd_H24, AYe_H24, AYb3_H24, AYw3_H24       &
      & , AYx3_H24, Amue, AMXM3, AMWM3, AMGM3, AMBM3, mq2_H24, ml2_H24          &
      & , mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24, mHw32, mHg32       &
      & , mHb32, mHx32, mHxb32, Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5, g1f)

    NGHb3 = 3._dp
    NGHg3 = 3._dp 
    NGHw3 = 3._dp 
    NGHx3 = 3._dp 
    NGHxb3 = 3._dp 
    ThresholdCrossed = 3

   Else

    Call GToParameters(g1, gauge_H24, Ye_H24, Yd_H24, Yu_H24, Mi_H24, AYe_H24 &
      & , AYd_H24, AYu_H24, Me2_H24, Ml2_H24, Md2_H24, Mq2_H24, Mu2_H24       &
      & , Mh2_H24, mu_H24, Amue)

   !-----------------------------------------------------
   ! 1 -> 2
   ! adding shifts to gauge couplings
   !-----------------------------------------------------
    If (TwoLoopRGE) Then

     gauge_h24(1) = gauge_h24(1) * (1._dp - oo16pi2 * gauge_h24(1)**2           &
                  &                  * 5._dp/2._dp*Log(MassMXM3(1)/MWM30(1,1)) )
     gauge_h24(2) = gauge_h24(2) * (1._dp - oo16pi2 * gauge_h24(2)**2           &
                &                     *( 1.5_dp *Log(MassMXM3(1)/MWM30(1,1))    &
                &                      + 2._dp *Log(MassMWM3(1)/MWM30(1,1))) )
     gauge_h24(3) = gauge_h24(3) * (1._dp - oo16pi2 * gauge_h24(3)**2           &
                &                       * (Log(MassMXM3(1)/MWM30(1,1))          &
                &                         + 3._dp*Log(MassMGM3(1)/MWM30(1,1)) ) )
    End If

    MWM3 = MWM3running(1,:,:)

    NGHb3 = 1._dp
    NGHg3 = 1._dp 
    NGHw3 = 1._dp 
    NGHx3 = 1._dp 
    NGHxb3 = 1._dp 
    ThresholdCrossed = 1 

    If (i_in.Eq.1) Then ! a first initialisation, later take saved values
     AMBM3 = 0._dp
     AMGM3 = 0._dp
     AMWM3 = 0._dp
     AMXM3 = 0._dp
     mHw32_s = 0._dp
     mHg32_s = 0._dp
     mHb32_s = 0._dp
     mHx32_s = 0._dp
     mHxb32_s = 0._dp
     AYb3_s = 0._dp
     AYw3_s = 0._dp
     AYx3_s = 0._dp
    End If

    Call ParametersToG555(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24  &
     & , Ye_H24, Yb30_H24(1,:,:),Yw30_H24(1,:,:),Yx30_H24(1,:,:), mu_H24, MXM3  &
     & , MWM3, MGM3, MBM3, AYu_H24, AYd_H24, AYe_H24, AYb3_s(1,:,:)             &
     & , AYw3_s(1,:,:), AYx3_s(1,:,:), Amue, AMXM3, AMWM3, AMGM3, AMBM3         &
     & , mq2_H24, ml2_H24, mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24    &
     & , mHw32_s(1,:,:), mHg32_s(1,:,:), mHb32_s(1,:,:), mHx32_s(1,:,:)         &
     & , mHxb32_s(1,:,:), Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5, g1f)

    NGHb3 = 3._dp

    m_lo = MWM30(1,1)
    If (MWM30(1,1).Ne.MWM30(2,2)) Then
     m_hi = MWM30(2,2)
     tz = Log(m_lo / m_hi)
     dt = - tz / 50._dp  
     Call odeint(g2f, 555, tz, 0._dp, delta, dt, 0._dp, rge555, kont)

     Call GToParameters555(g1f, gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24  &
      & , Yd_H24, Ye_H24, Yb30_H24(1,:,:),Yw30_H24(1,:,:),Yx30_H24(1,:,:)      &
      & , mu_H24, MXM3, MWM3, MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24, AYb3_s(1,:,:) &
      & , AYw3_s(1,:,:), AYx3_s(1,:,:), Amue, AMXM3, AMWM3, AMGM3, AMBM3       &
      & , mq2_H24, ml2_H24, mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24  &
      & , mHw32_s(1,:,:), mHg32_s(1,:,:), mHb32_s(1,:,:), mHx32_s(1,:,:)       &
      & , mHxb32_s(1,:,:), Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5)

     m_lo = m_hi
    End If

    Call ParametersToG555(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24        &
      & , Yd_H24, Ye_H24, Yb30_H24(1,:,:),Yw30_H24(1,:,:),Yx30_H24(1,:,:)      &
      & , mu_H24, MXM3, MWM3, MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24, AYb3_s(1,:,:) &
      & , AYw3_s(1,:,:), AYx3_s(1,:,:), Amue, AMXM3, AMWM3, AMGM3, AMBM3       &
      & , mq2_H24, ml2_H24, mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24  &
      & , mHw32_s(1,:,:), mHg32_s(1,:,:), mHb32_s(1,:,:), mHx32_s(1,:,:)       &
      & , mHxb32_s(1,:,:), Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5, g1f)

!    Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24 &
!       &  ,Ye_H24,Yb30_H24(2,:,:),Yw30_H24(2,:,:),Yx30_H24(2,:,:),g1f)

    !---------
    ! 2 -> 3
    !---------
    If (g1f(1).Gt.g1f(2)) Then ! already above the GUT scale
    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
     Call Set_Decoupling_Heavy_States(.True.)

     tz = Log(m_lo/1.e15_dp)
     dt = - tz / 50._dp
     Call odeintC(g1f, 555, tz, 0._dp, delta, dt, 0._dp, rge555, t_out, kont)
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1f(1)**2+g1f(2)**2) )
      g1f(1) = gGUT
      g1f(2) = gGUT
      If (StrictUnification) g1f(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If

     Call ParametersToG555(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24       &
      & , Yd_H24, Ye_H24, Yb30_H24(2,:,:),Yw30_H24(2,:,:),Yx30_H24(2,:,:)      &
      & , mu_H24, MXM3, MWM3, MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24, AYb3_s(2,:,:) &
      & , AYw3_s(2,:,:), AYx3_s(2,:,:), Amue, AMXM3, AMWM3, AMGM3, AMBM3       &
      & , mq2_H24, ml2_H24, mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24  &
      & , mHw32_s(2,:,:), mHg32_s(2,:,:), mHb32_s(2,:,:), mHx32_s(2,:,:)       &
      & , mHxb32_s(2,:,:), Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5, g1f)

     Call GToParameters555(g1f, gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24  &
      & , Yd_H24, Ye_H24, Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:)      &
      & , mu_H24, MXM3, MWM3, MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24, AYb3_s(3,:,:) &
      & , AYw3_s(3,:,:), AYx3_s(3,:,:), Amue, AMXM3, AMWM3, AMGM3, AMBM3       &
      & , mq2_H24, ml2_H24, mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24  &
      & , mHw32_s(3,:,:), mHg32_s(3,:,:), mHb32_s(3,:,:), mHx32_s(3,:,:)       &
      & , mHxb32_s(3,:,:), Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5)

!     Call GToParameters111(g1f,gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24   &
!       & ,Yd_H24,Ye_H24,Yb30_H24(2,:,:),Yw30_H24(2,:,:),Yx30_H24(2,:,:))
!    Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24 &
!       & ,Ye_H24,Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:),g1f)

     NGHb3 = 3._dp
     NGHg3 = 3._dp 
     NGHw3 = 3._dp 
     NGHx3 = 3._dp 
     NGHxb3 = 3._dp 
     ThresholdCrossed = 3 

    Else   
     MWM3 = MWM3running(2,:,:)
    !-----------------------------------------------------
    ! adding shifts to gauge couplings
    !-----------------------------------------------------
     If (TwoLoopRGE) Then
      gauge_h24(1) = gauge_h24(1) * (1._dp - oo16pi2 * gauge_h24(1)**2         &
                &                     * 5._dp/2._dp*Log(MassMXM3(2)/MWM30(2,2)) )
      gauge_h24(2) = gauge_h24(2) * (1._dp - oo16pi2 * gauge_h24(2)**2         &
                &                       *( 1.5_dp *Log(MassMXM3(2)/MWM30(2,2)) &
                &                   + 2._dp *Log(MassMWM3(2)/MWM30(2,2))) )
      gauge_h24(3) = gauge_h24(3) * (1._dp - oo16pi2 * gauge_h24(3)**2         &
                &                       * (Log(MassMXM3(2)/MWM30(2,2)) &
                &                         + 3._dp*Log(MassMGM3(2)/MWM30(2,2)) ) )

     End If

     Call ParametersToG555(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24       &
      & , Yd_H24, Ye_H24, Yb30_H24(2,:,:),Yw30_H24(2,:,:),Yx30_H24(2,:,:)      &
      & , mu_H24, MXM3, MWM3, MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24, AYb3_s(2,:,:) &
      & , AYw3_s(2,:,:), AYx3_s(2,:,:), Amue, AMXM3, AMWM3, AMGM3, AMBM3       &
      & , mq2_H24, ml2_H24, mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24  &
      & , mHw32_s(2,:,:), mHg32_s(2,:,:), mHb32_s(2,:,:), mHx32_s(2,:,:)       &
      & , mHxb32_s(2,:,:), Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5, g1f)

!     Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24&
!        &  ,Ye_H24,Yb30_H24(2,:,:),Yw30_H24(2,:,:),Yx30_H24(2,:,:),g1f)

     NGHb3 = 2._dp
     NGHg3 = 2._dp 
     NGHw3 = 2._dp 
     NGHx3 = 2._dp 
     NGHxb3 = 2._dp 
     ThresholdCrossed = 2 

     m_lo = MWM30(2,2)
     If (MWM30(2,2).Ne.MWM30(3,3)) Then
      m_hi = MWM30(3,3)
      tz = Log(m_lo / m_hi)
      dt = - tz / 50._dp  
      Call odeint(g1f, 555, tz, 0._dp, delta, dt, 0._dp, rge555, kont)

     Call GToParameters555(g1f, gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24  &
      & , Yd_H24, Ye_H24, Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:)      &
      & , mu_H24, MXM3, MWM3, MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24, AYb3_s(3,:,:) &
      & , AYw3_s(3,:,:), AYx3_s(3,:,:), Amue, AMXM3, AMWM3, AMGM3, AMBM3       &
      & , mq2_H24, ml2_H24, mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24  &
      & , mHw32_s(3,:,:), mHg32_s(3,:,:), mHb32_s(3,:,:), mHx32_s(3,:,:)       &
      & , mHxb32_s(3,:,:), Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5)

!      Call GToParameters111(g1f,gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24  &
!         &  ,Yd_H24,Ye_H24,Yb30_H24(2,:,:),Yw30_H24(2,:,:),Yx30_H24(2,:,:))

      m_lo = m_hi
     End If

     If (Ynu_at_MR3) Then
      Call FermionMass(Ye_H24,1._dp,EigMGM3,ZER,ZEL,kont)
      Ye_H24 = Matmul( Matmul(Conjg(ZER), Ye_H24), Transpose(Conjg(ZEL)))

      Yb30_H24(3,:,:) = Yb3_H24_GUT
      Yw30_H24(3,:,:) = Yb3_H24_GUT
      Yx30_H24(3,:,:) = Yb3_H24_GUT
      AYe_H24 = Matmul( Matmul(Conjg(ZER), AYe_H24), Transpose(Conjg(ZEL)))
      ME2_H24 = Matmul( Matmul(ZER, ME2_H24), Transpose(Conjg(ZER)))
      ML2_H24 = Matmul( Matmul(ZEL, ML2_H24), Transpose(Conjg(ZEL)))
     End If

     Call ParametersToG555(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24       &
      & , Yd_H24, Ye_H24, Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:)      &
      & , mu_H24, MXM3, MWM3, MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24, AYb3_s(3,:,:) &
      & , AYw3_s(3,:,:), AYx3_s(3,:,:), Amue, AMXM3, AMWM3, AMGM3, AMBM3       &
      & , mq2_H24, ml2_H24, mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24  &
      & , mHw32_s(3,:,:), mHg32_s(3,:,:), mHb32_s(3,:,:), mHx32_s(3,:,:)       &
      & , mHxb32_s(3,:,:), Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5, g1f)

!     Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24&
!         &  ,Ye_H24,Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:),g1f)

     If (g1f(1).Gt.g1f(2)) Then ! already above the GUT scale
    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
      Call Set_Decoupling_Heavy_States(.True.)

      tz = Log(m_lo/1.e15_dp)
      dt = - tz / 50._dp
      Call odeintC(g1f, 555, tz, 0._dp, delta, dt, 0._dp, rge555, t_out, kont)
      If (kont.Eq.0) Then
       FoundUnification = .True.
       mGUT = 1.e15_dp * Exp(t_out)
       gGUT = Sqrt( 0.5_dp * (g1f(1)**2+g1f(2)**2) )
       g1f(1) = gGUT
       g1f(2) = gGUT
       If (StrictUnification) g1f(3) = gGUT
      Else
       Iname = Iname - 1
       Return
      End If

     Call ParametersToG555(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24       &
      & , Yd_H24, Ye_H24, Yb30_H24(2,:,:),Yw30_H24(2,:,:),Yx30_H24(2,:,:)      &
      & , mu_H24, MXM3, MWM3, MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24, AYb3_s(2,:,:) &
      & , AYw3_s(2,:,:), AYx3_s(2,:,:), Amue, AMXM3, AMWM3, AMGM3, AMBM3       &
      & , mq2_H24, ml2_H24, mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24  &
      & , mHw32_s(2,:,:), mHg32_s(2,:,:), mHb32_s(2,:,:), mHx32_s(2,:,:)       &
      & , mHxb32_s(2,:,:), Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5, g1f)

     Call GToParameters555(g1f, gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24  &
      & , Yd_H24, Ye_H24, Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:)      &
      & , mu_H24, MXM3, MWM3, MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24, AYb3_s(3,:,:) &
      & , AYw3_s(3,:,:), AYx3_s(3,:,:), Amue, AMXM3, AMWM3, AMGM3, AMBM3       &
      & , mq2_H24, ml2_H24, mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24  &
      & , mHw32_s(3,:,:), mHg32_s(3,:,:), mHb32_s(3,:,:), mHx32_s(3,:,:)       &
      & , mHxb32_s(3,:,:), Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5)

!      Call GToParameters111(g1f,gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24  &
!        &  ,Yd_H24,Ye_H24,Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:))
!
!     Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24,Yd_H24&
!       & ,Ye_H24,Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:),g1f)


      NGHb3 = 3._dp
      NGHg3 = 3._dp 
      NGHw3 = 3._dp 
      NGHx3 = 3._dp 
      NGHxb3 = 3._dp 
      ThresholdCrossed = 3 

     Else

     !-----------------------------------------------------
     ! 3 -> GUT
     !-----------------------------------------------------

      MWM3 = MWM3running(3,:,:)

      NGHb3 = 3._dp
      NGHg3 = 3._dp 
      NGHw3 = 3._dp 
      NGHx3 = 3._dp 
      NGHxb3 = 3._dp 
      ThresholdCrossed = 3 

     !-----------------------------------------------------
     ! adding shifts to gauge couplings
     !-----------------------------------------------------
      If (TwoLoopRGE) Then
       gauge_h24(1) = gauge_h24(1) * (1._dp - oo16pi2 * gauge_h24(1)**2        &
                &                     * 5._dp/2._dp*Log(MassMXM3(3)/MWM30(3,3)) )
       gauge_h24(2) = gauge_h24(2) * (1._dp - oo16pi2 * gauge_h24(2)**2        &
                &                       *( 1.5_dp *Log(MassMXM3(3)/MWM30(3,3)) &
                &                    + 2._dp *Log(MassMWM3(3)/MWM30(3,3))) )
       gauge_h24(3) = gauge_h24(3) * (1._dp - oo16pi2 * gauge_h24(3)**2        &
                &                       * (Log(MassMXM3(3)/MWM30(3,3)) &
                &                         + 3._dp*Log(MassMGM3(3)/MWM30(3,3)) ) )
      
      End If

     Call ParametersToG555(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24       &
      & , Yd_H24, Ye_H24, Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:)      &
      & , mu_H24, MXM3, MWM3, MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24, AYb3_s(3,:,:) &
      & , AYw3_s(3,:,:), AYx3_s(3,:,:), Amue, AMXM3, AMWM3, AMGM3, AMBM3       &
      & , mq2_H24, ml2_H24, mH2_H24(1), mH2_H24(2), md2_H24, mu2_H24, me2_H24  &
      & , mHw32_s(3,:,:), mHg32_s(3,:,:), mHb32_s(3,:,:), mHx32_s(3,:,:)       &
      & , mHxb32_s(3,:,:), Mi_H24(1), Mi_H24(2), Mi_H24(3), MnuL5, g1f)

!      Call ParametersToG111(gauge_H24(1),gauge_H24(2),gauge_H24(3),Yu_H24     &
!      & ,Yd_H24,Ye_H24,Yb30_H24(3,:,:),Yw30_H24(3,:,:),Yx30_H24(3,:,:),g1f)

      If (UseFixedGUTScale) Then
       tz = Log(m_lo/GUT_scale)
       mGUT = GUT_scale
       dt = - tz / 50._dp
       Call odeint(g1f, 555, tz, 0._dp, delta, dt, 0._dp, rge555, kont)
       If (kont.Ne.0) Then
        Iname = Iname -1
        Return
       End If

      Else

       If (g1f(1).Lt.g1f(2)) Then ! I am still below GUT scale
        tz = Log(m_lo/1.e18_dp)
        dt = - tz / 50._dp
        Call odeintB(g1f, 555, tz, 0._dp, delta, dt, 0._dp, rge555, t_out, kont) 
        If (kont.Eq.0) Then
         FoundUnification = .True.
         mGUT = 1.e18_dp * Exp(t_out)
         gGUT = Sqrt( 0.5_dp * (g1f(1)**2+g1f(2)**2) )
         g1f(1) = gGUT
         g1f(2) = gGUT
         If (StrictUnification) g1f(3) = gGUT
        Else
         Iname = Iname - 1
         Return
        End If

       Else If (g1f(1).Eq.g1f(2)) Then ! I am at the GUT scale, very unlikely
                                    ! but possible
        FoundUnification = .True.
        mGUT = 1.e15_dp * Exp(t_out)
        gGUT = g1f(1)
        If (StrictUnification) g1f(3) = gGUT

       Else ! I have already crossed the GUT scale
        tz = Log(m_lo/1.e15_dp)
        dt = - tz / 50._dp
        Call odeintC(g1f, 555, tz, 0._dp, delta, dt, 0._dp, rge555, t_out, kont)
        If (kont.Eq.0) Then
         FoundUnification = .True.
         mGUT = 1.e15_dp * Exp(t_out)
         gGUT = Sqrt( 0.5_dp * (g1f(1)**2+g1f(2)**2) )
         g1f(1:2) = gGUT
         If (StrictUnification) g1f(3) = gGUT
        Else
         Iname = Iname - 1
         Return
        End If
       End If

      End If

     End If
    End If


   End If


  Else If (HighScaleModel.Eq.'SEESAW_II_SARAH') Then

   If (g1(1).Gt.g1(2)) Then ! already above the GUT scale
    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
    If (UseFixedGUTScale) Then
     tz = Log(m_lo/GUT_scale)
     mGUT = GUT_scale
     dt = - tz / 50._dp
     Call odeint(g1, 213, tz, 0._dp, delta, dt, 0._dp, rge213, kont)
     If (kont.Ne.0) Then
      Iname = Iname -1
      Return
     End If

    Else ! .not.UseFixedGUTScale
     tz = Log(m_lo/1.e15_dp)
     dt = - tz / 50._dp
     Call odeintC(g1, 213, tz, 0._dp, delta, dt, 0._dp, rge213, t_out, kont)
     If (kont.Eq.0) Then
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = Sqrt( 0.5_dp * (g1(1)**2+g1(2)**2) )
      g1(1) = gGUT
      g1(2) = gGUT
      If (StrictUnification) g1(3) = gGUT
     Else
      Iname = Iname - 1
      Return
     End If

    End If ! UseFixedGUTScale

    Call GToParameters(g1, gauge_H15, Ye_H15, Yd_H15, Yu_H15, Mi_H15, AYe_H15 &
      & , AYd_H15, AYu_H15, Me2_H15, Ml2_H15, Md2_H15, Mq2_H15, Mu2_H15       &
      & , Mh2_H15, mu_H15, Amue)

    g1_h15 = gauge_H15(1)
    g2_h15 = gauge_H15(2)
    g3_h15 = gauge_H15(3)
    MassB_h15 = Mi_H15(1)
    MassWB_h15 = Mi_H15(2)
    MassG_h15 = Mi_H15(3)

    Call ParametersToG365(g1_h15,g2_h15,g3_h15,Yu_h15,Yd_h15,Ye_h15,Yt0_h15     &
       & ,Ys0_h15,Yz0_h15, Lambda10,Lambda20,mu_h15,MTM,MZM,MSM,AYu_h15,AYd_h15 &
       & ,AYe_h15,AYt_h15,AYs_h15,AYz_h15,ALambda1,ALambda2,Amue,AMTM,AMZM,AMSM &
       & ,mq2_h15,ml2_h15,mHd2_h15, mHu2_h15,md2_h15,mu2_h15,me2_h15,mt2_H15    &
       & ,mtb2_H15,ms2_H15,msb2_H15,mz2_H15,mzb2_H15,MassB_h15,MassWB_h15       &
       & ,MassG_h15,MnuL5,g1g)

!    Call GToCouplings(g1, gauge_h15, Ye_h15, Yd_h15, Yu_h15)
!    M15 = M_H3(1)

!    Call ParametersToG115(gauge_h15(1),gauge_h15(2),gauge_h15(3), &
!    & Yu_h15,Yd_h15,Ye_h15,Yt0_h15,Ys0_h15,Yz0_h15,Lambda10,Lambda20,g1g)
 
   Else ! still below the GUT scale


    !---------------------------------------------------
    ! heavy states are to be included
    !---------------------------------------------------
    Call GToParameters(g1, gauge_H15, Ye_H15, Yd_H15, Yu_H15, Mi_H15, AYe_H15 &
      & , AYd_H15, AYu_H15, Me2_H15, Ml2_H15, Md2_H15, Mq2_H15, Mu2_H15       &
      & , Mh2_H15, mu_H15, Amue)

    If (TwoLoopRGE) Then
     !-----------------------------------------------------
     ! adding shifts to gauge couplings
     !-----------------------------------------------------
     gauge_h15(1) = gauge_h15(1) * (1._dp - 0.2_dp * oo16pi2 * gauge_h15(1)**2 &
                 &                       * (8._dp * Log(MSM/MTM_GUT)           &
                 &                         + 9._dp * Log(MTM0/MTM_GUT)         &
                 &                         + 0.5_dp*Log(MZM/MTM_GUT) ) )
     gauge_h15(2) = gauge_h15(2) * (1._dp - oo16pi2 * gauge_h15(2)**2          &
                 &    * (2._dp *Log(MTM0/MTM_GUT) + 1.5_dp *Log(MZM/MTM_GUT) ))
     gauge_h15(3) = gauge_h15(3) * (1._dp - oo16pi2 * gauge_h15(3)**2          &
                 &    * (2.5_dp*Log(MSM/MTM_GUT) + Log(MZM/MTM_GUT) ) )
    End If


    If (i_in.Eq.1) Then
     Lambda10 = ZeroC
     Lambda20 = ZeroC
     ALambda1 = ZeroC
     ALambda2 = ZeroC
     AMTM = ZeroC
     AMZM = ZeroC
     AMSM = ZeroC
     mt2_H15 = ZeroC
     mtb2_H15 = ZeroC
     ms2_H15 = ZeroC
     msb2_H15 = ZeroC
     mz2_H15 = ZeroC
     mzb2_H15 = ZeroC
    End If

    If (Ynu_at_MR3) Then
     Call FermionMass(Ye_H15,1._dp,EigMGM3,ZER,ZEL,kont)
     Ye_H15 = Matmul( Matmul(Conjg(ZER), Ye_H15), Transpose(Conjg(ZEL)))

     YT0_H15 = Y_T_0
     Ys0_H15 = Y_T_0
     Yz0_H15 = Y_T_0
     AYe_H15 = Matmul( Matmul(Conjg(ZER), AYe_H15), Transpose(Conjg(ZEL)))
     ME2_H15 = Matmul( Matmul(ZER, ME2_H15), Transpose(Conjg(ZER)))
     ML2_H15 = Matmul( Matmul(ZEL, ML2_H15), Transpose(Conjg(ZEL)))
    End If

    g1_h15 = gauge_H15(1)
    g2_h15 = gauge_H15(2)
    g3_h15 = gauge_H15(3)
    MassB_h15 = Mi_H15(1)
    MassWB_h15 = Mi_H15(2)
    MassG_h15 = Mi_H15(3)

    Call ParametersToG365(g1_h15,g2_h15,g3_h15,Yu_h15,Yd_h15,Ye_h15,Yt0_h15     &
       & ,Ys0_h15,Yz0_h15, Lambda10,Lambda20,mu_h15,MTM,MZM,MSM,AYu_h15,AYd_h15 &
       & ,AYe_h15,AYt_h15,AYs_h15,AYz_h15,ALambda1,ALambda2,Amue,AMTM,AMZM,AMSM &
       & ,mq2_h15,ml2_h15,mHd2_h15, mHu2_h15,md2_h15,mu2_h15,me2_h15,mt2_H15    &
       & ,mtb2_H15,ms2_H15,msb2_H15,mz2_H15,mzb2_H15,MassB_h15,MassWB_h15       &
       & ,MassG_h15,MnuL5,g1g)

!    Call ParametersToG115(gauge_h15(1),gauge_h15(2),gauge_h15(3),Yu_h15,Yd_h15 &
!       & ,Ye_h15,Yt0_h15,Ys0_h15,Yz0_h15,Lambda10,Lambda20,g1g)

    If (UseFixedGUTScale) Then
     tz = Log(m_lo/GUT_scale)
     mGUT = GUT_scale
     dt = - tz / 50._dp
     Call odeint(g1g, 365, tz, 0._dp, delta, dt, 0._dp, rge365, kont)
     If (kont.Ne.0) Then
      Iname = Iname -1
      Return
     End If

    Else

     If (g1g(1).Lt.g1g(2)) Then ! I am still below GUT scale
      tz = Log(m_lo/1.e18_dp)
      dt = - tz / 50._dp
      Call odeintB(g1g, 365, tz, 0._dp, delta, dt, 0._dp, rge365, t_out, kont) 

      If (kont.Eq.0) Then
       FoundUnification = .True.
       mGUT = 1.e18_dp * Exp(t_out)
       gGUT = Sqrt( 0.5_dp * (g1g(1)**2+g1g(2)**2) )
       g1g(1) = gGUT
       g1g(2) = gGUT
       If (StrictUnification) g1g(3) = gGUT
      Else
       Iname = Iname - 1
       Return
      End If

     Else If (g1g(1).Eq.g1g(2)) Then ! I am at the GUT scale, very unlikely
                                    ! but possible
      FoundUnification = .True.
      mGUT = 1.e15_dp * Exp(t_out)
      gGUT = g1g(1)
      If (StrictUnification) g1g(3) = gGUT

     End If
 
    End If


   End If ! check if below or above GUT scale
# endif SEESAWIII
! Florian Staub Seesaw II+III

  End If

  If ((.Not.UseFixedGUTScale).And.(.Not.FoundUnification)) Then
   Write (ErrCan,*) 'SUGRA: no unification found'
   SugraErrors(1) = .True.
   kont = -417
   Call AddError(417)
   Iname = Iname - 1
   Return
  End If

  !------------------------------------
  ! Saving parameters at high scale
  !------------------------------------
  mGUT_Save = mGUT

  !---------------------------------------
  ! boundary condition at the high scale
  !---------------------------------------
  If (HighScaleModel.Eq.'SUGRA_NuR') Then
   Call BoundaryHS2(g1a,g2a)

  Else If (HighScaleModel.Eq.'SUGRA_NuR1') Then
   Call BoundaryHS2(g1b,g2b)

  Else If ((HighScaleModel.Eq.'SEESAW_II').And.Fifteen_plet) Then
    Call BoundaryHS2(g1d,g2d)

! Florian Staub Seesaw II+III
# ifdef SEESAWIII
 Else If (HighScaleModel.Eq.'SEESAW_III_3G') Then
   Call BoundaryHS2(g1f,g2f)

  Else If (HighScaleModel.Eq.'SEESAW_II_SARAH') Then
   Call BoundaryHS2(g1g,g2g)
# endif SEESAWIII
! Florian Staub Seesaw II+III
  Else
   Call BoundaryHS2(g1,g2)
  End If
  
  !--------------------------------------
  ! running down to the electroweak scale
  !--------------------------------------
  If (HighScaleModel.Eq.'SUGRA_NuR') Then

   !------------------------
   ! m_GUT -> m_nuR_3
   !------------------------
   m_hi = mGUT
   m_lo = MNuR(3)
   If (Abs(m_lo).Lt.Abs(m_hi)) Then
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)
    m_hi = m_lo
   Endif

   Call GToParameters3(g2a, gauge_mR, y_l_mR(3,:,:), y_nu_mR(3,:,:)          &
      & , y_d_mR(3,:,:), y_u_mR(3,:,:), Mi_mR, A_l_mR(3,:,:), A_nu_mR(3,:,:) &
      & , A_d_mR(3,:,:), A_u_mR(3,:,:), M2_E_mR(3,:,:), M2_L_mR(3,:,:)       &
      & , M2_R_mR(3,:,:), M2_D_mR(3,:,:), M2_Q_mR(3,:,:), M2_U_mR(3,:,:)     &
      & , M2_H_mR, mu_mR, B_mR, MnuL5)

   If (Ynu_at_MR3) Y_nu_mR(3,:,:) = Y_nu_0 ! to enhance numerical stability
                                           ! and speed

   Do i1=1,3
    Do i2=1,3
     MnuL5(i1,i2) = - Y_nu_mR(3,3,i1) * Y_nu_mR(3,3,i2) / MNuR(3)
    End Do
   End Do
   Y_nu_mR(2,:,:) = Y_nu_mR(3,:,:)
   Y_nu_mR(2,3,:) = 0._dp
   A_nu_mR(2,:,:) = A_nu_mR(3,:,:)
   A_nu_mR(2,3,:) = 0._dp
   M2_R_mR(2,:,:) = M2_R_mR(3,:,:)
   M2_R_mR(2,3,:) = 0._dp
   M2_R_mR(2,:,3) = 0._dp
   !------------------------
   ! m_nuR_3 -> m_nuR_2
   !------------------------
   Call ParametersToG3(gauge_mR, y_l_mR(3,:,:), y_nu_mR(2,:,:), y_d_mR(3,:,:) &
      & , y_u_mR(3,:,:), Mi_mR, A_l_mR(3,:,:), A_nu_mR(2,:,:), A_d_mR(3,:,:)  &
      & , A_u_mR(3,:,:), M2_E_mR(3,:,:), M2_L_mR(3,:,:), M2_R_mR(2,:,:)       &
      & , M2_D_mR(3,:,:), M2_Q_mR(3,:,:), M2_U_mR(3,:,:), M2_H_mR, mu_mR      &
      & , B_mR, MnuL5, g2a)

   m_lo = MNuR(2)
   If (Abs(m_lo).Lt.Abs(m_hi)) Then
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)
    m_hi = m_lo
   Endif

   Call GToParameters3(g2a, gauge_mR, y_l_mR(2,:,:), y_nu_mR(2,:,:)          &
      & , y_d_mR(2,:,:), y_u_mR(2,:,:), Mi_mR, A_l_mR(2,:,:), A_nu_mR(2,:,:) &
      & , A_d_mR(2,:,:), A_u_mR(2,:,:), M2_E_mR(2,:,:), M2_L_mR(2,:,:)       &
      & , M2_R_mR(2,:,:), M2_D_mR(2,:,:), M2_Q_mR(2,:,:), M2_U_mR(2,:,:)     &
      & , M2_H_mR, mu_mR, B_mR, MnuL5)

   Do i1=1,3
    Do i2=1,3
     MnuL5(i1,i2) = MnuL5(i1,i2) - Y_nu_mR(2,2,i1) * Y_nu_mR(2,2,i2) / MNuR(2)
    End Do
   End Do
   Y_nu_mR(1,:,:) = Y_nu_mR(2,:,:)
   Y_nu_mR(1,2,:) = 0._dp
   A_nu_mR(1,:,:) = A_nu_mR(2,:,:)
   A_nu_mR(1,2,:) = 0._dp
   M2_R_mR(1,:,:) = M2_R_mR(2,:,:)
   M2_R_mR(1,2,:) = 0._dp
   M2_R_mR(1,:,2) = 0._dp
   !------------------------
   ! m_nuR_2 -> m_nuR_1
   !------------------------
   Call ParametersToG3(gauge_mR, y_l_mR(2,:,:), y_nu_mR(1,:,:), y_d_mR(2,:,:) &
      & , y_u_mR(2,:,:), Mi_mR, A_l_mR(2,:,:), A_nu_mR(1,:,:), A_d_mR(2,:,:)  &
      & , A_u_mR(2,:,:), M2_E_mR(2,:,:), M2_L_mR(2,:,:), M2_R_mR(1,:,:)       &
      & , M2_D_mR(2,:,:), M2_Q_mR(2,:,:), M2_U_mR(2,:,:), M2_H_mR, mu_mR      &
      & , B_mR, MnuL5, g2a)

   m_lo = MNuR(1)
   If (Abs(m_lo).Lt.Abs(m_hi)) Then
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)
    m_hi = m_lo
   Endif

   Call GToParameters3(g2a, gauge_mR, y_l_mR(1,:,:), y_nu_mR(1,:,:)          &
      & , y_d_mR(1,:,:), y_u_mR(1,:,:), Mi_mR, A_l_mR(1,:,:), A_nu_mR(1,:,:) &
      & , A_d_mR(1,:,:), A_u_mR(1,:,:), M2_E_mR(1,:,:), M2_L_mR(1,:,:)       &
      & , M2_R_mR(1,:,:), M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:)     &
      & , M2_H_mR, mu_mR, B_mR, MnuL5)

   Do i1=1,3
    Do i2=1,3
     MnuL5(i1,i2) = MnuL5(i1,i2) - Y_nu_mR(1,1,i1) * Y_nu_mR(1,1,i2) / MNuR(1)
    End Do
   End Do
   !------------------------
   ! m_nuR_1 -> Q_EWSB
   !------------------------
   Call ParametersToG3(gauge_mR, y_l_mR(1,:,:), Zero33C, y_d_mR(1,:,:)     &
        & , y_u_mR(1,:,:), Mi_mR, A_l_mR(1,:,:), Zero33C, A_d_mR(1,:,:)    &
        & , A_u_mR(1,:,:), M2_E_mR(1,:,:), M2_L_mR(1,:,:), Zero33C         &
        & , M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:), M2_H_mR, mu_mR &
        & , B_mR, MnuL5, g2a)

   mudim = GetRenormalizationScale()
   mudim = Max(mudim, mZ2)
   tz = 0.5_dp * Log(m_hi**2/mudim)
   dt = - tz / 100._dp

   Call odeint(g2a, 285, tz, 0._dp, delta, dt, 0._dp, rge285, kont)

   Call GToParameters3(g2a, gauge, y_l, Y_nu, y_d, y_u, Mi, A_l, A_nu, A_d &
      & , A_u, M2_E, M2_L, M2_R, M2_D, M2_Q, M2_U, M2_H, mu, B, MnuL5)

   Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u, M2_E, M2_L  &
      & , M2_D, M2_Q, M2_U, M2_H, mu, B, g2)

  Else If (HighScaleModel.Eq.'SUGRA_NuR1') Then
   mudim = MNuR(1)**2
   tz = 0.5_dp * Log(mudim/mGUT**2)
   dt = tz / 50._dp
   Call odeint(g2a, 267, 0._dp, tz, delta, dt, 0._dp, rge267, kont)

   Call GToParameters2(g2a, gauge_mR, y_l_mR(1,:,:), y_nu_mR(1,:,:)          &
      & , y_d_mR(1,:,:), y_u_mR(1,:,:), Mi_mR, A_l_mR(1,:,:), A_nu_mR(1,:,:) &
      & , A_d_mR(1,:,:), A_u_mR(1,:,:), M2_E_mR(1,:,:), M2_L_mR(1,:,:)       &
      & , M2_R_mR(1,:,:), M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:)     &
      & , M2_H_mR, mu_mR, B_mR)

   Call ParametersToG(gauge_mR, y_l_mR(1,:,:), y_d_mR(1,:,:), y_u_mR(1,:,:)  &
      & , Mi_mR, A_l_mR(1,:,:), A_d_mR(1,:,:), A_u_mR(1,:,:), M2_E_mR(1,:,:) &
      & , M2_L_mR(1,:,:), M2_D_mR(1,:,:), M2_Q_mR(1,:,:), M2_U_mR(1,:,:)     &
      & , M2_H_mR, mu_mR, B_mR, g2)

   mudim = GetRenormalizationScale()
   mudim = Max(mudim, mZ2)

   tz = 0.5_dp * Log(mudim/mNuR(1)**2)
   dt = tz / 100._dp
   Call odeint(g2, 213, 0._dp, tz, delta, dt, 0._dp, rge213, kont)

  Else If ((HighScaleModel.Eq.'SEESAW_II').And.Fifteen_plet) Then

   If ( (oo4pi*Maxval(g2d(1:115)**2)).Gt.0.5_dp) Then
    Write(ErrCan,*) "Non perturbative regime at M_GUT"
    If (ErrorLevel.Ge.2) Call TerminateProgram
    Write(errcan,*) " "
    kont = -418
    Call AddError(418)
    Iname = Iname - 1
    Return
   End If
   
   !-------------------------
   ! run only if m_H < m_GUT
   !-------------------------
   If (m_H3(1).Lt.mGUT) Then
    !---------------------------------------------------
    ! heavy states are to be included
    !---------------------------------------------------
    Call Set_Decoupling_Heavy_States(.False.)

    mudim = M_H3(1)**2
    tz = 0.5_dp * Log(mudim/mGUT**2)
    dt = tz / 50._dp
    Call odeint(g2d, 356 , 0._dp, tz, delta, dt, 0._dp, rge356, kont)
    m_lo = M_H3(1)
    If ( (oo4pi*Maxval(g2d(1:115)**2)).Gt.0.5_dp) Then
     Write(ErrCan,*) "Non perturbative regime at M_H3"
     If (ErrorLevel.Ge.2) Call TerminateProgram
     Write(errcan,*) " "
     kont = -419
     Call AddError(419)
     Iname = Iname - 1
     Return
    End If

   Else
     m_lo = mGUT
   End If

   Call GToParameters5(g2d, gauge_mH3, y_l_mH3, y_T_mH3, y_d_mH3, y_u_mH3     &
          & , y_Z_mH3, y_S_mH3, lam12_mH3(1), lam12_mH3(2), Mi_mH3, A_l_mH3   &
          & , A_T_mH3, A_d_mH3, A_u_mH3, A_Z_mH3, A_S_mH3, Alam12_MH3(1)      &
          & , Alam12_MH3(2), M2_E_mH3, M2_L_mH3, M2_D_mH3, M2_Q_mH3, M2_U_mH3 &
          & , M2_H_mH3, M2_T_mH3, M2_Z_mH3, M2_S_mH3, MT15_mH3, MZ15_mH3      &
          & , MS15_mH3, mu_mH3, B_mH3, MnuL5)

   MnuL5 = - lam12_MH3(2) * y_T_mH3 / M_H3(1)

   Delta_b_1 = 0._dp ! decoupling the Higgs triplets
   Delta_b_2 = 0._dp ! decoupling the Higgs triplets

   If ((m_H3(1).Lt.mGUT).And.TwoLoopRGE) Then
   !-----------------------------------------------------
   ! adding shifts to gauge couplings
   !-----------------------------------------------------
    gauge_mH3(1) = gauge_mH3(1) / (1._dp - oo16pi2 * gauge_mH3(1)**2           &
                &                       * (8._dp/3._dp*Log(MS15_mH3/MT15_mH3) &
                &                         + Log(MZ15_mH3/MT15_mH3) /6._dp ) )
    gauge_mH3(2) = gauge_mH3(2) / (1._dp - oo16pi2 * gauge_mH3(2)**2           &
                &                       * 1.5_dp *Log(MZ15_mH3/MT15_mH3) )
    gauge_mH3(3) = gauge_mH3(3) / (1._dp - oo16pi2 * gauge_mH3(3)**2           &
                &                       * (2.5_dp*Log(MS15_mH3/MT15_mH3) &
                &                         + Log(MZ15_mH3/MT15_mH3) ) )
    Mi_mH3(1) = Mi_mH3(1) / (1._dp - oo16pi2 * gauge_mH3(1)**2           &
                &                       * (8._dp/3._dp*Log(MS15_mH3/MT15_mH3) &
                &                         + Log(MZ15_mH3/MT15_mH3) /6._dp ) )
    Mi_mH3(2) = Mi_mH3(2) / (1._dp - oo16pi2 * gauge_mH3(2)**2           &
                &                       * 1.5_dp *Log(MZ15_mH3/MT15_mH3) )
    Mi_mH3(3) = Mi_mH3(3) / (1._dp - oo16pi2 * gauge_mH3(3)**2           &
                &                       * (2.5_dp*Log(MS15_mH3/MT15_mH3) &
                &                         + Log(MZ15_mH3/MT15_mH3) ) )
   End If

    !---------------------------------------------------
    ! heavy states are not to be included
    !---------------------------------------------------
    Call Set_Decoupling_Heavy_States(.True.)

   Call ParametersToG4(gauge_mH3, y_l_mH3, Zero33C, y_d_mH3, y_u_mH3   &
          & , ZeroC, ZeroC, Mi_mH3, A_l_mH3, Zero33C, A_d_mH3 &
          & , A_u_mH3, ZeroC, ZeroC, M2_E_mH3, M2_L_mH3     &
          & , M2_D_mH3, M2_Q_mH3, M2_U_mH3, M2_H_mH3, ZeroR2, mu_mH3      &
          & , B_mH3, MnuL5, g2c)

   mudim = GetRenormalizationScale()
   mudim = Max(mudim, mZ2)

   tz = 0.5_dp * Log(mudim/m_lo**2)
   dt = tz / 100._dp
   Call odeint(g2c, 277, 0._dp, tz, delta, dt, 0._dp, rge277, kont)

   Call GToParameters4(g2c, gauge, y_l, y_T, y_d, y_u, lam12(1), lam12(2), Mi &
          & , A_l, A_T, A_d, A_u, Alam12(1), Alam12(2), M2_E, M2_L, M2_D      &
          & , M2_Q, M2_U, M2_H, M2_T, mu, B, MnuL5)
   Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u, M2_E, M2_L  &
      & , M2_D, M2_Q, M2_U, M2_H, mu, B, g2)

! Florian Staub Seesaw II+III
# ifdef SEESAWIII
  Else If (HighScaleModel.Eq.'SEESAW_III_3G') Then 
   !---------------
   ! GUT -> 3
   !---------------
   m_hi = mGUT
   If (Abs(MWM30(3,3)).Lt.Abs(m_hi)) Then
    m_lo = Abs(MWM30(3,3))
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2f, 573, tz, 0._dp, delta, dt, 0._dp, rge555, kont)
    m_hi = m_lo
   Endif

   Call GToParameters555(g2f,g1_H24,g2_H24,g3_H24,Yu_H24,Yd_H24,Ye_H24,Yb3_H24 &
      & ,Yw3_H24,Yx3_H24,mu_H24,MXM3,MWM3,MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24    &
      & ,AYb3_H24,AYw3_H24,AYx3_H24,Amue,AMXM3,AMWM3,AMGM3,AMBM3, mq2_H24      &
      & ,ml2_H24,mHd2_H24,mHu2_H24,md2_H24,mu2_H24,me2_H24,mHw32,mHg32,mHb32   &
      & ,mHx32,mHxb32,MassB_H24,MassWB_H24,MassG_H24,MnuL5)

   MBM3Running(3,:,:)=MBM3
   MGM3Running(3,:,:)=MGM3
   MWM3Running(3,:,:)=MWM3
   MXM3Running(3,:,:)=MXM3

   NGHb3 = 2._dp
   NGHg3 = 2._dp 
   NGHw3 = 2._dp 
   NGHx3 = 2._dp 
   NGHxb3 = 2._dp 
   ThresholdCrossed = 2

   Call FermionMass(MWM3,sqrt2,EigMWM3,RotWL, RotWR,kont)  
   Call FermionMass(MXM3,sqrt2,EigMXM3,RotXL, RotXR,kont)
   Call FermionMass(MGM3,sqrt2,EigMGM3,RotGL, RotGR,kont)
   Call FermionMass(MBM3,sqrt2,EigMBM3,RotBL, RotBR,kont)

   Yb3_H24 = Matmul(RotBL,Yb3_H24)
   Yw3_H24 = Matmul(RotWL,Yw3_H24)
   Yx3_H24 = Matmul(RotXL,Yx3_H24)
   AYb3_H24 = Matmul(RotBL,AYb3_H24)
   AYw3_H24 = Matmul(RotwL,AYb3_H24)
   AYx3_H24 = Matmul(RotxL,AYb3_H24)
   MWM3 = Matmul(Conjg(Transpose(RotWR)),Matmul(MWM3,RotWL))
   MGM3 = Matmul(Conjg(Transpose(RotGR)),Matmul(MGM3,RotGL))
   MXM3 = Matmul(Conjg(Transpose(RotXR)),Matmul(MXM3,RotXL))
   MBM3 = Matmul(Conjg(Transpose(RotBR)),Matmul(MBM3,RotBL))
   mHx32 = Matmul(Conjg(Transpose(RotXL)),Matmul(mHx32,RotXL)) 
   mHxb32 = Matmul(Conjg(Transpose(RotXR)),Matmul(mHxb32,RotXR))  
   mHg32 = Matmul(Conjg(Transpose(RotGL)),Matmul(mHg32,RotGL)) 
   mHb32 = Matmul(Conjg(Transpose(RotBL)),Matmul(mHb32,RotBL))  
   mHw32 = Matmul(Conjg(Transpose(RotWL)),Matmul(mHw32,RotWL))

   MassMWM3(3) = Maxval(EigMWM3)
   MassMXM3(3) = Maxval(EigMXM3)
   MassMGM3(3) = Maxval(EigMGM3)
   MassMBM3(3) = Maxval(EigMBM3)

   If (MassMWM3(3).Lt.Abs(mGUT)) Then
    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = 0.5_dp*Yw3_H24(3,i1)*Yw3_H24(3,i2)/MassMWM3(3) + &
                   & 0.3_dp*Yb3_H24(3,i1)*Yb3_H24(3,i2)/MassMBM3(3) 
     End Do
    End Do
    If (TwoLoopRGE) Then
     g1_h24 = g1_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                &                     * 5._dp/2._dp*Log(MassMXM3(3)/MWM30(3,3)) )
     g2_h24 = g2_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                &                      *( 1.5_dp *Log(MassMXM3(3)/MWM30(3,3)) + &
                &                      2._dp *Log(MassMWM3(3)/MWM30(3,3))) )
     g3_h24 = g3_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                       * (Log(MassMXM3(3)/MWM30(3,3)) &
                &                         + 3._dp*Log(MassMGM3(3)/MWM30(3,3)) ) )

     MassB_h24 = MassB_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                &                     * 5._dp/2._dp*Log(MassMXM3(3)/MWM30(3,3)) )
     MassWB_h24 = MassWB_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                &                      *( 1.5_dp *Log(MassMXM3(3)/MWM30(3,3)) + &
                &                     2._dp *Log(MassMWM3(3)/MWM30(3,3))) )
     MassG_h24 = MassG_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                       * (Log(MassMXM3(3)/MWM30(3,3)) &
                &                         + 3._dp*Log(MassMGM3(3)/MWM30(3,3)) ) )
    End If

   Else 

    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = 0.5_dp * Yw3_H24(3,i1)*Yw3_H24(3,i2)/MWM30(3,3) + &
                   & 0.3_dp * Yb3_H24(3,i1)*Yb3_H24(3,i2)/MWM30(3,3) 
     End Do
    End Do

   End If

   Yb30_H24(3,:,:) = Yb3_H24
   Yw30_H24(3,:,:) = Yw3_H24
   Yx30_H24(3,:,:) = Yx3_H24
   Yb3_H24(3,:) = 0._dp
   Yw3_H24(3,:) = 0._dp
   Yx3_H24(3,:) = 0._dp
   AYb3_H24(3,:) = 0._dp
   AYw3_H24(3,:) = 0._dp
   AYx3_H24(3,:) = 0._dp
   MWM3(:,3) = 0._dp
   MWM3(3,:) = 0._dp
   MGM3(:,3) = 0._dp
   MGM3(3,:) = 0._dp
   MXM3(:,3) = 0._dp
   MXM3(3,:) = 0._dp
   MBM3(:,3) = 0._dp
   MBM3(3,:) = 0._dp

   mHx32(3,:) = 0._dp 
   mHx32(:,3) = 0._dp 
   mHxb32(3,:) = 0._dp 
   mHxb32(:,3) = 0._dp 
   mHg32(3,:) = 0._dp 
   mHg32(:,3) = 0._dp 
   mHb32(3,:) = 0._dp 
   mHb32(:,3) = 0._dp 
   mHw32(3,:) = 0._dp 
   mHw32(:,3) = 0._dp 

   Call ParametersToG555(g1_H24,g2_H24,g3_H24,Yu_H24,Yd_H24,Ye_H24,Yb3_H24   &
      & ,Yw3_H24,Yx3_H24,mu_H24,MXM3,MWM3,MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24  &
      & ,AYb3_H24,AYw3_H24,AYx3_H24,Amue,AMXM3,AMWM3,AMGM3,AMBM3, mq2_H24    &
      & ,ml2_H24,mHd2_H24,mHu2_H24,md2_H24,mu2_H24,me2_H24,mHw32,mHg32,mHb32 &
      & ,mHx32,mHxb32,MassB_H24,MassWB_H24,MassG_H24,MnuL5,g2f)

   !------------
   ! 3 -> 2
   !------------
   If (Abs(MWM30(2,2)).Lt.Abs(m_hi)) Then
    m_lo = Abs(MWM30(2,2))
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2f, 573, tz, 0._dp, delta, dt, 0._dp, rge555, kont)
    m_hi = m_lo
   Endif

   Call GToParameters555(g2f,g1_H24,g2_H24,g3_H24,Yu_H24,Yd_H24,Ye_H24,Yb3_H24 &
      & ,Yw3_H24,Yx3_H24,mu_H24,MXM3,MWM3,MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24    &
      & ,AYb3_H24,AYw3_H24,AYx3_H24,Amue,AMXM3,AMWM3,AMGM3,AMBM3, mq2_H24      &
      & ,ml2_H24,mHd2_H24,mHu2_H24,md2_H24,mu2_H24,me2_H24,mHw32,mHg32,mHb32   &
      & ,mHx32,mHxb32,MassB_H24,MassWB_H24,MassG_H24,MnuL5)

   MBM3Running(2,:,:)=MBM3
   MGM3Running(2,:,:)=MGM3
   MWM3Running(2,:,:)=MWM3
   MXM3Running(2,:,:)=MXM3

   NGHb3 = 1._dp
   NGHg3 = 1._dp 
   NGHw3 = 1._dp 
   NGHx3 = 1._dp 
   NGHxb3 = 1._dp 
   ThresholdCrossed = 1

   RotWR = 0._dp
   RotWL = 0._dp
   RotXR = 0._dp
   RotXL = 0._dp
   RotGR = 0._dp
   RotGL = 0._dp
   RotBR = 0._dp
   RotBL = 0._dp
   EigMWM3 = 0._dp
   EigMGM3 = 0._dp
   EigMXM3 = 0._dp
   EigMBM3 = 0._dp

   Call EigenSystem(Matmul(Transpose(Conjg(MWM3(1:2,1:2))),MWM3(1:2,1:2)) &
    & , EigMWM3(1:2),RotWR(1:2,1:2), ierr, test2)
   Call EigenSystem(Matmul(MWM3(1:2,1:2),Transpose(Conjg(MWM3(1:2,1:2)))) &
    & , EigMWM3(1:2),RotWL(1:2,1:2), ierr, test2)

   Call EigenSystem(Matmul(Transpose(Conjg(MBM3(1:2,1:2))),MBM3(1:2,1:2)) &
    & , EigMBM3(1:2),RotBR(1:2,1:2), ierr, test2)
   Call EigenSystem(Matmul(MBM3(1:2,1:2),Transpose(Conjg(MBM3(1:2,1:2)))) &
    & , EigMBM3(1:2),RotBL(1:2,1:2), ierr, test2)

   Call EigenSystem(Matmul(Transpose(Conjg(MGM3(1:2,1:2))),MGM3(1:2,1:2)) &
    & , EigMGM3(1:2),RotGR(1:2,1:2), ierr, test2)
   Call EigenSystem(Matmul(MGM3(1:2,1:2),Transpose(Conjg(MGM3(1:2,1:2)))) &
    & , EigMGM3(1:2),RotGL(1:2,1:2), ierr, test2)

   Call EigenSystem(Matmul(Transpose(Conjg(MXM3(1:2,1:2))),MXM3(1:2,1:2)) &
    & , EigMXM3(1:2),RotXR(1:2,1:2), ierr, test2)
   Call EigenSystem(Matmul(MXM3(1:2,1:2),Transpose(Conjg(MXM3(1:2,1:2)))) &
    & , EigMXM3(1:2),RotXL(1:2,1:2), ierr, test2)


   Yb3_H24 = Matmul(RotBL,Yb3_H24)
   Yw3_H24 = Matmul(RotWL,Yw3_H24)
   Yx3_H24 = Matmul(RotXL,Yx3_H24)
   AYb3_H24 = Matmul(RotBL,AYb3_H24)
   AYw3_H24 = Matmul(RotwL,AYb3_H24)
   AYx3_H24 = Matmul(RotxL,AYb3_H24)
   MWM3 = Matmul(Conjg(Transpose(RotWR)),Matmul(MWM3,RotWL))
   MGM3 = Matmul(Conjg(Transpose(RotGR)),Matmul(MGM3,RotGL))
   MXM3 = Matmul(Conjg(Transpose(RotXR)),Matmul(MXM3,RotXL))
   MBM3 = Matmul(Conjg(Transpose(RotBR)),Matmul(MBM3,RotBL))
   mHx32 = Matmul(Conjg(Transpose(RotXL)),Matmul(mHx32,RotXL)) 
   mHxb32 = Matmul(Conjg(Transpose(RotXR)),Matmul(mHxb32,RotXR))  
   mHg32 = Matmul(Conjg(Transpose(RotGL)),Matmul(mHg32,RotGL)) 
   mHb32 = Matmul(Conjg(Transpose(RotBL)),Matmul(mHb32,RotBL))  
   mHw32 = Matmul(Conjg(Transpose(RotWL)),Matmul(mHw32,RotWL))
  
   MassMWM3(2) = Sqrt(EigMWM3(2))
   MassMXM3(2) = Sqrt(EigMXM3(2))
   MassMGM3(2) = Sqrt(EigMGM3(2))
   MassMBM3(2) = Sqrt(EigMBM3(2))

   If (MassMWM3(2).Lt.Abs(mGUT)) Then

    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = MnuL5(i1,i2) &
                & + 0.5_dp * Yw3_H24(2,i1)*Yw3_H24(2,i2)/MassMWM3(2)  &
                & + 0.3_dp * Yb3_H24(2,i1)*Yb3_H24(2,i2)/MassMBM3(2) 
     End Do
    End Do

    If (TwoLoopRGE) Then
     g1_h24 = g1_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                &                   * 5._dp/2._dp*Log( MassMXM3(2)/MWM30(2,2)) )
     g2_h24 = g2_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                &                     *( 1.5_dp *Log(MassMXM3(2)/MWM30(2,2)) + &
                &                      2._dp *Log(MassMWM3(2)/MWM30(2,2))) )
     g3_h24 = g3_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                     * (Log( MassMXM3(2)/MWM30(2,2)) &
                &                       + 3._dp*Log( MassMGM3(2)/MWM30(2,2)) ) )

     MassB_h24 = MassB_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                 &                   * 5._dp/2._dp*Log( MassMXM3(2)/MWM30(2,2)) )
     MassWB_h24 = MassWB_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                 &                     *( 1.5_dp *Log(MassMXM3(2)/MWM30(2,2)) + &
                 &                       2._dp *Log(MassMWM3(2)/MWM30(2,2))) )
     MassG_h24 = MassG_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                      * (Log( MassMXM3(2)/MWM30(2,2)) &
                &                        + 3._dp*Log( MassMGM3(2)/MWM30(2,2)) ) )
    End If

   Else 

    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = MnuL5(i1,i2)   &
                 & + 0.5_dp * Yw3_H24(2,i1)*Yw3_H24(2,i2)/mGut  &
                 & + 0.3_dp * Yb3_H24(2,i1)*Yb3_H24(2,i2)/mGut 
     End Do
    End Do

   End If

   Yb30_H24(2,:,:) = Yb3_H24
   Yw30_H24(2,:,:) = Yw3_H24
   Yx30_H24(2,:,:) = Yx3_H24
   Yb3_H24(2,:) = 0._dp
   Yw3_H24(2,:) = 0._dp
   Yx3_H24(2,:) = 0._dp
   AYb3_H24(2,:) = 0._dp
   AYw3_H24(2,:) = 0._dp
   AYx3_H24(2,:) = 0._dp
   MWM3(:,2) = 0._dp
   MWM3(2,:) = 0._dp
   MGM3(:,2) = 0._dp
   MGM3(2,:) = 0._dp
   MXM3(:,2) = 0._dp
   MXM3(2,:) = 0._dp
   MBM3(:,2) = 0._dp
   MBM3(2,:) = 0._dp

   mHx32(2,:) = 0._dp 
   mHx32(:,2) = 0._dp 
   mHxb32(2,:) = 0._dp 
   mHxb32(:,2) = 0._dp 
   mHg32(2,:) = 0._dp 
   mHg32(:,2) = 0._dp 
   mHb32(2,:) = 0._dp 
   mHb32(:,2) = 0._dp 
   mHw32(2,:) = 0._dp 
   mHw32(:,2) = 0._dp 

   Call ParametersToG555(g1_H24,g2_H24,g3_H24,Yu_H24,Yd_H24,Ye_H24,Yb3_H24 &
      & ,Yw3_H24,Yx3_H24,mu_H24,MXM3,MWM3,MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24    &
      & ,AYb3_H24,AYw3_H24,AYx3_H24,Amue,AMXM3,AMWM3,AMGM3,AMBM3, mq2_H24      &
      & ,ml2_H24,mHd2_H24,mHu2_H24,md2_H24,mu2_H24,me2_H24,mHw32,mHg32,mHb32   &
      & ,mHx32,mHxb32,MassB_H24,MassWB_H24,MassG_H24,MnuL5,g2f)

   !------------
   ! 2 -> 1
   !------------

   If (Abs(MWM30(1,1)).Lt.Abs(m_hi)) Then
    m_lo = MWM30(1,1)
    tz = Log(Abs(m_hi/m_lo))
    dt = - tz / 50._dp
    Call odeint(g2f, 573, tz, 0._dp, delta, dt, 0._dp, rge555, kont)
    m_hi = m_lo
   Endif

   Call GToParameters555(g2f,g1_H24,g2_H24,g3_H24,Yu_H24,Yd_H24,Ye_H24,Yb3_H24 &
      & ,Yw3_H24,Yx3_H24,mu_H24,MXM3,MWM3,MGM3,MBM3,AYu_H24,AYd_H24,AYe_H24    &
      & ,AYb3_H24,AYw3_H24,AYx3_H24,Amue,AMXM3,AMWM3,AMGM3,AMBM3, mq2_H24      &
      & ,ml2_H24,mHd2_H24,mHu2_H24,md2_H24,mu2_H24,me2_H24,mHw32,mHg32,mHb32   &
      & ,mHx32,mHxb32,MassB_H24,MassWB_H24,MassG_H24,MnuL5)

   MBM3Running(1,:,:)=MBM3
   MGM3Running(1,:,:)=MGM3
   MWM3Running(1,:,:)=MWM3
   MXM3Running(1,:,:)=MXM3

   NGHb3 = 0._dp
   NGHg3 = 0._dp 
   NGHw3 = 0._dp 
   NGHx3 = 0._dp 
   NGHxb3 = 0._dp 
   ThresholdCrossed = 0

   MassMWM3(1) = Abs(MWM3(1,1))
   MassMXM3(1) = Abs(MXM3(1,1))
   MassMGM3(1) = Abs(MGM3(1,1))
   MassMBM3(1) = Abs(MBM3(1,1))

   Yb30_H24(1,:,:) = Yb3_H24
   Yw30_H24(1,:,:) = Yw3_H24
   Yx30_H24(1,:,:) = Yx3_H24
 
   gauge_h24(1) = g1_h24
   gauge_h24(2) = g2_h24
   gauge_h24(3) = g3_h24
   Mi_h24(1) = MassB_h24
   Mi_h24(2) = MassWB_h24
   Mi_h24(3) = MassG_h24
   mh2_h24(1) = mHd2_h24
   mh2_h24(2) = mHu2_h24

   If (MassMWM3(1).Lt.Abs(mGUT)) Then

    If (TwoLoopRGE) Then
     g1_h24 = g1_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                &                   * 5._dp/12._dp*Log( MassMXM3(1)/MWM30(1,1)) )
     g2_h24 = g2_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                &                      *( 0.5_dp *Log(MassMXM3(1)/MWM30(1,1)) + &
                &                       2._dp *Log(MassMWM3(1)/MWM30(1,1))) )
     g3_h24 = g3_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                      * (0.5_dp*Log( MassMXM3(1)/MWM30(1,1)) &
                &                        + 3._dp*Log( MassMGM3(1)/MWM30(1,1)) ) )

     MassB_h24 = MassB_h24 * (1._dp + oo16pi2 * g1_h24**2           &
                &                   * 5._dp/12._dp*Log( MassMXM3(1)/MWM30(1,1)) )
     MassWB_h24 = MassWB_h24 * (1._dp + oo16pi2 * g2_h24**2           &
                 &                    *( 0.5_dp *Log(MassMXM3(1)/MWM30(1,1)) + &
                 &                       2._dp *Log(MassMWM3(1)/MWM30(1,1))) )
     MassG_h24 = MassG_h24 * (1._dp + oo16pi2 * g3_h24**2           &
                &                      * (0.5_dp*Log( MassMXM3(1)/MWM30(1,1)) &
                &                        + 3._dp*Log( MassMGM3(1)/MWM30(1,1)) ) )
    End If
 
    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = MnuL5(i1,i2) &
                & + 0.5_dp * Yw3_H24(1,i1)*Yw3_H24(1,i2)/ MassMWM3(1)  &
                & + 0.3_dp * Yb3_H24(1,i1)*Yb3_H24(1,i2)/ MassMBM3(1) 
     End Do
    End Do

   Else 

    Do i1=1,3
     Do i2=1,3
      MnuL5(i1,i2) = MnuL5(i1,i2)   &
                & + 0.5_dp * Yw3_H24(1,i1)*Yw3_H24(1,i2)/Abs(MWM30(1,1))  &
                & + 0.3_dp * Yb3_H24(1,i1)*Yb3_H24(1,i2)/Abs(MWM30(1,1)) 
     End Do
    End Do

   End If


  Call ParametersToG4(gauge_H24, ye_H24, Zero33C, yd_H24, yu_H24   &
          & , ZeroC, ZeroC, Mi_H24, AYe_h24, Zero33C, AYd_h24 &
          & , AYu_h24, ZeroC, ZeroC, Me2_h24, ML2_h24     &
          & , MD2_h24, MQ2_h24, MU2_h24, MH2_h24, ZeroR2, mu_h24      &
          & , Amue, MnuL5, g2c)

  mudim = GetRenormalizationScale()
  mudim = Max(mudim, mZ2)

  tz = 0.5_dp * Log(mudim/m_hi**2)
  dt = tz / 100._dp
  Call odeint(g2c, 277, 0._dp, tz, delta, dt, 0._dp, rge277, kont)

  Call GToParameters4(g2c, gauge, y_l, y_T, y_d, y_u, lam12(1), lam12(2), Mi &
          & , A_l, A_T, A_d, A_u, Alam12(1), Alam12(2), M2_E, M2_L, M2_D      &
          & , M2_Q, M2_U, M2_H, M2_T, mu, B, MnuL5)

  Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u, M2_E, M2_L  &
      & , M2_D, M2_Q, M2_U, M2_H, mu, B, g2)


 Else If (HighScaleModel.Eq.'SEESAW_II_SARAH') Then

   If ( (oo4pi*Maxval(g2g(1:115)**2)).Gt.0.5_dp) Then
    Write(ErrCan,*) "Non perturbative regime at M_GUT"
    If (ErrorLevel.Ge.2) Call TerminateProgram
    Write(errcan,*) " "
    kont = -418
    Call AddError(418)
    Iname = Iname - 1
    Return
   End If

   !-------------------------
   ! run only if m_H < m_GUT
   !-------------------------

   If (Abs(MTM_GUT).Lt.mGUT) Then
    !---------------------------------------------------
    ! heavy states are to be included
    !---------------------------------------------------
    

    mudim = Abs(MTM_GUT)**2
    tz = 0.5_dp * Log(mudim/mGUT**2)
    dt = tz / 50._dp
    Call odeint(g2g, 365 , 0._dp, tz, delta, dt, 0._dp, rge365, kont)
    m_lo = Abs(MTM_GUT)

    If ( (oo4pi*Maxval(g2g(1:115)**2)).Gt.0.5_dp) Then
     Write(ErrCan,*) "Non perturbative regime at M_H3"
     If (ErrorLevel.Ge.2) Call TerminateProgram
     Write(errcan,*) " "
     kont = -419
     Call AddError(419)
     Iname = Iname - 1
     Return
    End If

   Else
     m_lo = mGUT
   End If


   Call GToParameters365(g2g,g1_h15,g2_h15,g3_h15,Yu_h15,Yd_h15,Ye_h15,Yt0_h15 &
      & ,Ys0_h15,Yz0_h15, Lambda10,Lambda20,mu_h15,MTM,MZM,MSM,AYu_h15,AYd_h15 &
      & ,AYe_h15,AYt_h15,AYs_h15,AYz_h15,ALambda1,ALambda2,Amue,AMTM,AMZM,AMSM &
      & ,mq2_h15,ml2_h15,mHd2_h15, mHu2_h15,md2_h15,mu2_h15,me2_h15,mt2_H15    &
      & ,mtb2_H15,ms2_H15,msb2_H15,mz2_H15,mzb2_H15,MassB_h15,MassWB_h15       &
      & ,MassG_h15,MnuL5)

   gauge_h15(1) = g1_h15
   gauge_h15(2) = g2_h15
   gauge_h15(3) = g3_h15
   Mi_h15(1) = MassB_h15
   Mi_h15(2) = MassWB_h15
   Mi_h15(3) = MassG_h15
   mh2_h15(1) = mHd2_h15
   mh2_h15(2) = mHu2_h15

   MnuL5 = - Lambda20 * YT0_h15 / MTM


   If ((Abs(MTM_GUT).Lt.mGUT).And.TwoLoopRGE) Then
   !-----------------------------------------------------
   ! adding shifts to gauge couplings and gaugino masses
   !-----------------------------------------------------
    gauge_h15(1) = gauge_h15(1) / (1._dp - 0.2_dp * oo16pi2 * gauge_h15(1)**2 &
                &                       * ( 8._dp * Log(Abs(MSM/MTM_GUT))     &
                &                         + 9._dp * Log(Abs(MTM/MTM_GUT))     &
                &                         + 0.5_dp * Log(Abs(MZM/MTM_GUT)) ) )
    gauge_h15(2) = gauge_h15(2) / (1._dp - oo16pi2 * gauge_h15(2)**2          &
              & *(2._dp*Log(Abs(MTM/MTM_GUT))+ 1.5_dp *Log(Abs(MZM/MTM_GUT)) ))
    gauge_h15(3) = gauge_h15(3) / (1._dp - oo16pi2 * gauge_h15(3)**2          &
                &   * (2.5_dp*Log(Abs(MSM/MTM_GUT)) + Log(Abs(MZM/MTM_GUT)) ) )
    Mi_h15(1) = Mi_h15(1) / (1._dp - 0.2_dp * oo16pi2 * gauge_h15(1)**2       &
                &                       * ( 8._dp * Log(Abs(MSM/MTM_GUT))     &
                &                         + 9._dp * Log(Abs(MTM/MTM_GUT))     &
                &                         + 0.5_dp * Log(Abs(MZM/MTM_GUT)) ) )
    Mi_h15(2) = Mi_h15(2) / (1._dp - oo16pi2 * gauge_h15(2)**2                &
             & *(2._dp *Log(Abs(MTM/MTM_GUT))+ 1.5_dp *Log(Abs(MZM/MTM_GUT)) ))
    Mi_h15(3) = Mi_h15(3) / (1._dp - oo16pi2 * gauge_h15(3)**2           &
                &   * (2.5_dp*Log(Abs(MSM/MTM_GUT)) + Log(Abs(MZM/MTM_GUT)) ) )
   End If

   !---------------------------------------------------
   ! heavy states are not to be included
   !---------------------------------------------------

   Call ParametersToG4(gauge_h15, ye_h15, Zero33C, yd_h15, yu_h15   &
          & , ZeroC, ZeroC, Mi_h15, AYe_h15, Zero33C, AYd_h15 &
          & , AYu_h15, ZeroC, ZeroC, ME2_h15, ML2_h15     &
          & , MD2_h15, MQ2_h15, MU2_h15, MH2_h15, ZeroR2, mu_h15      &
          & , Amue, MnuL5, g2c)

   mudim = GetRenormalizationScale()
   mudim = Max(mudim, mZ2)

   tz = 0.5_dp * Log(mudim/m_lo**2)
   dt = tz / 100._dp
   Call odeint(g2c, 277, 0._dp, tz, delta, dt, 0._dp, rge277, kont)



   Call GToParameters4(g2c, gauge, y_l, y_T, y_d, y_u, lam12(1), lam12(2), Mi &
          & , A_l, A_T, A_d, A_u, Alam12(1), Alam12(2), M2_E, M2_L, M2_D      &
          & , M2_Q, M2_U, M2_H, M2_T, mu, B, MnuL5)
   Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u, M2_E, M2_L  &
      & , M2_D, M2_Q, M2_U, M2_H, mu, B, g2)

# endif SEESAWIII
! Florian Staub Seesaw II+III

  Else

   mudim = GetRenormalizationScale()
   mudim = Max(mudim, mZ2)

   tz = 0.5_dp * Log(mudim/mGUT_save**2)
   dt = tz / 100._dp

   Call odeint(g2, 213, 0._dp, tz, delta, dt, 0._dp, rge213, kont)

  End If

  Iname = Iname - 1


  900 Format(a20,e15.6)
  910 Format(a15,3e15.6)

 End Subroutine RunRGE_2_boundaries

 Subroutine RunRGEup(g2, mGUT, Qvec, g_out, g_out2, kont)
 !-----------------------------------------------------------------------
 ! Uses Runge-Kutta method to integrate RGE's from M_Z to M_GUT
 ! Written by Werner Porod, 17.04.02
 ! 17.04.02 : - including the mSugra case and a simplified version
 !              for varying m_t
 !            - setting flag for writing
 ! 01.02.03: remodelling everything, because this routine should only
 !           serve as possiblity to run up RGEs, structure is build such
 !           that right handed neutrinos can be included
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: g2(:), mGUT
  Real(dp), Intent(out) :: g_out(:,:), g_out2(:,:), Qvec(:)
  
  Logical :: check
  Integer:: i1, len1, len2, len3, len4
  Real(dp) :: g2a(213), w2(213,3), g2b(267), w2b(267,3), tz, dt, t, mudim

  Iname = Iname + 1
  NameOfUnit(Iname) = 'RunRGEup'

  !-------------------------------------------
  ! first a check if everthing is consistent
  !-------------------------------------------
  len1 = Size( g2 )
  len2 = Size( Qvec )
  len3 = Size( g_out, dim=1 )
  len4 = Size( g_out, dim=2 )

  check = .True.  
  If (len1.Ne.len4) check = .False. ! size of parameter vector do not conincide
  If (len2.Ne.len3) check = .False. ! number of steps do not conincide
  If (.Not.check) Then
   Write(ErrCan,*) "Problem in routine "//NameOfUnit(Iname)
   Write(ErrCan,*) "dimension of g2   :",len1
   Write(ErrCan,*) "dimension of Qvec :",len2
   Write(ErrCan,*) "dimension of g2   :",len1
   If (ErrorLevel.Ge.0) Call TerminateProgram
  End If   

  !---------------------------------------------
  ! data for running: mudim is low energy scale
  !---------------------------------------------
  mudim = GetRenormalizationScale()
  mudim = Sqrt(mudim)
  tz = Log(mudim/mGUT)
  dt = - tz / Real(len2-1, dp)
  g2a = g2
  !------------------
  ! now the running
  !------------------
  g_out(1,:) = g2a
  g_out2 = 0._dp
  If (HighScaleModel.Eq.'SUGRA_NuR') Then
   tz = Log(mudim/mNuR(3))
   dt = - tz / Real(len2-18, dp)
   Do i1=0,len2-19
    t =  tz + dt * i1
    Qvec(i1+1) = mNuR(3) * Exp(t)
    Call rkstp(dt,t,g2a,rge213,w2)
    g_out(2+i1,:) = g2a
   End Do
    Qvec(i1+1) = mNuR(3)
   Call GToParameters(g2a, gauge_mR, y_l_mR, y_d_mR, y_u_mR, Mi_mR &
          & , A_l_mR, A_d_mR, A_u_mR, M2_E_mR, M2_L_mR     &
          & , M2_D_mR, M2_Q_mR, M2_U_mR, M2_H_mR, mu_mR, B_mR)

   Call ParameterstoG2( gauge_mR, y_l_mR, y_nu_mR, y_d_mR, y_u_mR, Mi_mR &
          & , A_l_mR, A_nu_mR, A_d_mR, A_u_mR, M2_E_mR, M2_L_mR, M2_R_mR     &
          & , M2_D_mR, M2_Q_mR, M2_U_mR, M2_H_mR, mu_mR, B_mR, g2b)


   g_out2(1+i1,:) = g2b
   tz = Log(MnuR(3)/mGUT)
   dt = - tz / Real(17, dp)
   Do i1=len2-18,len2-2
    t =  tz + dt * (i1 - len2 + 18)
    Qvec(i1+1) = Mgut * Exp(t)
    Call rkstp(dt,t,g2a,rge213,w2)
    g_out(2+i1,:) = g2a
    t =  tz + dt * (i1 - len2 + 18)
    Call rkstp(dt,t,g2b,rge267,w2b)
    g_out2(2+i1,:) = g2b

   End Do

  Else

   Do i1=0,len2-2
    t =  tz + dt * i1
    Qvec(i1+1) = Mgut * Exp(t)
    Call rkstp(dt,t,g2a,rge213,w2)
    g_out(2+i1,:) = g2a
   End Do
  End If

  Qvec(len2) = Mgut

 Iname = Iname - 1
 
 End Subroutine RunRGEup


 Subroutine RunRGEup2(lenR, lenQ, g2, Qi, Qf, Qvec, g_out)
 !-----------------------------------------------------------------------
 ! Uses Runge-Kutta method to integrate RGE's from Qi to Qf
 ! input:
 !   lenR .. length of data vector
 !   lenQ .. length of Q vector 
 !   g2 .... initial data
 !   Qi .... starting scale
 !   Qf .... final scale
 ! output:
 !   Qvec .. vector containing the scales
 !   g_out.. vector containing the data at the various scale
 ! Written by Werner Porod, 17.04.02
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: lenR, lenQ
  Real(dp), Intent(in) :: g2(lenR), Qi, Qf
  Real(dp), Intent(out) :: g_out(lenR,lenQ), Qvec(lenQ)
  
  Integer:: i1
  Real(dp) :: g2a(lenR), w2(lenR,3), tz, dt, t

  Iname = Iname + 1
  NameOfUnit(Iname) = 'RunRGEup2'

  tz = Log(Qi/Qf)
  dt = - tz / Real(lenQ-1, dp)
  g2a = g2
  !------------------
  ! now the running
  !------------------
  g_out(:,1) = g2a
  Qvec(1) = Qi
  Do i1=0,lenQ-2 
   t =  tz + dt * i1
   Qvec(i1+1) = Qf * Exp(t)
   If (lenR.Eq.180) Call rkstp(dt,t,g2a,rge_SU5,w2)
   If (lenR.Eq.213) Call rkstp(dt,t,g2a,rge213,w2)
   If (lenR.Eq.267) Call rkstp(dt,t,g2a,rge267,w2)
   If (lenR.Eq.285) Call rkstp(dt,t,g2a,rge285,w2)
   g_out(:,2+i1) = g2a
  End Do

  Qvec(lenQ) = Qf

 Iname = Iname - 1
 
 End Subroutine RunRGEup2



 Logical Function SetCheckSugraDetails(V1, V2, V3, V4, V5)
 !----------------------------------------------------------------------------
 ! Sets the variable CheckSugraDetails which controls writing of the details
 ! during the run. Default is .False. In case that one of the entries .True.
 ! the following action is performed:
 ! CheckSugraDetails(1) -> write high scale values for gauge and Yukawa
 !                         couplings to channel 10 for each iteration
 ! CheckSugraDetails(2) -> write low scale values for gauge and Yukawa
 !                         couplings to channel 10 for each iteration
 ! CheckSugraDetails(3) -> write high scale values for soft SUSY parameters
 !                         to channel 10 for each iteration
 ! CheckSugraDetails(4) -> write low scale values for soft SUSY parameters
 !                         to channel 10 for each iteration
 ! CheckSugraDetails(5) -> write inital low scale values for gauge and Yukawa
 !                         couplings to channel 10 for each iteration
 ! Written by Werner Porod, 24..09.01
 !----------------------------------------------------------------------------
 Implicit None
  Logical, Intent(in) :: V1, V2, V3, V4, V5
   SetCheckSugraDetails = .False.
   CheckSugraDetails(1) = V1
   CheckSugraDetails(2) = V2
   CheckSugraDetails(3) = V3
   CheckSugraDetails(4) = V4
   CheckSugraDetails(5) = V5
   SetCheckSugraDetails = .True.
 End Function SetCheckSugraDetails


 Subroutine SetGUTScale(scale)
 Implicit None
  Real(dp), Intent(in) :: scale

  If (scale.Lt.0._dp) Then
   UseFixedGUTScale = .False.
  Else
   UseFixedGUTScale = .True.
   GUT_scale = scale
  End If

 End Subroutine SetGUTScale


 Logical Function SetHighScaleModel(model)
  Implicit None
  Character(len=*), Intent(in) :: model

  SetHighScaleModel = .False.
  HighScaleModel = model
  SetHighScaleModel = .True.
 End  Function SetHighScaleModel

 Subroutine SetRGEScale(scale)
 Implicit None
  Real(dp), Intent(in) :: scale

  Real(dp) :: old_scale

  If (scale.Lt.0._dp) Then
   UseFixedScale = .False.
  Else
   UseFixedScale = .True.
   old_scale = SetRenormalizationScale(scale)
  End If

 End Subroutine SetRGEScale

 Logical Function SetStrictUnification(V1)
 !-----------------------------------------------------------------------
 ! Sets the parameter StrictUnification, which enforces g_3 = g_1 = g_2
 ! at the high scale, default is .false.
 ! written by Werner Porod, 24.09.01
 !-----------------------------------------------------------------------
 Implicit None
  Logical, Intent(in) :: V1
  SetStrictUnification = .False.
  StrictUnification = V1
  SetStrictUnification = .True.
 End Function SetStrictUnification


 Integer Function SetYukawaScheme(V1)
 !-----------------------------------------------------------------------
 ! Sets the parameter YukawaScheme, which controls wheter the top (=1) or the
 ! down (=2) Yukawa couplings stay diagonal at the low scale 
 ! written by Werner Porod, 20.11.01
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: V1
  SetYukawaScheme = YukawaScheme
  YukawaScheme = V1
 End Function SetYukawaScheme

 Subroutine Sugra(delta, vevSM, mC, U, V, mN, N, mS0, mS02, RS0, mP0, mP02,RP0&
    & , mSpm, mSpm2, RSpm, mDsquark, mDsquark2, RDsquark, mUsquark, mUsquark2 &
    & , RUsquark, mSlepton, mSlepton2, RSlepton, mSneutrino, mSneutrino2      &
    & , RSneutrino, mGlu, phase_glu, gauge, uL_L, uL_R, uD_L, uD_R, uU_L      &
    & , uU_R, Y_l, Y_d  &
    & , Y_u, Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B     &
    & , mGUT, kont, WriteComment, niter)
 !-----------------------------------------------------------------------
 ! Computes RGE's of the the SUSY parameters 
 ! Uses Runge-Kutta method to integrate RGE's from M_Z to M_GUT
 ! and back, putting in correct thresholds. For the first iteration
 ! only the first 6 couplings are included and a common threshold
 ! is used.
 ! written by Werner Porod, 4.8.1999
 ! 24.09.01: portation to f90
 !           taking masses as input to see whether situation improves
 !           in principal they can be zero, but this in general requires
 !           one additional run. therefore a good educated guess is
 !           useful. The parameters might be generated with the help
 !           of the routine FirstGuess which are then plugged into
 !           TreeMasses
 !         - Contrary to the old version I put now the calculation of the
 !           tree level masses to the routine LoopMasses
 ! 16.11.01: including logical variable WriteComment to get a better control
 !           on the various steps
 ! 16.09.02: the electroweak scale is now put here, except it has already
 !           been fixed in the main program
 ! 05.01.04: calculating now gauge and Yukawa couplings at m_Z using
 !           running masses, more Yukawa couplings from the previous run
 !           are taken as starting point for the next run
 !-----------------------------------------------------------------------
 Implicit None

  Logical, Intent(in) ::  WriteComment
  Integer, Intent(in) :: niter
  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: delta
  Real(dp), Intent(inout) :: vevSM(2)
  Real(dp), Intent(out) :: mGUT
  Real(dp), Intent(inout) :: mC(:), mN(:), mS0(:), mP0(:), mSpm(:) &
    & , mUsquark(:), mDsquark(:), mSlepton(:), mSneutrino(:)       &
    & , mUsquark2(:), mDsquark2(:), mSlepton2(:), mSneutrino2(:)   &
    & , mS02(:), mP02(:), mSpm2(:), RP0(:,:), RS0(:,:), mglu
  Complex(dp), Intent(inout) :: U(:,:), V(:,:), N(:,:), RSpm(:,:)      &
    & , RDsquark(:,:), RUsquark(:,:), RSlepton(:,:), RSneutrino(:,:)   &
    & , phase_Glu
  Real(dp), Intent(inout) :: gauge(3), M2_H(2)
  Complex(dp), Dimension(3,3), Intent(inout) :: Y_l, Y_d, Y_u, A_l, A_d, A_u &
          & , M2_E, M2_L, M2_D, M2_Q, M2_U, uL_L, uL_R, uU_L, uU_R, uD_L, uD_R
  Complex(dp), Intent(inout) :: Mi(3), mu, B

  Real(dp) :: deltag0, tanb, g0(213), t1, t2, mZ2_run, mW2_run, g1(57) &
    & , mc2(2), mn2(4), mudim, tz, dt, mc2_T(2), mn2_T(4), vev
  Integer :: j, n_C, n_N, n_S0, n_P0, n_Spm, n_Su, n_Sd, n_Sl &
    & , n_Sn, i1, n_tot
  Real(dp) :: mC_T(2), mN_T(4), mS0_T(2), mP0_T(2), mSpm_T(2)              &
    & , mUsquark_T(6), mDsquark_T(6), mSlepton_T(6), mSneutrino_T(3)       &
    & , mUsquark2_T(6), mDsquark2_T(6), mSlepton2_T(6), mSneutrino2_T(3)   &
    & , mS02_T(2), mP02_T(2), mSpm2_T(2), RP0_T(2,2), RS0_T(2,2), mglu_T   &
    & , mf_nu(3), Abs_Mu2, Abs_Mu, g58(58), g2(213), tanbQ
  Complex(dp) :: U_T(2,2), V_T(2,2), N_T(4,4), RSpm_T(2,2)                   &
    & , RDsquark_T(6,6), RUsquark_T(6,6), RSlepton_T(6,6), RSneutrino_T(3,3) &
    & , phase_Glu_T, Unu(3,3), mu_T, B_T
  Complex(dp), Dimension(3,3) :: Ad_s, Au_s, M2D_s, M2Q_s, M2U_s, Al_s, M2E_s &
    & , M2L_s

  Logical :: FoundResult
  Real(dp), Allocatable :: mass_new(:), mass_old(:), diff_m(:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Sugra'
  !-----------------
  ! Inititalization
  !-----------------
  kont = 0 
  FoundResult = .False.
  tanb = vevSM(2) / vevSM(1)
  !--------------------------------------------------------------------
  ! saving masses for checking + creation of corresponding variables
  !-------------------------------------------------------------------
  n_C = Size( mC )
  n_N = Size( mN )
  n_S0 = Size( mS0 )
  n_P0 = Size( mP0 )
  n_Spm = Size( mSpm )
  n_Su = Size( mUsquark )
  n_Sd = Size( mDsquark )
  n_Sl = Size( mSlepton )
  n_Sn = Size( mSneutrino )
  n_tot = 1 + n_C + n_N + n_S0 + n_P0 + n_Spm + n_Su + n_Sd + n_Sl + n_Sn
  Allocate( mass_old( n_tot ), mass_new( n_tot ), diff_m( n_tot) )
  mass_old(1) = mGlu
  n_tot = 1
  mass_old(n_tot+1:n_tot+n_C) = mC
  n_tot = n_tot + n_C
  mass_old(n_tot+1:n_tot+n_N) = mN
  n_tot = n_tot + n_N
  mass_old(n_tot+1:n_tot+n_S0) = mS0
  n_tot = n_tot + n_S0
  mass_old(n_tot+1:n_tot+n_P0) = mP0
  n_tot = n_tot + n_P0
  mass_old(n_tot+1:n_tot+n_Spm) = mSpm
  n_tot = n_tot + n_Spm
  mass_old(n_tot+1:n_tot+n_Su) = mUsquark
  n_tot = n_tot + n_Su
  mass_old(n_tot+1:n_tot+n_Sd) = mDsquark
  n_tot = n_tot + n_Sd
  mass_old(n_tot+1:n_tot+n_Sl) = mSlepton
  n_tot = n_tot + n_Sl
  mass_old(n_tot+1:n_tot+n_Sn) = mSneutrino

  !-----------------------------------------------------------------
  ! first setting of renormalization scale, if it is not yet fixed 
  ! somewhere else
  ! I take here the geometric mean of the stop masses
  !-----------------------------------------------------------------
  If (.Not.UseFixedScale) Then
   If (UseNewScale) Then
    mudim = CalcScale_from_Stops(Real(M2_U(3,3),dp), Real(M2_Q(3,3),dp)   &
               & , Y_u(3,3), A_u(3,3), vevSM, mu, gauge(2), gauge(1) )
    If (mudim.Lt.mZ**2) mudim = mZ**2
   Else
    mudim = Max(mZ**2 &
         &    , product_stop_masses(mUSquark, RUSquark, GenerationMixing))
   End If
   Call SetRGEScale(mudim)
   UseFixedScale = .False.
  End If
  !-----------------------------------------------------
  ! running of RGEs
  ! iterate entire process
  !-----------------------------------------------------
  A_l_mZ = A_l
  A_d_mZ = A_d
  A_u_mZ = A_u
  M2_E_mZ = M2_E
  M2_L_mZ = M2_L
  M2_D_mZ = M2_D
  M2_Q_mZ = M2_D
  M2_U_mZ = M2_U
  Mi_mZ = Mi
  mu_mZ = mu
  Y_l_mZ = Y_l
  Y_d_mZ = Y_d
  Y_u_mZ = Y_u
  tanb_mZ = tanb ! starting point, will be changed later
  Do j=1,niter
   !-------------------------------------
   ! boundary condition at the EW-scale
   !-------------------------------------
   If (WriteComment) Write(*,*) "Sugra",j
   Call Cpu_time(t1)

  !--------------------------------------------------------------
  ! neutrino mixings, if the dim-5 operator has non-zero entries
  !--------------------------------------------------------------
  If (Maxval(Abs(MnuL5)).Gt.0._dp) Then
   Call NeutrinoMasses(MnuL5, mf_nu, Unu, kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Deallocate(mass_old, mass_new, diff_m  )
    Return
   End If
  Else 
   Unu = id3C
  End If

  If (UseNewBoundaryEW) Then
   If (GenerationMixing) Then
    Call Switch_to_superCKM(Y_d_mZ, Y_u_mZ, A_d_mZ, A_u_mZ, M2_D_mZ, M2_Q_mZ, M2_U_mZ &
               & , Ad_s, Au_s, M2D_s, M2Q_s, M2U_s, .False.,ckm_out=uU_L )
    If (Maxval(Abs(MatNu)).Gt.0._dp) Then
     Call Switch_from_superPMNS(Y_l_mZ, A_l_mZ, M2_E_mZ, M2_L_mZ, MatNu &
               & , Al_s, M2E_s, M2L_s, .False. )
    Else
     Call Switch_from_superPMNS(Y_l_mZ, id3C, A_l_mZ, M2_E_mZ, M2_L_mZ &
               & , Al_s, M2E_s, M2L_s, .False. )
    End If
   Else
    Al_s = A_l
    Ad_s = A_d
    Au_s = A_u
    M2E_s = M2_E_mZ 
    M2L_s = M2_L_mZ
    M2D_s = M2_D_mZ
    M2Q_s = M2_Q_mZ
    M2U_s = M2_U_mZ
   End If
    Call BoundaryEW(j, mZ, tanb_mZ, Mi, Al_s, Ad_s, Au_s    &
     & , M2E_s, M2L_s, M2D_s, M2Q_s, M2U_s, mu_mZ, mP0(2)   &
     & , 0.1_dp*delta, GenerationMixing, .True., mZ2_run, mW2_run, g1, kont)

    uU_R = id3C
    uD_L = id3C
    uD_R = id3C
    uL_L = id3C
    uL_R = id3C
   Else 
    Call BoundaryEW(j, vevSM, mC, U, V, mN, N, mS02, RS0, mP02, RP0, mSpm, mSpm2&
     & , RSpm, mDsquark, mDsquark2, RDsquark, mUsquark, mUsquark2, RUsquark     &
     & , mSlepton, mSlepton2, RSlepton, mSneutrino2, RSneutrino, uU_L, uU_R     &
     & , uD_L, uD_R, uL_L, uL_R, Unu, mGlu, phase_glu, mZ2_run, mW2_run         &
     & , 0.1_dp*delta, g1, kont)
   End If

   If (kont.Ne.0) Then
    Iname = Iname - 1
    Deallocate(mass_old, mass_new, diff_m  )
    Return
   End If
   Call Cpu_time(t2)
   If (WriteComment) Write(*,*) "BoundaryEW",t2-t1
   !---------------------------------------------------
   ! now the running, 
   !---------------------------------------------------
   If ( (Sum(in_extpar(0,:)).Gt.0) .And. &
      & ((Sum(in_extpar(1:3,:))+Sum(in_extpar(21:24,:)) &
      &  +Sum(in_extpar(26,:))).Gt.0) ) Then  ! some low scale parameters 

    Call GToCouplings(g1, gauge, Y_l, Y_d, Y_u)
    Call ParametersToG(gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
          & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g2)
    
    mudim = Sqrt(GetRenormalizationScale())

    tz = Log(mudim/mZ)
    dt = tz / 50._dp

    If (tanb_in_at_Q) Then
     tanbQ = tanb_Q
    Else
     g58(1:57) = g1
     g58(58) = Log( tanb )

     Call odeint(g58, 58, 0._dp, tz,  0.1_dp * delta , dt, 0._dp, rge58, kont)
     tanbQ = Exp( g58(58) )
     tanb_Q = tanbQ
    End If

    tz = Log(mudim/mZ)
    dt = tz / 50._dp
    Call odeint(g2, 213, 0._dp, tz, 0.1_dp * delta, dt, 0._dp, rge213, kont)

    Call GToParameters(g2, gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
                  & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B)

    If (Sum(in_extpar(1,:)).Gt.0) Mi(1) = Cmplx(r_extpar(1),i_extpar(1),dp)
    If (Sum(in_extpar(2,:)).Gt.0) Mi(2) = Cmplx(r_extpar(2),i_extpar(2),dp)
    If (Sum(in_extpar(3,:)).Gt.0) Mi(3) = Cmplx(r_extpar(3),i_extpar(3),dp)

    If (Sum(in_extpar(21,:)).Gt.0) M2_H(1) = r_extpar(21)
    If (Sum(in_extpar(22,:)).Gt.0) M2_H(2) = r_extpar(22)
    If (Sum(in_extpar(23,:)).Gt.0) mu = Cmplx(r_extpar(23),i_extpar(23),dp)
    gauge(1) = Sqrt(3._dp/5._dp) * gauge(1)

    If (Sum(in_extpar(24,:)).Gt.0) &
     & Call Calcualte_MH2(delta, tanbQ, mudim, gauge, Y_l, Y_d, Y_u, Mi, A_l    &
                        & , A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, mu, mP02(2) &
                        & , .False., M2_H, B, kont)

    If (Sum(in_extpar(26,:)).Gt.0) &
     & Call Calcualte_MH2(delta, tanbQ, mudim, gauge, Y_l, Y_d, Y_u, Mi, A_l    &
                       & , A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, mu, P0(2)%m2 &
                       & , .True., M2_H, B, kont) 

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Deallocate(mass_old, mass_new, diff_m  )
     Return
    End If

    gauge(1) = Sqrt(5._dp/3._dp) * gauge(1)
    Call ParametersToG(gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
          & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g2)

    Call RunRGE_2_boundaries(kont, j+1, 0.1_dp*delta, mudim, g2, g0, mGUT)
    Call Cpu_time(t1)
    If (WriteComment) Write(*,*) "RunRGE2",t1-t2
   Else
    Call RunRGE(kont, 0.1_dp*delta, vevSM, g1, g0, mGUT)
    Call Cpu_time(t1)
    If (WriteComment) Write(*,*) "RunRGE",t1-t2
   End If

   If (kont.Ne.0) Then
    Iname = Iname - 1
    Deallocate(mass_old, mass_new, diff_m  )
    Return
   End If

   Call GToParameters(g0, gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
                  & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B)

   !----------------------------------------------------------------
   ! the RGE paper defines the Yukawas transposed to my conventions
   ! renormalize g_1 to g_Y
   !----------------------------------------------------------------
   Y_u = Transpose(Y_u)
   Y_d = Transpose(Y_d)
   Y_l = Transpose(Y_l)
   A_u = Transpose(A_u)
   A_d = Transpose(A_d)
   A_l = Transpose(A_l)
   gauge(1) = Sqrt(3._dp / 5._dp ) * gauge(1)

   If (Sum(in_extpar(26,:)).Gt.0) Then
    mp0(2) = P0(2)%m
    mp02(2) = P0(2)%m2
    Call LoopMassesMSSM_2(0.1_dp*delta, tanb_mZ, tanbQ, gauge, Y_l, Y_d, Y_u &
       & , Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, mu               &
       & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP02, RP0, mSpm        &
       & , mSpm2, RSpm, mDsquark, mDsquark2, RDsquark, mUsquark, mUsquark2   &
       & , RUsquark, mSlepton, mSlepton2, RSlepton, mSneutrino, mSneutrino2  &
       & , RSneutrino, mGlu, phase_glu, M2_H, B, kont)

    mu_loop = mu
    B_loop = B
    Call Cpu_time(t2)
    If (WriteComment) Write(*,*) "LoopMassesMSSM_2",t2-t1

   Else If (Sum(in_extpar(24,:)).Gt.0) Then
    Call LoopMassesMSSM_3(tanb_mZ, tanbQ, gauge, Y_l, Y_d, Y_u, Mi, A_l     &
       & , A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, mu, B, 0.1_dp*delta      &
       & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP0, mP02, RP0, mSpm  &
       & , mSpm2, RSpm, mDsquark, mDsquark2, RDsquark, mUsquark, mUsquark2  &
       & , RUsquark, mSlepton, mSlepton2, RSlepton, mSneutrino, mSneutrino2 &
       & , RSneutrino, mGlu, phase_glu, M2_H, kont)
    mu_loop = mu
    B_loop = B
    Call Cpu_time(t2)
    If (WriteComment) Write(*,*) "LoopMassesMSSM_3",t2-t1

   Else
    Call LoopMassesMSSM(0.1_dp*delta, tanb_mZ, tanb, gauge, Y_l, Y_d, Y_u &
     & , Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, phase_mu  &
     & , mu, B, j, uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                     &
     & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP0, mP02, RP0, mSpm  &
     & , mSpm2, RSpm, mDsquark, mDsquark2, RDsquark, mUsquark, mUsquark2  &
     & , RUsquark, mSlepton, mSlepton2, RSlepton, mSneutrino, mSneutrino2 &
     & , RSneutrino, mGlu, phase_glu, kont)
    tanb_Q = tanb

    g0(210) = Real(mu,dp)
    g0(211) = Aimag(mu)
    g0(212) = Real(B,dp)
    g0(213) = Aimag(B)
    Call Cpu_time(t2)
    If (WriteComment) Write(*,*) "LoopMassesMSSM",t2-t1
   End If

   If (kont.Ne.0) Then
    Iname = Iname - 1
    Deallocate(mass_old, mass_new, diff_m  )
    Return
   End If
   
   !-------------------
   ! comparing masses
   !-------------------
   mass_new(1) = mGlu
   n_tot = 1
   mass_new(n_tot+1:n_tot+n_C) = mC
   n_tot = n_tot + n_C
   mass_new(n_tot+1:n_tot+n_N) = mN
   n_tot = n_tot + n_N
   mass_new(n_tot+1:n_tot+n_S0) = mS0
   n_tot = n_tot + n_S0
   mass_new(n_tot+1:n_tot+n_P0) = mP0
   n_tot = n_tot + n_P0
   mass_new(n_tot+1:n_tot+n_Spm) = mSpm
   n_tot = n_tot + n_Spm
   mass_new(n_tot+1:n_tot+n_Su) = mUsquark
   n_tot = n_tot + n_Su
   mass_new(n_tot+1:n_tot+n_Sd) = mDsquark
   n_tot = n_tot + n_Sd
   mass_new(n_tot+1:n_tot+n_Sl) = mSlepton
   n_tot = n_tot + n_Sl
   mass_new(n_tot+1:n_tot+n_Sn) = mSneutrino

   diff_m = Abs(mass_new - mass_old)
   Where (Abs(mass_old).Gt.0._dp) diff_m = diff_m / Abs(mass_old)
 
   deltag0 = Maxval( diff_m )

   If (WriteComment) Write(*,*) "Sugra,Comparing",deltag0

   If ((deltag0.Lt.delta).And.(j.Gt.1)) Then ! require at least two iterations
    FoundResult = .True.
    B = B_loop
    mu = mu_loop
    Exit
   Else
    mass_old = mass_new
    !----------------------------------------------------------------
    ! recalculating massses at tree level; this is needed as input
    ! for BoundaryEW
    ! first running down to m_Z
    !----------------------------------------------------------------
    If (WriteComment) Write(*,*) "Sugra, Tree level masses",deltag0
    mudim = GetRenormalizationScale()
    tz = 0.5_dp * Log(mZ2/mudim)
    dt = tz / 100._dp

    !---------------------------------------------------------------
    ! recalculating scale for later use here, as I am using the 
    ! tree-level stop masses as default
    !---------------------------------------------------------------
    If (.Not.UseFixedScale) Then
     If (UseNewScale) Then
      mudim = CalcScale_from_Stops(Real(M2_U(3,3),dp), Real(M2_Q(3,3),dp)   &
                 & , Y_u(3,3), A_u(3,3), vevSM, mu, gauge(2), gauge(1) )
      If (mudim.Lt.mZ**2) mudim = mZ**2
     Else
      mudim = Max(mZ**2 &
           &    , product_stop_masses(mUSquark, RUSquark, GenerationMixing))
     End If
     Call SetRGEScale(mudim)
     UseFixedScale = .False.
    End If
    !-----------------------------------------------------------------
    ! and now the running, using the previous scale for consistency
    !-----------------------------------------------------------------
    Call odeint(g0, 213, 0._dp, tz, 0.1_dp * delta, dt, 0._dp, rge213, kont)

    Call GToParameters(g0, gauge_mZ, Y_l_mZ, Y_d_mZ, Y_u_mZ, Mi_mZ, A_l_mZ &
           & , A_d_mZ, A_u_mZ, M2_E_mZ, M2_L_mZ, M2_D_mZ, M2_Q_mZ, M2_U_mZ &
           & , M2_H_mZ, mu_mZ, B_mZ)
    Y_u_mZ = Transpose(Y_u_mZ)
    Y_d_mZ = Transpose(Y_d_mZ)
    Y_l_mZ = Transpose(Y_l_mZ)
    A_u_mZ = Transpose(A_u_mZ)
    A_d_mZ = Transpose(A_d_mZ)
    A_l_mZ = Transpose(A_l_mZ)
    gauge_mZ(1) = Sqrt(3._dp / 5._dp ) * gauge_mZ(1)

    vev = 2._dp * mW / gauge_mZ(2)
    vevSM(1) = vev / Sqrt(1._dp + tanb_mZ**2)
    vevSM(2) = tanb_mZ * vevSM(1)

    Abs_Mu2 = (M2_H_mZ(2) * tanb_MZ**2 - M2_H_mZ(1) ) / (1._dp - tanb_MZ**2) &
          & - 0.5_dp * mZ2
    If (Abs_Mu2.Lt.0._dp) Then
     Write (ErrCan,*) 'Warning from subroutine Sugra, |mu|^2 < 0 at m_Z'
     kont = -415
     Call AddError(415)

     Iname = Iname - 1
     Deallocate(mass_old, mass_new, diff_m  )
     Return
    End If
    Abs_Mu = Sqrt( Abs_Mu2 )
    mu_T = Abs_mu * phase_mu
    B_T = (M2_H_mZ(1) + M2_H_mZ(2) + 2._dp *  Abs_Mu2) * tanb_MZ / (1+tanb_MZ**2)

    Call TreeMassesMSSM2(gauge_mZ(1), gauge_mZ(2), vevSM, Mi_mZ(1), Mi_mZ(2)   &
     & , Mi_mZ(3), mu_T, B_T, tanb_mZ, M2_E_mZ, M2_L_mZ, A_l_mZ, Y_l_mZ        &
     & , M2_D_mZ, M2_U_mZ, M2_Q_mZ, A_d_mZ, A_u_mZ, Y_d_mZ, Y_u_mZ, uU_L, uU_R &
     & , uD_L, uD_R, uL_L, uL_R, mGlu_T, Phase_Glu_T, mC_T, mC2_T, U_T, V_T    &
     & , mN_T, mN2_T, N_T, mSneutrino_T, mSneutrino2_T, Rsneutrino_T           &
     & , mSlepton_T, mSlepton2_T, RSlepton_T, mDSquark_T, mDSquark2_T          &
     & , RDSquark_T, mUSquark_T, mUSquark2_T, RUSquark_T, mP0_T, mP02_T, RP0_T &
     & , mS0_T, mS02_T, RS0_T, mSpm_T, mSpm2_T, RSpm_T, mZ2_run, mW2_run       &
     & , GenerationMixing, kont, .False., .False.)

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Deallocate(mass_old, mass_new, diff_m  )
     Return
    End If

    !-----------------------------------------------------------------------
    ! the following variables need to be defined if m^2_A0 < 0,
    ! otherwise the NAG compiler creates problems
    !-----------------------------------------------------------------------
    If (mP02_T(2).Lt.0._dp) Then
     mspm2_t=0
     ms02_t = 0
    End If
    !-------------------------------------------------------------------
    ! checking if at tree level all masses squared are positiv and above
    ! approximately 0.9*m_Z at m_Z (the latter is for numerical stability),
    ! if not, the on-shell masses will be used 
    !-------------------------------------------------------------------
    If (Min(Minval(mUsquark2_T), Minval(mDSquark2_T), Minval(mSlepton2_T)   &
     &    , Minval(mSneutrino2_T), Minval(mS02_T), Minval(mP02_T)           &
     &    , Minval(mSpm2_T)).Gt. 5000._dp ) Then ! more than ~mZ^2/2
     mGlu = mGlu_T
     Phase_Glu = Phase_Glu_T
     mC = mC_T
     mC2 = mC2_T
     U = U_T
     V = V_T
     mN = mN_T
     mN2 = mN2_T
     N = N_T
     mSneutrino = mSneutrino_T
     mSneutrino2 = mSneutrino2_T
     Rsneutrino = Rsneutrino_T
     mSlepton = mSlepton_T
     mSlepton2 = mSlepton2_T
     RSlepton = RSlepton_T
     mDSquark = mDSquark_T
     mDSquark2 = mDSquark2_T
     RDSquark = RDSquark_T
     mUSquark = mUSquark_T
     mUSquark2 = mUSquark2_T
     RUSquark = RUSquark_T
     mP0 = mP0_T
     mP02 = mP02_T
     RP0 = RP0_T
     mS0 = mS0_T
     mS02 = mS02_T
     RS0 = RS0_T
     mSpm = mSpm_T
     mSpm2 = mSpm2_T
     RSpm = RSpm_T
     YukScen = 1
    Else
     YukScen = 2
     kont = 0
    End If
   End If

  End Do

  If (.Not. FoundResult ) Then
   Write (ErrCan,*) 'Warning from subroutine Sugra, no convergence has'
   Write (ErrCan,*) 'has been found after',niter,' iterations.'
   Write (ErrCan,*) 'required delta',delta 
   Write (ErrCan,*) 'found delta',deltag0
   kont = -412
   Call AddError(412)
  End If

  Deallocate( mass_old, mass_new, diff_m )

  Iname = Iname - 1

 End Subroutine Sugra

End Module SugraRuns
