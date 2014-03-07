Module NMSSM_tools

!---------------------------------
! loading necessary modules
!---------------------------------
Use Control
Use LoopFunctions 
Use RGEs
Use StandardModel
Use LoopCouplings
Use EplusEminusProduction
Use Model_Data
Use MSSM_Data
Use NMSSM_Data
Use BranchingRatios
Use LoopMasses
Use SugraRuns

Contains


 Subroutine Model_NMSSM(m32, Grav_fac, F_GMSB, Ecms, Pm, Pp, ISR, Beam      &
   & , SigSup , SigSdown, SigSle, SigSn, SigC, SigChi0, SigS0, SigSP, SigHp &
   & , kont)
 Implicit None

  !---------------------------------------
  ! input / output
  !---------------------------------------
  Integer, Intent (inout) :: kont
  !-------------------------------------------------------------------
  ! are needed in case of light gravitinos, e.g. in GMSB like models
  !-------------------------------------------------------------------
  Real(dp), Intent(in) :: m32, Grav_fac, F_GMSB
  !--------------------------------
  ! cross section calculation
  !--------------------------------
  Real(dp), Intent(in) :: Ecms(:), Pm(:), Pp(:)
  Logical, Intent(in) :: ISR(:), Beam(:)
  Real(dp), Intent(out) :: SigSup(:,:,:) , SigSdown(:,:,:), SigSle(:,:,:) &
         & , SigSn(:,:,:), SigC(:,:,:), SigChi0(:,:,:), SigS0(:,:)        &
         & , SigSP(:,:,:), SigHp(:,:,:)

  !---------------------------------------
  ! local variables
  !---------------------------------------
  Real(dp) :: gp, g, vev, sinW2, gauge_mZ(3), ht &
         & , hb, htau, Scale_Q, tz, dt, g8(8), g9(9), g6(6)
  Real(dp) :: delta = 1.e-5_dp, logQ, alphas_DR
  Complex(dp) :: mu_in
  Integer :: i1 
  Real(dp), Parameter :: &
    & as2loop = 1._dp / 24._dp + 2011._dp * oo32Pi2 / 12._dp          &
    &         + Log2 / 12._dp - oo8Pi2 * Zeta3                        &
    & , log2loop_a = 123._dp * oo32Pi2, log2loop_b = 33._dp * oo32Pi2
 !-------------------------------
 ! widths and branching ratios
 !-------------------------------
  Real(dp) ::  epsI, deltaM, ratioWoM
  Logical :: CalcTBD
  !--------------------------------
  ! cross section calculation
  !--------------------------------
  Integer :: p_max
  Complex(dp) :: Ylp(3,3)
 !------------------------------------------------------
 ! internal information on particle identities
 !------------------------------------------------------
  Integer :: id_gl, id_ph, id_grav
  Integer, Dimension(1) :: id_Z, id_W
  Integer, Dimension(3) :: id_nu, id_l, id_d, id_u


  Iname = Iname + 1 
  NameOfUnit(Iname) = "Model_NMSSM"

  Call Initialize_NMSSM(GenerationMixing, id_gl, id_ph, id_Z, id_W &
               & , id_nu, id_l, id_d, id_u, id_grav)
  kont = 0

  sinW2 = 1._dp - mW2/mZ2 
  
  alpha_mZ = Alpha_MSbar(mZ, mW)
  gauge_mZ(1) = Sqrt(4._dp * pi * alpha_mZ/(1._dp - sinW2))
  gauge_mZ(2) = Sqrt(4._dp * pi * alpha_mZ/sinW2) 
  alphas_DR =  AlphaS_mZ / (1._dp - oo4pi *AlphaS_mZ )
  gauge_mZ(3) = Sqrt(4._dp * pi * alphas_DR)

  tanb_mZ = tanb

  vev = Sqrt(mZ2 * (1._dp - sinW2) *sinW2 / (pi * alpha_mZ)) 
  vevSM(1) = vev / Sqrt(1._dp + tanb**2)
  vevSM(2) = tanb * vevSM(1)
! wird in InputOutput definiert zu vP = sqrt2 * mu_eff/lambda

  Y_l = 0._dp
  Y_d = 0._dp
  Y_u = 0._dp
  A_l = 0._dp
  A_d = 0._dp 
  A_u = 0._dp 

  Scale_Q = Sqrt(GetRenormalizationScale())

  htau = Sqrt2 * mf_l_mZ(3) / vevSM(1)
  hb = Sqrt2 * mf_d_mZ(3) / vevSM(1)
  ht = Sqrt2 * mf_u_mZ(3) / vevSM(2)
  logQ = 2._dp * Log(mZ/mf_u(3))
  ht = ht * (1._dp - alphas_DR / (3._dp*pi) * (5._dp +3._dp * LogQ  &
       &   + (as2loop + log2loop_a * logQ                      &
       &                         + log2loop_b * logQ**2) *4*pi*alphas_DR))

  !-------------------------------------------------------------
  ! running up within MSSM to get first approximation at Q_EWSB
  !-------------------------------------------------------------
   g6(1:3) = gauge_mZ
   g6(1) = Sqrt(5._dp/3._dp) * g6(1)
   g6(4) = htau
   g6(5) = hb
   g6(6) = ht

   tz = Log(Scale_Q/mZ)
   dt = tz / 50._dp
   Call odeint(g6, 6, 0._dp, tz, 0.1_dp*delta, dt, 0._dp, rge6, kont)
   !-----------------------------------------------------------------------
   ! adding NMSSM couplings, run down to m_Z to get the NMSSM couplings at 
   ! m_Z, repeating this a second time including tan(beta)
   !-----------------------------------------------------------------------
   g8(1:6) = g6
   g8(7) = h0
   g8(8) = 0.5_dp * lam ! using convention of hep-ph/9505326

   Call odeint(g8, 8, tz, 0._dp, 0.1_dp*delta, - dt, 0._dp, rge8_NMSSM, kont)

   g8(1:3) = gauge_mZ
   g8(1) = Sqrt(5._dp/3._dp) * g8(1) 
   g8(4) = htau
   g8(5) = hb
   g8(6) = ht
   Call odeint(g8, 8, 0._dp, tz, 0.1_dp*delta, dt, 0._dp, rge8_NMSSM, kont)

   g8(7) = h0
   g8(8) = 0.5_dp * lam ! using convention of hep-ph/9505326

   Call odeint(g8, 8, tz, 0._dp, 0.1_dp*delta, - dt, 0._dp, rge8_NMSSM, kont)

   g9(1:3) = gauge_mZ
   g9(1) = Sqrt(5._dp/3._dp) * g9(1)
   g9(4) = htau
   g9(5) = hb
   g9(6) = ht
   g9(7:8) = g8(7:8) 
   g9(9) = Log(tanb)
   Call odeint(g9, 9, 0._dp, tz, 0.1_dp*delta, dt, 0._dp, rge9_NMSSM, kont)

   gp = Sqrt(3._dp/5._dp) * g9(1)
   g = g9(2)
   gauge(2:3) = g9(2:3)
   gauge(1) = gp
   Y_l(3,3) = g9(4)
   Y_d(3,3) = g9(5)
   Y_u(3,3) = g9(6)
 
   A_u(3,3) = Y_u(3,3) * AoY_u(3,3)
   A_d(3,3) = Y_d(3,3) * AoY_d(3,3)
   A_l(3,3) = Y_l(3,3) * AoY_l(3,3)

   tanb_Q = Exp(g9(9))

   vevSM(1) = vev / Sqrt(1._dp + tanb_Q**2)
   vevSM(2) = tanb_Q * vevSM(1)
 
   mu_in = 0._dp

   If (GenerationMixing) Then
    i1 = GetYukawaScheme()
    If (i1.Eq.1) Then
     Y_u = Matmul(Transpose(CKM),Y_u) 
    Else
     Y_d = Matmul(Conjg(CKM),Y_d) 
    End If
   End If

 ! calculates all SusyMasses in the NMSSM at tree level except for
 ! neutral Higgs bosons where the 1-loop effective potential is used
 ! false = Higgs bosons treelevel 
 ! true = Higgs bosons 1loop effective potential

 If (.Not.External_Spectrum) &
  & Call TreeMassesNMSSM(gp, g, vevSM, vP, Mi(1), Mi(2), Mi(3), mu_in, B, h0  &
            &, lam, A_h0, A_lam, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q    &
            &, A_d, A_u, Y_d, Y_u, Glu%m, PhaseGlu, ChiPm%m, ChiPm%m2, U, V, Chi05%m      &
            &, Chi05%m2, N5, Sneut%m, Sneut%m2, Rsneut, Slepton%m, Slepton%m2      &
            &, RSlepton, Sdown%m, Sdown%m2, RSdown, Sup%m, Sup%m2, RSup        &
            &, P03%m, P03%m2, RP03, S03%m, S03%m2, RS03, Spm%m, Spm%m2 &
            &, RSpm, GenerationMixing, kont, .True.)
  Glu%m2 = Glu%m**2

  If (kont.Ne.0) Then
   Iname = Iname - 1
   Return 
  End If

    If (GenerationMixing) Then
     Call FermionMass(Y_l,vevSM(1),mf_l,uL_L,uL_R,kont)
     Call FermionMass(Y_d,vevSM(1),mf_d,uD_L,uD_R,kont)
     Call FermionMass(Y_u,vevSM(2),mf_u,uU_L,uU_R,kont)
    Else
     uL_L = id3C
     uL_R = id3C
     uD_L = id3C
     uD_R = id3C
     uU_L = id3C
     uU_R = id3C
    End If

   !----------------------------------------------
   ! reorder state identification if necessary
   !----------------------------------------------
   If (.Not.GenerationMixing) Then
    Call Swap_Order_Sf(RSlepton(1,1), Slepton(1)%id, Slepton(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSlepton(3,3), Slepton(3)%id, Slepton(4)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSdown(1,1), Sdown(1)%id, Sdown(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSdown(3,3), Sdown(3)%id, Sdown(4)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSup(1,1), Sup(1)%id, Sup(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSup(3,3), Sup(3)%id, Sup(4)%id, id_p, c_name)    
   End If
  !----------------------------------------
  ! decays of SUSY and Higgs particles
  !----------------------------------------
  If ((L_BR).And.(kont.Eq.0)) Then
   ! relative precision for the calculation of phase space integrals
   epsI = 1.e-5_dp
   deltaM = 1.e-3_dp ! if mass/phase space is smaller, than mass is treated as 0
   CalcTBD = .False. ! if .True. than calculation of 3-body decays is enforced
   ratioWoM = 0._dp ! 1.e-4_dp

   Call CalculateBR(n_nu, id_nu, n_l, id_l, n_d, id_d, n_u, id_u, n_Z, id_Z    &
     & , n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_g, n_s0, n_p0, n_Spm &
     & , id_grav, id_gl, id_ph, gauge, Glu, PhaseGlu, Chipm, U, V, Chi05, N5   &
     & , Sneut, RSneut, Slepton, RSlepton, Sup, RSup, Sdown, RSdown, uL_L      &
     & , uL_R, uD_L, uD_R, uU_L, uU_R, S03, RS03, P03, RP03, Spm, RSpm, epsI   &
     & , deltaM, CalcTBD, ratioWoM, Y_d, A_d, Y_l, A_l, Y_u, A_u, mu, h0, A_h0 &
     & , lam, A_lam, vevSM, vP, F_Gmsb, m32, grav_fac)

  End If

  !-------------------------------------------------------------
  ! cross section for e+ e- annihilation
  ! The following input quantities can be specified either in the file 
  ! CrossSections.in or in the file LesHouches.in:
  !                    Ecms .... c.m.s enerergy in GeV
  !                   Pm ...... degree of longitudinal polarization of incoming
  !                             electron
  !                   Pp ...... degree of longitudinal polarization of incoming
  !                             positron
  !                   ISR ..... if .TRUE. then the effect of initial state
  !                             radiation will be included
  ! In the case that none of these files exist, the following
  ! default values are used: Ecms = 500 GeV, Pm = Pp = 0, ISR = .TRUE.
  !-------------------------------------------------------------
  If ((L_CS).And.(kont.Eq.0)) Then
    Ylp = Y_l / gauge(2)
    p_max = Size(Pm)

    Do i1=1,p_max
     If (Ecms(i1).Eq.0._dp) Exit
     Call CalculateCrossSections(Ecms(i1), Pm(i1), Pp(i1), ISR(i1), Beam(i1)   &
            & , "Tesla800", Sup%m, RSup, mf_u, Sdown%m, RSdown, mf_d, glu%m    &
            & , SigSup(i1,:,:), SigSdown(i1,:,:), Slepton%m, RSlepton, Ylp     &
            & , Sneut%m, RSneut, SigSle(i1,:,:), SigSn(i1,:,:), ChiPm%m, U, V  &
            & , Chi05%m, N5, SigC(i1,:,:), SigChi0(i1,1:5,1:5), S03%m, RS03    &
            & , vevSM, P03%m, RP03, Spm%m, RSpm, SigS0(i1,1:3)                 &
            & , SigSP(i1,1:3,1:2), SigHp(i1,1,1) )
    End Do
   End If

  Iname = Iname - 1
 

 Contains

  Subroutine Swap_Order_Sf(Rij, id1, id2, id_p, names)
  Implicit None
   Complex(dp), Intent(in) :: Rij
   Integer, Intent(in) :: id1, id2
   Integer, Intent(inout) :: id_p(:)
   Character(len=12), Intent(inout) :: names(:)

   Integer :: k(2)
   Character(len=12) :: nam(2)
   !--------------------------------------------------------------------
   ! changes name and particle id, in case that the lighter of two
   ! sfermions is a right-sfermion
   !--------------------------------------------------------------------
   If (Abs(Rij).Lt.0.5_dp) Then
    k = id_p(id1:id1+1)
    id_p(id1:id1+1) = id_p(id2:id2+1)
    id_p(id2:id2+1) = k

    nam = names(id1:id1+1)
    names(id1:id1+1) = names(id2:id2+1)
    names(id2:id2+1) = nam
   End If

  End Subroutine Swap_Order_Sf

 End Subroutine Model_NMSSM

End Module NMSSM_tools
