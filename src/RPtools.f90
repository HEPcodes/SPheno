Module RPtools

Use Control
Use StandardModel
Use LoopCouplings
Use EplusEminusProduction
Use MSSM_data
Use RP_data
Use Experiment
Use BranchingRatios
Use LoopMasses
Use SugraRuns

 Interface Fit_Neutrino_Data
  Module Procedure Fit_Neutrino_Data_sp
 End Interface

 Interface ParMin
  Module Procedure ParMin_sp
 End Interface

Contains


 Subroutine Calculate_Bi(mu, eps, vevL, vevSM, gp, g, M2L, B)
 !----------------------------------------------------------------------
 ! calculates the R-parity violating B-parameters [GeV^2]
 ! written by Werner Porod, 18.08.05
 !----------------------------------------------------------------------
  Implicit None
   Complex(dp), Intent(in) :: mu, eps(3) ! bilinear parameters of Superpotential
   Complex(dp), Intent(in) :: M2L(3,3)   ! left slepton masses
   Real(dp), Intent(in) :: gp, g         ! U(1) & SU(2) gauge couplings
   Real(dp), Intent(in) :: vevSM(2)      ! (v_d, v_u)
   Real(dp), Intent(in) :: vevL(3)       ! (v_1, v_2, v_3)
   Complex(dp), Intent(out) :: B(3)      ! B-parameters [GeV^2]

   Integer :: i1, i2
   Real(dp) :: DTerm
   Complex(dp) :: SumC

   DTerm = 0.125_dp * (g**2+gp**2)  &
       & * ( vevSM(1)**2 - vevSM(2)**2 + Dot_product(vevL,vevL))

   Do i1=1,3
    sumC = - (M2L(i1,i1)+Dterm) * vevL(i1) + vevSM(1) * Conjg(mu) * eps(i1)
    Do i2=1,3
     sumC = sumC - Conjg( eps(i1) ) * eps(i2) * vevL(i2)
     If (i1.Ne.i2) sumC = sumC - vevL(i2) * M2_L(i2,i1)
    End Do
    B(i1) = sumC / vevSM(2) 
  End Do
   
 End Subroutine Calculate_Bi


 Subroutine Calculate_RP_Parameters(delta, eps, vevL, B_i, Lambda, RSpm, RS0 &
                                 & , RP0, U, V, mN7, N7, kont)
 Implicit None
  Real(dp), Intent(in) :: delta     ! numerical precision for RGE running
  Integer , Intent (inout) :: kont  ! checks if everyhint is o.k. = 0

  Real(dp), Intent(out) :: vevL(3)      ! sneutrino vevs
  Complex(dp), Intent(out) :: eps(3)    ! superpotential parameters epsilon_i
  Complex(dp), Intent(out) :: B_i(3)    ! soft SUSY parameters B_i [GeV^2]
  Real(dp), Intent(out) :: Lambda(3)    ! RP lambda vector = mu v_i + eps_i v_d
  Complex(dp), Intent(out) :: RSpm(8,8) ! charged scalar mixing matrix
  Real(dp), Intent(out) :: RS0(5,5)     ! neutral scalar mixing matrix
  Real(dp), Intent(out) :: RP0(5,5)     ! neutral pseudoscalar mixing matrix
  Complex(dp), Intent(out) :: U(5,5), V(5,5) ! chargino mixing matrices
  Real(dp), Intent(out) :: mN7(7)    ! loop corrected neutrino/neutralino masses
  Complex(dp), Intent(out) :: N7(7,7)   ! neutralino/neutrino mixing matrix

  Integer :: count, isol
  Real(dp) :: g0(213), mudim, tz, dt, gauge_mZ(3), M2H_mZ(3)
  Complex(dp), Dimension(3,3) :: y_l_mZ,  y_d_mZ, y_u_mZ , Al_mZ, Ad_mZ &
           &, Au_mZ,M2E_mZ, M2L_mZ, M2D_mZ, M2Q_mZ, M2U_mZ
  Complex(dp) :: Mi_mZ(3), B_mZ, B_4(4), bi(4)
 ! neutrino constraints
  Real(dp) :: Lam_Sq, m_sq, m2_atm_rp, m2_sol_rp, tan2_sol, tan2_atm &
   & , Ue32, tan2_sol_opt, tan2_atm_opt, epsT12, epsT22, Ue3_opt
  Logical :: check
 ! RP-masses + mixing
  Real(dp) :: mGlu, mC(5), mC2(5), mSdown(6), mSdown2(6), mSup(6), mSup2(6)   &
    & , mP0(5), mP02(5), mS0(5), mS02(5), mSpm(8), mSpm2(8) &
    & , mN(7), mN1L(7), mN2(7), mN1L2(7)
  Complex(dp) :: PhaseGlu, RSdown(6,6), RSup(6,6), N(7,7), N1L(7,7), lam(3), epst(3)
!  Logical, Save :: WriteOut = .False.
  Logical, Save :: WriteOut = .True.

  !------------------------------------------
  ! check if neutrino data are consistent
  !------------------------------------------
  If ((m2_atm.Lt.m2_atm_min).Or.(m2_atm.Gt.m2_atm_max)) &
      & m2_atm = 0.5_dp * (m2_atm_min +m2_atm_max)
  If ((m2_sol.Lt.m2_sol_min).Or.(m2_sol.Gt.m2_sol_max)) &
      & m2_sol = 0.5_dp * (m2_sol_min +m2_sol_max)
  Ue3_opt = 0.5_dp * (Ue32_max + Ue32_min)

  !-------------------------------------------------
  ! evolve the parameters down to m_Z if necessary 
  !-------------------------------------------------
  mudim = GetRenormalizationScale()
  Call ParametersToG(gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
                    &,M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g0)

  If (mudim.Ne.mZ2) Then
 
   tz = 0.5_dp * Log(mZ2/mudim)
   dt = tz / 100._dp
   g0(1) = Sqrt(5._dp / 3._dp ) * g0(1)
 
   Call odeint(g0, 213, 0._dp, tz, delta, dt, 0._dp, rge213, kont)
   g0(1) = Sqrt(3._dp / 5._dp ) * g0(1)
 
  End If
  
  Call GToParameters(g0,gauge_mZ, y_l_mZ,  y_d_mZ, y_u_mZ, Mi_mZ, Al_mZ, Ad_mZ &
          &, Au_mZ,M2E_mZ, M2L_mZ, M2D_mZ, M2Q_mZ, M2U_mZ, M2H_mZ, mu_mZ, B_mZ)

  kont = 0

  !--------------------------------------------
  ! try to find a consist set of RP-parameters
  !--------------------------------------------
  count = 0
  !-------------------------------------------------
  ! first set of masses to start with
  ! needed to get a first estimate of RP parameters
  !-------------------------------------------------
  tan2_atm_opt = 0.5_dp * (tan2_atm_min + tan2_atm_max)
  tan2_sol_opt = 0.5_dp * (tan2_sol_min + tan2_sol_max)

  eps = 1.e-4_dp * Abs(mu_mZ)
  eps(2) = -eps(2)

  bi(1) = mu_mZ
  bi(2:4) = eps
  b_4 = B_mZ

  Lam_sq = 4._dp * Sqrt(m2_atm)                                               &
         & * Real(Mi_mZ(1)*Mi_mZ(2)*mu_mZ**2                                  &
         &    -0.5*vevSM(1)*vevSM(2)*Mi_mZ(2)*mu_mZ*(g0(1)**2+g0(2)**2),dp)   &
         &       / (g0(1)**2*Real(Mi_mZ(2),dp) + g0(2)**2*Real(Mi_mZ(1),dp))
  Lam_Sq = Abs(Lam_Sq)
  Lambda(2:3) = Sqrt( 0.5_dp *(1._dp-Ue3_opt) *Lam_sq)
  Lambda(1) = Sqrt(Lam_sq-2._dp*Lambda(2)**2)

  vevL = (Lambda - eps*vevSM(1)) / mu_mZ
  !--------------------------------------
  ! tree level masses
  !--------------------------------------
  Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4) )

  Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                    &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                 &
                    &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ   &
                    &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N              &
                    &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup             &
                    &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm      &
                    &, GenerationMixing, kont)
  !----------------------------
  ! approx. squared solar mass
  !-----------------------------
   lam = Lambda
   Call CalculateEpsTilde_from_Eps(epsT, eps, Lam)

   epsT12 = Abs(epsT(1))**2
   epsT22 = Abs(epsT(2))**2

   m_Sq = 2._dp * oo16pi2 * (epsT12+epsT22) / Abs(mu_mZ)**2                   &
        &       * ( 3._dp*mf_d_mZ(3)*Rsdown(5,5)*Rsdown(5,6)*Y_d(3,3)**2      &
        &                *Log(mSdown2(6)/mSdown2(5))                          &
        &         + mf_l_mZ(3)* Rslepton(5,5)*Rslepton(5,6)*Y_l(3,3)**2       &
        &                     * Log(Slepton(6)%m2/Slepton(5)%m2))

   m_Sq = m_Sq**2

   If (m_Sq.Lt.m2_sol_min) eps = (m2_sol/m_sq)**0.25_dp * eps
   If (m_Sq.Gt.m2_sol_max) eps = (m2_sol/m_sq)**0.25_dp * eps

   bi(2:4) = eps
   vevL = (Lambda - eps*vevSM(1)) / mu_mZ 
   !-------------------------------------------------------------
   ! recalculating B_i to fulfill tadpoles
   !-------------------------------------------------------------
   Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4) )

  !--------------------------------------
  !recalculation of tree level masses
  !--------------------------------------
   Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                    &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                    &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                    &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                    &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                    &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                    &, GenerationMixing, kont)
  !-----------------------------------------------------
  ! the iteration to get a consistent set of parameters
  !-----------------------------------------------------
  isol=100
  Do count = 0,isol
   If (WriteOut) Write(Errcan,*) " "
   If (WriteOut) Write(Errcan,*) "Step",count

   Call NeutralinoMass_Loop_RP(g0(1), g0(2), Y_d_mZ, Y_l_mZ, Y_u_mZ, vevSM     &
          & , vevL, Mi_mZ(1), Mi_mZ(2), mu_mZ, eps, mC, mC2, U, V, mSup2, RSup &
          & , mSdown2, RSdown, mS02, RS0, mP02, RP0, mSpm2, RSpm, uD_L, uD_R   &
          & , uU_L, uU_R, mN, mN2, N, mN1L, mN1L2, N1L, kont, .False.)

   m2_sol_rp = mN1L2(2)-mN1L2(1)
   m2_atm_rp = mN1L2(3)-mN1L2(1)
   If (WriteOut) Then
    Write(Errcan,*) "m^2_atm,m^2_sol",m2_atm_rp,m2_sol_rp
    Write(Errcan,*) "           tree",mN2(3)-mN2(1),mN2(2)-mN2(1)
    Write(errcan,*) m2_atm_min,m2_atm_max
    Write(Errcan,*) (m2_atm_rp.Lt.m2_atm_min),(m2_atm_rp.Gt.m2_atm_max) &
                &   ,(m2_sol_rp.Lt.m2_sol_min),(m2_sol_rp.Gt.m2_sol_max)
   End If
   !------------------------------------------------
   ! checking experimental data, first the masses
   !------------------------------------------------
   check = .True. 
   If ((m2_atm_rp.Lt.m2_atm_min).Or.(m2_atm_rp.Gt.m2_atm_max)) Then
    check = .False.

    Lambda = (m2_atm/m2_atm_rp)**0.25_dp * Lambda
   End If

   If ((m2_sol_rp.Lt.m2_sol_min).Or.(m2_sol_rp.Gt.m2_sol_max)) Then
    check = .False.
    eps = (m2_sol/m2_sol_rp)**0.25_dp * eps
    bi(2:4) = eps
   End If

   If (.Not.check) Then
    vevL = (Lambda - eps*vevSM(1)) / mu_mZ 
    Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4))
    Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                     &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                     &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                     &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                     &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                     &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                     &, GenerationMixing, kont)
    Cycle  !stop here and try the next iteration
   End If
   !-------------------------------------
   ! and now atmospheric and solar angle
   !-------------------------------------
   tan2_atm = Abs(N1L(3,6)/N1L(3,7))**2
   tan2_sol = Abs(N1L(2,5)/N1L(1,5))**2

   If ((tan2_atm.Lt.tan2_atm_min).Or.(tan2_atm.Gt.tan2_atm_max)) Then
    check = .False.

    Lam_sq = Lambda(2)**2 + Lambda(3)**2
    Lambda(2) = Sqrt(tan2_atm_opt/(tan2_atm+tan2_atm_opt)) * Lambda(2)
    Lambda(3) = Sqrt(tan2_atm/(tan2_atm+tan2_atm_opt)) * Lambda(3)
    Lam_sq = Lam_sq / (Lambda(2)**2 + Lambda(3)**2)
    Lambda(2:3) =  Lambda(2:3) * Sqrt(Lam_sq) 

   End If

   If ((tan2_sol.Lt.tan2_sol_min).Or.(tan2_sol.Gt.tan2_sol_max)) Then
    check = .False.

    lam = Cmplx(lambda,0._dp,dp)
    Call CalculateEpsTilde_from_Eps(epsT,eps,lam)
    epsT12 = Abs(epsT(1))**2 + Abs(epsT(2))**2

    epsT(1) = Sqrt(tan2_sol_opt/(tan2_sol+tan2_sol_opt)) * epsT(1)
    epsT(2) = Sqrt(tan2_sol/(tan2_sol+tan2_sol_opt)) * epsT(2)
 
    epsT12 = epsT12 / (Abs(epsT(1))**2 + Abs(epsT(2))**2)
    epsT(1:2) = epsT(1:2) * Sqrt(epsT12)

    Call CalculateEps_from_EpsTilde(eps,epsT,lam)

    bi(2:4) = eps
   End If

   If (.Not.check) Then
    vevL = (Lambda - eps*vevSM(1)) / mu_mZ 
    Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4))
    Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                     &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                     &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                     &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                     &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                     &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                     &, GenerationMixing, kont)
    Cycle
   End If
   !-----------------------------------
   ! and now the reactor constraint
   !-----------------------------------
   Ue32 = Abs(N1l(3,5))**2

   If (Ue32.Gt.Ue32_max) Then
    check = .False.
    Lambda(1) = Lambda(1) / 2._dp
    Lam_Sq = Dot_product(Lambda, lambda)

    vevL = (Lambda - eps*vevSM(1)) / mu_mZ 
    Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4))
   End If
   !-----------------------------------
   ! leave loop if everything is fine
   !-----------------------------------
   If (check) Exit
   !---------------------------------------------------
   ! else recalculate tree level masses
   !---------------------------------------------------
   Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                    &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                    &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                    &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                    &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                    &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                    &, GenerationMixing, kont)

  End Do

  If(.Not.check)Then 
    Write(ErrCan,*)'Error, no solution found',isol
    kont = -isol
!    Write(*,*)'Error, no solution found',isol
!  Else
!    Write(*,*)'solution found after:',count,Sqrt(Sum(Abs(eps)**2))/Real(mu_mZ,dp) 
  End If

  mN7 = mN1L
  N7 = N1L
  B_i = B_4(2:4)

 End Subroutine Calculate_RP_Parameters

!---------------------------------------------------
! now the complex version
!---------------------------------------------------

 Subroutine Calculate_RP_ParametersC(delta, eps, vevL, B_i, Lambda, RSpm, RS0 &
                                 & , RP0, U, V, mN7, N7, kont, phases_eps)
 Implicit None
  Real(dp), Intent(in) :: delta     ! numerical precision for RGE running
  Integer , Intent (inout) :: kont  ! checks if everyhint is o.k. = 0
  Complex(dp), Intent(in), Optional :: phases_eps(3) ! in case the phases are given

  Real(dp), Intent(out) :: vevL(3)      ! sneutrino vevs
  Complex(dp), Intent(out) :: eps(3)    ! superpotential parameters epsilon_i
  Complex(dp), Intent(out) :: B_i(3)    ! soft SUSY parameters B_i [GeV^2]
  Complex(dp), Intent(out) :: Lambda(3) ! RP lambda vector = mu v_i + eps_i v_d
  Complex(dp), Intent(out) :: RSpm(8,8) ! charged scalar mixing matrix
  Real(dp), Intent(out) :: RS0(5,5)     ! neutral scalar mixing matrix
  Real(dp), Intent(out) :: RP0(5,5)     ! neutral pseudoscalar mixing matrix
  Complex(dp), Intent(out) :: U(5,5), V(5,5) ! chargino mixing matrices
  Real(dp), Intent(out) :: mN7(7)    ! loop corrected neutrino/neutralino masses
  Complex(dp), Intent(out) :: N7(7,7)   ! neutralino/neutrino mixing matrix

  Integer :: count, isol
  Real(dp) :: g0(213), mudim, tz, dt, gauge_mZ(3), M2H_mZ(3)
  Complex(dp), Dimension(3,3) :: y_l_mZ,  y_d_mZ, y_u_mZ , Al_mZ, Ad_mZ &
           &, Au_mZ,M2E_mZ, M2L_mZ, M2D_mZ, M2Q_mZ, M2U_mZ
  Complex(dp) :: Mi_mZ(3), B_mZ, B_4(4), bi(4)
 ! neutrino constraints
  Real(dp) :: Lam_Sq, m_sq, m2_atm_rp, m2_sol_rp, tan2_sol, tan2_atm &
   & , Ue32, tan2_sol_opt, tan2_atm_opt, epsT12, epsT22, Ue3_opt
  Logical :: check
 ! RP-masses + mixing
  Real(dp) :: mGlu, mC(5), mC2(5), mSdown(6), mSdown2(6), mSup(6), mSup2(6)   &
    & , mP0(5), mP02(5), mS0(5), mS02(5), mSpm(8), mSpm2(8) &
    & , mN(7), mN1L(7), mN2(7), mN1L2(7)
  Complex(dp) :: PhaseGlu, RSdown(6,6), RSup(6,6), N(7,7), N1L(7,7), lam(3), epst(3)
  Logical, Save :: WriteOut = .False.
  Logical :: l_phases

  Iname = Iname + 1
  NameOfUnit(Iname) = "Calculate_RP_ParametersC"

  !------------------------------------------
  ! check if neutrino data are consistent
  !------------------------------------------
  If ((m2_atm.Lt.m2_atm_min).Or.(m2_atm.Gt.m2_atm_max)) &
      & m2_atm = 0.5_dp * (m2_atm_min +m2_atm_max)
  If ((m2_sol.Lt.m2_sol_min).Or.(m2_sol.Gt.m2_sol_max)) &
      & m2_sol = 0.5_dp * (m2_sol_min +m2_sol_max)
  Ue3_opt = 0.5_dp * (Ue32_max + Ue32_min)
 !-------------------------------------------------
 ! evolve the parameters down to m_Z if necessary 
 !-------------------------------------------------
  mudim = GetRenormalizationScale()
  Call ParametersToG(gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
                    &,M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g0)

  If (mudim.Ne.mZ2) Then
 
   tz = 0.5_dp * Log(mZ2/mudim)
   dt = tz / 100._dp
   g0(1) = Sqrt(5._dp / 3._dp ) * g0(1)
 
   Call odeint(g0, 213, 0._dp, tz, delta, dt, 0._dp, rge213, kont)
   g0(1) = Sqrt(3._dp / 5._dp ) * g0(1)
 
  End If
  
  Call GToParameters(g0,gauge_mZ, y_l_mZ,  y_d_mZ, y_u_mZ, Mi_mZ, Al_mZ, Ad_mZ &
          &, Au_mZ,M2E_mZ, M2L_mZ, M2D_mZ, M2Q_mZ, M2U_mZ, M2H_mZ, mu_mZ, B_mZ)

  kont = 0

  !------------------------------------------
  ! check if there are phases
  !------------------------------------------
  If (Present(phases_eps).Or.(Aimag(mu_mz).Ne.0._dp)) Then
   l_phases = .True.  ! in this case we need the procedure
  Else
   l_phases = .False.  ! real parameters, tested procedure works
  End If

  !--------------------------------------------
  ! try to find a consist set of RP-parameters
  !--------------------------------------------
  count = 0
  !-------------------------------------------------
  ! first set of masses to start with
  ! needed to get a first estimate of RP parameters
  !-------------------------------------------------
  tan2_atm_opt = 0.5_dp * (tan2_atm_min + tan2_atm_max)
  tan2_sol_opt = 0.5_dp * (tan2_sol_min + tan2_sol_max)

  eps = 1.e-4_dp * Abs(mu_mZ)
  If (Present(phases_eps)) Then
   eps = eps * phases_eps
  Else
   eps(2) = -eps(2)
  End If

  bi(1) = mu_mZ
  bi(2:4) = eps
  b_4 = B_mZ

  Lam_sq = 4._dp * Sqrt(m2_atm)                                               &
         & * Real(Mi_mZ(1)*Mi_mZ(2)*mu_mZ**2                                  &
         &    -0.5*vevSM(1)*vevSM(2)*Mi_mZ(2)*mu_mZ*(g0(1)**2+g0(2)**2),dp)   &
         &       / (g0(1)**2*Real(Mi_mZ(2),dp) + g0(2)**2*Real(Mi_mZ(1),dp))
  Lam_Sq = Abs(Lam_Sq)
  
  If (l_phases) Then
   Call my_vi(mu_mZ, vevSM(1), eps(1), Ue3_opt * Lam_sq, vevL(1), Lambda(1) &
      &      , kont)
   If (kont.Ne.0._dp) Then ! let's take are more risky choice
    Ue3_opt = Ue32_max
    Call my_vi(mu_mZ, vevSM(1), eps(1), Ue3_opt * Lam_sq, vevL(1), Lambda(1) &
      &      , kont)
   End If
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   ! abuse of notation
   Lam_sq = 0.5_dp * (1._dp-Ue3_opt) * Lam_sq
   Call my_vi(mu_mZ, vevSM(1), eps(2), Lam_sq, vevL(2), Lambda(2), kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   Call my_vi(mu_mZ, vevSM(1), eps(3), Lam_sq, vevL(3), Lambda(3), kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If

  Else ! in this case the old procedure
   Lambda(2:3) = Sqrt( 0.5_dp *(1._dp-Ue3_opt) *Lam_sq)
   Lambda(1) = Sqrt(Lam_sq-2._dp*Abs(Lambda(2))**2)
   vevL = Lambda - eps*vevSM(1) / mu_mZ
  End If

  !--------------------------------------
  ! tree level masses
  !--------------------------------------
  Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4) )

  Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                    &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                 &
                    &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ   &
                    &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N              &
                    &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup             &
                    &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm      &
                    &, GenerationMixing, kont)
  !----------------------------
  ! approx. squared solar mass
  !-----------------------------
   Call CalculateEpsTilde_from_Eps(epsT, eps, Lambda)

   epsT12 = Abs(epsT(1))**2
   epsT22 = Abs(epsT(2))**2

   m_Sq = 2._dp * oo16pi2 * (epsT12+epsT22) / Abs(mu_mZ)**2                   &
        &       * ( 3._dp*mf_d_mZ(3)*Rsdown(5,5)*Rsdown(5,6)*Y_d(3,3)**2      &
        &                *Log(mSdown2(6)/mSdown2(5))                          &
        &         + mf_l_mZ(3)* Rslepton(5,5)*Rslepton(5,6)*Y_l(3,3)**2       &
        &                     * Log(Slepton(6)%m2/Slepton(5)%m2))

   m_Sq = m_Sq**2

   If (m_Sq.Lt.m2_sol_min) eps = (m2_sol/m_sq)**0.25_dp * eps
   If (m_Sq.Gt.m2_sol_max) eps = (m2_sol/m_sq)**0.25_dp * eps

   bi(2:4) = eps

   Lam_sq = Dot_product(Lambda,Lambda)
   Call my_vi(mu_mZ, vevSM(1), eps(1), Ue3_opt * Lam_sq, vevL(1), Lambda(1) &
      &      , kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   ! abuse of notation
   Lam_sq = 0.5_dp * (1._dp-Ue3_opt) * Lam_sq
   Call my_vi(mu_mZ, vevSM(1), eps(2), Lam_sq, vevL(2), Lambda(2), kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   Call my_vi(mu_mZ, vevSM(1), eps(3), Lam_sq, vevL(3), Lambda(3), kont)
   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If

   !-------------------------------------------------------------
   ! recalculating B_i to fulfill tadpoles
   !-------------------------------------------------------------
   Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4) )

  !--------------------------------------
  !recalculation of tree level masses
  !--------------------------------------
   Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                    &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                    &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                    &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                    &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                    &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                    &, GenerationMixing, kont)
  !-----------------------------------------------------
  ! the iteration to get a consistent set of parameters
  !-----------------------------------------------------
  isol=100
  Do count = 0,isol
   If (WriteOut) Write(Errcan,*) " "
   If (WriteOut) Write(Errcan,*) "Step",count

   Call NeutralinoMass_Loop_RP(g0(1), g0(2), Y_d_mZ, Y_l_mZ, Y_u_mZ, vevSM     &
          & , vevL, Mi_mZ(1), Mi_mZ(2), mu_mZ, eps, mC, mC2, U, V, mSup2, RSup &
          & , mSdown2, RSdown, mS02, RS0, mP02, RP0, mSpm2, RSpm, uD_L, uD_R   &
          & , uU_L, uU_R, mN, mN2, N, mN1L, mN1L2, N1L, kont, .False.)

   m2_sol_rp = mN1L2(2)-mN1L2(1)
   m2_atm_rp = mN1L2(3)-mN1L2(1)
   If (WriteOut) Then
    Write(Errcan,*) "m^2_atm,m^2_sol",m2_atm_rp,m2_sol_rp
    Write(Errcan,*) "               ",mN2(3)-mN2(1),mN2(2)-mN2(1)
    Write(errcan,*) m2_atm_min,m2_atm_max
    Write(Errcan,*) (m2_atm_rp.Lt.m2_atm_min),(m2_atm_rp.Gt.m2_atm_max) &
                &   ,(m2_sol_rp.Lt.m2_sol_min),(m2_sol_rp.Gt.m2_sol_max)
   End If
   !------------------------------------------------
   ! checking experimental data, first the masses
   !------------------------------------------------
   check = .True. 
   If ((m2_atm_rp.Lt.m2_atm_min).Or.(m2_atm_rp.Gt.m2_atm_max)) Then
    check = .False.

    Lambda = (m2_atm/m2_atm_rp)**0.25_dp * Lambda
   End If

   If ((m2_sol_rp.Lt.m2_sol_min).Or.(m2_sol_rp.Gt.m2_sol_max)) Then
    check = .False.
    eps = (m2_sol/m2_sol_rp)**0.25_dp * eps
    bi(2:4) = eps
   End If

   If (.Not.check) Then
    If (l_phases) Then
     Lam_sq = Dot_product(Lambda,Lambda)
     Call my_vi(mu_mZ, vevSM(1), eps(1), Ue3_opt * Lam_sq, vevL(1), Lambda(1) &
       &      , kont)
     If (kont.Ne.0) Then
      Iname = Iname - 1
      Return
     End If
     ! abuse of notation
     Lam_sq = 0.5_dp * (1._dp-Ue3_opt) * Lam_sq
     Call my_vi(mu_mZ, vevSM(1), eps(2), Lam_sq, vevL(2), Lambda(2), kont)
     If (kont.Ne.0) Then
      Iname = Iname - 1
      Return
     End If
     Call my_vi(mu_mZ, vevSM(1), eps(3), Lam_sq, vevL(3), Lambda(3), kont)
     If (kont.Ne.0) Then
      Iname = Iname - 1
      Return
     End If
    Else ! old procedure
     Lam_sq = Dot_product(Lambda,Lambda)
     vevL = (Lambda - eps*vevSM(1)) / mu_mZ 
    End If

    Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4))
    Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                     &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                     &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                     &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                     &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                     &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                     &, GenerationMixing, kont)
    Cycle  !stop here and try the next iteration
   End If
   !-------------------------------------
   ! and now atmospheric and solar angle
   !-------------------------------------
   tan2_atm = Abs(N1L(3,6)/N1L(3,7))**2
   tan2_sol = Abs(N1L(2,5)/N1L(1,5))**2

   If ((tan2_atm.Lt.tan2_atm_min).Or.(tan2_atm.Gt.tan2_atm_max)) Then
    check = .False.

    Lam_sq = Abs(Lambda(2))**2 + Abs(Lambda(3))**2
    Lambda(2) = Sqrt(tan2_atm_opt/(tan2_atm+tan2_atm_opt)) * Lambda(2)
    Lambda(3) = Sqrt(tan2_atm/(tan2_atm+tan2_atm_opt)) * Lambda(3)
    Lam_sq = Lam_sq / (Abs(Lambda(2))**2 + Abs(Lambda(3))**2)
    Lambda(2:3) =  Lambda(2:3) * Sqrt(Lam_sq) 

   End If

   If ((tan2_sol.Lt.tan2_sol_min).Or.(tan2_sol.Gt.tan2_sol_max)) Then
    check = .False.

    Call CalculateEpsTilde_from_Eps(epsT,eps,lambda)
    epsT12 = Abs(epsT(1))**2 + Abs(epsT(2))**2

    epsT(1) = Sqrt(tan2_sol_opt/(tan2_sol+tan2_sol_opt)) * epsT(1)
    epsT(2) = Sqrt(tan2_sol/(tan2_sol+tan2_sol_opt)) * epsT(2)
 
    epsT12 = epsT12 / (Abs(epsT(1))**2 + Abs(epsT(2))**2)
    epsT(1:2) = epsT(1:2) * Sqrt(epsT12)

    Call CalculateEps_from_EpsTilde(eps,epsT,lambda)

    bi(2:4) = eps
   End If

   If (.Not.check) Then
    If (l_phases) Then
     Lam_sq = Dot_product(Lambda,Lambda)
     Call my_vi(mu_mZ, vevSM(1), eps(1), Ue3_opt * Lam_sq, vevL(1), Lambda(1) &
       &      , kont)
     If (kont.Ne.0) Then
      Iname = Iname - 1
      Return
     End If
     ! abuse of notation
     Lam_sq = 0.5_dp * (1._dp-Ue3_opt) * Lam_sq
     Call my_vi(mu_mZ, vevSM(1), eps(2), Lam_sq, vevL(2), Lambda(2), kont)
     If (kont.Ne.0) Then
      Iname = Iname - 1
      Return
     End If
     Call my_vi(mu_mZ, vevSM(1), eps(3), Lam_sq, vevL(3), Lambda(3), kont)
     If (kont.Ne.0) Then
      Iname = Iname - 1
      Return
     End If
    Else ! old procedure
     Lam_sq = Dot_product(Lambda,Lambda)
     vevL = (Lambda - eps*vevSM(1)) / mu_mZ 
    End If

    Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4))
    Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                     &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                     &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                     &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                     &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                     &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                     &, GenerationMixing, kont)
    Cycle
   End If
   !-----------------------------------
   ! and now the reactor constraint
   !-----------------------------------
   Ue32 = Abs(N1l(3,5))**2

   If (Ue32.Gt.Ue32_max) Then
    check = .False.
    Lambda(1) = Lambda(1) / 2._dp

    If (l_phases) Then
     Lam_sq = Dot_product(Lambda,Lambda)
     Call my_vi(mu_mZ, vevSM(1), eps(1), Ue3_opt * Lam_sq, vevL(1), Lambda(1) &
       &      , kont)
     If (kont.Ne.0) Then
      Iname = Iname - 1
      Return
     End If
     ! abuse of notation
     Lam_sq = 0.5_dp * (1._dp-Ue3_opt) * Lam_sq
     Call my_vi(mu_mZ, vevSM(1), eps(2), Lam_sq, vevL(2), Lambda(2), kont)
     If (kont.Ne.0) Then
      Iname = Iname - 1
      Return
     End If
     Call my_vi(mu_mZ, vevSM(1), eps(3), Lam_sq, vevL(3), Lambda(3), kont)
     If (kont.Ne.0) Then
      Iname = Iname - 1
      Return
     End If
     Lam_sq = Dot_product(Lambda,Lambda)
    Else ! old procedure
     Lam_sq = Dot_product(Lambda,Lambda)
     vevL = (Lambda - eps*vevSM(1)) / mu_mZ 
    End If

    Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4))
   End If
   !-----------------------------------
   ! leave loop if everything is fine
   !-----------------------------------
   If (check) Then
    Iname = Iname - 1
    Exit
   End If
   !---------------------------------------------------
   ! else recalculate tree level masses
   !---------------------------------------------------
   Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                    &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                    &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                    &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                    &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                    &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                    &, GenerationMixing, kont)

  End Do

  If(.Not.check)Then 
    Write(ErrCan,*)'Error, no solution found',isol
    kont = -isol
!    Write(*,*)'Error, no solution found',isol
!  Else
!    Write(*,*)'solution found after:',count,Sqrt(Sum(Abs(eps)**2))/Real(mu_mZ,dp) 
  End If

  mN7 = mN1L
  N7 = N1L
  B_i = B_4(2:4)
  Iname = Iname - 1

 Contains

  Subroutine my_vi(mu, vd, eps, lam_sq, v_i, lam_i, kont)
  Implicit None
   Complex(dp), Intent(in) :: mu, eps
   Real(dp), Intent(in) :: vd, lam_sq
   Integer, Intent(inout) :: kont
  
   Real(dp), Intent(out) :: v_i
   Complex(dp), Intent(out) :: lam_i

   ! local variables
   Real(dp) :: alpha, beta2
 
   ! initialisation
   kont = 0
   v_i = 0._dp
   lam_i = 0._dp

   alpha = vd * ( Real(eps,dp)*Real(mu,dp) + Aimag(eps)*Aimag(mu)) / Abs(mu)**2
   beta2 = (lam_sq - vd**2 * Abs(eps)**2) / Abs(mu)**2 - alpha**2

   If (beta2.Gt.0._dp) Then ! there is a solution
    If (alpha.Lt.0._dp) Then ! taking arbitrarily this one
     v_i = -alpha - Sqrt(beta2)
    Else
     v_i = -alpha + Sqrt(beta2)
    End If

   Else ! there is no solution
    kont = -1201
   End If

  End Subroutine my_vi

 End Subroutine Calculate_RP_ParametersC

 Subroutine CalculateEpsTilde_from_Eps(epsT, eps, Lambda)
 Implicit None
  Complex(dp), Intent(in), Dimension(3) :: eps, Lambda 
  Complex(dp), Intent(out), Dimension(3) :: epsT

  Real(dp) :: AbsLam2, Abs23, AbsLam22, Abs232


  Abs232 = Abs(Lambda(2))**2 + Abs(Lambda(3))**2
  AbsLam22 = Abs232 + Abs(Lambda(1))**2

  Abs23 = Sqrt(Abs232)
  AbsLam2 = Sqrt(AbsLam22)

  epsT(1) = ( eps(1) * Abs232                                          &
          & - Conjg(Lambda(1)) * (Lambda(2)*eps(2)+Lambda(3)*eps(3))   &
          & ) / (Abs23*AbsLam2)
  epsT(2) = (Conjg(Lambda(3))*eps(2)-Conjg(Lambda(2))*eps(3)) / Abs23
  epsT(3) = (Lambda(1)*eps(1)+Lambda(2)*eps(2)+Lambda(3)*eps(3)) / AbsLam2

 End Subroutine CalculateEpsTilde_from_Eps

 Subroutine CalculateEps_from_EpsTilde(eps, epsT, Lambda)
 Implicit None
  Complex(dp), Intent(in), Dimension(3) :: epsT, Lambda 
  Complex(dp), Intent(out), Dimension(3) :: eps

  Real(dp) :: AbsLam2, Abs23, AbsLam22, Abs232

  Abs232 = Abs(Lambda(2))**2 + Abs(Lambda(3))**2
  AbsLam22 = Abs232 + Abs(Lambda(1))**2

  Abs23 = Sqrt(Abs232)
  AbsLam2 = Sqrt(AbsLam22)

  eps(1) = ( epsT(3) * Conjg(Lambda(1)) + epsT(1) * Abs23  )  / AbsLam2

  eps(2) = ( Conjg(Lambda(2)) *( epsT(3) * Abs23 -epsT(1) * Lambda(1))   &
         & + epsT(2) * Lambda(3) * AbsLam2 ) / (Abs23 * AbsLam2)

  eps(3) = ( Conjg(Lambda(3)) * ( epsT(3) * Abs23 - epsT(1) * Lambda(1)) &
         & - epsT(2) * Lambda(2) * AbsLam2 ) / (Abs23 * AbsLam2)

 End Subroutine CalculateEps_from_EpsTilde

  
 Subroutine Fit_Neutrino_Data_sp(Nfit, m2_atm_min, m2_atm_max, tan2_atm_min   &
       & , tan2_atm_max, m2_sol_min, m2_sol_max, tan2_sol_min, tan2_sol_max   &
       & , Ue32_max, g, gp, rh0, hpncs, lmbd, vevSM, vR, vS, vP, mPhi, MR, mu &
       & , M1, M2, hnu, vL)
  !--------------------------------------------------------------------------
  ! adjustes h_nu and v_L such that neutrino data are most likely fullfilled
  ! based on the routine FitNeu by Martin Hirsch
  !     options:
  !     Nfit:   Comment:
  !       <0    returns without doing anything, added by Werner
  !       0     No fit at all, Eps and Lam Random numbers, see below
  !       1     (Lam=Atm,Eps=Sol,sign condition = yes),
  !              28% success rate with current NuData
  !
  !       2     (Lam=Atm,Eps=Sol,sign condition = no),
  !              12% success rate with current NuDat
  !
  !       3     (Lam=Sol,Eps=Atm,sign condition = yes),
  !              28% success rate with current NuData
  !
  !       4     (Lam=Sol,Eps=Atm,sign condition = no),
  !               9% success rate with current NuData
  !
  !     Note: Success rate can be enhanced by narrowing the random regions,
  !           but at the expense of not covering uniformly the (currently)
  !           allowed neutrino parameter space
  !--------------------------------------------------------------------------
  Implicit none
   Integer, intent(in) :: Nfit
   Real(dp), intent(in) :: m2_atm_min, m2_atm_max, tan2_atm_min, tan2_atm_max &
       & , m2_sol_min, m2_sol_max, tan2_sol_min, tan2_sol_max, Ue32_max       &
       & , g, gp, rh0, hpncs, lmbd, vevSM(2), vR, vS, vP, mPhi, MR
   Complex(dp), Intent(in) :: M1, M2, mu
   Real(dp), Dimension(3), Intent(out) :: hnu, vL

   Real(dp) :: epssq, lambda, eps(3), lam(3), mptot, mRtot, mphot, x1, x2, x3 &
       & , det7, c1, c2, c3, vu, vd, sgnc
   Real(dp) :: vec8(8), vec13(13) ! vectors of random numbers
   integer :: i1
   !----------------------------------------------
   ! no adjustment if Nfit < 0
   !----------------------------------------------
   if (Nfit.lt.0) return 

   Iname = Iname + 1
   NameOfUnit(Iname) = "Fit_Neutrino_Data_sp"

   !--------------------------------------------
   ! initialization
   !--------------------------------------------
   epssq = 0._dp
   lambda = 0._dp
   eps = 0._dp
   lam = 0._dp

   If (Nfit.gt.0) then
    vd = vevSM(1)
    vu = vevSM(2)
    mptot=mPhi + lmbd * vP * oosqrt2
    mRtot=MR+hpncs*vP * oosqrt2
    mphot=g**2*m1+gp**2*m2
    x1=8.0_dp*m1*m2*mu
    x1=x1*(mptot*mRtot*mu-hpncs**2*mu*vR*vS+rh0**2*mRtot*vd*vu)
    x2=4.0_dp*mu*vd*(mptot*mRtot-hpncs**2*vR*vS)*vu
    x3=rh0**2*mRtot*(vd**2+vu**2)**2
    det7=mRtot*(x1-mphot*(x2+x3))/8.0_dp
  !      det4=mphot/(4.0_dp*m1*m2*mu**2-2.0_dp*mphot*mu*vd*vu)
  !
  !1.2  === c1,c2 and c3 ===
  !
    c1=-hpncs**2*vR*vS*mu+mptot*mRtot*mu
    c1=(c1+rh0**2*mRtot*vd*vu)*mphot*mRtot/4.0_dp/mu
    c2=rh0*mphot*mRtot*(rh0*mRtot+hpncs*mu)*vu*(vu**2-vd**2)/8.0_dp/mu
    c3=(rh0*mRtot+hpncs*mu)**2*vu**2*(2.0_dp*m1*m2*mu-mphot*vd*vu)
    c3=c3/4.0_dp/mu
    c1=-c1/det7
    c2=-c2/det7
    c3=-c3/det7
   end if

   Select Case(Nfit)
   Case(0)
    Call Random_Number(vec8)
    lambda=10.0_dp**(1.0_dp - 4.0_dp*vec8(1))
    epssq =10.0_dp**(0.0_dp - 6.0_dp*vec8(2))
    lam = vec8(3:5)
    lam = lam * lambda / Sqrt(Dot_product(lam,lam)) 
    eps = vec8(6:8)
    eps = eps * Sqrt(epssq/Dot_product(eps,eps)) 

   Case(1)
    Call Random_Number(vec13)
    lambda = Sqrt(sqrt(m2_atm_min+(m2_atm_max-m2_atm_min)*vec13(1))/abs(c1))
    epssq =  abs(sqrt(m2_sol_min+(m2_sol_max-m2_sol_min)*vec13(2))/abs(c3))

    lam(1) = sqrt(Ue32_max * vec13(3))
    lam(2) = Sqrt(tan2_atm_min + (tan2_atm_max - tan2_atm_min) * vec13(4))
    lam(3) = 1._dp
    lam = lam * lambda / Sqrt(Dot_product(lam,lam))

    eps(1) = sqrt(tan2_sol_min + (tan2_sol_max - tan2_sol_min) * vec13(5))
    eps(2) = 1._dp
    eps(3) = sqrt(tan2_sol_min + 0.8_dp * (tan2_sol_max - tan2_sol_min) * vec13(6))
    eps = eps * Sqrt(epssq/Dot_product(eps,eps))

    Do i1=1,3
     If (vec13(6+i1).Lt.0.5_dp) lam(i1) = - lam(i1)
     If (vec13(9+i1).Lt.0.5_dp) eps(i1) = - eps(i1)
    end do

    sgnc=(eps(2)/eps(3))*(lam(2)/lam(3))
    If (sgnc.Gt.0._dp) lam(2) = - lam(2)

   Case(3)
    Call Random_Number(vec13)
    lambda = Sqrt(sqrt(m2_sol_min+(m2_sol_max-m2_sol_min)*vec13(2))/Abs(c1))
    epssq = Abs(sqrt(m2_atm_min+(m2_atm_max-m2_atm_min)*vec13(1))/c3)

    lam(1) = sqrt(tan2_sol_min + (tan2_sol_max - tan2_sol_min) * vec13(5))
    lam(2) = 1._dp
    lam(3) = sqrt(tan2_sol_min + 0.8_dp * (tan2_sol_max - tan2_sol_min) * vec13(6))
    lam = lam * lambda / Sqrt(Dot_product(lam,lam))

    eps(1) = sqrt(Ue32_max * vec13(3))
    eps(2) = sqrt(tan2_atm_min + (tan2_atm_max - tan2_atm_min) * vec13(4))
    eps(3) = 1._dp
    eps = eps * Sqrt(epssq/Dot_product(eps,eps))

    Do i1=1,3
     If (vec13(6+i1).Lt.0.5_dp) lam(i1) = - lam(i1)
     If (vec13(9+i1).Lt.0.5_dp) eps(i1) = - eps(i1)
    end do

    sgnc=(eps(2)/eps(3))*(lam(2)/lam(3))
    If (sgnc.Gt.0._dp) lam(2) = - lam(2)

   Case default
    Write(ErrCan,*) "Problem in routine "//NameOfUnit(Iname)
    Write(ErrCan,*) "Option Nfit=",Nfit," does not exist."
    If (ErrorLevel.Gt.-2) call TerminateProgram
   end select
    vL = (lam - eps*vd) / mu
    hnu = sqrt2 * eps / vR

    Iname = Iname - 1
 
 End Subroutine Fit_Neutrino_Data_sp


 Subroutine Model_bilinear_Rparity(add_Rparity, HighScaleModel, delta, epsI   &
       & , deltaM, ratioWoM, m32, grav_fac, CalcTBD, Ecms, Pm, Pp, ISR, Beam  &
       & , SigSup , SigSdown, SigC, SigChi0, SigS0, SigSP, SigHp, M_GUT, kont)
  Implicit None

  !-----------------------------------------------------------
  ! input / ouput
  !-----------------------------------------------------------
  Logical, Intent(in) :: add_Rparity
  Character(len=15), Intent(inout) :: HighScaleModel
  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: delta   ! required precision for spectrum calculation
  Real(dp), Intent(out) :: M_GUT  ! scale of SUSY boundary conditions,
                                  ! usually the GUT scale
 !-------------------------------
 ! widths and branching ratios
 !-------------------------------
  Real(dp), Intent(in) ::  m32, grav_fac, epsI, deltaM, ratioWoM
  Logical, Intent(in) :: CalcTBD
 !--------------------------------
 ! cross section calculation
 !--------------------------------
  Real(dp), Intent(in) :: Ecms(:), Pm(:), Pp(:)
  Logical, Intent(in) :: ISR(:), Beam(:)
  Real(dp), Intent(out) :: SigSup(:,:,:) , SigSdown(:,:,:), SigC(:,:,:)   &
     & , SigChi0(:,:,:), SigS0(:,:), SigSP(:,:,:), SigHp(:,:,:)

  !-----------------------------------------------------------
  ! local variables
  !-----------------------------------------------------------
  Integer :: i1, i_min(3)
  Real(dp) :: sinW2, gp, g,vev
  Real(dp) :: g0(213), mudim, tz, dt, gauge_mZ(3), M2H_mZ(3), mN1L(7), mN1L2(7)
  Real(dp) :: mGlu_T, mSdown_T(6), mSdown2_T(6), mSup_T(6), mSup2_T(6) &
      & , abs_mu2, mz2_run
  Complex(dp), Dimension(3,3) :: y_l_mZ,  y_d_mZ, y_u_mZ , Al_mZ, Ad_mZ &
           &, Au_mZ,M2E_mZ, M2L_mZ, M2D_mZ, M2Q_mZ, M2U_mZ
  Complex(dp) :: Mi_mZ(3), B_mZ, bi(4), BiEpsi(4), N1L(7,7), PhaseGlu_T &
           & , RSdown_T(6,6), RSup_T(6,6)
  !--------------------------------------
  ! for pMSSM model
  !--------------------------------------
  Real(dp) :: mZ2_t, mW2_t, Scale_Q, g1(57), g58(58), g214(214) &
       & , mN2(4), mC2(2), g2(213), mP0_T(2), mP02_T(2), RP0_T(2,2)
  Real(dp) :: mC_T(2), mN_T(4), mS0_T(2), mSpm_T(2)  &
    & , mSlepton_T(6), mSneut_T(3), mSlepton2_T(6) &
    & , mSneut2_T(3), mS02_T(2), mSpm2_T(2), RS0_T(2,2), mC2_T(2), mN2_T(4)  &
    & , mass_new(32), mass_old(32), diff_m(32)
  Complex(dp) :: U_T(2,2), V_T(2,2), N_T(4,4), RSpm_T(2,2)  &
    & , RSlepton_T(6,6), RSneut_T(3,3), mu_T, B_T, mu_save
  Integer :: i2, i3
  Logical :: Converge, UseFixedScale

  !--------------------------------
  ! cross section calculation
  !--------------------------------
  Integer :: p_max

 !------------------------------------------------------
 ! internal information on particle identities
 !------------------------------------------------------
  Integer :: id_gl, id_ph, id_grav
  Integer, Dimension(1) :: id_Z, id_W
  Integer, Dimension(3) :: id_nu, id_l, id_d, id_u

  Iname = Iname + 1
  NameOfUnit(Iname) = "Model_bilinear_Rparity"

  kont = 0

  Call Initialize_RPexplicit(GenerationMixing, id_gl, id_ph, id_Z, id_W &
               & , id_nu, id_l, id_d, id_u, id_grav)

  If (add_Rparity) Then ! calculate parameters from a high scale model
                       ! assuming conserved R-parity

  !-------------------------------------------------------------------------
  ! Iterative procedure to get the 1-loop masses and mixing matrices;
  ! delta is the relative precision required for the masses after two runs
  ! The procedure is stopped after the i-th iteration 
  ! if delta > |m(i) - m(i-1)|/m(i) for all SUSY. In the case that more than
  ! 20 iterations are necessary the iteration loop is left and a warning is 
  ! given.
  !--------------------------------------------------------------------------
   ! In the SPA convention the the renormalization scale is fixed with 1 TeV
   If (SPA_Convention) Call SetRGEScale(1.e3_dp**2) ! 
   Scale_Q = Sqrt(GetRenormalizationScale())

   If (Scale_Q.Eq.1._dp) Then ! in this case no scale has been set and one
    UseFixedScale = .False.   ! has to use the square root of the stop masses
   Else
    UseFixedScale = .True.
   End If

  If (HighScaleModel.Eq."pMSSM") Then
   sinW2 = 1._dp - mW2/mZ2
   alpha_mZ = Alpha_MSbar(mZ, mW)
   gp = Sqrt( 4._dp*pi*alpha_mZ/(1._dp-sinW2) )
   g = Sqrt( 4._dp*pi*alpha_mZ/sinW2)

   vev =  2._dp * mW / g

   vevSM(1) = vev / Sqrt(1._dp + tanb**2)
   vevSM(2) = tanb * vevSM(1)

   Y_l = 0._dp
   Y_d = 0._dp
   Y_u = 0._dp

   Do i1=1,3
    Y_l(i1,i1) = sqrt2 * mf_l_mZ(i1) / vevSM(1)
    Y_d(i1,i1) = sqrt2 * mf_d_mZ(i1) / vevSM(1)
    Y_u(i1,i1) = sqrt2 * mf_u_mZ(i1) / vevSM(2)
   End Do
   If (GenerationMixing) Then

    i1 = GetYukawaScheme()
    If (i1.Eq.1) Then
     Y_u = Matmul(Transpose(CKM),Y_u) 
    Else
     Y_d = Matmul(Conjg(CKM),Y_d) 
    End If
    Call Switch_from_superCKM(Y_d, Y_u, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm &
              & , M2U_sckm, A_d, A_u, M2_D, M2_Q, M2_U, .False. )
    If (Maxval(Abs(MatNu)).Gt.0._dp) Then
     Call Switch_from_superPMNS(Y_l, MatNu, Al_pmns, M2E_pmns, M2L_pmns &
               & , A_l, M2_E, M2_L, .False. )
    Else
     Call Switch_from_superPMNS(Y_l, id3C, Al_pmns, M2E_pmns, M2L_pmns &
               & , A_l, M2_E, M2_L, .False. )
    End If

   Else
    A_d = Ad_sckm
    A_u = Au_sckm
    M2_D = M2D_sckm
    M2_U = M2U_sckm
    M2_Q = M2Q_sckm
    A_l = Al_pmns
    M2_E = M2E_pmns
    M2_L = M2L_pmns

   End If

   Call TreeMasses(gp, g, vevSM, Mi(1), Mi(2), Mi(3), mu, B                &
        &, tanb, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d               &
        &, A_u, Y_d, Y_u, Glu%m, PhaseGlu, ChiPm%m, ChiPm%m2, U, V, Chi0%m &
        &, Chi0%m2, N, Sneut%m, Sneut%m2, Rsneut, Slepton%m, Slepton%m2    &
        &, RSlepton, Sdown%m, Sdown%m2, RSdown, Sup%m, Sup%m2, RSup        &
        &, mP0_T, mP02_T, RP0, S0%m, S0%m2, RS0, Spm%m, Spm%m2, RSpm       &
        &, GenerationMixing, kont, .False.) ! tree-level Higgs mass

   P0(1)%m = mZ
   P0(1)%m2 = mZ2
   mP0 = P0%m
   mP02 = P0%m2
   
   mSpm = Spm%m    ! short cut to solve a problem with arrays in type variables
   mass_old(1) = Abs(glu%m)
   mass_old(2:3) = Abs(ChiPm%m)
   mass_old(4:7) = Abs(Chi0%m)
   mass_old(8:9) = S0%m
   mass_old(10) = mP0_T(2)
   mass_old(11) = mSpm(2)
   mass_old(12:17) = Sup%m   
   mass_old(18:23) = Sdown%m    
   mass_old(24:29) = Slepton%m  
   mass_old(30:32) = Sneut%m 
   ! calculate first gauge and Yukawa in DR-scheme at m_Z
   Call BoundaryEW(1, vevSM, ChiPm%m, U, V, Chi0%m, N, S0%m2, RS0, mP02, RP0 &
    & , Spm%m, Spm%m2, RSpm, Sdown%m, Sdown%m2, RSdown, Sup%m, Sup%m2, RSup  &
    & , Slepton%m, Slepton%m2, RSlepton, Sneut%m2, RSneut, uU_L, uU_R        &
    & , uD_L, uD_R, uL_L, uL_R, id3C, Glu%m, PhaseGlu, mZ2_t, mW2_t          &
    & , delta, g1, kont)

   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If

   converge = .False.
   mP0(1) = mZ
   mP02(1) = mZ2
   If (.Not.UseFixedScale) Then
    mSup = Sup%m
    If (GenerationMixing) Then
     Scale_Q = 1._dp
     Do i1=1,6
      If ((Abs(Rsup(i1,3))**2+Abs(Rsup(i1,6))**2).Gt.0.6_dp) &
           &  Scale_Q = Scale_Q*mSup(i1)
     End Do
     Scale_Q = Sqrt(Scale_Q)
    Else
     Scale_Q = Sqrt(mSup(5)*mSup(6))
    End If
    tz = SetRenormalizationScale(scale_Q**2)
   End If

   Do i1=1,n_run
    tz = Log(Scale_Q/mZ)
    dt = tz / 50._dp
    kont = 0
    If (tanb_in_at_Q) Then
     Call odeint(g1, 57, 0._dp, tz, 0.1_dp*delta, dt, 0._dp, rge57, kont)
    Else  ! include tanb(beta) in running
     g58(1:57) = g1
     g58(58) = Log(tanb_mZ)
     Call odeint(g58, 58, 0._dp, tz, 0.1_dp*delta, dt, 0._dp, rge58, kont)
     g1 = g58(1:57)
     tanb_Q = Exp(g58(58))
    End If

    !----------------------------------------------
    ! evolve first parameters up to Q
    !----------------------------------------------
    Call GToCouplings(g1, gauge, Y_l, Y_d, Y_u)
    gauge(1) = Sqrt(3._dp/5._dp) * gauge(1)

    If (GenerationMixing) Then
     Y_l = Transpose(Y_l)
     Y_d = Transpose(Y_d)
     Y_u = Transpose(Y_u)
     Call Switch_from_superCKM(Y_d, Y_u, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm &
               & , M2U_sckm, A_d, A_u, M2_D, M2_Q, M2_U, .False. )
     If (Maxval(Abs(MatNu)).Gt.0._dp) Then
      Call Switch_from_superPMNS(Y_l, MatNu, Al_pmns, M2E_pmns, M2L_pmns &
                & , A_l, M2_E, M2_L, .False. )
     Else
      Call Switch_from_superPMNS(Y_l, id3C, Al_pmns, M2E_pmns, M2L_pmns &
                & , A_l, M2_E, M2_L, .False. )
     End If
    Else ! .not. GenerationMixing
     Do i2=1,3
      If (.Not.l_Al) A_l(i2,i2) = AoY_l(i2,i2) * Y_l(i2,i2)
      If (.Not.l_Ad) A_d(i2,i2) = AoY_d(i2,i2) * Y_d(i2,i2)
      If (.Not.l_Au) A_u(i2,i2) = AoY_u(i2,i2) * Y_u(i2,i2)
      Do i3=1,3
       If (i3.Ne.i2) Then
        A_l(i3,i2) = ZeroC
        A_d(i3,i2) = ZeroC
        A_u(i3,i2) = ZeroC
       End If
      End Do
     End Do
    End If

    kont = 0
    Call LoopMassesMSSM_2(delta, tanb_mZ, tanb, gauge, Y_l, Y_d, Y_u, Mi  &
       & , A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, mu                &
       & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP02, RP0           &
       & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2        &
       & , RSup, mSlepton, mSlepton2, RSlepton, mSneut, mSneut2           &
       & , RSneut, mGlu, PhaseGlu, M2_H, B, kont)

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If

    mass_new(1) = Abs(mglu)
    mass_new(2:3) = Abs(mC)
    mass_new(4:7) = Abs(mN)
    mass_new(8:9) = mS0
    mass_new(10) = mP0(2)  
    mass_new(11) = mSpm(2)  
    mass_new(12:17) = mSup   
    mass_new(18:23) = mSdown   
    mass_new(24:29) = mSlepton 
    mass_new(30:32) = mSneut
    diff_m = Abs(mass_new - mass_old)
    Where (mass_old.Ne.0._dp) diff_m = diff_m / mass_old

    If (Maxval(diff_m).Lt.delta) Then
     Glu%m = mglu
     ChiPm%m = mC
     ChiPm%m2 = mC2
     Chi0%m = mN
     Chi0%m2 = mN2
     S0%m = mS0
     S0%m2 = mS02
     P0%m = mP0
     P0%m2 = mP02
     Spm%m = mSpm
     Spm%m2 = mSpm2
     Sup%m = mSup
     Sup%m2 = mSup2
     Sdown%m = mSdown
     Sdown%m2 = mSdown2
     Slepton%m = mSlepton
     Slepton%m2 = mSlepton2
     Sneut%m = mSneut
     Sneut%m2 = mSneut2
     converge = .True.
     Exit
    Else
     mass_old = mass_new
    End If

    If (i1.Eq.101) Then
     Write(*,*) "Problem with accuracy (pMSSM)",diff_m,delta,i1
     Exit
    End If

    gauge(1) = Sqrt(5._dp/3._dp) * gauge(1)

    Y_l = Transpose(Y_l) 
    Y_d = Transpose(Y_d) 
    Y_u = Transpose(Y_u) 
    A_l = Transpose(A_l) 
    A_d = Transpose(A_d) 
    A_u = Transpose(A_u) 
  
    Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u       &
                  & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g2)
    !------------------------------
    ! now the running
    !------------------------------
    If (tanb_in_at_Q) Then ! include tanb(beta) in running
     g214(1:213) = g2
     g214(214) = Log(tanb_Q)
     Call odeint(g214, 214, tz, 0._dp, 0.1_dp*delta, dt, 0._dp, rge214, kont)
     g2 = g214(1:213)

     tanb_mZ = Exp(g214(214))

    Else
     Call odeint(g2, 213, tz, 0._dp, 0.1_dp*delta, dt, 0._dp, rge213, kont)

    End If

    Call GToParameters(g2, gauge_mZ, Y_l_mZ, Y_d_mZ, Y_u_mZ, Mi_mZ, A_l_mZ &
       & , A_d_mZ, A_u_mZ, M2_E_mZ, M2_L_mZ, M2_D_mZ, M2_Q_mZ, M2_U_mZ     &
       & , M2_H_mZ, mu_mZ, B_mZ)
    Y_l = Transpose(Y_l)
    Y_d = Transpose(Y_d)
    Y_u = Transpose(Y_u)
    A_l = Transpose(A_l)
    A_d = Transpose(A_d)
    A_u = Transpose(A_u)

    Y_l_mZ = Transpose(Y_l_mZ) 
    Y_d_mZ = Transpose(Y_d_mZ) 
    Y_u_mZ = Transpose(Y_u_mZ) 
    A_l_mZ = Transpose(A_l_mZ) 
    A_d_mZ = Transpose(A_d_mZ) 
    a_u_mZ = Transpose(A_u_mZ)
    !-----------------------------------------
    ! use consistently running parameters
    !-----------------------------------------
    g2(1) = Sqrt(3._dp/5._dp) * g2(1)
    sinW2 = g2(1)**2 / (g2(1)**2 +g2(2)**2)
    vev =  2._dp * Sqrt(mZ2_t / (g2(1)**2 +g2(2)**2))
    vevSM(1) = vev / Sqrt(1._dp + tanb_mZ**2)
    vevSM(2) = tanb_mZ * vevSM(1)
 
    If (.Not.UseFixedScale) Then
     If (GenerationMixing) Then
      Scale_Q = 1._dp
      Do i2=1,6
       If ((Abs(Rsup(i2,3))**2+Abs(Rsup(i2,6))**2).Gt.0.6_dp) Scale_Q = Scale_Q*mSup(i2)
      End Do
      Scale_Q = Sqrt(Scale_Q)
     Else
      Scale_Q = Sqrt(mSup(5)*mSup(6))
     End If
     tz = SetRenormalizationScale(scale_Q**2)
    End If
    !-----------------------------------
    ! mu and B parameters at tree-level
    !-----------------------------------
    Abs_Mu2 = (M2_H_mZ(2) * tanb_mZ**2 - M2_H(1) ) / (1._dp - tanb_mZ**2) &
          & - 0.5_dp * mZ2_t
    If (Abs_Mu2.Le.0) Then
     kont = -12345
     Return
    End If

    mu_T = Sqrt(Abs_mu2) * phase_mu
    B_T = (M2_H_mZ(1) + M2_H_mZ(2) + 2._dp *  Abs_Mu2) * tanb_mZ / (1+tanb_mZ**2)

    Call TreeMasses(g2(1), g2(2), vevSM, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3), mu_T   &
      & , B_T , tanb_mZ, M2_E_mZ, M2_L_mZ, A_l_mZ, Y_l_mZ, M2_D_mZ, M2_U_mZ   &
      & , M2_Q_mZ, A_d_mZ, A_u_mZ, Y_d_mZ, Y_u_mZ, mGlu_T, PhaseGlu, mC_T     &
      & , mC2_T, U_T, V_T, mN_T, mN2_T, N_T, mSneut_T, mSneut2_T, Rsneut_T    &
      & , mSlepton_T, mSlepton2_T, RSlepton_T, mSdown_T, mSdown2_T, RSdown_T  &
      & , mSup_T, mSup2_T, RSup_T, mP0_T, mP02_T, RP0_T, mS0_T, mS02_T, RS0_T &
      & , mSpm_T, mSpm2_T, RSpm_T, GenerationMixing, kont, .False.)

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
     
    If (Min(Minval(mSup2_T), Minval(mSdown2_T), Minval(mSlepton2_T)   &
       &   , Minval(mSneut2_T), Minval(mS02_T), Minval(mP02_T)       &
       &   , Minval(mSpm2_T)).Gt. 0._dp ) Then
     mC = mC_T
     mC2 = mC2_T
     U = U_T
     V = V_T
     mN = mN_T
     mN2 = mN2_T
     N = N_T
     mSneut = mSneut_T
     mSneut2 = mSneut2_T
     Rsneut = Rsneut_T
     mSlepton = mSlepton_T
     mSlepton2 = mSlepton2_T
     RSlepton = RSlepton_T
     mSDown = mSDown_T
     mSDown2 = mSDown2_T
     RSDown = RSDown_T
     mSup = mSup_T
     mSup2 = mSup2_T
     RSup = RSup_T
     mS0 = mS0_T
     mS02 = mS02_T
     RS0 = RS0_T
     mSpm = mSpm_T
     mSpm2 = mSpm2_T
     RSpm = RSpm_T
     YukScen = 1 ! using running masses for boundary conditions at mZ
    Else
     YukScen = 2 ! using pole masses for boundary conditions at mZ
     mglu_T = mglu
     mP0_T = mP0
     mP02_T = mP02
     RP0_T = RP0     
    End If

    kont = 0
    Call BoundaryEW(i1+1,vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T &
      & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup   &
      & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut, uU_L, uU_R      &
      & , uD_L, uD_R , uL_L, uL_R, id3C, mGlu_T, PhaseGlu, mZ2_t, mW2_t   &
      & , delta, g1, kont)

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do

   If ((kont.Eq.0).And.(.Not.converge)) Then
    Write(ErrCan,*) 'Problem in subroutine CalculateSpectrum, model ' &
                              //Trim(HighScaleModel)//'!!'
    Write(ErrCan,*) "After",n_run,"iterations no convergence found"
    kont = -1200
   End If

   else ! Highscale model

    Call FirstGuess(phase_mu, tanb, Mi, M2_E, M2_L, A_l, M2_D, M2_Q, M2_U, A_d &
           & , A_u, mu, B, M2_H, gp, g, Y_l, Y_d, Y_u, vevSM, mP02, mP0, kont)
    If (kont.Ne.0) then
     Iname = Iname - 1
     Return
    End if
    Call TreeMasses(gp, g, vevSM, Mi(1), Mi(2), Mi(3), mu, B               &
        &, tanb, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d               &
        &, A_u, Y_d, Y_u, Glu%m, PhaseGlu, ChiPm%m, ChiPm%m2, U, V, Chi0%m &
        &, Chi0%m2, N, Sneut%m, Sneut%m2, Rsneut, Slepton%m, Slepton%m2    &
        &, RSlepton, Sdown%m, Sdown%m2, RSdown, Sup%m, Sup%m2, RSup        &
        &, mP0_T, mP02_T, RP0, S0%m, S0%m2, RS0, Spm%m, Spm%m2, RSpm       &
        &, GenerationMixing, kont, .False.) ! tree-level Higgs mass
    If (kont.Ne.0) then
     Iname = Iname - 1
     Return
    End if

    Call Sugra(delta, vevSM, ChiPm%m, U, V, Chi0%m, N, S0%m, S0%m2, RS0   &
     & , P0%m, P0%m2, RP0, Spm%m, Spm%m2, RSpm, Sdown%m, Sdown%m2, RSdown &
     & , Sup%m, Sup%m2, RSup, Slepton%m, Slepton%m2, RSlepton, Sneut%m    &
     & , Sneut%m2, RSneut, Glu%m, PhaseGlu, gauge, uL_L, uL_R, uD_L, uD_R &
     & , uU_L, uU_R, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D   &
     & , M2_Q, M2_U, M2_H, mu, B, m_GUT, kont, WriteOut, n_run)

   end if

   If (kont.Ne.0) then
    Iname = Iname - 1
    Return
   End if

   If (l_fit_RP_parameters) Then

    Call Calculate_RP_Parameters(delta_mass, eps, vevL, Beps, Lam_ex, RSpm8 &
                     & , RS05, RP05, U5, V5, Chi07%m, N7, kont)

   
    If (kont.Ne.0) then
     Iname = Iname - 1
     Return
    End if

   Else ! .not.l_fit_RP_parameters
   !-------------------------------------------------
   ! evolve the parameters down to m_Z if necessary 
   !-------------------------------------------------
    mudim = GetRenormalizationScale()
    Call ParametersToG(gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
                    &,M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g0)

    If (mudim.Ne.mZ2) Then
 
     tz = 0.5_dp * Log(mZ2/mudim)
     dt = tz / 100._dp
     g0(1) = Sqrt(5._dp / 3._dp ) * g0(1)
 
     Call odeint(g0, 213, 0._dp, tz, delta, dt, 0._dp, rge213, kont)
     g0(1) = Sqrt(3._dp / 5._dp ) * g0(1)
 
    End If
  
    Call GToParameters(g0, gauge_mZ, y_l_mZ,  y_d_mZ, y_u_mZ, Mi_mZ, Al_mZ  &
          & , Ad_mZ, Au_mZ,M2E_mZ, M2L_mZ, M2D_mZ, M2Q_mZ, M2U_mZ, M2H_mZ   &
          & , mu_mZ, B_mZ)

    mZ2_run = (gp**2+g**2)*0.25_dp*(vevSM(1)**2+vevSM(2)**2)
    abs_mu2 = (M2H_mZ(2) * tanb_mz**2 - M2H_mZ(1) ) / (1._dp-tanb_mZ**2)  &
          & - 0.5_dp * mZ2_run
    mu_mZ = phase_mu * Sqrt(abs_mu2)

    B_mZ = (M2H_mZ(1) + M2H_mZ(2) + 2._dp * Abs_Mu2) * tanb_mz / (1+tanb_mz**2)

    Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, Beps)
    BiEpsi(1) = B_mZ
    BiEpsi(2:4) = Beps
    bi(1) = mu_mZ
    bi(2:4) = eps
 
    Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2)        &
           & , Mi_mZ(3), bi, BiEpsi, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ           &
           & , M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ          &
           & , mGlu_T, PhaseGlu_T, ChiPm5%m, ChiPm5%m2, U5, V5, Chi07%m      &
           & , Chi07%m2, N7, mSdown_T, mSdown2_T, RSdown_T, mSup_T, mSup2_T  &
           & , RSup_T, P05%m, P05%m2, RP05, S05%m, S05%m2, RS05              &
           & , Spm8%m, Spm8%m2, RSpm8, GenerationMixing, kont)

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
  
    Call NeutralinoMass_Loop_RP(g0(1), g0(2), Y_d_mZ, Y_l_mZ, Y_u_mZ, vevSM    &
         & , vevL, Mi_mZ(1), Mi_mZ(2), mu_mZ, eps, ChiPm5%m, ChiPm5%m2, U5, V5 &
         & , mSup2_T, RSup_T, mSdown2_T, RSdown_T, S05%m2, RS05, P05%m2, RP05  &
         & , Spm8%m2, RSpm8, uD_L, uD_R, uU_L, uU_R, Chi07%m, Chi07%m2, N7     &
         & , mN1L, mN1L2, N1L, kont, .False.)

    Chi07%m = mN1L
    Chi07%m2 = mN1L2
    N7 = N1L

   End If  ! l_fit_RP_parameters

   !-------------------------------------------------------------------
   ! replacing the tree level masses by the 1-loop masses of the MSSM
   ! assuming that the effect of the RP parameter is tiny
   ! this is justified in the region of parameter space where neutrino
   ! data are correctely explained
   !-------------------------------------------------------------------
    ChiPm5(1:3)%m = mf_l
    ChiPm5(4:5)%m = ChiPm%m
    Chi07(4:7)%m = Chi0%m

    P05(1)%m = mZ
    Do i1=2,5
     If ( (RP05(i1,1)**2 +RP05(i1,2)**2).Gt.0.5_dp) P05(i1)%m = P0(2)%m
     If ( RP05(i1,3)**2.Gt.0.5_dp) P05(i1)%m = Sneut(1)%m
     If ( RP05(i1,4)**2.Gt.0.5_dp) P05(i1)%m = Sneut(2)%m
     If ( RP05(i1,5)**2.Gt.0.5_dp) P05(i1)%m = Sneut(3)%m
    End Do
    P05%m2 = P05%m**2

    i_min = 0
    Do i1=1,5
     If ( (RS05(i1,1)**2 +RS05(i1,2)**2).Gt.0.5_dp) Then
      S05(i1)%m = S0(1 + i_min(1) )%m
      i_min(1) = i_min(1) + 1
     End If
     If ( RS05(i1,3)**2.Gt.0.5_dp) S05(i1)%m = Sneut(1)%m
     If ( RS05(i1,4)**2.Gt.0.5_dp) S05(i1)%m = Sneut(2)%m
     If ( RS05(i1,5)**2.Gt.0.5_dp) S05(i1)%m = Sneut(3)%m
    End Do
    S05%m2 = S05%m**2

    Spm8(1)%m = mW
    i_min = 0
    Do i1=2,8
     If ( (Abs(RSpm8(i1,1))**2 +Abs(RSpm8(i1,2))**2).Gt.0.5_dp) &
        &  Spm8(i1)%m = Spm(2)%m
     If ( (Abs(RSpm8(i1,3))**2 +Abs(RSpm8(i1,6))**2).Gt.0.5_dp) Then
      Spm8(i1)%m = Slepton(1 + i_min(1) )%m
      i_min(1) = i_min(1) + 1
     End If
     If ( (Abs(RSpm8(i1,4))**2 +Abs(RSpm8(i1,7))**2).Gt.0.5_dp) Then
      Spm8(i1)%m = Slepton(3 + i_min(2) )%m
      i_min(2) = i_min(2) + 1
     End If
     If ( (Abs(RSpm8(i1,5))**2 +Abs(RSpm8(i1,8))**2).Gt.0.5_dp) Then
      Spm8(i1)%m = Slepton(5 + i_min(3) )%m
      i_min(3) = i_min(3) + 1
     End If
    End Do
    Spm8%m2 = Spm8%m**2

   HighScaleModel = "RPexplicit"

  Else ! .not.add_Rparity
   sinW2 = 1._dp - mW2/mZ2

   alpha_mZ = Alpha_MSbar(mZ, mW)
   gauge(1) = Sqrt( 4._dp*pi*alpha_mZ/(1._dp-sinW2) )
   gauge(2) = Sqrt( 4._dp*pi*alpha_mZ/sinW2)
   gauge(3) = Sqrt( 4._dp*pi*alphas_mZ)
   gp = gauge(1)
   g = gauge(2)

   vev =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
   vevSM(1) = vev / Sqrt(1._dp + tanb**2)
   vevSM(2) = tanb * vevSM(1)

   Y_l = 0._dp
   Y_d = 0._dp
   Y_u = 0._dp
   A_l = 0._dp
   A_d = 0._dp
   A_u = 0._dp

   Do i1=1,3
    Y_l(i1,i1) = sqrt2 * mf_l_mZ(i1) / vevSM(1)
    A_l(i1,i1) = AoY_l(i1,i1) *  Y_l(i1,i1)
    Y_d(i1,i1) = sqrt2 * mf_d_mZ(i1) / vevSM(1)
    A_d(i1,i1) = AoY_d(i1,i1) *  Y_d(i1,i1)
    Y_u(i1,i1) = sqrt2 * mf_u_mZ(i1) / vevSM(2)
    A_u(i1,i1) = AoY_u(i1,i1) *  Y_u(i1,i1)
   End Do
   If (GenerationMixing) Then
    i1 = GetYukawaScheme()
    If (i1.Eq.1) Then
     Y_u = Matmul(Transpose(CKM),Y_u) 
    Else
     Y_d = Matmul(Conjg(CKM),Y_d) 
    End If
   End If
   Call Calculate_Bi(mu, eps, vevL, vevSM, gp, g, M2L_pmns, BiEpsi(2:4))
   bi(1) = mu
   bi(2:4) = eps
   BiEpsi(1) = B * mu
!   BiEpsi(2:4) = Beps * eps 
   Beps = BiEpsi(2:4) / eps

   Call TreeMassesEps3(gp, g, vevSM, vevL, Mi(1), Mi(2), Mi(3), bi, BiEpsi     &
      & , M2E_pmns, M2L_pmns, A_l, Y_l, M2D_sckm, M2U_sckm, M2Q_sckm, A_d, A_u &
      & , Y_d, Y_u, Glu%m, PhaseGlu, ChiPm5%m, ChiPm5%m2, U5, V5, Chi07%m      &
      & , Chi07%m2, N7, Sdown%m, Sdown%m2, RSdown, Sup%m, Sup%m2, RSup, P05%m  &
      & , P05%m2, RP05, S05%m, S05%m2, RS05, Spm8%m, Spm8%m2, RSpm8            &
      & , GenerationMixing, kont)

   If (kont.Ne.0) then
    Iname = Iname - 1
    Return
   End if

  End If ! add_Rparity

   !----------------------------------------------
   ! reorder state identification if necessary
   !----------------------------------------------
   If (.Not.GenerationMixing) Then
    Call Swap_Order_Sf(RSdown(1,1), Sdown(1)%id, Sdown(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSdown(3,3), Sdown(3)%id, Sdown(4)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSup(1,1), Sup(1)%id, Sup(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSup(3,3), Sup(3)%id, Sup(4)%id, id_p, c_name)    
   End If
  !-----------------------------------------------------------------------
  ! as there is a large hierachy between m_nu and m_chi, there might have
  ! been numerical problems. However, in all known cases this has been
  ! the mass of the lightest neutrino, which is of no importance for the
  ! following. Therefore, the error flag is reset to zero. 
  ! This is of course dangerous and has to be improved.
  !-----------------------------------------------------------------------
  If (kont.Eq.-14) kont = 0

  
  If ((L_BR).And.(kont.Eq.0)) Then

   Call CalculateBR(0, id_nu, 0, id_l, 3, id_d, 3, id_u, 1, id_Z, 1, id_W, 0   &
    & , 0, 6, 6, 7, 5, 1, 5, 5, 8, id_grav, id_gl, id_ph, gauge, Glu, PhaseGlu &
    & , ChiPm5, U5, V5, Chi07, N7, Sup, RSup, Sdown, RSdown, uD_L, uD_R, uU_L  &
    & , uU_R, S05, RS05, P05, RP05, Spm8, RSpm8, epsI, deltaM, CalcTBD         &
    & , ratioWoM, Y_d, A_d, Y_l, A_l, Y_u, A_u, mu, eps, vevSM, vevL, F_Gmsb   &
    & , m32, grav_fac)

  End If


 !---------------------------------------------------------------------------
 ! Calculation of the cross sections in e+ e- annihilation provided L_Cs is
 ! set .TRUE. (default) and that the routine Sugra has finished
 ! correctly (kont.eq.0) 
 ! The following input quantities can be specified in the file 
 ! CrossSections.in: Ecms .... c.m.s enerergy in GeV
 !                   Pm ...... degree of longitudinal polarization of incoming
 !                             electron
 !                   Pp ...... degree of longitudinal polarization of incoming
 !                             positron
 !                   ISR ..... if .TRUE. then the effect of initial state
 !                             radiation will be included
 ! In the case that the file CrossSections.in does not exist, the following
 ! default values are used: Ecms = 500 GeV, Pm = Pp = 0, ISR = .TRUE.
 !----------------------------------------------------------------------------
! only partially usable and useful
  If ((L_CSrp).And.(kont.Eq.0)) Then
   p_max = Size(Pm)
   SigSup = 0._dp
   SigSdown = 0._dp
   SigC = 0._dp
   SigChi0  = 0._dp
   SigS0 = 0._dp
   SigSP = 0._dp
   SigHp = 0._dp
   Do i1=1,p_max
    If (Ecms(i1).Eq.0._dp) Exit
    Call CalculateCrossSections(Ecms(i1), Pm(i1), Pp(i1), ISR(i1), Beam(i1)    &
           & , "Tesla800", Sup%m, RSup, mf_u, Sdown%m, RSdown, mf_d, Glu%m     &
           & , SigSup(i1,:,:), SigSdown(i1,:,:), ChiPm5%m, U5, V5, Chi07%m, N7 &
           & , SigC(i1,:,:), SigChi0(i1,:,:), S05%m, RS05, vevSM, vevL         &
           & , P05%m, RP05, Spm8%m, RSpm8, SigS0(i1,:), SigSP(i1,:,:)          &
           & , SigHp(i1,:,:) )
   End Do

  End If

  lam_ex = mu_mZ * vevL + vevSM(1) * eps
  mf_nu = Chi07(1:3)%m

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


 End Subroutine Model_bilinear_Rparity

 Subroutine Model_trilinearL_Rparity(delta, epsI, deltaM, ratioWoM    &
       & , m32, grav_fac, CalcTBD, Ecms, Pm, Pp, ISR, Beam            &
       & , SigSup , SigSdown, SigC, SigChi0, SigS0, SigSP, SigHp, kont)
  Implicit None

  !-----------------------------------------------------------
  ! input / ouput
  !-----------------------------------------------------------
  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: delta   ! required precision for spectrum calculation
 !-------------------------------
 ! widths and branching ratios
 !-------------------------------
  Real(dp), Intent(in) ::  m32, grav_fac, epsI, deltaM, ratioWoM
  Logical, Intent(in) :: CalcTBD
 !--------------------------------
 ! cross section calculation
 !--------------------------------
  Real(dp), Intent(in) :: Ecms(:), Pm(:), Pp(:)
  Logical, Intent(in) :: ISR(:), Beam(:)
  Real(dp), Intent(out) :: SigSup(:,:,:) , SigSdown(:,:,:), SigC(:,:,:)   &
     & , SigChi0(:,:,:), SigS0(:,:), SigSP(:,:,:), SigHp(:,:,:)

  !-----------------------------------------------------------
  ! local variables
  !-----------------------------------------------------------
  Integer :: i1
  Real(dp) :: sinW2, gp, g,vev
  Real(dp) :: mN1L(7), mN1L2(7)
  Real(dp) :: mGlu_T, mSdown_T(6), mSdown2_T(6), mSup_T(6), mSup2_T(6)
  Complex(dp) :: bi(4), BiEpsi(4), N1L(7,7), RSdown_T(6,6), RSup_T(6,6)

 !------------------------------------------------------
 ! internal information on particle identities
 !------------------------------------------------------
  Integer :: id_gl, id_ph, id_grav
  Integer, Dimension(1) :: id_Z, id_W
  Integer, Dimension(3) :: id_nu, id_l, id_d, id_u

  Iname = Iname + 1
  NameOfUnit(Iname) = "Model_trilinearL_Rparity"

  kont = 0

  Call Initialize_RPexplicit(GenerationMixing, id_gl, id_ph, id_Z, id_W &
               & , id_nu, id_l, id_d, id_u, id_grav)

  sinW2 = 1._dp - mW2/mZ2

  alpha_mZ = Alpha_MSbar(mZ, mW)
  gauge(1) = Sqrt( 4._dp*pi*alpha_mZ/(1._dp-sinW2) )
  gauge(2) = Sqrt( 4._dp*pi*alpha_mZ/sinW2)
  gauge(3) = Sqrt( 4._dp*pi*alphas_mZ)
  gp = gauge(1)
  g = gauge(2)

  vev =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
  vevSM(1) = vev / Sqrt(1._dp + tanb**2)
  vevSM(2) = tanb * vevSM(1)

  Y_l = 0._dp
  Y_d = 0._dp
  Y_u = 0._dp

  Do i1=1,3
   Y_l(i1,i1) = sqrt2 * mf_l_mZ(i1) / vevSM(1)
   Y_d(i1,i1) = sqrt2 * mf_d_mZ(i1) / vevSM(1)
   Y_u(i1,i1) = sqrt2 * mf_u_mZ(i1) / vevSM(2)
  End Do
  If (GenerationMixing) Then
   i1 = GetYukawaScheme()
   If (i1.Eq.1) Then
    Y_u = Matmul(Transpose(CKM),Y_u) 
   Else
    Y_d = Matmul(Conjg(CKM),Y_d) 
   End If
  End If
  !-------------------------------------
  ! needs to be changed
  !-------------------------------------
   M2_E = M2E_pmns
   M2_L = M2L_pmns
   A_L = Al_pmns
   M2_D = M2D_sckm
   M2_Q = M2Q_sckm
   M2_U = M2U_sckm
   A_d = Ad_sckm
   A_u = Au_sckm
  !-------------------------------------
  ! the above needs to be changed
  !-------------------------------------
  Call Calculate_Bi(mu, eps, vevL, vevSM, gp, g, M2_L, Beps)
  bi(1) = mu
  bi(2:4) = eps
  BiEpsi(1) = B
  BiEpsi(2:4) = Beps
  !---------------------------------------------------------
  ! for completeness also the softs for the Higgs doublets
  !---------------------------------------------------------
  M2_H = - Abs(mu)**2 - 0.5_dp * mZ2 * (1._dp-tanb**2)/(1._dp+tanb**2)
  M2_H(1) = M2_H(1) + Real(B,dp)*tanb
  M2_H(2) = M2_H(2) + Real(B,dp)/tanb
 
  Call TreeMassesLam3(gp, g, vevSM, vevL, Mi(1), Mi(2), Mi(3), bi, BiEpsi     &
      & , M2_E, M2_L, A_l, Y_l, Rp_lam, RP_Alam, M2_D, M2_U, M2_Q, A_d, A_u   &
      & , Y_d, Y_u, Rp_lamp, RP_Alamp                                         &
      & , Glu%m, PhaseGlu, ChiPm5%m, ChiPm5%m2, U5, V5, Chi07%m, Chi07%m2, N7 &
      & , Sdown%m, Sdown%m2, RSdown, Sup%m, Sup%m2, RSup, P05%m, P05%m2, RP05 &
      & , S05%m, S05%m2, RS05, Spm8%m, Spm8%m2, RSpm8, GenerationMixing, kont)

  If (kont.Ne.0) then
   Iname = Iname - 1
   Return
  End if

  Call SetRGEscale(mZ2)  
  Call NeutralinoMass_Loop_RPtri(gp, g, Y_d, Y_l, Y_u, RP_lam, RP_lamp        &
         & , vevSM, vevL, Mi(1), Mi(2), mu, eps, ChiPm5%m, ChiPm5%m2, U5, V5  &
         & , Sup%m2, RSup, Sdown%m2, RSdown, S05%m2, RS05, P05%m2, RP05 &
         & , Spm8%m2, RSpm8, uD_L, uD_R, uU_L, uU_R, Chi07%m, Chi07%m2, N7    &
         & , mN1L, mN1L2, N1L, kont, .False.)

  If (kont.Ne.0) then
   Iname = Iname - 1
   Return
  End if

   Chi07%m = mN1L
   Chi07%m2 = mN1L2
   N7 = N1L
 
   !----------------------------------------------
   ! reorder state identification if necessary
   !----------------------------------------------
   If (.Not.GenerationMixing) Then
    Call Swap_Order_Sf(RSdown(1,1), Sdown(1)%id, Sdown(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSdown(3,3), Sdown(3)%id, Sdown(4)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSup(1,1), Sup(1)%id, Sup(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSup(3,3), Sup(3)%id, Sup(4)%id, id_p, c_name)    
   End If
  !-----------------------------------------------------------------------
  ! as there is a large hierachy between m_nu and m_chi, there might have
  ! been numerical problems. However, in all known cases this has been
  ! the mass of the lightest neutrino, which is of no importance for the
  ! following. Therefore, the error flag is reset to zero. 
  ! This is of course dangerous and has to be improved.
  !-----------------------------------------------------------------------
  If (kont.Eq.-14) kont = 0

  
  If ((L_BR).And.(kont.Eq.0)) Then

   Call CalculateBR(0, id_nu, 0, id_l, 3, id_d, 3, id_u, 1, id_Z, 1, id_W, 0   &
    & , 0, 6, 6, 7, 5, 1, 5, 5, 8, id_grav, id_gl, id_ph, gauge, Glu, PhaseGlu &
    & , ChiPm5, U5, V5, Chi07, N7, Sup, RSup, Sdown, RSdown, uD_L, uD_R, uU_L  &
    & , uU_R, S05, RS05, P05, RP05, Spm8, RSpm8, epsI, deltaM, CalcTBD         &
    & , ratioWoM, Y_d, A_d, Y_l, A_l, Y_u, A_u, mu, eps, RP_lam, RP_lamp       &
    & , vevSM, vevL, F_Gmsb, m32, grav_fac)

  End If


 !---------------------------------------------------------------------------
 ! Calculation of the cross sections in e+ e- annihilation provided L_Cs is
 ! set .TRUE. (default) and that the routine Sugra has finished
 ! correctly (kont.eq.0) 
 ! The following input quantities can be specified in the file 
 ! CrossSections.in: Ecms .... c.m.s enerergy in GeV
 !                   Pm ...... degree of longitudinal polarization of incoming
 !                             electron
 !                   Pp ...... degree of longitudinal polarization of incoming
 !                             positron
 !                   ISR ..... if .TRUE. then the effect of initial state
 !                             radiation will be included
 ! In the case that the file CrossSections.in does not exist, the following
 ! default values are used: Ecms = 500 GeV, Pm = Pp = 0, ISR = .TRUE.
 !----------------------------------------------------------------------------
! disabled for the moment being,
  If ((L_CS).And.(kont.Eq.0)) Then
  SigSup = 0._dp
  SigSdown = 0._dp
  SigC = 0._dp
  SigChi0  = 0._dp
  SigS0 = 0._dp
  SigSP = 0._dp
  SigHp = 0._dp
  Write(ErrCan,*) "Warning, no e+ e- cross sections in case of"
  Write(ErrCan,*) "trilinear RP violation!"
  End If

  lam_ex = mu_mZ * vevL + vevSM(1) * eps
  mf_nu = Chi07(1:3)%m

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


 End Subroutine Model_trilinearL_Rparity

 Subroutine ParMin_sp(g, gp, hn, rh0, h, lmbd, Ahn, Ah, RAh0, Almbd, vevSM, vL &
      & , vR, vS, vP, mu, MPhi, MR, B, bmr, BMphi, mH2, mL2, mR2, mS2, mP2)
 !------------------------------------------------------------------------
 ! this routine solves the tad-pol equations in terms of the soft susy masses
 ! squared in the spontaneous model.
 ! taken from J.Romao and adjusted for SPheno
 !------------------------------------------------------------------------
 Implicit None

   Real(dp), Intent(in) :: g, gp, rh0, h, lmbd, Ah, RAh0, Almbd, vevSM(2) &
      & , vL(3) , vR, vS, vP, MPhi, MR, bmr, BMphi, hn(3), Ahn(3)
   Complex(dp), Intent(in) :: B, mu
   Real(dp), Intent(out) :: mH2(2), mL2(3), mR2, mS2, mP2

  Real(dp) :: g2, gp2, vd, vu, t1, t2, t3

  Iname = Iname + 1
  NameOfUnit(Iname) = "ParMin_sp"

  g2=g*g
  gp2=gp*gp
  vd = vevSM(1)
  vu = vevSM(2)
  
  mH2(1)=(B*vu+rh0*vu*( (lmbd*vp**2)/4._dp+(h*vR*vS)/2._dp        &
    & -(rh0*vd*vu)/2._dp+(Mphi*vp)/sqrt2)+(RAh0*vp*vu)/sqrt2           &
    & +(gp2*vd*(-vd**2/2._dp+vu**2/2._dp-vL(1)**2/2._dp-vL(2)**2/2._dp  &
    & -vL(3)**2/2._dp))/4._dp-(- mu -(rh0*vp)/sqrt2)*sqrt2*(-(rh0*vd   &
    & *vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1)*vL(1))/sqrt2         &
    & +(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3))/sqrt2))/sqrt2)-(g2     &
    & *(2*vd*sqrt2*(vd**2/2._dp+vL(1)**2/2._dp+vL(2)**2/2._dp+vL(3)**2 &
    & /2._dp)-vd*sqrt2*(vd**2/2._dp+vu**2/2._dp+vL(1)**2/2._dp+vL(2)**2 &
    & /2._dp+vL(3)**2/2._dp)))/(4._dp*sqrt2))/vd

  t1=1/vu
  t2=B*vd-(vR**2*vu*hn(1)**2)/2._dp-(vR**2*vu*hn(2)**2)/2._dp&
    &-(vR**2*vu*hn(3)**2)/2._dp+rh0*vd*( (lmbd*vp**2)/4._dp+(h&
    &*vR*vS)/2._dp-(rh0*vd*vu)/2._dp+(Mphi*vp)/sqrt2)
  t3=(RAh0*vd*vp)/sqrt2-(- mu -(rh0*vp)/sqrt2)*(-(rh0*vp*vu) &
    /2._dp-( mu *vu)/sqrt2)*sqrt2-vR*((Ahn(1)*vL(1))/sqrt2&
    &+(Ahn(2)*vL(2))/sqrt2+(Ahn(3)*vL(3))/sqrt2)-sqrt2*((h&
    &*vp*vS)/2._dp+(MR*vS)/sqrt2+(vu*hn(1)*vL(1))/2._dp+(vu&
    &*hn(2)*vL(2))/2._dp+(vu*hn(3)*vL(3))/2._dp)*((hn(1)*vL(1)) &
    /sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3))/sqrt2)-(gp2&
    &*vu*(-vd**2/2._dp+vu**2/2._dp-vL(1)**2/2._dp-vL(2)**2/2._dp&
    &-vL(3)**2/2._dp))/4._dp-(g2*(vu**3*sqrt2-vu*sqrt2*(vd**2/2._dp&
    &+vu**2/2._dp+vL(1)**2/2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)))/(4._dp*sqrt2)
  mH2(2)=t1*(t2+t3)
  mL2(1)=(-((vR*vu*Ahn(1))/sqrt2)-vu*hn(1)*((h*vp&
    &*vS)/2._dp+(MR*vS)/sqrt2+(vu*hn(1)*vL(1))/2._dp+(vu*hn(2) &
    *vL(2))/2._dp+(vu*hn(3)*vL(3))/2._dp)+(gp2*vL(1)*(-vd**2/2._dp&
    &+vu**2/2._dp-vL(1)**2/2._dp-vL(2)**2/2._dp-vL(3)**2/2._dp))/4._dp-vR&
    &*hn(1)*(-(rh0*vd*vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1) &
    *vL(1))/sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3)) &
    /sqrt2))/sqrt2)-(g2*(2*sqrt2*vL(1)*(vd**2/2._dp+vL(1)**2 &
    /2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)-sqrt2*vL(1)*(vd**2/2._dp&
    &+vu**2/2._dp+vL(1)**2/2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)))/(4._dp*sqrt2))/vL(1)
  mL2(2)=(-((vR*vu*Ahn(2))/sqrt2)-vu*hn(2)*((h*vp&
    &*vS)/2._dp+(MR*vS)/sqrt2+(vu*hn(1)*vL(1))/2._dp+(vu*hn(2) &
    *vL(2))/2._dp+(vu*hn(3)*vL(3))/2._dp)+(gp2*vL(2)*(-vd**2/2._dp&
    &+vu**2/2._dp-vL(1)**2/2._dp-vL(2)**2/2._dp-vL(3)**2/2._dp))/4._dp-vR&
    &*hn(2)*(-(rh0*vd*vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1) &
    *vL(1))/sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3)) &
    /sqrt2))/sqrt2)-(g2*(2*sqrt2*vL(2)*(vd**2/2._dp+vL(1)**2 &
    /2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)-sqrt2*vL(2)*(vd**2/2._dp&
    &+vu**2/2._dp+vL(1)**2/2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)))/(4._dp*sqrt2))/vL(2)
  mL2(3)=(-((vR*vu*Ahn(3))/sqrt2)-vu*hn(3)*((h*vp&
    &*vS)/2._dp+(MR*vS)/sqrt2+(vu*hn(1)*vL(1))/2._dp+(vu*hn(2) &
    *vL(2))/2._dp+(vu*hn(3)*vL(3))/2._dp)+(gp2*vL(3)*(-vd**2/2._dp&
    &+vu**2/2._dp-vL(1)**2/2._dp-vL(2)**2/2._dp-vL(3)**2/2._dp))/4._dp-vR&
    &*hn(3)*(-(rh0*vd*vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1) &
    *vL(1))/sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3)) &
    /sqrt2))/sqrt2)-(g2*(2*sqrt2*vL(3)*(vd**2/2._dp+vL(1)**2 &
    /2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)-sqrt2*vL(3)*(vd**2/2._dp&
    &+vu**2/2._dp+vL(1)**2/2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)))/(4._dp*sqrt2))/vL(3)

  mP2=(-(BMphi*vp)-h*vR*((h*vp*vR)/2._dp+(MR*vR) &
    /sqrt2)+rh0*vu*(-(rh0*vp*vu)/2._dp-( mu *vu)/sqrt2)-(Almbd&
    &*vp**2)/(2._dp*sqrt2)-(Ah*vR*vS)/sqrt2+(RAh0*vd*vu)/sqrt2&
    &-(Mphi+(lmbd*vp)/sqrt2)*( (lmbd*vp**2) &
    /4._dp+(h*vR*vS)/2._dp-(rh0*vd*vu)/2._dp+(Mphi*vp)/sqrt2)*sqrt2&
    &-h*vS*((h*vp*vS)/2._dp+(MR*vS)/sqrt2+(vu*hn(1)*vL(1))/2._dp&
    &+(vu*hn(2)*vL(2))/2._dp+(vu*hn(3)*vL(3))/2._dp)+rh0*vd*(-(rh0&
    &*vd*vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1)*vL(1))/sqrt2&
    &+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3))/sqrt2))/sqrt2))/vp
  mR2=(-(BMR*vR)-h*vR*( (lmbd*vp**2)/4._dp+(h*vR&
    &*vS)/2._dp-(rh0*vd*vu)/2._dp+(Mphi*vp)/sqrt2)-(Ah*vp*vR)/sqrt2&
    &-(MR+(h*vp)/sqrt2)*sqrt2*((h*vp*vS)/2._dp+(MR*vS)/sqrt2&
    &+(vu*hn(1)*vL(1))/2._dp+(vu*hn(2)*vL(2))/2._dp+(vu*hn(3)&
    & *vL(3))/2._dp))/vS
  mS2=(-(BMR*vS)-(vR*vu**2*hn(1)**2)/2._dp-(vR*vu**2&
    &*hn(2)**2)/2._dp-(vR*vu**2*hn(3)**2)/2._dp-h*vS*( (lmbd&
    &*vp**2)/4._dp+(h*vR*vS)/2._dp-(rh0*vd*vu)/2._dp+(Mphi*vp)/sqrt2) &
    -(Ah*vp*vS)/sqrt2-(MR+(h*vp)/sqrt2)*((h*vp*vR)/2._dp+(MR*vR) &
    /sqrt2)*sqrt2-(vu*Ahn(1)*vL(1))/sqrt2-(vu*Ahn(2) &
    *vL(2))/sqrt2-(vu*Ahn(3)*vL(3))/sqrt2-sqrt2*((hn(1) &
    *vL(1))/sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3)) &
    /sqrt2)*(-(rh0*vd*vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1) &
    *vL(1))/sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3)) /sqrt2))/sqrt2))/vR

  Iname = Iname - 1

 End Subroutine ParMin_sp


End Module RPtools
