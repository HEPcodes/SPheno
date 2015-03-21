!------------------------------------------------------------------------
! This is the main program for the package SPheno, version 3.0beta
! SPheno stands for S(upersymmetric) Pheno(menology) 
! A documentation is provided in the file Spheno.ps
! Please send any comments and/or bugs to porod@physik.uni-wuerzburg.de
! In case of a bug:
!   - please check that the integer in the first line of Control.in
!     is 0 or larger. If it is negative, please set it to 0. This
!     enables the writing of several warnings. Afterwards rerun the 
!     program.
!   - please send the file Messages.out  and all input files (*.in) to
!     the email above, if possible with a copy of the screen output 
! written by Werner Porod
!-----------------------------------------------------------------
Program SPheno

!---------------------------------
! loading necessary modules
!---------------------------------
Use Control
Use Mathematics
Use LoopFunctions
Use RGEs
Use StandardModel
Use LHC_observables
Use EplusEminusProduction
Use Model_Data
Use MSSM_Data
Use BranchingRatios
Use LoopMasses
Use LowEnergy
Use SugraRuns
Use NMSSM_tools 
Use RPtools
Use InputOutput

! Use f90_unix_proc
 Implicit None
 !-------------------------------
 ! widths and branching ratios
 !-------------------------------
  Real(dp) ::  epsI, deltaM
 !--------------------------------
 ! cross section calculation
 !--------------------------------
 Integer, Parameter :: p_max=100
 Complex(dp) :: Ylp(3,3)
 Real(dp) :: Ecms(p_max), Pm(p_max), Pp(p_max), SigSup(p_max,6,6)        &
         & , SigSdown(p_max,6,6), SigSle(p_max,6,6), SigSn(p_max,3,3)    &
         & , SigC(p_max,5,5), SigChi0(p_max,7,7), SigS0(p_max,5)         &
         & , SigSP(p_max,5,4), SigHp(p_max,7,7)
 Logical :: ISR(p_max)=.False., Beam(p_max)=.False.
 !----------------------------------
 ! low energy constraints
 !----------------------------------
 Real(dp) :: BRbtosgamma, Bs_ll(3), BrBToSLL, BR_Bu_TauNu, a_e, a_mu, a_tau    &
   & , d_e, d_mu, d_tau, BrMutoEGamma, BrTautoEGamma, BrTautoMuGamma           &
   & , BrMu3e, BrTau3e, BrTau3Mu, BR_Z_e_mu, BR_Z_e_tau, BR_Z_mu_tau, BtoSNuNu &
   & , epsK, DeltaMK2, K0toPi0NuNu, KptoPipNuNu, Bd_ll(3), R_Bu_TauNu
 Complex(dp) :: DeltaMBd, DeltaMBs
 !------------------------------------
 ! auxiliary parameters
 !------------------------------------
 Logical :: CalcTBD
 Integer :: kont, i1
 Real(dp) :: m_Gut, ratioWoM, Q_in
 !------------------------------------------------------
 ! internal information on particle identities
 !------------------------------------------------------
  Integer :: id_gl, id_ph, id_grav
  Integer, Dimension(1) :: id_Z, id_W
  Integer, Dimension(3) :: id_nu, id_l, id_d, id_u
 !------------------------------------------------------
 ! Input and output files
 !------------------------------------------------------
  Character(len=60) :: inputFileName, outputFileName

 Iname = 1
 NameOfUnit(1) = "SPheno3"
 !--------------------------------------------------------------------------
 ! set all parameters and low energy observables to zero
 !--------------------------------------------------------------------------
 Call Set_All_Parameters_0()

 BRbtosgamma = 0._dp
 BToSNuNu = 0._dp
 BrBToSLL = 0._dp
 DeltaMBd = 0._dp
 DeltaMBs = 0._dp
 Bs_ll = 0._dp
 Bd_ll = 0._dp
 BR_Bu_TauNu = 0._dp
 R_Bu_TauNu = 0._dp

 epsK = 0._dp
 DeltaMK2 = 0._dp
 K0toPi0NuNu = 0._dp
 KptoPipNuNu = 0._dp

 a_e = 0._dp
 a_mu = 0._dp
 a_tau = 0._dp
 d_e = 0._dp
 d_mu = 0._dp
 d_tau = 0._dp
 BrMutoEGamma = 0._dp
 BrTautoEGamma = 0._dp
 BrTautoMuGamma = 0._dp
 BrMu3e = 0._dp
 BrTau3e = 0._dp
 BrTau3Mu = 0._dp
 BR_Z_e_mu = 0._dp
 BR_Z_e_tau = 0._dp
 BR_Z_mu_tau = 0._dp
 rho_parameter = 0._dp
 mf_nu = 0
 !--------------------------------------------------------------------------
 ! This routine call routines to
 !   - initializ the error system
 !   - to calculate the constants which are required for the 
 !     calculation of the loop functions
 !   - to get the standard model parameters
 !   - to get SUSY parameters
 ! The following steps are performed to get the parameters and flags
 !   (i) The file LesHouches.in exists containing all necessary information.
 !       In this case the remaining steps are skipped
 !   (ii) Check if Control.in exists to set the error level and
 !        to check if widths and cross sections shall be calculated
 !   (iii) Check if StandardModel.in exists to change the default SM values
 !   (iv) Read the information concerning the SUSY model from the file
 !        HighScale.in
 ! Note that either the file LesHouches.in or the file HighScale.in
 ! must exist.
 !--------------------------------------------------------------------------
 Call ReadingData(kont)
 !---------------------------------------------
 ! parameters for branching ratio calculations
 !---------------------------------------------
  ! relative precision for the calculation of phase space integrals
  epsI = 1.e-5_dp
  deltaM = 1.e-3_dp ! if mass/phase space is smaller, than mass is treated as 0
  CalcTBD = .False. ! if .True. than calculation of 3-body decays is enforced
  ratioWoM = 0._dp ! 1.e-4_dp

 If ((HighScaleModel.Eq."NMSSM").And.(kont.Eq.0)) Then ! NMSSM model

  Call Model_NMSSM(m32, Grav_fac, F_GMSB, Ecms, Pm, Pp, ISR, Beam           &
   & , SigSup , SigSdown, SigSle, SigSn, SigC, SigChi0, SigS0, SigSP, SigHp &
   & , kont)

 Else If (((HighScaleModel2.Eq."RPexplicit").Or.(Add_Rparity)).And.(kont.Eq.0)) Then ! bilinear RP

  Call Model_bilinear_Rparity(add_Rparity, HighScaleModel, delta_mass, epsI     &
       & , deltaM, ratioWoM, m32, grav_fac, CalcTBD, Ecms, Pm, Pp, ISR, Beam    &
       & , SigSup , SigSdown, SigC, SigChi0, SigS0, SigSP, SigHp, M_GUT, kont)
  
 Else If (RP_trilinear) Then ! trilinear RP
  HighScaleModel = "RPexplicit"
  Call Model_trilinearL_Rparity(delta_mass, epsI, deltaM, ratioWoM    &
       & , m32, grav_fac, CalcTBD, Ecms, Pm, Pp, ISR, Beam            &
       & , SigSup , SigSdown, SigC, SigChi0, SigS0, SigSP, SigHp, kont)
  
 Else If (kont.Eq.0) Then  ! models with MSSM particle content
                           ! at the electroweak scale 
  Call Initialize_MSSM(GenerationMixing, id_gl, id_ph, id_Z, id_W      &
                    & , id_nu, id_l, id_d, id_u, id_grav)
  !---------------------------------------------------------------------------
  ! calculation of the spectrum, the following parameters can be changed
  ! with the help of the SLHA input file LesHouches.in or the file Control.in
  ! 
  !---------------------------------------------------------------------------
   Call CalculateSpectrum(n_run, delta_mass, WriteOut, kont, tanb           &
    & , vevSM, ChiPm%m, U, V, Chi0%m, N, S0%m, S0%m2, RS0, P0%m, P0%m2, RP0 &
    & , Spm%m, Spm%m2, RSpm, Sdown%m, Sdown%m2, RSdown, Sup%m, Sup%m2, RSup &
    & , Slepton%m, Slepton%m2, RSlepton, Sneut%m, Sneut%m2, RSneut, Glu%m   &
    & , PhaseGlu, gauge, uL_L, uL_R, uD_L, uD_R, uU_L, uU_R, Y_l, Y_d, Y_u  &
    & , Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, m_GUT)

   Q_in = Sqrt(GetRenormalizationScale())

   If (kont.Eq.0) Then
    Chi0%m2 = Chi0%m**2
    ChiPm%m2 = ChiPm%m**2
    Glu%m2 = Glu%m**2

    Call Low_Energy_Constraints_MSSM(Q_in, gauge, Y_l, Y_d, Y_u, A_l, A_d, A_u &
     & , Mi, mu, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, B, tanb_Q, P0%m2, S0%m2   &
     & , Spm%m2, CKM, kont, GenerationMixing, rho_parameter, DeltaMBd          &
     & , DeltaMBs, BRbtosgamma, Bs_ll, Bd_ll, BrBToSLL, BtoSNuNu, BR_Bu_TauNu  &
     & , R_Bu_TauNu, epsK, DeltaMK2, K0toPi0NuNu, KptoPipNuNu                  &
     & , a_e, a_mu, a_tau, d_e, d_mu, d_tau, BrMutoEGamma, BrTautoEGamma       &
     & , BrTautoMuGamma, BrMu3e, BrTau3e, BrTau3Mu, BR_Z_e_mu, BR_Z_e_tau      &
     & , BR_Z_mu_tau)
   End If

   !----------------------------------------------
   ! reorder state identification if necessary
   !----------------------------------------------
   If ((.Not.GenerationMixing).And.(kont.Eq.0)) Then
    Call Swap_Order_Sf(RSlepton(1,1), Slepton(1)%id, Slepton(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSlepton(3,3), Slepton(3)%id, Slepton(4)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSdown(1,1), Sdown(1)%id, Sdown(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSdown(3,3), Sdown(3)%id, Sdown(4)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSup(1,1), Sup(1)%id, Sup(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSup(3,3), Sup(3)%id, Sup(4)%id, id_p, c_name)    
   End If
 !-------------------------------------------------------------------
 ! Calculation of the branching ratios and widths provided L_BR is
 ! set .TRUE. (default) and that the routine Sugra has finished
 ! correctly (kont.eq.0) 
 !-------------------------------------------------------------------
  If ((L_BR).And.(kont.Eq.0)) Then
   Call CalculateBR(n_nu, id_nu, n_l, id_l, n_d, id_d, n_u, id_u, n_Z, id_Z     &
     & , n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_g, n_s0, n_p0, n_Spm  &
     & , id_grav, id_gl, id_ph, gauge, Glu, PhaseGlu, ChiPM, U, V, Chi0, N      &
     & , Sneut, RSneut, Slepton, RSlepton, Sup, RSup, Sdown, RSdown, uL_L, uL_R &
     & , uD_L, uD_R, uU_L, uU_R, S0, RS0, P0, RP0, Spm, RSpm, epsI, deltaM      &
     & , CalcTBD, ratioWoM, Y_d, A_d, Y_l, A_l, Y_u, A_u, mu, vevSM, F_Gmsb     &
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

  If ((L_CS).And.(kont.Eq.0)) Then
   Ylp = Y_l / gauge(2)
   Do i1=1,p_max
    If (Ecms(i1).Eq.0._dp) Exit

    Call CalculateCrossSections(Ecms(i1), Pm(i1), Pp(i1), ISR(i1), Beam(i1)  &
           & , "Tesla800", Sup%m, RSup, mf_u, Sdown%m, RSdown, mf_d, Glu%m      &
           & , SigSup(i1,:,:), SigSdown(i1,:,:), Slepton%m, RSlepton, Ylp     &
           & , Sneut%m, RSneut, SigSle(i1,:,:), SigSn(i1,:,:), ChiPm%m, U, V       &
           & , Chi0%m, N, SigC(i1,1:2,1:2), SigChi0(i1,1:4,1:4), S0%m, RS0, vevSM &
           & , P0%m, RP0, Spm%m, RSpm, SigS0(i1,1:2), SigSP(i1,1:2,1)          &
           & , SigHp(i1,1,1) )
   End Do

  End If

 End If

 If (kont.Ne.0) Call WriteNumberOfErrors(ErrCan)

 Call LesHouches_Out(67, kont, id_p, c_name, HighScaleModel, M_GUT         &
      & , BRbtosgamma, Bs_ll, Bd_ll, DeltaMBd, DeltaMBs, BrBToSLL      &
      & , BtoSNuNu, BR_Bu_TauNu, R_Bu_TauNu                                &
      & , a_e, a_mu, a_tau, d_e, d_mu, d_tau, BrMuToEGamma, BrTauToEGamma  &
      & , BrTauToMuGamma, BrMu3e, BrTau3e, BrTau3Mu, BR_Z_e_mu, BR_Z_e_tau &
      & , BR_Z_mu_tau, epsK, DeltaMK2, K0toPi0NuNu, KptoPipNuNu             &
      & , Rho_parameter, Ecms, Pm, Pp, ISR, SigSup, SigSdown, SigSle       &
      & , SigSn, SigChi0, SigC, SigS0, SigSP, SigHp, f_name=Trim(outputFileName))
 !------------------------------------------------------------
 ! programs like micrOmegas do not yet use flavour mixing, in
 ! this case a modified SLHA file is needed
 !------------------------------------------------------------
 If ((Write_SLHA1).And.(kont.Eq.0)) Call WriteSPhenoOutputLHA1(35, M_GUT)

 Call closing() ! closes the files
 If ((kont.Ne.0).And.Non_Zero_Exit) Stop 99

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

 Subroutine CalculateSpectrum(n_run, delta, WriteOut, kont, tanb, vevSM     &
     & , mC, U, V, mN, N, mS0, mS02, RS0, mP0, mP02, RP0, mSpm, mSpm2, RSpm &
     & , mSdown, mSdown2, RSdown, mSup, mSup2, RSup, mSlepton, mSlepton2    &
     & , RSlepton, mSneut, mSneut2, RSneut, mGlu, PhaseGlu, gauge           &
     & , uL_L, uL_R, uD_L, uD_R, uU_L, uU_R, Y_l, Y_d, Y_u                  &
     & , Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B       &
     & , m_GUT)
 Implicit None
  Integer, Intent(in) :: n_run
  Logical, Intent(in) :: WriteOut
  Real(dp), Intent(in) :: delta
  Integer, Intent(inout) :: kont
  Real(dp), Intent(inout) :: gauge(3), M2_H(2), tanb, vevSM(2), mP0(2) &
     & , mP02(2), RP0(2,2), mGlu, mSpm(2), mSpm2(2)
  Complex(dp), Intent(inout) :: Y_l(3,3), Y_d(3,3), Y_u(3,3), mu, B, Mi(3)  &
     & , A_l(3,3), A_d(3,3), A_u(3,3), M2_E(3,3), M2_L(3,3), M2_D(3,3)      &
     & , M2_U(3,3), M2_Q(3,3), PhaseGlu, RSpm(2,2)
  Real(dp), Intent(inout) :: mS0(2), mS02(2), RS0(2,2) &
     & , mC(2), mN(4), mSneut(3), mSneut2(3), mSlepton(6), mSlepton2(6) &
     & , mSdown(6), mSdown2(6), mSup(6), mSup2(6), m_GUT
  Complex(dp), Intent(inout) :: U(2,2), V(2,2), N(4,4), Rsneut(3,3)  &
     & , RSlepton(6,6), RSdown(6,6), RSup(6,6)
  Complex(dp), Intent(inout), Dimension(3,3) :: uL_L, uL_R, uD_L, uD_R, uU_L &
     & , uU_R

  Real(dp) :: mZ2_t, mW2_t, Scale_Q, g, gp, dt, tz, g1(57), g58(58), g214(214) &
       & , mN2(4), mC2(2), g2(213), mglu_T, mP0_T(2), mP02_T(2), RP0_T(2,2)
  Real(dp) :: mC_T(2), mN_T(4), mS0_T(2), mSpm_T(2), mSup_T(6), mSdown_T(6)  &
    & , mSlepton_T(6), mSneut_T(3), mSup2_T(6), mSdown2_T(6), mSlepton2_T(6) &
    & , mSneut2_T(3), mS02_T(2), mSpm2_T(2), RS0_T(2,2), mC2_T(2), mN2_T(4)  &
    & , mass_new(32), mass_old(32), diff_m(32), sinW2, mf3(3), vev, Abs_mu2
  Complex(dp) :: U_T(2,2), V_T(2,2), N_T(4,4), RSpm_T(2,2), RSdown_T(6,6)  &
    & , RSup_T(6,6), RSlepton_T(6,6), RSneut_T(3,3), mu_T, B_T
  Complex(dp), Dimension(3,3) :: Al_s, M2E_s, M2L_s, Ad_s, Au_s, M2D_s, M2Q_s &
    & , M2U_s
  Integer :: i1, i2, i3, ierr
  Logical :: Converge, UseFixedScale

  !------------------------------------------------------------------
  ! Performing a first, very rough calculation of the parameters
  ! using 1-loop RGEs and tree-level boundary conditions for gauge and
  ! Yukawa couplings at m_Z for the parameters.
  ! These parameters are used to get a first set of tree-level masses
  ! which are needed for the first of computation of the SUSY thresholds
  ! to gauge and Yukawa couplings. In case of the general MSSM the
  ! parameters are already given.
  !------------------------------------------------------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "CalculateSpectrum"

  kont = 0
  If ( (HighScaleModel(1:4).Eq."MSSM").Or.(HighScaleModel.Eq."pMSSM") )  Then
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

  Else ! high scale model

   If (Len(Trim(Old_data)).Ne.0) Then
    Call FirstGuess(phase_mu, tanb, Mi, M2_E, M2_L, A_l, M2_D   &
           & , M2_Q, M2_U, A_d, A_u, mu, B, M2_H, gp, g, Y_l  &
           & , Y_d, Y_u, vevSM, mP02, mP0, kont, .True., delta, Trim(Old_data) )
   Else 
    Call FirstGuess(phase_mu, tanb, Mi, M2_E, M2_L, A_l, M2_D   &
           & , M2_Q, M2_U, A_d, A_u, mu, B, M2_H, gp, g, Y_l  &
           & , Y_d, Y_u, vevSM, mP02, mP0, kont)
   End If
  End If

  If (External_Spectrum) Then
   Iname = Iname - 1
   Return ! using the externaly given spectrum
  End If

  If (HighScaleModel.Eq."MSSMtree") Then 
   If (GenerationMixing) Then
    Call Switch_from_superCKM(Y_d, Y_u, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm &
              & , M2U_sckm, A_d, A_u, M2_D, M2_Q, M2_U, .False. )
    If (Maxval(Abs(MatNu)).Gt.0._dp) Then
     Call Switch_from_superPMNS(Y_l, MatNu, Al_pmns, M2E_pmns, M2L_pmns &
               & , A_l, M2_E, M2_L, .False. )
    Else
     Call Switch_from_superPMNS(Y_l, id3C, Al_pmns, M2E_pmns, M2L_pmns &
               & , A_l, M2_E, M2_L, .False. )
    End If
   End If
   !------------------------------------------------------------------
   ! last flag implies that Higgs masses at 1-loop effective potential
   !------------------------------------------------------------------
   Call TreeMasses(gp, g, vevSM, Mi(1), Mi(2), Mi(3), mu, B        &
        &, tanb, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d       &
        &, A_u, Y_d, Y_u, mGlu_T, PhaseGlu, mC, mC2, U, V, mN      &
        &, mN2, N, mSneut, mSneut2, Rsneut, mSlepton, mSlepton2    &
        &, RSlepton, mSdown, mSdown2, RSdown, mSup, mSup2, RSup    &
        &, mP0_T, mP02_T, RP0_T, mS0, mS02, RS0, mSpm, mSpm2, RSpm &
        &, GenerationMixing, kont, .True.) 

   tanb_Q = tanb
   mA2_Q = mP02(2)
   vev_Q = Sqrt(vevSM(1)**2 + vevSM(2)**2)
   rho_parameter = 0

  Else 

   Call TreeMasses(gp, g, vevSM, Mi(1), Mi(2), Mi(3), mu, B        &
        &, tanb, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d       &
        &, A_u, Y_d, Y_u, mGlu_T, PhaseGlu, mC, mC2, U, V, mN      &
        &, mN2, N, mSneut, mSneut2, Rsneut, mSlepton, mSlepton2    &
        &, RSlepton, mSdown, mSdown2, RSdown, mSup, mSup2, RSup    &
        &, mP0_T, mP02_T, RP0_T, mS0, mS02, RS0, mSpm, mSpm2, RSpm &
        &, GenerationMixing, kont, .False.) ! tree-level Higgs mass

   mass_old(1) = Abs(mglu_T)
   mass_old(2:3) = Abs(mC)
   mass_old(4:7) = Abs(mN)
   mass_old(8:9) = mS0
   mass_old(10) = mP0_T(2)  
   mass_old(11) = mSpm(2)  
   mass_old(12:17) = mSup   
   mass_old(18:23) = mSdown   
   mass_old(24:29) = mSlepton 
   mass_old(30:32) = mSneut
  End If

  If (kont.Ne.0) Then
   Iname = Iname - 1
   Return
  End If

   ! In the SPA convention the the renormalization scale is fixed with 1 TeV
  If (SPA_Convention) Call SetRGEScale(1.e3_dp**2)
  Scale_Q = Sqrt(GetRenormalizationScale())
  If (Scale_Q.Eq.1._dp) Then ! in this case no scale has been set and one
   UseFixedScale = .False.   ! has to use the square root of the stop masses
  Else
   UseFixedScale = .True.
  End If

  If (HighScaleModel.Eq."MSSM") Then
   ! MSSM parameters and masses at loop level, all parameters are given at
   ! scale Q, mu and m_A(tree) serve as input in the Higgs sector
   ! calculate first gauge and Yukawa in DR-scheme at m_Z
   kont = 0
   If (UseNewBoundaryEW) Then
    Call BoundaryEW(1, mZ, tanb, Mi, Al_pmns, Ad_sckm, Au_sckm, M2E_pmns &
    & , M2L_pmns, M2D_sckm, M2Q_sckm, M2U_sckm, mu, mP0_T(2), delta      &
    & , GenerationMixing, .True., mZ2_t, mW2_t, g1, kont)
   Else
    Call BoundaryEW(1, vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T &
     & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup  &
     & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut, uU_L, uU_R     &
     & , uD_L, uD_R, uL_L, uL_R, id3C, mGlu_T, PhaseGlu, mZ2_t, mW2_t   &
     & , delta, g1, kont)
   End If

   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If

   converge = .False.
   If (.Not.UseFixedScale) Then
    If (UseNewScale) Then ! note, that is actually the scale squared
     Scale_Q = CalcScale_from_Stops(Real(M2_U(3,3),dp), Real(M2_Q(3,3),dp)   &
               & , Y_u(3,3), A_u(3,3), vevSM, mu, gauge(2), gauge(1) )
    Else
     Scale_Q = product_stop_masses(mSup, RSup, GenerationMixing)
    End If
    tz = SetRenormalizationScale(scale_Q)
    Scale_Q = Sqrt(Scale_Q)
   End If

   Do i1=1,n_run
    tz = Log(Scale_Q/mZ)
    dt = tz / 50._dp
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

     If (.Not.l_Al) A_l = AoY_l * Y_l
     If (.Not.l_Ad) A_d = AoY_d * Y_d
     If (.Not.l_Au) A_u = AoY_u * Y_u

     Call FermionMass(Y_u, sqrt2, mf3, uU_L, uU_R, ierr)
     Call FermionMass(Y_d, sqrt2, mf3, uD_L, uD_R, ierr)

     M2_D = Matmul( Matmul( Transpose(Conjg(uD_R)), M2D_sckm), uD_R)
     M2_U = Matmul( Matmul( Transpose(Conjg(uU_R)), M2U_sckm), uU_R)
     M2_Q = Matmul( Matmul( Transpose(Conjg(uU_L)), M2Q_sckm), uU_L)

    Else

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

    Call LoopMassesMSSM_3(tanb_mZ, tanb_Q, gauge, Y_l, Y_d, Y_u, Mi, A_l  &
       & , A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, mu, B, 0.1_dp*delta    &
       & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP0, mP02, RP0      &
       & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2        &
       & , RSup, mSlepton, mSlepton2, RSlepton, mSneut, mSneut2           &
       & , RSneut, mGlu, PhaseGlu, M2_H, kont)

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

    If (Maxval(diff_m).Lt.0.1_dp*delta) Then
     converge = .True.
     Exit
    Else
     mass_old = mass_new
    End If

    If (.Not.UseFixedScale) Then
     If (UseNewScale) Then ! note, that is actually the scale squared
      Scale_Q = CalcScale_from_Stops(Real(M2_U(3,3),dp), Real(M2_Q(3,3),dp)   &
                & , Y_u(3,3), A_u(3,3), vevSM, mu, gauge(2), gauge(1) )
     Else
      Scale_Q = product_stop_masses(mSup, RSup, GenerationMixing)
     End If
     tz = SetRenormalizationScale(scale_Q)
     Scale_Q = Sqrt(Scale_Q)
    End If

    !----------------------------------------------
    ! evolve parameters down to m_Z
    !----------------------------------------------
    gauge(1) = Sqrt(5._dp/3._dp) * gauge(1)
    Y_l = Transpose(Y_l)
    Y_d = Transpose(Y_d)
    Y_u = Transpose(Y_u)
    A_l = Transpose(A_l)
    A_d = Transpose(A_d)
    A_u = Transpose(A_u)
    Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u       &
                  & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g2)

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
        & , A_d_mZ, A_u_mZ, M2_E_mZ, M2_L_mZ, M2_D_mZ, M2_Q_mZ, M2_U_mZ    &
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
    A_u_mZ = Transpose(A_u_mZ)
    !-----------------------------------------
    ! use consistently running parameters
    !-----------------------------------------
    g2(1) = Sqrt(3._dp/5._dp) * g2(1)
    sinW2 = g2(1)**2 / (g2(1)**2 +g2(2)**2)
    vev =  2._dp * Sqrt(mZ2_t / (g2(1)**2 +g2(2)**2))
    vevSM(1) = vev / Sqrt(1._dp + tanb_mZ**2)
    vevSM(2) = tanb_mZ * vevSM(1)
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

    Call TreeMasses(g2(1), g2(2), vevSM, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3), mu_T  &
      & , B_T, tanb_mZ, M2_E_mZ, M2_L_mZ, A_l_mZ, Y_l_mZ, M2_D_mZ, M2_U_mZ   &
      & , M2_Q_mZ, A_d_mZ, A_u_mZ, Y_d_mZ, Y_u_mZ, mGlu_T, PhaseGlu, mC_T    &
      & , mC2_T, U_T, V_T, mN_T, mN2_T, N_T, mSneut_T, mSneut2_T, Rsneut_T   &
      & , mSlepton_T, mSlepton2_T, RSlepton_T, mSdown_T, mSdown2_T, RSdown_T &
      & , mSup_T, mSup2_T, RSup_T, mP0_T, mP02_T, RP0_T, mS0_T, mS02_T       &
      & , RS0_T, mSpm_T, mSpm2_T, RSpm_T, GenerationMixing, kont, .False.)

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
     mGlu_T = mGlu
     mP0_T = mP0
     mP02_T = mP02
     RP0_T = RP0
     YukScen = 2 ! using pole masses for boundary conditions at mZ
    End If

    kont = 0
    If (UseNewBoundaryEW) Then
     If (GenerationMixing) then
      Call Switch_to_superCKM(Y_d_mZ, Y_u_mZ, A_d_mZ, A_u_mZ, M2_D_mZ, M2_Q_mZ &
       & , M2_U_mZ, Ad_s, Au_s, M2D_s, M2Q_s, M2U_s, .False.,ckm_out=uU_L )
      If (Maxval(Abs(MatNu)).Gt.0._dp) Then
       Call Switch_from_superPMNS(Y_l_mZ, A_l_mZ, M2_E_mZ, M2_L_mZ, MatNu &
                 & , Al_s, M2E_s, M2L_s, .False. )
      Else
       Call Switch_from_superPMNS(Y_l_mZ, id3C, A_l_mZ, M2_E_mZ, M2_L_mZ &
                & , Al_s, M2E_s, M2L_s, .False. )
      End If
     Else
      Ad_s = A_d_mZ
      Au_s = A_u_mZ
      Al_s = A_l_mZ
      M2D_s = M2_D_mZ
      M2Q_s = M2_Q_mZ
      M2U_s = M2_U_mZ
      M2E_s = M2_E_mZ
      M2L_s = M2_L_mZ
     End If

     Call BoundaryEW(i1+1, mZ, tanb_mZ, Mi_mZ, Al_s, Ad_s, Au_s &
      & , M2E_s, M2L_s, M2D_s, M2Q_s, M2U_s, mu_mZ, mP0_T(2)    &
      & , 0.1_dp*delta, GenerationMixing, .True., mZ2_t, mW2_t, g1, kont)

     uU_R = id3C
     uD_L = id3C
     uD_R = id3C
     uL_L = id3C
     uL_R = id3C
    Else 
     Call BoundaryEW(i1+1,vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T &
       & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup   &
       & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut, uU_L, uU_R      &
       & , uD_L, uD_R , uL_L, uL_R, id3C, mGlu_T, PhaseGlu, mZ2_t, mW2_t   &
       & , delta, g1, kont)
    End If

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If

   End Do

   If ((kont.Eq.0).And.(.Not.converge)) Then
    Write(ErrCan,*) 'Problem in subroutine CalculateSpectrum!!'
    Write(ErrCan,*) "After",n_run,"iterations no convergence found"
    kont = -1200
   End If

  Else If (HighScaleModel.Eq."MSSM1") Then
   ! MSSM parameters and masses at loop level, all parameters are given at
   ! scale Q, m^2_(H_i) serve as input in the Higgs sector
   ! calculate first gauge and Yukawa in DR-scheme at m_Z
   kont = 0
   If (UseNewBoundaryEW) Then
    Call BoundaryEW(1, mZ, tanb, Mi, Al_pmns, Ad_sckm, Au_sckm, M2E_pmns &
    & , M2L_pmns, M2D_sckm, M2Q_sckm, M2U_sckm, mu, mP0_T(2), delta      &
    & , GenerationMixing, .True., mZ2_t, mW2_t, g1, kont)
   Else
    Call BoundaryEW(1, vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T, mSpm &
     & , mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
     & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut, uU_L, uU_R           &
     & , uD_L, uD_R, uL_L, uL_R, id3C, mGlu_T, PhaseGlu, mZ2_t, mW2_t         &
     & , delta, g1, kont)
   End If

   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If

   converge = .False.
   tanb_mZ = tanb ! first gues, justified because tanb runs weakly

   If (.Not.UseFixedScale) Then
    If (UseNewScale) Then ! note, that is actually the scale squared
     Scale_Q = CalcScale_from_Stops(Real(M2_U(3,3),dp), Real(M2_Q(3,3),dp)   &
               & , Y_u(3,3), A_u(3,3), vevSM, mu, gauge(2), gauge(1) )
    Else
     Scale_Q = product_stop_masses(mSup, RSup, GenerationMixing)
    End If
    tz = SetRenormalizationScale(scale_Q)
    Scale_Q = Sqrt(Scale_Q)
   End If

   Do i1=1,n_run
    !----------------------------------------------
    ! evolve first parameters up to Q
    !----------------------------------------------
    tz = Log(Scale_Q/mZ)
    dt = tz / 50._dp
    If (tanb_in_at_Q) Then
     Call odeint(g1, 57, 0._dp, tz, 0.1_dp*delta, dt, 0._dp, rge57, kont)
    Else  ! include tanb(beta) in running
     g58(1:57) = g1
     g58(58) = Log(tanb_mZ)
     Call odeint(g58, 58, 0._dp, tz, 0.1_dp*delta, dt, 0._dp, rge58, kont)
     g1 = g58(1:57)
     tanb_Q = Exp(g58(58))
    End If

    Call GToCouplings(g1, gauge, Y_l, Y_d, Y_u)
    gauge(1) = Sqrt(3._dp/5._dp) * gauge(1)

    If (GenerationMixing) Then
     Y_l = Transpose(Y_l) 
     Y_d = Transpose(Y_d) 
     Y_u = Transpose(Y_u) 
     If (At_save.Ne.ZeroC) Au_sckm(3,3) = At_save * Y_u(3,3)
     If (Ab_save.Ne.ZeroC) Ad_sckm(3,3) = Ab_save * Y_d(3,3)
     If (Atau_save.Ne.ZeroC) Al_pmns(3,3) = Atau_save * Y_l(3,3)
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

    Call LoopMassesMSSM(delta, tanb_mZ, tanb_Q, gauge, Y_l, Y_d, Y_u, Mi   &
       & , A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, phase_mu, mu &
       & , B, i1, uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                       &
       & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP0, mP02, RP0       &
       & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2         &
       & , RSup, mSlepton, mSlepton2, RSlepton, mSneut, mSneut2            &
       & , RSneut, mGlu, PhaseGlu, kont)

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

    If (Maxval(diff_m).Lt.0.1_dp*delta) Then
     converge = .True.
     Exit
    Else
     mass_old = mass_new
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
    A_u_mZ = Transpose(A_u_mZ)
    !-----------------------------------------
    ! use consistently running parameters
    !-----------------------------------------
    g2(1) = Sqrt(3._dp/5._dp) * g2(1)
    sinW2 = g2(1)**2 / (g2(1)**2 +g2(2)**2)
    vev =  2._dp * Sqrt(mZ2_t / (g2(1)**2 +g2(2)**2))
    vevSM(1) = vev / Sqrt(1._dp + tanb_mZ**2)
    vevSM(2) = tanb_mZ * vevSM(1)
    
    If (.Not.UseFixedScale) Then
     If (UseNewScale) Then ! note, that is actually the scale squared
      Scale_Q = CalcScale_from_Stops(Real(M2_U(3,3),dp), Real(M2_Q(3,3),dp)   &
                & , Y_u(3,3), A_u(3,3), vevSM, mu, gauge(2), gauge(1) )
     Else
      Scale_Q = product_stop_masses(mSup, RSup, GenerationMixing)
     End If
     tz = SetRenormalizationScale(scale_Q)
     Scale_Q = Sqrt(Scale_Q)
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

    Call TreeMassesMSSM2(g2(1), g2(2), vevSM, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
      & , mu_T, B_T , tanb_mZ, M2_E_mZ, M2_L_mZ, A_l_mZ, Y_l_mZ, M2_D_mZ   &
      & , M2_U_mZ, M2_Q_mZ, A_d_mZ, A_u_mZ, Y_d_mZ, Y_u_mZ, uU_L, uU_R     &
      & , uD_L, uD_R, uL_L, uL_R, mGlu_T, PhaseGlu, mC_T, mC2_T, U_T, V_T  &
      & , mN_T, mN2_T, N_T, mSneut_T, mSneut2_T, Rsneut_T, mSlepton_T      &
      & , mSlepton2_T, RSlepton_T, mSdown_T, mSdown2_T, RSdown_T, mSup_T   &
      & , mSup2_T, RSup_T, mP0_T, mP02_T, RP0_T, mS0_T, mS02_T, RS0_T      &
      & , mSpm_T, mSpm2_T, RSpm_T, mZ2_t, mW2_t, GenerationMixing, kont    &
      & , .False., .False.)

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If

    If (Min(Minval(mSup2_T), Minval(mSdown2_T), Minval(mSlepton2_T)  &
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
    If (UseNewBoundaryEW) Then
     If (GenerationMixing) then
      Call Switch_to_superCKM(Y_d_mZ, Y_u_mZ, A_d_mZ, A_u_mZ, M2_D_mZ, M2_Q_mZ &
       & , M2_U_mZ, Ad_s, Au_s, M2D_s, M2Q_s, M2U_s, .False.,ckm_out=uU_L )
      If (Maxval(Abs(MatNu)).Gt.0._dp) Then
       Call Switch_from_superPMNS(Y_l_mZ, A_l_mZ, M2_E_mZ, M2_L_mZ, MatNu &
                 & , Al_s, M2E_s, M2L_s, .False. )
      Else
       Call Switch_from_superPMNS(Y_l_mZ, id3C, A_l_mZ, M2_E_mZ, M2_L_mZ &
                 & , Al_s, M2E_s, M2L_s, .False. )
      End If
     Else
      Ad_s = A_d_mZ
      Au_s = A_u_mZ
      Al_s = A_l_mZ
      M2D_s = M2_D_mZ
      M2Q_s = M2_Q_mZ
      M2U_s = M2_U_mZ
      M2E_s = M2_E_mZ
      M2L_s = M2_L_mZ
     End If

     Call BoundaryEW(i1+1, mZ, tanb_mZ, Mi_mZ, Al_s, Ad_s, Au_s &
      & , M2E_s, M2L_s, M2D_s, M2Q_s, M2U_s, mu_mZ, mP0_T(2)    &
      & , 0.1_dp*delta, GenerationMixing, .True., mZ2_t, mW2_t, g1, kont)

     uU_R = id3C
     uD_L = id3C
     uD_R = id3C
     uL_L = id3C
     uL_R = id3C
    Else 
     Call BoundaryEW(i1+1,vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T &
       & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup   &
       & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut, uU_L, uU_R      &
       & , uD_L, uD_R , uL_L, uL_R, id3C, mGlu_T, PhaseGlu, mZ2_t, mW2_t   &
       & , delta, g1, kont)
    End If

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

  Else If (HighScaleModel.Eq."pMSSM") Then
   ! calculate first gauge and Yukawa in DR-scheme at m_Z

   If (UseNewBoundaryEW) Then
    Call BoundaryEW(1, mZ, tanb, Mi, Al_pmns, Ad_sckm, Au_sckm, M2E_pmns &
    & , M2L_pmns, M2D_sckm, M2Q_sckm, M2U_sckm, mu, mP0_T(2), delta      &
    & , GenerationMixing, .True., mZ2_t, mW2_t, g1, kont)
   Else
    Call BoundaryEW(1, vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T, mSpm &
     & , mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
     & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut, uU_L, uU_R           &
     & , uD_L, uD_R, uL_L, uL_R, id3C, mGlu_T, PhaseGlu, mZ2_t, mW2_t         &
     & , delta, g1, kont)
   End If

   If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If

   converge = .False.
   mP0(1) = mZ
   mP02(1) = mZ2
   If (.Not.UseFixedScale) Then
    If (UseNewScale) Then ! note, that is actually the scale squared
     Scale_Q = CalcScale_from_Stops(Real(M2_U(3,3),dp), Real(M2_Q(3,3),dp)   &
               & , Y_u(3,3), A_u(3,3), vevSM, mu, gauge(2), gauge(1) )
    Else
     Scale_Q = product_stop_masses(mSup, RSup, GenerationMixing)
    End If
    tz = SetRenormalizationScale(scale_Q)
    Scale_Q = Sqrt(Scale_Q)
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
               & , M2U_sckm, A_d, A_u, M2_D, M2_Q, M2_U, .False.)

     If (Maxval(Abs(MatNu)).Gt.0._dp) Then
      Call Switch_from_superPMNS(Y_l, MatNu, Al_pmns, M2E_pmns, M2L_pmns &
                & , A_l, M2_E, M2_L, .False. )
     Else
      Call Switch_from_superPMNS(Y_l, id3C, Al_pmns, M2E_pmns, M2L_pmns &
                & , A_l, M2_E, M2_L, .False. )
     End If
    Else ! .not. GenerationMixing
     A_l = ZeroC
     A_d = ZeroC
     A_u = ZeroC
     Do i2=1,3
      If (.Not.l_Al) A_l(i2,i2) = AoY_l(i2,i2) * Y_l(i2,i2)
      If (.Not.l_Ad) A_d(i2,i2) = AoY_d(i2,i2) * Y_d(i2,i2)
      If (.Not.l_Au) A_u(i2,i2) = AoY_u(i2,i2) * Y_u(i2,i2)
     End Do
    End If

    kont = 0
    Call LoopMassesMSSM_2(delta, tanb_mZ, tanb_Q, gauge, Y_l, Y_d, Y_u &
       & , Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, mu         &
       & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP02, RP0        &
       & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2     &
       & , RSup, mSlepton, mSlepton2, RSlepton, mSneut, mSneut2        &
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
    gauge_MZ(1) = Sqrt(3._dp/5._dp) * gauge_MZ(1)
    sinW2 = g2(1)**2 / (g2(1)**2 +g2(2)**2)
    vev =  2._dp * Sqrt(mZ2_t / (g2(1)**2 +g2(2)**2))
    vevSM(1) = vev / Sqrt(1._dp + tanb_mZ**2)
    vevSM(2) = tanb_mZ * vevSM(1)
 
    If (.Not.UseFixedScale) Then
     If (UseNewScale) Then ! note, that is actually the scale squared
      Scale_Q = CalcScale_from_Stops(Real(M2_U(3,3),dp), Real(M2_Q(3,3),dp)   &
                & , Y_u(3,3), A_u(3,3), vevSM, mu, gauge(2), gauge(1) )
     Else
      Scale_Q = product_stop_masses(mSup, RSup, GenerationMixing)
     End If
     tz = SetRenormalizationScale(scale_Q)
     Scale_Q = Sqrt(Scale_Q)
    End If
    !-----------------------------------
    ! mu and B parameters at tree-level
    !-----------------------------------
    Abs_Mu2 = (M2_H_mZ(2) * tanb_mZ**2 - M2_H_mZ(1) ) / (1._dp - tanb_mZ**2) &
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
    If (UseNewBoundaryEW) Then
     If (GenerationMixing) then
      Call Switch_to_superCKM(Y_d_mZ, Y_u_mZ, A_d_mZ, A_u_mZ, M2_D_mZ, M2_Q_mZ &
       & , M2_U_mZ, Ad_s, Au_s, M2D_s, M2Q_s, M2U_s, .False.,ckm_out=uU_L )
      If (Maxval(Abs(MatNu)).Gt.0._dp) Then
       Call Switch_from_superPMNS(Y_l_mZ, A_l_mZ, M2_E_mZ, M2_L_mZ, MatNu &
                 & , Al_s, M2E_s, M2L_s, .False. )
      Else
       Call Switch_from_superPMNS(Y_l_mZ, id3C, A_l_mZ, M2_E_mZ, M2_L_mZ &
                 & , Al_s, M2E_s, M2L_s, .False. )
      End If
     Else
      Ad_s = A_d_mZ
      Au_s = A_u_mZ
      Al_s = A_l_mZ
      M2D_s = M2_D_mZ
      M2Q_s = M2_Q_mZ
      M2U_s = M2_U_mZ
      M2E_s = M2_E_mZ
      M2L_s = M2_L_mZ
     End If

     Call BoundaryEW(i1+1, mZ, tanb_mZ, Mi_mZ, Al_s, Ad_s, Au_s &
!      & , M2E_s, M2L_s, M2D_s, M2Q_s, M2U_s, mu_T, mP0_T(2)    &
      & , M2E_s, M2L_s, M2D_s, M2Q_s, M2U_s, mu, mP0(2)    &
      & , 0.1_dp*delta, GenerationMixing, .True., mZ2_t, mW2_t, g1, kont)

     uU_R = id3C
     uD_L = id3C
     uD_R = id3C
     uL_L = id3C
     uL_R = id3C
    Else 
     Call BoundaryEW(i1+1,vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T &
      & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup   &
      & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut, uU_L, uU_R      &
      & , uD_L, uD_R , uL_L, uL_R, id3C, mGlu_T, PhaseGlu, mZ2_t, mW2_t   &
      & , delta, g1, kont)

    End If

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

  Else If (HighScaleModel.Eq."MSSMtree") Then
   mP0 = mP0_T
   mP02 = mP02_T
   RP0 = RP0_T
   mglu = mGlu_T
   If (GenerationMixing) Then

    Call FermionMass(Y_l,vevSM(1),mf_l_mZ,uL_L,uL_R,kont)
    Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevSM, mf_d_mZ, uD_L, uD_R &
                                     & , mf_u_mZ, uU_L, uU_R)

   Else
    uL_L = id3C
    uL_R = id3C
    uD_L = id3C
    uD_R = id3C
    uU_L = id3C
    uU_R = id3C
   End If
   
  Else ! high scale models
   mP0 = mP0_T
   mP02 = mP02_T
   RP0 = RP0_T
   mglu = mGlu_T
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

   Call Sugra(delta, vevSM, mC, U, V, mN, N, mS0, mS02, RS0 &
     & , mP0, mP02, RP0, mSpm, mSpm2, RSpm, mSdown, mSdown2 &
     & , RSdown, mSup, mSup2, RSup, mSlepton, mSlepton2     &
     & , RSlepton, mSneut, mSneut2, RSneut, mGlu, PhaseGlu  &
     & , gauge, uL_L, uL_R, uD_L, uD_R, uU_L, uU_R, Y_l     &
     & , Y_d, Y_u, Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D      &
     & , M2_Q, M2_U, M2_H, mu, B, m_GUT, kont, WriteOut, n_run)
  End If

  !--------------------------------------------------------------
  ! Saving the masses in an extra vector in case that the
  ! output for at the mass scale of a specific is required
  !--------------------------------------------------------------
  Q_PDG_out(1) = mGlu
  Q_PDG_out(2:3) = Abs(mC)
  Q_PDG_out(4:7) = Abs(mN)
  Q_PDG_out(8:9) = mS0
  Q_PDG_out(10) = mP0(2)
  Q_PDG_out(11) = mSpm(2)
  Q_PDG_out(12:17) = mSup
  Q_PDG_out(18:23) = mSdown
  Q_PDG_out(24:29) = mSlepton
  Q_PDG_out(30:32) = mSneut
  Q_PDG_out(36) = mf_u(3) 

  Iname = Iname - 1

 End Subroutine CalculateSpectrum


 Subroutine ReadingData(kont)
 !--------------------------------------------------------
 ! reading the input, all routines used can be found in
 ! InputOutput.f90 
 !--------------------------------------------------------
 Implicit None
  Integer, Intent(out) :: kont

  Logical :: file_exists=.False., file_exists2=.False.

  Call get_command_argument(1,inputFileName)
!  inputFileName="LesHouches.in"
  If (Len_trim(inputFileName) == 0) Then
   inputFileName="LesHouches.in"
  Else
   inputFileName=Trim(inputFileName)
   Inquire(file=Trim(inputFileName),exist=file_exists)
   If (.Not.file_exists) Then
    Write(*,*) "Input file ",Trim(inputFileName)," not found!"
    Write(*,*) "Trying to find LesHouches.in"
    Write(ErrCan,*) "Input file ",Trim(inputFileName)," not found!"
    Write(ErrCan,*) "Trying with LesHouches.in"
    inputFileName="LesHouches.in"
   End If
  End If

  Call get_command_argument(2,outputFileName)
!  outputFileName="SPheno.spc"
  If (Len_trim(outputFileName) == 0) Then
   outputFileName="SPheno.spc"
  Else
   outputFileName=Trim(outputFileName)
  End If

  kont = -123456

 !------------------------------------------------------------------
 ! Checked, if the file LesHouches.in exists
 !------------------------------------------------------------------
  Inquire(file="LesHouches.in",exist=file_exists2)
 !---------------------------------------
 !   if yes, use the data from this file
 !---------------------------------------
  If (file_exists) Then
   kont = 0

   Call LesHouches_Input(inputFileName, kont, HighScaleModel &
      &                , Ecms, Pm, Pp, ISR, F_GMSB)

  Else If (file_exists2) Then
   kont = 0

   Call LesHouches_Input(inputFileName, kont, HighScaleModel &
      &                , Ecms, Pm, Pp, ISR, F_GMSB)

  Else

   Write(*,*) "No input file has been found. Please provide one."
   Write(ErrCan,*) "No input file has been found. Please provide one."
   Call TerminateProgram 

  End If

 End Subroutine ReadingData

End Program SPheno

